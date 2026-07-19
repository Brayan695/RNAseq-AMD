# The CNC-VAE model (tf.keras version): a single-hidden-layer variational
# autoencoder over the concatenated clinical+gene-expression vector ("CNC" =
# Concatenated iNputs Cnc-vae). Trains an encoder (input -> latent) and decoder
# (latent -> reconstructed input) jointly; downstream scripts use only the
# trained encoder's latent output as a compact per-patient representation.
#
# NOTE: the PyTorch port in ../../CNC-VAE/models/cncvae.py implements the exact
# same architecture/math without needing the _VAELossLayer workaround below -
# if you're debugging Keras-specific issues, it can help to compare against
# that simpler version.
import numpy as np
import tensorflow as tf
from tensorflow.keras import backend as K
from tensorflow.keras import optimizers
from tensorflow.keras.layers import BatchNormalization as BN, Dense, Input, Lambda, Dropout, Layer
from tensorflow.keras.losses import mse as mean_squared_error
from tensorflow.keras.models import Model

from models.common import mmd, sampling, kl_regu


class _VAELossLayer(Layer):
    """Adds the VAE loss (reconstruction + beta * distance) as a model loss.

    Keras 3 forbids calling backend/TF ops directly on symbolic tensors outside
    of a Layer's call() - so the loss combination has to live in here rather
    than inline in build_model(), which is how the original Keras-2 code did it.
    This layer is otherwise a pass-through (returns `outputs` unchanged); its
    only job is to compute the loss as a side effect via self.add_loss(...).
    """

    def __init__(self, beta, distance_type, **kwargs):
        super().__init__(**kwargs)
        self.beta = beta
        self.distance_type = distance_type

    def call(self, inputs):
        x, outputs, z_mean, z_log_sigma, z = inputs
        if self.distance_type == 'mmd':
            true_samples = K.random_normal(tf.shape(z))  # samples from the target N(0, I)
            distance = mmd(true_samples, z)
        else:
            distance = kl_regu(z_mean, z_log_sigma)

        reconstruction_loss = mean_squared_error(x, outputs)
        vae_loss = K.mean(reconstruction_loss + self.beta * distance)
        self.add_loss(vae_loss)  # registers this as the model's training loss
        return outputs


class CNCVAE:
    """CNC-VAE for Clin+mRNA integration (AMD cohort has no CNA arm, unlike METABRIC)."""

    def __init__(self, args):
        self.args = args
        self.vae = None       # full encoder+decoder model, used for training
        self.encoder = None   # encoder-only model, used for predict()/downstream use

    def build_model(self):
        """Construct the encoder, decoder, and combined VAE model, then compile
        it with the Adam optimizer. Must be called before train()/predict().
        Reads its hyperparameters off self.args (input_size, ds, ls, act,
        dropout, beta, distance)."""
        np.random.seed(42)
        tf.random.set_seed(42)  # fixed seed -> reproducible weight initialization

        # ------------ Input -----------------
        inputs = Input(shape=(self.args.input_size,), name='concat_input')

        # ------------ Encoding Layer -----------------
        # Single hidden layer: Dense -> BatchNorm (no activation-layer dropout
        # in the encoder, only in the decoder below).
        x = Dense(self.args.ds, activation=self.args.act)(inputs)
        x = BN()(x)

        # ------------ Embedding Layer --------------
        # Two parallel Dense heads turn the hidden representation into the mean
        # and log-variance of a Gaussian over the latent space (what makes this
        # a *variational* autoencoder). z_log_sigma is zero-initialized so
        # every sample starts training with unit variance.
        z_mean = Dense(self.args.ls, name='z_mean')(x)
        z_log_sigma = Dense(self.args.ls, name='z_log_sigma', kernel_initializer='zeros')(x)
        # Lambda wraps the sampling() function (see models/common.py) as a
        # proper Keras layer, since Keras 3 doesn't allow raw TF ops outside a
        # layer's call() - this is also why sampling() needs to be registered
        # with @keras.saving.register_keras_serializable.
        z = Lambda(sampling, output_shape=(self.args.ls,), name='z')([z_mean, z_log_sigma])

        self.encoder = Model(inputs, [z_mean, z_log_sigma, z], name='encoder')
        self.encoder.summary()

        # ------------ Decoding Layer -----------------
        # Mirror of the encoder: latent -> Dense -> BatchNorm -> activation ->
        # Dropout -> back out to input_size. No output activation (e.g. no
        # sigmoid), since the reconstruction target isn't bounded to [0, 1]
        # and the loss is plain MSE.
        latent_inputs = Input(shape=(self.args.ls,), name='z_sampling')
        x = Dense(self.args.ds, activation=self.args.act)(latent_inputs)
        x = BN()(x)
        x = Dropout(self.args.dropout)(x)

        # ------------ Out -----------------------
        concat_out = Dense(self.args.input_size)(x)

        decoder = Model(latent_inputs, concat_out, name='decoder')
        decoder.summary()

        # Wire encoder -> decoder into one end-to-end model, using the SAMPLED
        # z (index [2] of the encoder's outputs) rather than z_mean, so the
        # decoder is trained against the stochastic latent code.
        outputs = decoder(self.encoder(inputs)[2])
        outputs = _VAELossLayer(self.args.beta, self.args.distance, name='vae_loss')(
            [inputs, outputs, z_mean, z_log_sigma, z]
        )
        self.vae = Model(inputs, outputs, name='vae_mlp')

        adam = optimizers.Adam(learning_rate=0.001, beta_1=0.9, beta_2=0.999)
        # No explicit `loss=` here - the loss was already registered via
        # self.add_loss() inside _VAELossLayer above.
        self.vae.compile(optimizer=adam)
        self.vae.summary()

    def train(self, s1_train, s2_train, s1_test=None, s2_test=None):
        """Train the VAE. s1_*/s2_* are the two input sources (clinical, gene
        expression) as separate arrays - concatenated here into one input
        matrix, since CNC-VAE trains on the combined vector. Passing
        s1_test/s2_test additionally reports validation loss each epoch (no
        effect on training - the model never sees the test set during fitting)."""
        train = np.concatenate((s1_train, s2_train), axis=-1)
        validation_data = None
        if s1_test is not None and s2_test is not None:
            test = np.concatenate((s1_test, s2_test), axis=-1)
            validation_data = (test, None)  # None target - the loss is self-supervised (reconstruction)

        self.vae.fit(
            train,
            epochs=self.args.epochs,
            batch_size=self.args.bs,
            shuffle=True,
            validation_data=validation_data,
        )
        if self.args.save_model:
            self.vae.save_weights(self.args.model_out)

    def predict(self, s1_data, s2_data):
        """Encode (clinical, gene expression) into the latent embedding.
        Returns z_mean (index [0] of the encoder's outputs) - the deterministic
        center of the latent distribution, not a random sample - since that's
        what downstream classifiers/plots should use for reproducible results."""
        return self.encoder.predict(
            np.concatenate((s1_data, s2_data), axis=1), batch_size=self.args.bs
        )[0]
