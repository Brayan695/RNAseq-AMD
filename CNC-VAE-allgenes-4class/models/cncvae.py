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
    """

    def __init__(self, beta, distance_type, **kwargs):
        super().__init__(**kwargs)
        self.beta = beta
        self.distance_type = distance_type

    def call(self, inputs):
        x, outputs, z_mean, z_log_sigma, z = inputs
        if self.distance_type == 'mmd':
            true_samples = K.random_normal(tf.shape(z))
            distance = mmd(true_samples, z)
        else:
            distance = kl_regu(z_mean, z_log_sigma)

        reconstruction_loss = mean_squared_error(x, outputs)
        vae_loss = K.mean(reconstruction_loss + self.beta * distance)
        self.add_loss(vae_loss)
        return outputs


class CNCVAE:
    """CNC-VAE for Clin+mRNA integration (AMD cohort has no CNA arm, unlike METABRIC)."""

    def __init__(self, args):
        self.args = args
        self.vae = None
        self.encoder = None

    def build_model(self):
        np.random.seed(42)
        tf.random.set_seed(42)

        # ------------ Input -----------------
        inputs = Input(shape=(self.args.input_size,), name='concat_input')

        # ------------ Encoding Layer -----------------
        x = Dense(self.args.ds, activation=self.args.act)(inputs)
        x = BN()(x)

        # ------------ Embedding Layer --------------
        z_mean = Dense(self.args.ls, name='z_mean')(x)
        z_log_sigma = Dense(self.args.ls, name='z_log_sigma', kernel_initializer='zeros')(x)
        z = Lambda(sampling, output_shape=(self.args.ls,), name='z')([z_mean, z_log_sigma])

        self.encoder = Model(inputs, [z_mean, z_log_sigma, z], name='encoder')
        self.encoder.summary()

        # ------------ Decoding Layer -----------------
        latent_inputs = Input(shape=(self.args.ls,), name='z_sampling')
        x = Dense(self.args.ds, activation=self.args.act)(latent_inputs)
        x = BN()(x)
        x = Dropout(self.args.dropout)(x)

        # ------------ Out -----------------------
        concat_out = Dense(self.args.input_size)(x)

        decoder = Model(latent_inputs, concat_out, name='decoder')
        decoder.summary()

        outputs = decoder(self.encoder(inputs)[2])
        outputs = _VAELossLayer(self.args.beta, self.args.distance, name='vae_loss')(
            [inputs, outputs, z_mean, z_log_sigma, z]
        )
        self.vae = Model(inputs, outputs, name='vae_mlp')

        adam = optimizers.Adam(learning_rate=0.001, beta_1=0.9, beta_2=0.999)
        self.vae.compile(optimizer=adam)
        self.vae.summary()

    def train(self, s1_train, s2_train, s1_test=None, s2_test=None):
        train = np.concatenate((s1_train, s2_train), axis=-1)
        validation_data = None
        if s1_test is not None and s2_test is not None:
            test = np.concatenate((s1_test, s2_test), axis=-1)
            validation_data = (test, None)

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
        return self.encoder.predict(
            np.concatenate((s1_data, s2_data), axis=1), batch_size=self.args.bs
        )[0]
