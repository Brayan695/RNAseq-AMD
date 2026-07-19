# Small math building blocks shared by the tf.keras VAE model: the two
# alternative "keep the latent space close to a standard normal prior"
# regularizers (MMD or KL divergence - pick one via --distance), plus the
# reparameterization trick used to sample a latent vector z in a way that's
# still differentiable. (See CNC-VAE/models/common.py for the PyTorch
# equivalent with the same math, if comparing the two implementations.)
import keras
import tensorflow as tf
from tensorflow.keras import backend as K


def sse(true, pred):
    """Sum of squared errors per sample. Not currently used by CNCVAE (which
    uses plain MSE instead) - kept for parity with the original repo."""
    return K.sum(K.square(true - pred), axis=1)


def bce(true, pred):
    """Binary cross-entropy per sample, summed over features. Not currently
    used by CNCVAE - kept for parity with the original repo."""
    return K.sum(K.binary_crossentropy(true, pred), axis=1)


def compute_kernel(x, y):
    """Gaussian/RBF kernel with bandwidth 2*dim (standard InfoVAE convention,
    matching the CNC_VAE_AMD_V3 reference notebook).

    Pairwise similarity between every row of x and every row of y:
    compute_kernel(x,y)[i,j] = exp(-||x_i - y_j||^2 / bandwidth). Used as a
    building block for MMD below.
    """
    x_size = tf.shape(x)[0]
    y_size = tf.shape(y)[0]
    dim = tf.shape(x)[1]
    tiled_x = tf.tile(tf.reshape(x, [x_size, 1, dim]), [1, y_size, 1])
    tiled_y = tf.tile(tf.reshape(y, [1, y_size, dim]), [x_size, 1, 1])
    bandwidth = 2.0 * tf.cast(dim, tf.float32)
    return tf.exp(-tf.reduce_mean(tf.square(tiled_x - tiled_y), axis=2) / bandwidth)


def mmd(x, y):
    """Maximum Mean Discrepancy between samples x and y: near 0 when x and y
    look like they came from the same distribution, grows as they differ.
    MMD^2 estimator = mean within-x similarity + mean within-y similarity
    - 2 * mean cross similarity."""
    x_kernel = compute_kernel(x, x)
    y_kernel = compute_kernel(y, y)
    xy_kernel = compute_kernel(x, y)
    return tf.reduce_mean(x_kernel) + tf.reduce_mean(y_kernel) - 2 * tf.reduce_mean(xy_kernel)


# Keras 3 needs custom functions used inside a Lambda layer to be explicitly
# registered, so a saved model (.keras file) can be reloaded later - without
# this decorator, loading a saved encoder would fail with "Could not locate
# function 'sampling'". Importing this module (e.g. `import models.common`)
# is what actually runs this registration; a script that only calls
# keras.models.load_model() without ever importing models.common will still fail.
@keras.saving.register_keras_serializable(package='cncvae')
def sampling(args):
    """Reparameterization trick: sample from an isotropic unit Gaussian.

    Draws z ~ N(z_mean, exp(z_log_var)) in a way that's still differentiable
    w.r.t. z_mean/z_log_var (by sampling independent noise `epsilon` and
    combining it via simple arithmetic, rather than sampling z directly).

    Arguments:
        args (tensor): mean and log-variance of Q(z|X)
    Returns:
        z (tensor): sampled latent vector
    """
    z_mean, z_log_var = args
    batch = tf.shape(z_mean)[0]
    dim = tf.shape(z_mean)[1]
    epsilon = tf.random.normal(shape=(batch, dim))
    return z_mean + tf.exp(0.5 * z_log_var) * epsilon


def kl_regu(z_mean, z_log_sigma):
    """Per-sample KL( N(z_mean, exp(z_log_sigma)) || N(0, I) ) - the standard
    VAE regularization term, closed-form since both distributions are Gaussian.
    Penalizes the encoder for producing latent distributions that stray far
    from a standard normal prior."""
    kl_loss = 1 + z_log_sigma - K.square(z_mean) - K.exp(z_log_sigma)
    kl_loss = K.sum(kl_loss, axis=-1)  # sum over latent dims -> one value per sample
    kl_loss *= -0.5
    return kl_loss
