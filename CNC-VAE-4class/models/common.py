import tensorflow as tf
from tensorflow.keras import backend as K


def sse(true, pred):
    return K.sum(K.square(true - pred), axis=1)


def bce(true, pred):
    return K.sum(K.binary_crossentropy(true, pred), axis=1)


def compute_kernel(x, y):
    """Gaussian/RBF kernel with bandwidth 2*dim (standard InfoVAE convention,
    matching the CNC_VAE_AMD_V3 reference notebook)."""
    x_size = tf.shape(x)[0]
    y_size = tf.shape(y)[0]
    dim = tf.shape(x)[1]
    tiled_x = tf.tile(tf.reshape(x, [x_size, 1, dim]), [1, y_size, 1])
    tiled_y = tf.tile(tf.reshape(y, [1, y_size, dim]), [x_size, 1, 1])
    bandwidth = 2.0 * tf.cast(dim, tf.float32)
    return tf.exp(-tf.reduce_mean(tf.square(tiled_x - tiled_y), axis=2) / bandwidth)


def mmd(x, y):
    x_kernel = compute_kernel(x, x)
    y_kernel = compute_kernel(y, y)
    xy_kernel = compute_kernel(x, y)
    return tf.reduce_mean(x_kernel) + tf.reduce_mean(y_kernel) - 2 * tf.reduce_mean(xy_kernel)


def sampling(args):
    """Reparameterization trick: sample from an isotropic unit Gaussian.

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
    kl_loss = 1 + z_log_sigma - K.square(z_mean) - K.exp(z_log_sigma)
    kl_loss = K.sum(kl_loss, axis=-1)
    kl_loss *= -0.5
    return kl_loss
