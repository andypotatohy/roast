"""Ported from lib/multiaxial/utils.py.

Only the loss/metric functions actually referenced by the bundled models'
Keras ``custom_objects`` at load time are ported here (the training-only
utilities in utils.py -- data generators, model-building functions,
callbacks -- are out of scope: roast_py only needs to run inference with
the already-trained weights, not retrain them).
"""

from __future__ import annotations

import tensorflow as tf


def dice_coef(y_true, y_pred):
    smooth = 1e-6
    y_true_f = tf.keras.backend.flatten(y_true)
    y_pred_f = tf.keras.backend.flatten(y_pred)
    intersection = tf.reduce_sum(y_true_f * y_pred_f)
    return (2.0 * intersection + smooth) / (
        tf.reduce_sum(y_true_f**2) + tf.reduce_sum(y_pred_f**2) + smooth
    )


def Generalised_dice_coef_multilabel7(y_true, y_pred, numLabels=7):
    dice = 0
    for index in range(numLabels):
        dice -= dice_coef(y_true[:, :, :, index], y_pred[:, :, :, index])
    return numLabels + dice


def _dice_coef_multilabel_bin(index):
    def f(y_true, y_pred):
        numerator = 2 * tf.math.reduce_sum(y_true[:, :, :, index] * y_pred[:, :, :, index])
        denominator = tf.math.reduce_sum(y_true[:, :, :, index] + y_pred[:, :, :, index])
        return numerator / denominator

    return f


dice_coef_multilabel_bin0 = _dice_coef_multilabel_bin(0)
dice_coef_multilabel_bin1 = _dice_coef_multilabel_bin(1)
dice_coef_multilabel_bin2 = _dice_coef_multilabel_bin(2)
dice_coef_multilabel_bin3 = _dice_coef_multilabel_bin(3)
dice_coef_multilabel_bin4 = _dice_coef_multilabel_bin(4)
dice_coef_multilabel_bin5 = _dice_coef_multilabel_bin(5)
dice_coef_multilabel_bin6 = _dice_coef_multilabel_bin(6)

CUSTOM_OBJECTS = {
    "Generalised_dice_coef_multilabel7": Generalised_dice_coef_multilabel7,
    "dice_coef_multilabel_bin0": dice_coef_multilabel_bin0,
    "dice_coef_multilabel_bin1": dice_coef_multilabel_bin1,
    "dice_coef_multilabel_bin2": dice_coef_multilabel_bin2,
    "dice_coef_multilabel_bin3": dice_coef_multilabel_bin3,
    "dice_coef_multilabel_bin4": dice_coef_multilabel_bin4,
    "dice_coef_multilabel_bin5": dice_coef_multilabel_bin5,
    "dice_coef_multilabel_bin6": dice_coef_multilabel_bin6,
}
