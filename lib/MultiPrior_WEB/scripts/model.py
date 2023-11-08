# -*- coding: utf-8 -*-
"""
Created on Thu Sep 28 15:21:22 2023

@author: Andrew
"""

import keras
import numpy as np
import keras.backend as K

def dice_coef(y_true, y_pred):
    smooth = 1e-9
    y_true_f = K.flatten(y_true)
    y_pred_f = K.flatten(y_pred)
    intersection = K.sum(y_true_f * y_pred_f)
    return (2. * intersection + smooth) / (K.sum(y_true_f) + K.sum(y_pred_f) + smooth)


def Generalised_dice_coef_multilabel7(y_true, y_pred, numLabels=7):
    dice=0
    for index in range(numLabels):
        dice -= dice_coef(y_true[:,:,:,:,index], y_pred[:,:,:,:,index])
    return numLabels + dice

def dice_coef_multilabel0(y_true, y_pred):
    index = 0
    dice = dice_coef(y_true[:,:,:,:,index], y_pred[:,:,:,:,index])
    return dice
def dice_coef_multilabel1(y_true, y_pred):
    index = 1
    dice = dice_coef(y_true[:,:,:,:,index], y_pred[:,:,:,:,index])
    return dice
def dice_coef_multilabel2(y_true, y_pred):
    index = 2
    dice = dice_coef(y_true[:,:,:,:,index], y_pred[:,:,:,:,index])
    return dice
def dice_coef_multilabel3(y_true, y_pred):
    index = 3
    dice = dice_coef(y_true[:,:,:,:,index], y_pred[:,:,:,:,index])
    return dice
def dice_coef_multilabel4(y_true, y_pred):
    index = 4
    dice = dice_coef(y_true[:,:,:,:,index], y_pred[:,:,:,:,index])
    return dice
def dice_coef_multilabel5(y_true, y_pred):
    index = 5
    dice = dice_coef(y_true[:,:,:,:,index], y_pred[:,:,:,:,index])
    return dice
def dice_coef_multilabel6(y_true, y_pred):
    index = 6
    dice = dice_coef(y_true[:,:,:,:,index], y_pred[:,:,:,:,index])
    return dice

my_custom_objects = {'Generalised_dice_coef_multilabel7':Generalised_dice_coef_multilabel7,
                     'dice_coef_multilabel0':dice_coef_multilabel0,
                     'dice_coef_multilabel1':dice_coef_multilabel1,
                     'dice_coef_multilabel2':dice_coef_multilabel2,
                     'dice_coef_multilabel3':dice_coef_multilabel3,
                     'dice_coef_multilabel4':dice_coef_multilabel4,
                     'dice_coef_multilabel5':dice_coef_multilabel5,
                     'dice_coef_multilabel6':dice_coef_multilabel6}

model_path = 'C:/Users/Andrew/Documents/MultiPriors_WEB/models/best_model_Tensorflow2.h5'.replace('/','\\')

model = keras.models.load_model(model_path, custom_objects=my_custom_objects)


