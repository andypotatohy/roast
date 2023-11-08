#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Dec  8 10:26:01 2017

@author: lukas
"""
from keras.models import Model
from keras.layers.core import Dense, Activation, Dropout
from keras.layers.convolutional import Conv3D
from keras.initializers import he_normal, Orthogonal
from keras.layers.normalization import BatchNormalization
from keras.regularizers import l1_l2
from keras.layers import Input, Flatten, Reshape, Permute
from keras.layers.merge import Concatenate
from keras.layers import MaxPooling3D
from keras.layers import AveragePooling3D
from keras.layers.convolutional import Cropping3D
from keras.layers import UpSampling3D
from keras.layers import concatenate
from keras.layers.advanced_activations import LeakyReLU
from keras.utils import print_summary
from keras import regularizers
from keras.optimizers import RMSprop
import keras.backend as K
from keras.optimizers import Adam
import numpy as np
from keras.activations import softmax


def dice_coef(y_true, y_pred):
    smooth = 1e-9
    y_true_f = K.flatten(y_true)
    y_pred_f = K.flatten(y_pred)
    intersection = K.sum(y_true_f * y_pred_f)
    return (2. * intersection + smooth) / (K.sum(y_true_f**2) + K.sum(y_pred_f**2) + smooth)

def Generalised_dice_coef_multilabel2(y_true, y_pred, numLabels=2):
    dice=0
    for index in range(numLabels):
        dice -= dice_coef(y_true[:,:,:,:,index], y_pred[:,:,:,:,index])
    return numLabels + dice

def Generalised_dice_coef_multilabel6(y_true, y_pred, numLabels=6):
    dice=0
    for index in range(numLabels):
        dice -= dice_coef(y_true[:,:,:,:,index], y_pred[:,:,:,:,index])
    return numLabels + dice
    
def Generalised_dice_coef_multilabel7(y_true, y_pred, numLabels=7):
    dice=0
    for index in range(numLabels):
        dice -= dice_coef(y_true[:,:,:,:,index], y_pred[:,:,:,:,index])
    return numLabels + dice

def w_dice_coef(y_true, y_pred, PENALTY):
    smooth = 1e-9
    y_true_f = K.flatten(y_true)
    y_pred_f = K.flatten(y_pred)
    intersection = K.sum(y_true_f * y_pred_f) * PENALTY
    return (2. * intersection + smooth) / (K.sum(y_true_f) + K.sum(y_pred_f) + smooth)

def w_dice_coef_multilabel2(y_true, y_pred, numLabels=2):
                                    
    PENALTY = np.array([[ 1,    0],
                        [-1,  1]], dtype='float32')
                          
    #PENALTY[PENALTY < 0] = -5                      
                             
    dice = []
    for index in range(numLabels):
        dice_class = []
        for j in range(numLabels):
            wDice = w_dice_coef(y_true[:,:,:,:,index], y_pred[:,:,:,:,j], PENALTY[index,j])
            dice_class.append(wDice)
        dice.append(K.sum(dice_class)) 
        
    final_dice = K.sum(dice)
    return numLabels - final_dice

def w_dice_coef_multilabel6(y_true, y_pred, numLabels=6):

    PENALTY = np.array([   [ 1, -1, -1, -1,  0,  0],
                           [-1,  1,  0,  0, -1, -1],
                           [-1,  0,  1,  0, -1, -1],
                           [-1,  0,  0,  1, -1, -1],
                           [ 0, -1, -1, -1,  1,  0],
                           [ 0, -1, -1, -1,  0,  1]], dtype='float32')
                          
    PENALTY[PENALTY < 0] = -1                   

    #PENALTY = np.eye(6,6, dtype='float32')
                             
    dice = []
    for index in range(numLabels):
        dice_class = []
        for j in range(numLabels):
            wDice = w_dice_coef(y_true[:,:,:,:,index], y_pred[:,:,:,:,j], PENALTY[index,j])
            dice_class.append(wDice)
        dice.append(K.sum(dice_class)) 
        
    final_dice = K.sum(dice)
    return numLabels - final_dice

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

    

class MultiPriors_Model():
    
    def __init__(self, output_classes, num_channels, L2, dropout, learning_rate, optimizer_decay, loss_function):
        
        self.output_classes = output_classes
        self.conv_features = [30, 30, 30, 30, 50, 50, 50, 50]   
        self.fc_features = [100,150]
        self.num_channels = num_channels
        self.L2 = L2
        self.dropout = dropout
        self.learning_rate = learning_rate
        self.optimizer_decay = optimizer_decay
        self.loss_function = loss_function

    def createModel(self):
        '''Creates model architecture
        Input: Data input dimensions, eventually architecture specifications parsed from a config file? (activations, costFunction, hyperparameters (nr layers), dropout....)
        Output: Keras Model'''
    
        highRes      = Input((None, None, None, self.num_channels))    # Dim = 25^3 (from a 57^3 cube, cropped 16 per side)
        lowRes       = Input((None, None, None, self.num_channels))    # Dim = 19^3 (from a 57^3 cube downsampled by three)
        #############   Normal pathway   ##################  

        # This input is now already cropped     = -32 
        x1 = highRes        
        
        # 25 --> 9  =  -16
        for kernels in self.conv_features:   
            x1        = Conv3D(filters = kernels, 
                               kernel_size = (3,3,3), 
                               kernel_initializer=Orthogonal(),
                               kernel_regularizer=regularizers.l2(self.L2))(x1)
            x1        = LeakyReLU()(x1)
            x1        = BatchNormalization()(x1)   
            
        #############   Downsampled pathway   ##################   

        # This input is already downsampled       = /3  
        x2 = lowRes        

        # in total (x/3 - 16)*3  =  -66
        # -22
        for kernels in self.conv_features:   
            x2        = Conv3D(filters = kernels, 
                               kernel_size = (3,3,3), 
                               kernel_initializer=Orthogonal(),
                               kernel_regularizer=regularizers.l2(self.L2))(x2)
            x2        = LeakyReLU()(x2)
            x2        = BatchNormalization()(x2)   
        
        #x2        = UpSampling3D(size=(3,3,3))(x2)
        
        #############   Fully connected layers   ################## 

        tpm = Input((None,None,None,6))
        
        x        = concatenate([x1, x2, tpm])

        for feature in self.fc_features:    
            x        = Conv3D(filters = feature, 
                               kernel_size = (1,1,1), 
                               kernel_initializer=Orthogonal(),
                               kernel_regularizer=regularizers.l2(self.L2))(x)
            x        = LeakyReLU()(x)
            x        = BatchNormalization()(x)   

        x        = Conv3D(filters = self.output_classes, 
                   kernel_size = (1,1,1), 
                   kernel_initializer=Orthogonal(),
                   kernel_regularizer=regularizers.l2(self.L2))(x)
        x        = Activation(softmax)(x)
        
        model     = Model(inputs=[highRes, lowRes, tpm], outputs=x)       
        
        if self.loss_function == 'Multinomial':
            model.compile(loss='categorical_crossentropy', optimizer=Adam(lr=self.learning_rate), metrics=[dice_coef_multilabel0,dice_coef_multilabel1,dice_coef_multilabel2,dice_coef_multilabel3, dice_coef_multilabel4,dice_coef_multilabel5,dice_coef_multilabel6])
        elif self.loss_function == 'Dice2':
            model.compile(loss=Generalised_dice_coef_multilabel2, optimizer=Adam(lr=self.learning_rate), metrics=[dice_coef_multilabel0,dice_coef_multilabel1])
        elif self.loss_function == 'Dice6':
            model.compile(loss=dice_coef_multilabel6, optimizer=Adam(lr=self.learning_rate), metrics=[dice_coef_multilabel0,dice_coef_multilabel1,dice_coef_multilabel2,dice_coef_multilabel3, dice_coef_multilabel4,dice_coef_multilabel5])
        elif self.loss_function == 'wDice6':
            model.compile(loss=w_dice_coef_multilabel6, optimizer=Adam(lr=self.learning_rate), metrics=[dice_coef_multilabel0,dice_coef_multilabel1,dice_coef_multilabel2,dice_coef_multilabel3, dice_coef_multilabel4,dice_coef_multilabel5]) 
        elif self.loss_function == 'wDice2':
            model.compile(loss=w_dice_coef_multilabel2, optimizer=Adam(lr=self.learning_rate), metrics=[dice_coef_multilabel0,dice_coef_multilabel1]) 
        elif self.loss_function == 'Dice7':
            model.compile(loss=Generalised_dice_coef_multilabel7, optimizer=Adam(lr=self.learning_rate), metrics=[dice_coef_multilabel0,dice_coef_multilabel1,dice_coef_multilabel2,dice_coef_multilabel3, dice_coef_multilabel4,dice_coef_multilabel5,dice_coef_multilabel6 ])
        return model


#dm = MultiPriors_Model(7, 1, 0.001, [0], 0.01, 0, 'Dice7' )
#model = dm.createModel()            
#model.summary()  
#from keras.utils import plot_model
#plot_model(model, to_file='/home/hirsch/Documents/projects/brainSegmentation/DeepPriors/multiscale_shapes_noUpsampling.png', show_shapes=True)    
#
#dpatch = 42*3
#
#X = np.random.randn(1,dpatch,dpatch,dpatch,1)
#TPM = np.random.randn(1,dpatch-48,dpatch-48,dpatch-48,6)
#
#from skimage.transform import resize
#highRes = X[:,16:-16,16:-16,16:-16,:]
#highRes.shape
#lowRes = resize(X, output_shape=highRes.shape, anti_aliasing=True)
#lowRes.shape
#
#yhat = model.predict([highRes,lowRes, TPM])
#yhat.shape

