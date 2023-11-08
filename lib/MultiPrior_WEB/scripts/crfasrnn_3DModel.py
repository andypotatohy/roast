# -*- coding: utf-8 -*-
"""
Created on Thu Aug 30 15:26:01 2018

@author: hirsch
"""


from keras.models import Model
from keras.layers.core import Dense, Activation, Dropout
from keras.layers.convolutional import Conv3D
from keras.initializers import he_normal
from keras.initializers import Orthogonal
from keras.layers.normalization import BatchNormalization
from keras.regularizers import l1_l2
from keras.layers import Input, Flatten, Reshape, Permute
from keras.layers.merge import Concatenate
from keras.layers import MaxPooling3D
from keras.layers import AveragePooling3D
from keras.layers.convolutional import Cropping3D
from keras.layers import UpSampling3D
from keras.layers import concatenate
from keras.layers.advanced_activations import PReLU
from keras.utils import print_summary
from keras import regularizers
from keras.optimizers import RMSprop
from keras.optimizers import Adam
import numpy as np
from keras.activations import softmax
import keras.backend as K
from keras.layers import Lambda
from keras.layers import Multiply

import sys
sys.path.append('/home/hirsch/Documents/projects/brainSegmentation/crfasrnn_keras-master/src/')          
from crfrnn_layer import CrfRnnLayer   



def dice_coef(y_true, y_pred):
    smooth = 1e-9
    y_true_f = K.flatten(y_true)
    y_pred_f = K.flatten(y_pred)
    intersection = K.sum(y_true_f * y_pred_f)
    return (2. * intersection + smooth) / (K.sum(y_true_f) + K.sum(y_pred_f) + smooth)

def Generalised_dice_coef_multilabel2(y_true, y_pred, numLabels=2):
    dice=0
    for index in range(numLabels):
        dice -= dice_coef(y_true[:,:,:,:,index], y_pred[:,:,:,:,index])
    return numLabels + dice

def dice_coef_multilabel6(y_true, y_pred, numLabels=6):
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
                                    
    #PENALTY = np.array([   [ 1, -1, -1, -1,  0,  0],
    #                       [-1,  1,  0,  0, -1, -1],
    #                       [-1,  0,  1,  0, -1, -1],
    #                       [-1,  0,  0,  1,  0, -1],
    #                       [ 0, -1, -1,  0,  1,  0],
    #                       [ 0, -1, -1, -1,  0,  1]], dtype='float32')

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


dpatch = 51
output_classes = 6
conv_features = [10,10,10,10,10,10,10,10]#[30,30,40,40,50,50,60,60,60,60]
fc_features = [50,100, output_classes]
d_factor = 3  # downsampling factor = stride in downsampling pathway
num_channels = 1
L2 = 0.001

learning_rate = 0.0002

#w_initializer=w_initializer, # initialization of layer parameters? Needed here?
#w_regularizer=w_regularizer,
#b_initializer=b_initializer, # initialization of bias parameters? Needed here?
#b_regularizer=b_regularizer,
#acti_func=acti_func


mod1      = Input((dpatch,dpatch,dpatch, num_channels))

#############   Normal pathway   ##################  

x1        = Cropping3D(cropping = ((13,13),(13,13),(13,13)), input_shape=(dpatch,dpatch,dpatch, num_channels))(mod1)
x1        = Conv3D(filters = conv_features[0], 
                   kernel_size = (3,3,3), 
                   #kernel_initializer=he_normal(seed=seed),
                   kernel_initializer=Orthogonal(),
                   kernel_regularizer=regularizers.l2(L2))(x1)
x1        = BatchNormalization()(x1)
#x1        = Activation('relu')(x1)
x1        = PReLU()(x1)
#x1        = BatchNormalization()(x1)



for feature in conv_features[1:]:  
    x1        = Conv3D(filters = feature, 
                       kernel_size = (3,3,3), 
                       #kernel_initializer=he_normal(seed=seed),
                       kernel_initializer=Orthogonal(),
                       kernel_regularizer=regularizers.l2(L2))(x1)
    x1        = BatchNormalization()(x1)
    #x1        = Activation('relu')(x1)
    x1        = PReLU()(x1)
    #x1        = BatchNormalization()(x1)
    
#############   Downsampled pathway   ##################   
#x2        = MaxPooling3D(pool_size=(d_factor,d_factor,d_factor), padding="same")(mod1)

x2        = AveragePooling3D(pool_size=(d_factor,d_factor,d_factor), padding="same")(mod1)

x3        = Conv3D(filters = 5, 
                   kernel_size = (7,7,7), 
                   padding="same",
                   strides=(3,3,3),
                   kernel_initializer=Orthogonal(),
                   )(mod1)
x3        = BatchNormalization()(x3)
x3        = PReLU()(x3)

x2 = concatenate([x2,x3])

x2        = Conv3D(filters = conv_features[0], 
                   kernel_size = (3,3,3), 
                   #kernel_initializer=he_normal(seed=seed),
                   #dilation_rate = (17,17,17),
                   kernel_initializer=Orthogonal(),
                   kernel_regularizer=regularizers.l2(L2))(x2)
x2        = BatchNormalization()(x2)
#x2        = Activation('relu')(x2)
x2        = PReLU()(x2)
#x2        = BatchNormalization()(x2)

for feature in (conv_features[1:]):    
    x2        = Conv3D(filters = feature, 
                       kernel_size = (3,3,3), 
                       #kernel_initializer=he_normal(seed=seed),
                       kernel_initializer=Orthogonal(),
                       kernel_regularizer=regularizers.l2(L2))(x2)
    x2        = BatchNormalization()(x2)
    #x2        = Activation('relu')(x2)
    x2        = PReLU()(x2)
    #x2        = BatchNormalization()(x2)

x2        = UpSampling3D(size=(9,9,9))(x2)

#############   Fully connected layers   ################## 

x        = concatenate([x1,x2])
#x        = Reshape(target_shape = (1, 60))(x)   # do I need this?
'''x        = Flatten()(x)
x        = Dense(units = fc_features[0], activation = 'elu')(x)
x        = Dropout(rate = 0.5)(x)
x        = Dense(units = fc_features[1], activation = 'elu')(x)
x        = Dropout(rate = 0.5)(x)    
x        = Dense(units = fc_features[2], activation = 'softmax', name = 'softmax')(x)'''

#   Fully convolutional variant

#x        = Dropout(rate = dropout[0])(x)
x        = Conv3D(filters = fc_features[0], 
                   kernel_size = (1,1,1), 
                   #kernel_initializer=he_normal(seed=seed),
                   kernel_initializer=Orthogonal(),
                   kernel_regularizer=regularizers.l2(L2))(x)
x        = BatchNormalization()(x)
#x        = Activation('relu')(x)
x        = PReLU()(x)
#x        = BatchNormalization()(x)


x        = Conv3D(filters = fc_features[1], 
                   kernel_size = (1,1,1), 
                   #kernel_initializer=he_normal(seed=seed),
                   kernel_initializer=Orthogonal(),
                   kernel_regularizer=regularizers.l2(L2))(x)
x        = BatchNormalization()(x)
#x        = Activation('relu')(x)
x        = PReLU()(x)
#x        = BatchNormalization()(x)

#x        = Flatten()(x)
      
x        = Conv3D(filters = output_classes, 
                   kernel_size = (1,1,1), 
                   #kernel_initializer=he_normal(seed=seed),
                   kernel_initializer=Orthogonal(),
                   kernel_regularizer=regularizers.l2(L2))(x)
#x        = BatchNormalization()(x)

# NO ACTIVATION (LINEAR ACTIVATION), THIS IS JUST THE LAST LAYER BEFORE SOFTMAX. hERE WE GET USUAL LOGITS.

# Skipping this bandfilter and going straigth to the softmax makes everything pointless (no nonlinearity besides softmax), and pushes performance to the floor.
x        = Conv3D(filters = fc_features[2], 
           kernel_size = (1,1,1), 
           #kernel_initializer=he_normal(seed=seed),
           kernel_initializer=Orthogonal(),
           kernel_regularizer=regularizers.l2(L2))(x)
x        = BatchNormalization()(x)
x        = PReLU()(x)
       
       


img_output = Cropping3D(cropping = (((dpatch-9)/2,(dpatch-9)/2),((dpatch-9)/2,(dpatch-9)/2),((dpatch-9)/2,(dpatch-9)/2)), input_shape=(dpatch,dpatch,dpatch, num_channels))(mod1)


F = Lambda(lambda x, i: x[:,:,i,:], output_shape=(9,9,), name="lambda_layer")  # Define your lambda layer
  # Update extra arguments to F


F.arguments = {'i': 0}
y0 = CrfRnnLayer(image_dims=(9,9,),
                 num_classes=6,
                 theta_alpha=160.,
                 theta_beta=3.,
                 theta_gamma=3.,
                 num_iterations=4,
                 name='crfrnn0')([F(x), F(img_output)])

F.arguments = {'i': 1}
y1 = CrfRnnLayer(image_dims=(9,9,),
                 num_classes=6,
                 theta_alpha=160.,
                 theta_beta=3.,
                 theta_gamma=3.,
                 num_iterations=4,
                 name='crfrnn1')([F(x), F(img_output)])

F.arguments = {'i': 2}          
y2 = CrfRnnLayer(image_dims=(9,9,),
                 num_classes=6,
                 theta_alpha=160.,
                 theta_beta=3.,
                 theta_gamma=3.,
                 num_iterations=4,
                 name='crfrnn2')([F(x), F(img_output)])

F.arguments = {'i': 3}          
y3 = CrfRnnLayer(image_dims=(9,9,),
                 num_classes=6,
                 theta_alpha=160.,
                 theta_beta=3.,
                 theta_gamma=3.,
                 num_iterations=4,
                 name='crfrnn3')([F(x), F(img_output)])
F.arguments = {'i': 4}          
y4 = CrfRnnLayer(image_dims=(9,9,),
                 num_classes=6,
                 theta_alpha=160.,
                 theta_beta=3.,
                 theta_gamma=3.,
                 num_iterations=4,
                 name='crfrnn4')([F(x), F(img_output)])
F.arguments = {'i': 5}          
y5 = CrfRnnLayer(image_dims=(9,9,),
                 num_classes=6,
                 theta_alpha=160.,
                 theta_beta=3.,
                 theta_gamma=3.,
                 num_iterations=4,
                 name='crfrnn5')([F(x), F(img_output)])
F.arguments = {'i': 6}          
y6 = CrfRnnLayer(image_dims=(9,9,),
                 num_classes=6,
                 theta_alpha=160.,
                 theta_beta=3.,
                 theta_gamma=3.,
                 num_iterations=4,
                 name='crfrnn6')([F(x), F(img_output)])
F.arguments = {'i': 7}          
y7 = CrfRnnLayer(image_dims=(9,9,),
                 num_classes=6,
                 theta_alpha=160.,
                 theta_beta=3.,
                 theta_gamma=3.,
                 num_iterations=4,
                 name='crfrnn7')([F(x), F(img_output)])
F.arguments = {'i': 8}          
y8 = CrfRnnLayer(image_dims=(9,9,),
                 num_classes=6,
                 theta_alpha=160.,
                 theta_beta=3.,
                 theta_gamma=3.,
                 num_iterations=4,
                 name='crfrnn8')([F(x), F(img_output)])

  
from keras.backend import stack


S = Lambda(lambda x: stack(x,axis=3), output_shape=(9,9,9,6,))
output = S([y0,y1,y2,y3,y4,y5,y6,y7,y8])
  
model     = Model(inputs=[mod1], outputs=output)
model.summary()
print_summary(model, positions=[.33, .6, .67,1])


model.compile(loss=dice_coef_multilabel6, optimizer=Adam(lr=learning_rate), metrics=[dice_coef_multilabel0,dice_coef_multilabel1,dice_coef_multilabel2,dice_coef_multilabel3, dice_coef_multilabel4,dice_coef_multilabel5])






X = np.zeros((5,51,51,51,1))

Y = model.predict(X)
X.shape
Y.shape


model.train_on_batch(X,Y)




X = np.zeros((1,51,51,51,1))

Y = model.predict(X)
X.shape
Y.shape


model.train_on_batch(X,Y)





valbatch.shape
vallabels.shape

x1 = valbatch[0:2]
x1.shape

y1 = vallabels[0:2]
y1.shape

model.train_on_batch(x1,y1)
    

    
Y2 = model.predict(x1)
Y2.shape