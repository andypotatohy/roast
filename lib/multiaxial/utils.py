#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun  5 12:20:18 2024

@author: deeperthought
"""


import numpy as np
import tensorflow as tf
from tensorflow.keras import backend as K


from tensorflow.keras import Input, Model
from tensorflow.keras.layers import Conv2D, MaxPooling2D, UpSampling2D, Activation, BatchNormalization, Conv2DTranspose, concatenate
from tensorflow.keras.optimizers import Adam
from tensorflow.keras import regularizers




def dice_loss(y_true, y_pred):
  numerator = 2 * tf.math.reduce_sum(y_true * y_pred)
  denominator = tf.math.reduce_sum(y_true + y_pred)
  return 1 - numerator / denominator


def dice_coef(y_true, y_pred):
    smooth = 1e-6
    y_true_f = K.flatten(y_true)
    y_pred_f = K.flatten(y_pred)
    intersection = tf.reduce_sum(y_true_f * y_pred_f)
    return (2. * intersection + smooth) / (tf.reduce_sum(y_true_f**2) + tf.reduce_sum(y_pred_f**2) + smooth)


def Generalised_dice_coef_multilabel7(y_true, y_pred, numLabels=7):
    """This is the loss function to MINIMIZE. A perfect overlap returns 0. Total disagreement returns numeLabels"""
    dice=0
    for index in range(numLabels):
        dice -= dice_coef(y_true[:,:,:,index], y_pred[:,:,:,index])
        
    return numLabels + dice

@tf.function
def generalized_dice_coef_multilabel7_present_class(y_true, y_pred):
    """This is the loss function to MINIMIZE. A perfect overlap returns 0. Total disagreement returns numLabels"""
    y_true = tf.cast(y_true, tf.float32)
    y_pred = tf.cast(y_pred, tf.float32)
    "If the class is absent and the prediction of the class is absent, it should return 1. Perfect agreement!"
    dice = 7
    empty = tf.constant(1, tf.float32)
    for index in range(7):
        flag1 = tf.math.reduce_sum(y_true[:,:,:,index])
        #flag2 = tf.math.reduce_sum(y_pred[:,:,:,index])
        if flag1 == 0 :#and flag2 == 0:
            dice -= empty
        else:
            dice -= dice_coef(y_true[:,:,:,index], y_pred[:,:,:,index])
               
    return dice

# def Generalised_dice_coef_neighbors(y_true, y_pred, numLabels=7):
#     """This is the loss function to MINIMIZE. A perfect overlap returns 0. Total disagreement returns numeLabels"""
#     dice=0
#     for index in range(numLabels):
#         dice -= dice_coef(y_true[:,:,:,index], y_pred[:,:,:,index])
        
#       # Get the dimensions of the predicted volume
#     depth, height, width = y_pred.shape[1:]
    
#     brain_classes = [2,3]
#     nonbrain_classes = [0,6]
    
#     # Create a mask for neighboring voxels
#     neighbor_mask = tf.zeros_like(y_pred)
#     for i in range(1, depth - 1):
#       for j in range(1, height - 1):
#         for k in range(1, width - 1):
#           if y_pred[0, i, j, k] in brain_classes and y_pred[0, i, j - 1, k] in nonbrain_classes:
#             neighbor_mask[0, i, j, k] = 1
#           elif y_pred[0, i, j, k] in brain_classes and y_pred[0, i, j + 1, k] in nonbrain_classes:
#             neighbor_mask[0, i, j, k] = 1
#     #       # Add similar checks for other neighboring voxels (e.g., up, down, diagonal)
    
#     # Calculate the penalty based on the number of voxels in the mask
#     dice += tf.reduce_sum(neighbor_mask)    
    
#     return numLabels + dice

def dice_coef_multilabel_bin0(y_true, y_pred):
  numerator = 2 * tf.math.reduce_sum(y_true[:,:,:,0] * y_pred[:,:,:,0])
  denominator = tf.math.reduce_sum(y_true[:,:,:,0] + y_pred[:,:,:,0])
  return numerator / denominator

def dice_coef_multilabel_bin1(y_true, y_pred):
  numerator = 2 * tf.math.reduce_sum(y_true[:,:,:,1] * y_pred[:,:,:,1])
  denominator = tf.math.reduce_sum(y_true[:,:,:,1] + y_pred[:,:,:,1])
  return numerator / denominator

def dice_coef_multilabel_bin2(y_true, y_pred):
    
  numerator = 2 * tf.math.reduce_sum(y_true[:,:,:,2] * y_pred[:,:,:,2])
  denominator = tf.math.reduce_sum(y_true[:,:,:,2] + y_pred[:,:,:,2])
  return numerator / denominator

def dice_coef_multilabel_bin3(y_true, y_pred):
  numerator = 2 * tf.math.reduce_sum(y_true[:,:,:,3] * y_pred[:,:,:,3])
  denominator = tf.math.reduce_sum(y_true[:,:,:,3] + y_pred[:,:,:,3])
  return numerator / denominator

def dice_coef_multilabel_bin4(y_true, y_pred):
  numerator = 2 * tf.math.reduce_sum(y_true[:,:,:,4] * y_pred[:,:,:,4])
  denominator = tf.math.reduce_sum(y_true[:,:,:,4] + y_pred[:,:,:,4])
  return numerator / denominator

def dice_coef_multilabel_bin5(y_true, y_pred):
  numerator = 2 * tf.math.reduce_sum(y_true[:,:,:,5] * y_pred[:,:,:,5])
  denominator = tf.math.reduce_sum(y_true[:,:,:,5] + y_pred[:,:,:,5])
  return numerator / denominator

def dice_coef_multilabel_bin6(y_true, y_pred):
  numerator = 2 * tf.math.reduce_sum(y_true[:,:,:,6] * y_pred[:,:,:,6])
  denominator = tf.math.reduce_sum(y_true[:,:,:,6] + y_pred[:,:,:,6])
  return numerator / denominator



def myDiceMetric(tissueType):
    def dice_metric(y_true, y_pred):
        # I want a discrete dice score. So taking the argmax and then retriving score per class
        y_pred = tf.math.argmax(y_pred, -1) # (batch, x,y)
        y_true = tf.math.argmax(y_true, -1)
        
        y_pred_mask = tf.cast(y_pred == tissueType, tf.float32)
        y_true_mask = tf.cast(y_true == tissueType, tf.float32)
        
        numerator = 2 * tf.math.reduce_sum(y_true_mask * y_pred_mask)
        denominator = tf.math.reduce_sum(y_true_mask + y_pred_mask)
        return (numerator + tf.keras.backend.epsilon())/ (denominator + tf.keras.backend.epsilon())
    return dice_metric


def dice_metric0(y_true, y_pred):
    # I want a discrete dice score. So taking the argmax and then retriving score per class
    y_pred = tf.math.argmax(y_pred, -1) # (batch, x,y)
    y_true = tf.math.argmax(y_true, -1)
    
    y_pred_mask = tf.cast(y_pred == 0, tf.float32)
    y_true_mask = tf.cast(y_true == 0, tf.float32)
    
    numerator = 2 * tf.math.reduce_sum(y_true_mask * y_pred_mask)
    denominator = tf.math.reduce_sum(y_true_mask + y_pred_mask)
    return (numerator + tf.keras.backend.epsilon())/ (denominator + tf.keras.backend.epsilon())


def dice_metric1(y_true, y_pred):
    # I want a discrete dice score. So taking the argmax and then retriving score per class
    y_pred = tf.math.argmax(y_pred, -1) # (batch, x,y)
    y_true = tf.math.argmax(y_true, -1)
    
    y_pred_mask = tf.cast(y_pred == 1, tf.float32)
    y_true_mask = tf.cast(y_true == 1, tf.float32)
    
    numerator = 2 * tf.math.reduce_sum(y_true_mask * y_pred_mask)
    denominator = tf.math.reduce_sum(y_true_mask + y_pred_mask)
    return (numerator + tf.keras.backend.epsilon())/ (denominator + tf.keras.backend.epsilon())


def dice_metric2(y_true, y_pred):
    # I want a discrete dice score. So taking the argmax and then retriving score per class
    y_pred = tf.math.argmax(y_pred, -1) # (batch, x,y)
    y_true = tf.math.argmax(y_true, -1)
    
    y_pred_mask = tf.cast(y_pred == 2, tf.float32)
    y_true_mask = tf.cast(y_true == 2, tf.float32)
    
    numerator = 2 * tf.math.reduce_sum(y_true_mask * y_pred_mask)
    denominator = tf.math.reduce_sum(y_true_mask + y_pred_mask)
    return (numerator + tf.keras.backend.epsilon())/ (denominator + tf.keras.backend.epsilon())


def dice_metric3(y_true, y_pred):
    # I want a discrete dice score. So taking the argmax and then retriving score per class
    y_pred = tf.math.argmax(y_pred, -1) # (batch, x,y)
    y_true = tf.math.argmax(y_true, -1)
    
    y_pred_mask = tf.cast(y_pred == 3, tf.float32)
    y_true_mask = tf.cast(y_true == 3, tf.float32)
    
    numerator = 2 * tf.math.reduce_sum(y_true_mask * y_pred_mask)
    denominator = tf.math.reduce_sum(y_true_mask + y_pred_mask)
    return (numerator + tf.keras.backend.epsilon())/ (denominator + tf.keras.backend.epsilon())



def dice_metric4(y_true, y_pred):
    # I want a discrete dice score. So taking the argmax and then retriving score per class
    y_pred = tf.math.argmax(y_pred, -1) # (batch, x,y)
    y_true = tf.math.argmax(y_true, -1)
    
    y_pred_mask = tf.cast(y_pred == 4, tf.float32)
    y_true_mask = tf.cast(y_true == 4, tf.float32)
    
    numerator = 2 * tf.math.reduce_sum(y_true_mask * y_pred_mask)
    denominator = tf.math.reduce_sum(y_true_mask + y_pred_mask)
    return (numerator + tf.keras.backend.epsilon())/ (denominator + tf.keras.backend.epsilon())



def dice_metric5(y_true, y_pred):
    # I want a discrete dice score. So taking the argmax and then retriving score per class
    y_pred = tf.math.argmax(y_pred, -1) # (batch, x,y)
    y_true = tf.math.argmax(y_true, -1)
    
    y_pred_mask = tf.cast(y_pred == 5, tf.float32)
    y_true_mask = tf.cast(y_true == 5, tf.float32)
    
    numerator = 2 * tf.math.reduce_sum(y_true_mask * y_pred_mask)
    denominator = tf.math.reduce_sum(y_true_mask + y_pred_mask)
    return (numerator + tf.keras.backend.epsilon())/ (denominator + tf.keras.backend.epsilon())



def dice_metric6(y_true, y_pred):
    # I want a discrete dice score. So taking the argmax and then retriving score per class
    y_pred = tf.math.argmax(y_pred, -1) # (batch, x,y)
    y_true = tf.math.argmax(y_true, -1)
    
    y_pred_mask = tf.cast(y_pred == 6, tf.float32)
    y_true_mask = tf.cast(y_true == 6, tf.float32)
    
    numerator = 2 * tf.math.reduce_sum(y_true_mask * y_pred_mask)
    denominator = tf.math.reduce_sum(y_true_mask + y_pred_mask)
    return (numerator + tf.keras.backend.epsilon())/ (denominator + tf.keras.backend.epsilon())

#%%



def double_conv_block(x, n_filters):
   # Conv2D then ReLU activation
   x = tf.keras.layers.Conv2D(n_filters, 3, padding = "same", activation = "relu", kernel_initializer = "he_normal")(x)
   # Conv2D then ReLU activation
   x = tf.keras.layers.Conv2D(n_filters, 3, padding = "same", activation = "relu", kernel_initializer = "he_normal")(x)
   return x

def downsample_block(x, n_filters):
   f = double_conv_block(x, n_filters)
   p = tf.keras.layers.MaxPool2D(2)(f)
   p = tf.keras.layers.Dropout(0)(p)
   return f, p

def upsample_block(x, conv_features, n_filters):
   # upsample
   x = tf.keras.layers.Conv2DTranspose(n_filters, 3, 2, padding="same")(x)
   # concatenate
   x = tf.keras.layers.concatenate([x, conv_features])
   # dropout
   x = tf.keras.layers.Dropout(0)(x)
   # Conv2D twice with ReLU activation
   x = double_conv_block(x, n_filters)
   return x



def create_convolution_block(input_layer, n_filters, kernel=(3, 3), padding='same', strides=(1, 1), L2=0, use_batch_norm=True):

    layer = Conv2D(n_filters, kernel, padding=padding, strides=strides, kernel_regularizer=regularizers.l2(L2))(input_layer)
    if use_batch_norm:
        layer = BatchNormalization()(layer)

    return Activation('relu')(layer)


def get_up_convolution(n_filters, pool_size=(2,2), kernel_size=(2,2), strides=(2, 2),
                       deconvolution=True, bilinear_upsampling=False, L2=0):
    if deconvolution:
        if bilinear_upsampling:
            return Conv2DTranspose(filters=n_filters, kernel_size=(3,3),
                                   strides=strides, trainable=False)#, kernel_initializer=make_bilinear_filter_5D(shape=(3,3,3,n_filters,n_filters)), trainable=False)
        else:
            return Conv2DTranspose(filters=n_filters, kernel_size=kernel_size,
                                   strides=strides, kernel_regularizer=regularizers.l2(L2))            
    else:
        return UpSampling2D(size=pool_size)

def my_init(shape, dtype=None):
    return K.random_normal(shape, dtype=dtype)


def slice_number_embedding(slice_number, embedding_dim):
  # One-hot encode slice number (optional)
  # slice_one_hot = tf.one_hot(slice_number, depth=max_slices)
  embedding = tf.keras.layers.Embedding(256, embedding_dim)(slice_number) 
  reshaped_embedding = tf.keras.layers.Reshape((1,1, embedding_dim))(embedding)

  return reshaped_embedding


def UNet_v0_2DTumorSegmenter(input_shape =  (256, 256,1), pool_size=(2, 2),initial_learning_rate=1e-5, deconvolution=True,
                      depth=6, n_base_filters=32, activation_name="softmax", L2=0, use_batch_norm=False, add_spatial_prior=False, 
                      embedding=False, embedding_dimensions=8 ):
        """ Simple version, padding 'same' on every layer, output size is equal to input size. Has border artifacts and checkerboard artifacts """
        inputs = Input(input_shape)
        
        if add_spatial_prior:
            positional_encoding_input = Input((256,256,4))
            # if embedding:
                # positional_encoding = slice_number_embedding(positional_encoding_input, embedding_dim=embedding_dimensions)
            # else:
                # positional_encoding = tf.keras.layers.Lambda(lambda x: tf.divide(x, 256.))(positional_encodin_input)

                # positional_encoding = tf.keras.layers.Reshape((1,256,256))(positional_encoding_input)
                
        
        levels = list()
        #current_layer = Conv2D(n_base_filters, (1, 1))(inputs)  # ???? needed??  Not even a nonlinear activation!!!
        current_layer = inputs
        positional_encoding = positional_encoding_input
    
        # add levels with max pooling
        for layer_depth in range(depth):
            layer1 = create_convolution_block(input_layer=current_layer, kernel=(3,3), n_filters=n_base_filters*(layer_depth+1), padding='same', L2=L2, use_batch_norm=use_batch_norm)
            layer2 = create_convolution_block(input_layer=layer1, kernel=(3,3),  n_filters=n_base_filters*(layer_depth+1), padding='same', L2=L2, use_batch_norm=use_batch_norm)
            if layer_depth < depth - 1:
                current_layer = MaxPooling2D(pool_size=(2,2))(layer2)
                levels.append([layer1, layer2, current_layer])
            else:
                current_layer = layer2
                levels.append([layer1, layer2])

        for layer_depth in range(depth-2, -1, -1):
            
            up_convolution = get_up_convolution(pool_size=(2,2), deconvolution=deconvolution, n_filters=n_base_filters*(layer_depth+1), L2=L2)(current_layer)

            concat = concatenate([up_convolution, levels[layer_depth][1]] , axis=-1)
            current_layer = create_convolution_block(n_filters=n_base_filters*(layer_depth+1),kernel=(3,3), input_layer=concat, padding='same', L2=L2, use_batch_norm=use_batch_norm)
            current_layer = create_convolution_block(n_filters=n_base_filters*(layer_depth+1),kernel=(3,3), input_layer=current_layer, padding='same', L2=L2, use_batch_norm=use_batch_norm)



        if add_spatial_prior:

            current_layer = tf.keras.layers.concatenate([current_layer, positional_encoding])
            current_layer = Conv2D(128, (3, 3),  padding='same', activation='relu')(current_layer)
            current_layer = tf.keras.layers.concatenate([current_layer, positional_encoding])
            current_layer = Conv2D(256, (1, 1), activation='relu')(current_layer)
            current_layer = tf.keras.layers.concatenate([current_layer,positional_encoding])
            current_layer = Conv2D(512, (1, 1), activation='relu')(current_layer)
            
            # current_layer = tf.keras.layers.concatenate([current_layer, tf.tile(positional_encoding, multiples=[1,256, 256, 1])])
            # current_layer = Conv2D(256, (1, 1), activation='relu')(current_layer)
            # current_layer = tf.keras.layers.concatenate([current_layer, tf.tile(positional_encoding, multiples=[1,256, 256, 1])])
            # current_layer = Conv2D(512, (1, 1), activation='relu')(current_layer)
            
        else:
    
            current_layer = Conv2D(256, (1, 1), activation='relu')(current_layer)
            current_layer = Conv2D(512, (1, 1), activation='relu')(current_layer)
        
        final_convolution = Conv2D(7, (1, 1))(current_layer)   

        act = Activation(activation_name)(final_convolution)
        
        if add_spatial_prior:

            model = Model(inputs=[inputs, positional_encoding_input], outputs=act)
    
        else:
            model = Model(inputs=[inputs], outputs=act)

        
        model.compile(loss=Generalised_dice_coef_multilabel7, 
                  optimizer=Adam(lr=initial_learning_rate), 
                  metrics=[dice_coef_multilabel_bin0, 
                           dice_coef_multilabel_bin1, 
                           dice_coef_multilabel_bin2,
                           dice_coef_multilabel_bin3,
                           dice_coef_multilabel_bin4,
                           dice_coef_multilabel_bin5,
                           dice_coef_multilabel_bin6])


        return model
    


def UNet_v0_2DTumorSegmenter_V2(input_shape =  (256, 256,1), pool_size=(2, 2),initial_learning_rate=1e-5, deconvolution=True,
                      depth=6, n_base_filters=32, activation_name="softmax", L2=0, use_batch_norm=False, add_spatial_prior=False):
        """ Simple version, padding 'same' on every layer, output size is equal to input size. Has border artifacts and checkerboard artifacts """
        inputs = Input(input_shape)
        
        if add_spatial_prior:
            positional_encoding_input = Input((256,256,3))
            positional_encoding = positional_encoding_input
        
        levels = list()
        #current_layer = Conv2D(n_base_filters, (1, 1))(inputs)  # ???? needed??  Not even a nonlinear activation!!!
        current_layer = inputs
        
    
        # add levels with max pooling
        for layer_depth in range(depth):
            layer1 = create_convolution_block(input_layer=current_layer, kernel=(3,3), n_filters=n_base_filters*(layer_depth+1), padding='same', L2=L2, use_batch_norm=use_batch_norm)
            layer2 = create_convolution_block(input_layer=layer1, kernel=(3,3),  n_filters=n_base_filters*(layer_depth+1), padding='same', L2=L2, use_batch_norm=use_batch_norm)
            if layer_depth < depth - 1:
                current_layer = MaxPooling2D(pool_size=(2,2))(layer2)
                levels.append([layer1, layer2, current_layer])
            else:
                current_layer = layer2
                levels.append([layer1, layer2])

        for layer_depth in range(depth-2, -1, -1):
            
            up_convolution = get_up_convolution(pool_size=(2,2), deconvolution=deconvolution, n_filters=n_base_filters*(layer_depth+1), L2=L2)(current_layer)

            concat = concatenate([up_convolution, levels[layer_depth][1]] , axis=-1)
            current_layer = create_convolution_block(n_filters=n_base_filters*(layer_depth+1),kernel=(3,3), input_layer=concat, padding='same', L2=L2, use_batch_norm=use_batch_norm)
            current_layer = create_convolution_block(n_filters=n_base_filters*(layer_depth+1),kernel=(3,3), input_layer=current_layer, padding='same', L2=L2, use_batch_norm=use_batch_norm)



        if add_spatial_prior:

            current_layer = tf.keras.layers.concatenate([current_layer, positional_encoding])
            current_layer = Conv2D(128, (3, 3),  padding='same', activation='relu')(current_layer)
            current_layer = tf.keras.layers.concatenate([current_layer, positional_encoding])
            current_layer = Conv2D(256, (1, 1), activation='relu')(current_layer)
            current_layer = tf.keras.layers.concatenate([current_layer,positional_encoding])
            current_layer = Conv2D(512, (1, 1), activation='relu')(current_layer)
            

        else:
    
            current_layer = Conv2D(256, (1, 1), activation='relu')(current_layer)
            current_layer = Conv2D(512, (1, 1), activation='relu')(current_layer)
        
        final_convolution = Conv2D(7, (1, 1))(current_layer)   

        act = Activation(activation_name)(final_convolution)
        
        if add_spatial_prior:

            model = Model(inputs=[inputs, positional_encoding_input], outputs=act)
    
        else:
            model = Model(inputs=[inputs], outputs=act)

        
        model.compile(loss=Generalised_dice_coef_multilabel7, 
                  optimizer=Adam(lr=initial_learning_rate), 
                  metrics=[dice_coef_multilabel_bin0, 
                           dice_coef_multilabel_bin1, 
                           dice_coef_multilabel_bin2,
                           dice_coef_multilabel_bin3,
                           dice_coef_multilabel_bin4,
                           dice_coef_multilabel_bin5,
                           dice_coef_multilabel_bin6])


        return model


#%%

def compute_positional_matrix(batch_size, image_size):
        
    z_grid = np.full((image_size,image_size), np.arange(image_size))
    
    half = image_size//2
    
    center = 0
    x_grid, y_grid = np.meshgrid(np.arange(-half,half), np.arange(-half,half))
    x_grid = x_grid/image_size
    y_grid = y_grid/image_size
    
    distances = np.sqrt(((x_grid - center) ** 2) + ((y_grid - center) ** 2))
    center = distances/np.max(distances)
    
    
    center = np.repeat(center[ np.newaxis,:], batch_size, axis=0), 
    x_grid = np.repeat(x_grid[ np.newaxis,:], batch_size, axis=0), 
    y_grid = np.repeat(y_grid[ np.newaxis,:], batch_size, axis=0)

    return np.stack([z_grid,center,x_grid,y_grid],-1)



class DataGenerator(tf.keras.utils.Sequence): # inheriting from Sequence allows for multiprocessing functionalities

    def __init__(self, list_IDs, batch_size=4, dim=(256,256,1), n_channels=3,
                 n_classes=2, shuffledata=True, data_path='', labels_path='', 
                 do_augmentation=True, use_slice_location = True, debug=False):
        
        self.dim = dim
        self.batch_size = batch_size
        self.list_IDs = list_IDs
        self.n_channels = n_channels
        self.n_classes = n_classes
        self.shuffledata = shuffledata
        self.on_epoch_end()
        self.data_path = data_path
        self.labels_path = labels_path
        self.seed = 0
        self.do_augmentation = do_augmentation
        self.debug = debug
        self.use_slice_location = use_slice_location
        
        self.augmentor = tf.keras.preprocessing.image.ImageDataGenerator(
                    rotation_range=0,
                    shear_range=0.0,
                    zoom_range=0.0,
                    horizontal_flip=False,
                    vertical_flip=False,
                    fill_mode='nearest',  # ????
                    
                )
        
        self.augmentor_mask = tf.keras.preprocessing.image.ImageDataGenerator( 
                    rotation_range=0,
                    shear_range=0.0,
                    zoom_range=0.0,
                    horizontal_flip=False,
                    vertical_flip=False,
                    fill_mode='constant'
                    
#                    preprocessing_function = lambda x: np.where(x>0, 1, 0).astype('float32')
                    
#                    preprocessing_function = lambda x:np.where(x < 0.5, 0, 
#                                                                (np.where(x < 1.5, 1, 
#                                                                (np.where(x < 2.5, 2, 
#                                                                (np.where(x < 3.5, 3, 
#                                                                (np.where(x < 4.5, 4, 
#                                                                (np.where(x < 5.5, 5, 6))))))))))).astype('float32') 
                    )
        
        
        
        self.center_matrix, self.x_grid, self.y_grid = self._compute_positional_matrix(self.batch_size)
        
    def _compute_positional_matrix(self, batch_size):
        
       
        
        center = 0
        x_grid, y_grid = np.meshgrid(np.arange(-128,128), np.arange(-128,128))
        x_grid = x_grid/256
        y_grid = y_grid/256
        
        distances = np.sqrt(((x_grid - center) ** 2) + ((y_grid - center) ** 2))
        distance_from_center = distances/np.max(distances)

        
        return  np.repeat(distance_from_center[ np.newaxis,:], batch_size, axis=0), np.repeat(x_grid[ np.newaxis,:], batch_size, axis=0), np.repeat(y_grid[ np.newaxis,:], batch_size, axis=0)

    def __len__(self):
        'Denotes the number of batches per epoch'
        return int(np.floor(len(self.list_IDs) / self.batch_size))

    def __getitem__(self, index):
        'Generate one batch of data'
        # Generate indexes of the batch
        indexes = self.indexes[index*self.batch_size:(index+1)*self.batch_size]

        # Find list of IDs
        list_IDs_temp = [self.list_IDs[k] for k in indexes]

        # Generate data
        X, y, pos = self.__data_generation(list_IDs_temp)
        
        # if self.do_augmentation:
        #     for i in range(self.batch_size):
        #         X[i] *= np.random.uniform(low=0.9, high=1.1, size=1)                                                     
                    
        y = tf.keras.utils.to_categorical(y, num_classes=self.n_classes) 

        if self.use_slice_location:
            positional_encoding = np.stack([pos,self.center_matrix, self.x_grid, self.y_grid],-1)

        if self.debug:
            
            return [X, positional_encoding], y, list_IDs_temp
        else:
            if self.use_slice_location:

                return [X, positional_encoding], y

            else:
                return X, y

    def on_epoch_end(self):
        'Updates indexes after each epoch'
        self.indexes = np.arange(len(self.list_IDs))
        if self.shuffledata == True:
            np.random.shuffle(self.indexes)


    def __data_generation(self, list_IDs_temp):
        'Generates data containing batch_size samples' # X : (n_samples, *dim, n_channels)
        # Initialization
        X = np.empty((self.batch_size, self.dim[0], self.dim[1], self.n_channels))#, dtype='float32')
        y = np.zeros((self.batch_size, self.dim[0], self.dim[1]))
        
        if self.use_slice_location:
            positional_encoding_vector = np.empty((self.batch_size,256,256))
            
        # Generate data
        # orientations = ['sagittal','coronal','axial']
        for i, ID in enumerate(list_IDs_temp):          
            # X[i,:,:,0] = np.load(self.data_path.replace('sagittal',orientations[i%3]) + ID)   # Here we add the path. ID can be the path
            # y[i] = np.load(self.labels_path.replace('sagittal',orientations[i%3])  + ID)

            X[i,:,:,0] = np.load(self.data_path + ID)   # Here we add the path. ID can be the path
            y[i] = np.load(self.labels_path  + ID)

            
            SLICE = int(ID.split('_')[-1].split('.npy')[0].replace('slice',''))

            if self.use_slice_location:
                positional_encoding_vector[i] =  np.full((256, 256), SLICE/256.) 

        #y = y-1
        y = np.expand_dims(y, -1) #ImageDataGenerator needs arrays of rank 4. But wants channel dim = 1,3,4
        
        if self.do_augmentation:
            X_gen = self.augmentor.flow(X, batch_size=self.batch_size, shuffle=False, seed=self.seed)
            y_gen = self.augmentor_mask.flow(y, batch_size=self.batch_size, shuffle=False, seed=self.seed)


            return next(X_gen), next(y_gen)
        else:
            return X,y,positional_encoding_vector


#%%  DATA GENERATOR 2



class DataGenerator2(tf.keras.utils.Sequence): # inheriting from Sequence allows for multiprocessing functionalities

    def __init__(self, list_IDs, batch_size=4, dim=(256,256,1), n_channels=3,
                 n_classes=2, shuffledata=True, data_path='', labels_path='', coords_path=None, 
                 do_augmentation=True, axis=None, rotation_range=5, shear_range=5, debug=False):
        
        self.dim = dim
        self.batch_size = batch_size
        self.list_IDs = list_IDs
        self.n_channels = n_channels
        self.n_classes = n_classes
        self.shuffledata = shuffledata
        self.on_epoch_end()
        self.data_path = data_path
        self.labels_path = labels_path
        self.coords_path = coords_path
        self.seed = 0
        self.do_augmentation = do_augmentation
        self.debug = debug
        self.axis = axis
        self.rotation_range = rotation_range
        self.shear_range = shear_range
        
        # self.augmentor = tf.keras.preprocessing.image.ImageDataGenerator(
        #             rotation_range=self.rotation_range,
        #             shear_range=self.shear_range,
        #             fill_mode='constant', 
        #             interpolation_order=1,
        #         )
        
        # self.augmentor_mask = tf.keras.preprocessing.image.ImageDataGenerator( 
        #             rotation_range=self.rotation_range,
        #             shear_range=self.shear_range,
        #             fill_mode='constant', cval=0,
        #             interpolation_order=0,
        # )
        
        
        # self.augmentor_coordinates = tf.keras.preprocessing.image.ImageDataGenerator(
        #             rotation_range=self.rotation_range,
        #             shear_range=self.shear_range,
        #             interpolation_order=1,
        #             fill_mode='nearest',  # MRI, segmentation = constant.   Coordinates = nearest.  
                    
        #         )
        
   
    def __len__(self):
        'Denotes the number of batches per epoch'
        return int(np.floor(len(self.list_IDs) / self.batch_size))

    def __getitem__(self, index):
        'Generate one batch of data'
        # Generate indexes of the batch
        indexes = self.indexes[index*self.batch_size:(index+1)*self.batch_size]

        # Find list of IDs
        list_IDs_temp = [self.list_IDs[k] for k in indexes]

        # Generate data
        X, y = self.__data_generation(list_IDs_temp)
                
        y = tf.keras.utils.to_categorical(y, num_classes=self.n_classes) 

        if self.debug:
            return X, y, list_IDs_temp 
        else:
            return X, y

    def on_epoch_end(self):
        'Updates indexes after each epoch'
        self.indexes = np.arange(len(self.list_IDs))
        if self.shuffledata == True:
            np.random.shuffle(self.indexes)


    def __data_generation(self, list_IDs_temp):
        'Generates data containing batch_size samples' # X : (n_samples, *dim, n_channels)
        # Initialization
        X = np.empty((self.batch_size, self.dim[0], self.dim[1], self.n_channels))#, dtype='float32')
        y = np.zeros((self.batch_size, self.dim[0], self.dim[1]))
        
        if self.coords_path is not None:
            positional_encoding_vector = np.empty((self.batch_size,256,256, 3))
            

        for i, ID in enumerate(list_IDs_temp):          
            if self.axis == 'all':
                random_axis = np.random.choice(['sagittal','axial','coronal'], size=1)[0]
                X[i,:,:,0] = np.load(self.data_path.replace('sagittal',random_axis) + ID)   # Here we add the path. ID can be the path
                y[i] = np.load(self.labels_path.replace('sagittal',random_axis)  + ID)
            else:
                X[i,:,:,0] = np.load(self.data_path + ID)   # Here we add the path. ID can be the path
                y[i] = np.load(self.labels_path  + ID)

            if self.coords_path is not None:
                if self.axis == 'all':
                    positional_encoding_vector[i] =  np.load(self.coords_path.replace('sagittal',random_axis) + ID) /256.
                else:
                    positional_encoding_vector[i] =  np.load(self.coords_path + ID) /256.

        y = np.expand_dims(y, -1) #ImageDataGenerator needs arrays of rank 4. But wants channel dim = 1,3,4

        if self.do_augmentation:
            X_gen = self.augmentor.flow(X, batch_size=self.batch_size, shuffle=False, seed=self.seed)
            y_gen = self.augmentor_mask.flow(y, batch_size=self.batch_size, shuffle=False, seed=self.seed)

            if self.coords_path is not None:
                positional_encoding_vector_gen = self.augmentor_coordinates.flow(positional_encoding_vector, batch_size=self.batch_size, shuffle=False, seed=self.seed)
                return [next(X_gen), next(positional_encoding_vector_gen)], next(y_gen), 
                return [next(X_gen), next(positional_encoding_vector_gen)], next(y_gen),
            else:
                return next(X_gen), next(y_gen)
        else:
            if self.coords_path is not None:
                return [X,positional_encoding_vector],y
            else:
                return X,y


#%%

class MyHistory(tf.keras.callbacks.Callback):
    def __init__(self, OUT, NAME, loss=[], val_loss=[], 
                 dice_coef_multilabel_bin0=[], dice_coef_multilabel_bin1=[], 
                 dice_coef_multilabel_bin2=[], dice_coef_multilabel_bin3=[], 
                 dice_coef_multilabel_bin4=[], dice_coef_multilabel_bin5=[], dice_coef_multilabel_bin6=[], 
                 val_dice_coef_multilabel_bin0=[], val_dice_coef_multilabel_bin1=[], 
                 val_dice_coef_multilabel_bin2=[], val_dice_coef_multilabel_bin3=[], 
                 val_dice_coef_multilabel_bin4=[], val_dice_coef_multilabel_bin5=[], val_dice_coef_multilabel_bin6=[]):
        
        self.OUT = OUT
        self.NAME = NAME  
        
        self.loss = loss
        self.dice_coef_multilabel_bin0 = dice_coef_multilabel_bin0
        self.dice_coef_multilabel_bin1 = dice_coef_multilabel_bin1
        self.dice_coef_multilabel_bin2 = dice_coef_multilabel_bin2
        self.dice_coef_multilabel_bin3 = dice_coef_multilabel_bin3
        self.dice_coef_multilabel_bin4 = dice_coef_multilabel_bin4
        self.dice_coef_multilabel_bin5 = dice_coef_multilabel_bin5
        self.dice_coef_multilabel_bin6 = dice_coef_multilabel_bin6

        self.val_loss = val_loss
        self.val_dice_coef_multilabel_bin0 = val_dice_coef_multilabel_bin0
        self.val_dice_coef_multilabel_bin1 = val_dice_coef_multilabel_bin1
        self.val_dice_coef_multilabel_bin2 = val_dice_coef_multilabel_bin2
        self.val_dice_coef_multilabel_bin3 = val_dice_coef_multilabel_bin3
        self.val_dice_coef_multilabel_bin4 = val_dice_coef_multilabel_bin4
        self.val_dice_coef_multilabel_bin5 = val_dice_coef_multilabel_bin5
        self.val_dice_coef_multilabel_bin6 = val_dice_coef_multilabel_bin6        
        
    def on_train_begin(self, logs={}):

        self.loss = []
        self.dice_coef_multilabel_bin0 = []
        self.dice_coef_multilabel_bin1 = []
        self.dice_coef_multilabel_bin2 = []
        self.dice_coef_multilabel_bin3 = []
        self.dice_coef_multilabel_bin4 = []
        self.dice_coef_multilabel_bin5 = []
        self.dice_coef_multilabel_bin6 = []

        self.val_loss = []
        self.val_dice_coef_multilabel_bin0 = []
        self.val_dice_coef_multilabel_bin1 = []
        self.val_dice_coef_multilabel_bin2 = []
        self.val_dice_coef_multilabel_bin3 = []
        self.val_dice_coef_multilabel_bin4 = []
        self.val_dice_coef_multilabel_bin5 = []
        self.val_dice_coef_multilabel_bin6 = []

    # def on_batch_end(self, batch, logs={}):
    #     self.loss.append(logs.get('loss'))
    #     self.dice_coef_multilabel_bin1.append(logs.get('dice_coef_multilabel_bin1'))
        
    def on_epoch_end(self, epoch, logs={}):

        self.loss.append(logs.get('loss'))
        self.dice_coef_multilabel_bin0.append(logs.get('dice_coef_multilabel_bin0'))
        self.dice_coef_multilabel_bin1.append(logs.get('dice_coef_multilabel_bin1'))
        self.dice_coef_multilabel_bin2.append(logs.get('dice_coef_multilabel_bin2'))
        self.dice_coef_multilabel_bin3.append(logs.get('dice_coef_multilabel_bin3'))
        self.dice_coef_multilabel_bin4.append(logs.get('dice_coef_multilabel_bin4'))
        self.dice_coef_multilabel_bin5.append(logs.get('dice_coef_multilabel_bin5'))
        self.dice_coef_multilabel_bin6.append(logs.get('dice_coef_multilabel_bin6'))


        self.val_loss.append(logs.get('val_loss'))
        self.val_dice_coef_multilabel_bin0.append(logs.get('val_dice_coef_multilabel_bin0'))
        self.val_dice_coef_multilabel_bin1.append(logs.get('val_dice_coef_multilabel_bin1'))
        self.val_dice_coef_multilabel_bin2.append(logs.get('val_dice_coef_multilabel_bin2'))
        self.val_dice_coef_multilabel_bin3.append(logs.get('val_dice_coef_multilabel_bin3'))
        self.val_dice_coef_multilabel_bin4.append(logs.get('val_dice_coef_multilabel_bin4'))
        self.val_dice_coef_multilabel_bin5.append(logs.get('val_dice_coef_multilabel_bin5'))
        self.val_dice_coef_multilabel_bin6.append(logs.get('val_dice_coef_multilabel_bin6'))

        # self.loss.append(logs.get('loss'))
        # self.dice_coef_multilabel_bin0.append(logs.get('dice_metric0'))
        # self.dice_coef_multilabel_bin1.append(logs.get('dice_metric1'))
        # self.dice_coef_multilabel_bin2.append(logs.get('dice_metric2'))
        # self.dice_coef_multilabel_bin3.append(logs.get('dice_metric3'))
        # self.dice_coef_multilabel_bin4.append(logs.get('dice_metric4'))
        # self.dice_coef_multilabel_bin5.append(logs.get('dice_metric5'))
        # self.dice_coef_multilabel_bin6.append(logs.get('dice_metric6'))


        # self.val_loss.append(logs.get('val_loss'))
        # self.val_dice_coef_multilabel_bin0.append(logs.get('val_dice_metric0'))
        # self.val_dice_coef_multilabel_bin1.append(logs.get('val_dice_metric1'))
        # self.val_dice_coef_multilabel_bin2.append(logs.get('val_dice_metric2'))
        # self.val_dice_coef_multilabel_bin3.append(logs.get('val_dice_metric3'))
        # self.val_dice_coef_multilabel_bin4.append(logs.get('val_dice_metric4'))
        # self.val_dice_coef_multilabel_bin5.append(logs.get('val_dice_metric5'))
        # self.val_dice_coef_multilabel_bin6.append(logs.get('val_dice_metric6'))


        plt.figure(figsize=(15,5))
        
        plt.subplot(1,3,1); 
        plt.title('Loss')
        plt.plot(self.loss)
        plt.plot(self.val_loss)
        plt.legend(['Train','Val'])
        plt.grid()
        
        plt.subplot(1,3,2); plt.title('Dice')
        plt.plot(self.dice_coef_multilabel_bin0, label='Bg')
        plt.plot(self.dice_coef_multilabel_bin1, label='Air')
        plt.plot(self.dice_coef_multilabel_bin2, label='GM')
        plt.plot(self.dice_coef_multilabel_bin3, label='WM')
        plt.plot(self.dice_coef_multilabel_bin4, label='CSF')
        plt.plot(self.dice_coef_multilabel_bin5, label='Bone')
        plt.plot(self.dice_coef_multilabel_bin6, label='Skin')
        plt.ylim([0,1])
        plt.legend(loc='upper left', bbox_to_anchor=(0.95, 0.5), ncol=1, fancybox=True, shadow=True)
        plt.grid()
        
        
        plt.subplot(1,3,3); plt.title('Dice Val')
        plt.plot(self.val_dice_coef_multilabel_bin0, label='Bg')
        plt.plot(self.val_dice_coef_multilabel_bin1, label='Air')
        plt.plot(self.val_dice_coef_multilabel_bin2, label='GM')
        plt.plot(self.val_dice_coef_multilabel_bin3, label='WM')
        plt.plot(self.val_dice_coef_multilabel_bin4, label='CSF')
        plt.plot(self.val_dice_coef_multilabel_bin5, label='Bone')
        plt.plot(self.val_dice_coef_multilabel_bin6, label='Skin')
        #plt.legend(loc='upper left', bbox_to_anchor=(0.95, 0.5), ncol=1, fancybox=True, shadow=True)
        plt.grid()
        plt.ylim([0,1])
        
        plt.tight_layout()   
        
        plt.savefig(self.OUT + self.NAME + '/Training_curves.png')
        plt.close()

def freeze_layers(model):
    model_type = type(model) 

    for i in model.layers:
        i.trainable = False
        if type(i) == model_type:
            freeze_layers(i)
    return model


def save_model_and_weights(model, NAME, FOLDER):    
    model_to_save = tf.keras.models.clone_model(model)
    model_to_save.set_weights(model.get_weights())
    model_to_save = freeze_layers(model_to_save)
    model_to_save.save_weights(FOLDER + NAME + '_weights.h5')
    model_to_save.save(FOLDER + NAME + '.h5')       

class my_model_checkpoint(tf.keras.callbacks.Callback):
    
    def __init__(self, MODEL_PATH, MODEL_NAME, val_loss = [999]):
        self.MODEL_PATH = MODEL_PATH
        self.MODEL_NAME = MODEL_NAME    
        self.val_loss = val_loss

    def on_epoch_end(self, epoch, logs={}):
        min_val_loss = min(self.val_loss)
        current_val_loss = logs.get('val_loss')
        self.val_loss.append(current_val_loss)
        print('Min loss so far: {}, new loss: {}'.format(min_val_loss, current_val_loss))
        if current_val_loss < min_val_loss :
            print('New best model! Epoch: {}'.format(epoch))
            save_model_and_weights(self.model, self.MODEL_NAME, self.MODEL_PATH)
        else :
            save_model_and_weights(self.model, '/last_model', self.MODEL_PATH)            




#%%
def segment_MRI(img, coords, model_sagittal=None, model_axial=None, model_coronal=None,model_layer=None):
    model_segmentation_sagittal = None
    model_segmentation_coronal = None
    model_segmentation_axial = None

    print('SEGMENTING...will take a while...')
    # Predict using sagittal model
    if model_sagittal is not None:
        print('working on the sagittal slices...')
        yhat_sagittal = model_sagittal.predict([np.expand_dims(img, -1), coords], batch_size=1,verbose=2) if coords is not None else model_sagittal.predict(np.expand_dims(img, -1), batch_size=1,verbose=2)
            
    # Predict using coronal model
    if model_coronal is not None:
        print('working on the coronal slices...')
        img_coronal = np.swapaxes(img, 0, 1)
        coords_coronal = np.swapaxes(coords, 0, 1) if coords is not None else None
        yhat_coronal = model_coronal.predict([np.expand_dims(img_coronal, -1), coords_coronal], batch_size=1,verbose=2) if coords is not None else model_coronal.predict(np.expand_dims(img_coronal, -1), batch_size=1,verbose=2)
        yhat_coronal = np.swapaxes(yhat_coronal, 0, 1)  # Align to original orientation
            
    # Predict using axial model
    if model_axial is not None:
        print('working on the axial slices...')
        img_axial = np.swapaxes(np.swapaxes(img, 1, 2), 0, 1)
        coords_axial = np.swapaxes(np.swapaxes(coords, 1, 2), 0, 1) if coords is not None else None
        yhat_axial = model_axial.predict([np.expand_dims(img_axial, -1), coords_axial], batch_size=1,verbose=2) if coords is not None else model_axial.predict(np.expand_dims(img_axial, -1), batch_size=1,verbose=2)
        yhat_axial = np.swapaxes(np.swapaxes(yhat_axial, 0, 1), 1, 2)  # Align to original orientation
        
    #print('done 3 views!')
    # Prepare input data for the new model using probability maps
    try:
        if model_layer is not None:
            print('now combining results from the 3 views...')
            X_test = np.concatenate([np.expand_dims(img / np.percentile(img, 95), -1), yhat_sagittal, yhat_coronal, yhat_axial], axis=-1)
            # X_test = np.concatenate([np.expand_dims(img, -1), yhat_sagittal, yhat_coronal, yhat_axial], axis=-1)
            # X_test = np.concatenate([yhat_sagittal, yhat_coronal, yhat_axial], axis=-1)
            yhat_new_model = model_layer.predict(np.expand_dims(X_test, 0),verbose=2)
            model_segmentation_layer = np.argmax(yhat_new_model[0], -1).astype(np.int32)
        else:
	        model_segmentation_layer = None
    except:
        print('failed')
    #print('check here?')
        
    # Return all predictions
    return model_segmentation_layer
 


 

