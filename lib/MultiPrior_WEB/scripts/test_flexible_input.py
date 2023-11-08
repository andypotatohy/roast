# -*- coding: utf-8 -*-
"""
Created on Mon Jan 28 14:50:13 2019

@author: hirsch
"""


import os
import sys
import nibabel as nib
import numpy as np
import time
import random
from numpy.random import seed
import keras 
from keras import backend as K
from keras.utils import to_categorical
from tensorflow import set_random_seed
from sklearn import metrics
import matplotlib.pyplot as plt
from shutil import copyfile
from random import shuffle
from random import sample
from keras.callbacks import ModelCheckpoint
from matplotlib.patches import Rectangle
import pandas as pd
from keras.utils import print_summary

from lib import *

configFile = '/home/hirsch/Documents/projects/brainSegmentation/DeepPriors/configFiles/configFiles_stroke/configFile0_BIG_CNN_TPM_Dice_loss.py'
workingDir ='/home/hirsch/Documents/projects/brainSegmentation/DeepPriors/'

print(configFile)
path = '/'.join(configFile.split('/')[:-1])
print(path)
configFileName = configFile.split('/')[-1][:-3]   
sys.path.append(path)
cfg = __import__(configFileName)

cfg.TPM_channel = workingDir + cfg.TPM_channel
cfg.trainChannels = [workingDir + x for x in cfg.trainChannels]
cfg.trainLabels = workingDir +cfg.trainLabels 
cfg.testChannels = [workingDir + x for x in cfg.testChannels]
cfg.testLabels = workingDir + cfg.testLabels
cfg.validationChannels = [workingDir + x for x in cfg.validationChannels]
cfg.validationLabels = workingDir +cfg.validationLabels

cfg.num_channels = 1

cfg.dpatch = [61,61,61]

from BIG_multiscale_CNN_TPM_flexibleInput import BIG_multiscale_CNN_TPM_flexibleInput
dm = BIG_multiscale_CNN_TPM_flexibleInput(cfg.dpatch, cfg.output_classes, cfg.num_channels, cfg.L2, cfg.dropout, cfg.learning_rate, cfg.optimizer_decay, cfg.loss_function)
model = dm.createModel()            
print_summary(model, positions=[.33, .8, .9,1])

from keras.utils import plot_model
plot_model(model, to_file='/home/hirsch/Documents/projects/Breast_segmentation/DeepPriors_package/multiscale_TPM.png', show_shapes=True)

# Dummy data
x = np.random.rand(1,80,350,350,3)

x = np.random.rand(1,61,61,61,1)
x.shape
TPM = np.random.rand(1,9,9,9,4)
yhat = model.predict([x,TPM])
yhat.shape

dpatch = 181
tpm_size = dpatch - 52

x = np.random.rand(1,dpatch,dpatch,dpatch,1)
x.shape
TPM = np.random.rand(1,tpm_size,tpm_size,tpm_size,4)
yhat = model.predict([x,TPM])
yhat.shape

## Real MRI


nii = nib.load('/home/hirsch/Documents/projects/MSKCC/FINAL_TRAIN_DATA/test/highres/rrMSKCC_16-328_1_00113_20031116_T1_post_rsc_stand_stand_ras.nii')
X1 = nii.get_data()

nii2 = nib.load('/home/hirsch/Documents/projects/MSKCC/FINAL_TRAIN_DATA/test/lowres/rrMSKCC_16-328_1_00113_20031116_T1_post_rsc_stand_stand_ras.nii')
X2 = nii2.get_data()

X1.shape
X2.shape

padding = (0,0,0)

output_shape = X1.shape
X2.shape

X = np.zeros((X1.shape[0] + padding[0],X1.shape[1] + padding[1],X1.shape[2] + padding[2]))
XX = np.zeros((X1.shape[0] + padding[0],X1.shape[1] + padding[1],X1.shape[2] + padding[2]))


X[(X.shape[0] - X1.shape[0])/2 : X.shape[0] - (X.shape[0] - X1.shape[0])/2 ,
  (X.shape[1] - X1.shape[1])/2 : X.shape[1] - (X.shape[1] - X1.shape[1])/2 ,
  (X.shape[2] - X1.shape[2])/2 : X.shape[2] - (X.shape[2] - X1.shape[2])/2 ] = X1


XX[(XX.shape[0] - X2.shape[0])/2 : XX.shape[0] - (XX.shape[0] - X2.shape[0])/2 -1,
   (XX.shape[1] - X2.shape[1])/2 : XX.shape[1] - (XX.shape[1] - X2.shape[1])/2 -1,
   (XX.shape[2] - X2.shape[2])/2 : XX.shape[2] - (XX.shape[2] - X2.shape[2])/2 -1] = X2



X = X.reshape((1,) + X.shape + (1,) )
XX = XX.reshape((1,) + XX.shape + (1,) )

X.shape
XX.shape

start = time.time()
yhat = model.predict([X,XX])
end = time.time()
print(end - start)
yhat.shape


Y_nii = nib.load('/home/hirsch/Documents/projects/MSKCC/FINAL_TRAIN_DATA/test/highres/rrMSKCC_16-328_1_00113_20031116_label_ras.nii')
Y = Y_nii.get_data()
Y.shape
n_values = np.max(Y) +1
y_target = np.eye(n_values)[Y]
y_target = y_target.reshape((1,) + y_target.shape)

model.fit([X,XX],y_target)

y = np.argmax(yhat, axis=4)
#y = yhat[:,:,:,:,1]
y.shape
y_out = np.zeros((1,) + output_shape)
y_out[:,(y_out.shape[1] -y.shape[1])/2:y_out.shape[1] - (y_out.shape[1] -y.shape[1])/2,
        (y_out.shape[2] -y.shape[2])/2:y_out.shape[2] - (y_out.shape[2] -y.shape[2])/2,
        (y_out.shape[3] -y.shape[3])/2:y_out.shape[3] - (y_out.shape[3] -y.shape[3])/2] = y


y_out = y_out.reshape(output_shape)
y_out.shape
img = nib.Nifti1Image(y_out, nii.affine)
nib.save(img,'/home/hirsch/Documents/projects/Breast_segmentation/DeepPriors_package/example_FullInput.nii')