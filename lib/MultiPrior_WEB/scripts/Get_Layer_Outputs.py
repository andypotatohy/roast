# -*- coding: utf-8 -*-
"""
Created on Mon Aug 26 11:39:29 2019

@author: hirsch
"""


import nibabel as nib
import numpy as np
import sys
sys.path.append('/home/hirsch/Documents/projects/brainSegmentation/DeepPriors/scripts/')
from MultiPriors_noDownsampling import Generalised_dice_coef_multilabel7, dice_coef_multilabel0,dice_coef_multilabel1,dice_coef_multilabel2,dice_coef_multilabel3,dice_coef_multilabel4,dice_coef_multilabel5,dice_coef_multilabel6
from keras.models import load_model  
import os
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 16})


SESSION_PATH = '/home/hirsch/Documents/projects/brainSegmentation/DeepPriors/training_sessions/MultiPriors_noDownsampling_fullHeadSegmentation_configFile_Aphasic_Stroke_MultiPriors_noDownsampling_2019-08-23_2109'

MODEL_PATH = SESSION_PATH + '/models/noDownsampling_fullHeadSegmentation_configFile_Aphasic_Stroke_MultiPriors_noDownsampling_2019-08-23_2109.log_epoch91.h5'

X          = nib.load('/home/hirsch/Documents/projects/stroke_heads/AphasicStrokeTrial/GU008/GU008_Avg_for_dixon_big.nii').get_data()

TPM        = nib.load('/home/hirsch/Documents/projects/stroke_heads/AphasicStrokeTrial/GU008/indiTPM.nii').get_data()

save_plot = True

subject_name = 'GU008'

LAYER_NAME_LIST = ['conv3d_3','conv3d_8','conv3d_10','conv3d_16','up_sampling3d_1','conv3d_17','conv3d_18','conv3d_19']


def percentile95_normalizeMRI(data):
    p95 = np.percentile(data,95)
    data = data/p95
    return(data)    



my_custom_objects = {'Generalised_dice_coef_multilabel7':Generalised_dice_coef_multilabel7,
                                 'dice_coef_multilabel0':dice_coef_multilabel0,
                                 'dice_coef_multilabel1':dice_coef_multilabel1,
                                 'dice_coef_multilabel2':dice_coef_multilabel2,
                                 'dice_coef_multilabel3':dice_coef_multilabel3,
                                 'dice_coef_multilabel4':dice_coef_multilabel4,
                                 'dice_coef_multilabel5':dice_coef_multilabel5,
                                 'dice_coef_multilabel6':dice_coef_multilabel6}
                                 
model = load_model(MODEL_PATH, custom_objects = my_custom_objects )

BORDER = 0
AXIAL_BORDER = 0
CORONAL_BORDER = 81

X = np.pad(X, ((30,30),(30,30),(0,0)), 'minimum')

pad = []
for ch in range(TPM.shape[3]):
  pad.append(np.pad(TPM[:,:,:,ch],((30,30),(30,30),(0,0)),'minimum'))
TPM = np.stack(pad, axis=3)


# Make shape divisible by 3:
X = X[BORDER:X.shape[0] - X.shape[0]%3 -BORDER, 
      AXIAL_BORDER:X.shape[1] - X.shape[1]%3-AXIAL_BORDER, 
      CORONAL_BORDER:X.shape[2] - X.shape[2]%3-CORONAL_BORDER] 
      
     
TPM = TPM[BORDER:TPM.shape[0] - TPM.shape[0]%3-BORDER, 
          AXIAL_BORDER:TPM.shape[1] - TPM.shape[1]%3-AXIAL_BORDER, 
          CORONAL_BORDER:TPM.shape[2] - TPM.shape[2]%3-CORONAL_BORDER] 

X = percentile95_normalizeMRI(X)

X = X.reshape((1,) + X.shape + (1,))

#X = np.random.randn(1,150,210,210,1)
TPM = TPM[24:-24,24:-24,24:-24,:]
TPM = TPM.reshape((1,) + TPM.shape)

X.shape
TPM.shape

from skimage.transform import resize
lowRes = resize(X, output_shape=[X.shape[0], X.shape[1]/3, X.shape[2]/3, X.shape[3]/3, X.shape[4]], anti_aliasing=True)
lowRes.shape
highRes = X[:,16:-16,16:-16,16:-16,:]
highRes.shape
yhat = model.predict([highRes,lowRes, TPM])
yhat.shape
out = np.argmax(yhat, axis=4)

plt.figure(figsize=(18,10))
plt.subplot(2,2,1)
plt.xticks([])
plt.yticks([])
plt.xlabel('HighRes')
plt.imshow(np.rot90(highRes[0,:,:,highRes.shape[3]/2,0]), cmap = 'gray')

plt.subplot(2,2,2)
plt.xlabel('LowRes')
plt.xticks([])
plt.yticks([])
plt.imshow(np.rot90(lowRes[0,:,:,lowRes.shape[3]/2,0]), cmap = 'gray')

plt.subplot(2,2,3)
plt.xlabel('TPM Channel')
plt.xticks([])
plt.yticks([])
plt.imshow(np.rot90(TPM[0,:,:,TPM.shape[3]/2,0]), cmap = 'gray')

plt.subplot(2,2,4)
plt.xlabel('Model Output')
plt.xticks([])
plt.yticks([])
plt.imshow(np.rot90(out[0,:,:,out.shape[3]/2]))
plt.tight_layout()

if save_plot:
  if not os.path.exists(SESSION_PATH + '/Input_Output_{}.png'.format(subject_name)):  
    plt.savefig(SESSION_PATH + '/Input_Output_{}.png'.format(subject_name))
    plt.close()

from keract import get_activations

for layer_name in LAYER_NAME_LIST:

  print(layer_name)
  x = [highRes,lowRes, TPM]
  
  activations = get_activations(model, x, layer_name)
  
  key = activations.keys()[0]
  activations[key].shape
    
  plt.figure(figsize=(18,10))
  plt.suptitle(layer_name)
  for j in range(1,min(activations[key].shape[-1], 11)):
    plt.subplot(2,5,j)
    plt.axis('off')
    plt.imshow(np.rot90(activations[key][0,:,:,activations[key].shape[3]/2,j-1]), cmap='gray')
  
  plt.tight_layout()
  
  if save_plot:
    plt.savefig(SESSION_PATH + '/feature_maps_{}_{}.png'.format(layer_name, subject_name))
    plt.close()
  
