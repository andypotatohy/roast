# -*- coding: utf-8 -*-
"""
Created on Wed Oct  4 09:00:56 2023

@author: Andrew
"""

import tensorflow as tf
config = tf.compat.v1.ConfigProto()
config.gpu_options.allow_growth = True
config.gpu_options.visible_device_list="0"
tf.compat.v1.keras.backend.set_session(tf.compat.v1.Session(config=config))


import os
import sys
wd = os.getcwd()
subj = sys.argv[2]

############################## Load dataset #############################  
 
segmentChannels = subj.replace('\\','/')  
TPM_channel = os.path.splitext(subj)[0] + '_indiTPM.nii.gz'
 
segmentLabels = ''

output_classes = 7
    
#-------------------------------------------------------------------------------------------------------------

# Parameters 

######################################### MODEL PARAMETERS
# Models : 'CNN_TPM' , 'DeepMedic', 'BIG_multiscale_CNN_TPM_flexible', 'BIG_singleScale_CNN_TPM'
model = 'MultiPriors_noDownsampling'
#dpatch=61
segmentation_dpatch = 51*3
path_to_model = (wd +'/lib/multipriors/multipriors_best_model.h5').replace('\\','/')
session =  'Renamed_Session'

########################################### TEST PARAMETERS
quick_segmentation = True
test_subjects = 19
n_fullSegmentations = test_subjects
size_test_minibatches = 1
saveSegmentation = True

import numpy as np
penalty_MATRIX = np.eye(output_classes,output_classes,dtype='float32')

comments = ''

