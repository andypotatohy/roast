# -*- coding: utf-8 -*-
"""
Created on Fri Aug  9 16:52:09 2019

@author: hirsch
"""

import nibabel as nib
import numpy as np

label = open('/home/hirsch/Documents/projects/brainSegmentation/DeepPriors/CV_folds/aphasic_stroke/all_LABEL.txt').readlines()
mri  = open('/home/hirsch/Documents/projects/brainSegmentation/DeepPriors/CV_folds/aphasic_stroke/all_MRI.txt').readlines()
tpm  = open('/home/hirsch/Documents/projects/brainSegmentation/DeepPriors/CV_folds/aphasic_stroke/all_TPMs.txt').readlines()

bad = []
for i in range(len(label)):
  gt_img_shape = nib.load(label[i][:-1]).shape  
  mri_img_shape = nib.load(mri[i][:-1]).shape  
  tpm_img_shape = nib.load(tpm[i][:-1]).shape  
  if gt_img_shape != mri_img_shape:
    bad.append(label[i])
  elif mri_img_shape != tpm_img_shape[:-1]:
    bad.append(label[i])
