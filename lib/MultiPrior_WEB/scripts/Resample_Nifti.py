# -*- coding: utf-8 -*-
"""
Created on Tue Aug 27 14:35:08 2019

@author: hirsch
"""

import nibabel as nib
import os
from skimage.transform import resize
import numpy as np
import matplotlib.pyplot as plt
import nilearn


TARGET_AFFINE = np.diag((1,1,1))

FILES_PATH = '/home/hirsch/Documents/projects/test-retest_Maclaren/testretest_data/test-retest_data/brain_data/Subject_2/'

#####################################################################################################

myfiles = os.listdir(FILES_PATH)

myfiles = [x for x in myfiles if x.endswith('.nii')]

if not os.path.exists(FILES_PATH + 'isotropic/'):
  os.mkdir(FILES_PATH + 'isotropic/')

myfiles.sort()

for NIFTI in myfiles:
  print(NIFTI)
  nii = nib.load(FILES_PATH + NIFTI)
  img_res = nilearn.image.resample_img(nii, target_affine=TARGET_AFFINE, interpolation='continuous')
  nib.save(img_res, FILES_PATH + 'isotropic/' + NIFTI.split('.nii')[0] + '_iso.nii.gz')
  
