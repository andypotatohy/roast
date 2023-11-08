#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon May 28 15:35:32 2018

@author: lukas
"""

import nibabel as nib
import numpy as np

nii = nib.load('/home/lukas/Documents/projects/strokeHeads/processed/logTPM/TPM_padded.nii')
d = np.array(nii.get_data(),dtype='float32')
d.dtype
np.min(d[d>0])

d[d<1e-9] = 1e-9  # Values smaller than e-10 start generating too much floating point artifacts.
np.sum(d==0)
np.min(d)
dlog = np.log(d)

img = nib.Nifti1Image(dlog, nii.affine)
nib.save(img, '/home/lukas/Documents/projects/strokeHeads/processed/logTPM/logTPM.nii')
