# -*- coding: utf-8 -*-
"""
Created on Mon Feb 11 17:31:52 2019

@author: hirsch
"""

import nibabel as nib
import numpy as np
import os

TEXT_FILE_MRIs = '/home/hirsch/Documents/projects/stroke_heads/MODEL_INPUT/MRI/ALL_MRI.txt'

STAND_ADDRESS = '/home/hirsch/Documents/projects/stroke_heads/MODEL_INPUT/MRI/STAND/'


def rescale(img, min_new = 0, max_new = 255):
    " Rescale image. Default is int8 [0-255] "
    return ((img - img.min()) * (float(max_new - min_new) / float(img.max() - img.min()))) + min_new
    
def update_sum_N(nii, m, N):
    img = nib.load(nii)
    data = img.get_data()
    m = m + np.sum(data)
    N = N + reduce(lambda x,y: x*y, data.shape)
    return m, N    
    
def update_square_diffs(nii, mean, diffs):
    img = nib.load(nii)
    data = img.get_data()
    difs = np.array((data-mean)**2, dtype='uint64')
    diffs = diffs + np.sum(difs)
    return diffs
    
def get_overall_mean_std_training_set(MRIs):
    m = 0
    N = 0
    fi = open(MRIs)
    lines = fi.readlines()
    for nii in lines:
        m, N = update_sum_N(nii[:-1], m, N)
    mean = float(m)/N
    diffs = np.array(0, dtype='uint64')
    for nii in lines:
        diffs = update_square_diffs(nii[:-1], mean, diffs)
    s = np.sqrt(diffs/N)
    return mean, s

def normalizeMRI(data, mean, std):
    if (mean == 0) and (std == 1):
        mean = np.mean(data)
        std = np.std(data)
    data1 = (data - mean)/std
    return(data1)
    
def standardize_MRIs(MRIs, OUT_FILE, mean=None, std=None):
    if (mean == None) & (std == None):
        mean, std = get_overall_mean_std_training_set(MRIs)
    mris = open(MRIs)
    lines = mris.readlines()
    mris.close()
    fi = open(OUT_FILE + 'stand_padded_rescaled_coreg_MRIs.txt', 'a')
    for nii in lines:
        nii = nii[0:-1]
        img = nib.load(nii)
        data = img.get_data()
        aff = img.affine
        data = rescale(data)
        data_norm = normalizeMRI(data, mean, std)
        img_out = nib.Nifti1Image(data_norm, aff)
        img_out_path = OUT_FILE + nii.split('.')[0].split('/')[-1] + '.nii'
        nib.save(img_out, img_out_path)   
        fi.write("{}\n".format(img_out))
    fi.close()
    fi = open(OUT_FILE + 'MEAN_STD.txt', 'a')
    fi.write("{}".format([mean, std]))
    fi.close()



#%%  NORMALIZATION To Standard Normal Distribution

print('Starting normalization.')

if not os.path.exists(STAND_ADDRESS):
    os.mkdir(STAND_ADDRESS)
    
standardize_MRIs(TEXT_FILE_MRIs, STAND_ADDRESS, mean=0, std=1)
print('End standardization.')
