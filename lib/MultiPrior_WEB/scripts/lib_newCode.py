#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 12 11:35:14 2017

@author: lukas
"""
import os
from shutil import copy
import sys
import nibabel as nib
import numpy as np
np.set_printoptions(precision=3)
import time
import random
from scipy.interpolate import interp1d
import json 
import matplotlib.pyplot as plt
#from random import sample
import keras 
import keras.backend as K
from keras.utils import to_categorical
#from keras.callbacks import ModelCheckpoint
import pandas as pd
from keras.models import load_model   
from skimage.transform import resize
from keras.utils import multi_gpu_model
from numpy.random import seed
from tensorflow import set_random_seed

seed(1)
set_random_seed(2)
        
import multiprocessing        
from multiprocessing import Pool, Process
#from multiprocessing.dummy import Pool as ThreadPool        
        
#import tensorflow as tf
#config = tf.ConfigProto()
#
#config = tf.ConfigProto(device_count = {'GPU': 1})
#config.gpu_options.allow_growth = True
#
#tf.keras.backend.set_session(tf.Session(config=config))

from keras_radam import RAdam

################################## AUXILIARY FUNCTIONS #####################################

  
def flip_random(patches, labels, TPM_patches, coords, proportion_to_flip=0.5):

  # Why only flipping malignants???
  #malignant_indexes = np.argwhere(np.sum(np.sum(labels, axis=-1), axis=-1) > 0)[:,0]
  #indx_toflip = np.random.choice(malignant_indexes, int(len(malignant_indexes)*proportion_to_flip), replace=False)
  
  indx_toflip = np.random.choice(range(patches.shape[0]), int(patches.shape[0]*proportion_to_flip), replace=False)
  axis = np.random.choice(range(0,3),size=len(indx_toflip))
  for i in range(len(indx_toflip)):
    if axis[i] == 0:
      # SAGITTAL FLIP (no flipping of labels, TPM and coords as they are 2D)
      for ch in range(patches.shape[-1]):
          patches[indx_toflip[i],:,:,:,ch] = patches[indx_toflip[i],::-1,:,:,ch]
    elif axis[i] == 1:
      # AXIAL? FLIP
      for ch in range(patches.shape[-1]):
          patches[indx_toflip[i],:,:,:,ch] = patches[indx_toflip[i],:,::-1,:,ch]
      labels[indx_toflip[i],0,:] = np.flip(labels[indx_toflip[i],0,:],0)
      if len(TPM_patches) != 0:
          TPM_patches[indx_toflip[i],0,:] = np.flip(TPM_patches[indx_toflip[i],0,:],0)
      if len(coords) != 0:
          for ch in range(coords.shape[-1]):
              coords[indx_toflip[i],0,:,:,ch] = np.flip(coords[indx_toflip[i],0,:,:,ch],0)
    elif axis[i] == 2:
      # CORONAL? FLIP
      for ch in range(patches.shape[-1]):
          patches[indx_toflip[i],:,:,:,ch] = patches[indx_toflip[i],:,:,::-1,ch]
      labels[indx_toflip[i],0,:] = np.flip(labels[indx_toflip[i],0,:],1)        
      if len(TPM_patches) != 0:
          TPM_patches[indx_toflip[i],0,:] = np.flip(TPM_patches[indx_toflip[i],0,:],1)        
      if len(coords) != 0:
          for ch in range(coords.shape[-1]):
              coords[indx_toflip[i],0,:,:,ch] = np.flip(coords[indx_toflip[i],0,:,:,ch],1)          
              
  return patches, labels , TPM_patches
  
def subtract_MRI(t1post,t1pre):
    t1post_flat = t1post.flatten()
    t1pre_flat = t1pre.flatten()
    w,b = np.polyfit(x=t1post_flat, y=t1pre_flat, deg=1)
    t1post_flat_transformed = w*t1post_flat + b      
    sub = t1post_flat_transformed - t1pre_flat
    sub = sub.reshape(t1post.shape)  
    return sub

def percentile95_normalizeMRI(data, p95=0):
    if p95 == 0:
        p95 = np.percentile(data,95)
    data1 = data/p95
    return(data1)


def getVarFromFile(filename):
    import imp
    print('import using {}'.format(filename))
    f = open(filename)
    global cfg
    cfg = imp.load_source('cfg', '', f)
    f.close()


def file_len(fname):
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i + 1

def normalizeMRI(data, mean=0, std=1):
    if (mean == 0) and (std == 1):
        mean = np.mean(data)
        std = np.std(data)
    data1 = (data - mean)/std
    return(data1)


def classesInSample(minibatch_labels, output_classes):
	label_numbers = []
	print(minibatch_labels.shape)	
	minibatch_labels = np.argmax(minibatch_labels, axis=3)
	for c in range(output_classes):
		label_numbers.append(np.sum(minibatch_labels == c))
	#return label_numbers
	return np.sum(minibatch_labels, axis=-1)
 
################################## SAMPLING FUNCTIONS #####################################


def generateRandomIndexesSubjects(n_subjects, total_subjects):
    indexSubjects = random.sample(xrange(total_subjects), n_subjects)
    return indexSubjects

def getSubjectChannels(subjectIndexes, channel):
    "With the channels (any modality) and the indexes of the selected subjects, return the addresses of the subjects channels"
    fp = open(channel)
    # read file, per subject index extract patches given the indexesPatch
    lines = fp.readlines()
    selectedSubjects = [lines[i][:-1] for i in subjectIndexes]
    fp.close()
    return selectedSubjects


def getSubjectShapes_parallelization(args):

    subjectChannel = args[0]
    resolution = args[1]
    proxy_img = nib.load(subjectChannel)
    res = proxy_img.header['pixdim'][1:4]
    shape = proxy_img.shape
    
    if resolution == 'high':
    
        if res[1] > 0.6:    
            target_res = [res[0],res[1]/2.,res[2]/2.]
            out_shape = np.floor([float(s)*r1/r2 for s,r1,r2 in zip(shape, res, target_res)])
        else:
            out_shape = shape
            
    elif resolution == 'low':
    
        if res[1] < 0.6:    
            target_res = [res[0],res[1]*2.,res[2]*2.]
            out_shape = np.floor([float(s)*r1/r2 for s,r1,r2 in zip(shape, res, target_res)])
        else:
            out_shape = shape
            
            
    return out_shape    

def getSubjectShapes(subjectIndexes, n_patches, channelList):
    # Need to open every nifty file and get the shapes
    fp = open(channelList)
    # read file, per subject index extract patches given the indexesPatch
    lines = fp.readlines()
    selectedSubjects = [lines[i] for i in subjectIndexes]
    fp.close()
    shapes = []
    # Get shapes of all subjects to sample from. Can be a separate function (cause apparently I am needing this everywhere)
    for subjectChannel in selectedSubjects:
        subjectChannel = str(subjectChannel)[:-1]
        proxy_img = nib.load(subjectChannel)
        res = proxy_img.header['pixdim'][1:4]
        shape = proxy_img.shape
        if res[1] > 0.6:    
            target_res = [res[0],res[1]/2.,res[2]/2.]
            out_shape = np.floor([float(s)*r1/r2 for s,r1,r2 in zip(shape, res, target_res)])
        else:
            out_shape = shape
        shapes.append(out_shape)
    return shapes      



def generateVoxelIndexes_wrapper_parallelization(DATA_INPUT):
    subjectIndexes = DATA_INPUT[0]
    target_shape = DATA_INPUT[1]
    patches_per_subject = DATA_INPUT[2]
    channels = DATA_INPUT[3]
    channel_mri = DATA_INPUT[4]     
    dpatch = DATA_INPUT[5]  
    n_patches = DATA_INPUT[6]    
    samplingMethod = DATA_INPUT[7] 
    output_classes = DATA_INPUT[8]  
    percentile_voxel_intensity_sample_benigns = DATA_INPUT[9]
    CV_FOLDS_ARRAYS_PATH = DATA_INPUT[10]
    return generateVoxelIndexes_parallel(subjectIndexes,CV_FOLDS_ARRAYS_PATH, target_shape, patches_per_subject, 
                                         dpatch, n_patches, channels, channel_mri, samplingMethod, output_classes, 
                                         percentile_voxel_intensity_sample_benigns, allForegroundVoxels = "", verbose=False)


def truncated_normal(mean, stddev, minval, maxval):
    return np.clip(np.random.normal(mean, stddev), minval, maxval)

def generateVoxelIndexes_parallel(subjectIndexes,CV_FOLDS_ARRAYS_PATH, target_shape, patches_per_subject, dpatch, n_patches, channels, 
                                  channel_mri, samplingMethod, output_classes, percentile_voxel_intensity_sample_benigns, 
                                  allForegroundVoxels = "", verbose=False):
    allVoxelIndexes = {} #{a:None for a in subjectIndexes} 
    #--------------------------------------------------------------------------------------------------------------------------------------------------
    if samplingMethod == 0:
        "Sample voxels from random locations in a SLICE. This allows sampling from any class present in the SLICE."
        scanVoxels = []

        # only one 2D slice is labeled with the breast mask, so on the sagittal dimension there is no freedom of choice.
        
        if 'BENIGN' not in channels:      
            data_label = nib.load(channels).get_data()
            if np.sum(data_label) == 0:
                print('Found malignant scan with empty segmentation: {}'.format(channels))
                # Assign slice to pick a border slice..
                mySlice = 3
            else:
                #print('Getting foreground voxels..')
                foreground_voxels = np.argwhere(data_label>0)                
                #mySlice = foreground_voxels[0][0]  
                mySlice = foreground_voxels[:,0]
            #background_voxels = np.argwhere(data_label[mySlice] == 0)
            for _ in range(patches_per_subject):
                x = random.choice(mySlice)#mySlice 
                y = random.choice(xrange(dpatch[1]/2,int(target_shape[1])-(dpatch[1]/2)+1))
                z = random.choice(xrange(dpatch[2]/2,int(target_shape[2])-(dpatch[2]/2)+1))
                scanVoxels.append([x,y,z])
        else:
            for _ in range(patches_per_subject):
                #x = int(truncated_normal(mean=target_shape[0]/2, stddev=target_shape[0]/5, minval=0, maxval=target_shape[0])) # --> Make more probably to pick slices from the MIDDLE. 
                x = random.choice(xrange(dpatch[0]/2,int(target_shape[0])-(dpatch[0]/2)+1))   
                y = random.choice(xrange(dpatch[1]/2,int(target_shape[1])-(dpatch[1]/2)+1))
                z = random.choice(xrange(dpatch[2]/2,int(target_shape[2])-(dpatch[2]/2)+1))
                scanVoxels.append([x,y,z])            
            
        allVoxelIndexes[subjectIndexes] = scanVoxels
    #-------------------------------------------------------------------------------------------------------------
    elif samplingMethod == 2:    
        "When feeding a breast Mask, get voxels that are within the breast mask"
        assert os.path.exists(channels), 'Generating voxel-index for samplig: ERROR: path doesnt exist {}'.format(channels)           
        BM_PATH = '/home/deeperthought/kirby_MSK/BreastMasks/alignedNii-Aug2019/'  # LUKAS. Change. Include from config file
        #print('Getting BreastMask Scans from: {}'.format(BM_PATH))
                
        #Check if previously stored arrays indicating sampling locations:
        exam = channel_mri.split('/')[-2]
        side = channel_mri.split('T1_')[-1][0]
        scan_ID = exam + '_' + side
        if 'BENIGN' in channels:
            benign_scan = True            
            voxel_locations_array = CV_FOLDS_ARRAYS_PATH + scan_ID + '_Background_BreastMask_voxel_locations.npz'
        else:
            benign_scan = False
            voxel_locations_array = CV_FOLDS_ARRAYS_PATH + scan_ID + '_Foreground_voxel_locations.npz'        
        
        if os.path.exists(voxel_locations_array):
            #print('Found previously stored voxel locations..')
            candidate_voxels_for_sampling = np.load(voxel_locations_array, allow_pickle=True)
            candidate_voxels_for_sampling = candidate_voxels_for_sampling[candidate_voxels_for_sampling.keys()[0]]
            scanVoxels = candidate_voxels_for_sampling[random.sample(xrange(0,len(candidate_voxels_for_sampling)), min(len(candidate_voxels_for_sampling),patches_per_subject))]
        else:
            # No previously stored voxel coordinates for candidate sampling ############
            ############################################################################    
            bV = 0
            fg = 0      
            nifti_label = nib.load(channels)
            data_label = nifti_label.get_data()
            if np.sum(data_label) > 0:
                # Target label contains segmentation of a tumor. Scan is malignant.
                fg = getForegroundBackgroundVoxels(nifti_label, data_label, target_shape) # This function returns only foreground voxels based on labels.
                if len(fg) == 0:
                  print('resize function removed the foreground voxels...')
                  print(channels)
                  print('Target shape = {}'.format(target_shape))
                  print('Original shape = {}'.format(data_label.shape))
                  sys.exit(0)
                np.savez_compressed(CV_FOLDS_ARRAYS_PATH + scan_ID + '_Foreground_voxel_locations',fg)
                scanVoxels = fg[random.sample(xrange(0,len(fg)), min(len(fg),patches_per_subject))]
            else:  
                # Only getting non-tumor voxels from benign scans:    
                breastMask = [x for x in os.listdir(BM_PATH) if exam in x and side in x][0]
                bm_nii = nib.load(BM_PATH + breastMask).get_data()
                bV = np.argwhere(bm_nii > 0)
                np.savez_compressed(CV_FOLDS_ARRAYS_PATH + scan_ID + '_Background_BreastMask_voxel_locations', bV)                    
                # Half from the intensity-based sampling:
                scanVoxels = bV[random.sample(xrange(0,len(bV)), min( len(bV), patches_per_subject) )] 
            del fg
            del bV                  
        
    #--------------------------------------------------------------------------------------------------------------------------------------------------
    elif samplingMethod == 1:
        "Only for binary classes. Sample only foreground voxels when present. Sample background voxels only from scans that have NO foreground voxels."
        assert os.path.exists(channels), 'Generating voxel-index for samplig: ERROR: path doesnt exist {}'.format(channels)           
        #Check if previously stored arrays indicating sampling locations:
        exam = channel_mri.split('/')[-2]
        side = channel_mri.split('T1_')[-1][0]
        scan_ID = exam + '_' + side
        if 'BENIGN' in channels:
            benign_scan = True
            voxel_locations_array = CV_FOLDS_ARRAYS_PATH + scan_ID + '_Background_voxel_locations_{}_percentile.npz'.format(percentile_voxel_intensity_sample_benigns)
        else:
            benign_scan = False
            voxel_locations_array = CV_FOLDS_ARRAYS_PATH + scan_ID + '_Foreground_voxel_locations.npz'
            
        if os.path.exists(voxel_locations_array):
            #print('Found previously stored voxel locations..')
            candidate_voxels_for_sampling = np.load(voxel_locations_array, allow_pickle=True)
            candidate_voxels_for_sampling = candidate_voxels_for_sampling[candidate_voxels_for_sampling.keys()[0]]
            if benign_scan:
                #print('Scan benign. Sampling half from voxel-intensity > {}'.format(percentile_voxel_intensity_sample_benigns))
                # Half from the intensity-based sampling:
                scanVoxels = candidate_voxels_for_sampling[random.sample(xrange(0,len(candidate_voxels_for_sampling)), 
                                                                         min( len(candidate_voxels_for_sampling), patches_per_subject/2) )].tolist()   
                # Half from random locations:
                for _ in range(patches_per_subject/2):
                    x = int(truncated_normal(mean=target_shape[0]/2, stddev=target_shape[0]/5, minval=0, maxval=target_shape[0]))
                    y = random.choice(xrange(dpatch[1]/2,int(target_shape[1])-(dpatch[1]/2)+1))
                    z = random.choice(xrange(dpatch[2]/2,int(target_shape[2])-(dpatch[2]/2)+1))
                    scanVoxels.append([x,y,z])      
            else:
                #print('Scan malignant, sampling from labeled region.')
                scanVoxels = candidate_voxels_for_sampling[random.sample(xrange(0,len(candidate_voxels_for_sampling)), min(len(candidate_voxels_for_sampling),patches_per_subject))]
               
                
        else:
            # No previously stored voxel coordinates for candidate sampling ############
            ############################################################################    
            bV = 0
            fg = 0      
            nifti_label = nib.load(channels)
            data_label = nifti_label.get_data()
            if np.sum(data_label) > 0:
                # Target label contains segmentation of a tumor. Scan is malignant.
                fg = getForegroundBackgroundVoxels(nifti_label, data_label, target_shape) # This function returns only foreground voxels based on labels.
                if len(fg) == 0:
                  print('resize function removed the foreground voxels...')
                  print(channels)
                  print('Target shape = {}'.format(target_shape))
                  print('Original shape = {}'.format(data_label.shape))
                  sys.exit(0)
                np.savez_compressed(CV_FOLDS_ARRAYS_PATH + scan_ID + '_Foreground_voxel_locations',fg)
                scanVoxels = fg[random.sample(xrange(0,len(fg)), min(len(fg),patches_per_subject))]
            else:  
                # Only getting non-tumor voxels from benign scans:
                if percentile_voxel_intensity_sample_benigns > 0:
                    # Sample high-intensity voxels but also random from the scan.
                    bV = getBodyVoxels(channel_mri, percentile_voxel_intensity_sample_benigns)
                    np.savez_compressed(CV_FOLDS_ARRAYS_PATH + scan_ID + '_Background_voxel_locations_{}_percentile'.format(percentile_voxel_intensity_sample_benigns),bV)                    
                    # Half from the intensity-based sampling:
                    scanVoxels = bV[random.sample(xrange(0,len(bV)), min( len(bV), patches_per_subject/2) )].tolist()   
                    # Half from random locations:
                    for _ in range(patches_per_subject/2):
                        x = int(truncated_normal(mean=target_shape[0]/2, stddev=target_shape[0]/5, minval=0, maxval=target_shape[0])) 
                        y = random.choice(xrange(dpatch[1]/2,int(target_shape[1])-(dpatch[1]/2)+1))
                        z = random.choice(xrange(dpatch[2]/2,int(target_shape[2])-(dpatch[2]/2)+1))
                        scanVoxels.append([x,y,z])          
                else:
                    scanVoxels = []
                    for _ in range(patches_per_subject):
                        x = int(truncated_normal(mean=target_shape[0]/2, stddev=target_shape[0]/5, minval=0, maxval=target_shape[0]))#random.choice(xrange(dpatch[0]/2,int(target_shape[0])-(dpatch[0]/2)+1))  
                        y = random.choice(xrange(dpatch[1]/2,int(target_shape[1])-(dpatch[1]/2)+1))
                        z = random.choice(xrange(dpatch[2]/2,int(target_shape[2])-(dpatch[2]/2)+1))
                        scanVoxels.append([x,y,z])
            del fg
            del bV   
    
    if samplingMethod == 3:
        "Sample voxels from random locations."
        scanVoxels = []

        # only one 2D slice is labeled with the breast mask, so on the sagittal dimension there is no freedom of choice.
        for _ in range(patches_per_subject):
            #x = int(truncated_normal(mean=target_shape[0]/2, stddev=target_shape[0]/5, minval=0, maxval=target_shape[0])) # --> Make more probably to pick slices from the MIDDLE. 
            x = random.choice(xrange(dpatch[0]/2,int(target_shape[0])-(dpatch[0]/2)+1))   
            y = random.choice(xrange(dpatch[1]/2,int(target_shape[1])-(dpatch[1]/2)+1))
            z = random.choice(xrange(dpatch[2]/2,int(target_shape[2])-(dpatch[2]/2)+1))
            scanVoxels.append([x,y,z])            
            
        allVoxelIndexes[subjectIndexes] = scanVoxels        
    #--------------------------------------------------------------------------------------------------------------------------------------------------        
    allVoxelIndexes[subjectIndexes] = scanVoxels
    del scanVoxels
    return allVoxelIndexes    


def getAllForegroundClassesVoxels(groundTruthChannel, dpatch, output_classes):
    '''Get vector of voxel coordinates for all voxel values for all freground classes'''
    "e.g. groundTruthChannel = '/home/hirsch/Documents/projects/ATLASdataset/native_part2/c0011/c0011s0006t01/c0011s0006t01_LesionSmooth_Binary.nii.gz'"
    "NOTE: img in MRICRON starts at (1,1,1) and this function starts at (0,0,0), so points do not match when comparing in MRICRON. Add 1 to all dimensions to match in mricron. Function works properly though"
    img = nib.load(groundTruthChannel)
    data = np.array(img.dataobj[dpatch[0]/2:img.shape[0]-(dpatch[0]/2)+1, 
                                dpatch[1]/2:img.shape[1]-(dpatch[1]/2)+1, 
                                dpatch[2]/2:img.shape[2]-(dpatch[2]/2)+1],dtype='int16') # Get a cropped image, to avoid CENTRAL foreground voxels that are too near to the border. These will still be included, but not as central voxels. As long as they are in the 9x9x9 volume (-+ 4 voxels from the central, on a segment size of 25x25x25) they will still be included in the training.
    img.uncache()    
    voxels = []
    for c in range(1,output_classes):
        coords = np.argwhere(data==c)
        coords = [sum(x) for x in zip(coords , [x/2 for x in dpatch])]
        voxels.append(coords)
    return voxels  # This is a List! Use totuple() to convert if this makes any trouble
            
def getSubjectsToSample(channelList, subjectIndexes):
    "Actually returns channel of the subjects to sample"
    fp = open(channelList)
    lines = fp.readlines()
    subjects = [lines[i] for i in subjectIndexes]
    fp.close()
    return subjects



def extractLabels_parallelization_wrapper(DATA_INPUT_EXTRACT_LABELS_PATCH):
    subject_label_channel = DATA_INPUT_EXTRACT_LABELS_PATCH[0]
    voxelCoordinates = DATA_INPUT_EXTRACT_LABELS_PATCH[1]
    output_dpatch = DATA_INPUT_EXTRACT_LABELS_PATCH[2]
    output_shape = DATA_INPUT_EXTRACT_LABELS_PATCH[3]    
    return extractLabels_parallelization(subject_label_channel, voxelCoordinates, output_dpatch, output_shape)

def extractLabels_parallelization(subject_label_channel, voxelCoordinates, output_dpatch, output_shape):
    if len(voxelCoordinates) == 0:
      print('Within extractLabels_parallelization: \nERROR: len(voxelCoordinates) == 0!, subject_label_channel = {}'.format(subject_label_channel))
      sys.exit(0)
    labels = []       
    subject = str(subject_label_channel)[:-1]
    proxy_label = nib.load(subject)
    label_data = np.array(proxy_label.get_data(),dtype='int8')

    label_data = resize(label_data, order=0, output_shape=output_shape, preserve_range=True, anti_aliasing=True, mode='reflect') 
    #label_data[label_data > 0] = 1
    #DEBUG
    #out = nib.Nifti1Image(label_data, np.diag([1,1,1,0]))
    #nib.save(out, '/home/deeperthought/Projects/MultiPriors_MSKCC/DEBUG/' + 'Label_' + subject.split('MSKCC')[-1].replace('/','_label_'))     
     
    label_padded = np.pad(label_data,((0,60),(0,100),(0,100)),'constant')  # need to pad for segmentation with huge patches that go outside (only the end - ascending coordinates) boundaries. Scale stays the same, as the origin is not modified. 
    if np.sum(label_data) == 0:
      for j in range(len(voxelCoordinates)):
        labels.append(np.zeros((output_dpatch[0],output_dpatch[1],output_dpatch[2]),dtype='int8'))
    else:
      for j in range(len(voxelCoordinates)):
        D1,D2,D3 = voxelCoordinates[j]
        labels.append(label_padded[D1-output_dpatch[0]/2:D1+(output_dpatch[0]/2)+output_dpatch[0]%2,
                                   D2-output_dpatch[1]/2:D2+(output_dpatch[1]/2)+output_dpatch[1]%2,
                                   D3-output_dpatch[2]/2:D3+(output_dpatch[2]/2)+output_dpatch[2]%2])
    proxy_label.uncache()
    del label_data
    return labels


def extractLabels(groundTruthChannel_list, subjectIndexes, voxelCoordinates, output_dpatch, shapes):
    #print('extracting labels from ' + str(len(subjectIndexes))+ ' subjects.')    
    subjects = getSubjectsToSample(groundTruthChannel_list,subjectIndexes)
    labels = []       
    for i in range(len(subjects)):
        subject = str(subjects[i])[:-1]
        #print('extracting labels from subject index [{}] with path : {}'.format(subjectIndexes[i],subject))
        proxy_label = nib.load(subject)
        label_data = np.array(proxy_label.get_data(),dtype='int8')
        # WHERE AM I RESIZING THE SEGMENTATION LABEL???        
        label_data = resize(label_data, shapes[i], order=0, preserve_range=True, anti_aliasing=True)        
        #DEBUG
        #out = nib.Nifti1Image(label_data, np.diag([1,1,1,0]))
        #nib.save(out, '/home/deeperthought/Projects/MultiPriors_MSKCC/DEBUG/' + 'Label_' + subject.split('MSKCC')[-1].replace('/','_label_'))
        label_padded = np.pad(label_data,((0,60),(0,100),(0,100)),'constant')  # need to pad for segmentation with huge patches that go outside (only the end - ascending coordinates) boundaries. Scale stays the same, as the origin is not modified. 
        if np.sum(label_data) == 0:
          for j in range(len(voxelCoordinates[i])):
            labels.append(np.zeros((output_dpatch[0],output_dpatch[1],output_dpatch[2]),dtype='int8'))
        else:
          for j in range(len(voxelCoordinates[i])):
            D1,D2,D3 = voxelCoordinates[i][j]
            #print('Extracting labels from \n subject {} with shape {} and coords {},{},{}'.format(subjects[i], label_data.shape ,D1,D2,D3))
            labels.append(label_padded[D1-output_dpatch[0]/2:D1+(output_dpatch[0]/2)+output_dpatch[0]%2,
                                       D2-output_dpatch[1]/2:D2+(output_dpatch[1]/2)+output_dpatch[1]%2,
                                       D3-output_dpatch[2]/2:D3+(output_dpatch[2]/2)+output_dpatch[2]%2])
            #if len(labels[-1])==0:
            #  labels[-1] = np.zeros((9,9),dtype='int8')
        proxy_label.uncache()
        del label_data
    return labels

#shapes = [shape]
#voxelCoordinates = [voxelCoordinates]
   
def extractCoordinates(shapes, voxelCoordinates, output_dpatch):
    """ Given a list of voxel coordinates, it returns the absolute location coordinates for a given patch size (output 1x9x9) """
    #print('extracting coordinates from ' + str(len(subjectIndexes))+ ' subjects.')
    #subjects = getSubjectsToSample(channel, subjectIndexes)
    
    all_coordinates = []
    for i in xrange(len(shapes)):
        #subject = str(subjects[i])[:-1]
        #img = nib.load(subject)
        img_shape = shapes[i]
        for j in xrange(len(voxelCoordinates[i])):     
            D1,D2,D3 = voxelCoordinates[i][j]
            #all_coordinates.append(get_Coordinates_from_target_patch(img.shape,D1,D2,D3))                 
            all_coordinates.append(get_Coordinates_from_target_patch(img_shape,D1,D2,D3, output_dpatch))                    

        #img.uncache()
    return np.array(all_coordinates)    


def get_Coordinates_from_target_patch(img_shape,D1,D2,D3, output_dpatch) :

    x_ = range(D1-(output_dpatch[0]//2),D1+((output_dpatch[0]//2)+output_dpatch[0]%2))
    y_ = range(D2-(output_dpatch[1]//2),D2+((output_dpatch[1]//2)+output_dpatch[1]%2))
    z_ = range(D3-(output_dpatch[2]//2),D3+((output_dpatch[2]//2)+output_dpatch[2]%2))
    
    x_norm = np.array(x_)/float(img_shape[0])  
    y_norm = np.array(y_)/float(img_shape[1])  
    z_norm = np.array(z_)/float(img_shape[2])  
    
    x, y, z = np.meshgrid(x_norm, y_norm, z_norm, indexing='ij')    
    coords = np.stack([x,y,z], axis=-1)
    return coords
    
       
def get_patches_per_subject( n_patches, n_subjects):
    patches_per_subject = [n_patches/n_subjects]*n_subjects
    randomAdd = random.sample(range(0,len(patches_per_subject)),k=n_patches%n_subjects)
    randomAdd.sort()
    for index in randomAdd:
        patches_per_subject[index] = patches_per_subject[index] + 1
    return patches_per_subject



def extract_TPM_patches_parallelization_wrapper(TPM_INPUT_DATA):
    TPM_channel = TPM_INPUT_DATA[0]
    subjectIndexes = TPM_INPUT_DATA[1]
    voxelCoordinates_subject = TPM_INPUT_DATA[2]
    output_dpatch = TPM_INPUT_DATA[3]
    output_shape = TPM_INPUT_DATA[4]
    return extract_TPM_patches_parallelization(TPM_channel, subjectIndexes, voxelCoordinates_subject, output_dpatch, output_shape)

def extract_TPM_patches_parallelization(TPM_channel, subjectIndexes, voxelCoordinates_subject, output_dpatch, output_shape):    
    vol = np.zeros((len(voxelCoordinates_subject),output_dpatch[0],output_dpatch[1],output_dpatch[2]),dtype='float32')
    proxy_label = nib.load(TPM_channel)
    TPM_data = np.array(proxy_label.get_data())#,dtype='float32')  
    #print('Resizing TPM from {} to {}'.format(TPM_data.shape, output_shape))
    TPM_data = resize(TPM_data, output_shape, order=1, preserve_range=True, anti_aliasing=True, mode='reflect') 
    
    #DEBUG
    #out = nib.Nifti1Image(TPM_data, np.diag([1,1,1,0]))
    #nib.save(out, '/home/deeperthought/Projects/MultiPriors_MSKCC/DEBUG/' + str(subjectIndexes) + '_TPM.nii' )          
    
    padding_border = 100    # why so large?              
    TPM_data = np.pad(TPM_data, padding_border,'minimum')
    for j in range(len(voxelCoordinates_subject)):     
        D1,D2,D3 = voxelCoordinates_subject[j]
        # Scale these coordinates to fit the shape of the Breast Tumor TPM! 
        D1 = D1 + padding_border
        D2 = D2 + padding_border
        D3 = D3 + padding_border        

        try:
          vol[j] = TPM_data[D1-output_dpatch[0]/2:D1+(output_dpatch[0]/2)+output_dpatch[0]%2,
                                       D2-output_dpatch[1]/2:D2+(output_dpatch[1]/2)+output_dpatch[1]%2,
                                       D3-output_dpatch[2]/2:D3+(output_dpatch[2]/2)+output_dpatch[2]%2]
        except:
          print('FAILED TO EXTRACT TPM PATCH')
          print('Coordinates: {}'.format([D1,D2,D3]))
          print('From subject index {}'.format(subjectIndexes))
          print('TPM shape after padding: {}'.format(TPM_data.shape))
          
    proxy_label.uncache()
    del proxy_label
    del TPM_data
    return vol

    
def extract_TPM_patches(TPM_channel, subjectIndexes, voxelCoordinates, output_dpatch, shapes):
    print('extracting TPM patches from ' + str(len(subjectIndexes))+ ' subjects.')          
    n_patches = 0
    k = 0
    for i in range(len(voxelCoordinates)):
        n_patches += len(voxelCoordinates[i])
    vol = np.zeros((n_patches,output_dpatch[0],output_dpatch[1],output_dpatch[2]),dtype='float32')             # CHANGE THIS TO CHANNELS OF TPM!!
    print('TPM patches shape: {}'.format(vol.shape))
    for i in range(len(subjectIndexes)):            
        proxy_label = nib.load(TPM_channel)
        TPM_data = np.array(proxy_label.get_data(),dtype='float32')  
        # 'align' the TPM to the input image
        TPM_data = resize(TPM_data, shapes[i], order=1, preserve_range=True, anti_aliasing=True)
        
        #DEBUG
        #out = nib.Nifti1Image(TPM_data, np.diag([1,1,1,0]))
        #nib.save(out, '/home/deeperthought/Projects/MultiPriors_MSKCC/DEBUG/' + str(subjectIndexes) + '_TPM.nii' )         
                
        
        padding_border = 100                  
        label_data_padded = np.pad(TPM_data[:,:,:],
                                   ((padding_border,padding_border), 
                                    (padding_border,padding_border), 
                                    (padding_border,padding_border)),
                                   'minimum')
        for j in range(len(voxelCoordinates[i])):     
            D1,D2,D3 = voxelCoordinates[i][j]
            # Scale these coordinates to fit the shape of the Breast Tumor TPM! 
            D1 = D1 + padding_border
            D2 = D2 + padding_border
            D3 = D3 + padding_border        

            try:
              vol[k,:,:,:] = label_data_padded[D1-output_dpatch[0]/2:D1+(output_dpatch[0]/2)+output_dpatch[0]%2,
              D2-output_dpatch[1]/2:D2+(output_dpatch[1]/2)+output_dpatch[1]%2,
              D3-output_dpatch[2]/2:D3+(output_dpatch[2]/2)+output_dpatch[2]%2]
            except:
              print('FAILED TO EXTRACT TPM PATCH')
              print('Coordinates: {}'.format([D1,D2,D3]))
              print('From subject index {}'.format(subjectIndexes[i]))
              print('TPM shape after padding: {}'.format(label_data_padded.shape))
                                               
            k=k+1
        proxy_label.uncache()
        
    return vol

# Need one function that takes in a list of subject indexes and outputs patch volumes.
# Need one wrapper function that extracts all the arguments needed out of a single list/matrix.
#
#for i in range(len(DATA_INPUT_EXTRACT_IMAGE_PATCH)):
#    print(i)
#    DATA_INPUT = DATA_INPUT_EXTRACT_IMAGE_PATCH[i]
#    pathchs = extractImagePatch_parallelization(MRI_PATH, channel, subjectIndex, subject_channel_voxelCoordinates, output_shape, dpatch, 
#                                           intensity_normalization_method)
#    pathchs.shape

def extractImagePatch_parallelization_wrapper(DATA_INPUT):
  # This has to have all voxel coordinates for each subject... 
  subjectIndex = DATA_INPUT[0]
  channel = DATA_INPUT[1]
  dpatch = DATA_INPUT[2]
  subject_channel_voxelCoordinates = DATA_INPUT[3]
  output_shape = DATA_INPUT[4]
  intensity_normalization_method = DATA_INPUT[5]  #integer 1,2,3
  MRI_PATH = DATA_INPUT[6]   # string
  return extractImagePatch_parallelization(MRI_PATH, channel, subjectIndex, subject_channel_voxelCoordinates, output_shape, dpatch, 
                                           intensity_normalization_method)


def normalize_using_LUT(img, LUT_cluster):
    img[img < 0] = 0 
    data = img.reshape(np.prod(img.shape))
    img_new = LUT_cluster(data)        
    nifti_new_normalized = img_new.reshape(img.shape)
    return nifti_new_normalized


def extractImagePatch_parallelization(MRI_PATH, channel, subjectIndex, subject_channel_voxelCoordinates, output_shape, dpatch, 
                                      intensity_normalization_method, preprocess_image_data=True,fullSegmentationPhase=False):   
    subject_channel = getSubjectsToSample(channel, [subjectIndex])
#    if 'slope' in subject_channel[0]:
#        t1post_channel = subject_channel[0][:-1].replace('slope1','02_01').replace('slope2','02_01')
#        p95 = np.percentile(nib.load(t1post_channel).get_data(),95)
#    else:
    p95 = 0
    n_patches = len(subject_channel_voxelCoordinates)
    subject = str(subject_channel[0])[:-1]   
    
    scanID = subject_channel[0].split('/')[-2] + '_' + subject_channel[0].split('T1_')[-1][0]   
    exam = scanID[:-2]
    side = scanID[-1]    
    scans = os.listdir(MRI_PATH + exam)
    
    T1_pre_nii_path = [MRI_PATH + exam + '/' + nifti for nifti in scans if 'T1_' in nifti and nifti.endswith('01_01.nii.gz') and nifti[3] == side][0]
    T1_post_nii_path = [MRI_PATH + exam + '/' + nifti for nifti in scans if 'T1_' in nifti and nifti.endswith('02_01.nii.gz') and nifti[3] == side][0]    
    # NOT LOADING SLOPE1 FOR NOW. WHEN USING OLD MODEL, ADD OPTION TO ADD SLOPE1 AND SLOPE2. ITS JUST 2 MORE LINES HERE. 
    # ALSO INSTEAD OF DOING ALL PREPROCESSING PER IMAGE, DO A LOOP OVER ALL IMAGES, JUST LIKE BELOW WHEN EXTRACTING THE PATCHES. JUST MOVE THAT LOOP UPWARDS.    
    
    proxy_img_T1pre = nib.load(T1_pre_nii_path)            
    img_data_T1pre = np.array(proxy_img_T1pre.get_data(),dtype='float32')
    proxy_img_T1post = nib.load(T1_post_nii_path)            
    img_data_T1post = np.array(proxy_img_T1post.get_data(),dtype='float32')

    img_data_T1pre[img_data_T1pre < 0] = 0
    img_data_T1post[img_data_T1post < 0] = 0
    
    if preprocess_image_data:   
    
        
      if np.array(img_data_T1pre.shape != output_shape).any():
        #print('Resizing training data: \nInput_shape = {}, \nOutput_shape = {}. \nSubject = {}'.format(img_data.shape, output_shape, subject))
        img_data_T1pre = resize(img_data_T1pre, output_shape=output_shape, preserve_range=True, anti_aliasing=True, mode='reflect')
        img_data_T1post = resize(img_data_T1post, output_shape=output_shape, preserve_range=True, anti_aliasing=True, mode='reflect')
        
        
      if np.any(np.isnan(img_data_T1pre)):
        print('Nans found in scan {}'.format(subject))
        print('Nans replace by value: {}'.format(np.nanmin(img_data_T1pre)))
        img_data_T1pre[np.isnan(img_data_T1pre)] = np.nanmin(img_data_T1pre)
      if np.any(np.isnan(img_data_T1post)):
        print('Nans found in scan {}'.format(subject))
        print('Nans replace by value: {}'.format(np.nanmin(img_data_T1post)))
        img_data_T1post[np.isnan(img_data_T1post)] = np.nanmin(img_data_T1post)
      
      if intensity_normalization_method == 1:
        #print('Normalizing intensities using histogram matching')
        # ALL THIS HARD-CODED PATHS NEED TO GO
        keys = pd.read_csv('/home/deeperthought/Projects/Intensity_Normalization_stuff/scanID_cluster_keys.csv')          
        if scanID in list(keys['scanID']):
            cluster = keys.loc[keys['scanID'] == scanID, 'cluster_ID'].values[0]            
        else:
            print('{} : scan not in dataframe!'.format(scanID))
            return np.random.normal(size=(n_patches,dpatch[0],dpatch[1],dpatch[2], 3))
            #sys.exit(0)
        
        #print('Loading LUT..')
        LUT_PATH = '/home/deeperthought/Projects/Intensity_Normalization_stuff/LUT/'
        LUT_arr = json.load(open(LUT_PATH + 'LUT_cluster_{}.json'.format(cluster), 'r'))
        x = [tup[0] for tup in LUT_arr]
        y = [tup[1] for tup in LUT_arr]
        
        # Extrapolation. I need to re-make the LUTs. With all 70k scans!!
        x[-1] = 50*x[-1]
        y[-1] = 50*y[-1]
        
        LUT_cluster = interp1d(x,y)
        #print('Normalizing scan..')  
        img_data_T1pre = normalize_using_LUT(img_data_T1pre, LUT_cluster)
        img_data_T1post = normalize_using_LUT(img_data_T1post, LUT_cluster)      
          
      elif intensity_normalization_method == 2:
        img_data_T1pre = percentile95_normalizeMRI(img_data_T1pre, p95)
        img_data_T1post = percentile95_normalizeMRI(img_data_T1post, p95)

      elif intensity_normalization_method == 3:
        img_data_T1pre = normalizeMRI(img_data_T1pre)      
        img_data_T1post = normalizeMRI(img_data_T1post)    
        
      if not np.isfinite(img_data_T1pre).all():
        print('Normalization: Nans found in scan {}'.format(subject))
        print('Nans replace by value: {}'.format(np.nanmin(img_data_T1pre)))
        img_data_T1pre[ ~ np.isfinite(img_data_T1pre)] = np.nanmin(img_data_T1pre)
      if not np.isfinite(img_data_T1post).all():
        print('Normalization: Nans found in scan {}'.format(subject))
        print('Nans replace by value: {}'.format(np.nanmin(img_data_T1post)))
        img_data_T1post[ ~ np.isfinite(img_data_T1post)] = np.nanmin(img_data_T1post)
        
    
    slope1 = img_data_T1post - img_data_T1pre
    #nii_out = nib.Nifti1Image(slope1, nii.affine)
    #nib.save(nii_out, MRI_PATH + exam + '/T1_right_slope1_Normalized.nii')                 
        
    subject_all_channels = [img_data_T1pre, img_data_T1post, slope1]

    vol = np.zeros((n_patches,dpatch[0],dpatch[1],dpatch[2], len(subject_all_channels)),dtype='float32') 


    for mri_modality in range(len(subject_all_channels)):
        img_data = subject_all_channels[mri_modality]
        if fullSegmentationPhase:      
            if np.max(dpatch) > 200:  # This is the case with the full-image U-Net_v0. If we pad too big, this takes a lot of time and unnecessary resources.
                padding_border = 10#np.max(dpatch)#np.max(dpatch)/2 + 10#550
            else:
                padding_border = np.max(dpatch)
        else:
            padding_border = np.max(dpatch)/2 + 10
        # Padding needs to be larger than dpatch/2. During training all patches are centered within the image so here its enough.
        # But during testing, we need to sample over all center patches, which means that we can go outside the original image boundaries.
        # Example: image size = 7, center patch size = 3
        # Image : [0123456] center-patches: [012][345][678] . Last patch was centered on 7, just to capture the 6 on the border. 
        #print('Padding image..')
        img_data_padded = np.pad(img_data, padding_border,'reflect')    
        
        for j in range(n_patches):      
            D1,D2,D3 = subject_channel_voxelCoordinates[j]           
            D1 = D1 + padding_border
            D2 = D2 + padding_border
            D3 = D3 + padding_border
            try:
              vol[j,:,:,:,mri_modality] = img_data_padded[D1-(dpatch[0]/2):D1+(dpatch[0]/2)+dpatch[0]%2,
                                             D2-(dpatch[1]/2):D2+(dpatch[1]/2)+dpatch[1]%2,
                                             D3-(dpatch[2]/2):D3+(dpatch[2]/2)+dpatch[2]%2]
            except:
              print('Failed to extract image data into shape... This is: \n{}, \nimg_data_padded.shape = {}, \nCoords = {}, \nCoords+Padding = {}'.format(subject_channel,img_data_padded.shape, subject_channel_voxelCoordinates[j] , [D1,D2,D3] ))
              sys.exit(0)
              
              
    proxy_img_T1pre.uncache()
    proxy_img_T1post.uncache()
    del img_data
    del img_data_padded
    return vol


#channel = testChannels[0]    
#subjectIndexes = subjectIndex 
#dpatch = segmentation_dpatch
#voxelCoordinates = [minibatch_voxelCoordinates]

#
#def extractImagePatch(channel, subjectIndexes, patches, voxelCoordinates, dpatch, debug=False, preprocess_image_data=True):
#    subjects = getSubjectsToSample(channel, subjectIndexes)
#    n_patches = 0   
#    # Replace this thing. No need to compute. Have this information in list patches_per_subject!
#    for i in range(len(voxelCoordinates)):
#        n_patches += len(voxelCoordinates[i])
#    #print('Starting extraction of {} patches from {} subjects.'.format(n_patches,len(voxelCoordinates)))
#    vol = np.ones((n_patches,dpatch[0],dpatch[1],dpatch[2]),dtype='float32')
#    k = 0
#    
#    for i in range(len(subjectIndexes)):   
#        #if i%20==0:
#        #  print('{}%'.format(round(i*100./len(voxelCoordinates),2)))
#        subject = str(subjects[i])[:-1]
#        #print('Subject with path: {}'.format(subject))
#        proxy_img = nib.load(subject)            
#        img_data = np.array(proxy_img.get_data(),dtype='float32')
#
#        if preprocess_image_data:
#          # Change resolution:
#          res = proxy_img.header['pixdim'][1:4]
#          shape = img_data.shape
#          if res[1] > 0.6: 
#              'Upsampling image..'
#              target_res = [res[0],res[1]/2.,res[2]/2.]
#              out_shape = np.floor([float(s)*r1/r2 for s,r1,r2 in zip(shape, res, target_res)])
#              img_data = resize(img_data, output_shape=out_shape, preserve_range=True, anti_aliasing=True)
#          #else:
#              #out_shape = shape
#          if np.any(np.isnan(img_data)):
#            print('Nans found in scan {}'.format(subject))
#            print('Nans replace by value: {}'.format(np.nanmin(img_data)))
#            img_data[np.isnan(img_data)] = np.nanmin(img_data)
#          
#          #shape = img_data.shape
#          #print('Shape after normalization: {}'.format(shape))
#          # Standardize image.
#          #for i in range(img_data.shape[0]):
#          #  img_data[i,:,:] = percentile95_normalizeMRI(img_data[i,:,:])
#
#          # NEED TO ADD CONDITION SET FROM CONFIG FILE !! BUG
#          img_data = percentile95_normalizeMRI(img_data)
#          #img_data = normalizeMRI(img_data)
#
#        padding_border = np.max(dpatch)#np.max(dpatch)/2 + 10#550
#        # Padding needs to be larger than dpatch/2. During training all patches are centered within the image so here its enough.
#        # But during testing, we need to sample over all center patches, which means that we can go outside the original image boundaries.
#        # Example: image size = 7, center patch size = 3
#        # Image : [0123456] center-patches: [012][345][678] . Last patch was centered on 7, just to capture the 6 on the border. 
#        img_data_padded = np.pad(img_data, 
#                                 padding_border,
#                                 'reflect')
#        
#        # Loop over voxelCoordinates tuples of subject i
#        for j in range(len(voxelCoordinates[i])):   
#            #print(voxelCoordinates[i][j] )     
#            D1,D2,D3 = voxelCoordinates[i][j]           
#
#            D1 = D1 + padding_border#dpatch[0]/2
#            D2 = D2 + padding_border#dpatch[1]/2
#            D3 = D3 + padding_border#dpatch[2]/2
#
#            vol[k,:,:,:] = img_data_padded[D1-(dpatch[0]/2):D1+(dpatch[0]/2)+dpatch[0]%2,
#                                           D2-(dpatch[1]/2):D2+(dpatch[1]/2)+dpatch[1]%2,
#                                           D3-(dpatch[2]/2):D3+(dpatch[2]/2)+dpatch[2]%2]
#
#            k = k+1  
#        
#        proxy_img.uncache()
#        del img_data
#        if debug: print('extracted [' + str(len(voxelCoordinates[i])) + '] patches from subject ' + str(i) +'/'+ str(len(subjectIndexes)) +  ' with index [' + str(subjectIndexes[i]) + ']')        
#    #print('In this batch found {} Bad Coordinates \n'.format(badCoords))
#    #print('From subject(s): {}'.format(list(set(badCoords_subj))))
#    #raw_input("Press Enter to continue...")
#    return vol



#trainChannels = cfg.trainChannels
#trainLabels = cfg.trainLabels
#TPM_channel = cfg.TPM_channel
#n_subjects = 5#cfg.n_subjects
#n_patches = 100#cfg.n_patches
#dpatch = cfg.dpatch
#output_classes = cfg.output_classes
#samplingMethod = cfg.samplingMethod_train
#use_coordinates = cfg.use_coordinates
#balanced_sample_subjects = cfg.balanced_sample_subjects
#verbose=False
#debug=False
#proportion_malignants_to_sample = cfg.proportion_malignants_to_sample_val
#percentile_voxel_intensity_sample_benigns = cfg.percentile_voxel_intensity_sample_benigns
#intensity_normalization_method = cfg.intensity_normalization_method
#proportion_to_flip = cfg.proportion_to_flip
#procnum = 'TRAINING'
#resolution = cfg.resolution
#model_patch_reduction = cfg.model_patch_reduction
#MRI_PATH = cfg.MRI_PATH

def sampleTrainData_daemon(return_dict, procnum, resolution, trainChannels,CV_FOLDS_ARRAYS_PATH, trainLabels, TPM_channel, n_patches, 
                           n_subjects, dpatch, output_classes, samplingMethod, use_coordinates, proportion_malignants_to_sample, 
                           percentile_voxel_intensity_sample_benigns, data_augmentation, proportion_to_flip, 
                           intensity_normalization_method, MRI_PATH, model_patch_reduction, model_crop, balanced_sample_subjects=True, 
                           verbose=False, debug=False, using_unet=True):
    
    print('Setting number of channels manually to 3, during testing new code (only passing T1pre and T1post, with new normalization method generating slope1 on the fly)')
    num_channels = 3 #len(trainChannels)
    
    output_dpatch = dpatch[0] - model_patch_reduction[0], dpatch[1] - model_patch_reduction[1], dpatch[2] - model_patch_reduction[2]
    patches_per_subject = get_patches_per_subject( n_patches, n_subjects)    
    labelsFile = open(trainLabels).readlines()    
    total_subjects = len(labelsFile)

    if balanced_sample_subjects:     
      proportion_malignants = int(np.ceil(n_subjects*proportion_malignants_to_sample))
      malignant_subjects_index = [labelsFile.index(x) for x in labelsFile if not 'BENIGN' in x]
      benign_subjects_index = list(set(range(total_subjects)) - set(malignant_subjects_index))
      subjectIndexes = random.sample(malignant_subjects_index, min(len(malignant_subjects_index), proportion_malignants))
      print('{} : sampling {} malignants from partition'.format(procnum,len(subjectIndexes)))
      try:
        subjectIndexes.extend(random.sample(benign_subjects_index, n_subjects - len(subjectIndexes)))
      except:
        # if not enough or only malignants in set.
        subjectIndexes.extend(random.sample(malignant_subjects_index, n_subjects - len(subjectIndexes)))
      random.shuffle(subjectIndexes)
    else:
      print('{} : Extracting data from randomly selected subjects.. [breast mask model]'.format(procnum))  
      subjectIndexes = generateRandomIndexesSubjects(n_subjects, total_subjects)  
    
    #------------- Parallelization ----------------
    channel_mri = getSubjectChannels(subjectIndexes, trainChannels[0]) 
    #GET_SHAPES_INPUT = zip(channel_mri, [CV_FOLDS_SHAPES_PATH]*len(channel_mri))
    pool = Pool(min(multiprocessing.cpu_count()/4,19))#mp.cpu_count() -1)
    time1 = time.time()
    shapes = pool.map(getSubjectShapes_parallelization, zip(channel_mri,[resolution]*len(channel_mri)))
    print('{} : Getting scan shapes took {} s'.format(procnum, round(time.time() - time1,2)))

    ############ Generating Voxel Coordinates For training ##############    
    print('{} : ------------ Generating List of Voxel Indexes for sampling  ------------'.format(procnum))
    #------------ Parallelization --------------------
    channels = getSubjectChannels(subjectIndexes, trainLabels)
    channel_mri = getSubjectChannels(subjectIndexes, trainChannels[0]) 
    
    # INPUT SHOULD BE OUTPUT PATCH, NOT DPATCH!!! WE PUT CONSTRAIN ON BORDER OF OUTPUT NOT INPUT!! 
    
    DATA_INPUT_WRAPPER = zip(subjectIndexes, shapes, patches_per_subject, channels, channel_mri, [dpatch] * len(patches_per_subject), 
                             [n_patches]*len(patches_per_subject), [samplingMethod]*len(patches_per_subject), [output_classes]*len(patches_per_subject), 
                             [percentile_voxel_intensity_sample_benigns]*len(patches_per_subject), [CV_FOLDS_ARRAYS_PATH]*len(patches_per_subject))
    DATA_INPUT_WRAPPER = list(DATA_INPUT_WRAPPER)
    #pool = Pool(multiprocessing.cpu_count() - 1)#mp.cpu_count() -1)
    time1 = time.time()
    try:      
      voxelCoordinates_full = pool.map(generateVoxelIndexes_wrapper_parallelization, DATA_INPUT_WRAPPER)
    except IndexError:
      print('...IndexError ')
      print(DATA_INPUT_WRAPPER)
    print('{} : Generating List of Voxel Indexes for sampling took {} s'.format(procnum, round(time.time() - time1,2)))
    #pool.close()
    #pool.join()  
    subjectIndexes_check = []
    voxelCoordinates = []
    for i in range(len(voxelCoordinates_full)):
      subjectIndexes_check.append(voxelCoordinates_full[i].keys()[0])
      voxelCoordinates.extend(voxelCoordinates_full[i].values())
    assert subjectIndexes_check == subjectIndexes, 'Subject Indexes got out of order through multiprocessing...'
    del subjectIndexes_check
    del voxelCoordinates_full
    
    if debug:
        affine = np.diag((1,1,1,0))
        data = np.zeros((65,400,400))
        for coord_set in voxelCoordinates:
            for coord in coord_set:
                D1,D2,D3 = coord
                data[D1-4:D1+5,D2-4:D2+5,D3-4:D3+5] += 1
        img = nib.Nifti1Image(data, affine)
        nib.save(img,'/home/hirsch/Documents/projects/Breast_segmentation/DeepPriors_package/debug_voxel_sampling_location.nii'.format(n_subjects, n_patches))

    real_n_patches = 0
    for i in range(len(voxelCoordinates)):
        if len(voxelCoordinates[i]) == 0:
          print('Empty voxelCoordinates for i = {}, channel_mri[i] = {}, \n\n DATA_INPUT_WRAPPER[i] = {}, \n\n patches_per_subject = {},\n patches_per_subject[i] = {}'.format(i, channel_mri[i], DATA_INPUT_WRAPPER[i], patches_per_subject, patches_per_subject[i]))
          sys.exit(0)
        real_n_patches += len(voxelCoordinates[i])        

    ############## Parallelization Extract Image Patches ####################
    print('{} : ------  Extracting {} image patches from {} subjects, for each of {} channels --------'.format(procnum, real_n_patches,len(voxelCoordinates), len(trainChannels) ))
    patches = np.zeros((real_n_patches,dpatch[0],dpatch[1],dpatch[2],num_channels),dtype='float32')       
    #for i in range(num_channels):
        #print('{} : Extracting image patches from channel: {}'.format(procnum, trainChannels[i]))  


    DATA_INPUT_EXTRACT_IMAGE_PATCH = list(zip(subjectIndexes, [trainChannels[0]]*len(subjectIndexes), [dpatch] * len(subjectIndexes), 
                                         voxelCoordinates, shapes, [intensity_normalization_method]*len(subjectIndexes), 
                                         [MRI_PATH]*len(subjectIndexes)))
    
    #pool = Pool(multiprocessing.cpu_count() -1) #-1 )
    time1 = time.time()
    channel_patches = pool.map(extractImagePatch_parallelization_wrapper, DATA_INPUT_EXTRACT_IMAGE_PATCH)
    print('{} : Extracting image patches took {} s'.format(procnum, round(time.time() - time1,2))) 
    
    start = 0
    for ii in range(len(channel_patches)): # iterates over number of subjects
      patches[start:start+len(channel_patches[ii])] = channel_patches[ii]
      start = start+len(channel_patches[ii])        
    del channel_patches
        
    print('{} : ------  Extracting {} target-label patches from {} subjects --------'.format(procnum, real_n_patches,len(voxelCoordinates) ))

    ############# Parallelization Extract Label Patches ###################### 
    subjects_label_channels = getSubjectsToSample(trainLabels,subjectIndexes)          
    DATA_INPUT_EXTRACT_LABELS_PATCH = zip(subjects_label_channels, voxelCoordinates, [output_dpatch] * len(subjectIndexes), shapes)
    DATA_INPUT_EXTRACT_LABELS_PATCH = list(DATA_INPUT_EXTRACT_LABELS_PATCH)
    #pool = Pool(multiprocessing.cpu_count() -1) #-1 )
    time1 = time.time()
    labels_list_unflattened = pool.map(extractLabels_parallelization_wrapper, DATA_INPUT_EXTRACT_LABELS_PATCH)
    print('{} : Extracting target label patches took {} s'.format(procnum, round(time.time() - time1,2)))   
    #pool.close()
    #pool.join() 
    labels = np.zeros(([real_n_patches] + list(output_dpatch) ), dtype='int8')
    start = 0
    try:
      for ii in range(len(labels_list_unflattened)):
        labels[start:start+len(labels_list_unflattened[ii])] = labels_list_unflattened[ii]
        start = start+len(labels_list_unflattened[ii])        
    except ValueError:
      print('ValueError... ii = {}, start={}, len(labels_list_unflattened[ii]) = {},  subjects_label_channels[ii]: {}'.format(ii, start,
            len(labels_list_unflattened[ii]), subjects_label_channels[ii] ))
    del labels_list_unflattened    

    labels = np.array(labels,dtype='int8')
    labels_list = np.array(labels)
    labels = np.array(to_categorical(labels.astype(int),output_classes),dtype='int8')
    #if(samplingMethod == 2):
    #    patches = patches[0:len(labels)]  # when using equal sampling (samplingMethod 2), because some classes have very few voxels in a head, there are fewer patches as intended. Patches is initialized as the maximamum value, so needs to get cut to match labels.

    ############ Location Coordinates ####################
    if use_coordinates:
      all_coordinates = extractCoordinates(shapes, voxelCoordinates, output_dpatch)
      if debug:
        y_coords = all_coordinates[:,2,:,:]
        plt.imshow(y_coords[1])
        center_y = y_coords[:,5,1]
        plt.hist(center_y,  200)
        plt.xlabel('Normalized Y coordinate')

    else:
      all_coordinates = []
      
    ############ TPM Patches #############################  
    if len(TPM_channel) > 0:
      print('{} : ------  Extracting {} TPM patches from {} subjects --------'.format(procnum, real_n_patches,len(voxelCoordinates) ))      
      TPM_INPUT_DATA = zip([TPM_channel]*len(subjectIndexes), subjectIndexes, voxelCoordinates, [output_dpatch] * len(subjectIndexes), shapes)
      TPM_INPUT_DATA = list(TPM_INPUT_DATA)
      time1 = time.time()
      TPM_patches_unflattened = pool.map(extract_TPM_patches_parallelization_wrapper, TPM_INPUT_DATA)
      print('{} : Extracting TPM patches took {} s'.format(procnum, round(time.time() - time1,2))) 
      TPM_patches = np.zeros(([real_n_patches] + list(output_dpatch)),dtype='float32')
      start = 0
      for ii in range(len(TPM_patches_unflattened)):
        TPM_patches[start:start+len(TPM_patches_unflattened[ii])] = TPM_patches_unflattened[ii]
        start = start+len(TPM_patches_unflattened[ii])        
      del TPM_patches_unflattened              
         
    else:
      TPM_patches = []  
      
    pool.close()
    pool.join() 
      
    if data_augmentation:
      print('{} : Data augmentation: Randomly flipping {}% patches..'.format(procnum, proportion_to_flip*100))
      patches, labels, TPM_patches = flip_random(patches, labels, TPM_patches, all_coordinates, proportion_to_flip)  # WHY THIS DOESNT RETURN ALL_COORDINATES TOO???
       
    if debug:
        patches.shape
        display_index = 5
        #y_coords = all_coordinates[:,2,:,:]
        labels_img = np.array(labels_list,dtype='int8')
        plt.figure(figsize=(12,8))
        plt.subplot(231)
        plt.imshow(patches[display_index,5,:,:,0], cmap='gray')
        plt.title('T1post')
        plt.subplot(232)
        plt.imshow(patches[display_index,5,:,:,1], cmap='gray')
        plt.title('T1pre')
        plt.subplot(233)
        plt.imshow(patches[display_index,5,:,:,2], cmap='gray')
        plt.title('T2')
        plt.subplot(234)
        plt.imshow(labels_img[display_index,0], cmap='gray')
        plt.title('Target label')
        plt.subplot(235)
        #plt.imshow(y_coords[display_index])
        #plt.title('Y-Coords')        
        plt.savefig('/home/andy/projects/mskProj/DeepPriors_package/One_Patch_example_MSKCC_16-328_1_09687_20100716-l_{}.png'.format(display_index))
         
    print('{} : Checking for NaN in image patches..'.format(procnum))    
    if np.any(np.isnan(patches)):
        print('{} : nan found in the input data batch for training..'.format(procnum))
        print(patches[np.isnan(patches)].shape)
        patches[np.isnan(patches)] = 0.0
    assert not np.any(np.isnan(patches)), 'STILL NANs!'          
    if np.any(~ np.isfinite(patches)):
        patches[~ np.isfinite(patches)] = 0.0
    assert np.all(np.isfinite(patches)), 'STILL Non-Finite Values!'
    #print('{} : Number of class 0 samples in whole batch: {}'.format(procnum, np.sum(labels[:,0])))
    #print('{} : Number of class 1 samples in whole batch: {}'.format(procnum, np.sum(labels[:,1])))
    #### SHUFFLE ####
    print('{} : Shuffling data..'.format(procnum))
    shuffleOrder = np.arange(patches.shape[0])
    np.random.shuffle(shuffleOrder)
    patches = patches[shuffleOrder]
    labels = labels[shuffleOrder]  
    if len(all_coordinates) > 0:
        all_coordinates = all_coordinates[shuffleOrder]
    if len(TPM_patches) > 0:
        TPM_patches = TPM_patches[shuffleOrder]        
    #--------------------------------------------------------------------------
    # Preprocess patches according to model needs... (BreastMask vs Segmenter...)
    # Resize Context Patch
    if using_unet:
        context = []
    else:
        context = np.array(patches[:,:,:,:,0],'float')
        context = resize(image=context, order=1, 
                             output_shape=(context.shape[0],context.shape[1],context.shape[2]/3,context.shape[3]/3), 
                             anti_aliasing=True, preserve_range=True )    
    # Crop Detail Patch
    if model_crop > 2:
        patches = np.array(patches[:,:,model_crop/2:-model_crop/2,model_crop/2:-model_crop/2,:])
    print('{} : ------- Finished data sampling -------'.format(procnum))    
    return_dict[procnum] = context, patches, labels, all_coordinates, TPM_patches


  
def sampleTestData(TPM_channel, testChannels, testLabels, subjectIndex, output_classes, output_dpatch, shape, use_coordinates):
   
    xend = output_dpatch[0] * int(round(float(shape[0])/output_dpatch[0] + 0.5)) 
    if shape[1] == output_dpatch[1]:
        yend = output_dpatch[1]
    else:
        yend = output_dpatch[1] * int(round(float(shape[1])/output_dpatch[1] + 0.5)) 
    if shape[2] == output_dpatch[2]:
        zend = output_dpatch[2]
    else:           
        zend = output_dpatch[2] * int(round(float(shape[2])/output_dpatch[2] + 0.5))
    voxelCoordinates = []
    # Remember in python the end is not included! Last voxel will be the prior-to-last in list.
    # It is ok if the center voxel is outside the image, PROPER padding will take care of that (if outside the image size, then we need larger than dpatch/2 padding)
    
    for x in range(output_dpatch[0]/2,xend,output_dpatch[0]): 
        for y in range(output_dpatch[1]/2,yend,output_dpatch[1]):
            for z in range(output_dpatch[2]/2,zend,output_dpatch[2]):
                voxelCoordinates.append([x,y,z])
    
    if len(TPM_channel) > 0:
      TPM_patches = extract_TPM_patches(TPM_channel, subjectIndex, [voxelCoordinates], output_dpatch, [shape])
    else:
      TPM_patches = []
    if len(testLabels) > 0:
      labels = np.array(extractLabels(testLabels, subjectIndex, [voxelCoordinates], output_dpatch,[shape]))
      labels = to_categorical(labels.astype(int),output_classes)
    else:
      labels = []
    if use_coordinates:
      spatial_coordinates = extractCoordinates([shape], [voxelCoordinates], output_dpatch) 
    else:
      spatial_coordinates = []
    #print("Finished extracting " + str(n_patches) + " patches, from "  + str(n_subjects) + " subjects and " + str(num_channels) + " channels. Timing: " + str(round(end-start,2)) + "s")
    return TPM_patches, labels, voxelCoordinates, spatial_coordinates, shape        



def getForegroundBackgroundVoxels(nifti_label, data_label, target_shape):
    "NOTE: img in MRICRON starts at (1,1,1) and this function starts at (0,0,0), so points do not match when comparing in MRICRON. Add 1 to all dimensions to match in mricron. Function works properly though"
    shape = nifti_label.shape
    if np.array(shape != target_shape).any():
      data = resize(data_label, order=0, output_shape=target_shape, preserve_range=True, anti_aliasing=True, mode='reflect') 
    else:
      data = data_label
    if np.sum(data) == 0:
      data = resize(data_label,order=1, output_shape=target_shape, preserve_range=True, anti_aliasing=True, mode='reflect') 
      data[data > 0] = 1
    nifti_label.uncache()    
    foregroundVoxels = np.argwhere(data>0)
    return foregroundVoxels
      
def getBodyVoxels(channel, percentile_voxel_intensity_sample_benigns):
    '''Get vector of voxel coordinates for all voxel values > 0'''
    "e.g. groundTruthChannel = '/home/hirsch/Documents/projects/ATLASdataset/native_part2/c0011/c0011s0006t01/c0011s0006t01_LesionSmooth_Binary.nii.gz'"
    "NOTE: img in MRICRON starts at (1,1,1) and this function starts at (0,0,0), so points do not match when comparing in MRICRON. Add 1 to all dimensions to match in mricron. Function works properly though"
    img = nib.load(channel)
    data = img.get_data()  
    res = img.header['pixdim'][1:4]
    shape = img.shape 
    if res[1] > 0.6:    
      target_res = [res[0],res[1]/2.,res[2]/2.]
      out_shape = np.floor([float(s)*r1/r2 for s,r1,r2 in zip(shape, res, target_res)])
      data = resize(data,order=1, output_shape=out_shape, preserve_range=True, anti_aliasing=True, mode='reflect')
    img.uncache()    
    bodyVoxels = np.argwhere(data > np.percentile(data, percentile_voxel_intensity_sample_benigns))
    return bodyVoxels

####################################### METRIC FUNCTIONS #################################################

def weighted_generalized_dice_completeImages(img1,img2,penalty_MATRIX):
    classes = np.array(range(0,len(penalty_MATRIX)), dtype='int8')   
    dice = []
    
    for i in classes:
        dice_2 = []
        #DICE = 2*np.sum(np.multiply(img1==i,img2==i))/float(np.sum(img1==i)+np.sum(img2==i))
        for j in classes:
            wDice = 2*np.sum(np.multiply(img1==i,img2==j) * penalty_MATRIX[i,j] )/float(np.sum(img1==i)+np.sum(img2==j))
            dice_2.append(wDice)
        dice.append(np.sum(dice_2)) 
    return np.sum(dice)/len(classes), [round(x,2) for x in dice]

def dice_completeImages(img1,img2):
    return(2*np.sum(np.multiply(img1>0,img2>0))/float(np.sum(img1>0)+np.sum(img2>0)))
       
def generalized_dice_completeImages(img1,img2):
    assert img1.shape == img2.shape, 'Images of different size!'
    #assert (np.unique(img1) == np.unique(img2)).all(), 'Images have different classes!'
    classes = np.array(np.unique(img1), dtype='int8')   
    if len(classes) < len(np.array(np.unique(img2), dtype='int8')   ):
      classes = np.array(np.unique(img2), dtype='int8')   
    dice = []
    for i in classes:
        dice.append(2*np.sum(np.multiply(img1==i,img2==i))/float(np.sum(img1==i)+np.sum(img2==i)))   
    return np.sum(dice)/len(classes), [round(x,2) for x in dice]

def dice_coef_multilabel_bin0(y_true, y_pred):
    index = 0
    dice = dice_coef(y_true[:,:,:,:,index], K.round(y_pred[:,:,:,:,index]))
    return dice

def dice_coef_multilabel_bin1(y_true, y_pred):
    index = 1
    dice = dice_coef(y_true[:,:,:,:,index], K.round(y_pred[:,:,:,:,index]))
    return dice

################################## DOCUMENTATION FUNCTIONS ################################################


def plot_training(session,losses, metrics,val_performance,full_segm_DICE, smooth=50, loss_name = ['Multiclass Dice'], class_names = ['Air','GM','WM','CSF','Bone','Skin']):

    losses_df = pd.DataFrame(losses)
    losses_df.columns=loss_name
    
    losses_mv_avg = losses_df.rolling(smooth,center=False).mean()
    metrics_df = pd.DataFrame(metrics)
    metrics_df.columns = class_names
    color_dict = {'Air':'black','GM':'blue','WM':'green','CSF':'yellow','Bone':'orange','Skin':'red'}
    metrics_mv_avg = metrics_df.rolling(smooth,center=False).mean()
    
    n_plots = 2 + np.sum([int(x) for x in [2*(len(val_performance) > 0), len(full_segm_DICE) > 0]])
            
    f, axarr = plt.subplots(n_plots, sharex=False, figsize=(8,10))
    losses_mv_avg.plot(ax=axarr[0])
    axarr[0].set_title(session)
    metrics_mv_avg.plot(ax=axarr[1], color=[color_dict.get(x, '#333333') for x in metrics_mv_avg.columns])
    #axarr[1].plot(metrics_mv_avg)
    #axarr[1].set_title('Single Class Dice Loss')
    axarr[1].set_xlabel('Training Iterations')
    axarr[1].legend(loc='upper left')
       
    if len(val_performance) > 0  :
    
        loss_val = [x[0] for x in val_performance]
        metrics_val = [x[1:len(x)] for x in val_performance]
        
        loss_val_df = pd.DataFrame(loss_val)
        loss_val_df.columns=loss_name
        #loss_val_df = loss_val_df.rolling(smooth,center=False).mean()
        metrics_val_df = pd.DataFrame(metrics_val)
        metrics_val_df.columns = class_names
        #metrics_val_df = metrics_val_df.rolling(smooth,center=False).mean()
        loss_val_df.plot(ax=axarr[2])
        #axarr[2].set_title(loss_name[0])
        metrics_val_df.plot(ax=axarr[3], color=[color_dict.get(x, '#333333') for x in metrics_mv_avg.columns])
        #axarr[1].plot(metrics_mv_avg)
        #axarr[3].set_title('Single Class Dice Loss')
        #axarr[3].set_xlabel('Training Iterations')
        
        axarr[3].legend(loc='upper left')
    
    if len(full_segm_DICE) > 0:
        
        full_segm_DICE = pd.DataFrame(full_segm_DICE)
        full_segm_DICE.columns=['Full Segmentation DICE']
        full_segm_DICE.plot(ax=axarr[n_plots-1],style='-o',color='green')
        axarr[n_plots-1].legend(loc='lower right')
        

def my_logger(string, logfile, print_out=True):
    f = open(logfile,'a')
    f.write('\n' + str(string))
    f.close()
    if print_out:
        print(string)
    

def start_training_session_logger(logfile,threshold_EARLY_STOP, TPM_channel, load_model,saveSegmentation,path_to_model,model,dropout, trainChannels, trainLabels, validationChannels, validationLabels, testChannels, testLabels, num_iter, epochs, n_patches, n_patches_val, n_subjects, samplingMethod_train, size_minibatches, n_full_segmentations, epochs_for_fullSegmentation, size_test_minibatches):
    my_logger('#######################################  NEW TRAINING SESSION  #######################################', logfile)    
    my_logger(trainChannels, logfile)
    my_logger(trainLabels, logfile)
    my_logger(validationChannels, logfile)        
    my_logger(validationLabels, logfile)  
    my_logger(testChannels, logfile) 
    my_logger(testLabels, logfile)
    my_logger('TPM channel (if given):', logfile)
    my_logger(TPM_channel, logfile)
    my_logger('Session parameters: ', logfile)
    my_logger('[num_iter, epochs, n_patches, n_patches_val, n_subjects, samplingMethod_train, size_minibatches, n_full_segmentations, epochs_for_fullSegmentation, size_test_minibatches]', logfile)
    my_logger([num_iter, epochs, n_patches, n_patches_val, n_subjects, samplingMethod_train, size_minibatches, n_full_segmentations, epochs_for_fullSegmentation, size_test_minibatches], logfile)
    my_logger('Dropout for last two fully connected layers: ' + str(dropout), logfile)
    my_logger('Model loss function: ' + str(model.loss), logfile)
    my_logger('Model number of parameters: ' + str(model.count_params()), logfile)
    my_logger('Optimizer used: ' +  str(model.optimizer.from_config), logfile)
    my_logger('Optimizer parameters: ' + str(model.optimizer.get_config()), logfile)
    my_logger('Save full head segmentation of subjects: ' + str(saveSegmentation), logfile)
    my_logger('EARLY STOP Threshold last 3 epochs: ' + str(threshold_EARLY_STOP), logfile)
    if load_model:
        my_logger("USING PREVIOUSLY SAVED MODEL -  Model retrieved from: " + path_to_model, logfile)


class LossHistory_multiDice6(keras.callbacks.Callback):
    def on_train_begin(self, logs={}):
        self.losses = []
        self.dice = []
        self.metrics = []

    def on_batch_end(self, batch, logs={}):
        self.dice = []
        self.losses.append(logs.get('loss'))
        self.dice.append(logs.get('dice_coef_multilabel0'))
        self.dice.append(logs.get('dice_coef_multilabel1'))
        self.dice.append(logs.get('dice_coef_multilabel2'))
        self.dice.append(logs.get('dice_coef_multilabel3'))
        self.dice.append(logs.get('dice_coef_multilabel4'))
        self.dice.append(logs.get('dice_coef_multilabel5'))
        self.metrics.append(self.dice)



class LossHistory_multiDice2(keras.callbacks.Callback):
    def on_train_begin(self, logs={}):
        self.losses = []
        self.dice = []
        self.metrics = []

    def on_batch_end(self, batch, logs={}):
        self.dice = []
        self.losses.append(logs.get('loss'))
        self.dice.append(logs.get('dice_coef_multilabel_bin0'))
        self.dice.append(logs.get('dice_coef_multilabel_bin1'))
        self.metrics.append(self.dice)

        
        
################################### SEGMENTATION FUNCTIONS ##################################################


def fullSegmentation_Flexible(wd, penalty_MATRIX, dice_compare, dsc, model, testChannels, testLabels,TPM_channel, subjectIndex, output_classes,
                              dpatch, size_minibatches,logfile, epoch, use_coordinates, saveSegmentation = False, full_evaluation = False,debug=False):    

    subjectIndex = [subjectIndex]
    subject_channels = []
    for modality in testChannels:
      channelsFile = open(modality,"r")   
      ch = channelsFile.readlines()
      subject_channels.append(ch[subjectIndex[0]][:-1])#
      channelsFile.close()
    my_logger('Segmenting subject with channels: ' + str(subject_channels), logfile)  
    
    subID = '_'.join(subject_channels[0].split('/')[-2:])
    
    images = []
    for channel in subject_channels: 
      proxy_img = nib.load(channel)
      X = proxy_img.get_data()
      images.append(X)
    X = np.stack(images, axis=3)

    shape = X.shape
    original_shape = shape

    tpm_nii = nib.load(TPM_channel)      
    TPM_data = tpm_nii.get_data()  
    # Because the TPM has a different size, have to resize it to have the same size as the input MRI.
    # Actually I need to do this too for extracting patches. Else the center is right, but the dimensions are wrong!
    # Need to resize the TPM for each MRI input! 
    TPM_data.shape
    TPM_data = resize(TPM_data, original_shape[:-1], order=1, preserve_range=True, anti_aliasing=True)

    Res = X[:,:,:,0] + TPM_data*10
    img = nib.Nifti1Image(Res, proxy_img.affine)
    nib.save(img, '/home/andy/projects/mskProj/DeepPriors_package/test_t1-TPM.nii')
    
    affine = proxy_img.affine
    if shape[0]*shape[1]*shape[2] > 59*211*211:  # if shape exceeeds 55*261*261
      # Set boundaries for maximum allowed shape
      a = np.max([0,(shape[0] - 59)])/2   
      b = np.max([0,(shape[1] - 211)])/2
      c = np.max([0,(shape[2] - 211)])/2    
      X = X[a:shape[0]-a,:,:,:]
      X = X[:,b:shape[1]-b,:,:]
      X = X[:,:,c:shape[1]-c,:]
  
      TPM_data = TPM_data[a:shape[0]-a,:,:]
      TPM_data = TPM_data[:,b:shape[1]-b,:]
      TPM_data = TPM_data[:,:,c:shape[1]-c]
           
    shape = X.shape
    X = X.reshape( (1,) + shape)
    TPM_data = TPM_data[6:TPM_data.shape[0]-6, 33:TPM_data.shape[1]-33, 33:TPM_data.shape[2]-33]
      
    if use_coordinates:
      coords_shape = [X.shape[1] - 12 ,X.shape[2] - 66, X.shape[3] - 66]
      y_coords = np.tile(np.array([range(6,coords_shape[2]+6)]).transpose(), (1,coords_shape[1]))
      y_coords = y_coords/float(X.shape[2])
      #z_coords = np.tile(range(33,coords_shape[1]+33), (coords_shape[2],1))
      #z_coords = z_coords/float(X.shape[3])
      y_coords = np.repeat(y_coords[np.newaxis, :,: ], coords_shape[0], axis=0)    
      #z_coords = np.repeat(z_coords[np.newaxis, :,: ], coords_shape[0], axis=0)   
      y_coords = y_coords.reshape((1,) + y_coords.shape + (1,))
      #z_coords = z_coords.reshape((1,)+ z_coords.shape + (1,))  
      #print('X shape: {} , y_coords shape {} , z_coords shape {}'.format(X_padded.shape, y_coords.shape, z_coords.shape))
      
    T1post = X[:,:,:,:,0].reshape(X[:,:,:,:,0].shape + (1,))
    T1pre = X[:,:,:,:,1].reshape(X[:,:,:,:,1].shape + (1,))
    T2 = X[:,:,:,:,2].reshape(X[:,:,:,:,2].shape + (1,))
    TPM_data = TPM_data.reshape((1,) + TPM_data.shape + (1,))
    
    if debug:
      img = nib.Nifti1Image(y_coords[0,:,:,:,0], np.diag((1,1,1,0)))
      nib.save(img,'/home/andy/projects/mskProj/DeepPriors_package/y_coords_fullSegm.nii' )
      img = nib.Nifti1Image(T1post[0,:,:,:,0], np.diag((1,1,1,0)))
      nib.save(img,'/home/andy/projects/mskProj/DeepPriors_package/T1post_fullSegm.nii' )
      img = nib.Nifti1Image(T1pre[0,:,:,:,0], np.diag((1,1,1,0)))
      nib.save(img,'/home/andy/projects/mskProj/DeepPriors_package/T1pre_fullSegm.nii' )
      img = nib.Nifti1Image(T2[0,:,:,:,0], np.diag((1,1,1,0)))
      nib.save(img,'/home/andy/projects/mskProj/DeepPriors_package/T2_fullSegm.nii' )
          
      
    yhat = model.predict([T1post, T1pre, T2, TPM_data])

    #y = np.argmax(yhat, axis=4)   # For classification output
    y = yhat[:,:,:,:,1]            # For logits for class 2
    print('y shape: {}'.format(y.shape))
    y = y.reshape(y.shape[1],y.shape[2],y.shape[3])
    #y = y.reshape(shape[0]-24,shape[1]-78,shape[2]-78)
    

    y_out = np.zeros((original_shape[0],original_shape[1],original_shape[2]))

    try:
        y_out[  :, abs(original_shape[1] -y.shape[1])/2:original_shape[1] - abs(original_shape[1] -y.shape[1])/2,
                   abs(original_shape[2] -y.shape[2])/2:original_shape[2] - abs(original_shape[2] -y.shape[2])/2] = y[abs(original_shape[0] -y.shape[0])/2:y.shape[0] - abs(original_shape[0] -y.shape[0])/2,:,:]
    
    except:
        y_out[abs(original_shape[0] -y.shape[0])/2:original_shape[0] - abs(original_shape[0] -y.shape[0])/2, 
              abs(original_shape[1] -y.shape[1])/2:original_shape[1] - abs(original_shape[1] -y.shape[1])/2,
              abs(original_shape[2] -y.shape[2])/2:original_shape[2] - abs(original_shape[2] -y.shape[2])/2] = y
    
    img = nib.Nifti1Image(y_out, affine)
    segmentationName = '/predictions/' + subID + str(epoch)
    output = wd +'/' + segmentationName + '.nii'
    nib.save(img, output)
    #my_logger('Saved segmentation of subject at: ' + output, logfile)
      
   
#penalty_MATRIX = cfg.penalty_MATRIX
#TPM_channel = cfg.TPM_channel
#output_classes = cfg.output_classes
#segmentation_dpatch = cfg.segmentation_dpatch
#size_test_minibatches =  cfg.size_test_minibatches
#output_probability =  cfg.output_probability
#saveSegmentation =  cfg.saveSegmentation    
#size_minibatches = cfg.size_test_minibatches
#subjectIndex = 0
#use_coordinates= cfg.use_coordinates
#epoch = 0
#intensity_normalization_method = cfg.intensity_normalization_method    
#model_patch_reduction = cfg.model_patch_reduction
#model_crop = cfg.model_crop
#resolution = cfg.resolution
#testChannels = cfg.testChannels
#testLabels = cfg.testLabels
#
#testChannels = cfg.segmentChannels
#testLabels = cfg.segmentLabels
#intensity_normalization_method = cfg.intensity_normalization_method
#MRI_PATH = cfg.MRI_PATH
#OUTPUT_PATH = cfg.OUTPUT_PATH

def fullSegmentation(wd, penalty_MATRIX, resolution, OUTPUT_PATH, TPM_channel, dice_compare, dsc, smooth_dice_scores,foreground_percent_list, model, 
                     testChannels, testLabels, subjectIndex, output_classes, segmentation_dpatch, size_minibatches,output_probability, use_coordinates, 
                     intensity_normalization_method, MRI_PATH, model_patch_reduction, model_crop, epoch, using_breastMaskModel=False, MASK_BREAST=False,
                     using_Unet=False, using_unet_breastMask=False):    
    
    output_dpatch = segmentation_dpatch[0] - model_patch_reduction[0], segmentation_dpatch[1] - model_patch_reduction[1], segmentation_dpatch[2] - model_patch_reduction[2]
    if len(testLabels) == 0:
        dice_compare = False         

    subjectIndex = [subjectIndex]
    num_channels = len(testChannels)
    firstChannelFile = open(testChannels[0],"r")   
    ch = firstChannelFile.readlines()
    subjectGTchannel = ch[subjectIndex[0]][:-1]
    subID = subjectGTchannel.split('/')[-2] + '_' + subjectGTchannel.split('/')[-1].split('.nii')[0]
    print('SEGMENTATION : Segmenting subject: ' + str(subID))  
    segmentationName =  subID + '_epoch' + str(epoch)
    if len(OUTPUT_PATH) > 0:
        output = OUTPUT_PATH + '/' + segmentationName + '.nii.gz' 
    else:
        output = wd + '/predictions/' + segmentationName + '.nii.gz'
    if os.path.exists(output):
      print('SEGMENTATION : Segmentation already done. Skip')
      return [None]*5
    
    firstChannelFile.close()      
    proxy_img = nib.load(subjectGTchannel)
    shape = proxy_img.shape
    affine = proxy_img.affine      
    res = proxy_img.header['pixdim'][1:4]

    if resolution == 'high':
        if res[1] > 0.6:    
          target_res = [res[0],res[1]/2.,res[2]/2.]
          shape = [int(x) for x in np.floor([float(s)*r1/r2 for s,r1,r2 in zip(shape, res, target_res)])]
        else:
          target_res = res
    elif resolution == 'low':
        if res[1] < 0.6:    
          target_res = [res[0],res[1]*2.,res[2]*2.]
          shape = [int(x) for x in np.floor([float(s)*r1/r2 for s,r1,r2 in zip(shape, res, target_res)])]
        else:
          target_res = res          
          
    print('SEGMENTATION : Sampling data..')  
    TPM_patches, labels, voxelCoordinates, spatial_coordinates, shape = sampleTestData(TPM_channel, testChannels, testLabels, subjectIndex, output_classes, 
                                                                                       output_dpatch, shape, use_coordinates)    
    affine = np.diag(list(target_res) + [0])        
    n_minibatches = 0 # min(0,len(voxelCoordinates)/size_minibatches) 
    total_number_of_patches = (len(voxelCoordinates)-n_minibatches*size_minibatches)  
    
    #########################################################################
    print('SEGMENTATION : Extracting {} image patches..'.format(total_number_of_patches))

    patches = extractImagePatch_parallelization(MRI_PATH, testChannels[0], subjectIndex[0], voxelCoordinates, shape, segmentation_dpatch, 
                                                intensity_normalization_method, fullSegmentationPhase=True)    


    print('SEGMENTATION : Finished sampling data.')
#    if debug_coords:
#        patches = np.ones(patches.shape)
#        spatial_coordinates = np.zeros(spatial_coordinates.shape)
    
    INPUT_DATA = []  
    
    # NEED TO ADAPT FOR BREAST MASK MODEL: Inputs: Context (no resizing, 13,75,75) and spatial coordinates.
    if using_breastMaskModel:
        INPUT_DATA.append(patches[:,:,:,:,0].reshape(patches[:,:,:,:,0].shape + (1,)))  
        INPUT_DATA.append(spatial_coordinates)    
  
    elif using_Unet:  # This means the model is the U-Net
        if using_unet_breastMask:
            patches = patches[:,:,:,:,0]
            patches = patches.reshape(patches.shape + (1,))
        INPUT_DATA.append(patches)
#        if len(TPM_patches) > 0:
#            INPUT_DATA.append(TPM_patches[:,:,:,:].reshape(TPM_patches[:,:,:,:].shape + (1,)))   
        if len(spatial_coordinates) > 0:
            INPUT_DATA.append(spatial_coordinates)            
            
        
    else:
        # Context
        context = np.array(patches[:,:,:,:,0],'float')
        context = resize(image=context, order=1, 
                             output_shape=(context.shape[0],context.shape[1],context.shape[2]/3,context.shape[3]/3), 
                             anti_aliasing=True, preserve_range=True )
        INPUT_DATA.append(context.reshape(context.shape + (1,)))        
        
        for jj in range(patches.shape[-1]):
            INPUT_DATA.append(patches[:,:,model_crop/2:-model_crop/2,model_crop/2:-model_crop/2,jj].reshape(patches[:,:,model_crop/2:-model_crop/2,model_crop/2:-model_crop/2,jj].shape + (1,)))  
        if len(TPM_patches) > 0:
            INPUT_DATA.append(TPM_patches[:,:,:,:].reshape(TPM_patches[:,:,:,:].shape + (1,)))   
        if len(spatial_coordinates) > 0:
            INPUT_DATA.append(spatial_coordinates)    
    
    print("SEGMENTATION : Finished preprocessing data for segmentation.")
    #########################################################################
      
    prediction = model.predict(INPUT_DATA, verbose=1, batch_size=size_minibatches)
     
    ##########  Output binary ############          
    indexes = []
    class_pred = np.argmax(prediction, axis=4)
    indexes.extend(class_pred)     
           
    head = np.ones(shape, dtype=np.float32)  # same size as input head, start index for segmentation start at 26,26,26, rest filled with zeros....
    i = 0
    for x,y,z in voxelCoordinates:
        patch_shape = head[x-output_dpatch[0]/2:min(x+(output_dpatch[0]/2+output_dpatch[0]%2), shape[0]),
                           y-output_dpatch[1]/2:min(y+(output_dpatch[1]/2+output_dpatch[1]%2), shape[1]),
                           z-output_dpatch[2]/2:min(z+(output_dpatch[2]/2+output_dpatch[2]%2), shape[2])].shape
        #print(np.array(indexes[i])[0:patch_shape[0], 0:patch_shape[1],0:patch_shape[2]])
        head[x-output_dpatch[0]/2:min(x+(output_dpatch[0]/2+output_dpatch[0]%2), shape[0]),
             y-output_dpatch[1]/2:min(y+(output_dpatch[1]/2+output_dpatch[1]%2), shape[1]),
             z-output_dpatch[2]/2:min(z+(output_dpatch[2]/2+output_dpatch[2]%2), shape[2])] = np.array(indexes[i])[0:patch_shape[0], 
                                                                                                                   0:patch_shape[1],
                                                                                                                   0:patch_shape[2]]
        i = i+1
    #img_binary = nib.Nifti1Image(head, affine)
    img_binary = head

#    if(saveSegmentation):
#        nib.save(img_binary, output)
#        my_logger('Saved segmentation of subject at: ' + output, logfile)

    ##########  Output probabilities ############
    if output_probability:
        indexes = []        
        class_pred = prediction[:,:,:,:,1]
        indexes.extend(class_pred)     
               
        head = np.ones(shape, dtype=np.float32)  # same size as input head, start index for segmentation start at 26,26,26, rest filled with zeros....
        i = 0
        for x,y,z in voxelCoordinates:
            patch_shape = head[x-output_dpatch[0]/2:min(x+(output_dpatch[0]/2+output_dpatch[0]%2), shape[0]),
                               y-output_dpatch[1]/2:min(y+(output_dpatch[1]/2+output_dpatch[1]%2), shape[1]),
                               z-output_dpatch[2]/2:min(z+(output_dpatch[2]/2+output_dpatch[2]%2), shape[2])].shape
            #print(np.array(indexes[i])[0:patch_shape[0], 0:patch_shape[1],0:patch_shape[2]])
            head[x-output_dpatch[0]/2:min(x+(output_dpatch[0]/2+output_dpatch[0]%2), shape[0]),
                 y-output_dpatch[1]/2:min(y+(output_dpatch[1]/2+output_dpatch[1]%2), shape[1]),
                 z-output_dpatch[2]/2:min(z+(output_dpatch[2]/2+output_dpatch[2]%2), shape[2])] = np.array(indexes[i])[0:patch_shape[0], 
                                                                                                                       0:patch_shape[1],
                                                                                                                       0:patch_shape[2]]
            i = i+1
        #img_probs = nib.Nifti1Image(head, affine)
        img_probs = head

    foreground_percent = 0 
    score_smooth = 0    
    if dice_compare:
      LABEL_CHANNEL = open(testLabels).readlines()[subjectIndex[0]][:-1]
      print('SEGMENTATION : Comparing with Ground-Truth label: {}'.format(LABEL_CHANNEL))
      if 'BENIGN' in LABEL_CHANNEL:
        dice_compare = False
        label_data = np.zeros((img_binary.shape))
        foreground_percent = np.sum(img_binary)/float(head.size)    
        foreground_percent_list.append(foreground_percent)    
      elif np.sum(nib.load(LABEL_CHANNEL).get_data()) == 0:   
        print('SEGMENTATION : Empty label!')
        label_data = np.zeros((img_binary.shape))
      else:
        label_data = nib.load(LABEL_CHANNEL).get_data()    
        print('SEGMENTATION : Label shape = {}'.format(label_data.shape))
        print('SEGMENTATION : Segmentation shape = {}'.format(shape))
        
        print('label data sum = {}'.format(np.sum(label_data)))
        
        if label_data.shape != shape:  
          print('SEGMENTATION : Resizing.')
          label_data = resize(label_data, output_shape=shape, preserve_range=True, anti_aliasing=True, order=0)
          print('SEGMENTATION: np.sum after resizing: {}'.format(np.sum(label_data)))
          
        if np.any(np.isnan(label_data)):
          label_data[np.isnan(label_data)] = 0
        # Get only segmented slice
        try:  
#          if MASK_BREAST:  
#              breast_masks_path = '/home/deeperthought/kirby_MSK/BreastMasks/alignedNii-Aug2019/'
#              print('SEGMENTATION : Post-processing segmentation with provided breast mask from ' + breast_masks_path)
#              exam = output.split('/')[-1].split('_T1')[0]
#              side = output.split('T1_')[-1].split('_post')[0]
#              bm = [x for x in os.listdir(breast_masks_path) if exam in x and side[0] in x]
#              bm = nib.load(breast_masks_path + bm[0]).get_data()
#              img_probs *= bm
#              img_binary *= bm  
              
          slice_of_interest = np.argwhere(label_data>0)[0][0]
          print('SLICE : {}'.format(slice_of_interest))
          label_data = label_data[slice_of_interest, :, :]
          img_binary_2d = img_binary[slice_of_interest,:,:]
          img_probs_2d = img_probs[slice_of_interest,:,:]

        except:
          print('SEGMENTATION : Target label removed on resizing, subject: {}'.format(LABEL_CHANNEL))
          sys.exit(0)
          dice_compare = False
                
        score = generalized_dice_completeImages(img_binary_2d, label_data)
        score_smooth = Generalised_dice_coef_multilabel2(label_data, img_probs_2d)/2.
        dsc.append(score[0])
        smooth_dice_scores.append(score_smooth)
        print(dsc[-1])
        print('per class dice score: {}'.format(score[1]))
        print('mean DCS so far:' + str(np.mean(dsc)))
        print('smooth_dice score: {}'.format(score_smooth))
        print('mean SMOOTH_DCS so far:' + str(np.mean(smooth_dice_scores)))
    
    if output_probability:
        img_probs = nib.Nifti1Image(img_probs, affine)
        return img_probs, output, dsc, score_smooth,foreground_percent
    else:
        img_binary = nib.Nifti1Image(img_binary, affine)
        return img_binary, output, dsc, score_smooth,foreground_percent
    
def dice_coef(y_true_f, y_pred_f):
    smooth = 1e-6
    #y_true_f = y_true.reshape(np.prod(y_true.shape))
    #y_pred_f = y_pred.reshape(np.prod(y_pred.shape))
    intersection = np.sum(y_true_f * y_pred_f)
    return (2. * intersection + smooth) / (np.sum(y_true_f**2) + np.sum(y_pred_f**2) + smooth)

def Generalised_dice_coef_multilabel2(label_data, img_probs, numLabels=2):
    y_true = np.zeros((np.prod(label_data.shape), 2))
    y_true[:,1] = label_data.reshape(np.prod(label_data.shape))
    y_true[:,0] = 1 - y_true[:,1]
    y_pred = np.zeros((np.prod(img_probs.shape), 2))
    y_pred[:,1] = img_probs.reshape(np.prod(img_probs.shape))
    y_pred[:,0] = 1 - y_pred[:,1]
    
    dice=0
    for index in range(numLabels):
        dice -= dice_coef(y_true[:,index], y_pred[:,index])
    return dice*-1     

#######################################################################################################################################
#configFile = '/home/deeperthought/Projects/MultiPriors_MSKCC/configFiles/UNet_v0_TumorSegmenter.py'
#workingDir = '/home/deeperthought/Projects/MultiPriors_MSKCC/'

def segment(configFile,workingDir):

    scripts_path = configFile.split('configFiles')[0] + 'scripts'
    sys.path.append(scripts_path)
    
    path = '/'.join(configFile.split('/')[:-1])
    configFileName = configFile.split('/')[-1][:-3]   
    sys.path.append(path)
    cfg = __import__(configFileName)
           
    epoch = 0 #int(cfg.path_to_model.split('.')[-2][cfg.path_to_model.split('.')[-2].find('epoch') + 5 : ]) + 1
        
    os.chdir(workingDir + '/training_sessions/')
    session = cfg.session
    wd = workingDir + '/training_sessions/' +session
    print('\n CURRENTLY IN SESSION {} \n'.format(session))
    if not os.path.exists(wd):    
        os.mkdir(session)
        os.mkdir(session + '/models')
        os.mkdir(session + '/predictions')
    os.chdir(wd)
    
    logfile = 'segmentations.log'
    dice_compare = cfg.dice_compare
    #cfg.TPM_channel = workingDir + cfg.TPM_channel
    cfg.segmentChannels = [workingDir + x for x in cfg.segmentChannels]
    if len(cfg.segmentLabels) > 0:
      cfg.segmentLabels = workingDir + cfg.segmentLabels 
      dice_compare = True
    if len(cfg.TPM_channel) != 0:
      cfg.TPM_channel = workingDir + cfg.TPM_channel
      dice_compare = False
    if cfg.output_classes == 6:
    	try:
    	    from MultiPriors_MSKCC_MultiScale import dice_coef_multilabel6, dice_coef_multilabel0,dice_coef_multilabel1,dice_coef_multilabel2,dice_coef_multilabel3,dice_coef_multilabel4,dice_coef_multilabel5
    	    my_custom_objects = {'dice_coef_multilabel6':dice_coef_multilabel6,
    				     'dice_coef_multilabel0':dice_coef_multilabel0,
    				     'dice_coef_multilabel1':dice_coef_multilabel1,
    				     'dice_coef_multilabel2':dice_coef_multilabel2,
    				     'dice_coef_multilabel3':dice_coef_multilabel3,
    				     'dice_coef_multilabel4':dice_coef_multilabel4,
    				     'dice_coef_multilabel5':dice_coef_multilabel5}
    		#custom_metrics =[dice_coef_multilabel6,dice_coef_multilabel0,dice_coef_multilabel1,dice_coef_multilabel2,dice_coef_multilabel3,dice_coef_multilabel4,dice_coef_multilabel5]
    		#my_custom_objects = dict(zip(np.sort(my_custom_objects.keys()), custom_metrics))
    	    model = load_model(cfg.path_to_model, custom_objects = my_custom_objects )
    	except:
    	    from MultiPriors_MSKCC_MultiScale import w_dice_coef_multilabel6, dice_coef_multilabel0,dice_coef_multilabel1,dice_coef_multilabel2,dice_coef_multilabel3,dice_coef_multilabel4,dice_coef_multilabel5
    	    my_custom_objects = {'w_dice_coef_multilabel6':w_dice_coef_multilabel6,
    					     'dice_coef_multilabel0':dice_coef_multilabel0,
    					     'dice_coef_multilabel1':dice_coef_multilabel1,
    					     'dice_coef_multilabel2':dice_coef_multilabel2,
    					     'dice_coef_multilabel3':dice_coef_multilabel3,
    					     'dice_coef_multilabel4':dice_coef_multilabel4,
    					     'dice_coef_multilabel5':dice_coef_multilabel5}
    elif cfg.output_classes == 2:
        try:
            from MultiPriors_Models_Collection import Generalised_dice_coef_multilabel2, dice_coef_multilabel_bin0,dice_coef_multilabel_bin1
            my_custom_objects = {'Generalised_dice_coef_multilabel2':Generalised_dice_coef_multilabel2,
                                 'dice_coef_multilabel_bin0':dice_coef_multilabel_bin0,
                                 'dice_coef_multilabel_bin1':dice_coef_multilabel_bin1}
        except:
            from MultiPriors_MSKCC_MultiScale import w_dice_coef_multilabel2, dice_coef_multilabel0,dice_coef_multilabel1
            my_custom_objects = {'w_dice_coef_multilabel2':w_dice_coef_multilabel2,
				     'dice_coef_multilabel0':dice_coef_multilabel0,
				     'dice_coef_multilabel1':dice_coef_multilabel1}
        model = load_model(cfg.path_to_model, custom_objects = my_custom_objects )

    full_segm_DICE = 0
    full_segm_SMOOTH_DICE = 0
    np.set_printoptions(precision=3)

    print("------------------------------------------------------")
    print("                 WHOLE SCAN SEGMENTATION")
    print("------------------------------------------------------")
    dsc = []
    smooth_dice_scores = []
    foreground_percent_list = []
    epoch_foreground_percent = []
            
#    if 'BreastSegmentor' in cfg.model:
#        using_breastMaskModel = True
#        print'USING BREASTMASK MODEL'
#    else:
#        using_breastMaskModel = False
    with open(cfg.segmentChannels[0]) as vl:
        n_segmentSubjects = len(vl.readlines())
    if cfg.test_subjects > n_segmentSubjects:
        print("Given number of subjects for test set (" + str(cfg.test_subjects) +") is larger than the amount of \
        subjects in test set (" +str(n_segmentSubjects)+ ")")
        cfg.test_subjects = n_segmentSubjects
        print('Using {} number of test subjects'.format(n_segmentSubjects))
    if len(cfg.list_subjects_fullSegmentation) == 0:
      list_subjects_fullSegmentation = range(cfg.test_subjects)
    else:
      list_subjects_fullSegmentation = cfg.list_subjects_fullSegmentation
    for subjectIndex in list_subjects_fullSegmentation: 
        t_segment = time.time()
        if cfg.full_segmentation_patches:
                        
            segmentation_img, output, dsc, score_smooth,foreground_percent = fullSegmentation(wd, cfg.penalty_MATRIX, cfg.resolution, cfg.OUTPUT_PATH, cfg.TPM_channel, dice_compare, dsc, smooth_dice_scores, 
                                                                                              foreground_percent_list, model, cfg.segmentChannels, cfg.segmentLabels, subjectIndex,
                                                                                              cfg.output_classes, cfg.segmentation_dpatch, cfg.size_test_minibatches,cfg.output_probability, 
                                                                                              cfg.use_coordinates, cfg.intensity_normalization_method, cfg.MRI_PATH, cfg.model_patch_reduction, cfg.model_crop, 
                                                                                              epoch, using_Unet=True, using_unet_breastMask=cfg.using_unet_breastMask)
            
        else:
            fullSegmentation_Flexible(wd, cfg.penalty_MATRIX, dice_compare, dsc, model, cfg.segmentChannels, cfg.segmentLabels, cfg.TPM_channel, subjectIndex, \
            cfg.output_classes, cfg.dpatch, cfg.size_test_minibatches, logfile, epoch, cfg.use_coordinates,cfg.saveSegmentation)

 
        if dice_compare:
            my_logger('--------------- TEST EVALUATION ---------------', logfile)
            my_logger('          Full segmentation evaluation of subject' + str(subjectIndex), logfile)
            my_logger('foreground_percent {}'.format(foreground_percent), logfile)
            my_logger('SMOOTH_DCS ' + str(score_smooth),logfile)   
            my_logger('DCS ' + str(dsc[-1]),logfile)
        if cfg.saveSegmentation:
            if segmentation_img != None:
                if cfg.save_as_nifti:
                    if len(cfg.OUTPUT_PATH) > 0:
                        output = cfg.OUTPUT_PATH + '/' + output.split('/')[-1]
                    nib.save(segmentation_img, output)
                else:
                    np.savez_compressed(output.replace('.nii.gz',''), segmentation_img.get_data())
            
       #if (dice_compare & len(dsc)>0):
        #  my_logger('DCS ' + str(dsc[-1]),logfile)
        print('Segmentation of subject took {} s'.format(time.time()-t_segment))
    my_logger('         FULL SEGMENTATION SUMMARY STATISTICS ', logfile)
    
    full_segm_DICE.append(np.mean(dsc))   
    full_segm_SMOOTH_DICE.append(np.mean(smooth_dice_scores))   
    my_logger('Overall DCS:   ' + str(full_segm_DICE[-1]),logfile)
    my_logger('Overall SMOOTH_DCS:   ' + str(full_segm_SMOOTH_DICE[-1]),logfile)
    epoch_foreground_percent.append(np.mean(foreground_percent_list))            
    my_logger('Epoch_foreground_percent {}'.format(epoch_foreground_percent[-1]), logfile)

############################# MODEL TRAINING AND VALIDATION FUNCTIONS ############################################

#batch = return_dict['TRAINING'][0]
#labels = return_dict['TRAINING'][1]
#TPM_patches = return_dict['TRAINING'][3] 
#coords = return_dict['TRAINING'][2]
#size_minibatches = cfg.size_minibatches

def train_validate_model_on_batch(model_name, model,context,batch,labels,coords,TPM_patches,size_minibatches,history,losses,metrics,output_classes,logfile=0, TRAINING_FLAG=True, using_unet_breastMask=False, verbose=False):
    batch_performance = []   
    INPUT_DATA = []
    
    if 'UNet' in model_name:
        if using_unet_breastMask:
            batch = batch[:,:,:,:,0]
            batch = batch.reshape(batch.shape + (1,))
        INPUT_DATA.append(batch)
        if len(coords) > 0:
          #coords = coords.reshape( (coords.shape[0],) + (1,) + coords.shape[1:] )  
          INPUT_DATA.append(coords)
       
    else:    
        # Context
        INPUT_DATA.append(context.reshape(context.shape + (1,)))    
        
        for jj in range(batch.shape[-1]):
          INPUT_DATA.append(batch[:,:,:,:,jj].reshape(batch[:,:,:,:,jj].shape + (1,)))
          
        if len(TPM_patches) > 0:
          INPUT_DATA.append(TPM_patches[:,:,:,:].reshape(TPM_patches[:,:,:,:].shape + (1,)))   
    
        if len(coords) > 0:
          #coords = coords.reshape( (coords.shape[0],) + (1,) + coords.shape[1:] )  
          INPUT_DATA.append(coords)

    ######### TRAINING ###########
    if TRAINING_FLAG:
        print('Training..')
        model.fit(INPUT_DATA, labels, verbose = 1, callbacks = [history], batch_size = size_minibatches)

        if verbose:
          freq = classesInSample(labels, output_classes)
          print("Sampled following number of classes in training MINIBATCH: " + str(freq))
      
        if logfile != 0:
            output_results = zip(['Train cost and metrics     ']*len(history.losses), history.losses, history.metrics)
            for line in output_results:
                my_logger(' '.join(map(str, line)),logfile,print_out=False)
            
    ######### VALIDATION ###########
    else:
        print('Validation..')
        batch_performance.append(model.evaluate(INPUT_DATA, labels, verbose=1, batch_size = size_minibatches))
        
    del batch
    del labels
    if TRAINING_FLAG:
        return history.losses, history.metrics
    else:    
        val_performance = np.mean(batch_performance, 0)
        my_logger('Validation cost and accuracy ' + str(val_performance),logfile)                    
        return list(val_performance)   


################################ MAIN TRAINING FUNCTION ###########################################
#configFile = '/home/deeperthought/Projects/MultiPriors_MSKCC/configFiles/configFile_UNet_3D_v4_STROKE_DATA.py'
#workingDir = '/home/deeperthought/Projects/MultiPriors_MSKCC/'

def train_test_model(configFile, workingDir):
    print(configFile)
    path = '/'.join(configFile.split('/')[:-1])
    print(path)
    configFileName = configFile.split('/')[-1][:-3]   
    sys.path.append(path)
    #sys.path.append(path.replace('configFiles','scripts'))
    cfg = __import__(configFileName)
    if len(cfg.TPM_channel) != 0:
      cfg.TPM_channel = workingDir + cfg.TPM_channel
    cfg.trainChannels = [workingDir + x for x in cfg.trainChannels]
    cfg.trainLabels = workingDir +cfg.trainLabels 
    cfg.testChannels = [workingDir + x for x in cfg.testChannels]
    cfg.testLabels = workingDir + cfg.testLabels
    cfg.validationChannels = [workingDir + x for x in cfg.validationChannels]
    cfg.validationLabels = workingDir +cfg.validationLabels
        
    cfg.dataset  # Use this to change extraction code for brain data...
    
    if cfg.load_model == False:

        if cfg.model == 'MultiPriors_v0':
            from MultiPriors_Models_Collection import MultiPriors_v0
            mp = MultiPriors_v0(cfg.output_classes, cfg.num_channels, cfg.L2, cfg.dropout, cfg.learning_rate, cfg.optimizer_decay, cfg.loss_function)
            model = mp.createModel()            
            model.summary()                             

        elif cfg.model == 'MultiPriors_v1':
            from MultiPriors_Models_Collection import MultiPriors_v1
            mp = MultiPriors_v1(cfg.output_classes, cfg.num_channels, cfg.L2, cfg.dropout, cfg.learning_rate, cfg.optimizer_decay, cfg.loss_function)
            model = mp.createModel()            
            model.summary()           

        elif cfg.model == 'MultiPriors_v2':
            from MultiPriors_Models_Collection import MultiPriors_v2
            mp = MultiPriors_v2(cfg.output_classes, cfg.num_channels, cfg.L2, cfg.dropout, cfg.learning_rate, cfg.optimizer_decay, cfg.loss_function)
            model = mp.createModel()            
            model.summary()        

        elif cfg.model == 'MultiPriors_v2_Big':
            from MultiPriors_Models_Collection import MultiPriors_v2_Big
            mp = MultiPriors_v2_Big(cfg.output_classes, cfg.num_channels, cfg.L2, cfg.dropout, cfg.learning_rate, cfg.optimizer_decay, cfg.loss_function)
            model = mp.createModel()            
            model.summary()      

        elif cfg.model == 'MultiPriors_v2_Big_BreastMask':
            from MultiPriors_Models_Collection import MultiPriors_v2_Big_BreastMask
            mp = MultiPriors_v2_Big_BreastMask(cfg.output_classes, cfg.num_channels, cfg.L2, cfg.dropout, cfg.learning_rate, cfg.optimizer_decay, cfg.loss_function)
            model = mp.createModel()            
            model.summary()                  
            

        elif cfg.model == 'MultiPriors_TEST':
            from MultiPriors_Models_Collection import MultiPriors_v2_ContextOutput
            mp = MultiPriors_v2_ContextOutput(cfg.output_classes, cfg.num_channels, cfg.L2, cfg.dropout, cfg.learning_rate, cfg.optimizer_decay, cfg.loss_function)
            model = mp.createModel()            
            model.summary()     

        elif cfg.model == 'BreastSegmentor_v0':
            from BreastMask_Models_Collection import BreastSegmentor_v0
            mp = BreastSegmentor_v0(cfg.output_classes, cfg.num_channels, cfg.L2, cfg.dropout, cfg.learning_rate, cfg.optimizer_decay, cfg.loss_function)
            model = mp.createModel()            
            model.summary()  
            
        elif cfg.model == 'BreastSegmentor_v1':
            from BreastMask_Models_Collection import BreastSegmentor_v1
            mp = BreastSegmentor_v1(cfg.output_classes, cfg.num_channels, cfg.L2, cfg.dropout, cfg.learning_rate, cfg.optimizer_decay, cfg.loss_function)
            model = mp.createModel()            
            model.summary()              

        elif cfg.model == 'UNet_3D_v0':
            from Unet_3D_Class import Unet_3D
            model = Unet_3D().create_model0((12,76,76,2), pool_size=(2, 2, 2), n_labels=2, initial_learning_rate=0.00001, deconvolution=True, depth=4, n_base_filters=32, include_label_wise_dice_coefficients=False, batch_normalization=True, activation_name="softmax")     
            model.summary()          
        elif cfg.model == 'UNet_3D_v1':
            from Unet_3D_Class import Unet_3D
            model = Unet_3D().create_model1((12,76,76,2), pool_size=(2, 2, 2), n_labels=2, initial_learning_rate=0.00001, deconvolution=True, depth=4, n_base_filters=32, include_label_wise_dice_coefficients=False, batch_normalization=True, activation_name="softmax")     
            model.summary()     
            
        elif cfg.model == 'UNet_3D_v4':
            from Unet_3D_Class import UNet_v4
            model = UNet_v4(input_shape=(None,None,None,3), pool_size=(2, 2, 2), n_labels=2, initial_learning_rate=cfg.learning_rate, deconvolution=True,
                  depth=4, n_base_filters=48, include_label_wise_dice_coefficients=False, 
                  batch_normalization=True, activation_name="softmax", bilinear_upsampling=True)
            model.summary()                



        elif cfg.model == 'UNet_v0_BreastMask':
            from Unet_3D_Class import UNet_v0_BreastMask
            model = UNet_v0_BreastMask(input_shape =  (3, 256,256,1), pool_size=(1, 2, 2), n_labels=1, initial_learning_rate=0.0001, deconvolution=False,
                      depth=4, n_base_filters=32, include_label_wise_dice_coefficients=False, 
                      batch_normalization=False, add_spatial_coordinates=cfg.use_coordinates)
            model.summary()    
            
        elif cfg.model == 'UNet_v4_BreastMask':
            from Unet_3D_Class import UNet_v4_BreastMask
            model = UNet_v4_BreastMask(input_shape=(19,75,75,1), pool_size=(2, 2, 2), n_labels=2, initial_learning_rate=0.00001, deconvolution=True,
                  depth=4, n_base_filters=32, include_label_wise_dice_coefficients=False, 
                  batch_normalization=True, activation_name="softmax", bilinear_upsampling=True)
            model.summary()    

        elif cfg.model == 'UNet_v4_TPM':
            from Unet_3D_Class import UNet_v4_TPM
            model = UNet_v4_TPM(input_shape=(19,75,75,4), pool_size=(2, 2, 2), n_labels=2, initial_learning_rate=0.00001, deconvolution=True,
                  depth=4, n_base_filters=32, include_label_wise_dice_coefficients=False, 
                  batch_normalization=True, activation_name="softmax", bilinear_upsampling=True)
            model.summary()    

        elif cfg.model == 'UNet_v0_TumorSegmenter':
            from Unet_3D_Class import UNet_v0_TumorSegmenter
            model = UNet_v0_TumorSegmenter(input_shape =  (3, 512, 512,4), pool_size=(1, 2, 2), n_labels=1, initial_learning_rate=0.00001, deconvolution=False,
                      depth=6, n_base_filters=16, include_label_wise_dice_coefficients=False, 
                      batch_normalization=True, activation_name="softmax")
            model.summary()  

            
        else: 
            print('ERROR: No model selected.')
            return 0          
        
        
        if cfg.merge_breastMask_model:

            from keras.models import load_model  
            from MultiPriors_Models_Collection import Generalised_dice_coef_multilabel2, dice_coef_multilabel0,dice_coef_multilabel1
            my_custom_objects = {'Generalised_dice_coef_multilabel2':Generalised_dice_coef_multilabel2,
                                 'dice_coef_multilabel0':dice_coef_multilabel0,
                                 'dice_coef_multilabel1':dice_coef_multilabel1}            
            bm_model = load_model(cfg.path_to_breastMask_model, custom_objects = my_custom_objects )            
            preTrained_convLayers = [x for x in bm_model.layers if 'Conv3D' in str(x)]
            preTrained_batchNormLayers = [x for x in bm_model.layers if 'BatchNorm' in str(x)]
            if cfg.Context_parameters_trainable:
                # If fine-tuning allowed, then skip bottleneck of last layer of the breastMask model 
                preTrained_convLayers = preTrained_convLayers[:-1]
                
            #preTrained_convLayers = [x for x in bm_model.layers if 'T1post_Context' in x.name]
            #preTrained_batchNormLayers = [x for x in bm_model.layers if 'BatchNorm' in x.name]   
            newModel_convLayers = [x for x in model.layers if 'T1post_Context' in x.name]           
            newModel_batchNormLayers = [x for x in model.layers if 'BatchNorm' in x.name]            
            
            assert len(preTrained_convLayers ) == len(newModel_convLayers), 'Models have incompatible architecture..'
            assert len(preTrained_batchNormLayers) == len(newModel_batchNormLayers), 'Models have incompatible architecture..'
            print('Transfering weights from breastMask model {} '.format(cfg.path_to_breastMask_model))       
#            for i in range(len(newModel_convLayers)):
#                print('Equal layer {}: {}'.format(newModel_convLayers[i].name, (model.get_layer(newModel_convLayers[i].name).get_weights()[0] == preTrained_convLayers[i].get_weights()[0]).all()))
            for i in range(len(newModel_convLayers)):
                print('Weight transfer of layer : {}'.format(newModel_convLayers[i].name))
                model.get_layer(newModel_convLayers[i].name).set_weights(preTrained_convLayers[i].get_weights())        
                model.get_layer(newModel_convLayers[i].name).trainable = cfg.Context_parameters_trainable
            for i in range(len(newModel_batchNormLayers)):
                print('Weight transfer of layer : {}'.format(newModel_batchNormLayers[i].name))
                model.get_layer(newModel_batchNormLayers[i].name).set_weights(preTrained_batchNormLayers[i].get_weights())        
                
            # Need to re-compile when changing the TRAINABLE attribute:   
            #model = multi_gpu_model(model, gpus=4)    
            from MultiPriors_Models_Collection import Generalised_dice_coef_multilabel2, dice_coef_multilabel_bin0, dice_coef_multilabel_bin1
            model.compile(loss=Generalised_dice_coef_multilabel2, optimizer=keras.optimizers.adam(lr=cfg.learning_rate), 
                          metrics=['acc', dice_coef_multilabel_bin0, dice_coef_multilabel_bin1])
            model.summary()     
                    
        start_epoch = 0
        os.chdir(workingDir + '/training_sessions/')
        session = cfg.model + '_' + cfg.dataset + '_' + configFileName + '_' + time.strftime("%Y-%m-%d_%H%M") 
        wd = workingDir + '/training_sessions/' +session
        if not os.path.exists(wd):    
            os.mkdir(session)
            os.mkdir(session + '/models')
            os.mkdir(session + '/predictions')
        os.chdir(wd) 

        CV_FOLDS_ARRAYS_PATH = '/home/deeperthought/Projects/MultiPriors_MSKCC/CV_folds/CV_alignedNii-Aug2019_actual-F4-training/arrays/'#'/'.join(cfg.trainChannels[0].split('/')[:-1]) + '/arrays/'
        if not os.path.exists(CV_FOLDS_ARRAYS_PATH):
            os.mkdir(CV_FOLDS_ARRAYS_PATH)
        
        #copy(workingDir + configFile[1:], wd)
        copy(configFile, wd)
        logfile = session +'.log'            
        print(model.summary())
        val_performance = []
        from keras.utils import plot_model
        plot_model(model, to_file=wd+'/multiscale_TPM.png', show_shapes=True)
        with open(wd+'/model_summary.txt','w') as fh:
            # Pass the file handle in as a lambda function to make it callable
            model.summary(print_fn=lambda x: fh.write(x + '\n'))
        if len(cfg.comments) > 0:
            f = open('Comments.txt','w')
            f.write(str(cfg.comments))
            f.close()
        
    elif cfg.load_model == True:
        from keras.models import load_model  
        if cfg.loss_function == 'Dice6':
            from MultiPriors_Models_Collection import dice_coef_multilabel6, dice_coef_multilabel0,dice_coef_multilabel1,dice_coef_multilabel2,dice_coef_multilabel3,dice_coef_multilabel4,dice_coef_multilabel5
            my_custom_objects = {'dice_coef_multilabel6':dice_coef_multilabel6,
                                 'dice_coef_multilabel0':dice_coef_multilabel0,
                                 'dice_coef_multilabel1':dice_coef_multilabel1,
                                 'dice_coef_multilabel2':dice_coef_multilabel2,
                                 'dice_coef_multilabel3':dice_coef_multilabel3,
                                 'dice_coef_multilabel4':dice_coef_multilabel4,
                                 'dice_coef_multilabel5':dice_coef_multilabel5}
        elif cfg.loss_function == 'wDice6':
            from MultiPriors_Models_Collection import w_dice_coef_multilabel6, dice_coef_multilabel0,dice_coef_multilabel1,dice_coef_multilabel2,dice_coef_multilabel3,dice_coef_multilabel4,dice_coef_multilabel5
            my_custom_objects = {'w_dice_coef_multilabel6':w_dice_coef_multilabel6,
                                 'dice_coef_multilabel0':dice_coef_multilabel0,
                                 'dice_coef_multilabel1':dice_coef_multilabel1,
                                 'dice_coef_multilabel2':dice_coef_multilabel2,
                                 'dice_coef_multilabel3':dice_coef_multilabel3,
                                 'dice_coef_multilabel4':dice_coef_multilabel4,
                                 'dice_coef_multilabel5':dice_coef_multilabel5}
        elif cfg.loss_function == 'Dice':
            from MultiPriors_Models_Collection import Generalised_dice_coef_multilabel2, dice_coef_multilabel_bin0,dice_coef_multilabel_bin1
            my_custom_objects = {'Generalised_dice_coef_multilabel2':Generalised_dice_coef_multilabel2,
                                 'dice_coef_multilabel_bin0':dice_coef_multilabel_bin0,
                                 'dice_coef_multilabel_bin1':dice_coef_multilabel_bin1}
            
        elif cfg.loss_function == 'ReLU_Dice':
            from MultiPriors_Models_Collection import Generalised_dice_coef_multilabel2, dice_coef_multilabel0,dice_coef_multilabel1
            my_custom_objects = {'Generalised_dice_coef_multilabel2':Generalised_dice_coef_multilabel2,
                                 'dice_coef_multilabel0':dice_coef_multilabel0,
                                 'dice_coef_multilabel1':dice_coef_multilabel1,
                                 'RAdam':RAdam}            
        elif cfg.loss_function == 'wDice2':
            from DM_MSKCC_Atrous_model import w_dice_coef_multilabel2, dice_coef_multilabel0,dice_coef_multilabel1
            my_custom_objects = {'w_dice_coef_multilabel2':w_dice_coef_multilabel2,
                                 'dice_coef_multilabel0':dice_coef_multilabel0,
                                 'dice_coef_multilabel1':dice_coef_multilabel1}
        model = load_model(cfg.path_to_model, custom_objects = my_custom_objects )
        print('LOADED MODEL FROM SESSION {}'.format(cfg.session))
        session = cfg.session
        start_epoch = 0 #int(cfg.path_to_model.split('.')[-2][cfg.path_to_model.split('.')[-2].find('epoch') + 5 : ]) + 1
        #cfg.epochs_for_fullSegmentation = range(start_epoch+1, cfg.epochs)
        os.chdir(workingDir + '/training_sessions/')
        wd = workingDir + '/training_sessions/' +session
        if not os.path.exists(wd):    
            os.mkdir(session)
            os.mkdir(session + '/models')
            os.mkdir(session + '/predictions')
        os.chdir(wd)
        logfile = session +'.log'
        CV_FOLDS_ARRAYS_PATH = '/home/deeperthought/Projects/MultiPriors_MSKCC/CV_folds/CV_alignedNii-Aug2019_actual-F4-training/arrays/'#'/'.join(cfg.trainChannels[0].split('/')[:-1]) + '/arrays/'
        if not os.path.exists(CV_FOLDS_ARRAYS_PATH):
            os.mkdir(CV_FOLDS_ARRAYS_PATH)
    #################################################################################################
    #                                                                                               #
    #                                         START SESSION                                         #
    #                                                                                               #
    #################################################################################################
    
    val_performance = []
    full_segm_DICE = []
    full_segm_SMOOTH_DICE = []
    epoch_foreground_percent = []
    losses = []
    metrics = []
    EARLY_STOP = False      



    start_training_session_logger(logfile, cfg.threshold_EARLY_STOP, cfg.TPM_channel, cfg.load_model, cfg.saveSegmentation, cfg.path_to_model, model, \
        cfg.dropout, cfg.trainChannels, cfg.trainLabels, cfg.validationChannels, cfg.validationLabels, \
        cfg.testChannels, cfg.testLabels, cfg.num_iter, cfg.epochs, cfg.n_patches, cfg.n_patches_val, cfg.n_subjects, cfg.samplingMethod_train, \
        cfg.size_minibatches, cfg.n_full_segmentations, cfg.epochs_for_fullSegmentation, cfg.size_test_minibatches)
    # Callback history    
    if cfg.output_classes == 2:
        history = LossHistory_multiDice2() 
    elif cfg.output_classes == 6:
        history = LossHistory_multiDice6()

    ############## START SAMPLING DATA FOR FIRST TRAINING AND VALIDATION STEP ###################### 
    manager = multiprocessing.Manager()
    return_dict = manager.dict()
    train_data_process = Process(target=sampleTrainData_daemon, args=(return_dict, 'TRAINING', cfg.resolution, cfg.trainChannels, CV_FOLDS_ARRAYS_PATH, cfg.trainLabels,cfg.TPM_channel, 
                                                                      cfg.n_patches, cfg.n_subjects, cfg.dpatch, cfg.output_classes, 
                                                                      cfg.samplingMethod_train, cfg.use_coordinates, cfg.proportion_malignants_to_sample_train, 
                                                                      cfg.percentile_voxel_intensity_sample_benigns,  cfg.data_augmentation, cfg.proportion_to_flip, 
                                                                      cfg.intensity_normalization_method, cfg.MRI_PATH, cfg.model_patch_reduction, cfg.model_crop, cfg.balanced_sample_subjects))
    train_data_process.start()
    val_data_process = Process(target=sampleTrainData_daemon, args=(return_dict, 'VALIDATION', cfg.resolution, cfg.validationChannels, CV_FOLDS_ARRAYS_PATH, cfg.validationLabels, cfg.TPM_channel, cfg.n_patches_val, 
                                                            cfg.n_subjects_val, cfg.dpatch, cfg.output_classes, cfg.samplingMethod_val, cfg.use_coordinates, 
                                                            cfg.proportion_malignants_to_sample_val, cfg.percentile_voxel_intensity_sample_benigns, 
                                                            cfg.data_augmentation, cfg.proportion_to_flip, cfg.intensity_normalization_method, cfg.MRI_PATH, 
                                                            cfg.model_patch_reduction, cfg.model_crop, cfg.balanced_sample_subjects))

    val_data_process.start()

    for epoch in xrange(start_epoch,cfg.epochs):
      t1 = time.time()
      my_logger("######################################################",logfile)
      my_logger("                   TRAINING EPOCH " + str(epoch) + "/" + str(cfg.epochs),logfile)
      my_logger("######################################################",logfile)
   
      ####################### FULL HEAD SEGMENTATION ##############################
                
      if epoch in cfg.epochs_for_fullSegmentation:
        my_logger("------------------------------------------------------", logfile)
        my_logger("                 FULL HEAD SEGMENTATION", logfile)
        my_logger("------------------------------------------------------", logfile)
        dice_compare = True
        dsc = []
        smooth_dice_scores = []
        foreground_percent_list = []
        subjectIndex = 0

        with open(cfg.validationLabels) as vl:
            n_valSubjects = len(vl.readlines())
        if cfg.test_subjects > n_valSubjects:
            print("Given number of subjects for test set (" + str(cfg.test_subjects) +") is larger than the amount of \
            subjects in test set (" +str(n_valSubjects)+ ")")
            cfg.test_subjects = n_valSubjects
            cfg.n_full_segmentations = n_valSubjects
            print('Using {} number of test subjects'.format(n_valSubjects))
        if len(cfg.list_subjects_fullSegmentation) == 0:
            #list_subjects_fullSegmentation = sample(range(cfg.test_subjects), cfg.n_full_segmentations)
            if cfg.balanced_sample_subjects:
                proportion_malignants = int(np.ceil(cfg.n_full_segmentations * cfg.proportion_malignants_fullSegmentation))
                labelsFile = [x[:-1] for x in open(cfg.testLabels).readlines()]
                malignant_subjects_index = [labelsFile.index(x) for x in labelsFile if not 'BENIGN' in x]
                benign_subjects_index = list(set(range(len(labelsFile))) - set(malignant_subjects_index))
                list_subjects_fullSegmentation = random.sample(malignant_subjects_index, min(len(malignant_subjects_index), proportion_malignants))
                print('sampling {} malignants from partition'.format(len(list_subjects_fullSegmentation)))
                list_subjects_fullSegmentation.extend(random.sample(benign_subjects_index, cfg.n_full_segmentations - len(list_subjects_fullSegmentation)))
                random.shuffle(list_subjects_fullSegmentation)
            else:
                list_subjects_fullSegmentation = random.sample(xrange(cfg.test_subjects ), cfg.n_full_segmentations)

        else:
            list_subjects_fullSegmentation = cfg.list_subjects_fullSegmentation
               
        for subjectIndex in list_subjects_fullSegmentation: 
            t_segment = time.time()
            if cfg.full_segmentation_patches:
                segmentation_img, output, dsc, score_smooth,foreground_percent = fullSegmentation(wd, cfg.penalty_MATRIX,cfg.resolution, cfg.OUTPUT_PATH, cfg.TPM_channel, 
                                                                                                  dice_compare, dsc, smooth_dice_scores, foreground_percent_list, 
                                                                                                  model, cfg.testChannels, cfg.testLabels, subjectIndex, 
                                                                                                  cfg.output_classes, cfg.segmentation_dpatch, cfg.size_test_minibatches,
                                                                                                  cfg.output_probability, cfg.use_coordinates, cfg.intensity_normalization_method, cfg.MRI_PATH,
                                                                                                  cfg.model_patch_reduction, cfg.model_crop, epoch, using_Unet=True, 
                                                                                                  using_unet_breastMask=cfg.using_unet_breastMask)
            else:
                fullSegmentation_Flexible(wd, cfg.penalty_MATRIX, dice_compare, dsc, model, cfg.testChannels, cfg.testLabels, cfg.TPM_channel, subjectIndex, \
                cfg.output_classes, cfg.dpatch, cfg.size_test_minibatches, logfile, epoch, cfg.use_coordinates,cfg.saveSegmentation)

            my_logger('--------------- TEST EVALUATION ---------------', logfile)
            my_logger('          Full segmentation evaluation of subject' + str(subjectIndex), logfile)
            my_logger('foreground_percent {}'.format(foreground_percent), logfile)
            my_logger('SMOOTH_DCS ' + str(score_smooth),logfile)  
            if dsc != None:
                if len(dsc) > 0:
                    my_logger('DCS ' + str(dsc[-1]),logfile)
            if cfg.saveSegmentation:
                nib.save(segmentation_img, output)
                my_logger('Saved segmentation of subject at: ' + output, logfile)    
                
           #if (dice_compare & len(dsc)>0):
            #  my_logger('DCS ' + str(dsc[-1]),logfile)
            print('Segmentation of subject took {} s'.format(time.time()-t_segment))
        if dsc != None:
            my_logger('         FULL SEGMENTATION SUMMARY STATISTICS ', logfile)
            full_segm_DICE.append(np.mean(dsc))   
            full_segm_SMOOTH_DICE.append(np.mean(smooth_dice_scores))   
            my_logger('Overall DCS:   ' + str(full_segm_DICE[-1]),logfile)
            my_logger('Overall SMOOTH_DCS:   ' + str(full_segm_SMOOTH_DICE[-1]),logfile)
            epoch_foreground_percent.append(np.mean(foreground_percent_list))            
            my_logger('Epoch_foreground_percent {}'.format(epoch_foreground_percent[-1]), logfile)

        # Function to define if STOP flag goes to True or not, based on difference between last three or two segmentations.
        if len(full_segm_DICE) > 5:                        
            if np.max(np.abs(np.diff([full_segm_DICE[-3], full_segm_DICE[-2], full_segm_DICE[-1]] ))) < cfg.threshold_EARLY_STOP:
                EARLY_STOP = True
            #elif np.max(np.abs(np.diff([full_segm_DICE[-5],full_segm_DICE[-4],full_segm_DICE[-3], full_segm_DICE[-2], full_segm_DICE[-1]] ))) < 0.03:
            #    EARLY_STOP = True

        if len(full_segm_DICE) == 1:
            my_logger('###### SAVING TRAINED MODEL AT : ' + wd + '/models/best_model.h5', logfile)
            model.save(wd+'/models/best_model.h5')            
        # Save model if best results achieved
        if len(full_segm_DICE) > 1:
          if np.max(full_segm_DICE[:-1]) <= full_segm_DICE[-1]:   # If last results is larger than any other before, then this is the best current model.
            #my_logger('###### SAVING TRAINED MODEL AT : ' + wd +'/Output/models/'+logfile[12:]+'.h5', logfile)
            #model.save(wd+'/models/'+logfile[12:]+'_epoch' + str(epoch) + '.h5')
            my_logger('###### SAVING TRAINED MODEL AT : ' + wd + '/models/best_model.h5', logfile)
            model.save(wd+'/models/best_model.h5')

        if EARLY_STOP:
          my_logger('Convergence criterium met. Stopping training.',logfile)
          break           
      #################################################################################################
      #                                                                                               #
      #                               Training and Validation                                         #
      #                                                                                               #
      #################################################################################################
      if cfg.sample_intensity_based:
#          if epoch > 20:
#            cfg.percentile_voxel_intensity_sample_benigns = 99   # This will sample only from VERY high intensity areas          
#            my_logger('Reached epoch {}, changing sampling benign voxels to voxel intensity > {}'.format(epoch, cfg.percentile_voxel_intensity_sample_benigns ), logfile)             

          if epoch > 5:
            cfg.percentile_voxel_intensity_sample_benigns = 90   # This will sample only from high intensity areas
            my_logger('Reached epoch {}, changing sampling benign voxels to voxel intensity > {}'.format(epoch, cfg.percentile_voxel_intensity_sample_benigns ), logfile)          
       
#          elif epoch > 5:
#            cfg.percentile_voxel_intensity_sample_benigns = 90   # This will sample only from breast
#            my_logger('Reached epoch {}, changing sampling benign voxels to voxel intensity > {}'.format(epoch, cfg.percentile_voxel_intensity_sample_benigns ), logfile)

      ####################################### PARALLELIZED TRAINING ITERATIONS ########################################
    
      for i in range(0, cfg.num_iter):
        print("###################################################### ")          
        my_logger("                   Batch " + str(i+1) + "/" + str(cfg.num_iter) ,logfile)
        print("------------------------------------------------------ ")          
                    
        ####################### VALIDATION ON BATCHES ############################      
        print('\n###################### VALIDATION ####################')                                      
        with open(cfg.validationLabels) as vl:
            n_valSubjects = len(vl.readlines())
        if cfg.n_subjects_val > n_valSubjects:
            print("Given number of subjects for test set (" + str(cfg.n_subjects_val) +") is larger than the amount of subjects in test set (" +str(n_valSubjects)+ ")")
            cfg.n_subjects_val = n_valSubjects
            print('Using {} number of test subjects'.format(n_valSubjects))
    
        #--------- VALIDATION ---------
        # Make sure data sampling finishes    
        val_data_process.join()    
        # Start next data sampling process
        val_data_process = Process(target=sampleTrainData_daemon, args=(return_dict, 'VALIDATION',cfg.resolution, cfg.validationChannels, CV_FOLDS_ARRAYS_PATH, cfg.validationLabels, cfg.TPM_channel, cfg.n_patches_val, 
                                                                            cfg.n_subjects_val, cfg.dpatch, cfg.output_classes, cfg.samplingMethod_val, cfg.use_coordinates, 
                                                                            cfg.proportion_malignants_to_sample_val, cfg.percentile_voxel_intensity_sample_benigns, 
                                                                            cfg.data_augmentation, cfg.proportion_to_flip, cfg.intensity_normalization_method, cfg.MRI_PATH,
                                                                            cfg.model_patch_reduction, cfg.model_crop, cfg.balanced_sample_subjects))
        val_data_process.start()
        # Validate model
        val_performance.append(train_validate_model_on_batch(cfg.model, model, return_dict['VALIDATION'][0],  return_dict['VALIDATION'][1], return_dict['VALIDATION'][2], 
                                                             return_dict['VALIDATION'][3], return_dict['VALIDATION'][4], cfg.size_minibatches_val, history, losses,  metrics, 
                                                               cfg.output_classes, logfile, TRAINING_FLAG=False, using_unet_breastMask=cfg.using_unet_breastMask)) 
        #--------- TRAINING ---------
        print('\n###################### TRAINING ####################')                                                   
        with open(cfg.trainLabels) as vl:
            n_trainSubjects = len(vl.readlines())                
        if cfg.n_subjects > n_trainSubjects:
            print("Given number of subjects for test set (" + str(cfg.n_subjects) +") is larger than the amount of \
            subjects in test set (" +str(n_trainSubjects)+ ")")
            cfg.n_subjects = n_trainSubjects
            print('Using {} number of test subjects'.format(n_trainSubjects))
        print('sampling {} patches'.format(cfg.n_patches))    
        # Make sure data sampling finishes
        train_data_process.join()    

        # Start next data sampling process
        train_data_process = Process(target=sampleTrainData_daemon, args=(return_dict, 'TRAINING', cfg.resolution,cfg.trainChannels, CV_FOLDS_ARRAYS_PATH, cfg.trainLabels,cfg.TPM_channel, 
                                                                      cfg.n_patches, cfg.n_subjects, cfg.dpatch, cfg.output_classes, 
                                                                      cfg.samplingMethod_train, cfg.use_coordinates, cfg.proportion_malignants_to_sample_train, 
                                                                      cfg.percentile_voxel_intensity_sample_benigns,  cfg.data_augmentation, cfg.proportion_to_flip, 
                                                                      cfg.intensity_normalization_method, cfg.MRI_PATH, cfg.model_patch_reduction, cfg.model_crop, cfg.balanced_sample_subjects))
        train_data_process.start()
        
               
        context = return_dict['TRAINING'][0]
        patches = return_dict['TRAINING'][1]
        target_labels = return_dict['TRAINING'][2]
        spatial_coordinates = return_dict['TRAINING'][3]
        TPM_patches = return_dict['TRAINING'][4] 
        
        train_data_process.join()    
        train_data_process = Process(target=sampleTrainData_daemon, args=(return_dict, 'TRAINING', cfg.resolution, cfg.trainChannels, CV_FOLDS_ARRAYS_PATH, cfg.trainLabels,cfg.TPM_channel, 
                                                                          cfg.n_patches, cfg.n_subjects, cfg.dpatch, cfg.output_classes, 
                                                                          cfg.samplingMethod_train, cfg.use_coordinates, cfg.proportion_malignants_to_sample_train, 
                                                                          cfg.percentile_voxel_intensity_sample_benigns,  cfg.data_augmentation, cfg.proportion_to_flip, 
                                                                          cfg.intensity_normalization_method, cfg.MRI_PATH, cfg.model_patch_reduction, cfg.balanced_sample_subjects))       
        train_data_process.start() 
        
        context = np.concatenate([context, return_dict['TRAINING'][0]])
        patches = np.concatenate([patches, return_dict['TRAINING'][1]])
        target_labels = np.concatenate([target_labels, return_dict['TRAINING'][2]])
        spatial_coordinates = np.concatenate([spatial_coordinates, return_dict['TRAINING'][3]])     
        TPM_patches = np.concatenate([TPM_patches, return_dict['TRAINING'][4]])
#        
#        train_data_process = Process(target=sampleTrainData_daemon, args=(return_dict, 'TRAINING', cfg.resolution,cfg.trainChannels, CV_FOLDS_ARRAYS_PATH, cfg.trainLabels,cfg.TPM_channel, 
#                                                                      cfg.n_patches, cfg.n_subjects, cfg.dpatch, cfg.output_classes, 
#                                                                      cfg.samplingMethod_train, cfg.use_coordinates, cfg.proportion_malignants_to_sample_train, 
#                                                                      cfg.percentile_voxel_intensity_sample_benigns,  cfg.data_augmentation, cfg.proportion_to_flip, 
#                                                                      cfg.intensity_normalization_method, cfg.MRI_PATH, cfg.model_patch_reduction, cfg.model_crop, cfg.balanced_sample_subjects))
#        train_data_process.start() 
       
        print('patches.shape = {}'.format(patches.shape))
        print('target_labels.shape = {}'.format(target_labels.shape))
        #print('spatial_coordinates.shape = {}'.format(spatial_coordinates.shape))

        # Train model
        epoch_loss, epoch_metrics = train_validate_model_on_batch(cfg.model, model,context,patches,
                                                              target_labels, spatial_coordinates, TPM_patches,
                                                              cfg.size_minibatches,history,losses,metrics,cfg.output_classes, logfile, using_unet_breastMask=cfg.using_unet_breastMask)  
        # For large datasets, save model after every weight updates
        print('Saving model..')
        model.save(wd+'/models/last_model.h5')        
        try:
          global_loss = np.concatenate([np.load(wd + '/LOSS.npy', allow_pickle=True), epoch_loss])
          global_metrics = np.concatenate([np.load(wd + '/METRICS.npy', allow_pickle=True), epoch_metrics]) 
          np.save(wd + '/LOSS.npy', global_loss)
          np.save(wd + '/METRICS.npy', global_metrics)
        except:
          np.save(wd + '/LOSS.npy', epoch_loss)
          np.save(wd + '/METRICS.npy', epoch_metrics)              


      my_logger('Total training this epoch took ' + str(round(time.time()-t1,2)) + ' seconds',logfile)
      #############################################################################################################################

#    plot_training(session,losses, metrics, val_performance, full_segm_DICE, smooth=20, loss_name = [cfg.loss_function], class_names = [str(x) for x in range(cfg.output_classes)])
#    plt.savefig(wd + '/' + session + '.png')
#    plt.close()




