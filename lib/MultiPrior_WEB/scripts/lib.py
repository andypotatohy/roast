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
import time
import random
from numpy.random import seed
import keras 
from keras import backend as K
from keras.utils import to_categorical
import tensorflow 

#from sklearn import metrics
import matplotlib.pyplot as plt
from shutil import copyfile
from random import shuffle
from random import sample
import pandas as pd

import multiprocessing        
from multiprocessing import Pool
from multiprocessing.dummy import Pool as ThreadPool      

from skimage.transform import resize

seed(1)
tensorflow.random.set_seed(2)

#import tensorflow as tf
#config = tf.ConfigProto()
#config.gpu_options.allow_growth = True
#tf.keras.backend.set_session(tf.Session(config=config))


############################## METRICS FUNCTIONS ##############################################


def weighted_generalized_dice_completeImages(img1,img2,penalty_MATRIX):
    assert img1.shape == img2.shape, 'Images of different size!'
    classes = np.array(np.unique(img1), dtype='int8')   
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
    dice = []
    for i in classes:
        dice.append(2*np.sum(np.multiply(img1==i,img2==i))/float(np.sum(img1==i)+np.sum(img2==i)))   
    return np.sum(dice)/len(classes), [round(x,2) for x in dice]
    
    
def percentile95_normalizeMRI(data):
    p95 = np.percentile(data,95)
    data = data/p95
    return(data)    
    
def normalizeMRI(data, mean=0, std=1):
    if (mean == 0) and (std == 1):
        mean = np.mean(data)
        std = np.std(data)
    data = (data - mean)/std
    return(data)

################################ SAMPLING FUNCTIONS ####################################################

def flip_random(patches, labels, TPM_patches, proportion_to_flip=0.5):
  indx_toflip = np.random.choice(np.arange(0,patches.shape[0]), int(patches.shape[0]*proportion_to_flip), replace=False)
  axis = np.random.choice(range(0,3),size=len(indx_toflip))
  for i in range(len(indx_toflip)):
    if axis[i] == 0:    
      # SAGITTAL FLIP
      patches[indx_toflip[i],:,:,:,0] = patches[indx_toflip[i],::-1,:,:,0]
      labels[indx_toflip[i],:,:,:] = labels[indx_toflip[i],::-1,:,:]
      if len(TPM_patches) > 0:
        for ii in range(0,TPM_patches.shape[-1]):
          TPM_patches[indx_toflip[i],:,:,:,ii] = TPM_patches[indx_toflip[i],::-1,:,:,ii]
    elif axis[i] == 1:        
      # AXIAL FLIP      
      patches[indx_toflip[i],:,:,:,0] = patches[indx_toflip[i],:,::-1,:,0]
      labels[indx_toflip[i],:,:,:] = labels[indx_toflip[i],:,::-1,:]
      if len(TPM_patches) > 0:
        for ii in range(0,TPM_patches.shape[-1]):
          TPM_patches[indx_toflip[i],:,:,:,ii] = TPM_patches[indx_toflip[i],:,::-1,:,ii]
    elif axis[i] == 2:      
      # CORONAL FLIP
      patches[indx_toflip[i],:,:,:,0] = patches[indx_toflip[i],:,:,::-1,0]
      labels[indx_toflip[i],:,:,:] = labels[indx_toflip[i],:,:,::-1]
      if len(TPM_patches) > 0:
        for ii in range(0,TPM_patches.shape[-1]):
          TPM_patches[indx_toflip[i],:,:,:,ii] = TPM_patches[indx_toflip[i],:,:,::-1,ii]
        
    return patches, labels, TPM_patches    

def generateRandomIndexesSubjects(n_subjects, total_subjects):
    indexSubjects = random.sample(range(total_subjects), n_subjects)
    return indexSubjects

def getSubjectChannels(subjectIndexes, channel):
    "With the channels (any modality) and the indexes of the selected subjects, return the addresses of the subjects channels"
    fp = open(channel)
    # read file, per subject index extract patches given the indexesPatch
    lines = fp.readlines()
    selectedSubjects = [lines[i][:-1] for i in subjectIndexes]
    fp.close()
    return selectedSubjects

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
        shapes.append(proxy_img.shape)
    return shapes      

def getSubjectShapes_parallelization(subjectChannel):
    # Need to open every nifty file and get the shapes
    proxy_img = nib.load(subjectChannel)
    shape = proxy_img.shape
    proxy_img.uncache()
    del proxy_img
    return shape


def generateVoxelIndexes(subjectIndexes, shapes, patches_per_subject, dpatch, n_patches, groundTruthChannel_list, samplingMethod, output_classes, allForegroundVoxels = ""):
    "Alternative improved function of the same named one above."
    "Here extract the channels from the subject indexes, and loop over them. Then in second loop extract as many needed voxel coordinates per subject."
    methods =["Random sampling","Equal sampling background/foreground","Equal sampling all classes center voxel","Equal sampling background/foreground with exhaustive foreground samples"]
    print("Generating voxel indexes with method: " + methods[samplingMethod])
    channels = getSubjectChannels(subjectIndexes, groundTruthChannel_list)
    allVoxelIndexes = []    
    if samplingMethod == 0:
        for i in range(0, len(shapes)):
            voxelIndexesSubj = []
            #loop over voxels per subject
            for j in range(0,patches_per_subject[i]):
                # unform sampling
                # central voxel of final 9x9x9 label cube
                voxelIndexesSubj.append((np.random.randint(4, shapes[i][0]-5),np.random.randint(4, shapes[i][1]-5),np.random.randint(4, shapes[i][2]-5)))
               
               
#                voxelIndexesSubj.append((np.random.randint(0+dpatch/2, shapes[i][0]-(dpatch/2)-1),
#                                         np.random.randint(0+dpatch/2, shapes[i][1]-(dpatch/2)-1),
#                                         np.random.randint(0+dpatch/2, shapes[i][2]-(dpatch/2)-1)))
            allVoxelIndexes.append(voxelIndexesSubj)            
        random.shuffle(allVoxelIndexes[i])
        return allVoxelIndexes    
    elif samplingMethod == 1:
        "This samples equally background/foreground. Assumption that foreground is very seldom: Only foreground voxels are sampled, and background voxels are just random samples which are then proofed against foreground ones"
        "Still need to proof for repetition. Although unlikely and uncommon"
        for i in range(0,len(channels)): 
            voxelIndexesSubj = []
            backgroundVoxels = []           
            fg = getForegroundBackgroundVoxels(channels[i], dpatch) # This function returns only foreground voxels
            #print("Extracting foreground voxels " + str(patches_per_subject[i]/2 + patches_per_subject[i]%2) +' from ' + str(len(fg)) + " from channel " + str(channels[i]) + " with index " + str(i))
            foregroundVoxels = fg[random.sample(range(0,len(fg)), min(len(fg),patches_per_subject[i]/2 + patches_per_subject[i]%2))].tolist()            
            # get random voxel coordinates
            for j in range(0,patches_per_subject[i]/2):
                backgroundVoxels.append((np.random.randint(0+dpatch/2, shapes[i][0]-(dpatch/2)-1),np.random.randint(0+dpatch/2, shapes[i][1]-(dpatch/2)-1),np.random.randint(0+dpatch/2, shapes[i][2]-(dpatch/2)-1)))                
            #Replace the ones that by chance are foreground voxels (not so many in tumor data)
            while any([e for e in foregroundVoxels if e in backgroundVoxels]):
                ix = [e for e in foregroundVoxels if e in backgroundVoxels]
                for index in ix:
                    newVoxel = [np.random.randint(dpatch/2, shapes[i][0]-(dpatch/2)-1),np.random.randint(dpatch/2, shapes[i][1]-(dpatch/2)-1),np.random.randint(dpatch/2, shapes[i][2]-(dpatch/2)-1)]
                    backgroundVoxels[backgroundVoxels.index(index)] = newVoxel
            allVoxelIndexes.append(foregroundVoxels + backgroundVoxels)
            random.shuffle(allVoxelIndexes[i])
        return allVoxelIndexes    
    elif samplingMethod == 2:
        "sample from each class equally"        
        "use function getAllForegroundClassesVoxels to get coordinates from all classes (not including background)"
        for i in range(0,len(channels)):  # iteration over subjects
            voxelIndexesSubj = []
            backgroundVoxels = []
            fg = getAllForegroundClassesVoxels(channels[i], dpatch, output_classes) # This function returns only foreground voxels            
            # WATCH OUT, FG IS A LIST OF LISTS. fIRST DIMENSION IS THE CLASS, SECOND IS THE LIST OF VOXELS OF THAT CLASS
            foregroundVoxels = []
            patches_to_sample = [patches_per_subject[i]/output_classes] * output_classes #
            extra = random.sample(range(output_classes),1)
            patches_to_sample[extra[0]] = patches_to_sample[extra[0]] + patches_per_subject[i]%output_classes            
            for c in range(0, output_classes-1):
                foregroundVoxels.extend(fg[c][random.sample(range(0,len(fg[c])), min(patches_to_sample[c],len(fg[c])))].tolist())
            # get random voxel coordinates
            for j in range(0,patches_per_subject[i]/output_classes):
                backgroundVoxels.append([np.random.randint(0+dpatch/2, shapes[i][0]-(dpatch/2)-1),np.random.randint(0+dpatch/2, shapes[i][1]-(dpatch/2)-1),np.random.randint(0+dpatch/2, shapes[i][2]-(dpatch/2)-1)])
                #backgroundVoxels.extend([0])
            #Replace the ones that by chance are foreground voxels (not so many in tumor data)
            while any([e for e in foregroundVoxels if e in backgroundVoxels]):
                ix = [e for e in foregroundVoxels if e in backgroundVoxels]
                for index in ix:
                    newVoxel = [np.random.randint(dpatch/2, shapes[i][0]-(dpatch/2)-1),np.random.randint(dpatch/2, shapes[i][1]-(dpatch/2)-1),np.random.randint(dpatch/2, shapes[i][2]-(dpatch/2)-1)]
                    backgroundVoxels[backgroundVoxels.index(index)] = newVoxel
            allVoxelIndexes.append(foregroundVoxels + backgroundVoxels)
            random.shuffle(allVoxelIndexes[i])
        return allVoxelIndexes
    
def getAllForegroundClassesVoxels(groundTruthChannel, dpatch, output_classes):
    '''Get vector of voxel coordinates for all voxel values for all freground classes'''
    "e.g. groundTruthChannel = '/home/hirsch/Documents/projects/ATLASdataset/native_part2/c0011/c0011s0006t01/c0011s0006t01_LesionSmooth_Binary.nii.gz'"
    "NOTE: img in MRICRON starts at (1,1,1) and this function starts at (0,0,0), so points do not match when comparing in MRICRON. Add 1 to all dimensions to match in mricron. Function works properly though"
    img = nib.load(groundTruthChannel)
    data = np.array(img.dataobj[dpatch/2:img.shape[0]-(dpatch/2)-1, dpatch/2:img.shape[1]-(dpatch/2)-1, dpatch/2:img.shape[2]-(dpatch/2)-1],dtype='int16') # Get a cropped image, to avoid CENTRAL foreground voxels that are too near to the border. These will still be included, but not as central voxels. As long as they are in the 9x9x9 volume (-+ 4 voxels from the central, on a segment size of 25x25x25) they will still be included in the training.
    img.uncache()    
    voxels = []
    for c in range(1,output_classes):
        coords = np.argwhere(data==c)
        coords = coords + dpatch/2
        voxels.append(coords)
    return voxels  # This is a List! Use totuple() to convert if this makes any trouble
            
def getSubjectsToSample(channelList, subjectIndexes):
    "Actually returns channel of the subjects to sample"
    fp = open(channelList)
    lines = fp.readlines()
    subjects = [lines[i] for i in subjectIndexes]
    fp.close()
    return subjects

def extractLabels(groundTruthChannel_list, subjectIndexes, voxelCoordinates, output_dpatch):
    print('extracting labels from ' + str(len(subjectIndexes))+ ' subjects.')
    subjects = getSubjectsToSample(groundTruthChannel_list,subjectIndexes)
    labels = []
    for i in range(len(subjects)):
        subject = str(subjects[i])[:-1]
        proxy_label = nib.load(subject)
        label_data = np.array(proxy_label.get_fdata(),dtype='int8')
        if 0 not in label_data:
          label_data -= 1
        label_data = np.pad(label_data[:,:,:],((0, 2*output_dpatch), (0, 2*output_dpatch), (0, 2*output_dpatch)),'minimum')        
        for j in range(len(voxelCoordinates[i])):     
            D1,D2,D3 = voxelCoordinates[i][j]
            labels.append(label_data[D1-output_dpatch/2:D1+(output_dpatch/2)+1,D2-output_dpatch/2:D2+(output_dpatch/2)+1,D3-output_dpatch/2:D3+(output_dpatch/2)+1])
        proxy_label.uncache()
        del label_data
    return labels    


def extract_TPM_patches(TPM_channel, subjectIndexes, voxelCoordinates, output_dpatch):
    print('extracting TPM patches from ' + str(len(subjectIndexes))+ ' subjects.')
    n_patches = 0
    k = 0
    if not TPM_channel.endswith('.nii'):
        print('TPM channel is not a nifti file. Handling as a text file for Individual TPMs.')
        subjects_TPM = [x[:-1] for x in getSubjectsToSample(TPM_channel,subjectIndexes)]
    else:
        subjects_TPM = [TPM_channel]*len(subjectIndexes)

    for i in range(len(voxelCoordinates)):
        n_patches += len(voxelCoordinates[i])
    vol = np.zeros((n_patches,output_dpatch,output_dpatch,output_dpatch,6),dtype='float32')            
    for i in range(len(subjectIndexes)):  # same as iterating over subjects_TPM
        
        proxy_label = nib.load(subjects_TPM[i])
        label_data = np.array(proxy_label.get_fdata(),dtype='float32')

        # Take the log
        #label_data[label_data < 1e-7] = 1e-7
        #label_data = np.log(label_data)    

        pad = []
        for ch in range(label_data.shape[3]):
          pad.append(np.pad(label_data[:,:,:,ch],((0, 2*output_dpatch), (0, 2*output_dpatch), (0, 2*output_dpatch)),'minimum'))
        label_data = np.stack(pad, axis=3)
        del pad
        for j in range(0,len(voxelCoordinates[i])):     
            D1,D2,D3 = voxelCoordinates[i][j]
            
            vol[k,:,:,:,:] = label_data[D1-output_dpatch//2:D1+(output_dpatch//2)+1,
                                        D2-output_dpatch//2:D2+(output_dpatch//2)+1,
                                        D3-output_dpatch//2:D3+(output_dpatch//2)+1]
            

            k=k+1
        proxy_label.uncache()
        del label_data
    return vol    

def extractImagePatch(channel, subjectIndexes, patches, voxelCoordinates, dpatch, debug=False):
    subjects = getSubjectsToSample(channel, subjectIndexes)
    n_patches = 0   
    # Replace this thing. No need to compute. Have this information in list patches_per_subject!
    for i in range(len(voxelCoordinates)):
        n_patches += len(voxelCoordinates[i])
    #print('Starting extraction of {} patches from {} subjects.'.format(n_patches,len(voxelCoordinates)))
    vol = np.ones((n_patches,dpatch,dpatch,dpatch),dtype='float32')
    k = 0   
    for i in range(len(subjectIndexes)):
        #if i%20==0:
        #  print('{}%'.format(round(i*100./len(voxelCoordinates),2)))
        subject = str(subjects[i])[:-1]
        #print('Subject with path: {}'.format(subject))
        proxy_img = nib.load(subject)            
        img_data = np.array(proxy_img.get_fdata(),dtype='float32')
        #if percentile_normalization:
        img_data = percentile95_normalizeMRI(img_data)
        #else:
        #img_data = normalizeMRI(img_data)      
        if not np.isfinite(img_data).all():
        #if np.any(np.isnan(img_data)):
          print('Normalization: Nans found in scan {}'.format(subject))
          print('Nans replace by value: {}'.format(0.0))
          img_data[np.isfinite(img_data)] = 0.0
              
        npad = dpatch + 25
        # Need to pad. When training I extract border patches that go outside the image in any direction. When segmenting with large patches we get large patches on end-borders.
        # This thus depends on the size of images patches we are extracting, so we pad using that number.
        img_data_padded = np.pad(img_data,npad,'reflect')#'minimum')
        # Loop over voxelCoordinates tuples of subject i
        for j in range(len(voxelCoordinates[i])):      
            D1,D2,D3 = voxelCoordinates[i][j]           
            D1 += npad
            D2 += npad
            D3 += npad
            vol[k,:,:,:] = img_data_padded[D1-(dpatch//2):D1+(dpatch//2)+dpatch%2,
                                           D2-(dpatch//2):D2+(dpatch//2)+dpatch%2,
                                           D3-(dpatch//2):D3+(dpatch//2)+dpatch%2]
            k = k+1          
        proxy_img.uncache()
        del img_data
        if debug: print('extracted [' + str(len(voxelCoordinates[i])) + '] patches from subject ' + str(i) +'/'+ str(len(subjectIndexes)) +  ' with index [' + str(subjectIndexes[i]) + ']')        
    #print('In this batch found {} Bad Coordinates \n'.format(badCoords))
    #print('From subject(s): {}'.format(list(set(badCoords_subj))))
    #raw_input("Press Enter to continue...")
    return vol
    
    

#TPM_channel = cfg.TPM_channel
#trainChannels = cfg.trainChannels
#trainLabels = cfg.trainLabels
#validationChannels = cfg.validationChannels
#validationLabels = cfg.validationLabels 
#n_patches = cfg.n_patches
#n_subjects = cfg.n_subjects
#dpatch = cfg.dpatch
#output_classes = cfg.output_classes
#samplingMethod = cfg.samplingMethod_train
#model_name = cfg.model
#segmentation_dpatch = cfg.segmentation_dpatch
#logfile = None
#testChannels = cfg.validationChannels
#testLabels = cfg.validationLabels

def sampleTrainData(TPM_channel, trainChannels, trainLabels, n_patches, n_subjects, dpatch, output_classes, samplingMethod, model_name, data_augmentation=True, logfile = None, proportion_to_flip=0.5):
    '''output is a batch containing n-patches and their labels'''
    '''main function, called in the training process'''  
  
    if ('Baseline' in model_name) or (' ' in model_name) or ('DeepMedic' in model_name) or ('Vanilla' in model_name):
      output_dpatch = dpatch - 48
    else:
      output_dpatch = dpatch - 52    
    
    num_channels = len(trainChannels)
    start = time.time()
    patches_per_subject = get_patches_per_subject( n_patches, n_subjects)   
    labelsFile = open(trainLabels,"r")    
    total_subjects = file_len(trainLabels) # this was trainLabels before... was counting length of the string...???
    labelsFile.close()        
    subjectIndexes = generateRandomIndexesSubjects(n_subjects, total_subjects) 
  
    #------------ Parallelization Get Subject Shapes -------------------
    channel_mri = getSubjectChannels(subjectIndexes, trainChannels[0]) 
    pool = Pool(multiprocessing.cpu_count() - 1)#mp.cpu_count() -1)
    time1 = time.time()
    shapes = pool.map(getSubjectShapes_parallelization, channel_mri)
    print('Getting scan shapes took {} s'.format(round(time.time() - time1,2)))    
    pool.close()
    pool.join()     
    #--------------------------------------------------------------------
    
    voxelCoordinates = generateVoxelIndexes(subjectIndexes, shapes, patches_per_subject, dpatch, n_patches, trainLabels, samplingMethod, output_classes)
    # Get real number of patches to sample (as counted by the voxelCoordinates extracted, which is <= n_patches, as some classes are sparse)
    real_n_patches = 0
    for i in range(len(voxelCoordinates)):
        real_n_patches += len(voxelCoordinates[i])
    patches = np.zeros((real_n_patches,dpatch,dpatch,dpatch,num_channels),dtype='float32')
    for i in range(0,len(trainChannels)):
        patches[:,:,:,:,i] = extractImagePatch(trainChannels[i], subjectIndexes, patches, voxelCoordinates, dpatch, debug=False)             
    labels = np.array(extractLabels(trainLabels, subjectIndexes, voxelCoordinates, output_dpatch),dtype='int8')

    if len(TPM_channel) > 0:
      TPM_patches = np.array(extract_TPM_patches(TPM_channel, subjectIndexes, voxelCoordinates, output_dpatch),dtype='float32')    
    else:
      TPM_patches = []
    if(samplingMethod == 2):
        patches = patches[0:len(labels)]  # when using equal sampling (samplingMethod 2), because some classes have very few voxels in a head, there are fewer patches as intended. Patches is initialized as the maximamum value, so needs to get cut to match labels.
    end = time.time()
        
    if data_augmentation:
       print('Data augmentation: Randomly flipping {}% patches..'.format(proportion_to_flip*100))
       patches, labels, TPM_patches = flip_random(patches, labels, TPM_patches, proportion_to_flip)     

    labels = np.array(to_categorical(labels.astype(int),output_classes),dtype='int8')      
    
    del voxelCoordinates, subjectIndexes, patches_per_subject    
    
    if logfile is not None:
      my_logger("Finished extracting " + str(real_n_patches) + " patches, from "  + str(n_subjects) + " subjects and " + str(num_channels) + " channels. Timing: " + str(round(end-start,2)) + "s", logfile)
    return patches, labels, TPM_patches
       

def sampleTestData(TPM_channel, testChannels, testLabels, subjectIndex, output_classes, output_dpatch,logfile):
    labelsFile = open(testChannels[0],"r")   
    ch = labelsFile.readlines()
    subjectGTchannel = ch[subjectIndex[0]][:-1]
    my_logger('Segmenting subject with channel: ' + str(subjectGTchannel), logfile)
    labelsFile.close()      
    proxy_img = nib.load(subjectGTchannel)
    shape = proxy_img.shape
    affine = proxy_img.affine
    xend = output_dpatch * int(round(float(shape[0])/output_dpatch + 0.5)) # shape[0] + output_dpatch#-output_dpatch/2#5  #--> Size of final segmentation volume + 1
    yend = output_dpatch * int(round(float(shape[1])/output_dpatch + 0.5)) #shape[1] + output_dpatch#-output_dpatch/2+1+25#5
    zend = output_dpatch * int(round(float(shape[2])/output_dpatch + 0.5))#shape[2] + output_dpatch#-output_dpatch/2+1+25#5
    voxelCoordinates = []
    for x in range(output_dpatch//2,xend,output_dpatch):
        for y in range(output_dpatch//2,yend,output_dpatch):
            for z in range(output_dpatch//2,zend,output_dpatch):
                voxelCoordinates.append([x,y,z])
    if len(TPM_channel) > 0:
      TPM_patches = extract_TPM_patches(TPM_channel, subjectIndex, [voxelCoordinates], output_dpatch)
    else:
      TPM_patches = []
    if len(testLabels) > 0:
      labels = np.array(extractLabels(testLabels, subjectIndex, [voxelCoordinates], output_dpatch))
      labels = to_categorical(labels.astype(int),output_classes)
    else:
      labels = []
    #print("Finished extracting " + str(n_patches) + " patches, from "  + str(n_subjects) + " subjects and " + str(num_channels) + " channels. Timing: " + str(round(end-start,2)) + "s")
    return TPM_patches, labels, voxelCoordinates, shape, affine    
    
    
def get_patches_per_subject( n_patches, n_subjects):
    patches_per_subject = [n_patches/n_subjects]*n_subjects
    randomAdd = random.sample(range(len(patches_per_subject)),k=n_patches%n_subjects)
    randomAdd.sort()
    for index in randomAdd:
        patches_per_subject[index] = patches_per_subject[index] + 1
    return patches_per_subject
    
def getForegroundBackgroundVoxels(groundTruthChannel, dpatch):
    '''Get vector of voxel coordinates for all voxel values > 0'''
    "e.g. groundTruthChannel = '/home/hirsch/Documents/projects/ATLASdataset/native_part2/c0011/c0011s0006t01/c0011s0006t01_LesionSmooth_Binary.nii.gz'"
    "NOTE: img in MRICRON starts at (1,1,1) and this function starts at (0,0,0), so points do not match when comparing in MRICRON. Add 1 to all dimensions to match in mricron. Function works properly though"
    img = nib.load(groundTruthChannel)
    data = np.array(img.dataobj[dpatch/2:img.shape[0]-(dpatch/2)-1, dpatch/2:img.shape[1]-(dpatch/2)-1, dpatch/2:img.shape[2]-(dpatch/2)-1],dtype='int16') # Get a cropped image, to avoid CENTRAL foreground voxels that are too near to the border. These will still be included, but not as central voxels. As long as they are in the 9x9x9 volume (-+ 4 voxels from the central, on a segment size of 25x25x25) they will still be included in the training.
    img.uncache()    
    foregroundVoxels = np.argwhere(data>0)
    foregroundVoxels = foregroundVoxels + dpatch/2 # need to add this, as the cropped image starts again at (0,0,0)
    #backgroundVoxels = np.argwhere(data==0)
    return foregroundVoxels#, backgroundVoxels  # This is a List! Use totuple() to convert if this makes any trouble

#################################### DOCUMENTATION FUNCTIONS ####################################################

def my_logger(string, logfile):
    f = open(logfile,'a')
    f.write('\n' + str(string))
    f.close()
    print(string)    

def start_training_session_logger(logfile,threshold_EARLY_STOP, TPM_channel, load_model,saveSegmentation,path_to_model,model,dropout, trainChannels, trainLabels, validationChannels, validationLabels, num_iter, epochs, n_patches, n_patches_val, n_subjects, samplingMethod_train, size_minibatches, n_fullSegmentations, epochs_for_fullSegmentation, size_test_minibatches):
    my_logger('#######################################  NEW TRAINING SESSION  #######################################', logfile)    
    my_logger(trainChannels, logfile)
    my_logger(trainLabels, logfile)
    my_logger(validationChannels, logfile)        
    my_logger(validationLabels, logfile)  
    my_logger('TPM channel (if given):', logfile)
    my_logger(TPM_channel, logfile)
    my_logger('Session parameters: ', logfile)
    my_logger('[num_iter, epochs, n_patches, n_patches_val, n_subjects, samplingMethod_train, size_minibatches, n_fullSegmentations, epochs_for_fullSegmentation, size_test_minibatches]', logfile)
    my_logger([num_iter, epochs, n_patches, n_patches_val, n_subjects, samplingMethod_train, size_minibatches, n_fullSegmentations, epochs_for_fullSegmentation, size_test_minibatches], logfile)
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

class LossHistory_multiDice7(keras.callbacks.Callback):
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
        self.dice.append(logs.get('dice_coef_multilabel6'))
        self.metrics.append(self.dice)


class LossHistory_multiDice2(keras.callbacks.Callback):
    def on_train_begin(self, logs={}):
        self.losses = []
        self.dice = []
        self.metrics = []

    def on_batch_end(self, batch, logs={}):
        self.dice = []
        self.losses.append(logs.get('loss'))
        self.dice.append(logs.get('dice_coef_multilabel0'))
        self.dice.append(logs.get('dice_coef_multilabel1'))
        self.metrics.append(self.dice)   
   
########################################### MODEL TRAINING FUNCTIONS ##########################################################33   
   
def validation_on_batch(model, valbatch, TPM_patches, size_minibatches_val, vallabels, output_classes, logfile, model_name):
    batch_performance = []
    start = 0
    n_minibatches = len(valbatch)/size_minibatches_val
    for j in range(n_minibatches):
        print("validation on minibatch " +str(j+1)+ "/" + str(n_minibatches))
        end = start + size_minibatches_val
        minivalbatch = valbatch[start:end,:,:,:,:]    
        minivalbatch_labels = vallabels[start:end,:,:,:,:]    
        if len(TPM_patches) > 0:
          TPM_patch = TPM_patches[start:end,:,:,:,:]  
          if 'noDownsampling' in model_name:
            output_shape=[minivalbatch.shape[0], minivalbatch.shape[1]/3, minivalbatch.shape[2]/3, minivalbatch.shape[3]/3, minivalbatch.shape[4]]
            lowRes = np.zeros(output_shape)
            for iii in range(minivalbatch.shape[0]):
              lowRes[iii] = resize(minivalbatch[iii], output_shape[1:], anti_aliasing=True, preserve_range=True)
            highRes = minivalbatch[:,16:-16,16:-16,16:-16,:]
            #print('highRes shape, {}, lowRes.shape {}, TPM.shape {}'.format(highRes.shape, lowRes.shape, TPM_patch.shape))  
            batch_performance.append(model.evaluate([highRes, lowRes, TPM_patch], minivalbatch_labels, verbose=0))
          else:
            highRes = minivalbatch[:,16:-16,16:-16,16:-16,:]
            batch_performance.append(model.evaluate([highRes, TPM_patch], minivalbatch_labels, verbose=0))
        else:
          if 'noDownsampling' in model_name:
            output_shape=[minivalbatch.shape[0], minivalbatch.shape[1]/3, minivalbatch.shape[2]/3, minivalbatch.shape[3]/3, minivalbatch.shape[4]]
            lowRes = np.zeros(output_shape)
            for iii in range(minivalbatch.shape[0]):
              lowRes[iii] = resize(minivalbatch[iii], output_shape[1:], anti_aliasing=True, preserve_range=True)
            highRes = minivalbatch[:,16:-16,16:-16,16:-16,:]
            batch_performance.append(model.evaluate([highRes, lowRes], minivalbatch_labels, verbose=0))
          else:
            highRes = minivalbatch[:,16:-16,16:-16,16:-16,:]
            batch_performance.append(model.evaluate([highRes], minivalbatch_labels, verbose=0))
    val_performance = np.mean(batch_performance, 0)
    my_logger('Validation cost and accuracy ' + str(val_performance),logfile)  
    del valbatch
    del vallabels     
    return list(val_performance)

def train_model_on_batch(model,batch,TPM_patches,labels,size_minibatches,history,losses,metrics,logfile, output_classes, model_name):
    start = 0
    n_minibatches = len(batch)/size_minibatches
    for j in range(n_minibatches):
        print("training on minibatch " +str(j+1)+ "/" + str(n_minibatches))
        end = start + size_minibatches
        minibatch = batch[start:end,:,:,:,:]             
        minibatch_labels = labels[start:end,:,:,:,:]   
        if len(TPM_patches) > 0:
          TPM_patch = TPM_patches[start:end,:,:,:,:]    
          if 'noDownsampling' in model_name:
            output_shape=[minibatch.shape[0], minibatch.shape[1]/3, minibatch.shape[2]/3, minibatch.shape[3]/3, minibatch.shape[4]]
            lowRes = np.zeros(output_shape)
            for iii in range(minibatch.shape[0]):
              lowRes[iii] = resize(minibatch[iii], output_shape[1:], anti_aliasing=True, preserve_range=True)
            highRes = minibatch[:,16:-16,16:-16,16:-16,:]
            model.fit([highRes,lowRes, TPM_patch], minibatch_labels, verbose = 0, callbacks = [history])
          else:  
            highRes = minibatch[:,16:-16,16:-16,16:-16,:]
            model.fit([highRes, TPM_patch], minibatch_labels, verbose = 0, callbacks = [history])
        else:
          if 'noDownsampling' in model_name:
            output_shape=[minibatch.shape[0], minibatch.shape[1]/3, minibatch.shape[2]/3, minibatch.shape[3]/3, minibatch.shape[4]]
            lowRes = np.zeros(output_shape)
            for iii in range(minibatch.shape[0]):
              lowRes[iii] = resize(minibatch[iii], output_shape[1:], anti_aliasing=True, preserve_range=True)
            highRes = minibatch[:,16:-16,16:-16,16:-16,:]
            model.fit([highRes,lowRes], minibatch_labels, verbose = 0, callbacks = [history])
          else:  
            highRes = minibatch[:,16:-16,16:-16,16:-16,:]
            model.fit([highRes], minibatch_labels, verbose = 0, callbacks = [history])
        losses.extend(history.losses)
        metrics.extend(history.metrics)
        start = end
        my_logger('Train cost and metrics     ' + str(losses[-1]) + ' ' + str(metrics[-1]),logfile)
    del batch
    del labels
    
    
    
def fullHeadSegmentation(wd, penalty_MATRIX, TPM_channel, dice_compare, dsc, model, model_name, testChannels, testLabels, subjectIndex, output_classes, segmentation_dpatch, size_minibatches,logfile, epoch, saveSegmentation = False):    

    if ('Baseline' in model_name) or ('MultiPriors' in model_name) or ('DeepMedic' in model_name) or ('Vanilla' in model_name):
      output_dpatch = segmentation_dpatch - 48
    else:
      output_dpatch = segmentation_dpatch - 52    

    if len(testLabels) == 0:
        dice_compare = False    
    subjectIndex = [subjectIndex]
    flairCh = getSubjectsToSample(testChannels[0], subjectIndex)
    subID = flairCh[0].split('.')[0].split('/')[-1][:7]
    segmentationName = 'predictions/' + subID + '_segmentation_epoch' + str(epoch)
    output = wd +'/' + segmentationName + '.nii.gz'
    if os.path.exists(output):
      print('Segmentation already found in: {}\nSkip..'.format(output))
      return

    # Open subject MRI to extract patches dynamically
    num_channels = len(testChannels)
    firstChannelFile = open(testChannels[0],"r")   
    ch = firstChannelFile.readlines()
    subjectGTchannel = ch[subjectIndex[0]][:-1]
    my_logger('Segmenting subject with channel: ' + str(subjectGTchannel), logfile)
    firstChannelFile.close()      
    proxy_img = nib.load(subjectGTchannel)
    shape = proxy_img.shape
    affine = proxy_img.affine
    
    TPM_patches, labels, voxelCoordinates, shape, affine = sampleTestData(TPM_channel, testChannels, testLabels, subjectIndex, output_classes, output_dpatch,logfile)
    
    print("Extracted image patches for full head segmentation.")
    print('Len TPM patches: {}'.format(len(TPM_patches)))
    print('Shape: {}'.format(shape))
    start = 0
    n_minibatches = len(voxelCoordinates)//size_minibatches
    indexes = []
    #print('Number of patches: {}'.format(n_minibatches))
    for j in range(0,n_minibatches):
        print("Segmenting patch " +str(j)+ "/" + str(n_minibatches))
        end = start + size_minibatches
        patches = np.zeros((size_minibatches,segmentation_dpatch,segmentation_dpatch,segmentation_dpatch,num_channels),dtype='float32')
        minibatch_voxelCoordinates = voxelCoordinates[start:end]
        for i in range(0,len(testChannels)):
            patches[:,:,:,:,i] = extractImagePatch(testChannels[i], subjectIndex, patches, [minibatch_voxelCoordinates], segmentation_dpatch, debug=False)  
        if len(TPM_patches) > 0:                                
          TPM_batch = TPM_patches[start:end,:,:,:,:]  
          if 'noDownsampling' in model_name:
            output_shape=[patches.shape[0], patches.shape[1]//3, patches.shape[2]//3, patches.shape[3]//3, patches.shape[4]]
            lowRes = np.zeros(output_shape)
            for iii in range(patches.shape[0]):
              lowRes[iii] = resize(patches[iii], output_shape[1:], anti_aliasing=True, preserve_range=True)
            highRes = patches[:,16:-16,16:-16,16:-16,:]
            #print('highRes shape, {}, lowRes.shape {}, TPM.shape {}'.format(highRes.shape, lowRes.shape, TPM_batch.shape))  
            prediction = model.predict([highRes, lowRes, TPM_batch], verbose=0)
            
          else:
            highRes = patches[:,16:-16,16:-16,16:-16,:]
            prediction = model.predict([highRes,TPM_batch], verbose=0)
        else:
          if 'noDownsampling' in model_name:
            output_shape=[patches.shape[0], patches.shape[1]/3, patches.shape[2]/3, patches.shape[3]/3, patches.shape[4]]
            lowRes = np.zeros(output_shape)
            for iii in range(patches.shape[0]):
              lowRes[iii] = resize(patches[iii], output_shape[1:], anti_aliasing=True, preserve_range=True)
            highRes = patches[:,16:-16,16:-16,16:-16,:]
            prediction = model.predict([highRes, lowRes], verbose=0)
          else:  
            highRes = patches[:,16:-16,16:-16,16:-16,:]
            prediction = model.predict([highRes], verbose=0)

        class_pred = np.argmax(prediction, axis=4)
        indexes.extend(class_pred)        
        start = end
    #last one
    
    size_last_minibatch = (len(voxelCoordinates)-n_minibatches*size_minibatches)
    if size_last_minibatch > 0:
      print("Segmenting LAST patch " +str(j)+ "/" + str(n_minibatches))
      end = start + size_minibatches
      patches = np.zeros((size_minibatches,segmentation_dpatch,segmentation_dpatch,segmentation_dpatch,num_channels),dtype='float32')
      minibatch_voxelCoordinates = voxelCoordinates[start:end]
      for i in range(0,len(testChannels)):
          patches[:,:,:,:,i] = extractImagePatch(testChannels[i], subjectIndex, patches, [minibatch_voxelCoordinates], segmentation_dpatch, debug=False)  
      if len(TPM_patches) > 0:                                
        TPM_batch = TPM_patches[start:end,:,:,:,:]  
        if 'noDownsampling' in model_name:
          output_shape=[patches.shape[0], patches.shape[1]/3, patches.shape[2]/3, patches.shape[3]/3, patches.shape[4]]
          lowRes = np.zeros(output_shape)
          for iii in range(patches.shape[0]):
            lowRes[iii] = resize(patches[iii], output_shape[1:], anti_aliasing=True, preserve_range=True)
          highRes = patches[:,16:-16,16:-16,16:-16,:]
          #print('highRes shape, {}, lowRes.shape {}, TPM.shape {}'.format(highRes.shape, lowRes.shape, TPM_batch.shape))  
          prediction = model.predict([highRes, lowRes, TPM_batch], verbose=0)
          
        else:
          highRes = patches[:,16:-16,16:-16,16:-16,:]
          prediction = model.predict([highRes,TPM_batch], verbose=0)
      else:
        if 'noDownsampling' in model_name:
          output_shape=[patches.shape[0], patches.shape[1]/3, patches.shape[2]/3, patches.shape[3]/3, patches.shape[4]]
          lowRes = np.zeros(output_shape)
          for iii in range(patches.shape[0]):
            lowRes[iii] = resize(patches[iii], output_shape[1:], anti_aliasing=True, preserve_range=True)
          highRes = patches[:,16:-16,16:-16,16:-16,:]
          prediction = model.predict([highRes, lowRes], verbose=0)
        else:  
          highRes = patches[:,16:-16,16:-16,16:-16,:]
          prediction = model.predict([highRes], verbose=0)

      class_pred = np.argmax(prediction, axis=4)
      indexes.extend(class_pred)    
      
      
      
#      print("Segmenting patch " +str(j)+ "/" + str(n_minibatches))
#      end = start + size_last_minibatch
#      patches = np.zeros((size_last_minibatch,segmentation_dpatch,segmentation_dpatch,segmentation_dpatch,num_channels),dtype='float32')
#      minibatch_voxelCoordinates = voxelCoordinates[start:end]
#      for i in range(0,len(testChannels)):
#          patches[:,:,:,:,i] = extractImagePatch(testChannels[i], subjectIndex, patches, [minibatch_voxelCoordinates], segmentation_dpatch, debug=False)            
#
#      if len(TPM_patches) > 0:          
#        TPM_batch = TPM_patches[start:end,:,:,:,:]  
#        prediction = model.predict([patches,TPM_batch], verbose=0)
#      else:
#
#        prediction = model.predict([patches], verbose=0)
#
#      class_pred = np.argmax(prediction, axis=4)
#      indexes.extend(class_pred)     
      del patches  
     

    head = np.zeros(shape, dtype=np.int16)  
    i = 0
    for x,y,z in voxelCoordinates:
        
        patch_shape = head[x-output_dpatch//2:min(x+output_dpatch//2+output_dpatch%2, shape[0]),
                           y-output_dpatch//2:min(y+output_dpatch//2+output_dpatch%2, shape[1]),
                           z-output_dpatch//2:min(z+output_dpatch//2+output_dpatch%2, shape[2])].shape
        
        head[x-output_dpatch//2:min(x+output_dpatch//2+output_dpatch%2, shape[0]),
             y-output_dpatch//2:min(y+output_dpatch//2+output_dpatch%2, shape[1]),
             z-output_dpatch//2:min(z+output_dpatch//2+output_dpatch%2, shape[2])] = np.array(indexes[i])[0:patch_shape[0], 0:patch_shape[1],0:patch_shape[2]]
        i = i+1
    img = nib.Nifti1Image(head, affine)
        
    if dice_compare:           
        labelsFile = open(testLabels,"r")   
        ch = labelsFile.readlines()
        subjectGTchannel = ch[subjectIndex[0]][:-1]
        GT = nib.load(subjectGTchannel).get_fdata()  
        if 0 not in GT:
          GT -= 1
        score = weighted_generalized_dice_completeImages(GT, img.get_fdata(), penalty_MATRIX)
        #score = generalized_dice_completeImages(img.get_fdata(), GT.get_fdata())
        dsc.append(score[0])
        print(dsc[-1])
        print('per class dice score: {}'.format(score[1]))
        print('mean DCS so far:' + str(np.mean(dsc)))
        
    if(saveSegmentation):

        nib.save(img, output)
        print('Saved segmentation of subject at: ' + output)
        response = input('Do you want to add path to save results (yes or no)? ')
        if response.lower() == 'yes':
            output2 = input('Path: ') + '\\' + subID + '_segmentation_epoch' + str(epoch) + '.nii.gz'
            nib.save(img, output2)
            print('Saved segmentation of subject at: ' + output2)
        my_logger('Saved segmentation of subject at: ' + output, logfile)
       
      



def fullHeadSegmentation_Flexible(wd, penalty_MATRIX, TPM_channel, dice_compare, dsc, model, testChannels, testLabels, subjectIndex, output_classes, dpatch, size_minibatches,logfile, epoch, saveSegmentation = False, full_evaluation = False):    
    
    subjectIndex = [subjectIndex]
    subject_channels = []
    for modality in testChannels:
      channelsFile = open(testChannels[0],"r")   
      ch = channelsFile.readlines()
      subject_channels.append(ch[subjectIndex[0]][:-1])#
      channelsFile.close()
    my_logger('Segmenting subject with channels: ' + str(subject_channels), logfile)  
    
    subID = subject_channels[0].split('/')[-1][:7]
    
    images = []
    for channel in subject_channels: 
      proxy_img = nib.load(channel)
      X = proxy_img.get_fdata()
      images.append(X)

    # load TPM and add to stack of inputs.
    tpm_nii = nib.load(TPM_channel)      
    TPM_data = tpm_nii.get_fdata()
    
    X = np.stack(images, axis=3)
    print(X.shape)
    shape = X.shape
    original_shape = shape
      
    
    if shape[0]*shape[1]*shape[2] > 182**3:  
      boundary = 275
      # Set boundaries for maximum allowed shape
      a = np.max([0,(shape[0] - boundary)])/2   
      b = np.max([0,(shape[1] - boundary)])/2
      c = np.max([0,(shape[2] - boundary)])/2    
      X = X[a:shape[0]-a,:,:,:]
      X = X[:,b:shape[1]-b,:,:]
      X = X[:,:,(c+1):shape[2]-c,:]

      TPM_data = TPM_data[26:X.shape[0]-26,:,:,:]
      TPM_data = TPM_data[:,26:X.shape[1]-26,:,:]
      TPM_data = TPM_data[:,:,26:X.shape[2]-26,:]     
      
    shape = X.shape

    X = X.reshape(((1,)+ X.shape))
    TPM_data = TPM_data.reshape((1,)+TPM_data.shape)
    
    
    # Padding
    X_pad = []
    for m in range(X.shape[4]):      
      #print('padding modality {}'.format(m))
      X_pad.append(np.pad(X[0,:,:,:,m],10,'constant'))
    TPM_pad = []
    for m in range(TPM_data.shape[4]):  
      TPM_pad.append(np.pad(TPM_data[0,:,:,:,m], 10, 'constant'))
  
    X_padded = np.stack(X_pad,axis=3)
    X_padded = X_padded.reshape(((1,) + X_padded.shape))
    
    TPM_data_pad = np.stack(TPM_pad,axis=3)    
    TPM_data_pad = TPM_data_pad.reshape(((1,) + TPM_data_pad.shape))
    print('X shape: {}'.format(X_padded.shape))
    print('TPM shape: {}'.format(TPM_data_pad.shape))
    time_segment = time.time()
    yhat = model.predict([X_padded, TPM_data_pad])
    print('Segmentation done. Time: {}'.format(time.time() - time_segment))
    y = np.argmax(yhat, axis=4)   # For classification output
    #y = yhat[:,:,:,:,1]            # For logits for class 2
    print('VALUES IN Y: {}'.format(np.unique(y)))
    
    y = y.reshape(y.shape[1],y.shape[2],y.shape[3])
    #y = y.reshape(shape[0]-24,shape[1]-78,shape[2]-78)
    print('segmentation output shape: {}'.format(y.shape))

    y_out = np.zeros((original_shape[0],original_shape[1],original_shape[2]))


    y_out[abs(original_shape[0] -y.shape[0])/2:original_shape[0] - abs(original_shape[0] -y.shape[0])/2,
          abs(original_shape[1] -y.shape[1])/2:original_shape[1] - abs(original_shape[1] -y.shape[1])/2,
          abs(original_shape[2] -y.shape[2])/2:original_shape[2] - abs(original_shape[2] -y.shape[2])/2 -1] = y
          
    print('VALUES IN Y_OUT: {}'.format(np.unique(y_out)))          
          
    labelsFile = open(testLabels,"r")   
    ch = labelsFile.readlines()
    subjectGTchannel = ch[subjectIndex[0]][:-1]
    GT = nib.load(subjectGTchannel)      
    print('Subject GT channel: {}'.format(subjectGTchannel))       
    y_out = np.array(y_out, dtype='int16')
    img = nib.Nifti1Image(y_out, GT.affine)


    
    for att in img.header:
      img.header[att] = GT.header[att]
    print('VALUES IN IM : {}'.format(np.unique(img.get_fdata())))      

    
    segmentationName = '/predictions/' + subID + 'segmented_epoch_' + str(epoch)
    output = wd +'/' + segmentationName + '.nii'
    nib.save(img, output)
    #my_logger('Saved segmentation of subject at: ' + output, logfile)
    if dice_compare:           
        labelsFile = open(testLabels,"r")   
        ch = labelsFile.readlines()
        subjectGTchannel = ch[subjectIndex[0]][:-1]
        GT = nib.load(subjectGTchannel)  
        score = weighted_generalized_dice_completeImages(GT.get_fdata(), img.get_fdata(), penalty_MATRIX)
        #score = generalized_dice_completeImages(img.get_fdata(), GT.get_fdata())
        dsc.append(score[0])
        print(dsc[-1])
        print('per class dice score: {}'.format(score[1]))
        print('mean DCS so far:' + str(np.mean(dsc)))


############################## create or load model ###########################################

def make_model(model_configFile):
    
    path = '/'.join(model_configFile.split('/')[:-1])
    model_configFileName = model_configFile.split('/')[-1][:-3]   
    sys.path.append(path)
    cfg = __import__(model_configFileName)  
    if cfg.model == 'CNN_TPM':
        from multiscale_CNN_TPM import multiscale_CNN_TPM            
        dm = multiscale_CNN_TPM(cfg.dpatch, cfg.output_classes, cfg.num_channels, cfg.L2, cfg.dropout, cfg.learning_rate, cfg.optimizer_decay, cfg.loss_function)
        model = dm.createModel()          
    elif cfg.model == 'DeepMedic':
        from DeepMedic_model import DeepMedic
        dm = DeepMedic(cfg.dpatch, cfg.output_classes, cfg.num_channels, cfg.L2, cfg.dropout, cfg.learning_rate, cfg.optimizer_decay, cfg.loss_function)
        model = dm.createModel()
    else: 
        print('ERROR: No model selected.')
        return 0
    print(model.summary())
    from keras.utils import plot_model
    plot_model(model, to_file=cfg.workingDir+'/models/'+cfg.model +'.png', show_shapes=True)
    model_path = cfg.workingDir+'/models/'+cfg.model +'.h5'
    model.save(model_path)
    print('Saved model at {}'.format(model_path))


############################## RUNNING TRAINING AND SEGMENTATION FUNCTIONS ###########################################

#configFile = '/home/hirsch/Documents/projects/brainSegmentation/DeepPriors/configFiles/configFiles_stroke/configFile_Aphasic_Stroke_MultiPriors_noDownsampling.py'
#workingDir = '/home/hirsch/Documents/projects/brainSegmentation/DeepPriors/'

def train_test_model(configFile, workingDir):
    # import configuration file and create working environment
    print(configFile)
    path = '/'.join(configFile.split('/')[:-1])
    print(path)
    configFileName = configFile.split('/')[-1][:-3]   
    sys.path.append(path)
    cfg = __import__(configFileName)
    
    if len(cfg.TPM_channel) > 0:
      cfg.TPM_channel = workingDir + cfg.TPM_channel
    if len(cfg.train_TPM_channel) > 0:
      cfg.train_TPM_channel = workingDir + cfg.train_TPM_channel
      cfg.TPM_channel = cfg.train_TPM_channel  
    if len(cfg.validation_TPM_channel) > 0:
      cfg.validation_TPM_channel = workingDir + cfg.validation_TPM_channel
      
    cfg.trainChannels = [workingDir + x for x in cfg.trainChannels]
    cfg.trainLabels = workingDir +cfg.trainLabels 

    cfg.validationChannels = [workingDir + x for x in cfg.validationChannels]
    if len(cfg.validationLabels) != 0:
      cfg.validationLabels = workingDir +cfg.validationLabels
    dice_compare = cfg.dice_compare  
    # Create or load CNN model
    if cfg.load_model == False:
        if cfg.model == 'CNN_TPM':
            from multiscale_CNN_TPM import multiscale_CNN_TPM            
            dm = multiscale_CNN_TPM(cfg.dpatch, cfg.output_classes, cfg.num_channels, cfg.L2, cfg.dropout, cfg.learning_rate, cfg.optimizer_decay, cfg.loss_function)
            model = dm.createModel()                
        elif cfg.model == 'DeepMedic':
            from DeepMedic_model import DeepMedic
            dm = DeepMedic(cfg.dpatch, cfg.output_classes, cfg.num_channels, cfg.L2, cfg.dropout, cfg.learning_rate, cfg.optimizer_decay, cfg.loss_function)
            model = dm.createModel()            
        elif cfg.model == 'BIG_multiscale_CNN_TPM':
            from BIG_multiscale_CNN_TPM import BIG_multiscale_CNN_TPM
            dm = BIG_multiscale_CNN_TPM(cfg.dpatch, cfg.output_classes, cfg.num_channels, cfg.L2, cfg.dropout, cfg.learning_rate, cfg.optimizer_decay, cfg.loss_function)
            model = dm.createModel()   
        elif cfg.model == 'BIG_multiscale_CNN_TPM_flexible':
            from BIG_multiscale_CNN_TPM_flexibleInput import BIG_multiscale_CNN_TPM_flexibleInput
            dm = BIG_multiscale_CNN_TPM_flexibleInput(cfg.dpatch, cfg.output_classes, cfg.num_channels, cfg.L2, cfg.dropout, cfg.learning_rate, cfg.optimizer_decay, cfg.loss_function)
            model = dm.createModel()   
        elif cfg.model == 'BIG_singleScale_CNN_TPM':
            from BIG_singleScale_CNN_TPM import BIG_singleScale_CNN_TPM
            dm = BIG_singleScale_CNN_TPM(cfg.dpatch, cfg.output_classes, cfg.num_channels, cfg.L2, cfg.dropout, cfg.learning_rate, cfg.optimizer_decay, cfg.loss_function)
            model = dm.createModel()               
        elif cfg.model == 'DeepPriors_CRFasRNN':
            from DeepPriors_CRFasRNN import DeepPriors_CRFasRNN
            dm = DeepPriors_CRFasRNN(cfg.dpatch, cfg.output_classes, cfg.num_channels, cfg.L2, cfg.dropout, cfg.learning_rate, cfg.optimizer_decay, cfg.loss_function)
            model = dm.createModel()
        elif cfg.model == 'DeepMedic_TPM':
            from DeepMedic_TPM import DeepMedic
            dm = DeepMedic(cfg.output_classes, cfg.num_channels, cfg.L2, cfg.dropout, cfg.learning_rate, cfg.optimizer_decay, cfg.loss_function)
            model = dm.createModel()            
            model.summary()               
        elif cfg.model == 'MultiPriors_noDownsampling':
            from MultiPriors_noDownsampling import MultiPriors_Model
            dm = MultiPriors_Model(cfg.output_classes, cfg.num_channels, cfg.L2, cfg.dropout, cfg.learning_rate, cfg.optimizer_decay, cfg.loss_function)
            model = dm.createModel()            
            model.summary()                  
        elif cfg.model == 'DeepMedic_noDownsampling':
            from DeepMedic_noDownsampling import MultiPriors_Model
            dm = MultiPriors_Model(cfg.output_classes, cfg.num_channels, cfg.L2, cfg.dropout, cfg.learning_rate, cfg.optimizer_decay, cfg.loss_function)
            model = dm.createModel()            
            model.summary()                        
        elif cfg.model == 'MultiPriors':
            from MultiPriors import MultiPriors_Model
            dm = MultiPriors_Model(cfg.output_classes, cfg.num_channels, cfg.L2, cfg.dropout, cfg.learning_rate, cfg.optimizer_decay, cfg.loss_function)
            model = dm.createModel()            
            model.summary()          
        elif cfg.model == 'MultiPriors_noTPM':
            from MultiPriors_noTPM import MultiPriors_noTPM_Model
            dm = MultiPriors_noTPM_Model(cfg.output_classes, cfg.num_channels, cfg.L2, cfg.dropout, cfg.learning_rate, cfg.optimizer_decay, cfg.loss_function)
            model = dm.createModel()            
            model.summary()                      
        elif cfg.model == 'Vanilla':
            from Vanilla import DeepMedic
            dm = DeepMedic(cfg.output_classes, cfg.num_channels, cfg.L2, cfg.dropout, cfg.learning_rate, cfg.optimizer_decay, cfg.loss_function)
            model = dm.createModel()            
            model.summary()      
        elif cfg.model == 'Baseline':
            from Baseline import Baseline_model
            dm = Baseline_model(cfg.output_classes, cfg.num_channels, cfg.L2, cfg.dropout, cfg.learning_rate, cfg.optimizer_decay, cfg.loss_function)
            model = dm.createModel()            
            model.summary()                
        elif cfg.model == 'Baseline_TPM':
            from Baseline_TPM import Baseline_model
            dm = Baseline_model(cfg.output_classes, cfg.num_channels, cfg.L2, cfg.dropout, cfg.learning_rate, cfg.optimizer_decay, cfg.loss_function)
            model = dm.createModel()            
            model.summary()                        
        elif cfg.model == 'Vanilla_TPM':
            from Vanilla_TPM import DeepMedic
            dm = DeepMedic(cfg.output_classes, cfg.num_channels, cfg.L2, cfg.dropout, cfg.learning_rate, cfg.optimizer_decay, cfg.loss_function)
            model = dm.createModel()            
            model.summary()                  
            
        else: 
            print('ERROR: No model selected.')
            return 0       
  
        start_epoch = 0
        os.chdir(workingDir + '/training_sessions/')
        session = cfg.model + '_' + cfg.dataset + '_' + configFileName + '_' + time.strftime("%Y-%m-%d_%H%M") 
        wd = workingDir + '/training_sessions/' +session
        if not os.path.exists(wd):    
            os.mkdir(session)
            os.mkdir(session + '/models')
            os.mkdir(session + '/predictions')
        os.chdir(wd)  
        copy(workingDir + configFile[1:], wd)
        logfile = session +'.log'
        print(model.summary())
        val_performance = []
        full_segm_DICE = []
        best_results_epoch = 0
        from keras.utils import plot_model
        plot_model(model, to_file=wd+'/Model_Architecture.png', show_shapes=True)
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
            from multiscale_CNN_TPM import dice_coef_multilabel6, dice_coef_multilabel0,dice_coef_multilabel1,dice_coef_multilabel2,dice_coef_multilabel3,dice_coef_multilabel4,dice_coef_multilabel5
            my_custom_objects = {'dice_coef_multilabel6':dice_coef_multilabel6,
                                 'dice_coef_multilabel0':dice_coef_multilabel0,
                                 'dice_coef_multilabel1':dice_coef_multilabel1,
                                 'dice_coef_multilabel2':dice_coef_multilabel2,
                                 'dice_coef_multilabel3':dice_coef_multilabel3,
                                 'dice_coef_multilabel4':dice_coef_multilabel4,
                                 'dice_coef_multilabel5':dice_coef_multilabel5}
            #custom_metrics = [dice_coef_multilabel6,dice_coef_multilabel0,dice_coef_multilabel1,dice_coef_multilabel2,dice_coef_multilabel3,dice_coef_multilabel4,dice_coef_multilabel5]
        elif cfg.loss_function == 'Dice7':
            from BIG_multiscale_CNN_TPM import Generalised_dice_coef_multilabel7, dice_coef_multilabel6, dice_coef_multilabel0,dice_coef_multilabel1,dice_coef_multilabel2,dice_coef_multilabel3,dice_coef_multilabel4,dice_coef_multilabel5
            my_custom_objects = {'Generalised_dice_coef_multilabel7':Generalised_dice_coef_multilabel7,
                                 'dice_coef_multilabel0':dice_coef_multilabel0,
                                 'dice_coef_multilabel1':dice_coef_multilabel1,
                                 'dice_coef_multilabel2':dice_coef_multilabel2,
                                 'dice_coef_multilabel3':dice_coef_multilabel3,
                                 'dice_coef_multilabel4':dice_coef_multilabel4,
                                 'dice_coef_multilabel5':dice_coef_multilabel5,
                                 'dice_coef_multilabel6':dice_coef_multilabel6}
        elif cfg.loss_function == 'wDice6':
            from multiscale_CNN_TPM import w_dice_coef_multilabel6, dice_coef_multilabel0,dice_coef_multilabel1,dice_coef_multilabel2,dice_coef_multilabel3,dice_coef_multilabel4,dice_coef_multilabel5
            my_custom_objects = {'w_dice_coef_multilabel6':w_dice_coef_multilabel6,
                                 'dice_coef_multilabel0':dice_coef_multilabel0,
                                 'dice_coef_multilabel1':dice_coef_multilabel1,
                                 'dice_coef_multilabel2':dice_coef_multilabel2,
                                 'dice_coef_multilabel3':dice_coef_multilabel3,
                                 'dice_coef_multilabel4':dice_coef_multilabel4,
                                 'dice_coef_multilabel5':dice_coef_multilabel5}
        elif cfg.loss_function == 'Dice2':
            from multiscale_CNN_TPM import Generalised_dice_coef_multilabel2, dice_coef_multilabel0,dice_coef_multilabel1
            my_custom_objects = {'Generalised_dice_coef_multilabel2':Generalised_dice_coef_multilabel2,
                                 'dice_coef_multilabel0':dice_coef_multilabel0,
                                 'dice_coef_multilabel1':dice_coef_multilabel1}
        elif cfg.loss_function == 'wDice2':
            from multiscale_CNN_TPM import w_dice_coef_multilabel2, dice_coef_multilabel0,dice_coef_multilabel1
            my_custom_objects = {'w_dice_coef_multilabel2':w_dice_coef_multilabel2,
                                 'dice_coef_multilabel0':dice_coef_multilabel0,
                                 'dice_coef_multilabel1':dice_coef_multilabel1} 
        elif cfg.loss_function == 'Multinomial':
            from multiscale_CNN_TPM import dice_coef_multilabel0,dice_coef_multilabel1,dice_coef_multilabel2,dice_coef_multilabel3,dice_coef_multilabel4,dice_coef_multilabel5,dice_coef_multilabel6
            my_custom_objects = {'dice_coef_multilabel0':dice_coef_multilabel0,
                                 'dice_coef_multilabel1':dice_coef_multilabel1,
                                 'dice_coef_multilabel2':dice_coef_multilabel2,
                                 'dice_coef_multilabel3':dice_coef_multilabel3,
                                 'dice_coef_multilabel4':dice_coef_multilabel4,
                                 'dice_coef_multilabel5':dice_coef_multilabel5,
                                 'dice_coef_multilabel6':dice_coef_multilabel6} 

        model = load_model(cfg.path_to_model, custom_objects = my_custom_objects )
        print('LOADED MODEL FROM SESSION {}'.format(cfg.session))
        session = cfg.session
        start_epoch = int(cfg.path_to_model.split('.')[-2][cfg.path_to_model.split('.')[-2].find('epoch') + 5 : ]) + 1
        cfg.epochs_for_fullSegmentation = range(start_epoch+1, cfg.epochs)
        os.chdir(workingDir + '/training_sessions/')
        wd = workingDir + '/training_sessions/' +session
        dice_file = wd + '/Dice_scores.txt'
        if os.path.exists(dice_file):
        	full_segm_DICE = [float(x[:-1]) for x in open(dice_file).readlines()]
	        best_results_epoch = full_segm_DICE.index(np.max(full_segm_DICE))
        else:
            full_segm_DICE = []
            results_epoch = 0
        if not os.path.exists(wd):    
            os.mkdir(session)
            os.mkdir(session + '/models')
            os.mkdir(session + '/predictions')
        os.chdir(wd)   
        logfile = session +'.log'
    #################################################################################################
    #                                                                                               #
    #                                         START SESSION                                         #
    #                                                                                               #
    #################################################################################################
    # OUTCOMMENTED SO I CAN KEEP USING SAME TRAINING DATA FOR SAME MODEL.
        
    val_performance = []
    full_segm_DICE = []
    losses = []
    metrics = []
    np.set_printoptions(precision=3)

    start_training_session_logger(logfile, cfg.threshold_EARLY_STOP, cfg.TPM_channel, cfg.load_model, cfg.saveSegmentation, cfg.path_to_model, model, \
        cfg.dropout, cfg.trainChannels, cfg.trainLabels, cfg.validationChannels, cfg.validationLabels, \
        cfg.num_iter, cfg.epochs, cfg.n_patches, cfg.n_patches_val, cfg.n_subjects, cfg.samplingMethod_train, \
        cfg.size_minibatches, cfg.n_fullSegmentations, cfg.epochs_for_fullSegmentation, cfg.size_test_minibatches)
    # Callback history    
    if cfg.output_classes == 2:
        history = LossHistory_multiDice2() 
    elif cfg.output_classes == 6:
        history = LossHistory_multiDice6()
    elif cfg.output_classes == 7:
        history = LossHistory_multiDice7()
    for epoch in range(start_epoch,cfg.epochs+1):
      t1 = time.time()
      my_logger("######################################################",logfile)
      my_logger("                   TRAINING EPOCH " + str(epoch) + "/" + str(cfg.epochs),logfile)
      my_logger("######################################################",logfile)
              
      ####################### FULL HEAD SEGMENTATION ##############################
      
      if epoch in cfg.epochs_for_fullSegmentation:
        my_logger("------------------------------------------------------", logfile)
        my_logger("                 FULL HEAD SEGMENTATION", logfile)
        my_logger("------------------------------------------------------", logfile)
        
        if len(cfg.validationLabels) > 0:
          dice_compare = True
        dsc = []
        subjectIndex = 0
        with open(cfg.validationChannels[0]) as vl:
          n_valSubjects = len(vl.readlines())
          print('Using {} val subjects'.format(n_valSubjects))
        if cfg.test_subjects > n_valSubjects:
          print("Given number of subjects for test set (" + str(cfg.test_subjects) +") is larger than the amount of \
          subjects in test set (" +str(n_valSubjects)+ ")")
          cfg.test_subjects = n_valSubjects
          cfg.n_fullSegmentations = min(n_valSubjects, cfg.n_fullSegmentations)
        list_subjects_fullSegmentation = sample(range(cfg.test_subjects), cfg.n_fullSegmentations)
        for subjectIndex in list_subjects_fullSegmentation: 
          time_segmentation = time.time()
          fullHeadSegmentation(wd, cfg.penalty_MATRIX, cfg.validation_TPM_channel, dice_compare, dsc, model, cfg.model, cfg.validationChannels, cfg.validationLabels, subjectIndex, \
          cfg.output_classes, cfg.segmentation_dpatch, cfg.size_test_minibatches, logfile, epoch, cfg.saveSegmentation)
          print('Segmentation took {}'.format(round(time.time() - time_segmentation, 3)))
          my_logger('--------------- TEST EVALUATION ---------------', logfile)
          my_logger('          Full segmentation evaluation of subject' + str(subjectIndex), logfile)
          if dice_compare: 
            my_logger('DCS ' + str(dsc[-1]),logfile)
        my_logger('         FULL SEGMENTATION SUMMARY STATISTICS ', logfile)
        full_segm_DICE.append(np.mean(dsc))
        f = open(wd + '/Dice_scores.txt', 'w+')
        for d in full_segm_DICE:
          f.write(str(d))
          f.write('\n')
        f.close()
        my_logger('Overall DCS:   ' + str(full_segm_DICE[-1]),logfile)
        # Function to define if STOP flag goes to True or not, based on difference between last three or two segmentations.
        if len(full_segm_DICE) > 5:                        
          if np.max(np.abs(np.diff(full_segm_DICE[-4:] ))) < cfg.threshold_EARLY_STOP:
            my_logger('Convergence criterium met: No absolute change in DICE score over {} . Stopping training.'.format(cfg.threshold_EARLY_STOP),logfile)
            break    
        my_logger('###### SAVING TRAINED MODEL AT : ' + wd +'/Output/models/'+logfile[12:]+'.h5', logfile)
        model.save(wd+'/models/'+logfile[12:]+'_epoch' + str(epoch) + '.h5')           
        my_logger('Total training this epoch took ' + str(round(time.time()-t1,2)) + ' seconds',logfile)           
      # Save model if best results achieved
#      if len(full_segm_DICE) > 0:
#        if np.max(full_segm_DICE) == full_segm_DICE[-1]:
#          my_logger('###### SAVING TRAINED MODEL AT : ' + wd +'/Output/models/'+logfile[12:]+'.h5', logfile)
#          model.save(wd+'/models/'+logfile[12:]+'BestModel.h5')
#          best_results_epoch = epoch
#      if epoch == best_results_epoch + cfg.early_stop_patience:
#          my_logger('Convergence criterium met: No DICE score performance increment in the last {} epochs. Stopping training.'.format(cfg.early_stop_patience),logfile)
#          break      
      if dice_compare == False:
        my_logger('###### SAVING TRAINED MODEL AT : ' + wd +'/Output/models/'+logfile[12:]+'.h5', logfile)
        model.save(wd+'/models/'+logfile[12:]+'_epoch' + str(epoch) + '.h5')
            
      #################################################################################################
      #                                                                                               #
      #                               Training and Validation                                         #
      #                                                                                               #
      #################################################################################################   
      for i in range(0, cfg.num_iter):
        my_logger("                   Batch " + str(i+1) + "/" + str(cfg.num_iter) ,logfile)
        my_logger("###################################################### ",logfile)
        if not cfg.quickmode:
          ####################### VALIDATION ON BATCHES ############################                      
          with open(cfg.validationChannels[0]) as vl:
              n_valSubjects = len(vl.readlines())
          if cfg.n_subjects_val > n_valSubjects:
              print("Given number of subjects for validation set (" + str(cfg.n_subjects_val) +") is larger than the amount of \
              subjects in validation set (" +str(n_valSubjects)+ ")")
              cfg.n_subjects_val = n_valSubjects
              print('Using {} number of validation subjects'.format(n_valSubjects))
              
          valbatch, vallabels, val_TPM_patches = sampleTrainData(cfg.validation_TPM_channel, cfg.validationChannels, cfg.validationLabels, cfg.n_patches_val, cfg.n_subjects_val, \
          cfg.dpatch, cfg.output_classes, cfg.samplingMethod_val,cfg.model,cfg.data_augmentation, logfile, cfg.proportion_to_flip)
           
          val_performance.append(validation_on_batch(model, valbatch, val_TPM_patches, cfg.size_minibatches_val, vallabels, cfg.output_classes, logfile, cfg.model))               
          del valbatch, vallabels, val_TPM_patches
                      
        ####################### TRAINING ON BATCHES ##############################
        with open(cfg.trainLabels) as vl:
            n_trainSubjects = len(vl.readlines())                
        if cfg.n_subjects > n_trainSubjects:
            print("Given number of subjects for training set (" + str(cfg.n_subjects) +") is larger than the amount of \
            subjects in training set (" +str(n_trainSubjects)+ ")")
            cfg.n_subjects = n_trainSubjects
            print('Using {} number of training subjects'.format(n_trainSubjects))

        batch, labels, TPM_patches = sampleTrainData(cfg.TPM_channel, cfg.trainChannels,cfg.trainLabels, cfg.n_patches, cfg.n_subjects, cfg.dpatch, cfg.output_classes, \
        cfg.samplingMethod_train,cfg.model,cfg.data_augmentation, logfile, cfg.proportion_to_flip)
        assert not np.any(np.isnan(batch)), 'nan found in the input data!'          
        
        shuffleOrder = np.arange(batch.shape[0])
        np.random.shuffle(shuffleOrder)
        batch = batch[shuffleOrder]
        labels = labels[shuffleOrder]      
        if len(TPM_patches) > 0:
          TPM_patches = TPM_patches[shuffleOrder]
          
        del shuffleOrder  
        train_model_on_batch(model,batch,TPM_patches,labels,cfg.size_minibatches,history,losses,metrics,logfile,cfg.output_classes, cfg.model)  
        del batch, labels, TPM_patches 



def segment(configFile,workingDir):

    path = '/'.join(configFile.split('/')[:-1])
    configFileName = configFile.split('/')[-1][:-3]   
    sys.path.append(path)
    cfg = __import__(configFileName)
           
    start_epoch = int(cfg.path_to_model.split('.')[-2][cfg.path_to_model.split('.')[-2].find('epoch') + 5 : ]) + 1
    workingDir = workingDir.replace('\\','/') #for windows
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

    if len(cfg.TPM_channel) > 0:
      cfg.TPM_channel = workingDir + cfg.TPM_channel
    cfg.segmentChannels = [workingDir + x for x in cfg.segmentChannels]
    if len(cfg.segmentLabels) > 0:

        cfg.segmentLabels = workingDir + cfg.segmentLabels
        dice_compare = True
    else:
        dice_compare = False
    from keras.models import load_model   
    if cfg.output_classes == 6:
        try:
    	    from multiscale_CNN_TPM import Generalised_dice_coef_multilabel6, dice_coef_multilabel0,dice_coef_multilabel1,dice_coef_multilabel2,dice_coef_multilabel3,dice_coef_multilabel4,dice_coef_multilabel5
    	    my_custom_objects = {'Generalised_dice_coef_multilabel6':Generalised_dice_coef_multilabel6,
    				     'dice_coef_multilabel0':dice_coef_multilabel0,
    				     'dice_coef_multilabel1':dice_coef_multilabel1,
    				     'dice_coef_multilabel2':dice_coef_multilabel2,
    				     'dice_coef_multilabel3':dice_coef_multilabel3,
    				     'dice_coef_multilabel4':dice_coef_multilabel4,
    				     'dice_coef_multilabel5':dice_coef_multilabel5}
    		#custom_metrics =[dice_coef_multilabel6,dice_coef_multilabel0,dice_coef_multilabel1,dice_coef_multilabel2,dice_coef_multilabel3,dice_coef_multilabel4,dice_coef_multilabel5]
    		#my_custom_objects = dict(zip(np.sort(my_custom_objects.keys()), custom_metrics))
    
        except:
    	    from multiscale_CNN_TPM import w_dice_coef_multilabel6, dice_coef_multilabel0,dice_coef_multilabel1,dice_coef_multilabel2,dice_coef_multilabel3,dice_coef_multilabel4,dice_coef_multilabel5
    	    my_custom_objects = {'w_dice_coef_multilabel6':w_dice_coef_multilabel6,
    					     'dice_coef_multilabel0':dice_coef_multilabel0,
    					     'dice_coef_multilabel1':dice_coef_multilabel1,
    					     'dice_coef_multilabel2':dice_coef_multilabel2,
    					     'dice_coef_multilabel3':dice_coef_multilabel3,
    					     'dice_coef_multilabel4':dice_coef_multilabel4,
    					     'dice_coef_multilabel5':dice_coef_multilabel5}
        model = load_model(cfg.path_to_model, custom_objects = my_custom_objects )

    elif cfg.output_classes == 2:
        try:
            from multiscale_CNN_TPM import Generalised_dice_coef_multilabel2, dice_coef_multilabel0,dice_coef_multilabel1
            my_custom_objects = {'Generalised_dice_coef_multilabel2':Generalised_dice_coef_multilabel2,
				     'dice_coef_multilabel0':dice_coef_multilabel0,
				     'dice_coef_multilabel1':dice_coef_multilabel1}
        except:
            from multiscale_CNN_TPM import w_dice_coef_multilabel2, dice_coef_multilabel0,dice_coef_multilabel1
            my_custom_objects = {'w_dice_coef_multilabel2':w_dice_coef_multilabel2,
				     'dice_coef_multilabel0':dice_coef_multilabel0,
				     'dice_coef_multilabel1':dice_coef_multilabel1}
        model = load_model(cfg.path_to_model, custom_objects = my_custom_objects )

    elif cfg.output_classes == 7:
        try:
    	    from BIG_multiscale_CNN_TPM import Generalised_dice_coef_multilabel7, dice_coef_multilabel0,dice_coef_multilabel1,dice_coef_multilabel2,dice_coef_multilabel3,dice_coef_multilabel4,dice_coef_multilabel5, dice_coef_multilabel6
    	    my_custom_objects = {'Generalised_dice_coef_multilabel7':Generalised_dice_coef_multilabel7,
    				     'dice_coef_multilabel0':dice_coef_multilabel0,
    				     'dice_coef_multilabel1':dice_coef_multilabel1,
    				     'dice_coef_multilabel2':dice_coef_multilabel2,
    				     'dice_coef_multilabel3':dice_coef_multilabel3,
    				     'dice_coef_multilabel4':dice_coef_multilabel4,
    				     'dice_coef_multilabel5':dice_coef_multilabel5,
    				     'dice_coef_multilabel6':dice_coef_multilabel6}
        except:
    	    from multiscale_CNN_TPM import Generalised_dice_coef_multilabel7, dice_coef_multilabel0,dice_coef_multilabel1,dice_coef_multilabel2,dice_coef_multilabel3,dice_coef_multilabel4,dice_coef_multilabel5, dice_coef_multilabel6
    	    my_custom_objects = {'Generalised_dice_coef_multilabel7':Generalised_dice_coef_multilabel7,
    				     'dice_coef_multilabel0':dice_coef_multilabel0,
    				     'dice_coef_multilabel1':dice_coef_multilabel1,
    				     'dice_coef_multilabel2':dice_coef_multilabel2,
    				     'dice_coef_multilabel3':dice_coef_multilabel3,
    				     'dice_coef_multilabel4':dice_coef_multilabel4,
    				     'dice_coef_multilabel5':dice_coef_multilabel5,
    				     'dice_coef_multilabel6':dice_coef_multilabel6}

    print(my_custom_objects)
    model = load_model(cfg.path_to_model, custom_objects = my_custom_objects )


    full_segm_DICE = []
    np.set_printoptions(precision=3)

    print("------------------------------------------------------")
    print("                 FULL HEAD SEGMENTATION")
    print("------------------------------------------------------")

    dsc = []
    subjectIndex = 0
    epoch = start_epoch
    with open(cfg.segmentChannels[0]) as vl:
        n_segmentSubjects = len(vl.readlines())
    if cfg.test_subjects > n_segmentSubjects:
        print("Given number of subjects for test set (" + str(cfg.test_subjects) +") is larger than the amount of \
        subjects in test set (" +str(n_segmentSubjects)+ ")")
        cfg.test_subjects = n_segmentSubjects
    print('Using {} number of test subjects'.format(n_segmentSubjects))
    list_subjects_fullSegmentation = range(cfg.test_subjects)
    for subjectIndex in list_subjects_fullSegmentation:        
        fullHeadSegmentation(wd, cfg.penalty_MATRIX, cfg.TPM_channel, dice_compare, dsc, model, cfg.model, cfg.segmentChannels, cfg.segmentLabels, subjectIndex, \
        cfg.output_classes,cfg.segmentation_dpatch, cfg.size_test_minibatches, logfile, epoch, cfg.saveSegmentation)
        my_logger('--------------- TEST EVALUATION ---------------', logfile)
        my_logger('          Full segmentation evaluation of subject' + str(subjectIndex), logfile)
        if dice_compare: my_logger('DCS ' + str(dsc[-1]),logfile)
            
    my_logger('         FULL SEGMENTATION SUMMARY STATISTICS ', logfile)
    full_segm_DICE.append(np.mean(dsc))
    my_logger('Overall DCS:   ' + str(full_segm_DICE[-1]),logfile)
    

############################## AUXILIARY FUNCTIONS ############################################

def file_len(fname):
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i + 1