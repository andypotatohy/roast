## -*- coding: utf-8 -*-
"""
Created on Mon Oct  2 14:59:07 2023

@author: Andrew
"""
import os
import sys
import nibabel as nib
import numpy as np
import random
from numpy.random import seed
from keras.utils import to_categorical
import tensorflow 
from skimage.transform import resize

seed(1)
tensorflow.random.set_seed(2)

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
    
def percentile95_normalizeMRI(data):
    p95 = np.percentile(data,95)
    data = data/p95
    return(data)    
    


################################ SAMPLING FUNCTIONS ####################################################

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

    for i in range(len(voxelCoordinates)):
        n_patches += len(voxelCoordinates[i])
    vol = np.zeros((n_patches,output_dpatch,output_dpatch,output_dpatch,6),dtype='float32')            
    for i in range(len(subjectIndexes)):  # same as iterating over subjects_TPM
        
        proxy_label = nib.load(TPM_channel)
        label_data = np.array(proxy_label.get_fdata(),dtype='float32')

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
   # subjects = getSubjectsToSample(channel, subjectIndexes)
    subjects = channel
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
       # subject = str(subjects[i])[:-1]
        #print('Subject with path: {}'.format(subject))
        proxy_img = nib.load(subjects)            
        img_data = np.array(proxy_img.get_fdata(),dtype='float32')
        #if percentile_normalization:
        img_data = percentile95_normalizeMRI(img_data)
        #else:
        #img_data = normalizeMRI(img_data)      
        if not np.isfinite(img_data).all():
        #if np.any(np.isnan(img_data)):
          print('Normalization: Nans found in scan {}'.format(subjects))
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
    
def sampleTestData(TPM_channel, testChannels, testLabels, subjectIndex, output_classes, output_dpatch,logfile):
   
    proxy_img = nib.load(testChannels[0])
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


def fullHeadSegmentation(wd,subj, penalty_MATRIX, TPM_channel, dice_compare, dsc, model, model_name, testChannels, testLabels, subjectIndex, output_classes, segmentation_dpatch, size_minibatches,logfile, epoch, saveSegmentation = False):    

    if ('Baseline' in model_name) or ('MultiPriors' in model_name) or ('DeepMedic' in model_name) or ('Vanilla' in model_name):
      output_dpatch = segmentation_dpatch - 48
    else:
      output_dpatch = segmentation_dpatch - 52    

    if len(testLabels) == 0:
        dice_compare = False    
    subjectIndex = [subjectIndex]
   # flairCh = getSubjectsToSample(testChannels[0], subjectIndex)
    flairCh = testChannels[0]
    #subID = flairCh[0].split('.')[0].split('/')[-1][:7]
    subID = flairCh.split('/')[-1].split('.')[0]
    #segmentationName = 'predictions/' + subID + '_segmentation_epoch' + str(epoch)
    #output = wd +'/' + segmentationName + '.nii.gz'
    #output = subj + subID + '_segmentation_epoch' + str(epoch)+ '.nii.gz'
    output = subj + subID + '_masks.nii.gz'
    if os.path.exists(output):
      print('Segmentation already found in: {}\nSkip..'.format(output))
      return

    # Open subject MRI to extract patches dynamically
    num_channels = 1
    proxy_img = nib.load(testChannels[0])
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
      
        dsc.append(score[0])
        print(dsc[-1])
        print('per class dice score: {}'.format(score[1]))
        print('mean DCS so far:' + str(np.mean(dsc)))
        
    if(saveSegmentation):

        nib.save(img, output)
        print('Saved segmentation of subject at: ' + output)
    
       
############################## create or load model ###########################################



def segment(configFile,workingDir,subj):
  
    configFile = configFile.replace('\\','/')
    path = '/'.join(configFile.split('/')[:-1])+ '/scripts'
    configFileName = configFile.split('/')[-1][:-3]   
    sys.path.append(path)
    subj = os.path.dirname(subj)+'/'
    cfg = __import__(configFileName)
           
    start_epoch = int(cfg.path_to_model.split('.')[-2][cfg.path_to_model.split('.')[-2].find('epoch') + 5 : ]) + 1
    wd = workingDir.replace('\\','/') #for windows
    logfile = 'segmentations.log'

    if len(cfg.TPM_channel) > 0:
        cfg.TPM_channel = cfg.TPM_channel
    cfg.segmentChannels = [cfg.segmentChannels] 
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
          
    fullHeadSegmentation(wd, subj, cfg.penalty_MATRIX, cfg.TPM_channel, dice_compare, dsc, model, cfg.model, cfg.segmentChannels, cfg.segmentLabels, subjectIndex, \
    cfg.output_classes,cfg.segmentation_dpatch, cfg.size_test_minibatches, logfile, epoch, cfg.saveSegmentation)
      
    

############################## AUXILIARY FUNCTIONS ############################################

def file_len(fname):
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i + 1

