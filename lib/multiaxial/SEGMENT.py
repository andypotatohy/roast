# -*- coding: utf-8 -*-
"""
Created on Wed Oct  4 09:57:23 2023

@author: Andrew
"""
import os
import sys
import nibabel as nib
import numpy as np
import tensorflow as tf
GPU = 0
gpus = tf.config.experimental.list_physical_devices('GPU')
if gpus:
  # Restrict TensorFlow to only use the first GPU
  try:
    tf.config.experimental.set_visible_devices(gpus[GPU], 'GPU')
    tf.config.experimental.set_memory_growth(gpus[GPU], True)
  except RuntimeError as e:
    # Visible devices must be set at program startup
    print(e)
print("Num GPUs Available: ", len(tf.config.list_physical_devices('GPU')))

def main():
    # print("Starting main function")
    if len(sys.argv) < 2:
        print('Please include an input image:')
        print('>>> python SEGMENT.py </subject_value> \n')
        sys.exit()

    working_dir = os.getcwd().replace('\\', '/')
    # Extract subj from command-line arguments
    subj = sys.argv[1].replace('\\', '/')
    # print(f"Subject: {subj}")
    
    
    if not os.path.isabs(subj):
        subj = os.path.join(working_dir, subj)
    OUTPUT_PATH = os.path.join(os.path.dirname(subj), '')
    scan_path = subj
    # print(f"Output path: {OUTPUT_PATH}")
    # print(f"Scan path: {scan_path}")
    SAGITTAL_MODEL_SESSION_PATH = 'sagittal_model.h5'
    AXIAL_MODEL_SESSION_PATH = 'axial_model.h5'
    CORONAL_MODEL_SESSION_PATH = 'coronal_model.h5'
    CONSENSUS_LAYER_PATH = 'consensus_layer.h5'

    # print(f"Sagittal model session path: {SAGITTAL_MODEL_SESSION_PATH}")
    # print(f"Axial model session path: {AXIAL_MODEL_SESSION_PATH}")
    # print(f"Coronal model session path: {CORONAL_MODEL_SESSION_PATH}")
    # print(f"Consensus layer path: {CONSENSUS_LAYER_PATH}")
    
    
    segmentation_path = None 
    anterior_commissure = None 
    save_segmentation = False
    os.chdir(os.path.dirname(os.path.abspath(__file__)))
    # os.chdir(r'C:\Users\Andrew\Documents\roast_multiaxial\lib\multiaxial')
    from preprocessing_lib import preprocess_head_MRI, reshape_back_to_original
    from utils import segment_MRI, Generalised_dice_coef_multilabel7, dice_coef_multilabel_bin0, dice_coef_multilabel_bin1, dice_coef_multilabel_bin2, dice_coef_multilabel_bin3, dice_coef_multilabel_bin4, dice_coef_multilabel_bin5, dice_coef_multilabel_bin6

    my_custom_objects = {
        'Generalised_dice_coef_multilabel7': Generalised_dice_coef_multilabel7,
        'dice_coef_multilabel_bin0': dice_coef_multilabel_bin0,
        'dice_coef_multilabel_bin1': dice_coef_multilabel_bin1,
        'dice_coef_multilabel_bin2': dice_coef_multilabel_bin2,
        'dice_coef_multilabel_bin3': dice_coef_multilabel_bin3,
        'dice_coef_multilabel_bin4': dice_coef_multilabel_bin4,
        'dice_coef_multilabel_bin5': dice_coef_multilabel_bin5,
        'dice_coef_multilabel_bin6': dice_coef_multilabel_bin6
    }

    if SAGITTAL_MODEL_SESSION_PATH is not None:
        print(f"Loading sagittal model from {SAGITTAL_MODEL_SESSION_PATH}")
        model_sagittal = tf.keras.models.load_model(SAGITTAL_MODEL_SESSION_PATH, custom_objects=my_custom_objects)
    if AXIAL_MODEL_SESSION_PATH is not None:
        print(f"Loading axial model from {AXIAL_MODEL_SESSION_PATH}")
        model_axial = tf.keras.models.load_model(AXIAL_MODEL_SESSION_PATH, custom_objects=my_custom_objects)
    if CORONAL_MODEL_SESSION_PATH is not None:
        print(f"Loading coronal model from {CORONAL_MODEL_SESSION_PATH}")
        model_coronal = tf.keras.models.load_model(CORONAL_MODEL_SESSION_PATH, custom_objects=my_custom_objects)

    # model_coronal.summary()
    print(f"Loading consensus model from {CONSENSUS_LAYER_PATH}")
    consensus_model = tf.keras.models.load_model(CONSENSUS_LAYER_PATH)
    # consensus_model.summary()

    print(f"Loading scan from {scan_path}")
    nii = nib.load(scan_path)
    if segmentation_path is not None:
        print(f"Loading segmentation from {segmentation_path}")
        nii_seg = nib.load(segmentation_path)
    else:
        nii_seg = None

    # print(f"Scan shape: {nii.shape}")

    subject = scan_path.split('/')[-1].replace('.nii', '')
    # subject = scan_path.split('/')[-2]
    # if subject.startswith('r'):
    #    subject = subject[1:]
    # print(f"Subject ID: {subject}")
     
    nii_out, nii_seg_out, coords, anterior_commissure, reconstruction_parms = preprocess_head_MRI(
        nii, nii_seg, anterior_commissure=anterior_commissure, keep_parameters_for_reconstruction=True)

    segmentation = segment_MRI(nii_out.get_fdata(), coords, model_sagittal, model_axial, model_coronal, consensus_model)
    # print("Reshaping segmentation back to original")
    nii_model_seg_reconstructed = reshape_back_to_original(segmentation, nii, reconstruction_parms, resample_order=0).get_fdata()

    img_ROAST = np.array(nii_model_seg_reconstructed, dtype='uint8')
  
    # Define mapping from MP to ROAST labels
    mapping = {1: 6, 3: 1, 4: 3, 5: 4, 6: 5}

    # Apply the mapping
    for mp_label, roast_label in mapping.items():
        img_ROAST[nii_model_seg_reconstructed == mp_label] = roast_label
        
    nii_out_pred = nib.Nifti1Image(img_ROAST, nii.affine)
    img_multiaxial = OUTPUT_PATH + subject + '_multiaxial_masks.nii'
    nib.save(nii_out_pred,  img_multiaxial)
    print("Segmentation saved successfully")

if __name__ == "__main__":
    main()
