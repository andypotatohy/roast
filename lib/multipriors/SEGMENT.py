# -*- coding: utf-8 -*-
"""
Created on Wed Oct  4 09:57:23 2023

@author: Andrew
"""
import os
import sys
import nibabel as nib
import numpy as np

def main():
    if len(sys.argv) < 2:
        print('Please include a model configuration file:')
        print('>>> python SEGMENT2.py </path_to/segmentation_config.py> </subject_value> \n')
        sys.exit()

    working_dir = os.getcwd().replace('\\', '/')

    # Extract subj from command-line arguments
    subj = sys.argv[2].replace('\\', '/') if len(sys.argv) >= 3 else None

    from scripts.lib_seg import segment
    config_file = sys.argv[1].replace('\\', '/')
    segment(config_file, working_dir, subj)

    # Call process_segmentation after segment to avoid the 'NameError'
    process_segmentation(subj)

def process_segmentation(subj):
    # Extract the name and extension from the subject file
    name, _ = os.path.splitext(subj)
    segmentation_file = name + '_masks.nii.gz'

    # Check if the segmentation file exists
    if os.path.exists(segmentation_file):
        print("Found segmentation file:", segmentation_file)

        # Load the segmentation file
        nii = nib.load(segmentation_file)
        img_MP = nii.get_fdata()

        # Convert to uint8
        img_ROAST = np.array(img_MP, dtype='uint8')

        # Define mapping from MP to ROAST labels
        mapping = {1: 6, 3: 1, 4: 3, 5: 4, 6: 5}

        # Apply the mapping
        for mp_label, roast_label in mapping.items():
            img_ROAST[img_MP == mp_label] = roast_label

        output_filename = name + '_masks.nii'
        nib.save(nib.Nifti1Image(img_ROAST, nii.affine), output_filename)
        print("Processed segmentation file saved as:", output_filename)
    else:
        print("No segmentation file with the specified ending found in the directory.")

if __name__ == "__main__":
    main()
