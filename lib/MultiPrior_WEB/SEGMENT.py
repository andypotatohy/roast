# -*- coding: utf-8 -*-
"""
Created on Wed Oct  4 09:57:23 2023

@author: Andrew
"""
if __name__ == "__main__":
    import warnings
    warnings.filterwarnings("ignore")
    import os
    import sys
    sys.path.append(os.getcwd())
    sys.path.append(os.getcwd() + '/scripts')
 
    if len(sys.argv) < 2:
        print('Please include a model configuration file:')
        print('>>> python SEGMENT2.py </path_to/segmentation_config.py> </subject_value> \n')
        sys.exit()

    workingDir = os.getcwd().replace('\\','/')

    # Extract subj from command-line arguments
    if len(sys.argv) >= 3:
        subj = sys.argv[2].replace('\\','/')
          
    else:
        subj = None  # Provide a default value or handle the case when subj is not provided

    from scripts.lib_seg import segment
    configFile = sys.argv[1].replace('\\','/')
    segment(configFile, workingDir, subj)

import os
import nibabel as nib
import numpy as np

path = os.path.dirname(subj)

# List all files in the specified directory
all_files = os.listdir(path)

# Find the segmentation file with the ending "segmentation_epoch101.nii.gz"
segmentation_file = next((os.path.join(path, x) for x in all_files if x.endswith("segmentation_epoch101.nii.gz")), None).replace('\\','/')

if segmentation_file is not None:
    # The 'segmentation_file' variable contains the path to the segmentation file
    print("Found segmentation file:", segmentation_file)

    nii = nib.load(segmentation_file)

    #name = segmentation_file.split('_')[0] + '_1mm_MultiPriors_Segmentation.nii'
    name = subj.split('.')[0] + '_T1orT2_masks_MultiPrior_Segmentation.nii'
    print(name)

    if not os.path.exists(os.path.join(path, name)):
        img_MP = nii.get_fdata()

        np.unique(img_MP)

        # ROAST
        # 0 - Bkg
        # 1 - WM
        # 2 - GM
        # 3 - CSF
        # 4 - Bone
        # 5 - skin
        # 6 - sinus (air)

        # MP
        # 0 - Bkg
        # 1 - Sinus
        # 2 - GM
        # 3 - WM
        # 4 - CSF
        # 5 - Bone
        # 6 - Skin

        img_ROAST = np.array(img_MP, dtype='uint8')

        # MP --> ROAST
        img_ROAST[img_MP == 1] = 6  # Sinus (MP) --> Sinus (ROAST)
        img_ROAST[img_MP == 3] = 1
        img_ROAST[img_MP == 4] = 3
        img_ROAST[img_MP == 5] = 4
        img_ROAST[img_MP == 6] = 5

        np.unique(img_ROAST)

        nib.save(nib.Nifti1Image(img_ROAST, nii.affine), os.path.join(path, name))
        os.remove(segmentation_file)
    else:
        print("MultiPriors_Segmentation.nii already exists in the output directory. Skipping processing.")
else:
    print("No segmentation file with the specified ending found in the directory.")
