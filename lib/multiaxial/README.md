# multiaxial_brain_segmenter

Main script: Segment_MRI.py

# INPUT

SUBJECT_PATH:
Path to T1w MRI to segment in format nifti file (.nii )

SEGMENTATION_PATH:
[optional] segmentation of MRI to use as ground truth

OUTPUT_PATH:
Path to folder to which output segmentations will be saved (it will create folder if nonexistent)

SAGITTAL_MODEL_SESSION_PATH:
Path to sagittal segmenter model (included in repository), in format .h5

AXIAL_MODEL_SESSION_PATH:
Path to axial segmenter model (included in repository), in format .h5

CORONAL_MODEL_SESSION_PATH:
Path to coronal segmenter model (included in repository), in format .h5

CONSENSUS_LAYER_PATH:
Path to layer for merging output probabilities from 3 uniaxial models (included in repository), in format .h5

# OUTPUT

In OUTPUT_PATH you will find 4 files: Individual segmentations along each axis by each respective model, and the consensus segmentation made by majority vote. All segmentations will be in nifti format (.nii) with the header of the original MRI.

Segmented tissues and values:

0 - Background

1 - Air (Sinus cavities)

2 - Gray Matter

3 - White Matter

4 - Cerebrospinal Fluid

5 - Bone (Skull)

6 - Skin


Tested on Tensorflow 2 (2.0.0) + Python 3 and Tensorflow 1 (1.14) + Python 2

