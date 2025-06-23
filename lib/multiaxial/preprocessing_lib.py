#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct  2 15:01:29 2024

@author: deeperthought
"""


import nibabel as nib
import numpy as np

from skimage.transform import resize

from nibabel.orientations import axcodes2ornt
from nibabel.orientations import ornt_transform
print("hello")

def reorient( nii, orientation) -> nib.Nifti1Image:
    """Reorients a nifti image to specified orientation. Orientation string or tuple
    must consist of "R" or "L", "A" or "P", and "I" or "S" in any order."""
    orig_ornt = nib.io_orientation(nii.affine)
    targ_ornt = axcodes2ornt(orientation)
    transform = ornt_transform(orig_ornt, targ_ornt)
    reoriented_nii = nii.as_reoriented(transform)
    return reoriented_nii


def create_coordinate_matrix(shape, anterior_commissure):
    x, y, z = shape
    meshgrid = np.meshgrid(np.linspace(0, x - 1, x), np.linspace(0, y - 1, y), np.linspace(0, z - 1, z), indexing='ij')
    coordinates = np.stack(meshgrid, axis=-1) - np.array(anterior_commissure)
    matrix_with_ones = np.concatenate([coordinates, np.ones((coordinates.shape[0], coordinates.shape[1], coordinates.shape[2], 1))], axis=-1)

    return matrix_with_ones



def preprocess_head_MRI(nii: nib.Nifti1Image, nii_seg: nib.Nifti1Image = None, anterior_commissure: tuple = None, keep_parameters_for_reconstruction: bool = False):
    """
    Preprocesses a head MRI image.
    
    Args:
      nii: A NIfTI image object representing the MRI scan.
      nii_seg: (Optional) A NIfTI image object representing the segmentation mask with 7 classes.
      anterior_commissure: (Optional) A tuple representing the coordinates of the anterior commissure. If None, it uses the center of the image.
    
    Returns:
      nii_out: The preprocessed MRI image as a NIfTI object.
      nii_seg_out: (Optional) The preprocessed segmentation mask as a NIfTI object, if provided.
      coords: A NumPy array containing the coordinates of the anterior commissure relative to the preprocessed image.
      anterior_commissure: The updated coordinates of the anterior commissure after preprocessing.
    """    
    if nii_seg is not None:
        assert nii.shape == nii_seg.shape
    if anterior_commissure is None:
        print('No anterior commissure location given.. centering to center of image..')
        anterior_commissure = [nii.shape[0]//2, nii.shape[1]//2, nii.shape[2]//2]
    else:
        print(f'anterior commissure given: {anterior_commissure}')
    orientation = nib.aff2axcodes(nii.affine)
    
    if ''.join(orientation) != 'RAS':
        print(f'Image orientation : {orientation}. Changing to RAS..')
        nii = reorient(nii, "RAS")
        if nii_seg is not None:
            nii_seg = reorient(nii_seg, "RAS")
   
    if nii_seg is not None:
        img_seg = nii_seg.get_fdata()
    else:
        print('Segmentation image not given')
                
    ############### ISOTROPIC #######################
    
    res = nii.header['pixdim'][1:4]
    img = nii.get_fdata()
    new_shape = np.array(np.array(nii.shape)*res, dtype='int')
    if np.any(np.array(nii.shape) != new_shape):
        img = resize(img, new_shape, anti_aliasing=True, preserve_range=True)
        if nii_seg is not None: img_seg = resize(img_seg, new_shape, order=0, anti_aliasing=True, preserve_range=True)
    
    
    nii_new = nib.Nifti1Image(nii.get_fdata().copy(), nii.affine, nii.header)
    
    nii_new.affine[0][0] = 1.
    nii_new.affine[1][1] = 1.
    nii_new.affine[2][2] = 1.
    
    nii_new.header['pixdim'][1:4] = np.diag(nii_new.affine)[0:3]
    
    ############### Crop/Pad to make shape 256,256,256 ###############
    
    d1, d2, d3 = new_shape
    start = None
    end = None    
    
    if d1 < 256:
        pad1 = 256-d1
        img = np.pad(img, ((pad1//2, pad1//2+pad1%2),(0,0),(0,0)))
        if nii_seg is not None: img_seg = np.pad(img_seg, ((pad1//2, pad1//2+pad1%2),(0,0),(0,0)))

        anterior_commissure[0] += pad1//2
        
    
    if d2 > 256: 
        crop2 = d2-256
        img = img[:,crop2//2:-(crop2//2+crop2%2)]
        if nii_seg is not None: img_seg = img_seg[:,crop2//2:-(crop2//2+crop2%2)]
        anterior_commissure[1] -= crop2//2
            
    elif d2 < 256:
        pad2 = 256-d2
        img = np.pad(img, ((0,0),(pad2//2, pad2//2+pad2%2),(0,0)))
        if nii_seg is not None: img_seg = np.pad(img_seg, ((0,0),(pad2//2, pad2//2+pad2%2),(0,0)))
        anterior_commissure[1] += pad2//2
        
    if d3 > 256: 

        #--- head start
        proj = np.max(img,(0,1))
        proj[proj < np.percentile(proj, 50) ] = 0
        proj[proj > 0] = 1
        end = np.max(np.argwhere(proj == 1))
        end = np.min([end + 20, d3]) # leave some space above the head
        start = end-256

        if start < 0:
            crop3 = d3 - 256
            img = img[:,:,crop3:]
            if nii_seg is not None: img_seg = img_seg[:,:,crop3:]
            anterior_commissure[2] -= crop3
         
        else:
            img = img[:,:,start:end]
            if nii_seg is not None: img_seg = img_seg[:,:,start:end]
            anterior_commissure[2] -= start

    elif d3 < 256:

        pad3 = 256-d3
        img = np.pad(img, ((0,0),(0,0),(pad3//2, pad3//2+pad3%2)))
        if nii_seg is not None: img_seg = np.pad(img_seg, ((0,0),(0,0),(pad3//2, pad3//2+pad3%2)))
        anterior_commissure[2] += pad3//2

    coords = create_coordinate_matrix(img.shape, anterior_commissure)        
    
    # Intensity normalization
    p95 = np.percentile(img, 95)
    img = img/p95
    
    coords = coords[:,:,:,:3]
    coords = coords/256.
    
                
    if nii_seg is not None:
        if keep_parameters_for_reconstruction:
            reconstruction_parms = d1,d2,d3,start,end
            return nib.Nifti1Image(img, nii_new.affine), nib.Nifti1Image(np.array(img_seg, dtype='int8'), nii_new.affine), np.array(coords, dtype='int16'), np.array(anterior_commissure, dtype='int'), reconstruction_parms
        else:
            return nib.Nifti1Image(img, nii_new.affine), nib.Nifti1Image(np.array(img_seg, dtype='int8'), nii_new.affine), np.array(coords, dtype='int16'), np.array(anterior_commissure, dtype='int')

    else:
        if keep_parameters_for_reconstruction:
            reconstruction_parms = d1,d2,d3,start,end
            return nib.Nifti1Image(img, nii_new.affine), None, np.array(coords, dtype='int16'), np.array(anterior_commissure, dtype='int'), reconstruction_parms
        else:
            return nib.Nifti1Image(img, nii_new.affine), None, np.array(coords, dtype='int16'), np.array(anterior_commissure, dtype='int')
        


def reshape_back_to_original(img, nii_original, reconstruction_parms, resample_order=1):
    d1,d2,d3, start, end = reconstruction_parms

    #f'{nii_out.shape}     --> (crop/pad)-->  {d1,d2,d3} --> (resample) --> {nii.shape} --> (reorient)'

    # pad or crop z axis
    if d3 > 256:     
        if start < 0:
            pad3 = d3 - 256
            img = np.pad(img, ((0,0),(0,0),(pad3,0)))#img[:,:,crop3:]
        else:
            pad_end = d3-end
            img = np.pad(img, ((0,0),(0,0),(start,pad_end))) # img[:,:,start:end]           
    elif d3 < 256:
        crop3 = 256-d3
        if crop3 == 1:
            img = img[:,:,crop3//2:]
    
        else:
            img = img[:,:,crop3//2:-crop3//2+crop3%2]
    
    # pad or crop y axis
    if d2 > 256:    
        pad2 = d2-256
        img = np.pad(img, ((0,0),(pad2//2,pad2//2+pad2%2),(0,0)))#img[:,:,crop3:]
            
    elif d2 < 256:
        crop2 = 256-d2
        if crop2 == 1:
            img = img[:,crop2//2:,:]
        
        else:
            img = img[:,crop2//2:-crop2//2+crop2%2,:]
        
    # pad or crop x axis
    if d1 < 256:
        crop1 = 256-d1
        img = img[crop1//2:-crop1//2+crop1%2,:,:]
        
        
    # resample to original resolution
    img = resize(img, nii_original.shape, anti_aliasing=True, preserve_range=True, order=resample_order)
    
    nii_out = nib.Nifti1Image(np.array(img, dtype='int16'), nii_original.affine)
    return nii_out
