function [Affine,mri2mni] = getAffine(spmOut,landmarks)
%
% Computes the affine transformation from subject MRI space to MNI space
% using provided anatomical landmarks and known coordinates in the MNI152
% template.
%
% INPUTS:
%   spmOut     - Full path to the subject’s MRI NIfTI file (e.g., 'subj.nii')
%   landmarks  - 14×3 matrix of landmark coordinates in subject space
%                Order must match the MNI landmark set in the code:
%                  [Nasion, Inion, Right Ear, Left Ear, Scalp Center,
%                   and 9 central electrode locations]
%
% PROCESS:
%   - Loads affine matrices of subject MRI and MNI152_T1_1mm_ras.nii
%   - Uses pseudo-inverse to calculate an affine matrix mapping subject
%     landmarks to MNI landmarks
%   - Constructs a final 4×4 affine matrix mapping subject voxel space
%     directly to MNI space
%
% OUTPUTS:
%   Affine   - 4×4 matrix converting subject voxel coordinates to MNI voxel space
%   mri2mni  - 4×4 matrix mapping subject world coordinates to MNI world coordinates
%
% NOTES:
%   - The landmarks are expected in millimeter coordinates (world space).
%   - The function assumes a specific reference image in the path:
%       'example\MNI152_T1_1mm_ras.nii'
%
% DEPENDENCIES:
%   - Requires SPM installed and on the MATLAB path (`spm_vol` used)
%
% EXAMPLE USAGE:
%   [Affine, mri2mni] = getAffine('subj.nii', landmarks_from_mri);
%
% See also: spm_vol, spm_get_space
%
% (c) Andrew Birnbaum, Parra Lab at CCNY  
% 
% June 2025

% MNI_RAS head 
% May need to run MNI through ROAST to get RAS version
tpm = spm_vol('example\MNI152_T1_1mm_ras.nii');
affine_TPM = tpm.mat;

mri = spm_vol(spmOut);
affine_mri = mri.mat;

% Based on MNI_RAS head 
landmarkInTPM =    [92   210    31 % nasion
    92    12    33 % inion
   169    94    23 % right
    16    94    23 % left
    92   111    76 % center of scalp
    % 9 center electrodes
    92     8    75
    92    17   112
    92    37   147
    92    71   170
    92   112   172
    92   150   161
    92   183   141
    92   206   110
    92   213    73];

landmarkInTPM = [landmarkInTPM,ones(14,1)];
landmarks2 = [landmarks,ones(14,1)];
A = transpose(landmarkInTPM) / transpose((landmarks2));

Affine  = affine_TPM *A*inv(affine_mri );
Affine(end,:) = [0 0 0 1];
mri2mni = Affine*affine_mri;
