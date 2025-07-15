function alignHeader2mni(T1,T2,multiaxial,mniAffine)
% alignHeader2mni(T1,T2,spmOut,segOut)
% 
% Set the header info such that the voxel-to-world mapping maps directly
% into the MNI space.
%  
% (c) Andrew Birnbaum, Parra Lab at CCNY
%    Yu (Andy) Huang
% June 2025

[dirname,subject,ext] = fileparts(T1);
if isempty(dirname), dirname = pwd; end
% segOutName = '_multiaxial_manual'; %this is just for debugging 
if ~exist('mniAffine','var')
    if multiaxial
        %NiftyReg
       segOutName = '_multiaxial';
       mniAffine = runNiftyReg(T1);
    else
        %SPM
       seg8 = load([dirname filesep subject '_T1orT2_seg8.mat']);
       segOutName = '_SPM';
       mniAffine = seg8.Affine;
    end
end

if ~exist([dirname filesep subject segOutName '_MNI' ext],'file')
    disp('Aligning T1 MRI''s header to the MNI space ...');
    [Vol,info] = update_affine(T1,mniAffine);
    niftiwrite(Vol,[dirname filesep subject segOutName '_MNI' ext],info);
else
    disp('T1 MRI''s header already aligned to the MNI space. Skipping ...');
end

if ~isempty(T2)
    [dirname2,t2name,ext2] = fileparts(T2);
    if isempty(dirname2), dirname2 = pwd; end
    if ~exist([dirname2 filesep t2name segOutName '_MNI' ext2],'file')
        disp('Aligning T2 MRI''s header to the MNI space ...');
        [Vol,info] = update_affine(T2,mniAffine);
        niftiwrite(Vol,[dirname2 filesep t2name segOutName '_MNI' ext2],info);
    else
        disp('T2 MRI''s header already aligned to the MNI space. Skipping ...');
    end
end

if ~exist([dirname filesep subject '_T1orT2' segOutName '_masks_MNI.nii'],'file')
    disp('Aligning segmentation''s header to the MNI space ...');
    [Vol,info] = update_affine([dirname filesep subject '_T1orT2' segOutName '_masks.nii'],mniAffine);
    niftiwrite(Vol,[dirname filesep subject '_T1orT2' segOutName '_masks_MNI.nii'],info);
else
    disp('Segmentation''s header already aligned to the MNI space. Skipping ...');
end

function [Vol,info] = update_affine(subj,mniAffine)
info=niftiinfo(subj);
Vol = niftiread(subj);
Affine = info.Transform.T;
newAffine = mniAffine*Affine';
info.Transform.T = newAffine';
%Makes sure both MRI and Segmentation are the same form
info.TransformName = 'Sform';