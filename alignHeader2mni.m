function alignHeader2mni(T1,T2,spmOut,segOut)
% alignHeader2mni(T1,T2,spmOut,segOut)
% 
% Set the header info such that the voxel-to-world mapping maps directly
% into the MNI space.
% 
% (c) Andrew Birnbaum, Parra Lab at CCNY
%     Yu (Andy) Huang
% April 2024

[dirname,t1name,ext] = fileparts(T1);
if isempty(dirname), dirname = pwd; end

[~,spmOutName] = fileparts(spmOut);
seg8 = load([dirname filesep spmOutName '_seg8.mat']);
seg8Affine = seg8.Affine;

if ~exist([dirname filesep t1name '_MNI' ext],'file')
    disp('Aligning T1 MRI''s header to the MNI space ...');
    [Vol,info] = update_affine(T1,seg8Affine);
    niftiwrite(Vol,[dirname filesep t1name '_MNI' ext],info);
else
    disp('T1 MRI''s header already aligned to the MNI space. Skipping ...');
end

if ~isempty(T2)
    [dirname2,t2name,ext2] = fileparts(T2);
    if isempty(dirname2), dirname2 = pwd; end
    if ~exist([dirname2 filesep t2name '_MNI' ext2],'file')
        disp('Aligning T2 MRI''s header to the MNI space ...');
        [Vol,info] = update_affine(T2,seg8Affine);
        niftiwrite(Vol,[dirname2 filesep t2name '_MNI' ext2],info);
    else
        disp('T2 MRI''s header already aligned to the MNI space. Skipping ...');
    end
end

[~,segOutName] = fileparts(segOut);
if ~exist([dirname filesep segOutName '_masks_MNI.nii'],'file')
    disp('Aligning segmentation''s header to the MNI space ...');
    [Vol,info] = update_affine([dirname filesep segOutName '_masks.nii'],seg8Affine);
    niftiwrite(Vol,[dirname filesep segOutName '_masks_MNI.nii'],info);
else
    disp('Segmentation''s header already aligned to the MNI space. Skipping ...');
end

function [Vol,info] = update_affine(subj,seg8Affine)
info=niftiinfo(subj);
Vol = niftiread(subj);
Affine = info.Transform.T;
newAffine = seg8Affine*Affine';
info.Transform.T = newAffine';
%Makes sure both MRI and Segmentation are the same form
info.TransformName = 'Sform';