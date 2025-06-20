function [landmarks,mri2mni] = getLandmarks(spmOut,multiaxial)

% landmarks = getLandmarks(spmOut)
%
% Mapping landmarks in the TPM (hard-coded) to the individual head, by
% using the mapping computed during the segmentation process in SPM.
% 
% Landmarks follow the order of: nasion, inion, right, left, front neck,
% and back neck.
% 
% (c) Yu (Andy) Huang, Parra Lab at CCNY
% yhuang16@citymail.cuny.edu
% April 2018 
% (c) Andrew Birnbaum, Parra Lab at CCNY 
% May 2025

if multiaxial
    % mapping landmarks from TPM to individual MRI % ANDREW 2025-05-11
    % NIFTYREG 
    tpm = spm_vol('eTPM.nii');
    affine_TPM = tpm.mat;
    mri = spm_vol(spmOut);
    affine_mri = mri.mat;
    Affine = runNiftyReg(spmOut);
    tpm2mri = inv(affine_mri)*inv(Affine)*affine_TPM;
else
    % SPM 
    [dirname,spmOutName] = fileparts(spmOut);
    if isempty(dirname), dirname = pwd; end
    load([dirname filesep spmOutName '_T1orT2_seg8.mat'],'image','tpm','Affine');
    tpm2mri = inv(image(1).mat)*inv(Affine)*tpm(1).mat;
    mri2mni = Affine*image(1).mat;

end
% mapping landmarks from TPM to individual MRI % ANDY 2017-05-17

landmarkInTPM = [61 139 98; % nasion
                 61 9 100; % inion
                 11 62 93; % right
                 111 63 93;]; % left; note here because eTPM is LAS orientation
                 % 61 113 7; % front_neck
                 % 61 7 20]; % back_neck

landmarks = zeros(size(landmarkInTPM,1),size(landmarkInTPM,2));
for i=1:size(landmarkInTPM,1)

    temp = tpm2mri * [landmarkInTPM(i,:) 1]';
    landmarks(i,:) = round(temp(1:3)');

end

if ~exist('mri2mni', 'var')
% Affine = inv(affine_mri*tpm2mri*inv(affine_TPM)); already have Affine but
% this is to check that they are equal
mri2mni = Affine*affine_mri;
end


end