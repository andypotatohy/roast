function landmarks = getLandmarks(P,T2)
% landmarks = getLandmarks(P,T2)
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

[dirname,baseFilename] = fileparts(P);
if isempty(T2)
    baseFilename = [baseFilename '_T1orT2'];
else
    baseFilename = [baseFilename '_T1andT2'];
end
load(fullfile(dirname,[baseFilename '_seg8.mat']),'image','tpm','Affine');
tpm2mri = inv(image(1).mat)*inv(Affine)*tpm(1).mat;
% mapping landmarks from TPM to individual MRI % ANDY 2017-05-17

landmarkInTPM = [61 139 98; % nasion
                 61 9 100; % inion
                 11 62 93; % right
                 111 63 93; % left; note here because eTPM is LAS orientation
                 61 113 7; % front_neck
                 61 7 20]; % back_neck

landmarks = zeros(size(landmarkInTPM,1),size(landmarkInTPM,2));
for i=1:size(landmarkInTPM,1)
    
    temp = tpm2mri * [landmarkInTPM(i,:) 1]';
    landmarks(i,:) = round(temp(1:3)');
    
end