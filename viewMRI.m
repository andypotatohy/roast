function viewMRI(T1,T2,mri2mni)
% viewMRI(T1,T2,mri2mni)
% 
% (c) Yu (Andy) Huang
% yhuang16@citymail.cuny.edu
% July 2025

t1Data = load_untouch_nii(T1);
sliceshow(t1Data.img,[],'gray',[],[],'MRI: Click anywhere to navigate.',[],mri2mni); drawnow

if ~isempty(T2)
    t2Data = load_untouch_nii(T2);
    sliceshow(t2Data.img,[],'gray',[],[],'MRI: T2. Click anywhere to navigate.',[],mri2mni); drawnow
end