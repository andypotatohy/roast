function viewSeg(segOut,mri2mni)
% viewSeg(segOut,mri2mni)
% 
% (c) Yu (Andy) Huang
% yhuang16@citymail.cuny.edu
% July 2025

[dirname,segOutName] = fileparts(segOut);
data = load_untouch_nii([dirname filesep segOutName '_masks.nii']);
allMaskShow = data.img;

% More anatomical looking colormap
color_map = [
    0, 0, 0;                   % background: black
    1, 1, 1;                   % white matter: white
    0.7, 0.7, 0.7;             % gray matter: gray
    105/255, 175/255, 255/255; % CSF: blue
    241/255, 214/255, 145/255; % bone
    177/255, 122/255, 101/255; % skin
    0.6863, 0.8824, 0.6863;    % air cavities
    %0, 0, 139/255;             % gel&elec
];

sliceshow(allMaskShow,[],color_map,[],'Tissue index','Segmentation. Click anywhere to navigate.',[],mri2mni); drawnow