function segTouchup(P,T2,isSmooth,conn)
% segTouchup(P,T2,isSmooth,conn)
%
% This function combines mysegment() and autoPatching() in ROAST version
% 2.1 and earlier. It performs automated clean up on the output of SPM8
% New Segmentation and applies auto patch on the remaining holes.
%
% See Huang et al 2013 (DOI: 10.1088/1741-2560/10/6/066004), and Huang et
% al 2017 (DOI: https://doi.org/10.1101/217331) for details.
% 
%
% (c) Yu Huang (Andy), Chris Thomas, Soterix Medical Inc., June 2017
% (c) Yu Huang (Andy), Lucas C. Parra, Neural Engineering Lab, City College of New York, August 2011
% (c) Zhewei Jiang, Jacek Dmochowski, Neural Engineering Lab, City College of New York, July 2010
% Send bugs to yhuang16@citymail.cuny.edu

if nargin < 3 || isempty(isSmooth)
    isSmooth = 1; % smooth the tissue by default
end

if nargin < 4
    conn = 18;
end

disp('loading data...')
% cd(dirname)
[dirname,baseFilename] = fileparts(P);
if isempty(dirname), dirname = pwd; end

if isempty(T2)
    baseFilename = [baseFilename '_T1orT2'];
else
    baseFilename = [baseFilename '_T1andT2'];
end

gray = load_untouch_nii([dirname filesep 'c1' baseFilename '.nii']);
white = load_untouch_nii([dirname filesep 'c2' baseFilename '.nii']);
csf = load_untouch_nii([dirname filesep 'c3' baseFilename '.nii']);
bone = load_untouch_nii([dirname filesep 'c4' baseFilename '.nii']);
skin = load_untouch_nii([dirname filesep 'c5' baseFilename '.nii']);
air = load_untouch_nii([dirname filesep 'c6' baseFilename '.nii']);
% load the masks

gray_temp = gray.img; white_temp = white.img; csf_temp = csf.img;
bone_temp = bone.img; skin_temp = skin.img; air_temp = air.img;

if isSmooth
    disp('smoothing masks...')
%     disp('smoothing GM...')
    smt_fil = fspecial('gaussian', 5, 0.2);
    for i = 1:size(gray_temp,3)
        gray_temp(:,:,i) = imfilter(gray_temp(:,:,i),smt_fil);
    end
    
%     disp('smoothing WM...')
    smt_fil = fspecial('gaussian', 5, 0.1);
    for i = 1:size(white_temp,3)
        white_temp(:,:,i) = imfilter(white_temp(:,:,i),smt_fil);
    end
    
%     disp('smoothing CSF...')
    smt_fil = fspecial('gaussian', 5, 0.1);
    for i = 1:size(csf_temp,3)
        csf_temp(:,:,i) = imfilter(csf_temp(:,:,i),smt_fil);
    end
    
%     disp('smoothing bone...')
    smt_fil = fspecial('gaussian', 5, 0.4);
    for i = 1:size(bone_temp,3)
        bone_temp(:,:,i) = imfilter(bone_temp(:,:,i),smt_fil);
    end
    
%     disp('smoothing skin...')
    smt_fil = fspecial('gaussian', 5, 1);
    for i = 1:size(skin_temp,3)
        skin_temp(:,:,i) = imfilter(skin_temp(:,:,i),smt_fil);
    end
    
%     disp('smoothing air...')
    smt_fil = fspecial('gaussian', 5, 1);
    for i = 1:size(air_temp,3)
        air_temp(:,:,i) = imfilter(air_temp(:,:,i),smt_fil);
    end
end
% Smooth the mask if it's not clean enough

disp('creating binary masks...')
[empt_temp,gray_temp,white_temp,csf_temp,bone_temp,skin_temp,air_temp]...
    = binaryMaskGenerate(gray_temp,white_temp,csf_temp,bone_temp,skin_temp,air_temp);
% Create binary mask for each tissue class

disp('fixing CSF continuity...')
se=ones(3,3,3);
dcsf=imdilate(csf_temp, se);
dbone=imdilate(bone_temp, se);
contin=(empt_temp&dcsf)|(dbone&gray_temp);
csf_temp=csf_temp|contin;

[~,csf_temp,bone_temp,gray_temp]...
    = binaryMaskGenerate(csf_temp,bone_temp,gray_temp); % NOTE: no removal, then no empty voxel generated
% Fix csf continuity

disp('removing disconnected voxels...')
% disp('removing disconnected voxels for GM...')
siz = sizeOfObject(gray_temp,conn);
thres = siz(2)+1; % may not be robust, but will be tested in future data % ANDY 2017-05-17
gray_temp = bwareaopen(gray_temp,thres,conn);

% disp('removing disconnected voxels for WM...')
siz = sizeOfObject(white_temp,conn);
thres = siz(2)+1; % may not be robust, but will be tested in future data % ANDY 2017-05-17
white_temp = bwareaopen(white_temp,thres,conn);

% disp('removing disconnected voxels for CSF...')
siz = sizeOfObject(csf_temp,conn);
thres = siz(4)+1;
csf_temp = bwareaopen(csf_temp,thres,conn);

% disp('removing disconnected voxels for bone...')
siz = sizeOfObject(bone_temp,conn);
thres = siz(2)+1; % this is aggressive, removing all floating structures including those in the spine.
                  % But maybe it does not matter much, as we don't care much about the lower part of the head
                  % ANDY 2017-05-17
                  % If we do this less aggressively, ie, those wrong
                  % floatings will be kept, then the model may not be
                  % solved without manually removing them
% thres = 300; % this is hard to generalize robustly
bone_temp = bwareaopen(bone_temp,thres,conn);

% disp('removing disconnected voxels for skin...')
siz = sizeOfObject(skin_temp,conn);
thres = siz(2)+1;
skin_temp = bwareaopen(skin_temp,thres,conn);

% disp('removing disconnected voxels for air...')
thres = 20;
air_temp = bwareaopen(air_temp,thres);
% Identify the disconnected voxels

disp('generating and labeling empty voxels...')
empt_temp = binaryMaskGenerate(gray_temp,white_temp,csf_temp,bone_temp,skin_temp,air_temp);
% Generate unassigned voxels (empty voxels)

% for j = 1:2 % usually all empty voxels will be labelled in two loops
while any(empt_temp(:))
    gray_fil = uint8(gray_temp)*255; white_fil = uint8(white_temp)*255;
    csf_fil = uint8(csf_temp)*255; bone_fil = uint8(bone_temp)*255;
    skin_fil = uint8(skin_temp)*255; air_fil = uint8(air_temp)*255;
    sigma = [1 1 1 1 1 1];
    smt_fil = fspecial('gaussian', 5, sigma(1));
    for i = 1:size(gray_fil,3)
        img = gray_fil(:,:,i);
        gray_fil(:,:,i) = imfilter(img,smt_fil);
    end
    smt_fil = fspecial('gaussian', 5, sigma(2));
    for i = 1:size(white_fil,3)
        img = white_fil(:,:,i);
        white_fil(:,:,i) = imfilter(img,smt_fil);
    end
    smt_fil = fspecial('gaussian', 5, sigma(3));
    for i = 1:size(csf_fil,3)
        img = csf_fil(:,:,i);
        csf_fil(:,:,i) = imfilter(img,smt_fil);
    end
    smt_fil = fspecial('gaussian', 5, sigma(4));
    for i = 1:size(bone_fil,3)
        img = bone_fil(:,:,i);
        bone_fil(:,:,i) = imfilter(img,smt_fil);
    end
    smt_fil = fspecial('gaussian', 5, sigma(5));
    for i = 1:size(skin_fil,3)
        img = skin_fil(:,:,i);
        skin_fil(:,:,i) = imfilter(img,smt_fil);
    end
    smt_fil = fspecial('gaussian', 5, sigma(6));
    for i = 1:size(air_fil,3)
        img = air_fil(:,:,i);
        air_fil(:,:,i) = imfilter(img,smt_fil);
    end
    
    [~,air_fil,gray_fil,white_fil,csf_fil,bone_fil,skin_fil]...
        = binaryMaskGenerate(air_fil,gray_fil,white_fil,csf_fil,bone_fil,skin_fil);

    gray_temp = (empt_temp&gray_fil)|gray_temp;
    white_temp = (empt_temp&white_fil)|white_temp;
    csf_temp = (empt_temp&csf_fil)|csf_temp;
    bone_temp = (empt_temp&bone_fil)|bone_temp;
    skin_temp = (empt_temp&skin_fil)|skin_temp;
    air_temp = (empt_temp&air_fil)|air_temp;
    
    empt_temp = xor(empt_temp,((empt_temp&gray_fil) | (empt_temp&white_fil) | (empt_temp&csf_fil)...
        | (empt_temp&bone_fil) | (empt_temp&skin_fil) | (empt_temp&air_fil))); % update empty voxels
end
% Relabel each empty voxel to its nearest tissue type
% The Gaussian filter is used to calculate distances, and max operation
% relabels each empty voxel based on the distances.

disp('removing outside air...')
temp = xor(air_temp, ones(size(air_temp)));
temp = imerode(imfill(imclose(temp,ones(10,10,10)),'holes'),ones(12,12,12));
air_temp = air_temp & temp;
% Remove the air outside the head.
% This won't remove all the outside air if there is empty space outside the
% head. To avoid empty space outside the head, make sure the 2nd image (if
% available) used for segmentation does not have empty space (e.g. Dixon
% image from stroke data has empty space) % ANDY 2017-05-17

load([dirname filesep baseFilename '_rmask.mat'],'holes_vol','eyes_vol','WMexclude_vol');
% for exceptions in auto touchup

% assign labels to tissues in this order: white,gray,csf,bone,skin,air
% note white-gray order differs from SPM outputs. Changed after ROAST V2.1
allMask = zeros(size(white_temp));
allMask(white_temp) = 1;
allMask(gray_temp) = 2;
allMask(csf_temp) = 3;
allMask(bone_temp) = 4;
allMask(skin_temp) = 5;
allMask(air_temp) = 6;

% Create matrix of interested neighbor directions
[x1,x2,x3] = ndgrid(-1:1);
perms = [x1(:) x2(:) x3(:)];
vect_g = zeros(18,3);
h = 1;
for i = 1:length(perms)
    %     test_val = abs(perms(i,1)) + abs(perms(i,2)) + abs(perms(i,3));
    test_val = norm(perms(i,:),1);
    if test_val == 1 || test_val == 2 %to omit corners (3D vertices) and center of 3x3 area
        vect_g(h,:) = perms(i,:);
        h = h+1;
    end
end
clear i h

% % ==================Patching gray matter=======================
% Rule: convert white matter pixels to gray matter if they touch CSF/bone/skin/air
% ind_maskOfInterest = find(allMask==1);
% handle local specifics for WM here (there is no gray matter in brain stem and around ventricles)
disp('patching gray matter...');
WMall = allMask==1;
WMkeep = WMall & imdilate(WMexclude_vol>1e-4,ones(10,10,10)); % this part of WM won't be touched
% image dilation needed to make it work, sometimes dilated a lot
% dilation needs free para, ugly
ind_maskOfInterest = find(xor(WMall,WMkeep));
[i,j,k] = ind2sub(size(allMask),ind_maskOfInterest);

isInside = i>1 & i<size(allMask,1) & j>1 & j<size(allMask,2) & k>1 & k<size(allMask,3);
% only process those pixels that are not on the image edges
ind_maskOfInterest = ind_maskOfInterest(isInside);
i = i(isInside); j = j(isInside); k = k(isInside);

wrongNeighborInd = 3:6; % Andy 2017-09-13 added more wrong neighbors other than CSF
for ii = 1:length(vect_g) %cycle through all interested neighbors
    ind_maskOfIntestNeighbor = sub2ind(size(allMask),i+vect_g(ii,1),j+vect_g(ii,2),k+vect_g(ii,3));
    for t = 1:length(wrongNeighborInd)
        allMask(ind_maskOfInterest(find(allMask(ind_maskOfIntestNeighbor)==wrongNeighborInd(t)))) = 2;
        % now operation on the pixel itself, not its neighbor
    end
end

% % ==================Patching CSF=======================
% Rule: If GM or WM touches bone/skin/air, then convert the wrong
% neighbors to CSF
disp('patching CSF...');
ind_maskOfInterest = find(allMask==1 | allMask==2);
[i,j,k] = ind2sub(size(allMask),ind_maskOfInterest);

isInside = i>1 & i<size(allMask,1) & j>1 & j<size(allMask,2) & k>1 & k<size(allMask,3);
% only process those pixels that are not on the image edges
% ind_maskOfInterest = ind_maskOfInterest(isInside);
i = i(isInside); j = j(isInside); k = k(isInside);

wrongNeighborInd = 4:6;
for ii = 1:length(vect_g) %cycle through all interested neighbors
    ind_maskOfIntestNeighbor = sub2ind(size(allMask),i+vect_g(ii,1),j+vect_g(ii,2),k+vect_g(ii,3));
    for t = 1:length(wrongNeighborInd)
        allMask(ind_maskOfIntestNeighbor(find(allMask(ind_maskOfIntestNeighbor)==wrongNeighborInd(t)))) = 3;
    end
end

% % ==================Patching bone=======================
% Rule: If CSF touches skin/air, then convert the wrong
% neighbors to bone
disp('patching bone...');
ind_maskOfInterest = find(allMask==3);

% exclude the eyeballs
mask_eyes = eyes_vol>1e-4; %graythresh(eyes_vol(:));
% because tiny = 1e-5; in spm_load_mask.m
iseyes = mask_eyes(ind_maskOfInterest);
ind_maskOfInterest = ind_maskOfInterest(~iseyes);

[i,j,k] = ind2sub(size(allMask),ind_maskOfInterest);

isInside = i>1 & i<size(allMask,1) & j>1 & j<size(allMask,2) & k>1 & k<size(allMask,3);
% only process those pixels that are not on the image edges
% ind_maskOfInterest = ind_maskOfInterest(isInside);
i = i(isInside); j = j(isInside); k = k(isInside);

wrongNeighborInd = 5:6;
for ii = 1:length(vect_g) %cycle through all interested neighbors
    ind_maskOfIntestNeighbor = sub2ind(size(allMask),i+vect_g(ii,1),j+vect_g(ii,2),k+vect_g(ii,3));
    for t = 1:length(wrongNeighborInd)
        allMask(ind_maskOfIntestNeighbor(find(allMask(ind_maskOfIntestNeighbor)==wrongNeighborInd(t)))) = 4;
    end
end
% handle local specifics for bone here
allMask(find(allMask==4 & holes_vol>1e-4)) = 5;
% because tiny = 1e-5; in spm_load_mask.m

% save out
% disp('Saving...');
white.img = uint8(allMask);
white.hdr.dime.scl_slope=1; % so that display of NIFTI will not alter the data
% In SPM results, this is 1/255, then all uint8 data will be displayed
% in the range of [0 1] % ANDY 2018-06-04
white.fileprefix = [dirname filesep baseFilename '_masks'];
white.hdr.hist.descrip = 'tissue masks';
save_untouch_nii(white,[dirname filesep baseFilename '_masks.nii']);