function autoPatching(P,T2)
% autoPatching(P,T2)
%
% Auto-patching to fill in holes that cannot be fully handled by mysegment()
%
% (c) Yu (Andy) Huang, Parra Lab at CCNY
% yhuang16@citymail.cuny.edu
% (c) Chris Thomas, Soterix Medical Inc.
% June 2017

[dirname,baseFilename] = fileparts(P);
if isempty(dirname), dirname = pwd; end

if isempty(T2)
    baseFilename = [baseFilename '_T1orT2'];
else
    baseFilename = [baseFilename '_T1andT2'];
end

load([dirname filesep baseFilename '_rmask.mat'],'holes_vol','eyes_vol','WMexclude_vol');
maskName = {'gray','white','csf','bone','skin','air'};

for i=1:length(maskName)
    
    data = load_untouch_nii([dirname filesep baseFilename '_mask_' maskName{i} '.nii']);
    img = data.img;
    
    if i==1
        [dim1,dim2,dim3] = size(img);
        allMask = zeros(dim1,dim2,dim3);
    end
    allMask(img==255) = i;
    
end

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
% ind_maskOfInterest = find(allMask==2);
% handle local specifics for WM here (there is no gray matter in brain stem and around ventricles)
disp('patching gray matter...');
WMall = allMask==2;
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
        allMask(ind_maskOfInterest(find(allMask(ind_maskOfIntestNeighbor)==wrongNeighborInd(t)))) = 1;
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
for i=1:length(maskName)
    data.img = uint8(allMask == i)*255;
    save_untouch_nii(data,[dirname filesep baseFilename '_mask_' maskName{i} '.nii']);
end