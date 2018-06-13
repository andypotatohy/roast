function [node,elem,face] = meshByIso2mesh(P1,P2,T2,opt,uniTag)
% [node,elem,face] = meshByIso2mesh(P1,P2,T2,opt,uniTag)
%
% Generate volumetric tetrahedral mesh using iso2mesh toolbox
% http://iso2mesh.sourceforge.net/cgi-bin/index.cgi?Download
%
% (c) Yu (Andy) Huang, Parra Lab at CCNY
% yhuang16@citymail.cuny.edu
% October 2017

[dirname,baseFilename] = fileparts(P1);
if isempty(dirname), dirname = pwd; end
[~,baseFilenameRSPD] = fileparts(P2);
if isempty(T2)
    baseFilenameRSPD = [baseFilenameRSPD '_T1orT2'];
else
    baseFilenameRSPD = [baseFilenameRSPD '_T1andT2'];
end

data = load_untouch_nii([dirname filesep baseFilenameRSPD '_masks.nii']);
allMask = data.img;
allMaskShow = data.img;
numOfTissue = 6; % hard coded across ROAST.  max(allMask(:));
% data = load_untouch_nii([dirname filesep baseFilename '_' uniTag '_mask_gel.nii']);
% allMask(data.img==255) = 7;
% data = load_untouch_nii([dirname filesep baseFilename '_' uniTag '_mask_elec.nii']);
% allMask(data.img==255) = 8;

data = load_untouch_nii([dirname filesep baseFilename '_' uniTag '_mask_gel.nii']);
numOfGel = max(data.img(:));
for i=1:numOfGel
    allMask(data.img==i) = numOfTissue + i;
end
allMaskShow(data.img>0) = numOfTissue + 1;
data = load_untouch_nii([dirname filesep baseFilename '_' uniTag '_mask_elec.nii']);
numOfElec = max(data.img(:));
for i=1:numOfElec
    allMask(data.img==i) = numOfTissue + numOfGel + i;
end
allMaskShow(data.img>0) = numOfTissue + 2;

% sliceshow(allMask,[],[],[],'Tissue index','Segmentation. Click anywhere to navigate.')
sliceshow(allMaskShow,[],[],[],'Tissue index','Segmentation. Click anywhere to navigate.')
drawnow

% allMask = uint8(allMask);

% opt.radbound = 5; % default 6, maximum surface element size
% opt.angbound = 30; % default 30, miminum angle of a surface triangle
% opt.distbound = 0.4; % default 0.5, maximum distance
% % between the center of the surface bounding circle and center of the element bounding sphere
% opt.reratio = 3; % default 3, maximum radius-edge ratio
% maxvol = 10; %100; % target maximum tetrahedral elem volume

[node,elem,face] = cgalv2m(allMask,opt,opt.maxvol);
node(:,1:3) = node(:,1:3) + 0.5; % then voxel space

% figure;
% % plotmesh(node(:,1:3),face,elem)
%
% % visualize tissue by tissue
% for i=1:length(maskName)
%     indElem = find(elem(:,5) == i);
%     indFace = find(face(:,4) == i);
%     plotmesh(node(:,1:3),face(indFace,:),elem(indElem,:))
%     title(maskName{i})
%     pause
% end

disp('saving mesh...')
% maskName = {'WHITE','GRAY','CSF','BONE','SKIN','AIR','GEL','ELEC'};
maskName = cell(1,numOfTissue+numOfGel+numOfElec);
maskName(1:numOfTissue) = {'WHITE','GRAY','CSF','BONE','SKIN','AIR'};
for i=1:numOfGel, maskName{numOfTissue+i} = ['GEL' num2str(i)]; end
for i=1:numOfElec, maskName{numOfTissue+numOfGel+i} = ['ELEC' num2str(i)]; end
savemsh(node(:,1:3),elem,[dirname filesep baseFilename '_' uniTag '.msh'],maskName);
save([dirname filesep baseFilename '_' uniTag '.mat'],'node','elem','face');