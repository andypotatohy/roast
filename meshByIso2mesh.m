function [node,elem,face] = meshByIso2mesh(subj,spmOut,segOut,opt,hdrInfo,uniTag)
% [node,elem,face] = meshByIso2mesh(subj,spmOut,segOut,opt,hdrInfo,uniTag)
%
% Generate volumetric tetrahedral mesh using iso2mesh toolbox
% http://iso2mesh.sourceforge.net/cgi-bin/index.cgi?Download
%
% (c) Yu (Andy) Huang, Parra Lab at CCNY
% yhuang16@citymail.cuny.edu
% October 2017

[dirname,subjName] = fileparts(subj);
if isempty(dirname), dirname = pwd; end

[~,spmOutName] = fileparts(spmOut);
mappingFile = [dirname filesep spmOutName '_seg8.mat'];
if ~exist(mappingFile,'file')
    error(['Mapping file ' mappingFile ' from SPM not found. Please check if you run through SPM segmentation in ROAST.']);
else
    load(mappingFile,'image','Affine');
    mri2mni = Affine*image(1).mat;
    % mapping from MRI voxel space to MNI space
end

[~,segOutName] = fileparts(segOut);
data = load_untouch_nii([dirname filesep segOutName '_masks.nii']);
allMask = data.img;
allMaskShow = data.img;
numOfTissue = 6; % hard coded across ROAST.  max(allMask(:));
% data = load_untouch_nii([dirname filesep subjName '_' uniTag '_mask_gel.nii']);
% allMask(data.img==255) = 7;
% data = load_untouch_nii([dirname filesep subjName '_' uniTag '_mask_elec.nii']);
% allMask(data.img==255) = 8;

data = load_untouch_nii([dirname filesep subjName '_' uniTag '_mask_gel.nii']);
numOfGel = max(data.img(:));
for i=1:numOfGel
    allMask(data.img==i) = numOfTissue + i;
end
allMaskShow(data.img>0) = numOfTissue + 1;
data = load_untouch_nii([dirname filesep subjName '_' uniTag '_mask_elec.nii']);
numOfElec = max(data.img(:));
for i=1:numOfElec
    allMask(data.img==i) = numOfTissue + numOfGel + i;
end
allMaskShow(data.img>0) = numOfTissue + 2;

% sliceshow(allMask,[],[],[],'Tissue index','Segmentation. Click anywhere to navigate.')
sliceshow(allMaskShow,[],[],[],'Tissue index','Segmentation. Click anywhere to navigate.',[],mri2mni,[])
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

for i=1:3, node(:,i) = node(:,i)*hdrInfo.pixdim(i); end
% Put mesh coordinates into pseudo-world space (voxel space but scaled properly
% using the scaling factors in the header) to avoid mistakes in
% solving. Putting coordinates into pure-world coordinates causes other
% complications. Units of coordinates are mm here. No need to convert into
% meter as voltage output from solver is mV.
% ANDY 2019-03-13

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
savemsh(node(:,1:3),elem,[dirname filesep subjName '_' uniTag '.msh'],maskName);
save([dirname filesep subjName '_' uniTag '.mat'],'node','elem','face');