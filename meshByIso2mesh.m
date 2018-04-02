function [node,elem,face] = meshByIso2mesh(P,opt,uniTag)
% [node,elem,face] = meshByIso2mesh(P,opt,uniTag)
%
% generate volumetric tetrahedral mesh using iso2mesh toolbox
% http://iso2mesh.sourceforge.net/cgi-bin/index.cgi?Download
%
% (c) Yu (Andy) Huang, Parra Lab at CCNY
% October 2017

[dirname,baseFilename] = fileparts(P);
if isempty(dirname), dirname = pwd; end

maskName = {'white','gray','csf','bone','skin','air','gel','elec'};

for i=1:length(maskName)
    
    if strcmp(maskName{i},'gel') || strcmp(maskName{i},'elec') 
        data = load_untouch_nii([dirname filesep baseFilename '_' uniTag '_mask_' maskName{i} '.nii']);
    else
        data = load_untouch_nii([dirname filesep baseFilename '_mask_' maskName{i} '.nii']);
    end
    img = data.img;
    
    if i==1, [dim1,dim2,dim3] = size(img); allMask = zeros(dim1,dim2,dim3); end
    
    allMask(img==255) = i;
end

sliceshow(allMask,[],[],[],'Tissue index','Segmentation. Click anywhere to navigate.')
drawnow

allMask = uint8(allMask);

% opt.radbound = 5; % default 6, maximum surface element size
% opt.angbound = 30; % default 30, miminum angle of a surface triangle
% opt.distbound = 0.4; % default 0.5, maximum distance
% % between the center of the surface bounding circle and center of the element bounding sphere
% opt.reratio = 3; % default 3, maximum radius-edge ratio
% maxvol = 10; %100; % target maximum tetrahedral elem volume

[node,elem,face] = cgalv2m(allMask,opt,maxvol);
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
maskName = {'WHITE','GRAY','CSF','BONE','SKIN','AIR','GEL','ELEC'};
savemsh(node(:,1:3),elem,[dirname filesep baseFilename '_' uniTag '.msh'],maskName);
save([dirname filesep baseFilename '_' uniTag '.mat'],'node','elem','face');