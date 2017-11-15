function visualizeRes(P,node,elem,face,vol_all,ef_mag,inCurrent,uniTag,showAll)
% visualizeRes(P,node,elem,face,vol_all,ef_mag,inCurrent,uniTag,showAll)
%
% Display the simulation results (in voxel space)
%
% (c) Yu (Andy) Huang, Parra Lab at CCNY
% October 2017

[dirname,baseFilename] = fileparts(P);
if isempty(dirname), dirname = pwd; end

if showAll
    
    disp('showing MRI and segmentations...');
    
    ref = deblank(P);
    data = load_untouch_nii(ref); sliceshow(data.img,[],'gray',[],[],'MRI: T1. Click anywhere to navigate.'); drawnow
    
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
    
end

% node = node + 0.5; already done right after mesh

% scrsz = get(groot,'ScreenSize');

disp('generating 3D renderings...')
indFace = face(find(face(:,4) == 2),1:3);
indElem = elem(find(elem(:,5) == 2),1:4);

totInCurMag = sum(abs(inCurrent))/2;

fid = fopen([dirname filesep baseFilename '_' uniTag '_v.pos']);
fgetl(fid);
C = textscan(fid,'%d %f');
fclose(fid);

C{2} = C{2} - min(C{2}); % re-reference the voltage

% dataShow = [node(C{1},1:3) C{2}];
color = nan(size(node,1),1);
color(C{1}) = C{2};
dataShow = [node(:,1:3) color];

figName = ['Voltage in Simulation: ' uniTag];
figure('Name',[figName '. Move your mouse to rotate.'],'NumberTitle','off');
set(gcf,'color','w');
colormap(jet);
plotmesh(dataShow,indFace,indElem,'LineStyle','none');
hold on; plotmesh(dataShow,face(find(face(:,4) == 8),1:3),elem(find(elem(:,5) == 8),1:4),'LineStyle','none');
axis off; rotate3d on;
% set(hp1,'SpecularColorReflectance',0,'SpecularExponent',50);
caxis([min(dataShow(unique(indElem(:)),4)) max(dataShow(unique(indElem(:)),4))]);
lightangle(-90,45)
lightangle(90,45)
lightangle(-90,-45)
hc1 = colorbar; set(hc1,'FontSize',18,'YAxisLocation','right');
title(hc1,'Voltage (mV)','FontSize',18);
drawnow

fid = fopen([dirname filesep baseFilename '_' uniTag '_e.pos']);
fgetl(fid);
C = textscan(fid,'%d %f %f %f');
fclose(fid);

C_ef_mag = sqrt(C{2}.^2+C{3}.^2+C{4}.^2);
% dataShow = [node(C{1},1:3), C_ef_mag];
color = nan(size(node,1),1);
color(C{1}) = C_ef_mag;
dataShow = [node(:,1:3) color];

figName = ['Electric field in Simulation: ' uniTag];
figure('Name',[figName '. Move your mouse to rotate.'],'NumberTitle','off');
set(gcf,'color','w');
colormap(jet);
plotmesh(dataShow,indFace,indElem,'LineStyle','none');
hold on; plotmesh(dataShow,face(find(face(:,4) == 8),1:3),elem(find(elem(:,5) == 8),1:4),'LineStyle','none');
axis off; rotate3d on;
% set(hp2,'SpecularColorReflectance',0,'SpecularExponent',50);
caxis([0 0.3*totInCurMag]);
lightangle(-90,45)
lightangle(90,45)
lightangle(-90,-45)
hc2 = colorbar; set(hc2,'FontSize',18,'YAxisLocation','right');
title(hc2,'Electric field (V/m)','FontSize',18);
drawnow

disp('generating slice views...');

data = load_untouch_nii([dirname filesep baseFilename '_mask_gray.nii']); gray = data.img;
data = load_untouch_nii([dirname filesep baseFilename '_mask_white.nii']); white = data.img;
brain = gray | white;
nan_mask_brain = nan(size(brain));
nan_mask_brain(find(brain)) = 1;

figName = ['Voltage in Simulation: ' uniTag];
sliceshow(vol_all.*nan_mask_brain,[],[],[],'Voltage (mV)',[figName '. Click anywhere to navigate.']); drawnow

figName = ['Electric field in Simulation: ' uniTag];
sliceshow(ef_mag.*nan_mask_brain,[],[],[0 0.3*totInCurMag],'Electric field (V/m)',[figName '. Click anywhere to navigate.']); drawnow
