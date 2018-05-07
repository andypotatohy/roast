function visualizeRes(P1,P2,T2,node,elem,face,vol_all,ef_mag,inCurrent,label_elec,hdrInfo,uniTag,showAll)
% visualizeRes(P1,P2,T2,node,elem,face,vol_all,ef_mag,inCurrent,label_elec,hdrInfo,uniTag,showAll)
%
% Display the simulation results. The 3D rendering is displayed in the
% world space, while the slice view is done in the voxel space.
%
% (c) Yu (Andy) Huang, Parra Lab at CCNY
% yhuang16@citymail.cuny.edu
% April 2018

[dirname,baseFilename] = fileparts(P1);
if isempty(dirname), dirname = pwd; end
[~,baseFilenameRSPD] = fileparts(P2);
if isempty(T2)
    baseFilenameRSPD = [baseFilenameRSPD '_T1orT2'];
else
    baseFilenameRSPD = [baseFilenameRSPD '_T1andT2'];
end

if showAll
    
    if ~strcmp(baseFilename,'nyhead')
        disp('showing MRI and segmentations...');
        data = load_untouch_nii(P2); sliceshow(data.img,[],'gray',[],[],'MRI: Click anywhere to navigate.'); drawnow
        
        if ~isempty(T2) %T2 specified
            data = load_untouch_nii(T2);
            sliceshow(data.img,[],'gray',[],[],'MRI: T2. Click anywhere to navigate.'); drawnow
        end
    else
        disp('NEW YORK HEAD selected, there is NO MRI for it to show.')
    end
    
    maskName = {'white','gray','csf','bone','skin','air','gel','elec'};
    
    for i=1:length(maskName)
        
        if strcmp(maskName{i},'gel') || strcmp(maskName{i},'elec')
            data = load_untouch_nii([dirname filesep baseFilename '_' uniTag '_mask_' maskName{i} '.nii']);
        else
            data = load_untouch_nii([dirname filesep baseFilenameRSPD '_mask_' maskName{i} '.nii']);
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
voxCoord = [node(:,1:3) ones(size(node,1),1)];
worldCoord = (hdrInfo.v2w*voxCoord')';
% do the 3D rendering in world space, to avoid confusion in left-right;
% sliceshow below is still in voxel space though
node(:,1:3) = worldCoord(:,1:3);

indNode_grayFace = face(find(face(:,4) == 2),1:3);
indNode_grayElm = elem(find(elem(:,5) == 2),1:4);
indNode_elecFace = face(find(face(:,4) == 8),1:3);
indNode_elecElm = elem(find(elem(:,5) == 8),1:4);

% node(:,1:3) = sms(node(:,1:3),indNode_grayFace);
% % smooth the surface that's to be displayed
% % just for display, the output data is not smoothed
% % very slow if mesh is big

% totInCurMag = sum(abs(inCurrent))/2;
totInCurMag = max(abs(inCurrent)); % maybe a user option?
inCurrentRange = [min(inCurrent) max(inCurrent)];

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
plotmesh(dataShow,indNode_grayFace,indNode_grayElm,'LineStyle','none');
dataShowRange = [min(dataShow(unique(indNode_grayElm(:)),4)) max(dataShow(unique(indNode_grayElm(:)),4))];
dataShowForElec = interp1(inCurrentRange,dataShowRange,inCurrent);
for i=1:length(inCurrent)
    indNodeTemp = indNode_elecElm(find(label_elec==i),:);
    dataShow(unique(indNodeTemp(:)),4) = dataShowForElec(i);
end % to show injected current intensities properly
hold on;
plotmesh(dataShow,indNode_elecFace,indNode_elecElm,'LineStyle','none');
axis off; rotate3d on;
% set(hp1,'SpecularColorReflectance',0,'SpecularExponent',50);
caxis(dataShowRange);
lightangle(-90,45)
lightangle(90,45)
lightangle(-90,-45)
hc1 = colorbar; set(hc1,'FontSize',18,'YAxisLocation','right');
title(hc1,'Voltage (mV)','FontSize',18);
a1 = gca;
a2 = axes('Color','none','Position',get(a1,'Position'),'XLim',get(a1,'XLim'),'YLim',get(a1,'YLim'),'ZLim',get(a1,'ZLim'));
axis off;
hc2 = colorbar; set(hc2,'FontSize',18,'YAxisLocation','right','Location','westoutside');
title(hc2,'Injected current (mA)','FontSize',18);
caxis(inCurrentRange);
axes(a1);
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
plotmesh(dataShow,indNode_grayFace,indNode_grayElm,'LineStyle','none');
dataShowRange = [0 0.3*totInCurMag];
dataShowForElec = interp1(inCurrentRange,dataShowRange,inCurrent);
for i=1:length(inCurrent)
    indNodeTemp = indNode_elecElm(find(label_elec==i),:);
    dataShow(unique(indNodeTemp(:)),4) = dataShowForElec(i);
end % to show injected current intensities properly
hold on;
plotmesh(dataShow,indNode_elecFace,indNode_elecElm,'LineStyle','none');
axis off; rotate3d on;
% set(hp2,'SpecularColorReflectance',0,'SpecularExponent',50);
caxis(dataShowRange);
lightangle(-90,45)
lightangle(90,45)
lightangle(-90,-45)
hc1 = colorbar; set(hc1,'FontSize',18,'YAxisLocation','right');
title(hc1,'Electric field (V/m)','FontSize',18);
a1 = gca;
a2 = axes('Color','none','Position',get(a1,'Position'),'XLim',get(a1,'XLim'),'YLim',get(a1,'YLim'),'ZLim',get(a1,'ZLim'));
axis off;
hc2 = colorbar; set(hc2,'FontSize',18,'YAxisLocation','right','Location','westoutside');
title(hc2,'Injected current (mA)','FontSize',18);
caxis(inCurrentRange);
axes(a1);
drawnow

disp('generating slice views...');

data = load_untouch_nii([dirname filesep baseFilenameRSPD '_mask_gray.nii']); gray = data.img;
data = load_untouch_nii([dirname filesep baseFilenameRSPD '_mask_white.nii']); white = data.img;
brain = gray | white;
nan_mask_brain = nan(size(brain));
nan_mask_brain(find(brain)) = 1;

cm = colormap(jet(512)); cm = [1 1 1;cm];
figName = ['Voltage in Simulation: ' uniTag];
sliceshow(vol_all.*nan_mask_brain,[],cm,[],'Voltage (mV)',[figName '. Click anywhere to navigate.']); drawnow

figName = ['Electric field in Simulation: ' uniTag];
sliceshow(ef_mag.*nan_mask_brain,[],cm,[0 0.3*totInCurMag],'Electric field (V/m)',[figName '. Click anywhere to navigate.']); drawnow

disp('=======================ALL DONE=======================');