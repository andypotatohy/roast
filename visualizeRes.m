function visualizeRes(subj,subjRasRSPD,segOut,mri2mni,T2,node,elem,face,inCurrent,imgHdr,uniTag,showAll,varargin)
% visualizeRes(subj,subjRasRSPD,segOut,T2,node,elem,face,inCurrent,imgHdr,uniTag,showAll,varargin)
%
% Display the simulation results. The 3D rendering is displayed in the
% world space, while the slice view is done in the voxel space.
%
% (c) Yu (Andy) Huang, Parra Lab at CCNY
% yhuang16@citymail.cuny.edu
% April 2018
% August 2019 callable by roast_target()
%
% (c) May 2024, sliceshow improved by Gavin Hsu and Andrew Birnbaum

if ndims(varargin{1})==3
    isRoast = 1;
    vol_all = varargin{1}; ef_mag = varargin{2}; ef_all = varargin{3};
else
    isRoast = 0;
    C = varargin{1}; ef_mag = varargin{2}; ef_all = varargin{3}; targetCoord = varargin{4};
end

[dirname,subjName] = fileparts(subj);
if isempty(dirname), dirname = pwd; end

if showAll    
    if ~strcmp(subjName,'nyhead')
        disp('showing MRI and segmentations...');
        data = load_untouch_nii(subjRasRSPD); sliceshow(data.img,[],'gray',[],[],'MRI: Click anywhere to navigate.',[],mri2mni); drawnow
        
        if ~isempty(T2) %T2 specified
            data = load_untouch_nii(T2);
            sliceshow(data.img,[],'gray',[],[],'MRI: T2. Click anywhere to navigate.',[],mri2mni); drawnow
        end
    else
        disp('NEW YORK HEAD selected, there is NO MRI for it to show.')
    end    
end

[~,segOutName] = fileparts(segOut);
masks = load_untouch_nii([dirname filesep segOutName '_masks.nii']);
allMask = masks.img;
numOfTissue = 6; % hard coded across ROAST.  max(allMask(:));
if isRoast
    gel = load_untouch_nii([dirname filesep subjName '_' uniTag '_mask_gel.nii']);
    numOfGel = max(gel.img(:));
    elec = load_untouch_nii([dirname filesep subjName '_' uniTag '_mask_elec.nii']);
    % numOfElec = max(elec.img(:));
else
    numOfGel = length(inCurrent);
    indMonElec = find(abs(inCurrent)>1e-3); % this is not perfect
end

if showAll
%     if isRoast,
%     viewElectrodes(allMask,elec.img>0,gel.img>0,landmarks,imgHdr,uniTag);
%     end % NEED TO THINK ABOUT HOW TO DO THIS MORE ELEGANTLY

    % More anatomical looking colormap
    color_map = [
        0, 0, 0;                   % background: black
        1, 1, 1;                   % white matter: white
        0.7, 0.7, 0.7;             % gray matter: gray
        105/255, 175/255, 255/255; % CSF: blue
        241/255, 214/255, 145/255; % bone
        177/255, 122/255, 101/255; % skin
        0.6863, 0.8824, 0.6863;    % air cavities
        0, 0, 139/255;             % gel&elec
        ];

    allMaskShow = masks.img;
    allMaskShow(gel.img>0) = numOfTissue + 1;
    allMaskShow(elec.img>0) = numOfTissue + 2;
    sliceshow(allMaskShow,[],color_map,[],'Tissue index','Segmentation. Click anywhere to navigate.',[],mri2mni)
    drawnow
end

% node = node + 0.5; already done right after mesh

% scrsz = get(groot,'ScreenSize');

disp('generating 3D renderings...')

for i=1:3, node(:,i) = node(:,i)/imgHdr.mat(i,i); end
% convert pseudo-world coordinates back to voxel coordinates so that the
% following conversion to pure-world space is meaningful

voxCoord = [node(:,1:3) ones(size(node,1),1)];
worldCoord = (imgHdr.mat*voxCoord')';
% do the 3D rendering in world space, to avoid confusion in left-right;
% sliceshow below is still in voxel space though
node(:,1:3) = worldCoord(:,1:3);

indNode_grayFace = face(find(face(:,4) == 2),1:3);
indNode_grayElm = elem(find(elem(:,5) == 2),1:4);

% node(:,1:3) = sms(node(:,1:3),indNode_grayFace);
% % smooth the surface that's to be displayed
% % just for display, the output data is not smoothed
% % very slow if mesh is big

if isRoast
    
    indNode_elecFace = face(find(face(:,4) > numOfTissue+numOfGel),1:3);
    indNode_elecElm = elem(find(elem(:,5) > numOfTissue+numOfGel),1:4);
    inCurrentRange = [min(inCurrent) max(inCurrent)];
    
    fid = fopen([dirname filesep subjName '_' uniTag '_v.pos']);
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
        %     indNodeTemp = indNode_elecElm(find(label_elec==i),:);
        %     dataShow(unique(indNodeTemp(:)),4) = dataShowForElec(i);
        dataShow(unique(elem(find(elem(:,5) == numOfTissue+numOfGel+i),1:4)),4) = dataShowForElec(i);
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
    
    fid = fopen([dirname filesep subjName '_' uniTag '_e.pos']);
    fgetl(fid);
    C = textscan(fid,'%d %f %f %f');
    fclose(fid);
    
    % dataShow = [node(C{1},1:3), C_ef_mag];
    color = nan(size(node,1),1);
    color(C{1}) = sqrt(C{2}.^2+C{3}.^2+C{4}.^2);
    dataShow = [node(:,1:3) color];
    
    figName = ['Electric field in Simulation: ' uniTag];
    figure('Name',[figName '. Move your mouse to rotate.'],'NumberTitle','off');
    set(gcf,'color','w');
    colormap(jet);
    plotmesh(dataShow,indNode_grayFace,indNode_grayElm,'LineStyle','none');
%     dataShowVal = dataShow(unique(indNode_grayElm(:)),4);
    dataShowRange = [min(dataShow(unique(indNode_grayElm(:)),4)) prctile(dataShow(unique(indNode_grayElm(:)),4),95)];
    dataShowForElec = interp1(inCurrentRange,dataShowRange,inCurrent);
    for i=1:length(inCurrent)
        %     indNodeTemp = indNode_elecElm(find(label_elec==i),:);
        %     dataShow(unique(indNodeTemp(:)),4) = dataShowForElec(i);
        dataShow(unique(elem(find(elem(:,5) == numOfTissue+numOfGel+i),1:4)),4) = dataShowForElec(i);
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
    
else
    
    indNode_elecFace = face(ismember(face(:,4),numOfTissue+numOfGel+indMonElec),1:3);
    indNode_elecElm = elem(ismember(elem(:,5),numOfTissue+numOfGel+indMonElec),1:4);
    inCurrentRange = [min(inCurrent) max(inCurrent)];
    
    color = nan(size(node,1),1);
    color(C(:,1)) = sqrt(sum(C(:,2:4).^2,2));
    dataShow = [node(:,1:3) color];
    
    figName = ['Electric field in Targeting: ' uniTag];
    figure('Name',[figName '. Move your mouse to rotate.'],'NumberTitle','off');
    set(gcf,'color','w');
    colormap(jet);
    plotmesh(dataShow,indNode_grayFace,indNode_grayElm,'LineStyle','none');
%     dataShowVal = dataShow(unique(indNode_grayElm(:)),4);
    dataShowRange = [min(dataShow(unique(indNode_grayElm(:)),4)) prctile(dataShow(unique(indNode_grayElm(:)),4),95)];
    dataShowForElec = interp1(inCurrentRange,dataShowRange,inCurrent(indMonElec));
    for i=1:length(indMonElec)
        %     indNodeTemp = indNode_elecElm(find(label_elec==i),:);
        %     dataShow(unique(indNodeTemp(:)),4) = dataShowForElec(i);
        dataShow(unique(elem(find(elem(:,5) == numOfTissue+numOfGel+indMonElec(i)),1:4)),4) = dataShowForElec(i);
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
    
end

disp('generating slice views...');

brain = (allMask==1 | allMask==2);
nan_mask_brain = nan(size(brain));
nan_mask_brain(find(brain)) = 1;

cm = colormap(jet(2^11)); cm = [1 1 1;cm];
bbox = brainCrop(segOut);

if isRoast
    figName = ['Voltage in Simulation: ' uniTag];
    sliceshow(vol_all.*nan_mask_brain,round(mean(bbox)),cm,[],'Voltage (mV)',[figName '. Click anywhere to navigate.'],[],mri2mni,bbox); drawnow
end

for i=1:size(ef_all,4), ef_all(:,:,:,i) = ef_all(:,:,:,i).*nan_mask_brain; end
ef_mag = ef_mag.*nan_mask_brain;
dataShowVal = ef_mag(~isnan(ef_mag(:)));
if isRoast
    figName = ['Electric field in Simulation: ' uniTag];
    sliceshow(ef_mag,round(mean(bbox)),cm,[min(dataShowVal) prctile(dataShowVal,95)],'Electric field (V/m)',[figName '. Click anywhere to navigate.'],ef_all,mri2mni,bbox); drawnow
else
    
    for i=1:size(targetCoord,1)
        figName = ['Electric field at Target ' num2str(i) ' in Targeting: ' uniTag];
        sliceshow(ef_mag,targetCoord(i,:),cm,[min(dataShowVal) prctile(dataShowVal,95)],'Electric field (V/m)',[figName '. Click anywhere to navigate.'],ef_all,mri2mni,bbox); drawnow
    end
end