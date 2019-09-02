function reviewRes(subj,simTag,tissue,fastRender,tarTag)
% reviewRes(subj,simTag,tissue,fastRender,tarTag)
%
% A simpler interface to visualize the simulations/targetings that are already done.
% Users do not have to enter all the parameters as they would have in
% roast() or roast_target() main function. Instead, they just need to enter the path
% to the input MRI and the unique simulation/targeting tag for that MRI.
%
% Please refer to the README.md on the github repo for better formated
% documentations: https://github.com/andypotatohy/roast
% 
% If you use ROAST in your research, please cite these:
% 
% Huang, Y., Datta, A., Bikson, M., Parra, L.C., Realistic vOlumetric-Approach
% to Simulate Transcranial Electric Stimulation -- ROAST -- a fully automated
% open-source pipeline, Journal of Neural Engineering, Vol. 16, No. 5, 2019 (prefered reference)
% 
% Huang, Y., Datta, A., Bikson, M., Parra, L.C., ROAST: an open-source,
% fully-automated, Realistic vOlumetric-Approach-based Simulator for TES,
% Proceedings of the 40th Annual International Conference of the IEEE Engineering
% in Medicine and Biology Society, Honolulu, HI, July 2018
% 
% If you use New York head to run simulation, please also cite the following:
% Huang, Y., Parra, L.C., Haufe, S.,2016. The New York Head - A precise
% standardized volume conductor model for EEG source localization and tES
% targeting, NeuroImage,140, 150-162
% 
% If you also use the targeting feature (`roast_target`), please cite these:
% 
% Dmochowski, J.P., Datta, A., Bikson, M., Su, Y., Parra, L.C., Optimized 
% multi-electrode stimulation increases focality and intensity at target,
% Journal of Neural Engineering 8 (4), 046011, 2011
% 
% Dmochowski, J.P., Datta, A., Huang, Y., Richardson, J.D., Bikson, M.,
% Fridriksson, J., Parra, L.C., Targeted transcranial direct current stimulation 
% for rehabilitation after stroke, NeuroImage, 75, 12-19, 2013
% 
% Huang, Y., Thomas, C., Datta, A., Parra, L.C., Optimized tDCS for Targeting
% Multiple Brain Regions: An Integrated Implementation. Proceedings of the 40th
% Annual International Conference of the IEEE Engineering in Medicine and Biology
% Society, Honolulu, HI, July 2018, 3545-3548
% 
% ROAST was supported by NIH through grants R01MH111896, R01MH111439, 
% R01NS095123, R44NS092144, R41NS076123, and by Soterix Medical Inc.
% 
% General Public License version 3 or later. See LICENSE.md for details.
% 
% This software uses free packages from the Internet, except Matlab, which
% is a proprietary software by the MathWorks. You need a valid Matlab license
% to run this software.
% 
% ROAST is considered as an "aggregate" rather than "derived work", based on
% the definitions in GPL FAQ. The ROAST license only applies to the scripts,
% documentation and the individual MRI data under example/ folder in this 
% package and excludes those programs stored in the lib/ directory. The software 
% under lib/ follow their respective licenses. This software is only intended
% for non-commercial use.
% 
% (c) Yu (Andy) Huang, Parra Lab at CCNY
% yhuang16@citymail.cuny.edu
% September 2019

fprintf('\n\n');
disp('=============================================================')
disp('ROAST is an aggregated work by Yu (Andy) Huang licensed under')
disp('General Public License version 3 or later. It''s supported by')
disp('both NIH grants and Soterix Medical Inc.')
disp('=============================================================')

addpath(genpath([fileparts(which(mfilename)) filesep 'lib/']));

% check subject name
if nargin<1 || isempty(subj)
    subj = 'example/MNI152_T1_1mm.nii';
end

% check simulation tag
if nargin<2 || isempty(simTag)
    error(['Please provide a valid simulation tag for Subject ' subj]);
end

if strcmpi(subj,'nyhead')
    subj = 'example/nyhead.nii';
end

% check tissue to be visualized
if nargin<3 || isempty(tissue)
    tissue = 'brain';
end

switch lower(tissue)
    case 'white'
        indSurfShow = 1;
        indSliceShow = 1;
    case 'gray'
        indSurfShow = 2;
        indSliceShow = 2;
    case 'csf'
        indSurfShow = 3;
        indSliceShow = 3;
    case 'bone'
        indSurfShow = 4;
        indSliceShow = 4;
    case 'skin'
        indSurfShow = 5;
        indSliceShow = 5;
    case 'air'
        indSurfShow = 6;
        indSliceShow = 6;
    case 'brain'
        indSurfShow = 2;
        indSliceShow = 1:2;
    case 'all'
        indSurfShow = 5;
        indSliceShow = 1:6;
    otherwise
        error('Supported tissues to be displayed are: ''white'', ''gray'', ''CSF'', ''bone'', ''skin'', ''air'', ''brain'' and ''all''.');
end

% check if do fast rendering
if nargin<4 || isempty(fastRender)
    fastRender = 1; % no smoothing on surface
end

[dirname,baseFilename,ext] = fileparts(subj);
optionFile = [dirname filesep baseFilename '_' simTag '_roastOptions.mat'];
if ~exist(optionFile,'file')
    error(['Option file not found. Simulation ' simTag ' may never be run. Please run it first.']);
else
    load(optionFile,'opt');
    optRoast = opt;
end

% check if visualizing results from roast() or roast_target()
if ~strcmp(optRoast.configTxt,'leadFieldGeneration')
    
    isRoast = 1;
    if exist('tarTag','var')
        warning(['Simulation ' simTag ' was not run for generating the lead field for subject ' subj ', so targeting tag ' tarTag ' will be ignored.']);
    end
    
    disp(['Showing results for Simulation ' simTag ' ...']);
    
else
    
    isRoast = 0;
    
    if nargin<5 || isempty(tarTag)
        error(['Simulation ' simTag ' was run for generating the lead field for subject ' subj '. reviewRes() will visualize the results from roast_target(), but no targeting tag was provided.']);
    end
    
    optionFile = [dirname filesep baseFilename '_' tarTag '_targetOptions.mat'];
    if ~exist(optionFile,'file')
        error(['Option file not found. Targeting ' tarTag ' may never be run. Please run it first.']);
    else
        load(optionFile,'opt');
        optTarget = opt;
    end
    
    disp(['Showing results for Targeting ' tarTag ' ...']);
    
    resFile = [dirname filesep baseFilename '_' tarTag '_targetResult.mat'];
    if ~exist(resFile,'file')
        error(['Result file ' resFile ' not found. Check if you run through targeting under tag ' tarTag '.']);
    else
        load(resFile,'r');
    end
    
end

% to locate related files (e.g. MRI header, *_seg8 mapping, tissue masks)
if optRoast.isNonRAS
    subjRas = [dirname filesep baseFilename '_ras' ext];
else
    subjRas = subj;
end

if optRoast.resamp
    [dirname2,baseFilename2,ext2] = fileparts(subjRas);
    subjRasRS = [dirname filesep baseFilename2 '_1mm' ext];
else
    subjRasRS = subjRas;
end

if optRoast.zeroPad>0
    [dirname2,baseFilename2,ext2] = fileparts(subjRasRS);
    subjRasRSPD = [dirname2 filesep baseFilename2 '_padded' num2str(optRoast.zeroPad) ext2];
    %     subjRasRSPD = ['example/nyhead_padded' num2str(paddingAmt) '.nii'];
else
    subjRasRSPD = subjRasRS;
end

[~,baseFilenameRasRSPD] = fileparts(subjRasRSPD);
if isempty(optRoast.T2)
    mappingFile = [dirname filesep baseFilenameRasRSPD '_T1orT2_seg8.mat'];
else
    mappingFile = [dirname filesep baseFilenameRasRSPD '_T1andT2_seg8.mat'];
end
if ~exist(mappingFile,'file')
    error(['Mapping file ' mappingFile ' from SPM not found. Please check if you run through SPM segmentation in ROAST.']);
else
    load(mappingFile,'image','Affine');
    mri2mni = Affine*image(1).mat;
    % mapping from MRI voxel space to MNI space
end

if isRoast
    
    lp = strfind(optRoast.configTxt,'(');
    rp = strfind(optRoast.configTxt,')');
    
    inCurrent = zeros(length(lp),1);
    for i=1:length(lp)
        inCurrent(i) = str2num(optRoast.configTxt(lp(i)+1:rp(i)-4));
    end
    
    if ~strcmp(baseFilename,'nyhead')
        
        disp('showing MRI and segmentations...');
        if ~exist(subjRasRSPD,'file')
            error(['The subject MRI you provided ' subjRasRSPD ' does not exist. Check if you run through resampling or zero-padding if you tried to do that.']);
        else
            data = load_untouch_nii(subjRasRSPD); sliceshow(data.img,[],'gray',[],[],'MRI: Click anywhere to navigate.',[],mri2mni); drawnow
        end
        
        if ~isempty(optRoast.T2) %T2 specified
            if ~exist(optRoast.T2,'file')
                error(['T2 file ' optRoast.T2 ' does not exist. You used that to run ROAST but maybe later deleted it.']);
            else
                data = load_untouch_nii(optRoast.T2);
                sliceshow(data.img,[],'gray',[],[],'MRI: T2. Click anywhere to navigate.',[],mri2mni); drawnow
            end
        end
    else
        disp('NEW YORK HEAD selected, there is NO MRI for it to show.')
    end
    
else
    
    fid = fopen('./elec72.loc'); C = textscan(fid,'%d %f %f %s'); fclose(fid);
    elecName = C{4}; for i=1:length(elecName), elecName{i} = strrep(elecName{i},'.',''); end
    elecPara = struct('capType','1010');
    
    [~,indInUsrInput] = elecPreproc(subj,elecName,elecPara);
    inCurrent = r.mon(indInUsrInput);
    
end

if isempty(optRoast.T2)
    masksFile = [dirname filesep baseFilenameRasRSPD '_T1orT2_masks.nii'];
else
    masksFile = [dirname filesep baseFilenameRasRSPD '_T1andT2_masks.nii'];
end
if ~exist(masksFile,'file')
    error(['Segmentation masks ' masksFile ' not found. Check if you run through MRI segmentation.']);
else
    masks = load_untouch_nii(masksFile);
end

numOfTissue = 6; % hard coded across ROAST.  max(allMask(:));

if isRoast
    
    gelMask = [dirname filesep baseFilename '_' simTag '_mask_gel.nii'];
    if ~exist(gelMask,'file')
        error(['Gel mask ' gelMask ' not found. Check if you run through electrode placement.']);
    else
        gel = load_untouch_nii(gelMask);
        numOfGel = max(gel.img(:));
    end
    elecMask = [dirname filesep baseFilename '_' simTag '_mask_elec.nii'];
    if ~exist(elecMask,'file')
        error(['Electrode mask ' elecMask ' not found. Check if you run through electrode placement.']);
    else
        elec = load_untouch_nii(elecMask);
        % numOfElec = max(elec.img(:));
    end
    
    allMaskShow = masks.img;
    allMaskShow(gel.img>0) = numOfTissue + 1;
    allMaskShow(elec.img>0) = numOfTissue + 2;
    sliceshow(allMaskShow,[],[],[],'Tissue index','Segmentation. Click anywhere to navigate.',[],mri2mni)
    drawnow
    
else
    
    numOfGel = length(inCurrent);
    indMonElec = find(abs(inCurrent)>1e-3); % this is not perfect
    
    cm = colormap(jet(64));
    if strcmpi(optTarget.optType,'max-l1') || strcmpi(optTarget.optType,'max-l1per')
        cm(3:62,:) = ones(60,3);
    end
    figure('Name',['Montage in Targeting: ' tarTag],'NumberTitle','off');
    mytopoplot(r.mon,'./elec72.loc','numcontour',0,'plotrad',0.9,'shading','flat','gridscale',1000,'whitebk','off','colormap',cm);
    hc = colorbar; set(hc,'FontSize',18,'YAxisLocation','right');
    title(hc,'Injected current (mA)','FontSize',18);
    caxis([min(r.mon) max(r.mon)]);
    drawnow
    
    disp('Electrodes used are:')
    disp(strrep(r.montageTxt,', ',newline));
    
end

% node = node + 0.5; already done right after mesh

disp('generating 3D renderings...')

meshFile = [dirname filesep baseFilename '_' simTag '.mat'];
if ~exist(meshFile,'file')
    error(['Mesh file ' meshFile ' not found. Check if you run through meshing.']);
else
    load(meshFile,'node','elem','face');
end

indNode_showFace = face(find(face(:,4) == indSurfShow),1:3);
indNode_showElm = elem(find(elem(:,5) == indSurfShow),1:4);

if ~fastRender
    node(:,1:3) = sms(node(:,1:3),indNode_showFace);
    % smooth the surface that's to be displayed (just for display, the output data is not smoothed)
    % very slow if mesh is big
end

hdrFile = [dirname filesep baseFilenameRasRSPD '_header.mat'];
if ~exist(hdrFile,'file')
    error(['Header file ' hdrFile ' not found. Check if you run through electrode placement.']);
else
    load(hdrFile,'hdrInfo');
end

for i=1:3, node(:,i) = node(:,i)/hdrInfo.pixdim(i); end
% convert pseudo-world coordinates back to voxel coordinates so that the
% following conversion to pure-world space is meaningful
voxCoord = [node(:,1:3) ones(size(node,1),1)];
worldCoord = (hdrInfo.v2w*voxCoord')';
% do the 3D rendering in world space, to avoid confusion in left-right;
% sliceshow below is still in voxel space though
node(:,1:3) = worldCoord(:,1:3);

if isRoast
    
    indNode_elecFace = face(find(face(:,4) > numOfTissue+numOfGel),1:3);
    indNode_elecElm = elem(find(elem(:,5) > numOfTissue+numOfGel),1:4);
    
    inCurrentRange = [min(inCurrent) max(inCurrent)];
    
    volFile = [dirname filesep baseFilename '_' simTag '_v.pos'];
    if ~exist(volFile,'file')
        error(['Solution file ' volFile ' not found. Check if you run through solving.']);
    else
        fid = fopen(volFile);
        fgetl(fid);
        C = textscan(fid,'%d %f');
        fclose(fid);
    end
    
    C{2} = C{2} - min(C{2}); % re-reference the voltage
    
    % dataShow = [node(C{1},1:3) C{2}];
    color = nan(size(node,1),1);
    color(C{1}) = C{2};
    dataShow = [node(:,1:3) color];
    
    figName = ['Voltage in Simulation: ' simTag];
    figure('Name',[figName '. Move your mouse to rotate.'],'NumberTitle','off');
    set(gcf,'color','w');
    colormap(jet);
    plotmesh(dataShow,indNode_showFace,indNode_showElm,'LineStyle','none');
    dataShowRange = [min(dataShow(unique(indNode_showElm(:)),4)) max(dataShow(unique(indNode_showElm(:)),4))];
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
    
    efFile = [dirname filesep baseFilename '_' simTag '_e.pos'];
    if ~exist(efFile,'file')
        error(['Solution file ' efFile ' not found. Check if you run through solving.']);
    else
        fid = fopen(efFile);
        fgetl(fid);
        C = textscan(fid,'%d %f %f %f');
        fclose(fid);
    end
    
    C_ef_mag = sqrt(C{2}.^2+C{3}.^2+C{4}.^2);
    % dataShow = [node(C{1},1:3), C_ef_mag];
    color = nan(size(node,1),1);
    color(C{1}) = C_ef_mag;
    dataShow = [node(:,1:3) color];
    
    figName = ['Electric field in Simulation: ' simTag];
    figure('Name',[figName '. Move your mouse to rotate.'],'NumberTitle','off');
    set(gcf,'color','w');
    colormap(jet);
    plotmesh(dataShow,indNode_showFace,indNode_showElm,'LineStyle','none');
%     dataShowVal = dataShow(unique(indNode_showElm(:)),4);
    dataShowRange = [min(dataShow(unique(indNode_showElm(:)),4)) prctile(dataShow(unique(indNode_showElm(:)),4),95)];
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
    C = r.xopt;
    color(C(:,1)) = sqrt(sum(C(:,2:4).^2,2));
    dataShow = [node(:,1:3) color];
    
    figName = ['Electric field in Targeting: ' tarTag];
    figure('Name',[figName '. Move your mouse to rotate.'],'NumberTitle','off');
    set(gcf,'color','w');
    colormap(jet);
    plotmesh(dataShow,indNode_showFace,indNode_showElm,'LineStyle','none');
%     dataShowVal = dataShow(unique(indNode_showElm(:)),4);
    dataShowRange = [min(dataShow(unique(indNode_showElm(:)),4)) prctile(dataShow(unique(indNode_showElm(:)),4),95)];
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

allMask = masks.img;
mask = zeros(size(allMask));
for i=1:length(indSliceShow)
    mask = (mask | allMask==indSliceShow(i));
end
nan_mask = nan(size(mask));
nan_mask(find(mask)) = 1;

cm = colormap(jet(2^11)); cm = [1 1 1;cm];

if isRoast
    
    resFile = [dirname filesep baseFilename '_' simTag '_roastResult.mat'];
    if ~exist(resFile,'file')
        error(['Result file ' resFile ' not found. Check if you run through post processing after solving.']);
    else
        load(resFile,'vol_all','ef_mag','ef_all');
    end
    
    figName = ['Voltage in Simulation: ' simTag];
    sliceshow(vol_all.*nan_mask,[],cm,[],'Voltage (mV)',[figName '. Click anywhere to navigate.'],[],mri2mni); drawnow
    
    figName = ['Electric field in Simulation: ' simTag];
    for i=1:size(ef_all,4), ef_all(:,:,:,i) = ef_all(:,:,:,i).*nan_mask; end
    ef_mag = ef_mag.*nan_mask;
    dataShowVal = ef_mag(~isnan(ef_mag(:)));
    sliceshow(ef_mag,[],cm,[min(dataShowVal) prctile(dataShowVal,95)],'Electric field (V/m)',[figName '. Click anywhere to navigate.'],ef_all,mri2mni); drawnow
    
else
    
    for i=1:size(r.ef_all,4), r.ef_all(:,:,:,i) = r.ef_all(:,:,:,i).*nan_mask; end
    r.ef_mag = r.ef_mag.*nan_mask;
    dataShowVal = r.ef_mag(~isnan(r.ef_mag(:)));
    for i=1:size(r.targetCoord,1)
        figName = ['Electric field at Target ' num2str(i) ' in Targeting: ' tarTag];
        sliceshow(r.ef_mag,r.targetCoord(i,:),cm,[min(dataShowVal) prctile(dataShowVal,95)],'Electric field (V/m)',[figName '. Click anywhere to navigate.'],r.ef_all,mri2mni); drawnow
    end
end