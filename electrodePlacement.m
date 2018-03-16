function [rnge_elec,rnge_gel] = electrodePlacement(P,elecNeeded,elecPara,uniTag)
% help text

[dirname,baseFilename] = fileparts(P);
if isempty(dirname), dirname = pwd; end

%% options
capType = elecPara.capType;
% elecType = elecPara.elecType;
% elecSize = elecPara.elecSize;
% elecOri = elecPara.elecOri;
% legacy = elecPara.legacy;
doPredefined = elecPara.doPredefined;
doNeck = elecPara.doNeck;
doCustom = elecPara.doCustom;

%% cap options
if doPredefined
    switch capType
        case {'1020','1010','1005'}
            load('cap1005FullWithExtra.mat','capInfo');
            elecPool_P = capInfo{1};
            %         elec_template = cell2mat(capInfo(2:4));
            isBiosemi = 0;
            %         isCustomizedCap = 0;
            [isPredefined,indPredefined]=ismember(elecNeeded,elecPool_P);
            indP = indPredefined(isPredefined);
        case {'biosemi','Biosemi','bioSemi','BioSemi','BIOSEMI'}
            load('capBioSemiFullWithExtra.mat','capInfo');
            elecPool_P = capInfo{1};
            %         elec_template = cell2mat(capInfo(2:4));
            isBiosemi = 1;
            %         isCustomizedCap = 0;
            [isPredefined,indPredefined]=ismember(elecNeeded,elecPool_P);
            indP = indPredefined(isPredefined);
    end
% else
%     elecPool_P = '';
end

if doNeck
    elecPool_N = {'Nk1';'Nk2';'Nk3';'Nk4'};
    [isNeck,indNeck]=ismember(elecNeeded,elecPool_N);
    indN = indNeck(isNeck);
% else
%     elecPool_N = '';
end

if doCustom
    fid = fopen([dirname filesep baseFilename '_customLocations']);
    capInfo_C = textscan(fid,'%s %f %f %f');
    fclose(fid);
    elecPool_C = capInfo_C{1};
    elecLoc_C = cell2mat(capInfo_C(2:4));
    [isCustom,indCustom]=ismember(elecNeeded,elecPool_C);
    indC = indCustom(isCustom);
    elecLoc_C = elecLoc_C(indC,:);
% else
%     elecPool_C = '';
end

% doPredefined = 0;
% doNeck = 0;
% doCustom = 0;
% for i=1:length(elecNeeded)
%     if ismember(elecNeeded{i},elecPool)
%         doPredefined = 1;
%     end
%     if ismember(elecNeeded{i},{'Nk1';'Nk2';'Nk3';'Nk4'})
%         doNeck = 1;
%     end
%     if ~isempty(strfind(elecNeeded{i},'custom')) || ~isempty(strfind(elecNeeded{i},'Custom'))
%         doCustom = 1;
%         fid = fopen([dirname filesep 'customLocations']);
%         capInfoCustom = textscan(fid,'%s %f %f %f');
%         fclose(fid);
%         elecPoolCustom = capInfoCustom{1};
%         elecLocCustom = cell2mat(capInfoCustom(2:4));
%     end
% end

% isCustomElec = zeros(length(elecNeeded),1);
% for i=1:length(elecNeeded)
%     if ~isempty(strfind(elecNeeded{i},'custom')) || ~isempty(strfind(elecNeeded{i},'Custom'))
%         isCustomElec(i) = 1;
%     end
% end
% if all(isCustomElec)
%     pureCustom = 1;
% else
%     pureCustom = 0;
%
% switch capType
%     case {'1020','1010','1005'}
%         load('cap1005FullWithExtra.mat','capInfo');
%         elec = capInfo{1};
% %         elec_template = cell2mat(capInfo(2:4));
%         isBiosemi = 0;
% %         isCustomizedCap = 0;
%     case 'biosemi'
%         load('capBioSemiFullWithExtra.mat','capInfo');
%         elec = capInfo{1};
% %         elec_template = cell2mat(capInfo(2:4));
%         isBiosemi = 1;
% %         isCustomizedCap = 0;
% %     case 'customized'
% % %         read user-provided electrode coordniate text file
% %         elec = capInfo{1};
% %         elec_template = cell2mat(capInfo(2:4));
% %         isBiosemi = 0;
% %         isCustomizedCap = 1;
% end
% end

% %% electrode shape options
% switch elecType
%     case 'disc'
%         radius = elecSize(1);
%         height = elecSize(2);
%     case 'pad'
%         padSize = elecSize; % [length width thickness]
%         padOri = elecOri;
%     case 'ring'
%         radiusInner = elecSize(1);
%         radiusOutter = elecSize(2);
%         height = elecSize(3);
% end

%% can be any non-ras head (to be consistent with user-provided coordinates)
landmarks_original = getLandmarks(P);

[perm,iperm,isFlipInner,isFlipOutter] = how2getRAS(landmarks_original);

template = load_untouch_nii([dirname filesep baseFilename '_mask_skin.nii']); % Load the scalp mask
% template is used for saving the results with the same header info as the input
scalp_original = template.img;

scalp = changeOrientationVolume(scalp_original,perm,isFlipInner);

if doPredefined || doNeck
    landmarks = changeOrientationPointCloud(landmarks_original,perm,isFlipInner,size(scalp));
end
if doCustom
    elecLoc_C = changeOrientationPointCloud(elecLoc_C,perm,isFlipInner,size(scalp));
end

scalp_surface = mask2EdgePointCloud(scalp,'erode',ones(3,3,3));

%% fit cap position on the individual's head
if doPredefined
   [electrode_coord_P,center_P]= fitCap2individual(scalp,scalp_surface,landmarks,capInfo,indP,isBiosemi);
%     if ~isBiosemi
%         
%         if ~exist([dirname filesep baseFilename '_fittedCap1005.mat'],'file')
%             [electrode_coord_P,center_P]= fitCap2individual(scalp,scalp_surface,landmarks,capInfo,isBiosemi);
%             save([dirname filesep baseFilename '_fittedCap1005.mat'],'electrode_coord_P','center_P');
%         else
%             load([dirname filesep baseFilename '_fittedCap1005.mat'],'electrode_coord_P','center_P');
%         end
%         
%     else
%         
%         if ~exist([dirname filesep baseFilename '_fittedCapBioSemi.mat'],'file')
%             [electrode_coord_P,center_P]= fitCap2individual(scalp,scalp_surface,landmarks,capInfo,isBiosemi);
%             save([dirname filesep baseFilename '_fittedCapBioSemi.mat'],'electrode_coord_P','center_P');
%         else
%             load([dirname filesep baseFilename '_fittedCapBioSemi.mat'],'electrode_coord_P','center_P');
%         end
%         
%     end
else
    electrode_coord_P = []; center_P = [];
end

if doNeck
    if any(landmarks(5:6,3)<=0)
        error('MRI does not cover the neck, so cannot place electrodes on the neck.');
    else
        [electrode_coord_N,center_N] = placeNeckElec(scalp,scalp_surface,landmarks,indN);
    end
else
    electrode_coord_N = []; center_N = [];
end

if doCustom
    [~,indOnScalpSurf] = map2Points(elecLoc_C,scalp_surface,'closest');
    electrode_coord_C = scalp_surface(indOnScalpSurf,:);
    center_C = mean(scalp_surface);
else
    electrode_coord_C = []; center_C = [];
end

%% head clean up for placing electrodes
[scalp_clean,scalp_filled] = cleanScalp(scalp,scalp_surface);
scalp_clean_surface = mask2EdgePointCloud(scalp_clean,'erode',ones(3,3,3));

disp('calculating gel amount for each electrode...')
if doPredefined
    %     elec_range_P = zeros(size(electrode_coord_P,1),100);
    [~,indOnScalpSurf] = project2ClosestSurfacePoints(electrode_coord_P,scalp_clean_surface,center_P);
    elec_range_P = indOnScalpSurf(1:100,:);
    %     for i=1:size(elec_range_P,1)
    %         elec_range_P(i,:) = indOnScalpSurf(1:100,i);
    % Get some points on the scalp surface that are close to the exact
    % location of each electrode for the calculation of local normal vector
    % for each electrode in the following step
    %     end
else
    elec_range_P = [];
end

if doNeck
    % Get local scalp points for neck electrodes
    %     elec_range_N = zeros(size(electrode_coord_N,1),100);
    [~,indOnScalpSurf] = project2ClosestSurfacePoints(electrode_coord_N,scalp_clean_surface,center_N);
    elec_range_N = indOnScalpSurf(1:100,:);
    %     for i=1:size(elec_range_N,1)
    %         elec_range_N(i,:) = indOnScalpSurf(1:100,i);
    %     end
else
    elec_range_N = [];
end

if doCustom
    [~,elec_range_C] = map2Points(electrode_coord_C,scalp_clean_surface,'closer',100);
    % [~,elec_range_C] = map2Points(electrode_coord_C,scalp_surface,'closer',100);
else
    elec_range_C = [];
end

% elecPool = cat(1,elecPool_P,elecPool_N,elecPool_C);
electrode_coord = cat(1,electrode_coord_P,electrode_coord_N,electrode_coord_C);
% electrode_center = cat(1,repmat(center_P,size(electrode_coord_P,1),1),...
%     repmat(center_N,size(electrode_coord_N,1),1),repmat(center_C,size(electrode_coord_C,1),1));
elec_range = cat(1,elec_range_P',elec_range_N',elec_range_C');

% [~,indElecNeeded] = ismember(elecNeeded,elecPool);
% doElec = zeros(size(electrode_coord,1),1);
% doElec(indElecNeeded) = 1;

%% placing and model the electrodes
[elec_C,gel_C] = placeAndModelElectrodes(electrode_coord,elec_range,scalp_clean_surface,scalp_filled,elecPara);

%% generate final results (elec and gel masks, and their coordinate ranges)
disp('constructing electrode and gel volume to be exported...')
for i = 1:length(elec_C)
    if ~isempty(elec_C{i})
        elec_C{i} = changeOrientationPointCloud(elec_C{i},iperm,isFlipOutter,size(scalp_original));
        gel_C{i} = changeOrientationPointCloud(gel_C{i},iperm,isFlipOutter,size(scalp_original));
    end
end

[volume_elec,rnge_elec] = generateElecMask(elec_C,size(scalp_original));
[volume_gel,rnge_gel] = generateElecMask(gel_C,size(scalp_original));

disp('final clean-up...')
volume_gel = xor(volume_gel,volume_gel & scalp_original); % remove the gel that goes into the scalp
volume_gel = xor(volume_gel,volume_gel & volume_elec); % remove the gel that overlap with the electrode
bone = load_untouch_nii([dirname filesep baseFilename '_mask_bone.nii']); volume_bone = bone.img;
volume_gel = xor(volume_gel,volume_gel & volume_bone); % remove the gel that gets into the bone

% disp('saving the results...')
template.img = uint8(volume_elec)*255;
save_untouch_nii(template,[dirname filesep baseFilename '_' uniTag '_mask_elec.nii']);
template.img = uint8(volume_gel)*255;
save_untouch_nii(template,[dirname filesep baseFilename '_' uniTag '_mask_gel.nii']);

for i=1:length(rnge_elec)
    if ~isempty(rnge_elec{i})
        rnge_elec{i} = rnge_elec{i}.*repmat(template.hdr.dime.pixdim(2:4),2,1);
        rnge_gel{i} = rnge_gel{i}.*repmat(template.hdr.dime.pixdim(2:4),2,1);
    end
end % use NIFTI header info to convert range info into world coordinates for subsequent electrode labeling

save([dirname filesep baseFilename '_' uniTag '_rnge.mat'],'rnge_elec','rnge_gel');