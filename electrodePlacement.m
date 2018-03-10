function [rnge_elec,rnge_gel] = electrodePlacement(P,elecNeeded,elecPara,uniTag)
% help text

%% options
capType = elecPara.capType;
elecType = elecPara.elecType;
elecSize = elecPara.elecSize;
elecOri = elecPara.elecOri;
legacy = elecPara.legacy;

%% cap options
switch capType
    case {'1020','1010','1005'}
        load('cap1005FullWithExtra.mat','capInfo');
        elecPool = capInfo{1};
%         elec_template = cell2mat(capInfo(2:4));
        isBiosemi = 0;
%         isCustomizedCap = 0;
    case 'biosemi'
        load('capBioSemiFullWithExtra.mat','capInfo');
        elecPool = capInfo{1};
%         elec_template = cell2mat(capInfo(2:4));
        isBiosemi = 1;
%         isCustomizedCap = 0;
end

doPredefined = 0;
doNeck = 0;
doCustom = 0;
for i=1:length(elecNeeded)
    if ismember(elecNeeded{i},elecPool)
        doPredefined = 1;
    end
    if ismember(elecNeeded{i},{'Nk1';'Nk2';'Nk3';'Nk4'})
        doNeck = 1;
    end
    if ~isempty(strfind(elecNeeded{i},'custom')) || ~isempty(strfind(elecNeeded{i},'Custom'))
        doCustom = 1;
        fid = fopen([dirname filesep 'customLocations']);
        capInfoCustom = textscan(fid,'%s %f %f %f');
        fclose(fid);
        elecPoolCustom = capInfoCustom{1};
        elecLocCustom = cell2mat(capInfoCustom(2:4));
    end
end
        
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

[perm,iperm,isFlip] = how2getRAS(landmarks_original);

[dirname,baseFilename] = fileparts(P);
if isempty(dirname), dirname = pwd; end
template = load_untouch_nii([dirname filesep baseFilename '_mask_skin.nii']); % Load the scalp mask
% template is used for saving the results with the same header info as the input
scalp_original = template.img;

scalp = changeOrientationVolume(scalp_original,perm,isFlip);

if doPredefined || doNeck
    landmarks = changeOrientationPointCloud(landmarks_original,perm,isFlip,size(scalp));
end
if doCustom
    elecLocCustom = changeOrientationPointCloud(elecLocCustom,perm,isFlip,size(scalp));
end

scalp_surface = mask2EdgePointCloud(scalp,'erode',ones(3,3,3));

%% fit cap position on the individual's head
if doPredefined
    [electrode_coord,center]= fitCap2individual(scalp,scalp_surface,landmarks,capInfo,isBiosemi);
end

if doNeck
    if any(landmarks(5:6,3)<=0)
        error('MRI does not cover the neck, so cannot place electrodes on the neck.');
    else
        [neck_coord,neck_center] = placeNeckElec(scalp,scalp_surface,landmarks);
    end
end

if doCustom
    [~,indOnScalpSurf] = map2Points(elecLocCustom,scalp_surface,'closest');
    electrode_coord = scalp_surface(indOnScalpSurf,:);
end


%% head clean up for placing electrodes
[scalp_clean,scalp_filled] = cleanScalp(scalp,scalp_surface);
scalp_clean_surface = mask2EdgePointCloud(scalp_clean,'erode',ones(3,3,3));

disp('calculating gel amount for each electrode...')
if doPredefined
    elec_range = zeros(100,size(electrode_coord,1));
    [~,indOnScalpSurf] = project2ClosestSurfacePoints(electrode_coord,scalp_clean_surface,center);
    for i=1:size(elec_range,1)
        elec_range(:,i) = indOnScalpSurf(1:100,i);
        % Get some points on the scalp surface that are close to the exact
        % location of each electrode for the calculation of local normal vector
        % for each electrode in the following step
    end
end

if doNeck
    % Get local scalp points for neck electrodes    
    neck_elec_range = zeros(100,size(neck_coord,1));
    [~,indOnScalpSurf] = project2ClosestSurfacePoints(neck_coord,scalp_clean_surface,neck_center);
    for i=1:size(neck_elec_range,1)
        neck_elec_range(:,i) = indOnScalpSurf(1:100,i);
    end
end
    
if doCustom    
    [~,elec_range] = map2Points(electrode_coord,scalp_clean_surface,'closer',100);
end

% if isBiosemi THIS DOES NOT MATTER ANY MORE FOR ROAST
%     aidElec = [CPz FCz AFz POz];
%     elecToPlace = setdiff(1:size(electrode_coord,1),aidElec);
%     electrode_coord = electrode_coord(elecToPlace,:);
%     elec_range = elec_range(:,elecToPlace);
% end

if any(ismember(elecNeeded,{'Nk1','Nk2','Nk3','Nk4'}))
    elec = cat(1,elec,{'Nk1','Nk2','Nk3','Nk4'});
    electrode_coord = cat(1,electrode_coord,neck_coord);
    elec_range = cat(2,elec_range,neck_elec_range);
end

[~,indElecNeeded] = ismember(elecNeeded,elec);
doElec = zeros(size(electrode_coord,1),1);
doElec(indElecNeeded) = 1;

%% placing and model the electrodes
[elec_C,gel_C] = placeAndModelElectrodes(electrode_coord,elec_range,scalp_clean_surface,scalp_filled,doElec,elecPara);

%% generate final results (elec and gel masks, and their coordinate ranges)
disp('constructing electrode and gel volume to be exported...')
for i = 1:length(elec_C)
    if ~isempty(elec_C{i})
        elec_C{i} = changeOrientationPointCloud(elec_C{i},iperm,isFlip,size(scalp_original));
        gel_C{i} = changeOrientationPointCloud(gel_C{i},iperm,isFlip,size(scalp_original));
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