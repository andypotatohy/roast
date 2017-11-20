function vargout = electrodePlacement(P,elecNeeded,elecPara,uniTag)
% help text

%% input check
capType = elecPara.capType;
elecType = elecPara.elecType;
elecDim = elecPara.elecDim;
elecOri = elecPara.elecOri;

%% cap options
switch capType
    case {'1020','1010','1005'}
        load('D_1005sysFULLextraAdded.mat');
        elec = capInfo{1};
        elec_template = cell2mat(capInfo(2:4));
        isBiosemi = 0;
        isCustomizedCap = 0;
    case 'biosemi'
        load('D_biosemiFULLextraAdded.mat');
        elec = capInfo{1};
        elec_template = cell2mat(capInfo(2:4));
        isBiosemi = 1;
        isCustomizedCap = 0;
    case 'customized'
%         read user-provided electrode coordniate text file
        elec = capInfo{1};
        elec_template = cell2mat(capInfo(2:4));
        isBiosemi = 0;
        isCustomizedCap = 1;
end

%% electrode shape options
switch elecType
    case 'disc'
        radius = ;
        height = ;
    case 'pad'
        padDim = [70 50 3]; % [length width thickness]
        padOri = 'ap';
    case 'ring'
        radiusInner = ;
        radiusOutter = ;
        height = ;
end

%% can be any non-ras head (to be consistent with user-provided coordinates)
landmarks_original = getLandmarks(P);

[perm,iperm,isFlip] = how2getRAS(landmarks_original);

[dirname,baseFilename] = fileparts(P);
if isempty(dirname), dirname = pwd; end
template = load_untouch_nii([dirname filesep baseFilename '_mask_skin.nii']); % Load the scalp mask
% template is used for saving the results with the same header info as the input
scalp_original = template.img;

[scalp,landmarks,sizRAS] = changeOrientationIn(perm,isFlip,scalp_original,landmarks_original);
if isCustomizedCap
    [~,elec_template] = changeOrientationIn(perm,isFlip,scalp_original,elec_template);
end

scalp_surface = mask2EdgePointCloud(scalp,'erode',ones(3,3,3));

%% fit cap position on the individual's head
if ~isCustomizedCap
    electrode_coord = fitCap2individual(scalp,scalp_surface,landmarks,capInfo,isBiosemi,elecNeeded);
else
    electrode_coord = map2ClosestPoints(elec_template,scalp_surface);
end

%% head clean up for placing electrodes
[scalp_clean,scalp_filled] = cleanScalp(scalp);
scalp_surface2 = mask2EdgePointCloud(scalp_clean,'erode',ones(3,3,3));

disp('calculating gel amount for each electrode...')
vec1 = repmat(center,size(electrode_coord,1),1)-electrode_coord;
% vectors connecting center to each electrode
vec2 = repmat(center,size(scalp_surface2,1),1)-scalp_surface2;
% vectors connecting center to each point on scalp surface (after clean-up)
elec_range = zeros(100,size(vec1,1));
for j=1:size(vec1,1)
    temp = dot(repmat(vec1(j,:),size(vec2,1),1),vec2,2)./(repmat(norm(vec1(j,:)),size(vec2,1),1).*sqrt(sum(vec2.^2,2)));
    [~,intemp] = sort(temp,'descend');
    elec_range(:,j) = intemp(1:100);
    % Get some points on the scalp surface that are close to the exact
    % location of each electrode for the calculation of local normal vector
    % for each electrode in the following step
end

if ~isempty(back_neck)
    vec1 = repmat(neckCenter,size(neck_coord,1),1)-neck_coord;
    vec2 = repmat(neckCenter,size(scalp_surface2,1),1)-scalp_surface2;
    neck_elec_range = zeros(100,size(vec1,1));
    for j=1:size(vec1,1)
        temp = dot(repmat(vec1(j,:),size(vec2,1),1),vec2,2)./(repmat(norm(vec1(j,:)),size(vec2,1),1).*sqrt(sum(vec2.^2,2)));
        [~,intemp] = sort(temp,'descend');
        neck_elec_range(:,j) = intemp(1:100);
    end
    electrode_coord = cat(1,electrode_coord,neck_coord);
    elec_range = cat(2,elec_range,neck_elec_range);
end
% Get local scalp points for neck electrodes

if isBiosemi
    aidElec = [CPz FCz AFz POz];
    elecToPlace = setdiff(1:size(electrode_coord,1),aidElec);
    electrode_coord = electrode_coord(elecToPlace,:);
    elec_range = elec_range(:,elecToPlace);
end

[~,indElecNeeded] = ismember(elecNeeded,elec);
doElec = zeros(size(electrode_coord,1),1);
doElec(indElecNeeded) = 1;

%% placing and model the electrodes

%% range generation

disp('constructing electrode and gel volume to be exported...')
for i = 1:size(electrode_coord,1)
    if doElec(i)
        for j=1:size(gel_C{i},1)
            if gel_C{i}(j,1)>0 && gel_C{i}(j,1)<=Nx && gel_C{i}(j,2)>0 && gel_C{i}(j,2)<=Ny && gel_C{i}(j,3)>0 && gel_C{i}(j,3)<=Nz
                volume_gel(gel_C{i}(j,1), gel_C{i}(j,2), gel_C{i}(j,3)) = 1;
            end
        end
        for j=1:size(elec_C{i},1)
            if elec_C{i}(j,1)>0 && elec_C{i}(j,1)<=Nx && elec_C{i}(j,2)>0 && elec_C{i}(j,2)<=Ny && elec_C{i}(j,3)>0 && elec_C{i}(j,3)<=Nz
                volume_elec(elec_C{i}(j,1), elec_C{i}(j,2), elec_C{i}(j,3)) = 1;
            end
        end
    end
end
volume_gel = permute(volume_gel,iperm); volume_elec = permute(volume_elec,iperm);
% permute back electrode and gel coordinates to match the original orientation of the head
% Construct electrode and gel volume

disp('final clean-up...')
scalp = template.img;
volume_gel = xor(volume_gel,volume_gel & scalp); % remove the gel that goes into the scalp
volume_gel = xor(volume_gel,volume_gel & volume_elec); % remove the gel that overlap with the electrode
bone = load_untouch_nii([dirname filesep baseFilename '_mask_bone.nii']); volume_bone = bone.img;
volume_gel = xor(volume_gel,volume_gel & volume_bone); % remove the gel that gets into the bone

% disp('exporting the results...')
template.img = im2uint8(volume_gel); save_untouch_nii(template,[dirname filesep baseFilename '_' uniTag '_mask_gel.nii']);
template.img = im2uint8(volume_elec); save_untouch_nii(template,[dirname filesep baseFilename '_' uniTag '_mask_elec.nii']);
% Save the results
% disp('DONE! (electrodes and gel were exported as mask_elec.nii and mask_gel.nii)')