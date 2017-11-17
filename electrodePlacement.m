function vargout = electrodePlacement(vargin)
% help text

%% input check

%% coord options
switch coordType
    case {'1020','1010','1005'}
%         or just detect by using electrode names
load('D_1005sysFULLextraAdded.mat')
    case 'biosemi'
%         or just detect by using electrode names
load('D_biosemiFULLextraAdded.mat')
    case 'customized'
end

%% shape options
switch shape
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

%% assume head is in RAS from the very beginning (using mricron to convert)
%% later remove this assumption (to be consistent with user-provided coordinates)
template = load_untouch_nii([dirname filesep baseFilename '_mask_skin.nii']); % Load the scalp mask
% template is used for saving the results with the same header info as the input
scalp = template.img;
[Nx, Ny, Nz]=size(scalp); % size of head in RAS orientation

scalp_surface = mask2EdgePointCloud(scalp);

%% fit cap position on the individual's head, if not customized coodinates

%% mapping landmarks
[dirname,baseFilename] = fileparts(P);
if isempty(dirname), dirname = pwd; end

load([dirname filesep baseFilename '_seg8.mat'],'image','tpm','Affine');
tpm2mri = inv(image(1).mat)*inv(Affine)*tpm(1).mat;
% mapping landmarks from TPM to individual MRI % ANDY 2017-05-17
landmark = [61 139 98; % nasion
    61 9 100; % inion
    11 62 93; % right
    111 63 93; % left; note here because eTPM is LAS orientation
    61 113 7; % front_neck
    61 7 20]; % back_neck
temp = tpm2mri * [landmark(1,:) 1]'; nasion = round(temp(1:3)');
temp = tpm2mri * [landmark(2,:) 1]'; inion = round(temp(1:3)');
temp = tpm2mri * [landmark(3,:) 1]'; right = round(temp(1:3)');
temp = tpm2mri * [landmark(4,:) 1]'; left = round(temp(1:3)');
temp = tpm2mri * [landmark(5,:) 1]'; front_neck = round(temp(1:3)');

electrode_coord = fitCap2individual(scalp,nasion,inion,right,left);

%% head clean up for placing electrodes
scalp_surface_clean = cleanScalp(scalp_surface);
scalp_surface2 = mask2EdgePointCloud(scalp_surface_clean);

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