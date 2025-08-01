function [elec,gel] = electrodePlacement(subj,template,imgHdr,landmarks,elecNeeded,options,uniTag)
% [elec,gel] = electrodePlacement(subj,template,imgHdr,landmarks,elecNeeded,options,uniTag)
%
% Place electrodes on the scalp surface. options.elecPara contains all the options
% info for each electrode. Enforced RAS in the first step starting from ROAST v3.0
% 
% (c) Yu (Andy) Huang, Parra Lab at CCNY
% yhuang16@citymail.cuny.edu
% April 2018

[dirname,subjName] = fileparts(subj);
if isempty(dirname), dirname = pwd; end

% scalp = template.img==5;
scalp=template.img>0; % fill in scalp first to avoid complication % ANDY 2024-03-12
scalp_surface = mask2EdgePointCloud(scalp,'erode',ones(3,3,3));

elecPara = options.elecPara;

indP = elecPara(1).indP;
indN = elecPara(1).indN;
indC = elecPara(1).indC;

%% fit cap position on the individual's head
if ~isempty(indP)
   switch lower(elecPara(1).capType)
       case {'1020','1010','1005'}
           load('./cap1005FullWithExtra.mat','capInfo');
           isBiosemi = 0;
           isEGI = 0;
       case 'biosemi'
           load('./capBioSemiFullWithExtra.mat','capInfo');
           isBiosemi = 1;
           isEGI = 0;
       case 'egi'
           load('./capEGIfull.mat','capInfo');
           isBiosemi = 0;
           isEGI = 1;
   end
   [electrode_coord_P,center_P]= fitCap2individual(scalp,scalp_surface,landmarks,imgHdr,capInfo,indP,isBiosemi,isEGI);
else
    electrode_coord_P = []; center_P = [];
end

if ~isempty(indN)
    if any(landmarks(5:6,3)<=0)
        error('MRI does not cover the neck, so cannot place electrodes on the neck. Consider using ''zeroPadding'' option to extend the input MRI.');
    else
        [electrode_coord_N,center_N] = placeNeckElec(scalp,scalp_surface,landmarks,indN);
    end
else
    electrode_coord_N = []; center_N = [];
end

if ~isempty(indC)
    fid = fopen([dirname filesep subjName '_customLocations']);
    capInfo_C = textscan(fid,'%s %f %f %f');
    fclose(fid);
    elecLoc_C = cell2mat(capInfo_C(2:4));
    % update customized elec coords according to options of isNonRAS, resamp, and zeroPad
    if options.isNonRAS
        [elecLoc_C,perm] = convertToRASpointCloud(subj,elecLoc_C);
    else
        perm = [1 2 3];
    end
    if options.resamp
        data = load_untouch_nii(subj);
        temp = data.hdr.dime.pixdim(2:4);
        temp = temp(perm);
        elecLoc_C = elecLoc_C.*repmat(temp,size(elecLoc_C,1),1);
    end
    if options.zeroPad>0
        elecLoc_C = elecLoc_C + options.zeroPad;
    end
    
    elecLoc_C = elecLoc_C(indC,:);
    [~,indOnScalpSurf] = map2Points(elecLoc_C,scalp_surface,'closest');
    electrode_coord_C = scalp_surface(indOnScalpSurf,:);
else
    electrode_coord_C = [];
end

%% head clean up for placing electrodes
[scalp_clean,scalp_filled] = cleanScalp(scalp,scalp_surface);
temp1 = scalp_filled(:,:,[1 end]);
temp2 = scalp_filled(:,[1 end],:);
temp3 = scalp_filled([1 end],:,:);
if any([temp1(:);temp2(:);temp3(:)])
    warning('Scalp touches image boundary. Electrodes may go out of image boundary. ROAST can continue but results may not be accurate. It is recommended that you expand the input MRI by specifying the ''zeroPadding'' option.');
end    
scalp_clean_surface = mask2EdgePointCloud(scalp_clean,'erode',ones(3,3,3));

disp('calculating gel amount for each electrode...')
if ~isempty(indP)
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

if ~isempty(indN)
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

if ~isempty(indC)
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

%% placing and model the electrodes
resolution = mean([imgHdr(1).mat(1,1),imgHdr(1).mat(2,2),imgHdr(1).mat(3,3)]);
% mean() here to handle anisotropic resolution; ugly. Maybe just
% resample MRI to isotropic in the very beginning?
[elec_C,gel_C] = placeAndModelElectrodes(electrode_coord,elec_range,scalp_clean_surface,scalp_filled,elecNeeded,elecPara,resolution,0,uniTag);

%% generate final results (elec and gel masks, and their coordinate ranges)
disp('constructing electrode and gel volume to be exported...')
% for i = 1:length(elec_C)
%     if ~isempty(elec_C{i})
%         elec_C{i} = changeOrientationPointCloud(elec_C{i},iperm,isFlipOutter,size(scalp_original));
%         gel_C{i} = changeOrientationPointCloud(gel_C{i},iperm,isFlipOutter,size(scalp_original));
%     end
% end

% [volume_elec,volume_elecLabel] = generateElecMask(elec_C,size(scalp_original),elecNeeded,1);
% [volume_gel,volume_gelLabel] = generateElecMask(gel_C,size(scalp_original),elecNeeded,0);
% volume_elec_C = generateElecMask(elec_C,size(scalp_original),elecNeeded,1);
% volume_gel_C = generateElecMask(gel_C,size(scalp_original),elecNeeded,0);
volume_elec_C = generateElecMask(elec_C,size(scalp),elecNeeded,1);
volume_gel_C = generateElecMask(gel_C,size(scalp),elecNeeded,0);

disp('final clean-up...')
volume_elec = volume_elec_C>0;
volume_gel = volume_gel_C>0;
volume_gel = xor(volume_gel,volume_gel & volume_elec); % remove the gel that overlaps with the electrode
for i=1:6
    volume_tissue = template.img==i;
    volume_gel = xor(volume_gel,volume_gel & volume_tissue);
end % remove the gel that goes into other tissue masks

disp('saving placed electrodes...')
elec = template;
elec.fileprefix = [dirname filesep subjName '_' uniTag '_mask_elec'];
elec.hdr.hist.descrip = 'electrode mask';
% elec.img = uint8(volume_elec)*255;
elec.img = uint8(volume_elec_C.*volume_elec);
save_untouch_nii(elec,[dirname filesep subjName '_' uniTag '_mask_elec.nii']);

gel = template;
gel.fileprefix = [dirname filesep subjName '_' uniTag '_mask_gel'];
gel.hdr.hist.descrip = 'gel mask';
% gel.img = uint8(volume_gel)*255;
gel.img = uint8(volume_gel_C.*volume_gel);
save_untouch_nii(gel,[dirname filesep subjName '_' uniTag '_mask_gel.nii']);

% save([dirname filesep subjName '_' uniTag '_labelVol.mat'],'volume_elecLabel','volume_gelLabel');