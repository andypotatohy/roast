function [rnge_elec,rnge_gel] = electrode_placement(P,model,elec_radius,gel_height,elec_height,elecNeeded,uniTag)
% [rnge_elec,rnge_gel] = electrode_placement(P,model,elec_radius,gel_height,elec_height,elecNeeded,uniTag)
%
% Automated electrode placement on scalp surface. It can automatically
% place electrodes on the scalp surface following standard international
% 10-10 system.
%
% P: the T1 image file, giving directory of the data;
%
% elec_radius,gel_height,elec_height: the radius of each electrode, the
% height of each electrode and the underlying gel. Default values: 6 mm
% radius for each electrode with a height of 2 mm, and 2 mm height for gel.
%
% Output: electrode and gel masks, i.e., mask_elec.nii and mask_gel.nii
%
% rnge_elec,rnge_gel: for labelling electrodes later.
% 
% See Huang et al 2013 (DOI: 10.1088/1741-2560/10/6/066004) for details.
%
% (c) Yu Huang (Andy), May 2017
% (c) Yu Huang (Andy), September 2012
% (c) Yu Huang (Andy), March 2011
% (c) Yuzhuo Su (Suzy), May 2010
% The Neural Engineering Lab, Dept. of Biomedical Engineering, City College of New York
% Send bugs to yhuang16@citymail.cuny.edu

if nargin < 3 || isempty(elec_radius)
    elec_radius = 6; % default electrode radius
end

if nargin < 4 || isempty(gel_height)
    gel_height = 2; % default gel height
end

if nargin < 5 || isempty(elec_height)
    elec_height = 2; % default electrode height
end

switch model
    case '1020'
        % something for future products
    case '1010' % the default option
        capType = 'D_74.mat';% FUTURE: consider using D_72 where TP9/TP10 are removed (otherwise complications on ears)
        %         numOfElec = 74;
    case '1010X' % extended 10-10 system, i.e., 93-electrode model
        capType = 'D_89.mat';
        %         numOfElec = 93;
    case '1005'
        % something for future products
    case 'biosemi256'
        % something for future products
    case 'biosemi256X' % extended BioSemi256 system, i.e., 336-electrode model
        capType = 'D_336.mat';
        %         numOfElec = 336;
end

disp('loading data...')
% cd(dirname)
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
if ~isempty(strfind(model,'X'))
    % use this 'X' to tell the program if neck electrodes are needed to be placed
    temp = tpm2mri * [landmark(6,:) 1]'; back_neck = round(temp(1:3)');
else
    back_neck = [];
end

template = load_untouch_nii([dirname filesep baseFilename '_mask_skin.nii']); % Load the scalp mask
% template is used for saving the results with the same header info as the input
scalp = template.img;

disp('adjusting the orientation of the head into RAS orientation...')
e1 = right-left; e1 = e1/norm(e1);
e2 = nasion-inion; e2 = e2/norm(e2);
e3 = nasion-front_neck; e3 = e3/norm(e3); % detect the orientation of the head based on the anatomical landmarks

[~,perm1] = max(abs(e1)); [~,perm2] = max(abs(e2)); [~,perm3] = max(abs(e3));
isFlip = [sign(e1(perm1)) sign(e2(perm2)) sign(e3(perm3))]; % detect if the head is flipped or not in each direction compared to RAS system

perm = [perm1,perm2,perm3]; % permutation order into RAS
[~,iperm] = sort(perm); % inverse permutation order back to original orientation of the head
scalp = permute(scalp,perm); % permute the head into RAS
[Nx, Ny, Nz]=size(scalp); % size of head in RAS orientation

nasion = nasion(perm); inion = inion(perm); right = right(perm); left = left(perm);
front_neck = front_neck(perm);
if ~isempty(back_neck), back_neck = back_neck(perm); end % permute the landmarks accordingly

for i = 1:length(isFlip)
    if isFlip(i) < 0
        scalp = flipdim(scalp,i); % flip the head in flipped direction
        switch i
            case 1
                temp = Nx;
            case 2
                temp = Ny;
            case 3
                temp = Nz;
        end
        nasion(i) = temp - nasion(i) + 1; inion(i) = temp - inion(i) + 1;
        right(i) = temp - right(i) + 1; left(i) = temp - left(i) + 1;
        front_neck(i) = temp - front_neck(i) + 1;
        if ~isempty(back_neck), back_neck(i) = temp - back_neck(i) + 1; end % flip the landmarks accordingly
        % Note the +1
    end
end

disp('initial calibration of the head...')
s = right-left; s = s/norm(s);
c = nasion-inion; c = c/norm(c);
a = cross(s,c); a = a/norm(a); % vectors to be used in affine transform

scalp_edge = scalp-imerode(scalp,ones(3,3,3)); % Get the edge of scalp
inde = find(scalp_edge==255);
scalp_surface = zeros(length(inde),3);
[scalp_surface(:,1),scalp_surface(:,2),scalp_surface(:,3)] = ind2sub(size(scalp_edge),inde);
% Obtain scalp surface points as a matrix

disp('measuring head size...')
L = norm(inion-nasion); % Distance between nasion and inion
line_center = (inion+nasion)/2; % Midpoint between nasion and inion

centralSag = round(line_center(1));
img_c = squeeze(scalp(centralSag,:,:)); % The central sagittal slice
indc = find(sum(img_c)>0,1,'first');
indr1 = find(img_c(:,indc)>0,1,'first');
indr2 = find(img_c(:,indc)>0,1,'last');
img_c(indr1:indr2,indc) = 255;
% Make sure the central sagittal slice can be closed completely in order to detect the edge correctly
se = 0;
isFilled = 0;
[ytemp,ztemp] = ind2sub(size(img_c),find(img_c));
centroid = round(mean([ytemp ztemp]));
while ~all(isFilled)
    se = se+8;
    im_test = imfill(imclose(img_c,ones(se,se)),'holes');
    isFilled = im_test(centroid(1),centroid(2));
    % isFilled = [im_test(round(Ny/2),round(Nz/2));
    %     im_test(round(Ny/2)+5,round(Nz/2));
    %     im_test(round(Ny/2),round(Nz/2)+5)];
end
% bw_c = edge(imopen(imfill(imclose(img_c,ones(8,8)),'holes'),ones(3,3)));
bw_c = edge(imopen(im_test,ones(3,3)));
% Get the edge of the central sagittal slice
[r_c,c_c] = find(bw_c==1);

indxinion = find(inion(3)==c_c,1,'first');
indxnasion = find(nasion(3)==c_c,1,'last');
[~,I] = max(c_c);
temp_right_up = find((c_c>=c_c(indxinion))&(r_c<(r_c(I))));
[~,I_up] = sort(c_c(temp_right_up));
temp_right_up = temp_right_up(I_up);
temp_right_down = find((c_c>=c_c(indxnasion))&(r_c>=(r_c(I))));
[~,I_down] = sort(c_c(temp_right_down),'descend');
temp_right_down = temp_right_down(I_down);
index = [temp_right_up; temp_right_down];
% Preparation for the calculation of the distance between nasion and inion
% along scalp surface using Natural Cubic Spline

[bx,by,finalbreaks]=ncs2dapprox(r_c(index),c_c(index));
% Approximation of 2-D Data by Natural Cubic Spline
% http://www.mathworks.co.jp/matlabcentral/fileexchange/7617
t = finalbreaks';
pp1= spline(t,[bx,by]');
range = linspace(1,finalbreaks(end),finalbreaks(end));
yizi = ppval(pp1,range);
yi=yizi(1,:)';
zi=yizi(2,:)';

distance_all = sum(sqrt(diff(yi).^2+diff(zi).^2));
% Calculate the distance between nasion and inion along the scalp surface

disp('wearing the cap...')
load(fullfile(fileparts(which(mfilename)),capType),'D'); %EasyCap electrode coordinates

elec = cell(length(D),1);
for i =1:length(D),
    elec{i} = D(i).labels;
    if strcmp(elec{i},'Oz'),Oz = i;end
    if strcmp(elec{i},'POz'),POz = i;end
    if strcmp(elec{i},'Pz'),Pz = i;end
    if strcmp(elec{i},'CPz'),CPz = i;end
    if strcmp(elec{i},'Cz'),Cz = i;end
    if strcmp(elec{i},'FCz'),FCz = i;end
    if strcmp(elec{i},'Fz'),Fz = i;end
    if strcmp(elec{i},'AFz'),AFz = i;end
    if strcmp(elec{i},'Fpz'),Fpz = i;end
end
% Find the electrodes along the central sagittal line

isBiosemi = 0;
if ~exist('Oz','var') % this is the case of placing a BioSemi Cap
    isBiosemi = 1;
    for i =1:length(D),
        elec{i} = D(i).labels;
        if strcmp(elec{i},'A19'),Oz = i;end
        if strcmp(elec{i},'POz'),POz = i;end
        if strcmp(elec{i},'A6'),Pz = i;end
        if strcmp(elec{i},'CPz'),CPz = i;end
        if strcmp(elec{i},'A1'),Cz = i;end
        if strcmp(elec{i},'FCz'),FCz = i;end
        if strcmp(elec{i},'E17'),Fz = i;end
        if strcmp(elec{i},'AFz'),AFz = i;end
        if strcmp(elec{i},'E12'),Fpz = i;end
    end
    
    if exist('Oz','var') && ~exist('POz','var')
        errMsg = 'BioSemi cap detected...';
        errMsg = [errMsg 'If you are placing a BioSemi cap, please add the following electrodes to help the program optimally place all the electrodes:\n'];
        errMsg = [errMsg 'electrode name \t inclination \t azimuth \n'];
        errMsg = [errMsg 'POz \t\t 69 \t\t -90 \n'];
        errMsg = [errMsg 'CPz \t\t 23 \t\t -90 \n'];
        errMsg = [errMsg 'FCz \t\t 23 \t\t 90 \n'];
        errMsg = [errMsg 'AFz \t\t 69 \t\t 90 \n'];
        error('error:convert',errMsg);
    end
end % Find the electrodes along the central sagittal line

if ~exist('Oz','var')
    error('Unrecognized electrode cap!');
end

elec_template = ones(4,length(D));
for i = 1:length(D), elec_template(1:3,i) = [D(i).X;D(i).Y;D(i).Z]; end
% Electrode coordinates from EasyCap

theta = 23;
alpha = ((360-10*theta)/2)*(pi/180);
h = (L/2)*(1/tan(alpha));
% For the calculation of the center of electrode coordinates

disp('adjust the cap for optimized position...this will take a while...')
factor = 1:-0.05:0.5; % Adjusting factor
CENTER = zeros(length(factor),3);
ELEC_COORD = zeros(length(D),3,length(factor));
F = zeros(length(factor),1);
for n = 1:length(factor)
    %     fprintf('Iteration No. %d...\n', n)
    center = line_center + h*factor(n)*a; % Adjust the center
    CENTER(n,:) = center; % buffer
    scale = 150;
    shift = center';
    
    affine = scale * [s' c' a' shift/scale;0 0 0 1/scale]; % Affine transform matrix
    
    elec_adjusted = elec_template;
    elec_adjusted(3,:) = elec_template(3,:)*factor(n);
    % Adjust the z-coordinate correspondingly
    elec_transformed = affine * elec_adjusted;
    % Affine transform the EasyCap coordinate to an approximate position for each electrode outside of the scalp surface
    
    vec1 = repmat(center,size(elec_transformed,2),1)-elec_transformed(1:3,:)';
    % vectors connecting center to the approximate position for each electrode
    vec2 = repmat(center,size(scalp_surface,1),1)-scalp_surface;
    % vectors connecting center to each point on scalp surface
    idx = zeros(size(vec1,1),1);
    for j=1:size(vec1,1)
        temp = dot(repmat(vec1(j,:),size(vec2,1),1),vec2,2)./(repmat(norm(vec1(j,:)),size(vec2,1),1).*sqrt(sum(vec2.^2,2)));
        [sorttemp,intemp] = sort(temp,'descend');
        testPts = scalp_surface(intemp(sorttemp> max(sorttemp)*0.99993),:);
        vecT = repmat(center,size(testPts,1),1)-testPts;
        dist = sqrt(sum(vecT.^2,2));
        idx(j) = intemp(find(dist==max(dist),1,'first'));
        % Find the only point on the outer surface of the scalp for each electrode, i.e., the exact coordinates for each electrode on the scalp surface
    end
    
    elec_interp = scalp_surface(idx,:); % exact coordinates for each electrode
    ELEC_COORD(:,:,n) = elec_interp; % buffer
    
    center_points = [inion;elec_interp(Oz,:);elec_interp(POz,:);elec_interp(Pz,:);...
        elec_interp(CPz,:);elec_interp(Cz,:);elec_interp(FCz,:);elec_interp(Fz,:);...
        elec_interp(AFz,:);elec_interp(Fpz,:);nasion]; % coordinates for electrodes on central sagittal line
    center_fit = [centralSag*ones(length(yi),1) yi zi]; % coordinates for each point on central sagittal line
    
    indxsave = 1;
    distance = zeros(size(center_points,1)-1,1);
    for ii = 2:size(center_points,1)
        electemp = repmat(center_points(ii,:),[size(center_fit,1) 1]);
        [~,indx] = min(sqrt(sum((electemp - center_fit).^2,2)));
        distance(ii-1) = sum(sqrt(diff(yi(indxsave:indx)).^2+diff(zi(indxsave:indx)).^2));
        % Calculate the distance between every adjacent pair of electrodes on central sagittal line
        indxsave = indx;
    end
    F(n) = sum(abs(distance-distance_all/10)); % the total error compared to ideal 10-10 system, this error needs to be minimized
    % Optimize the location of the center, to make the distance between each adjacent electrode on central sagittal line equal to
    % 1/10 of the total distance from nasion to inion
end

[~,index] = min(F);
electrode_coord = ELEC_COORD(:,:,index); % exact coordinate for each electrode projected on the scalp surface
center = CENTER(index,:); % center of electrode coordinates

if ~isempty(back_neck)
    neckCenter = (front_neck+back_neck)/2;
    vec1 = [neckCenter-front_neck;neckCenter-back_neck;-1 0 0;1 0 0];
    vec2 = repmat(neckCenter,size(scalp_surface,1),1)- scalp_surface;
    idx = zeros(size(vec1,1),1);
    for j=1:size(vec1,1)
        temp = dot(repmat(vec1(j,:),size(vec2,1),1),vec2,2)./(repmat(norm(vec1(j,:)),size(vec2,1),1).*sqrt(sum(vec2.^2,2)));
        [sorttemp,intemp] = sort(temp,'descend');
        testPts = scalp_surface(intemp(sorttemp> max(sorttemp)*0.99993),:);
        vecT = repmat(neckCenter,size(testPts,1),1)-testPts;
        dist = sqrt(sum(vecT.^2,2));
        idx(j) = intemp(find(dist==max(dist),1,'first'));
    end
    neck_coord = scalp_surface(idx,:);
end
% Place neck electrodes if needed

disp('cleaning the hair for gel injection...')
scalp(min(scalp_surface(:,1)):max(scalp_surface(:,1)),min(scalp_surface(:,2)):max(scalp_surface(:,2)),min(scalp_surface(:,3))) = 255;
% force to close the most bottom slice, because some MRIs with
% limited FOV just cut off in the middle of the face
se = 0;
% isFilled = zeros(4,1);
isFilled = 0;
[xtemp,ytemp,ztemp] = ind2sub(size(scalp),find(scalp));
centroid = round(mean([xtemp ytemp ztemp]));
while ~all(isFilled)
    se = se+10;
    img1 = imfill(imclose(scalp,ones(se,se,se)),'holes');
    isFilled = img1(centroid(1),centroid(2),centroid(3));
    %     isFilled = [img1(round(Nx/2),round(Ny/2),round(Nz/2));
    %         img1(round(Nx/2)+5,round(Ny/2),round(Nz/2));
    %         img1(round(Nx/2),round(Ny/2)+5,round(Nz/2));
    %         img1(round(Nx/2),round(Ny/2),round(Nz/2)+5)]; % Make sure the scalp is closed and filled completely
end

se = 0;
isOpen = 1;
while any(isOpen)
    se = se+10;
    if se>30
        img1(:,:,Nz) = 0; img1(:,1,:) = 0; img1(:,Ny,:) = 0; img1(1,:,:) = 0; img1(Nx,:,:) = 0;
        % force to make the image boundaries to be "opened"
    end
    img2 = imopen(img1,ones(se,se,se));
    imTemp = squeeze(img2(:,:,Nz)); isOpen = imTemp(:);
    imTemp = squeeze(img2(:,1,:)); isOpen = [isOpen;imTemp(:)];
    imTemp = squeeze(img2(1,:,:)); isOpen = [isOpen;imTemp(:)];
    imTemp = squeeze(img2(Nx,:,:)); isOpen = [isOpen;imTemp(:)]; % Make sure the scalp is "opened" to have a smooth outer surface
end
% Scalp clean-up: for calculation of local normal vector for each electrode

img_edge = img2-imerode(img2,ones(3,3,3)); % Get the edge of scalp (after clean-up)
inde = find(img_edge==255);
scalp_surface2 = zeros(length(inde),3);
[scalp_surface2(:,1),scalp_surface2(:,2),scalp_surface2(:,3)] = ind2sub(size(img_edge),inde);
% Obtain scalp surface points as a matrix (after clean-up)

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

% disp(['placing electrodes...(in total ' num2str(size(electrode_coord,1)) ' electrodes)'])
disp('placing electrodes...')
NOP = 500; verSamp = 10;
r = 0.05:0.05:elec_radius; % parameters used for modeling of electrodes and gel
% gel_height = 2; elec_height = 2; % heights of electrodes and gel
volume_gel = zeros(size(scalp)); volume_elec = zeros(size(scalp)); % electrode and gel volume to be output

% figure;hold on;plot3(scalp_surface(:,1),scalp_surface(:,2),scalp_surface(:,3),'y.');
gel_C = cell(1,size(electrode_coord,1)); elec_C = cell(1,size(electrode_coord,1));
% buffer for coordinates of each electrode and gel point
for i = 1:size(electrode_coord,1)
    if doElec(i)
        %     fprintf('%d... ',i);
        %     if mod(i,10) == 0, fprintf('\n'); end
        lcl = scalp_surface2(elec_range(:,i),:); % local scalp surface for each electrode
        %     plot3(lcl(:,1),lcl(:,2),lcl(:,3),'b.');
        [U,D] = eig(cov(lcl)); [~,ind] = min(diag(D));
        nv = U(:,ind)'; normal = nv/norm(nv); % Local normal for each electrode
        
        gel_out = electrode_coord(i,:) +  2*gel_height*normal;
        electrode = gel_out + elec_height*normal;
        gel_in = gel_out - 4*gel_height*normal; % coordinates of the boundaries of gel and electrode
        if norm(center - gel_out) < norm(center - electrode_coord(i,:))
            normal = -normal;
            gel_out = electrode_coord(i,:) +  2*gel_height*normal;
            electrode = gel_out + elec_height*normal;
            gel_in = gel_out - 4*gel_height*normal;
        end % make sure the normal is pointing out
        
        gel_X = zeros(length(r)*verSamp*4,NOP); gel_Y = zeros(length(r)*verSamp*4,NOP); gel_Z = zeros(length(r)*verSamp*4,NOP);
        elec_X = zeros(length(r)*verSamp,NOP); elec_Y = zeros(length(r)*verSamp,NOP); elec_Z = zeros(length(r)*verSamp,NOP);
        for j = 1:length(r)
            [gel_X(((j-1)*verSamp*4+1):verSamp*4*j,:), gel_Y(((j-1)*verSamp*4+1):verSamp*4*j,:), gel_Z(((j-1)*verSamp*4+1):verSamp*4*j,:)] = cylinder2P(ones(verSamp*4)*r(j),NOP,gel_in,gel_out);
            [elec_X(((j-1)*verSamp+1):verSamp*j,:), elec_Y(((j-1)*verSamp+1):verSamp*j,:), elec_Z(((j-1)*verSamp+1):verSamp*j,:)] = cylinder2P(ones(verSamp)*r(j),NOP,gel_out,electrode);
        end % Use cylinders to model electrodes and gel, and calculate the coordinates of the points that make up the cylinder
        
        gel_coor = floor([gel_X(:) gel_Y(:) gel_Z(:)]);
        gel_coor = unique(gel_coor,'rows');
        elec_coor = floor([elec_X(:) elec_Y(:) elec_Z(:)]);
        elec_coor = unique(elec_coor,'rows'); % clean-up of the coordinates
        
%         plot3(elec_coor(:,1),elec_coor(:,2),elec_coor(:,3),'.b');
%         plot3(gel_coor(:,1),gel_coor(:,2),gel_coor(:,3),'.m');
        
        gel_C{i} = gel_coor; elec_C{i} = elec_coor; % buffer for coordinates of each electrode and gel point
        %     fprintf('%d out of %d electrodes placed...\n',i,size(electrode_coord,1));
    end
end
% xlabel('x');ylabel('y');zlabel('z'); view([270 0]);
% hold off; % Place electrodes and visualize the results

% fprintf('\n');
disp('post-processing...')
for i = 1:length(isFlip)
    if isFlip(i) < 0
        switch i
            case 1
                temp = Nx;
            case 2
                temp = Ny;
            case 3
                temp = Nz;
        end
        for j = 1:size(electrode_coord,1)
            if doElec(j)
                elec_C{j}(:,i) = temp - elec_C{j}(:,i) + 1;
                gel_C{j}(:,i) = temp - gel_C{j}(:,i) + 1;
            end
        end
    end % flip back electrode and gel coordinates to match the original orientation of the head
end

elec_C_final = cell(1,size(electrode_coord,1));
rnge = cell(1,size(electrode_coord,1));
% isOut = zeros(length(rnge),1);
for i = 1:size(electrode_coord,1)
    if doElec(i)
        elec_C_final{i} = zeros(size(elec_C{i},1),3);
        rng = zeros(2,3);
        for j=1:size(elec_C{i},1)
            if elec_C{i}(j,1)>0 && elec_C{i}(j,1)<=Nx && elec_C{i}(j,2)>0 && elec_C{i}(j,2)<=Ny && elec_C{i}(j,3)>0 && elec_C{i}(j,3)<=Nz
                elec_C_final{i}(j,1) = elec_C{i}(j,iperm(1));
                elec_C_final{i}(j,2) = elec_C{i}(j,iperm(2));
                elec_C_final{i}(j,3) = elec_C{i}(j,iperm(3));
            end % permute back electrode coordinates to match the original orientation of the head
        end
        if all(sum(elec_C_final{i},2))
            rng(1,:) = max(elec_C_final{i});
            rng(2,:) = min(elec_C_final{i});
            rnge{i} = rng;
        elseif any(sum(elec_C_final{i},2)) % in case of some electrode voxels go out of image boundary then elec_C_final has 0
            elec_C_temp = elec_C_final{i}(sum(elec_C_final{i},2)>0,:);
            rng(1,:) = max(elec_C_temp);
            rng(2,:) = min(elec_C_temp);
            rnge{i} = rng;
        else rnge{i} = []; % this is the case that the electrode completely goes out of image boundary, thus this electrode actually does not exist
            error(['Electrode ' elec{i} ' goes out of image boundary! The program cannot continue. Please expand the MRI by adding empty slices on the boundaries.']);
            %                 isOut(i) = 1;
        end
    end
end
% if any(isOut)
%     error('Some of the electrodes go out of image boundary, the program cannot continue. Please expand the MRI by adding empty slices on the boundaries.');
% end
for i=1:length(rnge)
    if ~isempty(rnge{i})
        rnge{i} = rnge{i}.*repmat(template.hdr.dime.pixdim(2:4),2,1);
    end
end % use NIFTI header info to convert range info into world coordinates for subsequent electrode labeling ANDY 2014-08-12
rnge_elec = rnge;
% save([dirname filesep 'rnge_elec.mat'],'rnge_elec');

gel_C_final = cell(1,size(electrode_coord,1));
for i = 1:size(electrode_coord,1)
    if doElec(i)
        gel_C_final{i} = zeros(size(gel_C{i},1),3);
        for j=1:size(gel_C{i},1)
            if gel_C{i}(j,1)>0 && gel_C{i}(j,1)<=Nx && gel_C{i}(j,2)>0 && gel_C{i}(j,2)<=Ny && gel_C{i}(j,3)>0 && gel_C{i}(j,3)<=Nz
                gel_C_final{i}(j,1) = gel_C{i}(j,iperm(1));
                gel_C_final{i}(j,2) = gel_C{i}(j,iperm(2));
                gel_C_final{i}(j,3) = gel_C{i}(j,iperm(3));
            end % permute back gel coordinates to match the original orientation of the head
        end
    end
end % Get ready for the generation of the range of coordinates for gel in the following codes
% ANDY 2014-03-04

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

rnge = cell(1,size(electrode_coord,1));
for i = 1:size(electrode_coord,1)
    if doElec(i)
        isGel = zeros(size(gel_C_final{i},1),1);
        for j=1:length(isGel)
            if all(gel_C_final{i}(j,:)) && volume_gel(gel_C_final{i}(j,1),gel_C_final{i}(j,2),gel_C_final{i}(j,3))==1
                isGel(j)=1;
            end
        end
        gel_C_final{i} = gel_C_final{i}(isGel==1,:); % remove those gel coordinates that go into scalp/electrode/bone
        rng = zeros(2,3);
        if all(sum(gel_C_final{i},2))
            rng(1,:) = max(gel_C_final{i});
            rng(2,:) = min(gel_C_final{i});
            rnge{i} = rng;
        elseif any(sum(gel_C_final{i},2)) % in case of some gel voxels go out of image boundary then gel_C_final has 0
            gel_C_temp = gel_C_final{i}(sum(gel_C_final{i},2)>0,:);
            rng(1,:) = max(gel_C_temp);
            rng(2,:) = min(gel_C_temp);
            rnge{i} = rng;
        else rnge{i} = []; % this is the case that the gel completely goes out of image boundary, thus this gel actually does not exist
            error(['Electrode ' elec{i} ' goes out of image boundary! The program cannot continue. Please expand the MRI by adding empty slices on the boundaries.']);
        end
    end
end
for i=1:length(rnge)
    if ~isempty(rnge{i})
        rnge{i} = rnge{i}.*repmat(template.hdr.dime.pixdim(2:4),2,1);
    end
end % use NIFTI header info to convert range info into world coordinates for subsequent gel labeling ANDY 2014-08-12
rnge_gel = rnge;
% save([dirname filesep 'rnge_gel.mat'],'rnge_gel');
% ANDY 2014-03-04
save([dirname filesep baseFilename '_' uniTag '_rnge.mat'],'rnge_elec','rnge_gel');

% disp('exporting the results...')
template.img = im2uint8(volume_gel); save_untouch_nii(template,[dirname filesep baseFilename '_' uniTag '_mask_gel.nii']);
template.img = im2uint8(volume_elec); save_untouch_nii(template,[dirname filesep baseFilename '_' uniTag '_mask_elec.nii']);
% Save the results
% disp('DONE! (electrodes and gel were exported as mask_elec.nii and mask_gel.nii)')