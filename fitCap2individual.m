function electrode_coord = fitCap2individual(scalp,scalp_surface,landmarks,capInfo,isBiosemi,elecNeeded)

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
end
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

s = right-left; s = s/norm(s);
c = nasion-inion; c = c/norm(c);
a = cross(s,c); a = a/norm(a); % vectors to be used in affine transform

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
