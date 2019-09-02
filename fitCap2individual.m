function [electrode_coord,center]= fitCap2individual(scalp,scalp_surface,landmarks,P2,capInfo,indNeed,isBiosemi,isEGI)
% [electrode_coord,center]= fitCap2individual(scalp,scalp_surface,landmarks,P2,capInfo,indNeed,isBiosemi,isEGI)
%
% Place the electrodes with pre-defined coordinates in the standard EEG
% system (e.g., 10/05, BioSemi, or EGI system).
% 
% Landmarks follow the order of: nasion, inion, right, left, front neck,
% and back neck.
%
% Note for the EGI layout, just directly project the coordinates onto the
% scalp surface without adjusting.
% 
% (c) Yu (Andy) Huang, Parra Lab at CCNY
% yhuang16@citymail.cuny.edu
% April 2018
% August 2019 added EGI layout

nasion = landmarks(1,:);
inion = landmarks(2,:);
right = landmarks(3,:);
left = landmarks(4,:);
% front_neck = landmarks(5,:);
% back_neck = landmarks(6,:);

disp('measuring head size...')
L = norm(inion-nasion); % Distance between nasion and inion
line_center = (inion+nasion)/2; % Midpoint between nasion and inion

if ~isEGI
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
end

disp('wearing the cap...')
elec = capInfo{1};
if ~isEGI
    if ~isBiosemi
        centralElec = {'Oz';'POz';'Pz';'CPz';'Cz';'FCz';'Fz';'AFz';'Fpz'};
    else
        centralElec = {'A19';'POzAid';'A6';'CPzAid';'A1';'FCzAid';'E17';'AFzAid';'E12'};
    end
    [~,indCentralElec] = ismember(centralElec,elec);
else
    indCentralElec = [];
end

indFit = cat(1,indCentralElec,indNeed); % only fit those elec specified by users (to save time)
elec_template = cell2mat(capInfo(2:4));
elec_template = elec_template(indFit,:);
if isempty(strfind(P2,'example/nyhead'))
    data = load_untouch_nii(P2);
    elec_template = elec_template./repmat(data.hdr.dime.pixdim(2:4),length(indFit),1);
    % account for MRI resolution (so can do non-1mm, anisotropic MRI accurately)
else
    elec_template = elec_template./repmat([0.5 0.5 0.5],length(indFit),1);
end

theta = 23;
alpha = ((360-10*theta)/2)*(pi/180);
h = (L/2)*(1/tan(alpha));
% For the calculation of the center of electrode coordinates
% This applies to 1010 and BioSemi systems, but only roughly applies to EGI system

s = right-left; s = s/norm(s);
c = nasion-inion; c = c/norm(c);
a = cross(s,c); a = a/norm(a); % vectors to be used in affine transform

if ~isEGI
    disp('adjust the cap for optimized position...this will take a while...')
    factor = 1:-0.05:0.5; % Adjusting factor
    F = zeros(length(factor),1);
else
    factor = 1;
end
CENTER = zeros(length(factor),3);
ELEC_COORD = zeros(size(elec_template,1),3,length(factor));
for n = 1:length(factor)
    if ~isEGI, fprintf('Iteration No. %d...\n', n); end
    center = line_center + h*factor(n)*a; % Adjust the center
    CENTER(n,:) = center; % buffer
    scale = round(max(size(scalp))/2);
    shift = center';
    
    affine = scale * [s' c' a' shift/scale;0 0 0 1/scale]; % Affine transform matrix
    
    elec_adjusted = [elec_template';ones(1,size(elec_template,1))];
    elec_adjusted(3,:) = elec_adjusted(3,:)*factor(n);
    % Adjust the z-coordinate correspondingly
    elec_transformed = affine * elec_adjusted;
    elec_transformed = elec_transformed(1:3,:)';
    % Affine transform the EasyCap coordinate to an approximate position for each electrode outside of the scalp surface
    
    idx = zeros(size(elec_transformed,1),1);
    [cosineAngle,indOnScalpSurf] = project2ClosestSurfacePoints(elec_transformed,scalp_surface,center);
    for i = 1:length(idx)
%         testPts = scalp_surface(indOnScalpSurf(cosineAngle(:,i) > max(cosineAngle(:,i))*0.99993,i),:);
        testPts = scalp_surface(indOnScalpSurf(cosineAngle(:,i) > prctile(cosineAngle(:,i),99.99),i),:);
        [~,indFarthestOnTestPts] = map2Points(center,testPts,'farthest');
        idx(i) = indOnScalpSurf(indFarthestOnTestPts,i);
        % Find the only point on the outer surface of the scalp for each electrode,
        % i.e., the exact coordinates for each electrode on the scalp surface
    end
    
    elec_interp = scalp_surface(idx,:); % exact coordinates for each electrode
    ELEC_COORD(:,:,n) = elec_interp; % buffer
    
    if ~isEGI % only adjust the cap for non-EGI systems
        %     center_points = [inion;elec_interp(indCentralElec,:);nasion];
        center_points = [inion;elec_interp(1:length(indCentralElec),:);nasion];
        % coordinates for electrodes on central sagittal line
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
end

if ~isEGI, [~,index] = min(F); else, index = 1; end
% electrode_coord = ELEC_COORD(:,:,index); % exact coordinate for each electrode projected on the scalp surface
electrode_coord = ELEC_COORD(length(indCentralElec)+1:end,:,index); % exact coordinate for each electrode projected on the scalp surface
center = CENTER(index,:); % center of electrode coordinates