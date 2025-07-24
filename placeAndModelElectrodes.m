function [elec_allCoord,gel_allCoord] = placeAndModelElectrodes(elecLoc,elecRange,scalpCleanSurf,scalpFilled,elecPlacing,elecPara,res,isDebug,uniTag)
% [elec_allCoord,gel_allCoord] = placeAndModelElectrodes(elecLoc,elecRange,scalpCleanSurf,scalpFilled,elecPlacing,elecPara,res,isDebug,uniTag)
% 
% Place and generate the point cloud for each placed electrode and gel.
% 
% (c) Yu (Andy) Huang, Parra Lab at CCNY
% yhuang16@citymail.cuny.edu
% April 2018

disp('placing electrodes...')

if length(elecPara)==1
    numOfElec = size(elecLoc,1);
    temp = repmat(elecPara,numOfElec,1);
    if size(elecPara.elecSize,1)>1
        for i=1:length(temp), temp(i).elecSize = temp(i).elecSize(i,:); end
    end
    if size(elecPara.elecOri,1)>1
        for i=1:length(temp), temp(i).elecOri = temp(i).elecOri(i,:); end
    end
    elecPara = temp;
end

padH = zeros(length(elecPara),1);
for i=1:length(elecPara)
    if strcmpi(elecPara(i).elecType,'pad')
        padH(i) = elecPara(i).elecSize(3)/res;
    end
end
indPad = find(padH>0);
if ~isempty(indPad)
    ind2allPH = zeros(length(elecPara),1);
    [allPH,~,ind2allPH(indPad)] = unique(padH(indPad));
    gel_layer = cell(length(allPH),1);
    elec_layer = cell(length(allPH),1);
    for i=1:length(allPH)
        strel = ones(round(allPH(i)),round(allPH(i)),round(allPH(i))); % this is not perfect yet
        [gel_layer{i},scalpDilated] = mask2EdgePointCloud(scalpFilled,'dilate',strel);
        elec_layer{i} = mask2EdgePointCloud(scalpDilated,'dilate',strel);
        % Get the layer of electrode/gel to intersect with placed pad
    end
end

[Nx, Ny, Nz] = size(scalpFilled); % size of head in RAS orientation
scalpFilled(:,:,1) = 0; scalpFilled(:,:,Nz) = 0; scalpFilled(:,1,:) = 0; scalpFilled(:,Ny,:) = 0; scalpFilled(1,:,:) = 0; scalpFilled(Nx,:,:) = 0;

if isDebug
    figName = ['Electrode placement in Simulation: ' uniTag];
    figure('Name',[figName '. Move your mouse to rotate.'],'NumberTitle','off');
    set(gcf,'color','w');
    plot3(scalpCleanSurf(:,1),scalpCleanSurf(:,2),scalpCleanSurf(:,3),'y.');
    hold on;
end
elec_allCoord = cell(size(elecLoc,1),1); gel_allCoord = cell(size(elecLoc,1),1);
% buffer for coordinates of each electrode and gel point
for i = 1:length(elecPara) % size(elecLoc,1)    
    lcl = scalpCleanSurf(elecRange(i,:),:); % local scalp surface for each electrode
    %     plot3(lcl(:,1),lcl(:,2),lcl(:,3),'b.');
    %         localCentroid = mean(lcl);
    [U,D] = eig(cov(lcl)); [~,ind] = min(diag(D));
    nv = U(:,ind)'; normal = nv/norm(nv); % Local normal for each electrode
    
    lenTry=1;
    testPointIn = round(elecLoc(i,:) - lenTry*normal);
    testPointOut = round(elecLoc(i,:) + lenTry*normal);
    while all(min([testPointIn;testPointOut])>0) && all(max([testPointIn;testPointOut])<=size(scalpFilled)) && ...
            ~xor(scalpFilled(testPointIn(1),testPointIn(2),testPointIn(3)),scalpFilled(testPointOut(1),testPointOut(2),testPointOut(3)))
        lenTry = lenTry+1;
        testPointIn = round(elecLoc(i,:) - lenTry*normal);
        testPointOut = round(elecLoc(i,:) + lenTry*normal);
    end
    if all(min([testPointIn;testPointOut])>0) && all(max([testPointIn;testPointOut])<=size(scalpFilled)) && ...
            scalpFilled(testPointIn(1),testPointIn(2),testPointIn(3))==0
        normal = -normal;
    end % make sure the normal is pointing out
    
    switch lower(elecPara(i).elecType)
        case 'pad'
            if ~isempty(strfind(lower(elecPlacing{i}),'custom')) || ismember(lower(elecPlacing{i}),{'nk1';'nk2';'nk3';'nk4'})
                fprintf('placing electrode %s (%d out of %d)...\n',elecPlacing{i},i,size(elecLoc,1));
            else
                fprintf('placing electrode %s (%d out of %d) in the %s layout...\n',elecPlacing{i},i,size(elecLoc,1),lower(elecPara(i).capType));
            end
            
            pad_length = elecPara(i).elecSize(1)/res;
            pad_width = elecPara(i).elecSize(2)/res;
            
            dimTry = mean([pad_length pad_width]);
            % bigger electrode needs bigger dimTry (needs to establish a better relation)
            
            if ischar(elecPara(i).elecOri)
                switch lower(elecPara(i).elecOri)
                    case 'lr'
                        elecPara(i).elecOri = [1 0 0];
                    case 'ap'
                        elecPara(i).elecOri = [0 1 0];
                    case 'si'
                        elecPara(i).elecOri = [0 0 1];
                end
            end
            
            padOriShort = cross(normal,elecPara(i).elecOri); padOriShort = padOriShort/norm(padOriShort);
            padOriLong = cross(normal,padOriShort); padOriLong = padOriLong/norm(padOriLong);
            
            %     startPt = elecLoc(i,:) + normal * dimTry/4;
            den = 2; % 2 points per pixel
            pad_coor = drawCuboid(elecLoc(i,:),[pad_length pad_width dimTry],padOriLong,padOriShort,normal,den);
            
            pad_coor = unique(round(pad_coor),'rows'); % clean-up of the coordinates
            
            %     plot3(pad_coor(:,1),pad_coor(:,2),pad_coor(:,3),'.m');
            
            gel_coor = intersect(pad_coor,gel_layer{ind2allPH(i)},'rows');
            elec_coor = intersect(pad_coor,elec_layer{ind2allPH(i)},'rows');
            
            if isDebug
                plot3(elec_coor(:,1),elec_coor(:,2),elec_coor(:,3),'.b');
                plot3(gel_coor(:,1),gel_coor(:,2),gel_coor(:,3),'.m');
            end
            
            gel_allCoord{i} = gel_coor; elec_allCoord{i} = elec_coor; % buffer for coordinates of each electrode and gel point
            
        case 'disc'
            if ~isempty(strfind(lower(elecPlacing{i}),'custom')) || ismember(lower(elecPlacing{i}),{'nk1';'nk2';'nk3';'nk4'})
                fprintf('placing electrode %s (%d out of %d)...\n',elecPlacing{i},i,size(elecLoc,1));
            else
                fprintf('placing electrode %s (%d out of %d) in the %s layout...\n',elecPlacing{i},i,size(elecLoc,1),lower(elecPara(i).capType));
            end
            
            disc_radius = elecPara(i).elecSize(1)/res;
            disc_height = elecPara(i).elecSize(2)/res;
            
            dimTry = disc_radius;
            
            gel_out = elecLoc(i,:) +  2*disc_height*normal;
            electrode = gel_out + disc_height*normal;
            gel_in = gel_out - dimTry*normal; % coordinates of the boundaries of gel and electrode
            
            den = 2; % 2 points per pixel
            gel_coor = drawCylinder(0,disc_radius,gel_in,gel_out,den);
            elec_coor = drawCylinder(0,disc_radius,gel_out,electrode,den);
            % Use cylinders to model electrodes and gel, and calculate the coordinates of the points that make up the cylinder
            
            gel_coor = unique(round(gel_coor),'rows');
            elec_coor = unique(round(elec_coor),'rows'); % clean-up of the coordinates
            
            if isDebug
                plot3(elec_coor(:,1),elec_coor(:,2),elec_coor(:,3),'.b');
                plot3(gel_coor(:,1),gel_coor(:,2),gel_coor(:,3),'.m');
            end
            
            gel_allCoord{i} = gel_coor; elec_allCoord{i} = elec_coor; % buffer for coordinates of each electrode and gel point
            
        case 'ring'
            if ~isempty(strfind(lower(elecPlacing{i}),'custom')) || ismember(lower(elecPlacing{i}),{'nk1';'nk2';'nk3';'nk4'})
                fprintf('placing electrode %s (%d out of %d)...\n',elecPlacing{i},i,size(elecLoc,1));
            else
                fprintf('placing electrode %s (%d out of %d) in the %s layout...\n',elecPlacing{i},i,size(elecLoc,1),lower(elecPara(i).capType));
            end
            
            ring_radiusIn = elecPara(i).elecSize(1)/res;
            ring_radiusOut = elecPara(i).elecSize(2)/res;
            ring_height = elecPara(i).elecSize(3)/res;
            
            dimTry = mean([ring_radiusOut ring_radiusIn]);
            
            gel_out = elecLoc(i,:) +  2*ring_height*normal;
            electrode = gel_out + ring_height*normal;
            gel_in = gel_out - dimTry*normal; % coordinates of the boundaries of gel and electrode
            
            den = 2; % 2 points per pixel
            gel_coor = drawCylinder(ring_radiusIn,ring_radiusOut,gel_in,gel_out,den);
            elec_coor = drawCylinder(ring_radiusIn,ring_radiusOut,gel_out,electrode,den);
            % Use cylinders to model electrodes and gel, and calculate the coordinates of the points that make up the cylinder
            
            gel_coor = unique(round(gel_coor),'rows');
            elec_coor = unique(round(elec_coor),'rows'); % clean-up of the coordinates
            
            if isDebug
                plot3(elec_coor(:,1),elec_coor(:,2),elec_coor(:,3),'.b');
                plot3(gel_coor(:,1),gel_coor(:,2),gel_coor(:,3),'.m');
            end
            
            gel_allCoord{i} = gel_coor; elec_allCoord{i} = elec_coor; % buffer for coordinates of each electrode and gel point
    end
end
if isDebug
    xlabel('x');ylabel('y');zlabel('z'); axis equal;
    hold off; rotate3d on; % Place electrodes and visualize the results
end