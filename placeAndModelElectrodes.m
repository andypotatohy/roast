function [elec_allCoord,gel_allCoord] = placeAndModelElectrodes(elecLoc,elecRange,scalpCleanSurf,scalpFilled,elecPara)

% disp('placing electrodes...')

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

doPad = 0;
for i=1:length(elecPara)
    if any(strcmp(elecPara(i).elecType,{'pad','Pad','PAD'}))
        doPad = 1;
    end
end

if doPad % do this one time or for each electrode?
    pad_height = elecSize(3);
    scalpDilated = imdilate(scalpFilled,ones(pad_height,pad_height,pad_height));
    dilatedScalp = scalpDilated - scalpFilled; % scalpFilled: filled scalp
    inde = find(dilatedScalp==255);
    gel_layer = zeros(length(inde),3);
    [gel_layer(:,1),gel_layer(:,2),gel_layer(:,3)] = ind2sub(size(dilatedScalp),inde);
    % Get the layer of gel to intersect with placed pad
    
    scalpDilatedMore = imdilate(scalpDilated,ones(pad_height,pad_height,pad_height));
    dilatedScalp = scalpDilatedMore - scalpDilated; % scalpFilled: filled scalp
    inde = find(dilatedScalp==255);
    elec_layer = zeros(length(inde),3);
    [elec_layer(:,1),elec_layer(:,2),elec_layer(:,3)] = ind2sub(size(dilatedScalp),inde);
    % Get the layer of electrode to intersect with placed pad
end

[Nx, Ny, Nz] = size(scalpFilled); % size of head in RAS orientation
scalpFilled(:,:,1) = 0; scalpFilled(:,:,Nz) = 0; scalpFilled(:,1,:) = 0; scalpFilled(:,Ny,:) = 0; scalpFilled(1,:,:) = 0; scalpFilled(Nx,:,:) = 0;

% DIST = zeros(size(elecLoc,1),2);
figure;hold on;plot3(scalpCleanSurf(:,1),scalpCleanSurf(:,2),scalpCleanSurf(:,3),'y.');
elec_allCoord = cell(size(elecLoc,1),1); gel_allCoord = cell(size(elecLoc,1),1);
% buffer for coordinates of each electrode and gel point
for i = 1:size(elecLoc,1)
    %     if doPlace(i)
    
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
    
    switch elecType
        case {'pad','Pad','PAD'}
            
            pad_length = elecSize(1);
            pad_width = elecSize(2);
            pad_height = elecSize(3);
            
            dimTry = 10*pad_height; % should at least be 4*pad_height, so that one direction can cover 2 pad_height
            %                 pad_out = elecLoc(i,:) +  dimTry*normal;
            %                 pad_in = elecLoc(i,:) -  dimTry*normal;
            % coordinates of the boundaries of pad
            %                 if norm(pad_out-localCentroid) <= norm(pad_in-localCentroid)
            %                     normal = -normal;
            %                 end % make sure the normal is pointing out
            
            if ischar(elecOri)
                switch elecOri
                    case 'lr'
                        elecOri = [1 0 0];
                    case 'ap'
                        elecOri = [0 1 0];
                    case 'si'
                        elecOri = [0 0 1];
                end
            end
            
            padOriShort = cross(normal,elecOri);
            padOriLong = cross(normal,padOriShort);
            
            %     startPt = elecLoc(i,:) + normal * dimTry/4;
            den = 5;
            pad_coor = drawCuboid(elecLoc(i,:),[pad_length pad_width dimTry],padOriLong,padOriShort,normal,den);
            
            %     gel_X = zeros(length(r)*verSamp*4,NOP); gel_Y = zeros(length(r)*verSamp*4,NOP); gel_Z = zeros(length(r)*verSamp*4,NOP);
            %     elec_X = zeros(length(r)*verSamp,NOP); elec_Y = zeros(length(r)*verSamp,NOP); elec_Z = zeros(length(r)*verSamp,NOP);
            %     for j = 1:length(r)
            %         [gel_X(((j-1)*verSamp*4+1):verSamp*4*j,:), gel_Y(((j-1)*verSamp*4+1):verSamp*4*j,:), gel_Z(((j-1)*verSamp*4+1):verSamp*4*j,:)] = cylinder2P(ones(verSamp*4)*r(j),NOP,gel_in,gel_out);
            %         [elec_X(((j-1)*verSamp+1):verSamp*j,:), elec_Y(((j-1)*verSamp+1):verSamp*j,:), elec_Z(((j-1)*verSamp+1):verSamp*j,:)] = cylinder2P(ones(verSamp)*r(j),NOP,gel_out,electrode);
            %     end % Use cylinders to model electrodes and gel, and calculate the coordinates of the points that make up the cylinder
            %
            %     gel_coor = floor([gel_X(:) gel_Y(:) gel_Z(:)]);
            pad_coor = round(pad_coor);
            pad_coor = unique(pad_coor,'rows');
            %     elec_coor = floor([elec_X(:) elec_Y(:) elec_Z(:)]);
            %     elec_coor = unique(elec_coor,'rows'); % clean-up of the coordinates
            
            %     plot3(pad_coor(:,1),pad_coor(:,2),pad_coor(:,3),'.m');
            
            gel_coor = intersect(pad_coor,gel_layer,'rows');
            elec_coor = intersect(pad_coor,elec_layer,'rows');
            plot3(elec_coor(:,1),elec_coor(:,2),elec_coor(:,3),'.b');
            plot3(gel_coor(:,1),gel_coor(:,2),gel_coor(:,3),'.m');
            
            gel_allCoord{i} = gel_coor; elec_allCoord{i} = elec_coor; % buffer for coordinates of each electrode and gel point
            %     fprintf('%d out of %d electrodes placed...\n',i,size(elecLoc,1));
            
        case {'disc','Disc','DISC'}
            
            disc_radius = elecSize(1);
            disc_height = elecSize(2);
            
            gel_out = elecLoc(i,:) +  2*disc_height*normal;
            electrode = gel_out + disc_height*normal;
            gel_in = gel_out - 4*disc_height*normal; % coordinates of the boundaries of gel and electrode
            
            % %               if norm(gel_out-localCentroid) <= norm(gel_in-localCentroid)
            %                 if all(min([testPointIn;testPointOut])>0) && all(max([testPointIn;testPointOut])<=size(scalpFilled)) && ...
            %                         scalpFilled(testPointIn(1),testPointIn(2),testPointIn(3))==0
            %                     normal = -normal;
            %                     gel_out = elecLoc(i,:) +  2*disc_height*normal;
            %                     electrode = gel_out + disc_height*normal;
            %                     gel_in = gel_out - 4*disc_height*normal;
            %                 end % make sure the normal is pointing out
            
            %                 DIST(i,:) = [norm(gel_out-localCentroid) norm(gel_in-localCentroid)];
            
            NOP = 500; verSamp = 10;
            r = 0.05:0.05:disc_radius; % parameters used for modeling of electrodes and gel
            
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
            
            plot3(elec_coor(:,1),elec_coor(:,2),elec_coor(:,3),'.b');
            plot3(gel_coor(:,1),gel_coor(:,2),gel_coor(:,3),'.m');
            
            gel_allCoord{i} = gel_coor; elec_allCoord{i} = elec_coor; % buffer for coordinates of each electrode and gel point
            fprintf('%d out of %d electrodes placed...\n',i,size(elecLoc,1));
            
        case {'ring','Ring','RING'}
            
            ring_radiusIn = elecSize(1);
            ring_radiusOut = elecSize(2);
            ring_height = elecSize(3);
            
            gel_out = elecLoc(i,:) +  2*ring_height*normal;
            electrode = gel_out + ring_height*normal;
            gel_in = gel_out - 4*ring_height*normal; % coordinates of the boundaries of gel and electrode
            %                 if norm(gel_out-localCentroid) <= norm(gel_in-localCentroid)
            %                     normal = -normal;
            %                     gel_out = elecLoc(i,:) +  2*ring_height*normal;
            %                     electrode = gel_out + ring_height*normal;
            %                     gel_in = gel_out - 4*ring_height*normal;
            %                 end % make sure the normal is pointing out
            
            NOP = 500; verSamp = 10;
            r = ring_radiusIn:0.05:ring_radiusOut; % parameters used for modeling of electrodes and gel
            
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
            
            plot3(elec_coor(:,1),elec_coor(:,2),elec_coor(:,3),'.b');
            plot3(gel_coor(:,1),gel_coor(:,2),gel_coor(:,3),'.m');
            
            gel_allCoord{i} = gel_coor; elec_allCoord{i} = elec_coor; % buffer for coordinates of each electrode and gel point
            %     fprintf('%d out of %d electrodes placed...\n',i,size(elecLoc,1));
            
    end
    %     end
end
xlabel('x');ylabel('y');zlabel('z'); view([270 0]);
hold off; % Place electrodes and visualize the results