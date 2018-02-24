function [elec_allCoord,gel_allCoord] = placeAndModelElectrodes(elecLoc,elecRange,scalpCleanSurf,scalpFilled,doPlace,elecPara)

% disp('placing electrodes...')

if ~isempty(elecPara.elecType,'pad')    
    se = padDim(3);
    scalpDilated = imdilate(scalpFilled,ones(se,se,se));
    dilatedScalp = scalpDilated - scalpFilled; % scalpFilled: filled scalp
    inde = find(dilatedScalp==255);
    gel_layer = zeros(length(inde),3);
    [gel_layer(:,1),gel_layer(:,2),gel_layer(:,3)] = ind2sub(size(dilatedScalp),inde);
    % Get the layer of gel to intersect with placed pad
    
    scalpDilatedMore = imdilate(scalpDilated,ones(se,se,se));
    dilatedScalp = scalpDilatedMore - scalpDilated; % scalpFilled: filled scalp
    inde = find(dilatedScalp==255);
    elec_layer = zeros(length(inde),3);
    [elec_layer(:,1),elec_layer(:,2),elec_layer(:,3)] = ind2sub(size(dilatedScalp),inde);
    % Get the layer of electrode to intersect with placed pad    
end

% figure;hold on;plot3(scalp_surface(:,1),scalp_surface(:,2),scalp_surface(:,3),'y.');
elec_allCoord = cell(1,size(elecLoc,1)); gel_allCoord = cell(1,size(elecLoc,1));
% buffer for coordinates of each electrode and gel point
for i = 1:size(elecLoc,1)
    if doPlace(i)
 
        lcl = scalpCleanSurf(elecRange(:,i),:); % local scalp surface for each electrode
        %     plot3(lcl(:,1),lcl(:,2),lcl(:,3),'b.');
        [U,D] = eig(cov(lcl)); [~,ind] = min(diag(D));
        nv = U(:,ind)'; normal = nv/norm(nv); % Local normal for each electrode
        
        switch elecPara.elecType
            case 'pad'
                
                dimTry = 10*padDim(3); % should at least be 4*padDim(3), so that one direction can cover 2 padDim(3)
                out = elecLoc(i,:) +  dimTry*normal;
                %     electrode = out + elec_height*normal;
                %     gel_in = out - 4*gel_height*normal; % coordinates of the boundaries of gel and electrode
                if norm(center - out) < norm(center - elecLoc(i,:))
                    % try computing this "center" inside this func: gravity center of the filled scalp?
                    normal = -normal;
                    %         out = elecLoc(i,:) +  dimTry*normal;
                    %         electrode = out + elec_height*normal;
                    %         gel_in = out - 4*gel_height*normal;
                end % make sure the normal is pointing out
                
                padOriShort = cross(normal,padOri);
                padOriLong = cross(normal,padOriShort);
                
                %     startPt = elecLoc(i,:) + normal * dimTry/4;
                den = 5;
                pad_coor = drawCuboid(elecLoc(i,:),[padDim(1:2) dimTry],padOriLong,padOriShort,normal,den);
                
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
                %     plot3(elec_coor(:,1),elec_coor(:,2),elec_coor(:,3),'.b');
                %     plot3(gel_coor(:,1),gel_coor(:,2),gel_coor(:,3),'.m');
                
                gel_allCoord{i} = gel_coor; elec_allCoord{i} = elec_coor; % buffer for coordinates of each electrode and gel point
                %     fprintf('%d out of %d electrodes placed...\n',i,size(elecLoc,1));
                
            case 'disc'
                
                gel_out = elecLoc(i,:) +  2*gel_height*normal;
                electrode = gel_out + elec_height*normal;
                gel_in = gel_out - 4*gel_height*normal; % coordinates of the boundaries of gel and electrode
                if norm(center - gel_out) < norm(center - elecLoc(i,:))
                    normal = -normal;
                    gel_out = elecLoc(i,:) +  2*gel_height*normal;
                    electrode = gel_out + elec_height*normal;
                    gel_in = gel_out - 4*gel_height*normal;
                end % make sure the normal is pointing out
                
                NOP = 500; verSamp = 10;
                r = 0.05:0.05:elec_radius; % parameters used for modeling of electrodes and gel
                
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
                
                gel_allCoord{i} = gel_coor; elec_allCoord{i} = elec_coor; % buffer for coordinates of each electrode and gel point
                %     fprintf('%d out of %d electrodes placed...\n',i,size(elecLoc,1));
                
            case 'ring'
                
                gel_out = elecLoc(i,:) +  2*gel_height*normal;
                electrode = gel_out + elec_height*normal;
                gel_in = gel_out - 4*gel_height*normal; % coordinates of the boundaries of gel and electrode
                if norm(center - gel_out) < norm(center - elecLoc(i,:))
                    normal = -normal;
                    gel_out = elecLoc(i,:) +  2*gel_height*normal;
                    electrode = gel_out + elec_height*normal;
                    gel_in = gel_out - 4*gel_height*normal;
                end % make sure the normal is pointing out
                
                NOP = 500; verSamp = 10;
                r = radiusInner:0.05:radiusOutter; % parameters used for modeling of electrodes and gel
                
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
                
                gel_allCoord{i} = gel_coor; elec_allCoord{i} = elec_coor; % buffer for coordinates of each electrode and gel point
                %     fprintf('%d out of %d electrodes placed...\n',i,size(elecLoc,1));
                
        end
    end
end
% xlabel('x');ylabel('y');zlabel('z'); view([270 0]);
% hold off; % Place electrodes and visualize the results