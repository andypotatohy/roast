function [scalpOut,scalpMid] = cleanScalp(scalpIn)

disp('cleaning the hair for gel injection...')
scalpIn(min(scalpIn(:,1)):max(scalpIn(:,1)),min(scalpIn(:,2)):max(scalpIn(:,2)),min(scalpIn(:,3))) = 255;
% force to close the most bottom slice, because some MRIs with
% limited FOV just cut off in the middle of the face
se = 0;
isFilled = 0;
[xtemp,ytemp,ztemp] = ind2sub(size(scalpIn),find(scalpIn));
centroid = round(mean([xtemp ytemp ztemp]));
while ~all(isFilled)
    se = se+10;
    scalpMid = imfill(imclose(scalpIn,ones(se,se,se)),'holes');
    isFilled = scalpMid(centroid(1),centroid(2),centroid(3));
    % Make sure the scalp is closed and filled completely
end

se = 0;
isOpen = 1;
while any(isOpen)
    se = se+10;
    if se>30
        scalpMid(:,:,Nz) = 0; scalpMid(:,1,:) = 0; scalpMid(:,Ny,:) = 0; scalpMid(1,:,:) = 0; scalpMid(Nx,:,:) = 0;
        % force to make the image boundaries to be "opened"
    end
    scalpOut = imopen(scalpMid,ones(se,se,se));
    imgTemp = squeeze(scalpOut(:,:,Nz)); isOpen = imgTemp(:);
    imgTemp = squeeze(scalpOut(:,1,:)); isOpen = [isOpen;imgTemp(:)];
    imgTemp = squeeze(scalpOut(1,:,:)); isOpen = [isOpen;imgTemp(:)];
    imgTemp = squeeze(scalpOut(Nx,:,:)); isOpen = [isOpen;imgTemp(:)]; % Make sure the scalp is "opened" to have a smooth outer surface
end
% Scalp clean-up: for calculation of local normal vector for each electrode