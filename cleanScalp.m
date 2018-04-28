function [scalpOut,scalpMid] = cleanScalp(scalpIn,scalpSurf)
% [scalpOut,scalpMid] = cleanScalp(scalpIn,scalpSurf)
% 
% Morphological operations on the scalp mask to make it ready to place
% electrodes.
% 
% (c) Yu (Andy) Huang, Parra Lab at CCNY
% yhuang16@citymail.cuny.edu
% April 2018

disp('cleaning the hair for gel injection...')
scalpIn(min(scalpSurf(:,1)):max(scalpSurf(:,1)),min(scalpSurf(:,2)):max(scalpSurf(:,2)),min(scalpSurf(:,3))) = 255;
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

[Nx, Ny, Nz] = size(scalpIn); % size of head in RAS orientation
se = 0;
isOpen = 1;
while any(isOpen)
    se = se+10;
    if se>30
        scalpMid(:,:,1) = 0; scalpMid(:,:,Nz) = 0; scalpMid(:,1,:) = 0; scalpMid(:,Ny,:) = 0; scalpMid(1,:,:) = 0; scalpMid(Nx,:,:) = 0;
        % force to make the image boundaries to be "opened"
    end
    scalpOut = imopen(scalpMid,ones(se,se,se));
    imgTemp = squeeze(scalpOut(:,:,Nz)); isOpen = imgTemp(:);
    imgTemp = squeeze(scalpOut(:,1,:)); isOpen = [isOpen;imgTemp(:)];
    imgTemp = squeeze(scalpOut(1,:,:)); isOpen = [isOpen;imgTemp(:)];
    imgTemp = squeeze(scalpOut(Nx,:,:)); isOpen = [isOpen;imgTemp(:)]; % Make sure the scalp is "opened" to have a smooth outer surface
end
% Scalp clean-up: for calculation of local normal vector for each electrode

% scalpMid(:,:,1) = 0; scalpMid(:,:,Nz) = 0; scalpMid(:,1,:) = 0; scalpMid(:,Ny,:) = 0; scalpMid(1,:,:) = 0; scalpMid(Nx,:,:) = 0;
% scalpOut=zeros(size(scalpMid)); smt_fil = fspecial('gaussian', 5, 2);
% for i = 1:size(scalpMid,3)
% scalpOut(:,:,i) = imfilter(scalpMid(:,:,i),smt_fil);
% end
% scalpOut = uint8(scalpOut>=round(255/2))*255;