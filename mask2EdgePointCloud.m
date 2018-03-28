function [edgePtCloud,maskProc] = mask2EdgePointCloud(mask,method,se)

switch method
    case 'erode'
        maskProc = imerode(mask,se);
        mask_edge = mask-maskProc; % Get the edge of mask
    case 'dilate'
        maskProc = imdilate(mask,se);
        mask_edge = maskProc-mask; % Get the edge of mask
    otherwise
        error('Please enter either erode or dilate for the Method.');
end

inde = find(mask_edge==max(mask_edge(:)));

edgePtCloud = zeros(length(inde),3);
[edgePtCloud(:,1),edgePtCloud(:,2),edgePtCloud(:,3)] = ind2sub(size(mask_edge),inde);
% Obtain mask surface points as a matrix