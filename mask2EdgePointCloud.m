function [edgePtCloud,maskProc] = mask2EdgePointCloud(mask,method,se)
% [edgePtCloud,maskProc] = mask2EdgePointCloud(mask,method,se)
%
% Detect the edge of a 3D mask and convert the edge into point cloud.
% 
% (c) Yu (Andy) Huang, Parra Lab at CCNY
% yhuang16@citymail.cuny.edu
% April 2018

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
% Obtain mask surface points as a point cloud