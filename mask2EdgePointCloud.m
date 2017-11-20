function edgePtCloud = mask2EdgePointCloud(mask,method,se)

switch method
    case 'erode'
        mask_edge = mask-imerode(mask,se); % Get the edge of mask
    case 'dilate'
        mask_edge = imdilate(mask,se)-mask; % Get the edge of mask
    otherwise
        error('Please enter either erode or dilate for the Method.');
end

inde = find(mask_edge==max(mask_edge(:)));

edgePtCloud = zeros(length(inde),3);
[edgePtCloud(:,1),edgePtCloud(:,2),edgePtCloud(:,3)] = ind2sub(size(mask_edge),inde);
% Obtain mask surface points as a matrix