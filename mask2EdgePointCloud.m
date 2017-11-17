function edgePtCloud = mask2EdgePointCloud(mask)

mask_edge = mask-imerode(mask,ones(3,3,3)); % Get the edge of mask

inde = find(mask_edge==max(mask_edge(:)));

edgePtCloud = zeros(length(inde),3);
[edgePtCloud(:,1),edgePtCloud(:,2),edgePtCloud(:,3)] = ind2sub(size(mask_edge),inde);
% Obtain mask surface points as a matrix