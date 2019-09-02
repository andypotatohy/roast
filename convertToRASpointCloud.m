function [data,permOrder] = convertToRASpointCloud(mri,data)
% [data,permOrder] = convertToRASpointCloud(mri,data)
%
% Flip and permute a 3D point cloud, using info read from the MRI header.
%
% (c) Yu (Andy) Huang, Parra Lab at CCNY
% yhuang16@citymail.cuny.edu
% September 2019

nii = load_untouch_nii(mri);
M1 = [nii.hdr.hist.srow_x; nii.hdr.hist.srow_y; nii.hdr.hist.srow_z; 0 0 0 1];

M_orient = M1(1:3,1:3);

[~,oriInd] = max(abs(M_orient));
[~,permOrder] = sort(oriInd); % permutation order to RAS system
flipTag = [sign(M_orient(oriInd(1),1)) sign(M_orient(oriInd(2),2)) sign(M_orient(oriInd(3),3))];
% detect if the head is flipped in each direction compared to RAS system

% flipping
for j=1:length(flipTag)
    if flipTag(j)<0
        data(:,j) = nii.hdr.dime.dim(j+1) - data(:,j) + 1; % Note the +1
    end
end

% permute
for i=1:size(data,1)
    temp = data(i,:);
    data(i,:) = temp(permOrder);
end