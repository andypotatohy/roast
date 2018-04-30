function data = changeOrientationPointCloud(data,perm,flipTag,coordRange)
% data = changeOrientationPointCloud(data,perm,flipTag,coordRange)
% 
% Permute and flip a 3D point cloud, using permutation info in "perm" and
% flipping info in "flipTag".
% 
% (c) Yu (Andy) Huang, Parra Lab at CCNY
% yhuang16@citymail.cuny.edu
% April 2018

if ~(ndims(data)==2 && ~isempty(find(size(data)==3)))
    error('Unsupported data format. Please supply only 3D point cloud.');
end

if find(size(data)==3)==1
    data = data'; % make sure for point cloud, each row in data is for a point in the 3D space
end

for i=1:size(data,1)
    temp = data(i,:);
    data(i,:) = temp(perm);
end % permute

for j = 1:length(flipTag)
    if flipTag(j) < 0
        for i=1:size(data,1)
            temp = data(i,:);
            data(i,j) = coordRange(j) - temp(j) + 1; % Note the +1
        end
    end
end % flip