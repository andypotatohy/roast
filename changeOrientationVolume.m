function data = changeOrientationVolume(data,perm,flipTag)
% data = changeOrientationVolume(data,perm,flipTag)
% 
% Permute and flip a 3D matrix, using permutation info in "perm" and
% flipping info in "flipTag".
% 
% (c) Yu (Andy) Huang, Parra Lab at CCNY
% yhuang16@citymail.cuny.edu
% April 2018

if ndims(data)~=3
    error('Unsupported data format. Please supply only 3D volume.');
end

data = permute(data,perm); % permute

for j = 1:length(flipTag)
    if flipTag(j) < 0
        data = flipdim(data,j); % flip
    end
end