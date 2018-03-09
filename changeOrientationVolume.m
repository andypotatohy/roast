function data = changeOrientationVolume(data,perm,flipTag)

if ndims(data)~=3
    error('Unsupported data format. Please supply only 3D volume.');
end

data = permute(data,perm); % permute

for j = 1:length(flipTag)
    if flipTag(j) < 0
        data = flipdim(data,j); % flip
    end
end