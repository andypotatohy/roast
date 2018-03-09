function data = changeOrientation(data,perm,flipTag)

if ndims(data)==1 || ndims(data)>3
    error('Unsupported data format. Please supply only 3D volume or point cloud.');
end

if ndims(data)==2 && isempty(find(size(data)==3))
    error('For point cloud, one dimension must be 3');
end

if ndims(data)==2 && find(size(data)==3)==1
    data = data'; % make sure for point cloud, each row in data is for a point in the 3D space
end

switch ndims(data)
    case 2 % point cloud
        
        for i=1:size(data,1)
            temp = data(i,:);
            if ~isempty(temp), data(i,:) = temp(perm); end
        end % permute
        
        for j = 1:length(flipTag)
            if flipTag(j) < 0
                for i=1:size(data,1)
                    temp = data(i,:);
                    if ~isempty(temp)
                        data(i,j) = sizeRAS(j) - temp(j) + 1; % Note the +1
                    end
                end
            end
        end % flip
        
    case 3 % 3D volume
        
        data = permute(data,perm); % permute
        
        for j = 1:length(flipTag)
            if flipTag(j) < 0
                data = flipdim(data,j); % flip
            end
        end
end