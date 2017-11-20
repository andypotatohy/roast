function [mask,pointCloud,sizeRAS] = changeOrientationIn(permIn,flipTag,mask,pointCloud)

mask = permute(mask,permIn); % permute mask into RAS

for i=1:size(pointCloud,1)
    temp = pointCloud(i,:);
    if ~isempty(temp), pointCloud(i,:) = temp(permIn); end
end % permute points accordingly

sizeRAS = size(mask); % size of head in RAS orientation

for j = 1:length(flipTag)
    if flipTag(j) < 0
        mask = flipdim(mask,j); % flip mask in flipped direction
        
        for i=1:size(pointCloud,1)
            temp = pointCloud(i,:);
            if ~isempty(temp)
                pointCloud(i,j) = sizeRAS(j) - temp(j) + 1; % Note the +1
            end
        end
        
    end
end % flip points accordingly