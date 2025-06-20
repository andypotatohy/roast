function coord = drawCuboid(center,dim,longAxis,shortAxis,normalAxis,density)
% coord = drawCuboid(center,dim,longAxis,shortAxis,normalAxis,density)
% 
% Generate a point cloud of a cuboid.
% 
% (c) Yu (Andy) Huang, Parra Lab at CCNY
% yhuang16@citymail.cuny.edu
% April 2018

% dim = dim/2; % [length width thickness]

corner = center - longAxis*dim(1)/2 - shortAxis*dim(2)/2 - normalAxis*dim(3)/2;

numSamp_shortAxis = round(dim(2)*density+1);
numSamp_normalAxis = round(dim(3)*density+1);

coord = cell(numSamp_shortAxis*numSamp_normalAxis,1);
i=1;
for s = 1:numSamp_shortAxis
    for n = 1:numSamp_normalAxis
        start = corner + 1/density*shortAxis*(s-1) + 1/density*normalAxis*(n-1);
        coord{i} = drawLine(start,longAxis,dim(1),density);
        i = i+1;
    end
end

coord = cell2mat(coord);