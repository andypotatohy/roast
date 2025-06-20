function coord = drawCylinder(innerRadius,outterRadius,top,bottom,density)
% coord = drawCylinder(innerRadius,outterRadius,top,bottom,density)
%
% Generate a point cloud of a cylinder.
%
% (c) Yu (Andy) Huang, Parra Lab at CCNY
% yhuang16@citymail.cuny.edu
% May 2018

r = (innerRadius+1/density):1/density:outterRadius;
cylinderHeight = norm(top-bottom);

coord = cell(length(r),1);
for i = 1:length(r)
    [x,y,z] = cylinder2P(r(i)*ones(round(cylinderHeight*density)),round(2*pi*r(i)*density),top,bottom);
    coord{i} = single([x(:) y(:) z(:)]);
end

coord = cell2mat(coord);