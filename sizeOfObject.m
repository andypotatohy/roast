function [size_descend,ind] = sizeOfObject(img, conn)

% [size_descend,ind] = sizeOfObject(img, conn)
%
% Compute the size of each connected component in binary image and sort the
% size in descending order.
% Input: img--2D or 3D binary image;
%        conn--connectivity defined by the user. It can have the following
%        value:
%                  Value     Meaning
%                   Two-dimensional connectivities
%                    4        4-connected neighborhood
%                    8        8-connected neighborhood
%                   Three-dimensional connectivities
%                    6        6-connected neighborhood
%                    18       18-connected neighborhood
%                    26       26-connected neighborhood
%         The default connectivity is 8 for 2D, and 26 for 3D.
% Output: size_descend--size of each connected component in descending
%         order
%         ind--index of sorted connected component
%
% (c) Yu Huang (Andy), May 2011
% The Neural Engineering Lab, Dept. of Biomedical Engineering, City College of New York
% Send bugs to yhuang16@citymail.cuny.edu

if nargin < 2
    if ndims(img) == 2
        conn = 8;
    elseif ndims(img) == 3
        conn = 26;
    end
end

CC = bwconncomp(img,conn);
size_obj = zeros(1,CC.NumObjects);
for k = 1:CC.NumObjects
    size_obj(k) = length(CC.PixelIdxList{k});
end
[size_descend,ind] = sort(size_obj,'descend');