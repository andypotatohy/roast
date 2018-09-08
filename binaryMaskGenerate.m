function [empt,mask1,mask2,mask3,mask4,mask5,mask6] = binaryMaskGenerate(data1,data2,data3,data4,data5,data6)

% [empt,mask1,mask2,mask3,mask4,mask5,mask6] = binaryMaskGenerate(data1,data2,data3,data4,data5,data6)
%
% Use max operation to generate binary mask for each tissue type from the
% results generated from SPM. It can also be used to update each mask
% after one mask is processed.
% Input: data1~data6: 6 tissue types, at least need 1 tissue type;
% Output: empt: empty mask represents those empty voxels which do not
% belong to any tissue type; mask1~mask6: 6 masks corresponding to
% data1~data6.
%
% (c) Yu Huang (Andy), May 2011
% The Neural Engineering Lab, Dept. of Biomedical Engineering, City College of New York
% Send bugs to yhuang16@citymail.cuny.edu

if nargin < 2
    data2 = [];
end
if nargin < 3
    data3 = [];
end
if nargin < 4
    data4 = [];
end
if nargin < 5
    data5 = [];
end
if nargin < 6
    data6 = [];
end

data = [];
data(:,:,:,1) = zeros(size(data1));
data(:,:,:,2) = data1;
if ~isempty(data2)
    data(:,:,:,3) = data2;
end
if ~isempty(data3)
    data(:,:,:,4) = data3;
end
if ~isempty(data4)
    data(:,:,:,5) = data4;
end
if ~isempty(data5)
    data(:,:,:,6) = data5;
end
if ~isempty(data6)
    data(:,:,:,7) = data6;
end

[~,maxind] = max(data,[],4);

empt = (maxind==1);
mask1 = (maxind==2);
if ~isempty(data2)
    mask2 = (maxind==3);
end
if ~isempty(data3)
    mask3 = (maxind==4);
end
if ~isempty(data4)
    mask4 = (maxind==5);
end
if ~isempty(data5)
    mask5 = (maxind==6);
end
if ~isempty(data6)
    mask6 = (maxind==7);
end