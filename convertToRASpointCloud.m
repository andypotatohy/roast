function [data,permOrder] = convertToRASpointCloud(mri,data)
% [data,permOrder] = convertToRASpointCloud(mri,data)
%
% Flip and permute a 3D point cloud, using info read from the MRI header.
%
% (c) Yu (Andy) Huang, Parra Lab at CCNY
% yhuang16@citymail.cuny.edu
% September 2019
% June 2023 added support for qform_code

nii = load_untouch_nii(mri);

% Below code is following doc on https://nifti.nimh.nih.gov/pub/dist/src/niftilib/nifti1.h
if nii.hdr.hist.sform_code>0 % METHOD 3 as in doc 
    M1 = [nii.hdr.hist.srow_x; nii.hdr.hist.srow_y; nii.hdr.hist.srow_z; 0 0 0 1];
elseif nii.hdr.hist.sform_code==0 && nii.hdr.hist.qform_code>0 % METHOD 2 as in doc
    b = nii.hdr.hist.quatern_b; c = nii.hdr.hist.quatern_c; d = nii.hdr.hist.quatern_d;
    a=sqrt(1-b^2-c^2-d^2);

    M1 = [(a*a+b*b-c*c-d*d)*nii.hdr.dime.pixdim(2) (2*b*c-2*a*d)*nii.hdr.dime.pixdim(3)      (2*b*d+2*a*c)*nii.hdr.dime.pixdim(4)*nii.hdr.dime.pixdim(1)      nii.hdr.hist.qoffset_x;
          (2*b*c+2*a*d)*nii.hdr.dime.pixdim(2)     (a*a+c*c-b*b-d*d)*nii.hdr.dime.pixdim(3)  (2*c*d-2*a*b)*nii.hdr.dime.pixdim(4)*nii.hdr.dime.pixdim(1)      nii.hdr.hist.qoffset_y;
          (2*b*d-2*a*c)*nii.hdr.dime.pixdim(2)     (2*c*d+2*a*b)*nii.hdr.dime.pixdim(3)      (a*a+d*d-c*c-b*b)*nii.hdr.dime.pixdim(4)*nii.hdr.dime.pixdim(1)  nii.hdr.hist.qoffset_z;
          0                                        0                                         0                                                                1];
    nii.hdr.hist.sform_code = nii.hdr.hist.qform_code;
    nii.hdr.hist.qform_code = 0;
else
    error('Unknown MRI header.');
end

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