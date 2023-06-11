function [mriRAS,isNonRAS] = convertToRAS(mri)
% [mriRAS,isNonRAS] = convertToRAS(mri)
%
% Check if input T1 MRI is in non-RAS orientation, if yes, re-orient it to
% RAS.
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

[dirname,baseFilename,ext] = fileparts(mri);
if any(permOrder~=[1 2 3]) || any(flipTag<0)
    
    isNonRAS = 1;
    mriRAS = [dirname filesep baseFilename '_ras' ext];
    
    if ~exist(mriRAS,'file')
        
        warning(['Input MRI ' mri ' is not in RAS orientation. ROAST will re-orient it into RAS now...']);
        
        img = nii.img;
        
        siz = nii.hdr.dime.dim(2:4);
        
        resolution = nii.hdr.dime.pixdim(2:4);
        
        origin = round(inv(M1)*[0;0;0;1]+1);
        origin = origin(1:3); % voxel coordinates of the origin
        
        for i=1:3
            if flipTag(i)<0
                img = flipdim(img,i); % flipping the volume
                origin(i) = siz(i) - origin(i) + 1; % flipping the origin voxel coordinates
                M_orient(oriInd(i),i) = abs(M_orient(oriInd(i),i)); % update header
            end
        end
        
        img = permute(img,permOrder); % permute the volume
        nii.img = img;
        
        origin = origin(permOrder); % permute the origin coordinates
        
        nii.hdr.dime.dim(2:4) = siz(permOrder); % update header
        
        nii.hdr.dime.pixdim(1) = abs(nii.hdr.dime.pixdim(1)); % update header
        nii.hdr.dime.pixdim(2:4) = resolution(permOrder); % update header
        
        M_orient = M_orient(:,permOrder);
        nii.hdr.hist.srow_x(1:3) = M_orient(1,:);
        nii.hdr.hist.srow_y(1:3) = M_orient(2,:);
        nii.hdr.hist.srow_z(1:3) = M_orient(3,:); % update header
        
        nii.hdr.hist.srow_x(4) = -origin(1)*nii.hdr.hist.srow_x(1)-origin(2)*nii.hdr.hist.srow_x(2)-origin(3)*nii.hdr.hist.srow_x(3);
        nii.hdr.hist.srow_y(4) = -origin(1)*nii.hdr.hist.srow_y(1)-origin(2)*nii.hdr.hist.srow_y(2)-origin(3)*nii.hdr.hist.srow_y(3);
        nii.hdr.hist.srow_z(4) = -origin(1)*nii.hdr.hist.srow_z(1)-origin(2)*nii.hdr.hist.srow_z(2)-origin(3)*nii.hdr.hist.srow_z(3);
        nii.hdr.hist.qoffset_x = nii.hdr.hist.srow_x(4);
        nii.hdr.hist.qoffset_y = nii.hdr.hist.srow_y(4);
        nii.hdr.hist.qoffset_z = nii.hdr.hist.srow_z(4); % update header
        
        save_untouch_nii(nii,mriRAS);
        
        disp([mri ' is now in RAS orientation and saved as:']);
        disp(mriRAS);
        disp('It''ll be used as the input for ROAST.');
        
    else
        
        warning(['Input MRI ' mri ' is not in RAS orientation, and has been re-oriented into RAS and saved as ' mriRAS '. ROAST will use that file as the input.']);
        
    end
    
else
    
    isNonRAS = 0;
    mriRAS = mri;
    
end