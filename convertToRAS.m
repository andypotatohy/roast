function convertToRAS(mri,subject)
% convertToRAS(mri,subject)
% 
% Check if input T1 MRI is in non-RAS orientation, if yes, re-orient it to
% RAS.
% 
% (c) Yu (Andy) Huang, Parra Lab at CCNY
% yhuang16@citymail.cuny.edu
% August 2019

M1 = [mri.hdr.hist.srow_x; mri.hdr.hist.srow_y; mri.hdr.hist.srow_z; 0 0 0 1];

M_orient = M1(1:3,1:3);

[~,oriInd] = max(abs(M_orient));
[~,permOrder] = sort(oriInd); % permutation order to RAS system
flipTag = [sign(M_orient(oriInd(1),1)) sign(M_orient(oriInd(2),2)) sign(M_orient(oriInd(3),3))];
% detect if the head is flipped in each direction compared to RAS system

if any(permOrder~=[1 2 3]) || any(flipTag<0)
    
    [pth,nam,ext] = fileparts(subject);
    subjectBak = [pth filesep nam '_backup' ext];
    warning(['Input MRI ' subject ' is not in RAS orientation. ROAST will re-orient it into RAS.']);
    disp(['Re-orienting ' subject ' ...']);

    copyfile(subject,subjectBak); % backup the original
    
    img = mri.img;
    
    siz = mri.hdr.dime.dim(2:4);
    
    resolution = mri.hdr.dime.pixdim(2:4);
    
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
    mri.img = img;
    
    origin = origin(permOrder); % permute the origin coordinates
    
    mri.hdr.dime.dim(2:4) = siz(permOrder); % update header
    
    mri.hdr.dime.pixdim(1) = abs(mri.hdr.dime.pixdim(1)); % update header
    mri.hdr.dime.pixdim(2:4) = resolution(permOrder); % update header
    
    M_orient = M_orient(:,permOrder);
    mri.hdr.hist.srow_x(1:3) = M_orient(1,:);
    mri.hdr.hist.srow_y(1:3) = M_orient(2,:);
    mri.hdr.hist.srow_z(1:3) = M_orient(3,:); % update header
    
    mri.hdr.hist.srow_x(4) = -origin(1)*mri.hdr.hist.srow_x(1)-origin(2)*mri.hdr.hist.srow_x(2)-origin(3)*mri.hdr.hist.srow_x(3);
    mri.hdr.hist.srow_y(4) = -origin(1)*mri.hdr.hist.srow_y(1)-origin(2)*mri.hdr.hist.srow_y(2)-origin(3)*mri.hdr.hist.srow_y(3);
    mri.hdr.hist.srow_z(4) = -origin(1)*mri.hdr.hist.srow_z(1)-origin(2)*mri.hdr.hist.srow_z(2)-origin(3)*mri.hdr.hist.srow_z(3);
    mri.hdr.hist.qoffset_x = mri.hdr.hist.srow_x(4);
    mri.hdr.hist.qoffset_y = mri.hdr.hist.srow_y(4);
    mri.hdr.hist.qoffset_z = mri.hdr.hist.srow_z(4); % update header
    
    save_untouch_nii(mri,subject);
    
    disp([subject ' is now in RAS orientation, original MRI is backed up as ' subjectBak]);
end