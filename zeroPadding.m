function zeroPadding(mri)
% zeroPadding(mri)
% 
% Adding empty slices to the MRI to avoid complication when placing
% electrodes at locations close to the image boundaries
% 
% (c) Yu (Andy) Huang, Parra Lab at CCNY
% yhuang16@citymail.cuny.edu
% April 2018

nii = load_untouch_nii(mri);
newSize = size(nii.img) + [20 20 20];

img = zeros(newSize);
img(11:end-10,11:end-10,11:end-10) = nii.img;
nii.img = img;
nii.hdr.dime.dim(2:4) = newSize;

origin = inv([nii.hdr.hist.srow_x;nii.hdr.hist.srow_y;nii.hdr.hist.srow_z;0 0 0 1])*[0;0;0;1];
origin = origin(1:3) + 10;

nii.hdr.hist.srow_x(4) = -origin(1)*nii.hdr.hist.srow_x(1)-origin(2)*nii.hdr.hist.srow_x(2)-origin(3)*nii.hdr.hist.srow_x(3);
nii.hdr.hist.srow_y(4) = -origin(1)*nii.hdr.hist.srow_y(1)-origin(2)*nii.hdr.hist.srow_y(2)-origin(3)*nii.hdr.hist.srow_y(3);
nii.hdr.hist.srow_z(4) = -origin(1)*nii.hdr.hist.srow_z(1)-origin(2)*nii.hdr.hist.srow_z(2)-origin(3)*nii.hdr.hist.srow_z(3);

nii.hdr.hist.qoffset_x = nii.hdr.hist.srow_x(4);
nii.hdr.hist.qoffset_y = nii.hdr.hist.srow_y(4);
nii.hdr.hist.qoffset_z = nii.hdr.hist.srow_z(4);

[dirname,baseFilename,ext] = fileparts(mri);
save_untouch_nii(nii,[dirname filesep baseFilename 'Padded' ext]);

disp([mri ' has been zero-padded, and is saved as:']);
disp([dirname filesep baseFilename 'Padded' ext]);
disp('Please use this file as the input for ROAST main function.');