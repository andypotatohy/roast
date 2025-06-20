function mriPD = zeroPadding(mri,padNum)
% mriPD = zeroPadding(mri,padNum)
%
% Adding empty slices to the MRI to avoid complication when placing
% electrodes at locations close to the image boundaries
%
% (c) Yu (Andy) Huang, Parra Lab at CCNY
% yhuang16@citymail.cuny.edu
% April 2018

[dirname,baseFilename,ext] = fileparts(mri);
if isempty(dirname), dirname = pwd; end

if ~isempty(strfind(mri,'_padded'))
    warning([mri ' has already been zero-padded. Nothing will happen here. If you meant to add empty slices on the MRI, please provide MRI name without the _padded suffix.']);
    mriPD = mri;
    return;
end

if ~isempty(strfind(mri,'example/nyhead_'))
    mriPD = ['example/nyhead_padded' num2str(padNum) baseFilename(7:end) ext];
else
    mriPD = [dirname filesep baseFilename '_padded' num2str(padNum) ext];
end

if exist(mriPD,'file')
    
    warning([mri ' has already been zero-padded by ' num2str(padNum) ' empty slices in the six directions and saved as ' mriPD '. ROAST will use that file as the input.']);
    
else
    
    disp(['Padding ' mri ' by ' num2str(padNum) ' empty slices in all the six directions...']);
    nii = load_untouch_nii(mri);
    newSize = size(nii.img) + 2*[padNum padNum padNum];
    
    img = zeros(newSize);
    img(padNum+1:end-padNum,padNum+1:end-padNum,padNum+1:end-padNum) = nii.img;
    nii.img = img;
    nii.hdr.dime.dim(2:4) = newSize;
    
    origin = inv([nii.hdr.hist.srow_x;nii.hdr.hist.srow_y;nii.hdr.hist.srow_z;0 0 0 1])*[0;0;0;1];
    origin = origin(1:3) + padNum;
    
    nii.hdr.hist.srow_x(4) = -origin(1)*nii.hdr.hist.srow_x(1)-origin(2)*nii.hdr.hist.srow_x(2)-origin(3)*nii.hdr.hist.srow_x(3);
    nii.hdr.hist.srow_y(4) = -origin(1)*nii.hdr.hist.srow_y(1)-origin(2)*nii.hdr.hist.srow_y(2)-origin(3)*nii.hdr.hist.srow_y(3);
    nii.hdr.hist.srow_z(4) = -origin(1)*nii.hdr.hist.srow_z(1)-origin(2)*nii.hdr.hist.srow_z(2)-origin(3)*nii.hdr.hist.srow_z(3);
    
    nii.hdr.hist.qoffset_x = nii.hdr.hist.srow_x(4);
    nii.hdr.hist.qoffset_y = nii.hdr.hist.srow_y(4);
    nii.hdr.hist.qoffset_z = nii.hdr.hist.srow_z(4);
    
    save_untouch_nii(nii,mriPD);
    
    disp([mri ' has been zero-padded, and is saved as:']);
    disp(mriPD);
    disp('It''ll be used as the input for ROAST.');
    
end