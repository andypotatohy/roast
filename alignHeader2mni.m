function alignHeader2mni(T1,T2,segOut,mri2mni)
% alignHeader2mni(T1,T2,segOut,mri2mni)
%
% Set the header info such that the voxel-to-world mapping maps directly
% into the MNI space.
%
% (c) Andrew Birnbaum, Parra Lab at CCNY
%     Yu (Andy) Huang
% April 2024

[dirname,segOutName] = fileparts(segOut);
if isempty(dirname), dirname = pwd; end

if exist(T1,'file')
%     disp('Aligning T1 MRI''s header to the MNI space ...');
    update_affine(T1,mri2mni,segOutName);
end

if ~isempty(T2)
    if exist(T2,'file')
%         disp('Aligning T2 MRI''s header to the MNI space ...');
        [~,t2name] = fileparts(T2);
        update_affine(T2,mri2mni,t2name);
    end
end

if exist([dirname filesep segOutName '_masks.nii'],'file')
%     disp('Aligning segmentation''s header to the MNI space ...');
    update_affine([dirname filesep segOutName '_masks.nii'],mri2mni,[segOutName '_masks']);
end

function update_affine(mri,v2w,outName)
nii = load_untouch_nii(mri);
nii.hdr.hist.srow_x = v2w(1,:);
nii.hdr.hist.srow_y = v2w(2,:);
nii.hdr.hist.srow_z = v2w(3,:);

nii.hdr.hist.qoffset_x = nii.hdr.hist.srow_x(4);
nii.hdr.hist.qoffset_y = nii.hdr.hist.srow_y(4);
nii.hdr.hist.qoffset_z = nii.hdr.hist.srow_z(4);
nii.hdr.hist.qform_code = 0;
nii.hdr.hist.sform_code = 1; % force it to be sform

[dirname,~,ext] = fileparts(mri);
if isempty(dirname), dirname = pwd; end
save_untouch_nii(nii,[dirname filesep outName '_MNI' ext]);