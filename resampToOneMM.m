function resampToOneMM(mri)
% resampToOneMM(mri)
% 
% Resample the MRI into 1 mm isotropic resolution if it has non-1 mm
% resolution.
% 
% (c) Yu (Andy) Huang, Parra Lab at CCNY
% yhuang16@citymail.cuny.edu
% April 2018

V = spm_vol(mri);

pixdim = [V.mat(1,1) V.mat(2,2) V.mat(3,3)];

if all(pixdim==1)
    
    warning(['The MRI ' mri ' already has a 1 mm isotropic resolution. No need to resample.']);
    return;
    
else
    
    [dirname,baseFilename,ext] = fileparts(mri);
    
    if exist([dirname filesep 'r' baseFilename ext],'file')
        
        disp([mri ' has already been resampled to 1mm resolution and saved as ' dirname filesep 'r' baseFilename ext]);
        disp('Please use that file as the input for ROAST main function.');
        return;
        
    else
        
        disp(['Resampling ' mri ' to 1 mm isotropic resolution...']);
        
        voxsiz = [1 1 1]; % new voxel size in mm
        bb = spm_get_bbox(V);
        VV(1:2) = V;
        VV(1).mat = spm_matrix([bb(1,:) 0 0 0 voxsiz])*spm_matrix([-1 -1 -1]);
        VV(1).dim = ceil(VV(1).mat \ [bb(2,:) 1]' - 0.1)';
        VV(1).dim = VV(1).dim(1:3);
        spm_reslice(VV,struct('mean',false,'which',1,'interp',7)); % 1 for linear
        
        disp([mri ' has been resampled to 1 mm isotropic resolution, and is saved as:']);
        disp([dirname filesep 'r' baseFilename ext]);
        disp('Please use this file as the input for ROAST main function.');
        
    end
    
end