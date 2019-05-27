function srcR = realignT2(src,ref)

refData = load_untouch_nii(ref);
srcData = load_untouch_nii(src);
if any(refData.hdr.dime.dim(1:5)~=srcData.hdr.dime.dim(1:5)) || ...
        any(refData.hdr.dime.pixdim(1:4)~=srcData.hdr.dime.pixdim(1:4)) || ...
        any(refData.hdr.hist.srow_x~=srcData.hdr.hist.srow_x) || ...
        any(refData.hdr.hist.srow_y~=srcData.hdr.hist.srow_y) || ...
        any(refData.hdr.hist.srow_z~=srcData.hdr.hist.srow_z)
    
    [~,nam] = fileparts(ref);
    [pth,~,ext] = fileparts(src);
        
    srcR = fullfile(pth ,[nam '_T2_aligned' ext]);

    if exist(srcR,'file')
        srcData = load_untouch_nii(srcR);
        if all(refData.hdr.dime.dim(1:5)==srcData.hdr.dime.dim(1:5)) && ...
                all(refData.hdr.dime.pixdim(1:4)==srcData.hdr.dime.pixdim(1:4)) && ...
                all(refData.hdr.hist.srow_x==srcData.hdr.hist.srow_x) && ...
                all(refData.hdr.hist.srow_y==srcData.hdr.hist.srow_y) && ...
                all(refData.hdr.hist.srow_z==srcData.hdr.hist.srow_z)
            warning(['Original T2 ' src ' is not aligned to T1 image space, but an aligned T2 has been found (' srcR '). ROAST will use that file as the input.']);
            return;
        end
    end
    
    warning(['T2 image ' src ' is not aligned to T1 image space. ROAST will align it now.']);
    disp('Aligning T2 to T1...');
    
    % make a copy as SPM will overwrite the header in source image
    srcForRealign = fullfile(pth ,[nam '_T2' ext]);
    copyfile(src,srcForRealign);
    
    if exist('matlabbatch','var')
        clear matlabbatch;
    end
    
    spm_jobman('initcfg');
    
    matlabbatch{1}.spm.spatial.coreg.estwrite.ref = {[ref ',1']};
    matlabbatch{1}.spm.spatial.coreg.estwrite.source = {[srcForRealign ',1']};
    matlabbatch{1}.spm.spatial.coreg.estwrite.other = {''};
    matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.cost_fun = 'nmi';
    matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.sep = [4 2];
    matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
    matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.fwhm = [7 7];
    matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.interp = 7; % 1 for linear, 7 is highest degree (most accurate, slowest). Use 7 here to be consistent with resampToOneMM()
    matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.wrap = [0 0 0];
    matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.mask = 0;
    matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.prefix = '_aligned'; % keep using 'prefix' as it's bad to hack SPM variable names
    
    spm_jobman('run',matlabbatch);
    
    disp(['Aligned T2 is saved as ' srcR ', and will be used by ROAST.']);
    
else
    
    srcR = src;
    
end