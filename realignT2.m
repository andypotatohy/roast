function srcR = realignT2(src,ref)

refData = load_untouch_nii(ref);
srcData = load_untouch_nii(src);
if any(size(refData.img)~=size(srcData.img))
    
    [pth,nam,ext] = fileparts(src);
    srcR = [pth filesep nam '_aligned' ext];

    if exist(srcR,'file')
        srcData = load_untouch_nii(srcR);
        if all(size(refData.img)==size(srcData.img))
            warning(['Original T2 ' src ' is not aligned to T1 image space, but an aligned T2 has been found (' srcR '). ROAST will use that file as the input.']);
            return;
        end
    end
    
    warning(['T2 image ' src ' is not aligned to T1 image space. ROAST will align it now.']);
    
    if exist('matlabbatch','var')
        clear matlabbatch;
    end
    
    disp('Aligning T2 to T1...');
    
    spm_jobman('initcfg');
    
    matlabbatch{1}.spm.spatial.coreg.estwrite.ref = {[ref ',1']};
    matlabbatch{1}.spm.spatial.coreg.estwrite.source = {[src ',1']};
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