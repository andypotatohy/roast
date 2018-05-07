function srcR = realignT2(src,ref)

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
matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.interp = 1;
matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.wrap = [0 0 0];
matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.mask = 0;
matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.suffix = '_aligned';

spm_jobman('run',matlabbatch);

[pth,nam,ext] = fileparts(src);
srcR = [pth filesep nam '_aligned' ext];
disp(['Done. Aligned T2 is saved as ' srcR ', and will be used by ROAST.']);