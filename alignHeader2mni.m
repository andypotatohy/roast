function alignHeader2mni(subjRasRSPD)
    [dirname,baseFilename,ext] = fileparts(subjRasRSPD);
    seg8 = load([dirname,'/',baseFilename, '_seg8.mat']);
    seg8Affine = seg8.Affine;
    if exist(fullfile([dirname,'/',baseFilename,'_MNI',ext]), 'file')
    %disp('File exists. Skipping Conversion...')
    else
    disp('Converting MRI to MNI...')
  
    parts = strsplit(baseFilename, '_');
    newParts = parts(1:end-2);
    resultString = strjoin(newParts, '_');
    mri = [dirname filesep resultString ext];
    [Vol,info] = update_affine(mri,seg8Affine);
    niftiwrite(Vol,[dirname,'/',baseFilename,'_MNI',ext],info)
    disp([dirname,'/',baseFilename,'_MNI',ext])
    end
   
    if (exist(fullfile([dirname,'/',baseFilename,'_masks.nii']), 'file') &&...
       ~exist(fullfile([dirname,'/',baseFilename,'_masks_MNI.nii']), 'file'))
        segmentation = [dirname,'/',baseFilename, '_masks.nii'];
        [Vol,info] = update_affine(segmentation,seg8Affine);
        niftiwrite(Vol,[dirname,'/',baseFilename,'_masks_MNI.nii'],info)
        disp([dirname,'\',baseFilename,'_masks_MNI.nii'])   
        disp('Converting Segmentation to MNI...')
    else
        %disp('File exists. Skipping Conversion...')
    end 
end

function [Vol,info] = update_affine(subj,seg8Affine)
            info=niftiinfo(subj);
            Vol = niftiread(subj);
            Affine = info.Transform.T;
            newAffine = seg8Affine*Affine';
            info.Transform.T = newAffine';
            %Makes sure both MRI and Segmentation are the same form
            info.TransformName = 'Sform';
end