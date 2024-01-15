function multipriors(subj)

[dirname,baseFilename] = fileparts(subj);
if ~exist(fullfile(dirname, [baseFilename '_T1orT2_masks_MultiPriors_Segmentation.nii']), 'file')

    WARP_indiTPM(subj)
   
    SEGMENT(subj) 
else 
    disp('MultiPriors Segmentation File: ')
    disp([dirname, [baseFilename '_T1orT2_masks_MultiPriors_Segmentation.nii']])
    disp('Already Exists. Skipping Segmentation...')
end
end
