function MultiPrior(subj)

[dirname,baseFilename] = fileparts(subj);
if ~exist(fullfile(dirname, [baseFilename '_T1orT2_masks_MultiPrior_Segmentation.nii']), 'file')

    WARP_indiTPM(subj)
    
    cd ..
    cd ..
    
    SEGMENT(subj) 
else 
    disp('MultiPrior Segmentation File: ')
    disp([dirname, [baseFilename '_T1orT2_masks_MultiPrior_Segmentation.nii']])
    disp('Already Exists. Skipping Segmentation...')
end
end
