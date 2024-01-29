function Multipriors(subj)

[dirname,baseFilename] = fileparts(subj);
if ~exist(fullfile(dirname, [baseFilename '_masks.nii']), 'file')

    WARP_indiTPM(subj)
   
    SEGMENT(subj) 
else 
    disp('MultiPriors Segmentation File: ')
    disp([dirname, [baseFilename '_masks.nii']])
    disp('Already Exists. Skipping Segmentation...')
end
end
