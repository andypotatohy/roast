function renameFiles(dirname, baseFilenameRasRSPD, multiprior)   
    if multiprior
        disp('Renaming MultiPrior Segmentation ...')
        % Check if the MultiPrior segmentation file exists
        if exist(fullfile(dirname, [baseFilenameRasRSPD '_T1orT2_masks_MultiPrior_Segmentation.nii']), 'file')
            disp(['File exists. Renamamed as ',[baseFilenameRasRSPD '_T1orT2_masks.nii']])
            % Copy and rename the MultiPrior segmentation file to T1_orT2_masks.nii
            copyfile(fullfile(dirname, [baseFilenameRasRSPD '_T1orT2_masks_MultiPrior_Segmentation.nii']), ...
                fullfile(dirname, [baseFilenameRasRSPD '_T1orT2_masks.nii']));
        end
    else
        disp('Renaming Roast Segmentation ...')
        % Check if T1orT2_masks_Roast_Segmentation file exists
        if exist(fullfile(dirname, [baseFilenameRasRSPD '_T1orT2_masks_Roast_Segmentation.nii']), 'file')
            disp(['File exists. Renamed as ',[baseFilenameRasRSPD '_T1orT2_masks.nii']])
            % Copy and rename the T1orT2_masks file to T1orT2_masks.nii
            copyfile(fullfile(dirname, [baseFilenameRasRSPD '_T1orT2_masks_Roast_Segmentation.nii']), ...
                fullfile(dirname, [baseFilenameRasRSPD '_T1orT2_masks.nii']));
        end
        % Check if T1andT2_masks_Roast_Segmentation file exists
        if exist(fullfile(dirname, [baseFilenameRasRSPD '_T1andT2_masks_Roast_Segmentation.nii']), 'file')
            disp(['File exists. Renamed as ',[baseFilenameRasRSPD '_T1andT2_masks.nii']])
            % Rename the T1andT2_masks file to T1andT2_masks.nii
            movefile(fullfile(dirname, [baseFilenameRasRSPD '_T1andT2_masks_Roast_Segmentation.nii']), ...
                fullfile(dirname, [baseFilenameRasRSPD '_T1andT2_masks.nii']));
        end
    end
end
