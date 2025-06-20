function viewElectrodes(subj, elec, gel, landmarks)
%
% 3D Visualization of skin, brain, electrodes, gel, and anatomical landmarks
% from segmented MRI data. Displays a rendered volume with overlaid surfaces
% and points for interactive inspection of electrode placement.
%
% INPUTS:
%   subj       - Path to the original NIfTI segmentation (.nii) file.
%                The function assumes the segmented mask is in a file named 
%                [subjName '_masks.nii'] in the same directory.
%   elec       - 3D matrix of the same size as segmentation, where electrode
%                regions have non-zero values.
%   gel        - 3D matrix of the same size as segmentation, where gel 
%                regions have non-zero values.
%   landmarks  - (Optional) Nx3 array of anatomical landmark coordinates. If
%                provided, the first four are labeled and displayed.
%
% USAGE NOTES:
% - This function assumes the tissue labels in the segmentation:
%     Label 5 = Skin
%     Label 2 = Brain
% - Landmarks will be displayed in red with text labels for Nasion, Inion,
%   Left Ear, and Right Ear.
%
% OUTPUT:
%   None. Launches an interactive 3D viewer with overlaid anatomy and electrodes.
%
% See also: getLandmarksManual, checkLandmarks
%
% (c) Andrew Birnbaum, Parra Lab at CCNY  
% 
% June 2025

    if nargin < 4
        landmarks = [];
    end

    % Load NIfTI segmentation and metadata
    [dirname, subjName, ext] = fileparts(subj);
    segOut = fullfile(dirname, [subjName, '_masks', ext]);
    nii_info = niftiinfo(segOut);
    nii_mask = niftiread(segOut);
    nii_elec = elec;
    nii_gel = gel;
    % nii_mask = flip(nii_mask, 2);      % Flip x-axis (left-right)
    % nii_elec = flip(nii_elec, 2);
    % nii_gel  = flip(nii_gel, 2);

    voxel_size = nii_info.PixelDimensions;

    % Clean segmentation: only keep relevant labels
    mask_skin  = double(nii_mask == 5); % Skin
    mask_brain = double(nii_mask == 2); % Brain
    electrodes = double(nii_elec > 0);
    gel = double(nii_gel > 0);

    % Smoothing
    mask_skin = imgaussfilt3(mask_skin, 1);
    mask_brain = imgaussfilt3(mask_brain, 1);
    electrodes = imgaussfilt3(electrodes, 1);
    gel = imgaussfilt3(gel, 1);

    % Create figure
    figure('Name', '3D Segmentation Viewer', ...
           'NumberTitle', 'off', ...
           'Position', [100, 100, 1200, 800], ...
           'Color', 'white');
    hold on;
    daspect(1 ./ voxel_size);

    % Plot skin (semi-transparent)
    if any(mask_skin(:))
        p1 = patch(isosurface(mask_skin, 0.5));
        p1.FaceColor = [229/255, 181/255, 161/255]; % light skin
        p1.EdgeColor = 'none';
        p1.FaceAlpha = 0.2;
    end

    % Plot brain (pink)
    if any(mask_brain(:))
        p2 = patch(isosurface(mask_brain, 0.5));
        p2.FaceColor = [1, 0.6, 0.8]; % pink
        p2.EdgeColor = 'none';
        p2.FaceAlpha = 1;
    end

    % Plot electrodes (blue)
    if any(electrodes(:))
        p3 = patch(isosurface(electrodes, 0.5));
        p3.FaceColor = 'blue';
        p3.EdgeColor = 'none';
        p3.FaceAlpha = .8;
    end

    % Plot gel (green)
    if any(gel(:))
        p4 = patch(isosurface(gel, 0.5));
        p4.FaceColor = 'green';
        p4.EdgeColor = 'none';
        p4.FaceAlpha = .8;
    end

    % Plot landmarks if provided
    if ~isempty(landmarks)
        keepIdx = [1, 2, 3, 4];
        scatter3(landmarks(keepIdx, 2), landmarks(keepIdx, 1), landmarks(keepIdx, 3), ...
            200, 'red', 'filled');
        labels = {'     Nasion', '     Inion', '     Left Ear', '     Right Ear'};
        for i = 1:length(keepIdx)
            text(landmarks(keepIdx(i), 2), landmarks(keepIdx(i), 1), landmarks(keepIdx(i), 3), ...
                labels{i}, 'FontSize', 14, 'Color', 'red', 'FontWeight', 'bold');
        end
    end

    view(3);
    axis off;
    grid off;
    light('Position', [-1, 0, 0], 'Style', 'infinite');
    light('Position', [1, 0, 1], 'Style', 'infinite');
    lighting phong;
    rotate3d on;
    title('3D View: Skin, Brain, Electrodes, and Landmarks');
    % Save figure (optional — change path as needed)
    saveas(gcf, fullfile(dirname, [subjName '_3DView.fig']))
end

