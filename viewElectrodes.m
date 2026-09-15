function viewElectrodes(mask,elec,gel,landmarks,imgHdr,uniTag)
% viewElectrodes(mask,elec,gel,landmarks,imgHdr,uniTag)
%
% 3D Visualization of skin, brain, electrodes, gel, and anatomical landmarks
% from segmented MRI data. Displays a rendered volume with overlaid surfaces
% and points for interactive inspection of electrode placement.
%
% USAGE NOTES:
% - This function assumes the tissue labels in the segmentation:
%     Label 5 = Skin
%     Label 2 = Brain
% - Landmarks will be displayed in red with text labels for Nasion, Inion,
%   Right Ear, and Left Ear.
%
% See also: getLandmarksManual, checkLandmarks
%
% (c) Andrew Birnbaum, Parra Lab at CCNY
%     Yu (Andy) Huang
% June 2025

% nii_mask = flip(nii_mask, 2);      % Flip x-axis (left-right)
% nii_elec = flip(nii_elec, 2);
% nii_gel  = flip(nii_gel, 2);

mask_skin = imgaussfilt3(single(mask.img == 5), 1);
mask_brain = imgaussfilt3(single(mask.img == 2), 1);

elec = imgaussfilt3(single(elec.img>0), 1);
gel = imgaussfilt3(single(gel.img>0), 1);

% Create figure or draw into the ROAST GUI when it is active.
embeddedTabs = [];
embeddedTarget = [];
if isappdata(0, 'ROAST_GUI_PLOT_TARGET')
    embeddedTarget = getappdata(0, 'ROAST_GUI_PLOT_TARGET');
end
if isappdata(0, 'ROAST_GUI_SLICE_TABS')
    embeddedTabs = getappdata(0, 'ROAST_GUI_SLICE_TABS');
    if ~ishandle(embeddedTabs)
        embeddedTabs = [];
    end
end

if ~isempty(embeddedTarget)
    [panel, fh] = roastGuiAddPlot('Electrodes', '3d');
    if isempty(panel)
        embeddedTarget = [];
    else
        ax = axes('Parent', panel, 'Units', 'normalized', 'Position', [0.02 0.02 0.96 0.94]);
    end
end

if isempty(embeddedTarget) && ~isempty(embeddedTabs)
    tab = uitab(embeddedTabs, 'Title', 'Electrodes');
    setappdata(tab, 'ROAST_GUI_PLOT_MODE', '3d');
    panel = uipanel(tab, 'Units', 'normalized', 'Position', [0 0 1 1], ...
        'BorderType', 'none', 'BackgroundColor', 'white');
    ax = axes('Parent', panel, 'Units', 'normalized', 'Position', [0.02 0.02 0.96 0.92]);
    set(embeddedTabs, 'SelectedTab', tab);
    fh = ancestor(panel, 'figure');
elseif isempty(embeddedTarget)
    fh=figure('Name', '3D Viewer. Please rotate and inspect.', ...
           'NumberTitle', 'off', ...
           'Position', [100, 100, 1200, 800], ...
           'Color', 'white');
    ax = axes('Parent', fh);
end
isEmbedded = ~isempty(embeddedTarget) || ~isempty(embeddedTabs);
hold(ax, 'on');
daspect(ax, 1 ./ [imgHdr(1).mat(1,1),imgHdr(1).mat(2,2),imgHdr(1).mat(3,3)]);

% Plot skin (semi-transparent)
if any(mask_skin(:))
    p1 = patch(ax, isosurface(mask_skin, 0.5,'noshare'));
    p1.FaceColor = [229/255, 181/255, 161/255]; % light skin
    p1.EdgeColor = 'none';
    p1.FaceAlpha = 0.2;
end

% Plot brain (pink)
if any(mask_brain(:))
    p2 = patch(ax, isosurface(mask_brain, 0.5,'noshare'));
    p2.FaceColor = [1, 0.6, 0.8]; % pink
    p2.EdgeColor = 'none';
    p2.FaceAlpha = 1;
end

% Plot electrodes (blue)
if any(elec(:))
    p3 = patch(ax, isosurface(elec, 0.5,'noshare'));
    p3.FaceColor = 'blue';
    p3.EdgeColor = 'none';
    p3.FaceAlpha = .8;
end

% Plot gel (green)
if any(gel(:))
    p4 = patch(ax, isosurface(gel, 0.5,'noshare'));
    p4.FaceColor = 'green';
    p4.EdgeColor = 'none';
    p4.FaceAlpha = .8;
end

% Plot landmarks if provided
if ~isempty(landmarks)
    keepIdx = [1, 2, 3, 4];
    scatter3(ax, landmarks(keepIdx, 2), landmarks(keepIdx, 1), landmarks(keepIdx, 3), ...
        200, 'red', 'filled');
    labels = {'     Nasion', '     Inion', '     Right Ear', '     Left Ear'};
    for i = 1:length(keepIdx)
        text(ax, landmarks(keepIdx(i), 2), landmarks(keepIdx(i), 1), landmarks(keepIdx(i), 3), ...
            labels{i}, 'FontSize', 14, 'Color', 'red', 'FontWeight', 'bold');
    end
end

view(ax, 3);
axis(ax, 'ij'); % use axis ij, so that we can LR flip the axis from patch command, without flipping the data or hacking the order of labels
axis(ax, 'off');
grid(ax, 'off');
light(ax, 'Position', [-1, 0, 0], 'Style', 'infinite');
light(ax, 'Position', [1, 0, 1], 'Style', 'infinite');
lighting(ax, 'phong');
if ~isEmbedded
    title(ax, ['Electrode placement in Simulation: ' uniTag]);
end
if isEmbedded
    rotate3d(fh, 'on');
    try
        ax.Interactions = [rotateInteraction zoomInteraction dataTipInteraction];
    catch
        try
            enableDefaultInteractivity(ax);
        catch
        end
    end
else
    rotate3d(fh, 'on');
end
if ~isEmbedded
    movegui(fh,'center')
end
drawnow
%     % Save figure (optional — change path as needed)
%     saveas(gcf, fullfile(dirname, [subjName '_3DView.fig']));
