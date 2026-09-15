function landmarks = checkLandmarks(mask,landmarks)
% landmarks = checkLandmarks(mask,landmarks)
%
% GUI for visual inspection and optional modification of anatomical landmarks
% (Nasion, Inion, Right Ear, Left Ear) on a segmented MRI head.
%
% This function loads the segmentation, renders a 3D visualization of the skin and
% brain, and overlays the current landmark locations. The user can confirm the
% existing selections or (re)launch the manual selection tool to modify them.
%
% USAGE NOTES:
% - Segmentation labels: 5 = skin, 2 = brain.
% - Automatically plots landmarks and allows interaction in 3D space.
% - Includes buttons to confirm or modify landmark selection.
%
% See also: getLandmarksManual
%
% (c) Andrew Birnbaum, Parra Lab at CCNY
%     Yu (Andy) Huang
% June 2025

% % Initialize updated flag only if it is not already set
% if ~exist('updated', 'var') || isempty(updated)
%     updated = false;
% end

% smooth=0;

disp('=============    OPENING MANUAL GUI ...    =============')

% Load NIfTI segmentation
skin = single(mask.img == 5);
brain = single(mask.img == 2);
% Apply Gaussian smoothing to the segmentation masks
% if smooth
skin = imgaussfilt3(skin, 1);
% end

embedded = false;
action = 'confirm';
if isappdata(0, 'ROAST_GUI_PLOT_TARGET')
    [panel, fig] = roastGuiAddPlot('Check landmarks', 'interactive');
    if ~isempty(panel)
        embedded = true;
        set(panel, 'BackgroundColor', 'white');
    end
end

if ~embedded
    fig = figure('Name', '3D Viewer. Please rotate and inspect the landmarks.', ...
                 'NumberTitle', 'off', ...
                 'Position', [100, 100, 1200, 800]);
    panel = uipanel('Parent', fig, 'Position', [0, 0, 1, 1]);
end

% Create new axes inside the panel
if embedded
    imgAx = axes('Parent', panel, 'Position', [0.02, 0.10, 0.96, 0.86]);
else
    imgAx = axes('Parent', panel, 'Position', [0, 0, 1, 1]);
end
axis(imgAx, 'off');
axes(imgAx);

% Render skin layer (start fresh)
p = patch(isosurface(skin, 0.5,'noshare'));
p.FaceColor = '#E5B5A1';
p.EdgeColor = 'none';
alpha(p, 0.3);  % Make the brain semi-transparent
hold on;

p_brain = patch(isosurface(brain, 0.5,'noshare'));
p_brain.FaceColor = [1, 0.6, 0.8];  % Light blue color
p_brain.EdgeColor = 'none';

% Indices of points to keep (All four landmarks)
keepIdx = [1, 2, 3, 4];
% Scatter plot for landmarks

% if smooth
% scatter3(smoothLandmarks(keepIdx, 2),smoothLandmarks(keepIdx, 1), smoothLandmarks(keepIdx, 3), 200, 'red', 'filled');
% else
scatter3(landmarks(keepIdx, 2), landmarks(keepIdx, 1), landmarks(keepIdx, 3), 200, 'red', 'filled');
% end
% Labels for points
labels = {'     Nasion', '     Inion', '     Right Ear', '     Left Ear'};

for i = 1:length(keepIdx)
%     if smooth
%     text(smoothLandmarks(keepIdx(i), 2),smoothLandmarks(keepIdx(i), 1), smoothLandmarks(keepIdx(i), 3), labels{i}, ...
%          'FontSize', 14, 'Color', 'red', 'FontWeight', 'bold');
%     else
    text(landmarks(keepIdx(i), 2), landmarks(keepIdx(i), 1), landmarks(keepIdx(i), 3), labels{i}, ...
     'FontSize', 14, 'Color', 'red', 'FontWeight', 'bold');
%     end
end

% Improve 3D visualization
view(3);
axis equal;
axis ij; % use axis ij, so that we can LR flip the axis from patch command, without flipping the data or hacking the order of labels

% Set two lights, one on the left and one on the right
light('Position', [-1, 0, 0], 'Style', 'infinite'); % Back light
light('Position', [1, 0, 1], 'Style', 'infinite');  % Front light
lighting phong;

% Enable 3D rotation interaction
rotate3d(fig, 'on');
try
    imgAx.Interactions = [rotateInteraction zoomInteraction dataTipInteraction];
catch
end

% Update title
title('Selecting Landmarks ...');

% Define selectedPoints correctly
selectedPoints = landmarks(keepIdx, :);
disp('Current voxel coordinates of landmarks (nasion, inion, right ear, left ear):');
disp(selectedPoints);
disp('3D interaction enabled. Rotate and inspect the landmarks.');

if embedded
    uicontrol(panel, 'Style', 'pushbutton', 'String', 'Confirm', ...
              'Units', 'normalized', ...
              'Position', [0.78 0.025 0.09 0.055], ...
              'FontSize', 12, ...
              'Callback', @(~,~) confirmSelection());
    uicontrol(panel, 'Style', 'pushbutton', 'String', 'Modify', ...
              'Units', 'normalized', ...
              'Position', [0.88 0.025 0.09 0.055], ...
              'FontSize', 12, ...
              'Callback', @(~,~) modifyLandmarks());
else
    uicontrol(panel, 'Style', 'pushbutton', 'String', 'Confirm', ...
              'Units', 'pixels', ...
              'Position', [1000 10 80 40], ...
              'FontSize', 14, ...
              'Callback', @(~,~) confirmSelection());
    uicontrol(panel, 'Style', 'pushbutton', 'String', 'Modify', ...
              'Units', 'pixels', ...
              'Position', [1100 10 80 40], ...
              'FontSize', 14, ...
              'Callback', @(~,~) modifyLandmarks());
end

if embedded
    roastGuiShowPlot(panel);
else
    drawnow;
end

% Wait for the figure to close before returning landmarks
uiwait(fig);  % Wait here until the user confirms or modifies
if isRoastGuiCancelRequested() || (embedded && ~ishandle(fig))
    error('ROAST run cancelled from ROAST Launcher.');
end
if embedded && ishandle(panel)
    roastGuiRemovePlot(panel);
end
if strcmp(action, 'modify')
    updatedLandmarks = getLandmarksManual(mask);
    landmarks(keepIdx, :) = reorderLandmarks(updatedLandmarks);
    landmarks = checkLandmarks(mask,landmarks);
end

% Callback function for confirming landmarks
function confirmSelection()
    uiresume(fig);
    if ~embedded && ishandle(fig)
        close(fig);
    end
end

% Nested function to call segmentation_viewer_nii
function modifyLandmarks()
    action = 'modify';
    uiresume(fig);
    if ~embedded && ishandle(fig)
        close(fig);
    end
end

function reorderedLandmarks = reorderLandmarks(inputMatrix)
    % Define the new row order
    newOrder = [1, 5, 2, 3];
    % Reorder the rows of the matrix according to the new order
    reorderedLandmarks = inputMatrix(newOrder, :);
end

function tf = isRoastGuiCancelRequested()
tf = isappdata(0, 'ROAST_GUI_CANCEL_REQUESTED');
end

% Only return landmarks after the figure is closed (confirmed or modified)
% The function will return after uiwait completes
end
