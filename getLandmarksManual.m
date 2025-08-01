function [landmarks, smoothLandmarks] = getLandmarksManual(mask)
% Interactive GUI for manual selection of five anatomical landmarks 
% (Nasion, Right Ear, Left Ear, and two Inion points) from a segmented MRI.
% Visualizes the 3D skin and skull surfaces and guides the user step-by-step.
%
% INPUT:
%   filepath - path to NIfTI (.nii) segmentation file with labeled tissues:
%              label 5 = skin, label 4 = skull.
%
% OUTPUT:
%   landmarks        - Nx3 array of selected voxel coordinates (original).
%   smoothLandmarks  - Nx3 array of corresponding smoothed coordinates.
%
% Developed for clinical tDCS modeling applications, where robust and accurate
% anatomical landmark selection is essential for electrode placement.
%
% (c) Andrew Birnbaum, Parra Lab at CCNY  
%
% June 2025

disp('=============    LOADING GUI FOR MODIFY ...    =============')

% Load NIfTI segmentation
skin = single(mask.img == 5);  % Skin layer
skull = single(mask.img == 4); % Skull layer
[rows, cols, slices] = size(mask.img);

% Apply Gaussian smoothing to the segmentation masks
skinSmooth = imgaussfilt3(skin, 1);
skullSmooth = imgaussfilt3(skull, 1);

% Get screen size
screenSize = get(0, 'ScreenSize');
figWidth = 1200;  % Increased width
figHeight = 800;  % Increased height
figX = (screenSize(3) - figWidth) / 2;
figY = (screenSize(4) - figHeight) / 2;

% Create figure and axes
fig = figure('Name', 'Selecting Landmarks ...', ...
             'NumberTitle', 'off', ...
             'Position', [figX, figY, figWidth, figHeight], 'MenuBar', 'none');

ax = axes('Parent', fig, 'NextPlot', 'add', 'DataAspectRatio', [1 1 1], 'Position', [0.05, 0.1, 0.7, 0.8]);
hold(ax, 'on');
axis(ax, 'off');
grid off;

% Render segmentation (Skin first)
skinPatch = patch(isosurface(skinSmooth, 0.5,'noshare'));
skinPatch.FaceColor = '#E5B5A1'; % Skin is brown
skinPatch.EdgeColor = 'none';

% Render skull (after skin)
skullPatch = patch(isosurface(skullSmooth, 0.5,'noshare')); 
skullPatch.FaceColor = '#F1D691'; % Skull is light yellow
skullPatch.EdgeColor = 'none';

axis(ax,'ij'); % now the data visualized genuinely follow the real orientation, so that the code below will be clean

view(90, 0); % Initial view
title('Select the Nasion: The Point in Between the Eyes','FontSize', 16);

% Apply lighting
updateLighting([1, 0, 0]); % Initial front light

% Adjust axes limits
xlim([1 cols]); ylim([1 rows]); zlim([1 slices]);

% Store selected points
selectedPoint = [];
selectedPointSmooth = [];
selectedPointSkull = [];
selectionPhase = 0;  % Start at 0 for Nasion
lastClickedPoint = [];
lastClickedPointSmooth  = [];
lastClickedPointSkull = [];
scatterHandle = []; % Handle for the scatter plot
% Button to submit selection
submitButton = uicontrol('Style', 'pushbutton', 'String', 'Submit', ...
                         'Position', [570 5 100 40], ...
                         'FontSize', 16, ...
                         'Callback', @submitSelection, ...
                         'Enable', 'off');

panel = uipanel('Parent', fig, 'Position', [0.7, 0.4, 0.3, 0.5], 'BorderType', 'none'); 

% Create a centered title ABOVE the panel
titleText = uicontrol('Parent', fig, 'Style', 'text', ...
    'String', 'Example Selection', 'Units', 'normalized', ...
    'Position', [0.7, 0.87, 0.3, 0.05], ... % Positioned above the panel
    'FontSize', 16, 'FontWeight', 'bold', ...
    'HorizontalAlignment', 'center');

imgAx = axes('Parent', panel, 'Position', [0, 0, 1, 1]); % Keeps image full inside panel
axis(imgAx, 'off');
updateInstructions('lib/screenshots/Nasion_New.png');

% Callback for mouse click
set(fig, 'WindowButtonDownFcn', @onMouseClick);

function onMouseClick(~, ~)
    clickPoint = get(ax, 'CurrentPoint');
    clickPoint = round(clickPoint(1, :));

    if ismember(selectionPhase, [0, 4])  % Nasion and Inion
        x = clickPoint(2);
        z = clickPoint(3);
        if x >= 1 && x <= rows && z >= 1 && z <= slices
            % fprintf('Selected X = %d, Z = %d\n', x, z);
            yPoints = 1:cols;

            % Get both smooth and original values
            skinValues = skin(x, yPoints, z);
            skullValues = skull(x, yPoints, z);
            skinSmoothValues = skinSmooth(x, yPoints, z);

            % Identify zero & non-zero points
            zeroPoints = yPoints((skinValues == 0)&(skullValues == 0));
            zeroPointsSmooth = yPoints((skinSmoothValues == 0)&(skullValues == 0));
            nonZeroPoints = yPoints(skinValues ~= 0);
            nonZeroSkullPoints = yPoints(skullValues ~= 0);
            nonZeroPointsSmooth = yPoints(skinSmoothValues ~= 0);

            if isempty(nonZeroPoints)
                disp('Click is out of bounds. Please click following the example.');
                return;
            end

            if ~isempty(nonZeroPointsSmooth)
                if selectionPhase == 0  % Nasion (max y)
                    ySmooth = max(nonZeroPointsSmooth);
                else  % Inion (min y)
                    ySmooth = min(nonZeroPointsSmooth);
                    ySkull = min(nonZeroSkullPoints);
                    % Find the first background (zero) point after the skull
                    if ~isempty(nonZeroSkullPoints) && ~isempty(zeroPointsSmooth)
                     % Search for the first zero after the skull point
                     yBackground = find(zeroPointsSmooth < ySkull, 1, 'last'); 
                     ySmooth = zeroPointsSmooth(yBackground)+1;
                    end
                end
            end


            if ~isempty(nonZeroPoints)
                if selectionPhase == 0  % Nasion (max y)
                    y = max(nonZeroPoints);
                else  % Inion (min y)
                    y = min(nonZeroPoints);
                    ySkull = min(nonZeroSkullPoints);
                    % Find the first background (zero) point after the skull
                    if ~isempty(nonZeroSkullPoints) && ~isempty(zeroPoints)
                     % Search for the first zero after the skull point
                     yBackground = find(zeroPoints < ySkull, 1, 'last'); 
                     y  = zeroPoints(yBackground)+1;
                    end
        
                end
            end
            

            % Remove previous scatter plot if it exists
            if ~isempty(scatterHandle) && isvalid(scatterHandle)
                delete(scatterHandle);
            end

            % Plot using smoothed values
            scatterHandle = scatter3(ySmooth, x, z, 100, 'blue', 'filled');

            % Store both smooth and non-smooth selections
            lastClickedPoint = [x, y, z];  
            lastClickedPointSmooth = [x, ySmooth, z];

            disp("Clicked at:"); disp(lastClickedPoint);
%             disp("Original:"); disp(lastClickedPoint);
%             disp("Smoothed:"); disp(lastClickedPointSmooth);

            % Enable submission
            set(submitButton, 'Enable', 'on');
        else
            disp('Click is out of bounds. Please click following the example.');
        end

    elseif ismember(selectionPhase, [1, 2, 3])  % Right Ear and Left Ear
        y = clickPoint(1);
        z = clickPoint(3);
        if y >= 1 && y <= cols && z >= 1 && z <= slices
            % fprintf('Selected Y = %d, Z = %d\n', y, z);
            xPoints = 1:rows;

            % Get both smooth and original values
            skinValues = skin(xPoints, y, z);
            skinSmoothValues = skinSmooth(xPoints, y, z);
            skullValues = skull(xPoints, y, z);

            % Identify non-zero points for both
            nonZeroPoints = xPoints(skinValues ~= 0);
            nonZeroPointsSmooth = xPoints(skinSmoothValues ~= 0);
            nonZeroSkullPoints = xPoints(skullValues ~= 0);

            if isempty(nonZeroPoints)
                disp('Click is out of bounds. Please click following the example.');
                return;
            end

            if ~isempty(nonZeroPointsSmooth)
                if selectionPhase == 1  % Right Ear (max x)
                    xSmooth = max(nonZeroPointsSmooth);
                else  % Left Ear (min x)
                    xSmooth = min(nonZeroPointsSmooth);
                end
            end

            if ~isempty(nonZeroPoints)
                if selectionPhase == 1  % Right Ear (max x)
                    x = max(nonZeroPoints);
                else  % Left Ear (min x)
                    x = min(nonZeroPoints);
                end
            end
            
            % Remove previous scatter plot if it exists
            if ~isempty(scatterHandle) && isvalid(scatterHandle)
                delete(scatterHandle);
            end

            % Remove previous grid lines if they exist
            previousGridLines = findall(ax, 'Type', 'Surface');
            if ~isempty(previousGridLines)
                delete(previousGridLines);
            end

            % Plot using smoothed values
            scatterHandle = scatter3(y, xSmooth, z, 100, 'blue', 'filled');

            % Store both smooth and non-smooth selections
            lastClickedPoint = [x, y, z];  
            lastClickedPointSmooth = [xSmooth, y, z];

            disp("Clicked at:"); disp(lastClickedPoint);
%             disp("Original:"); disp(lastClickedPoint);
%             disp("Smoothed:"); disp(lastClickedPointSmooth);

            if selectionPhase == 3
                drawGridLines(z);
            end

            set(submitButton, 'Enable', 'on');
        else
            disp('Click is out of bounds. Please click following the example.');
        end
    end
end

function drawGridLines(z)
    % Get axis limits
    xLimits = xlim(ax);
    yLimits = ylim(ax);

    % Draw grid lines only along the axes, no diagonals
    hold on;
    % Horizontal line at the selected point
    plot3(xLimits(1):xLimits(2),ones(xLimits(2),1),ones(xLimits(2),1)*z,'--k','linewidth',2); 
    plot3(ones(yLimits(2),1),yLimits(1):yLimits(2),ones(yLimits(2),1)*z,'--k','linewidth',2);
    % this is much faster than the surf command
end

function submitSelection(~, ~)
    if ~isempty(lastClickedPoint)
        selectedPoint = [selectedPoint; lastClickedPoint];  % Store point
        selectedPointSmooth = [selectedPointSmooth; lastClickedPointSmooth];
        lastClickedPoint = [];  % Reset last clicked point
        lastClickedPointSmooth = [];
             % Remove the scatter plot after submission
        if exist('scatterHandle', 'var') && isvalid(scatterHandle)
            delete(scatterHandle);
            scatterHandle = [];  % Reset handle
        end
        set(submitButton, 'Enable', 'off');  % Disable button until next selection

        % Move to next phase
        switch selectionPhase
            case 0
                selectionPhase = 1;
                view(0, 0);
                updateLighting([1, -1, 1]);
                title('Select the Right Ear: In the Front Middle');
                updateInstructions('lib/screenshots/Right_Ear_New.png');
            case 1
                selectionPhase = 2;
                view(180, 0);
                updateLighting([-1, 0, 1]);
                title('Select the Left Ear: In the Front Middle');
                updateInstructions('lib/screenshots/Left_Ear_New.png');
            case 2
                selectionPhase = 3;
                if exist('skinPatch', 'var')
                    delete(skinPatch);  % Delete the skin patch
                end
                view(180, 0);
                updateLighting([-1, 0, 1]);
                title('Select the Inion: The Skull Begins to Slope Inwards');
                updateInstructions('lib/screenshots/Inion1_New.png');
            case 3
                selectionPhase = 4;
                view(-90, 0);
                updateLighting([-1, -1, 1]);
                title('Select the Inion: In the Middle of the Skull');
                updateInstructions('lib/screenshots/Inion2_New.png');
            case 4
                selectionPhase = 5;
                set(submitButton, 'Enable', 'off');
                close(gcf);

        end
    end
end

function updateLighting(position)
    delete(findall(gcf, 'Type', 'light')); % Clear previous lights
    light('Position', position, 'Style', 'local');
    camlight;
    lighting gouraud; % Smoother lighting
end

function updateInstructions(imagePath)
    img = imread(imagePath); % Load the new image
    imshow(img, 'Parent', imgAx); % Display the new image
end

uiwait(fig); % Wait for the figure to close before returning landmarks
landmarks = selectedPoint;
smoothLandmarks = selectedPointSmooth;
end
