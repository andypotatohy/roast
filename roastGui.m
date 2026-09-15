function roastGui()
% roastGui()
%
% Graphical launcher for ROAST. The GUI builds and runs a normal roast(...)
% command so command-line and GUI behavior remain equivalent.

rootDir = fileparts(which(mfilename));
if ~contains(path, [rootDir filesep 'lib' filesep])
    addpath(genpath([rootDir filesep 'lib']));
end

ui = struct();
ui.roastRunning = false;
colors = warmLabTheme();

figW = 1400;
figH = 960;
fig = figure('Name', 'ROAST Launcher', ...
    'NumberTitle', 'off', ...
    'MenuBar', 'none', ...
    'ToolBar', 'none', ...
    'Units', 'pixels', ...
    'Position', centerFigure(figW, figH), ...
    'Resize', 'on', ...
    'Color', colors.window, ...
    'CloseRequestFcn', @closeGui);

main = uipanel(fig, 'Units', 'pixels', 'Position', [20 18 figW-40 figH-36], ...
    'BorderType', 'none', 'BackgroundColor', colors.window);

ui.logoPanel = uipanel(main, 'Units', 'pixels', 'Position', [4 figH-128 132 78], ...
    'BorderType', 'line', 'BackgroundColor', [1 1 1], 'Tag', 'ROASTLogoPanel');
logoAx = axes('Parent', ui.logoPanel, 'Units', 'pixels', 'Position', [0 0 132 78]);
drawRoastLogo(logoAx, colors);
ui.headerTitle = addText(main, 'ROAST Launcher', [152 figH-96 260 28], 18, 'bold', colors.window, colors.text);
ui.headerSubtitle = addText(main, 'Set up a simulation, review the command, run ROAST, and browse generated plots.', ...
    [416 figH-92 680 20], 10, 'normal', colors.window, colors.muted);

ui.tabs = uitabgroup(main, 'Units', 'pixels', 'Position', [4 48 figW-48 figH-200]);
ui.setupTab = uitab(ui.tabs, 'Title', 'Setup');
ui.plotsTab = uitab(ui.tabs, 'Title', 'Plots');
set(ui.setupTab, 'BackgroundColor', colors.window);
set(ui.plotsTab, 'BackgroundColor', colors.window);

leftX = 12;
rightX = 778;
leftW = 740;
rightW = 520;

buildSubjectPanel(ui.setupTab, leftX, 628, leftW, 96);
buildMontagePanel(ui.setupTab, leftX, 326, leftW, 280);
buildProcessingPanel(ui.setupTab, leftX, 158, leftW, 146);
buildAdvancedPanel(ui.setupTab, rightX, 376, rightW, 348);
buildCommandPanel(ui.setupTab, rightX, 80, rightW, 274);
buildResultsPanel(ui.plotsTab, 12, 16, figW-72, figH-252);
ui.statusText = addText(main, 'Ready.', [8 16 780 22], 10, 'bold', colors.window, colors.ok);

set(ui.stimRadio, 'Value', 1);
set(ui.leadFieldRadio, 'Value', 0);
set(ui.targetRadio, 'Value', 0);
advancedChanged();
modeChanged();
applyThemeColors();
refreshCommand();

    function applyThemeColors()
        if ishandle(fig)
            set(fig, 'Color', colors.window);
        end
        if ishandle(main)
            set(main, 'BackgroundColor', colors.window);
        end
        if ishandle(logoAx)
            drawRoastLogo(logoAx, colors);
        end
        delete(findobj(main, 'Tag', 'ROASTHeaderThemeBand'));
        setIfHandle(ui.setupTab, 'BackgroundColor', colors.window);
        setIfHandle(ui.plotsTab, 'BackgroundColor', colors.window);
        tabs = findall(fig, 'Type', 'uitab');
        for ii = 1:numel(tabs)
            if isequal(tabs(ii), ui.setupTab) || isequal(tabs(ii), ui.plotsTab)
                setIfHandle(tabs(ii), 'BackgroundColor', colors.window);
            else
                setIfHandle(tabs(ii), 'BackgroundColor', colors.panelAlt);
            end
        end
        setIfHandle(ui.headerTitle, 'BackgroundColor', colors.header, 'ForegroundColor', colors.headerText);
        setIfHandle(ui.headerSubtitle, 'BackgroundColor', colors.header, 'ForegroundColor', colors.headerMuted);
        setIfHandle(ui.statusText, 'BackgroundColor', colors.window);

        panels = findall(fig, 'Type', 'uipanel');
        for ii = 1:numel(panels)
            if strcmp(get(panels(ii), 'Tag'), 'ROASTLogoPanel')
                setIfHandle(panels(ii), 'BackgroundColor', [1 1 1]);
            elseif strcmp(get(panels(ii), 'BorderType'), 'none')
                setIfHandle(panels(ii), 'BackgroundColor', colors.window, 'ForegroundColor', colors.text);
            else
                setIfHandle(panels(ii), 'BackgroundColor', colors.panel, 'ForegroundColor', colors.text);
            end
        end
        headerBand = uipanel(main, 'Units', 'pixels', 'Position', headerBandPosition(), ...
            'BorderType', 'none', 'BackgroundColor', colors.header, 'Tag', 'ROASTHeaderThemeBand');
        layoutHeaderContent();
        uistack(headerBand, 'bottom');
        setIfHandle(ui.logoPanel, 'BackgroundColor', [1 1 1]);
        uistack(ui.logoPanel, 'top');
        setIfHandle(ui.headerTitle, 'BackgroundColor', colors.header, 'ForegroundColor', colors.headerText);
        setIfHandle(ui.headerSubtitle, 'BackgroundColor', colors.header, 'ForegroundColor', colors.headerMuted);
        uistack(ui.headerTitle, 'top');
        uistack(ui.headerSubtitle, 'top');
        setIfHandle(ui.embeddedPlotCanvas, 'BackgroundColor', [1 1 1]);
        plotPanels = getPlotCanvasAppdata('ROAST_GUI_PLOT_PANELS', gobjects(0));
        plotPanels = plotPanels(ishandle(plotPanels));
        for ii = 1:numel(plotPanels)
            setIfHandle(plotPanels(ii), 'BackgroundColor', [1 1 1]);
        end

        controls = findall(fig, 'Type', 'uicontrol');
        for ii = 1:numel(controls)
            style = get(controls(ii), 'Style');
            switch lower(style)
                case 'text'
                    parentBg = parentBackground(controls(ii), colors.panel);
                    setIfHandle(controls(ii), 'BackgroundColor', parentBg, 'ForegroundColor', colors.label);
                case {'edit', 'listbox', 'popupmenu'}
                    setIfHandle(controls(ii), 'BackgroundColor', colors.input, 'ForegroundColor', colors.text);
                    if strcmp(get(controls(ii), 'Tag'), 'ROASTPlotList')
                        setIfHandle(controls(ii), 'BackgroundColor', colors.plotListBg, ...
                            'ForegroundColor', colors.text, 'FontSize', 10);
                    end
                case {'checkbox', 'radiobutton'}
                    setIfHandle(controls(ii), 'BackgroundColor', parentBackground(controls(ii), colors.panel), ...
                        'ForegroundColor', colors.text);
                case 'pushbutton'
                    setIfHandle(controls(ii), 'BackgroundColor', colors.button, 'ForegroundColor', colors.text);
            end
        end
        setIfHandle(ui.headerTitle, 'BackgroundColor', colors.header, 'ForegroundColor', colors.headerText);
        setIfHandle(ui.headerSubtitle, 'BackgroundColor', colors.header, 'ForegroundColor', colors.headerMuted);
        setIfHandle(ui.currentText, 'ForegroundColor', colors.ok);
        setIfHandle(ui.modeText, 'ForegroundColor', colors.muted);

        tables = findall(fig, 'Type', 'uitable');
        for ii = 1:numel(tables)
            setIfHandle(tables(ii), 'BackgroundColor', [colors.input; colors.inputAlt], ...
                'ForegroundColor', colors.text);
        end
        drawnow limitrate;
    end

    function pos = headerBandPosition()
        tabPos = get(ui.tabs, 'Position');
        mainPos = get(main, 'Position');
        y = tabPos(2) + tabPos(4);
        h = max(116, mainPos(4) - y);
        borderAllowance = 4;
        pos = [tabPos(1) y mainPos(3) - tabPos(1) * 2 + borderAllowance h];
    end

    function layoutHeaderContent()
        headerPos = headerBandPosition();
        logoW = 132;
        logoH = 78;
        logoX = headerPos(1) + 8;
        logoY = headerPos(2) + round((headerPos(4) - logoH) / 2);
        setIfHandle(ui.logoPanel, 'Position', [logoX logoY logoW logoH]);
        if ishandle(logoAx)
            set(logoAx, 'Position', [0 0 logoW logoH]);
        end

        textY = headerPos(2) + round(headerPos(4) / 2);
        setIfHandle(ui.headerTitle, 'Position', [logoX + logoW + 20 textY - 12 280 32]);
        setIfHandle(ui.headerSubtitle, 'Position', [logoX + logoW + 310 textY - 6 700 24]);
    end

    function buildSubjectPanel(parent, x, y, w, h)
        panel = addPanel(parent, '1  Subject', [x y w h], colors.panel);
        addText(panel, 'MRI', [18 50 88 20], 10, 'bold', colors.panel, [0.12 0.12 0.12]);
        ui.subjectEdit = uicontrol(panel, 'Style', 'edit', 'String', '', ...
            'HorizontalAlignment', 'left', 'Units', 'pixels', 'Position', [94 48 w-304 26], ...
            'BackgroundColor', [1 1 1], 'Callback', @refreshCommand);
        uicontrol(panel, 'Style', 'pushbutton', 'String', 'Browse...', ...
            'Units', 'pixels', 'Position', [w-198 47 84 28], 'Callback', @browseSubject);
        uicontrol(panel, 'Style', 'pushbutton', 'String', 'MNI', ...
            'Units', 'pixels', 'Position', [w-106 47 42 28], 'Callback', @useDefaultMni);
        uicontrol(panel, 'Style', 'pushbutton', 'String', 'NY', ...
            'Units', 'pixels', 'Position', [w-56 47 36 28], 'Callback', @useNyHead);

        addText(panel, 'T2', [18 16 88 20], 10, 'bold', colors.panel, [0.12 0.12 0.12]);
        ui.t2Edit = uicontrol(panel, 'Style', 'edit', 'String', '', ...
            'HorizontalAlignment', 'left', 'Units', 'pixels', 'Position', [94 14 w-304 26], ...
            'BackgroundColor', [1 1 1], 'Callback', @refreshCommand);
        uicontrol(panel, 'Style', 'pushbutton', 'String', 'Browse...', ...
            'Units', 'pixels', 'Position', [w-198 13 84 28], 'Callback', @browseT2);
        uicontrol(panel, 'Style', 'pushbutton', 'String', 'Clear', ...
            'Units', 'pixels', 'Position', [w-106 13 86 28], 'Callback', @clearT2);
    end

    function buildMontagePanel(parent, x, y, w, h)
        panel = addPanel(parent, '2  Montage', [x y w h], colors.panel);
        ui.modeGroup = uibuttongroup(panel, 'Units', 'pixels', 'Position', [18 h-60 430 38], ...
            'BorderType', 'none', 'BackgroundColor', colors.panel, ...
            'SelectionChangedFcn', @modeChanged);
        ui.stimRadio = uicontrol(ui.modeGroup, 'Style', 'radiobutton', 'String', 'Stimulation', ...
            'Units', 'pixels', 'Position', [0 8 110 24], 'BackgroundColor', colors.panel);
        ui.leadFieldRadio = uicontrol(ui.modeGroup, 'Style', 'radiobutton', 'String', 'Lead field', ...
            'Units', 'pixels', 'Position', [128 8 108 24], 'BackgroundColor', colors.panel);
        ui.targetRadio = uicontrol(ui.modeGroup, 'Style', 'radiobutton', 'String', 'Targeting', ...
            'Units', 'pixels', 'Position', [252 8 108 24], 'BackgroundColor', colors.panel);

        ui.recipeTable = uitable(panel, 'Units', 'pixels', 'Position', [18 56 w-210 h-124], ...
            'Data', {'Fp1', 1; 'P4', -1}, ...
            'ColumnName', {'Electrode', 'mA'}, ...
            'ColumnEditable', [true true], ...
            'ColumnFormat', {'char', 'numeric'}, ...
            'ColumnWidth', {300 100}, ...
            'RowName', [], ...
            'CellEditCallback', @refreshCommand);

        uicontrol(panel, 'Style', 'pushbutton', 'String', 'Add', ...
            'Units', 'pixels', 'Position', [w-158 h-128 96 32], 'Callback', @addRecipeRow);
        uicontrol(panel, 'Style', 'pushbutton', 'String', 'Remove', ...
            'Units', 'pixels', 'Position', [w-158 h-170 96 32], 'Callback', @removeRecipeRow);
        uicontrol(panel, 'Style', 'pushbutton', 'String', 'Default', ...
            'Units', 'pixels', 'Position', [w-158 h-212 96 32], 'Callback', @defaultRecipe);

        ui.currentText = addText(panel, 'Total current: 0 mA', [18 17 210 20], ...
            10, 'bold', colors.panel, colors.ok);
        ui.modeText = addText(panel, '', [232 17 w-260 20], 10, 'normal', colors.panel, colors.muted);
    end

    function buildProcessingPanel(parent, x, y, w, h)
        panel = addPanel(parent, '3  Processing', [x y w h], colors.panel);
        addText(panel, 'Segmentation', [18 82 104 20], 10, 'bold', colors.panel, [0.12 0.12 0.12]);
        ui.segGroup = uibuttongroup(panel, 'Units', 'pixels', 'Position', [126 72 228 36], ...
            'BorderType', 'none', 'BackgroundColor', colors.panel, ...
            'SelectionChangedFcn', @refreshCommand);
        ui.spmRadio = uicontrol(ui.segGroup, 'Style', 'radiobutton', 'String', 'SPM', ...
            'Units', 'pixels', 'Position', [0 8 76 22], 'Value', 1, 'BackgroundColor', colors.panel);
        ui.multiaxialRadio = uicontrol(ui.segGroup, 'Style', 'radiobutton', 'String', 'Multiaxial', ...
            'Units', 'pixels', 'Position', [82 8 112 22], 'BackgroundColor', colors.panel);

        ui.manualGuiCheck = uicontrol(panel, 'Style', 'checkbox', 'String', 'Manual landmarks', ...
            'Units', 'pixels', 'Position', [398 80 150 24], 'BackgroundColor', colors.panel, ...
            'Callback', @refreshCommand);
        ui.resamplingCheck = uicontrol(panel, 'Style', 'checkbox', 'String', 'Resample to 1 mm', ...
            'Units', 'pixels', 'Position', [18 34 150 24], 'BackgroundColor', colors.panel, ...
            'Callback', @refreshCommand);
        ui.zeroPaddingLabel = addText(panel, 'Zero padding', [210 36 92 20], 10, 'normal', colors.panel, [0.12 0.12 0.12]);
        ui.zeroPaddingEdit = uicontrol(panel, 'Style', 'edit', 'String', '', ...
            'HorizontalAlignment', 'left', 'Units', 'pixels', 'Position', [306 34 66 25], ...
            'BackgroundColor', [1 1 1], 'Callback', @refreshCommand);
        ui.simTagLabel = addText(panel, 'Simulation tag', [398 36 100 20], 10, 'normal', colors.panel, [0.12 0.12 0.12]);
        ui.simTagEdit = uicontrol(panel, 'Style', 'edit', 'String', '', ...
            'HorizontalAlignment', 'left', 'Units', 'pixels', 'Position', [500 34 w-528 25], ...
            'BackgroundColor', [1 1 1], 'Callback', @refreshCommand);
    end

    function buildAdvancedPanel(parent, x, y, w, h)
        panel = addPanel(parent, 'Advanced Options', [x y w h], colors.panel);
        tabs = uitabgroup(panel, 'Units', 'pixels', 'Position', [12 14 w-24 h-42]);
        ui.advancedTabs = tabs;

        elecTab = uitab(tabs, 'Title', 'Electrodes');
        set(elecTab, 'BackgroundColor', colors.panelAlt);
        addText(elecTab, 'Cap', [18 244 74 20], 10, 'bold', colors.panelAlt, [0.12 0.12 0.12]);
        ui.capPopup = uicontrol(elecTab, 'Style', 'popupmenu', ...
            'String', {'1010', '1020', '1005', 'BioSemi', 'EGI', 'custom file'}, ...
            'Units', 'pixels', 'Position', [98 242 180 26], 'Callback', @capPopupChanged);
        ui.capCustomEdit = uicontrol(elecTab, 'Style', 'edit', 'String', '', ...
            'HorizontalAlignment', 'left', 'Units', 'pixels', 'Position', [18 208 w-60 26], ...
            'Enable', 'off', 'BackgroundColor', [1 1 1], 'Callback', @refreshCommand);

        addText(elecTab, 'Type', [18 168 74 20], 10, 'bold', colors.panelAlt, [0.12 0.12 0.12]);
        ui.elecTypePopup = uicontrol(elecTab, 'Style', 'popupmenu', ...
            'String', {'disc', 'pad', 'ring', 'custom MATLAB'}, ...
            'Units', 'pixels', 'Position', [98 166 180 26], 'Callback', @elecTypePopupChanged);
        ui.elecTypeCustomEdit = uicontrol(elecTab, 'Style', 'edit', 'String', '', ...
            'HorizontalAlignment', 'left', 'Units', 'pixels', 'Position', [18 132 w-60 26], ...
            'Enable', 'off', 'BackgroundColor', [1 1 1], 'Callback', @refreshCommand);

        addText(elecTab, 'Size', [18 92 74 20], 10, 'bold', colors.panelAlt, [0.12 0.12 0.12]);
        ui.elecSizeEdit = uicontrol(elecTab, 'Style', 'edit', 'String', '', ...
            'HorizontalAlignment', 'left', 'Units', 'pixels', 'Position', [98 90 w-140 26], ...
            'BackgroundColor', [1 1 1], 'Callback', @refreshCommand);
        addText(elecTab, 'Orientation', [18 54 74 20], 10, 'bold', colors.panelAlt, [0.12 0.12 0.12]);
        ui.elecOriEdit = uicontrol(elecTab, 'Style', 'edit', 'String', '', ...
            'HorizontalAlignment', 'left', 'Units', 'pixels', 'Position', [98 52 w-140 26], ...
            'BackgroundColor', [1 1 1], 'Callback', @refreshCommand);

        meshTab = uitab(tabs, 'Title', 'Mesh');
        set(meshTab, 'BackgroundColor', colors.panelAlt);
        ui.meshCheck = uicontrol(meshTab, 'Style', 'checkbox', 'String', 'Custom mesh options', ...
            'Units', 'pixels', 'Position', [18 246 190 24], 'BackgroundColor', colors.panelAlt, ...
            'Callback', @advancedChanged);
        meshNames = {'radbound', 'angbound', 'distbound', 'reratio', 'maxvol'};
        meshDefaults = {'5', '30', '0.3', '3', '10'};
        ui.meshEdits = gobjects(1, numel(meshNames));
        for ii = 1:numel(meshNames)
            y0 = 222 - (ii-1)*38;
            addText(meshTab, meshNames{ii}, [18 y0 90 20], 10, 'normal', colors.panelAlt, [0.12 0.12 0.12]);
            ui.meshEdits(ii) = uicontrol(meshTab, 'Style', 'edit', 'String', meshDefaults{ii}, ...
                'HorizontalAlignment', 'left', 'Units', 'pixels', 'Position', [118 y0-2 90 25], ...
                'BackgroundColor', [1 1 1], 'Callback', @refreshCommand);
        end

        condTab = uitab(tabs, 'Title', 'Conductivity');
        set(condTab, 'BackgroundColor', colors.panelAlt);
        ui.conductivityCheck = uicontrol(condTab, 'Style', 'checkbox', 'String', 'Custom conductivities', ...
            'Units', 'pixels', 'Position', [18 246 190 24], 'BackgroundColor', colors.panelAlt, ...
            'Callback', @advancedChanged);
        condNames = {'white', 'gray', 'csf', 'bone', 'skin', 'air', 'gel', 'electrode'};
        condDefaults = {'0.126', '0.276', '1.65', '0.01', '0.465', '2.5e-14', '0.3', '5.9e7'};
        ui.condEdits = gobjects(1, numel(condNames));
        for ii = 1:numel(condNames)
            col = floor((ii-1)/4);
            row = mod(ii-1, 4);
            x0 = 18 + col*160;
            y0 = 214 - row*42;
            addText(condTab, condNames{ii}, [x0 y0 82 20], 10, 'normal', colors.panelAlt, [0.12 0.12 0.12]);
            ui.condEdits(ii) = uicontrol(condTab, 'Style', 'edit', 'String', condDefaults{ii}, ...
                'HorizontalAlignment', 'left', 'Units', 'pixels', 'Position', [x0 y0-24 118 25], ...
                'BackgroundColor', [1 1 1], 'Callback', @refreshCommand);
        end

        extraTab = uitab(tabs, 'Title', 'Extra');
        set(extraTab, 'BackgroundColor', colors.panelAlt);
        addText(extraTab, 'Additional name-value pairs', [18 244 220 20], ...
            10, 'bold', colors.panelAlt, [0.12 0.12 0.12]);
        ui.extraOptionsEdit = uicontrol(extraTab, 'Style', 'edit', 'Max', 3, 'Min', 0, ...
            'String', '', 'HorizontalAlignment', 'left', 'Units', 'pixels', ...
            'Position', [18 64 w-60 170], 'BackgroundColor', [1 1 1], ...
            'Callback', @refreshCommand);

        targetTab = uitab(tabs, 'Title', 'Targeting');
        ui.targetTab = targetTab;
        set(targetTab, 'BackgroundColor', colors.panelAlt);
        ui.targetControls = gobjects(0);
        ui.targetControls(end+1) = addText(targetTab, 'Lead-field tag', [18 244 104 20], ...
            10, 'bold', colors.panelAlt, [0.12 0.12 0.12]);
        ui.leadFieldTagEdit = uicontrol(targetTab, 'Style', 'edit', 'String', '', ...
            'HorizontalAlignment', 'left', 'Units', 'pixels', 'Position', [130 242 w-172 26], ...
            'BackgroundColor', [1 1 1], 'Callback', @refreshCommand);
        ui.targetControls(end+1) = ui.leadFieldTagEdit;

        ui.targetControls(end+1) = addText(targetTab, 'Target coord', [18 208 104 20], ...
            10, 'bold', colors.panelAlt, [0.12 0.12 0.12]);
        ui.targetCoordEdit = uicontrol(targetTab, 'Style', 'edit', 'String', '[-48 -8 50]', ...
            'HorizontalAlignment', 'left', 'Units', 'pixels', 'Position', [130 206 w-244 26], ...
            'BackgroundColor', [1 1 1], 'Callback', @refreshCommand);
        ui.coordTypePopup = uicontrol(targetTab, 'Style', 'popupmenu', 'String', {'mni', 'voxel'}, ...
            'Units', 'pixels', 'Position', [w-106 206 64 26], 'Callback', @refreshCommand);
        ui.targetControls(end+1:end+2) = [ui.targetCoordEdit ui.coordTypePopup];

        ui.targetControls(end+1) = addText(targetTab, 'Target tag', [18 172 104 20], ...
            10, 'bold', colors.panelAlt, [0.12 0.12 0.12]);
        ui.targetingTagEdit = uicontrol(targetTab, 'Style', 'edit', 'String', '', ...
            'HorizontalAlignment', 'left', 'Units', 'pixels', 'Position', [130 170 w-172 26], ...
            'BackgroundColor', [1 1 1], 'Callback', @refreshCommand);
        ui.targetControls(end+1) = ui.targetingTagEdit;

        ui.targetControls(end+1) = addText(targetTab, 'Opt type', [18 136 104 20], ...
            10, 'bold', colors.panelAlt, [0.12 0.12 0.12]);
        ui.optTypePopup = uicontrol(targetTab, 'Style', 'popupmenu', ...
            'String', {'max-l1', 'max-l1per', 'wls-l1', 'wls-l1per', 'lcmv-l1', 'lcmv-l1per', 'unconstrained-wls', 'unconstrained-lcmv'}, ...
            'Units', 'pixels', 'Position', [130 134 w-172 26], 'Callback', @refreshCommand);
        ui.targetControls(end+1) = ui.optTypePopup;

        ui.targetControls(end+1) = addText(targetTab, 'Orient', [18 100 104 20], ...
            10, 'bold', colors.panelAlt, [0.12 0.12 0.12]);
        ui.targetOrientEdit = uicontrol(targetTab, 'Style', 'edit', 'String', '', ...
            'HorizontalAlignment', 'left', 'Units', 'pixels', 'Position', [130 98 w-172 26], ...
            'BackgroundColor', [1 1 1], 'Callback', @refreshCommand);
        ui.targetControls(end+1) = ui.targetOrientEdit;

        ui.targetControls(end+1) = addText(targetTab, 'Intensity', [18 64 70 20], ...
            10, 'bold', colors.panelAlt, [0.12 0.12 0.12]);
        ui.targetIntensityEdit = uicontrol(targetTab, 'Style', 'edit', 'String', '', ...
            'HorizontalAlignment', 'left', 'Units', 'pixels', 'Position', [88 62 54 26], ...
            'BackgroundColor', [1 1 1], 'Callback', @refreshCommand);
        ui.targetControls(end+1) = ui.targetIntensityEdit;

        ui.targetControls(end+1) = addText(targetTab, 'Elec #', [158 64 54 20], ...
            10, 'bold', colors.panelAlt, [0.12 0.12 0.12]);
        ui.elecNumEdit = uicontrol(targetTab, 'Style', 'edit', 'String', '', ...
            'HorizontalAlignment', 'left', 'Units', 'pixels', 'Position', [214 62 42 26], ...
            'BackgroundColor', [1 1 1], 'Callback', @refreshCommand);
        ui.targetControls(end+1) = ui.elecNumEdit;

        ui.targetControls(end+1) = addText(targetTab, 'Radius', [18 28 70 20], ...
            10, 'bold', colors.panelAlt, [0.12 0.12 0.12]);
        ui.targetRadiusEdit = uicontrol(targetTab, 'Style', 'edit', 'String', '', ...
            'HorizontalAlignment', 'left', 'Units', 'pixels', 'Position', [88 26 54 26], ...
            'BackgroundColor', [1 1 1], 'Callback', @refreshCommand);
        ui.targetControls(end+1) = ui.targetRadiusEdit;

        ui.targetControls(end+1) = addText(targetTab, 'k', [158 28 54 20], ...
            10, 'bold', colors.panelAlt, [0.12 0.12 0.12]);
        ui.targetKEdit = uicontrol(targetTab, 'Style', 'edit', 'String', '', ...
            'HorizontalAlignment', 'left', 'Units', 'pixels', 'Position', [214 26 42 26], ...
            'BackgroundColor', [1 1 1], 'Callback', @refreshCommand);
        ui.targetControls(end+1) = ui.targetKEdit;
    end

    function buildResultsPanel(parent, x, y, w, h)
        panel = addPanel(parent, 'ROAST Plots', [x y w h], colors.panel);
        ui.embeddedPlotPanel = panel;
        ui.resultFigures = gobjects(0);
        ui.resultIndex = 0;
        ui.embeddedPlotIndex = 0;
        ui.embeddedPlotList = uicontrol(panel, 'Style', 'listbox', 'String', {'Plots will appear here'}, ...
            'Units', 'pixels', 'Position', [18 18 230 h-56], ...
            'BackgroundColor', colors.plotListBg, 'ForegroundColor', colors.text, ...
            'FontSize', 10, 'Tag', 'ROASTPlotList', 'Callback', @embeddedPlotListChanged);
        ui.embeddedPlotTitle = uicontrol(panel, 'Style', 'text', 'String', '', ...
            'Units', 'pixels', 'Position', [266 h-58 w-284 24], ...
            'FontSize', 12, 'FontWeight', 'bold', 'HorizontalAlignment', 'center', ...
            'BackgroundColor', colors.panel, 'ForegroundColor', colors.text, ...
            'Tag', 'ROASTPlotTitle');
        ui.embeddedPlotCanvas = uipanel(panel, 'Units', 'pixels', 'Position', [266 18 w-284 h-82], ...
            'BorderType', 'line', 'BackgroundColor', [1 1 1], 'Tag', 'ROASTPlotCanvas');
        addText(ui.embeddedPlotCanvas, 'Run ROAST to show plots here.', [32 round((h-82)/2) 300 22], ...
            11, 'normal', [1 1 1], colors.muted);
    end

    function buildCommandPanel(parent, x, y, w, h)
        panel = addPanel(parent, '4  Review and Run', [x y w h], colors.panel);
        ui.commandEdit = uicontrol(panel, 'Style', 'edit', 'Max', 8, 'Min', 0, ...
            'String', '', 'HorizontalAlignment', 'left', 'Units', 'pixels', ...
            'Position', [14 56 w-28 h-96], 'Enable', 'on', ...
            'BackgroundColor', [1 1 1], 'FontName', fixedWidthFont(), ...
            'Tag', 'ROASTCommandEdit');
        buttonW = 92;
        buttonGap = 14;
        buttonStart = round((w - (buttonW*3 + buttonGap*2)) / 2);
        uicontrol(panel, 'Style', 'pushbutton', 'String', 'Copy', ...
            'Units', 'pixels', 'Position', [buttonStart 16 buttonW 30], 'Callback', @copyCommand);
        uicontrol(panel, 'Style', 'pushbutton', 'String', 'Run ROAST', ...
            'Units', 'pixels', 'Position', [buttonStart+buttonW+buttonGap 16 buttonW 30], 'FontWeight', 'bold', ...
            'Callback', @runRoastFromGui);
        uicontrol(panel, 'Style', 'pushbutton', 'String', 'Close', ...
            'Units', 'pixels', 'Position', [buttonStart+2*(buttonW+buttonGap) 16 buttonW 30], 'Callback', @(~, ~) close(fig));
    end

    function browseSubject(~, ~)
        [fn, pn] = uigetfile({'*.nii;*.hdr;*.img', 'MRI files (*.nii, *.hdr, *.img)'; '*.*', 'All files'}, ...
            'Choose subject MRI');
        if isequal(fn, 0), return; end
        set(ui.subjectEdit, 'String', fullfile(pn, fn));
        refreshCommand();
    end

    function browseT2(~, ~)
        [fn, pn] = uigetfile({'*.nii;*.hdr;*.img', 'MRI files (*.nii, *.hdr, *.img)'; '*.*', 'All files'}, ...
            'Choose T2 MRI');
        if isequal(fn, 0), return; end
        set(ui.t2Edit, 'String', fullfile(pn, fn));
        refreshCommand();
    end

    function clearT2(~, ~)
        set(ui.t2Edit, 'String', '');
        refreshCommand();
    end

    function useDefaultMni(~, ~)
        set(ui.subjectEdit, 'String', '');
        refreshCommand();
    end

    function useNyHead(~, ~)
        set(ui.subjectEdit, 'String', 'nyhead');
        refreshCommand();
    end

    function modeChanged(~, ~)
        runMode = currentRunMode();
        isLeadField = strcmpi(runMode, 'Lead field');
        isTargeting = strcmpi(runMode, 'Targeting');
        if isLeadField
            set(ui.recipeTable, 'Enable', 'off');
            set(ui.modeText, 'String', 'Generates data for roast_target.');
        elseif isTargeting
            set(ui.recipeTable, 'Enable', 'off');
            set(ui.modeText, 'String', 'Runs roast_target using the lead-field tag.');
            set(ui.advancedTabs, 'SelectedTab', ui.targetTab);
        else
            set(ui.recipeTable, 'Enable', 'on');
            set(ui.modeText, 'String', 'Currents must balance to 0 mA.');
        end
        setProcessingControlsEnabled(~isTargeting);
        setTargetControlsEnabled(isTargeting);
        refreshCommand();
    end

    function addRecipeRow(~, ~)
        data = get(ui.recipeTable, 'Data');
        data(end+1, :) = {'', 0};
        set(ui.recipeTable, 'Data', data);
        refreshCommand();
    end

    function removeRecipeRow(~, ~)
        data = get(ui.recipeTable, 'Data');
        if size(data, 1) > 1
            data(end, :) = [];
            set(ui.recipeTable, 'Data', data);
        end
        refreshCommand();
    end

    function defaultRecipe(~, ~)
        set(ui.recipeTable, 'Data', {'Fp1', 1; 'P4', -1});
        refreshCommand();
    end

    function capPopupChanged(~, ~)
        items = get(ui.capPopup, 'String');
        custom = strcmp(items{get(ui.capPopup, 'Value')}, 'custom file');
        set(ui.capCustomEdit, 'Enable', onOff(custom));
        refreshCommand();
    end

    function elecTypePopupChanged(~, ~)
        items = get(ui.elecTypePopup, 'String');
        custom = strcmp(items{get(ui.elecTypePopup, 'Value')}, 'custom MATLAB');
        set(ui.elecTypeCustomEdit, 'Enable', onOff(custom));
        refreshCommand();
    end

    function advancedChanged(~, ~)
        set(ui.meshEdits, 'Enable', onOff(get(ui.meshCheck, 'Value')));
        set(ui.condEdits, 'Enable', onOff(get(ui.conductivityCheck, 'Value')));
        if isfield(ui, 'targetControls')
            setTargetControlsEnabled(strcmpi(currentRunMode(), 'Targeting'));
        end
        refreshCommand();
    end

    function refreshCommand(~, ~)
        try
            [cmd, status, isError] = buildCommand(false);
            set(ui.commandEdit, 'String', cmd);
            if isError
                set(ui.statusText, 'ForegroundColor', colors.error, 'String', status);
            else
                set(ui.statusText, 'ForegroundColor', colors.ok, 'String', status);
            end
        catch ME
            set(ui.commandEdit, 'String', '');
            set(ui.statusText, 'ForegroundColor', colors.error, 'String', ME.message);
        end
    end

    function copyCommand(~, ~)
        cmd = getCommandString();
        if isempty(cmd), return; end
        clipboard('copy', cmd);
        set(ui.statusText, 'ForegroundColor', colors.ok, 'String', 'Command copied.');
    end

    function runRoastFromGui(~, ~)
        plotTimer = [];
        try
            ui.roastRunning = true;
            rmappdataIfPresent('ROAST_GUI_CANCEL_REQUESTED');
            cmd = strtrim(getCommandString());
            if isempty(cmd)
                [cmd, ~, ~] = buildCommand(true);
                set(ui.commandEdit, 'String', cmd);
            end
            set(ui.statusText, 'ForegroundColor', colors.muted, 'String', 'Running ROAST. See Command Window for progress.');
            ui.resultFigures = gobjects(0);
            ui.resultIndex = 0;
            clearEmbeddedPlotTabs();
            [plotTimer, progressTarget] = showPlotLoading();
            set(ui.tabs, 'SelectedTab', ui.plotsTab);
            drawnow;
            beforeFigures = findall(0, 'Type', 'figure');
            setappdata(0, 'ROAST_GUI_PLOT_TARGET', plotTarget());
            setappdata(0, 'ROAST_GUI_PROGRESS_TARGET', progressTarget);
            eval(cmd);
            stopAndDeleteTimer(plotTimer);
            ui.roastRunning = false;
            rmappdataIfPresent('ROAST_GUI_PLOT_TARGET');
            rmappdataIfPresent('ROAST_GUI_PROGRESS_TARGET');
            if isappdata(0, 'ROAST_GUI_CANCEL_REQUESTED') || ~isGuiOpen()
                rmappdataIfPresent('ROAST_GUI_CANCEL_REQUESTED');
                return;
            end
            afterFigures = findall(0, 'Type', 'figure');
            newFigures = setdiff(afterFigures, [beforeFigures; fig]);
            collectFigures(newFigures, true);
            finishPlotLoading();
            embeddedCount = embeddedPlotCount();
            set(ui.statusText, 'ForegroundColor', colors.ok, ...
                'String', sprintf('ROAST finished. %d plot(s) available.', embeddedCount));
        catch ME
            stopAndDeleteTimer(plotTimer);
            ui.roastRunning = false;
            rmappdataIfPresent('ROAST_GUI_PLOT_TARGET');
            rmappdataIfPresent('ROAST_GUI_PROGRESS_TARGET');
            wasCanceled = isappdata(0, 'ROAST_GUI_CANCEL_REQUESTED') || contains(ME.message, 'cancelled from ROAST Launcher');
            rmappdataIfPresent('ROAST_GUI_CANCEL_REQUESTED');
            if ~isGuiOpen()
                return;
            end
            finishPlotLoading();
            collectOpenFigures();
            if wasCanceled
                clearEmbeddedPlotTabs();
                set(ui.statusText, 'ForegroundColor', colors.muted, 'String', 'ROAST run cancelled.');
            else
                set(ui.statusText, 'ForegroundColor', colors.error, 'String', ME.message);
                errordlg(ME.message, 'ROAST Launcher');
            end
        end
    end

    function closeGui(~, ~)
        if ui.roastRunning
            setappdata(0, 'ROAST_GUI_CANCEL_REQUESTED', true);
            rmappdataIfPresent('ROAST_GUI_PLOT_TARGET');
            rmappdataIfPresent('ROAST_GUI_PROGRESS_TARGET');
        end
        if ishandle(fig)
            delete(fig);
        end
    end

    function collectOpenFigures(varargin)
        if ~isGuiOpen()
            return;
        end
        figs = setdiff(findall(0, 'Type', 'figure'), fig);
        collectFigures(figs, true);
        set(ui.statusText, 'ForegroundColor', colors.ok, ...
            'String', sprintf('%d plot(s) available.', embeddedPlotCount()));
    end

    function collectRunFigures(beforeFigures)
        if ~isGuiOpen()
            return;
        end
        figs = setdiff(findall(0, 'Type', 'figure'), [beforeFigures; fig]);
        collectFigures(figs, false);
        if ~isempty(ui.resultFigures)
            set(ui.statusText, 'ForegroundColor', colors.muted, ...
                'String', 'Running ROAST. Preparing plots...');
        end
        drawnow limitrate;
    end

    function collectFigures(figHandles, replaceExisting)
        if nargin < 2
            replaceExisting = false;
        end
        figHandles = figHandles(ishandle(figHandles));
        if isempty(figHandles)
            if replaceExisting || isempty(ui.resultFigures)
                ui.resultFigures = gobjects(0);
                ui.resultIndex = 0;
            end
            return;
        end
        if replaceExisting
            merged = figHandles(:);
        else
            merged = [ui.resultFigures(:); figHandles(:)];
        end
        merged = merged(ishandle(merged));
        if isempty(merged)
            return;
        end
        if replaceExisting
            ui.resultFigures = uniqueGraphicsHandles(figHandles);
        else
            ui.resultFigures = uniqueGraphicsHandles(merged);
        end
        if ui.resultIndex < 1 || ui.resultIndex > numel(ui.resultFigures)
            ui.resultIndex = 1;
        end
    end

    function plotListChanged(~, ~)
        if isempty(ui.resultFigures), return; end
        ui.resultIndex = get(ui.resultList, 'Value');
        renderPlotPreview();
    end

    function previousPlot(~, ~)
        if isempty(ui.resultFigures), return; end
        ui.resultIndex = max(1, ui.resultIndex - 1);
        set(ui.resultList, 'Value', ui.resultIndex);
        renderPlotPreview();
    end

    function nextPlot(~, ~)
        if isempty(ui.resultFigures), return; end
        ui.resultIndex = min(numel(ui.resultFigures), ui.resultIndex + 1);
        set(ui.resultList, 'Value', ui.resultIndex);
        renderPlotPreview();
    end

    function focusPlot(~, ~)
        if isempty(ui.resultFigures), return; end
        target = ui.resultFigures(ui.resultIndex);
        if ishandle(target)
            figure(target);
            set(ui.statusText, 'ForegroundColor', colors.ok, 'String', 'Focused selected plot window.');
        end
    end

    function renderPlotPreview()
        if isempty(ui.resultFigures) || ui.resultIndex < 1 || ui.resultIndex > numel(ui.resultFigures)
            clearPreview('No plots collected yet.');
            return;
        end
        target = ui.resultFigures(ui.resultIndex);
        if ~ishandle(target)
            collectOpenFigures();
            return;
        end
        delete(allchild(ui.previewPanel));
        sourceAxes = findall(target, 'Type', 'axes');
        sourceAxes = sourceAxes(arrayfun(@(h) strcmp(get(h, 'Visible'), 'on'), sourceAxes));
        if isempty(sourceAxes)
            clearPreview('Open the selected plot for interaction.');
            return;
        end
        try
            ax = copyobj(sourceAxes(1), ui.previewPanel);
            set(ax, 'Units', 'normalized', 'Position', [0.07 0.11 0.86 0.80]);
            title(ax, '', 'Interpreter', 'none');
            rotate3d(fig, 'off');
            try
                enableDefaultInteractivity(ax);
            catch
            end
        catch
            clearPreview('Preview unavailable. Use Focus.');
        end
    end

    function clearPreview(message)
        delete(allchild(ui.previewPanel));
        axes('Parent', ui.previewPanel, 'Units', 'normalized', 'Position', [0.08 0.12 0.84 0.76], ...
            'XTick', [], 'YTick', [], 'Box', 'on', 'Color', [0.985 0.99 0.995]);
        addText(ui.previewPanel, message, [34 42 300 20], ...
            9, 'normal', [0.985 0.99 0.995], colors.muted);
    end

    function clearEmbeddedPlotTabs()
        if ~isGuiOpen() || ~isfield(ui, 'embeddedPlotCanvas') || ~ishandle(ui.embeddedPlotCanvas)
            return;
        end
        deleteLoadingOverlay();
        setPlotControlsVisible(true);
        panels = getPlotCanvasAppdata('ROAST_GUI_PLOT_PANELS', gobjects(0));
        if ~isempty(panels)
            panelsToDelete = panels(ishandle(panels));
            if ~isempty(panelsToDelete)
                delete(panelsToDelete);
            end
        end
        if isappdata(ui.embeddedPlotCanvas, 'ROAST_GUI_PLOT_PANELS')
            rmappdata(ui.embeddedPlotCanvas, 'ROAST_GUI_PLOT_PANELS');
        end
        if isappdata(ui.embeddedPlotCanvas, 'ROAST_GUI_PLOT_NAMES')
            rmappdata(ui.embeddedPlotCanvas, 'ROAST_GUI_PLOT_NAMES');
        end
        ui.embeddedPlotIndex = 0;
        set(ui.embeddedPlotList, 'String', {'Plots will appear here'}, 'Value', 1);
        set(ui.embeddedPlotTitle, 'String', '');
        delete(allchild(ui.embeddedPlotCanvas));
        addText(ui.embeddedPlotCanvas, 'Run ROAST to show plots here.', [32 340 300 22], ...
            11, 'normal', [1 1 1], colors.muted);
    end

    function [t, progressTarget] = showPlotLoading()
        progressTarget = struct();
        if ~isGuiOpen()
            t = [];
            return;
        end
        deleteLoadingOverlay();
        delete(allchild(ui.embeddedPlotCanvas));
        setPlotControlsVisible(false);
        set(ui.embeddedPlotTitle, 'String', '');
        loadingPanel = uipanel(ui.embeddedPlotPanel, 'Units', 'pixels', 'Position', [18 18 0 0], ...
            'BorderType', 'line', 'BackgroundColor', [0.985 0.99 0.995], 'Tag', 'ROASTPlotLoading');
        resizeLoadingOverlay(loadingPanel);
        setappdata(ui.embeddedPlotPanel, 'ROAST_GUI_LOADING_PANEL', loadingPanel);
        loadingLogoAx = axes('Parent', loadingPanel, 'Units', 'normalized', ...
            'Position', [0.455 0.63 0.09 0.11]);
        drawRoastLogo(loadingLogoAx, colors);
        uicontrol(loadingPanel, 'Style', 'text', 'Units', 'normalized', ...
            'Position', [0.28 0.56 0.44 0.06], 'String', 'Running ROAST', ...
            'FontSize', 16, 'FontWeight', 'bold', 'HorizontalAlignment', 'center', ...
            'BackgroundColor', [0.985 0.99 0.995], 'ForegroundColor', colors.text);
        statusLabel = uicontrol(loadingPanel, 'Style', 'text', 'Units', 'normalized', ...
            'Position', [0.24 0.50 0.52 0.04], ...
            'String', 'Starting ROAST...', ...
            'FontSize', 11, 'HorizontalAlignment', 'center', ...
            'BackgroundColor', [0.985 0.99 0.995], 'ForegroundColor', colors.muted);
        progressLabel = uicontrol(loadingPanel, 'Style', 'text', 'Units', 'normalized', ...
            'Position', [0.24 0.405 0.52 0.035], 'String', 'Waiting for step output', ...
            'FontSize', 9, 'HorizontalAlignment', 'center', ...
            'BackgroundColor', [0.985 0.99 0.995], 'ForegroundColor', colors.muted, ...
            'Tag', 'ROASTPlotProgressText');
        uicontrol(loadingPanel, 'Style', 'text', 'Units', 'normalized', ...
            'Position', [0.18 0.33 0.64 0.04], ...
            'String', 'Inputs     Segmentation     Touchup     Electrodes     Mesh     Solve     Results', ...
            'FontSize', 9, 'HorizontalAlignment', 'center', ...
            'BackgroundColor', [0.985 0.99 0.995], 'ForegroundColor', [0.30 0.36 0.40]);
        barOuter = uipanel(loadingPanel, 'Units', 'normalized', 'Position', [0.30 0.45 0.40 0.024], ...
            'BorderType', 'line', 'BackgroundColor', [1 1 1]);
        barFill = uipanel(barOuter, 'Units', 'normalized', 'Position', [0.02 0.18 0.02 0.64], ...
            'BorderType', 'none', 'BackgroundColor', [0.34 0.48 0.62], 'Tag', 'ROASTPlotLoadingFill');
        setappdata(loadingPanel, 'ROAST_GUI_LOADING_PHASE', 0);
        setappdata(loadingPanel, 'ROAST_GUI_PROGRESS_CURRENT', 0);
        setappdata(loadingPanel, 'ROAST_GUI_PROGRESS_TARGET', 0.14);
        setappdata(loadingPanel, 'ROAST_GUI_PROGRESS_TOTAL_STEPS', 6);
        setappdata(loadingPanel, 'ROAST_GUI_PROGRESS_STEP_INDEX', 1);
        setappdata(loadingPanel, 'ROAST_GUI_PROGRESS_STEPS_LEFT', 6);
        set(ui.embeddedPlotList, 'String', {'Loading plots...'}, 'Value', 1);
        progressTarget = struct('Figure', fig, 'LoadingPanel', loadingPanel, ...
            'Title', ui.embeddedPlotTitle, 'List', ui.embeddedPlotList, ...
            'Status', statusLabel, 'Progress', progressLabel, 'BarFill', barFill);
        try
            t = timer('ExecutionMode', 'fixedSpacing', 'Period', 0.08, ...
                'BusyMode', 'drop', 'TimerFcn', @(~, ~) advancePlotLoadingBar(loadingPanel, barFill));
            start(t);
        catch
            t = [];
        end
    end

    function advancePlotLoadingBar(loadingPanel, barFill)
        if ~ishandle(loadingPanel) || ~ishandle(barFill)
            return;
        end
        phase = getappdata(loadingPanel, 'ROAST_GUI_LOADING_PHASE');
        if isempty(phase)
            phase = 0;
        end
        phase = phase + 0.04;
        if phase > 1
            phase = 0;
        end
        setappdata(loadingPanel, 'ROAST_GUI_LOADING_PHASE', phase);
        current = getappdata(loadingPanel, 'ROAST_GUI_PROGRESS_CURRENT');
        target = getappdata(loadingPanel, 'ROAST_GUI_PROGRESS_TARGET');
        if isempty(current), current = 0; end
        if isempty(target), target = min(0.96, current + 0.08); end
        current = min(target, current + max(0.0015, (target - current) * 0.025));
        setappdata(loadingPanel, 'ROAST_GUI_PROGRESS_CURRENT', current);
        set(barFill, 'Position', [0.02 0.18 max(0.02, 0.96 * current) 0.64]);
        progressLabel = findobj(loadingPanel, 'Tag', 'ROASTPlotProgressText');
        if ~isempty(progressLabel) && ishandle(progressLabel(1))
            totalSteps = getappdata(loadingPanel, 'ROAST_GUI_PROGRESS_TOTAL_STEPS');
            stepIndex = getappdata(loadingPanel, 'ROAST_GUI_PROGRESS_STEP_INDEX');
            stepsLeft = getappdata(loadingPanel, 'ROAST_GUI_PROGRESS_STEPS_LEFT');
            if isempty(totalSteps), totalSteps = 6; end
            if isempty(stepIndex), stepIndex = 1; end
            if isempty(stepsLeft), stepsLeft = max(0, totalSteps - stepIndex); end
            set(progressLabel(1), 'String', sprintf('%d%% complete - step %d of %d - %d step(s) left', ...
                round(current * 100), stepIndex, totalSteps, stepsLeft));
        end
        drawnow limitrate;
    end

    function finishPlotLoading()
        if ~isGuiOpen() || ~isfield(ui, 'embeddedPlotCanvas') || ~ishandle(ui.embeddedPlotCanvas)
            return;
        end
        panels = getPlotCanvasAppdata('ROAST_GUI_PLOT_PANELS', gobjects(0));
        panels = panels(ishandle(panels));
        names = getPlotCanvasAppdata('ROAST_GUI_PLOT_NAMES', {});
        deleteLoadingOverlay();
        setPlotControlsVisible(true);
        if isempty(panels)
            clearEmbeddedPlotTabs();
            return;
        end
        canvasChildren = allchild(ui.embeddedPlotCanvas);
        staleChildren = setdiff(canvasChildren, panels);
        delete(staleChildren(ishandle(staleChildren)));
        if isempty(names)
            names = arrayfun(@(ii) sprintf('Plot %d', ii), 1:numel(panels), 'UniformOutput', false)';
        end
        for ii = 1:numel(panels)
            set(panels(ii), 'Visible', 'off');
        end
        ui.embeddedPlotIndex = 1;
        set(ui.embeddedPlotList, 'String', names, 'Value', 1);
        selectEmbeddedPlot(1);
        drawnow;
    end

    function setPlotControlsVisible(isVisible)
        state = onOff(isVisible);
        if isfield(ui, 'embeddedPlotList') && ishandle(ui.embeddedPlotList)
            set(ui.embeddedPlotList, 'Visible', state);
        end
        if isfield(ui, 'embeddedPlotTitle') && ishandle(ui.embeddedPlotTitle)
            set(ui.embeddedPlotTitle, 'Visible', state);
        end
        if isfield(ui, 'embeddedPlotCanvas') && ishandle(ui.embeddedPlotCanvas)
            set(ui.embeddedPlotCanvas, 'Visible', state);
        end
    end

    function resizeLoadingOverlay(loadingPanel)
        if ~isfield(ui, 'embeddedPlotPanel') || ~ishandle(ui.embeddedPlotPanel) || ~ishandle(loadingPanel)
            return;
        end
        parentPos = get(ui.embeddedPlotPanel, 'Position');
        set(loadingPanel, 'Position', [18 18 max(10, parentPos(3)-36) max(10, parentPos(4)-56)]);
    end

    function deleteLoadingOverlay()
        loadingPanels = gobjects(0);
        if isfield(ui, 'embeddedPlotPanel') && ishandle(ui.embeddedPlotPanel)
            loadingPanels = [loadingPanels; findobj(ui.embeddedPlotPanel, 'Tag', 'ROASTPlotLoading')];
            if isappdata(ui.embeddedPlotPanel, 'ROAST_GUI_LOADING_PANEL')
                panel = getappdata(ui.embeddedPlotPanel, 'ROAST_GUI_LOADING_PANEL');
                if ishandle(panel)
                    loadingPanels = [loadingPanels; panel];
                end
                rmappdata(ui.embeddedPlotPanel, 'ROAST_GUI_LOADING_PANEL');
            end
        end
        if isfield(ui, 'embeddedPlotCanvas') && ishandle(ui.embeddedPlotCanvas)
            loadingPanels = [loadingPanels; findobj(ui.embeddedPlotCanvas, 'Tag', 'ROASTPlotLoading')];
        end
        loadingPanels = uniqueGraphicsHandles(loadingPanels);
        delete(loadingPanels(ishandle(loadingPanels)));
    end

    function n = embeddedPlotCount()
        if ~isGuiOpen() || ~isfield(ui, 'embeddedPlotCanvas') || ~ishandle(ui.embeddedPlotCanvas)
            n = 0;
            return;
        end
        panels = getPlotCanvasAppdata('ROAST_GUI_PLOT_PANELS', gobjects(0));
        n = numel(panels(ishandle(panels)));
    end

    function target = plotTarget()
        target = struct();
        target.Figure = fig;
        target.Canvas = ui.embeddedPlotCanvas;
        target.List = ui.embeddedPlotList;
        target.Title = ui.embeddedPlotTitle;
        target.Loading = true;
    end

    function embeddedPlotListChanged(~, ~)
        names = get(ui.embeddedPlotList, 'String');
        panels = getPlotCanvasAppdata('ROAST_GUI_PLOT_PANELS', gobjects(0));
        if isempty(panels) || (iscell(names) && isscalar(names) && strcmp(names{1}, 'Plots will appear here'))
            return;
        end
        ui.embeddedPlotIndex = get(ui.embeddedPlotList, 'Value');
        selectEmbeddedPlot(ui.embeddedPlotIndex);
    end

    function selectEmbeddedPlot(index)
        panels = getPlotCanvasAppdata('ROAST_GUI_PLOT_PANELS', gobjects(0));
        panels = panels(ishandle(panels));
        if isempty(panels)
            return;
        end
        index = max(1, min(index, numel(panels)));
        for ii = 1:numel(panels)
            set(panels(ii), 'Visible', onOff(ii == index));
        end
        ui.embeddedPlotIndex = index;
        set(ui.embeddedPlotList, 'Value', index);
        names = getPlotCanvasAppdata('ROAST_GUI_PLOT_NAMES', {});
        if numel(names) >= index
            set(ui.embeddedPlotTitle, 'String', names{index});
        end
        mode = 'slice';
        if isappdata(panels(index), 'ROAST_GUI_PLOT_MODE')
            mode = getappdata(panels(index), 'ROAST_GUI_PLOT_MODE');
        end
        if isappdata(panels(index), 'ROAST_GUI_PANEL_COLORMAP')
            panelMap = getappdata(panels(index), 'ROAST_GUI_PANEL_COLORMAP');
            selectedAxes = findall(panels(index), 'Type', 'axes');
            for ii = 1:numel(selectedAxes)
                try
                    colormap(selectedAxes(ii), panelMap);
                catch
                end
            end
        end
        rotate3d(fig, 'off');
        selectedAxes = findall(panels(index), 'Type', 'axes');
        if strcmp(mode, '3d')
            rotate3d(fig, 'on');
            for ii = 1:numel(selectedAxes)
                enableEmbedded3DInteractions(selectedAxes(ii));
            end
        else
            for ii = 1:numel(selectedAxes)
                try
                    disableDefaultInteractivity(selectedAxes(ii));
                catch
                end
                try
                    selectedAxes(ii).Interactions = [];
                catch
                end
            end
        end
    end

    function value = getPlotCanvasAppdata(name, defaultValue)
        if isGuiOpen() && isfield(ui, 'embeddedPlotCanvas') && ishandle(ui.embeddedPlotCanvas) && isappdata(ui.embeddedPlotCanvas, name)
            value = getappdata(ui.embeddedPlotCanvas, name);
        else
            value = defaultValue;
        end
    end

    function tf = isGuiOpen()
        tf = ishandle(fig);
    end

    function plotTabChanged(~, event)
        mode = 'overview';
        selectedTab = [];
        if nargin > 1 && isobject(event) && isprop(event, 'NewValue')
            selectedTab = event.NewValue;
        elseif nargin > 1 && isstruct(event) && isfield(event, 'NewValue')
            selectedTab = event.NewValue;
        elseif isfield(ui, 'plotTabs') && ishandle(ui.plotTabs)
            selectedTab = ui.plotTabs.SelectedTab;
        end
        if ~isempty(selectedTab) && isappdata(selectedTab, 'ROAST_GUI_PLOT_MODE')
            mode = getappdata(selectedTab, 'ROAST_GUI_PLOT_MODE');
        elseif isappdata(ui.plotTabs.SelectedTab, 'ROAST_GUI_PLOT_MODE')
            mode = getappdata(ui.plotTabs.SelectedTab, 'ROAST_GUI_PLOT_MODE');
        end
        rotate3d(fig, 'off');
        if strcmp(mode, '3d')
            selectedAxes = findall(selectedTab, 'Type', 'axes');
            for ii = 1:numel(selectedAxes)
                try
                    enableDefaultInteractivity(selectedAxes(ii));
                catch
                end
            end
        end
    end

    function enableEmbedded3DInteractions(ax)
        try
            ax.Interactions = [rotateInteraction zoomInteraction dataTipInteraction];
        catch
            try
                enableDefaultInteractivity(ax);
            catch
            end
        end
    end

    function stopAndDeleteTimer(t)
        if isempty(t) || ~isvalid(t)
            return;
        end
        if strcmp(get(t, 'Running'), 'on')
            stop(t);
        end
        delete(t);
    end

    function rmappdataIfPresent(name)
        if isappdata(0, name)
            rmappdata(0, name);
        end
    end

    function [cmd, status, isError] = buildCommand(validateRun)
        isError = false;
        subj = strtrim(get(ui.subjectEdit, 'String'));
        t2 = strtrim(get(ui.t2Edit, 'String'));
        runMode = currentRunMode();
        isLeadField = strcmpi(runMode, 'Lead field');
        isTargeting = strcmpi(runMode, 'Targeting');

        if isempty(subj)
            subjExpr = '[]';
        else
            subjExpr = matlabString(subj);
        end

        if isTargeting
            cmd = buildTargetCommand(subjExpr, validateRun);
            status = 'Ready.';
            return;
        end

        if isLeadField
            recipeExpr = '''leadField''';
        else
            recipeExpr = recipeToExpression(validateRun);
        end

        args = {subjExpr, recipeExpr};
        if ~isLeadField
            args = appendElectrodeOptions(args);
        end

        if ~isempty(t2)
            if strcmpi(get(get(ui.segGroup, 'SelectedObject'), 'String'), 'Multiaxial')
                msg = 'Multiaxial cannot be used with T2.';
                if validateRun, error(msg); end
                isError = true;
                status = msg;
            end
            args(end+1:end+2) = {'''T2''', matlabString(t2)};
        end

        if strcmpi(get(get(ui.segGroup, 'SelectedObject'), 'String'), 'Multiaxial')
            args(end+1:end+2) = {'''multiaxial''', '''on'''};
        end
        if get(ui.manualGuiCheck, 'Value')
            args(end+1:end+2) = {'''manualGui''', '''on'''};
        end
        if get(ui.resamplingCheck, 'Value')
            args(end+1:end+2) = {'''resampling''', '''on'''};
        end

        padding = strtrim(get(ui.zeroPaddingEdit, 'String'));
        if ~isempty(padding)
            paddingValue = str2double(padding);
            if isnan(paddingValue) || paddingValue <= 0 || mod(paddingValue, 1) ~= 0
                error('Zero padding must be a positive integer.');
            end
            args(end+1:end+2) = {'''zeroPadding''', num2str(paddingValue)};
        end

        tag = strtrim(get(ui.simTagEdit, 'String'));
        if ~isempty(tag)
            args(end+1:end+2) = {'''simulationTag''', matlabString(tag)};
        end

        if get(ui.meshCheck, 'Value')
            args(end+1:end+2) = {'''meshOptions''', meshOptionsExpression()};
        end
        if get(ui.conductivityCheck, 'Value')
            args(end+1:end+2) = {'''conductivities''', conductivityExpression()};
        end

        extra = strtrim(getMultilineString(ui.extraOptionsEdit));
        if ~isempty(extra)
            args{end+1} = extra;
        end

        cmd = ['roast(' strjoin(args, ', ') ');'];
        if ~exist('status', 'var')
            [status, isError] = currentBalanceText();
        end
    end

    function cmd = buildTargetCommand(subjExpr, validateRun)
        leadFieldTag = strtrim(get(ui.leadFieldTagEdit, 'String'));
        if isempty(leadFieldTag)
            if validateRun
                error('Targeting mode needs the lead-field simulation tag.');
            end
            simTagExpr = '''leadFieldTag''';
        else
            simTagExpr = matlabString(leadFieldTag);
        end

        targetExpr = strtrim(get(ui.targetCoordEdit, 'String'));
        if isempty(targetExpr)
            targetExpr = '[]';
        end

        coordItems = get(ui.coordTypePopup, 'String');
        coordType = coordItems{get(ui.coordTypePopup, 'Value')};
        args = {subjExpr, simTagExpr, targetExpr};
        if ~strcmpi(coordType, 'mni')
            args(end+1:end+2) = {'''coordType''', matlabString(coordType)};
        end

        optItems = get(ui.optTypePopup, 'String');
        optType = optItems{get(ui.optTypePopup, 'Value')};
        if ~strcmpi(optType, 'max-l1')
            args(end+1:end+2) = {'''optType''', matlabString(optType)};
        end

        orientExpr = strtrim(get(ui.targetOrientEdit, 'String'));
        if ~isempty(orientExpr)
            if isOptionKeyword(orientExpr)
                orientExpr = matlabString(orientExpr);
            end
            args(end+1:end+2) = {'''orient''', orientExpr};
        end

        intensityExpr = strtrim(get(ui.targetIntensityEdit, 'String'));
        if ~isempty(intensityExpr)
            args(end+1:end+2) = {'''desiredIntensity''', intensityExpr};
        end

        elecNumExpr = strtrim(get(ui.elecNumEdit, 'String'));
        if ~isempty(elecNumExpr)
            args(end+1:end+2) = {'''elecNum''', elecNumExpr};
        end

        radiusExpr = strtrim(get(ui.targetRadiusEdit, 'String'));
        if ~isempty(radiusExpr)
            args(end+1:end+2) = {'''targetRadius''', radiusExpr};
        end

        kExpr = strtrim(get(ui.targetKEdit, 'String'));
        if ~isempty(kExpr)
            args(end+1:end+2) = {'''k''', kExpr};
        end

        targetTag = strtrim(get(ui.targetingTagEdit, 'String'));
        if ~isempty(targetTag)
            args(end+1:end+2) = {'''targetingTag''', matlabString(targetTag)};
        end

        extra = strtrim(getMultilineString(ui.extraOptionsEdit));
        if ~isempty(extra)
            args{end+1} = extra;
        end

        cmd = ['roast_target(' strjoin(args, ', ') ');'];
    end

    function args = appendElectrodeOptions(args)
        capItems = get(ui.capPopup, 'String');
        capValue = capItems{get(ui.capPopup, 'Value')};
        if strcmp(capValue, 'custom file')
            capValue = strtrim(get(ui.capCustomEdit, 'String'));
            if isempty(capValue), error('Custom cap file is selected but empty.'); end
        end
        if ~strcmpi(capValue, '1010')
            args(end+1:end+2) = {'''capType''', matlabString(capValue)};
        end

        typeItems = get(ui.elecTypePopup, 'String');
        typeValue = typeItems{get(ui.elecTypePopup, 'Value')};
        if strcmp(typeValue, 'custom MATLAB')
            typeExpr = strtrim(get(ui.elecTypeCustomEdit, 'String'));
            if isempty(typeExpr), error('Custom electrode type is selected but empty.'); end
            args(end+1:end+2) = {'''elecType''', typeExpr};
        elseif ~strcmpi(typeValue, 'disc')
            args(end+1:end+2) = {'''elecType''', matlabString(typeValue)};
        end

        sizeExpr = strtrim(get(ui.elecSizeEdit, 'String'));
        if ~isempty(sizeExpr)
            args(end+1:end+2) = {'''elecSize''', sizeExpr};
        end
        oriExpr = strtrim(get(ui.elecOriEdit, 'String'));
        if ~isempty(oriExpr)
            if any(strcmpi(oriExpr, {'lr', 'ap', 'si'}))
                oriExpr = matlabString(oriExpr);
            end
            args(end+1:end+2) = {'''elecOri''', oriExpr};
        end
    end

    function recipeExpr = recipeToExpression(validateRun)
        data = get(ui.recipeTable, 'Data');
        parts = {};
        total = 0;
        for ii = 1:size(data, 1)
            name = strtrim(char(data{ii, 1}));
            current = data{ii, 2};
            if ischar(current), current = str2double(current); end
            if isempty(name)
                if validateRun, error('Recipe contains an empty electrode name.'); end
                continue;
            end
            if isempty(current) || ~isnumeric(current) || isnan(current)
                error('Recipe current for %s must be numeric.', name);
            end
            parts{end+1} = matlabString(name); %#ok<AGROW>
            parts{end+1} = sprintf('%.15g', current); %#ok<AGROW>
            total = total + current;
        end
        if isempty(parts)
            recipeExpr = '[]';
        else
            recipeExpr = ['{' strjoin(parts, ', ') '}'];
        end
        if validateRun && abs(total) > eps
            error('Recipe currents must sum to 0 mA. Current sum is %.15g mA.', total);
        end
    end

    function expr = meshOptionsExpression()
        names = {'radbound', 'angbound', 'distbound', 'reratio', 'maxvol'};
        parts = cell(1, numel(names)*2);
        for ii = 1:numel(names)
            value = strtrim(get(ui.meshEdits(ii), 'String'));
            if isempty(value), error('Mesh option %s is empty.', names{ii}); end
            parts{ii*2-1} = matlabString(names{ii});
            parts{ii*2} = value;
        end
        expr = ['struct(' strjoin(parts, ', ') ')'];
    end

    function expr = conductivityExpression()
        names = {'white', 'gray', 'csf', 'bone', 'skin', 'air', 'gel', 'electrode'};
        parts = cell(1, numel(names)*2);
        for ii = 1:numel(names)
            value = strtrim(get(ui.condEdits(ii), 'String'));
            if isempty(value), error('Conductivity %s is empty.', names{ii}); end
            parts{ii*2-1} = matlabString(names{ii});
            parts{ii*2} = value;
        end
        expr = ['struct(' strjoin(parts, ', ') ')'];
    end

    function [text, isError] = currentBalanceText()
        runMode = currentRunMode();
        isLeadField = strcmpi(runMode, 'Lead field');
        isTargeting = strcmpi(runMode, 'Targeting');
        total = 0;
        data = get(ui.recipeTable, 'Data');
        for ii = 1:size(data, 1)
            current = data{ii, 2};
            if ischar(current), current = str2double(current); end
            if isnumeric(current) && ~isnan(current)
                total = total + current;
            end
        end
        if isTargeting
            set(ui.currentText, 'ForegroundColor', colors.muted, 'String', 'Targeting mode');
            text = 'Ready.';
            isError = false;
        elseif isLeadField
            set(ui.currentText, 'ForegroundColor', colors.muted, 'String', 'Lead field mode');
            text = 'Ready.';
            isError = false;
        elseif abs(total) > eps
            set(ui.currentText, 'ForegroundColor', colors.error, ...
                'String', sprintf('Total current: %.15g mA', total));
            text = sprintf('Currents must sum to 0 mA. Current sum: %.15g mA.', total);
            isError = true;
        else
            set(ui.currentText, 'ForegroundColor', colors.ok, ...
                'String', sprintf('Total current: %.15g mA', total));
            text = 'Ready.';
            isError = false;
        end
    end

    function cmd = getCommandString()
        raw = get(ui.commandEdit, 'String');
        if iscell(raw)
            cmd = strjoin(raw, newline);
        else
            cmd = raw;
        end
    end

    function mode = currentRunMode()
        selected = get(ui.modeGroup, 'SelectedObject');
        if isempty(selected)
            mode = 'Stimulation';
        else
            mode = get(selected, 'String');
        end
    end

    function setTargetControlsEnabled(enabled)
        if ~isfield(ui, 'targetControls') || isempty(ui.targetControls)
            return;
        end
        controls = ui.targetControls(ishandle(ui.targetControls));
        set(controls, 'Enable', onOff(enabled));
    end

    function setProcessingControlsEnabled(enabled)
        controls = [ui.spmRadio ui.multiaxialRadio ui.manualGuiCheck ui.resamplingCheck ...
            ui.zeroPaddingLabel ui.zeroPaddingEdit ui.simTagLabel ui.simTagEdit];
        controls = controls(ishandle(controls));
        set(controls, 'Enable', onOff(enabled));
    end
end

function names = figureNames(figHandles)
names = cell(numel(figHandles), 1);
for ii = 1:numel(figHandles)
    if ishandle(figHandles(ii))
        name = get(figHandles(ii), 'Name');
        if isempty(name)
            name = sprintf('Figure %d', double(figHandles(ii)));
        end
        names{ii} = name;
    else
        names{ii} = sprintf('Closed plot %d', ii);
    end
end
end

function handlesOut = uniqueGraphicsHandles(handlesIn)
handlesIn = handlesIn(:);
handlesOut = gobjects(0);
for ii = 1:numel(handlesIn)
    if ~ishandle(handlesIn(ii))
        continue;
    end
    isDuplicate = false;
    for jj = 1:numel(handlesOut)
        if isequal(handlesIn(ii), handlesOut(jj))
            isDuplicate = true;
            break;
        end
    end
    if ~isDuplicate
        handlesOut(end+1, 1) = handlesIn(ii); %#ok<AGROW>
    end
end
end

function drawRoastLogo(ax, ~)
cla(ax);
hold(ax, 'on');
axis(ax, 'equal');
axis(ax, 'off');
set(ax, 'Color', [1 1 1], 'XLim', [0 180], 'YLim', [0 112]);

assetPath = fullfile(fileparts(which('roastGui')), 'roastLogo.png');
if exist(assetPath, 'file')
    [img, ~, alpha] = imread(assetPath);
    img = im2double(img);
    if ~isempty(alpha)
        alpha = im2double(alpha);
        if ismatrix(alpha)
            alpha = repmat(alpha, 1, 1, size(img, 3));
        end
        img = img .* alpha + (1 - alpha);
    end
    image(ax, [0 180], [112 0], img);
    set(ax, 'YDir', 'normal');
    rectangle(ax, 'Position', [0.5 0.5 179 111], 'EdgeColor', [0.86 0.88 0.90], 'LineWidth', 1);
    hold(ax, 'off');
    return;
end

ink = [0 0 0];
softInk = [0.34 0.34 0.34];
lw = 1.6;
thin = 0.85;

rectangle(ax, 'Position', [8 3 164 106], 'Curvature', 0.035, ...
    'FaceColor', [1 1 1], 'EdgeColor', [0.84 0.86 0.88], 'LineWidth', 1.0);

% Outer oven body and upper control box.
line(ax, [38 142 142 38 38], [14 14 80 80 14], 'Color', ink, 'LineWidth', lw);
line(ax, [54 126 126 54 54], [80 80 106 106 80], 'Color', ink, 'LineWidth', lw);

% Handle and side caps.
line(ax, [20 160], [72 72], 'Color', ink, 'LineWidth', lw);
line(ax, [20 160], [64 64], 'Color', ink, 'LineWidth', lw);
line(ax, [20 20 16 16 20], [64 72 72 64 64], 'Color', ink, 'LineWidth', lw);
line(ax, [160 160 164 164 160], [64 72 72 64 64], 'Color', ink, 'LineWidth', lw);

% Feet.
line(ax, [50 50 62 62], [14 4 4 14], 'Color', ink, 'LineWidth', lw);
line(ax, [118 118 130 130], [14 4 4 14], 'Color', ink, 'LineWidth', lw);

% Front glass and lower panel.
line(ax, [48 132 132 48 48], [28 28 57 57 28], 'Color', ink, 'LineWidth', lw);
line(ax, [48 132], [25 25], 'Color', ink, 'LineWidth', lw);
line(ax, [90 90], [14 25], 'Color', softInk, 'LineWidth', thin);

% Knobs.
theta = linspace(0, 2*pi, 120);
knobX = [70 90 110];
for ii = 1:numel(knobX)
    plot(ax, knobX(ii) + 8*cos(theta), 94 + 8*sin(theta), 'Color', ink, 'LineWidth', thin);
    rectangle(ax, 'Position', [knobX(ii)-10 91 20 6], 'EdgeColor', ink, 'LineWidth', thin);
end
line(ax, [62 78], [94 94], 'Color', ink, 'LineWidth', thin);
line(ax, [82 98], [86 102], 'Color', ink, 'LineWidth', thin);
line(ax, [102 118], [94 94], 'Color', ink, 'LineWidth', thin);

% Brain silhouette and folds.
brainX = [68 71 76 83 91 101 109 115 119 118 113 105 95 86 78 71 67 68];
brainY = [42 49 54 57 58 56 53 48 42 36 31 28 27 28 30 34 38 42];
patch(ax, brainX, brainY, [1 1 1], 'EdgeColor', ink, 'LineWidth', thin);
plot(ax, [72 78 82 88 94 100 106 112], [43 47 45 50 49 52 48 44], 'Color', softInk, 'LineWidth', 0.65);
plot(ax, [76 82 89 96 104], [37 40 39 42 39], 'Color', softInk, 'LineWidth', 0.65);
plot(ax, [82 84 88 91], [54 49 46 41], 'Color', softInk, 'LineWidth', 0.65);
plot(ax, [99 99 102 104], [54 49 45 40], 'Color', softInk, 'LineWidth', 0.65);
plot(ax, [73 78 84 89], [34 31 30 28], 'Color', softInk, 'LineWidth', 0.65);
plot(ax, [70 76 83], [38 36 36], 'Color', softInk, 'LineWidth', 0.65);

% Cerebellum hatch.
for yy = 29:2:37
    plot(ax, [71 84], [yy yy+2], 'Color', softInk, 'LineWidth', 0.55);
end

rectangle(ax, 'Position', [0.5 0.5 179 111], 'EdgeColor', [0.86 0.88 0.90], 'LineWidth', 1);
hold(ax, 'off');
end

function panel = addPanel(parent, titleText, position, bgColor)
panel = uipanel(parent, 'Title', titleText, 'Units', 'pixels', ...
    'Position', position, 'BackgroundColor', bgColor, ...
    'FontWeight', 'bold', 'BorderType', 'line');
end

function label = addText(parent, textValue, position, fontSize, fontWeight, bgColor, fgColor)
label = uicontrol(parent, 'Style', 'text', 'String', textValue, ...
    'Units', 'pixels', 'Position', position, ...
    'HorizontalAlignment', 'left', 'FontSize', fontSize, ...
    'FontWeight', fontWeight, 'BackgroundColor', bgColor, ...
    'ForegroundColor', fgColor);
end

function colors = warmLabTheme()
window = [0.98 0.94 0.88];
input = [1.00 0.995 0.98];
colors = struct('window', window, ...
    'panel', [1.00 0.985 0.955], ...
    'panelAlt', [0.96 0.91 0.84], ...
    'button', [0.94 0.86 0.76], ...
    'input', input, ...
    'inputAlt', input * 0.96 + window * 0.04, ...
    'plotBg', [1 1 1], ...
    'plotListBg', [1.00 0.975 0.925], ...
    'accent', [0.58 0.23 0.10], ...
    'text', [0.16 0.10 0.06], ...
    'brand', [0.77 0.16 0.10], ...
    'label', [0.16 0.10 0.06], ...
    'muted', [0.46 0.34 0.24], ...
    'error', [0.68 0.08 0.04], ...
    'ok', [0.22 0.42 0.18], ...
    'header', [0.60 0.28 0.12], ...
    'headerText', [1.00 0.96 0.90], ...
    'headerMuted', [0.94 0.80 0.68]);
end

function setIfHandle(handle, varargin)
if ~isempty(handle) && ishandle(handle)
    try
        set(handle, varargin{:});
    catch
    end
end
end

function bg = parentBackground(handle, defaultBg)
bg = defaultBg;
try
    parent = get(handle, 'Parent');
    if ishandle(parent) && isprop(parent, 'BackgroundColor')
        bg = get(parent, 'BackgroundColor');
    end
catch
end
end

function value = matlabString(textValue)
value = ['''' strrep(textValue, '''', '''''') ''''];
end

function value = onOff(tf)
if tf
    value = 'on';
else
    value = 'off';
end
end

function textValue = getMultilineString(control)
raw = get(control, 'String');
if iscell(raw)
    textValue = strjoin(raw, ' ');
else
    textValue = raw;
end
end

function tf = isOptionKeyword(value)
keywords = {'radial-in','radial-out','right','left','anterior','posterior', ...
    'right-anterior','right-posterior','left-anterior','left-posterior','optimal'};
tf = any(strcmpi(value, keywords));
end

function fontName = fixedWidthFont()
if ismac
    fontName = 'Menlo';
elseif ispc
    fontName = 'Consolas';
else
    fontName = 'Monospaced';
end
end

function position = centerFigure(width, height)
screenSize = get(0, 'ScreenSize');
left = max(20, round((screenSize(3) - width) / 2));
bottom = max(40, round((screenSize(4) - height) / 2));
position = [left bottom width height];
end
