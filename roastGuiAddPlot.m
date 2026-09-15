function [panel, figHandle] = roastGuiAddPlot(plotName, plotMode)
% Add a plot panel to the ROAST Launcher plot list.

target = getappdata(0, 'ROAST_GUI_PLOT_TARGET');
if ~isstruct(target) || ~isfield(target, 'Canvas') || ~ishandle(target.Canvas)
    panel = [];
    figHandle = [];
    return;
end

if nargin < 2 || isempty(plotMode)
    plotMode = 'slice';
end
if nargin < 1 || isempty(plotName)
    plotName = 'Plot';
end

canvas = target.Canvas;
figHandle = target.Figure;
plotList = target.List;
isInteractive = any(strcmpi(plotMode, {'interactive', 'landmarks'}));
isLoading = isfield(target, 'Loading') && target.Loading;
displayMode = plotMode;
if isInteractive
    displayMode = '3d';
end

panels = getCanvasAppdata(canvas, 'ROAST_GUI_PLOT_PANELS', gobjects(0));
names = getCanvasAppdata(canvas, 'ROAST_GUI_PLOT_NAMES', {});

if isInteractive && ~isLoading
    parentPanel = get(canvas, 'Parent');
    loadingPanels = findobj(parentPanel, 'Tag', 'ROASTPlotLoading');
    set(loadingPanels(ishandle(loadingPanels)), 'Visible', 'off');
    set(canvas, 'Visible', 'on');
    if ishandle(plotList)
        set(plotList, 'Visible', 'on');
    end
    if isfield(target, 'Title') && ishandle(target.Title)
        set(target.Title, 'Visible', 'on');
    end
end

if isempty(panels)
    if ~isLoading
        delete(allchild(canvas));
    end
else
    panels = panels(ishandle(panels));
    for ii = 1:numel(panels)
        set(panels(ii), 'Visible', 'off');
    end
end

panel = uipanel(canvas, 'Units', 'normalized', 'Position', [0 0 1 1], ...
    'BorderType', 'none', 'BackgroundColor', 'white', 'Visible', onOff(~isLoading));
setappdata(panel, 'ROAST_GUI_PLOT_MODE', displayMode);

panels(end+1, 1) = panel;
names{end+1, 1} = char(plotName);
setappdata(canvas, 'ROAST_GUI_PLOT_PANELS', panels);
setappdata(canvas, 'ROAST_GUI_PLOT_NAMES', names);

if ishandle(plotList) && ~isLoading
    set(plotList, 'String', names, 'Value', numel(names));
end
if isfield(target, 'Title') && ishandle(target.Title) && ~isLoading
    set(target.Title, 'String', char(plotName));
end
if ishandle(figHandle)
    rotate3d(figHandle, 'off');
end

function value = onOff(tf)
if tf
    value = 'on';
else
    value = 'off';
end
end
end

function value = getCanvasAppdata(handle, name, defaultValue)
if isappdata(handle, name)
    value = getappdata(handle, name);
else
    value = defaultValue;
end
end
