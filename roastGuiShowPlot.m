function roastGuiShowPlot(panel)
%ROASTGUISHOWPLOT Reveal a finished embedded plot in the ROAST Launcher.

if nargin < 1 || ~ishandle(panel)
    return;
end

target = getappdata(0, 'ROAST_GUI_PLOT_TARGET');
if ~isstruct(target) || ~isfield(target, 'Canvas') || ~ishandle(target.Canvas)
    return;
end

canvas = target.Canvas;
panels = getCanvasAppdata(canvas, 'ROAST_GUI_PLOT_PANELS', gobjects(0));
names = getCanvasAppdata(canvas, 'ROAST_GUI_PLOT_NAMES', {});
panels = panels(ishandle(panels));
[found, index] = ismember(panel, panels);
if ~found
    return;
end

parentPanel = get(canvas, 'Parent');
loadingPanels = findobj(parentPanel, 'Tag', 'ROASTPlotLoading');
set(loadingPanels(ishandle(loadingPanels)), 'Visible', 'off');
set(panel, 'BackgroundColor', [1 1 1]);
set(canvas, 'Visible', 'on');

for ii = 1:numel(panels)
    set(panels(ii), 'Visible', onOff(ii == index));
end

if isfield(target, 'List') && ishandle(target.List)
    set(target.List, 'Visible', 'on', 'String', names, 'Value', index);
end

if isfield(target, 'Title') && ishandle(target.Title)
    set(target.Title, 'Visible', 'on');
    if numel(names) >= index
        set(target.Title, 'String', names{index});
    end
end

drawnow;
end

function value = getCanvasAppdata(handle, name, defaultValue)
if isappdata(handle, name)
    value = getappdata(handle, name);
else
    value = defaultValue;
end
end

function value = onOff(tf)
if tf
    value = 'on';
else
    value = 'off';
end
end
