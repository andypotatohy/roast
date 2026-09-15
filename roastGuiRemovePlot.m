function roastGuiRemovePlot(panel)
% Remove a temporary embedded plot from the ROAST Launcher plot list.

if nargin < 1 || ~ishandle(panel)
    return;
end

target = getappdata(0, 'ROAST_GUI_PLOT_TARGET');
if ~isstruct(target) || ~isfield(target, 'Canvas') || ~ishandle(target.Canvas)
    if ishandle(panel)
        delete(panel);
    end
    return;
end

canvas = target.Canvas;
panels = getCanvasAppdata(canvas, 'ROAST_GUI_PLOT_PANELS', gobjects(0));
names = getCanvasAppdata(canvas, 'ROAST_GUI_PLOT_NAMES', {});
panels = panels(ishandle(panels));
keep = panels ~= panel;

if ishandle(panel)
    delete(panel);
end

panels = panels(keep);
if numel(names) >= numel(keep)
    names = names(keep);
else
    names = names(1:min(numel(names), numel(panels)));
end

setappdata(canvas, 'ROAST_GUI_PLOT_PANELS', panels);
setappdata(canvas, 'ROAST_GUI_PLOT_NAMES', names);

if isfield(target, 'List') && ishandle(target.List)
    if isempty(names)
        set(target.List, 'String', {'Loading plots...'}, 'Value', 1);
    else
        set(target.List, 'String', names, 'Value', min(get(target.List, 'Value'), numel(names)));
    end
end

if isfield(target, 'Title') && ishandle(target.Title)
    if isempty(names)
        set(target.Title, 'String', 'Preparing plots');
    else
        currentIndex = 1;
        if isfield(target, 'List') && ishandle(target.List)
            currentIndex = get(target.List, 'Value');
        end
        currentIndex = max(1, min(currentIndex, numel(names)));
        set(target.Title, 'String', names{currentIndex});
    end
end

if isfield(target, 'Loading') && target.Loading && isempty(names)
    parentPanel = get(canvas, 'Parent');
    loadingPanels = findobj(parentPanel, 'Tag', 'ROASTPlotLoading');
    if ~isempty(loadingPanels)
        set(canvas, 'Visible', 'off');
        if isfield(target, 'List') && ishandle(target.List)
            set(target.List, 'Visible', 'off');
        end
        if isfield(target, 'Title') && ishandle(target.Title)
            set(target.Title, 'Visible', 'off');
        end
        set(loadingPanels(ishandle(loadingPanels)), 'Visible', 'on');
    end
end

function value = getCanvasAppdata(handle, name, defaultValue)
if isappdata(handle, name)
    value = getappdata(handle, name);
else
    value = defaultValue;
end
end
end
