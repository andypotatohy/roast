function roastGuiProgress(step, totalSteps, message)
%ROASTGUIPROGRESS Update ROAST Launcher progress when it is active.

if isappdata(0, 'ROAST_GUI_CANCEL_REQUESTED')
    error('ROAST run cancelled from ROAST Launcher.');
end

if ~isappdata(0, 'ROAST_GUI_PROGRESS_TARGET')
    return;
end

target = getappdata(0, 'ROAST_GUI_PROGRESS_TARGET');
if ~isstruct(target) || ~isfield(target, 'Figure') || ~ishandle(target.Figure)
    return;
end

fraction = max(0, min(1, step / totalSteps));
stepIndex = min(totalSteps, max(0, floor(step)));
stepCeiling = min(0.98, (stepIndex + 0.92) / totalSteps);
if step >= totalSteps
    stepCeiling = 1;
end
percent = round(fraction * 100);
stepsLeft = max(0, totalSteps - ceil(step));

displayFraction = fraction;
if isfield(target, 'LoadingPanel') && ishandle(target.LoadingPanel)
    current = getappdata(target.LoadingPanel, 'ROAST_GUI_PROGRESS_CURRENT');
    if isempty(current)
        current = 0;
    end
    current = max(current, fraction);
    displayFraction = current;
    setappdata(target.LoadingPanel, 'ROAST_GUI_PROGRESS_CURRENT', current);
    setappdata(target.LoadingPanel, 'ROAST_GUI_PROGRESS_TARGET', max(current, stepCeiling));
    setappdata(target.LoadingPanel, 'ROAST_GUI_PROGRESS_TOTAL_STEPS', totalSteps);
    setappdata(target.LoadingPanel, 'ROAST_GUI_PROGRESS_STEP_INDEX', ...
        min(totalSteps, max(1, ceil(step))));
    setappdata(target.LoadingPanel, 'ROAST_GUI_PROGRESS_STEPS_LEFT', stepsLeft);
end

if isfield(target, 'BarFill') && ishandle(target.BarFill)
    width = max(0.02, 0.96 * displayFraction);
    set(target.BarFill, 'Position', [0.02 0.18 width 0.64]);
end

if isfield(target, 'Title') && ishandle(target.Title)
    set(target.Title, 'String', 'Running ROAST');
end

if isfield(target, 'Status') && ishandle(target.Status)
    set(target.Status, 'String', message);
end

if isfield(target, 'Progress') && ishandle(target.Progress)
    set(target.Progress, 'String', ...
        sprintf('%d%% complete - step %d of %d - %d step(s) left', ...
        percent, min(totalSteps, max(1, ceil(step))), totalSteps, stepsLeft));
end

if isfield(target, 'List') && ishandle(target.List)
    set(target.List, 'String', {sprintf('%d%% complete', percent); message}, 'Value', 1);
end

drawnow limitrate;

if isappdata(0, 'ROAST_GUI_CANCEL_REQUESTED')
    error('ROAST run cancelled from ROAST Launcher.');
end
end
