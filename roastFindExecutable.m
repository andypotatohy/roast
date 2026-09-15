function executablePath = roastFindExecutable(envVar, executableName, bundledPath, extraPaths)
% roastFindExecutable Locate a platform executable used by ROAST.
%
% Priority:
%   1. An explicit environment variable path.
%   2. A bundled executable path, when provided.
%   3. Extra candidate paths.
%   4. The user's PATH.

if nargin < 3
    bundledPath = '';
end
if nargin < 4
    extraPaths = {};
end

explicitPath = strtrim(getenv(envVar));
if ~isempty(explicitPath)
    if exist(explicitPath, 'file')
        executablePath = explicitPath;
        return;
    end
    error('%s is set to "%s", but that file does not exist.', envVar, explicitPath);
end

candidates = {};
if ~isempty(bundledPath)
    candidates{end+1} = bundledPath;
end
for i = 1:numel(extraPaths)
    candidates{end+1} = extraPaths{i};
end

for i = 1:numel(candidates)
    if exist(candidates{i}, 'file')
        executablePath = candidates{i};
        return;
    end
end

if ispc
    lookupCmd = ['where ' executableName];
else
    lookupCmd = ['command -v ' executableName];
end
[status, output] = system(lookupCmd);
if status == 0
    matches = strsplit(strtrim(output), newline);
    if ~isempty(matches) && exist(matches{1}, 'file')
        executablePath = matches{1};
        return;
    end
end

error(['Could not find "%s". Install a native version for this platform, ' ...
       'or set %s to the full executable path.'], executableName, envVar);
end
