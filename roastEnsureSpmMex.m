function roastEnsureSpmMex()
% roastEnsureSpmMex Ensure bundled SPM has native Apple Silicon segmentation MEX files.

if ~strcmp(computer('arch'), 'maca64')
    return;
end

roastDir = fileparts(which('roast'));
spmDir = fullfile(roastDir, 'lib', 'spm12');
mexSuffix = ['.' mexext];

requiredDirs = { ...
    spmDir, ...
    fullfile(spmDir, '@file_array', 'private') ...
};

missing = {};
for i = 1:numel(requiredDirs)
    intelMex = dir(fullfile(requiredDirs{i}, '*.mexmaci64'));
    for j = 1:numel(intelMex)
        [~, mexName] = fileparts(intelMex(j).name);
        targetPath = fullfile(requiredDirs{i}, [mexName mexSuffix]);
        if ~exist(targetPath, 'file')
            missing{end+1} = targetPath; %#ok<AGROW>
        end
    end
end

if isempty(missing)
    return;
end

sourceDir = strtrim(getenv('ROAST_SPM12_MEX_DIR'));
if isempty(sourceDir)
    sourceDir = strtrim(getenv('ROAST_SPM12_DIR'));
end
if ~isempty(sourceDir)
    copyMissingMex(sourceDir, spmDir, missing, mexSuffix);
    rehash toolboxcache;
    missing = missingExisting(missing);
    if isempty(missing)
        return;
    end
end

fprintf('\n');
fprintf('ROAST is running in native Apple Silicon MATLAB (%s), but bundled SPM12\n', computer('arch'));
fprintf('does not include Apple Silicon MEX files (%s).\n\n', mexSuffix);
fprintf('First missing file:\n  %s\n\n', missing{1});
fprintf('This ROAST package is expected to include Apple Silicon SPM MEX files. The\n');
fprintf('installation appears incomplete or the files were removed.\n\n');
fprintf('Fix options:\n');
fprintf('  1. Restore the bundled lib/spm12/*.mexmaca64 files from a complete ROAST package.\n');
fprintf('  2. Install or download an SPM12 build that contains *.mexmaca64 files, then set:\n');
fprintf('       setenv(''ROAST_SPM12_MEX_DIR'', ''/path/to/spm12'')\n');
fprintf('     before running roast. ROAST will copy matching MEX files into lib/spm12.\n');
fprintf('  3. Run Intel MATLAB under Rosetta, which can use the bundled *.mexmaci64 files.\n');
fprintf('  4. Use Multiaxial segmentation if you do not need SPM segmentation:\n');
fprintf('       roast(''example/subject1.nii'', [], ''multiaxial'', ''on'')\n\n');
error('Missing Apple Silicon SPM12 MEX files.');
end

function copyMissingMex(sourceDir, spmDir, missing, mexSuffix)
for i = 1:numel(missing)
    [targetDir, targetName] = fileparts(missing{i});
    relDir = strrep(targetDir, spmDir, '');
    if ~isempty(relDir) && relDir(1) == filesep
        relDir = relDir(2:end);
    end

    candidates = { ...
        fullfile(sourceDir, relDir, [targetName mexSuffix]), ...
        fullfile(sourceDir, [targetName mexSuffix]) ...
    };

    for j = 1:numel(candidates)
        if exist(candidates{j}, 'file')
            copyfile(candidates{j}, missing{i}, 'f');
            break;
        end
    end
end
end

function missing = missingExisting(pathsToCheck)
missing = {};
for i = 1:numel(pathsToCheck)
    if ~exist(pathsToCheck{i}, 'file')
        missing{end+1} = pathsToCheck{i}; %#ok<AGROW>
    end
end
end
