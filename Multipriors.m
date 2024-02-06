function Multipriors(subj)

%Installs miniconda and multipriors enviornment
str = computer('arch');

switch str
    case 'win64'
        if ~exist('multipriorsEnv', 'dir')
            system([strrep(pwd, '\', '/') '/lib/multipriors/setupWindows.bat']);
        else
            disp('Files already exist on Windows. Skipping setup.');
        end

        % Specify the path to libiomp5md.dll (duplicate)
        dllPath = [strrep(pwd, '\', '/') '/lib/multipriors/multipriorsENV/Library/bin/libiomp5md.dll'];
        if exist(dllPath, 'file')
            delete(dllPath);
            disp('Duplicate libiomp5md.dll has been deleted.');
        end

    case 'glnxa64'
        setupScript = [strrep(pwd, '\', '/') '/lib/multipriors/setupLinux.sh'];

        % Grant execute permissions to the script
        system(['chmod +x ' setupScript]);

        if ~exist('multipriorsEnvLinux', 'dir')
            system([strrep(pwd, '\', '/') '/lib/multipriors/setupLinux.sh']);
        else
            disp('Files already exist on Linux. Skipping setup.');
        end

    case 'maci64'
        setupScript = [strrep(pwd, '\', '/') '/lib/multipriors/setupMac.sh'];

        % Grant execute permissions to the script
        system(['chmod +x ' setupScript]);

        if ~exist('multipriorsEnvMac', 'dir')
            system([strrep(pwd, '\', '/') '/lib/multipriors/setupMac.sh']);
        else
            disp('Files already exist on Mac. Skipping setup.');
        end

    otherwise
        error('Unsupported operating system!');
end

[dirname,baseFilename] = fileparts(subj);
if ~exist(fullfile(dirname, [baseFilename '_masks.nii']), 'file')

    WARP_indiTPM(subj)
   
    SEGMENT(subj) 
else 
    disp('MultiPriors Segmentation File: ')
    disp([dirname, [baseFilename '_masks.nii']])
    disp('Already Exists. Skipping Segmentation...')
end
end
