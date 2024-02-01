function Multipriors(subj)

str = computer('arch');

switch str
    case 'win64'
        if ~exist('multipriorsEnv', 'dir')
            system('setupWindows.bat');
        else
            disp('Files already exist on Windows. Skipping setup.');
        end
    case 'glnxa64'
        if ~exist('multipriorsEnvLinux', 'dir')
            system('setupLinux.sh');
        else
            disp('Files already exist on Linux. Skipping setup.');
        end
    case 'maci64'
        if ~exist('multipriorsEnvMac', 'dir')
            system('setupMac.sh');
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
