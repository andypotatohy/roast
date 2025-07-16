function runMultiaxial(T1)
% runMultiaxial(T1)
% 
% This calls Multiaxial segmentation that is based on a deep CNN.
% 
% See Hirsch et al 2021 (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8209627/)
% for details.
%
% (c) Andrew Birnbaum, Parra Lab at CCNY
%     Yu (Andy) Huang
% April 2024

% Install miniconda and Multiaxial enviornment if it doesn't exist
% And prepare to run Multiaxial python script
str = computer('arch');
switch str
    case 'win64'
        if ~exist('multiaxialEnv', 'dir')
            system([pwd '\lib\multiaxial\setupWindows.bat']);
        else
            disp('Enviornment already exist on Windows. Skipping setup...');
        end

        % Specify the path to libiomp5md.dll (duplicate)
        dllPath = [pwd '\lib\multiaxial\multiaxialENV\Library\bin\libiomp5md.dll'];
        if exist(dllPath, 'file')
            delete(dllPath);
           % disp('Duplicate libiomp5md.dll has been deleted');
        end

        pythonExecutable = [pwd '\lib\multiaxial\multiaxialEnv\python.exe'];

    case 'glnxa64'
        setupScript = [pwd '/lib/multiaxial/setupLinux.sh'];

        % Grant execute permissions to the script
        system(['chmod +x ' setupScript]);

        if ~exist('multiaxialEnvLinux', 'dir')
            system([pwd '/lib/multiaxial/setupLinux.sh']);
        else
            disp('Enviornment already exist on Linux. Skipping setup...');
        end

        pythonExecutable = [pwd '/lib/multiaxial/multiaxialEnvLinux/bin/python3'];

    case 'maci64'
        setupScript = [pwd '/lib/multiaxial/setupMac.sh'];

        % Grant execute permissions to the script
        system(['chmod +x ' setupScript]);

        if ~exist('multiaxialEnvMac', 'dir')
            system([pwd '/lib/multiaxial/setupMac.sh']);
        else
            disp('Enviornment already exist on Mac. Skipping setup...');
        end

        pythonExecutable = [pwd '/lib/multiaxial/multiaxialEnvMac/bin/python3'];

    otherwise
        error('Unsupported operating system!');
end

pythonScript = [pwd '/lib/multiaxial/SEGMENT.py'];

T1 = ['"' T1 '"']; % To allow spaces in subject path

% Construct the command to execute the Python script
command = [pythonExecutable ' ' pythonScript ' ' T1];

% Execute the Python script
status = system(command);

if status == 0
    disp('Python script executed successfully.');
else
    error('Error running Python script.');
end