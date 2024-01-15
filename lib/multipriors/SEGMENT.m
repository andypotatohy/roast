function SEGMENT(subj) 

    % Change the current working directory to the specified working directory  
    str = computer('arch');
    switch str
    case 'win64'
         pythonExecutable = [strrep(pwd, '\', '/') '/lib/multipriors/multipriorsEnv/python.exe'];  % Relative path to the Python executable
    case 'glnxa64'
         pythonExecutable = [strrep(pwd, '\', '/') '/lib/multipriors/multipriorsEnvLinux/bin/python3'];  % Relative path to the Python executable
    case 'maci64'
         pythonExecutable = [strrep(pwd, '\', '/') '/lib/multipriors/multipriorsEnvMac/bin/python3'];  % Relative path to the Python executable
    otherwise
        error('Unsupported operating system!');
    end
   
    pythonScript = [strrep(pwd, '\', '/') '/lib/multipriors/SEGMENT.py'];  % Relative path to the Python script
    configFile = [strrep(pwd, '\', '/')  '/lib/multipriors/Segmentation_config.py'];  % Relative path to the config file
    subj = ['"' subj '"']; % To allow spaces in subject path

    % Construct the command to execute the Python script
    command = [pythonExecutable ' ' pythonScript ' ' configFile ' ' subj]; % aka command line 

    % Execute the Python script
    status = system(command);

    if status == 0
        disp('Python script executed successfully.');
    else
        error('Error running Python script.');
    end
end
