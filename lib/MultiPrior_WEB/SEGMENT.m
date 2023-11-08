function SEGMENT(subj) 

    % Change the current working directory to the specified working directory  
    pythonExecutable = [strrep(pwd, '\', '/') '/lib/MultiPriors_WEB/MultiPriors_env/python.exe'];  % Relative path to the Python executable
    pythonScript = [strrep(pwd, '\', '/') '/lib/MultiPriors_WEB/SEGMENT.py'];  % Relative path to the Python script
    configFile = [strrep(pwd, '\', '/')  '/lib/MultiPriors_WEB/Segmentation_config.py'];  % Relative path to the config file
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
