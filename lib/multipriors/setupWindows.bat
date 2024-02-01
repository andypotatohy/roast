@echo off
set CONDA_INSTALLER=Miniconda3-latest-Windows-x86_64.exe
set ENV_NAME=multipriors_env
set YAML_FILE=multipriorsEnv.yml

:: Download and install Miniconda
curl -O https://repo.anaconda.com/miniconda/%CONDA_INSTALLER%
start /wait %CONDA_INSTALLER% /InstallationType=JustMe /AddToPath=0 /RegisterPython=0 /S
del %CONDA_INSTALLER%

:: Activate conda base environment (this is required for conda command to work)
call "%USERPROFILE%\miniconda3\Scripts\activate"

:: Create Conda environment
conda env create -f %YAML_FILE% --prefix .\%ENV_NAME%

:: Pause to see the output (optional, for debugging)
pause
