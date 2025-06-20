@echo off
set CONDA_INSTALLER=Miniconda3-latest-Windows-x86_64.exe
set ENV_NAME=%~dp0multiaxialEnv
set YAML_FILE=%~dp0multiaxialEnv.yml

:: Download and install Miniconda
curl -O https://repo.anaconda.com/miniconda/%CONDA_INSTALLER%
start /wait %CONDA_INSTALLER% /InstallationType=JustMe /AddToPath=0 /RegisterPython=0 /S
del %CONDA_INSTALLER%

:: Activate conda base environment (this is required for conda command to work)
call "%USERPROFILE%\miniconda3\Scripts\activate"

:: Create Conda environment
conda env create -f %YAML_FILE% --prefix %ENV_NAME%

:: Get the path of the conda environment
for /f "tokens=*" %%i in ('conda info --json ^| find "prefix"') do set "CONDA_ENV_PATH=%%i"
set CONDA_ENV_PATH=%CONDA_ENV_PATH:~14,-1%
echo Conda environment path: %CONDA_ENV_PATH%

:: Pause to see the output (optional, for debugging)
pause