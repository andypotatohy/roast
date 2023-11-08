@echo off
:: Command to create a Conda environment
conda env create -f MultiPriors_env.yml --prefix .\MultiPriors_env_test

:: Pause to see the output (optional, for debugging)
pause
