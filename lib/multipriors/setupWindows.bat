@echo off
:: Command to create a Conda environment

conda env create -f  multipriorsEnv.yml --prefix ./multipriors_env
:: Pause to see the output (optional, for debugging)
pause
