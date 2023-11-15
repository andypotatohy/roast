@echo off
:: Command to create a Conda environment

conda env create -f  MultiPriors_env.yml --prefix ./roast-3.0/lib/MultiPriors_WEB/MultiPriors_env
:: Pause to see the output (optional, for debugging)
pause
