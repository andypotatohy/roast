#!/usr/bin/env bash

CONDA_INSTALLER=Miniconda3-latest-MacOSX-x86_64.sh
ENV_NAME=multipriorsEnvMac
YAML_FILE=multipriorsEnvMac.yml

# Download and install Miniconda
curl -O https://repo.anaconda.com/miniconda/${CONDA_INSTALLER}
bash ${CONDA_INSTALLER} -b -p $HOME/miniconda3
rm ${CONDA_INSTALLER}

# Activate conda base environment (this is required for conda command to work)
source $HOME/miniconda3/etc/profile.d/conda.sh
conda activate base

# Create Conda environment
conda env create -f ${YAML_FILE} --prefix ./${ENV_NAME}

# Display activation instructions for the new environment
echo "To activate the environment, run: conda activate ./${ENV_NAME}"

# Optional: Pause to see the output (for debugging)
read -p "Press enter to exit..."
