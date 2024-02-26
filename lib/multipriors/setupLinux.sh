#!/bin/bash

CONDA_INSTALLER=Miniconda3-latest-Linux-x86_64.sh
CONDA_INSTALL_PATH="$(pwd)/lib/multipriors/miniconda"
ENV_NAME="$(pwd)/lib/multipriors/multipriorsEnvLinux"
YAML_FILE="$(pwd)/lib/multipriors/multipriorsEnvLinux.yml"

# Download and install Miniconda
curl -O https://repo.anaconda.com/miniconda/${CONDA_INSTALLER}
bash ${CONDA_INSTALLER} -b -p ${CONDA_INSTALL_PATH}
rm ${CONDA_INSTALLER}

# Activate conda base environment (this is required for conda command to work)
source ${CONDA_INSTALL_PATH}/etc/profile.d/conda.sh
conda activate base

# Create Conda environment
conda env create -f ${YAML_FILE} --prefix ${ENV_NAME}

# Remove Miniconda installation directory
rm -rf ${CONDA_INSTALL_PATH}

# Display activation instructions for the new environment
echo "To activate the environment, run: conda activate ${CONDA_INSTALL_PATH}/${ENV_NAME}"