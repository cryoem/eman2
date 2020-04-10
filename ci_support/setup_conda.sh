#!/usr/bin/env bash

set -xe

# Download and install Miniconda
export MINICONDA_URL="https://repo.continuum.io/miniconda"

curl -L -O "${MINICONDA_URL}/${MINICONDA_FILE}"
bash $MINICONDA_FILE -b

# Configure conda
source ${HOME}/miniconda3/bin/activate root
conda config --set show_channel_urls true
conda config --set auto_update_conda False

conda install conda conda-build cmake=3.14 -c defaults --yes
