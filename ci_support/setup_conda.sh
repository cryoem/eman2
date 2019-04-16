#!/usr/bin/env bash

set -xe

# Download and install Miniconda
export MINICONDA_URL="https://repo.continuum.io/miniconda"

curl -L -O "${MINICONDA_URL}/${MINICONDA_FILE}"
bash $MINICONDA_FILE -b

# Configure conda
source ${HOME}/miniconda2/bin/activate root
conda config --set show_channel_urls true

conda update conda -c defaults --yes --quiet
conda install "cmake=3.9" -c defaults --yes --quiet
