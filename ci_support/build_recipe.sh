#!/usr/bin/env bash

set -x

# Download and install Miniconda
export MINICONDA_URL="https://repo.continuum.io/miniconda"

curl -L -O "${MINICONDA_URL}/${MINICONDA_FILE}"
bash $MINICONDA_FILE -b

# Configure conda
source ${HOME}/miniconda2/bin/activate root
conda config --set show_channel_urls true

export CPU_COUNT=2

conda install conda-build -c defaults --yes --quiet

conda build recipes/eman -c cryoem -c defaults -c conda-forge --numpy 1.8
