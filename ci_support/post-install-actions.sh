#!/usr/bin/env bash

set -xe

source ${PREFIX}/bin/activate

mkdir ${PREFIX}/install_logs

conda info -a          > ${PREFIX}/install_logs/info_log.txt 2>&1
conda list             > ${PREFIX}/install_logs/list_log.txt 2>&1
conda list --explicit >> ${PREFIX}/install_logs/list_log.txt 2>&1

SP_DIR=$(python -c "import site; print(site.getsitepackages()[0])")

if [ -d site-packages ]; then
    cp -av site-packages/* "${SP_DIR}"
    rm -rv site-packages
fi > ${PREFIX}/install_logs/install_log.txt 2>&1

conda install eman-deps=14.1 -c cryoem -c defaults -c conda-forge -y | tee -a ${PREFIX}/install_logs/install_log.txt 2>&1

cat <<EOF

INSTALLATION IS NOW COMPLETE

Please, go to http://blake.bcm.edu/emanwiki/EMAN2/Install/BinaryInstallAnaconda
for detailed installation instructions, testing and troubleshooting information.
If this installation is on a Linux cluster,
you will require additional steps before installation is complete!

EOF
