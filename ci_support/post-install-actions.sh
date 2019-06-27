#!/usr/bin/env bash

set -xe

source ${PREFIX}/bin/activate

mkdir ${PREFIX}/install_logs

conda info -a          > ${PREFIX}/install_logs/info_log.txt 2>&1
conda list             > ${PREFIX}/install_logs/list_log.txt 2>&1
conda list --explicit >> ${PREFIX}/install_logs/list_log.txt 2>&1

conda install --force-reinstall conda=4.6.14 conda-build=3.17.8 pytz backports backports.functools_lru_cache filelock tqdm -y
conda install eman-deps=14.1 -c cryoem -c defaults -c conda-forge -y | tee ${PREFIX}/install_logs/install_log.txt 2>&1

cat <<EOF

INSTALLATION IS NOW COMPLETE

Please, go to http://blake.bcm.edu/emanwiki/EMAN2/Install/BinaryInstallAnaconda
for detailed installation instructions, testing and troubleshooting information.
If this installation is on a Linux cluster,
you will require additional steps before installation is complete!

EOF
