#!/usr/bin/env bash

set -xe

source ${PREFIX}/bin/activate

conda config --env --set auto_update_conda False

mkdir ${PREFIX}/install_logs

conda info -a            | tee ${PREFIX}/install_logs/info_log.txt 2>&1
conda list               | tee ${PREFIX}/install_logs/list_log.txt 2>&1
conda list --explicit | tee -a ${PREFIX}/install_logs/list_log.txt 2>&1

case ${EMAN_INSTALL_DONT_UPDATE_DEPS:-} in
    0|"")
        conda install eman-deps=27.0 -c cryoem -c defaults -c conda-forge -y | tee -a ${PREFIX}/install_logs/install_log.txt 2>&1
        ;;
    *)
        echo "WARNING: Skipping installation of dependencies per user request..."
        ;;
esac

cat <<EOF

INSTALLATION IS NOW COMPLETE

Please, go to http://blake.bcm.edu/emanwiki/EMAN2/Install/BinaryInstallAnaconda
for detailed installation instructions, testing and troubleshooting information.
If this installation is on a Linux cluster,
you will require additional steps before installation is complete!

EOF
