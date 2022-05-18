#!/usr/bin/env bash

set -xeo pipefail

REPO_ROOT="$( cd "$( dirname "$0" )/.." >/dev/null && pwd )"

docker info

# In order for the conda-build process in the container to write to the mounted
# volumes, we need to run with the same id as the host machine, which is
# normally the owner of the mounted volumes, or at least has write permission
export HOST_USER_ID=$(id -u)
# Check if docker-machine is being used (normally on OSX) and get the uid from
# the VM
if hash docker-machine 2> /dev/null && docker-machine active > /dev/null; then
    export HOST_USER_ID=$(docker-machine ssh $(docker-machine active) id -u)
fi

DOCKER_IMAGE="quay.io/condaforge/linux-anvil-cos7-x86_64"

# Allow people to specify extra default arguments to `docker run` (e.g. `--rm`)
DOCKER_RUN_ARGS="${CONDA_FORGE_DOCKER_RUN_ARGS}"
if [ -z "${CI}" ]; then
    DOCKER_RUN_ARGS="-it ${DOCKER_RUN_ARGS}"
fi

docker pull "${DOCKER_IMAGE}"
docker run ${DOCKER_RUN_ARGS} \
           -v "${REPO_ROOT}":/home/conda/repo_root:rw,z,delegated \
           -e HOST_USER_ID \
           -e GIT_BRANCH \
           -e CI \
           -e CPU_COUNT \
           -e BUILD_OUTPUT_ID \
           "${DOCKER_IMAGE}" \
           bash \
           "/home/conda/repo_root/ci_support/build_steps.sh"
