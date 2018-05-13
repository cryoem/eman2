#!/usr/bin/env bash

set -xe

GITHASH=$(python -c "from EMAN2_meta import GITHASH; print(GITHASH)")

if [ ! -z "$JENKINS_HOME" ];then
    THIS_COMMIT_HASH="$GIT_COMMIT_SHORT"
elif [ ${CIRCLECI} ];then
    THIS_COMMIT_HASH="$CIRCLE_SHA1"
elif [ ${TRAVIS} ];then
    THIS_COMMIT_HASH="$TRAVIS_COMMIT"
else
    THIS_COMMIT_HASH="       "
fi

GITHASH=${GITHASH:0:7}
THIS_COMMIT_HASH=${THIS_COMMIT_HASH:0:7}

test "${GITHASH}" == "${THIS_COMMIT_HASH}"
