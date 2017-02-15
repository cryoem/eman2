#!/bin/bash

# Update PR refs for testing.
if [[ -n "${CIRCLE_PR_NUMBER}" ]]; then
    FETCH_REFS="${FETCH_REFS} +refs/pull/${CIRCLE_PR_NUMBER}/head:pr/${CIRCLE_PR_NUMBER}/head"
    FETCH_REFS="${FETCH_REFS} +refs/pull/${CIRCLE_PR_NUMBER}/merge:pr/${CIRCLE_PR_NUMBER}/merge"

# Retrieve the refs.
    git fetch -u origin ${FETCH_REFS}

# Checkout the PR merge ref.
    git checkout -qf "pr/${CIRCLE_PR_NUMBER}/merge"

# Check for merge conflicts.
    git branch --merged | grep "pr/${CIRCLE_PR_NUMBER}/head" > /dev/null
fi
