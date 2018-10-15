#!/bin/sh

set -ev

# Run all debug tests
if [ "$TRAVIS_PULL_REQUEST" != "false" ]; then
  ./run_tests.py -c
fi
