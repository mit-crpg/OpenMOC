import sys

## Run all tests.

## work on directory stuff

## This is not ready yet.

## Determine if slow version is wanted.
## TODO: command line args
long_version = False

## Load benchmark tests from inside sample-input

## C5G7
if long_version:
    sys.path.append('/Users/kenausis/Desktop/OpenMOC/sample-input/benchmarks/c5g7')
    import testc5g7

## Homog One Group
sys.path.append('/Users/kenausis/Desktop/OpenMOC/sample-input/benchmarks/homogeneous-one-group')
import testH1G

## Homog Two Groups
sys.path.append('/Users/kenausis/Desktop/OpenMOC/sample-input/benchmarks/homogeneous-two-group')
import testH2G

## LRA
sys.path.append('/Users/kenausis/Desktop/OpenMOC/sample-input/benchmarks/LRA')
#import testLRA

## Romano
sys.path.append('/Users/kenausis/Desktop/OpenMOC/sample-input/benchmarks/romano')
#import testromano


## Load regression tests

sys.path.append('/Users/kenausis/Desktop/OpenMOC/sample-input/')
import testregression


## Load tests from within openmoc
sys.path.append('/Users/kenausis/Desktop/OpenMOC/openmoc/')

import testlog
#import testmaterializeall
import testoptions
