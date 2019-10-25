#!/bin/bash
set -ex

# Install mpi4py and MPI-wrapped compilers for MPI configurations
if [[ $MPI == 'y' ]]; then
    $HOME/miniconda/bin/conda config --add channels conda-forge
    $HOME/miniconda/bin/conda install mpi4py openmpi-mpicc openmpi-mpicxx
    export LD_PRELOAD="$HOME/miniconda/envs/test-environment/lib/libmpi_cxx.so"
else
# Get regular C++11-compatible compiler for non-MPI configurations
    sudo apt-get install -qq gcc-6 g++-6
    sudo rm /usr/bin/gcc
    sudo rm /usr/bin/g++
    sudo ln -s /usr/bin/gcc-6 /usr/bin/gcc
    sudo ln -s /usr/bin/g++-6 /usr/bin/g++
fi
