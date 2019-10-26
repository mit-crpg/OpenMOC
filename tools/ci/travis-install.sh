#!/bin/bash
set -ex

# Handle Python third-party packages
if [ -f $HOME/miniconda ]; then
    echo "Using cached Miniconda"
else
    wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh
    bash miniconda.sh -b -p $HOME/miniconda
    export PATH="$HOME/miniconda/bin:$PATH"
    hash -r
    conda config --set always_yes yes --set changeps1 no
  # conda update -q --force conda
    conda info -a
    conda create -q -n test-environment python=$TRAVIS_PYTHON_VERSION numpy scipy h5py pandas matplotlib pillow
    source activate test-environment
    sudo apt-get install swig
  # For coverage testing of C++ and Python source files
    pip install coverage
    pip install cpp-coveralls
    pip install coveralls
fi

# Install openmc for Superhomogeneization test
if [-f "openmc" ]; then
    echo "Using cached OpenMC"
else
    git clone --depth=50 --branch=master https://github.com/openmc-dev/openmc.git openmc
fi 

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
