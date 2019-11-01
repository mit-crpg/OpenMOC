#!/bin/bash
set -ex

# Handle Python third-party packages
if [ -d "$HOME/miniconda" ] && [ -d "$HOME/miniconda/bin" ]; then
    echo "Using cached Miniconda"
    export PATH="$HOME/miniconda/bin:$PATH"
else
    wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh
    rm -rf $HOME/miniconda
    bash miniconda.sh -b -p $HOME/miniconda
    export PATH="$HOME/miniconda/bin:$PATH"
    hash -r
    conda config --set always_yes yes --set changeps1 no
    conda config --add channels conda-forge
  # conda update -q --force conda
    conda info -a
    conda create -q -n test-environment python=$TRAVIS_PYTHON_VERSION numpy scipy h5py pandas matplotlib pillow uncertainties pandas
fi
source activate test-environment
sudo apt-get install swig

# For coverage testing of C++ and Python source files
pip install coverage
pip install cpp-coveralls
pip install coveralls

# Install openmc for Superhomogeneization test
if [ -d "$HOME/openmc" ]  && [ -d "$HOME/openmc/openmc" ]; then
    echo "Using cached OpenMC"
else
    git clone --depth=50 --branch=develop https://github.com/openmc-dev/openmc.git $HOME/openmc
fi
export PYTHONPATH="$HOME/openmc/:$PYTHONPATH"

# Install mpi4py and MPI-wrapped compilers for MPI configurations
if [[ $MPI == 'y' ]]; then
    conda config --add channels conda-forge
    conda install -y mpi4py openmpi-mpicc openmpi-mpicxx
    export LD_PRELOAD="$HOME/miniconda/envs/test-environment/lib/libmpi_cxx.so"
# Get regular C++11-compatible compiler for non-MPI configurations
else
    sudo apt-get install -qq gcc-6 g++-6
    sudo rm /usr/bin/gcc
    sudo rm /usr/bin/g++
    sudo ln -s /usr/bin/gcc-6 /usr/bin/gcc
    sudo ln -s /usr/bin/g++-6 /usr/bin/g++
fi
