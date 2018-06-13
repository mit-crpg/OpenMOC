#!/bin/sh
#######################################
# Specify nodes, processors per node
# and maximum running time
#######################################

#PBS -N xxxNAMExxx
#PBS -l nodes=xxxNODESxxx:ppn=xxxPPNxxx
#PBS -l walltime=xxxWTxxx
#PBS -l mem=xxxMEMxxxGB

######################################
# Enter directory and set PATH
######################################

cd $PBS_O_WORKDIR
PATH=$PBS_O_PATH

######################################
# Run OpenMOC - MAKE YOUR CHANGES HERE
######################################

module load gcc
module load mpich2/gnu
python xxxSTRIPxxx
mpiexec -machinefile xxxMACHINExxx ./xxxEXECUTExxx
