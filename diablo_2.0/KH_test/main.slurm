#!/bin/bash
# Job name:
#SBATCH --job-name=KH_test
#
# Account:
#SBATCH --account=co_stratflows
#
# Partition:
#SBATCH --partition=savio4_htc
#
# Number of tasks needed:
#SBATCH --ntasks=4
#
# Processors per task:
#SBATCH --cpus-per-task=1
#
# Wall clock limit:
#SBATCH --time=08:00:00
#
## Command(s) to run (example):
module load gcc/11.4.0 openmpi/4.1.6 hdf5/1.14.3 fftw/2.1.5
mpirun -n 4 ./diablo >output.dat


