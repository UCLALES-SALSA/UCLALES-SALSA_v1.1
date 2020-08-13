#!/bin/bash -l
#SBATCH --job-name=DYCOMS1ppb
#SBATCH --ntasks=144
#SBATCH --time=48:00:00
#SBATCH --output=output.txt
#SBATCH --error=error.txt
#SBATCH --account="project_2002155"
#SBATCH --partition="fmi"
module load StdEnv
module load netcdf
module load parallel-netcdf
module load netcdf-fortran
export MPICH_ENV_DISPLAY=1
cd /fmi/scratch/project_2002155/kudzotsa/DYCOMS/UCLALES-SALSA_gases_3ppb/bin
srun -n 144 ./les.mpi | tee uclales-salsa.output

set -ex

