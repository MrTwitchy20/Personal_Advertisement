#!/bin/bash
##SBATCH -n 9			# 12 CPU cores
##SBATCH -t 0:10:00		# time
##SBATCH -p debug		# partition name (debug, compute)
##SBATCH -J FiniteDifference	# sensible job name

# load the correct modules
module load gcc openmpi Python/3.7.6-gnu intel/18.0.5

# launch the code
make
mpirun -n 9 ./2D_Rect_Fixed_Boundary_SS 50 50 .001 50000 outputfile 
