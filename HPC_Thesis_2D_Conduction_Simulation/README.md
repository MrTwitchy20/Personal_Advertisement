# Thesis_Repo
author: 	F.OSuibhne
date:		23.09.21

This repository will serve as a cloud based backup of all code based works to date related to my thesis. It will provide a means of tracking changes,
maintaining file continuity between the multiple devices I use and also provide my supervisor access to track progress and comment/review content.

Running Instructions:
	1. Load intel compiller "module load intel"
	
	2. Load openMPI on the cluster using command "module load openmpi"
	
	3. Go to directory with Makefile and use the command "make" to compile relevent code
	
	4. Items can then be run individually using command:
		"mpirun -np <number of procs> ./<executable> <M> <N> <Tolerance> <Max_no.Iteration> <filename>"
	   executable: Two executables are available 2D_Annulus_Fixed_Boundary or 2D_Rect_Fixed_Boundary_SS
	   M: Number of Rows desired
	   N: Number of Columns desired
	   Tolerance: The convergence criteria, maximum allowed difference between any two nodes in consecutive itterations of a numerical solution.
	   Max_no.Iteration: This is the maximum number of grid solves to be attempted before the program must exit (secondary exit criteria) 
	   filename: The name of the desired output binary file name which the programs will print grids to
	
	5. Items can be scheduled for running on a cluster with slurm using the Sbatch.sh file.
	   This file must be edited each time to specify the desired execution information as well as SBATCH commands such as number of processors,
	   allocated execution time, cluster partition and desired job name. Once this is done the command "sbatch Sbatch.sh" is used to run the job.
	
	6. Load python using the command "module load "Python"

	7. Upon completion the executed code should have produced a file of the specified name. This is not currently user readable. To convert the file from 
	   binary to a .txt with decimal representation the command "python Binary_to_decimal_txt.py <M> <N> <filename>"
	   M: Number of Rows used originally in to run the simulation
	   N: Number of Columns used originally to run the simulation
	   filename: Name of the output binary file from the original simulation
	   This will generate a new file with the original "filename" and an appended "_decimal.txt" which can now be opened and read to show the grid structure
	
	8. The generated output .txt grid files can now be used to generate visual graphic representations and compared against similar simulations and analytical solutions
	   opening Mathlab and importing the .txt files into the three available Mathlab programs "Analytical_poisson", "Cylindrical_plot" and "Test_multiple_procs_match". 
	   Both the Analytical_poisson & Cylindrical_plot are hardcoded to compare a specific grid size/ types and analytical solutions and require input files matching these
	   variables.  
Notes:	 
- The binary to decimal converter can be configured to print the desired level of precision and adjusted for Int types instead of double
- Two aditional executables exist 2D_Annulus_Transient_Fixed_Boundary & 2D_Annulus_Convective however neither are in complete and working state suitable for running.
- A number of parameters such as grid boundary conditions and grid dimensions are hardcoded into the C files and must be changed and recompiled for simulating different conditions.
- A file named Dear_Diary_Thesis.txt exists which was used to keep track of progress, problems encountered and general project related comments over the duration of the Thesis.   
