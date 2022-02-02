/**
 * \filename	main.c
 * \brief	This file contains the main functions required for executing 2d rectangular simulations
 * 		of conduction to steady state conditions. In particular the file is set up for a very 
 * 		specific boundary condition problem which is obtained numerically for comparison with 
 * 		an analytical solution. The user can specify the desired matrix size, and convergence criteria
 * 		however other parameters like actual x,y length dimensions or boundary temperature values must be
 * 		changed within this file.	
 * \author	F.OSuibhne
 * \date	15.08.21
 * \Usage	Command line input function "mpirun -np <n> ./<executable> <M> <N> <Tolerance> <Max_No._Iterations> <filename>"
 **/

#include<stdio.h>
#include<stdlib.h>  //< for test use of malloc
#include<mpi.h>
#include<string.h>
#include<math.h>
#include"jacobi_Itteration.h"
#include"memory_alloc.h"
#include"decomp2d.h"
#include"Parallel_write_to_file.h"

int main(int argc, char **argv){
	double x_dim = 1;		//< (m) length
	double y_dim = 1;		
	int M = 50; 			//< Grid rows
	int N = 50;  			//< Grid columns

	double ** U1_grid = NULL;	//< Grid Memory Pointer
	double ** U2_grid = NULL;
	double ** F_grid = NULL;	//< Pointer to spacial Internal Heat Generation grid values 	
	double Tolerance = 1E-3;	//< Default Tolerance for convergence criteria	
	double Num_Iteration = 20;	//< Desired Max No.Itterations

	int myid;			//< MPI processor ID
	int Total_nprocs;		//< Number of allocated MPI procs
	int s_x = 0, s_y = 0, e_x = 0, e_y = 0;
	char fname[200];
 
	int X_procs = 0;		//< Num Procs split over columns
	int Y_procs = 0;		//< Num Procs split over rows
	
	MPI_Comm cartcomm;

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &myid);
	MPI_Comm_size(MPI_COMM_WORLD,&Total_nprocs);
	
	double start_time = MPI_Wtime(); //< Obtain an initial run time

	//_________ Read In from CommandLine _______// 
  	if( argc != 6 ){		//< Read in on all processors
  		char errormsg[200];
		sprintf(errormsg,"Failure Insufficient command line arguments > mpirun -np <n> ./<executable> <M> <N> <Tollerance> <Max_no._Iteration> <filename>, occured at File: %s, Line: %d\n", __FILE__, __LINE__);
		perror(errormsg);	
		MPI_Abort(MPI_COMM_WORLD, 2);
   	}
 
	sscanf(argv[1], "%d", &M);	//< No Error checking exists on input values to ensure non negatives etc.	
	sscanf(argv[2], "%d", &N);	//  Obviously this is very bad practice and should be changed in next itteration //		
	sscanf(argv[3], "%lf", &Tolerance);	 
	sscanf(argv[4], "%lf", &Num_Iteration);	
	sprintf(fname, "%s", argv[5]);
 	printf("(rank: %d) M: %d, N: %d, Tol: %lf, Max_no._Itterations: %lf, file: %s \n",myid, M, N, Tolerance, Num_Iteration, fname);


	double dx = x_dim/((double)N-1);
	double dy = y_dim/((double)M-1);
	int grid_size = N;		//< Nessisary for Parallel print
		
	
	//_______ Decomposing Procs over X,Y axis _______//
	int Proc_Decomp_Error = MPI_Rect_Axis_Decomp(Total_nprocs, &X_procs, &Y_procs, N, M);	//< Splits processors over X,Y dimensions	
	if(Proc_Decomp_Error == 1){	//< Error check Proc Decomp
		char errormsg[150];
		sprintf(errormsg,"Failed to arrange processors over X,Y axis in File: %s, Line: %d\n", __FILE__, __LINE__);
		perror(errormsg);
		exit(EXIT_FAILURE);
	}

	int Decomp_check = MPI_Rec_Grid_Decomp2d(Total_nprocs, X_procs, Y_procs, N, M, myid, &s_x, &e_x, &s_y, &e_y);	//< Decomposes grid and allocates to processors
	if(Decomp_check == 1){		//< Error check Decomp2d
		char errormsg[150];
		sprintf(errormsg,"Failed to Decompose matrix in File: %s, Line: %d\n", __FILE__, __LINE__);
		perror(errormsg);
		exit(EXIT_FAILURE);
	}
	
	
	//________ CartComm _____//
	int ndims = 2;
	int dims[2] = {Y_procs,X_procs};	//< Note inverted order, to get Cart_create/Cart_shift
						//  to order as desired.
	int periods[2]={0,0};			//< Neither axis are periodic (no wrap around)
	int reorder = 0;
	int nbrleft, nbrright, nbrup, nbrdown;
	
	MPI_Cart_create(MPI_COMM_WORLD, ndims, dims, periods, reorder, &cartcomm);
	MPI_Cart_shift(cartcomm, 1, 1, &nbrleft, &nbrright);	//< Right shift info to neighbours	
	MPI_Cart_shift(cartcomm, 0, 1, &nbrup, &nbrdown);	//< Right shift info to neighbours	
	
	int Num_x = e_x-s_x+1;		//< quantity of x-axis nodes
  	int Num_y = e_y-s_y+1;		//< quantity of y-axis nodes
  	int MY_BUFSIZE = (Num_x+2)*(Num_y+2);	//< Required for Parallel Print to file	**changed for ghost (DELETE LATER): can this be isolated to function its used in ? 

	int sweep_sx = 1;
	int sweep_ex = Num_x;
	int sweep_sy = 1;
	int sweep_ey = Num_y;
	
	if(s_x == 1){
		nbrleft = MPI_PROC_NULL;	//< Specifying grid boundarys as non-neighbours	
		sweep_sx = 2;			//< Preventing boundary conditions perform sweep calc/overwrite
		//printf("myid=%d, sweep_sx=%d\n", myid, sweep_sx);
	}
	if(e_x == N){
		nbrright = MPI_PROC_NULL;
		sweep_ex = Num_x-1;
		//printf("myid=%d, sweep_ex=%d\n", myid, sweep_ex);
	}	
	if(s_y == 1){
		nbrup = MPI_PROC_NULL;		
		sweep_sy = 2;
		//printf("myid=%d, sweep_sy=%d\n", myid, sweep_sy);
	}
	if(e_y == M){
		nbrdown = MPI_PROC_NULL;
		sweep_ey = Num_y-1;
		//printf("myid=%d, sweep_ey=%d\n", myid, sweep_ey);

	}
	//_______ End CartComm _______//	


	// printf("Im proc %d of %d sharing s_x = %d e_x = %d, s_y = %d, e_y =%d, nbrleft = %d, nbrright= %d, nbrup= %d, nbrdown = %d \n", myid, Total_nprocs, s_x, e_x, s_y, e_y, nbrleft, nbrright, nbrup, nbrdown);	
	//MPI_Barrier(MPI_COMM_WORLD);	//< Usefull in combination with above print statement

	alloc_2D_rect_matrix(&U1_grid, Num_y+2, Num_x+2); 	//< Dynamic memory allocation with ghost row/columns		
	alloc_2D_rect_matrix(&U2_grid, Num_y+2, Num_x+2); 			
	alloc_2D_rect_matrix(&F_grid, Num_y+2, Num_x+2); 	//(DELETE LATER): Ghost columns are unessisary for F however its easier to keep them		
	
	//Label_2D_matrix_positions(U1_grid, Num_y+2, Num_x+2); 	//< Numbers matrix positions
	//Label_2D_matrix_positions(U2_grid, Num_y+2, Num_x+2); 	//< Numbers matrix positions
	//print_rect_matrix(U1_grid,Num_y, Num_y);	  		//< prints matrix 
		
	Initialize_2D_Grid_Values(U1_grid, Num_y+2, Num_x+2, 0.0); 	//< sets grid values to specified value 
	Initialize_2D_Grid_Values(U2_grid, Num_y+2, Num_x+2, 0.0);
	Initialize_2D_Grid_Values(F_grid, Num_y+2, Num_x+2, 0.0);
	
	double Tmax=100;	//< Specific to the Analytical problem we are trying to simulate numerically
	double T1=0;
	
	Initialize_2D_Grid_Boundries(U1_grid, T1, T1, T1, T1, s_x, e_x, s_y, e_y, M, N, Num_y, Num_x); //< Setting Grid Boundary conditions
	Initialize_2D_Grid_Boundries(U2_grid, T1, T1, T1, T1, s_x, e_x, s_y, e_y, M, N, Num_y, Num_x);	
		
	// Initializing non-uniform Bottom Boundary Condition 
	if(e_y == M){	//< Specific to to Analytical Problem we are trying to simulate numerically
		for(int i=1; i< Num_x+1; i++){
			U1_grid[Num_y][i] = T1 + Tmax*sin((M_PI*((double)s_x+i-2)*dx)/x_dim); 
			U2_grid[Num_y][i] = T1 + Tmax*sin((M_PI*((double)s_x+i-2)*dx)/x_dim);  
		}
	}

	//______ Data Exchange/Jacobi ______//
	
	const int Max_itt = (int)(Num_Iteration/2);	//< Specify some global Max_itt
	double Local_grid_diff = 0.0;
	double Global_grid_diff = 0.0;
	double Convergence = 0.0;

	//___ Datatype for exchanging columns __//
	MPI_Datatype column_type; 	
	MPI_Type_vector(Num_y, 1, Num_x+2, MPI_DOUBLE, &column_type);	
	MPI_Type_commit(&column_type);	//< Required type for non-contiguous access of columns
	
	int Itteration=0;
	for(Itteration = 0; Itteration<Max_itt; Itteration++ ){
	
		Exchange_Data_2D(Num_x, Num_y, U1_grid, MPI_COMM_WORLD, s_x, e_x, s_y, e_y, nbrleft, nbrright, nbrup, nbrdown, column_type);
	
		Sweep_Solve_2D(U1_grid, U2_grid, F_grid, sweep_sx, sweep_ex, sweep_sy, sweep_ey, dx, dy); //< Sweep function assumes convergance test	
		
		Exchange_Data_2D(Num_x, Num_y, U2_grid, MPI_COMM_WORLD, s_x, e_x, s_y, e_y, nbrleft, nbrright, nbrup, nbrdown, column_type);
	
		Sweep_Solve_2D(U2_grid, U1_grid, F_grid, sweep_sx, sweep_ex, sweep_sy, sweep_ey, dx, dy); //< Sweep function assumes convergance test

		//____ Convergence Test ____//
		Local_grid_diff = Calc_grid_diff_pointwise(U1_grid, U2_grid, Num_x, Num_y, Tolerance);		//< Obtain local sum of Indicators
  		MPI_Allreduce(&Local_grid_diff, &Global_grid_diff, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);	//< Obtain global sum of Indicator
		Convergence = Global_grid_diff; //< Value Indicator for Entire Grid 

	
		if(Convergence == 0){
			
			if(myid == 0){
				printf("Simulation has converged to the specified tollerance %lf after %d sweeps of the grid, \n", Tolerance, 2*Itteration);
			}

			break;	//< Breaks Loop if Convergence criteria met //		
		}
	
	}
	
	//______ All Processors Print Own Section to File ______//	
	Parallel_Print_Double_To_File_wBoundary( fname, Num_x, Num_y, MY_BUFSIZE, *(&U1_grid[1])+1, grid_size, myid, s_x, e_x, s_y, e_y);
		
	//_____ Free dynamically allocated memory _____// 		
	MPI_Type_free(&column_type);	//< Datatype for exchange
	free(*(&U1_grid[0]));
	free(*(&U1_grid));	
 	free(*(&U2_grid[0]));
	free(*(&U2_grid));	
 	free(*(&F_grid[0]));
	free(*(&F_grid));	

  	double end_time = MPI_Wtime();
	if(myid == 0){
		printf("Convergence = %lf to a desired Tollerance = %lf, Itteration %d \n", Convergence, Tolerance, Itteration*2);	
		printf("Time elapsed = %lfs \n", end_time-start_time);
	}

	MPI_Finalize();
	return(0);
}
