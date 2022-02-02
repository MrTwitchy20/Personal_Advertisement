/**
 * \filename	main.c
 * \brief	This file contains the main function for executing 2d rectangular simulations	
 * \author	F.OSuibhne
 * \date	05.09.21
 * \notes	- Command line input function
 *
 *
 **/

#include<stdio.h>
#include<stdlib.h>  //< for test use of malloc
#include<mpi.h>
#include<string.h>
#include<math.h>

void Parallel_Print_Double_To_File_wBoundary( const char *, const int , const int, const int, double *, const int, const int, const int, const int, const int, const int);
void Exchange_Data_2D(int, int,  double * const * const, MPI_Comm, int, int, int, int,int, int, int, int, MPI_Datatype);
//void Sweep_Solve_1D_SteadyState_Annulus_FixedBound(double *const *const, double *const *const, double *const *const, int, int, int, int, int, int, int, int, double, double);
void Sweep_Solve_1D_Transient_Annulus_Heat_Gen(double *const *const, double *const *const, double *const *const, const double, const double, const double, const double, const double, const double, const double, const double, const double, const int);
void Transient_1D_Stability_Check( const double, const double, const double);	

double Calc_grid_diff(double *const *const, double *const *const, int , int );
	
int alloc_2D_rect_matrix(double ***, const int, const int);
void print_rect_matrix(double *const *const, const int, const int);
void Initialize_2D_Grid_Values(double * const *const, const int, const int, const double);
void Label_2D_matrix_positions(double *const *const, const int, const int);

int MPI_Rect_Axis_Decomp(const int, int *, int *, const int, const int);
int MPI_Rec_Grid_Decomp2d(const int, const int, const int, const int, const int, const int, int *,int *,int *,int *);

void Initialize_2D_Grid_Boundries(double * const * const, const double, const double, const double, const double, const int, const int, const int, const int, const int, const int, const int, const int);	

int main(int argc, char **argv){
	
	//___ Initialise Grid Properties ___//
	int M = 10; // Grid rows 
	int N = 15;  // Grid columns 
	// Pure Copper: Heat transfer textbook third edition john H. Lienhard IV
	const double k = 437.8;
	const double alpha= .0001157; // m^2/s
	const double r_i = .002; //< m, -> 2mm	cross section of an annulus, 18650 cell (image from journal of power sources 472 2020)
	const double r_o = .009; //< m, -> 9mm
	const double delta_r = (r_o-r_i)/((double)N-1); 	//< step size in r Ozisik p465
	const double delta_t = .001; // seconds  	
	
	Transient_1D_Stability_Check(alpha, delta_t, delta_r);	

	double ** U1_grid = NULL;	// Grid Memory Pointer
	double ** U2_grid = NULL;
	double ** F_grid = NULL;	// Pointer to function value grid 	
	double Tolerance = 1E-6;	//< Set Tolerance for convergence criteria	

	int myid;		// MPI processor ID
	int Total_nprocs;	// Number of allocated MPI procs
	int s_x = 0, s_y = 0, e_x = 0, e_y = 0;
	int grid_size = N;	// Nessisary for Parallel print
	
	char fname[200];
  	int fname_len;
 
	int X_procs = 0;	// Num Procs split over columns
	int Y_procs = 0;	// Num Procs split over rows
	
	MPI_Comm cartcomm;

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &myid);
	MPI_Comm_size(MPI_COMM_WORLD,&Total_nprocs);


	//_________ Read In from CommandLine _______//
	if( argc != 2 ){	// Read in on all processors
  		fprintf(stderr,"Error: usage: %s <output filename>\n",argv[0]);
   		MPI_Abort(MPI_COMM_WORLD, 2);
   	}
	sprintf(fname, "%s", argv[1]);		
	
	//_______ Decomposing Procs over X,Y axis _______//
	int Proc_Decomp_Error = MPI_Rect_Axis_Decomp(Total_nprocs, &X_procs, &Y_procs, N, M);	//< Splits processors over X,Y dimensions	
	if(Proc_Decomp_Error == 1){	//Error check Proc Decomp
		char errormsg[150];
		sprintf(errormsg,"Failed to arrange processors over X,Y axis in File: %s, Line: %d\n", __FILE__, __LINE__);
		perror(errormsg);
		exit(EXIT_FAILURE);
	}

	int Decomp_check = MPI_Rec_Grid_Decomp2d(Total_nprocs, X_procs, Y_procs, N, M, myid, &s_x, &e_x, &s_y, &e_y);	//< Decomposes grid and allocates to processors
	if(Decomp_check == 1){	//Error check Decomp2d
		char errormsg[150];
		sprintf(errormsg,"Failed to Decompose matrix in File: %s, Line: %d\n", __FILE__, __LINE__);
		perror(errormsg);
		exit(EXIT_FAILURE);
	}
	
	
	//________ CartComm _____//
	int ndims = 2;
	int dims[2] = {Y_procs,X_procs};	//< Note inverted order, to get Cart_create/Cart_shift
						//  to order as desired.
	int periods[2]={1,0};			//< Periodic in Y dimension
	int reorder = 0;
	int nbrleft, nbrright, nbrup, nbrdown;
	
	MPI_Cart_create(MPI_COMM_WORLD, ndims, dims, periods, reorder, &cartcomm);
	MPI_Cart_shift(cartcomm, 1, 1, &nbrleft, &nbrright);	//< Right shift info to neighbours	
	MPI_Cart_shift(cartcomm, 0, 1, &nbrup, &nbrdown);	//< Right shift info to neighbours	
	
	int Num_x = e_x-s_x+1;	//< quantity of x-axis nodes
  	int Num_y = e_y-s_y+1;	//< quantity of y-axis nodes
  	int MY_BUFSIZE = (Num_x+2)*(Num_y+2);	//< Required for Parallel Print to file	**changed for ghost (DELETE LATER): can this be isolated to function its used in ? 

	int sweep_sx = 1;
	int sweep_ex = Num_x;
	int sweep_sy = 1;
	int sweep_ey = Num_y;
	
	if(s_x == 1){
		nbrleft = MPI_PROC_NULL;	//< Specifying grid boundarys as non-neighbours	
		sweep_sx = 2;			//< Preventing boundary conditions perform sweep calc/overwrite
		printf("myid=%d, sweep_sx=%d\n", myid, sweep_sx);
	}
	if(e_x == N){
		nbrright = MPI_PROC_NULL;
		sweep_ex = Num_x-1;
		printf("myid=%d, sweep_ex=%d\n", myid, sweep_ex);
	}	
	/*
	if(s_y == 1){
		nbrup = MPI_PROC_NULL;	//< Specifying grid boundarys as non-neighbours	
		sweep_sy = 2;			//< Preventing boundary conditions perform sweep calc/overwrite
		printf("myid=%d, sweep_sy=%d\n", myid, sweep_sy);
	}
	if(e_y == M){
		nbrdown = MPI_PROC_NULL;
		sweep_ey = Num_y-1;
		printf("myid=%d, sweep_ey=%d\n", myid, sweep_ey);
	}	*/
	
	
  	alloc_2D_rect_matrix(&U1_grid, Num_y+2, Num_x+2); 	//< Dynamic memory allocation with ghost row/columns		
	alloc_2D_rect_matrix(&U2_grid, Num_y+2, Num_x+2); 			
	alloc_2D_rect_matrix(&F_grid, Num_y+2, Num_x+2); 	
	
	Initialize_2D_Grid_Values(U1_grid, Num_y+2, Num_x+2, 0.0); 	//< sets grid values to "myid"
	Initialize_2D_Grid_Values(U2_grid, Num_y+2, Num_x+2, 0.0);
	Initialize_2D_Grid_Values(F_grid, Num_y+2, Num_x+2, 0.0);
	Initialize_2D_Grid_Boundries(U1_grid, 1.0, 2.0, 0.0, 0.0, s_x, e_x, s_y, e_y, M, N, Num_y, Num_x); //< Setting Grid Boundary conditions
	Initialize_2D_Grid_Boundries(U2_grid, 1.0, 2.0, 0.0, 0.0, s_x, e_x, s_y, e_y, M, N, Num_y, Num_x); //< No Y boundarys
	
	double start_time = MPI_Wtime();

	//______ Data Exchange/Jacobi ______//
	
	const int Max_itt = 2000000; //< seconds, Specify value for termination
	double Local_grid_diff = 0.0;
	double Global_grid_diff = 0.0;
	double RMS_diff = 0.0;

	//___ Datatype for exchanging columns __//
	MPI_Datatype column_type; 	
	MPI_Type_vector(Num_y, 1, Num_x+2, MPI_DOUBLE, &column_type);	
	MPI_Type_commit(&column_type);	//< Required type for non-contiguous access of columns

	for(int Itteration = 0; Itteration<Max_itt; Itteration++ ){	//< Secondary criteria for termination
		Exchange_Data_2D(Num_x, Num_y, U1_grid, MPI_COMM_WORLD, s_x, e_x, s_y, e_y, nbrleft, nbrright, nbrup, nbrdown, column_type);
		Sweep_Solve_1D_Transient_Annulus_Heat_Gen(U1_grid, U2_grid, F_grid, alpha, sweep_sx, sweep_ex, sweep_sy, sweep_ey, delta_r, delta_t, k, r_i,s_x);

		Exchange_Data_2D(Num_x, Num_y, U2_grid, MPI_COMM_WORLD, s_x, e_x, s_y, e_y, nbrleft, nbrright, nbrup, nbrdown, column_type);
		Sweep_Solve_1D_Transient_Annulus_Heat_Gen(U2_grid, U1_grid, F_grid, alpha, sweep_sx, sweep_ex, sweep_sy, sweep_ey, delta_r, delta_t, k, r_i, s_x);
	}

	double end_time = MPI_Wtime();
	if(myid == 0){
		printf("Time elapsed = %lfs \n", end_time-start_time);
	}

	Parallel_Print_Double_To_File_wBoundary( fname, Num_x, Num_y, MY_BUFSIZE, *(&U1_grid[1])+1, grid_size, myid, s_x, e_x, s_y, e_y);
	
	//_____ Free dynamically allocated memory _____// 		
	MPI_Type_free(&column_type);	
	free(*(&U1_grid[0]));
	free(*(&U1_grid));	
 	free(*(&U2_grid[0]));
	free(*(&U2_grid));	
 	free(*(&F_grid[0]));
	free(*(&F_grid));	
	
  	MPI_Finalize();
	return(0);
}
