/**
  * \filename	Jacobi_Itteration.c
  * \brief	This file contains the required functions to perform Jacobi Itteration
  * \notes	
  * \todo	Change file name	
  * \author	F.OSuibhne
  * \date	19-08-21
  **/

#include<mpi.h>
#include<stdio.h>
#include<string.h>	//< required for error messages
#include<stdlib.h>

/*
 * \brief	function is responcible for 2D data exchange between processors,
 * 		Left, Right, Up, Down neighbours all exchange ghost columns/rows
 *		with one another in a non-blocking fashion. These are all placed under the same tag.
 * \todo	Should include error checking if possible
 * \note	- Non-ghost row/column of sender to neighbours ghost row/column
 * 		- Message tags all identical.
 * 
 * \visual	Initial Primary Matrix with Boundary values L,R,T,B: 
 *     		|C|_T___T___T___T_|C|	
 *              |L| 0   1   2   3 |R| 
 *              |L| 4   5   6   7 |R| 
 *              |L| 8   9  10  11 |R| 
 *              |L|12__13__14__15_|R| 
 *		|C| B   B   B   B |C|
 *	  	
 *	  	Matrix section on each processor before exchange:
 *	  	  Proc 0    Proc 1     Proc 2     Proc 3
 * 		|C T T 0| |0 T T C| |0  0  0 0| |0  0  0 0|
 * 		|L 0 1 0| |0 2 3 R| |L  8  9 0| |0 10 11 R|
 *		|L 4 5 0| |0 6 7 R| |L 12 13 0| |0 14 15 R|
 *              |0 0 0 0| |0 0 0 0| |C  B  B 0| |0  B  B C|
 *		
 *		Matrix section on each processor after exchange:
 *		  Proc 0    Proc 1     Proc 2     Proc 3
 * 		|C T T 0| |0  T  T C| |0  4  5  0| | 0  6  7 0|
 * 		|L 0 1 2| |1  2  3 R| |L  8  9 10| | 9 10 11 R|
 * 		|L 4 5 6| |5  6  7 R| |L 12 13 14| |13 14 15 R|
 *              |0 8 9 0| |0 10 11 0| |C  B  B  0| | 0  B  B C|
 */ 
void Exchange_Data_2D(int nx, int ny, double *const *const A, MPI_Comm comm, int s_x, int e_x, int s_y, int e_y, int nbrleft, int nbrright, int nbrup, int nbrdown, MPI_Datatype column_type){

	MPI_Request requests[4];	//< declare request for each direction

	MPI_Isend(&A[1][1], 1, column_type, nbrleft, 0, comm, &requests[2]);		//< Send left column to left neighbour
	MPI_Isend(&A[1][nx], 1, column_type, nbrright, 0, comm, &requests[3]);		//< Send right column to right neighbour

	MPI_Irecv(&A[1][nx+1], 1, column_type, nbrright, 0, comm, &requests[2]);	//< Receive left column from right neighbour into my right ghost
	MPI_Irecv(&A[1][0], 1, column_type, nbrleft, 0, comm, &requests[3]);		//< Receive right column from left neighbour into my left ghost	

	MPI_Isend(&A[1][1], nx, MPI_DOUBLE, nbrup, 0, comm, &requests[0]);		//< Send top row up to neibhour above
	MPI_Isend(&A[ny][1], nx, MPI_DOUBLE, nbrdown, 0, comm, &requests[1]);		//< Send bottom row down to neighbour beneath
	
	MPI_Irecv(&A[ny+1][1], nx, MPI_DOUBLE, nbrdown, 0, comm, &requests[0]);		//< Receive top row from neighbour below into my bottom ghost	
	MPI_Irecv(&A[0][1], nx, MPI_DOUBLE, nbrup, 0, comm, &requests[1]);		//< Receive top bottom row from neighbour above into my top ghost

	MPI_Waitall(4, requests, MPI_STATUS_IGNORE);	 	
}



/*
 * \brief	Function is responsible for sweeping across whole array performing a single calculation 
 * 		itteration of poisson.  		
 * \note	- Calculation occurs locally on each proc
 * 		- Three grids are passed in, U_1 & U_2 hold previous position values and new position values respectivly
 * 		  as the new values are soley based on the previous itteration seperate grids are required to prevent mixing
 * 	 	  of new and old itteration data. F grid can be used for spacialy varying heat generation but is currently 
 * 	 	  unused in this version. It requires devision by a k thermal conductivity term not passed in however would be
 * 		  better to apply this division once externally to the itteration loop as it is constant. 
 *		- Should be accessing data with unit stride 
 *		- Accessing Data U[y][x] hence U[i][j]
 * \visual	
 * 			     U_1[i][j+1]
 *				  |
 *				  |
 *		U_1[i-1][j]----U_1[i][j]----U_1[i+1][j]
 *				  |
 *				  |
 *			     U_1[i][j-1]
 *		
 */
void Sweep_Solve_2D(double *const *const U_1, double *const *const U_2, double *const *const F, int sweep_sx, int sweep_ex, int sweep_sy, int sweep_ey, int dx, int dy){	
				
  	for(int i=sweep_sy; i<=sweep_ey; i++){	//< Sweep ignores ghost columns //
    		for(int j= sweep_sx; j<= sweep_ex; j++){ 
      			U_2[i][j] = 0.25 * (U_1[i-1][j] + U_1[i+1][j] + U_1[i][j+1] + U_1[i][j-1] - dx*dy*F[i][j]); 
			// Adapted Finite difference formula  
    		}
  	}
}


/*
 * \brief	This function performs a grid sweep solve across an Annulus shaped cross section with 2D conductive heat transfer
 * 		in both the radial r direction (cartesian x equivalent) and tangentially uppercase phi (cartesian y equivalent)
 * \note	Not subject to stability criterion which comes from time dependancy.
 * \todo	-delta_r, delta_phi should also be constants passed in indicating step size presumable relative to an input r,phi dimension as opposed to M,N as we have
 * 		- If input Matrix size is M*N, there are also input radius required, it needs to be determined whether boundary conditions will utilize ghosts or simply
 * 		require updating to shift in left and right to use actual grid points as opposed to ficticious as is currently in place.
 * 		If deciding to shift in consider what happens as main for loop must be changed to prevent same grid position attempting to apply finite diffrence.
 * \notes	- Accessing data with unit stride
 * 		- Formula: Page 465 of Osizik Heat Conduction 
 */
void Sweep_Solve_1D_SteadyState_Annulus_Convective(double *const *const U_1, double *const *const U_2, double *const *const F, int sweep_sx, int sweep_ex, int sweep_sy, int sweep_ey,  int nx, int ny, int N, int s_x, int e_x, double r_i, double delta_r, double h, double k, double alpha ){
	// phi = 2pi in a complete revolution so 2pi/num_y = phi
	
	//____ Material Properties _____//
	// Pure Copper: Heat transfer textbook third edition john H. Lienhard IV 
//	const double k = 398;		//< W/m.degC	
	//double roh = 8954;	//< Kg/m^3
	//double C_p = 384; 	//< kJ/kg.degC
//	const double h = 13.14;	//< W/m^2K Convective heat transfer coefficient (Journal of Engineering Technology and Applied Sciences2019)
//	const double alpha = .0001157	//< m^2/s  	

	//double r_i = .002; //< m -> 4mm	//< 2D cross section of 18650 cell (image from journal of power sources 472 2020)
	//double r_o = .009; //< m -> 9mm
		
//	const double delta_r = (r_o-r_i)/((double)N-1); 	//< step size in r Note: p465 uses 0-M()
	//int delta_phi = (2*PI)/M;	//< phi is periodic with 2PI radians 	
	const double Int_heat_gen = 0.0; //(alpha*delta_t)/k; //DELETE LATER: need to check this function, also may be able to add delta_r^2 to save computation 
	 
	//There is an assumption that a>delta term p466

	const double T_Inf_ri = 20; //< Ambient temperature at inner radius
	const double T_Inf_ro = 20; //< Ambient temperature at outer radius

	const double Beta_0 = 1+(1-(1/(2*(r_i/delta_r))))*((delta_r*h)/k);	//< DELETE LATER Is delta an Int, need to Cast ?  	
	const double Gamma_0 = (1-1/(2*(r_i/delta_r)))*(delta_r/k)*(h*T_Inf_ri);
	const double Beta_N = 1+(1+(1/(2*((r_i/delta_r)+(N-1)))))*((delta_r*h)/k);	//< DELETE Later same as above + N ?	
	const double Gamma_N = (1+1/(2*((r_i/delta_r)+(N-1))))*(delta_r/k)*(h*T_Inf_ro); 

	if(s_x == 1){	//< Only processors with a Leftmost boundary
		for(int i=0; i<= ny; i++ ){
			U_2[i][1] = (2*U_1[i][2]+2*Gamma_0+delta_r*delta_r*Int_heat_gen)/(2*Beta_0*U_1[i][1]);	//< Update left column convective boundary   
			//^ See page 466, 12-82a	
		}
	}	

  	for(int i=sweep_sy; i<=sweep_ey; i++){	//< Sweep ignores ghost columns//
    		for(int j=sweep_sx; j<= sweep_ex; j++){  
      			U_2[i][j] = 0.5*((1-1/(2*((r_i/delta_r)+(double)s_x+j-1)))*U_1[i][j-1] + (1+1/(2*((r_i/(delta_r))+(double)s_x+j-1)))*U_1[i][j+1] + delta_r*delta_r*Int_heat_gen); 
			//^ See page 465 for formula, a->r_i and j-> global grid so (s_x+j), in book i=0->M, in our case j starts at 1 so to replicate say 0->N-1
		}		
  	}


	if(e_x == N){	//< Only processors with a Rightmost boundary	
		for(int i=0; i<= ny; i++){
			U_2[i][nx] = (2*U_1[i][nx-1]+2*Gamma_N+delta_r*delta_r*Int_heat_gen)/(2*Beta_N*U_1[i][nx]);		//< Update right column convective boundary
			//^ See page 466,12-82b
		}
	}
}


/*
 * \brief	This function performs a 2D grid sweep solve across an Annulus shaped cross section with 1D heat conduction
 * 		in the radial r direction (cartesian x equivalent). This is done with fixed left and right boundary condition
 * 		and does not involve time dependancy. 
 * \notes	- It may seem as though this function should be capable of 2D heat transfer however it is not. This would require
 * 		  a different calculation to take place and so here only 1D heat conduction is simulated.
 * 		- Accessing data with unit stride
 * 		- Formula: Page 465 of Osizik Heat Conduction 
 */
void Sweep_Solve_1D_SteadyState_Annulus_FixedBound(double *const *const U_1, double *const *const U_2, double *const *const F, int sweep_sx, int sweep_ex, int sweep_sy, int sweep_ey,  int nx, int ny, int N, int s_x, double r_i, double delta_r){

	double Int_heat_gen = 0.0; //(alpha*delta_t)/k; //DELETE LATER: need to check this function, also may be able to add delta_r^2 to save computation 
  	for(int i=sweep_sy; i<=sweep_ey; i++){	//< Sweep ignores ghost columns//
    		for(int j=sweep_sx; j<= sweep_ex; j++){    
			U_2[i][j] = 0.5*((1-1/(2*((r_i/delta_r)+(double)s_x+j-1)))*U_1[i][j-1] + (1+1/(2*((r_i/(delta_r))+(double)s_x+j-1)))*U_1[i][j+1] + delta_r*delta_r*Int_heat_gen); 				//^ See page 465 for formula, a->r_i and j-> global grid so (s_x+j), in book i=0->M, in our case j starts at 1 so to replicate say 0->N-1
		}		
  	}
}



/*
 *  \brief	This function checks to verify that the parameters selected for simulation with the
 * 		Time dependant explicit method meet the criteria for numerical stability. In essence, 
 * 		the finite difference solution solves one node at a time preventing imediate propogation 
 * 		of information through the grid.
 * 		If this propogation rate is slower than the rate at which heat transfer is occuring then the 
 * 		simulation cant keep up and the system becomes unstable and amplify error.	
 * \notes	- There is a maximum allowable time step to conform with stability requirements
 *		- (alpha*delta_t/delta_x^2)+(alpha*delta_+t/delta_y^2)<=1/2	
 * \params	Thermal diffusivity -> alpha = k/(C_p*roh) material dependant property
 *		delta_t -> step size in the time dimension
 *		delta_x & delta_y -> step sizes in respective spacial dimensions
 */

void Transient_1D_Stability_Check(double alpha, double delta_t, double delta_x){
	double r_x = (alpha*delta_t)/(delta_x*delta_x); 
	
	if(r_x > .5){	//< Stability check
       		char errormsg[150];
		sprintf(errormsg,"Error: The selected simulation parameters does not meet the stability criteria, occured in File: %s, Line: %d\n", __FILE__, __LINE__);
		perror(errormsg);
		exit(EXIT_FAILURE);	//< Critical failure, these values would lead to instability and waste computation time //     
	}
}

/*
 *  \brief	This function checks to verify that the parameters selected for simulation with the
 * 		Time dependant explicit method meet the criteria for numerical stability. In essence, 
 * 		the finite difference solution solves one node at a time preventing imediate propogation 
 * 		of information through the grid.
 * 		If this propogation rate is slower than the rate at which heat transfer is occuring then the 
 * 		simulation cant keep up and the system becomes unstable and amplify error.	
 * \notes	- There is a maximum allowable time step to conform with stability requirements
 *		- (alpha*delta_t/delta_x^2)+(alpha*delta_+t/delta_y^2)<=1/2	
 * \params	Thermal diffusivity -> alpha = k/(C_p*roh) material dependant property
 *		delta_t -> step size in the time dimension
 *		delta_x & delta_y -> step sizes in respective spacial dimensions
 */

void Transient_2D_Stability_Check(double alpha, double delta_t, double delta_x, double delta_y){
	double r_x = (alpha*delta_t)/(delta_x*delta_x); 
	double r_y = (alpha*delta_t)/(delta_y*delta_y);  
	
	if(r_x > .5 && r_y > .5){	//< Stability check
       		char errormsg[150];
		sprintf(errormsg,"Error: The selected simulation parameters does not meet the stability criteria, occured in File: %s, Line: %d\n", __FILE__, __LINE__);
		perror(errormsg);
		exit(EXIT_FAILURE);	//< Critical failure, these values would lead to instability and waste computation time //     
	}
}


/*
 * \brief	This function performs a 2D grid sweep solve across an Annulus shaped cross section with 1D transient conduction in a radial direction and internal heat generation.
 * \todo	- f(x,y,t) does this term exist (assumed zero?) Must check and understand where this term has gone, significant memory saving if its not needed
 *		- Int_heat_gen is a constant and should be calculated once and passed in
 *		- delta_x, delta_y should also be constants passed in indicating step size presumable relative to an input x,y dimension as opposed to M,N as we have
 *		- need a stability check (outside of this function)
 * \notes	- Accessing data with unit stride 
 * \todo	- This function is incomplete	
 */

void Sweep_Solve_1D_Transient_Annulus_Heat_Gen(double *const *const U_1, double *const *const U_2, double *const *const F, const double alpha, const int sweep_sx, const int sweep_ex, const int sweep_sy, const int sweep_ey, const double delta_r, const double delta_t, const double k, const double r_i, const int s_x){
	
	double Int_heat_gen = 0.0;//10*delta_t/k; //< heat generation is a spacial factor but is uniform in homogenous material...
  	for(int i=sweep_sy; i<=sweep_ey; i++){	//< Sweep ignores ghost columns//
    		for(int j=sweep_sx; j<= sweep_ex; j++){  
			U_2[i][j] = U_1[i][j]+.0001*.0001157*((1-1/(2*((r_i/delta_r)+(double)s_x+j-1)))*U_1[i][j-1] + (1+1/(2*((r_i/(delta_r))+(double)s_x+j-1)))*U_1[i][j+1] + delta_r*delta_r*Int_heat_gen); 
			
		//U_2[i][j] = U_1[i][j]+ alpha(((U_1[i][j-1]-2*U_1[i][j]+U_1[i][j+1])/(delta_r*delta_r))+(1/(r_i+((double)s_x-1+j-1)*delta_r))*((U_1[i][j+1]-U_1[i][j-1])/(2*delta_r)) + (1/k)*Int_heat_gen);
      		//U_2[i][j] = ((1-1/(2*((r_i/delta_r)+(double)s_x+j-1)))*U_1[i][j-1] + (1+1/(2*((r_i/(delta_r))+(double)s_x+j-1)))*U_1[i][j+1] + delta_r*delta_r*Int_heat_gen); 		
	//U_2[i][j]= U_1[i][j] + alpha*((1-(1/(2*((r_i/delta_r)+(double)s_x-1+j-1))))*U_1[i][j-1] -2*U_1[i][j]+ (1+(1/(2*((r_i/delta_r)+(double)s_x-1+j-1))))*U_1[i][j+1]+ (delta_r*delta_r*Int_heat_gen)/k);
		}		
  	}
}


void Sweep_Solve_1D_Transient_Annulus_FixedBound(double *const *const U_1, double *const *const U_2, double *const *const F, int sweep_sx, int sweep_ex, int sweep_sy, int sweep_ey,  int nx, int ny, int N, int s_x, double r_i, double delta_r){

	const double alpha= .0001157; // m^2/s
	const double delta_t = 1; // seconds  	
	

	double Int_heat_gen = 0.0; //(alpha*delta_t)/k; //DELETE LATER: need to check this function, also may be able to add delta_r^2 to save computation 
  	for(int i=sweep_sy; i<=sweep_ey; i++){	//< Sweep ignores ghost columns//
    		for(int j=sweep_sx; j<= sweep_ex; j++){    
			U_2[i][j] = U_1[i][j] +delta_t*alpha*((1-1/(2*((r_i/delta_r)+(double)s_x+j-1)))*U_1[i][j-1] + (1+1/(2*((r_i/(delta_r))+(double)s_x+j-1)))*U_1[i][j+1] + delta_r*delta_r*Int_heat_gen); 				//^ See page 465 for formula, a->r_i and j-> global grid so (s_x+j), in book i=0->M, in our case j starts at 1 so to replicate say 0->N-1
		}		
  	}
}



/*
 * \brief	Function measures the difference between two sequential itterations of a grid.
 * 		This can be used to an error measure to indicate when a method has converged
 * 		to within a desired tollerance. The function excludes a root|mean calculation
 * 		to allow local proc calculation of sum of differences before being combined
 * 		with other processor results to obtain global solution.
 * \notes	- The idea was to use this to calculate the norm of a matrix however would
 * 		  only work for a square matrix, calculation of rectangular matrix norms
 * 		  could be acheived with manipulation of this functiuon but not hold
 * 		  universal applicability and would require greater amounts of inter-processor
 * 		  communication.
 * 		- Should be accessing data with unit stride
 * 		- Can include a Root\mean calculation externally.
 * 		- Sweep over grids consideres only true grid values and ignores ghosts
 */
double Calc_grid_diff(double *const *const U_1, double *const *const U_2, int nx, int ny){
	double Diff = 0.0;
	double Sum_sq_diff = 0.0;	
	
  	for(int i=1; i<=ny; i++){			
    		for(int j=1; j<= nx; j++){  
			Diff = U_1[i][j] - U_2[i][j]; 	//< Calc change //      			
			Sum_sq_diff += Diff*Diff;	//< Sum of squared difference // 
    		}
  	}
	return(Sum_sq_diff);
}

/*
 * \brief	Function measures the difference between two sequential itterations of a grid and squares
 *		the difference to get a positive value which it divides by a similiarly squared Tollerance value.
 *		This will provide either 0 or some positive integer value for all grid values. If a matrix has converged
 *		to a desired tollerance at every single node all processors would return 0 for Sum_diff.
 * 		The function allows a local proc calculation, the result of which can be combined with other processor
 * 		results to obtain global solution.
 *
 * \notes	- Should be accessing data with unit stride
 * 		- This function was originally considered with the math.h's floor() command however was beleived to be more time intensive
 * 		  and to prevent vectorization. This would significantly impact the time to solve a problem as convergence
 * 		  needs to be checked at each itteration of the grid and at each point.
 * 		- Returns 0 if grids have converged and any other value otherwise
 * 		- The proposed solution may be suseptible to overflow, where the cast double is a value greater than the maximum allowable
 * 		  stored in an Integer. It is significant to note however that the numerical value of Sum_sq_diff is irrelevent as it 
 * 		  communicates whether all nodes have met a convergence criteria with an indicator of 0 as converged or any other number 
 * 		  as not-converged, which can include any positive/negative numbers that can result from undefined behaviour. 
 * \usage	- Include a check external to this function if value == 0 converged and use MPI_Reduce() to obtain global value.  
 */ 
double Calc_grid_diff_pointwise(double *const *const U_1, double *const *const U_2, const int nx, const int ny, double Tollerance){
	double Diff = 0.0;
	double Convergance = 0.0;
	Tollerance *= Tollerance;			//< Squared because Diff will be squared // 	
	
  	for(int i=1; i<=ny; i++){			
    		for(int j=1; j<= nx; j++){  
			Diff = U_1[i][j] - U_2[i][j]; 	//< Calc change //      			
			Diff = (Diff*Diff)/Tollerance;
			Convergance += (int)Diff;	//< Sum of squared difference // 
    		}
  	}
	return(Convergance);

}	
