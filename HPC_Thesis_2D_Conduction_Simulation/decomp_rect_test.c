


#include<mpi.h>
#include<stdio.h>
#include<stdlib.h>
#include<string.h>

int MPI_Rect_Decomp(const int, int *, int *, const int, const int);
int MPI_Rec_Decomp2d(const int, const int, const int, const int, const int, const int, int *,int *,int *,int *);


int main(int argc, char **argv){

	int Total_No_procs;	// = 6;
	int Grid_dim_X = 30;
	int Grid_dim_Y = 10;
	int X_procs = 0;
	int Y_procs = 0;	
	
	int myid; 
	int s_x=0, s_y=0, e_x=0, e_y =0;
 
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &myid);
	MPI_Comm_size(MPI_COMM_WORLD,&Total_No_procs);
		 
	int test = MPI_Rect_Decomp(Total_No_procs, &X_procs, &Y_procs, Grid_dim_X, Grid_dim_Y);

//	printf("%d procs, %d columns , %d rows \n ", Total_No_procs, Grid_dim_X, Grid_dim_Y);
//	printf("%d procs in X dim : %d Procs in Y dim \n", X_procs, Y_procs);
	
	int Decomp_check = MPI_Rec_Decomp2d(Total_No_procs, X_procs, Y_procs, Grid_dim_X, Grid_dim_Y, myid, &s_x, &e_x, &s_y, &e_y);

	if(Decomp_check == 1){	//Error check Decomp2d
		char errormsg[150];
		sprintf(errormsg,"Failed to Decompose matrix in File: %s, Line: %d\n", __FILE__, __LINE__);
		perror(errormsg);
		exit(EXIT_FAILURE);
	}
	printf("Im proc %d of %d sharing s_x = %d e_x = %d, s_y = %d, e_y =%d \n", myid, Total_No_procs, s_x, e_x, s_y, e_y);	
	
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Finalize();
	return(0);
}
