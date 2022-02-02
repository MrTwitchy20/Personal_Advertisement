/**
 * \filename 	decomp2test.c
 * \brief	File is simply a test for demonstrating successfull decomposition of a 2d grid by printing possessed
 * 		grid start and end x,y positions for each relevent node/proc.	
 * \author	F.OSuibhne
 * \date	27.07.21 
 * \ToDo	- Should filename be read in by all procs? instead of distributed
 * 		- Requires input commandline inpute "filename"  
**/


#include<mpi.h>
#include<stdio.h>
#include<stdlib.h>
#include<string.h>

int MPI_Decomp2d(int, int, int, int *, int *, int *, int *);

int main(int argc, char **argv){
	int myid,nprocs;
	int s_x=0, s_y=0, e_x=0, e_y =0;
	int grid_size=10;
	char fname[200];
  	//int fname_len;
 
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &myid);
	MPI_Comm_size(MPI_COMM_WORLD,&nprocs);
	
	int Decomp_check = MPI_Decomp2d(grid_size, nprocs, myid, &s_x, &e_x, &s_y, &e_y);	
	if(Decomp_check == 1){	//Error check Decomp2d
		char errormsg[150];
		sprintf(errormsg,"Failed to Decompose matrix in File: %s, Line: %d\n", __FILE__, __LINE__);
		perror(errormsg);
		exit(EXIT_FAILURE);
	}
	printf("Im proc %d of %d sharing s_x = %d e_x = %d, s_y = %d, e_y =%d \n", myid, nprocs, s_x, e_x, s_y, e_y);	
	MPI_Barrier(MPI_COMM_WORLD);
/* 
  //___ Read in filename from commandline and distribute to all procs____//
  //NOTE: REQUIRES ERROR CHECKING INCASE NO_FILENAME ENTERED ALSO WHAT FILENAMES ARE NOT ACCEPTABLE? //
  if(myid == 0){
    	if( argc != 2 ){
    		fprintf(stderr,"Error: usage: %s <output filename>\n",argv[0]);
    		MPI_Abort(MPI_COMM_WORLD, 2);
    	}
    sprintf(fname, "%s", argv[1]);
    fname_len = strlen(fname);    //set filename
  }

  MPI_Bcast(&fname_len, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(fname, (fname_len+1), MPI_CHAR, 0, MPI_COMM_WORLD);
*/ 	
	if(argc != 2){
		fprintf(stderr,"Error: usage: %s <output filename>\n",argv[0]);
  		MPI_Abort(MPI_COMM_WORLD, 2);
  	}

  
  printf("rank: %d): fname: %s\n",myid,fname);
  //MPI_Barrier(MPI_COMM_WORLD);
  MPI_Finalize();

	return(0);
}
