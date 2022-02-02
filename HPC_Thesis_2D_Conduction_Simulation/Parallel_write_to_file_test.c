/**
 * \filename 	Parallel_write_to_file_test.c
 * \brief	Contains a main function which utilises the Parallel_write_to_file to demonstrate successfull function.	
 * \author	F.OSuibhne
 * \date	27.07.21 
 * \Note	- Requires input command line input "filename"  
**/


#include<mpi.h>
#include<stdio.h>
#include<stdlib.h>
#include<string.h>

int MPI_Decomp2d(int, int, int, int *, int *, int *, int *);
//int alloc_2D_rect_matrix(double ***, const int, const int);
void Parallel_Print_Int_To_File(const  char *, const int , const int, const int, const int *, const int, const int, const int, const int, const int, const int);
void Parallel_Print_Double_To_File_wBoundary( const char *, const int , const int, const int, const double *, const int, const int, const int, const int, const int, const int);
 

int main(int argc, char **argv){
	
	int myid,nprocs;
	int s_x=0, s_y=0, e_x=0, e_y =0;
	int grid_size=10;
	char fname[200];
  	int fname_len;
 
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
	printf("Im proc %d of %d sharing s_x = %d e_x = %d, s_y = %d, s_y =%d \n", myid, nprocs, s_x, e_x, s_y, e_y);	
	MPI_Barrier(MPI_COMM_WORLD);

  int Num_x = e_x-s_x+1;
  int Num_y = e_y-s_y+1;
  int MY_BUFSIZE = (Num_x+2)*(Num_y+2);
  //int buf_int[MY_BUFSIZE];	//< wont be how array is allocated
  double buf_double[MY_BUFSIZE];
	//double ** buf_double = NULL;
	//alloc_2D_rect_matrix(&buf_double, Num_x, Num_y); //testing with square matrix

  //___ Read in filename from commandline and distribute to all procs____//
  //NOTE: REQUIRES ERROR CHECKING IN CASE NO_FILENAME ENTERED ALSO WHAT FILENAMES ARE NOT ACCEPTABLE? //
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
  printf("rank: %d): fname: %s\n",myid,fname);
  MPI_Barrier(MPI_COMM_WORLD);

  int count=1;
  //____________ Fill Buffer ________//  
  for(int i=0; i < MY_BUFSIZE; i++){
  	// buf_int[i] = (int)myid;   //used to test with multiple procs
        buf_double[i] = (double)i; //(double)myid;
  	printf("%.2lf \t", buf_double[i]); 
	
	if(count==12){
		printf("\n");
		count=0;
	} 
 	count+=1;
  }
  MPI_Barrier(MPI_COMM_WORLD);
  //____________ Print _____________ //
  //Parallel_Print_Int_To_File( fname, Num_x, Num_y, MY_BUFSIZE, buf_int, grid_size, myid, s_x, e_x, s_y, e_y);  
  Parallel_Print_Double_To_File_wBoundary( fname, Num_x, Num_y, MY_BUFSIZE, &buf_double[Num_x+3], grid_size, myid, s_x, e_x, s_y, e_y);
 
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Finalize();
	
	//free(*(&buf_double[0]));	//free allocated 2D memory arrays
	//free(*(buf_double));
	
	return(0);
}
