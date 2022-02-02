/**
 * \filename	Dynamic_write_to_filei.c
 * \breif	MPI_write_to_file accepts a single pointer (double A[]), this file will use dynamic 2D arrys of type double pointer (double A[][])
 * \author	F.OSuibhne
 * \date	09.08.21
 * \version	V0.1
 **/

#include<stdio.h>
#include<stdlib.h>  //for test use of malloc
#include<mpi.h>
#include<string.h>

//int MPI_Decomp2d(int, int, int, int *, int *, int *, int *);
void Parallel_Print_Double_To_File( const char *, const int , const int, const int, double *, const int, const int, const int, const int, const int, const int);

int alloc_2D_rect_matrix(double ***, const int, const int);
void print_rect_matrix(double *const *const, const int, const int);
void Initialize_2D_Grid_Values(double * const *const, const int, const int, const double);
void Label_2D_matrix_positions(double *const *const, const int, const int);

int MPI_Rect_Decomp(const int, int *, int *, const int, const int);
int MPI_Rec_Decomp2d(const int, const int, const int, const int, const int, const int, int *,int *,int *,int *);




int main(int argc, char **argv){
	int M = 10; //rows
	int N = 6; //columns
	double ** Grid=NULL;
	
	int myid,nprocs;
	int s_x=0, s_y=0, e_x=0, e_y =0;
	int grid_size = N;
	char fname[200];
  	int fname_len;
 
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &myid);
	MPI_Comm_size(MPI_COMM_WORLD,&nprocs);
	
//	int Decomp_check = MPI_Decomp2d(grid_size, nprocs, myid, &s_x, &e_x, &s_y, &e_y);	
//__________ Changed __________//
	int X_procs = 0;
	int Y_procs = 0;	
		 
	int test = MPI_Rect_Decomp(nprocs, &X_procs, &Y_procs, N, M);	
	int Decomp_check = MPI_Rec_Decomp2d(nprocs, X_procs, Y_procs, N, M, myid, &s_x, &e_x, &s_y, &e_y);
//__________ End of Changed ____//


	if(Decomp_check == 1){	//Error check Decomp2d
		char errormsg[150];
		sprintf(errormsg,"Failed to Decompose matrix in File: %s, Line: %d\n", __FILE__, __LINE__);
		perror(errormsg);
		exit(EXIT_FAILURE);
	}
	printf("Im proc %d of %d sharing s_x = %d e_x = %d, s_y = %d, e_y =%d \n", myid, nprocs, s_x, e_x, s_y, e_y);	
	MPI_Barrier(MPI_COMM_WORLD);

  int Num_x = e_x-s_x+1;
  //int M=Num_x; //rows
  int Num_y = e_y-s_y+1;
  //int N=Num_y; //columns
  int MY_BUFSIZE = Num_x*Num_y;

	//__ Deallocation_test __//
	alloc_2D_rect_matrix(&Grid, M, N);// allocate memmory
	Label_2D_matrix_positions(Grid,M,N); //Numbers matrix positions
	print_rect_matrix(Grid,M,N);	  // prints matrix
	
//  double buf_double[MY_BUFSIZE];
	for(int i=0; i<M; i++){
		for(int j=0; j<N; j++){
			Grid[i][j]= (double) myid;
		}
	}

	Initialize_2D_Grid_Values(Grid, M, N, (double) myid);
	print_rect_matrix(Grid,M,N);	  // prints matrix

  	if( argc != 2 ){
  		fprintf(stderr,"Error: usage: %s <output filename>\n",argv[0]);
   		MPI_Abort(MPI_COMM_WORLD, 2);
   	}

    	sprintf(fname, "%s", argv[1]); 
  	printf("rank: %d): fname: %s\n",myid,fname);
  	MPI_Barrier(MPI_COMM_WORLD);
 
	Parallel_Print_Double_To_File( fname, Num_x, Num_y, MY_BUFSIZE, *(&Grid[0]), grid_size, myid, s_x, e_x, s_y, e_y);
 		
	free(*(&Grid[0]));
	free(*(&Grid));	
 
	MPI_Barrier(MPI_COMM_WORLD);
  	MPI_Finalize();
	
	return(0);
}


/**
 * \brief   Function allocates a 2D rectangular block of memory of size MxN using an initial 
 *          pointer input as the starting adress.
 * \note    - do we want to use malloc here or alloc, its desireable to have contiguous memory ?
 *          - What cluster/compillers will this code run on, gcc or intel ?
 * \todo    - Memory allocated will need to be free'd externally to this function
 */
int alloc_2D_rect_matrix(double ***A, const int M, const int N){

    *A = malloc(M*sizeof(double*));             //< allocates an array of pointers pointed to by A //
    double *data = malloc(M*N*sizeof(double));  //< allocates the 2D matrix of memory // 
    if (*A == NULL || data == NULL){            //< Malloc error checking //
        perror("Error allocating memory in alloc_2D_matrix");
        return(-1);                             //< should include further error checking in main
                                                // this is a fatal error, should close program //     
    }

    for (int i=0; i<M; i++){	//one of these should be N not M..
        (*A)[i] = (data + N*i); //< Setting pointer array to point to start of columns??//
    }

    return(0);  //< Indicates sucessfull allocation //  
}


/*
 * \brief   Function used to set the data values of a 2D matrix to 
            their respective position values from 0->1-M*N
 * \notes      
 * \illustration    -e.g. for an M=3, N=4 matrix:
 *                      | 0  1   2   3  |
 *                      | 4  5   6   7  |
 *                      | 8  9  10  11  |
 */
void Label_2D_matrix_positions(double *const *const A, const int M, const int N){
    int count = 0;  
    for (int i=0; i<M; i++){
        for(int j=0; j<N; j++){
            A[i][j] = count++;
        }
    }
}

/*
 * \brief   Function Initialises 2D Grid data to a specified starting value. 
 */
void Initialize_2D_Grid_Values(double * const *const A, const int M, const int N, const double Initial_val){
    for (int i=0; i<M; i++){    //< rows //
        for(int j=0; j<N; j++){ //< columns  //
            A[i][j] = Initial_val;
        }
    }
}

/*
 *\breif    Function prints a 2D rectangular matrix to command line
 */
void print_rect_matrix(double *const *const B, const int x, const int y){
    
    for(int i=0; i<x; i++){
        printf("|\t");
        for(int j=0; j<y; j++){
            if(j==0){   //< If start of row place //
                printf("%0.1lf",B[i][j]);
            }
            else if(j==y-1){
                printf("%0.1lf", B[i][j]);   //< print value at matrix co-ordinate //
            }
            else{   
                printf("%0.1lf", B[i][j]);   //< print value at matrix co-ordinate //
            }
            printf("\t");
        }
        printf("|\n");
    }
    printf("\n");
}

/*
 * \ brief	Function utilises MPI IO commands to allow each node/proc to print
 *		their respective section of the grid to a binary file in parallel.  	
 * \ notes	- Function requires use of openmpi commands
 *		- MPI_COMM_WORLD is a global communicator and as such cant be passed in with 
 *		  "Function(MPI_Comm MPI_COMM_WORLD)". Not sure of usage case were specified 
 *		  communicators would need to print so will ignore however could present as an issue
 *		  if desired later to specify custom comm. Safest bet might be to force a new
 *		  MPI_COMM == MPI_COMM_WORLD
 *	
 *		Benefits: 
 *		- Saving to a single file by all procs similtaniously avoiding alternative solutions which
 *		  would essentially serialize the process wasting significant resource time as nodes sit idle
 *		  waiting for their turn to print or communicate to a master node for printing.
 * 		- This is a more scalable solution, it would allow for larger grid sizes. Nodes can maxamise
 *		  memory usage without the need to reserve memory for communication buffers.
 *		
 *		Disadvantage:
 *		- The output file resulting from this function is not immediatly user readable and requires a further
 *		  process to interpret/reformat the data.
 *
 * \ ToDo	- MPI_File_set_view() using "natve", this should be changed as will cause issues if run on disimilar device.
 *		- Should an MPI_Barrier be placed at end of function ? //< I suspect not nessisary
 *		- Include error checking 
 */
void Parallel_Print_Double_To_File(const char *fname, const int Num_x, const int Num_y, const int MY_BUFSIZE, double *buf, const int grid_size, const int myid, const int s_x, const int e_x, const int s_y, const int e_y){

  MPI_Aint lb, extent; //is this still correct or should it be a double
  MPI_Datatype etype, filetype, contig;
  MPI_Offset disp;
  MPI_File fh;
 
  //fh is refered to as the file handle I imagine it as a file pointer
  MPI_File_open(MPI_COMM_WORLD, fname, MPI_MODE_CREATE | MPI_MODE_RDWR,
  MPI_INFO_NULL, &fh);	//< All procs open file and have a filepointer //
			// MPI_MODE_WRONLY-- write only //
			// MPI_MODE_APPEND //
			// Should have fatal errors, however FIle_Open_Error=?? //
 
  MPI_Type_contiguous(Num_x, MPI_DOUBLE, &contig); //< Proc dependant datatype is created //
  
  extent = grid_size*sizeof(double);       //< contigous datatypes are spaced apart by grid_size //
  lb = 0;       			//< "New lower bound of datatype (address integer)"
  MPI_Type_create_resized(contig, lb, extent, &filetype);  //< contiguous types, no offset, of 6 doubles ?
  MPI_Type_commit(&filetype);

  disp = (s_y-1)*extent + (s_x-1)*sizeof(double); // displacement your offset from 0 + row displacement
 
  etype = MPI_DOUBLE;
  MPI_File_set_view(fh, disp, etype, filetype, "native", MPI_INFO_NULL);	//< SHOULD WE BE STATING NATIVE OR LOWER ELISION?
  MPI_File_write(fh, buf, MY_BUFSIZE, MPI_DOUBLE, MPI_STATUS_IGNORE);	//< writes proc specific buffer to file //
  MPI_File_close(&fh);
}
