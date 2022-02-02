/**
 * \file        memory_alloc.c
 * \brief       Thesis:
 * \note        - Assumes use of gcc compillers with use of "malloc"
 * \author      F.OSuibhne
 * \date        30.06.21
 * \version     0.1
 **/


#include<stdio.h>
#include<string.h>
#include<stdlib.h>  

/*
void print_sq_matrix(double *const *const, const int);
int alloc_2D_sq_matrix(double ***, const int);
int alloc_2D_rect_matrix(double ***, const int, const int);
void print_rect_matrix(double *const *const, const int, const int);
void Label_2D_matrix_positions(double *const *const, const int, const int);
void Initialize_2D_Grid_Values(double * const *const A, const int M, const int N, const int );

int main(void){

    int N = 6;  //< row or column ?//
    int M = 15;
    double **A = NULL; //< Declare unknown values as NULL//

    alloc_2D_sq_matrix(&A, N);   //< Passing Address of A //
    
    // setting values (for testing only) //
    int count = 0;  
    for (int i=0; i<N; i++){
        for(int j=0; j<N; j++){
            A[i][j] = count++;
        }
    }

    print_sq_matrix(A, N);  //< print Matrix A of dimension N //
    
    free(*(&A[0]));
    free(*(&A)); 
    
    alloc_2D_rect_matrix(&A,M,N);
    
    // setting values (for testing only) //
    //count = 0;  
    //for (int i=0; i<M; i++){
    //    for(int j=0; j<N; j++){
    //       A[i][j] = count++;
    //    }
    //}
    
    //Label_2D_matrix_positions(A,M,N);
    Initialize_2D_Grid_Values(A,M,N);
    print_rect_matrix(A,M,N);

    free(*(&A[0]));
    free(*(&A)); 

    return(0);
}
*/

/*
 * \brief   Function used to set the data values of a 2D matrix to 
 *          their respective position values from 0->1-M*N
 * \notes      
 * \illustration    -e.g. for an M=4, N=3 matrix:
 *                      | 0  1   2   3  |
 *                      | 4  5   6   7  |
 *                      | 8  9  10  11  |
 */
void Label_2D_matrix_positions(double * const *const A, const int M, const int N){
    int count = 0;  
    for (int i=0; i<M; i++){
        for(int j=0; j<N; j++){
            A[i][j] = count++;
        }
    }
}


/*
 * \breif	Function initialises the grid boundaries distributed over multiple
 * 		processors. The grid boundaries are set to the edges of the Primary
 * 		Matrix which is held in portions on each processor, these are not the
 * 		same as the ghost columns that each processor will generate 
 * \notes	- The corner nodes of Primary matrix where two boundaries meet should 
 * 		  be irrelevent to calculation and in writing these boundaries the order
 * 		  of if statements in this function will dictate the arbitrary value these
 * 		  hold.
 * 		- ERROR: Due to the sequential manner these grids are allocated a corner nodes are overwritten
 * 		based on order of if statements, this is not typically an issue however when it is desired to 
 * 		utilize a periodic boundary these corner nodes become valid boundary nodes and will cause failure
 * 		based on the order of call 
 *
 * \illustration    An M=4, N=3 Primary Matrix:
 * 		   
 * 		     Primary Matrix        Proc 0    Proc 1     Proc 2     Proc 3
 *                  | 0   1   2   3 | -> |0 0 0 0| |0 0 0 0| |0  0  0 0| |0  0  0 0|
 *                  | 4   5   6   7 | -> |0 0 1 0| |0 2 3 0| |0  8  9 0| |0 10 11 0|
 *                  | 8   9  10  11 | -> |0 4 5 0| |0 6 7 0| |0 12 13 0| |0 14 15 0|
 *                  |12  13  14  15 | -> |0 0 0 0| |0 0 0 0| |0  0  0 0| |0  0  0 0|
 *					 
 *		     Proc 0     Proc 1    Proc 2      Proc 3	  Primary Matrix	  
 *		 -> |0 0 0 0| |0 0 0 0| |0  0 0 0| |0  0  0 0| -> |  0 U  U  3 |
 *               -> |0 0 U 0| |0 U 3 0| |0  L 9 0| |0 10  R 0| -> |  L 5  6  R |
 *               -> |0 L 5 0| |0 6 R 0| |0 12 D 0| |0  D 15 0| -> |  L 9 10  R |
 *               -> |0 0 0 0| |0 0 0 0| |0  0 0 0| |0  0  0 0| -> | 12 D  D 15 |
 */
void Initialize_2D_Grid_Boundries(double * const * const A, const double left, const double right, const double up, const double down, const int s_x, const int e_x, const int s_y, const int e_y, const int M, const int N, const int dim_y, const int dim_x){	


	if(s_y == 1){	//< Set top boundary
		for(int i=1; i< dim_x+1; i++){
			A[1][i] = up;
		}
	}

	if(e_y == M){	//< Set bottom boundary
		for(int i=1; i< dim_x+1; i++){
			A[dim_y][i] = down;
		}
	}

	if(s_x == 1){ 	//< Set top row	
		for(int i=1; i < dim_y+1; i++ ){
			A[i][1] = left; 
		}
	}
	
	if(e_x == N){	//< Set right boundary
		for(int i=1; i< dim_y+1; i++){
			A[i][dim_x] = right;
		}
	}
/*
	if(s_y == 1){	//< Set top boundary
		for(int i=1; i< dim_x+1; i++){
			A[1][i] = up;
		}
	}

	if(e_y == M){	//< Set bottom boundary
		for(int i=1; i< dim_x+1; i++){
			A[dim_y][i] = down;
		}
	}*/
}


/*
 * \breif	Function initialises the grid boundaries distributed over multiple
 * 		processors into the outermost ghost column/rows. The grid boundaries
 * 		are set to the edges of the Primary Matrix which is held in portions
 * 		on each processor, these are ghosts which are not exchanged.  
 * \notes	- The corner nodes of Primary matrix where two boundaries meet should 
 * 		  be irrelevent to calculation and in writing these boundaries the order
 * 		  of if statements in this function will dictate the arbitrary value these
 * 		  hold. 
 *		- This is Legacy code and should not be used: It can be used for troubleshooting
 *		however it was initially used to store outermost grid boundaries in unused outermost 
 *		ghost columns/rows. This was changed however in order to accomodate a print function that
 *		could include boundary nodes.
 * \illustration    An M=4, N=3 Primary Matrix:
 * 		   
 * 		     Primary Matrix        Proc 0    Proc 1     Proc 2     Proc 3
 *                  | 0   1   2   3 | -> |0 0 0 0| |0 0 0 0| |0  0  0 0| |0  0  0 0|
 *                  | 4   5   6   7 | -> |0 0 1 0| |0 2 3 0| |0  8  9 0| |0 10 11 0|
 *                  | 8   9  10  11 | -> |0 4 5 0| |0 6 7 0| |0 12 13 0| |0 14 15 0|
 *                  |12  13  14  15 | -> |0 0 0 0| |0 0 0 0| |0  0  0 0| |0  0  0 0|
 *					 
 *		     Proc 0     Proc 1    Proc 2      Proc 3	  Primary Matrix	  
 *		 -> |0 0 0 0| |0 0 0 0| |0  0 0 0| |0  0  0 0| -> |  0 U  U  3 |
 *               -> |0 0 U 0| |0 U 3 0| |0  L 9 0| |0 10  R 0| -> |  L 5  6  R |
 *               -> |0 L 5 0| |0 6 R 0| |0 12 D 0| |0  D 15 0| -> |  L 9 10  R |
 *               -> |0 0 0 0| |0 0 0 0| |0  0 0 0| |0  0  0 0| -> | 12 D  D 15 |
 */
void Initialize_2D_Grid_Ghost_Boundries(double * const * const A, const double left, const double right, const double up, const double down, const int s_x, const int e_x, const int s_y, const int e_y, const int M, const int N, const int dim_y, const int dim_x){	

	if(s_x == 1){ 	//< Set top row	
		for(int i=0; i < dim_y+2; i++ ){
			A[i][0] = left; 
		}
	}
	
	if(e_x == N){	//< Set right boundary
		for(int i=0; i< dim_y+2; i++){
			A[i][dim_x+1] = right;
		}
	}

	if(s_y == 1){	//< Set top boundary
		for(int i=0; i< dim_x+2; i++){
			A[0][i] = up;
		}
	}

	if(e_y == M){	//< Set bottom boundary
		for(int i=0; i< dim_x+2; i++){
			A[dim_y+1][i] = down;
		}
	}
}


/*
 * \brief   	Function Initialises 2D Grid data to a specified starting value. 
 */
void Initialize_2D_Grid_Values(double * const *const A, const int M, const int N, const double Initial_val){
	for (int i=0; i<M; i++){    //< rows //
        	for(int j=0; j<N; j++){ //< columns //
            	A[i][j] = Initial_val;
       		}
    	}
}


/**
 * \brief   Function allocates a 2D square block of memory of size NxN using the starting address pointer input
 * \note    - Use of alloc here instead could initialise to zero 
 * 	    - Will require alteration depending on use of Intel or gcc compiller
 * 	    - Returns -1 on failure and - on success.
 * \todo    - Memory allocated will need to be free'd externally to this function
 *
 * \usage   A pointer to a pointer can be made externally to function call
 * 	    "double **A=NULL;" and address of pointer passed into function "&A".
 * 	    This should enable calling or referal to any node in matrix using
 * 	    "A[M][N]". Once matrix is no longer needed two calls must be made to prevent 
 * 	    memory leak "free(*(&A[0]));" followed by "free(*(&A));"	
 **/
int alloc_2D_sq_matrix(double ***A, const int N){

	*A = malloc(N*sizeof(double*)); 
	double *data = malloc(N*N*sizeof(double));  

  	if(*A == NULL || data == NULL){   //< Malloc error checking //
       		char errormsg[150];
		sprintf(errormsg,"Error allocating memory in alloc_2d_matrix, occured in File: %s, Line: %d\n", __FILE__, __LINE__);
		perror(errormsg);
		exit(EXIT_FAILURE);	//< This is a critical failure//
		return(-1); // Unessisary      
    	}

	for(int i=0; i<N; i++){
        	(*A)[i] = (data + N*i);
    	}
   	return(0);  //< Indicate sucessfull allocation //  
}


/**
 * \brief   Function allocates a 2D rectangular block of memory of size MxN using an initial 
 *          pointer input as the starting adress.
 * \note    - Use of alloc here instead could initialise to zero
 *          - Will require alteration depending on use of intel or gcc compiller
 *          - Returns -1 on failure and 0 on success. 
 * \todo    - Memory allocated will need to be free'd externally to this function
 *
 *\usage   A pointer to a pointer can be made externally to function call
 * 	    "double **A=NULL;" and address of pointer passed into function "&A".
 * 	    This should enable calling or referal to any node in matrix using
 * 	    "A[M][N]". Once matrix is no longer needed two calls must be made to prevent 
 * 	    memory leak "free(*(&A[0]));" followed by "free(*(&A));"	
 **/
int alloc_2D_rect_matrix(double ***A, const int M, const int N){

	*A = malloc(M*sizeof(double*));             //< allocates an array of pointers pointed to by A //	
	double *data = malloc(M*N*sizeof(double));  //< allocates the 2D matrix of memory // 	
	
	if(*A == NULL || data == NULL){            //< Malloc error checking //
       		char errormsg[150];
		sprintf(errormsg,"Error allocating memory in alloc_2d_matrix, occured in File: %s, Line: %d\n", __FILE__, __LINE__);
		perror(errormsg);
		exit(EXIT_FAILURE);	//< This is a critical failure//
		return(-1); // Unessisary        
	}

    	for (int i=0; i<M; i++){
        	(*A)[i] = (data + N*i); //< Setting pointer array to point to start of columns??//
    	}
	return(0);  //< Indicates sucessfull allocation //  
}


/*
 *\breif    Function prints a formatted 2D square matrix to command line
 */
void print_sq_matrix(double *const *const B, const int dim){
    for(int i=0; i<dim; i++){
        printf("|\t");
        for(int j=0; j<dim; j++){
            if(j==0){   //< If start of row place //
                printf("%0.2lf",B[i][j]);
            }
            else if(j==dim-1){
                printf("%0.2lf", B[i][j]);   //< print value at matrix co-ordinate //
            }
            else{   
                printf("%0.2lf", B[i][j]);   //< print value at matrix co-ordinate //
            }
            printf("\t");
        }
        printf("|\n");
    }
    printf("\n");
}


/*
 *\breif    Function prints a formatted 2D rectangular matrix to command line
 */
void print_rect_matrix(double *const *const B, const int M, const int N){
    
    for(int i=0; i<M; i++){
        printf("|\t");
        for(int j=0; j<N; j++){
            if(j==0){   //< If start of row place //
                printf("%0.2lf",B[i][j]);
            }
            else if(j==N-1){
                printf("%0.2lf", B[i][j]);   //< print value at matrix co-ordinate //
            }
            else{   
                printf("%0.2lf", B[i][j]);   //< print value at matrix co-ordinate //
            }
            printf("\t");
        }
        printf("|\n");
    }
    printf("\n");
}
