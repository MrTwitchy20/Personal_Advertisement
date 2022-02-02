/**
 * \file        memory_alloc_test.c
 * \brief       Thesis:
 * \note        - Assumes use of gcc compillers with use of "malloc"
 * \author      F.OSuibhne
 * \date        30.06.21
 * \version     0.1
 **/

#include<stdio.h>
#include<string.h>
#include<stdlib.h>  //for test use of malloc

void print_sq_matrix(double *const *const, const int);
int alloc_2D_sq_matrix(double ***, const int);
int alloc_2D_rect_matrix(double ***, const int, const int);
void print_rect_matrix(double *const *const, const int, const int);
void Label_2D_matrix_positions(double *const *const, const int, const int);
void Initialize_2D_Grid_Values(double * const *const A, const int M, const int N, const int);

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
    /*
    // setting values (for testing only) //
    count = 0;  
    for (int i=0; i<M; i++){
        for(int j=0; j<N; j++){
            A[i][j] = count++;
        }
    }
    */
    //Label_2D_matrix_positions(A,M,N);
    Initialize_2D_Grid_Values(A,M,N,0.0);
    print_rect_matrix(A,M,N);

    free(*(&A[0]));
    free(*(&A)); 

    return(0);
}
