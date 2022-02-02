/**
 * \filename	memory_alloc.h
 * \brief	Header file for memory_alloc.c
 * \author	F.OSuibhne
 * \date	06.09.21
 */

#ifndef memory_alloc_h
#define memory_alloc_h
void Label_2D_matrix_positions(double * const *const, const int, const int);
void Initialize_2D_Grid_Boundries(double * const * const, const double, const double, const double, const double, const int, const int, const int, const int, const int, const int, const int, const int);	
void Initialize_2D_Grid_Ghost_Boundries(double * const * const, const double, const double, const double, const double, const int, const int, const int, const int, const int, const int, const int, const int);	
void Initialize_2D_Grid_Values(double * const *const, const int, const int, const double);
int alloc_2D_sq_matrix(double ***, const int);
int alloc_2D_rect_matrix(double ***, const int, const int);
void print_sq_matrix(double *const *const, const int);
void print_rect_matrix(double *const *const, const int, const int);   
#endif 
