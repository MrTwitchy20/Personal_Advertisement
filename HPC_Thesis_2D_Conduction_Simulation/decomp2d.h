/**
 * \filename	decomp2d.h
 * \brief	Header file for decomp2d.c
 * \author	F.OSuibhne
 * \date	07.09.21
 */

#ifndef decomp2d_h 
#define decomp2d_h
int MPI_Rect_Axis_Decomp(const int, int *, int *, const int, const int);
int MPI_Decomp1d(const int, const int, const int, int *, int *);
int MPI_Decomp2d(const int, const int, const int, int *, int *, int *, int *);
int MPI_Rec_Grid_Decomp2d(const int, const int, const int, const int, const int, const int, int *, int *, int *, int *);

#endif 
