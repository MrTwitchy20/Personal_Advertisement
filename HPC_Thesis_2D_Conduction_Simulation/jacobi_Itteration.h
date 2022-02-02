/**
 * \filename	jacobi_Itteration.h
 * \brief	Header file for jacobi_Itteration.c
 * \author	F.OSuibhe
 * \date	06.09.21
 */

#ifndef Jacobi_Itteration_h
#define Jaccobi_Itteration_h

void Exchange_Data_2D(int, int, double *const *const, MPI_Comm, int, int, int, int, int, int, int, int, MPI_Datatype);
//void Sweep_Solve_2D(double *const *const, double *const *const, double *const *const, int, int);
void Sweep_Solve_2D(double *const *const, double *const *const, double *const *const, int, int, int, int, int, int);
void Sweep_Solve_1D_SteadyState_Annulus_Convective(double *const *const, double *const *const, double *const *const, int, int, int, int, int, int, int, int, int, double, double, double, double, double);
void Sweep_Solve_1D_SteadyState_Annulus_FixedBound(double *const *const, double *const *const, double *const *const, int, int, int, int, int, int, int, int, double, double);
void Transient_1D_Stability_Check(double, double, double);
void Transient_2D_Stability_Check(double, double, double, double);
void Sweep_Solve_1D_Transient_Annulus_Heat_Gen(double *const *const, double *const *const, double *const *const, const double, const int, const int, const int, const int, const double, const double, const double, const double, const int);
void Sweep_Solve_1D_Transient_Annulus_FixedBound(double *const *const, double *const *const, double *const *const, int, int, int, int,  int, int, int, int, double, double);
double Calc_grid_diff(double *const *const, double *const *const, int, int);
double Calc_grid_diff_pointwise(double *const *const, double *const *const, const int, const int, double);
			
#endif 
