/**
 * \filename	Parallel_write_to_file.h
 * \brief	Header file for Parallel_write_to_file.c
 * \author	F.OSuibhne
 * \date	07.09.21
 */

#ifndef Parallel_write_to_file_h
#define Parallel_write_to_file_h
void Parallel_Print_Int_To_File(const char *, const int, const int, const int, const int, const int, const int, const int, const int, const int, const int);
void Parallel_Print_Double_To_File(const char *, const int, const int, const int, const double, const int, const int, const int, const int, const int, const int);
void Parallel_Print_Double_To_File_wBoundary(const char *, const int, const int, const int, const double *, const int, const int, const int, const int, const int, const int);

#endif

// Parallel_Print_Double_To_Fule_wBoundary(); previously declared as Buff[] -> "double *" in main function,  
