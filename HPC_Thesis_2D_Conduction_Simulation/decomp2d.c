/**
 *\filename	decomp2d.c
 *\breif	This file contains functions for distributing the number of processors
 *		in balanced manner in 2 dimensions. It also contains further functions for 
 *		decomposing or distributing a grid across the processors given to each dimension.
 *\author	F.OSuibhne
 *\date		12-04-21
 */

#include<stdio.h>
#include<math.h>

/*
 * \brief	This Function is run once by a single proc and takes input dimensions of a rectangular grid to be decomposed, using
 * 		ratios of X,Y dimensions it determines how best to distribute allocated quantity of processors assigned
 * 		to run the job.
 * 
 * \note	- Not to be mixed up with MPI_Rect_Grid_Decomp2d();
 * 		- The function is trying to optimise the prescribed number of processors over a fixed grid size, this
 * 		  does not nessisarily translate to the optimised number of processors however. There is an oversimplication here
 * 		  as factors such as memory access may favour one axis over another in performance or a smaller no. procs.
 * 		  In short this evenly divides the work to be done however does not ensure each proc can perform the
 * 		  task designated equivalently fast to one another or to fewer.
 * 		- The function wont allow distribution of prime numbers (with exception to 1) of procs as these are only 1D distributable.
 * 		- It checks to ensure sufficient data is available for distribution across proposed qty of procs
 * 		- The returned value indicates an error if 0 and sucess if 1.
 * 		- Writing this function was unessisary as there exists an MPI function called MPI_Dims_create(); which 
 * 		  likely features greater flexibility see the following link: https://www.open-mpi.org/doc/v3.1/man3/MPI_Dims_create.3.php
 * 		  which provides a good description of its function. It would be better practice to use the existing function. 
 *
 * \todo	- A nice addition could include a suggestor that will indicate the next smallest number of procs that can be distributed,
 * 		this shouldnt be difficult simply cause function to loop back around and provide it with reduced No_procs.	
 */
int MPI_Rect_Axis_Decomp(const int No_procs, int * No_procs_X, int * No_procs_Y, const int X_dim, const int Y_dim){
	
	if( No_procs > X_dim*Y_dim ){		//< Number of procs allocated cant be less than the qty of data. //
		printf(" Error: There is insufficient data to distribute this job across %d Processors! \n ", No_procs);	
		return(1);
	}	
	
	double Grid_ratio = (double)Y_dim/(double)X_dim;
	int Largest_multiple = No_procs;	//< initialised to largest number iteself (always true for if statement)
	double result_double = 0;
	int result_int = 0;
	int Activate = 1;			//< Acts as a switch for an if statement
	double Current_ratio = 1.0;
	double Proposed_ratio = 1.0;
	int Best_result = 1;
	
	//printf("Grid ratio %lf \n", Grid_ratio);

	for(int factor=2; factor <= Largest_multiple && factor < No_procs; factor++){	// stops when the largest multiple of that number is obtained, or in event its a prime stops before
											// it divides into istelf (No_procs)
		result_double = No_procs/ (double)factor; 
		result_int = No_procs/factor;
		//printf("result_double= %lf result_int= %lf \n", result_double, (double)result_int);
		
		if(result_double == (double) result_int){ 	//< Comparing the rational part of same number to determine if result_double is rational
							  	// and inferring that the current factor is a multiple of No_procs.
								// Note: A better solution would have been to use modulus or floor/ceil commands
 
			if( Activate ){ 			//< Run once only on first true factor //
				Largest_multiple = result_int;	//< Sets the largest possible rational multiple for loop
								// also indicating non-prime 
				Activate = 0; 			//< Lock if Statment closed
			}		
			
			Proposed_ratio = factor/result_double;
			if( fabs(Proposed_ratio - Grid_ratio) <= fabs(Current_ratio - Grid_ratio)){ //< Ratio of new factor better approximates the ratio than previous keep it
				Current_ratio = Proposed_ratio;
				Best_result = result_int;
			}		
		}
	}		
	
	if( Largest_multiple == No_procs && No_procs!= 1){ 		//< No other multiples found
		printf("Error: No. of Procs cant be a prime number as its not possible to decompose in 2D plane \n");
		return(1);
	}	
	else{
		*No_procs_X = Best_result;		//< Grid_ratio is initialised so as to obtain No_procs for X axis. 
		*No_procs_Y = No_procs/Best_result;	//< As No_procs_X is now a known factor, remaining procs are Y axis.
		return(0);	//< Exit success //
	}
}


/**
 * \brief	Function takes an array size n and a number of processes
 * 		p to evenly divide and distribute the array between 
 * 		processes differing by at most 1 array positon. 
 * 		This acheives a load ballanced decomposition of 
 * 		an input array and returns respective start s and 
 * 		end e array position values for each
 *		respective process myid.
 * \notes	- Assumes myid values are between 0 - (p-1)
 *		- Returns 0 on success or 1 on failure
 *		- Returns failure if there are more processes
 *		  than vector positions as this would not be balanced.
 *		- Output incompatability s=1 => array[0]  
 */
int MPI_Decomp1d(const int n, const int p, const int myid, int *s, int *e){

	if(myid+1>p){		//< Error checking myid is within range//
		printf("myid %d is greater than the number of processors %d \n", myid,p);
		return(1);	//< Error //	
	}
	
	int whole = n/p; 	//< determines number of even divisions //
	int remainder = n-whole*p; //< calculates excess n vector positions to add to n processes //

	if(remainder == 0){	// array is perfectly devisable by processes
		*s = (myid*whole) + 1;
		*e = *s + whole-1;
	}
	else if(n<p){		//< More processors than array positions//
		printf("%d more processors than can be utilised for array \n",(p-n));
		return(1);	//< Error//
	}
	else{			//< Array divides unevenly into processes //
	
		if(myid>remainder-1){	//< Myid process above processes with additional array positions 
			*s = (remainder)*(whole+1)+1+(myid-(remainder))*whole;	//< Processes with extra vector value + processes without //	
			*e = *s+whole-1;
		}		
		else{		//< Myid sits within processes which get additional vector prosition for balancing//
			*s = (myid*(whole+1))+1;
			*e = *s + (whole+1)-1; //< (whole+1) counts additional vector position //
		}
	}
	return(0);		//< success //
}


/**
 * \brief	Function takes a grid size n and number of
 *		processes p to evenly divide and distribute
 *		the matrix between processos. This acheives a load 
 *		ballanced decomposition of an input matrix and 
 *		returns respective start s and end e array position
 *		values for each respective dimension & process myid.
 * 
 * \notes	- Assumes a square grid n*n
 *		- Assumes myid values are between 0 - (p-1)
 *		- Returns 0 on success or 1 on failure
 *		- Returns failure if more processors than vector
 *		  positions as this would not be balanced.
 *		- Output incompatability s=1 => array[0] 
 * \illustration
 *		- The distribution of an 2D x,y matrix with
 *		  a gridsize=5, across 9 processors
 *
 * 		 ____________________    proc:
 * 	row 1	|__1__2_|__3__4_|__5_|__  3
 * 	row 2 	|  6  7 |  8  9 | 10 |
 * 	row 3	|_11_12_|_13_14_|_15_|__  6
 *	row 4	| 16 17 | 18 19 | 20 |
 *	row 5	|_21_22_|_23_24_|_25_|__  9
 *		|	|	|    |
 *	column:	   1  2    3  4    5
 *	proc:	    7       8      9  
 *
 *		This code should find for myid=4 is designated centre of grid
 *		and s_x=3, e_x=4, s_y=3, e_y=4 
 *
 */
int MPI_Decomp2d(int n, int p, int myid, int *s_x, int *e_x, int *s_y, int *e_y){
	
	if( p - sqrt(p)*sqrt(p) != 0){			//< Error checking //
		printf("myid: %d, number of processors %d entered is not a power of 2 \n",myid,p);
		return(1);
	}
	
	int p_xy = sqrt(p);	//< Divide no. procs evenly over two axis //  
	int row, column;	
	if(myid+1 <= p_xy){	
		row = 1;
		column = ((myid+1)-(p_xy*(row-1)));	//< calculates column myid is in //
	
	}
	else{
		row = floor((myid+1)/p_xy);		//< calculates the row myid is in //
		if((myid+1)%p_xy !=0){
			row=row+1;
		}
		column = ((myid+1)-(p_xy*(row-1)));
	}

	if(myid+1>p){					//< Error checking myid is within range//
		printf("myid %d is greater than the number of processors %d \n", myid,p);
		return(1);				//< Error //	
	}


	//_________ X Dimension __________//

	int whole_x = n/p_xy;				//< determines number of even divisions //
	int remainder_x = n-whole_x*p_xy; 		//< calculates excess n vector positions to add to n processes //

	if(remainder_x == 0){				//< array is perfectly devisable by processes//
		*s_x = ((column-1)*whole_x) + 1;
		*e_x = *s_x + whole_x-1;
	}
	else if(n<p_xy){				//< More processors than array positions//
		printf("%d more processors than can be utilised for x dimension array \n",(p_xy-n));
		return(1);				//< Error//
	}
	else{						//< Array divides unevenly into processes //
	
		if(column>remainder_x){			//< Myid process above processes with additional array positions 
			*s_x = (remainder_x)*(whole_x+1)+1+((column-1)-(remainder_x))*whole_x;	//< Processes with extra vector value + processes without //	
			*e_x = *s_x+whole_x-1;
		}		
		else{					//< Myid sits within processes which get additional vector prosition for balancing//
			*s_x = ((column-1)*(whole_x+1))+1;
			*e_x = *s_x + (whole_x+1)-1; 	//< (whole+1) counts additional vector position //
		}
	}
	
	
	
	//_________ Y Dimension _____________//
	
	int whole_y = n/p_xy; 				//< determines number of even divisions //
	int remainder_y = n-whole_y*p_xy; 		//< calculates excess n vector positions to add to n processes //

	if(remainder_y == 0){				//< array is perfectly devisable by processes
		*s_y = ((row-1)*whole_y)+ 1;
		*e_y = *s_y + whole_y -1;
	}
	else if(n<p_xy){				//< More processors than array positions//
		printf("%d more processors than can be utilised for y dimension array \n",(p_xy-n));
		return(1);				//< Error//
	}
	else{						//< Array divides unevenly into processes //
	
		if(row>remainder_y){			//< Myid process above processes with additional array positions 
			*s_y = (remainder_y)*(whole_y+1)+1+((row-1)-(remainder_y))*whole_y;	//< Processes with extra vector value + processes without //	
			*e_y = *s_y+whole_y-1;
		}		
		else{					//< Myid sits within processes which get additional vector prosition for balancing//
			*s_y = ((row-1)*(whole_y+1))+1;
			*e_y = *s_y + (whole_y+1)-1; 	//< (whole+1) counts additional vector position //
		}
	}

	return(0);					//< Success //
	
}


/**
 * \brief	Function takes a grid size MxN and number of
 *		processes No_procs to evenly divide and distribute
 *		the matrix between processos. This acheives a load 
 *		ballanced decomposition of an input matrix and 
 *		returns respective start s and end e array position
 *		values for each respective dimension & process myid.
 * 
 * \notes	- This does not require a square grid or an even
 * 		  distribution of processors in x and y dimension
 *		- Assumes myid values are between 0 - (No_proc-1)
 *		- Returns 0 on success or 1 on failure
 *		- Returns failure if more processors than vector
 *		  positions as this would not be balanced.
 *		- Output incompatability s=1 => array[0] 
 * \illustration
 *		- The distribution of a 2D x,y matrix with
 *		  a gridsize=5, across 9 processors
 *
 * 		 ____________________    proc:
 * 	row 1	|__1__2_|__3__4_|__5_|__  3
 * 	row 2 	|  6  7 |  8  9 | 10 |
 * 	row 3	|_11_12_|_13_14_|_15_|__  6
 *	row 4	| 16 17 | 18 19 | 20 |
 *	row 5	|_21_22_|_23_24_|_25_|__  9
 *		|	|	|    |
 *	column:	   1  2    3  4    5
 *	proc:	    7       8      9  
 *
 *		This code should find for myid=4 is designated centre of grid
 *		and s_x=3, e_x=4, s_y=3, e_y=4 
 *
 */
int MPI_Rec_Grid_Decomp2d(const int No_procs, const int No_procs_X, const int No_procs_Y, const int X_dim, const int Y_dim, const int myid, int *s_x, int *e_x, int *s_y, int *e_y){

	// ______ Designate Proc Numbers to Grids _____//
	int row, column;	
	if(myid+1 <= No_procs_X){
		row = 1;
		column = ((myid+1)-(No_procs_X*(row-1)));	//< calculates column myid is in
	}
	else{				
		row = floor((myid+1)/No_procs_X);	//< calculates the row myid is in
		if((myid+1)%No_procs_X != 0){		//< if divides evenly	
			row=row+1;
		}
		column = ((myid+1)-(No_procs_X*(row-1)));	 
	}
		if(myid+1 > No_procs){			//< Error checking myid is within range//	
		printf("myid %d is greater than the number of processors %d \n", myid, No_procs);
		return(1);				//< Error //	
	}

	//_________ X Dimension __________//
	int whole_x = X_dim/No_procs_X; 		//< determines number of even divisions //
	int remainder_x = X_dim-whole_x*No_procs_X;	//< calculates excess n vector positions to add to n processes //

	if(remainder_x == 0){				//< Array is perfectly devisable by processes
		*s_x = ((column-1)*whole_x) + 1;
		*e_x = *s_x + whole_x-1;
	}
	else if(X_dim<No_procs_X){			//< More processors than array positions//
		printf("%d more processors than can be utilised for x dimension array \n",(No_procs_X-X_dim));
		return(1);				//< Error//
	}
	else{						//< Array divides unevenly into processes //
	
		if(column>remainder_x){			//< Myid process above processes with additional array positions 
			*s_x = (remainder_x)*(whole_x+1)+1+((column-1)-(remainder_x))*whole_x;	//< Processes with extra vector value + processes without //	
			*e_x = *s_x+whole_x-1;
		}		
		else{					//< Myid sits within processes which get additional vector prosition for balancing//
			*s_x = ((column-1)*(whole_x+1))+1;
			*e_x = *s_x + (whole_x+1)-1;	//< (whole+1) counts additional vector position //
		}
	}
		
	//_________ Y Dimension _____________//
	int whole_y = Y_dim/No_procs_Y;			//< determines number of even divisions //
	int remainder_y = Y_dim-whole_y*No_procs_Y; 	//< calculates excess n vector positions to add to n processes //

	if(remainder_y == 0){				//< array is perfectly devisable by processes //
		*s_y = ((row-1)*whole_y)+ 1;
		*e_y = *s_y + whole_y -1;
	}
	else if(Y_dim<No_procs_Y){			//< More processors than array positions//
		printf("%d more processors than can be utilised for y dimension array \n",(No_procs_Y-Y_dim));
		return(1);				//< Error//
	}
	else{						//< Array divides unevenly into processes //
	
		if(row>remainder_y){			//< Myid process above processes with additional array positions 
			*s_y = (remainder_y)*(whole_y+1)+1+((row-1)-(remainder_y))*whole_y;	//< Processes with extra vector value + processes without //	
			*e_y = *s_y+whole_y-1;
		}		
		else{					//< Myid sits within processes which get additional vector prosition for balancing//
			*s_y = ((row-1)*(whole_y+1))+1;
			*e_y = *s_y + (whole_y+1)-1; 	//< (whole+1) counts additional vector position //
		}
	}

	return(0); 					//< Success //
}


