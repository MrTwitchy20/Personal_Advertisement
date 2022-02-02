/**
 * \filename 	Parallel_write_to_file.c
 * \brief	Contains functions for printing int/double Grids to a file in Parallel using fileviews and MPI IO functions
 * \author	F.OSuibhne
 * \date	27.07.21 
**/

#include<mpi.h>
#include<stdio.h>

/*
 * \brief	Function utilises MPI IO commands to allow each node/proc to print
 *		their respective section of the grid to a binary file in parallel as integers.  	
 * \notes	- Function requires use of openmpi commands
 *		- This is unable to ignore ghost columns/rows
 *		- MPI_COMM_WORLD is a global communicator and as such cant be passed in with 
 *		  "Function(MPI_Comm MPI_COMM_WORLD)".
 *		- This could present an issue if it is desired to specify individual communicators 
 *		  to print. A solution could be force a new MPI_COMM == MPI_COMM_WORLD
 *		Benefits: 
 *		- Saving to a single file by all procs similtaniously avoiding alternative solutions which
 *		  would essentially serialize the process wasting significant resource time as nodes sit idle
 *		  waiting for their turn to print or communicate to a master node for printing.
 * 		- This is a more scalable solution, it would allow for larger grid sizes. Nodes can maxamise
 *		  memory usage without the need to reserve memory for communication buffers.
 *		- This is particularly nessiary for a Transient usage case were multiple prints may be needed throughout
 *		  the solution 
 *		Disadvantage:
 *		- The output file resulting from this function is not immediatly user readable and requires a further
 *		  process to interpret/reformat the data.
 *
 * \ToDo	- MPI_File_set_view() using "natve", Might this cause issues if run on disimilar devices?
 */
void Parallel_Print_Int_To_File(const char *fname, const int Num_x, const int Num_y, const int MY_BUFSIZE, const int buf[], const int grid_size, const int myid, const int s_x, const int e_x, const int s_y, const int e_y){

  	MPI_Aint lb, extent;
  	MPI_Datatype etype, filetype, contig;
  	MPI_Offset disp;
  	MPI_File fh;
 
 	MPI_File_open(MPI_COMM_WORLD, fname, MPI_MODE_CREATE | MPI_MODE_RDWR,
  	MPI_INFO_NULL, &fh);	//< All procs open file and have a filepointer //
 
 	MPI_Type_contiguous(Num_x, MPI_INT, &contig); //< Proc dependant datatype is created //
  
  	extent = grid_size*sizeof(int);       //< contigous datatypes are spaced apart by grid_size //
  	lb = 0;       			//< "New lower bound of datatype (address integer)"
  	MPI_Type_create_resized(contig, lb, extent, &filetype);  //< contiguous types, no offset, of 6 doubles ?
  	MPI_Type_commit(&filetype);

  	disp = (s_y-1)*extent + (s_x-1)*sizeof(int); //<  Changed 19.07.21 and seems to work, 2 should be number in MPI_TYPE_contiguous above !!
  	//disp=0; //< need somthing to consider rows
 
 	etype = MPI_INT;
  	MPI_File_set_view(fh, disp, etype, filetype, "native", MPI_INFO_NULL);	//< SHOULD WE BE STATING NATIVE OR LOWER ELISION?
  	MPI_File_write(fh, buf, MY_BUFSIZE, MPI_INT, MPI_STATUS_IGNORE);	//< writes proc specific buffer to file //
  	MPI_File_close(&fh);

	MPI_Type_free(&filetype);	//< Must free datatype to prevent memory leak
}



/*
 * \brief	Function utilises MPI IO commands to allow each node/proc to print
 *		their respective section of the grid to a binary file in parallel as doubles.
 *		 
 * \notes	- Function requires use of openmpi commands
 * 		- This is unable to ignore ghost columns/rows
 *		- MPI_COMM_WORLD is a global communicator and as such cant be passed in with 
 *		  "Function(MPI_Comm MPI_COMM_WORLD)".
 *		- This could present an issue if it is desired to specify individual communicators 
 *		  to print. A solution could be force a new MPI_COMM == MPI_COMM_WORLD
 *		Benefits: 
 *		- Saving to a single file by all procs similtaniously avoiding alternative solutions which
 *		  would essentially serialize the process wasting significant resource time as nodes sit idle
 *		  waiting for their turn to print or communicate to a master node for printing.
 * 		- This is a more scalable solution, it would allow for larger grid sizes. Nodes can maxamise
 *		  memory usage without the need to reserve memory for communication buffers.
 *		- This is particularly nessiary for a Transient usage case were multiple prints may be needed throughout
 *		  the solution 
 *		Disadvantage:
 *		- The output file resulting from this function is not immediatly user readable and requires a further
 *		  process to interpret/reformat the data.
 *
 * \ToDo	- MPI_File_set_view() using "natve", Might this cause issues if run on disimilar devices?
 */
void Parallel_Print_Double_To_File(const char *fname, const int Num_x, const int Num_y, const int MY_BUFSIZE, const double buf[], const int grid_size, const int myid, const int s_x, const int e_x, const int s_y, const int e_y){

	MPI_Aint lb, extent; 
  	MPI_Datatype etype, filetype, contig;
  	MPI_Offset disp;
  	MPI_File fh;	//< file handle, can be imagined as a pointer or curser within a file //
  	
	MPI_File_open(MPI_COMM_WORLD, fname, MPI_MODE_CREATE | MPI_MODE_RDWR,
  	MPI_INFO_NULL, &fh);	//< All procs open file and have a filepointer //
				// MPI_MODE_WRONLY-- write only //
 
 	MPI_Type_contiguous(Num_x, MPI_DOUBLE, &contig); //< Proc dependant datatype is created//
  
 	extent = grid_size*sizeof(double);      //< contigous datatypes are spaced apart by grid_size //
  	lb = 0;       				//< Lower bound of datatype (address integer)"
  	MPI_Type_create_resized(contig, lb, extent, &filetype);  //< contiguous types, no offset, of 6 doubles ?
  	MPI_Type_commit(&filetype);

 	disp = (s_y-1)*extent + (s_x-1)*sizeof(double); //< displacement your offset from 0 + row displacement
 
  	etype = MPI_DOUBLE;
  	MPI_File_set_view(fh, disp, etype, filetype, "native", MPI_INFO_NULL);	//< Ensures each proc only see's file areas relevent to it//	
  	MPI_File_write(fh, buf, MY_BUFSIZE, MPI_DOUBLE, MPI_STATUS_IGNORE);	//< writes proc specific buffer to file //
  	MPI_File_close(&fh);

	MPI_Type_free(&filetype);	//< Must free datatype to prevent memory leak
}

/*
 * \brief	Function utilises MPI IO commands to allow each node/proc to print
 *		their respective section of the grid to a binary file in parallel as doubles.
 *		This is the file which will be used primarily as it allows printin of boundary conditions
 *		and ignores all the ghost columns/rows.
 *				 
 * \notes	- Function requires use of openmpi commands
 * 		- This is unable to ignore ghost columns/rows
 *		- MPI_COMM_WORLD is a global communicator and as such cant be passed in with 
 *		  "Function(MPI_Comm MPI_COMM_WORLD)".
 *		- This could present an issue if it is desired to specify individual communicators 
 *		  to print. A solution could be force a new MPI_COMM == MPI_COMM_WORLD
 *		Benefits: 
 *		- Saving to a single file by all procs similtaniously avoiding alternative solutions which
 *		  would essentially serialize the process wasting significant resource time as nodes sit idle
 *		  waiting for their turn to print or communicate to a master node for printing.
 * 		- This is a more scalable solution, it would allow for larger grid sizes. Nodes can maxamise
 *		  memory usage without the need to reserve memory for communication buffers.
 *		- This is particularly nessiary for a Transient usage case were multiple prints may be needed throughout
 *		  the solution 
 *		Disadvantage:
 *		- The output file resulting from this function is not immediatly user readable and requires a further
 *		  process to interpret/reformat the data.
 * \usage	When calling function with ghost columns input buff[] should be e.g. "(&Grid[1])+1", this skips the first row
 * \ToDo	- MPI_File_set_view() using "natve", Might this cause issues if run on disimilar devices?
 */
void Parallel_Print_Double_To_File_wBoundary(const char *fname, const int Num_x, const int Num_y, const int MY_BUFSIZE, const double buf[], const int grid_size, const int myid, const int s_x, const int e_x, const int s_y, const int e_y){

	MPI_Aint lb, extent; 
  	MPI_Datatype etype, filetype, contig;
  	MPI_Offset disp;
  	MPI_File fh;	//< File handle, can be imagined as a curser or poointer within a file	
 	int Error_check;
		
  	Error_check = MPI_File_open(MPI_COMM_WORLD, fname, MPI_MODE_CREATE | MPI_MODE_RDWR,
  	MPI_INFO_NULL, &fh);	//< All procs open file and have a filepointer //
				// MPI_MODE_WRONLY //
				// MPI_MODE_APPEND //
	if(Error_check != 0){
		char errormsg[150];
		sprintf(errormsg,"Failure occured! MPI_File_open(); returned an error (%d) with proc= %d at File: %s, Line: %d\n",Error_check, myid, __FILE__, __LINE__);
		perror(errormsg);
		// Doesnt envoke a fatal error, maybe some procs are successfull and will just be missing chunk//
	}	

  	MPI_Type_contiguous(Num_x, MPI_DOUBLE, &contig);	//< Proc dependant datatype is created // 
  	extent = grid_size*sizeof(double);    			//< contigous datatypes spaced apart by grid_size //
  	lb = 0;       						//< lower bound of datatype // 
  	MPI_Type_create_resized(contig, lb, extent, &filetype); //< contiguous types, no offset
  	MPI_Type_commit(&filetype);

 	MPI_Aint g_lb, g_extent;
  	MPI_Datatype contig_ghost, from_proc; 
  
  	MPI_Type_contiguous(Num_x, MPI_DOUBLE, &contig_ghost);			//< Required to send from buffer avoiding ghost columns
  	g_lb = 1*sizeof(double); 						//< Offset to skip first ghost column 
  	g_extent = (Num_x+2)*sizeof(double);					//< Extent of local processor grid 
  	MPI_Type_create_resized(contig_ghost, g_lb, g_extent, &from_proc);	//< Resized to ignore ghost columns
  	MPI_Type_commit(&from_proc);
		
  
 	disp = (s_y-1)*extent + (s_x-1)*sizeof(double); 			//< file pointer offset processor specific //
  	etype = MPI_DOUBLE;
  	MPI_File_set_view(fh, disp, etype, filetype, "native", MPI_INFO_NULL);	//< Hiding all non-proc specifid data slots //
  	MPI_File_write(fh, buf, Num_y, from_proc, MPI_STATUS_IGNORE);		//< Writes proc specific buffer to file //
  	
	Error_check = MPI_File_close(&fh);
  	if(Error_check != 0){
		char errormsg[150];
		sprintf(errormsg,"Failure occured! MPI_File_close(); returned an error (%d) with proc= %d at File: %s, Line: %d\n",Error_check, myid, __FILE__, __LINE__);
		perror(errormsg);
		// Doesnt envoke a fatal error //
  	}
	
	MPI_Type_free(&filetype);	//< Must free committed datatypes to prevent memory leaks
	MPI_Type_free(&from_proc);	
}

