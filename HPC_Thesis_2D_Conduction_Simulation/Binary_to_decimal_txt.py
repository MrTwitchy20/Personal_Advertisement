#/**
# * \file       Binary_to_decimal_txt.py
# * \brief      Code converts the output of an MPI I/O binary file to decimal string representation in a .txt
# *		The binary file is a non-avoidable format from MPI as is the only option available 
# *             where it is desired to have the output simulation results written in parallel by multiple cores.
# *		The binary format is not easily human readable particularly where double precision numbers
# *		must be read in hence the need for a readable/ manipulatable format.
# *		This file can be made to suit (int or double) precison.
# * \author     F.OSuibhne
# * \date       05-07-21
# * \version    2.0v
# * \note       - The binary file can be read/converted to Little-Endian binary from command line using: 
# *		  "xxd -a -b filename > filename.txt" //xxd is a software for reading binary files.
# *               Could also use "chmod +x Binary_to_decimal_txt.py" followed by ./Binary_to_decimal_txt.py m n filename
# *             - python has a cvs module available which could also visulaise with excell if needed as opposed to .txt output
# *		- The binary on chuck is stored in Little-Endian format and the below code is written with this in mind,
# *		  this may potentially pose an issue/error where run with procs utilizing different formats.
# *		- The function will overwrite an existing file if one exists with the same output filename.
# *		- Format types for representing the double prescision number = 1.00:
# *		-----------------------------------------------------------------------------------------------------------
# *		| Padded Hexidecimal -> | 			     0000 0000 0000 f03f				  |
# *		| Hexbyte -> 		| 			\x00\x00\x00\x00\x00\x00\xf0? 				  |
# *		| Little-Endian binary->| 	0000000000000000000000000000000000000000000000001111000000111111          |
# *		| Big-Endian ->		| 	0011111111110000000000000000000000000000000000000000000000000000 	  |
# *		| Big-Endian binary ->	| 00111111 11110000 00000000 00000000 00000000 00000000 00000000 00000000         | 
# *             | 			| sign |  Exponent   |                         Mantissa                           |
# *		| IEE 64-bit breakdown  |   0  | 01111111111 | 0000 00000000 00000000 00000000 00000000 00000000 00000000 |
# *		----------------------------------------------------------------------------------------------------------- 
# * \usage	Obtain raw binary file output from a MPI_File_write(); then simply run "python3 Binary_to_decimal_txt.py m n filename"
# *		where m = rows, n = columns and filename is user specified. this will cause a filename_binary.txt to be produced and saved
# *		in the directory it is run from. This should contain the grid values seperated by spaces and newlines.
# **/

import sys  # Required to take command line arguments

#print(sys.argv,"\n",len(sys.argv)) #< echos user input back to commandline
if len(sys.argv) < 4:
	print("Error! Executable called with fewer arguments than expected \nCommand should be of the form: 'python3 Binary_to_decimal_txt.py m n filename'")
if len(sys.argv) > 4:
	print("Error! Executable called with more arguments than expected \nCommand should be of the form: 'python3 Binary_to_decimal_txt.py m n filename'")

m = int(sys.argv[1])   #< rows of data (nodes)
n = int(sys.argv[2])   #< columns of data (nodes)
# Num_Nodes = m*n #< qty of data to be read
# print("Num_Nodes:", Num_Nodes)

file_in = str(sys.argv[3])  # input filename
file_out = str(sys.argv[3]) + str("_decimal.txt") # generates own output filename
#print("input_file:", file_in) 
print("output_file:", file_out)

import struct   # contains unpack function

file_out_ptr = open(file_out,"w")   #< open to write to file
file_in_ptr = open(file_in,'rb')    #< open to read in bits from raw binary file (not .txt)

with file_in_ptr:
    for i in range (0,m):           #< loop over rows
        for j in range (0,n):       #< loop over columns
            byte = file_in_ptr.read(8)  #< Read in 8 bits at a time (double)
            #byte = file_in_ptr.read(4)   #< Read in 4 bits at a time (int)
            #print(byte) #< prints raw binary
            #print(format(struct.unpack('d',byte)[0],'.8f'),file=file_out_ptr)  #< format is used to print precision, print used to output number as opposed to strings to file
                                                                                # write only writes strings, we could manipulate bype in this way first and still write string if we desired
                                                                                # print(x,file=filename)    prints to the file instead of command line
            byte=str(format(struct.unpack('d',byte)[0],'.5f'))  #< Converts binary->Decimal String. unpacking 8bits at a time, creates a turple of one elememt "[0]" 
								# sets byte to string of that element to a set precision level ".2f"
								# NOTE: Could pose issues when reading back from file to form a heatmap)
                                                                # unpacking 8 bits at a time, creates a turple of one element "[0]" sets byte to string of element
            # byte = str(format(struct.unpack('i',byte)[0]))    #< As previous but for but int            
            # print(byte) #< prints decimal converted
            file_out_ptr.write(byte)				#< writes to file one "int" or "double" at a time
            file_out_ptr.write(" ")
        file_out_ptr.write("\n")

file_out_ptr.close() # close files
file_in_ptr.close()

