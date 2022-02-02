CC = gcc

MPICC = mpicc

CFLAGS = -g -Wall

OPTIMIZE = -xSSE4.2 -O3 

LDFLAGS = -lm

EXECS = 2D_Annulus_Transient_Fixed_Boundary 2D_Annulus_Fixed_Boundary 2D_Annulus_Convective 2D_Rect_Fixed_Boundary_SS decomp2test Parallel_write_to_file_test 

all: $(EXECS)

decomp2test: decomp2test.o decomp2d.o 
	$(MPICC) decomp2d.c $(CFLAGS) -o $@ $< $(LDFLAGS)

Parallel_write_to_file_test: Parallel_write_to_file_test.o Parallel_write_to_file.o decomp2d.o
	$(MPICC) Parallel_write_to_file.c decomp2d.c $(CFLAGS) -o $@ $< $(LDFLAGS)

2D_Rect_Fixed_Boundary_SS: 2D_Rect_Fixed_Boundary_SS.o Parallel_write_to_file.o decomp2d.o memory_alloc.o jacobi_Itteration.o
	$(MPICC) Parallel_write_to_file.c decomp2d.c memory_alloc.c jacobi_Itteration.c $(CFLAGS) -o $@ $< $(LDFLAGS) $(OPTIMIZE)

2D_Annulus_Fixed_Boundary: 2D_Annulus_Fixed_Boundary.o Parallel_write_to_file.o decomp2d.o memory_alloc.o jacobi_Itteration.o
	$(MPICC) Parallel_write_to_file.c decomp2d.c memory_alloc.c jacobi_Itteration.c $(CFLAGS) -o $@ $< $(LDFLAGS) $(OPTIMIZE)

2D_Annulus_Convective: 2D_Annulus_Convective.o Parallel_write_to_file.o decomp2d.o memory_alloc.o jacobi_Itteration.o
	$(MPICC) Parallel_write_to_file.c decomp2d.c memory_alloc.c jacobi_Itteration.c $(CFLAGS) -o $@ $< $(LDFLAGS)

2D_Annulus_Transient_Fixed_Boundary: 2D_Annulus_Transient_Fixed_Boundary.o Parallel_write_to_file.o decomp2d.o memory_alloc.o jacobi_Itteration.o
	$(MPICC) Parallel_write_to_file.c decomp2d.c memory_alloc.c jacobi_Itteration.c $(CFLAGS) -o $@ $< $(LDFLAGS)


.PHONY: clean

clean:
	$(RM) *.o $(EXECS)


