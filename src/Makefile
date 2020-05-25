# Intel Fotran setting
#MPIFC=mpiifort
#FC=ifort

# GNU Fotran setting
MPIFC=mpif90 -f90=gfortran
FC=gfortran

# Release mode flags
#LFLAGS=  -g -fpe0 -fpe-all=0 -check -debug -traceback
#CFLAGS= -c -g -fpe-all=0 -fpe0 -check -debug -traceback

#LFLAGS=   -fopenmp
#CFLAGS= -c   -fopenmp

# Release mode flags
LFLAGS= -fopenmp
CFLAGS= -c -fopenmp

# Intel Fortran setting

#Object files
OBJS   = velocityGrid.o physicalGrid.o mpiParams.o flow.o \
		 solver.o main.o gaussHermite.o parameters.o \
		 string_utility_module.o

#Link into an excutable
mpi_dvm2d: $(OBJS)
	$(MPIFC) $(LFLAGS) $(OBJS) -o mpi_dvm2d

#Compile the modules
main.o : main.F90 flow.o mpiParams.o solver.o parameters.o
	$(MPIFC) $(CFLAGS) main.F90
parameters.o : parameters.F90 velocityGrid.o physicalGrid.o \
	mpiParams.o solver.o flow.o string_utility_module.o 
	$(MPIFC) $(CFLAGS) parameters.F90
string_utility_module.o : string_utility_module.F90
	$(MPIFC) $(CFLAGS) string_utility_module.F90
velocityGrid.o : velocityGrid.F90 gaussHermite.o
	$(MPIFC) $(CFLAGS) velocityGrid.F90
physicalGrid.o : physicalGrid.F90
	$(MPIFC) $(CFLAGS) physicalGrid.F90
flow.o : flow.F90 velocityGrid.o physicalGrid.o mpiParams.o
	$(MPIFC) $(CFLAGS) flow.F90
mpiParams.o : mpiParams.F90 physicalGrid.o
	$(MPIFC) $(CFLAGS) mpiParams.F90
solver.o : solver.F90 velocityGrid.o physicalGrid.o flow.o mpiParams.o
	$(MPIFC) $(CFLAGS) solver.F90
gaussHermite.o : gaussHermite.F90
	$(MPIFC) $(CFLAGS) gaussHermite.F90

# Build options
all: mpi_dvm2d

clean:
	rm -f ./*.o ./*.mod ./mpi_dvm2d
