SHELL = /bin/sh

# GNU Fotran setting
MPIFC=mpiifort
FC=ifort
#LFLAGS=  -g -fpe0 -fpe-all=0 -check -debug -traceback
#CFLAGS= -c -g -fpe-all=0 -fpe0 -check -debug -traceback
LFLAGS=  -O3 -fopenmp
CFLAGS= -c -O3 -fopenmp

# Intel Fortran setting


#Object files
OBJS   = velocityGrid.o physicalGrid.o mpiParams.o flow.o solver.o main.o

#Link into an excutable
dvm.x: $(OBJS); $(MPIFC) $(LFLAGS) $(OBJS) -o dvm.x

#Compile the modules
main.o           : main.F90 flow.F90 mpiParams.F90 solver.F90; $(MPIFC) $(CFLAGS) main.F90
velocityGrid.o   : velocityGrid.F90; $(MPIFC) $(CFLAGS) velocityGrid.F90
physicalGrid.o   : physicalGrid.F90; $(MPIFC) $(CFLAGS) physicalGrid.F90
flow.o           : flow.F90 velocityGrid.F90 physicalGrid.F90 mpiParams.F90; $(MPIFC) $(CFLAGS) flow.F90
mpiParams.o      : mpiParams.F90 physicalGrid.F90; $(MPIFC) $(CFLAGS) mpiParams.F90
solver.o         : solver.F90 velocityGrid.F90 physicalGrid.F90 flow.F90 mpiParams.F90; $(MPIFC) $(CFLAGS) solver.F90

# Build options
all: dvm.x

clean:
	rm -f ./*.o ./*.mod ./dvm.x
