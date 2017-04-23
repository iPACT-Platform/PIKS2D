SHELL = /bin/sh

# GNU Fotran setting
#MPIFC=mpiifort
#FC=ifort
MPIFC=ftn
FC=ftn
#LFLAGS=  -g -fpe0 -fpe-all=0 -check -debug -traceback
#CFLAGS= -c -g -fpe-all=0 -fpe0 -check -debug -traceback
#LFLAGS=   -fopenmp
#CFLAGS= -c   -fopenmp
LFLAGS= -g
CFLAGS= -c -g

# Intel Fortran setting


#Object files
OBJS   = velocityGrid.o physicalGrid.o mpiParams.o flow.o solver.o main.o gaussHermite.o parameters.o

#Link into an excutable
dvm.x: $(OBJS); $(MPIFC) $(LFLAGS) $(OBJS) -o dvm.x

#Compile the modules
main.o           : main.F90 flow.o mpiParams.o solver.o parameters.o; $(MPIFC) $(CFLAGS) main.F90
parameters.o     : velocityGrid.o physicalGrid.o mpiParams.o solver.o flow.o parameters.F90; $(MPIFC) $(CFLAGS) parameters.F90
velocityGrid.o   : gaussHermite.o velocityGrid.F90; $(MPIFC) $(CFLAGS) velocityGrid.F90
physicalGrid.o   : physicalGrid.F90; $(MPIFC) $(CFLAGS) physicalGrid.F90
flow.o           : flow.F90 velocityGrid.o physicalGrid.o mpiParams.o; $(MPIFC) $(CFLAGS) flow.F90
mpiParams.o      : mpiParams.F90 physicalGrid.o; $(MPIFC) $(CFLAGS) mpiParams.F90
solver.o         : solver.F90 velocityGrid.o physicalGrid.o flow.o mpiParams.o; $(MPIFC) $(CFLAGS) solver.F90
gaussHermite.o   : gaussHermite.F90; $(MPIFC) $(CFLAGS) gaussHermite.F90

# Build options
all: dvm.x

clean:
	rm -f ./*.o ./*.mod ./dvm.x
