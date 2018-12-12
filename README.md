### Short cuts (alias)

Put add this line on your ~/.bashrc

```bash
alias dvm='export PRJ=${HOME}/mpi_dvm; source ${PRJ}/conf/bashrc cylinders'
```

Any time your login and needs to work on the project, type `dvm`. 
Then the aliases and related environments are set for you. 
For the compelete list of environment and aliases see `conf/bashrc`

### How to compile for the first time and the re-compilation

Supposing your are current in the project directory (the one with src and config direcories in it),
Copy the reference configure file for cmake to the current direcotry and name it as `CMakeListst.txt`.

```bash
cp conf/CMakeLists.txt.example CMakeListst.txt
```

#### For the Release mode compilation

```bash
mkdir build
cd build
cmake ..
make -j 6
```

This will do the Release mode compilation. 
Anytime you make some modification on the source code, and you want to re-compile, just cd to the `build` directory, and make:

```
cd $BUILD
make -j 6
```

or simplely,

```bash
make -C $BUILD
```

or even more simply using the alias we have set,

```bash
mb
```
which means make for build. The actural command mb issues can be seen by `alias mb`. 

#### For the Debug mode compilation
For Debug build, similarly, start from the project directory

```bash
mkdir debug
cmake .. -DCMAKE_BUILD_TYPE=Debug
make
```
Anytime you make some modifications on the source code, and you want to re-compile for the Debug mode, just cd to the `debug` directory,
```
cd $DEBUG
make -j 6
```

or

```bash
make -C $DEBUG
```

or using the alias

```bash
md
```

### How to run it

Copy the demo directory to your working directory, and modify the `para.in` file in it, and then replace the image file.

Run the solver using the folowing way

```bash
OMP_NUM_THREADS=12 mpirun -np 2 $BUILD/dvm.x

```
which will use 2 MPI processes and 12 OpenMP threads for each MPI process.

All parameters are set throught the `para.in`.  Below is an example of `para.in` file

```
&physicalNml
imageFileName = 'example.dat'
Nx = 500,
Ny = 500,
xPadRatio = 0,
wallExtOrder = 1/

&velocityNml
Nc_fundamental = 4,
halfRange = T /

&mpiNml
mpi_xdim = 4,
mpi_ydim = 4/

&solverNml
maxStep = 10000000,
chkConvergeStep = 100,
saveStep = 5000000,
saveFormat = 1,
eps = 1.0d-8/

&flowNml
allKnStr = '1.0d-1,0.5d-2,0.2d-2'
pressDrop = 1.0d-3,
accom = 1.d0/
```

where
* line starting with & means a parameter group, do not change
* `imageFileName` is the file name of the geometry flag (0 for fluid node, and 1 for solid node) file in ascii format. 
   Node flags are seprated with a space in each line
* `Nx` and `Ny` is the grid size in physical space, should be consisent to the grid flag file.
* `xPadRatio` is the total relative fulid padding length in the inlet and outlet. Default value is zero and can be omit if no fluid layers is added.
* `wallExtOrder` is  the accuracy order of extropolation at the wall nodes. can be [1|2|3]. 1 means first order; 2 means second order; 3 means use first order when the wall is one fluid node off the communication bounary, and use second order otherwise. For stability, recommand to use 3.
* `Nc_fundamental` is number of discrete point in half of an axis.  Here 4 means a total of $(2x4)^3$ for 3D problmes. Valida values are  2, 4, 6, 8, 12, 16.
* `halfRange` can be [T|F]; Use half-range Gauss-Hermite or not, eighr T or F.
* `mpi_ydim` and `mpi_xdim` : Use a Cartesian even domain decomposition, these are the numbers of subdomains in X and Y direction
* `maxSteps` is the maximum of the number of iterations.
* `chkConvergeStep` the frequency to check the convergence of the permeability
* `saveStep` the frequency to save the flow fields.
* `saveFormat` the format to save data, can be [1|2|3]. 1 means in .pvti format, 2 means in tecplot format, 3 means in vtk format.
* `eps` is the convergence crieterion.
* `allKnStr` lists all the Knudsen nubmers want to be simulated.
* `pressDrop` is the non-dimensional pressure drop between the inelt and outlet, no need to change.
* `accom` is the accomodation coefficient on the solid wall.

**Note**: Putting higer Kn numbers ahead in allKnStr is better than the reverse, as for high Kn number, the solver converges fast, and the converged flow fields are used as the initialization value of the next lower Kn case. If use multiple Kn nubmers, the sovler will create a output directory for each of the Kn. A `Results.dat` will be generated and the screent output will be logged, which contains the permeability values for each of the cases, and their convergence history.


