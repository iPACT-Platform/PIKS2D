

### How to compile 

Suppose that your are currently in the project directory (the one with `src` and `config` directories in it) and `openmpi` is available

```bash
cd src
make
```
This will create executive file `piks2d` inside the `src`  directory.
### How to run simulation

Copy the `tutorials` folder, which consists of a few subfolder for some example cases, to your working directory. In each subfolder of `tutorials` folder, there are an image file and a `para.in`  file. To adapt your simulation setup, you can change the values of parameters in the `para.in` file and replace the image file.

From the directory of tutorial case, run the solver using the following command

```bash
OMP_NUM_THREADS=8 mpirun -np 2 ../../src/piks2d

```
which will use, for example, 2 MPI processes and 8 OpenMP threads for each MPI process.

### How to interpret the simulation parameters
All parameters are set in the `para.in`.  Below is an example of `para.in` file in the directory `tutorials/Sierpinski_carpet/`
```
&physicalNml
imageFileName = 'carpetL3.dat'
Nx = 540,
Ny = 540,
wallExtOrder = 1/

&velocityNml
Nc_fundamental = 4,
halfRange = T /

&mpiNml
mpi_xdim = 2,
mpi_ydim = 1/

&solverNml
maxStep = 100,
chkConvergeStep = 10,
saveStep = 50,
saveFormat = 1,
eps = 1.0d-8/

&flowNml
allKnStr = '1.0d-1,5.0d-2,2.0d-2',
pressDrop = 1.0d-3,
accom = 1.d0/
```
where
* line starting with & means a parameter group, do not change!
* `imageFileName` is the file name of the image file, i.e. geometry flag (0 for fluid node, and 1 for solid node) file in ASCII format. 
   The coordinate origin is set at bottom left corner in rectangular domain of the image. The positive direction in x-, y-direction is to the right and to the top, respectively. Each line and column in the image file corresponds to a fixed y- and x-coordinate, respectively. Nodes are separated by a space " " in x-direction.
* `Nx` and `Ny` is the number of nodes in x- and y-direction in the spatial space, should be consistent to the image file.
* `wallExtOrder` is  the accuracy order of extrapolation at the wall nodes. can be [1|2|3]. 1 means first order; 2 means second order; 3 means use first order when the wall is one fluid node off the communication boundary, and use second order otherwise. 
* `Nc_fundamental` is number of discrete point in half of an axis.  Here 4 means a total of $(2x4)^2$ for 2D problems. Available values are  2, 4, 6, 8, 12, 16.
* `halfRange` can be [T|F]: choose half-range Gauss-Hermite or not, either T or F.
* `mpi_ydim` and `mpi_xdim` : these are the numbers of subdomains in x- and y-direction in a frame of Cartesian even domain decomposition. The number of MPI processes corresponds to the multiple of `mpi_ydim` and `mpi_xdim`.  
* `maxSteps` is the maximum number of iterations.
* `chkConvergeStep` the frequency to check the convergence of the permeability.
* `saveStep` the frequency to save the flow fields.
* `saveFormat` the format to save data, can be [1|2|3]. 1 means in .pvti format, 2 means in Tecplot format, 3 means in vtk format.
* `eps` is value for the convergence criterion. The simulation will stop either the `eps` or the `maxSteps` is achieved, whichever comes first. 
* `allKnStr` lists all the Knudsen numbers of interest for simulations.
* `pressDrop` is the non-dimensional pressure drop between the inlet and outlet. With the current linearized formulation and nondimensionalization, the dimensionless results are independent of positive`pressDrop`.  
* `accom` is the accommodation coefficient in the diffuse-specular model for the gas-solid interaction.

**Note**: 
* If multiple Kn numbers are used in `allKnStr`, the solver will create a output directory for each of the Kn. 
*  Putting higher Kn numbers ahead in `allKnStr` is better than the reverse, as for high Kn number, the solver converges fast, and the converged flow fields are used as the initialization value of the next lower Kn case. 

### How to find the simulation results
All simulation results are stored in the working directory
*  A `Results.dat` file logs the screen output, which contains the permeability values for all of the cases, and their convergence history. Below is and example of `Results.dat` file in the directory `tutorials/Sierpinski_carpet/` when the all simulations finish. From the left to right, columns for Knudsen number, dimensionless mass flow rate, dimensionless permeability, residual and iteration number are represented. 
```
   1.000000E-01   9.433632E-02   5.322357E-03   9.056637E-02     10
   1.000000E-01   5.542642E-02   3.127101E-03   4.124594E-02     20
   1.000000E-01   4.634453E-02   2.614710E-03   1.638548E-02     30
   1.000000E-01   4.257757E-02   2.402182E-03   8.128167E-03     40
   1.000000E-01   4.088889E-02   2.306909E-03   3.966130E-03     50
   1.000000E-01   4.014601E-02   2.264996E-03   1.816818E-03     60
   1.000000E-01   3.981451E-02   2.246293E-03   8.257390E-04     70
   1.000000E-01   3.966996E-02   2.238138E-03   3.630695E-04     80
   1.000000E-01   3.960545E-02   2.234498E-03   1.626111E-04     90
   1.000000E-01   3.957665E-02   2.232873E-03   7.272301E-05    100
   5.000000E-02   4.040217E-02   1.139724E-03   2.085884E-03     10
   5.000000E-02   4.040392E-02   1.139774E-03   4.330736E-06     20
   5.000000E-02   4.043641E-02   1.140690E-03   8.042186E-05     30
   5.000000E-02   4.046309E-02   1.141443E-03   6.597477E-05     40
   5.000000E-02   4.047311E-02   1.141725E-03   2.476907E-05     50
   5.000000E-02   4.048055E-02   1.141935E-03   1.838429E-05     60
   5.000000E-02   4.048407E-02   1.142035E-03   8.694014E-06     70
   5.000000E-02   4.048514E-02   1.142065E-03   2.639421E-06     80
   5.000000E-02   4.048566E-02   1.142079E-03   1.284225E-06     90
   5.000000E-02   4.048530E-02   1.142069E-03   8.908546E-07    100
   2.000000E-02   4.787465E-02   5.402076E-04   1.825194E-02     10
   2.000000E-02   4.781709E-02   5.395581E-04   1.202381E-04     20
   2.000000E-02   4.760965E-02   5.372174E-04   4.338165E-04     30
   2.000000E-02   4.755085E-02   5.365538E-04   1.235168E-04     40
   2.000000E-02   4.751590E-02   5.361595E-04   7.349885E-05     50
   2.000000E-02   4.749186E-02   5.358883E-04   5.058723E-05     60
   2.000000E-02   4.747645E-02   5.357144E-04   3.244868E-05     70
   2.000000E-02   4.746720E-02   5.356100E-04   1.948711E-05     80
   2.000000E-02   4.746133E-02   5.355438E-04   1.236042E-05     90
   2.000000E-02   4.745750E-02   5.355006E-04   8.062526E-06    100
```
* Flow field files, e.g. `Field.pvti`, are stored in   Kn _i_ subfolder for each value of Knudsen number. The _i_ is the value of Knudsen number, e.g. Kn1.0d-1 is the subfolder of the working directory corresponds to the case Knudsen number of 0.1. 
