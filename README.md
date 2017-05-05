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

```bash
mkdir -p run/yourCase
cd run/yourCase
cp $CONF/para.in .
# do some adjusting in para.in
# then prepare the iamge data file in the current case direcory
# run the release mode exe.:
mpirun -np 2 $BUILD/dvm.x 
# or run the debug mode exe.:
mpirun -np 2 $DEBUG/dmv.x
```
