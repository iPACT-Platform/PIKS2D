 ```
$$$$$$$\ $$$$$$\ $$\   $$\  $$$$$$\   $$$$$$\  $$$$$$$\  
$$  __$$\\_$$  _|$$ | $$  |$$  __$$\ $$  __$$\ $$  __$$\ 
$$ |  $$ | $$ |  $$ |$$  / $$ /  \__|\__/  $$ |$$ |  $$ |
$$$$$$$  | $$ |  $$$$$  /  \$$$$$$\   $$$$$$  |$$ |  $$ |
$$  ____/  $$ |  $$  $$<    \____$$\ $$  ____/ $$ |  $$ |
$$ |       $$ |  $$ |\$$\  $$\   $$ |$$ |      $$ |  $$ |
$$ |     $$$$$$\ $$ | \$$\ \$$$$$$  |$$$$$$$$\ $$$$$$$  |
\__|     \______|\__|  \__| \______/ \________|\_______/                                                   
                                                             
     Parallel Image-based Kinetic Solver (2D)
```

## What is PIKS2D?

PIKS2D is an open-source 2D parallel pore-scale rarefied gas-flow simulator. 
It solves the linearized gas-kinetic equation on a uniform cartesian grid using the Discrete
Velocity Method (DVM). It can be run in parallel both OpenMP and MPI. The porous
structure can be arbitrary complex and and is input into the simulator as an ASCII
based binary images of '0' (fluid) and '1' (solid). An example simulation is shown below using the solver.

<p align="center"><a href="https://ibb.co/v4jqQsy"><img src="https://i.ibb.co/0yG2FB4/dvm-Kn5e-4-U.png" alt="dvm-Kn5e-4-U" border="0" width="480"></a></p>


## How do I use PIKS2D?

See the [wiki](https://github.com/iPACT-Platform/PIKS2D/wiki) pages for installation and tutorials.
See the reference below for the theory and numerical method. The full-text PDFs are provided in the directory `reference`.

## How do I cite PIKS2D?

* Minh Tuan Ho, Lianhua Zhu, Lei Wu, Peng Wang, Zhaoli Guo, Zhi-Hui Li, and Yonghao Zhang. “A Multi-Level Parallel Solver for Rarefied Gas Flows in Porous Media.” Computer Physics Communications 234 (January 1, 2019): 14–25. [DOI: 10.1016/j.cpc.2018.08.009](https://doi.org/10.1016/j.cpc.2018.08.009).
* Minh Tuan Ho, Lianhua Zhu, Lei Wu, Peng Wang, Zhaoli Guo, Jingsheng Ma, and Yonghao Zhang. “Pore-Scale Simulations of Rarefied Gas Flows in Ultra-Tight Porous Media.” Fuel 249 (August 1, 2019): 341–51. [DOI: 10.1016/j.fuel.2019.03.106](https://doi.org/10.1016/j.fuel.2019.03.106).


## License

The PIKS2D is licensed under the MIT license, see the file `LICENSE`.

## Who has funded PIKS2D
Development of PIKS2D received funding from Engineering and Physical Sciences Research Council, European Union’s Horizon 2020 Marie Skłodowska-Curie Individual Fellowship and the “Advanced Hybrid Method for Post-Scale Simulation of Shale Gas Flows”, a global partnership project grant funded by [KFUPM](http://www.kfupm.edu.sa).
