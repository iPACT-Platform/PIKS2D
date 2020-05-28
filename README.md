 ```
 _______  ___   ___   _  _______  _______  ______  
|       ||   | |   | | ||       ||       ||      | 
|    _  ||   | |   |_| ||  _____||____   ||  _    |
|   |_| ||   | |      _|| |_____  ____|  || | |   |
|    ___||   | |     |_ |_____  || ______|| |_|   |
|   |    |   | |    _  | _____| || |_____ |       |
|___|    |___| |___| |_||_______||_______||______| 

     Parallel Image-based Kinetic Solver (2D)
```

## What is PIKS2D?

PIKS is an open-source 2D parallel pore-scale rarefied gas-flow simulator. 
It solves the linearized gas-kinetic equation on a uniform cartesian grid using the Discrete
Velocity Method (DVM). It can be run in parallel both OpenMP and MPI. The porous
structure can be arbitrary complex and and is input into the simulator as an ASCII
based binary images of '0' (fluid) and '1' (solid).

## How do I use PIKS2D?

See the [wiki](https://github.com/iPACT-Platform/PIKS2D/wiki) pages for installation and tutorials.
See the reference below for the theory and numerical method.

## How do I cite PIKS2D?

* Minh Tuan Ho, Lianhua Zhu, Lei Wu, Peng Wang, Zhaoli Guo, Zhi-Hui Li, and Yonghao Zhang. “A Multi-Level Parallel Solver for Rarefied Gas Flows in Porous Media.” Computer Physics Communications 234 (January 1, 2019): 14–25. [DOI: 10.1016/j.cpc.2018.08.009](https://doi.org/10.1016/j.cpc.2018.08.009).
* Minh Tuan Ho, Lianhua Zhu, Lei Wu, Peng Wang, Zhaoli Guo, Jingsheng Ma, and Yonghao Zhang. “Pore-Scale Simulations of Rarefied Gas Flows in Ultra-Tight Porous Media.” Fuel 249 (August 1, 2019): 341–51. [DOI: 10.1016/j.fuel.2019.03.106](https://doi.org/10.1016/j.fuel.2019.03.106).


## LICENSE

The PIKS2D is licensed under the MIT license, see the file `LICENSE`.

## Who is funding PIKS2D
Development of PIKS2D is supported by the Engineering and Physical Sciences Research Council, European Union’s Horizon 2020 Research and Innovation Programme under the Marie Skłodowska-Curie Action.