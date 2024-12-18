# Installation

Let us run throught the installation of ADEpostLES

## Requirements

- **Fortran compiler**: recomended gfortran
- **Cmake**: recomended to just download and install the [binary](https://cmake.org/download/)
- **MPI**: tested with Open MPI > 4.1.4. Installation guides: [1](https://webpages.charlotte.edu/abw/coit-grid01.uncc.edu/ParallelProgSoftware/Software/OpenMPIInstall.pdf),[2](https://docs.open-mpi.org/en/main/installing-open-mpi/quickstart.html)
- **NetCDF**: Netcdf fortran, recomended installation using [conda](https://anaconda.org/conda-forge/netcdf-fortran)  
  **NOTE** You might need to export the paths to your netcdf-fortran instalations, you can find them using `nf-config`. See the following example, cmake compilation first looks for these environment variables.

```
export NETCDF_INCLUDE=/home/user/miniconda3/envs/fort/include
export NETCDF_LIB=/home/user/miniconda3/envs/fort/lib
```

```
#!/bin/sh

module load 2023r1-gcc11
module load openmpi/4.1.4
module load cmake/3.24.3
module load netcdf-fortran/4.6.0
```

## Compilation

1. Clone the github repo
   `git clone https://github.com/lorenzodonadio/adepostles.git`
2. Create a build directory
