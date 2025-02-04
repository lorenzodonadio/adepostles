# Config file sample:

Sample _namoptions_ file:

```
&RUN
runtime    =  700
dtmax      =  1
ladaptivedt  =  .true.
lanisotrop = .true.
rkmethod = 2
field_load_chunk_size = 2
output_save_interval = 5
courant_limit = 0.8
vonneumann_limit = 2

field_dump_path = "fielddump_30s.nc"
sources_prefix = "tracer_inp/tracer_inp_001_001.nc"
outputfile_path = "must_conc_30s_periodictest.nc"
/

&BOUNDARY
xboundary = 12
yboundary = 12
/

&IBM
lapplyibm = .true.
ibm_input_file = "ibm.inp.001"
/
```


## RUN Section
Parameters controlling the core runtime behavior and input/output settings.

#### Time Control
- `runtime`: duration of the simulation is seconds
- `dtmax` : max timestep of the simulation in seconds
- `ladaptivedt`: (optional default true) switch to use adaptive timestep, if tuned off then dtmax is the timestep for the simulation.
- `courant_limit`: (optional) CFL limit, adviced around 0.3 for complex geometries (defaults depend on RK scheme)
- `vonneumann_limit`: (optional) von neumann limit, usually not the limiting factor for the timestep but adviced around 1 (defaults depend on RK scheme),
- `output_save_interval`: number of seconds between result saves (default 20)

#### Simulation Parameters
- `lanisotrop`: (default true), use anisotrophic diffusion scheme (moddiff.f90)
- `rkmethod`: whick Runge-Kutta method to use: 1. euler, 2.heun, 3.rk3, 4. rk4. (default 1)
- `field_load_chunk_size`: integer that controls how many time steps of the wind field are loaded; how big is the chunk. The larger it is, the faster the simulation (less i/o) but the more memory requirements.

#### File Configuration
- `field_dump_path`: path for the LES output .nc file ([example](/run/#example-files))
- `sources_prefix`: file or prefix of the tracer source .nc file ([example](/run/#example-files))
- `outputfile_path`: path where the results will be saved. This filename will be suffixed with .001, .002, etc for the processor number in case parallel simulations are run.

## BOUNDARY Section

possible values:  11: Neumann 1st Order, 12: Neumann 2nd Order, 2: Periodic. It is not recomended to use 1st order BC

- `xboundary`: boundary condition for east-west direction
- `yboundary`: boundary condition for north-south direction


## IBM Section
Immersed Boundary Method configuration.

- `lapplyibm`: boolean to use ibm or not
- `ibm_input_file`: path to the ibm file
