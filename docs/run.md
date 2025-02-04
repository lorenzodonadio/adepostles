# Run

After you have [compiled](installation.md) Adepostles, you are almost ready to run your first simulation.

It is adviced to copy the `main` executable to a separate directory where you want to run the simulation not to clutter build/bin

## Input files

The executable takes only **one config file** as input, this is usually called _namoptions_ in DALES. but can have any name. A complete description of _namoptions_ can be found [here](namoptions.md).

Within _namoptions_ we must specify 3 file paths:

**1. field_dump_path**  
path to the fielddump, the output of the LES simulation merged into **one single NetCDF file**.
The file must contain:

- time dimension
- x,y,z coordinates of the domain
- 4D fields (time,z,y,x): wind velocity (u,v,w), eddy viscosity (ekh)
- 2D vertical profiles (time, z) of density: rhobf

**2. sources_prefix**

_sources_prefix_ specifies the source file for the concentration, which has both initial conditions and time varying sources.

For runs with **one single core** then _source_prefix_ must be exactly the path to the source file and end with .nc

For **runs in parallel** ,i.e with multiple source files, then only the prefix to the NetCDF file should be provided, and the program will automatically append `_001.nc` for the first core and `_002.nc` for the second core and so on. So you must prepare source files that end exactly that way up to the number of cores you are using for the run.

**3. ibm_input_file**

The IBM file is **not NetCDF** only because DALES still takes a text file as input for the IBM and it is important that the same file use for IBM in the LES is used for Adepostles.

## Example files

**field_dump**
Here is an example of the result of `ncdump -h` of a fielddump with the expected format:

```
dimensions:
        time = 1480 ;
        zt = 20 ;
        zm = 20 ;
        xt = 128 ;
        xm = 128 ;
        yt = 128 ;
        ym = 128 ;
variables:
        float time(time) ;
        float zt(zt) ;
        float zm(zm) ;
        float xt(xt) ;
        float xm(xm) ;
        float yt(yt) ;
        float ym(ym) ;
        float u(time, zt, yt, xm) ;
        float v(time, zt, ym, xt) ;
        float w(time, zm, yt, xt) ;
        float ekh(time, zt, yt, xt) ;
        float rhof(time, zt) ;
        float rhobf(time, zt) ;
        float rhobh(time, zm) ;
        float presh(time, zt) ;
```

**sources**
Example of netcdf structure for a source file, that has 2 separate groups (set of tracers) each with only one source

init is left empty so it means 0 everywhere

one could imagine multiple sources per tracer

some utilities to create these source files: TODO insert link

```
netcdf tracer_inp_001_001 {
dimensions:
        time = 10 ;
        x = 256 ;
        y = 256 ;
        z = 40 ;

// global attributes:
                :description = "NetCDF file for tracer properties" ;
                :history = "Created with T_tracer structure equivalent" ;
                :numtracers = 2LL ;
                :totnumsources = 2LL ;

group: s0 {
  variables:
        float init(z, y, x) ;
        float source_1(time) ;
                source_1:x = 8LL ;
                source_1:y = 30LL ;
                source_1:z = 2LL ;

  // group attributes:
                :tracname = "s0" ;
                :traclong = "s0_x_8_y_30" ;
                :unit = "kg/kg" ;
                :molar_mass = -999. ;
                :lemis = "true" ;
                :numsources = 1LL ;
  } // group s0

group: s1 {
  variables:
        float init(z, y, x) ;
        float source_1(time) ;
                source_1:x = 8LL ;
                source_1:y = 45LL ;
                source_1:z = 2LL ;

  // group attributes:
                :tracname = "s1" ;
                :traclong = "s1_x_8_y_45" ;
                :unit = "kg/kg" ;
                :molar_mass = -999. ;
                :lemis = "true" ;
                :numsources = 1LL ;
  } // group s1

}
```
