# Run

After you have [compiled](installation.md) Adepostles, you are almost ready to run your first simulation.

It is adviced to copy the `main` executable to a separate directory where you want to run the simulation not to clutter build/bin

## Input files

The executable takes only **one config file** as input, this is usually called _namoptions_ in DALES. but can have any name. A complete description of _namoptions_ can be found [here](namoptions.md).

Within _namoptions_ we must specify 3 file paths:

1. **field_dump_path**  
   path to the fielddump, the output of the LES simulation merged into **one single NetCDF file**.
   The file must contain:

- time dimension
- x,y,z coordinates of the domain
- 4D fields (time,z,y,x): wind velocity (u,v,w), eddy viscosity (ekh)
- 2D vertical profiles (time, z) of density: rhobf

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

2. **sources_prefix**
3. **ibm_input_file**
