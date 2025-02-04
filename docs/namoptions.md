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
