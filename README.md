# adepostles

full docs [here](https://lorenzodonadio.github.io/adepostles/)

The **Advection Diffusion Equation Post Large Eddy Simulation** (adepostles) code solves tracer transport as a post processing step to LES simulations.

It was initially developped to treat output from [DALES](https://github.com/dalesteam/dales), but the idea is to be able to post process a broad set of LES outputs, as long as they conform to NetCDF.

## Features

- IBM
- MPI parallelization
- NetCDF input/output format

> NOTE:<br>
> The parallelization goal is to allow multiple initial conditions and source profiles to be calculated in parallel for the same LES simulation output
