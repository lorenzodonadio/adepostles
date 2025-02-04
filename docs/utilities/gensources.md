# Tracer Generation Utility for NetCDF Files

## Overview

This utility is designed to generate tracer source profiles for atmospheric or environmental modeling, specifically creating NetCDF files with point source emissions across multiple spatial and temporal dimensions.

## Components

### 1. `sourcegen.py`
A Python module that provides core functionality for creating NetCDF tracer profiles:
- Creates NetCDF files with tracer information
- Supports adding multiple tracers
- Allows definition of point sources for each tracer
- Manages dimensions and metadata for tracer profiles

### 2. `must_tracergen.py` (example)
The main script for generating tracer input files, **modify this to match your needs**:
- Generates multiple tracer profiles based on predefined configurations
- Supports parallel processing with configurable parameters
- Creates multiple NetCDF files with point source emissions

### 3. `gensource.sh`
A shell script for running the tracer generation on a high-performance computing (HPC) environment:
- Uses SLURM job scheduling
- Configures environment and dependencies
- Executes the tracer generation script

## Configuration Parameters

Key configuration parameters in `must_tracergen.py`:
- `tracers_per_core`: Number of tracers per computational core (default: 3)
- `nproc`: Number of processors (default: 2)
- `nruns`: Number of simulation runs (default: 11)
- `z_idx_src`: Vertical index for point sources (default: 2)
- `step_src_mat`: Source profile matrix defining emission values over time
- `point_sources`: Array of point source coordinates

## Usage

### Prerequisites
- Python 3.x
- Required libraries: 
  - `numpy`
  - `netCDF4`
  - `tqdm`

### Running the Script

```bash
python must_tracergen.py <fielddump_netcdf_path> --output_dir <output_directory>
```

Example:
```bash
python must_tracergen.py fielddump_10s.nc --output_dir ./tracers/
```

### HPC Deployment

The provided `gensource.sh` script can be used to submit the job to a SLURM-managed cluster:

```bash
sbatch gensource.sh
```

## Output

The script generates multiple NetCDF files in the specified output directory:
- Filename format: `tracer_inp_XXX_YYY.nc`
  - `XXX`: Run number (001-011)
  - `YYY`: Processor number (001-002)

Each file contains:
- Tracer metadata
- Point source emission profiles
- Spatial and temporal emission information

## Customization

To modify the tracer generation:
1. Adjust `point_sources` array to change source locations
2. Modify `step_src_mat` to alter emission profiles
3. Change `tracers_per_core`, `nproc`, and `nruns` to match your simulation requirements

## Note

This utility is specifically designed for atmospheric or environmental modeling workflows that require detailed tracer source definition across multiple dimensions.