# NetCDF Field Merger and Profile Adder

## Overview

This utility provides a comprehensive solution for merging distributed NetCDF files from parallel computational simulations and adding profile variables to the merged dataset.

## Components

### 1. `fieldmerge.py`
A Python script with two primary functionalities:
- Merging distributed NetCDF files from multiple processing tiles
- Adding profile variables to an existing NetCDF dataset

### 2. `merge.sh`
A shell script for submitting the field merging job to a SLURM-managed high-performance computing (HPC) environment.

## Key Classes and Functions

### `PostProcessField` Class
Manages the merging of NetCDF files from multiple processing tiles.

#### Methods
- `__init__(file_path)`: Initializes the merger by discovering and organizing NetCDF files
- `extract_field(var_name, time_start, time_end)`: Extracts variable data across all tiles for a specified time range
- `save_to_single_netcdf(output_file, chunk_size)`: Merges all discovered NetCDF files into a single file

### Main Functions
- `add_profiles_to_dataset(profile_path, dataset_path)`: Adds profile variables to an existing NetCDF file

## Usage

### Command-line Interface

Two primary sub-commands are available:

1. Merge NetCDF Files
```bash
python fieldmerge.py merge --input_dir /path/to/netcdf/files \
                           --output_file completefielddump.nc \
                           --chunk_size 50 \
                           --profile_file optional_profiles.nc
```

2. Add Profiles to Existing Dataset
```bash
python fieldmerge.py add_profiles --profile_file profiles.nc \
                                  --dataset_file existing_dataset.nc
```

### HPC Deployment

The `merge.sh` script provides a SLURM job configuration for running the merge process:

```bash
sbatch merge.sh
```

## Configuration Parameters

### Merge Command
- `--input_dir`: Directory containing distributed NetCDF files
- `--output_file`: Path for the merged output file
- `--chunk_size`: Number of time steps to process in each iteration (default: 50)
- `--profile_file`: Optional profile file to merge with the dataset

### Add Profiles Command
- `--profile_file`: Path to the profile NetCDF file
- `--dataset_file`: Path to the existing NetCDF dataset

## Requirements

- Python 3.x
- Libraries:
  - `numpy`
  - `netCDF4`
  - `tqdm`

## Limitations and Considerations

- Assumes a specific file naming convention for distributed NetCDF files
- Requires matching time dimensions between profile and main dataset
- Memory usage depends on `chunk_size` and total dataset size

## Notes

- Supports merging files from multiple processing tiles
- Allows optional addition of profile variables post-merging

## Example Workflow

1. Run parallel simulation generating multiple NetCDF files
2. Use `fieldmerge.py merge` to combine files into a single dataset
3. Optionally add profile variables using the same script


# NetCDF Field Subsampling Utility

## Overview

This utility provides a flexible tool for subsampling large NetCDF files, allowing users to extract time-series data at different sampling rates while preserving the original dataset's structure and attributes.

## Components

### 1. `fieldsubsample.py`
A Python script for subsampling NetCDF files with the following key features:
- Extracting time series data at specified sampling intervals
- Skipping initial non-physical time steps
- Batch processing of large datasets
- Preserving original variable attributes and global metadata

### 2. `subsample.sh`
A SLURM job script for submitting the field subsampling process to a high-performance computing (HPC) environment.

## Key Functions

### `subsample_netcdf()`
Primary function for subsampling NetCDF files

#### Parameters
- `input_file` (str): Path to the source NetCDF file
- `output_dir` (str, optional): Directory to save subsampled files (default: current directory)
- `skip_first` (int, optional): Number of initial time steps to skip (default: 120)
- `sampling_rates` (list, optional): List of sampling rates in seconds (default: [5, 10, 30, 60])
- `nbatches` (int, optional): Number of batches to process data (default: 10)

## Usage

### Command-line Interface

```bash
python fieldsubsample.py input_file.nc \
    --output_dir ./output \
    --skip_first 432 \
    --sampling_rates 1,5,10,30,60 \
    --nbatches 10
```

### HPC Deployment

Use the provided `subsample.sh` SLURM script:

```bash
sbatch subsample.sh
```

## Configuration Options

### Command-line Arguments
- `input_file`: Required input NetCDF file path
- `--output_dir`: Output directory for subsampled files
- `--skip_first`: Number of initial time steps to skip
- `--sampling_rates`: Comma-separated sampling rates (in seconds)
- `--nbatches`: Number of data processing batches

## Workflow

1. Load the original NetCDF file
2. Skip specified initial time steps
3. Create subsampled files at different sampling rates
4. Preserve original variable attributes and metadata
5. Save output files with rate-specific naming (e.g., `fielddump_10s.nc`)

## Requirements

- Python 3.x
- Libraries:
  - `numpy`
  - `netCDF4`
  - `tqdm`

## Key Features

- Flexible sampling rate configuration
- Batch processing for large datasets
- Preserves original data structure
- Handles time-dependent and time-independent variables
- Supports multiple output files for different sampling intervals

## Limitations and Considerations

- Assumes uniform time step in source data
- Memory usage depends on dataset size and batch configuration
- Requires selection of skip_first and sampling rates
