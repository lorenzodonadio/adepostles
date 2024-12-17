# Welcome to ADEpostLES Docs

The Advection Diffusion Equation Post Large Eddy Simulation (adepostles) code solves tracer transport as a post processing step to LES simulations.

It was initially developped to treat output from [DALES](https://github.com/dalesteam/dales), but the idea is to be able to post process a broad set of LES outputs, as long as they conform to NetCDF.

## Commands

- `mkdocs new [dir-name]` - Create a new project.
- `mkdocs serve` - Start the live-reloading docs server.
- `mkdocs build` - Build the documentation site.
- `mkdocs -h` - Print help message and exit.

## Project layout

    mkdocs.yml    # The configuration file.
    docs/
        index.md  # The documentation homepage.
        ...       # Other markdown pages, images and other files.
