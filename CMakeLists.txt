# Minimum version of CMake required
cmake_minimum_required(VERSION 3.15)

# Project name and language
project(adeSolver LANGUAGES Fortran)

# Set the output directories for the build
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)

# Set compiler flags for Fortran
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -O3 -Wall -Wextra -std=f2008 -fPIC -fcheck=all -fbacktrace -fmax-errors=3 -g")


# # Check for HDF5 environment variables first
# if(DEFINED ENV{HDF5_INCLUDE} AND DEFINED ENV{HDF5_LIB})
#     message(STATUS "Using HDF5 paths from environment variables.")
#     set(HDF5_INCLUDE_DIR $ENV{HDF5_INCLUDE})
#     set(HDF5_LIB_DIR $ENV{HDF5_LIB})
#     set(HDF5_LIBRARIES "-L${HDF5_LIB_DIR} -lhdf5_fortran -lhdf5")
#     include_directories(${HDF5_INCLUDE_DIR})
# else()
#     message(STATUS "Using find_package to locate HDF5.")
#     find_package(HDF5 REQUIRED COMPONENTS Fortran)
#     if(HDF5_FOUND)
#         message(STATUS "Found HDF5: ${HDF5_LIBRARIES}")
#         include_directories(${HDF5_INCLUDE_DIRS})
#     endif()
# endif()


# Check if the environment variables for NetCDF are set, otherwise use nf-config
if(DEFINED ENV{NETCDF_INCLUDE} AND DEFINED ENV{NETCDF_LIB})
    message(STATUS "Using NetCDF paths from environment variables.")
    set(NC_INCLUDE_DIR $ENV{NETCDF_INCLUDE})
    set(NC_LIB_DIR $ENV{NETCDF_LIB})
    set(FLIBS_NC "-L${NC_LIB_DIR} -lnetcdff -lnetcdf")
else()
    message(STATUS "Using nf-config to get NetCDF paths.")
    find_program(NF_CONFIG_EXECUTABLE nf-config)
    if (NF_CONFIG_EXECUTABLE)
        execute_process(COMMAND ${NF_CONFIG_EXECUTABLE} --fflags OUTPUT_VARIABLE FFLAGS_NC OUTPUT_STRIP_TRAILING_WHITESPACE)
        execute_process(COMMAND ${NF_CONFIG_EXECUTABLE} --flibs OUTPUT_VARIABLE FLIBS_NC OUTPUT_STRIP_TRAILING_WHITESPACE)
        execute_process(COMMAND ${NF_CONFIG_EXECUTABLE} --includedir OUTPUT_VARIABLE NC_INCLUDE_DIR OUTPUT_STRIP_TRAILING_WHITESPACE)
        execute_process(COMMAND ${NF_CONFIG_EXECUTABLE} --prefix OUTPUT_VARIABLE NC_LIB_DIR OUTPUT_STRIP_TRAILING_WHITESPACE)
        set(NC_LIB_DIR "${NC_LIB_DIR}/lib")  # Ensure the library path is correct
    else()
        message(FATAL_ERROR "nf-config not found. Make sure NetCDF-Fortran is installed or set NETCDF_INCLUDE and NETCDF_LIB environment variables.")
    endif()
endif()

# Include the NetCDF-Fortran paths
include_directories(${NC_INCLUDE_DIR})
link_directories(${NC_LIB_DIR})


# Update the NetCDF linking to include HDF5
# set(FLIBS_NC "${FLIBS_NC} ${HDF5_LIBRARIES}")
message("FLIBS_NC ${FLIBS_NC}")

# Add the source directory where the .f90 files are located
add_subdirectory(src)

# Retrieve the Git version information
find_package(Git QUIET)
if (GIT_FOUND)
message("GIT_VERSION ${GIT_VERSION}")
execute_process(
    COMMAND ${GIT_EXECUTABLE} describe --tags --dirty
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    OUTPUT_VARIABLE GIT_VERSION
    OUTPUT_STRIP_TRAILING_WHITESPACE
    )
    set(GIT_VERSION ${GIT_VERSION})
    # Pass the Git version information to the src directory
    endif()
if (GIT_FOUND)
endif()

if(EXISTS GIT_VERSION)
    set(GIT_VERSION "unknown")
endif()
        
message("GIT_VERSION ${GIT_VERSION}")
    