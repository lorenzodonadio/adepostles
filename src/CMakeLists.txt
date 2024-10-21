# Set the executable name with versioning
set(EXECUTABLE_NAME "main-${GIT_VERSION}")

# Gather all the Fortran source files in the src directory
file(GLOB SRC_FILES *.f90)

# Create the executable from the Fortran source files
add_executable(${EXECUTABLE_NAME} ${SRC_FILES})

# Link necessary libraries for NetCDF
target_link_libraries(${EXECUTABLE_NAME} ${FLIBS_NC})
