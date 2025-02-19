module config
   implicit none

   ! < File unit numbers
   integer, parameter :: ifconfig = 10      ! Unit number for the configuration file
   integer, parameter :: ifinputibm  = 20      ! Unit number for input file
   integer, parameter :: ifoutput = 30      ! Unit number for output file
   integer, parameter :: ifnamopt = 40      ! Unit number for optional namelist files

   ! Configuration variables for &RUN namelist
   ! /RUN/

   character(len=256) :: field_dump_path = ''  ! Must be specified
   character(len=256) :: sources_prefix = ''  ! Must be specified
   character(len=256) :: outputfile_path = 'concentration_out.nc'  ! Must be specified
   character(len=256) :: ibm_input_file = ''   ! Must be specified
   character(len=256) :: config_filename


   logical :: ladaptivedt = .false.          ! Sensible default
   logical :: lanisotrop = .true.
   logical :: lperiodic_field_pad = .true.  !< if to use periodic BC for reconstructing the field ghost cells, if set to false it uses interpolation

   integer :: rkmethod = 1                   !< 1. euler, 2.heun, 3.rk3, 4. rk4
   integer :: runtime = 25                   !< seconds
   integer :: field_load_chunk_size = 100    ! Sensible default

   real :: dtmax = 1.                      !< seconds
   real :: output_save_interval = 20. !< how many seconds between save of concentration output

   real :: courant_limit = -1 !< CFL criterion, default is -1 so we use the hard coded value in time_integrate.f90
   real :: vonneumann_limit = -1  !< Von Neuman criterion, default is -1 so we use the hard coded value in time_integrate.f90

   ! Configuration variables for &IBM namelist
   ! /IBM/
   logical :: lapplyibm = .false.                 ! Sensible default


   !other flags

   ! /BOUNDARY/
   ! 11: Neumann 1st Order, 12: Neumann 2nd Order, 2: Periodic
   integer :: xboundary = 1
   integer :: yboundary = 2
contains

   subroutine read_config_file()
      integer :: ios
      ! Namelists for the configuration
      namelist /RUN/ runtime, dtmax,output_save_interval,rkmethod, ladaptivedt, lanisotrop,lperiodic_field_pad,&
         field_load_chunk_size, field_dump_path,sources_prefix, outputfile_path, courant_limit, vonneumann_limit

      namelist /IBM/ lapplyibm, ibm_input_file

      namelist /BOUNDARY/ xboundary, yboundary

      if (command_argument_count() >=1) then
         call get_command_argument(1,config_filename)
      end if

      write(*,*) 'Config path:', trim(config_filename)
      ! Open the configuration file and read the namelists using the named unit number ifconfig
      open(unit=ifconfig, file=config_filename, status='old', action='read', iostat=ios)
      if (ios /= 0) then
         write(*,*) 'Error: Unable to open configuration file ', trim(config_filename)
         stop 'Execution halted due to missing configuration file.'
      end if

      ! Read the &RUN namelist
      read(ifconfig, nml=RUN, iostat=ios)
      if (ios /= 0) then
         write(*,*) 'Error:', ios, ' In file: ', trim(config_filename)
         stop 'Execution halted due to errors in the &RUN namelist.'
      end if

      ! Read the &IBM namelist
      ! rewind(ifconfig)  ! Rewind to ensure the file pointer is at the beginning for a fresh read not needed here
      read(ifconfig, nml=IBM, iostat=ios)
      if (ios /= 0) then
         write(*,*) 'Error:', ios, ' In file: ', trim(config_filename)
         stop 'Execution halted due to errors in the &IBM namelist.'
      end if

      ! Close the configuration file
      close(ifconfig)

      ! Check required parameters and set defaults where necessary
      call validate_parameters()
   end subroutine read_config_file

   subroutine validate_parameters()
      ! Ensure all necessary parameters have been set and provide default values or errors if required

      ! Validate &RUN namelist parameters

      if (trim(field_dump_path) == '') then
         write(*,*) 'Error: Parameter field_dump_path is required in the &RUN namelist but was not provided.'
         stop 'Execution halted due to missing field_dump_path.'
      end if

      ! Validate &IBM namelist parameters
      if (lapplyibm .and. trim(ibm_input_file) == '') then
         write(*,*) 'Error: Parameter ibm_input_file is required in the &IBM namelist but was not provided.'
         stop 'Execution halted due to missing ibm_input_file.'
      end if

      if (rkmethod > 4 .or. rkmethod < 1) then
         write(*,*) 'Parameter rkmethod not supported: only 1 = euler, 2 = heun, 3 = rk3, 4 = rk4'
         stop 'Execution halted due to wrong rkmethod'
      end if
      write(*,*) 'Configuration loaded successfully:'
      write(*,*) 'runtime:', runtime
      write(*,*) 'dtmax:', dtmax
      write(*,*) 'ladaptivedt:', ladaptivedt
      write(*,*) 'rkmethod:', rkmethod
      write(*,*) 'field_load_chunk_size:', field_load_chunk_size
      write(*,*) 'field_dump_path:', trim(field_dump_path)
      write(*,*) 'lapplyibm:', lapplyibm
      write(*,*) 'ibm_input_file:', trim(ibm_input_file)
   end subroutine validate_parameters

end module config
