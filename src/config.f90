module config
   implicit none

   ! < File unit numbers
   integer, parameter :: ifconfig = 10      ! Unit number for the configuration file
   integer, parameter :: ifinput  = 20      ! Unit number for input file
   integer, parameter :: ifoutput = 30      ! Unit number for output file
   integer, parameter :: ifnamopt = 40      ! Unit number for optional namelist files

   ! Configuration variables for &RUN namelist
   integer :: iexpnr = -1                    ! Default value indicates it must be specified
   integer :: runtime = 25                   ! Sensible default
   integer :: dtmax = 1                      ! Sensible default
   logical :: ladaptivedt = .false.          ! Sensible default
   logical :: lanisotrop = .true.
   integer :: integration_scheme = 1         ! Sensible default
   integer :: field_load_chunk_size = 100    ! Sensible default
   character(len=256) :: field_dump_path = ''  ! Must be specified

   ! Configuration variables for &IBM namelist
   logical :: libm = .false.                 ! Sensible default
   character(len=256) :: ibm_input_file = ''   ! Must be specified

contains

   subroutine read_config_file()
      character(len=256) :: config_filename
      integer :: ios
      ! Namelists for the configuration
      namelist /RUN/ iexpnr, runtime, dtmax, ladaptivedt,lanisotrop, integration_scheme, field_load_chunk_size, field_dump_path
      namelist /IBM/ libm, ibm_input_file


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
      if (iexpnr == -1) then
         write(*,*) 'Error: Parameter iexpnr is required in the &RUN namelist but was not provided.'
         stop 'Execution halted due to missing iexpnr.'
      end if

      if (trim(field_dump_path) == '') then
         write(*,*) 'Error: Parameter field_dump_path is required in the &RUN namelist but was not provided.'
         stop 'Execution halted due to missing field_dump_path.'
      end if

      ! Validate &IBM namelist parameters
      if (libm .and. trim(ibm_input_file) == '') then
         write(*,*) 'Error: Parameter ibm_input_file is required in the &IBM namelist but was not provided.'
         stop 'Execution halted due to missing ibm_input_file.'
      end if

      write(*,*) 'Configuration loaded successfully:'
      write(*,*) 'iexpnr:', iexpnr
      write(*,*) 'runtime:', runtime
      write(*,*) 'dtmax:', dtmax
      write(*,*) 'ladaptivedt:', ladaptivedt
      write(*,*) 'integration_scheme:', integration_scheme
      write(*,*) 'field_load_chunk_size:', field_load_chunk_size
      write(*,*) 'field_dump_path:', trim(field_dump_path)
      write(*,*) 'libm:', libm
      write(*,*) 'ibm_input_file:', trim(ibm_input_file)
   end subroutine validate_parameters

end module config
