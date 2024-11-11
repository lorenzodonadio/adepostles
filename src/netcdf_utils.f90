module netcdf_utils
   use netcdf
   use iso_fortran_env, only: real32
   implicit none

contains
   subroutine nchandle_error(status, error_message)
      integer, intent(in) :: status
      character(len=*), optional, intent(in) :: error_message

      if (status /= nf90_noerr) then
         write(*,*) trim(nf90_strerror(status))
         if (present(error_message)) then
            write(*,*) error_message
         end if
         stop 'Execution stopped due to NetCDF error.'
      end if
   end subroutine nchandle_error

   subroutine get_dimension_size(ncid, dim_name, dimid, dim_size)
      integer, intent(in) :: ncid
      character(len=*), intent(in) :: dim_name
      integer, intent(out) :: dimid, dim_size
      integer :: retval
      character(len=256) :: name

      retval = nf90_inq_dimid(ncid, trim(dim_name), dimid)
      call nchandle_error(retval, 'Error: Could not get dimension ID for '//trim(dim_name))

      retval = nf90_inquire_dimension(ncid, dimid, name, dim_size)
      call nchandle_error(retval, 'Error: Could not get dimension size for '//trim(dim_name))
   end subroutine get_dimension_size

   subroutine get_1d_variable(ncid, var_name, varid, var_data)
      integer, intent(in) :: ncid
      character(len=*), intent(in) :: var_name
      integer, intent(out) :: varid
      real(real32), intent(out) :: var_data(:)
      integer :: retval
      ! integer :: retval,var_type

      retval = nf90_inq_varid(ncid, trim(var_name), varid)
      call nchandle_error(retval, 'Error: Could not get variable ID for '//trim(var_name))

      retval = nf90_get_var(ncid, varid, var_data)
      call nchandle_error(retval, 'Error: Could not read data for variable '//trim(var_name))
   end subroutine get_1d_variable

   subroutine get_profile_chunk(ncid, var_name, varid, var_data, chunk_number)
      use config, only: field_load_chunk_size

      integer, intent(in) :: ncid, chunk_number
      character(len=*), intent(in) :: var_name
      integer, intent(out) :: varid
      real(real32), intent(out) :: var_data(:,:)
      integer :: retval
      integer :: start(4)

      ! For debugging purposes
      ! call get_variable_dimensions(ncid,var_name)

      retval = nf90_inq_varid(ncid, trim(var_name), varid)
      call nchandle_error(retval, 'Error: Could not get variable ID for '//trim(var_name))
      start = (/1, 1, 1, field_load_chunk_size*(chunk_number-1)+1 /)

      ! write(*,*) 'Start', start
      ! write(*,*) 'var_data shape', shape(var_data)

      retval = nf90_get_var(ncid, varid, var_data, start = start)
      call nchandle_error(retval, 'Error: Could not read data range for variable '//trim(var_name))
   end subroutine get_profile_chunk
   subroutine get_field_chunk(ncid, var_name, varid, var_data, chunk_number)
      use config, only: field_load_chunk_size

      integer, intent(in) :: ncid, chunk_number
      character(len=*), intent(in) :: var_name
      integer, intent(out) :: varid
      real(real32), intent(out) :: var_data(:,:,:,:)
      integer :: retval
      integer :: start(4)

      ! For debugging purposes
      ! call get_variable_dimensions(ncid,var_name)

      retval = nf90_inq_varid(ncid, trim(var_name), varid)
      call nchandle_error(retval, 'Error: Could not get variable ID for '//trim(var_name))
      start = (/1, 1, 1, field_load_chunk_size*(chunk_number-1)+1 /)

      ! write(*,*) 'Start', start
      ! write(*,*) 'var_data shape', shape(var_data)

      retval = nf90_get_var(ncid, varid, var_data, start = start)
      call nchandle_error(retval, 'Error: Could not read data range for variable '//trim(var_name))
   end subroutine get_field_chunk

   subroutine get_variable_dimensions(ncid, var_name)
      use iso_fortran_env, only: int32
      integer, intent(in) :: ncid
      character(len=*), intent(in) :: var_name
      integer :: retval, varid, num_dims, dim_ids(nf90_max_var_dims)
      integer :: dim_sizes(nf90_max_var_dims), i
      character(len=128) :: dim_name

      ! Get the variable ID
      retval = nf90_inq_varid(ncid, trim(var_name), varid)
      call nchandle_error(retval, 'Error: Could not get variable ID for '//trim(var_name))

      ! Get the number of dimensions and their IDs for the specified variable
      retval = nf90_inquire_variable(ncid, varid, ndims=num_dims, dimids=dim_ids)
      call nchandle_error(retval, 'Error: Could not inquire variable dimensions for '//trim(var_name))

      ! Loop through each dimension to get its size and name
      write(*,*) 'Variable: ', trim(var_name), ' has ', num_dims, ' dimensions:'
      do i = 1, num_dims
         retval = nf90_inquire_dimension(ncid, dim_ids(i), dim_name, dim_sizes(i))
         call nchandle_error(retval, 'Error: Could not inquire size of dimension '//trim(dim_name))
         write(*,*) 'Dimension ', i, ': ', trim(dim_name), ' Size: ', dim_sizes(i)
      end do
   end subroutine get_variable_dimensions

   subroutine create_concentration_file_nc(ncid, filename, x, y, z, nsv)
      integer, intent(out) :: ncid
      character(len=256), intent(in) :: filename
      real, intent(in) :: x(:), y(:), z(:)
      integer, intent(in) :: nsv  ! Number of scalar fields (pollutants)
      integer :: retval, time_dim, x_dim, y_dim, z_dim, nsv_dim, varid, xid, yid, zid

      ! Create the NetCDF file
      retval = nf90_create(trim(filename), nf90_clobber, ncid)
      call nchandle_error(retval, 'Error: Could not create NetCDF file '//trim(filename))

      ! Define dimensions
      retval = nf90_def_dim(ncid, "x", size(x), x_dim)
      call nchandle_error(retval, 'Error: Could not define x dimension')

      retval = nf90_def_dim(ncid, "y", size(y), y_dim)
      call nchandle_error(retval, 'Error: Could not define y dimension')

      retval = nf90_def_dim(ncid, "z", size(z), z_dim)
      call nchandle_error(retval, 'Error: Could not define z dimension')

      retval = nf90_def_dim(ncid, "time", nf90_unlimited, time_dim)
      call nchandle_error(retval, 'Error: Could not define time dimension')

      retval = nf90_def_dim(ncid, "nsv", nsv, nsv_dim)  ! Define the nsv dimension
      call nchandle_error(retval, 'Error: Could not define nsv dimension')

      ! Define variables for the dimensions x, y, z, and time
      retval = nf90_def_var(ncid, "x", nf90_real, x_dim, xid)
      call nchandle_error(retval, 'Error: Could not define variable x')

      retval = nf90_def_var(ncid, "y", nf90_real, y_dim, yid)
      call nchandle_error(retval, 'Error: Could not define variable y')

      retval = nf90_def_var(ncid, "z", nf90_real, z_dim, zid)
      call nchandle_error(retval, 'Error: Could not define variable z')

      retval = nf90_def_var(ncid, "time", nf90_real, time_dim, varid)
      call nchandle_error(retval, 'Error: Could not define variable time')

      ! Define the 5D concentration variable with time as the first dimension: c(time, x, y, z, nsv)
      !   retval = nf90_def_var(ncid, "c", nf90_real, (/time_dim, x_dim, y_dim, z_dim, nsv_dim/), varid)
      retval = nf90_def_var(ncid, "c", nf90_real, (/x_dim, y_dim, z_dim, nsv_dim, time_dim/), varid)
      call nchandle_error(retval, 'Error: Could not define concentration variable c')

      ! End definition mode
      retval = nf90_enddef(ncid)
      call nchandle_error(retval, 'Error: Could not end definition mode in NetCDF file')

      ! Write x, y, and z to the corresponding variables x, y, z
      retval = nf90_put_var(ncid, xid, x)
      call nchandle_error(retval, 'Error: Could not write values to variable x')

      retval = nf90_put_var(ncid, yid, y)
      call nchandle_error(retval, 'Error: Could not write values to variable y')

      retval = nf90_put_var(ncid, zid, z)
      call nchandle_error(retval, 'Error: Could not write values to variable z')

   end subroutine create_concentration_file_nc

   subroutine write_concentration_nc(ncid, c, rsts)
      integer, intent(in) :: ncid
      real, intent(in) :: rsts
      real, intent(in) :: c(:,:,:,:)  ! 4D array with dimensions (x, y, z, nsv)
      integer :: retval, time_size, varid, time_dimid
      integer :: start(5)
      character(len=72) :: varname

      ! Get the dimension ID for "time"
      retval = nf90_inq_dimid(ncid, "time", time_dimid)
      call nchandle_error(retval, 'Error: Could not get dimension ID for time')

      ! Get the current size of the "time" dimension to determine the next index
      retval = nf90_inquire_dimension(ncid, time_dimid, varname, len=time_size)
      call nchandle_error(retval, 'Error: Could not get the current size of time dimension')

      ! Set up the start indices for writing to c(x, y, z, nsv, time)
      start = (/1, 1, 1, 1,time_size + 1/)

      ! Write the concentration field at the current time index for all nsv
      retval = nf90_inq_varid(ncid, 'c', varid)
      call nchandle_error(retval, 'Error: Could not get variable ID for c')

      retval = nf90_put_var(ncid, varid, c, start = start)
      call nchandle_error(retval, 'Error: Could not write variable c')

      ! Write the current time step to the time variable
      retval = nf90_inq_varid(ncid, 'time', varid)
      call nchandle_error(retval, 'Error: Could not get variable ID for time')

      retval = nf90_put_var(ncid, varid, rsts, start = (/time_size + 1/))
      call nchandle_error(retval, 'Error: Could not write time step to variable time')

   end subroutine write_concentration_nc



end module netcdf_utils
! !!!-------------- WRITING CONCENTRATION ------------

! !TODO: add an extra input, integer nsv, that means the number of scalar fields, so now we will write the concentration
! ! to a 5D field: x,y,time,nsv, so that we can handle many different pollutants in the same simulation

! subroutine create_concentration_file_nc(ncid,filename,x,y,z)
!    integer, intent(out) :: ncid
!    character(len=256), intent(in) :: filename
!    real, intent(in) :: x(:),y(:),z(:)
!    integer :: retval, time_dim, x_dim, y_dim, z_dim, varid,xid,yid,zid

!    ! Create the NetCDF file
!    retval = nf90_create(trim(filename), nf90_clobber, ncid)
!    call nchandle_error(retval, 'Error: Could not create NetCDF file '//trim(filename))

!    ! Define dimensions

!    retval = nf90_def_dim(ncid, "x", size(x), x_dim)
!    call nchandle_error(retval, 'Error: Could not define x dimension')

!    retval = nf90_def_dim(ncid, "y", size(y), y_dim)
!    call nchandle_error(retval, 'Error: Could not define y dimension')

!    retval = nf90_def_dim(ncid, "z", size(z), z_dim)
!    call nchandle_error(retval, 'Error: Could not define z dimension')

!    retval = nf90_def_dim(ncid, "time", nf90_unlimited, time_dim)
!    call nchandle_error(retval, 'Error: Could not define time dimension')

!    !define actual variables containing the values of the dimensions since they might not be unifornmly distributed
!    retval = nf90_def_var(ncid, "x", nf90_real, x_dim , xid)
!    call nchandle_error(retval, 'Error: Could not define variable x')
!    retval = nf90_def_var(ncid, "y", nf90_real, y_dim , yid)
!    call nchandle_error(retval, 'Error: Could not define variable y')
!    retval = nf90_def_var(ncid, "z", nf90_real, z_dim , zid)
!    call nchandle_error(retval, 'Error: Could not define variable z')
!    retval = nf90_def_var(ncid, "time", nf90_real, time_dim , varid)
!    call nchandle_error(retval, 'Error: Could not define variable time')
!    ! Define the variable with dimensions, concentration only for now: c(x, y, z,time)
!    retval = nf90_def_var(ncid, "c", nf90_real, (/x_dim, y_dim, z_dim, time_dim/) , varid)
!    call nchandle_error(retval, 'Error: Could not define variable c')

!    ! End definition mode
!    retval = nf90_enddef(ncid)
!    call nchandle_error(retval, 'Error: Could not end definition mode in NetCDF file')

!    ! Write x, y, and z to the corresponding variables x, y, z
!    retval = nf90_put_var(ncid, xid, x)
!    call nchandle_error(retval, 'Error: Could not write values to variable x')
!    retval = nf90_put_var(ncid, yid, y)
!    call nchandle_error(retval, 'Error: Could not write values to variable y')
!    retval = nf90_put_var(ncid, zid, z)
!    call nchandle_error(retval, 'Error: Could not write values to variable z')

! end subroutine create_concentration_file_nc

! subroutine write_concentration_nc(ncid,c,rsts)
!    ! use modglobal, only: rsts  ! real simulation time in seconds
!    ! use modfields, only: c0    ! 3D field to be written
!    integer, intent(in) :: ncid
!    real, intent(in) :: rsts
!    real, intent(in) :: c(:,:,:)
!    integer :: retval, time_size, varid, time_dimid
!    integer :: start(4)
!    character(len=72) :: varname
!    ! Get the dimension ID for "time"
!    retval = nf90_inq_dimid(ncid, "time", time_dimid)
!    call nchandle_error(retval, 'Error: Could not get dimension ID for time')

!    ! Get the current size of the "time" dimension to determine the next index
!    retval = nf90_inquire_dimension(ncid, time_dimid, varname, len=time_size)
!    call nchandle_error(retval, 'Error: Could not get the current size of time dimension')
!    ! Write concentration field c0 at the current time index
!    start = (/1, 1, 1, time_size + 1/)
!    ! start = (/0, 0, 0, time_size/)
!    retval = nf90_inq_varid(ncid, 'c', varid)
!    call nchandle_error(retval, 'Error: Could not get variable ID for c')

!    retval = nf90_put_var(ncid, varid, c, start = start)
!    call nchandle_error(retval, 'Error: Could not write variable c0')

!    ! Write the current time step to the time variable
!    retval = nf90_inq_varid(ncid, 'time', varid)
!    call nchandle_error(retval, 'Error: Could not get variable ID for time')

!    retval = nf90_put_var(ncid, varid, rsts, start = (/time_size + 1/))
!    call nchandle_error(retval, 'Error: Could not write time step to variable time')
! end subroutine write_concentration_nc
