module netcdf_loader
   use netcdf
   use iso_fortran_env, only: real32
   implicit none

contains

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

   subroutine get_and_read_variable(ncid, var_name, varid, var_data)
      integer, intent(in) :: ncid
      character(len=*), intent(in) :: var_name
      integer, intent(out) :: varid
      real(real32), intent(out) :: var_data(:)
      integer :: retval
      ! integer :: retval,var_type

      retval = nf90_inq_varid(ncid, trim(var_name), varid)
      call nchandle_error(retval, 'Error: Could not get variable ID for '//trim(var_name))

      ! retval = nf90_inquire_variable(ncid, varid, xtype=var_type)
      ! call nchandle_error(retval, 'Error: Could not inquire variable type for '//trim(var_name))

      ! ! Check the variable type and read the data based on the type
      ! select case (var_type)
      !  case (nf90_float)
      !    print *, 'FLOAT'
      !  case (nf90_double)
      !    print *, 'DOUBLE'
      !  case default
      !    print *, "UNKNOWN", var_type
      ! end select

      retval = nf90_get_var(ncid, varid, var_data)
      call nchandle_error(retval, 'Error: Could not read data for variable '//trim(var_name))
      print *, 'READ: ', var_name
   end subroutine get_and_read_variable

   subroutine get_and_read_variable_chunk(ncid, var_name, varid, var_data, chunk_number)
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
   end subroutine get_and_read_variable_chunk

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




end module netcdf_loader
