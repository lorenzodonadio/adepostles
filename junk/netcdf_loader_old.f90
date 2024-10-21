module netcdf_loader
   use netcdf
   use iso_fortran_env, only: real32
   use config, only: field_load_chunk_size
   implicit none

   real(real32), allocatable :: zt(:), zm(:), xt(:), xm(:), yt(:), ym(:),rtime(:)

contains

   subroutine load_static_variables(filename)
      character(len=*), intent(in) :: filename
      integer :: ncid, retval

      ! Dimension IDs and sizes
      integer :: zt_dim, zm_dim, xt_dim, xm_dim, yt_dim, ym_dim,time_dim
      integer :: zt_size, zm_size, xt_size, xm_size, yt_size, ym_size, time_size

      ! Variable IDs
      integer :: zt_varid, zm_varid, xt_varid, xm_varid, yt_varid, ym_varid,time_varid

      ! Static variable arrays (dynamic allocation)

      ! Open the NetCDF file
      retval = nf90_open(filename, NF90_NOWRITE, ncid)
      call nchandle_error(retval, 'Error opening file: '//trim(filename))

      ! Get dimension sizes for static variables
      call get_dimension_size(ncid, 'zt', zt_dim, zt_size)
      call get_dimension_size(ncid, 'zm', zm_dim, zm_size)
      call get_dimension_size(ncid, 'xt', xt_dim, xt_size)
      call get_dimension_size(ncid, 'xm', xm_dim, xm_size)
      call get_dimension_size(ncid, 'yt', yt_dim, yt_size)
      call get_dimension_size(ncid, 'ym', ym_dim, ym_size)
      call get_dimension_size(ncid, 'time', time_dim, time_size)

      ! Allocate memory for the static variables
      allocate(zt(zt_size), zm(zm_size), xt(xt_size), xm(xm_size), yt(yt_size), ym(ym_size),rtime(time_size))
      ! Debug output to verify allocation
      ! Read static variables from the file
      call get_and_read_variable(ncid, 'zt', zt_varid, zt)
      call get_and_read_variable(ncid, 'zm', zm_varid, zm)
      call get_and_read_variable(ncid, 'xt', xt_varid, xt)
      call get_and_read_variable(ncid, 'xm', xm_varid, xm)
      call get_and_read_variable(ncid, 'yt', yt_varid, yt)
      call get_and_read_variable(ncid, 'ym', ym_varid, ym)
      call get_and_read_variable(ncid, 'time', time_varid, rtime)

      ! Close the NetCDF file
      retval = nf90_close(ncid)
      call nchandle_error(retval, 'Error closing file: '//trim(filename))

      print *, 'Static variables loaded successfully from: ', trim(filename)
      call validate_time_size

   end subroutine load_static_variables

   subroutine load_time_dependent_variables(filename, time_start, time_end)
      character(len=*), intent(in) :: filename
      integer, intent(inout) :: time_start, time_end
      integer :: ncid, retval

      ! Dimension IDs and sizes
      integer :: time_dim, time_size

      ! Variable IDs for time-dependent variables
      integer :: time_varid, u_varid, v_varid, w_varid, ekh_varid

      ! Time-dependent variable arrays (dynamic allocation)
      real(real32), allocatable :: time(:)
      real(real32), allocatable :: u(:,:,:,:), v(:,:,:,:), w(:,:,:,:), ekh(:,:,:,:)

      ! Open the NetCDF file
      retval = nf90_open(filename, NF90_NOWRITE, ncid)
      call nchandle_error(retval, 'Error opening file: '//trim(filename))

      ! Get dimension size for the time dimension
      call get_dimension_size(ncid, 'time', time_dim, time_size)

      ! Check and adjust time range boundaries
      if (time_start < 1) time_start = 1
      if (time_end > time_size) time_end = time_size

      ! Allocate memory for the time range
      !   allocate(time(time_start:time_end))
      !   allocate(u(time_end - time_start + 1, :, :, :))
      !   allocate(v(time_end - time_start + 1, :, :, :))
      !   allocate(w(time_end - time_start + 1, :, :, :))
      !   allocate(ekh(time_end - time_start + 1, :, :, :))

      !   ! Read time-dependent variables from the file
      !   call get_and_read_variable_range(ncid, 'time', time_varid, time, time_start, time_end)
      !   call get_and_read_variable_range(ncid, 'u', u_varid, u, time_start, time_end)
      !   call get_and_read_variable_range(ncid, 'v', v_varid, v, time_start, time_end)
      !   call get_and_read_variable_range(ncid, 'w', w_varid, w, time_start, time_end)
      !   call get_and_read_variable_range(ncid, 'ekh', ekh_varid, ekh, time_start, time_end)

      !   ! Close the NetCDF file
      !   retval = nf90_close(ncid)
      !   if (retval /= nf90_noerr) then
      !      print *, 'Error closing file: ', trim(filename)
      !      call handle_netcdf_error(retval)
      !   end if

      print *, 'Time-dependent variables loaded successfully from ', trim(filename), ' for range: ', time_start, ' to ', time_end

      ! Close the NetCDF file
      retval = nf90_close(ncid)
      call nchandle_error(retval, 'Error closing file: '//trim(filename))
   end subroutine load_time_dependent_variables

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
      integer :: retval,var_type

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

   subroutine get_and_read_variable_range(ncid, var_name, varid, var_data, start_index, end_index)
      integer, intent(in) :: ncid, start_index, end_index
      character(len=*), intent(in) :: var_name
      integer, intent(out) :: varid
      real(real32), allocatable, intent(out) :: var_data(:,:,:,:)
      integer :: retval
      integer :: start(4), count(4)

      retval = nf90_inq_varid(ncid, trim(var_name), varid)
      call nchandle_error(retval, 'Error: Could not get variable ID for '//trim(var_name))

      start = (/ start_index - 1, 0, 0, 0 /)
      count = (/ end_index - start_index + 1, -1, -1, -1 /)

      retval = nf90_get_var(ncid, varid, var_data, start = start, count = count)
      call nchandle_error(retval, 'Error: Could not read data range for variable '//trim(var_name))
   end subroutine get_and_read_variable_range

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

   subroutine validate_time_size()
      if (mod(size(rtime), field_load_chunk_size) /= 0) then
         print *, 'time_size:', size(rtime)
         print *, 'field_load_chunk_size :', field_load_chunk_size
         stop 'Time size must be a multiple of Field load chunk, please adjust the options'
      endif
   end subroutine validate_time_size
end module netcdf_loader
