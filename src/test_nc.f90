module test_nc
   use netcdf
   use iso_fortran_env, only: real32
   implicit none
contains
   subroutine do_test_nc(filename)
      character(len=256), intent(in) :: filename   ! Replace with your filename
      integer :: ncid, varid, ret
      real(real32), allocatable :: test_data(:)
      ! Try to open the file
      ret = nf90_open(filename, NF90_NOWRITE, ncid)
      if (ret /= nf90_noerr) then
         write(*,*) 'Error opening file:', trim(nf90_strerror(ret))
         stop
      end if

      ! Try to read zt specifically
      ret = nf90_inq_varid(ncid, "zt", varid)
      if (ret /= nf90_noerr) then
         write(*,*) 'Error getting zt variable:', trim(nf90_strerror(ret))
         stop
      end if

      ! Allocate and read
      allocate(test_data(20))  ! Since zt dimension is 20
      ret = nf90_get_var(ncid, varid, test_data)
      if (ret /= nf90_noerr) then
         write(*,*) 'Error reading zt data:', trim(nf90_strerror(ret))
      else
         write(*,*) 'Successfully read zt. First value:', test_data(1)
      end if

      ret = nf90_close(ncid)
   end subroutine do_test_nc
end module test_nc
