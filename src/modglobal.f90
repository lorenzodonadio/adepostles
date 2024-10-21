module modglobal
   use, intrinsic :: iso_fortran_env
   use netcdf
   use netcdf_loader, only :nchandle_error, get_dimension_size, get_and_read_variable
   use config, only: field_load_chunk_size
   use modfields, only: u,v,w,ekh
   implicit none
   integer :: total_chunks
   ! dimensions zm = zf, zt = zh in dales, so why not unify how it writes to nc. idk
   real(real32), allocatable :: zt(:), zm(:), xt(:), xm(:), yt(:), ym(:),rtime(:)
   ! fields

   !dimensions
   integer ::  imax, jmax, kmax, time_size
   integer ::  i1,i2,ih,j1,j2,jh,k1

   real :: dx              !<  grid spacing in x-direction
   real :: dy              !<  grid spacing in y-direction
   ! real :: dz              !<  grid spacing in z-direction
   real :: dxi             !<  1/dx
   real :: dyi             !<  1/dy
   ! real :: dzi             !<  1/dz
   ! real :: dxiq            !<  1/(dx*4)
   ! real :: dyiq            !<  1/(dy*4)
   ! real :: dziq            !<  1/(dz*4)
   real :: dxi5            !<  1/(2*dx)
   real :: dyi5            !<  1/(2*dy)
   ! real :: dzi5            !<  1/(2*dz)
   real :: dx2i            !<  (1/dx)**2
   real :: dy2i            !<  (1/dy)**2

   real, allocatable :: dzf(:)         !<  thickness of full level
   real, allocatable :: dzh(:)         !<  thickness of half level
   real, allocatable :: delta(:)       !<  (dx*dy*dz)**(1/3)
   real, allocatable :: deltai(:)       !<  (dx*dy*dz)**(-1/3)  or dzf**-1 for anisotropic diffusion

contains
   subroutine init_global()
      i1 = imax +1
      j1 = jmax +1
      i2 = imax +2
      j2 = jmax +2

      k1 = kmax +1

      dx = xm(2)-xm(1)
      dy = ym(2)-ym(1)

      dxi = 1/dx
      dyi = 1/dy

      dxi5 = 0.5*dxi
      dyi5 = 0.5*dyi

      dx2i = dxi**2
      dy2i = dyi**2

      allocate(dzf(k1),dzh(k1),delta(k1),deltai(k1))
      ! zm = zf, zt = zh
      do  k=1,kmax
         dzf(k) = zt(k+1) - zt(k)
      end do
      dzf(k1) = dzf(kmax)

      dzh(1) = 2*zm(1)
      do k=2,k1
         dzh(k) = zm(k) - zm(k-1)
      end do

      do k=1,k1

         delta(k) = (dx*dy*dzf(k))**(1./3.)
         deltai(k) = 1./delta(k)     !can be overruled in modsubgrid in case anisotropic diffusion is applied
      end do
   end subroutine init_global

   subroutine load_dimensions(filename)
      character(len=*), intent(in) :: filename
      integer :: ncid, retval

      ! Dimension IDs and sizes
      integer :: zt_dim, zm_dim, xt_dim, xm_dim, yt_dim, ym_dim,time_dim
      integer :: zt_size, zm_size, xt_size, xm_size, yt_size, ym_size

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

      if (xt_size /= xm_size) stop 'xt and xm dimensions must have the same size'
      if (yt_size /= ym_size) stop 'yt and ym dimensions must have the same size'
      if (zt_size /= zm_size) stop 'zt and zm dimensions must have the same size'
      ! Set global vars exported by module
      imax = xt_size
      jmax = yt_size
      kmax = zt_size
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

      print *, 'Dimensions loaded successfully from: ', trim(filename)

   end subroutine load_dimensions

   subroutine load_fields_chunk(filename, chunk_number)
      character(len=*), intent(in) :: filename
      integer, intent(in) :: chunk_number
      integer :: ncid, retval


      ! Variable IDs for time-dependent variables
      integer :: u_varid, v_varid, w_varid, ekh_varid

      ! Time-dependent variable arrays (dynamic allocation)

      ! Open the NetCDF file
      retval = nf90_open(filename, NF90_NOWRITE, ncid)
      call nchandle_error(retval, 'Error opening file: '//trim(filename))

      ! Read time-dependent variables from the file
      call get_and_read_variable_chunk(ncid, 'u', u_varid, u, chunk_number)
      call get_and_read_variable_chunk(ncid, 'v', v_varid, v, chunk_number)
      call get_and_read_variable_chunk(ncid, 'w', w_varid, w, chunk_number)
      call get_and_read_variable_chunk(ncid, 'ekh', ekh_varid, ekh, chunk_number)

      ! Close the NetCDF file
      retval = nf90_close(ncid)
      call nchandle_error(retval, 'Error closing file: '//trim(filename))

      ! print *, 'Fields variables loaded successfully from ', trim(filename), 'for chunk: ', chunk_number,'/',total_chunks
      print *, 'Fields loaded - chunk: ', chunk_number,'/',total_chunks
   end subroutine load_fields_chunk


   subroutine validate_time_size()
      if (mod(size(rtime), field_load_chunk_size) /= 0) then
         print *, 'time_size:', size(rtime)
         print *, 'field_load_chunk_size :', field_load_chunk_size
         stop 'Time size must be a multiple of Field load chunk, please adjust the options'
      endif

      total_chunks = size(rtime)/field_load_chunk_size
   end subroutine validate_time_size

end module modglobal
