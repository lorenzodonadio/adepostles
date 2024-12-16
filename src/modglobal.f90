module modglobal
   use, intrinsic :: iso_fortran_env
   implicit none
   ! dimensions zm = zf, zt = zh in dales, so why not unify how it writes to nc. idk
   real(real32), allocatable :: zt(:), zm(:), xt(:), xm(:), yt(:), ym(:)
   real(real32), allocatable :: rtime(:)
   ! time loop variables
   real,parameter :: tres = 0.001 !< to convert simtime to seconds simtime*tres
   integer(int64) :: simtime = 0 !< simulation time in ms
   integer(int64) :: maxtime
   integer        :: dt = 50 !< Delta Time in ms
   integer        :: current_chunk = 0 !< from 0 to total_chunks, 0 means no fields have been read yet
   integer        :: total_chunks !< calculated as time_size/field_load_chunk_size, must be int otherwise program stops
   real(real32)   :: rsts = 0. !< real simulation time seconds
   real(real32)   :: next_chunk_load_time = -1. !< read immediatly please!
   real(real32)   :: next_save !<time for the next save

   !dimensions
   integer ::  imax, jmax, kmax, time_size
   integer ::  i1,i2,j1,j2,k1
   integer ::  ih = 2
   integer ::  jh = 2
   integer ::  kh = 1
   integer ::  nsv = 0       !< Number of scalar fields, loaded in modtracers

   !grid data
   logical :: luniformz = .false.
   real :: dx              !<  grid spacing in x-direction
   real :: dy              !<  grid spacing in y-direction
   real :: dz              !<  grid spacing in z-direction
   real :: dxi             !<  1/dx
   real :: dyi             !<  1/dy
   real :: dzi             !<  1/dz
   ! real :: dxiq            !<  1/(dx*4)
   ! real :: dyiq            !<  1/(dy*4)
   ! real :: dziq            !<  1/(dz*4)
   real :: dxi5            !<  1/(2*dx)
   real :: dyi5            !<  1/(2*dy)
   real :: dzi5            !<  1/(2*dz)
   real :: dx2i            !<  (1/dx)**2
   real :: dy2i            !<  (1/dy)**2
   real :: dz2i      !< (1/dz)**2
   real, allocatable :: dzf(:)         !<  thickness of full level
   real, allocatable :: dzh(:)         !<  thickness of half level
   real, allocatable :: delta(:)       !<  (dx*dy*dz)**(1/3)
   real, allocatable :: deltai(:)       !<  (dx*dy*dz)**(-1/3)  or dzf**-1 for anisotropic diffusion

contains
   subroutine init_global()
      use config, only: runtime,dtmax,output_save_interval
      integer :: k

      if (rtime(time_size)<runtime) then
         runtime = int(rtime(time_size))
         write(*,*) "selected runtime exceeded time span of provided data, setting runtime to: ", runtime, " s"
      endif
      maxtime = int(runtime/tres)

      dt = min(int(dtmax/tres),dt) !dtmax comes from namoptions and its in seconds

      write (*,*) "init dt: ", dt, "ms"
      next_save = output_save_interval

      i1 = imax +1
      j1 = jmax +1
      i2 = imax +2
      j2 = jmax +2

      k1 = kmax +1

      dx = xm(2)-xm(1)
      dy = ym(2)-ym(1)
      dxi = 1./dx
      dyi = 1./dy

      write(*,*) "dx,dy,dxi,dyi", dx,dy,dxi,dyi

      dxi5 = 0.5*dxi
      dyi5 = 0.5*dyi

      dx2i = dxi**2
      dy2i = dyi**2

      allocate(dzf(k1),dzh(k1),delta(k1),deltai(k1))

      ! zm = zf, zt = zh
      do  k=1,kmax-1
         dzf(k) = zm(k+1) - zm(k)
      end do
      dzf(kmax) = dzf(kmax-1) ! outofbounds error if not
      dzf(k1) = dzf(kmax)

      dzh(1) = 2*zt(1)
      do k=2,kmax
         dzh(k) = zt(k) - zt(k-1)
      end do

      dzh(k1) = dzh(kmax) !avoid outofbounds error

      luniformz = all(dzh .eq. dzh(1))

      if (luniformz) then
         dz = dzh(1)
         dzi = 1./dz
         dz2i = dzi**2
         dzi5 = 0.5*dzi

         write(*,*) 'Working with uniform z grid:'
      endif
      ! dzh(1) = 2*zm(1)
      ! do k=2,kmax
      !    dzh(k) = zm(k) - zm(k-1)
      ! end do

      ! dzh(k1) = dzh(kmax) !avoid outofbounds error

      do k=1,k1
         delta(k) = (dx*dy*dzf(k))**(1./3.)
         deltai(k) = 1./delta(k)     !can be overruled in modsubgrid in case anisotropic diffusion is applied
      end do
   end subroutine init_global

   subroutine load_dimensions()
      use netcdf_utils, only :nchandle_error, get_dimension_size, get_1d_variable
      use config, only: field_dump_path
      use netcdf

      integer :: ncid, retval

      ! Dimension IDs and sizes
      integer :: zt_dim, zm_dim, xt_dim, xm_dim, yt_dim, ym_dim,time_dim
      integer :: zt_size, zm_size, xt_size, xm_size, yt_size, ym_size

      ! Variable IDs
      integer :: zt_varid, zm_varid, xt_varid, xm_varid, yt_varid, ym_varid,time_varid

      ! Static variable arrays (dynamic allocation)

      ! Open the NetCDF file
      retval = nf90_open(field_dump_path, NF90_NOWRITE, ncid)
      call nchandle_error(retval, 'Error opening file: '//trim(field_dump_path))

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
      call get_1d_variable(ncid, 'zt', zt_varid, zt)
      call get_1d_variable(ncid, 'zm', zm_varid, zm)
      call get_1d_variable(ncid, 'xt', xt_varid, xt)
      call get_1d_variable(ncid, 'xm', xm_varid, xm)
      call get_1d_variable(ncid, 'yt', yt_varid, yt)
      call get_1d_variable(ncid, 'ym', ym_varid, ym)
      call get_1d_variable(ncid, 'time', time_varid, rtime)

      ! Close the NetCDF file
      retval = nf90_close(ncid)
      call nchandle_error(retval, 'Error closing file: '//trim(field_dump_path))

      print *, 'Dimensions loaded successfully from: ', trim(field_dump_path)

   end subroutine load_dimensions

   subroutine validate_time_size()
      use config, only: field_load_chunk_size

      if (mod(size(rtime), field_load_chunk_size) /= 0) then
         print *, 'time_size:', size(rtime)
         print *, 'field_load_chunk_size :', field_load_chunk_size
         stop 'Time size must be a multiple of Field load chunk, please adjust the options'
      endif

      total_chunks = size(rtime)/field_load_chunk_size
   end subroutine validate_time_size

end module modglobal
