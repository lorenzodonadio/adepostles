module modfields
   use iso_fortran_env
   implicit none

   real(real32), allocatable :: ekh(:,:,:,:) !< k-coefficient for eddy diffusivity ekh(xt,yt,zt,time) m2/s
   real(real32), allocatable :: u(:,:,:,:)   !< u(xm,yt,zt,time) m/s
   real(real32), allocatable :: v(:,:,:,:)   !< v(xt,ym,zt,time) m/s
   real(real32), allocatable :: w(:,:,:,:)   !< w(xt,yt,zm,time) m/s
   ! at different time stamps
   real(real32), allocatable :: c0(:,:,:)   !< Concentration
   real(real32), allocatable :: cp(:,:,:)   !< Concentration
   real(real32), allocatable :: cm(:,:,:)   !< Concentration

   real(real32), allocatable :: ekh0(:,:,:) !< ekhm at time t
   real(real32), allocatable :: u0(:,:,:)   !< u at time t m/s
   real(real32), allocatable :: v0(:,:,:)   !< v at time t m/s
   real(real32), allocatable :: w0(:,:,:)   !< w at time t m/s

   real(real32), allocatable :: ekhp(:,:,:) !< ekhm at time t+1
   real(real32), allocatable :: up(:,:,:)   !< u at time t+1 m/s
   real(real32), allocatable :: vp(:,:,:)   !< v at time t+1 m/s
   real(real32), allocatable :: wp(:,:,:)   !< w at time t+1 m/s

   real(real32), allocatable :: ekhm(:,:,:) !< ekhm at time t-1
   real(real32), allocatable :: um(:,:,:)   !< u at time t-1m/s
   real(real32), allocatable :: vm(:,:,:)   !< v at time t-1 m/s
   real(real32), allocatable :: wm(:,:,:)   !< w at time t-1m/s

   real(real32), allocatable :: rhobf(:)   !< density full level
   real(real32), allocatable :: rhobh(:)   !< density half level

contains

   subroutine allocate_fields()
      use modglobal, only:i1,ih,j1,jh,k1
      use config, only: field_load_chunk_size

      allocate(rhobf(k1))
      allocate(rhobh(k1))

      allocate(ekh  (2-ih:i1+ih,2-jh:j1+jh,k1,field_load_chunk_size))
      allocate(u    (2-ih:i1+ih,2-jh:j1+jh,k1,field_load_chunk_size))
      allocate(v    (2-ih:i1+ih,2-jh:j1+jh,k1,field_load_chunk_size))
      allocate(w    (2-ih:i1+ih,2-jh:j1+jh,k1,field_load_chunk_size))

      ! TODO make this parallelizable, so that we can handle multiple ADE solutions with the same u,v,w,ekh,rho, thats the goal
      ! maybe we just add a fourth dimension and thats it?
      allocate(c0   (2-ih:i1+ih,2-jh:j1+jh,k1))
      allocate(cp   (2-ih:i1+ih,2-jh:j1+jh,k1))
      allocate(cm   (2-ih:i1+ih,2-jh:j1+jh,k1))

      allocate(ekh0 (2-ih:i1+ih,2-jh:j1+jh,k1))
      allocate(u0   (2-ih:i1+ih,2-jh:j1+jh,k1))
      allocate(v0   (2-ih:i1+ih,2-jh:j1+jh,k1))
      allocate(w0   (2-ih:i1+ih,2-jh:j1+jh,k1))

      allocate(ekhp (2-ih:i1+ih,2-jh:j1+jh,k1))
      allocate(up   (2-ih:i1+ih,2-jh:j1+jh,k1))
      allocate(vp   (2-ih:i1+ih,2-jh:j1+jh,k1))
      allocate(wp   (2-ih:i1+ih,2-jh:j1+jh,k1))

      allocate(ekhm (2-ih:i1+ih,2-jh:j1+jh,k1))
      allocate(um   (2-ih:i1+ih,2-jh:j1+jh,k1))
      allocate(vm   (2-ih:i1+ih,2-jh:j1+jh,k1))
      allocate(wm   (2-ih:i1+ih,2-jh:j1+jh,k1))

   end subroutine allocate_fields


   subroutine load_fields_chunk(filename, chunk_number)
      use netcdf
      use netcdf_loader, only :nchandle_error, get_and_read_variable_chunk
      use modglobal,only: total_chunks
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

end module modfields
