module modfields
   use iso_fortran_env
   implicit none

   integer :: nti = 0 !< next time index, index of closest fields measurement in future time
   real :: dtfield !< delta t between the two currecly active fields
   !this makes life easier
   real(real32), allocatable :: chunktime(:) !< a time vector similar to rtime but that only holds the times for the loaded chunk of fields

   real(real32), allocatable :: rhobf_chunk(:,:)   !< density full level (zf,time)
   real(real32), allocatable :: rhobh_chunk(:,:)   !< density half level (zh,time)
   real(real32), allocatable :: rhobf(:)   !< density full level (zf)
   real(real32), allocatable :: rhobfi(:)  !< inverse density full level (zf)
   real(real32), allocatable :: rhobh(:)   !< density half level (zh)
   integer,allocatable :: ksfc(:,:) !< Kmin of the surface level, for IBM, default is 1, allocated in modfields but ibm_init fills it according to the surface heigh

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

contains
   subroutine interpolate_fields_to_simtime()
      !! populates u0,v0,w0,and ek0 to match simtime
      use config, only: field_load_chunk_size
      use modglobal, only: rsts
      real  :: d !< distance to the next point in -1,0 space
      !this happens only once per simulation, maybe solve with an initial field, DALES does not save the sim at time = 0

      ! if the simtime is greater than the time at that index, increment the index
      ! the index can not be greater than field_load_chunk_size, this is done by resetting to 1 every load fields
      if (rsts > chunktime(nti)) then
         nti = nti+1

         um = up
         vm = vp
         wm = wp
         ekhm = ekhp

         up = u(:,:,:,nti)
         vp = v(:,:,:,nti)
         wp = w(:,:,:,nti)
         ekhp = ekh(:,:,:,nti)

         dtfield = chunktime(nti)-chunktime(nti-1)
      endif

      ! for now simple linear interpolation, TODO improve this
      ! the time from simulation to the next point divided by dt of the field
      d = (chunktime(nti) - rsts) / dtfield
      u0 = up*(1-d)+um*d
      v0 = vp*(1-d)+vm*d
      w0 = wp*(1-d)+wm*d
      ekh0 = ekhp*(1-d)+ekhm*d
      !TODO make a better handling profiles interpolation at nti ==1
      if (nti == 1) then
         ! to calculate the distance use a combination of old and new dtfield, only needed if irregulat timestep and even so i doubt it
         ! d = (chunktime(nti) - rsts)/(0.5*dtfield+0.5*(chunktime(nti+1)-chunktime(nti)))
         rhobf = rhobf_chunk(:,1)
         rhobh = rhobh_chunk(:,1)
      else
         rhobf = rhobf_chunk(:,nti)*(1-d)+rhobf_chunk(:,nti-1)*d
         rhobh = rhobh_chunk(:,nti)*(1-d)+rhobh_chunk(:,nti-1)*d
      endif

      rhobfi = 1./rhobf
      ! write(*,*) 'simtime: ', rsts,' next time: ', chunktime(nti), ' delta: ',chunktime(nti) - rsts
   end subroutine interpolate_fields_to_simtime

   subroutine lerp_field_backwards(f_in,d,tidx, f_out)
      !f(t) ~ f(0) + (dtsim/dtfield)*( f(-1) - f(0)) for interval -1,0, or nti-1 to nti
      real, intent(in) :: f_in(:,:,:,:)
      real, intent(in) :: d !< distance "
      integer, intent(in) :: tidx
      real, intent(out) ::  f_out(:,:,:)
      f_out = f_in(:,:,:,tidx)*(1-d)+f_in(:,:,:,tidx-1)*d

   end subroutine lerp_field_backwards

   subroutine allocate_fields()
      !< Allocates and sets to zero, all fields and profiles
      use modglobal, only:i1,ih,j1,jh,k1
      use config, only: field_load_chunk_size

      allocate(chunktime(field_load_chunk_size))

      allocate(rhobf_chunk(k1,field_load_chunk_size))
      allocate(rhobh_chunk(k1,field_load_chunk_size))

      allocate(rhobf(k1))
      allocate(rhobfi(k1))
      allocate(rhobh(k1))
      allocate(ksfc     (2-ih:i1+ih,2-jh:j1+jh))
      ksfc = 1 !this is the default value, if not IBM, then all k-loops are executed from k=ksfc=1

      allocate(ekh  (2-ih:i1+ih,2-jh:j1+jh,k1,field_load_chunk_size))
      allocate(u    (2-ih:i1+ih,2-jh:j1+jh,k1,field_load_chunk_size))
      allocate(v    (2-ih:i1+ih,2-jh:j1+jh,k1,field_load_chunk_size))
      allocate(w    (2-ih:i1+ih,2-jh:j1+jh,k1,field_load_chunk_size))

      allocate(ekh0 (2-ih:i1+ih,2-jh:j1+jh,k1))
      allocate(u0   (2-ih:i1+ih,2-jh:j1+jh,k1))
      allocate(v0   (2-ih:i1+ih,2-jh:j1+jh,k1))
      allocate(w0   (2-ih:i1+ih,2-jh:j1+jh,k1))

      allocate(ekhm (2-ih:i1+ih,2-jh:j1+jh,k1))
      allocate(um   (2-ih:i1+ih,2-jh:j1+jh,k1))
      allocate(vm   (2-ih:i1+ih,2-jh:j1+jh,k1))
      allocate(wm   (2-ih:i1+ih,2-jh:j1+jh,k1))

      allocate(ekhp (2-ih:i1+ih,2-jh:j1+jh,k1))
      allocate(up   (2-ih:i1+ih,2-jh:j1+jh,k1))
      allocate(vp   (2-ih:i1+ih,2-jh:j1+jh,k1))
      allocate(wp   (2-ih:i1+ih,2-jh:j1+jh,k1))

      !   c0 = 0.
      !   cp = 0.
      !   cm = 0.

      ! c0(10:12,20:22,3:5) = 5.
      ! c0(40:44,70:72,3:5) = 5.
      ! c0(10:12,5:100,2:3) = 5.
      ! c0(50:52,60:100,3:4) = 5.
   end subroutine allocate_fields

   subroutine init_interp_fields()
      use config, only: field_dump_path,field_load_chunk_size
      use modglobal, only :rtime
      call load_fields_chunk(field_dump_path,1)

      up = u(:,:,:,1)
      vp = v(:,:,:,1)
      wp = w(:,:,:,1)
      ekhp = ekh(:,:,:,1)

      um = up
      vm = vp
      wm = wp
      ekhm = ekhp

      dtfield = rtime(1)
   end subroutine init_interp_fields

   subroutine load_fields_intimeloop()
      use modglobal, only: current_chunk, next_chunk_load_time, rtime, rsts
      use config, only: field_dump_path,field_load_chunk_size

      if (rsts > next_chunk_load_time) then
         current_chunk = current_chunk + 1
         next_chunk_load_time = rtime(current_chunk*field_load_chunk_size)
         call load_fields_chunk(field_dump_path,current_chunk)

         chunktime  = rtime(1+(current_chunk-1)*field_load_chunk_size:current_chunk*field_load_chunk_size)

         nti = 1 ! reset the index here every loaded chunk
      endif

   end subroutine load_fields_intimeloop

!    subroutine load_fields_chunk(filename, chunk_number)
!       use netcdf
!       use mpi_f08
!       use modglobal, only: i1, j1, k1
!       use netcdf_utils, only: nchandle_error, get_field_chunk, get_profile_chunk
!       use modglobal, only: total_chunks
!       implicit none

!       integer :: dimids(4), dim_len
!       integer :: i, ndims, dimsizes(4)

!       character(len=*), intent(in) :: filename
!       integer, intent(in) :: chunk_number
!       integer :: ncid, retval, varid, my_id, comm_size

!       ! Get MPI rank and communicator size
!       call MPI_Comm_rank(MPI_COMM_WORLD, my_id, retval)
!       call MPI_Comm_size(MPI_COMM_WORLD, comm_size, retval)

!       if (my_id == 0) then
!          ! Only root process loads data from file
!          if (chunk_number > total_chunks .or. chunk_number < 1) then
!             write(*,*) 'Error loading chunks, invalid chunk number, exiting'
!             write(*,*) 'Chunk number: ', chunk_number, " - Total chunks: ", total_chunks
!             return
!          endif

!          ! Open the NetCDF file
!          retval = nf90_open(filename, NF90_NOWRITE, ncid)
!          call nchandle_error(retval, 'Error opening file: '//trim(filename))
!          ! Get the variable ID for 'u'
!          retval = nf90_inq_varid(ncid, 'u', varid)
!          call nchandle_error(retval, 'Error: Could not get variable ID for u')
!          call nchandle_error(retval, 'Error: Could not inquire variable u')

!          ! Inquire variable dimensions and print them
!          retval = nf90_inquire_variable(ncid, varid, ndims=ndims, dimids=dimids)

!          ! Loop through each dimension and get its size
!          do i = 1, ndims
!             retval = nf90_inquire_dimension(ncid, dimids(i), len=dim_len)
!             call nchandle_error(retval, 'Error: Could not inquire dimension size')
!             dimsizes(i) = dim_len
!          end do


!          ! Print the shape of 'u' as stored in the NetCDF file
!          print *, 'Expected shape of u in NetCDF file: (', dimsizes(1:ndims), ')'
!          print *, 'Expected shape U', shape(u(2:i1,2:j1,2:k1,:)), ')'
!          ! Load fields
!          call get_field_chunk(ncid, 'u',   varid, u  (2:i1,2:j1,2:k1,:), chunk_number)
!          call get_field_chunk(ncid, 'v',   varid, v  (2:i1,2:j1,2:k1,:), chunk_number)
!          call get_field_chunk(ncid, 'w',   varid, w  (2:i1,2:j1,2:k1,:), chunk_number)
!          call get_field_chunk(ncid, 'ekh', varid, ekh(2:i1,2:j1,2:k1,:), chunk_number)

!          ! Load profiles
!          call get_profile_chunk(ncid, 'rhobh', varid, rhobh_chunk, chunk_number)
!          call get_profile_chunk(ncid, 'rhobf', varid, rhobf_chunk, chunk_number)

!          ! Close the NetCDF file
!          retval = nf90_close(ncid)
!          call nchandle_error(retval, 'Error closing file: '//trim(filename))

!          ! Pad the fields before broadcasting (only on root)
!          call pad_field(u)
!          call pad_field(v)
!          call pad_field(w)
!          call pad_field(ekh)
!          call pad_profile(rhobf_chunk)
!          call pad_profile(rhobh_chunk)
!       endif

!       call MPI_Barrier(MPI_COMM_WORLD)
!       ! Broadcast the data to all other processes
!       call MPI_Bcast(u, size(u), MPI_REAL, 0, MPI_COMM_WORLD, retval)
!       call MPI_Bcast(v, size(v), MPI_REAL, 0, MPI_COMM_WORLD, retval)
!       call MPI_Bcast(w, size(w), MPI_REAL, 0, MPI_COMM_WORLD, retval)
!       call MPI_Bcast(ekh, size(ekh), MPI_REAL, 0, MPI_COMM_WORLD, retval)
!       call MPI_Bcast(rhobf_chunk, size(rhobf_chunk), MPI_REAL, 0, MPI_COMM_WORLD, retval)
!       call MPI_Bcast(rhobh_chunk, size(rhobh_chunk), MPI_REAL, 0, MPI_COMM_WORLD, retval)

!       ! Check for errors in MPI_Bcast
!       if (retval /= MPI_SUCCESS) then
!          print *, 'Error broadcasting data from root'
!          call MPI_Abort(MPI_COMM_WORLD, retval)
!       endif

!       print *, 'Loaded fields chunk: ', chunk_number, '/', total_chunks, ' on process ', my_id
!    end subroutine load_fields_chunk

   subroutine load_fields_chunk(filename, chunk_number)
      !< Loads fields: u,v,w,ekh (3d+time) and vertical profiles rhobf & rhobh (1d+time) for the specific time specified by the chunk_number
      use netcdf
      use mpi_f08
      use modglobal, only:i1,j1,kmax
      use netcdf_utils, only :nchandle_error, get_field_chunk,get_profile_chunk
      use modglobal,only: total_chunks
      use modmpi, only: is_root
      character(len=*), intent(in) :: filename
      integer, intent(in) :: chunk_number
      integer :: ncid, retval, varid

      if (is_root) then
         if (chunk_number > total_chunks .or. chunk_number < 1) then
            write(*,*) 'Error loading chunks, invalid chunk number, exiting'
            write(*,*) 'Chunk number: ',chunk_number, " - Total chunks: ", total_chunks
            return
         endif

         ! Time-dependent variable arrays (dynamic allocation)

         ! Open the NetCDF file
         retval = nf90_open(filename, NF90_NOWRITE, ncid)
         call nchandle_error(retval, 'Error opening file: '//trim(filename))


         ! Read time-dependent variables from the file
         ! kmax instead of k1 because of number of ghost cells
         ! For debuggind purposes,TODO implement a debug flag?
         ! print *, shape(u)
         ! print *, shape(u(2:i1,2:j1,2:kmax,:))
         ! print *, 'Vertical U slice: u(5,5,:,1)'
         ! print *, u(5,5,:,1)
         ! print *, 'Horizontal U slice'
         ! print *, u(:,5,5,1)
         ! Load fields
         call get_field_chunk(ncid, 'u',   varid, u  (2:i1,2:j1,2:kmax,:), chunk_number)
         call get_field_chunk(ncid, 'v',   varid, v  (2:i1,2:j1,2:kmax,:), chunk_number)
         call get_field_chunk(ncid, 'w',   varid, w  (2:i1,2:j1,2:kmax,:), chunk_number)
         call get_field_chunk(ncid, 'ekh', varid, ekh(2:i1,2:j1,2:kmax,:), chunk_number)
         ! Load profiles
         call get_profile_chunk(ncid, 'rhobh', varid, rhobh_chunk(2:kmax,:), chunk_number)
         call get_profile_chunk(ncid, 'rhobf', varid, rhobf_chunk(2:kmax,:), chunk_number)

         ! Close the NetCDF file
         retval = nf90_close(ncid)
         call nchandle_error(retval, 'Error closing file: '//trim(filename))

         ! print *, 'Fields variables loaded successfully from ', trim(filename), 'for chunk: ', chunk_number,'/',total_chunks
         call pad_field(u)
         call pad_field(v)
         call pad_field(w)
         call pad_field(ekh)
         call pad_profile(rhobf_chunk)
         call pad_profile(rhobh_chunk)

         print *, 'Loaded fields chunk: ', chunk_number,'/',total_chunks
      endif

      call MPI_Barrier(MPI_COMM_WORLD)
      ! Broadcast the data to all other processes
      call MPI_Bcast(u, size(u), MPI_REAL, 0, MPI_COMM_WORLD, retval)
      call MPI_Bcast(v, size(v), MPI_REAL, 0, MPI_COMM_WORLD, retval)
      call MPI_Bcast(w, size(w), MPI_REAL, 0, MPI_COMM_WORLD, retval)
      call MPI_Bcast(ekh, size(ekh), MPI_REAL, 0, MPI_COMM_WORLD, retval)
      call MPI_Bcast(rhobf_chunk, size(rhobf_chunk), MPI_REAL, 0, MPI_COMM_WORLD, retval)
      call MPI_Bcast(rhobh_chunk, size(rhobh_chunk), MPI_REAL, 0, MPI_COMM_WORLD, retval)

      ! Check for errors in MPI_Bcast
      if (retval /= MPI_SUCCESS) then
         print *, 'Error broadcasting data from root'
         call MPI_Abort(MPI_COMM_WORLD, retval)
      endif

   end subroutine load_fields_chunk


   subroutine pad_profile(p)
      use modglobal,only:k1,kh,kmax,dzh
      real(real32), intent(inout) :: p(:,:) !< some profile (z,time)
      p(1:kh,:) = (1.+1./dzh(1))*p(1+kh:2*kh,:)-(1/dzh(1))*p(2+kh:1+2*kh,:)
      p(k1-kh:k1,:) = (1.+1./dzh(kmax))*p(k1-2*kh:k1-kh,:) - (1/dzh(kmax))*p(kmax-2*kh:kmax-kh,:)
   end subroutine pad_profile

   subroutine pad_field(f)
      !! fills ghost cells form fields, reconstructing boundary conditions
      !! Options: 1. Periodic BC (default), 2. Neumann BC (copying points), 3. Second derivative to Zero (interpolation)
      !! f: some 4d dimensional field (x,y,z,time)
      use modglobal,only:i1,ih,j1,jh,k1,kh,imax, jmax, kmax,dzh,dxi,dyi
      use config, only: lperiodic_field_pad
      real(real32), intent(inout) :: f(:,:,:,:)
      !! treat horizontal (X,y) boundaries
      !! periodic bc
      if (lperiodic_field_pad) then
         !West from East
         f(1:ih,:,:,:) = f(i1:imax+ih,:,:,:)
         !East from west
         f(i1+ih:imax+2*ih,:,:,:) = f(1+ih:2*ih,:,:,:)
         !South form north
         f(:,1:jh,:,:) = f(:,j1:jmax+jh,:,:)
         !North from south
         f(:,j1+jh:jmax+2*jh,:,:) = f(:,1+jh:2*jh,:,:)
      else
         !! lerp 2nd order aprox of derivative
         !! X direction interp
         !! f(-1) ~~ (1+1.5/dx)*f0-(2/dx)*f(1)+(1/2dx)*f(2)
         f(ih,:,:,:) = (1+1.5*dxi)*f(ih+1,:,:,:) -2*dxi*f(ih+2,:,:,:)+0.5*dxi*f(ih+3,:,:,:)
         f(ih-1,:,:,:) = (1+1.5*dxi)*f(ih,:,:,:) -2*dxi*f(ih+1,:,:,:)+0.5*dxi*f(ih+2,:,:,:)
         !! f(1) ~~ (1+1.5/dx)*f0-(2/dx)*f(-1)+(1/2dx)*f(-2)
         f(i1+ih,:,:,:) = (1+1.5*dxi)*f(i1+ih-1,:,:,:) -2*dxi*f(i1+ih-2,:,:,:)+0.5*dxi*f(i1+ih-3,:,:,:)
         f(i1+ih+1,:,:,:) = (1+1.5*dxi)*f(i1+ih,:,:,:) -2*dxi*f(i1+ih-1,:,:,:)+0.5*dxi*f(i1+ih-2,:,:,:)
         !! Y direction interp
         f(:,jh,:,:) = (1+1.5*dyi)*f(:,jh+1,:,:) -2*dyi*f(:,jh+2,:,:)+0.5*dyi*f(:,jh+3,:,:)
         f(:,jh-1,:,:) = (1+1.5*dyi)*f(:,jh,:,:) -2*dyi*f(:,jh+1,:,:)+0.5*dyi*f(:,jh+2,:,:)
         !! f(1) ~~ (1+1.5/dx)*f0-(2/dx)*f(-1)+(1/2dx)*f(-2)
         f(:,j1+jh,:,:) = (1+1.5*dyi)*f(:,j1+jh-1,:,:) -2*dyi*f(:,j1+jh-2,:,:)+0.5*dyi*f(:,j1+jh-3,:,:)
         f(:,j1+jh+1,:,:) = (1+1.5*dyi)*f(:,j1+jh,:,:) -2*dyi*f(:,j1+jh-2,:,:)+0.5*dyi*f(:,j1+jh-2,:,:)
      endif

      !! Z field is always just extended (maybe better to lerp it?) yes
      !! extend the fields
      ! f(:,:,1:kh,:) = f(:,:,1+kh:2*kh,:)
      ! f(:,:,k1-kh:k1,:) = f(:,:,k1-2*kh:k1-kh,:)
      !! lerp
      f(:,:,1:kh,:) = (1.+1./dzh(1))*f(:,:,1+kh:2*kh,:)-(1/dzh(1))*f(:,:,2+kh:1+2*kh,:)
      f(:,:,k1-kh:k1,:) = (1.+1./dzh(kmax))*f(:,:,k1-2*kh:k1-kh,:) - (1/dzh(kmax))*f(:,:,kmax-2*kh:kmax-kh,:)

   end subroutine pad_field


end module modfields


