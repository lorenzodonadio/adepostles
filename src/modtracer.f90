module modtracer
   use iso_fortran_env, only: real32
   use modtracer_type, only: T_tracer_collection
   use modglobal, only: nsv
   use netcdf
   implicit none

   integer :: output_ncid !< where to write the output concentration

   !    integer :: x_dim, y_dim, z_dim, time_dim  ! Dimension IDs
   !    integer :: x_size, y_size, z_size, time_size
   !    character(len=80) :: group_name

   ! Arrays for storing dimensions and variable data
   type(T_tracer_collection) :: tracers
   real, allocatable :: cm(:,:,:,:)
   real, allocatable :: c0(:,:,:,:)
   real, allocatable :: cp(:,:,:,:)
   real, allocatable :: source(:,:)
   integer, allocatable :: source_meta(:,:) !> contains i,j,k,tracernumber,sourcenumber,idx, corresponding to the source row that holds the actual time series data
   integer :: src_t_idx = 1
   
contains
   !       allocate(c0   (2-ih:i1+ih,2-jh:j1+jh,k1))
   ! allocate(cp   (2-ih:i1+ih,2-jh:j1+jh,k1))
   ! allocate(cm   (2-ih:i1+ih,2-jh:j1+jh,k1))
   subroutine apply_source()
      use modglobal, only: rtime,rsts
      integer :: nsrc,i,j,k,sv

      if (rsts > rtime(src_t_idx)) src_t_idx = src_t_idx + 1

      do nsrc = 1, size(source_meta,1)
         i = source_meta(nsrc,1)
         j = source_meta(nsrc,2)
         k = source_meta(nsrc,3)
         sv = source_meta(nsrc,4)
         ! timerow = sou
         c0(i,j,k,sv) = source(nsrc,src_t_idx)

      end do
      ! if (rkstep == 1) then
      ! c0(48:50,80:82,3) = 1
      !   c0(51,100:104,3,2) = 1
      !   c0(51,62:66,3,1) = 1
      !   c0(51,4:8,3,1) = 1
      ! endif

   end subroutine apply_source

   subroutine load_tracer_init_and_sources()
      use netcdf_utils, only: nchandle_error
      use modutils, only: str_ends_with
      use modglobal, only: i1,ih,j1,jh,k1,time_size,nsv
      use config, only: sources_prefix
      use mpi_f08
      integer :: varid, src_varid, numsources,totnumsources, retval
      integer :: i, itsv, itsrc, src_counter                ! Loop counters
      integer :: src_x,src_y,src_z     !grid
      integer :: ncid, numgrps,groupid !netcdf stuff
      integer :: my_id !mpi stuff
      integer,allocatable :: grpids(:)
      character(len=80) :: varname
      character(len=256) :: sources_path
      character(len=256) :: group_name
      !   character(len=256) ::desc
      !   character(len=256) ::hist
      !   DEBUG VARS
      !   integer :: dimids(3), dim_size,i
      !   character(len=16) ::dimname

      ! Get the MPI rank (ID) of the current process
      call MPI_Comm_rank(MPI_COMM_WORLD, my_id, retval)
      if (retval /= MPI_SUCCESS) then
         print *, 'Error: Could not get MPI rank'
         call MPI_Abort(MPI_COMM_WORLD, retval)
      endif

      ! Format the sources_path based on my_id
      sources_path = trim(sources_prefix)

      if (.not.str_ends_with(sources_path,'.nc')) then
         write(sources_path, '(A,"_",I3.3,".nc")') trim(sources_prefix), my_id + 1
      endif

      print *,'Sources File: ', sources_path
      retval = nf90_open(sources_path, nf90_nowrite, ncid)
      call nchandle_error(retval, 'Error: Could not open NetCDF file'//sources_path)
      retval = nf90_get_att(ncid, nf90_global, 'numtracers', nsv)
      !   retval = nf90_get_att(ncid, nf90_global, 'decription', desc)
      !   retval = nf90_get_att(ncid, nf90_global, 'history', hist)
      call nchandle_error(retval, 'Error: Could not read global attribute "numtracers"')
      allocate(grpids(nsv))
      retval = nf90_get_att(ncid, nf90_global, 'totnumsources', totnumsources)
      call nchandle_error(retval, 'Error: Could not read global attribute "totnumsources"')
      retval = nf90_inq_grps(ncid, numgrps, grpids)
      call nchandle_error(retval, 'Error: Could not Inquire groups')

      if (nsv /= numgrps) stop 'numtracers attribute does not equal the number of groups in sources_path'

      allocate(tracers%data(nsv))
      write(*,*) "Tracers data shape: ",shape(tracers%data(nsv))
      
      !   write(*,*) numgrps, grpids
      ! Allocate arrays for tracers and source variables
      allocate(cm(2-ih:i1+ih,2-jh:j1+jh,k1,nsv))
      allocate(c0(2-ih:i1+ih,2-jh:j1+jh,k1,nsv))
      allocate(cp(2-ih:i1+ih,2-jh:j1+jh,k1,nsv))

      c0 = 0.
      cp = 0.
      cm = 0.

      allocate(source(totnumsources,time_size))
      allocate(source_meta(totnumsources,6)) ! this contains i,j,k,
      !   allocate(source(time_size))

      ! Read global attribute `nsv`

      ! Loop over each tracer group based on `nsv`
      !   write(*,*) shape(c0)
      src_counter = 0

      do itsv = 1, nsv
         !  ! Generate group name based on a known pattern (e.g., "tracer_1", "tracer_2", etc.)
         !  write(group_name, '(A,I0)') 'tracer_', itsv
         !  write(*,*) group_name
         !  ! Access the group by name
         !  retval = nf90_inq_ncid(ncid, trim(group_name), groupid)
         !  call nchandle_error(retval, 'Error: Could not access group '//trim(group_name))
         groupid = grpids(itsv)
         print *, 'Reading group:', groupid
         
         ! Retrieve `numsources` attribute within each group
         retval = nf90_get_att(groupid, nf90_global, 'numsources', numsources)
         call nchandle_error(retval, 'Error: Could not retrieve "numsources" attribute in group ')

         ! Read the 'init' variable in the current group
         retval = nf90_inq_varid(groupid, 'init', varid)
         call nchandle_error(retval, 'Error: Could not read "init" variable in group ')

         call read_tracer_attributes(groupid,tracers%data(itsv))
         tracers%data(itsv)%trac_idx = itsv

         call tracers%data(itsv)%print_properties()

         !!!!!!!!!!!!!!!! DEBUG
         !  ! Inquire the dimensions of `init` in the NetCDF file
         !  retval = nf90_inquire_variable(groupid, varid, dimids=dimids)
         !  call nchandle_error(retval, 'Error: Could not inquire dimensions of "init" variable in group ')

         !  ! Print dimension sizes for each dimension of init in the file

         !  do i = 1, 3
         !     retval = nf90_inquire_dimension(groupid, dimids(i),dimname, len = dim_size)
         !     call nchandle_error(retval, 'Error: Could not inquire dimension size')
         !     write(*, *) 'File init dimension ', i, dimname,':', dim_size
         !  end do

         !  write(*, *) 'Expected init dimensions (X, Y, Z): ', (2-ih+i1+ih), (2-jh+j1+jh), k1

         !!!!!!!!!!!!!

         retval = nf90_get_var(groupid, varid, c0(2:i1,2:j1,1:k1-1,itsv))
         call nchandle_error(retval, 'Error: Could not read "init" variable in group ')

         ! Loop over each source variable in this group based on `numsources`
         do itsrc = 1, numsources
            src_counter = src_counter + 1
            ! Generate source variable name, e.g., "source_1", "source_2", etc.
            write(varname, '(A,I0)') 'source_', itsrc
            retval = nf90_inq_varid(groupid, trim(varname), src_varid)
            call nchandle_error(retval, 'Error: Could not read variable '//trim(varname)//' in group ')
            ! function nf90_inquire_attribute(ncid, varid, name, xtype, len, attnum)

            retval = nf90_get_att(groupid, src_varid, 'x',src_x)
            call nchandle_error(retval, 'Error: Could not read attribute x')
            retval = nf90_get_att(groupid, src_varid, 'y',src_y)
            call nchandle_error(retval, 'Error: Could not read attribute y')
            retval = nf90_get_att(groupid, src_varid, 'z',src_z)
            call nchandle_error(retval, 'Error: Could not read attribute z')
            retval = nf90_get_var(groupid, src_varid, source(src_counter,:))
            call nchandle_error(retval, 'Error: Could not read variable '//trim(varname)//' in group ')
            source_meta(src_counter,:) = (/src_x,src_y,src_z,itsv,itsrc,src_counter/)

            print *, 'Read source variable:', trim(varname), 'in group'
         end do
      end do
      ! Close the NetCDF file
      retval = nf90_close(ncid)
      call nchandle_error(retval, 'Error: Could not close NetCDF file')




      ! The initial conditions should apply also to the "previous iteration", so no time gradients
      ! Also the ghost cells are filled with the values of the initial conditions at the edges
      do i=1,ih
         c0(2-i,:,:,:) = c0(2-ih+1,:,:,:) ! West side
         c0(i1+i,:,:,:) = c0(i1,:,:,:)    ! east boundary
         c0(:,2-i,:,:) = c0(:,2-jh+2,:,:) ! south boundary
         c0(:,j1+i,:,:) = c0(:,j1,:,:)    ! north boundary
      end do

      cm = c0

   end subroutine load_tracer_init_and_sources

   subroutine read_tracer_attributes(groupid, tracer)
      use modtracer_type , only: T_tracer
      use netcdf_utils, only: nchandle_error
      use netcdf
      implicit none
      integer, intent(in) :: groupid
      type(T_tracer), intent(inout) :: tracer
      integer :: retval
      character(len=64) :: tracname
      character(len=64) :: group_name
      character(len=64) :: traclong
      character(len=16) :: unit
      real(real32) :: molar_mass
      integer :: lemis_int

      retval = nf90_inq_grpname(groupid, group_name)
      print *, 'Group name: ', group_name
      call nchandle_error(retval, 'Error: Could not inquire group name')
      retval = nf90_get_att(groupid, nf90_global, 'tracname', tracname)

      if (trim(group_name) /= trim(tracname)) then
         write(*,*) "ERROR: Group name: ", trim(group_name), "Must Equal  Tracname attribute: ", trim(tracname)
         stop "Mismatch between group name and tracname attribute."
      end if
      retval = nf90_get_att(groupid, nf90_global, 'traclong', traclong)
      retval = nf90_get_att(groupid, nf90_global, 'unit', unit)
      retval = nf90_get_att(groupid, nf90_global, 'molar_mass', molar_mass)
      retval = nf90_get_att(groupid, nf90_global, 'lemis', lemis_int)

      tracer%tracname = trim(tracname)
      tracer%traclong = trim(traclong)
      tracer%unit = trim(unit)
      tracer%molar_mass = molar_mass
      tracer%lemis = lemis_int /= 0
      
   end subroutine read_tracer_attributes


   subroutine init_concentration_output_nc()
      use modglobal, only: xt,yt,zt,nsv
      use config, only: outputfile_path
      use netcdf_utils, only: create_concentration_file_nc
      use modmpi,only: my_id
      character(len=256) :: outputfile
      write(outputfile, '(A,".",I3.3)') trim(outputfile_path), my_id+1
      !   print *, outputfile
      call create_concentration_file_nc(output_ncid,outputfile,xt,yt,zt,nsv)
   end subroutine init_concentration_output_nc

   subroutine write_concentration_timeloop()
      use modglobal, only: rsts,next_save,i1,j1,kmax,dt
      use config, only: output_save_interval
      use netcdf_utils, only: write_concentration_nc
      use modmpi,only: is_root


      if (rsts>=next_save) then
         call write_concentration_nc(output_ncid, c0(2:i1,2:j1,1:kmax,:), rsts)
         next_save = next_save + output_save_interval

         if (is_root) then
            write (*,*) "Saved concentration Field at: ", rsts
            write(*,*)  "Max value:", maxval(abs(c0))
            write(*,*)  "Min value:", minval(c0)
            write(*,*)  "dt (ms):", dt
         endif
         if (maxval(abs(c0)) > 1e3) then
            call close_concentration_nc
            stop "CONCENTRTION VALUE EXCEEDED 1e3"
         endif
      endif

   end subroutine write_concentration_timeloop

   subroutine close_concentration_nc()
      use netcdf
      use netcdf_utils, only: nchandle_error
      integer :: retval
      ! Close the NetCDF file
      retval = nf90_close(output_ncid)
      call nchandle_error(retval, 'Error: Could not close the NetCDF file')

   end subroutine close_concentration_nc
end module modtracer
