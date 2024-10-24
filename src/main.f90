program adepostles
   use config, only: read_config_file,runtime,field_load_chunk_size
   use modglobal, only: load_dimensions,validate_time_size,init_global,maxtime,simtime , ih,i1,jh,j1,k1,dt,tres
   use modfields,only:allocate_fields, load_fields_intimeloop,init_interp_fields,interpolate_fields_to_simtime,c0,cp,u0,v0
   use time_integrate, only: time_step
   use advec_kappa, only: advecc_kappa
   use moddiff,only: init_diff,diffc
   implicit none
   ! Read the configuration from the namelist file
   call read_config_file

   ! Proceed with using the configuration variables as needed
   ! call do_test_nc(field_dump_path)
   call load_dimensions()
   call init_global !this needs to run after loading the dimensions, maybe we put it inside the subroutine? but for now here to be explicit
   call validate_time_size
   call allocate_fields
   call init_interp_fields
   call init_diff
   do while (simtime < maxtime)
      call load_fields_intimeloop
      call interpolate_fields_to_simtime
      write(*,*) 'u0: '
      write(*,*) u0(9:13,10,3)
      write(*,*) 'v0: '
      write(*,*) v0(9:13,10,3)
      write(*,*) 'c0: '
      write(*,*) c0(9:13,10,3)
      call advecc_kappa(c0,cp)
      write(*,*) 'cp adv: '
      write(*,*) cp(9:13,10,3)
      call diffc(c0,cp)
      !Time integration is the last step
      ! c0(9:13,9:13,1:3) = 5.
      write(*,*) 'cp diff: '
      write(*,*) cp(9:13,10,3)
      c0 = c0 + dt*tres*cp
      cp = 0.
      call time_step
   end do
end program adepostles
