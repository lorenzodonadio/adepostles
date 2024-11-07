program adepostles
   use config, only: read_config_file,runtime,field_load_chunk_size
   use modglobal, only: load_dimensions,validate_time_size,init_global,maxtime,simtime
   use modfields,only:allocate_fields, load_fields_intimeloop,init_interp_fields,interpolate_fields_to_simtime,&
      c0,cp,init_concentration_output_nc,write_concentration_timeloop,close_concentration_nc
   use time_integrate, only: time_step,allocate_k
   use advec_kappa, only: advecc_kappa
   use advec_upw, only: advecc_upw
   use moddiff,only: init_diff,diffc
   use modboundary, only: apply_bc,apply_source
   use modibm, only : init_ibm,apply_ibm
   implicit none

   ! Read the configuration from the namelist file
   call read_config_file

   ! Proceed with using the configuration variables as needed
   ! call do_test_nc(field_dump_path)
   call load_dimensions
   call init_global !this needs to run after loading the dimensions, maybe we put it inside the subroutine? but for now here to be explicit
   call validate_time_size
   call allocate_fields
   call allocate_k
   call init_interp_fields
   call init_diff
   call init_ibm
   call init_concentration_output_nc

   do while (simtime < maxtime)

      call write_concentration_timeloop
      call load_fields_intimeloop
      call interpolate_fields_to_simtime

      call diffc(c0,cp)
      ! call advecc_kappa(c0,cp)
      call advecc_upw(c0,cp)

      call time_step
      call apply_source
      call apply_ibm
      call apply_bc

   end do


   call close_concentration_nc

   write (*,*) "Simulation Finished Succesfully at: ", simtime, "ms"

end program adepostles
