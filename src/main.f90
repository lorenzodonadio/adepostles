program adepostles
   use config, only: read_config_file,runtime,field_load_chunk_size
   use modglobal, only: load_dimensions,validate_time_size,init_global,maxtime,simtime
   use modfields,only:allocate_fields, load_fields_intimeloop,init_interp_fields,interpolate_fields_to_simtime
   use time_integrate, only: time_step,allocate_k
   use modadvect, only: apply_advection
   use moddiff,only: init_diff,apply_diff
   use modboundary, only: apply_bc
   use modibm, only : init_ibm,apply_ibm
   use modtracer, only: load_tracer_init_and_sources,init_concentration_output_nc,&
      write_concentration_timeloop,close_concentration_nc,apply_source
   use modmpi, only: init_mpi,exit_mpi
   implicit none

   ! Read the configuration from the namelist file
   call init_mpi
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
   call load_tracer_init_and_sources
   call init_concentration_output_nc
!    stop 'Init completed succesfully'

   do while (simtime < maxtime)

      call write_concentration_timeloop
      call load_fields_intimeloop
      call interpolate_fields_to_simtime

      call apply_source

      call apply_diff
      ! call advecc_kappa(c0,cp)
      call apply_advection

      call time_step
      call apply_ibm
      call apply_bc
      !   stop 'ACTUALLY ALL GOOD :)'

   end do


   call close_concentration_nc
   call exit_mpi

   write (*,*) "Simulation Finished Succesfully at: ", simtime, "ms"

end program adepostles
