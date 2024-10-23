program adepostles
   use config, only: read_config_file,runtime,field_load_chunk_size
   use modglobal, only: load_dimensions,validate_time_size,init_global,maxtime,simtime
   use modfields,only:allocate_fields, load_fields_intimeloop,init_interp_fields,interpolate_fields_to_simtime , u0
   use time_integrate, only: time_step
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
   do while (simtime < maxtime)
      call load_fields_intimeloop
      call interpolate_fields_to_simtime
      !Time integration is the last step
      call time_step
   end do
end program adepostles
