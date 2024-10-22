program adepostles
   use config, only: field_dump_path, read_config_file
   use modglobal, only:current_chunk, load_dimensions,validate_time_size,init_global
   use modfields,only:allocate_fields,load_fields_chunk
   implicit none
   ! Read the configuration from the namelist file
   call read_config_file

   print *, field_dump_path
   ! Proceed with using the configuration variables as needed
   ! call do_test_nc(field_dump_path)
   call load_dimensions(field_dump_path)
   call init_global() !this needs to run after loading the dimensions, maybe we put it inside the subroutine? but for now here to be explicit
   call validate_time_size
   call allocate_fields
   call load_fields_chunk(field_dump_path,current_chunk)
end program adepostles
