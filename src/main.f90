program adepostles
   use config, only: field_dump_path, read_config_file
   use test_nc, only: do_test_nc
   use modglobal, only: total_chunks,load_dimensions,validate_time_size
   use modfields,only:allocate_fields,load_fields_chunk
   implicit none
   integer :: chunk
   ! Read the configuration from the namelist file
   call read_config_file

   print *, field_dump_path
   ! Proceed with using the configuration variables as needed
   call do_test_nc(field_dump_path)
   call load_dimensions(field_dump_path)
   call validate_time_size
   call allocate_fields
   call load_fields_chunk(field_dump_path,1)
end program adepostles
