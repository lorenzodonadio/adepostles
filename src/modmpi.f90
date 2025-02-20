module modmpi
   use mpi_f08
   implicit none
   integer :: ierr, my_id, mpi_size
   logical :: is_root

   real     :: CPU_program    !end time
   real     :: CPU_program0   !start time
contains
   subroutine init_mpi()
      call MPI_Init(ierr)
      call MPI_Comm_rank(MPI_COMM_WORLD, my_id, ierr)
      call MPI_Comm_size(MPI_COMM_WORLD, mpi_size, ierr)

      is_root = my_id == 0

      if(myid==0)then
         CPU_program0 = MPI_Wtime()
         write(*,*) 'MPI mpi_size: ', mpi_size
       end if
   end subroutine init_mpi

   subroutine exit_mpi()
      
      if(myid==0)then
         CPU_program = MPI_Wtime() - CPU_program0
         write(6,*)'TOTAL wall time = ', CPU_program
      end if

      call MPI_Finalize(ierr)
   end subroutine exit_mpi
end module modmpi
