module modmpi
   use mpi_f08
   implicit none
   integer :: ierr, my_id, mpi_size
   logical :: is_root
contains
   subroutine init_mpi()
      call MPI_Init(ierr)
      call MPI_Comm_rank(MPI_COMM_WORLD, my_id, ierr)
      call MPI_Comm_size(MPI_COMM_WORLD, mpi_size, ierr)

      is_root = my_id == 0
      print *, 'Hello from rank ', my_id, ' of ', mpi_size

   end subroutine init_mpi

   subroutine exit_mpi()
      call MPI_Finalize(ierr)

   end subroutine exit_mpi
end module modmpi
