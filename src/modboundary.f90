module modboundary
   implicit none

contains

   subroutine apply_bc()
      use time_integrate, only: rkstep
      use modtracer, only: c0,cm
      use config, only: rkmethod
      use modglobal, only: nsv
      integer :: i
      do i = 1, nsv

         ! if (rkstep == 1) then
         call x_bc(c0(:,:,:,i))
         call y_bc(c0(:,:,:,i))
         call z_bc(c0(:,:,:,i))
      end do

      ! if (rkmethod>1) then
      !    call x_bc(cm)
      !    call y_bc(cm)
      !    call z_bc(cm)
      ! endif
      ! endif

   end subroutine apply_bc

   subroutine z_bc(field)
      ! Neumann BC for top and bottom, generalized to accept an array
      use modibm, only: groundlibm
      use modglobal, only: k1, kh, kmax, i1,ih,j1,jh,k1
      real, intent(inout) :: field(2-ih:i1+ih,2-jh:j1+jh,k1)
      integer :: i, j

      field(:,:,k1) = field(:,:,k1-1)  ! Top level, mirror

      ! Only apply boundary conditions where there are no buildings above
      do i = 1, i1
         do j = 1, j1
            if (.not.groundlibm(i,j)) then
               field(i,j,1) = field(i,j,2)  ! Ground level, mirror
            endif
         end do
      end do
   end subroutine z_bc

   subroutine x_bc(field)
      ! Neumann BC for west and east boundaries, generalized to accept an array
      use modglobal, only: i1,ih,j1,jh,k1
      real, intent(inout) :: field(2-ih:i1+ih,2-jh:j1+jh,k1)
      integer :: i

      do i = 1, ih
         ! WEST
         field(ih-i,:,:) = field(ih-i+1,:,:)
         ! EAST
         field(i1+i,:,:) = field(i1+i-1,:,:)
      end do
   end subroutine x_bc

   subroutine y_bc(field)
      ! Neumann BC for north and south boundaries, generalized to accept an array
      use modglobal, only: i1,ih,j1,jh,k1
      real, intent(inout) :: field(2-ih:i1+ih,2-jh:j1+jh,k1)
      integer :: j

      do j = 1, jh
         ! NORTH
         field(:,jh-j,:) = field(:,jh-j+1,:)
         ! SOUTH
         field(:,j1+j,:) = field(:,j1+j-1,:)
      end do
   end subroutine y_bc

end module modboundary
