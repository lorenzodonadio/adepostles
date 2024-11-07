module modboundary
   implicit none

contains
   subroutine apply_source()
      use time_integrate, only: rkstep
      use modfields, only: c0

      ! if (rkstep == 1) then
      ! c0(48:50,80:82,3) = 1
      c0(51,100:104,3) = 1
      c0(51,62:66,3) = 1
      c0(51,4:8,3) = 1
      ! endif

   end subroutine apply_source

   subroutine apply_bc()
      use time_integrate, only: rkstep
      use modfields, only: c0,cm
      use config, only: rkmethod
      ! if (rkstep == 1) then
      call x_bc(c0)
      call y_bc(c0)
      call z_bc(c0)

      ! if (rkmethod>1) then
      !    call x_bc(cm)
      !    call y_bc(cm)
      !    call z_bc(cm)
      ! endif
      ! endif

   end subroutine apply_bc

   subroutine z_bc(field)
      ! Neumann BC for top and bottom, generalized to accept an array
      use modibm, only: libm
      use modglobal, only: k1, kh, kmax, i1,ih,j1,jh,k1
      real, intent(inout) :: field(2-ih:i1+ih,2-jh:j1+jh,k1)
      integer :: i, j

      field(:,:,k1) = field(:,:,k1-1)  ! Top level, mirror

      ! Only apply boundary conditions where there are no buildings above
      do i = 1, i1
         do j = 1, j1
            if (.not.libm(i,j,1)) then
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
