module modboundary
   implicit none

   real :: r6_11 = 0.5454545454545454
   real :: r3_11 = 0.2727272727272727
   real :: r24_33 = 0.7272727272727273
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


      ! field(:,:,k1) = field(:,:,k1-1)  ! Top level, mirror
      ! 2nd ORder
      field(:,:,k1) = r6_11*field(:,:,k1-2)+r24_33*field(:,:,k1-2)-r3_11*field(:,:,k1-3)

      ! Only apply boundary conditions where there are no buildings above
      do i = 1, i1
         do j = 1, j1
            if (.not.groundlibm(i,j)) then
               ! field(i,j,1) = field(i,j,2)  ! Ground level, mirror

               field(i,j,1) = r6_11*field(i,j,2)+r24_33*field(i,j,3)-r3_11*field(i,j,4)
            endif
         end do
      end do
   end subroutine z_bc

   subroutine x_bc(field)
      ! BC for west and east boundaries, generalized to accept an array
      use modglobal, only: i1,ih,j1,jh,k1
      use config, only: xboundary
      real, intent(inout) :: field(2-ih:i1+ih,2-jh:j1+jh,k1)
      integer :: i

      select case(xboundary)
       case(11)  ! NEUMANN 1st Order
         do i = 1, ih
            field(ih-i,:,:) = field(ih-i+1,:,:) ! WEST
            field(i1+i,:,:) = field(i1+i-1,:,:) ! EAST
         end do

       case(12)  ! NEUMANN 2nd Order
         do i = 1, ih
            ! WEST
            field(ih-i,:,:) = r6_11*field(ih-i+1,:,:)+r24_33*field(ih-i+2,:,:)-r3_11*field(ih-i+3,:,:)
            ! EAST
            field(i1+i,:,:) = r6_11*field(i1+i-1,:,:)+r24_33*field(i1+i-2,:,:)-r3_11*field(i1+i-3,:,:)
         end do
       case(2)  ! PERIODIC
         field(ih-2:ih,:,:) = field(i1-2:i1,:,:)
         field(i1:i1+2,:,:) = field(ih-2:ih,:,:)
      end select

   end subroutine x_bc

   subroutine y_bc(field)
      ! Neumann BC for north and south boundaries, generalized to accept an array
      use modglobal, only: i1,ih,j1,jh,k1
      use config, only: yboundary

      real, intent(inout) :: field(2-ih:i1+ih,2-jh:j1+jh,k1)
      integer :: j

      select case(yboundary)
       case(11)  ! NEUMANN 1st Order
         do j = 1, jh
            field(:,jh-j,:) = field(:,jh-j+1,:) ! NORTH
            field(:,j1+j,:) = field(:,j1+j-1,:) ! SOUTH
         end do
       case(12)  ! NEUMANN 2nd Order
         do j = 1, jh
            ! NORTH
            field(:,jh-j,:) = r6_11*field(:,jh-j+1,:)+r24_33*field(:,jh-j+2,:)-r3_11*field(:,jh-j+3,:)
            ! SOUTH
            field(:,jh+j,:) = r6_11*field(:,jh+j-1,:)+r24_33*field(:,jh+j-2,:)-r3_11*field(:,jh+j-3,:)
         end do
       case(2)  ! PERIODIC
         field(:,jh-2:jh,:) = field(:,j1-2:j1,:)
         field(:,j1:j1+2,:) = field(:,jh-2:jh,:)
      endselect

   end subroutine y_bc
end module modboundary
