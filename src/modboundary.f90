module modboundary
   implicit none

contains
   subroutine apply_source()
      use time_integrate, only: rkstep
      use modfields, only: c0

      if (rkstep == 1) then
         c0(48:50,80:82,3) = 5
      endif
      
   end subroutine apply_source

   subroutine apply_bc()
      call x_bc
      call y_bc
      call z_bc
   end subroutine apply_bc

   subroutine z_bc()
      ! Neumann BC for top and bottom
      ! TODO maybe make this general and accept an array as input?
      ! (2-ih:i1+ih,2-jh:j1+jh,k1))
      use modglobal,only: k1,kh,kmax,i1,j1
      use modfields, only: c0
      use modibm, only: libm
      integer :: i,j
      c0(:,:,k1) = c0(:,:,k1-1) !top level just mirror

      ! only apply boundary conditions where there are no building above
      do i = 1, i1
         do j = 1, j1
            if (.not.libm(i,j,1)) then
               c0(i,j,1) = c0(i,j,2) !ground level just mirror
            endif
         end do
      end do

   end subroutine z_bc

   subroutine x_bc()
      ! Neumann BC for top and bottom
      ! TODO maybe make this general and accept an array as input?
      ! (2-ih:i1+ih,2-jh:j1+jh,k1))
      use modglobal,only: ih,i1
      use modfields, only: c0

      integer :: i
      ! we have 2 ghost cells actually
      do i = 1, ih
         !  WEST
         c0(ih-i,:,:) = c0(ih-i+1,:,:)
         !  EAST
         c0(i1+i,:,:) = c0(ih+i-1,:,:)
      end do

   end subroutine x_bc

   subroutine y_bc()
      ! Neumann BC for top and bottom
      ! TODO maybe make this general and accept an array as input?
      ! (2-ih:i1+ih,2-jh:j1+jh,k1))
      use modglobal,only: jh,j1
      use modfields, only: c0

      integer :: j
      ! we have 2 ghost cells actually
      do j = 1, jh
         !  WEST
         c0(:,jh-j,:) = c0(:,jh-j+1,:)
         !  EAST
         c0(:,j1+j,:) = c0(:,jh+j-1,:)
      end do
   end subroutine y_bc


end module modboundary
