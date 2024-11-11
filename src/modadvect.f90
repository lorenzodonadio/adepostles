!> \file modadvect.f90
!!  Does advection with a 1st order upwind scheme.
!! \par Revision list
!! \par Authors
!  This file is part of DALES.
!
! DALES is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 3 of the License, or
! (at your option) any later version.
!
! DALES is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.
!
!  Copyright 1993-2009 Delft University of Technology, Wageningen University, Utrecht University, KNMI
!

module modadvect
   use iso_fortran_env, only:real32
contains
   !> Advection at cell center
   subroutine apply_advection()
      use modtracer, only: c0,cp
      use modglobal, only: nsv
      integer :: i
      do i = 1, nsv
         call advecc_upw(c0(:,:,:,i),cp(:,:,:,i))
      end do

   end subroutine apply_advection

   subroutine advecc_upw(a_in,a_out)

      use modglobal, only : i1,ih,j1,jh,k1,kmax,dxi,dyi,dzi
      use modfields, only : u0, v0, w0, rhobf
      implicit none

      real(real32), dimension(2-ih:i1+ih,2-jh:j1+jh,k1), intent(in)  :: a_in !< Input: the cell centered field
      real(real32), dimension(2-ih:i1+ih,2-jh:j1+jh,k1), intent(inout) :: a_out !< Output: the tendency

      real(real32), dimension(2-ih:i1+ih,2-jh:j1+jh,k1) :: put
      real(real32), dimension(2-ih:i1+ih,2-jh:j1+jh,k1) :: rho_a_in
      integer :: i,j,k


      do k=1,k1
         do j=2-jh,j1+jh
            do i=2-ih,i1+ih
               rho_a_in(i,j,k)=rhobf(k)*a_in(i,j,k)
            end do
         end do
      end do

      do k=1,k1
         do j=2,j1
            do i=2,i1+1
               if( u0(i,j,k) > 0 ) then
                  put(i,j,k) = rho_a_in(i-1,j,k)
               else
                  put(i,j,k) = rho_a_in(i,j,k)
               endif
            enddo
         enddo
      enddo

      do k=1,k1
         do j=2,j1
            do i=2,i1
               a_out(i,j,k) = a_out(i,j,k) - &
                  (1./rhobf(k))*(u0(i+1,j,k)*put(i+1,j,k)-u0(i,j,k)*put(i,j,k))*dxi
            enddo
         enddo
      enddo

      do k=1,k1
         do j=2,j1+1
            do i=2,i1
               if( v0(i,j,k) > 0 ) then
                  put(i,j,k) = rho_a_in(i,j-1,k)
               else
                  put(i,j,k) = rho_a_in(i,j,k)
               endif
            enddo
         enddo
      enddo
      do k=1,k1
         do j=2,j1
            do i=2,i1
               a_out(i,j,k) = a_out(i,j,k) - &
                  (1./rhobf(k))*(v0(i,j+1,k)*put(i,j+1,k)-v0(i,j,k)*put(i,j,k))*dyi
            enddo
         enddo
      enddo

      put(2:i1,2:j1, 1) = 0
      put(2:i1,2:j1,k1) = 0
      do k=2,kmax
         do j=2,j1
            do i=2,i1
               if( w0(i,j,k) > 0 ) then
                  put(i,j,k) = rho_a_in(i,j,k-1)
               else
                  put(i,j,k) = rho_a_in(i,j,k)
               endif
            enddo
         enddo
      enddo
      do k=1,kmax
         do j=2,j1
            do i=2,i1
               a_out(i,j,k) = a_out(i,j,k) - &
                  (1./rhobf(k))*(w0(i,j,k+1)*put(i,j,k+1)-w0(i,j,k)*put(i,j,k))*dzi
            enddo
         enddo
      enddo

   end subroutine advecc_upw

   subroutine advecc_kappa(a_in,a_out)
      use modglobal, only : i1,i2,ih,j1,j2,jh,k1,kmax,dxi,dyi,dzf
      use modfields, only : u0, v0, w0, rhobf
      implicit none
      real(real32), dimension(2-ih:i1+ih,2-jh:j1+jh,k1), intent(in) :: a_in
      real(real32), dimension(2-ih:i1+ih,2-jh:j1+jh,k1), intent(inout) :: a_out

      real      d1,d2,cf
      real :: d1m, d2m, d1p, cfm, cfp, work
      integer   i,j,k,k_low,k_high

      ! find the lowest and highest k level with a non-zero value in a_in
      k_low = -1
      do k=1,k1
         if (any(a_in(:,:,k) /= 0.)) then
            k_low = k
            exit
         endif
      enddo
      if (k_low == -1) then
         ! a_in == zero
         return
      endif

      do k=k1,1,-1
         if (any(a_in(:,:,k) /= 0.)) then
            k_high = k
            exit
         endif
      enddo

      if (k_low <= 3) then
         ! vertical advection from layer 1 to 2, special case. k=2
         do j=2,j1
            do i=2,i1 ! YES
               d1m = 0
               d2m = rhobf(2) * a_in(i,j,2) - rhobf(1) * a_in(i,j,1)
               cfm = rhobf(1) * a_in(i,j,1)
               d1p = rhobf(2) * a_in(i,j,2) - rhobf(3) * a_in(i,j,3)
               cfp = rhobf(2) * a_in(i,j,2)

               if (w0(i,j,2) > 0) then
                  d1 = d1m
                  d2 = d2m
                  cf = cfm
               else
                  d1 = d1p
                  d2 = -d2m
                  cf = cfp
               end if

               work = cf + &
                  min(abs(d1), abs(d2), abs((d1/6.0) + (d2/3.0))) * &
                  (sign(0.5, d1) + sign(0.5, d2))

               work = work * w0(i,j,2)
               a_out(i,j,1) = a_out(i,j,1) - (1./(rhobf(1)*dzf(1)))*work
               a_out(i,j,2) = a_out(i,j,2) + (1./(rhobf(2)*dzf(2)))*work
            end do
         end do
      end if

      !do k=1,kmax
      do k= max(k_low-1,1), min(k_high+2, kmax) ! loop accesses k-2, k-1, k, k+1
         do j=2,j1
            do i=2,i2 ! YES
               d2m = a_in(i  ,j,k)-a_in(i-1,j,k)
               d1m = a_in(i-1,j,k)-a_in(i-2,j,k)
               d1p = a_in(i  ,j,k)-a_in(i+1,j,k)
               cfm = a_in(i-1,j,k)
               cfp = a_in(i  ,j,k)

               if (u0(i,j,k) > 0) then
                  d1 = d1m
                  d2 = d2m
                  cf = cfm
               else
                  d1 = d1p
                  d2 = -d2m
                  cf = cfp
               end if

               work = cf + &
                  min(abs(d1), abs(d2), abs((d1/6.0) + (d2/3.0))) * &
                  (sign(0.5, d1) + sign(0.5, d2))

               work = work * u0(i,j,k) * dxi
               a_out(i-1,j,k) = a_out(i-1,j,k) - work
               a_out(i,j,k)   = a_out(i,j,k)   + work
            end do
         end do
         !  end do

         !  do k=1,kmax
         do j=2,j2
            do i=2,i1 ! YES
               d1m = a_in(i,j-1,k)-a_in(i,j-2,k)
               d1p = a_in(i,j  ,k)-a_in(i,j+1,k)
               d2m = a_in(i,j  ,k)-a_in(i,j-1,k)
               cfm = a_in(i,j-1,k)
               cfp = a_in(i,j  ,k)

               if (v0(i,j,k) > 0) then
                  d1 = d1m
                  d2 = d2m
                  cf = cfm
               else
                  d1 = d1p
                  d2 = -d2m
                  cf = cfp
               end if

               work = cf + &
                  min(abs(d1), abs(d2), abs((d1/6.0) + (d2/3.0))) * &
                  (sign(0.5, d1) + sign(0.5, d2))

               work = work * v0(i,j,k) * dyi
               a_out(i,j-1,k) = a_out(i,j-1,k) - work
               a_out(i,j,k)   = a_out(i,j,k)   + work
            end do
         end do
         !  end do

         !  do k=3,kmax
         if (k >= 3) then
            do j=2,j1
               do i=2,i1 ! YES
                  d1m = rhobf(k-1) * a_in(i,j,k-1) - rhobf(k-2) * a_in(i,j,k-2)
                  d2m = rhobf(k)   * a_in(i,j,k  ) - rhobf(k-1) * a_in(i,j,k-1)
                  d1p = rhobf(k)   * a_in(i,j,k  ) - rhobf(k+1) * a_in(i,j,k+1)
                  cfm = rhobf(k-1) * a_in(i,j,k-1)
                  cfp = rhobf(k)   * a_in(i,j,k  )

                  if (w0(i,j,k) > 0) then
                     d1 = d1m
                     d2 = d2m
                     cf = cfm
                  else
                     d1 = d1p
                     d2 = -d2m
                     cf = cfp
                  end if

                  work = cf + &
                     min(abs(d1), abs(d2), abs((d1/6.0) + (d2/3.0))) * &
                     (sign(0.5, d1) + sign(0.5, d2))

                  work = work * w0(i,j,k)
                  a_out(i,j,k-1) = a_out(i,j,k-1) - (1./(rhobf(k-1)*dzf(k-1)))*work
                  a_out(i,j,k)   = a_out(i,j,k)   + (1./(rhobf(k)  *dzf(k)  ))*work
               end do
            end do
         end if
      end do

   end subroutine advecc_kappa
end module modadvect
