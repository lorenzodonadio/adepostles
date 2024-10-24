module moddiff
   use iso_fortran_env
   implicit none
   real,allocatable :: anis_fac(:)
contains

   subroutine init_diff()
      use modglobal, only : ih,i1,jh,j1,k1,dx,dy,dzf
      use config, only: lanisotrop

      implicit none

      integer   :: k
      allocate(anis_fac(k1))
      if(lanisotrop) then
         ! Anisotropic diffusion scheme  https://doi.org/10.1029/2022MS003095
         ! length scale in TKE equation is delta z (private communication with Marat)
         if ((dx.ne.dy)) stop "The anisotropic diffusion assumes dx=dy." !what if not ? cant we have dxhor = (dx+dy/2)
         do k = 1,k1
            ! deltai   (k) = 1./dzf(k)       !overrules deltai (k) = 1/delta(k) as defined in initglobal
            anis_fac (k) = (dx/dzf(k))**2  !assumes dx=dy. is used to enhance horizontal diffusion
         end do
      else
         do k = 1,k1
            anis_fac (k) = 1.   !horizontal = vertical diffusion
         end do

      endif
   end subroutine init_diff

!    subroutine diffc (a_in,a_out,flux)
   subroutine diffc (a_in,a_out)

      use modglobal, only : i1,ih,i2,j1,jh,j2,k1,kmax,dx2i,dzf,dy2i,dzh
      use modfields, only : rhobf,rhobh,ekh0 !ksfc
      implicit none

      real(real32), intent(in)    :: a_in(2-ih:i1+ih,2-jh:j1+jh,k1)
      real(real32), intent(inout) :: a_out(2-ih:i1+ih,2-jh:j1+jh,k1)
      !   real, intent(in)    :: flux (i2,j2)

      integer i,j,k,jm,jp,km,kp,kmin


      do j=2,j1
         jp=j+1
         jm=j-1

         do i=2,i1
            ! kmin = ksfc(i,j)
!cibm       do k=2,kmax
            ! do k=kmin+1,kmax
            do k=2,kmax
               kp=k+1
               km=k-1

               a_out(i,j,k) = a_out(i,j,k) &
                  +  0.5 * ( &
                  ( (ekh0(i+1,j,k)+ekh0(i,j,k))*(a_in(i+1,j,k)-a_in(i,j,k)) &
                  -(ekh0(i,j,k)+ekh0(i-1,j,k))*(a_in(i,j,k)-a_in(i-1,j,k)))*dx2i * anis_fac(k) &
                  + &
                  ( (ekh0(i,jp,k)+ekh0(i,j,k)) *(a_in(i,jp,k)-a_in(i,j,k)) &
                  -(ekh0(i,j,k)+ekh0(i,jm,k)) *(a_in(i,j,k)-a_in(i,jm,k)) )*dy2i * anis_fac(k) &
                  + &
                  ( rhobh(kp)/rhobf(k) * (dzf(kp)*ekh0(i,j,k) + dzf(k)*ekh0(i,j,kp)) &
                  *  (a_in(i,j,kp)-a_in(i,j,k)) / dzh(kp)**2 &
                  - &
                  rhobh(k)/rhobf(k) * (dzf(km)*ekh0(i,j,k) + dzf(k)*ekh0(i,j,km)) &
                  *  (a_in(i,j,k)-a_in(i,j,km)) / dzh(k)**2           )/dzf(k) &
                  )

            end do
         end do
      end do

      do j=2,j1
         do i=2,i1
            ! kmin = ksfc(i,j)
            kmin = 1
            a_out(i,j,kmin) = a_out(i,j,kmin) &
               + 0.5 * ( &
               ( (ekh0(i+1,j,kmin)+ekh0(i,j,kmin))*(a_in(i+1,j,kmin)-a_in(i,j,kmin)) &
               -(ekh0(i,j,kmin)+ekh0(i-1,j,kmin))*(a_in(i,j,kmin)-a_in(i-1,j,kmin)) )*dx2i * anis_fac(kmin) &
               + &
               ( (ekh0(i,j+1,kmin)+ekh0(i,j,kmin))*(a_in(i,j+1,kmin)-a_in(i,j,kmin)) &
               -(ekh0(i,j,kmin)+ekh0(i,j-1,kmin))*(a_in(i,j,kmin)-a_in(i,j-1,kmin)) )*dy2i * anis_fac(kmin) &
               + &
               ( rhobh(kmin+1)/rhobf(kmin) * (dzf(kmin+1)*ekh0(i,j,kmin) + dzf(kmin)*ekh0(i,j,kmin+1)) &
               *  (a_in(i,j,kmin+1)-a_in(i,j,kmin)) / dzh(kmin+1)**2 )/dzf(kmin))
            !    + &
            !    ( rhobh(kmin+1)/rhobf(kmin) * (dzf(kmin+1)*ekh0(i,j,kmin) + dzf(kmin)*ekh0(i,j,kmin+1)) &
            !    *  (a_in(i,j,kmin+1)-a_in(i,j,kmin)) / dzh(kmin+1)**2 &
            !    + rhobh(kmin)/rhobf(kmin)*flux(i,j) *2.                        )/dzf(kmin) &
            !    )

         end do
      end do

   end subroutine diffc

end module moddiff
