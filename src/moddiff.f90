module moddiff
   implicit none

contains

   subroutine init_diff()
      if(lanisotrop) then
         ! Anisotropic diffusion scheme  https://doi.org/10.1029/2022MS003095
         ! length scale in TKE equation is delta z (private communication with Marat)
         if ((dx.ne.dy) .and. myid == 0) stop "The anisotropic diffusion assumes dx=dy."
         do k = 1,k1
            deltai   (k) = 1./dzf(k)       !overrules deltai (k) = 1/delta(k) as defined in initglobal
            anis_fac (k) = (dx/dzf(k))**2  !assumes dx=dy. is used to enhance horizontal diffusion
         end do
      else
         do k = 1,k1
            anis_fac (k) = 1.   !horizontal = vertical diffusion
         end do

      endif
   end subroutine init_diff

   subroutine diffc (a_in,a_out,flux)

      use modglobal, only : i1,ih,i2,j1,jh,j2,k1,kmax,dx2i,dzf,dy2i,dzh
      use modfields, only : rhobf,rhobh,ksfc
      implicit none

      real(field_r), intent(in)    :: a_in(2-ih:i1+ih,2-jh:j1+jh,k1)
      real(field_r), intent(inout) :: a_out(2-ih:i1+ih,2-jh:j1+jh,k1)
      real, intent(in)    :: flux (i2,j2)

      integer i,j,k,jm,jp,km,kp,kmin


      do j=2,j1
         jp=j+1
         jm=j-1

         do i=2,i1
            kmin = ksfc(i,j)
!cibm       do k=2,kmax
            do k=kmin+1,kmax
               kp=k+1
               km=k-1

               a_out(i,j,k) = a_out(i,j,k) &
                  +  0.5 * ( &
                  ( (ekh(i+1,j,k)+ekh(i,j,k))*(a_in(i+1,j,k)-a_in(i,j,k)) &
                  -(ekh(i,j,k)+ekh(i-1,j,k))*(a_in(i,j,k)-a_in(i-1,j,k)))*dx2i * anis_fac(k) &
                  + &
                  ( (ekh(i,jp,k)+ekh(i,j,k)) *(a_in(i,jp,k)-a_in(i,j,k)) &
                  -(ekh(i,j,k)+ekh(i,jm,k)) *(a_in(i,j,k)-a_in(i,jm,k)) )*dy2i * anis_fac(k) &
                  + &
                  ( rhobh(kp)/rhobf(k) * (dzf(kp)*ekh(i,j,k) + dzf(k)*ekh(i,j,kp)) &
                  *  (a_in(i,j,kp)-a_in(i,j,k)) / dzh(kp)**2 &
                  - &
                  rhobh(k)/rhobf(k) * (dzf(km)*ekh(i,j,k) + dzf(k)*ekh(i,j,km)) &
                  *  (a_in(i,j,k)-a_in(i,j,km)) / dzh(k)**2           )/dzf(k) &
                  )

            end do
         end do
      end do

      do j=2,j1
         do i=2,i1
            kmin = ksfc(i,j)
            a_out(i,j,kmin) = a_out(i,j,kmin) &
               + 0.5 * ( &
               ( (ekh(i+1,j,kmin)+ekh(i,j,kmin))*(a_in(i+1,j,kmin)-a_in(i,j,kmin)) &
               -(ekh(i,j,kmin)+ekh(i-1,j,kmin))*(a_in(i,j,kmin)-a_in(i-1,j,kmin)) )*dx2i * anis_fac(kmin) &
               + &
               ( (ekh(i,j+1,kmin)+ekh(i,j,kmin))*(a_in(i,j+1,kmin)-a_in(i,j,kmin)) &
               -(ekh(i,j,kmin)+ekh(i,j-1,kmin))*(a_in(i,j,kmin)-a_in(i,j-1,kmin)) )*dy2i * anis_fac(kmin) &
               + &
               ( rhobh(kmin+1)/rhobf(kmin) * (dzf(kmin+1)*ekh(i,j,kmin) + dzf(kmin)*ekh(i,j,kmin+1)) &
               *  (a_in(i,j,kmin+1)-a_in(i,j,kmin)) / dzh(kmin+1)**2 &
               + rhobh(kmin)/rhobf(kmin)*flux(i,j) *2.                        )/dzf(kmin) &
               )

         end do
      end do

   end subroutine diffc

end module moddiff
