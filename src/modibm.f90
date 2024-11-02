module modibm
   use iso_fortran_env
   implicit none

   !< Normal immersed boundary layers for incoming x,y,z-velocities
   logical, allocatable :: libm(:,:,:)!, lnorm_x(:,:,:), lnorm_y(:,:,:), lnorm_z(:,:,:)
!    integer,allocatable :: inorm_west(:,:),inorm_east(:,:),inorm_south(:,:),inorm_north(:,:),inorm_top(:,:)
   integer,allocatable :: inorm_ibm(:,:)
contains

   subroutine init_ibm()
      use config, only: lapplyibm,ibm_input_file,ifinputibm
      use modglobal, only: i1,ih,j1,jh,k1,kh,kmax,zm,zt !zm = zh omgg
      real(real32), allocatable :: bc_height(:,:)     !< Height of immersed boundary at grid pos x,y
      integer       :: i, j, k, nsize
      character(128) :: readstring

      if (.not. lapplyibm) return


      allocate(bc_height (i1,j1))
      allocate(libm(2-ih:i1+ih,2-jh:j1+jh,k1))

      nsize = 0 !allboundaries
      libm = .false.
      !! READ FILE
      write(6,*) 'Reading inputfile in modibm'
      open (ifinputibm,file=ibm_input_file)
      do k=1,7
         read (ifinputibm,'(a100)') readstring
         !  read (ifinputibm,*) readstring
         write (6,*) readstring
      enddo

      do j=j1,2,-1  !cstep
         do i=2,i1
            read(ifinputibm,'(F6.1)') bc_height(i,j)
         enddo
      end do
      close(ifinputibm)

      bc_height(1,:)=bc_height(i1,:)
      bc_height(:,1)=bc_height(:,j1)
      write(6,*) 'Succesfully read inputfile in modibm'
      !   write(*,*) bc_height

      do i=2,i1
         do j=2,j1
            do k=1,kmax
               if (zt(k) <= bc_height(i,j)) then  !obstacle height is above mid point of vertical grid
                  libm (i,j,k) = .true.
                  !   TODO make use of ksfc, mainly to improve diffusion loop
                  !   ksfc (i,j)   = k + 1   !half (flux) level
                  !   if (ksfc(i,j).gt.kibm_maxl) then
                  !      kibm_maxl = k
                  !   endif
                  ! write (6,*) 'libm',i+myidx*imax,j+myidy*jmax,i,j,k,libm(i,j,k),bc_height(i+myidx*imax,j+myidy*jmax),zh(ksfc(i,j))
                  !   write (6,*) i,j,k
               endif
            end do
         end do
      end do

      deallocate (bc_height)


      !cstep            0     X      X     X      0       ,building position
      !cstep           i-2   i-1     i    i+1    i+2
      !cstep  libm      F     T      T     T      F
      !cstep  lnorm_x   F     T      F     F      T  , the true points refer to u-positions on the
      !                                                grid box (on its left)
      do k=2,kmax
         do j=2,j1
            do i=2,i1
               !WEST
               if (libm(i,j,k) .and. .not.libm(i-1,j,k)) then
                  nsize = nsize + 1
                  !EAST
               elseif(libm(i,j,k) .and. .not.libm(i+1,j,k)) then
                  nsize = nsize + 1
                  !SOUTH
               elseif (libm(i,j,k) .and. .not.libm(i,j-1,k)) then
                  nsize = nsize + 1
                  !NORTH
               elseif (libm(i,j,k) .and. .not.libm(i,j+1,k)) then
                  nsize = nsize + 1
                  !TOP
               elseif (libm(i,j,k) .and. .not.libm(i,j,k+1)) then
                  nsize = nsize + 1
               endif
            end do
         end do
      end do

      write(*,*) 'IBM DIMENSIONS: ', nsize

      allocate(inorm_ibm(3,nsize))
      !   do k=1,kmax
      !reuse those ints
      nsize = 1
      do k=2,kmax
         do j=2,j1
            do i=2,i1
               !WEST
               if (libm(i,j,k) .and. .not.libm(i-1,j,k)) then
                  inorm_ibm(:,nsize) = (/i,j,k/)
                  nsize = nsize + 1
                  !EAST
               elseif (libm(i,j,k) .and. .not.libm(i+1,j,k)) then
                  inorm_ibm(:,nsize) = (/i,j,k/)
                  nsize = nsize + 1
                  !SOUTH
               elseif (libm(i,j,k) .and. .not.libm(i,j-1,k)) then
                  inorm_ibm(:,nsize) = (/i,j,k/)
                  nsize = nsize + 1
                  !NORTH
               elseif (libm(i,j,k) .and. .not.libm(i,j+1,k)) then
                  inorm_ibm(:,nsize) = (/i,j,k/)
                  nsize = nsize + 1
                  !TOP
               elseif (libm(i,j,k) .and. .not.libm(i,j,k+1)) then
                  inorm_ibm(:,nsize) = (/i,j,k/)
                  nsize = nsize + 1
               endif
            end do
         end do
      end do


      write(*,*) 'Succesfully found normal layers in all directions'
      ! write(*,*) 'ALL BOUNDARIES', inorm_ibm(:,:20)
   end subroutine init_ibm

   subroutine apply_ibm()
      use modglobal, only: i1,j1,kmax
      use modfields, only: c0
      integer       :: i, j, k,nboundary
      real          :: cwest,ceast,csouth,cnorth,ctop,csum
      do nboundary = 1, size(inorm_ibm,2)
         i = inorm_ibm(1,nboundary)
         j = inorm_ibm(2,nboundary)
         k = inorm_ibm(3,nboundary)

         if (c0(i,j,k) > 1e-6) then
            ! write(*,*) i,j,k,c0(i,j,k)
            csum = c0(i+1,j,k)+c0(i-1,j,k)+c0(i,j+1,k)+c0(i,j-1,k)+c0(i,j,k+1)
            c0(i+1,j,k) = c0(i+1,j,k)*c0(i,j,k)/csum
            c0(i-1,j,k) = c0(i-1,j,k)*c0(i,j,k)/csum
            c0(i,j+1,k) = c0(i,j+1,k)*c0(i,j,k)/csum
            c0(i,j-1,k) = c0(i,j-1,k)*c0(i,j,k)/csum
            c0(i,j,k+1) = c0(i,j,k+1)*c0(i,j,k)/csum
         endif
      end do

      do k=2,kmax
         do j=2,j1
            do i=2,i1
               if (libm(i,j,k)) c0(i,j,k) = 0.
            end do
         end do
      end do

   end subroutine apply_ibm
end module modibm

!subroutine with all the boundaries, but i think we dont even need that to be soo complex
!    subroutine init_ibm()
!       use config, only: lapplyibm,ibm_input_file,ifinputibm
!       use modglobal, only: i1,ih,j1,jh,k1,kh,kmax,zm,zt !zm = zh omgg
!       real(real32), allocatable :: bc_height(:,:)     !< Height of immersed boundary at grid pos x,y
!       integer       :: i, j, k, ws, es,ss,ns,ts !, ierr,ii,jj,kk,n  !cstep , kmin
!       character(128) :: readstring

!       if (.not. lapplyibm) return


!       allocate(bc_height (i1,j1))
!       allocate(libm(2-ih:i1+ih,2-jh:j1+jh,k1))

!       ws = 0 !west
!       es = 0 !east
!       ss = 0 !south
!       ns = 0 !north
!       ts = 0 !top
!       libm = .false.
!       !! READ FILE
!       write(6,*) 'Reading inputfile in modibm'
!       open (ifinputibm,file=ibm_input_file)
!       do k=1,7
!          read (ifinputibm,'(a100)') readstring
!          !  read (ifinputibm,*) readstring
!          write (6,*) readstring
!       enddo

!       do j=j1,2,-1  !cstep
!          do i=2,i1
!             read(ifinputibm,'(F6.1)') bc_height(i,j)
!          enddo
!       end do
!       close(ifinputibm)

!       bc_height(1,:)=bc_height(i1,:)
!       bc_height(:,1)=bc_height(:,j1)
!       write(6,*) 'Succesfully read inputfile in modibm'
!       !   write(*,*) bc_height

!       do i=2,i1
!          do j=2,j1
!             do k=1,kmax
!                if (zt(k) <= bc_height(i,j)) then  !obstacle height is above mid point of vertical grid
!                   libm (i,j,k) = .true.
!                   !   TODO make use of ksfc, mainly to improve diffusion loop
!                   !   ksfc (i,j)   = k + 1   !half (flux) level
!                   !   if (ksfc(i,j).gt.kibm_maxl) then
!                   !      kibm_maxl = k
!                   !   endif
!                   ! write (6,*) 'libm',i+myidx*imax,j+myidy*jmax,i,j,k,libm(i,j,k),bc_height(i+myidx*imax,j+myidy*jmax),zh(ksfc(i,j))
!                   !   write (6,*) i,j,k
!                endif
!             end do
!          end do
!       end do

!       deallocate (bc_height)


!       !cstep            0     X      X     X      0       ,building position
!       !cstep           i-2   i-1     i    i+1    i+2
!       !cstep  libm      F     T      T     T      F
!       !cstep  lnorm_x   F     T      F     F      T  , the true points refer to u-positions on the
!       !                                                grid box (on its left)
!       do k=2,kmax
!          do j=2,j1
!             do i=2,i1
!                !WEST
!                if (libm(i,j,k) .and. .not.libm(i-1,j,k)) ws = ws + 1
!                !EAST
!                if (libm(i,j,k) .and. .not.libm(i+1,j,k)) es = es +1
!                !SOUTH
!                if (libm(i,j,k) .and. .not.libm(i,j-1,k)) ss = ss +1
!                !NORTH
!                if (libm(i,j,k) .and. .not.libm(i,j+1,k)) ns = ns +1
!                !TOP
!                if (libm(i,j,k) .and. .not.libm(i,j,k+1)) ts = ts +1
!             end do
!          end do
!       end do

!       write(*,*) 'IBM DIMENSIONS W E S N T'
!       write(*,*) ws,es,ss,ns,ts

!       if (ws /= es .or. ss /= ns) stop "IBM NOT ALLOWED IN BOUNDARIES"

!       allocate(inorm_west(3,ws),inorm_east(3,es),inorm_south(3,ss),inorm_north(3,ns),inorm_top(3,ts))
!       !   do k=1,kmax
!       !reuse those ints
!       ws = 1 !west
!       es = 1 !east
!       ss = 1 !south
!       ns = 1 !north
!       ts = 1 !top

!       do k=2,kmax
!          do j=2,j1
!             do i=2,i1
!                !WEST
!                if (libm(i,j,k) .and. .not.libm(i-1,j,k)) then
!                   inorm_west(:,ws) = (/i,j,k/)
!                   ws = ws +1
!                endif
!                !EAST
!                if (libm(i,j,k) .and. .not.libm(i+1,j,k)) then
!                   inorm_east(:,es) = (/i,j,k/)
!                   es = es +1
!                endif
!                !SOUTH
!                if (libm(i,j,k) .and. .not.libm(i,j-1,k)) then
!                   inorm_south(:,ss) = (/i,j,k/)
!                   ss = ss +1
!                endif
!                !NORTH
!                if (libm(i,j,k) .and. .not.libm(i,j+1,k)) then
!                   inorm_north(:,ns) = (/i,j,k/)
!                   ns = ns +1
!                endif
!                !TOP
!                if (libm(i,j,k) .and. .not.libm(i,j,k+1)) then
!                   inorm_top(:,ts) = (/i,j,k/)
!                   ts = ts +1
!                endif
!             end do
!          end do
!       end do


!       write(*,*) 'Succesfully found normal layers in all directions'
!       write(*,*) 'WEST', inorm_west
!       write(*,*) 'EAST', inorm_east
!    end subroutine init_ibm


! subroutine with lnorm, but i dont want that actually
!    subroutine init_ibm()
!       use config, only: lapplyibm,ibm_input_file,ifinputibm
!       use modglobal, only: i1,ih,j1,jh,k1,kh,kmax,zm,zt !zm = zh omgg
!       real(real32), allocatable :: bc_height(:,:)     !< Height of immersed boundary at grid pos x,y
!       integer       :: i, j, k !, ierr,ii,jj,kk,n  !cstep , kmin
!       character(128) :: readstring

!       if (.not. lapplyibm) return
!       allocate(bc_height (i1,j1))
!       allocate(libm(2-ih:i1+ih,2-jh:j1+jh,k1))

!       allocate(lnorm_x (2-ih:i1+ih,2-jh:j1+jh,k1))
!       allocate(lnorm_y (2-ih:i1+ih,2-jh:j1+jh,k1))
!       allocate(lnorm_z (2-ih:i1+ih,2-jh:j1+jh,k1))
!       !! READ FILE
!       write(6,*) 'Reading inputfile in modibm'
!       open (ifinputibm,file=ibm_input_file)
!       do k=1,7
!          read (ifinputibm,'(a100)') readstring
!          !  read (ifinputibm,*) readstring
!          write (6,*) readstring
!       enddo

!       do j=j1,2,-1  !cstep
!          do i=2,i1
!             read(ifinputibm,'(F6.1)') bc_height(i,j)
!          enddo
!       end do
!       close(ifinputibm)

!       bc_height(1,:)=bc_height(i1,:)
!       bc_height(:,1)=bc_height(:,j1)
!       write(6,*) 'Succesfully read inputfile in modibm'
!       !   write(*,*) bc_height

!       do i=2,i1
!          do j=2,j1
!             do k=1,kmax
!                if (zt(k) <= bc_height(i,j)) then  !obstacle height is above mid point of vertical grid
!                   libm (i,j,k) = .true.
!                   !   TODO make use of ksfc, mainly to improve diffusion loop
!                   !   ksfc (i,j)   = k + 1   !half (flux) level
!                   !   if (ksfc(i,j).gt.kibm_maxl) then
!                   !      kibm_maxl = k
!                   !   endif
!                   !   write (6,*) 'libm',i+myidx*imax,j+myidy*jmax,i,j,k,libm(i,j,k),bc_height(i+myidx*imax,j+myidy*jmax),zh(ksfc(i,j))
!                endif
!             end do

!          end do
!       end do

!       deallocate (bc_height)

!       !   do k=1,kmax
!       do k=2,k1
!          do j=2,j1
!             do i=2,i1
!                !libm   ipos=i+myidx*imax
!                !libm   jpos=j+myidy*jmax
!                !libm   if (.not. (limmersed_boundary(ipos,jpos,k)==limmersed_boundary(ipos-1,jpos,k))) then

!                !cstep            0     X      X     X      0       ,building position
!                !cstep           i-2   i-1     i    i+1    i+2
!                !cstep  libm      F     T      T     T      F
!                !cstep  lnorm_x   F     T      F     F      T  , the true points refer to u-positions on the
!                !                                                grid box (on its left)

!                if (libm(i,j,k).neqv.libm(i-1,j,k)) then
!                   lnorm_x(i,j,k)=.true.  !cstep a wall at position i with its normal pointing in the x-direction
!                endif

!                if (libm(i,j,k).neqv.libm(i,j-1,k)) then
!                   lnorm_y(i,j,k)=.true.
!                endif

!                if (libm(i,j,k).neqv.libm(i,j,k-1)) then
!                   lnorm_z(i,j,k)=.true.
!                endif
!             end do
!          end do
!       end do

!       write(*,*) 'Succesfully found normal layers in all directions'
!    end subroutine init_ibm
