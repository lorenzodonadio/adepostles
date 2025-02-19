module modchem
   use modglobal, only: nsv
   use netcdf
   implicit none

   ! TODO decide if this is to come from user input or not. I think like this is betteer
   ! Molar weights g/mol
   real :: M_NO2 = 46.005 !> Molar weight of NO2
   real :: M_NO  = 30.006 !> Molar weight of NO
   real :: M_O3  = 47.997 !> Molar weight of O3

   ! Reaction constants
   ! JNO2​​ aprox 0.005−0.01, varies with sun... so maybe a way to input this also?
   real :: jNO2 = 0.005
   real :: kNO_03 = 3e-12
   logical :: lusechem = .false.
   logical :: lnullcycle = .false.
   ! Indices of tracers
   integer :: io3, ino2, ino
contains

   subroutine init_chem
      use modtracer, only: tracers
      use config, only: config_filename,ifconfig

      integer :: ios
      namelist /CHEMISTRY/ lusechem, lnullcycle, jNO2, kNO_03

      open(unit=ifconfig, file=config_filename, status='old', action='read', iostat=ios)
      read(ifconfig, nml=CHEMISTRY, iostat=ios)

      if (.not.lusechem) return

      if (lnullcycle) then
         ino2 = tracers%find_index_by_name('no2')
         ino = tracers%find_index_by_name('no')
         io3 = tracers%find_index_by_name('o3')

         if (ino2 == -1) stop "no2 needed for nullcycle chemistry"
         if (ino == -1)  stop "no needed for nullcycle chemistry"
         if (io3 == -1)  stop "o3 needed for nullcycle chemistry"
      endif

      write(*,*) 'Initialized chemistry: ', lusechem
   endsubroutine

   subroutine apply_chem
      if (.not.lusechem) return

      if (lnullcycle) then
         call nullcycle
      endif
   end subroutine apply_chem

   subroutine nullcycle
      use modglobal,  only : dt
      use modtracer,  only : c0
      use time_integrate, only: rkstep
      use modglobal, only:i1,ih,j1,jh,k1
      implicit none
      real , dimension(2-ih:i1+ih,2-jh:j1+jh,k1) :: tNO, tNO2, tO3,den
      real :: rdts
      ! real, dimension(ib-ihc:ie+ihc,jb-jhc:je+jhc,kb:ke+khc)     :: tNO, tNO2, tO3
      rdts = dt*1e-3 ! dt is in miliseconds so we convert it to seconds here

      if (.not. rkstep == 1)  return

      ! convert into mol/m^3
      tNO = c0(:,:,:,ino)/M_NO
      tNO2= c0(:,:,:,ino2)/M_NO2
      tO3 = c0(:,:,:,io3)/M_O3
      den =  1. + ( ( tNO + tO3 ) * kNO_03 + jNO2 ) * rdts !Denominator
      ! fully implicit chemistry
      c0(:,:,:,ino) = M_NO * ( tNO + ( rdts * (- kNO_03 * tNO * tO3 + jNO2 * tNO2) ) / den)

      c0(:,:,:,ino2) = M_NO2 * ( tNO2 - ( rdts * (- kNO_03 * tNO * tO3 + jNO2 * tNO2) ) / den)

      c0(:,:,:,io3) = M_O3 * ( tO3 + ( rdts * (- kNO_03 * tNO * tO3 + jNO2 * tNO2) ) / den)

      ! alternative method in Zhong 2017 eqn. 7, analytical solution!

   end subroutine nullcycle

endmodule modchem


! module modchem
!    use modglobal, only: nsv
!    use netcdf
!    implicit none

!    ! Molar weights g/mol
!    real :: M_NO2 = 46.005
!    real :: M_NO  = 30.006
!    real :: M_O3  = 47.997
!    real :: M_CO  = 28.010
!    real :: M_CO2 = 44.009
!    real :: M_VOC = 72.000 ! Generic VOC molar weight
!    real :: M_HO  = 17.008 ! Hydroxyl radical (OH)
   
!    ! Reaction constants
!    real :: jNO2 = 0.005
!    real :: kNO_O3 = 3e-12
!    real :: kCO_OH = 2.4e-13 ! CO + OH -> CO2
!    real :: kVOC_OH = 1.0e-11 ! VOC + OH -> Products

!    logical :: lusechem = .false.
!    logical :: lnullcycle = .false.
!    logical :: lCO_OH = .false.
!    logical :: lVOC_OH = .false.

!    ! Indices of tracers
!    integer :: io3, ino2, ino, ico, ico2, ivoc, iho

! contains

!    subroutine init_chem
!       use modtracer, only: tracers
!       use config, only: config_filename, ifconfig
!       integer :: ios
!       namelist /CHEMISTRY/ lusechem, lnullcycle, lCO_OH, lVOC_OH, jNO2, kNO_O3, kCO_OH, kVOC_OH

!       open(unit=ifconfig, file=config_filename, status='old', action='read', iostat=ios)
!       read(ifconfig, nml=CHEMISTRY, iostat=ios)
!       if (.not. lusechem) return

!       ino2 = tracers%find_index_by_name('no2')
!       ino = tracers%find_index_by_name('no')
!       io3 = tracers%find_index_by_name('o3')
!       ico = tracers%find_index_by_name('co')
!       ico2 = tracers%find_index_by_name('co2')
!       ivoc = tracers%find_index_by_name('voc')
!       iho = tracers%find_index_by_name('ho')

!       if (ino2 == -1 .or. ino == -1 .or. io3 == -1) stop "NOx species needed for chemistry"
!       if (lCO_OH .and. (ico == -1 .or. ico2 == -1 .or. iho == -1)) stop "CO-OH chemistry requires CO, CO2, and OH"
!       if (lVOC_OH .and. (ivoc == -1 .or. iho == -1)) stop "VOC-OH chemistry requires VOC and OH"
      
!       write(*,*) 'Initialized chemistry: ', lusechem
!    end subroutine init_chem

!    subroutine apply_chem
!       if (.not. lusechem) return

!       if (lnullcycle) call nullcycle
!       if (lCO_OH) call co_oh_reaction
!       if (lVOC_OH) call voc_oh_reaction
!    end subroutine apply_chem

!    subroutine co_oh_reaction
!       use modglobal, only: dt
!       use modtracer, only: c0
!       use time_integrate, only: rkstep
!       implicit none
!       real, dimension(:,:,:), allocatable :: tCO, tCO2, tOH, den
!       real :: rdts

!       if (rkstep /= 1) return
!       rdts = dt * 1e-3

!       allocate(tCO, tCO2, tOH, den)
!       tCO = c0(:,:,:,ico) / M_CO
!       tCO2 = c0(:,:,:,ico2) / M_CO2
!       tOH = c0(:,:,:,iho) / M_HO

!       den = 1.0 + kCO_OH * tCO * tOH * rdts

!       c0(:,:,:,ico) = M_CO * (tCO - kCO_OH * tCO * tOH * rdts / den)
!       c0(:,:,:,ico2) = M_CO2 * (tCO2 + kCO_OH * tCO * tOH * rdts / den)

!       deallocate(tCO, tCO2, tOH, den)
!    end subroutine co_oh_reaction

!    subroutine voc_oh_reaction
!       use modglobal, only: dt
!       use modtracer, only: c0
!       use time_integrate, only: rkstep
!       implicit none
!       real, dimension(:,:,:), allocatable :: tVOC, tOH, den
!       real :: rdts

!       if (rkstep /= 1) return
!       rdts = dt * 1e-3

!       allocate(tVOC, tOH, den)
!       tVOC = c0(:,:,:,ivoc) / M_VOC
!       tOH = c0(:,:,:,iho) / M_HO

!       den = 1.0 + kVOC_OH * tVOC * tOH * rdts

!       c0(:,:,:,ivoc) = M_VOC * (tVOC - kVOC_OH * tVOC * tOH * rdts / den)

!       deallocate(tVOC, tOH, den)
!    end subroutine voc_oh_reaction

! endmodule modchem
