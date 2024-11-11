module time_integrate
   use iso_fortran_env, only: real32
   implicit none
   integer, save :: rkstep = 1   ! Current RK substep

   real(real32), allocatable :: rk1(:,:,:,:)
   real(real32), allocatable :: rk2(:,:,:,:)
   real(real32), allocatable :: rk3(:,:,:,:)

contains

   subroutine allocate_k()
      use modglobal, only:i1,ih,j1,jh,k1,nsv
      use config, only: field_load_chunk_size,rkmethod

      if (rkmethod > 2) allocate(rk1(2-ih:i1+ih,2-jh:j1+jh,k1,nsv)) !2nd order or higher may need k1
      if (rkmethod > 2) allocate(rk2(2-ih:i1+ih,2-jh:j1+jh,k1,nsv)) !3nd order or higher may need k2

      if (rkmethod > 3) allocate(rk3(2-ih:i1+ih,2-jh:j1+jh,k1,nsv)) !4nd order or higher may need k3


   end subroutine allocate_k

   subroutine increase_simtime(in_dt)
      use modglobal, only: simtime, tres, rsts
      integer, intent(in) :: in_dt

      simtime = simtime + in_dt     ! microseconds
      rsts = real(simtime)*tres  ! seconds

   end subroutine increase_simtime

   subroutine time_step()
      use config, only: ladaptivedt,rkmethod
      use modglobal, only: simtime, dt, tres, rsts
      use modtracer, only:  c0
      ! Check if we need adaptive timestep (only at start of RK steps)
      if (rkstep == 1 .and. ladaptivedt) call adaptative_dt

      select case(rkmethod)
       case(1)
         call euler_step()
       case(2)
         call rk2_step()
       case(3)
         call rk3_step()
       case(4)
         call rk4_step()
      end select

      ! Update time only after completing all substeps
      if (rkstep == 1) then
         ! ! Optional: Add positivity preservation
         where (c0 < 0.0)
            c0 = 0.0
         end where
         ! actually every rk method should increase the dt in their respective steps
         ! call increase_simtime(dt)
         ! simtime = simtime + dt     ! microseconds
         ! rsts = real(simtime)*tres  ! seconds
      endif
   end subroutine time_step

   subroutine euler_step()
      use modtracer, only: c0, cp
      use modglobal, only: dt, tres

      ! Euler method is single step
      c0 = c0 + dt*tres*cp
      cp = 0.
      call increase_simtime(dt)
   end subroutine euler_step

   subroutine rk2_step()
      use modtracer, only: cm, c0, cp
      use modglobal, only: dt, tres

      select case(rkstep)
       case(1)
         ! First step: Store initial state in cm and update c0
         cm = c0  ! Store initial state
         c0 = c0 + dt*tres*cp
         ! cp = 0.
         rkstep = 2
         call increase_simtime(dt)

       case(2)
         ! Second step: Complete RK2 update using cm as initial state
         c0 = cm + dt*tres*cp ! cp here contains k2
         cp = 0.
         rkstep = 1
      end select
   end subroutine rk2_step

   subroutine rk3_step()
      !Ralston method
      ! 0   |   0    0     0
      ! 1/2 |   1/2  0     0
      ! 3/4 |   0    3/4   0
      ! __________________
      !     |   2/9  1/3   4/9

      use modtracer, only: cm, c0, cp
      use modglobal, only: dt, tres
      real :: rdt
      rdt = dt*tres

      select case(rkstep)
       case(1)
         ! First stage
         cm = c0            ! Store initial state
         rk1 = cp
         c0 = cm + 0.5*rdt*rk1
         rkstep = 2
         call increase_simtime(int(dt/2))
       case(2)
         rk2 = cp
         c0 = cm + 0.75*rdt*rk2
         rkstep = 3
         call increase_simtime(int(dt/4))
       case(3)
         ! Third stage
         ! k3 = cp
         ! c0 = cm + 2/9 * k1 + (1/3)*k3 + 4/9 * k3
         c0 = cm + 0.2222222222*rdt*rk1 + 0.3333333333*rdt*rk2 + 0.4444444444*rdt*cp
         ! c0 = cm + (2.0/9.0)*rdt*rk1 + (1.0/3.0)*rdt*rk2 + (4.0/9.0)*rdt*cp
         cp = 0.
         rkstep = 1
         call increase_simtime(int(dt/4))
      end select
   end subroutine rk3_step

   subroutine rk4_step()
      ! 3/8 Runge-Kutta Method (Fourth-Order)
      ! Butcher tableau:
      !   0   |   0       0      0      0
      !  1/3  |   1/3     0      0      0
      !  2/3  |  -1/3     1      0      0
      !   1   |   1      -1      1      0
      ! ------------------------------------
      !       |  1/8     3/8    3/8    1/8

      use modtracer, only: cm, c0, cp
      use modglobal, only: dt, tres
      real :: rdt
      integer :: dt_div_3
      dt_div_3 = int(dt/3)
      rdt = dt * tres

      select case(rkstep)
       case(1)
         cm = c0                   ! Store the initial state
         rk1 = cp                  ! k1 = cp
         c0 = cm + 0.3333333333 * rdt * rk1
         rkstep = 2
         call increase_simtime(dt_div_3)

       case(2)
         rk2 = cp                  ! k2 = cp after recalculation
         c0 = cm - 0.3333333333* rdt * rk1 + rdt * rk2
         rkstep = 3
         call increase_simtime(dt_div_3)

       case(3)
         rk3 = cp                  ! k3 = cp after recalculation
         c0 = cm + rdt * rk1 - rdt * rk2 + rdt * rk3
         rkstep = 4
         call increase_simtime(dt_div_3)

       case(4)
         ! k4 = cp after recalculation
         c0 = cm + 0.125 * rdt * rk1 + 0.375 * rdt * rk2 + 0.375 * rdt * rk3 + 0.125 * rdt * cp
         cp = 0.                   ! Final reset of cp after 3/8 RK4 step
         rkstep = 1                ! Reset rkstep for the next cycle

      end select
   end subroutine rk4_step

   subroutine adaptative_dt()
      use config, only: rkmethod
      use modglobal, only: dxi,dyi,dzi,dx2i,dy2i,dz2i,dzh,luniformz,i1,j1,k1,dt,tres
      use modfields, only: u0,v0,w0,ekh0
      real :: cfl,vnm,cfl_limit,vn_limit
      real :: mydzi,mydz2i
      integer :: maxdt

      ! Set stability limits based on RK method
      select case(rkmethod)
       case(1)  ! RK1
         ! cfl_limit = 1.0
         cfl_limit = 0.9
         vn_limit = 0.5
       case(2)  ! RK2
         cfl_limit = 0.95
         vn_limit = 1.0
       case(3)  ! RK3
         ! cfl_limit = 1.73  ! sqrt(3) and a bit less
         cfl_limit = 0.9
         vn_limit = 1.36
       case(4)  ! RK4
         ! cfl_limit = 2.  ! sqrt(8) = 2.82 and a bit less
         cfl_limit = 0.9  ! sqrt(8) and a bit less
         vn_limit = 2.78
      end select

      if (luniformz) then
         mydzi = dzi
         mydz2i = dz2i
      else
         mydzi = size(dzh)/sum(dzh)
         mydz2i = mydzi**2
      endif

      ! Calculate CFL and VN numbers (now multiplied by their limits)
      cfl = cfl_limit/maxval((u0(2:i1,2:j1,2:k1)*dxi+v0(2:i1,2:j1,2:k1)*dyi+w0(2:i1,2:j1,2:k1)*mydzi))
      vnm = vn_limit*maxval((ekh0(2:i1,2:j1,2:k1)*(dx2i + dy2i + mydz2i)))

      maxdt = int(min(cfl,vnm)/tres)
      if (dt<maxdt) then
         dt = int(0.85*maxdt+0.15*dt)
      else
         dt = maxdt
      endif

      ! write (*,*) "cfl: ", cfl, "vnm: ", vnm , "New dt: ", dt*tres , 's'
   end subroutine adaptative_dt
end module time_integrate
