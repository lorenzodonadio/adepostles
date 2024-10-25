module time_integrate
   implicit none
   integer, save :: rkstep = 1   ! Current RK substep

contains
   subroutine time_step()
      use config, only: ladaptivedt,rkmethod
      use modglobal, only: simtime, dt, tres, rsts

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
         simtime = simtime + dt     ! microseconds
         rsts = real(simtime)*tres  ! seconds
      endif
   end subroutine time_step

   subroutine euler_step()
      use modfields, only: c0, cp
      use modglobal, only: dt, tres

      ! Euler method is single step
      c0 = c0 + dt*tres*cp
      cp = 0.
      rkstep = 1
   end subroutine euler_step

   subroutine rk2_step()
      use modfields, only: cm, c0, cp
      use modglobal, only: dt, tres

      select case(rkstep)
       case(1)
         ! First step: Store initial state in cm and update c0
         cm = c0  ! Store initial state
         c0 = c0 + dt*tres*cp
         cp = 0.
         rkstep = 2

       case(2)
         ! Second step: Complete RK2 update using cm as initial state
         c0 = cm + 0.5*dt*tres*(cp + cp)  ! cp here contains k2
         cp = 0.
         rkstep = 1
      end select
   end subroutine rk2_step

   subroutine rk3_step()
      use modfields, only: cm, c0, cp
      use modglobal, only: dt, tres
      real :: rk3coef

      rk3coef = dt*tres/(4. - real(rkstep))

      if (rkstep /= 3) then
         ! Steps 1 and 2
         c0 = cm + rk3coef * cp
         cp = 0.
         rkstep = rkstep + 1
      else
         ! Step 3 - store result in both c0 and cm
         cm = cm + rk3coef * cp
         c0 = cm
         cp = 0.
         rkstep = 1
      endif
   end subroutine rk3_step

   subroutine rk4_step()
      use modfields, only: cm, c0, cp,rk4_acc
      use modglobal, only: dt, tres

      select case(rkstep)
       case(1)
         ! First step
         cm = c0           ! Store initial state
         rk4_acc = cp     ! Store k1
         c0 = cm + 0.5*dt*tres*cp  ! Prepare for k2
         cp = 0.
         rkstep = 2

       case(2)
         ! Second step - now cp contains k2
         rk4_acc = rk4_acc + 2.0*cp  ! Accumulate 2*k2
         c0 = cm + 0.5*dt*tres*cp  ! Prepare for k3
         cp = 0.
         rkstep = 3

       case(3)
         ! Third step - now cp contains k3
         rk4_acc = rk4_acc + 2.0*cp  ! Accumulate 2*k3
         c0 = cm + dt*tres*cp  ! Prepare for k4
         cp = 0.
         rkstep = 4

       case(4)
         ! Fourth step - now cp contains k4
         ! Final update using accumulated terms
         c0 = cm + (dt*tres/6.0)*(rk4_acc + cp)
         cp = 0.
         rkstep = 1
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
         cfl_limit = 1.0
         vn_limit = 0.5
       case(2)  ! RK2
         cfl_limit = 1.0
         vn_limit = 1.0
       case(3)  ! RK3
         cfl_limit = 1.73  ! sqrt(3) and a bit less
         vn_limit = 1.36
       case(4)  ! RK4
         cfl_limit = 2.82  ! sqrt(8) and a bit less
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
         dt = int(0.7*maxdt+0.3*dt)
      else
         dt = int(0.95*maxdt)
      endif

      write (*,*) "cfl: ", cfl, "vnm: ", vnm , "New dt: ", dt*tres , 's'
   end subroutine adaptative_dt
end module time_integrate
