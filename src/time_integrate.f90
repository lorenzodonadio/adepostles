module time_integrate

   implicit none
   integer,parameter :: rkmethod = 1 ! 1,2,3,4
contains
   subroutine time_step()
      use config, only: ladaptivedt
      use modglobal, only: simtime,dt,tres,rsts
      use modfields, only: cm,c0,cp
      !Time integration is the last step
      if (ladaptivedt) call adaptative_dt
      simtime = simtime + dt
      rsts = real(simtime)*tres
      c0 = c0 + dt*tres*cp
      cp = 0.
   end subroutine time_step

   subroutine adaptative_dt()
      use modglobal, only: dxi,dyi,dzi,dx2i,dy2i,dz2i,dzh,luniformz,i1,j1,k1,dt,tres
      use modfields, only: u0,v0,w0,ekh0
      real :: cfl,vnm
      real :: mydzi,mydz2i !< variables for the non uniform case, im not proud
      integer :: maxdt
      !!TODO make tis better for non uniform z
      if (luniformz) then
         mydzi = dzi
         mydz2i = dz2i
      else
         mydzi = size(dzh)/sum(dzh)
         mydz2i = mydzi**2
      endif

      cfl = 1./maxval((u0(2:i1,2:j1,2:k1)*dxi+v0(2:i1,2:j1,2:k1)*dyi+w0(2:i1,2:j1,2:k1)*mydzi))
      vnm = 0.5*maxval((ekh0(2:i1,2:j1,2:k1)*(dx2i + dy2i + mydz2i)))
      write (*,*) "cfl: ", cfl, "vnm: ", vnm
      maxdt = int(min(cfl,vnm)/tres)
      if (dt<maxdt) then
         dt = int(0.7*maxdt+0.3*dt)
      else
         dt = int(0.95*maxdt)
      endif

      write (*,*) "New Dt: ", dt
   end subroutine adaptative_dt

   !!! STABILITY ANALYSIS:
   !! ADVECTION CFL critetion:
   !! since we use explicit scheme Cmax = 1
   !! dt(u/dx+v/dy+w/dz) < 1 ==> dt < 1./(u/dx+v/dy+w/dz)
   !! DIFFUSION Von Neumann stability Analysis:
   !! since we use explicit scheme Cmax = 1
   !! dt < 1./(2*ekh*(1/dx^2 + 1/dy^2 + 1/dz^2)
end module time_integrate
