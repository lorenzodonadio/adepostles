module time_integrate

   implicit none

contains
   subroutine time_step()
      use modglobal, only: simtime,dt,tres,rsts
      !Time integration is the last step
      simtime = simtime + dt
      rsts = real(simtime)*tres
   end subroutine time_step
end module time_integrate
