module modfields
   use iso_fortran_env
   use modglobal, only:i1,ih,j1,jh,k1
   use config, only: field_load_chunk_size

   implicit none
   real(real32), allocatable :: ekh(:,:,:,:) !< k-coefficient for eddy diffusivity ekh(xt,yt,zt,time) m2/s
   real(real32), allocatable :: u(:,:,:,:)   !< u(xm,yt,zt,time) m/s
   real(real32), allocatable :: v(:,:,:,:)   !< v(xt,ym,zt,time) m/s
   real(real32), allocatable :: w(:,:,:,:)   !< w(xt,yt,zm,time) m/s
   ! at different time stamps
   real(real32), allocatable :: c0(:,:,:)   !< Concentration
   real(real32), allocatable :: cp(:,:,:)   !< Concentration
   real(real32), allocatable :: cm(:,:,:)   !< Concentration

   real(real32), allocatable :: ekh0(:,:,:) !< ekhm at time t
   real(real32), allocatable :: u0(:,:,:)   !< u at time t m/s
   real(real32), allocatable :: v0(:,:,:)   !< v at time t m/s
   real(real32), allocatable :: w0(:,:,:)   !< w at time t m/s

   real(real32), allocatable :: ekhp(:,:,:) !< ekhm at time t+1
   real(real32), allocatable :: up(:,:,:)   !< u at time t+1 m/s
   real(real32), allocatable :: vp(:,:,:)   !< v at time t+1 m/s
   real(real32), allocatable :: wp(:,:,:)   !< w at time t+1 m/s

   real(real32), allocatable :: ekhm(:,:,:) !< ekhm at time t-1
   real(real32), allocatable :: um(:,:,:)   !< u at time t-1m/s
   real(real32), allocatable :: vm(:,:,:)   !< v at time t-1 m/s
   real(real32), allocatable :: wm(:,:,:)   !< w at time t-1m/s

contains

   subroutine allocate_fields()
      !< k-coefficient for eddy diffusivity ekh(xt,yt,zt,time) m2/s
      !< u(xm,yt,zt,time) m/s
      !< v(xt,ym,zt,time) m/s
      !< w(xt,yt,zm,time) m/s
      allocate(ekh  (2-ih:i1+ih,2-jh:j1+jh,k1,field_load_chunk_size))
      allocate(u    (2-ih:i1+ih,2-jh:j1+jh,k1,field_load_chunk_size))
      allocate(v    (2-ih:i1+ih,2-jh:j1+jh,k1,field_load_chunk_size))
      allocate(w    (2-ih:i1+ih,2-jh:j1+jh,k1,field_load_chunk_size))

      ! TODO make this parallelizable, so that we can handle multiple ADE solutions with the same u,v,w,ekh,rho, thats the goal
      ! maybe we just add a fourth dimension and thats it?
      allocate(c0   (2-ih:i1+ih,2-jh:j1+jh,k1))
      allocate(cp   (2-ih:i1+ih,2-jh:j1+jh,k1))
      allocate(cm   (2-ih:i1+ih,2-jh:j1+jh,k1))

      allocate(ekh0 (2-ih:i1+ih,2-jh:j1+jh,k1))
      allocate(u0   (2-ih:i1+ih,2-jh:j1+jh,k1))
      allocate(v0   (2-ih:i1+ih,2-jh:j1+jh,k1))
      allocate(w0   (2-ih:i1+ih,2-jh:j1+jh,k1))

      allocate(ekh0 (2-ih:i1+ih,2-jh:j1+jh,k1))
      allocate(u0   (2-ih:i1+ih,2-jh:j1+jh,k1))
      allocate(v0   (2-ih:i1+ih,2-jh:j1+jh,k1))
      allocate(w0   (2-ih:i1+ih,2-jh:j1+jh,k1))

      allocate(ekh0 (2-ih:i1+ih,2-jh:j1+jh,k1))
      allocate(u0   (2-ih:i1+ih,2-jh:j1+jh,k1))
      allocate(v0   (2-ih:i1+ih,2-jh:j1+jh,k1))
      allocate(w0   (2-ih:i1+ih,2-jh:j1+jh,k1))

   end subroutine allocate_fields

end module modfields
