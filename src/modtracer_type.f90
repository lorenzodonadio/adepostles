module modtracer_type

   use iso_fortran_env, only: real32

   implicit none

   type T_tracer
      ! Fixed tracer properties
      ! Tracer name
      character(len=16) :: tracname="tracname"
      ! Tracer long name
      character(len=64) :: traclong="dummy_long_name"
      ! Tracer unit
      character(len=16) :: unit="dummy_unit"
      ! Moleculare mass of tracer (g mol-1)
      real(real32)     :: molar_mass=-999.

      integer           :: trac_idx=-1
      ! Boolean if tracer is emitted
      logical           :: lemis=.false.
      ! Boolean if tracer is reactive
      logical           :: lreact=.false.
      ! Boolean if tracer is deposited
      logical           :: ldep=.false.
      ! Boolean if in A-gs
      logical           :: lags=.false.
      ! Boolean if in cloud microphysics
      logical           :: lmicro=.false.
      ! ! Static tracer properties:
      ! real :: diffusivity

   contains
      procedure :: print_properties => tracer_print_properties
   end type T_tracer

   ! New Type to manage tracers
   type T_tracer_collection
      type(T_tracer), allocatable :: data(:)
   contains
      procedure :: find_tracer_by_name
      procedure :: find_index_by_name
   end type T_tracer_collection

contains

   subroutine tracer_print_properties(self)

      class(T_tracer), intent(in) :: self

      write(*,*) "Tracer: ", self%tracname
      write(*,*) "  long name  : ", trim(self%traclong)
      write(*,*) "  unit       : ", trim(self%unit)
      write(*,*) "  molar_mass : ", self%molar_mass
      write(*,*) "  trac_idx      : ", self%trac_idx
      write(*,*) "  lemis      : ", self%lemis
      write(*,*) "  lreact     : ", self%lreact
      write(*,*) "  ldep       : ", self%ldep
      write(*,*) "  lags       : ", self%lags
      write(*,*) "  lmicro     : ", self%lmicro

   end subroutine tracer_print_properties

   function find_tracer_by_name(self, name) result(tracer)
      class(T_tracer_collection), intent(in) :: self
      character(len=*), intent(in) :: name
      type(T_tracer) :: tracer
      integer :: i

      do i = 1, size(self%data)
         if (trim(self%data(i)%tracname) == trim(name)) then
            tracer = self%data(i)
            return
         end if
      end do
   end function find_tracer_by_name

   function find_index_by_name(self, name) result(index)
      class(T_tracer_collection), intent(in) :: self
      character(len=*), intent(in) :: name
      type(T_tracer) :: tracer
      integer :: index
      tracer = self%find_tracer_by_name(name)
      index = tracer%trac_idx
   end function find_index_by_name
end module modtracer_type

