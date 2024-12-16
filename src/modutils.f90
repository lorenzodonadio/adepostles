module modutils
   implicit none

contains

   logical function str_ends_with(string, suffix)
      implicit none
      character(len=*), intent(in) :: string
      character(len=*), intent(in) :: suffix
      integer :: string_len, suffix_len

      string_len = len_trim(string)
      suffix_len = len_trim(suffix)

      if (string_len >= suffix_len) then
         str_ends_with = string(string_len - suffix_len + 1:string_len) == suffix
      else
         str_ends_with = .false.
      end if
   end function str_ends_with

end module modutils
