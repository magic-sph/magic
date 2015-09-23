module cutils

   use, intrinsic :: iso_c_binding
 
   implicit none
 
   interface 
      !-------------------------------------------------------
      subroutine print_address(str,arr) BIND(c)
 
         use, intrinsic :: ISO_C_BINDING
         character, intent(in) :: str(*)
         complex(C_DOUBLE_COMPLEX) :: arr(*)
 
      end subroutine print_address
      !-------------------------------------------------------
   end interface
 
 
   interface 
      !-------------------------------------------------------
      subroutine print_cache_info_dcmplx(str,arr) BIND(c)
 
         use, intrinsic :: ISO_C_BINDING
         character, intent(in) :: str(*)
         complex(C_DOUBLE_COMPLEX) :: arr(*)
 
      end subroutine print_cache_info_dcmplx
      !-------------------------------------------------------
      subroutine print_cache_info_dreal(str,arr) BIND(c)
 
         use, intrinsic :: ISO_C_BINDING
         character, intent(in) :: str(*)
         real(C_DOUBLE) :: arr(*)
 
      end subroutine print_cache_info_dreal
      !-------------------------------------------------------
      subroutine print_cache_info_integer(str,ptr) BIND(c)
 
        use, intrinsic :: ISO_C_BINDING
        character, intent(in) :: str(*)
        integer(C_INT) :: ptr
 
      end subroutine print_cache_info_integer
      !-------------------------------------------------------
   end interface

end module cutils
