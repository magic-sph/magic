MODULE cutils
  USE, intrinsic :: iso_c_binding
  IMPLICIT NONE

  INTERFACE 
     SUBROUTINE print_address(str,arr) BIND(c)
       USE,intrinsic :: ISO_C_BINDING
       CHARACTER,DIMENSION(*),INTENT(IN) :: str
       COMPLEX(C_DOUBLE_COMPLEX),DIMENSION(*) :: arr
     END SUBROUTINE print_address
  END INTERFACE


  INTERFACE 
     !void print_cache_info_dcmplx(char *str, DOUBLE _Complex *arr) {
     SUBROUTINE print_cache_info_dcmplx(str,arr) BIND(c)
       USE,intrinsic :: ISO_C_BINDING
       CHARACTER,DIMENSION(*),INTENT(IN) :: str
       COMPLEX(C_DOUBLE_COMPLEX),DIMENSION(*) :: arr
     END SUBROUTINE print_cache_info_dcmplx
     !void print_cache_info_dreal(char *str, double *arr) {
     SUBROUTINE print_cache_info_dreal(str,arr) BIND(c)
       USE,intrinsic :: ISO_C_BINDING
       CHARACTER,DIMENSION(*),INTENT(IN) :: str
       real(C_DOUBLE),DIMENSION(*) :: arr
     END SUBROUTINE print_cache_info_dreal
     SUBROUTINE print_cache_info_integer(str,ptr) BIND(c)
       USE,intrinsic :: ISO_C_BINDING
       CHARACTER,DIMENSION(*),INTENT(IN) :: str
       integer(C_INT) :: ptr
     END SUBROUTINE print_cache_info_integer
  END INTERFACE

END MODULE cutils
