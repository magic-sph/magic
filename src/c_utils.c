#include <stdio.h>

void print_address(char *str,double _Complex *arr) {
  printf("address of %s is %p\n",str,arr);
}

void print_cache_info_void(char *str, void *thisptr) {
  unsigned long tag,index;
  unsigned long ptr;
  ptr=(unsigned long)thisptr;

  tag=ptr>>15;
  ptr=(unsigned long)thisptr-(tag<<15);
  index=ptr>>6;
  printf("%20s: tag = %15lu, set index = %3lu\n",str,tag,index);
}

void print_cache_info_dcmplx(char *str, double _Complex *arr) {
  print_cache_info_void(str,(void*)arr);
}

void print_cache_info_dreal(char *str, double *arr) {
  print_cache_info_void(str,(void*)arr);
}

void print_cache_info_integer(char *str, int *arr) {
  print_cache_info_void(str,(void*)arr);
}

