file(GLOB srcfiles "*.f90")
#set(red_sources ${srcfiles})
list(REMOVE_ITEM srcfiles "${SRC}/algebra_lapack.f90")
list(REMOVE_ITEM srcfiles "${SRC}/algebra.f90")
list(REMOVE_ITEM srcfiles "${SRC}/algebra_hipfort.f90")
list(REMOVE_ITEM srcfiles "${SRC}/cosine_transform_odd.f90")
list(REMOVE_ITEM srcfiles "${SRC}/dct_fftw.f90")
list(REMOVE_ITEM srcfiles "${SRC}/dct_hipfort_hipfftw.f90")
list(REMOVE_ITEM srcfiles "${SRC}/hipfort_hipfftw.f90")
list(REMOVE_ITEM srcfiles "${SRC}/hipfort_check.f90")
list(REMOVE_ITEM srcfiles "${SRC}/hipfort_hipblas.f90")
list(REMOVE_ITEM srcfiles "${SRC}/hipfort_hipblas_auxiliary.f90")
list(REMOVE_ITEM srcfiles "${SRC}/hipfort_hipblas_enums.f90")
list(REMOVE_ITEM srcfiles "${SRC}/hipfort_hipsolver.f90")
list(REMOVE_ITEM srcfiles "${SRC}/hipfort_hipsolver_enums.f90")
list(REMOVE_ITEM srcfiles "${SRC}/algebra_hipfort.f90")
list(REMOVE_ITEM srcfiles "${SRC}/hipfort.f90")
list(REMOVE_ITEM srcfiles "${SRC}/hipfort_auxiliary.f90")
list(REMOVE_ITEM srcfiles "${SRC}/hipfort_cuda_errors.f90")
list(REMOVE_ITEM srcfiles "${SRC}/hipfort_enums.f90")
list(REMOVE_ITEM srcfiles "${SRC}/hipfort_hipmalloc.f90")
list(REMOVE_ITEM srcfiles "${SRC}/hipfort_hipmemcpy.f90")
list(REMOVE_ITEM srcfiles "${SRC}/hipfort_types.f90")
list(REMOVE_ITEM srcfiles "${SRC}/hipfort_hipsparse.f90")
list(REMOVE_ITEM srcfiles "${SRC}/hipfort.f90")

#set(red_sources ${srcfiles})
list(REMOVE_ITEM srcfiles "${SRC}/fft_fftw.f90")
list(REMOVE_ITEM srcfiles "${SRC}/fft.f90")
list(REMOVE_ITEM srcfiles "${SRC}/magic.f90")

set_source_files_properties(mpimod.f90 PROPERTIES COMPILE_FLAGS ${FFLAG_NO_STD})

if(USE_MPI MATCHES no)
   # Dirty hack to make the code switching endianess at reading stage
   # Allows to read MPI-IO checkpoints when running the serial version of the
   # code
   set_source_files_properties(readCheckPoints.f90 PROPERTIES COMPILE_FLAGS ${FFLAG_NO_STD})
endif()

##################
# File selection #
##################

#Hipfort FFTW/BLAS/LAPACK/SPARSE
if(USE_GPU MATCHES yes AND HIPFORT MATCHES HIPFORT-NOTFOUND)
   set(hipAlgebra hipfort_enums.f90 hipfort_hipfftw.f90 hipfort_check.f90 hipfort_hipblas.f90 hipfort_hipblas_auxiliary.f90 hipfort_hipblas_enums.f90 hipfort_hipsolver.f90 hipfort_hipsolver_enums.f90 hipfort_hipsparse_enums.f90 hipfort_hipsparse.f90 hipfort.f90)
   set(hip hipfort_types.f90 hipfort_auxiliary.f90 hipfort_cuda_errors.f90 hipfort_check.f90 hipfort.f90 hipfort_hipmalloc.f90 hipfort_hipmemcpy.f90)
endif()

# Lapack
if(USE_LAPACKLIB MATCHES JW)
   if(USE_GPU MATCHES yes)
      if(NOT HIPFORT MATCHES HIPFORT-NOTFOUND)
         set(lapack algebra.f90 algebra_hipfort.f90)
      else()
         set(lapack algebra.f90 ${hipAlgebra} algebra_hipfort.f90)
      endif()
   else()
      set(lapack algebra.f90)
   endif(USE_GPU)
elseif(USE_LAPACKLIB MATCHES LAPACK OR USE_LAPACKLIB MATCHES MKL OR USE_LAPACKLIB MATCHES LIBFLAME OR USE_LAPACKLIB MATCHES LIBSCI)
   if(USE_GPU MATCHES yes)
      if(NOT HIPFORT MATCHES HIPFORT-NOTFOUND)
         set(lapack algebra_lapack.f90 algebra_hipfort.f90)
      else()
         set(lapack algebra_lapack.f90 ${hipAlgebra} algebra_hipfort.f90)
      endif()
   else()
      set(lapack algebra_lapack.f90)
   endif()
endif()

#DCT
if(USE_DCTLIB MATCHES JW)
   if(USE_GPU MATCHES yes)
      if(NOT HIPFORT MATCHES HIPFORT-NOTFOUND)
         set(dct cosine_transform_odd.f90 dct_hipfort_hipfftw.f90)
      else()
         set(dct cosine_transform_odd.f90 ${hipAlgebra} dct_hipfort_hipfftw.f90)
      endif()
   else()
      set(dct cosine_transform_odd.f90)
   endif()
else()
   if(USE_GPU MATCHES yes)
      if(NOT HIPFORT MATCHES HIPFORT-NOTFOUND)
         set(dct dct_fftw.f90 dct_hipfort_hipfftw.f90)
      else()
         set(dct dct_fftw.f90 ${hipAlgebra} dct_hipfort_hipfftw.f90)
      endif()
   else()
      set(dct dct_fftw.f90)
   endif()
endif()

# FFT
if(USE_FFTLIB MATCHES JW)
   set(fft fft.f90)
elseif(USE_FFTLIB MATCHES MKL)
   set(fft fft_fftw.f90)
elseif(USE_FFTLIB MATCHES FFTW)
   set(fft fft_fftw.f90)
endif()

#SHTNS
if(USE_SHTNS MATCHES yes AND NOT SHTNS MATCHES SHTNS-NOTFOUND)
   list(REMOVE_ITEM srcfiles "${SRC}/shtransforms.f90")
   list(REMOVE_ITEM srcfiles "${SRC}/sht_native.f90")
   if(USE_GPU MATCHES yes)
      if(NOT HIPFORT MATCHES HIPFORT-NOTFOUND)
         set(whichsht shtns.f90)
      else()
         set(whichsht ${hip} shtns.f90)
      endif()
   else()
      set(whichsht shtns.f90)
   endif()
else()
   list(REMOVE_ITEM srcfiles "${SRC}/shtns.f90")
   set(whichsht sht_native.f90)
endif()

add_executable(magic.exe magic.f90 ${srcfiles} ${fft} ${dct} ${lapack} ${whichsht})
set_property(TARGET ${EXEC} PROPERTY LINKER_LANGUAGE Fortran)

########
# Link #
########

# Link MPI
if(USE_MPI MATCHES yes)
   if(MPI_Fortran_FOUND AND NOT MPI_WRAPPER_FOUND)
      target_link_libraries(${EXEC} ${MPI_Fortran_LIBRARIES})
   endif()
endif()

# Link MKL
if(USE_LAPACKLIB MATCHES MKL AND USE_FFTLIB MATCHES MKL)
   add_library(lapack95 ${MKLROOT}/include/lapack.f90)
   target_link_libraries(lapack95 ${MKLCORE} ${MKLLAPACK} ${MKLSEQ} ${MKL64})
   target_link_libraries(${EXEC} lapack95)
elseif(USE_LAPACKLIB MATCHES MKL AND USE_FFTLIB MATCHES JW)
   add_library(lapack95 ${MKLROOT}/include/lapack.f90)
   target_link_libraries(lapack95 ${MKLCORE} ${MKLLAPACK} ${MKLSEQ} ${MKL64})
   target_link_libraries(${EXEC} lapack95)
elseif(USE_LAPACKLIB MATCHES MKL AND USE_FFTLIB MATCHES FFTW)
   add_library(lapack95 ${MKLROOT}/include/lapack.f90)
   target_link_libraries(lapack95 ${MKLCORE} ${MKLLAPACK} ${MKLSEQ} ${MKL64})
   target_link_libraries(${EXEC} lapack95)
elseif(USE_LAPACKLIB MATCHES JW AND USE_FFTLIB MATCHES MKL)
   target_link_libraries(mkl_fft ${MKLCORE} ${MKLSEQ} ${MKL64})
   target_link_libraries(${EXEC} mkl_fft)
elseif(USE_LAPACKLIB MATCHES LIBFLAME AND USE_FFTLIB MATCHES MKL)
   target_link_libraries(mkl_fft ${MKLCORE} ${MKLSEQ} ${MKL64})
   target_link_libraries(${EXEC} mkl_fft)
endif()

# Link LAPACK
if(USE_LAPACKLIB MATCHES LAPACK)
   target_link_libraries(${EXEC} ${LAPACK_LIBRARIES})
endif()

# Link HIPFORT
if(USE_GPU MATCHES yes)
   if(NOT HIPFORT MATCHES HIPFORT-NOTFOUND)
      target_link_libraries(${EXEC} -lhipfft)
      target_link_libraries(${EXEC} -lrocfft)
      target_link_libraries(${EXEC} -lhipsolver)
      target_link_libraries(${EXEC} -lrocsolver)
      target_link_libraries(${EXEC} -lhipblas)
      target_link_libraries(${EXEC} ${HIPFORT})
      #target_link_libraries(${EXEC} -lamdhip64)
   else()
      target_link_libraries(${EXEC} -lhipfft)
      target_link_libraries(${EXEC} -lhipsolver)
      target_link_libraries(${EXEC} -lhipblas)
      #target_link_libraries(${EXEC} -lamdhip64)
   endif()
else()
   target_link_libraries(${EXEC} -lhipfft)	   
endif()

# Link LIBFLAME and BLIS
if(USE_LAPACKLIB MATCHES LIBFLAME)
   target_link_libraries (${EXEC} ${FLAME} ${BLIS})
#elseif(USE_LAPACKLIB MATCHES LIBSCI)
   #target_link_libraries (${EXEC} ${SCI})
endif()

if(USE_SHTNS MATCHES yes AND NOT SHTNS MATCHES SHTNS-NOTFOUND)
   # check shtns backend for ffts (mkl or fftw)
   FILE(WRITE "${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeTmp/testSHTNS.f90"
   "
   program TestFFTW
    call shtns_verbose(0)
   end program TestFFTW
   ")

   try_run(SHTNS_RUN SHTNS_LINK ${CMAKE_BINARY_DIR}
      ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeTmp/testSHTNS.f90
      COMPILE_DEFINITIONS ${CMAKE_REQUIRED_DEFINITIONS}
      LINK_LIBRARIES ${SHTNS} ${MKLCORE} ${MKLSEQ} ${MKL64}
      COMPILE_OUTPUT_VARIABLE OUTPUT)

   if(SHTNS_LINK)
      set(SHTNS_DEPS ${MKLCORE} ${MKLSEQ} ${MKL64})
   endif()

   if(NOT SHTNS_LINK)
      set(SAFE_CMAKE_LIBS "${LIBS}")
      file(WRITE "${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeTmp/testSHTNS.f90"
      "
      program TestFFTW
      call shtns_verbose(0)
      end program TestFFTW
      ")
      try_run(SHTNS_RUN SHTNS_LINK ${CMAKE_BINARY_DIR}
         ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeTmp/testSHTNS.f90
         COMPILE_DEFINITIONS ${CMAKE_REQUIRED_DEFINITIONS}
         LINK_LIBRARIES ${SHTNS} ${FFTW3_OMP} ${FFTW3}
         COMPILE_OUTPUT_VARIABLE OUTPUT)
      if(SHTNS_LINK)
         set(SHTNS_DEPS ${FFTW3_OMP} ${FFTW3})
      endif(SHTNS_LINK)
   endif(NOT SHTNS_LINK)

   set(LIBS "${LIBS} ${SHTNS_DEPS}")
   target_link_libraries(${EXEC} ${SHTNS})
   target_link_libraries(${EXEC} ${SHTNS_DEPS})
else()
   if(USE_FFTLIB MATCHES FFTW)
      if(USE_OMP)
         target_link_libraries(${EXEC} ${FFTW3_OMP} ${FFTW3})
      else()
         target_link_libraries(${EXEC} ${FFTW3})
      endif()
   endif()
endif()
