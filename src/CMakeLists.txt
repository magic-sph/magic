file(GLOB srcfiles "*.f90")
#set(red_sources ${srcfiles})
list(REMOVE_ITEM srcfiles "${SRC}/algebra_lapack.f90")
list(REMOVE_ITEM srcfiles "${SRC}/algebra.f90")
list(REMOVE_ITEM srcfiles "${SRC}/cosine_transform_odd.f90")
list(REMOVE_ITEM srcfiles "${SRC}/dct_fftw.f90")
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

# Lapack
if(USE_LAPACKLIB MATCHES JW)
   set(lapack algebra.f90)
elseif(USE_LAPACKLIB MATCHES LAPACK OR USE_LAPACKLIB MATCHES MKL OR USE_LAPACKLIB MATCHES LIBFLAME OR USE_LAPACKLIB MATCHES LIBSCI)
   set(lapack algebra_lapack.f90)
endif()

#DCT
if(USE_DCTLIB MATCHES JW)
   set(dct cosine_transform_odd.f90)
else()
   set(dct dct_fftw.f90)
endif()

# FFT
if(USE_FFTLIB MATCHES JW)
   set(fft fft.f90)
elseif(USE_FFTLIB MATCHES MKL)
   set(fft fft_fftw.f90)
elseif(USE_FFTLIB MATCHES FFTW)
   set(fft fft_fftw.f90)
endif()

if(USE_SHTNS MATCHES yes AND NOT SHTNS MATCHES SHTNS-NOTFOUND)
   list(REMOVE_ITEM srcfiles "${SRC}/shtransforms.f90")
   list(REMOVE_ITEM srcfiles "${SRC}/sht_native.f90")
   set(whichsht shtns.f90)
else()
   list(REMOVE_ITEM srcfiles "${SRC}/shtns.f90")
   set(whichsht sht_native.f90)
endif()

add_executable(magic.exe magic.f90 ${srcfiles} ${fft} ${dct} ${lapack} ${wichsht})
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
