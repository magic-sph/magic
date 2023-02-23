cmake_minimum_required (VERSION 3.0)
project (magic)
enable_language (Fortran C)
cmake_policy(SET CMP0054 NEW)

# generate a Git version string
execute_process(COMMAND git describe --abbrev=7 --always --tags OUTPUT_VARIABLE GIT_VERSION OUTPUT_STRIP_TRAILING_WHITESPACE)

# build date
execute_process(COMMAND date "+%Y/%m/%d %H:%M:%S" OUTPUT_VARIABLE BUILD_DATE OUTPUT_STRIP_TRAILING_WHITESPACE)

set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")
set(CMAKE_Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR}/modules)
# Location of the binary output
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})

set(SRC ${CMAKE_SOURCE_DIR}/src)
set(EXEC "magic.exe")

# Input options
set(PRECISION   "dble" CACHE STRING "Precision of the calculations")
set(ARCH        "64"   CACHE STRING "Architecture of the platform: 32/64")
set(OUT_PREC    "sngl" CACHE STRING "Precision of the binary outputs")
set(USE_FFTLIB  "MKL"  CACHE STRING "Library for the FFTs")
set(USE_DCTLIB  "MKL"  CACHE STRING "Library for the DCTs")
set(USE_MPI     "yes"  CACHE STRING "Use MPI parallelization")
set(USE_OMP     "yes"  CACHE STRING "Use OpenMP")
set(USE_LAPACKLIB "MKL"  CACHE STRING "Use the Math Kernel Library for lapack calls")
set(USE_PRECOND "yes"  CACHE STRING "Use matrix preconditioning")
set(USE_SHTNS   "no"   CACHE STRING "Use libshtns for SH projection")
set(USE_GPU     "no"   CACHE STRING "Use GPU")
find_package(SSE)

# Specific compiler options
if(CMAKE_Fortran_COMPILER_ID MATCHES Intel)
   set(FFLAG_ARCH      " ")
   set(FFLAG_STD       "-stand f08 -diag-disable 7025")
   set(FFLAG_NO_STD    "-nostand")
   set(FFLAG_TRACEBACK "-traceback")
   set(FFLAG_CONVERT   "-convert big")
   set(FFLAG_PREPROCESSOR "-fpp")
   set(FFLAG_WARN_ALL  "-warn all,nointerfaces")
   set(FFLAG_CHECK_ALL "-check all")
   set(FFLAG_FPE       "-fpe0 -fp-stack-check -ftrapuv")
   set(FFLAG_OPT     "-xHost")
endif()
if(CMAKE_Fortran_COMPILER_ID MATCHES GNU)
   set(FFLAG_ARCH      "-m${ARCH}")
   set(FFLAG_STD       "-std=f2008")
   set(FFLAG_NO_STD    "-std=gnu")
   set(FFLAG_TRACEBACK "-fbacktrace")
   set(FFLAG_CONVERT   "-fconvert=big-endian")
   set(FFLAG_PREPROCESSOR "-cpp")
   set(FFLAG_WARN_ALL  "-Wall -Wno-maybe-uninitialized")
   set(FFLAG_CHECK_ALL "-fcheck=all")
   set(FFLAG_FPE       "-ffpe-trap=invalid,zero,overflow")
   set(FFLAG_OPT      "-march=native")
endif()
if(CMAKE_Fortran_COMPILER_ID MATCHES Cray)
   set(FFLAG_ARCH      " ")
   set(FFLAG_STD       "-en")
   set(FFLAG_NO_STD    " ")
   set(FFLAG_TRACEBACK " ")
   set(FFLAG_CONVERT   "-h byteswapio")
   set(FFLAG_PREPROCESSOR "-eZz")
   set(FFLAG_WARN_ALL  " ")
   set(FFLAG_CHECK_ALL "-G 0 -R b")
   set(FFLAG_FPE       "-K trap=fp")
   set(FFLAG_OPT       " ")
   #set(CMAKE_EXE_LINKER_FLAGS "-h dynamic")
   # Remove -rdynamic option (workaround for test machine):
   # https://cmake.org/pipermail/cmake/2014-November/059154.html
   set(CMAKE_SHARED_LIBRARY_LINK_Fortran_FLAGS "")
endif()
if(CMAKE_Fortran_COMPILER_ID MATCHES PGI)
   set(FFLAG_ARCH      "-m${ARCH}")
   set(FFLAG_STD       " ")
   set(FFLAG_NO_STD    " ")
   set(FFLAG_TRACEBACK "-traceback")
   set(FFLAG_CONVERT   "-Mbyteswapio")
   set(FFLAG_PREPROCESSOR "-Mpreprocess")
   set(FFLAG_WARN_ALL  "-Wall -Wno-maybe-uninitialized")
   set(FFLAG_CHECK_ALL "-fcheck=all")
   set(FFLAG_FPE       "-Ktrap=divz,inv,ovf")
   set(FFLAG_OPT      "-fastsse")
endif()

# Matrix preconditioning
if (USE_PRECOND MATCHES yes)
   add_definitions(-DWITH_PRECOND_S -DWITH_PRECOND_Z -DWITH_PRECOND_BJ
                   -DWITH_PRECOND_Z10)
endif()

# Git version
add_definitions(-DGIT_VERSION="${GIT_VERSION}")
# Build date
add_definitions(-DBUILD_DATE="${BUILD_DATE}")
# Precision of the calculations
add_definitions(-Dsngl=1 -Ddble=2)
add_definitions(-DDEFAULT_PRECISION=${PRECISION})
add_definitions(-DDEFAULT_OUTPUT_PRECISION=${OUT_PREC})
message(STATUS "Set architecture to '${ARCH}'")
message(STATUS "Set precision to '${PRECISION}'")
message(STATUS "Set output precision to '${OUT_PREC}'")

#find_package(MKL)

# Start with the list of MPI wrappers. Acknowledgement to http://www.openflipper.org/svnrepo/CoMISo/trunk/CoMISo/cmake/FindMPI.cmake
# Generic MPI compilers
set(_MPI_C_COMPILER_NAMES                  "mpicc|mpcc|mpicc_r|mpcc_r")
set(_MPI_Fortran_COMPILER_NAMES            "mpif90|mpf90|mpif90_r|mpf90_r")
# GNU compiler names
set(_MPI_GNU_C_COMPILER_NAMES              "mpigcc|mpgcc|mpigcc_r|mpgcc_r")
set(_MPI_GNU_Fortran_COMPILER_NAMES        "mpigfortran|mpgfortran|mpigfortran_r|mpgfortran_r")
# Intel MPI compiler names
set(_MPI_Intel_C_COMPILER_NAMES            "mpiicc")
set(_MPI_Intel_Fortran_COMPILER_NAMES      "mpiifort")
# PGI compiler names
set(_MPI_PGI_C_COMPILER_NAMES              "mpipgcc|mppgcc")
set(_MPI_PGI_Fortran_COMPILER_NAMES        "mpipgf|mppgf")
# Cray compiler names
set(_MPI_Cray_C_COMPILER_NAMES             "cc")
set(_MPI_Cray_Fortran_COMPILER_NAMES       "ftn")
# XLC MPI Compiler names
set(_MPI_XL_C_COMPILER_NAMES               "mpxlc|mpxlc_r|mpixlc|mpixlc_r")
set(_MPI_XL_Fortran_COMPILER_NAMES         "mpxlf|mpixlf")

# Find CC, CXX and Fortran MPI wrappers
string(REGEX MATCH "${_MPI_C_COMPILER_NAMES}|${_MPI_GNU_C_COMPILER_NAMES}|${_MPI_Intel_C_COMPILER_NAMES}|${_MPI_PGI_C_COMPILER_NAMES}|${_MPI_XL_C_COMPILER_NAMES}|${_MPI_Cray_C_COMPILER_NAMES}" MPICC_COMPILER "${CMAKE_C_COMPILER}")
string(REGEX MATCH "${_MPI_Fortran_COMPILER_NAMES}|${_MPI_GNU_Fortran_COMPILER_NAMES}|${_MPI_Intel_Fortran_COMPILER_NAMES}|${_MPI_PGI_Fortran_COMPILER_NAMES}|${_MPI_XL_Fortran_COMPILER_NAMES}|${_MPI_Cray_Fortran_COMPILER_NAMES}" MPIFORTRAN_COMPILER "${CMAKE_Fortran_COMPILER}")

# Check whether compiler_options() is available
file(WRITE "${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeTmp/testCompOpt.f90"
"
program TestCompOpt
use iso_fortran_env
print*, compiler_options()
print*, compiler_version()
end program TestCompOpt
")
try_compile(COMP_OPT_AV ${CMAKE_BINARY_DIR}
   ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeTmp/testCompOpt.f90
   COMPILE_DEFINITIONS ${CMAKE_REQUIRED_DEFINITIONS}
   OUTPUT_VARIABLE OUTPUT)
if(COMP_OPT_AV)
   add_definitions(-DCOMP_OPT)
   message(STATUS "Fortran 2008 compiler_options(), compiler_version() available")
else()
   message(STATUS "Fortran 2008 compiler_options(), compiler_version() not available")
endif()

# Look for MPI and add the required flags if found
if(USE_MPI MATCHES yes)
   if("${MPICC_COMPILER}" STREQUAL "" AND "${MPIFORTRAN_COMPILER}" STREQUAL "")
      message(STATUS "Could not find MPI wrappers in CMAKE_C_COMPILER or CMAKE_Fortran_COMPILER. Trying to find MPI libs by default")
      set(MPI_WRAPPER_FOUND FALSE CACHE BOOL "Did not find the MPI Wrappers")

      find_package(MPI)

      if(MPI_Fortran_FOUND)
         message(STATUS "Use MPI")
         add_definitions(-DWITH_MPI)
         add_definitions(${MPI_Fortran_COMPILE_FLAGS})
         link_directories(${MPI_LIBRARY_DIRS})
         include_directories(${MPI_Fortran_INCLUDE_PATH})
      else()
         message(STATUS "MPI was not found: serial job")
      endif()
   else()
      # Setting wrappers and mpi found to true
      set(MPI_WRAPPER_FOUND TRUE CACHE BOOL "found the MPI Wrappers")
      set(MPI_FOUND TRUE CACHE BOOL "Found the MPI library")
      message(STATUS "Use MPI")
         add_definitions(-DWITH_MPI)
   endif()
else()
      message(STATUS "No MPI: serial job")
endif()

# Check whether openMP can be used
if (USE_MPI MATCHES yes)
   if(MPI_WRAPPER_FOUND OR MPI_Fortran_FOUND)
      if(USE_OMP MATCHES yes)
         find_package(OpenMP_Fortran)
         message(STATUS "Use OpenMP")
         add_definitions(-DWITHOMP)
         set(CMAKE_Fortran_FLAGS ${CMAKE_Fortran_FLAGS} ${OpenMP_Fortran_FLAGS})
      else()
         message(STATUS "No OpenMP")
      endif()
   endif()
endif()

# Detect the MKL library and find the libraries
if(USE_LAPACKLIB MATCHES MKL OR USE_FFTLIB MATCHES MKL)
   if (EXISTS $ENV{MKLROOT})
      set (MKLROOT "$ENV{MKLROOT}")
      message(STATUS "MKL was found: ${MKLROOT}")
      find_library(MKLCORE mkl_core HINTS ${MKLROOT}/lib/intel64)
      find_library(MKLLAPACK mkl_lapack95_lp64 HINTS ${MKLROOT}/lib/intel64)
      if (USE_SHTNS MATCHES yes AND USE_OMP MATCHES yes)
         find_library(MKLSEQ mkl_gnu_thread HINTS ${MKLROOT}/lib/intel64)
      else()
         find_library(MKLSEQ mkl_sequential HINTS ${MKLROOT}/lib/intel64)
      endif()
      if(CMAKE_Fortran_COMPILER_ID MATCHES Intel)
         find_library(MKL64 mkl_intel_lp64 HINTS ${MKLROOT}/lib/intel64)
      endif()
      if(CMAKE_Fortran_COMPILER_ID MATCHES GNU)
         find_library(MKL64 mkl_gf_lp64 HINTS ${MKLROOT}/lib/intel64)
      endif()
      link_directories(${MKLROOT}/lib/intel64)
      include_directories(${MKLROOT}/include ${MKLROOT}/include/fftw)
   else()
      message(STATUS "MKL was not found")
      set(USE_LAPACKLIB "LIBFLAME")
      set(USE_FFTLIB "FFTW")
   endif()
endif()


# Detect BLIS and LIBFLAME
if(USE_LAPACKLIB MATCHES LIBFLAME)
   find_library(FLAME NAME flame
       PATHS $ENV{HOME}/local/lib
             $ENV{HOME}/local/lib64)
   find_library(BLIS NAME blis
       PATHS $ENV{HOME}/local/lib
             $ENV{HOME}/local/lib64)
   if ( FLAME MATCHES FLAME-NOTFOUND )
      message(STATUS "LIBFLAME was not found")
      set(USE_LAPACKLIB "LIBSCI")
   else()
      message(STATUS "BLIS: '${BLIS}'")
      message(STATUS "LIBFLAME: '${FLAME}'")
      add_definitions(-DWITH_LIBFLAME)
   endif()
endif()

# Detect LIBSCI (Cray)
if(USE_LAPACKLIB MATCHES LIBSCI)
   if(CMAKE_Fortran_COMPILER_ID MATCHES Cray)
      message(STATUS "LIBSCI: 'built-in with ftn'")
   else()
      find_library(SCI NAME sci_cray_mp
         PATHS $ENV{CRAY_LIBSCI_PREFIX}/lib
               $ENV{HOME}/local/lib
               $ENV{HOME}/local/lib64)
      if ( SCI MATCHES SCI-NOTFOUND )
         message(STATUS "LIBSCI was not found")
         set(USE_LAPACKLIB "LAPACK")
      else()
         message(STATUS "LIBSCI: '${SCI}'")
      endif()
   endif()
endif()

# Detect LAPACK
if(USE_LAPACKLIB MATCHES LAPACK)
   find_package(LAPACK)
   if (LAPACK_FOUND)
      message(STATUS "LAPACK was found")
      set(USE_LAPACKLIB "LAPACK")
   else()
      message(STATUS "LAPACK was not found")
      set(USE_LAPACKLIB "JW")
   endif()
endif()


find_library(FFTW3 NAME fftw3
    PATHS $ENV{HOME}/local/lib
          $ENV{HOME}/.local/lib
          $ENV{HOME}/local/lib64
          /usr/local/fftw/3.3.4_openmp/lib
          $ENV{FFTW_ROOT}/lib
          $ENV{FFTW3_ROOT}/lib
          $ENV{FFTW_HOME}/lib                 #Hydra
          $ENV{FFTWDIR})                      #GWDG


if(USE_OMP)
   find_library(FFTW3_OMP NAME fftw3_omp
       PATHS $ENV{HOME}/local/lib
             $ENV{HOME}/local/lib64
             /usr/local/fftw/3.3.4_openmp/lib
             $ENV{FFTW_ROOT}/lib
             $ENV{FFTW3_ROOT}/lib
             $ENV{FFTW_HOME}/lib              #Hydra
             $ENV{FFTWDIR})                   #GWDG
   message(STATUS "FFTW3_OMP: '${FFTW3_OMP}'")
endif()


if( USE_FFTLIB MATCHES FFTW )
   include_directories(/usr/include)
   if ( FFTW3-NOTFOUND )
      message(STATUS "FFTW3 was not found")
      set(USE_FFTLIB "JW")
   else()
      message(STATUS "FFTW3: '${FFTW3}'")
      set(USE_FFTLIB "FFTW")
   endif()
endif()

if ( USE_FFTLIB MATCHES MKL )
   set(USE_DCTLIB "MKL")
elseif( USE_FFTLIB MATCHES FFTW )
   set(USE_DCTLIB "FFTW")
else()
   set(USE_DCTLIB "JW")
endif()

message(STATUS "Use '${USE_FFTLIB}' for the FFTs")
message(STATUS "Use '${USE_DCTLIB}' for the DCTs")
message(STATUS "Use '${USE_LAPACKLIB}' for the LU factorisations")

if(APPLE)
    include_directories(/usr/local/include)
endif()

include_directories(/usr/include)  #For Ubuntu-like systems

# Detect SHTNS
if(USE_SHTNS MATCHES yes)
   set(ROCM_PATH $ENV{ROCM_PATH})
   if(USE_GPU)
      find_library(SHTNS shtns_cuda HINTS $ENV{HOME}/local/lib)
   else()
      if(USE_OMP)
         find_library(SHTNS shtns_cuda HINTS $ENV{HOME}/local/lib) #Need for batched version
      else()
         find_library(SHTNS shtns HINTS $ENV{HOME}/local/lib)
      endif(USE_OMP)
   endif(USE_GPU)
   if(NOT SHTNS MATCHES SHTNS-NOTFOUND)
      add_definitions(-DWITH_SHTNS)
      include_directories($ENV{HOME}/local/include)
      message(STATUS "SHTNS: ${SHTNS}")
   else()
      message(STATUS "Failed to find SHTNS")
   endif()
else()
   message(STATUS "Use SHTNS: no")
endif()

# Check whether SHT_ALLOW_PADDING is available
if(USE_SHTNS MATCHES yes AND NOT SHTNS MATCHES SHTNS-NOTFOUND)
   file(WRITE "${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeTmp/testPad.f90"
   "
   program TestPad
   use iso_c_binding
   implicit none
   include 'shtns.f03'
   print*, SHT_ALLOW_PADDING
   end program TestPad
   ")
   try_compile(SHT_PAD ${CMAKE_BINARY_DIR}
      ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeTmp/testPad.f90
      CMAKE_FLAGS "-DINCLUDE_DIRECTORIES=$ENV{HOME}/local/include"
      COMPILE_DEFINITIONS ${CMAKE_REQUIRED_DEFINITIONS}
      OUTPUT_VARIABLE OUTPUT)
   if(SHT_PAD)
      add_definitions(-DSHT_PADDING)
      message(STATUS "theta padding allowed with this version of SHTns")
   else()
      message(STATUS "No theta padding allowed with this version of SHTns")
   endif()
endif()

if(USE_GPU MATCHES yes)
   add_definitions(-DWITH_OMP_GPU)
   message(STATUS $ENV{ROCM_PATH})
   #Detect Hipfort (lib name is hipfort-amdgcn)
   find_library(HIPFORT hipfort-amdgcn $ENV{HIPFORT_PATH}/lib)
   if(NOT HIPFORT MATCHES HIPFORT-NOTFOUND)
      #Test if HIPFORT installation
      file(WRITE "${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeTmp/testHipfort.f90"
      "
      program TestHipfort
      use omp_lib
      use iso_c_binding
      use hipfort_hipfft
      use hipfort_check
      print*, \"Hipfort Hipfft test\"
      end program TestHipfort
      ")
      try_compile(HIPFORT_WORK ${CMAKE_BINARY_DIR}
         ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeTmp/testHipfort.f90
         CMAKE_FLAGS "-DINCLUDE_DIRECTORIES=$ENV{HIPFORT_PATH}/include/hipfort/amdgcn"
         COMPILE_DEFINITIONS ${CMAKE_REQUIRED_DEFINITIONS}
         LINK_LIBRARIES "-L${ROCM_PATH}/hipfft/lib/ -lhipfft -L${ROCM_PATH}/rocfft/lib/ -lrocfft -L$ENV{HIPFORT_PATH}/lib -lhipfort-amdgcn"
         OUTPUT_VARIABLE OUTPUT)
      if(HIPFORT_WORK)
         add_definitions(-DWITH_HIPFORT)
         include_directories("$ENV{HIPFORT_PATH}/include/hipfort/amdgcn")
         message(STATUS "HIPFORT: ${HIPFORT}")
      else()
         message(STATUS ${OUTPUT})
         message(STATUS "Failed to test HIPFORT installation")
      endif()
   else()
           message(STATUS "Failed to find HIPFORT")
   endif()
endif()

# Final Fortran flags (after OpenMP/MPI to avoid the conflict
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${FFLAG_ARCH} ${FFLAG_STD} -g ${FFLAG_TRACEBACK} ${FFLAG_CONVERT} ${FFLAG_PREPROCESSOR}")

if(USE_GPU MATCHES yes)
   set(CMAKE_Fortran_FLAGS_DEBUG  "-h keepfiles -hmsgs -hlist=a -O0 -h acc_model=auto_async_none:no_fast_addr:deep_copy ${FFLAG_CHECK_ALL} ${FFLAG_WARN_ALL} ${FFLAG_FPE}")
   set(CMAKE_Fortran_FLAGS_RELEASE  "-h keepfiles -hmsgs -hlist=a -O2 -h acc_model=auto_async_none:no_fast_addr:deep_copy ${FFLAG_OPT}")
else()
   set(CMAKE_Fortran_FLAGS_DEBUG  "-O0 ${FFLAG_CHECK_ALL} ${FFLAG_WARN_ALL} ${FFLAG_FPE}")
   set(CMAKE_Fortran_FLAGS_RELEASE  "-O2 ${FFLAG_OPT}")
endif()

# Default to Release build type
if(NOT CMAKE_BUILD_TYPE)
   set(CMAKE_BUILD_TYPE "Release")
endif()

# Print FLAGS
message(STATUS "Compilation flags: ${CMAKE_Fortran_FLAGS}")
if(CMAKE_BUILD_TYPE MATCHES "Release")
   message(STATUS "Optimisation flags: ${CMAKE_Fortran_FLAGS_RELEASE}")
elseif(CMAKE_BUILD_TYPE MATCHES "Debug")
   message(STATUS "Optimisation flags: ${CMAKE_Fortran_FLAGS_DEBUG}")
endif()

add_subdirectory(${SRC})
