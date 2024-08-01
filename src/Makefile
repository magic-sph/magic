# COMPILER: intel, gnu, portland
COMPILER=intel
# precision of the calculations: sngl/dble
PRECISION=dble
# precision for binary outputs: sngl/dble
OUT_PREC=sngl
# production run: yes/no (no= -O0, yes = -03)
PRODRUN=yes
#USE_MPI: yes/no
USE_MPI=yes
USE_OMP=yes
# DEBUG: no/noopt/check/all/marmot
DEBUG=no
# Pre-conditioning of the matrices: yes/no
USE_PRECOND=yes
# USE_FFTLIB can be MKL, FFTW or JW
USE_FFTLIB=MKL
#USE_DCTLIB can be MKL, FFTW or JW
USE_DCTLIB=MKL
#USE_LAPACKLIB can be MKL, LAPACK or JW
USE_LAPACKLIB=MKL
#Use shtns for Legendre/Fourier transforms
USE_SHTNS=no

USE_PERFLIB=none

ifeq ($(USE_MPI),no)
	USE_OMP=no
endif

SHTNS_INC = $(HOME)/local/include
ifeq ($(USE_OMP),yes)
	FFLAG_SHTNS = -L$(HOME)/local/lib -lshtns_omp -L$(MKLROOT)/lib/intel64 -lmkl_intel_lp64 -lmkl_gnu_thread -lmkl_core -lrt -lm -Wl,-rpath,$(MKLROOT)/lib/intel64
else
	FFLAG_SHTNS = -L$(HOME)/local/lib -lshtns -L$(MKLROOT)/lib/intel64 -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lrt -lm -Wl,-rpath,$(MKLROOT)/lib/intel64
endif


# In case MKL is not found, switch back to built-in libraries
ifndef MKLROOT
	USE_FFTLIB=FFTW
endif
ifeq ($(USE_LAPACKLIB),MKL)
ifndef MKLROOT
	USE_LAPACKLIB=JW
endif
endif

ifeq ($(USE_PRECOND),yes)
	PREPROC +=	-DWITH_PRECOND_S \
			-DWITH_PRECOND_Z \
			-DWITH_PRECOND_BJ \
			-DWITH_PRECOND_Z10
endif

ifeq ($(COMPILER),intel)
	FFLAG_TRACEBACK = -traceback
	FFLAG_ARCH_OPT = -xHost
	FFLAG_OPENMP = -qopenmp
	FFLAG_OPENMP_STUBS = -openmp-stubs
	FFLAG_CHECK_ALL = -check all
	FFLAG_FPE = -fpe0 -fp-stack-check -ftrapuv
	FFLAG_WARN_ALL  = -warn all,nointerfaces
	FFLAG_PREPROCESSOR = -fpp
	FFLAG_CONVERT = -convert big
	FFLAG_PROFILING = -prof-gen=srcpos
	FFLAG_REPORT = -vec-report6
	FFLAG_STD = -stand f08  -diag-disable 7025
	ifneq ($(PRODRUN),yes)
		FFLAG_STD += -fp-model source
	endif
	COMP_FC      = ifort
	COMP_CC      = icc
	ifeq ($(DEBUG),marmot)
		COMP_MPFC = marmotf90
	else
		COMP_MPFC = mpiifort
		COMP_MPCC = mpiicc
	endif	

	FFLAG_MKL_LIBS = -L$(MKLROOT)/lib/intel64 -lmkl_intel_lp64 \
                         -lmkl_lapack95_lp64\
			 -lmkl_sequential -lmkl_core -lpthread -lm \
			 -Wl,-rpath,$(MKLROOT)/lib/intel64 

endif
ifeq ($(COMPILER),gnu)
	FFLAGS = -m64
	FFLAG_ARCH_OPT = -march=native
	FFLAG_TRACEBACK = -fbacktrace
	FFLAG_OPENMP = -fopenmp
	FFLAG_CONVERT = -fconvert=big-endian
	FFLAG_STD = -std=f2008
	FFLAG_PREPROCESSOR = -cpp
	FFLAG_WARN_ALL = -Wall -Wno-maybe-uninitialized
	FFLAG_CHECK_ALL = -fcheck=all
	FFLAG_FPE = -ffpe-trap=invalid,zero,overflow 
	COMP_FC = gfortran
	COMP_CC = gcc
	COMP_MPFC = mpif90
	COMP_MPCC = mpicc

	FFLAG_MKL_LIBS = -L$(MKLROOT)/lib/intel64 -lmkl_gf_lp64\
                         -lmkl_lapack95_lp64\
			 -lmkl_sequential -lmkl_core -lpthread -lm \
			 -Wl,-rpath,$(MKLROOT)/lib/intel64
endif
ifeq ($(COMPILER),portland)
	FFLAGS = -m64
	FFLAG_ARCH_OPT = -fastsse
	FFLAG_TRACEBACK = -traceback
	FFLAG_OPENMP = -mp
	FFLAG_CONVERT = -Mbyteswapio
	FFLAG_PREPROCESSOR = -Mpreprocess
	FFLAG_WARN_ALL = -Wall -Wno-maybe-uninitialized
	FFLAG_CHECK_ALL = -fcheck=all
	FFLAG_FPE = -Ktrap=divz,inv,ovf
	COMP_FC = pgf95
	COMP_CC = pgcc
	COMP_MPFC = mpif90
	COMP_MPCC = mpicc
	FFLAG_MKL_LIBS = -L$(MKLROOT)/lib/intel64 -lmkl_gf_lp64\
                         -lmkl_lapack95_lp64\
			 -lmkl_sequential -lmkl_core -lpthread -lm \
			 -Wl,-rpath,$(MKLROOT)/lib/intel64
endif


# the following path is only used, if the mpi.mod file is
# not found in the standard location. It is used for example
# on hydra, if one compiles with the gnu compilers.
#MPI_INCPATH=/usr/lib64/mpi/gcc/openmpi/include
#MPI_INCPATH=/opt/ibmhpc/pecurrent/mpich2/gnu/include64

# generate a Git version string
GIT_VERSION := $(shell git describe --abbrev=7 --always --tags)
# get the date of the build
BUILD_DATE := $(shell date "+%Y/%m/%d %H:%M:%S")

PREPROC += -Dsngl=1 -Ddble=2 -DGIT_VERSION=\"$(GIT_VERSION)\" #-DBUILD_DATE=\"$(BUILD_DATE)\"
ifeq ($(USE_MPI),yes)
	PREPROC += -DWITH_MPI
	FC = $(COMP_MPFC)
	CC = $(COMP_MPCC)
else
	FC = $(COMP_FC)
	CC = $(COMP_CC)
endif

ifeq ($(USE_PERFLIB),perf)
     PREPROC += -DWITHPERF
#     INCPATHS += -I$(PERFLIB_HOME)/include
#     INCPATHS += -I/afs/.ipp/common/soft/perflib/2.2/intel-12.1.6/include
     LIBS += -L$(PERFLIB_HOME)/lib -looperf -Wl,-rpath,$(PERFLIB_HOME)/lib
endif
ifeq ($(USE_PERFLIB),scalasca)
	FC:= scalasca -instrument $(FC)
endif
ifeq ($(USE_PERFLIB),itac)
	FFLAGS += -tcollect -tcollect-filter itac_filter.list
endif
ifeq ($(USE_PERFLIB),coverage)
	FFLAGS += $(FFLAG_PROFILING)
endif
ifeq ($(USE_PERFLIB),mpitrace)
	LIBS += -L$(IHPCT_BASE)/lib64 -lmpitrace
endif

ifeq ($(USE_PERFLIB),likwid)
	PREPROC += -DWITH_LIKWID
	INCPATHS += -I$(LIKWID_HOME)/include
	LIBS += -L$(LIKWID_HOME)/lib -llikwid
endif

ifeq ($(USE_FFTLIB),MKL)
	INCPATHS += -I$(MKLROOT)/include
	INCPATHS += -I$(MKLROOT)/include/fftw
	LIBS +=  $(FFLAG_MKL_LIBS)
endif

ifeq ($(USE_LAPACKLIB),MKL)
	INCPATHS += -I$(MKLROOT)/include
	LIBS += $(FFLAG_MKL_LIBS)
endif

ifeq ($(USE_LAPACKLIB),LAPACK)
	LIBS += -llapack
endif

ifeq ($(USE_SHTNS),yes)
	PREPROC += -DWITH_SHTNS
	INCPATHS += -I$(SHTNS_INC)
	LIBS += $(FFLAG_SHTNS)
endif

ifeq ($(PRECISION),sngl)
	PREPROC += -DDEFAULT_PRECISION=sngl
else ifeq ($(PRECISION),dble)
	PREPROC += -DDEFAULT_PRECISION=dble
endif

ifeq ($(OUT_PREC),sngl)
	PREPROC += -DDEFAULT_OUTPUT_PRECISION=sngl
else ifeq ($(OUT_PREC),dble)
	PREPROC += -DDEFAULT_OUTPUT_PRECISION=dble
endif

#FFLAGS += -heap-arrays 1024
FFLAGS += $(FFLAG_STD)
FFLAGS += -g $(FFLAG_TRACEBACK)
ifeq ($(PRODRUN),yes)
	DEBUG=no
	OPT += -O3 $(FFLAG_ARCH_OPT)
else
	OPT += -O 
endif
ifeq ($(USE_OMP),yes)
	FFLAGS += $(FFLAG_OPENMP)
	PREPROC += -DWITHOMP
endif
ifeq ($(USE_OMP),stubs)
	FFLAGS += $(FFLAG_OPENMP_STUBS)
endif

FFLAGS_NOCHECK=$(FFLAGS)
ifeq ($(DEBUG),noopt)
	OPT = -O0
else ifeq ($(DEBUG),check)
	OPT = -O0
	PREPROC += -DWITH_DEBUG
	FFLAGS += $(FFLAG_CHECK_ALL) $(FFLAG_WARN_ALL)
else ifeq ($(DEBUG),all)
	OPT = -O0
	FFLAGS += $(FFLAG_CHECK_ALL) $(FFLAG_WARN_ALL) $(FFLAG_FPE)
else ifeq ($(DEBUG),marmot)
	LIBS += -L$(MARMOT_HOME)/lib -Wl,-rpath,$(MARMOT_HOME)/lib
endif

INCPATHS += -I.

FFLAGS += $(FFLAG_CONVERT) 
FFLAGS  += $(FFLAG_PREPROCESSOR) $(OPT)
LDFLAGS =

FFLAGS_NOCHECK = $(filter-out $(FFLAG_CHECK_ALL),$(FFLAGS))

OUT     = magic.exe
SOURCES:= $(wildcard *.f90)
RED_SOURCES = $(filter-out fft% algebra% sht% dct_fftw% cosine_transform_odd% mpimod%,$(SOURCES))
RED_SOURCES += fft_fac.o

ifeq ($(USE_FFTLIB),JW)
	RED_SOURCES += fft.f90
	FFT_OBJS += fft.o
else ifeq ($(USE_FFTLIB),FFTW)
	RED_SOURCES += fft_fftw.f90
	FFT_OBJS += fft_fftw.o
else ifeq ($(USE_FFTLIB),MKL)
	RED_SOURCES += fft_fftw.f90
	FFT_OBJS += fft_fftw.o
endif

ifeq ($(USE_LAPACKLIB),MKL)
	RED_SOURCES += algebra_lapack.f90
	LAPACK_OBJS += algebra_lapack.o
else ifeq ($(USE_LAPACKLIB),LAPACK)
	RED_SOURCES += algebra_lapack.f90
	LAPACK_OBJS += algebra_lapack.o
else ifeq ($(USE_LAPACKLIB),LIBFLAME)
	RED_SOURCES += algebra_lapack.f90
	LAPACK_OBJS += algebra_lapack.o
else ifeq ($(USE_LAPACKLIB),JW)
	RED_SOURCES += algebra.f90
	LAPACK_OBJS += algebra.o
endif

ifeq ($(USE_DCTLIB),MKL)
	RED_SOURCES += dct_fftw.f90
	DCT_OBJS += dct_fftw.o
else ifeq ($(USE_DCTLIB),FFTW)
	RED_SOURCES += dct_fftw.f90
	DCT_OBJS += dct_fftw.o
else ifeq ($(USE_DCTLIB),JW)
	RED_SOURCES += cosine_transform_odd.f90
	DCT_OBJS += cosine_transform_odd.o
endif

ifeq ($(USE_SHTNS),yes)
	RED_SOURCES += shtns.f90
	SHT_OBJS += shtns.o
else
	RED_SOURCES += sht_native.f90 shtransforms.f90
	SHT_OBJS += sht_native.o
endif

OBJS   := $(addsuffix .o, $(basename $(RED_SOURCES)))
OBJS += truncation.o

.SUFFIXES: 

#ifeq ($(COMPILER),gnu)
#ifeq ($(USE_MPI),yes)
#all: mpi.mod
#	make $(OUT)
#endif
#endif

$(OUT): $(OBJS)
	$(FC) $(FFLAGS) -o $@ $^ $(LIBS)

%.o: %.f90
	$(FC) -c $(FFLAGS) $(PREPROC) $(INCPATHS) -o $@ $<

startFields.o: startFields.f90
	$(FC) -c $(FFLAGS_NOCHECK) $(PREPROC) $(INCPATHS) -o $@ $<

timing.o: timing.f90
	$(FC) -c $(FFLAGS_NOCHECK) $(PREPROC) $(INCPATHS) -o $@ $<

%.o: %.c
	$(CC) -c $(CFLAGS) $(PREPROC) $(INCPATHS) -o $@ $<

# cancelling the implicit rule for .F files
#%.f: %.F

precision_mod.o       : mpimod.o
truncation.o          : precision_mod.o logic.o useful.o
signals.o             : precision_mod.o char_manip.o output_data.o\
                        parallel.o
magnetic_energy.o     : truncation.o blocking.o horizontal.o\
                        phys_param.o num_param.o constants.o radial.o\
                        logic.o output_data.o movie.o special.o\
                        integration.o useful.o  mem_alloc.o
finite_differences.o  : truncation.o useful.o constants.o parallel.o\
                        radial_scheme.o mem_alloc.o
chebyshev.o           : mem_alloc.o constants.o blocking.o radial_scheme.o\
                        chebyshev_polynoms.o $(DCT_OBJS) useful.o
radial_scheme.o       : precision_mod.o
useful.o              : logic.o output_data.o parallel.o precision_mod.o\
                        constants.o
LMmapping.o           : mem_alloc.o truncation.o
mem_alloc.o           : parallel.o precision_mod.o output_data.o
out_TO.o              : truncation.o horizontal.o num_param.o logic.o\
                        phys_param.o radial.o parallel.o integration.o\
                        constants.o output_data.o TO.o useful.o mem_alloc.o
dtB_arrays.o          : truncation.o precision_mod.o constants.o\
                        mem_alloc.o
courant.o             : truncation.o blocking.o num_param.o\
                        phys_param.o radial.o horizontal.o logic.o\
                        useful.o radial_data.o parallel.o
nonlinear_bcs.o       : truncation.o horizontal.o $(FFT_OBJS)\
                        $(SHT_OBJS) blocking.o phys_param.o\
                        radial.o radial_data.o
fields_average.o      : truncation.o output_data.o blocking.o time_schemes.o\
                        logic.o constants.o radial.o horizontal.o\
                        radial_derivatives.o radial_derivatives_even.o\
                        spectra.o kinetic_energy.o out_graph_file.o\
                        magnetic_energy.o $(SHT_OBJS)\
                        $(FFT_OBJS) out_coeff.o storeCheckPoints.o\
                        mem_alloc.o
get_nl.o              : truncation.o horizontal.o logic.o time_schemes.o\
                        phys_param.o radial.o blocking.o mem_alloc.o
out_movie_file.o      : truncation.o output_data.o blocking.o\
                        radial.o horizontal.o movie.o fields.o\
                        num_param.o $(FFT_OBJS) logic.o outMisc.o\
                        out_dtB_frame.o parallel.o communications.o
storeCheckPoints.o    : truncation.o phys_param.o num_param.o logic.o\
                        output_data.o init_fields.o dt_fieldsLast.o\
                        communications.o parallel.o char_manip.o time_array.o\
                        time_schemes.o
updateS.o             : truncation.o blocking.o horizontal.o time_array.o\
                        phys_param.o radial.o logic.o output_data.o\
                        init_fields.o constants.o mem_alloc.o time_schemes.o\
                        matrices.o $(DCT_OBJS) radial_derivatives.o matrices.o
updateXI.o            : truncation.o blocking.o horizontal.o time_array.o\
                        phys_param.o radial.o logic.o output_data.o\
                        init_fields.o constants.o mem_alloc.o time_schemes.o\
                        matrices.o $(DCT_OBJS) radial_derivatives.o
power.o               : truncation.o blocking.o horizontal.o\
                        phys_param.o num_param.o radial.o logic.o\
                        output_data.o outRot.o integration.o useful.o\
                        mem_alloc.o mean_sd.o
integration.o         : $(DCT_OBJS) precision_mod.o constants.o chebyshev.o\
                        communications.o
out_coeff.o           : logic.o precision_mod.o parallel.o blocking.o\
                        truncation.o communications.o\
                        mem_alloc.o radial.o output_data.o phys_param.o
magic.o               : truncation.o num_param.o parallel.o logic.o\
                        output_data.o timing.o constants.o precision_mod.o\
                        useful.o startFields.o preCalculations.o time_schemes.o\
                        step_time.o kinetic_energy.o mem_alloc.o\
			magnetic_energy.o horizontal.o fields_average.o \
			outGeos.o spectra.o out_TO.o output_data.o outMisc.o probes.o\
                        Namelists.o blocking.o out_coeff.o
step_time.o           : truncation.o phys_param.o updateB.o updateZ.o\
                        num_param.o radial_data.o updateWP.o time_schemes.o\
                        logic.o output_data.o output.o movie.o updateS.o\
                        timing.o courant.o signals.o updateXI.o updatePHI.o\
                        radialLoop.o nonlinear_bcs.o updateWPS.o\
                        LMLoop.o dt_fieldsLast.o useful.o probes.o
phys_param.o          : precision_mod.o
char_manip.o          : precision_mod.o
constants.o           : precision_mod.o
outMisc.o             : truncation.o num_param.o phys_param.o logic.o\
                        output_data.o blocking.o constants.o radial.o\
                        horizontal.o integration.o startFields.o
chebyshev_polynoms.o  : constants.o logic.o precision_mod.o num_param.o
rIteration.o          : precision_mod.o time_schemes.o truncation.o radial_data.o
rIter.o               : precision_mod.o parallel.o truncation.o logic.o radial_data.o\
                        radial.o constants.o get_nl.o get_td.o TO.o dtB.o\
                        out_graph_file.o nonlinear_bcs.o out_movie_file.o\
                        courant.o outRot.o fields.o time_schemes.o phys_param.o\
                        ${SHT_OBJS} outGeos.o
radialLoop.o          : truncation.o time_schemes.o precision_mod.o mem_alloc.o\
                        radial_data.o rIteration.o rIter.o
LMLoop.o              : truncation.o blocking.o parallel.o time_array.o\
                        logic.o radial_data.o updateWPS.o\
                        output_data.o fields.o dt_fieldsLast.o\
                        updateS.o updateB.o mem_alloc.o updateXI.o\
                        updateZ.o updateWP.o  debugging.o
mpi_transpose.o       : precision_mod.o truncation.o radial_data.o\
                        mem_alloc.o parallel.o blocking.o
debugging.o           : precision_mod.o precision_mod.o
timing.o              : parallel.o precision_mod.o
radial_derivatives.o  : $(DCT_OBJS) constants.o precision_mod.o chebyshev.o
output.o              : truncation.o blocking.o phys_param.o\
                        num_param.o logic.o output_data.o radial.o\
                        horizontal.o constants.o fields.o radial_spectra.o\
                        out_movie_file.o outGeos.o time_schemes.o\
                        outRot.o kinetic_energy.o magnetic_energy.o\
                        spectra.o fields_average.o power.o getDlm.o\
                        useful.o out_TO.o outMisc.o out_graph_file.o\
                        storeCheckPoints.o integration.o parallel.o \
                        communications.o out_coeff.o outPar.o dtB.o\
                        store_movie_IC.o RMS.o
communications.o      : blocking.o parallel.o truncation.o\
			radial_data.o mem_alloc.o mpi_transpose.o
outPar.o              : truncation.o blocking.o phys_param.o\
                        logic.o radial.o output_data.o radial_data.o\
                        horizontal.o num_param.o integration.o constants.o mem_alloc.o
spectra.o             : truncation.o blocking.o num_param.o mean_sd.o\
                        horizontal.o phys_param.o radial.o logic.o\
                        output_data.o useful.o integration.o mem_alloc.o
RMS.o                 : truncation.o parallel.o blocking.o num_param.o\
                        phys_param.o radial.o constants.o logic.o mean_sd.o\
                        output_data.o radial_derivatives.o useful.o\
                        RMS_helpers.o integration.o communications.o\
                        horizontal.o dtB.o chebyshev_polynoms.o mem_alloc.o\
                        finite_differences.o radial_scheme.o chebyshev.o
updateWP.o            : truncation.o blocking.o num_param.o integration.o\
                        phys_param.o radial.o horizontal.o logic.o\
                        RMS.o $(DCT_OBJS) radial_derivatives.o\
                        RMS_helpers.o mem_alloc.o fields.o matrices.o\
                        time_array.o time_schemes.o
updateWPS.o           : truncation.o blocking.o num_param.o time_array.o\
                        phys_param.o radial.o horizontal.o logic.o\
                        $(DCT_OBJS) radial_derivatives.o\
                        mem_alloc.o init_fields.o RMS.o time_schemes.o
updatePHI.o           : truncation.o blocking.o num_param.o time_array.o\
                        phys_param.o radial.o horizontal.o logic.o\
                        $(DCT_OBJS) radial_derivatives.o init_fields.o\
                        mem_alloc.o time_schemes.o
get_td.o              : truncation.o blocking.o horizontal.o\
                        phys_param.o num_param.o radial.o logic.o\
                        RMS.o RMS_helpers.o fields.o mem_alloc.o
store_movie_IC.o      : truncation.o blocking.o logic.o radial_data.o\
                        movie.o radial.o horizontal.o\
                        $(SHT_OBJS) $(FFT_OBJS) out_movie_file.o phys_param.o
radial_spectra.o      : truncation.o blocking.o horizontal.o radial.o\
                        num_param.o output_data.o logic.o useful.o\
                        char_manip.o radial_data.o LMmapping.o constants.o
fft_fac.o             : constants.o precision_mod.o mem_alloc.o
cosine_transform_even.o: truncation.o fft_fac.o constants.o useful.o\
                        mem_alloc.o
readCheckPoints.o     : truncation.o blocking.o phys_param.o time_array.o\
                        init_fields.o constants.o radial.o time_schemes.o\
                        logic.o $(DCT_OBJS) finite_differences.o\
                        chebyshev.o radial_scheme.o num_param.o updateB.o\
			updateWP.o updateS.o updateZ.o updateXI.o updateWPS.o
radial_derivatives_even.o : $(DCT_OBJS) precision_mod.o constants.o\
                        cosine_transform_even.o
updateZ.o             : truncation.o blocking.o phys_param.o time_array.o\
                        num_param.o constants.o radial.o horizontal.o\
                        logic.o init_fields.o RMS.o TO.o time_schemes.o\
                        matrices.o $(DCT_OBJS)\
                        radial_derivatives.o outRot.o RMS_helpers.o \
                        communications.o mem_alloc.o special.o
out_graph_file.o      : truncation.o horizontal.o phys_param.o parallel.o\
                        num_param.o radial.o logic.o blocking.o\
                        output_data.o $(FFT_OBJS) $(SHT_OBJS)
startFields.o         : truncation.o phys_param.o num_param.o constants.o\
                        radial.o init_fields.o blocking.o logic.o output_data.o\
                        dt_fieldsLast.o special.o readCheckPoints.o updateWPS.o\
                        fields.o updateS.o updateZ.o updateWP.o updateB.o timing.o\
                        radial_derivatives.o useful.o communications.o updateXI.o\
                        time_schemes.o
getDlm.o              : truncation.o blocking.o num_param.o constants.o\
                        horizontal.o radial.o constants.o\
                        useful.o integration.o
out_dtB_frame.o       : truncation.o blocking.o radial.o logic.o\
                        horizontal.o dtB.o radial_derivatives.o\
                        radial_derivatives_even.o $(FFT_OBJS) constants.o
updateB.o             : truncation.o blocking.o horizontal.o fields.o\
                        phys_param.o num_param.o init_fields.o time_array.o\
                        logic.o radial.o special.o constants.o time_schemes.o\
                        RMS.o matrices.o radial_derivatives.o\
                        $(DCT_OBJS) radial_derivatives_even.o\
                        RMS_helpers.o mem_alloc.o
outGeos.o             : truncation.o horizontal.o radial.o constants.o\
                        output_data.o logic.o radial_data.o parallel.o\
                        communications.o integration.o mem_alloc.o logic.o\
                        movie.o
preCalculations.o     : truncation.o phys_param.o num_param.o constants.o\
                        radial.o horizontal.o init_fields.o time_schemes.o\
                        blocking.o logic.o output_data.o\
                        useful.o integration.o parallel.o special.o
RMS_helpers.o         : truncation.o blocking.o radial.o constants.o\
                        horizontal.o  useful.o integration.o\
                        LMmapping.o 
kinetic_energy.o      : truncation.o blocking.o num_param.o radial.o\
                        horizontal.o output_data.o logic.o constants.o\
                        phys_param.o useful.o integration.o\
                        communications.o mem_alloc.o
outRot.o              : truncation.o blocking.o num_param.o\
                        phys_param.o logic.o radial.o output_data.o\
                        constants.o integration.o special.o horizontal.o
blocking.o            : truncation.o LMmapping.o logic.o parallel.o\
                        output_data.o useful.o precision_mod.o constants.o\
                        mem_alloc.o
num_param.o           : truncation.o precision_mod.o mem_alloc.o timing.o
radial.o              : truncation.o radial_data.o $(LAPACK_OBJS)\
                        phys_param.o logic.o precision_mod.o finite_differences.o\
                        output_data.o $(DCT_OBJS)\
                        chebyshev_polynoms.o radial_derivatives.o\
			cosine_transform_even.o mem_alloc.o useful.o radial_scheme.o\
                        chebyshev.o finite_differences.o
output_data.o         : precision_mod.o mem_alloc.o
horizontal.o          : truncation.o phys_param.o num_param.o precision_mod.o\
                        radial.o logic.o blocking.o plms.o $(FFT_OBJS)\
                        mem_alloc.o
RMS.o	              : truncation.o blocking.o mem_alloc.o
init_fields.o         : truncation.o blocking.o radial.o horizontal.o\
                        constants.o logic.o $(FFT_OBJS) $(SHT_OBJS)\
                        useful.o phys_param.o mpi_transpose.o\
                        $(DCT_OBJS) mem_alloc.o parallel.o 
mean_sd.o             : precision_mod.o mem_alloc.o matrices.o communications.o
movie.o               : truncation.o parallel.o radial_data.o\
                        output_data.o logic.o radial.o mem_alloc.o\
                        horizontal.o char_manip.o useful.o
dtB.o                 : truncation.o radial_data.o parallel.o \
                        $(FFT_OBJS) radial.o logic.o\
                        blocking.o horizontal.o phys_param.o\
                        mpi_transpose.o $(SHT_OBJS)\
                        radial_spectra.o radial_derivatives.o mem_alloc.o
TO.o                  : truncation.o blocking.o horizontal.o\
                        radial_data.o phys_param.o radial.o logic.o\
                        $(SHT_OBJS) mem_alloc.o
$(LAPACK_OBJS)        : precision_mod.o constants.o useful.o
time_schemes.o        : time_array.o logic.o output_data.o precision_mod.o
time_array.o          : precision_mod.o constants.o mem_alloc.o
dirk_schemes.o        : parallel.o precision_mod.o num_param.o constants.o\
                        mem_alloc.o useful.o logic.o output_data.o time_schemes.o\
                        time_array.o
multistep_schemes.o   : parallel.o precision_mod.o num_param.o constants.o\
                        mem_alloc.o useful.o logic.o output_data.o time_schemes.o\
                        time_array.o
dt_fieldsLast.o       : truncation.o blocking.o precision_mod.o time_array.o\
                        logic.o mem_alloc.o constants.o
fields.o              : truncation.o blocking.o radial_data.o communications.o\
                        precision_mod.o	mem_alloc.o logic.o
special.o             : precision_mod.o mem_alloc.o truncation.o
probes.o              : parallel.o precision_mod.o truncation.o radial.o num_param.o\
                        blocking.o horizontal.o output_data.o constants.o
parallel_solver.o     : precision_mod.o parallel.o radial_data.o mem_alloc.o\
                        constants.o blocking.o truncation.o
radial_data.o         : truncation.o parallel.o logic.o
matrices.o            : precision_mod.o mem_alloc.o $(LAPACK_OBJS)
plms.o                : constants.o precision_mod.o useful.o
$(FFT_OBJS)           : truncation.o blocking.o useful.o constants.o
Namelists.o           : constants.o truncation.o phys_param.o num_param.o logic.o \
                        init_fields.o output_data.o parallel.o movie.o\
                        TO.o special.o radial.o char_manip.o probes.o\
                        time_schemes.o dirk_schemes.o multistep_schemes.o
$(DCT_OBJS)           : precision_mod.o fft_fac.o constants.o useful.o mem_alloc.o

ifeq ($(USE_SHTNS),yes)
shtns.o               : truncation.o blocking.o horizontal.o radial.o  parallel.o\
                        precision_mod.o constants.o blocking.o
rIter.o               : shtns.o
magic.o               : shtns.o
else
shtransforms.o        : precision_mod.o mem_alloc.o truncation.o blocking.o horizontal.o\
                        plms.o constants.o ${FFT_OBJS} parallel.o
sht_native.o          : precision_mod.o blocking.o constants.o truncation.o\
                        horizontal.o parallel.o shtransforms.o
rIter.o               : sht_native.o
magic.o               : sht_native.o
endif

ifeq ($(USE_LAPACKLIB),MKL)
$(LAPACK_OBJS)        : lapack95.mod
lapack95.mod          : $(MKLROOT)/include/lapack.f90
	$(FC) -c $(FFLAGS) $(INCPATHS) $<
endif

ifeq ($(USE_DCTLIB),MKL)
$(DCT_OBJS)           :blocking.o
else ifeq ($(USE_DCTLIB),FFTW)
$(DCT_OBJS)           :blocking.o
endif


.PHONY: clean echo
clean:
	@rm -f *.o $(OUT) *.mod
	-rm -f *__genmod.*

echo:
	@echo $(RED_SOURCES)
	@echo $(OBJS)
