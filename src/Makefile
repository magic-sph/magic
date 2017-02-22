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
# USE_FFTLIB can be MKL or JW
USE_FFTLIB=MKL
#USE_LAPACKLIB can be MKL, LAPACK or JW
USE_LAPACKLIB=MKL
#Use hdf5 for input/output of the restart files: yes/no
USE_HDF5=no
#Use shtns for Legendre/Fourier transforms
USE_SHTNS=no

USE_PERFLIB=none

ifdef HDF5_HOME
	HDF5LIB=$(HDF5_HOME)/lib 
	HDF5ICL=$(HDF5_HOME)/include
endif
ifdef HDF5DIR
	HDF5LIB=$(HDF5DIR) 
	HDF5ICL=$(HDF5INCLUDE)
endif

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
	USE_FFTLIB=JW
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
			-DWITH_PRECOND_S0 \
			-DWITH_PRECOND_Z10
endif

ifeq ($(COMPILER),intel)
	FFLAG_TRACEBACK = -traceback
	FFLAG_ARCH_OPT = -xHost
	FFLAG_OPENMP = -openmp
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

	FFLAG_HDF5_LIBS = -L$(HDF5LIB) -lhdf5hl_fortran -lhdf5 \
                          -Wl,-rpath,$(HDF5LIB)

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

PREPROC += -Dsngl=1 -Ddble=2
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

ifeq ($(USE_HDF5), yes)
	PREPROC += -DWITH_HDF5
	INCPATHS += -I$(HDF5ICL)
	LIBS += $(FFLAG_HDF5_LIBS)
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
RED_SOURCES = $(filter-out fft% algebra% shtns% rIterThetaBlocking_shtns% hdf5Helpers% mpimod%,$(SOURCES))
RED_SOURCES += fft_fac.o

ifeq ($(USE_FFTLIB),JW)
	RED_SOURCES += fft.f90
	FFT_OBJS += fft.o
else ifeq ($(USE_FFTLIB),MKL)
	RED_SOURCES += fft_mkl.f90
	FFT_OBJS += fft_mkl.o
endif

ifeq ($(USE_LAPACKLIB),MKL)
	RED_SOURCES += algebra_mkl.f90
	LAPACK_OBJS += algebra_mkl.o
else ifeq ($(USE_LAPACKLIB),LAPACK)
	RED_SOURCES += algebra_lapack.f90
	LAPACK_OBJS += algebra_lapack.o
else ifeq ($(USE_LAPACKLIB),JW)
	RED_SOURCES += algebra.f90
	LAPACK_OBJS += algebra.o
endif

ifeq ($(USE_HDF5),yes)
	RED_SOURCES += hdf5Helpers.f90
endif

ifeq ($(USE_SHTNS),yes)
	RED_SOURCES += shtns.f90 rIterThetaBlocking_shtns.f90
endif

OBJS   := $(addsuffix .o, $(basename $(RED_SOURCES)))
OBJS += truncation.o
OBJS += c_utils.o

.SUFFIXES: 

ifeq ($(COMPILER),gnu)
ifeq ($(USE_MPI),yes)
all: mpi.mod
	make $(OUT)
endif
endif

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

truncation.o          : precision_mod.o logic.o
magnetic_energy.o     : truncation.o blocking.o horizontal.o\
                        phys_param.o num_param.o constants.o radial.o\
                        logic.o output_data.o movie.o special.o\
                        integration.o useful.o LMLoop_data.o mem_alloc.o
finite_differences.o  : truncation.o useful.o constants.o parallel.o\
                        radial_scheme.o mem_alloc.o
chebyshev.o           : mem_alloc.o constants.o LMLoop_data.o radial_scheme.o\
                        chebyshev_polynoms.o cosine_transform_odd.o useful.o
radial_scheme.o       : precision_mod.o
useful.o              : logic.o output_data.o parallel.o precision_mod.o\
                        constants.o
LMmapping.o           : mem_alloc.o truncation.o
mem_alloc.o           : parallel.o precision_mod.o output_data.o
out_TO.o              : truncation.o horizontal.o num_param.o\
                        phys_param.o radial.o logic.o blocking.o\
                        constants.o output_data.o TO.o legendre_grid_to_spec.o\
                        cosine_transform_odd.o chebInt.o plms.o TO_helpers.o\
                        char_manip.o integration.o useful.o mem_alloc.o\
                        communications.o
TO_helpers.o          : truncation.o blocking.o horizontal.o precision_mod.o\
                        constants.o
TO_arrays.o           : truncation.o precision_mod.o mem_alloc.o
dtB_arrays.o          : truncation.o precision_mod.o constants.o\
                        mem_alloc.o
courant.o             : truncation.o blocking.o num_param.o\
                        phys_param.o radial.o horizontal.o logic.o\
                        useful.o radial_data.o parallel.o
nonlinear_bcs.o       : truncation.o horizontal.o $(FFT_OBJS)\
                        legendre_grid_to_spec.o blocking.o phys_param.o\
                        radial.o radial_data.o
fields_average.o      : truncation.o output_data.o blocking.o\
                        logic.o constants.o radial.o horizontal.o\
                        radial_derivatives.o radial_derivatives_even.o\
                        spectra.o kinetic_energy.o out_graph_file.o\
                        magnetic_energy.o legendre_spec_to_grid.o\
                        $(FFT_OBJS) out_coeff.o\
                        legendre_helpers.o mem_alloc.o
get_nl.o              : truncation.o horizontal.o logic.o\
                        phys_param.o radial.o blocking.o mem_alloc.o
out_movie_file.o      : truncation.o output_data.o blocking.o\
                        radial.o horizontal.o movie.o fields.o\
                        num_param.o $(FFT_OBJS) logic.o\
                        out_dtB_frame.o
storeCheckPoints.o    : truncation.o phys_param.o num_param.o logic.o\
                        output_data.o init_fields.o dt_fieldsLast.o 
ifeq ($(USE_HDF5),yes)
readCheckPoints.o     : hdf5Helpers.o
storeCheckPoints.o    : hdf5Helpers.o LMLoop_data.o
hdf5Helpers.o         : blocking.o LMLoop_data.o precision_mod.o
endif
updateS.o             : truncation.o blocking.o horizontal.o\
                        phys_param.o radial.o logic.o output_data.o\
                        init_fields.o constants.o mem_alloc.o\
                        $(LAPACK_OBJS) cosine_transform_odd.o\
                        radial_derivatives.o LMLoop_data.o
updateXI.o            : truncation.o blocking.o horizontal.o\
                        phys_param.o radial.o logic.o output_data.o\
                        init_fields.o constants.o mem_alloc.o\
                        $(LAPACK_OBJS) cosine_transform_odd.o\
                        radial_derivatives.o LMLoop_data.o
power.o               : truncation.o blocking.o horizontal.o\
                        phys_param.o num_param.o radial.o logic.o\
                        output_data.o outRot.o integration.o useful.o\
                        mem_alloc.o legendre_spec_to_grid.o
integration.o         : cosine_transform_odd.o precision_mod.o constants.o
out_coeff.o           : logic.o precision_mod.o parallel.o blocking.o\
                        truncation.o communications.o LMLoop_data.o\
                        mem_alloc.o radial.o output_data.o phys_param.o
magic.o               : truncation.o num_param.o parallel.o logic.o\
                        output_data.o timing.o constants.o precision_mod.o\
                        useful.o startFields.o preCalculations.o\
                        step_time.o kinetic_energy.o mem_alloc.o\
			magnetic_energy.o horizontal.o fields_average.o \
			outGeos.o spectra.o out_TO.o output_data.o outMisc.o probes.o\
                        Namelists.o blocking.o out_coeff.o
step_time.o           : truncation.o blocking.o phys_param.o\
                        num_param.o radial_data.o mem_alloc.o\
                        logic.o output_data.o output.o movie.o\
                        char_manip.o timing.o courant.o\
                        radialLoop.o nonlinear_bcs.o\
                        LMLoop.o dt_fieldsLast.o \
			c_utils.o useful.o probes.o
phys_param.o          : precision_mod.o
char_manip.o          : precision_mod.o
constants.o           : precision_mod.o
outMisc.o             : truncation.o num_param.o phys_param.o logic.o\
                        output_data.o blocking.o constants.o radial.o\
                        horizontal.o legendre_spec_to_grid.o\
                        integration.o startFields.o
chebyshev_polynoms.o  : constants.o logic.o precision_mod.o
rIteration.o          : precision_mod.o
radialLoop.o          : truncation.o blocking.o parallel.o\
                        phys_param.o logic.o output_data.o\
                        constants.o legendre_spec_to_grid.o\
                        get_td.o out_graph_file.o radial_data.o rIteration.o \
			rIterThetaBlocking.o rIterThetaBlocking_OpenMP.o \
			rIterThetaBlocking_seq.o mem_alloc.o
rIterThetaBlocking.o  : rIteration.o truncation.o blocking.o \
			logic.o radial.o radial_data.o\
			$(FFT_OBJS) legendre_spec_to_grid.o get_nl.o\
                        legendre_helpers.o get_td.o nonlinear_bcs.o\
                        legendre_grid_to_spec.o  out_graph_file.o
rIterThetaBlocking_OpenMP.o : rIterThetaBlocking.o TO.o dtB.o\
                        out_movie_file.o outRot.o nl_special_calc.o courant.o\
			TO_arrays.o dtB_arrays.o probes.o
rIterThetaBlocking_seq.o : rIterThetaBlocking.o TO.o dtB.o\
                        out_movie_file.o outRot.o nl_special_calc.o\
			TO_arrays.o dtB_arrays.o probes.o
LMLoop.o              : truncation.o blocking.o parallel.o\
                        logic.o radial_data.o updateWPS.o\
                        output_data.o fields.o dt_fieldsLast.o\
                        updateS.o updateB.o mem_alloc.o updateXI.o\
                        updateZ.o updateWP.o timing.o LMLoop_data.o\
                        communications.o debugging.o updateWPT.o
debugging.o           : precision_mod.o
timing.o              : parallel.o precision_mod.o
LMLoop_data.o         : parallel.o blocking.o logic.o
legendre_spec_to_grid.o: truncation.o blocking.o horizontal.o logic.o\
                         constants.o legendre_helpers.o
radial_derivatives.o  : cosine_transform_odd.o constants.o precision_mod.o
legendre_helpers.o    : truncation.o num_param.o radial.o horizontal.o\
                        logic.o constants.o blocking.o special.o\
                        logic.o blocking.o fields.o TO.o mem_alloc.o
output.o              : truncation.o blocking.o phys_param.o\
                        num_param.o logic.o output_data.o radial.o\
                        horizontal.o constants.o fields.o radial_spectra.o\
                        out_movie_file.o outGeos.o\
                        outRot.o kinetic_energy.o magnetic_energy.o\
                        spectra.o fields_average.o power.o getDlm.o\
                        useful.o out_TO.o outMisc.o out_graph_file.o\
                        char_manip.o storeCheckPoints.o\
                        integration.o outOmega.o parallel.o \
                        communications.o out_coeff.o outPar.o dtB.o\
                        store_movie_IC.o RMS.o
communications.o      : blocking.o parallel.o LMLoop_data.o truncation.o\
			radial_data.o mem_alloc.o
outPar.o              : truncation.o blocking.o phys_param.o\
                        logic.o radial.o output_data.o radial_data.o\
                        horizontal.o fields.o legendre_spec_to_grid.o\
                        num_param.o integration.o constants.o mem_alloc.o
spectra.o             : truncation.o blocking.o num_param.o\
                        horizontal.o phys_param.o radial.o logic.o\
                        output_data.o useful.o LMLoop_data.o\
                        integration.o mem_alloc.o
RMS.o                 : truncation.o parallel.o blocking.o num_param.o\
                        phys_param.o radial.o constants.o logic.o\
                        output_data.o radial_derivatives.o useful.o\
                        RMS_helpers.o integration.o communications.o\
                        horizontal.o dtB.o chebyshev_polynoms.o mem_alloc.o\
                        finite_differences.o radial_scheme.o chebyshev.o
updateWP.o            : truncation.o blocking.o num_param.o integration.o\
                        phys_param.o radial.o horizontal.o logic.o\
                        RMS.o cosine_transform_odd.o radial_derivatives.o\
                        RMS_helpers.o LMLoop_data.o mem_alloc.o fields.o
updateWPT.o           : truncation.o blocking.o num_param.o fields.o\
                        phys_param.o radial.o horizontal.o logic.o\
                        cosine_transform_odd.o radial_derivatives.o\
                        LMLoop_data.o mem_alloc.o init_fields.o RMS.o
updateWPS.o           : truncation.o blocking.o num_param.o\
                        phys_param.o radial.o horizontal.o logic.o\
                        cosine_transform_odd.o radial_derivatives.o\
                        LMLoop_data.o mem_alloc.o init_fields.o RMS.o
chebInt.o             : chebyshev_polynoms.o cosine_transform_odd.o\
                        radial_derivatives.o
get_td.o              : truncation.o blocking.o horizontal.o\
                        phys_param.o num_param.o radial.o logic.o\
                        RMS.o RMS_helpers.o legendre_helpers.o fields.o\
                        cutils_iface.o mem_alloc.o
store_movie_IC.o      : truncation.o blocking.o logic.o radial_data.o\
                        movie.o radial.o horizontal.o\
                        legendre_helpers.o legendre_spec_to_grid.o\
                        $(FFT_OBJS) out_movie_file.o phys_param.o
radial_spectra.o      : truncation.o blocking.o horizontal.o radial.o\
                        num_param.o output_data.o logic.o useful.o\
                        char_manip.o radial_data.o LMmapping.o constants.o
nl_special_calc.o     : truncation.o blocking.o horizontal.o\
                        radial.o legendre_grid_to_spec.o radial_data.o\
                        phys_param.o constants.o
outOmega.o            : truncation.o blocking.o logic.o\
                        output_data.o radial.o cosine_transform_odd.o plms.o
fft_fac.o             : constants.o precision_mod.o
cosine_transform_even.o: truncation.o fft_fac.o constants.o useful.o\
                        mem_alloc.o
cosine_transform_odd.o : truncation.o fft_fac.o constants.o useful.o\
                        mem_alloc.o
readCheckPoints.o     : truncation.o blocking.o phys_param.o\
                        init_fields.o constants.o radial.o\
                        logic.o cosine_transform_odd.o finite_differences.o\
                        chebyshev.o radial_scheme.o
radial_derivatives_even.o : cosine_transform_odd.o precision_mod.o constants.o\
                        cosine_transform_even.o
updateZ.o             : truncation.o blocking.o phys_param.o\
                        num_param.o constants.o radial.o horizontal.o\
                        logic.o init_fields.o RMS.o TO.o\
                        $(LAPACK_OBJS) cosine_transform_odd.o\
                        radial_derivatives.o outRot.o RMS_helpers.o \
                        LMLoop_data.o communications.o mem_alloc.o special.o
out_graph_file.o      : truncation.o horizontal.o phys_param.o parallel.o\
                        num_param.o radial.o logic.o blocking.o\
                        output_data.o legendre_helpers.o $(FFT_OBJS)\
                        legendre_spec_to_grid.o
startFields.o         : truncation.o phys_param.o num_param.o constants.o\
                        radial.o init_fields.o blocking.o logic.o\
                        dt_fieldsLast.o special.o readCheckPoints.o\
                        fields.o radial_derivatives_even.o\
                        radial_derivatives.o useful.o communications.o
legendre_grid_to_spec.o : truncation.o horizontal.o blocking.o 
getDlm.o              : truncation.o blocking.o num_param.o constants.o\
                        horizontal.o radial.o constants.o LMLoop_data.o\
                        useful.o integration.o
out_dtB_frame.o       : truncation.o blocking.o radial.o logic.o\
                        horizontal.o dtB.o radial_derivatives.o\
                        radial_derivatives_even.o $(FFT_OBJS) constants.o
updateB.o             : truncation.o blocking.o horizontal.o fields.o\
                        phys_param.o num_param.o init_fields.o\
                        logic.o radial.o special.o constants.o\
                        RMS.o $(LAPACK_OBJS) radial_derivatives.o\
                        cosine_transform_odd.o radial_derivatives_even.o\
                        RMS_helpers.o LMLoop_data.o mem_alloc.o
outGeos.o             : truncation.o phys_param.o num_param.o\
                        horizontal.o radial.o blocking.o constants.o\
                        output_data.o logic.o cosine_transform_odd.o chebInt.o\
                        plms.o $(FFT_OBJS) LMLoop_data.o communications.o\
                        mem_alloc.o TO_helpers.o
preCalculations.o     : truncation.o phys_param.o num_param.o constants.o\
                        radial.o horizontal.o init_fields.o\
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
num_param.o           : truncation.o precision_mod.o mem_alloc.o
radial.o              : truncation.o radial_data.o $(LAPACK_OBJS)\
                        phys_param.o logic.o precision_mod.o finite_differences.o\
                        output_data.o cosine_transform_odd.o\
                        chebyshev_polynoms.o radial_derivatives.o\
			cosine_transform_even.o mem_alloc.o useful.o radial_scheme.o\
                        chebyshev.o finite_differences.o
output_data.o         : precision_mod.o
horizontal.o          : truncation.o phys_param.o num_param.o precision_mod.o\
                        radial.o logic.o blocking.o plms.o $(FFT_OBJS)\
                        mem_alloc.o
RMS.o	              : truncation.o blocking.o mem_alloc.o
init_fields.o         : truncation.o blocking.o radial.o horizontal.o\
                        constants.o logic.o $(FFT_OBJS) legendre_grid_to_spec.o\
                        useful.o phys_param.o LMLoop_data.o\
                        cosine_transform_odd.o mem_alloc.o
movie.o               : truncation.o parallel.o radial_data.o\
                        output_data.o logic.o radial.o mem_alloc.o\
                        horizontal.o char_manip.o useful.o
dtB.o                 : truncation.o radial_data.o parallel.o \
                        LMLoop_data.o $(FFT_OBJS) radial.o logic.o\
                        blocking.o horizontal.o phys_param.o\
                        communications.o legendre_grid_to_spec.o\
                        radial_spectra.o radial_derivatives.o mem_alloc.o
TO.o                  : truncation.o blocking.o horizontal.o\
                        radial_data.o phys_param.o radial.o logic.o\
                        legendre_grid_to_spec.o mem_alloc.o
$(LAPACK_OBJS)        : precision_mod.o constants.o
dt_fieldsLast.o       : truncation.o LMLoop_data.o precision_mod.o\
                        mem_alloc.o
fields.o              : truncation.o LMLoop_data.o radial_data.o\
                        precision_mod.o	mem_alloc.o logic.o
special.o             : precision_mod.o mem_alloc.o truncation.o
probes.o              : parallel.o precision_mod.o truncation.o radial.o num_param.o\
                        blocking.o horizontal.o output_data.o constants.o
radial_data.o         : truncation.o parallel.o logic.o
plms.o                : constants.o precision_mod.o
$(FFT_OBJS)           : truncation.o blocking.o useful.o constants.o
Namelists.o           : constants.o truncation.o phys_param.o num_param.o logic.o \
			init_fields.o output_data.o parallel.o \
			TO.o special.o radial.o char_manip.o probes.o

ifeq ($(USE_MPI),yes)
mpi.mod: $(MPI_INCPATH)/mpif.h
	$(FC) -c -o mpi.o -I$(MPI_INCPATH) mpimod.f90
endif

ifeq ($(USE_SHTNS),yes)
rIterThetaBlocking_shtns.o : rIterThetaBlocking.o TO.o dtB.o shtns.o\
                        out_movie_file.o outRot.o nl_special_calc.o courant.o\
			TO_arrays.o dtB_arrays.o probes.o
shtns.o               : truncation.o blocking.o horizontal.o radial.o  parallel.o\
                        precision_mod.o constants.o
radialLoop.o          : rIterThetaBlocking_shtns.o
magic.o               : shtns.o
endif

ifeq ($(USE_LAPACKLIB),MKL)
$(LAPACK_OBJS)        : lapack95.mod
lapack95.mod          : $(MKLROOT)/include/lapack.f90
	$(FC) -c $(FFLAGS) $(INCPATHS) $<
endif


.PHONY: clean echo
clean:
	@rm -f *.o $(OUT) *.mod
	-rm -f *__genmod.*

echo:
	@echo $(RED_SOURCES)
	@echo $(OBJS)
