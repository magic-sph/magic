# COMPILER can be intel, gnu or amd
COMPILER=intel

USE_PERFLIB=none
PRODRUN=yes
OPENMP=yes
# DEBUG = no noopt check all marmot
DEBUG=no
USE_MPI=yes

# USE_FFTLIB can be MKL or JW
USE_FFTLIB=JW
USE_MKL=no

# preconditioning of the matrices
# switches are WITH_PRECOND_S,WITH_PRECOND_BJ,
# WITH_PRECOND_Z
# PRECOND_WP is always set
# to switch the preconditioner off and get the old ill-conditioned matrices
# just comment the following statement
PREPROC += 	-DWITH_PRECOND_S \
		-DWITH_PRECOND_Z \
		-DWITH_PRECOND_BJ \
		-DWITH_PRECOND_S0 \
		-DWITH_PRECOND_Z10

ifeq ($(COMPILER),intel)
	FFLAG_TRACEBACK = -traceback
#	FFLAG_ARCH_OPT = -xsse3
	FFLAG_ARCH_OPT = -xAVX
	FFLAG_OPENMP = -openmp
	FFLAG_OPENMP_STUBS = -openmp-stubs
	FFLAG_CHECK_ALL = -check all
	FFLAG_FPE = -fpe0 -fp-stack-check -ftrapuv
	FFLAG_WARN_ALL  = -warn all,nointerfaces
	FFLAG_PREPROCESSOR = -fpp
	FFLAG_CONVERT = -convert big
	FFLAG_PROFILING = -prof-gen=srcpos
	FFLAG_REPORT = -vec-report6
	FFLAG_STD = -stand f03  -diag-disable 7025
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
	FFLAG_TRACEBACK = -fbacktrace
	FFLAG_OPENMP = -fopenmp
	FFLAG_CONVERT = -fconvert=big-endian
	FFLAG_STD = -std=f2003
	FFLAG_PREPROCESSOR = -cpp
	FFLAG_WARN_ALL = -Wall -Wno-maybe-uninitialized
	FFLAG_CHECK_ALL = -fcheck=all
	FFLAG_FPE = -ffpe-trap=invalid,zero,overflow 
	COMP_FC = gfortran
	COMP_CC = gcc
	COMP_MPFC = mpif90
	COMP_MPCC = mpicc

	FFLAG_MKL_LIBS = -L$(MKLROOT)/lib/intel64 -lmkl_gf_lp64 -lmkl_sequential -lmkl_core -lpthread -lm \
		-Wl,-rpath,$(MKLROOT)/lib/intel64
endif
ifeq ($(COMPILER),amd)
	FFLAGS = -m64
	FFLAG_TRACEBACK = 
	FFLAG_ARCH_OPT = -march=barcelona
	FFLAG_OPENMP = -mp
	FFLAG_CONVERT = -convert big_endian
	FFLAG_PREPROCESSOR = -cpp
	FFLAG_WARN_ALL = -Wall
	FC = openf95

	FFLAG_MKL_LIBS = -L$(MKLROOT)/lib/intel64 -lmkl_gf_lp64 -lmkl_sequential -lmkl_core -lpthread -lm \
		-Wl,-rpath,$(MKLROOT)/lib/intel64
endif

# the following path is only used, if the mpi.mod file is
# not found in the standard location. It is used for example
# on hydra, if one compiles with the gnu compilers.
MPI_INCPATH=/opt/ibmhpc/pecurrent/mpich2/gnu/include64
MPI_INCPATH=/usr/include

PREPROC += -DJW=1 -DMKL=2
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

ifeq ($(USE_FFTLIB),JW)
	PREPROC += -DFFTLIB=JW
endif
ifeq ($(USE_FFTLIB),MKL)
	PREPROC += -DFFTLIB=MKL
	INCPATHS += -I$(MKLROOT)/include
	LIBS +=  $(FFLAG_MKL_LIBS)
endif

ifeq ($(USE_MKL),yes)
	PREPROC += -DWITH_MKL_LU
	INCPATHS += -I$(MKLROOT)/include
	LIBS += $(FFLAG_MKL_LIBS)
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
ifeq ($(OPENMP),yes)
	FFLAGS += $(FFLAG_OPENMP)
	PREPROC += -DWITHOMP
endif
ifeq ($(OPENMP),stubs)
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
SOURCES:= $(wildcard [ms]*.F90)
RED_SOURCES = $(filter-out m_fft% s_graphOut_mpi% s_graphOut_IC_mpi% m_store_mpi% mpimod%,$(SOURCES))
#RED_SOURCES = $(filter-out m_fft% s_graphOut_mpi%,$(SOURCES))

ifeq ($(USE_FFTLIB),JW)
	RED_SOURCES += m_fftJW.F90
	FFT_OBJS += m_fftJW.o
endif
ifeq ($(USE_FFTLIB),MKL)
	RED_SOURCES += m_fft_mkl.F90
	FFT_OBJS += m_fft_mkl.o
endif

ifeq ($(USE_MPI),yes)
	RED_SOURCES += s_graphOut_mpi.F90 m_store_mpi.F90
endif

OBJS   := $(addsuffix .o, $(basename $(RED_SOURCES)))
OBJS += truncation.o
OBJS += c_utils.o

.SUFFIXES: 

ifeq ($(COMPILER),gnu)
all: mpi.mod
	make $(OUT)
endif

$(OUT): $(OBJS)
	$(FC) $(FFLAGS) -o $@ $^ $(LIBS)

%.o: %.F90
	$(FC) -c $(FFLAGS) $(PREPROC) $(INCPATHS) -o $@ $<

s_getStartFields.o: s_getStartFields.F90
	$(FC) -c $(FFLAGS_NOCHECK) $(PREPROC) $(INCPATHS) -o $@ $<

m_timing.o: m_timing.F90
	$(FC) -c $(FFLAGS_NOCHECK) $(PREPROC) $(INCPATHS) -o $@ $<

%.o: %.c
	$(CC) -c $(CFLAGS) $(PREPROC) $(INCPATHS) -o $@ $<

# cancelling the implicit rule for .F files
#%.f: %.F

s_get_B_surface.o       : truncation.o m_blocking.o m_radial.o m_horizontal.o\
                          $(FFT_OBJS) 
m_magnetic_energy.o     : truncation.o m_blocking.o m_horizontal.o\
                          m_phys_param.o m_num_param.o m_const.o m_radial.o\
                          m_logic.o m_output_data.o m_movie.o m_Bext.o\
                          m_integration.o m_usefull.o
s_legtfAS.o             : truncation.o m_horizontal.o 
s_getTO.o               : truncation.o m_blocking.o m_horizontal.o\
                          m_num_param.o m_phys_param.o m_radial.o m_logic.o\
                          m_output_data.o m_TO.o s_legtfAS2.o 
s_init_costf2.o         : m_usefull.o 
s_legtfAS2.o            : truncation.o m_horizontal.o 
s_preCalcTimes.o        : truncation.o m_logic.o m_num_param.o m_output_data.o\
                          m_init_fields.o s_get_hit_times.o 
m_outTO.o               : truncation.o m_horizontal.o m_num_param.o\
                          m_phys_param.o m_radial.o m_logic.o m_blocking.o\
                          m_const.o m_output_data.o m_TO.o s_legtfAS.o s_legtfAS2.o\
                          s_costf1.o s_chebInt.o m_plms.o s_getPAStr.o\
                          s_getAStr.o s_get_PAS.o m_char_manip.o m_integration.o
s_courant.o             : truncation.o m_blocking.o m_num_param.o\
                          m_phys_param.o m_radial.o m_horizontal.o m_logic.o 
s_get_br_v_bcs.o        : truncation.o m_horizontal.o $(FFT_OBJS)\
                          s_legTF2.o 
m_fields_average.o      : truncation.o m_num_param.o m_output_data.o m_blocking.o\
                          m_logic.o m_const.o m_radial.o m_horizontal.o\
                          s_get_drNS.o s_get_ddrNS_even.o s_get_drNS_even.o\
                          s_spectrum.o s_spectrumC.o m_kinetic_energy.o\
                          m_magnetic_energy.o s_graphOut.o s_legPrep.o s_legTF.o\
                          $(FFT_OBJS) m_write_special.o s_storePotW.o 
ifeq ($(USE_MPI),yes)
m_fields_average.o	: s_graphOut_mpi.o
endif
s_storePotW.o           : truncation.o m_radial.o m_output_data.o m_phys_param.o\
                          m_horizontal.o m_logic.o m_char_manip.o s_costf1.o 
s_hInt2dPol.o           : truncation.o m_blocking.o m_radial.o m_horizontal.o\
                          m_usefull.o
s_get_nl.o              : truncation.o m_horizontal.o m_num_param.o\
                          m_phys_param.o m_radial.o m_logic.o m_blocking.o \
			  m_leg_helper.o
s_store_movie_frame.o   : truncation.o m_output_data.o m_blocking.o\
                          m_output_data.o m_radial.o m_horizontal.o m_movie.o\
                          s_get_B_surface.o m_output_data.o m_movie.o\
                          m_horizontal.o m_num_param.o \
                          m_output_data.o m_movie.o m_horizontal.o\
                          m_phys_param.o s_get_fl.o\
                          s_get_sl.o m_output_data.o m_movie.o\
                          m_horizontal.o \
                          m_output_data.o m_movie.o m_horizontal.o
s_store.o               : truncation.o m_blocking.o m_phys_param.o\
                          m_num_param.o m_radial.o m_logic.o m_output_data.o\
                          m_init_fields.o m_dt_fieldsLast.o 
m_updateS.o             : truncation.o m_blocking.o m_horizontal.o\
                          m_phys_param.o m_radial.o m_logic.o m_output_data.o\
                          m_init_fields.o m_const.o m_mat.o s_get_dr.o\
                          s_get_s0Mat.o m_algebra.o s_get_sMat.o s_costf1.o\
                          s_get_ddr.o m_LMLoop_data.o
s_getTOfinish.o         : truncation.o m_blocking.o m_horizontal.o\
                          m_num_param.o m_phys_param.o m_radial.o m_logic.o\
                          m_output_data.o m_TO.o 
m_power.o               : truncation.o m_blocking.o m_horizontal.o\
                          m_phys_param.o m_num_param.o m_radial.o m_logic.o\
                          m_output_data.o m_movie.o s_get_viscous_torque.o\
                          m_integration.o m_usefull.o
m_integration.o         : s_costf1.o 
s_get_dtBLM.o           : truncation.o m_horizontal.o m_num_param.o\
                          m_phys_param.o m_radial.o m_logic.o m_blocking.o\
                          $(FFT_OBJS) s_legTF3.o s_legTF2.o 
s_mapDataIC.o           : truncation.o m_radial.o m_init_fields.o m_blocking.o\
                          s_mapDataICR.o 
s_init_costf1.o         : m_usefull.o 
s_get_dr.o              : s_costf1.o s_get_dcheb.o 
s_get_b_nl_bcs.o        : truncation.o m_blocking.o m_horizontal.o\
                          m_phys_param.o 
s_safeOut.o             : m_logic.o
s_initB.o               : truncation.o m_blocking.o m_phys_param.o\
                          m_num_param.o m_radial.o m_horizontal.o m_logic.o\
                          m_init_fields.o m_Grenoble.o m_const.o s_j_cond.o\
                          m_usefull.o
m_outPV3.o              : truncation.o m_horizontal.o m_num_param.o\
                          m_phys_param.o m_radial.o m_logic.o m_blocking.o\
                          m_output_data.o s_costf1.o m_plms.o s_getPAStr.o\
                          s_getPVptr.o 
m_write_special.o       : m_logic.o
s_graph_write.o         : truncation.o 
s_lm2pt.o               : truncation.o m_horizontal.o m_blocking.o m_radial.o\
                          $(FFT_OBJS) 
magic.o                 : truncation.o m_num_param.o m_parallel.o m_logic.o\
                          m_output_data.o m_timing.o s_openFiles.o\
                          s_logWrite.o s_preCalc.o s_getStartFields.o\
                          s_preCalcTimes.o s_writeInfo.o s_checkMovie.o\
                          m_step_time.o s_closeFiles.o m_kinetic_energy.o \
			  m_magnetic_energy.o m_horizontal.o m_fields_average.o \
			  m_Egeos.o m_spectrum_average.o m_spectrumC_average.o \
			  m_outTO.o m_output_data.o m_outPV3.o m_Namelists.o m_blocking.o
m_step_time.o           : truncation.o m_blocking.o m_phys_param.o\
                          m_num_param.o m_init_fields.o m_const.o m_radial.o\
                          m_horizontal.o m_logic.o m_output_data.o m_output.o m_movie.o\
                          m_char_manip.o m_timing.o s_logWrite.o\
                          s_safeOut.o m_radialLoop.o s_get_b_nl_bcs.o\
                          s_check_time_hits.o m_LMLoop.o\
			  m_radialLoop.o m_dt_fieldsLast.o \
			  c_utils.o
s_get_PAS.o             : truncation.o m_blocking.o m_num_param.o m_radial.o\
                          m_horizontal.o m_logic.o 
s_get_wpMat.o           : truncation.o m_blocking.o m_phys_param.o\
                          m_num_param.o m_radial.o m_logic.o m_algebra.o 
s_outMisc.o             : truncation.o m_num_param.o m_phys_param.o m_logic.o\
                          m_output_data.o m_blocking.o m_const.o m_radial.o\
                          m_horizontal.o s_lmAS2pt.o m_Egeos.o m_integration.o
m_chebyshev_polynoms.o  : m_const.o
s_spectrumC.o           : truncation.o m_blocking.o m_num_param.o\
                          m_horizontal.o m_phys_param.o m_radial.o m_logic.o\
                          m_output_data.o m_const.o m_integration.o
m_radialLoop.o          : truncation.o m_blocking.o m_parallel.o\
                          m_phys_param.o m_num_param.o m_logic.o m_output_data.o\
                          m_radial.o m_horizontal.o m_const.o m_TO.o \
                          s_legPrepG.o m_legendre_trafo.o\
                          $(FFT_OBJS) s_v_rigid_boundary.o \
                          s_legTF3.o s_legTF2.o s_legTF1.o s_get_br_v_bcs.o\
                          s_get_td.o s_graphOut.o m_radial_data.o m_rIteration.o \
			  m_rIterThetaBlocking.o m_rIterThetaBlocking_OpenMP.o \
			  m_rIterThetaBlocking_seq.o
m_rIterThetaBlocking.o  : m_rIteration.o truncation.o m_blocking.o \
			  m_logic.o m_radial.o m_radial_data.o m_output_data.o \
			  $(FFT_OBJS) m_legendre_trafo.o s_get_nl.o m_leg_helper.o \
			  m_grid_space_arrays.o m_nonlinear_lm.o
m_rIterThetaBlocking_OpenMP.o : m_rIterThetaBlocking.o
m_rIterThetaBlocking_seq.o : m_rIterThetaBlocking.o
m_grid_space_arrays.o   : m_blocking.o truncation.o
ifeq ($(USE_MPI),yes)
m_radialLoop.o          : s_graphOut_mpi.o
endif
m_LMLoop.o              : truncation.o m_blocking.o m_parallel.o\
                          m_phys_param.o m_logic.o m_radial.o m_horizontal.o\
                          m_output_data.o m_fields.o m_dt_fieldsLast.o m_mat.o\
			  m_updateS.o m_updateB.o \
                          m_updateZ.o m_updateWP.o m_timing.o m_LMLoop_data.o\
			  m_communications.o m_debugging.o
m_timing.o		: m_parallel.o
m_LMLoop_data.o		: m_parallel.o m_blocking.o
s_get_ddrNS_even.o      : s_costf1.o s_get_ddcheb_even.o s_costf2.o 
s_get_movie_type.o      : truncation.o m_output_data.o m_phys_param.o m_radial.o\
                          m_logic.o m_movie.o m_horizontal.o m_char_manip.o
s_dt_courant.o          : truncation.o m_output_data.o m_logic.o  m_radial_data.o
m_legendre_trafo.o      : truncation.o m_blocking.o m_horizontal.o m_logic.o\
                          m_const.o m_leg_helper.o
s_get_s0Mat.o           : truncation.o m_phys_param.o m_num_param.o m_radial.o\
                          m_algebra.o 
s_get_dddr.o            : s_costf1.o s_get_dddcheb.o 
s_openFiles.o           : truncation.o m_logic.o m_output_data.o m_init_fields.o\
                          m_parallel.o m_char_manip.o
s_writeInfo.o           : truncation.o m_phys_param.o m_num_param.o m_logic.o\
                          m_init_fields.o m_output_data.o m_const.o 
m_spectrumC_average.o   : truncation.o m_blocking.o m_num_param.o\
                          m_horizontal.o m_phys_param.o m_radial.o m_logic.o\
                          m_output_data.o m_const.o s_safeOut.o\
                          m_integration.o m_usefull.o
s_get_dtB.o             : truncation.o m_blocking.o m_radial.o m_horizontal.o\
                          m_logic.o 
s_zeroRms.o             : truncation.o m_blocking.o m_RMS.o 
s_get_RAS.o             : truncation.o m_blocking.o m_num_param.o m_radial.o\
                          m_horizontal.o m_logic.o 
s_legPrepG.o            : truncation.o m_num_param.o m_radial.o m_horizontal.o\
                          m_logic.o m_const.o m_blocking.o m_Grenoble.o\
                          m_fields.o m_TO.o m_leg_helper.o
m_output.o              : truncation.o m_blocking.o m_phys_param.o\
                          m_num_param.o m_logic.o m_output_data.o m_radial.o\
                          m_horizontal.o m_const.o m_fields.o s_rBrSpec.o\
                          s_rBpSpec.o s_logWrite.o s_write_movie_frame.o\
                          s_write_rot.o m_kinetic_energy.o m_magnetic_energy.o\
                          m_spectrum_average.o m_spectrumC_average.o\
                          m_fields_average.o m_power.o m_getDlm.o\
                          s_get_u_square.o s_safeOut.o m_outTO.o s_outMisc.o\
                          s_dtVrms.o s_zeroRms.o m_char_manip.o s_store.o\
                          m_integration.o m_outPV3.o m_outOmega.o m_parallel.o \
			  m_communications.o m_write_special.o m_outPar.o\
                          m_outPerpPar.o
m_communications.o	: m_blocking.o m_parallel.o m_LMLoop_data.o truncation.o
ifeq ($(USE_MPI),yes)
m_output.o		: m_store_mpi.o
endif
m_outPar.o              : truncation.o m_blocking.o m_phys_param.o\
                          m_logic.o m_radial.o m_output_data.o m_radial_data.o\
                          m_horizontal.o m_const.o m_fields.o s_lmAS2pt.o
m_outPerpPar.o          : truncation.o m_blocking.o m_horizontal.o m_const.o\
                          m_radial.o m_output_data.o s_lmAS2pt.o m_radial_data.o\
                          m_parallel.o m_integration.o m_logic.o m_num_param.o
s_cheb_x_map_e_xr.o     : m_logic.o 
s_mapDataICR.o          : truncation.o m_blocking.o m_radial.o s_init_costf1.o 
m_spectrum_average.o    : truncation.o m_blocking.o m_num_param.o\
                          m_horizontal.o m_phys_param.o m_radial.o m_logic.o\
                          m_output_data.o s_safeOut.o m_usefull.o m_LMLoop_data.o\
                          m_integration.o
s_dtVrms.o              : truncation.o m_parallel.o m_blocking.o\
                          m_phys_param.o m_radial.o m_const.o m_logic.o\
                          m_output_data.o m_RMS.o s_init_rNB.o s_get_drNS_R.o\
                          s_hInt2dPol.o m_integration.o m_communications.o
s_get_drNS.o            : s_costf1.o s_get_dcheb.o 
m_updateWP.o            : truncation.o m_blocking.o m_num_param.o\
                          m_phys_param.o m_radial.o m_horizontal.o m_logic.o\
                          m_mat.o m_RMS.o s_get_wpMat.o s_costf1.o\
                          s_get_dddr.o s_get_dr.o s_hInt2Pol.o m_LMLoop_data.o
s_get_sMat.o            : truncation.o m_phys_param.o m_num_param.o m_radial.o\
                          m_algebra.o 
s_get_z10Mat.o          : truncation.o m_phys_param.o m_num_param.o m_const.o\
                          m_radial.o m_logic.o m_algebra.o 
s_chebInt.o             : s_cheb_x_map_e.o s_init_costf1.o s_get_dcheb.o 
s_get_dr_even.o         : s_costf1.o s_get_dcheb_even.o s_costf2.o 
s_parallel.o            : m_parallel.o 
s_get_td.o              : truncation.o m_blocking.o m_horizontal.o\
                          m_phys_param.o m_num_param.o m_radial.o m_logic.o\
                          m_RMS.o s_hInt2Pol.o s_hInt2Tor.o m_leg_helper.o \
			  m_nonlinear_lm.o
m_nonlinear_lm.o	: m_cutils_iface.o
s_get_fl.o              : truncation.o m_blocking.o m_radial.o m_fields.o\
                          m_horizontal.o m_logic.o m_Grenoble.o
s_store_movie_frame_IC.o: truncation.o m_blocking.o m_logic.o m_output_data.o\
                          m_movie.o m_radial.o m_horizontal.o\
                          m_phys_param.o s_legPrep_IC.o s_legTF.o\
                          $(FFT_OBJS) s_get_fl.o 
s_legTF2.o              : truncation.o m_horizontal.o m_blocking.o 
s_rBrSpec.o             : truncation.o m_blocking.o m_horizontal.o m_radial.o\
                          m_num_param.o m_phys_param.o m_output_data.o m_logic.o\
                          m_char_manip.o
s_get_ddr.o             : s_costf1.o s_get_ddcheb.o 
s_check_time_hits.o     : truncation.o m_output_data.o m_num_param.o m_logic.o 
s_graphOut_IC.o         : truncation.o m_blocking.o m_num_param.o m_radial.o\
                          m_horizontal.o m_logic.o m_output_data.o s_legPrep_IC.o\
                          s_legTF.o $(FFT_OBJS) s_graph_write.o 
ifeq ($(USE_MPI),yes)
s_graphOut_IC.o		: s_graphOut_mpi.o
endif
s_get_Bpol.o            : truncation.o m_blocking.o m_num_param.o m_radial.o\
                          m_horizontal.o m_logic.o $(FFT_OBJS) 
s_getHelLM.o            : truncation.o m_blocking.o m_horizontal.o\
                          m_num_param.o m_phys_param.o m_radial.o m_logic.o\
                          m_output_data.o s_legtfAS2.o 
s_get_nlBLayers.o       : truncation.o m_blocking.o m_horizontal.o\
                          m_radial.o s_legtfAS2.o s_legtfAS.o
s_get_fluxes.o          : truncation.o m_blocking.o m_horizontal.o\
                          m_radial.o s_legtfAS2.o s_legtfAS.o m_const.o
s_get_perpPar.o         : truncation.o m_blocking.o m_horizontal.o\
                          m_radial.o s_legtfAS2.o
s_hInt2Pol.o            : truncation.o m_blocking.o m_const.o m_radial.o\
                          m_horizontal.o m_usefull.o
s_write_movie_frame.o   : truncation.o m_blocking.o m_radial.o m_horizontal.o\
                          m_phys_param.o m_num_param.o m_output_data.o m_logic.o\
                          m_movie.o s_write_dtB_frame.o 
m_outOmega.o            : truncation.o m_blocking.o m_logic.o m_const.o\
                          m_output_data.o m_radial.o s_costf1.o m_plms.o
s_storePot.o            : truncation.o m_radial.o m_output_data.o m_phys_param.o\
                          m_horizontal.o m_logic.o m_char_manip.o s_costf1.o 
s_costf2.o              : s_fft_fac.o truncation.o
s_checkMovie.o          : truncation.o m_logic.o m_output_data.o m_parallel.o\
                          s_get_movie_type.o s_logWrite.o 
s_v_rigid_boundary.o    : truncation.o m_blocking.o m_radial.o m_horizontal.o 
s_readStartFields.o     : truncation.o m_blocking.o m_phys_param.o\
                          m_init_fields.o m_const.o m_logic.o m_radial.o\
                          m_output_data.o s_mapData.o s_mapDataIC.o 
s_get_ddr_even.o        : s_costf1.o s_get_ddcheb_even.o s_costf2.o 
s_legTF3.o              : truncation.o m_horizontal.o m_blocking.o 
s_hInt2Tor.o            : truncation.o m_blocking.o m_radial.o m_horizontal.o 
m_updateZ.o             : truncation.o m_blocking.o m_phys_param.o\
                          m_num_param.o m_const.o m_radial.o m_horizontal.o\
                          m_logic.o m_init_fields.o m_mat.o m_RMS.o m_TO.o\
                          s_get_z10Mat.o m_algebra.o s_get_zMat.o s_costf1.o\
                          s_get_ddr.o s_get_angular_moment.o s_hInt2Tor.o \
			  m_LMLoop_data.o m_communications.o
s_closeFiles.o          : truncation.o m_logic.o m_output_data.o m_parallel.o 
s_getPVptr.o            : truncation.o m_blocking.o m_radial.o m_horizontal.o\
                          m_const.o $(FFT_OBJS) 
s_initS.o               : truncation.o m_blocking.o m_phys_param.o\
                          m_num_param.o m_radial.o m_const.o m_horizontal.o\
                          m_init_fields.o s_s_cond.o $(FFT_OBJS) s_legTF1.o\
                          m_algebra.o m_usefull.o
s_get_drNS_even.o       : s_costf1.o s_get_dcheb_even.o s_costf2.o 
s_initV.o               : truncation.o m_blocking.o m_phys_param.o\
                          m_num_param.o m_radial.o m_horizontal.o\
                          m_init_fields.o m_const.o m_logic.o $(FFT_OBJS)\
                          s_legTF1.o m_usefull.o
s_dtBrms.o              : truncation.o m_parallel.o m_blocking.o m_num_param.o\
                          m_output_data.o m_radial.o m_horizontal.o m_const.o\
                          m_logic.o m_phys_param.o m_dtB.o m_RMS.o\
                          s_get_drNS.o s_get_PolTorRms.o s_hInt2dPol.o\
                          s_get_RAS.o s_get_PASLM.o m_integration.o \
			  m_communications.o
s_s_cond.o              : truncation.o m_horizontal.o m_phys_param.o\
                          m_radial.o m_mat.o m_init_fields.o\
                          m_algebra.o s_costf1.o 
s_graphOut.o            : truncation.o m_horizontal.o m_phys_param.o\
                          m_num_param.o m_radial.o m_logic.o m_blocking.o\
                          m_output_data.o s_graph_write.o 
s_graphOut_mpi.o        : truncation.o m_horizontal.o m_phys_param.o\
                          m_num_param.o m_radial.o m_logic.o m_blocking.o\
                          m_output_data.o s_graph_write.o m_parallel.o
s_get_angular_moment.o  : truncation.o m_blocking.o m_const.o m_radial.o\
                          m_integration.o
s_mapData.o             : truncation.o m_init_fields.o m_radial.o m_blocking.o\
                          m_logic.o s_mapDataR.o 
s_getStartFields.o      : truncation.o m_phys_param.o m_num_param.o m_const.o\
                          m_radial.o m_horizontal.o m_init_fields.o\
                          m_blocking.o m_logic.o m_output_data.o m_fields.o\
                          m_dt_fieldsLast.o m_Grenoble.o s_readStartFields.o\
                          s_logWrite.o s_s_cond.o s_get_dr.o s_initB.o\
                          s_get_ddr.o s_get_ddr_even.o m_usefull.o m_communications.o
s_rBpSpec.o             : truncation.o m_blocking.o m_horizontal.o m_radial.o\
                          m_num_param.o m_phys_param.o m_output_data.o m_logic.o\
                          m_char_manip.o m_usefull.o
s_legTF1.o              : truncation.o m_horizontal.o m_blocking.o 
s_legTF.o		: truncation.o
s_logWrite.o            : truncation.o m_output_data.o m_parallel.o m_logic.o 
s_mapDataR.o            : truncation.o m_blocking.o m_horizontal.o m_radial.o\
                          s_init_costf1.o 
s_lmAS2pt.o             : truncation.o m_blocking.o m_horizontal.o m_radial.o 
s_get_dtBLMfinish.o     : truncation.o m_horizontal.o m_phys_param.o\
                          m_radial.o m_logic.o m_blocking.o m_dtB.o\
                          s_get_drNS.o s_rBrSpec.o s_rBpSpec.o 
m_getDlm.o              : truncation.o m_blocking.o m_num_param.o\
                          m_horizontal.o m_phys_param.o m_radial.o m_const.o\
                          m_usefull.o m_integration.o
s_costf1.o              : s_fft_fac.o truncation.o
s_write_dtB_frame.o     : truncation.o m_output_data.o m_blocking.o m_radial.o\
                          m_horizontal.o m_dtB.o s_get_dtB.o s_get_drNS.o\
                          s_get_Bpol.o s_get_Btor.o s_lm2pt.o s_get_drNS_even.o 
s_spectrum.o            : truncation.o m_blocking.o m_num_param.o\
                          m_horizontal.o m_phys_param.o m_radial.o m_logic.o\
                          m_output_data.o  m_usefull.o m_integration.o
s_get_PolTorRms.o       : truncation.o m_blocking.o m_const.o m_radial.o\
                          m_horizontal.o m_usefull.o m_integration.o
m_updateB.o             : truncation.o m_blocking.o m_horizontal.o\
                          m_phys_param.o m_num_param.o m_init_fields.o\
                          m_logic.o m_radial.o m_Bext.o m_const.o m_mat.o\
                          m_RMS.o s_get_drNS.o s_get_bMat.o m_algebra.o\
                          s_costf1.o s_get_ddr.o s_get_ddr_even.o s_hInt2Pol.o\
                          s_hInt2Tor.o m_LMLoop_data.o
s_getTOnext.o           : truncation.o m_blocking.o m_horizontal.o\
                          m_num_param.o m_phys_param.o m_radial.o m_logic.o\
                          m_output_data.o m_TO.o 
m_Egeos.o               : truncation.o m_phys_param.o m_num_param.o\
                          m_horizontal.o m_radial.o m_blocking.o m_const.o\
                          m_output_data.o m_logic.o s_costf1.o s_chebInt.o\
                          m_plms.o s_getDVptr.o m_LMLoop_data.o
s_get_Btor.o            : truncation.o m_blocking.o m_num_param.o m_radial.o\
                          m_horizontal.o m_logic.o $(FFT_OBJS) 
s_preCalc.o             : truncation.o m_phys_param.o m_num_param.o m_const.o\
                          m_radial.o m_horizontal.o m_init_fields.o\
                          m_blocking.o m_logic.o m_output_data.o s_get_dr.o\
                          s_logWrite.o m_integration.o
s_get_sl.o              : truncation.o m_blocking.o m_radial.o m_fields.o\
                          m_horizontal.o m_logic.o 
s_init_rNB.o            : s_cheb_x_map_e.o s_init_costf1.o s_get_dr.o 
s_get_u_square.o        : truncation.o m_blocking.o m_num_param.o m_radial.o\
                          m_horizontal.o m_output_data.o m_logic.o m_const.o\
                          m_phys_param.o m_integration.o
s_get_zMat.o            : truncation.o m_blocking.o m_phys_param.o\
                          m_num_param.o m_radial.o m_logic.o m_algebra.o 
s_get_lorentz_torque.o  : truncation.o m_radial.o m_horizontal.o m_Grenoble.o\
                          m_phys_param.o 
s_j_cond.o              : truncation.o m_horizontal.o m_phys_param.o\
                          m_radial.o m_init_fields.o m_logic.o m_mat.o\
                          m_algebra.o s_costf1.o 
s_get_PASLM.o           : truncation.o m_blocking.o m_num_param.o m_radial.o\
                          m_horizontal.o m_logic.o 
m_kinetic_energy.o      : truncation.o m_blocking.o m_num_param.o m_radial.o\
                          m_horizontal.o m_output_data.o m_logic.o m_const.o\
                          m_phys_param.o m_usefull.o m_integration.o\
                          m_communications.o
s_get_bMat.o            : truncation.o m_blocking.o m_phys_param.o\
                          m_num_param.o m_radial.o m_logic.o m_init_fields.o\
                          m_Bext.o m_algebra.o 
s_get_dH_dtBLM.o        : truncation.o m_blocking.o m_horizontal.o\
                          m_phys_param.o m_num_param.o m_radial.o m_logic.o\
                          m_dtB.o 
s_getB0.o               : truncation.o m_blocking.o m_horizontal.o m_radial.o\
                          m_num_param.o m_const.o m_phys_param.o m_output_data.o\
                          m_logic.o m_Grenoble.o m_usefull.o m_integration.o
s_getDVptr.o            : truncation.o m_blocking.o m_radial.o m_horizontal.o\
                          m_const.o $(FFT_OBJS) 
s_write_rot.o           : truncation.o m_blocking.o m_num_param.o\
                          m_phys_param.o m_logic.o m_radial.o m_output_data.o\
                          m_const.o m_horizontal.o s_get_viscous_torque.o\
                          s_get_angular_moment.o m_integration.o
m_blocking.o            : truncation.o m_LMmapping.o m_logic.o m_parallel.o m_output_data.o
m_num_param.o           : truncation.o
m_radial.o              : truncation.o m_radial_data.o m_mat.o m_algebra.o\
			  m_num_param.o m_phys_param.o m_logic.o\
                          m_output_data.o s_cheb_x_map_e_xr.o m_init_fields.o\
                          m_chebyshev_polynoms.o s_init_costf1.o s_get_dr.o\
                          s_get_chebs_even.o s_init_costf2.o
m_output_data.o         : truncation.o m_movie.o
m_horizontal.o          : truncation.o m_phys_param.o m_num_param.o\
                          m_radial.o m_logic.o m_output_data.o\
                          m_blocking.o m_plms.o $(FFT_OBJS)
m_RMS.o	                : truncation.o m_blocking.o
m_Grenoble.o            : truncation.o
m_init_fields.o         : truncation.o
m_movie.o               : truncation.o m_parallel.o m_radial_data.o
m_dtB.o                 : truncation.o m_radial_data.o m_parallel.o \
			  m_LMLoop_data.o
m_TO.o                  : truncation.o
m_mat.o                 : truncation.o
m_dt_fieldsLast.o       : truncation.o m_LMLoop_data.o
m_fields.o              : truncation.o m_LMLoop_data.o m_radial_data.o
m_radial_data.o		: truncation.o m_parallel.o m_logic.o

m_fftJW.o               : truncation.o m_blocking.o m_usefull.o
m_fft_mkl.o             : truncation.o m_blocking.o
s_legTF.o               : truncation.o
m_Namelists.o           : m_const.o truncation.o m_phys_param.o m_num_param.o m_logic.o \
			  m_init_fields.o m_output_data.o m_parallel.o \
			  m_TO.o m_Grenoble.o m_radial.o m_Bext.o m_char_manip.o
m_store_mpi.o		: truncation.o m_output_data.o m_parallel.o m_logic.o m_dt_fieldsLast.o

mpi.mod: $(MPI_INCPATH)/mpif.h
	$(FC) -c -o mpi.o -I$(MPI_INCPATH) mpimod.F90


ifeq ($(USE_MKL),yes)
s_initS.o		: lapack95.mod
m_radial.o		: lapack95.mod
s_s_cond.o		: lapack95.mod
s_j_cond.o		: lapack95.mod
s_get_bMat.o		: lapack95.mod
m_updateB.o		: lapack95.mod
s_get_z10Mat.o		: lapack95.mod
s_get_zMat.o		: lapack95.mod
m_updateZ.o		: lapack95.mod
s_get_wpMat.o		: lapack95.mod
m_updateWP.o		: lapack95.mod
m_updateS.o		: lapack95.mod
s_get_s0Mat.o		: lapack95.mod
s_get_sMat.o		: lapack95.mod
lapack95.mod		: $(MKLROOT)/include/lapack.f90
	$(FC) -c $(FFLAGS) $(INCPATHS) $<
endif


.PHONY: clean echo
clean:
	@rm -f *.o $(OUT) *.mod
	-rm -f *__genmod.*

echo:
	@echo $(RED_SOURCES)
	@echo $(OBJS)

