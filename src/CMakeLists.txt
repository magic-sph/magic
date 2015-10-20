add_library(truncation truncation.f90)
add_library(num_param num_param.f90)
add_library(precision_mod precision_mod.f90)
set_source_files_properties(mpimod.f90 PROPERTIES COMPILE_FLAGS ${FFLAG_NO_STD})
add_library(mpimod mpimod.f90)
add_library(Bext Bext.f90)
add_library(debugging debugging.f90)
add_library(phys_param phys_param.f90)
add_library(char_manip char_manip.f90)
add_library(constants constants.f90)
add_library(parallel parallel.f90)
add_library(logic logic.f90)
add_library(useful useful.f90)
add_library(LMmapping LMmapping.f90)
add_library(blocking blocking.f90)
add_library(RMS RMS.f90)
add_library(RMS_helpers RMS_helpers.f90)
add_library(Grenoble Grenoble.f90)
add_library(timing timing.f90)
add_library(radial_data radial_data.f90)
add_library(plms plms.f90)
add_library(LMLoop_data LMLoop_data.f90)
add_library(radial radial.f90)
add_library(horizontal horizontal.f90)
add_library(mat mat.f90)
add_library(fields fields.f90)
add_library(dt_fieldsLast dt_fieldsLast.f90)
add_library(radial_derivatives radial_derivatives.f90)
add_library(radial_derivatives_even radial_derivatives_even.f90)
add_library(chebyshev_polynoms chebyshev_polynoms.f90)
add_library(init_costf init_costf.f90)
add_library(cosine_transform cosine_transform.f90)
add_library(chebInt chebInt.f90)
add_library(fft_fac fft_fac.f90)
add_library(Namelists Namelists.f90)
add_library(TO TO.f90)
add_library(TO_helpers TO_helpers.f90)
add_library(legendre_grid_to_spec legendre_grid_to_spec.f90)
add_library(legendre_helpers legendre_helpers.f90)
add_library(legendre_spec_to_grid legendre_spec_to_grid.f90)
add_library(init_fields init_fields.f90)
add_library(preCalculations preCalculations.f90)
add_library(startFields startFields.f90)
add_library(readCheckPoints readCheckPoints.f90)
add_library(storeCheckPoints storeCheckPoints.f90)
add_library(integration integration.f90)
add_library(movie movie.f90)
add_library(communications communications.f90)
add_library(courant courant.f90)
add_library(get_nl get_nl.f90)
add_library(radialLoop radialLoop.f90)
add_library(rIteration rIteration.f90)
add_library(rIterThetaBlocking rIterThetaBlocking.f90)
add_library(rIterThetaBlocking_seq rIterThetaBlocking_seq.f90)
add_library(rIterThetaBlocking_OpenMP rIterThetaBlocking_OpenMP.f90)
add_library(nonlinear_bcs nonlinear_bcs.f90)
add_library(dtB dtB.f90)

# Outputs
add_library(nl_special_calc nl_special_calc.f90)
add_library(fields_average fields_average.f90)
add_library(output_data output_data.f90)
add_library(outOmega outOmega.f90)
add_library(outRot outRot.f90)
add_library(outMisc outMisc.f90)
add_library(out_coeff out_coeff.f90)
add_library(store_pot store_pot.f90)
add_library(out_movie_file out_movie_file.f90)
add_library(store_movie_IC store_movie_IC.f90)
add_library(out_graph_file out_graph_file.f90)
add_library(out_dtB_frame out_dtB_frame.f90)
add_library(out_RMS out_RMS.f90)
add_library(outPV3 outPV3.f90)
add_library(outPar outPar.f90)
add_library(out_TO out_TO.f90)
add_library(Egeos Egeos.f90)
add_library(kinetic_energy kinetic_energy.f90)
add_library(magnetic_energy magnetic_energy.f90)
add_library(getDlm getDlm.f90)
add_library(power power.f90)
add_library(spectra spectra.f90)
add_library(output output.f90)
add_library(radial_spectra radial_spectra.f90)

# Time advance
add_library(get_td get_td.f90)
add_library(step_time step_time.f90)
add_library(LMLoop LMLoop.f90)
add_library(updateS updateS.f90)
add_library(updateB updateB.f90)
add_library(updateZ updateZ.f90)
add_library(updateWP updateWP.f90)

add_library(cutils_iface cutils_iface.f90)
add_library(c_utils c_utils.c)


# Lapack
if(USE_MKL MATCHES yes)
   add_library(lapack algebra_mkl.f90)
else()
   add_library(lapack algebra.f90)
endif()

# FFTs
if(USE_FFTLIB MATCHES JW)
   add_library(fft fft.f90)
endif()
if(USE_FFTLIB MATCHES MKL)
   add_library(fft fft_mkl.f90)
endif()

# HDF5
if(USE_HDF5 MATCHES yes)
   add_library(hdf5Helpers hdf5Helpers.f90)
   target_link_libraries(hdf5Helpers blocking LMLoop_data precision_mod)
   target_link_libraries(readCheckPoints hdf5Helpers)
   target_link_libraries(storeCheckPoints hdf5Helpers LMLoop_data)
endif()

target_link_libraries(parallel mpimod)
target_link_libraries(precision_mod mpimod)
target_link_libraries(phys_param precision_mod)
target_link_libraries(num_param truncation precision_mod)
target_link_libraries(char_manip precision_mod)
target_link_libraries(Bext precision_mod)
target_link_libraries(debugging precision_mod)
target_link_libraries(constants precision_mod)
target_link_libraries(mat truncation precision_mod)
target_link_libraries(fields truncation radial_data LMLoop_data precision_mod)
target_link_libraries(dt_fieldsLast truncation LMLoop_data precision_mod)
target_link_libraries(output_data logic parallel char_manip precision_mod)
target_link_libraries(useful logic output_data parallel precision_mod)
target_link_libraries(blocking truncation LMmapping logic parallel output_data
                      useful precision_mod)
target_link_libraries(plms constants precision_mod)
target_link_libraries(RMS blocking truncation)
target_link_libraries(Grenoble truncation precision_mod)
target_link_libraries(timing parallel precision_mod)
target_link_libraries(radial_data truncation parallel logic)
target_link_libraries(plms constants precision_mod)
target_link_libraries(LMLoop_data parallel blocking logic)
target_link_libraries(radial truncation radial_data parallel precision_mod)
target_link_libraries(useful logic output_data parallel constants precision_mod)
target_link_libraries(horizontal truncation phys_param num_param precision_mod
                      radial logic blocking plms fft)
target_link_libraries(radial truncation radial_data mat lapack phys_param num_param 
                      precision_mod output_data init_costf cosine_transform chebyshev_polynoms
                      radial_derivatives)
target_link_libraries(init_costf useful constants precision_mod)
target_link_libraries(cosine_transform truncation fft_fac constants)
target_link_libraries(fft_fac constants precision_mod)
target_link_libraries(chebyshev_polynoms constants logic precision_mod)
target_link_libraries(radial_derivatives cosine_transform constants precision_mod)
target_link_libraries(radial_derivatives_even cosine_transform constants precision_mod)
target_link_libraries(TO_helpers truncation blocking horizontal precision_mod constants)
target_link_libraries(Namelists truncation constants phys_param num_param logic
                      init_fields output_data parallel TO Grenoble radial Bext char_manip
                      movie)
target_link_libraries(TO truncation blocking horizontal radial_data phys_param radial
                      logic legendre_grid_to_spec)
target_link_libraries(legendre_grid_to_spec truncation blocking horizontal)
target_link_libraries(legendre_spec_to_grid truncation blocking horizontal logic constants
                      legendre_helpers)
target_link_libraries(init_fields truncation blocking radial horizontal constants logic
                      fft legendre_grid_to_spec useful phys_param LMLoop_data
                      cosine_transform)
target_link_libraries(preCalculations truncation phys_param num_param constants radial
                      horizontal init_fields blocking logic output_data useful integration
                      parallel Bext)
target_link_libraries(integration cosine_transform precision_mod constants)
target_link_libraries(startFields truncation phys_param num_param constants radial
                      init_fields blocking logic dt_fieldsLast Grenoble readCheckPoints
                      fields radial_derivatives_even radial_derivatives useful communications)
target_link_libraries(movie truncation parallel radial_data output_data logic radial
                      horizontal char_manip useful)
target_link_libraries(radialLoop truncation blocking parallel phys_param logic output_data
                      constants legendre_spec_to_grid fft get_td out_graph_file radial_data 
                      rIteration rIterThetaBlocking rIterThetaBlocking_OpenMP
                      rIterThetaBlocking_seq)
target_link_libraries(rIteration precision_mod)
target_link_libraries(rIterThetaBlocking rIteration truncation blocking logic radial
                      radial_data fft legendre_spec_to_grid get_nl
                      legendre_helpers get_td nonlinear_bcs
                      legendre_grid_to_spec  out_graph_file)
target_link_libraries(rIterThetaBlocking_OpenMP rIterThetaBlocking TO dtB
                      out_movie_file outRot nl_special_calc courant)
target_link_libraries(rIterThetaBlocking_seq rIterThetaBlocking TO dtB
                      out_movie_file outRot nl_special_calc)
target_link_libraries(legendre_helpers truncation num_param radial horizontal
                      logic constants blocking Grenoble fields TO)
target_link_libraries(outRot truncation blocking num_param phys_param logic radial 
                      output_data constants integration Grenoble horizontal)
target_link_libraries(communications blocking parallel LMLoop_data truncation
                      radial_data)
target_link_libraries(get_nl truncation horizontal logic phys_param radial blocking)
target_link_libraries(nl_special_calc truncation blocking horizontal radial 
                      legendre_grid_to_spec radial_data phys_param constants)
target_link_libraries(courant truncation blocking num_param phys_param radial
                      horizontal logic useful radial_data parallel)
target_link_libraries(dtB truncation radial_data parallel LMLoop_data fft radial logic
                      blocking horizontal phys_param communications legendre_grid_to_spec
                      radial_spectra radial_derivatives)
target_link_libraries(radial_spectra truncation blocking horizontal radial num_param 
                      output_data logic useful char_manip radial_data LMmapping constants)
target_link_libraries(out_graph_file truncation horizontal phys_param parallel
                      num_param radial logic blocking output_data legendre_helpers fft
                      legendre_spec_to_grid)
target_link_libraries(outOmega truncation blocking logic output_data radial
                      cosine_transform plms)
target_link_libraries(chebInt chebyshev_polynoms init_costf radial_derivatives
                      constants)
target_link_libraries(readCheckPoints truncation blocking phys_param init_fields constants 
                      radial logic init_costf cosine_transform)
target_link_libraries(storeCheckPoints phys_param num_param logic output_data init_fields
                      dt_fieldsLast)
target_link_libraries(out_movie_file truncation output_data blocking radial horizontal 
                      movie fields num_param fft logic out_dtB_frame)
target_link_libraries(out_dtB_frame truncation blocking radial logic horizontal dtB 
                      radial_derivatives radial_derivatives_even fft constants)
target_link_libraries(out_coeff logic precision_mod)
target_link_libraries(outMisc truncation num_param phys_param logic output_data blocking
                      constants radial horizontal legendre_spec_to_grid Egeos integration)
target_link_libraries(power truncation blocking horizontal phys_param num_param radial 
                      logic output_data outRot integration useful)
target_link_libraries(nonlinear_bcs truncation horizontal fft legendre_grid_to_spec
                      blocking phys_param radial radial_data)
target_link_libraries(out_TO truncation horizontal num_param phys_param radial logic 
                      blocking constants output_data TO legendre_grid_to_spec
                      cosine_transform chebInt plms TO_helpers char_manip integration 
                      useful)
target_link_libraries(magnetic_energy truncation blocking horizontal phys_param num_param
                      constants radial logic output_data movie Bext integration useful)
target_link_libraries(outPV3 truncation horizontal cosine_transform phys_param radial
                      logic blocking output_data plms TO_helpers constants fft)
target_link_libraries(Egeos truncation phys_param num_param horizontal radial blocking
                      constants output_data logic cosine_transform chebInt plms fft
                      LMLoop_data communications)
target_link_libraries(output truncation blocking phys_param num_param logic output_data
                      radial horizontal constants fields radial_spectra out_movie_file
                      outRot kinetic_energy magnetic_energy spectra fields_average power
                      getDlm useful out_TO outMisc out_graph_file out_RMS char_manip
                      storeCheckPoints integration outPV3 outOmega parallel
                      communications out_coeff outPar dtB store_movie_IC RMS store_pot)
target_link_libraries(outPar truncation blocking phys_param logic radial output_data
                      radial_data horizontal fields legendre_spec_to_grid
                      num_param integration constants)
target_link_libraries(spectra truncation blocking num_param horizontal phys_param radial
                      logic output_data useful LMLoop_data integration)
target_link_libraries(out_RMS truncation parallel blocking num_param phys_param radial
                      constants logic output_data RMS radial_derivatives RMS_helpers
                      integration communications horizontal dtB)
target_link_libraries(RMS_helpers truncation blocking radial constants horizontal useful
                      integration LMmapping chebyshev_polynoms init_costf radial_derivatives)
target_link_libraries(kinetic_energy truncation blocking num_param radial horizontal
                      output_data logic constants phys_param useful integration
                      communications)
target_link_libraries(getDlm truncation blocking num_param constants horizontal radial 
                      constants LMLoop_data useful integration)
target_link_libraries(fields_average truncation output_data blocking logic constants 
                      radial horizontal radial_derivatives radial_derivatives_even
                      spectra kinetic_energy out_graph_file magnetic_energy
                      legendre_spec_to_grid fft out_coeff store_pot legendre_helpers)
target_link_libraries(store_movie_IC truncation blocking logic radial_data movie radial
                      horizontal legendre_helpers legendre_spec_to_grid fft
                      out_movie_file phys_param)
target_link_libraries(store_pot truncation radial output_data phys_param logic
                      cosine_transform LMLoop_data parallel communications)

target_link_libraries(updateB truncation blocking horizontal phys_param num_param
                      init_fields logic radial Bext constants mat RMS lapack
                      radial_derivatives cosine_transform radial_derivatives_even
                      RMS_helpers LMLoop_data)
target_link_libraries(updateZ truncation blocking phys_param num_param constants radial
                      horizontal logic init_fields mat RMS TO lapack cosine_transform
                      radial_derivatives outRot RMS_helpers LMLoop_data communications)
target_link_libraries(updateS truncation blocking horizontal phys_param radial logic 
                      output_data init_fields constants mat lapack cosine_transform
                      radial_derivatives LMLoop_data)
target_link_libraries(updateWP truncation blocking num_param phys_param radial 
                      horizontal logic mat RMS cosine_transform radial_derivatives
                      RMS_helpers LMLoop_data)

target_link_libraries(LMLoop truncation blocking parallel logic radial_data
                      output_data fields dt_fieldsLast mat updateS updateB
                      updateZ updateWP timing LMLoop_data communications debugging)
target_link_libraries(step_time truncation blocking phys_param num_param radial_data
                      logic output_data output movie char_manip timing courant
                      radialLoop nonlinear_bcs LMLoop dt_fieldsLast useful c_utils)
target_link_libraries(get_td truncation blocking horizontal phys_param num_param radial 
                      logic RMS RMS_helpers legendre_helpers fields cutils_iface)

target_link_libraries(lapack constants precision_mod)
target_link_libraries(fft truncation blocking useful constants)



# Main executable
add_executable(${EXEC} magic.f90 )
set_property(TARGET ${EXEC} PROPERTY LINKER_LANGUAGE Fortran)
target_link_libraries(${EXEC} truncation num_param parallel logic output_data timing 
                      constants precision_mod useful startFields preCalculations
                      step_time kinetic_energy magnetic_energy horizontal fields_average
                      Egeos spectra out_TO output_data outPV3 Namelists blocking)


# Link MPI
if(USE_MPI MATCHES yes)
   if(MPI_Fortran_FOUND)
      target_link_libraries(${EXEC} ${MPI_Fortran_LIBRARIES})
   endif()
endif()

# Link MKL
if(USE_MKL MATCHES yes AND USE_FFTLIB MATCHES MKL)
   add_library(lapack95 ${MKLROOT}/include/lapack.f90)
   target_link_libraries(lapack95 ${MKLCORE} ${MKLLAPACK} ${MKLSEQ} ${MKL64})
   target_link_libraries(lapack lapack95)
elseif(USE_MKL MATCHES yes AND USE_FFTLIB MATCHES JW)
   add_library(lapack95 ${MKLROOT}/include/lapack.f90)
   target_link_libraries(lapack95 ${MKLCORE} ${MKLLAPACK} ${MKLSEQ} ${MKL64})
   target_link_libraries(lapack lapack95)
elseif(USE_MKL MATCHES no AND USE_FFTLIB MATCHES MKL)
   target_link_libraries(${EXEC} ${MKLCORE} ${MKLLAPACK} ${MKLSEQ} ${MKL64})
endif()

# Link HDF5
if(USE_HDF5 MATCHES yes)
   if (HDF5_FOUND AND HDF5_IS_PARALLEL)
      target_link_libraries(${EXEC} ${HDF5_Fortran_LIBRARIES})
   endif()
endif()
