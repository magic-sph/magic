!$Id$
module logic
   !-------------------------------------------------------
   !  Module containing the logicals that control the run
   !-------------------------------------------------------

   implicit none
 
   logical :: l_update_v,l_update_b,l_update_s
   logical :: l_mag,l_conv,l_mag_kin,l_SRIC,l_SRMA
   logical :: l_heat,l_heat_nl
   logical :: l_mag_nl,l_conv_nl,l_mag_LF,l_corr
   logical :: l_rot_ic,l_rot_ma
   logical :: l_z10mat,l_cond_ic,l_cond_ma
   logical :: l_time_hits 
   logical :: l_average
   logical :: l_movie,l_movie_oc,l_movie_ic
   logical :: l_save_out
   logical :: l_true_time
   logical :: l_cmb_field
   logical :: l_dt_cmb_field
   logical :: l_storeBpot
   logical :: l_storeVpot
   logical :: l_storeTpot
   logical :: l_storePot
   logical :: l_r_field
   logical :: l_r_fieldT
   logical :: l_b_nl_cmb,l_b_nl_icb
   logical :: l_correct_AMe,l_correct_AMz
   logical :: l_HTmovie,l_HT
   logical :: l_dtBmovie,l_dtB
   logical :: l_store_frame
   logical :: l_non_rot
   logical :: l_rMagSpec,l_DTrMagSpec
   logical :: l_TO,l_TOmovie 
   logical :: l_hel
   logical :: l_anel,l_isothermal
   logical :: l_anelastic_liquid
   logical :: l_AM
   logical :: l_power
   logical :: l_drift
   logical :: l_iner
   logical :: l_runTimeLimit
   logical :: l_RMS,l_RMStest
   logical :: l_par
   logical :: l_corrMov
   logical :: l_PV
   logical :: l_newmap
   logical :: l_viscBcCalc 
   logical :: l_fluxProfs
   logical :: l_perpPar
   logical :: l_LCR
   logical :: lVerbose

end module logic
