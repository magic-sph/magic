!$Id$
!********************************************************************
!  Module containing logicals that control the run
!********************************************************************
        
!   !------------ This is release 1 level 1  --------------!
!   !------------ Created on 1/17/02  by JW. --------------!
MODULE logic
  IMPLICIT NONE
  LOGICAL :: l_update_v,l_update_b,l_update_s
  LOGICAL :: l_mag,l_conv,l_mag_kin,l_SRIC,l_SRMA
  LOGICAL :: l_heat,l_heat_nl
  LOGICAL :: l_mag_nl,l_conv_nl,l_mag_LF,l_corr
  LOGICAL :: l_rot_ic,l_rot_ma
  LOGICAL :: l_z10mat,l_cond_ic,l_cond_ma
  LOGICAL :: l_time_hits 
  LOGICAL :: l_average
  LOGICAL :: l_movie,l_movie_oc,l_movie_ic
  LOGICAL :: l_save_out
  LOGICAL :: l_true_time
  LOGICAL :: l_cmb_field
  LOGICAL :: l_dt_cmb_field
  LOGICAL :: l_storeBpot
  LOGICAL :: l_storeVpot
  LOGICAL :: l_storeTpot
  LOGICAL :: l_storePot
  LOGICAL :: l_r_field
  LOGICAL :: l_b_nl_cmb,l_b_nl_icb
  LOGICAL :: l_correct_AMe,l_correct_AMz
  LOGICAL :: l_HTmovie,l_HT
  LOGICAL :: l_dtBmovie,l_dtB
  LOGICAL :: l_store_frame
  LOGICAL :: l_non_rot
  LOGICAL :: l_rMagSpec,l_DTrMagSpec
  LOGICAL :: l_TO,l_TOmovie 
  LOGICAL :: l_hel
  LOGICAL :: l_anel,l_isothermal,l_interior_model
  LOGICAL :: l_AM
  LOGICAL :: l_power
  LOGICAL :: l_drift
  LOGICAL :: l_iner
  LOGICAL :: l_runTimeLimit
  LOGICAL :: l_RMS,l_RMStest
  LOGICAL :: l_par
  LOGICAL :: l_corrMov
  LOGICAL :: l_PV
  LOGICAL :: l_newmap
  LOGICAL :: l_plotmap
  LOGICAL :: l_prms
  LOGICAL :: l_viscBcCalc 

  !COMMON/logic/l_update_v,l_update_b,l_update_s,                  &
  !     &               l_mag,l_conv,l_mag_nl,l_conv_nl,l_mag_kin,         &
  !     &               l_heat,l_heat_nl,l_SRIC,l_SRMA,l_corr,             &
  !     &               l_mag_LF,l_rot_ic,l_rot_ma,l_z10mat,               &
  !     &               l_cond_ic,l_cond_ma,l_time_hits,                   &
  !     &               l_average,l_save_out,l_non_rot,l_anel,             &
  !     &               l_interior_model,l_isothermal,                     &
  !     &               l_movie,l_movie_oc,l_movie_ic,l_HTmovie,           &
  !     &               l_dtBmovie,l_dtB,l_HT,l_true_time,l_cmb_field,     &
  !     &               l_dt_cmb_field,l_r_field,                          &
  !     &               l_b_nl_cmb,l_b_nl_icb,                             &
  !     &               l_correct_AMe,l_correct_AMz,l_store_frame,         &
  !     &               l_rMagSpec,l_DTrMagSpec,                           &
  !     &               l_TO,l_TOmovie,l_hel,l_AM,l_power,                 &
  !     &               l_drift,l_runTimeLimit,l_storeBpot,l_storeVpot,    &
  !     &               l_storeTpot,l_storePot,l_RMS,l_RMStest,            &
  !     &               l_par,l_corrMov,l_newmap,l_plotmap,l_prms,         &
  !     &               l_PV,l_iner,lviscBcCalc
END MODULE logic
