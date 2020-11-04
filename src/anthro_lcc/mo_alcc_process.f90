!>
!! @brief Contains the routines for the alcc processes
!!
!! @author
!!  Rainer Schneck, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Rainer Schneck (2019-01-08)
!!
MODULE mo_alcc_process
#ifndef __NO_JSBACH__

   USE mo_kind,      ONLY: wp

   IMPLICIT NONE
   PRIVATE

   PUBLIC :: relocate_vars_to_new_pft_fractions

   CHARACTER(len=*), PARAMETER :: modname = 'mo_alcc_process'

 CONTAINS


  ! Calculations to relocate carbon to new pft fractions from alcc
  ! R: Note, this must be done for all vars that differ between pfts!
  !    We could also call it for each variable separately instead of naming each variable.
  !    We could also split this Subroutine for each process separately.
  !    As soon as the general scheme for C and N is implemented, this procedure
  !    has to be adapted.
  ELEMENTAL SUBROUTINE relocate_vars_to_new_pft_fractions( &
      & d_fract_grow_pfts,          &
      & f_old_to_new_grow_pfts,     &
      & f_d_to_new_grow_pfts,       &
      !
      ! CARBON
      & c_green_lost_avg,           &
      & c_woods_lost_avg,           &
      & c_reserve_lost_avg,         &
      & c_crop_harvest_lost_avg,    &
      ! size class 1
      & c_acid_ag1_lost_avg,        &
      & c_water_ag1_lost_avg,       &
      & c_ethanol_ag1_lost_avg,     &
      & c_nonsoluble_ag1_lost_avg,  &
      & c_acid_bg1_lost_avg,        &
      & c_water_bg1_lost_avg,       &
      & c_ethanol_bg1_lost_avg,     &
      & c_nonsoluble_bg1_lost_avg,  &
      & c_humus_1_lost_avg,         &
      ! size class 2
      & c_acid_ag2_lost_avg,        &
      & c_water_ag2_lost_avg,       &
      & c_ethanol_ag2_lost_avg,     &
      & c_nonsoluble_ag2_lost_avg,  &
      & c_acid_bg2_lost_avg,        &
      & c_water_bg2_lost_avg,       &
      & c_ethanol_bg2_lost_avg,     &
      & c_nonsoluble_bg2_lost_avg,  &
      & c_humus_2_lost_avg,         &
      & c_bg_sum_lost_avg,          &
      & c_ag_sum_1_lost_avg,        &
      & c_litter_green_ag_lost_avg, &
      & c_litter_green_bg_lost_avg, &
      & c_litter_wood_ag_lost_avg,  &
      & c_litter_wood_bg_lost_avg,  &
      ! rest
      & NPP_pot_sum_lost_avg,          &
      & GPP_sum_lost_avg,              &
      & alpha_lost_avg,                &
      & alpha_sum_lost_avg,            &
      & alpha_yDayMean_lost_avg,       &
      & soil_respiration_lost_avg,     &
      & NPP_flux_correction_lost_avg,  &
      & NPP_pot_yDayMean_lost_avg,     &
      & NPP_act_yDayMean_lost_avg,     &
      & excess_NPP_lost_avg,           &
      & GPP_yDayMean_lost_avg,         &
      & root_exudates_lost_avg,        &
      & cflux_c_green_2_herb_lost_avg,           &
      & cflux_herb_2_littergreen_lost_avg,       &
      & cflux_herb_2_atm_lost_avg,               &
      & cflux_dist_greenreserve_2_soil_lost_avg, &
      & cflux_dist_woods_2_soil_lost_avg,        &
      & cflux_fire_all_2_atm_lost_avg,           &
      & current_max_green_lost_avg,              &
      & veg_carbon_at_max_green_lost_avg,        &
      !
      ! ASSIMI
      & gross_assimilation_ca_lost_avg,    &
      & gross_assimilation_lost_avg,       &
      & dark_respiration_ca_lost_avg,      &
      & dark_respiration_lost_avg,         &
      & NPP_pot_rate_ca_lost_avg,          &
      & canopy_cond_lost_avg,    &
      & cover_fract_pot_lost_avg,          &
      ! for each cl:
      & assimi_scaling_fact_cl_lost_avg,   &
      & canopy_cond_cl_lost_avg,           &
      & gross_assimilation_cl_lost_avg,    &
      & dark_respiration_cl_lost_avg,      &
      & carbox_rate_max_cl_lost_avg,       &
      & e_transport_rate_max_cl_lost_avg,  &
      & carbox_rate_cl_lost_avg,           &
      & e_transport_rate_cl_lost_avg,      &
      & CO2_conc_leaf_cl_lost_avg,         &
      !
      ! DISTURB
      & cconservation_fire_lost_avg, &
      & cconservation_wind_lost_avg, &
      & burned_fract_lost_avg,       &
      & damaged_fract_lost_avg,      &
      !
      ! PHENO
      & fract_fpc_lost_avg,                    &
      !& fract_fpc_max_lost_avg,                &
      & previous_day_NPP_pot_rate_ca_lost_avg, &
      & cconservation_allom_lost_avg,          &
      & day_NPP_sum_lost_avg,                  &
      !
      ! RAD
      !& sw_srf_net_lost_avg,    &
      !& swvis_srf_net_lost_avg, &
      !& swnir_srf_net_lost_avg, &
      !& lw_srf_net_lost_avg,    &
      !& rad_srf_net_lost_avg,   &
      !& alb_vis_lost_avg,       &
      !& alb_nir_lost_avg,       &
      !& alb_vis_soil_lost_avg,  &
      !& alb_nir_soil_lost_avg,  &
      !& alb_vis_lnd_lost_avg,   &
      !& alb_nir_lnd_lost_avg,   &
      !& alb_vis_can_lost_avg,   &
      !& alb_nir_can_lost_avg,   &
      & faPAR_lost_avg,         &
      & soil_reflectivity_par_lost_avg, &
      ! & apar_nacc_lost_avg,     &
      ! for each cl:
      & faPAR_cl_lost_avg,             &
      & apar_per_lai_cl_lost_avg,      &
      ! & apar_per_lai_cl_nacc_lost_avg, &
      ! & apar_per_lai_cl_tavg_lost_avg, &
      !
      ! TURB
      fact_q_air_lost_avg, &
      fact_qsat_srf_lost_avg,       &
      fact_qsat_trans_srf_lost_avg, &
      !
      !
      ! CARBON
      & c_green_, c_woods_, c_reserve_, c_crop_harvest_,              &
      & c_acid_ag1_, c_water_ag1_, c_ethanol_ag1_, c_nonsoluble_ag1_, &
      & c_acid_bg1_, c_water_bg1_, c_ethanol_bg1_, c_nonsoluble_bg1_, &
      & c_acid_ag2_, c_water_ag2_, c_ethanol_ag2_, c_nonsoluble_ag2_, &
      & c_acid_bg2_, c_water_bg2_, c_ethanol_bg2_, c_nonsoluble_bg2_, &
      & c_humus_1_, c_humus_2_,                                       &
      & c_bg_sum_, c_ag_sum_1_,                                       &
      & c_litter_green_ag_, c_litter_green_bg_,                       &
      & c_litter_wood_ag_, c_litter_wood_bg_,                         &
      ! rest
      & NPP_pot_sum_,          &
      & GPP_sum_,              &
      & alpha_,                &
      & alpha_sum_,            &
      & alpha_yDayMean_,       &
      & soil_respiration_,     &
      & NPP_flux_correction_,  &
      & NPP_pot_yDayMean_,     &
      & NPP_act_yDayMean_,     &
      & excess_NPP_,           &
      & GPP_yDayMean_,         &
      & root_exudates_,        &
      & cflux_c_green_2_herb_,           &
      & cflux_herb_2_littergreen_,       &
      & cflux_herb_2_atm_,               &
      & cflux_dist_greenreserve_2_soil_, &
      & cflux_dist_woods_2_soil_,        &
      & cflux_fire_all_2_atm_,           &
      & current_max_green_,              &
      & veg_carbon_at_max_green_,        &
      !
      !ASSIMI
      & gross_assimilation_ca,    &
      & gross_assimilation,       &
      & dark_respiration_ca,      &
      & dark_respiration,         &
      & NPP_pot_rate_ca,          &
      & canopy_cond,    &
      & cover_fract_pot,          &
      ! for each cl:
      & assimi_scaling_fact_cl,   &
      & canopy_cond_cl,           &
      & gross_assimilation_cl,    &
      & dark_respiration_cl,      &
      & carbox_rate_max_cl,       &
      & e_transport_rate_max_cl,  &
      & carbox_rate_cl,           &
      & e_transport_rate_cl,      &
      & CO2_conc_leaf_cl,         &
      !
      ! DISTURB
      & cconservation_fire, &
      & cconservation_wind, &
      & burned_fract,       &
      & damaged_fract,      &
      !
      ! PHENO
      & fract_fpc,                    &
      !& fract_fpc_max,                &
      & previous_day_NPP_pot_rate_ca, &
      & cconservation_allom,          &
      & day_NPP_sum,                  &
      !
      ! RAD
      !& sw_srf_net,    &
      !& swvis_srf_net, &
      !& swnir_srf_net, &
      !& lw_srf_net,    &
      !& rad_srf_net,   &
      !& alb_vis,       &
      !& alb_nir,       &
      !& alb_vis_soil,  &
      !& alb_nir_soil , &
      !& alb_vis_lnd,   &
      !& alb_nir_lnd,   &
      !& alb_vis_can,   &
      !& alb_nir_can,   &
      & faPAR,         &
      & soil_reflectivity_par, &
      ! & apar_nacc,     &
      ! for each cl:
      & faPAR_cl,             &
      & apar_per_lai_cl,      &
      ! & apar_per_lai_cl_nacc, &
      ! & apar_per_lai_cl_tavg, &
      !
      ! TURB
      fact_q_air, &
      fact_qsat_srf,       &
      fact_qsat_trans_srf  &
      & )



    ! Arguments
    REAL(wp), INTENT(in) ::     &
      & d_fract_grow_pfts,      &
      & f_old_to_new_grow_pfts, &
      & f_d_to_new_grow_pfts
      !

    REAL(wp), INTENT(in), OPTIONAL :: &
      ! CARBON
      & c_green_lost_avg,           &
      & c_woods_lost_avg,           &
      & c_reserve_lost_avg,         &
      & c_crop_harvest_lost_avg,    &
      ! size class 1
      & c_acid_ag1_lost_avg,        &
      & c_water_ag1_lost_avg,       &
      & c_ethanol_ag1_lost_avg,     &
      & c_nonsoluble_ag1_lost_avg,  &
      & c_acid_bg1_lost_avg,        &
      & c_water_bg1_lost_avg,       &
      & c_ethanol_bg1_lost_avg,     &
      & c_nonsoluble_bg1_lost_avg,  &
      & c_humus_1_lost_avg,         &
      ! size class 2
      & c_acid_ag2_lost_avg,        &
      & c_water_ag2_lost_avg,       &
      & c_ethanol_ag2_lost_avg,     &
      & c_nonsoluble_ag2_lost_avg,  &
      & c_acid_bg2_lost_avg,        &
      & c_water_bg2_lost_avg,       &
      & c_ethanol_bg2_lost_avg,     &
      & c_nonsoluble_bg2_lost_avg,  &
      & c_humus_2_lost_avg,         &
      & c_bg_sum_lost_avg,          &
      & c_ag_sum_1_lost_avg,        &
      & c_litter_green_ag_lost_avg, &
      & c_litter_green_bg_lost_avg, &
      & c_litter_wood_ag_lost_avg,  &
      & c_litter_wood_bg_lost_avg,  &
      ! rest
      & NPP_pot_sum_lost_avg,          &
      & GPP_sum_lost_avg,              &
      & alpha_lost_avg,                &
      & alpha_sum_lost_avg,            &
      & alpha_yDayMean_lost_avg,       &
      & soil_respiration_lost_avg,     &
      & NPP_flux_correction_lost_avg,  &
      & NPP_pot_yDayMean_lost_avg,     &
      & NPP_act_yDayMean_lost_avg,     &
      & excess_NPP_lost_avg,           &
      & GPP_yDayMean_lost_avg,         &
      & root_exudates_lost_avg,        &
      & cflux_c_green_2_herb_lost_avg,           &
      & cflux_herb_2_littergreen_lost_avg,       &
      & cflux_herb_2_atm_lost_avg,               &
      & cflux_dist_greenreserve_2_soil_lost_avg, &
      & cflux_dist_woods_2_soil_lost_avg,        &
      & cflux_fire_all_2_atm_lost_avg,           &
      & current_max_green_lost_avg,              &
      & veg_carbon_at_max_green_lost_avg,        &
      !
      ! ASSIMI
      & gross_assimilation_ca_lost_avg,    &
      & gross_assimilation_lost_avg,       &
      & dark_respiration_ca_lost_avg,      &
      & dark_respiration_lost_avg,         &
      & NPP_pot_rate_ca_lost_avg,          &
      & canopy_cond_lost_avg,    &
      & cover_fract_pot_lost_avg,          &
      ! for each cl:
      & assimi_scaling_fact_cl_lost_avg,   &
      & canopy_cond_cl_lost_avg,           &
      & gross_assimilation_cl_lost_avg,    &
      & dark_respiration_cl_lost_avg,      &
      & carbox_rate_max_cl_lost_avg,       &
      & e_transport_rate_max_cl_lost_avg,  &
      & carbox_rate_cl_lost_avg,           &
      & e_transport_rate_cl_lost_avg,      &
      & CO2_conc_leaf_cl_lost_avg,         &
      !
      ! DISTURB
      & cconservation_fire_lost_avg, &
      & cconservation_wind_lost_avg, &
      & burned_fract_lost_avg,       &
      & damaged_fract_lost_avg,      &
      !
      ! PHENO
      & fract_fpc_lost_avg,                    &
      !& fract_fpc_max_lost_avg,                &
      & previous_day_NPP_pot_rate_ca_lost_avg, &
      & cconservation_allom_lost_avg,          &
      & day_NPP_sum_lost_avg,                  &
      !
      ! RAD
      !& sw_srf_net_lost_avg,    &
      !& swvis_srf_net_lost_avg, &
      !& swnir_srf_net_lost_avg, &
      !& lw_srf_net_lost_avg,    &
      !& rad_srf_net_lost_avg,   &
      !& alb_vis_lost_avg,       &
      !& alb_nir_lost_avg,       &
      !& alb_vis_soil_lost_avg,  &
      !& alb_nir_soil_lost_avg,  &
      !& alb_vis_lnd_lost_avg,   &
      !& alb_nir_lnd_lost_avg,   &
      !& alb_vis_can_lost_avg,   &
      !& alb_nir_can_lost_avg,   &
      & faPAR_lost_avg,         &
      & soil_reflectivity_par_lost_avg, &
      ! & apar_nacc_lost_avg,     &
      ! for each cl:
      & faPAR_cl_lost_avg,             &
      & apar_per_lai_cl_lost_avg,      &
      ! & apar_per_lai_cl_nacc_lost_avg, &
      ! & apar_per_lai_cl_tavg_lost_avg, &
      !
      ! TURB
      fact_q_air_lost_avg, &
      fact_qsat_srf_lost_avg,       &
      fact_qsat_trans_srf_lost_avg

    REAL(wp), INTENT(inout), OPTIONAL  :: &
      ! Note: _ is to indicate that these vars represent _ca or_ta vars
      & c_green_,          &
      & c_woods_,          &
      & c_reserve_,        &
      & c_crop_harvest_,   &
      ! size class 1
      & c_acid_ag1_,       &
      & c_water_ag1_,      &
      & c_ethanol_ag1_,    &
      & c_nonsoluble_ag1_, &
      & c_acid_bg1_,       &
      & c_water_bg1_,      &
      & c_ethanol_bg1_,    &
      & c_nonsoluble_bg1_, &
      ! size class 2
      & c_acid_ag2_,       &
      & c_water_ag2_,      &
      & c_ethanol_ag2_,    &
      & c_nonsoluble_ag2_, &
      & c_acid_bg2_,       &
      & c_water_bg2_,      &
      & c_ethanol_bg2_,    &
      & c_nonsoluble_bg2_, &
      !
      & c_humus_1_,         &
      & c_humus_2_,         &
      & c_bg_sum_,          &
      & c_ag_sum_1_,        &
      & c_litter_green_ag_, &
      & c_litter_green_bg_, &
      & c_litter_wood_ag_,  &
      & c_litter_wood_bg_,  &
      !
      & NPP_pot_sum_,          &
      & GPP_sum_,              &
      & alpha_,                &
      & alpha_sum_,            &
      & alpha_yDayMean_,       &
      & soil_respiration_,     &
      & NPP_flux_correction_,  &
      & NPP_pot_yDayMean_,     &
      & NPP_act_yDayMean_,     &
      & excess_NPP_,           &
      & GPP_yDayMean_,         &
      & root_exudates_,        &
      & cflux_c_green_2_herb_,           &
      & cflux_herb_2_littergreen_,       &
      & cflux_herb_2_atm_,               &
      & cflux_dist_greenreserve_2_soil_, &
      & cflux_dist_woods_2_soil_,        &
      & cflux_fire_all_2_atm_,           &
      & current_max_green_,              &
      & veg_carbon_at_max_green_,        &
      !
      !ASSIMI
      & gross_assimilation_ca,    &
      & gross_assimilation,       &
      & dark_respiration_ca,      &
      & dark_respiration,         &
      & NPP_pot_rate_ca,          &
      & canopy_cond,    &
      & cover_fract_pot,          &
      ! for each cl:
      & assimi_scaling_fact_cl,   &
      & canopy_cond_cl,           &
      & gross_assimilation_cl,    &
      & dark_respiration_cl,      &
      & carbox_rate_max_cl,       &
      & e_transport_rate_max_cl,  &
      & carbox_rate_cl,           &
      & e_transport_rate_cl,      &
      & CO2_conc_leaf_cl,         &
      !
      ! DISTURB
      & cconservation_fire, &
      & cconservation_wind, &
      & burned_fract,       &
      & damaged_fract,      &
      !
      ! PHENO
      & fract_fpc,                    &
      !& fract_fpc_max,                &
      & previous_day_NPP_pot_rate_ca, &
      & cconservation_allom,          &
      & day_NPP_sum,                  &
      !
      ! RAD
      !& sw_srf_net,    &
      !& swvis_srf_net, &
      !& swnir_srf_net, &
      !& lw_srf_net,    &
      !& rad_srf_net,   &
      !& alb_vis,       &
      !& alb_nir,       &
      !& alb_vis_soil,  &
      !& alb_nir_soil , &
      !& alb_vis_lnd,   &
      !& alb_nir_lnd,   &
      !& alb_vis_can,   &
      !& alb_nir_can,   &
      & faPAR,         &
      & soil_reflectivity_par, &
      ! & apar_nacc,     &
      ! for each cl:
      & faPAR_cl,             &
      & apar_per_lai_cl,      &
      ! & apar_per_lai_cl_nacc, &
      ! & apar_per_lai_cl_tavg, &
      !
      ! TURB
      fact_q_air, &
      fact_qsat_srf,       &
      fact_qsat_trans_srf


    IF (d_fract_grow_pfts > 0._wp) THEN ! Change only for growing pfts existing values of variable
      ! CARBON
      IF (PRESENT(c_green_))          c_green_   = c_green_   * f_old_to_new_grow_pfts   +  c_green_lost_avg    * f_d_to_new_grow_pfts
      IF (PRESENT(c_woods_))          c_woods_   = c_woods_   * f_old_to_new_grow_pfts   +  c_woods_lost_avg    * f_d_to_new_grow_pfts
      IF (PRESENT(c_reserve_))        c_reserve_ = c_reserve_ * f_old_to_new_grow_pfts   +  c_reserve_lost_avg  * f_d_to_new_grow_pfts

      IF (PRESENT(c_crop_harvest_))   c_crop_harvest_ = c_crop_harvest_    * f_old_to_new_grow_pfts   &
                                 & +  c_crop_harvest_lost_avg    * f_d_to_new_grow_pfts
      IF (PRESENT(c_acid_ag1_))       c_acid_ag1_ = c_acid_ag1_            * f_old_to_new_grow_pfts   &
                                 & +  c_acid_ag1_lost_avg        * f_d_to_new_grow_pfts
      IF (PRESENT(c_water_ag1_))      c_water_ag1_ = c_water_ag1_          * f_old_to_new_grow_pfts   &
                                 & +  c_water_ag1_lost_avg       * f_d_to_new_grow_pfts
      IF (PRESENT(c_ethanol_ag1_))    c_ethanol_ag1_ = c_ethanol_ag1_      * f_old_to_new_grow_pfts   &
                                 & +  c_ethanol_ag1_lost_avg     * f_d_to_new_grow_pfts
      IF (PRESENT(c_nonsoluble_ag1_)) c_nonsoluble_ag1_ = c_nonsoluble_ag1_ * f_old_to_new_grow_pfts  &
                                 & +  c_nonsoluble_ag1_lost_avg  * f_d_to_new_grow_pfts
      IF (PRESENT(c_acid_bg1_))       c_acid_bg1_ = c_acid_bg1_            * f_old_to_new_grow_pfts   &
                                 & +  c_acid_bg1_lost_avg        * f_d_to_new_grow_pfts
      IF (PRESENT(c_water_bg1_))      c_water_bg1_ = c_water_bg1_          * f_old_to_new_grow_pfts   &
                                 & +  c_water_bg1_lost_avg       * f_d_to_new_grow_pfts
      IF (PRESENT(c_ethanol_bg1_))    c_ethanol_bg1_ = c_ethanol_bg1_      * f_old_to_new_grow_pfts   &
                                 & +  c_ethanol_bg1_lost_avg     * f_d_to_new_grow_pfts
      IF (PRESENT(c_nonsoluble_bg1_)) c_nonsoluble_bg1_ = c_nonsoluble_bg1_ * f_old_to_new_grow_pfts  &
                                 & +  c_nonsoluble_bg1_lost_avg  * f_d_to_new_grow_pfts
      IF (PRESENT(c_acid_ag2_))       c_acid_ag2_ = c_acid_ag2_            * f_old_to_new_grow_pfts   &
                                 & +  c_acid_ag2_lost_avg        * f_d_to_new_grow_pfts
      IF (PRESENT(c_water_ag2_))      c_water_ag2_ = c_water_ag2_          * f_old_to_new_grow_pfts   &
                                 & +  c_water_ag2_lost_avg       * f_d_to_new_grow_pfts
      IF (PRESENT(c_ethanol_ag2_))    c_ethanol_ag2_ = c_ethanol_ag2_      * f_old_to_new_grow_pfts   &
                                 & +  c_ethanol_ag2_lost_avg     * f_d_to_new_grow_pfts
      IF (PRESENT(c_nonsoluble_ag2_)) c_nonsoluble_ag2_ = c_nonsoluble_ag2_ * f_old_to_new_grow_pfts  &
                                 & +  c_nonsoluble_ag2_lost_avg  * f_d_to_new_grow_pfts
      IF (PRESENT(c_acid_bg2_))       c_acid_bg2_ = c_acid_bg2_             * f_old_to_new_grow_pfts  &
                                 & +  c_acid_bg2_lost_avg        * f_d_to_new_grow_pfts
      IF (PRESENT(c_water_bg2_))      c_water_bg2_ = c_water_bg2_           * f_old_to_new_grow_pfts  &
                                 & +  c_water_bg2_lost_avg       * f_d_to_new_grow_pfts
      IF (PRESENT(c_ethanol_bg2_))    c_ethanol_bg2_ = c_ethanol_bg2_       * f_old_to_new_grow_pfts  &
                                 & +  c_ethanol_bg2_lost_avg     * f_d_to_new_grow_pfts
      IF (PRESENT(c_nonsoluble_bg2_)) c_nonsoluble_bg2_ = c_nonsoluble_bg2_ * f_old_to_new_grow_pfts  &
                                 & +  c_nonsoluble_bg2_lost_avg  * f_d_to_new_grow_pfts
      IF (PRESENT(c_humus_1_))        c_humus_1_ = c_humus_1_               * f_old_to_new_grow_pfts  &
                                 & +  c_humus_1_lost_avg         * f_d_to_new_grow_pfts
      IF (PRESENT(c_humus_2_))        c_humus_2_ = c_humus_2_               * f_old_to_new_grow_pfts  &
                                 & +  c_humus_2_lost_avg         * f_d_to_new_grow_pfts
      IF (PRESENT(c_bg_sum_))         c_bg_sum_ = c_bg_sum_                 * f_old_to_new_grow_pfts  &
                                 & +  c_bg_sum_lost_avg          * f_d_to_new_grow_pfts
      IF (PRESENT(c_ag_sum_1_))       c_ag_sum_1_ = c_ag_sum_1_             * f_old_to_new_grow_pfts  &
                                 & +  c_ag_sum_1_lost_avg        * f_d_to_new_grow_pfts
      IF (PRESENT(c_litter_green_ag_)) c_litter_green_ag_ = c_litter_green_ag_ * f_old_to_new_grow_pfts &
                                 & +  c_litter_green_ag_lost_avg * f_d_to_new_grow_pfts
      IF (PRESENT(c_litter_green_bg_)) c_litter_green_bg_ = c_litter_green_bg_ * f_old_to_new_grow_pfts &
                                 & +  c_litter_green_bg_lost_avg * f_d_to_new_grow_pfts
      IF (PRESENT(c_litter_wood_ag_)) c_litter_wood_ag_ = c_litter_wood_ag_ * f_old_to_new_grow_pfts  &
                                 & +  c_litter_wood_ag_lost_avg  * f_d_to_new_grow_pfts
      IF (PRESENT(c_litter_wood_bg_)) c_litter_wood_bg_ = c_litter_wood_bg_ * f_old_to_new_grow_pfts  &
                                 & +  c_litter_wood_bg_lost_avg  * f_d_to_new_grow_pfts
      IF (PRESENT(NPP_pot_sum_))     NPP_pot_sum_ = NPP_pot_sum_             * f_old_to_new_grow_pfts &
                                 & +  NPP_pot_sum_lost_avg       * f_d_to_new_grow_pfts
      IF (PRESENT(GPP_sum_))         GPP_sum_ = GPP_sum_                     * f_old_to_new_grow_pfts &
                                 & +  GPP_sum_lost_avg           * f_d_to_new_grow_pfts
      IF (PRESENT(alpha_))            alpha_ = alpha_                        * f_old_to_new_grow_pfts &
                                 & +  alpha_lost_avg             * f_d_to_new_grow_pfts
      IF (PRESENT(alpha_sum_))        alpha_sum_ = alpha_sum_                * f_old_to_new_grow_pfts &
                                 & +  alpha_sum_lost_avg         * f_d_to_new_grow_pfts
      IF (PRESENT(alpha_yDayMean_)) alpha_yDayMean_ = alpha_yDayMean_        * f_old_to_new_grow_pfts &
                                 & +  alpha_yDayMean_lost_avg    * f_d_to_new_grow_pfts
      IF (PRESENT(soil_respiration_)) soil_respiration_ = soil_respiration_  * f_old_to_new_grow_pfts &
                                 & +  soil_respiration_lost_avg  * f_d_to_new_grow_pfts
      IF (PRESENT(NPP_flux_correction_)) NPP_flux_correction_ = NPP_flux_correction_ * f_old_to_new_grow_pfts &
                                 & +  NPP_flux_correction_lost_avg * f_d_to_new_grow_pfts
      IF (PRESENT(NPP_pot_yDayMean_)) NPP_pot_yDayMean_ = NPP_pot_yDayMean_  * f_old_to_new_grow_pfts &
                                 & +  NPP_pot_yDayMean_lost_avg  * f_d_to_new_grow_pfts
      IF (PRESENT(NPP_act_yDayMean_)) NPP_act_yDayMean_ = NPP_act_yDayMean_  * f_old_to_new_grow_pfts &
                                 & +  NPP_act_yDayMean_lost_avg  * f_d_to_new_grow_pfts
      IF (PRESENT(excess_NPP_))       excess_NPP_ = excess_NPP_              * f_old_to_new_grow_pfts &
                                 & +  excess_NPP_lost_avg        * f_d_to_new_grow_pfts
      IF (PRESENT(GPP_yDayMean_))     GPP_yDayMean_ = GPP_yDayMean_          * f_old_to_new_grow_pfts &
                                 & +  GPP_yDayMean_lost_avg      * f_d_to_new_grow_pfts
      IF (PRESENT(root_exudates_)) root_exudates_ = root_exudates_           * f_old_to_new_grow_pfts &
                                 & +  root_exudates_lost_avg     * f_d_to_new_grow_pfts
      IF (PRESENT(cflux_c_green_2_herb_)) cflux_c_green_2_herb_ = cflux_c_green_2_herb_ * f_old_to_new_grow_pfts             &
                                 & +  cflux_c_green_2_herb_lost_avg * f_d_to_new_grow_pfts
      IF (PRESENT(cflux_herb_2_littergreen_)) cflux_herb_2_littergreen_ = cflux_herb_2_littergreen_ * f_old_to_new_grow_pfts &
                                 & +  cflux_herb_2_littergreen_lost_avg * f_d_to_new_grow_pfts
      IF (PRESENT(cflux_herb_2_atm_)) cflux_herb_2_atm_ = cflux_herb_2_atm_  * f_old_to_new_grow_pfts                        &
                                 & +  cflux_herb_2_atm_lost_avg  * f_d_to_new_grow_pfts
      IF (PRESENT(cflux_dist_greenreserve_2_soil_)) cflux_dist_greenreserve_2_soil_ = cflux_dist_greenreserve_2_soil_ * f_old_to_new_grow_pfts &
                                 & +  cflux_dist_greenreserve_2_soil_lost_avg * f_d_to_new_grow_pfts
      IF (PRESENT(cflux_dist_woods_2_soil_)) cflux_dist_woods_2_soil_ = cflux_dist_woods_2_soil_ * f_old_to_new_grow_pfts    &
                                 & +  cflux_dist_woods_2_soil_lost_avg * f_d_to_new_grow_pfts
      IF (PRESENT(cflux_fire_all_2_atm_)) cflux_fire_all_2_atm_ = cflux_fire_all_2_atm_ * f_old_to_new_grow_pfts             &
                                 & +  cflux_fire_all_2_atm_lost_avg * f_d_to_new_grow_pfts
      IF (PRESENT(current_max_green_)) current_max_green_ = current_max_green_ * f_old_to_new_grow_pfts                      &
                                 & +  current_max_green_lost_avg  * f_d_to_new_grow_pfts
      IF (PRESENT(veg_carbon_at_max_green_)) veg_carbon_at_max_green_ = veg_carbon_at_max_green_ * f_old_to_new_grow_pfts    &
                                 & +  veg_carbon_at_max_green_lost_avg * f_d_to_new_grow_pfts

      ! ASSIMI
      IF (PRESENT(gross_assimilation_ca))   gross_assimilation_ca = gross_assimilation_ca    * f_old_to_new_grow_pfts   &
                                 & +  gross_assimilation_ca_lost_avg    * f_d_to_new_grow_pfts
      IF (PRESENT(gross_assimilation))   gross_assimilation = gross_assimilation    * f_old_to_new_grow_pfts   &
                                 & +  gross_assimilation_lost_avg    * f_d_to_new_grow_pfts
      IF (PRESENT(dark_respiration_ca))   dark_respiration_ca = dark_respiration_ca    * f_old_to_new_grow_pfts   &
                                 & +  dark_respiration_ca_lost_avg    * f_d_to_new_grow_pfts
      IF (PRESENT(dark_respiration))   dark_respiration = dark_respiration    * f_old_to_new_grow_pfts   &
                                 & +  dark_respiration_lost_avg    * f_d_to_new_grow_pfts
      IF (PRESENT(NPP_pot_rate_ca))   NPP_pot_rate_ca = NPP_pot_rate_ca    * f_old_to_new_grow_pfts   &
                                 & +  NPP_pot_rate_ca_lost_avg    * f_d_to_new_grow_pfts
      IF (PRESENT(canopy_cond))   canopy_cond = canopy_cond    * f_old_to_new_grow_pfts   &
                                 & +  canopy_cond_lost_avg    * f_d_to_new_grow_pfts
      IF (PRESENT(cover_fract_pot))   cover_fract_pot = cover_fract_pot    * f_old_to_new_grow_pfts   &
                                 & +  cover_fract_pot_lost_avg    * f_d_to_new_grow_pfts
      IF (PRESENT(assimi_scaling_fact_cl))   assimi_scaling_fact_cl =assimi_scaling_fact_cl    * f_old_to_new_grow_pfts   &
                                 & + assimi_scaling_fact_cl_lost_avg    * f_d_to_new_grow_pfts
      IF (PRESENT(canopy_cond_cl))   canopy_cond_cl = canopy_cond_cl    * f_old_to_new_grow_pfts   &
                                 & +  canopy_cond_cl_lost_avg    * f_d_to_new_grow_pfts
      IF (PRESENT(gross_assimilation_cl))   gross_assimilation_cl = gross_assimilation_cl    * f_old_to_new_grow_pfts   &
                                 & +  gross_assimilation_cl_lost_avg    * f_d_to_new_grow_pfts
      IF (PRESENT(dark_respiration_cl))   dark_respiration_cl = dark_respiration_cl    * f_old_to_new_grow_pfts   &
                                 & +  dark_respiration_cl_lost_avg    * f_d_to_new_grow_pfts
      IF (PRESENT(carbox_rate_max_cl))   carbox_rate_max_cl = carbox_rate_max_cl    * f_old_to_new_grow_pfts   &
                                 & +  carbox_rate_max_cl_lost_avg    * f_d_to_new_grow_pfts
      IF (PRESENT(e_transport_rate_max_cl))   e_transport_rate_max_cl = e_transport_rate_max_cl    * f_old_to_new_grow_pfts   &
                                 & +  e_transport_rate_max_cl_lost_avg    * f_d_to_new_grow_pfts
      IF (PRESENT(carbox_rate_cl))   carbox_rate_cl = carbox_rate_cl    * f_old_to_new_grow_pfts   &
                                 & +  carbox_rate_cl_lost_avg    * f_d_to_new_grow_pfts
      IF (PRESENT(e_transport_rate_cl))   e_transport_rate_cl = e_transport_rate_cl    * f_old_to_new_grow_pfts   &
                                 & +  e_transport_rate_cl_lost_avg    * f_d_to_new_grow_pfts
      IF (PRESENT(CO2_conc_leaf_cl))   CO2_conc_leaf_cl = CO2_conc_leaf_cl    * f_old_to_new_grow_pfts   &
                                 & +  CO2_conc_leaf_cl_lost_avg    * f_d_to_new_grow_pfts

      ! DISTURB
      IF (PRESENT(cconservation_fire))   cconservation_fire = cconservation_fire    * f_old_to_new_grow_pfts   &
                                 & +  cconservation_fire_lost_avg    * f_d_to_new_grow_pfts
      IF (PRESENT(cconservation_wind))   cconservation_wind = cconservation_wind    * f_old_to_new_grow_pfts   &
                                 & +  cconservation_wind_lost_avg    * f_d_to_new_grow_pfts
      IF (PRESENT(burned_fract))   burned_fract = burned_fract    * f_old_to_new_grow_pfts   &
                                 & +  burned_fract_lost_avg    * f_d_to_new_grow_pfts
      IF (PRESENT(damaged_fract))   damaged_fract = damaged_fract    * f_old_to_new_grow_pfts   &
                                 & +  damaged_fract_lost_avg    * f_d_to_new_grow_pfts

      ! PHENO
      IF (PRESENT(fract_fpc))   fract_fpc =  fract_fpc   * f_old_to_new_grow_pfts   &
                                      & +  fract_fpc_lost_avg    * f_d_to_new_grow_pfts
      !IF (PRESENT(fract_fpc_max))  fract_fpc_max  = fract_fpc_max    * f_old_to_new_grow_pfts   &
      !                                & +  fract_fpc_max_lost_avg    * f_d_to_new_grow_pfts
      IF (PRESENT(previous_day_NPP_pot_rate_ca))   previous_day_NPP_pot_rate_ca =  previous_day_NPP_pot_rate_ca   * f_old_to_new_grow_pfts   &
                                      & +  previous_day_NPP_pot_rate_ca_lost_avg    * f_d_to_new_grow_pfts
      IF (PRESENT(cconservation_allom))   cconservation_allom = cconservation_allom    * f_old_to_new_grow_pfts   &
                                      & +  cconservation_allom_lost_avg    * f_d_to_new_grow_pfts
      IF (PRESENT(day_NPP_sum))   day_NPP_sum = day_NPP_sum    * f_old_to_new_grow_pfts   &
                                      & +  day_NPP_sum_lost_avg    * f_d_to_new_grow_pfts

      ! RAD
!       IF (PRESENT(sw_srf_net))  sw_srf_net  = sw_srf_net    * f_old_to_new_grow_pfts   &
!                                       & +  sw_srf_net_lost_avg    * f_d_to_new_grow_pfts
!       IF (PRESENT(swvis_srf_net))   swvis_srf_net = swvis_srf_net    * f_old_to_new_grow_pfts   &
!                                       & +  swvis_srf_net_lost_avg    * f_d_to_new_grow_pfts
!       IF (PRESENT(swnir_srf_net))   swnir_srf_net = swnir_srf_net    * f_old_to_new_grow_pfts   &
!                                       & +  swnir_srf_net_lost_avg    * f_d_to_new_grow_pfts
!       IF (PRESENT(lw_srf_net))   lw_srf_net = lw_srf_net    * f_old_to_new_grow_pfts   &
!                                       & +  lw_srf_net_lost_avg    * f_d_to_new_grow_pfts
!       IF (PRESENT(rad_srf_net))   rad_srf_net = rad_srf_net    * f_old_to_new_grow_pfts   &
!                                       & +  rad_srf_net_lost_avg    * f_d_to_new_grow_pfts
!       IF (PRESENT(alb_vis))   alb_vis = alb_vis    * f_old_to_new_grow_pfts   &
!                                       & +  alb_vis_lost_avg    * f_d_to_new_grow_pfts
!       IF (PRESENT(alb_nir))   alb_nir = alb_nir    * f_old_to_new_grow_pfts   &
!                                       & +  alb_nir_lost_avg    * f_d_to_new_grow_pfts
!       IF (PRESENT(alb_vis_soil))   alb_vis_soil = alb_vis_soil    * f_old_to_new_grow_pfts   &
!                                       & +  alb_vis_soil_lost_avg    * f_d_to_new_grow_pfts
!       IF (PRESENT(alb_nir_soil))   alb_nir_soil = alb_nir_soil    * f_old_to_new_grow_pfts   &
!                                       & + alb_nir_soil_lost_avg    * f_d_to_new_grow_pfts
!       IF (PRESENT(alb_vis_lnd))   alb_vis_lnd = alb_vis_lnd     * f_old_to_new_grow_pfts   &
!                                       & + alb_vis_lnd_lost_avg    * f_d_to_new_grow_pfts
!       IF (PRESENT(alb_nir_lnd))   alb_nir_lnd = alb_nir_lnd    * f_old_to_new_grow_pfts   &
!                                       & +  alb_nir_lnd_lost_avg    * f_d_to_new_grow_pfts
!       IF (PRESENT(alb_vis_can))  alb_vis_can  = alb_vis_can     * f_old_to_new_grow_pfts   &
!                                       & +  alb_vis_can_lost_avg    * f_d_to_new_grow_pfts
!       IF (PRESENT(alb_nir_can))   alb_nir_can = alb_nir_can    * f_old_to_new_grow_pfts   &
!                                       & +  alb_nir_can_lost_avg    * f_d_to_new_grow_pfts
      IF (PRESENT(faPAR))   faPAR = faPAR    * f_old_to_new_grow_pfts   &
                                      & +  faPAR_lost_avg    * f_d_to_new_grow_pfts
      IF (PRESENT(soil_reflectivity_par))   soil_reflectivity_par = soil_reflectivity_par    * f_old_to_new_grow_pfts   &
                                      & +  soil_reflectivity_par_lost_avg    * f_d_to_new_grow_pfts
      ! IF (PRESENT(apar_nacc))   apar_nacc = apar_nacc    * f_old_to_new_grow_pfts   &
      !                                 & +  apar_nacc_lost_avg    * f_d_to_new_grow_pfts
      IF (PRESENT(faPAR_cl))   faPAR_cl = faPAR_cl    * f_old_to_new_grow_pfts   &
                                      & +  faPAR_cl_lost_avg    * f_d_to_new_grow_pfts
      IF (PRESENT(apar_per_lai_cl))  apar_per_lai_cl  = apar_per_lai_cl    * f_old_to_new_grow_pfts   &
                                      & +  apar_per_lai_cl_lost_avg    * f_d_to_new_grow_pfts
      ! IF (PRESENT(apar_per_lai_cl_nacc))   apar_per_lai_cl_nacc = apar_per_lai_cl_nacc    * f_old_to_new_grow_pfts   &
      !                                 & +  apar_per_lai_cl_nacc_lost_avg    * f_d_to_new_grow_pfts
      ! IF (PRESENT(apar_per_lai_cl_tavg))   apar_per_lai_cl_tavg = apar_per_lai_cl_tavg    * f_old_to_new_grow_pfts   &
      !                                 & +  apar_per_lai_cl_tavg_lost_avg    * f_d_to_new_grow_pfts

      ! TURB
      IF (PRESENT(fact_q_air))  fact_q_air = fact_q_air    * f_old_to_new_grow_pfts   &
                                      & +  fact_q_air_lost_avg    * f_d_to_new_grow_pfts
      IF (PRESENT(fact_qsat_srf))  fact_qsat_srf = fact_qsat_srf    * f_old_to_new_grow_pfts   &
                                      & +  fact_qsat_srf_lost_avg    * f_d_to_new_grow_pfts
      IF (PRESENT(fact_qsat_trans_srf))  fact_qsat_trans_srf = fact_qsat_trans_srf    * f_old_to_new_grow_pfts   &
                                      & +  fact_qsat_trans_srf_lost_avg    * f_d_to_new_grow_pfts

    END IF

  END SUBROUTINE

#endif
END MODULE mo_alcc_process
