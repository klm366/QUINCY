!>
!! @brief Contains the routines for the disturb processes
!!
!! @author
!!  Rainer Schneck, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Rainer Schneck (2015-07-29)
!!
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/license .
!!

MODULE mo_disturb_process
#ifndef __NO_JSBACH__

  USE mo_kind,      ONLY: wp
  USE mo_exception, ONLY: finish

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: burned_fract_jsbach, broken_woody_fract_jsbach, get_relative_humidity_air
  PUBLIC :: burned_fract_dl, get_relative_humidity_air_v2    ! <Rackhun Son>

  CHARACTER(len=*), PARAMETER :: modname = 'mo_disturb_process'

CONTAINS

  PURE ELEMENTAL SUBROUTINE burned_fract_jsbach( &
    & fire_rel_hum_threshold,                    & ! in
    & fire_litter_threshold,                     & ! in
    & fire_minimum,                              & ! in
    & fire_tau,                                  & ! in
    & q_rel_air_climbuf,                         & ! in
    & fuel,                                      & ! in
    & burned_fract                               & ! inout
    & )

    ! Input Arguments
    REAL(wp), INTENT(in) ::     &
      & fire_rel_hum_threshold, &
      & fire_litter_threshold,  &
      & fire_minimum,           &
      & fire_tau,               &
      & q_rel_air_climbuf,      &
      & fuel

    ! Output Arguments
    REAL(wp), INTENT(inout) :: burned_fract

    ! Locals
    REAL(wp) :: delta_time_yr
                !fire_rel_hum_threshold, & ! maximal relative humidity for fire
                !fire_litter_threshold,  & ! minimal amount of litter [mol(C)/m^2(grid box)] for fire
                !fire_minimum_woody,     & ! minimal fraction of act_fpc of woody PFT to be burned each year
                !fire_minimum_grass,     & ! minimal fraction of act_fpc of grass PFT to be burned each year
                !fire_tau_woody,         & ! return period of fire for woody PFT [year] at 0% relative humidity
                !fire_tau_grass,         & ! return period of fire for grass PFT [year] at 0% relative humidity

    delta_time_yr = 1._wp / 365._wp

    burned_fract = fire_minimum * delta_time_yr

    IF (       q_rel_air_climbuf < fire_rel_hum_threshold   &
         & .AND. fuel > fire_litter_threshold               &
       ) THEN
      burned_fract = burned_fract                                                                        &
        &           + (fire_rel_hum_threshold - q_rel_air_climbuf) / (fire_tau * fire_rel_hum_threshold) &
        &             * delta_time_yr
    END IF

  END SUBROUTINE burned_fract_jsbach

  !!! DL based Fire model  <Rackhun Son>
  SUBROUTINE burned_fract_dl( &
    & t_air,       t_clim_o,      t_clim_m,       & ! in
    & q_rel_air,   q_rel_clim_o,  q_rel_clim_m,   & ! in
    & q_rel_air_climbuf, fire_rel_hum_threshold,  & ! in
    & q_air,       q_clim_o,      q_clim_m,       & ! in
    & wind_10m,                                   & ! in
    & rain,                                       & ! in
    & swl1,        swl1_clim_o,   swl1_clim_m,    & ! in
    & swl2,        swl2_clim_o,   swl2_clim_m,    & ! in
    & swl3,        swl3_clim_o,   swl3_clim_m,    & ! in
    & swl4,        swl4_clim_o,   swl4_clim_m,    & ! in
    & lai,         lai_clim_o,    lai_clim_m,     & ! in
    & modified_pfts_fract,                        & ! in
    !& fuel,        fuel_yearly,                   & ! in
    & fuel, fire_litter_threshold,                & ! in
    !& population_density,                         & ! in
    & lightning_climatology,                      & ! in
    & elevation,                                  & ! in
    & slope,                                      & ! in
    & roughness,                                  & ! in
    !& topographic_position_idx,                   & ! in
    !& terrain_ruggedness_idx,                     & ! in
    !& vector_ruggedness_measure,                  & ! in
    !& profile_curvature,                          & ! in
    !& tangential_curvature,                       & ! in

    & population_density,                         & ! in
    & gdp,                                        & ! in
    & hdi,                                        & ! in
    & road_density,                               & ! in
    & land_use_state_c3ann,                       & ! in
    & land_use_state_c3nfx,                       & ! in
    & land_use_state_c3per,                       & ! in
    & land_use_state_c4ann,                       & ! in
    & land_use_state_c4per,                       & ! in
    & land_use_state_pastr,                       & ! in
    & land_use_state_primf,                       & ! in
    & land_use_state_primn,                       & ! in
    & land_use_state_range,                       & ! in
    & land_use_state_secdf,                       & ! in
    & land_use_state_secdn,                       & ! in
    & land_use_state_secma,                       & ! in
    & land_use_state_secmb,                       & ! in
    & land_use_state_urban,                       & ! in

    & nc,                                         & ! in
    & tilename,                                   & ! in

    & d_climate_input,     d_climate_hidden,      & ! in
    & d_vegetation_input,  d_vegetation_hidden,   & ! in
    !& d_lnt_input,         d_lnt_hidden,          & ! in
    & d_ant_input,         d_ant_hidden,          & ! in
    & d_output,                                   & ! in
    & d_emb_input,         d_emb_hidden,          & ! in
    & d_emb_output,                               & ! in

    & climate_i2i_weight, climate_i2i_bias,       & ! in
    & climate_h2i_weight, climate_h2i_bias,       & ! in
    & climate_i2f_weight, climate_i2f_bias,       & ! in
    & climate_h2f_weight, climate_h2f_bias,       & ! in
    & climate_i2g_weight, climate_i2g_bias,       & ! in
    & climate_h2g_weight, climate_h2g_bias,       & ! in
    & climate_i2o_weight, climate_i2o_bias,       & ! in
    & climate_h2o_weight, climate_h2o_bias,       & ! in
    & climate_out_weight, climate_out_bias,       & ! in

    & vegetation_i2i_weight, vegetation_i2i_bias, & ! in
    & vegetation_h2i_weight, vegetation_h2i_bias, & ! in
    & vegetation_i2f_weight, vegetation_i2f_bias, & ! in
    & vegetation_h2f_weight, vegetation_h2f_bias, & ! in
    & vegetation_i2g_weight, vegetation_i2g_bias, & ! in
    & vegetation_h2g_weight, vegetation_h2g_bias, & ! in
    & vegetation_i2o_weight, vegetation_i2o_bias, & ! in
    & vegetation_h2o_weight, vegetation_h2o_bias, & ! in
    & vegetation_out_weight, vegetation_out_bias, & ! in

    !& lig_topo_layer1_weight,  lig_topo_layer1_bias,  & ! in
    !& lig_topo_layer2_weight,  lig_topo_layer2_bias,  & ! in

    & anthropogenic_layer1_weight,  anthropogenic_layer1_bias,  & ! in
    & anthropogenic_layer2_weight,  anthropogenic_layer2_bias,  & ! in

    & climate_emb_layer1_weight,    climate_emb_layer1_bias,      & ! in
    & climate_emb_layer2_weight,    climate_emb_layer2_bias,      & ! in
    & vegetation_emb_layer1_weight, vegetation_emb_layer1_bias,   & ! in
    & vegetation_emb_layer2_weight, vegetation_emb_layer2_bias,   & ! in

    & climate_min_max,               & ! in 
    & vegetation_min_max,            & ! in
    & embedding_min_max,             & ! in
    !& lig_topo_min_max,              & ! in
    & anthropogenic_min_max,         & ! in

    & climate_lstm_hidden,                        & ! inout
    & climate_lstm_cell,                          & ! inout
    & vegetation_lstm_hidden,                     & ! inout
    & vegetation_lstm_cell,                       & ! inout
    & burned_fract                                & ! inout
    & )

    ! Input Arguments
    REAL(wp), INTENT(in) ::                             &
      & t_air(:),     t_clim_o(:),     t_clim_m(:),     &
      & q_rel_air(:), q_rel_clim_o(:), q_rel_clim_m(:), &
      & q_rel_air_climbuf(:), fire_rel_hum_threshold,   &
      & q_air(:),     q_clim_o(:),     q_clim_m(:),     &
      & wind_10m(:),                                    &
      & rain(:),                                        &
      & swl1(:),      swl1_clim_o(:),  swl1_clim_m(:),  &
      & swl2(:),      swl2_clim_o(:),  swl2_clim_m(:),  &
      & swl3(:),      swl3_clim_o(:),  swl3_clim_m(:),  &
      & swl4(:),      swl4_clim_o(:),  swl4_clim_m(:),  &
      & lai(:),       lai_clim_o(:),   lai_clim_m(:),   &
      & modified_pfts_fract(:,:),                       &
      !& fuel(:),       fuel_yearly(:),                  &
      & fuel(:), fire_litter_threshold,                 &
      & lightning_climatology(:),                       &
      & elevation(:),                                   &
      & slope(:),                                       &
      & roughness(:),                                   &
      !& topographic_position_idx(:),                    &
      !& terrain_ruggedness_idx(:),                      &
      !& vector_ruggedness_measure(:),                   &
      !& profile_curvature(:),                           &
      !& tangential_curvature(:),                        &

      & population_density(:),                          &
      & gdp(:),                                               &
      & hdi(:),                                               &
      & road_density(:),                                      &
      & land_use_state_c3ann(:),                              &
      & land_use_state_c3nfx(:),                              &
      & land_use_state_c3per(:),                              &
      & land_use_state_c4ann(:),                              &
      & land_use_state_c4per(:),                              &
      & land_use_state_pastr(:),                              &
      & land_use_state_primf(:),                              &
      & land_use_state_primn(:),                              &
      & land_use_state_range(:),                              &
      & land_use_state_secdf(:),                              &
      & land_use_state_secdn(:),                              &
      & land_use_state_secma(:),                              &
      & land_use_state_secmb(:),                              &
      & land_use_state_urban(:),                              &

      & climate_i2i_weight(:,:),    climate_i2i_bias(:),      &
      & climate_h2i_weight(:,:),    climate_h2i_bias(:),      &
      & climate_i2f_weight(:,:),    climate_i2f_bias(:),      &
      & climate_h2f_weight(:,:),    climate_h2f_bias(:),      &
      & climate_i2g_weight(:,:),    climate_i2g_bias(:),      &
      & climate_h2g_weight(:,:),    climate_h2g_bias(:),      &
      & climate_i2o_weight(:,:),    climate_i2o_bias(:),      &
      & climate_h2o_weight(:,:),    climate_h2o_bias(:),      &
      & climate_out_weight(:,:),    climate_out_bias(:),      &
      
      & vegetation_i2i_weight(:,:), vegetation_i2i_bias(:),   &
      & vegetation_h2i_weight(:,:), vegetation_h2i_bias(:),   &
      & vegetation_i2f_weight(:,:), vegetation_i2f_bias(:),   &
      & vegetation_h2f_weight(:,:), vegetation_h2f_bias(:),   &
      & vegetation_i2g_weight(:,:), vegetation_i2g_bias(:),   &
      & vegetation_h2g_weight(:,:), vegetation_h2g_bias(:),   &
      & vegetation_i2o_weight(:,:), vegetation_i2o_bias(:),   &
      & vegetation_h2o_weight(:,:), vegetation_h2o_bias(:),   &
      & vegetation_out_weight(:,:), vegetation_out_bias(:),   &

      !& lig_topo_layer1_weight(:,:),              &
      !& lig_topo_layer1_bias(:),                  &
      !& lig_topo_layer2_weight(:,:),              &
      !& lig_topo_layer2_bias(:),                  &

      & anthropogenic_layer1_weight(:,:),         &
      & anthropogenic_layer1_bias(:),             &
      & anthropogenic_layer2_weight(:,:),         &
      & anthropogenic_layer2_bias(:),             &

      & climate_emb_layer1_weight(:,:), climate_emb_layer1_bias(:),       &
      & climate_emb_layer2_weight(:,:), climate_emb_layer2_bias(:),       &

      & vegetation_emb_layer1_weight(:,:), vegetation_emb_layer1_bias(:), &
      & vegetation_emb_layer2_weight(:,:), vegetation_emb_layer2_bias(:), &

      & climate_min_max(:,:),              &
      & vegetation_min_max(:,:),           &
      & embedding_min_max(:,:),            &
      !& lig_topo_min_max(:,:),             &
      & anthropogenic_min_max(:,:)

    INTEGER, INTENT(in) ::                          &
      & nc,                                         & ! in
      & d_climate_input,     d_climate_hidden,      & ! in
      & d_vegetation_input,  d_vegetation_hidden,   & ! in
      !& d_lnt_input,         d_lnt_hidden,          & ! in
      & d_ant_input,         d_ant_hidden,          & ! in
      & d_output,                                   & ! in
      & d_emb_input,         d_emb_hidden,          & ! in
      & d_emb_output                                  ! in

    CHARACTER(len=*), INTENT(in)::                  &
      & tilename                                      ! in

    ! Output Arguments
    REAL(wp), INTENT(inout) ::                      &
      & climate_lstm_hidden(:,:),                   & ! inout
      & climate_lstm_cell(:,:),                     & ! inout
      & vegetation_lstm_hidden(:,:),                & ! inout
      & vegetation_lstm_cell(:,:),                  & ! inout
      & burned_fract(:)                               ! inout

    ! Locals
    REAL(wp) ::                                          &
      climate_variables(nc, d_climate_input),            &
      climate_lstm_i(nc, d_climate_hidden),              &
      climate_lstm_f(nc, d_climate_hidden),              &
      climate_lstm_g(nc, d_climate_hidden),              &
      climate_lstm_o(nc, d_climate_hidden),              &
      climate_lstm_output(nc, d_output*d_emb_output),    &
      climate_lstm_output_(nc, d_output),                &
      climate_score(nc, d_emb_output),                   &

      vegetation_variables(nc, d_vegetation_input),      &
      vegetation_lstm_i(nc, d_vegetation_hidden),        &
      vegetation_lstm_f(nc, d_vegetation_hidden),        &
      vegetation_lstm_g(nc, d_vegetation_hidden),        &
      vegetation_lstm_o(nc, d_vegetation_hidden),        &
      vegetation_lstm_output(nc, d_output*d_emb_output), &
      vegetation_lstm_output_(nc, d_output),             &
      vegetation_score(nc, d_emb_output),                &

      !lig_topo_variables(nc, d_lnt_input),               &
      !lig_topo_network(nc, d_lnt_hidden),                &
      !lig_topo_network_act(nc, d_lnt_hidden),            &
      !lig_topo_network_output(nc, d_output),             &

      anthro_variables(nc, d_ant_input),                 &
      anthro_network(nc, d_ant_hidden),                  &
      anthro_network_act(nc, d_ant_hidden),              &
      anthro_network_output(nc, d_output),               &

      emb_variables(nc, d_emb_input),                    &
      emb_network(nc, d_emb_hidden),                     &
      emb_network_act(nc, d_emb_hidden),                 &
      climate_emb_output(nc, d_emb_output),              &
      vegetation_emb_output(nc, d_emb_output),           &

      no_fire_cover(nc),                                 &
      fuel_(nc),                                         &
      q_rel_air_(nc)

    REAL(wp) ::  scale
    INTEGER  ::  pft_idx

    !climate_variables(:,1) = (t_air                 - t_air_minmax(1))                     &
    !climate_variables(:,1) = (t_air-t_clim_m+t_clim_o - t_air_minmax(1))             &
    !                       / (t_air_minmax(2)       - t_air_minmax(1))
    !climate_variables(:,2) = (t_air-t_clim_m        - t_air_ano_minmax(1))           &
    !                       / (t_air_ano_minmax(2)   - t_air_ano_minmax(1)) 
    !climate_variables(:,3) = (q_rel_air             - q_rel_minmax(1))                     &
    !climate_variables(:,3) = (q_rel_air-q_rel_clim_m+q_rel_clim_o - q_rel_minmax(1)) &
    !                       / (q_rel_minmax(2)       - q_rel_minmax(1)) 
    !climate_variables(:,4) = (q_rel_air-q_rel_clim_m- q_rel_ano_minmax(1))           &
    !                       / (q_rel_ano_minmax(2)   - q_rel_ano_minmax(1)) 
    !climate_variables(:,5) = (q_air                 - q_minmax(1))                         &
    !climate_variables(:,5) = (q_air-q_clim_m+q_clim_o - q_minmax(1))                 &
    !                       / (q_minmax(2)           - q_minmax(1)) 
    !climate_variables(:,6) = (q_air-q_clim_m        - q_ano_minmax(1))               &
    !                       / (q_ano_minmax(2)       - q_ano_minmax(1)) 
    !climate_variables(:,7) = (wind_10m              - wind_10m_minmax(1))            &
    !                       / (wind_10m_minmax(2)    - wind_10m_minmax(1)) 
    !climate_variables(:,8) = (rain*86400            - rain_minmax(1))                &
    !                       / (rain_minmax(2)        - rain_minmax(1)) 
    !climate_variables(:,9) = (swl1                  - swl1_minmax(1))                      &
    !climate_variables(:,9) = (swl1-swl1_clim_m+swl1_clim_o - swl1_minmax(1))         &
    !                       / (swl1_minmax(2)        - swl1_minmax(1)) 
    !climate_variables(:,10)= (swl1-swl1_clim_m      - swl1_ano_minmax(1))            &
    !                       / (swl1_ano_minmax(2)    - swl1_ano_minmax(1)) 
    !climate_variables(:,11)= (swl2                  - swl2_minmax(1))                      &
    !climate_variables(:,11)= (swl2-swl2_clim_m+swl2_clim_o - swl2_minmax(1))         &
    !                       / (swl2_minmax(2)        - swl2_minmax(1))
    !climate_variables(:,12)= (swl2-swl2_clim_m      - swl2_ano_minmax(1))            &
    !                       / (swl2_ano_minmax(2)    - swl2_ano_minmax(1))  
    !climate_variables(:,13)= (swl3                  - swl3_minmax(1))                      &
    !climate_variables(:,13)= (swl3-swl3_clim_m+swl3_clim_o - swl3_minmax(1))         &
    !                       / (swl3_minmax(2)        - swl3_minmax(1))
    !climate_variables(:,14)= (swl3-swl3_clim_m      - swl3_ano_minmax(1))            &
    !                       / (swl3_ano_minmax(2)    - swl3_ano_minmax(1))  
    !climate_variables(:,15)= (swl4                  - swl4_minmax(1))                      &
    !climate_variables(:,15)= (swl4-swl4_clim_m+swl4_clim_o - swl4_minmax(1))         &
    !                       / (swl4_minmax(2)        - swl4_minmax(1))
    !climate_variables(:,16)= (swl4-swl4_clim_m      - swl4_ano_minmax(1))            &
    !                       / (swl4_ano_minmax(2)    - swl4_ano_minmax(1)) 
    
    !vegetation_variables(:,1) = (lai-lai_clim_m+lai_clim_o - lai_minmax(1))          &
    !                          / (lai_minmax(2)    - lai_minmax(1))
    !vegetation_variables(:,2) = (lai-lai_clim_m   - lai_ano_minmax(1))               &
    !                          / (lai_ano_minmax(2)- lai_ano_minmax(1))

    !static_variables(:,1) = (population_density             - population_density_minmax(1))    &
    !                      / (population_density_minmax(2)   - population_density_minmax(1))
    !static_variables(:,2) = (lightning_climatology          - lightning_climatology_minmax(1)) &
    !                      / (lightning_climatology_minmax(2)- lightning_climatology_minmax(1))
    !static_variables(:,3) = (elevation                      - elevation_minmax(1))             &
    !                      / (elevation_minmax(2)            - elevation_minmax(1))  

    !emb_variables(:,1) = (t_clim_o             - t_air_clim_minmax(1)) &
    !                   / (t_air_clim_minmax(2) - t_air_clim_minmax(1)) 
    !emb_variables(:,2) = (q_rel_clim_o         - q_rel_clim_minmax(1)) &
    !                   / (q_rel_clim_minmax(2) - q_rel_clim_minmax(1))
    !emb_variables(:,3) = (q_clim_o             - q_clim_minmax(1))     &
    !                   / (q_clim_minmax(2)     - q_clim_minmax(1))
    !emb_variables(:,4) = (swl1_clim_o          - swl1_clim_minmax(1))  &
    !                   / (swl1_clim_minmax(2)  - swl1_clim_minmax(1)) 
    !emb_variables(:,5) = (swl2_clim_o          - swl2_clim_minmax(1))  &
    !                   / (swl2_clim_minmax(2)  - swl2_clim_minmax(1)) 
    !emb_variables(:,6) = (swl3_clim_o          - swl3_clim_minmax(1))  &
    !                   / (swl3_clim_minmax(2)  - swl3_clim_minmax(1)) 
    !emb_variables(:,7) = (swl4_clim_o          - swl4_clim_minmax(1))  &
    !                   / (swl4_clim_minmax(2)  - swl4_clim_minmax(1))
    !emb_variables(:,8) = (lai_clim_o           - lai_clim_minmax(1))   &
    !                   / (lai_clim_minmax(2)   - lai_clim_minmax(1)) 
    !emb_variables(:,9) = (fuel                 - agb_minmax(1))        &
    !                   / (agb_minmax(2)        - agb_minmax(1)) 
    !emb_variables(:,9) = (fuel_yearly          - agb_minmax(1))        &
    !                   / (agb_minmax(2)        - agb_minmax(1)) 
    !emb_variables(:,10)= modified_pfts_fract(:,1)
    !emb_variables(:,10)= (modified_pfts_fract(:,1) - pft1_minmax(1))   &
    !                   / (pft1_minmax(2)           - pft1_minmax(1))
    !emb_variables(:,11)= modified_pfts_fract(:,2)
    !emb_variables(:,11)= (modified_pfts_fract(:,2) - pft2_minmax(1))   &
    !                   / (pft2_minmax(2)           - pft2_minmax(1))
    !emb_variables(:,12)= modified_pfts_fract(:,3)
    !emb_variables(:,12)= (modified_pfts_fract(:,3) - pft3_minmax(1))   &
    !                   / (pft3_minmax(2)           - pft3_minmax(1))
    !emb_variables(:,13)= modified_pfts_fract(:,4)
    !emb_variables(:,13)= (modified_pfts_fract(:,4) - pft4_minmax(1))   &
    !                   / (pft4_minmax(2)           - pft4_minmax(1))
    !emb_variables(:,14)= modified_pfts_fract(:,5)
    !emb_variables(:,14)= (modified_pfts_fract(:,5) - pft5_minmax(1))   &
    !                   / (pft5_minmax(2)           - pft5_minmax(1))
    !emb_variables(:,15)= modified_pfts_fract(:,6)
    !emb_variables(:,15)= (modified_pfts_fract(:,6) - pft6_minmax(1))   &
    !                   / (pft6_minmax(2)           - pft6_minmax(1))
    !emb_variables(:,16)= modified_pfts_fract(:,7)
    !emb_variables(:,16)= (modified_pfts_fract(:,7) - pft7_minmax(1))   &
    !                   / (pft7_minmax(2)           - pft7_minmax(1))
    !emb_variables(:,17)= modified_pfts_fract(:,8)
    !emb_variables(:,17)= (modified_pfts_fract(:,8) - pft8_minmax(1))   &
    !                   / (pft8_minmax(2)           - pft8_minmax(1))
    !emb_variables(:,18)= modified_pfts_fract(:,9)
    !emb_variables(:,18)= (modified_pfts_fract(:,9) - pft9_minmax(1))   &
    !                   / (pft9_minmax(2)           - pft9_minmax(1))


    climate_variables(:,1) = t_air - t_clim_m + t_clim_o
    climate_variables(:,2) = t_air - t_clim_m
    climate_variables(:,3) = q_rel_air - q_rel_clim_m + q_rel_clim_o
    climate_variables(:,4) = q_rel_air - q_rel_clim_m
    climate_variables(:,5) = q_air - q_clim_m + q_clim_o
    climate_variables(:,6) = q_air - q_clim_m
    climate_variables(:,7) = wind_10m 
    climate_variables(:,8) = rain*86400
    !climate_variables(:,1) = t_air
    !climate_variables(:,2) = t_air - t_clim_m
    !climate_variables(:,3) = q_rel_air
    !climate_variables(:,4) = q_rel_air - q_rel_clim_m
    !climate_variables(:,5) = q_air
    !climate_variables(:,6) = q_air - q_clim_m
    !climate_variables(:,7) = wind_10m 
    !climate_variables(:,8) = rain*86400
    !climate_variables(:,9) = swl1 - swl1_clim_m + swl1_clim_o
    !climate_variables(:,10)= swl1 - swl1_clim_m
    !climate_variables(:,11)= swl2 - swl2_clim_m + swl2_clim_o
    !climate_variables(:,12)= swl2 - swl2_clim_m
    !climate_variables(:,13)= swl3 - swl3_clim_m + swl3_clim_o
    !climate_variables(:,14)= swl3 - swl3_clim_m
    !climate_variables(:,15)= swl4 - swl4_clim_m + swl4_clim_o
    !climate_variables(:,16)= swl4 - swl4_clim_m

    !climate_variables(:,1) = t_air
    !climate_variables(:,2) = t_air - t_clim_m
    !climate_variables(:,3) = q_rel_air
    !climate_variables(:,4) = q_rel_air - q_rel_clim_m
    !climate_variables(:,5) = q_air
    !climate_variables(:,6) = q_air - q_clim_m
    !climate_variables(:,7) = wind_10m 
    !climate_variables(:,8) = rain*86400
    !climate_variables(:,9) = swl1
    !climate_variables(:,10)= swl1 - swl1_clim_m
    !climate_variables(:,11)= swl2
    !climate_variables(:,12)= swl2 - swl2_clim_m
    !climate_variables(:,13)= swl3
    !climate_variables(:,14)= swl3 - swl3_clim_m
    !climate_variables(:,15)= swl4
    !climate_variables(:,16)= swl4 - swl4_clim_m

    vegetation_variables(:,1) = lai  - lai_clim_m  + lai_clim_o
    vegetation_variables(:,2) = lai  - lai_clim_m
    vegetation_variables(:,3) = swl1 - swl1_clim_m + swl1_clim_o
    vegetation_variables(:,4) = swl1 - swl1_clim_m
    vegetation_variables(:,5) = swl2 - swl2_clim_m + swl2_clim_o
    vegetation_variables(:,6) = swl2 - swl2_clim_m
    vegetation_variables(:,7) = swl3 - swl3_clim_m + swl3_clim_o
    vegetation_variables(:,8) = swl3 - swl3_clim_m
    vegetation_variables(:,9) = swl4 - swl4_clim_m + swl4_clim_o
    vegetation_variables(:,10)= swl4 - swl4_clim_m
    !vegetation_variables(:,1) = lai
    !vegetation_variables(:,2) = lai  - lai_clim_m
    !vegetation_variables(:,3) = swl1
    !vegetation_variables(:,4) = swl1 - swl1_clim_m
    !vegetation_variables(:,5) = swl2
    !vegetation_variables(:,6) = swl2 - swl2_clim_m
    !vegetation_variables(:,7) = swl3
    !vegetation_variables(:,8) = swl3 - swl3_clim_m
    !vegetation_variables(:,9) = swl4
    !vegetation_variables(:,10)= swl4 - swl4_clim_m

    !lig_topo_variables(:,1) = lightning_climatology
    !lig_topo_variables(:,2) = elevation
    !lig_topo_variables(:,3) = slope
    !lig_topo_variables(:,4) = roughness
    !lig_topo_variables(:,5) = topographic_position_idx
    !lig_topo_variables(:,6) = terrain_ruggedness_idx
    !lig_topo_variables(:,7) = vector_ruggedness_measure
    !lig_topo_variables(:,8) = profile_curvature
    !lig_topo_variables(:,9) = tangential_curvature

    anthro_variables(:,1) = population_density
    anthro_variables(:,2) = gdp
    anthro_variables(:,3) = hdi
    anthro_variables(:,4) = road_density
    anthro_variables(:,5) = land_use_state_c3ann
    anthro_variables(:,6) = land_use_state_c3nfx
    anthro_variables(:,7) = land_use_state_c3per
    anthro_variables(:,8) = land_use_state_c4ann
    anthro_variables(:,9) = land_use_state_c4per
    anthro_variables(:,10)= land_use_state_pastr
    anthro_variables(:,11)= land_use_state_primf
    anthro_variables(:,12)= land_use_state_primn
    anthro_variables(:,13)= land_use_state_range
    anthro_variables(:,14)= land_use_state_secdf
    anthro_variables(:,15)= land_use_state_secdn
    anthro_variables(:,16)= land_use_state_secma
    anthro_variables(:,17)= land_use_state_secmb
    anthro_variables(:,18)= land_use_state_urban

    emb_variables(:,1) = t_clim_o
    emb_variables(:,2) = q_rel_clim_o
    emb_variables(:,3) = q_clim_o
    emb_variables(:,4) = swl1_clim_o
    emb_variables(:,5) = swl2_clim_o
    emb_variables(:,6) = swl3_clim_o 
    emb_variables(:,7) = swl4_clim_o
    emb_variables(:,8) = lai_clim_o
    !emb_variables(:,1) = t_clim_m
    !emb_variables(:,2) = q_rel_clim_m
    !emb_variables(:,3) = q_clim_m
    !emb_variables(:,4) = swl1_clim_m
    !emb_variables(:,5) = swl2_clim_m
    !emb_variables(:,6) = swl3_clim_m 
    !emb_variables(:,7) = swl4_clim_m
    !emb_variables(:,8) = lai_clim_m
    emb_variables(:,9) = fuel
    emb_variables(:,10)= lightning_climatology
    emb_variables(:,11)= elevation
    emb_variables(:,12)= slope
    emb_variables(:,13)= roughness
    !emb_variables(:,10)= modified_pfts_fract(:,1)
    !emb_variables(:,11)= modified_pfts_fract(:,2)
    !emb_variables(:,12)= modified_pfts_fract(:,3)
    !emb_variables(:,13)= modified_pfts_fract(:,4)
    !emb_variables(:,14)= modified_pfts_fract(:,5)
    !emb_variables(:,15)= modified_pfts_fract(:,6)
    !emb_variables(:,16)= modified_pfts_fract(:,7)
    !emb_variables(:,17)= modified_pfts_fract(:,8)
    !emb_variables(:,18)= modified_pfts_fract(:,9)
    emb_variables(:,14)= modified_pfts_fract(:,1)
    emb_variables(:,15)= modified_pfts_fract(:,2)
    emb_variables(:,16)= modified_pfts_fract(:,3)
    emb_variables(:,17)= modified_pfts_fract(:,4)
    emb_variables(:,18)= modified_pfts_fract(:,5)
    emb_variables(:,19)= modified_pfts_fract(:,6)
    emb_variables(:,20)= modified_pfts_fract(:,7)
    emb_variables(:,21)= modified_pfts_fract(:,8)
    emb_variables(:,22)= modified_pfts_fract(:,9)

    !where (climate_variables(:,3) < 0.)  climate_variables(:,3)  = 0.
    !where (climate_variables(:,5) < 0.)  climate_variables(:,5)  = 0.
    !where (climate_variables(:,7) < 0.)  climate_variables(:,7)  = 0.
    !where (climate_variables(:,8) < 0.)  climate_variables(:,8)  = 0.
    !where (climate_variables(:,9) < 0.)  climate_variables(:,9)  = 0.
    !where (climate_variables(:,11) < 0.) climate_variables(:,11) = 0.
    !where (climate_variables(:,13) < 0.) climate_variables(:,13) = 0.
    !where (climate_variables(:,15) < 0.) climate_variables(:,15) = 0.
    !where (vegetation_variables(:,1) < 0.) vegetation_variables(:,1) = 0.

    ! Normalize
    climate_variables(:,:)    = (climate_variables(:,:)                    - spread(climate_min_max(:,1),1,nc))    &
                              / (spread(climate_min_max(:,2),1,nc)         - spread(climate_min_max(:,1),1,nc))
    vegetation_variables(:,:) = (vegetation_variables(:,:)                 - spread(vegetation_min_max(:,1),1,nc)) &
                              / (spread(vegetation_min_max(:,2),1,nc)      - spread(vegetation_min_max(:,1),1,nc))
    !emb_variables(:,1:9)      = (emb_variables(:,1:9)                      - spread(embedding_min_max(:,1),1,nc))  &
    !                          / (spread(embedding_min_max(:,2),1,nc)       - spread(embedding_min_max(:,1),1,nc))
    emb_variables(:,1:13)     = (emb_variables(:,1:13)                     - spread(embedding_min_max(:,1),1,nc))  &
                              / (spread(embedding_min_max(:,2),1,nc)       - spread(embedding_min_max(:,1),1,nc))
    !lig_topo_variables(:,:)   = (lig_topo_variables(:,:)                   - spread(lig_topo_min_max(:,1),1,nc))   &
    !                          / (spread(lig_topo_min_max(:,2),1,nc)        - spread(lig_topo_min_max(:,1),1,nc))
    anthro_variables(:,1:4)   = (anthro_variables(:,1:4)                   - spread(anthropogenic_min_max(1:4,1),1,nc))   &
                              / (spread(anthropogenic_min_max(1:4,2),1,nc) - spread(anthropogenic_min_max(1:4,1),1,nc))
    
    !if (maxval(emb_variables(:,1))>1.1 .or. minval(emb_variables(:,1))<-0.1) then
    !  print*, '1', minval(emb_variables(:,1)), maxval(emb_variables(:,1))
    !endif

    ! climate-LSTM model
    climate_lstm_i = matmul(climate_variables,   climate_i2i_weight) + spread(climate_i2i_bias, 1, nc) &
                   + matmul(climate_lstm_hidden, climate_h2i_weight) + spread(climate_h2i_bias, 1, nc)
    climate_lstm_f = matmul(climate_variables,   climate_i2f_weight) + spread(climate_i2f_bias, 1, nc) &
                   + matmul(climate_lstm_hidden, climate_h2f_weight) + spread(climate_h2f_bias, 1, nc)
    climate_lstm_g = matmul(climate_variables,   climate_i2g_weight) + spread(climate_i2g_bias, 1, nc) &
                   + matmul(climate_lstm_hidden, climate_h2g_weight) + spread(climate_h2g_bias, 1, nc)
    climate_lstm_o = matmul(climate_variables,   climate_i2o_weight) + spread(climate_i2o_bias, 1, nc) &
                   + matmul(climate_lstm_hidden, climate_h2o_weight) + spread(climate_h2o_bias, 1, nc)
         
    climate_lstm_i = 1/(1+exp(-climate_lstm_i))
    climate_lstm_f = 1/(1+exp(-climate_lstm_f))
    climate_lstm_g = tanh(climate_lstm_g)
    climate_lstm_o = 1/(1+exp(-climate_lstm_o))

    climate_lstm_cell   = climate_lstm_f*climate_lstm_cell + climate_lstm_i*climate_lstm_g 
    climate_lstm_hidden = climate_lstm_o*tanh(climate_lstm_cell)
    climate_lstm_output = matmul(climate_lstm_hidden, climate_out_weight) + spread(climate_out_bias, 1, nc)


    ! climate emb NN model
    emb_network = matmul(emb_variables, climate_emb_layer1_weight) + spread(climate_emb_layer1_bias, 1, nc)
    emb_network_act = log(1+exp(emb_network))
    where (emb_network > 20) emb_network_act = emb_network
    climate_emb_output = matmul(emb_network_act, climate_emb_layer2_weight) + spread(climate_emb_layer2_bias, 1, nc)
    where (climate_emb_output > 700._wp) climate_emb_output = 700._wp

    climate_score = exp(climate_emb_output)
    climate_score = climate_score / spread(sum(climate_score,dim=2), 2, d_emb_output)

    climate_lstm_output_ = sum(                                      &
      reshape(climate_lstm_output, (/nc, d_output, d_emb_output/))   &
      *spread(climate_score, 2, d_output)                            &
      , dim=3)


    ! vegetation-LSTM model
    vegetation_lstm_i = matmul(vegetation_variables,   vegetation_i2i_weight) + spread(vegetation_i2i_bias, 1, nc) &
                      + matmul(vegetation_lstm_hidden, vegetation_h2i_weight) + spread(vegetation_h2i_bias, 1, nc)
    vegetation_lstm_f = matmul(vegetation_variables,   vegetation_i2f_weight) + spread(vegetation_i2f_bias, 1, nc) &
                      + matmul(vegetation_lstm_hidden, vegetation_h2f_weight) + spread(vegetation_h2f_bias, 1, nc)
    vegetation_lstm_g = matmul(vegetation_variables,   vegetation_i2g_weight) + spread(vegetation_i2g_bias, 1, nc) &
                      + matmul(vegetation_lstm_hidden, vegetation_h2g_weight) + spread(vegetation_h2g_bias, 1, nc)
    vegetation_lstm_o = matmul(vegetation_variables,   vegetation_i2o_weight) + spread(vegetation_i2o_bias, 1, nc) &
                      + matmul(vegetation_lstm_hidden, vegetation_h2o_weight) + spread(vegetation_h2o_bias, 1, nc)

    vegetation_lstm_i = 1/(1+exp(-vegetation_lstm_i))
    vegetation_lstm_f = 1/(1+exp(-vegetation_lstm_f))
    vegetation_lstm_g = tanh(vegetation_lstm_g)
    vegetation_lstm_o = 1/(1+exp(-vegetation_lstm_o))

    vegetation_lstm_cell   = vegetation_lstm_f*vegetation_lstm_cell + vegetation_lstm_i*vegetation_lstm_g
    vegetation_lstm_hidden = vegetation_lstm_o*tanh(vegetation_lstm_cell)
    vegetation_lstm_output = matmul(vegetation_lstm_hidden, vegetation_out_weight) + spread(vegetation_out_bias, 1, nc)


    ! vegetation emb NN model
    emb_network = matmul(emb_variables, vegetation_emb_layer1_weight) + spread(vegetation_emb_layer1_bias, 1, nc)
    emb_network_act = log(1+exp(emb_network))
    where (emb_network > 20) emb_network_act = emb_network
    vegetation_emb_output = matmul(emb_network_act, vegetation_emb_layer2_weight) + spread(vegetation_emb_layer2_bias, 1, nc)
    where (vegetation_emb_output > 700._wp) vegetation_emb_output = 700._wp

    vegetation_score = exp(vegetation_emb_output)
    vegetation_score = vegetation_score / spread(sum(vegetation_score,dim=2), 2, d_emb_output)

    vegetation_lstm_output_ = sum(                                      &
      reshape(vegetation_lstm_output, (/nc, d_output, d_emb_output/))   &
      *spread(vegetation_score, 2, d_output)                            &
      , dim=3)
    

    ! lightning & topography NN model
    !lig_topo_network = matmul(lig_topo_variables, lig_topo_layer1_weight) + spread(lig_topo_layer1_bias, 1, nc)
    !lig_topo_network_act = log(1+exp(lig_topo_network))
    !where (lig_topo_network > 20) lig_topo_network_act = lig_topo_network
    !lig_topo_network_output = matmul(lig_topo_network_act, lig_topo_layer2_weight) + spread(lig_topo_layer2_bias, 1, nc)


    ! anthropogenic NN model
    anthro_network = matmul(anthro_variables, anthropogenic_layer1_weight) + spread(anthropogenic_layer1_bias, 1, nc)
    anthro_network_act = log(1+exp(anthro_network))
    where (anthro_network > 20) anthro_network_act = anthro_network
    anthro_network_output = matmul(anthro_network_act, anthropogenic_layer2_weight) + spread(anthropogenic_layer2_bias, 1, nc)


    ! output
    !IF (tilename == 'pft01') THEN 
    !  pft_idx = 2
    !  scale = 1._wp
    !ELSE IF (tilename == 'pft02') THEN
    !  pft_idx = 3
    !  scale = 1._wp
    !ELSE IF (tilename == 'pft03') THEN
    !  pft_idx = 4
    !  scale = 1._wp
    !ELSE IF (tilename == 'pft04') THEN
    !  pft_idx = 5
    !  scale = 1._wp
    !ELSE IF (tilename == 'pft05') THEN
    !  pft_idx = 6
    !  scale = 1._wp
    !ELSE IF (tilename == 'pft06') THEN
    !  pft_idx = 7
    !  scale = 1._wp
    !ELSE IF (tilename == 'pft07') THEN
    !  pft_idx = 8
    !  scale = 1._wp/2._wp
    !ELSE IF (tilename == 'pft08') THEN
    !  pft_idx = 8
    !  scale = 1._wp/2._wp
    !ELSE
    !  print*, 'new land cover type'
    !endif
    
    !burned_fract = climate_lstm_output_(:,pft_idx)                &
    !             * vegetation_lstm_output_(:,pft_idx)             &
    !             * lig_topo_network_output(:,pft_idx)             &
    !             * modified_pfts_fract(:,pft_idx)                 &
    !             * scale/100.

    !burned_fract = burned_fract                                   &
    !             + climate_lstm_output_(:,1)                      &
    !             * vegetation_lstm_output_(:,1)                   &
    !             * lig_topo_network_output(:,1)                   &
    !             * modified_pfts_fract(:,1)                       &
    !             * scale/100./8.
                 !* 8._wp / 100

    burned_fract = sum(climate_lstm_output_         &
                     * vegetation_lstm_output_      &
                     !* lig_topo_network_output      &
                     * anthro_network_output        &
                     * modified_pfts_fract(:,1:8), dim=2)

    where (burned_fract < 0.) burned_fract = 0.
    !where (burned_fract > 1.) burned_fract = 1.
    !where (burned_fract > 100.) burned_fract = 100.
    burned_fract(:) = burned_fract(:)/100.

    no_fire_cover(:) = modified_pfts_fract(:,9)+modified_pfts_fract(:,1)

    where (modified_pfts_fract(:,9) > 0.3) no_fire_cover = 1.
    !where (modified_pfts_fract(:,1) > 0.2) no_fire_cover = 1.

    !burned_fract(:) = burned_fract(:)*(1.-modified_pfts_fract(:,9)-modified_pfts_fract(:,1))
    burned_fract(:) = burned_fract(:)*(1.-no_fire_cover(:))

    !print*, maxval(modified_pfts_fract(:,9)+modified_pfts_fract(:,1)), maxval(modified_pfts_fract(:,9))

    fuel_(:) = emb_variables(:,9)
    !where (fuel_ < 0.)   fuel_ = 0.
    !where (fuel_ > 1.)   fuel_ = 1.
    !where (fuel < fire_litter_threshold/2) fuel_ = 0.
    !burned_fract = burned_fract*(fuel_(:)**2)  ! v1
    burned_fract(:) = burned_fract(:)*fuel_(:)  ! v1
    !burned_fract(:) = burned_fract(:)*emb_variables(:,9)  ! v1


    !q_rel_air_(:) = q_rel_air(:) - q_rel_clim_m(:) + q_rel_clim_o(:)
    !q_rel_air_(:) = q_rel_air_(:)/100.
    q_rel_air_(:) = q_rel_air(:)
    where (q_rel_air_ < 0.) q_rel_air_ = 0.
    where (q_rel_air_ > 100.) q_rel_air_ = 100.
    !where (q_rel_air_climbuf > fire_rel_hum_threshold) q_rel_air_ = 100.
    !where (q_rel_air_ > fire_rel_hum_threshold) q_rel_air_ = 100.
    

    !burned_fract = burned_fract*((1.-q_rel_air_(:)/100.)**2)  ! v1
    burned_fract(:) = burned_fract(:)*(1.-q_rel_air_(:)/100.)  ! v1
    !q_rel_air_(:) = 1._wp / (1._wp + exp(-20._wp*((q_rel_air_(:)/100.)-0.5)))
    !burned_fract = burned_fract*(1.-q_rel_air_(:))
    
    !burned_fract(:) = burned_fract(:)/100.
    !where (burned_fract < 0.) burned_fract = 0.
    where (burned_fract > 1.) burned_fract = 1.

    !if (maxval(burned_fract) .gt. 1e-3) then
    !  print*, maxval(burned_fract)*100, minval(burned_fract)*100
    !endif

    !if (maxval(burned_fract) .gt. 0) then
    !  print*, maxval(burned_fract)
    !endif

    !print*, MAXVAL(burned_fract), MAXVAL(burned_fract)  

    ! ratio for bared & snow area
    !burned_fract = burned_fract*(1-modified_pfts_fract(:,9)-modified_pfts_fract(:,1))

    ! ratio for dry (fuel)
    !q_rel_air_(:) = q_rel_air(:)/100._wp
    !q_rel_air_(:) = q_rel_air(:) - q_rel_clim_m(:) + q_rel_clim_o(:)
    !q_rel_air_(:) = q_rel_air_(:)/100._wp
    !where (q_rel_air_ < 0._wp) q_rel_air_ = 0._wp
    !where (q_rel_air_ > 1._wp) q_rel_air_ = 1._wp
    !q_rel_air_(:) = 1._wp / (1._wp + exp(-20._wp*(q_rel_air_(:)-0.5)))
    !q_rel_air_(:) = 1._wp / (1._wp + exp(-5._wp*(q_rel_air_(:)-0.5)))
    !q_rel_air_(:) = 1._wp / (1._wp + exp(-10._wp*(q_rel_air_(:)-0.5)))
    !q_rel_air_(:) = 1._wp / (1._wp + exp(-15._wp*(q_rel_air_(:)-0.5)))
    !q_rel_air_(:) = 1._wp / (1._wp + exp(-25._wp*(q_rel_air_(:)-0.5)))
    
  END SUBROUTINE burned_fract_dl

  ! Calculate windbreak for woody types
  PURE ELEMENTAL SUBROUTINE broken_woody_fract_jsbach( &
    & wind_threshold,                                  & ! in
    & wind_damage_scale,                               & ! in
    & cover_fract_pot,                                 & ! in
    & prev_day_max_wind_10m,                           & ! in
    & max_wind_10m,                                    & ! in
    & damaged_fract                                    & ! inout
    & )

    ! Arguments
    REAL(wp), INTENT(in) ::    &
      & wind_threshold,        &
      & wind_damage_scale,     &
      & cover_fract_pot,       &
      & prev_day_max_wind_10m, &
      & max_wind_10m

    REAL(wp), INTENT(inout) :: &
      & damaged_fract

    ! Locals
    REAL(wp) :: delta_time_yr

    delta_time_yr = 1._wp / 365._wp

    IF (      cover_fract_pot > EPSILON(1._wp)                      &
        .AND. prev_day_max_wind_10m > max_wind_10m * wind_threshold &
       ) THEN
      damaged_fract = MIN( 1._wp - EPSILON(1._wp) / cover_fract_pot,                                          &
        &                 delta_time_yr * wind_damage_scale * prev_day_max_wind_10m ** 3._wp / max_wind_10m   &
        &               )
    ENDIF

  END SUBROUTINE broken_woody_fract_jsbach

  FUNCTION get_relative_humidity_air(nc, air_moisture, air_temperature, air_press) RESULT(relative_humidity_air)

    USE mo_physical_constants_iface, ONLY: rd  ! gas constant for dry air [J/(K*kg)]
    USE mo_physical_constants_iface, ONLY: rv  ! gas constant for water vapor [J/(K*kg)]
    !USE mo_jsb4_forcing,             ONLY: sat_specific_humidity ! Function

    INTEGER, INTENT(in)  :: nc
    REAL(wp), INTENT(in) :: &
      air_moisture(:),      & ! Specific humidity of air [kg kg-1]
      air_temperature(:),   & ! Temperature of air [K]
      air_press(:)            ! Pressure of air [Pa]

    REAL(wp) :: air_qsat(SIZE(air_moisture)), &
                relative_humidity_air(SIZE(air_moisture))

    CHARACTER(len=*), PARAMETER :: routine = modname//':get_relative_humidity_air'

    CALL sat_specific_humidity( &
      & nc,                     & ! in
      & air_temperature,        & ! in
      & air_press,              & ! in
      & air_qsat                & ! out
      & )

        IF (ANY( air_qsat(:) > 0.99_wp*HUGE(1._wp) )) THEN
           CALL finish(TRIM(routine), 'lookup table overflow.')
        END IF

    relative_humidity_air(:) = 100._wp * air_moisture(:)                                       &
                               * ( (1._wp - rd / rv)  * air_qsat(:)  + (rd / rv) )             &
                               / ( ((1._wp - rd / rv) * air_moisture(:) + (rd / rv)) * air_qsat(:) )

  END FUNCTION get_relative_humidity_air

  FUNCTION get_relative_humidity_air_v2(nc, air_moisture, air_temperature, air_press) RESULT(relative_humidity_air)

    INTEGER, INTENT(in)  :: nc
    REAL(wp), INTENT(in) :: &
      air_moisture(:),      & ! Specific humidity of air [kg kg-1]
      air_temperature(:),   & ! Temperature of air [K]
      air_press(:)            ! Pressure of air [Pa]

    REAL(wp) :: relative_humidity_air(SIZE(air_moisture))

    CHARACTER(len=*), PARAMETER :: routine = modname//':get_relative_humidity_air_v2'

    relative_humidity_air(:) = 0.263 * air_press(:) * air_moisture(:)                          &
                             / ( EXP((17.67*(air_temperature(:)-273.16))/(air_temperature(:)-29.65)) ) 

  END FUNCTION get_relative_humidity_air_v2

  ! R: TBD: I copied this subroutine from MO_JSB4_FORCING, because
  !    this mo is not compiled with icon/echam.
  !    Reiner has to decide if we change config or the place of this procedure.
  !    I made some changes because the dimensions (vector lengths) were unnecessary
  !    cryptic. We could even skip "nc" if we would change the DO loop below...
  SUBROUTINE sat_specific_humidity(nc, temp, pressure, qsat)
    ! Returns saturation specific humidity for given temperature and pressure (at some atmospheric level or at surface)
    ! Uses Eq. 2.27 of Fundamentals of Atmospheric Modeling for the saturated case, but the saturation vapor pressure
    ! of water over a liquid surface resp. ice (see pp 32-34 of Ref.) is computed as in ECHAM5 (mo_convect_tables)
    USE mo_exception,                ONLY: message_text
    USE mo_jsb_convect_tables_iface, ONLY: tlucua, jptlucu1, jptlucu2
    USE mo_jsb_physical_constants,   ONLY: rvd1

    INTEGER,  INTENT(IN) :: nc
    REAL(wp), INTENT(IN) :: temp(nc)      ! Air temperature at level [K]
    REAL(wp), INTENT(IN) :: pressure(nc)  ! Pressure at level [Pa]
    REAL(wp), INTENT(OUT):: qsat(nc)

    CHARACTER(len=*), PARAMETER :: routine = modname//':sat_specific_humidity'

    INTEGER  :: i
    INTEGER  :: it(nc)
    REAL(wp) :: tluc(nc)

    !Assertion: the expected temperature range is indirectly defined by the tlucua bounds -> consistency check
    DO i = 1,nc
      IF ((temp(i) .LT. REAL(LBOUND(tlucua,dim=1)/1000._wp)) .OR. (temp(i) .GT. REAL(UBOUND(tlucua,dim=1)/1000._wp))) THEN
        WRITE (message_text,*)                                                                                      &
          & 'Inconsistency with tlucua bounds (would give a seg. fault or a floating invalid operation exception)', &
          & '- temperature was: ', temp(i), 'K. One thing to check: consistency of land-sea mask and forcing.'
        CALL finish(TRIM(routine), TRIM(message_text))
      ENDIF
    ENDDO

    it = NINT(temp*1000._wp)
    tluc = tlucua(it)

    WHERE (it >= jptlucu1 .AND. it <= jptlucu2)
      qsat = tluc / (pressure - rvd1*tluc)
    ELSEWHERE
      qsat = HUGE(1._wp)
    END WHERE

  END SUBROUTINE sat_specific_humidity

#endif
END MODULE mo_disturb_process
