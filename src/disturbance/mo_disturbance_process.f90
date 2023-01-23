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
    !& q_rel_air_climbuf, fire_rel_hum_threshold,  & ! in
    & q_air,       q_clim_o,      q_clim_m,       & ! in
    & wind_10m,                                   & ! in
    & rain,                                       & ! in
    & swl1,        swl1_clim_o,   swl1_clim_m,    & ! in
    & swl2,        swl2_clim_o,   swl2_clim_m,    & ! in
    & swl3,        swl3_clim_o,   swl3_clim_m,    & ! in
    & swl4,        swl4_clim_o,   swl4_clim_m,    & ! in
    & lai,         lai_clim_o,    lai_clim_m,     & ! in
    !& t_air,        t_clim,                       & ! in
    !& q_rel_air,    q_rel_clim,                   & ! in
    !& q_air,        q_clim,                       & ! in
    !& wind_10m,                                   & ! in
    !& rain,                                       & ! in
    !& swl1,        swl1_clim,                     & ! in
    !& swl2,        swl2_clim,                     & ! in
    !& swl3,        swl3_clim,                     & ! in
    !& swl4,        swl4_clim,                     & ! in
    !& lai,         lai_clim,                      & ! in
    & modified_pfts_fract,                        & ! in
    & fuel,                                       & ! in
    !& fire_litter_threshold,                     & ! in
    !& population_density,                        & ! in
    & lightning_climatology,                      & ! in
    & elevation,                                  & ! in
    & slope,                                      & ! in
    & roughness,                                  & ! in

    & population,                                 & ! in
    & gdp,                                        & ! in
    & hdi,                                        & ! in
    & road_density,                               & ! in
    & land_use_c3ann,                             & ! in
    & land_use_c3nfx,                             & ! in
    & land_use_c3per,                             & ! in
    & land_use_c4ann,                             & ! in
    & land_use_c4per,                             & ! in
    & land_use_pastr,                             & ! in
    & land_use_primf,                             & ! in
    & land_use_primn,                             & ! in
    & land_use_range,                             & ! in
    & land_use_secdf,                             & ! in
    & land_use_secdn,                             & ! in
    & land_use_secma,                             & ! in
    & land_use_secmb,                             & ! in
    & land_use_urban,                             & ! in

    & nc,                                         & ! in
    & tilename,                                   & ! in

    & d_weather_input,     d_weather_hidden,      & ! in
    & d_land_input,        d_land_hidden,         & ! in
    & d_land_input_n_pft,                         & ! in
    !& d_lnt_input,         d_lnt_hidden,         & ! in
    & d_anthro_input,      d_anthro_hidden,       & ! in
    & d_output,            d_region,              & ! in
    !& d_emb_input,         d_emb_hidden,          & ! in
    !& d_emb_output,                               & ! in

    & weather_i2i_weight, weather_i2i_bias,       & ! in
    & weather_h2i_weight, weather_h2i_bias,       & ! in
    & weather_i2f_weight, weather_i2f_bias,       & ! in
    & weather_h2f_weight, weather_h2f_bias,       & ! in
    & weather_i2g_weight, weather_i2g_bias,       & ! in
    & weather_h2g_weight, weather_h2g_bias,       & ! in
    & weather_i2o_weight, weather_i2o_bias,       & ! in
    & weather_h2o_weight, weather_h2o_bias,       & ! in
    & weather_out_weight, weather_out_bias,       & ! in

    & land_i2i_weight,    land_i2i_bias,          & ! in
    & land_h2i_weight,    land_h2i_bias,          & ! in
    & land_i2f_weight,    land_i2f_bias,          & ! in
    & land_h2f_weight,    land_h2f_bias,          & ! in
    & land_i2g_weight,    land_i2g_bias,          & ! in
    & land_h2g_weight,    land_h2g_bias,          & ! in
    & land_i2o_weight,    land_i2o_bias,          & ! in
    & land_h2o_weight,    land_h2o_bias,          & ! in
    & land_out_weight,    land_out_bias,          & ! in

    & anthropogenic_l1_weight, anthropogenic_l1_bias,  & ! in
    & anthropogenic_l2_weight, anthropogenic_l2_bias,  & ! in

    !& lig_topo_layer1_weight,  lig_topo_layer1_bias,  & ! in
    !& lig_topo_layer2_weight,  lig_topo_layer2_bias,  & ! in

    !& climate_emb_layer1_weight,    climate_emb_layer1_bias,      & ! in
    !& climate_emb_layer2_weight,    climate_emb_layer2_bias,      & ! in
    !& vegetation_emb_layer1_weight, vegetation_emb_layer1_bias,   & ! in
    !& vegetation_emb_layer2_weight, vegetation_emb_layer2_bias,   & ! in
 
    & w_minmax,                                   & ! in 
    & l_minmax,                                   & ! in
    & a_minmax,                                   & ! in
    !& embedding_min_max,                         & ! in
    !& lig_topo_min_max,                          & ! in

    & gfed_region,                                & ! in

    & weather_lstm_hidden,                        & ! inout
    & weather_lstm_cell,                          & ! inout
    & land_lstm_hidden,                           & ! inout
    & land_lstm_cell,                             & ! inout
    & burned_fract                                & ! inout
    & )

    ! Input Arguments
    REAL(wp), INTENT(in) ::                             &
      & t_air(:),     t_clim_o(:),     t_clim_m(:),     &
      & q_rel_air(:), q_rel_clim_o(:), q_rel_clim_m(:), &
      !& q_rel_air_climbuf(:), fire_rel_hum_threshold,   &
      & q_air(:),     q_clim_o(:),     q_clim_m(:),     &
      & wind_10m(:),                                    &
      & rain(:),                                        &
      & swl1(:),      swl1_clim_o(:),  swl1_clim_m(:),  &
      & swl2(:),      swl2_clim_o(:),  swl2_clim_m(:),  &
      & swl3(:),      swl3_clim_o(:),  swl3_clim_m(:),  &
      & swl4(:),      swl4_clim_o(:),  swl4_clim_m(:),  &
      & lai(:),       lai_clim_o(:),   lai_clim_m(:),   &
      !& t_air(:),     t_clim(:),                       &
      !& q_rel_air(:), q_rel_clim(:),                   &
      !& q_air(:),     q_clim(:),                       &
      !& wind_10m(:),                                   &
      !& rain(:),                                       &
      !& swl1(:),      swl1_clim(:),                    &
      !& swl2(:),      swl2_clim(:),                    &
      !& swl3(:),      swl3_clim(:),                    &
      !& swl4(:),      swl4_clim(:),                    &
      !& lai(:),       lai_clim(:),                     &
      & modified_pfts_fract(:,:),                       &
      !& fuel(:),       fuel_yearly(:),                 &
      & fuel(:),                                        &
      !& fire_litter_threshold,                         &
      & lightning_climatology(:),                       &
      & elevation(:),                                   &
      & slope(:),                                       &
      & roughness(:),                                   &

      & population(:),                                  &
      & gdp(:),                                         &
      & hdi(:),                                         &
      & road_density(:),                                &
      & land_use_c3ann(:),                              &
      & land_use_c3nfx(:),                              &
      & land_use_c3per(:),                              &
      & land_use_c4ann(:),                              &
      & land_use_c4per(:),                              &
      & land_use_pastr(:),                              &
      & land_use_primf(:),                              &
      & land_use_primn(:),                              &
      & land_use_range(:),                              &
      & land_use_secdf(:),                              &
      & land_use_secdn(:),                              &
      & land_use_secma(:),                              &
      & land_use_secmb(:),                              &
      & land_use_urban(:),                              &

      & weather_i2i_weight(:,:),    weather_i2i_bias(:),      &
      & weather_h2i_weight(:,:),    weather_h2i_bias(:),      &
      & weather_i2f_weight(:,:),    weather_i2f_bias(:),      &
      & weather_h2f_weight(:,:),    weather_h2f_bias(:),      &
      & weather_i2g_weight(:,:),    weather_i2g_bias(:),      &
      & weather_h2g_weight(:,:),    weather_h2g_bias(:),      &
      & weather_i2o_weight(:,:),    weather_i2o_bias(:),      &
      & weather_h2o_weight(:,:),    weather_h2o_bias(:),      &
      & weather_out_weight(:,:),    weather_out_bias(:),      &
      
      & land_i2i_weight(:,:),       land_i2i_bias(:),         &
      & land_h2i_weight(:,:),       land_h2i_bias(:),         &
      & land_i2f_weight(:,:),       land_i2f_bias(:),         &
      & land_h2f_weight(:,:),       land_h2f_bias(:),         &
      & land_i2g_weight(:,:),       land_i2g_bias(:),         &
      & land_h2g_weight(:,:),       land_h2g_bias(:),         &
      & land_i2o_weight(:,:),       land_i2o_bias(:),         &
      & land_h2o_weight(:,:),       land_h2o_bias(:),         &
      & land_out_weight(:,:),       land_out_bias(:),         &

      & anthropogenic_l1_weight(:,:),                         &
      & anthropogenic_l1_bias(:),                             &
      & anthropogenic_l2_weight(:,:),                         &
      & anthropogenic_l2_bias(:),                             &

      !& lig_topo_layer1_weight(:,:),              &
      !& lig_topo_layer1_bias(:),                  &
      !& lig_topo_layer2_weight(:,:),              &
      !& lig_topo_layer2_bias(:),                  &

      !& climate_emb_layer1_weight(:,:), climate_emb_layer1_bias(:),       &
      !& climate_emb_layer2_weight(:,:), climate_emb_layer2_bias(:),       &

      !& vegetation_emb_layer1_weight(:,:), vegetation_emb_layer1_bias(:), &
      !& vegetation_emb_layer2_weight(:,:), vegetation_emb_layer2_bias(:), &

      & w_minmax(:,:),                                       &
      & l_minmax(:,:),                                       &
      & a_minmax(:,:),                                       &
      !& embedding_min_max(:,:),            &
      !& lig_topo_min_max(:,:),             &
      & gfed_region(:)
      

    INTEGER, INTENT(in) ::                          &
      & nc,                                         & ! in
      & d_weather_input,     d_weather_hidden,      & ! in
      & d_land_input,        d_land_hidden,         & ! in
      & d_land_input_n_pft,                         & ! in
      & d_anthro_input,      d_anthro_hidden,       & ! in
      & d_output,            d_region                 ! in
      !& d_lnt_input,         d_lnt_hidden,         & ! in
      !& d_emb_input,         d_emb_hidden,         & ! in
      !& d_emb_output                                 ! in

    CHARACTER(len=*), INTENT(in)::                  &
      & tilename                                      ! in

    ! Output Arguments
    REAL(wp), INTENT(inout) ::                      &
      & weather_lstm_hidden(:,:),                   & ! inout
      & weather_lstm_cell(:,:),                     & ! inout
      & land_lstm_hidden(:,:),                      & ! inout
      & land_lstm_cell(:,:),                        & ! inout
      & burned_fract(:)                               ! inout

    ! Locals
    REAL(wp) ::                                          &
      weather_inputs(nc, d_weather_input),               &
      weather_lstm_i(nc, d_weather_hidden),              &
      weather_lstm_f(nc, d_weather_hidden),              &
      weather_lstm_g(nc, d_weather_hidden),              &
      weather_lstm_o(nc, d_weather_hidden),              &
      weather_lstm_outputs(nc, d_output),                &
      weather_lstm_outputs_sum(nc, d_output),            &
      !weather_lstm_output_(nc, d_output),               &
      !climate_score(nc, d_emb_output),                  &

      !& regional_w_i2i_w(d_weather_input, d_weather_hidden), &
      !& regional_w_i2i_b(d_weather_hidden),                  &
      !& regional_w_h2i_w(d_weather_hidden,d_weather_hidden), &
      !& regional_w_h2i_b(d_weather_hidden),                  &
      !& regional_w_i2f_w(d_weather_input, d_weather_hidden), &
      !& regional_w_i2f_b(d_weather_hidden),                  &
      !& regional_w_h2f_w(d_weather_hidden,d_weather_hidden), &
      !& regional_w_h2f_b(d_weather_hidden),                  &
      !& regional_w_i2g_w(d_weather_input, d_weather_hidden), &
      !& regional_w_i2g_b(d_weather_hidden),                  &
      !& regional_w_h2g_w(d_weather_hidden,d_weather_hidden), &
      !& regional_w_h2g_b(d_weather_hidden),                  &
      !& regional_w_i2o_w(d_weather_input, d_weather_hidden), &
      !& regional_w_i2o_b(d_weather_hidden),                  &
      !& regional_w_h2o_w(d_weather_hidden,d_weather_hidden), &
      !& regional_w_h2o_b(d_weather_hidden),                  &
      !& regional_w_out_w(d_weather_hidden, d_output),        &
      !& regional_w_out_b(d_output),                          &

      land_inputs(nc, d_land_input_n_pft),               &
      land_lstm_i(nc, d_land_hidden),                    &
      land_lstm_f(nc, d_land_hidden),                    &
      land_lstm_g(nc, d_land_hidden),                    &
      land_lstm_o(nc, d_land_hidden),                    &
      land_lstm_outputs(nc, d_output),                   &
      land_lstm_outputs_sum(nc, d_output),               &
      !vegetation_lstm_output_(nc, d_output),            &
      !vegetation_score(nc, d_emb_output),               &

      !& regional_l_i2i_w(d_land_input, d_land_hidden),   &
      !& regional_l_i2i_b(d_land_hidden),                 &
      !& regional_l_h2i_w(d_land_hidden,d_land_hidden),   &
      !& regional_l_h2i_b(d_land_hidden),                 &
      !& regional_l_i2f_w(d_land_input, d_land_hidden),   &
      !& regional_l_i2f_b(d_land_hidden),                 &
      !& regional_l_h2f_w(d_land_hidden,d_land_hidden),   &
      !& regional_l_h2f_b(d_land_hidden),                 &
      !& regional_l_i2g_w(d_land_input, d_land_hidden),   &
      !& regional_l_i2g_b(d_land_hidden),                 &
      !& regional_l_h2g_w(d_land_hidden,d_land_hidden),   &
      !& regional_l_h2g_b(d_land_hidden),                 &
      !& regional_l_i2o_w(d_land_input, d_land_hidden),   &
      !& regional_l_i2o_b(d_land_hidden),                 &
      !& regional_l_h2o_w(d_land_hidden,d_land_hidden),   &
      !& regional_l_h2o_b(d_land_hidden),                 &
      !& regional_l_out_w(d_land_hidden, d_output),       &
      !& regional_l_out_b(d_output),                      &

      anthro_inputs(nc, d_anthro_input),                 &
      anthro_nn(nc, d_anthro_hidden),                    &
      anthro_nn_act(nc, d_anthro_hidden),                &
      anthro_nn_outputs(nc, d_output),                   &
      anthro_nn_outputs_sum(nc, d_output),               &

      & regional_a_l1_w(d_anthro_input, d_anthro_hidden), &
      & regional_a_l1_b(d_anthro_hidden),                 &
      & reginoal_a_l2_w(d_anthro_hidden, d_output),       &
      & regional_a_l2_b(d_output),                        &

      !lig_topo_variables(nc, d_lnt_input),               &
      !lig_topo_network(nc, d_lnt_hidden),                &
      !lig_topo_network_act(nc, d_lnt_hidden),            &
      !lig_topo_network_output(nc, d_output),             &

      !emb_variables(nc, d_emb_input),                    &
      !emb_network(nc, d_emb_hidden),                     &
      !emb_network_act(nc, d_emb_hidden),                 &
      !climate_emb_output(nc, d_emb_output),              &
      !vegetation_emb_output(nc, d_emb_output),           &

      !no_fire(nc),                                       &
      dum1(nc), dum2(nc), dum3(nc), dum4(nc), &
      num1(nc), num2(nc), num3(nc), num4(nc), &
      !fuel_constraint(nc, d_region),                      &
      fuel_constraint(nc),                      &
      q_rel_air_constraint(nc)


    REAL(wp) ::  scale
    INTEGER  ::  pft_idx, r_i, w_i, l_i, a_i, nc_i
    INTEGER  ::  is, ie, hs, he, os, oe


    !climate_variables(:,1) = t_air - t_clim_m + t_clim_o
    !climate_variables(:,2) = t_air - t_clim_m
    !climate_variables(:,3) = q_rel_air - q_rel_clim_m + q_rel_clim_o
    !climate_variables(:,4) = q_rel_air - q_rel_clim_m
    !climate_variables(:,5) = q_air - q_clim_m + q_clim_o
    !climate_variables(:,6) = q_air - q_clim_m
    !climate_variables(:,7) = wind_10m 
    !climate_variables(:,8) = rain*86400

    !vegetation_variables(:,1) = lai  - lai_clim_m  + lai_clim_o
    !vegetation_variables(:,2) = lai  - lai_clim_m
    !vegetation_variables(:,3) = swl1 - swl1_clim_m + swl1_clim_o
    !vegetation_variables(:,4) = swl1 - swl1_clim_m
    !vegetation_variables(:,5) = swl2 - swl2_clim_m + swl2_clim_o
    !vegetation_variables(:,6) = swl2 - swl2_clim_m
    !vegetation_variables(:,7) = swl3 - swl3_clim_m + swl3_clim_o
    !vegetation_variables(:,8) = swl3 - swl3_clim_m
    !vegetation_variables(:,9) = swl4 - swl4_clim_m + swl4_clim_o
    !vegetation_variables(:,10)= swl4 - swl4_clim_m

    !emb_variables(:,1) = t_clim_o
    !emb_variables(:,2) = q_rel_clim_o
    !emb_variables(:,3) = q_clim_o
    !emb_variables(:,4) = swl1_clim_o
    !emb_variables(:,5) = swl2_clim_o
    !emb_variables(:,6) = swl3_clim_o 
    !emb_variables(:,7) = swl4_clim_o
    !emb_variables(:,8) = lai_clim_o
    !emb_variables(:,1) = t_clim_m
    !emb_variables(:,2) = q_rel_clim_m
    !emb_variables(:,3) = q_clim_m
    !emb_variables(:,4) = swl1_clim_m
    !emb_variables(:,5) = swl2_clim_m
    !emb_variables(:,6) = swl3_clim_m 
    !emb_variables(:,7) = swl4_clim_m
    !emb_variables(:,8) = lai_clim_m
    !emb_variables(:,9) = fuel
    !emb_variables(:,10)= lightning_climatology
    !emb_variables(:,11)= elevation
    !emb_variables(:,12)= slope
    !emb_variables(:,13)= roughness
    !emb_variables(:,14)= modified_pfts_fract(:,1)
    !emb_variables(:,15)= modified_pfts_fract(:,2)
    !emb_variables(:,16)= modified_pfts_fract(:,3)
    !emb_variables(:,17)= modified_pfts_fract(:,4)
    !emb_variables(:,18)= modified_pfts_fract(:,5)
    !emb_variables(:,19)= modified_pfts_fract(:,6)
    !emb_variables(:,20)= modified_pfts_fract(:,7)
    !emb_variables(:,21)= modified_pfts_fract(:,8)
    !emb_variables(:,22)= modified_pfts_fract(:,9)

    !anthro_variables(:,1) = population_density
    !anthro_variables(:,2) = gdp
    !anthro_variables(:,3) = hdi
    !anthro_variables(:,4) = road_density
    !anthro_variables(:,5) = land_use_state_c3ann
    !anthro_variables(:,6) = land_use_state_c3nfx
    !anthro_variables(:,7) = land_use_state_c3per
    !anthro_variables(:,8) = land_use_state_c4ann
    !anthro_variables(:,9) = land_use_state_c4per
    !anthro_variables(:,10)= land_use_state_pastr
    !anthro_variables(:,11)= land_use_state_primf
    !anthro_variables(:,12)= land_use_state_primn
    !anthro_variables(:,13)= land_use_state_range
    !anthro_variables(:,14)= land_use_state_secdf
    !anthro_variables(:,15)= land_use_state_secdn
    !anthro_variables(:,16)= land_use_state_secma
    !anthro_variables(:,17)= land_use_state_secmb
    !anthro_variables(:,18)= land_use_state_urban

    ! initialize
    weather_lstm_outputs_sum(:,:) = 0.
    land_lstm_outputs_sum(:,:)    = 0.
    anthro_nn_outputs_sum(:,:)    = 0.
    !fuel_constraint(:)          = 0.
    fuel_constraint(:)            = 0.

    DO r_i=1,d_region
      w_i = (r_i-1)*d_weather_input
      l_i = (r_i-1)*d_land_input
      a_i = (r_i-1)*d_anthro_input

      weather_inputs(:,:) = 0.
      land_inputs(:,:)    = 0.
      anthro_inputs(:,:)  = 0.

    DO nc_i=1,nc
      ! reginally normalized input
      !print*, gfed_region(nc_i), int(gfed_region(nc_i))

      IF (r_i .eq. int(gfed_region(nc_i))) THEN
        weather_inputs(nc_i,1) = (t_air(nc_i)                   - w_minmax(1+w_i,1)) &
                               / (w_minmax(1+w_i,2)             - w_minmax(1+w_i,1)+1e-8)
        weather_inputs(nc_i,2) = (t_air(nc_i)-t_clim_m(nc_i)    - w_minmax(2+w_i,1)) &
                               / (w_minmax(2+w_i,2)             - w_minmax(2+w_i,1)+1e-8)
        weather_inputs(nc_i,3) = (q_rel_air(nc_i)               - w_minmax(3+w_i,1)) &
                               / (w_minmax(3+w_i,2)             - w_minmax(3+w_i,1)+1e-8)
        weather_inputs(nc_i,4) = (q_rel_air(nc_i)-q_rel_clim_m(nc_i) - w_minmax(4+w_i,1)) &
                               / (w_minmax(4+w_i,2)                  - w_minmax(4+w_i,1)+1e-8)
        weather_inputs(nc_i,5) = (q_air(nc_i)                   - w_minmax(5+w_i,1)) &
                               / (w_minmax(5+w_i,2)             - w_minmax(5+w_i,1)+1e-8)
        weather_inputs(nc_i,6) = (q_air(nc_i)-q_clim_m(nc_i)    - w_minmax(6+w_i,1)) &
                               / (w_minmax(6+w_i,2)             - w_minmax(6+w_i,1)+1e-8)
        weather_inputs(nc_i,7) = (wind_10m(nc_i)              - w_minmax(7+w_i,1)) &
                               / (w_minmax(7+w_i,2)           - w_minmax(7+w_i,1)+1e-8)
        weather_inputs(nc_i,8) = (rain(nc_i)*86400            - w_minmax(8+w_i,1)) &
                               / (w_minmax(8+w_i,2)           - w_minmax(8+w_i,1)+1e-8)
        weather_inputs(nc_i,9) = (lightning_climatology(nc_i) - w_minmax(9+w_i,1)) &
                               / (w_minmax(9+w_i,2)           - w_minmax(9+w_i,1)+1e-8)

        land_inputs(nc_i,1) = (lai(nc_i)                    - l_minmax(1+l_i,1))  &
                            / (l_minmax(1+l_i,2)            - l_minmax(1+l_i,1)+1e-8)
        land_inputs(nc_i,2) = (lai(nc_i)-lai_clim_m(nc_i)   - l_minmax(2+l_i,1))  &
                            / (l_minmax(2+l_i,2)            - l_minmax(2+l_i,1)+1e-8)
        land_inputs(nc_i,3) = (swl1(nc_i)                   - l_minmax(3+l_i,1))  &
                            / (l_minmax(3+l_i,2)            - l_minmax(3+l_i,1)+1e-8)
        land_inputs(nc_i,4) = (swl1(nc_i)-swl1_clim_m(nc_i) - l_minmax(4+l_i,1))  &
                            / (l_minmax(4+l_i,2)            - l_minmax(4+l_i,1)+1e-8)
        land_inputs(nc_i,5) = (swl2(nc_i)                   - l_minmax(5+l_i,1))  &
                            / (l_minmax(5+l_i,2)            - l_minmax(5+l_i,1)+1e-8)
        land_inputs(nc_i,6) = (swl2(nc_i)-swl2_clim_m(nc_i) - l_minmax(6+l_i,1))  &
                            / (l_minmax(6+l_i,2)            - l_minmax(6+l_i,1)+1e-8)
        land_inputs(nc_i,7) = (swl3(nc_i)                   - l_minmax(7+l_i,1))  &
                            / (l_minmax(7+l_i,2)            - l_minmax(7+l_i,1)+1e-8)
        land_inputs(nc_i,8) = (swl3(nc_i)-swl3_clim_m(nc_i) - l_minmax(8+l_i,1))  &
                            / (l_minmax(8+l_i,2)            - l_minmax(8+l_i,1)+1e-8)
        land_inputs(nc_i,9) = (swl4(nc_i)                   - l_minmax(9+l_i,1))  &
                            / (l_minmax(9+l_i,2)            - l_minmax(9+l_i,1)+1e-8)
        land_inputs(nc_i,10)= (swl4(nc_i)-swl4_clim_m(nc_i) - l_minmax(10+l_i,1)) &
                            / (l_minmax(10+l_i,2)           - l_minmax(10+l_i,1)+1e-8) 
        land_inputs(nc_i,11)= (fuel(nc_i)                 - l_minmax(11+l_i,1)) &
                            / (l_minmax(11+l_i,2)         - l_minmax(11+l_i,1)+1e-8)
        land_inputs(nc_i,12)= (elevation(nc_i)            - l_minmax(12+l_i,1)) &
                            / (l_minmax(12+l_i,2)         - l_minmax(12+l_i,1)+1e-8)
        land_inputs(nc_i,13)= (slope(nc_i)                - l_minmax(13+l_i,1)) &
                            / (l_minmax(13+l_i,2)         - l_minmax(13+l_i,1)+1e-8)
        land_inputs(nc_i,14)= (roughness(nc_i)            - l_minmax(14+l_i,1)) &
                            / (l_minmax(14+l_i,2)         - l_minmax(14+l_i,1)+1e-8)
        
        land_inputs(nc_i,15)= modified_pfts_fract(nc_i,1)
        land_inputs(nc_i,16)= modified_pfts_fract(nc_i,2)
        land_inputs(nc_i,17)= modified_pfts_fract(nc_i,3)
        land_inputs(nc_i,18)= modified_pfts_fract(nc_i,4)
        land_inputs(nc_i,19)= modified_pfts_fract(nc_i,5)
        land_inputs(nc_i,20)= modified_pfts_fract(nc_i,6)
        land_inputs(nc_i,21)= modified_pfts_fract(nc_i,7)
        land_inputs(nc_i,22)= modified_pfts_fract(nc_i,8)
        land_inputs(nc_i,23)= modified_pfts_fract(nc_i,9)

        anthro_inputs(nc_i,1) = (LOG(population(nc_i)+1)    - a_minmax(1+a_i,1))  &
                              / (a_minmax(1+a_i,2)          - a_minmax(1+a_i,1)+1e-8)
        anthro_inputs(nc_i,2) = (gdp(nc_i)                  - a_minmax(2+a_i,1))  &
                              / (a_minmax(2+a_i,2)          - a_minmax(2+a_i,1)+1e-8)
        anthro_inputs(nc_i,3) = (hdi(nc_i)                  - a_minmax(3+a_i,1))  &
                              / (a_minmax(3+a_i,2)          - a_minmax(3+a_i,1)+1e-8)
        anthro_inputs(nc_i,4) = (LOG(road_density(nc_i)+1)  - a_minmax(4+a_i,1))  &
                              / (a_minmax(4+a_i,2)          - a_minmax(4+a_i,1)+1e-8)
        anthro_inputs(nc_i,5) = (land_use_c3ann(nc_i)       - a_minmax(5+a_i,1))  &
                              / (a_minmax(5+a_i,2)          - a_minmax(5+a_i,1)+1e-8)
        anthro_inputs(nc_i,6) = (land_use_c3nfx(nc_i)       - a_minmax(6+a_i,1))  &
                              / (a_minmax(6+a_i,2)          - a_minmax(6+a_i,1)+1e-8)
        anthro_inputs(nc_i,7) = (land_use_c3per(nc_i)       - a_minmax(7+a_i,1))  &
                              / (a_minmax(7+a_i,2)          - a_minmax(7+a_i,1)+1e-8)
        anthro_inputs(nc_i,8) = (land_use_c4ann(nc_i)       - a_minmax(8+a_i,1))  &
                              / (a_minmax(8+a_i,2)          - a_minmax(8+a_i,1)+1e-8)
        anthro_inputs(nc_i,9) = (land_use_c4per(nc_i)       - a_minmax(9+a_i,1))  &
                              / (a_minmax(9+a_i,2)          - a_minmax(9+a_i,1)+1e-8)
        anthro_inputs(nc_i,10)= (land_use_pastr(nc_i)       - a_minmax(10+a_i,1)) &
                              / (a_minmax(10+a_i,2)         - a_minmax(10+a_i,1)+1e-8)
        anthro_inputs(nc_i,11)= (land_use_primf(nc_i)       - a_minmax(11+a_i,1)) &
                              / (a_minmax(11+a_i,2)         - a_minmax(11+a_i,1)+1e-8)
        anthro_inputs(nc_i,12)= (land_use_primn(nc_i)       - a_minmax(12+a_i,1)) &
                              / (a_minmax(12+a_i,2)         - a_minmax(12+a_i,1)+1e-8)
        anthro_inputs(nc_i,13)= (land_use_range(nc_i)       - a_minmax(13+a_i,1)) &
                              / (a_minmax(13+a_i,2)         - a_minmax(13+a_i,1)+1e-8)
        anthro_inputs(nc_i,14)= (land_use_secdf(nc_i)       - a_minmax(14+a_i,1)) &
                              / (a_minmax(14+a_i,2)         - a_minmax(14+a_i,1)+1e-8)
        anthro_inputs(nc_i,15)= (land_use_secdn(nc_i)       - a_minmax(15+a_i,1)) &
                              / (a_minmax(15+a_i,2)         - a_minmax(15+a_i,1)+1e-8)
        anthro_inputs(nc_i,16)= (land_use_secma(nc_i)       - a_minmax(16+a_i,1)) &
                              / (a_minmax(16+a_i,2)         - a_minmax(16+a_i,1)+1e-8)
        anthro_inputs(nc_i,17)= (land_use_secmb(nc_i)       - a_minmax(17+a_i,1)) &
                              / (a_minmax(17+a_i,2)         - a_minmax(17+a_i,1)+1e-8)
        anthro_inputs(nc_i,18)= (land_use_urban(nc_i)       - a_minmax(18+a_i,1)) &
                              / (a_minmax(18+a_i,2)         - a_minmax(18+a_i,1)+1e-8)

        !fuel_constraint(nc_i, r_i) = land_inputs(nc_i,11)
        fuel_constraint(nc_i) = land_inputs(nc_i,11)
        
        !if (isnan(land_inputs(nc_i,11))) then
        !  print*, r_i
        !endif

        where (isnan(anthro_inputs)) anthro_inputs = 0.

      ENDIF
    ENDDO   ! nc iteration

      os=1+(r_i-1)*d_output
      oe=r_i*d_output

      ! weather-LSTM model
      hs=1+(r_i-1)*d_weather_hidden
      he=r_i*d_weather_hidden
      
      weather_lstm_i = matmul(weather_inputs, weather_i2i_weight(:,hs:he)) + spread(weather_i2i_bias(hs:he),1,nc) &
       + matmul(weather_lstm_hidden(:,hs:he), weather_h2i_weight(:,hs:he)) + spread(weather_h2i_bias(hs:he),1,nc)
      weather_lstm_f = matmul(weather_inputs, weather_i2f_weight(:,hs:he)) + spread(weather_i2f_bias(hs:he),1,nc) &
       + matmul(weather_lstm_hidden(:,hs:he), weather_h2f_weight(:,hs:he)) + spread(weather_h2f_bias(hs:he),1,nc)
      weather_lstm_g = matmul(weather_inputs, weather_i2g_weight(:,hs:he)) + spread(weather_i2g_bias(hs:he),1,nc) &
       + matmul(weather_lstm_hidden(:,hs:he), weather_h2g_weight(:,hs:he)) + spread(weather_h2g_bias(hs:he),1,nc)
      weather_lstm_o = matmul(weather_inputs, weather_i2o_weight(:,hs:he)) + spread(weather_i2o_bias(hs:he),1,nc) &
       + matmul(weather_lstm_hidden(:,hs:he), weather_h2o_weight(:,hs:he)) + spread(weather_h2o_bias(hs:he),1,nc)

      weather_lstm_i = 1/(1+exp(-weather_lstm_i))
      weather_lstm_f = 1/(1+exp(-weather_lstm_f))
      weather_lstm_g = tanh(weather_lstm_g)
      weather_lstm_o = 1/(1+exp(-weather_lstm_o))

      weather_lstm_cell(:,hs:he)   = weather_lstm_f * weather_lstm_cell(:,hs:he) &
                                   + weather_lstm_i * weather_lstm_g
      weather_lstm_hidden(:,hs:he) = weather_lstm_o * tanh(weather_lstm_cell(:,hs:he))

      weather_lstm_outputs = matmul(weather_lstm_hidden(:,hs:he), weather_out_weight(:,os:oe)) &
                           + spread(weather_out_bias(os:oe),1,nc)

      ! make 0 for other regions
      where (spread(gfed_region, 2, d_output) /= float(r_i)) weather_lstm_outputs = 0.
      weather_lstm_outputs_sum(:,:) = weather_lstm_outputs_sum(:,:) + weather_lstm_outputs(:,:)


      ! land-LSTM model
      hs=1+(r_i-1)*d_land_hidden
      he=r_i*d_land_hidden
      
      land_lstm_i =    matmul(land_inputs, land_i2i_weight(:,hs:he)) + spread(land_i2i_bias(hs:he),1,nc) &
       + matmul(land_lstm_hidden(:,hs:he), land_h2i_weight(:,hs:he)) + spread(land_h2i_bias(hs:he),1,nc)
      land_lstm_f =    matmul(land_inputs, land_i2f_weight(:,hs:he)) + spread(land_i2f_bias(hs:he),1,nc) &
       + matmul(land_lstm_hidden(:,hs:he), land_h2f_weight(:,hs:he)) + spread(land_h2f_bias(hs:he),1,nc)
      land_lstm_g =    matmul(land_inputs, land_i2g_weight(:,hs:he)) + spread(land_i2g_bias(hs:he),1,nc) &
       + matmul(land_lstm_hidden(:,hs:he), land_h2g_weight(:,hs:he)) + spread(land_h2g_bias(hs:he),1,nc)
      land_lstm_o =    matmul(land_inputs, land_i2o_weight(:,hs:he)) + spread(land_i2o_bias(hs:he),1,nc) &
       + matmul(land_lstm_hidden(:,hs:he), land_h2o_weight(:,hs:he)) + spread(land_h2o_bias(hs:he),1,nc)


      land_lstm_i = 1/(1+exp(-land_lstm_i))
      land_lstm_f = 1/(1+exp(-land_lstm_f))
      land_lstm_g = tanh(land_lstm_g)
      land_lstm_o = 1/(1+exp(-land_lstm_o))

      land_lstm_cell(:,hs:he)   = land_lstm_f * land_lstm_cell(:,hs:he) &
                                + land_lstm_i * land_lstm_g
      land_lstm_hidden(:,hs:he) = land_lstm_o * tanh(land_lstm_cell(:,hs:he))

      land_lstm_outputs = matmul(land_lstm_hidden(:,hs:he), land_out_weight(:,os:oe)) &
                        + spread(land_out_bias(os:oe),1,nc)

      ! make 0 for other regions
      where (spread(gfed_region, 2, d_output) /= float(r_i)) land_lstm_outputs = 0.
      land_lstm_outputs_sum(:,:) = land_lstm_outputs_sum(:,:) + land_lstm_outputs(:,:)


      ! anthropogenic NN model
      hs=1+(r_i-1)*d_anthro_hidden
      he=r_i*d_anthro_hidden

      anthro_nn = matmul(anthro_inputs, anthropogenic_l1_weight(:,hs:he)) &
                + spread(anthropogenic_l1_bias(hs:he),1,nc)
      anthro_nn_act = log(1+exp(anthro_nn))
      where (anthro_nn > 20) anthro_nn_act = anthro_nn
      anthro_nn_outputs = matmul(anthro_nn_act, anthropogenic_l2_weight(:,os:oe)) &
                        + spread(anthropogenic_l2_bias(os:oe),1,nc)

      ! make 0 for other regions
      where (spread(gfed_region, 2, d_output) /= float(r_i)) anthro_nn_outputs = 0.
      anthro_nn_outputs_sum(:,:) = anthro_nn_outputs_sum(:,:) + anthro_nn_outputs(:,:)

    ENDDO   ! r (region) iteration

    !dum1 = sum(weather_lstm_outputs_sum, dim=2)
    !dum2 = sum(land_lstm_outputs_sum, dim=2)
    !dum3 = sum(anthro_nn_outputs_sum, dim=2)
    !dum4 = sum(fuel_constraint, dim=2)
    !dum4 = fuel_constraint

    !where(gfed_region /= 4) dum1 = 0.
    !where(gfed_region /= 4) dum2 = 0.
    !where(gfed_region /= 4) dum3 = 0.
    !where(gfed_region /= 4) dum4 = 0.

    !num1 = gfed_region
    !num2 = gfed_region
    !num3 = gfed_region
    !num4 = gfed_region

    !where(gfed_region /= 4) num1 = 0.
    !where(gfed_region /= 4) num2 = 0.
    !where(gfed_region /= 4) num3 = 0.
    !where(gfed_region /= 4) num4 = 0.
    !where(gfed_region == 4) num1 = 1.
    !where(gfed_region == 4) num2 = 1.
    !where(gfed_region == 4) num3 = 1.
    !where(gfed_region == 4) num4 = 1.

    !if (ISNAN(sum(dum4))) then
    !  print*, fuel
    !endif

    !if (sum(num4) .ne. 0) then
    !  print*, sum(dum1)/(sum(num1)+1e-7), sum(dum2)/(sum(num2)+1e-7), sum(dum3)/(sum(num3)+1e-7), sum(dum4)/(sum(num4)+1e-7)
    !endif
    !if (sum(fuel_constraint) /= 0.) then
    !  print*, sum(fuel_constraint)/(sum(num4)+1e-7)
    !endif

    burned_fract = sum(weather_lstm_outputs_sum     &
                     * land_lstm_outputs_sum        &
                     * anthro_nn_outputs_sum        &
                     * modified_pfts_fract(:,1:8), dim=2)

    burned_fract(:) = burned_fract(:)/100.
    burned_fract(:) = burned_fract(:)*(1.-modified_pfts_fract(:,9)-modified_pfts_fract(:,1))

    !burned_fract(:) = burned_fract(:)*sum(fuel_constraint, dim=2)
    burned_fract(:) = burned_fract(:)*fuel_constraint

    !if (maxval(fuel_constraint) .gt. 0.1) then
    !  print*, minval(fuel_constraint), maxval(fuel_constraint), maxval(fuel)
    !endif

    q_rel_air_constraint(:) = q_rel_air(:)/100.
    where (q_rel_air_constraint < 0) q_rel_air_constraint = 0.
    where (q_rel_air_constraint > 1) q_rel_air_constraint = 1.

    !q_rel_air_constraint(:) = (q_rel_air(:)-q_rel_clim_m(:)+q_rel_clim_o(:))/100.
    q_rel_air_constraint(:) = 1/(1+exp(-20*(q_rel_air_constraint(:)-0.5)))
    !where (q_rel_air > 50) q_rel_air_constraint = 1.
    where (q_rel_air > 60) q_rel_air_constraint = 1.
    !where (q_rel_air > 70) q_rel_air_constraint = 1.
    burned_fract(:) = burned_fract(:)*(1-q_rel_air_constraint(:))

    !dum1 = burned_fract
    !where(gfed_region /= 4) dum1 = 0.
    !if (sum(num1) .ne. 0) then
    !  print*, sum(dum1)/(sum(num1)+1e-8) *100*365
    !endif

    where (burned_fract < 0.) burned_fract = 0.
    where (burned_fract > 1.) burned_fract = 1.


    
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
