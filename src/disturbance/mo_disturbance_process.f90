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
    & t_air,       t_clim_m,                & ! in
    & q_rel_air,   q_rel_clim_m,            & ! in
    & q_air,       q_clim_m,                & ! in
    & wind_10m,                             & ! in
    & rain,                                 & ! in
    & swl1,        swl1_clim_m,             & ! in
    & swl2,        swl2_clim_m,             & ! in
    & swl3,        swl3_clim_m,             & ! in
    & swl4,        swl4_clim_m,             & ! in
    & swl5,        swl5_clim_m,             & ! in
    & lai,         lai_clim_m,              & ! in
    !& modified_pfts_fract,                  & ! in
    & pfts_fract,                           & ! in
    & fuel,                                 & ! in
    & lightning_climatology,                & ! in
    & elevation,                            & ! in
    !& slope,                                & ! in
    !& roughness,                            & ! in

    & population,                           & ! in
    & gdp,                                  & ! in
    & hdi,                                  & ! in
    !& road_density,                         & ! in
    !& land_use_c3ann,                       & ! in
    !& land_use_c3nfx,                       & ! in
    !& land_use_c3per,                       & ! in
    !& land_use_c4ann,                       & ! in
    !& land_use_c4per,                       & ! in
    !& land_use_pastr,                       & ! in
    !& land_use_primf,                       & ! in
    !& land_use_primn,                       & ! in
    !& land_use_range,                       & ! in
    !& land_use_secdf,                       & ! in
    !& land_use_secdn,                       & ! in
    !& land_use_secma,                       & ! in
    !& land_use_secmb,                       & ! in
    !& land_use_urban,                       & ! in

    & fba_smooth,                           & ! in

    & nc,                                   & ! in
    & tilename,                             & ! in

    & d_input,     d_hidden,                & ! in
    & d_output,    d_region,                & ! in

    & layer1_weight, layer1_bias,            & ! in
    & layer2_weight, layer2_bias,            & ! in
    & layer3_weight, layer3_bias,            & ! in
    & layer4_weight, layer4_bias,            & ! in
    & layer5_weight, layer5_bias,            & ! in
    & layer6_weight, layer6_bias,            & ! in

    & layer1_weight_update, layer1_bias_update,            & ! in
    & layer2_weight_update, layer2_bias_update,            & ! in
    & layer3_weight_update, layer3_bias_update,            & ! in
    & layer4_weight_update, layer4_bias_update,            & ! in
    & layer5_weight_update, layer5_bias_update,            & ! in
    & layer6_weight_update, layer6_bias_update,            & ! in

    & mW1,           mb1,                    & ! in
    & mW2,           mb2,                    & ! in
    & mW3,           mb3,                    & ! in
    & mW4,           mb4,                    & ! in
    & mW5,           mb5,                    & ! in
    & mW6,           mb6,                    & ! in

    & vW1,           vb1,                    & ! in
    & vW2,           vb2,                    & ! in
    & vW3,           vb3,                    & ! in
    & vW4,           vb4,                    & ! in
    & vW5,           vb5,                    & ! in
    & vW6,           vb6,                    & ! in

    & iteration,                             & ! in
    & is_training,                           & ! in
    & is_saving,                             & ! in
    & is_infer,                              & ! in
    & memory_idx,                            & ! in

    & V1_memory,  V2_memory,  V3_memory,     & ! inout
    & V4_memory,  V5_memory,  V6_memory,     & ! inout
    & V7_memory,  V8_memory,  V9_memory,     & ! inout
    & V10_memory, V11_memory, V12_memory,    & ! inout
    & V13_memory, V14_memory, V15_memory,    & ! inout
    & V16_memory, V17_memory, V18_memory,    & ! inout
    & V19_memory, V20_memory, V21_memory,    & ! inout
    & V22_memory, V23_memory, V24_memory,    & ! inout
    & V25_memory, V26_memory, V27_memory,    & ! inout
    & V28_memory, V29_memory, V30_memory,    & ! inout
    & V31_memory, V32_memory, V33_memory,    & ! inout
    & V34_memory, V35_memory, V36_memory,    & ! inout
    & V37_memory, V38_memory, V39_memory,    & ! inout
    & RHC_memory, Y_memory,                  & ! inout

    & minmax,                               & ! in 
    & gfed_region,                          & ! in

    & burned_fract,                         & ! inout
    & burned_fract_infer,                   & ! inout
    & burned_fract_test                     & ! inout
    & )

    !USE mo_global_variable

    ! Input Arguments
    REAL(wp), INTENT(in) ::                         &
      & t_air(:),          t_clim_m(:),             &
      & q_rel_air(:),  q_rel_clim_m(:),             &
      & q_air(:),          q_clim_m(:),             &
      & wind_10m(:),                                &
      & rain(:),                                    &
      & swl1(:),        swl1_clim_m(:),             &
      & swl2(:),        swl2_clim_m(:),             &
      & swl3(:),        swl3_clim_m(:),             &
      & swl4(:),        swl4_clim_m(:),             &
      & swl5(:),        swl5_clim_m(:),             &
      & lai(:),         lai_clim_m(:),              &

      !& modified_pfts_fract(:,:),                   &
      & pfts_fract(:,:),                            &
      & fuel(:),                                    &
      & lightning_climatology(:),                   &
      & elevation(:),                               &
      !& slope(:),                                   &
      !& roughness(:),                               &

      & population(:),                              &
      & gdp(:),                                     &
      & hdi(:),                                     &
      !& road_density(:),                            &
      !& land_use_c3ann(:),                          &
      !& land_use_c3nfx(:),                          &
      !& land_use_c3per(:),                          &
      !& land_use_c4ann(:),                          &
      !& land_use_c4per(:),                          &
      !& land_use_pastr(:),                          &
      !& land_use_primf(:),                          &
      !& land_use_primn(:),                          &
      !& land_use_range(:),                          &
      !& land_use_secdf(:),                          &
      !& land_use_secdn(:),                          &
      !& land_use_secma(:),                          &
      !& land_use_secmb(:),                          &
      !& land_use_urban(:),                          &

      & layer1_weight(:,:),    layer1_bias(:),      &
      & layer2_weight(:,:),    layer2_bias(:),      &
      & layer3_weight(:,:),    layer3_bias(:),      &
      & layer4_weight(:,:),    layer4_bias(:),      &
      & layer5_weight(:,:),    layer5_bias(:),      &
      & layer6_weight(:,:),    layer6_bias(:),      &

      & fba_smooth(:),                              &
      & minmax(:,:),                                &
      & gfed_region(:)
      

    INTEGER, INTENT(in) ::                          &
      & nc,                                         & ! in
      & d_input,             d_hidden,              & ! in
      & d_output,            d_region,              & ! in
      & iteration,           is_training,           & ! in
      & is_saving,           is_infer,              & ! in
      & memory_idx                                    ! in

    CHARACTER(len=*), INTENT(in)::                  &
      & tilename                                      ! in

    ! Output Arguments
    REAL(wp), INTENT(inout) ::                      &
      & layer1_weight_update(:,:),    layer1_bias_update(:),      &  ! inout (training)
      & layer2_weight_update(:,:),    layer2_bias_update(:),      &  ! inout (training)
      & layer3_weight_update(:,:),    layer3_bias_update(:),      &  ! inout (training)
      & layer4_weight_update(:,:),    layer4_bias_update(:),      &  ! inout (training)
      & layer5_weight_update(:,:),    layer5_bias_update(:),      &  ! inout (training)
      & layer6_weight_update(:,:),    layer6_bias_update(:),      &  ! inout (training)

      & mW1(:,:),              mb1(:),              &  ! inout (training)
      & mW2(:,:),              mb2(:),              &  ! inout (training)
      & mW3(:,:),              mb3(:),              &  ! inout (training)
      & mW4(:,:),              mb4(:),              &  ! inout (training)
      & mW5(:,:),              mb5(:),              &  ! inout (training)
      & mW6(:,:),              mb6(:),              &  ! inout (training)

      & vW1(:,:),              vb1(:),              &  ! inout (training)
      & vW2(:,:),              vb2(:),              &  ! inout (training)
      & vW3(:,:),              vb3(:),              &  ! inout (training)
      & vW4(:,:),              vb4(:),              &  ! inout (training)
      & vW5(:,:),              vb5(:),              &  ! inout (training)
      & vW6(:,:),              vb6(:),              &  ! inout (training)

      & V1_memory(:,:),  V2_memory(:,:),  V3_memory(:,:),     & ! inout (training)
      & V4_memory(:,:),  V5_memory(:,:),  V6_memory(:,:),     & ! inout (training)
      & V7_memory(:,:),  V8_memory(:,:),  V9_memory(:,:),     & ! inout (training)
      & V10_memory(:,:), V11_memory(:,:), V12_memory(:,:),    & ! inout (training)
      & V13_memory(:,:), V14_memory(:,:), V15_memory(:,:),    & ! inout (training)
      & V16_memory(:,:), V17_memory(:,:), V18_memory(:,:),    & ! inout (training)
      & V19_memory(:,:), V20_memory(:,:), V21_memory(:,:),    & ! inout (training)
      & V22_memory(:,:), V23_memory(:,:), V24_memory(:,:),    & ! inout (training)
      & V25_memory(:,:), V26_memory(:,:), V27_memory(:,:),    & ! inout (training)
      & V28_memory(:,:), V29_memory(:,:), V30_memory(:,:),    & ! inout (training)
      & V31_memory(:,:), V32_memory(:,:), V33_memory(:,:),    & ! inout (training)
      & V34_memory(:,:), V35_memory(:,:), V36_memory(:,:),    & ! inout (training)
      & V37_memory(:,:), V38_memory(:,:), V39_memory(:,:),    & ! inout (training)
      & RHC_memory(:,:), Y_memory(:,:),                       & ! inout (training)

      & burned_fract(:),                                      & ! inout
      & burned_fract_infer(:),                                & ! inout
      & burned_fract_test(:)                                    ! inout

    ! Locals
    INTEGER, PARAMETER :: bs = 1024
    !INTEGER, PARAMETER :: bs = 512
    REAL(wp) ::                                                   &
      & nn_inputs    (nc, d_input),                               &
      & nn_layer1_out(nc, d_hidden), nn_layer1_act(nc, d_hidden), &
      & nn_layer2_out(nc, d_hidden), nn_layer2_act(nc, d_hidden), &
      & nn_layer3_out(nc, d_hidden), nn_layer3_act(nc, d_hidden), &
      & nn_layer4_out(nc, d_hidden), nn_layer4_act(nc, d_hidden), &
      & nn_layer5_out(nc, d_hidden), nn_layer5_act(nc, d_hidden), &
      & nn_layer6_out(nc, d_output), nn_layer6_act(nc, d_output), &
      & nn_outputs   (nc, d_output),                              &

      & nn_layer1_mask(nc, d_hidden),                             &  ! dropout
      & nn_layer2_mask(nc, d_hidden),                             &
      & nn_layer3_mask(nc, d_hidden),                             &
      & nn_layer4_mask(nc, d_hidden),                             &
      & nn_layer5_mask(nc, d_hidden),                             &
      & nn_layer6_mask(nc, d_output),                             &

      & tr_inputs    (bs, d_input),                               &  ! <training>
      & tr_layer1_out(bs, d_hidden), tr_layer1_act(bs, d_hidden), &
      & tr_layer2_out(bs, d_hidden), tr_layer2_act(bs, d_hidden), &
      & tr_layer3_out(bs, d_hidden), tr_layer3_act(bs, d_hidden), &
      & tr_layer4_out(bs, d_hidden), tr_layer4_act(bs, d_hidden), &
      & tr_layer5_out(bs, d_hidden), tr_layer5_act(bs, d_hidden), &
      & tr_layer6_out(bs, d_output), tr_layer6_act(bs, d_output), &
      & tr_outputs   (bs, d_output),                              &

      & tr_layer1_mask(bs, d_hidden),                             &  ! <training> dropout
      & tr_layer2_mask(bs, d_hidden),                             &
      & tr_layer3_mask(bs, d_hidden),                             &
      & tr_layer4_mask(bs, d_hidden),                             &
      & tr_layer5_mask(bs, d_hidden),                             &
      & tr_layer6_mask(bs, d_output),                             &

      & dW1(d_input,  d_hidden),     db1(d_hidden),               &  ! derivatives
      & dW2(d_hidden, d_hidden),     db2(d_hidden),               &
      & dW3(d_hidden, d_hidden),     db3(d_hidden),               &
      & dW4(d_hidden, d_hidden),     db4(d_hidden),               &
      & dW5(d_hidden, d_hidden),     db5(d_hidden),               &
      & dW6(d_hidden, d_output),     db6(d_output)

                 
    INTEGER  :: i, j, k, idx
    real(wp) :: lr_rate, keep_prob, loss, u
    real(wp) :: q_rel_air_constraint(nc)
    real(wp) :: tr_q_rel_air_constraint(bs), tr_target(bs)

    ! initialize
    nn_inputs(:,1) = (t_air(:)                     - minmax(1,1))       &
                   / (minmax(1,2)                  - minmax(1,1)+1e-8)
    nn_inputs(:,2) = (t_air(:)-t_clim_m(:)         - minmax(2,1))       &
                   / (minmax(2,2)                  - minmax(2,1)+1e-8)
    nn_inputs(:,3) = (q_rel_air(:)                 - minmax(3,1))       &
                   / (minmax(3,2)                  - minmax(3,1)+1e-8)
    nn_inputs(:,4) = (q_rel_air(:)-q_rel_clim_m(:) - minmax(4,1))       &
                   / (minmax(4,2)                  - minmax(4,1)+1e-8)
    nn_inputs(:,5) = (q_air(:)                     - minmax(5,1))       &
                   / (minmax(5,2)                  - minmax(5,1)+1e-8)
    nn_inputs(:,6) = (q_air(:)-q_clim_m(:)         - minmax(6,1))       &
                   / (minmax(6,2)                  - minmax(6,1)+1e-8)
    nn_inputs(:,7) = (wind_10m(:)                  - minmax(7,1))       &
                   / (minmax(7,2)                  - minmax(7,1)+1e-8)
    nn_inputs(:,8) = (rain(:)                      - minmax(8,1))       &
                   / (minmax(8,2)                  - minmax(8,1)+1e-8)
    nn_inputs(:,9) = (swl1(:)                      - minmax(9,1))       &
                   / (minmax(9,2)                  - minmax(9,1)+1e-8)
    nn_inputs(:,10)= (swl1(:)-swl1_clim_m(:)       - minmax(10,1))      &
                   / (minmax(10,2)                 - minmax(10,1)+1e-8)
    nn_inputs(:,11)= (swl2(:)                      - minmax(11,1))      &
                   / (minmax(11,2)                 - minmax(11,1)+1e-8)
    nn_inputs(:,12)= (swl2(:)-swl2_clim_m(:)       - minmax(12,1))      &
                   / (minmax(12,2)                 - minmax(12,1)+1e-8)
    nn_inputs(:,13)= (swl3(:)                      - minmax(13,1))      &
                   / (minmax(13,2)                 - minmax(13,1)+1e-8)
    nn_inputs(:,14)= (swl3(:)-swl3_clim_m(:)       - minmax(14,1))      &
                   / (minmax(14,2)                 - minmax(14,1)+1e-8)
    nn_inputs(:,15)= (swl4(:)                      - minmax(15,1))      &
                   / (minmax(15,2)                 - minmax(15,1)+1e-8)
    nn_inputs(:,16)= (swl4(:)-swl4_clim_m(:)       - minmax(16,1))      &
                   / (minmax(16,2)                 - minmax(16,1)+1e-8) 
    nn_inputs(:,17)= (swl5(:)                      - minmax(17,1))      &
                   / (minmax(17,2)                 - minmax(17,1)+1e-8)
    nn_inputs(:,18)= (swl5(:)-swl5_clim_m(:)       - minmax(18,1))      &
                   / (minmax(18,2)                 - minmax(18,1)+1e-8) 

    nn_inputs(:,19)= (lai(:)                       - minmax(19,1))      &
                   / (minmax(19,2)                 - minmax(19,1)+1e-8)
    nn_inputs(:,20)= (lai(:)-lai_clim_m(:)         - minmax(20,1))      &
                   / (minmax(20,2)                 - minmax(20,1)+1e-8)
    nn_inputs(:,21)= (fuel(:)                      - minmax(21,1))      &
                   / (minmax(21,2)                 - minmax(21,1)+1e-8)
    nn_inputs(:,22)= (lightning_climatology(:)     - minmax(22,1))      &
                   / (minmax(22,2)                 - minmax(22,1)+1e-8)
    nn_inputs(:,23)= (LOG(gdp(:)+1._wp)            - minmax(23,1))      &
                   / (minmax(23,2)                 - minmax(23,1)+1e-8)
    nn_inputs(:,24)= (hdi(:)                       - minmax(24,1))      &
                   / (minmax(24,2)                 - minmax(24,1)+1e-8)
    nn_inputs(:,25)= (LOG(population(:)+1._wp)     - minmax(25,1))      &
                   / (minmax(25,2)                 - minmax(25,1)+1e-8)
    nn_inputs(:,26)= (elevation(:)                 - minmax(26,1))      &
                   / (minmax(26,2)                 - minmax(26,1)+1e-8) 

    nn_inputs(:,27:39)= pfts_fract(:,1:13)
     
    ! calculate FBA
    keep_prob = 1.
    nn_layer1_mask(:,:) = 1.
    nn_layer2_mask(:,:) = 1.
    nn_layer3_mask(:,:) = 1.
    nn_layer4_mask(:,:) = 1.
    nn_layer5_mask(:,:) = 1.
    nn_layer6_mask(:,:) = 1.

    ! RH constraint
    q_rel_air_constraint(:) = q_rel_air(:)
    CALL modified_sigmoid(q_rel_air_constraint, nc)

    CALL forward(                                    &
      nn_inputs,     nn_outputs,                     & 
      layer1_weight, layer1_bias, 'leaky_relu',      & 
      nn_layer1_out, nn_layer1_act, nn_layer1_mask,  &
      layer2_weight, layer2_bias, 'leaky_relu',      &
      nn_layer2_out, nn_layer2_act, nn_layer2_mask,  &
      layer3_weight, layer3_bias, 'leaky_relu',      &
      nn_layer3_out, nn_layer3_act, nn_layer3_mask,  &
      layer4_weight, layer4_bias, 'leaky_relu',      &
      nn_layer4_out, nn_layer4_act, nn_layer4_mask,  &
      layer5_weight, layer5_bias, 'leaky_relu',      &
      nn_layer5_out, nn_layer5_act, nn_layer5_mask,  &
      layer6_weight, layer6_bias, 'linear',          &
      nn_layer6_out, nn_layer6_act, nn_layer6_mask,  &
      nc, d_input, d_hidden, d_output, keep_prob,    &
      q_rel_air_constraint)

    burned_fract_infer = nn_outputs(:,1)/100.
    WHERE (burned_fract_infer < 0.) burned_fract_infer = 0.
    WHERE (burned_fract_infer > 1.) burned_fract_infer = 1.

    CALL forward(                                    &
      nn_inputs,     nn_outputs,                     & 
      layer1_weight_update, layer1_bias_update, 'leaky_relu',      & 
      nn_layer1_out, nn_layer1_act, nn_layer1_mask,  &
      layer2_weight_update, layer2_bias_update, 'leaky_relu',      &
      nn_layer2_out, nn_layer2_act, nn_layer2_mask,  &
      layer3_weight_update, layer3_bias_update, 'leaky_relu',      &
      nn_layer3_out, nn_layer3_act, nn_layer3_mask,  &
      layer4_weight_update, layer4_bias_update, 'leaky_relu',      &
      nn_layer4_out, nn_layer4_act, nn_layer4_mask,  &
      layer5_weight_update, layer5_bias_update, 'leaky_relu',      &
      nn_layer5_out, nn_layer5_act, nn_layer5_mask,  &
      layer6_weight_update, layer6_bias_update, 'linear',          &
      nn_layer6_out, nn_layer6_act, nn_layer6_mask,  &
      nc, d_input, d_hidden, d_output, keep_prob,    &
      q_rel_air_constraint)

    burned_fract_test = nn_outputs(:,1)/100.
    WHERE (burned_fract_test < 0.) burned_fract_test = 0.
    WHERE (burned_fract_test > 1.) burned_fract_test = 1.

    IF (is_infer .eq. 0) THEN
      burned_fract(:) = fba_smooth(:)
    ELSE
      burned_fract(:) = burned_fract_test(:)
    ENDIF

    ! saving in memory
    IF (is_saving .eq. 1) THEN
      V1_memory (:,memory_idx) = nn_inputs(:,1)
      V2_memory (:,memory_idx) = nn_inputs(:,2)
      V3_memory (:,memory_idx) = nn_inputs(:,3)
      V4_memory (:,memory_idx) = nn_inputs(:,4)
      V5_memory (:,memory_idx) = nn_inputs(:,5)
      V6_memory (:,memory_idx) = nn_inputs(:,6)
      V7_memory (:,memory_idx) = nn_inputs(:,7)
      V8_memory (:,memory_idx) = nn_inputs(:,8)
      V9_memory (:,memory_idx) = nn_inputs(:,9)
      V10_memory(:,memory_idx) = nn_inputs(:,10)
      V11_memory(:,memory_idx) = nn_inputs(:,11)
      V12_memory(:,memory_idx) = nn_inputs(:,12)
      V13_memory(:,memory_idx) = nn_inputs(:,13)
      V14_memory(:,memory_idx) = nn_inputs(:,14)
      V15_memory(:,memory_idx) = nn_inputs(:,15)
      V16_memory(:,memory_idx) = nn_inputs(:,16)
      V17_memory(:,memory_idx) = nn_inputs(:,17)
      V18_memory(:,memory_idx) = nn_inputs(:,18)
      V19_memory(:,memory_idx) = nn_inputs(:,19)
      V20_memory(:,memory_idx) = nn_inputs(:,20)
      V21_memory(:,memory_idx) = nn_inputs(:,21)
      V22_memory(:,memory_idx) = nn_inputs(:,22)
      V23_memory(:,memory_idx) = nn_inputs(:,23)
      V24_memory(:,memory_idx) = nn_inputs(:,24)
      V25_memory(:,memory_idx) = nn_inputs(:,25)
      V26_memory(:,memory_idx) = nn_inputs(:,26)
      V27_memory(:,memory_idx) = nn_inputs(:,27)
      V28_memory(:,memory_idx) = nn_inputs(:,28)
      V29_memory(:,memory_idx) = nn_inputs(:,29)
      V30_memory(:,memory_idx) = nn_inputs(:,30)
      V31_memory(:,memory_idx) = nn_inputs(:,31)
      V32_memory(:,memory_idx) = nn_inputs(:,32)
      V33_memory(:,memory_idx) = nn_inputs(:,33)
      V34_memory(:,memory_idx) = nn_inputs(:,34)
      V35_memory(:,memory_idx) = nn_inputs(:,35)
      V36_memory(:,memory_idx) = nn_inputs(:,36)
      V37_memory(:,memory_idx) = nn_inputs(:,37)
      V38_memory(:,memory_idx) = nn_inputs(:,38)
      V39_memory(:,memory_idx) = nn_inputs(:,39)
      RHC_memory(:,memory_idx) = q_rel_air_constraint(:)
      Y_memory  (:,memory_idx) = fba_smooth(:)
    ENDIF

    ! training
    IF (is_training .eq. 1) THEN

      ! random training batch from memory
      k=0
      DO i = 1,nc
      DO j = 1,INT(bs/nc)
        k = k+1

        !IF (k .eq. 1) THEN
        !  idx = memory_idx
        !ELSE
        call RANDOM_NUMBER(u)
        idx = 1 + INT(FLOOR((730)*u))
        !idx = 1 + INT(FLOOR((1095)*u))
        !ENDIF
        
        tr_inputs(k,1)  = V1_memory(i,idx)
        tr_inputs(k,2)  = V2_memory(i,idx)
        tr_inputs(k,3)  = V3_memory(i,idx)
        tr_inputs(k,4)  = V4_memory(i,idx)
        tr_inputs(k,5)  = V5_memory(i,idx)
        tr_inputs(k,6)  = V6_memory(i,idx)
        tr_inputs(k,7)  = V7_memory(i,idx)
        tr_inputs(k,8)  = V8_memory(i,idx)
        tr_inputs(k,9)  = V9_memory(i,idx)
        tr_inputs(k,10) = V10_memory(i,idx)
        tr_inputs(k,11) = V11_memory(i,idx)
        tr_inputs(k,12) = V12_memory(i,idx)
        tr_inputs(k,13) = V13_memory(i,idx)
        tr_inputs(k,14) = V14_memory(i,idx)
        tr_inputs(k,15) = V15_memory(i,idx)
        tr_inputs(k,16) = V16_memory(i,idx)
        tr_inputs(k,17) = V17_memory(i,idx)
        tr_inputs(k,18) = V18_memory(i,idx)
        tr_inputs(k,19) = V19_memory(i,idx)
        tr_inputs(k,20) = V20_memory(i,idx)
        tr_inputs(k,21) = V21_memory(i,idx)
        tr_inputs(k,22) = V22_memory(i,idx)
        tr_inputs(k,23) = V23_memory(i,idx)
        tr_inputs(k,24) = V24_memory(i,idx)
        tr_inputs(k,25) = V25_memory(i,idx)
        tr_inputs(k,26) = V26_memory(i,idx)
        tr_inputs(k,27) = V27_memory(i,idx)
        tr_inputs(k,28) = V28_memory(i,idx)
        tr_inputs(k,29) = V29_memory(i,idx)
        tr_inputs(k,30) = V30_memory(i,idx)
        tr_inputs(k,31) = V31_memory(i,idx)
        tr_inputs(k,32) = V32_memory(i,idx)
        tr_inputs(k,33) = V33_memory(i,idx)
        tr_inputs(k,34) = V34_memory(i,idx)
        tr_inputs(k,35) = V35_memory(i,idx)
        tr_inputs(k,36) = V36_memory(i,idx)
        tr_inputs(k,37) = V37_memory(i,idx)
        tr_inputs(k,38) = V38_memory(i,idx)
        tr_inputs(k,39) = V39_memory(i,idx)
        tr_q_rel_air_constraint(k) = RHC_memory(i,idx)
        tr_target(k) = Y_memory(i,idx)
      ENDDO
      ENDDO

      ! generate dropout mask
      keep_prob = 0.9
      CALL generate_dropout_mask(tr_layer1_mask, keep_prob, bs, d_hidden)
      CALL generate_dropout_mask(tr_layer2_mask, keep_prob, bs, d_hidden)
      CALL generate_dropout_mask(tr_layer3_mask, keep_prob, bs, d_hidden)
      CALL generate_dropout_mask(tr_layer4_mask, keep_prob, bs, d_hidden)
      CALL generate_dropout_mask(tr_layer5_mask, keep_prob, bs, d_hidden)
      tr_layer6_mask(:,:) = 1.


      ! Forward propagation
      CALL forward(                                    &
        tr_inputs,     tr_outputs,                     & 
        layer1_weight_update, layer1_bias_update, 'leaky_relu',      & 
        tr_layer1_out, tr_layer1_act, tr_layer1_mask,  &
        layer2_weight_update, layer2_bias_update, 'leaky_relu',      &
        tr_layer2_out, tr_layer2_act, tr_layer2_mask,  &
        layer3_weight_update, layer3_bias_update, 'leaky_relu',      &
        tr_layer3_out, tr_layer3_act, tr_layer3_mask,  &
        layer4_weight_update, layer4_bias_update, 'leaky_relu',      &
        tr_layer4_out, tr_layer4_act, tr_layer4_mask,  &
        layer5_weight_update, layer5_bias_update, 'leaky_relu',      &
        tr_layer5_out, tr_layer5_act, tr_layer5_mask,  &
        layer6_weight_update, layer6_bias_update, 'linear',          &
        tr_layer6_out, tr_layer6_act, tr_layer6_mask,  &
        bs, d_input, d_hidden, d_output, keep_prob,    &
        tr_q_rel_air_constraint)

      ! clear gradients
      dW1(:,:) = 0.;  db1(:) = 0.
      dW2(:,:) = 0.;  db2(:) = 0.
      dW3(:,:) = 0.;  db3(:) = 0.
      dW4(:,:) = 0.;  db4(:) = 0.
      dW5(:,:) = 0.;  db5(:) = 0.
      dW6(:,:) = 0.;  db6(:) = 0.

      ! Backward propagation
      CALL backward(                                                               &
        tr_inputs, tr_outputs, tr_target,                                          &
        layer1_weight_update, 'leaky_relu', tr_layer1_out, tr_layer1_act, tr_layer1_mask, &
        layer2_weight_update, 'leaky_relu', tr_layer2_out, tr_layer2_act, tr_layer2_mask, &
        layer3_weight_update, 'leaky_relu', tr_layer3_out, tr_layer3_act, tr_layer3_mask, &
        layer4_weight_update, 'leaky_relu', tr_layer4_out, tr_layer4_act, tr_layer4_mask, &
        layer5_weight_update, 'leaky_relu', tr_layer5_out, tr_layer5_act, tr_layer5_mask, &
        layer6_weight_update, 'linear',     tr_layer6_out, tr_layer6_act, tr_layer6_mask, &
        dW1, dW2, dW3, dW4, dW5, dW6,                                              &
        db1, db2, db3, db4, db5, db6,                                              &
        bs, d_input, d_hidden, d_output, keep_prob,                                &
        tr_q_rel_air_constraint)

      lr_rate = 1e-5
      CALL update(                                     &
        dW1, dW2, dW3, dW4, dW5, dW6,                  &
        db1, db2, db3, db4, db5, db6,                  &
        mW1, mW2, mW3, mW4, mW5, mW6,                  &
        mb1, mb2, mb3, mb4, mb5, mb6,                  &
        vW1, vW2, vW3, vW4, vW5, vW6,                  &
        vb1, vb2, vb3, vb4, vb5, vb6,                  &
        layer1_weight_update,  layer2_weight_update,  layer3_weight_update, &
        layer4_weight_update,  layer5_weight_update,  layer6_weight_update, &
        layer1_bias_update,    layer2_bias_update,    layer3_bias_update,   &
        layer4_bias_update,    layer5_bias_update,    layer6_bias_update,   &
        lr_rate, iteration, d_input, d_hidden, d_output)

      !layer1_weight_global(:,:) = layer1_weight_global(:,:)+1
      !layer1_weight_global(:,:) = layer1_weight_update(:,:)
      !layer2_weight_global(:,:) = layer2_weight_update(:,:)
      !layer3_weight_global(:,:) = layer3_weight_update(:,:)
      !layer4_weight_global(:,:) = layer4_weight_update(:,:)
      !layer5_weight_global(:,:) = layer5_weight_update(:,:)
      !layer6_weight_global(:,:) = layer6_weight_update(:,:)

      !layer1_bias_global(:) = layer1_bias_update(:)
      !layer2_bias_global(:) = layer2_bias_update(:)
      !layer3_bias_global(:) = layer3_bias_update(:)
      !layer4_bias_global(:) = layer4_bias_update(:)
      !layer5_bias_global(:) = layer5_bias_update(:)
      !layer6_bias_global(:) = layer6_bias_update(:)

      ! control tower
      !layer1_weight(10,20) = layer1_weight_global(10,20)
      !layer1_weight(10,20) = layer1_weight(10,20) + 1._wp
      !layer1_weight_global(10,20) = layer1_weight(10,20)
      !print*, layer1_weight_global(10,20)


    ENDIF

    
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

  ! NN-training
  SUBROUTINE generate_dropout_mask(x, keep_prob, d1, d2)
      IMPLICIT NONE
      REAL(wp), INTENT(IN)    :: keep_prob
      INTEGER,  INTENT(IN)    :: d1, d2
      REAL(wp), INTENT(INOUT) :: x(d1,d2)

      CALL RANDOM_NUMBER(x)
      WHERE (x < keep_prob)
        x = 1.
      ELSE WHERE
        x = 0.
      END WHERE
      END SUBROUTINE generate_dropout_mask


  SUBROUTINE ReLU(x, y, d1, d2)
    IMPLICIT NONE
    INTEGER,  INTENT(IN)    :: d1, d2
    REAL(wp), INTENT(IN)    :: x(d1,d2)
    REAL(wp), INTENT(inout) :: y(d1,d2)

    y(:,:) = 0.
    WHERE(x > 0.) y = x
    END SUBROUTINE ReLU


  SUBROUTINE leaky_ReLU(x, y, d1, d2)
    IMPLICIT NONE
    INTEGER,  INTENT(IN)    :: d1, d2
    REAL(wp), INTENT(IN)    :: x(d1,d2)
    REAL(wp), INTENT(INOUT) :: y(d1,d2)

    y(:,:) = x(:,:)
    WHERE(x < 0.) y = 0.01_wp*x
    END SUBROUTINE leaky_ReLU


  SUBROUTINE modified_sigmoid(x, d1)
      IMPLICIT NONE
      INTEGER,  INTENT(IN) :: d1
      REAL(wp), INTENT(INOUT) :: x(d1)

      x = x/100._wp
      WHERE(x < 0.) x = 0._wp
      WHERE(x > 1.) x = 1._wp
      x = 1._wp/(1._wp+EXP(-20._wp*(x-0.5_wp)))
      END SUBROUTINE modified_sigmoid


  SUBROUTINE single_forward_propagation(    &
    & input, weight, bias, activation_func, &
    & batch_size, d_in, d_out, keep_prob,   &
    & Z, output, mask)

    IMPLICIT NONE
    INTEGER,  INTENT(IN) :: batch_size, d_in, d_out
    REAL(wp), INTENT(IN) :: keep_prob
    REAL(wp), INTENT(IN) :: input(batch_size,d_in)
    REAL(wp), INTENT(IN) :: weight(d_in,d_out), bias(d_out)
    CHARACTER(len=*), INTENT(IN) :: activation_func
    REAL(wp), INTENT(IN)    :: mask  (batch_size,d_out)
    REAL(wp), INTENT(INOUT) :: Z     (batch_size,d_out)
    REAL(wp), INTENT(INOUT) :: output(batch_size,d_out)

    Z = matmul(input, weight) + spread(bias,1,batch_size)

    IF (activation_func .eq. 'relu') THEN
      CALL ReLU(Z, output, batch_size, d_out)
    ELSE if (activation_func .eq. 'leaky_relu') THEN
      CALL leaky_ReLU(Z, output, batch_size, d_out)
    ELSE
      output = Z
    ENDIF

    output = output*mask/keep_prob

    END SUBROUTINE single_forward_propagation


    SUBROUTINE forward(                                   &
      input, output,                                      &
      W1, b1, act1, layer1_out, layer1_act, layer1_mask,  &
      W2, b2, act2, layer2_out, layer2_act, layer2_mask,  &
      W3, b3, act3, layer3_out, layer3_act, layer3_mask,  &
      W4, b4, act4, layer4_out, layer4_act, layer4_mask,  &
      W5, b5, act5, layer5_out, layer5_act, layer5_mask,  &
      W6, b6, act6, layer6_out, layer6_act, layer6_mask,  &
      batch_size, d_in, d_hid, d_out, keep_prob,          &
      rh_constraint                                       )

    IMPLICIT NONE
    INTEGER,  INTENT(IN) :: batch_size, d_in, d_hid, d_out
    REAL(wp), INTENT(IN) :: keep_prob
    REAL(wp), INTENT(IN) :: input(batch_size,d_in)

    REAL(wp), INTENT(IN) :: W1(d_in, d_hid), b1(d_hid)
    REAL(wp), INTENT(IN) :: W2(d_hid,d_hid), b2(d_hid)
    REAL(wp), INTENT(IN) :: W3(d_hid,d_hid), b3(d_hid)
    REAL(wp), INTENT(IN) :: W4(d_hid,d_hid), b4(d_hid)
    REAL(wp), INTENT(IN) :: W5(d_hid,d_hid), b5(d_hid)
    REAL(wp), INTENT(IN) :: W6(d_hid,d_out), b6(d_out)

    CHARACTER(len=*),  INTENT(IN) :: act1, act2, act3
    CHARACTER(len=*),  INTENT(IN) :: act4, act5, act6

    REAL(wp), INTENT(IN) :: layer1_mask(batch_size, d_hid)
    REAL(wp), INTENT(IN) :: layer2_mask(batch_size, d_hid)
    REAL(wp), INTENT(IN) :: layer3_mask(batch_size, d_hid)
    REAL(wp), INTENT(IN) :: layer4_mask(batch_size, d_hid)
    REAL(wp), INTENT(IN) :: layer5_mask(batch_size, d_hid)
    REAL(wp), INTENT(IN) :: layer6_mask(batch_size, d_out)

    REAL(wp), INTENT(INOUT) :: layer1_out(batch_size, d_hid)
    REAL(wp), INTENT(INOUT) :: layer1_act(batch_size, d_hid)
    REAL(wp), INTENT(INOUT) :: layer2_out(batch_size, d_hid)
    REAL(wp), INTENT(INOUT) :: layer2_act(batch_size, d_hid)
    REAL(wp), INTENT(INOUT) :: layer3_out(batch_size, d_hid)
    REAL(wp), INTENT(INOUT) :: layer3_act(batch_size, d_hid)
    REAL(wp), INTENT(INOUT) :: layer4_out(batch_size, d_hid)
    REAL(wp), INTENT(INOUT) :: layer4_act(batch_size, d_hid)
    REAL(wp), INTENT(INOUT) :: layer5_out(batch_size, d_hid)
    REAL(wp), INTENT(INOUT) :: layer5_act(batch_size, d_hid)
    REAL(wp), INTENT(INOUT) :: layer6_out(batch_size, d_out)
    REAL(wp), INTENT(INOUT) :: layer6_act(batch_size, d_out)
    REAL(wp), INTENT(INOUT) :: output    (batch_size, d_out)
    REAL(wp), INTENT(IN)    :: rh_constraint(batch_size)


    CALL single_forward_propagation(input,      W1, b1, act1, &
      & batch_size, d_in,  d_hid,  keep_prob,                 &
      & layer1_out, layer1_act, layer1_mask                   )
    CALL single_forward_propagation(layer1_act, W2, b2, act2, &
      & batch_size, d_hid, d_hid,  keep_prob,                 &
      & layer2_out, layer2_act, layer2_mask                   )
    CALL single_forward_propagation(layer2_act, W3, b3, act3, &
      & batch_size, d_hid, d_hid,  keep_prob,                 &
      & layer3_out, layer3_act, layer3_mask                   )
    CALL single_forward_propagation(layer3_act, W4, b4, act4, &
      & batch_size, d_hid, d_hid,  keep_prob,                 &
      & layer4_out, layer4_act, layer4_mask                   )
    CALL single_forward_propagation(layer4_act, W5, b5, act5, &
      & batch_size, d_hid, d_hid,  keep_prob,                 &
      & layer5_out, layer5_act, layer5_mask                   )
    CALL single_forward_propagation(layer5_act, W6, b6, act6, &
      & batch_size, d_hid, d_out,  1._wp,                     &
      & layer6_out, layer6_act, layer6_mask                   )

    output(:,:) = layer6_act(:,:)

    ! constraints
    output(:,1) = output(:,1)*(1._wp-input(:,27)-input(:,39))
    output(:,1) = output(:,1)*input(:,21)
    output(:,1) = output(:,1)*(1._wp-rh_constraint)

    END SUBROUTINE forward


    SUBROUTINE dReLU(dy, z, dz, d1, d2)
      IMPLICIT NONE
      INTEGER,  INTENT(IN)    :: d1, d2
      REAL(wp), INTENT(IN)    :: dy(d1,d2)
      REAL(wp), INTENT(IN)    :: z (d1,d2)
      REAL(wp), INTENT(INOUT) :: dz(d1,d2)

      dz(:,:) = 0._wp
      WHERE(z > 0.) dz = dy
      END SUBROUTINE dReLU


    SUBROUTINE dleaky_ReLU(dy, z, dz, d1, d2)
      IMPLICIT NONE
      INTEGER,  INTENT(IN)    :: d1, d2
      REAL(wp), INTENT(IN)    :: dy(d1,d2)
      REAL(wp), INTENT(IN)    ::  z(d1,d2)
      REAL(wp), INTENT(INOUT) :: dz(d1,d2)

      dz(:,:) = dy(:,:)
      WHERE(z < 0.) dz = 0.01_wp*dy
      END SUBROUTINE dleaky_ReLU


    SUBROUTINE single_backward_propagation(        &
        & derivative_Y, mask, Z, X,                &
        & weight, activation_func,                 &
        & batch_size, d_in, d_out, keep_prob,      &
        & derivative_X, derivative_W, derivative_b )

      IMPLICIT NONE
      INTEGER,  INTENT(in) :: batch_size, d_in, d_out
      REAL(wp), INTENT(in) :: keep_prob
      REAL(wp), INTENT(in) :: derivative_Y(batch_size,d_out)
      REAL(wp), INTENT(in) :: mask        (batch_size,d_out)
      REAL(wp), INTENT(in) :: Z           (batch_size,d_out)
      REAL(wp), INTENT(in) :: X           (batch_size,d_in)
      REAL(wp), INTENT(in) :: weight      (d_in,d_out)
      CHARACTER(len=*),  INTENT(in) :: activation_func

      REAL(wp), INTENT(inout) :: derivative_X(batch_size,d_in)
      REAL(wp), INTENT(inout) :: derivative_W(d_in,d_out)
      REAL(wp), INTENT(inout) :: derivative_b(d_out)

      ! local
      REAL(wp) :: derivative_Y_(batch_size,d_out)
      REAL(wp) :: derivative_X_(batch_size,d_in)
      REAL(wp) :: derivative_W_(d_in,d_out)
      REAL(wp) :: derivative_b_(d_out)
      REAL(wp) :: derivative_Z_(batch_size,d_out)


      derivative_Y_ = derivative_Y*mask/keep_prob   ! dropout backprop

      IF (activation_func .eq. 'relu') THEN
        CALL dReLU(derivative_Y_, Z, derivative_Z_, batch_size, d_out)
      ELSE IF (activation_func .eq. 'leaky_relu') THEN
        CALL dleaky_ReLU(derivative_Y_, Z, derivative_Z_, batch_size, d_out)
      ELSE
        derivative_Z_ = derivative_Y_
      ENDIF  

      derivative_W_      = matmul(transpose(X), derivative_Z_)

      derivative_W_(:,:) = derivative_W_(:,:)/float(batch_size)
      derivative_W (:,:) = derivative_W(:,:)+derivative_W_(:,:)
      WHERE(derivative_W >  1e-1) derivative_W =  1e-1
      WHERE(derivative_W < -1e-1) derivative_W = -1e-1

      derivative_b_    = SUM(derivative_Z_, dim=1)
      derivative_b_(:) = derivative_b_(:)/float(batch_size)
      derivative_b (:) = derivative_b(:)+derivative_b_(:)
      WHERE(derivative_b >  1e-1) derivative_b =  1e-1
      WHERE(derivative_b < -1e-1) derivative_b = -1e-1

      derivative_X_ = matmul(derivative_Z_, transpose(weight))

      derivative_X(:,:) = derivative_X_(:,:)

      END SUBROUTINE single_backward_propagation


    SUBROUTINE backward(                                  &
          input, output, target_y,                        &
          W1, act1, layer1_out, layer1_act, layer1_mask,  &
          W2, act2, layer2_out, layer2_act, layer2_mask,  &
          W3, act3, layer3_out, layer3_act, layer3_mask,  &
          W4, act4, layer4_out, layer4_act, layer4_mask,  &
          W5, act5, layer5_out, layer5_act, layer5_mask,  &
          W6, act6, layer6_out, layer6_act, layer6_mask,  &
          dW1, dW2, dW3, dW4, dW5, dW6,                   &
          db1, db2, db3, db4, db5, db6,                   &
          batch_size, d_in, d_hid, d_out, keep_prob,      &
          rh_constraint                                   )

      IMPLICIT NONE
      INTEGER,  INTENT(IN) :: batch_size, d_in, d_hid, d_out
      REAL(wp), INTENT(IN) :: keep_prob
      REAL(wp), INTENT(IN) :: input   (batch_size,d_in)
      REAL(wp), INTENT(IN) :: output  (batch_size,d_out)
      REAL(wp), INTENT(IN) :: target_y(batch_size)

      REAL(wp), INTENT(IN) :: W1(d_in, d_hid)
      REAL(wp), INTENT(IN) :: W2(d_hid,d_hid)
      REAL(wp), INTENT(IN) :: W3(d_hid,d_hid)
      REAL(wp), INTENT(IN) :: W4(d_hid,d_hid)
      REAL(wp), INTENT(IN) :: W5(d_hid,d_hid)
      REAL(wp), INTENT(IN) :: W6(d_hid,d_out)
      character(len=*),  INTENT(IN) :: act1, act2, act3
      character(len=*),  INTENT(IN) :: act4, act5, act6

      REAL(wp), INTENT(IN) :: layer1_out(batch_size, d_hid)
      REAL(wp), INTENT(IN) :: layer1_act(batch_size, d_hid)
      REAL(wp), INTENT(IN) :: layer2_out(batch_size, d_hid)
      REAL(wp), INTENT(IN) :: layer2_act(batch_size, d_hid)
      REAL(wp), INTENT(IN) :: layer3_out(batch_size, d_hid)
      REAL(wp), INTENT(IN) :: layer3_act(batch_size, d_hid)
      REAL(wp), INTENT(IN) :: layer4_out(batch_size, d_hid)
      REAL(wp), INTENT(IN) :: layer4_act(batch_size, d_hid)
      REAL(wp), INTENT(IN) :: layer5_out(batch_size, d_hid)
      REAL(wp), INTENT(IN) :: layer5_act(batch_size, d_hid)
      REAL(wp), INTENT(IN) :: layer6_out(batch_size, d_out)
      REAL(wp), INTENT(IN) :: layer6_act(batch_size, d_out)

      REAL(wp), INTENT(IN) :: layer1_mask(batch_size, d_hid)
      REAL(wp), INTENT(IN) :: layer2_mask(batch_size, d_hid)
      REAL(wp), INTENT(IN) :: layer3_mask(batch_size, d_hid)
      REAL(wp), INTENT(IN) :: layer4_mask(batch_size, d_hid)
      REAL(wp), INTENT(IN) :: layer5_mask(batch_size, d_hid)
      REAL(wp), INTENT(IN) :: layer6_mask(batch_size, d_out)

      REAL(wp), INTENT(INOUT) :: dW1(d_in, d_hid), db1(d_hid)
      REAL(wp), INTENT(INOUT) :: dW2(d_hid,d_hid), db2(d_hid)
      REAL(wp), INTENT(INOUT) :: dW3(d_hid,d_hid), db3(d_hid)
      REAL(wp), INTENT(INOUT) :: dW4(d_hid,d_hid), db4(d_hid)
      REAL(wp), INTENT(INOUT) :: dW5(d_hid,d_hid), db5(d_hid)
      REAL(wp), INTENT(INOUT) :: dW6(d_hid,d_out), db6(d_out)
      REAL(wp), INTENT(IN)    :: rh_constraint(batch_size)

      ! local
      REAL(wp) :: dlayer1_act_(batch_size, d_hid)
      REAL(wp) :: dlayer2_act_(batch_size, d_hid)
      REAL(wp) :: dlayer3_act_(batch_size, d_hid)
      REAL(wp) :: dlayer4_act_(batch_size, d_hid)
      REAL(wp) :: dlayer5_act_(batch_size, d_hid)
      REAL(wp) :: dlayer6_act_(batch_size, d_out)
      REAL(wp) :: dinput_     (batch_size, d_in )

      dlayer6_act_(:,1) = 2.*(output(:,1)-target_y(:)*100)    ! MSE

      ! constraints
      dlayer6_act_(:,1) = dlayer6_act_(:,1)*(1._wp-input(:,27)-input(:,39))
      dlayer6_act_(:,1) = dlayer6_act_(:,1)*input(:,21)
      dlayer6_act_(:,1) = dlayer6_act_(:,1)*(1._wp-rh_constraint(:))

      CALL single_backward_propagation(                      &
        & dlayer6_act_, layer6_mask, layer6_out, layer5_act, &
        & W6, act6, batch_size, d_hid, d_out, 1._wp,         &
        & dlayer5_act_, dW6, db6                             )
      CALL single_backward_propagation(                      &
        & dlayer5_act_, layer5_mask, layer5_out, layer4_act, &
        & W5, act5, batch_size, d_hid, d_hid, keep_prob,     &
        & dlayer4_act_, dW5, db5                             )
      CALL single_backward_propagation(                      &
        & dlayer4_act_, layer4_mask, layer4_out, layer3_act, &
        & W4, act4, batch_size, d_hid, d_hid, keep_prob,     &
        & dlayer3_act_, dW4, db4                             )
      CALL single_backward_propagation(                      &
        & dlayer3_act_, layer3_mask, layer3_out, layer2_act, &
        & W3, act3, batch_size, d_hid, d_hid, keep_prob,     &
        & dlayer2_act_, dW3, db3                             )
      CALL single_backward_propagation(                      &
        & dlayer2_act_, layer2_mask, layer2_out, layer1_act, &
        & W2, act2, batch_size, d_hid, d_hid, keep_prob,     &
        & dlayer1_act_, dW2, db2                             )
      CALL single_backward_propagation(                      &
        & dlayer1_act_, layer1_mask, layer1_out, input,      &
        & W1, act1, batch_size, d_in,  d_hid, keep_prob,     &
        & dinput_     , dW1, db1                             )

      END SUBROUTINE backward


    SUBROUTINE Adam_optimizer_1d(                    &
          weights, gradients, momentum1, momentum2,  &
          learning_rate, n_iter, d1                  )

      IMPLICIT NONE
      INTEGER,  INTENT(IN)    :: n_iter, d1
      REAL(wp), INTENT(IN)    :: learning_rate
      REAL(wp), INTENT(IN)    :: gradients(d1)
      REAL(wp), INTENT(INOUT) :: momentum1(d1)
      REAL(wp), INTENT(INOUT) :: momentum2(d1)
      REAL(wp), INTENT(INOUT) :: weights  (d1)

      ! hyper-parameters
      REAL(wp) :: beta1, beta2
      beta1 = 0.9_wp     !0.9
      beta2 = 0.999_wp   !0.999

      momentum1(:)=beta1*momentum1(:)+(1._wp-beta1)*gradients(:)
      momentum2(:)=beta2*momentum2(:)+(1._wp-beta2)*(gradients(:)**2)

      !weights(:) = weights(:)-learning_rate                     &
      !           * (momentum1(:)  /(1-beta1**n_iter))           &
      !           / (((momentum2(:)/(1-beta2**n_iter))**0.5)+1e-7)
      weights(:) = weights(:)-learning_rate                                   &
                 * (momentum1(:)  /((1._wp-beta1**n_iter)+1e-8))              &
                 / (((momentum2(:)/((1._wp-beta2**n_iter)+1e-8))**0.5)+1e-8)


      END SUBROUTINE Adam_optimizer_1d


    SUBROUTINE Adam_optimizer_2d(                    &
          weights, gradients, momentum1, momentum2,  &
          learning_rate, n_iter, d1, d2              )

      IMPLICIT NONE
      INTEGER,  INTENT(IN)    :: n_iter, d1, d2
      REAL(wp), INTENT(IN)    :: learning_rate
      REAL(wp), INTENT(IN)    :: gradients(d1,d2)
      REAL(wp), INTENT(INOUT) :: momentum1(d1,d2)
      REAL(wp), INTENT(INOUT) :: momentum2(d1,d2)
      REAL(wp), INTENT(INOUT) :: weights  (d1,d2)

      ! hyper-parameters
      REAL(wp) :: beta1, beta2
      beta1 = 0.9_wp     !0.9
      beta2 = 0.999_wp   !0.999

      momentum1(:,:)=beta1*momentum1(:,:)+(1._wp-beta1)*gradients(:,:)
      momentum2(:,:)=beta2*momentum2(:,:)+(1._wp-beta2)*(gradients(:,:)**2)

      weights(:,:) = weights(:,:)-learning_rate                              &
                   * (momentum1(:,:)  /((1._wp-beta1**n_iter)+1e-8))         &
                   / (((momentum2(:,:)/((1._wp-beta2**n_iter)+1e-8))**0.5    &
                      )+1e-8)

      END SUBROUTINE Adam_optimizer_2d


    SUBROUTINE update(                       &
          dW1, dW2, dW3, dW4, dW5, dW6,      &
          db1, db2, db3, db4, db5, db6,      &
          mW1, mW2, mW3, mW4, mW5, mW6,      &
          mb1, mb2, mb3, mb4, mb5, mb6,      &
          vW1, vW2, vW3, vW4, vW5, vW6,      &
          vb1, vb2, vb3, vb4, vb5, vb6,      &
          W1,  W2,  W3,  W4,  W5,  W6,       &
          b1,  b2,  b3,  b4,  b5,  b6,       &
          learning_rate, iteration, d_in, d_hid, d_out)

      REAL(wp), INTENT(IN) :: learning_rate
      INTEGER,  INTENT(IN) :: iteration, d_in, d_hid, d_out

      REAL(wp), INTENT(IN) :: dW1(d_in, d_hid), db1(d_hid)
      REAL(wp), INTENT(IN) :: dW2(d_hid,d_hid), db2(d_hid)
      REAL(wp), INTENT(IN) :: dW3(d_hid,d_hid), db3(d_hid)
      REAL(wp), INTENT(IN) :: dW4(d_hid,d_hid), db4(d_hid)
      REAL(wp), INTENT(IN) :: dW5(d_hid,d_hid), db5(d_hid)
      REAL(wp), INTENT(IN) :: dW6(d_hid,d_out), db6(d_out)

      REAL(wp), INTENT(INOUT) :: mW1(d_in, d_hid), mb1(d_hid)
      REAL(wp), INTENT(INOUT) :: mW2(d_hid,d_hid), mb2(d_hid)
      REAL(wp), INTENT(INOUT) :: mW3(d_hid,d_hid), mb3(d_hid)
      REAL(wp), INTENT(INOUT) :: mW4(d_hid,d_hid), mb4(d_hid)
      REAL(wp), INTENT(INOUT) :: mW5(d_hid,d_hid), mb5(d_hid)
      REAL(wp), INTENT(INOUT) :: mW6(d_hid,d_out), mb6(d_out)

      REAL(wp), INTENT(INOUT) :: vW1(d_in, d_hid), vb1(d_hid)
      REAL(wp), INTENT(INOUT) :: vW2(d_hid,d_hid), vb2(d_hid)
      REAL(wp), INTENT(INOUT) :: vW3(d_hid,d_hid), vb3(d_hid)
      REAL(wp), INTENT(INOUT) :: vW4(d_hid,d_hid), vb4(d_hid)
      REAL(wp), INTENT(INOUT) :: vW5(d_hid,d_hid), vb5(d_hid)
      REAL(wp), INTENT(INOUT) :: vW6(d_hid,d_out), vb6(d_out)

      REAL(wp), INTENT(INOUT) :: W1(d_in, d_hid), b1(d_hid)
      REAL(wp), INTENT(INOUT) :: W2(d_hid,d_hid), b2(d_hid)
      REAL(wp), INTENT(INOUT) :: W3(d_hid,d_hid), b3(d_hid)
      REAL(wp), INTENT(INOUT) :: W4(d_hid,d_hid), b4(d_hid)
      REAL(wp), INTENT(INOUT) :: W5(d_hid,d_hid), b5(d_hid)
      REAL(wp), INTENT(INOUT) :: W6(d_hid,d_out), b6(d_out)

      CALL Adam_optimizer_2d(W1, dW1, mW1, vW1,  &
        & learning_rate, iteration, d_in,  d_hid )
      CALL Adam_optimizer_2d(W2, dW2, mW2, vW2,  &
        & learning_rate, iteration, d_hid, d_hid )
      CALL Adam_optimizer_2d(W3, dW3, mW3, vW3,  &
        & learning_rate, iteration, d_hid, d_hid )
      CALL Adam_optimizer_2d(W4, dW4, mW4, vW4,  &
        & learning_rate, iteration, d_hid, d_hid )
      CALL Adam_optimizer_2d(W5, dW5, mW5, vW5,  &
        & learning_rate, iteration, d_hid, d_hid )
      CALL Adam_optimizer_2d(W6, dW6, mW6, vW6,  &
        & learning_rate, iteration, d_hid, d_out )

      CALL Adam_optimizer_1d(b1, db1, mb1, vb1,  &
        & learning_rate, iteration, d_hid        )
      CALL Adam_optimizer_1d(b2, db2, mb2, vb2,  &
        & learning_rate, iteration, d_hid        )
      CALL Adam_optimizer_1d(b3, db3, mb3, vb3,  &
        & learning_rate, iteration, d_hid        )
      CALL Adam_optimizer_1d(b4, db4, mb4, vb4,  &
        & learning_rate, iteration, d_hid        )
      CALL Adam_optimizer_1d(b5, db5, mb5, vb5,  &
        & learning_rate, iteration, d_hid        )
      CALL Adam_optimizer_1d(b6, db6, mb6, vb6,  &
        & learning_rate, iteration, d_out        )

      END SUBROUTINE update

#endif
END MODULE mo_disturb_process
