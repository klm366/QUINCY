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
  PUBLIC :: burned_fract_dl, get_relative_humidity_air_v2, Koppen_Geiger_classification_lv2    ! <Rackhun Son>

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
    & land_use_c3ann,                       & ! in
    & land_use_c3nfx,                       & ! in
    & land_use_c3per,                       & ! in
    & land_use_c4ann,                       & ! in
    & land_use_c4per,                       & ! in
    & land_use_pastr,                       & ! in
    & land_use_primf,                       & ! in
    & land_use_primn,                       & ! in
    & land_use_range,                       & ! in
    & land_use_secdf,                       & ! in
    & land_use_secdn,                       & ! in
    & land_use_secma,                       & ! in
    & land_use_secmb,                       & ! in
    & land_use_urban,                       & ! in

    !& t_air_annual,                         & ! in
    !& rain_annual,                          & ! in
    & koppen,                                & ! inout

    & fba_smooth,                           & ! in

    & nc,                                   & ! in
    & iblk,                                 & ! in
    & tilename,                             & ! in
    & max_pft,                              & ! in

    & d_input,     d_hidden,                & ! in
    & d_output,    d_region,                & ! in
    & d_class,                              & ! in

    & layer1_weight, layer1_bias,            & ! in
    & layer2_weight, layer2_bias,            & ! in
    & layer3_weight, layer3_bias,            & ! in
    & layer4_weight, layer4_bias,            & ! in
    !& layer5_weight, layer5_bias,            & ! in
    !& layer6_weight, layer6_bias,            & ! in

    & layer1_weight_update, layer1_bias_update,            & ! in
    & layer2_weight_update, layer2_bias_update,            & ! in
    & layer3_weight_update, layer3_bias_update,            & ! in
    & layer4_weight_update, layer4_bias_update,            & ! in
    !& layer5_weight_update, layer5_bias_update,            & ! in
    !& layer6_weight_update, layer6_bias_update,            & ! in

    & mW1,           mb1,                    & ! in
    & mW2,           mb2,                    & ! in
    & mW3,           mb3,                    & ! in
    & mW4,           mb4,                    & ! in
    !& mW5,           mb5,                    & ! in
    !& mW6,           mb6,                    & ! in

    & vW1,           vb1,                    & ! in
    & vW2,           vb2,                    & ! in
    & vW3,           vb3,                    & ! in
    & vW4,           vb4,                    & ! in
    !& vW5,           vb5,                    & ! in
    !& vW6,           vb6,                    & ! in

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
    & V40_memory, V41_memory, V42_memory,    & ! inout
    & V43_memory, V44_memory, V45_memory,    & ! inout
    & V46_memory, V47_memory, V48_memory,    & ! inout
    & V49_memory, V50_memory, V51_memory,    & ! inout
    & V52_memory, V53_memory,                & ! inout
    & RHC_memory, Y_memory,   R_memory,      & ! inout
    & Koppen_memory,                         & ! inout

    & minmax,                                & ! in 
    & gfed_region,                           & ! in

    !& koppen,                                & ! inout
    & burned_fract,                          & ! inout
    & burned_fract_infer,                    & ! inout
    & burned_fract_test,                     & ! inout
    & weight)

    USE mo_global_variable

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
      & land_use_c3ann(:),                          &
      & land_use_c3nfx(:),                          &
      & land_use_c3per(:),                          &
      & land_use_c4ann(:),                          &
      & land_use_c4per(:),                          &
      & land_use_pastr(:),                          &
      & land_use_primf(:),                          &
      & land_use_primn(:),                          &
      & land_use_range(:),                          &
      & land_use_secdf(:),                          &
      & land_use_secdn(:),                          &
      & land_use_secma(:),                          &
      & land_use_secmb(:),                          &
      & land_use_urban(:),                          &

      !& t_air_annual(:,:),     rain_annual(:,:),    &
      & koppen(:),                                  &

      & layer1_weight(:,:),    layer1_bias(:),      &
      & layer2_weight(:,:),    layer2_bias(:),      &
      & layer3_weight(:,:),    layer3_bias(:),      &
      & layer4_weight(:,:),    layer4_bias(:),      &

      & fba_smooth(:),                              &
      & minmax(:,:),                                &
      & gfed_region(:)
      

    INTEGER, INTENT(in) ::                          &
      & nc,                  iblk,                  & ! in
      & d_input,             d_hidden,              & ! in
      & d_output,            d_region,              & ! in
      & d_class,                                    & ! in
      & iteration,           is_training,           & ! in
      & is_saving,           is_infer,              & ! in
      & memory_idx                                    ! in

    CHARACTER(len=*), INTENT(in)::                  &
      & tilename,            max_pft                  ! in

    ! Output Arguments
    REAL(wp), INTENT(inout) ::                      &
      & layer1_weight_update(:,:),    layer1_bias_update(:),      &  ! inout (training)
      & layer2_weight_update(:,:),    layer2_bias_update(:),      &  ! inout (training)
      & layer3_weight_update(:,:),    layer3_bias_update(:),      &  ! inout (training)
      & layer4_weight_update(:,:),    layer4_bias_update(:),      &  ! inout (training)

      & mW1(:,:),              mb1(:),              &  ! inout (training)
      & mW2(:,:),              mb2(:),              &  ! inout (training)
      & mW3(:,:),              mb3(:),              &  ! inout (training)
      & mW4(:,:),              mb4(:),              &  ! inout (training)

      & vW1(:,:),              vb1(:),              &  ! inout (training)
      & vW2(:,:),              vb2(:),              &  ! inout (training)
      & vW3(:,:),              vb3(:),              &  ! inout (training)
      & vW4(:,:),              vb4(:),              &  ! inout (training)

      & V1_memory (:,:),  V2_memory(:,:),  V3_memory(:,:),   & ! inout (training)
      & V4_memory (:,:),  V5_memory(:,:),  V6_memory(:,:),   & ! inout (training)
      & V7_memory (:,:),  V8_memory(:,:),  V9_memory(:,:),   & ! inout (training)
      & V10_memory(:,:), V11_memory(:,:), V12_memory(:,:),   & ! inout (training)
      & V13_memory(:,:), V14_memory(:,:), V15_memory(:,:),   & ! inout (training)
      & V16_memory(:,:), V17_memory(:,:), V18_memory(:,:),   & ! inout (training)
      & V19_memory(:,:), V20_memory(:,:), V21_memory(:,:),   & ! inout (training)
      & V22_memory(:,:), V23_memory(:,:), V24_memory(:,:),   & ! inout (training)
      & V25_memory(:,:), V26_memory(:,:), V27_memory(:,:),   & ! inout (training)
      & V28_memory(:,:), V29_memory(:,:), V30_memory(:,:),   & ! inout (training)
      & V31_memory(:,:), V32_memory(:,:), V33_memory(:,:),   & ! inout (training)
      & V34_memory(:,:), V35_memory(:,:), V36_memory(:,:),   & ! inout (training)
      & V37_memory(:,:), V38_memory(:,:), V39_memory(:,:),   & ! inout (training)
      & V40_memory(:,:), V41_memory(:,:), V42_memory(:,:),   & ! inout (training)
      & V43_memory(:,:), V44_memory(:,:), V45_memory(:,:),   & ! inout (training)
      & V46_memory(:,:), V47_memory(:,:), V48_memory(:,:),   & ! inout (training)
      & V49_memory(:,:), V50_memory(:,:), V51_memory(:,:),   & ! inout (training)
      & V52_memory(:,:), V53_memory(:,:),                    & ! inout (training)
      & RHC_memory(:,:), Y_memory  (:,:),   R_memory(:,:),   & ! inout (training)
      & Koppen_memory(:,:),                                  & ! inout (training)

      !& koppen(:),                                           & ! inout
      & burned_fract(:),                                     & ! inout
      & burned_fract_infer(:),                               & ! inout
      & burned_fract_test(:),                                & ! inout
      & weight

    ! Locals
    !INTEGER, PARAMETER :: bs = 1024
    INTEGER, PARAMETER :: bs = 1280
    REAL(wp) ::                                                   &
      & nn_inputs    (nc, d_input),                               &
      & nn_layer1_out(nc, d_hidden*d_class), nn_layer1_act(nc, d_hidden*d_class), &
      & nn_layer2_out(nc, d_hidden*d_class), nn_layer2_act(nc, d_hidden*d_class), &
      & nn_layer3_out(nc, d_hidden*d_class), nn_layer3_act(nc, d_hidden*d_class), &
      & nn_layer4_out(nc, d_output*d_class), nn_layer4_act(nc, d_output*d_class), &
      & nn_outputs   (nc, d_output),                              &

      & nn_layer1_mask(nc, d_hidden*d_class),                     &  ! dropout
      & nn_layer2_mask(nc, d_hidden*d_class),                     &
      & nn_layer3_mask(nc, d_hidden*d_class),                     &
      & nn_layer4_mask(nc, d_output*d_class),                     &

      !& tr_inputs    (bs, d_input),                               &  ! <training>
      !& tr_layer1_out(bs, d_hidden*d_class), tr_layer1_act(bs, d_hidden*d_class), &
      !& tr_layer2_out(bs, d_hidden*d_class), tr_layer2_act(bs, d_hidden*d_class), &
      !& tr_layer3_out(bs, d_hidden*d_class), tr_layer3_act(bs, d_hidden*d_class), &
      !& tr_layer4_out(bs, d_hidden*d_class), tr_layer4_act(bs, d_hidden*d_class), &
      !& tr_layer5_out(bs, d_hidden*d_class), tr_layer5_act(bs, d_hidden*d_class), &
      !& tr_layer6_out(bs, d_output*d_class), tr_layer6_act(bs, d_output*d_class), &
      !& tr_outputs   (bs, d_output),                              &

      !& tr_layer1_mask(bs, d_hidden*d_class),                     &  ! <training> dropout
      !& tr_layer2_mask(bs, d_hidden*d_class),                     &
      !& tr_layer3_mask(bs, d_hidden*d_class),                     &
      !& tr_layer4_mask(bs, d_hidden*d_class),                     &
      !& tr_layer5_mask(bs, d_hidden*d_class),                     &
      !& tr_layer6_mask(bs, d_output*d_class),                     &

      & tr_inputs    (bs, d_input),                               &  ! <training>
      & tr_layer1_out(1, d_hidden), tr_layer1_act(1, d_hidden),   &
      & tr_layer2_out(1, d_hidden), tr_layer2_act(1, d_hidden),   &
      & tr_layer3_out(1, d_hidden), tr_layer3_act(1, d_hidden),   &
      & tr_layer4_out(1, d_output), tr_layer4_act(1, d_output),   &
      & tr_outputs   (1, d_output),                               &

      & tr_layer1_mask(1, d_hidden),                              &  ! <training> dropout
      & tr_layer2_mask(1, d_hidden),                              &
      & tr_layer3_mask(1, d_hidden),                              &
      & tr_layer4_mask(1, d_output),                              &

      & dW1(d_input,  d_hidden*d_class),     db1(d_hidden*d_class),               &  ! derivatives
      & dW2(d_hidden, d_hidden*d_class),     db2(d_hidden*d_class),               &
      & dW3(d_hidden, d_hidden*d_class),     db3(d_hidden*d_class),               &
      & dW4(d_hidden, d_output*d_class),     db4(d_output*d_class)
                 
    INTEGER  :: i, j, k, idx, pid, getpid, c, min_idx, max_idx
    real(wp) :: lr_rate, keep_prob, loss, u, r
    real(wp) :: q_rel_air_constraint(nc)
    real(wp) :: tr_q_rel_air_constraint(bs), tr_target(bs), tr_weights(bs), tr_koppen(bs)
    real(wp) :: c_num(d_class)

    ! initialize
    DO i =1,nc
      c = INT(koppen(i))
      min_idx = 2*(c-1)+1
      max_idx = 2*(c-1)+2
      nn_inputs(i,1) = (t_air(i)                     - minmax(1,min_idx))       &
                     / (minmax(1,max_idx)            - minmax(1,min_idx)+1e-8)
      nn_inputs(i,2) = (t_air(i)-t_clim_m(i)         - minmax(2,min_idx))       &
                     / (minmax(2,max_idx)            - minmax(2,min_idx)+1e-8)
      nn_inputs(i,3) = (q_rel_air(i)                 - minmax(3,min_idx))       &
                     / (minmax(3,max_idx)            - minmax(3,min_idx)+1e-8)
      nn_inputs(i,4) = (q_rel_air(i)-q_rel_clim_m(i) - minmax(4,min_idx))       &
                     / (minmax(4,max_idx)            - minmax(4,min_idx)+1e-8)
      nn_inputs(i,5) = (q_air(i)                     - minmax(5,min_idx))       &
                     / (minmax(5,max_idx)            - minmax(5,min_idx)+1e-8)
      nn_inputs(i,6) = (q_air(i)-q_clim_m(i)         - minmax(6,min_idx))       &
                     / (minmax(6,max_idx)            - minmax(6,min_idx)+1e-8)
      nn_inputs(i,7) = (wind_10m(i)                  - minmax(7,min_idx))       &
                     / (minmax(7,max_idx)            - minmax(7,min_idx)+1e-8)
      nn_inputs(i,8) = (rain(i)                      - minmax(8,min_idx))       &
                     / (minmax(8,max_idx)            - minmax(8,min_idx)+1e-8)
      nn_inputs(i,9) = (swl1(i)                      - minmax(9,min_idx))       &
                     / (minmax(9,max_idx)            - minmax(9,min_idx)+1e-8)
      nn_inputs(i,10)= (swl1(i)-swl1_clim_m(i)       - minmax(10,min_idx))      &
                     / (minmax(10,max_idx)           - minmax(10,min_idx)+1e-8)
      nn_inputs(i,11)= (swl2(i)                      - minmax(11,min_idx))      &
                     / (minmax(11,max_idx)           - minmax(11,min_idx)+1e-8)
      nn_inputs(i,12)= (swl2(i)-swl2_clim_m(i)       - minmax(12,min_idx))      &
                     / (minmax(12,max_idx)           - minmax(12,min_idx)+1e-8)
      nn_inputs(i,13)= (swl3(i)                      - minmax(13,min_idx))      &
                     / (minmax(13,max_idx)           - minmax(13,min_idx)+1e-8)
      nn_inputs(i,14)= (swl3(i)-swl3_clim_m(i)       - minmax(14,min_idx))      &
                     / (minmax(14,max_idx)           - minmax(14,min_idx)+1e-8)
      nn_inputs(i,15)= (swl4(i)                      - minmax(15,min_idx))      &
                     / (minmax(15,max_idx)           - minmax(15,min_idx)+1e-8)
      nn_inputs(i,16)= (swl4(i)-swl4_clim_m(i)       - minmax(16,min_idx))      &
                     / (minmax(16,max_idx)           - minmax(16,min_idx)+1e-8) 
      nn_inputs(i,17)= (swl5(i)                      - minmax(17,min_idx))      &
                     / (minmax(17,max_idx)           - minmax(17,min_idx)+1e-8)
      nn_inputs(i,18)= (swl5(i)-swl5_clim_m(i)       - minmax(18,min_idx))      &
                     / (minmax(18,max_idx)           - minmax(18,min_idx)+1e-8) 

      nn_inputs(i,19)= (lai(i)                       - minmax(19,min_idx))      &
                     / (minmax(19,max_idx)           - minmax(19,min_idx)+1e-8)
      nn_inputs(i,20)= (lai(i)-lai_clim_m(i)         - minmax(20,min_idx))      &
                     / (minmax(20,max_idx)           - minmax(20,min_idx)+1e-8)
      nn_inputs(i,21)= (fuel(i)                      - minmax(21,min_idx))      &
                     / (minmax(21,max_idx)           - minmax(21,min_idx)+1e-8)
      nn_inputs(i,22)= (lightning_climatology(i)     - minmax(22,min_idx))      &
                     / (minmax(22,max_idx)           - minmax(22,min_idx)+1e-8)
      nn_inputs(i,23)= (LOG(gdp(i)+1._wp)            - minmax(23,min_idx))      &
                     / (minmax(23,max_idx)           - minmax(23,min_idx)+1e-8)
      nn_inputs(i,24)= (hdi(i)                       - minmax(24,min_idx))      &
                     / (minmax(24,max_idx)           - minmax(24,min_idx)+1e-8)
      nn_inputs(i,25)= (LOG(population(i)+1._wp)     - minmax(25,min_idx))      &
                     / (minmax(25,max_idx)           - minmax(25,min_idx)+1e-8)
      nn_inputs(i,26)= (elevation(i)                 - minmax(26,min_idx))      &
                     / (minmax(26,max_idx)           - minmax(26,min_idx)+1e-8)
    ENDDO 

    nn_inputs(:,27:39)= pfts_fract(:,1:13)

    nn_inputs(:,40) = land_use_c3ann(:)
    nn_inputs(:,41) = land_use_c3nfx(:)
    nn_inputs(:,42) = land_use_c3per(:)
    nn_inputs(:,43) = land_use_c4ann(:)
    nn_inputs(:,44) = land_use_c4per(:)
    nn_inputs(:,45) = land_use_pastr(:)
    nn_inputs(:,46) = land_use_primf(:)
    nn_inputs(:,47) = land_use_primn(:)
    nn_inputs(:,48) = land_use_range(:)
    nn_inputs(:,49) = land_use_secdf(:)
    nn_inputs(:,50) = land_use_secdn(:)
    nn_inputs(:,51) = land_use_secma(:)
    nn_inputs(:,52) = land_use_secmb(:)
    nn_inputs(:,53) = land_use_urban(:)

    ! calculate FBA
    keep_prob = 1.
    nn_layer1_mask(:,:) = 1._wp
    nn_layer2_mask(:,:) = 1._wp
    nn_layer3_mask(:,:) = 1._wp
    nn_layer4_mask(:,:) = 1._wp

    ! RH constraint
    q_rel_air_constraint(:) = q_rel_air(:)
    CALL modified_sigmoid(q_rel_air_constraint, nc)
    WHERE (q_rel_air_constraint > 0.6) q_rel_air_constraint = 1._wp

    DO c = 1,d_class
      CALL forward(                                                   &
        nn_inputs,     nn_outputs(:,c:c),                             & 
        layer1_weight (:,d_hidden*(c-1)+1:d_hidden*c),                &
        layer1_bias     (d_hidden*(c-1)+1:d_hidden*c), 'leaky_relu',  & 
        nn_layer1_out (:,d_hidden*(c-1)+1:d_hidden*c),                &
        nn_layer1_act (:,d_hidden*(c-1)+1:d_hidden*c),                &
        nn_layer1_mask(:,d_hidden*(c-1)+1:d_hidden*c),                &

        layer2_weight (:,d_hidden*(c-1)+1:d_hidden*c),                &
        layer2_bias     (d_hidden*(c-1)+1:d_hidden*c), 'leaky_relu',  & 
        nn_layer2_out (:,d_hidden*(c-1)+1:d_hidden*c),                &
        nn_layer2_act (:,d_hidden*(c-1)+1:d_hidden*c),                &
        nn_layer2_mask(:,d_hidden*(c-1)+1:d_hidden*c),                &

        layer3_weight (:,d_hidden*(c-1)+1:d_hidden*c),                &
        layer3_bias     (d_hidden*(c-1)+1:d_hidden*c), 'leaky_relu',  & 
        nn_layer3_out (:,d_hidden*(c-1)+1:d_hidden*c),                &
        nn_layer3_act (:,d_hidden*(c-1)+1:d_hidden*c),                &
        nn_layer3_mask(:,d_hidden*(c-1)+1:d_hidden*c),                &

        layer4_weight (:,d_output*(c-1)+1:d_output*c),                &
        layer4_bias     (d_output*(c-1)+1:d_output*c), 'linear',      & 
        nn_layer4_out (:,d_output*(c-1)+1:d_output*c),                &
        nn_layer4_act (:,d_output*(c-1)+1:d_output*c),                &
        nn_layer4_mask(:,d_output*(c-1)+1:d_output*c),                &

        nc, d_input, d_hidden, d_output, keep_prob,                   &
        q_rel_air_constraint)
    ENDDO

    ! Koppen-Geiger classification
    !CALL Koppen_Geiger_classification(nc, t_air_annual,  rain_annual, koppen)

    DO i = 1, nc
      c = INT(koppen(i))
      burned_fract_infer(i) = nn_outputs(i,c)/100.
      if (c .eq. 0) Then
        print* ,'wrong'
      ENDIF
    ENDDO
    WHERE (burned_fract_infer < 0.) burned_fract_infer = 0.
    WHERE (burned_fract_infer > 1.) burned_fract_infer = 1.


    DO c = 1,d_class
      CALL forward(                                                        &
        nn_inputs,            nn_outputs(:,c:c),                           & 
        layer1_weight_update(:,d_hidden*(c-1)+1:d_hidden*c),               &
        layer1_bias_update    (d_hidden*(c-1)+1:d_hidden*c), 'leaky_relu', & 
        nn_layer1_out (:,d_hidden*(c-1)+1:d_hidden*c),                     &
        nn_layer1_act (:,d_hidden*(c-1)+1:d_hidden*c),                     &
        nn_layer1_mask(:,d_hidden*(c-1)+1:d_hidden*c),                     &

        layer2_weight_update(:,d_hidden*(c-1)+1:d_hidden*c),               &
        layer2_bias_update    (d_hidden*(c-1)+1:d_hidden*c), 'leaky_relu', &
        nn_layer2_out (:,d_hidden*(c-1)+1:d_hidden*c),                     &
        nn_layer2_act (:,d_hidden*(c-1)+1:d_hidden*c),                     &
        nn_layer2_mask(:,d_hidden*(c-1)+1:d_hidden*c),                     &

        layer3_weight_update(:,d_hidden*(c-1)+1:d_hidden*c),               &
        layer3_bias_update    (d_hidden*(c-1)+1:d_hidden*c), 'leaky_relu', &
        nn_layer3_out (:,d_hidden*(c-1)+1:d_hidden*c),                     &
        nn_layer3_act (:,d_hidden*(c-1)+1:d_hidden*c),                     &
        nn_layer3_mask(:,d_hidden*(c-1)+1:d_hidden*c),                     &

        layer4_weight_update(:,d_output*(c-1)+1:d_output*c),               &
        layer4_bias_update    (d_output*(c-1)+1:d_output*c), 'linear',     &
        nn_layer4_out (:,d_output*(c-1)+1:d_output*c),                     &
        nn_layer4_act (:,d_output*(c-1)+1:d_output*c),                     &
        nn_layer4_mask(:,d_output*(c-1)+1:d_output*c),                     &

        nc, d_input, d_hidden, d_output, keep_prob,                        &
        q_rel_air_constraint)
    ENDDO

    DO i = 1, nc
      c = INT(koppen(i))
      burned_fract_test(i) = nn_outputs(i,c)/100.
    ENDDO
    WHERE (burned_fract_test < 0.) burned_fract_test = 0.
    WHERE (burned_fract_test > 1.) burned_fract_test = 1.

    !IF (is_infer .eq. 0) THEN
    !  burned_fract(:) = fba_smooth(:)
    !ELSE
    !  burned_fract(:) = burned_fract_test(:)
    !ENDIF
    burned_fract(:) = burned_fract_test(:)


    ! saving in memory
    IF ((is_saving .eq. 1) .and. (max_pft .eq. tilename)) THEN
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
      V40_memory(:,memory_idx) = nn_inputs(:,40)
      V41_memory(:,memory_idx) = nn_inputs(:,41)
      V42_memory(:,memory_idx) = nn_inputs(:,42)
      V43_memory(:,memory_idx) = nn_inputs(:,43)
      V44_memory(:,memory_idx) = nn_inputs(:,44)
      V45_memory(:,memory_idx) = nn_inputs(:,45)
      V46_memory(:,memory_idx) = nn_inputs(:,46)
      V47_memory(:,memory_idx) = nn_inputs(:,47)
      V48_memory(:,memory_idx) = nn_inputs(:,48)
      V49_memory(:,memory_idx) = nn_inputs(:,49)
      V50_memory(:,memory_idx) = nn_inputs(:,50)
      V51_memory(:,memory_idx) = nn_inputs(:,51)
      V52_memory(:,memory_idx) = nn_inputs(:,52)
      V53_memory(:,memory_idx) = nn_inputs(:,53)
      RHC_memory(:,memory_idx) = q_rel_air_constraint(:)
      Y_memory  (:,memory_idx) = fba_smooth(:)
      R_memory  (:,memory_idx) = gfed_region(:)
      Koppen_memory(:,memory_idx) = koppen(:)

    ENDIF

    pid = getpid()
    pid = mod(pid, 8)

    IF ((is_training .eq. 1) .and. (pid .eq. 0) .and. (iblk .eq. first_iblk) .and. (max_pft .eq. tilename)) THEN
      print*, 'train'
      print*, layer1_weight_update(10,128*2+20), iblk, first_iblk
      CALL load_global_batch()
      tr_inputs(:,:) = batch_global(:,1:53)
      tr_q_rel_air_constraint(:) = batch_global(:,54)
      tr_target(:) = batch_global(:,55)
      tr_koppen(:) = batch_global(:,57)

      DO i=1,bs
      !  IF (((batch_global(i,56) .eq. 8) .or. (batch_global(i,56) .eq. 9) .or. (batch_global(i,56) .eq. 14)) &
      !    .and. (tr_target(i) .gt. 0) &
      !  IF ((batch_global(i,56) .eq. 8) .and. (tr_target(i) .gt. 0)) THEN
      !  IF ((batch_global(i,56) .eq. 8) .or. (batch_global(i,56) .eq. 9) .or. &
      !  IF ((batch_global(i,42) .eq. 8) .or. (batch_global(i,42) .eq. 9)) THEN!  .or. &
      !      (batch_global(i,42) .eq. 7) .or. (batch_global(i,42) .eq. 11) .or. &
      !    ) THEN
      !    tr_weights(i) = 4._wp

        IF ((batch_global(i,56) .eq. 8) .and. (batch_global(i,57) .eq. 4) .and. (batch_global(i,55) .gt. 0)) THEN
          tr_weights(i) = 20._wp
        ELSEIF ((batch_global(i,56) .eq. 9) .and. (batch_global(i,57) .eq. 3) .and. (batch_global(i,55) .gt. 0)) THEN
          tr_weights(i) = 20._wp
        ELSEIF ((batch_global(i,56) .eq. 5) .and. (batch_global(i,57) .eq. 3) .and. (batch_global(i,55) .eq. 0)) THEN
          tr_weights(i) = 10._wp
        ELSEIF ((batch_global(i,57) .eq. 5) .and. (batch_global(i,55) .eq. 0)) THEN
          tr_weights(i) = 10._wp
        ELSEIF ((batch_global(i,57) .eq. 6) .and. (batch_global(i,55) .eq. 0)) THEN
          tr_weights(i) = 10._wp
        ELSE
          tr_weights(i) = 1._wp
        ENDIF
      ENDDO
      !tr_weights(:) = 1._wp

      ! clear gradients
      dW1(:,:) = 0._wp;  db1(:) = 0._wp
      dW2(:,:) = 0._wp;  db2(:) = 0._wp
      dW3(:,:) = 0._wp;  db3(:) = 0._wp
      dW4(:,:) = 0._wp;  db4(:) = 0._wp

      c_num(:) = 1e-8
      DO i = 1,bs

        c = INT(tr_koppen(i))
        r = batch_global(i,56)

        !IF (c .ne. 0) THEN
        IF ((c .ne. 0) .and. (r .ge. 1._wp)) THEN
        !IF ((c .ne. 0) .and. ((r .eq. 8) .or. (r .eq. 9) .or. (r .eq. 5))) THEN
        !IF ((r .ge. 1) .and. ((c .ge. 4) .and. (c .le. 6))) THEN
          c_num(c) = c_num(c) + 1._wp

          ! generate dropout mask
          keep_prob = 0.7
          CALL generate_dropout_mask(tr_layer1_mask, keep_prob, 1, d_hidden)
          CALL generate_dropout_mask(tr_layer2_mask, keep_prob, 1, d_hidden)
          CALL generate_dropout_mask(tr_layer3_mask, keep_prob, 1, d_hidden)
          tr_layer4_mask(:,:) = 1._wp

          ! Forward propagation
          CALL forward(                                                 &
            tr_inputs(i:i,:),     tr_outputs(1:1,:),                    & 
            layer1_weight_update(:,d_hidden*(c-1)+1:d_hidden*c),        &
            layer1_bias_update    (d_hidden*(c-1)+1:d_hidden*c),        &
            'leaky_relu', tr_layer1_out, tr_layer1_act, tr_layer1_mask, &

            layer2_weight_update(:,d_hidden*(c-1)+1:d_hidden*c),        &
            layer2_bias_update    (d_hidden*(c-1)+1:d_hidden*c),        &
            'leaky_relu', tr_layer2_out, tr_layer2_act, tr_layer2_mask, &

            layer3_weight_update(:,d_hidden*(c-1)+1:d_hidden*c),        &
            layer3_bias_update    (d_hidden*(c-1)+1:d_hidden*c),        &
            'leaky_relu', tr_layer3_out, tr_layer3_act, tr_layer3_mask, &

            layer4_weight_update(:,d_output*(c-1)+1:d_output*c),        &
            layer4_bias_update    (d_output*(c-1)+1:d_output*c),        &
            'linear',     tr_layer4_out, tr_layer4_act, tr_layer4_mask, &

            1, d_input, d_hidden, d_output, keep_prob,                  &
            tr_q_rel_air_constraint(i:i))

          ! Backward propagation
          CALL backward(                                                &
            tr_inputs(i:i,:), tr_outputs(1:1,:), tr_target(i:i),        &
            layer1_weight_update(:,d_hidden*(c-1)+1:d_hidden*c),        &
            'leaky_relu', tr_layer1_out, tr_layer1_act, tr_layer1_mask, &

            layer2_weight_update(:,d_hidden*(c-1)+1:d_hidden*c),        &
            'leaky_relu', tr_layer2_out, tr_layer2_act, tr_layer2_mask, &

            layer3_weight_update(:,d_hidden*(c-1)+1:d_hidden*c),        &
            'leaky_relu', tr_layer3_out, tr_layer3_act, tr_layer3_mask, &

            layer4_weight_update(:,d_output*(c-1)+1:d_output*c),        &
            'linear',     tr_layer4_out, tr_layer4_act, tr_layer4_mask, &

            dW1(:,d_hidden*(c-1)+1:d_hidden*c), dW2(:,d_hidden*(c-1)+1:d_hidden*c), &
            dW3(:,d_hidden*(c-1)+1:d_hidden*c), dW4(:,d_output*(c-1)+1:d_output*c), &

            db1(d_hidden*(c-1)+1:d_hidden*c),   db2(d_hidden*(c-1)+1:d_hidden*c), &
            db3(d_hidden*(c-1)+1:d_hidden*c),   db4(d_output*(c-1)+1:d_output*c), & 

            1, d_input, d_hidden, d_output, keep_prob,                            &
            tr_q_rel_air_constraint(i:i), tr_weights(i:i))

        ENDIF
      ENDDO


      lr_rate = 1e-4
      DO c = 1, d_class
        dW1(:,d_hidden*(c-1)+1:d_hidden*c) = dW1(:,d_hidden*(c-1)+1:d_hidden*c) / c_num(c)
        dW2(:,d_hidden*(c-1)+1:d_hidden*c) = dW2(:,d_hidden*(c-1)+1:d_hidden*c) / c_num(c)
        dW3(:,d_hidden*(c-1)+1:d_hidden*c) = dW3(:,d_hidden*(c-1)+1:d_hidden*c) / c_num(c)
        dW4(:,d_output*(c-1)+1:d_output*c) = dW4(:,d_output*(c-1)+1:d_output*c) / c_num(c)

        db1(d_hidden*(c-1)+1:d_hidden*c) = db1(d_hidden*(c-1)+1:d_hidden*c) / c_num(c)
        db2(d_hidden*(c-1)+1:d_hidden*c) = db2(d_hidden*(c-1)+1:d_hidden*c) / c_num(c)
        db3(d_hidden*(c-1)+1:d_hidden*c) = db3(d_hidden*(c-1)+1:d_hidden*c) / c_num(c)
        db4(d_output*(c-1)+1:d_output*c) = db4(d_output*(c-1)+1:d_output*c) / c_num(c)

        CALL update(                                                              &
          dW1(:,d_hidden*(c-1)+1:d_hidden*c), dW2(:,d_hidden*(c-1)+1:d_hidden*c), &
          dW3(:,d_hidden*(c-1)+1:d_hidden*c), dW4(:,d_output*(c-1)+1:d_output*c), &
          db1  (d_hidden*(c-1)+1:d_hidden*c), db2  (d_hidden*(c-1)+1:d_hidden*c), &
          db3  (d_hidden*(c-1)+1:d_hidden*c), db4  (d_output*(c-1)+1:d_output*c), &

          mW1(:,d_hidden*(c-1)+1:d_hidden*c), mW2(:,d_hidden*(c-1)+1:d_hidden*c), &
          mW3(:,d_hidden*(c-1)+1:d_hidden*c), mW4(:,d_output*(c-1)+1:d_output*c), &
          mb1  (d_hidden*(c-1)+1:d_hidden*c), mb2  (d_hidden*(c-1)+1:d_hidden*c), &
          mb3  (d_hidden*(c-1)+1:d_hidden*c), mb4  (d_output*(c-1)+1:d_output*c), &

          vW1(:,d_hidden*(c-1)+1:d_hidden*c), vW2(:,d_hidden*(c-1)+1:d_hidden*c), &
          vW3(:,d_hidden*(c-1)+1:d_hidden*c), vW4(:,d_output*(c-1)+1:d_output*c), &
          vb1  (d_hidden*(c-1)+1:d_hidden*c), vb2  (d_hidden*(c-1)+1:d_hidden*c), &
          vb3  (d_hidden*(c-1)+1:d_hidden*c), vb4  (d_output*(c-1)+1:d_output*c), &

          layer1_weight_update(:,d_hidden*(c-1)+1:d_hidden*c),                    &
          layer2_weight_update(:,d_hidden*(c-1)+1:d_hidden*c),                    &
          layer3_weight_update(:,d_hidden*(c-1)+1:d_hidden*c),                    &
          layer4_weight_update(:,d_output*(c-1)+1:d_output*c),                    &

          layer1_bias_update    (d_hidden*(c-1)+1:d_hidden*c),                    &
          layer2_bias_update    (d_hidden*(c-1)+1:d_hidden*c),                    &
          layer3_bias_update    (d_hidden*(c-1)+1:d_hidden*c),                    &
          layer4_bias_update    (d_output*(c-1)+1:d_output*c),                    &

          lr_rate, iteration, d_input, d_hidden, d_output)

      ENDDO
      
      CALL update_trainable_parameters(               &
        & layer1_weight_update,  layer1_bias_update,  &
        & layer2_weight_update,  layer2_bias_update,  &
        & layer3_weight_update,  layer3_bias_update,  &
        & layer4_weight_update,  layer4_bias_update   )


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
  SUBROUTINE Koppen_Geiger_classification_lv2(nc, t_air_annual, rain_annual, lat, koppen)
    IMPLICIT NONE
    INTEGER,  INTENT(IN)    :: nc
    REAL(wp), INTENT(IN)    :: t_air_annual(nc, 366), rain_annual(nc, 366), lat(nc)
    REAL(wp), INTENT(INOUT) :: koppen(nc)

    ! local
    REAL(wp) :: t_monthly(12), p_monthly(12)
    REAL(wp) :: t_cold, t_hot, t_ann, p_ann, p_summer, p_winter, p_threshold
    REAL(wp) :: is_boreal, t_ondjfm, t_amjjas, p_threshold1, p_threshold2, p_threshold3
    REAL(wp) :: p_dry, p_sdry, p_swet, p_wdry, p_wwet
    REAL(wp) :: Am_or_w, Bw_or_s, Cw_or_f
    INTEGER  :: i, Aw_or_s

    ! Koppen Geiger classification
    ! Tropical
    !   1: Rainforest (Af)
    !   2: Monsoon (Am)
    !   3: Savannah (Aw)

    ! 2: Arid (B)
    ! 3: Temperate (C)
    ! 4: Cold (D)
    ! 5: Polar (E)

    DO i = 1, nc
      t_monthly(1)  = SUM(t_air_annual(i, 92:121))/30._wp - 273.15   ! Apr
      t_monthly(2)  = SUM(t_air_annual(i,122:152))/31._wp - 273.15   ! May
      t_monthly(3)  = SUM(t_air_annual(i,153:182))/30._wp - 273.15   ! Jun
      t_monthly(4)  = SUM(t_air_annual(i,183:213))/31._wp - 273.15   ! Jul
      t_monthly(5)  = SUM(t_air_annual(i,214:244))/31._wp - 273.15   ! Aug
      t_monthly(6)  = SUM(t_air_annual(i,245:274))/30._wp - 273.15   ! Sep
      t_monthly(7)  = SUM(t_air_annual(i,275:305))/31._wp - 273.15   ! Oct
      t_monthly(8)  = SUM(t_air_annual(i,306:335))/30._wp - 273.15   ! Nov
      t_monthly(9)  = SUM(t_air_annual(i,336:366))/31._wp - 273.15   ! Dec
      t_monthly(10) = SUM(t_air_annual(i,  1: 31))/31._wp - 273.15   ! Jan
      t_monthly(11) = SUM(t_air_annual(i, 32: 60))/29._wp - 273.15   ! Feb
      t_monthly(12) = SUM(t_air_annual(i, 61: 91))/31._wp - 273.15   ! Mar

      p_monthly(1)  = SUM(rain_annual(i, 92:121)) * 86400            ! Apr
      p_monthly(2)  = SUM(rain_annual(i,122:152)) * 86400            ! May
      p_monthly(3)  = SUM(rain_annual(i,153:182)) * 86400            ! Jun
      p_monthly(4)  = SUM(rain_annual(i,183:213)) * 86400            ! Jul
      p_monthly(5)  = SUM(rain_annual(i,214:244)) * 86400            ! Aug
      p_monthly(6)  = SUM(rain_annual(i,245:274)) * 86400            ! Sep
      p_monthly(7)  = SUM(rain_annual(i,275:305)) * 86400            ! Oct
      p_monthly(8)  = SUM(rain_annual(i,306:335)) * 86400            ! Nov
      p_monthly(9)  = SUM(rain_annual(i,336:366)) * 86400            ! Dec
      p_monthly(10) = SUM(rain_annual(i,  1: 31)) * 86400            ! Jan
      p_monthly(11) = SUM(rain_annual(i, 32: 60)) * 86400            ! Feb
      p_monthly(12) = SUM(rain_annual(i, 61: 91)) * 86400            ! Mar

      t_cold   = MINVAL(t_monthly)
      t_hot    = MAXVAL(t_monthly)
      t_ann    = SUM(t_monthly)/12._wp
      p_ann    = SUM(p_monthly)

      p_dry    = MINVAL(p_monthly)

      !IF (lat(i) .ge. 0) THEN 
      !  p_summer = SUM(p_monthly(1:6))
      !  p_winter = SUM(p_monthly(7:12))
      !  p_sdry   = MINVAL(p_monthly(1:6))
      !  p_swet   = MAXVAL(p_monthly(1:6))
      !  p_wdry   = MINVAL(p_monthly(7:12))
      !  p_wwet   = MAXVAL(p_monthly(7:12))
      !ELSE
      !  p_summer = SUM(p_monthly(7:12))
      !  p_winter = SUM(p_monthly(1:6))
      !  p_sdry   = MINVAL(p_monthly(7:12))
      !  p_swet   = MAXVAL(p_monthly(7:12))
      !  p_wdry   = MINVAL(p_monthly(1:6))
      !  p_wwet   = MAXVAL(p_monthly(1:6))
      !ENDIF

      !IF (p_summer .gt. SUM(p_monthly)*0.7) THEN
      !  p_threshold = 2._wp*t_ann + 28._wp  
      !ELSEIF (p_winter .gt. SUM(p_monthly)*0.7) THEN
      !  p_threshold = 2._wp*t_ann
      !ELSE
      !  p_threshold = 2._wp*t_ann + 14._wp
      !ENDIF

      t_ondjfm = SUM(t_monthly(7:12))
      t_amjjas = SUM(t_monthly(1:6))
      IF (t_amjjas .gt. t_ondjfm) THEN
        is_boreal = 1
      ELSE
        is_boreal = 0
      ENDIF

      IF (is_boreal .eq. 1) THEN
        p_summer = SUM(p_monthly(1:6))
        p_winter = SUM(p_monthly(7:12))
        p_sdry   = MINVAL(p_monthly(1:6))
        p_swet   = MAXVAL(p_monthly(1:6))
        p_wdry   = MINVAL(p_monthly(7:12))
        p_wwet   = MAXVAL(p_monthly(7:12))
      ELSE
        p_summer = SUM(p_monthly(7:12))
        p_winter = SUM(p_monthly(1:6))
        p_sdry   = MINVAL(p_monthly(7:12))
        p_swet   = MAXVAL(p_monthly(7:12))
        p_wdry   = MINVAL(p_monthly(1:6))
        p_wwet   = MAXVAL(p_monthly(1:6))
      ENDIF

      IF (p_winter .gt. SUM(p_monthly)*0.7) THEN
        p_threshold1 = 2._wp*t_ann
      ELSE
        p_threshold1 = 0.
      ENDIF

      IF (p_summer .gt. SUM(p_monthly)*0.7) THEN
        p_threshold2 = 2._wp*t_ann + 28._wp
      ELSE
        p_threshold2 = 0.
      ENDIF

      IF ((p_winter .le. SUM(p_monthly)*0.7) .and. (p_summer .le. SUM(p_monthly)*0.7)) THEN
        p_threshold3 = 2._wp*t_ann + 14._wp
      ELSE
        p_threshold3 = 0.
      ENDIF

      p_threshold = p_threshold1+p_threshold2+p_threshold3


      !!! classification
      ! Tropical
      IF (t_cold .ge. 18) THEN
        Am_or_w = 100 - p_ann/25
        IF(p_dry .ge. 60) THEN
          koppen(i) = 1._wp
        ELSEIF (p_dry .ge. Am_or_w) THEN
          koppen(i) = 2._wp
        ELSE
          Aw_or_s = MINLOC(p_monthly, dim=1)
          IF (Aw_or_s .le. 6) THEN
            koppen(i) = 3._wp
          ELSE
            koppen(i) = 4._wp
          ENDIF
        ENDIF

      ! Temperate
      ELSEIF ((t_hot .ge. 10) .and. (t_cold .gt. 0) .and. (t_cold .lt. 18)) THEN
        Cw_or_f = p_swet/10
        IF ((p_sdry .lt. 40) .and. (p_sdry .lt. p_wwet/3)) THEN
          koppen(i) = 7._wp
        ELSEIF (p_wdry .lt. Cw_or_f) THEN
          koppen(i) = 8._wp
        ELSE
          koppen(i) = 9._wp
        ENDIF

      ! Cold
      ELSEIF ((t_hot .ge. 10) .and. (t_cold .le. 0)) THEN
        IF ((p_sdry .lt. 40) .and. (p_sdry .lt. p_wwet/3)) THEN
          koppen(i) = 10._wp
        ELSEIF (p_wdry .lt. p_swet/10) THEN
          koppen(i) = 11._wp
        ELSE
          koppen(i) = 12._wp
        ENDIF

      ! Polar 
      ELSEIF (t_hot .lt. 10) THEN
        !IF (t_hot .gt. 0) THEN
        koppen(i) = 13._wp
        !ELSE
        !  koppen(i) = 13._wp
        !ENDIF
      ENDIF

      ! Arid
      IF (p_ann .lt. 10*p_threshold) THEN
        Bw_or_s = 5*p_threshold
        IF (p_ann .lt. Bw_or_s) THEN
          Koppen(i) = 5._wp
        ELSE
          Koppen(i) = 6._wp
        ENDIF
      ENDIF

    ENDDO
        
  END SUBROUTINE Koppen_Geiger_classification_lv2


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
    batch_size, d_in, d_hid, d_out, keep_prob,          &
    rh_constraint                                       )

    IMPLICIT NONE
    INTEGER,  INTENT(IN) :: batch_size, d_in, d_hid, d_out
    REAL(wp), INTENT(IN) :: keep_prob
    REAL(wp), INTENT(IN) :: input(batch_size,d_in)

    REAL(wp), INTENT(IN) :: W1(d_in, d_hid), b1(d_hid)
    REAL(wp), INTENT(IN) :: W2(d_hid,d_hid), b2(d_hid)
    REAL(wp), INTENT(IN) :: W3(d_hid,d_hid), b3(d_hid)
    REAL(wp), INTENT(IN) :: W4(d_hid,d_out), b4(d_out)

    CHARACTER(len=*),  INTENT(IN) :: act1, act2, act3, act4

    REAL(wp), INTENT(IN) :: layer1_mask(batch_size, d_hid)
    REAL(wp), INTENT(IN) :: layer2_mask(batch_size, d_hid)
    REAL(wp), INTENT(IN) :: layer3_mask(batch_size, d_hid)
    REAL(wp), INTENT(IN) :: layer4_mask(batch_size, d_out)

    REAL(wp), INTENT(INOUT) :: layer1_out(batch_size, d_hid)
    REAL(wp), INTENT(INOUT) :: layer1_act(batch_size, d_hid)
    REAL(wp), INTENT(INOUT) :: layer2_out(batch_size, d_hid)
    REAL(wp), INTENT(INOUT) :: layer2_act(batch_size, d_hid)
    REAL(wp), INTENT(INOUT) :: layer3_out(batch_size, d_hid)
    REAL(wp), INTENT(INOUT) :: layer3_act(batch_size, d_hid)
    REAL(wp), INTENT(INOUT) :: layer4_out(batch_size, d_out)
    REAL(wp), INTENT(INOUT) :: layer4_act(batch_size, d_out)
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
      & batch_size, d_hid, d_out,  1._wp,                     &
      & layer4_out, layer4_act, layer4_mask                   )

    output(:,:) = layer4_act(:,:)

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
    !WHERE(derivative_W >  1e-1) derivative_W =  1e-1
    !WHERE(derivative_W < -1e-1) derivative_W = -1e-1

    derivative_b_    = SUM(derivative_Z_, dim=1)
    derivative_b_(:) = derivative_b_(:)/float(batch_size)
    derivative_b (:) = derivative_b(:)+derivative_b_(:)
    !WHERE(derivative_b >  1e-1) derivative_b =  1e-1
    !WHERE(derivative_b < -1e-1) derivative_b = -1e-1

    derivative_X_ = matmul(derivative_Z_, transpose(weight))
    derivative_X(:,:) = derivative_X_(:,:)

  END SUBROUTINE single_backward_propagation


  SUBROUTINE backward(                                  &
    input, output, target_y,                        &
    W1, act1, layer1_out, layer1_act, layer1_mask,  &
    W2, act2, layer2_out, layer2_act, layer2_mask,  &
    W3, act3, layer3_out, layer3_act, layer3_mask,  &
    W4, act4, layer4_out, layer4_act, layer4_mask,  &
    dW1, dW2, dW3, dW4, db1, db2, db3, db4,         &
    batch_size, d_in, d_hid, d_out, keep_prob,      &
    rh_constraint, weights                          )

    IMPLICIT NONE
    INTEGER,  INTENT(IN) :: batch_size, d_in, d_hid, d_out
    REAL(wp), INTENT(IN) :: keep_prob
    REAL(wp), INTENT(IN) :: input   (batch_size,d_in)
    REAL(wp), INTENT(IN) :: output  (batch_size,d_out)
    REAL(wp), INTENT(IN) :: target_y(batch_size)

    REAL(wp), INTENT(IN) :: W1(d_in, d_hid)
    REAL(wp), INTENT(IN) :: W2(d_hid,d_hid)
    REAL(wp), INTENT(IN) :: W3(d_hid,d_hid)
    REAL(wp), INTENT(IN) :: W4(d_hid,d_out)

    character(len=*),  INTENT(IN) :: act1, act2, act3, act4

    REAL(wp), INTENT(IN) :: layer1_out(batch_size, d_hid)
    REAL(wp), INTENT(IN) :: layer1_act(batch_size, d_hid)
    REAL(wp), INTENT(IN) :: layer2_out(batch_size, d_hid)
    REAL(wp), INTENT(IN) :: layer2_act(batch_size, d_hid)
    REAL(wp), INTENT(IN) :: layer3_out(batch_size, d_hid)
    REAL(wp), INTENT(IN) :: layer3_act(batch_size, d_hid)
    REAL(wp), INTENT(IN) :: layer4_out(batch_size, d_out)
    REAL(wp), INTENT(IN) :: layer4_act(batch_size, d_out)

    REAL(wp), INTENT(IN) :: layer1_mask(batch_size, d_hid)
    REAL(wp), INTENT(IN) :: layer2_mask(batch_size, d_hid)
    REAL(wp), INTENT(IN) :: layer3_mask(batch_size, d_hid)
    REAL(wp), INTENT(IN) :: layer4_mask(batch_size, d_out)

    REAL(wp), INTENT(INOUT) :: dW1(d_in, d_hid), db1(d_hid)
    REAL(wp), INTENT(INOUT) :: dW2(d_hid,d_hid), db2(d_hid)
    REAL(wp), INTENT(INOUT) :: dW3(d_hid,d_hid), db3(d_hid)
    REAL(wp), INTENT(INOUT) :: dW4(d_hid,d_out), db4(d_out)

    REAL(wp), INTENT(IN)    :: rh_constraint(batch_size)
    REAL(wp), INTENT(IN)    :: weights(batch_size)

    ! local
    REAL(wp) :: dlayer1_act_(batch_size, d_hid)
    REAL(wp) :: dlayer2_act_(batch_size, d_hid)
    REAL(wp) :: dlayer3_act_(batch_size, d_hid)
    REAL(wp) :: dlayer4_act_(batch_size, d_out)
    REAL(wp) :: dinput_     (batch_size, d_in )

    !dlayer6_act_(:,1) = 2.*(output(:,1)-target_y(:)*100)    ! MSE
    dlayer4_act_(:,1) = weights(:)*2.*(output(:,1)-target_y(:)*100)    ! MSE
    !dlayer4_act_(:,1) = weights(:)*2.*(output(:,1)-target_y(:)*100)/(target_y(:)*100)    ! MSE
    !dlayer6_act_(:,1) = weights(:)*100.*0.5/SQRT(ABS(output(:,1)-target_y(:)))    ! MSE

    ! constraints
    dlayer4_act_(:,1) = dlayer4_act_(:,1)*(1._wp-input(:,27)-input(:,39))
    dlayer4_act_(:,1) = dlayer4_act_(:,1)*input(:,21)
    dlayer4_act_(:,1) = dlayer4_act_(:,1)*(1._wp-rh_constraint(:))

    CALL single_backward_propagation(                      &
      & dlayer4_act_, layer4_mask, layer4_out, layer3_act, &
      & W4, act4, batch_size, d_hid, d_out, 1._wp,         &
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
    REAL(wp), INTENT(INOUT) :: gradients(d1)
    REAL(wp), INTENT(INOUT) :: momentum1(d1)
    REAL(wp), INTENT(INOUT) :: momentum2(d1)
    REAL(wp), INTENT(INOUT) :: weights  (d1)

    ! hyper-parameters
    REAL(wp) :: beta1, beta2
    beta1 = 0.9_wp     !0.9
    beta2 = 0.999_wp   !0.999

    WHERE(gradients >  1e-1) gradients =  1e-1
    WHERE(gradients < -1e-1) gradients = -1e-1

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
    REAL(wp), INTENT(INOUT) :: gradients(d1,d2)
    REAL(wp), INTENT(INOUT) :: momentum1(d1,d2)
    REAL(wp), INTENT(INOUT) :: momentum2(d1,d2)
    REAL(wp), INTENT(INOUT) :: weights  (d1,d2)

    ! hyper-parameters
    REAL(wp) :: beta1, beta2
    beta1 = 0.9_wp     !0.9
    beta2 = 0.999_wp   !0.999

    WHERE(gradients >  1e-1) gradients =  1e-1
    WHERE(gradients < -1e-1) gradients = -1e-1

    momentum1(:,:)=beta1*momentum1(:,:)+(1._wp-beta1)*gradients(:,:)
    momentum2(:,:)=beta2*momentum2(:,:)+(1._wp-beta2)*(gradients(:,:)**2)

    weights(:,:) = weights(:,:)-learning_rate                              &
                 * (momentum1(:,:)  /((1._wp-beta1**n_iter)+1e-8))         &
                 / (((momentum2(:,:)/((1._wp-beta2**n_iter)+1e-8))**0.5    &
                    )+1e-8)

  END SUBROUTINE Adam_optimizer_2d


  SUBROUTINE update(                        &
    dW1, dW2, dW3, dW4, db1, db2, db3, db4, &
    mW1, mW2, mW3, mW4, mb1, mb2, mb3, mb4, &
    vW1, vW2, vW3, vW4, vb1, vb2, vb3, vb4, &
    W1,  W2,  W3,  W4,  b1,  b2,  b3,  b4,  &
    learning_rate, iteration, d_in, d_hid, d_out)

    REAL(wp), INTENT(IN) :: learning_rate
    INTEGER,  INTENT(IN) :: iteration, d_in, d_hid, d_out

    REAL(wp), INTENT(INOUT) :: dW1(d_in, d_hid), db1(d_hid)
    REAL(wp), INTENT(INOUT) :: dW2(d_hid,d_hid), db2(d_hid)
    REAL(wp), INTENT(INOUT) :: dW3(d_hid,d_hid), db3(d_hid)
    REAL(wp), INTENT(INOUT) :: dW4(d_hid,d_out), db4(d_out)

    REAL(wp), INTENT(INOUT) :: mW1(d_in, d_hid), mb1(d_hid)
    REAL(wp), INTENT(INOUT) :: mW2(d_hid,d_hid), mb2(d_hid)
    REAL(wp), INTENT(INOUT) :: mW3(d_hid,d_hid), mb3(d_hid)
    REAL(wp), INTENT(INOUT) :: mW4(d_hid,d_out), mb4(d_out)

    REAL(wp), INTENT(INOUT) :: vW1(d_in, d_hid), vb1(d_hid)
    REAL(wp), INTENT(INOUT) :: vW2(d_hid,d_hid), vb2(d_hid)
    REAL(wp), INTENT(INOUT) :: vW3(d_hid,d_hid), vb3(d_hid)
    REAL(wp), INTENT(INOUT) :: vW4(d_hid,d_out), vb4(d_out)

    REAL(wp), INTENT(INOUT) :: W1(d_in, d_hid), b1(d_hid)
    REAL(wp), INTENT(INOUT) :: W2(d_hid,d_hid), b2(d_hid)
    REAL(wp), INTENT(INOUT) :: W3(d_hid,d_hid), b3(d_hid)
    REAL(wp), INTENT(INOUT) :: W4(d_hid,d_out), b4(d_out)

    CALL Adam_optimizer_2d(W1, dW1, mW1, vW1,  &
      & learning_rate, iteration, d_in,  d_hid )
    CALL Adam_optimizer_2d(W2, dW2, mW2, vW2,  &
      & learning_rate, iteration, d_hid, d_hid )
    CALL Adam_optimizer_2d(W3, dW3, mW3, vW3,  &
      & learning_rate, iteration, d_hid, d_hid )
    CALL Adam_optimizer_2d(W4, dW4, mW4, vW4,  &
      & learning_rate, iteration, d_hid, d_out )

    CALL Adam_optimizer_1d(b1, db1, mb1, vb1,  &
      & learning_rate, iteration, d_hid        )
    CALL Adam_optimizer_1d(b2, db2, mb2, vb2,  &
      & learning_rate, iteration, d_hid        )
    CALL Adam_optimizer_1d(b3, db3, mb3, vb3,  &
      & learning_rate, iteration, d_hid        )
    CALL Adam_optimizer_1d(b4, db4, mb4, vb4,  &
      & learning_rate, iteration, d_out        )

  END SUBROUTINE update

#endif
END MODULE mo_disturb_process
