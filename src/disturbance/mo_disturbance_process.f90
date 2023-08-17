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

    & pfts_fract,                           & ! in
    & fuel,                                 & ! in
    & lightning_climatology,                & ! in
    & elevation,                            & ! in

    & population,                           & ! in
    & gdp,                                  & ! in
    & hdi,                                  & ! in

    & fba_smooth,                           & ! in
    & nc,                                   & ! in
    & iblk,                                 & ! in
    & tilename,                             & ! in

    & d_input,     d_hidden,                & ! in
    & d_output,    d_region,                & ! in

    & layer1_weight, layer1_bias,            & ! in
    & layer2_weight, layer2_bias,            & ! in
    & layer3_weight, layer3_bias,            & ! in
    & layer4_weight, layer4_bias,            & ! in
    & layer5_weight, layer5_bias,            & ! in
    & layer6_weight, layer6_bias,            & ! in

    & minmax,                                & ! in 
    & gfed_region,                           & ! in

    & burned_fract                           ) ! inout

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

      & pfts_fract(:,:),                            &
      & fuel(:),                                    &
      & lightning_climatology(:),                   &
      & elevation(:),                               &

      & population(:),                              &
      & gdp(:),                                     &
      & hdi(:),                                     &

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
      & nc,                  iblk,                  & ! in
      & d_input,             d_hidden,              & ! in
      & d_output,            d_region

    CHARACTER(len=*), INTENT(in)::                  &
      & tilename

    ! Output Arguments
    REAL(wp), INTENT(inout) :: burned_fract(:)

    ! Locals
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
      & nn_layer6_mask(nc, d_output)

                 
    INTEGER  :: i, j, k, idx
    real(wp) :: lr_rate, keep_prob
    real(wp) :: q_rel_air_constraint(nc)

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

    burned_fract = nn_outputs(:,1)/100.
    WHERE (burned_fract < 0.) burned_fract = 0.
    WHERE (burned_fract > 1.) burned_fract = 1.
    !burned_fract(:) = fba_smooth(:)

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

#endif
END MODULE mo_disturb_process
