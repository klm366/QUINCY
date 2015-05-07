!>
!! @brief Contains the interfaces to the soil processes
!!
!! @author
!!  Reiner Schnur, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Reiner Schnur (2013-04-04)
!!
MODULE mo_soil_interface
#ifndef __NO_JSBACH__

  USE mo_jsb_model_class, ONLY: t_jsb_model
  USE mo_jsb_tile_class,  ONLY: t_jsb_tile
  USE mo_jsb_base,        ONLY: get_model
  USE mo_jsb_control,     ONLY: debug_on
  USE mo_kind,            ONLY: wp
  USE mo_exception,       ONLY: message, finish
  USE mo_util,            ONLY: one_of

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: update_soil_1, update_soil_2

  CHARACTER(len=*), PARAMETER :: modname = 'mo_soil_interface'

CONTAINS

  ! TBD: At the moment, the hydrology for soil and canopy is handled together in the soil component (as in echam6).
  ! This will be changed later.
  SUBROUTINE update_soil_1( tile, iblk, ics, ice, nc, dtime, steplen)

    USE mo_soil_process, ONLY: update_surf_down, update_soiltemp, update_surf_up

    TYPE(t_jsb_tile),  POINTER    :: tile
    INTEGER,           INTENT(in) :: iblk, ics, ice, nc
    REAL(wp),          INTENT(in) :: dtime, steplen

    ! Local variables
    !
    TYPE(t_jsb_model), POINTER    :: model

    ! Pointers to variables in memory
    REAL(wp), POINTER ::       &
      & fract          (:),    &
      & t_srf_unfilt   (:),    &
      & tte_corr       (:),    &
      & t_air          (:),    &
      & q_air          (:),    &
      & wind_10m       (:),    &
      & rain           (:),    &
      & snow           (:),    &
      & wsr_soil       (:),    &
      & wsr_srf        (:),    &
      & wfract_srf     (:),    &
      & wsn_soil       (:),    &
      & wsn_srf        (:),    &
      & sfract_srf     (:),    &
      & wsn_can        (:),    &
      & lai            (:),    &
      & fract_fpc_max  (:),    &
      & hcap_grnd      (:),    &
      & hcap_grnd_old  (:),    &
      & evapotrans     (:),    &
      & trans          (:),    &
      & evapopot       (:),    &
      & vhcap          (:),    &
      & hdiff          (:),    &
      & t_soil_acoef   (:,:),  &
      & t_soil_bcoef   (:,:),  &
      & t_soil         (:,:),  &
      & grnd_hflx      (:),    &
      & grnd_hflx_old  (:),    &
      & oro_stddev     (:),    &
      & runoff_acc     (:),    &
      & drainage_acc   (:),    &
      & ws             (:),    &
      & ws_l           (:,:),  &
      & ws_fc_l        (:,:),  &
      & ws_pwp_l       (:,:),  &
      & soil_depth     (:),    &
      & ws_root        (:),    &
      & ws_fc_root     (:),    &
      & root_depth     (:),    &
      & shc            (:),    &
      & mpot           (:),    &
      & bclapp         (:),    &
      & vpor           (:),    &
      & psdi           (:),    &
      & overflow       (:),    &
      & glacier        (:),    &
      & wsn_glacier    (:)


    ! Locally allocated vectors
    !
    REAL(wp) :: &
      wsr_max(nc), &
      evapotrans_no_snow_skin(nc), &
      glacier_precip_minus_evap(nc), &
      snow_acc(nc), &
      snow_melt_acc(nc), &
      glacier_runoff_acc(nc), &
      melt_water_excess(nc), &
      snow_melt(nc), &
      runoff_hd(nc), &
      drainage_hd(nc), &
      skin_res_evap(nc), &
      reduced_evap(nc)

    LOGICAL :: &
      is_present(nc), &
      is_glacier(nc)

    INTEGER :: jllog

    LOGICAL :: lstart
    LOGICAL :: l_tpe_open = .FALSE.
    LOGICAL :: l_tpe_closed = .FALSE.

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':update_soil_1'

    IF (one_of('soil', tile%components) < 1) RETURN

    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting update')

    model => get_model(tile%owner_model_id)

    IF (model%config%usecase == 'jsbach_lite_tpe_open')   l_tpe_open = .TRUE.
    IF (model%config%usecase == 'jsbach_lite_tpe_closed') l_tpe_closed = .TRUE.

    ! Get reference to variables for current block
    !
    fract                => tile%mem%fract                               (ics:ice,    iblk)
    IF (ALL(fract(:) < EPSILON(1._wp))) RETURN

    t_srf_unfilt   => tile%mem%srf  %energy_state   %t_srf_unfilt        (ics:ice,    iblk)
    tte_corr       => tile%mem%srf  %energy_state   %tte_corr            (ics:ice,    iblk)
    t_air          => tile%mem      %atm2land       %t_air               (ics:ice,    iblk)
    q_air          => tile%mem      %atm2land       %q_air               (ics:ice,    iblk)
    wind_10m       => tile%mem      %atm2land       %wind_10m            (ics:ice,    iblk)
    rain           => tile%mem      %atm2land       %rain                (ics:ice,    iblk)
    snow           => tile%mem      %atm2land       %snow                (ics:ice,    iblk)
    wsr_soil       => tile%mem%soil %water_state    %wsr_soil            (ics:ice,    iblk)
    wsr_srf        => tile%mem%srf  %water_state    %wsr_srf             (ics:ice,    iblk)
    wfract_srf     => tile%mem%srf  %water_state    %wfract_srf          (ics:ice,    iblk)
    wsn_soil       => tile%mem%soil %water_state    %wsn_soil            (ics:ice,    iblk)
    wsn_srf        => tile%mem%srf  %water_state    %wsn_srf             (ics:ice,    iblk)
    sfract_srf     => tile%mem%srf  %water_state    %sfract_srf          (ics:ice,    iblk)
    wsn_can        => tile%mem%veg  %water_state    %wsn_can             (ics:ice,    iblk)
    lai            => tile%mem%veg  %phys_state     %lai                 (ics:ice,    iblk)
    fract_fpc_max  => tile%mem%veg  %phys_state     %fract_fpc_max       (ics:ice,    iblk)
    hcap_grnd      => tile%mem%soil %energy_state   %hcap_grnd           (ics:ice,    iblk)
    hcap_grnd_old  => tile%mem%soil %energy_state   %hcap_grnd_old       (ics:ice,    iblk)
    evapotrans     => tile%mem%srf  %water_flux     %evapotrans          (ics:ice,    iblk)
    trans          => tile%mem%srf  %water_flux     %trans               (ics:ice,    iblk)
    evapopot       => tile%mem%srf  %water_flux     %evapopot            (ics:ice,    iblk)
    vhcap          => tile%mem%soil %phys_state     %vhcap               (ics:ice,    iblk)
    hdiff          => tile%mem%soil %phys_state     %hdiff               (ics:ice,    iblk)
    t_soil_acoef   => tile%mem%soil %energy_state   %t_soil_acoef        (ics:ice, :, iblk)
    t_soil_bcoef   => tile%mem%soil %energy_state   %t_soil_bcoef        (ics:ice, :, iblk)
    t_soil         => tile%mem%soil %energy_state   %t_soil              (ics:ice, :, iblk)
    grnd_hflx      => tile%mem%soil %energy_flux    %grnd_hflx           (ics:ice,    iblk)
    grnd_hflx_old  => tile%mem%soil %energy_flux    %grnd_hflx_old       (ics:ice,    iblk)
    oro_stddev     => tile%mem%srf  %phys_state     %oro_stddev          (ics:ice,    iblk)
    runoff_acc     => tile%mem%soil %water_flux     %runoff              (ics:ice,    iblk)
    drainage_acc   => tile%mem%soil %water_flux     %drainage            (ics:ice,    iblk)
    ws             => tile%mem%soil %water_state    %ws                  (ics:ice,    iblk)
    ws_l           => tile%mem%soil %water_state    %ws_l                (ics:ice, :, iblk)
    ws_fc_l        => tile%mem%soil %water_state    %ws_fc_l             (ics:ice, :, iblk)
    ws_pwp_l       => tile%mem%soil %water_state    %ws_pwp_l            (ics:ice, :, iblk)
    soil_depth     => tile%mem%soil %phys_state     %soil_depth          (ics:ice,    iblk)
    ws_root        => tile%mem%veg  %water_state    %ws_root             (ics:ice,    iblk)
    ws_fc_root     => tile%mem%veg  %water_state    %ws_fc_root          (ics:ice,    iblk)
    root_depth     => tile%mem%veg  %phys_state     %root_depth          (ics:ice,    iblk)
    shc            => tile%mem%soil %phys_state     %shc                 (ics:ice,    iblk)
    mpot           => tile%mem%soil %phys_state     %mpot                (ics:ice,    iblk)
    bclapp         => tile%mem%soil %phys_state     %bclapp              (ics:ice,    iblk)
    vpor           => tile%mem%soil %phys_state     %vpor                (ics:ice,    iblk)
    psdi           => tile%mem%soil %phys_state     %psdi                (ics:ice,    iblk)
    overflow       => tile%mem%soil %water_state    %overflow            (ics:ice,    iblk)
    glacier        => tile%mem%srf  %phys_state     %glacier             (ics:ice, iblk)
    wsn_glacier    => tile%mem%srf  %water_state    %wsn_glacier         (ics:ice, iblk)

    is_present(:)    = .TRUE.
    is_glacier(:)    = ( glacier(:) > 0.5_wp )

    IF (ALL(wsn_srf(:) < 0._wp)) THEN          ! Start of experiment
      lstart = .TRUE.
      wsn_srf(:) = 0._wp
    ELSE
      lstart = .FALSE.
    END IF

    ! Maximum capacity of skin reservoir (soil + canopy)
    wsr_max(:) = model%config%soil%wsr_capacity * (1._wp + lai(:) * fract_fpc_max(:))

    IF (lstart) THEN
      evapotrans_no_snow_skin(:) = 0._wp
      glacier_precip_minus_evap(:) = 0._wp
      snow_acc(:) = 0._wp
      snow_melt_acc(:) = 0._wp
      glacier_runoff_acc(:) = 0._wp
      melt_water_excess(:) = 0._wp
      tte_corr(:) = 0._wp
      snow_melt(:) = 0._wp
    END IF

    CALL update_surf_down(nc, dtime, lstart, &
      & q_air(:), t_srf_unfilt(:), wind_10m(:), t_air(:), wsr_soil(:), wfract_srf(:), &
      & wsr_max(:), wsn_soil(:), sfract_srf(:), wsn_can(:), wsn_glacier(:), &
      & hcap_grnd(:), evapotrans(:), evapopot(:), evapotrans_no_snow_skin(:), rain(:),  snow(:), &
      & is_present(:), is_glacier(:), glacier_precip_minus_evap(:), &
!!$      & snow_acc(:), snow_melt_acc(:), glacier_runoff_acc(:), &
      melt_water_excess(:), &
      & tte_corr(:), snow_melt(:) &
      & )

    ! Save old ground heat flux and heat capacity before updating
    ! The old flux and capacity (which have been used in the implicit calculation of the new surface temperature) are
    ! given back to the atmosphere so that the surface energy balance can be computed there
    grnd_hflx_old(:) = grnd_hflx(:)
    hcap_grnd_old(:) = hcap_grnd(:)

    CALL update_soiltemp(nc, model%config%soil%nsoil_water, dtime, lstart, &
      & t_srf_unfilt(:), wsn_soil(:), hdiff(:), vhcap(:), &
      & t_soil_acoef(:,:), t_soil_bcoef(:,:), t_soil(:,:), &
      & hcap_grnd(:), grnd_hflx(:), &
      & is_present(:), is_glacier(:))

print*, 'vg-in: runoff_acc: ', MINVAL(runoff_acc), MAXVAL(runoff_acc), SUM(runoff_acc)
!vg:  runoff_hd and drainage_hd no longer needed
    CALL update_surf_up(nc, dtime, model%config%soil%nsoil_water, &
      & wsr_soil(:), wfract_srf(:), wsr_max(:), &
      & ws_root(:), ws_fc_root(:), &
      & t_soil(:,1), sfract_srf(:), oro_stddev(:), &
      & evapotrans(:), evapopot(:), evapotrans_no_snow_skin(:), &
      & rain(:), is_present(:), is_glacier(:), &
      & runoff_acc(:), drainage_acc(:), runoff_hd(:), drainage_hd(:), &
      & melt_water_excess(:), skin_res_evap(:), &
      & model%config%soil%dz_water(1:model%config%soil%nsoil_water), &
      & ws_l(:,:), ws_fc_l(:,:), ws_pwp_l(:,:), root_depth(:), soil_depth(:), &
      & shc(:), mpot(:), bclapp(:), vpor(:), &
      & trans(:), reduced_evap(:), &
!!$      & SPREAD(0.318517_wp, NCOPIES=nc, DIM=1), & !Temporary
!!$      & SPREAD(0.185245_wp, NCOPIES=nc, DIM=1), & !Temporary
      & psdi(:), &
      & overflow(:), &
      & jllog, &
      & l_tpe_open, &
      & l_tpe_closed &
      & )
print*, 'vg-out: runoff_acc: ', MINVAL(runoff_acc), MAXVAL(runoff_acc), SUM(runoff_acc)

    ws(:) = SUM(ws_l(:,:), DIM=2)

    wsn_srf(:) = wsn_soil(:) + wsn_can(:)
    wsr_srf(:) = wsr_soil(:) ! + wsr_can(:)

  END SUBROUTINE update_soil_1

  SUBROUTINE update_soil_2( tile, iblk, ics, ice, nc)


    TYPE(t_jsb_tile),  POINTER    :: tile
    INTEGER,           INTENT(in) :: iblk, ics, ice, nc

    ! Local variables
    !
    TYPE(t_jsb_model), POINTER    :: model

    ! Pointers to variables in memory
    REAL(wp), POINTER ::       &
      & fract          (:)

    ! Locally allocated vectors
    !
!!$    REAL(wp) :: &

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':update_soil_2'

    IF (one_of('soil', tile%components) < 1) RETURN

    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting update')

    model => get_model(tile%owner_model_id)

    ! Get reference to variables for current block
    !
    fract                => tile%mem%fract                               (ics:ice,    iblk)
    IF (ALL(fract(:) < EPSILON(1._wp))) RETURN




  END SUBROUTINE update_soil_2

#endif
END MODULE mo_soil_interface
