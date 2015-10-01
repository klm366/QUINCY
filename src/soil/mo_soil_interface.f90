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

#include "dsl_jsbach.inc"

  USE mo_kind,            ONLY: wp
  USE mo_exception,       ONLY: message, finish

  USE mo_jsb_model_class, ONLY: t_jsb_model
  USE mo_jsb_class,       ONLY: Get_model
  USE mo_jsb_tile_class,  ONLY: t_jsb_tile
  USE mo_jsb_control,     ONLY: debug_on

  Use_processes SOIL, SRF, VEG, ATM2LAND
  Use_config(soil)
  USE mo_atmland_memory_class, ONLY: t_atm2land_memory
!!$  Use_memory(atm2land)
  Use_memory(soil)
  Use_memory(srf)
  Use_memory(veg)

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: update_soil_1, update_soil_2

  CHARACTER(len=*), PARAMETER :: modname = 'mo_soil_interface'

CONTAINS

  ! TBD: At the moment, the hydrology for soil and canopy is handled together in the soil process (as in echam6).
  ! This will be changed later.
  SUBROUTINE update_soil_1( tile, iblk, ics, ice, nc, dtime, steplen)

    USE mo_soil_process, ONLY: update_surf_down, update_soiltemp, update_surf_up

    TYPE(t_jsb_tile),  POINTER    :: tile
    INTEGER,           INTENT(in) :: iblk, ics, ice, nc
    REAL(wp),          INTENT(in) :: dtime, steplen

    ! Local variables
    !
    TYPE(t_jsb_model), POINTER    :: model
    Def_config(soil_conf, soil)
    Def_memory(soil_mem, soil)
    Def_memory(srf_mem, srf)
    Def_memory(veg_mem, veg)
    Def_memory(atm2land_mem, atm2land)

    ! Pointers to variables in memory
    Def_real2d_onChunk(fract)
    Def_real2d_onChunk(t_srf_unfilt)
    Def_real2d_onChunk(tte_corr)
    Def_real2d_onChunk(t_air)
    Def_real2d_onChunk(q_air)
    Def_real2d_onChunk(wind_10m)
    Def_real2d_onChunk(rain)
    Def_real2d_onChunk(snow)
    Def_real2d_onChunk(wsr_soil)
    Def_real2d_onChunk(wsr_srf)
    Def_real2d_onChunk(wfract_srf)
    Def_real2d_onChunk(wsn_soil)
    Def_real2d_onChunk(wsn_srf)
    Def_real2d_onChunk(sfract_srf)
    Def_real2d_onChunk(wsn_can)
    Def_real2d_onChunk(lai)
    Def_real2d_onChunk(fract_fpc_max)
    Def_real2d_onChunk(hcap_grnd)
    Def_real2d_onChunk(hcap_grnd_old)
    Def_real2d_onChunk(evapotrans)
    Def_real2d_onChunk(trans)
    Def_real2d_onChunk(evapopot)
    Def_real2d_onChunk(vhcap)
    Def_real2d_onChunk(hdiff)
    Def_real3d_onChunk(t_soil_acoef)
    Def_real3d_onChunk(t_soil_bcoef)
    Def_real3d_onChunk(t_soil)
    Def_real2d_onChunk(grnd_hflx)
    Def_real2d_onChunk(grnd_hflx_old)
    Def_real2d_onChunk(oro_stddev)
    Def_real2d_onChunk(runoff)
    Def_real2d_onChunk(drainage)
    Def_real2d_onChunk(pme_glacier)
    Def_real2d_onChunk(ws)
    Def_real3d_onChunk(ws_l)
    Def_real3d_onChunk(ws_fc_l)
    Def_real3d_onChunk(ws_pwp_l)
    Def_real2d_onChunk(soil_depth)
    Def_real2d_onChunk(ws_root)
    Def_real2d_onChunk(ws_fc_root)
    Def_real2d_onChunk(root_depth)
    Def_real2d_onChunk(shc)
    Def_real2d_onChunk(mpot)
    Def_real2d_onChunk(bclapp)
    Def_real2d_onChunk(vpor)
    Def_real2d_onChunk(psdi)
    Def_real2d_onChunk(overflow)
    Def_real2d_onChunk(glacier)
    Def_real2d_onChunk(wsn_glacier)


    ! Locally allocated vectors
    !
    REAL(wp) :: &
      wsr_max(nc), &
      evapotrans_no_snow_skin(nc), &
      snow_acc(nc), &
      snow_melt_acc(nc), &
      glacier_runoff_acc(nc), &
      melt_water_excess(nc), &
      snow_melt(nc), &
!vg       runoff_hd(nc), &
!vg       drainage_hd(nc), &
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

    IF (.NOT. tile%Is_process_active(SOIL)) RETURN

    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting update')

    model => Get_model(tile%owner_model_id)

!!$    soil_conf => Get_soil_config(model)
    Get_config(soil_conf, SOIL)

    IF (model%config%usecase == 'jsbach_lite_tpe_open')   l_tpe_open = .TRUE.
    IF (model%config%usecase == 'jsbach_lite_tpe_closed') l_tpe_closed = .TRUE.

    ! Get reference to variables for current block
    !
    fract                => tile%fract                               (ics:ice,    iblk)
    IF (ALL(fract(:) < EPSILON(1._wp))) RETURN

    Get_memory(soil_mem, SOIL)
    Get_memory(srf_mem, SRF)
    Get_memory(veg_mem, VEG)
    Get_memory(atm2land_mem, ATM2LAND)

    t_srf_unfilt   => srf_mem  %t_srf_unfilt%ptr        (ics:ice,    iblk)
    tte_corr       => srf_mem  %tte_corr%ptr            (ics:ice,    iblk)
    t_air          => atm2land_mem       %t_air%ptr               (ics:ice,    iblk)
    q_air          => atm2land_mem       %q_air%ptr               (ics:ice,    iblk)
    wind_10m       => atm2land_mem       %wind_10m%ptr            (ics:ice,    iblk)
    rain           => atm2land_mem       %rain%ptr                (ics:ice,    iblk)
    snow           => atm2land_mem       %snow%ptr                (ics:ice,    iblk)
    wsr_soil       => soil_mem %wsr_soil%ptr            (ics:ice,    iblk)
    wsr_srf        => srf_mem  %wsr_srf%ptr             (ics:ice,    iblk)
    wfract_srf     => srf_mem  %wfract_srf%ptr          (ics:ice,    iblk)
    wsn_soil       => soil_mem %wsn_soil%ptr            (ics:ice,    iblk)
    wsn_srf        => srf_mem  %wsn_srf%ptr             (ics:ice,    iblk)
    sfract_srf     => srf_mem  %sfract_srf%ptr          (ics:ice,    iblk)
    wsn_can        => veg_mem  %wsn_can%ptr             (ics:ice,    iblk)
    lai            => veg_mem  %lai%ptr                 (ics:ice,    iblk)
    fract_fpc_max  => veg_mem  %fract_fpc_max%ptr       (ics:ice,    iblk)
    hcap_grnd      => soil_mem %hcap_grnd%ptr           (ics:ice,    iblk)
    hcap_grnd_old  => soil_mem %hcap_grnd_old%ptr       (ics:ice,    iblk)
    evapotrans     => srf_mem  %evapotrans%ptr          (ics:ice,    iblk)
    trans          => srf_mem  %trans%ptr               (ics:ice,    iblk)
    evapopot       => srf_mem  %evapopot%ptr            (ics:ice,    iblk)
    vhcap          => soil_mem %vhcap%ptr               (ics:ice,    iblk)
    hdiff          => soil_mem %hdiff%ptr               (ics:ice,    iblk)
    t_soil_acoef   => soil_mem %t_soil_acoef%ptr        (ics:ice, :, iblk)
    t_soil_bcoef   => soil_mem %t_soil_bcoef%ptr        (ics:ice, :, iblk)
    t_soil         => soil_mem %t_soil%ptr              (ics:ice, :, iblk)
    grnd_hflx      => soil_mem %grnd_hflx%ptr           (ics:ice,    iblk)
    grnd_hflx_old  => soil_mem %grnd_hflx_old%ptr       (ics:ice,    iblk)
    oro_stddev     => srf_mem  %oro_stddev%ptr          (ics:ice,    iblk)
    runoff         => soil_mem %runoff%ptr              (ics:ice,    iblk)
    drainage       => soil_mem %drainage%ptr            (ics:ice,    iblk)
    pme_glacier    => soil_mem %pme_glacier%ptr         (ics:ice,    iblk)
    ws             => soil_mem %ws%ptr                  (ics:ice,    iblk)
    ws_l           => soil_mem %ws_l%ptr                (ics:ice, :, iblk)
    ws_fc_l        => soil_mem %ws_fc_l%ptr             (ics:ice, :, iblk)
    ws_pwp_l       => soil_mem %ws_pwp_l%ptr            (ics:ice, :, iblk)
    soil_depth     => soil_mem %soil_depth%ptr          (ics:ice,    iblk)
    ws_root        => veg_mem  %ws_root%ptr             (ics:ice,    iblk)
    ws_fc_root     => veg_mem  %ws_fc_root%ptr          (ics:ice,    iblk)
    root_depth     => veg_mem  %root_depth%ptr          (ics:ice,    iblk)
    shc            => soil_mem %shc%ptr                 (ics:ice,    iblk)
    mpot           => soil_mem %mpot%ptr                (ics:ice,    iblk)
    bclapp         => soil_mem %bclapp%ptr              (ics:ice,    iblk)
    vpor           => soil_mem %vpor%ptr                (ics:ice,    iblk)
    psdi           => soil_mem %psdi%ptr                (ics:ice,    iblk)
    overflow       => soil_mem %overflow%ptr            (ics:ice,    iblk)
    glacier        => srf_mem  %glacier%ptr             (ics:ice, iblk)
    wsn_glacier    => srf_mem  %wsn_glacier%ptr         (ics:ice, iblk)

    is_present(:)    = .TRUE.
    is_glacier(:)    = ( glacier(:) > 0.5_wp )

    IF (ALL(wsn_srf(:) < 0._wp)) THEN          ! Start of experiment
      lstart = .TRUE.
      wsn_srf(:) = 0._wp
    ELSE
      lstart = .FALSE.
    END IF

    ! Maximum capacity of skin reservoir (soil + canopy)
    wsr_max(:) = soil_conf%wsr_capacity * (1._wp + lai(:) * fract_fpc_max(:))

    IF (lstart) THEN
      evapotrans_no_snow_skin(:) = 0._wp
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
      & is_present(:), is_glacier(:), pme_glacier(:), &
!!$      & snow_acc(:), snow_melt_acc(:), glacier_runoff_acc(:), &
      melt_water_excess(:), &
      & tte_corr(:), snow_melt(:) &
      & )

    ! Save old ground heat flux and heat capacity before updating
    ! The old flux and capacity (which have been used in the implicit calculation of the new surface temperature) are
    ! given back to the atmosphere so that the surface energy balance can be computed there
    grnd_hflx_old(:) = grnd_hflx(:)
    hcap_grnd_old(:) = hcap_grnd(:)

    CALL update_soiltemp(nc, soil_conf%nsoil_water, dtime, lstart, &
      & t_srf_unfilt(:), wsn_soil(:), hdiff(:), vhcap(:), &
      & t_soil_acoef(:,:), t_soil_bcoef(:,:), t_soil(:,:), &
      & hcap_grnd(:), grnd_hflx(:), &
      & is_present(:), is_glacier(:))

!vg:  runoff_hd and drainage_hd no longer needed
!vg:  accumulation does not make sense (no reset)
    CALL update_surf_up(nc, dtime, soil_conf%nsoil_water, &
      & wsr_soil(:), wfract_srf(:), wsr_max(:), &
      & ws_root(:), ws_fc_root(:), &
      & t_soil(:,1), sfract_srf(:), oro_stddev(:), &
      & evapotrans(:), evapopot(:), evapotrans_no_snow_skin(:), &
      & rain(:), is_present(:), is_glacier(:), &
      & runoff(:), drainage(:), & !vg runoff_hd(:), drainage_hd(:), &
      & melt_water_excess(:), skin_res_evap(:), &
      & soil_conf%dz_water(1:soil_conf%nsoil_water), &
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

    IF (.NOT. tile%Is_process_active(SOIL)) RETURN

    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting update')

    model => Get_model(tile%owner_model_id)

    ! Get reference to variables for current block
    !
    fract                => tile%fract                               (ics:ice,    iblk)
    IF (ALL(fract(:) < EPSILON(1._wp))) RETURN




  END SUBROUTINE update_soil_2

#endif
END MODULE mo_soil_interface
