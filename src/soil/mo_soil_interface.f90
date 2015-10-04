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

  dsl$jsb Use_processes SOIL, SRF, VEG, ATM2LAND
  dsl$jsb Use_config(soil)
  USE mo_atmland_memory_class, ONLY: t_atm2land_memory
!!$  dsl$jsb Use_memory(atm2land)
  dsl$jsb Use_memory(soil)
  dsl$jsb Use_memory(srf)
  dsl$jsb Use_memory(veg)

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
    dsl$jsb Def_config(soil_conf, soil)
    dsl$jsb Def_memory(soil_mem, soil)
    dsl$jsb Def_memory(srf_mem, srf)
    dsl$jsb Def_memory(veg_mem, veg)
    dsl$jsb Def_memory(atm2land_mem, atm2land)

    ! Pointers to variables in memory
    dsl$jsb onChunk_Real2D :: fract
    dsl$jsb onChunk_Real2D :: t_srf_unfilt
    dsl$jsb onChunk_Real2D :: tte_corr
    dsl$jsb onChunk_Real2D :: t_air
    dsl$jsb onChunk_Real2D :: q_air
    dsl$jsb onChunk_Real2D :: wind_10m
    dsl$jsb onChunk_Real2D :: rain
    dsl$jsb onChunk_Real2D :: snow
    dsl$jsb onChunk_Real2D :: wsr_soil
    dsl$jsb onChunk_Real2D :: wsr_srf
    dsl$jsb onChunk_Real2D :: wfract_srf
    dsl$jsb onChunk_Real2D :: wsn_soil
    dsl$jsb onChunk_Real2D :: wsn_srf
    dsl$jsb onChunk_Real2D :: sfract_srf
    dsl$jsb onChunk_Real2D :: wsn_can
    dsl$jsb onChunk_Real2D :: lai
    dsl$jsb onChunk_Real2D :: fract_fpc_max
    dsl$jsb onChunk_Real2D :: hcap_grnd
    dsl$jsb onChunk_Real2D :: hcap_grnd_old
    dsl$jsb onChunk_Real2D :: evapotrans
    dsl$jsb onChunk_Real2D :: trans
    dsl$jsb onChunk_Real2D :: evapopot
    dsl$jsb onChunk_Real2D :: vhcap
    dsl$jsb onChunk_Real2D :: hdiff
    dsl$jsb onChunk_Real3D :: t_soil_acoef
    dsl$jsb onChunk_Real3D :: t_soil_bcoef
    dsl$jsb onChunk_Real3D :: t_soil
    dsl$jsb onChunk_Real2D :: grnd_hflx
    dsl$jsb onChunk_Real2D :: grnd_hflx_old
    dsl$jsb onChunk_Real2D :: oro_stddev
    dsl$jsb onChunk_Real2D :: runoff
    dsl$jsb onChunk_Real2D :: drainage
    dsl$jsb onChunk_Real2D :: pme_glacier
    dsl$jsb onChunk_Real2D :: ws
    dsl$jsb onChunk_Real3D :: ws_l
    dsl$jsb onChunk_Real3D :: ws_fc_l
    dsl$jsb onChunk_Real3D :: ws_pwp_l
    dsl$jsb onChunk_Real2D :: soil_depth
    dsl$jsb onChunk_Real2D :: ws_root
    dsl$jsb onChunk_Real2D :: ws_fc_root
    dsl$jsb onChunk_Real2D :: root_depth
    dsl$jsb onChunk_Real2D :: shc
    dsl$jsb onChunk_Real2D :: mpot
    dsl$jsb onChunk_Real2D :: bclapp
    dsl$jsb onChunk_Real2D :: vpor
    dsl$jsb onChunk_Real2D :: psdi
    dsl$jsb onChunk_Real2D :: overflow
    dsl$jsb onChunk_Real2D :: glacier
    dsl$jsb onChunk_Real2D :: wsn_glacier


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
    dsl$jsb Get_config(soil_conf, SOIL)

    IF (model%config%usecase == 'jsbach_lite_tpe_open')   l_tpe_open = .TRUE.
    IF (model%config%usecase == 'jsbach_lite_tpe_closed') l_tpe_closed = .TRUE.

    ! Get reference to variables for current block
    !
    fract                => tile%fract                               (ics:ice,    iblk)
    IF (ALL(fract(:) < EPSILON(1._wp))) RETURN

    dsl$jsb Get_memory(soil_mem, SOIL)
    dsl$jsb Get_memory(srf_mem, SRF)
    dsl$jsb Get_memory(veg_mem, VEG)
    dsl$jsb Get_memory(atm2land_mem, ATM2LAND)

    dsl$jsb Get_var2D_onChunk( atm2land_mem, t_air,         IN)
    dsl$jsb Get_var2D_onChunk( atm2land_mem, q_air,         IN)
    dsl$jsb Get_var2D_onChunk( atm2land_mem, wind_10m,      IN)
    dsl$jsb Get_var2D_onChunk( atm2land_mem, rain,          IN)
    dsl$jsb Get_var2D_onChunk( atm2land_mem, snow,          IN)

    dsl$jsb Get_var2D_onChunk( soil_mem,     wsr_soil,      OUT)
    dsl$jsb Get_var2D_onChunk( soil_mem,     wsn_soil,      OUT)
    dsl$jsb Get_var2D_onChunk( soil_mem,     hcap_grnd,     OUT)
    dsl$jsb Get_var2D_onChunk( soil_mem,     hcap_grnd_old, OUT)
    dsl$jsb Get_var2D_onChunk( soil_mem,     vhcap,         OUT)
    dsl$jsb Get_var2D_onChunk( soil_mem,     hdiff,         OUT)
    dsl$jsb Get_var3D_onChunk( soil_mem,     t_soil_acoef,  OUT)
    dsl$jsb Get_var3D_onChunk( soil_mem,     t_soil_bcoef,  OUT)
    dsl$jsb Get_var3D_onChunk( soil_mem,     t_soil,        OUT)
    dsl$jsb Get_var2D_onChunk( soil_mem,     grnd_hflx,     OUT)
    dsl$jsb Get_var2D_onChunk( soil_mem,     grnd_hflx_old, OUT)
    dsl$jsb Get_var2D_onChunk( soil_mem,     runoff,        OUT)
    dsl$jsb Get_var2D_onChunk( soil_mem,     drainage,      OUT)
    dsl$jsb Get_var2D_onChunk( soil_mem,     pme_glacier,   OUT)
    dsl$jsb Get_var2D_onChunk( soil_mem,     ws,            OUT)
    dsl$jsb Get_var3D_onChunk( soil_mem,     ws_l,          OUT)
    dsl$jsb Get_var3D_onChunk( soil_mem,     ws_fc_l,       OUT)
    dsl$jsb Get_var3D_onChunk( soil_mem,     ws_pwp_l,      OUT)
    dsl$jsb Get_var2D_onChunk( soil_mem,     soil_depth,    OUT)
    dsl$jsb Get_var2D_onChunk( soil_mem,     shc,           IN)
    dsl$jsb Get_var2D_onChunk( soil_mem,     mpot,          IN)
    dsl$jsb Get_var2D_onChunk( soil_mem,     bclapp,        IN)
    dsl$jsb Get_var2D_onChunk( soil_mem,     vpor,          IN)
    dsl$jsb Get_var2D_onChunk( soil_mem,     psdi,          IN)
    dsl$jsb Get_var2D_onChunk( soil_mem,     overflow,      OUT)

    dsl$jsb Get_var2D_onChunk( srf_mem,      t_srf_unfilt,  OUT)
    dsl$jsb Get_var2D_onChunk( srf_mem,      tte_corr,      OUT)
    dsl$jsb Get_var2D_onChunk( srf_mem,      wsr_srf,       OUT)
    dsl$jsb Get_var2D_onChunk( srf_mem,      wfract_srf,    OUT)
    dsl$jsb Get_var2D_onChunk( srf_mem,      wsn_srf,       OUT)
    dsl$jsb Get_var2D_onChunk( srf_mem,      sfract_srf,    OUT)
    dsl$jsb Get_var2D_onChunk( srf_mem,      evapotrans,    OUT)
    dsl$jsb Get_var2D_onChunk( srf_mem,      trans,         OUT)
    dsl$jsb Get_var2D_onChunk( srf_mem,      evapopot,      OUT)
    dsl$jsb Get_var2D_onChunk( srf_mem,      oro_stddev,    IN)
    dsl$jsb Get_var2D_onChunk( srf_mem,      glacier,       OUT)
    dsl$jsb Get_var2D_onChunk( srf_mem,      wsn_glacier,   OUT)

    dsl$jsb Get_var2D_onChunk( veg_mem,      wsn_can,       OUT)
    dsl$jsb Get_var2D_onChunk( veg_mem,      lai,           OUT)
    dsl$jsb Get_var2D_onChunk( veg_mem,      fract_fpc_max, OUT)
    dsl$jsb Get_var2D_onChunk( veg_mem,      ws_root,       OUT)
    dsl$jsb Get_var2D_onChunk( veg_mem,      ws_fc_root,    OUT)
    dsl$jsb Get_var2D_onChunk( veg_mem,      root_depth,    OUT)

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
