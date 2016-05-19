!>
!! @brief Contains the interfaces to the soil processes
!!
!! @author
!!  Reiner Schnur, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Reiner Schnur (2013-04-04)
!! added one layer soil model by Rainer Schneck (2016-05-19)
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

  dsl4jsb_Use_processes SOIL, SRF, VEG, ATM2LAND
  dsl4jsb_Use_config(SOIL)
  USE mo_atmland_memory_class, ONLY: t_atm2land_memory
!!$  dsl4jsb_Use_memory(ATM2LAND)
  dsl4jsb_Use_memory(SOIL)
  dsl4jsb_Use_memory(SRF)
  dsl4jsb_Use_memory(VEG)

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
    dsl4jsb_Def_config(SOIL)
    dsl4jsb_Def_memory(SOIL)
    dsl4jsb_Def_memory(SRF)
    dsl4jsb_Def_memory(VEG)
    dsl4jsb_Def_memory(ATM2LAND)

    ! Pointers to variables in memory
    dsl4jsb_Real2D_onChunk :: fract
    dsl4jsb_Real2D_onChunk :: t_srf_unfilt
    dsl4jsb_Real2D_onChunk :: tte_corr
    dsl4jsb_Real2D_onChunk :: t_air
    dsl4jsb_Real2D_onChunk :: q_air
    dsl4jsb_Real2D_onChunk :: wind_10m
    dsl4jsb_Real2D_onChunk :: rain
    dsl4jsb_Real2D_onChunk :: snow
    dsl4jsb_Real2D_onChunk :: wsr_soil
    dsl4jsb_Real2D_onChunk :: wsr_srf
    dsl4jsb_Real2D_onChunk :: wfract_srf
    dsl4jsb_Real2D_onChunk :: wsn_soil
    dsl4jsb_Real2D_onChunk :: wsn_srf
    dsl4jsb_Real2D_onChunk :: sfract_srf
    dsl4jsb_Real2D_onChunk :: wsn_can
    dsl4jsb_Real2D_onChunk :: lai
    dsl4jsb_Real2D_onChunk :: fract_fpc_max
    dsl4jsb_Real2D_onChunk :: hcap_grnd
    dsl4jsb_Real2D_onChunk :: hcap_grnd_old
    dsl4jsb_Real2D_onChunk :: evapotrans
    dsl4jsb_Real2D_onChunk :: trans
    dsl4jsb_Real2D_onChunk :: evapopot
    dsl4jsb_Real2D_onChunk :: vhcap
    dsl4jsb_Real2D_onChunk :: hdiff
    dsl4jsb_Real3D_onChunk :: t_soil_acoef
    dsl4jsb_Real3D_onChunk :: t_soil_bcoef
    dsl4jsb_Real3D_onChunk :: t_soil
    dsl4jsb_Real2D_onChunk :: grnd_hflx
    dsl4jsb_Real2D_onChunk :: grnd_hflx_old
    dsl4jsb_Real2D_onChunk :: oro_stddev
    dsl4jsb_Real2D_onChunk :: runoff
    dsl4jsb_Real2D_onChunk :: drainage
    dsl4jsb_Real2D_onChunk :: pme_glacier
    dsl4jsb_Real2D_onChunk :: ws
    dsl4jsb_Real3D_onChunk :: ws_l
    dsl4jsb_Real3D_onChunk :: ws_fc_l
    dsl4jsb_Real3D_onChunk :: ws_pwp_l
    dsl4jsb_Real2D_onChunk :: soil_depth
    dsl4jsb_Real2D_onChunk :: ws_root
    dsl4jsb_Real2D_onChunk :: ws_fc_root
    dsl4jsb_Real2D_onChunk :: root_depth
    dsl4jsb_Real2D_onChunk :: shc
    dsl4jsb_Real2D_onChunk :: mpot
    dsl4jsb_Real2D_onChunk :: bclapp
    dsl4jsb_Real2D_onChunk :: vpor
    dsl4jsb_Real2D_onChunk :: psdi
    dsl4jsb_Real2D_onChunk :: overflow
    dsl4jsb_Real2D_onChunk :: glacier
    dsl4jsb_Real2D_onChunk :: wsn_glacier


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
    dsl4jsb_Get_config(SOIL)

    IF (model%config%usecase == 'jsbach_lite_tpe_open')   l_tpe_open = .TRUE.
    IF (model%config%usecase == 'jsbach_lite_tpe_closed') l_tpe_closed = .TRUE.

    ! Get reference to variables for current block
    !
    fract                => tile%fract                               (ics:ice,    iblk)
    IF (ALL(fract(:) < EPSILON(1._wp))) RETURN

    dsl4jsb_Get_memory(SOIL)
    dsl4jsb_Get_memory(SRF)
    dsl4jsb_Get_memory(VEG)
    dsl4jsb_Get_memory(ATM2LAND)

    dsl4jsb_Get_var2D_onChunk(ATM2LAND, t_air,    IN)
    dsl4jsb_Get_var2D_onChunk(ATM2LAND, q_air,    IN)
    dsl4jsb_Get_var2D_onChunk(ATM2LAND, wind_10m, IN)
    dsl4jsb_Get_var2D_onChunk(ATM2LAND, rain,     IN)
    dsl4jsb_Get_var2D_onChunk(ATM2LAND, snow,     IN)

    dsl4jsb_Get_var2D_onChunk(SOIL, wsr_soil,      OUT)
    dsl4jsb_Get_var2D_onChunk(SOIL, wsn_soil,      OUT)
    dsl4jsb_Get_var2D_onChunk(SOIL, hcap_grnd,     OUT)
    dsl4jsb_Get_var2D_onChunk(SOIL, hcap_grnd_old, OUT)
    dsl4jsb_Get_var2D_onChunk(SOIL, vhcap,         OUT)
    dsl4jsb_Get_var2D_onChunk(SOIL, hdiff,         OUT)
    dsl4jsb_Get_var3D_onChunk(SOIL, t_soil_acoef,  OUT)
    dsl4jsb_Get_var3D_onChunk(SOIL, t_soil_bcoef,  OUT)
    dsl4jsb_Get_var3D_onChunk(SOIL, t_soil,        OUT)
    dsl4jsb_Get_var2D_onChunk(SOIL, grnd_hflx,     OUT)
    dsl4jsb_Get_var2D_onChunk(SOIL, grnd_hflx_old, OUT)
    dsl4jsb_Get_var2D_onChunk(SOIL, runoff,        OUT)
    dsl4jsb_Get_var2D_onChunk(SOIL, drainage,      OUT)
    dsl4jsb_Get_var2D_onChunk(SOIL, pme_glacier,   OUT)
    dsl4jsb_Get_var2D_onChunk(SOIL, ws,            OUT)
    dsl4jsb_Get_var3D_onChunk(SOIL, ws_l,          OUT)
    dsl4jsb_Get_var3D_onChunk(SOIL, ws_fc_l,       OUT)
    dsl4jsb_Get_var3D_onChunk(SOIL, ws_pwp_l,      OUT)
    dsl4jsb_Get_var2D_onChunk(SOIL, soil_depth,    OUT)
    dsl4jsb_Get_var2D_onChunk(SOIL, shc,           IN)
    dsl4jsb_Get_var2D_onChunk(SOIL, mpot,          IN)
    dsl4jsb_Get_var2D_onChunk(SOIL, bclapp,        IN)
    dsl4jsb_Get_var2D_onChunk(SOIL, vpor,          IN)
    dsl4jsb_Get_var2D_onChunk(SOIL, psdi,          IN)
    dsl4jsb_Get_var2D_onChunk(SOIL, overflow,      OUT)

    dsl4jsb_Get_var2D_onChunk(SRF, t_srf_unfilt, OUT)
    dsl4jsb_Get_var2D_onChunk(SRF, tte_corr,     OUT)
    dsl4jsb_Get_var2D_onChunk(SRF, wsr_srf,      OUT)
    dsl4jsb_Get_var2D_onChunk(SRF, wfract_srf,   OUT)
    dsl4jsb_Get_var2D_onChunk(SRF, wsn_srf,      OUT)
    dsl4jsb_Get_var2D_onChunk(SRF, sfract_srf,   OUT)
    dsl4jsb_Get_var2D_onChunk(SRF, evapotrans,   OUT)
    dsl4jsb_Get_var2D_onChunk(SRF, trans,        OUT)
    dsl4jsb_Get_var2D_onChunk(SRF, evapopot,     OUT)
    dsl4jsb_Get_var2D_onChunk(SRF, oro_stddev,   IN)
    dsl4jsb_Get_var2D_onChunk(SRF, glacier,      OUT)
    dsl4jsb_Get_var2D_onChunk(SRF, wsn_glacier,  OUT)

    dsl4jsb_Get_var2D_onChunk(VEG, wsn_can,       OUT)
    dsl4jsb_Get_var2D_onChunk(VEG, lai,           IN)
    dsl4jsb_Get_var2D_onChunk(VEG, fract_fpc_max, OUT)
    dsl4jsb_Get_var2D_onChunk(VEG, ws_root,       OUT)
    dsl4jsb_Get_var2D_onChunk(VEG, ws_fc_root,    OUT)
    dsl4jsb_Get_var2D_onChunk(VEG, root_depth,    OUT)

    is_present(:)    = .TRUE.
    is_glacier(:)    = ( glacier(:) > 0.5_wp )

    IF (ALL(wsn_srf(:) < 0._wp)) THEN          ! Start of experiment
      lstart = .TRUE.
      wsn_srf(:) = 0._wp
    ELSE
      lstart = .FALSE.
    END IF

    ! Maximum capacity of skin reservoir (soil + canopy)
    wsr_max(:) = dsl4jsb_Config(SOIL)%wsr_capacity * (1._wp + lai(:) * fract_fpc_max(:))

    IF (lstart) THEN
      evapotrans_no_snow_skin(:) = 0._wp
      snow_acc(:) = 0._wp
      snow_melt_acc(:) = 0._wp
      glacier_runoff_acc(:) = 0._wp
      melt_water_excess(:) = 0._wp
      tte_corr(:) = 0._wp
      snow_melt(:) = 0._wp
    END IF

    CALL update_surf_down( & ! R: Wasserfluxe in die Surface und Schneeschmelze/Sublimation
      & nc,              & ! In
      & dtime,           & ! In
      & lstart,          & ! In
      & q_air(:),        & ! In
      & t_srf_unfilt(:), & ! InOut
      & wind_10m(:),     & ! In
      & t_air(:),        & ! In
      & wsr_soil(:),     & ! InOut
      & wfract_srf(:),   & ! In
      & wsr_max(:),      & ! In
      & wsn_soil(:),     & ! InOut
      & sfract_srf(:),   & ! In
      & wsn_can(:),      & ! InOut
      & wsn_glacier(:),  & ! InOut
      & hcap_grnd(:),    & ! In
      & evapotrans(:),   & ! In
      & evapopot(:),     & ! In
      & evapotrans_no_snow_skin(:), & ! Out
      & rain(:),         & ! In
      & snow(:),         & ! In
      & is_present(:),   & ! In
      & is_glacier(:),   & ! In
      & pme_glacier(:),  & ! Out
!!$      & snow_acc(:), snow_melt_acc(:), glacier_runoff_acc(:), &
      & melt_water_excess(:), & ! Out
      & tte_corr(:),     & ! Out
      & snow_melt(:)     & ! Out
      & )

    ! Save old ground heat flux and heat capacity before updating
    ! The old flux and capacity (which have been used in the implicit calculation of the new surface temperature) are
    ! given back to the atmosphere so that the surface energy balance can be computed there
    grnd_hflx_old(:) = grnd_hflx(:)
    hcap_grnd_old(:) = hcap_grnd(:)

    CALL update_soiltemp(nc, dsl4jsb_Config(SOIL)%nsoil_water, dtime, lstart, &
      & t_srf_unfilt(:), wsn_soil(:), hdiff(:), vhcap(:), &
      & t_soil_acoef(:,:), t_soil_bcoef(:,:), t_soil(:,:), &
      & hcap_grnd(:), grnd_hflx(:), &
      & is_present(:), is_glacier(:))
      
    if (dsl4jsb_Config(SOIL)%nsoil_water > 1) THEN
      CALL update_surf_up(&
        & nc, dtime, dsl4jsb_Config(SOIL)%nsoil_water, & ! In (all)
        & wsr_soil(:),    & ! InOut. Water content in skin reservoir
!                           ! R: Note, this is possibly not fully correct, it was an average value over the tiles in JSBACH3 (SkinReservoire_tile)
!                           !    which was especially created and handled (avoids that it is the same as in the 1 layer soil model
!                           !    as it is here in JSBACH4).
        & wfract_srf(:),  & ! In
                            ! R: Note, this is possibly not fully correct, it was an average value over the tiles in JSBACH3 (WetSkin_frac)
!                           !    which was especially created (avoids that it is the same as in the 1 layer soil model
!                           !    as it is here in JSBACH4).
        & wsr_max(:),     & ! In
                            ! R: Note,  this is possibly not fully correct, it was an average value over the tiles in JSBACH3 (MaxSkinResCap)
!                           !    which was especially created (avoids that it is the same as in the 1 layer soil model
!                           !    as it is here in JSBACH4). Intend "in" when used in update_surf_down.
        & ws_root(:),     & ! InOut
        & ws_fc_root(:),  & ! In
        & t_soil(:,1),    & ! In
        & sfract_srf(:),  & ! In
        & oro_stddev(:),  & ! In
        & evapotrans(:),  & ! In. New in JSBACH4, but not used in the moment
        & evapopot(:),    & ! In
        & evapotrans_no_snow_skin(:), & ! InOut
        & rain(:),        & ! In
        & is_present(:),  & ! In
        & is_glacier(:),  & ! In
        & runoff(:),      & ! InOut
        & drainage(:),    & ! InOut
!vg     & runoff_hd(:),   &
!vg     & drainage_hd(:), &
        & melt_water_excess(:), & ! In
        & skin_res_evap(:), & ! Out
        & dsl4jsb_Config(SOIL)%dz_water(1:dsl4jsb_Config(SOIL)%nsoil_water), &  ! In
        & ws_l(:,:),      & ! InOut
        & ws_fc_l(:,:),   & ! In
        & ws_pwp_l(:,:),  & ! In
        & root_depth(:),  & ! In
        & soil_depth(:),  & ! In
        & shc(:),         & ! In
        & mpot(:),        & ! In
        & bclapp(:),      & ! In
        & vpor(:),        & ! In
        & trans(:),       & ! In
        & reduced_evap(:), & ! Out
!!$      & SPREAD(0.318517_wp, NCOPIES=nc, DIM=1), & !Temporary
!!$      & SPREAD(0.185245_wp, NCOPIES=nc, DIM=1), & !Temporary
        & psdi(:),        & ! In
        & overflow(:),    & ! InOut
        & jllog,          & ! Inout
        & l_tpe_open,     & ! In
        & l_tpe_closed    & ! In
        & )    
    ELSE
      CALL update_surf_up(&
        & nc, dtime, dsl4jsb_Config(SOIL)%nsoil_water, & ! In (all)
        & wsr_soil(:),    & ! InOut. Water content in skin reservoir
!                           ! R: Note, this is correct as it was NOT an average value over the tiles in JSBACH3 but soil%skin_reservoir
        & wfract_srf(:),  & ! In. R: Note, this is correct as it was NOT an average value over the tiles in JSBACH3 but wet_skin_fract
        & wsr_max(:),     & ! In. R: Note, this is correct as it was NOT an average value over the tiles in JSBACH3 but skin_reservoir_max
        & ws_root(:),     & ! InOut
        & ws_fc_root(:),  & ! In
        & t_soil(:,1),    & ! In
        & sfract_srf(:),  & ! In
        & oro_stddev(:),  & ! In
        & evapotrans(:),  & ! In
        & evapopot(:),    & ! In
        & evapotrans_no_snow_skin(:), & ! InOut
        & rain(:),        & ! In
        & is_present(:),  & ! In
        & is_glacier(:),  & ! In
        & runoff(:),      & ! InOut
        & drainage(:),    & ! InOut
!vg     & runoff_hd(:),   &
!vg     & drainage_hd(:), &
        & melt_water_excess(:), & ! In
        & skin_res_evap(:) & ! Out
        & )
    END IF
    
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
