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

  USE mo_jsb_control,     ONLY: debug_on

  USE mo_jsb_tile_class,     ONLY: t_jsb_tile_abstract
  USE mo_jsb_memory_class,   ONLY: t_jsb_memory, t_jsb_memory_p
  USE mo_jsb_config_class,   ONLY: t_jsb_config, t_jsb_config_p
  USE mo_jsb_process_class,  ONLY: t_jsb_process
  USE mo_jsb_task_class,     ONLY: t_jsb_process_task, t_jsb_task_options
  USE mo_jsb_lct_class,      ONLY: GLACIER_TYPE

  dsl4jsb_Use_processes SOIL, SEB, VEG, ATM2LAND, HYDRO
  dsl4jsb_Use_config(SOIL)
  USE mo_atmland_memory_class, ONLY: t_atm2land_memory
!!$  dsl4jsb_Use_memory(ATM2LAND)
  dsl4jsb_Use_memory(SOIL)
  dsl4jsb_Use_memory(SEB)
  dsl4jsb_Use_memory(VEG)
  dsl4jsb_Use_memory(HYDRO)

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: update_soil_1, update_soil_2
  PUBLIC :: t_soil_process

  TYPE, EXTENDS(t_jsb_process) :: t_soil_process
CONTAINS
  PROCEDURE :: Register_tasks => Register_soil_tasks
END TYPE t_soil_process

  CHARACTER(len=*), PARAMETER :: modname = 'mo_soil_interface'

CONTAINS

  SUBROUTINE Register_soil_tasks(this, model_id)

    CLASS(t_soil_process), INTENT(inout) :: this
    INTEGER, INTENT(in) :: model_id

    !CALL this%Register_task(tsk_radiation_surface(name="radiation_surface", process_id=SOIL))

  END SUBROUTINE Register_soil_tasks

  ! TBD: At the moment, the hydrology for soil and canopy is handled together in the soil process (as in echam6).
  ! This will be changed later.
  SUBROUTINE update_soil_1(tile, config, options)

    USE mo_soil_process, ONLY: update_surf_down, update_soiltemp, update_surf_up

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    CLASS(t_jsb_config_p),      INTENT(in)    :: config(:)
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    ! Local variables
    !
    dsl4jsb_Def_config(SOIL)
    dsl4jsb_Def_memory(SOIL)
    dsl4jsb_Def_memory(SEB)
    dsl4jsb_Def_memory(VEG)
    dsl4jsb_Def_memory(ATM2LAND)
    dsl4jsb_Def_memory(HYDRO)

    ! Pointers to variables in memory
    dsl4jsb_Real2D_onChunk :: fract
    dsl4jsb_Real2D_onChunk :: t_unfilt
    dsl4jsb_Real2D_onChunk :: tte_corr
    dsl4jsb_Real2D_onChunk :: t_air
    dsl4jsb_Real2D_onChunk :: q_air
    dsl4jsb_Real2D_onChunk :: wind_10m
    dsl4jsb_Real2D_onChunk :: rain
    dsl4jsb_Real2D_onChunk :: snow
    dsl4jsb_Real2D_onChunk :: w_skin_soil
    dsl4jsb_Real2D_onChunk :: w_skin_srf
    dsl4jsb_Real2D_onChunk :: fract_water_srf
    dsl4jsb_Real2D_onChunk :: w_snow_soil
    dsl4jsb_Real2D_onChunk :: w_snow_srf
    dsl4jsb_Real2D_onChunk :: fract_snow_srf
    dsl4jsb_Real2D_onChunk :: w_snow_can
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
    dsl4jsb_Real2D_onChunk :: w_soil_column
    dsl4jsb_Real3D_onChunk :: w_soil_sl
    dsl4jsb_Real3D_onChunk :: w_soil_fc_sl
    dsl4jsb_Real3D_onChunk :: w_soil_pwp_sl
    dsl4jsb_Real2D_onChunk :: soil_depth
    dsl4jsb_Real2D_onChunk :: w_soil_root
    dsl4jsb_Real2D_onChunk :: w_soil_root_fc
    dsl4jsb_Real2D_onChunk :: root_depth
    dsl4jsb_Real2D_onChunk :: shc
    dsl4jsb_Real2D_onChunk :: mpot
    dsl4jsb_Real2D_onChunk :: bclapp
    dsl4jsb_Real2D_onChunk :: vpor
    dsl4jsb_Real2D_onChunk :: psdi
    dsl4jsb_Real2D_onChunk :: w_soil_overflow
    dsl4jsb_Real2D_onChunk :: w_snow_glacier


    ! Locally allocated vectors
    !
    REAL(wp) :: &
      wsr_max(options%nc), &
      evapotrans_no_snow_skin(options%nc), &
      snow_acc(options%nc), &
      snow_melt_acc(options%nc), &
      glacier_runoff_acc(options%nc), &
      melt_water_excess(options%nc), &
      snow_melt(options%nc), &
!vg       runoff_hd(options%nc), &
!vg       drainage_hd(options%nc), &
      skin_res_evap(options%nc), &
      reduced_evap(options%nc)

    LOGICAL :: &
      is_present(options%nc), &
      is_glacier(options%nc)

    INTEGER :: jllog

    LOGICAL  :: lstart
    LOGICAL  :: l_tpe_open = .FALSE.
    LOGICAL  :: l_tpe_closed = .FALSE.
    INTEGER  :: iblk, ics, ice, nc
    REAL(wp) :: dtime, steplen

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':update_soil_1'

    iblk  = options%iblk
    ics   = options%ics
    ice   = options%ice
    nc    = options%nc
    dtime = options%dtime
    steplen = options%steplen

    IF (.NOT. tile%Is_process_active(SOIL)) RETURN

    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting update')

!!$    soil_conf => Get_soil_config(model)
    dsl4jsb_Get_config(SOIL)

!!$    IF (model%config%usecase == 'jsbach_lite_tpe_open')   l_tpe_open = .TRUE.
!!$    IF (model%config%usecase == 'jsbach_lite_tpe_closed') l_tpe_closed = .TRUE.

    ! Get reference to variables for current block
    !
    fract                => tile%fract                               (ics:ice,    iblk)
    IF (ALL(fract(:) < EPSILON(1._wp))) RETURN

    dsl4jsb_Get_memory(SOIL)
    dsl4jsb_Get_memory(SEB)
    dsl4jsb_Get_memory(VEG)
    dsl4jsb_Get_memory(ATM2LAND)
    dsl4jsb_Get_memory(HYDRO)

    dsl4jsb_Get_var2D_onChunk(ATM2LAND, t_air)         ! IN
    dsl4jsb_Get_var2D_onChunk(ATM2LAND, q_air)         ! IN
    dsl4jsb_Get_var2D_onChunk(ATM2LAND, wind_10m)      ! IN
    dsl4jsb_Get_var2D_onChunk(ATM2LAND, rain)          ! IN
    dsl4jsb_Get_var2D_onChunk(ATM2LAND, snow)          ! IN

    dsl4jsb_Get_var2D_onChunk(SOIL,     hcap_grnd)     ! OUT
    dsl4jsb_Get_var2D_onChunk(SOIL,     hcap_grnd_old) ! OUT
    dsl4jsb_Get_var2D_onChunk(SOIL,     vhcap)         ! OUT
    dsl4jsb_Get_var2D_onChunk(SOIL,     hdiff)         ! OUT
    dsl4jsb_Get_var3D_onChunk(SOIL,     t_soil_acoef)  ! OUT
    dsl4jsb_Get_var3D_onChunk(SOIL,     t_soil_bcoef)  ! OUT
    dsl4jsb_Get_var3D_onChunk(SOIL,     t_soil)        ! OUT
    dsl4jsb_Get_var2D_onChunk(SOIL,     grnd_hflx)     ! OUT
    dsl4jsb_Get_var2D_onChunk(SOIL,     grnd_hflx_old) ! OUT
    dsl4jsb_Get_var2D_onChunk(SOIL,     soil_depth)    ! OUT
    dsl4jsb_Get_var2D_onChunk(SOIL,     shc)           ! IN
    dsl4jsb_Get_var2D_onChunk(SOIL,     mpot)          ! IN
    dsl4jsb_Get_var2D_onChunk(SOIL,     bclapp)        ! IN
    dsl4jsb_Get_var2D_onChunk(SOIL,     vpor)          ! IN
    dsl4jsb_Get_var2D_onChunk(SOIL,     psdi)          ! IN

    dsl4jsb_Get_var3D_onChunk(HYDRO,    w_soil_sl)     ! OUT
    dsl4jsb_Get_var3D_onChunk(HYDRO,    w_soil_fc_sl)  ! OUT
    dsl4jsb_Get_var3D_onChunk(HYDRO,    w_soil_pwp_sl) ! OUT
    dsl4jsb_Get_var2D_onChunk(HYDRO,    w_skin_srf)    ! OUT
    dsl4jsb_Get_var2D_onChunk(HYDRO,    fract_water_srf) ! OUT
    dsl4jsb_Get_var2D_onChunk(HYDRO,    w_snow_srf)    ! OUT
    dsl4jsb_Get_var2D_onChunk(HYDRO,    fract_snow_srf) ! OUT
    dsl4jsb_Get_var2D_onChunk(HYDRO,    evapotrans)    ! OUT
    dsl4jsb_Get_var2D_onChunk(HYDRO,    trans)         ! OUT
    dsl4jsb_Get_var2D_onChunk(HYDRO,    evapopot)      ! OUT
    dsl4jsb_Get_var2D_onChunk(HYDRO,    w_snow_glacier) ! OUT
    dsl4jsb_Get_var2D_onChunk(HYDRO,    w_snow_can)    ! OUT
    dsl4jsb_Get_var2D_onChunk(HYDRO,    w_soil_root)   ! OUT
    dsl4jsb_Get_var2D_onChunk(HYDRO,    w_soil_root_fc) ! OUT
    dsl4jsb_Get_var2D_onChunk(HYDRO,    w_soil_column) ! OUT
    dsl4jsb_Get_var2D_onChunk(HYDRO,    w_skin_soil)   ! OUT
    dsl4jsb_Get_var2D_onChunk(HYDRO,    w_snow_soil)   ! OUT
    dsl4jsb_Get_var2D_onChunk(HYDRO,    runoff)        ! OUT
    dsl4jsb_Get_var2D_onChunk(HYDRO,    drainage)      ! OUT
    dsl4jsb_Get_var2D_onChunk(HYDRO,    pme_glacier)   ! OUT
    dsl4jsb_Get_var2D_onChunk(HYDRO,    w_soil_overflow) ! OUT

    dsl4jsb_Get_var2D_onChunk(SEB,      t_unfilt)      ! OUT
    dsl4jsb_Get_var2D_onChunk(HYDRO,    tte_corr)      ! OUT
    dsl4jsb_Get_var2D_onChunk(SOIL,     oro_stddev)    ! IN

    dsl4jsb_Get_var2D_onChunk(VEG, lai)                ! IN
    dsl4jsb_Get_var2D_onChunk(VEG, fract_fpc_max)      ! OUT
    dsl4jsb_Get_var2D_onChunk(VEG, root_depth)         ! OUT

    is_present(:)    = .TRUE.
    IF (tile%lct%id == GLACIER_TYPE) THEN
      is_glacier(:) = .TRUE.
    ELSE
      is_glacier(:) = .FALSE.
    END IF

    IF (ALL(w_snow_srf(:) < 0._wp)) THEN          ! Start of experiment
      lstart = .TRUE.
      w_snow_srf(:) = 0._wp
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

    CALL update_surf_down(nc, dtime, lstart, &
      & q_air(:), t_unfilt(:), wind_10m(:), t_air(:), w_skin_soil(:), fract_water_srf(:), &
      & wsr_max(:), w_snow_soil(:), fract_snow_srf(:), w_snow_can(:), w_snow_glacier(:), &
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

    CALL update_soiltemp(nc, dsl4jsb_Config(SOIL)%nsoil_water, dtime, lstart, &
      & t_unfilt(:), w_snow_soil(:), hdiff(:), vhcap(:), &
      & t_soil_acoef(:,:), t_soil_bcoef(:,:), t_soil(:,:), &
      & hcap_grnd(:), grnd_hflx(:), &
      & is_present(:), is_glacier(:))

!vg:  runoff_hd and drainage_hd no longer needed
!vg:  accumulation does not make sense (no reset)
    CALL update_surf_up(nc, dtime, dsl4jsb_Config(SOIL)%nsoil_water, &
      & w_skin_soil(:), fract_water_srf(:), wsr_max(:), &
      & w_soil_root(:), w_soil_root_fc(:), &
      & t_soil(:,1), fract_snow_srf(:), oro_stddev(:), &
      & evapotrans(:), evapopot(:), evapotrans_no_snow_skin(:), &
      & rain(:), is_present(:), is_glacier(:), &
      & runoff(:), drainage(:), & !vg runoff_hd(:), drainage_hd(:), &
      & melt_water_excess(:), skin_res_evap(:), &
      & dsl4jsb_Config(SOIL)%dz_water(1:dsl4jsb_Config(SOIL)%nsoil_water), &
      & w_soil_sl(:,:), w_soil_fc_sl(:,:), w_soil_pwp_sl(:,:), root_depth(:), soil_depth(:), &
      & shc(:), mpot(:), bclapp(:), vpor(:), &
      & trans(:), reduced_evap(:), &
!!$      & SPREAD(0.318517_wp, NCOPIES=nc, DIM=1), & !Temporary
!!$      & SPREAD(0.185245_wp, NCOPIES=nc, DIM=1), & !Temporary
      & psdi(:), &
      & w_soil_overflow(:), &
      & jllog, &
      & l_tpe_open, &
      & l_tpe_closed & ! IN
      & )

    w_soil_column(:) = SUM(w_soil_sl(:,:), DIM=2)

    w_snow_srf(:) = w_snow_soil(:) + w_snow_can(:)
    w_skin_srf(:) = w_skin_soil(:) ! + w_skin_can(:)

  END SUBROUTINE update_soil_1


  SUBROUTINE update_soil_2(tile, config, options)

    ! ---------------------------
    ! Variables

    ! Arguments
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    CLASS(t_jsb_config_p),      INTENT(in)    :: config(:)
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    ! Local variables


    ! Pointers to variables in memory
    REAL(wp), POINTER ::       &
      & fract          (:)

    ! Locally allocated vectors
    !
!!$    REAL(wp) :: &
    INTEGER :: iblk, ics, ice, nc

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':update_soil_2'

    ! ---------------------------
    ! Go

    iblk = options%iblk
    ics  = options%ics
    ice  = options%ice
    nc   = options%nc

    IF (.NOT. tile%Is_process_active(SOIL)) RETURN

    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting update')

    ! Get reference to variables for current block
    !
    fract                => tile%fract                               (ics:ice,    iblk)
    IF (ALL(fract(:) < EPSILON(1._wp))) RETURN

  END SUBROUTINE update_soil_2

    SUBROUTINE aggregate_soil(tile, config, options)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    CLASS(t_jsb_config_p),      INTENT(in)    :: config(:)
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

!!$    TYPE(t_jsb_model), POINTER :: model
!!$    CLASS(t_jsb_tile), POINTER :: tile

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':aggregate_soil'

    INTEGER :: iblk, ics, ice, nc

    iblk = options%iblk
    ics  = options%ics
    ice  = options%ice
    nc   = options%nc

  END SUBROUTINE aggregate_soil

#endif
END MODULE mo_soil_interface
