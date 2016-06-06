!>
!! @brief Contains the interfaces to the vegetation processes
!!
!! @author
!!  Reiner Schnur, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Reiner Schnur (2013-04-04)
!!
MODULE mo_veg_interface
#ifndef __NO_JSBACH__

#include "dsl_jsbach.inc"

  USE mo_jsb_control,        ONLY: debug_on
  USE mo_kind,               ONLY: wp
  USE mo_exception,          ONLY: message, finish
  USE mo_util,               ONLY: int2string

  USE mo_jsb_parallel_iface, ONLY: get_my_global_mpi_id

  USE mo_jsb_model_class,    ONLY: t_jsb_model
  USE mo_jsb_class,          ONLY: get_model
  USE mo_jsb_tile_class,     ONLY: t_jsb_tile_abstract
  USE mo_jsb_memory_class,   ONLY: t_jsb_memory, t_jsb_memory_p
  USE mo_jsb_config_class,   ONLY: t_jsb_config
  USE mo_jsb_process_class,  ONLY: t_jsb_process
  USE mo_jsb_task_class,     ONLY: t_jsb_process_task, t_jsb_task_options

  dsl4jsb_Use_processes VEG, SOIL, SEB, ATM2LAND, RAD, HYDRO !, ASSIMI
  dsl4jsb_Use_config(VEG)
  dsl4jsb_Use_config(SOIL)
  dsl4jsb_Use_config(RAD)

  dsl4jsb_Use_memory(atm2land)
  dsl4jsb_Use_memory(SOIL)
  dsl4jsb_Use_memory(SEB)
  dsl4jsb_Use_memory(VEG)
  dsl4jsb_Use_memory(RAD)
  dsl4jsb_Use_memory(HYDRO)
  !dsl4jsb_Use_memory(ASSIMI)

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: update_vegetation_2, update_vegetation_conductance
  PUBLIC :: Register_veg_tasks !,t_veg_process

!!$  TYPE, EXTENDS(t_jsb_process) :: t_veg_process
!!$  CONTAINS
!!$    PROCEDURE :: Register_tasks => Register_veg_tasks
!!$  END TYPE t_veg_process

  CHARACTER(len=*), PARAMETER :: modname = 'mo_veg_interface'

CONTAINS

  SUBROUTINE Register_veg_tasks(this, model_id)

    CLASS(t_jsb_process), INTENT(inout) :: this
    INTEGER,              INTENT(in)    :: model_id

    !CALL this%Register_task(tsk_radiation_surface(name="radiation_surface", process_id=VEG))

  END SUBROUTINE Register_veg_tasks

  SUBROUTINE update_vegetation_conductance(tile, options)
    ! ---------------------------
    ! Variables

    ! Used variables
    USE mo_veg_process, ONLY: get_canopy_conductance, get_water_stress_factor
    USE mo_veg_util,    ONLY: get_water_in_root_zone
    USE mo_phy_schemes, ONLY: qsat

    ! Arguments
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    dsl4jsb_Def_config(VEG)
    dsl4jsb_Def_config(SOIL)

    dsl4jsb_Def_memory(SOIL)
    dsl4jsb_Def_memory(SEB)
    dsl4jsb_Def_memory(VEG)
    dsl4jsb_Def_memory(ATM2LAND)
    dsl4jsb_Def_memory(HYDRO)
    !dsl4jsb_Def_memory(ASSIMI)

    dsl4jsb_Real2D_onChunk :: lai
    dsl4jsb_Real2D_onChunk :: swpar_srf_down
    dsl4jsb_Real2D_onChunk :: canopy_cond
    dsl4jsb_Real2D_onChunk :: canopy_cond_limited
    dsl4jsb_Real3D_onChunk :: soil_depth_l
    dsl4jsb_Real3D_onChunk :: root_depth_l
    !dsl4jsb_Real2D_onChunk :: w_soil_column   !< Water content in soil column
    dsl4jsb_Real3D_onChunk :: w_soil_sl            !< Water content in soil layers
    dsl4jsb_Real3D_onChunk :: w_soil_fc_sl         !< Water content of soil layers at field capacity
    dsl4jsb_Real3D_onChunk :: w_soil_pwp_sl        !< Water content of soil layers at permanent wilting point
    dsl4jsb_Real2D_onChunk :: w_soil_root         !< Water content in root zone
    dsl4jsb_Real2D_onChunk :: w_soil_root_fc      !< Water content in root zone at field capacity
    dsl4jsb_Real2D_onChunk :: w_soil_root_pwp     !< Water content in root zone at permanent wilting point
    dsl4jsb_Real2D_onChunk :: water_stress
    dsl4jsb_Real2D_onChunk :: t
    dsl4jsb_Real2D_onChunk :: q_air
    dsl4jsb_Real2D_onChunk :: press_srf

    INTEGER :: nsoil

    INTEGER :: p_pe
    INTEGER :: iblk, ics, ice, nc

    TYPE(t_jsb_model), POINTER :: model

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':update_vegetation_conductance'

    iblk = options%iblk
    ics  = options%ics
    ice  = options%ice
    nc   = options%nc

    IF (.NOT. tile%Is_process_active(VEG)) RETURN

!!$    p_pe = get_my_global_mpi_id()
!!$    IF (iblk == 1) CALL message(TRIM(routine), 'Starting update on PE '//int2string(p_pe), all_print=.TRUE.)
    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting update of vegetation_1 ')

    model => Get_model(tile%owner_model_id)

    dsl4jsb_Get_config(VEG)
    dsl4jsb_Get_config(SOIL)

    nsoil = dsl4jsb_Config(SOIL)%nsoil_water

    ! Get reference to variables for current block
    !
    dsl4jsb_Get_memory(VEG)
    dsl4jsb_Get_memory(SEB)
    dsl4jsb_Get_memory(SOIL)
    dsl4jsb_Get_memory(ATM2LAND)
    dsl4jsb_Get_memory(HYDRO)
    !dsl4jsb_Get_memory(ASSIMI)

    dsl4jsb_Get_var2D_onChunk(ATM2LAND, swpar_srf_down) ! IN
    dsl4jsb_Get_var2D_onChunk(ATM2LAND, q_air)          ! IN
    dsl4jsb_Get_var2D_onChunk(ATM2LAND, press_srf)      ! IN

    dsl4jsb_Get_var2D_onChunk(SEB,      t)              ! IN
    dsl4jsb_Get_var3D_onChunk(SOIL, soil_depth_l)       ! IN
    !dsl4jsb_Get_var2D_onChunk(HYDRO, w_soil_column)    ! OUT
    dsl4jsb_Get_var3D_onChunk(HYDRO, w_soil_sl)         ! IN
    dsl4jsb_Get_var3D_onChunk(HYDRO, w_soil_fc_sl)      ! IN
    dsl4jsb_Get_var3D_onChunk(HYDRO, w_soil_pwp_sl)     ! IN

    dsl4jsb_Get_var2D_onChunk(VEG,      lai)            ! IN
    dsl4jsb_Get_var3D_onChunk(VEG,      root_depth_l)   ! IN
    dsl4jsb_Get_var2D_onChunk(VEG,      canopy_cond)    ! OUT
    dsl4jsb_Get_var2D_onChunk(VEG,      canopy_cond_limited) ! OUT
    dsl4jsb_Get_var3D_onChunk(VEG,      root_depth_l)   ! OUT
    dsl4jsb_Get_var2D_onChunk(HYDRO,    w_soil_root)        ! OUT
    dsl4jsb_Get_var2D_onChunk(HYDRO,    w_soil_root_fc)     ! OUT
    dsl4jsb_Get_var2D_onChunk(HYDRO,    w_soil_root_pwp)    ! OUT
    dsl4jsb_Get_var2D_onChunk(HYDRO,    water_stress)   ! OUT



    ! Compute (max) canopy (stomatal) conductance with no water stress
    SELECT CASE(TRIM(dsl4jsb_Config(VEG)%scheme_stom_cond))
      CASE('echam5')  ! R: only if the assimilation process is not active!?
        canopy_cond(:) = get_canopy_conductance( lai(:), swpar_srf_down(:) )

!!$      CASE('bethy') ! R: the assimilation process is called from the main interface before this procedure and should not be called here!
!!$        ! TBD: Call bethy to compute unstressed canopy conductance

      CASE default
        CALL finish( TRIM(routine), 'Canopy conductance scheme '//TRIM(dsl4jsb_Config(VEG)%scheme_stom_cond)//&
                   & 'not supported.')
    END SELECT

    ! Get water stress

    ! Calculate water content and field capacity in the root zone
    w_soil_root   (:)  = get_water_in_root_zone(w_soil_sl     (:,:), soil_depth_l(:,:), root_depth_l(:,:))
    w_soil_root_fc(:)  = get_water_in_root_zone(w_soil_fc_sl  (:,:), soil_depth_l(:,:), root_depth_l(:,:))
    w_soil_root_pwp(:) = get_water_in_root_zone(w_soil_pwp_sl (:,:), soil_depth_l(:,:), root_depth_l(:,:))

    ! Calculate water stress factor
    ! Temporary: later use the distributed permanent wilting point from above and change get_water_stress_factor function
    water_stress(:) = get_water_stress_factor(w_soil_root(:), w_soil_root_fc(:), w_soil_root_pwp(:))

    ! Compute (actual) canopy (stomatal) conductance under water stress.
    ! R: diese canopy_cond (stressed) entspricht der canopy_conductance_limited die in JSBACH3 in update_soil gerechnet wurde.
    !    Sie wird auf jeden Fall benötigt unabhängig davon ob der Assimilations-Prozess an ist (wird dann als Input benötigt) oder nicht.
    canopy_cond_limited(:) = get_canopy_conductance(canopy_cond(:),            & ! in, unstressed canopy conductance
                                                    water_stress(:),           & ! in, water stress factor
                                                    q_air(:) >  qsat(t(:),press_srf(:)) & ! in, atmosphere saturated?
                                                    )

    ! Compute skin reservoir fraction on canopy

  END SUBROUTINE update_vegetation_conductance


  SUBROUTINE update_vegetation_2(tile, options)

    ! ---------------------------
    ! Variables

    ! Used variables
    ! sollen jetzt in dieses File hier kommen:
    !USE mo_veg_process, ONLY: canopy_albedo, get_foliage_projected_cover
    USE mo_veg_process, ONLY: get_foliage_projected_cover
    USE mo_jsb_time,    ONLY: get_time_interpolation_weights

    ! Arguments
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    !----------------------------
    ! Local variables
    REAL(wp), PARAMETER :: &
      albedo_snow_min = 0.4_wp,  & !< Minimum possible snow albedo   ! Temporary, will depend on vegetation type
      albedo_snow_max = 0.8_wp     !< Maximum possible snow albedo   ! Temporary


    INTEGER :: p_pe

    REAL(wp) :: wgt1, wgt2
    INTEGER  :: nmw1,nmw2

    ! Declare pointers of type process __conf and __mem
    dsl4jsb_Def_memory(VEG)
    dsl4jsb_Def_memory(SEB)
    dsl4jsb_Def_memory(SOIL)
    dsl4jsb_Def_memory(RAD)
    dsl4jsb_Def_memory(HYDRO)

    ! Declare pointers of type variable in memory
    dsl4jsb_Real2D_onChunk :: t
    dsl4jsb_Real2D_onChunk :: lai
    dsl4jsb_Real2D_onChunk :: fract_fpc
    dsl4jsb_Real2D_onChunk :: fract_fpc_max
    dsl4jsb_Real2D_onChunk :: fract_forest
!!$   dsl4jsb_Real2D_onChunk :: alb_vis_can
!!$   dsl4jsb_Real2D_onChunk :: alb_nir_can
!!$   dsl4jsb_Real2D_onChunk :: alb_vis_soil
!!$   dsl4jsb_Real2D_onChunk :: alb_nir_soil
    dsl4jsb_Real2D_onChunk :: alb_vis
    dsl4jsb_Real2D_onChunk :: alb_nir
    dsl4jsb_Real2D_onChunk :: alb_background
    dsl4jsb_Real2D_onChunk :: fract_snow_can
    dsl4jsb_Real2D_onChunk :: fract_snow_soil
    REAL(wp), POINTER :: lai_mon       (:,:)


    INTEGER  :: iblk, ics, ice, nc

    TYPE(t_jsb_model), POINTER :: model

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':update_vegetation_2'

    ! ---------------------------
    ! Go

    iblk = options%iblk
    ics  = options%ics
    ice  = options%ice
    nc   = options%nc

    IF (.NOT. tile%Is_process_active(VEG)) RETURN

    p_pe = get_my_global_mpi_id()
    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting update')

    model => Get_model(tile%owner_model_id)

    ! Get reference to variables for current block
    !
    dsl4jsb_Get_memory(VEG)
    dsl4jsb_Get_memory(SEB)
    dsl4jsb_Get_memory(SOIL)

    dsl4jsb_Get_var2D_onChunk(VEG,   lai)            ! OUT
    dsl4jsb_Get_var2D_onChunk(VEG,   fract_fpc)      ! OUT
    dsl4jsb_Get_var2D_onChunk(VEG,   fract_fpc_max)  ! OUT
    dsl4jsb_Get_var2D_onChunk(VEG,   fract_forest)   ! OUT
    dsl4jsb_Get_var2D_onChunk(HYDRO, fract_snow_can) ! OUT
!!$    alb_vis_can  => tile%mem%veg  %alb_vis_can%ptr     (ics:ice,    iblk)
!!$    alb_nir_can  => tile%mem%veg  %alb_nir_can%ptr     (ics:ice,    iblk)

!!$    alb_vis_soil => tile%mem%soil %alb_vis_soil%ptr    (ics:ice,    iblk)
!!$    alb_nir_soil => tile%mem%soil %alb_nir_soil%ptr    (ics:ice,    iblk)
    dsl4jsb_Get_var2D_onChunk(SEB,   t)              ! OUT
    dsl4jsb_Get_var2D_onChunk(RAD,   alb_vis)        ! OUT
    dsl4jsb_Get_var2D_onChunk(RAD,   alb_nir)        ! OUT
    dsl4jsb_Get_var2D_onChunk(RAD,   alb_background) ! OUT

    dsl4jsb_Get_var2D_onChunk(HYDRO, fract_snow_soil)! OUT

    lai_mon      => veg__mem%lai_mon (ics:ice,    iblk, :)

    ! Update leaf area index
    !
    ! At the moment by interpolating monthly climatological LAI
    CALL get_time_interpolation_weights(wgt1, wgt2, nmw1, nmw2)
    lai(:) = REAL(wgt1,wp) * lai_mon(:,nmw1+1) + REAL(wgt2,wp) * lai_mon(:,nmw2+1)

    ! Update foliage projected cover
    !
    ! TBD: Use clumpiness factor instead of 1._wp
    fract_fpc(:) = get_foliage_projected_cover(fract_fpc_max(:), lai(:), SPREAD(1._wp, DIM=1, NCOPIES=SIZE(lai,1)))

    ! Update albedo
    ! TBD: use newer albedo scheme
! R: entspricht meiner Funktion get_surface_albedo_simple. Wird bei mir aber vom Hauptinterface gecallt und soll auch nicht in update_vegetation. 
!     alb_vis(:) = canopy_albedo(t(:), fract_forest(:), fract_snow_soil(:), fract_snow_can(:), lai(:), &
!                    & SPREAD(albedo_snow_min, NCOPIES=nc, DIM=1), SPREAD(albedo_snow_max, NCOPIES=nc, DIM=1), alb_background(:))
!     alb_nir(:) = alb_vis(:)

  END SUBROUTINE update_vegetation_2


    SUBROUTINE aggregate_veg(tile, options)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    TYPE(t_jsb_model), POINTER :: model

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':aggregate_veg'

    INTEGER :: iblk, ics, ice, nc

    iblk = options%iblk
    ics  = options%ics
    ice  = options%ice
    nc   = options%nc

    model => Get_model(tile%owner_model_id)

  END SUBROUTINE aggregate_veg

#endif
END MODULE mo_veg_interface
