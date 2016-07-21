!>
!! @brief Contains the interfaces to the radiation processes
!!
!! @author
!!  Rainer Schneck, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Rainer Schneck (2015-04-14)
!!

MODULE mo_rad_interface
#ifndef __NO_JSBACH__

#include "dsl_jsbach.inc"
  ! -------------------------------------------------------------------------------------------------------
  ! Used variables of module

  ! Use of basic structures
  USE mo_jsb_control,         ONLY: debug_on
  USE mo_kind,                ONLY: wp
  USE mo_exception,           ONLY: message, finish

  USE mo_jsb_model_class,     ONLY: t_jsb_model
  USE mo_jsb_class,           ONLY: Get_model
  USE mo_jsb_tile_class,      ONLY: t_jsb_tile_abstract, t_jsb_aggregator
  USE mo_jsb_lct_class,       ONLY: LAND_TYPE, LAKE_TYPE, GLACIER_TYPE, BARE_TYPE, VEG_TYPE, Contains_lct
  USE mo_jsb_config_class,    ONLY: t_jsb_config, t_jsb_config_p
  USE mo_jsb_process_class,   ONLY: t_jsb_process
  USE mo_jsb_task_class,      ONLY: t_jsb_process_task, t_jsb_task_options
  !USE mo_jsb_aggregate_class, ONLY: t_jsb_aggregator

  ! Use of processes in this module
  dsl4jsb_Use_processes SEB, SOIL, VEG, PHENO, ATM2LAND, RAD, HYDRO, ASSIMI

  ! Use process configurations
  dsl4jsb_Use_config(SEB)
  dsl4jsb_Use_config(SOIL)
  dsl4jsb_Use_config(RAD)
  dsl4jsb_Use_config(ASSIMI)

  ! Use of process memories
  dsl4jsb_Use_memory(ATM2LAND)
  dsl4jsb_Use_memory(SOIL)
  dsl4jsb_Use_memory(SEB)
  dsl4jsb_Use_memory(VEG)
  dsl4jsb_Use_memory(PHENO)
  dsl4jsb_Use_memory(RAD)
  dsl4jsb_Use_memory(HYDRO)

  ! -------------------------------------------------------------------------------------------------------
  ! Module variables

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: Register_rad_tasks

  !> Type definition for radiation_surface task
  TYPE, EXTENDS(t_jsb_process_task) :: tsk_surface_radiation
  CONTAINS
    PROCEDURE, NOPASS :: Integrate => update_surface_radiation
    PROCEDURE, NOPASS :: Aggregate => aggregate_surface_radiation
  END TYPE tsk_surface_radiation

  INTERFACE tsk_surface_radiation
    PROCEDURE Create_task_surface_radiation
  END INTERFACE tsk_surface_radiation

  !> Type definition for radiation_par task
  TYPE, EXTENDS(t_jsb_process_task) :: tsk_radiation_par
  CONTAINS
    PROCEDURE, NOPASS :: Integrate => update_radiation_par
    ! R: Ich sehe keinen Sinn darin eine der Outputvariablen zu aggregieren, muss aber eine Prozedur bereitstellen
    PROCEDURE, NOPASS :: Aggregate => aggregate_radiation_par
  END TYPE tsk_radiation_par

  INTERFACE tsk_radiation_par
    PROCEDURE Create_task_radiation_par
  END INTERFACE tsk_radiation_par


  !> Type definition for albedo task
  TYPE, EXTENDS(t_jsb_process_task) :: tsk_albedo
  CONTAINS
    PROCEDURE, NOPASS :: Integrate => update_albedo
    PROCEDURE, NOPASS :: Aggregate => aggregate_albedo
  END TYPE tsk_albedo

  INTERFACE tsk_albedo
    PROCEDURE Create_task_albedo
  END INTERFACE tsk_albedo

  CHARACTER(len=*), PARAMETER :: modname = 'mo_rad_interface'

CONTAINS

  ! ================================================================================================================================
  !! Constructors for tasks

  ! -------------------------------------------------------------------------------------------------------
  !> Constructor for surface_radiation task
  !!
  !! @param[in]     model_id     Model id
  !! @return        return_ptr   Instance of process task "radiation"
  !!
  FUNCTION Create_task_surface_radiation(model_id) RESULT(return_ptr)

    INTEGER,                   INTENT(in) :: model_id
    CLASS(t_jsb_process_task), POINTER    :: return_ptr

    ALLOCATE(tsk_surface_radiation::return_ptr)

    return_ptr%name = 'surface_radiation'
    return_ptr%process_id = RAD
    return_ptr%owner_model_id = model_id

  END FUNCTION Create_task_surface_radiation

  ! -------------------------------------------------------------------------------------------------------
  !> Constructor for radiation_par task
  !!
  !! @param[in]     model_id     Model id
  !! @return        return_ptr   Instance of process task "radiation"
  !!
  FUNCTION Create_task_radiation_par(model_id) RESULT(return_ptr)

    INTEGER,                   INTENT(in) :: model_id
    CLASS(t_jsb_process_task), POINTER    :: return_ptr

    ALLOCATE(tsk_radiation_par::return_ptr)

    return_ptr%name = 'radiation_par'
    return_ptr%process_id = RAD
    return_ptr%owner_model_id = model_id

  END FUNCTION Create_task_radiation_par
  
  ! -------------------------------------------------------------------------------------------------------
  !> Constructor for albedo task
  !!
  !! @param[in]     model_id     Model id
  !! @return        return_ptr   Instance of process task "albedo"
  !!
  FUNCTION Create_task_albedo(model_id) RESULT(return_ptr)

    INTEGER,                   INTENT(in) :: model_id
    CLASS(t_jsb_process_task), POINTER    :: return_ptr

    ALLOCATE(tsk_albedo::return_ptr)

    return_ptr%name = 'albedo'
    return_ptr%process_id = RAD
    return_ptr%owner_model_id = model_id

  END FUNCTION Create_task_albedo

  ! ================================================================================================================================
  !> Register tasks for radiation process
  !!
  !! @param[in,out] this      Instance of <PROCESS_NAME_LOWER_CASE> process class
  !! @param[in]     model_id  Model id
  !!
  SUBROUTINE Register_rad_tasks(this, model_id)

    CLASS(t_jsb_process), INTENT(inout) :: this
    INTEGER, INTENT(in) :: model_id

    CALL this%Register_task(tsk_surface_radiation(model_id))
    CALL this%Register_task(tsk_radiation_par    (model_id))
    CALL this%Register_task(tsk_albedo           (model_id))

  END SUBROUTINE Register_rad_tasks

  ! ================================================================================================================================
  !>
  !! @brief Implementation to update the surface radiation
  !! Task "surface_radiation" calculates the net radiations at the surface.
  !!
  !! @param[in,out] tile    Tile for which routine is executed.
  !! @param[in]     options Additional run-time parameters.
  !!
  SUBROUTINE update_surface_radiation(tile, options)

    USE mo_rad_process, ONLY: calc_radiation_surface_net

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options), INTENT(in) :: options

    ! Declare pointers to process configs and memory
    dsl4jsb_Def_config(RAD)
    dsl4jsb_Def_memory(ATM2LAND)
    dsl4jsb_Def_memory(SEB)
    dsl4jsb_Def_memory(RAD)
    dsl4jsb_Def_memory(VEG)

    ! Declare pointers for variables in memory
    dsl4jsb_Real2D_onChunk :: &
      & swvis_srf_net,        &
      & swnir_srf_net,        &
      & sw_srf_net,           &
      & lw_srf_net,           &
      & rad_srf_net,          &
      & lw_srf_down,          &
      & swvis_srf_down,       &
      & swnir_srf_down,       &
      & alb_vis,              &
      & alb_nir,              &
      & fract_lice,           &
      & rad_net_lwtr,         &
      & rad_net_lice,         &
      & albedo_lwtr,          &
      & albedo_lice,          &
      & t,                    &
      & t_lwtr,               &
      & t_lice

    REAL(wp) ::                  &
      & lw_net_lwtr(options%nc), &
      & lw_net_lice(options%nc), &
      & sw_net_lwtr(options%nc), &
      & sw_net_lice(options%nc)

    INTEGER :: iblk, ics, ice, nc

    TYPE(t_jsb_model), POINTER :: model

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':update_surface_radiation'
    
    ! Get local variables from options argument
    iblk = options%iblk
    ics  = options%ics
    ice  = options%ice
    nc   = options%nc
    
    ! If process is not active on this tile, do nothing
    IF (.NOT. tile%Is_process_active(RAD)) RETURN
    
    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    model => Get_model(tile%owner_model_id)

    dsl4jsb_Get_memory(ATM2LAND)
    dsl4jsb_Get_memory(RAD)
    dsl4jsb_Get_memory(SEB)

    dsl4jsb_Get_var2D_onChunk(ATM2LAND, swvis_srf_down) ! in
    dsl4jsb_Get_var2D_onChunk(ATM2LAND, swnir_srf_down) ! in
    dsl4jsb_Get_var2D_onChunk(ATM2LAND, lw_srf_down)    ! in

    dsl4jsb_Get_var2D_onChunk(RAD,      sw_srf_net)     ! out
    dsl4jsb_Get_var2D_onChunk(RAD,      lw_srf_net)     ! out
    dsl4jsb_Get_var2D_onChunk(RAD,      rad_srf_net)    ! out

    ! ---------------------------
    ! Go

    IF (tile%lcts(1)%id == LAKE_TYPE) THEN

      ! Water
      dsl4jsb_Get_var2D_onChunk(SEB,      t_lwtr)       ! in
      dsl4jsb_Get_var2D_onChunk(ATM2LAND, albedo_lwtr)  ! in
      dsl4jsb_Get_var2D_onChunk(RAD,      rad_net_lwtr) ! out

      CALL calc_radiation_surface_net( &
        & swvis_srf_down(:),        & ! in
        & swnir_srf_down(:),        & ! in
        & albedo_lwtr(:),           & ! in
        & albedo_lwtr(:),           & ! in
        & lw_srf_down(:),           & ! in
        & t_lwtr(:),                & ! in
        & rad_net_lwtr(:),          & ! out
        !& swvis_srf_net(:),        & ! out
        !& swnir_srf_net(:),        & ! out
        & sw_net=sw_net_lwtr(:),    & ! out
        & lw_net=lw_net_lwtr(:)     & ! out
        & )

      ! Ice on lake
      dsl4jsb_Get_var2D_onChunk(SEB, t_lice)       ! in
      dsl4jsb_Get_var2D_onChunk(SEB, fract_lice)   ! in
      dsl4jsb_Get_var2D_onChunk(RAD, albedo_lice)  ! in
      dsl4jsb_Get_var2D_onChunk(RAD, rad_net_lice) ! out

      CALL calc_radiation_surface_net( &
        & swvis_srf_down(:),        & ! in
        & swnir_srf_down(:),        & ! in
        & albedo_lice(:),           & ! in
        & albedo_lice(:),           & ! in
        & lw_srf_down(:),           & ! in
        & t_lice(:),                & ! in
        & rad_net_lice(:),          & ! out
        !& swvis_srf_net(:),        & ! out
        !& swnir_srf_net(:),        & ! out
        & sw_net=sw_net_lice(:),    & ! out
        & lw_net=lw_net_lice(:)     & ! out
        & )

      rad_srf_net(:) = (1._wp - fract_lice(:)) * rad_net_lwtr(:) + fract_lice(:) * rad_net_lice(:)
      sw_srf_net (:) = (1._wp - fract_lice(:)) * sw_net_lwtr(:)  + fract_lice(:) * sw_net_lice(:)
      lw_srf_net (:) = (1._wp - fract_lice(:)) * lw_net_lwtr(:)  + fract_lice(:) * lw_net_lice(:)

    ELSE

      dsl4jsb_Get_var2D_onChunk(SEB,      t)              ! in
      dsl4jsb_Get_var2D_onChunk(RAD,      alb_vis)        ! in
      dsl4jsb_Get_var2D_onChunk(RAD,      alb_nir)        ! in
      dsl4jsb_Get_var2D_onChunk(RAD,      swvis_srf_net)  ! out
      dsl4jsb_Get_var2D_onChunk(RAD,      swnir_srf_net)  ! out

      CALL calc_radiation_surface_net( &
        & swvis_srf_down(:),    & ! in
        & swnir_srf_down(:),    & ! in
        & alb_vis(:),           & ! in
        & alb_nir(:),           & ! in
        & lw_srf_down(:),       & ! in
        & t(:),                 & ! in
        & rad_srf_net(:),       & ! out
        & swvis_srf_net(:),     & ! out
        & swnir_srf_net(:),     & ! out
        & sw_srf_net(:),        & ! out
        & lw_srf_net(:)         & ! out
        & )
    END IF

    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE update_surface_radiation

  ! -------------------------------------------------------------------------------------------------------
  !>
  !! Implementation of "aggregate" for task "surface_radiation"
  !!
  !! @param[in,out] tile    Tile for which routine is executed.
  !! @param[in]     options Additional run-time parameters.
  !!
  SUBROUTINE aggregate_surface_radiation(tile, options)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    dsl4jsb_Def_memory(RAD)

    CLASS(t_jsb_aggregator), POINTER :: weighted_by_fract

    TYPE(t_jsb_model), POINTER :: model

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':aggregate_surface_radiation'

    INTEGER :: iblk, ics, ice

    iblk = options%iblk
    ics  = options%ics
    ice  = options%ice

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    model => Get_model(tile%owner_model_id)

    dsl4jsb_Get_memory(RAD)

    weighted_by_fract => tile%Get_aggregator("weighted_by_fract")

    dsl4jsb_Aggregate_onChunk(RAD, rad_srf_net,    weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(RAD, lw_srf_net,     weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(RAD, sw_srf_net,     weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(RAD, swvis_srf_net,  weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(RAD, swnir_srf_net,  weighted_by_fract)

    dsl4jsb_Aggregate_onChunk(RAD, rad_net_lwtr,   weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(RAD, rad_net_lice,   weighted_by_fract)

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE aggregate_surface_radiation
  
  ! ================================================================================================================================
  !>
  !! @brief Implementation to update the photosynthetically active radiation
  !! Task "radiation_par" calculates PAR and APAR.
  !!
  !! @param[in,out] tile    Tile for which routine is executed.
  !! @param[in]     options Additional run-time parameters.
  !!
  SUBROUTINE update_radiation_par(tile, options)

    ! Use declarations
    USE mo_rad_process,    ONLY: calc_fapar_lai_cl, calc_par

    ! Arguments
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    ! Local variables
    INTEGER                           :: icanopy

    INTEGER                           :: iblk, ics, ice, nc
    REAL(wp)                          :: dtime

    REAL(wp)                          :: canopy_bound_lai_delta

    REAL(wp), DIMENSION(options%nc)   :: B4_layer_above ! R: cannot be scalar as it will be used as inout argument for elemental procedure
                                                        !    and thus has to have the same dimensions as the other out arguments.

    TYPE(t_jsb_model), POINTER :: model

    ! Declare pointers for process configuration and memory
    dsl4jsb_Def_config(RAD)
    dsl4jsb_Def_config(ASSIMI)
    dsl4jsb_Def_config(SEB)
    dsl4jsb_Def_memory(RAD)
    dsl4jsb_Def_memory(SEB)
    dsl4jsb_Def_memory(PHENO)
    dsl4jsb_Def_memory(ATM2LAND)

    INTEGER :: ncanopy               ! number of canopy layers
    !INTEGER, POINTER       :: acc_counter           ! counter for each delta_time step since the last variable output
    LOGICAL :: use_alb_veg_simple    !

    dsl4jsb_Real2D_onChunk :: par                   ! Downward PAR flux in [W/m^2], Photosynth. Active Radiation (400-700nm)
    dsl4jsb_Real2D_onChunk :: par_down_mol          ! Downward PAR flux in mol (photons)/(m^2 s)
    dsl4jsb_Real2D_onChunk :: par_down_mol_nacc     ! PAR accumulated [MOL PHOTONS/M^2 ]
    dsl4jsb_Real2D_onChunk :: par_down_mol_tavg     ! PAR [MOL PHOTONS/M^2 S] (time mean value)
    dsl4jsb_Real2D_onChunk :: soil_reflectivity_par ! Soil reflectivity within the PAR spectrum [-]
    dsl4jsb_Real2D_onChunk :: alb_vis_soil          ! Soil albedo in the visible range
    !dsl4jsb_Real2D_onChunk :: alb_background        ! USED ONLY FOR COMPARISON WITH JSBACH3
    dsl4jsb_Real2D_onChunk :: lai                   ! Leaf area index [-]
    dsl4jsb_Real2D_onChunk :: cos_zenith_angle      ! Angle of incoming radiation [-]
    dsl4jsb_Real2D_onChunk :: frac_par_direct       ! Fraction of direct (in contrast to diffuse) PAR [-]
    dsl4jsb_Real2D_onChunk :: frac_par_diffuse      ! Fraction of diffuse (in contrast to direct) PAR [-]
    dsl4jsb_Real2D_onChunk :: canopy_bound_lai
    dsl4jsb_Real2D_onChunk :: apar_nacc

    dsl4jsb_Real3D_onChunk :: faPAR_cl              ! Fraction of absorbed PAR per leaf area for each canopy layer
    dsl4jsb_Real3D_onChunk :: apar_per_lai_cl       ! Absorbed PAR of canopy layer
    dsl4jsb_Real3D_onChunk :: lai_cl                ! Leaf area index of canopy layer
    dsl4jsb_Real3D_onChunk :: apar_per_lai_cl_nacc  ! Absorbed PAR accumulated [MOL PHOTONS/M^2 ]
    dsl4jsb_Real3D_onChunk :: apar_per_lai_cl_tavg  ! Absorbed PAR [MOL PHOTONS/M^2 S] (time mean value)

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':update_radiation_par'
    
    ! Get local variables from options argument
    iblk  = options%iblk
    ics   = options%ics
    ice   = options%ice
    nc    = options%nc
    dtime = options%dtime ! time step in seconds e.g. for one day=86400

    IF (tile%lcts(1)%id /= VEG_TYPE) RETURN

    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')
    
    ! Get model
    model => Get_model(tile%owner_model_id)
    
    ! Get pointers to process configs and memory
    dsl4jsb_Get_config(RAD)
    dsl4jsb_Get_config(ASSIMI)
    dsl4jsb_Get_memory(RAD)
    dsl4jsb_Get_memory(PHENO)
    dsl4jsb_Get_memory(ATM2LAND)

    ! If process is not active on this tile, do nothing
    IF (.NOT. tile%Is_process_active(RAD)) RETURN

    ncanopy                     =  dsl4jsb_Config(ASSIMI)%ncanopy
    use_alb_veg_simple          =  dsl4jsb_Config(RAD)%use_alb_veg_simple
    canopy_bound_lai(0:ncanopy) => dsl4jsb_Config(ASSIMI)% canopy_bound_lai(0:ncanopy)
    
    ! Set pointers to variables in memory
    dsl4jsb_Get_var2D_onChunk(ATM2LAND, cos_zenith_angle)      ! in (received from the atmosphere)
    dsl4jsb_Get_var2D_onChunk(ATM2LAND, frac_par_diffuse)      ! in (received from the atmosphere)
    dsl4jsb_Get_var2D_onChunk(PHENO,    lai)                   ! in    
    dsl4jsb_Get_var2D_onChunk(RAD,      frac_par_direct)       ! in R: should be received from the atmosphere but it is NOT in the moment!
    dsl4jsb_Get_var2D_onChunk(RAD,      alb_vis_soil)          ! in (initialized from input file)
    !dsl4jsb_Get_var2D_onChunk(RAD,      alb_background)        ! in (initialized from input file)
    dsl4jsb_Get_var2D_onChunk(RAD,      par)                   ! inout
    dsl4jsb_Get_var2D_onChunk(RAD,      par_down_mol)          ! inout

    ! R: These variables could also be local as they are not used outside this procedure in the moment.
    !    However I want it in the std model output and thus it has to be in the memory
    dsl4jsb_Get_var2D_onChunk(RAD,      par_down_mol_nacc)     ! inout
    dsl4jsb_Get_var2D_onChunk(RAD,      par_down_mol_tavg)     ! inout
    dsl4jsb_Get_var2D_onChunk(RAD,      soil_reflectivity_par) ! inout
    dsl4jsb_Get_var3D_onChunk(RAD,      lai_cl)                ! inout
    dsl4jsb_Get_var3D_onChunk(RAD,      faPAR_cl)              ! inout
    ! used outside for assimilation:
    dsl4jsb_Get_var3D_onChunk(RAD,      apar_per_lai_cl)       ! inout
    dsl4jsb_Get_var3D_onChunk(RAD,      apar_per_lai_cl_nacc)  ! inout
    dsl4jsb_Get_var3D_onChunk(RAD,      apar_per_lai_cl_tavg)  ! inout
    dsl4jsb_Get_var2D_onChunk(RAD,      apar_nacc)             ! inout
    
    ! ---------------------------
    ! Go
    
    ! Initialize local variables
    par_down_mol           = 0._wp ! for now also in the memory but is only for current time step/"delta_time"
    soil_reflectivity_par  = 0._wp ! for now also in the memory but is only for current time step/"delta_time"
    apar_per_lai_cl        = 0._wp

    ! Compute direct fraction of PAR.
    ! R:  frac_par_direct = in ECHAM jsswdifpar, which is zsw_par_frac_diffuse = frac_diffuse_par  =sw_par_fract_direct
    frac_par_direct = 1._wp - frac_par_diffuse
                                                 
    B4_layer_above = 7777777 ! Used to transfer variable B4 from one canopy layer to the next. Startvalue only to indicate first layer
    DO icanopy=1,ncanopy
       ! Compute faPAR_cl (absorbed photsynthetic active radiation for each canopy layer) and lai_cl (leaf area index for each canopy layer)
       canopy_bound_lai_delta = canopy_bound_lai(icanopy) - canopy_bound_lai(icanopy-1)

       CALL calc_par(                 &
         & dtime,                     & ! in
         & ncanopy,                   & ! in
         & icanopy,                   & ! in
         & use_alb_veg_simple,        & ! in
         & par(:),                    & ! in
         & alb_vis_soil(:),           & ! in
         & lai(:),                    & ! in
         & cos_zenith_angle(:),       & ! in
         & frac_par_direct(:),        & ! in
         & canopy_bound_lai(icanopy), & ! in
         & canopy_bound_lai_delta,    & ! in
         & B4_layer_above,            & ! inout
         & par_down_mol(:),           & ! out, used outside for assimilation
         & soil_reflectivity_par(:),  & ! out, R: could also be local as it is not used outside in the moment,
                                        ! however I want it in the model output and thus it has to be updated by returning it
         & faPAR_cl(:,icanopy),       & ! out, R: could also be local as it is not used outside in the moment,
                                        ! however I want it in the model output and thus it has to be updated by returning it
         & lai_cl(:,icanopy),         & ! out, R: could also be local as it is not used outside in the moment,
                                        ! however I want it in the model output and thus it has to be updated by returning it
         & apar_per_lai_cl(:,icanopy) ) ! out, used outside for assimilation
    END DO

    ! R: Note, the following code is preliminary and should be replaced later by an laccu-procedure.
    !    For the correct integration over time as it is done now, this task must be called only once per timestep.
    ! Accumulate incoming par in mol(photons)/(m^2) over time
    par_down_mol_nacc(:) = par_down_mol_nacc(:) + par_down_mol(:) * dtime ! Change explanation:
                                                            ! "par/Epar" was unnecessarily double calculated. In JSBACH3
                                                            ! "CALL add(..., laccu=.TRUE.)" accumulates variable and sets it to zero
                                                            ! after writen out. As "laccu" does not exist in JSBACH4 (it could be
                                                            ! added) I added a workaround.
                                                            ! acc_counter, par_down_mol_nacc, apar_per_lai_cl_nacc is set to zero
                                                            ! when the output is written.
    ! R: This whole workaround should get abitary soon...
    !acc_counter = acc_counter + 1
    !par_down_mol_tavg = par_down_mol_nacc / dtime / acc_counter
    apar_per_lai_cl_nacc(:,:)  = apar_per_lai_cl_nacc(:,:) + apar_per_lai_cl(:,:) * dtime
    ! R: This whole workaround (acc_counter) should get abitary soon...
    !apar_per_lai_cl_tavg(:,icanopy)  = apar_per_lai_cl_nacc(:,icanopy) / dtime / acc_counter
    DO icanopy=1,ncanopy
      apar_nacc(:) = apar_nacc(:) + par_down_mol(:) * faPAR_cl(:,icanopy) * dtime
    END DO

!        ! Compute absorbed PAR per leaf area in canopy layer [units: (absorbed photons) / (m^2(leaf area) s)] from
!        ! par and fraction of absorbed PAR
!        apar_per_lai_cl(:,icanopy) = par_down_mol(:) * faPAR_cl(:,icanopy) / (MAX(lai_cl(:,icanopy),1.e-10_wp))
!        ! R: Folgendes fällt weg:
!        ! mask(1:nidx,itile) <--theLand%Surface%is_vegetation(kidx0:kidx1,1:ntiles)
!        ! apar_per_lai_cl(:,icanopy) = MERGE(apar_per_lai_cl(:,icanopy),0._dp,mask(1:nidx,itile))
!        ! Da damit die Wüstenflächen auf 0 gesetzt werden und JSBACH4 diese Prozedur dann
!        ! gar nicht erst aufruft für diese tile.
!        IF (waterLimitationFlag) THEN ! to make sure to accumulate only once per timestep
!           apar_per_lai_cl_nacc(:,icanopy)  = apar_per_lai_cl_nacc(:,icanopy) + apar_per_lai_cl(:,icanopy) * dtime
!           ! R: This whole workaround (acc_counter) should get abitary soon...
!           !apar_per_lai_cl_tavg(:,icanopy)  = apar_per_lai_cl_nacc(:,icanopy) / dtime / acc_counter
!           apar_nacc(:)                     = apar_nacc(:) + par_down_mol(:) * faPAR_cl(:,icanopy) * dtime
!        END IF

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE update_radiation_par

  ! -------------------------------------------------------------------------------------------------------
  !>
  !! Implementation of "aggregate" for task "radiation_par"
  !!
  !! @param[in,out] tile    Tile for which routine is executed.
  !! @param[in]     options Additional run-time parameters.
  !!
  SUBROUTINE aggregate_radiation_par(tile, options)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    TYPE(t_jsb_model), POINTER :: model

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':aggregate_radiation_par'

    INTEGER :: iblk, ics, ice, nc

    iblk = options%iblk
    ics  = options%ics
    ice  = options%ice
    nc   = options%nc

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    model => Get_model(tile%owner_model_id)

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE aggregate_radiation_par
  
  ! ================================================================================================================================
  !>
  !! @brief Implementation to update the albedo
  !! Task "albedo" calculates several albedos.
  !!
  !! @param[in,out] tile    Tile for which routine is executed.
  !! @param[in]     options Additional run-time parameters.
  !!
  SUBROUTINE update_albedo(tile, options)

    ! Use declarations
    USE mo_rad_process,    ONLY: get_surface_albedo_simple, calc_soil_albedo, calc_alb_lice, calc_glacier_albedo, &
                                 calc_sky_view_fractions, calc_snow_albedo, Merge_albedos_of_vegtile, &
                                 Has_minimal_radiation

    ! Arguments
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    ! Local variables
    TYPE(t_jsb_model), POINTER   :: model
    INTEGER                      :: iblk, ics, ice, nc
    REAL(wp)                     :: dtime

    LOGICAL                      :: &
                                    & l_day(options%nc), &
                                    &  ForestFlag  ! to be taken from the lctlib ("true" for forest landcover types)

    REAL(wp)                     :: StemArea             ! to be taken from the lctlib. (Area of stems and branches of woody plants.)

    CHARACTER(len=*), PARAMETER    :: routine = TRIM(modname)//':update_albedo'

    ! Declare pointers to process configuration and memory
    dsl4jsb_Def_config(RAD)
    dsl4jsb_Def_memory(RAD)
    dsl4jsb_Def_memory(SEB)
    dsl4jsb_Def_memory(ATM2LAND)
    !RMS dsl4jsb_Def_memory(CARBON)
    dsl4jsb_Def_memory(HYDRO)
    dsl4jsb_Def_memory(VEG)
    dsl4jsb_Def_memory(PHENO)

    ! Declare configuration variables
    CHARACTER(len=10) :: use_alb_soil_organic_C

    REAL(wp) :: &
      & albedo_age_weight, &    ! 0 < albedo_age_weight < 1: snow albedo is calcualted by linearly weighting
                                !                            the snow albedo resulting from both schemes
      & AlbedoCanopySnow
                                              
    LOGICAL ::           &
      & use_alb_veg_simple,        &
      & use_alb_soil_scheme,       &
      & use_alb_soil_litter,       &
      & use_alb_mineralsoil_const, &
      & use_alb_canopy
      
    ! Declare pointers to variables in memory
    dsl4jsb_Real2D_onChunk ::  &
      & cos_zenith_angle,      &
      & swvis_srf_down,        &
      & swnir_srf_down,        &
      & t,                     &
      & alb_vis,               &
      & alb_nir,               &
      & alb_vis_can,           &
      & alb_nir_can,           &
      & alb_vis_soil,          &
      & alb_nir_soil,          &
      & alb_vis_snow,          &
      & alb_nir_snow,          &
      & alb_background,        &
      & alb_vis_mineralsoil,   &
      & alb_nir_mineralsoil,   &
      & sky_view_fract,        &
      & sky_view_fract_stem,   &
      & snow_age,              &
      & lai,                   &
      & fract_forest,          &
      & fract_snow_can,        &
      & fract_snow_soil,       &
      & Cpool_litter_green_ag, &
      & Cpool_slow,            &
      & swvis_srf_net,         &
      & swnir_srf_net
    ! Lakes
    dsl4jsb_Real2D_onChunk ::  &
      & albedo_lwtr,           &
      & fract_lice,            &
      & t_lice,                &
      & w_snow_lice,           &
      & albedo_lice

    ! Set local variables from options argument
    iblk  = options%iblk
    ics   = options%ics
    ice   = options%ice
    nc    = options%nc
    dtime = options%dtime ! time step in seconds e.g. for one day=86400

    ! Set model
    model => Get_model(tile%owner_model_id)

    ! Set pointers to process configs and memory
    dsl4jsb_Get_config(RAD)
    dsl4jsb_Get_memory(RAD)
    dsl4jsb_Get_memory(ATM2LAND)
    dsl4jsb_Get_memory(SEB)
    !RMS dsl4jsb_Get_memory(CARBON)
    dsl4jsb_Get_memory(HYDRO)
    IF (Contains_lct(tile%lcts, VEG_TYPE)) THEN
      dsl4jsb_Get_memory(VEG)
      dsl4jsb_Get_memory(PHENO)
    END IF
    
    ! Set pointers to variables in config
    use_alb_veg_simple        = dsl4jsb_Config(RAD) %use_alb_veg_simple
    use_alb_soil_scheme       = dsl4jsb_Config(RAD) %use_alb_soil_scheme
    use_alb_mineralsoil_const = dsl4jsb_Config(RAD) %use_alb_mineralsoil_const
    use_alb_soil_litter       = dsl4jsb_Config(RAD) %use_alb_soil_litter
    use_alb_soil_organic_C    = dsl4jsb_Config(RAD) %use_alb_soil_organic_C
    use_alb_canopy            = dsl4jsb_Config(RAD) %use_alb_canopy
    albedo_age_weight         = dsl4jsb_Config(RAD) %albedo_age_weight
    AlbedoCanopySnow          = dsl4jsb_Config(RAD) %AlbedoCanopySnow

    ! Set pointers to variables in memory
    dsl4jsb_Get_var2D_onChunk(ATM2LAND, cos_zenith_angle)      ! in
    dsl4jsb_Get_var2D_onChunk(ATM2LAND, swvis_srf_down)        ! in
    dsl4jsb_Get_var2D_onChunk(ATM2LAND, swnir_srf_down)        ! in
    dsl4jsb_Get_var2D_onChunk(SEB,      t)                     ! in

    !RMS dsl4jsb_Get_var2D_onChunk(CARBON,   Cpool_litter_green_ag) ! in
    !RMS dsl4jsb_Get_var2D_onChunk(CARBON,   Cpool_slow)            ! in

    dsl4jsb_Get_var2D_onChunk(RAD,  swvis_srf_net)         ! in
    dsl4jsb_Get_var2D_onChunk(RAD,  swnir_srf_net)         ! in
    dsl4jsb_Get_var2D_onChunk(RAD,  alb_background)        ! in
    dsl4jsb_Get_var2D_onChunk(RAD,  alb_vis)               ! inout
    dsl4jsb_Get_var2D_onChunk(RAD,  alb_nir)               ! inout
    dsl4jsb_Get_var2D_onChunk(RAD,  alb_vis_snow)          ! inout
    dsl4jsb_Get_var2D_onChunk(RAD,  alb_nir_snow)          ! inout

    ! ---------------------------
    ! Go
    
    ! Set day-or-night-switch
    l_day(:) = Has_minimal_radiation(swvis_srf_down(:), swnir_srf_down(:))

    ! Calculate snow albedos. Used later only for bare and vegetation tile types
    ! Note, for vegetation tiles this snow albedos are used only for snow on soil
    IF ( (tile%lcts(1)%id == BARE_TYPE) .OR. (tile%lcts(1)%id == VEG_TYPE) ) THEN

      dsl4jsb_Get_var2D_onChunk(HYDRO, fract_snow_soil)  ! in
      dsl4jsb_Get_var2D_onChunk(RAD,   snow_age)         ! in

      CALL calc_snow_albedo(   &
        & l_day(:),            & ! in
        & fract_snow_soil(:),  & ! in
        & snow_age(:),         & ! in
        & cos_zenith_angle(:), & ! in
        & albedo_age_weight,   & ! in
        & t(:),                & ! in
        & alb_vis_snow(:),     & ! inout
        & alb_nir_snow(:)      & ! inout
        & )
    END IF


    IF (tile%lcts(1)%id == LAKE_TYPE) THEN
      dsl4jsb_Get_var2D_onChunk(SEB,      t_lice)       ! in
      dsl4jsb_Get_var2D_onChunk(SEB,      fract_lice)   ! in
      dsl4jsb_Get_var2D_onChunk(HYDRO,    w_snow_lice)  ! in
      dsl4jsb_Get_var2D_onChunk(ATM2LAND, albedo_lwtr)  ! in
      dsl4jsb_Get_var2D_onChunk(RAD,      albedo_lice)  ! inout

      CALL calc_alb_lice( &
        & l_day(:),       & ! in
        & t_lice(:),      & ! in
        & w_snow_lice(:), & ! in
        & albedo_lice(:)  & ! inout
        & )
      WHERE (l_day(:))
        alb_vis(:) = (1._wp - fract_lice(:)) * albedo_lwtr(:) + fract_lice(:) * albedo_lice(:)
        alb_nir(:) = alb_vis(:)
      END WHERE


    ELSE IF (tile%lcts(1)%id == GLACIER_TYPE) THEN
      CALL calc_glacier_albedo( &
        & l_day(:),             & ! in
        & t(:),                 & ! in       
        & alb_vis(:),           & ! inout
        & alb_nir(:)            & ! inout
        & )


    ELSE IF (tile%lcts(1)%id == BARE_TYPE) THEN
      dsl4jsb_Get_var2D_onChunk(RAD, alb_vis_mineralsoil)   ! inout
      dsl4jsb_Get_var2D_onChunk(RAD, alb_nir_mineralsoil)   ! inout
      dsl4jsb_Get_var2D_onChunk(RAD, alb_vis_soil)          ! inout
      dsl4jsb_Get_var2D_onChunk(RAD, alb_nir_soil)          ! inout

      IF (use_alb_soil_scheme) THEN
        CALL calc_soil_albedo(                        &
          & l_day(:),                                 & ! in
          & use_alb_mineralsoil_const,                & ! in
          & use_alb_soil_organic_C,                   & ! in
          & use_alb_soil_litter,                      & ! in
          & Cpool_litter_green_ag(:),                 & ! in
          & Cpool_slow(:),                            & ! in
          & dsl4jsb_Lctlib_param(specificLeafArea_C), & ! in
          & dsl4jsb_Lctlib_param(AlbedoLitterVIS),    & ! in
          & alb_vis_mineralsoil(:),                   & ! inout
          & alb_nir_mineralsoil(:),                   & ! inout
          & alb_vis_soil(:),                          & ! inout
          & alb_nir_soil(:)                           & ! inout
          & )
      ELSE
        ! Do nothing: alb_vis_soil and alb_nir_soil are already in memory initialited from bc_land_phy.nc file
        ! and represent the soil albedo including carbon in the soil as well as litter on the soil, but without
        ! the canopy above.
        IF (debug_on() .AND. iblk == 1) THEN
          CALL message(TRIM(routine), 'Soil albedo scheme not used for bare soil tiles but taken from bc_land_phy.nc file')
        END IF
      END IF

      ! Merge the albedo of soil with the snow on it for the overall bare soil tile albedo
      alb_vis = (fract_snow_soil(:) * alb_vis_snow(:)  + (1.0_wp - fract_snow_soil(:)) * alb_vis_soil(:))
      alb_nir = (fract_snow_soil(:) * alb_nir_snow(:)  + (1.0_wp - fract_snow_soil(:)) * alb_nir_soil(:))


    ELSE IF (tile%lcts(1)%id == VEG_TYPE) THEN

      dsl4jsb_Get_var2D_onChunk(PHENO,    lai)                   ! in
      dsl4jsb_Get_var2D_onChunk(VEG,      fract_forest)          ! in
      dsl4jsb_Get_var2D_onChunk(HYDRO,    fract_snow_can)        ! in
      dsl4jsb_Get_var2D_onChunk(HYDRO,    fract_snow_soil)       ! in
      dsl4jsb_Get_var2D_onChunk(RAD,      alb_vis_mineralsoil)   ! inout
      dsl4jsb_Get_var2D_onChunk(RAD,      alb_nir_mineralsoil)   ! inout      
      dsl4jsb_Get_var2D_onChunk(RAD,      alb_vis_can)           ! inout
      dsl4jsb_Get_var2D_onChunk(RAD,      alb_nir_can)           ! inout
      dsl4jsb_Get_var2D_onChunk(RAD,      alb_vis_soil)          ! inout
      dsl4jsb_Get_var2D_onChunk(RAD,      alb_nir_soil)          ! inout

      IF (use_alb_veg_simple) THEN
        alb_vis(:) = get_surface_albedo_simple( &
                       & t(:),                  &
                       & fract_forest(:),       &
                       & fract_snow_soil(:),    &
                       & fract_snow_can(:),     &
                       & lai(:),                &
                       & alb_background(:)      & ! Note, alb_background <-- bc_land_phy.nc
                       & )
        alb_nir(:) = alb_vis(:)
      ELSE

        ! Get canopy albedo
        IF (.NOT. use_alb_canopy) THEN                            ! TRUE:  nothing to do. The canopy albedo is already taken from
                                                                  !        the bc_land_phys.nc file during initialisation.
          alb_vis_can(:) = dsl4jsb_Lctlib_param(AlbedoCanopyVIS)  ! FALSE: The canopy albedo is taken from lctlib file.
          alb_nir_can(:) = dsl4jsb_Lctlib_param(AlbedoCanopyNIR)
        END IF

        ! Get bare soil albedo
        IF (use_alb_soil_scheme) THEN
          ! The albedo for bare soil under the canopy of a vegetation tile is the same as for the bare soil tile:
          CALL calc_soil_albedo(                        &
            & l_day(:),                                 & ! in
            & use_alb_mineralsoil_const,                & ! in
            & use_alb_soil_organic_C,                   & ! in
            & use_alb_soil_litter,                      & ! in
            & Cpool_litter_green_ag(:),                 & ! in
            & Cpool_slow(:),                            & ! in
            & dsl4jsb_Lctlib_param(specificLeafArea_C), & ! in
            & dsl4jsb_Lctlib_param(AlbedoLitterVIS),    & ! in
            & alb_vis_mineralsoil(:),                   & ! inout
            & alb_nir_mineralsoil(:),                   & ! inout
            & alb_vis_soil(:),                          & ! inout
            & alb_nir_soil(:)                           & ! inout
            & )
        ELSE
          ! Do nothing: alb_vis_soil and alb_nir_soil are already in memory initialited from the bc_land_phy.nc file
          ! and represent the soil albedo including carbon in the soil as well as litter on the soil, but without
          ! the canopy above.

          IF (debug_on() .AND. iblk == 1) THEN
            CALL message(TRIM(routine), 'Soil albedo scheme not used for veg tiles but taken from bc_land_phy.nc file')
          END IF
        END IF

        ! Get the fractions of the vegetation tile that is bare soil from the sky point of view
        ! Needed later for merging the contributing snow-, bare soil- and canopy-albedos for this vegetation tile
        dsl4jsb_Get_var2D_onChunk(RAD, sky_view_fract)       ! inout
        dsl4jsb_Get_var2D_onChunk(RAD, sky_view_fract_stem)  ! inout

        CALL calc_sky_view_fractions(       &
          & l_day(:),                       & ! in
          & lai(:),                         & ! in
          & cos_zenith_angle(:),            & ! in
          & dsl4jsb_Lctlib_param(StemArea), & ! in
          & sky_view_fract(:),              & ! out
          & sky_view_fract_stem(:)          & ! out
          & )

        CALL Merge_albedos_of_vegtile(        &
          & l_day(:),                         & ! in
          & dsl4jsb_Lctlib_param(ForestFlag), & ! in
          & sky_view_fract(:),                & ! in
          & sky_view_fract_stem(:),           & ! in
          & AlbedoCanopySnow,                 & ! in
          & fract_snow_soil(:),               & ! in
          & fract_snow_can(:),                & ! in
          & alb_vis_snow(:),                  & ! in
          & alb_nir_snow(:),                  & ! in
          & alb_vis_soil(:),                  & ! in
          & alb_nir_soil(:),                  & ! in
          & alb_vis_can(:),                   & ! in
          & alb_nir_can(:),                   & ! in
          & alb_vis(:),                       & ! inout
          & alb_nir(:)                        & ! inout
          & )
        
      END IF ! use_alb_veg_simple 
    END IF ! tile%lcts(1)%id == ...

  END SUBROUTINE update_albedo
  
  ! -------------------------------------------------------------------------------------------------------
  !>
  !! Implementation of "aggregate" for task "albedo"
  !!
  !! @param[in,out] tile    Tile for which routine is executed.
  !! @param[in]     options Additional run-time parameters.
  !!
  SUBROUTINE aggregate_albedo(tile, options)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    TYPE(t_jsb_model), POINTER :: model

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':aggregate_albedo'

    INTEGER :: iblk, ics, ice, nc

    iblk = options%iblk
    ics  = options%ics
    ice  = options%ice
    nc   = options%nc

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    model => Get_model(tile%owner_model_id)

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE aggregate_albedo

#endif
END MODULE mo_rad_interface
