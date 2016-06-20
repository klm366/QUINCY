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
  dsl4jsb_Use_processes SEB, SOIL, VEG, ATM2LAND, RAD, HYDRO

  ! Use process configurations
  dsl4jsb_Use_config(SEB)
  dsl4jsb_Use_config(SOIL)
  dsl4jsb_Use_config(RAD)

  ! Use of process memories
  dsl4jsb_Use_memory(ATM2LAND)
  dsl4jsb_Use_memory(SOIL)
  dsl4jsb_Use_memory(SEB)
  dsl4jsb_Use_memory(VEG)
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
  !> Constructor for radiation task
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

  ! -------------------------------------------------------------------------------------------------------
  !> Register tasks for radiation process
  !!
  !! @param[in,out] this      Instance of <PROCESS_NAME_LOWER_CASE> process class
  !! @param[in]     model_id  Model id
  !!
  SUBROUTINE Register_rad_tasks(this, model_id)

    CLASS(t_jsb_process), INTENT(inout) :: this
    INTEGER, INTENT(in) :: model_id

    CALL this%Register_task(tsk_surface_radiation(model_id))
    CALL this%Register_task(tsk_albedo           (model_id))

  END SUBROUTINE Register_rad_tasks

  ! ================================================================================================================================
  !! Surface radiation
  ! ================================================================================================================================

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

    ! Pointers to variables in memory
    dsl4jsb_Real2D_onChunk :: &
      & swvis_srf_net, &
      & swnir_srf_net, &
      & sw_srf_net,    &
      & lw_srf_net,    &
      & rad_srf_net,   &
      & lw_srf_down,   &
      & swvis_srf_down, &
      & swnir_srf_down, &
      & alb_vis,        &
      & alb_nir,        &
      & fract_lice,     &
      & rad_net_wtr,    &
      & rad_net_ice,    &
      & albedo_wtr,     &
      & albedo_ice,     &
      & t,              &
      & t_wtr,          &
      & t_ice

    REAL(wp) :: &
      & lw_srf_net_wtr(options%nc), &
      & lw_srf_net_ice(options%nc), &
      & sw_srf_net_wtr(options%nc), &
      & sw_srf_net_ice(options%nc)

    INTEGER :: iblk, ics, ice, nc

    TYPE(t_jsb_model), POINTER :: model

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':update_surface_radiation'

    iblk = options%iblk
    ics  = options%ics
    ice  = options%ice
    nc   = options%nc
    !
    !-----------------------------------------------------------------------
    !
    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    model => Get_model(tile%owner_model_id)

    dsl4jsb_Get_memory(ATM2LAND)
    dsl4jsb_Get_memory(RAD)
    dsl4jsb_Get_memory(SEB)

    dsl4jsb_Get_var2D_onChunk(ATM2LAND, swvis_srf_down) ! IN
    dsl4jsb_Get_var2D_onChunk(ATM2LAND, swnir_srf_down) ! IN
    dsl4jsb_Get_var2D_onChunk(ATM2LAND, lw_srf_down)    ! IN

    dsl4jsb_Get_var2D_onChunk(RAD,      sw_srf_net)     ! OUT
    dsl4jsb_Get_var2D_onChunk(RAD,      lw_srf_net)     ! OUT
    dsl4jsb_Get_var2D_onChunk(RAD,      rad_srf_net)    ! OUT

    IF (Contains_lct(tile%lcts, LAKE_TYPE)) THEN

      ! Water
      dsl4jsb_Get_var2D_onChunk(SEB,      t_wtr)       ! in
      dsl4jsb_Get_var2D_onChunk(ATM2LAND, albedo_wtr)  ! in
      dsl4jsb_Get_var2D_onChunk(RAD,      rad_net_wtr) ! out

      CALL calc_radiation_surface_net( &
        & swvis_srf_down(:),        & ! in
        & swnir_srf_down(:),        & ! in
        & albedo_wtr(:),            & ! in
        & albedo_wtr(:),            & ! in
        & lw_srf_down(:),           & ! in
        & t_wtr(:),                 & ! in
        & rad_net_wtr(:),           & ! out
        & lw_net=lw_srf_net_wtr(:), & ! out
        & sw_net=sw_srf_net_wtr(:)  & ! out
        & )

      ! Ice
      dsl4jsb_Get_var2D_onChunk(SEB, t_ice)        ! in
      dsl4jsb_Get_var2D_onChunk(SEB, fract_lice)   ! in
      dsl4jsb_Get_var2D_onChunk(RAD, albedo_ice)   ! in
      dsl4jsb_Get_var2D_onChunk(RAD, rad_net_ice)  ! out
      dsl4jsb_Get_var2D_onChunk(RAD, lw_srf_net)   ! out

      CALL calc_radiation_surface_net( &
        & swvis_srf_down(:),        & ! in
        & swnir_srf_down(:),        & ! in
        & albedo_ice(:),            & ! in
        & albedo_ice(:),            & ! in
        & lw_srf_down(:),           & ! in
        & t_ice(:),                 & ! in
        & rad_net_ice(:),           & ! out
        & lw_net=lw_srf_net_ice(:), & ! out
        & sw_net=sw_srf_net_ice(:)  & ! out
        & )

      rad_srf_net(:) = (1._wp - fract_lice(:)) * rad_net_wtr   (:) + fract_lice(:) * rad_net_ice   (:)
      lw_srf_net (:) = (1._wp - fract_lice(:)) * lw_srf_net_wtr(:) + fract_lice(:) * lw_srf_net_ice(:)
      sw_srf_net (:) = (1._wp - fract_lice(:)) * sw_srf_net_wtr(:) + fract_lice(:) * sw_srf_net_ice(:)

    ELSE

      dsl4jsb_Get_var2D_onChunk(SEB,      t)              ! IN
      dsl4jsb_Get_var2D_onChunk(RAD,      alb_vis)        ! IN
      dsl4jsb_Get_var2D_onChunk(RAD,      alb_nir)        ! IN
      dsl4jsb_Get_var2D_onChunk(RAD,      swvis_srf_net)  ! OUT
      dsl4jsb_Get_var2D_onChunk(RAD,      swnir_srf_net)  ! OUT

      CALL calc_radiation_surface_net( &
        & swvis_srf_down(:),    & ! input
        & swnir_srf_down(:),    & ! input
        & alb_vis(:),           & ! input
        & alb_nir(:),           & ! input
        & lw_srf_down(:),       & ! input
        & t(:),                 & ! input
        & rad_srf_net(:),       & ! output
        & swvis_srf_net(:),     & ! output
        & swnir_srf_net(:),     & ! output
        & sw_srf_net(:),        & ! output
        & lw_srf_net(:)         & ! output
        & )
    END IF

    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE update_surface_radiation

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

    dsl4jsb_Aggregate_onChunk(RAD, rad_net_wtr,    weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(RAD, rad_net_ice,    weighted_by_fract)

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE aggregate_surface_radiation

  ! ================================================================================================================================
  !! Photosynthetically active radiation
  ! ================================================================================================================================

  SUBROUTINE update_radiation_par(tile, options, waterLimitationFlag)

    ! ---------------------------
    ! Variables

    ! Used variables
    ! sollen jetzt in dieses File hier kommen:
    USE mo_rad_process,    ONLY: calc_fapar_lai_cl, calc_par

    ! Arguments
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options
    LOGICAL,                    INTENT(in)    :: waterLimitationFlag  ! Chooses between computation of water limited (true) and ...
                                                                      ! ... water unlimited (.FALSE.)
    !----------------------------
    ! Local variables

!    REAL(wp), POINTER, DIMENSION(:,:) :: apar_per_lai_cl   ! Absorbed PAR per canopy layer. If it can be removed from memory
                                                            ! it should be placed here.
!    REAL(wp), POINTER, DIMENSION(:)   :: par_down_mol      ! Downward PAR flux in mol (photons)/(m^2 s). If it can be removed
                                                            ! from memory it should be placed here.

    INTEGER                           :: icanopy

    INTEGER                           :: iblk, ics, ice, nc
    REAL(wp)                          :: dtime

    REAL(wp)                          :: canopy_bound_lai_delta

    REAL(wp), DIMENSION(options%nc)   :: B4_layer_above ! R: cannot be scalar as it will be used as inout argument for elemental procedure
                                                        !    and thus has to have the same dimensions as the other out arguments.

!    REAL(wp), DIMENSION(:)            :: frac_diffuse_par(nc), & ! in JSBACH3 received from atmosphere; however only used to calculate
                                                                 ! frac_par_direct, which is in JSBACH4 already in the "memory" => should be unnecessary and frac_par_direct should be received from the atmosphere...


    TYPE(t_jsb_model), POINTER :: model

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':update_par'

    dsl4jsb_Def_config(RAD)
    dsl4jsb_Def_memory(RAD)
    dsl4jsb_Def_config(SEB)
    dsl4jsb_Def_memory(SEB)
    dsl4jsb_Def_memory(VEG)
    dsl4jsb_Def_memory(ATM2LAND)

    ! Declare pointers of type variable in memory
    ! e.g. REAL(wp), POINTER, DIMENSION(:) variable
    INTEGER, POINTER       :: ncanopy               ! number of canopy layers
    INTEGER, POINTER       :: acc_counter           ! counter for each delta_time step since the last variable output
    LOGICAL, POINTER       :: use_alb_veg_simple    !

    dsl4jsb_Real2D_onChunk :: par                   ! Downward PAR flux in [W/m^2], Photosynth. Active Radiation (400-700nm)
    dsl4jsb_Real2D_onChunk :: par_down_mol          ! Downward PAR flux in mol (photons)/(m^2 s)
    dsl4jsb_Real2D_onChunk :: par_down_mol_nacc     ! PAR accumulated [MOL PHOTONS/M^2 ]
    dsl4jsb_Real2D_onChunk :: par_down_mol_tavg     ! PAR [MOL PHOTONS/M^2 S] (time mean value)
    dsl4jsb_Real2D_onChunk :: soil_reflectivity_par ! Soil reflectivity within the PAR spectrum [-]
    dsl4jsb_Real2D_onChunk :: alb_vis_soil          ! Soil albedo in the visible range
    dsl4jsb_Real2D_onChunk :: alb_background        ! USED ONLY FOR COMPARISON WITH JSBACH3
    dsl4jsb_Real2D_onChunk :: lai                   ! Leaf area index [-]
    dsl4jsb_Real2D_onChunk :: cos_zenith_angle      ! Angle of incoming radiation [-]
    dsl4jsb_Real2D_onChunk :: frac_par_direct       ! Fraction of direct (in contrast to diffuse) PAR [-]
    dsl4jsb_Real2D_onChunk :: canopy_bound_lai
    dsl4jsb_Real2D_onChunk :: apar_nacc

    dsl4jsb_Real3D_onChunk :: faPAR_cl              ! Fraction of absorbed PAR per leaf area for each canopy layer
    dsl4jsb_Real3D_onChunk :: apar_per_lai_cl       ! Absorbed PAR of canopy layer
    dsl4jsb_Real3D_onChunk :: lai_cl                ! Leaf area index of canopy layer
    dsl4jsb_Real3D_onChunk :: apar_per_lai_cl_nacc  ! Absorbed PAR accumulated [MOL PHOTONS/M^2 ]
    dsl4jsb_Real3D_onChunk :: apar_per_lai_cl_tavg  ! Absorbed PAR [MOL PHOTONS/M^2 S] (time mean value)

    ! ---------------------------
    ! Go

    iblk  = options%iblk
    ics   = options%ics
    ice   = options%ice
    nc    = options%nc
    dtime = options%dtime ! time step in seconds e.g. for one day=86400

    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting update of par')

    model => Get_model(tile%owner_model_id)

    ! Get reference to variables
    ! laut momentaner dsl:  CALL Get_process_config(config(PROC)%p, PROC__conf)
    ! müsste jetzt aber direkt aus Vektor config(:) genommen werden den ich oben daher als Argument übergebe ...
    dsl4jsb_Get_config(RAD)
    dsl4jsb_Get_memory(RAD)
    dsl4jsb_Get_memory(VEG)
    dsl4jsb_Get_memory(ATM2LAND)

    ! Set pointers
    !R: dsl fehlt Befehl um config Variablen zu laden...
    ! sollte jetzt z.B. ncanopy => config(4)%ncanopy sein ...
    !RMS ncanopy                  => model%config              %ncanopy
    !RMS canopy_bound_lai         => model%config              %canopy_bound_lai
    use_alb_veg_simple       => dsl4jsb_Config(RAD)       %use_alb_veg_simple

    !acc_counter              => tile%acc_counter

    dsl4jsb_Get_var2D_onChunk(RAD,      par)                   ! INOUT
    dsl4jsb_Get_var2D_onChunk(RAD,      par_down_mol)          ! INOUT

    ! R: could also be local as it is not used outside this procedure in the moment. However I want it in the model output and thus it has to be updated by returning it (as INOUT)
    dsl4jsb_Get_var2D_onChunk(RAD,      par_down_mol_nacc)     ! INOUT
    dsl4jsb_Get_var2D_onChunk(RAD,      par_down_mol_tavg)     ! INOUT
    dsl4jsb_Get_var2D_onChunk(RAD,      soil_reflectivity_par) ! INOUT
    dsl4jsb_Get_var2D_onChunk(RAD,      frac_par_direct)       ! IN
    dsl4jsb_Get_var2D_onChunk(RAD,      alb_vis_soil)          ! IN
    dsl4jsb_Get_var2D_onChunk(RAD,      alb_background)        ! IN
    dsl4jsb_Get_var2D_onChunk(VEG,      lai)                   ! INOUT
    dsl4jsb_Get_var2D_onChunk(ATM2LAND, cos_zenith_angle)      ! IN
    dsl4jsb_Get_var3D_onChunk(RAD,      lai_cl)                ! INOUT
    dsl4jsb_Get_var3D_onChunk(RAD,      faPAR_cl)              ! INOUT
    ! used outside for assimilation:
    dsl4jsb_Get_var3D_onChunk(RAD,      apar_per_lai_cl)       ! INOUT
    dsl4jsb_Get_var3D_onChunk(RAD,      apar_per_lai_cl_nacc)  ! INOUT
    dsl4jsb_Get_var3D_onChunk(RAD,      apar_per_lai_cl_tavg)  ! INOUT
    dsl4jsb_Get_var2D_onChunk(RAD,      apar_nacc)             ! INOUT

    ! initialize local arguments and variables
    par_down_mol           = 0._wp ! for now also in the memory but is only for current time step/"delta_time"
    soil_reflectivity_par  = 0._wp ! for now also in the memory but is only for current time step/"delta_time"
    !apar_per_lai_cl        = 0._wp ! local, only for current time step/dtime. Used outside for assimilation


    B4_layer_above = 7777777 ! Used to transfer variable B4 from one canopy layer to the next. Startvalue only to indicate first layer
    DO icanopy=1,ncanopy
       ! Compute faPAR_cl (absorbed photsynthetic active radiation for each canopy layer) and lai_cl (leaf area index for each canopy layer)
       ! from LAI, soil_reflectivity_par, cos_zenith_angle, frac_par_direct, canopy_boundaries_lai.
       canopy_bound_lai_delta = canopy_bound_lai(icanopy) - canopy_bound_lai(icanopy-1)

       CALL calc_par(dtime,                     & ! input
                     ncanopy,                   & ! input
                     icanopy,                   & ! input
                     waterLimitationFlag,       & ! input
                     use_alb_veg_simple,        & ! input
                     par(:),                    & ! input
                     alb_vis_soil(:),           & ! input
                     lai(:),                    & ! input
                     cos_zenith_angle(:),       & ! input
                     frac_par_direct(:),        & ! input
                     canopy_bound_lai(icanopy), & ! input
                     canopy_bound_lai_delta,    & ! input
                     B4_layer_above,            & ! in and out
                     par_down_mol,              & ! output, used outside for assimilation
                     soil_reflectivity_par,     & ! output, R: could also be local as it is not used outside in the moment,
                                                  ! however I want it in the model output and thus it has to be updated by returning it
                     faPAR_cl(:,icanopy),       & ! output, R: could also be local as it is not used outside in the moment,
                                                  ! however I want it in the model output and thus it has to be updated by returning it
                     lai_cl(:,icanopy),         & ! output, R: could also be local as it is not used outside in the moment,
                                                  ! however I want it in the model output and thus it has to be updated by returning it
                     apar_per_lai_cl(:,icanopy) ) ! output, used outside for assimilation
    END DO

    ! Accumulate incoming par in mol(photons)/(m^2) over time
    IF (.NOT. waterLimitationFlag) THEN ! to make sure to accumulate only once per timestep
       par_down_mol_nacc = par_down_mol_nacc + par_down_mol * dtime ! Change explanation:
                                                                          ! "par/Epar" was unnecessarily double calculated. In JSBACH3 
                                                                          ! "CALL add(..., laccu=.TRUE.)" accumulates variable and sets it to zero
                                                                          ! after writen out. As "laccu" does not exist in JSBACH4 (it could be
                                                                          ! added) I added a workaround.
                                                                          ! acc_counter, par_down_mol_nacc, apar_per_lai_cl_nacc is set to zero
                                                                          ! when the output is written.
       ! R: This whole workaround should get abitary soon...
       !acc_counter = acc_counter + 1
       !par_down_mol_tavg = par_down_mol_nacc / dtime / acc_counter
    END IF

    ! R: Compute direct fraction of PAR. Unneccessary as it should be taken directly from atmosphere (frac_par_direct = 1-wp - jsswdifpar)
    ! frac_par_direct = 1._wp - frac_par_diffuse ! Taken from mo_jsbach_interface.f90 in JSBACH3.
                                                 ! frac_par_direct =sw_par_fract_direct
                                                 ! frac_diffuse_par=zsw_par_frac_diffuse <- von ECHAM: jsswdifpar

    IF (waterLimitationFlag) THEN ! to make sure to accumulate only once per timestep
       apar_per_lai_cl_nacc  = apar_per_lai_cl_nacc + apar_per_lai_cl * dtime
       ! R: This whole workaround (acc_counter) should get abitary soon...
       !apar_per_lai_cl_tavg(:,icanopy)  = apar_per_lai_cl_nacc(:,icanopy) / dtime / acc_counter
       apar_nacc                     = apar_nacc + par_down_mol * faPAR_cl(:,icanopy) * dtime
    END IF

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

  END SUBROUTINE update_radiation_par

  ! ================================================================================================================================
  !! Albedo
  ! ================================================================================================================================

  SUBROUTINE update_albedo(tile, options)

    USE mo_rad_process,    ONLY: get_surface_albedo_simple, calc_soil_albedo, calc_albedo_ice, calc_glacier_albedo, &
                                 calc_sky_view_fractions, calc_snow_albedo, Merge_albedos, &
                                 Has_minimal_radiation
    USE mo_rad_constants,  ONLY: AlbedoCanopySnow_age, AlbedoCanopySnow_temp

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    ! Declare pointers to process configs and memory
    dsl4jsb_Def_config(RAD)
    dsl4jsb_Def_memory(RAD)
    dsl4jsb_Def_memory(SEB)
    dsl4jsb_Def_memory(ATM2LAND)
    !RMS dsl4jsb_Def_memory(CARBON)
    dsl4jsb_Def_memory(HYDRO)
    dsl4jsb_Def_memory(VEG)

    ! Declare pointers to variables in memory
    dsl4jsb_Real2D_onChunk :: &
      & cos_zenith_angle,     &
      & swvis_srf_down,       &
      & swnir_srf_down,       &
      & t,                    &
      & alb_vis,              &
      & alb_nir,              &
      & alb_vis_can,          &
      & alb_nir_can,          &
      & alb_vis_soil,         &
      & alb_nir_soil,         &
      & alb_vis_snow,         &
      & alb_nir_snow,         &
      & alb_background,       &
      & alb_vis_mineralsoil,  &
      & alb_nir_mineralsoil,  &
      & sky_view_fract,       &
      & sky_view_fract_stem,  &
      & snow_age,             &
      & lai,                  &
      & fract_forest,         &
      & fract_snow_can,       &
      & fract_snow_soil,      &
      & Cpool_litter_green_ag,&
      & Cpool_slow,           &
      & swvis_srf_net,        &
      & swnir_srf_net
    ! Lakes
    dsl4jsb_Real2D_onChunk :: &
      & albedo_wtr,  &
      & fract_lice,  &
      & t_ice,       &
      & w_snow_lice, &
      & albedo_ice

    !----------------------------
    ! Local variables

    LOGICAL,  POINTER :: &
      & use_alb_veg_simple,        &
      & use_alb_soil_scheme,       &
      & use_alb_soil_litter,       &
      & use_alb_mineralsoil_const, &
      & use_alb_canopy

    LOGICAL :: l_day(options%nc)

    CHARACTER(len=10), POINTER :: use_alb_soil_organic_C

    REAL(wp), POINTER :: albedo_age_weight    ! 0 < albedo_age_weight < 1: snow albedo is calcualted by linearly weighting
                                              !                            the snow albedo resulting from both schemes
    REAL(wp)          :: AlbedoCanopySnow


    INTEGER  :: iblk, ics, ice, nc
    REAL(wp) :: dtime

    TYPE(t_jsb_model), POINTER :: model

    CHARACTER(len=*), PARAMETER    :: routine = TRIM(modname)//':update_albedo'

    iblk  = options%iblk
    ics   = options%ics
    ice   = options%ice
    nc    = options%nc
    dtime = options%dtime ! time step in seconds e.g. for one day=86400

    model => Get_model(tile%owner_model_id)

    ! Get pointers to process configs and memory
    dsl4jsb_Get_config(RAD)
    dsl4jsb_Get_memory(RAD)
    dsl4jsb_Get_memory(ATM2LAND)
    dsl4jsb_Get_memory(SEB)
    !RMS dsl4jsb_Get_memory(CARBON)
    dsl4jsb_Get_memory(HYDRO)
    IF (Contains_lct(tile%lcts, VEG_TYPE)) THEN
      dsl4jsb_Get_memory(VEG)
    END IF

    use_alb_veg_simple        => dsl4jsb_Config(RAD)       %use_alb_veg_simple
    use_alb_soil_scheme       => dsl4jsb_Config(RAD)       %use_alb_soil_scheme
    use_alb_mineralsoil_const => dsl4jsb_Config(RAD)       %use_alb_mineralsoil_const
    use_alb_soil_litter       => dsl4jsb_Config(RAD)       %use_alb_soil_litter
    use_alb_soil_organic_C    => dsl4jsb_Config(RAD)       %use_alb_soil_organic_C
    use_alb_canopy            => dsl4jsb_Config(RAD)       %use_alb_canopy
    albedo_age_weight         => dsl4jsb_Config(RAD)       %albedo_age_weight

    ! Set pointers to variables in memory
    dsl4jsb_Get_var2D_onChunk(ATM2LAND, cos_zenith_angle)      ! IN
    dsl4jsb_Get_var2D_onChunk(ATM2LAND, swvis_srf_down)        ! IN
    dsl4jsb_Get_var2D_onChunk(ATM2LAND, swnir_srf_down)        ! IN


    dsl4jsb_Get_var2D_onChunk(SEB,      t)                     ! IN

    !RMS dsl4jsb_Get_var2D_onChunk(CARBON,   Cpool_litter_green_ag) ! IN
    !RMS dsl4jsb_Get_var2D_onChunk(CARBON,   Cpool_slow)            ! IN

    dsl4jsb_Get_var2D_onChunk(RAD,      swvis_srf_net)         ! IN
    dsl4jsb_Get_var2D_onChunk(RAD,      swnir_srf_net)         ! IN
    dsl4jsb_Get_var2D_onChunk(RAD,      alb_background)        ! IN

    dsl4jsb_Get_var2D_onChunk(RAD,      alb_vis)               ! INOUT
    dsl4jsb_Get_var2D_onChunk(RAD,      alb_nir)               ! INOUT
    dsl4jsb_Get_var2D_onChunk(RAD,      alb_vis_snow)          ! INOUT
    dsl4jsb_Get_var2D_onChunk(RAD,      alb_nir_snow)          ! INOUT

    l_day(:) = Has_minimal_radiation(swvis_srf_down(:), swnir_srf_down(:))

    IF (Contains_lct(tile%lcts, LAKE_TYPE)) THEN

      dsl4jsb_Get_var2D_onChunk(SEB,      t_ice)        ! in
      dsl4jsb_Get_var2D_onChunk(SEB,      fract_lice)   ! in
      dsl4jsb_Get_var2D_onChunk(HYDRO,    w_snow_lice)  ! in
      dsl4jsb_Get_var2D_onChunk(ATM2LAND, albedo_wtr)   ! in
      dsl4jsb_Get_var2D_onChunk(RAD,      albedo_ice)   ! out

      CALL calc_albedo_ice( &
        & t_ice(:),         &
        & w_snow_lice(:),   &
        & albedo_ice(:),    &
        & l_day(:)          &
        & )
      WHERE (l_day)
        alb_vis(:) = (1._wp - fract_lice(:)) * albedo_wtr(:) + fract_lice(:) * albedo_ice(:)
        alb_nir(:) = alb_vis(:)
      END WHERE

    ELSE IF (Contains_lct(tile%lcts, GLACIER_TYPE)) THEN

!        IF (use_alb_veg_simple) THEN
!           alb_vis(:) =albedo_snow_grnd aus Prozedur get_surface_albedo_simple
!           alb_nir(:) =albedo_snow_grnd aus Prozedur get_surface_albedo_simple
!           ! R: Die Prozedur get_surface_albedo_simple mit Glacier-maske auszustatten wäre simpel aber ich finde das unschön.
!                Die Frage ist ob die Prozedur wirklich benötigt wird und ob für Glacier dasselbe rauskommen muss...
!        ELSE

      CALL calc_glacier_albedo( &
        & t(:),           & ! in
        & alb_vis(:),     & ! inout
        & alb_nir(:),     & ! inout
        & l_day           & ! in
        & )

!        END IF
      !
    ELSE IF (Contains_lct(tile%lcts, BARE_TYPE)) THEN
!
      dsl4jsb_Get_var2D_onChunk(RAD,      alb_vis_mineralsoil)   ! INOUT
      dsl4jsb_Get_var2D_onChunk(RAD,      alb_nir_mineralsoil)   ! INOUT
      dsl4jsb_Get_var2D_onChunk(RAD,      alb_vis_soil)          ! INOUT
      dsl4jsb_Get_var2D_onChunk(RAD,      alb_nir_soil)          ! INOUT


!         CALL calc_soil_albedo(use_alb_mineralsoil_const, &
!                               use_alb_soil_organic_C,    &
!                               use_alb_soil_litter,       &
!                               Cpool_litter_green_ag(:),  &
!                               Cpool_slow(:),             &
!                               alb_vis_mineralsoil(:),    & ! Inout
!                               alb_nir_mineralsoil(:),    & ! Inout
!                               alb_vis_soil(:),           & ! Inout
!                               alb_nir_soil(:) )            ! Inout
!
!         alb_vis(:) =alb_vis_soil(:)
!         alb_nir(:) =alb_nir_soil(:)

    ELSE IF (Contains_lct(tile%lcts, VEG_TYPE)) THEN

      dsl4jsb_Get_var2D_onChunk(RAD,      alb_vis_mineralsoil)   ! INOUT
      dsl4jsb_Get_var2D_onChunk(RAD,      alb_nir_mineralsoil)   ! INOUT
      dsl4jsb_Get_var2D_onChunk(VEG,      lai)                   ! IN
      dsl4jsb_Get_var2D_onChunk(VEG,      fract_forest)          ! IN
      dsl4jsb_Get_var2D_onChunk(HYDRO,    fract_snow_can)        ! IN
      dsl4jsb_Get_var2D_onChunk(HYDRO,    fract_snow_soil)       ! IN
      dsl4jsb_Get_var2D_onChunk(RAD,      alb_vis_can)           ! INOUT
      dsl4jsb_Get_var2D_onChunk(RAD,      alb_nir_can)           ! INOUT
      dsl4jsb_Get_var2D_onChunk(RAD,      alb_vis_soil)          ! INOUT
      dsl4jsb_Get_var2D_onChunk(RAD,      alb_nir_soil)          ! INOUT



 !!!       IF (use_alb_veg_simple) THEN

           ! if tile is lowest child tile
              ! Einfache Berechnung der Albedo eines Vegetation-Tiles mit der Funktion canopy_albedo_echam5:
      alb_vis(:) = get_surface_albedo_simple( &
        & t(:), fract_forest(:), fract_snow_soil(:), fract_snow_can(:), &
        & lai(:), alb_background(:))
      alb_nir(:) = alb_vis(:)


! !         ELSE
! !            IF (.NOT. use_alb_canopy) THEN ! use_alb_canopy =TRUE: nothing to do. Canopy albedo is already taken from bc_land_phys.nc file during initialisation.
! !                                           !                =FALSE: canopy albedo is taken from lctlib file...
! ! !LCT              alb_vis_can = lctlib%AlbedoCanopyVIS(cover_type(i,itile))
! ! !LCT              alb_nir_can = lctlib%AlbedoCanopyNIR(cover_type(i,itile))
! !            END IF
! !
! !            IF (use_alb_soil_scheme) THEN
! !
! !               CALL calc_soil_albedo(use_alb_mineralsoil_const, & ! In
! !                                     use_alb_soil_organic_C,    & ! In
! !                                     use_alb_soil_litter,       & ! In
! !                                     Cpool_litter_green_ag(:),  & ! In
! !                                     Cpool_slow(:),             & ! In
! !                                     alb_vis_mineralsoil(:),    & ! InOut
! !                                     alb_nir_mineralsoil(:),    & ! InOut
! !                                     alb_vis_soil(:),           & ! InOut
! !                                     alb_nir_soil(:),           & ! InOut
! !                                     dsl4jsb_Lctlib_param(specificLeafArea_C), &
! !                                     dsl4jsb_Lctlib_param(AlbedoLitterVIS) &
! !                                     )
! !
! !            ELSE
! !               ! Do nothing: alb_vis_soil and alb_nir_soil are already in memory (e.g. from bc_land_phy.nc file)
! !               ! R: In JSBACH3: background_albedo_soil_vis(:,:) = background_albedo_soil_vis_base(:,:) <-- theLand%Soil%albedo_soil_nir 
! !               !    <-- from init_soil from jsbach.nc infile from variable albedo_soil_nir
! !
! !               IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), &
! !                       'Soil albedo scheme not used. Soil albedo read in from jsbach.nc file')
! !
! !            END IF ! from use_alb_soil_scheme

    END IF

    IF (Contains_lct(tile%lcts, BARE_TYPE) &
      & .OR. Contains_lct(tile%lcts, VEG_TYPE)) THEN

      dsl4jsb_Get_var2D_onChunk(HYDRO,    fract_snow_soil)       ! IN
      dsl4jsb_Get_var2D_onChunk(RAD,      snow_age)              ! IN

      CALL calc_snow_albedo( &
        fract_snow_soil(:),  & ! In
        snow_age(:),         & ! In
        cos_zenith_angle(:), & ! In
        albedo_age_weight,   & ! In
        t(:),                & ! In
        alb_vis_snow(:),     & ! Out
        alb_nir_snow(:)      & ! Out
        )

    END IF

    IF (Contains_lct(tile%lcts, VEG_TYPE)) THEN

      dsl4jsb_Get_var2D_onChunk(RAD,      sky_view_fract)        ! INOUT
      dsl4jsb_Get_var2D_onChunk(RAD,      sky_view_fract_stem)   ! INOUT


      CALL calc_sky_view_fractions(lai(:),                & ! In
        & cos_zenith_angle(:),   & ! In
        & sky_view_fract(:),     & ! Out
        & sky_view_fract_stem(:), & ! Out
        & dsl4jsb_Lctlib_param(StemArea) &
        & )

      ! R: All inputs are constants, not necessary to calculate this for each grid cell again. Therefore here
      !    and not in merge_albedos.               
      ! Weight between age and temperature impact on snow albedo for snow on canopy
      AlbedoCanopySnow = albedo_age_weight * AlbedoCanopySnow_age + (1._wp - albedo_age_weight) * AlbedoCanopySnow_temp

      CALL Merge_albedos( &
        & sky_view_fract(:),          & ! in
        & sky_view_fract_stem(:),     & ! in
        & AlbedoCanopySnow,           & ! in
        & fract_snow_soil(:),         & ! in
        & fract_snow_can(:),          & ! in
        & alb_vis_snow(:),            & ! in
        & alb_nir_snow(:),            & ! in
        & alb_vis_soil(:),            & ! in
        & alb_nir_soil(:),            & ! in
        & alb_vis_can(:),             & ! in
        & alb_nir_can(:),             & ! in
        & alb_vis(:),                 & ! out
        & alb_nir(:),                 & ! out
        & dsl4jsb_Lctlib_param(ForestFlag), & ! in
        & l_day(:)                    &
        & )

    END IF

! !         END IF ! from use_alb_veg_simple

  END SUBROUTINE update_albedo

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
