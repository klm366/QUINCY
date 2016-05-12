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
  USE mo_jsb_control,        ONLY: debug_on
  USE mo_kind,               ONLY: wp
  USE mo_exception,          ONLY: message, finish

  USE mo_jsb_model_class,    ONLY: t_jsb_model
  USE mo_jsb_class,          ONLY: Get_model
  USE mo_jsb_tile_class,     ONLY: t_jsb_tile_abstract
  USE mo_jsb_lct_class,      ONLY: LAND_TYPE, LAKE_TYPE, GLACIER_TYPE, BARE_TYPE, VEG_TYPE, Contains_lct
  USE mo_jsb_config_class,   ONLY: t_jsb_config, t_jsb_config_p
  USE mo_jsb_process_class,  ONLY: t_jsb_process
  USE mo_jsb_task_class,     ONLY: t_jsb_process_task, t_jsb_task_options

  ! Use of processes in this module
  dsl4jsb_Use_processes SEB, SOIL, VEG, ATM2LAND, RAD, HYDRO

  ! Use process configurations
  dsl4jsb_Use_config(SEB)
  dsl4jsb_Use_config(SOIL)
  dsl4jsb_Use_config(RAD)

  ! Use of process memories
  USE mo_atmland_memory_class, ONLY: t_atm2land_memory
!  dsl4jsb_Use_memory(ATM2LAND)
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
  TYPE, EXTENDS(t_jsb_process_task) :: tsk_radiation_surface
  CONTAINS
    PROCEDURE, NOPASS :: Integrate => update_radiation_surface
    PROCEDURE, NOPASS :: Aggregate => aggregate_radiation_surface
  END TYPE tsk_radiation_surface

  INTERFACE tsk_radiation_surface
    PROCEDURE Create_task_radiation_surface
  END INTERFACE tsk_radiation_surface

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
  FUNCTION Create_task_radiation_surface(model_id) RESULT(return_ptr)

    INTEGER,                   INTENT(in) :: model_id
    CLASS(t_jsb_process_task), POINTER    :: return_ptr

    ALLOCATE(tsk_radiation_surface::return_ptr)

    return_ptr%name = 'radiation_surface'
    return_ptr%process_id = RAD
    return_ptr%owner_model_id = model_id
! !!$    return_ptr%Integrate => update_radiation_surface
! !!$    return_ptr%Aggregate => aggregate_radiation_surface

  END FUNCTION Create_task_radiation_surface

  ! -------------------------------------------------------------------------------------------------------
  !> Register tasks for radiation process
  !!
  !! @param[in,out] this      Instance of <PROCESS_NAME_LOWER_CASE> process class
  !! @param[in]     model_id  Model id
  !!
  SUBROUTINE Register_rad_tasks(this, model_id)

    CLASS(t_jsb_process), INTENT(inout) :: this
    INTEGER, INTENT(in) :: model_id

    CALL this%Register_task(tsk_Radiation_surface(model_id))

  END SUBROUTINE Register_rad_tasks

  SUBROUTINE update_radiation_surface(tile, options)

    !-----------------------------------------------------------------------
    !  DECLARATIONS
    USE mo_rad_process, ONLY: calc_radiation_surface_net

    !-----------------------------------------------------------------------
    !  ARGUMENTS
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options), INTENT(in) :: options

    !-----------------------------------------------------------------------
    ! POINTERS TO VARIABLES IN MEMORY

    dsl4jsb_Def_config(RAD)
    dsl4jsb_Def_memory(ATM2LAND)
    dsl4jsb_Def_memory(SEB)
    dsl4jsb_Def_memory(RAD)
    dsl4jsb_Def_memory(VEG)

    ! Pointers to variables in memory
    dsl4jsb_Real2D_onChunk :: swvis_srf_net
    dsl4jsb_Real2D_onChunk :: swnir_srf_net
    dsl4jsb_Real2D_onChunk :: sw_srf_net
    dsl4jsb_Real2D_onChunk :: lw_srf_net
    dsl4jsb_Real2D_onChunk :: rad_srf_net
    dsl4jsb_Real2D_onChunk :: lw_srf_down
    dsl4jsb_Real2D_onChunk :: swvis_srf_down
    dsl4jsb_Real2D_onChunk :: swnir_srf_down
    dsl4jsb_Real2D_onChunk :: alb_vis
    dsl4jsb_Real2D_onChunk :: alb_nir
    dsl4jsb_Real2D_onChunk :: t

    INTEGER :: iblk, ics, ice, nc

    TYPE(t_jsb_model), POINTER :: model

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':update_radiation_surface'

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

    dsl4jsb_Get_var2D_onChunk(RAD,      swvis_srf_net)  ! OUT
    dsl4jsb_Get_var2D_onChunk(RAD,      swnir_srf_net)  ! OUT
    dsl4jsb_Get_var2D_onChunk(RAD,      sw_srf_net)     ! OUT
    dsl4jsb_Get_var2D_onChunk(RAD,      lw_srf_net)     ! OUT
    dsl4jsb_Get_var2D_onChunk(RAD,      rad_srf_net)    ! OUT
    dsl4jsb_Get_var2D_onChunk(SEB,      t)              ! IN
    dsl4jsb_Get_var2D_onChunk(RAD,      alb_vis)        ! IN
    dsl4jsb_Get_var2D_onChunk(RAD,      alb_nir)        ! IN

    IF (Contains_lct(tile%lcts, LAKE_TYPE)) THEN
    END IF
    CALL calc_radiation_surface_net( &
      & swvis_srf_down(:),    & ! input
      & swnir_srf_down(:),    & ! input
      & alb_vis(:),           & ! input
      & alb_nir(:),           & ! input
      & lw_srf_down(:),       & ! input
      & t(:),                 & ! input
      & swvis_srf_net(:),     & ! output
      & swnir_srf_net(:),     & ! output
      & sw_srf_net(:),        & ! output
      & lw_srf_net(:),        & ! output
      & rad_srf_net(:) )        ! output

    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE update_radiation_surface

  SUBROUTINE aggregate_radiation_surface(tile, options)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    TYPE(t_jsb_model), POINTER :: model

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':aggregate_radiation_surface'

    INTEGER :: iblk, ics, ice, nc

    iblk = options%iblk
    ics  = options%ics
    ice  = options%ice
    nc   = options%nc

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    model => Get_model(tile%owner_model_id)

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE aggregate_radiation_surface

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

!    REAL(wp), POINTER, DIMENSION(:,:) :: apar_per_lai_cl       ! Absorbed PAR per canopy layer. If it can be removed from memory it should be placed here.
!    REAL(wp), POINTER, DIMENSION(:)   :: par_down_mol  ! Downward PAR flux in mol (photons)/(m^2 s). If it can be removed from memory it should be placed here.

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

    ! Declare pointers of type process __conf and __mem
    ! TYPE(t_radiation_config), POINTER :: radiation__conf
    ! TYPE(t_radiation_memory), POINTER :: radiation__mem
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


  ! R: The principal sequence of calls for a process in the new JSBACH4 scheme should be similar for each tile. However, the schemes for
  !    the albedo can not be called for each tile in the same sequence of calls/procedures only with different tile specific constants/properties
  !    (like alb_vis_mineralsoil or snow_fract). Therefore it is still necessary to differentiate the albedo schemes for each tile.

  SUBROUTINE update_radiation_albedo( tile, options)

    ! ---------------------------
    ! Variables

    ! Used variables
    ! sollen jetzt in dieses File hier kommen:
    USE mo_rad_process,    ONLY: get_surface_albedo_simple, calc_soil_albedo, calc_glacier_albedo, &
                                       calc_sky_view_fractions, calc_snow_albedo, calc_merge_albedos
    USE mo_rad_constants,  ONLY: AlbedoCanopySnow_age, AlbedoCanopySnow_temp

    ! Arguments
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    !----------------------------
    ! Local variables

    REAL(wp)                        :: AlbedoCanopySnow


    ! Declare pointers of type process __conf and __mem
    ! TYPE(t_radiation_config), POINTER :: radiation__conf
    ! TYPE(t_radiation_memory), POINTER :: radiation__mem
    dsl4jsb_Def_config(RAD)
    dsl4jsb_Def_memory(RAD)
    dsl4jsb_Def_memory(SEB)
    dsl4jsb_Def_memory(ATM2LAND)
    !RMS dsl4jsb_Def_memory(CARBON)
    dsl4jsb_Def_memory(HYDRO)
    dsl4jsb_Def_memory(VEG)

    ! Declare pointers of type variable in memory
    ! e.g. REAL(wp), POINTER, DIMENSION(:) variable
    LOGICAL,  POINTER      :: use_alb_veg_simple,        &
                              use_alb_soil_scheme,       &
                              use_alb_soil_litter,       &
                              use_alb_mineralsoil_const, &
                              use_alb_canopy

    CHARACTER(len=10), POINTER :: use_alb_soil_organic_C

    REAL(wp), POINTER :: albedo_age_weight    ! 0 < albedo_age_weight < 1: snow albedo is calcualted by linearly weighting
                                              !                            the snow albedo resulting from both schemes

    INTEGER  :: iblk, ics, ice, nc
    REAL(wp) :: dtime

    dsl4jsb_Real2D_onChunk :: cos_zenith_angle
    dsl4jsb_Real2D_onChunk :: t
    dsl4jsb_Real2D_onChunk :: alb_vis
    dsl4jsb_Real2D_onChunk :: alb_nir
    dsl4jsb_Real2D_onChunk :: alb_vis_can
    dsl4jsb_Real2D_onChunk :: alb_nir_can
    dsl4jsb_Real2D_onChunk :: alb_vis_soil
    dsl4jsb_Real2D_onChunk :: alb_nir_soil
    dsl4jsb_Real2D_onChunk :: alb_vis_snow
    dsl4jsb_Real2D_onChunk :: alb_nir_snow
    dsl4jsb_Real2D_onChunk :: alb_background
    dsl4jsb_Real2D_onChunk :: alb_vis_mineralsoil
    dsl4jsb_Real2D_onChunk :: alb_nir_mineralsoil
    dsl4jsb_Real2D_onChunk :: sky_view_fract
    dsl4jsb_Real2D_onChunk :: sky_view_fract_stem
    dsl4jsb_Real2D_onChunk :: snow_age
    dsl4jsb_Real2D_onChunk :: lai
    dsl4jsb_Real2D_onChunk :: fract_forest
    dsl4jsb_Real2D_onChunk :: fract_snow_can
    dsl4jsb_Real2D_onChunk :: fract_snow_soil
    dsl4jsb_Real2D_onChunk :: Cpool_litter_green_ag
    dsl4jsb_Real2D_onChunk :: Cpool_slow
    dsl4jsb_Real2D_onChunk :: swvis_srf_net
    dsl4jsb_Real2D_onChunk :: swnir_srf_net

    TYPE(t_jsb_model), POINTER :: model

    CHARACTER(len=*), PARAMETER    :: routine = TRIM(modname)//':update_radiation_albedo'

    ! ---------------------------
    ! Go

    iblk  = options%iblk
    ics   = options%ics
    ice   = options%ice
    nc    = options%nc
    dtime = options%dtime ! time step in seconds e.g. for one day=86400

    model => Get_model(tile%owner_model_id)

    ! Get reference to variables
    ! laut momentaner dsl:  CALL Get_process_config(config(PROC)%p, PROC__conf)
    ! müsste jetzt aber direkt aus Vektor config(:) genommen werden den ich oben daher als Argument übergebe ...
    dsl4jsb_Get_config(RAD)
    dsl4jsb_Get_memory(RAD)
    dsl4jsb_Get_memory(ATM2LAND)
    dsl4jsb_Get_memory(SEB)
    !RMS dsl4jsb_Get_memory(CARBON)
    dsl4jsb_Get_memory(HYDRO)
    dsl4jsb_Get_memory(VEG)

    ! Set pointers
    use_alb_veg_simple        => dsl4jsb_Config(RAD)       %use_alb_veg_simple
    use_alb_soil_scheme       => dsl4jsb_Config(RAD)       %use_alb_soil_scheme
    use_alb_mineralsoil_const => dsl4jsb_Config(RAD)       %use_alb_mineralsoil_const
    use_alb_soil_litter       => dsl4jsb_Config(RAD)       %use_alb_soil_litter
    use_alb_soil_organic_C    => dsl4jsb_Config(RAD)       %use_alb_soil_organic_C
    use_alb_canopy            => dsl4jsb_Config(RAD)       %use_alb_canopy
    albedo_age_weight         => dsl4jsb_Config(RAD)       %albedo_age_weight

    dsl4jsb_Get_var2D_onChunk(ATM2LAND, cos_zenith_angle)      ! IN
    dsl4jsb_Get_var2D_onChunk(SEB,      t)                     ! IN
    dsl4jsb_Get_var2D_onChunk(VEG,      lai)                   ! IN
    dsl4jsb_Get_var2D_onChunk(VEG,      fract_forest)          ! IN
    dsl4jsb_Get_var2D_onChunk(HYDRO,    fract_snow_can)        ! IN
    dsl4jsb_Get_var2D_onChunk(HYDRO,    fract_snow_soil)       ! IN

    !RMS dsl4jsb_Get_var2D_onChunk(CARBON,   Cpool_litter_green_ag) ! IN
    !RMS dsl4jsb_Get_var2D_onChunk(CARBON,   Cpool_slow)            ! IN
    dsl4jsb_Get_var2D_onChunk(RAD,      swvis_srf_net)         ! IN
    dsl4jsb_Get_var2D_onChunk(RAD,      swnir_srf_net)         ! IN
    dsl4jsb_Get_var2D_onChunk(RAD,      snow_age)              ! IN
    dsl4jsb_Get_var2D_onChunk(RAD,      alb_background)        ! IN

    dsl4jsb_Get_var2D_onChunk(RAD,      alb_vis)               ! INOUT
    dsl4jsb_Get_var2D_onChunk(RAD,      alb_nir)               ! INOUT
    dsl4jsb_Get_var2D_onChunk(RAD,      alb_vis_snow)          ! INOUT
    dsl4jsb_Get_var2D_onChunk(RAD,      alb_nir_snow)          ! INOUT
    dsl4jsb_Get_var2D_onChunk(RAD,      alb_vis_can)           ! INOUT
    dsl4jsb_Get_var2D_onChunk(RAD,      alb_nir_can)           ! INOUT
    dsl4jsb_Get_var2D_onChunk(RAD,      alb_vis_soil)          ! INOUT
    dsl4jsb_Get_var2D_onChunk(RAD,      alb_nir_soil)          ! INOUT
    dsl4jsb_Get_var2D_onChunk(RAD,      alb_vis_mineralsoil)   ! INOUT
    dsl4jsb_Get_var2D_onChunk(RAD,      alb_nir_mineralsoil)   ! INOUT
    dsl4jsb_Get_var2D_onChunk(RAD,      sky_view_fract)        ! INOUT
    dsl4jsb_Get_var2D_onChunk(RAD,      sky_view_fract_stem)   ! INOUT

       ! if tile is lake:

       ! if tile is glacier:
!        IF (use_alb_veg_simple) THEN
!           alb_vis(:) =albedo_snow_grnd aus Prozedur get_surface_albedo_simple
!           alb_nir(:) =albedo_snow_grnd aus Prozedur get_surface_albedo_simple
!           ! R: Die Prozedur get_surface_albedo_simple mit Glacier-maske auszustatten wäre simpel aber ich finde das unschön.
!                Die Frage ist ob die Prozedur wirklich benötigt wird und ob für Glacier dasselbe rauskommen muss...
!        ELSE
!            CALL calc_glacier_albedo(t(:),             & ! In
!                                     swvis_srf_net(:), & ! In
!                                     swnir_srf_net(:), & ! In
!                                     alb_vis(:),       & ! Out
!                                     alb_nir(:) )        ! Out
!        END IF
!
!
!         ! if tile is bare soil:
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

        ! if tile is vegetation:
 !!!       IF (use_alb_veg_simple) THEN

           ! if tile is lowest child tile
              ! Einfache Berechnung der Albedo eines Vegetation-Tiles mit der Funktion canopy_albedo_echam5:
              alb_vis(:) = get_surface_albedo_simple(t(:), fract_forest(:), fract_snow_soil(:), fract_snow_can(:), &
                                                         lai(:), alb_background(:))
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
! !                                     alb_nir_soil(:) )            ! InOut
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

           CALL calc_snow_albedo(fract_snow_soil(:),  & ! In
                                 snow_age(:),         & ! In
                                 cos_zenith_angle(:), & ! In
                                 albedo_age_weight,   & ! In
                                 t(:),                & ! In
                                 alb_vis_snow(:),     & ! Out
                                 alb_nir_snow(:) )      ! Out

! !
! !            CALL calc_sky_view_fractions(lai(:),                & ! In
! !                                         cos_zenith_angle(:),   & ! In
! !                                         sky_view_fract(:),     & ! Out
! !                                         sky_view_fract_stem(:) ) ! Out
! !
! !            ! R: All inputs are constants, not necessary to calculate this for each grid cell again. Therefore here and not in merge_albedos.               
! !            ! Weight between age and temperature impact on snow albedo for snow on canopy
! !            AlbedoCanopySnow = albedo_age_weight * AlbedoCanopySnow_age + (1._wp - albedo_age_weight) * AlbedoCanopySnow_temp
! !
! !            CALL calc_merge_albedos(swvis_srf_net(:),       & ! In
! !                                    swnir_srf_net(:),       & ! In
! !                                    sky_view_fract(:),      & ! In
! !                                    sky_view_fract_stem(:), & ! In
! !                                    AlbedoCanopySnow,       & ! In
! !                                    fract_snow_soil(:),     & ! In
! !                                    fract_snow_can(:),      & ! In
! !                                    alb_vis_snow(:),        & ! In
! !                                    alb_nir_snow(:),        & ! In
! !                                    alb_vis_soil(:),        & ! In
! !                                    alb_nir_soil(:),        & ! In
! !                                    alb_vis_can(:),         & ! In
! !                                    alb_nir_can(:),         & ! In
! !                                    alb_vis_srf(:),         & ! Out
! !                                    alb_nir_srf(:) )          ! Out
! !
! !         END IF ! from use_alb_veg_simple

  END SUBROUTINE update_radiation_albedo

#endif
END MODULE mo_rad_interface
