!>
!! @brief Contains the interfaces to the assimi process
!!
!! @author
!!  Rainer Schneck, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Rainer Schneck (2015-05-30)
!!
MODULE mo_assimi_interface
#ifndef __NO_JSBACH__

#include "dsl_jsbach.inc"
  ! -------------------------------------------------------------------------------------------------------
  ! Used variables of module

  ! Use of basic structures
  USE mo_jsb_control,     ONLY: debug_on
  USE mo_kind,            ONLY: wp
  USE mo_exception,       ONLY: message

  USE mo_jsb_grid_class,     ONLY: t_jsb_grid
  USE mo_jsb_model_class,    ONLY: t_jsb_model
  USE mo_jsb_class,          ONLY: Get_model
  USE mo_jsb_tile_class,     ONLY: t_jsb_tile_abstract
  USE mo_jsb_config_class,   ONLY: t_jsb_config, t_jsb_config_p
  USE mo_jsb_process_class,  ONLY: t_jsb_process
  USE mo_jsb_task_class,     ONLY: t_jsb_process_task, t_jsb_task_options

  ! Use of processes in this module (Get_assimi_memory and Get_assimi_config)
  dsl4jsb_Use_processes ASSIMI, VEG, RAD, ATM2LAND

  ! Use of process configurations (t_assimi_config)
  dsl4jsb_Use_config(ASSIMI)

  ! Use of process memories (t_assimi_memory)
  dsl4jsb_Use_memory(ATM2LAND)
  dsl4jsb_Use_memory(ASSIMI)
  dsl4jsb_Use_memory(VEG)
  dsl4jsb_Use_memory(RAD)

  ! -------------------------------------------------------------------------------------------------------
  ! Module variables

  IMPLICIT NONE
  PRIVATE
  PUBLIC ::  Register_assimi_tasks !,t_assimi_process

  CHARACTER(len=*), PARAMETER :: modname = 'mo_assimi_interface'

  !> Type definition for assimi process
!   TYPE, EXTENDS(t_jsb_process) :: t_assimi_process
!   CONTAINS
!     PROCEDURE :: Register_tasks => Register_assimi_tasks        !< Registers tasks for assimi process
!   END TYPE t_assimi_process


  !> Type definition for assimilation_scaling_factors task
  TYPE, EXTENDS(t_jsb_process_task) ::   tsk_assimilation_scaling_factors
  CONTAINS
    PROCEDURE, NOPASS :: Integrate => update_assimilation_scaling_factors    !< Advances task computation for one timestep
    PROCEDURE, NOPASS :: Aggregate => aggregate_assimi !aggregate_assimilation_scaling_factors !< Aggregates computed task variables
  END TYPE   tsk_assimilation_scaling_factors

  !> Constructor interface for assimilation_scaling_factors task
  INTERFACE   tsk_assimilation_scaling_factors
    PROCEDURE Create_task_assimilation_scaling_factors         !< Constructor function for task
  END INTERFACE   tsk_assimilation_scaling_factors


  !> Type definition for assimilation water unlimited task
  TYPE, EXTENDS(t_jsb_process_task) :: tsk_assimilation_waterunlimited
  CONTAINS
    PROCEDURE, NOPASS :: Integrate => update_assimilation_waterunlimited    !< Advances task computation for one timestep
    PROCEDURE, NOPASS :: Aggregate => aggregate_assimi !aggregate_assimilation_waterunlimited !< Aggregates computed task variables
  END TYPE tsk_assimilation_waterunlimited

  !> Constructor interface for assimi task
  INTERFACE tsk_assimilation_waterunlimited
    PROCEDURE Create_task_assimilation_waterunlimited     !< Constructor function for task
  END INTERFACE tsk_assimilation_waterunlimited



  !> Type definition for assimilation water limited task
  TYPE, EXTENDS(t_jsb_process_task) :: tsk_assimilation_waterlimited
  CONTAINS
    PROCEDURE, NOPASS :: Integrate => update_assimilation_waterlimited    !< Advances task computation for one timestep
    PROCEDURE, NOPASS :: Aggregate => aggregate_assimi !aggregate_assimilation_waterlimited !< Aggregates computed task variables
  END TYPE tsk_assimilation_waterlimited

  !> Constructor interface for assimi task
  INTERFACE tsk_assimilation_waterlimited
    PROCEDURE Create_task_assimilation_waterlimited     !< Constructor function for task
  END INTERFACE tsk_assimilation_waterlimited



CONTAINS

  ! ================================================================================================================================
  !! Constructors for tasks

  ! -------------------------------------------------------------------------------------------------------
  !> Constructor for assimi task
  !!
  !! @param[in]     model_id     Model id
  !! @return        return_ptr   Instance of process task "assimi"
  !!
  FUNCTION Create_task_assimilation_scaling_factors(model_id) RESULT(return_ptr)

    INTEGER,                   INTENT(in) :: model_id
    CLASS(t_jsb_process_task), POINTER    :: return_ptr

    ALLOCATE(tsk_assimilation_scaling_factors::return_ptr)

    return_ptr%name = 'assimilation_scaling_factors'
    return_ptr%process_id = ASSIMI
    return_ptr%owner_model_id = model_id

  END FUNCTION Create_task_assimilation_scaling_factors

  ! -------------------------------------------------------------------------------------------------------
  !> Constructor for assimi task
  !!
  !! @param[in]     model_id     Model id
  !! @return        return_ptr   Instance of process task "assimi"
  !!
  FUNCTION Create_task_assimilation_waterunlimited(model_id) RESULT(return_ptr)

    INTEGER,                   INTENT(in) :: model_id
    CLASS(t_jsb_process_task), POINTER    :: return_ptr

    ALLOCATE(tsk_assimilation_waterunlimited::return_ptr)

    return_ptr%name = 'assimilation_waterunlimited'
    return_ptr%process_id = ASSIMI
    return_ptr%owner_model_id = model_id

  END FUNCTION Create_task_assimilation_waterunlimited

  ! -------------------------------------------------------------------------------------------------------
  !> Constructor for assimi task
  !!
  !! @param[in]     model_id     Model id
  !! @return        return_ptr   Instance of process task "assimi"
  !!
  FUNCTION Create_task_assimilation_waterlimited(model_id) RESULT(return_ptr)

    INTEGER,                   INTENT(in) :: model_id
    CLASS(t_jsb_process_task), POINTER    :: return_ptr

    ALLOCATE(tsk_assimilation_waterlimited::return_ptr)

    return_ptr%name = 'assimilation_waterlimited'
    return_ptr%process_id = ASSIMI
    return_ptr%owner_model_id = model_id

  END FUNCTION Create_task_assimilation_waterlimited

  ! -------------------------------------------------------------------------------------------------------
  !> Register tasks for assimi process
  !!
  !! @param[in,out] this      Instance of assimi process class
  !! @param[in]     model_id  Model id
  !!
  SUBROUTINE Register_assimi_tasks(this, model_id)

    CLASS(t_jsb_process), INTENT(inout) :: this
    INTEGER,                 INTENT(in) :: model_id

    CALL this%Register_task(tsk_assimilation_scaling_factors(model_id))
    CALL this%Register_task(tsk_assimilation_waterlimited(model_id))
    CALL this%Register_task(tsk_assimilation_waterunlimited(model_id))

  END SUBROUTINE Register_assimi_tasks

  ! ================================================================================================================================
  !>
  !! @brief Implementation to calculate the height scaling factors of plant enzyms for assimilation
  !! Task "assimi" calculates the carbon assimilation of plants.
  !!
  !! @param[in,out] tile    Tile for which routine is executed.
  !! @param[in]     options Additional run-time parameters.
  !!
  ! The scaling factors represent the enzymatic distribution of plants in their canopy layers. This correltates
  ! with the amount of nitrogen in these layers (therefore called nitrogen_scaling_factors in JSBACH3).
  ! Eqs. (106), (107) in Knorr (1997)
  ! Note that this routine should only be applied to special vegetation classes
  SUBROUTINE update_assimilation_scaling_factors(tile, options)

    ! Declarations
    !USE mo_orbit,                 ONLY: inquire_declination
     ! R: Returns declination calculated in last orbit call.
     !    Problem: FATAL ERROR in inquire_declination: Not Initialized
     !    This is because ICON did not call the calculation of orbital parameters in SUBROUTINE orbit_vsop87 or
     !    alternatively in SUBROUTINE orbit_kepler. Therefore this seems to be a problem solved with the development of ICON...
     !    Note: problem still exists in the 26.10.2015 version of icon-les-dev !
     !    However taken from mo_psrad_orbit works... perhaps precompiler command necessary here...
    USE mo_psrad_orbit,       ONLY: inquire_declination
    USE mo_jsb_grid,          ONLY: Get_grid
    USE mo_assimi_process,    ONLY: calc_assimilation_scaling_factors

    ! Arguments
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    ! Local variables
    TYPE(t_jsb_model), POINTER   :: model
    TYPE(t_jsb_grid),  POINTER   :: grid
    INTEGER                      :: iblk, ics, ice, nc
    REAL(wp)                     :: lat(options%nc)
    REAL(wp)                     :: declination ! Solar Declination angle
    INTEGER                      :: i
    INTEGER                      :: ncanopy          ! number of canopy layers
    REAL(wp),          POINTER   :: canopy_bound_lai(:)

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':assimilation_scaling_factors'

    ! Declare process configuration and memory Pointers
    dsl4jsb_Def_config(ASSIMI)
    dsl4jsb_Def_memory(ASSIMI)
    dsl4jsb_Def_memory(VEG)

    ! Declare pointers to variables in memory
    dsl4jsb_Real3D_onChunk :: assimi_scaling_fact_cl
    dsl4jsb_Real2D_onChunk :: lai

    ! Get local variables from options argument
    iblk    = options%iblk
    ics     = options%ics
    ice     = options%ice
    nc      = options%nc

    ! If process is not active on this tile, do nothing
    IF (.NOT. tile%Is_process_active(ASSIMI)) RETURN

    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    model => Get_model(tile%owner_model_id)

    ! Get process config
    dsl4jsb_Get_config(ASSIMI)

    ! Get process memories
    dsl4jsb_Get_memory(ASSIMI)
    dsl4jsb_Get_memory(VEG)

    ! Get process variables
    ncanopy          =  dsl4jsb_Config(ASSIMI)%ncanopy
    canopy_bound_lai(0:ncanopy) => dsl4jsb_Config(ASSIMI)%canopy_bound_lai(0:ncanopy)

    grid => Get_grid(model%grid_id)
    lat(:)  = grid%lat(ics:ice, iblk)                           ! in

    dsl4jsb_Get_var2D_onChunk(VEG, lai)                       ! in
    dsl4jsb_Get_var3D_onChunk(ASSIMI, assimi_scaling_fact_cl) ! out

    ! ---------------------------
    ! Go

    CALL inquire_declination(declination) ! R: do not put this in elemental procedure as it would not be pure...
                                          ! e.g. declination = 0.006918_wp

    IF (dsl4jsb_Lctlib_param(NitrogenScalingFlag)) THEN
      DO i=1,ncanopy
        CALL calc_assimilation_scaling_factors(declination,                 & ! input
                                               canopy_bound_lai(i-1),       & ! input
                                               lai(:),                      & ! input
                                               lat(:),                      & ! input
                                               assimi_scaling_fact_cl(:,i)  & ! output
                                               )
      END DO
    ELSE                                           ! to some vegetation classes nitrogen scaling is not applied
      assimi_scaling_fact_cl(:,:) = 1._wp
    END IF

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE update_assimilation_scaling_factors

  ! ================================================================================================================================
  !>
  !! @brief Implementation of assimilation for water unlimit conditions
  !! Task "assimi" calculates the carbon assimilation of plants.
  !!
  !! @param[in,out] tile    Tile for which routine is executed.
  !! @param[in]     options Additional run-time parameters.
  !!
  SUBROUTINE update_assimilation_waterunlimited(tile, options)

    USE mo_assimi_process,      ONLY: calc_assimilation_waterunlimited, calc_NPP_rate
    USE mo_bc_greenhouse_gases, ONLY: ghg_co2mmr ! CO2 mass mixing ratio of dry air in kg/kg. Could be named herer CO2_kg_mixing_ratio.
                                                 ! This is preliminary.
                                                 ! ghg_co2mmr [kg/kg] = zco2int [m³/m³]  * amco2 [g/mol] / amd [g/mol]
                                                 ! mass mixing ratio  = vol mixing ratio * MM CO2        / MM dry air
    USE mo_carbon_constants,    ONLY: molarMassDryAir, molarMassCO2

    ! Arguments
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    ! Local variables
    TYPE(t_jsb_model), POINTER                :: model
    REAL(wp) :: dtime
    INTEGER  :: iblk, ics, ice, nc
    INTEGER                             :: icanopy

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':update_assimilation_waterunlimited'

    ! Declare process memories
    dsl4jsb_Def_memory(ASSIMI)
    dsl4jsb_Def_memory(VEG)
    dsl4jsb_Def_memory(RAD)
    dsl4jsb_Def_memory(ATM2LAND)

    dsl4jsb_Def_config(ASSIMI)

    ! Declare pointers to variables in memory
    INTEGER, POINTER        :: ncanopy                 ! number of canopy layers

    !dsl4jsb_Real2D_onChunk :: C4flag ! Photosynthetic pathway (C3: 0; C4: 1)
                                      ! R: nur falls man C4flag doch als REAL haben wollte

    dsl4jsb_Real2D_onChunk :: t_air                 ! Atmosphere temperature (lowest layer) in Kelvin!
    dsl4jsb_Real2D_onChunk :: press_srf             ! Surface pressure
    dsl4jsb_Real2D_onChunk :: par_down_mol          ! Downward PAR flux in mol (photons)/(m^2 s)
    dsl4jsb_Real3D_onChunk :: apar_per_lai_cl       ! Input. Absorbed PAR of canopy layer.
    dsl4jsb_Real2D_onChunk :: canopy_cond           ! Output. Canopy conductance no limitation
    dsl4jsb_Real2D_onChunk :: canopy_cond_limited   ! Canopy conductance water limited

   ! R: As soon as I get CO2 from the atmosphere
    ! dsl4jsb_Real2D_onChunk :: CO2_kg_mixing_ratio ! CO2 mass mixing ratio [kgCO2/kgDryAir]
    dsl4jsb_Real2D_onChunk :: CO2_mol_mixing_ratio  ! Output. CO2 particle mixing ratio [molCO2/molDryAir]

    dsl4jsb_Real2D_onChunk :: gross_assimilation
    dsl4jsb_Real2D_onChunk :: dark_respiration
    dsl4jsb_Real2D_onChunk :: NPP_rate
    dsl4jsb_Real3D_onChunk :: assimi_scaling_fact_cl ! Input.
    dsl4jsb_Real3D_onChunk :: canopy_cond_cl         ! Output (for waterunlimited case)
    dsl4jsb_Real3D_onChunk :: gross_assimilation_cl  ! Output.
    dsl4jsb_Real3D_onChunk :: dark_respiration_cl    ! Output.
    dsl4jsb_Real3D_onChunk :: carbox_rate_max_cl     ! Output. Maximum carboxylation rate of canopy layer (=VCmax)
                                                         ! [Mol(CO2)/m^2(canopy)/s]. carbox_rate_max_leaf adapted to canopy layer
                                                         ! and actual temperature.
    dsl4jsb_Real3D_onChunk :: e_transport_rate_max_cl ! Output.
    dsl4jsb_Real3D_onChunk :: carbox_rate_cl          ! Output.
    dsl4jsb_Real3D_onChunk :: e_transport_rate_cl     ! Output.
    dsl4jsb_Real3D_onChunk :: CO2_conc_leaf_cl        ! CO2 concentration inside leaf  [mol(CO2)/mol(Air)]:
                                                      ! Input:  for waterLimitationFlag=.false..(Input identical for all layers)
                                                      ! Output: for waterLimitationFlag=.true.
    dsl4jsb_Real3D_onChunk :: lai_cl                  ! Leaf area index [-]


    ! Get local variables from options argument
    dtime   = options%dtime
    iblk    = options%iblk
    ics     = options%ics
    ice     = options%ice
    nc      = options%nc

    ! If process is not active on this tile, do nothing
    IF (.NOT. tile%Is_process_active(ASSIMI)) RETURN

    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    model => Get_model(tile%owner_model_id)

    dsl4jsb_Get_config(ASSIMI)

    ! Get process memories
    dsl4jsb_Get_memory(ASSIMI)
    dsl4jsb_Get_memory(VEG)
    dsl4jsb_Get_memory(RAD)
    dsl4jsb_Get_memory(ATM2LAND)

    ! Get process variables
    dsl4jsb_Get_var2D_onChunk(ATM2LAND, t_air)                   ! IN
    dsl4jsb_Get_var2D_onChunk(ATM2LAND, press_srf)               ! IN
    dsl4jsb_Get_var2D_onChunk(RAD,      par_down_mol)            ! IN
    dsl4jsb_Get_var3D_onChunk(RAD,      apar_per_lai_cl)         ! IN
    dsl4jsb_Get_var3D_onChunk(RAD,      lai_cl)                  ! IN
    dsl4jsb_Get_var2D_onChunk(VEG,      canopy_cond)             ! OUT

    dsl4jsb_Get_var2D_onChunk(ASSIMI,   CO2_mol_mixing_ratio)    ! OUT

    dsl4jsb_Get_var2D_onChunk(ASSIMI,   gross_assimilation)      ! OUT
    dsl4jsb_Get_var2D_onChunk(ASSIMI,   dark_respiration)        ! OUT
    dsl4jsb_Get_var2D_onChunk(ASSIMI,   NPP_rate)                ! OUT

    dsl4jsb_Get_var3D_onChunk(ASSIMI,   assimi_scaling_fact_cl)  ! OUT
    dsl4jsb_Get_var3D_onChunk(ASSIMI,   canopy_cond_cl)          ! OUT without water limitation.
                                                                 ! R: canopy_cond_cl will be overwritten within this procedure,
                                                                 !    nomatter what value comes in with it
    !dsl4jsb_Get_var3D_onChunk(ASSIMI,   gross_assimilation_cl)   ! INOUT
    dsl4jsb_Get_var3D_onChunk(ASSIMI,   dark_respiration_cl)     ! INOUT
    !dsl4jsb_Get_var3D_onChunk(ASSIMI,   carbox_rate_max_cl)      ! INOUT
    !dsl4jsb_Get_var3D_onChunk(ASSIMI,   e_transport_rate_max_cl) ! INOUT
    !dsl4jsb_Get_var3D_onChunk(ASSIMI,   e_transport_rate_cl)     ! INOUT
    dsl4jsb_Get_var3D_onChunk(ASSIMI,   carbox_rate_cl)          ! INOUT
    dsl4jsb_Get_var3D_onChunk(ASSIMI,   CO2_conc_leaf_cl)        ! INOUT for no water limitation case here


    ncanopy                  => dsl4jsb_Config(ASSIMI)%ncanopy

    ! Convert CO2 mass mixing ratio [kg/kg] to particle mixing ratio [mol/mol]
    ! Note: the conversion to CO2_mol_mixing_ratio is only used here therefore a separate function is unnecessary.

    CO2_mol_mixing_ratio = ghg_co2mmr * molarMassDryAir / molarMassCO2 ! e.g. CO2_mol_mixing_ratio = 0.280


    ! ---------------------------
    ! Go

    DO icanopy=1,ncanopy
       CALL calc_assimilation_waterunlimited(                      &
                               !waterLimitationFlag,                & ! in
                               dsl4jsb_Lctlib,                     & ! in
                               t_air(:),                           & ! in
                               press_srf(:),                       & ! in
                               par_down_mol(:),                    & ! in
                               apar_per_lai_cl(:,icanopy),         & ! in
                               CO2_mol_mixing_ratio(:),            & ! in
                               assimi_scaling_fact_cl(:,icanopy),  & ! in
                               CO2_conc_leaf_cl(:,icanopy),        & ! inout for no water limitation case here
                               canopy_cond_cl(:,icanopy),          & ! out   for no water limitation case here
                               !gross_assimilation_cl(:,icanopy),   & ! out
                               dark_respiration_cl(:,icanopy),     & ! out
                               !carbox_rate_max_cl(:,icanopy),      & ! out
                               !e_transport_rate_max_cl(:,icanopy), & ! out
                               carbox_rate_cl(:,icanopy)          & ! out
                               !e_transport_rate_cl(:,icanopy)      & ! out
                               )
    END DO


    ! Compute unlimited conductance for the whole canopy (all layers) considering the lai
    canopy_cond(:) = MAX(1.e-20_wp,SUM(canopy_cond_cl(:,:) * lai_cl(:,:),DIM=2))

    ! R: raus genommen, da wahrscheinlich nicht mehr von Interesse:
    ! ! Integrate over time
    ! canopy_cond_over_time_unlimited(:) = canopy_cond_over_time_unlimited(:) + canopy_cond(:) * dtime

    ! Same for gross_assimilation
    ! R: dies ist eine alternative Schreibweise die verwendet wurde da der PGI Kompiler probleme gemacht hat.
    !gross_assimilation(:) = 0._wp
    !DO icanopy=1,ncanopy
    !   gross_assimilation(:) = gross_assimilation(:) + gross_assimilation_cl(:,icanopy) * lai_cl(:,icanopy)
    !END DO
    ! R: habe ich raus genommen auch wenn das JSBACH3 entspricht. Macht keinen Sinn da gross_assimilation_cl oben nicht neu berechnet wird.
    !gross_assimilation(:) = SUM(gross_assimilation_cl(:,:) * lai_cl(:,:),DIM=2)

    ! Same for dark_respiration
    ! R: dies ist eine alternative Schreibweise die verwendet wurde da der PGI Kompiler probleme gemacht hat.
    !dark_respiration(:) = 0._wp
    !DO icanopy=1,ncanopy
    !   dark_respiration(:) =  dark_respiration(:) + dark_respiration_cl(:,icanopy) * lai_cl(:,icanopy)
    !END DO

    dark_respiration(:) = SUM(dark_respiration_cl(:,:) * lai_cl(:,:),DIM=2)

!    R: siehe mo_jsbach_interface.f90, Zeile 671
!    R: this call should be after the first call of update_assimilation and before the first
!       call of update_soil.
!       Diese Übergabe ist in JSBACH4 nicht mehr notwendig, da Variable canopy_cond nicht mehr lokal ist
!       canopy_conductance_unstressed(1:nidx,1:ntiles)=theLand%Vegetation%canopy_conductance_bethy(kidx0:kidx1,1:ntiles)


    CALL calc_NPP_rate(gross_assimilation, dark_respiration, NPP_rate)

    RETURN



    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE update_assimilation_waterunlimited

  ! ================================================================================================================================
  !>
  !! @brief Implementation of assimilation for water unlimit conditions
  !! Task "assimi" calculates the carbon assimilation of plants.
  !!
  !! @param[in,out] tile    Tile for which routine is executed.
  !! @param[in]     options Additional run-time parameters.
  !!
  SUBROUTINE update_assimilation_waterlimited(tile, options)

    USE mo_assimi_process,      ONLY: calc_assimilation_waterlimited, calc_NPP_rate
    USE mo_bc_greenhouse_gases, ONLY: ghg_co2mmr ! CO2 mass mixing ratio of dry air in kg/kg. Could be named herer CO2_kg_mixing_ratio.
                                                 ! This is preliminary.
                                                 ! ghg_co2mmr [kg/kg] = zco2int [m³/m³]  * amco2 [g/mol] / amd [g/mol]
                                                 ! mass mixing ratio  = vol mixing ratio * MM CO2        / MM dry air
    USE mo_carbon_constants,    ONLY: molarMassDryAir, molarMassCO2

    ! Arguments
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    ! Local variables
    TYPE(t_jsb_model), POINTER                :: model
    REAL(wp) :: dtime
    INTEGER  :: iblk, ics, ice, nc
    INTEGER                             :: icanopy

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':update_assimilation_waterunlimited'

    ! Declare process memories
    dsl4jsb_Def_memory(ASSIMI)
    dsl4jsb_Def_memory(VEG)
    dsl4jsb_Def_memory(RAD)
    dsl4jsb_Def_memory(ATM2LAND)

    dsl4jsb_Def_config(ASSIMI)

    ! Declare pointers to variables in memory
    INTEGER, POINTER        :: ncanopy                 ! number of canopy layers

    !dsl4jsb_Real2D_onChunk :: C4flag ! Photosynthetic pathway (C3: 0; C4: 1)
                                      ! R: nur falls man C4flag doch als REAL haben wollte

    dsl4jsb_Real2D_onChunk :: t_air                 ! Atmosphere temperature (lowest layer) in Kelvin!
    dsl4jsb_Real2D_onChunk :: press_srf             ! Surface pressure
    dsl4jsb_Real2D_onChunk :: par_down_mol          ! Downward PAR flux in mol (photons)/(m^2 s)
    dsl4jsb_Real3D_onChunk :: apar_per_lai_cl       ! Input. Absorbed PAR of canopy layer.
    dsl4jsb_Real2D_onChunk :: canopy_cond           ! Canopy conductance no limitation
    dsl4jsb_Real2D_onChunk :: canopy_cond_limited   ! Canopy conductance water limited
    dsl4jsb_Real2D_onChunk :: canopy_cond_over_time ! Canopy conductance integrated over time

   ! R: As soon as I get CO2 from the atmosphere
    ! dsl4jsb_Real2D_onChunk :: CO2_kg_mixing_ratio ! CO2 mass mixing ratio [kgCO2/kgDryAir]
    dsl4jsb_Real2D_onChunk :: CO2_mol_mixing_ratio  ! Output. CO2 particle mixing ratio [molCO2/molDryAir]

    dsl4jsb_Real2D_onChunk :: gross_assimilation
    dsl4jsb_Real2D_onChunk :: dark_respiration
    dsl4jsb_Real2D_onChunk :: NPP_rate
    dsl4jsb_Real3D_onChunk :: assimi_scaling_fact_cl ! Input.
    dsl4jsb_Real3D_onChunk :: canopy_cond_cl         ! InOut for water limited case here.
    dsl4jsb_Real3D_onChunk :: gross_assimilation_cl  ! Output.
    dsl4jsb_Real3D_onChunk :: dark_respiration_cl    ! Output.
    dsl4jsb_Real3D_onChunk :: carbox_rate_max_cl     ! Output. Maximum carboxylation rate of canopy layer (=VCmax)
                                                         ! [Mol(CO2)/m^2(canopy)/s]. carbox_rate_max_leaf adapted to canopy layer
                                                         ! and actual temperature.
    dsl4jsb_Real3D_onChunk :: e_transport_rate_max_cl ! Output.
    dsl4jsb_Real3D_onChunk :: carbox_rate_cl          ! Output.
    dsl4jsb_Real3D_onChunk :: e_transport_rate_cl     ! Output.
    dsl4jsb_Real3D_onChunk :: CO2_conc_leaf_cl        ! CO2 concentration inside leaf  [mol(CO2)/mol(Air)]:
                                                      ! Output only for the water limited case here
    dsl4jsb_Real3D_onChunk :: lai_cl                  ! Leaf area index [-]

    ! Get local variables from options argument
    dtime   = options%dtime
    iblk    = options%iblk
    ics     = options%ics
    ice     = options%ice
    nc      = options%nc

    ! If process is not active on this tile, do nothing
    IF (.NOT. tile%Is_process_active(ASSIMI)) RETURN

    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    model => Get_model(tile%owner_model_id)

    dsl4jsb_Get_config(ASSIMI)

    ! Get process memories
    dsl4jsb_Get_memory(ASSIMI)
    dsl4jsb_Get_memory(VEG)
    dsl4jsb_Get_memory(RAD)
    dsl4jsb_Get_memory(ATM2LAND)

    ! Get process variables
    dsl4jsb_Get_var2D_onChunk(ATM2LAND, t_air)                   ! IN
    dsl4jsb_Get_var2D_onChunk(ATM2LAND, press_srf)               ! IN
    dsl4jsb_Get_var2D_onChunk(RAD,      par_down_mol)            ! IN
    dsl4jsb_Get_var3D_onChunk(RAD,      apar_per_lai_cl)         ! IN
    dsl4jsb_Get_var3D_onChunk(RAD,      lai_cl)                  ! IN
    dsl4jsb_Get_var2D_onChunk(VEG,      canopy_cond)             ! OUT
    dsl4jsb_Get_var2D_onChunk(VEG,      canopy_cond_limited)     ! IN , only here with water limitation

    dsl4jsb_Get_var2D_onChunk(ASSIMI,   CO2_mol_mixing_ratio)    ! OUT

    dsl4jsb_Get_var2D_onChunk(ASSIMI,   gross_assimilation)      ! OUT
    dsl4jsb_Get_var2D_onChunk(ASSIMI,   dark_respiration)        ! OUT
    dsl4jsb_Get_var2D_onChunk(ASSIMI,   NPP_rate)                ! OUT
    dsl4jsb_Get_var2D_onChunk(ASSIMI,   canopy_cond_over_time)   ! INOUT

    dsl4jsb_Get_var3D_onChunk(ASSIMI,   assimi_scaling_fact_cl)  ! OUT
    dsl4jsb_Get_var3D_onChunk(ASSIMI,   canopy_cond_cl)          ! INOUT for water limitation case here
                                                                 ! R: canopy_cond_cl is here in the second call of update_assimilation
                                                                 !    for the water limitation case an input.
    dsl4jsb_Get_var3D_onChunk(ASSIMI,   gross_assimilation_cl)   ! INOUT
    dsl4jsb_Get_var3D_onChunk(ASSIMI,   dark_respiration_cl)     ! INOUT
    dsl4jsb_Get_var3D_onChunk(ASSIMI,   carbox_rate_max_cl)      ! INOUT
    dsl4jsb_Get_var3D_onChunk(ASSIMI,   e_transport_rate_max_cl) ! INOUT
    dsl4jsb_Get_var3D_onChunk(ASSIMI,   e_transport_rate_cl)     ! INOUT
    dsl4jsb_Get_var3D_onChunk(ASSIMI,   carbox_rate_cl)          ! INOUT
    dsl4jsb_Get_var3D_onChunk(ASSIMI,   CO2_conc_leaf_cl)        ! OUT


    ncanopy                  => dsl4jsb_Config(ASSIMI)%ncanopy          ! IN

    ! Convert CO2 mass mixing ratio [kg/kg] to particle mixing ratio [mol/mol]
    ! Note: the conversion to CO2_mol_mixing_ratio is only used here therefore a separate function is unnecessary.

    CO2_mol_mixing_ratio = ghg_co2mmr * molarMassDryAir / molarMassCO2 ! e.g. CO2_mol_mixing_ratio = 0.280


    ! ---------------------------
    ! Go

    ! R: Diese Berechung hier habe ich raus genommen, da CO2_conc_leaf_unlimited_acc weder für assimilation noch wo anderst verwendet wird.
    !    Wird nur ausgegeben. Braucht man wahrscheinlich nicht mehr, auch nicht im Output.
    ! bethy%CO2_conc_leaf_unlimited_acc(kidx0:kidx1,itile) = bethy%CO2_conc_leaf_unlimited_acc(kidx0:kidx1,itile) + &
    !                                                             CO2_conc_leaf_cl * dtime


    ! Compute unlimited conductance for the whole canopy (integrate over canopy layers)
    canopy_cond(:) = MAX(1.e-20_wp,SUM(canopy_cond_cl(:,:) * lai_cl(:,:),DIM=2))

    ! Scale unlimited conductance per leaf area with the ratio of limited to unlimited conductance for
    ! whole canopy, ie. compute new limited conductance per leaf area and canopy layer
    DO icanopy=1,ncanopy
       canopy_cond_cl(:,icanopy) = canopy_cond_cl(:,icanopy) * canopy_cond_limited(:) / canopy_cond(:) ! canopy_cond_limited comes from soil
    END DO

    DO icanopy=1,ncanopy
       CALL calc_assimilation_waterlimited(                      &
                             !waterLimitationFlag,                & ! in
                             dsl4jsb_Lctlib,                     & ! in
                             t_air(:),                           & ! in
                             press_srf(:),                       & ! in
                             par_down_mol(:),                    & ! in
                             apar_per_lai_cl(:,icanopy),         & ! in
                             CO2_mol_mixing_ratio(:),            & ! in
                             assimi_scaling_fact_cl(:,icanopy),  & ! in
                             CO2_conc_leaf_cl(:,icanopy),        & ! out only for the water limited case here
                             canopy_cond_cl(:,icanopy),          & ! in  only for the water limited case here
                             gross_assimilation_cl(:,icanopy),   & ! out
                             dark_respiration_cl(:,icanopy),     & ! out
                             carbox_rate_max_cl(:,icanopy),      & ! out
                             e_transport_rate_max_cl(:,icanopy), & ! out
                             carbox_rate_cl(:,icanopy),          & ! out
                             e_transport_rate_cl(:,icanopy)      & ! out
                             )
    END DO


    ! R: siehe mo_bethy.f90 Zeile 879
    ! Compute unlimited conductance for the whole canopy (all layers) considering the lai
    canopy_cond(:) = MAX(1.e-20_wp,SUM(canopy_cond_cl(:,:) * lai_cl(:,:),DIM=2))

    ! Integrate over time
    canopy_cond_over_time(:) = canopy_cond_over_time(:) + canopy_cond(:) * dtime

    ! Same for gross_assimilation
    ! R: dies ist eine alternative Schreibweise die verwendet wurde da der PGI Kompiler probleme gemacht hat.
    !gross_assimilation(:) = 0._wp
    !DO icanopy=1,ncanopy
    !   gross_assimilation(:) = gross_assimilation(:) + gross_assimilation_cl(:,icanopy) * lai_cl(:,icanopy)
    !END DO

    gross_assimilation(:) = SUM(gross_assimilation_cl(:,:) * lai_cl(:,:),DIM=2)

    ! Same for dark_respiration
    ! R: dies ist eine alternative Schreibweise die verwendet wurde da der PGI Kompiler probleme gemacht hat.
    !dark_respiration(:) = 0._wp
    !DO icanopy=1,ncanopy
    !   dark_respiration(:) =  dark_respiration(:) + dark_respiration_cl(:,icanopy) * lai_cl(:,icanopy)
    !END DO

    dark_respiration(:) = SUM(dark_respiration_cl(:,:) * lai_cl(:,:),DIM=2)


    CALL calc_NPP_rate(gross_assimilation, dark_respiration, NPP_rate)

    RETURN



    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')


  END SUBROUTINE update_assimilation_waterlimited


  ! -------------------------------------------------------------------------------------------------------
  !>
  !! Implementation of "aggregate" for task "assimi"
  !!
  !! @param[in,out] tile    Tile for which routine is executed.
  !! @param[in]     config  Vector of process configurations.
  !! @param[in]     options Additional run-time parameters.
  !!
  SUBROUTINE aggregate_assimi(tile, options)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options
    TYPE(t_jsb_model), POINTER                :: model

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':aggregate_assimi'

    ! Local variables
    INTEGER  :: iblk, ics, ice, nc
    !X if necessary: REAL(wp) :: dtime, steplen

    ! Get local variables from options argument
    iblk = options%iblk
    ics  = options%ics
    ice  = options%ice
    nc   = options%nc
    !X if necessary: dtime   = options%dtime
    !X if necessary: steplen = options%steplen

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    model => Get_model(tile%owner_model_id)

    !X Implementation: Start your process scheme here...
    !R: Das ist ein Provisorium bisher. Es muss noch festgelegt werden ob und wie die einzelnen assimi Tasks aggregiert werden sollen.

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE aggregate_assimi


#endif
END MODULE mo_assimi_interface
