!>
!! @brief Initialization of the the vegetation memory.
!!
!! @author
!!  Jan Engel, MPI-BGC Jena
!!
!! @par Revision History
!! First version by Jan Engel (20161213)
!!
MODULE mo_veg_init

#include "dsl_jsbach.inc"

  USE mo_kind,                    ONLY: wp
  USE mo_exception,               ONLY: message, finish, message_text

  USE mo_jsb_class,               ONLY: Get_model
  USE mo_base_memory,             ONLY: t_slm_tile, t_model, t_jsb_task_options     ! *slmtemp
  USE mo_jsb_lctlib_class,        ONLY: t_lctlib_element
  USE mo_jsb_process_class,       ONLY: A2L_, L2A_, SEB_, TURB_, SOIL_, HYDRO_, VEG_, HD_, RAD_, ASSIMI_, PHENO_, SB_
  USE mo_veg_config_class,        ONLY: veg_config

  ! Use of process memories
  dsl4jsb_Use_memory(VEG_)        !  USE mo_veg_memory_class,      ONLY: t_veg_memory

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: veg_init

  CHARACTER(len=*), PARAMETER :: modname = 'mo_veg_init'

CONTAINS

  !-----------------------------------------------------------------------------------------------------
  !> Intialize vegetation process (after memory has been set up)
  !! 
  !! differtiates between the various SLM models: land, plant, canopy, soil
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE veg_init(tile, options)

    USE mo_base_config_class,       ONLY: base_config

    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    TYPE(t_slm_tile),           INTENT(inout)     :: tile         !< data structure for dynamic phenoetation and soil variables
    TYPE(t_jsb_task_options),   INTENT(in)        :: options      !< model options
    ! ---------------------------
    ! 0.2 Local
    TYPE(t_model),          POINTER              :: model
    TYPE(t_lctlib_element), POINTER              :: lctlib     !< land-cover-type library - parameter across pft's
    INTEGER                                      :: iblk, ics, ice, nc
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':veg_init'
    ! ---------------------------
    ! 0.3 Declare Memory
    ! Declare process configuration and memory Pointers
!    dsl4jsb_Def_memory(VEG_)
    ! Declare pointers to variables in memory
    ! 
    ! Get local variables from options argument
    iblk    = options%iblk
    ics     = options%ics
    ice     = options%ice
    nc      = options%nc
    ! ---------------------------
    ! 0.4 Process Activity, Debug Option
!    IF (.NOT. tile%Is_process_active(VEG_)) RETURN
!    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')
    ! ---------------------------
    ! 0.5 Get Memory
    ! Get process config
!    dsl4jsb_Get_config(VEG_)
    ! Get process memories
!    dsl4jsb_Get_memory(VEG_)
    ! Get process variables (Set pointers to variables in memory)
!    model  => Get_model(tile%owner_model_id)
!    lctlib => model%lctlib(tile%lcts%lib_id)
    ! ---------------------------


    ! ------------------------------------------------------------------------------------------------------------
    ! Go science
    ! ------------------------------------------------------------------------------------------------------------


    !< 1.0 run vegetation init routine for an specific SLM model
    !! 
    !! SLM models: land, plant,soil (case 1) and canopy (case 2)
    SELECT CASE(TRIM(base_config%slm_model_name))
      CASE('land', 'plant', 'soil')
        CALL init_slm_plant(tile, options)

      CASE('canopy')
        CALL init_slm_canopy(tile, options)

      CASE('test_canopy')
        ! init is done within the run_test_canopy_model(), which is reasonable because this model runs for only one timestep

      CASE DEFAULT
        CALL finish('veg_init','Invalid SLM model name. Try either land, plant, soil, canopy or test_canopy.')
     END SELECT


  END SUBROUTINE veg_init


  !-----------------------------------------------------------------------------------------------------
  ! vegetation init routine specific to the land model
  ! 
  ! tbc ...
  !-----------------------------------------------------------------------------------------------------

  !-----------------------------------------------------------------------------------------------------
  ! vegetation init routine specific to the soil model
  ! 
  ! tbc ...
  !-----------------------------------------------------------------------------------------------------



  !-----------------------------------------------------------------------------------------------------
  !> vegetation init routine specific to the plant model
  !!
  !! currently used for: SLM land plant and soil
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE init_slm_plant(tile, options)

    USE mo_assimi_constants,                ONLY: CiCa_default_C3, CiCa_default_C4, ic3phot
    USE mo_atmland_constants
    USE mo_veg_constants,                   ONLY: fresp_growth
    USE mo_shared_isotope,                  ONLY: calc_fractionation 
    USE mo_assimi_update_canopy_fluxes,     ONLY: discrimination_ps 

    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    TYPE(t_slm_tile),           INTENT(inout)     :: tile         !< data structure for dynamic phenoetation and soil variables
    TYPE(t_jsb_task_options),   INTENT(in)        :: options      !< model options
    ! ---------------------------
    ! 0.2 Local
    TYPE(t_model),          POINTER              :: model
    TYPE(t_lctlib_element), POINTER              :: lctlib     !< land-cover-type library - parameter across pft's
    REAL(wp)                                     :: hlp1       ! needs to be Scalar as long as pools are Scalar !
    INTEGER                                      :: iblk, ics, ice, nc
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':init_slm_plant'
    ! ---------------------------
    ! 0.3 Declare Memory
    ! Declare process configuration and memory Pointers
    dsl4jsb_Def_memory(VEG_)
    ! Declare pointers to variables in memory
    dsl4jsb_Real2D_onChunk      :: target_cn_leaf                   ! VEG_
    dsl4jsb_Real2D_onChunk      :: target_np_leaf
    dsl4jsb_Real2D_onChunk      :: dens_ind
    dsl4jsb_Real2D_onChunk      :: dphi
    dsl4jsb_Real2D_onChunk      :: npp_talloc_mavg
    dsl4jsb_Real2D_onChunk      :: n_fixation_talloc_mavg
    dsl4jsb_Real2D_onChunk      :: transpiration_talloc_mavg
    dsl4jsb_Real2D_onChunk      :: unit_npp_talloc_mavg
    dsl4jsb_Real2D_onChunk      :: unit_transpiration_talloc_mavg
    dsl4jsb_Real2D_onChunk      :: unit_uptake_n_act_talloc_mavg
    dsl4jsb_Real2D_onChunk      :: unit_uptake_p_act_talloc_mavg
    dsl4jsb_Real2D_onChunk      :: dphi_talloc_mavg
    dsl4jsb_Real2D_onChunk      :: growth_cn_talloc_mavg
    dsl4jsb_Real2D_onChunk      :: growth_cp_talloc_mavg
    dsl4jsb_Real2D_onChunk      :: labile_carbon_talloc_mavg
    dsl4jsb_Real2D_onChunk      :: labile_nitrogen_talloc_mavg
    dsl4jsb_Real2D_onChunk      :: labile_phosphorus_talloc_mavg
    dsl4jsb_Real2D_onChunk      :: reserve_carbon_talloc_mavg
    dsl4jsb_Real2D_onChunk      :: reserve_nitrogen_talloc_mavg
    dsl4jsb_Real2D_onChunk      :: reserve_phosphorus_talloc_mavg
    dsl4jsb_Real2D_onChunk      :: t_air_tfrac_mavg
    dsl4jsb_Real2D_onChunk      :: press_srf_tfrac_mavg
    dsl4jsb_Real2D_onChunk      :: co2_mixing_ratio_tfrac_mavg
    dsl4jsb_Real2D_onChunk      :: ga_tfrac_mavg
    dsl4jsb_Real2D_onChunk      :: beta_sinklim_ps_tfrac_mavg
    dsl4jsb_Real2D_onChunk      :: t_air_tcnl_mavg
    dsl4jsb_Real2D_onChunk      :: press_srf_tcnl_mavg
    dsl4jsb_Real2D_onChunk      :: co2_mixing_ratio_tcnl_mavg
    dsl4jsb_Real2D_onChunk      :: beta_sinklim_ps_tcnl_mavg
    dsl4jsb_Real3D_onChunk      :: fleaf_sunlit_tfrac_mavg_cl
    dsl4jsb_Real3D_onChunk      :: fleaf_sunlit_tcnl_mavg_cl
!! *slmtemp    
    TYPE(t_veg_memory),      POINTER :: veg         ! local pointer to veg__mem for SLM usability
    ! Get local variables from options argument
    iblk    = options%iblk
    ics     = options%ics
    ice     = options%ice
    nc      = options%nc
    ! ---------------------------
    ! 0.4 Process Activity, Debug Option
!    IF (.NOT. tile%Is_process_active(VEG_)) RETURN
!    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')
    ! ---------------------------
    ! 0.5 Get Memory
    ! Get process config
!    dsl4jsb_Get_config(VEG_)
    ! Get process memories
    dsl4jsb_Get_memory(VEG_)
    ! Get process variables (Set pointers to variables in memory)
    model  => Get_model(tile%owner_model_id)
    lctlib => model%lctlib(tile%lcts%lib_id)
    ! ---------------------------
    dsl4jsb_Get_var2D_onChunk(VEG_, target_cn_leaf)
    dsl4jsb_Get_var2D_onChunk(VEG_, target_np_leaf)
    dsl4jsb_Get_var2D_onChunk(VEG_, dens_ind)
    dsl4jsb_Get_var2D_onChunk(VEG_, dphi)
    dsl4jsb_Get_var2D_onChunk(VEG_, npp_talloc_mavg)
    dsl4jsb_Get_var2D_onChunk(VEG_, n_fixation_talloc_mavg)
    dsl4jsb_Get_var2D_onChunk(VEG_, transpiration_talloc_mavg)
    dsl4jsb_Get_var2D_onChunk(VEG_, unit_npp_talloc_mavg)
    dsl4jsb_Get_var2D_onChunk(VEG_, unit_transpiration_talloc_mavg)
    dsl4jsb_Get_var2D_onChunk(VEG_, unit_uptake_n_act_talloc_mavg)
    dsl4jsb_Get_var2D_onChunk(VEG_, unit_uptake_p_act_talloc_mavg)
    dsl4jsb_Get_var2D_onChunk(VEG_, dphi_talloc_mavg)
    dsl4jsb_Get_var2D_onChunk(VEG_, growth_cn_talloc_mavg)
    dsl4jsb_Get_var2D_onChunk(VEG_, growth_cp_talloc_mavg)
    dsl4jsb_Get_var2D_onChunk(VEG_, labile_carbon_talloc_mavg)
    dsl4jsb_Get_var2D_onChunk(VEG_, labile_nitrogen_talloc_mavg)
    dsl4jsb_Get_var2D_onChunk(VEG_, labile_phosphorus_talloc_mavg)
    dsl4jsb_Get_var2D_onChunk(VEG_, reserve_carbon_talloc_mavg)
    dsl4jsb_Get_var2D_onChunk(VEG_, reserve_nitrogen_talloc_mavg)
    dsl4jsb_Get_var2D_onChunk(VEG_, reserve_phosphorus_talloc_mavg)
    dsl4jsb_Get_var2D_onChunk(VEG_, t_air_tfrac_mavg)
    dsl4jsb_Get_var2D_onChunk(VEG_, press_srf_tfrac_mavg)
    dsl4jsb_Get_var2D_onChunk(VEG_, co2_mixing_ratio_tfrac_mavg)
    dsl4jsb_Get_var2D_onChunk(VEG_, ga_tfrac_mavg)
    dsl4jsb_Get_var2D_onChunk(VEG_, beta_sinklim_ps_tfrac_mavg)
    dsl4jsb_Get_var2D_onChunk(VEG_, t_air_tcnl_mavg)
    dsl4jsb_Get_var2D_onChunk(VEG_, press_srf_tcnl_mavg)
    dsl4jsb_Get_var2D_onChunk(VEG_, co2_mixing_ratio_tcnl_mavg)
    dsl4jsb_Get_var2D_onChunk(VEG_, beta_sinklim_ps_tcnl_mavg)
    dsl4jsb_Get_var3D_onChunk(VEG_, fleaf_sunlit_tfrac_mavg_cl)
    dsl4jsb_Get_var3D_onChunk(VEG_, fleaf_sunlit_tcnl_mavg_cl)


!! *slmtemp
    veg     => veg__mem 


    ! ------------------------------------------------------------------------------------------------------------
    ! Go science
    ! ------------------------------------------------------------------------------------------------------------


      target_cn_leaf(:)             = lctlib%cn_leaf
      target_np_leaf(:)             = lctlib%np_leaf

      IF(veg_config%flag_tile_veg_dynamics)THEN
         veg%init_pool%carbon%ptr         = lctlib%seed_size * 10.0_wp
         dens_ind(:)                 = 0.0_wp 
      ELSE
         veg%init_pool%carbon%ptr         = 50.0_wp 
         dens_ind(:)                 = 1.0_wp 
      ENDIF

      veg%init_pool%nitrogen    = veg%init_pool%carbon    / (1._wp + fresp_growth) / lctlib%cn_leaf
      veg%init_pool%phosphorus  = veg%init_pool%nitrogen  / lctlib%np_leaf

      ! init_pool%carbon13
      IF(lctlib%ps_pathway == ic3phot)THEN
         hlp1 = discrimination_ps(CiCa_default_C3 * def_co2_mixing_ratio, def_co2_mixing_ratio, 'C13',lctlib%ps_pathway)
      ELSE
         hlp1 = discrimination_ps(CiCa_default_C4 * def_co2_mixing_ratio, def_co2_mixing_ratio, 'C13',lctlib%ps_pathway)
      ENDIF
      veg%init_pool%carbon13 = veg%init_pool%carbon * calc_fractionation(def_co2_mixing_ratio, def_co2_mixing_ratio_C13, hlp1)


      ! init_pool%carbon14
      IF(lctlib%ps_pathway == ic3phot)THEN
         hlp1 = discrimination_ps(CiCa_default_C3 * def_co2_mixing_ratio, def_co2_mixing_ratio, 'C14',lctlib%ps_pathway)
      ELSE
         hlp1 = discrimination_ps(CiCa_default_C4 * def_co2_mixing_ratio, def_co2_mixing_ratio, 'C14',lctlib%ps_pathway)
      ENDIF
      veg%init_pool%carbon14 = veg%init_pool%carbon * calc_fractionation(def_co2_mixing_ratio, def_co2_mixing_ratio_C14,hlp1)

      ! IF start bareground with vegetation dynamics switched on
      IF(veg_config%flag_tile_veg_dynamics)THEN
         veg%pool%seed_bed = veg%init_pool 
      ELSE
         veg%pool%reserve  = veg%init_pool 
      END IF

      ! plant variables
      dphi(:) = - lctlib%phi_leaf_min

      ! initialising moving averages
      ! allocation
      npp_talloc_mavg(:) = 0.0_wp
      n_fixation_talloc_mavg(:) = 0.0_wp
      transpiration_talloc_mavg(:) = 0.0_wp
      unit_npp_talloc_mavg(:) = 0.0_wp
      unit_transpiration_talloc_mavg(:) = 0.0_wp
      unit_uptake_n_act_talloc_mavg(:) = 0.0_wp
      unit_uptake_p_act_talloc_mavg(:) = 0.0_wp
      dphi_talloc_mavg(:) =  - lctlib%phi_leaf_min
      growth_cn_talloc_mavg(:) = 60.0_wp
      growth_cp_talloc_mavg(:) = 2000.0_wp
      labile_carbon_talloc_mavg(:) = 0.0_wp
      labile_nitrogen_talloc_mavg(:) = 0.0_wp
      labile_phosphorus_talloc_mavg(:) = 0.0_wp
      reserve_carbon_talloc_mavg(:) = 0.0_wp
      reserve_nitrogen_talloc_mavg(:) = 0.0_wp
      reserve_phosphorus_talloc_mavg(:) = 0.0_wp

      ! N fraction optimisation
      fleaf_sunlit_tfrac_mavg_cl(:,:) = 0.5_wp
      t_air_tfrac_mavg(:) = 283.15_wp
      press_srf_tfrac_mavg(:) = standard_press_srf
      co2_mixing_ratio_tfrac_mavg(:) = def_co2_mixing_ratio
      ga_tfrac_mavg(:) = 0.2_wp
      beta_sinklim_ps_tfrac_mavg(:) = 1.0_wp
      
      ! leaf N optimisation
      fleaf_sunlit_tcnl_mavg_cl(:,:) = 0.5_wp
      t_air_tcnl_mavg(:) = 283.15_wp
      press_srf_tcnl_mavg(:) = standard_press_srf
      co2_mixing_ratio_tcnl_mavg(:) = def_co2_mixing_ratio
      ga_tfrac_mavg(:) = 0.2_wp
      beta_sinklim_ps_tcnl_mavg(:) = 1.0_wp


  END SUBROUTINE init_slm_plant


  !-----------------------------------------------------------------------------------------------------
  !> init routine specific to the canopy model
  !!
  !! called from init_slm
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE init_slm_canopy(tile, options)

    USE mo_shared_constants,                ONLY: molar_mass_N

    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    TYPE(t_slm_tile),           INTENT(inout)     :: tile         !< data structure for dynamic phenoetation and soil variables
    TYPE(t_jsb_task_options),   INTENT(in)        :: options      !< model options
    ! ---------------------------
    ! 0.2 Local
    TYPE(t_model),          POINTER              :: model
    TYPE(t_lctlib_element), POINTER              :: lctlib     !< land-cover-type library - parameter across pft's
    REAL(wp)                                     :: ref_lai    !< LAI
    REAL(wp)                                     :: ref_nleaf  !< N leaf per unit leaf area in gN/m2
    INTEGER                                      :: iblk, ics, ice, nc
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':init_slm_canopy'
    ! ---------------------------
    ! 0.3 Declare Memory
    ! Declare process configuration and memory Pointers
    dsl4jsb_Def_memory(VEG_)
    ! Declare pointers to variables in memory
    dsl4jsb_Real2D_onChunk      :: lai              ! VEG_
!! *slmtemp    
    TYPE(t_veg_memory),      POINTER :: veg         ! local pointer to veg__mem for SLM usability
    ! Get local variables from options argument
    iblk    = options%iblk
    ics     = options%ics
    ice     = options%ice
    nc      = options%nc
    ! ---------------------------
    ! 0.4 Process Activity, Debug Option
!    IF (.NOT. tile%Is_process_active(VEG_)) RETURN
!    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')
    ! ---------------------------
    ! 0.5 Get Memory
    ! Get process config
!    dsl4jsb_Get_config(VEG_)
    ! Get process memories
    dsl4jsb_Get_memory(VEG_)
    ! Get process variables (Set pointers to variables in memory)
    model  => Get_model(tile%owner_model_id)
    lctlib => model%lctlib(tile%lcts%lib_id)
    ! ---------------------------
    dsl4jsb_Get_var2D_onChunk(VEG_, lai)


!! *slmtemp
    veg     => veg__mem 


    ! ------------------------------------------------------------------------------------------------------------
    ! Go science
    ! ------------------------------------------------------------------------------------------------------------


    !< 1.0 init local variables
    !! 
    ref_lai   = 6._wp
    ref_nleaf = 2.0_wp

    !< 2.0 Initialise tile
    !! 
    !! Initialisation of soil and vegetation, so far relevant
    lai(:) = ref_lai 

    IF(veg_config%pft_id==0) THEN
      veg%pool%leaf%nitrogen%ptr   = ref_nleaf / molar_mass_N * lai(1)
    ELSE
      veg%pool%leaf%nitrogen%ptr   = lai(1) / lctlib%sla / lctlib%cn_leaf 
      veg%pool%leaf%phosphorus = veg%pool%leaf%nitrogen / lctlib%np_leaf
    ENDIF

  END SUBROUTINE init_slm_canopy

















END MODULE mo_veg_init
