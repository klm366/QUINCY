!>
!! @brief Contains the routines for reading parameters from parameter_slm_run.list and replace parameter-values as defined in mo_PROC_constants
!!
!! @author
!!  Jan Engel, MPI-BGC Jena 
!!
!! @par Revision History
!! First version by Jan Engel (2016-10-04)
!!
MODULE mo_shared_set_parameters

  USE mo_kind,                ONLY: wp, eps8
  USE mo_exception,           ONLY: message_text, message, finish
  USE mo_jsb_io,              ONLY: missval  ! = -9.E+15_wp

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: parameters_slm_from_file                        ! it's an object
  PUBLIC :: set_parameter_lctlib_values_with_list_values    ! it's a subroutine

  !-----------------------------------------------------------------------------------------------------
  !> SLM process parameter values from file parameter_slm_run.list - replacing the values from mo_PROC_constants
  !! 
  !-----------------------------------------------------------------------------------------------------
  TYPE t_parameters_slm_from_file
    !LOGICAL            :: active   ! de/activate this TYPE (this is non-functional for now 20161004)
    CHARACTER(LEN=30)  :: filename = 'parameter_slm_run.list'                 ! filename 
    ! SHARED_CTL
    REAL(wp)           :: mavg_period_tlabile
    REAL(wp)           :: mavg_period_tphen
    REAL(wp)           :: mavg_period_tuptake
    REAL(wp)           :: mavg_period_tfrac
    REAL(wp)           :: mavg_period_tcnl
    REAL(wp)           :: mavg_period_talloc
    REAL(wp)           :: mavg_period_tvegdyn
    REAL(wp)           :: mavg_period_tacclim
    REAL(wp)           :: mavg_period_tmic
    ! VEGETATION_CTL
    REAL(wp)           :: delta_n_fraction
    REAL(wp)           :: delta_n_leaf
    REAL(wp)           :: cost_n_fixation_max
    REAL(wp)           :: cost_n_fixation_min
    REAL(wp)           :: vmax_bnf                        !< maximum rate of N fixation
    REAL(wp)           :: k_bnf                           !< half_saturating cost of the response of N fixation
    REAL(wp)           :: offset_bnf                      !< parameter used for bnf-scheme 'resistor'
    REAL(wp)           :: f_resp_acclim
    REAL(wp)           :: fmaint_rate_base
    ! SB_CTL
    REAL(wp)           :: floss_nmin
    REAL(wp)           :: tau_fast
    REAL(wp)           :: tau_slow
    REAL(wp)           :: tau_structural_litter
    REAL(wp)           :: microbial_cue
    REAL(wp)           :: microbial_cue_min
    REAL(wp)           :: microbial_cue_max
    REAL(wp)           :: k_diff_org
    ! additional jsm parameters
    REAL(wp)           :: vmax_mic_uptake_nh4, &    
                          km_mic_uptake_nh4, &      
                          vmax_mic_uptake_no3, &    
                          km_mic_uptake_no3, &      
                          f_mic_carrier_n, &        
                          vmax_plant_uptake_nh4, &  
                          km_plant_uptake_nh4, &    
                          vmax_plant_uptake_no3, &  
                          km_plant_uptake_no3, &    
                          f_plant_carrier_n, &      
                          vmax_mic_uptake_po4, &    
                          km_mic_uptake_po4, &      
                          vmax_plant_uptake_po4, &  
                          km_plant_uptake_po4, &    
                          k_mic_uptake_nh4, &       
                          k_plant_uptake_nh4, &     
                          k_mic_uptake_no3, &       
                          k_plant_uptake_no3, &     
                          k_mic_uptake_po4, &       
                          k_plant_uptake_po4, &     
                          f_mic_carrier_po4, &      
                          f_plant_carrier_po4, &    
                          qmax_po4_mineral, &       
                          km_adsorpt_mineral_po4 
  CONTAINS
    PROCEDURE :: read_param => read_parameter_slm_list
    PROCEDURE :: set_param  => set_parameter_slm_values_with_list_values
  END TYPE t_parameters_slm_from_file

  TYPE(t_parameters_slm_from_file), SAVE :: parameters_slm_from_file


  !-----------------------------------------------------------------------------------------------------
  !> lctlib parameter values from file parameter_slm_run.list - replacing the values from lctlib_slm_*.def file
  !! 
  !-----------------------------------------------------------------------------------------------------
  TYPE t_parameters_lctlib_from_file
    !LOGICAL            :: active   
    CHARACTER(LEN=30)   :: filename = 'parameter_slm_run.list'                 ! filename 

    REAL(wp)            :: sla
    REAL(wp)            :: sigma_vis
    REAL(wp)            :: sigma_nir
    REAL(wp)            :: omega_clumping
    REAL(wp)            :: crown_shape_factor
    REAL(wp)            :: cn_leaf
    REAL(wp)            :: np_leaf
    REAL(wp)            :: k0_fn_struc
    REAL(wp)            :: fn_oth_min
    REAL(wp)            :: t_jmax_opt
    REAL(wp)            :: t_jmax_omega
    REAL(wp)            :: g0
    REAL(wp)            :: g1
    REAL(wp)            :: gmin
    REAL(wp)            :: lambda_soil_ps
    REAL(wp)            :: weibull_k_soil_ps
    REAL(wp)            :: lambda_soil_gs
    REAL(wp)            :: weibull_k_soil_gs
    ! turnover times
    REAL(wp)            :: tau_leaf
    REAL(wp)            :: tau_fine_root
    REAL(wp)            :: tau_coarse_root
    REAL(wp)            :: tau_branch
    REAL(wp)            :: tau_sap_wood
    REAL(wp)            :: tau_fruit
    REAL(wp)            :: tau_seed_litter
    REAL(wp)            :: tau_seed_est
    REAL(wp)            :: tau_mycorrhiza
    ! N uptake parameters
    REAL(wp)            :: vmax_uptake_n
    REAL(wp)            :: vmax_uptake_p
    REAL(wp)            :: bnf_base
    ! Vegetation dynamics parameters
    REAL(wp)            :: lambda_est_light
    REAL(wp)            :: k_est_light
    REAL(wp)            :: seed_size
    ! Phenology parameters
    REAL(wp)            :: beta_soil_flush
    REAL(wp)            :: beta_soil_senescence
    REAL(wp)            :: gdd_t_air_threshold
    REAL(wp)            :: gdd_req_max
    REAL(wp)            :: k_gdd_dormance
    REAL(wp)            :: t_air_senescence
    REAL(wp)            :: min_leaf_age
    ! Allocation paramters
    REAL(wp)            :: frac_sapwood_branch
    REAL(wp)            :: wood_density
    REAL(wp)            :: k_latosa
    REAL(wp)            :: k_crtos
    REAL(wp)            :: k_rtos
    REAL(wp)            :: k2_fruit_alloc
    REAL(wp)            :: allom_k1
    REAL(wp)            :: allom_k2
    REAL(wp)            :: phi_leaf_min
    REAL(wp)            :: k_root
    REAL(wp)            :: k_sapwood
    REAL(wp)            :: c0_allom
    REAL(wp)            :: fstore_target
    ! Soil
    REAL(wp)            :: k_root_dist

  END TYPE t_parameters_lctlib_from_file

  TYPE(t_parameters_lctlib_from_file) :: param_pft1_from_file,  &
                                         param_pft2_from_file,  &
                                         param_pft3_from_file,  &
                                         param_pft4_from_file,  &
                                         param_pft5_from_file,  &
                                         param_pft6_from_file,  &
                                         param_pft7_from_file,  &
                                         param_pft8_from_file


  CHARACTER(len=*), PARAMETER :: modname = 'mo_shared_set_parameters'

CONTAINS


  !-----------------------------------------------------------------------------------------------------
  !> read SLM parameter values from file
  !! 
  !! 
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE read_parameter_slm_list(this)

    USE mo_namelist,           ONLY: open_nml, POSITIONED, position_nml, close_nml  ! read model options from file

    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    CLASS(t_parameters_slm_from_file), INTENT(inout) :: this
    ! ---------------------------
    ! 0.2 Local
    ! SHARED_CTL
    REAL(wp)           :: mavg_period_tlabile
    REAL(wp)           :: mavg_period_tphen
    REAL(wp)           :: mavg_period_tuptake
    REAL(wp)           :: mavg_period_tfrac
    REAL(wp)           :: mavg_period_tcnl
    REAL(wp)           :: mavg_period_talloc
    REAL(wp)           :: mavg_period_tvegdyn
    REAL(wp)           :: mavg_period_tacclim
    REAL(wp)           :: mavg_period_tmic
    ! VEGETATION_CTL
    REAL(wp)           :: delta_n_fraction
    REAL(wp)           :: delta_n_leaf
    REAL(wp)           :: cost_n_fixation_max
    REAL(wp)           :: cost_n_fixation_min
    REAL(wp)           :: vmax_bnf
    REAL(wp)           :: k_bnf
    REAL(wp)           :: offset_bnf
    REAL(wp)           :: f_resp_acclim
    REAL(wp)           :: fmaint_rate_base
    ! SB_CTL
    REAL(wp)           :: floss_nmin
    REAL(wp)           :: tau_fast
    REAL(wp)           :: tau_slow
    REAL(wp)           :: tau_structural_litter
    REAL(wp)           :: microbial_cue
    REAL(wp)           :: microbial_cue_min
    REAL(wp)           :: microbial_cue_max
    REAL(wp)           :: k_diff_org
    ! additional jsm parameters
    REAL(wp)           :: vmax_mic_uptake_nh4, &    
                          km_mic_uptake_nh4, &      
                          vmax_mic_uptake_no3, &    
                          km_mic_uptake_no3, &      
                          f_mic_carrier_n, &        
                          vmax_plant_uptake_nh4, &  
                          km_plant_uptake_nh4, &    
                          vmax_plant_uptake_no3, &  
                          km_plant_uptake_no3, &    
                          f_plant_carrier_n, &      
                          vmax_mic_uptake_po4, &    
                          km_mic_uptake_po4, &      
                          vmax_plant_uptake_po4, &  
                          km_plant_uptake_po4, &    
                          k_mic_uptake_nh4, &       
                          k_plant_uptake_nh4, &     
                          k_mic_uptake_no3, &       
                          k_plant_uptake_no3, &     
                          k_mic_uptake_po4, &       
                          k_plant_uptake_po4, &     
                          f_mic_carrier_po4, &      
                          f_plant_carrier_po4, &    
                          qmax_po4_mineral, &       
                          km_adsorpt_mineral_po4 

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':read_parameter_slm_list'

    ! Namelist definition of section parameters
    ! 
    NAMELIST /shared_ctl/         &  
        mavg_period_tlabile       , &
        mavg_period_tphen         , &
        mavg_period_tuptake       , &
        mavg_period_tfrac         , &
        mavg_period_tcnl          , &
        mavg_period_talloc        , &
        mavg_period_tvegdyn       , &
        mavg_period_tacclim       , &
        mavg_period_tmic

    NAMELIST /vegetation_ctl/     &  
        delta_n_fraction          , &
        delta_n_leaf              , &
        cost_n_fixation_max       , &
        cost_n_fixation_min       , &
        vmax_bnf                  , &
        k_bnf                     , &
        offset_bnf                , &
        f_resp_acclim             , &
        fmaint_rate_base

    NAMELIST /sb_ctl/             &
        floss_nmin                , &
        tau_fast                  , &
        tau_slow                  , &
        tau_structural_litter     , &
        microbial_cue             , & 
        microbial_cue_min         , &      
        microbial_cue_max         , &
        k_diff_org                , &
        vmax_mic_uptake_nh4       , &    
        km_mic_uptake_nh4         , &      
        vmax_mic_uptake_no3       , &    
        km_mic_uptake_no3         , &      
        f_mic_carrier_n           , &        
        vmax_plant_uptake_nh4     , &  
        km_plant_uptake_nh4       , &    
        vmax_plant_uptake_no3     , &  
        km_plant_uptake_no3       , &    
        f_plant_carrier_n         , &      
        vmax_mic_uptake_po4       , &    
        km_mic_uptake_po4         , &      
        vmax_plant_uptake_po4     , &  
        km_plant_uptake_po4       , &    
        k_mic_uptake_nh4          , &       
        k_plant_uptake_nh4        , &     
        k_mic_uptake_no3          , &       
        k_plant_uptake_no3        , &     
        k_mic_uptake_po4          , &       
        k_plant_uptake_po4        , &     
        f_mic_carrier_po4         , &      
        f_plant_carrier_po4       , &    
        qmax_po4_mineral          , &       
        km_adsorpt_mineral_po4 


    ! variables for reading model-options from namelist
    ! 
    INTEGER  :: nml_handler, nml_unit, istat

    CALL message(TRIM(routine), 'Starting to read parameter list')

    ! Set defaults
    ! 
    ! SHARED_CTL
    mavg_period_tlabile         = missval
    mavg_period_tphen           = missval
    mavg_period_tuptake         = missval
    mavg_period_tfrac           = missval
    mavg_period_tcnl            = missval
    mavg_period_talloc          = missval
    mavg_period_tvegdyn         = missval
    mavg_period_tacclim         = missval
    mavg_period_tmic            = missval
    ! VEGETATION_CTL
    delta_n_fraction            = missval
    delta_n_leaf                = missval
    cost_n_fixation_max         = missval
    cost_n_fixation_min         = missval
    vmax_bnf                    = missval
    k_bnf                       = missval
    offset_bnf                  = missval
    f_resp_acclim               = missval
    fmaint_rate_base            = missval
    ! SB_CTL
    floss_nmin                  = missval
    tau_fast                    = missval
    tau_slow                    = missval
    tau_structural_litter       = missval
    microbial_cue               = missval
    microbial_cue_min           = missval
    microbial_cue_max           = missval
    k_diff_org                  = missval
    vmax_mic_uptake_nh4         = missval    
    km_mic_uptake_nh4           = missval      
    vmax_mic_uptake_no3         = missval    
    km_mic_uptake_no3           = missval      
    f_mic_carrier_n             = missval        
    vmax_plant_uptake_nh4       = missval  
    km_plant_uptake_nh4         = missval    
    vmax_plant_uptake_no3       = missval  
    km_plant_uptake_no3         = missval    
    f_plant_carrier_n           = missval      
    vmax_mic_uptake_po4         = missval    
    km_mic_uptake_po4           = missval      
    vmax_plant_uptake_po4       = missval  
    km_plant_uptake_po4         = missval    
    k_mic_uptake_nh4            = missval       
    k_plant_uptake_nh4          = missval     
    k_mic_uptake_no3            = missval       
    k_plant_uptake_no3          = missval     
    k_mic_uptake_po4            = missval       
    k_plant_uptake_po4          = missval     
    f_mic_carrier_po4           = missval      
    f_plant_carrier_po4         = missval    
    qmax_po4_mineral            = missval       
    km_adsorpt_mineral_po4      = missval

    ! read the namelist
    ! 
    ! &SHARED_CTL
    nml_handler = open_nml(TRIM(this%filename))
    nml_unit = position_nml('SHARED_CTL', nml_handler, STATUS=istat)
    IF (istat == POSITIONED) READ(nml_unit, shared_ctl)
    CALL close_nml(nml_handler)
    ! &VEGETATION_CTL
    nml_handler = open_nml(TRIM(this%filename))
    nml_unit = position_nml('VEGETATION_CTL', nml_handler, STATUS=istat)
    IF (istat == POSITIONED) READ(nml_unit, vegetation_ctl)
    CALL close_nml(nml_handler)
    ! &SB_CTL
    nml_handler = open_nml(TRIM(this%filename))
    nml_unit = position_nml('SB_CTL', nml_handler, STATUS=istat)
    IF (istat == POSITIONED) READ(nml_unit, sb_ctl)
    CALL close_nml(nml_handler)

    ! pass values as read from file
    ! 
    ! SHARED_CTL
    this%mavg_period_tlabile         = mavg_period_tlabile
    this%mavg_period_tphen           = mavg_period_tphen
    this%mavg_period_tuptake         = mavg_period_tuptake
    this%mavg_period_tfrac           = mavg_period_tfrac
    this%mavg_period_tcnl            = mavg_period_tcnl
    this%mavg_period_talloc          = mavg_period_talloc
    this%mavg_period_tvegdyn         = mavg_period_tvegdyn
    this%mavg_period_tacclim         = mavg_period_tacclim
    this%mavg_period_tmic            = mavg_period_tmic
    ! VEGETATION_CTL
    this%delta_n_fraction            = delta_n_fraction
    this%delta_n_leaf                = delta_n_leaf
    this%cost_n_fixation_max         = cost_n_fixation_max
    this%cost_n_fixation_min         = cost_n_fixation_min
    this%vmax_bnf                    = vmax_bnf
    this%k_bnf                       = k_bnf
    this%offset_bnf                  = offset_bnf
    this%f_resp_acclim               = f_resp_acclim
    this%fmaint_rate_base            = fmaint_rate_base 
    ! SB_CTL
    this%floss_nmin                  = floss_nmin
    this%tau_fast                    = tau_fast
    this%tau_slow                    = tau_slow
    this%tau_structural_litter       = tau_structural_litter
    this%microbial_cue               = microbial_cue
    this%microbial_cue_min           = microbial_cue_min
    this%microbial_cue_max           = microbial_cue_max
    this%k_diff_org                  = k_diff_org
    this%vmax_mic_uptake_nh4         = vmax_mic_uptake_nh4       
    this%km_mic_uptake_nh4           = km_mic_uptake_nh4           
    this%vmax_mic_uptake_no3         = vmax_mic_uptake_no3       
    this%km_mic_uptake_no3           = km_mic_uptake_no3           
    this%f_mic_carrier_n             = f_mic_carrier_n               
    this%vmax_plant_uptake_nh4       = vmax_plant_uptake_nh4   
    this%km_plant_uptake_nh4         = km_plant_uptake_nh4       
    this%vmax_plant_uptake_no3       = vmax_plant_uptake_no3   
    this%km_plant_uptake_no3         = km_plant_uptake_no3       
    this%f_plant_carrier_n           = f_plant_carrier_n           
    this%vmax_mic_uptake_po4         = vmax_mic_uptake_po4       
    this%km_mic_uptake_po4           = km_mic_uptake_po4           
    this%vmax_plant_uptake_po4       = vmax_plant_uptake_po4   
    this%km_plant_uptake_po4         = km_plant_uptake_po4       
    this%k_mic_uptake_nh4            = k_mic_uptake_nh4             
    this%k_plant_uptake_nh4          = k_plant_uptake_nh4         
    this%k_mic_uptake_no3            = k_mic_uptake_no3             
    this%k_plant_uptake_no3          = k_plant_uptake_no3         
    this%k_mic_uptake_po4            = k_mic_uptake_po4             
    this%k_plant_uptake_po4          = k_plant_uptake_po4         
    this%f_mic_carrier_po4           = f_mic_carrier_po4           
    this%f_plant_carrier_po4         = f_plant_carrier_po4       
    this%qmax_po4_mineral            = qmax_po4_mineral             
    this%km_adsorpt_mineral_po4      = km_adsorpt_mineral_po4

  END SUBROUTINE read_parameter_slm_list


  !-----------------------------------------------------------------------------------------------------
  !> replace hard coded SLM parameter values with values from file
  !! 
  !! 
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE set_parameter_slm_values_with_list_values(this)

    USE mo_base_constants,                  ONLY: one_day, one_year, dtime
    !USE mo_assimi_constants
    !USE mo_atmland_constants
    !USE mo_rad_constants
    USE mo_shared_constants
    !USE mo_sse_constants
    USE mo_sb_constants
    USE mo_veg_constants

    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    CLASS(t_parameters_slm_from_file), INTENT(inout) :: this
    ! ---------------------------
    ! 0.2 Local
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':set_parameter_slm_values_with_list_values'


    !> 1.0 read parameter values from file
    !! 
    CALL this%read_param

    !> 2.0 replace parameter values from mo_PROC_constants for parameters provided by parameter_slm_run.list
    !! 
    ! SHARED_CTL
    IF (ABS(this%mavg_period_tlabile    - missval) > eps8)  mavg_period_tlabile   = this%mavg_period_tlabile
    IF (ABS(this%mavg_period_tphen      - missval) > eps8)  mavg_period_tphen     = this%mavg_period_tphen
    IF (ABS(this%mavg_period_tuptake    - missval) > eps8)  mavg_period_tuptake   = this%mavg_period_tuptake
    IF (ABS(this%mavg_period_tfrac      - missval) > eps8)  mavg_period_tfrac     = this%mavg_period_tfrac
    IF (ABS(this%mavg_period_tcnl       - missval) > eps8)  mavg_period_tcnl      = this%mavg_period_tcnl
    IF (ABS(this%mavg_period_talloc     - missval) > eps8)  mavg_period_talloc    = this%mavg_period_talloc
    IF (ABS(this%mavg_period_tvegdyn    - missval) > eps8)  mavg_period_tvegdyn   = this%mavg_period_tvegdyn
    IF (ABS(this%mavg_period_tacclim    - missval) > eps8)  mavg_period_tacclim   = this%mavg_period_tacclim
    IF (ABS(this%mavg_period_tmic       - missval) > eps8)  mavg_period_tmic      = this%mavg_period_tmic
    
    ! VEGETATION_CTL
    IF (ABS(this%delta_n_fraction       - missval) > eps8)  delta_n_fraction      = this%delta_n_fraction
    IF (ABS(this%delta_n_leaf           - missval) > eps8)  delta_n_leaf          = this%delta_n_leaf
    IF (ABS(this%cost_n_fixation_max    - missval) > eps8)  cost_n_fixation_max   = this%cost_n_fixation_max / &
                                                                                      (molar_mass_C/molar_mass_N)
    IF (ABS(this%cost_n_fixation_min    - missval) > eps8)  cost_n_fixation_min   = this%cost_n_fixation_min / &
                                                                                      (molar_mass_C/molar_mass_N)
    IF (ABS(this%vmax_bnf               - missval) > eps8)  vmax_bnf              = this%vmax_bnf * 0.0225_wp / &
                                                                                      (molar_mass_N/molar_mass_C) / &
                                                                                      (one_day*one_year) * 1.e6_wp
    IF (ABS(this%k_bnf                  - missval) > eps8)  k_bnf                 = this%k_bnf / &
                                                                                      (molar_mass_C/molar_Mass_N)
    IF (ABS(this%offset_bnf             - missval) > eps8)  offset_bnf            = this%offset_bnf * & 
                                                                                      (f_nacq_nh4 * transform_cost_nh4 + &
                                                                                       f_nacq_no3 * transform_cost_no3)

    IF (ABS(this%f_resp_acclim          - missval) > eps8)  f_resp_acclim         = this%f_resp_acclim
    IF (ABS(this%fmaint_rate_base       - missval) > eps8)  fmaint_rate_base      = this%fmaint_rate_base

    ! SB_CTL
    IF (ABS(this%floss_nmin             - missval) > eps8)  floss_nmin             = this%floss_nmin
    IF (ABS(this%tau_fast               - missval) > eps8)  tau_fast               = this%tau_fast    * one_year * one_day
    IF (ABS(this%tau_slow               - missval) > eps8)  tau_slow               = this%tau_slow    * one_year * one_day
    IF (ABS(this%tau_structural_litter  - missval) > eps8)  tau_structural_litter  = this%tau_structural_litter * &
                                                                                       one_year * one_day
    IF (ABS(this%microbial_cue          - missval) > eps8)  microbial_cue          = this%microbial_cue
    IF (ABS(this%microbial_cue_min      - missval) > eps8)  microbial_cue_min      = this%microbial_cue_min
    IF (ABS(this%microbial_cue_max      - missval) > eps8)  microbial_cue_max      = this%microbial_cue_max
    IF (ABS(this%k_diff_org             - missval) > eps8)  k_diff_org             = this%k_diff_org / one_year / one_day
    IF (ABS(this%vmax_mic_uptake_nh4    - missval) > eps8)  vmax_mic_uptake_nh4    = this%vmax_mic_uptake_nh4 / &
                                                                                       3600._wp * molar_mass_C
    IF (ABS(this%km_mic_uptake_nh4      - missval) > eps8)  km_mic_uptake_nh4      = this%km_mic_uptake_nh4 / molar_mass_N
    IF (ABS(this%vmax_mic_uptake_no3    - missval) > eps8)  vmax_mic_uptake_no3    = this%vmax_mic_uptake_no3 / &
                                                                                       3600._wp * molar_mass_C
    IF (ABS(this%km_mic_uptake_no3      - missval) > eps8)  km_mic_uptake_no3      = this%km_mic_uptake_no3 / molar_mass_N
    IF (ABS(this%f_mic_carrier_n        - missval) > eps8)  f_mic_carrier_n        = this%f_mic_carrier_n
    IF (ABS(this%vmax_plant_uptake_nh4  - missval) > eps8)  vmax_plant_uptake_nh4  = this%vmax_plant_uptake_nh4 / &
                                                                                       3600._wp * molar_mass_C
    IF (ABS(this%km_plant_uptake_nh4    - missval) > eps8)  km_plant_uptake_nh4    = this%km_plant_uptake_nh4 / molar_mass_N
    IF (ABS(this%vmax_plant_uptake_no3  - missval) > eps8)  vmax_plant_uptake_no3  = this%vmax_plant_uptake_no3 / &
                                                                                       3600._wp * molar_mass_C
    IF (ABS(this%km_plant_uptake_no3    - missval) > eps8)  km_plant_uptake_no3    = this%km_plant_uptake_no3 / molar_mass_N
    IF (ABS(this%f_plant_carrier_n      - missval) > eps8)  f_plant_carrier_n      = this%f_plant_carrier_n
    IF (ABS(this%vmax_mic_uptake_po4    - missval) > eps8)  vmax_mic_uptake_po4    = this%vmax_mic_uptake_po4 / &
                                                                                       3600._wp * molar_mass_C
    IF (ABS(this%km_mic_uptake_po4      - missval) > eps8)  km_mic_uptake_po4      = this%km_mic_uptake_po4 / molar_mass_P
    IF (ABS(this%vmax_plant_uptake_po4  - missval) > eps8)  vmax_plant_uptake_po4  = this%vmax_plant_uptake_po4 / &
                                                                                       3600._wp * molar_mass_C
    IF (ABS(this%km_plant_uptake_po4    - missval) > eps8)  km_plant_uptake_po4    = this%km_plant_uptake_po4 / molar_mass_P
    IF (ABS(this%k_mic_uptake_nh4       - missval) > eps8)  k_mic_uptake_nh4       = this%k_mic_uptake_nh4 / one_day
    IF (ABS(this%k_plant_uptake_nh4     - missval) > eps8)  k_plant_uptake_nh4     = this%k_plant_uptake_nh4 / one_day
    IF (ABS(this%k_mic_uptake_no3       - missval) > eps8)  k_mic_uptake_no3       = this%k_mic_uptake_no3 / one_day
    IF (ABS(this%k_plant_uptake_no3     - missval) > eps8)  k_plant_uptake_no3     = this%k_plant_uptake_no3 / one_day
    IF (ABS(this%k_mic_uptake_po4       - missval) > eps8)  k_mic_uptake_po4       = this%k_mic_uptake_po4 / one_day
    IF (ABS(this%k_plant_uptake_po4     - missval) > eps8)  k_plant_uptake_po4     = this%k_plant_uptake_po4 / one_day
    IF (ABS(this%f_mic_carrier_po4      - missval) > eps8)  f_mic_carrier_po4      = this%f_mic_carrier_po4
    IF (ABS(this%f_plant_carrier_po4    - missval) > eps8)  f_plant_carrier_po4    = this%f_plant_carrier_po4
    IF (ABS(this%qmax_po4_mineral       - missval) > eps8)  qmax_po4_mineral       = this%qmax_po4_mineral / molar_mass_P
    IF (ABS(this%km_adsorpt_mineral_po4 - missval) > eps8)  km_adsorpt_mineral_po4 = this%km_adsorpt_mineral_po4 / 1000._wp / &
                                                                                       molar_mass_P
  END SUBROUTINE set_parameter_slm_values_with_list_values


  !-----------------------------------------------------------------------------------------------------
  !> read lctlib parameter values from file
  !! 
  !! 
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE read_parameter_lctlib_list

    USE mo_namelist,           ONLY: open_nml, POSITIONED, position_nml, close_nml  ! read model options from file

    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    !  just working on the above defined (private) objects 'param_pft1_from_file' etc.
    ! ---------------------------
    ! 0.2 Local
    ! PFT 1
    REAL(wp)            :: sla_p1
    REAL(wp)            :: sigma_vis_p1
    REAL(wp)            :: sigma_nir_p1
    REAL(wp)            :: omega_clumping_p1
    REAL(wp)            :: crown_shape_factor_p1
    REAL(wp)            :: cn_leaf_p1
    REAL(wp)            :: np_leaf_p1
    REAL(wp)            :: k0_fn_struc_p1
    REAL(wp)            :: fn_oth_min_p1
    REAL(wp)            :: t_jmax_opt_p1
    REAL(wp)            :: t_jmax_omega_p1
    REAL(wp)            :: g0_p1
    REAL(wp)            :: g1_p1
    REAL(wp)            :: gmin_p1
    REAL(wp)            :: lambda_soil_ps_p1
    REAL(wp)            :: weibull_k_soil_ps_p1
    REAL(wp)            :: lambda_soil_gs_p1
    REAL(wp)            :: weibull_k_soil_gs_p1
    REAL(wp)            :: tau_leaf_p1
    REAL(wp)            :: tau_fine_root_p1
    REAL(wp)            :: tau_coarse_root_p1
    REAL(wp)            :: tau_branch_p1
    REAL(wp)            :: tau_sap_wood_p1
    REAL(wp)            :: tau_fruit_p1
    REAL(wp)            :: tau_seed_litter_p1
    REAL(wp)            :: tau_seed_est_p1
    REAL(wp)            :: tau_mycorrhiza_p1
    REAL(wp)            :: vmax_uptake_n_p1
    REAL(wp)            :: vmax_uptake_p_p1
    REAL(wp)            :: bnf_base_p1
    REAL(wp)            :: lambda_est_light_p1
    REAL(wp)            :: k_est_light_p1
    REAL(wp)            :: seed_size_p1
    REAL(wp)            :: beta_soil_flush_p1
    REAL(wp)            :: beta_soil_senescence_p1
    REAL(wp)            :: gdd_t_air_threshold_p1
    REAL(wp)            :: gdd_req_max_p1
    REAL(wp)            :: k_gdd_dormance_p1
    REAL(wp)            :: t_air_senescence_p1
    REAL(wp)            :: min_leaf_age_p1
    REAL(wp)            :: frac_sapwood_branch_p1
    REAL(wp)            :: wood_density_p1
    REAL(wp)            :: k_latosa_p1
    REAL(wp)            :: k_crtos_p1
    REAL(wp)            :: k_rtos_p1
    REAL(wp)            :: k2_fruit_alloc_p1
    REAL(wp)            :: allom_k1_p1
    REAL(wp)            :: allom_k2_p1
    REAL(wp)            :: phi_leaf_min_p1
    REAL(wp)            :: k_root_p1
    REAL(wp)            :: k_sapwood_p1
    REAL(wp)            :: c0_allom_p1
    REAL(wp)            :: fstore_target_p1
    REAL(wp)            :: k_root_dist_p1
    ! PFT 2
    REAL(wp)            :: sla_p2
    REAL(wp)            :: sigma_vis_p2
    REAL(wp)            :: sigma_nir_p2
    REAL(wp)            :: omega_clumping_p2
    REAL(wp)            :: crown_shape_factor_p2
    REAL(wp)            :: cn_leaf_p2
    REAL(wp)            :: np_leaf_p2
    REAL(wp)            :: k0_fn_struc_p2
    REAL(wp)            :: fn_oth_min_p2
    REAL(wp)            :: t_jmax_opt_p2
    REAL(wp)            :: t_jmax_omega_p2
    REAL(wp)            :: g0_p2
    REAL(wp)            :: g1_p2
    REAL(wp)            :: gmin_p2
    REAL(wp)            :: lambda_soil_ps_p2
    REAL(wp)            :: weibull_k_soil_ps_p2
    REAL(wp)            :: lambda_soil_gs_p2
    REAL(wp)            :: weibull_k_soil_gs_p2
    REAL(wp)            :: tau_leaf_p2
    REAL(wp)            :: tau_fine_root_p2
    REAL(wp)            :: tau_coarse_root_p2
    REAL(wp)            :: tau_branch_p2
    REAL(wp)            :: tau_sap_wood_p2
    REAL(wp)            :: tau_fruit_p2
    REAL(wp)            :: tau_seed_litter_p2
    REAL(wp)            :: tau_seed_est_p2
    REAL(wp)            :: tau_mycorrhiza_p2
    REAL(wp)            :: vmax_uptake_n_p2
    REAL(wp)            :: vmax_uptake_p_p2
    REAL(wp)            :: bnf_base_p2
    REAL(wp)            :: lambda_est_light_p2
    REAL(wp)            :: k_est_light_p2
    REAL(wp)            :: seed_size_p2
    REAL(wp)            :: beta_soil_flush_p2
    REAL(wp)            :: beta_soil_senescence_p2
    REAL(wp)            :: gdd_t_air_threshold_p2
    REAL(wp)            :: gdd_req_max_p2
    REAL(wp)            :: k_gdd_dormance_p2
    REAL(wp)            :: t_air_senescence_p2
    REAL(wp)            :: min_leaf_age_p2
    REAL(wp)            :: frac_sapwood_branch_p2
    REAL(wp)            :: wood_density_p2
    REAL(wp)            :: k_latosa_p2
    REAL(wp)            :: k_crtos_p2
    REAL(wp)            :: k_rtos_p2
    REAL(wp)            :: k2_fruit_alloc_p2
    REAL(wp)            :: allom_k1_p2
    REAL(wp)            :: allom_k2_p2
    REAL(wp)            :: phi_leaf_min_p2
    REAL(wp)            :: k_root_p2
    REAL(wp)            :: k_sapwood_p2
    REAL(wp)            :: c0_allom_p2
    REAL(wp)            :: fstore_target_p2
    REAL(wp)            :: k_root_dist_p2
    ! PFT 3
    REAL(wp)            :: sla_p3
    REAL(wp)            :: sigma_vis_p3
    REAL(wp)            :: sigma_nir_p3
    REAL(wp)            :: omega_clumping_p3
    REAL(wp)            :: crown_shape_factor_p3
    REAL(wp)            :: cn_leaf_p3
    REAL(wp)            :: np_leaf_p3
    REAL(wp)            :: k0_fn_struc_p3
    REAL(wp)            :: fn_oth_min_p3
    REAL(wp)            :: t_jmax_opt_p3
    REAL(wp)            :: t_jmax_omega_p3
    REAL(wp)            :: g0_p3
    REAL(wp)            :: g1_p3
    REAL(wp)            :: gmin_p3
    REAL(wp)            :: lambda_soil_ps_p3
    REAL(wp)            :: weibull_k_soil_ps_p3
    REAL(wp)            :: lambda_soil_gs_p3
    REAL(wp)            :: weibull_k_soil_gs_p3
    REAL(wp)            :: tau_leaf_p3
    REAL(wp)            :: tau_fine_root_p3
    REAL(wp)            :: tau_coarse_root_p3
    REAL(wp)            :: tau_branch_p3
    REAL(wp)            :: tau_sap_wood_p3
    REAL(wp)            :: tau_fruit_p3
    REAL(wp)            :: tau_seed_litter_p3
    REAL(wp)            :: tau_seed_est_p3
    REAL(wp)            :: tau_mycorrhiza_p3
    REAL(wp)            :: vmax_uptake_n_p3
    REAL(wp)            :: vmax_uptake_p_p3
    REAL(wp)            :: bnf_base_p3
    REAL(wp)            :: lambda_est_light_p3
    REAL(wp)            :: k_est_light_p3
    REAL(wp)            :: seed_size_p3
    REAL(wp)            :: beta_soil_flush_p3
    REAL(wp)            :: beta_soil_senescence_p3
    REAL(wp)            :: gdd_t_air_threshold_p3
    REAL(wp)            :: gdd_req_max_p3
    REAL(wp)            :: k_gdd_dormance_p3
    REAL(wp)            :: t_air_senescence_p3
    REAL(wp)            :: min_leaf_age_p3
    REAL(wp)            :: frac_sapwood_branch_p3
    REAL(wp)            :: wood_density_p3
    REAL(wp)            :: k_latosa_p3
    REAL(wp)            :: k_crtos_p3
    REAL(wp)            :: k_rtos_p3
    REAL(wp)            :: k2_fruit_alloc_p3
    REAL(wp)            :: allom_k1_p3
    REAL(wp)            :: allom_k2_p3
    REAL(wp)            :: phi_leaf_min_p3
    REAL(wp)            :: k_root_p3
    REAL(wp)            :: k_sapwood_p3
    REAL(wp)            :: c0_allom_p3
    REAL(wp)            :: fstore_target_p3
    REAL(wp)            :: k_root_dist_p3
    ! PFT 4
    REAL(wp)            :: sla_p4
    REAL(wp)            :: sigma_vis_p4
    REAL(wp)            :: sigma_nir_p4
    REAL(wp)            :: omega_clumping_p4
    REAL(wp)            :: crown_shape_factor_p4
    REAL(wp)            :: cn_leaf_p4
    REAL(wp)            :: np_leaf_p4
    REAL(wp)            :: k0_fn_struc_p4
    REAL(wp)            :: fn_oth_min_p4
    REAL(wp)            :: t_jmax_opt_p4
    REAL(wp)            :: t_jmax_omega_p4
    REAL(wp)            :: g0_p4
    REAL(wp)            :: g1_p4
    REAL(wp)            :: gmin_p4
    REAL(wp)            :: lambda_soil_ps_p4
    REAL(wp)            :: weibull_k_soil_ps_p4
    REAL(wp)            :: lambda_soil_gs_p4
    REAL(wp)            :: weibull_k_soil_gs_p4
    REAL(wp)            :: tau_leaf_p4
    REAL(wp)            :: tau_fine_root_p4
    REAL(wp)            :: tau_coarse_root_p4
    REAL(wp)            :: tau_branch_p4
    REAL(wp)            :: tau_sap_wood_p4
    REAL(wp)            :: tau_fruit_p4
    REAL(wp)            :: tau_seed_litter_p4
    REAL(wp)            :: tau_seed_est_p4
    REAL(wp)            :: tau_mycorrhiza_p4
    REAL(wp)            :: vmax_uptake_n_p4
    REAL(wp)            :: vmax_uptake_p_p4
    REAL(wp)            :: bnf_base_p4
    REAL(wp)            :: lambda_est_light_p4
    REAL(wp)            :: k_est_light_p4
    REAL(wp)            :: seed_size_p4
    REAL(wp)            :: beta_soil_flush_p4
    REAL(wp)            :: beta_soil_senescence_p4
    REAL(wp)            :: gdd_t_air_threshold_p4
    REAL(wp)            :: gdd_req_max_p4
    REAL(wp)            :: k_gdd_dormance_p4
    REAL(wp)            :: t_air_senescence_p4
    REAL(wp)            :: min_leaf_age_p4
    REAL(wp)            :: frac_sapwood_branch_p4
    REAL(wp)            :: wood_density_p4
    REAL(wp)            :: k_latosa_p4
    REAL(wp)            :: k_crtos_p4
    REAL(wp)            :: k_rtos_p4
    REAL(wp)            :: k2_fruit_alloc_p4
    REAL(wp)            :: allom_k1_p4
    REAL(wp)            :: allom_k2_p4
    REAL(wp)            :: phi_leaf_min_p4
    REAL(wp)            :: k_root_p4
    REAL(wp)            :: k_sapwood_p4
    REAL(wp)            :: c0_allom_p4
    REAL(wp)            :: fstore_target_p4
    REAL(wp)            :: k_root_dist_p4
    ! PFT 5
    REAL(wp)            :: sla_p5
    REAL(wp)            :: sigma_vis_p5
    REAL(wp)            :: sigma_nir_p5
    REAL(wp)            :: omega_clumping_p5
    REAL(wp)            :: crown_shape_factor_p5
    REAL(wp)            :: cn_leaf_p5
    REAL(wp)            :: np_leaf_p5
    REAL(wp)            :: k0_fn_struc_p5
    REAL(wp)            :: fn_oth_min_p5
    REAL(wp)            :: t_jmax_opt_p5
    REAL(wp)            :: t_jmax_omega_p5
    REAL(wp)            :: g0_p5
    REAL(wp)            :: g1_p5
    REAL(wp)            :: gmin_p5
    REAL(wp)            :: lambda_soil_ps_p5
    REAL(wp)            :: weibull_k_soil_ps_p5
    REAL(wp)            :: lambda_soil_gs_p5
    REAL(wp)            :: weibull_k_soil_gs_p5
    REAL(wp)            :: tau_leaf_p5
    REAL(wp)            :: tau_fine_root_p5
    REAL(wp)            :: tau_coarse_root_p5
    REAL(wp)            :: tau_branch_p5
    REAL(wp)            :: tau_sap_wood_p5
    REAL(wp)            :: tau_fruit_p5
    REAL(wp)            :: tau_seed_litter_p5
    REAL(wp)            :: tau_seed_est_p5
    REAL(wp)            :: tau_mycorrhiza_p5
    REAL(wp)            :: vmax_uptake_n_p5
    REAL(wp)            :: vmax_uptake_p_p5
    REAL(wp)            :: bnf_base_p5
    REAL(wp)            :: lambda_est_light_p5
    REAL(wp)            :: k_est_light_p5
    REAL(wp)            :: seed_size_p5
    REAL(wp)            :: beta_soil_flush_p5
    REAL(wp)            :: beta_soil_senescence_p5
    REAL(wp)            :: gdd_t_air_threshold_p5
    REAL(wp)            :: gdd_req_max_p5
    REAL(wp)            :: k_gdd_dormance_p5
    REAL(wp)            :: t_air_senescence_p5
    REAL(wp)            :: min_leaf_age_p5
    REAL(wp)            :: frac_sapwood_branch_p5
    REAL(wp)            :: wood_density_p5
    REAL(wp)            :: k_latosa_p5
    REAL(wp)            :: k_crtos_p5
    REAL(wp)            :: k_rtos_p5
    REAL(wp)            :: k2_fruit_alloc_p5
    REAL(wp)            :: allom_k1_p5
    REAL(wp)            :: allom_k2_p5
    REAL(wp)            :: phi_leaf_min_p5
    REAL(wp)            :: k_root_p5
    REAL(wp)            :: k_sapwood_p5
    REAL(wp)            :: c0_allom_p5
    REAL(wp)            :: fstore_target_p5
    REAL(wp)            :: k_root_dist_p5
    ! PFT 6
    REAL(wp)            :: sla_p6
    REAL(wp)            :: sigma_vis_p6
    REAL(wp)            :: sigma_nir_p6
    REAL(wp)            :: omega_clumping_p6
    REAL(wp)            :: crown_shape_factor_p6
    REAL(wp)            :: cn_leaf_p6
    REAL(wp)            :: np_leaf_p6
    REAL(wp)            :: k0_fn_struc_p6
    REAL(wp)            :: fn_oth_min_p6
    REAL(wp)            :: t_jmax_opt_p6
    REAL(wp)            :: t_jmax_omega_p6
    REAL(wp)            :: g0_p6
    REAL(wp)            :: g1_p6
    REAL(wp)            :: gmin_p6
    REAL(wp)            :: lambda_soil_ps_p6
    REAL(wp)            :: weibull_k_soil_ps_p6
    REAL(wp)            :: lambda_soil_gs_p6
    REAL(wp)            :: weibull_k_soil_gs_p6
    REAL(wp)            :: tau_leaf_p6
    REAL(wp)            :: tau_fine_root_p6
    REAL(wp)            :: tau_coarse_root_p6
    REAL(wp)            :: tau_branch_p6
    REAL(wp)            :: tau_sap_wood_p6
    REAL(wp)            :: tau_fruit_p6
    REAL(wp)            :: tau_seed_litter_p6
    REAL(wp)            :: tau_seed_est_p6
    REAL(wp)            :: tau_mycorrhiza_p6
    REAL(wp)            :: vmax_uptake_n_p6
    REAL(wp)            :: vmax_uptake_p_p6
    REAL(wp)            :: bnf_base_p6
    REAL(wp)            :: lambda_est_light_p6
    REAL(wp)            :: k_est_light_p6
    REAL(wp)            :: seed_size_p6
    REAL(wp)            :: beta_soil_flush_p6
    REAL(wp)            :: beta_soil_senescence_p6
    REAL(wp)            :: gdd_t_air_threshold_p6
    REAL(wp)            :: gdd_req_max_p6
    REAL(wp)            :: k_gdd_dormance_p6
    REAL(wp)            :: t_air_senescence_p6
    REAL(wp)            :: min_leaf_age_p6
    REAL(wp)            :: frac_sapwood_branch_p6
    REAL(wp)            :: wood_density_p6
    REAL(wp)            :: k_latosa_p6
    REAL(wp)            :: k_crtos_p6
    REAL(wp)            :: k_rtos_p6
    REAL(wp)            :: k2_fruit_alloc_p6
    REAL(wp)            :: allom_k1_p6
    REAL(wp)            :: allom_k2_p6
    REAL(wp)            :: phi_leaf_min_p6
    REAL(wp)            :: k_root_p6
    REAL(wp)            :: k_sapwood_p6
    REAL(wp)            :: c0_allom_p6
    REAL(wp)            :: fstore_target_p6
    REAL(wp)            :: k_root_dist_p6
    ! PFT 7
    REAL(wp)            :: sla_p7
    REAL(wp)            :: sigma_vis_p7
    REAL(wp)            :: sigma_nir_p7
    REAL(wp)            :: omega_clumping_p7
    REAL(wp)            :: crown_shape_factor_p7
    REAL(wp)            :: cn_leaf_p7
    REAL(wp)            :: np_leaf_p7
    REAL(wp)            :: k0_fn_struc_p7
    REAL(wp)            :: fn_oth_min_p7
    REAL(wp)            :: t_jmax_opt_p7
    REAL(wp)            :: t_jmax_omega_p7
    REAL(wp)            :: g0_p7
    REAL(wp)            :: g1_p7
    REAL(wp)            :: gmin_p7
    REAL(wp)            :: lambda_soil_ps_p7
    REAL(wp)            :: weibull_k_soil_ps_p7
    REAL(wp)            :: lambda_soil_gs_p7
    REAL(wp)            :: weibull_k_soil_gs_p7
    REAL(wp)            :: tau_leaf_p7
    REAL(wp)            :: tau_fine_root_p7
    REAL(wp)            :: tau_coarse_root_p7
    REAL(wp)            :: tau_branch_p7
    REAL(wp)            :: tau_sap_wood_p7
    REAL(wp)            :: tau_fruit_p7
    REAL(wp)            :: tau_seed_litter_p7
    REAL(wp)            :: tau_seed_est_p7
    REAL(wp)            :: tau_mycorrhiza_p7
    REAL(wp)            :: vmax_uptake_n_p7
    REAL(wp)            :: vmax_uptake_p_p7
    REAL(wp)            :: bnf_base_p7
    REAL(wp)            :: lambda_est_light_p7
    REAL(wp)            :: k_est_light_p7
    REAL(wp)            :: seed_size_p7
    REAL(wp)            :: beta_soil_flush_p7
    REAL(wp)            :: beta_soil_senescence_p7
    REAL(wp)            :: gdd_t_air_threshold_p7
    REAL(wp)            :: gdd_req_max_p7
    REAL(wp)            :: k_gdd_dormance_p7
    REAL(wp)            :: t_air_senescence_p7
    REAL(wp)            :: min_leaf_age_p7
    REAL(wp)            :: frac_sapwood_branch_p7
    REAL(wp)            :: wood_density_p7
    REAL(wp)            :: k_latosa_p7
    REAL(wp)            :: k_crtos_p7
    REAL(wp)            :: k_rtos_p7
    REAL(wp)            :: k2_fruit_alloc_p7
    REAL(wp)            :: allom_k1_p7
    REAL(wp)            :: allom_k2_p7
    REAL(wp)            :: phi_leaf_min_p7
    REAL(wp)            :: k_root_p7
    REAL(wp)            :: k_sapwood_p7
    REAL(wp)            :: c0_allom_p7
    REAL(wp)            :: fstore_target_p7
    REAL(wp)            :: k_root_dist_p7
    ! PFT 8
    REAL(wp)            :: sla_p8
    REAL(wp)            :: sigma_vis_p8
    REAL(wp)            :: sigma_nir_p8
    REAL(wp)            :: omega_clumping_p8
    REAL(wp)            :: crown_shape_factor_p8
    REAL(wp)            :: cn_leaf_p8
    REAL(wp)            :: np_leaf_p8
    REAL(wp)            :: k0_fn_struc_p8
    REAL(wp)            :: fn_oth_min_p8
    REAL(wp)            :: t_jmax_opt_p8
    REAL(wp)            :: t_jmax_omega_p8
    REAL(wp)            :: g0_p8
    REAL(wp)            :: g1_p8
    REAL(wp)            :: gmin_p8
    REAL(wp)            :: lambda_soil_ps_p8
    REAL(wp)            :: weibull_k_soil_ps_p8
    REAL(wp)            :: lambda_soil_gs_p8
    REAL(wp)            :: weibull_k_soil_gs_p8
    REAL(wp)            :: tau_leaf_p8
    REAL(wp)            :: tau_fine_root_p8
    REAL(wp)            :: tau_coarse_root_p8
    REAL(wp)            :: tau_branch_p8
    REAL(wp)            :: tau_sap_wood_p8
    REAL(wp)            :: tau_fruit_p8
    REAL(wp)            :: tau_seed_litter_p8
    REAL(wp)            :: tau_seed_est_p8
    REAL(wp)            :: tau_mycorrhiza_p8
    REAL(wp)            :: vmax_uptake_n_p8
    REAL(wp)            :: vmax_uptake_p_p8
    REAL(wp)            :: bnf_base_p8
    REAL(wp)            :: lambda_est_light_p8
    REAL(wp)            :: k_est_light_p8
    REAL(wp)            :: seed_size_p8
    REAL(wp)            :: beta_soil_flush_p8
    REAL(wp)            :: beta_soil_senescence_p8
    REAL(wp)            :: gdd_t_air_threshold_p8
    REAL(wp)            :: gdd_req_max_p8
    REAL(wp)            :: k_gdd_dormance_p8
    REAL(wp)            :: t_air_senescence_p8
    REAL(wp)            :: min_leaf_age_p8
    REAL(wp)            :: frac_sapwood_branch_p8
    REAL(wp)            :: wood_density_p8
    REAL(wp)            :: k_latosa_p8
    REAL(wp)            :: k_crtos_p8
    REAL(wp)            :: k_rtos_p8
    REAL(wp)            :: k2_fruit_alloc_p8
    REAL(wp)            :: allom_k1_p8
    REAL(wp)            :: allom_k2_p8
    REAL(wp)            :: phi_leaf_min_p8
    REAL(wp)            :: k_root_p8
    REAL(wp)            :: k_sapwood_p8
    REAL(wp)            :: c0_allom_p8
    REAL(wp)            :: fstore_target_p8
    REAL(wp)            :: k_root_dist_p8


    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':read_parameter_lctlib_list'

    ! Namelist definitions for each PFT
    ! 
    NAMELIST /LCTLIB_PFT1_CTL/  &
          sla_p1, &
          sigma_vis_p1, &
          sigma_nir_p1, &
          omega_clumping_p1, &
          crown_shape_factor_p1, &
          cn_leaf_p1, &
          np_leaf_p1, &
          k0_fn_struc_p1, &
          fn_oth_min_p1, &
          t_jmax_opt_p1, &
          t_jmax_omega_p1, &
          g0_p1, &
          g1_p1, &
          gmin_p1, &
          lambda_soil_ps_p1, &
          weibull_k_soil_ps_p1, &
          lambda_soil_gs_p1, &
          weibull_k_soil_gs_p1, &
          tau_leaf_p1, &
          tau_fine_root_p1, &
          tau_coarse_root_p1, &
          tau_branch_p1, &
          tau_sap_wood_p1, &
          tau_fruit_p1, &
          tau_seed_litter_p1, &
          tau_seed_est_p1, &
          tau_mycorrhiza_p1, &
          vmax_uptake_n_p1, &
          vmax_uptake_p_p1, &
          bnf_base_p1, &
          lambda_est_light_p1, &
          k_est_light_p1, &
          seed_size_p1, &
          beta_soil_flush_p1, &
          beta_soil_senescence_p1, &
          gdd_t_air_threshold_p1, &
          gdd_req_max_p1, &
          k_gdd_dormance_p1, &
          t_air_senescence_p1, &
          min_leaf_age_p1, &
          frac_sapwood_branch_p1, &
          wood_density_p1, &
          k_latosa_p1, &
          k_crtos_p1, &
          k_rtos_p1, &
          k2_fruit_alloc_p1, &
          allom_k1_p1, &
          allom_k2_p1, &
          phi_leaf_min_p1, &
          k_root_p1, &
          k_sapwood_p1, &
          c0_allom_p1, &
          fstore_target_p1, &
          k_root_dist_p1

    NAMELIST /LCTLIB_PFT2_CTL/  &
          sla_p2, &
          sigma_vis_p2, &
          sigma_nir_p2, &
          omega_clumping_p2, &
          crown_shape_factor_p2, &
          cn_leaf_p2, &
          np_leaf_p2, &
          k0_fn_struc_p2, &
          fn_oth_min_p2, &
          t_jmax_opt_p2, &
          t_jmax_omega_p2, &
          g0_p2, &
          g1_p2, &
          gmin_p2, &
          lambda_soil_ps_p2, &
          weibull_k_soil_ps_p2, &
          lambda_soil_gs_p2, &
          weibull_k_soil_gs_p2, &
          tau_leaf_p2, &
          tau_fine_root_p2, &
          tau_coarse_root_p2, &
          tau_branch_p2, &
          tau_sap_wood_p2, &
          tau_fruit_p2, &
          tau_seed_litter_p2, &
          tau_seed_est_p2, &
          tau_mycorrhiza_p2, &
          vmax_uptake_n_p2, &
          vmax_uptake_p_p2, &
          bnf_base_p2, &
          lambda_est_light_p2, &
          k_est_light_p2, &
          seed_size_p2, &
          beta_soil_flush_p2, &
          beta_soil_senescence_p2, &
          gdd_t_air_threshold_p2, &
          gdd_req_max_p2, &
          k_gdd_dormance_p2, &
          t_air_senescence_p2, &
          min_leaf_age_p2, &
          frac_sapwood_branch_p2, &
          wood_density_p2, &
          k_latosa_p2, &
          k_crtos_p2, &
          k_rtos_p2, &
          k2_fruit_alloc_p2, &
          allom_k1_p2, &
          allom_k2_p2, &
          phi_leaf_min_p2, &
          k_root_p2, &
          k_sapwood_p2, &
          c0_allom_p2, &
          fstore_target_p2, &
          k_root_dist_p2

    NAMELIST /LCTLIB_PFT3_CTL/  &
          sla_p3, &
          sigma_vis_p3, &
          sigma_nir_p3, &
          omega_clumping_p3, &
          crown_shape_factor_p3, &
          cn_leaf_p3, &
          np_leaf_p3, &
          k0_fn_struc_p3, &
          fn_oth_min_p3, &
          t_jmax_opt_p3, &
          t_jmax_omega_p3, &
          g0_p3, &
          g1_p3, &
          gmin_p3, &
          lambda_soil_ps_p3, &
          weibull_k_soil_ps_p3, &
          lambda_soil_gs_p3, &
          weibull_k_soil_gs_p3, &
          tau_leaf_p3, &
          tau_fine_root_p3, &
          tau_coarse_root_p3, &
          tau_branch_p3, &
          tau_sap_wood_p3, &
          tau_fruit_p3, &
          tau_seed_litter_p3, &
          tau_seed_est_p3, &
          tau_mycorrhiza_p3, &
          vmax_uptake_n_p3, &
          vmax_uptake_p_p3, &
          bnf_base_p3, &
          lambda_est_light_p3, &
          k_est_light_p3, &
          seed_size_p3, &
          beta_soil_flush_p3, &
          beta_soil_senescence_p3, &
          gdd_t_air_threshold_p3, &
          gdd_req_max_p3, &
          k_gdd_dormance_p3, &
          t_air_senescence_p3, &
          min_leaf_age_p3, &
          frac_sapwood_branch_p3, &
          wood_density_p3, &
          k_latosa_p3, &
          k_crtos_p3, &
          k_rtos_p3, &
          k2_fruit_alloc_p3, &
          allom_k1_p3, &
          allom_k2_p3, &
          phi_leaf_min_p3, &
          k_root_p3, &
          k_sapwood_p3, &
          c0_allom_p3, &
          fstore_target_p3, &
          k_root_dist_p3

    NAMELIST /LCTLIB_PFT4_CTL/  &
          sla_p4, &
          sigma_vis_p4, &
          sigma_nir_p4, &
          omega_clumping_p4, &
          crown_shape_factor_p4, &
          cn_leaf_p4, &
          np_leaf_p4, &
          k0_fn_struc_p4, &
          fn_oth_min_p4, &
          t_jmax_opt_p4, &
          t_jmax_omega_p4, &
          g0_p4, &
          g1_p4, &
          gmin_p4, &
          lambda_soil_ps_p4, &
          weibull_k_soil_ps_p4, &
          lambda_soil_gs_p4, &
          weibull_k_soil_gs_p4, &
          tau_leaf_p4, &
          tau_fine_root_p4, &
          tau_coarse_root_p4, &
          tau_branch_p4, &
          tau_sap_wood_p4, &
          tau_fruit_p4, &
          tau_seed_litter_p4, &
          tau_seed_est_p4, &
          tau_mycorrhiza_p4, &
          vmax_uptake_n_p4, &
          vmax_uptake_p_p4, &
          bnf_base_p4, &
          lambda_est_light_p4, &
          k_est_light_p4, &
          seed_size_p4, &
          beta_soil_flush_p4, &
          beta_soil_senescence_p4, &
          gdd_t_air_threshold_p4, &
          gdd_req_max_p4, &
          k_gdd_dormance_p4, &
          t_air_senescence_p4, &
          min_leaf_age_p4, &
          frac_sapwood_branch_p4, &
          wood_density_p4, &
          k_latosa_p4, &
          k_crtos_p4, &
          k_rtos_p4, &
          k2_fruit_alloc_p4, &
          allom_k1_p4, &
          allom_k2_p4, &
          phi_leaf_min_p4, &
          k_root_p4, &
          k_sapwood_p4, &
          c0_allom_p4, &
          fstore_target_p4, &
          k_root_dist_p4

    NAMELIST /LCTLIB_PFT5_CTL/  &
          sla_p5, &
          sigma_vis_p5, &
          sigma_nir_p5, &
          omega_clumping_p5, &
          crown_shape_factor_p5, &
          cn_leaf_p5, &
          np_leaf_p5, &
          k0_fn_struc_p5, &
          fn_oth_min_p5, &
          t_jmax_opt_p5, &
          t_jmax_omega_p5, &
          g0_p5, &
          g1_p5, &
          gmin_p5, &
          lambda_soil_ps_p5, &
          weibull_k_soil_ps_p5, &
          lambda_soil_gs_p5, &
          weibull_k_soil_gs_p5, &
          tau_leaf_p5, &
          tau_fine_root_p5, &
          tau_coarse_root_p5, &
          tau_branch_p5, &
          tau_sap_wood_p5, &
          tau_fruit_p5, &
          tau_seed_litter_p5, &
          tau_seed_est_p5, &
          tau_mycorrhiza_p5, &
          vmax_uptake_n_p5, &
          vmax_uptake_p_p5, &
          bnf_base_p5, &
          lambda_est_light_p5, &
          k_est_light_p5, &
          seed_size_p5, &
          beta_soil_flush_p5, &
          beta_soil_senescence_p5, &
          gdd_t_air_threshold_p5, &
          gdd_req_max_p5, &
          k_gdd_dormance_p5, &
          t_air_senescence_p5, &
          min_leaf_age_p5, &
          frac_sapwood_branch_p5, &
          wood_density_p5, &
          k_latosa_p5, &
          k_crtos_p5, &
          k_rtos_p5, &
          k2_fruit_alloc_p5, &
          allom_k1_p5, &
          allom_k2_p5, &
          phi_leaf_min_p5, &
          k_root_p5, &
          k_sapwood_p5, &
          c0_allom_p5, &
          fstore_target_p5, &
          k_root_dist_p5

    NAMELIST /LCTLIB_PFT6_CTL/  &
          sla_p6, &
          sigma_vis_p6, &
          sigma_nir_p6, &
          omega_clumping_p6, &
          crown_shape_factor_p6, &
          cn_leaf_p6, &
          np_leaf_p6, &
          k0_fn_struc_p6, &
          fn_oth_min_p6, &
          t_jmax_opt_p6, &
          t_jmax_omega_p6, &
          g0_p6, &
          g1_p6, &
          gmin_p6, &
          lambda_soil_ps_p6, &
          weibull_k_soil_ps_p6, &
          lambda_soil_gs_p6, &
          weibull_k_soil_gs_p6, &
          tau_leaf_p6, &
          tau_fine_root_p6, &
          tau_coarse_root_p6, &
          tau_branch_p6, &
          tau_sap_wood_p6, &
          tau_fruit_p6, &
          tau_seed_litter_p6, &
          tau_seed_est_p6, &
          tau_mycorrhiza_p6, &
          vmax_uptake_n_p6, &
          vmax_uptake_p_p6, &
          bnf_base_p6, &
          lambda_est_light_p6, &
          k_est_light_p6, &
          seed_size_p6, &
          beta_soil_flush_p6, &
          beta_soil_senescence_p6, &
          gdd_t_air_threshold_p6, &
          gdd_req_max_p6, &
          k_gdd_dormance_p6, &
          t_air_senescence_p6, &
          min_leaf_age_p6, &
          frac_sapwood_branch_p6, &
          wood_density_p6, &
          k_latosa_p6, &
          k_crtos_p6, &
          k_rtos_p6, &
          k2_fruit_alloc_p6, &
          allom_k1_p6, &
          allom_k2_p6, &
          phi_leaf_min_p6, &
          k_root_p6, &
          k_sapwood_p6, &
          c0_allom_p6, &
          fstore_target_p6, &
          k_root_dist_p6

    NAMELIST /LCTLIB_PFT7_CTL/  &
          sla_p7, &
          sigma_vis_p7, &
          sigma_nir_p7, &
          omega_clumping_p7, &
          crown_shape_factor_p7, &
          cn_leaf_p7, &
          np_leaf_p7, &
          k0_fn_struc_p7, &
          fn_oth_min_p7, &
          t_jmax_opt_p7, &
          t_jmax_omega_p7, &
          g0_p7, &
          g1_p7, &
          gmin_p7, &
          lambda_soil_ps_p7, &
          weibull_k_soil_ps_p7, &
          lambda_soil_gs_p7, &
          weibull_k_soil_gs_p7, &
          tau_leaf_p7, &
          tau_fine_root_p7, &
          tau_coarse_root_p7, &
          tau_branch_p7, &
          tau_sap_wood_p7, &
          tau_fruit_p7, &
          tau_seed_litter_p7, &
          tau_seed_est_p7, &
          tau_mycorrhiza_p7, &
          vmax_uptake_n_p7, &
          vmax_uptake_p_p7, &
          bnf_base_p7, &
          lambda_est_light_p7, &
          k_est_light_p7, &
          seed_size_p7, &
          beta_soil_flush_p7, &
          beta_soil_senescence_p7, &
          gdd_t_air_threshold_p7, &
          gdd_req_max_p7, &
          k_gdd_dormance_p7, &
          t_air_senescence_p7, &
          min_leaf_age_p7, &
          frac_sapwood_branch_p7, &
          wood_density_p7, &
          k_latosa_p7, &
          k_crtos_p7, &
          k_rtos_p7, &
          k2_fruit_alloc_p7, &
          allom_k1_p7, &
          allom_k2_p7, &
          phi_leaf_min_p7, &
          k_root_p7, &
          k_sapwood_p7, &
          c0_allom_p7, &
          fstore_target_p7, &
          k_root_dist_p7

    NAMELIST /LCTLIB_PFT8_CTL/  &
          sla_p8, &
          sigma_vis_p8, &
          sigma_nir_p8, &
          omega_clumping_p8, &
          crown_shape_factor_p8, &
          cn_leaf_p8, &
          np_leaf_p8, &
          k0_fn_struc_p8, &
          fn_oth_min_p8, &
          t_jmax_opt_p8, &
          t_jmax_omega_p8, &
          g0_p8, &
          g1_p8, &
          gmin_p8, &
          lambda_soil_ps_p8, &
          weibull_k_soil_ps_p8, &
          lambda_soil_gs_p8, &
          weibull_k_soil_gs_p8, &
          tau_leaf_p8, &
          tau_fine_root_p8, &
          tau_coarse_root_p8, &
          tau_branch_p8, &
          tau_sap_wood_p8, &
          tau_fruit_p8, &
          tau_seed_litter_p8, &
          tau_seed_est_p8, &
          tau_mycorrhiza_p8, &
          vmax_uptake_n_p8, &
          vmax_uptake_p_p8, &
          bnf_base_p8, &
          lambda_est_light_p8, &
          k_est_light_p8, &
          seed_size_p8, &
          beta_soil_flush_p8, &
          beta_soil_senescence_p8, &
          gdd_t_air_threshold_p8, &
          gdd_req_max_p8, &
          k_gdd_dormance_p8, &
          t_air_senescence_p8, &
          min_leaf_age_p8, &
          frac_sapwood_branch_p8, &
          wood_density_p8, &
          k_latosa_p8, &
          k_crtos_p8, &
          k_rtos_p8, &
          k2_fruit_alloc_p8, &
          allom_k1_p8, &
          allom_k2_p8, &
          phi_leaf_min_p8, &
          k_root_p8, &
          k_sapwood_p8, &
          c0_allom_p8, &
          fstore_target_p8, &
          k_root_dist_p8


    ! variables for reading model-options from namelist
    ! 
    INTEGER  :: nml_handler, nml_unit, istat

    CALL message(TRIM(routine), 'Starting to read parameter list')

    ! Set defaults
    ! 
    ! PFT1
    sla_p1      = missval
    sigma_vis_p1      = missval
    sigma_nir_p1      = missval
    omega_clumping_p1      = missval
    crown_shape_factor_p1      = missval
    cn_leaf_p1      = missval
    np_leaf_p1      = missval
    k0_fn_struc_p1      = missval
    fn_oth_min_p1      = missval
    t_jmax_opt_p1      = missval
    t_jmax_omega_p1      = missval
    g0_p1      = missval
    g1_p1      = missval
    gmin_p1      = missval
    lambda_soil_ps_p1      = missval
    weibull_k_soil_ps_p1      = missval
    lambda_soil_gs_p1      = missval
    weibull_k_soil_gs_p1      = missval
    tau_leaf_p1      = missval
    tau_fine_root_p1      = missval
    tau_coarse_root_p1      = missval
    tau_branch_p1      = missval
    tau_sap_wood_p1      = missval
    tau_fruit_p1      = missval
    tau_seed_litter_p1      = missval
    tau_seed_est_p1      = missval
    tau_mycorrhiza_p1      = missval
    vmax_uptake_n_p1      = missval
    vmax_uptake_p_p1      = missval
    bnf_base_p1      = missval
    lambda_est_light_p1      = missval
    k_est_light_p1      = missval
    seed_size_p1      = missval
    beta_soil_flush_p1      = missval
    beta_soil_senescence_p1      = missval
    gdd_t_air_threshold_p1      = missval
    gdd_req_max_p1      = missval
    k_gdd_dormance_p1      = missval
    t_air_senescence_p1      = missval
    min_leaf_age_p1      = missval
    frac_sapwood_branch_p1      = missval
    wood_density_p1      = missval
    k_latosa_p1      = missval
    k_crtos_p1      = missval
    k_rtos_p1      = missval
    k2_fruit_alloc_p1      = missval
    allom_k1_p1      = missval
    allom_k2_p1      = missval
    phi_leaf_min_p1      = missval
    k_root_p1      = missval
    k_sapwood_p1      = missval
    c0_allom_p1      = missval
    fstore_target_p1      = missval
    k_root_dist_p1      = missval
    ! PFT2
    sla_p2      = missval
    sigma_vis_p2      = missval
    sigma_nir_p2      = missval
    omega_clumping_p2      = missval
    crown_shape_factor_p2      = missval
    cn_leaf_p2      = missval
    np_leaf_p2      = missval
    k0_fn_struc_p2      = missval
    fn_oth_min_p2      = missval
    t_jmax_opt_p2      = missval
    t_jmax_omega_p2      = missval
    g0_p2      = missval
    g1_p2      = missval
    gmin_p2      = missval
    lambda_soil_ps_p2      = missval
    weibull_k_soil_ps_p2      = missval
    lambda_soil_gs_p2      = missval
    weibull_k_soil_gs_p2      = missval
    tau_leaf_p2      = missval
    tau_fine_root_p2      = missval
    tau_coarse_root_p2      = missval
    tau_branch_p2      = missval
    tau_sap_wood_p2      = missval
    tau_fruit_p2      = missval
    tau_seed_litter_p2      = missval
    tau_seed_est_p2      = missval
    tau_mycorrhiza_p2      = missval
    vmax_uptake_n_p2      = missval
    vmax_uptake_p_p2      = missval
    bnf_base_p2      = missval
    lambda_est_light_p2      = missval
    k_est_light_p2      = missval
    seed_size_p2      = missval
    beta_soil_flush_p2      = missval
    beta_soil_senescence_p2      = missval
    gdd_t_air_threshold_p2      = missval
    gdd_req_max_p2      = missval
    k_gdd_dormance_p2      = missval
    t_air_senescence_p2      = missval
    min_leaf_age_p2      = missval
    frac_sapwood_branch_p2      = missval
    wood_density_p2      = missval
    k_latosa_p2      = missval
    k_crtos_p2      = missval
    k_rtos_p2      = missval
    k2_fruit_alloc_p2      = missval
    allom_k1_p2      = missval
    allom_k2_p2      = missval
    phi_leaf_min_p2      = missval
    k_root_p2      = missval
    k_sapwood_p2      = missval
    c0_allom_p2      = missval
    fstore_target_p2      = missval
    k_root_dist_p2      = missval
    ! PFT3
    sla_p3      = missval
    sigma_vis_p3      = missval
    sigma_nir_p3      = missval
    omega_clumping_p3      = missval
    crown_shape_factor_p3      = missval
    cn_leaf_p3      = missval
    np_leaf_p3      = missval
    k0_fn_struc_p3      = missval
    fn_oth_min_p3      = missval
    t_jmax_opt_p3      = missval
    t_jmax_omega_p3      = missval
    g0_p3      = missval
    g1_p3      = missval
    gmin_p3      = missval
    lambda_soil_ps_p3      = missval
    weibull_k_soil_ps_p3      = missval
    lambda_soil_gs_p3      = missval
    weibull_k_soil_gs_p3      = missval
    tau_leaf_p3      = missval
    tau_fine_root_p3      = missval
    tau_coarse_root_p3      = missval
    tau_branch_p3      = missval
    tau_sap_wood_p3      = missval
    tau_fruit_p3      = missval
    tau_seed_litter_p3      = missval
    tau_seed_est_p3      = missval
    tau_mycorrhiza_p3      = missval
    vmax_uptake_n_p3      = missval
    vmax_uptake_p_p3      = missval
    bnf_base_p3      = missval
    lambda_est_light_p3      = missval
    k_est_light_p3      = missval
    seed_size_p3      = missval
    beta_soil_flush_p3      = missval
    beta_soil_senescence_p3      = missval
    gdd_t_air_threshold_p3      = missval
    gdd_req_max_p3      = missval
    k_gdd_dormance_p3      = missval
    t_air_senescence_p3      = missval
    min_leaf_age_p3      = missval
    frac_sapwood_branch_p3      = missval
    wood_density_p3      = missval
    k_latosa_p3      = missval
    k_crtos_p3      = missval
    k_rtos_p3      = missval
    k2_fruit_alloc_p3      = missval
    allom_k1_p3      = missval
    allom_k2_p3      = missval
    phi_leaf_min_p3      = missval
    k_root_p3      = missval
    k_sapwood_p3      = missval
    c0_allom_p3      = missval
    fstore_target_p3      = missval
    k_root_dist_p3      = missval
    ! PFT4
    sla_p4      = missval
    sigma_vis_p4      = missval
    sigma_nir_p4      = missval
    omega_clumping_p4      = missval
    crown_shape_factor_p4      = missval
    cn_leaf_p4      = missval
    np_leaf_p4      = missval
    k0_fn_struc_p4      = missval
    fn_oth_min_p4      = missval
    t_jmax_opt_p4      = missval
    t_jmax_omega_p4      = missval
    g0_p4      = missval
    g1_p4      = missval
    gmin_p4      = missval
    lambda_soil_ps_p4      = missval
    weibull_k_soil_ps_p4      = missval
    lambda_soil_gs_p4      = missval
    weibull_k_soil_gs_p4      = missval
    tau_leaf_p4      = missval
    tau_fine_root_p4      = missval
    tau_coarse_root_p4      = missval
    tau_branch_p4      = missval
    tau_sap_wood_p4      = missval
    tau_fruit_p4      = missval
    tau_seed_litter_p4      = missval
    tau_seed_est_p4      = missval
    tau_mycorrhiza_p4      = missval
    vmax_uptake_n_p4      = missval
    vmax_uptake_p_p4      = missval
    bnf_base_p4      = missval
    lambda_est_light_p4      = missval
    k_est_light_p4      = missval
    seed_size_p4      = missval
    beta_soil_flush_p4      = missval
    beta_soil_senescence_p4      = missval
    gdd_t_air_threshold_p4      = missval
    gdd_req_max_p4      = missval
    k_gdd_dormance_p4      = missval
    t_air_senescence_p4      = missval
    min_leaf_age_p4      = missval
    frac_sapwood_branch_p4      = missval
    wood_density_p4      = missval
    k_latosa_p4      = missval
    k_crtos_p4      = missval
    k_rtos_p4      = missval
    k2_fruit_alloc_p4      = missval
    allom_k1_p4      = missval
    allom_k2_p4      = missval
    phi_leaf_min_p4      = missval
    k_root_p4      = missval
    k_sapwood_p4      = missval
    c0_allom_p4      = missval
    fstore_target_p4      = missval
    k_root_dist_p4      = missval
    ! PFT5
    sla_p5      = missval
    sigma_vis_p5      = missval
    sigma_nir_p5      = missval
    omega_clumping_p5      = missval
    crown_shape_factor_p5      = missval
    cn_leaf_p5      = missval
    np_leaf_p5      = missval
    k0_fn_struc_p5      = missval
    fn_oth_min_p5      = missval
    t_jmax_opt_p5      = missval
    t_jmax_omega_p5      = missval
    g0_p5      = missval
    g1_p5      = missval
    gmin_p5      = missval
    lambda_soil_ps_p5      = missval
    weibull_k_soil_ps_p5      = missval
    lambda_soil_gs_p5      = missval
    weibull_k_soil_gs_p5      = missval
    tau_leaf_p5      = missval
    tau_fine_root_p5      = missval
    tau_coarse_root_p5      = missval
    tau_branch_p5      = missval
    tau_sap_wood_p5      = missval
    tau_fruit_p5      = missval
    tau_seed_litter_p5      = missval
    tau_seed_est_p5      = missval
    tau_mycorrhiza_p5      = missval
    vmax_uptake_n_p5      = missval
    vmax_uptake_p_p5      = missval
    bnf_base_p5      = missval
    lambda_est_light_p5      = missval
    k_est_light_p5      = missval
    seed_size_p5      = missval
    beta_soil_flush_p5      = missval
    beta_soil_senescence_p5      = missval
    gdd_t_air_threshold_p5      = missval
    gdd_req_max_p5      = missval
    k_gdd_dormance_p5      = missval
    t_air_senescence_p5      = missval
    min_leaf_age_p5      = missval
    frac_sapwood_branch_p5      = missval
    wood_density_p5      = missval
    k_latosa_p5      = missval
    k_crtos_p5      = missval
    k_rtos_p5      = missval
    k2_fruit_alloc_p5      = missval
    allom_k1_p5      = missval
    allom_k2_p5      = missval
    phi_leaf_min_p5      = missval
    k_root_p5      = missval
    k_sapwood_p5      = missval
    c0_allom_p5      = missval
    fstore_target_p5      = missval
    k_root_dist_p5      = missval
    ! PFT6
    sla_p6      = missval
    sigma_vis_p6      = missval
    sigma_nir_p6      = missval
    omega_clumping_p6      = missval
    crown_shape_factor_p6      = missval
    cn_leaf_p6      = missval
    np_leaf_p6      = missval
    k0_fn_struc_p6      = missval
    fn_oth_min_p6      = missval
    t_jmax_opt_p6      = missval
    t_jmax_omega_p6      = missval
    g0_p6      = missval
    g1_p6      = missval
    gmin_p6      = missval
    lambda_soil_ps_p6      = missval
    weibull_k_soil_ps_p6      = missval
    lambda_soil_gs_p6      = missval
    weibull_k_soil_gs_p6      = missval
    tau_leaf_p6      = missval
    tau_fine_root_p6      = missval
    tau_coarse_root_p6      = missval
    tau_branch_p6      = missval
    tau_sap_wood_p6      = missval
    tau_fruit_p6      = missval
    tau_seed_litter_p6      = missval
    tau_seed_est_p6      = missval
    tau_mycorrhiza_p6      = missval
    vmax_uptake_n_p6      = missval
    vmax_uptake_p_p6      = missval
    bnf_base_p6      = missval
    lambda_est_light_p6      = missval
    k_est_light_p6      = missval
    seed_size_p6      = missval
    beta_soil_flush_p6      = missval
    beta_soil_senescence_p6      = missval
    gdd_t_air_threshold_p6      = missval
    gdd_req_max_p6      = missval
    k_gdd_dormance_p6      = missval
    t_air_senescence_p6      = missval
    min_leaf_age_p6      = missval
    frac_sapwood_branch_p6      = missval
    wood_density_p6      = missval
    k_latosa_p6      = missval
    k_crtos_p6      = missval
    k_rtos_p6      = missval
    k2_fruit_alloc_p6      = missval
    allom_k1_p6      = missval
    allom_k2_p6      = missval
    phi_leaf_min_p6      = missval
    k_root_p6      = missval
    k_sapwood_p6      = missval
    c0_allom_p6      = missval
    fstore_target_p6      = missval
    k_root_dist_p6      = missval
    ! PFT7
    sla_p7      = missval
    sigma_vis_p7      = missval
    sigma_nir_p7      = missval
    omega_clumping_p7      = missval
    crown_shape_factor_p7      = missval
    cn_leaf_p7      = missval
    np_leaf_p7      = missval
    k0_fn_struc_p7      = missval
    fn_oth_min_p7      = missval
    t_jmax_opt_p7      = missval
    t_jmax_omega_p7      = missval
    g0_p7      = missval
    g1_p7      = missval
    gmin_p7      = missval
    lambda_soil_ps_p7      = missval
    weibull_k_soil_ps_p7      = missval
    lambda_soil_gs_p7      = missval
    weibull_k_soil_gs_p7      = missval
    tau_leaf_p7      = missval
    tau_fine_root_p7      = missval
    tau_coarse_root_p7      = missval
    tau_branch_p7      = missval
    tau_sap_wood_p7      = missval
    tau_fruit_p7      = missval
    tau_seed_litter_p7      = missval
    tau_seed_est_p7      = missval
    tau_mycorrhiza_p7      = missval
    vmax_uptake_n_p7      = missval
    vmax_uptake_p_p7      = missval
    bnf_base_p7      = missval
    lambda_est_light_p7      = missval
    k_est_light_p7      = missval
    seed_size_p7      = missval
    beta_soil_flush_p7      = missval
    beta_soil_senescence_p7      = missval
    gdd_t_air_threshold_p7      = missval
    gdd_req_max_p7      = missval
    k_gdd_dormance_p7      = missval
    t_air_senescence_p7      = missval
    min_leaf_age_p7      = missval
    frac_sapwood_branch_p7      = missval
    wood_density_p7      = missval
    k_latosa_p7      = missval
    k_crtos_p7      = missval
    k_rtos_p7      = missval
    k2_fruit_alloc_p7      = missval
    allom_k1_p7      = missval
    allom_k2_p7      = missval
    phi_leaf_min_p7      = missval
    k_root_p7      = missval
    k_sapwood_p7      = missval
    c0_allom_p7      = missval
    fstore_target_p7      = missval
    k_root_dist_p7      = missval
    ! PFT8
    sla_p8      = missval
    sigma_vis_p8      = missval
    sigma_nir_p8      = missval
    omega_clumping_p8      = missval
    crown_shape_factor_p8      = missval
    cn_leaf_p8      = missval
    np_leaf_p8      = missval
    k0_fn_struc_p8      = missval
    fn_oth_min_p8      = missval
    t_jmax_opt_p8      = missval
    t_jmax_omega_p8      = missval
    g0_p8      = missval
    g1_p8      = missval
    gmin_p8      = missval
    lambda_soil_ps_p8      = missval
    weibull_k_soil_ps_p8      = missval
    lambda_soil_gs_p8      = missval
    weibull_k_soil_gs_p8      = missval
    tau_leaf_p8      = missval
    tau_fine_root_p8      = missval
    tau_coarse_root_p8      = missval
    tau_branch_p8      = missval
    tau_sap_wood_p8      = missval
    tau_fruit_p8      = missval
    tau_seed_litter_p8      = missval
    tau_seed_est_p8      = missval
    tau_mycorrhiza_p8      = missval
    vmax_uptake_n_p8      = missval
    vmax_uptake_p_p8      = missval
    bnf_base_p8      = missval
    lambda_est_light_p8      = missval
    k_est_light_p8      = missval
    seed_size_p8      = missval
    beta_soil_flush_p8      = missval
    beta_soil_senescence_p8      = missval
    gdd_t_air_threshold_p8      = missval
    gdd_req_max_p8      = missval
    k_gdd_dormance_p8      = missval
    t_air_senescence_p8      = missval
    min_leaf_age_p8      = missval
    frac_sapwood_branch_p8      = missval
    wood_density_p8      = missval
    k_latosa_p8      = missval
    k_crtos_p8      = missval
    k_rtos_p8      = missval
    k2_fruit_alloc_p8      = missval
    allom_k1_p8      = missval
    allom_k2_p8      = missval
    phi_leaf_min_p8      = missval
    k_root_p8      = missval
    k_sapwood_p8      = missval
    c0_allom_p8      = missval
    fstore_target_p8      = missval
    k_root_dist_p8      = missval


    ! read the namelist (parameters of all PFT are stored in the same namelist file using param_pft1_from_file%filename here)
    ! 
    ! &LCTLIB_PFT1_CTL
    nml_handler = open_nml(TRIM(param_pft1_from_file%filename))
    nml_unit = position_nml('LCTLIB_PFT1_CTL', nml_handler, STATUS=istat)
    IF (istat == POSITIONED) READ(nml_unit, LCTLIB_PFT1_CTL)
    CALL close_nml(nml_handler)
    ! &LCTLIB_PFT2_CTL
    nml_handler = open_nml(TRIM(param_pft1_from_file%filename))
    nml_unit = position_nml('LCTLIB_PFT2_CTL', nml_handler, STATUS=istat)
    IF (istat == POSITIONED) READ(nml_unit, LCTLIB_PFT2_CTL)
    CALL close_nml(nml_handler)
    ! &LCTLIB_PFT3_CTL
    nml_handler = open_nml(TRIM(param_pft1_from_file%filename))
    nml_unit = position_nml('LCTLIB_PFT3_CTL', nml_handler, STATUS=istat)
    IF (istat == POSITIONED) READ(nml_unit, LCTLIB_PFT3_CTL)
    CALL close_nml(nml_handler)
    ! &LCTLIB_PFT4_CTL
    nml_handler = open_nml(TRIM(param_pft1_from_file%filename))
    nml_unit = position_nml('LCTLIB_PFT4_CTL', nml_handler, STATUS=istat)
    IF (istat == POSITIONED) READ(nml_unit, LCTLIB_PFT4_CTL)
    CALL close_nml(nml_handler)
    ! &LCTLIB_PFT5_CTL
    nml_handler = open_nml(TRIM(param_pft1_from_file%filename))
    nml_unit = position_nml('LCTLIB_PFT5_CTL', nml_handler, STATUS=istat)
    IF (istat == POSITIONED) READ(nml_unit, LCTLIB_PFT5_CTL)
    CALL close_nml(nml_handler)
    ! &LCTLIB_PFT6_CTL
    nml_handler = open_nml(TRIM(param_pft1_from_file%filename))
    nml_unit = position_nml('LCTLIB_PFT6_CTL', nml_handler, STATUS=istat)
    IF (istat == POSITIONED) READ(nml_unit, LCTLIB_PFT6_CTL)
    CALL close_nml(nml_handler)
    ! &LCTLIB_PFT7_CTL
    nml_handler = open_nml(TRIM(param_pft1_from_file%filename))
    nml_unit = position_nml('LCTLIB_PFT7_CTL', nml_handler, STATUS=istat)
    IF (istat == POSITIONED) READ(nml_unit, LCTLIB_PFT7_CTL)
    CALL close_nml(nml_handler)
    ! &LCTLIB_PFT8_CTL
    nml_handler = open_nml(TRIM(param_pft1_from_file%filename))
    nml_unit = position_nml('LCTLIB_PFT8_CTL', nml_handler, STATUS=istat)
    IF (istat == POSITIONED) READ(nml_unit, LCTLIB_PFT8_CTL)
    CALL close_nml(nml_handler)


    ! pass values as read from file
    ! 
    ! LCTLIB_PFT1_CTL
    param_pft1_from_file%sla   = sla_p1
    param_pft1_from_file%sigma_vis   = sigma_vis_p1
    param_pft1_from_file%sigma_nir   = sigma_nir_p1
    param_pft1_from_file%omega_clumping   = omega_clumping_p1
    param_pft1_from_file%crown_shape_factor   = crown_shape_factor_p1
    param_pft1_from_file%cn_leaf   = cn_leaf_p1
    param_pft1_from_file%np_leaf   = np_leaf_p1
    param_pft1_from_file%k0_fn_struc   = k0_fn_struc_p1
    param_pft1_from_file%fn_oth_min   = fn_oth_min_p1
    param_pft1_from_file%t_jmax_opt   = t_jmax_opt_p1
    param_pft1_from_file%t_jmax_omega   = t_jmax_omega_p1
    param_pft1_from_file%g0   = g0_p1
    param_pft1_from_file%g1   = g1_p1
    param_pft1_from_file%gmin   = gmin_p1
    param_pft1_from_file%lambda_soil_ps   = lambda_soil_ps_p1
    param_pft1_from_file%weibull_k_soil_ps   = weibull_k_soil_ps_p1
    param_pft1_from_file%lambda_soil_gs   = lambda_soil_gs_p1
    param_pft1_from_file%weibull_k_soil_gs   = weibull_k_soil_gs_p1
    param_pft1_from_file%tau_leaf   = tau_leaf_p1
    param_pft1_from_file%tau_fine_root   = tau_fine_root_p1
    param_pft1_from_file%tau_coarse_root   = tau_coarse_root_p1
    param_pft1_from_file%tau_branch   = tau_branch_p1
    param_pft1_from_file%tau_sap_wood   = tau_sap_wood_p1
    param_pft1_from_file%tau_fruit   = tau_fruit_p1
    param_pft1_from_file%tau_seed_litter   = tau_seed_litter_p1
    param_pft1_from_file%tau_seed_est   = tau_seed_est_p1
    param_pft1_from_file%tau_mycorrhiza   = tau_mycorrhiza_p1
    param_pft1_from_file%vmax_uptake_n   = vmax_uptake_n_p1
    param_pft1_from_file%vmax_uptake_p   = vmax_uptake_p_p1
    param_pft1_from_file%bnf_base   = bnf_base_p1
    param_pft1_from_file%lambda_est_light   = lambda_est_light_p1
    param_pft1_from_file%k_est_light   = k_est_light_p1
    param_pft1_from_file%seed_size   = seed_size_p1
    param_pft1_from_file%beta_soil_flush   = beta_soil_flush_p1
    param_pft1_from_file%beta_soil_senescence   = beta_soil_senescence_p1
    param_pft1_from_file%gdd_t_air_threshold   = gdd_t_air_threshold_p1
    param_pft1_from_file%gdd_req_max   = gdd_req_max_p1
    param_pft1_from_file%k_gdd_dormance   = k_gdd_dormance_p1
    param_pft1_from_file%t_air_senescence   = t_air_senescence_p1
    param_pft1_from_file%min_leaf_age   = min_leaf_age_p1
    param_pft1_from_file%frac_sapwood_branch   = frac_sapwood_branch_p1
    param_pft1_from_file%wood_density   = wood_density_p1
    param_pft1_from_file%k_latosa   = k_latosa_p1
    param_pft1_from_file%k_crtos   = k_crtos_p1
    param_pft1_from_file%k_rtos   = k_rtos_p1
    param_pft1_from_file%k2_fruit_alloc   = k2_fruit_alloc_p1
    param_pft1_from_file%allom_k1   = allom_k1_p1
    param_pft1_from_file%allom_k2   = allom_k2_p1
    param_pft1_from_file%phi_leaf_min   = phi_leaf_min_p1
    param_pft1_from_file%k_root   = k_root_p1
    param_pft1_from_file%k_sapwood   = k_sapwood_p1
    param_pft1_from_file%c0_allom   = c0_allom_p1
    param_pft1_from_file%fstore_target   = fstore_target_p1
    param_pft1_from_file%k_root_dist   = k_root_dist_p1
    ! LCTLIB_PFT2_CTL
    param_pft2_from_file%sla   = sla_p2
    param_pft2_from_file%sigma_vis   = sigma_vis_p2
    param_pft2_from_file%sigma_nir   = sigma_nir_p2
    param_pft2_from_file%omega_clumping   = omega_clumping_p2
    param_pft2_from_file%crown_shape_factor   = crown_shape_factor_p2
    param_pft2_from_file%cn_leaf   = cn_leaf_p2
    param_pft2_from_file%np_leaf   = np_leaf_p2
    param_pft2_from_file%k0_fn_struc   = k0_fn_struc_p2
    param_pft2_from_file%fn_oth_min   = fn_oth_min_p2
    param_pft2_from_file%t_jmax_opt   = t_jmax_opt_p2
    param_pft2_from_file%t_jmax_omega   = t_jmax_omega_p2
    param_pft2_from_file%g0   = g0_p2
    param_pft2_from_file%g1   = g1_p2
    param_pft2_from_file%gmin   = gmin_p2
    param_pft2_from_file%lambda_soil_ps   = lambda_soil_ps_p2
    param_pft2_from_file%weibull_k_soil_ps   = weibull_k_soil_ps_p2
    param_pft2_from_file%lambda_soil_gs   = lambda_soil_gs_p2
    param_pft2_from_file%weibull_k_soil_gs   = weibull_k_soil_gs_p2
    param_pft2_from_file%tau_leaf   = tau_leaf_p2
    param_pft2_from_file%tau_fine_root   = tau_fine_root_p2
    param_pft2_from_file%tau_coarse_root   = tau_coarse_root_p2
    param_pft2_from_file%tau_branch   = tau_branch_p2
    param_pft2_from_file%tau_sap_wood   = tau_sap_wood_p2
    param_pft2_from_file%tau_fruit   = tau_fruit_p2
    param_pft2_from_file%tau_seed_litter   = tau_seed_litter_p2
    param_pft2_from_file%tau_seed_est   = tau_seed_est_p2
    param_pft2_from_file%tau_mycorrhiza   = tau_mycorrhiza_p2
    param_pft2_from_file%vmax_uptake_n   = vmax_uptake_n_p2
    param_pft2_from_file%vmax_uptake_p   = vmax_uptake_p_p2
    param_pft2_from_file%bnf_base   = bnf_base_p2
    param_pft2_from_file%lambda_est_light   = lambda_est_light_p2
    param_pft2_from_file%k_est_light   = k_est_light_p2
    param_pft2_from_file%seed_size   = seed_size_p2
    param_pft2_from_file%beta_soil_flush   = beta_soil_flush_p2
    param_pft2_from_file%beta_soil_senescence   = beta_soil_senescence_p2
    param_pft2_from_file%gdd_t_air_threshold   = gdd_t_air_threshold_p2
    param_pft2_from_file%gdd_req_max   = gdd_req_max_p2
    param_pft2_from_file%k_gdd_dormance   = k_gdd_dormance_p2
    param_pft2_from_file%t_air_senescence   = t_air_senescence_p2
    param_pft2_from_file%min_leaf_age   = min_leaf_age_p2
    param_pft2_from_file%frac_sapwood_branch   = frac_sapwood_branch_p2
    param_pft2_from_file%wood_density   = wood_density_p2
    param_pft2_from_file%k_latosa   = k_latosa_p2
    param_pft2_from_file%k_crtos   = k_crtos_p2
    param_pft2_from_file%k_rtos   = k_rtos_p2
    param_pft2_from_file%k2_fruit_alloc   = k2_fruit_alloc_p2
    param_pft2_from_file%allom_k1   = allom_k1_p2
    param_pft2_from_file%allom_k2   = allom_k2_p2
    param_pft2_from_file%phi_leaf_min   = phi_leaf_min_p2
    param_pft2_from_file%k_root   = k_root_p2
    param_pft2_from_file%k_sapwood   = k_sapwood_p2
    param_pft2_from_file%c0_allom   = c0_allom_p2
    param_pft2_from_file%fstore_target   = fstore_target_p2
    param_pft2_from_file%k_root_dist   = k_root_dist_p2
    ! LCTLIB_PFT3_CTL
    param_pft3_from_file%sla   = sla_p3
    param_pft3_from_file%sigma_vis   = sigma_vis_p3
    param_pft3_from_file%sigma_nir   = sigma_nir_p3
    param_pft3_from_file%omega_clumping   = omega_clumping_p3
    param_pft3_from_file%crown_shape_factor   = crown_shape_factor_p3
    param_pft3_from_file%cn_leaf   = cn_leaf_p3
    param_pft3_from_file%np_leaf   = np_leaf_p3
    param_pft3_from_file%k0_fn_struc   = k0_fn_struc_p3
    param_pft3_from_file%fn_oth_min   = fn_oth_min_p3
    param_pft3_from_file%t_jmax_opt   = t_jmax_opt_p3
    param_pft3_from_file%t_jmax_omega   = t_jmax_omega_p3
    param_pft3_from_file%g0   = g0_p3
    param_pft3_from_file%g1   = g1_p3
    param_pft3_from_file%gmin   = gmin_p3
    param_pft3_from_file%lambda_soil_ps   = lambda_soil_ps_p3
    param_pft3_from_file%weibull_k_soil_ps   = weibull_k_soil_ps_p3
    param_pft3_from_file%lambda_soil_gs   = lambda_soil_gs_p3
    param_pft3_from_file%weibull_k_soil_gs   = weibull_k_soil_gs_p3
    param_pft3_from_file%tau_leaf   = tau_leaf_p3
    param_pft3_from_file%tau_fine_root   = tau_fine_root_p3
    param_pft3_from_file%tau_coarse_root   = tau_coarse_root_p3
    param_pft3_from_file%tau_branch   = tau_branch_p3
    param_pft3_from_file%tau_sap_wood   = tau_sap_wood_p3
    param_pft3_from_file%tau_fruit   = tau_fruit_p3
    param_pft3_from_file%tau_seed_litter   = tau_seed_litter_p3
    param_pft3_from_file%tau_seed_est   = tau_seed_est_p3
    param_pft3_from_file%tau_mycorrhiza   = tau_mycorrhiza_p3
    param_pft3_from_file%vmax_uptake_n   = vmax_uptake_n_p3
    param_pft3_from_file%vmax_uptake_p   = vmax_uptake_p_p3
    param_pft3_from_file%bnf_base   = bnf_base_p3
    param_pft3_from_file%lambda_est_light   = lambda_est_light_p3
    param_pft3_from_file%k_est_light   = k_est_light_p3
    param_pft3_from_file%seed_size   = seed_size_p3
    param_pft3_from_file%beta_soil_flush   = beta_soil_flush_p3
    param_pft3_from_file%beta_soil_senescence   = beta_soil_senescence_p3
    param_pft3_from_file%gdd_t_air_threshold   = gdd_t_air_threshold_p3
    param_pft3_from_file%gdd_req_max   = gdd_req_max_p3
    param_pft3_from_file%k_gdd_dormance   = k_gdd_dormance_p3
    param_pft3_from_file%t_air_senescence   = t_air_senescence_p3
    param_pft3_from_file%min_leaf_age   = min_leaf_age_p3
    param_pft3_from_file%frac_sapwood_branch   = frac_sapwood_branch_p3
    param_pft3_from_file%wood_density   = wood_density_p3
    param_pft3_from_file%k_latosa   = k_latosa_p3
    param_pft3_from_file%k_crtos   = k_crtos_p3
    param_pft3_from_file%k_rtos   = k_rtos_p3
    param_pft3_from_file%k2_fruit_alloc   = k2_fruit_alloc_p3
    param_pft3_from_file%allom_k1   = allom_k1_p3
    param_pft3_from_file%allom_k2   = allom_k2_p3
    param_pft3_from_file%phi_leaf_min   = phi_leaf_min_p3
    param_pft3_from_file%k_root   = k_root_p3
    param_pft3_from_file%k_sapwood   = k_sapwood_p3
    param_pft3_from_file%c0_allom   = c0_allom_p3
    param_pft3_from_file%fstore_target   = fstore_target_p3
    param_pft3_from_file%k_root_dist   = k_root_dist_p3
    ! LCTLIB_PFT4_CTL
    param_pft4_from_file%sla   = sla_p4
    param_pft4_from_file%sigma_vis   = sigma_vis_p4
    param_pft4_from_file%sigma_nir   = sigma_nir_p4
    param_pft4_from_file%omega_clumping   = omega_clumping_p4
    param_pft4_from_file%crown_shape_factor   = crown_shape_factor_p4
    param_pft4_from_file%cn_leaf   = cn_leaf_p4
    param_pft4_from_file%np_leaf   = np_leaf_p4
    param_pft4_from_file%k0_fn_struc   = k0_fn_struc_p4
    param_pft4_from_file%fn_oth_min   = fn_oth_min_p4
    param_pft4_from_file%t_jmax_opt   = t_jmax_opt_p4
    param_pft4_from_file%t_jmax_omega   = t_jmax_omega_p4
    param_pft4_from_file%g0   = g0_p4
    param_pft4_from_file%g1   = g1_p4
    param_pft4_from_file%gmin   = gmin_p4
    param_pft4_from_file%lambda_soil_ps   = lambda_soil_ps_p4
    param_pft4_from_file%weibull_k_soil_ps   = weibull_k_soil_ps_p4
    param_pft4_from_file%lambda_soil_gs   = lambda_soil_gs_p4
    param_pft4_from_file%weibull_k_soil_gs   = weibull_k_soil_gs_p4
    param_pft4_from_file%tau_leaf   = tau_leaf_p4
    param_pft4_from_file%tau_fine_root   = tau_fine_root_p4
    param_pft4_from_file%tau_coarse_root   = tau_coarse_root_p4
    param_pft4_from_file%tau_branch   = tau_branch_p4
    param_pft4_from_file%tau_sap_wood   = tau_sap_wood_p4
    param_pft4_from_file%tau_fruit   = tau_fruit_p4
    param_pft4_from_file%tau_seed_litter   = tau_seed_litter_p4
    param_pft4_from_file%tau_seed_est   = tau_seed_est_p4
    param_pft4_from_file%tau_mycorrhiza   = tau_mycorrhiza_p4
    param_pft4_from_file%vmax_uptake_n   = vmax_uptake_n_p4
    param_pft4_from_file%vmax_uptake_p   = vmax_uptake_p_p4
    param_pft4_from_file%bnf_base   = bnf_base_p4
    param_pft4_from_file%lambda_est_light   = lambda_est_light_p4
    param_pft4_from_file%k_est_light   = k_est_light_p4
    param_pft4_from_file%seed_size   = seed_size_p4
    param_pft4_from_file%beta_soil_flush   = beta_soil_flush_p4
    param_pft4_from_file%beta_soil_senescence   = beta_soil_senescence_p4
    param_pft4_from_file%gdd_t_air_threshold   = gdd_t_air_threshold_p4
    param_pft4_from_file%gdd_req_max   = gdd_req_max_p4
    param_pft4_from_file%k_gdd_dormance   = k_gdd_dormance_p4
    param_pft4_from_file%t_air_senescence   = t_air_senescence_p4
    param_pft4_from_file%min_leaf_age   = min_leaf_age_p4
    param_pft4_from_file%frac_sapwood_branch   = frac_sapwood_branch_p4
    param_pft4_from_file%wood_density   = wood_density_p4
    param_pft4_from_file%k_latosa   = k_latosa_p4
    param_pft4_from_file%k_crtos   = k_crtos_p4
    param_pft4_from_file%k_rtos   = k_rtos_p4
    param_pft4_from_file%k2_fruit_alloc   = k2_fruit_alloc_p4
    param_pft4_from_file%allom_k1   = allom_k1_p4
    param_pft4_from_file%allom_k2   = allom_k2_p4
    param_pft4_from_file%phi_leaf_min   = phi_leaf_min_p4
    param_pft4_from_file%k_root   = k_root_p4
    param_pft4_from_file%k_sapwood   = k_sapwood_p4
    param_pft4_from_file%c0_allom   = c0_allom_p4
    param_pft4_from_file%fstore_target   = fstore_target_p4
    param_pft4_from_file%k_root_dist   = k_root_dist_p4
    ! LCTLIB_PFT5_CTL
    param_pft5_from_file%sla   = sla_p5
    param_pft5_from_file%sigma_vis   = sigma_vis_p5
    param_pft5_from_file%sigma_nir   = sigma_nir_p5
    param_pft5_from_file%omega_clumping   = omega_clumping_p5
    param_pft5_from_file%crown_shape_factor   = crown_shape_factor_p5
    param_pft5_from_file%cn_leaf   = cn_leaf_p5
    param_pft5_from_file%np_leaf   = np_leaf_p5
    param_pft5_from_file%k0_fn_struc   = k0_fn_struc_p5
    param_pft5_from_file%fn_oth_min   = fn_oth_min_p5
    param_pft5_from_file%t_jmax_opt   = t_jmax_opt_p5
    param_pft5_from_file%t_jmax_omega   = t_jmax_omega_p5
    param_pft5_from_file%g0   = g0_p5
    param_pft5_from_file%g1   = g1_p5
    param_pft5_from_file%gmin   = gmin_p5
    param_pft5_from_file%lambda_soil_ps   = lambda_soil_ps_p5
    param_pft5_from_file%weibull_k_soil_ps   = weibull_k_soil_ps_p5
    param_pft5_from_file%lambda_soil_gs   = lambda_soil_gs_p5
    param_pft5_from_file%weibull_k_soil_gs   = weibull_k_soil_gs_p5
    param_pft5_from_file%tau_leaf   = tau_leaf_p5
    param_pft5_from_file%tau_fine_root   = tau_fine_root_p5
    param_pft5_from_file%tau_coarse_root   = tau_coarse_root_p5
    param_pft5_from_file%tau_branch   = tau_branch_p5
    param_pft5_from_file%tau_sap_wood   = tau_sap_wood_p5
    param_pft5_from_file%tau_fruit   = tau_fruit_p5
    param_pft5_from_file%tau_seed_litter   = tau_seed_litter_p5
    param_pft5_from_file%tau_seed_est   = tau_seed_est_p5
    param_pft5_from_file%tau_mycorrhiza   = tau_mycorrhiza_p5
    param_pft5_from_file%vmax_uptake_n   = vmax_uptake_n_p5
    param_pft5_from_file%vmax_uptake_p   = vmax_uptake_p_p5
    param_pft5_from_file%bnf_base   = bnf_base_p5
    param_pft5_from_file%lambda_est_light   = lambda_est_light_p5
    param_pft5_from_file%k_est_light   = k_est_light_p5
    param_pft5_from_file%seed_size   = seed_size_p5
    param_pft5_from_file%beta_soil_flush   = beta_soil_flush_p5
    param_pft5_from_file%beta_soil_senescence   = beta_soil_senescence_p5
    param_pft5_from_file%gdd_t_air_threshold   = gdd_t_air_threshold_p5
    param_pft5_from_file%gdd_req_max   = gdd_req_max_p5
    param_pft5_from_file%k_gdd_dormance   = k_gdd_dormance_p5
    param_pft5_from_file%t_air_senescence   = t_air_senescence_p5
    param_pft5_from_file%min_leaf_age   = min_leaf_age_p5
    param_pft5_from_file%frac_sapwood_branch   = frac_sapwood_branch_p5
    param_pft5_from_file%wood_density   = wood_density_p5
    param_pft5_from_file%k_latosa   = k_latosa_p5
    param_pft5_from_file%k_crtos   = k_crtos_p5
    param_pft5_from_file%k_rtos   = k_rtos_p5
    param_pft5_from_file%k2_fruit_alloc   = k2_fruit_alloc_p5
    param_pft5_from_file%allom_k1   = allom_k1_p5
    param_pft5_from_file%allom_k2   = allom_k2_p5
    param_pft5_from_file%phi_leaf_min   = phi_leaf_min_p5
    param_pft5_from_file%k_root   = k_root_p5
    param_pft5_from_file%k_sapwood   = k_sapwood_p5
    param_pft5_from_file%c0_allom   = c0_allom_p5
    param_pft5_from_file%fstore_target   = fstore_target_p5
    param_pft5_from_file%k_root_dist   = k_root_dist_p5
    ! LCTLIB_PFT6_CTL
    param_pft6_from_file%sla   = sla_p6
    param_pft6_from_file%sigma_vis   = sigma_vis_p6
    param_pft6_from_file%sigma_nir   = sigma_nir_p6
    param_pft6_from_file%omega_clumping   = omega_clumping_p6
    param_pft6_from_file%crown_shape_factor   = crown_shape_factor_p6
    param_pft6_from_file%cn_leaf   = cn_leaf_p6
    param_pft6_from_file%np_leaf   = np_leaf_p6
    param_pft6_from_file%k0_fn_struc   = k0_fn_struc_p6
    param_pft6_from_file%fn_oth_min   = fn_oth_min_p6
    param_pft6_from_file%t_jmax_opt   = t_jmax_opt_p6
    param_pft6_from_file%t_jmax_omega   = t_jmax_omega_p6
    param_pft6_from_file%g0   = g0_p6
    param_pft6_from_file%g1   = g1_p6
    param_pft6_from_file%gmin   = gmin_p6
    param_pft6_from_file%lambda_soil_ps   = lambda_soil_ps_p6
    param_pft6_from_file%weibull_k_soil_ps   = weibull_k_soil_ps_p6
    param_pft6_from_file%lambda_soil_gs   = lambda_soil_gs_p6
    param_pft6_from_file%weibull_k_soil_gs   = weibull_k_soil_gs_p6
    param_pft6_from_file%tau_leaf   = tau_leaf_p6
    param_pft6_from_file%tau_fine_root   = tau_fine_root_p6
    param_pft6_from_file%tau_coarse_root   = tau_coarse_root_p6
    param_pft6_from_file%tau_branch   = tau_branch_p6
    param_pft6_from_file%tau_sap_wood   = tau_sap_wood_p6
    param_pft6_from_file%tau_fruit   = tau_fruit_p6
    param_pft6_from_file%tau_seed_litter   = tau_seed_litter_p6
    param_pft6_from_file%tau_seed_est   = tau_seed_est_p6
    param_pft6_from_file%tau_mycorrhiza   = tau_mycorrhiza_p6
    param_pft6_from_file%vmax_uptake_n   = vmax_uptake_n_p6
    param_pft6_from_file%vmax_uptake_p   = vmax_uptake_p_p6
    param_pft6_from_file%bnf_base   = bnf_base_p6
    param_pft6_from_file%lambda_est_light   = lambda_est_light_p6
    param_pft6_from_file%k_est_light   = k_est_light_p6
    param_pft6_from_file%seed_size   = seed_size_p6
    param_pft6_from_file%beta_soil_flush   = beta_soil_flush_p6
    param_pft6_from_file%beta_soil_senescence   = beta_soil_senescence_p6
    param_pft6_from_file%gdd_t_air_threshold   = gdd_t_air_threshold_p6
    param_pft6_from_file%gdd_req_max   = gdd_req_max_p6
    param_pft6_from_file%k_gdd_dormance   = k_gdd_dormance_p6
    param_pft6_from_file%t_air_senescence   = t_air_senescence_p6
    param_pft6_from_file%min_leaf_age   = min_leaf_age_p6
    param_pft6_from_file%frac_sapwood_branch   = frac_sapwood_branch_p6
    param_pft6_from_file%wood_density   = wood_density_p6
    param_pft6_from_file%k_latosa   = k_latosa_p6
    param_pft6_from_file%k_crtos   = k_crtos_p6
    param_pft6_from_file%k_rtos   = k_rtos_p6
    param_pft6_from_file%k2_fruit_alloc   = k2_fruit_alloc_p6
    param_pft6_from_file%allom_k1   = allom_k1_p6
    param_pft6_from_file%allom_k2   = allom_k2_p6
    param_pft6_from_file%phi_leaf_min   = phi_leaf_min_p6
    param_pft6_from_file%k_root   = k_root_p6
    param_pft6_from_file%k_sapwood   = k_sapwood_p6
    param_pft6_from_file%c0_allom   = c0_allom_p6
    param_pft6_from_file%fstore_target   = fstore_target_p6
    param_pft6_from_file%k_root_dist   = k_root_dist_p6
    ! LCTLIB_PFT7_CTL
    param_pft7_from_file%sla   = sla_p7
    param_pft7_from_file%sigma_vis   = sigma_vis_p7
    param_pft7_from_file%sigma_nir   = sigma_nir_p7
    param_pft7_from_file%omega_clumping   = omega_clumping_p7
    param_pft7_from_file%crown_shape_factor   = crown_shape_factor_p7
    param_pft7_from_file%cn_leaf   = cn_leaf_p7
    param_pft7_from_file%np_leaf   = np_leaf_p7
    param_pft7_from_file%k0_fn_struc   = k0_fn_struc_p7
    param_pft7_from_file%fn_oth_min   = fn_oth_min_p7
    param_pft7_from_file%t_jmax_opt   = t_jmax_opt_p7
    param_pft7_from_file%t_jmax_omega   = t_jmax_omega_p7
    param_pft7_from_file%g0   = g0_p7
    param_pft7_from_file%g1   = g1_p7
    param_pft7_from_file%gmin   = gmin_p7
    param_pft7_from_file%lambda_soil_ps   = lambda_soil_ps_p7
    param_pft7_from_file%weibull_k_soil_ps   = weibull_k_soil_ps_p7
    param_pft7_from_file%lambda_soil_gs   = lambda_soil_gs_p7
    param_pft7_from_file%weibull_k_soil_gs   = weibull_k_soil_gs_p7
    param_pft7_from_file%tau_leaf   = tau_leaf_p7
    param_pft7_from_file%tau_fine_root   = tau_fine_root_p7
    param_pft7_from_file%tau_coarse_root   = tau_coarse_root_p7
    param_pft7_from_file%tau_branch   = tau_branch_p7
    param_pft7_from_file%tau_sap_wood   = tau_sap_wood_p7
    param_pft7_from_file%tau_fruit   = tau_fruit_p7
    param_pft7_from_file%tau_seed_litter   = tau_seed_litter_p7
    param_pft7_from_file%tau_seed_est   = tau_seed_est_p7
    param_pft7_from_file%tau_mycorrhiza   = tau_mycorrhiza_p7
    param_pft7_from_file%vmax_uptake_n   = vmax_uptake_n_p7
    param_pft7_from_file%vmax_uptake_p   = vmax_uptake_p_p7
    param_pft7_from_file%bnf_base   = bnf_base_p7
    param_pft7_from_file%lambda_est_light   = lambda_est_light_p7
    param_pft7_from_file%k_est_light   = k_est_light_p7
    param_pft7_from_file%seed_size   = seed_size_p7
    param_pft7_from_file%beta_soil_flush   = beta_soil_flush_p7
    param_pft7_from_file%beta_soil_senescence   = beta_soil_senescence_p7
    param_pft7_from_file%gdd_t_air_threshold   = gdd_t_air_threshold_p7
    param_pft7_from_file%gdd_req_max   = gdd_req_max_p7
    param_pft7_from_file%k_gdd_dormance   = k_gdd_dormance_p7
    param_pft7_from_file%t_air_senescence   = t_air_senescence_p7
    param_pft7_from_file%min_leaf_age   = min_leaf_age_p7
    param_pft7_from_file%frac_sapwood_branch   = frac_sapwood_branch_p7
    param_pft7_from_file%wood_density   = wood_density_p7
    param_pft7_from_file%k_latosa   = k_latosa_p7
    param_pft7_from_file%k_crtos   = k_crtos_p7
    param_pft7_from_file%k_rtos   = k_rtos_p7
    param_pft7_from_file%k2_fruit_alloc   = k2_fruit_alloc_p7
    param_pft7_from_file%allom_k1   = allom_k1_p7
    param_pft7_from_file%allom_k2   = allom_k2_p7
    param_pft7_from_file%phi_leaf_min   = phi_leaf_min_p7
    param_pft7_from_file%k_root   = k_root_p7
    param_pft7_from_file%k_sapwood   = k_sapwood_p7
    param_pft7_from_file%c0_allom   = c0_allom_p7
    param_pft7_from_file%fstore_target   = fstore_target_p7
    param_pft7_from_file%k_root_dist   = k_root_dist_p7
    ! LCTLIB_PFT8_CTL
    param_pft8_from_file%sla   = sla_p8
    param_pft8_from_file%sigma_vis   = sigma_vis_p8
    param_pft8_from_file%sigma_nir   = sigma_nir_p8
    param_pft8_from_file%omega_clumping   = omega_clumping_p8
    param_pft8_from_file%crown_shape_factor   = crown_shape_factor_p8
    param_pft8_from_file%cn_leaf   = cn_leaf_p8
    param_pft8_from_file%np_leaf   = np_leaf_p8
    param_pft8_from_file%k0_fn_struc   = k0_fn_struc_p8
    param_pft8_from_file%fn_oth_min   = fn_oth_min_p8
    param_pft8_from_file%t_jmax_opt   = t_jmax_opt_p8
    param_pft8_from_file%t_jmax_omega   = t_jmax_omega_p8
    param_pft8_from_file%g0   = g0_p8
    param_pft8_from_file%g1   = g1_p8
    param_pft8_from_file%gmin   = gmin_p8
    param_pft8_from_file%lambda_soil_ps   = lambda_soil_ps_p8
    param_pft8_from_file%weibull_k_soil_ps   = weibull_k_soil_ps_p8
    param_pft8_from_file%lambda_soil_gs   = lambda_soil_gs_p8
    param_pft8_from_file%weibull_k_soil_gs   = weibull_k_soil_gs_p8
    param_pft8_from_file%tau_leaf   = tau_leaf_p8
    param_pft8_from_file%tau_fine_root   = tau_fine_root_p8
    param_pft8_from_file%tau_coarse_root   = tau_coarse_root_p8
    param_pft8_from_file%tau_branch   = tau_branch_p8
    param_pft8_from_file%tau_sap_wood   = tau_sap_wood_p8
    param_pft8_from_file%tau_fruit   = tau_fruit_p8
    param_pft8_from_file%tau_seed_litter   = tau_seed_litter_p8
    param_pft8_from_file%tau_seed_est   = tau_seed_est_p8
    param_pft8_from_file%tau_mycorrhiza   = tau_mycorrhiza_p8
    param_pft8_from_file%vmax_uptake_n   = vmax_uptake_n_p8
    param_pft8_from_file%vmax_uptake_p   = vmax_uptake_p_p8
    param_pft8_from_file%bnf_base   = bnf_base_p8
    param_pft8_from_file%lambda_est_light   = lambda_est_light_p8
    param_pft8_from_file%k_est_light   = k_est_light_p8
    param_pft8_from_file%seed_size   = seed_size_p8
    param_pft8_from_file%beta_soil_flush   = beta_soil_flush_p8
    param_pft8_from_file%beta_soil_senescence   = beta_soil_senescence_p8
    param_pft8_from_file%gdd_t_air_threshold   = gdd_t_air_threshold_p8
    param_pft8_from_file%gdd_req_max   = gdd_req_max_p8
    param_pft8_from_file%k_gdd_dormance   = k_gdd_dormance_p8
    param_pft8_from_file%t_air_senescence   = t_air_senescence_p8
    param_pft8_from_file%min_leaf_age   = min_leaf_age_p8
    param_pft8_from_file%frac_sapwood_branch   = frac_sapwood_branch_p8
    param_pft8_from_file%wood_density   = wood_density_p8
    param_pft8_from_file%k_latosa   = k_latosa_p8
    param_pft8_from_file%k_crtos   = k_crtos_p8
    param_pft8_from_file%k_rtos   = k_rtos_p8
    param_pft8_from_file%k2_fruit_alloc   = k2_fruit_alloc_p8
    param_pft8_from_file%allom_k1   = allom_k1_p8
    param_pft8_from_file%allom_k2   = allom_k2_p8
    param_pft8_from_file%phi_leaf_min   = phi_leaf_min_p8
    param_pft8_from_file%k_root   = k_root_p8
    param_pft8_from_file%k_sapwood   = k_sapwood_p8
    param_pft8_from_file%c0_allom   = c0_allom_p8
    param_pft8_from_file%fstore_target   = fstore_target_p8
    param_pft8_from_file%k_root_dist   = k_root_dist_p8

  END SUBROUTINE read_parameter_lctlib_list


  !-----------------------------------------------------------------------------------------------------
  !> replace hard coded lctlib parameter values with values from file
  !! 
  !! 
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE set_parameter_lctlib_values_with_list_values(model)

    USE mo_jsb_model_class,                 ONLY: t_jsb_model

    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    TYPE(t_jsb_model), INTENT(inout) :: model
    ! ---------------------------
    ! 0.2 Local
    INTEGER                     :: pft_id 
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':set_parameter_lctlib_values_with_list_values'


    !> 1.0 read parameter lctlib values from file
    !! 
    CALL read_parameter_lctlib_list

    !> 2.0 replace parameter lctlib values for each PFT individually
    !! 
    !! lctlib(4)%sla is needed for calculation of k_rtos
    CALL replace_pft_specific_lctlib_parameter_values(model%lctlib(1), param_pft1_from_file, model%lctlib(4)%sla)
    CALL replace_pft_specific_lctlib_parameter_values(model%lctlib(2), param_pft2_from_file, model%lctlib(4)%sla)
    CALL replace_pft_specific_lctlib_parameter_values(model%lctlib(3), param_pft3_from_file, model%lctlib(4)%sla)
    CALL replace_pft_specific_lctlib_parameter_values(model%lctlib(4), param_pft4_from_file, model%lctlib(4)%sla)
    CALL replace_pft_specific_lctlib_parameter_values(model%lctlib(5), param_pft5_from_file, model%lctlib(4)%sla)
    CALL replace_pft_specific_lctlib_parameter_values(model%lctlib(6), param_pft6_from_file, model%lctlib(4)%sla)
    CALL replace_pft_specific_lctlib_parameter_values(model%lctlib(7), param_pft7_from_file, model%lctlib(4)%sla)
    CALL replace_pft_specific_lctlib_parameter_values(model%lctlib(8), param_pft8_from_file, model%lctlib(4)%sla)

  END SUBROUTINE set_parameter_lctlib_values_with_list_values


  !-----------------------------------------------------------------------------------------------------
  !> replace hard coded lctlib parameter values for one PFT with values from file
  !! 
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE replace_pft_specific_lctlib_parameter_values(lctlib, param_pft_from_file, pft4_sla)

    USE mo_jsb_lctlib_class,        ONLY: t_lctlib_element
    USE mo_base_constants,          ONLY: one_day, one_year
    USE mo_shared_constants,        ONLY: molar_mass_C, molar_mass_N, molar_mass_P, Dwv, Dco2, Tzero 
    USE mo_veg_constants,           ONLY: carbon_per_dryweight_leaf, sm2lm_grass, igrass, itree

    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    TYPE(t_lctlib_element),               INTENT(inout) :: lctlib
    TYPE(t_parameters_lctlib_from_file),  INTENT(in)    :: param_pft_from_file
    REAL(wp),                             INTENT(in)    :: pft4_sla
    ! ---------------------------
    ! 0.2 Local
    REAL(wp)                    :: pft4_sla_hlp
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':replace_pft_specific_lctlib_parameter_values'


    !> 0.9 init local variable
    !! 
    pft4_sla_hlp = pft4_sla


    !> 1.0 replace lctlib parameter-value for the particular PFT if a new value is provided
    !! 
! example:
!   IF (ABS(param_pft_from_file%frac_sapwood_branch       - missval) > eps8) &
!      lctlib%frac_sapwood_branch = param_pft_from_file%frac_sapwood_branch

    IF (ABS(param_pft_from_file%sla                       - missval) > eps8) THEN
       lctlib%sla = 1._wp / carbon_per_dryweight_leaf / 1000.0_wp * molar_mass_C * param_pft_from_file%sla
      ! change pft4_sla_hlp only of sla of PFT 4 was changed; see calculation of k_rtos below
      IF (lctlib%lctnumber == 4) pft4_sla_hlp = lctlib%sla
    ENDIF
    IF (ABS(param_pft_from_file%sigma_vis                 - missval) > eps8) &
       lctlib%sigma_vis = param_pft_from_file%sigma_vis
    IF (ABS(param_pft_from_file%sigma_nir                 - missval) > eps8) &
       lctlib%sigma_nir = param_pft_from_file%sigma_nir
    IF (ABS(param_pft_from_file%omega_clumping            - missval) > eps8) &
       lctlib%omega_clumping = param_pft_from_file%omega_clumping
    IF (ABS(param_pft_from_file%crown_shape_factor        - missval) > eps8) &
       lctlib%crown_shape_factor = param_pft_from_file%crown_shape_factor
    IF (ABS(param_pft_from_file%cn_leaf                   - missval) > eps8) &
       lctlib%cn_leaf = carbon_per_dryweight_leaf * 1000._wp * molar_mass_N / molar_mass_C / param_pft_from_file%cn_leaf
    IF (ABS(param_pft_from_file%np_leaf                   - missval) > eps8) &
       lctlib%np_leaf = carbon_per_dryweight_leaf * 1000._wp * molar_mass_P / molar_mass_C / lctlib%cn_leaf / & 
                        param_pft_from_file%np_leaf
    IF (ABS(param_pft_from_file%k0_fn_struc               - missval) > eps8) &
       lctlib%k0_fn_struc = param_pft_from_file%k0_fn_struc
    ! fn_oth_min might not be changed directly but via k0_fn_struc
    ! IF (ABS(param_pft_from_file%fn_oth_min                - missval) > eps8) &
    !    lctlib%fn_oth_min = param_pft_from_file%fn_oth_min
    IF (ABS(param_pft_from_file%k0_fn_struc               - missval) > eps8) &
       lctlib%fn_oth_min = param_pft_from_file%k0_fn_struc - 5.0_wp / molar_mass_n * 0.5_wp * 1.e3_wp * 71.4_wp * &
                           molar_mass_N / 1.e6_wp
    IF (ABS(param_pft_from_file%t_jmax_opt                - missval) > eps8) &
       lctlib%t_jmax_opt = param_pft_from_file%t_jmax_opt
    IF (ABS(param_pft_from_file%t_jmax_omega              - missval) > eps8) &
       lctlib%t_jmax_omega = param_pft_from_file%t_jmax_omega
    IF (ABS(param_pft_from_file%g0                        - missval) > eps8) &
       lctlib%g0 = param_pft_from_file%g0
    IF (ABS(param_pft_from_file%g1                        - missval) > eps8) &
       lctlib%g1 = param_pft_from_file%g1
    IF (ABS(param_pft_from_file%gmin                      - missval) > eps8) &
       lctlib%gmin = param_pft_from_file%gmin
    IF (ABS(param_pft_from_file%lambda_soil_ps            - missval) > eps8) &
       lctlib%lambda_soil_ps = param_pft_from_file%lambda_soil_ps
    IF (ABS(param_pft_from_file%weibull_k_soil_ps         - missval) > eps8) &
       lctlib%weibull_k_soil_ps = param_pft_from_file%weibull_k_soil_ps
    IF (ABS(param_pft_from_file%lambda_soil_gs            - missval) > eps8) &
       lctlib%lambda_soil_gs = param_pft_from_file%lambda_soil_gs
    IF (ABS(param_pft_from_file%weibull_k_soil_gs         - missval) > eps8) &
       lctlib%weibull_k_soil_gs = param_pft_from_file%weibull_k_soil_gs
    IF (ABS(param_pft_from_file%tau_leaf                  - missval) > eps8) &
       lctlib%tau_leaf = param_pft_from_file%tau_leaf / 12.0_wp
    IF (ABS(param_pft_from_file%tau_fine_root             - missval) > eps8) &
       lctlib%tau_fine_root = param_pft_from_file%tau_fine_root
    IF (ABS(param_pft_from_file%tau_coarse_root           - missval) > eps8) &
       lctlib%tau_coarse_root = param_pft_from_file%tau_coarse_root
    IF (ABS(param_pft_from_file%tau_branch                - missval) > eps8) &
       lctlib%tau_branch = param_pft_from_file%tau_branch
    IF (ABS(param_pft_from_file%tau_sap_wood              - missval) > eps8) &
       lctlib%tau_sap_wood = param_pft_from_file%tau_sap_wood
    IF (ABS(param_pft_from_file%tau_fruit                 - missval) > eps8) &
       lctlib%tau_fruit = param_pft_from_file%tau_fruit
    IF (ABS(param_pft_from_file%tau_seed_litter           - missval) > eps8) &
       lctlib%tau_seed_litter = param_pft_from_file%tau_seed_litter
    IF (ABS(param_pft_from_file%tau_seed_est              - missval) > eps8) &
       lctlib%tau_seed_est = param_pft_from_file%tau_seed_est
    IF (ABS(param_pft_from_file%tau_mycorrhiza            - missval) > eps8) &
       lctlib%tau_mycorrhiza = param_pft_from_file%tau_mycorrhiza
    IF (ABS(param_pft_from_file%vmax_uptake_n             - missval) > eps8) &
       lctlib%vmax_uptake_n = param_pft_from_file%vmax_uptake_n
    IF (ABS(param_pft_from_file%vmax_uptake_p             - missval) > eps8) &
       lctlib%vmax_uptake_p = param_pft_from_file%vmax_uptake_p
    IF (ABS(param_pft_from_file%bnf_base                  - missval) > eps8) &
       lctlib%bnf_base = param_pft_from_file%bnf_base / molar_mass_N / one_day / one_year * 1.e6_wp
    IF (ABS(param_pft_from_file%lambda_est_light          - missval) > eps8) &
       lctlib%lambda_est_light = param_pft_from_file%lambda_est_light
    IF (ABS(param_pft_from_file%k_est_light               - missval) > eps8) &
       lctlib%k_est_light = param_pft_from_file%k_est_light
    IF (ABS(param_pft_from_file%seed_size                 - missval) > eps8) &
       lctlib%seed_size = param_pft_from_file%seed_size
    IF (ABS(param_pft_from_file%beta_soil_flush           - missval) > eps8) &
       lctlib%beta_soil_flush = param_pft_from_file%beta_soil_flush
    IF (ABS(param_pft_from_file%beta_soil_senescence      - missval) > eps8) &
       lctlib%beta_soil_senescence = param_pft_from_file%beta_soil_senescence
    IF (ABS(param_pft_from_file%gdd_t_air_threshold       - missval) > eps8) &
       lctlib%gdd_t_air_threshold = param_pft_from_file%gdd_t_air_threshold + Tzero
    IF (ABS(param_pft_from_file%gdd_req_max               - missval) > eps8) &
       lctlib%gdd_req_max = param_pft_from_file%gdd_req_max
    IF (ABS(param_pft_from_file%k_gdd_dormance            - missval) > eps8) &
       lctlib%k_gdd_dormance = param_pft_from_file%k_gdd_dormance
    IF (ABS(param_pft_from_file%t_air_senescence          - missval) > eps8) &
       lctlib%t_air_senescence = param_pft_from_file%t_air_senescence + Tzero
    IF (ABS(param_pft_from_file%min_leaf_age              - missval) > eps8) &
       lctlib%min_leaf_age = param_pft_from_file%min_leaf_age
    IF (ABS(param_pft_from_file%frac_sapwood_branch       - missval) > eps8) &
       lctlib%frac_sapwood_branch = param_pft_from_file%frac_sapwood_branch
    IF (ABS(param_pft_from_file%wood_density              - missval) > eps8) &
       lctlib%wood_density = param_pft_from_file%wood_density * 1.e6_wp / molar_mass_C
    IF (ABS(param_pft_from_file%k_latosa                  - missval) > eps8) &
       lctlib%k_latosa = param_pft_from_file%k_latosa
    IF (ABS(param_pft_from_file%k_crtos                   - missval) > eps8) &
       lctlib%k_crtos = param_pft_from_file%k_crtos
    ! k_rtos might not be changed directly but via other lctlib parameters; code moved to below considering dependencies to other lctlib parameters
    ! IF (ABS(param_pft_from_file%k_rtos                    - missval) > eps8) &
    !    lctlib%k_rtos = param_pft_from_file%k_rtos
    IF (ABS(param_pft_from_file%k2_fruit_alloc            - missval) > eps8) &
       lctlib%k2_fruit_alloc = param_pft_from_file%k2_fruit_alloc
    IF (ABS(param_pft_from_file%allom_k1                  - missval) > eps8) &
       lctlib%allom_k1 = param_pft_from_file%allom_k1
    IF (ABS(param_pft_from_file%allom_k2                  - missval) > eps8) &
       lctlib%allom_k2 = param_pft_from_file%allom_k2
    IF (ABS(param_pft_from_file%phi_leaf_min              - missval) > eps8) &
       lctlib%phi_leaf_min = param_pft_from_file%phi_leaf_min
    IF (ABS(param_pft_from_file%k_root                    - missval) > eps8) &
       lctlib%k_root = param_pft_from_file%k_root * molar_mass_C * 1.e-10_wp
    IF (ABS(param_pft_from_file%k_sapwood                 - missval) > eps8) &
       lctlib%k_sapwood = 1.e-3_wp * param_pft_from_file%k_sapwood
    ! k_rtos
    IF (ABS(param_pft_from_file%wood_density              - missval) > eps8 .OR. &
        ABS(param_pft_from_file%k_latosa                  - missval) > eps8        ) THEN
      IF(lctlib%growthform /= igrass) THEN
         lctlib%k_rtos = 1.75_wp * (pft4_sla_hlp * lctlib%wood_density) / lctlib%k_latosa
      ELSE
         lctlib%k_rtos = 0.5_wp / sm2lm_grass
      ENDIF
    ENDIF
    ! c0_allom might not be changed directly but via other lctlib parameters:
    ! IF (ABS(param_pft_from_file%c0_allom                  - missval) > eps8) &
    !    lctlib%c0_allom = param_pft_from_file%c0_allom
    IF (ABS(param_pft_from_file%k_root                    - missval) > eps8 .OR. &
        ABS(param_pft_from_file%k_sapwood                 - missval) > eps8 .OR. &
        ABS(param_pft_from_file%tau_sap_wood              - missval) > eps8 .OR. &
        ABS(param_pft_from_file%tau_fine_root             - missval) > eps8 .OR. &
        ABS(param_pft_from_file%wood_density              - missval) > eps8        ) THEN
      IF(lctlib%growthform == itree) THEN
        lctlib%c0_allom = SQRT(lctlib%k_root / lctlib%k_sapwood * lctlib%tau_sap_wood / lctlib%tau_fine_root * &
                               lctlib%wood_density) * & 
                          lctlib%tau_fine_root / lctlib%tau_sap_wood
      ELSE
        lctlib%c0_allom = 0.0_wp
      ENDIF
    ENDIF
    IF (ABS(param_pft_from_file%fstore_target             - missval) > eps8) &
       lctlib%fstore_target = param_pft_from_file%fstore_target
    IF (ABS(param_pft_from_file%k_root_dist               - missval) > eps8) &
       lctlib%k_root_dist = param_pft_from_file%k_root_dist


  END SUBROUTINE replace_pft_specific_lctlib_parameter_values



END MODULE mo_shared_set_parameters
