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
    REAL(wp)            :: c0_allom

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
    REAL(wp)            :: c0_allom_p1
    ! PFT 2
    REAL(wp)            :: sla_p2
    REAL(wp)            :: sigma_vis_p2
    REAL(wp)            :: sigma_nir_p2
    REAL(wp)            :: c0_allom_p2
    ! PFT 3
    REAL(wp)            :: sla_p3
    REAL(wp)            :: sigma_vis_p3
    REAL(wp)            :: sigma_nir_p3
    REAL(wp)            :: c0_allom_p3
    ! PFT 4
    REAL(wp)            :: sla_p4
    REAL(wp)            :: sigma_vis_p4
    REAL(wp)            :: sigma_nir_p4
    REAL(wp)            :: c0_allom_p4        
    ! PFT 5
    REAL(wp)            :: sla_p5
    REAL(wp)            :: sigma_vis_p5
    REAL(wp)            :: sigma_nir_p5
    REAL(wp)            :: c0_allom_p5
    ! PFT 6
    REAL(wp)            :: sla_p6
    REAL(wp)            :: sigma_vis_p6
    REAL(wp)            :: sigma_nir_p6
    REAL(wp)            :: c0_allom_p6
    ! PFT 7
    REAL(wp)            :: sla_p7
    REAL(wp)            :: sigma_vis_p7
    REAL(wp)            :: sigma_nir_p7
    REAL(wp)            :: c0_allom_p7
    ! PFT 8
    REAL(wp)            :: sla_p8
    REAL(wp)            :: sigma_vis_p8
    REAL(wp)            :: sigma_nir_p8
    REAL(wp)            :: c0_allom_p8


    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':read_parameter_lctlib_list'

    ! Namelist definitions for each PFT
    ! 
    NAMELIST /LCTLIB_PFT1_CTL/      &  
          sla_p1                  , &
          sigma_vis_p1            , &
          sigma_nir_p1            , &
          c0_allom_p1             
    NAMELIST /LCTLIB_PFT2_CTL/      &  
          sla_p2                  , &
          sigma_vis_p2            , &
          sigma_nir_p2            , &
          c0_allom_p2             
    NAMELIST /LCTLIB_PFT3_CTL/      &  
          sla_p3                  , &
          sigma_vis_p3            , &
          sigma_nir_p3            , &
          c0_allom_p3             
    NAMELIST /LCTLIB_PFT4_CTL/      &  
          sla_p4                  , &
          sigma_vis_p4            , &
          sigma_nir_p4            , &
          c0_allom_p4             
    NAMELIST /LCTLIB_PFT5_CTL/      &  
          sla_p5                  , &
          sigma_vis_p5            , &
          sigma_nir_p5            , &
          c0_allom_p5             
    NAMELIST /LCTLIB_PFT6_CTL/      &  
          sla_p6                  , &
          sigma_vis_p6            , &
          sigma_nir_p6            , &
          c0_allom_p6             
    NAMELIST /LCTLIB_PFT7_CTL/      &  
          sla_p7                  , &
          sigma_vis_p7            , &
          sigma_nir_p7            , &
          c0_allom_p7             
    NAMELIST /LCTLIB_PFT8_CTL/      &  
          sla_p8                  , &
          sigma_vis_p8            , &
          sigma_nir_p8            , &
          c0_allom_p8             


    ! variables for reading model-options from namelist
    ! 
    INTEGER  :: nml_handler, nml_unit, istat

    CALL message(TRIM(routine), 'Starting to read parameter list')

    ! Set defaults
    ! 
    ! PFT 1
    sla_p1                    = missval
    sigma_vis_p1              = missval
    sigma_nir_p1              = missval
    c0_allom_p1               = missval
    ! PFT 2
    sla_p2                    = missval
    sigma_vis_p2              = missval
    sigma_nir_p2              = missval
    c0_allom_p2               = missval
    ! PFT 3
    sla_p3                    = missval
    sigma_vis_p3              = missval
    sigma_nir_p3              = missval
    c0_allom_p3               = missval
    ! PFT 4
    sla_p4                    = missval
    sigma_vis_p4              = missval
    sigma_nir_p4              = missval
    c0_allom_p4               = missval
    ! PFT 5
    sla_p5                    = missval
    sigma_vis_p5              = missval
    sigma_nir_p5              = missval
    c0_allom_p5               = missval
    ! PFT 6
    sla_p6                    = missval
    sigma_vis_p6              = missval
    sigma_nir_p6              = missval
    c0_allom_p6               = missval
    ! PFT 7
    sla_p7                    = missval
    sigma_vis_p7              = missval
    sigma_nir_p7              = missval
    c0_allom_p7               = missval
    ! PFT 8
    sla_p8                    = missval
    sigma_vis_p8              = missval
    sigma_nir_p8              = missval
    c0_allom_p8               = missval


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
    param_pft1_from_file%sla                    = sla_p1      
    param_pft1_from_file%sigma_vis              = sigma_vis_p1
    param_pft1_from_file%sigma_nir              = sigma_nir_p1
    param_pft1_from_file%c0_allom               = c0_allom_p1 
    ! LCTLIB_PFT2_CTL
    param_pft2_from_file%sla                    = sla_p2      
    param_pft2_from_file%sigma_vis              = sigma_vis_p2
    param_pft2_from_file%sigma_nir              = sigma_nir_p2
    param_pft2_from_file%c0_allom               = c0_allom_p2 
    ! LCTLIB_PFT3_CTL
    param_pft3_from_file%sla                    = sla_p3      
    param_pft3_from_file%sigma_vis              = sigma_vis_p3
    param_pft3_from_file%sigma_nir              = sigma_nir_p3
    param_pft3_from_file%c0_allom               = c0_allom_p3 
    ! LCTLIB_PFT4_CTL
    param_pft4_from_file%sla                    = sla_p4      
    param_pft4_from_file%sigma_vis              = sigma_vis_p4
    param_pft4_from_file%sigma_nir              = sigma_nir_p4
    param_pft4_from_file%c0_allom               = c0_allom_p4 
    ! LCTLIB_PFT5_CTL
    param_pft5_from_file%sla                    = sla_p5      
    param_pft5_from_file%sigma_vis              = sigma_vis_p5
    param_pft5_from_file%sigma_nir              = sigma_nir_p5
    param_pft5_from_file%c0_allom               = c0_allom_p5 
    ! LCTLIB_PFT6_CTL
    param_pft6_from_file%sla                    = sla_p6      
    param_pft6_from_file%sigma_vis              = sigma_vis_p6
    param_pft6_from_file%sigma_nir              = sigma_nir_p6
    param_pft6_from_file%c0_allom               = c0_allom_p6 
    ! LCTLIB_PFT7_CTL
    param_pft7_from_file%sla                    = sla_p7      
    param_pft7_from_file%sigma_vis              = sigma_vis_p7
    param_pft7_from_file%sigma_nir              = sigma_nir_p7
    param_pft7_from_file%c0_allom               = c0_allom_p7 
    ! LCTLIB_PFT8_CTL
    param_pft8_from_file%sla                    = sla_p8      
    param_pft8_from_file%sigma_vis              = sigma_vis_p8
    param_pft8_from_file%sigma_nir              = sigma_nir_p8
    param_pft8_from_file%c0_allom               = c0_allom_p8 

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
    CALL replace_pft_specific_lctlib_parameter_values(model%lctlib(1), param_pft1_from_file)
    CALL replace_pft_specific_lctlib_parameter_values(model%lctlib(2), param_pft2_from_file)
    CALL replace_pft_specific_lctlib_parameter_values(model%lctlib(3), param_pft3_from_file)
    CALL replace_pft_specific_lctlib_parameter_values(model%lctlib(4), param_pft4_from_file)
    CALL replace_pft_specific_lctlib_parameter_values(model%lctlib(5), param_pft5_from_file)
    CALL replace_pft_specific_lctlib_parameter_values(model%lctlib(6), param_pft6_from_file)
    CALL replace_pft_specific_lctlib_parameter_values(model%lctlib(7), param_pft7_from_file)
    CALL replace_pft_specific_lctlib_parameter_values(model%lctlib(8), param_pft8_from_file)

  END SUBROUTINE set_parameter_lctlib_values_with_list_values


  !-----------------------------------------------------------------------------------------------------
  !> replace hard coded lctlib parameter values for one PFT with values from file
  !! 
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE replace_pft_specific_lctlib_parameter_values(lctlib, param_pft_from_file)

    USE mo_jsb_lctlib_class,        ONLY: t_lctlib_element
    USE mo_base_constants,          ONLY: one_day, one_year
    USE mo_shared_constants,        ONLY: molar_mass_C, molar_mass_N, molar_mass_P, Dwv, Dco2, Tzero 
    USE mo_veg_constants,           ONLY: carbon_per_dryweight_leaf, sm2lm_grass, igrass, itree

    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    TYPE(t_lctlib_element),               INTENT(inout) :: lctlib
    TYPE(t_parameters_lctlib_from_file),  INTENT(in)    :: param_pft_from_file
    ! ---------------------------
    ! 0.2 Local
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':replace_pft_specific_lctlib_parameter_values'


    !> 1.0 replace lctlib parameter-value for the particular PFT if a new value is provided
    !! 
!   IF (ABS(param_pft_from_file%frac_sapwood_branch       - missval) > eps8) &
!      lctlib%frac_sapwood_branch = param_pft_from_file%frac_sapwood_branch

    IF (ABS(param_pft_from_file%sla                       - missval) > eps8) &
       lctlib%sla = 1._wp / carbon_per_dryweight_leaf / 1000.0_wp * molar_mass_C * param_pft_from_file%sla



    IF (ABS(param_pft_from_file%sigma_vis                 - missval) > eps8) &
       lctlib%sigma_vis = param_pft_from_file%sigma_vis
    IF (ABS(param_pft_from_file%sigma_nir                 - missval) > eps8) &
       lctlib%sigma_nir = param_pft_from_file%sigma_nir



    ! IF (ABS(param_pft_from_file%k_root                    - missval) > eps8 .OR. &
    !     ABS(param_pft_from_file%k_sapwood                 - missval) > eps8 .OR. &
    !     ABS(param_pft_from_file%tau_sap_wood              - missval) > eps8 .OR. &
    !     ABS(param_pft_from_file%tau_fine_root             - missval) > eps8 .OR. &
    !     ABS(param_pft_from_file%wood_density              - missval) > eps8        ) THEN
    !   IF(lctlib%growthform == itree) THEN
    !     lctlib%c0_allom = SQRT(lctlib%k_root / lctlib%k_sapwood * lctlib%tau_sap_wood / lctlib%tau_fine_root * &
    !                            lctlib%wood_density) * & 
    !                       lctlib%tau_fine_root / lctlib%tau_sap_wood
    !   ELSE
    !     lctlib%c0_allom = 0.0_wp
    !   ENDIF
    ! ENDIF



  END SUBROUTINE replace_pft_specific_lctlib_parameter_values



END MODULE mo_shared_set_parameters
