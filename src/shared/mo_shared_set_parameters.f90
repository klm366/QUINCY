!>
!! @brief Contains the routines for reading parameters from param.list and replace parameter-values as defined in mo_PROC_constants
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
  PUBLIC :: parameters_from_file

  !-----------------------------------------------------------------------------------------------------
  !> parameter values from file param.list - replacing the values from mo_PROC_constants
  !! 
  !! 
  !-----------------------------------------------------------------------------------------------------
  TYPE t_parameters_from_file
    !LOGICAL            :: active   ! de/activate this TYPE (this is non-functional for now 20161004)
    CHARACTER(LEN=10)  :: filename = 'param.list'                 ! filename 
    ! SHARED_CTL
    REAL(wp)           :: mavg_period_tlabile
    REAL(wp)           :: mavg_period_tphen
    REAL(wp)           :: mavg_period_tuptake
    REAL(wp)           :: mavg_period_tfrac
    REAL(wp)           :: mavg_period_tcnl
    REAL(wp)           :: mavg_period_talloc
    REAL(wp)           :: mavg_period_tvegdyn
    REAL(wp)           :: mavg_period_tacclim
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
    REAL(wp)           :: microbial_cue
    REAL(wp)           :: residue_cue
  CONTAINS
    PROCEDURE :: read_param => read_parameter_list
    PROCEDURE :: set_param => set_parameter_values_with_list_values
  END Type t_parameters_from_file

  TYPE(t_parameters_from_file), SAVE :: parameters_from_file


  CHARACTER(len=*), PARAMETER :: modname = 'mo_shared_set_parameters'

CONTAINS


  !-----------------------------------------------------------------------------------------------------
  !> read parameter values from file
  !! 
  !! 
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE Read_parameter_list(this)

    USE mo_namelist,           ONLY: open_nml, POSITIONED, position_nml, close_nml  ! read model options from file

    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    CLASS(t_parameters_from_file), INTENT(inout) :: this
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
    REAL(wp)           :: microbial_cue
    REAL(wp)           :: residue_cue

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':Read_parameter_list'

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
        mavg_period_tacclim

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
        microbial_cue             , &
        residue_cue

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
    microbial_cue               = missval
    residue_cue                 = missval

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
    this%microbial_cue               = microbial_cue
    this%residue_cue                 = residue_cue

  END SUBROUTINE Read_parameter_list


  !-----------------------------------------------------------------------------------------------------
  !> replace hard coded parameter values with values from file
  !! 
  !! 
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE set_parameter_values_with_list_values(this)

    USE mo_base_constants,                  ONLY: one_day, one_year, dtime
    !USE mo_assimi_constants
    !USE mo_atmland_constants
    !USE mo_rad_constants
    USE mo_shared_constants
    !USE mo_soil_constants
    USE mo_sb_constants
    USE mo_veg_constants


    IMPLICIT NONE
    ! 0.1 InOut
    CLASS(t_parameters_from_file), INTENT(inout) :: this
    ! 0.2 Local
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':set_parameter_values_with_list_values'


    !> 1.0 read parameter values from file
    !! 
    CALL this%read_param

    !> 2.0 replace parameter values from mo_PROC_constants for parameters provided by param.list
    !! 
    ! SHARED_CTL
    IF (ABS(this%mavg_period_tlabile  - missval) > eps8)  mavg_period_tlabile   = this%mavg_period_tlabile
    IF (ABS(this%mavg_period_tphen    - missval) > eps8)  mavg_period_tphen     = this%mavg_period_tphen
    IF (ABS(this%mavg_period_tuptake  - missval) > eps8)  mavg_period_tuptake   = this%mavg_period_tuptake
    IF (ABS(this%mavg_period_tfrac    - missval) > eps8)  mavg_period_tfrac     = this%mavg_period_tfrac
    IF (ABS(this%mavg_period_tcnl     - missval) > eps8)  mavg_period_tcnl      = this%mavg_period_tcnl
    IF (ABS(this%mavg_period_talloc   - missval) > eps8)  mavg_period_talloc    = this%mavg_period_talloc
    IF (ABS(this%mavg_period_tvegdyn  - missval) > eps8)  mavg_period_tvegdyn   = this%mavg_period_tvegdyn
    IF (ABS(this%mavg_period_tacclim  - missval) > eps8)  mavg_period_tacclim   = this%mavg_period_tacclim
    ! VEGETATION_CTL
    IF (ABS(this%delta_n_fraction     - missval) > eps8)  delta_n_fraction      = this%delta_n_fraction
    IF (ABS(this%delta_n_leaf         - missval) > eps8)  delta_n_leaf          = this%delta_n_leaf
    IF (ABS(this%cost_n_fixation_max  - missval) > eps8)  cost_n_fixation_max   = this%cost_n_fixation_max / &
                                                                                    (molar_mass_C/molar_mass_N)
    IF (ABS(this%cost_n_fixation_min  - missval) > eps8)  cost_n_fixation_min   = this%cost_n_fixation_min / &
                                                                                    (molar_mass_C/molar_mass_N)
    IF (ABS(this%vmax_bnf             - missval) > eps8)  vmax_bnf              = this%vmax_bnf / &
                                                                                    (molar_mass_N/molar_mass_C) / &
                                                                                    (one_day*one_year) *1.e6_wp
    IF (ABS(this%k_bnf                - missval) > eps8)  k_bnf                 = this%k_bnf / &
                                                                                    (molar_mass_C/molar_Mass_N)
    IF (ABS(this%offset_bnf           - missval) > eps8)  offset_bnf            = this%offset_bnf * 2.0_wp * & 
                                                                                    (transform_cost_nh4*transform_cost_no3)/ &
                                                                                    (transform_cost_nh4+transform_cost_no3)
    IF (ABS(this%f_resp_acclim        - missval) > eps8)  f_resp_acclim         = this%f_resp_acclim
    IF (ABS(this%fmaint_rate_base     - missval) > eps8)  fmaint_rate_base      = this%fmaint_rate_base

    ! SB_CTL
    IF (ABS(this%floss_nmin           - missval) > eps8)  floss_nmin            = this%floss_nmin
    IF (ABS(this%tau_fast             - missval) > eps8)  tau_fast              = this%tau_fast    * one_year * one_day
    IF (ABS(this%tau_slow             - missval) > eps8)  tau_slow              = this%tau_slow    * one_year * one_day
    IF (ABS(this%microbial_cue        - missval) > eps8)  microbial_cue         = this%microbial_cue
    IF (ABS(this%residue_cue          - missval) > eps8)  residue_cue           = this%residue_cue

  END SUBROUTINE set_parameter_values_with_list_values


END MODULE mo_shared_set_parameters
