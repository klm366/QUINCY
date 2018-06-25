!>
!! @brief Contains structures and methods for soil_biogeochemistry config
!!
!! @author
!!  Reiner Schnur, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Reiner Schnur (2013-04-04);
!! 2nd version by Jan Engel, MPI BGC (20170104) -> integration into SLM
!!
MODULE mo_sb_config_class
!#ifndef __NO_JSBACH__

  USE mo_exception,         ONLY: message_text, message, finish
  USE mo_io_units,          ONLY: filename_max
  USE mo_kind,              ONLY: wp
  USE mo_jsb_config_class,  ONLY: t_jsb_config

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: t_sb_config

  !-----------------------------------------------------------------------------------------------------
  !> configuration of the soil_biogeochemistry process, derived from t_jsb_config
  !! 
  !! currently it does mainly: reading parameters from namelist
  !-----------------------------------------------------------------------------------------------------
  TYPE, EXTENDS(t_jsb_config) :: t_sb_config
    !> options
    !! 
    CHARACTER(15)         :: sb_model_scheme
    CHARACTER(15)         :: sb_nloss_scheme
    CHARACTER(15)         :: sb_bnf_scheme
    CHARACTER(15)         :: sb_adsorp_scheme
    !> flags
    !! 
    LOGICAL               :: flag_sb_prescribe_nh4
    LOGICAL               :: flag_sb_prescribe_no3
    LOGICAL               :: flag_sb_prescribe_po4
    ! mycorrhiza
    LOGICAL               :: flag_mycorrhiza
    LOGICAL               :: flag_myco_inorganic
    LOGICAL               :: flag_myco_organic
    LOGICAL               :: flag_myco_priming
    ! JSM
    LOGICAL               :: flag_sb_jsm_transport
    LOGICAL               :: flag_sb_jsm_litter_input
    LOGICAL               :: flag_sb_jsm_OM_sorption
    ! ECA
    LOGICAL               :: flag_sb_ECA_upt_po4
    LOGICAL               :: flag_sb_CLE_po4
    LOGICAL               :: flag_sb_LK_po4
    !> parameter
    !! 
    INTEGER               :: usda_taxonomy_class
    INTEGER               :: nwrb_taxonomy_class
    INTEGER               :: lithology_class
    REAL(wp)              :: soil_ph                  !< soil pH - site specific from forcing info file
    REAL(wp)              :: soil_p_depth
    REAL(wp)              :: soil_p_labile
    REAL(wp)              :: soil_p_slow
    REAL(wp)              :: soil_p_occluded
    REAL(wp)              :: soil_p_primary
    REAL(wp)              :: qmax_org_fine_particle    !< maximum sorption capacity of OM to fine soil particle, mol C/ kg fine particle
   CONTAINS
     PROCEDURE :: Init => Init_sb_config
  END type t_sb_config

  CHARACTER(len=*), PARAMETER :: modname = 'mo_sb_config_class'

CONTAINS

  !-----------------------------------------------------------------------------------------------------
  !> configuration routine of t_sb_config
  !! 
  !! currently it does only: read parameters from namelist
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE Init_sb_config(config)

    USE mo_namelist,           ONLY: open_nml, POSITIONED, position_nml, close_nml  ! read model options from file
    USE mo_base_config_class,  ONLY: base_config

    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    CLASS(t_sb_config), INTENT(inout) :: config    !< config type for sb
    ! ---------------------------
    ! 0.2 Local
    ! variables for reading from namlist, identical to variable-name in namelist
    LOGICAL                     :: active
    CHARACTER(len=filename_max) :: ic_filename, bc_filename
    CHARACTER(15)               :: sb_model_scheme
    CHARACTER(15)               :: sb_nloss_scheme
    CHARACTER(15)               :: sb_bnf_scheme
    CHARACTER(15)               :: sb_adsorp_scheme
    LOGICAL                     :: flag_sb_prescribe_nh4
    LOGICAL                     :: flag_sb_prescribe_no3
    LOGICAL                     :: flag_sb_prescribe_po4
    LOGICAL                     :: flag_mycorrhiza
    LOGICAL                     :: flag_myco_inorganic
    LOGICAL                     :: flag_myco_organic
    LOGICAL                     :: flag_myco_priming
    LOGICAL                     :: flag_sb_jsm_transport
    LOGICAL                     :: flag_sb_jsm_litter_input
    LOGICAL                     :: flag_sb_jsm_OM_sorption
    LOGICAL                     :: flag_sb_ECA_upt_po4
    LOGICAL                     :: flag_sb_CLE_po4
    LOGICAL                     :: flag_sb_LK_po4
    INTEGER                     :: usda_taxonomy_class
    INTEGER                     :: nwrb_taxonomy_class
    INTEGER                     :: lithology_class
    REAL(wp)                    :: soil_ph
    REAL(wp)                    :: soil_p_depth,        &
                                   soil_p_labile,       &
                                   soil_p_slow,         &
                                   soil_p_occluded,     &
                                   soil_p_primary
    REAL(wp)                    :: qmax_org_fine_particle

    NAMELIST /soil_biogeochemistry_ctl/       & 
                      active                    , &
                      ic_filename               , &
                      bc_filename               , &
                      sb_model_scheme           , &
                      sb_nloss_scheme           , &
                      sb_bnf_scheme             , &
                      sb_adsorp_scheme          , &
                      flag_sb_prescribe_nh4     , &
                      flag_sb_prescribe_no3     , &
                      flag_sb_prescribe_po4     , &
                      flag_mycorrhiza           , &
                      flag_myco_inorganic       , &
                      flag_myco_organic         , &
                      flag_myco_priming         , & 
                      flag_sb_jsm_transport     , &
                      flag_sb_jsm_litter_input  , &
                      flag_sb_jsm_OM_sorption   , &
                      flag_sb_ECA_upt_po4       , &
                      flag_sb_CLE_po4           , &
                      flag_sb_LK_po4            , &    
                      usda_taxonomy_class       , &
                      nwrb_taxonomy_class       , &
                      lithology_class           , &
                      soil_ph                   , &
                      soil_p_depth              , &
                      soil_p_labile             , &
                      soil_p_slow               , &
                      soil_p_occluded           , &
                      soil_p_primary            , &
                      qmax_org_fine_particle

    ! variables for reading model-options from namelist
    INTEGER  :: nml_handler, nml_unit, istat

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':Init_sb_config'

    CALL message(TRIM(routine), 'Starting soil_biogeochemistry configuration')

    ! Set defaults
    active                      = .TRUE.
    ic_filename                 = 'ic_land_srf.nc'
    bc_filename                 = 'bc_land_srf.nc'
    sb_model_scheme             = "simple_0d"       ! simple_0d simple_1d jsm
    sb_nloss_scheme             = "dynamic"         ! none fixed dynamic 
    sb_bnf_scheme               = "standard"        ! none standard unlimited 
    sb_adsorp_scheme            = "eca_full"        ! eca_full; eca_part (only uptake, not for adsorption); eca_none
    flag_sb_prescribe_nh4       = .FALSE.
    flag_sb_prescribe_no3       = .FALSE.
    flag_sb_prescribe_po4       = .FALSE.
    flag_mycorrhiza             = .FALSE.
    flag_myco_inorganic         = .FALSE.
    flag_myco_organic           = .FALSE.
    flag_myco_priming           = .FALSE.
    flag_sb_jsm_transport       = .FALSE.
    flag_sb_jsm_litter_input    = .FALSE.
    flag_sb_jsm_OM_sorption     = .FALSE.
    flag_sb_ECA_upt_po4         = .TRUE.
    flag_sb_CLE_po4             = .TRUE.
    flag_sb_LK_po4              = .TRUE.
    usda_taxonomy_class         = 1
    nwrb_taxonomy_class         = 1
    lithology_class             = 1
    soil_ph                     = 6.5_wp 
    soil_p_depth                = 0.5_wp        ! these soil_p* values are valid for DE-Hainich ?!
    soil_p_labile               = 50._wp
    soil_p_slow                 = 70._wp
    soil_p_occluded             = 150._wp
    soil_p_primary              = 100._wp
    qmax_org_fine_particle      = 6.5537_wp

    ! read the namelist
    ! &SOIL_BIOGEOCHEMISTRY_CTL
    nml_handler = open_nml(TRIM(config%namelist_filename))

    nml_unit = position_nml('SOIL_BIOGEOCHEMISTRY_CTL', nml_handler, STATUS=istat)
    IF (istat == POSITIONED) READ(nml_unit, soil_biogeochemistry_ctl)

    config%active                    = active
    config%ic_filename               = ic_filename
    config%bc_filename               = bc_filename

    CALL close_nml(nml_handler)

    SELECT CASE (TRIM(sb_model_scheme))
      CASE ("simple_0d", "simple_1d")
        flag_sb_ECA_upt_po4         = .TRUE.
        flag_sb_CLE_po4             = .TRUE.
        flag_sb_LK_po4              = .TRUE.
        sb_adsorp_scheme            = "eca_full"        
      CASE ("jsm")
         IF(TRIM(base_config%slm_model_name) /= "soil") THEN
            WRITE(message_text,'(3a)') 'setting sb_model_scheme to ',TRIM(sb_model_scheme), & 
                                       ' currently only works for the standalone soil model'
            CALL finish("mo_sb_config_class:init_sb_config",message_text) 
         ENDIF 
      CASE DEFAULT
        WRITE(message_text,'(2a)') 'invalid sb_model_scheme',sb_model_scheme
        CALL finish("mo_sb_config_class:init_sb_config",message_text) 
    ENDSELECT

    ! pass values as read from file
    config%sb_model_scheme                 = sb_model_scheme
    config%sb_nloss_scheme                 = sb_nloss_scheme
    config%sb_bnf_scheme                   = sb_bnf_scheme
    config%sb_adsorp_scheme                = sb_adsorp_scheme
    config%flag_sb_prescribe_nh4           = flag_sb_prescribe_nh4
    config%flag_sb_prescribe_no3           = flag_sb_prescribe_no3
    config%flag_sb_prescribe_po4           = flag_sb_prescribe_po4
    config%flag_mycorrhiza                 = flag_mycorrhiza
    config%flag_myco_inorganic             = flag_myco_inorganic
    config%flag_myco_organic               = flag_myco_organic
    config%flag_myco_priming               = flag_myco_priming
    config%flag_sb_jsm_transport           = flag_sb_jsm_transport
    config%flag_sb_jsm_litter_input        = flag_sb_jsm_litter_input
    config%flag_sb_jsm_OM_sorption         = flag_sb_jsm_OM_sorption
    config%flag_sb_ECA_upt_po4             = flag_sb_ECA_upt_po4
    config%flag_sb_CLE_po4                 = flag_sb_CLE_po4
    config%flag_sb_LK_po4                  = flag_sb_LK_po4
    config%usda_taxonomy_class             = usda_taxonomy_class
    config%nwrb_taxonomy_class             = nwrb_taxonomy_class
    config%lithology_class                 = lithology_class
    config%soil_ph                         = soil_ph 
    config%soil_p_depth                    = soil_p_depth
    config%soil_p_labile                   = soil_p_labile 
    config%soil_p_slow                     = soil_p_slow 
    config%soil_p_occluded                 = soil_p_occluded 
    config%soil_p_primary                  = soil_p_primary    
    config%qmax_org_fine_particle          = qmax_org_fine_particle

  END SUBROUTINE Init_sb_config

!#endif
END MODULE mo_sb_config_class
