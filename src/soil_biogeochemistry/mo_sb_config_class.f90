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
    CHARACTER(15)         :: sb_model_scheme
    CHARACTER(15)         :: sb_nloss_scheme
    CHARACTER(15)         :: sb_bnf_scheme
    LOGICAL               :: flag_sb_prescribe_nh4
    LOGICAL               :: flag_sb_prescribe_no3
    LOGICAL               :: flag_sb_prescribe_po4
    LOGICAL               :: flag_mycorrhiza
    LOGICAL               :: flag_myco_inorganic
    LOGICAL               :: flag_myco_organic
    LOGICAL               :: flag_myco_priming
    LOGICAL               :: flag_sb_ECA_upt_po4
    LOGICAL               :: flag_sb_CLE_po4
    LOGICAL               :: flag_sb_LK_po4
    CHARACTER(15)         :: sb_adsorp_scheme
    INTEGER               :: usda_taxonomy_class
    INTEGER               :: nwrb_taxonomy_class
    REAL(wp)              :: soil_ph                  !< soil pH - site specific from forcing info file
    INTEGER               :: lithology_class
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

    ! 0.1 InOut
    CLASS(t_sb_config), INTENT(inout) :: config    !< config type for sb

    ! 0.2 Local
    ! variables for reading from namlist, identical to variable-name in namelist
    LOGICAL  :: active
    CHARACTER(len=filename_max) :: ic_filename, bc_filename
    CHARACTER(15)         :: sb_model_scheme
    CHARACTER(15)         :: sb_nloss_scheme
    CHARACTER(15)         :: sb_bnf_scheme
    LOGICAL               :: flag_sb_prescribe_nh4
    LOGICAL               :: flag_sb_prescribe_no3
    LOGICAL               :: flag_sb_prescribe_po4
    LOGICAL               :: flag_mycorrhiza
    LOGICAL               :: flag_myco_inorganic
    LOGICAL               :: flag_myco_organic
    LOGICAL               :: flag_myco_priming
    LOGICAL               :: flag_sb_ECA_upt_po4
    LOGICAL               :: flag_sb_CLE_po4
    LOGICAL               :: flag_sb_LK_po4
    CHARACTER(15)         :: sb_adsorp_scheme
    INTEGER               :: usda_taxonomy_class
    INTEGER               :: nwrb_taxonomy_class
    REAL(wp)              :: soil_ph
    INTEGER               :: lithology_class

    NAMELIST /soil_biogeochemistry_ctl/       & 
      active,                      &
      ic_filename, bc_filename,    &
                      sb_model_scheme      ,  &
                      sb_nloss_scheme      ,  &
                      sb_bnf_scheme        ,  &
                      flag_sb_prescribe_nh4,  &
                      flag_sb_prescribe_no3,  &
                      flag_sb_prescribe_po4,  &
                      flag_mycorrhiza      ,  &
                      flag_myco_inorganic  ,  &
                      flag_myco_organic    ,  &
                      flag_myco_priming    ,  & 
                      flag_sb_ECA_upt_po4  ,  &
                      flag_sb_CLE_po4      ,  &
                      flag_sb_LK_po4       ,  &    
                      sb_adsorp_scheme     ,  &
                      usda_taxonomy_class  ,  &
                      nwrb_taxonomy_class  ,  &
                      soil_ph              ,  &
                      lithology_class

    ! variables for reading model-options from namelist
    INTEGER  :: nml_handler, nml_unit, istat

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':Init_sb_config'

    CALL message(TRIM(routine), 'Starting soil_biogeochemistry configuration')

    ! Set defaults
    active          = .TRUE.
    ic_filename     = 'ic_land_srf.nc'
    bc_filename     = 'bc_land_srf.nc'
    sb_model_scheme             = "simple_0d"       ! simple_0d simple_1d jsm
    sb_nloss_scheme             = "fixed"           ! none fixed dynamic 
    sb_bnf_scheme               = "standard"        ! none standard unlimited 
    flag_sb_prescribe_nh4       = .FALSE.
    flag_sb_prescribe_no3       = .FALSE.
    flag_sb_prescribe_po4       = .FALSE.
    flag_mycorrhiza             = .FALSE.
    flag_myco_inorganic         = .FALSE.
    flag_myco_organic           = .FALSE.
    flag_myco_priming           = .FALSE.
    flag_sb_ECA_upt_po4         = .TRUE.
    flag_sb_CLE_po4             = .TRUE.
    flag_sb_LK_po4              = .TRUE.
    sb_adsorp_scheme            = "eca_full"         ! eca_full; eca_part (only uptake, not for adsorption); eca_none
    usda_taxonomy_class         = 1
    nwrb_taxonomy_class         = 1
    soil_ph                     = 6.5_wp
    lithology_class             = 1

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
    config%flag_sb_prescribe_nh4           = flag_sb_prescribe_nh4
    config%flag_sb_prescribe_no3           = flag_sb_prescribe_no3
    config%flag_sb_prescribe_po4           = flag_sb_prescribe_po4
    config%flag_mycorrhiza                 = flag_mycorrhiza
    config%flag_myco_inorganic             = flag_myco_inorganic
    config%flag_myco_organic               = flag_myco_organic
    config%flag_myco_priming               = flag_myco_priming
    config%flag_sb_ECA_upt_po4             = flag_sb_ECA_upt_po4
    config%flag_sb_CLE_po4                 = flag_sb_CLE_po4
    config%flag_sb_LK_po4                  = flag_sb_LK_po4
    config%sb_adsorp_scheme                = sb_adsorp_scheme
    config%usda_taxonomy_class             = usda_taxonomy_class
    config%nwrb_taxonomy_class             = nwrb_taxonomy_class
    config%soil_ph                         = soil_ph    
    config%lithology_class                 = lithology_class

  END SUBROUTINE Init_sb_config

!#endif
END MODULE mo_sb_config_class
