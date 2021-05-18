!>
!!  @par Copyright
!!    QUINCY terrestrial biosphere model with coupled carbon, nitrogen, and phosphorus cycles
!!    Copyright (C) 2019 Max Planck Institute for Biogeochemistry (MPI-BGC), Jena, Germany 
!!    Contributors: Caldararu S., Engel J., Kern M., Pallandt M., Schnur R., Thum T., Yu L., Zaehle S.
!!
!!    This file is part of the QUINCY model.
!!    The QUINCY model is free software: you can redistribute it and/or modify
!!    it under the terms of the GNU General Public License version 3,
!!    as published by the Free Software Foundation.
!!
!!    The QUINCY model is distributed in the hope that it will be useful,
!!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!!    GNU General Public License for more details.
!!
!!    You should have received a copy of the GNU General Public License
!!    along with the QUINCY model.  If not, see <https://www.gnu.org/licenses/>.
!! 
!! @brief Contains structures and methods for soil_biogeochemistry config
!!
!! @author
!!  Reiner Schnur, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Reiner Schnur (2013-04-04);
!! 2nd version by Jan Engel, MPI BGC (20170104) -> integration into QUINCY
!! 3rd version by Jan Engel, MPI BGC (20200619) -> integration into jsbach4
!!
MODULE mo_sb_config_class

  USE mo_kind,              ONLY: wp
  USE mo_exception,         ONLY: message_text, message, finish
  USE mo_io_units,          ONLY: filename_max
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
    LOGICAL               :: flag_mycorrhiza_org
    LOGICAL               :: flag_mycorrhiza_prim
    ! JSM
    LOGICAL               :: flag_sb_jsm_transport
    LOGICAL               :: flag_sb_jsm_litter_input
    LOGICAL               :: flag_sb_jsm_OM_sorption
    ! P sorpion
    LOGICAL               :: flag_sb_double_langmuir
    !> parameter
    !! 
    INTEGER               :: usda_taxonomy_class
    INTEGER               :: nwrb_taxonomy_class
    INTEGER               :: lithology_class
    REAL(wp)              :: soil_ph                  !< soil pH - site specific from forcing info file
    REAL(wp)              :: soil_p_depth             !< ... ?
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

    USE mo_jsb_namelist_iface, ONLY: open_nml, POSITIONED, position_nml, close_nml  ! read model options from file

    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    CLASS(t_sb_config), INTENT(inout) :: config    !< config type for sb
    ! ---------------------------
    ! 0.2 Local
    ! variables for reading from namlist, identical to variable-name in namelist
    LOGICAL                     :: active
    CHARACTER(len=filename_max) :: bc_filename, &
                                   ic_filename
    CHARACTER(15)               :: sb_model_scheme
    CHARACTER(15)               :: sb_nloss_scheme
    CHARACTER(15)               :: sb_bnf_scheme
    CHARACTER(15)               :: sb_adsorp_scheme
    LOGICAL                     :: flag_sb_prescribe_nh4
    LOGICAL                     :: flag_sb_prescribe_no3
    LOGICAL                     :: flag_sb_prescribe_po4
    LOGICAL                     :: flag_mycorrhiza
    LOGICAL                     :: flag_mycorrhiza_org
    LOGICAL                     :: flag_mycorrhiza_prim
    LOGICAL                     :: flag_sb_jsm_transport
    LOGICAL                     :: flag_sb_jsm_litter_input
    LOGICAL                     :: flag_sb_jsm_OM_sorption
    LOGICAL                     :: flag_sb_double_langmuir
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
                      bc_filename               , &
                      ic_filename               , &
                      sb_model_scheme           , &
                      sb_nloss_scheme           , &
                      sb_bnf_scheme             , &
                      sb_adsorp_scheme          , &
                      flag_sb_prescribe_nh4     , &
                      flag_sb_prescribe_no3     , &
                      flag_sb_prescribe_po4     , &
                      flag_mycorrhiza           , &
                      flag_mycorrhiza_org       , &
                      flag_mycorrhiza_prim      , & 
                      flag_sb_jsm_transport     , &
                      flag_sb_jsm_litter_input  , &
                      flag_sb_jsm_OM_sorption   , &
                      flag_sb_double_langmuir   , &
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
    bc_filename                 = 'bc_land_sb.nc'
    ic_filename                 = 'ic_land_sb.nc'
    sb_model_scheme             = "simple_0d"       ! simple_0d simple_1d jsm
    sb_nloss_scheme             = "dynamic"         ! none fixed dynamic 
    sb_bnf_scheme               = "standard"        ! none standard unlimited 
    sb_adsorp_scheme            = "eca_full"        ! eca_full; eca_part (only uptake, not for adsorption); eca_none
    flag_sb_prescribe_nh4       = .FALSE.
    flag_sb_prescribe_no3       = .FALSE.
    flag_sb_prescribe_po4       = .FALSE.
    flag_mycorrhiza             = .FALSE.
    flag_mycorrhiza_org         = .FALSE.
    flag_mycorrhiza_prim        = .FALSE.
    flag_sb_jsm_transport       = .FALSE.
    flag_sb_jsm_litter_input    = .FALSE.
    flag_sb_jsm_OM_sorption     = .FALSE.
    flag_sb_double_langmuir     = .FALSE.
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

    ! ! read the namelist
    ! ! &SOIL_BIOGEOCHEMISTRY_CTL
    ! nml_handler = open_nml(TRIM(config%namelist_filename))
    ! 
    ! nml_unit = position_nml('SOIL_BIOGEOCHEMISTRY_CTL', nml_handler, STATUS=istat)
    ! IF (istat == POSITIONED) READ(nml_unit, soil_biogeochemistry_ctl)
    ! 
    ! CALL close_nml(nml_handler)

    ! pass values as read from file
    config%active                          = active
    config%ic_filename                     = ic_filename
    config%bc_filename                     = bc_filename
    config%sb_model_scheme                 = sb_model_scheme
    config%sb_nloss_scheme                 = sb_nloss_scheme
    config%sb_bnf_scheme                   = sb_bnf_scheme
    config%sb_adsorp_scheme                = sb_adsorp_scheme
    config%flag_sb_prescribe_nh4           = flag_sb_prescribe_nh4
    config%flag_sb_prescribe_no3           = flag_sb_prescribe_no3
    config%flag_sb_prescribe_po4           = flag_sb_prescribe_po4
    config%flag_mycorrhiza                 = flag_mycorrhiza
    config%flag_mycorrhiza_org             = flag_mycorrhiza_org
    config%flag_mycorrhiza_prim            = flag_mycorrhiza_prim
    config%flag_sb_jsm_transport           = flag_sb_jsm_transport
    config%flag_sb_jsm_litter_input        = flag_sb_jsm_litter_input
    config%flag_sb_jsm_OM_sorption         = flag_sb_jsm_OM_sorption
    config%flag_sb_double_langmuir         = flag_sb_double_langmuir
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

END MODULE mo_sb_config_class
