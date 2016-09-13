!>
!! @brief Contains structures and methods for soil config
!!
!! @author
!!  Reiner Schnur, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Reiner Schnur (2013-04-04);
!! 2nd version by Jan Engel, MPI BGC (20160913) -> integration into SLM
!!
MODULE mo_soil_config_class
!#ifndef __NO_JSBACH__

  USE mo_kind,              ONLY: wp
  USE mo_exception,         ONLY: message_text, message, finish
  USE mo_io_units,          ONLY: filename_max
  USE mo_jsb_config_class,  ONLY: t_jsb_config

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: soil_config

  !-----------------------------------------------------------------------------------------------------
  !> configuration of the soil process, derived from t_jsb_config
  !! 
  !! currently it does mainly: reading parameters from namelist
  !-----------------------------------------------------------------------------------------------------
  TYPE, EXTENDS(t_jsb_config) :: t_soil_config
! 
! 

    CONTAINS
      PROCEDURE :: Init => Init_soil_config
  END Type t_soil_config

  TYPE(t_soil_config), SAVE :: soil_config

  CHARACTER(len=*), PARAMETER :: modname = 'mo_soil_config_class'

CONTAINS

  !-----------------------------------------------------------------------------------------------------
  !> configuration routine of t_soil_config
  !! 
  !! currently it does only: read parameters from namelist
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE Init_soil_config(config)

    USE mo_namelist,           ONLY: open_nml, POSITIONED, position_nml, close_nml  ! read model options from file
    ! 0.1 InOut
    CLASS(t_soil_config), INTENT(inout) :: config    !< config type for soil
    ! 0.2 Local
    ! variables for reading from namlist, identical to variable-name in namelist
    !LOGICAL      :: active
! 
! 


!    NAMELIST /soil_ctl/       & 
!      ,         &
!      

    ! variables for reading model-options from namelist
    INTEGER  :: nml_handler, nml_unit, istat

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':Init_soil_config'

    CALL message(TRIM(routine), 'Starting soil configuration')

    ! Set defaults
! 
! 


    ! read the namelist
    ! 
    ! &SOIL_CTL
!    nml_handler = open_nml(TRIM(config%namelist_filename))
!    nml_unit = position_nml('SOIL_CTL', nml_handler, STATUS=istat)
!    IF (istat == POSITIONED) READ(nml_unit, soil_ctl)
!    CALL close_nml(nml_handler)


    ! pass values as read from file
! veg examples: replace !   config%output_interval_pool          = output_interval_pool
!    config%output_interval_flux          = output_interval_flux
!    config%output_start_first_day_year   = output_start_first_day_year
!    config%output_end_last_day_year      = output_end_last_day_year
!    config%flag_calc_soil_moisture       = flag_calc_soil_moisture
!    config%flag_tile_veg_dynamics        = flag_tile_veg_dynamics  
!    config%flag_start_bareground         = flag_start_bareground


  END SUBROUTINE Init_soil_config

!#endif
END MODULE mo_soil_config_class
