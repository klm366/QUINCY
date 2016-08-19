!>
!! @brief Contains structures and methods for JSBACH model config
!!
!! @author
!!  Reiner Schnur, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Reiner Schnur (2013-04-04)
!! integration into SLM by Jan Engel, MPI BGC (20160401)
!!
MODULE mo_jsb_config_class
!#ifndef __NO_JSBACH__

  USE mo_io_units,          ONLY: filename_max
!  USE mo_exception,         ONLY: message

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: t_jsb_config!, t_jsb_config_p, t_jsb_model_config, new_model_config

!  ! model config
!  TYPE t_jsb_model_config
!  END TYPE t_jsb_model_config

  !-----------------------------------------------------------------------------------------------------
  !> abstract TYPE for configuration that cannot be initialized
  !! 
  !! is used by other process configurations, e.g. veg and base \n
  !! copied from JSBACH4 git 20160815
  !-----------------------------------------------------------------------------------------------------
  TYPE, ABSTRACT :: t_jsb_config
    LOGICAL :: active
    CHARACTER(len=filename_max) :: ic_filename
    CHARACTER(len=filename_max) :: bc_filename
    CHARACTER(LEN=filename_max) :: namelist_filename
  CONTAINS
    PROCEDURE (Init_config), DEFERRED, PASS(config) :: Init
  END TYPE t_jsb_config

  TYPE t_jsb_config_p
    CLASS(t_jsb_config), POINTER :: p
  END TYPE t_jsb_config_p

  ! https://de.wikibooks.org/wiki/Fortran:_Fortran_2003:_Zeiger#Abstraktes_Interface
  ! defines that the TYPE itself is passed to the procedure (i.e. t_*_config -> init_config subroutine) 
  ABSTRACT INTERFACE
    SUBROUTINE Init_config(config)
      IMPORT :: t_jsb_config    ! makes t_jsb_config accessible within this (abstract) interface
      CLASS(t_jsb_config), INTENT(inout) :: config
    END SUBROUTINE Init_config
  END INTERFACE


!  CHARACTER(len=*), PARAMETER :: modname = 'mo_jsb_config_class'

!CONTAINS

!  ! model config
!  FUNCTION new_model_config(namelist_filename) RESULT(config)

!    USE mo_jsb_namelist_iface, ONLY: open_nml, POSITIONED, position_nml, close_nml

!    CHARACTER(len=*), INTENT(in) :: namelist_filename
!    TYPE(t_jsb_model_config), POINTER :: config !R faende es besser wenn die config Pointer schon im Namen unterschieden wuerden, macht es einfacher zu grepen...hier und in "Subroutine init_base" m_config

!    TYPE(t_jsb_model_config), POINTER :: new

!    CHARACTER(len=filename_max) :: grid_filename
!    CHARACTER(len=30)           :: grid_name

!    NAMELIST /jsb_grid_nml/   &
!      grid_filename,          &
!      grid_name

!    CHARACTER(len=30)           :: usecase

!    NAMELIST /jsb_model_nml/  &
!      usecase

!    INTEGER :: nml_handler, nml_unit, istat

!    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':new_model_config'

!    CALL message(TRIM(routine), 'Starting model configuration from '//TRIM(namelist_filename))

!    ALLOCATE(new)

!    ! Set defaults
!    grid_filename = ''
!    grid_name     = ''
!    usecase       = ''

!    ! Read namelist
!    nml_handler = open_nml(TRIM(namelist_filename))

!    nml_unit = position_nml('jsb_grid_nml', nml_handler, STATUS=istat)
!    IF (istat == POSITIONED) READ(nml_unit, jsb_grid_nml)

!    new%grid_filename = grid_filename
!    new%grid_name     = grid_name

!    nml_unit =  position_nml('jsb_model_nml', nml_handler, STATUS=istat)
!    IF (istat == POSITIONED) READ(nml_unit, jsb_model_nml)

!    new%usecase       = usecase

!    CALL close_nml(nml_handler)

!    config => new

!  END FUNCTION new_model_config

!#endif
END MODULE mo_jsb_config_class
