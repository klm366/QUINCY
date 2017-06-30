!>
!! @brief Module to init, read and convert external forcing data
!!
!! @author
!!  Jan Engel, MPI-BGC Jena
!!
!! @par Revision History
!! First version by Jan Engel (2017-06-30): based on the git:JSBACH4 mo_jsb4_forcing module integrating only forcing_options and a basic version of Init_forcing()
!! 

MODULE mo_jsb4_forcing

#include "dsl_jsbach.inc"

  USE mo_kind

  IMPLICIT NONE


  PRIVATE
  PUBLIC :: Init_forcing
  PUBLIC :: forcing_options

  !-----------------------------------------------------------------------------------------------------
  !> Namelist parameters
  !! 
  !-----------------------------------------------------------------------------------------------------
  LOGICAL       :: flag_n_deposition


  !-----------------------------------------------------------------------------------------------------
  !> Structure for handling input data
  !! 
  !-----------------------------------------------------------------------------------------------------
  TYPE forcing_options_type
    LOGICAL       :: flag_n_deposition

  END TYPE forcing_options_type
  TYPE(forcing_options_type), SAVE :: forcing_options




  CHARACTER(len=*), PARAMETER :: modname = 'mo_jsb4_forcing'

CONTAINS

  !-----------------------------------------------------------------------------------------------------
  !> configuration routine of t_sb_config
  !! 
  !! currently it does only: read parameters from namelist
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE Init_forcing(model_id)

    USE mo_namelist,                        ONLY: open_nml, POSITIONED, position_nml, close_nml  ! read model options from file
    USE mo_exception,                       ONLY: message_text, message, finish
    USE mo_jsb_class,                       ONLY: Get_model
    USE mo_jsb_model_class,                 ONLY: t_jsb_model

    IMPLICIT NONE
    ! 0.1 InOut
    INTEGER, INTENT(in)                 :: model_id
    ! 0.2 Local
    TYPE(t_jsb_model), POINTER          :: model
    CHARACTER(len=*), PARAMETER         :: routine = TRIM(modname)//':Init_forcing'

    NAMELIST /jsb_forcing_ctl/       & 
                      flag_n_deposition

    ! variables for reading model-options from namelist
    INTEGER  :: nml_handler, nml_unit, istat

    CALL message(TRIM(routine), 'Starting Init_forcing')

    ! get model and namelist_filename
    model => Get_model(model_id)

    !> 1.0 Read namelist
    !! 

    ! Set defaults
    flag_n_deposition           = .FALSE.

    ! read the namelist
    ! 
    ! &JSB_FORCING_CTL
    nml_handler = open_nml(TRIM(model%namelist_filename))
    nml_unit = position_nml('JSB_FORCING_CTL', nml_handler, STATUS=istat)
    IF (istat == POSITIONED) READ(nml_unit, jsb_forcing_ctl)
    CALL close_nml(nml_handler)


    ! pass values as read from file
    forcing_options%flag_n_deposition       = flag_n_deposition


    !> 2.0 Modify .. Configure ..
    !! 




    !> 3.0 Init forcing input 
    !! 



  END SUBROUTINE Init_forcing
  
END MODULE mo_jsb4_forcing
