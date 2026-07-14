!>
!! @brief Initialization of the the flake memory.
!!
!! @author
!!  Reiner Schnur, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Reiner Schnur (2013-04-04)
!!
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/license .
!!

MODULE mo_flake_init
#ifndef __NO_JSBACH__
  
  USE mo_kind,              ONLY: wp
  USE mo_exception,         ONLY: message
  USE mo_jsb_control,       ONLY: debug_on

  USE mo_jsb_model_class,   ONLY: t_jsb_model
  USE mo_jsb_tile_class,    ONLY: t_jsb_tile_abstract
  USE mo_jsb_class,         ONLY: get_model
  !USE mo_jsb_io_netcdf,     ONLY: t_input_file !, jsb_netcdf_open_input

  dsl4jsb_Use_processes FLAKE_
  dsl4jsb_Use_config(FLAKE_)
  dsl4jsb_Use_memory(FLAKE_)

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: flake_init

  CHARACTER(len=*), PARAMETER :: modname = 'mo_flake_init'

CONTAINS

  !
  !> Intialize flake process (after memory has been set up)
  !
  SUBROUTINE flake_init(tile)
    USE mo_jsb_time,         ONLY: timestep_in_days
    USE mo_flake_constants, ONLY : tpl_T_r, &
                                  & tpl_T_f
    
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile

    CHARACTER(len=*), PARAMETER :: routine = modname//':flake_init'
    
    TYPE(t_jsb_model), POINTER :: model

    dsl4jsb_Real2D_onDomain :: &
      & c_t   , &
      & h_ice , &
      & h_ML  , &
      & h_snow, &
      & T_bot , &
      & T_lice , &
      & T_mnw , &
      !& T_sfc , &
      !& T_snow, &
      & T_wML, &
      & fetch

    dsl4jsb_Def_memory(FLAKE_)
    dsl4jsb_Def_config(FLAKE_)

    model => get_model(tile%owner_model_id)

    dsl4jsb_Get_memory(FLAKE_)
    dsl4jsb_Get_config(FLAKE_)

    IF (tile%contains_lake) THEN

      IF (.NOT. model%config%is_single_site) THEN
        CALL message(routine, 'WARNING: Initialization of the Flake process is currently using constant values from namelist')
      END IF

      dsl4jsb_Get_var2D_onDomain(FLAKE_, c_t)    ! OUT
      dsl4jsb_Get_var2D_onDomain(FLAKE_, h_ice)  ! OUT
      dsl4jsb_Get_var2D_onDomain(FLAKE_, h_ML)   ! OUT
      dsl4jsb_Get_var2D_onDomain(FLAKE_, h_snow) ! OUT
      dsl4jsb_Get_var2D_onDomain(FLAKE_, T_bot)  ! OUT
      dsl4jsb_Get_var2D_onDomain(FLAKE_, T_lice) ! OUT
      dsl4jsb_Get_var2D_onDomain(FLAKE_, T_mnw)  ! OUT
      ! dsl4jsb_Get_var2D_onDomain(FLAKE_, T_sfc)  ! OUT
      ! dsl4jsb_Get_var2D_onDomain(FLAKE_, T_snow) ! OUT
      dsl4jsb_Get_var2D_onDomain(FLAKE_, T_wML)  ! OUT
      dsl4jsb_Get_var2D_onDomain(FLAKE_, fetch)  ! OUT

      c_t(:,:)    = dsl4jsb_Config(FLAKE_)%c_t
      h_ice(:,:)  = dsl4jsb_Config(FLAKE_)%h_ice
      h_ML(:,:)   = dsl4jsb_Config(FLAKE_)%h_ML
      h_snow(:,:) = dsl4jsb_Config(FLAKE_)%h_snow
      T_bot(:,:)  = dsl4jsb_Config(FLAKE_)%T_bot
      T_lice(:,:) = dsl4jsb_Config(FLAKE_)%T_lice
      T_mnw(:,:)  = dsl4jsb_Config(FLAKE_)%T_mnw
      ! T_sfc       = dsl4jsb_Config(FLAKE_, T_sfc)
      ! T_snow      = dsl4jsb_Config(FLAKE_, T_snow)
      T_wML(:,:)  = dsl4jsb_Config(FLAKE_)%T_wML
      fetch(:,:)  = dsl4jsb_Config(FLAKE_)%fetch

    END IF

  END SUBROUTINE flake_init
  
#endif
END MODULE mo_flake_init
  