!>
!! @brief Initialization of the the lentic memory.
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

MODULE mo_lentic_init
#ifndef __NO_JSBACH__
  
  USE mo_kind,              ONLY: wp
  USE mo_exception,         ONLY: message
  USE mo_jsb_control,       ONLY: debug_on

  !USE mo_jsb_model_class,   ONLY: t_jsb_model
  USE mo_jsb_tile_class,    ONLY: t_jsb_tile_abstract
  !USE mo_jsb_class,         ONLY: get_model
  !USE mo_jsb_io_netcdf,     ONLY: t_input_file !, jsb_netcdf_open_input

  dsl4jsb_Use_processes LENTIC_
  dsl4jsb_Use_config(LENTIC_)
  dsl4jsb_Use_memory(LENTIC_)

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: lentic_init

  CHARACTER(len=*), PARAMETER :: modname = 'mo_lentic_init'

CONTAINS

  !
  !> Intialize lentic process (after memory has been set up)
  !
  SUBROUTINE lentic_init(tile)
    USE mo_jsb_time,         ONLY: timestep_in_days
    USE mo_lentic_constants, ONLY : tpl_T_r, &
                                  & tpl_T_f
    

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile

    CHARACTER(len=*), PARAMETER :: routine = modname//':lentic_init'

    dsl4jsb_Real2D_onDomain :: &
      & c_t   , &
      & h_ice , &
      & h_ML  , &
      & h_snow, &
      & T_bot , &
      & T_lice , &
      & T_mnw , &
      & T_sfc , &
      !& T_snow, &
      & T_wML


    dsl4jsb_Def_memory(LENTIC_)

    dsl4jsb_Get_memory(LENTIC_)


    IF (tile%contains_lake) THEN

      dsl4jsb_Get_var2D_onDomain(LENTIC_, c_t) ! OUT
      dsl4jsb_Get_var2D_onDomain(LENTIC_, h_ice) ! OUT
      dsl4jsb_Get_var2D_onDomain(LENTIC_, h_ML) ! OUT
      dsl4jsb_Get_var2D_onDomain(LENTIC_, h_snow) ! OUT
      dsl4jsb_Get_var2D_onDomain(LENTIC_, T_bot) ! OUT
      dsl4jsb_Get_var2D_onDomain(LENTIC_, T_lice) ! OUT
      dsl4jsb_Get_var2D_onDomain(LENTIC_, T_mnw) ! OUT
      dsl4jsb_Get_var2D_onDomain(LENTIC_, T_sfc) ! OUT
      ! dsl4jsb_Get_var2D_onDomain(LENTIC_, T_snow) ! OUT
      dsl4jsb_Get_var2D_onDomain(LENTIC_, T_wML) ! OUT

      c_t         = 0.5_wp 
      h_ice       = 1_wp
      h_ML        = 1_wp
      h_snow      = 0
      T_bot       = tpl_T_r
      T_lice       = tpl_T_f
      T_mnw       = tpl_T_f
      T_sfc       = tpl_T_f
      ! T_snow      = tpl_T_f
      T_wML       = tpl_T_f

      !$acc update device(F_pseudo_soil_temp, N_pseudo_soil_temp)
    END IF

  END SUBROUTINE lentic_init
  
#endif
END MODULE mo_lentic_init
  