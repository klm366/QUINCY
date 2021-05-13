!>
!! @brief Initialization of the cq memory.
!!
!! @author
!!  Julia Nabel, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Julia Nabel (2021-03-09)
!!
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/license .
!!

MODULE mo_cq_init
#ifndef __NO_JSBACH__

  USE mo_kind,              ONLY: wp
  USE mo_exception,         ONLY: message

  USE mo_jsb_model_class,   ONLY: t_jsb_model
  USE mo_jsb_class,         ONLY: get_model
  USE mo_jsb_tile_class,    ONLY: t_jsb_tile_abstract
  USE mo_jsb_io_netcdf,     ONLY: t_input_file, jsb_netcdf_open_input
  USE mo_jsb_time,          ONLY: get_time_dt
  USE mo_jsb_control,       ONLY: debug_on

  dsl4jsb_Use_processes CQ_
  dsl4jsb_Use_memory(CQ_)

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: cq_init

  CHARACTER(len=*), PARAMETER :: modname = 'mo_cq_init'

CONTAINS

  !
  !> Intialize cq process (after memory has been set up)
  !
  SUBROUTINE cq_init(tile)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile

    CHARACTER(len=*), PARAMETER :: routine = modname//':cq_init'

    dsl4jsb_Def_memory(CQ_)

    dsl4jsb_Real2D_onDomain :: a_dead_c_cq
    dsl4jsb_Real2D_onDomain :: a_veg_c_cq
    dsl4jsb_Real2D_onDomain :: another_dead_c_cq
    dsl4jsb_Real2D_onDomain :: an_implicit_scaling_cq

    dsl4jsb_Get_memory(CQ_)

    dsl4jsb_Get_var2D_onDomain(CQ_, a_dead_c_cq)
    dsl4jsb_Get_var2D_onDomain(CQ_, a_veg_c_cq)
    dsl4jsb_Get_var2D_onDomain(CQ_, another_dead_c_cq)
    dsl4jsb_Get_var2D_onDomain(CQ_, an_implicit_scaling_cq)

    !initialise with a constant depending on lct
    an_implicit_scaling_cq(:,:) = 10.0_wp * tile%lcts(1)%lib_id
    a_dead_c_cq(:,:) = tile%lcts(1)%lib_id
    a_veg_c_cq(:,:) = 5.0_wp * tile%lcts(1)%lib_id
    another_dead_c_cq(:,:) = 100.0_wp * tile%lcts(1)%lib_id

  END SUBROUTINE cq_init

#endif
END MODULE mo_cq_init
