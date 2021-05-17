!> Summary: cq (test conserved quantities) memory initialisation
!>
!>### Initialization of cq (test conserved quantities) variables
!>
!>
!> Copyright (c) 2021 Julia Nabel
!>
!>#### Revision History
!> First version and ongoing changes   by Julia Nabel (2021-03-09)
!>
!>#### License
!>
!> This code is subject to the MPI-M-ICON - License - Agreement in it's most recent form.
!> Please see [this link](https://code.mpimet.mpg.de/attachments/download/20888/MPI-M-ICONLizenzvertragV2.6.pdf) for details.
!> Where software is supplied by third parties, it is indicated in the headers of the routines.
!>
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

  ! ====================================================================================================== !
  !
  !> Initialize cq process (after memory has been set up)
  !
  SUBROUTINE cq_init(tile)

    ! -------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    ! -------------------------------------------------------------------------------------------------- !
    CHARACTER(len=*), PARAMETER :: routine = modname//':cq_init'

    dsl4jsb_Def_memory(CQ_)

    dsl4jsb_Real2D_onDomain :: a_dead_c_cq
    dsl4jsb_Real2D_onDomain :: a_veg_c_cq
    dsl4jsb_Real2D_onDomain :: another_dead_c_cq
    dsl4jsb_Real2D_onDomain :: an_implicit_scaling_cq
    ! -------------------------------------------------------------------------------------------------- !

    dsl4jsb_Get_memory(CQ_)

    dsl4jsb_Get_var2D_onDomain(CQ_, a_dead_c_cq)
    dsl4jsb_Get_var2D_onDomain(CQ_, a_veg_c_cq)
    dsl4jsb_Get_var2D_onDomain(CQ_, another_dead_c_cq)
    dsl4jsb_Get_var2D_onDomain(CQ_, an_implicit_scaling_cq)

    !>
    !> initialise with a constant depending on lct
    !>
    an_implicit_scaling_cq(:,:) = 10.0_wp * tile%lcts(1)%lib_id
    a_dead_c_cq(:,:) = tile%lcts(1)%lib_id
    a_veg_c_cq(:,:) = 5.0_wp * tile%lcts(1)%lib_id
    another_dead_c_cq(:,:) = 100.0_wp * tile%lcts(1)%lib_id

  END SUBROUTINE cq_init

#endif
END MODULE mo_cq_init
