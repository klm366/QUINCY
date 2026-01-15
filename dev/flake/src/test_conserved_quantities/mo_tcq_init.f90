!> Summary: tcq (test conserved quantities) memory initialisation
!>
!>### Initialization of tcq (test conserved quantities) variables
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
MODULE mo_tcq_init
#ifndef __NO_JSBACH__

  USE mo_kind,              ONLY: wp
  USE mo_jsb_tile_class,    ONLY: t_jsb_tile_abstract

  dsl4jsb_Use_processes TCQ_
  dsl4jsb_Use_memory(TCQ_)

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: tcq_init

  CHARACTER(len=*), PARAMETER :: modname = 'mo_tcq_init'

CONTAINS

  ! ====================================================================================================== !
  !
  !> Initialize tcq process (after memory has been set up)
  !
  SUBROUTINE tcq_init(tile)

    ! -------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    ! -------------------------------------------------------------------------------------------------- !
    CHARACTER(len=*), PARAMETER :: routine = modname//':tcq_init'

    dsl4jsb_Def_memory(TCQ_)

    dsl4jsb_Real2D_onDomain :: a_dead_c_tcq
    dsl4jsb_Real2D_onDomain :: a_veg_c_tcq
    dsl4jsb_Real2D_onDomain :: another_dead_c_tcq
    dsl4jsb_Real2D_onDomain :: an_implicit_scaling_tcq
    ! -------------------------------------------------------------------------------------------------- !

    dsl4jsb_Get_memory(TCQ_)

    dsl4jsb_Get_var2D_onDomain(TCQ_, a_dead_c_tcq)
    dsl4jsb_Get_var2D_onDomain(TCQ_, a_veg_c_tcq)
    dsl4jsb_Get_var2D_onDomain(TCQ_, another_dead_c_tcq)
    dsl4jsb_Get_var2D_onDomain(TCQ_, an_implicit_scaling_tcq)

    !>
    !> initialise with a constant depending on lct
    !>
    an_implicit_scaling_tcq(:,:) = 10.0_wp * tile%lcts(1)%lib_id
    a_dead_c_tcq(:,:) = tile%lcts(1)%lib_id
    a_veg_c_tcq(:,:) = 5.0_wp * tile%lcts(1)%lib_id
    another_dead_c_tcq(:,:) = 100.0_wp * tile%lcts(1)%lib_id

  END SUBROUTINE tcq_init

#endif
END MODULE mo_tcq_init
