!>
!! @brief Contains methods for initialization of JSBACH usecases.
!!
!! @author
!!  Reiner Schnur, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Reiner Schnur (2013-04-04)
!!
MODULE mo_jsb_model_usecases
#ifndef __NO_JSBACH__

  USE mo_jsb_model_class, ONLY: t_jsb_model
  USE mo_exception,       ONLY: message, finish
  USE mo_kind,            ONLY: wp

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: init_usecase

  CHARACTER(len=*), PARAMETER :: modname = 'mo_jsb_model_usecases'

CONTAINS

  SUBROUTINE init_usecase(model)

#ifdef HAVE_F2003
    TYPE(t_jsb_model), POINTER, INTENT(inout) :: model
#else
    TYPE(t_jsb_model), POINTER                :: model
#endif

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':init_usecase'

    CALL message(TRIM(routine), 'Initializing JSBACH usecase "'//TRIM(model%config%usecase)//'"')

    SELECT CASE (TRIM(model%config%usecase))

    CASE ('jsbach_lite')
      CALL init_usecase_lite(model)

    CASE ('jsbach_lite_tpe_open')
      CALL init_usecase_lite_tpe_open(model)

    CASE ('jsbach_lite_tpe_closed')
      CALL init_usecase_lite_tpe_closed(model)

    CASE DEFAULT
      CALL finish(TRIM(routine), 'JSBACH usecase not defined.')

    END SELECT

  END SUBROUTINE init_usecase

  SUBROUTINE init_usecase_lite(model)

    USE mo_jsb_tile_class,     ONLY: t_jsb_tile, new_tile
    USE mo_jsb_grid_class,     ONLY: get_nblks, get_nproma
    USE mo_jsb_memory,         ONLY: new_jsbach_memory

    CHARACTER(len=20) :: components(3)

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':init_usecase_lite'

#ifdef HAVE_F2003
    TYPE(t_jsb_model), POINTER, INTENT(inout) :: model
#else
    TYPE(t_jsb_model), POINTER                :: model
#endif

    TYPE(t_jsb_tile), POINTER :: tile
    INTEGER :: nblks, nproma

    nblks = get_nblks(model%grid)
    nproma = get_nproma(model%grid)

    ! Activate soil component
    model%config%soil%active = .TRUE.

    ! Create root of the tile structure representing the whole grid box
    components(1) = 'srf' ; components(2) = 'soil' ; components(3) = 'veg'
    model%land => new_tile(model%id, 'land', 'land', components)

    model%land%mem => new_jsbach_memory(TRIM(model%name), '', '', model%land)
    model%land%mem%fract    (:,:) = 1._wp
    model%land%mem%fract_max(:,:) = 1._wp

!!$    ! Create glacier tile
!!$    tile => new_tile(model%land, 'glacier', 'Glacier tile')
!!$
!!$    ! Create bare soil tile
!!$    components(1) = 'srf' ; components(2) = 'soil' ; components(3) = ''
!!$    tile => new_tile(model%land, 'bare', 'Bare soil tile', components)
!!$    tile%mem => new_jsbach_memory(TRIM(model%name), '', '_bare', tile)
!!$    tile%mem%fract     = 0._wp  ! Temporary (must later be initialized from input file)
!!$    tile%mem%fract_max = 0._wp
!!$
!!$    ! Create vegetated tile
!!$    components(1) = 'srf' ; components(2) = 'soil' ; components(3) = 'veg'
!!$    tile => new_tile(model%land, 'veg', 'Vegetated tile', components)
!!$    tile%mem => new_jsbach_memory(TRIM(model%name), '', '_veg', tile)
!!$    ! Note: %fract will be the foliar projective cover, i.e. the fraction of the parent tile that is actually covered by a
!!$    !       canopy (veg_ratio), depending on the leaf area index.
!!$    !       %fract_max is the fraction of the parent tile that can potentially be covered by a canopy (veg_ratio_max)
!!$    tile%mem%fract     = 1._wp
!!$    tile%mem%fract_max = 1._wp  ! Temporary (must later be initialized from input file)

  END SUBROUTINE init_usecase_lite

  SUBROUTINE init_usecase_lite_tpe_open(model)

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':init_usecase_lite_tpe_open'

#ifdef HAVE_F2003
    TYPE(t_jsb_model), POINTER, INTENT(inout) :: model
#else
    TYPE(t_jsb_model), POINTER                :: model
#endif

    CALL init_usecase_lite(model)

  END SUBROUTINE init_usecase_lite_tpe_open

  SUBROUTINE init_usecase_lite_tpe_closed(model)

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':init_usecase_lite_tpe_closed'

#ifdef HAVE_F2003
    TYPE(t_jsb_model), POINTER, INTENT(inout) :: model
#else
    TYPE(t_jsb_model), POINTER                :: model
#endif

    CALL init_usecase_lite(model)

  END SUBROUTINE init_usecase_lite_tpe_closed

#endif
END MODULE mo_jsb_model_usecases
