!
! This file was automatically generated from JSBACH source file:
!    src/phenology/mo_pheno_init.f90
! with md5 hash: e74cb3756f078fc832e0353dff5b0244
!
!>
!! @brief Initialization of the the phenology memory.
!!
!! @author
!!  Rainer Schneck, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Rainer Schneck (2015-09-29)
!! Julia Nabel (2018-01-10): l_forestRegrowth flag
!!
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/license .
!!

MODULE mo_pheno_init
#ifndef __NO_JSBACH__

  USE mo_kind,                ONLY: wp
  USE mo_exception,           ONLY: message
  USE mo_jsb_control,         ONLY: debug_on

  USE mo_jsb_model_class,     ONLY: t_jsb_model
  USE mo_jsb_class,           ONLY: get_model
  USE mo_jsb_tile_class,      ONLY: t_jsb_tile_abstract
  USE mo_jsb_io_netcdf,       ONLY: t_input_file, jsb_netcdf_open_input

  USE mo_jsb_process_class, ONLY: pheno_
  USE mo_pheno_config_class, ONLY: t_pheno_config
  USE mo_pheno_memory_class, ONLY: t_pheno_memory

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: pheno_init

  CHARACTER(len=*), PARAMETER :: modname = 'mo_pheno_init'

CONTAINS

  SUBROUTINE pheno_init(tile)

    ! Arguments
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile

    CHARACTER(len=*), PARAMETER :: routine = modname//':pheno_init'

    CALL pheno_init_bc(tile)
    CALL pheno_init_ic(tile)

  END SUBROUTINE pheno_init

  SUBROUTINE pheno_init_bc(tile)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile

    TYPE(t_jsb_model),       POINTER :: model

    TYPE(t_pheno_config), POINTER :: pheno__conf
    TYPE(t_pheno_memory), POINTER :: pheno__mem

    REAL(wp), POINTER ::  &
      & ptr_mon(:,:,:),   &
      & ptr_2D(:,:)

    TYPE(t_input_file) :: input_file

    CHARACTER(len=*), PARAMETER :: routine = modname//':pheno_init_bc'

    model => get_model(tile%owner_model_id)

    SELECT TYPE (cc => model%processes(pheno_)%p%config); TYPE IS (t_pheno_config) ; pheno__conf => cc ; END SELECT
    SELECT TYPE (mm => tile%mem(pheno_)%p); TYPE IS (t_pheno_memory); pheno__mem => mm ; END SELECT

    IF (debug_on()) CALL message(TRIM(routine), 'Setting boundary conditions of pheno memory for tile '// &
      &                          TRIM(tile%name)//' from '//TRIM(pheno__conf%bc_filename))

    input_file = jsb_netcdf_open_input(TRIM(pheno__conf%bc_filename), model%grid_id)

    ! Leaf area index climatology
    IF (debug_on()) CALL message(TRIM(routine), 'reading lai_clim ...')
    ptr_mon => input_file%Read_2d_time(    &
      & variable_name='lai_clim')
    ptr_mon(:,:,:) = MERGE(ptr_mon(:,:,:), 0._wp, ptr_mon(:,:,:) > 0._wp)
    pheno__mem%lai_mon(:,:,1:12) = ptr_mon(:,:,1:12)
    pheno__mem%lai_mon(:,:,0   ) = ptr_mon(:,:,12  )
    pheno__mem%lai_mon(:,:,13  ) = ptr_mon(:,:,1   )
    DEALLOCATE(ptr_mon)
    !$acc enter data attach(pheno__mem%lai_mon)
    !$acc update device(pheno__mem%lai_mon)

    ! Foliage projected cover climatology
    IF (debug_on()) CALL message(TRIM(routine), 'reading veg_fract ...')
    ptr_mon => input_file%Read_2d_time(    &
      & variable_name='veg_fract')
    ptr_mon(:,:,:) = MERGE(ptr_mon(:,:,:), 0._wp, ptr_mon(:,:,:) > 0._wp)
    pheno__mem%fract_fpc_mon(:,:,1:12) = ptr_mon(:,:,1:12)
    pheno__mem%fract_fpc_mon(:,:,0   ) = ptr_mon(:,:,12  )
    pheno__mem%fract_fpc_mon(:,:,13  ) = ptr_mon(:,:,1   )
    DEALLOCATE(ptr_mon)
    !$acc update device(pheno__mem%fract_fpc_mon)

    ! Forest fraction
    ! @todo: only used for use_alb_veg_simple=.TRUE.
    IF (debug_on()) CALL message(TRIM(routine), 'reading forest_fract ...')
    ptr_2D => input_file%Read_2d(      &
      & variable_name='forest_fract',  &
      & fill_array = pheno__mem%fract_forest%ptr)
    pheno__mem%fract_forest%ptr(:,:) = &
      MERGE(pheno__mem%fract_forest%ptr(:,:), 0._wp, pheno__mem%fract_forest%ptr(:,:) >= 0._wp)
    !$acc update device(pheno__mem%fract_forest%ptr(:,:))

    CALL input_file%Close()

    input_file = jsb_netcdf_open_input(TRIM(model%config%fract_filename), model%grid_id)     ! Temporary

    ! Maximum foliage projected cover
    IF (debug_on()) CALL message(TRIM(routine), 'reading veg_ratio_max ...')
    ptr_2D => input_file%Read_2d(      &
      & variable_name='veg_ratio_max', &
      & fill_array = pheno__mem%fract_fpc_max%ptr)

    pheno__mem%fract_fpc_max%ptr(:,:) = &
      MERGE(pheno__mem%fract_fpc_max%ptr(:,:), 0._wp, pheno__mem%fract_fpc_max%ptr(:,:) >= 0._wp)
    !$acc update device(pheno__mem%fract_fpc_max%ptr(:,:))

    CALL input_file%Close()

    ! avoid compiler warnings about pointers not being dereferenced
    IF (ASSOCIATED(ptr_2D)) CONTINUE
    NULLIFY(ptr_2D)

    IF (tile%lcts(1)%lib_id /= 0)  THEN ! only if the present tile is a pft
      !JN 10.01.18: for now: set the maxLai for non forest types to the one from the Lctlib
      !             later one might think about allometric relation ships in other multi annual plants, too?
      IF ( (pheno__conf%l_forestRegrowth) .AND. (.NOT. model%lctlib(tile%lcts(1)%lib_id)%ForestFlag) ) THEN
        pheno__mem%maxlai_allom%ptr(:,:) = model%lctlib(tile%lcts(1)%lib_id)%MaxLAI
        !$acc update device(pheno__mem%maxlai_allom%ptr(:,:))
      END IF

      !JN 10.01.18: with l_forestRegrowth the maxLAI can change over time for forest pfts
      IF ( (pheno__conf%l_forestRegrowth) .AND. (model%lctlib(tile%lcts(1)%lib_id)%ForestFlag) ) THEN
        pheno__mem%veg_fract_correction%ptr(:,:) =   &
          & 1.0_wp - exp(-pheno__mem%maxlai_allom%ptr(:,:) / model%lctlib(tile%lcts(1)%lib_id)%clumpinessFactor)
      !@todo: do not rely on the lib_id but use a flag instead?
      ELSE
        pheno__mem%veg_fract_correction%ptr(:,:) =   &
          & 1.0_wp - exp(-model%lctlib(tile%lcts(1)%lib_id)%MaxLAI / model%lctlib(tile%lcts(1)%lib_id)%clumpinessFactor)
      END IF
      !$acc update device(pheno__mem%veg_fract_correction%ptr(:,:))
    END IF

  END SUBROUTINE pheno_init_bc

  SUBROUTINE pheno_init_ic(tile)

    USE mo_pheno_process, ONLY: get_foliage_projected_cover
    USE mo_jsb_time,      ONLY: get_time_interpolation_weights

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile

    TYPE(t_jsb_model),       POINTER :: model

    TYPE(t_pheno_config), POINTER :: pheno__conf
    TYPE(t_pheno_memory), POINTER :: pheno__mem

    REAL(wp) :: wgt1, wgt2
    INTEGER  :: nmw1,nmw2

    CHARACTER(len=*), PARAMETER :: routine = modname//':pheno_init_ic'

    model => get_model(tile%owner_model_id)

    SELECT TYPE (cc => model%processes(pheno_)%p%config); TYPE IS (t_pheno_config) ; pheno__conf => cc ; END SELECT
    SELECT TYPE (mm => tile%mem(pheno_)%p); TYPE IS (t_pheno_memory); pheno__mem => mm ; END SELECT

    IF (debug_on()) CALL message(TRIM(routine), 'Setting initial conditions of pheno memory (lai) for tile '// &
      &                          TRIM(tile%name)//' from '//TRIM(pheno__conf%bc_filename))

    ! Set initial leaf area index from climatology (read in pheno_init_bc)
    CALL get_time_interpolation_weights(wgt1, wgt2, nmw1, nmw2)
    pheno__mem%lai%ptr(:,:) = &
      & wgt1 * pheno__mem%lai_mon(:,:,nmw1) + wgt2 * pheno__mem%lai_mon(:,:,nmw2)

    IF (model%config%l_compat401) THEN
      !$claw expand
      pheno__mem%fract_fpc%ptr(:,:) = get_foliage_projected_cover( &
        & pheno__mem%fract_fpc_max%ptr(:,:), pheno__mem%lai%ptr(:,:), 2._wp)
      !$acc update device(pheno__mem%fract_fpc%ptr(:,:))
    ELSE
      IF (       TRIM(pheno__conf%scheme) == 'climatology' &
          & .OR. tile%lcts(1)%lib_id == 0) THEN         ! lib_id=0: eneral vegetation tile, no specific PFT
        pheno__mem%fract_fpc%ptr(:,:) = &
          & wgt1 * pheno__mem%fract_fpc_mon(:,:,nmw1) + wgt2 * pheno__mem%fract_fpc_mon(:,:,nmw2)
      ELSE
        IF ( .NOT. pheno__conf%l_forestRegrowth )  THEN
          pheno__mem%lai%ptr(:,:) = MIN(pheno__mem%lai%ptr(:,:), model%lctlib(tile%lcts(1)%lib_id)%MaxLai)
        ELSE
          pheno__mem%lai%ptr(:,:) = &
            & MIN(pheno__mem%lai%ptr(:,:), pheno__mem%maxlai_allom%ptr(:,:))
        END IF
        !$claw expand
        pheno__mem%fract_fpc%ptr(:,:) = get_foliage_projected_cover(                &
          & pheno__mem%fract_fpc_max%ptr(:,:), pheno__mem%lai%ptr(:,:), &
          & model%lctlib(tile%lcts(1)%lib_id)%ClumpinessFactor)
      END IF
      !$acc update device(pheno__mem%fract_fpc%ptr(:,:))
    END IF

    !$acc update device(pheno__mem%lai%ptr(:,:))

  END SUBROUTINE pheno_init_ic

#endif
END MODULE mo_pheno_init
