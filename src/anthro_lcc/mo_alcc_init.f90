!>
!! @brief Initialization of the the alcc memory.
!!
!! @author
!!  Rainer Schneck, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Rainer Schneck (2018-11-28)
!! 
!!

MODULE mo_alcc_init
#ifndef __NO_JSBACH__

  USE mo_kind,                ONLY: wp
  USE mo_exception,           ONLY: message
  USE mo_jsb_control,         ONLY: debug_on

  USE mo_jsb_model_class,     ONLY: t_jsb_model
  USE mo_jsb_class,           ONLY: get_model
  USE mo_jsb_tile_class,      ONLY: t_jsb_tile_abstract
  USE mo_jsb_io_netcdf,       ONLY: t_input_file, jsb_netcdf_open_input

  dsl4jsb_Use_processes ALCC_
  dsl4jsb_Use_config(ALCC_)
  dsl4jsb_Use_memory(ALCC_)

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: alcc_init

  CHARACTER(len=*), PARAMETER :: modname = 'mo_alcc_init'

CONTAINS

  SUBROUTINE alcc_init(tile)

    ! Arguments
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile

    CHARACTER(len=*), PARAMETER :: routine = modname//':alcc_init'

    !CALL alcc_init_bc(tile)
    !CALL alcc_init_ic(tile)

  END SUBROUTINE alcc_init

!   SUBROUTINE alcc_init_bc(tile)
! 
!     CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
! 
!     TYPE(t_jsb_model),       POINTER :: model
! 
!     dsl4jsb_Def_config(ALCC_)
!     dsl4jsb_Def_memory(ALCC_)
! 
!     REAL(wp), POINTER ::  &
!       & ptr_mon(:,:,:),   &
!       & ptr_2D(:,:)
! 
!     TYPE(t_input_file) :: input_file
! 
!     CHARACTER(len=*), PARAMETER :: routine = modname//':alcc_init_bc'
! 
!     model => get_model(tile%owner_model_id)
! 
!     dsl4jsb_Get_config(ALCC_)
!     dsl4jsb_Get_memory(ALCC_)
! 
!     IF (debug_on()) CALL message(TRIM(routine), 'Setting boundary conditions of alcc memory for tile '// &
!       &                          TRIM(tile%name)//' from '//TRIM(dsl4jsb_Config(ALCC_)%bc_filename))
! 
! 
! 
!     ! Forest fraction
!     ! @todo: only used for use_alb_veg_simple=.TRUE.
!     IF (debug_on()) CALL message(TRIM(routine), 'reading forest_fract ...')
!     ptr_2D => input_file%Read_2d(      &
!       & variable_name='forest_fract',  &
!       & fill_array = dsl4jsb_var_ptr(ALCC_, fract_forest))
!     dsl4jsb_var2D_onDomain(ALCC_, fract_forest) = &
!       MERGE(dsl4jsb_var2D_onDomain(ALCC_, fract_forest), 0._wp, dsl4jsb_var2D_onDomain(ALCC_, fract_forest) >= 0._wp)
! 
!     CALL input_file%Close()
! 
!     input_file = jsb_netcdf_open_input(TRIM(model%config%fract_filename), model%grid_id)     ! Temporary
! 
!     ! Maximum foliage projected cover
!     IF (debug_on()) CALL message(TRIM(routine), 'reading veg_ratio_max ...')
!     ptr_2D => input_file%Read_2d(      &
!       & variable_name='veg_ratio_max', &
!       & fill_array = dsl4jsb_var_ptr(ALCC_, fract_fpc_max))
! 
!     dsl4jsb_var2D_onDomain(ALCC_, fract_fpc_max) = &
!       MERGE(dsl4jsb_var2D_onDomain(ALCC_, fract_fpc_max), 0._wp, dsl4jsb_var2D_onDomain(ALCC_, fract_fpc_max) >= 0._wp)
! 
!     CALL input_file%Close()
! 
!     ! avoid compiler warnings about pointers not being dereferenced
!     IF (ASSOCIATED(ptr_2D)) CONTINUE
!     NULLIFY(ptr_2D)
! 
!     IF (tile%lcts(1)%lib_id /= 0)  THEN ! only if the present tile is a pft
!       !JN 10.01.18: for now: set the maxLai for non forest types to the one from the Lctlib
!       !             later one might think about allometric relation ships in other multi annual plants, too?
!       IF ( (dsl4jsb_Config(ALCC_)%l_forestRegrowth) .AND. (.NOT. dsl4jsb_Lctlib_param(ForestFlag)) ) THEN
!         dsl4jsb_var2D_onDomain(ALCC_, maxLAI_allom) = dsl4jsb_Lctlib_param(MaxLAI)
!       END IF
! 
!       !JN 10.01.18: with l_forestRegrowth the maxLAI can change over time for forest pfts
!       IF ( (dsl4jsb_Config(ALCC_)%l_forestRegrowth) .AND. (dsl4jsb_Lctlib_param(ForestFlag)) ) THEN
!         dsl4jsb_var2D_onDomain(ALCC_, veg_fract_correction) =   &
!           & 1.0_wp - exp(-dsl4jsb_var2D_onDomain(ALCC_, maxLAI_allom) / dsl4jsb_Lctlib_param(clumpinessFactor))
!       !@todo: do not rely on the lib_id but use a flag instead?
!       ELSE
!         dsl4jsb_var2D_onDomain(ALCC_, veg_fract_correction) =   &
!           & 1.0_wp - exp(-dsl4jsb_Lctlib_param(MaxLAI) / dsl4jsb_Lctlib_param(clumpinessFactor))
!       END IF
!     END IF
! 
!   END SUBROUTINE alcc_init_bc
! 
!
!   SUBROUTINE alcc_init_ic(tile)
! 
!     USE mo_alcc_process, ONLY: get_foliage_projected_cover
!     USE mo_jsb_time,       ONLY: get_time_interpolation_weights
! 
!     CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
! 
!     TYPE(t_jsb_model),       POINTER :: model
! 
!     dsl4jsb_Def_config(ALCC_)
!     dsl4jsb_Def_memory(ALCC_)
! 
!     CHARACTER(len=*), PARAMETER :: routine = modname//':alcc_init_ic'
! 
!     model => get_model(tile%owner_model_id)
! 
!     dsl4jsb_Get_config(ALCC_)
!     dsl4jsb_Get_memory(ALCC_)
! 
!     IF (debug_on()) CALL message(TRIM(routine), 'Setting initial conditions of alcc memory for tile '// &
!       &                          TRIM(tile%name)//' from '//TRIM(dsl4jsb_Config(ALCC_)%ic_filename))
! 
! 
!     
!     IF ( TRIM(dsl4jsb_Config(ALCC_)%scheme) == 'maps' ) THEN
!       ! IF ( tile%lcts(1)%lib_id == 0) THEN         ! lib_id=0: general vegetation tile, no specific PFT
!       ! END IF
!       
!       ! Provide dummy values to indicate if before-fractions are really received from fract-initialization or restart file
!       dsl4jsb_var2D_onDomain(ALCC_, fract_pfts_before_alccstep) =7._wp
!     END IF
! 
!   END SUBROUTINE alcc_init_ic

#endif
END MODULE mo_alcc_init
