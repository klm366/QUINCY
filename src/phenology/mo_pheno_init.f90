!>
!! @brief Initialization of the the phenology memory.
!!
!! @author
!!  Rainer Schneck, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Rainer Schneck (2015-09-29)
!!



MODULE mo_pheno_init
#ifndef __NO_JSBACH__

  USE mo_kind,                ONLY: wp
  USE mo_exception,           ONLY: message

  USE mo_jsb_model_class,     ONLY: t_jsb_model  
  USE mo_jsb_class,           ONLY: get_model
  USE mo_jsb_tile_class,      ONLY: t_jsb_tile_abstract
  USE mo_jsb_io_netcdf,       ONLY: t_input_file, jsb_netcdf_open_input


  dsl4jsb_Use_processes PHENO_
  dsl4jsb_Use_config(PHENO_)
  dsl4jsb_Use_memory(PHENO_)
    
  IMPLICIT NONE
  PRIVATE
  
  PUBLIC :: pheno_init

  CHARACTER(len=*), PARAMETER :: modname = 'mo_pheno_init'

CONTAINS

  SUBROUTINE pheno_init(tile)
    
    ! Arguments
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    
    ! Local variables
    TYPE(t_jsb_model), POINTER  :: model

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':pheno_init'

    ! Declare pointers for process configuration and memory

    ! Declare pointers for variables in memory

    ! ---------------------------
    ! Go
    
    model => get_model(tile%owner_model_id)

    CALL pheno_init_bc(tile)
    CALL pheno_init_ic(tile)

  END SUBROUTINE pheno_init

  SUBROUTINE pheno_init_bc(tile)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile

    TYPE(t_jsb_model),       POINTER :: model

    dsl4jsb_Def_config(PHENO_)
    dsl4jsb_Def_memory(PHENO_)

    REAL(wp), POINTER ::  &
      & ptr_mon(:,:,:),   &
      & ptr_2D(:,:)

    TYPE(t_input_file) :: input_file

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':pheno_init_bc'

    model => get_model(tile%owner_model_id)

    dsl4jsb_Get_config(PHENO_)
    dsl4jsb_Get_memory(PHENO_)

    CALL message(TRIM(routine), 'Setting boundary conditions of pheno memory for tile '// &
      &                          TRIM(tile%name)//' from '//TRIM(dsl4jsb_Config(PHENO_)%bc_filename))

    input_file = jsb_netcdf_open_input(TRIM(dsl4jsb_Config(PHENO_)%bc_filename), model%grid_id)

    ! Leaf area index climatology
    ptr_mon => input_file%Read_2d_time(    &
      & variable_name='lai_clim')
    ptr_mon(:,:,:) = MERGE(ptr_mon(:,:,:), 0._wp, ptr_mon(:,:,:) > 0._wp)
    ! TODO Implement proper macros for getting lai_mon from mem
    dsl4jsb_memory(PHENO_)%lai_mon(:,:,1:12) = ptr_mon(:,:,1:12)
    dsl4jsb_memory(PHENO_)%lai_mon(:,:,0   ) = ptr_mon(:,:,12  )
    dsl4jsb_memory(PHENO_)%lai_mon(:,:,13  ) = ptr_mon(:,:,1   )
    DEALLOCATE(ptr_mon)
    ! Foliage projected cover climatology
    ptr_mon => input_file%Read_2d_time(    &
      & variable_name='veg_fract')
    ptr_mon(:,:,:) = MERGE(ptr_mon(:,:,:), 0._wp, ptr_mon(:,:,:) > 0._wp)
    ! TODO Implement proper macros for getting lai_mon from mem
    dsl4jsb_memory(PHENO_)%fract_fpc_mon(:,:,1:12) = ptr_mon(:,:,1:12)
    dsl4jsb_memory(PHENO_)%fract_fpc_mon(:,:,0   ) = ptr_mon(:,:,12  )
    dsl4jsb_memory(PHENO_)%fract_fpc_mon(:,:,13  ) = ptr_mon(:,:,1   )
    DEALLOCATE(ptr_mon)

    ! Forest fraction
    ! @todo: only used for use_alb_veg_simple=.TRUE.
    ptr_2D => input_file%Read_2d(        &
      & variable_name='forest_fract',            &
      & fill_array = dsl4jsb_var_ptr(PHENO_, fract_forest))
    dsl4jsb_var2D_onDomain(PHENO_, fract_forest) = &
      MERGE(dsl4jsb_var2D_onDomain(PHENO_, fract_forest), 0._wp, dsl4jsb_var2D_onDomain(PHENO_, fract_forest) >= 0._wp)

    CALL input_file%Close()

    input_file = jsb_netcdf_open_input('bc_land_frac.nc', model%grid_id)     ! Temporary

    ! Maximum foliage projected cover
    ptr_2D => input_file%Read_2d(        &
      & variable_name='veg_ratio_max',           &
      & fill_array = dsl4jsb_var_ptr(PHENO_, fract_fpc_max))
    dsl4jsb_var2D_onDomain(PHENO_, fract_fpc_max) = &
      MERGE(dsl4jsb_var2D_onDomain(PHENO_, fract_fpc_max), 0._wp, dsl4jsb_var2D_onDomain(PHENO_, fract_fpc_max) >= 0._wp)

    NULLIFY(ptr_2D)

    CALL input_file%Close()

  END SUBROUTINE pheno_init_bc

  SUBROUTINE pheno_init_ic(tile)

    USE mo_pheno_process, ONLY: get_foliage_projected_cover
    USE mo_jsb_time,      ONLY: get_time_interpolation_weights

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile

    TYPE(t_jsb_model),       POINTER :: model

    dsl4jsb_Def_config(PHENO_)
    dsl4jsb_Def_memory(PHENO_)

    REAL(wp) :: wgt1, wgt2
    INTEGER  :: nmw1,nmw2

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':pheno_init_ic'

    model => get_model(tile%owner_model_id)

    dsl4jsb_Get_config(PHENO_)
    dsl4jsb_Get_memory(PHENO_)

    CALL message(TRIM(routine), 'Setting initial conditions of pheno memory (lai) for tile '// &
      &                          TRIM(tile%name)//' from '//TRIM(dsl4jsb_Config(PHENO_)%bc_filename))

    ! Set initial leaf area index from climatology (read in pheno_init_bc)
    CALL get_time_interpolation_weights(wgt1, wgt2, nmw1, nmw2)
    dsl4jsb_var2D_onDomain(PHENO_, lai) = &
      & wgt1 * dsl4jsb_memory(PHENO_)%lai_mon(:,:,nmw1) + wgt2 * dsl4jsb_memory(PHENO_)%lai_mon(:,:,nmw2)

    IF (       TRIM(dsl4jsb_Config(PHENO_)%scheme) == 'climatology' &
        & .OR. tile%lcts(1)%lib_id == 0) THEN         ! lib_id=0: eneral vegetation tile, no specific PFT
      dsl4jsb_var2D_onDomain(PHENO_, fract_fpc) = &
        & wgt1 * dsl4jsb_memory(PHENO_)%fract_fpc_mon(:,:,nmw1) + wgt2 * dsl4jsb_memory(PHENO_)%fract_fpc_mon(:,:,nmw2)
    ELSE
      dsl4jsb_var2D_onDomain(PHENO_, lai) = MIN(dsl4jsb_var2D_onDomain(PHENO_, lai), dsl4jsb_Lctlib_param(MaxLai))
      dsl4jsb_var2D_onDomain(PHENO_, fract_fpc) = get_foliage_projected_cover( &
        & dsl4jsb_var2D_onDomain(PHENO_, fract_fpc_max), dsl4jsb_var2D_onDomain(PHENO_, lai), &
        & dsl4jsb_Lctlib_param(ClumpinessFactor))
    END IF

  END SUBROUTINE pheno_init_ic
  
#endif
END MODULE mo_pheno_init
