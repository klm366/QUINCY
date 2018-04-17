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
  USE mo_jsb_control,         ONLY: debug_on

  USE mo_jsb_model_class,     ONLY: t_jsb_model  
  USE mo_jsb_class,           ONLY: get_model
  USE mo_jsb_tile_class,      ONLY: t_jsb_tile_abstract
  USE mo_jsb_io_netcdf,       ONLY: t_input_file, jsb_netcdf_open_input

#ifdef HAVE_MOFFLEDAC
  USE mo_timer,             ONLY: timer_start, timer_stop, timer_input
  USE mo_input,             ONLY: InputVarGroupStart, InputVarGroupEnd, InputFileAdd, InputVarAdd, InputOptAdd, InputDimAdd, &
                                  InputGetData, InputVarReleaseMem, input_file_list, input_var_list
#endif

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
    dsl4jsb_Def_memory(PHENO_)
    dsl4jsb_Get_memory(PHENO_)
    !dsl4jsb_Def_config(PHENO_)
    !dsl4jsb_Get_config(PHENO_)

    ! Declare pointers for variables in memory

    ! ---------------------------
    ! Go

    model => get_model(tile%owner_model_id)

    CALL pheno_init_bc(tile)
    CALL pheno_init_ic(tile)

  END SUBROUTINE pheno_init

  SUBROUTINE pheno_init_bc(tile)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile

#ifdef HAVE_MOFFLEDAC
    TYPE (input_file_list), POINTER :: IFile
    TYPE (input_var_list),  POINTER :: var, var2
#endif

    TYPE(t_jsb_model),       POINTER :: model

    dsl4jsb_Def_config(PHENO_)
    dsl4jsb_Def_memory(PHENO_)

    REAL(wp), POINTER ::  &
      & ptr_mon(:,:,:),   &
#ifdef HAVE_MOFFLEDAC
      & ptr_mon2(:,:,:),  &
#endif
      & ptr_2D(:,:)

    TYPE(t_input_file) :: input_file

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':pheno_init_bc'

    model => get_model(tile%owner_model_id)

    dsl4jsb_Get_config(PHENO_)
    dsl4jsb_Get_memory(PHENO_)

    IF (debug_on()) CALL message(TRIM(routine), 'Setting boundary conditions of pheno memory for tile '// &
      &                          TRIM(tile%name)//' from '//TRIM(dsl4jsb_Config(PHENO_)%bc_filename))

#ifdef HAVE_MOFFLEDAC
    CALL InputDimAdd('clim_monthly_year',Extent=12)
    CALL InputFileAdd(TRIM(dsl4jsb_Config(PHENO_)%bc_filename),file_ref=IFile,dt_file=0)
    CALL InputVarGroupStart('pheno_'//TRIM(tile%name),IFile=IFile,dt_update=0)
    CALL InputVarAdd('nproma,nblock,clim_monthly_year',ptr_mon, 'lai_clim_' //TRIM(tile%name),var_ref=var)
    CALL InputVarAdd('nproma,nblock,clim_monthly_year',ptr_mon2,'veg_fract_'//TRIM(tile%name),var_ref=var2)
    CALL InputVarGroupEnd()
    CALL InputOptAdd(AFile=IFile,dim_name='time=>clim_monthly_year')
    CALL InputOptAdd(AVar=var, var_file_name='lai_clim', fill_val=0._wp,valid_range=(/0._wp,1.e300_wp/))
    CALL InputOptAdd(AVar=var2,var_file_name='veg_fract',fill_val=0._wp,valid_range=(/0._wp,1.e300_wp/))
    CALL InputOptAdd(var_name='pheno_fract_forest_'//TRIM(tile%name),var_file_name='forest_fract', &
                     file_name=TRIM(dsl4jsb_Config(PHENO_)%bc_filename),fill_val=0._wp,valid_range=(/0._wp,1.e300_wp/))
    CALL InputOptAdd(var_name='pheno_fract_fpc_max_'//TRIM(tile%name),var_file_name='veg_ratio_max', &
                     file_name=TRIM(model%config%fract_filename),fill_val=0._wp,valid_range=(/0._wp,1.e300_wp/))
    CALL InputGetData(msg=TRIM(routine))
    dsl4jsb_memory(PHENO_)%lai_mon      (:,:,1:12) = ptr_mon (:,:,1:12)
    dsl4jsb_memory(PHENO_)%lai_mon      (:,:,0   ) = ptr_mon (:,:,12  )
    dsl4jsb_memory(PHENO_)%lai_mon      (:,:,13  ) = ptr_mon (:,:,1   )
    dsl4jsb_memory(PHENO_)%fract_fpc_mon(:,:,1:12) = ptr_mon2(:,:,1:12)
    dsl4jsb_memory(PHENO_)%fract_fpc_mon(:,:,0   ) = ptr_mon2(:,:,12  )
    dsl4jsb_memory(PHENO_)%fract_fpc_mon(:,:,13  ) = ptr_mon2(:,:,1   )
    CALL InputVarReleaseMem(var_ref=var)
    CALL InputVarReleaseMem(var_ref=var2)
#ELSE
    input_file = jsb_netcdf_open_input(TRIM(dsl4jsb_Config(PHENO_)%bc_filename), model%grid_id)

    ! Leaf area index climatology
    IF (debug_on()) CALL message(TRIM(routine), 'reading lai_clim ...')
    ptr_mon => input_file%Read_2d_time(    &
      & variable_name='lai_clim')
    ptr_mon(:,:,:) = MERGE(ptr_mon(:,:,:), 0._wp, ptr_mon(:,:,:) > 0._wp)
    dsl4jsb_memory(PHENO_)%lai_mon(:,:,1:12) = ptr_mon(:,:,1:12)
    dsl4jsb_memory(PHENO_)%lai_mon(:,:,0   ) = ptr_mon(:,:,12  )
    dsl4jsb_memory(PHENO_)%lai_mon(:,:,13  ) = ptr_mon(:,:,1   )
    DEALLOCATE(ptr_mon)

    ! Foliage projected cover climatology
    IF (debug_on()) CALL message(TRIM(routine), 'reading veg_fract ...')
    ptr_mon => input_file%Read_2d_time(    &
      & variable_name='veg_fract')
    ptr_mon(:,:,:) = MERGE(ptr_mon(:,:,:), 0._wp, ptr_mon(:,:,:) > 0._wp)
    dsl4jsb_memory(PHENO_)%fract_fpc_mon(:,:,1:12) = ptr_mon(:,:,1:12)
    dsl4jsb_memory(PHENO_)%fract_fpc_mon(:,:,0   ) = ptr_mon(:,:,12  )
    dsl4jsb_memory(PHENO_)%fract_fpc_mon(:,:,13  ) = ptr_mon(:,:,1   )
    DEALLOCATE(ptr_mon)

    ! Forest fraction
    ! @todo: only used for use_alb_veg_simple=.TRUE.
    IF (debug_on()) CALL message(TRIM(routine), 'reading forest_fract ...')
    ptr_2D => input_file%Read_2d(        &
      & variable_name='forest_fract',            &
      & fill_array = dsl4jsb_var_ptr(PHENO_, fract_forest))
    dsl4jsb_var2D_onDomain(PHENO_, fract_forest) = &
      MERGE(dsl4jsb_var2D_onDomain(PHENO_, fract_forest), 0._wp, dsl4jsb_var2D_onDomain(PHENO_, fract_forest) >= 0._wp)

    CALL input_file%Close()

    input_file = jsb_netcdf_open_input(TRIM(model%config%fract_filename), model%grid_id)     ! Temporary

    ! Maximum foliage projected cover
    IF (debug_on()) CALL message(TRIM(routine), 'reading veg_ratio_max ...')
    ptr_2D => input_file%Read_2d(        &
      & variable_name='veg_ratio_max',           &
      & fill_array = dsl4jsb_var_ptr(PHENO_, fract_fpc_max))

    dsl4jsb_var2D_onDomain(PHENO_, fract_fpc_max) = &
      MERGE(dsl4jsb_var2D_onDomain(PHENO_, fract_fpc_max), 0._wp, dsl4jsb_var2D_onDomain(PHENO_, fract_fpc_max) >= 0._wp)

    NULLIFY(ptr_2D)

    CALL input_file%Close()

#endif

    ! R: veg_fract_correction is used in many modules in JS3 (e.g. in the albedo calculations). However, in JS4
    !    it is only necessary for the assimilation and for the carbon distribution of the tile (and in the
    !    disturbance process to calcultate fuel). As it is calculated for each tile from two corresponding values
    !    of the lctlib, it is not necesseary to calculate it for each time step again as in JSBACH3.
    ! Calculate veg_fract_correction
    IF ( tile%lcts(1)%lib_id /= 0 ) THEN ! only if the present tile is a pft!
      dsl4jsb_var2D_onDomain(PHENO_, veg_fract_correction) =   &
        & 1.0_wp - exp(-dsl4jsb_Lctlib_param(MaxLAI) / dsl4jsb_Lctlib_param(clumpinessFactor))
    END IF
    
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

    IF (debug_on()) CALL message(TRIM(routine), 'Setting initial conditions of pheno memory (lai) for tile '// &
      &                          TRIM(tile%name)//' from '//TRIM(dsl4jsb_Config(PHENO_)%bc_filename))

    ! Set initial leaf area index from climatology (read in pheno_init_bc)
    CALL get_time_interpolation_weights(wgt1, wgt2, nmw1, nmw2)
    dsl4jsb_var2D_onDomain(PHENO_, lai) = &
      & wgt1 * dsl4jsb_memory(PHENO_)%lai_mon(:,:,nmw1) + wgt2 * dsl4jsb_memory(PHENO_)%lai_mon(:,:,nmw2)

    IF (model%config%l_compat401) THEN
      dsl4jsb_var2D_onDomain(PHENO_, fract_fpc) = get_foliage_projected_cover( &
        & dsl4jsb_var2D_onDomain(PHENO_, fract_fpc_max), dsl4jsb_var2D_onDomain(PHENO_, lai), 2._wp)
    ELSE
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
    END IF

  END SUBROUTINE pheno_init_ic
  
#endif
END MODULE mo_pheno_init
