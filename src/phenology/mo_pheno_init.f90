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
  USE mo_exception,           ONLY: message, finish
  USE mo_jsb_control,         ONLY: debug_on

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

  !> Run phenology init
  !! 
  !! run jsb4 or/and quincy init routine
  !! 
  SUBROUTINE pheno_init(tile)

    ! Arguments
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile

    ! Local
    TYPE(t_jsb_model), POINTER :: model

    CHARACTER(len=*), PARAMETER :: routine = modname//':pheno_init'

    model => Get_model(tile%owner_model_id)

    ! select init routine depending on the "use_quincy" flag
    IF (.NOT. model%config%use_quincy) THEN
      CALL pheno_init_bc(tile)
      CALL pheno_init_ic(tile)
    ELSE IF (model%config%use_quincy) THEN     ! quincy
      CALL pheno_init_q_(tile)                 ! contains code parts of jsbach4 pheno init
    END IF

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

    CHARACTER(len=*), PARAMETER :: routine = modname//':pheno_init_bc'

    model => get_model(tile%owner_model_id)

    dsl4jsb_Get_config(PHENO_)
    dsl4jsb_Get_memory(PHENO_)

    IF (debug_on()) CALL message(TRIM(routine), 'Setting boundary conditions of pheno memory for tile '// &
      &                          TRIM(tile%name)//' from '//TRIM(dsl4jsb_Config(PHENO_)%bc_filename))

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
    !$acc enter data attach(dsl4jsb_memory(PHENO_)%lai_mon)
    !$acc update device(dsl4jsb_memory(PHENO_)%lai_mon)

    ! Foliage projected cover climatology
    IF (debug_on()) CALL message(TRIM(routine), 'reading veg_fract ...')
    ptr_mon => input_file%Read_2d_time(    &
      & variable_name='veg_fract')
    ptr_mon(:,:,:) = MERGE(ptr_mon(:,:,:), 0._wp, ptr_mon(:,:,:) > 0._wp)
    dsl4jsb_memory(PHENO_)%fract_fpc_mon(:,:,1:12) = ptr_mon(:,:,1:12)
    dsl4jsb_memory(PHENO_)%fract_fpc_mon(:,:,0   ) = ptr_mon(:,:,12  )
    dsl4jsb_memory(PHENO_)%fract_fpc_mon(:,:,13  ) = ptr_mon(:,:,1   )
    DEALLOCATE(ptr_mon)
    !$acc update device(dsl4jsb_memory(PHENO_)%fract_fpc_mon)

    ! Forest fraction
    ! @todo: only used for use_alb_veg_simple=.TRUE. and HYDRO_ init if l_organic=.TRUE.
    IF (debug_on()) CALL message(TRIM(routine), 'reading forest_fract ...')
    ptr_2D => input_file%Read_2d(      &
      & variable_name='forest_fract',  &
      & fill_array = dsl4jsb_var_ptr(PHENO_, fract_forest))
    dsl4jsb_var2D_onDomain(PHENO_, fract_forest) = &
      MERGE(dsl4jsb_var2D_onDomain(PHENO_, fract_forest), 0._wp, dsl4jsb_var2D_onDomain(PHENO_, fract_forest) >= 0._wp)
    !$acc update device(dsl4jsb_var2D_onDomain(PHENO_, fract_forest))

    CALL input_file%Close()

    input_file = jsb_netcdf_open_input(TRIM(model%config%fract_filename), model%grid_id)     ! Temporary

    ! Maximum foliage projected cover
    IF (debug_on()) CALL message(TRIM(routine), 'reading veg_ratio_max ...')
    ptr_2D => input_file%Read_2d(      &
      & variable_name='veg_ratio_max', &
      & fill_array = dsl4jsb_var_ptr(PHENO_, fract_fpc_max))

    dsl4jsb_var2D_onDomain(PHENO_, fract_fpc_max) = &
      MERGE(dsl4jsb_var2D_onDomain(PHENO_, fract_fpc_max), 0._wp, dsl4jsb_var2D_onDomain(PHENO_, fract_fpc_max) >= 0._wp)
    !$acc update device(dsl4jsb_var2D_onDomain(PHENO_, fract_fpc_max))

    CALL input_file%Close()

    ! avoid compiler warnings about pointers not being dereferenced
    IF (ASSOCIATED(ptr_2D)) CONTINUE
    NULLIFY(ptr_2D)

    IF (tile%lcts(1)%lib_id /= 0)  THEN ! only if the present tile is a pft
      !JN 10.01.18: for now: set the maxLai for non forest types to the one from the Lctlib
      !             later one might think about allometric relation ships in other multi annual plants, too?
      IF ( (dsl4jsb_Config(PHENO_)%l_forestRegrowth) .AND. (.NOT. dsl4jsb_Lctlib_param(ForestFlag)) ) THEN
        dsl4jsb_var2D_onDomain(PHENO_, maxLAI_allom) = dsl4jsb_Lctlib_param(MaxLAI)
        !$acc update device(dsl4jsb_var2D_onDomain(PHENO_, maxLAI_allom))
      END IF

      !JN 10.01.18: with l_forestRegrowth the maxLAI can change over time for forest pfts
      IF ( (dsl4jsb_Config(PHENO_)%l_forestRegrowth) .AND. (dsl4jsb_Lctlib_param(ForestFlag)) ) THEN
        dsl4jsb_var2D_onDomain(PHENO_, veg_fract_correction) =   &
          & 1.0_wp - exp(-dsl4jsb_var2D_onDomain(PHENO_, maxLAI_allom) / dsl4jsb_Lctlib_param(clumpinessFactor))
      !@todo: do not rely on the lib_id but use a flag instead?
      ELSE
        dsl4jsb_var2D_onDomain(PHENO_, veg_fract_correction) =   &
          & 1.0_wp - exp(-dsl4jsb_Lctlib_param(MaxLAI) / dsl4jsb_Lctlib_param(clumpinessFactor))
      END IF
  
     !-=-=-=-=--=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
      !PHS: JeDi does not need this vegetation correction so I set the correction fraction to one.
      !-=-=-=-=--=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=        
      IF(TRIM(dsl4jsb_Config(PHENO_)%scheme) == 'jedi_phenology') THEN
          dsl4jsb_var2D_onDomain(PHENO_, veg_fract_correction) =  1.0_wp 
      ELSE
          dsl4jsb_var2D_onDomain(PHENO_, veg_fract_correction) =   &
           & 1.0_wp - exp(-dsl4jsb_Lctlib_param(MaxLAI) / dsl4jsb_Lctlib_param(clumpinessFactor))
      END IF
      !-=-=-=-=--=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=      

      !$acc update device(dsl4jsb_var2D_onDomain(PHENO_, veg_fract_correction))

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

    CHARACTER(len=*), PARAMETER :: routine = modname//':pheno_init_ic'

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
      !$claw expand
      dsl4jsb_var2D_onDomain(PHENO_, fract_fpc) = get_foliage_projected_cover( &
        & dsl4jsb_var2D_onDomain(PHENO_, fract_fpc_max), dsl4jsb_var2D_onDomain(PHENO_, lai), 2._wp)
      !$acc update device(dsl4jsb_var2D_onDomain(PHENO_, fract_fpc))
    ELSE
      IF (       TRIM(dsl4jsb_Config(PHENO_)%scheme) == 'climatology' &
          & .OR. tile%lcts(1)%lib_id == 0) THEN         ! lib_id=0: eneral vegetation tile, no specific PFT
        dsl4jsb_var2D_onDomain(PHENO_, fract_fpc) = &
          & wgt1 * dsl4jsb_memory(PHENO_)%fract_fpc_mon(:,:,nmw1) + wgt2 * dsl4jsb_memory(PHENO_)%fract_fpc_mon(:,:,nmw2)
      ELSE
        IF ( .NOT. dsl4jsb_Config(PHENO_)%l_forestRegrowth )  THEN
        
           !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
           ! BCBCBCBCBCBCBCBCBC  in order to set the initial lai to zero
           !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
           IF(TRIM(dsl4jsb_Config(PHENO_)%scheme) == 'jedi_phenology') THEN
              dsl4jsb_var2D_onDomain(PHENO_, lai) = 0._wp
           ELSE
           !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-    
              dsl4jsb_var2D_onDomain(PHENO_, lai) = MIN(dsl4jsb_var2D_onDomain(PHENO_, lai), dsl4jsb_Lctlib_param(MaxLai))
           ENDIF
        
        
        ELSE
          dsl4jsb_var2D_onDomain(PHENO_, lai) = &
            & MIN(dsl4jsb_var2D_onDomain(PHENO_, lai), dsl4jsb_var2D_onDomain(PHENO_, maxLAI_allom))
        END IF
        !$claw expand
        dsl4jsb_var2D_onDomain(PHENO_, fract_fpc) = get_foliage_projected_cover(                &
          & dsl4jsb_var2D_onDomain(PHENO_, fract_fpc_max), dsl4jsb_var2D_onDomain(PHENO_, lai), &
          & dsl4jsb_Lctlib_param(ClumpinessFactor))
      END IF
      !$acc update device(dsl4jsb_var2D_onDomain(PHENO_, fract_fpc))
    END IF

    !$acc update device(dsl4jsb_var2D_onDomain(PHENO_, lai))

  END SUBROUTINE pheno_init_ic

  !-----------------------------------------------------------------------------------------------------
  !> Intialize pheno variables for boundary conditions (bc)
  !! 
  !! quincy
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE pheno_init_q_(tile)

    ! see also the USE statements in the module header
    USE mo_pheno_constants,       ONLY: ievergreen
    USE mo_jsb_impl_constants,    ONLY: false, true
    USE mo_pheno_process,         ONLY: get_foliage_projected_cover
    
    ! ---------------------------
    ! 0.1 InOut
    CLASS(t_jsb_tile_abstract), INTENT(inout)     :: tile         !< one tile with data structure for one lct
    ! ---------------------------
    ! 0.2 Local
    TYPE(t_jsb_model),       POINTER  :: model
    TYPE(t_input_file)                :: input_file
    REAL(wp),                POINTER  :: lai_clim(:,:,:)
    REAL(wp), POINTER                 :: ptr_2D(:,:) 
       
    CHARACTER(len=*), PARAMETER :: routine = modname//':pheno_init_q_'
    ! ---------------------------
    ! 0.3 Declare Memory
    ! Declare process configuration and memory Pointers
    dsl4jsb_Def_config(PHENO_)
    dsl4jsb_Def_memory(PHENO_)
    ! Declare pointers to variables in memory
    dsl4jsb_Real2D_onDomain    :: growing_season
    dsl4jsb_Real2D_onDomain    :: lai_max
    dsl4jsb_Real2D_onDomain    :: lai
    ! ---------------------------
    ! 0.4 Debug Option
    IF (debug_on()) CALL message(routine, 'Setting boundary conditions of pheno memory (quincy) for tile '// &
      &                          TRIM(tile%name))
    ! ---------------------------
    ! 0.5 Get Memory
    model => get_model(tile%owner_model_id)
    ! Get process config
    dsl4jsb_Get_config(PHENO_)
    ! Get process memories
    dsl4jsb_Get_memory(PHENO_)
    ! Get process variables (Set pointers to variables in memory)
    dsl4jsb_Get_var2D_onDomain(PHENO_, growing_season)
    dsl4jsb_Get_var2D_onDomain(PHENO_, lai_max)
    dsl4jsb_Get_var2D_onDomain(PHENO_, lai)

    ! ------------------------------------------------------------------------------------------------------------
    ! Go science
    ! ------------------------------------------------------------------------------------------------------------

    !> 1.0 original quincy code
    !!
    !! init growing_season, lai_max, lai
    IF (tile%lcts(1)%lib_id /= 0) THEN ! work with lctlib only if the present tile is a pft
      ! lai_max == jsbach4 lctlib parameter
      lai_max(:,:) = dsl4jsb_Lctlib_param(MaxLAI)
      ! lai init should also work with 0.0_wp (in QUINCY canopy it is initialized with 6.0_wp)
      lai(:,:)     = 1.0_wp
      ! growing season flag, set depending on whether this is a evergreen plant (TRUE) or not (FALSE)
      IF(dsl4jsb_Lctlib_param(phenology_type) == ievergreen) THEN
        growing_season(:,:) = true
      ELSE
        growing_season(:,:) = false
      ENDIF
    ELSE
      ! default for non-PFT tiles
      lai_max(:,:)        = 1.0_wp
      lai(:,:)            = 0.0_wp
      growing_season(:,:) = false
    ENDIF


    !> 2.0 jsbach4 "pheno init" code copied with minor modifications
    !!
    
    !> 2.1 get fract_fpc_max from bc file
    input_file = jsb_netcdf_open_input(TRIM(model%config%fract_filename), model%grid_id)     ! Temporary
    ! Maximum foliage projected cover
    IF (debug_on()) CALL message(TRIM(routine), 'reading veg_ratio_max ...')
    ptr_2D => input_file%Read_2d(      &
      & variable_name='veg_ratio_max', &
      & fill_array = dsl4jsb_var_ptr(PHENO_, fract_fpc_max))
    dsl4jsb_var2D_onDomain(PHENO_, fract_fpc_max) = &
      MERGE(dsl4jsb_var2D_onDomain(PHENO_, fract_fpc_max), 0._wp, dsl4jsb_var2D_onDomain(PHENO_, fract_fpc_max) >= 0._wp)
    !$acc update device(dsl4jsb_var2D_onDomain(PHENO_, fract_fpc_max))
    CALL input_file%Close()
    
    !> 2.2 calc initial fract_fpc 
    IF (tile%lcts(1)%lib_id /= 0) THEN ! work with lctlib only if the present tile is a pft
      !$claw expand
      dsl4jsb_var2D_onDomain(PHENO_, fract_fpc) = get_foliage_projected_cover(                &
        & dsl4jsb_var2D_onDomain(PHENO_, fract_fpc_max), dsl4jsb_var2D_onDomain(PHENO_, lai), &
        & dsl4jsb_Lctlib_param(ClumpinessFactor))
    ELSE
      ! other than PFT tiles (in jsbach4 used in case l_compat401=.TRUE.)
      !$claw expand
      dsl4jsb_var2D_onDomain(PHENO_, fract_fpc) = get_foliage_projected_cover( &
        & dsl4jsb_var2D_onDomain(PHENO_, fract_fpc_max), dsl4jsb_var2D_onDomain(PHENO_, lai), 2._wp)
    ENDIF
    
    !> 2.3 veg_fract_correction (a factor between 0 and 1, not modified at runtime, see pheno mem class for some docu)
    IF (tile%lcts(1)%lib_id /= 0) THEN ! work with lctlib only if the present tile is a pft
      dsl4jsb_var2D_onDomain(PHENO_, veg_fract_correction) =   &
        & 1.0_wp - exp(-dsl4jsb_Lctlib_param(MaxLAI) / dsl4jsb_Lctlib_param(clumpinessFactor))
    ENDIF

  END SUBROUTINE pheno_init_q_

#endif
END MODULE mo_pheno_init
