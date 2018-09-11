!>
!! @brief Initialization of the forest-age memory.
!!
!! @author
!!  Julia Nabel, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Julia Nabel (2017-11-21)
!!

MODULE mo_fage_init
#ifndef __NO_JSBACH__

  USE mo_kind,                ONLY: wp
  USE mo_exception,           ONLY: message, finish, message_text
  USE mo_jsb_control,         ONLY: debug_on

  USE mo_jsb_model_class,     ONLY: t_jsb_model
  USE mo_jsb_grid_class,      ONLY: t_jsb_vgrid, t_jsb_grid
  USE mo_jsb_class,           ONLY: get_model
  USE mo_jsb_grid,            ONLY: Get_grid
  USE mo_jsb_tile_class,      ONLY: t_jsb_tile_abstract
  USE mo_jsb_io_netcdf,       ONLY: t_input_file, jsb_netcdf_open_input

  dsl4jsb_Use_processes FAGE_
  dsl4jsb_Use_config(FAGE_)
  dsl4jsb_Use_memory(FAGE_)
    
  IMPLICIT NONE
  PRIVATE
  
  PUBLIC :: fage_init

  CHARACTER(len=*), PARAMETER :: modname = 'mo_fage_init'

CONTAINS

  SUBROUTINE fage_init(tile)
    
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    
    CHARACTER(len=*), PARAMETER :: routine = modname//':fage_init'
    
    CALL fage_init_bc(tile)
    CALL fage_init_ic(tile)

  END SUBROUTINE fage_init

  SUBROUTINE fage_init_bc(tile)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile

    dsl4jsb_Def_config(FAGE_)
    dsl4jsb_Def_memory(FAGE_)

    TYPE(t_input_file)          :: input_file
    TYPE(t_jsb_vgrid), POINTER  :: forest_age

    TYPE(t_jsb_model), POINTER  :: model

    CHARACTER(len=*), PARAMETER :: routine = modname//':fage_init_bc'

    INTEGER :: i, max_age, age_interval, ncohorts

    IF (.NOT. tile%Is_process_active(FAGE_)) RETURN

    model => get_model(tile%owner_model_id)

    dsl4jsb_Get_config(FAGE_)
    dsl4jsb_Get_memory(FAGE_)

    max_age = dsl4jsb_Config(FAGE_)%max_age
    ncohorts = dsl4jsb_Config(FAGE_)%ncohorts

    ! set/initialise age cohort bounds 
    ! - important: start with one ---> DO NOT CHANGE!  
    ! -- first age cohort needs to always only contain the first year!
    dsl4jsb_var3D_onDomain(FAGE_, cohort_max_age) = 1._wp 

    !TODO-JN if there is no dynamic merging of age cohorts, the cohort_max_age vector could be
    !        much more simple and set in the confic_class, etc
    IF (dsl4jsb_Config(FAGE_)%age_cohort_scheme == 'static-equal' ) THEN
      ! Equal intervals distributed equally over all but the first age cohort
      age_interval = int((max_age-1) / (ncohorts-1))
      DO i=1,ncohorts-1
        dsl4jsb_var_ptr(FAGE_, cohort_max_age) (:,i+1,:) &
          & = dsl4jsb_var_ptr(FAGE_, cohort_max_age) (:,i,:) + age_interval
      END DO
    ELSEIF (dsl4jsb_Config(FAGE_)%age_cohort_scheme == 'static-shorterYoung' ) THEN
      ! Equal intervals distributed unequally such that younger age cohorts have shorter steps
      ! (1, 1+n, 1+n+2*n, 1+n+2*n+3*n, ...)
      age_interval = 0
      DO i=1,ncohorts-1
        age_interval = age_interval + i
      END DO
      age_interval = int(max_age / age_interval)
      DO i=1,ncohorts-1
        dsl4jsb_var_ptr(FAGE_, cohort_max_age) (:,i+1,:) &
          & = dsl4jsb_var_ptr(FAGE_, cohort_max_age) (:,i,:) + (i * age_interval)
      END DO
    ELSE
      CALL finish(TRIM(routine), 'Currently only "static-shorterYoung" and "static-equal" implemented as forest ageing schemes')
    ENDIF

!TODO-JN: for now it would be sufficient to get it from iblk=1
!IF(tile%lcts(1)%lib_id .eq. 1) THEN
!print *, dsl4jsb_var_ptr(FAGE_, cohort_max_age)
!ENDIF

    !Last age cohort is not really restricted to a max age, but set to max_age for compliance with the fract_per_age field
    dsl4jsb_var_ptr(FAGE_, cohort_max_age) (:,ncohorts,:) = max_age

    !Assertion: 
    IF(ANY(dsl4jsb_var_ptr(FAGE_, cohort_max_age)(:,1,:) .NE. 1._wp)) THEN
      CALL finish(TRIM(routine), 'Assertion: lowest age cohort always needs to have max_age of 1!')
    ENDIF
    IF(ANY(dsl4jsb_var_ptr(FAGE_, cohort_max_age)(:,:,:) .GT. max_age)) THEN
      CALL finish(TRIM(routine), 'Assertion: no age cohort should have a max age larger than the total max_age!')
    ENDIF

    !TODO-JN other bcs?? Also from files?
!    dsl4jsb_var3D_onDomain(FM_, rotation_cycle) = REAL(dsl4jsb_Config(FM_)%global_rotation_cycle, wp)
!    IF (dsl4jsb_Config(FM_)%scheme == 'rotation_local' ) THEN 
!      IF (debug_on()) CALL message(TRIM(routine), 'Setting boundary conditions of fm memory (rotation_cycle) for tile '// &
!        &                          TRIM(tile%name)//' from '//TRIM(dsl4jsb_Config(FM_)%bc_filename))
!
!      input_file = jsb_netcdf_open_input(TRIM(dsl4jsb_Config(FM_)%bc_filename), model%grid_id)
!      ptr_2D => input_file%Read_2d(        &
!        & variable_name='rotation_cycle',           &
!        & fill_array = dsl4jsb_var_ptr(FM_, rotation_cycle))
!    ENDIF

  END SUBROUTINE fage_init_bc


  SUBROUTINE fage_init_ic(tile)

    USE mo_fage_util,      ONLY: get_mean_age

    !-----------------------------------------------------------------------
    !  ARGUMENTS
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile

    !-----------------------------------------------------------------------
    !  LOCAL VARIABLES
    TYPE(t_jsb_model), POINTER :: model
    TYPE(t_jsb_grid),  POINTER :: grid

!    REAL(wp), POINTER :: ptr_2D(:,:)
!    REAL(wp), POINTER :: ptr_3D(:,  :)

    TYPE(t_input_file) :: input_file
    TYPE(t_jsb_vgrid), POINTER :: forest_age

    dsl4jsb_Def_config(FAGE_)
    dsl4jsb_Def_memory(FAGE_)

    CHARACTER(len=*), PARAMETER :: routine = modname//':fage_init_ic'

    !-----------------------------------------------------------------------
    dsl4jsb_Real3D_onDomain :: fract_per_age


    model => get_model(tile%owner_model_id)
    grid => Get_grid(model%grid_id)

    dsl4jsb_Get_config(FAGE_)
    dsl4jsb_Get_memory(FAGE_)

    dsl4jsb_Get_var3D_onDomain(FAGE_, fract_per_age)

    !@todo: needs to be adapted when absolute areas are introduced!!!

    ! Distribution of age fractions - needs to comply with the fractions in the age cohorts (mo_jsb_tile - Init_tile)
    ! TODO-JN-COHORTS: options fromFile; allLast
    IF (dsl4jsb_Config(FAGE_)%age_cohort_area_fracts == 'allFirst') THEN
      ! all fracts to youngest age
      dsl4jsb_var3D_onDomain(FAGE_, fract_per_age) = 0._wp
      dsl4jsb_var_ptr(FAGE_, fract_per_age) (:,1,:) = 1._wp
    ELSE
      CALL finish(TRIM(routine), 'Currently only "allFirst" implemented for the distribution of age fractions and cohorts.')
    ENDIF

    ! input_file = jsb_netcdf_open_input(dsl4jsb_Config(FAGE_)%ic_filename, model%grid_id)
!    ptr_2D(:,:) = tile%fract(:,:)

!----------------------------------------------------------------------------------------
!TODO-JN would only be required when working with absolute areas?
!    ALLOCATE(ptr_2D(grid%nproma, grid%nblks))
!    CALL tile%Get_fraction(fract=ptr_2D(:,:))
!    dsl4jsb_var_ptr(FAGE_, fract_per_age) (:,1,:) = ptr_2D(:,:) !TODO-JN check: arrays seem to not be 0 based?!
!    DEALLOCATE(ptr_2D)
!----------------------------------------------------------------------------------------


    IF (tile%Has_children()) THEN

!      DO isoil=1,nsoil
!        ptr_3D => input_file%Read_2d_extdim(  &
!          & variable_name='layer_moist',                &
!          & start_extdim=isoil, end_extdim=isoil, extdim_name='soillev')
!        dsl4jsb_var_ptr(HYDRO_,w_soil_sl) (:,isoil,:) = &
!          & MERGE(ptr_3D(:,:,1), 0._wp, ptr_3D(:,:,1) >= 0._wp)
!      END DO

      dsl4jsb_var2D_onDomain(FAGE_, mean_age) = get_mean_age(dsl4jsb_var3D_onDomain(FAGE_, fract_per_age))
      !avgAge = util function - age weighted with area

    ELSE
      !Only one age -> for now initialise age with zero later e.g. from file
      dsl4jsb_var2D_onDomain(FAGE_, mean_age) = 0._wp
    ENDIF

  END SUBROUTINE fage_init_ic
  
#endif
END MODULE mo_fage_init
