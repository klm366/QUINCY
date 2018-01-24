!>
!! @brief Initialization of the the carbon memory.
!!
!! @author
!!  Rainer Schneck, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Rainer Schneck (2017-04-12)
!!
MODULE mo_carbon_init
#ifndef __NO_JSBACH__

  USE mo_kind,              ONLY: wp
  USE mo_exception,         ONLY: message

  USE mo_jsb_model_class,   ONLY: t_jsb_model
  USE mo_jsb_class,         ONLY: get_model
  USE mo_jsb_tile_class,    ONLY: t_jsb_tile_abstract
  USE mo_jsb_config_class,  ONLY: t_jsb_config_p
  USE mo_jsb_io_netcdf,     ONLY: t_input_file, jsb_netcdf_open_input
  USE mo_jsb_time,          ONLY: get_time_dt, is_time_experiment_start, is_time_restart

  dsl4jsb_Use_processes CARBON_
  dsl4jsb_Use_config(CARBON_)
  dsl4jsb_Use_memory(CARBON_)

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: carbon_init

  CHARACTER(len=*), PARAMETER :: modname = 'mo_carbon_init'

CONTAINS

  !
  !> Intialize carbon process (after memory has been set up)
  !
  SUBROUTINE carbon_init(tile)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_model), POINTER                :: model
    INTEGER                                   :: iblk, ics, ice, nc
    REAL(wp)                                  :: dtime

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':carbon_init'

    dsl4jsb_Def_memory(CARBON_)
    dsl4jsb_Def_config(CARBON_)

    dsl4jsb_Real3D_onChunk :: F_pseudo_temp
    dsl4jsb_Real3D_onChunk :: N_pseudo_temp
    dsl4jsb_Real3D_onChunk :: F_pseudo_precip
    dsl4jsb_Real3D_onChunk :: N_pseudo_precip

    model => get_model(tile%owner_model_id)

    dsl4jsb_Get_memory(CARBON_)
    dsl4jsb_Get_config(CARBON_)

    dtime = get_time_dt()

    !R: Subroutine carbon_init should be called only if experiment is started or restarted. Therefore this is arbitary:
    !   IF (dsl4jsb_Config(CARBON_)%read_cpools   .AND.  (is_time_experiment_start() .OR. is_time_restart()  ) THEN

    IF (dsl4jsb_Config(CARBON_)%read_cpools ) THEN ! This "IF" is necessary because we do not want the model to stop when
                                                   ! there is no ic_land_carbon.nc file.
      CALL carbon_init_ic(tile)
    END IF
    !CALL carbon_init_bc(tile)


    ! R: This should stand here because these variables are more technical constants (only depending on dtime)
    !    than a part of the carbon process parametrization. Furthermore they should be calculated only once when a simulation
    !    is started or restarted and not for each tile and for each time step again!
    !    It is not necessary to write and read them to/from the restart files.
    !    It is not necessary to write them to the output file.
    !    However, we could think about to kick them also out of the memory...
    !    AND we should make them 1D, later, as soon as we do not want them in the output for debugging!

     !dsl4jsb_Get_var2D_onChunk(CARBON_,   F_pseudo_temp)    ! for this we would need iblk, but we want the whole field!
     F_pseudo_temp   => carbon__mem%F_pseudo_temp%ptr(:,:)   ! out
     N_pseudo_temp   => carbon__mem%N_pseudo_temp%ptr(:,:)   ! out
     F_pseudo_precip => carbon__mem%F_pseudo_precip%ptr(:,:) ! out
     N_pseudo_precip => carbon__mem%N_pseudo_precip%ptr(:,:) ! out

    ! R: Note, I changed the "writing style" as compared to JSBACH3.
    ! F_pseudo_temp and N_pseudo_temp are constants used for updating the pseudo-15-day mean temperature.
    ! These constants define the contributions of t_air and pseudo_temp to the new pseudo_temp.
    F_pseudo_temp = EXP(-dtime / (15._wp * 86400._wp))   ! F_pseudo_temp decreases expotentially from 1 with higher
                                                         ! (dtime/15DaysInSeconds). For (dtime/15DaysInSeconds)
                                                         !  =1 it would be 0.368 and =0.5 it would be 0.607.
    N_pseudo_temp = 1._wp - F_pseudo_temp                ! Rises with the same amout as F_pseudo_temp falls.

    F_pseudo_precip = EXP(-dtime / (15._wp * 86400._wp))
    N_pseudo_precip = 1._wp - F_pseudo_precip

  END SUBROUTINE carbon_init

  SUBROUTINE carbon_init_bc(tile)
    !-----------------------------------------------------------------------
    !  ARGUMENTS
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile

    !-----------------------------------------------------------------------
    !  LOCAL VARIABLES
    TYPE(t_jsb_model),       POINTER :: model
    dsl4jsb_Def_memory(CARBON_)

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':carbon_init_bc'

    !-----------------------------------------------------------------------
    ! CONTENT

    model => get_model(tile%owner_model_id)

    ! Get carbon memory of the tile
    dsl4jsb_Get_memory(CARBON_)

    IF (.NOT. ASSOCIATED(carbon__mem)) RETURN

    CALL message(TRIM(routine), 'Setting  carbon boundary conditions for tile '//TRIM(tile%name))

    ! R: nothing to do in the moment ...

  END SUBROUTINE carbon_init_bc

  RECURSIVE SUBROUTINE carbon_init_ic(tile)

    !USE mo_carbon_constants,  ONLY: AlbedoCanopySnow_age, AlbedoCanopySnow_temp, AlbedoVisInitial, AlbedoNirInitial

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile

    TYPE(t_jsb_model),         POINTER :: model
    dsl4jsb_Def_config(CARBON_)
    dsl4jsb_Def_memory(CARBON_)

    !INTEGER :: icarbon
    REAL(wp), POINTER :: &
      & return_pointer  (:,  :) !, & !< temporary pointer
      ! & return_pointer3d(:,:,:)
    !LOGICAL :: &
    !  read_cpools

    TYPE(t_input_file) :: input_file
    !INTEGER :: i, n1, n2

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':carbon_init_ic'

    model => get_model(tile%owner_model_id)

    ! Get carbon memory and config  of the tile
    dsl4jsb_Get_memory(CARBON_)
    dsl4jsb_Get_config(CARBON_)

    IF (.NOT. dsl4jsb_Config(CARBON_)%active) RETURN

    CALL message(TRIM(routine), 'Reading carbon pools in memory for tile '//TRIM(tile%name))

    ! Initialize physical state
    !

    input_file = jsb_netcdf_open_input(TRIM(dsl4jsb_Config(CARBON_)%ic_filename), model%grid_id)

    ! TBD: adapt names of variables in input files to names in memory

    ! R: These variables should be 3D as they are tile dependent. However for THIS tile they would be 2D but 
    !    in the input file they are 3D!
    ! R: TBD: 3D Variablen einlesen. Dann wäre es dxsl4jsb_var3D_onDomain!
    !    Bisher wurden 3D variablen irgendwo über einlesen vom File initialisiert!
    return_pointer => input_file%Read_2d(         &
      & variable_name='Cpool_green',           &
      & fill_array = dsl4jsb_var_ptr(CARBON_,Cpool_green))
    dsl4jsb_var2D_onDomain(CARBON_,Cpool_green) = &
      MERGE(dsl4jsb_var2D_onDomain(CARBON_,Cpool_green), 0._wp, dsl4jsb_var2D_onDomain(CARBON_,Cpool_green) >= 0._wp)

    !
    return_pointer => input_file%Read_2d(         &
      & variable_name='Cpool_woods',             &
      & fill_array = dsl4jsb_var_ptr(CARBON_, Cpool_woods))
    dsl4jsb_var2D_onDomain(CARBON_, Cpool_woods) = &
      MERGE(dsl4jsb_var2D_onDomain(CARBON_,Cpool_woods), 0._wp, dsl4jsb_var2D_onDomain(CARBON_,Cpool_woods) >= 0._wp)

    !
    return_pointer => input_file%Read_2d(         &
      & variable_name='Cpool_reserve',          &
      & fill_array = dsl4jsb_var_ptr(CARBON_, Cpool_reserve))
    dsl4jsb_var2D_onDomain(CARBON_,Cpool_reserve) = &
      MERGE(dsl4jsb_var2D_onDomain(CARBON_,Cpool_reserve), 0._wp, dsl4jsb_var2D_onDomain(CARBON_,Cpool_reserve) >= 0._wp)

    !
    return_pointer => input_file%Read_2d(         &
      & variable_name='YCpool_acid_ag1',          &
      & fill_array = dsl4jsb_var_ptr(CARBON_, YCpool_acid_ag1))
    dsl4jsb_var2D_onDomain(CARBON_,YCpool_acid_ag1) = &
      MERGE(dsl4jsb_var2D_onDomain(CARBON_,YCpool_acid_ag1), 0._wp, dsl4jsb_var2D_onDomain(CARBON_,YCpool_acid_ag1) >= 0._wp)

    !
    return_pointer => input_file%Read_2d(         &
      & variable_name='YCpool_acid_ag1',          &
      & fill_array = dsl4jsb_var_ptr(CARBON_, YCpool_acid_ag1))
    dsl4jsb_var2D_onDomain(CARBON_,YCpool_acid_ag1) = &
      MERGE(dsl4jsb_var2D_onDomain(CARBON_,YCpool_acid_ag1), 0._wp, dsl4jsb_var2D_onDomain(CARBON_,YCpool_acid_ag1) >= 0._wp)
    !
    return_pointer => input_file%Read_2d(         &
      & variable_name='YCpool_ethanol_ag1',          &
      & fill_array = dsl4jsb_var_ptr(CARBON_, YCpool_ethanol_ag1))
    dsl4jsb_var2D_onDomain(CARBON_,YCpool_ethanol_ag1) = &
      MERGE(dsl4jsb_var2D_onDomain(CARBON_,YCpool_ethanol_ag1), 0._wp, dsl4jsb_var2D_onDomain(CARBON_,YCpool_ethanol_ag1) >= 0._wp)
    !
    return_pointer => input_file%Read_2d(         &
      & variable_name='YCpool_nonsoluble_ag1',          &
      & fill_array = dsl4jsb_var_ptr(CARBON_, YCpool_nonsoluble_ag1))
    dsl4jsb_var2D_onDomain(CARBON_,YCpool_nonsoluble_ag1) = &
      MERGE(dsl4jsb_var2D_onDomain(CARBON_,YCpool_nonsoluble_ag1), 0._wp, dsl4jsb_var2D_onDomain(CARBON_,YCpool_nonsoluble_ag1) >= 0._wp)
      
    !
    return_pointer => input_file%Read_2d(         &
      & variable_name='YCpool_acid_bg1',          &
      & fill_array = dsl4jsb_var_ptr(CARBON_, YCpool_acid_bg1))
    dsl4jsb_var2D_onDomain(CARBON_,YCpool_acid_bg1) = &
      MERGE(dsl4jsb_var2D_onDomain(CARBON_,YCpool_acid_bg1), 0._wp, dsl4jsb_var2D_onDomain(CARBON_,YCpool_acid_bg1) >= 0._wp)
    !
    return_pointer => input_file%Read_2d(         &
      & variable_name='YCpool_water_bg1',          &
      & fill_array = dsl4jsb_var_ptr(CARBON_, YCpool_water_bg1))
    dsl4jsb_var2D_onDomain(CARBON_,YCpool_water_bg1) = &
      MERGE(dsl4jsb_var2D_onDomain(CARBON_,YCpool_water_bg1), 0._wp, dsl4jsb_var2D_onDomain(CARBON_,YCpool_water_bg1) >= 0._wp)
    !
    return_pointer => input_file%Read_2d(         &
      & variable_name='YCpool_ethanol_bg1',          &
      & fill_array = dsl4jsb_var_ptr(CARBON_, YCpool_ethanol_bg1))
    dsl4jsb_var2D_onDomain(CARBON_,YCpool_ethanol_bg1) = &
      MERGE(dsl4jsb_var2D_onDomain(CARBON_,YCpool_ethanol_bg1), 0._wp, dsl4jsb_var2D_onDomain(CARBON_,YCpool_ethanol_bg1) >= 0._wp)
    !
    return_pointer => input_file%Read_2d(         &
      & variable_name='YCpool_nonsoluble_bg1',          &
      & fill_array = dsl4jsb_var_ptr(CARBON_, YCpool_nonsoluble_bg1))
    dsl4jsb_var2D_onDomain(CARBON_,YCpool_nonsoluble_bg1) = &
      MERGE(dsl4jsb_var2D_onDomain(CARBON_,YCpool_nonsoluble_bg1), 0._wp, dsl4jsb_var2D_onDomain(CARBON_,YCpool_nonsoluble_bg1) >= 0._wp)
    !
    return_pointer => input_file%Read_2d(         &
      & variable_name='YCpool_humus_1',          &
      & fill_array = dsl4jsb_var_ptr(CARBON_, YCpool_humus_1))
    dsl4jsb_var2D_onDomain(CARBON_,YCpool_humus_1) = &
      MERGE(dsl4jsb_var2D_onDomain(CARBON_,YCpool_humus_1), 0._wp, dsl4jsb_var2D_onDomain(CARBON_,YCpool_humus_1) >= 0._wp)
    !
    return_pointer => input_file%Read_2d(         &
      & variable_name='YCpool_acid_ag2',          &
      & fill_array = dsl4jsb_var_ptr(CARBON_, YCpool_acid_ag2))
    dsl4jsb_var2D_onDomain(CARBON_,YCpool_acid_ag2) = &
      MERGE(dsl4jsb_var2D_onDomain(CARBON_,YCpool_acid_ag2), 0._wp, dsl4jsb_var2D_onDomain(CARBON_,YCpool_acid_ag2) >= 0._wp)
    !
    return_pointer => input_file%Read_2d(         &
      & variable_name='YCpool_water_ag2',          &
      & fill_array = dsl4jsb_var_ptr(CARBON_, YCpool_water_ag2))
    dsl4jsb_var2D_onDomain(CARBON_,YCpool_water_ag2) = &
      MERGE(dsl4jsb_var2D_onDomain(CARBON_,YCpool_water_ag2), 0._wp, dsl4jsb_var2D_onDomain(CARBON_,YCpool_water_ag2) >= 0._wp)
    !
    return_pointer => input_file%Read_2d(         &
      & variable_name='YCpool_ethanol_ag2',          &
      & fill_array = dsl4jsb_var_ptr(CARBON_, YCpool_ethanol_ag2))
    dsl4jsb_var2D_onDomain(CARBON_,YCpool_ethanol_ag2) = &
      MERGE(dsl4jsb_var2D_onDomain(CARBON_,YCpool_ethanol_ag2), 0._wp, dsl4jsb_var2D_onDomain(CARBON_,YCpool_ethanol_ag2) >= 0._wp)
    !
    return_pointer => input_file%Read_2d(         &
      & variable_name='YCpool_nonsoluble_ag2',          &
      & fill_array = dsl4jsb_var_ptr(CARBON_, YCpool_nonsoluble_ag2))
    dsl4jsb_var2D_onDomain(CARBON_,YCpool_nonsoluble_ag2) = &
      MERGE(dsl4jsb_var2D_onDomain(CARBON_,YCpool_nonsoluble_ag2), 0._wp, dsl4jsb_var2D_onDomain(CARBON_,YCpool_nonsoluble_ag2) >= 0._wp)

    !
    return_pointer => input_file%Read_2d(         &
      & variable_name='YCpool_acid_bg2',          &
      & fill_array = dsl4jsb_var_ptr(CARBON_, YCpool_acid_bg2))
    dsl4jsb_var2D_onDomain(CARBON_,YCpool_acid_bg2) = &
      MERGE(dsl4jsb_var2D_onDomain(CARBON_,YCpool_acid_bg2), 0._wp, dsl4jsb_var2D_onDomain(CARBON_,YCpool_acid_bg2) >= 0._wp)
    !
    return_pointer => input_file%Read_2d(         &
      & variable_name='YCpool_water_bg2',          &
      & fill_array = dsl4jsb_var_ptr(CARBON_, YCpool_water_bg2))
    dsl4jsb_var2D_onDomain(CARBON_,YCpool_water_bg2) = &
      MERGE(dsl4jsb_var2D_onDomain(CARBON_,YCpool_water_bg2), 0._wp, dsl4jsb_var2D_onDomain(CARBON_,YCpool_water_bg2) >= 0._wp)
    !
    return_pointer => input_file%Read_2d(         &
      & variable_name='YCpool_ethanol_bg2',          &
      & fill_array = dsl4jsb_var_ptr(CARBON_, YCpool_ethanol_bg2))
    dsl4jsb_var2D_onDomain(CARBON_,YCpool_ethanol_bg2) = &
      MERGE(dsl4jsb_var2D_onDomain(CARBON_,YCpool_ethanol_bg2), 0._wp, dsl4jsb_var2D_onDomain(CARBON_,YCpool_ethanol_bg2) >= 0._wp)
    !
    return_pointer => input_file%Read_2d(         &
      & variable_name='YCpool_nonsoluble_bg2',          &
      & fill_array = dsl4jsb_var_ptr(CARBON_, YCpool_nonsoluble_bg2))
    dsl4jsb_var2D_onDomain(CARBON_,YCpool_nonsoluble_bg2) = &
      MERGE(dsl4jsb_var2D_onDomain(CARBON_,YCpool_nonsoluble_bg2), 0._wp, dsl4jsb_var2D_onDomain(CARBON_,YCpool_nonsoluble_bg2) >= 0._wp)
    !
    return_pointer => input_file%Read_2d(         &
      & variable_name='YCpool_humus_2',          &
      & fill_array = dsl4jsb_var_ptr(CARBON_, YCpool_humus_2))
    dsl4jsb_var2D_onDomain(CARBON_,YCpool_humus_2) = &
      MERGE(dsl4jsb_var2D_onDomain(CARBON_,YCpool_humus_2), 0._wp, dsl4jsb_var2D_onDomain(CARBON_,YCpool_humus_2) >= 0._wp)


    CALL input_file%Close()

  END SUBROUTINE carbon_init_ic

#endif
END MODULE mo_carbon_init
