!>
!! @brief Initialization of the fuel memory.
!!
!! @author
!!  Rainer Schneck, MPI-M Hamburg
!!  Stiig Wilkenskjeld, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Rainer Schneck (2017-06-20)
!! 2018-04-26, StW: Added Spitfire support
!!
MODULE mo_fuel_init
#ifndef __NO_JSBACH__

  USE mo_kind,              ONLY: wp
  USE mo_exception,         ONLY: message
  USE mo_jsb_control,       ONLY: debug_on
  USE mo_jsb_model_class,   ONLY: t_jsb_model
  USE mo_jsb_class,         ONLY: Get_model
  USE mo_jsb_tile_class,    ONLY: t_jsb_tile_abstract
  USE mo_jsb_config_class,  ONLY: t_jsb_config_p
  USE mo_jsb_io_netcdf,     ONLY: t_input_file, jsb_netcdf_open_input
  USE mo_jsb_time_iface,    ONLY: is_time_experiment_start, is_time_restart

  ! Use of processes in this module
  dsl4jsb_Use_processes FUEL_
  dsl4jsb_Use_processes DISTURB_

  ! Use process configurations
  dsl4jsb_Use_config(FUEL_)
  dsl4jsb_Use_config(DISTURB_)

  ! Use process memory
  dsl4jsb_Use_memory(FUEL_)

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: fuel_init

  CHARACTER(len=*), PARAMETER :: modname = 'mo_fuel_init'

CONTAINS

  !
  !> Intialize fuel process (after memory has been set up)
  !
  SUBROUTINE fuel_init(tile)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':fuel_init'

    dsl4jsb_Def_memory(FUEL_)
    dsl4jsb_Get_memory(FUEL_)

    CALL fuel_init_ic(tile)

    CALL fuel_init_bc(tile)

  END SUBROUTINE fuel_init

  SUBROUTINE fuel_init_bc(tile)

    USE mo_disturbance_constants, ONLY: FIRE_NONE, FIRE_JSB3, FIRE_SPITFIRE
#ifdef HAVE_MOFFLEDAC
    USE mo_input_strings, ONLY: GetKeyedString
    USE mo_input,         ONLY: InputOptAdd, moffledac
#endif

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile

    dsl4jsb_Def_config(FUEL_)
    dsl4jsb_Def_config(DISTURB_)

    TYPE(t_jsb_model),       POINTER :: model
    dsl4jsb_Def_memory(FUEL_)

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':fuel_init_bc'

    !-----------------------------------------------------------------------
    ! CONTENT

    model => get_model(tile%owner_model_id)

    dsl4jsb_Get_config(FUEL_)
    dsl4jsb_Get_config(DISTURB_)

    ! Get fuel memory of the tile
    dsl4jsb_Get_memory(FUEL_)

    IF (.NOT. tile%Is_process_active(FUEL_)) RETURN

    IF (debug_on()) CALL message(TRIM(routine), 'Setting fuel boundary conditions for tile '//TRIM(tile%name))

    SELECT CASE (dsl4jsb_Config(DISTURB_)%fire_algorithm)
      CASE (FIRE_NONE)
        ! By definition nothing to do
      CASE (FIRE_JSB3)
        ! Needs nothing
      CASE (FIRE_SPITFIRE)
#ifdef HAVE_MOFFLEDAC
        IF (.NOT. tile%Is_process_active(FUEL_)) RETURN
        IF (LEN_TRIM(dsl4jsb_Config(DISTURB_)%fire_spitfire_master_tile) == 0) THEN ! Should now be obsolete...
          CALL InputOptAdd(var_name='fuel_population_density_'//TRIM(tile%name), &
            file_name=TRIM(GetKeyedString('pop_dens',TRIM(dsl4jsb_Config(DISTURB_)%bc_filename))), &
            var_file_name='population_density',dt_update_unit='1m',valid_time='ACTUAL',reads_per=1,lstep_by_step='on')
          CALL InputOptAdd(var_name='fuel_lightning_frq_'//TRIM(tile%name), &
            file_name=TRIM(GetKeyedString('lightning',TRIM(dsl4jsb_Config(DISTURB_)%bc_filename))), &
            var_file_name='HRMC_COM_FR',dt_update_unit='1d',valid_time='MIDPOINT',action_cycle='all',reads_per=1)
          CALL InputOptAdd(var_name='fuel_lightning_c2g_'//TRIM(tile%name),var_file_name='c2g', &
            file_name=TRIM(GetKeyedString('lightning',TRIM(dsl4jsb_Config(DISTURB_)%bc_filename))), &
            dt_update_unit='2d',valid_time='MIDPOINT',action_cycle='all',def_val=-1._wp,action_var_miss='warn',reads_per=1)
          ! StW: Would make sense to make a_nd time dependent too, but for now we only have a static field
          CALL InputOptAdd(var_name='fuel_a_nd_'//TRIM(tile%name),var_file_name='a_nd', &
            file_name=TRIM(GetKeyedString('a_nd',TRIM(dsl4jsb_Config(DISTURB_)%bc_filename))), &
            dt_update=0,valid_time='MIDPOINT',action_cycle='all',def_val=1._wp,action_var_miss='warn')
          dsl4jsb_Config(DISTURB_)%fire_spitfire_master_tile = TRIM(tile%name)
        ENDIF
#endif
    END SELECT

  END SUBROUTINE fuel_init_bc

  RECURSIVE SUBROUTINE fuel_init_ic(tile)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile

    TYPE(t_jsb_model),         POINTER :: model

    dsl4jsb_Def_config(FUEL_)
    dsl4jsb_Def_memory(FUEL_)

    REAL(wp), POINTER :: &
      & return_pointer  (:,  :) !, & !< temporary pointer
      ! & return_pointer3d(:,:,:)

    TYPE(t_input_file) :: input_file

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':fuel_init_ic'

    model => get_model(tile%owner_model_id)

    ! Get fuel memory and config  of the tile
    dsl4jsb_Get_memory(FUEL_)
    dsl4jsb_Get_config(FUEL_)

    IF (.NOT. tile%Is_process_active(FUEL_)) RETURN

    CALL message(TRIM(routine), 'Setting fuel initial conditions for tile '//TRIM(tile%name))

    ! Initialize physical state
    !
!    input_file = jsb_netcdf_open_input(TRIM(dsl4jsb_Config(FUEL_)%ic_filename), model%grid_id)

!    CALL input_file%Close()

  END SUBROUTINE fuel_init_ic

#endif
END MODULE mo_fuel_init
