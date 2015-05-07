!>
!! @brief Initialization of the the hd memory.
!!
!! @author
!!  Veronika Gayler, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Veronika Gayler (2013-04-20)
!!
MODULE mo_hd_init
#ifndef __NO_JSBACH__

  USE mo_jsb_model_class,   ONLY: t_jsb_model
  USE mo_jsb_tile_class,    ONLY: t_jsb_tile
  USE mo_hd_memory_class,   ONLY: t_hd_memory, t_hd_param, t_hd_state
  USE mo_jsb_base,          ONLY: get_model
  USE mo_jsb_io_netcdf,     ONLY: t_input_file, jsb_netcdf_close_file, jsb_netcdf_open_input, &
                                  jsb_netcdf_read_2d, jsb_netcdf_read_2d_time
  USE mo_kind,              ONLY: wp
  USE mo_exception,         ONLY: message

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: hd_init_ic, hd_init_bc    ! initialization of initial and boundary conditions

  CHARACTER(len=*), PARAMETER :: modname = 'mo_hd_init'

CONTAINS

  SUBROUTINE hd_init_bc(tile)

    TYPE(t_jsb_tile),  INTENT(inout) :: tile

    TYPE(t_jsb_model), POINTER       :: model
    TYPE(t_hd_memory), POINTER       :: hd
    TYPE(t_hd_param),  POINTER       :: param
    TYPE(t_hd_state),  POINTER       :: state

    REAL(wp), POINTER :: &
      & return_pointer (:,  :)       !< temporary pointer

    TYPE(t_input_file) :: input_file

    CHARACTER(len=*),  PARAMETER     :: routine = TRIM(modname)//':hd_init_bc'

    model => get_model(tile%owner_model_id)

    ! Get hd memory of the tile (only the root tile is used)
    hd => tile%mem%hd
    IF (.NOT. ASSOCIATED(hd)) RETURN

    CALL message(TRIM(routine), 'Initializing hd memory for tile '//TRIM(tile%name)//' from '//TRIM(model%config%hd%bc_filename))

    ! Initialize boundary conditions
    !
    param => hd%hd_param

    input_file = jsb_netcdf_open_input('bc_hdpara.nc', model%id)     ! Temporary
!!$    input_file = jsb_netcdf_open_input(TRIM(model%config%hd%bc_filename), model%id)

    ! FLAG: HD model land sea mask
    return_pointer => jsb_netcdf_read_2d(         &
      & input_file,                               &
      & variable_name='FLAG',                     &
      & fill_array = param%hd_lsm)

    ! FDIR: river and base flow direction
    return_pointer => jsb_netcdf_read_2d(         &
      & input_file,                               &
      & variable_name='FDIR',                     &
      & fill_array = param%fdir)

    ! FILNEW: longitude of next grid cell downstream
    return_pointer => jsb_netcdf_read_2d(         &
      & input_file,                               &
      & variable_name='FILNEW',                   &
      & fill_array = param%filnew)

    ! FIBNEW: latitude of next grid cell downstream
    return_pointer => jsb_netcdf_read_2d(         &
      & input_file,                               &
      & variable_name='FIBNEW',                   &
      & fill_array = param%fibnew)

    ! ALF_K: retention constant of overland flow [day]
    return_pointer => jsb_netcdf_read_2d(         &
      & input_file,                               &
      & variable_name='ALF_K',                    &
      & fill_array = param%alf_k)

    ! ALF_N: number of overland flow reservoirs
    return_pointer => jsb_netcdf_read_2d(         &
      & input_file,                               &
      & variable_name='ALF_N',                    &
      & fill_array = param%alf_n)

    ! ARF_K: retention constant of river flow [day]
    return_pointer => jsb_netcdf_read_2d(         &
      & input_file,                               &
      & variable_name='ARF_K',                   &
      & fill_array = param%arf_k)

    ! ARF_N: number of river flow reservoirs
    return_pointer => jsb_netcdf_read_2d(         &
      & input_file,                               &
      & variable_name='ARF_N',                    &
      & fill_array = param%arf_n)

    ! AGF_K: retention constant for base flow [day]
    return_pointer => jsb_netcdf_read_2d(         &
      & input_file,                               &
      & variable_name='AGF_K',                    &
      & fill_array = param%agf_k)

    CALL jsb_netcdf_close_file(input_file)

  END SUBROUTINE hd_init_bc

  SUBROUTINE hd_init_ic(tile)


    TYPE(t_jsb_tile),  INTENT(inout) :: tile

    TYPE(t_jsb_model),       POINTER :: model
    TYPE(t_hd_memory),       POINTER :: hd
    TYPE(t_hd_param),        POINTER :: param
    TYPE(t_hd_state),        POINTER :: state

    REAL(wp), POINTER :: &
      & return_pointer (:,  :)       !< temporary pointer

    TYPE(t_input_file) :: input_file

    INTEGER  :: i
    REAL(wp) :: wgt1, wgt2
    INTEGER  :: nmw1,nmw2, n1, n2

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':hd_init_ic'

    model => get_model(tile%owner_model_id)

    ! Get hd memory of the tile
    hd => tile%mem%hd
    IF (.NOT. ASSOCIATED(hd)) RETURN

    ! read initial conditions

    input_file = jsb_netcdf_open_input('ic_hdstart.nc', model%id)     ! Temporary
    CALL message(TRIM(routine), 'Initializing hd memory for tile '//TRIM(tile%name)//' from '//TRIM(model%config%hd%ic_filename))

    ! Initialize boundary conditions
    !
    state => hd%hd_state

    ! FLFMEM
    return_pointer => jsb_netcdf_read_2d(         &
      & input_file,                               &
      & variable_name='FLFMEM',                   &
      & fill_array = state%flfmem)

    

  END SUBROUTINE hd_init_ic

#endif
END MODULE mo_hd_init
