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
!#ifndef __NO_JSBACH__
#if !defined(__NO_JSBACH__) && !defined(__NO_JSBACH_HD__)

#include "dsl_jsbach.inc"

  USE mo_kind,                ONLY: wp
  USE mo_exception,           ONLY: message, finish, message_text
  USE mo_jsb_model_class,     ONLY: t_jsb_model
  USE mo_jsb_tile_class,      ONLY: t_jsb_tile_abstract
  USE mo_jsb_config_class,    ONLY: t_jsb_config_p
  USE mo_jsb_class,           ONLY: get_model
  USE mo_jsb_io_netcdf,       ONLY: t_input_file, jsb_netcdf_open_input
  USE mo_jsb_grid_class,      ONLY: t_jsb_grid, get_blk_start, get_blk_end, get_col_start, get_col_end
  USE mo_jsb_grid,            ONLY: Get_grid
  USE mo_decomposition_tools, ONLY: get_local_index
  USE mo_sync,                ONLY: global_max
  USE mo_jsb_parallel_iface,  ONLY: get_my_global_mpi_id

  dsl4jsb_Use_processes HD_
  dsl4jsb_Use_config(HD_)
  dsl4jsb_Use_memory(HD_)

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: hd_init_ic, hd_init_bc    ! initialization of initial and boundary conditions

  CHARACTER(len=*), PARAMETER :: modname = 'mo_hd_init'

CONTAINS

  SUBROUTINE hd_init_bc(tile)

    CLASS(t_jsb_tile_abstract),  INTENT(inout) :: tile

    TYPE(t_jsb_model), POINTER       :: model
    TYPE(t_jsb_grid),  POINTER       :: grid

    dsl4jsb_Def_config(HD_)
    dsl4jsb_Def_memory(HD_)

    REAL(wp), POINTER ::         &
      & return_pointer   (:,:),  &   !< temporary pointer
      & return_pointer3d (:,:,:)     !< temporary pointer

    TYPE(t_input_file)  :: input_file

    INTEGER :: n, blk, blks, blke, jc, jcs, jce, idx 
    INTEGER :: nblks, nproma

    CHARACTER(len=*),  PARAMETER     :: routine = TRIM(modname)//':hd_init_bc'

    model => get_model(tile%owner_model_id)
    grid => Get_grid(model%grid_id)

    ! Get hd config
    dsl4jsb_Get_config(HD_)

    IF (.NOT. dsl4jsb_Config(HD_)%active) RETURN

    ! Get hd memory of the tile (only the root tile is used)
    dsl4jsb_Get_memory(HD_)

    IF (TRIM(dsl4jsb_Config(HD_)%routing_scheme) /= 'full') RETURN

    CALL message(TRIM(routine), 'Initializing hd memory for tile '//TRIM(tile%name)//&
                                ' from '//TRIM(dsl4jsb_Config(HD_)%bc_filename))

    ! Initialize boundary conditions
    !
    ! Number of base flow reservoirs
    dsl4jsb_var2D_onDomain(HD_, nres_baseflow) = 1._wp

    ! Maximum number of reservoirs
    hd__mem%nres_o_max = 1
    hd__mem%nres_b_max = 1
    hd__mem%nres_r_max = 5

    ! Maximum number of inflow neighbor cells
    hd__mem%nneigh = 6

    !CALL message(TRIM(routine), 'Initializing HD fractions for tile '//TRIM(tile%name)//&
    !  &                         ' from '//TRIM(dsl4jsb_Config(HD_)%fract_filneame))
    !input_file = jsb_netcdf_open_input('TRIM(dsl4jsb_Config(HD_)%fract_filneame', model%id)

    ! Initialize data from input file
    input_file = jsb_netcdf_open_input(TRIM(dsl4jsb_Config(HD_)%bc_filename), model%grid_id)

    ! HD model mask: -1 ocean, 0: ocean inflow cells, 1: land with outflow, 2: internal drainage
    return_pointer => input_file%Read_2d(         &
      & variable_name='MASK',                     &
      & fill_array = dsl4jsb_var2D_ptr(HD_,hd_mask))

    ! Retention constant of overland flow [day]
    return_pointer => input_file%Read_2d(         &
      & variable_name='ALF_K',                    &
      & fill_array = dsl4jsb_var2D_ptr(HD_,ret_overlflow))

    ! Number of overland flow reservoirs
    return_pointer => input_file%Read_2d(         &
      & variable_name='ALF_N',                    &
      & fill_array = dsl4jsb_var2D_ptr(HD_,nres_overlflow))

    ! Retention constant of river flow [day]
    return_pointer => input_file%Read_2d(         &
      & variable_name='ARF_K',                    &
      & fill_array = dsl4jsb_var2D_ptr(HD_,ret_riverflow))

    ! Number of river flow reservoirs
    return_pointer => input_file%Read_2d(         &
      & variable_name='ARF_N',                    &
      & fill_array = dsl4jsb_var2D_ptr(HD_,nres_riverflow))

    ! Retention constant for base flow [day]
    return_pointer => input_file%Read_2d(         &
      & variable_name='AGF_K',                    &
      & fill_array = dsl4jsb_var2D_ptr(HD_,ret_baseflow))

    ! Global index neighbouring grid cells upstream
    return_pointer3d => input_file%Read_2d_extdim(  &
      & variable_name='CELLS_UP',                   &
      & start_extdim=1, end_extdim=hd__mem%nneigh,    &
      & extdim_name='nneigh')

    ! find local idices of grid cells upstream
    nproma = grid%nproma
    nblks  = grid%nblks
    ALLOCATE(dsl4jsb_var3D_ptr(HD_,nidx_upstream) (nproma, hd__mem%nneigh, nblks))
    ALLOCATE(dsl4jsb_var3D_ptr(HD_,bidx_upstream) (nproma, hd__mem%nneigh, nblks))
    dsl4jsb_var3D_ptr(HD_,nidx_upstream) = -1
    dsl4jsb_var3D_ptr(HD_,bidx_upstream)  = -1

    blks = get_blk_start(grid)
    blke = get_blk_end(grid)
    DO blk = blks, blke
      jcs = get_col_start(blk, grid)
      jce = get_col_end  (blk, grid)
      DO n = 1, hd__mem%nneigh
        DO jc = jcs, jce
          IF (NINT(return_pointer3d(jc,blk,n)) > 0) THEN
            idx = get_local_index(grid%patch%cells%decomp_info%glb2loc_index, NINT(return_pointer3d(jc,blk,n)))
            !print*, 'vg: p_pe: ', p_pe, ' n: ', n, ' blk: ', blk, ' FDIR: ', NINT(return_pointer(n,blk)), 'idx: ',idx
            IF (idx ==  0) THEN
              WRITE (message_text,*) "index of grid cell upstream (", NINT(return_pointer3d(jc,blk,n)), ") not in global domain"
              CALL finish(routine, message_text)
            ELSE IF (idx == -1) THEN
              WRITE (message_text,*) "index of grid cell upstream (", NINT(return_pointer3d(jc,blk,n)), ") not in local domain"
              CALL finish(routine, message_text)
            END IF
            hd__mem%bidx_upstream%ptr(jc,n,blk) = (idx-1)/nproma + 1
            hd__mem%nidx_upstream%ptr(jc,n,blk) = idx - (hd__mem%bidx_upstream%ptr(jc,n,blk)-1)*nproma
          END IF
        END DO
      END DO
    END DO

    CALL input_file%Close()

  END SUBROUTINE hd_init_bc

  SUBROUTINE hd_init_ic(tile)

    CLASS(t_jsb_tile_abstract),  INTENT(inout) :: tile

    dsl4jsb_Def_config(HD_)
    dsl4jsb_Def_memory(HD_)

    REAL(wp), POINTER :: &
      & tile_fract(:,:), &
      & ptr(:,:)            !< temporary pointer

    TYPE(t_input_file) :: input_file

    TYPE(t_jsb_model), POINTER :: model
    TYPE(t_jsb_grid),  POINTER :: grid

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':hd_init_ic'

    model => get_model(tile%owner_model_id)
    grid => Get_grid(model%grid_id)

    ! Get hd config
    dsl4jsb_Get_config(HD_)

    IF (.NOT. dsl4jsb_Config(HD_)%active) RETURN

    ! Get hd memory of the tile (only the root tile is used)
    dsl4jsb_Get_memory(HD_)

    IF (TRIM(dsl4jsb_Config(HD_)%routing_scheme) /= 'full') RETURN

    ! read initial conditions

    CALL message(TRIM(routine), 'Initializing hd memory from '//TRIM(dsl4jsb_Config(HD_)%ic_filename))

!vg     input_file = jsb_netcdf_open_input(TRIM(model%config%hd%ic_filename), model%grid_id)

    ! Initialize boundary conditions
    !

!vg     ! Overlandflow reservoir
!vg     ptr => input_file%Read_2d_extdim(  &
!vg       & variable_name='FLFMEM',                   &
!vg       & start_extdim=1,                           &
!vg       & end_extdim=hd_mem%nres_o_max,             &
!vg       & fill_array = hd__mem%overlflow_res%ptr)
!vg 
!vg     ! Baseflow reservoir
!vg     ptr => input_file%Read_2d_extdim(  &
!vg       & variable_name='FGMEM',                    &
!vg       & start_extdim=1,                           &
!vg       & end_extdim=hd%hd_hd_mem%nres_b_max,        &
!vg       & fill_array = hd_mem%baseflow_res%ptr)
!vg 
!vg     ! Riverflow reservoir
!vg     ptr => input_file%Read_2d_extdim(  &
!vg       & variable_name='FRFMEM',                   &
!vg       & start_extdim=1,                           &
!vg       & end_extdim=hd%hd_hd_mem%nres_r_max,        &
!vg       & fill_array = hd_mem%riverflow_res%ptr)

    ALLOCATE(tile_fract(grid%nproma, grid%nblks))
    CALL tile%Get_fraction(fract=tile_fract(:,:))
    dsl4jsb_var3D_onDomain(HD_,overlflow_res) = 1000._wp * SPREAD(tile_fract(:,:), DIM=2, NCOPIES=hd__mem%nres_o_max)
    dsl4jsb_var3D_onDomain(HD_,baseflow_res)  = 5000._wp * SPREAD(tile_fract(:,:), DIM=2, NCOPIES=hd__mem%nres_b_max)
    dsl4jsb_var3D_onDomain(HD_,riverflow_res) = 1000._wp * SPREAD(tile_fract(:,:), DIM=2, NCOPIES=hd__mem%nres_r_max)

  END SUBROUTINE hd_init_ic

#endif
END MODULE mo_hd_init
