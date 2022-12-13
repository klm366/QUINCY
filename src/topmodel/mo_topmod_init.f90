!>
!! @brief Initialization of the the soil and snow memory.
!!
!! @author
!!  Stiig Wilkenskjeld, MPI-M Hamburg
!!
!! @par Revision History
!! Adopted from JSBACH3 code by Thomas Kleinen, 2022-12-06
!!
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/license .
!!

MODULE mo_top_init
#ifndef __NO_JSBACH__

  USE mo_kind,               ONLY: wp
  USE mo_exception,          ONLY: message
  USE mo_jsb_control,        ONLY: debug_on

  USE mo_jsb_model_class,    ONLY: t_jsb_model
  USE mo_jsb_grid_class,     ONLY: t_jsb_grid, t_jsb_vgrid
  USE mo_jsb_grid,           ONLY: Get_grid, Get_vgrid
  USE mo_jsb_tile_class,     ONLY: t_jsb_tile_abstract
  USE mo_jsb_class,          ONLY: get_model
  USE mo_jsb_io_netcdf,      ONLY: t_input_file, jsb_netcdf_open_input
  USE mo_jsb_io,             ONLY: missval
  USE mo_jsb_impl_constants, ONLY: ifs_nsoil, ifs_soil_depth
  USE mo_util,               ONLY: int2string

  dsl4jsb_Use_processes TOP_, SSE_, SEB_, HYDRO_
  dsl4jsb_Use_config(TOP_)
  dsl4jsb_Use_memory(TOP_)

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: topmodel_init !, init_soil_properties

  TYPE t_top_init_vars
!    REAL(wp), POINTER ::     &
  END TYPE t_top_init_vars

  TYPE(t_top_init_vars) :: top_init_vars

  CHARACTER(len=*), PARAMETER :: modname = 'mo_topmodel_init'

CONTAINS

  !
  !> Intialize soil and snow energy process (after memory has been set up)
  !
  SUBROUTINE topmodel_init(tile)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile

    TYPE(t_jsb_model), POINTER  :: model
    CHARACTER(len=*), PARAMETER :: routine = modname//':topmodel_init'
    
    model => Get_model(tile%owner_model_id)

    IF (.NOT. ASSOCIATED(tile%parent_tile)) THEN
      CALL topmodel_read_init_vars(tile)
    END IF

    CALL topmodel_init_bc(tile)
    CALL topmodel_init_ic(tile)

    IF (tile%Is_last_process_tile(SSE_)) THEN
      CALL topmodel_finalize_init_vars()
    END IF

  END SUBROUTINE topmodel_init

  SUBROUTINE topmodel_finalize_init_vars

!    DEALLOCATE(                     &
!      & sse_init_vars%fao,          &
!      & )

!   IF (ASSOCIATED(sse_init_vars%ifs_tskin)) DEALLOCATE(sse_init_vars%ifs_tskin)

  END SUBROUTINE topmodel_finalize_init_vars
      
  SUBROUTINE topmodel_read_init_vars(tile)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile

    dsl4jsb_Def_config(TOP_)

    REAL(wp), POINTER :: ptr_2D(:,:), ptr_3D(:,:,:) !< temporary pointers

    TYPE(t_input_file) :: input_file

    TYPE(t_jsb_model), POINTER :: model
    TYPE(t_jsb_grid),  POINTER :: grid
    TYPE(t_jsb_vgrid), POINTER :: soil_e
    INTEGER :: nproma, nblks, isoil

    CHARACTER(len=*), PARAMETER :: routine = modname//':topmodel_read_init_vars'

    model => get_model(tile%owner_model_id)
    grid  => get_grid(model%grid_id)
    nproma = grid%Get_nproma()
    nblks  = grid%Get_nblks()

!    ALLOCATE(                                          &
!      & sse_init_vars%fao         (nproma, nblks    ), &
!      & )

!    sse_init_vars%fao         (:,:)   = missval

    ! Get soil config
    dsl4jsb_Get_config(TOP_)

    IF (.NOT. tile%Is_process_active(TOP_)) RETURN

    IF (debug_on()) CALL message(TRIM(routine), 'Reading soil energy boundary conditions for tile '//TRIM(tile%name))

    input_file = jsb_netcdf_open_input(TRIM(dsl4jsb_Config(TOP_)%bc_filename), model%grid_id)

    ! Thermal diffusivity of soil
!    ptr_2D => input_file%Read_2d(                    &
!      & variable_name='fao')
!    ptr_2D(:,:) = MERGE(ptr_2D(:,:), 0._wp, ptr_2D(:,:) >= 0._wp)

    CALL input_file%Close()

    input_file = jsb_netcdf_open_input(dsl4jsb_Config(TOP_)%ic_filename, model%grid_id)

    CALL input_file%Close()

  END SUBROUTINE topmodel_read_init_vars

  SUBROUTINE topmodel_init_bc(tile)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile

    dsl4jsb_Def_config(TOP_)
    dsl4jsb_Def_memory(TOP_)

    INTEGER :: iblk

    TYPE(t_jsb_model), POINTER :: model

    CHARACTER(len=*), PARAMETER :: routine = modname//':topmodel_init_bc'

    model => get_model(tile%owner_model_id)

    ! Get soil config
    dsl4jsb_Get_config(TOP_)

    IF (.NOT. tile%Is_process_active(TOP_)) RETURN

    IF (debug_on()) CALL message(TRIM(routine), 'Setting top-model boundary conditions for tile '//TRIM(tile%name))

    ! Get soil memory of the tile
    dsl4jsb_Get_memory(TOP_)

    ! Volumetric heat capacity of dry soil [J/(m^3*K)]
!    IF (dsl4jsb_Config(SSE_)%l_heat_cap_map) THEN
!      dsl4jsb_var_ptr(SSE_, vol_heat_cap) = sse_init_vars%vol_heat_cap
!    ELSE
!      DO iblk=1,SIZE(sse_init_vars%fao,2)
!        dsl4jsb_var_ptr(SSE_,vol_heat_cap)(:,iblk) = fao_vol_hcap(NINT(sse_init_vars%fao(:,iblk)))
!      END DO
!    END IF

    !$acc data present(dsl4jsb_var2D_onDomain(SSE_,thermal_diffusivity), &
    !$acc dsl4jsb_var2D_onDomain(SSE_,vol_heat_cap), &
    !$acc update device(dsl4jsb_var2D_onDomain(SSE_,thermal_diffusivity), &
    !$acc dsl4jsb_var2D_onDomain(SSE_,vol_heat_cap), &
    !$acc end data

  END SUBROUTINE topmodel_init_bc

  SUBROUTINE topmodel_init_ic(tile)

!    USE mo_sse_process, ONLY: init_soil_temperature
!    USE mo_util,        ONLY: ifs2soil

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile

    TYPE(t_jsb_model), POINTER :: model
    TYPE(t_jsb_vgrid), POINTER :: soil_e

    dsl4jsb_Def_config(TOP_)
    dsl4jsb_Def_memory(TOP_)

    INTEGER :: iblk

    CHARACTER(len=*), PARAMETER :: routine = modname//':topmodel_init_ic'

    model => get_model(tile%owner_model_id)

    ! Get soil memory of the tile
    dsl4jsb_Get_memory(TOP_)

    ! Get soil config
    dsl4jsb_Get_config(TOP_)

    IF (.NOT. tile%Is_process_active(TOP_)) RETURN

    IF (debug_on()) CALL message(TRIM(routine), 'Initializing top-model  memory for tile '//TRIM(tile%name))

    ! Soil and surface temperatures
    !
!    IF (model%config%init_from_ifs) THEN
!
!      CALL ifs2soil(sse_init_vars%ifs_tsoil(:,:,:), ifs_soil_depth(:), dsl4jsb_var_ptr(SSE_,t_soil_sl)(:,:,:), soil_e%ubounds(:))
!
!    ELSE
!      DO iblk=1,SIZE(sse_init_vars%t_srf_mon, DIM=2)
!        CALL init_soil_temperature(                   &
!          & soil_e%levels(:),                         &
!          & sse_init_vars%t_srf_mon(:,iblk,1:12),     &
!          & dsl4jsb_var_ptr(SSE_,t_soil_sl)(:,:,iblk))
!      END DO
!    END IF

    !$acc data present(dsl4jsb_var3D_onDomain(SSE_,t_soil_sl))
    !$acc update device(dsl4jsb_var3D_onDomain(SSE_,t_soil_sl))
    !$acc end data

!    IF (TRIM(tile%name) == 'land') THEN
!      dsl4jsb_var2D_onDomain(SEB_,t)        = dsl4jsb_var_ptr(SSE_,t_soil_sl)(:,1,:)

      !$acc data present(dsl4jsb_var2D_onDomain(SEB_,t_unfilt), & 
      !$acc              dsl4jsb_var2D_onDomain(SEB_,t))
      !$acc update device(dsl4jsb_var2D_onDomain(SEB_,t_unfilt), & 
      !$acc               dsl4jsb_var2D_onDomain(SEB_,t))
      !$acc end data
!    END IF

  END SUBROUTINE topmodel_init_ic

#endif
END MODULE mo_top_init
