!>
!> Initialization of the methane production and transport module.
!>
!>
!> ICON-Land
!>
!> ---------------------------------------
!> Copyright (C) 2013-2024, MPI-M, MPI-BGC
!>
!> Contact: icon-model.org
!> Authors: AUTHORS.md
!> See LICENSES/ for license information
!> SPDX-License-Identifier: BSD-3-Clause
!> ---------------------------------------
!>
!>
!> author
!>  Stiig Wilkenskjeld, MPI-M Hamburg
!>
!> Revision History
!> Adopted from JSBACH3 code by Thomas Kleinen, 2023-03-17
!>

MODULE mo_methane_init
#ifndef __NO_JSBACH__

  USE mo_kind,                ONLY: lp
  USE mo_exception,           ONLY: message
  USE mo_jsb_control,         ONLY: debug_on
  USE mo_jsb_model_class,     ONLY: t_jsb_model
  USE mo_jsb_grid_class,      ONLY: t_jsb_grid, t_jsb_vgrid
  USE mo_jsb_grid,            ONLY: Get_grid, Get_vgrid
  USE mo_jsb_tile_class,      ONLY: t_jsb_tile_abstract
  USE mo_jsb_class,           ONLY: get_model
  USE mo_jsb_io_netcdf,       ONLY: t_input_file, jsb_netcdf_open_input
  USE mo_jsb_io_netcdf_iface, ONLY: nf !, nf_inq_dimid, nf_inq_dimlen, nf_get_vara_double, nf_inq_varndims, nf_inq_varid
  USE mo_jsb_parallel,        ONLY: my_process_is_stdio, p_io, mpi_comm, p_bcast
  USE mo_ch4_config_class,    ONLY: t_ch4_source_config
  USE mo_ch4_memory_class,    ONLY: t_timeseries
  USE mo_methane_constants,   ONLY: t_specie_const, CH4_consts, CO2_consts, O2_consts

  dsl4jsb_Use_processes HYDRO_, CH4_
  dsl4jsb_Use_config(CH4_)
  dsl4jsb_Use_memory(CH4_)

  IMPLICIT NONE
  PRIVATE
include "netcdf.inc"
  PUBLIC :: methane_init

  TYPE t_ch4_init_vars
    TYPE (t_timeseries) :: &
      & CH4,               &
      & CO2,               &
      &  O2
  END TYPE t_ch4_init_vars

  TYPE(t_ch4_init_vars) :: ch4_init_vars

  CHARACTER(len=*), PARAMETER :: modname = 'mo_methane_init'

CONTAINS

  !
  !> Intialize methane production and transport process (after memory has been set up)
  !
  SUBROUTINE methane_init(tile)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile

    TYPE(t_jsb_model), POINTER  :: model
    CHARACTER(len=*), PARAMETER :: routine = modname//':methane_init'

    model => Get_model(tile%owner_model_id)

    IF (.NOT. ASSOCIATED(tile%parent_tile)) THEN
      NULLIFY(ch4_init_vars%ch4%atm_vmr)
      NULLIFY(ch4_init_vars%co2%atm_vmr)
      NULLIFY(ch4_init_vars%o2%atm_vmr)
      CALL methane_read_init_vars(tile)
    END IF

    CALL methane_init_bc(tile)
    CALL methane_init_ic(tile)

    IF (tile%Is_last_process_tile(CH4_)) THEN
      CALL methane_finalize_init_vars()
    END IF

  END SUBROUTINE methane_init

  SUBROUTINE methane_finalize_init_vars

    IF (ASSOCIATED(ch4_init_vars%ch4%atm_vmr)) DEALLOCATE(ch4_init_vars%ch4%atm_vmr)
    IF (ASSOCIATED(ch4_init_vars%co2%atm_vmr)) DEALLOCATE(ch4_init_vars%co2%atm_vmr)
    IF (ASSOCIATED(ch4_init_vars%o2%atm_vmr )) DEALLOCATE(ch4_init_vars%o2%atm_vmr)

  END SUBROUTINE methane_finalize_init_vars

  SUBROUTINE methane_read_GHG_timeseries(input_file,source,varname,fyears,year_start,year_end,GHG_scale,GAS,ts)
    USE mo_methane_util, ONLY: convert_2_vmr

    TYPE(t_input_file),   INTENT(inout) :: input_file
    INTEGER,              INTENT(in)    :: source
    CHARACTER(len=*),     INTENT(in)    :: varname
    REAL(lp),             INTENT(in)    :: fyears(:)
    INTEGER,              INTENT(in)    :: year_start
    INTEGER,              INTENT(in)    :: year_end
    REAL(lp),             INTENT(in)    :: GHG_scale
    TYPE(t_specie_const), INTENT(in)    :: GAS
    TYPE(t_timeseries),   INTENT(inout) :: ts

    CHARACTER(len=*), PARAMETER :: routine = modname//':methane_read_GHG_timeseries'

    INTEGER :: i, nc_id, ndim
    INTEGER :: ind_start(1), ind_end(1), istart(7), isize(7)
    REAL(lp) :: tmp

    IF (my_process_is_stdio()) THEN
      ind_start = MINLOC(ABS(INT(fyears)-year_start))
      ind_end   = MINLOC(ABS(INT(fyears)-year_end))
      ts%atm_year_min = year_start
      ts%atm_year_max = year_end
      IF (INT(fyears(ind_start(1))) > year_start) ts%atm_year_min = INT(fyears(ind_start(1)))
      IF (INT(fyears(ind_end(1)))   < year_end)   ts%atm_year_max = INT(fyears(ind_end(1)))
      ! TODO: This expects the variable in the file to be a global scalar per year
      CALL nf(nf_inq_varid(input_file%file_id, varname, nc_id), routine)
      istart(:) = 1
      isize (:) = 1
      CALL nf(nf_inq_varndims(input_file%file_id,nc_id,ndim),routine)
      ALLOCATE(ts%atm_vmr(ts%atm_year_max-ts%atm_year_min+1))
      istart(ndim) = ind_start(1)
      isize (ndim) = SIZE(ts%atm_vmr)
      CALL nf(nf_get_vara_double(input_file%file_id,nc_id,istart(1:ndim),isize(1:ndim),ts%atm_vmr),routine)
      ts%atm_vmr(:) = ts%atm_vmr(:) * GHG_scale
      SELECT CASE (source)
        CASE (3) ! Already vmr - don't do anything
        CASE (4) ! mmr - convert
          DO i = 1, SIZE(ts%atm_vmr)
            tmp = ts%atm_vmr(i)
            CALL convert_2_vmr(GAS,tmp,ts%atm_vmr(i))
          ENDDO
        CASE (5) ! concentration - convert
          DO i = 1, SIZE(ts%atm_vmr)
            tmp = ts%atm_vmr(i)
            CALL convert_2_vmr(GAS,vmr=ts%atm_vmr(i),conc=tmp)
          ENDDO
      END SELECT
    ENDIF
    CALL p_bcast(ts%atm_year_min,p_io,mpi_comm)
    CALL p_bcast(ts%atm_year_max,p_io,mpi_comm)
    IF (.NOT. my_process_is_stdio()) ALLOCATE(ts%atm_vmr(ts%atm_year_max-ts%atm_year_min+1))
    CALL p_bcast(ts%atm_vmr,     p_io,mpi_comm)

  END SUBROUTINE methane_read_GHG_timeseries

  SUBROUTINE methane_copy_ts(init,GHG_atm)
    TYPE(t_timeseries), INTENT(in)    :: init
    TYPE(t_timeseries), INTENT(inout) :: GHG_atm

    CHARACTER(len=*), PARAMETER :: routine = modname//':methane_copy_ts'

    GHG_atm%atm_year_min = init%atm_year_min
    GHG_atm%atm_year_max = init%atm_year_max
    ALLOCATE(GHG_atm%atm_vmr(SIZE(init%atm_vmr)))
    GHG_atm%atm_vmr(:) = init%atm_vmr(:)
    !$acc data present(GHG_atm%atm_vmr)
    !$acc update device(GHG_atm%atm_vmr)
    !$acc end data

  END SUBROUTINE methane_copy_ts


  SUBROUTINE methane_read_init_vars(tile)
    USE mo_time_config,         ONLY: time_config
    USE mo_jsb_time,            ONLY: get_date_components
    USE mo_util,                ONLY: int2string, real2string

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile

    dsl4jsb_Def_config(CH4_)

    TYPE(t_jsb_model), POINTER :: model
    TYPE(t_jsb_grid),  POINTER :: grid
    TYPE(t_jsb_vgrid), POINTER :: soil_e
    INTEGER :: nproma, nblks, isoil

    TYPE(t_input_file) :: input_file
    INTEGER :: iblk
    INTEGER :: year_start, year_end, ind_start(1), ind_end(1), time_length
    INTEGER :: istart(7), isize(7), ndim
    INTEGER :: nc_id
    LOGICAL :: lfile_source
    REAL(lp), POINTER :: fyears(:)

    CHARACTER(len=*), PARAMETER :: routine = modname//':methane_read_init_vars'

    model => get_model(tile%owner_model_id)
    grid  => get_grid(model%grid_id)
    nproma = grid%Get_nproma()
    nblks  = grid%Get_nblks()

    ! Get methane config
    dsl4jsb_Get_config(CH4_)

    IF (.NOT. tile%Is_process_active(CH4_)) RETURN

    IF (debug_on() .OR. dsl4jsb_Config(CH4_)%dbg_process) &
      CALL message(TRIM(routine), 'Reading methane boundary conditions for tile '//TRIM(tile%name))

    lfile_source = (dsl4jsb_Config(CH4_)%CH4_source > 2) .OR. &
      &            (dsl4jsb_Config(CH4_)%CO2_source > 2) .OR. &
      &            (dsl4jsb_Config(CH4_)%O2_source  > 2)

    NULLIFY(fyears)
    IF (lfile_source) THEN
      CALL get_date_components(time_config%tc_startdate,year=year_start)
      CALL get_date_components(time_config%tc_stopdate, year=year_end)
      IF (my_process_is_stdio()) THEN
        CALL message(routine,'Getting greenhouse gas concentrations '// &
                             int2string(year_start)//' - '//int2string(year_end)//' from file')
        input_file = jsb_netcdf_open_input(TRIM(dsl4jsb_Config(CH4_)%bc_filename))
        ! TODO: This expects the "time" in the file to be 1d and given as years only (integer/floating point)
        CALL nf(nf_inq_dimid(input_file%file_id, 'time', nc_id), routine)
        CALL nf(nf_inq_dimlen(input_file%file_id, nc_id, time_length), routine)
        ALLOCATE(fyears(time_length)) ! Ok to keep on CPU
        CALL nf(nf_inq_varid(input_file%file_id, 'time', nc_id), routine)
        CALL nf(nf_get_vara_double(input_file%file_id,nc_id,(/1/),(/time_length/),fyears), routine)
      ENDIF
      IF (dsl4jsb_Config(CH4_)%CH4_source > 2) THEN
        CALL methane_read_GHG_timeseries(input_file,dsl4jsb_Config(CH4_)%CH4_source,dsl4jsb_Config(CH4_)%CH4_varname, &
                                         fyears,year_start,year_end,dsl4jsb_Config(CH4_)%CH4_scale,CH4_consts,ch4_init_vars%ch4)
      ENDIF
      IF (dsl4jsb_Config(CH4_)%CO2_source > 2) THEN
        CALL methane_read_GHG_timeseries(input_file,dsl4jsb_Config(CH4_)%CO2_source,dsl4jsb_Config(CH4_)%CO2_varname, &
                                         fyears,year_start,year_end,dsl4jsb_Config(CH4_)%CO2_scale,CO2_consts,ch4_init_vars%co2)
      ENDIF
      IF (dsl4jsb_Config(CH4_)%O2_source > 2) THEN
        CALL methane_read_GHG_timeseries(input_file,dsl4jsb_Config(CH4_)%O2_source,dsl4jsb_Config(CH4_)%O2_varname, &
                                         fyears,year_start,year_end,dsl4jsb_Config(CH4_)%O2_scale,O2_consts,ch4_init_vars%o2)
      ENDIF
      IF (my_process_is_stdio()) CALL input_file%Close()
    ENDIF
    IF (ASSOCIATED(fyears)) DEALLOCATE(fyears)

  END SUBROUTINE methane_read_init_vars

  SUBROUTINE methane_init_bc(tile)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile

    dsl4jsb_Def_config(CH4_)
    dsl4jsb_Def_memory(CH4_)

    TYPE(t_ch4_source_config), POINTER :: curr
    TYPE(t_jsb_model), POINTER :: model

    CHARACTER(len=*), PARAMETER :: routine = modname//':methane_init_bc'

    model => get_model(tile%owner_model_id)

    ! Get methane config
    dsl4jsb_Get_config(CH4_)
    dsl4jsb_Get_memory(CH4_)

    IF (.NOT. tile%Is_process_active(CH4_)) RETURN

    IF (debug_on() .OR. dsl4jsb_Config(CH4_)%dbg_process) &
      CALL message(TRIM(routine), 'Setting methane boundary conditions for tile '//TRIM(tile%name))

    ! Copy atmospheric concentrations from GHG-file if needed
    IF (dsl4jsb_Config(CH4_)%CH4_source > 2) CALL methane_copy_ts(ch4_init_vars%CH4,dsl4jsb_Memory(CH4_)%GHG_atm_CH4)
    IF (dsl4jsb_Config(CH4_)%CO2_source > 2) CALL methane_copy_ts(ch4_init_vars%CO2,dsl4jsb_Memory(CH4_)%GHG_atm_CO2)
    IF (dsl4jsb_Config(CH4_)%O2_source  > 2) CALL methane_copy_ts(ch4_init_vars%O2 ,dsl4jsb_Memory(CH4_)%GHG_atm_O2)

    ! Associate simple process data to atmospheric coupling if requested as temporary placeholders to avoid accidential allocation
    curr => dsl4jsb_Config(CH4_)%ch4_sources
    DO WHILE (ASSOCIATED(curr))
      SELECT CASE (curr%isource)
        CASE (5) ! Termites
          curr%dta => dsl4jsb_var2d_onDomain(CH4_,CH4_out_termite)
        CASE (6) ! Herbivores
          curr%dta => dsl4jsb_var2d_onDomain(CH4_,CH4_out_herbivore)
        CASE (7) ! Fires
          curr%dta => dsl4jsb_var2d_onDomain(CH4_,CH4_out_fire)
      END SELECT
      curr => curr%next
    ENDDO

  END SUBROUTINE methane_init_bc

  SUBROUTINE methane_init_ic(tile)

    USE mo_methane_util, ONLY: find_bedrock_layer, find_layer
    USE mo_methane_process, ONLY: determine_root_dist
    dsl4jsb_Use_memory(HYDRO_)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile

    TYPE(t_jsb_model), POINTER :: model
    TYPE(t_jsb_grid),  POINTER :: grid
    TYPE(t_jsb_vgrid), POINTER :: layer

    dsl4jsb_Def_config(CH4_)
    dsl4jsb_Def_memory(CH4_)
    dsl4jsb_Def_memory(HYDRO_)

    INTEGER :: nblks, i, j, nproma

    CHARACTER(len=*), PARAMETER :: routine = modname//':methane_init_ic'

    dsl4jsb_Real2D_onDomain :: soil_depth
    dsl4jsb_Real2D_onDomain :: mxlay_above_rock
    dsl4jsb_Real3D_onDomain :: root_dist_sl

    IF (.NOT. tile%Is_process_active(CH4_)) RETURN

    model => get_model(tile%owner_model_id)
    layer => Get_vgrid('soil_depth_water')
    grid  => get_grid(model%grid_id)
    nproma = grid%Get_nproma()
    nblks  = grid%Get_nblks()

    ! Get soil memory of the tile
    dsl4jsb_Get_memory(CH4_)
    dsl4jsb_Get_memory(HYDRO_)

    ! Get soil config
    dsl4jsb_Get_config(CH4_)

    IF (debug_on() .OR. dsl4jsb_Config(CH4_)%dbg_process) &
      CALL message(TRIM(routine), 'Initializing methane memory for tile '//TRIM(tile%name))

    ! Determine dz for soil layer midpoints
    dsl4jsb_memory(CH4_)%layer_midpoint_dz(1) = layer%dz(1) * 0.5_lp
    DO i=2,layer%n_levels
      dsl4jsb_memory(CH4_)%layer_midpoint_dz(i) = (layer%dz(i) + layer%dz(i-1)) * 0.5_lp
    ENDDO
    !$acc update device(dsl4jsb_memory(CH4_)%layer_midpoint_dz)

    ! Determine (constant) layer for in which bedrock is located
    dsl4jsb_Get_var2D_onDomain(HYDRO_, soil_depth)
    dsl4jsb_Get_var2D_onDomain(CH4_, mxlay_above_rock)
    dsl4jsb_Get_var3D_onDomain(CH4_, root_dist_sl)
    DO j=1,nblks
      DO i=1,nproma
        mxlay_above_rock(i,j) = REAL(find_bedrock_layer(layer%n_levels,layer%ubounds(:),soil_depth(i,j)),lp)
      ENDDO
    ENDDO
    !$acc update device(mxlay_above_rock)

    ! Determine distribution of roots in soil layers.
    ! NOTE: In case of dynamical roots, this should be moved into the time loop in update_methane
    !       and be changed from working on domain to chunk
    IF (ASSOCIATED(tile%parent_tile)) THEN
      IF (TRIM(tile%parent_tile%name(:3)) == 'veg') THEN
        DO j=1,nblks
          DO i=1,nproma
             CALL determine_root_dist(INT(mxlay_above_rock(i,j)+0.5_lp),layer%dz(:),layer%ubounds(:), &
               &                      dsl4jsb_Lctlib_param(root_dist_beta),                           &
               &                      root_dist_sl(i,:,j))
          ENDDO
        ENDDO
      ELSE
        root_dist_sl(:,:,:) = 0._lp
      ENDIF
    END IF
    !$acc update device(root_dist_sl)

  END SUBROUTINE methane_init_ic

#endif

END MODULE mo_methane_init
