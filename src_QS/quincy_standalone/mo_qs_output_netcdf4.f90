!> QUINCY netcdf4 output handling
!>
!> ICON-Land
!>
!> ---------------------------------------
!> Copyright (C) 2013-2026, MPI-M, MPI-BGC
!>
!> Contact: icon-model.org
!> Authors: AUTHORS.md
!> See LICENSES/ for license information
!> SPDX-License-Identifier: BSD-3-Clause
!> ---------------------------------------
!>
!> For more information on the QUINCY model see: <https://doi.org/10.17871/quincy-model-2019>
!>
!>#### Provides routines to specify, create and write quincy netcdf4 output
!>
!> Usage & Application of the netcdf4 output as provided by this module
!>
!> - To modify or add output files and/or variables you should only need to touch the routine `define_nc4_output_files`
!> - In that routine, new files can be added via the function call:
!>        `CALL add_output_file(p_output, NAME_PREFIX_FOR_FILE, p_current_file)`
!>   where 
!>   p_output and p_current_file are standardised pointers required in each call and
!>   NAME_PREFIX_FOR_FILE is a prefix for the file name 
!> - "prefix" since it gets the output period and if applicable 'spinup' appended.
!>
!> - After creating the new output file, different interface routines can be called to append variables to the file, 
!>   many of them having defaults for how variables are aggregated, which units variables have, etc.:
!>   - add_pool_variable                     add normal pool variable that has no layers
!>   - add_sl_pool_variable                  add pool variable that has soil layers
!>   - add_cl_pool_variable                  add pool variable that has canopy layers
!>   - add_flux_variable                     add flux from normal variable
!>   - add_sl_flux_variable                  add flux from variable that has soil layers
!>   - add_flux_sum_for_veg_bgcm_element     add flux calculated as sum of all veg bgcm compartments for a given element
!>   - add_element_of_sb_bgcm_compartment    add a given element for a given compartment of a sb bgcm
!>   - add_element_of_veg_bgcm_compartment   add a given element for a given compartment of a veg bgcm
!>   - add_all_elements_of_1L_veg_bgcm       add all elements for a given one-layered veg bgcm (i.e. veg bgcm without compartments)
!>   - add_all_elements_of_1L_sb_bgcm        add all elements for a given one-layered sb bgcm (i.e. sb bgcm without compartments)
!>   - add_all_parts_of_veg_bgcm_for_elem    add all compartments for a given element of a veg bgcm
!>   - add_all_parts_of_sb_bgcm_for_elem     add all compartments for a given element of an sb bgcm (reduced set for simple_1d)
!>   - ... each of these routines takes a couple of input arguments, whereby all require
!>     the standardised pointer p_current_file, a NAME_PREFIX for the name of the output variable, a unit 
!>     and a pointer to the variable in the quincy memory,
!>     depending on the output you want to add this could be a 2D or 3D variable, a veg or sb pool, or other
!>
MODULE mo_qs_output_netcdf4

  USE mo_kind,                ONLY: wp
  USE mo_exception,           ONLY: message_text, finish, message
  USE mo_util,                ONLY: int2string, real2string
  USE mo_jsb_model_class,     ONLY: t_jsb_model
  USE mo_jsb_config_class,    ONLY: t_jsb_model_config
  USE mo_jsb_class,           ONLY: Get_model
  USE mo_quincy_model_config, ONLY: Get_quincy_model_config_name
  USE mo_jsb_var_class,       ONLY: t_jsb_var, t_jsb_var_p, t_jsb_var_real2d, t_jsb_var_real3d
  USE mo_jsb_varlist,         ONLY: VARNAME_LEN
  USE mo_jsb_tile_class,      ONLY: t_jsb_tile_abstract
  USE mo_jsb_task_class,      ONLY: t_jsb_task_options

  USE mo_jsb_grid_class,    ONLY: t_jsb_vgrid
  USE mo_jsb_grid,          ONLY: Get_vgrid

  USE mo_lnd_bgcm_idx
  USE mo_lnd_bgcm_class,    ONLY: get_short_name_for_element_id, get_element_index, &
    &                             ELEM_C_ID, ELEM_N_ID, ELEM_P_ID, ELEM_C13_ID, ELEM_C14_ID, ELEM_N15_ID

  USE mo_jsb_process_class,  ONLY: A2L_, SPQ_, VEG_, RAD_, Q_ASSIMI_, Q_PHENO_, SB_, PHYD_, L2A_, Q_SYL_, Q_AGR_, HYDRO_, SSE_, &
    &                              SEB_
  USE mo_jsb_math_constants, ONLY: one_year, one_day
  USE mo_sb_memory_class,       ONLY: t_sb_bgcm_with_components, t_sb_bgcm_with_elements
  USE mo_veg_memory_class,      ONLY: t_veg_bgcm_with_components, t_veg_bgcm_with_elements

  USE mo_veg_config_class,      ONLY: get_number_of_veg_compartments
  USE mo_sb_config_class,       ONLY: get_number_of_sb_compartments

  USE mo_quincy_model_config,   ONLY: QLAND, QPLANT, QCANOPY, QSOIL
  USE mo_quincy_output_class,   ONLY: t_quincy_output, t_quincy_outputfile, t_quincy_out_var, &
    & SOIL_LAYERED, CANOPY_LAYERED,                                                           &
    & NO_FACT, USE_NDAT_ACC, USE_VEG_FLUX_FACT,                                               &
    & unit_flux, unit_pool, unit_soillayer_flux, unit_soillayer_pool, unit_veg_pool_flux,     &
    & unit_hydro_pool, unit_hydro_flux, unit_water_potential, unitless, unit_fraction,        &
    & unit_energy, unit_temperature, unit_speed, unit_size, unit_bulk_dens,                   &
    & unit_cn, unit_nc, unit_np, unit_pn, unit_cp, unit_assimi_rates, unit_rain,              &
    & unit_conductivity
  ! USE mo_phy_schemes,           ONLY: qsat_mixed

#ifdef __QUINCY_WITH_NETCDF__
  USE netcdf
#endif

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: init_quincy_nc4_output, define_nc4_output_files,                              &
    &       prepare_nc4_output_writing, close_nc4_output_files, deallocate_quincy_output, &
    &       write_netcdf4_quincy_buffered, accumulate_on_nc_acc_vars,                     &
    &       add_to_nc_buffer, reset_nc_acc_vars, assert_possible_output_interval

  INTEGER, PRIVATE :: ipoi = 1  !> the point in ics:ice of which data will be written to file; per default = 1 
  ! ... which is how its done in mo_qs_output_txt...

  CHARACTER(len=*), PARAMETER, PRIVATE :: name_var_time             = 'time'
  CHARACTER(len=*), PARAMETER, PRIVATE :: this_calendar_str         = '365_day'
  CHARACTER(len=*), PARAMETER, PRIVATE :: modname                   = 'mo_qs_output_netcdf4'
  CHARACTER(len=*), PARAMETER, PRIVATE :: no_prefix                 = ''

CONTAINS

  ! ======================================================================================================= !
  !>
  !> Initialise quincy output on the model with model_id 
  !>
  SUBROUTINE init_quincy_nc4_output(model_id)
    INTEGER, INTENT(IN) :: model_id !< ID of model instance
    ! ----------------------------------------------------------------------------------------------------- !
    TYPE(t_jsb_model),     POINTER   :: model
    TYPE(t_quincy_output), POINTER   :: quincy_output
    CHARACTER(len=*),      PARAMETER :: routine = TRIM(modname)//':init_quincy_nc4_output'
    ! ----------------------------------------------------------------------------------------------------- !
    model => Get_model(model_id)

    ALLOCATE(quincy_output)

    model%quincy_output => quincy_output

    CALL read_output_variable_list(model_id)

  END SUBROUTINE init_quincy_nc4_output

  ! ======================================================================================================= !
  !>
  !> Initialise quincy output on the model with model_id
  !>
  SUBROUTINE read_output_variable_list(model_id)
    INTEGER, INTENT(IN) :: model_id !< ID of model instance
    ! ----------------------------------------------------------------------------------------------------- !
    TYPE(t_jsb_model),     POINTER   :: model
    TYPE(t_quincy_output), POINTER   :: quincy_output
    CHARACTER(len=*),      PARAMETER :: routine = TRIM(modname)//':read_output_variable_list'
    ! ----------------------------------------------------------------------------------------------------- !
    CHARACTER(len=256)               :: filepath
    INTEGER                          :: nlines = 0
    INTEGER                          :: fid = 1
    INTEGER                          :: i = 0, ierr = 0
    CHARACTER(len=256)               :: dummy_line
    LOGICAL                          :: flag_file_exists

    model => Get_model(model_id)
    filepath = model%config%file_sel_output_variables

    PRINT *, ' '
    PRINT *, 'Netcdf output information: '

    IF ((TRIM(filepath) == "not_set") .OR. (TRIM(filepath) == "")) THEN
      PRINT *, 'Note: No selected output variable file specified.'// &
        & ' Using all available output variables.'
      RETURN
    ENDIF

    ! Check if file exists 
    INQUIRE(file=TRIM(filepath), exist=flag_file_exists)
    IF (flag_file_exists) THEN 
    ! Open file 
      OPEN(unit=fid, file=TRIM(filepath), action='read', iostat=ierr)

      ! Obtain number of lines in the input file.
      DO WHILE (ierr .EQ. 0)
        nlines = nlines + 1
        READ(fid, *, iostat=ierr) dummy_line
      END DO
      nlines = nlines - 1
      CLOSE(fid)

      ! Allocate string array according to the number of variables
      ALLOCATE(model%quincy_output%sel_output_variables(nlines))

      ! Populate array with file content
      OPEN(unit=fid, file=TRIM(filepath), action='read')
      REWIND(fid)
      DO i = 1, nlines
          READ(fid,'(A)') model%quincy_output%sel_output_variables(i)
      END DO
      CLOSE(fid)
      PRINT *, 'Found selected output variable file: '//TRIM(filepath)// &
        & ' containing '//int2string(nlines)// ' output variables.'
    ELSE
      PRINT *, 'Warning: Could not find the selected output variable file: '//TRIM(filepath)// &
        & '. Using all available output variables.'
    ENDIF     
  END SUBROUTINE read_output_variable_list

  ! ======================================================================================================= !
  !>
  !> Define and create output files -> This is currently THE place for the (and a lot of) manual work!
  !>
  SUBROUTINE define_nc4_output_files(model_id, tile)
    dsl4jsb_Use_memory(L2A_)
    dsl4jsb_Use_memory(A2L_)
    dsl4jsb_Use_memory(Q_ASSIMI_)
    dsl4jsb_Use_memory(Q_PHENO_)
    dsl4jsb_Use_memory(RAD_)
    dsl4jsb_Use_memory(SB_)
    dsl4jsb_Use_memory(SPQ_)
    dsl4jsb_Use_memory(HYDRO_)
    dsl4jsb_Use_memory(SSE_)
    dsl4jsb_Use_memory(SEB_) 
    dsl4jsb_Use_memory(VEG_)
    dsl4jsb_Use_memory(PHYD_)
    dsl4jsb_Use_memory(Q_SYL_)
    dsl4jsb_Use_config(SB_)
    dsl4jsb_Use_config(PHYD_)
    dsl4jsb_Use_config(Q_SYL_)
    dsl4jsb_Use_config(Q_AGR_)
    dsl4jsb_Use_config(VEG_)

    ! ----------------------------------------------------------------------------------------------------- !
    INTEGER,                    INTENT(IN) :: model_id   !< ID of model instance
    CLASS(t_jsb_tile_abstract), INTENT(IN) :: tile       !< the tile we are working on...
    ! ----------------------------------------------------------------------------------------------------- !
    TYPE(t_quincy_outputfile), POINTER  :: p_current_file  !< a pointer to a specific file 
    TYPE(t_quincy_output),     POINTER  :: p_output        !< a pointer to the quincy output object
    TYPE(t_jsb_model),         POINTER  :: model           !< the model
    TYPE(t_jsb_vgrid),         POINTER  :: vgrid_soil_w    !< Vertical grid
    INTEGER                             :: nsoil_w         !< number of soil layers in the vgrid
    INTEGER                             :: nr_of_veg_compartments !< number of veg compartments
    INTEGER                             :: nr_of_sb_compartments  !< number of sb compartments
    LOGICAL                             :: l_acc           !< if a variable is to be accumulated or not
    LOGICAL                             :: flag_mycorrhiza !< if using mycorrhiza (from sb config)
    CHARACTER(len=:), ALLOCATABLE       :: sb_model_scheme !< used sb model scheme (from sb config)
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':define_nc4_output_files'
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Def_memory(L2A_)
    dsl4jsb_Def_memory(A2L_)
    dsl4jsb_Def_memory(Q_ASSIMI_)
    dsl4jsb_Def_memory(Q_PHENO_)
    dsl4jsb_Def_memory(RAD_)
    dsl4jsb_Def_memory(SB_)
    dsl4jsb_Def_memory(SPQ_)
    dsl4jsb_Def_memory(HYDRO_)
    dsl4jsb_Def_memory(SSE_)
    dsl4jsb_Def_memory(SEB_) 
    dsl4jsb_Def_memory(VEG_)
    dsl4jsb_Def_memory(PHYD_)
    dsl4jsb_Def_memory(Q_SYL_)
    dsl4jsb_Def_config(SB_)
    dsl4jsb_Def_config(PHYD_)
    dsl4jsb_Def_config(Q_SYL_)
    dsl4jsb_Def_config(Q_AGR_)
    dsl4jsb_Def_config(VEG_)
    ! ----------------------------------------------------------------------------------------------------- !
    model         => Get_model(model_id)
    vgrid_soil_w  => Get_vgrid('soil_depth_water')
    nsoil_w       =  vgrid_soil_w%n_levels
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Get_memory(L2A_)
    dsl4jsb_Get_memory(A2L_)
    dsl4jsb_Get_memory(Q_ASSIMI_)
    dsl4jsb_Get_memory(Q_PHENO_)
    dsl4jsb_Get_memory(RAD_)
    dsl4jsb_Get_memory(SB_)
    dsl4jsb_Get_memory(SPQ_)
    dsl4jsb_Get_memory(HYDRO_)
    dsl4jsb_Get_memory(SSE_)
    dsl4jsb_Get_memory(SEB_)
    dsl4jsb_Get_memory(VEG_)
    dsl4jsb_Get_memory(PHYD_)
    dsl4jsb_Get_memory(Q_SYL_)
    dsl4jsb_Get_config(SB_)
    dsl4jsb_Get_config(PHYD_)
    dsl4jsb_Get_config(Q_SYL_)
    dsl4jsb_Get_config(Q_AGR_)
    dsl4jsb_Get_config(VEG_)

    !> 0.9 test availability of 'quincy_output' object
    !>
    IF (.NOT. ASSOCIATED(model%quincy_output)) THEN
      WRITE(message_text,*) 'Quincy output not initialised, please check!'
      CALL finish(TRIM(routine),message_text)
    ENDIF
    p_output => model%quincy_output

    ! ----------------------------------------------------------------------------------------------------- !
    ! NOTES...
    ! Note: 'no pool' flux vars were previously converted; now default: only averaged over time (USE_NDAT_ACC)
    ! Note: now molC molP molN -> mol (CMOR!)
    ! Note: Growth_X and Litterfall_X do not include the seed_bed here, instead these are passed as separate variables
    !      (i.e. different than in mo_qs_output_txt.f90)
    ! ----------------------------------------------------------------------------------------------------- !
    !>
    !> 1. assimilation output file (ASSIMI)
    !>
    SELECT CASE (model%config%qmodel_id)
    CASE (QLAND, QPLANT, QSOIL, QCANOPY)
      CALL add_output_file(p_output, 'Q_ASSIMI', p_current_file)
      CALL add_pool_variable(p_current_file, 'rootzone_soilwater_potential', dsl4jsb_var(Q_ASSIMI_,wtr_soil_root_pot), &
        &  unit_water_potential)
      SELECT CASE (model%config%qmodel_id)
      CASE (QLAND, QPLANT, QCANOPY)
        CALL add_cl_pool_variable(p_current_file, 'jmax', dsl4jsb_var(Q_ASSIMI_,jmax_cl), unit=unit_assimi_rates)
        CALL add_cl_pool_variable(p_current_file, 'vcmax', dsl4jsb_var(Q_ASSIMI_,vcmax_cl), unit=unit_assimi_rates)
        CALL add_cl_pool_variable(p_current_file, 'chlfl_yield', dsl4jsb_var(Q_ASSIMI_,chlfl_yield_cl), unitless)
        CALL add_flux_variable(p_current_file, 'gpp', dsl4jsb_var(Q_ASSIMI_,gross_assimilation))
        CALL add_flux_variable(p_current_file, 'gpp_c13', dsl4jsb_var(Q_ASSIMI_,gross_assimilation_C13))
        CALL add_flux_variable(p_current_file, 'gpp_c14', dsl4jsb_var(Q_ASSIMI_,gross_assimilation_C14))
        CALL add_flux_variable(p_current_file, 'gc', dsl4jsb_var(Q_ASSIMI_,canopy_cond), unit=unit_speed)
        CALL add_flux_variable(p_current_file, 'ga', dsl4jsb_var(Q_ASSIMI_,aerodyn_cond), unit=unit_speed)
        SELECT CASE (model%config%qmodel_id)
        CASE (QLAND, QPLANT)
          CALL add_pool_variable(p_current_file, 'beta_soil_gs', dsl4jsb_var(Q_ASSIMI_,beta_soil_gs), unitless)
          CALL add_pool_variable(p_current_file, 'beta_soil_ps', dsl4jsb_var(Q_ASSIMI_,beta_soil_ps), unitless)
        END SELECT ! ** if land or plant model
      END SELECT ! ** if land, plant or canopy model
    END SELECT ! ** if land, plant, soil or canopy model

    ! ----------------------------------------------------------------------------------------------------- !
    !>
    !> 2. atmland output file (a2l)
    !>
    SELECT CASE (model%config%qmodel_id)
    CASE (QLAND, QPLANT, QCANOPY)
      CALL add_output_file(p_output, 'A2L', p_current_file)
      CALL add_flux_variable(p_current_file, 'rain', dsl4jsb_var(A2L_,rain), unit=unit_rain)
      CALL add_flux_variable(p_current_file, 'snow', dsl4jsb_var(A2L_,snow), unit=unit_rain)
      CALL add_flux_variable(p_current_file, 'air_temperature', dsl4jsb_var(A2L_,t_air), unit=unit_temperature)
      CALL add_flux_variable(p_current_file, 'wind', dsl4jsb_var(A2L_,wind_air), unit=unit_speed)
      CALL add_flux_variable(p_current_file, 'swnir_srf_down', dsl4jsb_var(A2L_,swnir_srf_down), unit=unit_energy)
      CALL add_flux_variable(p_current_file, 'swvis_srf_down', dsl4jsb_var(A2L_,swvis_srf_down), unit=unit_energy)
      CALL add_flux_variable(p_current_file, 'swpar_srf_down', dsl4jsb_var(A2L_,swpar_srf_down), unit=unit_energy)
      CALL add_flux_variable(p_current_file, 'lw_down', dsl4jsb_var(A2L_,lw_srf_down), unit=unit_energy)
      ! @TODO: how to  vpd_acc = vpd_acc + (( qsat_mixed(t_air(ipoi),press_srf(ipoi)) - q_air(ipoi) ) / &
      !                         eps_vpd * press_srf(ipoi) / 1000.0_wp)
      ! @TODO SZ: make new variable, probalby best in ATMLand? -> #92
      CALL add_flux_variable(p_current_file, 'atm_co2', dsl4jsb_var(A2L_,co2_mixing_ratio), unit=unit_fraction)
      CALL add_flux_variable(p_current_file, 'atm_co2_c13', dsl4jsb_var(A2L_,co2_mixing_ratio_C13), unit=unit_fraction)
      ! @TODO: check unit
      CALL add_flux_variable(p_current_file, 'atm_co2_c14', dsl4jsb_var(A2L_,co2_mixing_ratio_C14), unit='scaled?')
      CALL add_flux_variable(p_current_file, 'nhx_deposition', dsl4jsb_var(A2L_,nhx_deposition))
      CALL add_flux_variable(p_current_file, 'noy_deposition', dsl4jsb_var(A2L_,noy_deposition))
      CALL add_flux_variable(p_current_file, 'nhx_n15_deposition', dsl4jsb_var(A2L_,nhx_n15_deposition))
      CALL add_flux_variable(p_current_file, 'noy_n15_deposition', dsl4jsb_var(A2L_,noy_n15_deposition))
      CALL add_flux_variable(p_current_file, 'p_deposition', dsl4jsb_var(A2L_,p_deposition))
    END SELECT ! ** if land, plant or canopy model

    ! ----------------------------------------------------------------------------------------------------- !
    !>
    !> 3. phenology output file (q_pheno)
    !>
    SELECT CASE (model%config%qmodel_id)
    CASE (QLAND, QPLANT)
      CALL add_output_file(p_output, 'Q_PHENO', p_current_file)
      CALL add_pool_variable(p_current_file, 'growing_season', dsl4jsb_var(Q_PHENO_,growing_season), unitless)
    END SELECT ! ** if land or plant model
    
    ! ----------------------------------------------------------------------------------------------------- !
    !>
    !> 4. radiation output file (rad)
    !>
    SELECT CASE (model%config%qmodel_id)
    CASE (QLAND, QPLANT, QCANOPY)
      CALL add_output_file(p_output, 'RAD', p_current_file)
      CALL add_flux_variable(p_current_file, 'r_net', dsl4jsb_var(RAD_,rad_srf_net), unit=unit_energy)
      ! Noon variables for longer term integrations (these make sensible daily + timescale means
      CALL add_pool_variable(p_current_file, 'albedo_noon', dsl4jsb_var(RAD_,albedo_noon), unit=unitless)
      CALL add_flux_variable(p_current_file, 'fapar_noon', dsl4jsb_var(RAD_,fapar_noon))
      ! Instantaneous albedo (probably makes sense only with time-step output because of the averaging problem
      CALL add_pool_variable(p_current_file, 'albedo_vis_lnd', dsl4jsb_var(RAD_,alb_vis_lnd), unitless)
      CALL add_pool_variable(p_current_file, 'albedo_nir_lnd', dsl4jsb_var(RAD_,alb_nir_lnd), unitless)
    END SELECT ! ** if land, plant or canopy model

    ! ----------------------------------------------------------------------------------------------------- !
    !>
    !> 5. soil biogeochemistry output file (sb)
    !>
    IF(model%config%flag_write_sb_output) THEN
      nr_of_sb_compartments = get_number_of_sb_compartments()
      sb_model_scheme = TRIM(dsl4jsb_Config(SB_)%sb_model_scheme)
      flag_mycorrhiza = dsl4jsb_Config(SB_)%flag_mycorrhiza

      SELECT CASE (model%config%qmodel_id)
      CASE (QLAND, QSOIL)
        CALL add_output_file(p_output, 'SB', p_current_file)
        !@TODO SZ: deltas in postprocessing -> #93

        ! ... total c in sb pool
        CALL add_sl_pool_variable(p_current_file, 'sb_total_c', dsl4jsb_var(SB_,sb_pool_total_c))
        CALL add_sl_pool_variable(p_current_file, 'sb_total_som_c', dsl4jsb_var(SB_,sb_pool_total_bg_soil_c))
        CALL add_pool_variable(p_current_file, 'sb_total_ag_litter_c', dsl4jsb_var(SB_,sb_pool_total_ag_litter_c))

        !--- previously soilpoolOrg and as 1D: soilpooldC13DC14dN15 / soilpoolC13C14N15
        CALL add_all_parts_of_sb_bgcm_for_elem(p_current_file, no_prefix, ELEM_C_ID, &
          & dsl4jsb_bgcm(SB_,sbbpool), model%config, nsoil_w, nr_of_sb_compartments, sb_model_scheme)
        CALL add_all_parts_of_sb_bgcm_for_elem(p_current_file, no_prefix, ELEM_N_ID, &
          & dsl4jsb_bgcm(SB_,sbbpool), model%config, nsoil_w, nr_of_sb_compartments, sb_model_scheme)
        CALL add_all_parts_of_sb_bgcm_for_elem(p_current_file, no_prefix, ELEM_P_ID, &
          & dsl4jsb_bgcm(SB_,sbbpool), model%config, nsoil_w, nr_of_sb_compartments, sb_model_scheme)
        CALL add_all_parts_of_sb_bgcm_for_elem(p_current_file, no_prefix, ELEM_C13_ID, &
          & dsl4jsb_bgcm(SB_,sbbpool), model%config, nsoil_w, nr_of_sb_compartments, sb_model_scheme)
        CALL add_all_parts_of_sb_bgcm_for_elem(p_current_file, no_prefix, ELEM_C14_ID, &
          & dsl4jsb_bgcm(SB_,sbbpool), model%config, nsoil_w, nr_of_sb_compartments, sb_model_scheme)
        CALL add_all_parts_of_sb_bgcm_for_elem(p_current_file, no_prefix, ELEM_N15_ID, &
          & dsl4jsb_bgcm(SB_,sbbpool), model%config, nsoil_w, nr_of_sb_compartments, sb_model_scheme)

        !--- previously soilpoolInOrg_L and as 1D: soilpoolInOrg
        CALL add_sl_pool_variable(p_current_file, 'mineral_assoc_nh4', dsl4jsb_var(SB_,nh4_assoc))
        CALL add_sl_pool_variable(p_current_file, 'mineral_nh4', dsl4jsb_var(SB_,nh4_solute))
        CALL add_sl_pool_variable(p_current_file, 'mineral_no3', dsl4jsb_var(SB_,no3_solute))
        CALL add_sl_pool_variable(p_current_file, 'mineral_p', dsl4jsb_var(SB_,po4_solute))
        CALL add_sl_pool_variable(p_current_file, 'noy', dsl4jsb_var(SB_,noy))
        CALL add_sl_pool_variable(p_current_file, 'n2o', dsl4jsb_var(SB_,n2o))
        CALL add_sl_pool_variable(p_current_file, 'n2', dsl4jsb_var(SB_,n2))
        CALL add_sl_pool_variable(p_current_file, 'n15_nh4_assoc', dsl4jsb_var(SB_,nh4_n15_assoc))
        CALL add_sl_pool_variable(p_current_file, 'n15_nh4_solute', dsl4jsb_var(SB_,nh4_n15_solute))
        CALL add_sl_pool_variable(p_current_file, 'n15_no3_solute', dsl4jsb_var(SB_,no3_n15_solute))
        CALL add_sl_pool_variable(p_current_file, 'po4_assoc_fast', dsl4jsb_var(SB_,po4_assoc_fast))
        CALL add_sl_pool_variable(p_current_file, 'po4_assoc_slow', dsl4jsb_var(SB_,po4_assoc_slow))
        CALL add_sl_pool_variable(p_current_file, 'po4_assoc_occl', dsl4jsb_var(SB_,po4_occluded))
        CALL add_sl_pool_variable(p_current_file, 'po4_solute', dsl4jsb_var(SB_,po4_solute))
        CALL add_sl_pool_variable(p_current_file, 'po4_primary', dsl4jsb_var(SB_,po4_primary))
        !--- previously soil_diagnostics_microbial

        IF(sb_model_scheme == "jsm") THEN
          ! NOTE: vars are only added as sl vars and not as averaged 1d vars (as in mo_qs_output_txt)
          CALL add_sl_pool_variable(p_current_file, 'cue_mic', &
            & dsl4jsb_var(SB_,microbial_cue_eff_tmic_mavg), l_write_1d=.FALSE., unit_sl=unitless)
          CALL add_sl_pool_variable(p_current_file, 'nue_mic', &
            & dsl4jsb_var(SB_,microbial_nue_eff_tmic_mavg), l_write_1d=.FALSE., unit_sl=unitless)
          CALL add_sl_pool_variable(p_current_file, 'pue_mic', &
            & dsl4jsb_var(SB_,microbial_pue_eff_tmic_mavg), l_write_1d=.FALSE., unit_sl=unitless)
          !@TODO -howto: veg_pool%fine_root%carbon(ipoi) * root_fraction_sl(ipoi,:) / soil_depth_sl(ipoi,:)
          CALL add_sl_pool_variable(p_current_file, 'enz_poly', &
            & dsl4jsb_var(SB_,enzyme_frac_poly), l_write_1d=.FALSE., unit_sl=unitless)
          CALL add_sl_pool_variable(p_current_file, 'enz_res', &
            & dsl4jsb_var(SB_,enzyme_frac_residue), l_write_1d=.FALSE., unit_sl=unitless)
          CALL add_sl_pool_variable(p_current_file, 'enz_poly_c', &
            & dsl4jsb_var(SB_,enzyme_frac_poly_c), l_write_1d=.FALSE., unit_sl=unitless)
          CALL add_sl_pool_variable(p_current_file, 'enz_poly_n', &
            & dsl4jsb_var(SB_,enzyme_frac_poly_n), l_write_1d=.FALSE., unit_sl=unitless)
          CALL add_sl_pool_variable(p_current_file, 'enz_poly_p', &
            & dsl4jsb_var(SB_,enzyme_frac_poly_p), l_write_1d=.FALSE., unit_sl=unitless)
          CALL add_sl_pool_variable(p_current_file, 'enz_poly_c_mavg', &
            & dsl4jsb_var(SB_,enzyme_frac_poly_c_mavg), l_write_1d=.FALSE., unit_sl=unitless)
          CALL add_sl_pool_variable(p_current_file, 'enz_poly_n_mavg', &
            & dsl4jsb_var(SB_,enzyme_frac_poly_n_mavg), l_write_1d=.FALSE., unit_sl=unitless)
          CALL add_sl_pool_variable(p_current_file, 'enz_poly_p_mavg', &
            & dsl4jsb_var(SB_,enzyme_frac_poly_p_mavg), l_write_1d=.FALSE., unit_sl=unitless)
          !@TODO: dom_cn_mavg and dom_cp_mavg are unitless in sb memory -> #100 
          !@TODO SZ: unit mol C mol-1 N
          CALL add_sl_pool_variable(p_current_file, 'dom_cn_mavg', dsl4jsb_var(SB_,dom_cn_mavg), l_write_1d=.FALSE., unit_sl=unit_cn)
          !@TODO SZ: unit mol N mol-1 P --> Julia: C!?
          CALL add_sl_pool_variable(p_current_file, 'dom_cp_mavg', dsl4jsb_var(SB_,dom_cp_mavg), l_write_1d=.FALSE., unit_sl=unit_cp)
        ENDIF ! sb_model_scheme == "jsm"

        !--- previously: soil_diagnostics_sorption_lay
        ! NOTE: vars are only added as sl vars and not as averaged 1d vars (as in mo_qs_output_txt)
        CALL add_sl_pool_variable(p_current_file, 'qmax_org', &
          & dsl4jsb_var(SB_,qmax_org), l_write_1d=.FALSE., unit_sl='kg C kg-1 mineral soil')
        CALL add_sl_pool_variable(p_current_file, 'rho_cor', & 
          & dsl4jsb_var(SB_,bulk_dens_corr_sl), l_write_1d=.FALSE., unit_sl=unit_bulk_dens)
        !@TODO SZ: question for Lin Yu, should be kg P kg-1 min soil, but I'm not sure 
        ! -- in sb memory several of these have "mol P m-3 soil" -> now mol/m3 --> is this correct? or not?
        CALL add_sl_pool_variable(p_current_file, 'qmax_po4', dsl4jsb_var(SB_,qmax_po4), l_write_1d=.FALSE.)
        CALL add_sl_pool_variable(p_current_file, 'mm_po4', dsl4jsb_var(SB_,km_adsorpt_po4_sl), l_write_1d=.FALSE.)
        !--- previously SPQ_ memory
        CALL add_sl_pool_variable(p_current_file, 'min_vol', dsl4jsb_var(SB_,volume_min_sl), l_write_1d=.FALSE., unit_sl=unitless)

        !*** fluxes
        !--- previously: soilfluxCdC13DC14 + soilfluxC_lay
        CALL add_sl_flux_variable(p_current_file, 'het_respiration', dsl4jsb_var(SB_,het_respiration), nsoil_w)
        CALL add_sl_flux_variable(p_current_file, 'het_respiration_c13', dsl4jsb_var(SB_,het_respiration_c13), nsoil_w)
        CALL add_sl_flux_variable(p_current_file, 'het_respiration_c14', dsl4jsb_var(SB_,het_respiration_c14), nsoil_w)

        IF(sb_model_scheme == "jsm") THEN
          CALL add_flux_variable(p_current_file, 'leaching_dom_C', dsl4jsb_var(SB_,leaching_dom_carbon))
          CALL add_element_of_sb_bgcm_compartment(p_current_file, 'transport_dom_C', &
            & dsl4jsb_bgcm(SB_,sbbtransport), ix_dom, ELEM_C_ID, model%config, nsoil_w)
          CALL add_sl_flux_variable(p_current_file, 'lateral_loss_dom_C', dsl4jsb_var(SB_,lateral_loss_dom_carbon_sl), nsoil_w)
          CALL add_sl_flux_variable(p_current_file, 'lateral_loss_dom_13C', dsl4jsb_var(SB_,lateral_loss_dom_carbon13_sl), nsoil_w)
          CALL add_sl_flux_variable(p_current_file, 'lateral_loss_dom_14C', dsl4jsb_var(SB_,lateral_loss_dom_carbon14_sl), nsoil_w)
        ENDIF ! sb_model_scheme == "jsm"

        !@TODO -- howto       root_respiration_lay_acc(1:nsoil_w)
        !@TODO SZ: needs a new variable in SB -> #92
        !  root_respiration_lay_acc(1:nsoil_w)    = root_respiration_lay_acc(1:nsoil_w) + &
        !                                            fmaint_rate_root(ipoi) * &
        !                                            (veg_pool%fine_root%nitrogen(ipoi) + veg_pool%coarse_root%nitrogen(ipoi)) * &
        !                                            root_fraction_sl(ipoi,:) / soil_depth_sl(ipoi,:) * dtime_one

        !--- previously: soilfluxN  [mol N /m2 /interval] and soilfluxN_lay  [mol N /m3 /interval]
        CALL add_sl_flux_variable(p_current_file, 'n_fix_asymb', dsl4jsb_var(SB_,asymb_n_fixation), nsoil_w)
        CALL add_sl_flux_variable(p_current_file, 'net_min_nh4', dsl4jsb_var(SB_,net_mineralisation_nh4), nsoil_w)
        CALL add_sl_flux_variable(p_current_file, 'net_min_no3', dsl4jsb_var(SB_,net_mineralisation_no3), nsoil_w)
        ! NOTE: in mo_qs_output_txt.f90 nitrification_acc was written as a sum of no3 noy and n2o 
        !       -> SZ simplify to only write out nitrification_n2o
        CALL add_sl_flux_variable(p_current_file, 'nitrification_n2o', dsl4jsb_var(SB_,nitrification_n2o), nsoil_w)
        ! NOTE: in mo_qs_output_txt.f90 denitrification_acc was written as a sum of no3 noy and n2o 
        !       -> SZ simplify to only write out denitrification_n2o
        CALL add_sl_flux_variable(p_current_file, 'denitrification_n2o', dsl4jsb_var(SB_,denitrification_n2o), nsoil_w)
        CALL add_flux_variable(p_current_file, 'emission_noy', dsl4jsb_var(SB_,emission_noy))
        CALL add_flux_variable(p_current_file, 'emission_n2o', dsl4jsb_var(SB_,emission_n2o))
        CALL add_flux_variable(p_current_file, 'emission_n2', dsl4jsb_var(SB_,emission_n2))
        CALL add_sl_flux_variable(p_current_file, 'plant_uptake_nh4', dsl4jsb_var(SB_,plant_uptake_nh4_sl), nsoil_w)
        ! MicUptake_nh4 and MicUptake_no3: in mo_qs_output_txt: only in soilfluxN_lay (i.e. only layered var not 1D var)
        CALL add_sl_flux_variable(p_current_file, 'mic_uptake_nh4', & 
          &  dsl4jsb_var(SB_,microbial_uptake_nh4_sl), nsoil_w, l_write_1d=.FALSE.)
        CALL add_sl_flux_variable(p_current_file, 'mic_uptake_no3', & 
          & dsl4jsb_var(SB_,microbial_uptake_no3_sl), nsoil_w, l_write_1d=.FALSE.)
        CALL add_sl_flux_variable(p_current_file, 'plant_uptake_no3', dsl4jsb_var(SB_,plant_uptake_no3_sl), nsoil_w)
        CALL add_flux_variable(p_current_file, 'leaching_nh4', dsl4jsb_var(SB_,leaching_nh4_solute))
        CALL add_flux_variable(p_current_file, 'leaching_no3', dsl4jsb_var(SB_,leaching_no3_solute))
        IF(sb_model_scheme == "jsm") THEN
          CALL add_flux_variable(p_current_file, 'leaching_dom_n', dsl4jsb_var(SB_,leaching_dom_nitrogen))
        ENDIF
        CALL add_sl_flux_variable(p_current_file, 'fast_exchange_nh4', &
          & dsl4jsb_var(SB_,fast_exchange_nh4), nsoil_w, l_write_2d=.FALSE.) ! NOTE: only 1D
        IF(flag_mycorrhiza) THEN
          CALL add_sl_flux_variable(p_current_file, 'myc_uptake_nh4', dsl4jsb_var(SB_,mycorrhiza_uptake_nh4_sl), nsoil_w)
          CALL add_sl_flux_variable(p_current_file, 'myc_uptake_no3', dsl4jsb_var(SB_,mycorrhiza_uptake_no3_sl), nsoil_w)
          CALL add_sl_flux_variable(p_current_file, 'myc_uptake_n_org', dsl4jsb_var(SB_,mycorrhiza_uptake_norg_sl), nsoil_w)
        ENDIF
        ! Transport_nh4 and Transport_assoc_nh4: in mo_qs_output_txt: only in soilfluxN_lay (i.e. only layered var not 1D var)
        CALL add_sl_flux_variable(p_current_file, 'transport_nh4', &
          & dsl4jsb_var(SB_,transport_nh4_solute), nsoil_w, l_write_1d=.FALSE.)
        CALL add_sl_flux_variable(p_current_file, 'lateral_loss_nh4', dsl4jsb_var(SB_,lateral_loss_nh4_solute_sl), nsoil_w)
        CALL add_sl_flux_variable(p_current_file, 'transport_assoc_nh4', &
          & dsl4jsb_var(SB_,transport_nh4_assoc), nsoil_w, l_write_1d=.FALSE.)
        ! Transport_no3: in mo_qs_output_txt: only in soilfluxN_lay (i.e. only layered var not 1D var)
        CALL add_sl_flux_variable(p_current_file, 'transport_no3', &
          & dsl4jsb_var(SB_,transport_no3_solute), nsoil_w, l_write_1d=.FALSE.)
        CALL add_sl_flux_variable(p_current_file, 'lateral_loss_no3', dsl4jsb_var(SB_,lateral_loss_no3_solute_sl), nsoil_w)
        IF(sb_model_scheme == "jsm") THEN
          CALL add_element_of_sb_bgcm_compartment(p_current_file, 'transport_dom_N', &
            & dsl4jsb_bgcm(SB_,sbbtransport), ix_dom, ELEM_N_ID, model%config, nsoil_w)
          CALL add_sl_flux_variable(p_current_file, 'lateral_loss_dom_N', dsl4jsb_var(SB_,lateral_loss_dom_nitrogen_sl), nsoil_w)
          ! Adsorpt_nh4 : in mo_qs_output_txt: only in soilfluxN_lay (i.e. only layered var not 1D var)
        ENDIF
        CALL add_sl_flux_variable(p_current_file, 'adsorpt_nh4', dsl4jsb_var(SB_,fast_adsorpt_nh4), nsoil_w, l_write_1d=.FALSE.)
        !--- previously: soilfluxdN
        CALL add_sl_flux_variable(p_current_file, 'n15_fix_asymb', dsl4jsb_var(SB_,asymb_n_fixation_n15), &
          & nsoil_w, l_write_2d=.FALSE.) ! only write 1D
        CALL add_flux_variable(p_current_file, 'emission_noy_n15', dsl4jsb_var(SB_,emission_noy_n15))
        CALL add_flux_variable(p_current_file, 'emission_n2o_n15', dsl4jsb_var(SB_,emission_n2o_n15))
        CALL add_flux_variable(p_current_file, 'emission_n2_n15', dsl4jsb_var(SB_,emission_n2_n15))
        CALL add_flux_variable(p_current_file, 'leaching_nh4_n15', dsl4jsb_var(SB_,leaching_nh4_n15_solute))
        CALL add_flux_variable(p_current_file, 'leaching_no3_n15', dsl4jsb_var(SB_,leaching_no3_n15_solute))
        CALL add_sl_flux_variable(p_current_file, 'lateral_loss_nh4_n15', dsl4jsb_var(SB_,lateral_loss_nh4_n15_solute_sl),nsoil_w)
        CALL add_sl_flux_variable(p_current_file, 'lateral_loss_no3_n15', dsl4jsb_var(SB_,lateral_loss_no3_n15_solute_sl),nsoil_w)
        !--- previously: soilfluxP [mol P /m2 /interval] and soilfluxP_lay  [mol P /m3 /interval]
        CALL add_flux_variable(p_current_file, 'leaching_po4', dsl4jsb_var(SB_,leaching_po4_solute))
        IF(sb_model_scheme == "jsm") THEN
          CALL add_flux_variable(p_current_file, 'leaching_dom_p', dsl4jsb_var(SB_,leaching_dom_phosphorus))
        ENDIF
        CALL add_sl_flux_variable(p_current_file, 'net_mineralisation_po4', dsl4jsb_var(SB_,net_mineralisation_po4), nsoil_w)
        CALL add_sl_flux_variable(p_current_file, 'plant_uptake_po4', dsl4jsb_var(SB_,plant_uptake_po4_sl), nsoil_w)
        CALL add_sl_flux_variable(p_current_file, 'weathering_po4', dsl4jsb_var(SB_,weathering_po4), nsoil_w)
        CALL add_sl_flux_variable(p_current_file, 'occlusion_po4', dsl4jsb_var(SB_,occlusion_po4), nsoil_w)
        CALL add_sl_flux_variable(p_current_file, 'slow_exchange_po4', dsl4jsb_var(SB_,slow_exchange_po4), nsoil_w)
        CALL add_sl_flux_variable(p_current_file, 'fast_exchange_po4', dsl4jsb_var(SB_,fast_exchange_po4), nsoil_w)
        CALL add_sl_flux_variable(p_current_file, 'biochem_mineralisation_po4', &
          & dsl4jsb_var(SB_,biochem_mineralisation_po4), nsoil_w)
        IF(sb_model_scheme == "jsm") THEN
          CALL add_sl_flux_variable(p_current_file, 'lateral_loss_dom_p', dsl4jsb_var(SB_,lateral_loss_dom_phosphorus_sl),nsoil_w)
        ENDIF
        CALL add_sl_flux_variable(p_current_file, 'lateral_loss_po4', dsl4jsb_var(SB_,lateral_loss_po4_solute_sl),nsoil_w)
        !--- previously: sb_massBalance
        !@TODO: what is the unit here? (in mo_qs_output_txt these were accumulated but not averaged)-> mumol/m2?!
        !@TODO SZ: mol m-3
        !@TODO: but these are (weighted) summes over the whole soil column (and over the output timestep)
        !       -> really not micro mol/m2 ? -> discuss...
        CALL add_sl_flux_variable(p_current_file, 'sb_delta_pool_c', dsl4jsb_elem_of_1L_bgcm(SB_,sbbdelta_pools,carbon), &
          & nsoil_w, out_fact_kind=NO_FACT, l_write_2d=.FALSE., unit_1d='?') ! only write 1D
        CALL add_sl_flux_variable(p_current_file, 'sb_delta_pool_n', dsl4jsb_elem_of_1L_bgcm(SB_,sbbdelta_pools,nitrogen), &
          & nsoil_w, out_fact_kind=NO_FACT, l_write_2d=.FALSE., unit_1d='?') ! only write 1D
        CALL add_sl_flux_variable(p_current_file, 'sb_delta_pool_p', dsl4jsb_elem_of_1L_bgcm(SB_,sbbdelta_pools,phosphorus), &
          & nsoil_w, out_fact_kind=NO_FACT, l_write_2d=.FALSE., unit_1d='?') ! only write 1D
        CALL add_sl_flux_variable(p_current_file, 'sb_delta_pool_c13', dsl4jsb_elem_of_1L_bgcm(SB_,sbbdelta_pools,carbon13), &
          & nsoil_w, out_fact_kind=NO_FACT, l_write_2d=.FALSE., unit_1d='?') ! only write 1D
        CALL add_sl_flux_variable(p_current_file, 'sb_delta_pool_c14', dsl4jsb_elem_of_1L_bgcm(SB_,sbbdelta_pools,carbon14), &
          & nsoil_w, out_fact_kind=NO_FACT, l_write_2d=.FALSE., unit_1d='?') ! only write 1D
        CALL add_sl_flux_variable(p_current_file, 'sb_delta_pool_n15', dsl4jsb_elem_of_1L_bgcm(SB_,sbbdelta_pools,nitrogen15), &
          & nsoil_w, out_fact_kind=NO_FACT, l_write_2d=.FALSE., unit_1d='?') ! only write 1D

        CALL add_flux_variable(p_current_file, 'sb_total_flux_c', dsl4jsb_var(SB_,total_flux_carbon))
        CALL add_flux_variable(p_current_file, 'sb_total_flux_n', dsl4jsb_var(SB_,total_flux_nitrogen))
        CALL add_flux_variable(p_current_file, 'sb_total_flux_p', dsl4jsb_var(SB_,total_flux_phosphorus))
        CALL add_flux_variable(p_current_file, 'sb_total_flux_c13', dsl4jsb_var(SB_,total_flux_carbon13))
        CALL add_flux_variable(p_current_file, 'sb_total_flux_c14', dsl4jsb_var(SB_,total_flux_carbon14))
        CALL add_flux_variable(p_current_file, 'sb_total_flux_n15', dsl4jsb_var(SB_,total_flux_nitrogen15))

        IF(flag_mycorrhiza .AND. model%config%qmodel_id == QLAND) THEN
          !--- previously: mycfluxCNP
          CALL add_all_elements_of_1L_sb_bgcm(p_current_file, 'myc_export', dsl4jsb_bgcm(SB_,sbbmyco_export), nsoil_w)
        ENDIF
      END SELECT ! ** if land or soil model
    END IF ! IF model%config%flag_write_sb_output
    
    ! ----------------------------------------------------------------------------------------------------- !
    !> soil phys output files
    !>
    SELECT CASE (model%config%qmodel_id)
    CASE (QLAND, QPLANT, QSOIL, QCANOPY)
      !>
      !> 6. soil phys quincy output file (spq)
      !>
      IF (.NOT. model%config%use_soil_phys_jsbach) THEN
        CALL add_output_file(p_output, 'SPQ', p_current_file)
        SELECT CASE (model%config%qmodel_id)
        CASE (QLAND, QPLANT, QCANOPY)
          CALL add_pool_variable(p_current_file, 'snow_height', dsl4jsb_var(SPQ_,snow_height), unit_hydro_pool)
          CALL add_flux_variable(p_current_file, 'evapopot', dsl4jsb_var(SPQ_,evapopot), unit = unit_hydro_flux)
          CALL add_flux_variable(p_current_file, 'evaporation', dsl4jsb_var(SPQ_,evaporation), unit = unit_hydro_flux)
          CALL add_flux_variable(p_current_file, 'interception', dsl4jsb_var(SPQ_,interception), unit = unit_hydro_flux)
          CALL add_flux_variable(p_current_file, 'srf_runoff', dsl4jsb_var(SPQ_,srf_runoff), unit = unit_hydro_flux)
          CALL add_flux_variable(p_current_file, 'gw_runoff', dsl4jsb_var(SPQ_,gw_runoff), unit = unit_hydro_flux)
          CALL add_flux_variable(p_current_file, 'drainage', dsl4jsb_var(SPQ_,drainage), unit = unit_hydro_flux)
          !--- previously: vegfluxEnergy 
          ! mo_qs_output_txt: vegfluxEnergy [Wm-2] [Wm-2] [Wm-2] [Wm-2] [--] [µmol m-2 s-1] [µmol m-2 s-1] [--] [--] [--]
          ! NOTE: Qle and Qh had a negative factor USE_DTIME_ONE_NEG_FACT; now also default (USE_NDAT_ACC) i.e. not negated
          CALL add_flux_variable(p_current_file, 'qh', dsl4jsb_var(SPQ_,sensible_heat_flx), unit=unit_energy)
          CALL add_flux_variable(p_current_file, 'qle', dsl4jsb_var(SPQ_,latent_heat_flx), unit=unit_energy)
          ! NOTE: ground heat flux had negated negation in mo_qs_output_txt; now default (USE_NDAT_ACC)
          CALL add_flux_variable(p_current_file, 'qg', dsl4jsb_var(SPQ_,ground_heat_flx), unit=unit_energy) 
        END SELECT ! ** if land, plant or canopy model
      !>
      !> 7. soil phys jsbach output file (hydro)
      !>
      ELSE 
        CALL add_output_file(p_output, 'HYDRO', p_current_file)
        CALL add_pool_variable(p_current_file, 'rootzone_soil_moisture', dsl4jsb_var(HYDRO_,wtr_rootzone), unit_hydro_pool)
        CALL add_sl_pool_variable(p_current_file, 'water_content_soil', dsl4jsb_var(HYDRO_,wtr_soil_sl), &
          &  unit=unit_hydro_pool, unit_sl=unit_hydro_pool)
        SELECT CASE (model%config%qmodel_id)
        CASE (QLAND, QPLANT, QCANOPY)
          CALL add_pool_variable(p_current_file, 'snow_height', dsl4jsb_var(HYDRO_,weq_snow_soil), unit_hydro_pool)
          CALL add_flux_variable(p_current_file, 'transpiration', dsl4jsb_var(HYDRO_,transpiration), unit = unit_hydro_flux)
          CALL add_flux_variable(p_current_file, 'evapopot', dsl4jsb_var(HYDRO_,evapopot), unit = unit_hydro_flux)
          CALL add_flux_variable(p_current_file, 'evaporation', dsl4jsb_var(HYDRO_,evapo_skin), unit = unit_hydro_flux)
          CALL add_flux_variable(p_current_file, 'srf_runoff', dsl4jsb_var(HYDRO_,runoff), unit = unit_hydro_flux)
          CALL add_flux_variable(p_current_file, 'gw_runoff', dsl4jsb_var(HYDRO_,drainage), unit = unit_hydro_flux)
        END SELECT ! ** if land, plant or canopy model
      END IF 
    END SELECT ! ** if land, plant, canopy or soil model

    ! ----------------------------------------------------------------------------------------------------- !
    !>
    !> 8. vegetation output file (veg)
    !>
    IF(model%config%flag_write_veg_output) THEN
      nr_of_veg_compartments = get_number_of_veg_compartments()

      SELECT CASE (model%config%qmodel_id)
      CASE (QLAND, QPLANT, QCANOPY, QSOIL)
        CALL add_output_file(p_output, 'VEG', p_current_file)

        SELECT CASE (model%config%qmodel_id)
        CASE (QLAND, QPLANT)
          !*** pools
          l_acc = .FALSE.  ! pool variables are not accumulated

          CALL add_flux_variable(p_current_file, 'net_biosphere_production', dsl4jsb_var(VEG_,net_biosphere_production))

          !--- previously: vegpoolC
          CALL add_all_parts_of_veg_bgcm_for_elem(p_current_file, 'vegpool_', ELEM_C_ID, nr_of_veg_compartments, l_acc,  &
            & dsl4jsb_bgcm(VEG_,vegbpool), model%config)
          CALL add_pool_variable(p_current_file, 'seed_bed_c',    dsl4jsb_elem_of_1L_bgcm(VEG_,seed_bed_pool,carbon))
          CALL add_pool_variable(p_current_file, 'total_veg_c',   dsl4jsb_var(VEG_, veg_pool_total_c), unit_pool)
          CALL add_pool_variable(p_current_file, 'total_veg_c13', dsl4jsb_var(VEG_, veg_pool_total_c13), unit_pool)
          CALL add_pool_variable(p_current_file, 'total_veg_c14', dsl4jsb_var(VEG_, veg_pool_total_c14), unit_pool)

          !--- previously: vegpoolC13C14 / vegpooldC13DC14
          CALL add_all_parts_of_veg_bgcm_for_elem(p_current_file, 'vegpool_', ELEM_C13_ID, nr_of_veg_compartments, l_acc, &
            & dsl4jsb_bgcm(VEG_,vegbpool), model%config)
          CALL add_pool_variable(p_current_file, 'seed_bed_c13', dsl4jsb_elem_of_1L_bgcm(VEG_,seed_bed_pool,carbon13))
          CALL add_all_parts_of_veg_bgcm_for_elem(p_current_file, 'vegpool_', ELEM_C14_ID, nr_of_veg_compartments, l_acc, &
            & dsl4jsb_bgcm(VEG_,vegbpool), model%config)
          CALL add_pool_variable(p_current_file, 'seed_bed_c14', dsl4jsb_elem_of_1L_bgcm(VEG_,seed_bed_pool,carbon14))

          !--- previously: vegpoolN
          CALL add_all_parts_of_veg_bgcm_for_elem(p_current_file, 'vegpool_', ELEM_N_ID, nr_of_veg_compartments, l_acc, &
            & dsl4jsb_bgcm(VEG_,vegbpool), model%config)
          CALL add_pool_variable(p_current_file, 'seed_bed_n', dsl4jsb_elem_of_1L_bgcm(VEG_,seed_bed_pool,nitrogen))
          CALL add_pool_variable(p_current_file, 'total_veg_n', dsl4jsb_var(VEG_, veg_pool_total_n), unit_pool)
          !--- previously: vegpoolP
          CALL add_all_parts_of_veg_bgcm_for_elem(p_current_file, 'vegpool_', ELEM_P_ID, nr_of_veg_compartments, l_acc, &
            & dsl4jsb_bgcm(VEG_,vegbpool), model%config)
          CALL add_pool_variable(p_current_file, 'seed_bed_p', dsl4jsb_elem_of_1L_bgcm(VEG_,seed_bed_pool,phosphorus))
          CALL add_pool_variable(p_current_file, 'total_veg_p', dsl4jsb_var(VEG_, veg_pool_total_p), unit_pool)
          !--- previously: vegpooldN15 / vegpoolN15
          CALL add_all_parts_of_veg_bgcm_for_elem(p_current_file, 'vegpool_', ELEM_N15_ID, nr_of_veg_compartments, l_acc, &
            & dsl4jsb_bgcm(VEG_,vegbpool), model%config)
          CALL add_pool_variable(p_current_file, 'seed_bed_n15', dsl4jsb_elem_of_1L_bgcm(VEG_,seed_bed_pool,nitrogen15))
          CALL add_pool_variable(p_current_file, 'total_veg_n15', dsl4jsb_var(VEG_, veg_pool_total_n15), unit_pool)
          ! ecopoolCNP -> does currently not work
          !@TODO SZ: generate new variables to have the total for organic vegetation, litter, soil. -> #96
          !     Ecosystem stores only in postprocessing -> #93
          ! WRITE(lun,'(a)') 'Time TotalSoilOrgC TotalSoilOrgN TotalSoilInOrgN EcoC EcoN EcoP EcoN15'
          !--- previously: veg_diagnostics
          !@TODO SZ: If we go for VEG_flux and VEG_state, then this should be in the state file -> #99
          CALL add_pool_variable(p_current_file, 'height', dsl4jsb_var(VEG_,height), unit_size)
          CALL add_pool_variable(p_current_file, 'diameter', dsl4jsb_var(VEG_,diameter), unit_size)
          CALL add_pool_variable(p_current_file, 'ind_density', dsl4jsb_var(VEG_,dens_ind), '# m-2')
          !@TODO: target_cn_leaf does not have a unit in mem%Add_var call in mo_veg_memory_class -> #100
          CALL add_pool_variable(p_current_file, 'cn_target', dsl4jsb_var(VEG_,target_cn_leaf), unit_cn)
          CALL add_pool_variable(p_current_file, 'cost_nup', dsl4jsb_var(VEG_,cost_n_uptake_root), unit_cn)
          CALL add_pool_variable(p_current_file, 'n_demand', dsl4jsb_var(VEG_,demand_uptake_n_tuptake_mavg), unitless)
          CALL add_pool_variable(p_current_file, 'p_demand', dsl4jsb_var(VEG_,demand_uptake_p_tuptake_mavg), unitless) 
          !@TODO SZ2: sinklim should be called beta_sinklim and added to the photosynthesis output stream -> #99
          ! !!! no, its a variable in the veg memory if it should go to assimi than the variable would need to go there
          CALL add_pool_variable(p_current_file, 'beta_sinklim', dsl4jsb_var(VEG_,beta_sinklim_ps_tacclim_mavg), unitless)
          ! according to #42 the following two variables should be kept:
          !SZ: this is an entry to the allocation calculation. Unit mol N mol-1 C, IMHO this can be removed. TBD @ QUINCY meeting
          !@TODO: check: growth_req_n_tuptake_mavg does not have a unit in mem%Add_var call in mo_veg_memory_class
          !       + upon declaration another unit is suggested "mol N m-2 s-1" => clarify! -> #100
          CALL add_pool_variable(p_current_file, 'growth_req_n_upt_mavg', dsl4jsb_var(VEG_,growth_req_n_tuptake_mavg), unit_nc)
          !@TODO: check: growth_req_p_tuptake_mavg does not have a unit in mem%Add_var call in mo_veg_memory_class
          !       + upon declaration another unit is suggested "mol P m-2 s-1" => clarify! -> #100
          !SZ: this is an entry to the allocation calculation. Unit mol P mol-1 N, IMHO this can be removed. TBD @ QUINCY meeting
          CALL add_pool_variable(p_current_file, 'growth_req_p_upt_mavg', dsl4jsb_var(VEG_,growth_req_p_tuptake_mavg), unit_pn)
          CALL add_pool_variable(p_current_file, 'growth_p_limit_based_on_n_mavg', &
            & dsl4jsb_var(VEG_,growth_p_limit_based_on_n_mavg), unit_pn)

          !*** fluxes
          l_acc = .TRUE.  ! flux variables are accumulated
          !--- previously: vegfluxN
          CALL add_flux_variable(p_current_file, 'uptake_nh4', dsl4jsb_var(VEG_,uptake_nh4))
          CALL add_flux_variable(p_current_file, 'uptake_no3', dsl4jsb_var(VEG_,uptake_no3))
          !@TODO: how to    RecyclingN (sum of recycling_leaf_n recycling_fine_root_n and recycling_heart_wood_n)
          CALL add_flux_variable(p_current_file, 'reserve_use_n', dsl4jsb_elem_of_1L_bgcm(VEG_,vegbreserve_use,nitrogen))
          CALL add_flux_variable(p_current_file, 'fixation_n', dsl4jsb_var(VEG_,n_fixation))
          CALL add_flux_variable(p_current_file, 'fixation_n15', dsl4jsb_var(VEG_,n_fixation_n15))
          CALL add_flux_variable(p_current_file, 'uptake_nh4_n15', dsl4jsb_var(VEG_,uptake_nh4_n15))
          CALL add_flux_variable(p_current_file, 'uptake_no3_n15', dsl4jsb_var(VEG_,uptake_no3_n15))
          CALL add_flux_sum_for_veg_bgcm_element(p_current_file, 'growth_n', ELEM_N_ID, &
            & dsl4jsb_bgcm(VEG_,vegbgrowth), model%config)
          CALL add_element_of_veg_bgcm_compartment(p_current_file, 'growth_sapwood_n15', &
            & dsl4jsb_bgcm(VEG_,vegbgrowth), ix_sap_wood, ELEM_N15_ID, model%config, USE_VEG_FLUX_FACT)
          CALL add_all_parts_of_veg_bgcm_for_elem(p_current_file, 'growth_', ELEM_N_ID, nr_of_veg_compartments, l_acc, &
            & dsl4jsb_bgcm(VEG_,vegbgrowth), model%config)
          CALL add_all_parts_of_veg_bgcm_for_elem(p_current_file, 'litterfall_', ELEM_N_ID, nr_of_veg_compartments, l_acc, &
            & dsl4jsb_bgcm(VEG_,vegblitterfall), model%config)
          CALL add_flux_variable(p_current_file, 'seed_bed_growth_n', dsl4jsb_elem_of_1L_bgcm(VEG_,seed_bed_growth,nitrogen))
          CALL add_flux_variable(p_current_file, 'seed_bed_litterfall_n', &
            & dsl4jsb_elem_of_1L_bgcm(VEG_,seed_bed_litterfall,nitrogen))
          !--- previously: vegfluxP
          CALL add_flux_variable(p_current_file, 'uptake_po4', dsl4jsb_var(VEG_,uptake_po4))
          !@TODO -- howto:
          !  recycling_p_acc      = recycling_p_acc       + (recycling_leaf_p(ipoi) + recycling_fine_root_p(ipoi) + &
          !     recycling_heart_wood_p(ipoi))  * dtime_one *1000000._wp/dtime
          CALL add_flux_variable(p_current_file, 'reserve_use_p', &
            & dsl4jsb_elem_of_1L_bgcm(VEG_,vegbreserve_use,phosphorus), USE_VEG_FLUX_FACT)
          CALL add_flux_sum_for_veg_bgcm_element(p_current_file, 'growth_p', ELEM_P_ID, &
            & dsl4jsb_bgcm(VEG_,vegbgrowth), model%config)

          CALL add_all_parts_of_veg_bgcm_for_elem(p_current_file, 'growth_', ELEM_P_ID, nr_of_veg_compartments, l_acc, &
            & dsl4jsb_bgcm(VEG_,vegbgrowth), model%config)
          CALL add_all_parts_of_veg_bgcm_for_elem(p_current_file, 'litterfall_', ELEM_P_ID, nr_of_veg_compartments, l_acc, &
            & dsl4jsb_bgcm(VEG_,vegblitterfall), model%config)
          CALL add_flux_variable(p_current_file, 'seed_bed_growth_p', dsl4jsb_elem_of_1L_bgcm(VEG_,seed_bed_growth,phosphorus))
          CALL add_flux_variable(p_current_file, 'seed_bed_litterfall_p', &
            & dsl4jsb_elem_of_1L_bgcm(VEG_,seed_bed_litterfall,phosphorus))

          IF ((model%lctlib(tile%lcts(1)%lib_id)%ForestFlag .AND. dsl4jsb_Config(Q_SYL_)%active) &
              & .OR. (model%lctlib(tile%lcts(1)%lib_id)%CropFlag .AND. dsl4jsb_Config(Q_AGR_)%active)) THEN
            CALL add_all_parts_of_veg_bgcm_for_elem(p_current_file, 'harvest_litter_', ELEM_C_ID, nr_of_veg_compartments, &
              & l_acc, dsl4jsb_bgcm(VEG_,harvest_litter), model%config)
            CALL add_all_parts_of_veg_bgcm_for_elem(p_current_file, 'harvest_litter_', ELEM_N_ID, nr_of_veg_compartments, &
              & l_acc, dsl4jsb_bgcm(VEG_,harvest_litter), model%config)
          END IF

          IF (dsl4jsb_Config(VEG_)%l_use_product_pools .AND. tile%lcts(1)%lib_id > 0) THEN
            CALL add_pool_variable(p_current_file, 'veg_products_total_c', dsl4jsb_var(VEG_,veg_products_total_c))
            CALL add_pool_variable(p_current_file, 'veg_products_total_n', dsl4jsb_var(VEG_,veg_products_total_n))

            IF ((model%lctlib(tile%lcts(1)%lib_id)%ForestFlag .AND. dsl4jsb_Config(Q_SYL_)%active) &
                & .OR. (model%lctlib(tile%lcts(1)%lib_id)%CropFlag .AND. dsl4jsb_Config(Q_AGR_)%active)) THEN

              CALL add_flux_variable(p_current_file, 'prod_decay_c', &
                & dsl4jsb_elem_of_1L_bgcm(VEG_,prod_decay,carbon), USE_VEG_FLUX_FACT)
              CALL add_flux_variable(p_current_file, 'prod_decay_n', &
                & dsl4jsb_elem_of_1L_bgcm(VEG_,prod_decay,nitrogen), USE_VEG_FLUX_FACT)
              CALL add_flux_variable(p_current_file, 'prod_decay_p', &
                & dsl4jsb_elem_of_1L_bgcm(VEG_,prod_decay,phosphorus), USE_VEG_FLUX_FACT)
              CALL add_flux_variable(p_current_file, 'prod_decay_c13', &
                & dsl4jsb_elem_of_1L_bgcm(VEG_,prod_decay,carbon13), USE_VEG_FLUX_FACT)
              CALL add_flux_variable(p_current_file, 'prod_decay_c14', &
                & dsl4jsb_elem_of_1L_bgcm(VEG_,prod_decay,carbon14), USE_VEG_FLUX_FACT)
              CALL add_flux_variable(p_current_file, 'prod_decay_n15', &
                & dsl4jsb_elem_of_1L_bgcm(VEG_,prod_decay,nitrogen15), USE_VEG_FLUX_FACT)

              CALL add_flux_variable(p_current_file, 'harvest_to_prod', &
                & dsl4jsb_elem_of_1L_bgcm(VEG_,harvest_to_prod,carbon), USE_VEG_FLUX_FACT)

              ! crop pp bgcm only available for crop tiles
              IF (model%lctlib(tile%lcts(1)%lib_id)%CropFlag) THEN
                CALL add_all_elements_of_1L_veg_bgcm_pool(p_current_file, 'pp_crop', dsl4jsb_bgcm(VEG_,veg_pp_crop))
              END IF

              ! wood pp bgcms only available for forest tiles
              IF (model%lctlib(tile%lcts(1)%lib_id)%ForestFlag) THEN
                CALL add_all_elements_of_1L_veg_bgcm_pool(p_current_file, 'pp_fuel', dsl4jsb_bgcm(VEG_,veg_pp_fuel))
                CALL add_all_elements_of_1L_veg_bgcm_pool(p_current_file, 'pp_paper', dsl4jsb_bgcm(VEG_,veg_pp_paper))
                CALL add_all_elements_of_1L_veg_bgcm_pool( &
                  & p_current_file, 'pp_fiberboard', dsl4jsb_bgcm(VEG_,veg_pp_fiberboard))
                CALL add_all_elements_of_1L_veg_bgcm_pool(p_current_file, 'pp_oirw', dsl4jsb_bgcm(VEG_,veg_pp_oirw))
                CALL add_all_elements_of_1L_veg_bgcm_pool(p_current_file, 'pp_pv', dsl4jsb_bgcm(VEG_,veg_pp_pv))
                CALL add_all_elements_of_1L_veg_bgcm_pool(p_current_file, 'pp_sawnwood', dsl4jsb_bgcm(VEG_,veg_pp_sawnwood))
              END IF
            END IF
          END IF
        END SELECT ! ** if land or plant model

        SELECT CASE (model%config%qmodel_id)
        CASE (QLAND, QPLANT, QCANOPY)
          CALL add_pool_variable(p_current_file, 'LAI', dsl4jsb_var(VEG_,lai), 'm2 m-2')
          !--- previously: veg_canopy
          CALL add_cl_pool_variable(p_current_file, 'f_et', dsl4jsb_var(VEG_,fn_et_cl), unit=unitless)
          CALL add_cl_pool_variable(p_current_file, 'f_rub', dsl4jsb_var(VEG_,fn_rub_cl), unit=unitless)
          CALL add_cl_pool_variable(p_current_file, 'f_chl', dsl4jsb_var(VEG_,fn_chl_cl), unit=unitless)
          CALL add_cl_pool_variable(p_current_file, 'f_pepc', dsl4jsb_var(VEG_,fn_pepc_cl), unit=unitless)
          CALL add_cl_pool_variable(p_current_file, 'f_oth', dsl4jsb_var(VEG_,fn_oth_cl), unit=unitless)
          CALL add_cl_pool_variable(p_current_file, 'leaf_nitrogen', dsl4jsb_var(VEG_,leaf_nitrogen_cl), unit="mmol m-2 LAI")
          !--- previously: vegfluxC
          l_acc = .TRUE.  ! flux variables are accumulated
          CALL add_flux_variable(p_current_file, 'growth_resp', dsl4jsb_var(VEG_,growth_respiration))
          CALL add_flux_variable(p_current_file, 'maint_resp', dsl4jsb_var(VEG_,maint_respiration))
          CALL add_flux_variable(p_current_file, 'n_transform_resp', dsl4jsb_var(VEG_,n_transform_respiration))
          CALL add_flux_variable(p_current_file, 'n_fixation_resp', dsl4jsb_var(VEG_,n_fixation_respiration))
          CALL add_flux_variable(p_current_file, 'npp', dsl4jsb_var(VEG_,npp))
          CALL add_flux_variable(p_current_file, 'leaf_herbivory', dsl4jsb_var(VEG_,herbivory_leaf_resp))
          CALL add_flux_variable(p_current_file, 'fruit_herbivory', dsl4jsb_var(VEG_,herbivory_fruit_resp))
          CALL add_flux_variable(p_current_file, 'reserve_use_c', & 
            & dsl4jsb_elem_of_1L_bgcm(VEG_,vegbreserve_use,carbon), USE_VEG_FLUX_FACT)
          CALL add_flux_variable(p_current_file, 'root_Exudation_c', &
            & dsl4jsb_elem_of_1L_bgcm(VEG_,vegbexudation,carbon), USE_VEG_FLUX_FACT)
          CALL add_flux_sum_for_veg_bgcm_element(p_current_file, 'growth_c', ELEM_C_ID, &
            & dsl4jsb_bgcm(VEG_,vegbgrowth), model%config)
          CALL add_all_parts_of_veg_bgcm_for_elem(p_current_file, 'growth_', ELEM_C_ID, nr_of_veg_compartments, l_acc, &
            & dsl4jsb_bgcm(VEG_,vegbgrowth), model%config)
          CALL add_all_parts_of_veg_bgcm_for_elem(p_current_file, 'litterfall_', ELEM_C_ID, nr_of_veg_compartments, l_acc, &
            & dsl4jsb_bgcm(VEG_,vegblitterfall), model%config)
          CALL add_flux_variable(p_current_file, 'seed_bed_growth_c', dsl4jsb_elem_of_1L_bgcm(VEG_,seed_bed_growth,carbon))
          CALL add_flux_variable(p_current_file, 'seed_bed_litterfall_c', &
            & dsl4jsb_elem_of_1L_bgcm(VEG_,seed_bed_litterfall,carbon))
          CALL add_flux_variable(p_current_file, 'growth_resp_c13', dsl4jsb_var(VEG_,growth_respiration_c13))
          CALL add_flux_variable(p_current_file, 'growth_resp_c14', dsl4jsb_var(VEG_,growth_respiration_c14))
          CALL add_flux_variable(p_current_file, 'maint_resp_c13', dsl4jsb_var(VEG_,maint_respiration_c13))
          CALL add_flux_variable(p_current_file, 'maint_resp_c14', dsl4jsb_var(VEG_,maint_respiration_c14))
          CALL add_element_of_veg_bgcm_compartment(p_current_file, 'growth_sapwood_c13', &
            & dsl4jsb_bgcm(VEG_,vegbgrowth), ix_sap_wood, ELEM_C13_ID, model%config, USE_VEG_FLUX_FACT)
          CALL add_element_of_veg_bgcm_compartment(p_current_file, 'growth_sapwood_c14', &
            & dsl4jsb_bgcm(VEG_,vegbgrowth), ix_sap_wood, ELEM_C14_ID, model%config, USE_VEG_FLUX_FACT)

        END SELECT ! ** if land, plant or canopy model

        SELECT CASE (model%config%qmodel_id)
        CASE (QLAND, QSOIL)
          !--- previously: soilfluxC / soilfluxCdC13DC14 / soilfluxC_lay
          CALL add_flux_variable(p_current_file, 'litterfall_c', dsl4jsb_var(VEG_,veg_litterfall_total_c))
          CALL add_flux_variable(p_current_file, 'litterfall_c13', dsl4jsb_var(VEG_,veg_litterfall_total_c13))
          CALL add_flux_variable(p_current_file, 'litterfall_c14', dsl4jsb_var(VEG_,veg_litterfall_total_c14))
          CALL add_flux_variable(p_current_file, 'litterfall_n', dsl4jsb_var(VEG_,veg_litterfall_total_n))
          CALL add_flux_variable(p_current_file, 'litterfall_p', dsl4jsb_var(VEG_,veg_litterfall_total_p))
          CALL add_flux_variable(p_current_file, 'litterfall_n15', dsl4jsb_var(VEG_,veg_litterfall_total_n15))
        END SELECT ! ** if land or soil model

        ! ** only if land: myc_diagnostics -> IF(flag_mycorrhiza)
        !@TODO: howto...
        !  myc_respiration_c_acc = myc_respiration_c_acc + SUM(myc_respiration(ipoi,:) * soil_depth_sl(ipoi,:)) * &
        !                          molar_mass_C * dtime_one * 1000000._wp / dtime
        !  f_cm_som_acc          = f_cm_som_acc + organic_cn * molar_mass_C / molar_mass_N * &
        !                          SUM(mycorrhiza_uptake_norg_sl(ipoi,:) * soil_depth_sl(ipoi,:)) * &
        !                          dtime_one * 1000000._wp / dtime
      END SELECT ! ** if land, plant, canopy or soil model
    END IF ! IF model%config%flag_write_veg_output
    
    ! ----------------------------------------------------------------------------------------------------- !
    !>
    !> 9. plant hydraulic output file (phyd)
    !>
    IF (dsl4jsb_Config(PHYD_)%use_plant_hydraulics) THEN
      SELECT CASE (model%config%qmodel_id)
      CASE (QLAND, QPLANT)
        CALL add_output_file(p_output, 'PHYD', p_current_file)
        CALL add_flux_variable(p_current_file, 'psi_leaf', dsl4jsb_var(PHYD_ , psi_leaf), unit=unit_water_potential)
        CALL add_flux_variable(p_current_file, 'psi_stem', dsl4jsb_var(PHYD_ , psi_stem), unit=unit_water_potential)
      END SELECT ! ** if land or plant model
    END IF

    ! ----------------------------------------------------------------------------------------------------- !
    !>
    !> 10. land to atmosphere output file (l2a)
    !>
    SELECT CASE (model%config%qmodel_id)
    CASE (QLAND)
      CALL add_output_file(p_output, 'L2A', p_current_file)
      CALL add_pool_variable(p_current_file, 'q_total_c', dsl4jsb_var(L2A_,q_total_c))
      CALL add_flux_variable(p_current_file, 'q_c_conservation_test', &
        &                    dsl4jsb_var(L2A_,q_c_conservation_test), unit=unit_veg_pool_flux)
    END SELECT ! ** if land mode

    ! ----------------------------------------------------------------------------------------------------- !
    !>
    !> 11. quincy sylviculture output file (q_syl)
    !>
    IF (dsl4jsb_Config(Q_SYL_)%active .AND. model%lctlib(tile%lcts(1)%lib_id)%ForestFlag) THEN
      SELECT CASE (model%config%qmodel_id)
      CASE (QLAND)
        CALL add_output_file(p_output, 'Q_SYL', p_current_file)
        CALL add_flux_variable(p_current_file, 'fract_forest_harvest', dsl4jsb_var(Q_SYL_ , fract_forest_harvest), unit=unit_fraction)
      END SELECT ! ** if land model
    END IF

    ! ----------------------------------------------------------------------------------------------------- !
    !>
    !> 12. soil and snow energy output file (sse)
    !>
    SELECT CASE (model%config%qmodel_id)
    CASE (QLAND, QPLANT, QSOIL, QCANOPY)
      CALL add_output_file(p_output, 'SSE', p_current_file)
      CALL add_sl_pool_variable(p_current_file, 'rho', dsl4jsb_var(SSE_,bulk_dens_sl), l_write_1d=.FALSE., unit_sl=unit_bulk_dens)
      CALL add_sl_pool_variable(p_current_file, 'temperature_soil', dsl4jsb_var(SSE_,t_soil_sl), l_avg_sl=.TRUE., &
        & l_write_2d=.FALSE., unit=unit_temperature)
      IF (model%config%use_soil_phys_jsbach) THEN
        SELECT CASE (model%config%qmodel_id)
        CASE (QLAND, QPLANT, QCANOPY)
          ! NOTE: ground heat flux had negated negation in mo_qs_output_txt; now default (USE_NDAT_ACC)
          CALL add_flux_variable(p_current_file, 'qg', dsl4jsb_var(SSE_,grnd_hflx), unit=unit_energy) 
        END SELECT ! ** if land, plant or canopy model
      END IF
    END SELECT ! ** if land, plant, soil or canopy model

    ! ----------------------------------------------------------------------------------------------------- !
    !>
    !> 13. surface energy balance output file (seb)
    !>
    IF (model%config%use_soil_phys_jsbach) THEN
      SELECT CASE (model%config%qmodel_id)
      CASE (QLAND, QPLANT, QCANOPY)
        CALL add_output_file(p_output, 'SEB', p_current_file)
        CALL add_flux_variable(p_current_file, 'qh', dsl4jsb_var(SEB_,sensible_hflx), unit=unit_energy)
        CALL add_flux_variable(p_current_file, 'qle', dsl4jsb_var(SEB_,latent_hflx), unit=unit_energy)
      END SELECT ! ** if land, plant or canopy model
    END IF

  END SUBROUTINE define_nc4_output_files

  ! ======================================================================================================= !
  !>
  !> Add a new output file to the quincy outputfile collection
  !>
  SUBROUTINE add_output_file(quincy_output, file_name_tag, current_output_file)
    TYPE(t_quincy_output),              INTENT(INOUT):: quincy_output       !< quincy output structure
    CHARACTER(len=*),                   INTENT(IN)   :: file_name_tag       !< name tag of the file to add
    TYPE(t_quincy_outputfile), POINTER, INTENT(OUT)  :: current_output_file !< pointer to current output file
    ! ----------------------------------------------------------------------------------------------------- !
    TYPE(t_quincy_outputfile), POINTER  :: new_quincy_output_file, this_quincy_output_file
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':add_output_file'
    ! ----------------------------------------------------------------------------------------------------- !
    ALLOCATE(new_quincy_output_file)

    IF (ASSOCIATED(quincy_output%current_output_file)) THEN
      ! Assert that so far there is no file with the same file name tag
      this_quincy_output_file => quincy_output%first_output_file 
      DO WHILE(ASSOCIATED(this_quincy_output_file))
        IF (this_quincy_output_file%file_name_tag .EQ. file_name_tag) THEN
          WRITE(message_text,*) 'Attention: there is already a file with this name tag '//TRIM(file_name_tag)//', please check!'
          CALL finish(TRIM(routine),message_text)
        ENDIF
        this_quincy_output_file => this_quincy_output_file%next_output_file
      ENDDO

      quincy_output%current_output_file%next_output_file => new_quincy_output_file
      new_quincy_output_file%prev_output_file => quincy_output%current_output_file
    ELSE
      quincy_output%first_output_file => new_quincy_output_file
    ENDIF
    quincy_output%current_output_file => new_quincy_output_file

    quincy_output%nr_of_output_files = quincy_output%nr_of_output_files + 1
    quincy_output%current_output_file%file_name_tag = TRIM(file_name_tag)

    !> Assign pointer of input variables to QUINCY output type
    quincy_output%current_output_file%sel_output_variables => quincy_output%sel_output_variables

    !> all quincy output files have a time variable
    CALL add_variable_to_quincy_output_file(quincy_output%current_output_file, name_var_time, "?", .FALSE.)

    current_output_file => new_quincy_output_file
  END SUBROUTINE add_output_file

  ! ======================================================================================================= !
  !>
  !> Simplified specific interface to add a normal pool variable via the general add_variable_to_quincy_output_file
  !>
  SUBROUTINE add_pool_variable(quincy_output_file, name_tag, quincy_var, unit)
    TYPE(t_quincy_outputfile),   INTENT(INOUT), TARGET :: quincy_output_file !< current quincy output file
    CHARACTER(len=*),            INTENT(IN)            :: name_tag           !< name or name tag of the variable to add
    TYPE(t_jsb_var_real2d),      INTENT(IN), TARGET    :: quincy_var         !< the actual variable
    CHARACTER(len=*),  OPTIONAL, INTENT(IN)            :: unit               !< unit of the variable to add (else use default)
    ! ----------------------------------------------------------------------------------------------------- !
    LOGICAL, PARAMETER            :: l_acc         = .FALSE.  !< pool variables are not accumulated
    INTEGER, PARAMETER            :: out_fact_kind = NO_FACT  !< id of multiplicative factor var to output
    CHARACTER(len=:), ALLOCATABLE :: this_unit
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':add_pool_variable'
    ! ----------------------------------------------------------------------------------------------------- !
    IF (PRESENT(unit)) THEN
      this_unit = unit
    ELSE
      this_unit = unit_pool
    ENDIF
    
    CALL add_variable_to_quincy_output_file(quincy_output_file, name_tag, this_unit, l_acc, out_fact_kind, quincy_var)
  END SUBROUTINE add_pool_variable

  ! ======================================================================================================= !
  !>
  !> Simplified specific interface to add a soil layered pool variable via the general add_variable_to_quincy_output_file
  !>
  SUBROUTINE add_sl_pool_variable(quincy_output_file, name_tag, q_var_3d, l_avg_sl, l_write_2d, l_write_1d, unit, unit_sl)
    TYPE(t_quincy_outputfile),   INTENT(INOUT), TARGET :: quincy_output_file !< current quincy output file
    CHARACTER(len=*),            INTENT(IN)            :: name_tag           !< name or name tag of the variable to add
    TYPE(t_jsb_var_real3d),      INTENT(IN), TARGET    :: q_var_3d           !< the 3d variable
    LOGICAL,           OPTIONAL, INTENT(IN)            :: l_avg_sl           !< if true: avg, else (default!) sum over soil layers
    LOGICAL,           OPTIONAL, INTENT(IN)            :: l_write_2d         !< default true - if false do not write per sl
    LOGICAL,           OPTIONAL, INTENT(IN)            :: l_write_1d         !< default true - if false do not write avg
    CHARACTER(len=*),  OPTIONAL, INTENT(IN)            :: unit               !< unit of 1D variable (else use default)
    CHARACTER(len=*),  OPTIONAL, INTENT(IN)            :: unit_sl            !< unit of 2D variable (else use default)
    ! ----------------------------------------------------------------------------------------------------- !
    LOGICAL,          PARAMETER   :: l_acc         = .FALSE.    !< pool variables are not accumulated
    INTEGER,          PARAMETER   :: out_fact_kind = NO_FACT    !< id of multiplicative factor var to output
    CHARACTER(len=:), ALLOCATABLE :: this_unit
    CHARACTER(len=:), ALLOCATABLE :: this_sl_unit
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':add_sl_pool_variable'
    ! ----------------------------------------------------------------------------------------------------- !
    IF (PRESENT(unit)) THEN
      this_unit = unit
    ELSE
      this_unit = unit_pool ! default unit for averaged sl pool variables
    ENDIF

    IF (PRESENT(unit_sl)) THEN
      this_sl_unit = unit_sl
    ELSE
      this_sl_unit = unit_soillayer_pool ! default unit for sl pool variables
    ENDIF
    
    CALL add_variable_to_quincy_output_file(quincy_output_file, name_tag, this_unit, l_acc, out_fact_kind, &
      & q_var_3d=q_var_3d, l_avg_sl=l_avg_sl, l_write_2d=l_write_2d, l_write_1d=l_write_1d, unit_sl=this_sl_unit)
  END SUBROUTINE add_sl_pool_variable

  ! ======================================================================================================= !
  !>
  !> Simplified specific interface to add a canopy layered pool variable via the general add_variable_to_quincy_output_file
  !>
  SUBROUTINE add_cl_pool_variable(quincy_output_file, name_tag, q_var_3d, unit)
    TYPE(t_quincy_outputfile),   INTENT(INOUT), TARGET :: quincy_output_file !< current quincy output file
    CHARACTER(len=*),            INTENT(IN)            :: name_tag           !< name or name tag of the variable to add
    TYPE(t_jsb_var_real3d),      INTENT(IN), TARGET    :: q_var_3d           !< the 3d variable
    CHARACTER(len=*),            INTENT(IN)            :: unit               !< unit of the variable to add
    ! ----------------------------------------------------------------------------------------------------- !
    LOGICAL, PARAMETER  :: l_acc         = .FALSE.          !< pool variables are not accumulated
    INTEGER, PARAMETER  :: var_type      = CANOPY_LAYERED   !< id of canopy layered variables...
    INTEGER, PARAMETER  :: out_fact_kind = NO_FACT          !< id of multiplicative factor var to output
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':add_cl_pool_variable'
    ! ----------------------------------------------------------------------------------------------------- !
    CALL add_variable_to_quincy_output_file(quincy_output_file, name_tag, unit, &
      & l_acc, out_fact_kind, q_var_3d=q_var_3d, var_type=var_type )
  END SUBROUTINE add_cl_pool_variable

  ! ======================================================================================================= !
  !>
  !> Simplified specific interface to add a normal flux variable via the general add_variable_to_quincy_output_file
  !>
  SUBROUTINE add_flux_variable(quincy_output_file, name_tag, quincy_var, out_fact_kind, unit)
    TYPE(t_quincy_outputfile),  INTENT(INOUT), TARGET :: quincy_output_file !< current quincy output file
    CHARACTER(len=*),           INTENT(IN)            :: name_tag           !< name or name tag of the variable to add
    TYPE(t_jsb_var_real2d),     INTENT(IN), TARGET    :: quincy_var         !< the actual variable
    INTEGER,          OPTIONAL, INTENT(IN)            :: out_fact_kind      !< id of multiplicative factor var to output (else use default)
    CHARACTER(len=*), OPTIONAL, INTENT(IN)            :: unit               !< unit of the variable to add (else use default)
    ! ----------------------------------------------------------------------------------------------------- !
    CHARACTER(len=:), ALLOCATABLE :: this_flux_unit
    INTEGER                       :: this_fact_kind
    LOGICAL, PARAMETER            :: l_acc             = .TRUE.   !< flux variables are accumulated
    INTEGER, PARAMETER            :: default_fact_kind = USE_NDAT_ACC 
      !< since flux variables are accumulated the default (i.e. if not passed) fact kind should be
      !< to calculate the average over the output interval
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':add_flux_variable'
    ! ----------------------------------------------------------------------------------------------------- !
    IF (PRESENT(out_fact_kind)) THEN
      this_fact_kind = out_fact_kind
    ELSE
      this_fact_kind = default_fact_kind
    ENDIF

    ! the unit depends on the factor
    IF (PRESENT(unit)) THEN
      this_flux_unit = unit
    ELSE
      IF (this_fact_kind .EQ. USE_VEG_FLUX_FACT) THEN
        this_flux_unit = unit_veg_pool_flux
      ELSE
        this_flux_unit = unit_flux
      ENDIF
    ENDIF

    CALL add_variable_to_quincy_output_file(quincy_output_file, name_tag, this_flux_unit, l_acc, this_fact_kind, quincy_var)
  END SUBROUTINE add_flux_variable

  ! ======================================================================================================= !
  !>
  !> Simplified specific interface to add a soil flux variable via the general add_variable_to_quincy_output_file
  !>
  SUBROUTINE add_sl_flux_variable(quincy_output_file, name_tag, q_var_3d, & 
      &                             nsoil_w, out_fact_kind, l_write_2d, l_write_1d, unit_1d, unit_sl)
    TYPE(t_quincy_outputfile),   INTENT(INOUT), TARGET :: quincy_output_file  !< current quincy output file
    CHARACTER(len=*),            INTENT(IN)            :: name_tag            !< name or name tag of the variable to add
    TYPE(t_jsb_var_real3d),      INTENT(IN), TARGET    :: q_var_3d            !< the soil layered variable
    INTEGER,                     INTENT(IN)            :: nsoil_w             !< number of soil layers
    INTEGER,           OPTIONAL, INTENT(IN)            :: out_fact_kind       !< id of multiplicative factor var to output
    LOGICAL,           OPTIONAL, INTENT(IN)            :: l_write_2d          !< default true - if false do not write per sl
    LOGICAL,           OPTIONAL, INTENT(IN)            :: l_write_1d          !< default true - if false do not write avg
    CHARACTER(len=*),  OPTIONAL, INTENT(IN)            :: unit_1d             !< unit of the 1D variable
    CHARACTER(len=*),  OPTIONAL, INTENT(IN)            :: unit_sl             !< unit of the 2D variable
    ! ----------------------------------------------------------------------------------------------------- !
    LOGICAL, PARAMETER            :: l_acc             = .TRUE.               !< flux variables are accumulated
    INTEGER, PARAMETER            :: default_fact_kind = USE_NDAT_ACC
      !< since flux variables are accumulated the default (i.e. if not passed) fact kind should be
      !< to calculate the average over the output interval
    INTEGER                       :: this_fact_kind
    CHARACTER(len=:), ALLOCATABLE :: this_1D_flux_unit
    CHARACTER(len=:), ALLOCATABLE :: this_2D_flux_unit
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':add_sl_flux_variable'
    ! ----------------------------------------------------------------------------------------------------- !
    IF (PRESENT(out_fact_kind)) THEN
      this_fact_kind = out_fact_kind
    ELSE
      this_fact_kind = default_fact_kind
    ENDIF

    IF (PRESENT(unit_1d)) THEN
      this_1D_flux_unit = unit_1d 
    ELSE
      this_1D_flux_unit = unit_flux
    ENDIF

    IF (PRESENT(unit_sl)) THEN
      this_2D_flux_unit = unit_sl 
    ELSE
      this_2D_flux_unit = unit_soillayer_flux
    ENDIF

    CALL add_variable_to_quincy_output_file(quincy_output_file, name_tag, this_1D_flux_unit, l_acc, this_fact_kind, &
      & q_var_3d=q_var_3d, nsoil_w=nsoil_w, l_write_2d=l_write_2d, l_write_1d=l_write_1d, unit_sl=this_2D_flux_unit)
  END SUBROUTINE add_sl_flux_variable

  ! ======================================================================================================= !
  !>
  !> Simplified specific interface to add the compartement sum for a certain element in a veg bgcm
  !> ... Note: only veg compartments! I.e. does not include seeds!
  !>
  SUBROUTINE add_flux_sum_for_veg_bgcm_element(quincy_output_file, name_tag, element_id, veg_bgcm, model_config)
    TYPE(t_quincy_outputfile),        INTENT(INOUT), TARGET :: quincy_output_file !< current quincy output file
    CHARACTER(len=*),                 INTENT(IN)            :: name_tag           !< name or name tag of the variable to add
    INTEGER,                          INTENT(IN)            :: element_id         !< ID of the element to sum components for
    TYPE(t_veg_bgcm_with_components), INTENT(IN), TARGET    :: veg_bgcm
    TYPE(t_jsb_model_config),         INTENT(IN)            :: model_config       !< the model config
      !< the veg bgcm whose compartments are to be summed for the fiven element
    ! ----------------------------------------------------------------------------------------------------- !
    CHARACTER(len=*), PARAMETER :: unit          = unit_veg_pool_flux !< unit of fluxes from pool variables
    LOGICAL, PARAMETER          :: l_acc         = .TRUE.             !< pool variables representing fluxes are accumulated
    INTEGER, PARAMETER          :: out_fact_kind = USE_VEG_FLUX_FACT 
      !< id of multiplicative factor used -- apparently this is the one always used when veg pool compartments are summed
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':add_flux_sum_for_veg_bgcm_element'
    ! ----------------------------------------------------------------------------------------------------- !
    CALL add_variable_to_quincy_output_file(quincy_output_file, name_tag, unit, l_acc, out_fact_kind, &
      &  element_id=element_id, veg_bgcm=veg_bgcm, model_config=model_config)
  END SUBROUTINE add_flux_sum_for_veg_bgcm_element

  ! ======================================================================================================= !
  !>
  !> Routine that really adds a new variable to the passed quincy output file 
  !>
  SUBROUTINE add_variable_to_quincy_output_file(quincy_output_file, &
      &                                         name_tag, &
      &                                         unit, &
      &                                         l_is_acc_var, &
      &                                         out_fact_kind, &
      &                                         quincy_var, &
      &                                         element_id, &
      &                                         veg_bgcm, &
      &                                         q_var_3d, &
      &                                         l_avg_sl, &
      &                                         l_write_2d, &
      &                                         l_write_1d, &
      &                                         nsoil_w, &
      &                                         var_type, &
      &                                         unit_sl, &
      &                                         model_config)
    TYPE(t_quincy_outputfile),                  INTENT(INOUT), TARGET :: quincy_output_file !< current quincy output file
    CHARACTER(len=*),                           INTENT(IN)            :: name_tag           !< name or name tag of the variable to add
    CHARACTER(len=*),                           INTENT(IN)            :: unit               !< unit of the variable to add
    LOGICAL,                                    INTENT(IN)            :: l_is_acc_var       !< indicates if the variable should be accumulated
    INTEGER,                          OPTIONAL, INTENT(IN)            :: out_fact_kind      !< id of multiplicative factor var to output
    TYPE(t_jsb_var_real2d),           OPTIONAL, INTENT(IN),    TARGET :: quincy_var         !< the actual variable
    INTEGER,                          OPTIONAL, INTENT(IN)            :: element_id         !< ID of this element or the element to add components for
    TYPE(t_veg_bgcm_with_components), OPTIONAL, INTENT(IN),    TARGET :: veg_bgcm           !< to be accumulated veg bgcm
      !< -> t_veg_bgcm_with_components could probably be abstracted to t_jsb_pool and treated with select type below, if required
    TYPE(t_jsb_var_real3d),           OPTIONAL, INTENT(IN),    TARGET :: q_var_3d           !< the actual variable in case of 2D output var (i.e., var with a vertical grid)
    LOGICAL,                          OPTIONAL, INTENT(IN)            :: l_avg_sl           !< only for sl vars: if true average, else (default!) sum
    LOGICAL,                          OPTIONAL, INTENT(IN)            :: l_write_2d         !< only for sl vars: default true - if false do not write per sl
    LOGICAL,                          OPTIONAL, INTENT(IN)            :: l_write_1d         !< only for sl vars: default true - if false do not write avg
    INTEGER,                          OPTIONAL, INTENT(IN)            :: nsoil_w            !< only for 2D output (q_var_3d): number of soil layers
    INTEGER,                          OPTIONAL, INTENT(IN)            :: var_type           !< only for 2D output (q_var_3d): type of layered variable
    CHARACTER(len=*),                 OPTIONAL, INTENT(IN)            :: unit_sl            !< only for 2D output: unit of the 2D variable
    TYPE(t_jsb_model_config),         OPTIONAL, INTENT(IN)            :: model_config       !< the model config
    ! ----------------------------------------------------------------------------------------------------- !
    TYPE(t_quincy_out_var), POINTER :: new_quincy_output_variable
    CHARACTER(len=100), ALLOCATABLE :: compartment_names(:)
    INTEGER                         :: nr_of_dims, size_container, this_elem_idx
    CHARACTER(len=*),     PARAMETER :: routine = TRIM(modname)//':add_variable_to_quincy_output_file'
    ! ----------------------------------------------------------------------------------------------------- !
    
    IF (ASSOCIATED (quincy_output_file%sel_output_variables)) THEN
      IF (.NOT. ANY (quincy_output_file%sel_output_variables .EQ. TRIM(name_tag))) THEN 
        ! Could not find the variable the selected output variable list. Do not add this variable.
        RETURN
      END IF 
    END IF
    
    ! create the new output variable
    ALLOCATE(new_quincy_output_variable)

    ! and sort it into the linked list
    IF (ASSOCIATED(quincy_output_file%current_var)) THEN
      quincy_output_file%current_var%next_var => new_quincy_output_variable
      new_quincy_output_variable%prev_var => quincy_output_file%current_var
    ELSE
      quincy_output_file%first_var => new_quincy_output_variable
    ENDIF
    quincy_output_file%current_var => new_quincy_output_variable

    ! check for the type of the output variable: t_jsb_var_real2d, t_jsb_var_real3d or veg_pool_struct
    !
    nr_of_dims = 1 ! default case for all but q_var_3d

    IF (PRESENT(quincy_var)) THEN
      ! i.e. a 'normal' t_jsb_var_real2d variable is passed
      ! There are a few cases where two variables should be added for the output...
      ALLOCATE(quincy_output_file%current_var%var_container(1))
      quincy_output_file%current_var%var_container(1)%p => quincy_var

      ! Assertion: no mixing of t_jsb_var_real2d, t_jsb_var_real3d and veg_bgcm
      IF (PRESENT(q_var_3d)) THEN
        WRITE(message_text,*) 'Attention: both, a 2D and a 3D variable were passed, which does not make sense, please check!'
        CALL finish(TRIM(routine),message_text)
      ENDIF
      IF (PRESENT(veg_bgcm)) THEN
        WRITE(message_text,*) 'Attention: both, a 2D and a veg bgcm were passed, which does not make sense, please check!'
        CALL finish(TRIM(routine),message_text)
      ENDIF
    ELSEIF (PRESENT(q_var_3d)) THEN
      ! i.e. a 'normal' t_jsb_var_real3d variable is passed
      nr_of_dims = 2
      ALLOCATE(quincy_output_file%current_var%var_container(1))
      quincy_output_file%current_var%var_container(1)%p => q_var_3d

      ! define vertical grid of var - soil_layer is default
      IF (PRESENT(var_type)) THEN
        quincy_output_file%current_var%var_type = var_type
      ELSE
        ! If no var_type but a q_var_3d variable is passed, then the soil layered type is assumed!
        quincy_output_file%current_var%var_type = SOIL_LAYERED
      ENDIF

      SELECT CASE (quincy_output_file%current_var%var_type)
      CASE (SOIL_LAYERED)
        quincy_output_file%has_soil_layers = .TRUE.
      CASE (CANOPY_LAYERED)
        quincy_output_file%has_canopy_layers = .TRUE.
      CASE DEFAULT
        WRITE(message_text,*) 'Currently only soil or canopy layered variables are implemented, please check! Unexpected type: ' &
          & // int2string(var_type) // ' passed for variable ' // TRIM(quincy_output_file%current_var%var_name)
        CALL finish(TRIM(routine),message_text)
      END SELECT
      ! Assertion: no mixing of t_jsb_var_real3d and veg_bgcm
      IF (PRESENT(veg_bgcm)) THEN
        WRITE(message_text,*) 'Attention: both, a 3D and a veg bgcm were passed, which does not make sense, please check!'
        CALL finish(TRIM(routine),message_text)
      ENDIF

    ELSEIF (PRESENT(veg_bgcm)) THEN
      ! i.e. a quincy veg pool struct is passed
      IF (.NOT. PRESENT(element_id)) THEN
        WRITE(message_text,*) 'Attention: when a to be accumulated veg bgcm is passed,' &
          & //' an element id is expected to also be passed as argument - but was not, please check!'
        CALL finish(TRIM(routine),message_text)
      ENDIF
      this_elem_idx = get_element_index(element_id, model_config%elements_index_map, routine)

      size_container = get_number_of_veg_compartments()
      ALLOCATE(compartment_names(size_container))
      ALLOCATE(quincy_output_file%current_var%var_container(size_container))

      CALL collect_veg_bgcm_pointers_for_element(this_elem_idx, veg_bgcm, compartment_names, &
          &                                      quincy_output_file%current_var%var_container)

      DEALLOCATE(compartment_names)
    ENDIF

    quincy_output_file%current_var%nr_of_dims = nr_of_dims
    IF (nr_of_dims .GT. quincy_output_file%nr_of_dims) THEN
      quincy_output_file%nr_of_dims = nr_of_dims
    ENDIF

    ! creates and additonal vector for acc/avg output variables (to store individual values over output period)
    quincy_output_file%current_var%l_is_acc_var = l_is_acc_var
    IF (.NOT. l_is_acc_var) THEN
      quincy_output_file%current_var%var_name = TRIM(name_tag)
    ELSE
      ! for accumulated variables a tag is appended depending on the output factor, i.e. if its an average or accumulated
      quincy_output_file%current_var%var_name = TRIM(name_tag)//'_acc'
      IF (PRESENT(out_fact_kind)) THEN
        IF (out_fact_kind .EQ. USE_NDAT_ACC) THEN
          quincy_output_file%current_var%var_name = TRIM(name_tag)//'_avg'
        ENDIF
      ENDIF
  
      IF (nr_of_dims .EQ. 2) THEN
        IF (quincy_output_file%current_var%var_type .NE. SOIL_LAYERED) THEN
          WRITE(message_text,*) 'In case of an acc var with two dimensions, a soil layered variable is expected.' &
            & //' This was not the case for ' //TRIM(name_tag)// '! Please check!'
          CALL finish(TRIM(routine),message_text)
        ENDIF
        IF (.NOT. PRESENT(nsoil_w)) THEN
          WRITE(message_text,*) 'In case of an acc var with two dimensions, ' &
            & //' the routine to add variables expects to get the length of the dimension (soil layers) passed as argument.' &
            & //' This was not the case for ' //TRIM(name_tag)// '! Please check!'
          CALL finish(TRIM(routine),message_text)
        ENDIF
        ALLOCATE(quincy_output_file%current_var%accumulated_2D(nsoil_w))
        quincy_output_file%current_var%accumulated_2D(:) = 0._wp
      ENDIF
    ENDIF

    quincy_output_file%current_var%unit = TRIM(unit)
    IF (PRESENT(unit_sl)) THEN
      quincy_output_file%current_var%unit_sl = unit_sl
    ELSE
      quincy_output_file%current_var%unit_sl = quincy_output_file%current_var%unit
    ENDIF

    IF (PRESENT(out_fact_kind)) THEN
      quincy_output_file%current_var%out_fact_kind = out_fact_kind
    ENDIF
    IF (PRESENT(l_avg_sl)) THEN
      quincy_output_file%current_var%l_avg_sl = l_avg_sl
    ENDIF
    IF (PRESENT(l_write_2d)) THEN
      quincy_output_file%current_var%l_write_2d = l_write_2d
    ENDIF
    IF (PRESENT(l_write_1d)) THEN
      quincy_output_file%current_var%l_write_1d = l_write_1d
    ENDIF

  END SUBROUTINE add_variable_to_quincy_output_file

  ! ======================================================================================================= !
  !>
  !> For given veg bgcm add all compartments of a certain element as new variables to the given outputfile
  !> NOTE: since seeds are no veg pool compartements anymore these are now added 'manually' as own variables
  !>
  SUBROUTINE add_all_parts_of_veg_bgcm_for_elem(quincy_output_file, var_name_prefix, this_element_id, &
      &                         nr_of_veg_compartments, l_is_acc_var, quincy_veg_bgcm_struct, model_config)
    TYPE(t_quincy_outputfile),        INTENT(INOUT), TARGET :: quincy_output_file     !< current quincy output file
    CHARACTER(len=*),                 INTENT(IN)            :: var_name_prefix        !< variable name prefix for the output 
    INTEGER,                          INTENT(IN)            :: this_element_id        !< ID of the element to add components for
    INTEGER,                          INTENT(IN)            :: nr_of_veg_compartments !< number of veg pool compartments
    LOGICAL,                          INTENT(IN)            :: l_is_acc_var           !< indicates if the variable should be accumulated
    TYPE(t_veg_bgcm_with_components), INTENT(IN), TARGET    :: quincy_veg_bgcm_struct !< pointer to the quincy structure
    TYPE(t_jsb_model_config),         INTENT(IN)            :: model_config           !< the model config
    ! ----------------------------------------------------------------------------------------------------- !
    CHARACTER(len=100)            :: compartment_names(nr_of_veg_compartments)
    TYPE(t_jsb_var_p), TARGET     :: compartment_pointers(nr_of_veg_compartments)
    CHARACTER(len=:), ALLOCATABLE :: unit
    CLASS(t_jsb_var), POINTER     :: this_var
    INTEGER                       :: ind
    INTEGER                       :: out_fact_kind               !< id of multiplicative factor var to output
    INTEGER                       :: this_elem_idx               !< index of the queried element
    CHARACTER(len=*), PARAMETER   :: routine = TRIM(modname)//':add_all_parts_of_veg_bgcm_for_elem'
    ! ----------------------------------------------------------------------------------------------------- !

    ! Assert that the model is running with this element + get its index
    this_elem_idx = get_element_index(this_element_id, model_config%elements_index_map, routine)

    ! if accumulated a flux compartement is assumed, else a pool compartement
    IF (l_is_acc_var) THEN
      unit = unit_flux
      out_fact_kind = USE_VEG_FLUX_FACT
    ELSE
      unit = unit_pool
      out_fact_kind = NO_FACT
    ENDIF

    CALL collect_veg_bgcm_pointers_for_element(this_elem_idx, quincy_veg_bgcm_struct, compartment_names, compartment_pointers)

    DO ind = 1,SIZE(compartment_pointers)
      this_var => compartment_pointers(ind)%p
      SELECT TYPE (this_var)
      TYPE IS (t_jsb_var_real2d)
        CALL add_variable_to_quincy_output_file(quincy_output_file, &
          & TRIM(var_name_prefix)//TRIM(compartment_names(ind))//'_'//TRIM(get_short_name_for_element_id(this_element_id)), &
          & unit, l_is_acc_var, out_fact_kind, this_var)
      END SELECT
    ENDDO

  END SUBROUTINE add_all_parts_of_veg_bgcm_for_elem

  ! ======================================================================================================= !
  !>
  !> For given sb bgcm add all compartments of a certain element as new variables to the passed quincy outputfile
  !>
  SUBROUTINE add_all_parts_of_sb_bgcm_for_elem(quincy_output_file, var_name_prefix, this_element_id, &
      &                            quincy_sb_bgcm, model_config, nsoil_w, nr_of_sb_compartments, sb_model_scheme, l_write_2d)
    USE mo_sb_config_class, ONLY: simple_1d_compartments, SB_PART_RESIDUE_IDX, SB_PART_MICROBIAL_IDX
    ! ----------------------------------------------------------------------------------------------------- !
    TYPE(t_quincy_outputfile),     INTENT(INOUT), TARGET :: quincy_output_file    !< current quincy output file
    CHARACTER(len=*),              INTENT(IN)            :: var_name_prefix       !< variable name prefix for the compartments to add
    INTEGER,                       INTENT(IN)            :: this_element_id       !< ID of the element to add components for
    TYPE(t_sb_bgcm_with_components), INTENT(IN), TARGET  :: quincy_sb_bgcm        !< pointer to the quincy structure
    TYPE(t_jsb_model_config),      INTENT(IN)            :: model_config          !< the model config
    INTEGER,                       INTENT(IN)            :: nsoil_w               !< number of soil layers
    INTEGER,                       INTENT(IN)            :: nr_of_sb_compartments !< number of compartments of the sb pool
    CHARACTER(len=*),              INTENT(IN)            :: sb_model_scheme       !< used sb model scheme
    LOGICAL, OPTIONAL,             INTENT(IN)            :: l_write_2d            !< write 1D and 2D (if true - default!) else only 1D
    ! ----------------------------------------------------------------------------------------------------- !
    CHARACTER(len=100)          :: compartment_names(nr_of_sb_compartments)
    TYPE(t_jsb_var_p), TARGET   :: compartment_pointers(nr_of_sb_compartments)
    CLASS(t_jsb_var), POINTER   :: this_var
    INTEGER                     :: ind
    CHARACTER(len=*), PARAMETER :: unit           = unit_pool           !< assumed default unit
    CHARACTER(len=*), PARAMETER :: unit_sl        = unit_soillayer_pool !< assumed default unit
    LOGICAL,          PARAMETER :: l_is_acc_var   = .FALSE.             !< it is currently not implemented to accumulate compartments of sb pools
    INTEGER,          PARAMETER :: out_fact_kind  = NO_FACT             !< id of multiplicative factor var to output
    INTEGER                     :: this_elem_idx                        !< index of the queried element
    CHARACTER(len=100)          :: this_compartment_name                !< compartment name (possibly adapted for simple_1d)
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':add_all_parts_of_sb_bgcm_for_elem'
    ! ----------------------------------------------------------------------------------------------------- !

    ! Assert that the model is running with this element + get its index
    this_elem_idx = get_element_index(this_element_id, model_config%elements_index_map, routine)

    CALL collect_sb_bgcm_pointers_for_element(this_elem_idx, quincy_sb_bgcm, compartment_names, compartment_pointers)

    DO ind = 1,SIZE(compartment_pointers)
      IF((TRIM(sb_model_scheme) == "simple_1d") .AND. .NOT. ANY(ind == simple_1d_compartments)) THEN
        ! For the simple_1d sb model scheme certain compartments are skipped (those not contained in simple_1d_compartments)
        CYCLE
      ENDIF
      this_var => compartment_pointers(ind)%p
      this_compartment_name = TRIM(compartment_names(ind))
      ! For the simple_1d sb model scheme two compartments are renamed for the output
      IF((TRIM(sb_model_scheme) == "simple_1d") .AND. (ind == SB_PART_RESIDUE_IDX)) THEN
        this_compartment_name = "slow"
      ENDIF
      IF((TRIM(sb_model_scheme) == "simple_1d") .AND. (ind == SB_PART_MICROBIAL_IDX)) THEN
        this_compartment_name = "fast"
      ENDIF

      SELECT TYPE (this_var)
      TYPE IS (t_jsb_var_real3d)
        CALL add_variable_to_quincy_output_file(quincy_output_file, &
          & TRIM(var_name_prefix)//TRIM(this_compartment_name)//'_'//TRIM(get_short_name_for_element_id(this_element_id)), &
          & unit, l_is_acc_var, out_fact_kind, q_var_3d=this_var, l_write_2d=l_write_2d, nsoil_w=nsoil_w, unit_sl=unit_sl)
      END SELECT
    ENDDO
  END SUBROUTINE add_all_parts_of_sb_bgcm_for_elem

  ! ======================================================================================================= !
  !>
  !> Add all elements for a passed one-layered veg pool bgcm (i.e. veg bgcm without compartments) to the passed quincy outputfile
  !>
  SUBROUTINE add_all_elements_of_1L_veg_bgcm_pool(quincy_output_file, var_name_prefix, this_veg_bgcm)
    TYPE(t_quincy_outputfile),      INTENT(INOUT), TARGET :: quincy_output_file   !< current quincy output file
    CHARACTER(len=*),               INTENT(IN)            :: var_name_prefix      !< variable name prefix 
    TYPE(t_veg_bgcm_with_elements), INTENT(IN), TARGET    :: this_veg_bgcm
    ! ----------------------------------------------------------------------------------------------------- !
    INTEGER                     :: idx
    INTEGER,          PARAMETER :: out_fact_kind  = NO_FACT
    CHARACTER(len=*), PARAMETER :: unit           = unit_pool
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':add_all_elements_of_1L_veg_bgcm_pool'
    LOGICAL,          PARAMETER :: l_is_acc_var   = .FALSE.
    ! ----------------------------------------------------------------------------------------------------- !

    ! Note: from the pool implementation I understand that the elements are added via Add_element
    !       to the pools element list and thus only the active elements are contained
    DO idx = 1, SIZE(this_veg_bgcm%element_list)
      SELECT TYPE (selected => this_veg_bgcm%element_list(idx)%p)
      CLASS IS (t_jsb_var_real2d)
        CALL add_variable_to_quincy_output_file(quincy_output_file, &
          & TRIM(var_name_prefix)//'_'//this_veg_bgcm%element_list(idx)%p%element_shortname, &
          & unit, l_is_acc_var, out_fact_kind, selected)
      END SELECT
    ENDDO
  END SUBROUTINE add_all_elements_of_1L_veg_bgcm_pool

  ! ======================================================================================================= !
  !>
  !> Add all elements for a passed one-layered sb bgcm (i.e. sb bgcm without compartments) to the passed quincy outputfile
  !>
  SUBROUTINE add_all_elements_of_1L_sb_bgcm(quincy_output_file, var_name_prefix, this_sb_bgcm, nsoil_w)
    TYPE(t_quincy_outputfile),     INTENT(INOUT), TARGET :: quincy_output_file   !< current quincy output file
    CHARACTER(len=*),              INTENT(IN)            :: var_name_prefix      !< variable name prefix 
    TYPE(t_sb_bgcm_with_elements), INTENT(IN), TARGET    :: this_sb_bgcm
    INTEGER,                       INTENT(IN)            :: nsoil_w              !< number of soil layers
    ! ----------------------------------------------------------------------------------------------------- !
    INTEGER                     :: idx
    INTEGER,          PARAMETER :: out_fact_kind  = USE_VEG_FLUX_FACT
    CHARACTER(len=*), PARAMETER :: unit           = unit_veg_pool_flux
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':add_all_elements_of_1L_sb_bgcm'
    ! ----------------------------------------------------------------------------------------------------- !

    ! Note: from the pool implementation I understand that the elements are added via Add_element
    !       to the pools element list and thus only the active elements are contained
    DO idx = 1, SIZE(this_sb_bgcm%element_list)
      SELECT TYPE (selected => this_sb_bgcm%element_list(idx)%p)
      CLASS IS (t_jsb_var_real3d)
        CALL add_sl_flux_variable(quincy_output_file,                                  &
          &  var_name_prefix//'_'//this_sb_bgcm%element_list(idx)%p%element_shortname, &
          & selected, nsoil_w, out_fact_kind=out_fact_kind, unit_1d=unit)
      END SELECT
    ENDDO
  END SUBROUTINE add_all_elements_of_1L_sb_bgcm

  ! ======================================================================================================= !
  !>
  !> Add given element for a passed compartment from the passed sb bgcm as soil flux variable to the passed quincy outputfile
  !>
  SUBROUTINE add_element_of_sb_bgcm_compartment(quincy_output_file, var_name_prefix, &
      &                                         quincy_sb_bgcm, this_compartment_ix, this_element_id, model_config, nsoil_w)
    TYPE(t_quincy_outputfile),       INTENT(INOUT), TARGET :: quincy_output_file  !< current quincy output file
    CHARACTER(len=*),                INTENT(IN)            :: var_name_prefix     !< variable name prefix
    TYPE(t_sb_bgcm_with_components), INTENT(IN), TARGET    :: quincy_sb_bgcm      !< pointer to the quincy structure
    INTEGER,                         INTENT(IN)            :: this_compartment_ix !< index of the compartment
    INTEGER,                         INTENT(IN)            :: this_element_id     !< ID of the element
    TYPE(t_jsb_model_config),        INTENT(IN)            :: model_config        !< the model config
    INTEGER,                         INTENT(IN)            :: nsoil_w             !< number of soil layers
    ! ----------------------------------------------------------------------------------------------------- !
    INTEGER :: this_elem_ix, this_part_ix
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':add_element_of_sb_bgcm_compartment'
    ! ----------------------------------------------------------------------------------------------------- !

    this_elem_ix = get_element_index(this_element_id, model_config%elements_index_map, routine)

    SELECT TYPE (selected => quincy_sb_bgcm%pool_list(this_compartment_ix)%p%element_list(this_elem_ix)%p)
    TYPE IS (t_jsb_var_real3d)
      CALL add_sl_flux_variable(quincy_output_file, var_name_prefix, selected, nsoil_w)
    END SELECT

  END SUBROUTINE add_element_of_sb_bgcm_compartment

  ! ======================================================================================================= !
  !>
  !> Add given element for a passed compartment from the passed veg bgcm as flux variable to the passed quincy outputfile
  !>
  SUBROUTINE add_element_of_veg_bgcm_compartment(quincy_output_file, var_name_prefix, &
    &                                           quincy_veg_bgcm, this_compartment_ix, this_element_id, model_config, out_fact_kind)
  TYPE(t_quincy_outputfile),       INTENT(INOUT), TARGET :: quincy_output_file  !< current quincy output file
  CHARACTER(len=*),                INTENT(IN)            :: var_name_prefix     !< variable name prefix
  TYPE(t_veg_bgcm_with_components), INTENT(IN), TARGET   :: quincy_veg_bgcm     !< pointer to the quincy structure
  INTEGER,                         INTENT(IN)            :: this_compartment_ix !< index of the compartment
  INTEGER,                         INTENT(IN)            :: this_element_id     !< ID of the element
  TYPE(t_jsb_model_config),        INTENT(IN)            :: model_config        !< the model config
  INTEGER,               OPTIONAL, INTENT(IN)            :: out_fact_kind       !< id of multiplicative factor var to output (else use default)
  ! ----------------------------------------------------------------------------------------------------- !
  INTEGER :: this_elem_ix, this_part_ix
  CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':add_element_of_veg_bgcm_compartment'
  ! ----------------------------------------------------------------------------------------------------- !

  this_elem_ix = get_element_index(this_element_id, model_config%elements_index_map, routine)

  SELECT TYPE (selected => quincy_veg_bgcm%pool_list(this_compartment_ix)%p%element_list(this_elem_ix)%p)
  TYPE IS (t_jsb_var_real2d)
    CALL add_flux_variable(quincy_output_file, var_name_prefix, selected, out_fact_kind)
  END SELECT

END SUBROUTINE add_element_of_veg_bgcm_compartment

  ! ======================================================================================================= !
  !>
  !> create nc4 output files and allocate buffers
  !>
  SUBROUTINE prepare_nc4_output_writing(tile, iblk, n_output_years, dtime, &
      &                                 output_interval_in_days, output_interval_char, add_tag, start_date_char)
    USE mo_io_units,  ONLY: filename_max
    dsl4jsb_Use_memory(HYDRO_)
    ! ----------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_tile_abstract), INTENT(IN) :: tile                    !< one tile with data structure for one lct
    INTEGER,                    INTENT(IN) :: iblk                    !< ... 1 in QUINCY ...
    INTEGER,                    INTENT(IN) :: n_output_years          !< number of years to be written into this file
    REAL(wp),                   INTENT(IN) :: dtime                   !< timestep length
    REAL(wp),                   INTENT(IN) :: output_interval_in_days !< interval in days (1.0 == one day exactly)
    CHARACTER(len=*),           INTENT(IN) :: output_interval_char    !< output interval tag for filename
    CHARACTER(len=*),           INTENT(IN) :: add_tag                 !< additional tag to add to the files -- e.g. 'spinup'
    CHARACTER(len=*),           INTENT(IN) :: start_date_char         !< start date of the simulation - used for the time axis
    ! ----------------------------------------------------------------------------------------------------- !
    TYPE(t_jsb_model),         POINTER     :: model
    TYPE(t_quincy_outputfile), POINTER     :: this_q_out_file
    TYPE(t_quincy_out_var),    POINTER     :: this_q_out_var
    TYPE(t_jsb_vgrid),         POINTER     :: vgrid_soil_w      !< Vertical grid
    TYPE(t_jsb_vgrid),         POINTER     :: vgrid_canopy      !< Vertical grid
    INTEGER                                :: nsoil_w, ncanopy, time_dimid, soil_level_dimid, canopy_level_dimid
    INTEGER                                :: ntimesteps        !< number of timesteps to be written into this file (time dim)
    INTEGER                                :: soil_depth_varid, buffer_size
    INTEGER                                :: dim_ids_1d(1), dim_ids_2d(2)
    REAL(wp),                  ALLOCATABLE :: soil_depth_sl(:)
    CHARACTER(len=filename_max)            :: this_file_name
    CHARACTER(len=30)                      :: runtime_date_as_str
    CHARACTER(len=*),          PARAMETER   :: file_name_ending        = ".nc"
    CHARACTER(len=*),          PARAMETER   :: UNITS                   = "units"
    CHARACTER(len=*),          PARAMETER   :: soil_depth_var_name     = "soil_depth"
    CHARACTER(len=*),          PARAMETER   :: soil_thickness_var_name = "soil_thickness"
    CHARACTER(len=*),          PARAMETER   :: soil_level_var_name     = "soil_level"
    CHARACTER(len=*),          PARAMETER   :: routine = TRIM(modname)//':prepare_nc4_output_writing'
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Def_memory(HYDRO_)
    ! ----------------------------------------------------------------------------------------------------- !
    model         => Get_model(tile%owner_model_id)
    vgrid_soil_w  => Get_vgrid('soil_depth_water')
    nsoil_w       =  vgrid_soil_w%n_levels
    vgrid_canopy  => Get_vgrid('q_canopy_layer')
    ncanopy       =  vgrid_canopy%n_levels
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Get_memory(HYDRO_)
    ! ----------------------------------------------------------------------------------------------------- !
    ALLOCATE(soil_depth_sl(nsoil_w))
    soil_depth_sl(:) = hydro__mem%soil_depth_sl%ptr(1,:,iblk)

    ! Assert expected output interval
    CALL assert_possible_output_interval(output_interval_in_days, output_interval_char, routine)
    buffer_size = NINT(one_year/output_interval_in_days)

    ! iterate over the linked list of quincy output files
    this_q_out_file => model%quincy_output%first_output_file
    DO WHILE(ASSOCIATED(this_q_out_file))
      this_file_name = TRIM(this_q_out_file%file_name_tag) //'_'// TRIM(add_tag) &
        &              // TRIM(output_interval_char) // file_name_ending

      ! Files that do not actually contain to be written variables do not need to be created and can be skipped
      ! (this would be the case when a file with selected output variables is used, but no variable of 
      !  this file is listed there -- in this case the file would only contain one variable - the time ...)
      IF (.NOT. ASSOCIATED(this_q_out_file%first_var%next_var)) THEN
        this_q_out_file => this_q_out_file%next_output_file
        CYCLE
      END IF

      CALL message(TRIM(routine), 'Setting up file '//TRIM(this_file_name))

#ifdef __QUINCY_WITH_NETCDF__
      !> 1. create each file in the linked list 
      !>
      CALL check_nc4_command_stat(nf90_create(TRIM(this_file_name), nf90_netcdf4, this_q_out_file%file_id), &
        &                         routine//':create:'//TRIM(this_file_name))
      ! create time dimension

      SELECT CASE (TRIM(output_interval_char))
      CASE ('yearly')
        ntimesteps = n_output_years
      CASE ('weekly')
        ntimesteps = n_output_years * 52
      CASE ('daily')
        ntimesteps = n_output_years * NINT(one_year)
      CASE ('timestep')
        ntimesteps = n_output_years * NINT(one_year * one_day / dtime)
      CASE DEFAULT
        CALL finish(TRIM(routine), ' '//TRIM(output_interval_char)//' .')
      END SELECT

      CALL check_nc4_command_stat(nf90_def_dim(this_q_out_file%file_id, 'time', ntimesteps, time_dimid), &
        &                         routine//':add:dim:time')
      ! Check for second dimension
      IF (this_q_out_file%nr_of_dims .EQ. 2) THEN
        IF (this_q_out_file%has_soil_layers) THEN
          CALL check_nc4_command_stat(nf90_def_dim(this_q_out_file%file_id, 'soil_layer', nsoil_w, soil_level_dimid), &
            &                         routine//':add:dim:soil_layer')
          dim_ids_2d(1) = soil_level_dimid
          IF (this_q_out_file%has_canopy_layers) THEN
            WRITE(message_text,*) 'Currently it is not implemented that a file has soil and canopy layered vars, but ' &
              & // TRIM(this_q_out_file%file_name_tag)// ' has?! Please check!'
          ENDIF
        ELSEIF (this_q_out_file%has_canopy_layers) THEN
          CALL check_nc4_command_stat(nf90_def_dim(this_q_out_file%file_id, 'q_canopy_layer', ncanopy, canopy_level_dimid), &
            &                         routine//':add:dim:q_canopy_layer')
          dim_ids_2d(1) = canopy_level_dimid
        ELSE
          WRITE(message_text,*) 'File needs to either be soil or canopy layered if it has 2 dimensions! But ' &
            & // TRIM(this_q_out_file%file_name_tag)// ' was neither nor. Please check!'
        ENDIF
        ! unlimited dimension must come last (i.e. here second) in the list of dimids
        dim_ids_2d(2) = time_dimid
      ELSEIF (this_q_out_file%nr_of_dims .NE. 1) THEN
        WRITE(message_text,*) 'Currently only files with 1 or 2 dimensions are implemented.'       &
          & //' However, number of dims was ' // int2string(this_q_out_file%nr_of_dims) // ' for ' &
          & // TRIM(this_q_out_file%file_name_tag)// '! Please check!'
        CALL finish(TRIM(routine),message_text)
      ENDIF
      dim_ids_1d(1) = time_dimid

      IF ((this_q_out_file%nr_of_dims .EQ. 2) .AND. (this_q_out_file%has_soil_layers)) THEN
        ! define soil layer variables
        CALL check_nc4_command_stat(nf90_def_var(this_q_out_file%file_id, TRIM(soil_depth_var_name), NF90_REAL, &
          &                         (/ soil_level_dimid /), soil_depth_varid), routine//':add:var:'//TRIM(soil_depth_var_name))
      ENDIF
      
      ! iterate over the linked list of variables for this file
      this_q_out_var => this_q_out_file%first_var
      DO WHILE(ASSOCIATED(this_q_out_var))
        CALL message(TRIM(routine), 'Add var: '// TRIM(this_q_out_var%var_name))

        !For time: set unit with respect to current interval
        IF (TRIM(this_q_out_var%var_name) .EQ. name_var_time) THEN
          this_q_out_var%unit = "seconds since "//TRIM(start_date_char)
        ENDIF

        !> 2. Define a variable for each variable in the list of this file
        !>
        IF (this_q_out_var%nr_of_dims .EQ. 1) THEN
          IF (TRIM(this_q_out_var%var_name) == name_var_time) THEN
            CALL check_nc4_command_stat(nf90_def_var(this_q_out_file%file_id, TRIM(this_q_out_var%var_name), NF90_DOUBLE, &
              &                         dim_ids_1d, this_q_out_var%var_id), routine//':add:var:'//TRIM(this_q_out_var%var_name))
          ELSE
            CALL check_nc4_command_stat(nf90_def_var(this_q_out_file%file_id, TRIM(this_q_out_var%var_name), NF90_REAL, &
              &                         dim_ids_1d, this_q_out_var%var_id), routine//':add:var:'//TRIM(this_q_out_var%var_name))
          ENDIF

        ELSEIF (this_q_out_var%nr_of_dims .EQ. 2) THEN
          IF (this_q_out_file%has_soil_layers) THEN
            this_q_out_var%dim2_len = nsoil_w
            ! Assertion
            IF ((.NOT. this_q_out_var%l_write_2d) .AND. (.NOT. this_q_out_var%l_write_1d)) THEN
              WRITE(message_text,*) 'Not writing 1d and not writing 2d makes no sense!' &
                & //' Problem occurred with var ' //TRIM(this_q_out_var%var_name)// '. Please check!'
              CALL finish(TRIM(routine),message_text)
            ENDIF

            !> 2.a if l_write_2d is true var is written per sl
            !>
            IF (this_q_out_var%l_write_2d) THEN
              CALL check_nc4_command_stat(nf90_def_var(this_q_out_file%file_id, TRIM(this_q_out_var%var_name)//"_L", NF90_REAL, &
                &             dim_ids_2d, this_q_out_var%var_add_2d_id), routine//':add:var:'//TRIM(this_q_out_var%var_name)//"_L")
              ! together with the unit
              CALL check_nc4_command_stat(nf90_put_att(this_q_out_file%file_id, this_q_out_var%var_add_2d_id, UNITS, &
               &                          this_q_out_var%unit_sl), routine//':add:unit')
            ENDIF
            !> 2.b if l_write_1d, var is written for the sum/avg
            !>
            IF (this_q_out_var%l_write_1d) THEN
              CALL check_nc4_command_stat(nf90_def_var(this_q_out_file%file_id, TRIM(this_q_out_var%var_name), NF90_REAL, &
                & dim_ids_1d, this_q_out_var%var_id), routine//':add:var:'//TRIM(this_q_out_var%var_name))
            ENDIF
          ELSEIF (this_q_out_file%has_canopy_layers) THEN
            this_q_out_var%dim2_len = ncanopy
            !> 2.c For canopy layered variables only the 2D variable is added
            !>
            CALL check_nc4_command_stat(nf90_def_var(this_q_out_file%file_id, TRIM(this_q_out_var%var_name), NF90_REAL, &
              &                         dim_ids_2d, this_q_out_var%var_id), routine//':add:var:'//TRIM(this_q_out_var%var_name))
          ELSE
            CALL finish(TRIM(routine), 'Neither soil layered nor canopy layered should have been caught already before!')
          ENDIF
        ELSE
          CALL finish(TRIM(routine), 'Unimplemented number of dimensions for variable: ' // int2string(this_q_out_file%nr_of_dims))
        ENDIF
        ! And its unit - exception: sl var for which no sum/avg should be written...
        IF (.NOT.(this_q_out_file%has_soil_layers .AND. .NOT.(this_q_out_var%l_write_1d))) THEN
          CALL check_nc4_command_stat(nf90_put_att(this_q_out_file%file_id, this_q_out_var%var_id, UNITS, &
            &                         this_q_out_var%unit), routine//':add:unit')
        ENDIF
        ! For time: add calendar
        IF (TRIM(this_q_out_var%var_name) .EQ. name_var_time) THEN
          IF (.NOT. one_year .EQ. 365._wp) THEN
            CALL finish(TRIM(routine), 'Attention, variable one_year was not as expected 365 days, but: ' // real2string(one_year))
          ENDIF
          CALL check_nc4_command_stat(nf90_put_att(this_q_out_file%file_id, this_q_out_var%var_id, 'calendar', &
            & TRIM(this_calendar_str)), routine//':add:calendar')
        ENDIF
        ! Allocate arrays for buffered writing
        IF (this_q_out_var%nr_of_dims == 1) THEN
          IF(ALLOCATED(this_q_out_var%buffer_2D)) DEALLOCATE(this_q_out_var%buffer_2D)
          ALLOCATE(this_q_out_var%buffer_2D(buffer_size))
          this_q_out_var%buffer_2D(:) = 0._wp
        ELSE IF(this_q_out_var%nr_of_dims == 2) THEN
          IF(ALLOCATED(this_q_out_var%buffer_3D)) DEALLOCATE(this_q_out_var%buffer_3D)
          ALLOCATE(this_q_out_var%buffer_3D(this_q_out_var%dim2_len, buffer_size))
          this_q_out_var%buffer_3D(:, :) = 0._wp
        ELSE
          WRITE(message_text,*) 'Currently only 1D and 2D variables, please check! Variable name: ' &
            & // TRIM(this_q_out_var%var_name) // ' of file with tag ' // TRIM(this_q_out_file%file_name_tag)
          CALL finish(TRIM(routine),message_text)
        END IF
        ! go to next var in the linked list
        this_q_out_var => this_q_out_var%next_var
      ENDDO

      ! End define mode for this file
      CALL check_nc4_command_stat(nf90_enddef(this_q_out_file%file_id), routine//':enddef:file' )

      IF ((this_q_out_file%nr_of_dims .EQ. 2) .AND. (this_q_out_file%has_soil_layers)) THEN
        ! Write soil layer variables
        CALL check_nc4_command_stat(nf90_put_var(this_q_out_file%file_id, soil_depth_varid, soil_depth_sl(:), &
          &                         count=(/ nsoil_w /)), routine//':write:soil_layer' )
      ENDIF

      ! Add more meta-data to the file
      CALL check_nc4_command_stat(nf90_put_att(this_q_out_file%file_id, NF90_GLOBAL, 'model_config', &
        &                         TRIM(Get_quincy_model_config_name(model%config%qmodel_id))), routine//'metadata:model_config')
      CALL check_nc4_command_stat(nf90_put_att(this_q_out_file%file_id, NF90_GLOBAL, 'branch', &
        &                         TRIM(model%config%git_branch)), routine//'metadata:branch')
      CALL check_nc4_command_stat(nf90_put_att(this_q_out_file%file_id, NF90_GLOBAL, 'commit', &
        &                         TRIM(model%config%git_commit_SHA)), routine//'metadata:commit')
      IF (model%config%code_has_changed_since_SHA) THEN
        CALL check_nc4_command_stat(nf90_put_att(this_q_out_file%file_id, NF90_GLOBAL, 'commit_status', 'modified'), &
          &                         routine//'metadata:commit_status')
      ENDIF
      CALL check_nc4_command_stat(nf90_put_att(this_q_out_file%file_id, NF90_GLOBAL, 'run_by', &
        &                         TRIM(model%config%simulation_run_by)), routine//'metadata:run_by')
      CALL FDATE(runtime_date_as_str)
      CALL check_nc4_command_stat(nf90_put_att(this_q_out_file%file_id, NF90_GLOBAL, 'date', &
        &                         TRIM(runtime_date_as_str)), routine//'metadata:date')
#endif
      ! Set record variable to first entry
      this_q_out_file%current_record = 1._wp

      ! Go to next file in the linked list
      this_q_out_file => this_q_out_file%next_output_file
    ENDDO

    DEALLOCATE(soil_depth_sl)
  END SUBROUTINE prepare_nc4_output_writing

  ! ======================================================================================================= !
  !>
  !> reset acc variables
  !>
  SUBROUTINE reset_nc_acc_vars(tile, options)
    dsl4jsb_Use_memory(SPQ_)
    ! ----------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_tile_abstract), INTENT(IN) :: tile              !< one tile with data structure for one lct
    TYPE(t_jsb_task_options),   INTENT(IN) :: options           !< model options
    ! ----------------------------------------------------------------------------------------------------- !
    TYPE(t_jsb_model),          POINTER     :: model
    TYPE(t_quincy_outputfile),  POINTER     :: this_q_out_file
    TYPE(t_quincy_out_var),     POINTER     :: this_q_out_var
    TYPE(t_jsb_vgrid),          POINTER     :: vgrid_soil_w     !< Vertical grid
    TYPE(t_jsb_vgrid),          POINTER     :: vgrid_canopy     !< Vertical grid
    INTEGER                                 :: nsoil_w          !< number of soil layers as used/defined by the SB_ process
    INTEGER                                 :: ncanopy          !< number of canopy layeres
    INTEGER                                 :: iblk

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':reset_nc_acc_vars'
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Def_memory(SPQ_)
    ! ----------------------------------------------------------------------------------------------------- !
    model         => Get_model(tile%owner_model_id)
    iblk          = options%iblk
    vgrid_soil_w  => Get_vgrid('soil_depth_water')
    nsoil_w       =  vgrid_soil_w%n_levels
    vgrid_canopy  => Get_vgrid('q_canopy_layer')
    ncanopy       =  vgrid_canopy%n_levels
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Get_memory(SPQ_)
    ! ----------------------------------------------------------------------------------------------------- !
    
    ! iterate over the linked list of quincy output files
    this_q_out_file => model%quincy_output%first_output_file
    DO WHILE(ASSOCIATED(this_q_out_file))

      ! Files that do not actually contain to be written variables have not been created and can be skipped
      IF (.NOT. ASSOCIATED(this_q_out_file%first_var%next_var)) THEN
        this_q_out_file => this_q_out_file%next_output_file
        CYCLE
      END IF

      this_q_out_var => this_q_out_file%first_var

      DO WHILE(ASSOCIATED(this_q_out_var))
        SELECT CASE (this_q_out_var%nr_of_dims)
        CASE (1)
          IF (this_q_out_var%l_is_acc_var) THEN
          this_q_out_var%accumulated = 0._wp
          ENDIF
        CASE (2)
          IF (this_q_out_var%l_is_acc_var) THEN
            this_q_out_var%accumulated_2D(:) = 0._wp
          ENDIF
        CASE DEFAULT
          WRITE(message_text,*) 'Currently only 1D and 2D variables, please check! Variable name: ' &
            & // TRIM(this_q_out_var%var_name) // ' of file with tag ' // TRIM(this_q_out_file%file_name_tag)
          CALL finish(TRIM(routine),message_text)
        END SELECT
        this_q_out_var => this_q_out_var%next_var
      ENDDO

      this_q_out_file => this_q_out_file%next_output_file
    ENDDO

    model%quincy_output%ndat_acc = 0._wp
  END SUBROUTINE reset_nc_acc_vars

  ! ======================================================================================================= !
  !>
  !> write buffered output
  !>
  SUBROUTINE write_netcdf4_quincy_buffered(tile, options, buffer_size)
    dsl4jsb_Use_memory(HYDRO_)
    ! ----------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_tile_abstract), INTENT(IN) :: tile              !< one tile with data structure for one lct
    TYPE(t_jsb_task_options),   INTENT(IN) :: options           !< model options
    INTEGER,                    INTENT(IN) :: buffer_size       !< size of the buffers to write to netcdf
    ! ----------------------------------------------------------------------------------------------------- !
    TYPE(t_jsb_model),             POINTER :: model
    TYPE(t_quincy_outputfile),     POINTER :: this_q_out_file
    TYPE(t_quincy_out_var),        POINTER :: this_q_out_var
    TYPE(t_jsb_vgrid),             POINTER :: vgrid_soil_w     !< Vertical grid
    TYPE(t_jsb_vgrid),             POINTER :: vgrid_canopy     !< Vertical grid
    INTEGER                                :: nsoil_w          !< number of soil layers as used/defined by the SB_ process
    INTEGER                                :: ncanopy          !< number of canopy layeres
    INTEGER                                :: iblk
    REAL(wp),                  ALLOCATABLE :: weighted_layer_value(:)
    REAL(wp),                  ALLOCATABLE :: soil_depth_sl(:)
    REAL(WP),                  ALLOCATABLE :: buffered_output_2D(:)

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':write_netcdf4_quincy_buffered'
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Def_memory(HYDRO_)
    ! ----------------------------------------------------------------------------------------------------- !
    model         => Get_model(tile%owner_model_id)
    iblk          = options%iblk
    vgrid_soil_w  => Get_vgrid('soil_depth_water')
    nsoil_w       =  vgrid_soil_w%n_levels
    vgrid_canopy  => Get_vgrid('q_canopy_layer')
    ncanopy       =  vgrid_canopy%n_levels
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Get_memory(HYDRO_)
    ! ----------------------------------------------------------------------------------------------------- !
    ALLOCATE(soil_depth_sl(nsoil_w))
    ALLOCATE(weighted_layer_value(buffer_size))

    soil_depth_sl(:) = hydro__mem%soil_depth_sl%ptr(1,:,iblk)
    
    ! iterate over the linked list of quincy output files
    this_q_out_file => model%quincy_output%first_output_file
    DO WHILE(ASSOCIATED(this_q_out_file))

      ! Files that do not actually contain to be written variables have not been created and can be skipped
      IF (.NOT. ASSOCIATED(this_q_out_file%first_var%next_var)) THEN
        this_q_out_file => this_q_out_file%next_output_file
        CYCLE
      END IF

      this_q_out_var => this_q_out_file%first_var
      ! First var is by convention time
      IF (.NOT. TRIM(this_q_out_var%var_name) .EQ. name_var_time) THEN
        WRITE(message_text,*) 'First var of all files is by convention'//name_var_time// ', please check! File with tag ' &
          & // TRIM(this_q_out_file%file_name_tag) // ' has as first var: ' // TRIM(this_q_out_var%var_name)
        CALL finish(TRIM(routine),message_text)
      ENDIF
      
      ! Write the time variable for each timestep
      CALL write_array_to_1d_var_on_file(this_q_out_var%buffer_2D, this_q_out_var%var_id, &
        & this_q_out_file%file_id, this_q_out_file%current_record, buffer_size)

      ! Move on to the next variable
      this_q_out_var => this_q_out_var%next_var

      DO WHILE(ASSOCIATED(this_q_out_var))
        SELECT CASE (this_q_out_var%nr_of_dims)
        CASE (1)
          buffered_output_2D = this_q_out_var%buffer_2D
          IF (this_q_out_var%l_write_1d) THEN
            ! Write the 1D buffer
            CALL write_array_to_1d_var_on_file(buffered_output_2D, this_q_out_var%var_id, &
              & this_q_out_file%file_id, this_q_out_file%current_record, buffer_size)

            this_q_out_var%accumulated = 0.0_wp
          ENDIF 

        CASE (2)
          IF (this_q_out_var%l_write_2d) THEN
            ! 2d variables can either be soil or canopy layered. The only difference here is 
            ! the netcdf variable id that is different for soil and canopy variables
            IF (this_q_out_var%var_type == SOIL_LAYERED) THEN
              CALL write_2d_array_to_2d_var_on_file(this_q_out_var%buffer_3D, this_q_out_var%var_add_2d_id, &
                & this_q_out_file%file_id, this_q_out_file%current_record, buffer_size, this_q_out_var%dim2_len)
            ELSEIF(this_q_out_var%var_type == CANOPY_LAYERED) THEN
              CALL write_2d_array_to_2d_var_on_file(this_q_out_var%buffer_3D, this_q_out_var%var_id, &
                & this_q_out_file%file_id, this_q_out_file%current_record, buffer_size, this_q_out_var%dim2_len)
            ENDIF
            ! Reset the accumulation variables if the output variable is of flux type
            IF (this_q_out_var%l_is_acc_var) THEN
              this_q_out_var%accumulated_2D = 0.0_wp
            ENDIF
          ENDIF

          IF (this_q_out_var%var_type .EQ. SOIL_LAYERED) THEN
            ! If demanded, write the 1d variable
            IF (this_q_out_var%l_write_1d) THEN
              ! either soil depth weighted layer sum or avg
              IF (this_q_out_var%l_avg_sl) THEN
                CALL calculate_soil_depth_weighted_layer_avg(this_q_out_var%buffer_3D, soil_depth_sl, weighted_layer_value)
              ELSE
                CALL calculate_soil_depth_weighted_layer_sum(this_q_out_var%buffer_3D, soil_depth_sl, weighted_layer_value)
              ENDIF
              CALL write_array_to_1d_var_on_file(weighted_layer_value(:), this_q_out_var%var_id, &
                &                                this_q_out_file%file_id, this_q_out_file%current_record, buffer_size)
            ENDIF
          ENDIF

        CASE DEFAULT
          WRITE(message_text,*) 'Currently only 1D and 2D variables are supported, please check! Variable name: ' &
            & // TRIM(this_q_out_var%var_name) // ' of file with tag ' // TRIM(this_q_out_file%file_name_tag)
          CALL finish(TRIM(routine),message_text)
        END SELECT

        this_q_out_var => this_q_out_var%next_var
      ENDDO
      this_q_out_file%current_record = this_q_out_file%current_record + buffer_size
      this_q_out_file => this_q_out_file%next_output_file
    ENDDO

    DEALLOCATE(soil_depth_sl, weighted_layer_value)
    model%quincy_output%ndat_acc = 0._wp
  END SUBROUTINE write_netcdf4_quincy_buffered

  ! ======================================================================================================= !
  !>
  !> accumulate quincy nc4 output acc vars
  !>
  SUBROUTINE accumulate_on_nc_acc_vars(tile, options, output_interval_in_days)
    CLASS(t_jsb_tile_abstract), INTENT(IN) :: tile                     !< one tile with data structure for one lct
    TYPE(t_jsb_task_options),   INTENT(IN) :: options                  !< model options
    REAL(wp),                   INTENT(IN) :: output_interval_in_days  !< interval in days (1.0 == one day exactly)
    ! ----------------------------------------------------------------------------------------------------- !
    TYPE(t_jsb_model),          POINTER :: model
    TYPE(t_quincy_outputfile),  POINTER :: this_q_out_file
    TYPE(t_quincy_out_var),     POINTER :: this_q_out_var
    CLASS(t_jsb_var),           POINTER :: this_var
    REAL(wp)                            :: dtime_one, f_veg_flux_output, this_out_fact
    INTEGER                             :: iblk, ind_var
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':accumulate_on_nc_acc_vars'
    ! ----------------------------------------------------------------------------------------------------- !
    model => Get_model(tile%owner_model_id)
    iblk  = options%iblk

    model%quincy_output%ndat_acc = model%quincy_output%ndat_acc + 1.0_wp

    IF (output_interval_in_days < 1.0_wp) THEN
      dtime_one = 1.0_wp                 ! interval = one timestep
    ELSE
      dtime_one = options%dtime / 1000000.0_wp   ! interval > one timestep
    END IF
    !@TODO: should this still be converted? or what should be done with the fluxes from pools?  
    f_veg_flux_output = dtime_one * 1000000._wp / options%dtime

    ! iterate over the linked list of quincy output files
    this_q_out_file => model%quincy_output%first_output_file
    DO WHILE(ASSOCIATED(this_q_out_file))

      ! Files that do not actually contain to be written variables have not been created and can be skipped
      IF (.NOT. ASSOCIATED(this_q_out_file%first_var%next_var)) THEN
        this_q_out_file => this_q_out_file%next_output_file
        CYCLE
      END IF

      this_q_out_var => this_q_out_file%first_var

      DO WHILE(ASSOCIATED(this_q_out_var))
        IF (this_q_out_var%l_is_acc_var) THEN
          SELECT CASE (this_q_out_var%out_fact_kind)
          CASE (NO_FACT, USE_NDAT_ACC)
            this_out_fact = 1.0_wp
          CASE (USE_VEG_FLUX_FACT)
            this_out_fact = f_veg_flux_output
          CASE DEFAULT
            CALL finish(TRIM(routine), 'Unimplemented output factor kind: ' &
              &         // int2string(this_q_out_var%out_fact_kind))
          END SELECT

          IF (ALLOCATED(this_q_out_var%var_container)) THEN
            DO ind_var = 1,SIZE(this_q_out_var%var_container)
              this_var => this_q_out_var%var_container(ind_var)%p
              SELECT CASE(this_q_out_var%nr_of_dims)
              CASE(1)
                this_q_out_var%accumulated = this_q_out_var%accumulated &
                  & + this_q_out_var%var_container(ind_var)%p%field(ipoi,iblk,1,1,1) * this_out_fact
              !> @TODO **ATTENTION** with newer versions of ICON-Land base modules this might need to be called differently!
              !>   for 2D variables:   this_q_out_var%var_container(ind_var)%p%ptr2d(ics:ice,iblk)
              !>   and for 3D variables: this_q_out_var%var_container(ind_var)%p%ptr3d(?TODO?)
              !> 
              CASE(2)
                this_q_out_var%accumulated_2D(:) = this_q_out_var%accumulated_2D(:) &
                  & + this_q_out_var%var_container(ind_var)%p%field(ipoi,:,iblk,1,1) * this_out_fact
              END SELECT
            ENDDO
          ELSE
            WRITE(message_text,*) 'No pointer to acc variable, please check! Variable name: ' &
              & // TRIM(this_q_out_var%var_name) // ' of file with tag ' // TRIM(this_q_out_file%file_name_tag)
            CALL finish(TRIM(routine),message_text)
          ENDIF
        ENDIF

        this_q_out_var => this_q_out_var%next_var
      ENDDO
      this_q_out_file => this_q_out_file%next_output_file
    ENDDO
  END SUBROUTINE accumulate_on_nc_acc_vars

  ! ======================================================================================================= !
  !>
  !> add current data to buffer
  !>
  SUBROUTINE add_to_nc_buffer(tile, options, nc_buffer_index, elapsed_seconds)
    CLASS(t_jsb_tile_abstract), INTENT(IN) :: tile             !< one tile with data structure for one lct
    TYPE(t_jsb_task_options),   INTENT(IN) :: options          !< model options
    INTEGER,                    INTENT(IN) :: nc_buffer_index  !< where to write in the buffer (current idx)
    REAL(wp),                   INTENT(IN) :: elapsed_seconds  !< time in seconds since start
    ! ----------------------------------------------------------------------------------------------------- !
    TYPE(t_jsb_model),          POINTER :: model
    TYPE(t_quincy_outputfile),  POINTER :: this_q_out_file
    TYPE(t_quincy_out_var),     POINTER :: this_q_out_var
    INTEGER                             :: iblk, ind_var
    REAL(wp)                            :: output_value
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':add_to_nc_buffer'
    ! ----------------------------------------------------------------------------------------------------- !
    model => Get_model(tile%owner_model_id)
    iblk  = options%iblk

    ! iterate over the linked list of quincy output files
    this_q_out_file => model%quincy_output%first_output_file
    DO WHILE(ASSOCIATED(this_q_out_file))

      ! Files that do not actually contain to be written variables have not been created and can be skipped
      IF (.NOT. ASSOCIATED(this_q_out_file%first_var%next_var)) THEN
        this_q_out_file => this_q_out_file%next_output_file
        CYCLE
      END IF

      ! The first variable is time
      this_q_out_var => this_q_out_file%first_var
      ! For time we just write the elapsed sconds into the buffer
      this_q_out_var%buffer_2D(nc_buffer_index) = elapsed_seconds
      
      ! We skip the first variable which is time
      this_q_out_var => this_q_out_var%next_var

      DO WHILE(ASSOCIATED(this_q_out_var))
        SELECT CASE (this_q_out_var%nr_of_dims)
        CASE (1)
          ! If it is a FLUX variable use the accumulated values
          IF(this_q_out_var%l_is_acc_var) THEN 
            IF ((this_q_out_var%out_fact_kind .EQ. USE_NDAT_ACC) &
                & .AND. (model%quincy_output%ndat_acc .GT. 1.0_wp)) THEN
              this_q_out_var%accumulated  = this_q_out_var%accumulated / model%quincy_output%ndat_acc
            ENDIF
            this_q_out_var%buffer_2D(nc_buffer_index) = this_q_out_var%accumulated

          ! If it is a POOL variable just update the buffer with the current status of the pool
          ELSE
            output_value = 0._wp
            IF (ALLOCATED(this_q_out_var%var_container)) THEN
              DO ind_var = 1,SIZE(this_q_out_var%var_container)
                output_value = output_value &
                  & + this_q_out_var%var_container(ind_var)%p%field(ipoi,iblk,1,1,1)
                !> TODO **ATTENTION** with newer versions of ICON-Land base modules this might need to be called differently!
                ! for 2D variables:
                !  & + this_q_out_var%var_container(ind_var)%p%ptr2d(ics:ice,iblk)
              ENDDO
            ELSE
              WRITE(message_text,*) 'No pointer to acc variable, please check! Variable name: ' &
                & // TRIM(this_q_out_var%var_name) // ' of file with tag ' // TRIM(this_q_out_file%file_name_tag)
              CALL finish(TRIM(routine),message_text)
            ENDIF
            this_q_out_var%buffer_2D(nc_buffer_index) = output_value
          ENDIF
          
        ! Two dimensional output data
        CASE (2)
          IF (this_q_out_var%l_is_acc_var) THEN
            IF (this_q_out_var%var_type .EQ. SOIL_LAYERED) THEN
              IF ((this_q_out_var%out_fact_kind .EQ. USE_NDAT_ACC) .AND. (model%quincy_output%ndat_acc .GT. 1.0_wp)) THEN
                this_q_out_var%buffer_3D(:, nc_buffer_index) = this_q_out_var%accumulated_2D(:) / model%quincy_output%ndat_acc
              ELSE
                this_q_out_var%buffer_3D(:, nc_buffer_index) = this_q_out_var%accumulated_2D(:)
              ENDIF
            ELSE
              CALL finish(TRIM(routine), 'Accumulating only for soil layered vars should have been caught already before!')
            ENDIF
          ELSE
            !TODO -> is %field(ipoi,:,iblk,1,1) still valid if recent changes in ICON-Land memory structure are merged back?
                  ! ? this_q_out_var%var_container(ind_var)%p%ptr3d(?)
            IF (this_q_out_var%var_type .EQ. SOIL_LAYERED) THEN
              this_q_out_var%buffer_3D(:, nc_buffer_index) = this_q_out_var%var_container(1)%p%field(ipoi,:,iblk,1,1)
            ELSEIF (this_q_out_var%var_type .EQ. CANOPY_LAYERED) THEN
              this_q_out_var%buffer_3D(:, nc_buffer_index) = this_q_out_var%var_container(1)%p%field(ipoi,:,iblk,1,1)
            ELSE
              CALL finish(TRIM(routine), 'Neither soil layered nor canopy layered should have been caught already before!')
            ENDIF
          ENDIF

        END SELECT
        this_q_out_var => this_q_out_var%next_var
      ENDDO
      this_q_out_file => this_q_out_file%next_output_file
    ENDDO
  END SUBROUTINE add_to_nc_buffer

  ! ======================================================================================================= !
  !>
  !> close nc4 output files
  !>
  SUBROUTINE close_nc4_output_files(model_id)
    INTEGER, INTENT(IN) :: model_id !< ID of model instance
    ! ----------------------------------------------------------------------------------------------------- !
    TYPE(t_jsb_model),         POINTER  :: model
    TYPE(t_quincy_outputfile), POINTER  :: this_q_out_file
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':close_nc4_output_files'
    ! ----------------------------------------------------------------------------------------------------- !
    model => Get_model(model_id)

#ifdef __QUINCY_WITH_NETCDF__
    ! iterate over the linked list of quincy output files
    this_q_out_file => model%quincy_output%first_output_file
    DO WHILE(ASSOCIATED(this_q_out_file))

      ! Files that do not actually contain to be written variables have not been created and can be skipped
      IF (.NOT. ASSOCIATED(this_q_out_file%first_var%next_var)) THEN
        this_q_out_file => this_q_out_file%next_output_file
        CYCLE
      END IF

      ! Close the current file. This frees up any internal netCDF resources
      ! associated with the file, and flushes any buffers.
      CALL check_nc4_command_stat(nf90_close(this_q_out_file%file_id), routine//'close:file')
      this_q_out_file%file_id = -1

      this_q_out_file => this_q_out_file%next_output_file
    ENDDO
#endif
  END SUBROUTINE close_nc4_output_files

  ! ======================================================================================================= !
  !>
  !> deallocate the content of the quincy output structure
  !>
  SUBROUTINE deallocate_quincy_output(model_id)
    INTEGER, INTENT(IN) :: model_id !< ID of model instance
    ! ----------------------------------------------------------------------------------------------------- !
    TYPE(t_jsb_model),         POINTER  :: model
    TYPE(t_quincy_outputfile), POINTER  :: this_q_out_file
    TYPE(t_quincy_out_var),    POINTER  :: this_q_var
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':deallocate_quincy_output'
    ! ----------------------------------------------------------------------------------------------------- !
    model => Get_model(model_id)

    ! iterate backwards over the linked list of quincy output files
    this_q_out_file => model%quincy_output%current_output_file
    DO WHILE(ASSOCIATED(this_q_out_file))

      ! same for the variables
      this_q_var => this_q_out_file%current_var
      DO WHILE(ASSOCIATED(this_q_var))

        ! Deallocate the stuff on this variable
        IF (ALLOCATED(this_q_var%var_container)) DEALLOCATE(this_q_var%var_container)
        IF (ALLOCATED(this_q_var%accumulated_2D)) DEALLOCATE(this_q_var%accumulated_2D)

        IF (ALLOCATED(this_q_var%buffer_2D)) DEALLOCATE(this_q_var%buffer_2D)
        IF (ALLOCATED(this_q_var%buffer_3D)) DEALLOCATE(this_q_var%buffer_3D)

        IF (ASSOCIATED(this_q_var%prev_var)) THEN
          ! Go to previous variable in the linked list
          this_q_var => this_q_var%prev_var
          ! and deallocate the previous last variable
          DEALLOCATE(this_q_var%next_var)
        ELSE
          DEALLOCATE(this_q_var)
        ENDIF
      ENDDO

      IF (ASSOCIATED(this_q_out_file%prev_output_file)) THEN
        ! Go to previous file in the linked list
        this_q_out_file => this_q_out_file%prev_output_file
        ! and deallocate the already closed current last file can be deallocated 
        DEALLOCATE(this_q_out_file%next_output_file)
      ELSE
        DEALLOCATE(this_q_out_file)
      ENDIF
    ENDDO

    IF (ASSOCIATED(model%quincy_output%sel_output_variables)) DEALLOCATE(model%quincy_output%sel_output_variables)

    DEALLOCATE(model%quincy_output)
  END SUBROUTINE deallocate_quincy_output

  ! ======================================================================================================= !
  !>
  !> write 2d var to nc4 output file
  !>
  SUBROUTINE write_array_to_1d_var_on_file(output_values, var_id, file_id, current_record, buffer_size)
    REAL(wp), INTENT(IN) :: output_values(:)   !< output value to write
    INTEGER,  INTENT(IN) :: var_id            !< ID of var to write to
    INTEGER,  INTENT(IN) :: file_id           !< ID of file to write to
    INTEGER,  INTENT(IN) :: current_record    !< where to write in the record dimension (here: time)
    INTEGER,  INTENT(IN) :: buffer_size       !< number of soil layers
    ! ----------------------------------------------------------------------------------------------------- !
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':write_array_to_1d_var_on_file'
    ! ----------------------------------------------------------------------------------------------------- !
#ifdef __QUINCY_WITH_NETCDF__
    CALL check_nc4_command_stat(nf90_put_var(file_id, var_id, output_values, start=(/current_record /), &
      &                         count=(/ buffer_size /)), routine//'put_var::id'//int2string(var_id))
#endif
  END SUBROUTINE write_array_to_1d_var_on_file

  ! ======================================================================================================= !
  !>
  !> write 3d var to nc4 output file
  !>
  SUBROUTINE write_2d_array_to_2d_var_on_file(output_values, var_id, file_id, current_record, buffer_size, dim2_len)
    REAL(wp), INTENT(IN) :: output_values(:, :)   !< output value to write
    INTEGER,  INTENT(IN) :: var_id            !< ID of var to write to
    INTEGER,  INTENT(IN) :: file_id           !< ID of file to write to
    INTEGER,  INTENT(IN) :: current_record    !< where to write in the record dimension (here: time)
    INTEGER,  INTENT(IN) :: buffer_size       !< number timesteps per buffer
    INTEGER,  INTENT(IN) :: dim2_len          !< Length of second dimension
    ! ----------------------------------------------------------------------------------------------------- !
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':write_2d_array_to_2d_var_on_file'
    ! ----------------------------------------------------------------------------------------------------- !
#ifdef __QUINCY_WITH_NETCDF__
    CALL check_nc4_command_stat(nf90_put_var(file_id, var_id, output_values, start=(/1, current_record/), &
      &                         count=(/dim2_len, buffer_size/)), routine//'put_var::id'//int2string(var_id))
#endif
  END SUBROUTINE write_2d_array_to_2d_var_on_file

  ! ======================================================================================================= !
  !>
  !> check status of nc4 command
  !>
  SUBROUTINE check_nc4_command_stat(status, routine_and_context)
    INTEGER,          INTENT(IN) :: status
    CHARACTER(len=*), INTENT(IN) :: routine_and_context
    ! ----------------------------------------------------------------------------------------------------- !
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':check_nc4_command_stat'
    ! ----------------------------------------------------------------------------------------------------- !
#ifdef __QUINCY_WITH_NETCDF__
    IF (status /= nf90_noerr) THEN
      CALL finish(TRIM(routine_and_context), TRIM(nf90_strerror(status)))
    ENDIF
#endif
  END SUBROUTINE check_nc4_command_stat

  ! ======================================================================================================= !
  !>
  !> Calculated the weighted soil depth sum of the passed layered variable
  !>
  SUBROUTINE calculate_soil_depth_weighted_layer_sum(layered_array, soil_depth_sl, weighted_sum)
    REAL(wp), INTENT(IN)  :: layered_array(:,:)
    REAL(wp), INTENT(IN)  :: soil_depth_sl(:)
    REAL(wp), INTENT(OUT) :: weighted_sum(:)
    ! ----------------------------------------------------------------------------------------------------- !
    INTEGER :: idx
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calculate_soil_depth_weighted_layer_sum'
    ! ----------------------------------------------------------------------------------------------------- !
    DO idx = 1,SIZE(weighted_sum)
      weighted_sum(idx) = SUM(layered_array(:,idx) * soil_depth_sl(:))
    ENDDO
  END SUBROUTINE calculate_soil_depth_weighted_layer_sum

  ! ======================================================================================================= !
  !>
  !> Calculated the weighted soil depth avg of the passed layered variable
  !>
  SUBROUTINE calculate_soil_depth_weighted_layer_avg(layered_array, soil_depth_sl, weighted_avg)
    REAL(wp), INTENT(IN)  :: layered_array(:,:)
    REAL(wp), INTENT(IN)  :: soil_depth_sl(:)
    REAL(wp), INTENT(OUT) :: weighted_avg(:)
    ! ----------------------------------------------------------------------------------------------------- !
    INTEGER  :: idx
    REAL(wp) :: soil_depth
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calculate_soil_depth_weighted_layer_avg'
    ! ----------------------------------------------------------------------------------------------------- !
    soil_depth = SUM(soil_depth_sl(:))
  
    DO idx = 1,SIZE(weighted_avg)
      weighted_avg(idx) = SUM(layered_array(:,idx) * soil_depth_sl(:)) / soil_depth
    ENDDO
  END SUBROUTINE calculate_soil_depth_weighted_layer_avg

  ! ======================================================================================================= !
  !>
  !> Collects pointers for a certain element from the compartments of a given veg bgcm
  !>
  SUBROUTINE collect_veg_bgcm_pointers_for_element(this_elem_idx, veg_bgcm, compartment_names, compartment_pointers)
    INTEGER,                                  INTENT(IN)    :: this_elem_idx           !< Index of the element to add
    TYPE(t_veg_bgcm_with_components), TARGET, INTENT(IN)    :: veg_bgcm                !< pointer to the veg bgcm
    CHARACTER(len=100),                       INTENT(INOUT) :: compartment_names(:)    !< vector with the names of the compartments
    TYPE(t_jsb_var_p),                TARGET, INTENT(INOUT) :: compartment_pointers(:) !< vector with the pointers to the queried element on the compartments
    ! ----------------------------------------------------------------------------------------------------- !
    INTEGER :: ind, nr_of_compartments, i_compartment
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':collect_veg_bgcm_pointers_for_element'
    ! ----------------------------------------------------------------------------------------------------- !
    ind = 1
    nr_of_compartments = SIZE(veg_bgcm%pool_list)

    DO i_compartment = 1, nr_of_compartments
      compartment_names(ind) = veg_bgcm%pool_list(i_compartment)%p%shortname
      compartment_pointers(ind)%p => veg_bgcm%pool_list(i_compartment)%p%element_list(this_elem_idx)%p
      ind = ind + 1
    END DO ! compartments

  END SUBROUTINE collect_veg_bgcm_pointers_for_element


  ! ======================================================================================================= !
  !>
  !> Collects pointers for a certain element from the compartments of a given sb bgcm
  !>
  SUBROUTINE collect_sb_bgcm_pointers_for_element(this_elem_idx, quincy_sb_bgcm_struct, compartment_names, compartment_pointers)
    INTEGER,                                 INTENT(IN)    :: this_elem_idx           !< index of the element to add to the output
    TYPE(t_sb_bgcm_with_components), TARGET, INTENT(IN)    :: quincy_sb_bgcm_struct   !< pointer to the sb bgcm
    CHARACTER(len=100),                      INTENT(INOUT) :: compartment_names(:)    !< vector with the names of all compartments
    TYPE(t_jsb_var_p),               TARGET, INTENT(INOUT) :: compartment_pointers(:) !< vector with the pointers to the queried element on the compartments
    ! ----------------------------------------------------------------------------------------------------- !
    INTEGER :: ind
    CHARACTER(len=*), PARAMETER   :: routine = TRIM(modname)//':collect_sb_bgcm_pointers_for_element'
    ! ----------------------------------------------------------------------------------------------------- !
    ! collect compartment names and pointer
    DO ind = 1,SIZE(quincy_sb_bgcm_struct%pool_list)
      compartment_names(ind) = quincy_sb_bgcm_struct%pool_list(ind)%p%shortname
      compartment_pointers(ind)%p => quincy_sb_bgcm_struct%pool_list(ind)%p%element_list(this_elem_idx)%p
    ENDDO
  END SUBROUTINE collect_sb_bgcm_pointers_for_element


  ! ======================================================================================================= !
  !>
  !> checks that netcdf output is possible with the given output interval 
  !>
  SUBROUTINE assert_possible_output_interval(output_interval_in_days, output_interval_char, caller, tag)
    REAL(wp),                   INTENT(IN) :: output_interval_in_days !< interval in days (1.0 == one day exactly)
    CHARACTER(len=*),           INTENT(IN) :: output_interval_char    !< output interval tag for filename
    CHARACTER(len=*),           INTENT(IN) :: caller                  !< calling routine
    CHARACTER(len=*), OPTIONAL, INTENT(IN) :: tag                     !< providing more context, e.g. spin-up
    ! ----------------------------------------------------------------------------------------------------- !
    LOGICAL                       :: is_expected_output_interval
    CHARACTER(len=:), ALLOCATABLE :: this_tag
    ! ----------------------------------------------------------------------------------------------------- !

    ! with netcdf all 'regular' timesteps are possible, i.e. yearly, weekly, daily and timestep output
    ! monthly output however is not possible
    is_expected_output_interval = (output_interval_in_days == one_year) &
      & .OR. (output_interval_in_days == 7.0_wp)                        &
      & .OR. ((output_interval_in_days <= 1.0_wp) .AND. (output_interval_in_days > 0.0_wp))

    IF(PRESENT(tag)) THEN
      this_tag = tag
    ELSE
      this_tag = ""
    ENDIF

    IF(.NOT. is_expected_output_interval) THEN
      WRITE(message_text,*) 'Unexpected '//TRIM(this_tag)//' output interval'                   &
        & //' - For nc these are currently yearly, weekly, daily or timestep. '                 &
        & //' However, interval was ' // real2string(output_interval_in_days) // ' with tag '   &
        & //TRIM(output_interval_char) // '! Please check!'
      CALL finish(TRIM(caller),message_text)
    ENDIF

  END SUBROUTINE assert_possible_output_interval

END MODULE mo_qs_output_netcdf4
