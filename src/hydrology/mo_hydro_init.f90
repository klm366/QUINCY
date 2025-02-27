!> Initialization of the the hydrology process
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
!>#### Initialize soil properties and initial conditions for the hydrology process
!>
!> Soil properties are set in [[hydro_init_bc]] by reading external parameters from the file given by the
!> [[t_hydro_config:bc_filename]] namelist parameter.
!> The properties are either read directly from the corresponding fields or, with
!> [[t_hydro_config:l_soil_texture]]=true, soil textures for sand, silt (loam), clay and organic
!> matter (peat) are read from the input file and used to derive the properties from a lookup table.
!>
!> Standard deviation of topography (and elevation for the standalone model) is read from the file given
!> by [[t_hydro_config:bc_sso_filename]].
!>
!> Initial conditions for soil and snow water are read in [[hydro_init_ic]] from the file given by
!> [[t_hydro_config:ic_filename]] if [[t_jsb_model_config:init_from_ifs]] == .FALSE. otherwise they are
!> derived from IFS initial data read from the file given by [[t_jsb_model_config:ifs_input_file]].
!>
!> In order to avoid reading the input data for each tile again they are read only once for the first tile
!> in [[hydro_read_init_vars]] and stored in variable [[mo_hydro_init:hydro_init_vars]].
!>
MODULE mo_hydro_init
#ifndef __NO_JSBACH__

  USE mo_kind,               ONLY: wp
  USE mo_exception,          ONLY: message, message_text, warning, finish
  USE mo_jsb_control,        ONLY: debug_on, jsbach_runs_standalone

  USE mo_jsb_model_class,    ONLY: t_jsb_model
  USE mo_jsb_grid_class,     ONLY: t_jsb_grid, t_jsb_vgrid
  USE mo_jsb_grid,           ONLY: Get_grid, Get_vgrid
  USE mo_jsb_tile_class,     ONLY: t_jsb_tile_abstract
  USE mo_jsb_class,          ONLY: get_model
  USE mo_jsb_io_netcdf,      ONLY: t_input_file, jsb_netcdf_open_input
  USE mo_jsb_io,             ONLY: missval
  USE mo_jsb_impl_constants, ONLY: ifs_nsoil, ifs_soil_depth
  USE mo_util,               ONLY: int2string, soil_init_from_texture

  dsl4jsb_Use_processes HYDRO_, PHENO_
  dsl4jsb_Use_config(HYDRO_)
  dsl4jsb_Use_memory(HYDRO_)
  dsl4jsb_Use_memory(PHENO_)

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: hydro_init, hydro_sanitize_state

  ! -------------------------------------------------------------------------------------------------- !
  ! Lookup table for deriving soil properties from soil textures

  ! TODO move the below parameters to constants / parameter module

  ! Parameters for sand, loam (here called silt), clay, peat (here called oc) are used from the TERRA
  ! land model in ICON (see src/lnd_phy_schemes/sfc_terra_data.f90 in ICON)

  ! Pore volume (called cporv in TERRA)
  REAL(wp), PARAMETER ::  porosity_sand = 0.364_wp      !< pore volume of sand [vol. fraction]
  REAL(wp), PARAMETER ::  porosity_silt = 0.455_wp      !< pore volume of loam [vol. fraction]
  REAL(wp), PARAMETER ::  porosity_clay = 0.507_wp      !< pore volume of clay [vol. fraction]
  REAL(wp), PARAMETER ::  porosity_oc   = 0.863_wp      !< pore volume of peat [vol. fraction]

  ! Field capacity (called cfcap in TERRA)
  REAL(wp), PARAMETER ::  field_cap_sand = 0.196_wp     !< field capacity of sand [vol. fraction]
  REAL(wp), PARAMETER ::  field_cap_silt = 0.340_wp     !< field capacity of loam [vol. fraction]
  REAL(wp), PARAMETER ::  field_cap_clay = 0.463_wp     !< field capacity of clay [vol. fraction]
  REAL(wp), PARAMETER ::  field_cap_oc   = 0.763_wp     !< field capacity of peat [vol. fraction]

  ! Hydraulic conductivity at saturation (called ckw0 in TERRA)
  REAL(wp), PARAMETER ::  hyd_cond_sand = 4.79e-5_wp    !< hydraulic conductivity at saturation of sand [m/s]
  REAL(wp), PARAMETER ::  hyd_cond_silt = 5.31e-6_wp    !< hydraulic conductivity at saturation of loam [m/s]
  REAL(wp), PARAMETER ::  hyd_cond_clay = 8.50e-8_wp    !< hydraulic conductivity at saturation of clay [m/s]
  REAL(wp), PARAMETER ::  hyd_cond_oc   = 5.80e-8_wp    !< hydraulic conductivity at saturation of peat [m/s]

  ! Wilting point (called cpwp in TERRA)
  REAL(wp), PARAMETER ::  wilt_sand = 0.042_wp          !< wilting point of sand [vol. fraction]
  REAL(wp), PARAMETER ::  wilt_silt = 0.11_wp           !< wilting point of loam [vol. fraction]
  REAL(wp), PARAMETER ::  wilt_clay = 0.257_wp          !< wilting point of clay [vol. fraction]
  REAL(wp), PARAMETER ::  wilt_oc   = 0.265_wp          !< wilting point of peat [vol. fraction]

  ! Pore size index
  REAL(wp), PARAMETER :: pore_size_index_sand = 0.35_wp !< pore size index of sand []
  REAL(wp), PARAMETER :: pore_size_index_silt = 0.2_wp  !< pore size index of loam []
  REAL(wp), PARAMETER :: pore_size_index_clay = 0.13_wp !< pore size index of clay []
  REAL(wp), PARAMETER :: pore_size_index_oc   = 0.65_wp !< pore size index of loam []

  ! Clapp & Hornberger exponent b (see Beringer et al. 2001)
  REAL(wp), PARAMETER ::  bclapp_sand = 3.39_wp   !< Clapp & Hornberger exponent b for sand (type 1 in Beringer et al. 2001)
  REAL(wp), PARAMETER ::  bclapp_silt = 4.98_wp   !< Clapp & Hornberger exponent b for loam (type 5 in Beringer et al. 2001)
  REAL(wp), PARAMETER ::  bclapp_clay = 10.38_wp  !< Clapp & Hornberger exponent b for clay (type 10 in Beringer et al. 2001)
  REAL(wp), PARAMETER ::  bclapp_oc   = 4._wp     !< Clapp & Hornberger exponent b for peat (type 12 in Beringer et al. 2001)

  ! Soil matric potential (see Beringer et al. 2001)
  REAL(wp), PARAMETER ::  matric_pot_sand = -0.04729_wp !< soil matric potential of sand (type 1 in Beringer at al. 2001)
  REAL(wp), PARAMETER ::  matric_pot_silt = -0.45425_wp !< soil matric potential of loam (type 5 in Beringer at al. 2001)
  REAL(wp), PARAMETER ::  matric_pot_clay = -0.633_wp   !< soil matric potential of clay (type 10 in Beringer at al. 2001)
  REAL(wp), PARAMETER ::  matric_pot_oc = -0.12_wp      !< soil matric potential of peat (type 12 in Beringer at al. 2001)

  REAL(wp), PARAMETER ::  hyd_cond_sat_profile = 0.432332_wp
                          !< factor to compensate profile of hyd_cond_sat with depth in TERRA

  !> Residual soil moisture (fraction of volume; Maidment, Handbook of Hydrology, 1993)
  REAL(wp), PARAMETER ::  wres_sand = 0.020_wp          !< residual soil moisture (type sand) [vol. fraction]
  REAL(wp), PARAMETER ::  wres_silt = 0.015_wp          !< residual soil moisture (type silt loam) [vol. fraction]
  REAL(wp), PARAMETER ::  wres_clay = 0.090_wp          !< residual soil moisture (type clay) [vol. fraction]
  REAL(wp), PARAMETER ::  wres_oc   = 0.150_wp          !< residual soil moisture (Letts et al., 2000)

  !
  ! -------------------------------------------------------------------------------------------------- !
  !
  !> Type to hold variables read from input files
  TYPE t_hydro_init_vars
    REAL(wp), POINTER ::                  &
      & elevation        (:,:) => NULL(), & !< geometric height of surface above sea level
      & oro_stddev       (:,:) => NULL(), & !< standard deviation of subgrid-scale orography
      & soil_depth       (:,:) => NULL(), & !< soil depth until bedrock
      & vol_porosity     (:,:) => NULL(), & !< volumetric soil porosity
      & hyd_cond_sat     (:,:) => NULL(), & !< saturated hydraulic conductivity
      & matric_pot       (:,:) => NULL(), & !< saturated matric potential
      & bclapp           (:,:) => NULL(), & !< Clapp and Hornberger exponent b
      & pore_size_index  (:,:) => NULL(), & !< pore size distribution index
      & vol_field_cap    (:,:) => NULL(), & !< volumetric soil field capacity
      & vol_p_wilt       (:,:) => NULL(), & !< volumetric wilting point
      & vol_wres         (:,:) => NULL(), & !< volumetric residual water content
      & root_depth       (:,:) => NULL(), & !< root depth
      & weq_snow_soil    (:,:) => NULL(), & !< initial water equivalent of snow on surface
      & fract_pond_max   (:,:) => NULL(), & !< maximum pond fraction
      & depth_pond_max   (:,:) => NULL(), & !< maximum pond depth
      & wtr_soil_sl    (:,:,:) => NULL(), & !< initial soil moisture within soil layers
      & fract_org_sl   (:,:,:) => NULL(), & !< fraction of soil organic matter within soil layers
      & fr_sand          (:,:) => NULL(), & !< fraction of sand for soil index
      & fr_silt          (:,:) => NULL(), & !< fraction of silt for soil index
      & fr_clay          (:,:) => NULL(), & !< fraction of clay for soil index
      & fr_sand_deep     (:,:) => NULL(), & !< fraction of sand for deep soil index
      & fr_silt_deep     (:,:) => NULL(), & !< fraction of silt for deep soil index
      & fr_clay_deep     (:,:) => NULL(), & !< fraction of clay for deep soil index
      & ifs_smi_sl     (:,:,:) => NULL(), & !< volumetric soil water (soil moisture index) from IFS
      & ifs_weq_snow     (:,:) => NULL()    !< water equivalent of snow on surface from IFS
  END TYPE t_hydro_init_vars

  TYPE(t_hydro_init_vars) :: hydro_init_vars                !< Module variable holding input fields

  CHARACTER(len=*), PARAMETER :: modname = 'mo_hydro_init'  !< Name of this module

CONTAINS

  ! ====================================================================================================== !
  !>
  !>#### Intialize the hydrology process
  !>
  SUBROUTINE hydro_init(tile)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile               !< current tile

    CHARACTER(len=*), PARAMETER :: routine = modname//':hydro_init' !< name of this routine

    !> - Read hydrology initial and boundary condition data <br>
    !>   The initial data is read only once - by the top tile - and is made avilable to the other tiles.
    !>
    ! Note: this requires that the init routines in mo_jsb_model_init:jsbach_init are called starting
    !       from the top tile!
    IF (.NOT. ASSOCIATED(tile%parent_tile)) THEN  ! This is the top tile.
      CALL hydro_read_init_vars(tile)
    END IF

    !> - Set soil properties (on all tiles)
    CALL hydro_init_bc(tile)

    !> - Set soil initial conditions (on all tiles)
    CALL hydro_init_ic(tile)

    !> - Finalize hydrology initialization after the last tile has been initialized
    IF (tile%Is_last_process_tile(HYDRO_)) THEN
      CALL hydro_finalize_init_vars()
    END IF

  END SUBROUTINE hydro_init

  ! ====================================================================================================== !
  !>
  !>#### Finalize hydrology initialization
  !>
  SUBROUTINE hydro_finalize_init_vars

    !> Deallocation of the variables that had only been needed during the initialization phase

    DEALLOCATE( &
      & hydro_init_vars%oro_stddev      ,       &
      & hydro_init_vars%soil_depth      ,       &
      & hydro_init_vars%vol_porosity    ,       &
      & hydro_init_vars%hyd_cond_sat    ,       &
      & hydro_init_vars%matric_pot      ,       &
      & hydro_init_vars%bclapp          ,       &
      & hydro_init_vars%pore_size_index ,       &
      & hydro_init_vars%vol_field_cap   ,       &
      & hydro_init_vars%vol_p_wilt      ,       &
      & hydro_init_vars%vol_wres        ,       &
      & hydro_init_vars%root_depth      ,       &
      & hydro_init_vars%fract_org_sl    ,       &
      & hydro_init_vars%fract_pond_max  ,       &
      & hydro_init_vars%depth_pond_max  ,       &
      & hydro_init_vars%weq_snow_soil)

    IF (ASSOCIATED(hydro_init_vars%elevation))    DEALLOCATE(hydro_init_vars%elevation)
    IF (ASSOCIATED(hydro_init_vars%fr_sand))      DEALLOCATE(hydro_init_vars%fr_sand)
    IF (ASSOCIATED(hydro_init_vars%fr_silt))      DEALLOCATE(hydro_init_vars%fr_silt)
    IF (ASSOCIATED(hydro_init_vars%fr_clay))      DEALLOCATE(hydro_init_vars%fr_clay)
    IF (ASSOCIATED(hydro_init_vars%fr_sand_deep)) DEALLOCATE(hydro_init_vars%fr_sand_deep)
    IF (ASSOCIATED(hydro_init_vars%fr_silt_deep)) DEALLOCATE(hydro_init_vars%fr_silt_deep)
    IF (ASSOCIATED(hydro_init_vars%fr_clay_deep)) DEALLOCATE(hydro_init_vars%fr_clay_deep)
    IF (ASSOCIATED(hydro_init_vars%wtr_soil_sl))  DEALLOCATE(hydro_init_vars%wtr_soil_sl)
    IF (ASSOCIATED(hydro_init_vars%ifs_smi_sl))   DEALLOCATE(hydro_init_vars%ifs_smi_sl)
    IF (ASSOCIATED(hydro_init_vars%ifs_weq_snow)) DEALLOCATE(hydro_init_vars%ifs_weq_snow)

  END SUBROUTINE hydro_finalize_init_vars

  ! ====================================================================================================== !
  !
  !>#### Read variables from initial files
  !>
  !> Soil parameters and initial conditions for the hydrology need to be set on all tiles on which
  !> the HYDRO_ process is active. Currently, the data is identical for all tiles. We thus only read it
  !> once by the top tile and store it in variable [[mo_hydro_init:hydro_init_vars]] of type
  !> [[t_hydro_init_vars]]. The remaining tiles are then using this variable to initialize their memories.
  !>
  SUBROUTINE hydro_read_init_vars(tile)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile      !< Current tile

    dsl4jsb_Def_config(HYDRO_)                  !< Configurable parameters of the HYDRO_ process

    REAL(wp), POINTER :: &
      & ptr_2D(:,  :),   &                      !< Pointer to 2d variable
      & ptr_3D(:,:,:)                           !< Pointer to 3d variable
    LOGICAL,  POINTER   :: is_in_domain(:,:)    !< True: the cell belongs to domain; False: halo cell

    TYPE(t_jsb_model), POINTER :: model         !< Current instance of the model
    TYPE(t_jsb_grid),  POINTER :: grid          !< Horizontal grid
    TYPE(t_jsb_vgrid), POINTER :: soil_w        !< Vertical grid used for hydrology
    TYPE(t_input_file)         :: input_file    !< File containing input data to be read
    INTEGER                    :: nproma, nblks !< Dimensions of the horizontal grid
    INTEGER                    :: nsoil         !< Number of vertical soil layers

    INTEGER :: isoil                            !< Looping index for soil layers

    CHARACTER(len=*), PARAMETER :: routine = modname//':hydro_read_init_vars'

    ! Assertion: do nothing if hydrology process is not active on this tile (should not happen)
    IF (.NOT. tile%Is_process_active(HYDRO_)) RETURN

    IF (debug_on()) CALL message(routine, 'Reading/setting hydrology init vars')

    model => Get_model(tile%owner_model_id)

    dsl4jsb_Get_config(HYDRO_)

    grid   => Get_grid(model%grid_id)
    nproma = grid%Get_nproma()
    nblks  = grid%Get_nblks()
    soil_w => Get_vgrid('soil_depth_water')
    nsoil = soil_w%n_levels
    is_in_domain => grid%patch%cells%decomp_info%owner_mask(:,:)

    ! -------------------------------------------------------------------------------------------------- !
    !>
    !> - Allocate memory for the input variables
    !>
    ALLOCATE( &
      & hydro_init_vars%oro_stddev      (nproma, nblks),       &
      & hydro_init_vars%soil_depth      (nproma, nblks),       &
      & hydro_init_vars%vol_porosity    (nproma, nblks),       &
      & hydro_init_vars%hyd_cond_sat    (nproma, nblks),       &
      & hydro_init_vars%matric_pot      (nproma, nblks),       &
      & hydro_init_vars%bclapp          (nproma, nblks),       &
      & hydro_init_vars%pore_size_index (nproma, nblks),       &
      & hydro_init_vars%vol_field_cap   (nproma, nblks),       &
      & hydro_init_vars%vol_p_wilt      (nproma, nblks),       &
      & hydro_init_vars%vol_wres        (nproma, nblks),       &
      & hydro_init_vars%root_depth      (nproma, nblks),       &
      & hydro_init_vars%fract_pond_max  (nproma, nblks),       &
      & hydro_init_vars%depth_pond_max  (nproma, nblks),       &
      & hydro_init_vars%weq_snow_soil   (nproma, nblks)        & ! TODO only needed for init_from_ifs == .FALSE.
      & )

      ! Initialize input variables with missing value
      hydro_init_vars%oro_stddev     (:,:)   = missval
      hydro_init_vars%soil_depth     (:,:)   = missval
      hydro_init_vars%vol_porosity   (:,:)   = missval
      hydro_init_vars%hyd_cond_sat   (:,:)   = missval
      hydro_init_vars%matric_pot     (:,:)   = missval
      hydro_init_vars%bclapp         (:,:)   = missval
      hydro_init_vars%pore_size_index(:,:)   = missval
      hydro_init_vars%vol_field_cap  (:,:)   = missval
      hydro_init_vars%vol_p_wilt     (:,:)   = missval
      hydro_init_vars%vol_wres       (:,:)   = missval
      hydro_init_vars%root_depth     (:,:)   = missval
      hydro_init_vars%fract_pond_max (:,:)   = missval
      hydro_init_vars%depth_pond_max (:,:)   = missval
      hydro_init_vars%weq_snow_soil  (:,:)   = missval     ! TODO only needed for init_from_ifs == .FALSE.

    !
    ! More of that
    !
    ! TODO only necessary for l_socmap == .TRUE.
    ALLOCATE(hydro_init_vars%fract_org_sl(nproma, nsoil, nblks))
    hydro_init_vars%fract_org_sl(:,:,:) = missval

    IF (jsbach_runs_standalone()) THEN
       ALLOCATE(hydro_init_vars%elevation(nproma,nblks))
       hydro_init_vars%elevation(:,:) = missval
    END IF
    IF (model%config%init_from_ifs) THEN
      ALLOCATE(hydro_init_vars%ifs_weq_snow(nproma,nblks))
      hydro_init_vars%ifs_weq_snow(:,:) = missval
    END IF
    IF (model%config%init_from_ifs .AND. .NOT. dsl4jsb_Config(HYDRO_)%l_read_initial_moist) THEN
      ALLOCATE(hydro_init_vars%ifs_smi_sl(nproma,ifs_nsoil,nblks))
      hydro_init_vars%ifs_smi_sl(:,:,:) = missval
    ELSE
      ALLOCATE(hydro_init_vars%wtr_soil_sl(nproma, nsoil, nblks))
      hydro_init_vars%wtr_soil_sl(:,:,:) = missval
    END IF
    IF (dsl4jsb_Config(HYDRO_)%l_soil_texture) THEN
      ALLOCATE( &
        & hydro_init_vars%fr_sand     (nproma, nblks), &
        & hydro_init_vars%fr_silt     (nproma, nblks), &
        & hydro_init_vars%fr_clay     (nproma, nblks), &
        & hydro_init_vars%fr_sand_deep(nproma, nblks), &
        & hydro_init_vars%fr_silt_deep(nproma, nblks), &
        & hydro_init_vars%fr_clay_deep(nproma, nblks)  &
        & )
      hydro_init_vars%fr_sand     (:,:) = missval
      hydro_init_vars%fr_silt     (:,:) = missval
      hydro_init_vars%fr_clay     (:,:) = missval
      hydro_init_vars%fr_sand_deep(:,:) = missval
      hydro_init_vars%fr_silt_deep(:,:) = missval
      hydro_init_vars%fr_clay_deep(:,:) = missval
    END IF ! l_soil_texture

    ! -------------------------------------------------------------------------------------------------- !
    !
    !> - Read soil properties from input files `bc_sso_filename` and `bc_filename`
    !

    ! Open input file - bc_sso
    input_file = jsb_netcdf_open_input(TRIM(dsl4jsb_Config(HYDRO_)%bc_sso_filename), model%grid_id)

    ! Orography ...
    IF (jsbach_runs_standalone()) THEN
       ptr_2D => input_file%Read_2d(  &
         & variable_name='elevation', &
         & fill_array = hydro_init_vars%elevation)
       hydro_init_vars%elevation = MERGE(ptr_2D, 0._wp, ptr_2D >= 0._wp)
    END IF

    ptr_2D => input_file%Read_2d(  &
      & variable_name='orostd',    &
      & fill_array = hydro_init_vars%oro_stddev)
    hydro_init_vars%oro_stddev = MERGE(ptr_2D, 0._wp, ptr_2D >= 0._wp)

    ! Read in maximum fraction and depth of local surface depressions
    IF (dsl4jsb_Config(HYDRO_)%l_ponds .AND. .NOT. tile%is_lake .AND. .NOT. tile%is_glacier) THEN
      IF (input_file%Has_var('surf_depr_fract') .AND. input_file%Has_var('surf_depr_depth')) THEN
        ! Use depression area fraction as maximum pond fraction
        ptr_2D => input_file%Read_2d(                                 &
          & variable_name='surf_depr_fract',                           &
          & fill_array=hydro_init_vars%fract_pond_max)
        ptr_2D(:,:) = MERGE(ptr_2D(:,:), 0._wp, ptr_2D(:,:) >= 0._wp)
        ! Use depression depth as maximum pond depth (Attention: Values are negative!)
        ptr_2D => input_file%Read_2d(                                 &
          & variable_name='surf_depr_depth',                           &
          & fill_array=hydro_init_vars%depth_pond_max)
        ptr_2D(:,:) = MERGE(ptr_2D(:,:), 0._wp, ptr_2D(:,:) <= 0._wp)
      ELSE
        CALL finish(TRIM(routine), '*** Error: surf_depr_fract and/or surf_depr_depth not found in bc_sso file. '// &
          & 'Pond scheme cannot run without these fields')
      END IF
    END IF

    ! Close bc_sso input file
    CALL input_file%Close()

    IF (tile%contains_soil) THEN ! for tiles containing soil

      ! Open input file - bc
      input_file = jsb_netcdf_open_input(TRIM(dsl4jsb_Config(HYDRO_)%bc_filename), model%grid_id)

      ! Total soil depth until bedrock
      ptr_2D => input_file%Read_2d(   &
        & variable_name='soil_depth', &
        & fill_array = hydro_init_vars%soil_depth)
      hydro_init_vars%soil_depth = MAX(0.1_wp, MERGE(ptr_2D, 0._wp, ptr_2D >= 0._wp))

      IF (.NOT. dsl4jsb_Config(HYDRO_)%l_soil_texture) THEN
        !>     - Reads maps if not using the soil texture for calculating the soil properties

        ! Volumetric soil porosity
        IF (dsl4jsb_Config(HYDRO_)%l_organic) THEN
          ptr_2D => input_file%Read_2d(                  &
            & fill_array = hydro_init_vars%vol_porosity, &
            & variable_name='soil_porosity_mineral')
        ELSE
          ptr_2D => input_file%Read_2d(                  &
            & fill_array = hydro_init_vars%vol_porosity, &
            & variable_name='soil_porosity')
        END IF
        ptr_2D(:,:) = MERGE(ptr_2D(:,:), 0.2_wp, ptr_2D(:,:) >= 0._wp)

        ! Saturated hydraulic conductivity
        IF (dsl4jsb_Config(HYDRO_)%l_organic) THEN
          ptr_2D => input_file%Read_2d(                  &
            & variable_name='hyd_cond_sat_mineral',      &
            & fill_array = hydro_init_vars%hyd_cond_sat)
        ELSE
          ptr_2D => input_file%Read_2d(                  &
            & variable_name='hyd_cond_sat',              &
            & fill_array = hydro_init_vars%hyd_cond_sat)
        END IF
        ptr_2D(:,:) = MERGE(ptr_2D(:,:), 4.e-6_wp, ptr_2D(:,:) >= 0._wp)

        ! Volumetric field capacity [m/m]
        IF (dsl4jsb_Config(HYDRO_)%l_organic) THEN
          ptr_2D => input_file%Read_2d(                  &
            & variable_name='soil_field_cap_mineral',    &
            & fill_array = hydro_init_vars%vol_field_cap)
        ELSE
          ptr_2D => input_file%Read_2d(                  &
            & variable_name='soil_field_cap',            &
            & fill_array = hydro_init_vars%vol_field_cap)
        END IF
        ptr_2D(:,:) = MERGE(ptr_2D(:,:), 0._wp, ptr_2D(:,:) >= 0._wp)

        ! Volumetric permanent wilting point
        IF (dsl4jsb_Config(HYDRO_)%l_organic) THEN
          ptr_2D => input_file%Read_2d(                  &
            & variable_name='wilting_point_mineral',     &
            & fill_array=hydro_init_vars%vol_p_wilt)
        ELSE
          ptr_2D => input_file%Read_2d(                  &
            & variable_name='wilting_point',             &
            & fill_array=hydro_init_vars%vol_p_wilt)
        END IF
        ptr_2D(:,:) = MERGE(ptr_2D(:,:), 0._wp, ptr_2D(:,:) >= 0._wp)

        ! Volumetric residual water content (if missing deduce from porosity)
        IF (input_file%Has_var('residual_water')) THEN
          IF (dsl4jsb_Config(HYDRO_)%l_organic) THEN
            ptr_2D => input_file%Read_2d(                &
              & variable_name='residual_water_mineral',  &
              & fill_array=hydro_init_vars%vol_wres)
          ELSE
            ptr_2D => input_file%Read_2d(                &
              & variable_name='residual_water',          &
              & fill_array=hydro_init_vars%vol_wres)
          END IF
          ptr_2D(:,:) = MERGE(ptr_2D(:,:), 0._wp, ptr_2D(:,:) >= 0._wp)
        ELSE
          WRITE (message_text,*) 'BC File does not contain data on residual water content. Will use approximation instead!'
          CALL message (TRIM(routine), message_text)
          hydro_init_vars%vol_wres = hydro_init_vars%vol_porosity * 0.2_wp
        END IF

        ! Soil matric potential
        IF (dsl4jsb_Config(HYDRO_)%l_organic) THEN
          ptr_2D => input_file%Read_2d(                  &
            & variable_name='moisture_pot_mineral',      &
            & fill_array = hydro_init_vars%matric_pot)
        ELSE
          ptr_2D => input_file%Read_2d(                  &
            & variable_name='moisture_pot',              &
            & fill_array = hydro_init_vars%matric_pot)
        END IF
        ! @todo There is still old bc data around where the matric potential is given with a positiv sign while
        !       it should be negativ from a soil physics perspective and is expected as such from the soil hydrology
        !       routines. Until this bc data is not used anymore, the following check is required to ensure the
        !       sign of the soil matric potential is (converted) correct.
        ! Exclude values outside of the domain
        WHERE(.NOT. is_in_domain(:,:))
          ptr_2D(:,:) = 0._wp
        ENDWHERE
        ! Exclude (standard) missing values
        WHERE(ABS(ptr_2D(:,:)) > 1.0e+20_wp)
          ptr_2D(:,:) = 0._wp
        ENDWHERE
        ! Check if all values are negative or positive (negative values are needed)
        !   and convert positive -> negative values
        ! TODO please revise if still is correct and necessary with current ic/bc files
        IF (ALL(ptr_2D(:,:) < 1.0e-10_wp)) THEN
          ! Correct: all values are negative
          ptr_2D(:,:) = MERGE(ptr_2D(:,:), -0.2_wp, ptr_2D(:,:) <= 0._wp)
        ELSE IF (ALL(ptr_2D(:,:) > -1.0e-10_wp)) THEN
          ! False: all values are positive --> convert to negative values
          ptr_2D(:,:) = MERGE(ptr_2D(:,:), 0.2_wp, ptr_2D(:,:) >= 0._wp) * (-1._wp)
        ELSE
          ! Totally wrong: dataset contains positive and negativ values
          WRITE (message_text,*) 'Found positive and negative values in soil moisture matric potential: ', &
            & MINVAL(ptr_2D(:,:)),' <-> ',MAXVAL(ptr_2D(:,:))
          CALL finish (TRIM(routine), message_text)
        END IF

        ! Exponent B in Clapp and Hornberger
        IF (dsl4jsb_Config(HYDRO_)%l_organic) THEN
          ptr_2D => input_file%Read_2d(                  &
            & variable_name='bclapp_mineral',            &
            & fill_array = hydro_init_vars%bclapp)
        ELSE
          ptr_2D => input_file%Read_2d(                  &
            & variable_name='bclapp',                    &
            & fill_array = hydro_init_vars%bclapp)
        END IF
        ptr_2D(:,:) = MERGE(ptr_2D(:,:), 6._wp, ptr_2D(:,:) >= 0._wp)

        ! Pore size distribution index
        IF (dsl4jsb_Config(HYDRO_)%l_organic) THEN
          ptr_2D => input_file%Read_2d(                  &
            & variable_name='pore_size_index_mineral',   &
            & fill_array = hydro_init_vars%pore_size_index)
        ELSE
          ptr_2D => input_file%Read_2d(                  &
            & variable_name='pore_size_index',           &
            & fill_array = hydro_init_vars%pore_size_index)
        END IF
        ptr_2D(:,:) = MERGE(ptr_2D(:,:), 0.25_wp, ptr_2D(:,:) >= 0.05_wp)

      ELSE   ! l_soil_texture
        !>     - reads soil textures if using the soil texture for calculating the soil properties

        ! Read fraction of sand for upper soil
        ptr_2D => input_file%Read_2d(                  &
          & variable_name='FR_SAND',                   &
          & fill_array=hydro_init_vars%fr_sand)
        ptr_2D(:,:) = MERGE(ptr_2D(:,:), 0._wp, ptr_2D(:,:) >= 0._wp)

        ! Read fraction of silt for upper soil
        ptr_2D => input_file%Read_2d(                  &
          & variable_name='FR_SILT',                   &
          & fill_array=hydro_init_vars%fr_silt)
        ptr_2D(:,:) = MERGE(ptr_2D(:,:), 0._wp, ptr_2D(:,:) >= 0._wp)

        ! Read fraction of clay for upper soil
        ptr_2D => input_file%Read_2d(                  &
          & variable_name='FR_CLAY',                   &
          & fill_array=hydro_init_vars%fr_clay)
        ptr_2D(:,:) = MERGE(ptr_2D(:,:), 0._wp, ptr_2D(:,:) >= 0._wp)

        ! Read fraction of sand for deep soil
        ptr_2D => input_file%Read_2d(                  &
          & variable_name='SUB_FR_SAND',               &
          & fill_array=hydro_init_vars%fr_sand_deep)
        ptr_2D(:,:) = MERGE(ptr_2D(:,:), 0._wp, ptr_2D(:,:) >= 0._wp)

        ! Read fraction of silt for deep soil
        ptr_2D => input_file%Read_2d(                  &
          & variable_name='SUB_FR_SILT',               &
          & fill_array=hydro_init_vars%fr_silt_deep)
        ptr_2D(:,:) = MERGE(ptr_2D(:,:), 0._wp, ptr_2D(:,:) >= 0._wp)

        ! Read fraction of clay for deep soil
        ptr_2D => input_file%Read_2d(                  &
          & variable_name='SUB_FR_CLAY',               &
          & fill_array=hydro_init_vars%fr_clay_deep)
        ptr_2D(:,:) = MERGE(ptr_2D(:,:), 0._wp, ptr_2D(:,:) >= 0._wp)

        ! Calculate volumetric soil porosity of mineral soil
        CALL soil_init_from_texture( &
          & porosity_sand, porosity_silt, porosity_clay, porosity_oc,                                                &
          & hydro_init_vars%fr_sand(:,:), hydro_init_vars%fr_silt(:,:), hydro_init_vars%fr_clay(:,:),                &
          & hydro_init_vars%fr_sand_deep(:,:), hydro_init_vars%fr_silt_deep(:,:), hydro_init_vars%fr_clay_deep(:,:), &
          & hydro_init_vars%vol_porosity(:,:))

        ! Calculate saturated hydraulic conductivity of mineral soil
        CALL soil_init_from_texture( &
          & hyd_cond_sand, hyd_cond_silt, hyd_cond_clay, hyd_cond_oc,                                                &
          & hydro_init_vars%fr_sand(:,:), hydro_init_vars%fr_silt(:,:), hydro_init_vars%fr_clay(:,:),                &
          & hydro_init_vars%fr_sand_deep(:,:), hydro_init_vars%fr_silt_deep(:,:), hydro_init_vars%fr_clay_deep(:,:), &
          & hydro_init_vars%hyd_cond_sat(:,:))
        hydro_init_vars%hyd_cond_sat(:,:) = hydro_init_vars%hyd_cond_sat(:,:) * hyd_cond_sat_profile

        ! Calculate volumetric field capacity [m/m] of mineral soil
        CALL soil_init_from_texture( &
          & field_cap_sand, field_cap_silt, field_cap_clay, field_cap_oc,                                            &
          & hydro_init_vars%fr_sand(:,:), hydro_init_vars%fr_silt(:,:), hydro_init_vars%fr_clay(:,:),                &
          & hydro_init_vars%fr_sand_deep(:,:), hydro_init_vars%fr_silt_deep(:,:), hydro_init_vars%fr_clay_deep(:,:), &
          & hydro_init_vars%vol_field_cap(:,:))

        ! Calculate plant wilting point of mineral soil
        CALL soil_init_from_texture( &
          & wilt_sand, wilt_silt, wilt_clay, wilt_oc,                                                                &
          & hydro_init_vars%fr_sand(:,:), hydro_init_vars%fr_silt(:,:), hydro_init_vars%fr_clay(:,:),                &
          & hydro_init_vars%fr_sand_deep(:,:), hydro_init_vars%fr_silt_deep(:,:), hydro_init_vars%fr_clay_deep(:,:), &
          & hydro_init_vars%vol_p_wilt(:,:))

        ! Calculate residual soil water content of mineral soil
        CALL soil_init_from_texture( &
           wres_sand, wres_silt, wres_clay, wres_oc,                                                                &
           hydro_init_vars%fr_sand(:,:), hydro_init_vars%fr_silt(:,:), hydro_init_vars%fr_clay(:,:),                &
           hydro_init_vars%fr_sand_deep(:,:), hydro_init_vars%fr_silt_deep(:,:), hydro_init_vars%fr_clay_deep(:,:), &
           hydro_init_vars%vol_wres(:,:))

        ! Calculate matric potential of mineral soil
        CALL soil_init_from_texture( &
           matric_pot_sand, matric_pot_silt, matric_pot_clay, matric_pot_oc,                                        &
           hydro_init_vars%fr_sand(:,:), hydro_init_vars%fr_silt(:,:), hydro_init_vars%fr_clay(:,:),                &
           hydro_init_vars%fr_sand_deep(:,:), hydro_init_vars%fr_silt_deep(:,:), hydro_init_vars%fr_clay_deep(:,:), &
           hydro_init_vars%matric_pot(:,:))

        ! Calculate pore_size_index of mineral soil
        CALL soil_init_from_texture( &
          & pore_size_index_sand, pore_size_index_silt, pore_size_index_clay, pore_size_index_oc,                    &
          & hydro_init_vars%fr_sand(:,:), hydro_init_vars%fr_silt(:,:), hydro_init_vars%fr_clay(:,:),                &
          & hydro_init_vars%fr_sand_deep(:,:), hydro_init_vars%fr_silt_deep(:,:), hydro_init_vars%fr_clay_deep(:,:), &
          & hydro_init_vars%pore_size_index(:,:))

        ! Calculate exponent b in Clapp & Hornberger of mineral soil
        CALL soil_init_from_texture( &
          & bclapp_sand, bclapp_silt, bclapp_clay, bclapp_oc,                                                        &
          & hydro_init_vars%fr_sand(:,:), hydro_init_vars%fr_silt(:,:), hydro_init_vars%fr_clay(:,:),                &
          & hydro_init_vars%fr_sand_deep(:,:), hydro_init_vars%fr_silt_deep(:,:), hydro_init_vars%fr_clay_deep(:,:), &
          & hydro_init_vars%bclapp(:,:))

      END IF ! l_soil_texture

      !>     - Read soil organic carbon fractions for each soil layer
      IF (dsl4jsb_Config(HYDRO_)%l_socmap) THEN
        DO isoil=1,nsoil
          ptr_3D => input_file%Read_2d_extdim( &
            & variable_name='fract_org_sl',    &
            & start_extdim=isoil, end_extdim=isoil, extdim_name='soillev')
          hydro_init_vars%fract_org_sl(:,isoil,:) = MERGE(ptr_3D(:,:,1), 0._wp, ptr_3D(:,:,1) >= 0._wp)
        END DO
      END IF

      ! Root depth [m]
      ! TODO root depth should become a function of PFTs eventually, instead of being a soil property
      ptr_2D => input_file%Read_2d(   &
        & variable_name='root_depth', &
        & fill_array = hydro_init_vars%root_depth)
      hydro_init_vars%root_depth = MERGE(ptr_2D, 0._wp, ptr_2D >= 0._wp)
      ! Make sure there is root depth equivalent to 0.2 m weq in all cells
      !   that contain field capacity data but no root depth
      WHERE (hydro_init_vars%vol_field_cap(:,:) > 1.0e-10_wp)
        hydro_init_vars%root_depth(:,:) = MAX(hydro_init_vars%root_depth(:,:), &
          & 0.2_wp / hydro_init_vars%vol_field_cap(:,:))
      ELSEWHERE
        hydro_init_vars%root_depth(:,:) = 0._wp
      ENDWHERE
      ! Constrain root depth because in the bc file it may be larger then soil depth
      hydro_init_vars%root_depth = MIN(hydro_init_vars%root_depth, hydro_init_vars%soil_depth)

      ! close bc input file
      CALL input_file%Close()

    END IF

    IF (tile%contains_soil) THEN ! for tiles containing soil

      IF (.NOT. model%config%init_from_ifs .OR. dsl4jsb_Config(HYDRO_)%l_read_initial_moist) THEN
        input_file = jsb_netcdf_open_input(dsl4jsb_Config(HYDRO_)%ic_filename, model%grid_id)
      END IF

      ! Initialize hydrology variables from ic file

      IF (model%config%init_from_ifs) THEN
        !>     - If [[t_jsb_model_config:init_from_ifs]] == .TRUE., read soil moisture index and snow water
        !>       equivalent from IFS initial file [[t_jsb_model_config:ifs_filename]]

        ptr_3D => model%config%ifs_input_file%Read_2d_time( &
          & variable_name='W_SNOW', start_time_step=1, end_time_step=1)
        hydro_init_vars%ifs_weq_snow(:,:) = ptr_3D(:,:,1)

      ELSE

        ! Initial snow depth
        IF (input_file%Has_var('snow')) THEN
          ptr_2D => input_file%Read_2d( &
            & variable_name='snow',     &
            & fill_array = hydro_init_vars%weq_snow_soil)
          hydro_init_vars%weq_snow_soil = MERGE(ptr_2D, 0._wp, ptr_2D >= 0._wp)
        ELSE
          hydro_init_vars%weq_snow_soil = 0._wp
          CALL message(TRIM(routine), 'Initializing snow to zero')
        END IF

      END IF

      IF (model%config%init_from_ifs .AND. .NOT. dsl4jsb_Config(HYDRO_)%l_read_initial_moist) THEN

        DO isoil=1,ifs_nsoil
          ptr_2D => model%config%ifs_input_file%Read_2d_1lev_1time( &
            & variable_name='SMIL'//int2string(isoil))
          hydro_init_vars%ifs_smi_sl(:,isoil,:) = ptr_2D(:,:)
        END DO

      ELSE

        DO isoil=1,nsoil
          ptr_3D => input_file%Read_2d_extdim( &
            & variable_name='layer_moist',     &
            & start_extdim=isoil, end_extdim=isoil, extdim_name='soillev')
          hydro_init_vars%wtr_soil_sl(:,isoil,:) = MERGE(ptr_3D(:,:,1), 0._wp, ptr_3D(:,:,1) >= 0._wp)
        END DO

      END IF ! init_from_ifs

      IF (.NOT. model%config%init_from_ifs .OR. dsl4jsb_Config(HYDRO_)%l_read_initial_moist) THEN
        CALL input_file%Close()
      END IF

    END IF ! contains_soil

  END SUBROUTINE hydro_read_init_vars

  ! ====================================================================================================== !
  !
  !> Set soil properties for hydrology process on tiles
  !>
  !> This subroutine uses the input variables already read from [[mo_hydro_init:hydro_init_vars]]
  !> in [[hydro_read_init_vars]].
  !
  SUBROUTINE hydro_init_bc(tile)

    USE mo_util,          ONLY: soil_depth_to_layers_2d
    USE mo_hydro_process, ONLY: calc_orographic_features

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile !< current tile

    dsl4jsb_Def_config(HYDRO_) !< configuration of `HYDRO_` process
    dsl4jsb_Def_memory(HYDRO_) !< memory of `HYDRO_` process

    TYPE(t_jsb_model), POINTER :: model  !< pointer to model
    TYPE(t_jsb_grid),  POINTER :: grid   !< pointer to horizontal grid
    TYPE(t_jsb_vgrid), POINTER :: soil_w !< pointer to vertical hydrology grid
    INTEGER                    :: nsoil  !< number of soil layers

    INTEGER :: i

    CHARACTER(len=*), PARAMETER :: routine = modname//':hydro_init_bc'

    ! Do nothing if hydrology process is not active on this tile
    IF (.NOT. tile%Is_process_active(HYDRO_)) RETURN

    IF (debug_on()) CALL message(routine, 'Setting  hydrology boundary conditions for tile '//TRIM(tile%name))

    model => Get_model(tile%owner_model_id)

    grid   => Get_grid(model%grid_id)
    soil_w => Get_vgrid('soil_depth_water')
    nsoil = soil_w%n_levels

    dsl4jsb_Get_config(HYDRO_)
    dsl4jsb_Get_memory(HYDRO_)

    ! -------------------------------------------------------------------------------------------------- !
    !>
    !> - Set orographic parameters
    !
    !>     - Elevation (for standalone model)
    IF (jsbach_runs_standalone()) dsl4jsb_var2D_onDomain(HYDRO_,elevation) = hydro_init_vars%elevation
    !>     - Standard deviation of orography
    dsl4jsb_var2D_onDomain(HYDRO_,oro_stddev) = hydro_init_vars%oro_stddev

    !$ACC UPDATE ASYNC(1) &
    !$ACC   DEVICE(dsl4jsb_var2D_onDomain(HYDRO_,elevation)) &
    !$ACC   DEVICE(dsl4jsb_var2D_onDomain(HYDRO_,oro_stddev))

    !>     - Steepness parameter for surface runoff (see [[arno_scheme(subroutine)]])
    DO i = 1, grid%Get_nblks()
      CALL calc_orographic_features(                    &
        ! in
        & grid%Get_nproma(),                            &
        & grid%nlat_g,                                  & ! (Effective) number of latitudes for steepness parameter
        & dsl4jsb_var_ptr(HYDRO_, oro_stddev ) (:,i),   &
        ! out
        & dsl4jsb_var_ptr(HYDRO_, steepness  ) (:,i)    &
        & )
    END DO

    IF (.NOT. tile%is_lake .AND. .NOT. tile%is_glacier) THEN
      ! Maximum depth and extent of surface water ponds
      IF (dsl4jsb_Config(HYDRO_)%l_ponds) THEN
        ! Note, surface depression depth is given as negative value in the data, but represent
        !   a storage on top of the soil --> needs to be converted to positive values!
        dsl4jsb_var_ptr(HYDRO_, fract_pond_max) = hydro_init_vars%fract_pond_max(:,:)
        dsl4jsb_var_ptr(HYDRO_, weq_pond_max)   = hydro_init_vars%depth_pond_max(:,:) * (-1._wp)
      ELSE
        ! Note, that ponds are implicitely disabled by setting the maximum allowed pond fraction
        !   and depth to zero. Thus, we can avoid a large number of l_ponds conditions for surface
        !   and soil hydrology processes.
        dsl4jsb_var_ptr(HYDRO_, fract_pond_max) = 0._wp
        dsl4jsb_var_ptr(HYDRO_, weq_pond_max)   = 0._wp
      END IF
    END IF

    IF (tile%contains_soil) THEN
      !>
      !> If `tile` contains soil:
      !>
      ! -------------------------------------------------------------------------------------------------- !
      !>
      !> - Set total and layer soil depths
      !
      !>     - Total soil depth until bedrock
      dsl4jsb_var2D_onDomain(HYDRO_, soil_depth) = hydro_init_vars%soil_depth
      !>     - Actual soil layer depths based on the fixed layer thicknesses from the vertical grid and
      !>       the total bedrock depth
      dsl4jsb_var3D_onDomain(HYDRO_, soil_depth_sl) = soil_depth_to_layers_2d( &
        & dsl4jsb_var2D_onDomain(HYDRO_, soil_depth), & ! Total soil depth until bedrock
        & soil_w%dz(:))                                 ! Fixed soil layer thicknesses from vertical grid
      !>
      !
      ! -------------------------------------------------------------------------------------------------- !
      !>
      !> - Set hydrologic properties
      !
      !>     - Volumetric soil porosity
      dsl4jsb_var_ptr(HYDRO_, vol_porosity) = hydro_init_vars%vol_porosity
      !>     - Saturated hydraulic conductivity
      dsl4jsb_var_ptr(HYDRO_, hyd_cond_sat) = hydro_init_vars%hyd_cond_sat

      !>     - Matric potential
      dsl4jsb_var_ptr(HYDRO_, matric_pot) = hydro_init_vars%matric_pot

      !>     - Exponent B in Clapp and Hornberger
      dsl4jsb_var_ptr(HYDRO_, bclapp) = hydro_init_vars%bclapp
      !>     - Pore size distribution index
      dsl4jsb_var_ptr(HYDRO_, pore_size_index) = hydro_init_vars%pore_size_index
      !>     - Volumetric field capacity
      dsl4jsb_var_ptr(HYDRO_, vol_field_cap) = hydro_init_vars%vol_field_cap
      !>     - Get water content at field capacity by multiplying volumetric field capacity with
      !>       soil layer thicknesses
      !>       _Note_: At this point a potential ice fraction is ignored.
      DO i=1,nsoil
        dsl4jsb_var_ptr(HYDRO_,wtr_soil_fc_sl) (:,i,:) = hydro_init_vars%vol_field_cap(:,:) * &
          & dsl4jsb_var_ptr(HYDRO_,soil_depth_sl) (:,i,:)
      END DO
      !>     - Volumetric permanent wilting point
      ! TODO: Do we need this for bare soil tiles?
      dsl4jsb_var_ptr(HYDRO_, vol_p_wilt) = hydro_init_vars%vol_p_wilt
      !>     - Get permanent wilting point by multiplying volumetric permanent wilting point
      !>       with soil layer thicknesses
      !>       _Note_: At this point a potential ice fraction is ignored.
      DO i=1,nsoil
        dsl4jsb_var_ptr(HYDRO_,wtr_soil_pwp_sl) (:,i,:) = hydro_init_vars%vol_p_wilt(:,:) * &
          & dsl4jsb_var_ptr(HYDRO_,soil_depth_sl) (:,i,:)
      END DO

      !>     - Volumetric residual water contant
      dsl4jsb_var_ptr(HYDRO_, vol_wres) = hydro_init_vars%vol_wres

      ! Get absolute residual water content by multiplying with soil layer thicknesses
      ! Note: At this point a potential ice fraction is ignored.
      DO i=1,nsoil
        dsl4jsb_var_ptr(HYDRO_,wtr_soil_res_sl) (:,i,:) = hydro_init_vars%vol_wres(:,:) * &
          & dsl4jsb_var_ptr(HYDRO_,soil_depth_sl) (:,i,:)
      END DO

      ! Volumetric residual water contant
      dsl4jsb_var_ptr(HYDRO_, vol_wres) = hydro_init_vars%vol_wres

      !>     - Get absolute residual water content by multiplying with soil layer thicknesses
      !>       _Note: At this point a potential ice fraction is ignored.
      DO i=1,nsoil
        dsl4jsb_var_ptr(HYDRO_,wtr_soil_res_sl) (:,i,:) = hydro_init_vars%vol_wres(:,:) * &
          & dsl4jsb_var_ptr(HYDRO_,soil_depth_sl) (:,i,:)
      END DO

      !>     - Volumetric residual water contant
      dsl4jsb_var_ptr(HYDRO_, vol_wres) = hydro_init_vars%vol_wres

      ! Get absolute residual water content by multiplying with soil layer thicknesses
      ! Note: At this point a potential ice fraction is ignored.
      DO i=1,nsoil
        dsl4jsb_var_ptr(HYDRO_,wtr_soil_res_sl) (:,i,:) = hydro_init_vars%vol_wres(:,:) * &
          & dsl4jsb_var_ptr(HYDRO_,soil_depth_sl) (:,i,:)
      END DO

      !>     - Get water holding capacity (water content at saturation) by multiplying volumetric porosity
      !>       with soil layer thicknesses
      !>       _Note_: At this point a potential ice fraction is ignored.
      DO i=1,nsoil
        dsl4jsb_var_ptr(HYDRO_,wtr_soil_sat_sl) (:,i,:) = dsl4jsb_var2D_onDomain(HYDRO_,vol_porosity) * &
          & dsl4jsb_var_ptr(HYDRO_,soil_depth_sl) (:,i,:)
      END DO

      !>     - Soil organic carbon fraction from map if [[t_hydro_config:l_socmap]] == .TRUE.
      IF (dsl4jsb_Config(HYDRO_)%l_socmap) THEN
        dsl4jsb_var_ptr(HYDRO_, fract_org_sl) = hydro_init_vars%fract_org_sl
      END IF

      !$ACC UPDATE ASYNC(1) &
      !$ACC   DEVICE(dsl4jsb_var_ptr       (HYDRO_, wtr_soil_sat_sl)) &
      !$ACC   DEVICE(dsl4jsb_var_ptr       (HYDRO_, wtr_soil_pwp_sl)) &
      !$ACC   DEVICE(dsl4jsb_var_ptr       (HYDRO_, wtr_soil_fc_sl)) &
      !$ACC   DEVICE(dsl4jsb_var_ptr       (HYDRO_, wtr_soil_res_sl)) &
      !$ACC   DEVICE(dsl4jsb_var_ptr       (HYDRO_, vol_field_cap)) &
      !$ACC   DEVICE(dsl4jsb_var_ptr       (HYDRO_, vol_p_wilt)) &
      !$ACC   DEVICE(dsl4jsb_var_ptr       (HYDRO_, vol_wres)) &
      !$ACC   DEVICE(dsl4jsb_var_ptr       (HYDRO_, pore_size_index)) &
      !$ACC   DEVICE(dsl4jsb_var_ptr       (HYDRO_, bclapp)) &
      !$ACC   DEVICE(dsl4jsb_var_ptr       (HYDRO_, matric_pot)) &
      !$ACC   DEVICE(dsl4jsb_var_ptr       (HYDRO_, hyd_cond_sat)) &
      !$ACC   DEVICE(dsl4jsb_var_ptr       (HYDRO_, vol_porosity)) &
      !$ACC   DEVICE(dsl4jsb_var_ptr       (HYDRO_, fract_org_sl)) &
      !$ACC   DEVICE(dsl4jsb_var_ptr       (HYDRO_, fract_pond_max)), &
      !$ACC   DEVICE(dsl4jsb_var_ptr       (HYDRO_, weq_pond_max)), &
      !$ACC   DEVICE(dsl4jsb_var3D_onDomain(HYDRO_, soil_depth_sl)) &
      !$ACC   DEVICE(dsl4jsb_var2D_onDomain(HYDRO_, soil_depth))

    END IF ! contains_soil

    IF (tile%contains_vegetation) THEN ! tile contains vegetation
      !>
      !> If `tile` contains vegetation:
      !>
      ! -------------------------------------------------------------------------------------------------- !
      !>
      !> - Set total and layer root depths
      !
      !>     - Total root depth
      dsl4jsb_var_ptr(HYDRO_, root_depth) = hydro_init_vars%root_depth
      !>     - Actual root depths based on the fixed layer thicknesses from the vertical grid and
      !>       the total root depth
      dsl4jsb_var3D_onDomain(HYDRO_, root_depth_sl) = soil_depth_to_layers_2d( &
        & dsl4jsb_var2D_onDomain(HYDRO_, root_depth), & ! Total rooting depth
        & soil_w%dz(:))                                 ! Soil layer thicknesses from vertical grid

      !> - Set maximum root zone soil moisture
      !>   _Note1: This corresponds to water content at field capacity.
      dsl4jsb_var2D_onDomain(HYDRO_, wpi_rootzone_max) = &
        & hydro_init_vars%vol_field_cap * hydro_init_vars%root_depth
      !>   _Note2_: Uses [[t_hydro_config:w_soil_limit]] > 0 as upper limit but this is disabled per default.
      IF (dsl4jsb_Config(HYDRO_)%w_soil_limit > 0._wp) THEN
        dsl4jsb_var2D_onDomain(HYDRO_, wpi_rootzone_max) = MIN( &
          & dsl4jsb_var2D_onDomain(HYDRO_, wpi_rootzone_max), dsl4jsb_Config(HYDRO_)%w_soil_limit)
      END IF

      !$ACC UPDATE ASYNC(1) &
      !$ACC   DEVICE(dsl4jsb_var2D_onDomain(HYDRO_, root_depth)) &
      !$ACC   DEVICE(dsl4jsb_var3D_onDomain(HYDRO_, root_depth_sl)) &
      !$ACC   DEVICE(dsl4jsb_var2D_onDomain(HYDRO_, wpi_rootzone_max))
    END IF

  END SUBROUTINE hydro_init_bc

  ! ====================================================================================================== !
  !
  !> Set initial conditions for hydrology process on tiles
  !>
  !> This subroutine uses the input variables stored in [[mo_hydro_init:hydro_init_vars]] that had been
  !> read in [[hydro_read_init_vars]].
  !
  SUBROUTINE hydro_init_ic(tile)

    USE mo_hydro_util,          ONLY: get_amount_in_rootzone
    USE mo_util,                ONLY: ifs2soil

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile !< current tile

    dsl4jsb_Def_config(HYDRO_) !< configuration of `HYDRO_` process
    dsl4jsb_Def_memory(HYDRO_) !< memory of `HYDRO_` process
    dsl4jsb_Def_memory(PHENO_) !< memory of `PHENO_` process

    TYPE(t_jsb_model), POINTER :: model   !< pointer to model
    TYPE(t_jsb_vgrid), POINTER :: soil_w  !< pointer to vertical hydrology grid

    INTEGER :: i

    CHARACTER(len=*), PARAMETER :: routine = modname//':hydro_init_ic'

    ! Do nothing if hydrology process is not active on this tile
    IF (.NOT. tile%Is_process_active(HYDRO_)) RETURN

    IF (debug_on()) CALL message(TRIM(routine), 'Initializing hydrology memory for tile '//TRIM(tile%name))

    model  => Get_model(tile%owner_model_id)
    soil_w => Get_vgrid('soil_depth_water')

    dsl4jsb_Get_config(HYDRO_)
    dsl4jsb_Get_memory(HYDRO_)
    IF (tile%contains_vegetation) THEN
      dsl4jsb_Get_memory(PHENO_)
    END IF

    !----------

    IF (tile%contains_soil) THEN ! for tiles containing soil
      !>
      !> If `tile` contains soil:
      !>
      ! -------------------------------------------------------------------------------------------------- !
      !>
      !> - Set initial soil water and snow water equivalent If [[t_jsb_model_config:init_from_ifs]] == .TRUE.
      !

      IF (model%config%init_from_ifs) THEN
      ! Initial snow depth
        dsl4jsb_var_ptr(HYDRO_, weq_snow_soil)(:,:) = hydro_init_vars%ifs_weq_snow(:,:)
        dsl4jsb_var2D_onDomain(HYDRO_, weq_snow) = dsl4jsb_var2D_onDomain(HYDRO_, weq_snow_soil)
      ELSE
        dsl4jsb_var_ptr(HYDRO_, weq_snow_soil) = hydro_init_vars%weq_snow_soil
        dsl4jsb_var2D_onDomain(HYDRO_, weq_snow) = dsl4jsb_var2D_onDomain(HYDRO_, weq_snow_soil)
      END IF

      IF (model%config%init_from_ifs .AND. .NOT. dsl4jsb_Config(HYDRO_)%l_read_initial_moist) THEN
        !>     - Interpolate soil moisture index from IFS soil layers to our soil layers
        !        Note: wtr_soil_sl is used here to temporarilly store the soil moisture index interpolated to
        !              our soil layers
        CALL ifs2soil(hydro_init_vars%ifs_smi_sl(:,:,:), ifs_soil_depth(:),        &
          &           dsl4jsb_var_ptr(HYDRO_,wtr_soil_sl)(:,:,:), soil_w%ubounds(:))
        !>     - Convert soil moisture index to water content using field capacity and wilting point
        WHERE (dsl4jsb_var_ptr(HYDRO_,wtr_soil_sl)(:,:,:) < -999._wp) ! Fill missing values
          dsl4jsb_var_ptr(HYDRO_,wtr_soil_sl)(:,:,:) = &
            &  0.5_wp * (dsl4jsb_var_ptr(HYDRO_,wtr_soil_pwp_sl)(:,:,:) + dsl4jsb_var_ptr(HYDRO_,wtr_soil_fc_sl) (:,:,:))
        ELSE WHERE
          dsl4jsb_var_ptr(HYDRO_,wtr_soil_sl)(:,:,:) = &
            & MIN(dsl4jsb_var_ptr(HYDRO_,wtr_soil_sat_sl)(:,:,:), &
            &     MAX(0._wp, &
            &         dsl4jsb_var_ptr(HYDRO_,wtr_soil_sl)(:,:,:) &
            &         * (dsl4jsb_var_ptr(HYDRO_,wtr_soil_fc_sl) (:,:,:) - dsl4jsb_var_ptr(HYDRO_,wtr_soil_pwp_sl)(:,:,:)) &
            &         + dsl4jsb_var_ptr(HYDRO_,wtr_soil_pwp_sl)(:,:,:) &
            &        ) &
            &    )
        END WHERE

      ELSE ! init_from_ifs
        !>
        !> - Set initial soil water and snow water equivalent If [[t_jsb_model_config:init_from_ifs]] == .FALSE.
        !>
        !>     - Initialize soil water and snow water equivalent with values from [[t_hydro_config:ic_filename]]
        dsl4jsb_var_ptr(HYDRO_,wtr_soil_sl) (:,:,:) = hydro_init_vars%wtr_soil_sl(:,:,:)

      END IF ! init_from_ifs

      ! Initial soil water overflow (used in terra planet experiments)
      dsl4jsb_var2D_onDomain(HYDRO_, tpe_overflow) = 0._wp

      !$ACC UPDATE ASYNC(1) &
      !$ACC   DEVICE(dsl4jsb_var2D_onDomain(HYDRO_, weq_snow))      &
      !$ACC   DEVICE(dsl4jsb_var2D_onDomain(HYDRO_, weq_snow_soil)) &
      !$ACC   DEVICE(dsl4jsb_var3D_onDomain(HYDRO_, wtr_soil_sl))   &
      !$ACC   DEVICE(dsl4jsb_var2D_onDomain(HYDRO_, tpe_overflow))

      ! -------------------------------------------------------------------------------------------------- !
      !>
      !> - If [[t_hydro_config:l_organic]] == .TRUE. but [[t_hydro_config:l_socmap]] == .FALSE., initialize
      !>   organic soil layer fractions `fract_org_sl` on vegetated tiles in upper layer from forest fraction
      !>   and set them to zero elsewhere.
      !>   _Note_: At the moment, `fract_org_sl` doesn't change over time, but it could.
      !>
      IF (dsl4jsb_Config(HYDRO_)%l_organic) THEN
        IF (.NOT. dsl4jsb_Config(HYDRO_)%l_socmap) THEN ! Not using map for soil organic carbon
          ! Set to zero
          dsl4jsb_var3D_onDomain(HYDRO_, fract_org_sl) = 0._wp
          !
          IF (tile%is_vegetation) THEN ! tile has vegetation
            IF (tile%lcts(1)%lib_id /= 0)  THEN ! lctlib is only defined for PFTs
              ! For PFT tiles
              ! Note: in current setup hydrology does not run on PFTs
              IF (dsl4jsb_Lctlib_param(ForestFlag)) THEN
                dsl4jsb_var_ptr(HYDRO_, fract_org_sl)(:,1,:) = dsl4jsb_var2D_onDomain(PHENO_, fract_fpc_max)
              END IF
            ELSE
              ! For general vegetated tiles that are not PFT
              dsl4jsb_var_ptr(HYDRO_, fract_org_sl)(:,1,:) = dsl4jsb_var2D_onDomain(PHENO_, fract_forest)
            END IF
          END IF ! is_vegetation
        END IF ! l_socmap

        ! Set to zero if we don't have soil data
        WHERE (dsl4jsb_var3D_onDomain(HYDRO_, wtr_soil_fc_sl) <= 1.0E-10_wp)
          dsl4jsb_var3D_onDomain(HYDRO_, fract_org_sl) = 0._wp
        ENDWHERE

        !$ACC UPDATE DEVICE(dsl4jsb_var3D_onDomain(HYDRO_, fract_org_sl)) ASYNC(1)

      END IF ! l_organic

      ! update soil properties
      ! modify saturated, field capacity and permanent
      !   wilting point capacities accordingly.
      !> - Call [[mo_hydro_init:init_soil_properties]] to update soil properties
      !>     - distribute 2d initial fields to soil layers
      !>     - aggregate contributions from mineral and organic content
      CALL init_soil_properties(tile)

      ! Make sure that soil water content is not larger than the maximum soil water capacity.
      ! Note: the volumetric porosity represents the maximum water capacity of the specific soil layer. Reduce
      ! soil water and ice if this limit is exceeded.
      IF (ANY(dsl4jsb_var3D_onDomain(HYDRO_, wtr_soil_sl)     + dsl4jsb_var3D_onDomain(HYDRO_, ice_soil_sl) > &
        &     dsl4jsb_var3D_onDomain(HYDRO_, vol_porosity_sl) * dsl4jsb_var3D_onDomain(HYDRO_, soil_depth_sl))) THEN
        ! First limit soil water to maximum soil water capacity
        dsl4jsb_var3D_onDomain(HYDRO_, wtr_soil_sl) = MAX(0._wp, &
          & MIN(dsl4jsb_var3D_onDomain(HYDRO_, wtr_soil_sl),     &
          &     dsl4jsb_var3D_onDomain(HYDRO_, vol_porosity_sl) * dsl4jsb_var3D_onDomain(HYDRO_, soil_depth_sl)))
        ! Then limit soil ice to the remaining capacity (without soil water)
        dsl4jsb_var3D_onDomain(HYDRO_, ice_soil_sl) = MAX(0._wp, &
          & MIN(dsl4jsb_var3D_onDomain(HYDRO_, ice_soil_sl),     &
          &     dsl4jsb_var3D_onDomain(HYDRO_, vol_porosity_sl) * dsl4jsb_var3D_onDomain(HYDRO_, soil_depth_sl) &
          &     - dsl4jsb_var3D_onDomain(HYDRO_, wtr_soil_sl)))
      END IF

      ! The following can change later depending on the organic layers and soil ice content.
      ! TODO At the time hydro_init is called ice_soil_sl hasn't been initialized, yet, i.e. it should be 0 ?
      dsl4jsb_var3D_onDomain(HYDRO_, wtr_soil_sat_sl) = MAX(0._wp,                                          &
        & dsl4jsb_var3D_onDomain(HYDRO_, vol_porosity_sl)  * dsl4jsb_var3D_onDomain(HYDRO_, soil_depth_sl)  &
        & - dsl4jsb_var3D_onDomain(HYDRO_, ice_soil_sl))
      dsl4jsb_var3D_onDomain(HYDRO_, wtr_soil_fc_sl) =  MAX(0._wp,                                          &
        & dsl4jsb_var3D_onDomain(HYDRO_, vol_field_cap_sl) * dsl4jsb_var3D_onDomain(HYDRO_, soil_depth_sl)  &
        & - dsl4jsb_var3D_onDomain(HYDRO_, ice_soil_sl))
      dsl4jsb_var3D_onDomain(HYDRO_, wtr_soil_pwp_sl) = MAX(0._wp,                                          &
        & dsl4jsb_var3D_onDomain(HYDRO_, vol_p_wilt_sl)    * dsl4jsb_var3D_onDomain(HYDRO_, soil_depth_sl)  &
        & - dsl4jsb_var3D_onDomain(HYDRO_, ice_soil_sl))
      dsl4jsb_var3D_onDomain(HYDRO_, wtr_soil_res_sl) = MAX(0._wp,                                          &
        & dsl4jsb_var3D_onDomain(HYDRO_, vol_wres_sl)      * dsl4jsb_var3D_onDomain(HYDRO_, soil_depth_sl)  &
        & - dsl4jsb_var3D_onDomain(HYDRO_, ice_soil_sl))

      !$ACC UPDATE ASYNC(1) &
      !$ACC   DEVICE(dsl4jsb_var3D_onDomain(HYDRO_, wtr_soil_sl))      &
      !$ACC   DEVICE(dsl4jsb_var3D_onDomain(HYDRO_, ice_soil_sl))      &
      !$ACC   DEVICE(dsl4jsb_var3D_onDomain(HYDRO_, wtr_soil_sat_sl))  &
      !$ACC   DEVICE(dsl4jsb_var3D_onDomain(HYDRO_, wtr_soil_fc_sl))   &
      !$ACC   DEVICE(dsl4jsb_var3D_onDomain(HYDRO_, wtr_soil_pwp_sl))  &
      !$ACC   DEVICE(dsl4jsb_var3D_onDomain(HYDRO_, wtr_soil_res_sl))

    END IF ! contains_soil

    IF (tile%contains_vegetation) THEN   ! tile contains vegetation
      !>
      !> If `tile` contains vegetation (i.e. soil and roots present):
      !>
      ! -------------------------------------------------------------------------------------------------- !
      !>
      !> - Re-compute maximum root zone moisture based on soil parameters modified by organic fractions
      !>
      !>   _Note: All plant related computations use wpi_rootzone_max reduced by ice and supercooled water content,
      !>          which represents the unfrozen part of the soil.
      !>          This assumes that plants retract/extents their roots immediately if the soil freezes/thaws.
      !
      DO i=1,SIZE(dsl4jsb_var_ptr(HYDRO_, wpi_rootzone_max), 2)
          CALL get_amount_in_rootzone( &
          & dsl4jsb_var_ptr(HYDRO_, wtr_soil_fc_sl)   (:,:,i), &
          & dsl4jsb_var_ptr(HYDRO_, soil_depth_sl)    (:,:,i), &
          & dsl4jsb_var_ptr(HYDRO_, root_depth_sl)    (:,:,i), &
          & dsl4jsb_var_ptr(HYDRO_, wpi_rootzone_max) (:,  i)  )
      END DO
      !$ACC UPDATE HOST(dsl4jsb_var_ptr(HYDRO_, wpi_rootzone_max)) ASYNC(1)
      !$ACC WAIT(1)
      ! -------------------------------------------------------------------------------------------------- !
      !>
      !> - Set maximum root zone moisture to 0.2 for glacier and ocean areas
      !>
      ! Set maximum root zone moisture to 0.2 for glacier and ocean areas (where wpi_rootzone_max == 0)
      dsl4jsb_var2D_onDomain(HYDRO_, wpi_rootzone_max) = MERGE(       &
          & dsl4jsb_var2D_onDomain(HYDRO_, wpi_rootzone_max), 0.2_wp, &
          & dsl4jsb_var2D_onDomain(HYDRO_, wpi_rootzone_max) > 0._wp )
      ! DOCU not used by default
      IF (dsl4jsb_Config(HYDRO_)%w_soil_limit > 0._wp) THEN
        dsl4jsb_var2D_onDomain(HYDRO_, wpi_rootzone_max) = MIN( &
        & dsl4jsb_var2D_onDomain(HYDRO_, wpi_rootzone_max), dsl4jsb_Config(HYDRO_)%w_soil_limit)
      END IF
      !$ACC UPDATE DEVICE(dsl4jsb_var_ptr(HYDRO_, wpi_rootzone_max)) ASYNC(1)

      ! Initialization of relative root zone soil moisture (only water, not ice).
      !   This is needed for the LoGro-P phenology, as phenology is updated before hydrology.

      ! wtr_rootzone calculated here is only used to calculate relative soil moisture for phenology.
      ! It is re-calculated in [[mo_hydro_interface:update_water_stress]].
      DO i=1,SIZE(dsl4jsb_var_ptr(HYDRO_, wtr_rootzone), 2)
          CALL get_amount_in_rootzone( &
          & dsl4jsb_var_ptr(HYDRO_, wtr_soil_sl)    (:,:,i), &
          & dsl4jsb_var_ptr(HYDRO_, soil_depth_sl)  (:,:,i), &
          & dsl4jsb_var_ptr(HYDRO_, root_depth_sl)  (:,:,i), &
          & dsl4jsb_var_ptr(HYDRO_, wtr_rootzone)   (:,  i)  )
          CALL get_amount_in_rootzone( &
          & dsl4jsb_var_ptr(HYDRO_, ice_soil_sl)    (:,:,i), &
          & dsl4jsb_var_ptr(HYDRO_, soil_depth_sl)  (:,:,i), &
          & dsl4jsb_var_ptr(HYDRO_, root_depth_sl)  (:,:,i), &
          & dsl4jsb_var_ptr(HYDRO_, ice_rootzone)   (:,  i)  )
      END DO
      !$ACC UPDATE ASYNC(1) &
      !$ACC HOST(dsl4jsb_var_ptr(HYDRO_, wtr_rootzone)) &
      !$ACC HOST(dsl4jsb_var_ptr(HYDRO_, ice_rootzone))
      !$ACC WAIT(1)

      ! Make sure wtr_rootzone is smaller than wpi_rootzone_max - ice_rootzone.
      ! Note: this corresponds to the maximum possible amount of liquid water in the root zone.
      WHERE (dsl4jsb_var2D_onDomain(HYDRO_, wtr_rootzone) > &
        & dsl4jsb_var2D_onDomain(HYDRO_, wpi_rootzone_max) - dsl4jsb_var2D_onDomain(HYDRO_, ice_rootzone))
        dsl4jsb_var2D_onDomain(HYDRO_, wtr_rootzone) = MAX(0._wp, &
        & dsl4jsb_var2D_onDomain(HYDRO_, wpi_rootzone_max) - dsl4jsb_var2D_onDomain(HYDRO_, ice_rootzone))
      END WHERE
      !$ACC UPDATE DEVICE(dsl4jsb_var_ptr(HYDRO_, wtr_rootzone)) ASYNC(1)

      ! Initialization of wtr_rootzone_rel
      !   This is re-calculated every time step in [[mo_hydro_interface:update_water_stress]].
      WHERE (dsl4jsb_var2D_onDomain(HYDRO_, wpi_rootzone_max) - dsl4jsb_var2D_onDomain(HYDRO_, ice_rootzone) > 0._wp)
        dsl4jsb_var2D_onDomain(HYDRO_, wtr_rootzone_rel) = &
          &  dsl4jsb_var2D_onDomain(HYDRO_, wtr_rootzone) / &
          & (dsl4jsb_var2D_onDomain(HYDRO_, wpi_rootzone_max) - dsl4jsb_var2D_onDomain(HYDRO_, ice_rootzone))
      ELSE WHERE
        dsl4jsb_var2D_onDomain(HYDRO_, wtr_rootzone_rel) = 0._wp
      END WHERE

      !$ACC UPDATE DEVICE(dsl4jsb_var_ptr(HYDRO_, wtr_rootzone_rel)) ASYNC(1)

    END IF ! Tile with vegetation

  END SUBROUTINE hydro_init_ic

  ! ====================================================================================================== !
  !
  SUBROUTINE init_soil_properties(tile)
    ! This routine is needed because soil snow and energy (SSE) processes are calculated before
    ! the soil properties are calculated the first time within the hydrology in update_soil_properties.

    USE mo_hydro_constants, ONLY: &
      ! Parameters for organic soil component in top layer
      & vol_porosity_org_top, hyd_cond_sat_org_top, bclapp_org_top, matric_pot_org_top, &
      & pore_size_index_org_top, vol_field_cap_org_top, vol_p_wilt_org_top, vol_wres_org_top, &
      ! Parameters for organic soil component in layers below top layer
      & vol_porosity_org_below, hyd_cond_sat_org_below, bclapp_org_below, matric_pot_org_below, &
      & pore_size_index_org_below, vol_field_cap_org_below, vol_p_wilt_org_below, vol_wres_org_below

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile

    dsl4jsb_Def_config(HYDRO_) !< configuration of `HYDRO_` process
    dsl4jsb_Def_memory(HYDRO_) !< memory of `HYDRO_` process

    ! Pointers to variables in memory
    dsl4jsb_Real2D_onDomain :: &
      & hyd_cond_sat,          &
      & vol_porosity,          &
      & bclapp,                &
      & matric_pot,            &
      & pore_size_index,       &
      & vol_field_cap,         &
      & vol_p_wilt,            &
      & vol_wres
    dsl4jsb_Real3D_onDomain :: &
      & fract_org_sl,          &
      & hyd_cond_sat_sl,       &
      & vol_porosity_sl,       &
      & bclapp_sl,             &
      & matric_pot_sl,         &
      & pore_size_index_sl,    &
      & vol_field_cap_sl,      &
      & vol_p_wilt_sl,         &
      & vol_wres_sl

    ! Locally allocated vectors
    !
    TYPE(t_jsb_model), POINTER :: model
    TYPE(t_jsb_vgrid), POINTER :: soil_w

    REAL(wp) ::  hyd_cond_sat_org, vol_porosity_org, bclapp_org, matric_pot_org, pore_size_index_org, &
      & vol_field_cap_org, vol_p_wilt_org, vol_wres_org

    INTEGER  :: nsoil, i

    CHARACTER(len=*), PARAMETER :: routine = modname//':init_soil_properties'

    ! Do nothing if hydrology process is not active on this tile
    IF (.NOT. tile%Is_process_active(HYDRO_)) RETURN

    ! Do nothing if this is a lake tile
    IF (tile%is_lake) RETURN

    model => Get_model(tile%owner_model_id)

    soil_w => Get_vgrid('soil_depth_water')
    nsoil = soil_w%n_levels

    dsl4jsb_Get_config(HYDRO_)
    dsl4jsb_Get_memory(HYDRO_)

    dsl4jsb_Get_var2D_onDomain(HYDRO_, hyd_cond_sat)       ! in
    dsl4jsb_Get_var2D_onDomain(HYDRO_, vol_porosity)       ! in
    dsl4jsb_Get_var2D_onDomain(HYDRO_, bclapp)             ! in
    dsl4jsb_Get_var2D_onDomain(HYDRO_, matric_pot)         ! in
    dsl4jsb_Get_var2D_onDomain(HYDRO_, pore_size_index)    ! in
    dsl4jsb_Get_var2D_onDomain(HYDRO_, vol_field_cap)      ! in
    dsl4jsb_Get_var2D_onDomain(HYDRO_, vol_p_wilt)         ! in
    dsl4jsb_Get_var2D_onDomain(HYDRO_, vol_wres)           ! in
    IF (dsl4jsb_Config(HYDRO_)%l_organic)   &
     &   dsl4jsb_Get_var3D_onDomain(HYDRO_, fract_org_sl)  ! in
    dsl4jsb_Get_var3D_onDomain(HYDRO_, hyd_cond_sat_sl)    ! out
    dsl4jsb_Get_var3D_onDomain(HYDRO_, vol_porosity_sl)    ! out
    dsl4jsb_Get_var3D_onDomain(HYDRO_, bclapp_sl)          ! out
    dsl4jsb_Get_var3D_onDomain(HYDRO_, matric_pot_sl)      ! out
    dsl4jsb_Get_var3D_onDomain(HYDRO_, pore_size_index_sl) ! out
    dsl4jsb_Get_var3D_onDomain(HYDRO_, vol_field_cap_sl)   ! out
    dsl4jsb_Get_var3D_onDomain(HYDRO_, vol_p_wilt_sl)      ! out
    dsl4jsb_Get_var3D_onDomain(HYDRO_, vol_wres_sl)        ! out

    ! -------------------------------------------------------------------------------------------------- !
    !
    !> Distribute 2d input parameters over soil layers
    !
    hyd_cond_sat_sl   (:,:,:) = SPREAD(hyd_cond_sat   (:,:), DIM=2, ncopies=nsoil)
    vol_porosity_sl   (:,:,:) = SPREAD(vol_porosity   (:,:), DIM=2, ncopies=nsoil)
    bclapp_sl         (:,:,:) = SPREAD(bclapp         (:,:), DIM=2, ncopies=nsoil)
    matric_pot_sl     (:,:,:) = SPREAD(matric_pot     (:,:), DIM=2, ncopies=nsoil)
    pore_size_index_sl(:,:,:) = SPREAD(pore_size_index(:,:), DIM=2, ncopies=nsoil)
    vol_field_cap_sl  (:,:,:) = SPREAD(vol_field_cap  (:,:), DIM=2, ncopies=nsoil)
    vol_p_wilt_sl     (:,:,:) = SPREAD(vol_p_wilt     (:,:), DIM=2, ncopies=nsoil)
    vol_wres_sl       (:,:,:) = SPREAD(vol_wres       (:,:), DIM=2, ncopies=nsoil)

    ! -------------------------------------------------------------------------------------------------- !
    !
    !> Update soil parameters with organic fractions for deep and top soil layers
    !  Attention: hyd_cond_sat_sl is calculated differently (correct) in update_soil_properties
    !
    IF (dsl4jsb_Config(HYDRO_)%l_organic) THEN
      DO i=1,nsoil
        IF (i == 1) THEN
          hyd_cond_sat_org    = hyd_cond_sat_org_top
          vol_porosity_org    = vol_porosity_org_top
          bclapp_org          = bclapp_org_top
          matric_pot_org      = matric_pot_org_top
          pore_size_index_org = pore_size_index_org_top
          vol_field_cap_org   = vol_field_cap_org_top
          vol_p_wilt_org      = vol_p_wilt_org_top
          vol_wres_org        = vol_wres_org_top
        ELSE
          hyd_cond_sat_org    = hyd_cond_sat_org_below
          vol_porosity_org    = vol_porosity_org_below
          bclapp_org          = bclapp_org_below
          matric_pot_org      = matric_pot_org_below
          pore_size_index_org = pore_size_index_org_below
          vol_field_cap_org   = vol_field_cap_org_below
          vol_p_wilt_org      = vol_p_wilt_org_below
          vol_wres_org        = vol_wres_org_below
        END IF

        ! DOCU
        hyd_cond_sat_sl   (:,i,:) = (1 - fract_org_sl(:,i,:)) * hyd_cond_sat_sl   (:,i,:) &
          & + fract_org_sl(:,i,:) * hyd_cond_sat_org
        vol_porosity_sl   (:,i,:) = (1 - fract_org_sl(:,i,:)) * vol_porosity_sl   (:,i,:) &
          & + fract_org_sl(:,i,:) * vol_porosity_org
        bclapp_sl         (:,i,:) = (1 - fract_org_sl(:,i,:)) * bclapp_sl         (:,i,:) &
          & + fract_org_sl(:,i,:) * bclapp_org
        matric_pot_sl     (:,i,:) = (1 - fract_org_sl(:,i,:)) * matric_pot_sl     (:,i,:) &
          & + fract_org_sl(:,i,:) * matric_pot_org
        pore_size_index_sl(:,i,:) = (1 - fract_org_sl(:,i,:)) * pore_size_index_sl(:,i,:) &
          & + fract_org_sl(:,i,:) * pore_size_index_org
        vol_field_cap_sl  (:,i,:) = (1 - fract_org_sl(:,i,:)) * vol_field_cap_sl  (:,i,:) &
          & + fract_org_sl(:,i,:) * vol_field_cap_org
        vol_p_wilt_sl     (:,i,:) = (1 - fract_org_sl(:,i,:)) * vol_p_wilt_sl     (:,i,:) &
          & + fract_org_sl(:,i,:) * vol_p_wilt_org
        vol_wres_sl       (:,i,:) = (1 - fract_org_sl(:,i,:)) * vol_wres_sl       (:,i,:) &
          & + fract_org_sl(:,i,:) * vol_wres_org
      END DO
    END IF

    !$ACC UPDATE ASYNC(1) &
    !$ACC   DEVICE(hyd_cond_sat_sl, vol_porosity_sl, bclapp_sl)         &
    !$ACC   DEVICE(matric_pot_sl, pore_size_index_sl, vol_field_cap_sl) &
    !$ACC   DEVICE(vol_p_wilt_sl, vol_wres_sl)

  END SUBROUTINE init_soil_properties

  ! ====================================================================================================== !
  !
  !> Sanitize hydro variables.
  !>
  !> Loading a model state from a file with lower than double precision can lead to model aborts because
  !> the state violates model constraints due to rounding errors in soil water content and root depth.
  !> This routine sanitizes the current state by limiting soil water and ice to the field capacity, and by
  !> recomputing the root depth.
  !>
  SUBROUTINE hydro_sanitize_state(tile)

    USE mo_util, ONLY: soil_depth_to_layers_2d

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile

    CHARACTER(len=*), PARAMETER :: routine = modname//':hydro_sanitize_state'

    REAL(wp), PARAMETER :: warn_level_wpi_soil_excess = 1.05_wp
                           !< If water content exceeds field capacity by more than this factor prior
                           !< to sanitization, a warning is issued.
    REAL(wp), PARAMETER :: zfcmin = 1.e-10_wp
                           !< minimum field capacity (or rather epsilon?)

    INTEGER :: maxloc_wpi_soil_excess(3)
    TYPE(t_jsb_vgrid), POINTER :: soil_w

    dsl4jsb_Def_memory(HYDRO_) !< memory of `HYDRO_` process

    dsl4jsb_Real2D_onDomain :: &
      & root_depth,            &
      & soil_depth

    dsl4jsb_Real3D_onDomain :: &
      & root_depth_sl,         &
      & soil_depth_sl,         &
      & vol_porosity_sl,       &
      & wtr_soil_sat_sl,       &
      & wtr_soil_sl,           &
      & ice_soil_sl

    IF (.NOT. tile%Is_process_active(HYDRO_)) RETURN

    dsl4jsb_Get_memory(HYDRO_)

    dsl4jsb_Get_var2D_onDomain(HYDRO_, root_depth)
    dsl4jsb_Get_var2D_onDomain(HYDRO_, soil_depth)

    dsl4jsb_Get_var3D_onDomain(HYDRO_, root_depth_sl)
    dsl4jsb_Get_var3D_onDomain(HYDRO_, soil_depth_sl)
    dsl4jsb_Get_var3D_onDomain(HYDRO_, vol_porosity_sl)
    dsl4jsb_Get_var3D_onDomain(HYDRO_, wtr_soil_sat_sl)
    dsl4jsb_Get_var3D_onDomain(HYDRO_, wtr_soil_sl)
    dsl4jsb_Get_var3D_onDomain(HYDRO_, ice_soil_sl)

    IF (tile%contains_soil) THEN

      ! Issue a warning if current soil water anywherer exceeds saturation capacity by more than a factor of
      ! `warn_level_wpi_soil_excess`.
      IF (ANY(wtr_soil_sl(:,:,:) > warn_level_wpi_soil_excess * wtr_soil_sat_sl(:,:,:))) THEN
        maxloc_wpi_soil_excess(:) = MAXLOC(wtr_soil_sl(:,:,:) - wtr_soil_sat_sl(:,:,:))
        WRITE (message_text,*) 'liquid water content above saturation capacity at ', maxloc_wpi_soil_excess(:), &
          & ': wtr_soil_sl = ',                                                                                 &
          & wtr_soil_sl    (maxloc_wpi_soil_excess(1), maxloc_wpi_soil_excess(2), maxloc_wpi_soil_excess(3)),   &
          & ', wtr_soil_sat_sl = ',                                                                             &
          & wtr_soil_sat_sl(maxloc_wpi_soil_excess(1), maxloc_wpi_soil_excess(2), maxloc_wpi_soil_excess(3))
        CALL warning(routine, message_text)
      END IF
      ! Limit soil water and soil ice to saturation amounts.
      IF (ANY((wtr_soil_sl(:,:,:) + ice_soil_sl(:,:,:))                &
        &  > MAX(vol_porosity_sl (:,:,:) * soil_depth_sl(:,:,:), 0._wp))) THEN
        ! Limit liquid soil water to saturation amount
        wtr_soil_sl(:,:,:) = MIN(wtr_soil_sl(:,:,:), wtr_soil_sat_sl(:,:,:))

        ! Limit ice content such that soil water + soil ice <= saturation.
        ice_soil_sl(:,:,:) = MAX(0._wp, MIN(ice_soil_sl(:,:,:), &
          &                             MAX(vol_porosity_sl (:,:,:) * soil_depth_sl(:,:,:), 0._wp) - wtr_soil_sl(:,:,:)))

        !$ACC UPDATE DEVICE(wtr_soil_sl, ice_soil_sl) ASYNC(1)
      END IF
    END IF

    IF (tile%contains_vegetation) THEN
      ! Due to rounding, the root depth in a layer might differ from the soil depth even though the layer is contained
      ! in the root zone. The might not add up to the total root depth for the same reason. Fix by recomputing.
      soil_w => Get_vgrid('soil_depth_water')
      root_depth(:,:) = MIN(root_depth(:,:), soil_depth(:,:))
      root_depth_sl(:,:,:) = soil_depth_to_layers_2d(root_depth(:,:), soil_w%dz(:))

      !$ACC UPDATE DEVICE(root_depth, root_depth_sl) ASYNC(1)
    END IF

  END SUBROUTINE hydro_sanitize_state
  !
  ! ====================================================================================================== !

#endif
END MODULE mo_hydro_init
