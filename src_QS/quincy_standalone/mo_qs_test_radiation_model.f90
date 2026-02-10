!> Routine for running the test-radiation QUINCY model
!>
!> QUINCY standalone
!>
!> ---------------------------------------
!> Copyright (C) 2013-2026, MPI-M, MPI-BGC
!>
!> Contact: https://doi.org/10.17871/quincy-model-2019
!> Authors: AUTHORS.md
!> See LICENSES/ for license information
!> SPDX-License-Identifier: BSD-3-Clause
!> ---------------------------------------
!>
MODULE mo_qs_test_radiation_model

  USE mo_kind,                    ONLY: wp

  USE mo_lnd_bgcm_idx
  USE mo_lnd_bgcm_store,          ONLY: t_lnd_bgcm_store
  USE mo_lnd_bgcm_store_class,    ONLY: VEG_BGCM_POOL_ID

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: run_test_radiation_model

  CHARACTER(len=*), PARAMETER :: modname = 'mo_qs_test_radiation_model'

CONTAINS

  ! ======================================================================================================= !
  !>Run the test_radiation model from the land model interface
  !>
  !>  run one timestep
  !> 
  !>  process inits: assimi_init, rad_init, veg_init
  !>
  SUBROUTINE run_test_radiation_model(tile, options)
    USE mo_jsb_class,                     ONLY: Get_model
    USE mo_jsb_tile_class,                ONLY: t_jsb_tile_abstract
    USE mo_jsb_task_class,                ONLY: t_jsb_task_options
    USE mo_jsb_model_class,               ONLY: t_jsb_model
    USE mo_jsb_lctlib_class,              ONLY: t_lctlib_element
    USE mo_jsb_process_class,             ONLY: Q_ASSIMI_, RAD_, VEG_, HYDRO_, TURB_, A2L_
    USE mo_jsb_grid_class,                ONLY: t_jsb_vgrid
    USE mo_jsb_grid,                      ONLY: Get_vgrid
    USE mo_atmland_constants,             ONLY: frac_vis_swtot_srf, frac_par_swvis_srf, eps_vpd, &
      &                                         standard_press_srf, def_co2_mixing_ratio
    USE mo_phy_schemes,                   ONLY: qsat_mixed
    USE mo_jsb_physical_constants,        ONLY: molar_mass_N, Tzero, grav
    USE mo_q_veg_canopy,                  ONLY: calc_canopy_layers
    USE mo_qs_rad_scope,                     ONLY: calc_radiative_transfer
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Use_config(Q_ASSIMI_)
    dsl4jsb_Use_memory(RAD_)
    dsl4jsb_Use_memory(VEG_)
    dsl4jsb_Use_memory(HYDRO_)
    dsl4jsb_Use_memory(TURB_)
    dsl4jsb_Use_memory(A2L_)
    ! ----------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_tile_abstract), INTENT(inout)     :: tile         !< one tile with data structure for one lct
    TYPE(t_jsb_task_options),   INTENT(in)        :: options      !< model options
    ! ----------------------------------------------------------------------------------------------------- !
    TYPE(t_jsb_model),      POINTER  :: model             !< the model
    TYPE(t_lnd_bgcm_store), POINTER  :: bgcm_store        !< the bgcm store of this tile
    TYPE(t_lctlib_element), POINTER  :: lctlib           !< land-cover-type library - parameter across pft's
    TYPE(t_jsb_vgrid),      POINTER  :: vgrid_canopy     !< Vertical grid
    TYPE(t_jsb_vgrid),      POINTER  :: vgrid_spectral   !< Vertical grid used to display the sprectal bands of light
    INTEGER                          :: ncanopy          !< number of canopy layers
    INTEGER                          :: nspec            !< number of spectral bands
    REAL(wp), DIMENSION(options%nc)  :: hlp1
    INTEGER             :: i, lun           !< looping, logical unit number
    INTEGER             :: ref_pft          !< PFT identifier
    INTEGER, PARAMETER  :: npft  = 8        !< number of PFTs simulated by the test_canopy model
    ! parameters (specific to the test_canopy model)
    REAL(wp),SAVE :: ref_t_air = 25.0_wp              ! temperature in deg C
    REAL(wp),SAVE :: ref_vpd = 1000._wp               ! vpd in Pa
    REAL(wp),SAVE :: ref_swdown = 652.173913_wp       ! PAR in W m-2
    REAL(wp),SAVE :: ref_fdiffuse = 0.2_wp            ! diffuse fraction
    REAL(wp),SAVE :: ref_cos_angle = 0.8_wp           ! cosine of zenith angle
    REAL(wp),SAVE :: ref_co2_mixing_ratio = 380.0_wp  ! co2 in ppm
    REAL(wp),SAVE :: ref_wind = 1.0_wp                ! wind in m/s
    REAL(wp),SAVE :: ref_ga = 0.4_wp                  ! aerodynamic conductance in m/s
    REAL(wp),SAVE :: ref_theta = 1._wp                ! soil water availability in fraction of AWC
    REAL(wp),SAVE :: ref_lai = 3._wp                  ! LAI
    REAL(wp),SAVE :: ref_nleaf = 2.0_wp               ! N leaf per unit leaf area in gN/m2
    REAL(wp)                :: lctlib_g1              !< set to g1_medlyn or g1_bberry depending on canopy_conductance_scheme
    REAL(wp)                :: alpha                  !< implicitness factor
    INTEGER                 :: ipoi, ichunk           ! number of the chunk: ipoi 1:nc / ichunk ics:ice
    INTEGER                 :: iblk, ics, ice, nc
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':run_test_radiation_model'
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Def_mt2L2D :: veg_pool_mt
    ! ----------------------------------------------------------------------------------------------------- !
    ! Declare process configuration and memory Pointers
    dsl4jsb_Def_config(Q_ASSIMI_)
    dsl4jsb_Def_memory(RAD_)
    dsl4jsb_Def_memory(VEG_)
    dsl4jsb_Def_memory(HYDRO_)
    dsl4jsb_Def_memory(TURB_)
    dsl4jsb_Def_memory(A2L_)
    ! ----------------------------------------------------------------------------------------------------- !
    ! RAD_
    dsl4jsb_Real3D_onChunk      :: ppfd_sunlit_tfrac_mavg_cl
    dsl4jsb_Real3D_onChunk      :: ppfd_shaded_tfrac_mavg_cl
    dsl4jsb_Real3D_onChunk      :: canopy_reflectance_sp
    ! HYDRO_
    dsl4jsb_Real2D_onChunk      :: wtr_rootzone_rel
    ! TURB_
    dsl4jsb_Real2D_onChunk      :: drag_srf_stdalone
    ! A2L_
    dsl4jsb_Real2D_onChunk      :: t_air
    dsl4jsb_Real2D_onChunk      :: q_air
    dsl4jsb_Real2D_onChunk      :: swpar_srf_down
    dsl4jsb_Real2D_onChunk      :: cos_zenith_angle
    dsl4jsb_Real2D_onChunk      :: co2_mixing_ratio
    dsl4jsb_Real2D_onChunk      :: wind_air
    dsl4jsb_Real2D_onChunk      :: wind_10m
    dsl4jsb_Real2D_onChunk      :: swvis_srf_down
    dsl4jsb_Real2D_onChunk      :: swnir_srf_down
    dsl4jsb_Real2D_onChunk      :: fract_par_diffuse
    dsl4jsb_Real2D_onChunk      :: press_srf
    ! VEG_
    dsl4jsb_Real2D_onChunk      :: lai
    dsl4jsb_Real2D_onChunk      :: t_jmax_opt_mavg
    dsl4jsb_Real3D_onChunk      :: leaf_nitrogen_cl
    dsl4jsb_Real3D_onChunk      :: fn_rub_cl
    dsl4jsb_Real3D_onChunk      :: fn_et_cl
    dsl4jsb_Real3D_onChunk      :: fn_pepc_cl
    dsl4jsb_Real3D_onChunk      :: fn_chl_cl
    dsl4jsb_Real3D_onChunk      :: fn_oth_cl
    dsl4jsb_Real3D_onChunk      :: fleaf_sunlit_tfrac_mavg_cl
    dsl4jsb_Real3D_onChunk      :: lai_cl
    dsl4jsb_Real3D_onChunk      :: cumm_lai_cl
    ! ----------------------------------------------------------------------------------------------------- !
    iblk    = options%iblk
    ics     = options%ics
    ice     = options%ice
    nc      = options%nc
    alpha   = options%alpha
    ! ----------------------------------------------------------------------------------------------------- !
    ! if pftid is set in the namelist use this pftid
    IF(tile%lcts(1)%lib_id>0)THEN
      ref_pft = tile%lcts(1)%lib_id
    ! else use pftid 1
    ELSE
      ref_pft = 0
      tile%lcts(1)%lib_id = 1
    ENDIF
    model           => Get_model(tile%owner_model_id)
    lctlib          => model%lctlib(tile%lcts(1)%lib_id)
    vgrid_canopy    => Get_vgrid('q_canopy_layer')
    vgrid_spectral  => Get_vgrid('spectral_bands')
    ncanopy         =  vgrid_canopy%n_levels
    nspec           =  vgrid_spectral%n_levels
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Get_config(Q_ASSIMI_)
    dsl4jsb_Get_memory(RAD_)
    dsl4jsb_Get_memory(VEG_)
    dsl4jsb_Get_memory(HYDRO_)
    dsl4jsb_Get_memory(TURB_)
    dsl4jsb_Get_memory(A2L_)
    ! ----------------------------------------------------------------------------------------------------- !
    bgcm_store => tile%bgcm_store
    dsl4jsb_Get_mt2L2D(VEG_BGCM_POOL_ID, veg_pool_mt)
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Get_var3D_onChunk(RAD_, ppfd_sunlit_tfrac_mavg_cl)      ! in
    dsl4jsb_Get_var3D_onChunk(RAD_, ppfd_shaded_tfrac_mavg_cl)      ! in
    dsl4jsb_Get_var3D_onChunk(RAD_, canopy_reflectance_sp)
    ! ---------------------------
    dsl4jsb_Get_var2D_onChunk(HYDRO_, wtr_rootzone_rel)
    ! ---------------------------
    dsl4jsb_Get_var2D_onChunk(TURB_, drag_srf_stdalone)
    ! ---------------------------
    dsl4jsb_Get_var2D_onChunk(A2L_, t_air)
    dsl4jsb_Get_var2D_onChunk(A2L_, q_air)
    dsl4jsb_Get_var2D_onChunk(A2L_, swpar_srf_down)
    dsl4jsb_Get_var2D_onChunk(A2L_, cos_zenith_angle)
    dsl4jsb_Get_var2D_onChunk(A2L_, co2_mixing_ratio)
    dsl4jsb_Get_var2D_onChunk(A2L_, wind_air)
    dsl4jsb_Get_var2D_onChunk(A2L_, wind_10m)
    dsl4jsb_Get_var2D_onChunk(A2L_, swvis_srf_down)
    dsl4jsb_Get_var2D_onChunk(A2L_, swnir_srf_down)
    dsl4jsb_Get_var2D_onChunk(A2L_, fract_par_diffuse)
    dsl4jsb_Get_var2D_onChunk(A2L_, press_srf)
    ! ---------------------------
    dsl4jsb_Get_var2D_onChunk(VEG_, lai)
    dsl4jsb_Get_var2D_onChunk(VEG_, t_jmax_opt_mavg)
    dsl4jsb_Get_var3D_onChunk(VEG_, leaf_nitrogen_cl)
    dsl4jsb_Get_var3D_onChunk(VEG_, fn_rub_cl)
    dsl4jsb_Get_var3D_onChunk(VEG_, fn_et_cl)
    dsl4jsb_Get_var3D_onChunk(VEG_, fn_pepc_cl)
    dsl4jsb_Get_var3D_onChunk(VEG_, fn_chl_cl)
    dsl4jsb_Get_var3D_onChunk(VEG_, fn_oth_cl)
    dsl4jsb_Get_var3D_onChunk(VEG_, fleaf_sunlit_tfrac_mavg_cl)     ! in
    dsl4jsb_Get_var3D_onChunk(VEG_, lai_cl)
    dsl4jsb_Get_var3D_onChunk(VEG_, cumm_lai_cl)
    ! ----------------------------------------------------------------------------------------------------- !

    !>0.9 define chunk for output: ipoi 1:nc / ichunk ics:ice
    !>
    ipoi   = 1
    ichunk = ics+ipoi-1

    !>1.0 Init
    !>

    !>  1.1 set g1 according to canopy_conductance_scheme
    !>
    SELECT CASE(TRIM(dsl4jsb_Config(Q_ASSIMI_)%canopy_conductance_scheme))
      CASE ("medlyn")
        lctlib_g1 = lctlib%g1_medlyn
      CASE ("ballberry")
        lctlib_g1 = lctlib%g1_bberry
    END SELECT

    !>  1.2 Initalisation of met-forcing which is assumed constant
    !>
    t_air(:)            = ref_t_air + TZero
    press_srf(:)        = standard_press_srf
    q_air(:)            = qsat_mixed(t_air(:),press_srf(:)) - ref_vpd / press_srf(:) * eps_vpd
    co2_mixing_ratio(:) = def_co2_mixing_ratio   ! def_co2_mixing_ratio is used instead of ref_co2_mixing_ratio
    swpar_srf_down(:)   = ref_swdown
    swvis_srf_down(:)   = swpar_srf_down(:) / frac_par_swvis_srf
    swnir_srf_down(:)   = swvis_srf_down(:) * (1._wp - frac_vis_swtot_srf) / frac_vis_swtot_srf
    fract_par_diffuse(:)= ref_fdiffuse
    wind_air(:)         = ref_wind
    wind_10m(:)         = wind_air(:)
    drag_srf_stdalone(:)     = ref_ga / wind_air(:) * alpha * grav * model%options%dtime
    cos_zenith_angle(:) = ref_cos_angle
    !latitude                  = 50.0_wp
    !longitude                 = 0.0_wp

    !> 1.3 Initialisation of soil and vegetation (so far relevant)
    !!
    !!
    ! soil water & lai 
    wtr_rootzone_rel(:) = ref_theta
    lai(:)   = ref_lai
    ! leaf N & P
    IF(ref_pft==0) THEN
      veg_pool_mt(ix_leaf, ixN, :) = ref_nleaf / molar_mass_N * lai(:)
    ELSE
      veg_pool_mt(ix_leaf, ixN, :) = lai(:) / lctlib%sla / lctlib%cn_leaf
      veg_pool_mt(ix_leaf, ixP, :) = veg_pool_mt(ix_leaf, ixN, :) / lctlib%np_leaf
    ENDIF
    ! canopy
    CALL calc_canopy_layers(nc                                                      , & ! in
                            ncanopy                                                 , &
                            model%options%dtime                                     , &
                            vgrid_canopy%dz(:)                                      , &
                            vgrid_canopy%lbounds(:)                                 , &
                            vgrid_canopy%ubounds(:)                                 , &
                            lctlib%ps_pathway                                       , &
                            lctlib%k0_fn_struc                                      , &
                            lctlib%fn_oth_min                                       , &
                            lctlib%sla                                              , &
                            lctlib%np_leaf                                          , &
                            lctlib%gmin                                             , &
                            lctlib%g0                                               , &
                            lctlib_g1                                               , &
                            lctlib%t_jmax_omega                                     , &
                            dsl4jsb_Config(Q_ASSIMI_)%flag_optimal_Nfraction          , &
                            TRIM(dsl4jsb_Config(Q_ASSIMI_)%canopy_conductance_scheme) , &
                            veg_pool_mt(ix_leaf, ixN, :)                            , &
                            lai(:)                                                  , & ! in
                            ppfd_sunlit_tfrac_mavg_cl(:,:), &                           ! inout
                            ppfd_shaded_tfrac_mavg_cl(:,:), & 
                            fleaf_sunlit_tfrac_mavg_cl(:,:), &
                            fn_rub_cl(:,:), &
                            fn_et_cl(:,:), &
                            fn_pepc_cl(:,:), &
                            fn_chl_cl(:,:), &
                            fn_oth_cl(:,:), &                                           ! inout
                            lai_cl(:,:), &                                              ! out
                            cumm_lai_cl(:,:), &
                            leaf_nitrogen_cl(:,:), &                                    ! out
                            inquire_n_fractions = .TRUE.)                               ! optional in
    ! t_jmax_opt_mavg
    t_jmax_opt_mavg(:) = 20._wp+0.35_wp*(t_air(:) - Tzero)


    !>2.0 radiative transfer calculations
    !>
    CALL calc_radiative_transfer( &
      &                          nc                                              , & ! in
      &                          ncanopy                                         , &
      &                          nspec                                           , &
      &                          t_air(:)                                        , &
      &                          canopy_reflectance_sp(:,:)                        ) ! inout

  END SUBROUTINE run_test_radiation_model

END MODULE mo_qs_test_radiation_model
