!> Routine for running the test-canopy QUINCY model
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
MODULE mo_qs_test_canopy_model

  USE mo_kind,                    ONLY: wp

  USE mo_lnd_bgcm_idx
  USE mo_lnd_bgcm_store,          ONLY: t_lnd_bgcm_store
  USE mo_lnd_bgcm_store_class,    ONLY: VEG_BGCM_POOL_ID

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: run_test_canopy_model

  CHARACTER(len=*), PARAMETER :: modname = 'mo_qs_test_canopy_model'

CONTAINS

  ! ======================================================================================================= !
  !>Run the test_canopy model from the land model interface
  !>
  !>  run one timestep
  !>
  SUBROUTINE run_test_canopy_model(tile, options)
    USE mo_jsb_class,              ONLY: Get_model
    USE mo_jsb_tile_class,         ONLY: t_jsb_tile_abstract
    USE mo_jsb_task_class,         ONLY: t_jsb_task_options
    USE mo_jsb_model_class,        ONLY: t_jsb_model
    USE mo_jsb_lctlib_class,       ONLY: t_lctlib_element
    USE mo_jsb_process_class,      ONLY: Q_ASSIMI_, RAD_, VEG_, TURB_, A2L_
    USE mo_jsb_grid_class,         ONLY: t_jsb_vgrid
    USE mo_jsb_grid,               ONLY: Get_vgrid
    USE mo_rad_parameters,         ONLY: rad2ppfd
    USE mo_atmland_constants,      ONLY: frac_vis_swtot_srf, frac_par_swvis_srf,eps_vpd, standard_press_srf, def_co2_mixing_ratio
    USE mo_phy_schemes,            ONLY: qsat_mixed
    USE mo_veg_util,               ONLY: reset_veg_fluxes
    USE mo_jsb_physical_constants, ONLY: molar_mass_N, Tzero, grav
    USE mo_q_veg_canopy            ! .. e.g. calc_canopy_layers()
    USE mo_rad_interface,          ONLY: update_radiation_par
    USE mo_q_assimi_interface,     ONLY: update_canopy_fluxes, update_reset_q_assimi_fluxes
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Use_config(Q_ASSIMI_)
    dsl4jsb_Use_memory(Q_ASSIMI_)
    dsl4jsb_Use_memory(RAD_)
    dsl4jsb_Use_memory(VEG_)
    dsl4jsb_Use_memory(TURB_)
    dsl4jsb_Use_memory(A2L_)
    ! ----------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_tile_abstract), INTENT(inout)     :: tile         !< one tile with data structure for one lct
    TYPE(t_jsb_task_options),   INTENT(in)        :: options      !< model options
    ! ----------------------------------------------------------------------------------------------------- !
    TYPE(t_jsb_model),      POINTER  :: model             !< the model
    TYPE(t_lnd_bgcm_store), POINTER  :: bgcm_store        !< the bgcm store of this tile
    TYPE(t_lctlib_element), POINTER  :: lctlib            !< land-cover-type library - parameter across pft's
    TYPE(t_jsb_vgrid),      POINTER  :: vgrid_canopy      !< Vertical grid
    INTEGER                          :: ncanopy           !< number of canopy layers
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
    REAL(wp),SAVE :: ref_soil_root_pot = 0._wp        ! soil water potential, ref for test_canopy model
    REAL(wp),SAVE :: ref_lai = 3._wp                  ! LAI
    REAL(wp),SAVE :: ref_nleaf = 2.0_wp               ! N leaf per unit leaf area in gN/m2
    REAL(wp)                :: lctlib_g1              !< set to g1_medlyn or g1_bberry depending on canopy_conductance_scheme
    REAL(wp)                :: alpha                  !< implicitness factor
    INTEGER                 :: ipoi, ichunk           ! number of the chunk: ipoi 1:nc / ichunk ics:ice
    INTEGER                 :: iblk, ics, ice, nc
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':run_test_canopy_model'
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Def_mt2L2D :: veg_pool_mt
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Def_config(Q_ASSIMI_)
    dsl4jsb_Def_memory(Q_ASSIMI_)
    dsl4jsb_Def_memory(RAD_)
    dsl4jsb_Def_memory(VEG_)
    dsl4jsb_Def_memory(TURB_)
    dsl4jsb_Def_memory(A2L_)
    ! ----------------------------------------------------------------------------------------------------- !
    ! Q_ASSIMI_
    dsl4jsb_Real2D_onChunk      :: net_assimilation 
    dsl4jsb_Real2D_onChunk      :: canopy_cond
    dsl4jsb_Real2D_onChunk      :: co2_conc_leaf
    dsl4jsb_Real2D_onChunk      :: soa_tsoa_mavg
    dsl4jsb_Real2D_onChunk      :: vpd
    dsl4jsb_Real2D_onChunk      :: wtr_soil_root_pot
    ! RAD_
    dsl4jsb_Real3D_onChunk      :: ppfd_sunlit_tfrac_mavg_cl
    dsl4jsb_Real3D_onChunk      :: ppfd_shaded_tfrac_mavg_cl
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
    ! TURB_
    dsl4jsb_Real2D_onChunk      :: drag_srf_stdalone
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
    model  => Get_model(tile%owner_model_id)
    lctlib        => model%lctlib(tile%lcts(1)%lib_id)
    vgrid_canopy  => Get_vgrid('q_canopy_layer')
    ncanopy       =  vgrid_canopy%n_levels
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Get_config(Q_ASSIMI_)
    dsl4jsb_Get_memory(Q_ASSIMI_)
    dsl4jsb_Get_memory(RAD_)
    dsl4jsb_Get_memory(VEG_)
    dsl4jsb_Get_memory(TURB_)
    dsl4jsb_Get_memory(A2L_)
    ! ----------------------------------------------------------------------------------------------------- !
    bgcm_store => tile%bgcm_store
    dsl4jsb_Get_mt2L2D(VEG_BGCM_POOL_ID, veg_pool_mt)
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Get_var2D_onChunk(Q_ASSIMI_, net_assimilation)
    dsl4jsb_Get_var2D_onChunk(Q_ASSIMI_, canopy_cond)
    dsl4jsb_Get_var2D_onChunk(Q_ASSIMI_, co2_conc_leaf)
    dsl4jsb_Get_var2D_onChunk(Q_ASSIMI_, soa_tsoa_mavg)
    dsl4jsb_Get_var2D_onChunk(Q_ASSIMI_, vpd)
    dsl4jsb_Get_var2D_onChunk(Q_ASSIMI_, wtr_soil_root_pot)
    ! ---------------------------
    dsl4jsb_Get_var3D_onChunk(RAD_, ppfd_sunlit_tfrac_mavg_cl)      ! in
    dsl4jsb_Get_var3D_onChunk(RAD_, ppfd_shaded_tfrac_mavg_cl)      ! in
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
    soa_tsoa_mavg(:)    = 7.0_wp        ! is initialized in assimi init with -10.0_wp as the model starts in winter
    !latitude                  = 50.0_wp
    !longitude                 = 0.0_wp

    !> 1.2 Initialisation of soil and vegetation (so far relevant)
    !!
    !!
    lai(:)   = ref_lai
    IF(ref_pft==0) THEN
      veg_pool_mt(ix_leaf, ixN, :) = ref_nleaf / molar_mass_N * lai(:)
    ELSE
      veg_pool_mt(ix_leaf, ixN, :) = lai(:) / lctlib%sla / lctlib%cn_leaf
      veg_pool_mt(ix_leaf, ixP, :) = veg_pool_mt(ix_leaf, ixN, :) / lctlib%np_leaf
    ENDIF
    CALL calc_canopy_layers(nc                                                        , & ! in
                            ncanopy                                                   , &
                            model%options%dtime                                       , &
                            vgrid_canopy%dz(:)                                        , &
                            vgrid_canopy%lbounds(:)                                   , &
                            vgrid_canopy%ubounds(:)                                   , &
                            lctlib%ps_pathway                                         , &
                            lctlib%k0_fn_struc                                        , &
                            lctlib%fn_oth_min                                         , &
                            lctlib%sla                                                , &
                            lctlib%np_leaf                                            , &
                            lctlib%gmin                                               , &
                            lctlib%g0                                                 , &
                            lctlib_g1                                                 , &
                            lctlib%t_jmax_omega                                       , &
                            dsl4jsb_Config(Q_ASSIMI_)%flag_optimal_Nfraction            , &
                            TRIM(dsl4jsb_Config(Q_ASSIMI_)%canopy_conductance_scheme)   , &
                            leaf_nitrogen       = veg_pool_mt(ix_leaf, ixN, :)        , &
                            lai                 = lai(:)                              , & ! in
                            ppfd_sunlit_cl      = ppfd_sunlit_tfrac_mavg_cl(:,:)      , & ! inout
                            ppfd_shaded_cl      = ppfd_shaded_tfrac_mavg_cl(:,:)      , & 
                            fleaf_sunlit_cl     = fleaf_sunlit_tfrac_mavg_cl(:,:)     , &
                            fn_rub_cl           = fn_rub_cl(:,:), &
                            fn_et_cl            = fn_et_cl(:,:), &
                            fn_pepc_cl          = fn_pepc_cl(:,:), &
                            fn_chl_cl           = fn_chl_cl(:,:), &
                            fn_oth_cl           = fn_oth_cl(:,:), &                       ! inout
                            lai_cl              = lai_cl(:,:), &                          ! out
                            cumm_lai_cl         = cumm_lai_cl(:,:), &
                            leaf_nitrogen_cl    = leaf_nitrogen_cl(:,:), &                ! out   
                            inquire_n_fractions = .TRUE.)                                 ! optional in
    t_jmax_opt_mavg(:) = 20._wp+0.35_wp*(t_air(:) - Tzero)
  
    !> 2.0 Sensitivity plots
    !!
    !!
  
    !> 2.1 PAR (0-1500 µmol m-2 s-1)
    !!
    !!
    lun = 10
    WRITE(lun,*) 'PAR A GC CI'
  
    swpar_srf_down(:) = 0.0_wp
    DO i = 1,16
  
       swpar_srf_down(:) = swpar_srf_down(:) + 100.0_wp / rad2ppfd
       swvis_srf_down(:) = swpar_srf_down(:) / frac_par_swvis_srf
       swnir_srf_down(:) = swvis_srf_down(:) * (1._wp - frac_vis_swtot_srf) / frac_vis_swtot_srf
  
       !! canopy calculation fluxes
       CALL reset_veg_fluxes(tile,options)             !< mo_veg_util
       CALL update_reset_q_assimi_fluxes(tile, model%options)
       CALL update_radiation_par(tile,options)
       CALL update_canopy_fluxes(tile,options)  !< mo_q_assimi_process
  
       WRITE(lun,'(4F16.7)') swpar_srf_down(ipoi) * rad2ppfd, &
             net_assimilation(ipoi), canopy_cond(ipoi) *1000._wp, &
             co2_conc_leaf(ipoi)/co2_mixing_ratio(ipoi) !,&
    ! net_assimilation(ipoi),
    ! co2_conc_leaf(ipoi)/co2_mixing_ratio(ipoi),veg_pool_mt(ix_leaf, ixN, ipoi) * molar_mass_N,&
    ! net_assimilation(ipoi)/(veg_pool_mt(ix_leaf, ixN, ipoi) * molar_mass_N)
  
    END DO
    swpar_srf_down(:) = ref_swdown
    swvis_srf_down(:) = swpar_srf_down(:) / frac_par_swvis_srf
    swnir_srf_down(:) = swvis_srf_down(:) * (1._wp - frac_vis_swtot_srf) / frac_vis_swtot_srf
  
    !> 2.2 Temperature (-5,60 degC)
    !!
    !!
    lun = lun + 1
    WRITE(lun,*) 'T A GC CI'
    t_air(:) = - 5.0_wp + Tzero
    hlp1(:)  = t_jmax_opt_mavg(:)
  
    DO i = 1,70
  
       t_air(:) = t_air(:) + 1.0_wp
       q_air(:) = 0.8_wp * qsat_mixed(t_air(:),press_srf(:))
       t_jmax_opt_mavg(:) = 20._wp+0.35_wp*(t_air(:) - Tzero)
  
       !! canopy calculation fluxes
       CALL reset_veg_fluxes(tile,options)             !< mo_veg_util
       CALL update_reset_q_assimi_fluxes(tile, model%options)
       CALL update_radiation_par(tile,options)
       CALL update_canopy_fluxes(tile,options)  !< mo_q_assimi_process
  
       WRITE(lun,'(5F16.7)') t_air(ipoi), &
             net_assimilation(ipoi), canopy_cond(ipoi) *1000._wp, &
             co2_conc_leaf(ipoi)/co2_mixing_ratio(ipoi)
  
    END DO
    t_air(:)    = ref_t_air + Tzero
    q_air(:)    = qsat_mixed(t_air(:),press_srf(:)) - ref_vpd / press_srf(:) * eps_vpd
    t_jmax_opt_mavg(:) = hlp1(:)
  
    !> 2.3 atmospheric CO2 (175 - 1000 ppm)
    !!
    !!
    lun = lun + 1
    WRITE(lun,*) 'CO2 A GC CI'
    co2_mixing_ratio(:) = 150._wp
    DO i = 1,34
       co2_mixing_ratio(:) = co2_mixing_ratio(:) + 25._wp
  
       !! canopy calculation fluxes
       CALL reset_veg_fluxes(tile,options)             !< mo_veg_util
       CALL update_reset_q_assimi_fluxes(tile, model%options)
       CALL update_radiation_par(tile,options)
       CALL update_canopy_fluxes(tile,options)  !< mo_q_assimi_process
  
       WRITE(lun,'(4F16.7)') co2_mixing_ratio(ipoi), &
             net_assimilation(ipoi), canopy_cond(ipoi) *1000._wp, &
             co2_conc_leaf(ipoi)/co2_mixing_ratio(ipoi)
  
    END DO
    co2_mixing_ratio(:) = ref_co2_mixing_ratio
  
    !> 2.4 VPD (0-3.3 kPa)
    !!
    !!
    lun = lun + 1
    WRITE(lun,*) 'VPD A GC CI'
    DO i = 1,30
       q_air(:) = qsat_mixed(t_air(:),press_srf(:)) - i * 100.0_wp / press_srf(:) * eps_vpd
  
       !! canopy calculation fluxes
       CALL reset_veg_fluxes(tile,options)             !< mo_veg_util
       CALL update_reset_q_assimi_fluxes(tile, model%options)
       CALL update_radiation_par(tile,options)
       CALL update_canopy_fluxes(tile,options)  !< mo_q_assimi_process
  
       vpd(:) = ( qsat_mixed(t_air(:),press_srf(:)) - q_air(:) ) / eps_vpd * press_srf(:) / 1000.0_wp
  
       WRITE(lun,'(4F16.7)') vpd(ipoi), &
             net_assimilation(ipoi), canopy_cond(ipoi) *1000._wp, &
             co2_conc_leaf(ipoi)/co2_mixing_ratio(ipoi)
  
    END DO
    q_air(:) = qsat_mixed(t_air(:),press_srf(:)) - ref_vpd / press_srf(:) * eps_vpd
    vpd(:) = ( qsat_mixed(t_air(:),press_srf(:)) - q_air(:) ) / eps_vpd * press_srf(:) / 1000.0_wp
  
    !> 2.5 soil water potential (0 to -4)
    !!
    !!
    lun = lun + 1
    WRITE(lun,*) 'SW_pot A GC CI'
    wtr_soil_root_pot(:) = ref_soil_root_pot
    DO i = 1,80
  
       wtr_soil_root_pot(:) = wtr_soil_root_pot(:) - 0.05_wp
  
       !! canopy calculation fluxes
       CALL reset_veg_fluxes(tile,options)             !< mo_veg_util
       CALL update_reset_q_assimi_fluxes(tile, model%options)
       CALL update_radiation_par(tile,options)
       CALL update_canopy_fluxes(tile,options)  !< mo_q_assimi_process
  
       WRITE(lun,'(4F16.7)') wtr_soil_root_pot(ipoi) *100._wp, &
             net_assimilation(ipoi), canopy_cond(ipoi) *1000._wp, &
             co2_conc_leaf(ipoi)/co2_mixing_ratio(ipoi)
  
    END DO
    wtr_soil_root_pot(:) = ref_soil_root_pot
  
    !> 2.6 LAI (increase from 0 to 11)
    !!
    !!
    lun = lun + 1
    WRITE(lun,*) 'LAI A GC CI'
    lai(:) = 0.0_wp
    DO i = 1,22
  
       lai(:)                       = lai(:) + 0.5_wp
       veg_pool_mt(ix_leaf, ixN, :) = ref_nleaf / molar_mass_N * lai(:)
       CALL calc_canopy_layers( nc                                             , &   ! in
                                ncanopy                                        , &
                                model%options%dtime                            , &
                                vgrid_canopy%dz(:)                             , &
                                vgrid_canopy%lbounds(:)                        , &
                                vgrid_canopy%ubounds(:)                        , &
                                lctlib%ps_pathway                              , &
                                lctlib%k0_fn_struc                             , &
                                lctlib%fn_oth_min                              , &
                                lctlib%sla                                     , &
                                lctlib%np_leaf                                 , &
                                lctlib%gmin                                    , &
                                lctlib%g0                                      , &
                                lctlib_g1                                      , &
                                lctlib%t_jmax_omega                            , &
                                dsl4jsb_Config(Q_ASSIMI_)%flag_optimal_Nfraction , &
                                TRIM(dsl4jsb_Config(Q_ASSIMI_)%canopy_conductance_scheme)       , &
                                leaf_nitrogen       = veg_pool_mt(ix_leaf, ixN, :)    , &
                                lai                 = lai(:)                          , & ! in
                                ppfd_sunlit_cl      = ppfd_sunlit_tfrac_mavg_cl(:,:)  , & ! inout
                                ppfd_shaded_cl      = ppfd_shaded_tfrac_mavg_cl(:,:)  , & 
                                fleaf_sunlit_cl     = fleaf_sunlit_tfrac_mavg_cl(:,:) , &
                                fn_rub_cl           = fn_rub_cl(:,:), &
                                fn_et_cl            = fn_et_cl(:,:), &
                                fn_pepc_cl          = fn_pepc_cl(:,:), &
                                fn_chl_cl           = fn_chl_cl(:,:), &
                                fn_oth_cl           = fn_oth_cl(:,:), &         ! inout
                                lai_cl              = lai_cl(:,:), &            ! out
                                cumm_lai_cl         = cumm_lai_cl(:,:), &
                                leaf_nitrogen_cl    = leaf_nitrogen_cl(:,:), &  ! out                                
                                inquire_n_fractions = .TRUE.)                     
  
       !! canopy calculation fluxes
       CALL reset_veg_fluxes(tile,options)             !< mo_veg_util
       CALL update_reset_q_assimi_fluxes(tile, model%options)
       CALL update_radiation_par(tile,options)
       CALL update_canopy_fluxes(tile,options)  !< mo_q_assimi_process
  
       WRITE(lun,'(4F16.7)') lai(ipoi), &
             net_assimilation(ipoi), canopy_cond(ipoi) *1000._wp, &
             co2_conc_leaf(ipoi)/co2_mixing_ratio(ipoi)
  
    END DO
    lai(:) = ref_lai
  
    !> 2.7 leaf N (0.5-3 g N / m2(LAI))
    !!
    !!
    lun = lun + 1
    WRITE(lun,*) 'Nleaf A GC CI'
    DO i = 1,24
  
       veg_pool_mt(ix_leaf, ixN, :) = 0.125_wp * (i) / molar_mass_N * lai(:)
       CALL calc_canopy_layers( nc                                             , &   ! in
                                ncanopy                                        , &
                                model%options%dtime                            , &
                                vgrid_canopy%dz(:)                             , &
                                vgrid_canopy%lbounds(:)                        , &
                                vgrid_canopy%ubounds(:)                        , &
                                lctlib%ps_pathway                              , &
                                lctlib%k0_fn_struc                             , &
                                lctlib%fn_oth_min                              , &
                                lctlib%sla                                     , &
                                lctlib%np_leaf                                 , &
                                lctlib%gmin                                    , &
                                lctlib%g0                                      , &
                                lctlib_g1                                      , &
                                lctlib%t_jmax_omega                            , &
                                dsl4jsb_Config(Q_ASSIMI_)%flag_optimal_Nfraction , &
                                TRIM(dsl4jsb_Config(Q_ASSIMI_)%canopy_conductance_scheme)       , &
                                leaf_nitrogen       = veg_pool_mt(ix_leaf, ixN, :)    , &
                                lai                 = lai(:)                          , & ! in
                                ppfd_sunlit_cl      = ppfd_sunlit_tfrac_mavg_cl(:,:)  , & ! inout
                                ppfd_shaded_cl      = ppfd_shaded_tfrac_mavg_cl(:,:)  , & 
                                fleaf_sunlit_cl     = fleaf_sunlit_tfrac_mavg_cl(:,:) , &
                                fn_rub_cl           = fn_rub_cl(:,:), &
                                fn_et_cl            = fn_et_cl(:,:), &
                                fn_pepc_cl          = fn_pepc_cl(:,:), &
                                fn_chl_cl           = fn_chl_cl(:,:), &
                                fn_oth_cl           = fn_oth_cl(:,:), &           ! inout
                                lai_cl              = lai_cl(:,:), &              ! out
                                cumm_lai_cl         = cumm_lai_cl(:,:), &
                                leaf_nitrogen_cl    = leaf_nitrogen_cl(:,:), &    ! out
                                inquire_n_fractions = .TRUE.)                     ! optional in
  
       !! canopy calculation fluxes
       CALL reset_veg_fluxes(tile,options)             !< mo_veg_util
       CALL update_reset_q_assimi_fluxes(tile, model%options)
       CALL update_radiation_par(tile,options)
       CALL update_canopy_fluxes(tile,options)  !< mo_q_assimi_process
  
       WRITE(lun,'(4F16.7)') veg_pool_mt(ix_leaf, ixN, ichunk) / lai(ipoi) * molar_mass_N, &
             net_assimilation(ipoi), canopy_cond(ipoi) * 1000._wp, &
             co2_conc_leaf(ipoi)/co2_mixing_ratio(ipoi)
  
    END DO
    veg_pool_mt(ix_leaf, ixN, :) = ref_nleaf / molar_mass_N * lai(:)
  
    !> 2.8 PFT
    !!
    !!
    lun = lun + 1
    WRITE(lun,*) 'PFT A GC CI'
    DO i = 1,npft
       veg_pool_mt(ix_leaf, ixN, :) = lai(:) /model%lctlib(i)%sla /model%lctlib(i)%cn_leaf
       veg_pool_mt(ix_leaf, ixP, :) = veg_pool_mt(ix_leaf, ixN, :) /model%lctlib(i)%np_leaf
       ! change lctlib id
       tile%lcts(1)%lib_id = i
       ! set pointer to the new PFT
       lctlib  => model%lctlib(tile%lcts(1)%lib_id)
       CALL calc_canopy_layers( nc                                             , &   ! in
                                ncanopy                                        , &
                                model%options%dtime                            , &
                                vgrid_canopy%dz(:)                             , &
                                vgrid_canopy%lbounds(:)                        , &
                                vgrid_canopy%ubounds(:)                        , &
                                lctlib%ps_pathway                              , &
                                lctlib%k0_fn_struc                             , &
                                lctlib%fn_oth_min                              , &
                                lctlib%sla                                     , &
                                lctlib%np_leaf                                 , &
                                lctlib%gmin                                    , &
                                lctlib%g0                                      , &
                                lctlib_g1                                      , &
                                lctlib%t_jmax_omega                            , &
                                dsl4jsb_Config(Q_ASSIMI_)%flag_optimal_Nfraction , &
                                TRIM(dsl4jsb_Config(Q_ASSIMI_)%canopy_conductance_scheme) , &
                                leaf_nitrogen       = veg_pool_mt(ix_leaf, ixN, :)      , &
                                lai                 = lai(:)                            , & ! in
                                ppfd_sunlit_cl      = ppfd_sunlit_tfrac_mavg_cl(:,:)    , & ! inout
                                ppfd_shaded_cl      = ppfd_shaded_tfrac_mavg_cl(:,:)    , & 
                                fleaf_sunlit_cl     = fleaf_sunlit_tfrac_mavg_cl(:,:)   , &
                                fn_rub_cl           = fn_rub_cl(:,:), &
                                fn_et_cl            = fn_et_cl(:,:), &
                                fn_pepc_cl          = fn_pepc_cl(:,:), &
                                fn_chl_cl           = fn_chl_cl(:,:), &
                                fn_oth_cl           = fn_oth_cl(:,:), &           ! inout
                                lai_cl              = lai_cl(:,:), &              ! out
                                cumm_lai_cl         = cumm_lai_cl(:,:), &
                                leaf_nitrogen_cl    = leaf_nitrogen_cl(:,:), &    ! out
                                inquire_n_fractions = .TRUE.)                     ! optional in
  
       !! canopy calculation fluxes
       CALL reset_veg_fluxes(tile,options)             !< mo_veg_util
       CALL update_reset_q_assimi_fluxes(tile, model%options)
       CALL update_radiation_par(tile,options)
       CALL update_canopy_fluxes(tile,options)  !< mo_q_assimi_process
  
       WRITE(lun,'(I1,3F16.7)') i, &
             net_assimilation(ipoi), canopy_cond(ipoi) *1000._wp, &
             co2_conc_leaf(ipoi)/co2_mixing_ratio(ipoi)
    END DO
  
  END SUBROUTINE run_test_canopy_model

END MODULE mo_qs_test_canopy_model
