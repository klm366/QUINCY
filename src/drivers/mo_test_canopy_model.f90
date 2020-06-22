!>
!!  @par Copyright
!!    QUINCY terrestrial biosphere model with coupled carbon, nitrogen, and phosphorus cycles
!!    Copyright (C) 2019 Max Planck Institute for Biogeochemistry (MPI-BGC), Jena, Germany 
!!    Contributors: Caldararu S., Engel J., Kern M., Pallandt M., Schnur R., Thum T., Yu L., Zaehle S.
!!
!!    This file is part of the QUINCY model.
!!    The QUINCY model is free software: you can redistribute it and/or modify
!!    it under the terms of the GNU General Public License version 3,
!!    as published by the Free Software Foundation.
!!
!!    The QUINCY model is distributed in the hope that it will be useful,
!!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!!    GNU General Public License for more details.
!!
!!    You should have received a copy of the GNU General Public License
!!    along with the QUINCY model.  If not, see <https://www.gnu.org/licenses/>.
!! 
!! @brief
!!  routine for running the test_canopy_model \n
!!  run QUINCY for one time step
!!
!! @author
!!  Sönke Zaehle, MPI-BGC Jena
!!
!! @par Revision History
!!  copy/pasted from mo_land_model_interface
!!
MODULE mo_test_canopy_model

#include "dsl_jsbach.inc"

  USE mo_kind

  IMPLICIT NONE


  PRIVATE
  PUBLIC :: run_test_canopy_model

  CHARACTER(len=*), PARAMETER :: modname = 'mo_test_canopy_model'

CONTAINS

  !-----------------------------------------------------------------------------------------------------
  !> Run the test_canopy model from the land model interface
  !!
  !! for one timestep
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE run_test_canopy_model(tile, options)

    USE mo_jsb_class,             ONLY: Get_model
    USE mo_jsb_tile_class,        ONLY: t_jsb_tile_abstract
    USE mo_jsb_task_class,        ONLY: t_jsb_task_options
    USE mo_jsb_model_class,       ONLY: t_jsb_model
    USE mo_jsb_lctlib_class,      ONLY: t_lctlib_element
    USE mo_jsb_process_class,     ONLY: A2L_, L2A_, SEB_, TURB_, SSE_, HYDRO_, VEG_, HD_, RAD_, ASSIMI_, PHENO_, SB_
    USE mo_jsb_grid_class,        ONLY: t_jsb_vgrid
    USE mo_jsb_grid,              ONLY: Get_vgrid
    USE mo_rad_constants,         ONLY: rad2ppfd
    USE mo_atmland_constants,     ONLY: frac_vis_swtot_srf, frac_par_swvis_srf,eps_vpd, standard_press_srf, def_co2_mixing_ratio
    USE mo_atmland_util,          ONLY: calc_spec_humidity_sat
    USE mo_veg_util,              ONLY: reset_veg_fluxes
    USE mo_jsb_phys_constants,    ONLY: molar_mass_N, Tzero, cvdifts, grav
    USE mo_veg_canopy             ! .. e.g. calc_canopy_layers()

    ! Use of process configurations (t_PROC_config)
    dsl4jsb_Use_config(ASSIMI_)    
    
    ! Use of process memories
    dsl4jsb_Use_memory(ASSIMI_)     !  USE mo_assimi_memory_class,   ONLY: t_assimi_memory
    dsl4jsb_Use_memory(RAD_)        !  USE mo_rad_memory_class,      ONLY: t_rad_memory
    dsl4jsb_Use_memory(VEG_)        !  USE mo_veg_memory_class,      ONLY: t_veg_memory
    dsl4jsb_Use_memory(SSE_)        !  USE mo_sse_memory_class,      ONLY: t_sse_memory

    ! Process interfaces:
    USE mo_atmland_interface                ! update_a2l
    USE mo_rad_interface                    ! update_radiation
    USE mo_assimi_interface                 ! update_canopy_fluxes, reset_assimi_fluxes
    ! get parameters defined in init_slm_test_canopy()

    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    CLASS(t_jsb_tile_abstract), INTENT(inout)     :: tile         !< one tile with data structure for one lct
    TYPE(t_jsb_task_options),   INTENT(in)        :: options      !< model options
    ! ---------------------------
    ! 0.2 Local
    INTEGER, PARAMETER                        :: npft  = 8         !< number of PFTs simulated by the test_canopy model
    TYPE(t_jsb_model),      POINTER           :: model
    TYPE(t_lctlib_element), POINTER           :: lctlib            !< land-cover-type library - parameter across pft's
    TYPE(t_jsb_vgrid),      POINTER           :: vgrid_canopy      ! Vertical grid
    INTEGER                                   :: ncanopy           ! number of canopy layers
    REAL(wp), DIMENSION(options%nc)           :: &
                  t_air           , &
                  press_srf       , &
                  q_air           , &
                  co2_mixing_ratio, &
                  rain            , &
                  snow           , &
                  swvis_srf_down , &
                  swnir_srf_down , &
                  swpar_srf_down , &
                  swpar_fdiffuse , &
                  wind_air       , &
                  wind_10m       , &
                  lw_srf_down    , &
                  drag_srf       , &
                  t_acoef        , &
                  t_bcoef        , &
                  q_acoef        , &
                  q_bcoef        , &
                  pch            , &
                  cos_zenith_angle, &
                  ga             , &
                  !latitude      , &
                  !longitude     , &
                  vpd
    REAL(wp),DIMENSION(options%nc)      :: hlp1
    INTEGER       :: i,lun
    INTEGER       :: ref_pft                          ! PFT identifier
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
    REAL(wp),ALLOCATABLE    :: canopy_layer_depth(:)  ! thickness of each single canopy layer [m2 m-2 LAI]
    INTEGER                 :: ipoi, ichunk           ! number of the chunk: ipoi 1:nc / ichunk ics:ice
    INTEGER                 :: iblk, ics, ice, nc
    REAL(wp)                :: lctlib_g1              !< set to g1_medlyn or g1_bberry depending on canopy_conductance_scheme
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':run_test_canopy_model'
    ! ---------------------------
    ! 0.3 Declare Memory
    ! Declare process configuration and memory Pointers
    dsl4jsb_Def_config(ASSIMI_)
    dsl4jsb_Def_memory(ASSIMI_)                       !   TYPE(t_assimi_memory), POINTER      :: assimi__mem
    dsl4jsb_Def_memory(RAD_)                          !   TYPE(t_rad_memory),    POINTER      :: rad__mem
    dsl4jsb_Def_memory(VEG_)                          !   TYPE(t_veg_memory),    POINTER      :: veg__mem
    dsl4jsb_Def_memory(SSE_)                         !   TYPE(t_sse_memory),   POINTER      :: sse__mem
    ! Declare pointers to variables in memory
    ! ASSIMI_
    dsl4jsb_Real2D_onChunk      :: net_assimilation 
    dsl4jsb_Real2D_onChunk      :: canopy_cond
    dsl4jsb_Real2D_onChunk      :: co2_conc_leaf
    dsl4jsb_Real2D_onChunk      :: soa_tsoa_mavg
    ! RAD_
    dsl4jsb_Real3D_onChunk      :: ppfd_sunlit_tfrac_mavg_cl
    dsl4jsb_Real3D_onChunk      :: ppfd_shaded_tfrac_mavg_cl
    ! SSE_
    dsl4jsb_Real2D_onChunk      :: w_soil_root_pot 
    ! VEG_
    dsl4jsb_Real2D_onChunk      :: lai 
    dsl4jsb_Real2D_onChunk      :: t_jmax_opt_mavg
    dsl4jsb_Real3D_onChunk      :: lai_cl
    dsl4jsb_Real3D_onChunk      :: cumm_lai_cl
    dsl4jsb_Real3D_onChunk      :: leaf_nitrogen_cl
    dsl4jsb_Real3D_onChunk      :: fn_rub_cl
    dsl4jsb_Real3D_onChunk      :: fn_et_cl
    dsl4jsb_Real3D_onChunk      :: fn_pepc_cl
    dsl4jsb_Real3D_onChunk      :: fn_chl_cl
    dsl4jsb_Real3D_onChunk      :: fn_oth_cl
    dsl4jsb_Real3D_onChunk      :: fleaf_sunlit_tfrac_mavg_cl
    ! Get local variables from options argument
    iblk    = options%iblk
    ics     = options%ics
    ice     = options%ice
    nc      = options%nc
    ! ---------------------------
    ! 0.4 Process Activity, Debug Option
!    IF (.NOT. tile%Is_process_active(xx_)) RETURN
!    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')
    ! ---------------------------
    ! 0.5 Get Memory
    ! Get process variables (Set pointers to variables in memory)
    model  => Get_model(tile%owner_model_id)
    ! if pftid is set in the namelist use this pftid
    IF(tile%lcts%lib_id>0)THEN
      ref_pft = tile%lcts%lib_id
    ! else use pftid 1
    ELSE
      ref_pft = 0
      tile%lcts%lib_id = 1
    ENDIF
    lctlib        => model%lctlib(tile%lcts%lib_id)
    vgrid_canopy  => Get_vgrid('canopy_layer')
    ncanopy       =  vgrid_canopy%n_levels
    ALLOCATE(canopy_layer_depth(ncanopy))
    canopy_layer_depth(:) = vgrid_canopy%dz(:)
    ! Get process config
    dsl4jsb_Get_config(ASSIMI_)
    ! Get process memories
    dsl4jsb_Get_memory(ASSIMI_)
    dsl4jsb_Get_memory(RAD_)
    dsl4jsb_Get_memory(VEG_)
    dsl4jsb_Get_memory(SSE_)
    ! ---------------------------
    dsl4jsb_Get_var2D_onChunk(ASSIMI_, net_assimilation)
    dsl4jsb_Get_var2D_onChunk(ASSIMI_, canopy_cond)
    dsl4jsb_Get_var2D_onChunk(ASSIMI_, co2_conc_leaf)
    dsl4jsb_Get_var2D_onChunk(ASSIMI_, soa_tsoa_mavg)
    ! ---------------------------
    dsl4jsb_Get_var3D_onChunk(RAD_, ppfd_sunlit_tfrac_mavg_cl)      ! in
    dsl4jsb_Get_var3D_onChunk(RAD_, ppfd_shaded_tfrac_mavg_cl)      ! in
    ! ---------------------------
    dsl4jsb_Get_var2D_onChunk(SSE_, w_soil_root_pot)
    ! ---------------------------
    dsl4jsb_Get_var2D_onChunk(VEG_, lai)
    dsl4jsb_Get_var2D_onChunk(VEG_, t_jmax_opt_mavg)
    dsl4jsb_Get_var3D_onChunk(VEG_, lai_cl)
    dsl4jsb_Get_var3D_onChunk(VEG_, cumm_lai_cl)
    dsl4jsb_Get_var3D_onChunk(VEG_, leaf_nitrogen_cl)
    dsl4jsb_Get_var3D_onChunk(VEG_, fn_rub_cl)
    dsl4jsb_Get_var3D_onChunk(VEG_, fn_et_cl)
    dsl4jsb_Get_var3D_onChunk(VEG_, fn_pepc_cl)
    dsl4jsb_Get_var3D_onChunk(VEG_, fn_chl_cl)
    dsl4jsb_Get_var3D_onChunk(VEG_, fn_oth_cl)
    dsl4jsb_Get_var3D_onChunk(VEG_, fleaf_sunlit_tfrac_mavg_cl)     ! in

    ! define chunk for output: ipoi 1:nc / ichunk ics:ice
    ipoi   = 1
    ichunk = ics+ipoi-1


    ! ------------------------------------------------------------------------------------------------------------
    ! Go science
    ! ------------------------------------------------------------------------------------------------------------


    !> 1.0 Init
    !!

    !> 1.1 set g1 according to canopy_conductance_scheme
    !! 
    SELECT CASE(TRIM(dsl4jsb_Config(ASSIMI_)%canopy_conductance_scheme))
      CASE ("medlyn")
        lctlib_g1 = lctlib%g1_medlyn
      CASE ("ballberry")
        lctlib_g1 = lctlib%g1_bberry
    END SELECT

    !> 1.2 Initalisation of met-forcing which is assumed constant
    !!
    !!
    t_air(:)            = ref_t_air + TZero
    press_srf(:)        = standard_press_srf
    q_air(:)            = calc_spec_humidity_sat(t_air(:),press_srf(:)) - ref_vpd / press_srf(:) * eps_vpd
    co2_mixing_ratio(:) = def_co2_mixing_ratio   ! def_co2_mixing_ratio is used instead of ref_co2_mixing_ratio
    rain(:)             = 0.0_wp
    snow(:)             = 0.0_wp
    swpar_srf_down(:)   = ref_swdown
    swvis_srf_down(:)   = swpar_srf_down(:) / frac_par_swvis_srf
    swnir_srf_down(:)   = swvis_srf_down(:) * (1._wp - frac_vis_swtot_srf) / frac_vis_swtot_srf
    swpar_fdiffuse(:)   = ref_fdiffuse
    wind_air(:)         = ref_wind
    wind_10m(:)         = wind_air(:)
    lw_srf_down(:)      = 0.0_wp
    ga(:)               = ref_ga 
    drag_srf(:)         = ga(:) / wind_air(:) * cvdifts * grav * model%options%dtime 
    t_acoef(:)          = 0.0_wp
    t_bcoef(:)          = 0.0_wp
    q_acoef(:)          = 0.0_wp
    q_bcoef(:)          = 0.0_wp
    pch(:)              = 0.0_wp
    cos_zenith_angle(:) = ref_cos_angle
    soa_tsoa_mavg(:)    = 7.0_wp        ! is initialized in assimi init with -10.0_wp as the model starts in winter
    !latitude                  = 50.0_wp
    !longitude                 = 0.0_wp


    CALL update_a2l( &
           tile                                            , &
           options                                         , &
           t_air                   = t_air(:)              , &
           q_air                   = q_air(:)              , &
           wind_air                = wind_air(:)           , &
           wind_10m                = wind_10m(:)           , &
           swvis_srf_down          = swvis_srf_down(:)     , &
           swnir_srf_down          = swnir_srf_down(:)     , &
           swpar_srf_down          = swpar_srf_down(:)     , &
           swpar_fdiffuse          = swpar_fdiffuse(:)     , &
           press_srf               = press_srf(:)          , &
           drag_srf                = drag_srf(:)           , &
           cos_zenith_angle        = cos_zenith_angle(:)   , &
           co2_mixing_ratio        = co2_mixing_ratio(:)   , &
           ga                      = ga(:)                 )

    !> 1.2 Initialisation of soil and vegetation (so far relevant)
    !!
    !!
    lai(:)   = ref_lai
    IF(ref_pft==0) THEN
      veg__mem%pool%leaf%nitrogen%ptr(:,iblk)     = ref_nleaf / molar_mass_N * lai(:)
    ELSE
      veg__mem%pool%leaf%nitrogen%ptr(:,iblk)     = lai(:) / lctlib%sla / lctlib%cn_leaf
      veg__mem%pool%leaf%phosphorus%ptr(:,iblk)   = veg__mem%pool%leaf%nitrogen%ptr(:,iblk) / lctlib%np_leaf
    ENDIF
    CALL calc_canopy_layers(nc                                             , & ! in
                            ncanopy                                        , &
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
                            dsl4jsb_Config(ASSIMI_)%flag_optimal_Nfraction , &
                            TRIM(dsl4jsb_Config(ASSIMI_)%canopy_conductance_scheme)       , & ! ASSIMI_ config
                            canopy_layer_depth(:)                                         , &
                            leaf_nitrogen       = veg__mem%pool%leaf%nitrogen%ptr(:,iblk) , &
                            lai                 = lai(:)                          , &
                            ppfd_sunlit_cl      = ppfd_sunlit_tfrac_mavg_cl(:,:)  , &
                            ppfd_shaded_cl      = ppfd_shaded_tfrac_mavg_cl(:,:)  , & 
                            fleaf_sunlit_cl     = fleaf_sunlit_tfrac_mavg_cl(:,:) , &
                            lai_cl              = lai_cl(:,:), &
                            cumm_lai_cl         = cumm_lai_cl(:,:), &
                            leaf_nitrogen_cl    = leaf_nitrogen_cl(:,:), &
                            fn_rub_cl           = fn_rub_cl(:,:), &
                            fn_et_cl            = fn_et_cl(:,:), &
                            fn_pepc_cl          = fn_pepc_cl(:,:), &
                            fn_chl_cl           = fn_chl_cl(:,:), &
                            fn_oth_cl           = fn_oth_cl(:,:), &
                            inquire_n_fractions = .TRUE.)
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
     CALL update_a2l(tile, &
           options            , &
           swvis_srf_down = swvis_srf_down(:), &
           swnir_srf_down = swnir_srf_down(:), &
           swpar_srf_down = swpar_srf_down(:))

     !! canopy calculation fluxes
     CALL reset_veg_fluxes(tile,options)             !< mo_veg_util
     CALL reset_assimi_fluxes(tile, model%options)   !< mo_assimi_util
     CALL update_radiation(tile,options)      !< mo_rad_update_radiation
     CALL update_canopy_fluxes(tile,options)  !< mo_assimi_update_canopy_fluxes

     WRITE(lun,'(4F16.7)') swpar_srf_down(ipoi) * rad2ppfd, &
           net_assimilation(ipoi), canopy_cond(ipoi) *1000._wp, &
           co2_conc_leaf(ipoi)/co2_mixing_ratio(ipoi) !,&
!           net_assimilation(ipoi),
!           co2_conc_leaf(ipoi)/co2_mixing_ratio(ipoi),veg__mem%pool%leaf%nitrogen%ptr(ipoi,iblk)*molar_mass_N,&
!           net_assimilation(ipoi)/(veg__mem%pool%leaf%nitrogen%ptr(ipoi,iblk)*molar_mass_N)

  END DO
  swpar_srf_down(:) = ref_swdown
  swvis_srf_down(:) = swpar_srf_down(:) / frac_par_swvis_srf
  swnir_srf_down(:) = swvis_srf_down(:) * (1._wp - frac_vis_swtot_srf) / frac_vis_swtot_srf
  CALL update_a2l(tile, &
           options            , &
           swvis_srf_down = swvis_srf_down, &
           swnir_srf_down = swnir_srf_down, &
           swpar_srf_down = swpar_srf_down)

  !> 2.2 Temperature (-5,60 °C)
  !!
  !!
  lun = lun + 1
  WRITE(lun,*) 'T A GC CI'
  t_air(:) = - 5.0_wp + Tzero
  hlp1(:)  = t_jmax_opt_mavg(:)

  DO i = 1,70

     t_air(:) = t_air(:) + 1.0_wp
     q_air(:) = 0.8_wp * calc_spec_humidity_sat(t_air(:),press_srf(:))
     t_jmax_opt_mavg(:) = 20._wp+0.35_wp*(t_air(:) - Tzero)
     CALL update_a2l(tile, &
           options            , &
           t_air = t_air, &
           q_air = q_air)

     !! canopy calculation fluxes
     CALL reset_veg_fluxes(tile,options)             !< mo_veg_util
     CALL reset_assimi_fluxes(tile, model%options)   !< mo_assimi_util
     CALL update_radiation(tile,options)      !< mo_rad_update_radiation
     CALL update_canopy_fluxes(tile,options)  !< mo_assimi_update_canopy_fluxes

     WRITE(lun,'(5F16.7)') t_air(ipoi), &
           net_assimilation(ipoi), canopy_cond(ipoi) *1000._wp, &
           co2_conc_leaf(ipoi)/co2_mixing_ratio(ipoi)

  END DO
  t_air(:)    = ref_t_air + Tzero
  q_air(:)    = calc_spec_humidity_sat(t_air(:),press_srf(:)) - ref_vpd / press_srf(:) * eps_vpd
  CALL update_a2l(tile, &
        options            , &
        t_air = t_air, &
        q_air = q_air)
  t_jmax_opt_mavg(:) = hlp1(:)

  !> 2.3 atmospheric CO2 (175 - 1000 ppm)
  !!
  !!
  lun = lun + 1
  WRITE(lun,*) 'CO2 A GC CI'
  co2_mixing_ratio(:) = 150._wp
  DO i = 1,34

     co2_mixing_ratio(:) = co2_mixing_ratio(:) + 25._wp
     CALL update_a2l(tile, &
           options            , &
           co2_mixing_ratio = co2_mixing_ratio)

     !! canopy calculation fluxes
     CALL reset_veg_fluxes(tile,options)             !< mo_veg_util
     CALL reset_assimi_fluxes(tile, model%options)   !< mo_assimi_util
     CALL update_radiation(tile,options)      !< mo_rad_update_radiation
     CALL update_canopy_fluxes(tile,options)  !< mo_assimi_update_canopy_fluxes

     WRITE(lun,'(4F16.7)') co2_mixing_ratio(ipoi), &
           net_assimilation(ipoi), canopy_cond(ipoi) *1000._wp, &
           co2_conc_leaf(ipoi)/co2_mixing_ratio(ipoi)

  END DO
  co2_mixing_ratio(:) = ref_co2_mixing_ratio
  CALL update_a2l(tile    , &
                       options , &
                       co2_mixing_ratio = co2_mixing_ratio)

  !> 2.4 VPD (0-3.3 kPa)
  !!
  !!
  lun = lun + 1
  WRITE(lun,*) 'VPD A GC CI'
  DO i = 1,30
     q_air(:) = calc_spec_humidity_sat(t_air(:),press_srf(:)) - i * 100.0_wp / press_srf(:) * eps_vpd
     CALL update_a2l(tile    , &
                          options , &
                          q_air = q_air)

     !! canopy calculation fluxes
     CALL reset_veg_fluxes(tile,options)             !< mo_veg_util
     CALL reset_assimi_fluxes(tile, model%options)   !< mo_assimi_util
     CALL update_radiation(tile,options)      !< mo_rad_update_radiation
     CALL update_canopy_fluxes(tile,options)  !< mo_assimi_update_canopy_fluxes

     vpd(:) = ( calc_spec_humidity_sat(t_air(:),press_srf(:)) - q_air(:) ) / eps_vpd * press_srf(:) / 1000.0_wp

     WRITE(lun,'(4F16.7)') vpd(ipoi), &
           net_assimilation(ipoi), canopy_cond(ipoi) *1000._wp, &
           co2_conc_leaf(ipoi)/co2_mixing_ratio(ipoi)

  END DO
  q_air(:) = calc_spec_humidity_sat(t_air(:),press_srf(:)) - ref_vpd / press_srf(:) * eps_vpd
  vpd(:) = ( calc_spec_humidity_sat(t_air(:),press_srf(:)) - q_air(:) ) / eps_vpd * press_srf(:) / 1000.0_wp
     CALL update_a2l(tile    , &
                          options , &
                          q_air = q_air)

  !> 2.5 soil water potential (0 to -4)
  !!
  !!
  lun = lun + 1
  WRITE(lun,*) 'SW_pot A GC CI'
  w_soil_root_pot(:) = ref_soil_root_pot
  DO i = 1,80

     w_soil_root_pot(:) = w_soil_root_pot(:) - 0.05_wp

     !! canopy calculation fluxes
     CALL reset_veg_fluxes(tile,options)             !< mo_veg_util
     CALL reset_assimi_fluxes(tile, model%options)   !< mo_assimi_util
     CALL update_radiation(tile,options)      !< mo_rad_update_radiation
     CALL update_canopy_fluxes(tile,options)  !< mo_assimi_update_canopy_fluxes

     WRITE(lun,'(4F16.7)') w_soil_root_pot(ipoi) *100._wp, &
           net_assimilation(ipoi), canopy_cond(ipoi) *1000._wp, &
           co2_conc_leaf(ipoi)/co2_mixing_ratio(ipoi)

  END DO
  w_soil_root_pot(:) = ref_soil_root_pot

  !> 2.6 LAI (increase from 0 to 11)
  !!
  !!
  lun = lun + 1
  WRITE(lun,*) 'LAI A GC CI'
  lai(:) = 0.0_wp
  DO i = 1,22

     lai(:) = lai(:) + 0.5_wp
     veg__mem%pool%leaf%nitrogen%ptr(:,iblk)=ref_nleaf/molar_mass_N*lai(:)
     CALL calc_canopy_layers( nc                                             , &   ! in
                              ncanopy                                        , &
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
                              dsl4jsb_Config(ASSIMI_)%flag_optimal_Nfraction , &
                              TRIM(dsl4jsb_Config(ASSIMI_)%canopy_conductance_scheme)       , & ! ASSIMI_ config
                              canopy_layer_depth(:)                                         , &
                              leaf_nitrogen       = veg__mem%pool%leaf%nitrogen%ptr(:,iblk) , &
                              lai                 = lai(:)                          , &
                              ppfd_sunlit_cl      = ppfd_sunlit_tfrac_mavg_cl(:,:)  , &
                              ppfd_shaded_cl      = ppfd_shaded_tfrac_mavg_cl(:,:)  , & 
                              fleaf_sunlit_cl     = fleaf_sunlit_tfrac_mavg_cl(:,:) , &
                              lai_cl              = lai_cl(:,:), &
                              cumm_lai_cl         = cumm_lai_cl(:,:), &
                              leaf_nitrogen_cl    = leaf_nitrogen_cl(:,:), &
                              fn_rub_cl           = fn_rub_cl(:,:), &
                              fn_et_cl            = fn_et_cl(:,:), &
                              fn_pepc_cl          = fn_pepc_cl(:,:), &
                              fn_chl_cl           = fn_chl_cl(:,:), &
                              fn_oth_cl           = fn_oth_cl(:,:), &
                              inquire_n_fractions = .TRUE.)

     !! canopy calculation fluxes
     CALL reset_veg_fluxes(tile,options)             !< mo_veg_util
     CALL reset_assimi_fluxes(tile, model%options)   !< mo_assimi_util
     CALL update_radiation(tile,options)      !< mo_rad_update_radiation
     CALL update_canopy_fluxes(tile,options)  !< mo_assimi_update_canopy_fluxes

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

     veg__mem%pool%leaf%nitrogen%ptr(:,iblk)=0.125_wp*(i)/molar_mass_N*lai(:)
     CALL calc_canopy_layers( nc                                             , &   ! in
                              ncanopy                                        , &
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
                              dsl4jsb_Config(ASSIMI_)%flag_optimal_Nfraction , &
                              TRIM(dsl4jsb_Config(ASSIMI_)%canopy_conductance_scheme)       , & ! ASSIMI_ config
                              canopy_layer_depth(:)                                         , &
                              leaf_nitrogen       = veg__mem%pool%leaf%nitrogen%ptr(:,iblk) , &
                              lai                 = lai(:)                          , &
                              ppfd_sunlit_cl      = ppfd_sunlit_tfrac_mavg_cl(:,:)  , &
                              ppfd_shaded_cl      = ppfd_shaded_tfrac_mavg_cl(:,:)  , & 
                              fleaf_sunlit_cl     = fleaf_sunlit_tfrac_mavg_cl(:,:) , &
                              lai_cl              = lai_cl(:,:), &
                              cumm_lai_cl         = cumm_lai_cl(:,:), &
                              leaf_nitrogen_cl    = leaf_nitrogen_cl(:,:), &
                              fn_rub_cl           = fn_rub_cl(:,:), &
                              fn_et_cl            = fn_et_cl(:,:), &
                              fn_pepc_cl          = fn_pepc_cl(:,:), &
                              fn_chl_cl           = fn_chl_cl(:,:), &
                              fn_oth_cl           = fn_oth_cl(:,:), &
                              inquire_n_fractions = .TRUE.)

     !! canopy calculation fluxes
     CALL reset_veg_fluxes(tile,options)             !< mo_veg_util
     CALL reset_assimi_fluxes(tile, model%options)   !< mo_assimi_util
     CALL update_radiation(tile,options)      !< mo_rad_update_radiation
     CALL update_canopy_fluxes(tile,options)  !< mo_assimi_update_canopy_fluxes

     WRITE(lun,'(4F16.7)') veg__mem%pool%leaf%nitrogen%ptr(ichunk,iblk)/lai(ipoi)*molar_mass_N, &
           net_assimilation(ipoi), canopy_cond(ipoi) *1000._wp, &
           co2_conc_leaf(ipoi)/co2_mixing_ratio(ipoi)

  END DO
  veg__mem%pool%leaf%nitrogen%ptr(:,iblk)=ref_nleaf/molar_mass_N*lai(:)

  !> 2.8 PFT
  !!
  !!
  lun = lun + 1
  WRITE(lun,*) 'PFT A GC CI'
  DO i = 1,npft
     veg__mem%pool%leaf%nitrogen%ptr(:,iblk)  =lai(:) /model%lctlib(i)%sla /model%lctlib(i)%cn_leaf
     veg__mem%pool%leaf%phosphorus%ptr(:,iblk)=veg__mem%pool%leaf%nitrogen%ptr(:,iblk) /model%lctlib(i)%np_leaf
     ! change lctlib id
     tile%lcts%lib_id = i
     ! set pointer to the new PFT
     lctlib  => model%lctlib(tile%lcts%lib_id)
     CALL calc_canopy_layers( nc                                             , &   ! in
                              ncanopy                                        , &
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
                              dsl4jsb_Config(ASSIMI_)%flag_optimal_Nfraction , &
                              TRIM(dsl4jsb_Config(ASSIMI_)%canopy_conductance_scheme)       , & ! ASSIMI_ config
                              canopy_layer_depth(:)                                         , &
                              leaf_nitrogen       = veg__mem%pool%leaf%nitrogen%ptr(:,iblk) , &
                              lai                 = lai(:)                          , &
                              ppfd_sunlit_cl      = ppfd_sunlit_tfrac_mavg_cl(:,:)  , &
                              ppfd_shaded_cl      = ppfd_shaded_tfrac_mavg_cl(:,:)  , & 
                              fleaf_sunlit_cl     = fleaf_sunlit_tfrac_mavg_cl(:,:) , &
                              lai_cl              = lai_cl(:,:), &
                              cumm_lai_cl         = cumm_lai_cl(:,:), &
                              leaf_nitrogen_cl    = leaf_nitrogen_cl(:,:), &
                              fn_rub_cl           = fn_rub_cl(:,:), &
                              fn_et_cl            = fn_et_cl(:,:), &
                              fn_pepc_cl          = fn_pepc_cl(:,:), &
                              fn_chl_cl           = fn_chl_cl(:,:), &
                              fn_oth_cl           = fn_oth_cl(:,:), &
                              inquire_n_fractions = .TRUE.)

     !! canopy calculation fluxes
     CALL reset_veg_fluxes(tile,options)             !< mo_veg_util
     CALL reset_assimi_fluxes(tile, model%options)   !< mo_assimi_util
     CALL update_radiation(tile,options)      !< mo_rad_update_radiation
     CALL update_canopy_fluxes(tile,options)  !< mo_assimi_update_canopy_fluxes

     WRITE(lun,'(I1,3F16.7)') i, &
           net_assimilation(ipoi), canopy_cond(ipoi) *1000._wp, &
           co2_conc_leaf(ipoi)/co2_mixing_ratio(ipoi)

  END DO
  !STOP
  END SUBROUTINE run_test_canopy_model


END MODULE mo_test_canopy_model
