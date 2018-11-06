!>
!! @brief Contains the memory class for the radiation process.
!!
!! @author
!!  Jan Engel, MPI-BGC Jena (based on the JSBACH4 version from Rainer Schneck, MPI-M Hamburg)
!!
!! @par Revision History
!! First version by Jan Engel (20160819)
!!
MODULE mo_rad_memory_class
#ifndef __NO_JSBACH__

  USE mo_kind,                   ONLY: wp
  USE mo_exception,              ONLY: message
  USE mo_jsb_class,              ONLY: Get_model
  USE mo_jsb_memory_class,       ONLY: t_jsb_memory
  USE mo_jsb_var_class,          ONLY: t_jsb_var_real2d, t_jsb_var_real3d


  IMPLICIT NONE
  PRIVATE
  PUBLIC :: t_rad_memory, max_no_of_vars

  INTEGER, PARAMETER :: max_no_of_vars = 50

  !-----------------------------------------------------------------------------------------------------
  !> Type definition for rad memory
  !! 
  !-----------------------------------------------------------------------------------------------------
  TYPE, EXTENDS(t_jsb_memory) :: t_rad_memory

    TYPE(t_jsb_var_real2d)  :: &
      & sw_srf_net,            & !< net shortwave radiation at surface [W m-2]
      & swvis_srf_net,         & !< net shortwave radiation at surface in visible range [W m-2]
      & swnir_srf_net,         & !< net shortwave radiation at surface in NIR range [W m-2]
      & lw_srf_net,            & !< net longwave radiation at surface [W m-2]
      & rad_srf_net              !< net radiation at surface [W m-2]

    TYPE(t_jsb_var_real2d)  :: & 
      & alb_vis,               & !< Albedo of the whole tile in the VIS (visible) range [unitless]
      & alb_nir,               & !< Albedo of the whole tile in the NIR (near-infrared) range [unitless]
      & alb_vis_snow,          & !< Albedo of snow in the VIS range [unitless]
      & alb_nir_snow             !< Albedo of snow in the NIR range [unitless]

    ! Additional variables for soil
    TYPE(t_jsb_var_real2d)  :: &
      & alb_vis_soil,          & !< Albedo in VIS range of soil [unitless]
      & alb_nir_soil,          & !< Albedo in VIS range of soil [unitless]
      & arad_vis_soil,         & !< radiation in VIS range absorbed by the soil [W m-2]
      & arad_nir_soil            !< radiation in NIR range absorbed by the soil [W m-2]

    ! Additional variables for vegetation
    TYPE(t_jsb_var_real2d)  ::      &
      & alb_vis_can,                & !< Albedo of canopy (without snow cover) in VIS range [unitless]
      & alb_nir_can,                & !< Albedo of canopy (without snow cover) in NIR range [unitless]
      & arad_vis_can,               & !< radiation in VIS range absorbed by the canopy [W m-2]
      & arad_nir_can,               & !< radiation in NIR range absorbed by the canopy [W m-2]
      & appfd,                      & !< absorbed Photosynthetically Active Photon Flux Density [µmol m-2 s-1]
      & rfr_ratio_boc,              & !< red to far-red ratio at the bottom of the canopy [unitless]
      & rfr_ratio_boc_tvegdyn_mavg    !< red to far-red ratio at the bottom of the canopy [unitless]

    TYPE(t_jsb_var_real3d)  :: &   ! DIMENSION(ncanopy)
                            ppfd_sunlit_cl              , & !< Photosynthetically Active Photon Flux Density on sunlit leaves [µmol m-2 s-1]
                            ppfd_shaded_cl              , & !< Photosynthetically Active Photon Flux Density on shaded leaves [µmol m-2 s-1]
                            ppfd_sunlit_daytime_cl      , & !< average PPFD on sunlit leaves of the previous day [µmol m-2 s-1]
                            ppfd_shaded_daytime_cl      , & !< average PPFD on shaded leaves of the previous day [µmol m-2 s-1]
                            ppfd_sunlit_daytime_dacc_cl , & !< daily accumulated PPFD on sunlit leaves [µmol m-2 d-1] 
                            ppfd_shaded_daytime_dacc_cl , & !< daily accumulated PPFD on shaded leaves of the previous day [µmol m-2 d-1]
                            ppfd_sunlit_tfrac_mavg_cl   , & !< average PPFD on sunlit leaves at tfrac-timescale [µmol m-2 s-1] 
                            ppfd_shaded_tfrac_mavg_cl   , & !< average PPFD on shaded leaves at tfrac-timescale [µmol m-2 s-1]
                            ppfd_sunlit_tcnl_mavg_cl    , & !< average PPFD on sunlit leaves at tcnl-timescale [µmol m-2 s-1]
                            ppfd_shaded_tcnl_mavg_cl        !< average PPFD on shaded leaves at tcnl-timescale [µmol m-2 s-1]

    TYPE(t_jsb_var_real3d)  :: &   ! DIMENSION(ncanopy)
                            arad_sunlit_vis_cl   , & !< absorbed radiation in visible band on sunlit leaves [W m-2]
                            arad_shaded_vis_cl   , & !< absorbed radiation in visible band on shaded leaves [W m-2]
                            fleaf_sunlit_vis_cl  , & !< fraction of leaves that are sunlit in VIS band [unitless]
                            arad_sunlit_nir_cl   , & !< absorbed radiation in NIR band on sunlit leaves [W m-2]
                            arad_shaded_nir_cl   , & !< absorbed radiation in NIR band on shaded leaves [W m-2]
                            fleaf_sunlit_nir_cl      !< fraction of leaves that are sunlit in NIR band [unitless]

  CONTAINS
    PROCEDURE :: Init => Init_rad_memory
  END TYPE t_rad_memory


  CHARACTER(len=*), PARAMETER :: modname = 'mo_rad_memory_class'


CONTAINS

  !-----------------------------------------------------------------------------------------------------
  !> initialize RAD_ process memory variables
  !!
  !! 
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE Init_rad_memory(mem, prefix, suffix, lct_ids, model_id)

!    USE mo_jsb_varlist,       ONLY: t_jsb_varlist, jsb_add_var
    USE mo_slm_varlist,       ONLY: jsb_add_var   ! replace by mo_jsb_varlist later    *slmtemp
    USE mo_jsb_io,            ONLY: grib_bits, t_cf, t_grib1, t_grib2   ! , TSTEP_CONSTANT, tables
    USE mo_jsb_grid_class,    ONLY: t_jsb_grid, t_jsb_vgrid
    USE mo_jsb_grid,          ONLY: Get_grid, Get_vgrid
    USE mo_jsb_model_class,   ONLY: t_jsb_model

    IMPLICIT NONE
    ! ---------------------------
    ! 0.1.InOut
    CLASS(t_rad_memory), INTENT(inout), TARGET  :: mem             !> rad memory
    CHARACTER(len=*),     INTENT(in)            :: prefix          !> process name
    CHARACTER(len=*),     INTENT(in)            :: suffix          !> tile name
    INTEGER,              INTENT(in)            :: lct_ids !(:)    !> in jsbach4: vector of land cover type IDs
    INTEGER,              INTENT(in)            :: model_id        !> model ID model\%id
    ! ---------------------------
    ! 0.2 Local
    TYPE(t_jsb_model), POINTER :: model
    TYPE(t_jsb_grid),  POINTER :: hgrid                        ! Horizontal grid
    TYPE(t_jsb_vgrid), POINTER :: vgrid_surface                ! Vertical grid
    TYPE(t_jsb_vgrid), POINTER :: vgrid_canopy                 ! Vertical grid
    INTEGER                    :: table
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':Init_rad_memory'
    ! ---------------------------
    ! 0.5 Get Memory
    model         => Get_model(model_id)
    table         = 1     ! in jsbach4 this line is: table = tables(1)     *slmtemp
!    hgrid         => Get_grid(mem%grid_id)
    vgrid_surface => Get_vgrid('surface')
    vgrid_canopy  => Get_vgrid('canopy_layer')



    !> 1.0  init memory
    !! 
    IF (.NOT. ASSOCIATED(mem%sw_srf_net%ptr)) THEN

      CALL jsb_add_var(mem, 'sw_srf_net', mem%sw_srf_net, &
            hgrid, vgrid_surface, &
            t_cf('sw_srf_net', 'W m-2', 'net shortwave radiation at surface'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'swvis_srf_net', mem%swvis_srf_net, &
            hgrid, vgrid_surface, &
            t_cf('swvis_srf_net', 'W m-2', 'net shortwave radiation at surface in visible range'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'swnir_srf_net', mem%swnir_srf_net, &
            hgrid, vgrid_surface, &
            t_cf('swnir_srf_net', 'W m-2', 'net shortwave radiation at surface in NIR range'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'lw_srf_net', mem%lw_srf_net, &
            hgrid, vgrid_surface, &
            t_cf('lw_srf_net', 'W m-2', 'net longwave radiation at surface'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'rad_srf_net', mem%rad_srf_net, &
            hgrid, vgrid_surface, &
            t_cf('rad_srf_net', 'W m-2', 'net radiation at surface'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'alb_vis', mem%alb_vis, &
            hgrid, vgrid_surface, &
            t_cf('alb_vis', 'unitless', 'Albedo of the whole tile in the VIS (visible) range'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'alb_nir', mem%alb_nir, &
            hgrid, vgrid_surface, &
            t_cf('alb_nir', 'unitless', 'Albedo of the whole tile in the NIR (near-infrared) range'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'alb_vis_snow', mem%alb_vis_snow, &
            hgrid, vgrid_surface, &
            t_cf('alb_vis_snow', 'unitless', 'Albedo of snow in the VIS range'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'alb_nir_snow', mem%alb_nir_snow, &
            hgrid, vgrid_surface, &
            t_cf('alb_nir_snow', 'unitless', 'Albedo of snow in the NIR range'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'alb_vis_soil', mem%alb_vis_soil, &
            hgrid, vgrid_surface, &
            t_cf('alb_vis_soil', 'unitless', 'Albedo in VIS range of soil'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'alb_nir_soil', mem%alb_nir_soil, &
            hgrid, vgrid_surface, &
            t_cf('alb_nir_soil', 'unitless', 'Albedo in VIS range of soil'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'arad_vis_soil', mem%arad_vis_soil, &
            hgrid, vgrid_surface, &
            t_cf('arad_vis_soil', 'W m-2', 'radiation in VIS range absorbed by the soil'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'arad_nir_soil', mem%arad_nir_soil, &
            hgrid, vgrid_surface, &
            t_cf('arad_nir_soil', 'W m-2', 'radiation in NIR range absorbed by the soil'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'alb_vis_can', mem%alb_vis_can, &
            hgrid, vgrid_surface, &
            t_cf('alb_vis_can', 'unitless', 'Albedo of canopy (without snow cover) in VIS range'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'alb_nir_can', mem%alb_nir_can, &
            hgrid, vgrid_surface, &
            t_cf('alb_nir_can', 'unitless', 'Albedo of canopy (without snow cover) in NIR range'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'arad_vis_can', mem%arad_vis_can, &
            hgrid, vgrid_surface, &
            t_cf('arad_vis_can', 'W m-2', 'radiation in VIS range absorbed by the canopy'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'arad_nir_can', mem%arad_nir_can, &
            hgrid, vgrid_surface, &
            t_cf('arad_nir_can', 'W m-2', 'radiation in NIR range absorbed by the canopy'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'appfd', mem%appfd, &
            hgrid, vgrid_surface, &
            t_cf('appfd', 'µmol m-2 s-1', 'absorbed Photosynthetically Active Photon Flux Density'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'rfr_ratio_boc', mem%rfr_ratio_boc, &
            hgrid, vgrid_surface, &
            t_cf('rfr_ratio_boc', 'unitless', 'red to far-red ratio at the bottom of the canopy'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'rfr_ratio_boc_tvegdyn_mavg', mem%rfr_ratio_boc_tvegdyn_mavg, &
            hgrid, vgrid_surface, &
            t_cf('rfr_ratio_boc_tvegdyn_mavg', 'unitless', 'red to far-red ratio at the bottom of the canopy'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'ppfd_sunlit_cl', mem%ppfd_sunlit_cl, &
            hgrid, vgrid_canopy, &
            t_cf('ppfd_sunlit_cl', 'µmol m-2 s-1', 'Photosynthetically Active Photon Flux Density on sunlit leaves'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'ppfd_shaded_cl', mem%ppfd_shaded_cl, &
            hgrid, vgrid_canopy, &
            t_cf('ppfd_shaded_cl', 'µmol m-2 s-1', 'Photosynthetically Active Photon Flux Density on shaded leaves'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'ppfd_sunlit_daytime_cl', mem%ppfd_sunlit_daytime_cl, &
            hgrid, vgrid_canopy, &
            t_cf('ppfd_sunlit_daytime_cl', 'µmol m-2 s-1', 'average PPFD on sunlit leaves of the previous day'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'ppfd_shaded_daytime_cl', mem%ppfd_shaded_daytime_cl, &
            hgrid, vgrid_canopy, &
            t_cf('ppfd_shaded_daytime_cl', 'µmol m-2 s-1', 'average PPFD on shaded leaves of the previous day'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'ppfd_sunlit_daytime_dacc_cl', mem%ppfd_sunlit_daytime_dacc_cl, &
            hgrid, vgrid_canopy, &
            t_cf('ppfd_sunlit_daytime_dacc_cl', 'µmol m-2 d-1', 'daily accumulated PPFD on sunlit leaves'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'ppfd_shaded_daytime_dacc_cl', mem%ppfd_shaded_daytime_dacc_cl, &
            hgrid, vgrid_canopy, &
            t_cf('ppfd_shaded_daytime_dacc_cl', 'µmol m-2 d-1', 'daily accumulated PPFD on shaded leaves of the previous day'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'ppfd_sunlit_tfrac_mavg_cl', mem%ppfd_sunlit_tfrac_mavg_cl, &
            hgrid, vgrid_canopy, &
            t_cf('ppfd_sunlit_tfrac_mavg_cl', 'µmol m-2 s-1', 'average PPFD on sunlit leaves at tfrac-timescale'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'ppfd_shaded_tfrac_mavg_cl', mem%ppfd_shaded_tfrac_mavg_cl, &
            hgrid, vgrid_canopy, &
            t_cf('ppfd_shaded_tfrac_mavg_cl', 'µmol m-2 s-1', 'average PPFD on shaded leaves at tfrac-timescale'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'ppfd_sunlit_tcnl_mavg_cl', mem%ppfd_sunlit_tcnl_mavg_cl, &
            hgrid, vgrid_canopy, &
            t_cf('ppfd_sunlit_tcnl_mavg_cl', 'µmol m-2 s-1', 'average PPFD on sunlit leaves at tcnl-timescale'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'ppfd_shaded_tcnl_mavg_cl', mem%ppfd_shaded_tcnl_mavg_cl, &
            hgrid, vgrid_canopy, &
            t_cf('ppfd_shaded_tcnl_mavg_cl', 'µmol m-2 s-1', 'average PPFD on shaded leaves at tcnl-timescale'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue

      CALL jsb_add_var(mem, 'arad_sunlit_vis_cl', mem%arad_sunlit_vis_cl, &
            hgrid, vgrid_canopy, &
            t_cf('arad_sunlit_vis_cl', 'W m-2', 'absorbed radiation in visible band on sunlit leaves'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue

      CALL jsb_add_var(mem, 'arad_shaded_vis_cl', mem%arad_shaded_vis_cl, &
            hgrid, vgrid_canopy, &
            t_cf('arad_shaded_vis_cl', 'W m-2', 'absorbed radiation in visible band on shaded leaves'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue

      CALL jsb_add_var(mem, 'fleaf_sunlit_vis_cl', mem%fleaf_sunlit_vis_cl, &
            hgrid, vgrid_canopy, &
            t_cf('fleaf_sunlit_vis_cl', 'unitless', 'fraction of leaves that are sunlit in VIS band'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue

      CALL jsb_add_var(mem, 'arad_sunlit_nir_cl', mem%arad_sunlit_nir_cl, &
            hgrid, vgrid_canopy, &
            t_cf('arad_sunlit_nir_cl', 'W m-2', 'absorbed radiation in NIR band on sunlit leaves'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue

      CALL jsb_add_var(mem, 'arad_shaded_nir_cl', mem%arad_shaded_nir_cl, &
            hgrid, vgrid_canopy, &
            t_cf('arad_shaded_nir_cl', 'W m-2', 'absorbed radiation in NIR band on shaded leaves'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue

      CALL jsb_add_var(mem, 'fleaf_sunlit_nir_cl', mem%fleaf_sunlit_nir_cl, &
            hgrid, vgrid_canopy, &
            t_cf('fleaf_sunlit_nir_cl', 'unitless', 'fraction of leaves that are sunlit in NIR band'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue


    END IF

  END SUBROUTINE Init_rad_memory

#endif
END MODULE mo_rad_memory_class
