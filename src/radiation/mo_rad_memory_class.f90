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

  USE mo_kind,                   ONLY: wp
  USE mo_exception,              ONLY: message
  USE mo_jsb_class,              ONLY: Get_model

  USE mo_jsb_memory_class,       ONLY: t_jsb_memory
  USE mo_jsb_var_class,          ONLY: t_jsb_var_real2d, t_jsb_var_real3d

  USE mo_slm_varlist,            ONLY: slm_add_var   ! replace by mo_jsb_varlist later   *slmtemp



  IMPLICIT NONE
  PRIVATE
  PUBLIC :: t_rad_memory, max_no_of_vars

  INTEGER, PARAMETER :: max_no_of_vars = 40   ! default JSBACH4

  !-----------------------------------------------------------------------------------------------------
  !> Type definition for rad memory
  !! 
  !-----------------------------------------------------------------------------------------------------
  TYPE, EXTENDS(t_jsb_memory) :: t_rad_memory

    TYPE(t_jsb_var_real2d)  :: &
      & sw_srf_net,            & !< net shortwave radiation at surface
      & swvis_srf_net,         & !< net shortwave radiation at surface in visible range
      & swnir_srf_net,         & !< net shortwave radiation at surface in NIR range
      & lw_srf_net,            & !< net longwave radiation at surface
      & rad_srf_net              !< net radiation at surface

    TYPE(t_jsb_var_real2d)  :: &        ! REAL(wp) 1D ; implemented for SLM use only, see mo_jsb_var_class 
      & alb_vis,               & !< Albedo of the whole tile in the VIS (visible) range
      & alb_nir,               & !< Albedo of the whole tile in the NIR (near-infrared) range
      & alb_vis_snow,          & !< Albedo of snow in the VIS range
      & alb_nir_snow             !< Albedo of snow in the NIR range                               

    ! Additional variables for soil
    TYPE(t_jsb_var_real2d)  :: &
      & alb_vis_soil,          & !< Albedo in VIS range of soil
                                 ! Initialized from default:
      & alb_nir_soil,          & !< Albedo in VIS range of soil
      & arad_vis_soil,         & !< radiation in VIS range absorbed by the soil
      & arad_nir_soil            !< radiation in NIR range absorbed by the soil

    ! Additional variables for vegetation
    TYPE(t_jsb_var_real2d)  :: &
      & alb_vis_can,           & !< Albedo of canopy (without snow cover) in VIS range.
      & alb_nir_can,           & !< Albedo of canopy (without snow cover) in NIR range. 
                                 !< If use_alb_canopy=.TRUE. taken from bc_land_phys.nc file, else from lctlib-File (AlbedoCanopyNIR)
      & arad_vis_can,          & !< radiation in VIS range absorbed by the canopy      
      & arad_nir_can             !< radiation in NIR range absorbed by the canopy

    TYPE(t_jsb_var_real2d)  :: &
                            appfd                           !< absorbed Photosynthetically Active Photon Flux Density

    TYPE(t_jsb_var_real3d)  :: &   ! DIMENSION(ncanopy)
                            ppfd_sunlit_cl              , & !< Photosynthetically Active Photon Flux Density on sunlit leaves (µmol m-2 s-1)
                            ppfd_shaded_cl              , & !< Photosynthetically Active Photon Flux Density on shaded leaves (µmol m-2 s-1)
                            ppfd_sunlit_daytime_cl      , & !<  average of the previous day
                            ppfd_shaded_daytime_cl      , & !<  average of the previous day
                            ppfd_sunlit_daytime_dacc_cl , & !<  value accumulated over current day
                            ppfd_shaded_daytime_dacc_cl , & !<  value accumulated over current day
                            ppfd_sunlit_tfrac_mavg_cl   , & !< Photosynthetically Active Photon Flux Density on sunlit leaves
                            ppfd_shaded_tfrac_mavg_cl   , & !< Photosynthetically Active Photon Flux Density on shaded leaves
                            ppfd_sunlit_tcnl_mavg_cl    , & !< Photosynthetically Active Photon Flux Density on sunlit leaves
                            ppfd_shaded_tcnl_mavg_cl        !< Photosynthetically Active Photon Flux Density on shaded leaves

    TYPE(t_jsb_var_real3d)  :: &   ! DIMENSION(ncanopy)
                            arad_sunlit_vis_cl   , &
                            arad_shaded_vis_cl   , &
                            fleaf_sunlit_vis_cl  , & 
                            arad_sunlit_nir_cl   , &
                            arad_shaded_nir_cl   , &
                            fleaf_sunlit_nir_cl

  CONTAINS
    PROCEDURE :: Init => Init_rad_memory
  END TYPE t_rad_memory


  CHARACTER(len=*), PARAMETER :: modname = 'mo_rad_memory_class'

CONTAINS


  !-----------------------------------------------------------------------------------------------------
  !> initialize memory (variables) for the process: radiation
  !!
  !! 
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE Init_rad_memory(mem, prefix, suffix, lct_ids, model_id)

!    USE mo_jsb_varlist,       ONLY: t_jsb_varlist, jsb_add_var
!    USE mo_jsb_io,            ONLY: grib_bits, t_cf, t_grib1, t_grib2, TSTEP_CONSTANT, tables
!    USE mo_jsb_grid_class,    ONLY: t_jsb_grid, t_jsb_vgrid
!    USE mo_jsb_grid,          ONLY: Get_grid, Get_vgrid
    USE mo_jsb_model_class,       ONLY: t_jsb_model

    ! 0.1.InOut
    CLASS(t_rad_memory), INTENT(inout), TARGET :: mem     !> radiation memory
    CHARACTER(len=*),    INTENT(in)    :: prefix          !> ?
    CHARACTER(len=*),    INTENT(in)    :: suffix          !> ?
    INTEGER,             INTENT(in)    :: lct_ids !(:)    !> in jsbach4: vector of land cover type IDs
    INTEGER,             INTENT(in)    :: model_id        !> model ID model\%id
    ! 0.2 Local
    TYPE(t_jsb_model), POINTER :: model
!    TYPE(t_jsb_grid),  POINTER :: hgrid                        ! Horizontal grid
!    TYPE(t_jsb_vgrid), POINTER :: surface, canopy_layer        ! Vertical grid
!    INTEGER :: table

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':Init_rad_memory'

    model => Get_model(model_id)

!    table = tables(1)

!    hgrid   => Get_grid(mem%grid_id)
!    surface => Get_vgrid('surface')
!    canopy_layer => Get_vgrid('canopy_layer')


    IF (.NOT. ALLOCATED(mem%sw_srf_net%ptr)) THEN

      CALL slm_add_var(mem, 'sw_srf_net', mem%sw_srf_net, &
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'swvis_srf_net', mem%swvis_srf_net, &
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'swnir_srf_net', mem%swnir_srf_net, &
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'lw_srf_net', mem%lw_srf_net, &
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'rad_srf_net', mem%rad_srf_net, &
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'alb_vis', mem%alb_vis, &
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'alb_nir', mem%alb_nir, &
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'alb_vis_snow', mem%alb_vis_snow, &
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'alb_nir_snow', mem%alb_nir_snow, &
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'alb_vis_soil', mem%alb_vis_soil, &
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'alb_nir_soil', mem%alb_nir_soil, &
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'arad_vis_soil', mem%arad_vis_soil, &
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'arad_nir_soil', mem%arad_nir_soil, &
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'alb_vis_can', mem%alb_vis_can, &
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'alb_nir_can', mem%alb_nir_can, &
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'arad_vis_can', mem%arad_vis_can, &
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'arad_nir_can', mem%arad_nir_can, &
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'appfd', mem%appfd, &
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'ppfd_sunlit_cl', mem%ppfd_sunlit_cl, &
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'ppfd_shaded_cl', mem%ppfd_shaded_cl, &
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'ppfd_sunlit_daytime_cl', mem%ppfd_sunlit_daytime_cl, &
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'ppfd_shaded_daytime_cl', mem%ppfd_shaded_daytime_cl, &
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'ppfd_sunlit_daytime_dacc_cl', mem%ppfd_sunlit_daytime_dacc_cl, &
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'ppfd_shaded_daytime_dacc_cl', mem%ppfd_shaded_daytime_dacc_cl, &
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'ppfd_sunlit_tfrac_mavg_cl', mem%ppfd_sunlit_tfrac_mavg_cl, &
            initval_r=150.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'ppfd_shaded_tfrac_mavg_cl', mem%ppfd_shaded_tfrac_mavg_cl, &
            initval_r=10.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'ppfd_sunlit_tcnl_mavg_cl', mem%ppfd_sunlit_tcnl_mavg_cl, &
            initval_r=150.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'ppfd_shaded_tcnl_mavg_cl', mem%ppfd_shaded_tcnl_mavg_cl, &
            initval_r=10.0_wp) ! initvalue

      CALL slm_add_var(mem, 'arad_sunlit_vis_cl', mem%arad_sunlit_vis_cl, &
            initval_r=0.0_wp) ! initvalue

      CALL slm_add_var(mem, 'arad_shaded_vis_cl', mem%arad_shaded_vis_cl, &
            initval_r=0.0_wp) ! initvalue

      CALL slm_add_var(mem, 'fleaf_sunlit_vis_cl', mem%fleaf_sunlit_vis_cl, &
            initval_r=0.0_wp) ! initvalue

      CALL slm_add_var(mem, 'arad_sunlit_nir_cl', mem%arad_sunlit_nir_cl, &
            initval_r=0.0_wp) ! initvalue

      CALL slm_add_var(mem, 'arad_shaded_nir_cl', mem%arad_shaded_nir_cl, &
            initval_r=0.0_wp) ! initvalue

      CALL slm_add_var(mem, 'fleaf_sunlit_nir_cl', mem%fleaf_sunlit_nir_cl, &
            initval_r=0.0_wp) ! initvalue


    END IF

  END SUBROUTINE Init_rad_memory

END MODULE mo_rad_memory_class
