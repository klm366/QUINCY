!>
!! @brief Contains the memory class for the vegetation processes
!!
!! @author
!!  Reiner Schnur, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Reiner Schnur (2013-04-04)
!!
MODULE mo_veg_memory_class
#ifndef __NO_JSBACH__

  USE mo_kind, ONLY: wp
  USE mo_jsb_memory_class, ONLY: t_jsb_memory
  USE mo_jsb_var_class,    ONLY: t_jsb_var_real2d, t_jsb_var_real3d

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: t_veg_memory, max_no_of_vars

  INTEGER, PARAMETER :: max_no_of_vars = 20

  TYPE, EXTENDS(t_jsb_memory) :: t_veg_memory
    !
    ! Parameters
    TYPE(t_jsb_var_real2d) :: &
      & lai            , &  !< Leaf area index
      & fract_fpc      , &  !< Foliage projected cover (fraction of ground area occupied by the vertical projection of foliage)
      & fract_fpc_max  , &  !< Maximum fraction of tile that can be covered by leaves (assuming infinite lai)
      & fract_forest   , &  !< Forest fraction
      & canopy_cond    , &  !< Canopy conductance (1/resistance)
      & skin_can_depth , &  !< Depth of skin reservoir on canopy
      & root_depth     , &  !< Rooting depth [m]
      & alb_vis_can    , &  !< Albedo of the canopy in the visible range
      & alb_nir_can         !< Albedo of the canopy in the NIR range
    TYPE(t_jsb_var_real3d) :: &
      & root_depth_l     !< Rooted depth per soil layer (until rooting depth) [m]
    REAL(wp), POINTER :: lai_mon(:,:,:) !< Monthly leaf area index (per year or climatology) <nproma x nblks x 14>    !

    ! Water state
    TYPE(t_jsb_var_real2d) :: &
      & wfract_can  , &  !< Wet (skin reservoir) fraction of canopy
      & wsr_can     , &  !< Water content in skin reservoir on canopy [m]
      & sfract_can  , &  !< Snow fraction of canopy
      & wsn_can     , &  !< Water content of snow reservoir onn canopy [m water equivalent]
      & ws_root     , &  !< Water content in root zone of the soil [m]
      & ws_fc_root  , &  !< Water content at field capacity in root zone of the soil [m]
      & ws_pwp_root , &  !< Water content at permanent wilting point in root zone of the soil [m]
      & water_stress     !< Water stress factor (1: no water stress, 0: infinite stress)
  CONTAINS
    PROCEDURE :: Init => Init_veg_memory
  END TYPE t_veg_memory

  CHARACTER(len=*), PARAMETER :: modname = 'mo_veg_memory_class'

CONTAINS

  SUBROUTINE Init_veg_memory(mem, prefix, suffix)

    USE mo_jsb_varlist,       ONLY: t_jsb_varlist, jsb_add_var
    USE mo_jsb_io,            ONLY: grib_bits, nc_bits, t_cf, t_grib1, t_grib2, &
                                    TSTEP_CONSTANT, ZAXIS_SURFACE, tables
    USE mo_jsb_grid_class,    ONLY: t_jsb_grid, t_jsb_vgrid
    USE mo_jsb_grid,          ONLY: Get_grid, Get_vgrid

    CLASS(t_veg_memory), INTENT(inout) :: mem
    CHARACTER(len=*),   INTENT(in)    :: prefix
    CHARACTER(len=*),   INTENT(in)    :: suffix

    TYPE(t_jsb_grid),  POINTER :: hgrid                        ! Horizontal grid
    TYPE(t_jsb_vgrid), POINTER :: surface, soil_w              ! Vertical grids
    INTEGER :: table

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':Init_veg_memory'

    table = tables(1)

    hgrid        => Get_grid(mem%grid_id)

    surface => get_vgrid('surface')
    soil_w  => get_vgrid('soil_depth_water')

      CALL jsb_add_var( mem, prefix//'lai'//suffix, mem%lai,                             &
        & hgrid, surface,                                                 &
        & t_cf('leaf_area_index', '', '', nc_bits),                                &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),             &
        & initval_r=0.0_wp )

      ! Note: extra dimensions (except levels) are not allowed in echam streams
!!$      CALL jsb_add_var( mem, prefix//'lai_mon'//suffix, mem%lai_mon,                     &
!!$        & hgrid, surface,                                                 &
!!$        & t_cf('monthly_leaf_area_index', '', '', nc_bits),                        &
!!$        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),             &
!!$        & extra_dim_n=14, initval_r=0.0_wp, loutput=.FALSE., lrestart=.FALSE. )
      ALLOCATE(mem%lai_mon(hgrid%nproma, hgrid%nblks, 14))

      CALL jsb_add_var( mem, prefix//'fract_fpc'//suffix, mem%fract_fpc,                 &
        & hgrid, surface,                                                 &
        & t_cf('foliage_projected_cover', '', '', nc_bits),                        &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),             &
        & initval_r=0.0_wp )

      CALL jsb_add_var( mem, prefix//'fract_fpc_max'//suffix, mem%fract_fpc_max,         &
        & hgrid, surface,                                                 &
        & t_cf('maximum_fpc', '', '', nc_bits),                                    &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),             &
        & initval_r=0.0_wp )

      CALL jsb_add_var( mem, prefix//'fract_forest'//suffix, mem%fract_forest,           &
        & hgrid, surface,                                                 &
        & t_cf('forest_fraction', '', '', nc_bits),                                &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),             &
        & initval_r=0.0_wp )

      CALL jsb_add_var( mem, prefix//'canopy_cond'//suffix, mem%canopy_cond,             &
        & hgrid, surface,                                                 &
        & t_cf('canopy_conductance', '', '', nc_bits),                             &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),             &
        & initval_r=0.0_wp )

      CALL jsb_add_var( mem, prefix//'skin_can_depth'//suffix, mem%skin_can_depth,       &
        & hgrid, surface,                                                 &
        & t_cf('canopy_skin_depth', 'm', '', nc_bits),                             &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),             &
        & initval_r=0.0_wp )

      CALL jsb_add_var( mem, prefix//'root_depth'//suffix, mem%root_depth,               &
        & hgrid, surface,                                                 &
        & t_cf('rooting_depth', '', '', nc_bits),                                  &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),             &
        & initval_r=0.0_wp, isteptype=TSTEP_CONSTANT )

      CALL jsb_add_var( mem, prefix//'root_depth_l'//suffix, mem%root_depth_l,           &
        & hgrid, soil_w,                                                                      &
        & t_cf('root_depth_layer', 'm', '', nc_bits),                              &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),             &
        & initval_r=0.0_wp, isteptype=TSTEP_CONSTANT )

      CALL jsb_add_var( mem, prefix//'alb_vis_can'//suffix, mem%alb_vis_can,             &
        & hgrid, surface,                                                 &
        & t_cf('canopy_albedo_visible', '', '', nc_bits),                          &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),             &
        & initval_r=0.0_wp )

      CALL jsb_add_var( mem, prefix//'alb_nir_can'//suffix, mem%alb_nir_can,             &
        & hgrid, surface,                                                 &
        & t_cf('canopy_albedo_nir', '', '', nc_bits),                              &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),             &
        & initval_r=0.0_wp )

      CALL jsb_add_var( mem, prefix//'wfract_can'//suffix, mem%wfract_can,                         &
        & hgrid, surface,                                                           &
        & t_cf('wet_fraction_canopy', '', '', nc_bits),                                      &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                       &
        & initval_r=0.0_wp )

      CALL jsb_add_var( mem, prefix//'wsr_can'//suffix, mem%wsr_can,                               &
        & hgrid, surface,                                                           &
        & t_cf('skin_reservoir_canopy', 'm', '', nc_bits),                                   &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                       &
        & initval_r=0.0_wp )

      CALL jsb_add_var( mem, prefix//'sfract_can'//suffix, mem%sfract_can,                         &
        & hgrid, surface,                                                           &
        & t_cf('snow_fraction_canopy', '', '', nc_bits),                                     &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                       &
        & initval_r=0.0_wp )

      CALL jsb_add_var( mem, prefix//'wsn_can'//suffix, mem%wsn_can,                               &
        & hgrid, surface,                                                           &
        & t_cf('snow_water_equivalent_canopy', 'm', '', nc_bits),                            &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                       &
        & initval_r=0.0_wp )

      CALL jsb_add_var( mem, prefix//'ws_root'//suffix, mem%ws_root,                               &
        & hgrid, surface,                                                           &
        & t_cf('soil_moisture_root', 'm', '', nc_bits),                                      &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                       &
        & initval_r=0.0_wp )

      CALL jsb_add_var( mem, prefix//'ws_fc_root'//suffix, mem%ws_fc_root,                         &
        & hgrid, surface,                                                           &
        & t_cf('field_capacity_root', 'm', '', nc_bits),                                     &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                       &
        & initval_r=0.0_wp )

      CALL jsb_add_var( mem, prefix//'ws_pwp_root'//suffix, mem%ws_pwp_root,                       &
        & hgrid, surface,                                                           &
        & t_cf('wilting_root', 'm', '', nc_bits),                                            &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                       &
        & initval_r=0.0_wp )

      CALL jsb_add_var( mem, prefix//'water_stress'//suffix, mem%water_stress,                     &
        & hgrid, surface,                                                           &
        & t_cf('water_stress_factor', '', '', nc_bits),                                      &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                       &
        & initval_r=0.0_wp )

  END SUBROUTINE Init_veg_memory

#endif
END MODULE mo_veg_memory_class
