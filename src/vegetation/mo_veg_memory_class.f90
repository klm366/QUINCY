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

  USE mo_jsb_model_class,  ONLY: t_jsb_model
  USE mo_jsb_class,        ONLY: Get_model
  USE mo_jsb_memory_class, ONLY: t_jsb_memory
  USE mo_jsb_var_class,    ONLY: t_jsb_var_real2d, t_jsb_var_real3d

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: t_veg_memory, max_no_of_vars

  INTEGER, PARAMETER :: max_no_of_vars = 20

  TYPE, EXTENDS(t_jsb_memory) :: t_veg_memory
    !
    ! Parameters
    TYPE(t_jsb_var_real2d) ::                      &
      & lai,                 &  !< Leaf area index
      & fract_fpc,           &  !< Foliage projected cover (fraction of ground area occupied by the vertical projection of foliage)
      & fract_fpc_max,       &  !< Maximum fraction of tile that can be covered by leaves (assuming infinite lai)
      & fract_forest,        &  !< Forest fraction
      & canopy_cond,         & ! Canopy conductance without water limit  [m/s] (mean value)
      & canopy_cond_limited, & ! Canopy conductance water limit [m/s] (mean value)
      & skin_can_depth,      &  !< Depth of skin reservoir on canopy
      & root_depth,          &  !< Rooting depth [m]
      & alb_vis_can,         &  !< Albedo of the canopy in the visible range
      & alb_nir_can             !< Albedo of the canopy in the NIR range
    TYPE(t_jsb_var_real3d) ::               &
      & root_depth_l     !< Rooted depth per soil layer (until rooting depth) [m]
    REAL(wp), POINTER :: lai_mon(:,:,:) !< Monthly leaf area index (per year or climatology) <nproma x nblks x 14>    !

  CONTAINS
    PROCEDURE :: Init => Init_veg_memory
  END TYPE t_veg_memory

  CHARACTER(len=*), PARAMETER :: modname = 'mo_veg_memory_class'

CONTAINS

  SUBROUTINE Init_veg_memory(mem, prefix, suffix, lct_ids, model_id)

    USE mo_jsb_varlist,       ONLY: t_jsb_varlist, jsb_add_var
    USE mo_jsb_io,            ONLY: grib_bits, t_cf, t_grib1, t_grib2, &
                                    TSTEP_CONSTANT, tables
    USE mo_jsb_grid_class,    ONLY: t_jsb_grid, t_jsb_vgrid
    USE mo_jsb_grid,          ONLY: Get_grid, Get_vgrid

    CLASS(t_veg_memory), INTENT(inout), TARGET :: mem
    CHARACTER(len=*),    INTENT(in)    :: prefix
    CHARACTER(len=*),    INTENT(in)    :: suffix
    INTEGER,             INTENT(in)    :: lct_ids(:)
    INTEGER,             INTENT(in)    :: model_id

    TYPE(t_jsb_model), POINTER :: model
    TYPE(t_jsb_grid),  POINTER :: hgrid                        ! Horizontal grid
    TYPE(t_jsb_vgrid), POINTER :: surface, soil_w              ! Vertical grids
    INTEGER :: table

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':Init_veg_memory'

    model => Get_model(model_id)

    table = tables(1)

    hgrid        => Get_grid(mem%grid_id)

    surface => get_vgrid('surface')
    soil_w  => get_vgrid('soil_depth_water')

    IF (ASSOCIATED(mem%lai%ptr)) RETURN

    CALL jsb_add_var( mem, 'lai', mem%lai,                             &
      & hgrid, surface,                                                 &
      & t_cf('leaf_area_index', '', ''),                                         &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),             &
      & prefix, suffix,                                                      &
      & initval_r=0.0_wp )

      ! Note: extra dimensions (except levels) are not allowed in echam streams
!!$      CALL jsb_add_var( mem, 'lai_mon', mem%lai_mon,                     &
!!$        & hgrid, surface,                                                 &
!!$        & t_cf('monthly_leaf_area_index', '', ''),                                 &
!!$        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),             &
!!$        & extra_dim_n=14, initval_r=0.0_wp, loutput=.FALSE., lrestart=.FALSE. )
    ALLOCATE(mem%lai_mon(hgrid%nproma, hgrid%nblks, 14))

    CALL jsb_add_var( mem, 'fract_fpc', mem%fract_fpc,                 &
      & hgrid, surface,                                                 &
      & t_cf('foliage_projected_cover', '', ''),                                 &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),             &
      & prefix, suffix,                                                      &
      & initval_r=0.0_wp )

    CALL jsb_add_var( mem, 'fract_fpc_max', mem%fract_fpc_max,         &
      & hgrid, surface,                                                 &
      & t_cf('maximum_fpc', '', ''),                                             &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),             &
      & prefix, suffix,                                                      &
      & initval_r=0.0_wp )

    CALL jsb_add_var( mem, 'fract_forest', mem%fract_forest,           &
      & hgrid, surface,                                                 &
      & t_cf('forest_fraction', '', ''),                                         &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),             &
      & prefix, suffix,                                                      &
      & initval_r=0.0_wp )

    CALL jsb_add_var( mem, 'canopy_cond', mem%canopy_cond,             &
      & hgrid, surface,                                                 &
      & t_cf('canopy_conductance', '', ''),                                      &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),             &
      & prefix, suffix,                                                      &
      & initval_r=0.0_wp )

    CALL jsb_add_var( mem, 'canopy_cond_limited', mem%canopy_cond_limited,     &
      & hgrid, surface,                                                 &
      & t_cf('canopy_conductance_limited', '', ''),                                      &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),             &
      & prefix, suffix,                                                      &
      & initval_r=0.0_wp )

    CALL jsb_add_var( mem, 'skin_can_depth', mem%skin_can_depth,       &
      & hgrid, surface,                                                 &
      & t_cf('canopy_skin_depth', 'm', ''),                                      &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),             &
      & prefix, suffix,                                                      &
      & initval_r=0.0_wp )

    CALL jsb_add_var( mem, 'root_depth', mem%root_depth,               &
      & hgrid, surface,                                                 &
      & t_cf('rooting_depth', '', ''),                                  &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),             &
      & prefix, suffix,                                                      &
      & initval_r=0.0_wp, isteptype=TSTEP_CONSTANT )

    CALL jsb_add_var( mem, 'root_depth_l', mem%root_depth_l,           &
      & hgrid, soil_w,                                                                      &
      & t_cf('root_depth_layer', 'm', ''),                              &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),             &
      & prefix, suffix,                                                      &
      & initval_r=0.0_wp, isteptype=TSTEP_CONSTANT )

    CALL jsb_add_var( mem, 'alb_vis_can', mem%alb_vis_can,             &
      & hgrid, surface,                                                 &
      & t_cf('canopy_albedo_visible', '', ''),                          &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),             &
      & prefix, suffix,                                                      &
      & initval_r=0.0_wp )

    CALL jsb_add_var( mem, 'alb_nir_can', mem%alb_nir_can,             &
      & hgrid, surface,                                                 &
      & t_cf('canopy_albedo_nir', '', ''),                              &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),             &
      & prefix, suffix,                                                      &
      & initval_r=0.0_wp )


  END SUBROUTINE Init_veg_memory

#endif
END MODULE mo_veg_memory_class
