!>
!! @brief Contains the memory class for the soil processes.
!!
!! @author
!!  Reiner Schnur, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Reiner Schnur (2013-04-04)
!!
MODULE mo_soil_memory_class
#ifndef __NO_JSBACH__

  USE mo_kind, ONLY: wp

  USE mo_jsb_model_class,  ONLY: t_jsb_model
  USE mo_jsb_class,        ONLY: Get_model
  USE mo_jsb_memory_class, ONLY: t_jsb_memory, t_jsb_param
  USE mo_jsb_var_class,    ONLY: t_jsb_var_real2d, t_jsb_var_real3d

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: t_soil_memory, max_no_of_vars
!!$            copy_soil_mem, add_soil_mem

  INTEGER, PARAMETER :: max_no_of_vars = 40

!!$  TYPE, EXTENDS(t_jsb_param) :: t_soil_param
!!$    TYPE(t_jsb_var_real3d) :: &
!!$      soil_depth_l               !< Depth of each soil layer that can be saturated with water (until bedrock) [m]
!!$    TYPE(t_jsb_var_real2d) :: &
!!$      soil_depth          ,    & !< Soil depth derived from textures (bedrock) (for hydrology) [m]
!!$      vpor                ,    & !< Volumetric soil porosity []
!!$      psdi                ,    & !< Soil pore size distribution index []
!!$      bclapp              ,    & !< Exponent B in Clapp and Hornberger
!!$      mpot                ,    & !< Matrix potential [m]
!!$      shc                 ,    & !< Saturated hydraulic conductivity [m/s]
!!$      vhcap               ,    & !< Volumetric heat capacity of dry soil [J m-3 K-1]
!!$      hdiff               ,    & !< Thermal diffusivity of soil [m2 s-1]
!!$      alb_vis_soil        ,    & !< Albedo of the soil in the visible range
!!$      alb_nir_soil               !< Albedo of the soil in the NIR range
!!$  END type t_soil_param

  TYPE, EXTENDS(t_jsb_memory) :: t_soil_memory
    !
    ! Parameters
    TYPE(t_jsb_var_real3d) :: &
      soil_depth_l               !< Depth of each soil layer that can be saturated with water (until bedrock) [m]
    TYPE(t_jsb_var_real2d) :: &
      soil_depth          ,    & !< Soil depth derived from textures (bedrock) (for hydrology) [m]
      vpor                ,    & !< Volumetric soil porosity []
      psdi                ,    & !< Soil pore size distribution index []
      bclapp              ,    & !< Exponent B in Clapp and Hornberger
      mpot                ,    & !< Matrix potential [m]
      shc                 ,    & !< Saturated hydraulic conductivity [m/s]
      vhcap               ,    & !< Volumetric heat capacity of dry soil [J m-3 K-1]
      hdiff               ,    & !< Thermal diffusivity of soil [m2 s-1]
      alb_vis_soil        ,    & !< Albedo of the soil in the visible range
      alb_nir_soil        ,    & !< Albedo of the soil in the NIR range
      oro_stddev                 !< Standard deviation of orography
    !
    ! Energy state
    TYPE(t_jsb_var_real3d) :: &
      t_soil, &
      t_soil_acoef  , &
      t_soil_bcoef
    TYPE(t_jsb_var_real2d) :: &
      hcap_grnd     , &  !< Heat capacity of the ground [J m-2 K-1]
      hcap_grnd_old      !< Heat capacity of the ground [J m-2 K-1] (old)
    !
    ! Energy fluxes
    TYPE(t_jsb_var_real2d) :: &
      grnd_hflx        ,  &  !< ground heat flux
      grnd_hflx_old          !< ground heat flux (old)

  CONTAINS
    PROCEDURE :: Init => Init_soil_memory
  END TYPE t_soil_memory

  CHARACTER(len=*), PARAMETER :: modname = 'mo_soil_memory_class'

CONTAINS

  SUBROUTINE Init_soil_memory(mem, prefix, suffix, lct_ids, model_id)

    !USE mo_jsb_process
    USE mo_jsb_config_class
    USE mo_soil_config_class
    USE mo_jsb_varlist,       ONLY: jsb_add_var
    USE mo_jsb_io,            ONLY: grib_bits, t_cf, t_grib1, t_grib2, &
                                    TSTEP_CONSTANT, tables
    USE mo_jsb_grid_class,    ONLY: t_jsb_grid, t_jsb_vgrid
    USE mo_jsb_grid,          ONLY: Get_grid, Get_vgrid

    CLASS(t_soil_memory), INTENT(inout) :: mem
    CHARACTER(len=*),     INTENT(in)    :: prefix
    CHARACTER(len=*),     INTENT(in)    :: suffix
    INTEGER,              INTENT(in)    :: lct_ids(:)
    INTEGER,              INTENT(in)    :: model_id

    TYPE(t_jsb_model), POINTER :: model
    TYPE(t_jsb_grid),  POINTER :: hgrid                        ! Horizontal grid
    TYPE(t_jsb_vgrid), POINTER :: surface, soil_e, soil_w      ! Vertical grids
    INTEGER :: table

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':Init_soil_memory'

    model => Get_model(model_id)

    table = tables(1)

    hgrid => Get_grid(mem%grid_id)

    surface => Get_vgrid('surface')
    soil_e  => Get_vgrid('soil_depth_energy')
    soil_w  => Get_vgrid('soil_depth_water')

    IF (ASSOCIATED(mem%soil_depth%ptr)) RETURN

    CALL jsb_add_var( mem, 'soil_depth', mem%soil_depth,                         &
      & hgrid, surface,                                                                    &
      & t_cf('soil_depth', 'm', 'Depth of soil until bedrock'),                            &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
      & prefix, suffix,                                                      &
      & initval_r=0.0_wp, isteptype=TSTEP_CONSTANT )

    CALL jsb_add_var( mem, 'soil_depth_l', mem%soil_depth_l,                     &
      & hgrid, soil_w,                                                                     &
      & t_cf('soil_depth_layer', 'm', ''),                                                 &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
      & prefix, suffix,                                                      &
      & initval_r=0.0_wp, isteptype=TSTEP_CONSTANT )

    CALL jsb_add_var( mem, 'soil_por', mem%vpor,                                 &
      & hgrid, surface,                                                                    &
      & t_cf('vol_soil_porosity', '', ''),                                                 &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
      & prefix, suffix,                                                      &
      & initval_r=0.0_wp, isteptype=TSTEP_CONSTANT )

    CALL jsb_add_var( mem, 'pore_size', mem%psdi,                                &
      & hgrid, surface,                                                                    &
      & t_cf('pore_size_dist_index', '', ''),                                              &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
      & prefix, suffix,                                                      &
      & initval_r=0.0_wp, isteptype=TSTEP_CONSTANT )

    CALL jsb_add_var( mem, 'bclapp', mem%bclapp,                                 &
      & hgrid, surface,                                                                    &
      & t_cf('exp_clapp_hornberger', '', ''),                                              &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
      & prefix, suffix,                                                      &
      & initval_r=0.0_wp, isteptype=TSTEP_CONSTANT )

    CALL jsb_add_var( mem, 'matrix_pot', mem%mpot,                               &
      & hgrid, surface,                                                                    &
      & t_cf('matrix_potential', 'm', ''),                                                 &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
      & prefix, suffix,                                                      &
      & initval_r=0.0_wp, isteptype=TSTEP_CONSTANT )

    CALL jsb_add_var( mem, 'sat_hyd_cond', mem%shc,                              &
      & hgrid, surface,                                                                    &
      & t_cf('sat_hydraulic_cond', 'm s-1', ''),                                           &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
      & prefix, suffix,                                                      &
      & initval_r=0.0_wp, isteptype=TSTEP_CONSTANT )

    CALL jsb_add_var( mem, 'vol_heat_cap', mem%vhcap,                            &
      & hgrid, surface,                                                                    &
      & t_cf('vol_heat_capacity_soil', 'J m-3 K-1', ''),                                   &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
      & prefix, suffix,                                                      &
      & initval_r=0.0_wp, isteptype=TSTEP_CONSTANT )

    CALL jsb_add_var( mem, 'heat_diff', mem%hdiff,                               &
      & hgrid, surface,                                                                    &
      & t_cf('heat_diffusivity_soil', 'm2 s-1', ''),                                       &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
      & prefix, suffix,                                                      &
      & initval_r=0.0_wp, isteptype=TSTEP_CONSTANT )

    CALL jsb_add_var( mem, 'alb_vis_soil', mem%alb_vis_soil,                     &
      & hgrid, surface,                                                                    &
      & t_cf('soil_albedo_visible ', '', ''),                                              &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
      & prefix, suffix,                                                      &
      & initval_r=0.0_wp )

    CALL jsb_add_var( mem, 'alb_nir_soil', mem%alb_nir_soil,                     &
      & hgrid, surface,                                                                    &
      & t_cf('soil_albedo_nir', '', ''),                                                   &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
      & prefix, suffix,                                                      &
      & initval_r=0.0_wp )

    CALL jsb_add_var( mem, 'oro_stddev', mem%oro_stddev,                     &
      & hgrid, surface,                                                                    &
      & t_cf('orography_standard_deviation', '', ''),                        &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
      & prefix, suffix,                                                      &
      & initval_r=0.0_wp, isteptype=TSTEP_CONSTANT )

    CALL jsb_add_var( mem, 't_soil', mem%t_soil,                                 &
      & hgrid, soil_e,                                                                     &
      & t_cf('soil_temperature', 'K', ''),                                                 &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
      & prefix, suffix,                                                      &
      & initval_r=0.0_wp )

    CALL jsb_add_var( mem, 't_soil_aceof', mem%t_soil_acoef,                     &
      & hgrid, soil_e,                                                                     &
      & t_cf('soil_temperature_acoef', '', ''),                                            &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
      & prefix, suffix,                                                      &
      & initval_r=0.0_wp )

    CALL jsb_add_var( mem, 't_soil_bceof', mem%t_soil_bcoef,                     &
      & hgrid, soil_e,                                                                     &
      & t_cf('soil_temperature_bcoef', '', ''),                                            &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
      & prefix, suffix,                                                      &
      & initval_r=0.0_wp )

    CALL jsb_add_var( mem, 'heat_cap_grnd', mem%hcap_grnd,                       &
      & hgrid, surface,                                                                    &
      & t_cf('heat_capacity_ground', 'J m-2 K-1', 'Ground heat capacity'),                 &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
      & prefix, suffix,                                                      &
      & initval_r=0.0_wp )

    CALL jsb_add_var( mem, 'heat_cap_grnd_old', mem%hcap_grnd_old,&
      & hgrid, surface,                                                                    &
      & t_cf('heat_capacity_ground_old', 'J m-2 K-1', 'Ground heat capacity (old)'),          &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
      & prefix, suffix,                                                      &
      & initval_r=0.0_wp )

    CALL jsb_add_var( mem, 'grnd_hflx', mem%grnd_hflx,                            &
      & hgrid, surface,                                                                    &
      & t_cf('grnd_hflx', '', 'Ground heat flux'),                                         &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
      & prefix, suffix,                                                      &
      & initval_r=0.0_wp )

    CALL jsb_add_var( mem, 'grnd_hflx_old', mem%grnd_hflx_old,                    &
      & hgrid, surface,                                                                    &
      & t_cf('grnd_hflx_old', '', 'Ground heat flux (old)'),                               &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
      & prefix, suffix,                                                      &
      & initval_r=0.0_wp )

    END SUBROUTINE Init_soil_memory

#endif
END MODULE mo_soil_memory_class
