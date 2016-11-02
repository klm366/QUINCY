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
      alb_nir_soil               !< Albedo of the soil in the NIR range
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
    ! Water state
    TYPE(t_jsb_var_real2d) :: &
      wfract_soil   ,   & !< Wet (skin reservoir) fraction of soil []
      wsr_soil      ,   & !< Water content in skin reservoir on soil [m]
      sfract_soil   ,   & !< Snow fraction of soil []
      wsn_soil      ,   & !< Water content of snow reservoir of soil [m water equivalent]
      overflow      ,   & !< Water content of reservoir for overflown soil water (terraplanet) [m] 
      ws                  !< Water content in the whole soil column [m]
    TYPE(t_jsb_var_real3d) :: &
      ws_l          ,   & !< Water content in soil layers [m]
      ws_sat_l      ,   & !< Water content of soil layers at saturation (derived from soil porosity) [m]
      ws_fc_l       ,   & !< Water content of soil layers at field capacity (derived from volumetric field capacity)[m]
      ws_pwp_l            !< Water content of soil layers at permanent wilting point (derived from vol. perm wilt point) [m]
    !
    ! Energy fluxes
    TYPE(t_jsb_var_real2d) :: &
      grnd_hflx        ,  &  !< ground heat flux
      grnd_hflx_old          !< ground heat flux (old)
    !
    ! Water fluxes
    TYPE(t_jsb_var_real2d) :: &
      & runoff         , & !< surface runoff
      & drainage       , & !< drainage
      & pme_glacier    , & !< P-E on glaciers
      & discharge      , & !< discharge (local)
      & discharge_ocean, & !< discharge to the ocean
      & internal_drain,  & !< internal drain
      & water_flux,      & !< total water flux
      & water_content,   & !< total water content
      & water_budget       !< Water budget error
  CONTAINS
    PROCEDURE :: Init => Init_soil_memory
  END TYPE t_soil_memory

  ! Note: Would be nicer to be able to overload the "=" operator, but that would be more
  ! complicated because we only want to copy parts of the arrays (ics:ice,iblk)
!!$  INTERFACE copy_soil_mem
!!$    MODULE PROCEDURE copy_t_soil_energy_state
!!$  END INTERFACE copy_soil_mem
!!$
!!$  INTERFACE add_soil_mem
!!$    MODULE PROCEDURE add_t_soil_energy_state
!!$  END INTERFACE add_soil_mem

  CHARACTER(len=*), PARAMETER :: modname = 'mo_soil_memory_class'

CONTAINS 

  SUBROUTINE Init_soil_memory(mem, prefix, suffix)

    USE mo_jsb_varlist,       ONLY: jsb_add_var
    USE mo_jsb_io,            ONLY: grib_bits, t_cf, t_grib1, t_grib2, &
                                    TSTEP_CONSTANT, tables
    USE mo_jsb_grid_class,    ONLY: t_jsb_grid, t_jsb_vgrid
    USE mo_jsb_grid,          ONLY: Get_grid, Get_vgrid

    CLASS(t_soil_memory), INTENT(inout) :: mem
    CHARACTER(len=*),     INTENT(in)    :: prefix
    CHARACTER(len=*),     INTENT(in)    :: suffix

    TYPE(t_jsb_grid),  POINTER :: hgrid                        ! Horizontal grid
    TYPE(t_jsb_vgrid), POINTER :: surface, soil_e, soil_w      ! Vertical grids
    INTEGER :: table

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':Init_soil_memory'

    table = tables(1)

    hgrid => Get_grid(mem%grid_id)

    surface => Get_vgrid('surface')
    soil_e  => Get_vgrid('soil_depth_energy')
    soil_w  => Get_vgrid('soil_depth_water')

      CALL jsb_add_var( mem, prefix//'soil_depth'//suffix, mem%soil_depth,                         &
        & hgrid, surface,                                                                    &
        & t_cf('soil_depth', 'm', 'Depth of soil until bedrock'),                            &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
        & initval_r=0.0_wp, isteptype=TSTEP_CONSTANT )

      CALL jsb_add_var( mem, prefix//'soil_depth_l'//suffix, mem%soil_depth_l,                     &
        & hgrid, soil_w,                                                                     &
        & t_cf('soil_depth_layer', 'm', ''),                                                 &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
        & initval_r=0.0_wp, isteptype=TSTEP_CONSTANT )

      CALL jsb_add_var( mem, prefix//'soil_por'//suffix, mem%vpor,                                 &
        & hgrid, surface,                                                                    &
        & t_cf('vol_soil_porosity', '', ''),                                                 &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
        & initval_r=0.0_wp, isteptype=TSTEP_CONSTANT )

      CALL jsb_add_var( mem, prefix//'pore_size'//suffix, mem%psdi,                                &
        & hgrid, surface,                                                                    &
        & t_cf('pore_size_dist_index', '', ''),                                              &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
        & initval_r=0.0_wp, isteptype=TSTEP_CONSTANT )

      CALL jsb_add_var( mem, prefix//'bclapp'//suffix, mem%bclapp,                                 &
        & hgrid, surface,                                                                    &
        & t_cf('exp_clapp_hornberger', '', ''),                                              &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
        & initval_r=0.0_wp, isteptype=TSTEP_CONSTANT )

      CALL jsb_add_var( mem, prefix//'matrix_pot'//suffix, mem%mpot,                               &
        & hgrid, surface,                                                                    &
        & t_cf('matrix_potential', 'm', ''),                                                 &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
        & initval_r=0.0_wp, isteptype=TSTEP_CONSTANT )

      CALL jsb_add_var( mem, prefix//'sat_hyd_cond'//suffix, mem%shc,                              &
        & hgrid, surface,                                                                    &
        & t_cf('sat_hydraulic_cond', 'm s-1', ''),                                           &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
        & initval_r=0.0_wp, isteptype=TSTEP_CONSTANT )

      CALL jsb_add_var( mem, prefix//'vol_heat_cap'//suffix, mem%vhcap,                            &
        & hgrid, surface,                                                                    &
        & t_cf('vol_heat_capacity_soil', 'J m-3 K-1', ''),                                   &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
        & initval_r=0.0_wp, isteptype=TSTEP_CONSTANT )

      CALL jsb_add_var( mem, prefix//'heat_diff'//suffix, mem%hdiff,                               &
        & hgrid, surface,                                                                    &
        & t_cf('heat_diffusivity_soil', 'm2 s-1', ''),                                       &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
        & initval_r=0.0_wp, isteptype=TSTEP_CONSTANT )

      CALL jsb_add_var( mem, prefix//'alb_vis_soil'//suffix, mem%alb_vis_soil,                     &
        & hgrid, surface,                                                                    &
        & t_cf('soil_albedo_visible ', '', ''),                                              &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
        & initval_r=0.0_wp )

      CALL jsb_add_var( mem, prefix//'alb_nir_soil'//suffix, mem%alb_nir_soil,                     &
        & hgrid, surface,                                                                    &
        & t_cf('soil_albedo_nir', '', ''),                                                   &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
        & initval_r=0.0_wp )

      CALL jsb_add_var( mem, prefix//'t_soil'//suffix, mem%t_soil,                                 &
        & hgrid, soil_e,                                                                     &
        & t_cf('soil_temperature', 'K', ''),                                                 &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
        & initval_r=0.0_wp )

      CALL jsb_add_var( mem, prefix//'t_soil_aceof'//suffix, mem%t_soil_acoef,                     &
        & hgrid, soil_e,                                                                     &
        & t_cf('soil_temperature_acoef', '', ''),                                            &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
        & initval_r=0.0_wp )

      CALL jsb_add_var( mem, prefix//'t_soil_bceof'//suffix, mem%t_soil_bcoef,                     &
        & hgrid, soil_e,                                                                     &
        & t_cf('soil_temperature_bcoef', '', ''),                                            &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
        & initval_r=0.0_wp )

      CALL jsb_add_var( mem, prefix//'heat_cap_grnd'//suffix, mem%hcap_grnd,                       &
        & hgrid, surface,                                                                    &
        & t_cf('heat_capacity_ground', 'J m-2 K-1', 'Ground heat capacity'),                 &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
        & initval_r=0.0_wp )

      CALL jsb_add_var( mem, prefix//'heat_cap_grnd_old'//suffix, mem%hcap_grnd_old,&
        & hgrid, surface,                                                                    &
        & t_cf('heat_capacity_ground_old', 'J m-2 K-1', 'Ground heat capacity (old)'),          &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
        & initval_r=0.0_wp )

      CALL jsb_add_var( mem, prefix//'grnd_hflx'//suffix, mem%grnd_hflx,                            &
        & hgrid, surface,                                                                    &
        & t_cf('grnd_hflx', '', 'Ground heat flux'),                                         &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
        & initval_r=0.0_wp )

      CALL jsb_add_var( mem, prefix//'grnd_hflx_old'//suffix, mem%grnd_hflx_old,                    &
        & hgrid, surface,                                                                    &
        & t_cf('grnd_hflx_old', '', 'Ground heat flux (old)'),                               &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
        & initval_r=0.0_wp )

      CALL jsb_add_var( mem, prefix//'wfract_soil'//suffix, mem%wfract_soil,                       &
        & hgrid, surface,                                                                    &
        & t_cf('wet_fraction_soil', '', ''),                                                 &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
        & initval_r=0.0_wp )

      CALL jsb_add_var( mem, prefix//'wsr_soil'//suffix, mem%wsr_soil,                             &
        & hgrid, surface,                                                                    &
        & t_cf('skin_reservoir_soil', 'm', ''),                                              &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
        & initval_r=0.0_wp )

      CALL jsb_add_var( mem, prefix//'sfract_soil'//suffix, mem%sfract_soil,                       &
        & hgrid, surface,                                                                    &
        & t_cf('snow_fraction_soil', '', ''),                                                &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
        & initval_r=0.0_wp )

      CALL jsb_add_var( mem, prefix//'wsn_soil'//suffix, mem%wsn_soil,                             &
        & hgrid, surface,                                                                    &
        & t_cf('snow_water_equivalent_soil', 'm', ''),                                       &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
        & initval_r=0.0_wp )

      CALL jsb_add_var( mem, prefix//'overflow'//suffix, mem%overflow,                             &
        & hgrid, surface,                                                                    &
        & t_cf('overflown_soil_moisture', 'm', ''),                                       &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
        & initval_r=0.0_wp )

      CALL jsb_add_var( mem, prefix//'ws_l'//suffix, mem%ws_l,                                     &
        & hgrid, soil_w,                                                                     &
        & t_cf('soil_moisture', 'm', ''),                                                    &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
        & initval_r=0.0_wp )

      CALL jsb_add_var( mem, prefix//'ws_fc_l'//suffix, mem%ws_fc_l,                               &
        & hgrid, soil_w,                                                                     &
        & t_cf('field_capacity', 'm', ''),                                                   &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
        & initval_r=0.0_wp )

      CALL jsb_add_var( mem, prefix//'ws_pwp_l'//suffix, mem%ws_pwp_l,                             &
        & hgrid, soil_w,                                                                     &
        & t_cf('permanent_wilting_point', 'm', ''),                                          &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
        & initval_r=0.0_wp )

      CALL jsb_add_var( mem, prefix//'ws_sat_l'//suffix, mem%ws_sat_l,                             &
        & hgrid, soil_w,                                                                     &
        & t_cf('water holding capacity', 'm', ''),                                           &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
        & initval_r=0.0_wp )

      CALL jsb_add_var( mem, prefix//'ws'//suffix, mem%ws,                                         &
        & hgrid, surface,                                                                    &
        & t_cf('soil_moisture', 'm', ''),                                                    &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
        & initval_r=0.0_wp )

      CALL jsb_add_var( mem, prefix//'runoff'//suffix, mem%runoff,                   &
        & hgrid, surface,                                                                    &
        & t_cf('surface_runoff', 'kg m-2 s-1', 'Surface runoff'),                            &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                 &
        & initval_r=0.0_wp )

      CALL jsb_add_var( mem, prefix//'drainage'//suffix, mem%drainage,               &
        & hgrid, surface,                                                                    &
        & t_cf('drainage', 'kg m-2 s-1', 'drainage'),                                        &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                 &
        & initval_r=0.0_wp )

      CALL jsb_add_var( mem, prefix//'pme_glacier'//suffix, mem%pme_glacier,         &
        & hgrid, surface,                                                                    &
        & t_cf('pme_glacier', 'kg m-2 s-1', 'P-E on glaciers'),                              &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                 &
        & initval_r=0.0_wp )

      CALL jsb_add_var( mem, prefix//'discharge'//suffix, mem%discharge,             &
        & hgrid, surface,                                                                    &
        & t_cf('discharge', 'm3 s-1', 'local discharge'),                                    &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                 &
        & initval_r=0.0_wp )

      CALL jsb_add_var( mem, prefix//'discharge_ocean'//suffix, mem%discharge_ocean, &
        & hgrid, surface,                                                                    &
        & t_cf('discharge_ocean', 'm3 s-1', 'discharge to the ocean'),                       &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                 &
        & initval_r=0.0_wp )

      CALL jsb_add_var( mem, prefix//'internal_drain'//suffix, mem%internal_drain,           &
        & hgrid, surface,                                                                    &
        & t_cf('internal_drain', 'm3 s-1', 'internal drainage'),                             &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                 &
        & initval_r=0.0_wp )

      CALL jsb_add_var( mem, prefix//'water_flux'//suffix, mem%water_flux,                   &
        & hgrid, surface,                                                                    &
        & t_cf('water_flux', 'm3 s-1', 'total water flux'),                                  &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                 &
        & initval_r=0.0_wp )

      CALL jsb_add_var( mem, prefix//'water_content'//suffix, mem%water_content,             &
        & hgrid, surface,                                                                    &
        & t_cf('water_content', 'm3', 'total water content'),                                &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                 &
        & initval_r=0.0_wp )

      CALL jsb_add_var( mem, prefix//'water_budget'//suffix, mem%water_budget,               &
        & hgrid, surface,                                                                    &
        & t_cf('water_budget', 'm3', 'water budget error'),                                  &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                 &
        & initval_r=0.0_wp )


    END SUBROUTINE Init_soil_memory

!!$  SUBROUTINE copy_t_soil_energy_state(in, out, iblk, ics, ice)
!!$
!!$    TYPE(t_soil_energy_state), INTENT(in)    :: in
!!$    TYPE(t_soil_energy_state), INTENT(inout) :: out
!!$    INTEGER,                   INTENT(in)    :: iblk, ics, ice
!!$
!!$    out%t_soil        (ics:ice, :, iblk) = in%t_soil        (ics:ice, :, iblk)
!!$    out%t_soil_acoef  (ics:ice, :, iblk) = in%t_soil_acoef  (ics:ice, :, iblk)
!!$    out%t_soil_bcoef  (ics:ice, :, iblk) = in%t_soil_bcoef  (ics:ice, :, iblk)
!!$
!!$  END SUBROUTINE copy_t_soil_energy_state
!!$
!!$  SUBROUTINE add_t_soil_energy_state(in1, in2, out, iblk, ics, ice)
!!$
!!$    TYPE(t_soil_energy_state), INTENT(in)    :: in1, in2
!!$    TYPE(t_soil_energy_state), INTENT(inout) :: out
!!$    INTEGER,                   INTENT(in)    :: iblk, ics, ice
!!$
!!$    out%t_soil        (ics:ice, :, iblk) = in1%t_soil        (ics:ice, :, iblk) + in2%t_soil        (ics:ice, :, iblk)
!!$    out%t_soil_acoef  (ics:ice, :, iblk) = in1%t_soil_acoef  (ics:ice, :, iblk) + in2%t_soil_acoef  (ics:ice, :, iblk)
!!$    out%t_soil_bcoef  (ics:ice, :, iblk) = in1%t_soil_bcoef  (ics:ice, :, iblk) + in2%t_soil_bcoef  (ics:ice, :, iblk)
!!$
!!$  END SUBROUTINE add_t_soil_energy_state

#endif
END MODULE mo_soil_memory_class
