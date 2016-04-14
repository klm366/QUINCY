!>
!! @brief Contains the memory class for radiation.
!!
!! @author
!!  Rainer Schneck, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Rainer Schneck (2015-05-04)
!!

MODULE mo_rad_memory_class
#ifndef __NO_JSBACH__

  USE mo_kind, ONLY: wp
  USE mo_util, ONLY: One_of

  USE mo_jsb_memory_class, ONLY: t_jsb_memory
  USE mo_jsb_var_class,    ONLY: t_jsb_var_real2d, t_jsb_var_real3d
  USE mo_jsb_lct_class,    ONLY: BARE_TYPE, VEG_TYPE, LAKE_TYPE

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: t_rad_memory, max_no_of_vars

  INTEGER, PARAMETER :: max_no_of_vars = 40

  TYPE, EXTENDS(t_jsb_memory) :: t_rad_memory

    TYPE(t_jsb_var_real2d) :: &
      & sw_srf_net,           & !< net shortwave radiation at surface
      & swvis_srf_net,        & !< net shortwave radiation at surface in visible range
      & swnir_srf_net,        & !< net shortwave radiation at surface in NIR range
      & lw_srf_net,           & !< net longwave radiation at surface
      & rad_srf_net             !< net radiation at surface

    TYPE(t_jsb_var_real2d) ::                     &
      & alb_background,       & !< Background albedo (without snow or canopy)
                                !! Note, alb_background is ONLY for FUNCTION vegetation_albedo_simple_echam5 ?and update_par?.
      & alb_vis,              & !< Albedo in the VIS (visible) range
      & alb_nir,              & !< Albedo in the NIR (near-infrared) range
      & alb_vis_snow,         & !< Albedo of snow in the VIS range
      & alb_nir_snow            !< Albedo of snow in the NIR range

    !
    ! Additional variables for soil
    TYPE(t_jsb_var_real2d) :: &
      & alb_vis_soil,         & !< Albedo in VIS range of soil (may incl. litter layer and carbon within soil)
                                !! If use_albedo_soil_scheme = .FALSE. it is not calculated but taken from
                                !! bc_land_phys.nc file in the initialization
      & alb_nir_soil,         & ! Same for NIR
      ! R: Albedo of mineral soil is needed when Freya scheme is used. Not in the bc_lnd_phy.nc file yet!
      & alb_vis_mineralsoil,  & !< Albedo of soil without carbon or litter.
                                !! from bc_land_phys.nc file, necessary when use_albedo_soil_scheme = .TRUE.
                                !! for SOC-dependend albedo calculation
      & alb_nir_mineralsoil     ! from bc_land_phys.nc file, necessary when use_albedo_soil_scheme = .TRUE. for
                                ! SOC dependend albedo calculation

    !
    ! Additional variables for vegetation
    TYPE(t_jsb_var_real2d) :: &
      & alb_vis_can,          & !< Albedo of canopy (without snow cover) in VIS range. Is never calculated:
                                !! If use_alb_canopy=.TRUE. taken from bc_land_phys.nc file, else from lctlib-File (AlbedoCanopyVIS)
      & alb_nir_can,          & !< Albedo of canopy (without snow cover) in NIR range. Is never calculated:
                                !< If use_alb_canopy=.TRUE. taken from bc_land_phys.nc file, else from lctlib-File (AlbedoCanopyNIR)
      & sky_view_fract,       & !< Fraction of "visible" bare ground below vegetation without accounting for stems
      & sky_view_fract_stem,  & !< Fraction of soil below vegetation which is covered by stems but not by LAI
      & snow_age,             &
      !
      & par,                  & ! sw_par_down from atmosphere in [W/m^2]
      & par_down_mol,         & ! Downward PAR flux in mol(photons)/(m^2 s)
      & faPAR,                & ! fraction of PAR absorbed by the canopy
      & par_down_mol_nacc,    & ! incoming PAR accumulated (over output timestep) in mol(photons)/(m^2)
      & par_down_mol_tavg,    & ! incoming PAR (mean value over output timestep) [MOL PHOTONS/M^2 S]
      & soil_reflectivity_par,& ! Soil reflectivity in the PAR region
      & frac_par_direct,      & ! Fraction of direct PAR
      & apar_nacc

    TYPE(t_jsb_var_real3d) :: & ! alt: (nland, ncanopy, ntiles)
      & lai_cl,               & ! LAI per canopy layer; preliminarily, should later be put to another process...
      & faPAR_cl,             & ! Absorbed PAR per canopy layer [fraction of irradiance]
      & apar_per_lai_cl,      & ! Absorbed PAR per canopy layer. [(absorbed photons) / (m^2(leaf area) s)]
      & apar_per_lai_cl_nacc, & ! PAR absorbed by the canopy accumulated (over output timestep) in mol(photons)/m^2
      & apar_per_lai_cl_tavg    ! PAR absorbed by the canopy (mean value over output timestep)  in mol(photons)/(m^2 * s)

    !
    ! Additional variables for LAKE lct_type
    TYPE(t_jsb_var_real2d) :: &
      & alb_wtr,              & !< Albedo of open (lake) water
      & alb_ice                 !< Albedo of ice on lake

  CONTAINS
    PROCEDURE :: Init => Init_rad_memory
  END TYPE t_rad_memory

  CHARACTER(len=*), PARAMETER :: modname = 'mo_rad_memory_class'

CONTAINS

  SUBROUTINE Init_rad_memory(mem, prefix, suffix, lct_ids)

    USE mo_jsb_varlist,       ONLY: jsb_add_var
    USE mo_jsb_io,            ONLY: grib_bits, t_cf, t_grib1, t_grib2, &
                                    TSTEP_CONSTANT, tables
    USE mo_jsb_grid_class,    ONLY: t_jsb_grid, t_jsb_vgrid
    USE mo_jsb_grid,          ONLY: Get_grid, Get_vgrid

    CLASS(t_rad_memory), INTENT(inout) :: mem
    CHARACTER(len=*),    INTENT(in)    :: prefix
    CHARACTER(len=*),    INTENT(in)    :: suffix
    INTEGER,             INTENT(in)    :: lct_ids(:)

    TYPE(t_jsb_grid),  POINTER :: hgrid                        ! Horizontal grid
    TYPE(t_jsb_vgrid), POINTER :: surface, canopy_layer        ! Vertical grid
    INTEGER :: table

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':Init_rad_memory'

    table = tables(1)

    hgrid        => Get_grid(mem%grid_id)
    surface      => Get_vgrid('surface')
!!$    canopy_layer => Get_vgrid('canopy_layer')

    IF (.NOT. ASSOCIATED(mem%sw_srf_net%ptr)) THEN

      CALL jsb_add_var( mem, 'sw_srf_net', mem%sw_srf_net,                        &
        & hgrid, surface,                                                                                 &
        & t_cf('sw_srf_net', '', 'Net surface shortwave radiation'),                             &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                              &
        & prefix, suffix,                                                           &
        & initval_r=0.0_wp )

      CALL jsb_add_var( mem, 'swvis_srf_net', mem%swvis_srf_net,                  &
        & hgrid, surface,                                                                                 &
        & t_cf('swvis_srf_net', '', 'Net surface shortwave radiation in VIS'),                   &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                              &
        & prefix, suffix,                                                           &
        & initval_r=0.0_wp )

      CALL jsb_add_var( mem, 'swnir_srf_net', mem%swnir_srf_net,                  &
        & hgrid, surface,                                                                                 &
        & t_cf('swnir_srf_net', '', 'Net surface shortwave radiation in NIR'),                   &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                              &
        & prefix, suffix,                                                           &
        & initval_r=0.0_wp )

      CALL jsb_add_var( mem, 'lw_srf_net', mem%lw_srf_net,                        &
        & hgrid, surface,                                                                                 &
        & t_cf('lw_srf_net', '', 'Net surface longwave radiation'),                              &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                              &
        & prefix, suffix,                                                           &
        & initval_r=0.0_wp )

      CALL jsb_add_var( mem, 'rad_srf_net', mem%rad_srf_net,                      &
        & hgrid, surface,                                                                                 &
        & t_cf('rad_srf_net', '', 'Net surface radiation'),                                      &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                              &
        & prefix, suffix,                                                           &
        & initval_r=0.0_wp )

      CALL jsb_add_var( mem, 'alb_background', mem%alb_background,  &
        & hgrid, surface,                                                           &
        & t_cf('surface_albedo_background', '', ''),                       &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),        &
        & prefix, suffix,                                                           &
        & initval_r=0.0_wp, isteptype=TSTEP_CONSTANT )

      CALL jsb_add_var( mem, 'alb_vis', mem%alb_vis,                                &
        & hgrid, surface,                                                           &
        & t_cf('gridcell_surface_albedo_visible', '', ''),                 &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),        &
        & prefix, suffix,                                                           &
        & initval_r=0.0_wp, isteptype=TSTEP_CONSTANT )

      CALL jsb_add_var( mem, 'alb_nir', mem%alb_nir,                                &
        & hgrid, surface,                                                           &
        & t_cf('gridcell_surface_albedo_nir', '', ''),                     &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),        &
        & prefix, suffix,                                                           &
        & initval_r=0.0_wp )

      CALL jsb_add_var( mem, 'alb_vis_snow', mem%alb_vis_snow,      &
        & hgrid, surface,                                                           &
        & t_cf('albedo_vis_snow', '', ''),                                 &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),        &
        & prefix, suffix,                                                           &
        & loutput=.TRUE., lrestart=.TRUE.)

      CALL jsb_add_var( mem, 'alb_nir_snow', mem%alb_nir_snow,      &
        & hgrid, surface,                                                           &
        & t_cf('albedo_nir_snow', '', ''),                                 &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),        &
        & prefix, suffix,                                                           &
        & loutput=.TRUE., lrestart=.TRUE.)

    END IF

    ! Additional variables for soil
    IF ( (One_of(BARE_TYPE, lct_ids(:)) > 0 .OR. One_of(VEG_TYPE, lct_ids(:)) > 0) &
      & .AND. .NOT. ASSOCIATED(mem%alb_vis_soil%ptr)) THEN

      CALL jsb_add_var( mem, 'alb_vis_soil', mem%alb_vis_soil,      &
        & hgrid, surface,                                                           &
        & t_cf('soil_albedo_with_carbon_and_litter', '', ''),              &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),        &
        & prefix, suffix,                                                           &
        & loutput=.TRUE., lrestart=.TRUE., initval_r=0.0_wp )

      CALL jsb_add_var( mem, 'alb_nir_soil', mem%alb_nir_soil,      &
        & hgrid, surface,                                                           &
        & t_cf('soil_albedo_with_carbon_and_litter', '', ''),              &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),        &
        & prefix, suffix,                                                           &
        & loutput=.TRUE., lrestart=.TRUE., initval_r=0.0_wp )

      CALL jsb_add_var( mem, 'alb_vis_mineralsoil', mem%alb_vis_mineralsoil, &
        & hgrid, surface,                                                                    &
        & t_cf('alb_vis_mineralsoil', '', ''),                                      &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                 &
        & prefix, suffix,                                                           &
        & loutput=.TRUE., lrestart=.TRUE.)

      CALL jsb_add_var( mem, 'alb_nir_mineralsoil', mem%alb_nir_mineralsoil, &
        & hgrid, surface,                                                                    &
        & t_cf('alb_nir_mineralsoil', '', ''),                                      &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                 &
        & prefix, suffix,                                                           &
        & loutput=.TRUE., lrestart=.TRUE.)

    END IF

    ! Additional variables for vegetation
    IF (One_of(VEG_TYPE, lct_ids(:)) > 0 .AND. .NOT. ASSOCIATED(mem%alb_vis_can%ptr)) THEN

      CALL jsb_add_var( mem, 'alb_vis_can', mem%alb_vis_can,        &
      & hgrid, surface,                                                           &
      & t_cf('albedo_vis_canopy', '', ''),                               &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),        &
      & prefix, suffix,                                                           &
      & initval_r=0.0_wp )

      CALL jsb_add_var( mem, 'alb_nir_can', mem%alb_nir_can,        &
        & hgrid, surface,                                                           &
        & t_cf('albedo_nir_canpy', '', ''),                                &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),        &
        & prefix, suffix,                                                           &
        & initval_r=0.0_wp )

      CALL jsb_add_var( mem, 'snow_age', mem%snow_age,       &
        & hgrid, surface,                                                             &
        & t_cf('snow_age', '', ''),                                          &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),          &
        & prefix, suffix,                                                           &
        & loutput=.TRUE., lrestart=.TRUE., initval_r=0.0_wp )

      CALL jsb_add_var( mem, 'sky_view_fract', mem%sky_view_fract, &
        & hgrid, surface,                                                                             &
        & t_cf('sky_view_fract', '', ''),                                               &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                          &
        & prefix, suffix,                                                           &
        & loutput=.TRUE., lrestart=.FALSE.)

      CALL jsb_add_var( mem, 'sky_view_fract_stem', mem%sky_view_fract_stem,       &
        & hgrid, surface,                                                             &
        & t_cf('sky_view_fract_stem', '', ''),                                          &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),          &
        & prefix, suffix,                                                           &
        & loutput=.TRUE., lrestart=.FALSE.)

      CALL jsb_add_var( mem, 'par', mem%par, &
        & hgrid, surface,                                                                &
        & t_cf('Downward PAR flux in [W/m^2]', '', ''),                   &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),             &
        & prefix, suffix,                                                           &
        & lrestart=.FALSE., initval_r=0.0_wp )

      CALL jsb_add_var( mem, 'par_down_mol', mem%par_down_mol, &
        & hgrid, surface,                                                                &
        & t_cf('Downward PAR flux in mol (photons)', '', ''),                   &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),             &
        & prefix, suffix,                                                           &
        & lrestart=.FALSE., initval_r=0.0_wp ) ! R: initval =0

      CALL jsb_add_var( mem, 'faPAR', mem%faPAR,   &
        & hgrid, surface,                                                    &
        & t_cf('Fraction of PAR absorbed by the canopy', '', ''),   &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix,                                                           &
        & lrestart=.FALSE., initval_r=0.0_wp ) ! initval taken from JSBACH3

      CALL jsb_add_var( mem, 'par_down_mol_nacc', mem%par_down_mol_nacc, &
        & hgrid, surface,                                                                        &
        & t_cf('Incoming PAR accumulated (over output timestep)', '', ''),              &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                     &
        & prefix, suffix,                                                           &
        & lrestart=.FALSE., initval_r=0.0_wp )

      CALL jsb_add_var( mem, 'par_down_mol_tavg', mem%par_down_mol_tavg, &
        & hgrid, surface,                                                                          &
        & t_cf('Incoming PAR (mean value over output timestep)', '', ''),                 &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                       &
        & prefix, suffix,                                                           &
        & lrestart=.FALSE., initval_r=0.0_wp ) ! R: initval =0

      CALL jsb_add_var( mem, 'soil_reflectivity_par', mem%soil_reflectivity_par, &
        & hgrid, surface,                                                             &
        & t_cf('Soil reflectivity in the PAR region', '', ''),               &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),          &
        & prefix, suffix,                                                           &
        & loutput=.FALSE., lrestart=.FALSE., initval_r=0.0_wp )

      CALL jsb_add_var( mem, 'frac_par_direct', mem%frac_par_direct, &
        & hgrid, surface,                                                             &
        & t_cf('Fraction of direct PAR', '', ''),                            &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),          &
        & prefix, suffix,                                                           &
        & loutput=.TRUE., lrestart=.FALSE., initval_r=0.5_wp )

!!$    CALL jsb_add_var( mem, 'lai_cl', mem%lai_cl, &
!!$      & hgrid, canopy_layer,                                                          &
!!$      & t_cf('LAI for each canopy layer', '', ''),                           &
!!$      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),            &
!!$      & prefix, suffix,                                                           &
!!$      & lrestart=.FALSE., initval_r=0.0_wp )
!!$
!!$    CALL jsb_add_var( mem, 'faPAR_cl', mem%faPAR_cl, &
!!$      & hgrid, canopy_layer,                                                                   &
!!$      & t_cf('Absorbed PAR per canopy layer [fraction of irradiance]', '', ''),  &
!!$      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                &
!!$      & prefix, suffix,                                                           &
!!$      & lrestart=.FALSE., initval_r=0.0_wp )
!!$
!!$    CALL jsb_add_var( mem, 'apar_per_lai_cl', mem%apar_per_lai_cl, &
!!$      & hgrid, canopy_layer,                                                               &
!!$      & t_cf('Absorbed PAR per canopy layer (LAI area).', '', ''),                      &
!!$      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),            &
!!$      & prefix, suffix,                                                           &
!!$      & loutput=.FALSE., lrestart=.FALSE., initval_r=0.0_wp )
!!$
!!$    CALL jsb_add_var( mem, 'apar_per_lai_cl_nacc', mem%apar_per_lai_cl_nacc, &
!!$      & hgrid, canopy_layer,                                                                             &
!!$      & t_cf('PAR absorbed by the canopy (accumulated over output timestep)', '', ''),     &
!!$      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                          &
!!$      & prefix, suffix,                                                           &
!!$      & lrestart=.FALSE., initval_r=0.0_wp )
!!$
!!$    CALL jsb_add_var( mem, 'apar_per_lai_cl_tavg', mem%apar_per_lai_cl_tavg, &
!!$      & hgrid, canopy_layer,                                                                               &
!!$      & t_cf('PAR absorbed by the canopy (mean value over output timestep)', '', ''),        &
!!$      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                            &
!!$      & prefix, suffix,                                                           &
!!$      & lrestart=.FALSE., initval_r=0.0_wp )
!!$
!!$    CALL jsb_add_var( mem, 'apar_nacc', mem%apar_nacc, &
!!$      & hgrid, canopy_layer,                                                                             &
!!$      & t_cf('PAR absorbed by the canopy (accumulated over output timestep)', '', ''),     &
!!$      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                          &
!!$      & prefix, suffix,                                                           &
!!$      & lrestart=.FALSE., initval_r=0.0_wp )

    END IF

    ! Additional variables for lakes
    IF (One_of(LAKE_TYPE, lct_ids(:)) > 0 .AND. .NOT. ASSOCIATED(mem%alb_wtr%ptr)) THEN

      CALL jsb_add_var( mem, 'alb_wtr', mem%alb_wtr,                                &
        & hgrid, surface,                                                           &
        & t_cf('albedo_wtr', '', 'Albedo of open water'),                           &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),        &
        & prefix, suffix,                                                           &
        & loutput=.TRUE., lrestart=.TRUE., initval_r=0.0_wp )

      CALL jsb_add_var( mem, 'alb_ice', mem%alb_ice,                                &
        & hgrid, surface,                                                           &
        & t_cf('albedo_ice', '', 'Albedo of ice on water'),                         &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),        &
        & prefix, suffix,                                                           &
        & loutput=.TRUE., lrestart=.TRUE., initval_r=0.0_wp )

    END IF

  END SUBROUTINE Init_rad_memory

#endif
END MODULE mo_rad_memory_class
