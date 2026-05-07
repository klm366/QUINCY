!>
!>  @par Copyright
!>    QUINCY terrestrial biosphere model with coupled carbon, nitrogen, and phosphorus cycles
!>    Copyright (C) 2019 Max Planck Institute for Biogeochemistry (MPI-BGC), Jena, Germany
!>    Contributors: Caldararu S., Engel J., Kern M., Pallandt M., Schnur R., Thum T., Yu L., Zaehle S.
!>
!>    This file is part of the QUINCY model.
!>    The QUINCY model is free software: you can redistribute it and/or modify
!>    it under the terms of the GNU General Public License version 3,
!>    as published by the Free Software Foundation.
!>
!>    The QUINCY model is distributed in the hope that it will be useful,
!>    but WITHOUT ANY WARRANTY; without even the implied warranty of
!>    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!>    GNU General Public License for more details.
!>
!>    You should have received a copy of the GNU General Public License
!>    along with the QUINCY model.  If not, see <https://www.gnu.org/licenses/>.
!>
!> @brief Contains the memory class for the dist_fire process.
!>
!> @author
!>  Albert Jornet Puig, University of Leipzig / MPI-BGC Jena
!>  Ana Bastos, University of Leipzig / MPI-BGC Jena
!>
!>
MODULE mo_q_dfire_fuel_memory_class
#ifndef __NO_JSBACH__

  USE mo_kind,                   ONLY: wp
  USE mo_exception,              ONLY: message
  USE mo_jsb_class,              ONLY: Get_model
  USE mo_jsb_memory_class,       ONLY: t_jsb_memory
  USE mo_jsb_var_class,          ONLY: t_jsb_var_real2d, t_jsb_var_real3d

  USE mo_q_dfire_constants,      ONLY: df_sf_vgrid_dead_fuel_pools_name, df_sf_vgrid_live_fuel_pools_name, &
                                       df_sf_vgrid_elements_name, df_sf_vgrid_species_name

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: t_q_dfire_fuel_memory, max_no_of_vars

  INTEGER, PARAMETER :: max_no_of_vars = 50

  ! ======================================================================================================= !
  !>
  !> Type definition for dfire memory
  !>
  TYPE, EXTENDS(t_jsb_memory) :: t_q_dfire_fuel_memory

    TYPE(t_jsb_var_real2d)  :: forest_fract        !< forest fraction [%]
    TYPE(t_jsb_var_real2d)  :: grass_fract         !< grass fraction [%]

    ! lctlib, weighted by tile fraction
    TYPE(t_jsb_var_real2d)  :: fuel_bulk_density          !< Fuel bulk density [kg m-3]
    TYPE(t_jsb_var_real2d)  :: moisture_extinction        !< Moisture of extinction threshold [0 to 1]
    TYPE(t_jsb_var_real2d)  :: fwd_wind_factor            !< Forward wind factor [0 to 1]

    TYPE(t_jsb_var_real3d)  :: dead_fuel_1hr_dfp          !< 1hr dead fuel (#nc, #ndeadpools) [g DM m-2]
    TYPE(t_jsb_var_real3d)  :: dead_fuel_10hr_dfp         !< 10hr dead fuel [g DM m-2], _dfp = dead fuel pools
    TYPE(t_jsb_var_real3d)  :: dead_fuel_100hr_dfp        !< 100hr dead fuel [g DM m-2], _dfp = dead fuel pools
    TYPE(t_jsb_var_real3d)  :: dead_fuel_1000hr_dfp       !< 1000hr dead fuel [g DM m-2], _dfp = dead fuel pools

    TYPE(t_jsb_var_real3d)  :: live_grass_fuel_1hr_lfp    !< 1hr live grass fuel (#nc, #nlivepools) [g DM m-2], _lfp = live fuel pools
    !< cured liveg grass is considered dead fuel
    TYPE(t_jsb_var_real3d)  :: cured_live_grass_fuel_1hr_dfp !< 1hr cured dead fuel (#nc, #nlivepools) [g DM m-2], _lfp = dead fuel pools

    TYPE(t_jsb_var_real3d)  :: live_tree_fuel_1hr_lfp     !< 1hr dead fuel (#nc, #nlivepools) [g DM m-2], _lfp = live fuel pools
    TYPE(t_jsb_var_real3d)  :: live_tree_fuel_10hr_lfp    !< 10hr dead fuel [g DM m-2], _lfp = live fuel pools
    TYPE(t_jsb_var_real3d)  :: live_tree_fuel_100hr_lfp   !< 100hr dead fuel [g DM m-2], _lfp = live fuel pools
    TYPE(t_jsb_var_real3d)  :: live_tree_fuel_1000hr_lfp  !< 1000hr dead fuel [g DM m-2], _lfp = live fuel pools

    TYPE(t_jsb_var_real2d)  :: pm_fire                    !< Tree probability of mortality of tree [- ts-1]
    TYPE(t_jsb_var_real2d)  :: pm_fire_output             !< Tree probability of mortality of tree [- s-1]
    TYPE(t_jsb_var_real2d)  :: pm_prop_ck                 !< Prob mortality: fraction of the crown scorch [-]
    TYPE(t_jsb_var_real2d)  :: pm_tau                     !< Probability of mortality due to cambial damage [-]
    TYPE(t_jsb_var_real2d)  :: tau_l                      !< Residence time of the fire [min]
    TYPE(t_jsb_var_real2d)  :: tau_c                      !< Critical time for cambial damage [min]
    TYPE(t_jsb_var_real2d)  :: bark_thickness             !< Bark thickness [cm]
    TYPE(t_jsb_var_real2d)  :: DBH_diag                   !< Diameter at breast height [cm]
    TYPE(t_jsb_var_real2d)  :: par1_diag                  !< Parameters for woody PFTs used in bark thickness calculation [-]
    TYPE(t_jsb_var_real2d)  :: par2_diag                  !< Parameters for woody PFTs used in bark thickness calculation [-]

    TYPE(t_jsb_var_real2d)  :: ck                         !< Fraction of the crown scorch [-]
    TYPE(t_jsb_var_real2d)  :: r_ck                       !< Resistance factor against crown damage [-]
    TYPE(t_jsb_var_real2d)  :: ck_p                       !< Parameter for woody PFTs used in Pm (CK) equation [-]
    TYPE(t_jsb_var_real2d)  :: scorch_height              !< Scorch height [m]
    TYPE(t_jsb_var_real2d)  :: F_scorch_height            !< PFT-parameter in crown scorch equation [-]
    TYPE(t_jsb_var_real2d)  :: crown_length               !< Crown length of woody PFT [m]
    TYPE(t_jsb_var_real2d)  :: tree_height                !< Tree height [m]

    TYPE(t_jsb_var_real2d)  :: lvg_moisture               !< Livegrass moisture content [0 to 2.5]%

    TYPE(t_jsb_var_real3d)  :: fuel_consumption_els       !< Fuel consumption per elements (#nc,#elements) [umol m-2 s-1]
    TYPE(t_jsb_var_real3d)  :: emissions_spc              !< Fire emissions per species (nc, #species) [g X m-2 s-1]
    TYPE(t_jsb_var_real3d)  :: emissions_els              !< Fire emissions per elements (nc, #elements) [umolX m-2 s-1]
    TYPE(t_jsb_var_real3d)  :: residue_els                !< Fire residue (nc, #elements) [umolX m-2 s-1]
    TYPE(t_jsb_var_real3d)  :: overflow_els               !< Emissions mass imbalance (nc, #elements) [umolX m-2 s-1]

    !< flux diagnostics
    TYPE(t_jsb_var_real3d)  :: veg_fire_litterfall_total_els   !< Total veg fire litterfall [umol m-2 s-1]
    TYPE(t_jsb_var_real3d)  :: veg_fire_emissions_total_els    !< Total veg fire emissions [umol m-2 s-1]
    TYPE(t_jsb_var_real3d)  :: sb_fire_emissions_total_els     !< Total sb fire emissions [umol m-2 s-1]
    TYPE(t_jsb_var_real3d)  :: sb_fire_litter_ash_total_els    !< Total sb fire litter ash [umol m-2 s-1]
  CONTAINS
    PROCEDURE :: Init => Init_q_dfire_fuel_memory
  END TYPE t_q_dfire_fuel_memory

  CHARACTER(len=*), PARAMETER :: modname = 'mo_q_dfire_fuel_memory_class'

CONTAINS

  ! ======================================================================================================= !
  !>
  !> initialize memory (variables) for the process: dist_fire
  !>
  SUBROUTINE Init_q_dfire_fuel_memory(mem, prefix, suffix, lct_ids, lib_id, model_id)
    USE mo_jsb_varlist,       ONLY: BASIC, MEDIUM, FULL
    USE mo_jsb_io,            ONLY: grib_bits, t_cf, t_grib1, t_grib2, tables
    USE mo_jsb_grid_class,    ONLY: t_jsb_grid, t_jsb_vgrid
    USE mo_jsb_grid,          ONLY: Get_grid, Get_vgrid
    USE mo_jsb_model_class,   ONLY: t_jsb_model

    ! ----------------------------------------------------------------------------------------------------- !
    CLASS(t_q_dfire_fuel_memory), INTENT(inout), TARGET  :: mem    !> dfire memory
    CHARACTER(len=*),     INTENT(in)            :: prefix          !> process name
    CHARACTER(len=*),     INTENT(in)            :: suffix          !> tile name
    INTEGER,              INTENT(in)            :: lct_ids(:)      !< Primary lct (1) and lcts of descendant tiles
    INTEGER,              INTENT(in)            :: lib_id          !< id of primary lct in lctlib
    INTEGER,              INTENT(in)            :: model_id        !> model ID model\%id
    ! ----------------------------------------------------------------------------------------------------- !
    TYPE(t_jsb_model), POINTER :: model
    TYPE(t_jsb_grid),  POINTER :: hgrid                            !> Horizontal grid
    TYPE(t_jsb_vgrid), POINTER :: vgrid_surface                    !> Vertical grid

    INTEGER                    :: table

    TYPE(t_jsb_vgrid), POINTER :: vgrid_dead_fuel_pools
    TYPE(t_jsb_vgrid), POINTER :: vgrid_live_fuel_pools
    TYPE(t_jsb_vgrid), POINTER :: vgrid_elements
    TYPE(t_jsb_vgrid), POINTER :: vgrid_species

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':Init_q_dfire_fuel_memory'
    ! ----------------------------------------------------------------------------------------------------- !
    model                 => Get_model(model_id)
    table                 =  tables(1)
    hgrid                 => Get_grid(model%grid_id)
    vgrid_surface         => Get_vgrid('surface')
    vgrid_dead_fuel_pools => Get_vgrid(df_sf_vgrid_dead_fuel_pools_name)
    vgrid_live_fuel_pools => Get_vgrid(df_sf_vgrid_live_fuel_pools_name)
    vgrid_elements        => Get_vgrid(df_sf_vgrid_elements_name)
    vgrid_species         => Get_vgrid(df_sf_vgrid_species_name)
    ! ----------------------------------------------------------------------------------------------------- !

    ! dead pools
    CALL mem%Add_var('dead_fuel_1hr_dfp', mem%dead_fuel_1hr_dfp, &
        & hgrid, vgrid_dead_fuel_pools, &
        & t_cf('dead_fuel_1hr_dfp', 'gDM m-2', '1 hour dead fuel pools'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = BASIC, &
        & loutput = .TRUE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

    CALL mem%Add_var('dead_fuel_10hr_dfp', mem%dead_fuel_10hr_dfp, &
        & hgrid, vgrid_dead_fuel_pools, &
        & t_cf('dead_fuel_10hr_dfp', 'gDM m-2', '10 hour dead fuel pools'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = BASIC, &
        & loutput = .TRUE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

    CALL mem%Add_var('dead_fuel_100hr_dfp', mem%dead_fuel_100hr_dfp, &
        & hgrid, vgrid_dead_fuel_pools, &
        & t_cf('dead_fuel_100hr_dfp', 'gDM m-2', '100 hour dead fuel pools'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = BASIC, &
        & loutput = .TRUE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

    CALL mem%Add_var('dead_fuel_1000hr_dfp', mem%dead_fuel_1000hr_dfp, &
        & hgrid, vgrid_dead_fuel_pools, &
        & t_cf('dead_fuel_1000hr_dfp', 'gDM m-2', '1000 hour dead fuel pools'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = BASIC, &
        & loutput = .TRUE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

    ! live grass fuel pools
    CALL mem%Add_var('live_grass_fuel_1hr_lfp', mem%live_grass_fuel_1hr_lfp, &
        & hgrid, vgrid_live_fuel_pools, &
        & t_cf('live_grass_fuel_1hr_lfp', 'gDM m-2', '1 hours live grass fuel pools'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = BASIC, &
        & loutput = .TRUE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

    ! live grass fuel pools
    CALL mem%Add_var('cured_live_grass_fuel_1hr_dfp', mem%cured_live_grass_fuel_1hr_dfp, &
        & hgrid, vgrid_dead_fuel_pools, &
        & t_cf('cured_live_grass_fuel_1hr_dfp', 'gDM m-2', '1 hours cured live grass fuel pools'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = BASIC, &
        & loutput = .TRUE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

    ! live tree fuel pools
    CALL mem%Add_var('live_tree_fuel_1hr_lfp', mem%live_tree_fuel_1hr_lfp, &
        & hgrid, vgrid_live_fuel_pools, &
        & t_cf('live_tree_fuel_1hr_lfp', 'gDM m-2', '1 hours live tree fuel pools'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = BASIC, &
        & loutput = .TRUE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)
    CALL mem%Add_var('live_tree_fuel_10hr_lfp', mem%live_tree_fuel_10hr_lfp, &
        & hgrid, vgrid_live_fuel_pools, &
        & t_cf('live_tree_fuel_10hr_lfp', 'gDM m-2', '10 hours live tree fuel pools'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = BASIC, &
        & loutput = .TRUE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)
    CALL mem%Add_var('live_tree_fuel_100hr_lfp', mem%live_tree_fuel_100hr_lfp, &
        & hgrid, vgrid_live_fuel_pools, &
        & t_cf('live_tree_fuel_100hr_lfp', 'gDM m-2', '100 hours live tree fuel pools'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = BASIC, &
        & loutput = .TRUE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)
    CALL mem%Add_var('live_tree_fuel_1000hr_lfp', mem%live_tree_fuel_1000hr_lfp, &
        & hgrid, vgrid_live_fuel_pools, &
        & t_cf('live_tree_fuel_1000hr_lfp', 'gDM m-2', '1000 hours live tree fuel pools'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = BASIC, &
        & loutput = .TRUE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

    CALL mem%Add_var('forest_fract', mem%forest_fract, &
        & hgrid, vgrid_surface, &
        & t_cf('forest_fract', '%', 'Forest fraction'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = BASIC, &
        & loutput = .TRUE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

    CALL mem%Add_var('grass_fract', mem%grass_fract, &
        & hgrid, vgrid_surface, &
        & t_cf('grass_fract', '%', 'Grass fraction'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = BASIC, &
        & loutput = .TRUE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

    CALL mem%Add_var('fuel_bulk_density', mem%fuel_bulk_density, &
        & hgrid, vgrid_surface, &
        & t_cf('fuel_bulk_density', 'kg m-3', 'Fuel bulk density'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = BASIC, &
        & loutput = .TRUE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

    CALL mem%Add_var('moisture_extinction', mem%moisture_extinction, &
        & hgrid, vgrid_surface, &
        & t_cf('moisture_extinction', '0 to 1', 'Moisture of extinction threshold'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = BASIC, &
        & loutput = .TRUE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

    CALL mem%Add_var('fwd_wind_factor', mem%fwd_wind_factor, &
        & hgrid, vgrid_surface, &
        & t_cf('fwd_wind_factor', '0 to 1', 'Forward wind factor'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = BASIC, &
        & loutput = .TRUE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

    CALL mem%Add_var('pm_fire', mem%pm_fire, &
        & hgrid, vgrid_surface, &
        & t_cf('pm_fire', '- ts-1', 'Tree probability of mortality by fire'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = BASIC, &
        & loutput = .TRUE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

    CALL mem%Add_var('pm_fire_output', mem%pm_fire_output, &
        & hgrid, vgrid_surface, &
        & t_cf('pm_fire_output', '- s-1', 'Tree probability of mortality by fire'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = BASIC, &
        & loutput = .TRUE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

    CALL mem%Add_var('pm_prop_ck', mem%pm_prop_ck, &
        & hgrid, vgrid_surface, &
        & t_cf('pm_prop_ck', '-', 'Proportion of the crown scorch'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = BASIC, &
        & loutput = .TRUE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

    CALL mem%Add_var('pm_tau', mem%pm_tau, &
        & hgrid, vgrid_surface, &
        & t_cf('pm_tau', '%', 'Probability of post-fire mortality due to cambial damage'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = BASIC, &
        & loutput = .TRUE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

    CALL mem%Add_var('tau_l', mem%tau_l, &
        & hgrid, vgrid_surface, &
        & t_cf('tau_l', 'min', 'Residence time of the fire'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = BASIC, &
        & loutput = .TRUE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

    CALL mem%Add_var('tau_c', mem%tau_c, &
        & hgrid, vgrid_surface, &
        & t_cf('tau_c', 'min', 'Critical time for cambial damage'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = BASIC, &
        & loutput = .TRUE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

    CALL mem%Add_var('par2_diag', mem%par2_diag, &
        & hgrid, vgrid_surface, &
        & t_cf('par2_diag', '-', 'Parameter 2 for woody PFTs used in bark thickness calculation'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = BASIC, &
        & loutput = .TRUE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

    CALL mem%Add_var('par1_diag', mem%par1_diag, &
        & hgrid, vgrid_surface, &
        & t_cf('par1_diag', '-', 'Parameter 1 for woody PFTs used in bark thickness calculation'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = BASIC, &
        & loutput = .TRUE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

    CALL mem%Add_var('bark_thickness', mem%bark_thickness, &
        & hgrid, vgrid_surface, &
        & t_cf('bark_thickness', 'cm', 'Bark thickness'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = BASIC, &
        & loutput = .TRUE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

    CALL mem%Add_var('DBH_diag', mem%DBH_diag, &
        & hgrid, vgrid_surface, &
        & t_cf('DBH_diag', 'cm', 'Diameter at breast height'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = BASIC, &
        & loutput = .TRUE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

    CALL mem%Add_var('ck', mem%ck, &
        & hgrid, vgrid_surface, &
        & t_cf('ck', '-', 'Fraction of the crown scorch'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = BASIC, &
        & loutput = .TRUE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

    CALL mem%Add_var('r_ck', mem%r_ck, &
        & hgrid, vgrid_surface, &
        & t_cf('r_ck', '-', 'Resistance factor against crown damage'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = BASIC, &
        & loutput = .TRUE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

    CALL mem%Add_var('ck_p', mem%ck_p, &
        & hgrid, vgrid_surface, &
        & t_cf('ck_p', '-', 'Parameter for woody PFTs used in Pm (CK) equation'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = BASIC, &
        & loutput = .TRUE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

    CALL mem%Add_var('scorch_height', mem%scorch_height, &
        & hgrid, vgrid_surface, &
        & t_cf('scorch_height', '-', 'Scorch height'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = BASIC, &
        & loutput = .TRUE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

    CALL mem%Add_var('F_scorch_height', mem%F_scorch_height, &
        & hgrid, vgrid_surface, &
        & t_cf('F_scorch_height', '-', 'PFT-parameter in crown scorch equation'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = BASIC, &
        & loutput = .TRUE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

    CALL mem%Add_var('crown_length', mem%crown_length, &
        & hgrid, vgrid_surface, &
        & t_cf('crown_length', 'm', 'Crown length of woody PFT'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = BASIC, &
        & loutput = .TRUE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

    CALL mem%Add_var('tree_height', mem%tree_height, &
        & hgrid, vgrid_surface, &
        & t_cf('tree_height', 'm', 'Tree height'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = BASIC, &
        & loutput = .TRUE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

    CALL mem%Add_var('lvg_moisture', mem%lvg_moisture, &
        & hgrid, vgrid_surface, &
        & t_cf('lvg_moisture', '%', 'Livegrass moisture content'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = BASIC, &
        & loutput = .TRUE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)


    CALL mem%Add_var('fuel_consumption_els', mem%fuel_consumption_els, &
        & hgrid, vgrid_elements, &
        & t_cf('fuel_consumption_els', 'umol X m-2 s-1', 'Fuel consumption'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = BASIC, &
        & loutput = .TRUE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

    CALL mem%Add_var('emissions_spc', mem%emissions_spc, &
        & hgrid, vgrid_species, &
        & t_cf('emissions_spc', 'g X m-2 s-1', 'Fuel emissions'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = BASIC, &
        & loutput = .TRUE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

    CALL mem%Add_var('emissions_els', mem%emissions_els, &
        & hgrid, vgrid_elements, &
        & t_cf('emissions_els', 'umol X m-2 s-1', 'Fuel emissions'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = BASIC, &
        & loutput = .TRUE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

    CALL mem%Add_var('residue_els', mem%residue_els, &
        & hgrid, vgrid_elements, &
        & t_cf('residue_els', 'umol X m-2 s-1', 'Residue'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = BASIC, &
        & loutput = .TRUE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

    CALL mem%Add_var('overflow_els', mem%overflow_els, &
        & hgrid, vgrid_elements, &
        & t_cf('overflow_els', 'umol X m-2 s-1', 'Emissions overflow'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = BASIC, &
        & loutput = .TRUE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

    CALL mem%Add_var('veg_fire_litterfall_total_els', mem%veg_fire_litterfall_total_els, &
        & hgrid, vgrid_elements, &
        & t_cf('veg_fire_litterfall_total_els', 'umol X m-2 s-1', 'Total veg fire litterfall'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = BASIC, &
        & loutput = .TRUE., &
        & lrestart = .FALSE., &
        & initval_r = 0.0_wp)

    CALL mem%Add_var('veg_fire_emissions_total_els', mem%veg_fire_emissions_total_els, &
        & hgrid, vgrid_elements, &
        & t_cf('veg_fire_emissions_total_els', 'umol X m-2 s-1', 'Total veg fire emissions'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = BASIC, &
        & loutput = .TRUE., &
        & lrestart = .FALSE., &
        & initval_r = 0.0_wp)

    CALL mem%Add_var('sb_fire_emissions_total_els', mem%sb_fire_emissions_total_els, &
        & hgrid, vgrid_elements, &
        & t_cf('sb_fire_emissions_total_els', 'umol X m-2 s-1', 'Total sb fire emissions'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = BASIC, &
        & loutput = .TRUE., &
        & lrestart = .FALSE., &
        & initval_r = 0.0_wp)

    CALL mem%Add_var('sb_fire_litter_ash_total_els', mem%sb_fire_litter_ash_total_els, &
        & hgrid, vgrid_elements, &
        & t_cf('sb_fire_litter_ash_total_els', 'umol X m-2 s-1', 'Total sb fire litter ash'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = BASIC, &
        & loutput = .TRUE., &
        & lrestart = .FALSE., &
        & initval_r = 0.0_wp)
  END SUBROUTINE Init_q_dfire_fuel_memory

#endif
END MODULE mo_q_dfire_fuel_memory_class
