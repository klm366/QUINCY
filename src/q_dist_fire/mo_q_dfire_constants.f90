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
!> @brief Contains constants for the dfire process
!>
!> @author
!>  Albert Jornet Puig, University of Leipzig / MPI-BGC Jena
!>  Ana Bastos, University of Leipzig / MPI-BGC Jena
!>
!>
MODULE mo_q_dfire_constants
#ifndef __NO_JSBACH__

  USE mo_kind,                 ONLY: wp
  USE mo_veg_constants,        ONLY: carbon_per_dryweight_leaf

  IMPLICIT NONE

  PUBLIC

  !>
  !> dfire parameters
  !>
  CHARACTER(len=*), PARAMETER, PRIVATE :: modname = 'mo_q_dfire_constants'

  !< sb pool to dead fuel ratios (polymoric is not included as it is fundamental)
  TYPE :: t_sf_litter_ratios
    REAL(wp) :: polymeric, dom_assoc, dom, soluable, residue, residue_assoc
  END TYPE t_sf_litter_ratios

  !> soil top layers to consider in spitfire
  INTEGER, PARAMETER :: df_sf_spq_layer_top_sl = 1

  !< Top Soil layers to account for moisture in fuel
  INTEGER, PARAMETER :: df_sf_sl_dims_nr=2   !< soil moisture layer dimension from w_soil_sl(dim1, dim2)
  INTEGER, PARAMETeR :: df_sf_sl_indices(df_sf_sl_dims_nr) = [1,2]!< Select soil moisture layer to use

  !< Fuel moisture
  INTEGER, PARAMETER :: df_sf_omegas_nr = 2      !< Total number of omegas = 2 (omega_o and omega_n_l)
  ENUM, BIND(C)
    ENUMERATOR ::  &
        df_sf_omega_o_idx = 1, &    !< Array index for omega_o (fuel moisture for 1, 10, 100 hours + grass)
        df_sf_omega_nl_idx  !< Array index for omega_n_l (fuel moisture for 1 hours + grass)
  END ENUM

  CHARACTER(*), PARAMETER :: df_sf_vgrid_dead_fuel_pools_name="vg_dd_fl_pls" ! vgrid_dead_fuel_pools
  CHARACTER(*), PARAMETER :: df_sf_vgrid_live_fuel_pools_name="vg_lv_fl_pls" ! vgrid_live_fuel_pools
  CHARACTER(*), PARAMETER :: df_sf_vgrid_fuel_classes_name="vg_fl_cls" ! vgrid_fuel_classes
  CHARACTER(*), PARAMETER :: df_sf_vgrid_fuel_moisture_name = "vg_fuel_moist" ! vgrid_fuel_moisture
  CHARACTER(*), PARAMETER :: df_sf_vgrid_species_name = "vg_spec"
  CHARACTER(*), PARAMETER :: df_sf_vgrid_elements_name = "vg_elements"
  CHARACTER(*), PARAMETER :: df_sf_vgrid_months_in_year_name = "vg_months"

  REAL(wp), PARAMETER  :: df_sf_max_density_human_pop = 16._wp    !< Maximum density_human_population (eq. 3)

  !< "Dead" fuel array indices from t_sf_fuel_class
  INTEGER, PUBLIC, PARAMETER :: sf_fuel_dead_fuel_nr = 12  !< woody, leaf, sap wood, labile, reserve, fruit, polymeric, dom assoc, dom, soluable, res, res assoc

  ENUM, BIND(C)
  !< Ensure how much to remove from each pool
    ENUMERATOR :: t_sf_fuel_dead_fuel_woody_ix = 1, &     !< from sb_pool_mt
                  t_sf_fuel_dead_fuel_grass_leaf_ix, &    !< from veg_pool_mt
                  t_sf_fuel_dead_fuel_grass_sap_ix, &     !< from veg_pool_mt
                  t_sf_fuel_dead_fuel_grass_labile_ix, &  !< from veg_pool_mt
                  t_sf_fuel_dead_fuel_grass_reserve_ix, & !< from veg_pool_mt
                  t_sf_fuel_dead_fuel_grass_fruit_ix, &   !< from veg_pool_mt
                  t_sf_fuel_dead_fuel_polymeric_ix, &     !< from sb_pool_mt ( Optional)
                  t_sf_fuel_dead_fuel_dom_assoc_ix, &     !< from sb_pool_mt ( Optional)
                  t_sf_fuel_dead_fuel_dom_ix, &           !< from sb_pool_mt ( Optional)
                  t_sf_fuel_dead_fuel_soluable_ix, &      !< from sb_pool_mt ( Optional)
                  t_sf_fuel_dead_fuel_res_ix, &           !< from sb_pool_mt ( Optional)
                  t_sf_fuel_dead_fuel_res_assoc_ix        !< from sb_pool_mt ( Optional)
  END ENUM

  INTEGER, PARAMETER :: ixs_sf_dead_livegrass_fuel_pools(5) = [t_sf_fuel_dead_fuel_grass_leaf_ix, &
                                                  t_sf_fuel_dead_fuel_grass_sap_ix, &
                                                  t_sf_fuel_dead_fuel_grass_labile_ix, &
                                                  t_sf_fuel_dead_fuel_grass_reserve_ix, &
                                                  t_sf_fuel_dead_fuel_grass_fruit_ix ]


  INTEGER, PARAMETER :: ixs_sf_dead_litter_fuel_pools(7) = [t_sf_fuel_dead_fuel_woody_ix, t_sf_fuel_dead_fuel_polymeric_ix, &
      t_sf_fuel_dead_fuel_dom_assoc_ix, t_sf_fuel_dead_fuel_dom_ix, t_sf_fuel_dead_fuel_soluable_ix, &
      t_sf_fuel_dead_fuel_res_ix, t_sf_fuel_dead_fuel_res_assoc_ix]

  INTEGER, PARAMETER :: df_sf_elements_nr = 6 !< C, C13, C14, N, N15, P

  !> Emissions
  INTEGER, PARAMETER :: df_sf_species_nr = 6 !< Thonicke, Table 1 for now we use only (CO2, CO, CH4, NH3, N2O, N2)
  ENUM, BIND(C)
    ENUMERATOR ::                   &
        df_sf_species_co2_idx = 1, & !< 44.01 g/mol, where C = 12.011, C_ratio=0.2729
        df_sf_species_co_idx, &      !< 28.01 g/mol, where C = 12.011, C_ratio=0.4288
        df_sf_species_ch4_idx, &     !< 16 g/mol, where C = 12.011, C_ratio=0.7507
        !df_sf_species_voc_idx, &    !< Note: discuss later
        !df_sf_species_tpm_idx, &
        !df_sf_species_nrox_idx, &
        df_sf_species_nh3_idx, &     !< 17.021 g/mol, where N = 14.0067, N_ratio=0.8229
        df_sf_species_n2o_idx, &     !< 44.013 g/mol, where N = 14.0067, N_ratio=0.3182
        df_sf_species_n2_idx         !< 28.0134 g/mol, where N = 14.0067, N_ratio=0.5
  END ENUM

  !< Conversion factor from carbon to biomass
  REAL(wp), SAVE :: df_sf_carbon_to_biomass_factor

  ! Kelvin vs Celsius conversion
  REAL(wp), PARAMETER :: df_sf_units_c_to_k = 273.15_wp               !< conversion from degrees Celsius to kelvin
  REAL(wp), PARAMETER :: df_sf_units_k_to_c = -df_sf_units_c_to_k     !< conversion from degrees kelvin to Celsius
  REAL(wp), PARAMETER :: df_sf_units_ft_to_cm = 30.48_wp              !< conversion from feet to cm

  ! Kg/m2/s to mm/day units conversion (precipitation)
  INTEGER, PARAMETER :: df_sk_units_kgm2s_to_mmday = 60*60*24         !<

  ! Grams to Kilograms conversions
  REAL(wp), PARAMETER :: df_sf_units_g_to_kg = 0.001                  !<

  ! Seconds to minutes conversions
  REAL(wp), PARAMETER :: df_sf_units_sec_to_min = 60.0_wp

  ! Seconds to minutes conversions
  REAL(wp), PARAMETER :: df_sf_units_m_to_cm = 100.0_wp

  ! km2 to ha conversion
  REAL(wp), PARAMETER :: df_sf_units_km2_to_ha = 100._wp
  REAL(wp), PARAMETER :: df_sf_units_m2_to_ha = 0.0001_wp

  REAL(wp), PARAMETER :: df_sf_calorific_heat_content = 18000._wp !< From Thonicke et al. (2010), Appendix A [kJ kg-1]

  REAL(wp), PARAMETER :: df_sf_units_m2_to_cm2 = 10000._wp  !< units conversion from m-2 to cm-2

  !< "Live" fuel array indices from t_sf_fuel_class, todo: rename to sf_live_fuel_pools_nr
  INTEGER, PARAMETER, PUBLIC :: sf_fuel_live_fuel_nr = 6 !< (veg_pool tree) leaf, sap wood, heart wood, labile, reserve, fruit

  ENUM, BIND(C)
    ENUMERATOR :: t_sf_fuel_live_fuel_leaf_ix = 1, &
                  t_sf_fuel_live_fuel_sap_wood_ix, &
                  t_sf_fuel_live_fuel_heart_wood_ix, &
                  t_sf_fuel_live_fuel_labile_ix, &
                  t_sf_fuel_live_fuel_reserve_ix, &
                  t_sf_fuel_live_fuel_fruit_ix
  END ENUM

  ! Fuel allocation from pools into dead/live fuel
  INTEGER, PARAMETER :: sf_fuel_classes_nr = 4 !< 1, 10, 100, 1000
  REAL(wp), PARAMETER :: sf_slow_fuel(sf_fuel_classes_nr) = [0.045_wp, 0.075_wp, 0.21_wp, 0.67_wp]
  REAL(wp), PARAMETER :: sf_fast_fuel(sf_fuel_classes_nr) = [1._wp, 0._wp, 0._wp, 0._wp]

  !< Methods to compute daily precipitation, as it is not done by the hosting model.
  ENUM, BIND(C)
    ENUMERATOR :: Q_DFIRE_PRECIP_METHOD_DAYTIME=10    !< Previous day daily precipitation, based on sun movement (sunset vs sunrise)
    ENUMERATOR :: Q_DFIRE_PRECIP_METHOD_MAVG=20       !< Pure moving average, so it is never approximatted.
    ENUMERATOR :: Q_DFIRE_PRECIP_METHOD_APROX_MAVG=30 !< Approximation to moving average used by Quincy. Not recommended as it leaks water over time, as it never reaches 0.
  END ENUM

  !< total number of fuel classes per type
  INTEGER, PARAMETER, PUBLIC :: df_sf_fc_full_nr = 4        !< Set for 1, 10, 100 and 1000 hours
  INTEGER, PARAMETER, PUBLIC :: df_sf_fc_spread_nr = 3      !< Set for 1, 10 and 100 hours (involved when fire spreads)
  INTEGER, PARAMETER, PUBLIC :: df_sf_fc_fast_nr = 1        !< Set for 1 hour
  ENUM, BIND(C)
    ENUMERATOR ::  &
        df_sf_fc_1hr_idx = 1, &    !< 1 hour fuel class array index
        df_sf_fc_10hr_idx, &       !< 10 hours fuel class array indices
        df_sf_fc_100hr_idx, &      !< 100 hours fuel class array indices
        df_sf_fc_1000hr_idx        !< 1000 hours fuel class array indices
  END ENUM

  !< Set fraction of fuel to allocate from sb_pool (dead fuel)
  REAL(wp) :: df_sf_dead_frac(df_sf_fc_full_nr, sf_fuel_dead_fuel_nr)
  !< Set fraction of fuel to allocate from veg_pool (live tree fuel)
  REAL(wp), TARGET :: df_sf_live_frac_tree(df_sf_fc_full_nr, sf_fuel_live_fuel_nr)
  !< Set fraction of fuel to allocate from veg_pool (live grass fuel)
  REAL(wp), TARGET :: df_sf_live_frac_grass(df_sf_fc_full_nr, sf_fuel_live_fuel_nr)

  !< Fuel class allocation for dead fuel (Pfeiffer et al. (2013) )
  REAL(wp), PARAMETER :: df_sf_fuel_frac_slow(df_sf_fc_full_nr) = [0.045_wp, 0.075_wp, 0.21_wp, 0.67_wp]

  !< Fuel class allocation for live fuel: 100% of live fuel is considered to be in the 1h class
  REAL(wp), PARAMETER, PUBLIC :: df_sf_fuel_frac_fast(df_sf_fc_full_nr) = [1._wp, 0._wp, 0._wp, 0._wp]

  !< Pfeiffer el al. (2013), Eq. A18 - Used to calculate fuel bulk density
  REAL(wp), PARAMETER, PUBLIC :: df_sf_fuel_class_factors(df_sf_fc_spread_nr) = [1._wp, 1.2_wp, 1.4_wp] !< For 1hr, 10hrs and 100 hrs

  !> Alpha fuel classes indices
  INTEGER, PARAMETER :: df_sf_alpha_nr = 3        !< 1hr, 10hr and 100hr
  ENUM, BIND(C)
    ENUMERATOR ::  &
            df_sf_alpha_1hr = 1, &  !< 1 hr alpha index
            df_sf_alpha_10hr, &     !< 10 hr alpha index
            df_sf_alpha_100hr       !< 100 hr alplha index
  END ENUM

  !< Inverse proportion to their surface-area-to-volume ratio [Celsius-2]
  !< Given in comments to Thonicke's et al. (2010) Eq.6
  REAL(wp), PARAMETER :: df_sf_alphas(df_sf_alpha_nr) = [0.001_wp, 0.0000542_wp, 0.0000149_wp]

  !< Nesterov threshold constants
  REAL(wp), PARAMETER :: df_sf_precip_trhd = 3.0_wp !< Precipitation threshold (Thonicke et al. (2010), below eq. 5)
  REAL(wp), PARAMETER :: df_sf_frost_trhd = 4.0_wp !< Frost threshold, from source code (Thonicke et al. (2010), https://github.com/PIK-LPJmL/LPJmL/blob/811b1691c26645968803fab65f9e463d64d5f40b/src/spitfire/update_nesterov.c#L25)
  REAL(wp), PARAMETER :: df_sf_oven_dry_part_dens = 513._wp !< Table A1, Thonicke et al. (2010) [kg m-3]

  !< Total number of sigmas in Thonicke et al. (2010) [cm-1]
  INTEGER, PARAMETER :: df_sf_sigma_tho_nr = 3
  REAL(wp), PARAMETER :: df_sigmas_tho(df_sf_sigma_tho_nr) = (/ 66._wp, 3.58_wp, 0.98_wp /)
  !< Oberhagemann et al. (2025)
  INTEGER, PARAMETER :: df_sf_sigma_dd_nr = 5 !< Total number of sigmas (1,10, 100, 1000 hours and live herbaceous)
  INTEGER, PARAMETER :: df_sf_sigma_dd_ros_nr = 4 !< Total number of sigmas (1,10, 100 hours and live herbaceous)
  INTEGER, PARAMETER :: df_sf_sigma_lv_nr = 2 !< Total number of sigmas (live grass, live woody)
  INTEGER, PARAMETER :: df_sf_ix_sigma_dd_1hr = 1
  INTEGER, PARAMETER :: df_sf_ix_sigma_dd_10hr = 2
  INTEGER, PARAMETER :: df_sf_ix_sigma_dd_100hr = 3
  INTEGER, PARAMETER :: df_sf_ix_sigma_dd_1000hr = 4
  INTEGER, PARAMETER :: df_sf_ix_sigma_dd_clvg = 5
  INTEGER, PARAMETER :: df_sf_ix_sigma_lv_lvg = 1
  INTEGER, PARAMETER :: df_sf_ix_sigma_lv_wdy = 1
  INTEGER, PARAMETER :: df_sf_ix_sigma_dd_ff(2) = [df_sf_ix_sigma_dd_1hr, df_sf_ix_sigma_dd_clvg] ! fine fuels
  INTEGER, PARAMETER :: df_sf_ix_sigma_dd_ros(4) = [df_sf_ix_sigma_dd_1hr, df_sf_ix_sigma_dd_10hr, df_sf_ix_sigma_dd_100hr, df_sf_ix_sigma_dd_clvg] ! fuel moistre involved fire rate of spread

  !< Surface area to volume ratio (Oberhagemann et al. (2025)) [cm-1]
  REAL(wp), PARAMETER :: df_sigmas_dead(df_sf_sigma_dd_ros_nr) = (/ 66._wp, 3.58_wp, 0.98_wp, 66._wp/) !< 1 hr, 10 hrs, 100 hrs, live herbaceous (to burn due to low moisture)
  !< Live woody is skipped from implementation in latest source code version (lpjml 5.7.10 - spitfire 1.9)
  !< However, live woody used value 1500 ft2/ft3 -> 42 m2/m3, loosely averaged from fuel models Andrews pl (2018)
  REAL(wp), PARAMETER :: df_sigmas_live(df_sf_sigma_lv_nr) = (/ 66._wp, 42._wp/) * df_sf_units_ft_to_cm !< (live herbaceous, live woody) [cm-1]

  REAL(wp), PARAMETER :: df_sf_mineral_fract = 0.055_wp !< Total mineral content of fuel (fraction, S_T in Thonicke et al. (2010), Table A1) [0 to 1]
  REAL(wp), PARAMETER :: df_sf_mineral_damp_coeff = 0.41739_wp !< Table A1, Thonicke et al. (2010) [0 to 1]
  REAL(wp), PARAMETER :: df_sf_fuel_heat_content = 18000._wp !<  Appendix A, Thonicke et al. (2010) [kJ kg-1]

  REAL(wp), PARAMETER :: df_sf_isurface_extinc_th = 50._wp !< Fires extinguishes below this threshold [kW m-1]
                                                           !< below eq. 15, Thonicke et al. (2010)

  REAL(wp), PARAMETER :: df_sf_fire_duration_thd = 30      !< Maximum fire duration allowed, capped to timestep length [min]
                                                           !< This is an adaptation from daily to 30 min TS as Thonicke et al. (2010) runs at a daily TS

  REAL(wp), PARAMETER :: df_sf_lign_cg = 0.20_wp           !< Lightning events, from cloud to ground flashes, Thonicke et al. (2010) section 2.2.1
  REAL(wp), PARAMETER :: df_sf_lign_eff = 0.04_wp          !< Lightning events, Efficiency to start fires from cloud to ground flashes, Thonicke et al. (2010) section 2.2.1

  ! todo: smarter way to do it?
  LOGICAL, SAVE :: burnable_tree_pfts(13)  ! TTTTTTFFFFFFF, default
  LOGICAL, SAVE :: burnable_grass_pfts(13) ! FFFFFFTTFFFFF, default
  LOGICAL, SAVE :: burnable_pfts(13) != [.TRUE., .TRUE., .TRUE., .TRUE., .TRUE., .TRUE., & ! tree pfts
                                        !      .TRUE., .TRUE., &                                 ! grass pfts
                                        !      .FALSE., .FALSE., &                               ! pasture pfts , 13 pfts only
                                        !      .FALSE., .FALSE., &                               ! crop pfts
                                        !      .FALSE.]                                          ! Bare soil

  TYPE(t_sf_litter_ratios), SAVE :: ratio_as_dead_fuel                                 !< sb_pool ratios to dead fuel pools [0 to 1]
CONTAINS

  ! ======================================================================================================= !
  !>
  !> Initialize constants used in the dfire process
  !>
  SUBROUTINE init_q_dfire_constants
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':init_q_dfire_constants'
    ! ----------------------------------------------------------------------------------------------------- !
    df_sf_carbon_to_biomass_factor    = 1._wp / carbon_per_dryweight_leaf

    !< Initialize live tree fuel fraction allocation
    df_sf_live_frac_tree = 0._wp !< default to 0
    df_sf_live_frac_tree(:,t_sf_fuel_live_fuel_leaf_ix) =       df_sf_fuel_frac_fast
    df_sf_live_frac_tree(:,t_sf_fuel_live_fuel_fruit_ix) =      df_sf_fuel_frac_fast
    df_sf_live_frac_tree(:,t_sf_fuel_live_fuel_sap_wood_ix) =   df_sf_fuel_frac_slow
    df_sf_live_frac_tree(:,t_sf_fuel_live_fuel_heart_wood_ix) = df_sf_fuel_frac_slow
    df_sf_live_frac_tree(:,t_sf_fuel_live_fuel_labile_ix) =     df_sf_fuel_frac_slow
    df_sf_live_frac_tree(:,t_sf_fuel_live_fuel_reserve_ix) =    df_sf_fuel_frac_slow

    !< Initialize live grass fuel fraction allocation
    df_sf_live_frac_grass = 0._wp !< default to 0
    df_sf_live_frac_grass(:,t_sf_fuel_live_fuel_leaf_ix) =      df_sf_fuel_frac_fast
    df_sf_live_frac_grass(:,t_sf_fuel_live_fuel_fruit_ix) =     df_sf_fuel_frac_fast
    df_sf_live_frac_grass(:,t_sf_fuel_live_fuel_sap_wood_ix) =  df_sf_fuel_frac_fast
    df_sf_live_frac_grass(:,t_sf_fuel_live_fuel_labile_ix) =    df_sf_fuel_frac_fast
    df_sf_live_frac_grass(:,t_sf_fuel_live_fuel_reserve_ix) =   df_sf_fuel_frac_fast

    !< Initialize dead fuel fraction allocation
    df_sf_dead_frac = 0._wp !< default to 0
    df_sf_dead_frac(:,t_sf_fuel_dead_fuel_woody_ix) =         df_sf_fuel_frac_slow
    df_sf_dead_frac(:,t_sf_fuel_dead_fuel_polymeric_ix) =     df_sf_fuel_frac_fast
    df_sf_dead_frac(:,t_sf_fuel_dead_fuel_dom_assoc_ix) =     df_sf_fuel_frac_fast
    df_sf_dead_frac(:,t_sf_fuel_dead_fuel_dom_ix) =           df_sf_fuel_frac_fast
    df_sf_dead_frac(:,t_sf_fuel_dead_fuel_soluable_ix) =      df_sf_fuel_frac_fast
    df_sf_dead_frac(:,t_sf_fuel_dead_fuel_res_ix) =           df_sf_fuel_frac_fast
    df_sf_dead_frac(:,t_sf_fuel_dead_fuel_res_assoc_ix) =     df_sf_fuel_frac_fast
    df_sf_dead_frac(:,t_sf_fuel_dead_fuel_grass_leaf_ix) =    df_sf_fuel_frac_fast
    df_sf_dead_frac(:,t_sf_fuel_dead_fuel_grass_sap_ix) =     df_sf_fuel_frac_fast
    df_sf_dead_frac(:,t_sf_fuel_dead_fuel_grass_labile_ix) =  df_sf_fuel_frac_fast
    df_sf_dead_frac(:,t_sf_fuel_dead_fuel_grass_reserve_ix) = df_sf_fuel_frac_fast
    df_sf_dead_frac(:,t_sf_fuel_dead_fuel_grass_fruit_ix) =   df_sf_fuel_frac_fast
  END SUBROUTINE init_q_dfire_constants

#endif
END MODULE mo_q_dfire_constants
