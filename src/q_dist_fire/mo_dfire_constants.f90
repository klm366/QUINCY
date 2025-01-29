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
MODULE mo_dfire_constants
#ifndef __NO_JSBACH__

  USE mo_kind,                 ONLY: wp
  USE mo_veg_constants,        ONLY: carbon_per_dryweight_leaf

  IMPLICIT NONE

  PUBLIC

  !>
  !> dfire parameters
  !>
  CHARACTER(len=*), PARAMETER, PRIVATE :: modname = 'mo_dfire_constants'

  !> soil top layers to consider in spitfire
  INTEGER, PARAMETER :: df_sf_spq_layer_top_sl = 1

  !< Top Soil layers to account for moisture in fuel
  INTEGER, PARAMETER :: df_sf_sl_dims_nr=2   !< soil moisture layer dimension from w_soil_sl(dim1, dim2)
  INTEGER, PARAMETeR :: df_sf_sl_indices(df_sf_sl_dims_nr) = [1,2]!< Select soil moisture layer to use

  CHARACTER(*), PARAMETER :: df_sf_vgrid_dead_fuel_pools_name="vg_dd_fl_pls" ! vgrid_dead_fuel_pools
  CHARACTER(*), PARAMETER :: df_sf_vgrid_live_fuel_pools_name="vg_lv_fl_pls" ! vgrid_live_fuel_pools

  REAL(wp), PARAMETER  :: df_sf_max_density_human_pop = 16._wp    !< Maximum density_human_population (eq. 3)

  INTEGER, PARAMETER :: sf_fuel_dead_fuel_nr = 7  !< (sb_pool) woody/leaf, (veg_pool grass) sap wood/labile/reserve/fruit/polymeric
  ENUM, BIND(C)
  !< Ensure how much to remove from each pool
    ENUMERATOR :: t_sf_fuel_dead_fuel_woody_ix = 1, &             !< from sb_pool_mt
                  t_sf_fuel_dead_fuel_grass_leaf_ix, &            !< from veg_pool_mt
                  t_sf_fuel_dead_fuel_grass_sap_ix, &             !< from veg_pool_mt
                  t_sf_fuel_dead_fuel_grass_labile_ix, &          !< from veg_pool_mt
                  t_sf_fuel_dead_fuel_grass_reserve_ix, &         !< from veg_pool_mt
                  t_sf_fuel_dead_fuel_grass_fruit_ix, &           !< from veg_pool_mt
                  t_sf_fuel_dead_fuel_polymeric_ix                !< from sb_pool_mt ( Optional)
  END ENUM

  !< Conversion factor from carbon to biomass
  REAL(wp), SAVE :: df_sf_carbon_to_biomass_factor

  ! Kelvin vs Celsius conversion
  REAL(wp), PARAMETER :: df_sf_units_c_to_k = 273.15_wp               !< conversion from degrees Celsius to kelvin
  REAL(wp), PARAMETER :: df_sf_units_k_to_c = -df_sf_units_c_to_k     !< conversion from degrees kelvin to Celsius

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

  !REAL(wp), PARAMETER :: df_sf_units_m2_to_cm2 = 10000._wp  !< units conversion from m-2 to cm-2

  !< "Live" fuel array indices from t_sf_fuel_class
  INTEGER, PARAMETER, PUBLIC :: sf_fuel_live_fuel_nr = 6 !< (veg_pool tree) leaf, sap wood, heart wood, labile, reserve, fruit

  ENUM, BIND(C)
    ENUMERATOR :: t_sf_fuel_live_fuel_leaf_ix = 1, &
                  t_sf_fuel_live_fuel_sap_wood_ix, &
                  t_sf_fuel_live_fuel_heart_wood_ix, &
                  t_sf_fuel_live_fuel_labile_ix, &
                  t_sf_fuel_live_fuel_reserve_ix, &
                  t_sf_fuel_live_fuel_fruit_ix
  END ENUM

  INTEGER, PARAMETER :: sf_fuel_classes_nr = 4 !< 1, 10, 100, 1000
  REAL(wp), PARAMETER :: sf_slow_fuel(sf_fuel_classes_nr) = [0.045_wp, 0.075_wp, 0.21_wp, 0.67_wp]
  REAL(wp), PARAMETER :: sf_fast_fuel(sf_fuel_classes_nr) = [1._wp, 0._wp, 0._wp, 0._wp]

  !< Methods to compute daily precipitation, as it is not done by the hosting model.
  ENUM, BIND(C)
    ENUMERATOR :: DFIRE_PRECIP_METHOD_DAYTIME=10    !< Previous day daily precipitation, based on sun movement (sunset vs sunrise)
    ENUMERATOR :: DFIRE_PRECIP_METHOD_MAVG=20       !< Pure moving average, so it is never approximatted.
    ENUMERATOR :: DFIRE_PRECIP_METHOD_APROX_MAVG=30 !< Approximation to moving average used by Quincy. Not recommended as it leaks water over time, as it never reaches 0.
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

  !< Decide at runtime which df_sf_live_frac_XXX is required (PFT dependent)
  !< Fraction allocation for live fuel (#nfuels, #live pools) [%]
  !REAL(wp), ALLOCATABLE, PUBLIC :: live_frac_allocation(:,:)

  ! todo: smarter way to do it? 
  LOGICAL, PARAMETER :: dfire_pft_mask(11) = [.TRUE., .TRUE., .TRUE., .TRUE., .TRUE., .TRUE., & ! tree pfts
                                              .TRUE., .TRUE., &                                 ! grass pfts
                                              .FALSE., .FALSE., &                               ! crop pfts
                                              .FALSE.]                                          ! Bar soil
CONTAINS

  ! ======================================================================================================= !
  !>
  !> Initialize constants used in the dfire process
  !>
  SUBROUTINE init_dfire_constants
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':init_dfire_constants'
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
    df_sf_dead_frac(:,t_sf_fuel_dead_fuel_grass_leaf_ix) =    df_sf_fuel_frac_fast
    df_sf_dead_frac(:,t_sf_fuel_dead_fuel_grass_sap_ix) =     df_sf_fuel_frac_fast
    df_sf_dead_frac(:,t_sf_fuel_dead_fuel_grass_labile_ix) =  df_sf_fuel_frac_fast
    df_sf_dead_frac(:,t_sf_fuel_dead_fuel_grass_reserve_ix) = df_sf_fuel_frac_fast
    df_sf_dead_frac(:,t_sf_fuel_dead_fuel_grass_fruit_ix) =   df_sf_fuel_frac_fast
  END SUBROUTINE init_dfire_constants

#endif
END MODULE mo_dfire_constants
