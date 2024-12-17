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

  IMPLICIT NONE

  PUBLIC

  !>
  !> dfire parameters
  !>
  CHARACTER(len=*), PARAMETER, PRIVATE :: modname = 'mo_dfire_constants'

  !< Top Soil layers to account for moisture in fuel
  INTEGER, PARAMETER :: df_sf_sl_dims_nr=2   !< soil moisture layer dimension from w_soil_sl(dim1, dim2)
  INTEGER, PARAMETeR :: df_sf_sl_indices(df_sf_sl_dims_nr) = [1,2]!< Select soil moisture layer to use

  CHARACTER(*), PARAMETER :: df_sf_vgrid_dead_fuel_pools_name="vg_dd_fl_pls" ! vgrid_dead_fuel_pools
  CHARACTER(*), PARAMETER :: df_sf_vgrid_live_fuel_pools_name="vg_lv_fl_pls" ! vgrid_live_fuel_pools

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
CONTAINS

  ! ======================================================================================================= !
  !>
  !> Initialize constants used in the dfire process
  !>
  SUBROUTINE init_dfire_constants
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':init_dfire_constants'
    ! ----------------------------------------------------------------------------------------------------- !
  END SUBROUTINE init_dfire_constants

#endif
END MODULE mo_dfire_constants
