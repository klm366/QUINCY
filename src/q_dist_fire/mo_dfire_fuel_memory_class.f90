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
MODULE mo_dfire_fuel_memory_class
#ifndef __NO_JSBACH__

  USE mo_kind,                   ONLY: wp
  USE mo_exception,              ONLY: message
  USE mo_jsb_class,              ONLY: Get_model
  USE mo_jsb_memory_class,       ONLY: t_jsb_memory
  USE mo_jsb_var_class,          ONLY: t_jsb_var_real2d, t_jsb_var_real3d

  USE mo_dfire_constants,        ONLY: df_sf_vgrid_dead_fuel_pools_name, df_sf_vgrid_live_fuel_pools_name

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: t_dfire_fuel_memory, max_no_of_vars

  INTEGER, PARAMETER :: max_no_of_vars = 20

  ! ======================================================================================================= !
  !>
  !> Type definition for dfire memory
  !>
  TYPE, EXTENDS(t_jsb_memory) :: t_dfire_fuel_memory
    TYPE(t_jsb_var_real3d)  :: dead_fuel_1hr       !< 1hr dead fuel (#nc, #ndeadpools) [gDM m-2]
    TYPE(t_jsb_var_real3d)  :: dead_fuel_10hr      !< 10hr dead fuel (#nc, #ndeadpools) [gDM m-2]
    TYPE(t_jsb_var_real3d)  :: dead_fuel_100hr     !< 100hr dead fuel (#nc, #ndeadpools) [gDM m-2]
    TYPE(t_jsb_var_real3d)  :: dead_fuel_1000hr    !< 1000hr dead fuel (#nc, #ndeadpools) [gDM m-2]

    TYPE(t_jsb_var_real3d)  :: live_fuel_1hr       !< 1hr live fuel (#nc, #nlivepools) [gDM m-2]
    TYPE(t_jsb_var_real3d)  :: live_fuel_10hr      !< 10hr live fuel (#nc, #nlivepools) [gDM m-2]
    TYPE(t_jsb_var_real3d)  :: live_fuel_100hr     !< 100hr live fuel (#nc, #nlivepools) [gDM m-2]
    TYPE(t_jsb_var_real3d)  :: live_fuel_1000hr    !< 1000hr live fuel (#nc, #nlivepools) [gDM m-2]

    TYPE(t_jsb_var_real2d)  :: forest_fract        !< forest fraction [%] 
    TYPE(t_jsb_var_real2d)  :: grass_fract         !< grass fraction [%]

    TYPE(t_jsb_var_real2d)  :: fuel_bulk_density   !< Fuel bulk density [kg m-3]
  CONTAINS
    PROCEDURE :: Init => Init_dfire_fuel_memory
  END TYPE t_dfire_fuel_memory

  CHARACTER(len=*), PARAMETER :: modname = 'mo_dfire_fuel_memory_class'

CONTAINS

  ! ======================================================================================================= !
  !>
  !> initialize memory (variables) for the process: dist_fire
  !>
  SUBROUTINE Init_dfire_fuel_memory(mem, prefix, suffix, lct_ids, model_id)
    USE mo_jsb_varlist,       ONLY: BASIC, MEDIUM, FULL
    USE mo_jsb_io,            ONLY: grib_bits, t_cf, t_grib1, t_grib2, tables
    USE mo_jsb_grid_class,    ONLY: t_jsb_grid, t_jsb_vgrid
    USE mo_jsb_grid,          ONLY: Get_grid, Get_vgrid
    USE mo_jsb_model_class,   ONLY: t_jsb_model

    ! ----------------------------------------------------------------------------------------------------- !
    CLASS(t_dfire_fuel_memory), INTENT(inout), TARGET  :: mem      !> dfire memory
    CHARACTER(len=*),     INTENT(in)            :: prefix          !> process name
    CHARACTER(len=*),     INTENT(in)            :: suffix          !> tile name
    INTEGER,              INTENT(in)            :: lct_ids(:)      !< Primary lct (1) and lcts of descendant tiles
    INTEGER,              INTENT(in)            :: model_id        !> model ID model\%id
    ! ----------------------------------------------------------------------------------------------------- !
    TYPE(t_jsb_model), POINTER :: model
    TYPE(t_jsb_grid),  POINTER :: hgrid                            !> Horizontal grid
    TYPE(t_jsb_vgrid), POINTER :: vgrid_surface                    !> Vertical grid

    INTEGER                    :: table

    TYPE(t_jsb_vgrid), POINTER :: vgrid_dead_fuel_pools
    TYPE(t_jsb_vgrid), POINTER :: vgrid_live_fuel_pools

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':Init_dfire_fuel_memory'
    ! ----------------------------------------------------------------------------------------------------- !
    model                 => Get_model(model_id)
    table                 =  tables(1)
    hgrid                 => Get_grid(model%grid_id)
    vgrid_surface         => Get_vgrid('surface')
    vgrid_dead_fuel_pools => Get_vgrid(df_sf_vgrid_dead_fuel_pools_name)
    vgrid_live_fuel_pools => Get_vgrid(df_sf_vgrid_live_fuel_pools_name)
    ! ----------------------------------------------------------------------------------------------------- !
    ! dead pools
    CALL mem%Add_var('dead_fuel_1hr', mem%dead_fuel_1hr, &
        & hgrid, vgrid_dead_fuel_pools, &
        & t_cf('dead_fuel_1hr', 'gDM m-2', '1hr dead fuel'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & loutput = .TRUE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

    CALL mem%Add_var('dead_fuel_10hr', mem%dead_fuel_10hr, &
        & hgrid, vgrid_dead_fuel_pools, &
        & t_cf('dead_fuel_10hr', 'gDM m-2', '10hr dead fuel'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & loutput = .TRUE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

    CALL mem%Add_var('dead_fuel_100hr', mem%dead_fuel_100hr, &
        & hgrid, vgrid_dead_fuel_pools, &
        & t_cf('dead_fuel_100hr', 'gDM m-2', '100hr dead fuel'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & loutput = .TRUE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

    CALL mem%Add_var('dead_fuel_1000hr', mem%dead_fuel_1000hr, &
        & hgrid, vgrid_dead_fuel_pools, &
        & t_cf('dead_fuel_1000hr', 'gDM m-2', '1000hr dead fuel'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & loutput = .TRUE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

    ! live pools
    CALL mem%Add_var('live_fuel_1hr', mem%live_fuel_1hr, &
        & hgrid, vgrid_live_fuel_pools, &
        & t_cf('live_fuel_1hr', 'gDM m-2', '1hr live fuel'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & loutput = .TRUE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

    CALL mem%Add_var('live_fuel_10hr', mem%live_fuel_10hr, &
        & hgrid, vgrid_live_fuel_pools, &
        & t_cf('live_fuel_10hr', 'gDM m-2', '10hr live fuel'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & loutput = .TRUE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

    CALL mem%Add_var('live_fuel_100hr', mem%live_fuel_100hr, &
        & hgrid, vgrid_live_fuel_pools, &
        & t_cf('live_fuel_100hr', 'gDM m-2', '100hr live fuel'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & loutput = .TRUE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

    CALL mem%Add_var('live_fuel_1000hr', mem%live_fuel_1000hr, &
        & hgrid, vgrid_live_fuel_pools, &
        & t_cf('live_fuel_1000hr', 'gDM m-2', '1000hr live fuel'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & loutput = .TRUE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

    CALL mem%Add_var('forest_fract', mem%forest_fract, &
        & hgrid, vgrid_surface, &
        & t_cf('forest_fract', '%', 'Forest fraction'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & loutput = .TRUE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

    CALL mem%Add_var('grass_fract', mem%grass_fract, &
        & hgrid, vgrid_surface, &
        & t_cf('grass_fract', '%', 'Grass fraction'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & loutput = .TRUE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

    CALL mem%Add_var('fuel_bulk_density', mem%fuel_bulk_density, &
        & hgrid, vgrid_surface, &
        & t_cf('fuel_bulk_density', 'kg m-3', 'Fuel bulk density'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & loutput = .TRUE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

  END SUBROUTINE Init_dfire_fuel_memory

#endif
END MODULE mo_dfire_fuel_memory_class
