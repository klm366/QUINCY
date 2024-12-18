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
MODULE mo_dfire_memory_class
#ifndef __NO_JSBACH__

  USE mo_kind,                   ONLY: wp
  USE mo_exception,              ONLY: message
  USE mo_jsb_class,              ONLY: Get_model
  USE mo_jsb_memory_class,       ONLY: t_jsb_memory
  USE mo_jsb_var_class,          ONLY: t_jsb_var_real2d, t_jsb_var_real3d

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: t_dfire_memory, max_no_of_vars

  INTEGER, PARAMETER :: max_no_of_vars = 80

  ! ======================================================================================================= !
  !>
  !> Type definition for dfire memory
  !>
  TYPE, EXTENDS(t_jsb_memory) :: t_dfire_memory
    !TYPE(t_jsb_var_real3d)               :: fuel_moisture      !< moisture in fuel [0 to 1]
    !TYPE(t_jsb_var_real2d)               :: nesterov_index     !< cumulative function of tmax and tdew [C2]

    TYPE(t_jsb_var_real2d)               :: t_air_max           !< Aproximation to daily max temp (from prev day) [k]
    TYPE(t_jsb_var_real2d)               :: t_air_min           !< Aproximation to daily min temp (from prev day) [k]

    TYPE(t_jsb_var_real2d)               :: pot_burn_area       !< Potential burnable area (forest + grass area only) [m2?]

    !< Precipitation
    !< Multiple methods used to aproximate daily precipitation. Not directly available into IQ.
    TYPE(t_jsb_var_real2d)               :: precip_qc_mavg      !< Precipitation moving average (quincy moving average) [kg m-2 s-1]
    TYPE(t_jsb_var_real2d)               :: precip_daytime      !< Aproximation to daily precipitation (previous daily data, aproximatted by sun movement) [kg m-2 s-1]
    TYPE(t_jsb_var_real2d)               :: precip_mavg         !< Precipitation moving average (Pure implementation of moving average) [kg m-2 s-1]
    TYPE(t_jsb_var_real2d)               :: precip_daily        !< Daily precipitation, it links to precip_daytime, precip_mavg or precip_qc_mavg [kg m-2 s-1]
    TYPE(t_jsb_var_real2d)               :: precip_daily_mmd    !< Daily precipitation, same as precip_daily but in spitfire expected units [mm day-1]
  CONTAINS
    PROCEDURE :: Init => Init_dfire_memory
  END TYPE t_dfire_memory

  CHARACTER(len=*), PARAMETER :: modname = 'mo_dfire_memory_class'

CONTAINS

  ! ======================================================================================================= !
  !>
  !> initialize memory (variables) for the process: dist_fire
  !>
  SUBROUTINE Init_dfire_memory(mem, prefix, suffix, lct_ids, model_id)
    USE mo_jsb_varlist,       ONLY: BASIC, MEDIUM, FULL
    USE mo_jsb_io,            ONLY: grib_bits, t_cf, t_grib1, t_grib2, tables
    USE mo_jsb_grid_class,    ONLY: t_jsb_grid, t_jsb_vgrid
    USE mo_jsb_grid,          ONLY: Get_grid, Get_vgrid
    USE mo_jsb_model_class,   ONLY: t_jsb_model
    ! ----------------------------------------------------------------------------------------------------- !
    CLASS(t_dfire_memory), INTENT(inout), TARGET  :: mem           !> dfire memory
    CHARACTER(len=*),     INTENT(in)            :: prefix          !> process name
    CHARACTER(len=*),     INTENT(in)            :: suffix          !> tile name
    INTEGER,              INTENT(in)            :: lct_ids(:)      !< Primary lct (1) and lcts of descendant tiles
    INTEGER,              INTENT(in)            :: model_id        !> model ID model\%id
    ! ----------------------------------------------------------------------------------------------------- !
    TYPE(t_jsb_model), POINTER :: model
    TYPE(t_jsb_grid),  POINTER :: hgrid                            !> Horizontal grid
    TYPE(t_jsb_vgrid), POINTER :: vgrid_surface                    !> Vertical grid
    INTEGER                    :: table
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':Init_dfire_memory'
    ! ----------------------------------------------------------------------------------------------------- !
    model         => Get_model(model_id)
    table         =  tables(1)
    hgrid         => Get_grid(model%grid_id)
    vgrid_surface         => Get_vgrid('surface')
    ! ----------------------------------------------------------------------------------------------------- !

    CALL mem%Add_var('t_air_max', mem%t_air_max, &
        & hgrid, vgrid_surface, &
        & t_cf('t_air_max', 'K', 'Max temperature of the previous day (daytime_counter)'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & loutput = .TRUE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

    CALL mem%Add_var('t_air_min', mem%t_air_min, &
        & hgrid, vgrid_surface, &
        & t_cf('t_air_min', 'K', 'Min temperature of the previous day (daytime_counter)'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & loutput = .TRUE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

    !< Precipitation related variables
    CALL mem%Add_var('precip_daily', mem%precip_daily, &
        & hgrid, vgrid_surface, &
        & t_cf('precip_daily', 'kg m-2 s-1', 'Daily precipitation (precip_mavg, precip_qc_mavg, precip_daytime'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & loutput = .TRUE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

    CALL mem%Add_var('precip_mavg', mem%precip_mavg, &
        & hgrid, vgrid_surface, &
        & t_cf('precip_mavg', 'kg m-2 s-1', 'Precipitation moving average using pure method'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & loutput = .TRUE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

    CALL mem%Add_var('precip_qc_mavg', mem%precip_qc_mavg, &
        & hgrid, vgrid_surface, &
        & t_cf('precip_qc_mavg', 'kg m-2 s-1', 'Precipitation moving average using quincy method'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & loutput = .TRUE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

    CALL mem%Add_var('precip_daily_mmd', mem%precip_daily_mmd, &
        & hgrid, vgrid_surface, &
        & t_cf('precip_daily_mmd', 'mm/day', 'Daily precipitation in Spitfire units'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & loutput = .TRUE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

    CALL mem%Add_var('precip_daytime', mem%precip_daytime, &
        & hgrid, vgrid_surface, &
        & t_cf('precip_daytime', 'kg m-2 s-1', 'Daily precipitation from previous day, related to sun movement'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & loutput = .TRUE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

    CALL mem%Add_var('pot_burn_area', mem%pot_burn_area, &
        & hgrid, vgrid_surface, &
        & t_cf('pot_burn_area', 'm2', 'Potential burnable area (forest and grass area)'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & loutput = .TRUE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

  END SUBROUTINE Init_dfire_memory

#endif
END MODULE mo_dfire_memory_class
