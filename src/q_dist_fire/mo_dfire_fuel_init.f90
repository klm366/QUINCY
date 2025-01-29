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
!> @brief Initialization of the the dfire fuel memory.
!>
!> @author
!>  Albert Jornet Puig, University of Leipzig / MPI-BGC Jena
!>  Ana Bastos, University of Leipzig / MPI-BGC Jena
!>
!> @par Revision History
!>
MODULE mo_dfire_fuel_init
#ifndef __NO_JSBACH__

  USE mo_kind,                    ONLY: wp
  USE mo_exception,               ONLY: message
  USE mo_jsb_control,             ONLY: debug_on

  USE mo_veg_constants,           ONLY: itree, igrass   !< lctlib%growthform

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: dfire_fuel_init

  CHARACTER(len=*), PARAMETER :: modname = 'mo_dfire_fuel_init'

CONTAINS

  ! ======================================================================================================= !
  !>
  !> Intialize disturbance fire process (after memory has been set up)
  !>
  SUBROUTINE dfire_fuel_init(tile)
    USE mo_jsb_class,               ONLY: get_model
    USE mo_jsb_lctlib_class,        ONLY: t_lctlib_element
    USE mo_jsb_process_class,       ONLY: DFIRE_FUEL_
    USE mo_jsb_tile_class,          ONLY: t_jsb_tile_abstract
    USE mo_jsb_model_class,         ONLY: t_jsb_model

    !dsl4jsb_Use_config(DFIRE_FUEL_)
    dsl4jsb_Use_memory(DFIRE_FUEL_)
    ! ----------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_tile_abstract), INTENT(inout)     :: tile         !< one tile with data structure for one lct
    ! ----------------------------------------------------------------------------------------------------- !
    TYPE(t_jsb_model),      POINTER              :: model
    TYPE(t_lctlib_element), POINTER              :: lctlib !< land-cover-type library - parameter across pfts
    INTEGER                                      :: iblk, ics, ice, nc
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':dfire_fuel_init'
    ! ----------------------------------------------------------------------------------------------------- !
    ! Declare process configuration and memory Pointers
    !dsl4jsb_Def_config(DFIRE_FUEL_)
    dsl4jsb_Def_memory(DFIRE_FUEL_)
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Real3D_onDomain      :: dead_fuel
    dsl4jsb_Real3D_onDomain      :: live_fuel

    dsl4jsb_Real2D_onDomain      :: forest_fract
    dsl4jsb_Real2D_onDomain      :: grass_fract

    dsl4jsb_Real2D_onDomain      :: fuel_bulk_density
    ! ----------------------------------------------------------------------------------------------------- !
    IF (.NOT. tile%Is_process_active(DFIRE_FUEL_)) RETURN
    IF (tile%lcts(1)%lib_id == 0) RETURN                !< run this init only if the present tile is a pft
    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')
    ! ----------------------------------------------------------------------------------------------------- !
    model  => get_model(tile%owner_model_id)
    lctlib => model%lctlib(tile%lcts(1)%lib_id)
    ! ----------------------------------------------------------------------------------------------------- !
    !dsl4jsb_Get_config(DFIRE_FUEL_)
    dsl4jsb_Get_memory(DFIRE_FUEL_)
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Get_var3D_onDomain(DFIRE_FUEL_, dead_fuel)
    dsl4jsb_Get_var3D_onDomain(DFIRE_FUEL_, live_fuel)

    dsl4jsb_Get_var2D_onDomain(DFIRE_FUEL_, forest_fract)
    dsl4jsb_Get_var2D_onDomain(DFIRE_FUEL_, grass_fract)

    dsl4jsb_Get_var2D_onDomain(DFIRE_FUEL_, fuel_bulk_density)
    ! ----------------------------------------------------------------------------------------------------- !
    WRITE(*,*) "Computing dfire_fuel_init..."
    dead_fuel = 0._wp
    live_fuel = 0._wp

    forest_fract = 0._wp
    grass_fract = 0._wp

    fuel_bulk_density = 0._wp
  END SUBROUTINE dfire_fuel_init

#endif
END MODULE mo_dfire_fuel_init
