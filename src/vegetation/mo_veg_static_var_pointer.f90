!>
!!  @par Copyright
!!    QUINCY terrestrial biosphere model with coupled carbon, nitrogen, and phosphorus cycles
!!    Copyright (C) 2019 Max Planck Institute for Biogeochemistry (MPI-BGC), Jena, Germany 
!!    Contributors: Caldararu S., Engel J., Kern M., Pallandt M., Schnur R., Thum T., Yu L., Zaehle S.
!!
!!    This file is part of the QUINCY model.
!!    The QUINCY model is free software: you can redistribute it and/or modify
!!    it under the terms of the GNU General Public License version 3,
!!    as published by the Free Software Foundation.
!!
!!    The QUINCY model is distributed in the hope that it will be useful,
!!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!!    GNU General Public License for more details.
!!
!!    You should have received a copy of the GNU General Public License
!!    along with the QUINCY model.  If not, see <https://www.gnu.org/licenses/>.
!! 
!! @brief declare and init objects of variables used over one timestep
!!
!! @author
!!  Jan Engel, MPI-BGC Jena
!!
MODULE mo_veg_static_var_pointer
#ifndef __NO_JSBACH__

  USE mo_jsb_control,           ONLY: debug_on
  USE mo_exception,             ONLY: message
  USE mo_jsb_process_class,     ONLY: VEG_
  dsl4jsb_Use_memory(VEG_)
  dsl4jsb_Use_pools(VEG_) t_veg_tissues_pool, t_veg_elements_pool

  IMPLICIT NONE

  !>
  !> declare objects that are part of t_veg_main_pool
  !>
  TYPE(t_veg_tissues_pool),   POINTER, SAVE :: veg_pool
  TYPE(t_veg_tissues_pool),   POINTER, SAVE :: veg_growth    
  TYPE(t_veg_tissues_pool),   POINTER, SAVE :: veg_litterfall
  TYPE(t_veg_tissues_pool),   POINTER, SAVE :: veg_total_flux
  TYPE(t_veg_tissues_pool),   POINTER, SAVE :: veg_falloc
  TYPE(t_veg_elements_pool),  POINTER, SAVE :: veg_total_biomass
  TYPE(t_veg_elements_pool),  POINTER, SAVE :: veg_exudation    
  TYPE(t_veg_elements_pool),  POINTER, SAVE :: veg_establishment
  TYPE(t_veg_elements_pool),  POINTER, SAVE :: veg_reserve_use  


  PRIVATE
  PUBLIC :: init_veg_static_var_pointers
  PUBLIC :: veg_pool, veg_growth, veg_litterfall, veg_total_flux, veg_falloc
  PUBLIC :: veg_total_biomass, veg_exudation, veg_establishment, veg_reserve_use

  CHARACTER(len=*), PARAMETER :: modname = 'mo_veg_static_var_pointer'

CONTAINS

  ! ======================================================================================================= !
  !>
  !> init objects with pointers to the current chunk of variables used over one timestep
  !> 
  SUBROUTINE init_veg_static_var_pointers(tile, options)
    USE mo_jsb_tile_class,        ONLY: t_jsb_tile_abstract
    USE mo_jsb_task_class,        ONLY: t_jsb_task_options
    ! ----------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_tile_abstract), INTENT(inout)     :: tile         !< one tile with data structure for one lct
    TYPE(t_jsb_task_options),   INTENT(in)        :: options      !< model options
    ! ----------------------------------------------------------------------------------------------------- !
    INTEGER                                 :: iblk, ics, ice, nc
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':init_veg_static_var_pointers'
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Def_memory(VEG_)
    ! ----------------------------------------------------------------------------------------------------- !
    iblk    = options%iblk
    ics     = options%ics
    ice     = options%ice
    nc      = options%nc
    ! ----------------------------------------------------------------------------------------------------- !
    IF (.NOT. tile%Is_process_calculated(VEG_)) RETURN
    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Get_memory(VEG_)

    
    !>
    !> 1.0 define objects
    !>
    veg_pool            => veg__mem%veg_main_pool%veg_pool
    veg_growth          => veg__mem%veg_main_pool%veg_growth
    veg_litterfall      => veg__mem%veg_main_pool%veg_litterfall
    veg_total_flux      => veg__mem%veg_main_pool%veg_total_flux
    veg_falloc          => veg__mem%veg_main_pool%veg_falloc
    veg_total_biomass   => veg__mem%veg_main_pool%veg_total_biomass
    veg_exudation       => veg__mem%veg_main_pool%veg_exudation
    veg_establishment   => veg__mem%veg_main_pool%veg_establishment
    veg_reserve_use     => veg__mem%veg_main_pool%veg_reserve_use

    !>
    !> 2.0 init objects - set pointers to chunk
    !>
    CALL veg_pool%Associate_chunk_var_ptr(ics, ice, iblk)
    CALL veg_growth%Associate_chunk_var_ptr(ics, ice, iblk)
    CALL veg_litterfall%Associate_chunk_var_ptr(ics, ice, iblk)
    CALL veg_total_flux%Associate_chunk_var_ptr(ics, ice, iblk)
    CALL veg_falloc%Associate_chunk_var_ptr(ics, ice, iblk)
    CALL veg_total_biomass%Associate_chunk_var_ptr(ics, ice, iblk)
    CALL veg_exudation%Associate_chunk_var_ptr(ics, ice, iblk)
    CALL veg_establishment%Associate_chunk_var_ptr(ics, ice, iblk)
    CALL veg_reserve_use%Associate_chunk_var_ptr(ics, ice, iblk)
  END SUBROUTINE init_veg_static_var_pointers

#endif
END MODULE mo_veg_static_var_pointer
