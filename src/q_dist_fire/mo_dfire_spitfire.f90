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
!> @brief module to determine the spitfire fire fuel
!>
!> Main task to compute the spitfire
!>
!> -update_dfire_spitfire: integrate Quincy with the spitfire process:
!>    - convert units from Quincy into spitfire
!>    - gather lctlib parameters
!>    - select precalculated daily precipitation method
!>    - Soil moisture content
!>    - Fraction of the vegetation
!>    - Call calc_spitfire (science)
!>    - Dead feedback
!>    - Live feedback
!>    - Save data into process memory
!>
!> -calc_spitfire: spitfire subroutine which calls
!>    - [mask] method (file)
!>    - [-] Initialize t_sf_spitfire, used to output data
!>    - [-] Calc ignition events (mo_difre_sf_ignition.f90)
!>    - [-] Calc fuel load (mo_dfire_sf_fuel_class.f90)
!>    - [-] Calc Nesterov index (mo_dfire_sf_nesterov.f90)
!>    - [-] Calc fuel moisture (mo_dfire_sf_fmoist_class.f90)
!>    - [-] Calc Fire danger index (-)
!>    - [dead_fuel-all- > 0 & FDI>0] Calc Rate of spread (mo_dfire_sf_ros_class.f90)
!>    - [dead_fuel-all- > 0 & FDI>0] Calc Area Burned and its fraction (-)
!>    - [dead_fuel-1,10,100- > 0 & FDI>0] Calc dead fraction of consumption (mo_dfire_sf_frac_consump_class.f90)
!>    - [dead_fuel-1,10,100- > 0 & FDI>0]  Calc dead fuel consumption (mo_dfire_sf_fuel_consumption.f90)
!>    - [dead_fuel-1,10,100- > 0 & FDI>0]  Calc fire intensity (-)
!>    - [... & I_surface>50] Calc tree mortality (mo_dfire_sf_damage.f90)
!>    - [... & I_surface>50]  Calc live fuel consumption (mo_dfire_sf_fuel_consumption.f90)
!>    - [... & I_surface>50]  Calc emissions (mo_dfire_sf_emissions_class.f90)
!>    - [... & I_surface>50]  Derive fuel consumption for N,P,C13,C14,N15
!>    - [... & I_surface>50]  Calc residue from emissions
!>
!> mo_dfire_spitfire.f90 subroutines:
!>    -sf_calc_fire_intensity: calculate fire intensity
!>    -calc_area_burned_and_fract: calculate area burned and its fraction
!>
!> @author
!>  Albert Jornet Puig, University of Leipzig / MPI-BGC Jena
!>  Ana Bastos, University of Leipzig / MPI-BGC Jena
!>
MODULE mo_dfire_spitfire
#ifndef __NO_JSBACH__

  USE mo_kind,                  ONLY: wp
  USE mo_jsb_control,           ONLY: debug_on
  USE mo_exception,             ONLY: message, finish
  USE mo_lnd_bgcm_store,        ONLY: t_lnd_bgcm_store
  !USE mo_lnd_bgcm_store_class,  ONLY: VEG_BGCM_POOL_ID, VEG_BGCM_LOSS_ID, VEG_BGCM_LITTERFALL_ID

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: update_dfire_spitfire


  CHARACTER(len=*), PARAMETER :: modname = 'mo_dfire_spitfire'

CONTAINS


  ! ======================================================================================================= !
  !>
  !> Compute spitfire for current time step
  !> It links Quincy and Spitfire
  !>
  SUBROUTINE update_dfire_spitfire(tile, options)
    USE mo_jsb_class,               ONLY: Get_model
    USE mo_jsb_tile_class,          ONLY: t_jsb_tile_abstract
    USE mo_jsb_task_class,          ONLY: t_jsb_task_options
    USE mo_jsb_model_class,         ONLY: t_jsb_model
    USE mo_jsb_lctlib_class,        ONLY: t_lctlib_element
    USE mo_jsb_process_class,       ONLY: VEG_, DFIRE_ !, SB_, DFIRE_, SPQ_, A2L_
    ! ----------------------------------------------------------------------------------------------------- !
    USE mo_jsb_grid,                ONLY: Get_vgrid
    USE mo_jsb_grid_class,          ONLY: t_jsb_vgrid
    USE mo_jsb_math_constants,      ONLY: one_day

    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Use_memory(DFIRE_)
    !dsl4jsb_Use_memory(SPQ_)
    !dsl4jsb_Use_memory(VEG_)
    !dsl4jsb_Use_memory(SB_)
    !dsl4jsb_Use_memory(A2L_)
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Use_config(DFIRE_)
    ! ----------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_tile_abstract), INTENT(inout)     :: tile         !< one tile with data structure for one lct
    TYPE(t_jsb_task_options),   INTENT(in)        :: options      !< model options
    ! ----------------------------------------------------------------------------------------------------- !
    TYPE(t_jsb_model),      POINTER       :: model
    TYPE(t_lnd_bgcm_store), POINTER       :: bgcm_store !< the bgcm store of this tile
    TYPE(t_lctlib_element), POINTER       :: lctlib     !< land-cover-type library - parameter across pft's
    INTEGER                               :: iblk, ics, ice, nc
    REAL(wp)                              :: dtime      !< timestep
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':update_fdire_sf_fuel'
    ! ----------------------------------------------------------------------------------------------------- !
    !dsl4jsb_Def_mt2L2D :: veg_pool_mt
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Def_config(DFIRE_)
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Def_memory(DFIRE_)
    !dsl4jsb_Def_memory(VEG_)
    !dsl4jsb_Def_memory(SB_)
    !dsl4jsb_Def_memory(SPQ_)
    !dsl4jsb_Def_memory(A2L_)
    ! ----------------------------------------------------------------------------------------------------- !
    ! Local vars here
    ! ----------------------------------------------------------------------------------------------------- !
    !dsl4jsb_Real2D_onChunk                :: e_n_ignitions    ! out, Total rate of ignitions events [ha-1 ts-1]
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Get_memory(DFIRE_)            ! dfire__mem => dfire_memory
    ! ----------------------------------------------------------------------------------------------------- !
    iblk    = options%iblk
    ics     = options%ics
    ice     = options%ice
    nc      = options%nc
    dtime   = options%dtime
    ! ----------------------------------------------------------------------------------------------------- !
    IF (.NOT. tile%Is_process_calculated(DFIRE_)) RETURN
    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')
    ! ----------------------------------------------------------------------------------------------------- !
    model  => Get_model(tile%owner_model_id)
    lctlib => model%lctlib(tile%lcts(1)%lib_id) 
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Get_config(DFIRE_)
    dsl4jsb_Get_memory(DFIRE_)
    !dsl4jsb_Get_memory(VEG_)
    !dsl4jsb_Get_memory(SB_)
    !dsl4jsb_Get_memory(SPQ_)
    !dsl4jsb_Get_memory(A2L_)
    ! ----------------------------------------------------------------------------------------------------- !

    !dsl4jsb_Get_var2D_onChunk(DFIRE_, e_n_ignitions)                 ! out
    ! ----------------------------------------------------------------------------------------------------- !
    !bgcm_store => tile%bgcm_store
    !dsl4jsb_Get_mt2L2D(VEG_BGCM_POOL_ID, veg_pool_mt)
    ! ----------------------------------------------------------------------------------------------------- !

    !< Init variables

    ! ----------------------------------------------------------------------------------------------------- !
    ! Go science
    ! ----------------------------------------------------------------------------------------------------- !

  END SUBROUTINE update_dfire_spitfire


#endif
END MODULE mo_dfire_spitfire
