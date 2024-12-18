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
  USE mo_lnd_bgcm_store_class,  ONLY: SB_BGCM_LOSS_ID
  USE mo_jsb_tile_class,        ONLY: t_jsb_tile_abstract

  USE mo_dfire_constants,       ONLY: df_sf_sl_dims_nr, df_sf_sl_indices !< Soil related constants: top soil moisture layers 

  USE mo_dfire_constants,       ONLY: t_sf_fuel_dead_fuel_woody_ix

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: update_dfire_spitfire


  CHARACTER(len=*), PARAMETER :: modname = 'mo_dfire_spitfire'


  TYPE t_sf_bgcm_2L3D_p
    REAL(wp), POINTER :: p(:,:,:,:) !< dims: (part,elem,nc,layer)
  END TYPE t_sf_bgcm_2L3D_p

CONTAINS


  ! ======================================================================================================= !
  !>
  !> Compute spitfire for current time step
  !> It links Quincy and Spitfire
  !>
  SUBROUTINE update_dfire_spitfire(tile, options)
    USE mo_jsb_class,               ONLY: Get_model
    USE mo_jsb_task_class,          ONLY: t_jsb_task_options
    USE mo_jsb_model_class,         ONLY: t_jsb_model
    USE mo_jsb_process_class,       ONLY: VEG_, DFIRE_, SPQ_, DFIRE_FUEL_
    ! ----------------------------------------------------------------------------------------------------- !
    USE mo_jsb_grid,                ONLY: Get_vgrid, Get_grid
    USE mo_jsb_grid_class,          ONLY: t_jsb_vgrid, t_jsb_grid
    USE mo_jsb_math_constants,      ONLY: one_day

    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Use_memory(DFIRE_)
    dsl4jsb_Use_memory(DFIRE_FUEL_)
    dsl4jsb_Use_memory(SPQ_)
    dsl4jsb_Use_memory(VEG_)
    !dsl4jsb_Use_memory(SB_)
    !dsl4jsb_Use_memory(A2L_)
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Use_config(DFIRE_)
    ! ----------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_tile_abstract), INTENT(inout)     :: tile         !< one tile with data structure for one lct
    TYPE(t_jsb_task_options),   INTENT(in)        :: options      !< model options
    ! ----------------------------------------------------------------------------------------------------- !
    TYPE(t_jsb_grid),       POINTER       :: grid
    TYPE(t_jsb_model),      POINTER       :: model
    TYPE(t_lnd_bgcm_store), POINTER       :: bgcm_store !< the bgcm store of this tile
    INTEGER                               :: iblk, ics, ice, nc, ic
    INTEGER                               :: no_children
    REAL(wp)                              :: dtime      !< timestep
    CHARACTER(len=*), PARAMETER           :: routine = TRIM(modname)//':update_fdire_sf_fuel'

    REAL(wp), POINTER                     :: fract_fpc(:,:) ! in
    REAL(wp), POINTER                     :: grid_cell_area(:) ! [m2?]
    REAL(wp), POINTER                     :: tile_fract(:) ! [%]
    REAL(wp), ALLOCATABLE                 :: tile_area(:) ! [?]

    REAL(wp)                              :: ab_w_soil_sl(options%nc,df_sf_sl_dims_nr)    !< above ground soil water [m]
    REAL(wp)                              :: ab_w_soil_fc_sl(options%nc,df_sf_sl_dims_nr) !< above ground soil field capacity [m]
    REAL(wp)                              :: ab_moisture_sl(options%nc,df_sf_sl_dims_nr)  !< above ground relative soil moisture layers [-]
    REAL(wp)                              :: ab_moisture(options%nc)                      !< top relative soil moisture layer [-]

    INTEGER                               :: nr_of_children
    TYPE(t_sf_bgcm_2L3D_p), ALLOCATABLE   :: sb_loss_stores(:) !< Reference to each sb_loss pft pools (#npfts)
    ! ----------------------------------------------------------------------------------------------------- !
    !dsl4jsb_Def_mt2L2D :: veg_pool_mt
    dsl4jsb_Def_mt2L3D :: sb_loss_mt
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Def_config(DFIRE_)
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Def_memory(DFIRE_)
    dsl4jsb_Def_memory(DFIRE_FUEL_)
    dsl4jsb_Def_memory(VEG_)
    !dsl4jsb_Def_memory(SB_)
    dsl4jsb_Def_memory(SPQ_)
    !dsl4jsb_Def_memory(A2L_)
    ! ----------------------------------------------------------------------------------------------------- !
    ! Local vars here
    ! ----------------------------------------------------------------------------------------------------- !
    !dsl4jsb_Real2D_onChunk                :: e_n_ignitions    ! out, Total rate of ignitions events [ha-1 ts-1]
    dsl4jsb_Real3D_onChunk                :: w_soil_sl        ! soil water content [m]
    dsl4jsb_Real3D_onChunk                :: w_soil_fc_sl     ! soil water field capacity [m]

    ! DFIRE_FUEL_ 3D
    ! Note: t_air_max/t_air_min belongs to the previous days due to implementation constraints
    dsl4jsb_Real2D_onChunk                :: t_air_min         ! Min temperature of the day (previous day) [k]
    dsl4jsb_Real2D_onChunk                :: t_air_max         ! Max temperature of the day (previous day) [k]

    ! DFIRE_FUEL_ 3D
    dsl4jsb_Real3D_onChunk                :: dead_fuel_1hr     ! [g DM m-2]
    dsl4jsb_Real3D_onChunk                :: dead_fuel_10hr    ! 
    dsl4jsb_Real3D_onChunk                :: dead_fuel_100hr   ! 
    dsl4jsb_Real3D_onChunk                :: dead_fuel_1000hr  ! 
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
    model => Get_model(tile%owner_model_id)
    grid => Get_grid(model%grid_id)
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Get_config(DFIRE_)
    dsl4jsb_Get_memory(DFIRE_)
    dsl4jsb_Get_memory(DFIRE_FUEL_)
    dsl4jsb_Get_memory(VEG_)
    dsl4jsb_Get_memory(SPQ_)
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Get_var3D_onChunk(SPQ_, w_soil_sl)                ! in
    dsl4jsb_Get_var3D_onChunk(SPQ_, w_soil_fc_sl)             ! in

    dsl4jsb_Get_var2D_onChunk(DFIRE_, t_air_min)              ! in [k]
    dsl4jsb_Get_var2D_onChunk(DFIRE_, t_air_max)              ! in

    dsl4jsb_Get_var3D_onChunk(DFIRE_FUEL_, dead_fuel_1hr)     ! in [g DM m-2]
    dsl4jsb_Get_var3D_onChunk(DFIRE_FUEL_, dead_fuel_10hr)    ! in
    dsl4jsb_Get_var3D_onChunk(DFIRE_FUEL_, dead_fuel_100hr)   ! in
    dsl4jsb_Get_var3D_onChunk(DFIRE_FUEL_, dead_fuel_1000hr)  ! in
    ! ----------------------------------------------------------------------------------------------------- !
    bgcm_store => tile%bgcm_store
    !dsl4jsb_Get_mt2L2D(VEG_BGCM_POOL_ID, veg_pool_mt)
    ! ----------------------------------------------------------------------------------------------------- !
!    dsl4jsb_Get_mt2L3D(SB_BGCM_LOSS_ID, sb_loss_mt)

    no_children = tile%Get_no_of_children()

    !< Init variables
    fract_fpc => veg__mem%collect_var_real2d(ics, ice, iblk, no_children, veg__mem%fract_fpc)

    ! gather children's bgcm stores from specific flux pool
    sb_loss_stores = Get_childrens_bgcm_stores_2l_3d(SB_BGCM_LOSS_ID, ics, ice, iblk, tile)

    ab_moisture_sl = 0._wp
    ab_moisture = 0._wp

    write(*,*) "update_spitfire:: fract_fpc shape?", SHAPE(fract_fpc)

    ! ----------------------------------------------------------------------------------------------------- !
    ! Go science
    ! ----------------------------------------------------------------------------------------------------- !
    grid_cell_area => grid % area(ics:ice, iblk)
    tile_fract => tile%fract(ics:ice, iblk)

    tile_area = grid_cell_area * tile_fract

    !< Select the first 2 layers from the soil moisture and calculate soil moisture content as a fraction of field capacity
    ab_w_soil_sl(:,:) = w_soil_sl(:, [df_sf_sl_indices])
    ab_w_soil_fc_sl(:,:) = w_soil_fc_sl(:, [df_sf_sl_indices])
    WHERE (ab_w_soil_fc_sl > 0._wp) ! should it happen?
      ab_moisture_sl(:,:) = ab_w_soil_sl(:,:) / ab_w_soil_fc_sl(:,:)
    END WHERE

    !< moisture average of above ground soil layers
    DO ic=1, nc
      IF (SUM(ab_moisture_sl(ic,:)) > 0._wp) THEN
        ab_moisture(ic) = sum(ab_moisture_sl(ic,:)) / size(ab_moisture_sl, dim=2) !< 2 = soil layers index
      END IF
    ENDDO

    write(*,*) "(spitfire) I am burning..., iblk, dead_fuel_X=", iblk, dead_fuel_1hr(1,t_sf_fuel_dead_fuel_woody_ix), &
                                                            dead_fuel_10hr(1,t_sf_fuel_dead_fuel_woody_ix), &
                                                            dead_fuel_100hr(1,t_sf_fuel_dead_fuel_woody_ix), &
                                                            dead_fuel_1000hr(1,t_sf_fuel_dead_fuel_woody_ix)
    WRITE(*,*) "(spitfire) t_air_min, t_air_max = ", t_air_min(1), t_air_max(1)
  END SUBROUTINE update_dfire_spitfire


  ! ======================================================================================================= !
  !>
  !> Gather children's bgcm stores from specific flux pool
  !>
  FUNCTION Get_childrens_bgcm_stores_2l_3d(pool_id, ics, ice, iblk, tile) RESULT (bgcm_stores)
    INTEGER, INTENT(in) :: pool_id
    INTEGER, INTENT(in) :: ics, ice, iblk
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile         !< one tile with data structure for one lct
    TYPE(t_sf_bgcm_2L3D_p), ALLOCATABLE   :: bgcm_stores(:) ! output

    INTEGER :: i_tile, no_children
    CLASS(t_jsb_tile_abstract),  POINTER  :: current_tile

    no_children = tile%Get_no_of_children()
    ! gather children's bgcm stores from specific flux pool
    ALLOCATE(bgcm_stores(no_children)) !< move to spitfire_init, own function
    current_tile => tile%Get_first_child_tile()
    i_tile = 0
    DO WHILE (ASSOCIATED(current_tile))
      i_tile = i_tile + 1
      !write(*,*) "update_spitfire:: name tile=", current_tile%name
      bgcm_stores(i_tile)%p => current_tile%bgcm_store%Get_matrix_from_store_2l_3d(pool_id, ics, ice, iblk, current_tile%name)
      current_tile => current_tile%Get_next_sibling_tile()
    ENDDO
  END FUNCTION Get_childrens_bgcm_stores_2l_3d

#endif
END MODULE mo_dfire_spitfire
