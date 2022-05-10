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
!! @brief Contains the routines for the Task update_spq_physics
!!
!! @author
!!  Sönke Zaehle, MPI-BGC Jena 
!!
!! @par Revision History
!! First version by Sönke Zaehle (20160104)
!!
MODULE mo_spq_update_physics
#ifndef __NO_JSBACH__

#ifdef __QUINCY_STANDALONE__
#include "dsl_jsbach.inc"
#endif

  USE mo_kind,                  ONLY: wp
  USE mo_jsb_control,           ONLY: debug_on
  USE mo_exception,             ONLY: message

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: update_spq_physics

  CHARACTER(len=*), PARAMETER :: modname = 'mo_spq_update_physics'

CONTAINS


  !-----------------------------------------------------------------------------------------------------
  ! Main Task
  ! 
  !-----------------------------------------------------------------------------------------------------
  !> update soil moisture and theta for the plant only model
  !! 
  !! two cases considered: \n a) dynamic water uptake calculation, \n b) prescribed values from namelist
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE update_spq_physics(tile, options)

    USE mo_jsb_class,             ONLY: Get_model
    USE mo_jsb_tile_class,        ONLY: t_jsb_tile_abstract
    USE mo_jsb_task_class,        ONLY: t_jsb_task_options
    USE mo_jsb_model_class,       ONLY: t_jsb_model
    USE mo_jsb_lctlib_class,      ONLY: t_lctlib_element    
    USE mo_jsb_process_class,     ONLY: A2L_, SEB_, TURB_, SPQ_, HYDRO_, VEG_, HD_, RAD_, ASSIMI_, PHENO_, SB_
    USE mo_jsb_grid_class,        ONLY: t_jsb_vgrid
    USE mo_jsb_grid,              ONLY: Get_vgrid

    ! Use of process configurations (t_PROC_config)
    ! ..

    ! Use of process memories
    dsl4jsb_Use_memory(A2L_) 
    dsl4jsb_Use_memory(ASSIMI_)     !  USE mo_assimi_memory_class,   ONLY: t_assimi_memory
    dsl4jsb_Use_memory(RAD_)        !  USE mo_rad_memory_class,      ONLY: t_rad_memory
    dsl4jsb_Use_memory(SPQ_)        !  USE mo_spq_memory_class,      ONLY: t_spq_memory
    dsl4jsb_Use_memory(VEG_)        !  USE mo_veg_memory_class,      ONLY: t_veg_memory
    dsl4jsb_Use_memory(PHENO_)

    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    CLASS(t_jsb_tile_abstract), INTENT(inout)     :: tile         !< one tile with data structure for one lct
    TYPE(t_jsb_task_options),   INTENT(in)        :: options      !< model options
    ! ---------------------------
    ! 0.2 Local
    TYPE(t_jsb_model),      POINTER       :: model
    TYPE(t_lctlib_element), POINTER       :: lctlib         !< land-cover-type library - parameter across pft's
    TYPE(t_jsb_vgrid),      POINTER       :: vgrid_soil_sb  !< Vertical grid
    INTEGER                               :: nsoil_sb       !< number of soil layers as used/defined by the SB_ process
    REAL(wp), DIMENSION(options%nc)       :: rel_hum
    INTEGER                               :: isoil
    INTEGER                               :: iblk, ics, ice, nc
    REAL(wp)                              :: steplen
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':update_spq_physics'
    ! ---------------------------
    ! 0.3 Declare Memory
    ! Declare process configuration and memory Pointers
    dsl4jsb_Def_memory(A2L_) 
    dsl4jsb_Def_memory(ASSIMI_)                     !     TYPE(t_assimi_memory),   POINTER :: assimi__mem
    dsl4jsb_Def_memory(RAD_)                        !     TYPE(t_rad_memory),      POINTER :: rad__mem
    dsl4jsb_Def_memory(SPQ_)                       !     TYPE(t_spq_memory),     POINTER :: spq__mem
    dsl4jsb_Def_memory(VEG_)                        !     TYPE(t_veg_memory),      POINTER :: veg__mem
    dsl4jsb_Def_memory(PHENO_)
    ! Declare pointers to variables in memory
    ! A2L_
    dsl4jsb_Real2D_onChunk                :: t_air
    dsl4jsb_Real2D_onChunk                :: q_air
    dsl4jsb_Real2D_onChunk                :: press_srf
    dsl4jsb_Real2D_onChunk                :: wind_10m
    dsl4jsb_Real2D_onChunk                :: lw_srf_down
    dsl4jsb_Real2D_onChunk                :: rain
    dsl4jsb_Real2D_onChunk                :: snow
    dsl4jsb_Real2D_onChunk                :: t_acoef 
    dsl4jsb_Real2D_onChunk                :: t_bcoef 
    dsl4jsb_Real2D_onChunk                :: q_acoef
    dsl4jsb_Real2D_onChunk                :: q_bcoef
    dsl4jsb_Real2D_onChunk                :: ga
    dsl4jsb_Real2D_onChunk                :: pch 
    dsl4jsb_Real2D_onChunk                :: drag_srf
    ! SPQ_ 2D
    dsl4jsb_Real2D_onChunk                :: fact_q_air 
    dsl4jsb_Real2D_onChunk                :: fact_qsat_srf
    dsl4jsb_Real2D_onChunk                :: t_srf_old
    dsl4jsb_Real2D_onChunk                :: zril_old
    dsl4jsb_Real2D_onChunk                :: root_depth
    dsl4jsb_Real2D_onChunk                :: w_soil_root_pwp
    dsl4jsb_Real2D_onChunk                :: w_soil_root_fc
    dsl4jsb_Real2D_onChunk                :: t_soil_root
    dsl4jsb_Real2D_onChunk                :: transpiration
    dsl4jsb_Real2D_onChunk                :: w_skin
    dsl4jsb_Real2D_onChunk                :: w_soil_root
    dsl4jsb_Real2D_onChunk                :: w_soil_root_theta
    dsl4jsb_Real2D_onChunk                :: w_soil_root_pot
    dsl4jsb_Real2D_onChunk                :: interception 
    dsl4jsb_Real2D_onChunk                :: pet 
    dsl4jsb_Real2D_onChunk                :: evaporation 
    dsl4jsb_Real2D_onChunk                :: runoff
    dsl4jsb_Real2D_onChunk                :: drainage
    dsl4jsb_Real2D_onChunk                :: drainage_fraction
    dsl4jsb_Real2D_onChunk                :: ground_heat_flx_old
    dsl4jsb_Real2D_onChunk                :: ground_heat_flx
    dsl4jsb_Real2D_onChunk                :: latent_heat_flx
    dsl4jsb_Real2D_onChunk                :: sensible_heat_flx
    ! SPQ_ 3D        
    dsl4jsb_Real3D_onChunk                :: soil_depth_sl
    dsl4jsb_Real3D_onChunk                :: soil_lay_depth_ubound_sl
    dsl4jsb_Real3D_onChunk                :: bulk_dens_sl
    dsl4jsb_Real3D_onChunk                :: heat_capa_sl
    dsl4jsb_Real3D_onChunk                :: therm_cond_sl
    dsl4jsb_Real3D_onChunk                :: w_soil_pwp_sl
    dsl4jsb_Real3D_onChunk                :: w_soil_fc_sl 
    dsl4jsb_Real3D_onChunk                :: w_soil_sat_sl 
    dsl4jsb_Real3D_onChunk                :: saxtonA 
    dsl4jsb_Real3D_onChunk                :: saxtonB 
    dsl4jsb_Real3D_onChunk                :: saxtonC 
    dsl4jsb_Real3D_onChunk                :: kdiff_sat_sl 
    dsl4jsb_Real3D_onChunk                :: w_soil_sl
    dsl4jsb_Real3D_onChunk                :: w_soil_pot_sl
    dsl4jsb_Real3D_onChunk                :: t_soil_sl
    dsl4jsb_Real3D_onChunk                :: percolation_sl
    dsl4jsb_Real3D_onChunk                :: drainage_sl
    ! ASSIMI_
    dsl4jsb_Real2D_onChunk                :: canopy_cond
    ! RAD_
    dsl4jsb_Real2D_onChunk                :: sw_srf_net
    dsl4jsb_Real2D_onChunk                :: net_radiation
    ! PHENO_ 2D
    dsl4jsb_Real2D_onChunk                :: lai 
    ! VEG_ 2D
    dsl4jsb_Real2D_onChunk                :: height
    ! VEG_ 3D
    dsl4jsb_Real3D_onChunk                :: root_fraction_sl
    ! Get local variables from options argument
    iblk    = options%iblk
    ics     = options%ics
    ice     = options%ice
    nc      = options%nc
    steplen = options%steplen
    ! ---------------------------
    ! 0.4 Process Activity, Debug Option
    IF (.NOT. tile%Is_process_calculated(SPQ_)) RETURN
    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')
    ! ---------------------------
    ! 0.5 Get Memory
    model         => Get_model(tile%owner_model_id)
    lctlib        => model%lctlib(tile%lcts(1)%lib_id)
    vgrid_soil_sb => Get_vgrid('soil_layer_sb')
    nsoil_sb      =  vgrid_soil_sb%n_levels
    ! Get process config
    ! ..
    ! Get process memories
    dsl4jsb_Get_memory(A2L_) 
    dsl4jsb_Get_memory(ASSIMI_)
    dsl4jsb_Get_memory(RAD_)
    dsl4jsb_Get_memory(SPQ_)
    dsl4jsb_Get_memory(VEG_)
    dsl4jsb_Get_memory(PHENO_)
    ! Get process variables (Set pointers to variables in memory)
    dsl4jsb_Get_var2D_onChunk(A2L_, t_air)                ! in
    dsl4jsb_Get_var2D_onChunk(A2L_, q_air)                ! in
    dsl4jsb_Get_var2D_onChunk(A2L_, press_srf)            ! in
    dsl4jsb_Get_var2D_onChunk(A2L_, wind_10m)             ! in
    dsl4jsb_Get_var2D_onChunk(A2L_, lw_srf_down)          ! in
    dsl4jsb_Get_var2D_onChunk(A2L_, rain)                 ! in
    dsl4jsb_Get_var2D_onChunk(A2L_, snow)                 ! in
    dsl4jsb_Get_var2D_onChunk(A2L_, t_acoef)              ! in
    dsl4jsb_Get_var2D_onChunk(A2L_, t_bcoef)              ! in
    dsl4jsb_Get_var2D_onChunk(A2L_, q_acoef)              ! in
    dsl4jsb_Get_var2D_onChunk(A2L_, q_bcoef)              ! in
    dsl4jsb_Get_var2D_onChunk(A2L_, ga)                   !    not used actually !
    dsl4jsb_Get_var2D_onChunk(A2L_, pch)                  ! inout
    dsl4jsb_Get_var2D_onChunk(A2L_, drag_srf)             ! inout
    ! SPQ_ 2D
    dsl4jsb_Get_var2D_onChunk(SPQ_, fact_q_air)           ! inout
    dsl4jsb_Get_var2D_onChunk(SPQ_, fact_qsat_srf)        ! inout
    dsl4jsb_Get_var2D_onChunk(SPQ_, t_srf_old)            ! inout
    dsl4jsb_Get_var2D_onChunk(SPQ_, zril_old)             ! inout
    dsl4jsb_Get_var2D_onChunk(SPQ_, root_depth)           ! in
    dsl4jsb_Get_var2D_onChunk(SPQ_, w_soil_root_pwp)      ! in
    dsl4jsb_Get_var2D_onChunk(SPQ_, w_soil_root_fc)       ! in
    dsl4jsb_Get_var2D_onChunk(SPQ_, transpiration)        ! inout
    dsl4jsb_Get_var2D_onChunk(SPQ_, w_soil_root_theta)    ! inout
    dsl4jsb_Get_var2D_onChunk(SPQ_, w_soil_root_pot)      ! inout
    dsl4jsb_Get_var2D_onChunk(SPQ_, w_skin)               ! inout
    dsl4jsb_Get_var2D_onChunk(SPQ_, w_soil_root)          ! inout
    dsl4jsb_Get_var2D_onChunk(SPQ_, t_soil_root)          ! out
    dsl4jsb_Get_var2D_onChunk(SPQ_, pet)                  ! inout
    dsl4jsb_Get_var2D_onChunk(SPQ_, interception)         ! inout
    dsl4jsb_Get_var2D_onChunk(SPQ_, evaporation)          ! inout
    dsl4jsb_Get_var2D_onChunk(SPQ_, runoff)               ! inout
    dsl4jsb_Get_var2D_onChunk(SPQ_, drainage)             ! inout
    dsl4jsb_Get_var2D_onChunk(SPQ_, drainage_fraction)    ! inout
    dsl4jsb_Get_var2D_onChunk(SPQ_, ground_heat_flx_old)  ! inout
    dsl4jsb_Get_var2D_onChunk(SPQ_, ground_heat_flx)      ! inout
    dsl4jsb_Get_var2D_onChunk(SPQ_, latent_heat_flx)      ! inout
    dsl4jsb_Get_var2D_onChunk(SPQ_, sensible_heat_flx)    ! inout
    ! SPQ_ 3D
    dsl4jsb_Get_var3D_onChunk(SPQ_, soil_depth_sl)            ! in
    dsl4jsb_Get_var3D_onChunk(SPQ_, soil_lay_depth_ubound_sl) ! in
    dsl4jsb_Get_var3D_onChunk(SPQ_, bulk_dens_sl)             ! in
    dsl4jsb_Get_var3D_onChunk(SPQ_, heat_capa_sl)             ! in
    dsl4jsb_Get_var3D_onChunk(SPQ_, therm_cond_sl)            ! in
    dsl4jsb_Get_var3D_onChunk(SPQ_, w_soil_pwp_sl)            ! in
    dsl4jsb_Get_var3D_onChunk(SPQ_, w_soil_fc_sl)             ! in
    dsl4jsb_Get_var3D_onChunk(SPQ_, w_soil_sat_sl)            ! in
    dsl4jsb_Get_var3D_onChunk(SPQ_, saxtonA)                  ! in
    dsl4jsb_Get_var3D_onChunk(SPQ_, saxtonB)                  ! in
    dsl4jsb_Get_var3D_onChunk(SPQ_, saxtonC)                  ! in
    dsl4jsb_Get_var3D_onChunk(SPQ_, kdiff_sat_sl)             ! in
    dsl4jsb_Get_var3D_onChunk(SPQ_, w_soil_sl)                ! inout
    dsl4jsb_Get_var3D_onChunk(SPQ_, w_soil_pot_sl)            ! inout
    dsl4jsb_Get_var3D_onChunk(SPQ_, percolation_sl)           ! inout
    dsl4jsb_Get_var3D_onChunk(SPQ_, drainage_sl)              ! inout
    dsl4jsb_Get_var3D_onChunk(SPQ_, t_soil_sl)                ! inout
    ! ---------------------------
    dsl4jsb_Get_var2D_onChunk(ASSIMI_, canopy_cond)       ! in
    ! ---------------------------
    dsl4jsb_Get_var2D_onChunk(RAD_, sw_srf_net)           ! in
    dsl4jsb_Get_var2D_onChunk(RAD_, net_radiation)        ! inout
    ! PHENO_ 2D
    dsl4jsb_Get_var2D_onChunk(PHENO_, lai)                ! in 
    ! VEG_ 2D
    dsl4jsb_Get_var2D_onChunk(VEG_, height)               ! in 
    ! VEG_ 3D
    dsl4jsb_Get_var3D_onChunk(VEG_, root_fraction_sl)     ! in


    ! ------------------------------------------------------------------------------------------------------------
    ! Go science
    ! ------------------------------------------------------------------------------------------------------------

    !> 1.0 Diagnose surface relative humidity for energy budget calculations
    !!
    rel_hum(:) = MIN(1.0_wp, (w_soil_sl(:,1) - w_soil_pwp_sl(:,1)) / (w_soil_fc_sl(:,1) - w_soil_pwp_sl(:,1)))

    !> 2.0 surface fractions for JSBACH energy balance
    !! 


    !> 2.0 Calculate surface energy budget
    !!
    CALL calc_surface_energy_balance(nc, nsoil_sb, &
                                     steplen, &
                                     sw_srf_net(:),lw_srf_down(:), & 
                                     t_air(:),t_acoef(:),t_bcoef(:), & 
                                     q_air(:),q_acoef(:),q_bcoef(:), & 
                                     press_srf(:),wind_10m(:), & 
                                     height(:),lai(:),w_skin(:),canopy_cond(:),rel_hum(:),& 
                                     soil_depth_sl(:,:),heat_capa_sl(:,:),therm_cond_sl(:,:),bulk_dens_sl(:,:), &
                                     fact_q_air(:), fact_qsat_srf(:), &                             ! inout
                                     ga(:),drag_srf(:),pch(:), &                                    ! inout
                                     ground_heat_flx_old(:),zril_old(:),t_srf_old(:), &             ! inout
                                     pet(:),evaporation(:),interception(:),transpiration(:), &      ! inout
                                     sensible_heat_flx(:),latent_heat_flx(:),ground_heat_flx(:), &  ! inout
                                     net_radiation(:),t_soil_sl(:,:))                               ! inout

    !> 3.0 Calculate soil hydrology processes
    !!
    CALL calc_soil_hydrology(nc, nsoil_sb, &
                             soil_lay_depth_ubound_sl(:,:), &                                         ! (upper bound > lower bound !)
                             lctlib%phi_leaf_min, &
                             rain(:),snow(:),interception(:),evaporation(:),transpiration(:), & 
                             w_soil_root_fc(:),w_soil_root_pwp(:),lai(:),root_depth(:), & 
                             soil_depth_sl(:,:),w_soil_sat_sl(:,:),w_soil_fc_sl(:,:),w_soil_pwp_sl(:,:), & 
                             saxtonA(:,:),saxtonB(:,:),saxtonC(:,:),kdiff_sat_sl(:,:),root_fraction_sl(:,:),&
                             w_skin(:),w_soil_root(:),w_soil_root_theta(:),w_soil_root_pot(:), &      ! inout
                             runoff(:),drainage(:),drainage_fraction(:), &                            ! inout
                             w_soil_sl(:,:),w_soil_pot_sl(:,:),drainage_sl(:,:),percolation_sl(:,:))  ! inout

    !> 4.0 diagnose root temperature
    !!
    t_soil_root(:) = 0.0_wp
    DO isoil = 1,nsoil_sb
       t_soil_root(:) = t_soil_root(:) + root_fraction_sl(:,isoil) * t_soil_sl(:,isoil)
    ENDDO

  END SUBROUTINE update_spq_physics

  !-----------------------------------------------------------------------------------------------------
  ! Sub Task to update_spq_physics
  !
  !-----------------------------------------------------------------------------------------------------
  !> Routine to provide an update of the soil temperature given the net radiation, surface drag, \n
  !!   surface wetness, and canopy conductance. 
  !!  
  !!
  !! Calculates surface energy fluxes and deep soil temperatures 
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE calc_surface_energy_balance(nc, nsoil_sb, &
                                         steplen, &
                                         sw_srf_net,lw_srf_down, & 
                                         t_air,t_acoef,t_bcoef, & 
                                         q_air,q_acoef,q_bcoef, & 
                                         press_srf,wind_10m, & 
                                         height,lai,w_skin,canopy_cond,rel_hum,& 
                                         soil_depth_sl,heat_capa_sl,therm_cond_sl,bulk_dens_sl, &
                                         fact_q_air,fact_qsat_srf, &                                                ! inout
                                         ga,drag_srf,pch, &                                                         ! inout
                                         ground_heat_flx_old,zril_old,t_srf_old, &                                  ! inout
                                         pet,evaporation,interception,transpiration, &                              ! inout
                                         sensible_heat_flx,latent_heat_flx,ground_heat_flx,net_radiation,t_soil_sl) ! inout

    USE mo_jsb_math_constants,        ONLY: dtime, eps8
    USE mo_atmland_constants,         ONLY: min_wind, max_wind 
    USE mo_atmland_util,              ONLY: calc_spec_humidity_sat
    USE mo_spq_constants,             ONLY: w_skin_max
    USE mo_jsb_physical_constants,    ONLY: r_gas_dryair, LatentHeatSublimation, LatentHeatVaporization, &
                                            stbo, zemiss_def, grav, cvdifts, cpd, cpvd1
    USE mo_phy_schemes,               ONLY: heat_transfer_coef, q_effective, surface_dry_static_energy
    

    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut 
    INTEGER,                                          INTENT(in)    :: nc, &                  !< dimensions
                                                                       nsoil_sb
    REAL(wp),                                         INTENT(in)    :: steplen                !< step length for this process
    REAL(wp), DIMENSION(nc),                          INTENT(in)    :: sw_srf_net, &          !< net shortwave radiation [W/m2] 
                                                                       lw_srf_down, &         !< longwave downward radiation [W/m2]
                                                                       t_air, &               !< air temperature [K]
                                                                       t_acoef, &             !< Richtmeyer-Morton Coefficient Ta
                                                                       t_bcoef, &             !< Richtmeyer-Morton Coefficient Tb
                                                                       q_air, &               !< air specific humidity [g/g]
                                                                       q_acoef, &             !< Richtmeyer-Morton Coefficient Qa
                                                                       q_bcoef, &             !< Richtmeyer-Morton Coefficient Qb
                                                                       press_srf, &           !< surface pressure [Pa]
                                                                       wind_10m, &            !< wind speed [m/2]
                                                                       height, &              !< height
                                                                       lai, &                 !< LAI
                                                                       w_skin, &              !< skin water reservoir on leaves [m]
                                                                       canopy_cond, &         !< canopy conductance for water [m/s]
                                                                       rel_hum                !< relative humidity for soil evaporation
    REAL(wp), DIMENSION(nc,nsoil_sb),                 INTENT(in)    :: soil_depth_sl, &       !< soil depth per layer [m]
                                                                       heat_capa_sl, &        !< heat capacity of the soil layer [J/kg]
                                                                       therm_cond_sl, &       !< thermal conductivity of the soil layer []
                                                                       bulk_dens_sl           !< bulk density of the soil layer [kg/m3]
    REAL(wp), DIMENSION(nc),                          INTENT(inout) :: ground_heat_flx_old, & !< ground heat flux from previous time step [W/m2]
                                                                       t_srf_old, &           !< previous timestep's surface temperature 
                                                                       zril_old, &            !< previous timestep's Reynolds number 
                                                                       fact_q_air, &          !< 
                                                                       fact_qsat_srf, &       !< 
                                                                       ga, &                  !< aerodynamic conductance [m/s] 
                                                                       drag_srf, &            !< surface drag coefficient (JSBACH-style!)
                                                                       pch, &                 !< heat exchange parameter (JSBACH-style!)
                                                                       pet, &                 !< potential evaporation rate [kg / m2 /s]
                                                                       evaporation, &         !< evaporation rate [kg / m2 /s]
                                                                       interception, &        !< interception loss rate [kg / m2 /s]
                                                                       transpiration, &       !< transpiration rate [kg / m2 /s]
                                                                       sensible_heat_flx, &   !< sensible heat flux [W/m2]
                                                                       latent_heat_flx, &     !< latent heat flux [W/m2]
                                                                       ground_heat_flx, &     !< ground heat flux [W/m2]
                                                                       net_radiation          !< net radiation [W/m2]
    REAL(wp), DIMENSION(nc,nsoil_sb),                 INTENT(inout) :: t_soil_sl              !< soil temperature [K]
    ! ---------------------------
    ! 0.2 Local
    REAL(wp), DIMENSION(nc)                       ::     q_sat        , &
                                                         q_sat_srf    , & 
                                                         heat_tcoef   , & 
                                                         wind_air     , & 
                                                         t_srf_new    , & 
                                                         hlp1 
    REAL(wp), DIMENSION(nc,nsoil_sb)              :: ground_heat_flx_sl              !< soil temperature [K]
    ! LOGICAL, DIMENSION(nc)                        :: info_catch_t_jump    ! for debug/ write out temperature jumps larger 15 K
    ! INTEGER, SAVE                                 :: timesteps_help = 0   ! for debug/ write out temperature jumps larger 15 K
    REAL(wp), DIMENSION(nc)                       :: dQdT,s_old,t2s_conv,s_star, q_sat_star
    REAL(wp), DIMENSION(nc)                       :: q_air_eff,q_sat_srf_eff,evapotranspiration
    REAL(wp), DIMENSION(nc)                       :: fract_snow,fract_water,fract_fpc,fact_qsat_trans_srf
    REAL(wp), DIMENSION(nc)                       :: evaporation_snow,evaporation_bare

    INTEGER                                       :: isoil, it   ! looping
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_surface_energy_balance'


    ! ------------------------------------------------------------------------------------------------------------
    ! Go science
    ! ------------------------------------------------------------------------------------------------------------


    !> 1.0 calculate surface characteristic given old surface state 
    !!    
    q_sat(:)      = calc_spec_humidity_sat(t_air(:), press_srf(:))
    q_sat_srf(:)  = calc_spec_humidity_sat(t_soil_sl(:,1), press_srf(:))
    wind_air(:)   = MIN(MAX(wind_10m(:), min_wind), max_wind)     

    !> 1.1 Calculate heat transfer coefficient
    !!
    heat_tcoef(:) = heat_transfer_coef(drag_srf(:), steplen)
    
    !> 1.2 net radiation flux at the surface
    !!
    net_radiation(:) = sw_srf_net(:) + lw_srf_down(:) - stbo * zemiss_def * t_soil_sl(:,1)**4._wp  
    pet(:)           = -1._wp * heat_tcoef(:) * (q_air(:) - q_sat_srf(:))

    !> 1.3 calculate snow and water-saturated fraction of tile
    !!     following JSBACH4 hydrology 
    hlp1(:) = 0.0_wp
    CALL calc_wskin_fractions_veg(dtime, &
                                  steplen, & 
                                  w_skin_max * lai(:), &
                                  0.01_wp, & ! assume almost flat terrain
                                  pet(:), &
                                  w_skin, &
                                  0.0_wp, &     ! jsb4: w_snow_soil
                                  0.0_wp, &     ! jsb4: w_snow_can
                                  hlp1(:), &          ! out
                                  fract_water(:), &   ! out
                                  fract_snow(:))      ! out
    fract_fpc = 1.0_wp - EXP(-0.5_wp * lai(:)) 

    !> 1.4 calculate humidity scaling factors for energy balance calculation 
    !!
    CALL update_humidity_scaling(nc, &
                                 fract_snow, fract_water, fract_fpc, &
                                 q_sat_srf, q_air, &
                                 wind_air, pch, &
                                 canopy_cond, rel_hum, &
                                 fact_q_air, fact_qsat_srf, fact_qsat_trans_srf)

    !> 2.0 JSBACH4: update_surface_energy_land

    !> 2.1 Compute sensitivity of surface saturation specific humidity to temperature
    !!
    dQdT(:) = (calc_spec_humidity_sat(t_soil_sl(:,1) + 0.001_wp, press_srf(:)) - q_sat_srf(:)) * 1000._wp

    !> 2.2 Old dry static energy
    !!
    s_old(:) = surface_dry_static_energy( &
      & t_soil_sl(:,1), &
      & q_effective(q_sat_srf(:), q_air(:), fact_qsat_srf(:), fact_q_air(:)))
    t2s_conv(:) = s_old(:) / t_soil_sl(:,1) 

    !> 2.3 solve energy balance calculation using implicit formulation (JSBACH4)
    !!
    hlp1 = calc_temp_srf_new_implicit(steplen, &
                                      t2s_conv(:), & 
                                      t_acoef(:), &
                                      t_bcoef(:), &
                                      q_acoef(:), &
                                      q_bcoef(:), & 
                                      s_old(:), &
                                      q_sat_srf(:), &
                                      dQdT(:), &
                                      net_radiation(:), &
                                      ground_heat_flx(:), & 
                                      heat_tcoef(:), &
                                      fact_q_air(:), &
                                      fact_qsat_srf(:), &
                                      fract_snow(:), & 
                                      heat_capa_sl(:,1) * bulk_dens_sl(:,1) * soil_depth_sl(:,1)) 

    ! time-filtering as in JSBACH3
    s_star(:)     = 1._wp / cvdifts * hlp1(:) + (1._wp - 1._wp / cvdifts) * s_old(:) 
    t_srf_new(:)  = s_star(:) / t2s_conv(:)
    q_sat_star(:) = q_sat_srf(:) + dQdT(:) * (s_star(:) - s_old(:)) / t2s_conv(:) 

    !> 3.0 calculate water fluxes
    !!   This is from JSBACH4 hydrology/mo_hydro_interface:update_evaporation
    q_air_eff(:)          = q_effective(0.0_wp, q_air(:), 1.0_wp, 0.0_wp)
    q_sat_srf_eff(:)      = q_effective(q_sat_star(:), q_air(:), fact_qsat_srf(:), fact_q_air(:))
    evapotranspiration(:) = heat_tcoef(:) * (q_air_eff(:) - q_sat_srf_eff(:)) 
    pet(:)                = heat_tcoef(:) * (q_air(:) - q_sat_star(:))
  
    !> 3.1 convert to positive fluxes for water module
    !! needed for QUINCY, but not JSBACH4
    transpiration(:)    = -1._wp * fact_qsat_trans_srf(:) * pet(:) 
    interception(:)     = -1._wp * (1._wp - fract_snow(:)) * fract_water(:) * pet(:) 
    evaporation_snow(:) = -1._wp * fract_snow(:) * pet(:) 
    ! no dew into stomates or onto canopy
    WHERE(pet(:) > 0.0_wp)
       transpiration(:)    = 0.0_wp
       interception(:)     = 0.0_wp
       evaporation_snow(:) = 0.0_wp
    ENDWHERE
    evaporation_bare(:) = -1._wp * evapotranspiration(:) - transpiration(:) - interception(:) - evaporation_snow(:)
    ! BPLARR @TODO, as long as there is no snow
    evaporation = evaporation_bare

    !> 4.0 calculate sensible and latent heat flux
    !!  from JSBACH4 srf_energy_bal/mo_seb_land:update_surface_fluxes_land

    ! Compute sensible heat flux
    sensible_heat_flx(:) = heat_tcoef(:) * (t_acoef(:) * s_old(:) + t_bcoef(:) - s_star(:))

    ! Compute latent heat flux
    latent_heat_flx(:) = LatentHeatVaporization * evapotranspiration(:) + & 
                         (LatentHeatSublimation - LatentHeatVaporization) * fract_snow(:) * pet(:)

    !> 5.0 update surface temperature and remember previous value
    !>
    !>   catching temperature changes to be not larger than 15 K (hardcoded)
    !>   
    ! timesteps_help = timesteps_help +1    ! uncomment for debug only (also see declaration above)
    ! info_catch_t_jump(:) = .FALSE.        ! uncomment for debug only (also see declaration above)
    t_srf_old(:) = t_soil_sl(:,1)
    ! catch 15 K    
    WHERE (ABS(t_soil_sl(:,1) - t_srf_new(:)) > 15.0_wp)
      ! info_catch_t_jump(:) = .TRUE.       ! uncomment for debug only (also see declaration above)
      WHERE (t_soil_sl(:,1) < t_srf_new(:))
        t_soil_sl(:,1) = t_soil_sl(:,1) + 15.0_wp
      ELSEWHERE
        t_soil_sl(:,1) = t_soil_sl(:,1) - 15.0_wp
      END WHERE
    ELSEWHERE
      t_soil_sl(:,1) = t_srf_new(:)
    END WHERE
    ! ! write info if one or more values are TRUE -> uncomment for debug only
    ! IF(ANY(info_catch_t_jump(:))) THEN
    !   write(960,'(I14,20F26.8)') timesteps_help, t_srf_old(:), t_srf_new(:) - t_srf_old(:)
    ! END IF

    !> 5.1 calculate ground heat flux and update deep surface temperatures 
    !!
    ground_heat_flx_sl(:,:) = 0.0_wp 
    ! loop layers 1:(max-1)
    DO isoil = 1,nsoil_sb-1
      WHERE(soil_depth_sl(:,isoil) > eps8)
           ground_heat_flx_sl(:,isoil) = therm_cond_sl(:,isoil) * & 
                                         (t_soil_sl(:,isoil+1)-t_soil_sl(:,isoil)) / & 
                                         ((soil_depth_sl(:,isoil)+soil_depth_sl(:,isoil+1))/2.0_wp)
      END WHERE 
    ENDDO 
    ground_heat_flx(:)     = ground_heat_flx_sl(:,1) 
    ground_heat_flx_old(:) = ground_heat_flx(:)
    ! loop layers 2:max
    DO isoil = 2,nsoil_sb
      WHERE(soil_depth_sl(:,isoil) > eps8)
         t_soil_sl(:,isoil) = t_soil_sl(:,isoil) + &
                              (ground_heat_flx_sl(:,isoil) - ground_heat_flx_sl(:,isoil-1)) * dtime / & 
                              (heat_capa_sl(:,isoil) * bulk_dens_sl(:,isoil) * soil_depth_sl(:,isoil)) 
      END WHERE 
    ENDDO 
 
  END SUBROUTINE calc_surface_energy_balance


  !-----------------------------------------------------------------------------------------------------
  ! Sub Task to update_spq_physics
  !
  !-----------------------------------------------------------------------------------------------------
  !> Routine to provide an update of the soil water given rates of rain and snow fall \n
  !!  as well as interception loss, evaporation and transpiration. 
  !!
  !! Calculates runoff, drainage and an updated state of soil water variables (theta, potential, etc.)
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE calc_soil_hydrology(nc, nsoil_sb, &
                                 soil_layer_w_ubounds, & 
                                 lctlib_phi_leaf_min, &
                                 rain,snow,interception,evaporation,transpiration, & 
                                 w_soil_root_fc,w_soil_root_pwp,lai,root_depth, & 
                                 soil_depth_sl,w_soil_sat_sl,w_soil_fc_sl,w_soil_pwp_sl, & 
                                 saxtonA,saxtonB,saxtonC,kdiff_sat_sl,root_fraction_sl,&
                                 w_skin,w_soil_root,w_soil_root_theta,w_soil_root_pot, &  ! inout
                                 runoff,drainage,drainage_fraction, &                     ! inout
                                 w_soil_sl,w_soil_pot_sl,drainage_sl,percolation_sl)      ! inout

    USE mo_jsb_math_constants,          ONLY: dtime, eps8
    USE mo_jsb_physical_constants,      ONLY: grav
    USE mo_assimi_update_canopy_fluxes, ONLY: calc_soil_moisture_stress
    USE mo_spq_constants,               ONLY: drain_frac, fdrain_runoff, InterceptionEfficiency, w_skin_max, wpot_min


    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut 
    INTEGER,                                          INTENT(in)    :: nc, &                   !< dimensions
                                                                       nsoil_sb
    REAL(wp), DIMENSION(nc,nsoil_sb),                 INTENT(in)    :: soil_layer_w_ubounds    !< upper bound of soil layers water (upper bound > lower bound !)
    REAL(wp),                                         INTENT(in)    :: lctlib_phi_leaf_min     !< lctlib paramter
    REAL(wp), DIMENSION(nc),                          INTENT(in)    :: rain, &                 !< rainfall [kg / m2 / s ]
                                                                       snow, &                 !< snowfall [kg / m2 / s ]
                                                                       interception, &         !< interception [kg / m2 / s ]
                                                                       evaporation, &          !< evaporation [kg / m2 / s ]
                                                                       transpiration, &        !< transpiration [kg / m2 / s ]
                                                                       w_soil_root_fc, &       !< root moisture at field capacity [m]
                                                                       w_soil_root_pwp, &      !< root moisture at wilting point [m]
                                                                       lai, &                  !< leaf area index
                                                                       root_depth              !< rooting depth [m]
    REAL(wp), DIMENSION(nc,nsoil_sb),                 INTENT(in)    :: soil_depth_sl, &        !< soil depth per layer [m]
                                                                       w_soil_sat_sl, &        !< soil moisture per layer at saturation [m]
                                                                       w_soil_fc_sl, &         !< soil moisture per layer at field capacity [m]
                                                                       w_soil_pwp_sl, &        !< soil moisture per layer at wilting point [m]
                                                                       kdiff_sat_sl, &         !< saturated water diffusivity [m/s] 
                                                                       saxtonA, &              !< pedotransfer parameter 
                                                                       saxtonB, &              !< pedotransfer parameter
                                                                       saxtonC, &              !< pedotransfer parameter
                                                                       root_fraction_sl        !< fraction of root mass per soil layer
    REAL(wp), DIMENSION(nc),                          INTENT(inout) :: w_skin, &               !< water in skin reservoir [m]
                                                                       w_soil_root, &          !< water in root zone [m]
                                                                       w_soil_root_theta, &    !< fractional water availability in root zone [-]
                                                                       w_soil_root_pot, &      !< soil water potential in root zone [MPa]
                                                                       runoff, &               !< surface runoff [kg / m2 / s]
                                                                       drainage, &             !< deep drainage [kg / m2 / s]
                                                                       drainage_fraction       !< fraction of soil column lost to drainage [1/s]
    REAL(wp), DIMENSION(nc,nsoil_sb),                 INTENT(inout) :: w_soil_sl, &            !< soil moisture per layer [m]
                                                                       w_soil_pot_sl, &        !< soil water potential per layer [MPa]
                                                                       drainage_sl, &          !< drainage per layer [kg / m2 / s]
                                                                       percolation_sl          !< fractional water loss [1000/s] 
    ! ---------------------------
    ! 0.2 Local
    REAL(wp), DIMENSION(nc)                         :: &            
                                                       infiltration , &
                                                       water_to_skin, &
                                                       water_to_soil, &
                                                       w_soil_root_old, &
                                                       hlp1, hlp2, hlp3
    REAL(wp),DIMENSION(nc,nsoil_sb)                 :: &   
                                                       ftranspiration , &
                                                       fevaporation   , &
                                                       w_soil_sl_old  , &
                                                       k_diff         , &
                                                       hlp1_2d
    INTEGER                                         :: isoil, ichunk 
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_soil_hydrology'


    ! ------------------------------------------------------------------------------------------------------------
    ! Go science
    ! ------------------------------------------------------------------------------------------------------------


    !> 0.9 safe old water content for conservation tests
    !!
    w_soil_root_old(:) = w_soil_root(:)
    w_soil_sl_old(:,:) = w_soil_sl(:,:)

    !> 1.0 Update skin reservoir and calculate interception loss
    !!
    WHERE((w_skin_max * lai(:)) > eps8)
       ! m / timestep
       water_to_skin(:) = MIN(InterceptionEfficiency * (rain(:) + snow(:)) / 1000.0_wp * dtime, & 
                              MAX(interception(:) * dtime / 1000.0_wp + w_skin_max * lai(:) - w_skin(:), 0._wp))
    ELSEWHERE
       water_to_skin(:) = 0.0_wp
    ENDWHERE
    water_to_soil(:) = MAX(rain(:) + snow(:) - water_to_skin(:) / dtime * 1000.0_wp, 0._wp)
    w_skin(:)        = MAX(w_skin(:) + water_to_skin(:) - interception(:) * dtime / 1000.0_wp, 0.0_wp)


    !> 2.0 update soil moisture
    !!
    fevaporation(:,:) = 0.0_wp
    fevaporation(:,1) = 1.0_wp

    !> 2.1 determine fraction of transpiration per layer
    !!
    ftranspiration(:,:) = 0.0_wp
    DO isoil = 1,nsoil_sb
       WHERE((w_soil_fc_sl(:,isoil) - w_soil_pwp_sl(:,isoil)) > eps8)
          hlp1(:) = calc_soil_moisture_stress(w_soil_pot_sl(:, isoil), lctlib_phi_leaf_min)
          ftranspiration(:,isoil) = root_fraction_sl(:,isoil) * hlp1(:) 
       END WHERE
    ENDDO
    DO ichunk = 1,nc
      ! here SUM covers all soil layers per point of the chunk
      IF( SUM(ftranspiration(ichunk,:)) > eps8 ) THEN 
         ftranspiration(ichunk,:) = ftranspiration(ichunk,:) / SUM(ftranspiration(ichunk,:))
      ELSE 
         ftranspiration(ichunk,:) = 0.0_wp
      ENDIF
    END DO

    !> 2.2 determine water balance of the top layer 
    !! 
    !! infiltration rate (limited to attaining field capacity), surface runoff and deep drainage
    !! a fraction of runoff generation is leaked below the first layer
    !! a fraction of infiltration is leaked below the first layer (preferential flow) 
    ! first guess
    infiltration(:) = MAX( MIN( (w_soil_fc_sl(:,1) - w_soil_sl(:,1))*1000._wp/dtime, water_to_soil(:)), 0.0_wp)
    runoff(:)       = MAX(water_to_soil(:) - infiltration(:), 0.0_wp)

    ! modify surface runoff and infiltration such that a constant fraction is actually leaked below the first layer
    drainage_sl(:,1) = drain_frac / soil_depth_sl(:,1) * infiltration(:) + fdrain_runoff * runoff(:)
    runoff(:)        = ( 1.0_wp - fdrain_runoff ) * runoff(:)  
    infiltration(:)  = ( 1.0_wp - drain_frac / soil_depth_sl(:,1) ) * infiltration(:)

    w_soil_sl(:,1)   = w_soil_sl(:,1) + &
                       (infiltration(:) - evaporation(:) * fevaporation(:,1) - transpiration(:) * ftranspiration(:,1)) / &
                       1000._wp * dtime
    WHERE(w_soil_sl(:,1) > w_soil_fc_sl(:,1))
         runoff(:)      = runoff(:) + (w_soil_sl(:,1) - w_soil_fc_sl(:,1))*1000._wp / dtime
         w_soil_sl(:,1) = w_soil_fc_sl(:,1) 
    END WHERE 

    !> 2.3 update deeper soil layers given infiltration, transpiration and evaporation
    !!
    DO isoil = 2,nsoil_sb

      ! drainage from this layer due to preferential flow of infiltering water flux
      drainage_sl(:,isoil) = drain_frac / soil_depth_sl(:,isoil) * drainage_sl(:,isoil-1)

      ! water balance of soil layer
      w_soil_sl(:,isoil)   = w_soil_sl(:,isoil) + & 
                             (( 1.0_wp - drain_frac / soil_depth_sl(:,isoil) ) * drainage_sl(:,isoil-1) - &
                              evaporation(:) * fevaporation(:,isoil) - &
                              transpiration(:) * ftranspiration(:,isoil)) / &
                             1000._wp * dtime

      ! drain surplus water to lower layer to keep water storage below field capacity
      WHERE(w_soil_sl(:,isoil) > w_soil_fc_sl(:,isoil))
        drainage_sl(:,isoil) = drainage_sl(:,isoil) + (w_soil_sl(:,isoil) - w_soil_fc_sl(:,isoil))*1000._wp / dtime
        w_soil_sl(:,isoil)   = w_soil_fc_sl(:,isoil) 
      END WHERE 

    ENDDO

    !> 2.4 update deeper soil layer by diffusion
    !!
    !! lowest value of w_soil_sl given minimum number allowed as outcome of theta**C
    hlp1_2d(:,:) = exp(log(eps8)/saxtonC(:,:)) * w_soil_sat_sl(:,:)
    WHERE(w_soil_sat_sl(:,:) > eps8 .AND. w_soil_sl(:,:) > hlp1_2d(:,:))
       k_diff(:,:) = kdiff_sat_sl(:,:) * (w_soil_sl(:,:) / w_soil_sat_sl(:,:))**saxtonC(:,:)
    ELSEWHERE
       k_diff(:,:) = 0.0_wp
    END WHERE          

    !! update soil water potential for diffusion calculation
    WHERE(soil_depth_sl(:,:) > eps8)
       ! Campbell applies only water content below field capacity
       hlp1_2d(:,:) = MIN(w_soil_fc_sl(:,:),w_soil_sl(:,:)) / soil_depth_sl(:,:)
       ! for numerical reasons, constrain to minimum water potential 
       hlp1_2d(:,:) = MAX(exp(log(wpot_min/saxtonA(:,:)) /saxtonB(:,:)), hlp1_2d(:,:))
    ELSEWHERE
       hlp1_2d(:,:) = exp(log(wpot_min/saxtonA(:,:)) /saxtonB(:,:))
    ENDWHERE
    w_soil_pot_sl(:,:) = saxtonA(:,:) * (hlp1_2d(:,:) ** saxtonB(:,:))

    DO isoil = 2,nsoil_sb
       ! weighted mean of diffusivity (m/s)
       hlp1(:) = (k_diff(:,isoil-1) * soil_depth_sl(:,isoil-1) + k_diff(:,isoil) * soil_depth_sl(:,isoil)) / &
                 (soil_depth_sl(:,isoil-1) + soil_depth_sl(:,isoil))
       ! flux resulting from the pressure gradient between the centres of the soil layer
       ! note that matrix potential is converted to pressure head MPa -> J/N = m
       hlp1(:) = hlp1(:) * dtime * & 
                 (w_soil_pot_sl(:,isoil-1) - w_soil_pot_sl(:,isoil)) * &
                 1000.0_wp / grav / & 
                 (0.5_wp * (soil_depth_sl(:,isoil-1) + soil_depth_sl(:,isoil)))

       ! limit diffusion rate to the rate that produces equal fractional soil moisture in the two layers
       hlp3(:) = soil_depth_sl(:,isoil-1)/soil_depth_sl(:,isoil)             
       hlp2(:) = (w_soil_sl(:,isoil-1) - w_soil_sl(:,isoil) * hlp3(:)) / (1._wp + hlp3(:))
       WHERE(hlp1(:) < 0.0_wp  .AND. hlp1(:) < hlp2(:)) hlp1(:) = hlp2(:)
       WHERE(hlp1(:) >= 0.0_wp .AND. hlp1(:) > hlp2(:)) hlp1(:) = hlp2(:)
       ! avoid soil water content exceeding soil water capacity
       hlp1(:) = MIN(MAX(hlp1(:), w_soil_sl(:,isoil-1) - w_soil_fc_sl(:,isoil-1)), &
                     w_soil_fc_sl(:,isoil)-w_soil_sl(:,isoil))

       ! update soil layers with diffusion
       w_soil_sl(:,isoil-1)   = w_soil_sl(:,isoil-1) - hlp1(:)
       w_soil_sl(:,isoil)     = w_soil_sl(:,isoil) + hlp1(:) 
       drainage_sl(:,isoil-1) = drainage_sl(:,isoil-1) + hlp1(:) * 1000._wp / dtime 
    
    ENDDO

    !> 2.5. diagnose deep drainage
    !!
    drainage(:) = drainage_sl(:,nsoil_sb) 
    WHERE(w_soil_sl_old(:,:) > eps8)
      percolation_sl(:,:) = MIN( MAX((drainage_sl(:,:)/1000_wp*dtime) / w_soil_sl_old(:,:), -1.0_wp), 1.0_wp) * &
                            1000._wp/dtime
    ELSEWHERE
      percolation_sl(:,:) = 0.0_wp
    ENDWHERE


    !> 3.0 final diagnostics
    !!

    !> 3.1 calculate rooting zone water content and theta for compatability 
    !!       vegetation water stress calculation.
    !!       If a soil layer is only partially filled with roots, the plant-accessible 
    !!       water is limited to the AWC of the root zone. All of the layer's water 
    !!       (fc-w) deficit is substracted from the root accessible AWC (fc-pwp)*root_depth
    !! 
    !! note: upper bound is the larger value compared to lower bound !
    !! 
    w_soil_root(:) = w_soil_sl(:,1) 
    DO isoil = 2,nsoil_sb
       WHERE(root_depth(:) > soil_layer_w_ubounds(:,isoil-1))
          WHERE(root_depth(:) >= soil_layer_w_ubounds(:,isoil))
            w_soil_root(:) = w_soil_root(:) + w_soil_sl(:,isoil) 
          ELSEWHERE
            ! fraction of the layer that has roots
            hlp1(:) = (root_depth(:) - soil_layer_w_ubounds(:,isoil-1)) / &
                      (soil_layer_w_ubounds(:,isoil) - soil_layer_w_ubounds(:,isoil-1))

            w_soil_root(:) = w_soil_root(:) + w_soil_fc_sl(:,isoil) * hlp1(:) - &
                             MIN((w_soil_fc_sl(:,isoil) - w_soil_sl(:,isoil)), w_soil_pwp_sl(:,isoil) * hlp1(:))
          ENDWHERE
       ENDWHERE
    ENDDO

    !> 3.2 calculate plant fractional water availability 
    !!
    w_soil_root_theta(:) = ( w_soil_root(:) - w_soil_root_pwp(:) ) / &
                           ( w_soil_root_fc(:) - w_soil_root_pwp(:) ) 
    w_soil_root_theta(:) = MAX(w_soil_root_theta(:),0.0_wp)

    !> 3.3 calculate root zone and layer soil water potential according to Campbell et al 1974 
    !! 
    hlp1(:) = MIN(w_soil_root_fc(:),w_soil_root(:)) / root_depth(:)
    ! for numerical reasons, constrain to minimum water potential 
    hlp1(:) = MAX(exp(log(wpot_min/saxtonA(:,1))/saxtonB(:,1)),hlp1(:)) 
    w_soil_root_pot(:) = saxtonA(:,1) * (hlp1(:) ** saxtonB(:,1))

    WHERE(soil_depth_sl(:,:) > eps8)
       ! Campbell applies only water content below field capacity
       hlp1_2d(:,:) = MIN(w_soil_fc_sl(:,:),w_soil_sl(:,:)) / soil_depth_sl(:,:)
       ! for numerical reasons, constrain to minimum water potential 
       hlp1_2d(:,:) = MAX(exp(log(wpot_min/saxtonA(:,:))/saxtonB(:,:)),hlp1_2d(:,:)) 
    ELSEWHERE
       hlp1_2d(:,:) = exp(log(wpot_min/saxtonA(:,:))/saxtonB(:,:)) 
    ENDWHERE
    w_soil_pot_sl(:,:) = saxtonA(:,:) * (hlp1_2d(:,:) ** saxtonB(:,:))

    !> 3.4 Diagnostic fraction of water mass lost from the soil colum per second
    !!
    WHERE(w_soil_root_old(:) > eps8)
       drainage_fraction(:) = drainage(:) / 1000._wp / w_soil_root_old(:)
    ELSEWHERE
       drainage_fraction(:) = 0.0_wp 
    ENDWHERE

    ! for debug only:
    ! write(898,'(20F20.10)'),water_to_soil,interception,evaporation,transpiration,drainage,runoff,& 
    !                        (water_to_soil-evaporation-transpiration-drainage-runoff)/1000*dtime,(SUM(w_soil_sl)-SUM(w_soil_sl_old)
  END SUBROUTINE calc_soil_hydrology 


  !-----------------------------------------------------------------------------------------------------
  ! Sub Task to calc_surface_energy_balance 
  !
  !-----------------------------------------------------------------------------------------------------
  !> Routine to calculate snow and wet fraction of a tile 
  !!  mimicks the JSBACH4 routine calc_wskin_fractions_veg in 
  !!   hydroloogy/mo_hyrology_interface:calc_wskin_fractions_veg
  !!
  !! Calculates scaling factors fract_snow and fract_water 
  !-----------------------------------------------------------------------------------------------------
 ELEMENTAL SUBROUTINE calc_wskin_fractions_veg( &
    & dtime,                                     & ! in
    & steplen,                                   & ! in
    & wsr_max,                                   & ! in
    & oro_stddev,                                & ! in
    & evapo_pot,                                 & ! in
    & w_skin,                                    & ! in
    & w_snow_soil,                               & ! in
    & w_snow_can,                                & ! in
    & fract_snow_can,                            & ! out
    & fract_water    ,                           & ! out
    & fract_snow_soil                            & ! out
    & )

    USE mo_jsb_math_constants,     ONLY: eps8
    USE mo_jsb_physical_constants, ONLY: rhoh2o

    REAL(wp), INTENT(in) :: &
      & dtime,              &
      & steplen
    REAL(wp), INTENT(in) :: &
      & oro_stddev,         &
      & evapo_pot,          &
      & wsr_max,            &
      & w_skin,             &
      & w_snow_soil,        &
      & w_snow_can
    REAL(wp), INTENT(out) :: &
      & fract_snow_can,      &
      & fract_water    ,     &
      & fract_snow_soil


    ! Snow fraction on canopy
    ! Wet skin reservoir fraction on soil and canopy (between 0 and 1)
    IF (wsr_max > eps8) THEN
      fract_snow_can = MIN(1._wp, w_snow_can / wsr_max)
      fract_water    = MIN(1._wp, w_skin     / wsr_max)
    ELSE
      fract_snow_can = 0._wp
      fract_water    = 0._wp
    END IF
    fract_snow_can = 0._wp
    fract_snow_soil = 0._wp

    ! Snow fraction on soil
    fract_snow_soil = MERGE(fract_snow_can, fract_snow_soil, &
                               fract_snow_soil < EPSILON(1._wp) .AND. fract_snow_can > EPSILON(1._wp))


    ! Modify snow cover if snow loss during the time step due to potential evaporation is larger than
    ! snow water content from soil and canopy; same for skin reservoir
    IF (fract_snow_soil > 0._wp) THEN
      ! @TODO Shouldn't one take rhoice here insteady of rhoh2o?
      fract_snow_soil = fract_snow_soil / MAX(1._wp, fract_snow_soil * evapo_pot * dtime &
        &                                            / (rhoh2o * (w_snow_soil + w_snow_can)) &
        &                                    )
    END IF
    IF (fract_water > eps8 ) THEN
      fract_water = fract_water    / MAX(1._wp, (1._wp - fract_snow_soil) * evapo_pot * dtime &
        &                                       / (rhoh2o * MAX(EPSILON(1._wp), w_skin)) &
        &                                    )
    END IF

  END SUBROUTINE calc_wskin_fractions_veg

  !-----------------------------------------------------------------------------------------------------
  ! Sub Task to calc_surface_energy_balance
  !
  !-----------------------------------------------------------------------------------------------------
  !> Routine to update humidity scaling for the energy balance calculation 
  !!  mimicks the JSBACH4 routine update_humidity_scaling in 
  !!   turbulence/mo_turbulence_interface:update_humidity_scaling  
  !!
  !! Calculates scaling factors fact_qair and fact_qsat_srf 
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE update_humidity_scaling(nc, &
                                     fract_snow, fract_water, fract_fpc, &
                                     q_sat_srf, q_air, &
                                     wind_air, pch, &
                                     canopy_cond, rel_hum, &
                                     fact_q_air, fact_qsat_srf, fact_qsat_trans_srf)  ! out

    USE mo_jsb_math_constants,       ONLY: eps8
    
    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut 
    INTEGER, INTENT(in)                  :: nc                     !< dimensions
    REAL(wp), DIMENSION(nc), INTENT(in)  :: fract_snow,  &         !< surface fraction covered with snow 
                                            fract_water, &         !< surface fraction covered with water
                                            fract_fpc,   &         !< surface fraction covered with vegetation
                                            q_sat_srf,   &         !< surface saturated specific humidity
                                            q_air,   &             !< near-surface air specific humidity
                                            wind_air, &            !< wind speed [m/2]
                                            pch, &                 !<
                                            canopy_cond, &         !< canopy conductance [m/2]
                                            rel_hum                !< surface relative humidity
    REAL(wp), DIMENSION(nc), INTENT(out) :: fact_q_air, &
                                            fact_qsat_srf, &
                                            fact_qsat_trans_srf
                     
    REAL(wp), DIMENSION(nc) :: fact_qsat_veg,fact_qsat_soil,fact_qsat_trans
    REAL(wp), DIMENSION(nc) :: fact_qair_veg,fact_qair_soil

    !>1.0 vegetation
    !>    
    WHERE(canopy_cond(:) > eps8 .AND.(q_sat_srf(:) - q_air(:)) > 0.0_wp)
       fact_qsat_veg(:) =  fract_snow(:) + (1._wp - fract_snow(:)) *   &
          (fract_water(:) + ( 1._wp - fract_water(:)) /              &
          (1._wp + pch(:) * wind_air(:) / canopy_cond(:)))
       fact_qsat_trans(:) = (1._wp - fract_snow(:)) * &
          (1._wp - fract_water(:)) /                  &
          (1._wp + pch(:) * wind_air(:) / canopy_cond(:))
    ELSEWHERE
       fact_qsat_veg(:)   = fract_snow(:) + (1._wp - fract_snow(:)) * fract_water(:)
       fact_qsat_trans(:) = 0._wp 
    ENDWHERE
    fact_qair_veg(:) = fact_qsat_veg(:)

    !>2.0 soil evaporation
    !>
    WHERE (q_sat_srf(:) * rel_hum(:) > q_air(:) .AND. rel_hum(:) > eps8) 
       fact_qsat_soil(:) = fract_snow(:) + (1._wp - fract_snow(:)) * &
         & (fract_water(:) + (1._wp - fract_water(:)) * rel_hum(:) ) 
       fact_qair_soil(:) = 1._wp
    ELSEWHERE
       fact_qsat_soil(:) = fract_snow(:) + (1._wp - fract_snow(:)) * fract_water(:)
       fact_qair_soil(:) = fact_qsat_soil(:) 
    END WHERE
    WHERE (q_air(:) > q_sat_srf(:))
       fact_qsat_soil(:) = 1._wp
       fact_qair_soil(:) = 1._wp
    END WHERE

    !>3.0 totals surface 
    !>
    fact_qsat_srf(:)        = fract_fpc(:) * fact_qsat_veg(:) + (1._wp - fract_fpc(:)) * fact_qsat_soil(:)
    fact_q_air(:)           = fract_fpc(:) * fact_qair_veg(:) + (1._wp - fract_fpc(:)) * fact_qair_soil(:)
    fact_qsat_trans_srf(:)  = fract_fpc(:) * fact_qsat_trans(:)

  END SUBROUTINE update_humidity_scaling 

  !-----------------------------------------------------------------------------------------------------
  ! Sub Task to calc_surface_energy_balance
  !
  !-----------------------------------------------------------------------------------------------------
  !> Routine to calculate new surface dry static energy using JSBACHs implicit formulation 
  !!  mimicks the routine JSBACH4 surface_temp_implicit 
  !!   srf_energy_bal/mo_seb_land:update_surface_energy_land:surface_temp_implicit
  !!
  !-----------------------------------------------------------------------------------------------------
  ELEMENTAL FUNCTION calc_temp_srf_new_implicit(steplen,t2s_conv,t_Acoef,t_Bcoef,q_Acoef,q_Bcoef, & 
                                                  s_old,qsat_srf_old,dQdT, & 
                                                  net_radiation,ground_heat_flx,heat_transfer_coef, & 
                                                  cair,csat,snow_cover_fract,heat_capa_srf) RESULT (s_new)

    USE mo_jsb_physical_constants,    ONLY: LatentHeatVaporization, LatentHeatSublimation, cvdifts, stbo, zemiss_def

    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 In
    REAL(wp), INTENT(in) :: steplen, & 
                            t2s_conv, & 
                            t_Acoef, & 
                            t_Bcoef, & 
                            q_Acoef, & 
                            q_Bcoef, & 
                            s_old, &
                            qsat_srf_old, &
                            dQdT, &
                            net_radiation, &
                            ground_heat_flx, & 
                            heat_transfer_coef, & 
                            cair, & 
                            csat, & 
                            snow_cover_fract, & 
                            heat_capa_srf
     ! out
     REAL(wp)            :: s_new 
     ! local
     REAL(wp) :: pdt, zicp, zca, zcs, zcolin, zcohfl, zcoind

     pdt = cvdifts * steplen

     zicp = 1._wp / t2s_conv 
 
     zca = LatentHeatSublimation * snow_cover_fract + LatentHeatVaporization * (cair - snow_cover_fract)
     zcs = LatentHeatSublimation * snow_cover_fract + LatentHeatVaporization * (csat - snow_cover_fract)

     zcolin = heat_capa_srf * zicp + & 
              pdt * (zicp * 4.0_wp * stbo * zemiss_def * ((zicp * s_old) ** 3.0_wp) - &
              heat_transfer_coef * ( zca * q_Acoef - zcs) * zicp * dQdT)

     zcohfl = -pdt * heat_transfer_coef * (t_Acoef - 1.0_wp)

     zcoind = pdt * (net_radiation + heat_transfer_coef * t_Bcoef + heat_transfer_coef * &
                 ((zca * q_Acoef - zcs) * qsat_srf_old + zca * q_Bcoef) + ground_heat_flx)

     s_new  = (zcolin * s_old + zcoind) / (zcolin + zcohfl) 

  END FUNCTION calc_temp_srf_new_implicit

#endif
END MODULE mo_spq_update_physics
