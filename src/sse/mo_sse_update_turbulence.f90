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
!! @brief Contains the routines for the Task update_soil_turbulence
!!
!! @author
!!  Sönke Zaehle, MPI-BGC Jena 
!!
!! @par Revision History
!! First version by Sönke Zaehle (20171019): copyied from JSBACH3, mo_jsbalone.f90 
!!
MODULE mo_sse_update_turbulence

#include "dsl_jsbach.inc"

  USE mo_kind,                  ONLY: wp
  USE mo_jsb_control,           ONLY: debug_on
  USE mo_exception,             ONLY: message

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: update_soil_turbulence, calc_filtered_drag_coefficient

  CHARACTER(len=*), PARAMETER :: modname = 'mo_sse_update_turbulence'

CONTAINS


  !-----------------------------------------------------------------------------------------------------
  ! Main Task
  ! 
  !-----------------------------------------------------------------------------------------------------
  !> update soil moisture and theta for the plant only model
  !! 
  !! two cases considered: \n a) dynamic water uptake calculation, \n b) prescribed values from namelist
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE update_soil_turbulence(tile, options)

    USE mo_jsb_class,             ONLY: Get_model
    USE mo_jsb_process_class,     ONLY: A2L_, L2A_, SEB_, TURB_, SSE_, HYDRO_, VEG_, HD_, RAD_, ASSIMI_, PHENO_, SB_
    USE mo_jsb_tile_class,        ONLY: t_jsb_tile_abstract
    USE mo_jsb_task_class,        ONLY: t_jsb_task_options
    USE mo_jsb_model_class,       ONLY: t_jsb_model
    USE mo_jsb_lctlib_class,      ONLY: t_lctlib_element    
    USE mo_jsb_math_constants,    ONLY: dtime
    USE mo_atmland_constants,     ONLY: min_wind, max_wind 
    USE mo_atmland_util,          ONLY: calc_spec_humidity_sat
    USE mo_jsb_phys_constants,    ONLY: grav, r_gas_vapour, r_gas_dryair, cpd, cpvd1, ckap, cvdifts, zkappa, vtmpc1
    USE mo_sse_constants,         ONLY: height_wind

    ! Use of process configurations (t_PROC_config)
    dsl4jsb_Use_config(SSE_)

    ! Use of process memories
    dsl4jsb_Use_memory(A2L_) 
    dsl4jsb_Use_memory(SSE_)        !  USE mo_sse_memory_class,     ONLY: t_sse_memory
    dsl4jsb_Use_memory(VEG_)        !  USE mo_veg_memory_class,      ONLY: t_veg_memory

    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    CLASS(t_jsb_tile_abstract), INTENT(inout)     :: tile         !< one tile with data structure for one lct
    TYPE(t_jsb_task_options),   INTENT(in)        :: options      !< model options
    ! ---------------------------
    ! 0.2 Local
    TYPE(t_jsb_model),      POINTER       :: model
    TYPE(t_lctlib_element), POINTER       :: lctlib     !< land-cover-type library - parameter across pft's
    REAL(wp), DIMENSION(options%nc)       :: ztvd, ztvir, q_sat_surf, ztvs, zg, zgh, zdu2, zril, zcdnl, zchln, & 
                                             zscfl, zucfl, zucfhl, zcons, zcfncl, zcfnchl, zcfml, zcfhl, &
                                             csat, cair, z0m, z0h, zchnl, zlev, wind_air
    REAL(wp)                              :: zcons8, zcons9, zcons11, zcons12
    INTEGER                               :: iblk, ics, ice, nc
    REAL(wp)                              :: steplen
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':update_soil_turbulence'
    ! ---------------------------
    ! 0.3 Declare Memory
    ! Declare process configuration and memory Pointers
    dsl4jsb_Def_config(SSE_)
    dsl4jsb_Def_memory(A2L_) 
    dsl4jsb_Def_memory(SSE_)                        !     TYPE(t_sse_memory),     POINTER :: sse__mem
    dsl4jsb_Def_memory(VEG_)                        !     TYPE(t_veg_memory),      POINTER :: veg__mem
    ! Declare pointers to variables in memory
    ! A2L_
    dsl4jsb_Real2D_onChunk                :: t_air
    dsl4jsb_Real2D_onChunk                :: q_air
    dsl4jsb_Real2D_onChunk                :: wind_10m 
    dsl4jsb_Real2D_onChunk                :: press_srf
    dsl4jsb_Real2D_onChunk                :: drag_srf
    dsl4jsb_Real2D_onChunk                :: ga 
    dsl4jsb_Real2D_onChunk                :: pch 
    dsl4jsb_Real2D_onChunk                :: t_acoef 
    dsl4jsb_Real2D_onChunk                :: t_bcoef
    dsl4jsb_Real2D_onChunk                :: q_acoef
    dsl4jsb_Real2D_onChunk                :: q_bcoef
    ! SSE_
    dsl4jsb_Real3D_onChunk                :: t_soil_sl
    ! VEG_
    dsl4jsb_Real2D_onChunk                :: height 
    dsl4jsb_Real2D_onChunk                :: lai 
    ! Get local variables from options argument
    iblk    = options%iblk
    ics     = options%ics
    ice     = options%ice
    nc      = options%nc
    steplen = options%steplen
    ! ---------------------------
    ! 0.4 Process Activity, Debug Option
    IF (.NOT. tile%Is_process_active(SSE_)) RETURN
    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')
    ! ---------------------------
    ! 0.5 Get Memory
    model  => Get_model(tile%owner_model_id)
    lctlib => model%lctlib(tile%lcts%lib_id)
    ! Get process config
    dsl4jsb_Get_config(SSE_)
    ! Get process memories
    dsl4jsb_Get_memory(A2L_) 
    dsl4jsb_Get_memory(SSE_)
    dsl4jsb_Get_memory(VEG_)
    ! Get process variables (Set pointers to variables in memory)
    dsl4jsb_Get_var2D_onChunk(A2L_, t_air)                ! in
    dsl4jsb_Get_var2D_onChunk(A2L_, q_air)                ! in
    dsl4jsb_Get_var2D_onChunk(A2L_, wind_10m)             ! in
    dsl4jsb_Get_var2D_onChunk(A2L_, press_srf)            ! in
    dsl4jsb_Get_var2D_onChunk(A2L_, drag_srf)             ! out
    dsl4jsb_Get_var2D_onChunk(A2L_, ga)                   ! out
    dsl4jsb_Get_var2D_onChunk(A2L_, pch)                  ! out
    dsl4jsb_Get_var2D_onChunk(A2L_, t_acoef)              ! out
    dsl4jsb_Get_var2D_onChunk(A2L_, t_bcoef)              ! out
    dsl4jsb_Get_var2D_onChunk(A2L_, q_acoef)              ! out
    dsl4jsb_Get_var2D_onChunk(A2L_, q_bcoef)              ! out
    ! ---------------------------
    dsl4jsb_Get_var3D_onChunk(SSE_, t_soil_sl)            ! in
    ! ---------------------------
    dsl4jsb_Get_var2D_onChunk(VEG_, height)               ! in



    ! ------------------------------------------------------------------------------------------------------------
    ! Go science
    ! ------------------------------------------------------------------------------------------------------------


    ! 0.9 init local variables and constants
    ! Some constants from echam/mo_surface_land.f90
    csat(:)     = 0.5_wp
    cair(:)     = 0.5_wp
    wind_air(:) = MIN(MAX(wind_10m(:),min_wind),max_wind) 
    zlev(:)     = height_wind + 3._wp/10._wp * height(:)
    z0m(:)      = MAX(height(:) / 10.0_wp,0.1_wp)
    z0h(:)      = z0m(:) / 10.0_wp

    ! geopotential of the surface layer (see echam's auxhybc.f90 & geopot.f90)
    ! If height_wind is set, then the measurement height + an offset from the average vegetation height
    ! is used. Otherwise, the code defaults to the half-level of ECHAM's lowest atmospheric layer
    zg(:)       = zlev(:) * grav 
    zgh(:)      = zg(:)                          ! height_humidity and height_temperature have to be equal
    zdu2(:)     = wind_air(:)**2._wp

    ! virtual potential air temperature (see mo_surface_boundary.f90)
    ! according to Saucier, WJ Principles of Meteoroligical Analyses
    ! tv = t * (1 + 0.61 * q * t)    ! virtual temperature
    ! td = t * ( 1000 / p_mb ) ^ R/cdp  ! potential temperature
    ! tvd = tair * (100000/p_pa)^zkappa * 1 + vtmpc1 * q_air) ! virtual potential temperature
    ztvd(:)       = t_air(:) * ( 100000._wp/press_srf(:))**zkappa * ( 1._wp + vtmpc1 * q_air(:) )
    ztvir(:)      = t_air(:) * ( 1._wp + vtmpc1 * q_air(:) )
    ! virtual potential surface temperature
    q_sat_surf(:) = calc_spec_humidity_sat(t_soil_sl(:,1),press_srf(:))
    ztvs(:)       = t_soil_sl(:,1) * ( 100000._wp/press_srf(:))**zkappa * &
                    ( 1._wp + vtmpc1 * ( csat(:) * q_sat_surf(:) + ( 1._wp - cair(:) ) * q_air(:)))

    ! Richardson number (dry, Brinkop & Roeckner 1995, Tellus)
    ! ztvd, ztvs are now virtual potential temperatures, changed by Thomas Raddatz 07.2014
    zril(:)       = zg(:) * ( ztvd(:) - ztvs(:) ) / ( zdu2(:) * (ztvd(:)+ztvs(:))/2._wp )

    CALL calc_drag_coefficient(steplen, &
                               t_air(:),q_air(:),wind_air(:),press_srf(:), & 
                               zg(:),zgh(:),z0m(:),z0h(:), &
                               zril(:), &                 ! inout
                               drag_srf(:),pch(:),ga(:))  ! out

    ! Richtmeyer-Morton coefficient for explicit coupling
    t_acoef(:) = 0.0_wp
    t_bcoef(:) = t_air(:) * cpd * (1._wp + cpvd1 * q_air(:)) + zgh 
    q_acoef(:) = 0.0_wp
    q_bcoef(:) = q_air(:)

  END SUBROUTINE update_soil_turbulence


  !-----------------------------------------------------------------------------------------------------
  ! Sub Task to mo_sse_update_physics:calc_surface_energy_balance
  !
  !-----------------------------------------------------------------------------------------------------
  !> Routine to provide an new estimate of the drag coefficent, given the old (t-1), current (t) and \n
  !! new (t+1) surface temperature and other surface conditions for t 
  !!
  !! Calculates time-filtered drag coefficient 
  !-----------------------------------------------------------------------------------------------------
  ELEMENTAL SUBROUTINE calc_filtered_drag_coefficient( &
                                            steplen, &
                                            t_air,q_air,wind_air,press_srf,t_srf,t_srf_new,&                         
                                            height,zril_old,drag_srf,pch,ga)

    USE mo_atmland_util,          ONLY: calc_spec_humidity_sat
    USE mo_jsb_phys_constants,    ONLY: grav, r_gas_vapour, r_gas_dryair, cpd, cpvd1, ckap, cvdifts, zkappa, vtmpc1
    USE mo_sse_constants,         ONLY: height_wind

    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    REAL(wp),                        INTENT(in)    :: steplen                !< step length for this process
    REAL(wp),                        INTENT(in)    :: t_air, &
                                                      q_air, &
                                                      wind_air, &
                                                      press_srf, &
                                                      t_srf, &
                                                      t_srf_new, &
                                                      height
    REAL(wp),                        INTENT(inout) :: zril_old 
    REAL(wp),                        INTENT(out)   :: drag_srf, &
                                                      pch, &
                                                      ga
    ! ---------------------------
    ! 0.2 Local
    REAL(wp)                              :: ztvd, ztvir, q_sat_surf, ztvs, ztvs_new, zg, zgh, zdu2, zril, &
                                             z0m, z0h, csat, cair, zchnl, zlev, zril_initial, zril_update
    REAL(wp)                              :: zcons8, zcons9, zcons11, zcons12
    REAL(wp)                              :: coef_tm1, coef_t, coef_tp1
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_filtered_drag_coefficient'


    ! ------------------------------------------------------------------------------------------------------------
    ! Go science
    ! ------------------------------------------------------------------------------------------------------------


    ! 0.9 init local variables and constants
    ! Some constants from echam/mo_surface_land.f90
    coef_tm1    = 0.5_wp
    coef_t      = 0.25_wp
    coef_tp1    = 0.25_wp
    csat        = 0.5_wp
    cair        = 0.5_wp
    zlev        = height_wind + 3._wp/10._wp * height
    z0m         = MAX(height / 10.0_wp,0.1_wp)
    z0h         = z0m / 10.0_wp

    ! geopotential of the surface layer (see echam's auxhybc.f90 & geopot.f90)
    ! If height_wind is set, then the measurement height + an offset from the average vegetation height
    ! is used. Otherwise, the code defaults to the half-level of ECHAM's lowest atmospheric layer
    zg   = zlev * grav 
    zgh  = zg                          ! height_humidity and height_temperature have to be equal
    zdu2 = wind_air**2._wp

    ! virtual potential air temperature (see mo_surface_boundary.f90)
    ! virtual potential air temperature (see mo_surface_boundary.f90)
    ! according to Saucier, WJ Principles of Meteoroligical Analyses
    ! tv = t * (1 + 0.61 * q * t)    ! virtual temperature
    ! td = t * ( 1000 / p_mb ) ^ R/cdp  ! potential temperature
    ! tvd = tair * (100000/p_pa)^zkappa * 1 + vtmpc1 * q_air) ! virtual potential temperature
    ztvd       = t_air * ( 100000._wp/press_srf)**zkappa * ( 1._wp + vtmpc1 * q_air )
    ztvir      = t_air * ( 1._wp + vtmpc1 * q_air )
    ! virtual potential surface temperature
    q_sat_surf = calc_spec_humidity_sat(t_srf,press_srf)
    ztvs       = t_srf * ( 100000._wp/press_srf)**zkappa * &
                 ( 1._wp + vtmpc1 * ( csat * q_sat_surf + ( 1._wp - cair ) * q_air))
    ! virtual potential surface temperature
    q_sat_surf = calc_spec_humidity_sat(t_srf_new,press_srf)
    ztvs_new   = t_srf_new * ( 100000._wp/press_srf)**zkappa * &
                 ( 1._wp + vtmpc1 * ( csat * q_sat_surf + ( 1._wp - cair ) * q_air))

    ! Richardson number (dry, Brinkop & Roeckner 1995, Tellus)
    ! ztvd, ztvs are now virtual potential temperatures, changed by Thomas Raddatz 07.2014
    zril_initial = MAX(-5._wp,MIN(5._wp,zg * ( ztvd - ztvs )     / ( zdu2 * (ztvd+ztvs)/2._wp )))
    zril_update  = MAX(-5._wp,MIN(5._wp,zg * ( ztvd - ztvs_new ) / ( zdu2 * (ztvd+ztvs_new)/2._wp )))

    ! Richrdson number needs to be filtered for offline simulations in order to maintain stability
    zril = coef_tm1 * zril_old     + &
           coef_t   * zril_initial + &
           coef_tp1 * zril_update

    zril_old = zril   

    CALL calc_drag_coefficient(steplen, &
                               t_air,q_air,wind_air,press_srf, & 
                               zg,zgh,z0m,z0h,zril, & 
                               drag_srf,pch,ga)            

  END SUBROUTINE calc_filtered_drag_coefficient


  !-----------------------------------------------------------------------------------------------------
  ! Sub Task to update_soil_turbulence  &  calc_filtered_drag_coefficient
  !
  !-----------------------------------------------------------------------------------------------------
  !> Routine to calculate drag coefficient for a given surface condition and turbulence 
  !!
  !! Calculates drag coefficient and aerodynamic conductance 
  !-----------------------------------------------------------------------------------------------------
  ELEMENTAL SUBROUTINE calc_drag_coefficient( &
                                   steplen, &
                                   t_air,q_air,wind_air,press_srf, & 
                                   zg,zgh,z0m,z0h,zril, & 
                                   drag_srf,pch,ga)            

    USE mo_phy_schemes,           ONLY: heat_transfer_coef
    USE mo_jsb_phys_constants,    ONLY: grav, r_gas_vapour, r_gas_dryair, ckap, cvdifts, vtmpc1

    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut 
    REAL(wp),                        INTENT(in)    :: steplen                !< step length for this process
    REAL(wp),                        INTENT(in)    :: t_air, &
                                                      q_air, &
                                                      wind_air, &
                                                      press_srf, &
                                                      zg, &
                                                      zgh, &
                                                      z0m, &
                                                      z0h, &
                                                      zril
    REAL(wp),                        INTENT(out)   :: drag_srf, &
                                                      pch, &
                                                      ga
    ! ---------------------------
    ! 0.2 Local                                                      
    REAL(wp)                              :: zcdnl, zchln, & 
                                             zscfl, zucfl, zucfhl, zcons, zcfncl, zcfnchl, zcfml, zcfhl, &
                                             zchnl,zdu2
    REAL(wp)                              :: zcons8,zcons9,zcons11, zcons12
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_drag_coefficient'


    ! ------------------------------------------------------------------------------------------------------------
    ! Go science
    ! ------------------------------------------------------------------------------------------------------------


    ! 0.9 init local variables and constants
    zcons8     = 2._wp * 5.0_wp 
    zcons9     = 3._wp * 5.0_wp 
    zcons11    = 3._wp * 5.0_wp * 5.0_wp
    zcons12    = cvdifts * steplen * grav / r_gas_dryair 
    zdu2       = wind_air**2._wp

    ! Neutral drag coefficient for momentum and heat
    zcdnl = (ckap / LOG(1._wp + zg / (grav *z0m )))**2._wp
    zchnl = ckap**2._wp / (LOG(1._wp + zg / (grav * z0m )) * &
                           LOG( ( grav * z0m + zgh ) / (grav * z0h )) )

    ! account for stable/unstable case: helper variables
    zscfl  = SQRT ( 1._wp + 5._wp * ABS(zril))
    zucfl  = 1._wp / (1._wp + zcons11 * zcdnl * SQRT(ABS(zril) * (1._wp  + zg / (grav * z0m))))
    zucfhl = 1._wp / (1._wp + zcons11 * zchnl * SQRT(ABS(zril) * (1._wp  + zgh / (grav * z0h))))

    ! ignoring cloud water correction (see mo_surface_land.f90)
    zcons    = zcons12 * press_srf / ( t_air * (1._wp + vtmpc1 * q_air))
    zcfncl   = zcons * SQRT(zdu2) * zcdnl
    zcfnchl  = zcons * SQRT(zdu2) * zchnl

    ! Stable / Unstable case
    IF ( zril > 0._wp ) THEN
       zcfml = zcfncl  / (1._wp + zcons8 * zril / zscfl)
       zcfhl = zcfnchl / (1._wp + zcons9 * zril * zscfl)
    ELSE
       zcfml = zcfncl  * (1._wp - zcons8 * zril * zucfl)
       zcfhl = zcfnchl * (1._wp - zcons9 * zril * zucfhl)
    END IF

    ! total surface drag (ECHAM)
    drag_srf = zcfhl
    pch      = zcfhl / zcfnchl * zchnl
    ! aerodynamic conductance (m/s)
    ga       = heat_transfer_coef(drag_srf, steplen) * & 
               wind_air

  END SUBROUTINE calc_drag_coefficient

END MODULE mo_sse_update_turbulence
