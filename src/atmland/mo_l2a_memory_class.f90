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
!! @brief Contains types for the atmosphere-land interface memory.
!!
!! @author
!!  Reiner Schnur, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Reiner Schnur (2013-04-04)
!! Integrated into SLM by Jan Engel (20160823) using JSBACH4:git
!!
MODULE mo_l2a_memory_class
#ifndef __NO_JSBACH__

  USE mo_kind,                   ONLY: wp
  USE mo_exception,              ONLY: message
  USE mo_jsb_class,              ONLY: Get_model
  USE mo_jsb_memory_class,       ONLY: t_jsb_memory
  USE mo_jsb_var_class,          ONLY: t_jsb_var_real2d, t_jsb_var_real3d              ! *slmjsb4

! jsbach4 needed for lakes
!  USE mo_util,              ONLY: One_of
!  USE mo_jsb_lct_class,     ONLY: LAKE_TYPE


  IMPLICIT NONE
  PRIVATE
  PUBLIC :: t_l2a_memory, max_no_of_vars

  INTEGER, PARAMETER :: max_no_of_vars = 30


  !-----------------------------------------------------------------------------------------------------
  !> Type definition for l2a memory
  !! 
  !-----------------------------------------------------------------------------------------------------
  TYPE, EXTENDS(t_jsb_memory) :: t_l2a_memory

    TYPE(t_jsb_var_real2d) :: &
                              t_srf                 , &  !> surface temperature [K]
                              t_eff_srf             , &  !> effective surface temperature [K]
                              qsat_srf              , &  !> surface saturated specific humidity [g g-1]
                              fact_q_air            , &  !> factor representing air at actual water content [unitless]
                              fact_qsat_srf         , &  !> factor representing air at saturating water content [unitless]
                              zh_srf                , &  !> roughnesslength for heat transfer [m]
                              zm_srf                , &  !> roughnesslength for momentum transfer [m]
                              alb_vis_dif           , &  !> surface albedo for diffuse radiation in the visible band [unitless]
                              alb_nir_dif           , &  !> surface albedo for diffuse radiation in the near infrared band [unitless]
                              alb_vis_dir           , &  !> surface albedo for direct radiation in the visible band [unitless]
                              alb_nir_dir                !> surface albedo for direct radiation in the near infrared band [unitless]

    ! added from SLM, not existent in jsbach4 20160823
    TYPE(t_jsb_var_real2d) :: &                              
                              sensible_heat_flx     , &  !> sensible heat flux [W m-2]
                              latent_heat_flx       , &  !> latent heat flux [W m-2]
                              ground_heat_flx       , &  !> ground heat flux [W m-2]
                              ground_heat_flx_old        !> ground heat flux of previous timestep [W m-2]

    ! For lakes
    TYPE(t_jsb_var_real2d) :: &                              
                              latent_hflx_wtr       , &  !> latent heat flux over water [W m-2]
                              sensible_hflx_wtr     , &  !> sensible heat flux over water [W m-2]
                              latent_hflx_ice       , &  !> latent heat flux over ice [W m-2]
                              sensible_hflx_ice          !> sensible heat flux over ice [W m-2]

  CONTAINS
    PROCEDURE :: Init => Init_l2a_memory
  END TYPE t_l2a_memory

  CHARACTER(len=*), PARAMETER :: modname = 'mo_l2a_memory_class'

CONTAINS

  !-----------------------------------------------------------------------------------------------------
  !> initialize l2a process memory variables
  !!
  !! 
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE Init_l2a_memory(mem, prefix, suffix, lct_ids, model_id) ! jsbach4

!    USE mo_jsb_varlist,       ONLY: t_jsb_varlist, jsb_add_var
    USE mo_slm_varlist,       ONLY: jsb_add_var   ! replace by mo_jsb_varlist later    *quincytemp
    USE mo_jsb_io,            ONLY: grib_bits, t_cf, t_grib1, t_grib2   ! , TSTEP_CONSTANT, tables
    USE mo_jsb_grid_class,    ONLY: t_jsb_grid, t_jsb_vgrid
    USE mo_jsb_grid,          ONLY: Get_grid, Get_vgrid
    USE mo_jsb_model_class,   ONLY: t_jsb_model

    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    CLASS(t_l2a_memory),  INTENT(inout), TARGET :: mem             !> l2a memory
    CHARACTER(len=*),     INTENT(in)            :: prefix          !> process name
    CHARACTER(len=*),     INTENT(in)            :: suffix          !> tile name
    INTEGER,              INTENT(in)            :: lct_ids !(:)    !> in jsbach4: vector of land cover type IDs
    INTEGER,              INTENT(in)            :: model_id        !> model ID model\%id
    ! ---------------------------
    ! 0.2 Local
    TYPE(t_jsb_model), POINTER :: model
    TYPE(t_jsb_grid),  POINTER :: hgrid                        ! Horizontal grid
    TYPE(t_jsb_vgrid), POINTER :: vgrid_surface                ! Vertical grid
    INTEGER                    :: table
    CHARACTER(len=*),PARAMETER :: routine = TRIM(modname)//':Init_a2l_memory'
    ! ---------------------------
    ! 0.5 Get Memory
    model         => Get_model(model_id)
    table         = 1     ! in jsbach4 this line is: table = tables(1)     *quincytemp
    hgrid         => Get_grid(model%grid_id)
    vgrid_surface => Get_vgrid('surface')


    !> 1.0  init memory
    !! 
    CALL jsb_add_var( mem, 't_srf', mem%t_srf, &
            hgrid, vgrid_surface, &
            t_cf('t_srf', 'K', 'surface temperature'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp )

    CALL jsb_add_var( mem, 't_eff_srf', mem%t_eff_srf, &
            hgrid, vgrid_surface, &
            t_cf('t_eff_srf', 'K', 'effective surface temperature'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp )

    CALL jsb_add_var( mem, 'qsat_srf', mem%qsat_srf, &
            hgrid, vgrid_surface, &
            t_cf('qsat_srf', 'g g-1', 'surface saturated specific humidity'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp )

    CALL jsb_add_var( mem, 'fact_q_air', mem%fact_q_air, &
            hgrid, vgrid_surface, &
            t_cf('fact_q_air', 'unitless', 'factor representing air at actual water content'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp )

    CALL jsb_add_var( mem, 'fact_qsat_srf', mem%fact_qsat_srf, &
            hgrid, vgrid_surface, &
            t_cf('fact_qsat_srf', 'unitless', 'factor representing air at saturating water content'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp )

    CALL jsb_add_var( mem, 'zh_srf', mem%zh_srf, &
            hgrid, vgrid_surface, &
            t_cf('zh_srf', 'm', 'roughnesslength for heat transfer'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp )

    CALL jsb_add_var( mem, 'zm_srf', mem%zm_srf, &
            hgrid, vgrid_surface, &
            t_cf('zm_srf', 'm', 'roughnesslength for momentum transfer'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp )

    CALL jsb_add_var( mem, 'alb_vis_dif', mem%alb_vis_dif, &
            hgrid, vgrid_surface, &
            t_cf('alb_vis_dif', 'unitless', 'surface albedo for diffuse radiation in the visible band'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp )

    CALL jsb_add_var( mem, 'alb_nir_dif', mem%alb_nir_dif, &
            hgrid, vgrid_surface, &
            t_cf('alb_nir_dif', 'unitless', 'surface albedo for diffuse radiation in the near infrared band'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp )

    CALL jsb_add_var( mem, 'alb_vis_dir', mem%alb_vis_dir, &
            hgrid, vgrid_surface, &
            t_cf('alb_vis_dir', 'unitless', 'surface albedo for direct radiation in the visible band'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp )

    CALL jsb_add_var( mem, 'alb_nir_dir', mem%alb_nir_dir, &
            hgrid, vgrid_surface, &
            t_cf('alb_nir_dir', 'unitless', 'surface albedo for direct radiation in the near infrared band'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp )

    CALL jsb_add_var( mem, 'sensible_heat_flx', mem%sensible_heat_flx, &
            hgrid, vgrid_surface, &
            t_cf('sensible_heat_flx', 'W m-2', 'sensible heat flux'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp )

    CALL jsb_add_var( mem, 'latent_heat_flx', mem%latent_heat_flx, &
            hgrid, vgrid_surface, &
            t_cf('latent_heat_flx', 'W m-2', 'latent heat flux'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp )

    CALL jsb_add_var( mem, 'ground_heat_flx', mem%ground_heat_flx, &
            hgrid, vgrid_surface, &
            t_cf('ground_heat_flx', 'W m-2', 'ground heat flux'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp )

    CALL jsb_add_var( mem, 'ground_heat_flx_old', mem%ground_heat_flx_old, &
            hgrid, vgrid_surface, &
            t_cf('ground_heat_flx_old', 'W m-2', 'ground heat flux of previous timestep'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp )

!    IF (One_of(LAKE_TYPE, lct_ids(:)) > 0 .AND. .NOT. ASSOCIATED(mem%latent_hflx_wtr%ptr)) THEN

      CALL jsb_add_var( mem, 'latent_hflx_wtr', mem%latent_hflx_wtr, &
            hgrid, vgrid_surface, &
            t_cf('latent_hflx_wtr', 'W m-2', 'latent heat flux over water'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp )

      CALL jsb_add_var( mem, 'sensible_hflx_wtr', mem%sensible_hflx_wtr, &
            hgrid, vgrid_surface, &
            t_cf('sensible_hflx_wtr', 'W m-2', 'sensible heat flux over water'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp )

      CALL jsb_add_var( mem, 'latent_hflx_ice', mem%latent_hflx_ice, &
            hgrid, vgrid_surface, &
            t_cf('latent_hflx_ice', 'W m-2', 'latent heat flux over ice'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp )

      CALL jsb_add_var( mem, 'sensible_hflx_ice', mem%sensible_hflx_ice, &
            hgrid, vgrid_surface, &
            t_cf('sensible_hflx_ice', 'W m-2', 'sensible heat flux over ice'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp )

!    END IF   ! IF lakes

  END SUBROUTINE Init_l2a_memory

#endif
END MODULE mo_l2a_memory_class
