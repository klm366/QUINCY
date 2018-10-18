!>
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
  USE mo_slm_varlist,            ONLY: jsb_add_var   ! replace by mo_jsb_varlist later  *slmtemp

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
                              t_srf                 , &  !> surface temperature (K) 
                              t_eff_srf             , &  !> effective surface temperature (K)
                              qsat_srf              , &  !> surface saturated specific humidity (g g-1)
                              fact_q_air            , &  !> factor representing air at actual water content ()
                              fact_qsat_srf         , &  !> factor representing air at saturating water content ()
                              zh_srf                , &  !> roughnesslength for heat transfer (m)
                              zm_srf                , &  !> roughnesslength for momentum transfer (m)
                              alb_vis_dif           , &  !> surface albedo for diffuse radiation in the visible band ()
                              alb_nir_dif           , &  !> surface albedo for diffuse radiation in the near infrared band ()
                              alb_vis_dir           , &  !> surface albedo for direct radiation in the visible band ()
                              alb_nir_dir                !> surface albedo for direct radiation in the near infrared band ()

    ! added from SLM, not existent in jsbach4 20160823
    TYPE(t_jsb_var_real2d) :: &                              
                              albedo                , &  !> total surface albedo ()   
                              sensible_heat_flx     , &  !> sensible heat flux (W m-2)   
                              latent_heat_flx       , &  !> latent heat flux (W m-2)
                              ground_heat_flx       , &  !> ground heat flux (W m-2)   
                              ground_heat_flx_old   , &  !> ground heat flux of previous timestep (W m-2)  
                              net_radiation         , &  !> surface net radiation (W m-2)   
                              vpd                        !> 2m air vapour pressure deficit (Pa)

    ! For lakes
    TYPE(t_jsb_var_real2d) :: &                              
                              latent_hflx_wtr       , &  !> latent heat flux over water (W m-2)
                              sensible_hflx_wtr     , &  !> sensible heat flux over water (W m-2)
                              latent_hflx_ice       , &  !> latent heat flux over ice (W m-2)
                              sensible_hflx_ice          !> sensible heat flux over ice (W m-2)

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

!    USE mo_jsb_varlist,       ONLY: jsb_add_var
    USE mo_jsb_io,            ONLY: grib_bits, t_cf, t_grib1, t_grib2   ! , tables
    USE mo_jsb_grid_class,    ONLY: t_jsb_grid, t_jsb_vgrid
!    USE mo_jsb_grid,          ONLY: Get_grid, Get_vgrid
    USE mo_jsb_model_class,   ONLY: t_jsb_model

    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    CLASS(t_l2a_memory),   INTENT(inout), TARGET :: mem
    CHARACTER(len=*),      INTENT(in)            :: prefix
    CHARACTER(len=*),      INTENT(in)            :: suffix
    INTEGER,               INTENT(in)            :: lct_ids !(:)    !> in jsbach4: vector of land cover type IDs
    INTEGER,               INTENT(in)            :: model_id
    ! ---------------------------
    ! 0.2 Local
    TYPE(t_jsb_model),      POINTER   :: model 
    TYPE(t_jsb_grid),       POINTER   :: hgrid 
    TYPE(t_jsb_vgrid),      POINTER   :: surface  
    INTEGER                           :: table
    CHARACTER(len=*),       PARAMETER :: routine = TRIM(modname)//':Init_a2l_memory'

    model => Get_model(model_id) 

!    table = tables(1)  
    table = 1

!    hgrid => Get_grid(mem%grid_id) 
!    surface => Get_vgrid('surface') 


    !! allocate variables
    !! 
    !! 
    CALL jsb_add_var( mem, 't_srf', mem%t_srf, &
            hgrid, surface, &
            t_cf('t_srf', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp )

    CALL jsb_add_var( mem, 't_eff_srf', mem%t_eff_srf, &
            hgrid, surface, &
            t_cf('t_eff_srf', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp )

    CALL jsb_add_var( mem, 'qsat_srf', mem%qsat_srf, &
            hgrid, surface, &
            t_cf('qsat_srf', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp )

    CALL jsb_add_var( mem, 'fact_q_air', mem%fact_q_air, &
            hgrid, surface, &
            t_cf('fact_q_air', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp )

    CALL jsb_add_var( mem, 'fact_qsat_srf', mem%fact_qsat_srf, &
            hgrid, surface, &
            t_cf('fact_qsat_srf', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp )

    CALL jsb_add_var( mem, 'zh_srf', mem%zh_srf, &
            hgrid, surface, &
            t_cf('zh_srf', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp )

    CALL jsb_add_var( mem, 'zm_srf', mem%zm_srf, &
            hgrid, surface, &
            t_cf('zm_srf', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp )

    CALL jsb_add_var( mem, 'alb_vis_dif', mem%alb_vis_dif, &
            hgrid, surface, &
            t_cf('alb_vis_dif', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp )

    CALL jsb_add_var( mem, 'alb_nir_dif', mem%alb_nir_dif, &
            hgrid, surface, &
            t_cf('alb_nir_dif', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp )

    CALL jsb_add_var( mem, 'alb_vis_dir', mem%alb_vis_dir, &
            hgrid, surface, &
            t_cf('alb_vis_dir', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp )

    CALL jsb_add_var( mem, 'alb_nir_dir', mem%alb_nir_dir, &
            hgrid, surface, &
            t_cf('alb_nir_dir', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp )

    CALL jsb_add_var( mem, 'albedo', mem%albedo, &
            hgrid, surface, &
            t_cf('albedo', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp )

    CALL jsb_add_var( mem, 'sensible_heat_flx', mem%sensible_heat_flx, &
            hgrid, surface, &
            t_cf('sensible_heat_flx', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp )

    CALL jsb_add_var( mem, 'latent_heat_flx', mem%latent_heat_flx, &
            hgrid, surface, &
            t_cf('latent_heat_flx', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp )

    CALL jsb_add_var( mem, 'ground_heat_flx', mem%ground_heat_flx, &
            hgrid, surface, &
            t_cf('ground_heat_flx', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp )

    CALL jsb_add_var( mem, 'ground_heat_flx_old', mem%ground_heat_flx_old, &
            hgrid, surface, &
            t_cf('ground_heat_flx_old', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp )

    CALL jsb_add_var( mem, 'net_radiation', mem%net_radiation, &
            hgrid, surface, &
            t_cf('net_radiation', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp )

    CALL jsb_add_var( mem, 'vpd', mem%vpd, &
            hgrid, surface, &
            t_cf('vpd', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp )

!    IF (One_of(LAKE_TYPE, lct_ids(:)) > 0 .AND. .NOT. ASSOCIATED(mem%latent_hflx_wtr%ptr)) THEN

      CALL jsb_add_var( mem, 'latent_hflx_wtr', mem%latent_hflx_wtr, &
            hgrid, surface, &
            t_cf('latent_hflx_wtr', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp )

      CALL jsb_add_var( mem, 'sensible_hflx_wtr', mem%sensible_hflx_wtr, &
            hgrid, surface, &
            t_cf('sensible_hflx_wtr', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp )

      CALL jsb_add_var( mem, 'latent_hflx_ice', mem%latent_hflx_ice, &
            hgrid, surface, &
            t_cf('latent_hflx_ice', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp )

      CALL jsb_add_var( mem, 'sensible_hflx_ice', mem%sensible_hflx_ice, &
            hgrid, surface, &
            t_cf('sensible_hflx_ice', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp )

!    END IF   ! IF lakes

  END SUBROUTINE Init_l2a_memory

#endif
END MODULE mo_l2a_memory_class
