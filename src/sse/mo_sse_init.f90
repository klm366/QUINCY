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
!! @brief Initialization of the the sse memory.
!!
!! @author
!!  Jan Engel, MPI-BGC Jena
!!
!! @par Revision History
!! First version for SLM by Jan Engel (20161208) based on the JBACH4 version by Reiner Schnur (2013-04-04)
!!
MODULE mo_sse_init
#ifndef __NO_JSBACH__

#include "dsl_jsbach.inc"

  USE mo_kind,                  ONLY: wp
  USE mo_exception,             ONLY: message_text, finish, message
  USE mo_jsb_control,           ONLY: debug_on

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: sse_init

  CHARACTER(len=*), PARAMETER :: modname = 'mo_sse_init'

CONTAINS

  !-----------------------------------------------------------------------------------------------------
  !> Intialize sse process (after memory has been set up)
  !! 
  !! 
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE sse_init(tile, options)

    USE mo_jsb_class,             ONLY: Get_model
    USE mo_jsb_tile_class,        ONLY: t_jsb_tile_abstract
    USE mo_jsb_task_class,        ONLY: t_jsb_task_options
    USE mo_jsb_model_class,       ONLY: t_jsb_model
    USE mo_jsb_lctlib_class,      ONLY: t_lctlib_element
    USE mo_jsb_process_class,     ONLY: A2L_, L2A_, SEB_, TURB_, SSE_, HYDRO_, VEG_, HD_, RAD_, ASSIMI_, PHENO_, SB_
    USE mo_jsb_grid_class,        ONLY: t_jsb_vgrid
    USE mo_jsb_grid,              ONLY: Get_vgrid
    USE mo_util,                  ONLY: soil_depth_to_layers_2d
    USE mo_jsb_math_constants,    ONLY: eps8
    USE mo_sse_constants,         ONLY: soil_heat_cap, soil_therm_cond, kdiff_sat_max, &
                                        k_pwp_s, k_pwp_c, k_pwp_sc, k_pwp_a, k_pwp_at, k_pwp_bt, &
                                        k_fc_s, k_fc_c, k_fc_sc, k_fc_a, k_fc_at, k_fc_bt, k_fc_ct, &
                                        k_sat_s, k_sat_c, k_sat_sc, k_sat_a, k_sat_at, k_sat_bt, k_sat_ct, k_sat_dt, &
                                        wpot_pwp, wpot_fc, wpot_fc

    ! Use of process configurations (t_PROC_config)
    dsl4jsb_Use_config(SSE_)
    dsl4jsb_Use_config(SB_)

    ! USE PROC_memory
    dsl4jsb_Use_memory(SSE_)
    dsl4jsb_Use_memory(VEG_)

    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    CLASS(t_jsb_tile_abstract), INTENT(inout)     :: tile         !< one tile with data structure for one lct
    TYPE(t_jsb_task_options),   INTENT(in)        :: options      !< model options
    ! ---------------------------
    ! 0.2 Local
    TYPE(t_jsb_model),      POINTER                  :: model
    TYPE(t_lctlib_element), POINTER                  :: lctlib        !< land-cover-type library - parameter across pft's
    TYPE(t_jsb_vgrid),      POINTER                  :: vgrid_soil_sb !< Vertical grid
    INTEGER                                          :: nsoil_sb      !< number of soil layers as used/defined by the SB_ process
    REAL(wp), ALLOCATABLE, DIMENSION(:,:)            :: hlp1, &
                                                        theta_pwp_sl, theta_fc_sl, theta_sat_sl
    REAL(wp), ALLOCATABLE, DIMENSION(:,:)            :: sand_sse_init_sl, silt_sse_init_sl, clay_sse_init_sl
    REAL(wp), DIMENSION(options%nc)                  :: soil_awc, &
                                                        hlp2
    REAL(wp)                                         :: site_longitude
    LOGICAL                                          :: flag_spp1685
    CHARACTER(len=3)                                 :: site_ID_spp1685
    INTEGER                                          :: isoil, ichunk
    INTEGER                                          :: iblk, ics, ice, nc
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':sse_init'
    ! ---------------------------
    ! 0.3 Declare Memory
    ! Declare process configuration and memory Pointers
    dsl4jsb_Def_config(SSE_)
    dsl4jsb_Def_config(SB_)
    dsl4jsb_Def_memory(SSE_)
    dsl4jsb_Def_memory(VEG_)
    ! Declare pointers to variables in memory
    ! VEG_
    dsl4jsb_Real3D_onChunk    :: root_fraction_sl
    ! SSE_
    dsl4jsb_Real2D_onChunk    :: soil_depth
    dsl4jsb_Real2D_onChunk    :: root_depth
    dsl4jsb_Real2D_onChunk    :: t_srf_old 
    dsl4jsb_Real2D_onChunk    :: w_skin
    dsl4jsb_Real2D_onChunk    :: w_soil_root
    dsl4jsb_Real2D_onChunk    :: w_soil_root_pwp
    dsl4jsb_Real2D_onChunk    :: w_soil_root_fc
    dsl4jsb_Real2D_onChunk    :: w_soil_root_theta
    dsl4jsb_Real2D_onChunk    :: w_soil_root_pot
    dsl4jsb_Real3D_onChunk    :: soil_depth_sl          ! (soil_depth_sl == soil_lay_width_sl)
    dsl4jsb_Real3D_onChunk    :: soil_lay_width_sl
    dsl4jsb_Real3D_onChunk    :: soil_lay_depth_center_sl
    dsl4jsb_Real3D_onChunk    :: soil_lay_depth_ubound_sl
    dsl4jsb_Real3D_onChunk    :: soil_lay_depth_lbound_sl
    dsl4jsb_Real3D_onChunk    :: bulk_dens_sl
    dsl4jsb_Real3D_onChunk    :: heat_capa_sl
    dsl4jsb_Real3D_onChunk    :: therm_cond_sl
    dsl4jsb_Real3D_onChunk    :: w_soil_pwp_sl
    dsl4jsb_Real3D_onChunk    :: w_soil_fc_sl 
    dsl4jsb_Real3D_onChunk    :: w_soil_sat_sl 
    dsl4jsb_Real3D_onChunk    :: w_soil_sl
    dsl4jsb_Real3D_onChunk    :: w_soil_pot_sl
    dsl4jsb_Real3D_onChunk    :: saxtonA 
    dsl4jsb_Real3D_onChunk    :: saxtonB 
    dsl4jsb_Real3D_onChunk    :: saxtonC 
    dsl4jsb_Real3D_onChunk    :: kdiff_sat_sl 
    dsl4jsb_Real3D_onChunk    :: sand_sl
    dsl4jsb_Real3D_onChunk    :: silt_sl
    dsl4jsb_Real3D_onChunk    :: clay_sl
    dsl4jsb_Real3D_onChunk    :: volume_min_sl
    dsl4jsb_Real3D_onChunk    :: t_soil_sl
    ! Get local variables from options argument
    iblk      = options%iblk
    ics       = options%ics
    ice       = options%ice
    nc        = options%nc
    ! ---------------------------
    ! 0.4 Process Activity, Debug Option
    IF (.NOT. tile%Is_process_active(SSE_)) RETURN
    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')
    ! ---------------------------
    ! 0.5 Get Memory
    model          => Get_model(tile%owner_model_id)
    lctlib         => model%lctlib(tile%lcts%lib_id)
    vgrid_soil_sb  => Get_vgrid('soil_layer_sb')
    nsoil_sb       =  vgrid_soil_sb%n_levels
    site_longitude = model%config%lon
    flag_spp1685   = model%config%flag_spp1685
    ! Get process config
    dsl4jsb_Get_config(SSE_)
    dsl4jsb_Get_config(SB_)
    ! Get process memories
    dsl4jsb_Get_memory(SSE_)
    dsl4jsb_Get_memory(VEG_)
    ! Allocate local allocatable variables (deallocate below at the end of this routine)    
    ALLOCATE(hlp1(nc, nsoil_sb))
    ALLOCATE(theta_pwp_sl(nc, nsoil_sb))
    ALLOCATE(theta_fc_sl(nc, nsoil_sb))
    ALLOCATE(theta_sat_sl(nc, nsoil_sb))
    ALLOCATE(sand_sse_init_sl(nc, nsoil_sb))
    ALLOCATE(clay_sse_init_sl(nc, nsoil_sb))
    ALLOCATE(silt_sse_init_sl(nc, nsoil_sb))
    ! Get process variables (Set pointers to variables in memory)
    dsl4jsb_Get_var3D_onChunk(VEG_,  root_fraction_sl)
    ! ---------------------------
    dsl4jsb_Get_var2D_onChunk(SSE_, soil_depth)
    dsl4jsb_Get_var2D_onChunk(SSE_, root_depth)
    dsl4jsb_Get_var2D_onChunk(SSE_, t_srf_old)
    dsl4jsb_Get_var2D_onChunk(SSE_, w_skin)
    dsl4jsb_Get_var2D_onChunk(SSE_, w_soil_root)
    dsl4jsb_Get_var2D_onChunk(SSE_, w_soil_root_pwp)
    dsl4jsb_Get_var2D_onChunk(SSE_, w_soil_root_fc)
    dsl4jsb_Get_var2D_onChunk(SSE_, w_soil_root_theta)
    dsl4jsb_Get_var2D_onChunk(SSE_, w_soil_root_pot)
    dsl4jsb_Get_var3D_onChunk(SSE_, soil_depth_sl)              ! out   (soil_depth_sl == soil_lay_width_sl)
    dsl4jsb_Get_var3D_onChunk(SSE_, soil_lay_width_sl)          ! out
    dsl4jsb_Get_var3D_onChunk(SSE_, soil_lay_depth_center_sl)   ! out
    dsl4jsb_Get_var3D_onChunk(SSE_, soil_lay_depth_ubound_sl)   ! out
    dsl4jsb_Get_var3D_onChunk(SSE_, soil_lay_depth_lbound_sl)   ! out
    dsl4jsb_Get_var3D_onChunk(SSE_, bulk_dens_sl)
    dsl4jsb_Get_var3D_onChunk(SSE_, heat_capa_sl)
    dsl4jsb_Get_var3D_onChunk(SSE_, therm_cond_sl)
    dsl4jsb_Get_var3D_onChunk(SSE_, w_soil_pwp_sl)
    dsl4jsb_Get_var3D_onChunk(SSE_, w_soil_fc_sl)
    dsl4jsb_Get_var3D_onChunk(SSE_, w_soil_sat_sl)
    dsl4jsb_Get_var3D_onChunk(SSE_, w_soil_sl)
    dsl4jsb_Get_var3D_onChunk(SSE_, w_soil_pot_sl)
    dsl4jsb_Get_var3D_onChunk(SSE_, saxtonA)
    dsl4jsb_Get_var3D_onChunk(SSE_, saxtonB)
    dsl4jsb_Get_var3D_onChunk(SSE_, saxtonC)
    dsl4jsb_Get_var3D_onChunk(SSE_, kdiff_sat_sl)
    dsl4jsb_Get_var3D_onChunk(SSE_, sand_sl)
    dsl4jsb_Get_var3D_onChunk(SSE_, silt_sl)
    dsl4jsb_Get_var3D_onChunk(SSE_, clay_sl)
    dsl4jsb_Get_var3D_onChunk(SSE_, volume_min_sl)
    dsl4jsb_Get_var3D_onChunk(SSE_, t_soil_sl)


    ! ------------------------------------------------------------------------------------------------------------
    ! Go science
    ! ------------------------------------------------------------------------------------------------------------

    !> 1.0 Total soil depth until bedrock
    ! QUINCY soil_depth for one site (read from namelist)
    soil_depth(:)                             = dsl4jsb_Config(SSE_)%soil_depth
    ! code from jsbach4 (hydro init)
    ! for all points of the domain soil_depth(:,:) is at least 0.1
    dsl4jsb_var2D_onDomain(SSE_, soil_depth)  = MAX(0.1_wp, &
                                                    MERGE(dsl4jsb_var2D_onDomain(SSE_, soil_depth), 0._wp, &
                                                          dsl4jsb_var2D_onDomain(SSE_, soil_depth) >= 0._wp))

    !> 2.0 calc actual soil layer depths based on the fixed layer thicknesses from namelist and the bedrock depth
    !!  note, soil_depth_to_layers_2d() returns an actual 3D variable
    dsl4jsb_var3D_onDomain(SSE_, soil_depth_sl) = soil_depth_to_layers_2d( &
                                                    dsl4jsb_var2D_onDomain(SSE_, soil_depth), &   ! Total soil depth until bedrock (from textures), 2D variable
                                                    vgrid_soil_sb%dz(:))                          ! Soil layer thicknesses from namelist
    ! pass these values to a SSE_ variable with a correct name
    dsl4jsb_var3D_onDomain(SSE_, soil_lay_width_sl) = dsl4jsb_var3D_onDomain(SSE_, soil_depth_sl)

    !> 2.1 calc more metrics of the soil layers
    !! i)  lower & upper boundary
    !! ii) depth at the center of each layer
    ! @TODO sse_init must init variables onDomain for jsbach4, write a function similar to soil_depth_to_layers_2d()
    DO isoil = 1,nsoil_sb
      ! lower & upper boundary
      IF (isoil == 1) THEN
        soil_lay_depth_lbound_sl(:,isoil) = 0.0_wp
        soil_lay_depth_ubound_sl(:,isoil) = soil_lay_width_sl(:,isoil)
      ELSE
        soil_lay_depth_lbound_sl(:,isoil) = soil_lay_depth_lbound_sl(:,isoil-1) + soil_lay_width_sl(:,isoil-1)
        soil_lay_depth_ubound_sl(:,isoil) = soil_lay_depth_ubound_sl(:,isoil-1) + soil_lay_width_sl(:,isoil)
      ENDIF                                
      ! soil-layer center
      soil_lay_depth_center_sl(:,isoil)   = (soil_lay_depth_lbound_sl(:,isoil) + soil_lay_depth_ubound_sl(:,isoil)) * 0.5_wp 
    ENDDO

    !> 3.0 check if last soil layer has a thickness larger zero
    !! QUINCY specific code: for jsb4 it would need to work on the domain
    DO ichunk  = 1,nc
      IF (soil_depth_sl(ichunk,nsoil_sb) < eps8) THEN
        WRITE(message_text,*) 'The last soil layer has a thickness < eps8. Increase soil_depth from ', soil_depth, &
                              ' to at least ', SUM(soil_lay_width_sl(ics,1:nsoil_sb-1)), ' but better ', &
                              SUM(soil_lay_width_sl(ics,1:nsoil_sb)), ' and re-start the simulation.'
        CALL finish(TRIM(routine), message_text)
      END IF  
    END DO
    
    !> 4.0 init soil-property variables
    !! 
    
    !> 4.1 init with one site-specific value for all soil layers
    !! 
    bulk_dens_sl(:,:)   = dsl4jsb_Config(SSE_)%bulk_density 
    heat_capa_sl(:,:)   = soil_heat_cap
    therm_cond_sl(:,:)  = soil_therm_cond
    sand_sl(:,:)        = dsl4jsb_Config(SSE_)%soil_sand
    silt_sl(:,:)        = dsl4jsb_Config(SSE_)%soil_silt
    clay_sl(:,:)        = dsl4jsb_Config(SSE_)%soil_clay
    t_srf_old(:)        = 273.15_wp
    volume_min_sl(:,:)  = 1._wp
    t_soil_sl(:,:)      = 273.15_wp

    !------------------------------------------------------------------------------------------------------------------------------
    ! SPP site-set
    !
    !> 4.1.1 init with site-specific and soil-layer specific texture (sand silt clay) information
    !! 
    !! works only with the SPP_1685 site-set (namelist flag flag_spp1685 is set depending on the site-set); 
    !! met-forcing was created with the usual scripts; 
    !! soil forcing for JSM was created using QUINCY land (current trunk) with simple_1d and this texture information; 
    !! hence, JSM can run simulations for the SPP_1685 sites must use the below soil-layer specific texture init
    !! 
    !! the below texture data is based on: Lang, F. et al. (2017 ) "Soil phosphorus supply controls P nutrition strategies 
    !!   of beech forest ecosystems in Central Europe". Biogeochemistry, doi:10.1007/s10533-017-0375-0. Table S2 in the
    !!   Supplementary material 2.
    !!   These texture data were extrapolated: We use the measured soil texture of the topmost 1m soil for all the five sites 
    !!   (presented in this table S2). We extrapolated them to the whole soil depth (9.5m) following the pattern of 
    !!   texture change against depth increase observed in the topmost 1m soil.
    !! 
    IF (flag_spp1685) THEN
      IF (dsl4jsb_Config(SB_)%sb_model_scheme == "simple_1d" .OR. dsl4jsb_Config(SB_)%sb_model_scheme == "jsm") THEN
        ! check if nsoil_sb == 15
        IF (nsoil_sb /= 15) THEN
          CALL finish(TRIM(routine),'nsoil_sb is not 15 | Simulations with the SPP_1685 site-set work only with 15 soil layers.')
        ENDIF
        
        ! get site ID from longitude
        IF (ABS(site_longitude - 9.75_wp) < eps8) THEN
          site_ID_spp1685 = "BBR"
        ELSEIF (ABS(site_longitude - 10.25_wp) < eps8) THEN
          site_ID_spp1685 = "LUE"
        ELSEIF (ABS(site_longitude - 10.75_wp) < eps8) THEN
          site_ID_spp1685 = "VES"
        ELSEIF (ABS(site_longitude - 12.75_wp) < eps8) THEN
          site_ID_spp1685 = "MIT"
        ELSEIF (ABS(site_longitude - 7.75_wp) < eps8) THEN
          site_ID_spp1685 = "CON"
        ELSE
          site_ID_spp1685 = "ERR"
        ENDIF

        ! set site specific soil texture values per soil layer
        SELECT CASE(site_ID_spp1685)
        CASE ("BBR")
          DO ichunk=1,nc
            bulk_dens_sl(ichunk,1:nsoil_sb)   = (/1000._wp,1080._wp,1200._wp,1320._wp,1450._wp, &
                                                 1520._wp,1600._wp,1650._wp,1700._wp,1750._wp, &
                                                 1800._wp, 1850._wp, 1900._wp, 1950._wp, 2000._wp/)
            sand_sl(ichunk,1:nsoil_sb)   = (/0.07_wp,0.09_wp,0.1_wp,0.1228_wp,0.1785_wp, &
                                            0.24_wp,0.2726_wp,0.3519_wp,0.4_wp,0.45_wp, &
                                            0.5_wp,0.6_wp,0.7_wp,0.8_wp,0.8_wp/)
            silt_sl(ichunk,1:nsoil_sb)   = (/0.56_wp,0.54_wp,0.55_wp,0.5386_wp,0.5184_wp, &
                                            0.5_wp,0.5042_wp,0.4868_wp,0.4_wp,0.37_wp, &
                                            0.35_wp,0.3_wp,0.22_wp,0.15_wp,0.15_wp/)
          ENDDO
        CASE ("LUE")
          DO ichunk=1,nc
            bulk_dens_sl(ichunk,1:nsoil_sb)   = (/1570._wp,1600._wp,1610._wp,1630._wp,1660._wp, &
                                                 1800._wp,2000._wp,2000._wp,2200._wp,2200._wp, &
                                                 2200._wp, 2200._wp, 2200._wp, 2200._wp, 2200._wp/)
            sand_sl(ichunk,1:nsoil_sb)   = (/0.75_wp,0.76_wp,0.81_wp,0.7543_wp,0.7448_wp, &
                                             0.713_wp,0.7626_wp,0.8185_wp,0.84_wp,0.85_wp, &
                                             0.85_wp,0.88_wp,0.88_wp,0.9_wp,0.92_wp/)
            silt_sl(ichunk,1:nsoil_sb)   = (/0.19_wp,0.19_wp,0.14_wp,0.1864_wp,0.1926_wp, &
                                             0.2085_wp,0.1753_wp,0.1177_wp,0.1_wp,0.1_wp, &
                                             0.1_wp,0.08_wp,0.08_wp,0.06_wp,0.05_wp/)
          ENDDO    
        CASE ("VES")
          DO ichunk=1,nc
            bulk_dens_sl(ichunk,1:nsoil_sb)   = (/1000._wp,1000._wp,1000._wp,1000._wp,1050._wp, &
                                                 1050._wp,1050._wp,1200._wp,1250._wp,1300._wp, &
                                                 1350._wp, 1400._wp, 1450._wp, 1500._wp, 1600._wp/)
            sand_sl(ichunk,1:nsoil_sb)   = (/0.30_wp,0.3245_wp,0.3746_wp,0.4057_wp,0.4494_wp, &
                                             0.4923_wp,0.5627_wp,0.6669_wp,0.71_wp,0.75_wp, &
                                             0.77_wp,0.79_wp,0.83_wp,0.86_wp,0.88_wp/)
            silt_sl(ichunk,1:nsoil_sb)   = (/0.4638_wp,0.47_wp,0.4562_wp,0.4579_wp,0.4386_wp, &
                                             0.4012_wp,0.3578_wp,0.2864_wp,0.25_wp,0.22_wp, &
                                             0.2_wp,0.18_wp,0.15_wp,0.12_wp,0.1_wp/)
          ENDDO
        CASE ("MIT")
          DO ichunk=1,nc
            bulk_dens_sl(ichunk,1:nsoil_sb)   = (/1000._wp,1000._wp,1000._wp,1000._wp,1050._wp, &
                                                 1050._wp,1100._wp,1100._wp,1400._wp,1450._wp, &
                                                 1450._wp, 1500._wp, 1550._wp, 1600._wp, 1650._wp/)
            sand_sl(ichunk,1:nsoil_sb)   = (/0.447_wp,0.5115_wp,0.5102_wp,0.52_wp,0.5642_wp, &
                                             0.57_wp,0.56_wp,0.5681_wp,0.60_wp,0.62_wp, &
                                             0.64_wp,0.65_wp,0.70_wp,0.70_wp,0.70_wp/)
            silt_sl(ichunk,1:nsoil_sb)   = (/0.3223_wp,0.3115_wp,0.3148_wp,0.32_wp,0.3195_wp, &
                                             0.34_wp,0.3526_wp,0.3601_wp,0.33_wp,0.32_wp, &
                                             0.31_wp,0.30_wp,0.26_wp,0.26_wp,0.27_wp/)
          ENDDO
        CASE ("CON")
          DO ichunk=1,nc
            bulk_dens_sl(ichunk,1:nsoil_sb)   = (/1000._wp,1000._wp,1000._wp,1000._wp,1560._wp, &
                                                 1760._wp,1760._wp,1760._wp,1800._wp,2000._wp, &
                                                 2200._wp, 2200._wp, 2450._wp, 2450._wp, 2450._wp/)
            sand_sl(ichunk,1:nsoil_sb)   = (/0.40_wp,0.3969_wp,0.39_wp,0.40_wp,0.41_wp, &
                                             0.4447_wp,0.4872_wp,0.5159_wp,0.55_wp,0.58_wp, &
                                             0.60_wp,0.62_wp,0.65_wp,0.70_wp,0.70_wp/)
            silt_sl(ichunk,1:nsoil_sb)   = (/0.33_wp,0.3454_wp,0.38_wp,0.40_wp,0.42_wp, &
                                             0.40_wp,0.3992_wp,0.3961_wp,0.37_wp,0.35_wp, &
                                             0.33_wp,0.31_wp,0.29_wp,0.25_wp,0.25_wp/)
          ENDDO
        CASE DEFAULT
          CALL finish(TRIM(routine),'This site seems to be not a part of the SPP_1685 site-set although flag_spp1685 is set T.')
        END SELECT
        
        ! calculate clay content from sand and silt
        clay_sl(:,:) = 1.0_wp
        clay_sl(:,:) = clay_sl(:,:) - (sand_sl(:,:) + silt_sl(:,:))
      ENDIF  ! simple_1d / jsm
    ENDIF   ! flag_spp1685
    !------------------------------------------------------------------------------------------------------------------------------

    !> 4.2 ensure a certain minimum proportion for sand, clay and silt
    !! this is to ensure applying the empirical functions by Saxton within the range of the texture they were designed for \n
    !! avoiding ranges where calculation would return negative theta_pwp_sl or relatively large kdiff_sat_sl values
    !! 
    !! Saxton equations with too small soil-texture fractions may lead to unrealistic percolation rates and soil water potentials
    !!
    !! for more details see ticket #560
    !!  
    !! using the modified *_sse_init_sl soil texture for "saxton functions here in sse_init" only
    sand_sse_init_sl(:,:) = sand_sl(:,:)
    clay_sse_init_sl(:,:) = clay_sl(:,:)
    silt_sse_init_sl(:,:) = silt_sl(:,:)
    CALL calc_soil_texture_saxton_compatible(nc, &
                                             nsoil_sb, &
                                             sand_sse_init_sl(:,:), &
                                             silt_sse_init_sl(:,:), &
                                             clay_sse_init_sl(:,:))
    
    !> 5.0 soil hydraulic properties, based on Saxton and Rawls 2006, SSSA
    !! 

    !> 5.1 fractional water holding capacity at permanent wilting point
    !! note that this assumes that the soil hydraulic properties are identical between the layers
    !!
    hlp1(:,:)         = k_pwp_s * sand_sse_init_sl(:,:) + & 
                        k_pwp_c * clay_sse_init_sl(:,:) + &
                        k_pwp_sc * sand_sse_init_sl(:,:) * clay_sse_init_sl(:,:) + &
                        k_pwp_a
    theta_pwp_sl(:,:) = k_pwp_at + (1.0_wp + k_pwp_bt) * hlp1(:,:) 

    !> 5.2 fractional water holding capacity at field capacity 
    !! 
    hlp1(:,:)         = k_fc_s * sand_sse_init_sl(:,:) + &
                        k_fc_c * clay_sse_init_sl(:,:) + &
                        k_fc_sc * sand_sse_init_sl(:,:) * clay_sse_init_sl(:,:) + &
                        k_fc_a
    theta_fc_sl(:,:)  = k_fc_at + (1.0_wp + k_fc_bt) * hlp1(:,:) + k_fc_ct * hlp1(:,:)**2._wp

    !> 5.3 fractional water holding capacity at saturation 
    !! 
    hlp1(:,:)         = k_sat_s * sand_sse_init_sl(:,:) + & 
                        k_sat_c * clay_sse_init_sl(:,:) + &
                        k_sat_sc * sand_sse_init_sl(:,:) * clay_sse_init_sl(:,:) + &
                        k_sat_a
    theta_sat_sl(:,:) = theta_fc_sl(:,:) + k_sat_at + (1.0_wp + k_sat_bt) * hlp1(:,:) - & 
                        k_sat_ct * sand_sse_init_sl(:,:) + k_sat_dt 
    
    !> 5.4 parameters describing the moisture-tension curve                 
    !! note that for convenience, B has a negative sign and A is in MPa 
    saxtonB(:,:)      = - (log(wpot_pwp) - log(wpot_fc))/(log(theta_fc_sl(:,:)) - log(theta_pwp_sl(:,:)))
    saxtonA(:,:)      = - (exp(log(wpot_fc) - saxtonB(:,:) * log(theta_fc_sl(:,:))))/1000.0_wp

    !> 5.5 parameters describing the moisture-diffusivity curve
    !! 
    kdiff_sat_sl(:,:) = kdiff_sat_max * (theta_sat_sl(:,:) - theta_fc_sl(:,:)) ** (3.0_wp + 1._wp/saxtonB(:,:))
    saxtonC(:,:)      = 3.0_wp - 2._wp * saxtonB(:,:)


    !> 6.0 rooting depth given by prescribed AWC
    !! soil water in the rooting zone [m] \n
    !! the init of 'root_fraction_sl(:,:)' (sec 6.3) depends on this init of the root_depth(:)
    !! 
    soil_awc(:)         = 0._wp
    w_soil_root_pwp(:)  = 0._wp
    w_soil_root_fc(:)   = 0._wp
    root_depth(:)       = 0._wp
    DO isoil=1,nsoil_sb
      hlp1(:,isoil)     = soil_lay_width_sl(:,isoil) * (theta_fc_sl(:,isoil) - theta_pwp_sl(:,isoil))
      soil_awc(:)       = soil_awc(:) + hlp1(:,isoil)
      WHERE (soil_awc(:) < (dsl4jsb_Config(SSE_)%soil_awc_prescribe / 1000._wp)) 
        root_depth(:)       = root_depth(:) + soil_lay_width_sl(:,isoil)
        w_soil_root_pwp(:)  = w_soil_root_pwp(:) + theta_pwp_sl(:,isoil) * soil_lay_width_sl(:,isoil)
        w_soil_root_fc(:)   = w_soil_root_fc(:) + theta_fc_sl(:,isoil) * soil_lay_width_sl(:,isoil)
      ELSEWHERE
        hlp2(:)             = MAX(0._wp, dsl4jsb_Config(SSE_)%soil_awc_prescribe / 1000._wp - soil_awc(:) + hlp1(:,isoil))
        root_depth(:)       = root_depth(:) + hlp2(:) / (theta_fc_sl(:,isoil) - theta_pwp_sl(:,isoil))
        w_soil_root_pwp(:)  = w_soil_root_pwp(:) + theta_pwp_sl(:,isoil) * hlp2(:) / (theta_fc_sl(:,isoil) - theta_pwp_sl(:,isoil))
        w_soil_root_fc(:)   = w_soil_root_fc(:) + theta_fc_sl(:,isoil) * hlp2(:) / (theta_fc_sl(:,isoil) - theta_pwp_sl(:,isoil))
      END WHERE
    END DO

    !> 6.1 water in the surface reservoir [m]
    !! 
    w_skin(:) = 0.0_wp

    !> 6.2 adjust initial water amount to prescribed theta 
    !! 
    w_soil_root(:)        = w_soil_root_pwp(:) + & 
                            dsl4jsb_Config(SSE_)%soil_theta_prescribe * (w_soil_root_fc(:) - w_soil_root_pwp(:)) 
    w_soil_root_theta(:)  = dsl4jsb_Config(SSE_)%soil_theta_prescribe
    w_soil_root_pot(:)    = saxtonA(:,1) * (MIN(w_soil_root_fc,w_soil_root(:)) / root_depth(:)) ** saxtonB(:,1) 

    !> 6.3 root fraction 
    !! 
    !! root fractions assuming an exponential distribution of roots \n
    !!   Cr_l = Cr_0 * exp(-lctlib%k_root_dist*depth) \n
    !! integral to bottom of root zone: 1/lctlib%k_root_dist * (1 - exp(-lctlib%k_root_dist*depth) \n
    !!   thus for each layer
    !!
    !! root_fraction_sl(:,:) is a variable of veg_mem and depends on the init of root_depth(:), see section 6.0 and ticket #358
    !! 
    !! note: upper bound is the larger value compared to lower bound !
    !! 
    ! layer 1 (vgrid_soil_sb)
    root_fraction_sl(:,1) = (1._wp - exp(-lctlib%k_root_dist * soil_lay_depth_ubound_sl(:,1))) / &
                            (1._wp - exp(-lctlib%k_root_dist * root_depth(:)))
    ! layer 2 (vgrid_soil_sb)
    DO isoil=2,nsoil_sb
       WHERE(root_depth(:) > soil_lay_depth_ubound_sl(:,isoil-1))
          WHERE(root_depth(:) > soil_lay_depth_ubound_sl(:,isoil))
             root_fraction_sl(:,isoil) = (exp(-lctlib%k_root_dist * soil_lay_depth_ubound_sl(:,isoil-1)) - &
                                          exp(-lctlib%k_root_dist * soil_lay_depth_ubound_sl(:,isoil))) / & 
                                         (1._wp - exp(-lctlib%k_root_dist * root_depth(:)))
          ELSEWHERE
             root_fraction_sl(:,isoil) = (exp(-lctlib%k_root_dist * soil_lay_depth_ubound_sl(:,isoil-1)) - &
                                          exp(-lctlib%k_root_dist * root_depth(:))) / & 
                                         (1._wp - exp(-lctlib%k_root_dist * root_depth(:)))
          ENDWHERE
       ELSEWHERE
          root_fraction_sl(:,isoil) = 0.0_wp 
       ENDWHERE
    ENDDO


    !> 7.0 soil water per layer depth
    !! 
    w_soil_pwp_sl(:,:)  = 0.0_wp
    w_soil_fc_sl(:,:)   = 0.0_wp
    w_soil_sat_sl(:,:)  = 0.0_wp
    DO isoil=1,nsoil_sb
       WHERE(soil_depth_sl(:,isoil) > 0.0_wp) 
          w_soil_pwp_sl(:,isoil) = theta_pwp_sl(:,isoil) * soil_depth_sl(:,isoil) 
          w_soil_fc_sl(:,isoil)  = theta_fc_sl(:,isoil) * soil_depth_sl(:,isoil) 
          w_soil_sat_sl(:,isoil) = theta_sat_sl(:,isoil) * soil_depth_sl(:,isoil) 
       ENDWHERE
    ENDDO
    
    !> 7.1 adjust initial water amount to prescribed theta
    !! 
    w_soil_sl(:,:) = w_soil_pwp_sl(:,:) + dsl4jsb_Config(SSE_)%soil_theta_prescribe * (w_soil_fc_sl(:,:) - w_soil_pwp_sl(:,:))

    !> 7.2 soil water potential per layer
    !! 
    WHERE(soil_depth_sl(:,:) > eps8)
       hlp1(:,:)          = MIN(w_soil_fc_sl(:,:),w_soil_sl(:,:)) / soil_depth_sl(:,:)
       w_soil_pot_sl(:,:) = saxtonA(:,:) * hlp1(:,:) ** saxtonB(:,:)
    ELSEWHERE
       w_soil_pot_sl(:,:) = 0.0_wp 
    ENDWHERE


  ! Deallocate local allocatable variables
  DEALLOCATE(hlp1)
  DEALLOCATE(theta_pwp_sl)
  DEALLOCATE(theta_fc_sl)
  DEALLOCATE(theta_sat_sl)

  END SUBROUTINE sse_init


  !-----------------------------------------------------------------------------------------------------
  !> ensure a certain minimum proportion for sand, clay and silt
  !! 
  !! Saxton equations with too small soil-texture fractions may lead to unrealistic percolation rates and soil water potentials
  !!
  !! using these modified soil texture for "saxton functions in sse_init" only
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE calc_soil_texture_saxton_compatible(nc, nsoil_sb, sand_sl, silt_sl, clay_sl)

    USE mo_jsb_math_constants,          ONLY: eps8
    USE mo_kind,                        ONLY: wp
    USE mo_exception,                   ONLY: message


    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut 
    INTEGER,                                          INTENT(in)    :: nc, &             !< dimensions
                                                                       nsoil_sb
    REAL(wp), DIMENSION(nc,nsoil_sb),                 INTENT(inout) :: sand_sl, &        !< sand proportion per layer
                                                                       silt_sl, &        !< silt proportion per layer
                                                                       clay_sl           !< clay proportion per layer
    ! ---------------------------
    ! 0.2 Local
    REAL(wp), DIMENSION(nc,nsoil_sb)                :: hlp_sand_sl, &
                                                       hlp_silt_sl, &
                                                       hlp_clay_sl
    REAL(wp), DIMENSION(nc,nsoil_sb)                :: sand_original_sl, &
                                                       silt_original_sl, &
                                                       clay_original_sl
    REAL(wp)                                        :: min_prop_soil_texture, &
                                                       max_prop_soil_texture
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_soil_texture_saxton_compatible'


    ! ------------------------------------------------------------------------------------------------------------
    ! Go science
    ! ------------------------------------------------------------------------------------------------------------


    !> 0.9 init local variables
    !!  
    min_prop_soil_texture = 0.2_wp
    max_prop_soil_texture = 1.0_wp - (2.0_wp * min_prop_soil_texture)
    hlp_sand_sl(:,:)      = 0.0_wp
    hlp_clay_sl(:,:)      = 0.0_wp
    hlp_silt_sl(:,:)      = 0.0_wp
    sand_original_sl(:,:) = sand_sl(:,:)
    clay_original_sl(:,:) = clay_sl(:,:)
    silt_original_sl(:,:) = silt_sl(:,:)

    !> 1.0 ensure a minimum proportion 'min_prop_soil_texture' for sand, clay and silt
    !! this is to ensure applying the empirical functions by Saxton within the range of the texture they were designed for \n
    !! 
    WHERE (sand_sl(:,:) < min_prop_soil_texture .AND. clay_sl(:,:) < min_prop_soil_texture) 
      sand_sl(:,:)     = min_prop_soil_texture
      clay_sl(:,:)     = min_prop_soil_texture
      silt_sl(:,:)     = max_prop_soil_texture
    ELSEWHERE (sand_sl(:,:) < min_prop_soil_texture .AND. silt_sl(:,:) < min_prop_soil_texture)
      sand_sl(:,:)     = min_prop_soil_texture
      silt_sl(:,:)     = min_prop_soil_texture
      clay_sl(:,:)     = max_prop_soil_texture
    ELSEWHERE (clay_sl(:,:) < min_prop_soil_texture .AND. silt_sl(:,:) < min_prop_soil_texture)
      clay_sl(:,:)     = min_prop_soil_texture
      silt_sl(:,:)     = min_prop_soil_texture
      sand_sl(:,:)     = max_prop_soil_texture
    ELSEWHERE (sand_sl(:,:) < min_prop_soil_texture)
      hlp_sand_sl(:,:) = min_prop_soil_texture - sand_sl(:,:)
      ! proportion of hlp_sand_sl(:,:) to be drawn from clay and silt
      hlp_clay_sl(:,:) = clay_sl(:,:) + (1.0_wp - (clay_sl(:,:) + silt_sl(:,:))) / 2.0_wp
      hlp_silt_sl(:,:) = silt_sl(:,:) + (1.0_wp - (clay_sl(:,:) + silt_sl(:,:))) / 2.0_wp
      ! where new clay proportion in soil would be less than "min_prop_soil_texture + 0.01" use only from silt
      WHERE     ( (clay_sl(:,:) - hlp_sand_sl(:,:) * hlp_clay_sl(:,:)) < (min_prop_soil_texture + 0.01) )
        silt_sl(:,:) = silt_sl(:,:) - hlp_sand_sl(:,:)
      ! where new silt proportion in soil would be less than "min_prop_soil_texture + 0.01" use only from clay
      ELSEWHERE ( (silt_sl(:,:) - hlp_sand_sl(:,:) * hlp_silt_sl(:,:)) < (min_prop_soil_texture + 0.01) )
        clay_sl(:,:) = clay_sl(:,:) - hlp_sand_sl(:,:)
      ELSEWHERE
        clay_sl(:,:) = clay_sl(:,:) - hlp_sand_sl(:,:) * hlp_clay_sl(:,:)
        silt_sl(:,:) = silt_sl(:,:) - hlp_sand_sl(:,:) * hlp_silt_sl(:,:)
      ENDWHERE
      sand_sl(:,:) = min_prop_soil_texture
    ELSEWHERE (clay_sl(:,:) < min_prop_soil_texture)
      hlp_clay_sl(:,:) = min_prop_soil_texture - clay_sl(:,:)
      ! proportion of hlp_clay_sl(:,:) to be drawn from sand and silt
      hlp_sand_sl(:,:) = sand_sl(:,:) + (1.0_wp - (sand_sl(:,:) + silt_sl(:,:))) / 2.0_wp
      hlp_silt_sl(:,:) = silt_sl(:,:) + (1.0_wp - (sand_sl(:,:) + silt_sl(:,:))) / 2.0_wp
      ! where new sand proportion in soil would be less than "min_prop_soil_texture + 0.01" use only from silt
      WHERE     ( (sand_sl(:,:) - hlp_clay_sl(:,:) * hlp_sand_sl(:,:)) < (min_prop_soil_texture + 0.01) )
        silt_sl(:,:) = silt_sl(:,:) - hlp_clay_sl(:,:)
      ! where new silt proportion in soil would be less than "min_prop_soil_texture + 0.01" use only from sand
      ELSEWHERE ( (silt_sl(:,:) - hlp_clay_sl(:,:) * hlp_silt_sl(:,:)) < (min_prop_soil_texture + 0.01) )
        sand_sl(:,:) = sand_sl(:,:) - hlp_clay_sl(:,:)
      ELSEWHERE
        sand_sl(:,:) = sand_sl(:,:) - hlp_clay_sl(:,:) * hlp_sand_sl(:,:)
        silt_sl(:,:) = silt_sl(:,:) - hlp_clay_sl(:,:) * hlp_silt_sl(:,:)
      ENDWHERE
      clay_sl(:,:) = min_prop_soil_texture      
    ELSEWHERE (silt_sl(:,:) < min_prop_soil_texture)
      hlp_silt_sl(:,:) = min_prop_soil_texture - silt_sl(:,:)
      ! proportion of hlp_silt_sl(:,:) to be drawn from sand and clay
      hlp_sand_sl(:,:) = sand_sl(:,:) + (1.0_wp - (sand_sl(:,:) + clay_sl(:,:))) / 2.0_wp
      hlp_clay_sl(:,:) = clay_sl(:,:) + (1.0_wp - (sand_sl(:,:) + clay_sl(:,:))) / 2.0_wp
      ! where new sand proportion in soil would be less than "min_prop_soil_texture + 0.01" use only from clay
      WHERE     ( (sand_sl(:,:) - hlp_silt_sl(:,:) * hlp_sand_sl(:,:)) < (min_prop_soil_texture + 0.01) )
        clay_sl(:,:) = clay_sl(:,:) - hlp_silt_sl(:,:)
      ! where new clay proportion in soil would be less than "min_prop_soil_texture + 0.01" use only from sand
      ELSEWHERE ( (clay_sl(:,:) - hlp_silt_sl(:,:) * hlp_clay_sl(:,:)) < (min_prop_soil_texture + 0.01) )
        sand_sl(:,:) = sand_sl(:,:) - hlp_silt_sl(:,:)
      ELSEWHERE
        sand_sl(:,:) = sand_sl(:,:) - hlp_silt_sl(:,:) * hlp_sand_sl(:,:)
        clay_sl(:,:) = clay_sl(:,:) - hlp_silt_sl(:,:) * hlp_clay_sl(:,:)
      ENDWHERE
      silt_sl(:,:) = min_prop_soil_texture
    ENDWHERE
    
    !> 1.1 ensure all adds up to exactly 1.0_wp
    !! 
    clay_sl(:,:) = 1.0_wp
    clay_sl(:,:) = clay_sl(:,:) - (sand_sl(:,:) + silt_sl(:,:))

    !> 1.2 output any changes in soil texture to stdout or log
    !! 
    IF (ABS(SUM(sand_sl(:,:) - sand_original_sl(:,:))) > eps8) THEN
      CALL message(TRIM(routine), ' Modified the soil texture sand input for saxton calculation.')
      print*,TRIM(routine), " Modified the soil texture sand input for saxton calculation: ", sand_sl(:,:) - sand_original_sl(:,:)
    ENDIF
    IF (ABS(SUM(clay_sl(:,:) - clay_original_sl(:,:))) > eps8) THEN
      CALL message(TRIM(routine), ' Modified the soil texture clay input for saxton calculation.')
      print*,TRIM(routine), " Modified the soil texture clay input for saxton calculation: ", clay_sl(:,:) - clay_original_sl(:,:)
    ENDIF
    IF (ABS(SUM(silt_sl(:,:) - silt_original_sl(:,:))) > eps8) THEN
      CALL message(TRIM(routine), ' Modified the soil texture silt input for saxton calculation.')
      print*,TRIM(routine), " Modified the soil texture silt input for saxton calculation: ", silt_sl(:,:) - silt_original_sl(:,:)
    ENDIF

  END SUBROUTINE calc_soil_texture_saxton_compatible


#endif
END MODULE mo_sse_init
