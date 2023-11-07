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
!! @brief Initialization of the the spq memory.
!!
!! @author
!!  Jan Engel, MPI-BGC Jena
!!
!! @par Revision History
!! First version for QUINCY by Jan Engel (20161208) based on the JBACH4 version by Reiner Schnur (2013-04-04)
!! re-write to work onDomain by Jan Engel MPI-BGC (2023-05-30)
!!
MODULE mo_spq_init
#ifndef __NO_JSBACH__

  USE mo_kind,                    ONLY: wp
  USE mo_exception,               ONLY: message_text, finish, message
  USE mo_jsb_control,             ONLY: debug_on
  USE mo_jsb_math_constants,      ONLY: eps8

#ifdef __QUINCY_STANDALONE__
  USE mo_qs_process_init_util,    ONLY: init_spq_soil_properties_spp1685_sites
#endif

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: spq_init

  CHARACTER(len=*), PARAMETER :: modname = 'mo_spq_init'

CONTAINS

  ! ======================================================================================================= !
  !>Intialize SPQ_ process
  !>
  SUBROUTINE spq_init(tile)
    USE mo_jsb_class,              ONLY: Get_model
    USE mo_jsb_tile_class,         ONLY: t_jsb_tile_abstract
    USE mo_jsb_model_class,        ONLY: t_jsb_model
    USE mo_quincy_model_config,    ONLY: QLAND, QPLANT, QSOIL, QCANOPY, Q_TEST_CANOPY, Q_TEST_RADIATION
    USE mo_jsb_lctlib_class,       ONLY: t_lctlib_element
    USE mo_jsb_process_class,      ONLY: SPQ_, VEG_
    USE mo_jsb_grid_class,         ONLY: t_jsb_grid, t_jsb_vgrid
    USE mo_jsb_grid,               ONLY: Get_grid, Get_vgrid
    USE mo_util,                   ONLY: soil_depth_to_layers_2d
    USE mo_jsb_physical_constants, ONLY: Tzero
    USE mo_spq_constants,          ONLY: soil_heat_cap, soil_therm_cond, kdiff_sat_max, &
      &                                  k_pwp_s, k_pwp_c, k_pwp_sc, k_pwp_a, k_pwp_at, k_pwp_bt, &
      &                                  k_fc_s, k_fc_c, k_fc_sc, k_fc_a, k_fc_at, k_fc_bt, k_fc_ct, &
      &                                  k_sat_s, k_sat_c, k_sat_sc, k_sat_a, k_sat_at, k_sat_bt, k_sat_ct, k_sat_dt, &
      &                                  wpot_pwp, wpot_fc, wpot_fc
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Use_config(SPQ_)
    dsl4jsb_Use_config(VEG_)
    dsl4jsb_Use_memory(SPQ_)
    dsl4jsb_Use_memory(VEG_)
    ! ----------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_tile_abstract), INTENT(inout)     :: tile         !< one tile with data structure for one lct
    ! ----------------------------------------------------------------------------------------------------- !
    TYPE(t_jsb_model),      POINTER         :: model              !< the model
    TYPE(t_lctlib_element), POINTER         :: lctlib             !< land-cover-type library - parameter across pft's
    TYPE(t_jsb_grid),       POINTER         :: hgrid              !< Horizontal grid
    TYPE(t_jsb_vgrid),      POINTER         :: vgrid_soil_sb      !< Vertical grid
    TYPE(t_jsb_vgrid),      POINTER         :: vgrid_snow_e       !< Vertical grid snow
    INTEGER                                 :: nsoil_sb           !< number of soil layers as used/defined by the SB_ process
    INTEGER                                 :: nsnow              !< number of snow layers
    REAL(wp), ALLOCATABLE, DIMENSION(:,:)   :: soil_awc
    REAL(wp), ALLOCATABLE, DIMENSION(:,:)   :: hlp2
    REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) :: hlp1
    REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) :: theta_pwp_sl
    REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) :: theta_fc_sl
    REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) :: theta_sat_sl
    REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) :: sand_spq_init_sl
    REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) :: silt_spq_init_sl
    REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) :: clay_spq_init_sl
    CHARACTER(len=3)                        :: site_ID_spp1685
    INTEGER                                 :: isoil, ic, iblk
    INTEGER                                 :: nproma, nblks
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':spq_init'
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Def_config(SPQ_)
    dsl4jsb_Def_config(VEG_)
    dsl4jsb_Def_memory(SPQ_)
    dsl4jsb_Def_memory(VEG_)
    ! ----------------------------------------------------------------------------------------------------- !
    ! VEG_ 3D
    dsl4jsb_Real3D_onDomain    :: root_fraction_sl
    ! SPQ_ 2D
    dsl4jsb_Real2D_onDomain    :: soil_depth
    dsl4jsb_Real2D_onDomain    :: root_depth
    dsl4jsb_Real2D_onDomain    :: t_srf_old 
    dsl4jsb_Real2D_onDomain    :: w_skin
    dsl4jsb_Real2D_onDomain    :: w_soil_root
    dsl4jsb_Real2D_onDomain    :: w_soil_root_pwp
    dsl4jsb_Real2D_onDomain    :: w_soil_root_fc
    dsl4jsb_Real2D_onDomain    :: w_soil_root_theta
    dsl4jsb_Real2D_onDomain    :: w_soil_root_pot
    ! SPQ_ 3D
    dsl4jsb_Real3D_onDomain    :: soil_depth_sl          ! (soil_depth_sl == soil_lay_width_sl)
    dsl4jsb_Real3D_onDomain    :: soil_lay_width_sl
    dsl4jsb_Real3D_onDomain    :: soil_lay_depth_center_sl
    dsl4jsb_Real3D_onDomain    :: soil_lay_depth_ubound_sl
    dsl4jsb_Real3D_onDomain    :: soil_lay_depth_lbound_sl
    dsl4jsb_Real3D_onDomain    :: bulk_dens_sl
    dsl4jsb_Real3D_onDomain    :: heat_capa_sl
    dsl4jsb_Real3D_onDomain    :: therm_cond_sl
    dsl4jsb_Real3D_onDomain    :: w_soil_pwp_sl
    dsl4jsb_Real3D_onDomain    :: w_soil_fc_sl 
    dsl4jsb_Real3D_onDomain    :: w_soil_sat_sl 
    dsl4jsb_Real3D_onDomain    :: w_soil_sl
    dsl4jsb_Real3D_onDomain    :: w_soil_pot_sl
    dsl4jsb_Real3D_onDomain    :: w_ice_sl
    dsl4jsb_Real3D_onDomain    :: saxtonA 
    dsl4jsb_Real3D_onDomain    :: saxtonB 
    dsl4jsb_Real3D_onDomain    :: saxtonC 
    dsl4jsb_Real3D_onDomain    :: kdiff_sat_sl 
    dsl4jsb_Real3D_onDomain    :: sand_sl
    dsl4jsb_Real3D_onDomain    :: silt_sl
    dsl4jsb_Real3D_onDomain    :: clay_sl
    dsl4jsb_Real3D_onDomain    :: volume_min_sl
    dsl4jsb_Real3D_onDomain    :: t_snow_snl
    dsl4jsb_Real3D_onDomain    :: t_soil_sl
    ! ----------------------------------------------------------------------------------------------------- !
    IF (.NOT. tile%Is_process_calculated(SPQ_)) RETURN
    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')
    IF (tile%lcts(1)%lib_id == 0) RETURN                !< run this init only if the present tile is a pft
    ! ----------------------------------------------------------------------------------------------------- !
    model          => Get_model(tile%owner_model_id)
    lctlib         => model%lctlib(tile%lcts(1)%lib_id)
    hgrid          => Get_grid(model%grid_id)
    vgrid_soil_sb  => Get_vgrid('soil_layer_sb')
    vgrid_snow_e   => Get_vgrid('snow_depth_energy')
    nblks          =  hgrid%nblks
    nproma         =  hgrid%nproma
    nsoil_sb       =  vgrid_soil_sb%n_levels 
    nsnow          =  vgrid_snow_e%n_levels
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Get_config(SPQ_)
    dsl4jsb_Get_config(VEG_)
    dsl4jsb_Get_memory(SPQ_)
    dsl4jsb_Get_memory(VEG_)
    ! ----------------------------------------------------------------------------------------------------- !
    ALLOCATE(soil_awc(nproma, nblks))
    ALLOCATE(hlp2(nproma, nblks))
    ALLOCATE(hlp1(nproma, nsoil_sb, nblks))
    ALLOCATE(theta_pwp_sl(nproma, nsoil_sb, nblks))
    ALLOCATE(theta_fc_sl(nproma, nsoil_sb, nblks))
    ALLOCATE(theta_sat_sl(nproma, nsoil_sb, nblks))
    ALLOCATE(sand_spq_init_sl(nproma, nsoil_sb, nblks))
    ALLOCATE(clay_spq_init_sl(nproma, nsoil_sb, nblks))
    ALLOCATE(silt_spq_init_sl(nproma, nsoil_sb, nblks))
    ! ----------------------------------------------------------------------------------------------------- !
    ! VEG_ 3D
    dsl4jsb_Get_var3D_onDomain(VEG_,  root_fraction_sl)
    ! SPQ_ 2D
    dsl4jsb_Get_var2D_onDomain(SPQ_, soil_depth)
    dsl4jsb_Get_var2D_onDomain(SPQ_, root_depth)
    dsl4jsb_Get_var2D_onDomain(SPQ_, t_srf_old)
    dsl4jsb_Get_var2D_onDomain(SPQ_, w_skin)
    dsl4jsb_Get_var2D_onDomain(SPQ_, w_soil_root)
    dsl4jsb_Get_var2D_onDomain(SPQ_, w_soil_root_pwp)
    dsl4jsb_Get_var2D_onDomain(SPQ_, w_soil_root_fc)
    dsl4jsb_Get_var2D_onDomain(SPQ_, w_soil_root_theta)
    dsl4jsb_Get_var2D_onDomain(SPQ_, w_soil_root_pot)
    ! SPQ_ 3D
    dsl4jsb_Get_var3D_onDomain(SPQ_, soil_depth_sl)              ! out   (soil_depth_sl == soil_lay_width_sl)
    dsl4jsb_Get_var3D_onDomain(SPQ_, soil_lay_width_sl)          ! out
    dsl4jsb_Get_var3D_onDomain(SPQ_, soil_lay_depth_center_sl)   ! out
    dsl4jsb_Get_var3D_onDomain(SPQ_, soil_lay_depth_ubound_sl)   ! out
    dsl4jsb_Get_var3D_onDomain(SPQ_, soil_lay_depth_lbound_sl)   ! out
    dsl4jsb_Get_var3D_onDomain(SPQ_, bulk_dens_sl)
    dsl4jsb_Get_var3D_onDomain(SPQ_, heat_capa_sl)
    dsl4jsb_Get_var3D_onDomain(SPQ_, therm_cond_sl)
    dsl4jsb_Get_var3D_onDomain(SPQ_, w_soil_pwp_sl)
    dsl4jsb_Get_var3D_onDomain(SPQ_, w_soil_fc_sl)
    dsl4jsb_Get_var3D_onDomain(SPQ_, w_soil_sat_sl)
    dsl4jsb_Get_var3D_onDomain(SPQ_, w_soil_sl)
    dsl4jsb_Get_var3D_onDomain(SPQ_, w_soil_pot_sl)
    dsl4jsb_Get_var3D_onDomain(SPQ_, w_ice_sl)
    dsl4jsb_Get_var3D_onDomain(SPQ_, saxtonA)
    dsl4jsb_Get_var3D_onDomain(SPQ_, saxtonB)
    dsl4jsb_Get_var3D_onDomain(SPQ_, saxtonC)
    dsl4jsb_Get_var3D_onDomain(SPQ_, kdiff_sat_sl)
    dsl4jsb_Get_var3D_onDomain(SPQ_, sand_sl)
    dsl4jsb_Get_var3D_onDomain(SPQ_, silt_sl)
    dsl4jsb_Get_var3D_onDomain(SPQ_, clay_sl)
    dsl4jsb_Get_var3D_onDomain(SPQ_, volume_min_sl)
    dsl4jsb_Get_var3D_onDomain(SPQ_, t_snow_snl)
    dsl4jsb_Get_var3D_onDomain(SPQ_, t_soil_sl)
    ! ----------------------------------------------------------------------------------------------------- !

    !>0.9 snow init
    !>
    t_snow_snl(:,:,:) = Tzero 
    
    !>1.0 Total soil depth until bedrock
    !>
    ! QUINCY soil_depth for one site (read from namelist)
    soil_depth(:,:) = dsl4jsb_Config(SPQ_)%soil_depth
    ! for all points of the domain soil_depth(:,:) is at least 0.1
    soil_depth(:,:) = MAX(0.1_wp, MERGE(soil_depth(:,:), 0._wp, soil_depth(:,:) >= 0._wp))

    !>2.0 calc actual soil layer depths based on the fixed layer thicknesses from namelist and the bedrock depth
    !>
    ! NOTE: soil_depth_to_layers_2d() returns a 3D variable
    soil_depth_sl(:,:,:) = soil_depth_to_layers_2d(soil_depth(:,:), &     ! Total soil depth until bedrock (from textures), 2D variable
      &                                            vgrid_soil_sb%dz(:))   ! Soil layer thicknesses from namelist
    ! pass these values to a SPQ_ variable with a correct name
    soil_lay_width_sl(:,:,:) = soil_depth_sl(:,:,:)

    !>  2.1 calc more metrics of the soil layers
    !>
    ! i)  lower & upper boundary
    ! ii) depth at the center of each layer
    DO iblk = 1,nblks
      DO isoil = 1,nsoil_sb
        ! lower & upper boundary
        IF (isoil == 1) THEN
          soil_lay_depth_lbound_sl(:, isoil, iblk) = 0.0_wp
          soil_lay_depth_ubound_sl(:, isoil, iblk) = soil_lay_width_sl(:, isoil, iblk)
        ELSE
          soil_lay_depth_lbound_sl(:, isoil, iblk) = soil_lay_depth_lbound_sl(:, isoil-1, iblk) &
            &                                        + soil_lay_width_sl(:, isoil-1, iblk)
          soil_lay_depth_ubound_sl(:, isoil, iblk) = soil_lay_depth_ubound_sl(:, isoil-1, iblk) &
            &                                        + soil_lay_width_sl(:, isoil, iblk)
        ENDIF                                
        ! soil-layer center
        soil_lay_depth_center_sl(:, isoil, iblk) = (soil_lay_depth_lbound_sl(:, isoil, iblk) &
          &                                        + soil_lay_depth_ubound_sl(:, isoil, iblk)) &
          &                                        * 0.5_wp 
      END DO
    END DO

    !>3.0 check if last soil layer has a thickness larger zero
    !>
    DO iblk = 1,nblks
      DO ic  = 1,nproma
        IF (soil_depth_sl(ic, nsoil_sb, iblk) < eps8) THEN
          WRITE(message_text,*) 'The last soil layer has a thickness < eps8. Increase soil_depth from ', soil_depth, &
            &                   ' to at least ', SUM(soil_lay_width_sl(ic, 1:nsoil_sb-1, iblk)), &
            &                   ' + eps4, or more, and re-start the simulation.'
          CALL finish(TRIM(routine), message_text)
        END IF  
      END DO
    END DO
    
    !>4.0 init soil-property variables
    !>
    
    !>  4.1 init with one site-specific value for all soil layers
    !>
    bulk_dens_sl(:,:,:)   = dsl4jsb_Config(SPQ_)%bulk_density 
    heat_capa_sl(:,:,:)   = soil_heat_cap
    therm_cond_sl(:,:,:)  = soil_therm_cond
    sand_sl(:,:,:)        = dsl4jsb_Config(SPQ_)%soil_sand
    silt_sl(:,:,:)        = dsl4jsb_Config(SPQ_)%soil_silt
    clay_sl(:,:,:)        = dsl4jsb_Config(SPQ_)%soil_clay
    volume_min_sl(:,:,:)  = 1._wp
    t_soil_sl(:,:,:)      = 273.15_wp
    t_srf_old(:,:)        = 273.15_wp

#ifdef __QUINCY_STANDALONE__
    !>    4.1.1 soil-layer specific init for sites of the "SPP 1685" site-set
    !>
    IF (model%config%flag_spp1685) THEN
      DO iblk = 1,nblks
        CALL init_spq_soil_properties_spp1685_sites( &
          & nproma, &                           ! in
          & nsoil_sb, &
          & model%config%lon, &                 ! in ('site_longitude')
          & bulk_dens_sl(:, :, iblk), &         ! out
          & sand_sl(:, :, iblk), &
          & silt_sl(:, :, iblk), &
          & clay_sl(:, :, iblk) )               ! out
      END DO
    END IF
#endif

    !>  4.2 ensure a certain minimum proportion for sand, clay and silt
    !>    needed for the Saxton equations
    !>
    ! this is to ensure applying the empirical functions by Saxton within the range of the texture they were designed for \n
    ! avoiding ranges where calculation would return negative theta_pwp_sl or relatively large kdiff_sat_sl values
    ! 
    ! Saxton equations with too small soil-texture fractions may lead to unrealistic percolation rates and soil water potentials
    !
    ! for more details see ticket #560
    !  
    ! using the modified *_spq_init_sl soil texture for "saxton functions here in spq_init" only
    sand_spq_init_sl(:,:,:) = sand_sl(:,:,:)
    clay_spq_init_sl(:,:,:) = clay_sl(:,:,:)
    silt_spq_init_sl(:,:,:) = silt_sl(:,:,:)
    DO iblk = 1,nblks
      CALL calc_soil_texture_saxton_compatible( &
        & nproma, &                         ! in
        & nsoil_sb, &                       ! in
        & sand_spq_init_sl(:, :, iblk), &   ! inout
        & silt_spq_init_sl(:, :, iblk), &   ! inout
        & clay_spq_init_sl(:, :, iblk) )    ! inout
    END DO
    
    !>5.0 soil hydraulic properties, based on Saxton and Rawls 2006, SSSA
    !>

    !>  5.1 fractional water holding capacity at permanent wilting point
    !>  NOTE that this assumes that the soil hydraulic properties are identical between the layers
    !>
    hlp1(:,:,:)         = k_pwp_s * sand_spq_init_sl(:,:,:) + & 
                          k_pwp_c * clay_spq_init_sl(:,:,:) + &
                          k_pwp_sc * sand_spq_init_sl(:,:,:) * clay_spq_init_sl(:,:,:) + &
                          k_pwp_a
    theta_pwp_sl(:,:,:) = k_pwp_at + (1.0_wp + k_pwp_bt) * hlp1(:,:,:) 

    !>  5.2 fractional water holding capacity at field capacity 
    !>
    hlp1(:,:,:)         = k_fc_s * sand_spq_init_sl(:,:,:) + &
                          k_fc_c * clay_spq_init_sl(:,:,:) + &
                          k_fc_sc * sand_spq_init_sl(:,:,:) * clay_spq_init_sl(:,:,:) + &
                          k_fc_a
    theta_fc_sl(:,:,:)  = k_fc_at + (1.0_wp + k_fc_bt) * hlp1(:,:,:) + k_fc_ct * hlp1(:,:,:) ** 2._wp

    !>  5.3 fractional water holding capacity at saturation 
    !>
    hlp1(:,:,:)         = k_sat_s * sand_spq_init_sl(:,:,:) + & 
                          k_sat_c * clay_spq_init_sl(:,:,:) + &
                          k_sat_sc * sand_spq_init_sl(:,:,:) * clay_spq_init_sl(:,:,:) + &
                          k_sat_a
    theta_sat_sl(:,:,:) = theta_fc_sl(:,:,:) + k_sat_at + (1.0_wp + k_sat_bt) * hlp1(:,:,:) - & 
                          k_sat_ct * sand_spq_init_sl(:,:,:) + k_sat_dt 
    
    !>  5.4 parameters describing the moisture-tension curve                 
    !>  note that for convenience, B has a negative sign and A is in MPa 
    saxtonB(:,:,:)      = - (LOG(wpot_pwp) - LOG(wpot_fc)) / (LOG(theta_fc_sl(:,:,:)) - LOG(theta_pwp_sl(:,:,:)))
    saxtonA(:,:,:)      = - (EXP( LOG(wpot_fc) - saxtonB(:,:,:) * LOG(theta_fc_sl(:,:,:)))) / 1000.0_wp

    !>  5.5 parameters describing the moisture-diffusivity curve
    !>
    kdiff_sat_sl(:,:,:) = kdiff_sat_max * (theta_sat_sl(:,:,:) - theta_fc_sl(:,:,:)) ** (3.0_wp + 1._wp / saxtonB(:,:,:))
    saxtonC(:,:,:)      = 3.0_wp - 2._wp * saxtonB(:,:,:)

    !>6.0 calc depth of roots
    !>  given by prescribed AWC (soil water in the rooting zone [m])
    !>
    !>  the init of 'root_fraction_sl(:,:,:)' (in sec 6.3) depends on this init of the root_depth(:)
    !> 
    soil_awc(:,:)         = 0._wp
    w_soil_root_pwp(:,:)  = 0._wp
    w_soil_root_fc(:,:)   = 0._wp
    root_depth(:,:)       = 0._wp
    DO iblk = 1,nblks
      DO isoil = 1,nsoil_sb
        hlp1(:, isoil, iblk) = soil_lay_width_sl(:, isoil, iblk) * (theta_fc_sl(:, isoil, iblk) - theta_pwp_sl(:, isoil, iblk))
        soil_awc(:, iblk) = soil_awc(:, iblk) + hlp1(:, isoil, iblk)
        WHERE (soil_awc(:, iblk) < (dsl4jsb_Config(SPQ_)%soil_awc_prescribe / 1000._wp)) 
          root_depth(:, iblk)      = root_depth(:, iblk) + soil_lay_width_sl(:, isoil, iblk)
          w_soil_root_pwp(:, iblk) = w_soil_root_pwp(:, iblk) + theta_pwp_sl(:, isoil, iblk) &
            &                        * soil_lay_width_sl(:, isoil, iblk)
          w_soil_root_fc(:, iblk)  = w_soil_root_fc(:, iblk) + theta_fc_sl(:, isoil, iblk) &
            &                        * soil_lay_width_sl(:, isoil, iblk)
        ELSEWHERE
          hlp2(:, iblk)            = MAX(0._wp, dsl4jsb_Config(SPQ_)%soil_awc_prescribe &
            &                        / 1000._wp - soil_awc(:, iblk) + hlp1(:, isoil, iblk))
          root_depth(:, iblk)      = root_depth(:, iblk) + hlp2(:, iblk) &
            &                        / (theta_fc_sl(:, isoil, iblk) - theta_pwp_sl(:, isoil, iblk))
          w_soil_root_pwp(:, iblk) = w_soil_root_pwp(:, iblk) + theta_pwp_sl(:, isoil, iblk) * hlp2(:, iblk) &
            &                        / (theta_fc_sl(:, isoil, iblk) - theta_pwp_sl(:, isoil, iblk))
          w_soil_root_fc(:, iblk)  = w_soil_root_fc(:, iblk) + theta_fc_sl(:, isoil, iblk) * hlp2(:, iblk) &
            &                        / (theta_fc_sl(:, isoil, iblk) - theta_pwp_sl(:, isoil, iblk))
        END WHERE
      END DO
    END DO

    !>  6.1 water in the surface reservoir [m]
    !>
    w_skin(:,:) = 0.0_wp

    !>  6.2 adjust initial water amount to prescribed theta 
    !>
    w_soil_root_theta(:,:)  = dsl4jsb_Config(SPQ_)%soil_theta_prescribe
    w_soil_root(:,:)        = w_soil_root_pwp(:,:) &
      &                     + dsl4jsb_Config(SPQ_)%soil_theta_prescribe * (w_soil_root_fc(:,:) - w_soil_root_pwp(:,:)) 
    w_soil_root_pot(:,:)    = saxtonA(:, 1, :) * (MIN(w_soil_root_fc,w_soil_root(:,:)) / root_depth(:,:)) ** saxtonB(:, 1, :) 

    !>  6.3 root fraction 
    !>    root_fraction_sl(:,:,:) is a variable of veg_mem and depends on the init of root_depth(:,:), see section 6.0
    !>
    ! by default use: fixed roots for QCANOPY, Q_TEST_CANOPY, Q_TEST_RADIATION
    !                 because these models are not growing roots, suffering water limitations with dyn roots
    SELECT CASE(model%config%qmodel_id)
    CASE(QCANOPY, Q_TEST_CANOPY, Q_TEST_RADIATION)
      dsl4jsb_Config(VEG_)%flag_dynamic_roots = .FALSE.
      CALL message(TRIM(routine), ' set "flag_dynamic_roots = .FALSE." (using the canopy, test_canopy or test_radiation model)')
      print*,      TRIM(routine), ' set "flag_dynamic_roots = .FALSE."', &
        &                         ' (which is the default for the canopy, test_canopy or test_radiation models)'
    END SELECT

    ! init root fraction
    IF (dsl4jsb_Config(VEG_)%flag_dynamic_roots) THEN
      !>    6.3.1 dynamic root fractions across soil layers (see mo_veg_growth)
      !>
      root_fraction_sl(:, 1, :)          = 1.0_wp
      root_fraction_sl(:, 2:nsoil_sb, :) = 0.0_wp
    ELSE
      !>    6.3.2 fixed root fractions across soil layers (no changes after init)
      !>
      ! root fractions assuming an exponential distribution of roots \n
      !   Cr_l = Cr_0 * exp(-lctlib%k_root_dist*depth) \n
      ! integral to bottom of root zone: 1/lctlib%k_root_dist * (1 - exp(-lctlib%k_root_dist*depth) \n
      !   thus for each layer
      !
      ! NOTE: upper bound of a layer is the larger value compared to lower bound !
      ! 
      ! layer 1 (vgrid_soil_sb)
      root_fraction_sl(:, 1, :) = (1._wp - EXP(-lctlib%k_root_dist * soil_lay_depth_ubound_sl(:, 1, :))) &
        &                         / (1._wp - EXP(-lctlib%k_root_dist * root_depth(:,:)))
      ! layers 2 to x (vgrid_soil_sb)
      DO iblk = 1,nblks
        DO isoil = 2,nsoil_sb
          WHERE(root_depth(:, iblk) > soil_lay_depth_ubound_sl(:, isoil-1, iblk))
            WHERE(root_depth(:, iblk) > soil_lay_depth_ubound_sl(:, isoil, iblk))
              root_fraction_sl(:, isoil, iblk) = (EXP(-lctlib%k_root_dist * soil_lay_depth_ubound_sl(:, isoil-1, iblk)) &
                &                                - EXP(-lctlib%k_root_dist * soil_lay_depth_ubound_sl(:, isoil, iblk))) &
                &                                / (1._wp - EXP(-lctlib%k_root_dist * root_depth(:, iblk)))
            ELSEWHERE
              root_fraction_sl(:, isoil, iblk) = (EXP(-lctlib%k_root_dist * soil_lay_depth_ubound_sl(:, isoil-1, iblk)) &
                &                                - EXP(-lctlib%k_root_dist * root_depth(:, iblk))) &
                &                                / (1._wp - EXP(-lctlib%k_root_dist * root_depth(:, iblk)))
            ENDWHERE
          ELSEWHERE
            root_fraction_sl(:, isoil, iblk) = 0.0_wp 
          ENDWHERE
        END DO
      END DO
    END IF  ! flag_dynamic_roots

    !>7.0 soil water and ice pwp, fc and sat per layer depth
    !>
    w_soil_pwp_sl(:,:,:)  = 0.0_wp
    w_soil_fc_sl(:,:,:)   = 0.0_wp
    w_soil_sat_sl(:,:,:)  = 0.0_wp
    DO iblk = 1,nblks
      DO isoil = 1,nsoil_sb
        WHERE(soil_depth_sl(:, isoil, iblk) > 0.0_wp) 
          w_soil_pwp_sl(:, isoil, iblk) = theta_pwp_sl(:, isoil, iblk) * soil_depth_sl(:, isoil, iblk) 
          w_soil_fc_sl(:, isoil, iblk)  = theta_fc_sl(:, isoil, iblk) * soil_depth_sl(:, isoil, iblk) 
          w_soil_sat_sl(:, isoil, iblk) = theta_sat_sl(:, isoil, iblk) * soil_depth_sl(:, isoil, iblk) 
        ENDWHERE
      END DO
    END DO
    
    !>  7.1 adjust initial water amount to prescribed theta
    !>  avoid water stress (i.e., dry soil) of plants at the first timestep
    !>
    w_soil_sl(:,:,:) = w_soil_pwp_sl(:,:,:) &
      &                + (dsl4jsb_Config(SPQ_)%soil_theta_prescribe * (w_soil_fc_sl(:,:,:) - w_soil_pwp_sl(:,:,:)))

    !>  7.2 soil water potential per layer
    !>
    WHERE(soil_depth_sl(:,:,:) > eps8)
      hlp1(:,:,:)          = MIN(w_soil_fc_sl(:,:,:), w_soil_sl(:,:,:)) / soil_depth_sl(:,:,:)
      w_soil_pot_sl(:,:,:) = saxtonA(:,:,:) * hlp1(:,:,:) ** saxtonB(:,:,:)
    ELSEWHERE
      w_soil_pot_sl(:,:,:) = 0.0_wp 
    ENDWHERE

    ! Deallocate local allocatable variables
    DEALLOCATE(soil_awc)
    DEALLOCATE(hlp2)
    DEALLOCATE(hlp1)
    DEALLOCATE(theta_pwp_sl)
    DEALLOCATE(theta_fc_sl)
    DEALLOCATE(theta_sat_sl)
    DEALLOCATE(sand_spq_init_sl)
    DEALLOCATE(clay_spq_init_sl)
    DEALLOCATE(silt_spq_init_sl)

  END SUBROUTINE spq_init

  ! ======================================================================================================= !
  !>ensure a certain minimum proportion for sand, clay and silt
  !> 
  !> Saxton equations with too small soil-texture fractions may lead to unrealistic percolation rates and soil water potentials
  !>
  !> using these modified soil texture for "saxton functions in spq_init" only
  !>
  SUBROUTINE calc_soil_texture_saxton_compatible( &
    & nproma, &
    & nsoil_sb, &
    & sand_sl, &
    & silt_sl, &
    & clay_sl )

    INTEGER,                                INTENT(in)    :: nproma        !< dimensions
    INTEGER,                                INTENT(in)    :: nsoil_sb      !< number of soil layers
    REAL(wp), DIMENSION(nproma, nsoil_sb),  INTENT(inout) :: sand_sl       !< sand proportion per layer
    REAL(wp), DIMENSION(nproma, nsoil_sb),  INTENT(inout) :: silt_sl       !< silt proportion per layer
    REAL(wp), DIMENSION(nproma, nsoil_sb),  INTENT(inout) :: clay_sl       !< clay proportion per layer
    ! ----------------------------------------------------------------------------------------------------- !
    REAL(wp), DIMENSION(nproma, nsoil_sb) :: hlp_sand_sl, &
                                             hlp_silt_sl, &
                                             hlp_clay_sl
    REAL(wp), DIMENSION(nproma, nsoil_sb) :: sand_original_sl, &
                                             silt_original_sl, &
                                             clay_original_sl
    REAL(wp)                              :: min_prop_soil_texture, &
                                             max_prop_soil_texture
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_soil_texture_saxton_compatible'
    ! ----------------------------------------------------------------------------------------------------- !

    !>0.9 init local variables
    !>
    min_prop_soil_texture = 0.2_wp
    max_prop_soil_texture = 1.0_wp - (2.0_wp * min_prop_soil_texture)
    hlp_sand_sl(:,:)      = 0.0_wp
    hlp_clay_sl(:,:)      = 0.0_wp
    hlp_silt_sl(:,:)      = 0.0_wp
    sand_original_sl(:,:) = sand_sl(:,:)
    clay_original_sl(:,:) = clay_sl(:,:)
    silt_original_sl(:,:) = silt_sl(:,:)

    !>1.0 ensure a minimum proportion 'min_prop_soil_texture' for sand, clay and silt
    !>  this is to ensure applying the empirical functions by Saxton within the range of the texture they were designed for
    !>
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
    
    !>  1.1 ensure all adds up to exactly 1.0_wp
    !>
    clay_sl(:,:) = 1.0_wp
    clay_sl(:,:) = clay_sl(:,:) - (sand_sl(:,:) + silt_sl(:,:))

    !>  1.2 output any changes in soil texture to stdout or log
    !>
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
END MODULE mo_spq_init
