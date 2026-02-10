!> Init routines specific to QUINCY standalone
!>
!> QUINCY standalone
!>
!> ---------------------------------------
!> Copyright (C) 2013-2026, MPI-M, MPI-BGC
!>
!> Contact: https://doi.org/10.17871/quincy-model-2019
!> Authors: AUTHORS.md
!> See LICENSES/ for license information
!> SPDX-License-Identifier: BSD-3-Clause
!> ---------------------------------------
!>
!>#### Provides routines related to init of processes specific to QUINCY standalone (site-level)
!>
MODULE mo_qs_process_init_util

  USE mo_kind,                    ONLY: wp
  USE mo_exception,               ONLY: finish
  USE mo_jsb_math_constants,      ONLY: eps8

  IMPLICIT NONE
  
  PRIVATE
  PUBLIC :: init_spq_soil_properties_spp1685_sites
  PUBLIC :: init_sb_soil_properties_spp1685_sites, init_sb_soil_p_pools_spp1685_sites

  CHARACTER(len=*), PARAMETER, PRIVATE :: modname = 'mo_qs_process_init_util'

CONTAINS

  ! ======================================================================================================= !
  !>init soil properties for sites of the 'SPP 1685' site-set
  !>
  !>  soil-layer specific texture (sand silt clay) init for each particular site
  !>
  SUBROUTINE init_spq_soil_properties_spp1685_sites( &
    & nproma, &
    & nsoil_sb, &
    & site_longitude, &
    & bulk_dens_sl, &
    & sand_sl, &
    & silt_sl, &
    & clay_sl )
  
    INTEGER,                                INTENT(in)  :: nproma          !< dimensions
    INTEGER,                                INTENT(in)  :: nsoil_sb        !< number of soil layers
    REAL(wp),                               INTENT(in)  :: site_longitude  !< lon coordinate to identify site of SPP site-set
    REAL(wp), DIMENSION(nproma, nsoil_sb),  INTENT(out) :: bulk_dens_sl    !< bulk density ...
    REAL(wp), DIMENSION(nproma, nsoil_sb),  INTENT(out) :: sand_sl         !< sand proportion per layer
    REAL(wp), DIMENSION(nproma, nsoil_sb),  INTENT(out) :: silt_sl         !< silt proportion per layer
    REAL(wp), DIMENSION(nproma, nsoil_sb),  INTENT(out) :: clay_sl         !< clay proportion per layer
    ! ----------------------------------------------------------------------------------------------------- !
    INTEGER                     :: ic
    CHARACTER(len=3)            :: site_ID_spp1685
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':init_spq_soil_properties_spp1685_sites'
    ! ----------------------------------------------------------------------------------------------------- !

    !>0.9 init output variables
    !>
    bulk_dens_sl(:,:) = 0.0_wp
    sand_sl(:,:)      = 0.0_wp
    silt_sl(:,:)      = 0.0_wp
    clay_sl(:,:)      = 0.0_wp

    !>1.0 check number of soil layers
    !>
    IF (nsoil_sb /= 15) THEN
      CALL finish(TRIM(routine),'nsoil_sb is not 15 | Simulations with the SPP_1685 site-set work only with 15 soil layers.')
    ENDIF
  
    !>2.0 get site ID based on longitude
    !>
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
      ! error
      site_ID_spp1685 = "ERR"
    ENDIF

    !>3.0 set soil properties per soil layer - soecific to the SPP site
    !>
    ! the below texture data is based on: Lang, F. et al. (2017 ) "Soil phosphorus supply controls P nutrition strategies 
    !   of beech forest ecosystems in Central Europe". Biogeochemistry, doi:10.1007/s10533-017-0375-0. Table S2 in the
    !   Supplementary material 2.
    !   These texture data were extrapolated: We use the measured soil texture of the topmost 1m soil for all the five sites 
    !   (presented in this table S2). We extrapolated them to the whole soil depth (9.5m) following the pattern of 
    !   texture change against depth increase observed in the topmost 1m soil.
    SELECT CASE(site_ID_spp1685)
    CASE ("BBR")
      DO ic = 1,nproma
        bulk_dens_sl(ic,:)   = (/1000._wp,1080._wp,1200._wp,1320._wp,1450._wp, &
                                1520._wp,1600._wp,1650._wp,1700._wp,1750._wp, &
                                1800._wp, 1850._wp, 1900._wp, 1950._wp, 2000._wp/)
        sand_sl(ic,:)   = (/0.07_wp,0.09_wp,0.1_wp,0.1228_wp,0.1785_wp, &
                           .24_wp,0.2726_wp,0.3519_wp,0.4_wp,0.45_wp, &
                           .5_wp,0.6_wp,0.7_wp,0.8_wp,0.8_wp/)
        silt_sl(ic,:)   = (/0.56_wp,0.54_wp,0.55_wp,0.5386_wp,0.5184_wp, &
                           .5_wp,0.5042_wp,0.4868_wp,0.4_wp,0.37_wp, &
                           .35_wp,0.3_wp,0.22_wp,0.15_wp,0.15_wp/)
      ENDDO
    CASE ("LUE")
      DO ic = 1,nproma
        bulk_dens_sl(ic,:)   = (/1570._wp,1600._wp,1610._wp,1630._wp,1660._wp, &
                                1800._wp,2000._wp,2000._wp,2200._wp,2200._wp, &
                                2200._wp, 2200._wp, 2200._wp, 2200._wp, 2200._wp/)
        sand_sl(ic,:)   = (/0.75_wp,0.76_wp,0.81_wp,0.7543_wp,0.7448_wp, &
                           0.713_wp,0.7626_wp,0.8185_wp,0.84_wp,0.85_wp, &
                           0.85_wp,0.88_wp,0.88_wp,0.9_wp,0.92_wp/)
        silt_sl(ic,:)   = (/0.19_wp,0.19_wp,0.14_wp,0.1864_wp,0.1926_wp, &
                           0.2085_wp,0.1753_wp,0.1177_wp,0.1_wp,0.1_wp, &
                           0.1_wp,0.08_wp,0.08_wp,0.06_wp,0.05_wp/)
      ENDDO    
    CASE ("VES")
      DO ic = 1,nproma
        bulk_dens_sl(ic,:)   = (/1000._wp,1000._wp,1000._wp,1000._wp,1050._wp, &
                                1050._wp,1050._wp,1200._wp,1250._wp,1300._wp, &
                                1350._wp, 1400._wp, 1450._wp, 1500._wp, 1600._wp/)
        sand_sl(ic,:)   = (/0.30_wp,0.3245_wp,0.3746_wp,0.4057_wp,0.4494_wp, &
                           0.4923_wp,0.5627_wp,0.6669_wp,0.71_wp,0.75_wp, &
                           0.77_wp,0.79_wp,0.83_wp,0.86_wp,0.88_wp/)
        silt_sl(ic,:)   = (/0.4638_wp,0.47_wp,0.4562_wp,0.4579_wp,0.4386_wp, &
                           0.4012_wp,0.3578_wp,0.2864_wp,0.25_wp,0.22_wp, &
                           0.2_wp,0.18_wp,0.15_wp,0.12_wp,0.1_wp/)
      ENDDO
    CASE ("MIT")
      DO ic = 1,nproma
        bulk_dens_sl(ic,:)   = (/1000._wp,1000._wp,1000._wp,1000._wp,1050._wp, &
                                1050._wp,1100._wp,1100._wp,1400._wp,1450._wp, &
                                1450._wp, 1500._wp, 1550._wp, 1600._wp, 1650._wp/)
        sand_sl(ic,:)   = (/0.447_wp,0.5115_wp,0.5102_wp,0.52_wp,0.5642_wp, &
                           0.57_wp,0.56_wp,0.5681_wp,0.60_wp,0.62_wp, &
                           0.64_wp,0.65_wp,0.70_wp,0.70_wp,0.70_wp/)
        silt_sl(ic,:)   = (/0.3223_wp,0.3115_wp,0.3148_wp,0.32_wp,0.3195_wp, &
                           0.34_wp,0.3526_wp,0.3601_wp,0.33_wp,0.32_wp, &
                           0.31_wp,0.30_wp,0.26_wp,0.26_wp,0.27_wp/)
      ENDDO
    CASE ("CON")
      DO ic = 1,nproma
        bulk_dens_sl(ic,:)   = (/1000._wp,1000._wp,1000._wp,1000._wp,1560._wp, &
                                1760._wp,1760._wp,1760._wp,1800._wp,2000._wp, &
                                2200._wp, 2200._wp, 2450._wp, 2450._wp, 2450._wp/)
        sand_sl(ic,:)   = (/0.40_wp,0.3969_wp,0.39_wp,0.40_wp,0.41_wp, &
                           0.4447_wp,0.4872_wp,0.5159_wp,0.55_wp,0.58_wp, &
                           0.60_wp,0.62_wp,0.65_wp,0.70_wp,0.70_wp/)
        silt_sl(ic,:)   = (/0.33_wp,0.3454_wp,0.38_wp,0.40_wp,0.42_wp, &
                           0.40_wp,0.3992_wp,0.3961_wp,0.37_wp,0.35_wp, &
                           0.33_wp,0.31_wp,0.29_wp,0.25_wp,0.25_wp/)
      ENDDO
    CASE DEFAULT
      CALL finish(TRIM(routine),'This site seems to be not a part of the SPP_1685 site-set although flag_spp1685 is set T.')
    END SELECT
  
    !>4.0 calculate clay content from sand and silt
    !>
    clay_sl(:,:) = 1.0_wp
    clay_sl(:,:) = clay_sl(:,:) - (sand_sl(:,:) + silt_sl(:,:))
  END SUBROUTINE init_spq_soil_properties_spp1685_sites

  ! ======================================================================================================= !
  !>init soil properties for SPP sites only
  !> 
  !>  Output: ph and Qmax_AlFe_cor prescribed for spp sites
  !>
  !>  Qmax_AlFe_cor is further modified with calc_Psorption_parameter()
  !>
  SUBROUTINE init_sb_soil_properties_spp1685_sites( &
    & nproma, &
    & nsoil_sb, &
    & site_ID_spp1685, &
    & ph_sl, &
    & Qmax_AlFe_cor)

    INTEGER,                                INTENT(in)    :: nproma           !< dimensions
    INTEGER,                                INTENT(in)    :: nsoil_sb         !< dimensions
    CHARACTER(len=*), OPTIONAL,             INTENT(in)    :: site_ID_spp1685  !< site ID for spp1685  
    REAL(wp), DIMENSION(nproma, nsoil_sb),  INTENT(out)   :: ph_sl            !< ph of soil layer [-]
    REAL(wp), DIMENSION(nproma, nsoil_sb),  INTENT(out)   :: Qmax_AlFe_cor    !< qmax correlation factor of Al/Fe content [-]  
    ! ----------------------------------------------------------------------------------------------------- !
    INTEGER                     :: ic    !< loop over nproma
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':init_sb_soil_properties_spp1685_sites'
    ! ----------------------------------------------------------------------------------------------------- !
  
    !>1.0 prescribe ph and Qmax_AlFe_cor from measurements
    !>
    SELECT CASE(site_ID_spp1685)
    CASE ("BBR")
      DO ic = 1,nproma
        ph_sl(ic, 1:nsoil_sb)         = (/3.84_wp, 4.14_wp, 4.45_wp, 4.66_wp, 4.85_wp, &
                                        4.92_wp,5.19_wp,5.64_wp,5.8_wp,6._wp, &
                                        6.2_wp, 6.5_wp, 6.5_wp, 6.5_wp, 6.5_wp/)
        Qmax_AlFe_cor(ic, 1:nsoil_sb) = (/3.04_wp, 3.60_wp, 3.60_wp, 3.295_wp, 2.92_wp, &
                                        2.75_wp,2.40_wp,2.04_wp,2._wp,2._wp, &
                                        2._wp, 2._wp, 2._wp, 2._wp, 2._wp/)
      ENDDO
    CASE ("LUE")
      DO ic = 1,nproma
        ph_sl(ic, 1:nsoil_sb)         = (/3.51_wp,3.78_wp,3.84_wp,4.14_wp,4.20_wp, &
                                        4.46_wp,4.45_wp,4.62_wp,4.63_wp,4.8_wp, &
                                        5._wp, 5.2_wp, 5.5_wp, 6._wp, 6.5_wp/)
        Qmax_AlFe_cor(ic, 1:nsoil_sb) = (/0.17_wp, 0.29_wp, 0.41_wp, 0.41_wp, 0.396_wp, &
                                        0.3_wp,0.289_wp,0.236_wp,0.2_wp,0.2_wp, &
                                        0.2_wp, 0.2_wp, 0.2_wp, 0.2_wp, 0.2_wp/)
      ENDDO
    CASE ("VES")
      DO ic = 1,nproma
        ph_sl(ic, 1:nsoil_sb)         = (/3.47_wp,3.91_wp,4.25_wp,4.46_wp,4.52_wp, &
                                        4.66_wp,4.78_wp,4.77_wp,5._wp,5.2_wp, &
                                        5.5_wp, 5.8_wp, 6.0_wp, 6.5_wp, 6.5_wp/)
        Qmax_AlFe_cor(ic, 1:nsoil_sb) = (/1.275_wp, 1.515_wp, 1.378_wp, 1.459_wp, 1.90_wp, &
                                        2.469_wp,3.428_wp,5.33_wp,6._wp,6._wp, &
                                        6._wp, 6._wp, 6._wp, 6._wp, 6._wp/)
      ENDDO
    CASE ("MIT")
      DO ic = 1,nproma
        ph_sl(ic, 1:nsoil_sb)         = (/3.64_wp,4.02_wp,4.53_wp,4.46_wp,4.46_wp, &
                                        4.58_wp,4.56_wp,4.64_wp,4.7_wp,4.8_wp, &
                                        5._wp, 5.2_wp, 5.5_wp, 5.8_wp, 6._wp/)
        Qmax_AlFe_cor(ic, 1:nsoil_sb) = (/1.326_wp, 1.472_wp, 1.573_wp, 1.5_wp, 1.426_wp, &
                                        1.3_wp,1.258_wp,1.2_wp,1.2_wp,1.2_wp, &
                                        1.1_wp, 1.1_wp, 1._wp, 1._wp, 1._wp/)
      ENDDO
    CASE ("CON")
      DO ic = 1,nproma
        ph_sl(ic, 1:nsoil_sb)         = (/4.03_wp,4.12_wp,4.32_wp,4.37_wp,4.42_wp, &
                                        4.54_wp,4.60_wp,4.61_wp,4.7_wp,4.8_wp, &
                                        5._wp, 5.2_wp, 5.5_wp, 5.8_wp, 6._wp/)
        Qmax_AlFe_cor(ic, 1:nsoil_sb) = (/1.9_wp, 1.891_wp, 1.87_wp, 1.8_wp, 1.73_wp, &
                                        1.588_wp,1.436_wp,1.24_wp,1.2_wp,1.2_wp, &
                                        1.1_wp, 1.1_wp, 1._wp, 1._wp, 1._wp/)
      ENDDO
    END SELECT
  END SUBROUTINE init_sb_soil_properties_spp1685_sites
  
  ! ======================================================================================================= !
  !>SPP site-set - initialise mineral soil P pools
  !> 
  !>
  SUBROUTINE init_sb_soil_p_pools_spp1685_sites( &
    & nproma, &
    & nsoil_sb, &
    & soil_lay_depth_center_sl, &
    & soil_p_labile, &
    & soil_p_slow, &
    & soil_p_occluded, &
    & soil_p_primary, &
    & sitename, &
    & p_labile_slow_ratio, &
    & p_labile, &
    & po4_assoc_slow, &
    & po4_occluded, &
    & po4_primary)
  
    USE mo_sb_constants        ! ...
    ! ----------------------------------------------------------------------------------------------------- !
    INTEGER                              , INTENT(in)     :: nproma                    !< dimensions
    INTEGER                              , INTENT(in)     :: nsoil_sb                  !< number of soil layers as used/defined by the SB_ process
    REAL(wp), DIMENSION(nproma, nsoil_sb), INTENT(in)     :: soil_lay_depth_center_sl  !< depth at the center of each soil layer
    REAL(wp), DIMENSION(nproma, nsoil_sb), INTENT(in)     :: p_labile_slow_ratio   !< ratio between labile and labile+slow P pools
    REAL(wp), DIMENSION(nproma, nsoil_sb), INTENT(in)     :: soil_p_labile, &      !< labile P from Yang et al. 2013, g/m2
                                                             soil_p_slow, &        !< secondary P from Yang et al. 2013, g/m2
                                                             soil_p_occluded, &    !< Occluded P from Yang et al. 2013, g/m2
                                                             soil_p_primary        !< Primary from Yang et al. 2013, g/m2
    CHARACTER(len=*)                     , INTENT(in)     :: sitename              !< spp site name: BBR or LUE (MIT, VES, CON)
    REAL(wp), DIMENSION(nproma, nsoil_sb), INTENT(inout)  :: p_labile, &           !< PO4 in solution per soil layer [mol / m3]
                                                             po4_assoc_slow, &     !< slow minerally associated PO4 pool per soil layer [mol/m3]
                                                             po4_occluded, &       !< occluded PO4 pool per soil layer [mol/m3]
                                                             po4_primary           !< primary PO4 pool per soil layer [mol/m3]
    ! ----------------------------------------------------------------------------------------------------- !
    INTEGER                                      :: isoil, ic
    REAL(wp)                                     :: hlp1, hlp2, hlp3, hlp4, hlp5, hlp6
    REAL(wp), DIMENSION(nproma)                  :: arr_hlp2, arr_hlp3, arr_hlp4, soil_p_layer
    REAL(wp), DIMENSION(nsoil_sb)                :: arr_hlp5
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':init_sb_soil_p_pools_spp1685_sites'
    ! ----------------------------------------------------------------------------------------------------- !
    
    !>1.0 set site-level soil properties
    !>
    SELECT CASE (sitename)
      CASE("BBR")
        arr_hlp5(:) = (/9.24_wp,36.92_wp,39.44_wp,40.36_wp,42.75_wp, &  
                       45.72_wp,32.45_wp,95.89_wp,124.14_wp,133.93_wp, &
                       100._wp,126._wp,122.18_wp,145._wp,90._wp/)       ! observed total SIP, mol P m-3
        hlp1            = soil_p_labile(1,1)     ! p_labile from Yang et al. 2013, g/m2
        hlp2            = soil_p_slow(1,1)       ! p_slow from Yang et al. 2013, g/m2
        hlp3            = soil_p_primary(1,1)    ! p_occluded from Yang et al. 2013, g/m2
        hlp4            = soil_p_occluded(1,1)   ! p_primary from Yang et al. 2013, g/m2
        soil_p_layer(:) = 6._wp             ! field measurement depth corresponded model soil layers
        ! alternative: use [mg P/kg soil] measurement to upscale the data from Yang et al. 2013
        !hlp5        = 1750._wp / 1000._wp ! g P / kg soil
        !arr_hlp5(:) = hlp5 * bulk_dens_corr_sl(1,:) / molar_mass_P ! SIP, mol P m-3
      CASE("LUE")
        arr_hlp5(:) = (/0.8867_wp,3.187_wp,2.662_wp,4.887_wp,4.792_wp, &  
                       2.987_wp,3.313_wp,6.0616_wp,7.56_wp,5.787_wp, &
                       5.932_wp,6.845_wp,8.294_wp,8.413_wp,8.532_wp/)   ! observed total SIP, mol P m-3
        hlp1            = soil_p_labile(1,1)     ! p_labile from Yang et al. 2013, g/m2
        hlp2            = soil_p_slow(1,1)       ! p_slow from Yang et al. 2013, g/m2
        hlp3            = soil_p_occluded(1,1)   ! p_occluded from Yang et al. 2013, g/m2
        hlp4            = soil_p_primary (1,1)   ! p_primary from Yang et al. 2013, g/m2
        soil_p_layer(:) = 6._wp             ! field measurement depth corresponded model soil layers
        ! alternative: use [mg P/kg soil] measurement to upscale the data from Yang et al. 2013
        !hlp5        = 120._wp / 1000._wp ! g P / kg soil
        !arr_hlp5(:) = hlp5 * bulk_dens_corr_sl(1,:) / molar_mass_P ! SIP, mol P m-3
      CASE("VES")
        arr_hlp5(:) = (/2.388_wp,8.509_wp,9.033_wp,10.811_wp,15.365_wp, &  
                       20.661_wp,22.721_wp,29.321_wp,30.952_wp,32.821_wp, &
                       40.673_wp,41.880_wp,41.111_wp,48.438_wp,54._wp/)   ! observed total SIP, mol P m-3
        hlp1            = soil_p_labile(1,1)     ! p_labile from Yang et al. 2013, g/m2
        hlp2            = soil_p_slow(1,1)       ! p_slow from Yang et al. 2013, g/m2
        hlp3            = soil_p_primary(1,1)    ! p_occluded from Yang et al. 2013, g/m2
        hlp4            = soil_p_occluded(1,1)   ! p_primary from Yang et al. 2013, g/m2
        soil_p_layer(:) = 6._wp             ! field measurement depth corresponded model soil layers
        ! alternative: use [mg P/kg soil] measurement to upscale the data from Yang et al. 2013
        !hlp5        = 900._wp / 1000._wp ! g P / kg soil
        !arr_hlp5(:) = hlp5 * bulk_dens_corr_sl(1,:) / molar_mass_P ! SIP, mol P m-3
      CASE("MIT")
        arr_hlp5(:) = (/3.723_wp,6.272_wp,7.556_wp,8.292_wp,6.93_wp, &  
                       13.731_wp,13.296_wp,16.511_wp,12.890_wp,15._wp, &
                       16.667_wp,16.875_wp,17.5_wp,15.104_wp,15.625_wp/)  ! observed total SIP, mol P m-3
        hlp1            = soil_p_labile(1,1)     ! p_labile from Yang et al. 2013, g/m2
        hlp2            = soil_p_slow(1,1)       ! p_slow from Yang et al. 2013, g/m2
        hlp3            = soil_p_primary(1,1)    ! p_occluded from Yang et al. 2013, g/m2
        hlp4            = soil_p_occluded(1,1)   ! p_primary from Yang et al. 2013, g/m2
        soil_p_layer(:) = 9._wp             ! field measurement depth corresponded model soil layers
        ! alternative: use [mg P/kg soil] measurement to upscale the data from Yang et al. 2013
        !hlp5        = 900._wp / 1000._wp ! g P / kg soil
        !arr_hlp5(:) = hlp5 * bulk_dens_corr_sl(1,:) / molar_mass_P ! SIP, mol P m-3
      CASE("CON")
        arr_hlp5(:) = (/1.862_wp,2.924_wp,4.027_wp,4.048_wp,4.069_wp, &  
                       7.167_wp,10.855_wp,7.719_wp,18.75_wp,18.33_wp, &
                       23.44_wp,25._wp,21.093_wp,21.875_wp,15.104_wp/)   ! observed total SIP, mol P m-3
        hlp1            = soil_p_labile(1,1)     ! p_labile from Yang et al. 2013, g/m2
        hlp2            = soil_p_slow(1,1)       ! p_slow from Yang et al. 2013, g/m2
        hlp3            = soil_p_primary(1,1)    ! p_occluded from Yang et al. 2013, g/m2
        hlp4            = soil_p_occluded(1,1)   ! p_primary from Yang et al. 2013, g/m2
        soil_p_layer(:) = 9._wp             ! field measurement depth corresponded model soil layers
        ! alternative: use [mg P/kg soil] measurement to upscale the data from Yang et al. 2013
        !hlp5        = 900._wp / 1000._wp ! g P / kg soil
        !arr_hlp5(:) = hlp5 * bulk_dens_corr_sl(1,:) / molar_mass_P ! SIP, mol P m-3
    END SELECT
  
    !>2.0 calculate the initial slow P, occluded P, primary P pool, and p_labile
    !>
    ! @TODO this is inconsistent with section 1.0 of this routine, as section 1.0 does not work for multiple points in a chunk
    ! based on the assumption that soil P pool are linearly corrected with soil weight
    ! the depth correction of initial P pool distributions are fitted against the observed pattern of a young cambisol soil 
    ! with basalt as parent material (BBR, Lang et al. 2017) that, 
    ! the proportion of: labile+sorbed P decreases with increasing depth, and primary P increases with increasing depth (Yu et al., 2020)
    DO isoil = 1, nsoil_sb
      ! primary corrected for depth, [g P] 
      arr_hlp4(:) = MIN(k2_p_pool_depth_corr * soil_lay_depth_center_sl(:,isoil) ** k4_p_pool_depth_corr, &
                    (k1_p_pool_depth_corr - k_p_pool_depth_corr_min_occlud)) / k1_p_pool_depth_corr
      ! labile + slow
      arr_hlp2(:) = k3_p_pool_depth_corr ** (k2_p_pool_depth_corr * soil_lay_depth_center_sl(:,isoil)) / k1_p_pool_depth_corr
      ! occluded corrected for depth
      arr_hlp3(:) = MAX((1.0_wp - arr_hlp4(:) - arr_hlp2(:)), k_p_pool_depth_corr_min_occlud / k1_p_pool_depth_corr) 
  
      WHERE ((arr_hlp4(:) + arr_hlp2(:) + arr_hlp3(:)) > 1._wp)
        arr_hlp4(:) = arr_hlp4(:) / (arr_hlp4(:) + arr_hlp2(:) + arr_hlp3(:))
        arr_hlp2(:) = arr_hlp2(:) / (arr_hlp4(:) + arr_hlp2(:) + arr_hlp3(:))
        arr_hlp3(:) = arr_hlp3(:) / (arr_hlp4(:) + arr_hlp2(:) + arr_hlp3(:))
      ENDWHERE
  
      WHERE (isoil < soil_p_layer(:)) 
        po4_assoc_slow(:,isoil) = arr_hlp5(isoil) * (hlp1 + hlp2) / (hlp1 + hlp2 + hlp3 + hlp4) * &
                                  (1._wp - p_labile_slow_ratio(:,isoil))
        po4_occluded(:,isoil)   = arr_hlp5(isoil) * hlp3 / (hlp1 + hlp2 + hlp3 + hlp4)
        po4_primary(:,isoil)    = arr_hlp5(isoil) * hlp4 / (hlp1 + hlp2 + hlp3 + hlp4)
        p_labile(:,isoil)       = arr_hlp5(isoil) * (hlp1 + hlp2) / (hlp1 + hlp2 + hlp3 + hlp4) * p_labile_slow_ratio(:,isoil)
      ELSEWHERE
        po4_assoc_slow(:,isoil) = arr_hlp5(isoil) * arr_hlp2(:) * (1._wp - p_labile_slow_ratio(:,isoil))
        po4_occluded(:,isoil)   = arr_hlp5(isoil) * arr_hlp3(:)
        po4_primary(:,isoil)    = arr_hlp5(isoil) * arr_hlp4(:)
        p_labile(:,isoil)       = arr_hlp5(isoil) * arr_hlp2(:) * p_labile_slow_ratio(:,isoil)
      ENDWHERE
    ENDDO
  END SUBROUTINE init_sb_soil_p_pools_spp1685_sites

END MODULE mo_qs_process_init_util
