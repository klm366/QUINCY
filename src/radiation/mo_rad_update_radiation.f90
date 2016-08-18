!>
!! @brief Contains the routines for the shortwave surface radiation budget (fapar & albedo) 
!!
!! @author
!!  Sönke Zaehle, MPI-BGC Jena 
!!
!! @par Revision History
!! First version by Sönke Zaehle (2015-12-04)
!!
!! 18-08-2016: update to include VIS And NIR albedo calculation for the canopy
!!
MODULE mo_rad_update_radiation

  USE mo_kind,          ONLY: wp, eps8
  USE mo_rad_constants

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: update_radiation, calc_canopy_layer_depth

CONTAINS

  !-----------------------------------------------------------------------------------------------------
  ! Main Task
  ! 
  !-----------------------------------------------------------------------------------------------------
  !> calculates surface radiation reflection and absorption
  !! 
  !! @todo currently no soil and snow albedo calculation in the update_radiation routine
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE update_radiation(tile)

    USE mo_base_memory,      ONLY: t_slm_tile
    USE mo_base_constants,   ONLY: ncanop 
    USE mo_veg_memory,       ONLY: t_slm_plant
    USE mo_atmland_memory
    USE mo_atmland_constants,ONLY: rad2par
   
    IMPLICIT NONE
    ! 0.1 InOut
    TYPE(t_slm_tile),INTENT(inout), TARGET :: tile    !< one tile
    ! 0.2 Local
    TYPE(t_slm_plant),POINTER              :: plant   ! plant of one tile
    REAL(wp),DIMENSION(ncanop)             :: hlp1_cl
    REAL(wp)                               :: hlp1

    !  bi-directional pointer
    plant => tile%plant

    !> 1.0 initalise output
    !!
    !! atm2land land2atm are global TYPES, see mo_atmland_memory
    land2atm%alb_nir_dif = 0.0_wp
    land2atm%alb_nir_dir = 0.0_wp
    land2atm%alb_vis_dif = 0.0_wp
    land2atm%alb_vis_dir = 0.0_wp
    land2atm%albedo = 0.0_wp
    land2atm%net_radiation = 0.0_wp
    plant%apar = 0.0_wp

    !> 2.0 calculate canopy PAR radiation budget. 
    !! Note that the assumption is made that PAR albedo is identical to VIS albedo
    !! Note that the input is converted from W m-2 to µmol m-2 s-1 to have correct units for apar
    CALL calc_canopy_sw_radiation_budget(atm2land%swpar_srf_down*rad2par, & 
                     atm2land%swpar_fdiffuse, atm2land%cos_zenith_angle, & 
                     plant%lai, sigma_par, kbl0_par, kdf0_par, &
                     plant%canopy_layer%par_sunlit(:), & 
                     plant%canopy_layer%par_shaded(:), &
                     plant%canopy_layer%fleaf_sunlit(:), &
                     plant%apar, &
                     land2atm%alb_vis_dir,land2atm%alb_vis_dif)

    !> 3.0 calculate canopy NIR radiation budget. 
    !! Note that the assumption is made that diffuse and direct NIR are identical 
    CALL calc_canopy_sw_radiation_budget(atm2land%swnir_srf_down, & 
                     0.0_wp, atm2land%cos_zenith_angle, & 
                     plant%lai, sigma_nir, kbl0_nir, kdf0_nir, &
                     hlp1_cl,hlp1_cl,hlp1_cl,hlp1, & 
                     land2atm%alb_nir_dir,land2atm%alb_nir_dif)

    !> 4. calculate albedo for total shortwave radiation
    !!
    IF(atm2land%swvis_srf_down + atm2land%swnir_srf_down>eps8)THEN 
       land2atm%albedo = ( atm2land%swvis_srf_down * & 
                    ( land2atm%alb_vis_dir * (1._wp - atm2land%swpar_fdiffuse) + & 
                      land2atm%alb_vis_dif * atm2land%swpar_fdiffuse ) + & 
                      atm2land%swnir_srf_down * land2atm%alb_nir_dir ) / & 
                    (atm2land%swvis_srf_down + atm2land%swnir_srf_down) 
    ENDIF 

    !> 5. diagnose net radiation (temporarily here for diagnostic purposes)
    !!
    land2atm%net_radiation = (1._wp - land2atm%albedo) * (atm2land%swvis_srf_down + atm2land%swnir_srf_down)


  END SUBROUTINE update_radiation

  !-----------------------------------------------------------------------------------------------------
  !> calculate shortwave radiation budget of the canopy 
  !! 
  !! @par info on routine
  !!  This routine has been rewritten from ADFs implementation of the radiation calculations
  !!  by Spitters, C J T. 1986. \n “Separating the Diffuse and Direct Component of Global Radiation
  !!  and Its Implications for Modeling Canopy Photosynthesis .2. Calculation of Canopy Photosynthesis.”
  !!  Agricultural and Forest Meteorology 38 (1-3): 231–42. \n
  !!  Equations numbers refer to the Spitters paper 
  !!  The routine has been extended to diagnose albedo in the visible and NIR range. The implicit assumption 
  !!  is made that albedo in the PAR and VIS range are identical
  !! 
  !! @par Spitters Eq. 12 is wrong
  !!  accordingly to Wang2003, "radf4" was replaced by "-kbl" when calculating "idrdra"
  !! @par: APAR is calculated for all canopy layers irrespective of whether or not there are leaves 
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE calc_canopy_sw_radiation_budget(sw_srf_down,sw_fdiffuse,cos_zenith_angle,& 
                        lai,sigma,kbl0,kdf0,arad_sunlit,arad_shaded,fleaf_sunlit,arad,& 
                        alb_dir,alb_dif)

    USE mo_base_constants,     ONLY: ncanop,canopy_layer_depth
    USE mo_atmland_constants,  ONLY: min_cos_zenith_angle

    IMPLICIT NONE
    ! 0.1 InOut
    REAL(wp), INTENT(in) ::    sw_srf_down             !< shortwave downward radiation flux 
    REAL(wp), INTENT(in) ::    sw_fdiffuse             !< fraction of diffuse radiation (--)
    REAL(wp), INTENT(in) ::    cos_zenith_angle        !< solar angle (radians)
    REAL(wp), INTENT(in) ::    lai                     !< leaf area index (m-2 m-2)
    REAL(wp), INTENT(in) ::    sigma                   !< single leaf scattering coefficient 
    REAL(wp), INTENT(in) ::    kbl0                    !< extinction coefficient over black leaves 
    REAL(wp), INTENT(in) ::    kdf0                    !< extinction coefficient for diffuse radiation 
    REAL(wp), DIMENSION(ncanop), INTENT(out) :: arad_sunlit     !< absorbed radiation on sunlit leaves 
    REAL(wp), DIMENSION(ncanop), INTENT(out) :: arad_shaded     !< absorbed radiation of shaded leaves 
    REAL(wp), DIMENSION(ncanop), INTENT(out) :: fleaf_sunlit    !< fraction of leaves that is sunlit (--)
    REAL(wp), INTENT(out) ::   arad                   !< absorbed radiation 
    REAL(wp), INTENT(out) ::   alb_dir                !< albedo for direct light 
    REAL(wp), INTENT(out) ::   alb_dif                !< albedo for diffuse light
    ! 0.2 Local
    REAL(wp) :: irad_dir,irad_dif    ! downwelling direct and diffuse radiation 
    REAL(wp) :: rho                  ! fraction reflected by green canopy (--)
    REAL(wp) :: kdf,kbl              ! extinction coefficients for radiation over black canopy
    REAL(wp) :: lc                   ! cummulative LAI over mid point of a canopy layer
    REAL(wp) :: idfa                 ! diffuse radiation penetrating foliage at lc 
    REAL(wp) :: idra                 ! direct radiation penetrating sunlit foliage at lc 
    REAL(wp) :: idrdra               ! direct radiation remaining direct at lc 
    REAL(wp) :: temp1,temp2,radf1,radf2,radf3,radf3b,radf4 ! dummy variables
    REAL(wp), DIMENSION(ncanop) :: canopy_layer_depth_act  ! actual depth of canopy layers (account for incomplete lower bound) 
    INTEGER  :: icanop
    REAL(wp) :: tmp_cos_zenith_angle     ! local variable of solar angle (radians)
    REAL(wp) :: arad_dif,arad_dir        ! canopy integrated diffuse and direct absorbed radiation
    REAL(wp) :: irad_bc_dif,irad_bc_dir  ! bottom of the canopy light intensity

    arad_dif = 0.0    
    arad_dir = 0.0    
    arad_sunlit(:) = 0.0_wp 
    arad_shaded(:) = 0.0_wp
    fleaf_sunlit(:) = 0.0_wp
    tmp_cos_zenith_angle = cos_zenith_angle

    canopy_layer_depth_act = calc_canopy_layer_depth(lai)

    IF(cos_zenith_angle < min_cos_zenith_angle) tmp_cos_zenith_angle = min_cos_zenith_angle  ! allows negative zenith angle

    ! only do calculations if daylight is available, either direct or diffuse (in case of negative cos_zenith_angle)
    IF((sw_srf_down > eps8) .AND. lai .GT. eps8) THEN  

       ! Downwelling direct radiation 
       irad_dir = (1.0_wp - sw_fdiffuse) * sw_srf_down
       ! Downwelling diffuse radiation 
       irad_dif = sw_fdiffuse * sw_srf_down

       !> Eq.1 constants for albedo calculations
       radf3b = (1.0_wp - sigma)
       temp1 = SQRT(radf3b)
       temp2 = 2.0_wp * (1.0_wp - temp1)/(1.0_wp + temp1)

       !> Eq.1  Fraction of radiation reflected by green closed canopy (fraction)
       rho = temp2/(1.0_wp + rho2sbeta * tmp_cos_zenith_angle)

       ! Extinction coefficient for radiation over black leaves in canopy
       kdf = kdf0_par * temp1  ! temp1 == SQRT(1.0_wp - sigma)
       kbl = kbl0_par * kdf/(0.8_wp * temp1 * tmp_cos_zenith_angle)

       !> Eq.10 & 11 Some factors for radiation calculations that change with time
       radf1 = (1.0_wp - rho) * irad_dif * kdf            ! eq 10
       radf2 = (1.0_wp - rho) * temp1 * irad_dir * kbl    ! eq 11
       !> Eq.12 absorbed direct radiation
       radf3 = radf3b * irad_dir * kbl                    ! eq 12
       radf4 = -temp1 * kbl

       ! loop over canopy layers. All fluxes are evaluated at the centre of the canopy layer
       lc = canopy_layer_depth(1)/2.0_wp
       DO icanop = 1,ncanop
  
          IF (radf2 <= eps8) THEN ! case of only diffuse radiation
             ! Eq.10 Mean radiation penetrating foliage at lc 
             arad_shaded(icanop) = radf1 * exp(-kdf * lc)     ! eq 10
          ELSE ! Case of direct and diffuse radiation
             ! Sunlit foliage at lc (fraction)
             ! kbl: Extinction coefficient for radiation over black leaves in canopy
             fleaf_sunlit(icanop) = exp(-kbl * lc)
             ! Eq.10 Diffuse radiation penetrating foliage at lc 
             idfa = radf1 * exp(-kdf * lc)       
             ! Direct radiation penetrating sunlit foliage at lc 
             idra = radf2 * exp(radf4 * lc)
             ! Eq.12 Mean direct radiation remaining direct at lc 
             ! radf4 replaced by -kbl according to Wang2003 (error in Spitters 1986 Eqn.12)
             idrdra = radf3 * exp(-kbl * lc)  
             ! Radiation penetrating shaded foliage at lc (umol/m2/s)
             ! Eq.13 Sum of diffuse radiation and diffuse component of direct flux
             arad_shaded(icanop) = idfa + max((idra - idrdra), 0.0_wp)
             ! Radiation penetrating sunlit foliage at lc (umol/m2/s)
             ! sum of diffuse and direct radiation
             arad_sunlit(icanop) = arad_shaded(icanop) + radf3
          ENDIF 

          ! Integrate absorbed radiation over actually existing canopy layers
          IF(canopy_layer_depth_act(icanop)>eps8)THEN
             arad_dif = arad_dif + arad_shaded(icanop) * (1._wp-fleaf_sunlit(icanop)) * & 
                                   canopy_layer_depth_act(icanop)
             arad_dir = arad_dir + arad_sunlit(icanop) * fleaf_sunlit(icanop) * & 
                                   canopy_layer_depth_act(icanop)
          ENDIF

          IF(icanop < ncanop) & 
            lc = lc + canopy_layer_depth(icanop)/2.0_wp + canopy_layer_depth(icanop+1)/2.0_wp

       ENDDO

       ! canopy integrated incoming radiation that is absorbed
       arad = arad_dir+arad_dif

       ! light intensity at bottom of canopy
       irad_bc_dif = (1._wp - rho ) * irad_dif * exp(-kdf * lai)
       irad_bc_dir = (1._wp - rho ) * irad_dir * exp(radf4 * lai)

       ! albedo of the canopy, inferred from the difference between top and bottom of the 
       ! canopy light intensities minus the amount of absorbed radiation
       ! impossible to calculate a pure diffuse albedo, since direct radiation becomes diffuse
       ! such that the diffuse albedo is negative 
       alb_dir = ((irad_dir + irad_dif) - (irad_bc_dir + irad_bc_dif) - (arad_dir + arad_dif) ) / & 
                  (irad_dir + irad_dif)
       alb_dif = alb_dir 

    ENDIF
  END SUBROUTINE calc_canopy_sw_radiation_budget

  !-----------------------------------------------------------------------------------------------------
  !> calculate number of actual canopy layers present
  !! 
  !! 
  !-----------------------------------------------------------------------------------------------------
  FUNCTION calc_canopy_layer_depth(lai) RESULT(canopy_layer_depth_act)

    USE mo_base_constants,      ONLY: ncanop, canopy_layer_depth

    IMPLICIT NONE
    ! 0.1 InOut
    REAL(wp),INTENT(in)        :: lai
    REAL(wp),DIMENSION(ncanop) :: canopy_layer_depth_act
    ! 0.2 Local
    REAL(wp)  :: lc 
    INTEGER   :: icanop

    ! initialise output
    canopy_layer_depth_act(:) = 0._wp

    ! calculate number of actual canopy layers present
    IF(lai .GT. eps8) THEN ! avoid any calculation is there are no leaves 
       lc = 0._wp
       DO icanop = 1,ncanop
          IF(lc+canopy_layer_depth(icanop).LT.lai) THEN 
             canopy_layer_depth_act(icanop) = canopy_layer_depth(icanop)
          ELSE IF(lc.LT.lai.AND.lc+canopy_layer_depth(icanop).GE.lai) THEN
             ! lai = lc + canopy_layer_depth_act
             canopy_layer_depth_act(icanop) = lai - lc 
         ENDIF
         lc = lc + canopy_layer_depth_act(icanop)
       ENDDO
       ! make sure all LAI is represented by letting last layer go to the end of the canopy
       IF(lc < lai ) canopy_layer_depth_act(ncanop) = canopy_layer_depth_act(ncanop) + lai - lc
    ENDIF
  END FUNCTION calc_canopy_layer_depth

END MODULE mo_rad_update_radiation
