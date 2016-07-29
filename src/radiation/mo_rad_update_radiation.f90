!>
!! @brief Contains the routines for plant radiation (fapar & albedo) 
!!
!! @author
!!  Sönke Zaehle, MPI-BGC Jena 
!!
!! @par Revision History
!! First version by Sönke Zaehle (2015-12-04)
!! 2nd version by Jan Engel (20160719): svn move vegetation/mo_veg_radiation.f90 radiation/mo_rad_update_radiation.f90
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
  !> calculate surface radiation reflection and absorption in PAR band
  !! 
  !! @todo currently no albedo calculation in the update_radiation routine
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE update_radiation(tile)

    USE mo_base_memory,      ONLY: t_slm_tile
    USE mo_veg_memory,       ONLY: t_slm_plant
    USE mo_atmland_memory,   ONLY: atm2land, land2atm

    IMPLICIT NONE
    ! 0.1 InOut
    TYPE(t_slm_tile),INTENT(inout), TARGET :: tile    !< one tile
    ! 0.2 Local
    TYPE(t_slm_plant),POINTER              :: plant   ! plant of one tile

    !  bi-directional pointer
    plant => tile%plant

    !> 1.0 initalise output
    !!
    !! atm2land land2atm are global TYPES, see mo_atmland_memory
    land2atm%alb_nir_dif = 0.0_wp
    land2atm%alb_nir_dir = 0.0_wp
    land2atm%alb_vis_dif = 0.0_wp
    land2atm%alb_vis_dir = 0.0_wp
 
    !> 2.0 call subroutine
    !!
    !! 
    CALL calc_fapar(atm2land%swpar_srf_down, atm2land%swpar_fdiffuse, atm2land%cos_zenith_angle, & 
                     plant%lai, &
                     plant%canopy_layer(:)%par_sunlit, & 
                     plant%canopy_layer(:)%par_shaded, &
                     plant%canopy_layer(:)%fleaf_sunlit)
  END SUBROUTINE update_radiation

  !-----------------------------------------------------------------------------------------------------
  !> calculate number of actual canopy layers present
  !! 
  !! @par info on routine
  !!  This routine has been rewritten from ADFs implementation of the radiation calculations
  !!  by Spitters, C J T. 1986. \n “Separating the Diffuse and Direct Component of Global Radiation
  !!  and Its Implications for Modeling Canopy Photosynthesis .2. Calculation of Canopy Photosynthesis.”
  !!  Agricultural and Forest Meteorology 38 (1-3): 231–42. \n
  !!  Equations numbers refer to the Spitters paper 
  !! 
  !! @par Spitters Eq. 12 is wrong
  !!  accordingly to Wang2003, "radf4" was replaced by "-kbl" when calculating "idrdra"
  !! @todo: link canopy light calculation with albedo and fapar calculation. OBS: light is calculated 
  !!        for all canopy layers irrespective of whether or not there are leaves 
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE calc_fapar(swpar_srf_down,swpar_fdiffuse,cos_zenith_angle,& 
                        lai,par_sunlit,par_shaded,fleaf_sunlit)

    USE mo_base_constants,     ONLY: ncanop,canopy_layer_depth
    USE mo_atmland_constants,  ONLY: sw2par, min_cos_zenith_angle

    IMPLICIT NONE
    ! 0.1 InOut
    REAL(wp), INTENT(in) ::    swpar_srf_down          !< shortwave downward radiation flux (W m-2)
    REAL(wp), INTENT(in) ::    swpar_fdiffuse          !< fraction of diffuse radiation (--)
    REAL(wp), INTENT(in) ::    cos_zenith_angle        !< solar angle (radians)
    REAL(wp), INTENT(in) ::    lai                     !< leaf area index (m-2 m-2)
    REAL(wp), DIMENSION(ncanop), INTENT(out) :: par_sunlit     !< PAR on sunlit leaves (µmol m-2 s-1)
    REAL(wp), DIMENSION(ncanop), INTENT(out) :: par_shaded     !< PAR on shaded leaves (µmol m-2 s-1)
    REAL(wp), DIMENSION(ncanop), INTENT(out) :: fleaf_sunlit   !< fraction of leaves that is sunlit (--)
    ! 0.2 Local
    REAL(wp) :: ipar_dr,ipar_df      ! downwelling direct and diffuse PAR (µmol m-2 s-1)
    REAL(wp) :: rho                  ! fraction of PAR reflected by green canopy (--)
    REAL(wp) :: kdf,kbl              ! extinction coefficients for PAR over black canopy
    REAL(wp) :: lc                   ! cummulative LAI over mid point of a canopy layer
    REAL(wp) :: idfa                 ! diffuse PAR penetrating foliage at lc (umol/m2/s)
    REAL(wp) :: idra                 ! direct PAR penetrating sunlit foliage at lc (umol/m2/s)
    REAL(wp) :: idrdra               ! direct PAR remaining direct at lc (umol/m2/s) 
    REAL(wp) :: temp1,temp2,radf1,radf2,radf3,radf3b,radf4   ! dummy variables
    REAL(wp), DIMENSION(ncanop) :: canopy_layer_depth_act  ! actual depth of canopy layers (account for incomplete lower bound) 
    INTEGER  :: icanop
    REAL(wp) :: tmp_cos_zenith_angle   ! local variable of solar angle (radians)

    par_sunlit(:) = 0.0_wp 
    par_shaded(:) = 0.0_wp
    fleaf_sunlit(:) = 0.0_wp
    tmp_cos_zenith_angle = cos_zenith_angle

    canopy_layer_depth_act = calc_canopy_layer_depth(lai)

    IF(cos_zenith_angle < min_cos_zenith_angle) tmp_cos_zenith_angle = min_cos_zenith_angle  ! allows negative zenith angle

    ! only do calculations if daylight is available, either direct or diffuse (in case of negative cos_zenith_angle)
    IF((swpar_srf_down > eps8) .AND. lai .GT. eps8) THEN  

       ! Downwelling direct PAR umol/m2/s
       ipar_dr = (1.0_wp - swpar_fdiffuse) * sw2par * swpar_srf_down
       ! Downwelling diffuse PAR umol/m2/s
       ipar_df = swpar_fdiffuse * sw2par * swpar_srf_down

       !> Eq.1 Constants for visible albedo calculations
       radf3b = (1.0_wp - sigma)
       temp1 = SQRT(radf3b)
       temp2 = 2.0_wp * (1.0_wp - temp1)/(1.0_wp + temp1)

       !> Eq.1  Fraction of PAR reflected by green closed canopy (fraction)
       rho = temp2/(1.0_wp + rho2sbeta * tmp_cos_zenith_angle)

       ! Extinction coefficient for PAR over black leaves in canopy
       kdf = kdf0 * temp1  ! temp1 == SQRT(1.0_wp - sigma)
       kbl = kbl0 * kdf/(0.8_wp * temp1 * tmp_cos_zenith_angle)

       !> Eq.10 & 11 Some factors for radiation calculations that change with time
       radf1 = (1.0_wp - rho) * ipar_df * kdf            ! eq 10
       radf2 = (1.0_wp - rho) * temp1 * ipar_dr * kbl    ! eq 11
       !> Eq.12 absorbed direct radiation
       radf3 = radf3b * ipar_dr * kbl                    ! eq 12
       radf4 = -temp1 * kbl

       ! loop over canopy layers. All fluxes are evaluated at the centre of the canopy layer
       lc = canopy_layer_depth(1)/2.0_wp
       DO icanop = 1,ncanop
  
          IF (radf2 <= eps8) THEN ! case of only diffuse radiation
             !> Eq.10 Mean PAR penetrating foliage at lc (umol/m2/s)
             par_shaded(icanop) = radf1 * exp(-kdf * lc)     ! eq 10
          ELSE ! Case of direct and diffuse radiation
             ! Sunlit foliage at lc (fraction)
             ! kbl: Extinction coefficient for PAR over black leaves in canopy
             fleaf_sunlit(icanop) = exp(-kbl * lc)
             !> Eq.10 Diffuse PAR penetrating foliage at lc (umol/m2/s)
             idfa = radf1 * exp(-kdf * lc)       
             ! Direct PAR penetrating sunlit foliage at lc (umol/m2/s)
             idra = radf2 * exp(radf4 * lc)
             !> Eq.12 corrected following Wang2003 Mean direct PAR remaining direct at lc (umol/m2/s)
             idrdra = radf3 * exp(-kbl * lc)  ! radf4 replaced by -kbl according to Wang2003 (error in Spitters 1986 Eqn.12)
             ! PAR penetrating shaded foliage at lc (umol/m2/s)
             !> Eq.13 Sum of diffuse radiation and diffuse component of direct flux
             par_shaded(icanop) = idfa + max((idra - idrdra), 0.0_wp)
             ! PAR penetrating sunlit foliage at lc (umol/m2/s)
             ! sum of diffuse and direct radiation
             par_sunlit(icanop) = par_shaded(icanop) + radf3
          ENDIF 

          IF(icanop < ncanop) & 
            lc = lc + canopy_layer_depth(icanop)/2.0_wp + canopy_layer_depth(icanop+1)/2.0_wp

       ENDDO

    ENDIF
  END SUBROUTINE calc_fapar 

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
