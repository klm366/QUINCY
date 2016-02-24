!! @brief: Subroutines for computing tile specific turbulent exchange coefficients (tiled_srf_exchange_coeff),
!!         updating surface temperatures (surface_temp_implicit_flux_agg)
!!         and calculating dse and Q at the lowest atmospheric lvl (update_atm).
!!
!! @author Authors of ECHAM code
!! @Philipp de Vrese, MPI-M
!!
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
MODULE mo_flux_aggregation
  USE mo_io_units,          ONLY: nerr,nout
  USE mo_kind,              ONLY: wp
  USE mo_echam_vdiff_params,ONLY: ckap, cb,cc
  USE mo_physical_constants,ONLY: grav, rd, cpd, cpv, rd_o_cpd, rv,         &
    &                             vtmpc1,p0ref
#ifdef __ICON__
  USE mo_echam_phy_config,  ONLY: phy_config => echam_phy_config
#else
  USE mo_time_control,      ONLY: lstart
#endif
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: surface_temp_implicit_flux_agg, tiled_srf_exchange_coeff, update_atm, tiled_srf_exchange_coeff_jsbalone
CONTAINS
  !>
  !! Computation of the exchange coefficients at the Earth's terrestial surface with respect to different surface tiles (subroutine "tiled_srf_exchanged_coeff").
  !! Note that
  !! - the subroutine "tiled_srf_exchange_coeff" is a reduced version of the subroutine "sfc_exchange_coeff" contained in the module "mo_turbulence_diag".
  !! - the subroutine is called for each tile indvidually.
  !! @par Revision History
  !! Separated from vdiff.f90 of ECHAM6 and re-organized by Hui Wan (2010-09).
  !!  updated to echam-6.3.01 by Monika Esch (2014-11)
  !!  adapted to the tiled land surface by Philipp de Vrese (2016-01)
  SUBROUTINE tiled_srf_exchange_coeff(klon,                              &! in
                               & pz0m, paz0lh,                           &! in
                               & pcsat    , pcair,                       &! in 
                               & ptsfc    , pcpt_sfc,                    &! in
                               & pqsat_sfc,                              &! in      
                               & ppsfc,                                  &! in                         
                               & wind_b   ,                              &! in
                               & ptm1_b   , pgeom1_b,                    &! in
                               & pqm1_b   , pqxm1_b,                     &! in
                               & pqsat_b  , plh_b,                       &! in
                               & ptheta_b , pthetav_b,                   &! in
                               & pthetal_b, paclc_b,                     &! in
                               & pcfm_sfc,                               &! out
                               & pcfh_sfc                                &! out
                               & )

!!!!! de Vrese (01.12.2015) testing output - start  

    INTEGER, INTENT(IN) :: klon

    ! "_b" denotes value at the bottom level (the klev-th full level)

    REAL(wp),INTENT(IN) :: wind_b   (klon)  !< wind speed
    REAL(wp),INTENT(IN) :: ptm1_b   (klon)  !< temperature
    REAL(wp),INTENT(IN) :: pgeom1_b (klon)  !< geopotential
    REAL(wp),INTENT(IN) :: pqm1_b   (klon)  !< specific humidity
    REAL(wp),INTENT(IN) :: pqxm1_b  (klon)  !< total concentration of hydrometeors
    REAL(wp),INTENT(IN) :: pqsat_b  (klon)  !< saturation specific humidity
    REAL(wp),INTENT(IN) :: plh_b    (klon)  !< latent heat
    REAL(wp),INTENT(IN) :: ptheta_b (klon)  !< potential temp.
    REAL(wp),INTENT(IN) :: pthetav_b(klon)  !< virtual potential temp.
    REAL(wp),INTENT(IN) :: pthetal_b(klon)  !< liquid water (?) pot. temp.
    REAL(wp),INTENT(IN) :: paclc_b  (klon)  !< cloud cover at lowest model level

    ! "_sfc" denotes value at surface

    REAL(wp),INTENT(IN) :: pqsat_sfc (klon) !< saturation specific humidity at surface
    REAL(wp),INTENT(IN) :: pcpt_sfc  (klon) !< dry static energy

    REAL(wp),INTENT(IN) :: pz0m     (klon) !< aerodynamic roughness length
    REAL(wp),INTENT(IN) :: paz0lh (klon)  !< roughness length for heat over land
    
    REAL(wp),INTENT(IN) :: ptsfc    (klon) !< temp. at surface
    REAL(wp),INTENT(IN) :: ppsfc    (klon) !< surface pressure

    REAL(wp),INTENT(IN) :: pcsat  (klon)  !< area fraction with wet land surface
    REAL(wp),INTENT(IN) :: pcair  (klon)  !< area fraction with wet land surface (air)


    REAL(wp),INTENT(OUT) :: pcfm_sfc  (klon) !< exchange coeff. of momentum,
                                                      !< for each type of surface
    REAL(wp),INTENT(OUT) :: pcfh_sfc  (klon) !< exchange coeff. of heat and vapour
                                                      !< for each type of surface
    ! Local variables

    REAL(wp) :: pri_sfc   (klon) !< moist Richardson number

    REAL(wp) :: pchn_sfc (klon) !<
    REAL(wp) :: pcdn_sfc (klon) !<
    REAL(wp) :: pcfnc_sfc(klon) !<


    REAL(wp) :: zdu2   (klon) !<
    REAL(wp) :: zcfnch (klon) !<
    REAL(wp) :: zwst   (klon) !<
    REAL(wp) :: zqts   (klon)


    REAL(wp) :: zrdrv, zrvrd
    REAL(wp) :: z2b, z3b, z3bc, zepdu2
    REAL(wp) :: zqtl, zqtmit, zdqt, zqsmit
    REAL(wp) :: ztmit, ztheta, zthetav, zthetamit, zthetavmit, zfux, zfox
    REAL(wp) :: zmult1, zmult2, zmult3, zmult4, zmult5
    REAL(wp) :: zdus1, zdus2, zbuoy, zalo
    REAL(wp) :: zdthetal, zdthv, zscf, zucf, zucfh
    
    INTEGER  :: js

    !-------------------
    ! Some constants
    !-------------------
    zrvrd      = rv/rd
    zrdrv      = rd/rv
    z2b        = 2._wp*cb
    z3b        = 3._wp*cb
    z3bc       = 3._wp*cb*cc
    zepdu2     = 1._wp

    !------------------------------------------------------------------------------
    ! The prefactor (= air density * Rd) that will be multiplied to the exchange
    ! coefficients when building the linear algebraic equations. The density here is
    ! computed using air temperature of the lowest model level at time step n-1.
    !------------------------------------------------------------------------------
    !pprfac_sfc(1:klon) =  ppsfc(1:klon)                                     &
    !                     & /( ptm1_b(1:klon)                                  &
    !                     &   *(1._wp+vtmpc1*pqm1_b(1:klon)-pqxm1_b(1:klon)) )

    !-------------------------------------------------------------
    ! COMPUTATION OF BASIC QUANTITIES: WIND SHEAR,
    ! RICHARDSON NUMBER,SQUARED MIXING LENGTHS, UNSTABLE
    ! AND STABLE CASE COMMON FACTORS AND NEUTRAL CASE
    ! COMMON PART OF THE DRAG COEFFICIENTS.
    !-------------------------------------------------------------

    DO js = 1, klon

        zqts(js) = pcsat(js) * pqsat_sfc(js) + (1._wp - pcair(js))*pqm1_b(js) ! q_total at land surface

        ztheta      = ptsfc(js)*(p0ref/ppsfc(js))**rd_o_cpd
        zthetav     = ztheta*(1._wp+vtmpc1*zqts(js))
               
        zqtl       = pqm1_b(js) + pqxm1_b(js)                          ! q_total at lowest model level

        zqtmit     = 0.5_wp*( zqtl + zqts(js) )                   ! q_total, vertical average

        zqsmit     = 0.5_wp*( pqsat_b  (js) + pqsat_sfc(js) )  ! qs
        ztmit      = 0.5_wp*( ptm1_b   (js) + ptsfc    (js) )  ! temp.
        zthetamit  = 0.5_wp*( ptheta_b (js) + ztheta  )        ! potential temp.
        zthetavmit = 0.5_wp*( pthetav_b(js) + zthetav )        ! virtual potential temp.

        zfux = plh_b(js)/(cpd*ztmit)
        zfox = plh_b(js)/(rd*ztmit)

        zmult1 = 1._wp+vtmpc1*zqtmit   ! A in clear sky
        zmult2 = zfux*zmult1-zrvrd
        zmult3 = zrdrv*zfox*zqsmit/(1._wp+zrdrv*zfox*zfux*zqsmit)
        zmult5 = zmult1-zmult2*zmult3  ! A in cloud
        zmult4 = zfux*zmult5-1._wp     ! D in cloud

        zdus1 = paclc_b(js)*zmult5 + (1._wp-paclc_b(js))*zmult1   ! A avg
        zdus2 = paclc_b(js)*zmult4 + (1._wp-paclc_b(js))*vtmpc1   ! D avg

        zdqt     = zqtl - zqts(js)                           ! d qt
        zdthetal = pthetal_b(js) - ztheta                         ! d theta_l
!        zdu2(js) = MAX(zepdu2,(pum1_b(js)**2+pvm1_b(js)**2))
        zdu2(js) = MAX(zepdu2,(wind_b(js)**2))

        zbuoy        = zdus1*zdthetal + zdus2*zthetamit*zdqt
        pri_sfc(js)  = pgeom1_b(js)*zbuoy/(zthetavmit*zdu2(js))

        zalo = LOG( 1._wp + pgeom1_b(js)/(grav*pz0m(js)) )  ! ln[ 1 + zL/z0m ]
        pcdn_sfc(js) = (ckap/zalo)**2
        pcfnc_sfc(js)= SQRT(zdu2(js))*pcdn_sfc(js)
 !  REMARK:
 !  compared to the original (precalc_land,ocean,ice) the factor zcons is missing
 !  this factor is included as "prefactor for the exchange coefficients" in
 !  subroutine matrix_setup_elim
        zdthv         = MAX(0._wp,(zthetav-pthetav_b(js)))
        zwst(js) = zdthv*SQRT(zdu2(js))/zthetavmit

        pchn_sfc  (js) = (ckap / LOG(1._wp+pgeom1_b(js)/(grav*paz0lh(js))))**2
        zcfnch(js) = SQRT(zdu2(js))*pchn_sfc(js)
    END DO ! 1:klon
    !-------------------------------------------------------------------------
    ! Compute the exchange coefficients for momentum, heat and water vapour,
    ! for each type of surface
    !-------------------------------------------------------------------------

!TODO:    preset values to zero
     pcfh_sfc(1:klon) = 0._wp
     pcfm_sfc(1:klon) = 0._wp

      ! stable case              pri_sfc(js) > 0.
        DO js = 1,klon
          IF ( pri_sfc(js) > 0._wp ) THEN
            zscf = SQRT(1._wp+pri_sfc(js))
            pcfm_sfc(js)  = pcfnc_sfc (js)/(1._wp+z2b*pri_sfc(js)/zscf)   ! 5.2? 5.5
            pcfh_sfc(js)  = zcfnch    (js)/(1._wp+z2b*pri_sfc(js)*zscf)   ! 5.2? 5.6
          ENDIF
        ENDDO

      ! unstable case           pri_sfc(js) <= 0.
        DO js = 1,klon
          IF ( pri_sfc(js) <= 0._wp ) THEN
            zucf =  SQRT( -pri_sfc(js)*(1._wp+ pgeom1_b(js)/(grav*pz0m(js))) ) ! sqrt in (5.4)
            zucf =  1._wp+z3bc*pcdn_sfc(js)*zucf                   ! denominator in (5.4)
            zucf =  1._wp/zucf
            zucfh=  SQRT( -pri_sfc(js)*(1._wp+ pgeom1_b(js)/(grav*paz0lh(js))) ) ! sqrt in (5.4)
            zucfh=  1._wp+z3bc*pchn_sfc(js)*zucfh                   ! denominator in (5.4)
            zucfh=  1._wp/zucfh
            pcfm_sfc(js) = pcfnc_sfc(js)*(1._wp-z2b*pri_sfc(js)*zucf)  ! (5.2), (5.4)
            pcfh_sfc(js) = zcfnch   (js)*(1._wp-z3b*pri_sfc(js)*zucfh) ! (5.2), (5.4)
          ENDIF
        ENDDO
  END SUBROUTINE tiled_srf_exchange_coeff

!!!!! de Vrese (15.02.2016) Calculating surface drag coefficients according to the formulations used in jsbalone - start
  SUBROUTINE tiled_srf_exchange_coeff_jsbalone(klon,                              &! in
                               & pz0m, paz0lh,                           &! in
                               & pcsat    , pcair,                       &! in 
                               & ptsfc    ,                              &! in
                               & pqsat_sfc,                              &! in      
                               & ppsfc    ,                              &! in                         
                               & wind_b,                                 &! in
                               & ptm1_b,                                 &! in
                               & pqm1_b,                                 &! in
                               & pcfm_sfc,                               &! out
                               & pcfh_sfc                                &! out
                               & )

!!!!! de Vrese (01.12.2015) testing output - start  

    INTEGER, INTENT(IN) :: klon

    ! "_b" denotes value at the bottom level (the klev-th full level)
    REAL(wp),INTENT(IN) :: wind_b   (klon)  !< wind speed
    REAL(wp),INTENT(IN) :: ptm1_b   (klon)  !< temperature
    REAL(wp),INTENT(IN) :: pqm1_b   (klon)  !< specific humidity

    ! "_sfc" denotes value at surface
    REAL(wp),INTENT(IN) :: pqsat_sfc (klon) !< saturation specific humidity at surface
    REAL(wp),INTENT(IN) :: pz0m     (klon) !< aerodynamic roughness length
    REAL(wp),INTENT(IN) :: paz0lh (klon)  !< roughness length for heat over land    
    REAL(wp),INTENT(IN) :: ptsfc    (klon) !< temp. at surface
    REAL(wp),INTENT(IN) :: ppsfc    (klon) !< surface pressure
    REAL(wp),INTENT(IN) :: pcsat  (klon)  !< area fraction with wet land surface
    REAL(wp),INTENT(IN) :: pcair  (klon)  !< area fraction with wet land surface (air)

    REAL(wp),INTENT(OUT) :: pcfm_sfc  (klon) !< exchange coeff. of momentum,
                                                      !< for each type of surface
    REAL(wp),INTENT(OUT) :: pcfh_sfc  (klon) !< exchange coeff. of heat and vapour
                                                      !< for each type of surface


    ! Local variables

    REAL(wp) :: zdu2   (klon) !<
    REAL(wp) :: zcfnch (klon) !<
    REAL(wp) :: zwst   (klon) !<
    REAL(wp) :: zqts   (klon)
    REAL(wp) :: ztvd   (klon)
    REAL(wp) :: ztvs   (klon)
    REAL(wp) :: zg     (klon)
    REAL(wp) :: zril   (klon)
    REAL(wp) :: zcdnl  (klon)
    REAL(wp) :: zchnl  (klon)
    REAL(wp) :: zscfl  (klon)
    REAL(wp) :: zucfl  (klon)
    REAL(wp) :: zucfhl (klon)
    REAL(wp) :: zcfncl (klon)
    REAL(wp) :: zcfnchl(klon)
    REAL(wp) :: zcfml  (klon)
    REAL(wp) :: zcfhl  (klon)

    REAL(wp) :: z2b, z3b, z3bc

    !-------------------
    ! Some constants
    !-------------------
    z2b        = 2._wp*cb
    z3b        = 3._wp*cb
    z3bc       = 3._wp*cb*cc
    !zcons12       = cvdifts * time_step_len * Gravity / GasConstantDryAir

 !------------------------------------------------------------------------------------
    ! Approximation of cdrag
    !------------------------------------------------------------------------------------
    ! squared wind shear ! minimum wind speed square from echam/mo_surface_land.f90:precalc_land
    zdu2(:) = MAX(wind_b(:)**2,1._wp)   

    ! virtual potential air temperature (see mo_surface_boundary.f90)
    ! according to Saucier, WJ Principles of Meteoroligical Analyses
    ! tv = t * (1 + 0.61 * q * t)    ! virtual temperature
    ! td = t * ( 1000 / p_mb ) ^ R/cdp  ! potential temperature
    ! tvd = tair * (100000/p_pa)^zkappa * 1 + vtmpc1 * qair) ! virtual potential temperature
 !------------------------------------------------------------------------------------
    ztvd(:) = ( ptm1_b(:) * ( p0ref/ppsfc(:))**rd_o_cpd ) * & 
              ( 1._wp + vtmpc1 * pqm1_b(:) )
              
    ! virtual potential surface temperature 
    ztvs = ptsfc(:) * ( p0ref/ppsfc(:))**rd_o_cpd * &
           ( 1._wp + vtmpc1 * ( pcsat(:) * pqsat_sfc(:) + ( 1._wp - pcair(:) ) * pqm1_b(:)))

    ! geopotential of the surface layer (see echam's auxhybc.f90 & geopot.f90)
    ! If HeightWind is set, then the measurement height + an offset from the average vegetation height
    ! is used. Otherwise, the code defaults to the half-level of ECHAM's lowest atmospheric layer
    !IF(HeightWind > 0._dp) THEN
    !   zg(:) = (0.75_dp * Drv_Var%veg_height + HeightWind) * Gravity
    !ELSE
    
    zg(:) = ztvd(:) *  rd * LOG(1._wp / ( 1._wp - 0.0025_wp ))
    !ENDIF
    ! Richardson number (dry, Brinkop & Roeckner 1995, Tellus)
    ! ztvd, ztvs are now virtual potential temperatures, changed by Thomas Raddatz 07.2014
    zril(:) = zg(:) * ( ztvd(:) - ztvs(:) ) / ( zdu2(:) * ztvd(:))
    zril(:) = MAX(MIN(zril(:),5._wp),-5._wp)
    
    ! Neutral drag coefficient for momentum and heat
    zcdnl(:) = (ckap / LOG(1._wp + zg(:) / (grav * pz0m(:)  )))**2
    zchnl(:) = (ckap / LOG(1._wp + zg(:) / (grav * paz0lh(:))))**2

    ! account for stable/unstable case: helper variables
    zscfl(:) = SQRT (  1._wp + ABS(zril(:)))
    
    zucfl(:) = 1._wp / (1._wp + z3bc * zcdnl(:) * SQRT(ABS(zril(:)) * (1._wp  &
            + zg(:) / (grav * pz0m(:)))))
            
    zucfhl(:) = 1._wp / (1._wp + z3bc * zchnl(:) * SQRT(ABS(zril(:)) * (1._wp  &
            + zg(:) / (grav * paz0lh(:)))))

    ! ignoring cloud water correction (see mo_surface_land.f90)
    
    !  REMARK:
    !  compared to the original (precalc_land,ocean,ice) the factor zcons is missing
    !  this factor is included as "prefactor for the exchange coefficients" in
    !  subroutine matrix_setup_elim
    !------------------------------------------------------------------------------------
    !zcons(:) = zcons12 * pressure(:) / ( ptm1_b(:) * (1._dp + vtmpc1 * pqm1_b(:)))
    !zcfncl(:)  = zcons(:) * SQRT(zdu2(:)) * zcdnl(:)
    !zcfnchl(:)  = zcons(:) * SQRT(zdu2(:)) * zchnl(:)
    !------------------------------------------------------------------------------------
    zcfncl(:)  = SQRT(zdu2(:)) * zcdnl(:)
    zcfnchl(:) = SQRT(zdu2(:)) * zchnl(:)

    ! Stable / Unstable case
    WHERE ( zril(:) .GT. 0._wp )
       zcfml(:) = zcfncl(:)  / (1._wp + z2b * zril(:) * zscfl(:))
       zcfhl(:) = zcfnchl(:) / (1._wp + z2b * zril(:) * zscfl(:))
    ELSEWHERE
       zcfml(:) = zcfncl(:)  * (1._wp - z2b * zril(:) * zucfl(:))
       zcfhl(:) = zcfnchl(:) * (1._wp - z3b * zril(:) * zucfhl(:))
    ENDWHERE
    
! Standard formulation in jsbalone is: 
!    Drv_Var%cdrag(1:kland) = zcfhl(:)
!    Drv_Var%zchl(1:kland)  = zcfhl(:) / zcfnchl(:) * zchnl(:)
! But here cdrag corresponds to pcfh_sfc (exchange coeff. of heat and vapour) as in mo_turbulence_diag.f90 hence:
    pcfm_sfc(:)  = zcfml(:)
    pcfh_sfc(:)  = zcfhl(:)
  END SUBROUTINE tiled_srf_exchange_coeff_jsbalone
!!!!! de Vrese (15.02.2016) Calculating surface drag coefficients according to the formulations used in jsbalone - end 

!!!!! de Vrese (01.12.2015) Calculating new state based on flux aggregation see (M. J. Best, A. Beljaars, J. Polcher, and P. Viterbo, 2004) - start
  SUBROUTINE surface_temp_implicit_flux_agg(klon,              &
    &            steplen,                                      &
    &            ntiles,                                       &
    &            cfract_srf,                                   &
    &            t2s_conv,                                     &
    &            t_acoef_flux, t_bcoef_flux,                   &
    &            q_acoef_flux, q_bcoef_flux,                   &
    &            t_bcoef, q_bcoef,                             &
    &            s_srf_old, qsat_srf_old, dQdT,                &
    &            rad_srf_net, grnd_hflx,                       &
    &            heat_tcoef, fact_q_air, fact_qsat_srf,        &
!!!!! de Vrese (02.23.2016) Calculating Heat consumption by melting snow and land ice [W/m**2] - start  
    &            snmel_hflx,                                   &
!!!!! de Vrese (02.23.2016) Calculating Heat consumption by melting snow and land ice [W/m**2] - start  
    &            sfract_srf, hcap_grnd,                        &
    &            s_srf_new,                                    &
    &            t_acoef_equiv_tile, q_acoef_equiv_tile,       &
    &            t_acoef_equiv, q_acoef_equiv                  &
    &            )

!!!!! ------------------- NOTES ----------------------- !!!!!


!------------------  Surface energy balance euqation is solved using the flux aggregation technique discribed in:  ------------------
!  M. J. Best, A. Beljaars, J. Polcher, and P. Viterbo, 2004:
!  A Proposed Structure for Coupling Tiled Surfaces with the Planetary Boundary Layer. J. Hydrometeor, 5, 1271–1278.


!#ifndef STANDALONE
!    USE mo_physc2, ONLY: cvdifts
!#endif

    USE mo_echam_vdiff_params, ONLY: cvdifts   ! TBD: change to 1. for standalone JSBACH
    USE mo_physical_constants, ONLY: &
    & cpd, &
    & zemiss_def, &
    & stbo, &
    & alv,  &
    & als




    IMPLICIT NONE
    INTEGER,  INTENT(in)        :: klon
    INTEGER,  INTENT(in)        :: ntiles
!    REAL(wp), INTENT(in),  dimension(klon,ntiles):: cfract_srf
    REAL(wp), INTENT(in)        :: steplen
    REAL(wp), INTENT(in),  dimension(klon)       :: t_acoef_flux, t_bcoef_flux, q_acoef_flux, q_bcoef_flux
    REAL(wp), INTENT(IN),  dimension(klon)       :: q_bcoef, t_bcoef
    REAL(wp), INTENT(in),  dimension(klon,ntiles):: t2s_conv, s_srf_old, qsat_srf_old, dQdT, heat_tcoef, fact_q_air, fact_qsat_srf
    REAL(wp), INTENT(in),  dimension(klon,ntiles):: cfract_srf, rad_srf_net, grnd_hflx, hcap_grnd, sfract_srf, snmel_hflx
    REAL(wp), INTENT(out), dimension(klon,ntiles):: s_srf_new
    REAL(wp), INTENT(OUT), dimension(klon,ntiles):: t_acoef_equiv_tile, q_acoef_equiv_tile
    REAL(wp), INTENT(OUT), dimension(klon)       :: t_acoef_equiv, q_acoef_equiv
!   LOCAL VARIABLES
    INTEGER :: i,j
!    REAL(wp) :: zicp
    REAL(wp)                        :: pdt, pemi, pboltz
    REAL(wp)                        :: ztpfac1, ztmst !, s_atm_new, q_atm_new
    REAL(wp), dimension(klon)       :: DH1_av, DH2_av, DH4_av, DE1_av, DE2_av, DE4_av, M1, M2
    REAL(wp), dimension(klon,ntiles):: zicp, lair, lsat
    REAL(wp), dimension(klon,ntiles):: CH1, CH2, CH3, CH4, CE1, CE2, CE3, CE4
    REAL(wp), dimension(klon,ntiles):: DS1, DS2, DS4, Z, DH1, DH2, DH4, DE1, DE2, DE4 

    REAL(wp), dimension(klon,ntiles):: qsat_srf_new
    REAL(wp), dimension(klon)       :: s_atm_new, q_atm_new

#ifdef STANDALONE
    REAL(wp), PARAMETER :: cvdifts = 1.0_wp
#endif  
!-------------------------------------------------------------------------------------
!   Constants
    ztpfac1 = cvdifts
    ztmst   = steplen
    pdt     = ztpfac1*ztmst                  ! zcons29 in 'old' vdiff
    pemi    = zemiss_def                     ! emissivity
    pboltz  = stbo

!************************************************************************************

    zicp(:,:) = 1._wp / t2s_conv(:,:)
  
    do j = 1, klon
    do i = 1, ntiles
      if (fact_q_air(j,i) == 0.0_wp) Then
        lair(j,i) = 0.0_wp
      else
        lair(j,i) = sfract_srf(j,i)/fact_q_air(j,i) * (als - alv) + alv
      end if 
      if (fact_qsat_srf(j,i)== 0.0_wp) Then
        lsat(j,i) = 0.0_wp
      else 
        lsat(j,i) = sfract_srf(j,i)/fact_qsat_srf(j,i) * (als - alv) + alv
      end if 
    end do
    end do
!------------------  Coupling schemeas described in Best et al. (2004)  ------------------

!Calculating coefficients to linearize sensible flux, moisture flux and surface energy balance
! H = SA * CH1 + QA * CH2 + SS * CH3 + CH4
! E = SA * CE1 + QA * CE2 + SS * CE3 + CE4
! WITH H  = SENSIBLE HEAT FLUX AND E = EVAPOTRANSPIRATION
! ,    SA = S_ATM  OR  QA = QATM : ATM. VALUE OF S OR Q
! AND  SS = SSURF  SURFACE VALUE OF S 
!

    CH1(:,:) =   heat_tcoef(:,:)
    CH2(:,:) =   0.0_wp
    CH3(:,:) = - heat_tcoef(:,:)
    CH4(:,:) =   0.0_wp
 
    CE1(:,:) =   0.0_wp
    CE2(:,:) =   heat_tcoef(:,:) * fact_q_air(:,:)
    CE3(:,:) = - heat_tcoef(:,:) * fact_qsat_srf(:,:) * dQdT(:,:) * zicp(:,:)
    CE4(:,:) = - heat_tcoef(:,:) * fact_qsat_srf(:,:) &
             &* (qsat_srf_old(:,:) - dQdT(:,:)*zicp(:,:)*s_srf_old(:,:))
 
! Solving the surface energy balance for SS and linearizing it so that 
! SS = DS1 * SA + DS2 * QA + DS4

    Z(:,:) = (hcap_grnd(:,:) * zicp(:,:) / pdt + 4.0_wp * pemi*pboltz*zicp(:,:) *(zicp(:,:) * s_srf_old(:,:))**3.0_wp&
           &- CH3(:,:) - lsat(:,:) * CE3(:,:))
 
    DS1(:,:) =  CH1(:,:) / Z(:,:)
    DS2(:,:) = (CH2(:,:) + lair(:,:) * CE2(:,:)) / Z(:,:)
    DS4(:,:) = (rad_srf_net(:,:) + 4.0_wp * pemi * pboltz *(zicp(:,:)*s_srf_old(:,:))**4.0_wp + CH4(:,:) &
!!!!! de Vrese (02.23.2016) Using Heat consumption by melting snow and land ice [W/m**2] - start
    !        & +  lsat(:,:)*CE4(:,:) + grnd_hflx(:,:) &
             & +  lsat(:,:)*CE4(:,:) + grnd_hflx(:,:) + snmel_hflx(:,:)&
!!!!! de Vrese (02.23.2016) Using Heat consumption by melting snow and land ice [W/m**2] - start  
             & +  s_srf_old(:,:) * zicp(:,:) * hcap_grnd(:,:) / pdt) / Z(:,:)
 
! writing the sensible heat flux eqaution and replacing SS with the coefficients DS1 - 4 and linearizing it so that 
! H = DH1 * SA + DH2 * QA + DH4
 
    DH1(:,:) = CH1(:,:) + CH3(:,:) * DS1(:,:)
    DH2(:,:) = CH2(:,:) + CH3(:,:) * DS2(:,:)
    DH4(:,:) = CH4(:,:) + CH3(:,:) * DS4(:,:)
 
! Summing up cfract_srfi * DH1i, DH2i ...
 
    DH1_av = 0._wp
    DH2_av = 0._wp
    DH4_av = 0._wp
 
    do j = 1, klon 
    DH1_av(j) = sum(DH1(j,:) * cfract_srf(j,:))
    DH2_av(j) = sum(DH2(j,:) * cfract_srf(j,:))
    DH4_av(j) = sum(DH4(j,:) * cfract_srf(j,:))
    end do
! writing the moisture flux eqaution and replacing SS with the coefficients DS1 - 4 and linearizing it so that 
! E = DE1 * SA + DE2 * QA + DE4
 
    DE1(:,:) = CE1(:,:) + CE3(:,:) * DS1(:,:)
    DE2(:,:) = CE2(:,:) + CE3(:,:) * DS2(:,:)
    DE4(:,:) = CE4(:,:) + CE3(:,:) * DS4(:,:)
 
! Summing up cfract_srfi * DE1i, DE2i ...
 
    DE1_av = 0._wp
    DE2_av = 0._wp
    DE4_av = 0._wp

    do j = 1, klon 
    DE1_av(j) = sum(DE1(j,:) * cfract_srf(j,:))
    DE2_av(j) = sum(DE2(j,:) * cfract_srf(j,:))
    DE4_av(j) = sum(DE4(j,:) * cfract_srf(j,:))
    end do
!  substituting H* with "(SA - BS) / AS" and E* by "( QA - BQ ) / AQ" yields:
!  (SA - BS) / AS = DH1_av * SA + DH2_av * QA + DH4_av and
!  (QA - BQ) / AQ = DE1_av * SA + DE2_av * QA + DE4_av
!  WITH H*  = AGGREGATED SENSIBLE HEAT FLUX AND E* = AGGREGATED EVAPOTRANSPIRATION
!  AND AS, BS, AQ AND BQ = COEFFICIENTS FROM THE RICHTMYER AND MORTON SCHEME:
!  QA = AQ * E* + BQ
!  SA = AS * H* + BS

    s_atm_new(:)     =  (DH2_av(:) * (DE4_av(:) + q_bcoef_flux(:) / q_acoef_flux(:))/(1._wp / q_acoef_flux(:) - DE2_av(:)) &
                     & + DH4_av(:)              + t_bcoef_flux(:) / t_acoef_flux(:))&
                     & /(1._wp / t_acoef_flux - DH1_av - DH2_av * DE1_av / (1._wp / q_acoef_flux - DE2_av))

    q_atm_new(:)     = (s_atm_new(:) * DE1_av(:) + DE4_av(:) + q_bcoef_flux(:)/q_acoef_flux(:))/(1.0_wp/q_acoef_flux(:) - DE2_av(:))
 
    do j = 1, klon
    s_srf_new(j,:)   = DS1(j,:) * s_atm_new(j) + DS2(j,:) *q_atm_new(j) + DS4(j,:)
    qsat_srf_new(j,:) = qsat_srf_old(j,:) + dQdT(j,:)*zicp(j,:)*s_srf_new(j,:) - dQdT(j,:)*zicp(j,:)*s_srf_old(j,:)
!------------------  An equivalent RM coefficient:  ------------------
! X_Acoeff = (XA - X_BCOEFF) / X_SURF
! THE COEFFICIENTS FOR INDIVIDUAL TILES ARE REQUIRED TO CALCULATE ATMOSPHERIC STATE IN UPDATE_SURFACE_FLUXES
    t_acoef_equiv_tile(j,:) = (s_atm_new(j)-t_bcoef(j))/s_srf_new(j,:)
    q_acoef_equiv_tile(j,:) = (q_atm_new(j)-q_bcoef(j))/qsat_srf_new(j,:)
! NOTE THAT THE GRID AGGREGATED COEFFICIENTS NEED TO BE PASSED BACK TO ICON TO CALCULATE STATE ON THE LOWEST LVL IN VDIFF_UP
    t_acoef_equiv(j) = (s_atm_new(j)-t_bcoef(j))/SUM(cfract_srf(j,:)*s_srf_new(j,:))
    q_acoef_equiv(j) = (q_atm_new(j)-q_bcoef(j))/SUM(cfract_srf(j,:)*qsat_srf_new(j,:))
    end do

  END SUBROUTINE surface_temp_implicit_flux_agg
!!!!! de Vrese (01.12.2015) Calculating new state based on flux aggregation see (M. J. Best, A. Beljaars, J. Polcher, and P. Viterbo, 2004) - end

!!!!! de Vrese (01.12.2015) Apdates state of the lowest atm model level, this routine may be needed, to recalculate a grid box mean state  - start
  SUBROUTINE update_atm(klon,              &
    &            ntiles,                                       &
    &            qsat_tile,                                    &
    &            s_srf_flux,                                   &
    &            t_acoef_flux, t_bcoef_flux,                   &
    &            q_acoef_flux, q_bcoef_flux,                   &
    &            heat_tcoef_tile,                              &
!    &            cfract_srf,                                   &
    &            fact_qsat_srf,                                &
    &            fact_q_air,                                   &
    &            q_air_new,                                    &
    &            s_air_new                                     &
    &            )

!!!!! ------------------- NOTES ----------------------- !!!!!


!------------------  DSE and Q at the lowest level are calculated by:  ------------------
!  Q_A = A_Q * E* + B_Q
!  S_A = A_S * H* + B_S
    IMPLICIT NONE
    INTEGER,  INTENT(in)                             ::  klon
    INTEGER,  INTENT(in)                             ::  ntiles
    REAL(wp), INTENT(in),            dimension(klon) ::  t_acoef_flux, t_bcoef_flux, &
                                                     &   q_acoef_flux, q_bcoef_flux
                                          
    REAL(wp), INTENT(in),     dimension(klon,ntiles) ::  qsat_tile, s_srf_flux, &
                                                     &   heat_tcoef_tile,       &
                                                     &   fact_qsat_srf, fact_q_air
!   REAL(wp), INTENT(in),  dimension(klon,ntiles)      :: cfract_srf

    REAL(wp), INTENT(out),           dimension(klon) ::  q_air_new, s_air_new
    
    INTEGER :: j
    REAL(wp),                 dimension(klon,ntiles) ::  cfract_srf

    cfract_srf=1._wp

    do j = 1, klon
    !!!!! Humidity
      q_air_new(j) =(q_bcoef_flux(j)  -q_acoef_flux(j)*SUM(cfract_srf(j,:)*heat_tcoef_tile(j,:)*fact_qsat_srf(j,:)*qsat_tile(j,:))) &
                               /(1._wp-q_acoef_flux(j)*SUM(cfract_srf(j,:)*heat_tcoef_tile(j,:)*fact_q_air(j,:)))
    !!!!! Dry static energy    
      s_air_new(j) =(t_bcoef_flux(j)  -t_acoef_flux(j)*SUM(cfract_srf(j,:)*heat_tcoef_tile(j,:)*s_srf_flux(j,:))) &
                               /(1._wp-t_acoef_flux(j)*SUM(cfract_srf(j,:)*heat_tcoef_tile(j,:)))
    end do

  END SUBROUTINE update_atm
!-------------
!!!!! de Vrese (01.12.2015) Apdates state of the lowest atm model level, this routine may be needed, to recalculate a grid box mean state  - end
END MODULE mo_flux_aggregation
