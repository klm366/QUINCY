!>
!! @brief Contains the routines for the srf processes
!!
!! @author
!!  Reiner Schnur, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Reiner Schnur (2013-04-04)
!!
MODULE mo_srf_process
#ifndef __NO_JSBACH__

  USE mo_kind,      ONLY: wp
  USE mo_exception, ONLY: message

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: surface_temp_implicit
!!!!! de Vrese (01.12.2015) Calculating new state based on flux aggregation see (M. J. Best, A. Beljaars, J. Polcher, and P. Viterbo, 2004) - start
  PUBLIC :: surface_temp_implicit_flux_agg
!!!!! de Vrese (01.12.2015) Calculating new state based on flux aggregation see (M. J. Best, A. Beljaars, J. Polcher, and P. Viterbo, 2004) - end
  CHARACTER(len=*), PARAMETER :: modname = 'mo_srf_process'

CONTAINS

  SUBROUTINE surface_temp_implicit(klon, steplen,             &
    &            pcp,                                         &
    &            pescoe, pfscoe, peqcoe, pfqcoe,              &
    &            psold, pqsold, pdqsold,                      &
    &            pnetrad, pgrdfl,                             &
    &            pcfh, pcair, pcsat, pfracsu, pgrdcap,        &
    &            psnew)

    USE mo_echam_vdiff_params, ONLY: cvdifts   ! TBD: change to 1. for standalone JSBACH
    USE mo_physical_constants, ONLY: &
      & cpd, &
      & vtmpc2, &
      & zemiss_def, &
      & stbo, &
      & alv,  &
      & als


    INTEGER,  INTENT(in)    :: klon
    REAL(wp), INTENT(in)    :: steplen
    REAL(wp), INTENT(in)    :: pcp(klon), pfscoe(klon), pescoe(klon), pfqcoe(klon), peqcoe(klon)
    REAL(wp), INTENT(in)    :: psold(klon), pqsold(klon), pdqsold(klon)
    REAL(wp), INTENT(in)    :: pnetrad(klon), pgrdfl(klon)
    REAL(wp), INTENT(in)    :: pcfh(klon), pcair(klon), pcsat(klon), pfracsu(klon)
    REAL(wp), INTENT(in)    :: pgrdcap(klon)
    REAL(wp), INTENT(out)   :: psnew(klon)
    !
    REAL(wp) :: zcolin(klon), zcohfl(klon), zcoind(klon), zicp(klon), zca(klon), zcs(klon)
    REAL(wp) :: pdt
    REAL(wp) :: pc16

    pdt     = cvdifts * steplen
    pc16    = cpd * vtmpc2

    !
    zicp(:) = 1._wp / pcp(:)
    !
    zca(:)    = als * pfracsu(:) +  alv * (pcair(:) - pfracsu(:))
    zcs(:)    = als * pfracsu(:) +  alv * (pcsat(:) - pfracsu(:))
    !
    zcolin(:) = pgrdcap(:)*zicp(:) +                                                              &
              &      pdt * (zicp(:) * 4._wp * zemiss_def * stbo * ((zicp(:) * psold(:))**3)  -    &
              &      pcfh(:) * (zca(:) * peqcoe(:) - zcs(:)) *                                    &
              &      zicp(:) * pdqsold(:))
    !
    zcohfl(:) = -pdt * pcfh(:) * (pescoe(:) - 1._wp)
    !
    zcoind(:) = pdt * (pnetrad(:) + pcfh(:) * pfscoe(:) +  pcfh(:) *                              &
              &      ((zca(:) * peqcoe(:) - zcs(:)) * pqsold(:) + zca(:) * pfqcoe(:)) + pgrdfl(:))
    !
    psnew(:)  = (zcolin(:) * psold(:) + zcoind(:)) / (zcolin(:) + zcohfl(:))

  END SUBROUTINE surface_temp_implicit

!!!!! de Vrese (01.12.2015) Calculating new state based on flux aggregation see (M. J. Best, A. Beljaars, J. Polcher, and P. Viterbo, 2004) - start
  SUBROUTINE surface_temp_implicit_flux_agg(klon,              &
    &            ntiles,                                       &
    &            steplen,                                      &
    &            t2s_conv,                                     &
    &            t_acoef_flux, t_bcoef_flux,                   &
    &            q_acoef_flux, q_bcoef_flux,                   &
    &            s_srf_old, qsat_srf_old, dQdT,                &
    &            rad_srf_net, grnd_hflx,                       &
    &            heat_tcoef, fact_q_air, fact_qsat_srf,        &
!    &            cfract_srf,                                   &
    &            sfract_srf, hcap_grnd,                        &
    &            s_srf_new)!, qsat_srf_new, s_atm_new, q_atm_new)



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
    REAL(wp), INTENT(in)       :: steplen
    REAL(wp), INTENT(in),  dimension(klon)       :: t_acoef_flux, t_bcoef_flux, q_acoef_flux, q_bcoef_flux
    REAL(wp), INTENT(in),  dimension(klon,ntiles):: t2s_conv, s_srf_old, qsat_srf_old, dQdT, heat_tcoef, fact_q_air, fact_qsat_srf
!    REAL(wp), INTENT(in),  dimension(klon,ntiles):: cfract_srf, rad_srf_net, grnd_hflx, hcap_grnd, sfract_srf
    REAL(wp), INTENT(in),  dimension(klon,ntiles):: rad_srf_net, grnd_hflx, hcap_grnd, sfract_srf
    REAL(wp), INTENT(out), dimension(klon,ntiles):: s_srf_new    
!    REAL(wp), INTENT(out), dimension(klon,ntiles):: s_srf_new, qsat_srf_new
!    REAL(wp), INTENT(out), dimension(klon)       :: s_atm_new, q_atm_new
!   LOCAL VARIABLES
    INTEGER :: i,j
!    REAL(wp) :: zicp
    REAL(wp)                        :: pdt, pemi, pboltz
    REAL(wp)                        :: ztpfac1, ztmst !, s_atm_new, q_atm_new
    REAL(wp), dimension(klon)       :: DH1_av, DH2_av, DH4_av, DE1_av, DE2_av, DE4_av, M1, M2
    REAL(wp), dimension(klon,ntiles):: zicp, lair, lsat
    REAL(wp), dimension(klon,ntiles):: CH1, CH2, CH3, CH4, CE1, CE2, CE3, CE4
    REAL(wp), dimension(klon,ntiles):: DS1, DS2, DS4, Z, DH1, DH2, DH4, DE1, DE2, DE4 

    REAL(wp), dimension(klon,ntiles):: cfract_srf
    REAL(wp), dimension(klon,ntiles):: qsat_srf_new
    REAL(wp), dimension(klon)       :: s_atm_new, q_atm_new

#ifdef STANDALONE
    REAL(wp), PARAMETER :: cvdifts = 1.0_wp
#endif


    cfract_srf=1._wp 
  
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
             & +  lsat(:,:)*CE4(:,:) + grnd_hflx(:,:) &
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
    end do
 
  END SUBROUTINE surface_temp_implicit_flux_agg

!!!!! de Vrese (01.12.2015) Calculating new state based on flux aggregation see (M. J. Best, A. Beljaars, J. Polcher, and P. Viterbo, 2004) - end

#endif
END MODULE mo_srf_process
