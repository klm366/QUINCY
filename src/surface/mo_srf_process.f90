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

  PUBLIC :: surface_temp_and_qsat_implicit

  CHARACTER(len=*), PARAMETER :: modname = 'mo_srf_process'

CONTAINS

  SUBROUTINE surface_temp_and_qsat_implicit(klon, steplen,              &
    &            pcp,                                         &
    &            pescoe, pfscoe, peqcoe, pfqcoe,                       &
    &            psold, pqsold, pdqsold,                               &
    &            pnetrad, pgrdfl,                                      &
    &            pcfh, pcair, pcsat, pfracsu, pgrdcap,                 &
    &            psnew, pqsnew)

    ! TBD: Use cvdifts from mo_echam_vdiff_params once the unified vdiff for ICON/ECHAM is ready
#ifdef __ICON__
    USE mo_echam_vdiff_params, ONLY: cvdifts   ! TBD: change to 1. for standalone JSBACH
#else
    USE mo_physc2, ONLY: cvdifts
#endif
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
    REAL(wp), INTENT(out)   :: psnew(klon), pqsnew(klon)
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
    zcolin(:) = pgrdcap(:)*zicp(:) +                                                           &
              &      pdt * (zicp(:) * 4._wp * zemiss_def * stbo * ((zicp(:) * psold(:))**3)  - &
              &      pcfh(:) * (zca(:) * peqcoe(:) - zcs(:) -                                  &
              &      zicp(:) * pc16 * psold(:) *                                               &
              &      (pcair(:) * peqcoe(:) - pcsat(:)))*                                       &
              &      zicp(:) * pdqsold(:))
    !
    zcohfl(:) = -pdt * pcfh(:) * (pescoe(:) - 1._wp)
    !
    zcoind(:) = pdt * (pnetrad(:) + pcfh(:) * pfscoe(:) +  pcfh(:)*                            &
              &      ((zca(:) * peqcoe(:) - zcs(:)) * pqsold(:) + zca(:) * pfqcoe(:)  -        &
              &      zicp(:) * pc16 * psold(:) *                                               &
              &      ((pcair(:) * peqcoe(:) - pcsat(:)) * pqsold(:) +                          &
              &      pcair(:) * pfqcoe(:))) + pgrdfl(:))
    !
    psnew(:)  = (zcolin(:) * psold(:) + zcoind(:)) / (zcolin(:) + zcohfl(:))
    pqsnew(:) = pqsold(:) + zicp(:) * pdqsold(:) * (psnew(:) - psold(:))

  END SUBROUTINE surface_temp_and_qsat_implicit

#endif
END MODULE mo_srf_process
