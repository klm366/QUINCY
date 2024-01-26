!>
!! @brief Contains the utilities for the topmodel wetland diagnostics
!!
!! @author
!!  Stiig Wilkenskjeld, MPI-M Hamburg
!!
!! @par Revision History
!! Adopted from JSBACH3 code by Thomas Kleinen, 2022-12-06
!!
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/license .
!!

MODULE mo_top_util
#ifndef __NO_JSBACH__

  USE mo_kind,          ONLY: wp
!  USE mo_exception,     ONLY: finish
!  USE mo_util,          ONLY: real2string

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: gamm_p, gamm_inv_p, gamm_q, gamm_ln
  PUBLIC :: time_relaxation

  CHARACTER(len=*), PARAMETER :: modname = 'mo_top_util'

CONTAINS

  !!
  ! @brief Update a moving average (e-fold relaxation) with a new value
  !        Either dt and t_relax or nstep must be present
  !
  !!
  PURE ELEMENTAL SUBROUTINE time_relaxation(inst,relax,dt,t_relax,nstep,linit)
    REAL(wp),           INTENT(in)    :: inst    ! Instantaneous value to include in relaxation [?]
    REAL(wp),           INTENT(inout) :: relax   ! Mean value to update [?]
    REAL(wp), OPTIONAL, INTENT(in)    :: dt      ! Time step between updates [s]
    REAL(wp), OPTIONAL, INTENT(in)    :: t_relax ! E-fold time of the relaxation [s]
    REAL(wp), OPTIONAL, INTENT(in)    :: nstep   ! Number of update steps in the e-fold time
    LOGICAL,  OPTIONAL, INTENT(in)    :: linit   ! Initialize averages if linit=.TRUE.

    REAL(wp) :: n

    IF (PRESENT(nstep)) THEN
      n = nstep
    ELSE IF (PRESENT(dt) .AND. PRESENT(t_relax)) THEN
      n = t_relax / dt
    ELSE
      n = 1._wp ! relax = inst
    ENDIF
    IF (PRESENT(linit)) THEN
      IF (linit) n = 1._wp ! relax = inst => initialize averaging
    ENDIF

    relax = (relax * (n - 1._wp) + inst) / n

  END SUBROUTINE time_relaxation

  !------------------------------------------------------------------------------------
  ! Helper routines -- solve gamma function. See Numerical Recipes, Chap. 6
  !------------------------------------------------------------------------------------

  !   Returns the incomplete gamma function P (a, x).
  PURE ELEMENTAL FUNCTION gamm_p(a,x)
    REAL(wp), INTENT(in) :: a, x
    REAL(wp)             :: gamm_p
    !   USES gcf, gser
    REAL(wp)             :: gammcf, gamser, gln

    CHARACTER(len=*), PARAMETER :: routine = modname//':gamm_p'

!    IF (x < 0._wp .OR. a < 0._wp) CALL finish(routine,'bad arguments: '//real2string(x)//' '//real2string(a))
    IF (x < a+1._wp) THEN
    !    Use the series representation.
      CALL gser(gamser,a,x,gln)
      gamm_p = gamser
    ELSE
    !      Use the continued fraction representation
      CALL gcf(gammcf,a,x,gln)
      gamm_p = 1.-gammcf
    !      and take its complement.
    END IF

  END FUNCTION gamm_p

  !------------------------------------------------------------------------------------
  PURE ELEMENTAL FUNCTION gamm_inv_p(p, a)

    REAL(wp), INTENT(IN) :: p, a
    REAL(wp)             :: gamm_inv_p

    INTEGER              :: j
    REAL(wp)             :: x, err, t, u, pp, lna1, afac, a1, gln
    REAL(wp), PARAMETER  :: EPS = 1.e-8_wp

    CHARACTER(len=*), PARAMETER :: routine = modname//':gamm_inv_p'

    a1 = a-1
    gln = gamm_ln(a)

!    IF (a <= 0._wp) CALL finish(routine,'a must be pos in gamm_inv_p')
    IF (p >= 1._wp) THEN
      gamm_inv_p = MAX(100._wp,a + 100._wp*sqrt(a))
      RETURN
    END IF
    IF (p <= 0._wp) THEN
      gamm_inv_p = 0.0_wp
      RETURN
    END IF
    IF (a > 1._wp) THEN
      lna1 = LOG(a1)
      afac = EXP(a1*(lna1-1._wp)-gln)
      IF (p < 0.5_wp) THEN
        pp=p
      ELSE
        pp = 1._wp-p
      END IF
      t = SQRT(-2.*LOG(pp))
      x = (2.30753_wp+t*0.27061_wp)/(1._wp+t*(0.99229_wp+t*0.04481_wp)) - t
      IF (p < 0.5) x = -x
      x = MAX(1.e-3_wp, a*(1._wp-1._wp/(9._wp*a)-x/(3._wp*SQRT(a)))**3._wp)
    ELSE
      t = 1.0_wp - a*(0.253_wp+a*0.12_wp);
      IF (p < t) THEN
        x = (p/t)**(1._wp/a)
      ELSE
        x = 1._wp-LOG(1._wp-(p-t)/(1._wp-t))
      END IF
    END IF
    DO j=1,12
      IF (x <= 0.0_wp) THEN
        gamm_inv_p = 0.0_wp
        RETURN
      END IF
      err = gamm_p(a,x) - p
      IF (a > 1._wp) THEN
        t = afac*EXP(-(x-a1)+a1*(LOG(x)-lna1))
      ELSE
        t = EXP(-x+a1*LOG(x)-gln)
      END IF
      u = err/t
      t = u/(1._wp-0.5_wp*MIN(1._wp,u*((a-1._wp)/x - 1._wp)))
      x = x - t
      IF (x <= 0._wp) THEN
        x = 0.5_wp*(x + t)
      END IF
      IF (ABS(t) < EPS*x ) EXIT
    ENDDO
    gamm_inv_p = x

  END FUNCTION gamm_inv_p

  !------------------------------------------------------------------------------------
  ! Returns the incomplete gamma function Q(a, x) ≡ 1 − P (a, x).
  PURE ELEMENTAL FUNCTION gamm_q(a,x)
    REAL(wp), INTENT(in) :: a,x
    REAL(wp) :: gamm_q
    ! USES gcf,gser
    REAL(wp) :: gammcf,gamser,gln

    CHARACTER(len=*), PARAMETER :: routine = modname//':gamm_q'

!    IF (x < 0._wp .OR. a < 0._wp) CALL finish(routine,'bad arguments')
    IF (x < a+1._wp) THEN
    !    Use the series representation
      CALL gser(gamser,a,x,gln)
      gamm_q = 1.-gamser
    !    and take its complement.
    ELSE
    !    Use the continued fraction representation.
      CALL gcf(gammcf,a,x,gln)
      gamm_q = gammcf
    END IF

  END FUNCTION gamm_q

  !------------------------------------------------------------------------------------
  !   Returns the incomplete gamma function P (a, x) evaluated by its series representation as
  !   gamser . Also RETURNs ln Γ(a) as gln.
  PURE ELEMENTAL SUBROUTINE gser(gamser,a,x,gln)
    REAL(wp), INTENT(inout) :: gamser, gln
    REAL(wp), INTENT(in)    :: a, x

    INTEGER,  PARAMETER     :: ITMAX = 100
    REAL(wp), PARAMETER     :: EPS = 3.e-7_wp
    !   USES topmod_gammln
    INTEGER  :: n
    REAL(wp) :: ap,del,total

    CHARACTER(len=*), PARAMETER :: routine = modname//':gser'

    gln = gamm_ln(a)
    IF (x <= 0._wp) THEN
!      IF (x < 0._wp) CALL finish(routine,'x < 0 in topmod_gser')
      gamser = 0._wp
      RETURN
    END IF
    ap = a
    total = 1._wp/a
    del = total
    DO n = 1,ITMAX
      ap = ap + 1._wp
      del = del * x / ap
      total = total + del
      IF (ABS(del) < ABS(total)*EPS) EXIT
    END DO
!    IF (n == ITMAX) CALL finish(routine,'a too large, ITMAX too small in gser')
    gamser = total*EXP(-x+a*LOG(x)-gln)

  END SUBROUTINE gser

  !------------------------------------------------------------------------------------
  !   Returns the incomplete gamma function Q(a, x) evaluated by its continued fraction repre-
  !   sentation as gammcf. Also RETURNs ln Γ(a) as gln.
  !   Parameters: ITMAX is the maximum allowed number of iterations; EPS is the relative accu-
  !   racy; FPMIN is a number near the smallest representable floating-point number.
  PURE ELEMENTAL SUBROUTINE gcf(gammcf,a,x,gln)
    REAL(wp), INTENT(inout) :: gammcf, gln
    REAL(wp), INTENT(in)    :: a,x

    REAL(wp), PARAMETER     :: EPS = 3.e-7
    REAL(wp), PARAMETER     :: FPMIN = 1.e-30
    INTEGER,  PARAMETER     :: ITMAX = 100
    !   USES topmod_gammln
    INTEGER  :: i
    REAL(wp) :: an, b, c, d, del, h

    CHARACTER(len=*), PARAMETER :: routine = modname//':gcf'

    gln = gamm_ln(a)
    b = x+1.-a
    !    Set up for evaluating continued fraction by modified
    c = 1./FPMIN
    !    Lentz’s method (§5.2) with b0  =  0.
    d = 1./b
    h = d
    DO i = 1,ITMAX
    !    Iterate to convergence.
      an = -i*(i-a)
      b = b+2.
      d = an*d+b
      IF (ABS(d) .LT. FPMIN) d = FPMIN
      c = b+an/c
      IF (ABS(c) .LT. FPMIN) c = FPMIN
      d = 1./d
      del = d*c

      h = h*del
      IF (ABS(del-1.) < EPS) exit
    END DO
!    IF (i  == ITMAX) CALL finish(routine,'a too large, ITMAX too small in gcf')

    gammcf = EXP(-x+a*LOG(x)-gln)*h
    !    Put factors in front.

  END SUBROUTINE gcf

  !------------------------------------------------------------------------------------
  PURE ELEMENTAL FUNCTION gamm_ln(xx)
    REAL(wp), INTENT(in) :: xx
    REAL(wp)             :: gamm_ln
    ! Returns the value ln[Γ(xx )] for xx > 0.

    INTEGER        :: j
    REAL(wp)       :: ser,tmp,x,y
    REAL(wp), PARAMETER :: cof(6) = (/76.18009172947146_wp,-86.50532032941677_wp, &
    24.01409824083091_wp,-1.231739572450155_wp,0.1208650973866179e-2_wp, &
    -0.5395239384953e-5_wp/)
    REAL(wp), PARAMETER :: stp = 2.5066282746310005_wp

    x = xx
    y = x
    tmp = x + 5.5_wp
    tmp = (x + 0.5_wp) * LOG(tmp) - tmp
    ser = 1.000000000190015_wp
    DO j = 1,6
      y = y + 1._wp
      ser = ser + cof(j) / y
    END DO
    gamm_ln = tmp + LOG(stp * ser / x)

  END FUNCTION gamm_ln

#endif
END MODULE mo_top_util
