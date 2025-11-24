!>
!! @brief Contains the utilities for the topmodel wetland diagnostics
!!
!! @author
!!  Stiig Wilkenskjeld, MPI-M Hamburg
!!
!! @par Revision History
!! Adopted from JSBACH3 code by Thomas Kleinen, 2022-12-06
!! Added interpolation functionality,           2024-06-06
!!
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/license .
!!

MODULE mo_top_util
#ifndef __NO_JSBACH__

  USE mo_kind,          ONLY: wp
  USE mo_exception,     ONLY: finish, message

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: gamm_p, gamm_inv_p, gamm_q, gamm_ln
  PUBLIC :: time_relaxation
  PUBLIC :: interpol_reset, interpol_done, interpol_new, interpol_update, wp, t_interp_var, t_interp_var_p

  CHARACTER(len=*), PARAMETER :: modname = 'mo_top_util'

  TYPE t_interp_var
    REAL(wp), POINTER :: curr(:)                 ! Current value of the interpolation
    REAL(wp), POINTER :: a(:), b(:), c(:)        ! Coefficents for quadratic polynomias
    REAL(wp), POINTER :: excess(:)               ! Current excess area
    REAL(wp), POINTER :: excess_max(:)           ! Maximum of excess area
    REAL(wp), POINTER :: add(:)                  ! Area to add to current interpolation
    INTEGER,  POINTER :: flg(:)                  ! Bitwise flags for points to update "add" as well as modified endpoints:
                                                 !   1 : add should be updated
                                                 !   2 : first endpoint has been modified
                                                 !   4 : second endpoint has been modified
    REAL(wp)          :: valid_range(2)          ! Lower and upper boundaries for output
    REAL(wp)          :: frac_use_full           ! Fraction of valid_range for which the base value should be used
    REAL(wp)          :: dt                      ! Time between interpolation time steps
    INTEGER           :: nintrp                  ! Number of total interpolation time steps
    INTEGER           :: xcurr                   ! Number of current interpolation time step
    LOGICAL           :: ext_var                 ! Is var field allocated externally?
  END TYPE t_interp_var

  TYPE t_interp_var_p
    TYPE(t_interp_var), POINTER :: p => NULL()
  END TYPE t_interp_var_p

CONTAINS

  !!
  ! @brief Constructs or resets an interpolation structure
  !
  !!
  FUNCTION interpol_reset(n,intrp,var)
  !$acc routine seq
    TYPE(t_interp_var),           POINTER :: interpol_reset
    INTEGER,                   INTENT(in) :: n
    TYPE(t_interp_var), OPTIONAL, POINTER :: intrp
    REAL(wp),           OPTIONAL, POINTER :: var(:)

    TYPE(t_interp_var), POINTER :: res
    REAL(wp) :: tmp
    INTEGER :: i

    NULLIFY(res)
    IF (PRESENT(intrp)) THEN
      res => intrp
    ENDIF
    IF (.NOT. ASSOCIATED(res)) THEN
      ALLOCATE(res)
      !$acc enter data pcreate(res)
      !$acc enter data attach(res)
      NULLIFY(res%curr,res%a,res%b,res%c,res%excess,res%excess_max,res%add,res%flg)
    ENDIF
    res%valid_range   = (/-HUGE(tmp),HUGE(tmp)/)
    res%frac_use_full = -1._wp
    res%dt            = 0._wp
    res%nintrp        = 0
    res%xcurr         = 0
    res%ext_var = .FALSE.
    IF (.NOT. ASSOCIATED(res%curr)) THEN
      IF (PRESENT(var)) THEN
        res%curr => var
        res%ext_var = .TRUE.
      ELSE
        ALLOCATE(res%curr(n))
        !$acc enter data pcreate(res%curr)
        !$acc enter data attach(res%curr)
      ENDIF
    ENDIF
    IF (.NOT. ASSOCIATED(res%a))          ALLOCATE(res%a(n))
    IF (.NOT. ASSOCIATED(res%b))          ALLOCATE(res%b(n))
    IF (.NOT. ASSOCIATED(res%c))          ALLOCATE(res%c(n))
    IF (.NOT. ASSOCIATED(res%excess))     ALLOCATE(res%excess(n))
    IF (.NOT. ASSOCIATED(res%excess_max)) ALLOCATE(res%excess_max(n))
    IF (.NOT. ASSOCIATED(res%add))        ALLOCATE(res%add(n))
    !$acc enter data pcreate(res%a,res%b,res%c,res%excess,res%excess_max,res%add)
    !$acc enter data attach(res%a,res%b,res%c,res%excess,res%excess_max,res%add)
    IF (.NOT. ASSOCIATED(res%flg)) THEN
      ALLOCATE(res%flg(n))
      DO i=1,n
        res%flg(i) = 0
      ENDDO
    ELSE
      DO i=1,n
        res%flg(i) = ISHFT(IAND(res%flg(i),4),-1) ! Preserve status of previous endpoint as status of current startpoint
      ENDDO
    ENDIF
    DO i=1,n
      res%excess    (i) = 0._wp
      res%excess_max(i) = 0._wp
      res%add       (i) = 0._wp
    ENDDO

    interpol_reset => res

  END FUNCTION interpol_reset

  !!
  ! @brief Deep destruction of an interpolation structure
  !
  !!
  SUBROUTINE interpol_done(intrp)
  !$acc routine seq
    TYPE(t_interp_var), POINTER :: intrp

    TYPE(t_interp_var), POINTER :: done

    done => intrp
    IF ((.NOT. done%ext_var) .AND. ASSOCIATED(done%curr)) DEALLOCATE(done%curr)
    IF (ASSOCIATED(done%a))          DEALLOCATE(done%a)
    IF (ASSOCIATED(done%b))          DEALLOCATE(done%b)
    IF (ASSOCIATED(done%c))          DEALLOCATE(done%c)
    IF (ASSOCIATED(done%excess))     DEALLOCATE(done%excess)
    IF (ASSOCIATED(done%excess_max)) DEALLOCATE(done%excess_max)
    IF (ASSOCIATED(done%add))        DEALLOCATE(done%add)
    IF (ASSOCIATED(done%flg))        DEALLOCATE(done%flg)


    DEALLOCATE(done)

  END SUBROUTINE interpol_done

  !!
  ! @brief Initializes a new interpolation
  !
  !!
  FUNCTION interpol_new(dt,nintrp,y1,y2,y3,valid_range,frac_use_full,intrp,var)
  !$acc routine seq
    TYPE(t_interp_var),              POINTER :: interpol_new
    REAL(wp),                     INTENT(in) :: dt
    INTEGER,                      INTENT(in) :: nintrp
    REAL(wp),                     INTENT(in) :: y1(:), y2(:), y3(:)
    REAL(wp),           OPTIONAL, INTENT(in) :: valid_range(2)
    REAL(wp),           OPTIONAL, INTENT(in) :: frac_use_full
    TYPE(t_interp_var), OPTIONAL,    POINTER :: intrp
    REAL(wp),           OPTIONAL,    POINTER :: var(:)

    TYPE(t_interp_var), POINTER :: res
    REAL(wp) :: ys, ye, dx, dy, span
    INTEGER :: i, n, flgs

    ! flgs: Bitwise
    !   1 : y1 < lower threshold
    !   2 : y1 > upper threshold
    !   4 : y2 < lower threshold
    !   8 : y2 > upper threshold
    !  16 : y3 < lower threshold
    !  32 : y3 > upper threshold

    n = SIZE(y1)
    res => interpol_reset(n,intrp,var)
    IF (PRESENT(valid_range))   res%valid_range(:) = valid_range(:)
    IF (PRESENT(frac_use_full)) res%frac_use_full  = frac_use_full
    res%dt     = dt
    res%nintrp = nintrp

    span = res%valid_range(2) - res%valid_range(1)

    DO i=1,n
      ys = (y1(i)+y2(i))/2.0_wp
      ye = (y2(i)+y3(i))/2.0_wp

      flgs = 0
      IF (frac_use_full >= 0._wp) THEN
        IF (((y1(i)-res%valid_range(1))/span) <= frac_use_full) flgs = IOR(flgs, 1)
        IF (((y2(i)-res%valid_range(1))/span) <= frac_use_full) flgs = IOR(flgs, 4)
        IF (((y3(i)-res%valid_range(1))/span) <= frac_use_full) flgs = IOR(flgs,16)
        IF (((res%valid_range(2)-y1(i))/span) <= frac_use_full) flgs = IOR(flgs, 2)
        IF (((res%valid_range(2)-y2(i))/span) <= frac_use_full) flgs = IOR(flgs, 8)
        IF (((res%valid_range(2)-y3(i))/span) <= frac_use_full) flgs = IOR(flgs,32)

        res%flg(i) = IOR(res%flg(i),ISHFT(flgs,4))

        IF ((IAND(res%flg(i),2)/=0) .OR. (MINVAL(ABS((/4,8,10/)-IAND(flgs,15)))==0)) THEN
          ys = y2(i)
          res%flg(i) = IOR(res%flg(i),2)
        ENDIF
        IF (MINVAL(ABS((/4,8,16,32,40/)-IAND(flgs,60)))==0) THEN
          ye = y3(i)
          res%flg(i) = IOR(res%flg(i),4)
        ENDIF
      ENDIF

      dy = ye - ys
      dx = REAL(nintrp,wp)

      res%b(i) = 6._wp * (y2(i) - ys - dy/3._wp) / dx
      res%a(i) = (dy - res%b(i) * dx) / (dx * dx)
      res%c(i) = ys

    ENDDO

    interpol_new => res

  END FUNCTION interpol_new

  !!
  ! @brief Go to next step of an interpolation
  !
  !!
  SUBROUTINE interpol_update(intrp)
  !$acc routine seq
    TYPE(t_interp_var), POINTER :: intrp

    REAL(wp) :: x, y, tmp, su
    INTEGER :: i, n, nleft

    IF (intrp%xcurr >= intrp%nintrp) RETURN

    intrp%xcurr = intrp%xcurr + 1
    x = REAL(intrp%xcurr,wp)*intrp%dt
    nleft = intrp%nintrp - intrp%xcurr
    n = SIZE(intrp%a)
    su = 0._wp
    DO i=1,n
      y = intrp%a(i)*x*x + intrp%b(i)*x + intrp%c(i) ! + intrp%add(i)
      su = su + y
      intrp%excess(i) = intrp%excess(i) - intrp%add(i) * intrp%dt
      tmp = y - intrp%valid_range(1)
      IF (tmp < 0._wp) THEN
        intrp%excess(i) = intrp%excess(i) + tmp * intrp%dt
        y = intrp%valid_range(1)
      ENDIF
      tmp = y - intrp%valid_range(2)
      IF (tmp > 0._wp) THEN
        intrp%excess(i) = intrp%excess(i) + tmp * intrp%dt
        y = intrp%valid_range(2)
      ENDIF
      IF ((ABS(intrp%excess(i)) == intrp%excess_max(i)) .AND. (intrp%excess_max(i) > 0)) intrp%flg(i) = IOR(intrp%flg(i),1)
      IF (intrp%excess_max(i) < ABS(intrp%excess(i))) &
          intrp%excess_max(i) = ABS(intrp%excess(i))
      IF (IAND(intrp%flg(i),1)/=0) intrp%add(i) = intrp%excess(i)/intrp%dt/REAL(nleft,wp)*2._wp
      intrp%curr(i) = y
    ENDDO

  END SUBROUTINE interpol_update


  !!
  ! @brief Update a moving average (e-fold relaxation) with a new value
  !        Either dt and t_relax or nstep must be present
  !
  !!
#ifndef _OPENACC
  PURE &
#endif
  SUBROUTINE time_relaxation(inst,relax,dt,t_relax,nstep,linit)
  !$acc routine seq
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
#ifndef _OPENACC
  PURE &
#endif
  FUNCTION gamm_p(a,x)
  !$acc routine seq
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
#ifndef _OPENACC
  PURE &
#endif
  FUNCTION gamm_inv_p(p, a)
  !$acc routine seq

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
  ! Returns the incomplete gamma function Q(a, x) = 1 - P (a, x).
#ifndef _OPENACC
  PURE &
#endif
  FUNCTION gamm_q(a,x)
  !$acc routine seq
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
  !   gamser . Also RETURNs ln gamma(a) as gln.
#ifndef _OPENACC
  PURE &
#endif
  SUBROUTINE gser(gamser,a,x,gln)
  !$acc routine seq
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
  !   sentation as gammcf. Also RETURNs ln gamma(a) as gln.
  !   Parameters: ITMAX is the maximum allowed number of iterations; EPS is the relative accu-
  !   racy; FPMIN is a number near the smallest representable floating-point number.
#ifndef _OPENACC
  PURE &
#endif
  !$acc routine seq
  SUBROUTINE gcf(gammcf,a,x,gln)
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
    !    Lentz's method (paragraph 5.2) with b0 = 0.
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
#ifndef _OPENACC
  PURE &
#endif
  FUNCTION gamm_ln(xx)
  !$acc routine seq
    REAL(wp), INTENT(in) :: xx
    REAL(wp)             :: gamm_ln
    ! Returns the value ln[Gamma(xx )] for xx > 0.

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
