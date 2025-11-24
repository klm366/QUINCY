MODULE mo_debug_stw
  USE mo_kind,           ONLY: wp
  USE mo_util_string,    ONLY: int2string, logical2string, real2string
  USE mo_mpi,            ONLY: p_barrier, p_sum, p_global_sum, p_pe
  USE mo_exception,      ONLY: message
  USE mtime,             ONLY: t_datetime => datetime ! Directly from mtime to avoid recursive deps
IMPLICIT NONE
PRIVATE

  interface i2s
    module procedure i2s0
    module procedure i2s1
  end interface i2s

  interface r2s
    module procedure r2s0
    module procedure r2s1
  end interface r2s

  interface l2s
    module procedure l2s0
    module procedure l2s1
  end interface l2s

PUBLIC :: p_barrier, p_sum, p_global_sum, p_pe
PUBLIC :: t_datetime
PUBLIC :: msg, l2s, i2s, r2s, dati2s
PUBLIC :: isnan, anynan, nnan, isinf, anyinf, ninf, isbad, anybad, nbad

CONTAINS

  SUBROUTINE msg(name, text, pe, all_print)

    CHARACTER(len=*), INTENT(IN) :: name
    CHARACTER(len=*), INTENT(IN) :: text
    INTEGER, INTENT(IN), OPTIONAL :: pe
    LOGICAL, INTENT(IN), OPTIONAL :: all_print

    IF (PRESENT(pe)) THEN
      IF (p_pe==pe) CALL message(name, text, all_print=.TRUE.)
    ELSE
      CALL message(name,text,all_print)
    ENDIF

  END SUBROUTINE msg

  FUNCTION dati2s(datetime) RESULT(res)
    CHARACTER(len=20) :: res
    TYPE(t_datetime), TARGET, INTENT(in) :: datetime

    WRITE(res,'(A,I4.4,5(A,I2.2))') ' ',datetime%date%year,'-',datetime%date%month, '-',datetime%date%day,' ', &
                                        datetime%time%hour,':',datetime%time%minute,':',datetime%time%second

  END FUNCTION dati2s

  LOGICAL FUNCTION isnan(r)
    REAL(wp), INTENT(in) :: r

    isnan = r /= r

  END FUNCTION isnan

  LOGICAL FUNCTION anynan(n,r)
    INTEGER,  INTENT(in) :: n
    REAL(wp), INTENT(in) :: r(*)

    INTEGER :: i

    anynan = .TRUE.
    DO i=1,n
      IF (isnan(r(i))) RETURN
    ENDDO
    anynan = .FALSE.

  END FUNCTION anynan

  INTEGER FUNCTION nnan(n,r)
    INTEGER,  INTENT(in) :: n
    REAL(wp), INTENT(in) :: r(*)

    INTEGER :: i, nn

    nn = 0
    DO i=1,n
      IF (isnan(r(i))) nn = nn + 1
    ENDDO
    nnan = nn

  END FUNCTION nnan

  LOGICAL FUNCTION isinf(r)
    REAL(wp), INTENT(in) :: r

    isinf = ABS(r)>HUGE(r)

  END FUNCTION isinf

  LOGICAL FUNCTION anyinf(n,r)
    INTEGER,  INTENT(in) :: n
    REAL(wp), INTENT(in) :: r(*)

    INTEGER :: i

    anyinf = .TRUE.
    DO i=1,n
      IF (isinf(r(i))) RETURN
    ENDDO
    anyinf = .FALSE.

  END FUNCTION anyinf

  INTEGER FUNCTION ninf(n,r)
    INTEGER,  INTENT(in) :: n
    REAL(wp), INTENT(in) :: r(*)

    INTEGER :: i, nn

    nn = 0
    DO i=1,n
      IF (isinf(r(i))) nn = nn + 1
    ENDDO
    ninf = nn

  END FUNCTION ninf

  LOGICAL FUNCTION isbad(r)
    REAL(wp), INTENT(in) :: r

    isbad = isinf(r) .OR. isnan(r)

  END FUNCTION isbad

  LOGICAL FUNCTION anybad(n,r)
    INTEGER,  INTENT(in) :: n
    REAL(wp), INTENT(in) :: r(*)

    anybad = anyinf(n,r) .OR. anynan(n,r)

  END FUNCTION anybad

  INTEGER FUNCTION nbad(n,r)
    INTEGER,  INTENT(in) :: n
    REAL(wp), INTENT(in) :: r(*)

    nbad = ninf(n,r) + nnan(n,r)

  END FUNCTION nbad

  FUNCTION i2s0(n, sep, opt_fmt)
    CHARACTER(len=:), ALLOCATABLE :: i2s0
    INTEGER, INTENT(in) :: n
    CHARACTER(len=1), INTENT(in), OPTIONAL :: sep
    CHARACTER(len=*), INTENT(in), OPTIONAL :: opt_fmt

    IF (PRESENT(sep)) THEN
      IF (sep=='A') THEN
        i2s0 = int2string(n, opt_fmt)
      ELSE
        i2s0 = sep//int2string(n, opt_fmt)
      ENDIF
    ELSE
      i2s0 = ' '//int2string(n, opt_fmt)
    ENDIF

  END FUNCTION i2s0

  FUNCTION i2s1(n, sep, opt_fmt)
    CHARACTER(len=:), ALLOCATABLE :: i2s1
    INTEGER, INTENT(in) :: n(:)
    CHARACTER(len=1), INTENT(in), OPTIONAL :: sep
    CHARACTER(len=*), INTENT(in), OPTIONAL :: opt_fmt

    CHARACTER(len=1) :: se
    CHARACTER(len=:), ALLOCATABLE :: tmp1, tmp2
    INTEGER :: i

    se = ' ';
    tmp1 = ''
    DO i=1,SIZE(n)
      tmp2 = TRIM(tmp1)//se//int2string(n(i), opt_fmt)
      tmp1 = tmp2
      IF (PRESENT(sep)) se = sep ! First separator is always a space
    ENDDO

    i2s1 = TRIM(tmp1)

  END FUNCTION i2s1

  FUNCTION r2s0(x, sep, opt_fmt)
    CHARACTER(len=:), ALLOCATABLE :: r2s0
    REAL(wp), INTENT(in) :: x
    CHARACTER(len=1), INTENT(in), OPTIONAL :: sep
    CHARACTER(len=*), INTENT(in), OPTIONAL :: opt_fmt

    IF (PRESENT(sep)) THEN
      IF (sep=='A') THEN
        r2s0 = real2string(x, opt_fmt)
      ELSE
        r2s0 = sep//real2string(x, opt_fmt)
      ENDIF
    ELSE
      r2s0 = ' '//real2string(x, opt_fmt)
    ENDIF

  END FUNCTION r2s0

  FUNCTION r2s1(x, opt_fmt)
    CHARACTER(len=:), ALLOCATABLE :: r2s1
    REAL(wp), INTENT(in) :: x(:)
    CHARACTER(len=*), INTENT(in), OPTIONAL :: opt_fmt

    CHARACTER(len=:), ALLOCATABLE :: tmp1, tmp2
    INTEGER :: i

    tmp1 = ''
    DO i=1,SIZE(x)
      tmp2 = TRIM(tmp1)//' '//real2string(x(i), opt_fmt)
      tmp1 = tmp2
    ENDDO

    r2s1 = TRIM(tmp1)

  END FUNCTION r2s1

  FUNCTION l2s0(sel, sep)
    CHARACTER(len=:), ALLOCATABLE :: l2s0
    LOGICAL, INTENT(in) :: sel
    CHARACTER(len=1), INTENT(in), OPTIONAL :: sep

    IF (PRESENT(sep)) THEN
      IF (sep=='A') THEN
        l2s0 = logical2string(sel)
      ELSE
        l2s0 = sep//logical2string(sel)
      ENDIF
    ELSE
      l2s0 = ' '//logical2string(sel)
    ENDIF

  END FUNCTION l2s0

  FUNCTION l2s1(sel)
    CHARACTER(len=:), ALLOCATABLE :: l2s1
    LOGICAL, INTENT(in) :: sel(:)

    CHARACTER(len=:), ALLOCATABLE :: tmp1, tmp2
    INTEGER :: i

    tmp1 = ''
    DO i=1,SIZE(sel)
      tmp2 = TRIM(tmp1)//' '//logical2string(sel(i))
      tmp1 = tmp2
    ENDDO

    l2s1 = TRIM(tmp1)

  END FUNCTION l2s1

END MODULE mo_debug_stw
