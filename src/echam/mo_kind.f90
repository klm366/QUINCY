!>
!! @brief set kind and minimum precision of calculation (taken from MPI-ESM 1.2)
!!
!! @author
!!  Sönke Zaehle, MPI-BGC Jena 
!!
!! @par Revision History
!! First version by Sönke Zaehle (2015-12-04)
!!
!! @par ECHAM
!!  this module originates from ECHAM
!!
MODULE mo_kind

  ! L. Kornblueh, MPI, August 2001, added working precision and comments

  IMPLICIT NONE

  ! Number model from which the SELECTED_*_KIND are requested:
  ! 
  !                   4 byte REAL      8 byte REAL
  !          CRAY:        -            precision =   13
  !                                    exponent  = 2465
  !          IEEE:    precision =  6   precision =   15
  !                   exponent  = 37   exponent  =  307
  ! 
  ! Most likely this are the only possible models.

  ! Floating point section:

  ! Floating point section:

  INTEGER, PARAMETER :: ps = 6
  INTEGER, PARAMETER :: rs = 37

  INTEGER, PARAMETER :: pd = 12
  INTEGER, PARAMETER :: rd = 307

  INTEGER, PARAMETER :: pi4 = 9
  INTEGER, PARAMETER :: pi8 = 14

  INTEGER, PARAMETER :: sp = SELECTED_REAL_KIND(ps,rs)
  INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(pd,rd)

  ! Floating point working precision

  INTEGER, PARAMETER :: wp = dp

  REAL(wp), PARAMETER :: eps12 = 0.000000000001_wp
  REAL(wp), PARAMETER :: eps8  = 0.00000001_wp
  REAL(wp), PARAMETER :: eps4  = 0.0001_wp

  ! Integer section

  INTEGER, PARAMETER :: i4 = SELECTED_INT_KIND(pi4)
  INTEGER, PARAMETER :: i8 = SELECTED_INT_KIND(pi8)
!#if defined (__64BIT__) || defined (__LP64__) || defined (__x86_64__) || defined (__x86_64)
  INTEGER, PARAMETER :: widx = i8
!#else
!  INTEGER, PARAMETER :: widx = i4
!#endif

CONTAINS

  SUBROUTINE print_kinds

    PRINT *, 'single precision : ', sp
    PRINT *, 'double precision : ', dp
    PRINT *, 'working precision: ', wp
    PRINT *, '4 byte integer   : ', i4
    PRINT *, '8 byte integer   : ', i8
    PRINT *, 'index integer    : ', widx

  END SUBROUTINE print_kinds

END MODULE mo_kind

