!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/license .
!!
!! @brief Def. max. filename length and standard I/O-units
!!
!! @par Revision History
!!  copied to QUINCY by Jan Engel MPI-BGC (2016-03-18)
!!
!! @par ECHAM
!!  this module originates from ECHAM, file is left untouched besides the preprocessor directives
!!
MODULE mo_io_units

  IMPLICIT NONE
  
  PUBLIC
  
  ! This paramter is taken from /usr/include/stdio.h (ANSI C standard). If 
  ! problems with filename length appear, check the before mentioned file.
  
  INTEGER, PARAMETER :: filename_max = 1024
  
  ! Standard I/O-units

! preprocessor directived commented to make gfortran work, Jan Engel, 20160404  
!#ifdef hpux
!  INTEGER, PARAMETER :: nerr  = 7     ! error output
!#else
  INTEGER, PARAMETER :: nerr  = 0     ! error output
!#endif
  INTEGER, PARAMETER :: nlog  = 1     ! standard log file unit
  INTEGER, PARAMETER :: nin   = 5     ! standard input
  INTEGER, PARAMETER :: nout  = 6     ! standard output  
  
  INTEGER, PARAMETER, PRIVATE :: none = -1  ! unit given back, when nothing 
                                            ! in the allowed range is available
  
END MODULE mo_io_units
