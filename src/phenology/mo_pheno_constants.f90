!>
!! @brief Contains constants for the phenology process
!!
!! @author
!!  Jan Engel, MPI-BGC Jena
!!
!! @par Revision History
!! First version by Jan Engel (20161214)
!!
MODULE mo_pheno_constants

  USE mo_kind,                 ONLY: wp


  IMPLICIT NONE

  PUBLIC


  !< Phenology constants (parameter)
  !! 
  REAL(wp), SAVE :: &
      my_pheno_const  ! replace!


CONTAINS

  !-----------------------------------------------------------------------------------------------------
  !> initialize constants used in the phenology process
  !! 
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE init_pheno_constants


    IMPLICIT NONE

    my_pheno_const = 0.1_wp   ! replace!


  END SUBROUTINE init_pheno_constants

END MODULE mo_pheno_constants
