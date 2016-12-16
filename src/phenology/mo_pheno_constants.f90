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
  REAL(wp), SAVE ::               &
      min_leaf_age                    !< minimum leaf age before senescence is permitted (days) 
      

  !< some useful numbers for easy identification of photosynthetic pathways and growth forms
  !! need to be parameters (i.e. constant across model runtime)
  INTEGER, PARAMETER ::           &
      ievergreen    = 1         , &
      iraingreen    = 2         , &
      isummergreen  = 3         , &
      iperennial    = 4


CONTAINS

  !-----------------------------------------------------------------------------------------------------
  !> initialize constants used in the phenology process
  !! 
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE init_pheno_constants

    IMPLICIT NONE


    min_leaf_age                = 10.0_wp 


  END SUBROUTINE init_pheno_constants

END MODULE mo_pheno_constants
