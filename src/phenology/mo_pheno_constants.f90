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

  !< some useful numbers for easy identification of photosynthetic pathways and growth forms
  !! need to be parameters (i.e. constant across model runtime)
  !! 
  INTEGER, PARAMETER ::           &
      ievergreen    = 1         , &
      isummergreen  = 2         , &
      iraingreen    = 3         , &
      iperennial    = 4

  CHARACTER(len=*), PARAMETER, PRIVATE :: modname = 'mo_pheno_constants'

CONTAINS

  !-----------------------------------------------------------------------------------------------------
  !> initialize constants used in the phenology process
  !! 
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE init_pheno_constants

    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    ! 
    ! ---------------------------
    ! 0.2 Local
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':init_pheno_constants'

    
    ! 
    ! currently doing nothing, but called in model_init anyway
    ! PFT specific parameters you may find in the lctlib file
    ! 


  END SUBROUTINE init_pheno_constants

END MODULE mo_pheno_constants
