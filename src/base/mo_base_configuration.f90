!!
!! @brief flags (switches) to en-/disable certain model functionalities
!!
!! @author
!!  Jan Engel, MPI-BGC Jena
!!
!! @par Revision History
!! First version by Jan Engel (2016-03-08)
!!
MODULE mo_base_configuration

  USE mo_kind, ONLY : wp

  IMPLICIT NONE

  TYPE ModelOptionsClass
    ! Output related:
    ! Set output timestep of pools and fluxes
    CHARACTER(8) :: output_interval_pool  = 'daily'   ! options: timestep daily weekly monthly yearly (monthly == 30 days)
    CHARACTER(8) :: output_interval_flux  = 'daily'   ! ( TRIM is used to cut empty character fields)
    ! Set time of output start/end (e.g. only output last 3 years)
    INTEGER :: output_start_date     = 1
    INTEGER :: output_end_date       = 100
    ! Optimality related:
    ! Set whether model runs optimising schemes:
    LOGICAL :: flag_optimal_Nfraction = .FALSE.
    LOGICAL :: flag_optimal_Nleaf     = .FALSE.
  END TYPE ModelOptionsClass

  TYPE ForcingOptionsClass
    ! Driver related:
    ! Set forcing to one of three options
    CHARACTER(8) :: option_forcing    = 'readfile'  ! options: constant readfile generate 
  END TYPE ForcingOptionsClass

  TYPE VegetationOptionsClass
    ! docu
    
  END TYPE VegetationOptionsClass


  TYPE OptionsClass
    ! collects all sub-model specific Config TYPEs and makes them globally accessible
    TYPE( ModelOptionsClass )       :: model
    TYPE( ForcingOptionsClass )     :: forcing
    TYPE( VegetationOptionsClass )  :: vegetation
  END TYPE OptionsClass

  TYPE( OptionsClass ), SAVE :: options


  PRIVATE							! all other definitions of TYPEs are local, but accessible via options%... 

  PUBLIC :: options			        ! only this main options TYPE is global

END MODULE mo_base_configuration
