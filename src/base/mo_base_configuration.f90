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

  TYPE VegetationConfigClass
    ! currently: values are set in plant_model_driver section 0.1 -> later: read from file
    ! Output related:
    ! Set output timestep of pools and fluxes
    CHARACTER(8) :: output_interval_pool  = 'timestep'   ! options: timestep daily weekly monthly yearly (monthly == 30 days)
    CHARACTER(8) :: output_interval_flux  = 'timestep'   ! ( TRIM is used to cut empty character fields)
    ! Set time of output start/end (e.g. only output last 3 years)
    INTEGER :: output_start_date     = 1
    INTEGER :: output_end_date       = 100
    ! Driver related:
    ! Set whether driving variables (temperature, light etc.) vary seasonally: 0 == NO, 1 == YES
    CHARACTER(8) :: option_forcing    = 'constant'  ! options: constant readfile generate 
    ! Optimality related:
    ! Set whether model runs optimising schemes:
    LOGICAL :: flag_optimal_Nfraction = .FALSE.
    LOGICAL :: flag_optimal_Nleaf     = .FALSE.
  END TYPE VegetationConfigClass


  TYPE ModelConfigClass
  ! collects all sub-model specific Config TYPEs and makes them globally accessible
    TYPE( VegetationConfigClass ) :: vegetation			! values set in plant_model_driver.f90 (later: read from file)

  END TYPE ModelConfigClass

  TYPE( ModelConfigClass ) :: model_config


  PRIVATE							! all other definitions of TYPEs are local, but accessible via lodel_config%... 

  PUBLIC :: model_config			! only this main config TYPE is global

END MODULE mo_base_configuration
