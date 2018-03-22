!>
!! @brief Contains constants for the harvest processes
!!
!! @author
!!  Julia Nabel, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Julia Nabel (2018-03-22)
!!
MODULE mo_harvest_constants
#ifndef __NO_JSBACH__

  IMPLICIT NONE
  PUBLIC

  INTEGER , PARAMETER, PUBLIC  :: HARVEST_OLDEST_TREES =  1 ! Harvest until empty or demand filled: oldest trees first

#endif
END MODULE mo_harvest_constants
