!!
!! @brief Contains the routines for tile parameters 
!!
!! @author
!!  Jan Engel, MPI-BGC Jena 
!!
!! @par Revision History
!! First version by Jan Engel (2016-03-10)
!!
MODULE mo_land_parameters

  USE mo_kind, ONLY: wp

  TYPE AccumulationIntervalTileClass 
    ! length of period over which the values are accumulated/ averaged
    ! see mo_base_memory for defintion of TYPEs that keep the starting/default value
    REAL(wp) :: t_air = 1.
    REAL(wp) :: q_air = 1.
    REAL(wp) :: press_surf = 1.
  END TYPE AccumulationIntervalTileClass

END MODULE mo_land_parameters
