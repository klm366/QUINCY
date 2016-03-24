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
    REAL(wp) :: t_air = 1._wp
    REAL(wp) :: q_air = 1._wp
    REAL(wp) :: press_surf = 1._wp

    REAL(wp) :: weekly     = 7._wp     ! 7 days; constant across the model
    REAL(wp) :: monthly    = 30._wp    ! 30 days; constant across the model
    REAL(wp) :: longterm   = 10950._wp ! 10950 days (30 years); constant across the model
  END TYPE AccumulationIntervalTileClass

END MODULE mo_land_parameters
