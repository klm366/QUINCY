!>
!! @brief Contains the memory class for the surface.
!!
!! @author
!!  Veronika Gayler, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Veronika Gayler (2015-05-06)
!!
MODULE mo_hd_memory_class
#ifndef __NO_JSBACH__

  USE mo_kind, ONLY: wp

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: t_hd_memory, t_hd_param, t_hd_state

  TYPE t_hd_param
    REAL(wp), POINTER ::     &
      fdir            (:,:), &  !< index of next grid cell in flow direction
      ret_overlflow   (:,:), &  !< retention constant of overland flow
      ret_baseflow    (:,:), &  !< retention constant of baseflow
      ret_riverflow   (:,:), &  !< retention constant of riverflow
      nres_overlflow  (:,:), &  !< number of overland flow reservoirs
      nres_baseflow   (:,:), &  !< number of baseflow reservoirs
      nres_riverflow  (:,:)     !< number of riverflow reservoirs
    INTEGER, POINTER ::      &
      nres_o_max,            &  !< maximum number of overland flow reservoirs
      nres_b_max,            &  !< maximum number of baseflow reservoirs
      nres_r_max                !< maximum number of riverflow reservoirs
  END TYPE t_hd_param

  TYPE t_hd_state
    REAL(wp), POINTER ::     &
      overlflow_res (:,:,:), &  !< content of overland flow reservoir cascade
      baseflow_res  (:,:,:), &  !< content of baseflow reservoir cascade
      riverflow_res (:,:,:)     !< content of groundflow reservoir cascade
  END TYPE t_hd_state

  TYPE t_hd_memory
    TYPE(t_hd_param),   POINTER :: hd_param => NULL()
    TYPE(t_hd_state),   POINTER :: hd_state => NULL()
  END TYPE t_hd_memory


  CHARACTER(len=*), PARAMETER :: modname = 'mo_hd_memory_class'

CONTAINS


#endif
END MODULE mo_hd_memory_class
