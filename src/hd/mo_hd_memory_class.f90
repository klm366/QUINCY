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
    REAL(wp), POINTER ::           &
      hd_lsm  (:,:),    & !< HD model land sea mask (1: land, 0: ocean)
      fdir    (:,:),    & !< river flow direction (-1: ocean, 0: ocean inflow, 5: cells without outflow)
      alf_k   (:,:),    & !< retention parameter k for overland flow
      alf_n   (:,:),    & !< number of reservoirs n for overland flow
      arf_k   (:,:),    & !< retention parameter k for river flow
      arf_n   (:,:),    & !< number of reservoirs n for river flow
      agf_k   (:,:),    & !< retention parameter k for base flow
      hd_area (:,:),    & !< grid cell area [m2]
      filnew  (:,:),    & !< longitude index of flow destination cell
      fibnew  (:,:)       !< latitude index of flow destination cell
  END TYPE t_hd_param

  TYPE t_hd_state
    REAL(wp), POINTER ::           &
      flfmem   (:,:),    & !< overland flow reservoir (fed from surface runoff)
      frfmem   (:,:),    & !< river flow reservoir (inflow from surrounding grid boxes)
      fgmem    (:,:),    & !< base flow reservoir (fed from soil drainage)
      finfl    (:,:)       !< 
  END TYPE t_hd_state

  TYPE t_hd_memory
    TYPE(t_hd_param),   POINTER :: hd_param => NULL()
    TYPE(t_hd_state),   POINTER :: hd_state => NULL()
  END TYPE t_hd_memory

!vg   ! Note: Would be nicer to be able to overload the "=" operator, but that would be more
!vg   ! complicated because we only want to copy parts of the arrays (ics:ice,iblk)
!vg   INTERFACE copy_hd_mem
!vg     MODULE PROCEDURE copy_t_hd_state
!vg   END INTERFACE copy_srf_mem
!vg 
!vg   INTERFACE add_hd_mem
!vg     MODULE PROCEDURE add_t_hd_state
!vg   END INTERFACE add_hd_mem

  CHARACTER(len=*), PARAMETER :: modname = 'mo_hd_memory_class'

CONTAINS

!vg   SUBROUTINE copy_t_srf_energy_state(in, out, iblk, ics, ice)
!vg 
!vg     TYPE(t_srf_energy_state), INTENT(in)  :: in
!vg     TYPE(t_srf_energy_state), INTENT(out) :: out
!vg     INTEGER,                   INTENT(in)  :: iblk, ics, ice
!vg 
!vg     out%t_srf         (ics:ice,    iblk) = in%t_srf         (ics:ice,    iblk)
!vg     out%t_srf_old     (ics:ice,    iblk) = in%t_srf_old     (ics:ice,    iblk)
!vg     out%t_srf_unfilt  (ics:ice,    iblk) = in%t_srf_unfilt  (ics:ice,    iblk)
!vg     out%t_eff_srf     (ics:ice,    iblk) = in%t_eff_srf     (ics:ice,    iblk)
!vg     out%s_srf         (ics:ice,    iblk) = in%s_srf         (ics:ice,    iblk)
!vg !!$    out%s_srf_old     (ics:ice,    iblk) = in%s_srf_old     (ics:ice,    iblk)
!vg     out%qsat_srf      (ics:ice,    iblk) = in%qsat_srf      (ics:ice,    iblk)
!vg !!$    out%qsat_srf_old  (ics:ice,    iblk) = in%qsat_srf_old  (ics:ice,    iblk)
!vg     out%tte_corr      (ics:ice,    iblk) = in%tte_corr      (ics:ice,    iblk)
!vg     out%fact_q_air    (ics:ice,    iblk) = in%fact_q_air    (ics:ice,    iblk)
!vg     out%fact_qsat_srf (ics:ice,    iblk) = in%fact_qsat_srf (ics:ice,    iblk)
!vg     out%fact_qsat_trans_srf (ics:ice,    iblk) = in%fact_qsat_trans_srf (ics:ice,    iblk)
!vg 
!vg   END SUBROUTINE copy_t_srf_energy_state

!vg   SUBROUTINE add_t_srf_energy_state(in1, in2, out, iblk, ics, ice)
!vg 
!vg     TYPE(t_srf_energy_state), INTENT(in)  :: in1, in2
!vg     TYPE(t_srf_energy_state), INTENT(out) :: out
!vg     INTEGER,                   INTENT(in)  :: iblk, ics, ice
!vg 
!vg     out%t_srf         (ics:ice,    iblk) = in1%t_srf         (ics:ice,    iblk) + in2%t_srf         (ics:ice,    iblk)
!vg     out%t_srf_old     (ics:ice,    iblk) = in1%t_srf_old     (ics:ice,    iblk) + in2%t_srf_old     (ics:ice,    iblk)
!vg     out%t_srf_unfilt  (ics:ice,    iblk) = in1%t_srf_unfilt  (ics:ice,    iblk) + in2%t_srf_unfilt  (ics:ice,    iblk)
!vg     out%t_eff_srf     (ics:ice,    iblk) = in1%t_eff_srf     (ics:ice,    iblk) + in2%t_eff_srf     (ics:ice,    iblk)
!vg     out%s_srf         (ics:ice,    iblk) = in1%s_srf         (ics:ice,    iblk) + in2%s_srf         (ics:ice,    iblk)
!vg !!$    out%s_srf_old     (ics:ice,    iblk) = in1%s_srf_old     (ics:ice,    iblk) + in2%s_srf_old     (ics:ice,    iblk)
!vg     out%qsat_srf      (ics:ice,    iblk) = in1%qsat_srf      (ics:ice,    iblk) + in2%qsat_srf      (ics:ice,    iblk)
!vg !!$    out%qsat_srf_old  (ics:ice,    iblk) = in1%qsat_srf_old  (ics:ice,    iblk) + in2%qsat_srf_old  (ics:ice,    iblk)
!vg     out%tte_corr      (ics:ice,    iblk) = in1%tte_corr      (ics:ice,    iblk) + in2%tte_corr      (ics:ice,    iblk)
!vg     out%fact_q_air    (ics:ice,    iblk) = in1%fact_q_air    (ics:ice,    iblk) + in2%fact_q_air    (ics:ice,    iblk)
!vg     out%fact_qsat_srf (ics:ice,    iblk) = in1%fact_qsat_srf (ics:ice,    iblk) + in2%fact_qsat_srf (ics:ice,    iblk)
!vg 
!vg   END SUBROUTINE add_t_srf_energy_state

#endif
END MODULE mo_hd_memory_class
