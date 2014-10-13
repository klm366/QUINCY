!>
!! @brief Contains the memory class for the surface.
!!
!! @author
!!  Reiner Schnur, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Reiner Schnur (2013-04-04)
!!
MODULE mo_srf_memory_class
#ifndef __NO_JSBACH__

  USE mo_kind, ONLY: wp

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: t_srf_memory, t_srf_phys_state, &
            t_srf_energy_state, t_srf_water_state, &
            t_srf_energy_flux,  t_srf_water_flux, &
            copy_srf_mem, add_srf_mem

  TYPE t_srf_phys_state
    REAL(wp), POINTER ::           &
      glacier      (:,:),    & !< Glacier mask
      oro_stddev   (:,:),    & !< Standard deviation of orography
      rough_m_srf  (:,:),    & !< Surface roughness length for momentum
      rough_h_srf  (:,:),    & !< Surface roughness length for heat
      albedo_srf   (:,:),    & !< (Background) albedo of tile
      alb_vis_srf  (:,:),    & !< Albedo in the visible range
      alb_nir_srf  (:,:)       !< Albedo in the NIR range
  END TYPE t_srf_phys_state

  TYPE t_srf_energy_state
    REAL(wp), POINTER :: &
      t_srf            (:,:), &
      t_srf_old        (:,:), &
      t_srf_unfilt     (:,:), &
      t_srf_unfilt_old (:,:), &
      t_rad_srf        (:,:), &
!!$      t_eff_srf        (:,:), &
      s_srf            (:,:), &  !< Dry static energy at surface (new timestep)
!!$      s_srf_old        (:,:), &  !< Dry static energy at surface (old timestep)
      qsat_srf         (:,:), &
!!$      qsat_srf_old     (:,:), &
      tte_corr         (:,:), &  !< Correction term of tte for snow melt
      fact_q_air       (:,:), &  !< cair in echam6
      fact_qsat_srf    (:,:), &  !< csat in echam6
      fact_qsat_trans_srf (:,:)  !< csat_transpiration in echam6/jsbach
  END TYPE t_srf_energy_state

  TYPE t_srf_water_state
    REAL(wp), POINTER :: &
      wfract_srf    (:,:),   & !< Wet (skin reservoir) fraction of tile []
      wsr_srf       (:,:),   & !< Water content in skin reservoir [m]
      sfract_srf    (:,:),   & !< Snow fraction of tile []
      wsn_srf       (:,:),   & !< Water content of snow reservoir [m water equivalent]
      wsn_glacier   (:,:)      !< Water content of snow on glacier [m water equivalent]
  END TYPE t_srf_water_state

  TYPE t_srf_energy_flux
    REAL(wp), POINTER :: &
      latent_hflx      (:,:),   & !< latent heat flux
      sensible_hflx    (:,:),   & !< sensible heat flux
      swrad_srf_down   (:,:),   & !< net shortwave radiation at surface
      lwrad_srf_net    (:,:),   & !< net longwave radiation at surface
      rad_srf_net      (:,:)      !< net radiation at surface
  END TYPE t_srf_energy_flux

  TYPE t_srf_water_flux
    REAL(wp), POINTER :: &
      evapopot         (:,:),   & !< Potential evaporation
      evapotrans       (:,:),   & !< Evapotranspiration
      trans            (:,:)      !< Transpiration
  END TYPE t_srf_water_flux

  TYPE t_srf_memory
    TYPE(t_srf_phys_state),   POINTER :: phys_state   => NULL()
    TYPE(t_srf_energy_state), POINTER :: energy_state => NULL()
    TYPE(t_srf_water_state),  POINTER :: water_state  => NULL()
    TYPE(t_srf_energy_flux),  POINTER :: energy_flux  => NULL()
    TYPE(t_srf_water_flux),   POINTER :: water_flux   => NULL()
  END TYPE t_srf_memory

  ! Note: Would be nicer to be able to overload the "=" operator, but that would be more
  ! complicated because we only want to copy parts of the arrays (ics:ice,iblk)
  INTERFACE copy_srf_mem
    MODULE PROCEDURE copy_t_srf_energy_state
  END INTERFACE copy_srf_mem

  INTERFACE add_srf_mem
    MODULE PROCEDURE add_t_srf_energy_state
  END INTERFACE add_srf_mem

  CHARACTER(len=*), PARAMETER :: modname = 'mo_srf_memory_class'

CONTAINS

  SUBROUTINE copy_t_srf_energy_state(in, out, iblk, ics, ice)

    TYPE(t_srf_energy_state), INTENT(in)  :: in
    TYPE(t_srf_energy_state), INTENT(out) :: out
    INTEGER,                   INTENT(in)  :: iblk, ics, ice

    out%t_srf         (ics:ice,    iblk) = in%t_srf         (ics:ice,    iblk)
    out%t_srf_old     (ics:ice,    iblk) = in%t_srf_old     (ics:ice,    iblk)
    out%t_srf_unfilt  (ics:ice,    iblk) = in%t_srf_unfilt  (ics:ice,    iblk)
    out%t_rad_srf     (ics:ice,    iblk) = in%t_rad_srf     (ics:ice,    iblk)
!!$    out%t_eff_srf     (ics:ice,    iblk) = in%t_eff_srf     (ics:ice,    iblk)
    out%s_srf         (ics:ice,    iblk) = in%s_srf         (ics:ice,    iblk)
!!$    out%s_srf_old     (ics:ice,    iblk) = in%s_srf_old     (ics:ice,    iblk)
    out%qsat_srf      (ics:ice,    iblk) = in%qsat_srf      (ics:ice,    iblk)
!!$    out%qsat_srf_old  (ics:ice,    iblk) = in%qsat_srf_old  (ics:ice,    iblk)
    out%tte_corr      (ics:ice,    iblk) = in%tte_corr      (ics:ice,    iblk)
    out%fact_q_air    (ics:ice,    iblk) = in%fact_q_air    (ics:ice,    iblk)
    out%fact_qsat_srf (ics:ice,    iblk) = in%fact_qsat_srf (ics:ice,    iblk)
    out%fact_qsat_trans_srf (ics:ice,    iblk) = in%fact_qsat_trans_srf (ics:ice,    iblk)

  END SUBROUTINE copy_t_srf_energy_state

  SUBROUTINE add_t_srf_energy_state(in1, in2, out, iblk, ics, ice)

    TYPE(t_srf_energy_state), INTENT(in)  :: in1, in2
    TYPE(t_srf_energy_state), INTENT(out) :: out
    INTEGER,                   INTENT(in)  :: iblk, ics, ice

    out%t_srf         (ics:ice,    iblk) = in1%t_srf         (ics:ice,    iblk) + in2%t_srf         (ics:ice,    iblk)
    out%t_srf_old     (ics:ice,    iblk) = in1%t_srf_old     (ics:ice,    iblk) + in2%t_srf_old     (ics:ice,    iblk)
    out%t_srf_unfilt  (ics:ice,    iblk) = in1%t_srf_unfilt  (ics:ice,    iblk) + in2%t_srf_unfilt  (ics:ice,    iblk)
    out%t_rad_srf     (ics:ice,    iblk) = in1%t_rad_srf     (ics:ice,    iblk) + in2%t_rad_srf     (ics:ice,    iblk)
!!$    out%t_eff_srf     (ics:ice,    iblk) = in1%t_eff_srf     (ics:ice,    iblk) + in2%t_eff_srf     (ics:ice,    iblk)
    out%s_srf         (ics:ice,    iblk) = in1%s_srf         (ics:ice,    iblk) + in2%s_srf         (ics:ice,    iblk)
!!$    out%s_srf_old     (ics:ice,    iblk) = in1%s_srf_old     (ics:ice,    iblk) + in2%s_srf_old     (ics:ice,    iblk)
    out%qsat_srf      (ics:ice,    iblk) = in1%qsat_srf      (ics:ice,    iblk) + in2%qsat_srf      (ics:ice,    iblk)
!!$    out%qsat_srf_old  (ics:ice,    iblk) = in1%qsat_srf_old  (ics:ice,    iblk) + in2%qsat_srf_old  (ics:ice,    iblk)
    out%tte_corr      (ics:ice,    iblk) = in1%tte_corr      (ics:ice,    iblk) + in2%tte_corr      (ics:ice,    iblk)
    out%fact_q_air    (ics:ice,    iblk) = in1%fact_q_air    (ics:ice,    iblk) + in2%fact_q_air    (ics:ice,    iblk)
    out%fact_qsat_srf (ics:ice,    iblk) = in1%fact_qsat_srf (ics:ice,    iblk) + in2%fact_qsat_srf (ics:ice,    iblk)
    out%fact_qsat_trans_srf (ics:ice,  iblk) = in1%fact_qsat_trans_srf (ics:ice,  iblk) + in2%fact_qsat_trans_srf (ics:ice,  iblk)

  END SUBROUTINE add_t_srf_energy_state

#endif
END MODULE mo_srf_memory_class
