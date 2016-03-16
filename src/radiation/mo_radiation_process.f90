!>
!! @brief Contains the routines for the radiation processes
!!
!! @author
!!  Rainer Schneck, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Rainer Schneck (2015-04-27)
!!
MODULE mo_radiation_process
#ifndef __NO_JSBACH__

   USE mo_kind,      ONLY: wp
   USE mo_exception, ONLY: message

   IMPLICIT NONE
   PRIVATE

   PUBLIC :: calc_radiation_surface_net

   CHARACTER(len=*), PARAMETER :: modname = 'mo_radiation_process'

CONTAINS

    ELEMENTAL SUBROUTINE calc_radiation_surface_net(swvis_srf_down,    & ! input
                                                    swnir_srf_down,    & ! input
                                                    alb_vis_srf,       & ! input
                                                    alb_nir_srf,       & ! input
                                                    lw_srf_down,       & ! input
                                                    t_srf,             & ! input
                                                    swvis_srf_net,     & ! output
                                                    swnir_srf_net,     & ! output
                                                    sw_srf_net,        & ! output
                                                    lw_srf_net,        & ! output
                                                    rad_srf_net )        ! output

      !-----------------------------------------------------------------------
      !  DECLARATIONS
      USE mo_phy_schemes,          ONLY: lwnet_from_lwdown

      !-----------------------------------------------------------------------
      !  ARGUMENTS
      REAL(wp),  INTENT(in) ::                 &
                               swvis_srf_down, &
                               swnir_srf_down, &
                               alb_vis_srf,    &
                               alb_nir_srf,    &
                               lw_srf_down,    &
                               t_srf

      REAL(wp),  INTENT(out) ::               &
                               swvis_srf_net, &
                               swnir_srf_net, &
                               sw_srf_net,    &
                               lw_srf_net,    &
                               rad_srf_net

      ! Compute net SW radiation from downward SW and albedo
      swvis_srf_net = swvis_srf_down * (1._wp - alb_vis_srf)
      swnir_srf_net = swnir_srf_down * (1._wp - alb_nir_srf)
      sw_srf_net    = swvis_srf_net + swnir_srf_net

      ! Compute LW net radiation from incoming and the thermal radiation from the surface
      lw_srf_net = lwnet_from_lwdown(lw_srf_down, t_srf)

      ! Compute net radiation
      rad_srf_net = sw_srf_net + lw_srf_net

   END SUBROUTINE calc_radiation_surface_net

#endif
END MODULE mo_radiation_process
