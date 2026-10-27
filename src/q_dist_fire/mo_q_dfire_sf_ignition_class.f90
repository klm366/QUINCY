!>
!> ICON-Land
!>
!> ---------------------------------------
!> Copyright (C) 2013-2026, MPI-M, MPI-BGC
!>
!> Contact: icon-model.org
!> Authors: AUTHORS.md
!> See LICENSES/ for license information
!> SPDX-License-Identifier: BSD-3-Clause
!> ---------------------------------------
!>
!>#### routine to calculate fire ignitions
!>
MODULE mo_q_dfire_sf_ignition_class
#ifndef __NO_JSBACH__

  USE mo_kind,                ONLY: lp
  USE mo_jsb_control,         ONLY: debug_on
  USE mo_exception,           ONLY: message

  USE mo_q_dfire_sf_constants, ONLY: df_sf_lign_cg, df_sf_lign_eff
  USE mo_q_dfire_sf_classes,  ONLY: t_sf_ignitions, t_sf_config

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: calc_sf_ignitions, calc_sf_ignitions_main

  CHARACTER(len=*), PARAMETER :: modname = 'mo_q_dfire_ignition_class'

CONTAINS

  !>
  !>  Wrapper for calc_sf_ignitions
  !>
  !> Call main ignition method and applies time scaling to ignition events.
  !>
  SUBROUTINE calc_sf_ignitions_main(sf_config, sf_ignitions)
    TYPE(t_sf_config), INTENT(in) :: sf_config
    TYPE(t_sf_ignitions), INTENT(inout) :: sf_ignitions

    CALL calc_sf_ignitions(sf_config%nc, sf_config%const_max_density_human_pop, &
                          sf_ignitions%pop_dens, &
                          sf_ignitions%a_nd, &
                          sf_ignitions%lightning_freq,  &
                          sf_ignitions%e_n_humans, sf_ignitions%e_n_lightning, sf_ignitions%e_n_ignitions)

    !< Daily to ts
    sf_ignitions%e_n_ignitions = sf_ignitions%e_n_ignitions * (sf_config%dtime / sf_config%one_day)
    sf_ignitions%e_n_humans = sf_ignitions%e_n_humans * (sf_config%dtime / sf_config%one_day)
    sf_ignitions%e_n_lightning = sf_ignitions%e_n_lightning * (sf_config%dtime / sf_config%one_day)
  END SUBROUTINE calc_sf_ignitions_main


  ! ======================================================================================================= !
  !>
  !> Calculation of total ignitions (lightning + human) based on Section 2.2.1 from Thonicke et al. (2010)
  !>
  SUBROUTINE calc_sf_ignitions(nc, max_dens_hum_pop, dens_hum_pop, a_nd, lightning_freq, e_n_humans, e_n_lightnings, e_n_ignitions)
    INTEGER, INTENT(in)    :: nc                     !< dimensions
    REAL(lp), INTENT(in)   :: dens_hum_pop(nc)       !< average human population density in the grid cell
                                                     !< [individuals km-2]
    REAL(lp), INTENT(in)   :: max_dens_hum_pop       !< max human population density allowed per gridcell
                                                     !< [individuals km-2]
    REAL(lp), INTENT(in)   :: lightning_freq(nc)     !< average number of flashes per day [ha-1 d-1]
    REAL(lp), INTENT(in)   :: a_nd(nc)               !< propensity to people to produce ignitions events
                                                     !< [ignitions individual-1 d-1]
    REAL(lp), INTENT(out)  :: e_n_humans(nc)         !< Expected number of human-caused ignitions [ha-1 ts-1]
    REAL(lp), INTENT(out)  :: e_n_lightnings(nc)     !<
    REAL(lp), INTENT(out)  :: e_n_ignitions(nc)      !<
    ! ----------------------------------------------------------------------------------------------------- !
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_sf_ignition'
    REAL(lp)             :: sel_dens_hum_pop(nc)     !< selected average human population density [indiv km-2]
    REAL(lp)             :: k_pd(nc)                 !< Spatial dependence of human-caused ignitions [-]
    INTEGER              :: ic

    ! ----------------------------------------------------------------------------------------------------- !
    ! Go science
    ! ----------------------------------------------------------------------------------------------------- !
    sel_dens_hum_pop = dens_hum_pop
    DO ic=1, nc
      IF (sel_dens_hum_pop(ic) > max_dens_hum_pop) sel_dens_hum_pop(ic) = max_dens_hum_pop
    ENDDO

    !< Calculate how likely a fire is about to ignite due to human population
    !< Section 2.2.1, equations 4 from Thonicke et al. (2010)
    k_pd = 30.0_lp * exp(-0.5_lp * sqrt(sel_dens_hum_pop))

    !< Calculate how likely a fire is about to ignite due to human population
    !< Section 2.2.1, equation 3 from Thonicke et al. (2010)
    e_n_humans = sel_dens_hum_pop * k_pd * a_nd / 100._lp

    e_n_lightnings = lightning_freq * df_sf_lign_cg * df_sf_lign_eff
    !< Section 2.2, last paragraph (below eq. 2 from Thonicke et al. (2010))
    e_n_ignitions = e_n_humans + e_n_lightnings
  END SUBROUTINE calc_sf_ignitions

#endif
END MODULE mo_q_dfire_sf_ignition_class
