!> QUINCY calculate vegetation dynamics
!>
!> ICON-Land
!>
!> ---------------------------------------
!> Copyright (C) 2013-2024, MPI-M, MPI-BGC
!>
!> Contact: icon-model.org
!> Authors: AUTHORS.md
!> See LICENSES/ for license information
!> SPDX-License-Identifier: BSD-3-Clause
!> ---------------------------------------
!>
!> For more information on the QUINCY model see: <https://doi.org/10.17871/quincy-model-2019>
!>
!> @brief module to determine the spitfire fire ignition
!>
!> @author
!>  Albert Jornet Puig, University of Leipzig / MPI-BGC Jena
!>  Ana Bastos, University of Leipzig / MPI-BGC Jena
!>
MODULE mo_q_dfire_sf_ignition_class
#ifndef __NO_JSBACH__

  USE mo_kind,                ONLY: wp
  USE mo_jsb_control,         ONLY: debug_on
  USE mo_exception,           ONLY: message

  USE mo_q_dfire_constants,   ONLY: df_sf_lign_cg, df_sf_lign_eff

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: t_q_sf_ignitions, calc_q_sf_ignition

  CHARACTER(len=*), PARAMETER :: modname = 'mo_q_dfire_ignition_class'


  TYPE t_q_sf_ignitions
    REAL(wp), ALLOCATABLE :: human(:)
    REAL(wp), ALLOCATABLE :: lightning(:)
    REAL(wp), ALLOCATABLE :: total(:)
  END TYPE t_q_sf_ignitions


  INTERFACE t_q_sf_ignitions
    MODULE PROCEDURE :: create_t_q_sf_ignitions
  END INTERFACE


CONTAINS


  FUNCTION create_t_q_sf_ignitions(nc) RESULT (sf_ignitions)
    INTEGER, INTENT(in) :: nc
    TYPE(t_q_sf_ignitions) :: sf_ignitions
    ALLOCATE(sf_ignitions%human(nc))
    ALLOCATE(sf_ignitions%lightning(nc))
    ALLOCATE(sf_ignitions%total(nc))

    sf_ignitions%human = 0._wp
    sf_ignitions%lightning = 0._wp
    sf_ignitions%total = 0._wp
  END FUNCTION create_t_q_sf_ignitions


  ! ======================================================================================================= !
  !>
  !> Calculation of total ignitions (lightning + human) based on Section 2.2.1 from Thonicke et al. (2010)
  !>
  FUNCTION calc_q_sf_ignition(nc, dens_hum_pop, a_nd, n_lightning_ig, max_dens_hum_pop) RESULT (sf_ignitions)
    INTEGER, INTENT(in)    :: nc                     !< dimensions
    REAL(wp), INTENT(in)   :: dens_hum_pop(nc)       !< average human population density in the grid cell
                                                     !< [individuals km-2]
    REAL(wp), INTENT(in)   :: max_dens_hum_pop       !< max human population density allowed per gridcell
                                                     !< [individuals km-2]
    REAL(wp), INTENT(in)   :: n_lightning_ig(nc)     !< average number of flashes per day [ha-1 d-1]
    REAL(wp), INTENT(in)   :: a_nd(nc)               !< propensity to people to produce ignitions events
                                                     !< [ignitions individual-1 d-1]
    TYPE(t_q_sf_ignitions) :: sf_ignitions           !< (output)
    ! ----------------------------------------------------------------------------------------------------- !
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_q_sf_ignition'
    REAL(wp)             :: sel_dens_hum_pop(nc)     !< selected average human population density [indiv km-2]
    REAL(wp)             :: k_pd(nc)                 !< Spatial dependence of human-caused ignitions [-]
    INTEGER              :: ic

    !< Initialize output to 0
    sf_ignitions = t_q_sf_ignitions(nc)

    ! ----------------------------------------------------------------------------------------------------- !
    ! Go science
    ! ----------------------------------------------------------------------------------------------------- !
    sel_dens_hum_pop = dens_hum_pop
    DO ic=1, nc
      IF (sel_dens_hum_pop(ic) > max_dens_hum_pop) sel_dens_hum_pop(ic) = max_dens_hum_pop
    ENDDO

    !< Calculate how likely a fire is about to ignite due to human population
    !< Section 2.2.1, equations 4 from Thonicke et al. (2010)
    k_pd = 30.0_wp * exp(-0.5_wp * sqrt(sel_dens_hum_pop))

    !< Calculate how likely a fire is about to ignite due to human population
    !< Section 2.2.1, equation 3 from Thonicke et al. (2010)
    sf_ignitions%human = sel_dens_hum_pop * k_pd * a_nd / 100._wp

    sf_ignitions%lightning = n_lightning_ig * df_sf_lign_cg * df_sf_lign_eff
    !< Section 2.2, last paragraph (below eq. 2 from Thonicke et al. (2010))
    sf_ignitions%total = sf_ignitions%human + sf_ignitions%lightning
  END FUNCTION calc_q_sf_ignition

#endif
END MODULE mo_q_dfire_sf_ignition_class
