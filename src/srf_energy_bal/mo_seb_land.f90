!>
!! @brief Contains the routines for the surface energy balance on LAND lct_type.
!!
!! @author
!!  Reiner Schnur, MPI-M Hamburg
!!
!! @par Revision History
!! First version in mo_srf_interface.f90 by Reiner Schnur (2013-04-04)
!! Separated out the surface energy balance for LAND part by Reiner Schnur (2016-03-22)
!!
MODULE mo_seb_land
#ifndef __NO_JSBACH__

#include "dsl_jsbach.inc"

  USE mo_kind,      ONLY: wp
  USE mo_exception, ONLY: message

  USE mo_jsb_tile_class,     ONLY: t_jsb_tile_abstract
  USE mo_jsb_config_class,   ONLY: t_jsb_config, t_jsb_config_p
  USE mo_jsb_task_class,     ONLY: t_jsb_task_options
  USE mo_jsb_lct_class,      ONLY: GLACIER_TYPE, Contains_lct
  USE mo_jsb_control,        ONLY: debug_on

  dsl4jsb_Use_processes SEB, RAD, ATM2LAND, SOIL, HYDRO
  dsl4jsb_Use_config(SEB)
  dsl4jsb_Use_config(HYDRO)

  USE mo_atmland_memory_class, ONLY: t_atm2land_memory
!!$  dsl4jsb_Use_memory(ATM2LAND)
  dsl4jsb_Use_memory(SEB)
  dsl4jsb_Use_memory(RAD)
  dsl4jsb_Use_memory(SOIL)
  dsl4jsb_Use_memory(HYDRO)

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: update_surface_energy_land, update_asselin_land

  CHARACTER(len=*), PARAMETER :: modname = 'mo_seb_land'

CONTAINS
  !
  ! ================================================================================================================================
  !
  SUBROUTINE update_surface_energy_land(tile, config, options)

    USE mo_phy_schemes,            ONLY: qsat, q_effective, surface_dry_static_energy, heat_transfer_coef
    USE mo_echam_vdiff_params,     ONLY: tpfac2, tpfac3
    USE mo_jsb_physical_constants, ONLY: tmelt
    USE mo_jsb_time,               ONLY: get_asselin_coef

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    CLASS(t_jsb_config_p),      INTENT(in)    :: config(:)
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    ! Local variables
    !
    dsl4jsb_Def_config(SEB)
    dsl4jsb_Def_config(HYDRO)
    dsl4jsb_Def_memory(SEB)
    dsl4jsb_Def_memory(RAD)
    dsl4jsb_Def_memory(SOIL)
    dsl4jsb_Def_memory(HYDRO)
    dsl4jsb_Def_memory(ATM2LAND)

    LOGICAL                       :: lstart

    ! Pointers to variables in memory
    dsl4jsb_Real2D_onChunk :: t
    dsl4jsb_Real2D_onChunk :: t_old
    dsl4jsb_Real2D_onChunk :: t_unfilt
    dsl4jsb_Real2D_onChunk :: t_unfilt_old
    dsl4jsb_Real2D_onChunk :: t_eff
    dsl4jsb_Real2D_onChunk :: qsat_star
    dsl4jsb_Real2D_onChunk :: s_star
    dsl4jsb_Real2D_onChunk :: heat_cap
    dsl4jsb_Real2D_onChunk :: forc_hflx
    dsl4jsb_Real2D_onChunk :: press_srf
    dsl4jsb_Real2D_onChunk :: q_air
    dsl4jsb_Real2D_onChunk :: fact_qsat_srf
    dsl4jsb_Real2D_onChunk :: fact_q_air
    dsl4jsb_Real2D_onChunk :: drag_srf
    dsl4jsb_Real2D_onChunk :: rad_srf_net
    dsl4jsb_Real2D_onChunk :: grnd_hflx
    dsl4jsb_Real2D_onChunk :: hcap_grnd
    dsl4jsb_Real2D_onChunk :: w_snow_srf
    dsl4jsb_Real2D_onChunk :: fract_snow_srf
    dsl4jsb_Real2D_onChunk :: t_acoef
    dsl4jsb_Real2D_onChunk :: t_bcoef
    dsl4jsb_Real2D_onChunk :: q_acoef
    dsl4jsb_Real2D_onChunk :: q_bcoef

    ! Locally allocated vectors
    !
    REAL(wp) ::             &
      qsat_srf_old(options%nc),     &
      s_old(options%nc),        &
      dQdT(options%nc),             & !< Sensitivity of saturated surface specific humidity to temperature
      t2s_conv(options%nc),         & !< Conversion factor from temperature to dry static energy (C_pd * (1+(delta-1)*q_v))
      heat_tcoef(options%nc)          !< Heat transfer coefficient (rho*C_h*|v|)
    LOGICAL :: &
      is_glacier(options%nc)

    INTEGER :: iblk, ics, ice, nc
    REAL(wp) :: steplen
    REAL(wp) :: eps

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':update_surface_energy_land'

    iblk = options%iblk
    ics  = options%ics
    ice  = options%ice
    nc   = options%nc
    steplen = options%steplen

    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    dsl4jsb_Get_config(SEB)
    dsl4jsb_Get_config(HYDRO)

    ! Get reference to variables for current block
    !
    dsl4jsb_Get_memory(SEB)
    dsl4jsb_Get_memory(RAD)
    dsl4jsb_Get_memory(SOIL)
    dsl4jsb_Get_memory(HYDRO)
    dsl4jsb_Get_memory(ATM2LAND)

    dsl4jsb_Get_var2D_onChunk(ATM2LAND, t_acoef)        ! IN
    dsl4jsb_Get_var2D_onChunk(ATM2LAND, t_bcoef)        ! IN
    dsl4jsb_Get_var2D_onChunk(ATM2LAND, q_acoef)        ! IN
    dsl4jsb_Get_var2D_onChunk(ATM2LAND, q_bcoef)        ! IN
    dsl4jsb_Get_var2D_onChunk(ATM2LAND, press_srf)      ! IN
    dsl4jsb_Get_var2D_onChunk(ATM2LAND, q_air)          ! IN
    dsl4jsb_Get_var2D_onChunk(ATM2LAND, drag_srf)       ! IN

    dsl4jsb_Get_var2D_onChunk(SEB,      t)              ! OUT
    dsl4jsb_Get_var2D_onChunk(SEB,      t_old)          ! OUT
    dsl4jsb_Get_var2D_onChunk(SEB,      t_unfilt)       ! OUT
    dsl4jsb_Get_var2D_onChunk(SEB,      t_unfilt_old)   ! OUT
    dsl4jsb_Get_var2D_onChunk(SEB,      t_eff)          ! OUT
    dsl4jsb_Get_var2D_onChunk(SEB,      qsat_star)      ! OUT
    dsl4jsb_Get_var2D_onChunk(SEB,      s_star)         ! OUT
    dsl4jsb_Get_var2D_onChunk(SEB,      heat_cap)       ! OUT
    dsl4jsb_Get_var2D_onChunk(SEB,      forc_hflx)      ! OUT
    dsl4jsb_Get_var2D_onChunk(SEB,      fact_q_air)     ! OUT
    dsl4jsb_Get_var2D_onChunk(SEB,      fact_qsat_srf)  ! OUT
    dsl4jsb_Get_var2D_onChunk(RAD,      rad_srf_net)    ! OUT
    dsl4jsb_Get_var2D_onChunk(HYDRO,    w_snow_srf)     ! OUT
    dsl4jsb_Get_var2D_onChunk(HYDRO,    fract_snow_srf) ! OUT

    dsl4jsb_Get_var2D_onChunk( SOIL,    grnd_hflx)      ! IN
    dsl4jsb_Get_var2D_onChunk( SOIL,    hcap_grnd)      ! IN

    IF (Contains_lct(tile%lcts, GLACIER_TYPE)) THEN
      is_glacier(:) = .TRUE.
    ELSE
      is_glacier(:) = .FALSE.
    END IF

    IF (ALL(heat_cap(:) < EPSILON(1._wp))) THEN            ! Start of experiment
      lstart = .TRUE.
    ELSE
      lstart = .FALSE.
    END IF

    ! Save old surface temperature and saturation spec. humidity
!!$    t_srf_unfilt_old(:) = t_srf_unfilt(:)
!!$    t_srf_old   (:) = t_srf(:)
    ! Temporary fix for temperatures over ocean getting too warm or too cold, leading to crash in convect_tables
    t_unfilt_old(:) = MAX(MIN(t_unfilt(:), 380._wp), 100._wp)
    t_old(:)        = MAX(MIN(t(:),        380._wp), 100._wp)
    qsat_srf_old(:) = qsat(t_old(:), press_srf(:))

    ! The heat capacity of the surface layer is currently taken as the heat capacity of the upper soil layer
    heat_cap(:) = hcap_grnd(:)

    ! The last term of the rhs of the energy balance equation is the conductive heat flux from the ground below
    forc_hflx(:) = grnd_hflx(:)

    ! Calculate surface saturation specific humidity

    ! Update surface temperature
    !
    ! 1. Compute sensitivity of surface saturation specific humidity to temperature
    dQdT(:) = ( qsat(t(:) + 0.001_wp, press_srf(:)) - qsat_srf_old(:) ) * 1000._wp

    ! 2. Old dry static energy
    s_old(:) = surface_dry_static_energy( t_old(:), &
                                            & q_effective(qsat_srf_old(:), q_air(:), fact_qsat_srf(:), fact_q_air(:)))

    t2s_conv(:) = s_old(:) / t_old(:)  !

    heat_tcoef(:) = heat_transfer_coef(drag_srf(:), steplen)

    ! Compute new surface temperature and saturation humidity
    IF (lstart) THEN
      s_star(:) = s_old(:)
    ELSE
      CALL surface_temp_implicit( nc, steplen, &
        t2s_conv(:), &
        t_acoef(:), t_bcoef(:), q_acoef(:), q_bcoef(:), &
        s_old(:), qsat_srf_old(:), dQdT(:), &
        rad_srf_net(:), forc_hflx(:), heat_tcoef(:), &
        fact_q_air(:), fact_qsat_srf(:), &
        fract_snow_srf(:), heat_cap(:), &
        & s_star(:))
    END IF

    ! New unfiltered surface temperature (\tilde(X)^(t+1)) (see below Eq. 3.3.1.5 in ECHAM3 manual, alpha = 1/tpfac2)
    t_unfilt(:) = ( tpfac2 * s_star(:) + tpfac3 * s_old(:) ) / t2s_conv(:)

    ! Correct new dry static energy for snowmelt
    WHERE (w_snow_srf(:) > dsl4jsb_Config(HYDRO)%wsn_critical .OR. is_glacier(:))
      s_star(:) = MIN(s_star(:) / t2s_conv(:), tmelt) * t2s_conv(:)
    END WHERE

    ! Compute new surface qsat (after correction of dry static energy for snow melt!)
    IF (lstart) THEN
      qsat_star(:) = qsat_srf_old(:)
    ELSE
      qsat_star(:) = qsat_srf_old(:) + dQdT(:) * (s_star(:) - s_old(:)) / t2s_conv(:)
    END IF

    ! 'Effective temperature' used in radheat of the atmosphere model (see Eq. 6.3 in ECHAM5 manual)
    t_eff(:) = (t_old(:)**3 * (4._wp * s_star(:)/t2s_conv(:) - 3._wp * t_old(:))) ** 0.25_wp

    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE update_surface_energy_land

  SUBROUTINE surface_temp_implicit(klon, steplen,             &
    &            pcp,                                         &
    &            pescoe, pfscoe, peqcoe, pfqcoe,              &
    &            psold, pqsold, pdqsold,                      &
    &            pnetrad, pgrdfl,                             &
    &            pcfh, pcair, pcsat, pfracsu, pgrdcap,        &
    &            psnew)

    USE mo_echam_vdiff_params,     ONLY: cvdifts   ! TBD: change to 1. for standalone JSBACH
    USE mo_jsb_physical_constants, ONLY: &
      & cpd, &
      & cpvd1, &
      & zemiss_def, &
      & stbo, &
      & alv,  &
      & als


    INTEGER,  INTENT(in)    :: klon
    REAL(wp), INTENT(in)    :: steplen
    REAL(wp), INTENT(in)    :: pcp(klon), pfscoe(klon), pescoe(klon), pfqcoe(klon), peqcoe(klon)
    REAL(wp), INTENT(in)    :: psold(klon), pqsold(klon), pdqsold(klon)
    REAL(wp), INTENT(in)    :: pnetrad(klon), pgrdfl(klon)
    REAL(wp), INTENT(in)    :: pcfh(klon), pcair(klon), pcsat(klon), pfracsu(klon)
    REAL(wp), INTENT(in)    :: pgrdcap(klon)
    REAL(wp), INTENT(out)   :: psnew(klon)
    !
    REAL(wp) :: zcolin(klon), zcohfl(klon), zcoind(klon), zicp(klon), zca(klon), zcs(klon)
    REAL(wp) :: pdt
    REAL(wp) :: pc16

    pdt     = cvdifts * steplen
    pc16    = cpd * cpvd1

    !
    zicp(:) = 1._wp / pcp(:)
    !
    zca(:)    = als * pfracsu(:) +  alv * (pcair(:) - pfracsu(:))
    zcs(:)    = als * pfracsu(:) +  alv * (pcsat(:) - pfracsu(:))
    !
    zcolin(:) = pgrdcap(:)*zicp(:) +                                                              &
              &      pdt * (zicp(:) * 4._wp * zemiss_def * stbo * ((zicp(:) * psold(:))**3)  -    &
              &      pcfh(:) * (zca(:) * peqcoe(:) - zcs(:)) *                                    &
              &      zicp(:) * pdqsold(:))
    !
    zcohfl(:) = -pdt * pcfh(:) * (pescoe(:) - 1._wp)
    !
    zcoind(:) = pdt * (pnetrad(:) + pcfh(:) * pfscoe(:) +  pcfh(:) *                              &
              &      ((zca(:) * peqcoe(:) - zcs(:)) * pqsold(:) + zca(:) * pfqcoe(:)) + pgrdfl(:))
    !
    psnew(:)  = (zcolin(:) * psold(:) + zcoind(:)) / (zcolin(:) + zcohfl(:))

  END SUBROUTINE surface_temp_implicit
  !
  ! ================================================================================================================================
  !
  SUBROUTINE update_asselin_land(tile, config, options)

    USE mo_jsb_time,               ONLY: get_asselin_coef

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    CLASS(t_jsb_config_p),      INTENT(in)    :: config(:)
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    ! Local variables
    !
    dsl4jsb_Def_config(SEB)
    dsl4jsb_Def_memory(SEB)

    ! Pointers to variables in memory
    dsl4jsb_Real2D_onChunk :: t
    dsl4jsb_Real2D_onChunk :: t_old
    dsl4jsb_Real2D_onChunk :: t_unfilt
    dsl4jsb_Real2D_onChunk :: t_unfilt_old

    ! Locally allocated variables
    !
    INTEGER :: iblk, ics, ice, nc
    REAL(wp) :: eps

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':update_asselin_land'

    iblk = options%iblk
    ics  = options%ics
    ice  = options%ice
    nc   = options%nc

    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    dsl4jsb_Get_config(SEB)

    ! Get reference to variables for current block
    !
    dsl4jsb_Get_memory(SEB)

    dsl4jsb_Get_var2D_onChunk(SEB,      t)              ! OUT
    dsl4jsb_Get_var2D_onChunk(SEB,      t_old)          ! OUT
    dsl4jsb_Get_var2D_onChunk(SEB,      t_unfilt)       ! OUT
    dsl4jsb_Get_var2D_onChunk(SEB,      t_unfilt_old)   ! OUT

    ! Asselin time filter, if applicable
    eps = get_asselin_coef()
    IF (eps > 0._wp) THEN
      t(:) = t_unfilt_old(:) + eps * (t_old(:) - 2._wp * t_unfilt_old(:) + t_unfilt(:))
    ELSE
      t(:) = t_unfilt(:)
    END IF

    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE update_asselin_land

#endif
END MODULE mo_seb_land
