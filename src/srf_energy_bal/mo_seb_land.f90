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

  USE mo_kind,      ONLY: wp
  USE mo_exception, ONLY: message

  USE mo_jsb_tile_class,     ONLY: t_jsb_tile_abstract
  USE mo_jsb_config_class,   ONLY: t_jsb_config, t_jsb_config_p
  USE mo_jsb_task_class,     ONLY: t_jsb_task_options

  dsl4jsb_Use_processes SEB, RAD, ATM2LAND
  dsl4jsb_Use_config(SEB)

  USE mo_atmland_memory_class, ONLY: t_atm2land_memory
!!$  dsl4jsb_Use_memory(ATM2LAND)
  dsl4jsb_Use_memory(SEB)
  dsl4jsb_Use_memory(RAD)

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: update_surface_energy_land

  CHARACTER(len=*), PARAMETER :: modname = 'mo_seb_land'

CONTAINS

  SUBROUTINE update_surface_energy_land(tile, config, options)

    USE mo_phy_schemes,            ONLY: qsat, q_effective, surface_dry_static_energy, heat_transfer_coef
    USE mo_echam_vdiff_params,     ONLY: tpfac2, tpfac3
    USE mo_jsb_physical_constants, ONLY: tmelt

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    CLASS(t_jsb_config_p),      INTENT(in)    :: config(:)
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    ! Local variables
    !
    dsl4jsb_Def_config(SRF)
    dsl4jsb_Def_memory(SRF)
    dsl4jsb_Def_memory(SOIL)
    dsl4jsb_Def_memory(ATM2LAND)

    LOGICAL                       :: lstart

    ! Pointers to variables in memory
    dsl4jsb_Real2D_onChunk :: t_srf
    dsl4jsb_Real2D_onChunk :: t_srf_old
    dsl4jsb_Real2D_onChunk :: t_srf_unfilt
    dsl4jsb_Real2D_onChunk :: t_srf_unfilt_old
    dsl4jsb_Real2D_onChunk :: t_eff_srf
    dsl4jsb_Real2D_onChunk :: qsat_srf
!!$   dsl4jsb_Real2D_onChunk :: qsat_srf_old
    dsl4jsb_Real2D_onChunk :: s
!!$   dsl4jsb_Real2D_onChunk :: s_srf_old
    dsl4jsb_Real2D_onChunk :: press_srf
    dsl4jsb_Real2D_onChunk :: q_air
    dsl4jsb_Real2D_onChunk :: fact_qsat_srf
    dsl4jsb_Real2D_onChunk :: fact_q_air
    dsl4jsb_Real2D_onChunk :: drag_srf
    dsl4jsb_Real2D_onChunk :: rad_srf_net
    dsl4jsb_Real2D_onChunk :: grnd_hflx
    dsl4jsb_Real2D_onChunk :: hcap_grnd
    dsl4jsb_Real2D_onChunk :: wsn_srf
    dsl4jsb_Real2D_onChunk :: sfract_srf
    dsl4jsb_Real2D_onChunk :: t_acoef
    dsl4jsb_Real2D_onChunk :: t_bcoef
    dsl4jsb_Real2D_onChunk :: q_acoef
    dsl4jsb_Real2D_onChunk :: q_bcoef
    dsl4jsb_Real2D_onChunk :: glacier

    ! Locally allocated vectors
    !
    REAL(wp) ::             &
      qsat_srf_old(options%nc),     &
      s_old(options%nc),        &
      dQdT(options%nc),             & !< Sensitivity of saturated surface specific humidity to temperature
      t2s_conv(options%nc),         & !< Conversion factor from temperature to dry static energy (C_pd * (1+(delta-1)*q_v))
      heat_tcoef(options%nc)          !< Heat transfer coefficient (rho*C_h*|v|)

    INTEGER :: iblk, ics, ice, nc
    REAL(wp) :: steplen

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':update_surface_energy_land'

    iblk = options%iblk
    ics  = options%ics
    ice  = options%ice
    nc   = options%nc
    steplen = options%steplen

    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    dsl4jsb_Get_config(SRF)

    ! Get reference to variables for current block
    !
    dsl4jsb_Get_memory(SRF)
    dsl4jsb_Get_memory(SOIL)
    dsl4jsb_Get_memory(ATM2LAND)

    dsl4jsb_Get_var2D_onChunk(ATM2LAND, t_acoef)       ! IN
    dsl4jsb_Get_var2D_onChunk(ATM2LAND, t_bcoef)       ! IN
    dsl4jsb_Get_var2D_onChunk(ATM2LAND, q_acoef)       ! IN
    dsl4jsb_Get_var2D_onChunk(ATM2LAND, q_bcoef)       ! IN
    dsl4jsb_Get_var2D_onChunk(ATM2LAND, press_srf)     ! IN
    dsl4jsb_Get_var2D_onChunk(ATM2LAND, q_air)         ! IN
    dsl4jsb_Get_var2D_onChunk(ATM2LAND, drag_srf)      ! IN

    dsl4jsb_Get_var2D_onChunk(SRF,      t)             ! OUT
    dsl4jsb_Get_var2D_onChunk(SRF,      t_old)         ! OUT
    dsl4jsb_Get_var2D_onChunk(SRF,      t_srf_unfilt)  ! OUT
    dsl4jsb_Get_var2D_onChunk(SRF,      t_srf_unfilt_old) ! OUT
    dsl4jsb_Get_var2D_onChunk(SRF,      t_eff_srf)     ! OUT
    dsl4jsb_Get_var2D_onChunk(SRF,      qsat_srf)      ! OUT
    dsl4jsb_Get_var2D_onChunk(SRF,      s)             ! OUT
    dsl4jsb_Get_var2D_onChunk(SRF,      fact_q_air)    ! OUT
    dsl4jsb_Get_var2D_onChunk(SRF,      fact_qsat_srf) ! OUT
    dsl4jsb_Get_var2D_onChunk(RAD,      rad_srf_net)   ! OUT
    dsl4jsb_Get_var2D_onChunk(SRF,      wsn_srf)       ! OUT
    dsl4jsb_Get_var2D_onChunk(SRF,      sfract_srf)    ! OUT
    dsl4jsb_Get_var2D_onChunk(SRF,      glacier)       ! OUT

    dsl4jsb_Get_var2D_onChunk( SOIL,    grnd_hflx)     ! OUT
    dsl4jsb_Get_var2D_onChunk( SOIL,    hcap_grnd)     ! OUT

    ! Save old surface temperature and saturation spec. humidity
!!$    t_srf_unfilt_old(:) = t_srf_unfilt(:)
!!$    t_srf_old   (:) = t_srf(:)
    ! Temporary fix for temperatures over ocean getting too warm or too cold, leading to crash in convect_tables
    t_srf_unfilt_old(:) = MAX(MIN(t_srf_unfilt(:), 380._wp), 100._wp)
    t_old(:)        = MAX(MIN(t(:),        380._wp), 100._wp)
    qsat_srf_old(:) = qsat(t_old(:), press_srf(:))

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
    IF (ALL(qsat_srf(:) < 0._wp)) THEN          ! Start of experiment
      lstart = .TRUE.
      s(:) = s_old(:)
    ELSE
      lstart = .FALSE.
      CALL surface_temp_implicit( nc, steplen, &
        t2s_conv(:), &
        t_acoef(:), t_bcoef(:), q_acoef(:), q_bcoef(:), &
        s_old(:), qsat_srf_old(:), dQdT(:), &
        rad_srf_net(:), grnd_hflx(:), heat_tcoef(:), &
        fact_q_air(:), fact_qsat_srf(:), &
        sfract_srf(:), hcap_grnd(:), &
        & s(:))
    END IF

    ! New unfiltered surface temperature (\tilde(X)^(t+1)) (see below Eq. 3.3.1.5 in ECHAM3 manual, alpha = 1/tpfac2)
    t_srf_unfilt(:) = ( tpfac2 * s(:) + tpfac3 * s_old(:) ) / t2s_conv(:)

    ! Correct new dry static energy for snowmelt
    WHERE (wsn_srf(:) > dsl4jsb_Config(SRF)%wsn_critical .OR. glacier(:) > 0.5_wp)
      s(:) = MIN(s(:) / t2s_conv(:), tmelt) * t2s_conv(:)
    END WHERE

    ! Compute new surface qsat (after correction of dry static energy for snow melt!)
    qsat_srf(:) = 0._wp
    IF (lstart) THEN
      qsat_srf(:) = qsat_srf_old(:)
    ELSE
      qsat_srf(:) = qsat_srf_old(:) + dQdT(:) * (s(:) - s_old(:)) / t2s_conv(:)
    END IF

    ! 'Effective temperature' used in radheat of the atmosphere model (see Eq. 6.3 in ECHAM5 manual)
    t_eff_srf(:) = (t_old(:)**3 * (4._wp * s(:)/t2s_conv(:) - 3._wp * t_old(:))) ** 0.25_wp

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

#endif
END MODULE mo_seb_land
