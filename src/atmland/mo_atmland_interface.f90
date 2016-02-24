!>
!! @brief Interface to run JSBACH for one time step
!!
!! @author
!!  Reiner Schnur, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Reiner Schnur (2013-04-04)
!!
MODULE mo_atmland_interface
#ifndef __NO_JSBACH__

#include "dsl_jsbach.inc"

  USE mo_exception,      ONLY: finish
  USE mo_kind,           ONLY: wp
  USE mo_jsb_tile_class, ONLY: t_jsb_tile

  dsl4jsb_Use_processes ATM2LAND, SRF, SOIL
  USE mo_atmland_memory_class, ONLY: t_atm2land_memory
!!$  dsl4jsb_Use_memory(ATM2LAND)
  dsl4jsb_Use_memory(SRF)
  dsl4jsb_Use_memory(SOIL)

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: update_atm2land, update_land2atm

  CHARACTER(len=*), PARAMETER :: modname = 'mo_atmland_interface'

CONTAINS

  SUBROUTINE update_atm2land(tile, iblk, ics, ice, &
    & t_air, q_air, rain, snow, wind_air, wind_10m, lw_srf_down, swvis_srf_down, swnir_srf_down, swpar_srf_down, &
    & press_srf, drag_srf, t_acoef, t_bcoef, q_acoef, q_bcoef,&
!!!!! de Vrese (01.12.2015) Richtmyer-Morton coeffcients for the flux aggregation - start
    &t_acoef_flux, t_bcoef_flux, q_acoef_flux, q_bcoef_flux,  &
!!!!! de Vrese (01.12.2015) Richtmyer-Morton coeffcients for the flux aggregation - end
!!!!! de Vrese (18.01.2016) Fields of variables on the lowest atmosperic model levels required for the Calculation of tile specific drag coefficients - star
    & geom_b,                           &! in, pgeom1 = geopotential above ground
    & qxm1_b,                           &! in, xlm1 + xim1 = total concentration of hydrometeors
    & aclc_b,                           &! in, cloud fraction
    & theta_b,                          &! in, potential temp.
    & thetav_b,                         &! in, virtual potential temp.
    & thetal_b,                         &! in, liquid water (?) pot. temp.
    & qsat_b,                           &! in, saturation specific humidity   
    & lh_b,                             &! in, latent heat  
    & fac_sfc,                          &! in factor: drag_srf = pfac_sfc * pcfh_tile
!!!!! de Vrese (18.01.2016) Fields of variables on the lowest atmosperic model levels required for the Calculation of tile specific drag coefficients - end
    &pch, cos_zenith_angle)


    TYPE(t_jsb_tile),   INTENT(inout) :: tile
    INTEGER,            INTENT(in)    :: iblk, ics, ice
    REAL(wp), OPTIONAL, INTENT(in)    :: &
      & t_air(:), q_air(:), rain(:), snow(:), wind_air(:), wind_10m(:), lw_srf_down(:), swvis_srf_down(:), swnir_srf_down(:), &
      & swpar_srf_down(:), press_srf(:), drag_srf(:), t_acoef(:), t_bcoef(:), q_acoef(:), q_bcoef(:), pch(:), cos_zenith_angle(:) &
!!!!! de Vrese (01.12.2015) Richtmyer-Morton coeffcients for the flux aggregation - start
      &, t_acoef_flux(:), t_bcoef_flux(:), q_acoef_flux(:), q_bcoef_flux(:),&
!!!!! de Vrese (01.12.2015) Richtmyer-Morton coeffcients for the flux aggregation - end
!!!!! de Vrese (18.01.2016) Fields of variables on the lowest atmosperic model levels required for the Calculation of tile specific drag coefficients - star
      & geom_b(:),                           &! in, pgeom1 = geopotential above ground
      & qxm1_b(:),                           &! in, xlm1 + xim1 = total concentration of hydrometeors
      & aclc_b(:),                           &! in, cloud fraction
      & theta_b(:),                          &! in, potential temp.
      & thetav_b(:),                         &! in, virtual potential temp.
      & thetal_b(:),                         &! in, liquid water (?) pot. temp.
      & qsat_b(:),                           &! in, saturation specific humidity   
      & lh_b(:),                             &! in, latent heat  
      & fac_sfc(:)                            ! in factor: drag_srf = pfac_sfc * pcfh_tile
!!!!! de Vrese (18.01.2016) Fields of variables on the lowest atmosperic model levels required for the Calculation of tile specific drag coefficients - end

    dsl4jsb_Def_memory(ATM2LAND)

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':update_atm2land'

    IF (ice-ics+1 /= SIZE(t_air,1)) CALL finish(TRIM(routine), 'Wrong dimensions')

    dsl4jsb_Get_memory(ATM2LAND)

    IF (PRESENT(t_air))            dsl4jsb_var2D_onChunk(ATM2LAND,t_air)            = t_air            (:)
    IF (PRESENT(q_air))            dsl4jsb_var2D_onChunk(ATM2LAND,q_air)            = q_air            (:)
    IF (PRESENT(rain))             dsl4jsb_var2D_onChunk(ATM2LAND,rain)             = rain             (:)
    IF (PRESENT(snow))             dsl4jsb_var2D_onChunk(ATM2LAND,snow)             = snow             (:)
    IF (PRESENT(wind_air))         dsl4jsb_var2D_onChunk(ATM2LAND,wind_air)         = wind_air         (:)
    IF (PRESENT(wind_10m))         dsl4jsb_var2D_onChunk(ATM2LAND,wind_10m)         = wind_10m         (:)
    IF (PRESENT(lw_srf_down))      dsl4jsb_var2D_onChunk(ATM2LAND,lw_srf_down)      = lw_srf_down      (:)
    IF (PRESENT(swvis_srf_down))   dsl4jsb_var2D_onChunk(ATM2LAND,swvis_srf_down)   = swvis_srf_down   (:)
    IF (PRESENT(swnir_srf_down))   dsl4jsb_var2D_onChunk(ATM2LAND,swnir_srf_down)   = swnir_srf_down   (:)
    IF (PRESENT(swpar_srf_down))   dsl4jsb_var2D_onChunk(ATM2LAND,swpar_srf_down)   = swpar_srf_down   (:)
    IF (PRESENT(press_srf))        dsl4jsb_var2D_onChunk(ATM2LAND,press_srf)        = press_srf        (:)
    IF (PRESENT(drag_srf))         dsl4jsb_var2D_onChunk(ATM2LAND,drag_srf)         = drag_srf         (:)
    IF (PRESENT(t_acoef))          dsl4jsb_var2D_onChunk(ATM2LAND,t_acoef)          = t_acoef          (:)
    IF (PRESENT(t_bcoef))          dsl4jsb_var2D_onChunk(ATM2LAND,t_bcoef)          = t_bcoef          (:)
    IF (PRESENT(q_acoef))          dsl4jsb_var2D_onChunk(ATM2LAND,q_acoef)          = q_acoef          (:)
    IF (PRESENT(q_bcoef))          dsl4jsb_var2D_onChunk(ATM2LAND,q_bcoef)          = q_bcoef          (:)
!!!!! de Vrese (01.12.2015) Richtmyer-Morton coeffcients for the flux aggregation - start
    IF (PRESENT(t_acoef_flux))     dsl4jsb_var2D_onChunk(ATM2LAND,t_acoef_flux)     = t_acoef_flux     (:)
    IF (PRESENT(t_bcoef_flux))     dsl4jsb_var2D_onChunk(ATM2LAND,t_bcoef_flux)     = t_bcoef_flux     (:)
    IF (PRESENT(q_acoef_flux))     dsl4jsb_var2D_onChunk(ATM2LAND,q_acoef_flux)     = q_acoef_flux     (:)
    IF (PRESENT(q_bcoef_flux))     dsl4jsb_var2D_onChunk(ATM2LAND,q_bcoef_flux)     = q_bcoef_flux     (:)
!!!!! de Vrese (01.12.2015) Richtmyer-Morton coeffcients for the flux aggregation - end
!!!!! de Vrese (18.01.2016) Fields of variables on the lowest atmosperic model levels required for the Calculation of tile specific drag coefficients - star
    IF (PRESENT(geom_b))           dsl4jsb_var2D_onChunk(ATM2LAND,geom_b)           = geom_b           (:)
    IF (PRESENT(qxm1_b))           dsl4jsb_var2D_onChunk(ATM2LAND,qxm1_b)           = qxm1_b           (:)
    IF (PRESENT(aclc_b))           dsl4jsb_var2D_onChunk(ATM2LAND,aclc_b)           = aclc_b           (:)   
    IF (PRESENT(theta_b))          dsl4jsb_var2D_onChunk(ATM2LAND,theta_b)          = theta_b          (:)
    IF (PRESENT(thetav_b))         dsl4jsb_var2D_onChunk(ATM2LAND,thetav_b)         = thetav_b         (:)
    IF (PRESENT(thetal_b))         dsl4jsb_var2D_onChunk(ATM2LAND,thetal_b)         = thetal_b         (:)
    IF (PRESENT(qsat_b))           dsl4jsb_var2D_onChunk(ATM2LAND,qsat_b)           = qsat_b           (:)
    IF (PRESENT(lh_b))             dsl4jsb_var2D_onChunk(ATM2LAND,lh_b)             = lh_b             (:)
    IF (PRESENT(fac_sfc))          dsl4jsb_var2D_onChunk(ATM2LAND,fac_sfc)          = fac_sfc          (:)
!!!!! de Vrese (18.01.2016) Fields of variables on the lowest atmosperic model levels required for the Calculation of tile specific drag coefficients - end
    IF (PRESENT(pch))              dsl4jsb_var2D_onChunk(ATM2LAND,pch)              = pch              (:)
    IF (PRESENT(cos_zenith_angle)) dsl4jsb_var2D_onChunk(ATM2LAND,cos_zenith_angle) = cos_zenith_angle (:)


  END SUBROUTINE update_atm2land

  SUBROUTINE update_land2atm(tile, iblk, ics, ice,                      &
    & t_srf, t_eff_srf, qsat_srf, s_srf,                                &
    & fact_q_air, fact_qsat_srf,                                        &
    & evapotrans, latent_hflx, sensible_hflx, grnd_hflx, grnd_hcap,     &
    & rough_h_srf, rough_m_srf, tte_corr,                               &
    & alb_vis_dir, alb_nir_dir, alb_vis_dif, alb_nir_dif,               &
!!!!! de Vrese (01.12.2015) Variable for testing - start
    & test_flux_1, &
    & test_flux_2, &
    & test_flux_3, &
    & test_flux_4, &
    & test_flux_5  &
!!!!! de Vrese (01.12.2015) Variable for testing - end
    &)

    TYPE(t_jsb_tile),   INTENT(in) :: tile
    INTEGER,            INTENT(in) :: iblk, ics, ice
    REAL(wp), OPTIONAL, INTENT(out) ::                                               &
      & t_srf(:), t_eff_srf(:), qsat_srf(:), s_srf(:),                               &
      & fact_q_air(:), fact_qsat_srf(:),                                             &
      & evapotrans(:), latent_hflx(:), sensible_hflx(:), grnd_hflx(:), grnd_hcap(:), &
      & rough_h_srf(:), rough_m_srf(:), tte_corr(:),                                 &
      & alb_vis_dir(:), alb_nir_dir(:), alb_vis_dif(:), alb_nir_dif(:),              &
!!!!! de Vrese (01.12.2015) Variable for testing - start
      & test_flux_1(:),                                                              &
      & test_flux_2(:),                                                              &
      & test_flux_3(:),                                                              &
      & test_flux_4(:),                                                              &
      & test_flux_5(:)
!!!!! de Vrese (01.12.2015) Variable for testing - end


    dsl4jsb_Def_memory(SRF)
    dsl4jsb_Def_memory(SOIL)

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':update_land2atm'

    IF (ice-ics+1 /= SIZE(t_srf,1)) CALL finish(TRIM(routine), 'Wrong dimensions')

    dsl4jsb_Get_memory(SRF)
    dsl4jsb_Get_memory(SOIL)

    IF (PRESENT(t_srf))           t_srf(:)            = dsl4jsb_var2D_onChunk(SRF,t_srf)
    IF (PRESENT(t_eff_srf))       t_eff_srf(:)        = dsl4jsb_var2D_onChunk(SRF,t_eff_srf)
    IF (PRESENT(qsat_srf))        qsat_srf(:)         = dsl4jsb_var2D_onChunk(SRF,qsat_srf)
    IF (PRESENT(s_srf))           s_srf(:)            = dsl4jsb_var2D_onChunk(SRF,s_srf)
    IF (PRESENT(fact_q_air))      fact_q_air(:)       = dsl4jsb_var2D_onChunk(SRF,fact_q_air)
    IF (PRESENT(fact_qsat_srf))   fact_qsat_srf(:)    = dsl4jsb_var2D_onChunk(SRF,fact_qsat_srf)
    IF (PRESENT(evapotrans))      evapotrans(:)       = dsl4jsb_var2D_onChunk(SRF,evapotrans)
    IF (PRESENT(latent_hflx))     latent_hflx(:)      = dsl4jsb_var2D_onChunk(SRF,latent_hflx)
    IF (PRESENT(sensible_hflx))   sensible_hflx(:)    = dsl4jsb_var2D_onChunk(SRF,sensible_hflx)
    IF (PRESENT(grnd_hflx))       grnd_hflx(:)        = dsl4jsb_var2D_onChunk(SOIL,grnd_hflx_old)
    IF (PRESENT(grnd_hcap))       grnd_hcap(:)        = dsl4jsb_var2D_onChunk(SOIL,hcap_grnd_old)
    IF (PRESENT(rough_h_srf))     rough_h_srf(:)      = dsl4jsb_var2D_onChunk(SRF,rough_h_srf)
    IF (PRESENT(rough_m_srf))     rough_m_srf(:)      = dsl4jsb_var2D_onChunk(SRF,rough_m_srf)
    IF (PRESENT(tte_corr))        tte_corr(:)         = dsl4jsb_var2D_onChunk(SRF,tte_corr)
    IF (PRESENT(alb_vis_dir))     alb_vis_dir(:)      = dsl4jsb_var2D_onChunk(SRF,alb_vis_srf)  ! Temporary
    IF (PRESENT(alb_vis_dif))     alb_vis_dif(:)      = dsl4jsb_var2D_onChunk(SRF,alb_vis_srf)  ! Temporary
    IF (PRESENT(alb_nir_dir))     alb_nir_dir(:)      = dsl4jsb_var2D_onChunk(SRF,alb_nir_srf)  ! Temporary
    IF (PRESENT(alb_nir_dif))     alb_nir_dif(:)      = dsl4jsb_var2D_onChunk(SRF,alb_nir_srf)  ! Temporary
!!!!! de Vrese (01.12.2015) Variable for testing - start
    IF (PRESENT(test_flux_1))     test_flux_1(:)      = dsl4jsb_var2D_onChunk(SRF,test_flux_1)
    IF (PRESENT(test_flux_2))     test_flux_2(:)      = dsl4jsb_var2D_onChunk(SRF,test_flux_2)
    IF (PRESENT(test_flux_3))     test_flux_3(:)      = dsl4jsb_var2D_onChunk(SRF,test_flux_3)
    IF (PRESENT(test_flux_4))     test_flux_4(:)      = dsl4jsb_var2D_onChunk(SRF,test_flux_4)
    IF (PRESENT(test_flux_5))     test_flux_5(:)      = dsl4jsb_var2D_onChunk(SRF,test_flux_5)
!!!!! de Vrese (01.12.2015) Variable for testing - end

  END SUBROUTINE update_land2atm

#endif
END MODULE mo_atmland_interface
