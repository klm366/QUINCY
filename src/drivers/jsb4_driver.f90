!>
!! @brief Contains the driver of the standalone jsbach.
!!
!! @author
!!  Julia Nabel, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Julia Nabel (2015-11-02)
!!
PROGRAM jsbalone_driver

!USE mo_jsbach_interface, ONLY: jsbach_interface

!! call copied from echam_dev:
!       CALL jsbach_interface ( jg, nblock, 1, kproma, pdtime, psteplen,                  & ! in
!         & t_air            = ptemp(1:kproma),                                           & ! in
!         & q_air            = pq(1:kproma),                                              & ! in
!         & rain             = prsfl(1:kproma) + prsfc(1:kproma),                         & ! in
!         & snow             = pssfl(1:kproma) + pssfc(1:kproma),                         & ! in
!         & wind_air         = SQRT(pu(1:kproma)**2 + pv(1:kproma)**2),                   & ! in
!         & wind_10m         = SQRT(pu(1:kproma)**2 + pv(1:kproma)**2),                   & ! in, temporary
!         & lw_srf_down      = zlw_down(1:kproma),                                        & ! in
!         & swvis_srf_down   = zswvis_down(1:kproma),                                     & ! in
!         & swnir_srf_down   = zswnir_down(1:kproma),                                     & ! in
!         & swpar_srf_down   = pswpar_down(1:kproma),                                     & ! in
!         & press_srf        = presi_old(1:kproma,klev+1),                                & ! in
!         & drag_srf         = pfac_sfc(1:kproma) * pcfh_tile(1:kproma,idx_lnd),          & ! in
!         & t_acoef          = zen_h(1:kproma, idx_lnd),                                  & ! in
!         & t_bcoef          = zfn_h(1:kproma, idx_lnd),                                  & ! in
!         & q_acoef          = zen_qv(1:kproma, idx_lnd),                                 & ! in
!         & q_bcoef          = zfn_qv(1:kproma, idx_lnd),                                 & ! in
!         & pch              = MERGE(pch_tile(1:kproma,idx_lnd),1._wp,lsm(1:kproma)>0.5_wp),  & ! in
!         & cos_zenith_angle = pcosmu0(1:kproma),                                         & ! in
!         & t_srf            = ztsfc_lnd(1:kproma),                     & ! out (T_s^(n+1))
!                                                                         !   surface temp (filtered, if Asselin)
!         & t_eff_srf        = ztsfc_lnd_eff(1:kproma),                 & ! out (T_s^eff)
!                                                                         !   surface temp (effective, for radheat)
!         & qsat_srf         = sat_surface_specific_humidity(1:kproma),                   & ! out
!         & s_srf            = dry_static_energy(1:kproma),             & ! out
!                                                                         !(s_s^star, for vertical diffusion scheme)
!         & fact_q_air       = pcair(1:kproma),                                           & ! out
!         & fact_qsat_srf    = pcsat(1:kproma),                                           & ! out
!         & evapotrans       = evapotranspiration(1:kproma),                              & ! out
!         & latent_hflx      = plhflx_tile(1:kproma, idx_lnd),                            & ! out
!         & sensible_hflx    = pshflx_tile(1:kproma, idx_lnd),                            & ! out
!         & grnd_hflx        = zgrnd_hflx(1:kproma, idx_lnd),                             & ! out
!         & grnd_hcap        = zgrnd_hcap(1:kproma, idx_lnd),                             & ! out
!         & rough_h_srf      = z0h_lnd(1:kproma),                                         & ! out
!         & rough_m_srf      = z0m_lnd(1:kproma),                                         & ! out
!         & tte_corr         = tte_corr(1:kproma),                                        & ! out
!         & alb_vis_dir      = albvisdir(1:kproma),                                       & ! out
!         & alb_nir_dir      = albnirdir(1:kproma),                                       & ! out
!         & alb_vis_dif      = albvisdif(1:kproma),                                       & ! out
!         & alb_nir_dif      = albnirdif(1:kproma)                                        & ! out
!         )

print *, 'dummy'


END PROGRAM jsbalone_driver
