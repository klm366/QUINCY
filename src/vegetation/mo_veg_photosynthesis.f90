!>
!! @brief Contains the routines for photosynthesis related processes 
!!
!! @author
!!  Sönke Zaehle, MPI-BGC Jena 
!!
!! @par Revision History
!! First version by Sönke Zaehle (2015-12-04)
!!
MODULE mo_veg_photosynthesis

  USE mo_kind, ONLY: wp
  USE mo_base_memory,    ONLY: t_slm_tile
  USE mo_base_constants, ONLY: ncanop, dt
  USE mo_veg_memory
  USE mo_base_memory,    ONLY: t_slm_soil
  USE mo_atmland_memory, ONLY: atm2land, land2atm
  USE mo_atmland_forcing
  USE mo_veg_parameters
  USE mo_veg_radiation,  ONLY: calc_canopy_layer_depth

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: update_canopy_fluxes,discrimination_ps

CONTAINS

  !-----------------------------------------------------------------------------------------------------
  !> update canopy fluxes, includes calculation of photosynthesis
  !! 
  !! integrates: \n
  !! CALL calc_canopy_n_distribution \n
  !! CALL calc_photosynthesis
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE update_canopy_fluxes(tile, veg_param)

   USE mo_shared_constants
   USE mo_shared_isotope

   IMPLICIT NONE
   ! 0.1 InOut
   TYPE(t_slm_veg_param), INTENT(in)        :: veg_param      !< list of vegetation parameters
   TYPE(t_slm_tile),INTENT(inout), TARGET   :: tile           !< one tile
   ! 0.2 Local
   TYPE(t_slm_soil),POINTER         :: soil                   ! soil variables per tile
   TYPE(t_slm_plant),POINTER        :: plant                  ! plant variables per tile
   REAL(wp),DIMENSION(ncanop)       :: agl,anl,rl,gsl,cil
   REAL(wp),DIMENSION(ncanop)       :: canopy_layer_depth_act
   INTEGER                          :: pftid                  ! plant functional type ID
   INTEGER                          :: icanop                 ! for looping across ncanopy layers
   REAL(wp)                         :: hlp1                   ! dummy variable 

   !  bi-directional pointer
   plant => tile%plant
   soil  => tile%soil

   ! initialise output
   agl(:) = 0.0_wp 
   anl(:) = 0.0_wp 
   rl(:)  = 0.0_wp
   gsl(:) = 0.0_wp
   cil(:) = 0.0_wp     

   pftid = plant%pftid

   ! calculate the moisture scaling factor for photosynthesis
   plant%beta_air = calc_atm_moisture_stress(atm2land%t_air,atm2land%q_air,atm2land%press_srf) 
   plant%beta_soil_ps = calc_soil_moisture_stress(soil%theta, & 
                       veg_param%lambda_soil_ps(pftid),veg_param%weibull_k_soil_ps(pftid))
   plant%beta_soil_gs = calc_soil_moisture_stress(soil%theta, &
                       veg_param%lambda_soil_gs(pftid),veg_param%weibull_k_soil_gs(pftid))

   ! determine actual canopy layers
   canopy_layer_depth_act = calc_canopy_layer_depth(plant%lai)

   ! update N distribution in the canopy
   CALL calc_canopy_n_distribution(plant%pool%leaf%nitrogen,plant%lai, & 
                                   canopy_layer_depth_act, &
                                   veg_param%ps_pathway(pftid), & 
                                   veg_param%nf(pftid), & 
                                   veg_param%k0_fnl_struc(pftid), & 
                                   plant%canopy_layer(:)%nl_tot, & 
                                   plant%canopy_layer(:)%fnl_chl, & 
                                   plant%canopy_layer(:)%fnl_et, & 
                                   plant%canopy_layer(:)%fnl_rub, & 
                                   plant%canopy_layer(:)%fnl_pepc)
   ! calculate photosynthesis per layer
   DO icanop = 1,ncanop
      IF(canopy_layer_depth_act(icanop).GT.0.0_wp)THEN

        CALL calc_photosynthesis(atm2land%t_air,atm2land%press_srf,atm2land%co2_mixing_ratio,atm2land%ga, & 
                                 plant%canopy_layer(icanop)%par_sunlit, &
                                 plant%canopy_layer(icanop)%par_shaded, & 
                                 plant%canopy_layer(icanop)%fleaf_sunlit, &
                                 plant%beta_air,plant%beta_soil_ps,plant%beta_soil_gs, &
                                 plant%canopy_layer(icanop)%fnl_chl, & 
                                 plant%canopy_layer(icanop)%fnl_et, & 
                                 plant%canopy_layer(icanop)%fnl_rub, & 
                                 plant%canopy_layer(icanop)%fnl_pepc, & 
                                 plant%canopy_layer(icanop)%nl_tot, & 
                                 veg_param%ps_pathway(pftid), &
                                 veg_param%tjmax_opt(pftid), &
                                 veg_param%tjmax_omega(pftid),veg_param%g0(pftid),veg_param%g1(pftid), &
                                 agl(icanop),anl(icanop),rl(icanop),gsl(icanop),cil(icanop))

     ENDIF
  ENDDO

  !> integrate to canopy level values
  !! 
  !! fluxes were set to zero before by set_zero_fluxes_of_plant_pools in the driver
!  plant%flux%gpp = 0.0_wp
!  plant%flux%an = 0.0_wp
!  plant%flux%maint_respiration_leaf = 0.0_wp
!  plant%flux%ci = 0.0_wp
  plant%gc = 0.0_wp 
  DO icanop = 1,ncanop
     plant%flux%gpp = plant%flux%gpp + agl(icanop) * canopy_layer_depth_act(icanop) 
     IF(agl(icanop) > eps8) THEN
        hlp1 = discrimination_ps(cil(icanop),atm2land%co2_mixing_ratio,'C13',veg_param%ps_pathway(pftid))
        plant%flux%gpp_C13 = plant%flux%gpp_C13 + agl(icanop) * canopy_layer_depth_act(icanop) * & 
                             calc_fractionation(atm2land%co2_mixing_ratio,atm2land%co2_mixing_ratio_C13,hlp1)
        hlp1 = discrimination_ps(cil(icanop),atm2land%co2_mixing_ratio,'C14',veg_param%ps_pathway(pftid))
        plant%flux%gpp_C14 = plant%flux%gpp_C14 + agl(icanop) * canopy_layer_depth_act(icanop) * & 
                             calc_fractionation(atm2land%co2_mixing_ratio,atm2land%co2_mixing_ratio_C14, & 
                             hlp1,scale_correction_C14)
     ENDIF
     plant%flux%an = plant%flux%an + anl(icanop) * canopy_layer_depth_act(icanop) 
     plant%flux%maint_respiration_leaf = plant%flux%maint_respiration_leaf + & 
            rl(icanop) * canopy_layer_depth_act(icanop) 
     plant%gc = plant%gc + gsl(icanop) * canopy_layer_depth_act(icanop)
     plant%flux%ci = plant%flux%ci + cil(icanop) * canopy_layer_depth_act(icanop)
  ENDDO

  plant%flux%ci = plant%flux%ci / sum(canopy_layer_depth_act)

  END SUBROUTINE update_canopy_fluxes

  !-----------------------------------------------------------------------------------------------------
  !> Dummy routine to produce the within canopy distribution of N of OCN
  !! 
  !! @todo replace by Silvia's allocation scheme and feed back information on foliar C:N to growth routine
  !! @par note 
  !!  it includes "nf" as input, i.e. Andrew's tuning constant
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE calc_canopy_n_distribution(totalN,lai,canopy_layer_depth_act,ps_pathway,nf,k0_fnl_struc,& 
                                      nl_tot,fnl_chl,fnl_et,fnl_rub,fnl_pepc)

  USE mo_veg_parameters, ONLY: kn,k0_fnl_chl_C3,k0_fnl_chl_C4,k1_fnl_chl_C3,k1_fnl_chl_C4,kfnl_chl, & 
                               jmax2n,vcmax2n,jmax2vcmax_C3,jmax2vcmax_C4,k1_fnl_struc,k0_fnl_pepc,ic4phot
 
  IMPLICIT NONE
  ! 0.1 InOut
  REAL(wp), INTENT(in)    :: totalN           !< total foliar N (mol / m2)
  REAL(wp), INTENT(in)    :: lai              !< total foliar area (m2 / m2)
  REAL(wp), DIMENSION(ncanop), INTENT(in) :: canopy_layer_depth_act !< depth of canopy layers (in units LAI)
  INTEGER, INTENT(in)     :: ps_pathway       !< photosynthetic pathway (C3|C4) 
  REAL(wp), INTENT(in)    :: nf               !< Andrew's tuning constant 
  REAL(wp), INTENT(in)    :: k0_fnl_struc     !< maximum fraction of N not used for photosynthesis 
  REAL(wp),DIMENSION(ncanop),INTENT(inout) :: nl_tot           !< total N content of the layer
  REAL(wp),DIMENSION(ncanop),INTENT(inout) :: fnl_chl          !< fraction of N in chlorophyll content per layer
  REAL(wp),DIMENSION(ncanop),INTENT(inout) :: fnl_et           !< fraction of N used for electron-transport per layer
  REAL(wp),DIMENSION(ncanop),INTENT(inout) :: fnl_rub          !< fraction of N used for Rubisco per layer
  REAL(wp),DIMENSION(ncanop),INTENT(inout) :: fnl_pepc         !< fraction of N used for PEPC per layer
  ! 0.2 Local
  REAL(wp) :: n0                              ! foliar N content of top-most layer
  REAL(wp) :: fnl_oth                         ! fraction of N not associated with PS
  REAL(wp) :: k0_fnl_chl,k1_fnl_chl           ! factors in chlorophyll profile calculation 
  REAL(wp) :: jmax2vcmax                      ! jmax to vcmax ratio
  REAL(wp) :: lc                              ! cummulative LAI
  INTEGER :: icanop

  ! calculate top of the canopy N 
  IF(lai > canopy_layer_depth_act(1)) THEN
     n0 = totalN * 1000._wp * kn / (1._wp - exp(-kn * lai))
  ELSE 
     n0 = totalN * 1000._wp
  ENDIF

  lc = canopy_layer_depth_act(1) / 2._wp
  IF(ps_pathway == ic4phot) THEN
      jmax2vcmax = jmax2vcmax_C4
      k0_fnl_chl = k0_fnl_chl_C4
      k1_fnl_chl = k1_fnl_chl_C4
      fnl_pepc = k0_fnl_pepc
  ELSE 
      jmax2vcmax = jmax2vcmax_C3
      k0_fnl_chl = k0_fnl_chl_C3
      k1_fnl_chl = k1_fnl_chl_C3
      fnl_pepc = 0.0_wp
  ENDIF

  DO icanop = 1,ncanop
     IF(canopy_layer_depth_act(icanop).GT.eps8) THEN
        nl_tot(icanop) = n0 * exp(-kn * lc) 
        fnl_oth = k0_fnl_struc - k1_fnl_struc * nl_tot(icanop)
        fnl_chl(icanop) = (k0_fnl_chl - k1_fnl_chl * exp(-kfnl_chl * lc)) / chl2n
        ! standard OCN approach: fraction in indenpent of depth
        !fnl_et(icanop) = 0.273_wp * nf     ! OBS: currently only for code consistency 
        !fnl_rub(icanop) = 0.128_wp * nf    ! OBS: currently only for code consistency
        ! fraction of electron-transport is dependent on N in electron transport versus N in light harvesting
        fnl_et(icanop) = (1._wp - (fnl_oth + fnl_chl(icanop) + fnl_pepc(icanop) )) / (1._wp + jmax2n/(jmax2vcmax * vcmax2n))
        ! fraction of Rubisco is determined as the remainder
        fnl_rub(icanop) = 1._wp - fnl_oth - fnl_chl(icanop) - fnl_et(icanop) - fnl_pepc(icanop)
        lc = lc + canopy_layer_depth_act(icanop) / 2._wp & 
             + canopy_layer_depth_act(icanop+1) / 2._wp
    ENDIF
  ENDDO
  END SUBROUTINE calc_canopy_n_distribution

  !-----------------------------------------------------------------------------------------------------
  !> Calculation of net photosynthesis and canopy conductance
  !! 
  !! according to Friend & Kiang 2005, J Climate \n
  !! accounts for temperature acclimation of JMAX according to Friend, 2010 J Exp Bot
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE calc_photosynthesis(temp,press_srf,co2_mixing_ratio,ga,par_sunlit,par_shaded,fleaf_sunlit,beta_air,beta_soil_ps, &
                          beta_soil_gs,fnl_chl,fnl_et,fnl_rub,fnl_pepc,nfmol,ps_pathway,tjmax_opt,tjmax_omega,g0,g1, &
                          agl,anl,rl,gsl,cil)  
    
  USE mo_shared_constants,     ONLY: r_gas
  USE mo_atmland_constants,  ONLY: Tzero
  USE mo_veg_parameters
  USE mo_veg_respiration,    ONLY: temperature_response_respiration

  IMPLICIT NONE

  ! 0.1 InOut
  REAL(wp),INTENT(IN) :: temp                    !< temperature used for photosynthesis calculation (K) 
  REAL(wp),INTENT(IN) :: press_srf               !< air pressure (Pa) 
  REAL(wp),INTENT(IN) :: ga                      !< aerodynamic conductance (m s-1) 
  REAL(wp),INTENT(IN) :: par_sunlit              !< PAR of sunlit foliage par per canopy layer (µmol m-2 s-1)
  REAL(wp),INTENT(IN) :: par_shaded              !< PAR of shaded foliage par per canopy layer (µmol m-2 s-1)
  REAL(wp),INTENT(IN) :: fleaf_sunlit            !< fraction of canopy layer that is sunlit 
  REAL(wp),INTENT(IN) :: co2_mixing_ratio        !< CO2 concentration of the free air (ppm)
  REAL(wp),INTENT(IN) :: beta_air                !< scaling factor to account for air humidity effects on conductance
  REAL(wp),INTENT(IN) :: beta_soil_ps            !< scaling factor to account for soil moisture constraints on conductance
  REAL(wp),INTENT(IN) :: beta_soil_gs            !< scaling factor to account for soil moisture constraints on conductance
  REAL(wp),INTENT(IN) :: fnl_et                  !< fraction of foliar N allocated to electron transport 
  REAL(wp),INTENT(IN) :: fnl_rub                 !< fraction of foliar N allocated to Rubisco
  REAL(wp),INTENT(IN) :: fnl_pepc                !< fraction of foliar N allocated to PEPC complex (C4 only)
  REAL(wp),INTENT(IN) :: fnl_chl                 !< fraction of foliar N allocated to Chlorophyll
  REAL(wp),INTENT(IN) :: nfmol                   !< total foliar N content
  INTEGER,INTENT(IN)  :: ps_pathway              !< photosynthetic pathway 
  REAL(wp),INTENT(IN) :: tjmax_opt               !< temperature optimum of electron transport
  REAL(wp),INTENT(IN) :: tjmax_omega             !< shape-parameter of temperature sensitivity of electron transport
  REAL(wp),INTENT(IN) :: g0                      !< offset of the stomatal conductance ~ An relationship 
  REAL(wp),INTENT(IN) :: g1                      !< slope of the stomatal conductance ~ An relationship
  REAL(wp),INTENT(OUT) :: agl                    !< gross photsynthesis of layer (µmol CO2 / m2 /s) 
  REAL(wp),INTENT(OUT) :: anl                    !< net photsynthesis of layer (µmol CO2 / m2 /s)
  REAL(wp),INTENT(OUT) :: rl                     !< net canopy layer respiration (µmol CO2 / m2 /s)
  REAL(wp),INTENT(OUT) :: gsl                    !< canopy conductance of layer (m / s) 
  REAL(wp),INTENT(OUT) :: cil                    !< leaf internal CO2 concentration of layer (ppm)
  ! 0.2 Local
  REAL(wp) :: kc,ko,km                ! Michaelis-Menten Kinetic parameters for Rubisco (Pa)
  REAL(wp) :: pcp                     ! Photosynthetic compensation point (Pa)
  REAL(wp) :: f_resp                  ! N-specific layer respiration rate (µmol CO2 / mmol N /s) 
  REAL(wp) :: nsat,nlim               ! foliar N that is light saturated / limited (mmol N)
  REAL(wp) :: n1,n2,n3                ! N specific rates of electron-transport/Rubisco/PEPc limited PS
  REAL(wp) :: m1,m2,m3                ! ci sensitivity of electron-transport/Rubisco/PEPc limited PS
  REAL(wp) :: msat                    ! N-specific light saturated carboxylation (umol[CO2]/mmol[N]/s)
  REAL(wp) :: nco                     ! co-limitation factor 
  REAL(wp) :: agsh,agsl               ! gross photosynthesis on sunlit/shaded leaves (µmol / m2 / s) 
  REAL(wp) :: fabs,fabsb              ! fractional light absorption on chloroplast
  REAL(wp) :: ppm2Pa                  ! factor to convert gas concentrations from ppm to Pa
  REAL(wp) :: rtck                    ! dry static energy
  REAL(wp) :: anl_old,danl            ! dummy variables needed for iteration loop of CI
  REAL(wp) :: agsf1,agsf2             ! dummy variables for PS calculation
  INTEGER :: it                       ! Iteration loop counter
  LOGICAL :: last_call                ! boolean to stop iteration at low conductance

  !> static energy (R_gas: J mol-1 K-1 -> rtck: J mol-1)
  rtck = r_gas * temp 
  !> conversion of ppm to Pa
  ppm2Pa = press_srf * 1.e-06
  
  !>  Net canopy layer respiration (umol[CO2]/m2/s)
  rl = fmaint_rate_base/1000._wp * beta_soil_ps * temperature_response_respiration(temp) * nfmol

  !> @par IF: initialise output for no-light condition
  IF(par_sunlit + par_shaded <= eps8) THEN
    agl = 0.0_wp 
    anl = agl - rl
    gsl = MAX(g0,g0 + g1 * anl * beta_air * beta_soil_gs / co2_mixing_ratio) * rtck / press_srf 
    cil = co2_mixing_ratio * ppm2Pa
  !> @par ELSE: calculate photosynthesis only if it's day
  ELSE
     !> Michaelis-Menten constant of Rubisco for CO2 (Pa)
     kc = exp(E0kc - E1kc/rtck) * ppm2Pa 
     !> Michaelis-Menten constant of Rubisco for O2 (kPa)
     ko = exp(E0ko - E1ko/rtck) * ppm2Pa 
     !> Combined Rubisco Km for CO2 (Pa); kc (1 + pO2/ko) is part of equation 5 in Kull and Kruijt 1998
     km = kc * (1.0_wp + pO2/ko)
     !> Photorespiratory compensation point (Pa[CO2])
     pcp = exp(E0pcp - E1pcp/rtck) * ppm2Pa 

     !> N specific rate of electron transport rate
     n1 = jmax2n * fnl_et * exp(-(((temp - Tzero - tjmax_opt)/tjmax_omega)**2)) * beta_soil_ps

     !> N specific rate of Rubisco turnover 
     n2 = vcmax2n * fnl_rub * exp(E0v - E1v/rtck) * beta_soil_ps

     !> Initial guess of leaf internal CO2 concentration
     SELECT CASE (ps_pathway)
     CASE (ic3phot)
        cil = co2_mixing_ratio * CiCa_default_C3 * ppm2Pa
     CASE (ic4phot)
        cil = co2_mixing_ratio * CiCa_default_C4 * ppm2Pa
     ENDSELECT

     !> iteration to find correct combination of cil, Anl and gsl
     anl_old = 1000.0_wp
     danl = 1.0_wp
     it = 1
     last_call = .FALSE.
     !> @par DO WHILE (it < ps_it_max .AND. danl > eps4 )
     DO WHILE (it < ps_it_max .AND. danl > eps4 ) 

        SELECT CASE (ps_pathway)
        CASE (ic3phot)
           !> ci sensitivity of Electron-transport limited photosynthesis
           m1 = cil/(cil + 2.0_wp * pcp)

           !> ci sensitivity of Rubisco-limited photosynthesis
           m2 = cil/(cil + km)
 
           !> N-specific light saturated carboxylation (umol[CO2]/mmol[N]/s)
           msat = MIN(m1 * n1, m2 * n2)

           !> nlim factor (part of equation 12 (Kull and Kruijt 1998))
           nco = msat/(alpha_i_C3 * ka * m1)

           !> A couple of factors to calculate photosynthesis
           !!  parts of equations 14-16 of Kull and Kruijt 1998
           agsf1 = MAX(1.0_wp - pcp / cil, 0.0_wp)
           agsf2 = agsf1 * alpha_i_C3 * m1
           agsf1 = agsf1 * msat

        CASE (ic4phot)
           !> N specific rate of PePc limited photosynthessis 
           n3 = pepc2n * fnl_pepc * 2**((temp - TZero - Tref_pepc )/ Tbase_pepc)

           !> ci sensitivity of Electron-transport limited photosynthesis
           m1 = ci_max/(ci_max + 2.0_wp * pcp)

           !> ci sensitivity of Rubisco-limited photosynthesis
           m2 = ci_max/(ci_max + km)
 
           !> ci sensitivity of PePc-limited photosynthesis
           m3 = cil / press_srf

           !> N-specific light saturated carboxylation (umol[CO2]/mmol[N]/s)
           msat = MIN(MIN(m1 * n1, m2 * n2), m3 * n3)
           !> nlim factor (part of equation 12 (Kull and Kruijt 1998))
           nco = msat/(alpha_i_C4 * ka * m1)

           !> A couple of factors to calculate photosynthesis
           !!  parts of equations 14-16 of Kull and Kruijt 1998
           agsf1 = MAX(1.0_wp - pcp / ci_max, 0.0_wp)
           agsf2 = agsf1 * alpha_i_C4 * m1
           agsf1 = agsf1 * msat
  
        END SELECT

        !> fabsb: fraction of possible PAR absorption in present light conditions
        fabsb = exp(-ka * chl2n * fnl_chl * nfmol)

        !> @par IF: Case of no direct radiation
        IF (fleaf_sunlit <= eps8) THEN
           !>  Theoretical cumulative N at co-limitation by Rubisco and
           !!  light absorption in sunlit foliage (mmol[N]/m2))
           nlim = -log(nco/(par_shaded * chl2n * fnl_chl + eps8))/(ka * chl2n * fnl_chl)
           nsat = MAX(MIN(nlim, nfmol), 0.0_wp)

           !> Relative PAR absorption in light-limited region (fraction)
           fabs = exp(-ka * chl2n * fnl_chl * nsat) - fabsb

           !> Rate of gross photosynthesis (umol/m2/s)
           !!  sum of photosynthesis in colimited, light saturated PS and PS in the 
           !!  light-limited region (assumed to scale linearily with light) 
           agl = agsf1 * nsat + agsf2 * par_shaded * fabs

        !> @par ELSE: Case of direct and diffuse radiation
        ELSE 
          !>  Theoretical cumulative N at co-limitation by Rubisco and
          !!  light absorption in sunlit foliage (mmol[N]/m2))
          nlim = -log(nco/(par_sunlit * chl2n * fnl_chl + eps8))/(ka * chl2n * fnl_chl)
          nsat = MAX(MIN(nlim, nfmol), 0.0_wp)

          !>  Relative PAR absorption in light-limited region of sunlit foliage (fraction)
          fabs = exp(-ka * chl2n * fnl_chl * nsat) - fabsb

          !> Rate of gross photosynthesis in sunlit foliage (umol/m2/s)
          agsl = agsf1 * nsat + agsf2 * par_sunlit * fabs

          !> Theoretical cumulative N at co-limitation by Rubisco and
          !!  light absorption in shaded foliage (mmol[N]/m2))
          nlim = -log(nco/(par_shaded * chl2n * fnl_chl + eps8))/(ka * chl2n * fnl_chl)
          nsat = MAX(MIN(nlim, nfmol), 0.0_wp)

          !> Relative PAR absorption in light-limited region of shaded foliage (fraction)
          fabs = exp(-ka * chl2n * fnl_chl * nsat) - fabsb

          !> Rate of gross photosynthesis in shaded foliage (umol/m2/s)
          agsh = agsf1 * nsat + agsf2 * par_shaded* fabs

          !> Mean gross photosynthesis over sunlit and shaded fractions (umol/m2/s)
          agl = fleaf_sunlit * agsl + (1.0_wp - fleaf_sunlit) * agsh

        END IF
        !> END IF

        !> Net layer photosynthesis (umol[CO2]/m2/s)
        anl = agl - rl

        !> calculate canopy conductance (m/s) following Ball & Berry
        gsl = MAX(g0, g0 + g1 * anl * beta_air * beta_soil_gs / ( co2_mixing_ratio )) * & 
              rtck / press_srf
        !> calculate canopy conductance (m/s) following Medlyn et al. 2011
        !! gsl = MAX(g0, 1.6_wp * ( 1._wp + g1 * beta_air * beta_soil_gs ) * anl / ( co2_mixing_ratio )) * rtck / press_srf

        !> New internal CO2 from balance of supply and demand
        !! ci = cs - An * r_tot (Friend 2005 eq. 7)
        !! Friend 2005 eq. 8: r_tot = ( Dh2oco2_air / gsl + Dh2oco2_turb / (speed * q_cdrag)) 
        !! Corrected to Bonan et al. 2008, Ecol. Clim: accounting for the molecular diffusion differences
        !! at the stomatal interface (1.6) and in the boundary layer, where turbulence plays a role (1.37)
        cil = co2_mixing_ratio * ppm2Pa - anl * rtck * 1e-06 * (Dwv2co2_air/gsl + Dwv2co2_turb/ga)
        cil = MAX (cil, pcp) 

        !> change in An since last iteration and save new anl 
        danl = abs(anl-anl_old)
        anl_old = anl

        !> break-off point for numerical stability at low conductance
        IF(last_call) THEN
           it = ps_it_max + 1
           cil = pcp 
        ENDIF
!        IF( gsl <= g0 * rtck / press_srf .OR. cil <= pcp ) THEN 
!           last_call = .TRUE.
!           cil = pcp 
!        ENDIF 
        it =  it + 1
     END DO
     !> END DO WHILE
  ENDIF
  !> END IF

  !> for output purposes, convert ci to ppm
  cil = cil / ppm2Pa

  END SUBROUTINE calc_photosynthesis

  !-----------------------------------------------------------------------------------------------------
  !> calculate atmospheric drought stress factor
  !!
  !! @par Literature:
  !! Ball & Berry: beta = q_air/q_sat \n
  !! NOT used, Medlyn et al. 2011: beta = MIN(5._wp,1._wp / SQRT(vpd + eps8))
  !-----------------------------------------------------------------------------------------------------
  FUNCTION calc_atm_moisture_stress(temp,q_air,press_srf) RESULT(beta)

     USE mo_atmland_constants,  ONLY: eps_vpd
     USE mo_atmland_util,       ONLY: calc_spec_humidity_sat

     IMPLICIT NONE
     ! 0.1 InOut
     REAL(wp),INTENT(in) :: temp               !< temperature in K
     REAL(wp),INTENT(in) :: q_air              !< specific air humidity (g g-1)
     REAL(wp),INTENT(in) :: press_srf          !< air pressure (Pa)
     REAL(wp) :: beta                          !< atmospheric moisture scalar for use in conductance calculation
     ! 0.2 Local
     REAL(wp) :: q_sat                         ! saturating specific humidity 
     REAL(wp) :: vpd                           ! vapour pressure deficit 

     ! saturating specific humidity [kg kg-1]
     q_sat = calc_spec_humidity_sat(temp,press_srf) 
     ! vapour pressure deficit [kPa] 
     vpd = MAX(q_sat - q_air,0.0_wp) * press_srf / eps_vpd / 1000._wp
     ! scalar for atmospheric humidity as relative humidity (fractions)
     beta = q_air/q_sat ! Ball & Berry
     !beta = MIN(5._wp,1._wp / SQRT(vpd + eps8)) ! Medlyn et al. 2011
  END FUNCTION calc_atm_moisture_stress

  !-----------------------------------------------------------------------------------------------------
  !> calculate soil moisture stress on photosynthesis or conductance given current soil moisture
  !! 
  !! @todo this should account for rooting and soil moisture profile in the long-term
  !-----------------------------------------------------------------------------------------------------
  FUNCTION calc_soil_moisture_stress(theta,lambda,weibull_k) RESULT(beta)

    IMPLICIT NONE
    ! 0.1 InOut
    REAL(wp),INTENT(in) :: theta       !< soil plant available water (%)
    REAL(wp),INTENT(in) :: lambda      !< slope parameter of Weibull function
    REAL(wp),INTENT(in) :: weibull_k   !< shape parameter of Weibull function
    REAL(wp) :: beta                   !< plant\%beta_soil_ps / plant\%beta_soil_gs

    !> 1.0 Moisture stress now as a Weibull function to avoid MIN/MAX statements
    !!
    beta = 1.0_wp - exp(-(lambda*theta/100._wp)**weibull_k)
  END FUNCTION calc_soil_moisture_stress

  !-----------------------------------------------------------------------------------------------------
  !> Calculate discrimination of C13 and C14 due to photosynthesis 
  !!
  !! Follows the orignial ideas of Farquhar, as described in 
  !! Drake (2014), Radiocarbon, 56, 29-38 
  !! 
  !! Input: leaf-level Ci and Ca, photosynthetic pathway (C3/C4) and name of the isotope
  !!
  !! Output: Discrimination (per mill) 
  !-----------------------------------------------------------------------------------------------------
  PURE ELEMENTAL FUNCTION discrimination_ps(ci,ca,isotope_name,ps_pathway) RESULT(discrimination)

    USE mo_veg_parameters
    IMPLICIT NONE
    ! 0.1 Inout 
    REAL(wp),INTENT(in) :: ci                     !< leaf internal CO2 concentration (ppm) 
    REAL(wp),INTENT(in) :: ca                     !< atmospheric CO2 concentration (ppm) 
    CHARACTER(3),INTENT(in) :: isotope_name       !< isotope C13 or C14
    INTEGER,INTENT(in) :: ps_pathway              !< photosynthetic pathway (C3 or C4)
    REAL(wp) :: discrimination                    !< isotopic discrimination (per mill)

    SELECT CASE(ps_pathway)
    CASE (ic3phot)
      SELECT CASE(isotope_name)
        CASE ('C13')
          discrimination = discr_ps_a_C13 + (discr_ps_b_C13 - discr_ps_a_C13) * ci/ca
        CASE ('C14')
          discrimination = discr_ps_a_C14 + (discr_ps_b_C14 - discr_ps_a_C14) * ci/ca
      END SELECT
    CASE (ic4phot)
      SELECT CASE(isotope_name)
        CASE ('C13')
          discrimination = discr_ps_a_C13 + (discr_ps_c_C13 + discr_ps_phi * discr_ps_b_C13 - discr_ps_a_C13) * ci/ca
        CASE ('C14')
          discrimination = discr_ps_a_C14 + (discr_ps_c_C14 + discr_ps_phi * discr_ps_b_C14 - discr_ps_a_C14) * ci/ca
      END SELECT
    END SELECT

  END FUNCTION discrimination_ps

END MODULE mo_veg_photosynthesis
