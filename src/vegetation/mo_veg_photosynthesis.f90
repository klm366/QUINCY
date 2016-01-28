!!
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
 USE mo_base_constants, ONLY: ncanop,ntile,dt
 USE mo_base_memory
 USE mo_veg_parameters
 USE mo_veg_radiation, ONLY: calc_canopy_layer_depth

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: update_canopy_fluxes

CONTAINS

SUBROUTINE update_canopy_fluxes(tile,t_air,q_air,press_srf,co2_mixing_ratio,ga,veg_param)

   IMPLICIT NONE
   ! 0.1 Input
   TYPE(VegetationParameterClass),INTENT(in) :: veg_param  ! list of vegetation parameters
   REAL(wp),INTENT(in) :: t_air
   REAL(wp),INTENT(in) :: q_air
   REAL(wp),INTENT(in) :: press_srf
   REAL(wp),INTENT(in) :: co2_mixing_ratio
   REAL(wp),INTENT(in) :: ga
   ! 0.2 In-Output
   TYPE(TileClass),DIMENSION(ntile) :: tile  ! root tile
   ! 0.3 Output
   ! 0.4 Local
   REAL(wp),DIMENSION(ncanop) :: agl,anl,rl,gsl,cil
   REAL(wp),DIMENSION(ncanop) :: canopy_layer_depth_act
   REAL(wp) :: beta_air
   REAL(wp) :: beta_soil_ps
   REAL(wp) :: beta_soil_gs
   INTEGER :: pftid
   INTEGER :: t,icanop

   DO t = 1,ntile

     ! initialse output
     agl(:) = 0.0_wp 
     anl(:) = 0.0_wp 
     rl(:) = 0.0_wp
     gsl(:) = 0.0_wp
     cil(:) = 0.0_wp     

     pftid = tile(t)%plant%pftid

     ! calculate the moisture scaling factor for photosynthesis
     beta_air = calc_atm_moisture_stress(t_air,q_air,press_srf) 
     beta_soil_ps = calc_soil_moisture_stress(tile(t)%soil%theta, & 
                         veg_param%lambda_soil_ps(pftid),veg_param%weibull_k_soil_ps(pftid))
     beta_soil_gs = calc_soil_moisture_stress(tile(t)%soil%theta, &
                         veg_param%lambda_soil_gs(pftid),veg_param%weibull_k_soil_gs(pftid))

     ! determine actual canopy layers
     canopy_layer_depth_act = calc_canopy_layer_depth(tile(t)%plant%lai)

     ! update N distribution in the canopy
     CALL calc_canopy_n_distribution(tile(t)%plant%pool%leaf%nitrogen,tile(t)%plant%lai, & 
                                     canopy_layer_depth_act, &
                                     veg_param%nf(pftid), & 
                                     tile(t)%plant%pool%foliarN%nl_tot, & 
                                     tile(t)%plant%pool%foliarN%fnl_chl, & 
                                     tile(t)%plant%pool%foliarN%fnl_et, & 
                                     tile(t)%plant%pool%foliarN%fnl_rub)
     ! calculate photosynthesis per layer
     DO icanop = 1,ncanop
        IF(tile(t)%plant%lai.GT.0.0_wp)THEN

          CALL calc_photosynthesis(t_air,press_srf,co2_mixing_ratio,ga, & 
                                    tile(t)%plant%par_sunlit(icanop),tile(t)%plant%par_shaded(icanop),& 
                                    tile(t)%plant%fleaf_sunlit(icanop), &
                                    beta_air,beta_soil_ps,beta_soil_gs, &
                                    tile(t)%plant%pool%foliarN%fnl_chl(icanop),tile(t)%plant%pool%foliarN%fnl_et(icanop),& 
                                    tile(t)%plant%pool%foliarN%fnl_rub(icanop),tile(t)%plant%pool%foliarN%nl_tot(icanop),& 
                                    veg_param%tjmax_opt(pftid), &
                                    veg_param%tjmax_omega(pftid),veg_param%g0(pftid),veg_param%g1(pftid), &
                                    agl(icanop),anl(icanop),rl(icanop),gsl(icanop),cil(icanop))

        ENDIF
     ENDDO

     ! integrate to canopy level values
     tile(t)%plant%flux%gpp = 0.0_wp
     tile(t)%plant%flux%gc = 0.0_wp
     tile(t)%plant%flux%ci = 0.0_wp
     DO icanop = 1,ncanop
        tile(t)%plant%flux%gpp = tile(t)%plant%flux%gpp + anl(icanop) * canopy_layer_depth_act(icanop) &  
                * dt / 1000000._wp
        tile(t)%plant%flux%gc = tile(t)%plant%flux%gc + gsl(icanop) * canopy_layer_depth_act(icanop)
        tile(t)%plant%flux%ci = tile(t)%plant%flux%ci + cil(icanop) * canopy_layer_depth_act(icanop)
     ENDDO

     tile(t)%plant%flux%ci = tile(t)%plant%flux%ci / sum(canopy_layer_depth_act)

  ENDDO

  END SUBROUTINE update_canopy_fluxes

SUBROUTINE calc_canopy_n_distribution(totalN,lai,canopy_layer_depth_act,nf,nl_tot,fnl_chl,fnl_et,fnl_rub)
  !
  ! Dummy routine to produce the within canopy distribution of N of OCN
  ! NB: Should be replaced by Silvia's allocation scheme and feed back information on foliar C:N to 
  !     growth routine
  ! 
  USE mo_veg_parameters, ONLY: kn,k0_fnl_chl,k1_fnl_chl,kfnl_chl
 
  IMPLICIT NONE
  ! 0.1 Input
  REAL(wp), INTENT(in)    :: totalN           ! total foliar N (mol / m2)
  REAL(wp), INTENT(in)    :: lai              ! total foliar area (m2 / m2)
  REAL(wp), DIMENSION(ncanop), INTENT(in) :: canopy_layer_depth_act ! depth of canopy layers (in units LAI)
  REAL(wp), INTENT(in)    :: nf               ! Andrew's tuning constant 
  ! 0.2 In- and Output
  REAL(wp),DIMENSION(ncanop),INTENT(inout) :: nl_tot           ! total N content of the layer
  REAL(wp),DIMENSION(ncanop),INTENT(inout) :: fnl_chl          ! fraction of N in chlorophyll content per layer
  REAL(wp),DIMENSION(ncanop),INTENT(inout) :: fnl_et           ! fraction of N used for electron-transport per layer
  REAL(wp),DIMENSION(ncanop),INTENT(inout) :: fnl_rub          ! fraction of N used for Rubisco per layer
  ! 0.4 Local
  REAL(wp) :: n0                              ! foliar N content of top-most layer
  REAL(wp) :: lc
  INTEGER :: icanop
  REAL(wp) :: structural_n = 34.0_wp          ! non-photosynthetic leaf N default (mmol m-2; Evans, 1988)

  ! calculate top of the canopy N 
  ! standard OCN approach, including the conversion to mmol N / m2(lai)
  IF(lai >= canopy_layer_depth_act(1)) THEN
     n0 = totalN * 1000._wp * kn / (1._wp - exp(-kn * lai))
  ELSE 
     n0 = totalN * 1000._wp
  ENDIF
  ! dummy QUINCY appraoch: N for photosynthesis minus structural N (to be a function of leaf-type and sla)
  ! n0 = MAX( totalN - structural_n * lai, 0.0_wp ) * kn / (1._wp - exp(-kn * lai))
  ! OBS, redefining vcmax2n and jmax2n here temporarily for QUINCY approach
  ! jmax2n=6.0_wp  ! electron-transport limited carboxylation rate per µmol CO2 / mmol N in ET
  ! vcmax2n=1.2_wp ! Rubisco limited carboxylation rate per µmol CO2 / mmol N in RUB

  lc = canopy_layer_depth_act(1) / 2._wp
  DO icanop = 1,ncanop
     IF(canopy_layer_depth_act(icanop).GT.eps8) THEN
        nl_tot(icanop) = n0 * exp(-kn * lc) 
        fnl_chl(icanop) = (k0_fnl_chl - k1_fnl_chl * exp(-kfnl_chl * lc)) / chl2n
        ! standard OCN approach: fraction in indenpent of depth
        fnl_et(icanop) = 0.1_wp * nf     ! OBS: currently only for code consistency set to 0.1*nf 
        fnl_rub(icanop) = 0.1_wp * nf    ! OBS: currently only for code consistency set to 0.1*nf
        ! dummy QUINCY appraoch: fraction of electron transport linked to 
        ! fraction of electron-transport is dependent on N in electron transport versus N in light harvesting
        !fnl_et(icanop) = fnl_chl(icanop) * 24.0_wp/39.8_wp 
        ! fraction of Rubisco is determined as the remainder
        !fnl_rub(icanop) = 1._wp - fnl_chl(icanop) - fnl_et(icanop)
        lc = lc + canopy_layer_depth_act(icanop) / 2._wp & 
             + canopy_layer_depth_act(icanop+1) / 2._wp
    ENDIF
  ENDDO

END SUBROUTINE calc_canopy_n_distribution

SUBROUTINE calc_photosynthesis(temp,press_srf,co2_mixing_ratio,ga,par_sunlit,par_shaded,fleaf_sunlit,beta_air,beta_soil_ps, &
                          beta_soil_gs,fnl_chl,fnl_et,fnl_rub,nfmol,tjmax_opt,tjmax_omega,g0,g1, &
                          agl,anl,rl,gsl,cil)  

  ! Calculation of net photosynthesis and canopy conductance according to Friend & Kiang 2005, J Climate
  ! accounts for temperature acclimation of JMAX according to Friend, 2010 J Exp Bot
    
  USE mo_land_constants
  USE mo_veg_parameters

  IMPLICIT NONE

  ! 0.1 Input
  REAL(wp) :: temp                    ! temperature used for photosynthesis calculation (K) 
  REAL(wp) :: press_srf               ! air pressure (Pa) 
  REAL(wp) :: ga                      ! aerodynamic conductance (m s-1) 
  REAL(wp) :: par_sunlit              ! PAR of sunlit foliage par per canopy layer (µmol m-2 s-1)
  REAL(wp) :: par_shaded              ! PAR of shaded foliage par per canopy layer (µmol m-2 s-1)
  REAL(wp) :: fleaf_sunlit            ! fraction of canopy layer that is sunlit 
  REAL(wp) :: co2_mixing_ratio        ! CO2 concentration of the free air (ppm)
  REAL(wp) :: beta_air                ! scaling factor to account for air humidity effects on conductance
  REAL(wp) :: beta_soil_ps            ! scaling factor to account for soil moisture constraints on conductance
  REAL(wp) :: beta_soil_gs            ! scaling factor to account for soil moisture constraints on conductance
  REAL(wp) :: fnl_et                  ! N specific rate of electron transport (µmol CO2 / mmol N /s )
  REAL(wp) :: fnl_rub                 ! N specific rate of carboylation (µmol CO2 / mmol N /s )
  REAL(wp) :: fnl_chl                 ! fractional chlorophyll content
  REAL(wp) :: nfmol                   ! foliar N content
  REAL(wp) :: tjmax_opt               ! temperature optimum of electron transport
  REAL(wp) :: tjmax_omega             ! shape-parameter of temperature sensitivity of electron transport
  REAL(wp) :: g0                      ! offset of the stomatal conductance ~ An relationship 
  REAL(wp) :: g1                      ! slope of the stomatal conductance ~ An relationship
  ! 0.3 Output
  REAL(wp),INTENT(OUT) :: agl         ! gross photsynthesis of layer (µmol CO2 / m2 /s) 
  REAL(wp),INTENT(OUT) :: anl         ! net photsynthesis of layer (µmol CO2 / m2 /s)
  REAL(wp),INTENT(OUT) :: rl          ! net canopy layer respiration (µmol CO2 / m2 /s)
  REAL(wp),INTENT(OUT) :: gsl         ! canopy conductance of layer (m / s) !! TO BE REDEFINED AS mmol / m2 / s??
  REAL(wp),INTENT(OUT) :: cil         ! internal leaf concentration of layer (ppm)
  ! 0.4 Local
  REAL(wp) :: kc,ko,km                ! Michaelis-Menten Kinetic parameters for Rubisco (Pa)
  REAL(wp) :: pcp                     ! Photosynthetic compensation point (Pa)
  REAL(wp) :: f_resp                  ! N-specific layer respiration rate (µmol CO2 / mmol N /s) 
  REAL(wp) :: kt3                     ! temperature scalar for PeP limited-photosynthesis
  REAL(wp) :: nsat,nlim               ! foliar N that is light saturated / limited (mmol N)
  REAL(wp) :: n1,n2                   ! N specific rates of electron-transport/Rubisco limited PS
  REAL(wp) :: m1,m2                   ! ci sensitivity of electron-transport/Rubisco limited PS
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

  ! static energy (R_gas: J mol-1 K-1 -> rtck: J mol-1)
  rtck = r_gas * temp 
  ! conversion of ppm to Pa
  ppm2Pa = press_srf * 1.e-06
  
  !  Net canopy layer respiration (umol[CO2]/m2/s)
  rl = rl2mmolN * beta_soil_ps * exp(E0kt - E1kt/rtck) * nfmol

  ! initialise output for no-light condition
  IF(par_sunlit + par_shaded <= eps8) THEN
    agl = 0.0_wp 
    anl = agl - rl
    gsl = MAX(g0,g0 + g1 * anl * beta_air * beta_soil_gs / co2_mixing_ratio) * rtck / press_srf 
    cil = co2_mixing_ratio * ppm2Pa
  ELSE !  calculate photosynthesis only if it's day

     ! Michaelis-Menten constant of Rubisco for CO2 (Pa)
     kc = exp(E0kc - E1kc/rtck) * ppm2Pa 
     ! Michaelis-Menten constant of Rubisco for O2 (kPa)
     ko = exp(E0ko - E1ko/rtck) * ppm2Pa 
     ! Combined Rubisco Km for CO2 (Pa); kc (1 + pO2/ko) is part of equation 5 in Kull and Kruijt 1998
     km = kc * (1.0_wp + pO2/ko)
     ! Photorespiratory compensation point (Pa[CO2])
     pcp = exp(E0pcp - E1pcp/rtck) * ppm2Pa 

     ! N specific rate of electron transport rate
     n1 = jmax2n * fnl_et * exp(-(((temp - TZero - tjmax_opt)/tjmax_omega)**2)) * beta_soil_ps

     ! N specific rate of Rubisco turnover 
     n2 = vcmax2n * fnl_rub * exp(E0v - E1v/rtck) * beta_soil_ps

     ! temperature function for PeP C4 photosynthesis 
     kt3 = 2**((temp - TZero - Tref_pep )/ Tbase_pep)

     ! Initial guess of leaf internal CO2 concentration
     cil = co2_mixing_ratio * CiCa_default_C3 

     ! iteration to find correct combination of cil, Anl and gsl
     anl_old = 1000.0_wp
     danl = 1.0_wp
     it = 1
     last_call = .FALSE.
     DO WHILE (it < ps_it_max .AND. danl > eps4 ) 

        ! ci sensitivity of Electron-transport limited photosynthesis
        m1 = cil/(cil + 2.0_wp * pcp)

        ! ci sensitivity of Rubisco-limited photosynthesis
        m2 = cil/(cil + km)
 
        ! N-specific light saturated carboxylation (umol[CO2]/mmol[N]/s)
        msat = MIN(m1 * n1, m2 * n2)

        ! nlim factor (part of equation 12 (Kull and Kruijt 1998))
        nco = msat/(alpha_i * ka * m1)

        ! A couple of factors to calculate photosynthesis
        !  parts of equations 14-16 of Kull and Kruijt 1998
        agsf1 = MAX(1.0_wp - pcp / cil, 0.0_wp)
        agsf2 = agsf1 * alpha_i * m1
        agsf1 = agsf1 * msat

        ! fabsb: fraction of possible PAR absorption in present light conditions
        fabsb = exp(-ka * chl2n * fnl_chl * nfmol)

        !  Case of no direct radiation
        IF (fleaf_sunlit <= eps8) THEN
           !  Theoretical cumulative N at co-limitation by Rubisco and
           !  light absorption in sunlit foliage (mmol[N]/m2))
           nlim = -log(nco/(par_shaded * chl2n * fnl_chl + eps8))/(ka * chl2n * fnl_chl)
           nsat = MAX(MIN(nlim, nfmol), 0.0_wp)

           ! Relative PAR absorption in light-limited region (fraction)
           fabs = exp(-ka * chl2n * fnl_chl * nsat) - fabsb

           ! Rate of gross photosynthesis (umol/m2/s)
           !  sum of photosynthesis in colimited, light saturated PS and PS in the 
           !  light-limited region (assumed to scale linearily with light) 
           agl = agsf1 * nsat + agsf2 * par_shaded * fabs

        ELSE !  Case of direct and diffuse radiation
          !  Theoretical cumulative N at co-limitation by Rubisco and
          !  light absorption in sunlit foliage (mmol[N]/m2))
          nlim = -log(nco/(par_sunlit * chl2n * fnl_chl + eps8))/(ka * chl2n * fnl_chl)
          nsat = MAX(MIN(nlim, nfmol), 0.0_wp)

          !  Relative PAR absorption in light-limited region of sunlit foliage (fraction)
          fabs = exp(-ka * chl2n * fnl_chl * nsat) - fabsb

          ! Rate of gross photosynthesis in sunlit foliage (umol/m2/s)
          agsl = agsf1 * nsat + agsf2 * par_sunlit * fabs

          ! Theoretical cumulative N at co-limitation by Rubisco and
          !  light absorption in shaded foliage (mmol[N]/m2))
          nlim = -log(nco/(par_shaded * chl2n * fnl_chl + eps8))/(ka * chl2n * fnl_chl)
          nsat = MAX(MIN(nlim, nfmol), 0.0_wp)

          !  Relative PAR absorption in light-limited region of shaded foliage (fraction)
          fabs = exp(-ka * chl2n * fnl_chl * nsat) - fabsb

          ! Rate of gross photosynthesis in shaded foliage (umol/m2/s)
          agsh = agsf1 * nsat + agsf2 * par_shaded* fabs

          ! Mean gross photosynthesis over sunlit and shaded fractions (umol/m2/s)
          agl = fleaf_sunlit * agsl + (1.0_wp - fleaf_sunlit) * agsh

        ENDIF
        
        ! Net layer photosynthesis (umol[CO2]/m2/s)
        anl = agl - rl

        ! calculate canopy conductance (m/s)
        gsl = MAX(g0, g0 + g1 * anl * beta_air * beta_soil_gs / ( co2_mixing_ratio )) * & 
              rtck / press_srf
        ! New internal CO2 from balance of supply and demand
        ! ci = cs - An * r_tot (Friend 2005 eq. 7)
        ! Friend 2005 eq. 8: r_tot = ( 1.42 / gsl + 1.60 / (speed * q_cdrag)) 
        cil = co2_mixing_ratio * ppm2Pa - anl * rtck * 1e-06 * (1.42_wp/gsl + 1.6_wp/ga)
        cil = MAX (cil, pcp) 

        ! change in An since last iteration and save new anl 
        danl = abs(anl-anl_old)
        anl_old = anl

        ! break-off point for numerical stability at low conductance
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

  ENDIF

  ! for output purposes, convert ci to ppm
  cil = cil / ppm2Pa

END SUBROUTINE calc_photosynthesis

  FUNCTION calc_atm_moisture_stress(temp,q_air,press_srf) RESULT(beta)
     ! Function to calculate atmospheric drought stress factor
     USE mo_kind, ONLY: wp
     USE mo_land_constants, ONLY: calc_spec_humidity_sat,eps_vpd
     IMPLICIT NONE
     ! 0.1 Input
     REAL(wp),INTENT(in) :: temp               ! temperature in K
     REAL(wp),INTENT(in) :: q_air              ! specific air humidity (g g-1)
     REAL(wp),INTENT(in) :: press_srf          ! air pressure (Pa)
     ! 0.3 Output
     REAL(wp) :: beta                          ! atmospheric moisture scalar for
                                               ! use in conductance calculation
     ! 0.4 Local
     REAL(wp) :: q_sat                         ! saturating specific humidity 

     ! saturating specific humidity [kg kg-1]
     q_sat = calc_spec_humidity_sat(temp,press_srf) 
     ! scalar for atmospheric humidity as relative humidity (fractions)
     beta = q_air/q_sat 

  END FUNCTION calc_atm_moisture_stress

  FUNCTION calc_soil_moisture_stress(theta,lambda,weibull_k) RESULT(beta)
     ! Function to calculate soil moisture stress on photosynthesis or conductance
     ! given current soil moisture
     ! Note: this should account for rooting and soil moisture profile in the long-term
     USE mo_kind, ONLY : wp
     IMPLICIT NONE
     ! 0.1 Input
     REAL(wp),INTENT(in) :: theta       ! soil plant available water (%)
     REAL(wp),INTENT(in) :: lambda      ! slope parameter of Weibull function
     REAL(wp),INTENT(in) :: weibull_k   ! shape parameter of Weibull function
     ! 0.3 Output
     REAL(wp) :: beta
     ! 1. Moisture stress now as a Weibull function to avoid MIN/MAX statements
     beta = 1.0_wp - exp(-(lambda*theta/100._wp)**weibull_k)
  END FUNCTION calc_soil_moisture_stress

END MODULE mo_veg_photosynthesis
