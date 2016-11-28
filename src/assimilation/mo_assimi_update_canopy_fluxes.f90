!>
!! @brief Contains the routines for the Task update_canopy_fluxes (photosynthesis related processes)
!!
!! @author
!!  Sönke Zaehle, MPI-BGC Jena 
!!
!! @par Revision History
!! First version by Sönke Zaehle (2015-12-04)
!! 2nd version by Jan Engel (20160719): svn move vegetation/mo_veg_photosynthesis.f90 assimilation/mo_assimi_update_canopy_fluxes.f90
!!
MODULE mo_assimi_update_canopy_fluxes

#include "dsl_jsbach.inc"

  USE mo_kind,                  ONLY: wp, eps4, eps8
  USE mo_base_memory,           ONLY: t_slm_tile, t_model
  USE mo_base_memory,           ONLY: t_slm_soil
  USE mo_jsb_lctlib_class,      ONLY: t_lctlib_element
  USE mo_jsb_process_class,     ONLY: ATM2LAND_, LAND2ATM_, SEB_, TURB_, SOIL_, HYDRO_, VEG_, HD_, RAD_, ASSIMI_, PHENO_
  USE mo_assimi_constants
  USE mo_assimi_config_class,   ONLY: assimi_config
  USE mo_veg_constants,         ONLY: itree, igrass, ic3phot, ic4phot, ievergreen, iraingreen, isummergreen, iperennial, &
                                      fmaint_rate_base 

  ! Use of process memories
  dsl4jsb_Use_memory(ATM2LAND_)   !  USE mo_atm2land_memory_class, ONLY: t_atm2land_memory
  dsl4jsb_Use_memory(ASSIMI_)     !  USE mo_assimi_memory_class,   ONLY: t_assimi_memory
  dsl4jsb_Use_memory(VEG_)        !  USE mo_veg_memory_class,      ONLY: t_veg_memory


  IMPLICIT NONE

  PRIVATE
  PUBLIC :: update_canopy_fluxes, discrimination_ps, calc_photosynthesis
  PUBLIC :: set_zero_vegetation_fluxes    ! called in the land model interface prior to all Tasks

CONTAINS

  !-----------------------------------------------------------------------------------------------------
  ! Main Task
  ! 
  !-----------------------------------------------------------------------------------------------------
  !> update canopy fluxes, includes calculation of photosynthesis
  !! 
  !! integrates: \n
  !! CALL calc_photosynthesis \n
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE update_canopy_fluxes(tile)

    USE mo_jsb_class,             ONLY: Get_model
    USE mo_shared_constants
    USE mo_shared_isotope
 
    IMPLICIT NONE
    ! 0.1 InOut
    TYPE(t_slm_tile),       INTENT(inout), TARGET :: tile      !< one tile  =>  *slmdev may remove the TARGET argument when soil  => tile%soil is not needed anymore 
    ! 0.2 Local
    TYPE(t_model),          POINTER        :: model
    TYPE(t_lctlib_element), POINTER        :: lctlib     !< land-cover-type library - parameter across pft's
    TYPE(t_slm_soil),       POINTER        :: soil       ! soil variables per tile
    INTEGER                                :: icanopy     ! for looping across ncanopy layers
    REAL(wp)                               :: hlp1       ! dummy variable 
    ! 0.3 Memory allocation
    ! Declare process configuration and memory Pointers
    dsl4jsb_Def_memory(ATM2LAND_)     !   TYPE(t_atm2land_memory), POINTER :: atm2land__mem
    dsl4jsb_Def_memory(ASSIMI_)       !   TYPE(t_assimi_memory), POINTER   :: assimi__mem
    dsl4jsb_Def_memory(VEG_)          !   TYPE(t_veg_memory), POINTER      :: veg__mem  
    TYPE(t_assimi_memory), POINTER   :: assimi                 ! local pointer to assimi__mem for SLM usability
    TYPE(t_veg_memory), POINTER      :: veg                    ! local pointer to veg__mem for SLM usability
    REAL(wp), POINTER :: t_air
    REAL(wp), POINTER :: q_air
    REAL(wp), POINTER :: press_srf
    REAL(wp), POINTER :: ga
    REAL(wp), POINTER :: co2_mixing_ratio
    REAL(wp), POINTER :: co2_mixing_ratio_C13
    REAL(wp), POINTER :: co2_mixing_ratio_C14
    dsl4jsb_Get_memory(ATM2LAND_)
    dsl4jsb_Get_memory(ASSIMI_)
    dsl4jsb_Get_memory(VEG_)
!    SELECT TYPE (mm => tile%mem(ATM2LAND_)%p) 
!      TYPE IS(t_atm2land_memory)
!        atm2land__mem => mm
!    END SELECT
!    SELECT TYPE (mm => tile%mem(ASSIMI_)%p)
!      TYPE IS(t_assimi_memory)
!        assimi__mem => mm
!    END SELECT
!    SELECT TYPE (mm => tile%mem(VEG_)%p)     
!      TYPE IS(t_veg_memory)
!        veg__mem => mm
!    END SELECT
    assimi  => assimi__mem
    veg     => veg__mem 
    t_air                => atm2land__mem%t_air%ptr
    q_air                => atm2land__mem%q_air%ptr
    press_srf            => atm2land__mem%press_srf%ptr
    ga                   => atm2land__mem%ga%ptr
    co2_mixing_ratio     => atm2land__mem%co2_mixing_ratio%ptr
    co2_mixing_ratio_C13 => atm2land__mem%co2_mixing_ratio_C13%ptr
    co2_mixing_ratio_C14 => atm2land__mem%co2_mixing_ratio_C14%ptr

    model  => Get_model(tile%owner_model_id)
    lctlib => model%lctlib(tile%lcts%lib_id)

   !  bi-directional pointer  *slmdev
   soil  => tile%soil

   ! set output to zero  
   ! *slmdev  necessary ? 
!   assimi%agl%ptr = 0.0_wp 
!   assimi%anl%ptr = 0.0_wp 
!   assimi%rl%ptr  = 0.0_wp
!   assimi%gsl%ptr = 0.0_wp
!   assimi%cil%ptr = 0.0_wp     
   veg%gc%ptr     = 0.0_wp 


   !> 1.0 calculate the moisture scaling factor for photosynthesis
   !! 
   veg%beta_air%ptr     = calc_atm_moisture_stress( t_air                   , &
                                                    q_air                   , &
                                                    press_srf) 
   veg%beta_soil_ps%ptr = calc_soil_moisture_stress(soil%theta              , & 
                                                    lctlib%lambda_soil_ps   , &
                                                    lctlib%weibull_k_soil_ps)
   veg%beta_soil_gs%ptr = calc_soil_moisture_stress(soil%theta              , &
                                                    lctlib%lambda_soil_gs   , &
                                                    lctlib%weibull_k_soil_gs)

   !> 2.0 calculate photosynthesis per layer
   !! 
   DO icanopy = 1,assimi_config%ncanopy
      IF(veg%lai_cl%ptr(icanopy,1) > eps8)THEN

        CALL calc_photosynthesis(t_air                            , &   ! in
                                 press_srf                        , &   ! in
                                 co2_mixing_ratio                 , &   ! in
                                 ga                               , &   ! in
                                 veg%par_sunlit_cl%ptr(icanopy,1)  , &   ! in
                                 veg%par_shaded_cl%ptr(icanopy,1)  , &   ! in
                                 veg%fleaf_sunlit_cl%ptr(icanopy,1), &   ! in
                                 veg%beta_air%ptr                 , &   ! in
                                 veg%beta_soil_ps%ptr             , &   ! in
                                 veg%beta_soil_gs%ptr             , &   ! in
                                 veg%fnl_chl_cl%ptr(icanopy,1)     , &   ! in
                                 veg%fnl_et_cl%ptr(icanopy,1)      , &   ! in
                                 veg%fnl_rub_cl%ptr(icanopy,1)     , &   ! in
                                 veg%fnl_pepc_cl%ptr(icanopy,1)    , &   ! in
                                 veg%nl_tot_cl%ptr(icanopy,1)      , &   ! in
                                 lctlib%ps_pathway                , &   ! in
                                 veg%t_jmax_opt_mavg%ptr          , &   ! in
                                 lctlib%t_jmax_omega              , &   ! in
                                 lctlib%g0                        , &   ! in
                                 lctlib%g1                        , &   ! in
                                 assimi%agl%ptr                   , &   ! out
                                 assimi%anl%ptr                   , &   ! out
                                 assimi%rl%ptr                    , &   ! out
                                 assimi%gsl%ptr                   , &   ! out
                                 assimi%cil%ptr                   , &   ! out
                                 hlp1, hlp1, hlp1)                      ! out

        veg%gpp = veg%gpp + assimi%agl%ptr * veg%lai_cl%ptr(icanopy,1)

        IF(assimi%agl%ptr > eps8) THEN
           hlp1          = discrimination_ps(assimi%cil%ptr     , &
                                             co2_mixing_ratio   , &
                                             'C13'              , &
                                             lctlib%ps_pathway)
           veg%gpp_C13 = veg%gpp_C13 + assimi%agl%ptr * veg%lai_cl%ptr(icanopy,1) * calc_fractionation(co2_mixing_ratio  , &
                                                                                                   co2_mixing_ratio_C13 , &
                                                                                                   hlp1)
           hlp1          = discrimination_ps(assimi%cil%ptr     , &
                                             co2_mixing_ratio   , &
                                             'C14'              , &
                                             lctlib%ps_pathway)
           veg%gpp_C14 = veg%gpp_C14 + assimi%agl%ptr * veg%lai_cl%ptr(icanopy,1) * calc_fractionation(co2_mixing_ratio  , &
                                                                                                   co2_mixing_ratio_C14 , & 
                                                                                                   hlp1                 , &
                                                                                                   scale_correction_C14)
        ENDIF

        veg%an                      = veg%an + assimi%anl%ptr * veg%lai_cl%ptr(icanopy,1) 
        veg%maint_respiration_leaf  = veg%maint_respiration_leaf + assimi%rl%ptr * veg%lai_cl%ptr(icanopy,1) 
        veg%gc                      = veg%gc + assimi%gsl%ptr * veg%lai_cl%ptr(icanopy,1)
        veg%ci                      = veg%ci + assimi%cil%ptr * veg%lai_cl%ptr(icanopy,1)

        ! get the net photosynthesis of the lowest canopy layer at this timestep
        IF(veg%lai_cl%ptr(icanopy,1) > 0.0_wp) veg%an_boc%ptr = assimi%anl%ptr
     ENDIF
  ENDDO

  ! ci
  IF(veg%lai%ptr > eps8) veg%ci = veg%ci / veg%lai

  ! t_jmax_opt
  IF(assimi_config%flag_t_jmax_acclimation) THEN 
     veg%t_jmax_opt%ptr = t_jmax_offset + t_jmax_slope * ( t_air - TZero ) 
     veg%t_jmax_opt%ptr = MIN(MAX(veg%t_jmax_opt%ptr,t_jmax_opt_min),t_jmax_opt_max)
  ELSE 
     veg%t_jmax_opt%ptr = lctlib%t_jmax_opt
  ENDIF

  END SUBROUTINE update_canopy_fluxes

  !-----------------------------------------------------------------------------------------------------
  ! Sub Task to update_canopy_fluxes
  ! 
  !-----------------------------------------------------------------------------------------------------
  !> Calculation of net photosynthesis and canopy conductance
  !! 
  !! according to Friend & Kiang 2005, J Climate \n
  !! accounts for temperature acclimation of JMAX according to Friend, 2010 J Exp Bot
  !-----------------------------------------------------------------------------------------------------
  ELEMENTAL SUBROUTINE calc_photosynthesis( temp                , &   ! in
                                            press_srf           , &   ! in
                                            co2_mixing_ratio    , &   ! in
                                            ga                  , &   ! in
                                            par_sunlit_cl       , &   ! in
                                            par_shaded_cl       , &   ! in
                                            fleaf_sunlit_cl     , &   ! in
                                            beta_air            , &   ! in
                                            beta_soil_ps        , &   ! in
                                            beta_soil_gs        , &   ! in
                                            fnl_chl_cl          , &   ! in
                                            fnl_et_cl           , &   ! in
                                            fnl_rub_cl          , &   ! in
                                            fnl_pepc_cl         , &   ! in
                                            nfmol               , &   ! in
                                            ps_pathway          , &   ! in
                                            t_jmax_opt          , &   ! in
                                            t_jmax_omega        , &   ! in
                                            g0                  , &   ! in
                                            g1                  , &   ! in
                                            agl                 , &   ! out
                                            anl                 , &   ! out
                                            rl                  , &   ! out
                                            gsl                 , &   ! out
                                            cil                 , &   ! out
                                            m_rub               , &   ! out
                                            m_et                , &   ! out
                                            m_pepc)                   ! out
    
  USE mo_shared_constants,    ONLY: r_gas, Tzero
  USE mo_assimi_constants
  USE mo_veg_respiration,     ONLY: temperature_response_respiration

  IMPLICIT NONE

  ! 0.1 InOut
  REAL(wp),INTENT(IN) :: temp                    !< temperature used for photosynthesis calculation (K) 
  REAL(wp),INTENT(IN) :: press_srf               !< air pressure (Pa) 
  REAL(wp),INTENT(IN) :: ga                      !< aerodynamic conductance (m s-1) 
  REAL(wp),INTENT(IN) :: par_sunlit_cl              !< PAR of sunlit foliage par per canopy layer (µmol m-2 s-1)
  REAL(wp),INTENT(IN) :: par_shaded_cl              !< PAR of shaded foliage par per canopy layer (µmol m-2 s-1)
  REAL(wp),INTENT(IN) :: fleaf_sunlit_cl            !< fraction of canopy layer that is sunlit 
  REAL(wp),INTENT(IN) :: co2_mixing_ratio        !< CO2 concentration of the free air (ppm)
  REAL(wp),INTENT(IN) :: beta_air                !< scaling factor to account for air humidity effects on conductance
  REAL(wp),INTENT(IN) :: beta_soil_ps            !< scaling factor to account for soil moisture constraints on conductance
  REAL(wp),INTENT(IN) :: beta_soil_gs            !< scaling factor to account for soil moisture constraints on conductance
  REAL(wp),INTENT(IN) :: fnl_et_cl                  !< fraction of foliar N allocated to electron transport 
  REAL(wp),INTENT(IN) :: fnl_rub_cl                 !< fraction of foliar N allocated to Rubisco
  REAL(wp),INTENT(IN) :: fnl_pepc_cl                !< fraction of foliar N allocated to PEPC complex (C4 only)
  REAL(wp),INTENT(IN) :: fnl_chl_cl                 !< fraction of foliar N allocated to Chlorophyll
  REAL(wp),INTENT(IN) :: nfmol                   !< total foliar N content
  INTEGER,INTENT(IN)  :: ps_pathway              !< photosynthetic pathway 
  REAL(wp),INTENT(IN) :: t_jmax_opt              !< temperature optimum of electron transport
  REAL(wp),INTENT(IN) :: t_jmax_omega             !< shape-parameter of temperature sensitivity of electron transport
  REAL(wp),INTENT(IN) :: g0                      !< offset of the stomatal conductance ~ An relationship 
  REAL(wp),INTENT(IN) :: g1                      !< slope of the stomatal conductance ~ An relationship
  REAL(wp),INTENT(OUT) :: agl                    !< gross photsynthesis of layer (µmol CO2 / m2 /s) 
  REAL(wp),INTENT(OUT) :: anl                    !< net photsynthesis of layer (µmol CO2 / m2 /s)
  REAL(wp),INTENT(OUT) :: rl                     !< net canopy layer respiration (µmol CO2 / m2 /s)
  REAL(wp),INTENT(OUT) :: gsl                    !< canopy conductance of layer (m / s) 
  REAL(wp),INTENT(OUT) :: cil                    !< leaf internal CO2 concentration of layer (ppm)
  REAL(wp),INTENT(OUT) :: m_rub                  !< CO2 limitation of Rubisco factor for calculating optimal N fraction
  REAL(wp),INTENT(OUT) :: m_et                   !< CO2 limitation of ET factor for calculating optimal N fraction
  REAL(wp),INTENT(OUT) :: m_pepc                 !< CO2 limitation for PEPC factor for calculating optimal N fraction
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
  IF(par_sunlit_cl + par_shaded_cl <= eps8) THEN
    agl    = 0.0_wp 
    anl    = agl - rl
    gsl    = MAX(g0,g0 + g1 * anl * beta_air * beta_soil_gs / co2_mixing_ratio) * rtck / press_srf 
    cil    = co2_mixing_ratio * ppm2Pa
    m_et   = 0.0_wp
    m_rub  = 0.0_wp
    m_pepc = 0.0_wp
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
     n1 = jmax2n * fnl_et_cl * exp(-(((temp - Tzero - t_jmax_opt)/t_jmax_omega)**2)) * beta_soil_ps

     !> N specific rate of Rubisco turnover 
     n2 = vcmax2n * fnl_rub_cl * exp(E0v - E1v/rtck) * beta_soil_ps

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
           nco = msat/(alpha_i * ka * m1)

           !> A couple of factors to calculate photosynthesis
           !!  parts of equations 14-16 of Kull and Kruijt 1998
           agsf1 = MAX(1.0_wp - pcp / cil, 0.0_wp)
           agsf2 = agsf1 * alpha_i * m1
           agsf1 = agsf1 * msat

        CASE (ic4phot)
           !> N specific rate of PePc limited photosynthessis 
           n3 = pepc2n * fnl_pepc_cl * 2**((temp - TZero - Tref_pepc )/ Tbase_pepc)

           !> ci sensitivity of Electron-transport limited photosynthesis
           m1 = ci_max/(ci_max + 2.0_wp * pcp)

           !> ci sensitivity of Rubisco-limited photosynthesis
           m2 = ci_max/(ci_max + km)
 
           !> ci sensitivity of PePc-limited photosynthesis
           m3 = cil / press_srf

           !> N-specific light saturated carboxylation (umol[CO2]/mmol[N]/s)
           msat = MIN(MIN(m1 * n1, m2 * n2), m3 * n3)
           !> nlim factor (part of equation 12 (Kull and Kruijt 1998))
           nco = msat/(alpha_i * ka * m1)

           !> A couple of factors to calculate photosynthesis
           !!  parts of equations 14-16 of Kull and Kruijt 1998
           agsf1 = MAX(1.0_wp - pcp / ci_max, 0.0_wp)
           agsf2 = agsf1 * alpha_i * m1
           agsf1 = agsf1 * msat
  
        END SELECT

        ! output factors for optimal fraction calculation
        m_et = m1 * n1 / fnl_et_cl
        m_rub = m2 * n2 / fnl_rub_cl
        IF (fnl_pepc_cl>eps8) THEN
           m_pepc = m3 * n3 / fnl_pepc_cl
        ELSE
           m_pepc = 0._wp
        ENDIF

        !> fabsb: fraction of possible PAR absorption in present light conditions
        fabsb = exp(-ka * chl2n * fnl_chl_cl * nfmol)

        !> @par IF: Case of no direct radiation
        IF (fleaf_sunlit_cl <= eps8) THEN
           !>  Theoretical cumulative N at co-limitation by Rubisco and
           !!  light absorption in sunlit foliage (mmol[N]/m2))
           nlim = -log(nco/(par_shaded_cl * chl2n * fnl_chl_cl + eps8))/(ka * chl2n * fnl_chl_cl)
           nsat = MAX(MIN(nlim, nfmol), 0.0_wp)

           !> Relative PAR absorption in light-limited region (fraction)
           fabs = exp(-ka * chl2n * fnl_chl_cl * nsat) - fabsb

           !> Rate of gross photosynthesis (umol/m2/s)
           !!  sum of photosynthesis in colimited, light saturated PS and PS in the 
           !!  light-limited region (assumed to scale linearily with light) 
           agl = agsf1 * nsat + agsf2 * par_shaded_cl * fabs

        !> @par ELSE: Case of direct and diffuse radiation
        ELSE 
          !>  Theoretical cumulative N at co-limitation by Rubisco and
          !!  light absorption in sunlit foliage (mmol[N]/m2))
          nlim = -log(nco/(par_sunlit_cl * chl2n * fnl_chl_cl + eps8))/(ka * chl2n * fnl_chl_cl)
          nsat = MAX(MIN(nlim, nfmol), 0.0_wp)

          !>  Relative PAR absorption in light-limited region of sunlit foliage (fraction)
          fabs = exp(-ka * chl2n * fnl_chl_cl * nsat) - fabsb

          !> Rate of gross photosynthesis in sunlit foliage (umol/m2/s)
          agsl = agsf1 * nsat + agsf2 * par_sunlit_cl * fabs

          !> Theoretical cumulative N at co-limitation by Rubisco and
          !!  light absorption in shaded foliage (mmol[N]/m2))
          nlim = -log(nco/(par_shaded_cl * chl2n * fnl_chl_cl + eps8))/(ka * chl2n * fnl_chl_cl)
          nsat = MAX(MIN(nlim, nfmol), 0.0_wp)

          !> Relative PAR absorption in light-limited region of shaded foliage (fraction)
          fabs = exp(-ka * chl2n * fnl_chl_cl * nsat) - fabsb

          !> Rate of gross photosynthesis in shaded foliage (umol/m2/s)
          agsh = agsf1 * nsat + agsf2 * par_shaded_cl* fabs

          !> Mean gross photosynthesis over sunlit and shaded fractions (umol/m2/s)
          agl = fleaf_sunlit_cl * agsl + (1.0_wp - fleaf_sunlit_cl) * agsh
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
        !! 
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
  ! Sub Task to update_canopy_fluxes
  ! 
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

     !> saturating specific humidity [kg kg-1]
     q_sat = calc_spec_humidity_sat(temp,press_srf) 

     !> vapour pressure deficit [kPa] 
     vpd = MAX(q_sat - q_air,0.0_wp) * press_srf / eps_vpd / 1000._wp

     !> scalar for atmospheric humidity as relative humidity (fractions)
     beta = q_air/q_sat ! Ball & Berry
     !beta = MIN(5._wp,1._wp / SQRT(vpd + eps8)) ! Medlyn et al. 2011
  END FUNCTION calc_atm_moisture_stress


  !-----------------------------------------------------------------------------------------------------
  ! Sub Task to update_canopy_fluxes
  ! 
  !-----------------------------------------------------------------------------------------------------
  !> calculate soil moisture stress on photosynthesis or conductance given current soil moisture
  !! 
  !! @todo (calc_soil_moisture_stress) this should account for rooting and soil moisture profile in the long-term
  !-----------------------------------------------------------------------------------------------------
  FUNCTION calc_soil_moisture_stress(theta,lambda,weibull_k) RESULT(beta)

    IMPLICIT NONE
    ! 0.1 InOut
    REAL(wp),INTENT(in) :: theta       !< soil plant available water (%)
    REAL(wp),INTENT(in) :: lambda      !< slope parameter of Weibull function
    REAL(wp),INTENT(in) :: weibull_k   !< shape parameter of Weibull function
    REAL(wp) :: beta                   !< veg\%beta_soil_ps / veg\%beta_soil_gs

    !> 1.0 Moisture stress now as a Weibull function to avoid MIN/MAX statements
    !!
    beta = 1.0_wp - exp(-(lambda*theta/100._wp)**weibull_k)
  END FUNCTION calc_soil_moisture_stress


  !-----------------------------------------------------------------------------------------------------
  ! Sub Task to update_canopy_fluxes
  ! 
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
  PURE ELEMENTAL FUNCTION discrimination_ps(ci, ca, isotope_name, ps_pathway) RESULT(discrimination)

    USE mo_assimi_constants
    IMPLICIT NONE
    ! 0.1 Inout 
    REAL(wp),     INTENT(in) :: ci                     !< leaf internal CO2 concentration (ppm) 
    REAL(wp),     INTENT(in) :: ca                     !< atmospheric CO2 concentration (ppm) 
    CHARACTER(3), INTENT(in) :: isotope_name           !< isotope C13 or C14
    INTEGER,      INTENT(in) :: ps_pathway             !< photosynthetic pathway (C3 or C4)
    REAL(wp)                 :: discrimination         !< isotopic discrimination (per mill)

    !> SELECT PS pathway C3 C4
    SELECT CASE(ps_pathway)
    CASE (ic3phot)
      !> SELECT isotope C13 C14
      SELECT CASE(isotope_name)
        CASE ('C13')
          discrimination = discr_ps_a_C13 + (discr_ps_b_C13 - discr_ps_a_C13) * ci/ca
        CASE ('C14')
          discrimination = discr_ps_a_C14 + (discr_ps_b_C14 - discr_ps_a_C14) * ci/ca
      END SELECT
    CASE (ic4phot)
      !> SELECT isotope C13 C14
      SELECT CASE(isotope_name)
        CASE ('C13')
          discrimination = discr_ps_a_C13 + (discr_ps_c_C13 + discr_ps_phi * discr_ps_b_C13 - discr_ps_a_C13) * ci/ca
        CASE ('C14')
          discrimination = discr_ps_a_C14 + (discr_ps_c_C14 + discr_ps_phi * discr_ps_b_C14 - discr_ps_a_C14) * ci/ca
      END SELECT
    END SELECT

  END FUNCTION discrimination_ps


  !-----------------------------------------------------------------------------------------------------
  !> set all vegetation fluxes to zero
  !! 
  !! called in the land model interface prior to all Tasks
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE set_zero_vegetation_fluxes(tile)

    IMPLICIT NONE
    ! 0.1 InOut
    TYPE(t_slm_tile),INTENT(inout)      :: tile           !< one tile
    ! 0.2 Local
    REAL(wp) :: flux_zero
    ! 0.3 Memory allocation
    TYPE(t_veg_memory), POINTER         :: veg__mem   ! dsl4jsb_Def_memory(VEG_)
    TYPE(t_veg_memory), POINTER         :: veg        ! local pointer to veg__mem for SLM usability
    SELECT TYPE (mm => tile%mem(VEG_)%p)              ! dsl4jsb_Get_memory(VEG_)
      TYPE IS(t_veg_memory)
        veg__mem => mm
    END SELECT
    veg => veg__mem 

    flux_zero = 0.0_wp


    ! t_slm_struct_veg_element
    CALL veg%reserve_use%init(flux_zero)
    CALL veg%establishment%init(flux_zero)
    ! t_slm_struct_veg_pool
    CALL veg%growth%init(flux_zero)
    CALL veg%litterfall%init(flux_zero)
    ! variables
    CALL veg%gpp%init(flux_zero)                           !< gross photosynthesis (µmol CO2 m-2 s-1)
    CALL veg%an%init(flux_zero)                            !< net photosynthesis (µmol CO2 m-2 s-1)
    CALL veg%gpp_c13%init(flux_zero)                       !< gross photosynthesis (µmol 13CO2 m-2 s-1)
    CALL veg%gpp_c14%init(flux_zero)                       !< gross photosynthesis (µmol 14CO2 m-2 s-1)
    CALL veg%maint_respiration_leaf%init(flux_zero)        !< foliar maintenance respiration (µmol CO2 m-2 s-1)
    CALL veg%maint_respiration_pot%init(flux_zero)         !< maintenance respiration (µmol CO2 m-2 s-1)
    CALL veg%maint_respiration%init(flux_zero)             !< maintenance respiration (µmol CO2 m-2 s-1)
    CALL veg%maint_respiration_c13%init(flux_zero)         !< maintenance respiration (µmol 13CO2 m-2 s-1)
    CALL veg%maint_respiration_c14%init(flux_zero)         !< maintenance respiration (µmol 14CO2 m-2 s-1)
    CALL veg%growth_respiration%init(flux_zero)            !< growth respiration (µmol CO2 m-2 s-1)
    CALL veg%growth_respiration_c13%init(flux_zero)        !< growth respiration (µmol 13CO2 m-2 s-1)
    CALL veg%growth_respiration_c14%init(flux_zero)        !< growth respiration (µmol 14CO2 m-2 s-1)
    CALL veg%n_transform_respiration%init(flux_zero)       !< respiration associated with N transformation (µmol CO2 m-2 s-1)
    CALL veg%n_fixation_respiration%init(flux_zero)        !< respiration associated with N fixation (µmol CO2 m-2 s-1)
    CALL veg%n_processing_respiration%init(flux_zero)      !< respiration associated with N fixation (µmol CO2 m-2 s-1)
    CALL veg%n_processing_respiration_c13%init(flux_zero)  !< respiration associated with N fixation (µmol 13CO2 m-2 s-1)
    CALL veg%n_processing_respiration_c14%init(flux_zero)  !< respiration associated with N fixation (µmol 14CO2 m-2 s-1)
    CALL veg%npp%init(flux_zero)                           !< net primary production (µmol CO2 m-2 s-1)
    CALL veg%npp_c13%init(flux_zero)                       !< net primary production (µmol 13CO2 m-2 s-1)
    CALL veg%npp_c14%init(flux_zero)                       !< net primary production (µmol 14CO2 m-2 s-1)
    CALL veg%uptake_n%init(flux_zero)                      !< ??
    CALL veg%uptake_n15%init(flux_zero)                    !< ??
    CALL veg%n_fixation%init(flux_zero)                    !< ??
    CALL veg%n_fixation_n15%init(flux_zero)                !< ??
    CALL veg%uptake_p%init(flux_zero)                      !< ??
    CALL veg%recycling_leaf_n%init(flux_zero)              !< net flux of N from leaf to labile (senescence and maintenance) (mol m-2 timestep-1)
    CALL veg%recycling_leaf_n15%init(flux_zero)            !< net flux of 15N from leaf to labile (senescence and maintenance) (mol m-2 timestep-1) 
    CALL veg%recycling_leaf_p%init(flux_zero)              !< net flux of P from leaf to labile (senescence and maintenance) (mol m-2 timestep-1)
    CALL veg%recycling_fine_root_n%init(flux_zero)         !< net flux of N from fine roots to labile (senescence and maintenance) (mol m-2 timestep-1)
    CALL veg%recycling_fine_root_n15%init(flux_zero)       !< net flux of 15N from fine roots to labile (senescence and maintenance) (mol m-2 timestep-1) 
    CALL veg%recycling_fine_root_p%init(flux_zero)         !< net flux of P from fine roots to labile (senescence and maintenance) (mol m-2 timestep-1)
    CALL veg%recycling_sap_wood_n%init(flux_zero)          !< net flux of N from sapwood to labile (senescence and maintenance) (mol m-2 timestep-1)
    CALL veg%recycling_sap_wood_n15%init(flux_zero)        !< net flux of 15N from sapwood to labile (senescence and maintenance) (mol m-2 timestep-1) 
    CALL veg%recycling_sap_wood_p%init(flux_zero)          !< net flux of P from sapwood to labile (senescence and maintenance) (mol m-2 timestep-1)
    CALL veg%transpiration%init(flux_zero)                 !< transpiration in kg m-2 s-1
    CALL veg%ci%init(flux_zero)                            !< ??
    CALL veg%delta_dens_ind%init(flux_zero)                !< ??
    CALL veg%unit_transpiration%init(flux_zero)            !< ??

  END SUBROUTINE set_zero_vegetation_fluxes

END MODULE mo_assimi_update_canopy_fluxes
