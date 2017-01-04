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
  USE mo_jsb_process_class,     ONLY: ATM2LAND_, LAND2ATM_, SEB_, TURB_, SOIL_, HYDRO_, VEG_, HD_, RAD_, ASSIMI_, PHENO_, SB_
  USE mo_assimi_constants,      ONLY: ic3phot, ic4phot, t_jmax_offset, t_jmax_slope, t_jmax_opt_min, t_jmax_opt_max
  USE mo_assimi_config_class,   ONLY: assimi_config
  USE mo_veg_constants,         ONLY: itree, igrass, fmaint_rate_base


  ! Use of process memories
  dsl4jsb_Use_memory(ATM2LAND_)   !  USE mo_atm2land_memory_class, ONLY: t_atm2land_memory
  dsl4jsb_Use_memory(ASSIMI_)     !  USE mo_assimi_memory_class,   ONLY: t_assimi_memory
  dsl4jsb_Use_memory(VEG_)        !  USE mo_veg_memory_class,      ONLY: t_veg_memory


  IMPLICIT NONE

  PRIVATE
  PUBLIC :: update_canopy_fluxes, discrimination_ps, calc_photosynthesis


  CHARACTER(len=*), PARAMETER :: modname = 'mo_assimi_update_canopy_fluxes'

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

    ! Use of process memories 
    dsl4jsb_Use_memory(RAD_)        !  USE mo_rad_memory_class,      ONLY: t_rad_memory
    dsl4jsb_Use_memory(SOIL_)       !  USE mo_soil_memory_class,     ONLY: t_soil_memory

    IMPLICIT NONE
    ! 0.1 InOut
    TYPE(t_slm_tile),       INTENT(inout), TARGET :: tile      !< one tile  =>  *slmdev may remove the TARGET argument when soil  => tile%soil is not needed anymore 
    ! 0.2 Local
    TYPE(t_model),          POINTER        :: model
    TYPE(t_lctlib_element), POINTER        :: lctlib     !< land-cover-type library - parameter across pft's
    INTEGER                                :: icanopy     ! for looping across ncanopy layers
    REAL(wp)                               :: hlp1       ! dummy variable 
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':update_canopy_fluxes'
    ! 0.3 Memory allocation
    ! Declare process configuration and memory Pointers
    dsl4jsb_Def_memory(ATM2LAND_)     !   TYPE(t_atm2land_memory), POINTER :: atm2land__mem
    dsl4jsb_Def_memory(ASSIMI_)       !   TYPE(t_assimi_memory), POINTER   :: assimi__mem
    dsl4jsb_Def_memory(VEG_)          !   TYPE(t_veg_memory), POINTER      :: veg__mem  
    dsl4jsb_Def_memory(RAD_)                          !     TYPE(t_rad_memory),      POINTER :: rad__mem
    dsl4jsb_Def_memory(SOIL_)                         !     TYPE(t_soil_memory),      POINTER :: soil__mem
    TYPE(t_assimi_memory), POINTER   :: assimi                 ! local pointer to assimi__mem for SLM usability
    TYPE(t_veg_memory), POINTER      :: veg                    ! local pointer to veg__mem for SLM usability
    TYPE(t_rad_memory),       POINTER   :: rad        ! local pointer to rad__mem for SLM usability
    TYPE(t_soil_memory),      POINTER   :: soil       ! local pointer to soil__mem for SLM usability
    REAL(wp), POINTER :: t_air
    REAL(wp), POINTER :: q_air
    REAL(wp), POINTER :: press_srf
    REAL(wp), POINTER :: ga
    REAL(wp), POINTER :: co2_mixing_ratio
    REAL(wp), POINTER :: co2_mixing_ratio_C13
    REAL(wp), POINTER :: co2_mixing_ratio_C14

    REAL(wp), POINTER :: gross_assimilation, gross_assimilation_cl
    REAL(wp), POINTER :: gross_assimilation_C13, gross_assimilation_C14
    REAL(wp), POINTER :: net_assimilation, net_assimilation_cl, net_assimilation_boc
    REAL(wp), POINTER :: maint_respiration_leaf, maint_respiration_leaf_cl
    REAL(wp), POINTER :: canopy_cond, canopy_cond_cl
    REAL(wp), POINTER :: co2_conc_leaf, co2_conc_leaf_cl
    REAL(wp), POINTER :: beta_air, beta_soil_ps, beta_soil_gs
    REAL(wp), POINTER :: t_jmax_opt 

    dsl4jsb_Get_memory(ATM2LAND_)
    dsl4jsb_Get_memory(ASSIMI_)
    dsl4jsb_Get_memory(VEG_)
    dsl4jsb_Get_memory(RAD_)
    dsl4jsb_Get_memory(SOIL_)

    veg     => veg__mem 
    rad     => rad__mem
    soil    => soil__mem

    t_air                => atm2land__mem%t_air%ptr
    q_air                => atm2land__mem%q_air%ptr
    press_srf            => atm2land__mem%press_srf%ptr
    ga                   => atm2land__mem%ga%ptr
    co2_mixing_ratio     => atm2land__mem%co2_mixing_ratio%ptr
    co2_mixing_ratio_C13 => atm2land__mem%co2_mixing_ratio_C13%ptr
    co2_mixing_ratio_C14 => atm2land__mem%co2_mixing_ratio_C14%ptr

    gross_assimilation => assimi__mem%gross_assimilation%ptr
    gross_assimilation_C13 => assimi__mem%gross_assimilation_C13%ptr
    gross_assimilation_C14 => assimi__mem%gross_assimilation_C14%ptr
    gross_assimilation_cl => assimi__mem%gross_assimilation_cl%ptr
    net_assimilation => assimi__mem%net_assimilation%ptr
    net_assimilation_cl => assimi__mem%net_assimilation_cl%ptr
    net_assimilation_boc => assimi__mem%net_assimilation_boc%ptr
    maint_respiration_leaf => assimi__mem%maint_respiration_leaf%ptr
    maint_respiration_leaf_cl => assimi__mem%maint_respiration_leaf_cl%ptr
    canopy_cond => assimi__mem%canopy_cond%ptr
    canopy_cond_cl => assimi__mem%canopy_cond_cl%ptr
    co2_conc_leaf => assimi__mem%co2_conc_leaf%ptr
    co2_conc_leaf_cl => assimi__mem%co2_conc_leaf_cl%ptr
    beta_air => assimi__mem%beta_air%ptr
    beta_soil_gs => assimi__mem%beta_soil_gs%ptr
    beta_soil_ps => assimi__mem%beta_soil_ps%ptr
    t_jmax_opt => assimi__mem%t_jmax_opt%ptr

    model  => Get_model(tile%owner_model_id)
    lctlib => model%lctlib(tile%lcts%lib_id)

   ! set output to zero  
   gross_assimilation     = 0.0_wp
   gross_assimilation_C13 = 0.0_wp
   gross_assimilation_C14 = 0.0_wp
   net_assimilation       = 0.0_wp
   maint_respiration_leaf        = 0.0_wp
   canopy_cond            = 0.0_wp
   co2_conc_leaf          = 0.0_wp

   !> 1.0 calculate the moisture scaling factor for photosynthesis
   !! 
   beta_air = calc_atm_moisture_stress(t_air, q_air, press_srf) 
   beta_soil_ps = calc_soil_moisture_stress(soil%w_soil_root_theta%ptr(1,1), &
                                            lctlib%lambda_soil_ps   , &
                                            lctlib%weibull_k_soil_ps)
   beta_soil_gs = calc_soil_moisture_stress(soil%w_soil_root_theta%ptr(1,1), &
                                            lctlib%lambda_soil_gs   , &
                                            lctlib%weibull_k_soil_gs)

   !> 2.0 calculate photosynthesis per layer
   !! 
   DO icanopy = 1,assimi_config%ncanopy
      IF(veg%lai_cl%ptr(icanopy,1) > eps8)THEN

        CALL calc_photosynthesis(t_air                            , & ! in
                                 press_srf                        , &
                                 co2_mixing_ratio                 , &
                                 ga                               , &
                                 veg%t_air_tacclim_mavg%ptr       , &
                                 rad%ppfd_sunlit_cl%ptr(icanopy,1), &
                                 rad%ppfd_shaded_cl%ptr(icanopy,1), &
                                 veg%fleaf_sunlit_cl%ptr(icanopy,1), &
                                 beta_air                 , &
                                 beta_soil_ps             , &
                                 veg%beta_sinklim_ps%ptr          , &
                                 beta_soil_gs             , &
                                 veg%fn_chl_cl%ptr(icanopy,1)     , &
                                 veg%fn_et_cl%ptr(icanopy,1)      , &
                                 veg%fn_rub_cl%ptr(icanopy,1)     , &
                                 veg%fn_pepc_cl%ptr(icanopy,1)    , &
                                 veg%leaf_nitrogen_cl%ptr(icanopy,1)      , &
                                 lctlib%ps_pathway                , &
                                 veg%t_jmax_opt_mavg%ptr          , &
                                 lctlib%t_jmax_omega              , &
                                 lctlib%g0                        , &
                                 lctlib%g1                        , &
                                 gross_assimilation_cl            , & ! out
                                 net_assimilation_cl              , &
                                 maint_respiration_leaf_cl               , &
                                 canopy_cond_cl                   , &
                                 co2_conc_leaf_cl                 , &
                                 hlp1, hlp1, hlp1)                   

        gross_assimilation = gross_assimilation + gross_assimilation_cl * veg%lai_cl%ptr(icanopy,1)

        IF(gross_assimilation_cl > eps8) THEN
           hlp1          = discrimination_ps(co2_conc_leaf_cl     , &
                                             co2_mixing_ratio   , &
                                             'C13'              , &
                                             lctlib%ps_pathway)
           gross_assimilation_C13 = gross_assimilation_C13  + gross_assimilation_cl * veg%lai_cl%ptr(icanopy,1) * calc_fractionation(co2_mixing_ratio  , &
                                                                                                   co2_mixing_ratio_C13 , &
                                                                                                   hlp1)
           hlp1          = discrimination_ps(co2_conc_leaf_cl     , &
                                             co2_mixing_ratio   , &
                                             'C14'              , &
                                             lctlib%ps_pathway)
           gross_assimilation_C14 = gross_assimilation_C14 + gross_assimilation_cl * veg%lai_cl%ptr(icanopy,1) * calc_fractionation(co2_mixing_ratio  , &
                                                                                                   co2_mixing_ratio_C14 , & 
                                                                                                   hlp1                 , &
                                                                                                   scale_correction_C14)
        ENDIF

        net_assimilation = net_assimilation + net_assimilation_cl * veg%lai_cl%ptr(icanopy,1) 
        maint_respiration_leaf  = maint_respiration_leaf + maint_respiration_leaf_cl * veg%lai_cl%ptr(icanopy,1) 
        canopy_cond      = canopy_cond + canopy_cond_cl * veg%lai_cl%ptr(icanopy,1)
        co2_conc_leaf    = co2_conc_leaf + co2_conc_leaf_cl * veg%lai_cl%ptr(icanopy,1)

        ! get the net photosynthesis of the lowest canopy layer at this timestep
        IF(veg%lai_cl%ptr(icanopy,1) > 0.0_wp) net_assimilation_boc = net_assimilation_cl 

     ENDIF
  ENDDO

  ! get mean canopy CO2 concentration in the leaf
  IF(veg%lai%ptr > eps8) co2_conc_leaf = co2_conc_leaf / veg%lai%ptr

  ! uptake optimum temperature for electron transport 
  IF(assimi_config%flag_t_jmax_acclimation) THEN 
     t_jmax_opt = t_jmax_offset + t_jmax_slope * ( t_air - TZero ) 
     t_jmax_opt = MIN(MAX(t_jmax_opt,t_jmax_opt_min),t_jmax_opt_max)
  ELSE 
     t_jmax_opt = lctlib%t_jmax_opt
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
                                            press_srf           , &
                                            co2_mixing_ratio    , &
                                            ga                  , &
                                            temp_acclim         , &
                                            ppfd_sunlit_cl       , &
                                            ppfd_shaded_cl       , &
                                            fleaf_sunlit_cl     , &
                                            beta_air            , &
                                            beta_soil_ps        , &
                                            beta_sinklim_ps     , &
                                            beta_soil_gs        , &
                                            fn_chl_cl          , &
                                            fn_et_cl           , &
                                            fn_rub_cl          , &
                                            fn_pepc_cl         , &
                                            leaf_nitrogen_cl               , &
                                            ps_pathway          , &
                                            t_jmax_opt          , &
                                            t_jmax_omega        , &
                                            g0                  , &
                                            g1                  , &
                                            ag_cl                 , &   ! out
                                            an_cl                 , &
                                            maint_resp_cl                  , &
                                            gs_cl                 , &
                                            ci_cl                 , &
                                            m_rub               , &
                                            m_et                , &
                                            m_pepc)                
    
  USE mo_shared_constants,    ONLY: r_gas, Tzero
  USE mo_assimi_constants
  USE mo_veg_respiration,     ONLY: temperature_response_respiration

  IMPLICIT NONE

  ! 0.1 InOut
  REAL(wp),INTENT(IN) :: temp                    !< temperature used for photosynthesis calculation (K) 
  REAL(wp),INTENT(IN) :: press_srf               !< air pressure (Pa) 
  REAL(wp),INTENT(IN) :: ga                      !< aerodynamic conductance (m s-1) 
  REAL(wp),INTENT(IN) :: temp_acclim             !< veg\%t_air_tacclim_mavg  average temperature for respiration acclimation
  REAL(wp),INTENT(IN) :: ppfd_sunlit_cl           !< PAR of sunlit foliage par per canopy layer (µmol m-2 s-1)
  REAL(wp),INTENT(IN) :: ppfd_shaded_cl           !< PAR of shaded foliage par per canopy layer (µmol m-2 s-1)
  REAL(wp),INTENT(IN) :: fleaf_sunlit_cl         !< fraction of canopy layer that is sunlit 
  REAL(wp),INTENT(IN) :: co2_mixing_ratio        !< CO2 concentration of the free air (ppm)
  REAL(wp),INTENT(IN) :: beta_air                !< scaling factor to account for air humidity effects on conductance
  REAL(wp),INTENT(IN) :: beta_soil_ps            !< scaling factor to account for soil moisture constraints on photosynthesis
  REAL(wp),INTENT(IN) :: beta_sinklim_ps         !< scaling factor to account for sink limitation constraints on photosynthesis 
  REAL(wp),INTENT(IN) :: beta_soil_gs            !< scaling factor to account for soil moisture constraints on conductance
  REAL(wp),INTENT(IN) :: fn_et_cl                !< fraction of foliar N allocated to electron transport 
  REAL(wp),INTENT(IN) :: fn_rub_cl               !< fraction of foliar N allocated to Rubisco
  REAL(wp),INTENT(IN) :: fn_pepc_cl              !< fraction of foliar N allocated to PEPC complex (C4 only)
  REAL(wp),INTENT(IN) :: fn_chl_cl               !< fraction of foliar N allocated to Chlorophyll
  REAL(wp),INTENT(IN) :: leaf_nitrogen_cl                !< total foliar N content
  INTEGER,INTENT(IN)  :: ps_pathway              !< photosynthetic pathway 
  REAL(wp),INTENT(IN) :: t_jmax_opt              !< temperature optimum of electron transport
  REAL(wp),INTENT(IN) :: t_jmax_omega            !< shape-parameter of temperature sensitivity of electron transport
  REAL(wp),INTENT(IN) :: g0                      !< offset of the stomatal conductance ~ An relationship 
  REAL(wp),INTENT(IN) :: g1                      !< slope of the stomatal conductance ~ An relationship
  REAL(wp),INTENT(OUT) :: ag_cl                  !< gross photsynthesis of layer (µmol CO2 / m2 /s) 
  REAL(wp),INTENT(OUT) :: an_cl                  !< net photsynthesis of layer (µmol CO2 / m2 /s)
  REAL(wp),INTENT(OUT) :: maint_resp_cl          !< net canopy layer respiration (µmol CO2 / m2 /s)
  REAL(wp),INTENT(OUT) :: gs_cl                  !< canopy conductance of layer (m / s) 
  REAL(wp),INTENT(OUT) :: ci_cl                  !< leaf internal CO2 concentration of layer (ppm)
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
  REAL(wp) :: an_cl_old,dan_cl        ! dummy variables needed for iteration loop of CI
  REAL(wp) :: agsf1,agsf2             ! dummy variables for PS calculation
  INTEGER :: it                       ! Iteration loop counter
  LOGICAL :: last_call                ! boolean to stop iteration at low conductance
  CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_photosynthesis'

  !> static energy (R_gas: J mol-1 K-1 -> rtck: J mol-1)
  rtck = r_gas * temp 
  !> conversion of ppm to Pa
  ppm2Pa = press_srf * 1.e-06
  !>  Net canopy layer respiration (umol[CO2]/m2/s)
  maint_resp_cl = fmaint_rate_base/1000._wp * beta_soil_ps * temperature_response_respiration(temp, temp_acclim) * leaf_nitrogen_cl


  !> @par IF: initialise output for no-light condition
  IF(ppfd_sunlit_cl + ppfd_shaded_cl <= eps8) THEN
    ag_cl    = 0.0_wp 
    an_cl    = ag_cl - maint_resp_cl
    gs_cl    = MAX(g0,g0 + g1 * an_cl * beta_air * beta_soil_gs / co2_mixing_ratio) * rtck / press_srf 
    ci_cl    = co2_mixing_ratio * ppm2Pa
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
     n1 = jmax2n * fn_et_cl * exp(-(((temp - Tzero - t_jmax_opt)/t_jmax_omega)**2)) & 
          * beta_soil_ps * beta_sinklim_ps

     !> N specific rate of Rubisco turnover 
     n2 = vcmax2n * fn_rub_cl * exp(E0v - E1v/rtck) & 
          * beta_soil_ps * beta_sinklim_ps

     !> Initial guess of leaf internal CO2 concentration
     SELECT CASE (ps_pathway)
     CASE (ic3phot)
        ci_cl = co2_mixing_ratio * CiCa_default_C3 * ppm2Pa
     CASE (ic4phot)
        ci_cl = co2_mixing_ratio * CiCa_default_C4 * ppm2Pa
     ENDSELECT

     !> iteration to find correct combination of ci_cl, Anl and gs_cl
     an_cl_old = 1000.0_wp
     dan_cl = 1.0_wp
     it = 1
     last_call = .FALSE.
     !> @par DO WHILE (it < ps_it_max .AND. dan_cl > eps4 )
     DO WHILE (it < ps_it_max .AND. dan_cl > eps4 ) 

        SELECT CASE (ps_pathway)
        CASE (ic3phot)
           !> ci sensitivity of Electron-transport limited photosynthesis
           m1 = ci_cl/(ci_cl + 2.0_wp * pcp)

           !> ci sensitivity of Rubisco-limited photosynthesis
           m2 = ci_cl/(ci_cl + km)
 
           !> N-specific light saturated carboxylation (umol[CO2]/mmol[N]/s)
           msat = MIN(m1 * n1, m2 * n2)

           !> nlim factor (part of equation 12 (Kull and Kruijt 1998))
           nco = msat/(alpha_i * ka * m1)

           !> A couple of factors to calculate photosynthesis
           !!  parts of equations 14-16 of Kull and Kruijt 1998
           agsf1 = MAX(1.0_wp - pcp / ci_cl, 0.0_wp)
           agsf2 = agsf1 * alpha_i * m1
           agsf1 = agsf1 * msat

        CASE (ic4phot)
           !> N specific rate of PePc limited photosynthessis 
           n3 = pepc2n * fn_pepc_cl * 2**((temp - TZero - Tref_pepc )/ Tbase_pepc)

           !> ci sensitivity of Electron-transport limited photosynthesis
           m1 = ci_max/(ci_max + 2.0_wp * pcp)

           !> ci sensitivity of Rubisco-limited photosynthesis
           m2 = ci_max/(ci_max + km)
 
           !> ci sensitivity of PePc-limited photosynthesis
           m3 = ci_cl / press_srf

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
        m_et = m1 * n1 / fn_et_cl
        m_rub = m2 * n2 / fn_rub_cl
        IF (fn_pepc_cl>eps8) THEN
           m_pepc = m3 * n3 / fn_pepc_cl
        ELSE
           m_pepc = 0._wp
        ENDIF

        !> fabsb: fraction of possible PAR absorption in present light conditions
        fabsb = exp(-ka * chl2n * fn_chl_cl * leaf_nitrogen_cl)

        !> @par IF: Case of no direct radiation
        IF (fleaf_sunlit_cl <= eps8) THEN
           !>  Theoretical cumulative N at co-limitation by Rubisco and
           !!  light absorption in sunlit foliage (mmol[N]/m2))
           nlim = -log(nco/(ppfd_shaded_cl * chl2n * fn_chl_cl + eps8))/(ka * chl2n * fn_chl_cl)
           nsat = MAX(MIN(nlim, leaf_nitrogen_cl), 0.0_wp)

           !> Relative PAR absorption in light-limited region (fraction)
           fabs = exp(-ka * chl2n * fn_chl_cl * nsat) - fabsb

           !> Rate of gross photosynthesis (umol/m2/s)
           !!  sum of photosynthesis in colimited, light saturated PS and PS in the 
           !!  light-limited region (assumed to scale linearily with light) 
           ag_cl = agsf1 * nsat + agsf2 * ppfd_shaded_cl * fabs

        !> @par ELSE: Case of direct and diffuse radiation
        ELSE 
          !>  Theoretical cumulative N at co-limitation by Rubisco and
          !!  light absorption in sunlit foliage (mmol[N]/m2))
          nlim = -log(nco/(ppfd_sunlit_cl * chl2n * fn_chl_cl + eps8))/(ka * chl2n * fn_chl_cl)
          nsat = MAX(MIN(nlim, leaf_nitrogen_cl), 0.0_wp)

          !>  Relative PAR absorption in light-limited region of sunlit foliage (fraction)
          fabs = exp(-ka * chl2n * fn_chl_cl * nsat) - fabsb

          !> Rate of gross photosynthesis in sunlit foliage (umol/m2/s)
          agsl = agsf1 * nsat + agsf2 * ppfd_sunlit_cl * fabs

          !> Theoretical cumulative N at co-limitation by Rubisco and
          !!  light absorption in shaded foliage (mmol[N]/m2))
          nlim = -log(nco/(ppfd_shaded_cl * chl2n * fn_chl_cl + eps8))/(ka * chl2n * fn_chl_cl)
          nsat = MAX(MIN(nlim, leaf_nitrogen_cl), 0.0_wp)

          !> Relative PAR absorption in light-limited region of shaded foliage (fraction)
          fabs = exp(-ka * chl2n * fn_chl_cl * nsat) - fabsb

          !> Rate of gross photosynthesis in shaded foliage (umol/m2/s)
          agsh = agsf1 * nsat + agsf2 * ppfd_shaded_cl* fabs

          !> Mean gross photosynthesis over sunlit and shaded fractions (umol/m2/s)
          ag_cl = fleaf_sunlit_cl * agsl + (1.0_wp - fleaf_sunlit_cl) * agsh
        END IF
        !> END IF

        !> Net layer photosynthesis (umol[CO2]/m2/s)
        an_cl = ag_cl - maint_resp_cl

        !> calculate canopy conductance (m/s) following Ball & Berry
        gs_cl = MAX(g0, g0 + g1 * an_cl * beta_air * beta_soil_gs / ( co2_mixing_ratio )) * & 
              rtck / press_srf
        !> calculate canopy conductance (m/s) following Medlyn et al. 2011
        !! gs_cl = MAX(g0, 1.6_wp * ( 1._wp + g1 * beta_air * beta_soil_gs ) * an_cl / ( co2_mixing_ratio )) * rtck / press_srf

        !> New internal CO2 from balance of supply and demand
        !! 
        !! ci = cs - An * r_tot (Friend 2005 eq. 7)
        !! Friend 2005 eq. 8: r_tot = ( Dh2oco2_air / gs_cl + Dh2oco2_turb / (speed * q_cdrag)) 
        !! Corrected to Bonan et al. 2008, Ecol. Clim: accounting for the molecular diffusion differences
        !! at the stomatal interface (1.6) and in the boundary layer, where turbulence plays a role (1.37)
        ci_cl = co2_mixing_ratio * ppm2Pa - an_cl * rtck * 1e-06 * (Dwv2co2_air/gs_cl + Dwv2co2_turb/ga)
        ci_cl = MAX (ci_cl, pcp) 

        !> change in An since last iteration and save new an_cl 
        dan_cl = abs(an_cl-an_cl_old)
        an_cl_old = an_cl

        !> break-off point for numerical stability at low conductance
        IF(last_call) THEN
           it = ps_it_max + 1
           ci_cl = pcp 
        ENDIF
!        IF( gs_cl <= g0 * rtck / press_srf .OR. ci_cl <= pcp ) THEN 
!           last_call = .TRUE.
!           ci_cl = pcp 
!        ENDIF 
        it =  it + 1
     END DO
     !> END DO WHILE
  ENDIF
  !> END IF

  !> for output purposes, convert ci to ppm
  ci_cl = ci_cl / ppm2Pa

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
     CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_atm_moisture_stress'

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
    REAL(wp),INTENT(in) :: theta       !< soil plant available water (fractions)
    REAL(wp),INTENT(in) :: lambda      !< slope parameter of Weibull function
    REAL(wp),INTENT(in) :: weibull_k   !< shape parameter of Weibull function
    REAL(wp) :: beta                   !< veg\%beta_soil_ps / veg\%beta_soil_gs
    ! 0.2 Local
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_soil_moisture_stress'

    !> 1.0 Moisture stress now as a Weibull function to avoid MIN/MAX statements
    !!
    beta = 1.0_wp - exp(-(lambda*theta)**weibull_k)
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
    ! 0.2 Local
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':discrimination_ps'

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

END MODULE mo_assimi_update_canopy_fluxes
