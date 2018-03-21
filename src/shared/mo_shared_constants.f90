!! @brief Contains general parameters required for land calculations 
!!
!! @author
!!  Sönke Zaehle, MPI-BGC Jena 
!!
!! @par Revision History
!! First version by Sönke Zaehle (2015-12-04)
!!
MODULE mo_shared_constants
 
  USE mo_kind,   ONLY: wp

  IMPLICIT NONE
  PUBLIC

  ! physical constants
  REAL(wp), SAVE :: &
    r_gas        ,  &           !< universal gas constant (J / mol / K) 
    r_gas_dryair ,  &           !< specific gas constant for dry air (J / kg / K) 
    r_gas_vapour ,  &           !< specific gas constant for water vapour (J / kg / K) 
    Dwv          ,  &           !< diffusion coefficient for water vapour in air (m2 s-1) 
    Dco2         ,  &           !< diffusion coefficient for CO2 in air (m2 s-1) 
    Tzero        ,  &           !< melting point in K -> conversion factor for degreeC
!ECHAM
    grav         ,  &           !< gravitational acceleration ( m / s2 )
    rhoh2o       ,  &           !< density of liquid water (kg/m3)
    stbo         ,  &           !< Stephan-Boltzmann constant             [W/(m2 K4)]
    zemiss_def   ,  &           !< Long-wave surface emissivity factor    []
    cpd          ,  &           !< dry air specific heat at constant pressure [J/K/kg] 
    cpv          ,  &           !< water vapour specific heat at constant pressure  [J/K/kg] 
    cpvd1        ,  &           !< auxiliary constants
    cvdifts      ,  &           !< ???
    ckap         ,  &           !< von Karmann constant for evaporation
    LatentHeatVaporization, &   !< Latent heat for vapourisation (J/kg) 
    LatentHeatSublimation       !< Latent heat for sublimation (J/kg) 
  ! molar masses
  REAL(wp), SAVE :: &
    molar_mass_C    ,  &        !< g C / mol C
    molar_mass_C12  ,  &        !< g C / mol C
    molar_mass_C13  ,  &        !< g C / mol C
    molar_mass_C14  ,  &        !< g C / mol C
    molar_mass_N    ,  &        !< g N / mol N
    molar_mass_N14  ,  &        !< g N / mol N
    molar_mass_N15  ,  &        !< g N / mol N
    molar_mass_P                !< g P / mol P
  ! isotope standards
  REAL(wp), SAVE ::            &
    PDB_STANDARD_C13           , &  !< 13C/12C PDB standard
    AAB_STANDARD_C14           , &  !< 0.95·0.238 Bq/gC is the absolute specific activity of 
                                    !! the radiocarbon standard (i.e. 95% of the activity of the OxA-I standard)
    ATM_STANDARD_N15           , &  !< 15N/14N atmospheric standard
    AvogadroC                  , &  !< Avogadro Number
    lambda_C14                 , &  !< decay constant of C14 (1/s)
    AvogadroC_times_lambda_C14 , &  !< Avogadro Number * decay constant of C14 to reduce the use of the very large/small AvogadroC/lambda_C14 
    scale_correction_c14            !< conversion from mol to pico-mol for C14 (for numerical precision) 


  !!----------------------------------------------------------------------------------------------------
  !! needed for mo_shared_update_time_averages:update_time_averages()
  !! 
  !! moving average period lengths (you can add new parameters here)
  !! the period may be defined in days per default (but different units are possible, see calc_mavg() below)
                                                                 ! generic memory time-scales:
    REAL(wp),SAVE :: mavg_period_weekly            = 7.0_wp        !< 7 days; constant across the model
    REAL(wp),SAVE :: mavg_period_monthly           = 30.0_wp       !< 30 days; constant across the model
    REAL(wp),SAVE :: mavg_period_yearly            = 365.0_wp      !< 365 days
                                                                 ! specific memory time-scales of/for:
    REAL(wp),SAVE :: mavg_period_tsoa              = 2.0_wp        !< memory time_scale for state of acclimation
    REAL(wp),SAVE :: mavg_period_tuptake           = 3.0_wp        !< nutrient uptake processes (demand, BNF) 
    REAL(wp),SAVE :: mavg_period_tlabile           = 7.0_wp        !< labile pool 
    REAL(wp),SAVE :: mavg_period_tphen             = 7.0_wp        !< phenological processes 
    REAL(wp),SAVE :: mavg_period_tfrac             = 10.0_wp       !< leaf N fractions
    REAL(wp),SAVE :: mavg_period_tcnl              = 20.0_wp       !< leaf stoichiometry 
    REAL(wp),SAVE :: mavg_period_talloc            = 20.0_wp       !< biomass allocation
    REAL(wp),SAVE :: mavg_period_tacclim           = 30.0_wp       !< respiration acclimation
    REAL(wp),SAVE :: mavg_period_tmic              = 30.0_wp       !< microbial community acclimation
    REAL(wp),SAVE :: mavg_period_tgrowth           = 365.0_wp      !< vegetation growth (?)
    REAL(wp),SAVE :: mavg_period_tvegdyn           = 365.0_wp      !< vegetation dynamics processes 

!    REAL(wp),SAVE :: mavg_period_xxx         = xxx.0_wp       !< xxx
!    REAL(wp),SAVE :: mavg_period_xxx         = xxx.0_wp       !< xxx
!    REAL(wp),SAVE :: mavg_period_xxx         = xxx.0_wp       !< xxx
  !!----------------------------------------------------------------------------------------------------

  CHARACTER(len=*), PARAMETER, PRIVATE :: modname = 'mo_shared_constants'

CONTAINS

  !-----------------------------------------------------------------------------------------------------
  !> initialize constants from shared
  !! 
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE init_shared_constants

    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    ! 
    ! ---------------------------
    ! 0.2 Local
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':init_shared_constants'


    ! physical constants
    r_gas        = 8.3144_wp    
    r_gas_dryair = 287.058_wp   
    r_gas_vapour = 461.51_wp   
    Dwv          = 0.0000242_wp 
    Dco2         = 0.0000151_wp 
    Tzero        = 273.15_wp
    grav         = 9.81_wp
    rhoh2o       = 1000._wp 
    stbo         = 5.6704E-8_wp
    zemiss_def   = 0.996_wp
    cpd          = 1004.64_wp       
    cpv          = 1869.46_wp   
    cpvd1        = cpv/cpd-1.0_wp   
    cvdifts      = 1.0_wp
    ckap         = 0.41_wp
    LatentHeatVaporization = 2.5008e6_wp
    LatentHeatSublimation  = 2.8345e6_wp
    ! molar masses
    molar_mass_C   = 12.0112_wp 
    molar_mass_C12 = 12.0_wp    
    molar_mass_C13 = 13.0_wp    
    molar_mass_C14 = 14.0_wp    
    molar_mass_N   = 14.0072_wp 
    molar_mass_N14 = 14.0_wp    
    molar_mass_N15 = 15.0_wp    
    molar_mass_P   = 30.997_wp  
    ! isotope standards
    PDB_STANDARD_C13            = 0.0112372_wp   
    AAB_STANDARD_C14            = 0.2261_wp      
    ATM_STANDARD_N15            = 0.0036765_wp !Robinson 2001, TREE
    AvogadroC                   = 6.022_wp*1.e23_wp     
    lambda_C14                  = 3.8332_wp*1.e-12_wp    
    AvogadroC_times_lambda_C14  = 6.022_wp * 3.8332_wp * 1.e11_wp  ! used in e.g. calc_delta_C14() & calc_mixing_ratio_C14C()
    scale_correction_c14        = 1.e12_wp                         !1000000000000._wp   

  END SUBROUTINE init_shared_constants

END MODULE mo_shared_constants
