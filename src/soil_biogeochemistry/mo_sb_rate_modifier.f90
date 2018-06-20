!>
!! @brief Contains the routines for biophysical rate modifiers
!!
!! @author
!!  Sönke Zaehle 
!!
!! @par Revision History
!! First version by Jan Engel (20170509): copyied from mo_sb_simple_model.f90
!!
MODULE mo_sb_rate_modifier

#include "dsl_jsbach.inc"

  USE mo_kind,            ONLY: wp, eps8, eps4
  USE mo_exception,       ONLY: message, message_text, finish
  USE mo_jsb_control,     ONLY: debug_on

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: update_sb_rate_modifier


  CHARACTER(len=*), PARAMETER :: modname = 'mo_sb_rate_modifier'

CONTAINS


  !-----------------------------------------------------------------------------------------------------
  ! Main Task
  ! 
  !-----------------------------------------------------------------------------------------------------
  !> biophysical rate modifier 
  !! 
  !! 
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE update_sb_rate_modifier(tile, options)

    USE mo_jsb_class,             ONLY: Get_model
    USE mo_jsb_process_class,     ONLY: A2L_, L2A_, SEB_, TURB_, SOIL_, HYDRO_, VEG_, HD_, RAD_, ASSIMI_, PHENO_, SB_
    USE mo_jsb_tile_class,        ONLY: t_jsb_tile_abstract
    USE mo_jsb_task_class,        ONLY: t_jsb_task_options
    USE mo_jsb_model_class,       ONLY: t_jsb_model
    USE mo_jsb_lctlib_class,      ONLY: t_lctlib_element
    USE mo_sb_constants

    ! Use of process memory
    dsl4jsb_Use_memory(SB_)         !  USE mo_sb_memory_class,       ONLY: t_sb_memory
    dsl4jsb_Use_memory(SOIL_)       !  USE mo_soil_memory_class,     ONLY: t_soil_memory

    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    CLASS(t_jsb_tile_abstract), INTENT(inout)     :: tile         !< one tile with data structure for one lct
    TYPE(t_jsb_task_options),   INTENT(in)        :: options      !< model options
    ! ---------------------------
    ! 0.2 Local
    TYPE(t_jsb_model),      POINTER       :: model
    TYPE(t_lctlib_element), POINTER       :: lctlib     !< land-cover-type library - parameter across pft's
    INTEGER                               :: iblk, ics, ice, nc
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':update_sb_rate_modifier'
    ! ---------------------------
    ! 0.3 Declare Memory
    ! Declare process configuration and memory Pointers
    dsl4jsb_Def_memory(SB_)                               ! TYPE(t_sb_memory),       POINTER :: sb__mem
    dsl4jsb_Def_memory(SOIL_)                             ! TYPE(t_soil_memory),     POINTER :: soil__mem
    ! Declare pointers to variables in memory
    ! SB_
    dsl4jsb_Real3D_onChunk :: rtm_decomposition 
    dsl4jsb_Real3D_onChunk :: rtm_depolymerisation 
    dsl4jsb_Real3D_onChunk :: rtm_mic_uptake
    dsl4jsb_Real3D_onChunk :: rtm_sorption
    dsl4jsb_Real3D_onChunk :: rtm_desorption
    dsl4jsb_Real3D_onChunk :: rtm_hsc
    dsl4jsb_Real3D_onChunk :: rtm_nitrification
    dsl4jsb_Real3D_onChunk :: rtm_denitrification
    dsl4jsb_Real3D_onChunk :: rtm_gasdiffusion
    dsl4jsb_Real3D_onChunk :: rmm_decomposition 
    dsl4jsb_Real3D_onChunk :: rmm_depolymerisation 
    dsl4jsb_Real3D_onChunk :: rmm_mic_uptake
    dsl4jsb_Real3D_onChunk :: rmm_sorption
    dsl4jsb_Real3D_onChunk :: rmm_desorption
    dsl4jsb_Real3D_onChunk :: rmm_hsc
    dsl4jsb_Real3D_onChunk :: rmm_nitrification
    dsl4jsb_Real3D_onChunk :: rmm_gasdiffusion
    dsl4jsb_Real3D_onChunk :: anaerobic_volume_fraction_sl 
    ! SOIL_
    dsl4jsb_Real3D_onChunk :: t_soil_sl
    dsl4jsb_Real3D_onChunk :: w_soil_sl
    dsl4jsb_Real3D_onChunk :: w_soil_fc_sl
    dsl4jsb_Real3D_onChunk :: w_soil_sat_sl

    ! Get local variables from options argument
    iblk  = options%iblk
    ics   = options%ics
    ice   = options%ice
    nc    = options%nc
    ! ---------------------------
    ! 0.4 Process Activity, Debug Option
    IF (.NOT. tile%Is_process_active(SB_)) RETURN
    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')
    ! ---------------------------
    ! 0.5 Get Memory
    ! Get process config
!    dsl4jsb_Get_config(SB_)
    ! Get process memories
    dsl4jsb_Get_memory(SB_)
    dsl4jsb_Get_memory(SOIL_)
    ! Get process variables (Set pointers to variables in memory)
    model  => Get_model(tile%owner_model_id)
    lctlib => model%lctlib(tile%lcts%lib_id) 
    ! ---------------------------
    dsl4jsb_Get_var3D_onChunk(SB_,       rtm_decomposition)           ! out
    dsl4jsb_Get_var3D_onChunk(SB_,       rtm_depolymerisation)        ! out
    dsl4jsb_Get_var3D_onChunk(SB_,       rtm_mic_uptake)              ! out
    dsl4jsb_Get_var3D_onChunk(SB_,       rtm_sorption)                ! out
    dsl4jsb_Get_var3D_onChunk(SB_,       rtm_desorption)              ! out
    dsl4jsb_Get_var3D_onChunk(SB_,       rtm_hsc)                     ! out
    dsl4jsb_Get_var3D_onChunk(SB_,       rtm_nitrification)           ! out
    dsl4jsb_Get_var3D_onChunk(SB_,       rtm_denitrification)         ! out
    dsl4jsb_Get_var3D_onChunk(SB_,       rtm_gasdiffusion)            ! out
    dsl4jsb_Get_var3D_onChunk(SB_,       rmm_decomposition)           ! out
    dsl4jsb_Get_var3D_onChunk(SB_,       rmm_depolymerisation)        ! out
    dsl4jsb_Get_var3D_onChunk(SB_,       rmm_mic_uptake)              ! out
    dsl4jsb_Get_var3D_onChunk(SB_,       rmm_sorption)                ! out
    dsl4jsb_Get_var3D_onChunk(SB_,       rmm_desorption)              ! out
    dsl4jsb_Get_var3D_onChunk(SB_,       rmm_hsc)                     ! out
    dsl4jsb_Get_var3D_onChunk(SB_,       rmm_nitrification)           ! out
    dsl4jsb_Get_var3D_onChunk(SB_,       rmm_gasdiffusion)            ! out
    dsl4jsb_Get_var3D_onChunk(SB_,       anaerobic_volume_fraction_sl)! out 
    ! ---------------------------
    dsl4jsb_Get_var3D_onChunk(SOIL_,     t_soil_sl)                   ! in
    dsl4jsb_Get_var3D_onChunk(SOIL_,     w_soil_sl)                   ! in
    dsl4jsb_Get_var3D_onChunk(SOIL_,     w_soil_fc_sl)                ! in
    dsl4jsb_Get_var3D_onChunk(SOIL_,     w_soil_sat_sl)               ! in



    ! ------------------------------------------------------------------------------------------------------------
    ! Go science
    ! ------------------------------------------------------------------------------------------------------------


    !> 1.0 calc modifiers

    rtm_depolymerisation = calc_temperature_rate_modifier(t_soil_sl,ea_depolymerisation) 
    rtm_mic_uptake       = calc_temperature_rate_modifier(t_soil_sl,ea_mic_uptake)
    rtm_sorption         = calc_temperature_rate_modifier(t_soil_sl,ea_sorption) 
    rtm_desorption       = calc_temperature_rate_modifier(t_soil_sl,ea_desorption) 
    rtm_hsc              = calc_temperature_rate_modifier(t_soil_sl,ea_hsc) 
     
    rmm_depolymerisation = calc_oxygen_rate_modifier(w_soil_sl,w_soil_fc_sl,kmr_depolymerisation) 
    rmm_mic_uptake       = calc_oxygen_rate_modifier(w_soil_sl,w_soil_fc_sl,kmr_mic_uptake)
    rmm_sorption         = calc_convert_bulk_to_liquid(w_soil_sl)
    rmm_desorption       = calc_convert_bulk_to_liquid(w_soil_sl) 
    rmm_hsc              = calc_substrate_rate_modifier(w_soil_sl,w_soil_fc_sl,kmr_hsc) 

    !> 2.0 simple soil model 
    !!
    rtm_decomposition    = calc_temperature_rate_modifier(t_soil_sl,ea_depolymerisation) 
    rmm_decomposition    = calc_oxygen_rate_modifier(w_soil_sl,w_soil_fc_sl,kmr_depolymerisation) * &
                           calc_moisture_modifier(w_soil_sl,w_soil_fc_sl)

    !> 3.0 nitrification-denitrification model
    !!
    rtm_nitrification    = calc_peaked_temperature_rate_modifier(t_soil_sl,ea_nitrification,ed_nitrification,t_opt_nitrification) 
    rtm_denitrification  = calc_temperature_rate_modifier(t_soil_sl,ea_denitrification) 
    rtm_gasdiffusion     = calc_temperature_rate_modifier(t_soil_sl,ea_gasdiffusion) 
    rmm_nitrification    = calc_moisture_modifier(w_soil_sl,w_soil_fc_sl) 
    rmm_gasdiffusion     = 1.0_wp - calc_moisture_modifier(w_soil_sl,w_soil_fc_sl) 
    anaerobic_volume_fraction_sl = calc_anaerobic_volume_fraction(w_soil_sl,w_soil_sat_sl,lambda_anvf,k_anvf) 

  END SUBROUTINE update_sb_rate_modifier


  !-----------------------------------------------------------------------------------------------------
  ! Functionality used by the SB module, derived from JSM 
  !
  !-----------------------------------------------------------------------------------------------------
  !> Calculates temperature response factor for various soil biogeochemical reactions
  !!  according to the Arrhenius equation as defined in
  !!  Wang, G., Post, W.M., Mayes, M.A., Frerichs, J.T., Sindhu, J., 2012.
  !!  Parameter estimation for models of ligninolytic and cellulolytic enzyme kinetics.
  !!  Soil Biology and Biochemistry 48, 28-38.
  !!  NB if sensitivity of zero is specified, 1 will be returned for all layers
  !!
  !! Soil temperatures (t_soil_sl) lower eps8 are caught by using eps8 instead of t_soil_sl for calculations, 
  !!  that is, t_soil_sl is NOT modified in case: t_soil_sl < eps8
  !! 
  !! Input: soil temperature, activiation energy 
  !!
  !! Output: rate modifier (unitless) 
  !-----------------------------------------------------------------------------------------------------
  ELEMENTAL FUNCTION calc_temperature_rate_modifier(t_soil_sl,Ea) RESULT(rate_modifier)

    USE mo_shared_constants,    ONLY: r_gas
    USE mo_sb_constants,        ONLY: t_ref_tresponse 

    IMPLICIT NONE
    ! ---------------------------
    !0.1 InOut
    REAL(wp), INTENT(in)  :: t_soil_sl     !< soil temperature [K]
    REAL(wp), INTENT(in)  :: Ea            !< activation energy [J mol-1]
    REAL(wp)              :: rate_modifier !< temperature rate modifier [unitless] 
    ! ---------------------------
    ! 0.2 Local
    REAL(wp)              :: hlp           ! either t_soil_sl or eps8
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_temperature_rate_modifier'

    ! temperature must not be negative in the below function
    IF(t_soil_sl < eps8) THEN
      hlp = eps8
    ELSE
      hlp = t_soil_sl
    END IF

    IF(Ea > eps8)THEN
      rate_modifier = exp(-Ea/r_gas*(1._wp/hlp - 1._wp/t_ref_tresponse))
    ELSE
      rate_modifier = 1._wp
    ENDIF   

  END FUNCTION calc_temperature_rate_modifier

  !-----------------------------------------------------------------------------------------------------
  ! Functionality used by the SB module
  !
  !-----------------------------------------------------------------------------------------------------
  !> Calculates temperature response factor for various soil biogeochemical reactions
  !!  according to the peaked Arrhenius equation as defined in
  !!  Medlyn et al. 2002, PCE, eq. 18
  !!  NB if sensitivity of zero is specified, 1 will be returned for all layers
  !!
  !! Input: soil temperature, (de-)activiation energy, optimum temperature
  !!
  !! Output: rate modifier (unitless) 
  !-----------------------------------------------------------------------------------------------------
  ELEMENTAL FUNCTION calc_peaked_temperature_rate_modifier(t_soil_sl,Ea,Ed,t_opt) RESULT(rate_modifier)

    USE mo_shared_constants,    ONLY: r_gas
    USE mo_sb_constants,        ONLY: t_ref_tresponse 

    IMPLICIT NONE
    ! ---------------------------
    !0.1 InOut
    REAL(wp), INTENT(in)  :: t_soil_sl     !< soil temperature [K]
    REAL(wp), INTENT(in)  :: Ea            !< activation energy [J mol-1]
    REAL(wp), INTENT(in)  :: Ed            !< deactivation energy [J mol-1]
    REAL(wp), INTENT(in)  :: t_opt         !< optimum temperature [K] 
    REAL(wp)              :: rate_modifier !< temperature rate modifier [unitless] 
    ! ---------------------------
    ! 0.2 Local
    REAL(wp)              :: hlp1, hlp2
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_temperature_rate_modifier'

    IF(Ea>eps8)THEN
       hlp1 = t_soil_sl - t_opt
       hlp2 = t_soil_sl * t_opt * r_gas
       rate_modifier = (Ed * exp(Ea * hlp1 / hlp2)) / (Ed - Ea * (1._wp - exp(Ed * hlp1 / hlp2))) 
    ELSE
       rate_modifier = 1._wp
    ENDIF   

  END FUNCTION calc_peaked_temperature_rate_modifier


  !-----------------------------------------------------------------------------------------------------
  ! Functionality used by the SB module, derived from JSM 
  !
  !-----------------------------------------------------------------------------------------------------
  !> Oxygen limitation of decomposition and uptake due to soil moisture response function 
  !!  according to Davidson et al's DAMM model (2012)
  !!  NB if sensitivity of zero is specified, 1 will be returned for all layers
  !!
  !! Input: soil moisture, sensitivity parameter 
  !!
  !! Output: rate modifier (unitless) 
  !-----------------------------------------------------------------------------------------------------
  ELEMENTAL FUNCTION calc_oxygen_rate_modifier(w_soil_sl,w_soil_fc_sl,param) RESULT(rate_modifier)

    USE mo_sb_constants,        ONLY: k1_afps 

    IMPLICIT NONE
    ! ---------------------------
    !0.1 InOut
    REAL(wp), INTENT(in)  :: w_soil_sl     !< soil moisture [m]
    REAL(wp), INTENT(in)  :: w_soil_fc_sl  !< soil moisture at field capacity [m]
    REAL(wp), INTENT(in)  :: param         !< 
    REAL(wp)              :: rate_modifier !< oxygen availability rate modifier [unitless] 
    ! ---------------------------
    ! 0.2 Local
    REAL(wp)              :: air_filled_pore_space
    REAL(wp)              :: hlp
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_oxygen_rate_modifier'

    IF(param>eps8) THEN
       IF(w_soil_fc_sl > eps8) THEN
          air_filled_pore_space = (w_soil_fc_sl - w_soil_sl) / w_soil_fc_sl
          IF(air_filled_pore_space < 0.005_wp) air_filled_pore_space = 0.005_wp
          rate_modifier = (air_filled_pore_space**k1_afps) / (param + (air_filled_pore_space**k1_afps))
       ELSE
          rate_modifier = 1._wp
       ENDIF
    ELSE 
       rate_modifier = 1._wp
    ENDIF

  END FUNCTION calc_oxygen_rate_modifier


  !-----------------------------------------------------------------------------------------------------
  ! Functionality used by the SB module, derived from JSM 
  !
  !-----------------------------------------------------------------------------------------------------
  !> Converts bulk density to liquid concentration 
  !!  NB if sensitivity of zero is specified, 1 will be returned for all layers
  !!
  !! Input: soil moisture, sensitivity parameter 
  !!
  !! Output: conversion modifier (unitless) 
  !-----------------------------------------------------------------------------------------------------
  ELEMENTAL FUNCTION calc_convert_bulk_to_liquid(w_soil_sl) RESULT(factor)

    IMPLICIT NONE
    ! ---------------------------
    !0.1 InOut
    REAL(wp), INTENT(in)  :: w_soil_sl           !< soil moisture [m]
    REAL(wp)              :: factor              !< conversion factor [unitless] 
    ! ---------------------------
    ! 0.2 Local
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_convert_bulk_to_liquid'

    IF(w_soil_sl > eps4) THEN
      factor = 1._wp / w_soil_sl 
    ELSE
       factor = 1._wp / eps4 
    ENDIF

  END FUNCTION calc_convert_bulk_to_liquid


  !-----------------------------------------------------------------------------------------------------
  ! Functionality used by the SB module, derived from JSM 
  !
  !-----------------------------------------------------------------------------------------------------
  !> Substrate limitation of decomposition and uptake due to soil moisture response function 
  !!  according to Davidson et al's DAMM model (2012)
  !!  NB if sensitivity of zero is specified, 1 will be returned for all layers
  !!
  !! Input: soil moisture, sensitivity parameter 
  !!
  !! Output: rate modifier (unitless) 
  !-----------------------------------------------------------------------------------------------------
  ELEMENTAL FUNCTION calc_substrate_rate_modifier(w_soil_sl,w_soil_fc_sl,param) RESULT(rate_modifier)

    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    REAL(wp), INTENT(in)  :: w_soil_sl           !< soil moisture [m]
    REAL(wp), INTENT(in)  :: w_soil_fc_sl        !< soil moisture at field capacity [m]
    REAL(wp), INTENT(in)  :: param               !< 
    REAL(wp)              :: rate_modifier !< oxygen availability rate modifier [unitless] 
    ! ---------------------------
    ! 0.2 Local
    REAL(wp)              :: hlp1 
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_substrate_rate_modifier'

    IF(param>eps8.AND.w_soil_fc_sl > eps8) THEN
       hlp1 = (w_soil_sl / w_soil_fc_sl)**3._wp
       IF(hlp1 > param) THEN
          rate_modifier = 1._wp / hlp1 
       ELSE 
          rate_modifier = 1._wp / param 
       ENDIF
    ELSE
      rate_modifier = 1._wp
    ENDIF

  END FUNCTION calc_substrate_rate_modifier

  !-----------------------------------------------------------------------------------------------------
  ! Functionality used by the SB module
  !
  !-----------------------------------------------------------------------------------------------------
  !> Anaerobic fraction of the soil layer  
  !!  according to Zaehle et al. 2011 
  !!  NB if sensitivity of zero is specified, 1 will be returned for all layers
  !!
  !! Input: soil moisture, sensitivity parameters 
  !!
  !! Output: rate modifier (unitless) 
  !-----------------------------------------------------------------------------------------------------
  ELEMENTAL FUNCTION calc_anaerobic_volume_fraction(w_soil_sl,w_soil_sat_sl,lambda,k) RESULT(rate_modifier)

    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    REAL(wp), INTENT(in)  :: w_soil_sl           !< soil moisture [m]
    REAL(wp), INTENT(in)  :: w_soil_sat_sl       !< soil moisture at saturation [m]
    REAL(wp), INTENT(in)  :: lambda              !< lambda in Weibull function 
    REAL(wp), INTENT(in)  :: k                   !< k in Weibull function
    REAL(wp)              :: rate_modifier       !< oxygen availability rate modifier [unitless] 
    ! ---------------------------
    ! 0.2 Local
    REAL(wp)              :: air_filled_pore_space
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_substrate_rate_modifier'

    air_filled_pore_space = (w_soil_sat_sl - w_soil_sl ) / w_soil_sat_sl
    IF(air_filled_pore_space > eps4) THEN
       rate_modifier = 1.0_wp - exp(-(lambda * (1.0_wp - air_filled_pore_space))**k) 
    ELSE 
       rate_modifier = 1._wp - exp(-(lambda * (1.0_wp - eps4))**k)
    ENDIF

  END FUNCTION calc_anaerobic_volume_fraction

  !-----------------------------------------------------------------------------------------------------
  ! Functionality used by the SB module
  !
  !-----------------------------------------------------------------------------------------------------
  !> soil moisture response function for the nitrification model ((FC-soil_moisture)/FC) 
  !!
  !! Input: soil moisture, sensitivity parameter 
  !!
  !! Output: rate modifier (unitless) 
  !-----------------------------------------------------------------------------------------------------
  ELEMENTAL FUNCTION calc_moisture_modifier(w_soil_sl,w_soil_fc_sl) RESULT(rate_modifier)

    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    REAL(wp), INTENT(in)  :: w_soil_sl           !< soil moisture [m]
    REAL(wp), INTENT(in)  :: w_soil_fc_sl        !< soil moisture at field capacity [m]
    REAL(wp)              :: rate_modifier       !< rate modifier [unitless] 
    ! ---------------------------
    ! 0.2 Local
    REAL(wp)              :: air_filled_pore_space
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_moisture_modifier'

    IF(w_soil_fc_sl > eps8) THEN
       air_filled_pore_space = (w_soil_fc_sl - w_soil_sl) / w_soil_fc_sl
       rate_modifier = 1.0_wp - air_filled_pore_space
    ELSE
       rate_modifier = 1._wp
    ENDIF

  END FUNCTION calc_moisture_modifier


END MODULE mo_sb_rate_modifier
