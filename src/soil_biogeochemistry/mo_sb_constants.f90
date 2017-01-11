!>
!! @brief Contains constants for the soil_biogeochemistry process
!!
!! @author
!!  Jan Engel, MPI-BGC Jena
!!
!! @par Revision History
!! First version by Jan Engel (20170104)
!!
MODULE mo_sb_constants

  USE mo_kind, ONLY: wp

  IMPLICIT NONE
  PUBLIC

  ! parameters for simple soil biogeochemical model
  REAL(wp), SAVE :: & 
      tau_litter, &        !< soil litter pool turnover time (years)
      floss_nmin, &        !< fractional loss of N per unit N mineralisation 
      fnit_nmin            !< fractional nitrification per unit N mineralisation 

  ! parameter for JSM implementation
  REAL(wp), SAVE :: &
      t_ref_tresponse, &      !< reference temperature for temperature response function [K]
      ea_depolymerisation, &  !< activation energy for depolymerisation [J mol-1]
      ea_mic_uptake, &        !< activation energy for microbial uptake [J mol-1]
      ea_sorption, &          !< activation energy for sorption to mineral surfaces [J mol-1]
      ea_desorption, &        !< activation energy for desorption to mineral surfaces [J mol-1]
      ea_hsc, &               !< activation energy for half-saturation constants [J mol-1]
      k1_afps, &              !< scaling factor in the air-filled pore space to modifier response
      kmr_depolymerisation, & !< scaling factor for the sensitivity of depolymerisation to oxygen limitation 
      kmr_mic_uptake, &       !< scaling factor for the sensitivity of microbial uptake to oxygen limitation 
      kmr_sorption, &         !< scaling factor for the sensitivity of sorption to moisture limitation 
      kmr_desorption, &       !< scaling factor for the sensitivity of desorption to moisture limitation 
      kmr_hsc                 !< scaling factor for the sensitivity of half-saturation constant to moisture limitation 
                             
  

  CHARACTER(len=*), PARAMETER, PRIVATE :: modname = 'mo_sb_constants'

CONTAINS

  !-----------------------------------------------------------------------------------------------------
  !> initialize constants from soil_biogeochemistry
  !! 
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE init_sb_constants

    IMPLICIT NONE
    ! 0.2 Local
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':init_sb_constants'

    ! parameters for simple soil biogeochemical model 
    tau_litter = 5.0_wp  ! years
    floss_nmin = 0.05_wp
    fnit_nmin  = 0.1_wp

    ! parameter for JSM implementation
    t_ref_tresponse = 293.15_wp   ! Wang, G., et al. 2012., Soil Biology and Biochemistry 48, 28-38.
    ea_depolymerisation = 53000._wp   ! B. Ahrens, pers. comm. 2016
    ea_mic_uptake       = 47000._wp   ! B. Ahrens, pers. comm. 2016
    ea_sorption         = 5000._wp   ! B. Ahrens, pers. comm. 2016
    ea_desorption       = 20000._wp   ! B. Ahrens, pers. comm. 2016
    ea_hsc              = 30000._wp   ! B. Ahrens, pers. comm. 2016
    k1_afps             = 4._wp / 3._wp ! B. Ahrens, pers. comm. 2016
    kmr_depolymerisation= 0.001_wp
    kmr_mic_uptake      = 0.001_wp
    kmr_sorption        = 0.50853_wp
    kmr_desorption      = 0.50853_wp
    kmr_hsc             = 0.0_wp


  END SUBROUTINE init_sb_constants

END MODULE mo_sb_constants
