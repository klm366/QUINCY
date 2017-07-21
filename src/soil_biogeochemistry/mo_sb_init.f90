!>
!! @brief Initialization of the the soil_biogeochemistry memory.
!!
!! @author
!!  Sönke Zaehle, MPI-BGC Jena
!!
!! @par Revision History
!! First version for SLM by Jan Engel (20161208) based on the JBACH4 version by Reiner Schnur (2013-04-04)
!!
MODULE mo_sb_init

#include "dsl_jsbach.inc"

  USE mo_kind,                  ONLY: wp
  USE mo_exception,             ONLY: finish,message

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: sb_init

  CHARACTER(len=*), PARAMETER :: modname = 'mo_sb_init'

CONTAINS

  !-----------------------------------------------------------------------------------------------------
  !> Intialize soil_biogeochemistry process (after memory has been set up)
  !! 
  !! 
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE sb_init(tile, options)

    USE mo_jsb_class,             ONLY: Get_model
    USE mo_base_memory,           ONLY: t_slm_tile, t_jsb_task_options
    USE mo_jsb_model_class,       ONLY: t_jsb_model
    USE mo_base_config_class,     ONLY: base_config  
    USE mo_jsb_lctlib_class,      ONLY: t_lctlib_element
    USE mo_jsb_process_class,     ONLY: A2L_, L2A_, SEB_, TURB_, SOIL_, HYDRO_, VEG_, HD_, RAD_, ASSIMI_, PHENO_, SB_
    USE mo_shared_isotope,        ONLY: calc_mixing_ratio_C13C12, calc_mixing_ratio_C14C, calc_mixing_ratio_N15N14
    USE mo_shared_constants,      ONLY: molar_mass_C, molar_mass_N, molar_mass_P
    USE mo_sb_jsm_processes
    USE mo_sb_jsm_transport
    USE mo_veg_constants,         ONLY: itree
    USE mo_sb_constants

    ! Use of process configurations (t_PROC_config)
    dsl4jsb_Use_config(SB_)
    dsl4jsb_Use_config(SOIL_)

    ! USE PROC_memory
    dsl4jsb_Use_memory(SB_)
    dsl4jsb_Use_memory(SOIL_)

    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    TYPE(t_slm_tile),         INTENT(inout)      :: tile       !< data structure for dynamic vegetation and soil variables
    TYPE(t_jsb_task_options), INTENT(in)         :: options      !< model options
    ! ---------------------------
    ! 0.2 Local
    TYPE(t_jsb_model),      POINTER              :: model
    TYPE(t_lctlib_element), POINTER              :: lctlib     !< land-cover-type library - parameter across pft's
    INTEGER                                      :: isoil
    INTEGER                                      :: iblk, ics, ice, nc
    REAL(wp)                                     :: hlp1,hlp2,hlp3,hlp4,hlp5,hlp6
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':sb_init'
    ! ---------------------------
    ! 0.3 Memory allocation
    ! Declare process configuration and memory Pointers
    dsl4jsb_Def_config(SB_)
    dsl4jsb_Def_config(SOIL_)
    dsl4jsb_Def_memory(SB_)
    dsl4jsb_Def_memory(SOIL_)
!! *slmtemp    
    TYPE(t_sb_memory),    POINTER                :: sb         ! local pointer to sb__mem for SLM usability
    TYPE(t_soil_memory),    POINTER              :: soil       ! local pointer to soil__mem for SLM usability
    ! Get local variables from options argument
    iblk    = options%iblk
    ics     = options%ics
    ice     = options%ice
    nc      = options%nc
    ! ---------------------------
    ! 0.4 Process Activity, Debug Option
!    IF (.NOT. tile%Is_process_active(xx_)) RETURN
!    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')
    ! ---------------------------
    ! 0.5 Get Memory
    model  => Get_model(tile%owner_model_id)
    lctlib => model%lctlib(tile%lcts%lib_id)
    ! Get process config
    dsl4jsb_Get_config(SB_)
    dsl4jsb_Get_config(SOIL_)
    ! Get process memories
    dsl4jsb_Get_memory(SB_)
    dsl4jsb_Get_memory(SOIL_)
    ! Get process variables (Set pointers to variables in memory)
    ! ...


!! *slmtemp
    sb      => sb__mem 
    soil    => soil__mem 


    ! ------------------------------------------------------------------------------------------------------------
    ! Go science
    ! ------------------------------------------------------------------------------------------------------------

    SELECT CASE(TRIM(dsl4jsb_Config(SB_)%sb_model_scheme))

       ! Initial values for Soenke's four pool soil biogeochemical model
       CASE("simple_0d", "simple_1d")

          hlp1 = 100._wp ! g C m-2 metabolic litter
          hlp2 = 100._wp ! g C m-2 structural litter
          hlp3 = 1000._wp ! g C m-2 woody litter
          hlp4 = 4500._wp ! g C m-2 fast SOM
          hlp5 = 4500._wp ! g C m-2 slow SOM
          hlp6 = 0.6_wp   ! fraction of production not roots 

          IF(TRIM(dsl4jsb_Config(SB_)%sb_model_scheme) == "simple_0d") THEN
             sb%pool%soluable_litter%carbon%ptr(ics:ice,1,iblk) = hlp1 / molar_mass_C / dsl4jsb_Config(SOIL_)%dz_water(1)
             sb%pool%polymeric_litter%carbon%ptr(ics:ice,1,iblk) = hlp2 / molar_mass_C / dsl4jsb_Config(SOIL_)%dz_water(1)
             IF(lctlib%growthform==itree) &
                sb%pool%woody_litter%carbon%ptr(ics:ice,1,iblk) = hlp3 / molar_mass_C / dsl4jsb_Config(SOIL_)%dz_water(1)
             sb%pool%microbial%carbon%ptr(ics:ice,1,iblk) = hlp4 / molar_mass_C / dsl4jsb_Config(SOIL_)%dz_water(1)
             sb%pool%residue%carbon%ptr(ics:ice,1,iblk) = hlp5 / molar_mass_C / dsl4jsb_Config(SOIL_)%dz_water(1)

             sb%nh4_solute%ptr(ics:ice,1,iblk) = 0.1_wp / dsl4jsb_Config(SOIL_)%dz_water(1)
             sb%no3_solute%ptr(ics:ice,1,iblk) = 0.1_wp / dsl4jsb_Config(SOIL_)%dz_water(1)
             sb%po4_solute%ptr(ics:ice,1,iblk) = 0.01_wp / dsl4jsb_Config(SOIL_)%dz_water(1) 
          ELSE 
             sb%pool%soluable_litter%carbon%ptr(ics:ice,1,iblk)  = hlp1 /  molar_mass_C / dsl4jsb_Config(SOIL_)%dz_water(1) * &
                                                                       ( hlp6 + (1.0_wp - hlp6) * soil%root_fraction_sl%ptr(ics:ice,1,iblk)) 
             sb%pool%polymeric_litter%carbon%ptr(ics:ice,1,iblk) = hlp2 / molar_mass_C / dsl4jsb_Config(SOIL_)%dz_water(1) * &
                                                                       ( hlp6 + (1.0_wp - hlp6) * soil%root_fraction_sl%ptr(ics:ice,1,iblk)) 
             IF(lctlib%growthform==itree) &
                sb%pool%woody_litter%carbon%ptr(ics:ice,1,iblk)  = hlp3 / molar_mass_C / dsl4jsb_Config(SOIL_)%dz_water(1)
             sb%pool%microbial%carbon%ptr(ics:ice,1,iblk)        = hlp4 / molar_mass_C / dsl4jsb_Config(SOIL_)%dz_water(1) * &
                                                                       ( hlp6 + (1.0_wp - hlp6) * soil%root_fraction_sl%ptr(ics:ice,1,iblk)) 
             sb%pool%residue%carbon%ptr(ics:ice,1,iblk)          = hlp5 / molar_mass_C / dsl4jsb_Config(SOIL_)%dz_water(1) * & 
                                                                       ( hlp6 + (1.0_wp - hlp6) * soil%root_fraction_sl%ptr(ics:ice,1,iblk)) 
             DO isoil=2,options%nsoil_sb
                sb%pool%soluable_litter%carbon%ptr(ics:ice,isoil,iblk)  = hlp1 /  molar_mass_C / dsl4jsb_Config(SOIL_)%dz_water(isoil) * &
                                                                          (1.0_wp - hlp6) * soil%root_fraction_sl%ptr(ics:ice,isoil,iblk) 
                sb%pool%polymeric_litter%carbon%ptr(ics:ice,isoil,iblk) = hlp2 / molar_mass_C / dsl4jsb_Config(SOIL_)%dz_water(isoil) * &
                                                                          (1.0_wp - hlp6) * soil%root_fraction_sl%ptr(ics:ice,isoil,iblk) 
                sb%pool%microbial%carbon%ptr(ics:ice,isoil,iblk)        = hlp4 / molar_mass_C / dsl4jsb_Config(SOIL_)%dz_water(isoil) * &
                                                                          (1.0_wp - hlp6) * soil%root_fraction_sl%ptr(ics:ice,isoil,iblk) 
                sb%pool%residue%carbon%ptr(ics:ice,isoil,iblk)          = hlp5 / molar_mass_C / dsl4jsb_Config(SOIL_)%dz_water(isoil) * & 
                                                                          (1.0_wp - hlp6) * soil%root_fraction_sl%ptr(ics:ice,isoil,iblk)

                sb%nh4_solute%ptr(ics:ice,isoil,iblk) = 0.1_wp
                sb%no3_solute%ptr(ics:ice,isoil,iblk) = 0.1_wp 
                sb%po4_solute%ptr(ics:ice,isoil,iblk) = 0.01_wp 
             ENDDO
          ENDIF
 
          ! metabolic litter all other elements
          sb%pool%soluable_litter%nitrogen%ptr(ics:ice,:,iblk)   = sb%pool%soluable_litter%carbon%ptr(ics:ice,:,iblk) / & 
                                                                    (25._wp * molar_mass_N / molar_mass_C)
          sb%pool%soluable_litter%phosphorus%ptr(ics:ice,:,iblk) = sb%pool%soluable_litter%nitrogen%ptr(ics:ice,:,iblk) / & 
                                                                    (14._wp * molar_mass_P / molar_mass_N) 
          sb%pool%soluable_litter%carbon13%ptr(ics:ice,:,iblk)   = sb%pool%soluable_litter%carbon%ptr(ics:ice,:,iblk) / & 
                                                                    (1._wp + 1._wp / calc_mixing_ratio_C13C12(-26.0_wp))
          sb%pool%soluable_litter%carbon14%ptr(ics:ice,:,iblk)   = sb%pool%soluable_litter%carbon%ptr(ics:ice,:,iblk) * & 
                                                                    calc_mixing_ratio_C14C(-26.0_wp,2.0_wp) 
          sb%pool%soluable_litter%nitrogen15%ptr(ics:ice,:,iblk) = sb%pool%soluable_litter%nitrogen%ptr(ics:ice,:,iblk) / & 
                                                                    (1._wp + 1._wp / calc_mixing_ratio_N15N14(0.0_wp))
          ! non-woody litter all other elements
          sb%pool%polymeric_litter%nitrogen%ptr(ics:ice,:,iblk)   = sb%pool%polymeric_litter%carbon%ptr(ics:ice,:,iblk) / & 
                                                                    (100._wp * molar_mass_N / molar_mass_C)
          sb%pool%polymeric_litter%phosphorus%ptr(ics:ice,:,iblk) = sb%pool%polymeric_litter%nitrogen%ptr(ics:ice,:,iblk) / & 
                                                                    (14._wp * molar_mass_P / molar_mass_N) 
          sb%pool%polymeric_litter%carbon13%ptr(ics:ice,:,iblk)   = sb%pool%polymeric_litter%carbon%ptr(ics:ice,:,iblk) / & 
                                                                    (1._wp + 1._wp / calc_mixing_ratio_C13C12(-26.0_wp))
          sb%pool%polymeric_litter%carbon14%ptr(ics:ice,:,iblk)   = sb%pool%polymeric_litter%carbon%ptr(ics:ice,:,iblk) * & 
                                                                    calc_mixing_ratio_C14C(-26.0_wp,2.0_wp) 
          sb%pool%polymeric_litter%nitrogen15%ptr(ics:ice,:,iblk) = sb%pool%polymeric_litter%nitrogen%ptr(ics:ice,:,iblk) / & 
                                                                    (1._wp + 1._wp / calc_mixing_ratio_N15N14(0.0_wp))
          ! woody litter all other elements
          sb%pool%woody_litter%nitrogen%ptr(ics:ice,:,iblk)   = sb%pool%woody_litter%carbon%ptr(ics:ice,:,iblk) / & 
                                                                (260._wp * molar_mass_N / molar_mass_C)
          sb%pool%woody_litter%phosphorus%ptr(ics:ice,:,iblk) = sb%pool%woody_litter%nitrogen%ptr(ics:ice,:,iblk) / & 
                                                                (14._wp * molar_mass_P / molar_mass_N) 
          sb%pool%woody_litter%carbon13%ptr(ics:ice,:,iblk)   = sb%pool%woody_litter%carbon%ptr(ics:ice,:,iblk) / & 
                                                                (1._wp + 1._wp / calc_mixing_ratio_C13C12(-26.0_wp))
          sb%pool%woody_litter%carbon14%ptr(ics:ice,:,iblk)   = sb%pool%woody_litter%carbon%ptr(ics:ice,:,iblk) * & 
                                                                calc_mixing_ratio_C14C(-26.0_wp,2.0_wp) 
          sb%pool%woody_litter%nitrogen15%ptr(ics:ice,:,iblk) = sb%pool%woody_litter%nitrogen%ptr(ics:ice,:,iblk) / & 
                                                                (1._wp + 1._wp / calc_mixing_ratio_N15N14(0.0_wp))

          ! fast pool (=microbial) all other elements
          sb%pool%microbial%nitrogen%ptr(ics:ice,:,iblk)   = sb%pool%microbial%carbon%ptr(ics:ice,:,iblk) / microbial_cn
          sb%pool%microbial%phosphorus%ptr(ics:ice,:,iblk) = sb%pool%microbial%nitrogen%ptr(ics:ice,:,iblk) / microbial_np 
          sb%pool%microbial%carbon13%ptr(ics:ice,:,iblk)   = sb%pool%microbial%carbon%ptr(ics:ice,:,iblk) / & 
                                                             (1._wp + 1._wp / calc_mixing_ratio_C13C12(-26.0_wp))
          sb%pool%microbial%carbon14%ptr(ics:ice,:,iblk)   = sb%pool%microbial%carbon%ptr(ics:ice,:,iblk) * & 
                                                             calc_mixing_ratio_C14C(-26.0_wp,2.0_wp) 
          sb%pool%microbial%nitrogen15%ptr(ics:ice,:,iblk) = sb%pool%microbial%nitrogen%ptr(ics:ice,:,iblk) / & 
                                                             (1._wp + 1._wp / calc_mixing_ratio_N15N14(0.0_wp))

          sb%pool%residue%nitrogen%ptr(ics:ice,:,iblk)   = sb%pool%residue%carbon%ptr(ics:ice,:,iblk) / microbial_cn
          sb%pool%residue%phosphorus%ptr(ics:ice,:,iblk) = sb%pool%residue%nitrogen%ptr(ics:ice,:,iblk) / microbial_np 
          sb%pool%residue%carbon13%ptr(ics:ice,:,iblk)   = sb%pool%residue%carbon%ptr(ics:ice,:,iblk) / & 
                                                             (1._wp + 1._wp / calc_mixing_ratio_C13C12(-26.0_wp))
          sb%pool%residue%carbon14%ptr(ics:ice,:,iblk)   = sb%pool%residue%carbon%ptr(ics:ice,:,iblk) * & 
                                                             calc_mixing_ratio_C14C(-26.0_wp,2.0_wp) 
          sb%pool%residue%nitrogen15%ptr(ics:ice,:,iblk) = sb%pool%residue%nitrogen%ptr(ics:ice,:,iblk) / & 
                                                             (1._wp + 1._wp / calc_mixing_ratio_N15N14(0.0_wp))

          sb%nh4_n15_solute%ptr(ics:ice,:,iblk) = sb%nh4_solute%ptr(ics:ice,:,iblk) / ( 1._wp + 1._wp / calc_mixing_ratio_N15N14(0.0_wp)) 
          sb%no3_n15_solute%ptr(ics:ice,:,iblk) = sb%no3_solute%ptr(ics:ice,:,iblk) / ( 1._wp + 1._wp / calc_mixing_ratio_N15N14(0.0_wp)) 

       CASE('jsm')
 
         ! set microbial biomass to minimum possible density
         sb%pool%microbial%carbon%ptr(ics:ice,:,iblk)     = min_mic_biomass
         sb%pool%microbial%nitrogen%ptr(ics:ice,:,iblk)   = sb%pool%microbial%carbon%ptr(ics:ice,:,iblk) / microbial_cn 
         sb%pool%microbial%phosphorus%ptr(ics:ice,:,iblk) = sb%pool%microbial%nitrogen%ptr(ics:ice,:,iblk) / microbial_np
         sb%pool%microbial%carbon13%ptr(ics:ice,:,iblk)   = sb%pool%microbial%carbon%ptr(ics:ice,:,iblk) / ( 1._wp + 1._wp / calc_mixing_ratio_C13C12(-26.0_wp))
         sb%pool%microbial%carbon14%ptr(ics:ice,:,iblk)   = sb%pool%microbial%carbon%ptr(ics:ice,:,iblk) * calc_mixing_ratio_C14C(-26.0_wp,2.0_wp) 
         sb%pool%microbial%nitrogen15%ptr(ics:ice,:,iblk) = sb%pool%microbial%nitrogen%ptr(ics:ice,:,iblk)  / ( 1._wp + 1._wp / calc_mixing_ratio_N15N14(0.0_wp))
     
         ! some initial values for testing
         sb%nh4_solute%ptr(ics:ice,:,iblk) = 1.0_wp 
         sb%nh4_n15_solute%ptr(ics:ice,:,iblk) = sb%nh4_solute%ptr(ics:ice,:,iblk) / ( 1._wp + 1._wp / calc_mixing_ratio_N15N14(0.0_wp)) 
         sb%no3_solute%ptr(ics:ice,:,iblk) = 1.0_wp 
         sb%no3_n15_solute%ptr(ics:ice,:,iblk) = sb%no3_solute%ptr(ics:ice,:,iblk) / ( 1._wp + 1._wp / calc_mixing_ratio_N15N14(0.0_wp)) 
         sb%po4_solute%ptr(ics:ice,:,iblk) = 1.0_wp 

         ! inital Qmax based on current soil carbon content
         sb%bulk_soil_carbon_sl%ptr(ics:ice,:,iblk) = calc_bulk_soil_carbon(options, sb%pool)
         sb%bulk_dens_corr_sl%ptr(ics:ice,:,iblk) = calc_bulk_density_correction(sb%bulk_soil_carbon_sl%ptr(ics:ice,:,iblk), & 
                                                                         soil%bulk_dens_sl%ptr(ics:ice,:,iblk))
         sb%qmax_org%ptr(ics:ice,:,iblk) = calc_qmax_org(sb%bulk_dens_corr_sl%ptr(ics:ice,:,iblk))
 
         ! initialise diffusion by bioturbation
         sb%k_bioturb%ptr(ics:ice,:,iblk) = calc_bioturbation_rate(options,sb%bulk_dens_corr_sl%ptr(ics:ice,:,iblk)) 

       CASE DEFAULT
         CALL finish(routine,'Invalid setting of sb_model_scheme.')
 
    END SELECT

  END SUBROUTINE sb_init

END MODULE mo_sb_init
