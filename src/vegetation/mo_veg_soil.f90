!>
!! @brief diverse subroutines and functions to mimick soil for the stand-alone plant-model
!!
!! @author
!!  Sönke Zaehle, MPI-BGC Jena
!!
!! @par Revision History
!! First version by Sönke Zaehle (2016-03-31)
!!
MODULE mo_veg_soil

#include "dsl_jsbach.inc"

  USE mo_kind,               ONLY: wp, eps8
  USE mo_jsb_lctlib_class,   ONLY: t_lctlib_element
  USE mo_veg_config_class,   ONLY: veg_config
  USE mo_veg_constants


  IMPLICIT NONE
  PRIVATE
  PUBLIC :: update_vegsoil_bgc

  CHARACTER(len=*), PARAMETER :: modname = 'mo_veg_soil'

CONTAINS

  !-----------------------------------------------------------------------------------------------------
  ! Main Task
  ! 
  !-----------------------------------------------------------------------------------------------------
  !> routine to simulate N and P uptake from soil in the absence of a proper soil model
  !! 
  !! @par alternative use
  !! supply N and P to maintain a certain target C:N:P of the labile pool
  !! 
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE update_vegsoil_bgc(tile)

    USE mo_jsb_class,             ONLY: Get_model
    USE mo_base_memory,           ONLY: t_slm_tile, t_model
    USE mo_base_constants,        ONLY: dtime, one_day
    USE mo_jsb_process_class,     ONLY: ATM2LAND_, LAND2ATM_, SEB_, TURB_, SOIL_, HYDRO_, VEG_, HD_, RAD_, ASSIMI_, PHENO_, SB_

    ! Use of process memories
    dsl4jsb_Use_memory(VEG_)        !  USE mo_veg_memory_class,      ONLY: t_veg_memory
    dsl4jsb_Use_memory(SB_)         !  USE mo_sb_memory_class,       ONLY: t_sb_memory

    IMPLICIT NONE
    ! 0.1 InOut
    TYPE(t_slm_tile),INTENT(inout), TARGET   :: tile      !< one tile
    ! 0.2 Local
    TYPE(t_model),          POINTER        :: model
    TYPE(t_lctlib_element), POINTER        :: lctlib      !< land-cover-type library - parameter across pft's
    REAL(wp) :: dcost                                     ! cost_n_uptake_root - cost_n_fixation (mol C / mol N)
    REAL(wp) :: n_demand                                  ! N demand of plants (mu mol N /m2 / s) 
    REAL(wp) :: f_nacq_bnf,f_nacq_nup                     ! fraction of N uptake as BFN or Nuptake
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':update_vegsoil_bgc'
    ! 0.3 Memory allocation
    ! Declare process configuration and memory Pointers
    dsl4jsb_Def_memory(VEG_)                          ! TYPE(t_veg_memory),      POINTER :: veg__mem
    dsl4jsb_Def_memory(SB_)                           ! TYPE(t_sb_memory),      POINTER :: sb__mem
    TYPE(t_veg_memory), POINTER             :: veg        ! local pointer to veg__mem for SLM usability
    TYPE(t_sb_memory), POINTER              :: sb         ! local pointer to sb__mem for SLM usability
    dsl4jsb_Get_memory(VEG_)
    dsl4jsb_Get_memory(SB_)

    veg => veg__mem 
    sb  => sb__mem 

    model  => Get_model(tile%owner_model_id)
    lctlib => model%lctlib(tile%lcts%lib_id)

    !> 1.0 Calculate an instantaneous scalar of nutrient demand
    !!     taken as difference of the the weekly growth requirement of the nutrient to the stoichiometry of the labile pool
    !!
    !! this probably needs some thought - f_n_demand should not grow too large, but go to zero, if the plant has too much N/P 
    IF(veg%pool%labile%carbon%ptr > eps8) THEN
       veg%f_n_demand%ptr = lctlib%cn_leaf * & 
           MAX(veg%growth_req_n_tuptake_mavg%ptr - veg%pool%labile%nitrogen%ptr / veg%pool%labile%carbon%ptr, 0.0_wp) 
    ELSE
       veg%f_n_demand%ptr = 1.0_wp
    ENDIF
    IF(veg%pool%labile%nitrogen%ptr > eps8) THEN
       veg%f_p_demand%ptr = lctlib%np_leaf * & 
           MAX(veg%growth_req_p_tuptake_mavg%ptr - veg%pool%labile%phosphorus%ptr / veg%pool%labile%nitrogen%ptr, 0.0_wp)
    ELSE
       veg%f_p_demand%ptr = 1.0_wp
    ENDIF

    !> 2.0 inorganic N uptake (first guess)
    !!
    !!
    SELECT CASE (veg_config%n_uptake_scheme)
      !> 2.1 Case 1  maintain prescribed C:N:P irrespective of nutrient uptake capacity == no nutrient limitation case 
      !!
      CASE ("prescribed")
         ! uptake required to maintain a presribed C:N ratio of the labile pool in µmol N / m2 / s
         veg%uptake_n%ptr = MAX(veg%pool%labile%carbon%ptr / lctlib%cn_leaf - & 
                                           veg%pool%labile%nitrogen%ptr, 0.0_wp) * 1.e6_wp / dtime
         ! uptake required to maintain a presribed N:P ratio of the labile pool in µmol P / m2 / s
         veg%uptake_p%ptr = MAX(veg%pool%labile%nitrogen%ptr / lctlib%np_leaf + & 
                                              veg%uptake_n%ptr * dtime / 1.e6 - &
                                              veg%pool%labile%phosphorus%ptr, 0.0_wp) * 1.e6_wp / dtime

      !> 2.2 Case 2. Calculate inorganic N uptake according to Michaelis-Mnten Kinetics  
      !!
      CASE ("standard","optimal","resistance")
         ! OCN's MM Kinetics
         veg%uptake_n%ptr = lctlib%vmax_uptake_n * &
                 veg%demand_uptake_n_tuptake_mavg%ptr * & 
                 sb%nh4_solute%ptr(1,1,1) * ( lctlib%soil_n_aff1 + & 
                   1._wp / ( lctlib%soil_n_aff2 + sb%nh4_solute%ptr(1,1,1) )) * & 
                 veg%pool%fine_root%carbon%ptr
         veg%uptake_p%ptr = lctlib%vmax_uptake_p * & 
                 veg%demand_uptake_p_tuptake_mavg%ptr * & 
                 sb%po4_solute%ptr(1,1,1) * ( lctlib%soil_p_aff1 + & 
                   1._wp / ( lctlib%soil_p_aff2 + sb%po4_solute%ptr(1,1,1) )) * & 
                 veg%pool%fine_root%carbon%ptr

    ENDSELECT

    !> 3.0 N fixation
    !!
    !!
    SELECT CASE (veg_config%n_uptake_scheme)
      !> 3.1 Case 1: no N fixation
      !!
      CASE ("prescribed") 
        veg%n_fixation%ptr = 0.0_wp
      CASE("standard")
        IF(sb%nh4_solute%ptr(1,1,1) < soil_n_threshold_bnf) veg%n_fixation%ptr = lctlib%bnf_base 
      !> 3.2 Case 2: optimal N fixation
      !! BNF calculation with Michaelis Menten function (Rastetter, 2001) (MM)
      CASE ("optimal")
        IF (veg%cost_n_fixation_tuptake_mavg%ptr <= veg%cost_n_uptake_root%ptr) THEN 
           dcost  = veg%cost_n_uptake_root%ptr - veg%cost_n_fixation_tuptake_mavg%ptr
           veg%n_fixation%ptr = veg%pool%fine_root%carbon%ptr * & 
                              vmax_bnf * dcost / (k_bnf + dcost)
        ELSE
           veg%n_fixation%ptr = 0.0_wp
        ENDIF
      !> 3.3 Case 3: Ohm's law analogy 
      !! BNF calculation after FUN2.0 (Brzostek, 2014) (RA)
      CASE ("resistance")
         ! current N deficit of growth assuming default labile turnover rate (µmol N / m2 / s)
         n_demand = veg%demand_uptake_n_tuptake_mavg%ptr / lctlib%cn_leaf * veg%pool%labile%carbon%ptr / tau_labile / one_day * 1.e6
         f_nacq_bnf = (veg%cost_n_uptake_root%ptr - (transform_cost_nh4+transform_cost_no3)/2.0_wp) / &
                      (veg%cost_n_uptake_root%ptr+veg%cost_n_fixation_tuptake_mavg%ptr)
         f_nacq_nup = 1.0_wp - f_nacq_bnf
  
         veg%n_fixation%ptr = f_nacq_bnf * n_demand 
         veg%uptake_n%ptr   = MIN(f_nacq_nup * n_demand, veg%uptake_n%ptr)

    ENDSELECT 

    !> 3.0 Calculate potential and actual N uptake per unit fine root carbon
    !!
    veg%unit_uptake_n_pot%ptr = lctlib%vmax_uptake_n * &
                                sb%nh4_solute%ptr(1,1,1) * ( lctlib%soil_n_aff1 + &
                                1._wp / ( lctlib%soil_n_aff2 + sb%nh4_solute%ptr(1,1,1) ))
    veg%unit_uptake_p_pot%ptr = lctlib%vmax_uptake_p * &
                                sb%po4_solute%ptr(1,1,1) * ( lctlib%soil_p_aff1 + &
                                1._wp / ( lctlib%soil_p_aff2 + sb%po4_solute%ptr(1,1,1) ))
    IF(veg%pool%fine_root%carbon%ptr > eps8) THEN
      veg%unit_uptake_n_act%ptr = veg%uptake_n%ptr / veg%pool%fine_root%carbon%ptr
      veg%unit_uptake_p_act%ptr = veg%uptake_p%ptr / veg%pool%fine_root%carbon%ptr
    ELSE
      veg%unit_uptake_n_act%ptr = veg%unit_uptake_n_pot%ptr * veg%demand_uptake_n_tuptake_mavg%ptr
      veg%unit_uptake_p_act%ptr = veg%unit_uptake_p_pot%ptr * veg%demand_uptake_p_tuptake_mavg%ptr
    ENDIF

  END SUBROUTINE update_vegsoil_bgc

END MODULE mo_veg_soil
