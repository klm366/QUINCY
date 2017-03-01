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
  SUBROUTINE update_vegsoil_bgc(tile, options)

    USE mo_jsb_class,             ONLY: Get_model
    USE mo_base_memory,           ONLY: t_slm_tile, t_model, t_jsb_task_options
    USE mo_base_constants,        ONLY: dtime, one_day
    USE mo_jsb_process_class,     ONLY: A2L_, L2A_, SEB_, TURB_, SOIL_, HYDRO_, VEG_, HD_, RAD_, ASSIMI_, PHENO_, SB_

    ! Use of process memories
    dsl4jsb_Use_memory(VEG_)        !  USE mo_veg_memory_class,      ONLY: t_veg_memory
    dsl4jsb_Use_memory(SB_)         !  USE mo_sb_memory_class,       ONLY: t_sb_memory
    dsl4jsb_Use_memory(SOIL_)       !  USE mo_soil_memory_class,     ONLY: t_soil_memory

    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    TYPE(t_slm_tile),           INTENT(inout)     :: tile         !< one tile
    TYPE(t_jsb_task_options),   INTENT(in)        :: options      !< model options
    ! ---------------------------
    ! 0.2 Local
    TYPE(t_model),          POINTER       :: model
    TYPE(t_lctlib_element), POINTER       :: lctlib     !< land-cover-type library - parameter across pft's
    REAL(wp),DIMENSION(options%nc)        :: dcost        , &  ! cost_n_uptake_root - cost_n_fixation (mol C / mol N)
                                             n_demand     , &  ! N demand of plants (mu mol N /m2 / s) 
                                             f_nacq_bnf   , &  ! fraction of N uptake as BFN
                                             f_nacq_nup   , &  ! fraction of N uptake as Nuptake
                                             scal
    INTEGER                               :: iblk, ics, ice, nc
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':update_vegsoil_bgc'
    ! ---------------------------
    ! 0.3 Declare Memory
    ! Declare process configuration and memory Pointers
    dsl4jsb_Def_memory(VEG_)                          ! TYPE(t_veg_memory),      POINTER :: veg__mem
    dsl4jsb_Def_memory(SB_)                           ! TYPE(t_sb_memory),       POINTER :: sb__mem
    dsl4jsb_Def_memory(SOIL_)                         ! TYPE(t_soil_memory),     POINTER :: soil__mem
    ! Declare pointers to variables in memory
    dsl4jsb_Real3D_onChunk      :: soil_depth_sl                 ! SOIL_
    dsl4jsb_Real2D_onChunk      :: f_n_demand                   ! VEG_
    dsl4jsb_Real2D_onChunk      :: growth_req_n_tuptake_mavg
    dsl4jsb_Real2D_onChunk      :: f_p_demand
    dsl4jsb_Real2D_onChunk      :: growth_req_p_tuptake_mavg
    dsl4jsb_Real2D_onChunk      :: uptake_nh4
    dsl4jsb_Real2D_onChunk      :: demand_uptake_n_tuptake_mavg
    dsl4jsb_Real2D_onChunk      :: uptake_no3
    dsl4jsb_Real2D_onChunk      :: uptake_p
    dsl4jsb_Real2D_onChunk      :: demand_uptake_p_tuptake_mavg
    dsl4jsb_Real2D_onChunk      :: n_fixation
    dsl4jsb_Real2D_onChunk      :: cost_n_fixation_tuptake_mavg
    dsl4jsb_Real2D_onChunk      :: cost_n_uptake_root
    dsl4jsb_Real2D_onChunk      :: unit_uptake_n_pot
    dsl4jsb_Real2D_onChunk      :: unit_uptake_p_pot
    dsl4jsb_Real2D_onChunk      :: unit_uptake_n_act
    dsl4jsb_Real2D_onChunk      :: unit_uptake_p_act
! 2 be removed *slmtemp
    TYPE(t_veg_memory), POINTER             :: veg        ! local pointer to veg__mem for SLM usability
    TYPE(t_sb_memory), POINTER              :: sb         ! local pointer to sb__mem for SLM usability
    ! Get local variables from options argument
    iblk    = options%iblk
    ics     = options%ics
    ice     = options%ice
    nc      = options%nc
    ! ---------------------------
    ! 0.4 Process Activity, Debug Option
!    IF (.NOT. tile%Is_process_active(VEG_)) RETURN
!    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')
    ! ---------------------------
    ! 0.5 Get Memory
    ! Get process config
!    dsl4jsb_Get_config(xx_)
    ! Get process memories
    dsl4jsb_Get_memory(VEG_)
    dsl4jsb_Get_memory(SB_)
    dsl4jsb_Get_memory(SOIL_)
    ! Get process variables (Set pointers to variables in memory)
    model  => Get_model(tile%owner_model_id)
    lctlib => model%lctlib(tile%lcts%lib_id) 
    ! ---------------------------
    dsl4jsb_Get_var3D_onChunk(SOIL_, soil_depth_sl)
    ! ---------------------------
    dsl4jsb_Get_var2D_onChunk(VEG_, f_n_demand)
    dsl4jsb_Get_var2D_onChunk(VEG_, growth_req_n_tuptake_mavg)
    dsl4jsb_Get_var2D_onChunk(VEG_, f_p_demand)
    dsl4jsb_Get_var2D_onChunk(VEG_, growth_req_p_tuptake_mavg)
    dsl4jsb_Get_var2D_onChunk(VEG_, uptake_nh4)
    dsl4jsb_Get_var2D_onChunk(VEG_, demand_uptake_n_tuptake_mavg)
    dsl4jsb_Get_var2D_onChunk(VEG_, uptake_no3)
    dsl4jsb_Get_var2D_onChunk(VEG_, uptake_p)
    dsl4jsb_Get_var2D_onChunk(VEG_, demand_uptake_p_tuptake_mavg)
    dsl4jsb_Get_var2D_onChunk(VEG_, n_fixation)
    dsl4jsb_Get_var2D_onChunk(VEG_, cost_n_fixation_tuptake_mavg)
    dsl4jsb_Get_var2D_onChunk(VEG_, cost_n_uptake_root)
    dsl4jsb_Get_var2D_onChunk(VEG_, unit_uptake_n_pot)
    dsl4jsb_Get_var2D_onChunk(VEG_, unit_uptake_p_pot)
    dsl4jsb_Get_var2D_onChunk(VEG_, unit_uptake_n_act)
    dsl4jsb_Get_var2D_onChunk(VEG_, unit_uptake_p_act)


! 2 be removed  *slmtemp
    veg  => veg__mem 
    sb   => sb__mem 


    ! ------------------------------------------------------------------------------------------------------------
    ! Go science
    ! ------------------------------------------------------------------------------------------------------------



    !> 1.0 Calculate an instantaneous scalar of nutrient demand
    !!     taken as difference of the the weekly growth requirement of the nutrient to the stoichiometry of the labile pool
    !!
    !! this probably needs some thought - f_n_demand should not grow too large, but go to zero, if the plant has too much N/P 
    IF(veg%pool%labile%carbon%ptr > eps8) THEN
       f_n_demand(:) = lctlib%cn_leaf * & 
           MAX(growth_req_n_tuptake_mavg(:) - veg%pool%labile%nitrogen%ptr / veg%pool%labile%carbon%ptr, 0.0_wp) 
    ELSE
       f_n_demand(:) = 1.0_wp
    ENDIF

    IF(veg%pool%labile%nitrogen%ptr > eps8) THEN
       f_p_demand(:) = lctlib%np_leaf * & 
           MAX(growth_req_p_tuptake_mavg(:) - veg%pool%labile%phosphorus%ptr / veg%pool%labile%nitrogen%ptr, 0.0_wp)
    ELSE
       f_p_demand(:) = 1.0_wp
    ENDIF

    !> 2.0 inorganic uptake (first guess)
    !!
    !!

    !> 2.1 Nitrogen uptake
    !!
    SELECT CASE (veg_config%n_uptake_scheme)
      !> 2.1.1 Case 1  maintain prescribed C:N irrespective of nutrient uptake capacity == no nutrient limitation case 
      !!
      CASE ("prescribed")
         ! uptake required to maintain a presribed C:N ratio of the labile pool in µmol N / m2 / s
         uptake_nh4(:) = MAX(veg%pool%labile%carbon%ptr / lctlib%cn_leaf - veg%pool%labile%nitrogen%ptr, 0.0_wp) * 1.e6_wp / dtime

      !> 2.1.2 Case 2. Calculate inorganic N uptake according to Michaelis-Mnten Kinetics  
      !!
      CASE ("standard","optimal","resistance")
         ! OCN's MM Kinetics
         uptake_nh4(:) = lctlib%vmax_uptake_n * &
                 demand_uptake_n_tuptake_mavg(:) * & 
                 sb%nh4_solute%ptr(:,1,1) * soil_depth_sl(:,1) * ( lctlib%soil_n_aff1 + & 
                   1._wp / ( lctlib%soil_n_aff2 + sb%nh4_solute%ptr(:,1,1) * soil_depth_sl(:,1) )) * & 
                   veg%pool%fine_root%carbon%ptr
         
         uptake_no3(:) = lctlib%vmax_uptake_n * &
                 demand_uptake_n_tuptake_mavg(:) * & 
                 sb%no3_solute%ptr(:,1,1) * soil_depth_sl(:,1) * ( lctlib%soil_n_aff1 + & 
                   1._wp / ( lctlib%soil_n_aff2 + sb%no3_solute%ptr(:,1,1) * soil_depth_sl(:,1) )) * & 
                 veg%pool%fine_root%carbon%ptr

    ENDSELECT

    !> 2.2 inorganic P uptake 
    !!
    SELECT CASE (veg_config%p_uptake_scheme)
      !> 2.1 Case 1  maintain prescribed N:P irrespective of nutrient uptake capacity == no nutrient limitation case 
      !!
      CASE ("prescribed")
         ! uptake required to maintain a presribed N:P ratio of the labile pool in µmol P / m2 / s
         uptake_p(:) = MAX(veg%pool%labile%nitrogen%ptr / lctlib%np_leaf + uptake_nh4(:) * dtime / 1.e6 - &
                           veg%pool%labile%phosphorus%ptr, 0.0_wp) * &
                          1.e6_wp / dtime

      !> 2.2 Case 2. Calculate inorganic P uptake according to Michaelis-Mnten Kinetics  
      !!
      CASE ("standard")
         ! OCN's MM Kinetics
         uptake_p(:) = lctlib%vmax_uptake_p * & 
                 demand_uptake_p_tuptake_mavg(:) * & 
                 sb%po4_solute%ptr(:,1,1) * soil_depth_sl(:,1) * ( lctlib%soil_p_aff1 + & 
                   1._wp / ( lctlib%soil_p_aff2 + sb%po4_solute%ptr(:,1,1) * soil_depth_sl(:,1) )) * & 
                 veg%pool%fine_root%carbon%ptr

    ENDSELECT

    !> 3.0 N fixation
    !!
    !!
    SELECT CASE (veg_config%n_uptake_scheme)
      !> 3.1 Case 1: no N fixation
      !!
      CASE ("prescribed") 
        n_fixation(:) = 0.0_wp

      CASE("standard")
        WHERE(sb%nh4_solute%ptr(:,1,1) * soil_depth_sl(:,1) < soil_n_threshold_bnf) n_fixation(:) = lctlib%bnf_base 

      !> 3.2 Case 2: optimal N fixation
      !! BNF calculation with Michaelis Menten function (Rastetter, 2001) (MM)
      CASE ("optimal")
        WHERE(cost_n_fixation_tuptake_mavg(:) <= cost_n_uptake_root(:))
           dcost(:)         = cost_n_uptake_root(:) - cost_n_fixation_tuptake_mavg(:)
           n_fixation(:) = veg%pool%fine_root%carbon%ptr * vmax_bnf * dcost(:) / (k_bnf + dcost(:))
        ELSEWHERE
           n_fixation(:) = 0.0_wp
        END WHERE

      !> 3.3 Case 3: Ohm's law analogy 
      !! BNF calculation after FUN2.0 (Brzostek, 2014) (RA)
      CASE ("resistance")
         ! current N deficit of growth assuming default labile turnover rate (µmol N / m2 / s)
!SZ: why is n-demand defined on a static leaf C:N and not the actual growth requirement for N?
         n_demand(:) = demand_uptake_n_tuptake_mavg(:) / lctlib%cn_leaf * &
                    veg%pool%labile%carbon%ptr / tau_labile / one_day * 1.e6

         f_nacq_bnf(:) = MAX((cost_n_uptake_root(:) - offset_bnf) / &
                      (cost_n_uptake_root(:) + cost_n_fixation_tuptake_mavg(:)), 0.0_wp)
         f_nacq_nup(:) = 1.0_wp - f_nacq_bnf(:)
  
         n_fixation(:) = f_nacq_bnf(:) * n_demand(:) 
         WHERE(uptake_nh4(:) + uptake_no3(:) > eps8 .AND. f_nacq_nup(:) * n_demand(:) < uptake_nh4(:) + uptake_no3(:) )
               scal(:)       = (f_nacq_nup(:) * n_demand(:)) / (uptake_nh4(:) + uptake_no3(:))
! this needs some thoughts...
               uptake_nh4(:) = uptake_nh4(:) * scal(:)
               uptake_no3(:) = uptake_no3(:) * scal(:)
         END WHERE
  
    ENDSELECT 

    !> 4.0 Calculate potential and actual N uptake per unit fine root carbon
    !!
    unit_uptake_n_pot(:) = lctlib%vmax_uptake_n * &
                                (sb%nh4_solute%ptr(:,1,1) * soil_depth_sl(:,1) * ( lctlib%soil_n_aff1 + &
                                 1._wp / ( lctlib%soil_n_aff2 + sb%nh4_solute%ptr(:,1,1) * soil_depth_sl(:,1))) + &
                                 sb%no3_solute%ptr(:,1,1) * soil_depth_sl(:,1) * ( lctlib%soil_n_aff1 + &
                                 1._wp / ( lctlib%soil_n_aff2 + sb%no3_solute%ptr(:,1,1) * soil_depth_sl(:,1))))
    unit_uptake_p_pot(:) = lctlib%vmax_uptake_p * &
                                sb%po4_solute%ptr(:,1,1) * soil_depth_sl(:,1) * ( lctlib%soil_p_aff1 + &
                                1._wp / ( lctlib%soil_p_aff2 + sb%po4_solute%ptr(:,1,1) * soil_depth_sl(:,1) ))

    IF(veg%pool%fine_root%carbon%ptr > eps8) THEN
      unit_uptake_n_act(:) = (uptake_nh4(:) + uptake_no3(:)) / veg%pool%fine_root%carbon%ptr
      unit_uptake_p_act(:) = uptake_p(:) / veg%pool%fine_root%carbon%ptr
    ELSE
      unit_uptake_n_act(:) = unit_uptake_n_pot(:) * demand_uptake_n_tuptake_mavg(:)
      unit_uptake_p_act(:) = unit_uptake_p_pot(:) * demand_uptake_p_tuptake_mavg(:)
    ENDIF

  END SUBROUTINE update_vegsoil_bgc

END MODULE mo_veg_soil
