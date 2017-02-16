!>
!! @brief Contains the science routines for the JSM version in SLM 
!!
!! @author
!!  Sönke Zaehle 
!!
!! @par Revision History
!! First version by Sönke Zaehle (20170124): adopted from JSM_interface and JSM_core
!!
MODULE mo_sb_jsm_processes

#include "dsl_jsbach.inc"

  USE mo_kind,            ONLY: wp, eps4, eps8
  USE mo_exception,       ONLY: message, finish
  USE mo_base_memory,     ONLY: t_jsb_task_options
  USE mo_base_constants,  ONLY: dtime
  USE mo_sb_memory_class, ONLY: t_slm_struct_sb_pool,t_slm_struct_sb_element

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: calc_transfer_litterfall_pools,calc_litter_turnover, calc_depolymerisation, & 
            calc_microbial_turnover,calc_microbial_growth, &
            calc_bulk_soil_carbon, calc_bulk_density_correction, calc_qmax_org, &
            calc_om_sorption, & 
            calc_jsm_pool_change, & 
            calc_temperature_rate_modifier,calc_oxygen_rate_modifier, & 
            calc_convert_bulk_to_liquid,calc_substrate_rate_modifier

  CHARACTER(len=*), PARAMETER :: modname = 'mo_sb_jsm_processes'

CONTAINS

  !-----------------------------------------------------------------------------------------------------
  ! Sub Task called from update_sb_pools
  !
  !-----------------------------------------------------------------------------------------------------
  !> Calculates the transfer of vegetation litter fall to JSM litter pools
  !! still requires some thought as to the stoichiometry of the litter and the fraction of soluable litter 
  !! provides the litterfall interface between vegetation and soil_biogeochemistry
  !!
  !! Input: litterfall from vegetation, root distribution, soil layer depth
  !!
  !! Output: increase in soluable, polymeric, and woody litter
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE calc_transfer_litterfall_pools(options,nsoil,soil_depth_l,root_fraction_l,litterfall, &
                                  formation_soluable_litter,formation_polymeric_litter,formation_woody_litter)
                                 

    USE mo_sb_constants
    USE mo_veg_memory_class,      ONLY: t_slm_struct_veg_pool

    IMPLICIT NONE

    ! 0.1 InOut
    TYPE(t_jsb_task_options),              INTENT(in)    :: options           !< model options
    INTEGER,                               INTENT(in)    :: nsoil             !< number of soil layers 
    REAL(wp), DIMENSION(options%nc,nsoil), INTENT(in)    :: soil_depth_l      !< depth of soil layer [m]
    REAL(wp), DIMENSION(options%nc,nsoil), INTENT(in)    :: root_fraction_l   !< fraction of roots per layer 
    TYPE(t_slm_struct_veg_pool),           INTENT(in)    :: litterfall        !< litter input from vegetation [mol/m2/time step] 
    TYPE(t_slm_struct_sb_element),         INTENT(inout) :: formation_soluable_litter
    TYPE(t_slm_struct_sb_element),         INTENT(inout) :: formation_polymeric_litter
    TYPE(t_slm_struct_sb_element),         INTENT(inout) :: formation_woody_litter
    ! 0.2 Local 
    INTEGER                                                 :: iblk, ics, ice, nc
    INTEGER                                                 :: isoil
    REAL(wp), DIMENSION(options%nc)                         :: hlp1
    REAL(wp), DIMENSION(options%nc,nsoil) :: hlp2
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_transfer_litterfall_pools'

    iblk  = options%iblk
    ics   = options%ics
    ice   = options%ice
    nc    = options%nc

    !> 1.0 pure soluable litter -> labile + reserve?
    !!
    formation_soluable_litter%carbon%ptr(ics:ice,1,iblk) = formation_soluable_litter%carbon%ptr(ics:ice,1,iblk) + &
        (litterfall%labile%carbon%ptr + litterfall%reserve%carbon%ptr) / soil_depth_l(:,1)
    formation_soluable_litter%nitrogen%ptr(ics:ice,1,iblk) = formation_soluable_litter%nitrogen%ptr(ics:ice,1,iblk) + &
        (litterfall%labile%nitrogen%ptr + litterfall%reserve%nitrogen%ptr) / soil_depth_l(:,1)
    formation_soluable_litter%phosphorus%ptr(ics:ice,1,iblk) = formation_soluable_litter%phosphorus%ptr(ics:ice,1,iblk) + &
        (litterfall%labile%phosphorus%ptr + litterfall%reserve%phosphorus%ptr) / soil_depth_l(:,1)
    formation_soluable_litter%carbon13%ptr(ics:ice,1,iblk) = formation_soluable_litter%carbon13%ptr(ics:ice,1,iblk) + &
        (litterfall%labile%carbon13%ptr + litterfall%reserve%carbon13%ptr) / soil_depth_l(:,1)
    formation_soluable_litter%carbon14%ptr(ics:ice,1,iblk) = formation_soluable_litter%carbon14%ptr(ics:ice,1,iblk) + &
        (litterfall%labile%carbon14%ptr + litterfall%reserve%carbon14%ptr) / soil_depth_l(:,1)
    formation_soluable_litter%nitrogen15%ptr(ics:ice,1,iblk) = formation_soluable_litter%nitrogen15%ptr(ics:ice,1,iblk) + &
        (litterfall%labile%nitrogen15%ptr + litterfall%reserve%nitrogen15%ptr) / soil_depth_l(:,1)

    !> 2.0 mixed soluable and polymeric stuff
    !! 
    !> 2.1 leaves, still needs a way to split nutrients between soluable and polymeric..., and calculate the fraction going DOM
    hlp1 = f_soluable_leaf 
    formation_soluable_litter%carbon%ptr(ics:ice,1,iblk) = formation_soluable_litter%carbon%ptr(ics:ice,1,iblk) + &
        hlp1 * litterfall%leaf%carbon%ptr / soil_depth_l(:,1)
    formation_polymeric_litter%carbon%ptr(ics:ice,1,iblk) = formation_polymeric_litter%carbon%ptr(ics:ice,1,iblk) + &
        (1._wp - hlp1) * litterfall%leaf%carbon%ptr / soil_depth_l(:,1)
    formation_soluable_litter%nitrogen%ptr(ics:ice,1,iblk) = formation_soluable_litter%nitrogen%ptr(ics:ice,1,iblk) + &
        hlp1 * litterfall%leaf%nitrogen%ptr / soil_depth_l(:,1)
    formation_polymeric_litter%nitrogen%ptr(ics:ice,1,iblk) = formation_polymeric_litter%nitrogen%ptr(ics:ice,1,iblk) + &
        (1._wp - hlp1) * litterfall%leaf%nitrogen%ptr / soil_depth_l(:,1)
    formation_soluable_litter%phosphorus%ptr(ics:ice,1,iblk) = formation_soluable_litter%phosphorus%ptr(ics:ice,1,iblk) + &
        hlp1 * litterfall%leaf%phosphorus%ptr / soil_depth_l(:,1)
    formation_polymeric_litter%phosphorus%ptr(ics:ice,1,iblk) = formation_polymeric_litter%phosphorus%ptr(ics:ice,1,iblk) + &
        (1._wp - hlp1) * litterfall%leaf%phosphorus%ptr / soil_depth_l(:,1)
    formation_soluable_litter%carbon13%ptr(ics:ice,1,iblk) = formation_soluable_litter%carbon13%ptr(ics:ice,1,iblk) + &
        hlp1 * litterfall%leaf%carbon13%ptr / soil_depth_l(:,1)
    formation_polymeric_litter%carbon13%ptr(ics:ice,1,iblk) = formation_polymeric_litter%carbon13%ptr(ics:ice,1,iblk) + &
        (1._wp - hlp1) * litterfall%leaf%carbon13%ptr / soil_depth_l(:,1)
    formation_soluable_litter%carbon14%ptr(ics:ice,1,iblk) = formation_soluable_litter%carbon14%ptr(ics:ice,1,iblk) + &
        hlp1 * litterfall%leaf%carbon14%ptr / soil_depth_l(:,1)
    formation_polymeric_litter%carbon14%ptr(ics:ice,1,iblk) = formation_polymeric_litter%carbon14%ptr(ics:ice,1,iblk) + &
        (1._wp - hlp1) * litterfall%leaf%carbon14%ptr / soil_depth_l(:,1)
    formation_soluable_litter%nitrogen15%ptr(ics:ice,1,iblk) = formation_soluable_litter%nitrogen15%ptr(ics:ice,1,iblk) + &
        hlp1 * litterfall%leaf%nitrogen15%ptr / soil_depth_l(:,1)
    formation_polymeric_litter%nitrogen15%ptr(ics:ice,1,iblk) = formation_polymeric_litter%nitrogen15%ptr(ics:ice,1,iblk) + &
        (1._wp - hlp1) * litterfall%leaf%nitrogen15%ptr / soil_depth_l(:,1)

    !> 2.2 fruit, still needs a way to split nutrients between soluable and polymeric..., and calculate the fraction going DOM
    !!
    hlp1 = f_soluable_fruit 
    formation_soluable_litter%carbon%ptr(ics:ice,1,iblk) = formation_soluable_litter%carbon%ptr(ics:ice,1,iblk) + &
        hlp1 * litterfall%fruit%carbon%ptr / soil_depth_l(:,1)
    formation_polymeric_litter%carbon%ptr(ics:ice,1,iblk) = formation_polymeric_litter%carbon%ptr(ics:ice,1,iblk) + &
        (1._wp - hlp1) * litterfall%fruit%carbon%ptr / soil_depth_l(:,1)
    formation_soluable_litter%nitrogen%ptr(ics:ice,1,iblk) = formation_soluable_litter%nitrogen%ptr(ics:ice,1,iblk) + &
        hlp1 * litterfall%fruit%nitrogen%ptr / soil_depth_l(:,1)
    formation_polymeric_litter%nitrogen%ptr(ics:ice,1,iblk) = formation_polymeric_litter%nitrogen%ptr(ics:ice,1,iblk) + &
        (1._wp - hlp1) * litterfall%fruit%nitrogen%ptr / soil_depth_l(:,1)
    formation_soluable_litter%phosphorus%ptr(ics:ice,1,iblk) = formation_soluable_litter%phosphorus%ptr(ics:ice,1,iblk) + &
        hlp1 * litterfall%fruit%phosphorus%ptr / soil_depth_l(:,1)
    formation_polymeric_litter%phosphorus%ptr(ics:ice,1,iblk) = formation_polymeric_litter%phosphorus%ptr(ics:ice,1,iblk) + &
        (1._wp - hlp1) * litterfall%fruit%phosphorus%ptr / soil_depth_l(:,1)
    formation_soluable_litter%carbon13%ptr(ics:ice,1,iblk) = formation_soluable_litter%carbon13%ptr(ics:ice,1,iblk) + &
        hlp1 * litterfall%fruit%carbon13%ptr / soil_depth_l(:,1)
    formation_polymeric_litter%carbon13%ptr(ics:ice,1,iblk) = formation_polymeric_litter%carbon13%ptr(ics:ice,1,iblk) + &
        (1._wp - hlp1) * litterfall%fruit%carbon13%ptr / soil_depth_l(:,1)
    formation_soluable_litter%carbon14%ptr(ics:ice,1,iblk) = formation_soluable_litter%carbon14%ptr(ics:ice,1,iblk) + &
        hlp1 * litterfall%fruit%carbon14%ptr / soil_depth_l(:,1)
    formation_polymeric_litter%carbon14%ptr(ics:ice,1,iblk) = formation_polymeric_litter%carbon14%ptr(ics:ice,1,iblk) + &
        (1._wp - hlp1) * litterfall%fruit%carbon14%ptr / soil_depth_l(:,1)
    formation_soluable_litter%nitrogen15%ptr(ics:ice,1,iblk) = formation_soluable_litter%nitrogen15%ptr(ics:ice,1,iblk) + &
        hlp1 * litterfall%fruit%nitrogen15%ptr / soil_depth_l(:,1)
    formation_polymeric_litter%nitrogen15%ptr(ics:ice,1,iblk) = formation_polymeric_litter%nitrogen15%ptr(ics:ice,1,iblk) + &
        (1._wp - hlp1) * litterfall%fruit%nitrogen15%ptr / soil_depth_l(:,1)

    !> 2.3 seed_bed, still needs a way to split nutrients between soluable and polymeric..., and calculate the fraction going DOM
    !!
    hlp1 = f_soluable_seed_bed 
    formation_soluable_litter%carbon%ptr(ics:ice,1,iblk) = formation_soluable_litter%carbon%ptr(ics:ice,1,iblk) + &
        hlp1 * litterfall%seed_bed%carbon%ptr / soil_depth_l(:,1)
    formation_polymeric_litter%carbon%ptr(ics:ice,1,iblk) = formation_polymeric_litter%carbon%ptr(ics:ice,1,iblk) + &
        (1._wp - hlp1) * litterfall%seed_bed%carbon%ptr / soil_depth_l(:,1)
    formation_soluable_litter%nitrogen%ptr(ics:ice,1,iblk) = formation_soluable_litter%nitrogen%ptr(ics:ice,1,iblk) + &
        hlp1 * litterfall%seed_bed%nitrogen%ptr / soil_depth_l(:,1)
    formation_polymeric_litter%nitrogen%ptr(ics:ice,1,iblk) = formation_polymeric_litter%nitrogen%ptr(ics:ice,1,iblk) + &
        (1._wp - hlp1) * litterfall%seed_bed%nitrogen%ptr / soil_depth_l(:,1)
    formation_soluable_litter%phosphorus%ptr(ics:ice,1,iblk) = formation_soluable_litter%phosphorus%ptr(ics:ice,1,iblk) + &
        hlp1 * litterfall%seed_bed%phosphorus%ptr / soil_depth_l(:,1)
    formation_polymeric_litter%phosphorus%ptr(ics:ice,1,iblk) = formation_polymeric_litter%phosphorus%ptr(ics:ice,1,iblk) + &
        (1._wp - hlp1) * litterfall%seed_bed%phosphorus%ptr / soil_depth_l(:,1)
    formation_soluable_litter%carbon13%ptr(ics:ice,1,iblk) = formation_soluable_litter%carbon13%ptr(ics:ice,1,iblk) + &
        hlp1 * litterfall%seed_bed%carbon13%ptr / soil_depth_l(:,1)
    formation_polymeric_litter%carbon13%ptr(ics:ice,1,iblk) = formation_polymeric_litter%carbon13%ptr(ics:ice,1,iblk) + &
        (1._wp - hlp1) * litterfall%seed_bed%carbon13%ptr / soil_depth_l(:,1)
    formation_soluable_litter%carbon14%ptr(ics:ice,1,iblk) = formation_soluable_litter%carbon14%ptr(ics:ice,1,iblk) + &
        hlp1 * litterfall%seed_bed%carbon14%ptr / soil_depth_l(:,1)
    formation_polymeric_litter%carbon14%ptr(ics:ice,1,iblk) = formation_polymeric_litter%carbon14%ptr(ics:ice,1,iblk) + &
        (1._wp - hlp1) * litterfall%seed_bed%carbon14%ptr / soil_depth_l(:,1)
    formation_soluable_litter%nitrogen15%ptr(ics:ice,1,iblk) = formation_soluable_litter%nitrogen15%ptr(ics:ice,1,iblk) + &
        hlp1 * litterfall%seed_bed%nitrogen15%ptr / soil_depth_l(:,1)
    formation_polymeric_litter%nitrogen15%ptr(ics:ice,1,iblk) = formation_polymeric_litter%nitrogen15%ptr(ics:ice,1,iblk) + &
        (1._wp - hlp1) * litterfall%seed_bed%nitrogen15%ptr / soil_depth_l(:,1)

    !> 3.0 roots
    !!
    DO isoil = 1,nsoil
       !> 3.1 fine roots
       !! still needs a way to split nutrients between soluable and polymeric..., and calculate the fraction going DOM
       hlp1 = f_soluable_fine_root 
       WHERE(soil_depth_l(:,isoil) > eps8)
          formation_soluable_litter%carbon%ptr(ics:ice,isoil,iblk) = formation_soluable_litter%carbon%ptr(ics:ice,isoil,iblk) + &
              hlp1 * root_fraction_l(:,isoil) * litterfall%fine_root%carbon%ptr / soil_depth_l(:,isoil)
          formation_polymeric_litter%carbon%ptr(ics:ice,isoil,iblk) = formation_polymeric_litter%carbon%ptr(ics:ice,isoil,iblk) + &
              (1._wp - hlp1) * root_fraction_l(:,isoil) * litterfall%fine_root%carbon%ptr / soil_depth_l(:,isoil)
          formation_soluable_litter%nitrogen%ptr(ics:ice,isoil,iblk) = formation_soluable_litter%nitrogen%ptr(ics:ice,isoil,iblk) + &
              hlp1 * root_fraction_l(:,isoil) * litterfall%fine_root%nitrogen%ptr / soil_depth_l(:,isoil)
          formation_polymeric_litter%nitrogen%ptr(ics:ice,isoil,iblk) = formation_polymeric_litter%nitrogen%ptr(ics:ice,isoil,iblk) + &
              (1._wp - hlp1) * root_fraction_l(:,isoil) * litterfall%fine_root%nitrogen%ptr / soil_depth_l(:,isoil)
          formation_soluable_litter%phosphorus%ptr(ics:ice,isoil,iblk) = formation_soluable_litter%phosphorus%ptr(ics:ice,isoil,iblk) + &
              hlp1 * root_fraction_l(:,isoil) * litterfall%fine_root%phosphorus%ptr / soil_depth_l(:,isoil)
          formation_polymeric_litter%phosphorus%ptr(ics:ice,isoil,iblk) = formation_polymeric_litter%phosphorus%ptr(ics:ice,isoil,iblk) + &
              (1._wp - hlp1) * root_fraction_l(:,isoil) * litterfall%fine_root%phosphorus%ptr / soil_depth_l(:,isoil)
          formation_soluable_litter%carbon13%ptr(ics:ice,isoil,iblk) = formation_soluable_litter%carbon13%ptr(ics:ice,isoil,iblk) + &
              hlp1 * root_fraction_l(:,isoil) * litterfall%fine_root%carbon13%ptr / soil_depth_l(:,isoil)
          formation_polymeric_litter%carbon13%ptr(ics:ice,isoil,iblk) = formation_polymeric_litter%carbon13%ptr(ics:ice,isoil,iblk) + &
              (1._wp - hlp1) * root_fraction_l(:,isoil) * litterfall%fine_root%carbon13%ptr / soil_depth_l(:,isoil)
          formation_soluable_litter%carbon14%ptr(ics:ice,isoil,iblk) = formation_soluable_litter%carbon14%ptr(ics:ice,isoil,iblk) + &
              hlp1 * root_fraction_l(:,isoil) * litterfall%fine_root%carbon14%ptr / soil_depth_l(:,isoil)
          formation_polymeric_litter%carbon14%ptr(ics:ice,isoil,iblk) = formation_polymeric_litter%carbon14%ptr(ics:ice,isoil,iblk) + &
              (1._wp - hlp1) * root_fraction_l(:,isoil) * litterfall%fine_root%carbon14%ptr / soil_depth_l(:,isoil)
          formation_soluable_litter%nitrogen15%ptr(ics:ice,isoil,iblk) = formation_soluable_litter%nitrogen15%ptr(ics:ice,isoil,iblk) + &
              hlp1 * root_fraction_l(:,isoil) * litterfall%fine_root%nitrogen15%ptr / soil_depth_l(:,isoil)
          formation_polymeric_litter%nitrogen15%ptr(ics:ice,isoil,iblk) = formation_polymeric_litter%nitrogen15%ptr(ics:ice,isoil,iblk) + &
              (1._wp - hlp1) * root_fraction_l(:,isoil) * litterfall%fine_root%nitrogen15%ptr / soil_depth_l(:,isoil)

          !> 3.2 woody litter below ground 
          !!
          formation_woody_litter%carbon%ptr(ics:ice,isoil,iblk) = formation_woody_litter%carbon%ptr(ics:ice,isoil,iblk) + & 
             root_fraction_l(:,isoil) * litterfall%coarse_root%carbon%ptr / soil_depth_l(:,isoil)
          formation_woody_litter%nitrogen%ptr(ics:ice,isoil,iblk) = formation_woody_litter%nitrogen%ptr(ics:ice,isoil,iblk) + & 
             root_fraction_l(:,isoil) * litterfall%coarse_root%nitrogen%ptr / soil_depth_l(:,isoil)
          formation_woody_litter%phosphorus%ptr(ics:ice,isoil,iblk) = formation_woody_litter%phosphorus%ptr(ics:ice,isoil,iblk) + & 
             root_fraction_l(:,isoil) * litterfall%coarse_root%phosphorus%ptr / soil_depth_l(:,isoil)
          formation_woody_litter%carbon13%ptr(ics:ice,isoil,iblk) = formation_woody_litter%carbon13%ptr(ics:ice,isoil,iblk) + & 
             root_fraction_l(:,isoil) * litterfall%coarse_root%carbon13%ptr / soil_depth_l(:,isoil)
          formation_woody_litter%carbon14%ptr(ics:ice,isoil,iblk) = formation_woody_litter%carbon14%ptr(ics:ice,isoil,iblk) + & 
             root_fraction_l(:,isoil) * litterfall%coarse_root%carbon14%ptr / soil_depth_l(:,isoil)
          formation_woody_litter%nitrogen15%ptr(ics:ice,isoil,iblk) = formation_woody_litter%nitrogen15%ptr(ics:ice,isoil,iblk) + & 
             root_fraction_l(:,isoil) * litterfall%coarse_root%nitrogen15%ptr / soil_depth_l(:,isoil)
       ENDWHERE
    ENDDO

    !> 4.0 woody litter above ground 
    !!
    formation_woody_litter%carbon%ptr(ics:ice,1,iblk) = formation_woody_litter%carbon%ptr(ics:ice,1,iblk) + & 
        (litterfall%sap_wood%carbon%ptr + litterfall%heart_wood%carbon%ptr) / soil_depth_l(:,1)
    formation_woody_litter%nitrogen%ptr(ics:ice,1,iblk) = formation_woody_litter%nitrogen%ptr(ics:ice,1,iblk) + & 
        (litterfall%sap_wood%nitrogen%ptr + litterfall%heart_wood%nitrogen%ptr) / soil_depth_l(:,1)
    formation_woody_litter%phosphorus%ptr(ics:ice,1,iblk) = formation_woody_litter%phosphorus%ptr(ics:ice,1,iblk) + & 
        (litterfall%sap_wood%phosphorus%ptr + litterfall%heart_wood%phosphorus%ptr) / soil_depth_l(:,1)
    formation_woody_litter%carbon13%ptr(ics:ice,1,iblk) = formation_woody_litter%carbon13%ptr(ics:ice,1,iblk) + & 
        (litterfall%sap_wood%carbon13%ptr + litterfall%heart_wood%carbon13%ptr) / soil_depth_l(:,1)
    formation_woody_litter%carbon14%ptr(ics:ice,1,iblk) = formation_woody_litter%carbon14%ptr(ics:ice,1,iblk) + & 
        (litterfall%sap_wood%carbon14%ptr + litterfall%heart_wood%carbon14%ptr) / soil_depth_l(:,1)
    formation_woody_litter%nitrogen15%ptr(ics:ice,1,iblk) = formation_woody_litter%nitrogen15%ptr(ics:ice,1,iblk) + & 
        (litterfall%sap_wood%nitrogen15%ptr + litterfall%heart_wood%nitrogen15%ptr) / soil_depth_l(:,1)

  END SUBROUTINE calc_transfer_litterfall_pools 

  !-----------------------------------------------------------------------------------------------------
  ! Sub Task called from update_sb_pools
  ! 
  !
  !-----------------------------------------------------------------------------------------------------
  !> calculate turnover of litter pools according to first order decay 
  !! provides a link between the soluable and woody litter pools introduced in SLM, and the
  !! DOM and polymeric pools of JSM
  !!
  !! Input: temperature and moisture dependence of litter decomposition, litter pool size, and its
  !!        turnover time at standard conditions 
  !! Output: transfer of litter to its sink pool (soluable -> DOM; woody -> polymeric litter) 
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE calc_litter_turnover(options,nsoil,tau_litter_pool,rtm_depolymerisation,rmm_depolymerisation, & 
                                  pool_litter,loss_litter,formation_sink)

    IMPLICIT NONE

    ! 0.1 InOut
    TYPE(t_jsb_task_options),              INTENT(in)    :: options               !< model options
    INTEGER,                               INTENT(in)    :: nsoil                 !< number of soil layers 
    REAL(wp),                              INTENT(in)    :: tau_litter_pool       !< turnover time of the litter pool [1/s]
    REAL(wp), DIMENSION(options%nc,nsoil), INTENT(in)    :: rtm_depolymerisation  !< temperature rate modifier [-]
    REAL(wp), DIMENSION(options%nc,nsoil), INTENT(in)    :: rmm_depolymerisation  !< moisture rate modifier [-]
    TYPE(t_slm_struct_sb_element),         INTENT(in)    :: pool_litter           !< litter pool to decompose [mol/m3]
    TYPE(t_slm_struct_sb_element),         INTENT(inout) :: loss_litter           !< reduction in litter pool [mol/m3/timestep]
    TYPE(t_slm_struct_sb_element),         INTENT(inout) :: formation_sink        !< formation of SOM pool [mol/m3/timestep]
    ! 0.2 Local 
    INTEGER                                                 :: iblk, ics, ice, nc
    REAL(wp), DIMENSION(options%nc,nsoil) :: hlp1,decay_rate
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_litter_turnover'

    iblk  = options%iblk
    ics   = options%ics
    ice   = options%ice
    nc    = options%nc

    !> 1.0 calulate decay per timestep
    !!
    decay_rate = 1._wp/tau_litter_pool * dtime

    !> 2.0 calculate source and sink terms
    !! assuming to change in composition or mass of the material to decompose
    !! 
    hlp1 = decay_rate * pool_litter%carbon%ptr(ics:ice,:,iblk)
    loss_litter%carbon%ptr(ics:ice,:,iblk)    = loss_litter%carbon%ptr(ics:ice,:,iblk) + hlp1
    formation_sink%carbon%ptr(ics:ice,:,iblk) = formation_sink%carbon%ptr(ics:ice,:,iblk) + hlp1 

    hlp1 = decay_rate * pool_litter%nitrogen%ptr(ics:ice,:,iblk)
    loss_litter%nitrogen%ptr(ics:ice,:,iblk)    = loss_litter%nitrogen%ptr(ics:ice,:,iblk) + hlp1
    formation_sink%nitrogen%ptr(ics:ice,:,iblk) = formation_sink%nitrogen%ptr(ics:ice,:,iblk) + hlp1 

    hlp1 = decay_rate * pool_litter%phosphorus%ptr(ics:ice,:,iblk)
    loss_litter%phosphorus%ptr(ics:ice,:,iblk)    = loss_litter%phosphorus%ptr(ics:ice,:,iblk) + hlp1
    formation_sink%phosphorus%ptr(ics:ice,:,iblk) = formation_sink%phosphorus%ptr(ics:ice,:,iblk) + hlp1 

    hlp1 = decay_rate * pool_litter%carbon13%ptr(ics:ice,:,iblk)
    loss_litter%carbon13%ptr(ics:ice,:,iblk)    = loss_litter%carbon13%ptr(ics:ice,:,iblk) + hlp1
    formation_sink%carbon13%ptr(ics:ice,:,iblk) = formation_sink%carbon13%ptr(ics:ice,:,iblk) + hlp1 

    hlp1 = decay_rate * pool_litter%carbon14%ptr(ics:ice,:,iblk)
    loss_litter%carbon14%ptr(ics:ice,:,iblk)    = loss_litter%carbon14%ptr(ics:ice,:,iblk) + hlp1
    formation_sink%carbon14%ptr(ics:ice,:,iblk) = formation_sink%carbon14%ptr(ics:ice,:,iblk) + hlp1 

    hlp1 = decay_rate * pool_litter%nitrogen15%ptr(ics:ice,:,iblk)
    loss_litter%nitrogen15%ptr(ics:ice,:,iblk)    = loss_litter%nitrogen15%ptr(ics:ice,:,iblk) + hlp1
    formation_sink%nitrogen15%ptr(ics:ice,:,iblk) = formation_sink%nitrogen15%ptr(ics:ice,:,iblk) + hlp1 

  END SUBROUTINE calc_litter_turnover

  !-----------------------------------------------------------------------------------------------------
  ! Sub Task called from update_sb_pools
  !
  !-----------------------------------------------------------------------------------------------------
  !> calculate turnover of the microbial pool 
  !! corresponds to the term pi * Cb in Ahrens et al. 2015, eq. 7 
  !! 
  !! Input: microbial biomass 
  !!
  !! Output: transfer of microbial biomass to its sink pools (soluable -> DOM; polymeric -> residue) 
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE calc_microbial_turnover(options,nsoil,pool_microbial,loss_microbial,formation_dom,formation_residue)

    USE mo_sb_constants,          ONLY: k_mort_microbial, min_mic_biomass, f_soluable_residue

    IMPLICIT NONE

    ! 0.1 InOut
    TYPE(t_jsb_task_options),       INTENT(in)    :: options            !< model options
    INTEGER,                        INTENT(in)    :: nsoil              !< number of soil layers 
    TYPE(t_slm_struct_sb_element),  INTENT(in)    :: pool_microbial     !< microbial biomass [mol/m3]
    TYPE(t_slm_struct_sb_element),  INTENT(inout) :: loss_microbial     !< microbial death [mol/m3/timestep]
    TYPE(t_slm_struct_sb_element),  INTENT(inout) :: formation_dom      !< DOM formation from microbial death [mol/m3/timestep]
    TYPE(t_slm_struct_sb_element),  INTENT(inout) :: formation_residue  !< residue formation from microbial death [mol/m3/timestep]
    ! 0.2 Local 
    INTEGER                               :: iblk, ics, ice, nc
    REAL(wp), DIMENSION(options%nc,nsoil) :: hlp1, mortality_rate
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_microbial_turnover'

    iblk  = options%iblk
    ics   = options%ics
    ice   = options%ice
    nc    = options%nc

    !> 1.0 calculate microbial mortality rate, limited to minimum microbial biomass density (-> dormant microbes)
    !!
    mortality_rate = k_mort_microbial * dtime
    hlp1 = mortality_rate * pool_microbial%carbon%ptr(ics:ice,:,iblk) 
    WHERE(pool_microbial%carbon%ptr(ics:ice,:,iblk) - hlp1 < min_mic_biomass)
      mortality_rate = 1._wp - (min_mic_biomass / pool_microbial%carbon%ptr(ics:ice,:,iblk))
    ENDWHERE

    !> 2.0 calculate matter fluxes associated with microbial death
    !! a fraction turns to DOM, the remainder becomes residue
    !! the assumption is made that the composition of this material is not changed
    !!@todo: define C:N:P stoichiometries of DOM and residues from microbial death 
    !!
    hlp1 = mortality_rate * pool_microbial%carbon%ptr(ics:ice,:,iblk)
    loss_microbial%carbon%ptr(ics:ice,:,iblk)    = loss_microbial%carbon%ptr(ics:ice,:,iblk) + hlp1
    formation_dom%carbon%ptr(ics:ice,:,iblk)     = formation_dom%carbon%ptr(ics:ice,:,iblk) + & 
                                                   f_soluable_residue * hlp1 
    formation_residue%carbon%ptr(ics:ice,:,iblk) = formation_dom%carbon%ptr(ics:ice,:,iblk) + & 
                                                   (1._wp - f_soluable_residue) * hlp1 

    hlp1 = mortality_rate * pool_microbial%nitrogen%ptr(ics:ice,:,iblk)
    loss_microbial%nitrogen%ptr(ics:ice,:,iblk)    = loss_microbial%nitrogen%ptr(ics:ice,:,iblk) + hlp1
    formation_dom%nitrogen%ptr(ics:ice,:,iblk)     = formation_dom%nitrogen%ptr(ics:ice,:,iblk) + & 
                                                     f_soluable_residue * hlp1 
    formation_residue%nitrogen%ptr(ics:ice,:,iblk) = formation_dom%nitrogen%ptr(ics:ice,:,iblk) + & 
                                                     (1._wp - f_soluable_residue) * hlp1 

    hlp1 = mortality_rate * pool_microbial%phosphorus%ptr(ics:ice,:,iblk)
    loss_microbial%phosphorus%ptr(ics:ice,:,iblk)    = loss_microbial%phosphorus%ptr(ics:ice,:,iblk) + hlp1
    formation_dom%phosphorus%ptr(ics:ice,:,iblk)     = formation_dom%phosphorus%ptr(ics:ice,:,iblk) + & 
                                                       f_soluable_residue * hlp1 
    formation_residue%phosphorus%ptr(ics:ice,:,iblk) = formation_dom%phosphorus%ptr(ics:ice,:,iblk) + & 
                                                       (1._wp - f_soluable_residue) * hlp1 

    hlp1 = mortality_rate * pool_microbial%carbon13%ptr(ics:ice,:,iblk)
    loss_microbial%carbon13%ptr(ics:ice,:,iblk)    = loss_microbial%carbon13%ptr(ics:ice,:,iblk) + hlp1
    formation_dom%carbon13%ptr(ics:ice,:,iblk)     = formation_dom%carbon13%ptr(ics:ice,:,iblk) + & 
                                                     f_soluable_residue * hlp1 
    formation_residue%carbon13%ptr(ics:ice,:,iblk) = formation_dom%carbon13%ptr(ics:ice,:,iblk) + & 
                                                     (1._wp - f_soluable_residue) * hlp1 

    hlp1 = mortality_rate * pool_microbial%carbon14%ptr(ics:ice,:,iblk)
    loss_microbial%carbon14%ptr(ics:ice,:,iblk)    = loss_microbial%carbon14%ptr(ics:ice,:,iblk) + hlp1
    formation_dom%carbon14%ptr(ics:ice,:,iblk)     = formation_dom%carbon14%ptr(ics:ice,:,iblk) + & 
                                                     f_soluable_residue * hlp1 
    formation_residue%carbon14%ptr(ics:ice,:,iblk) = formation_dom%carbon14%ptr(ics:ice,:,iblk) + & 
                                                     (1._wp - f_soluable_residue) * hlp1 

    hlp1 = mortality_rate * pool_microbial%nitrogen15%ptr(ics:ice,:,iblk)
    loss_microbial%nitrogen15%ptr(ics:ice,:,iblk)    = loss_microbial%nitrogen15%ptr(ics:ice,:,iblk) + hlp1
    formation_dom%nitrogen15%ptr(ics:ice,:,iblk)     = formation_dom%nitrogen15%ptr(ics:ice,:,iblk) + & 
                                                       f_soluable_residue * hlp1 
    formation_residue%nitrogen15%ptr(ics:ice,:,iblk) = formation_dom%nitrogen15%ptr(ics:ice,:,iblk) + & 
                                                       (1._wp - f_soluable_residue) * hlp1 

  END SUBROUTINE calc_microbial_turnover

  !-----------------------------------------------------------------------------------------------------
  ! Sub Task called from update_sb_pools
  !
  !-----------------------------------------------------------------------------------------------------
  !> calculate growth of the microbial pool 
  !! corresponds to Ahrens et al. 2015, eq. 2 
  !! 
  !! OBS: currently assuming no nutrient constraint and infinite stoichiometric flexibility of mics
  !!
  !! Input: DOM and microbial pools, temperature and moisture rate modifiers 
  !!
  !! Output: DOM consumption, microbial growth, hetertrophic respiration and net mineralisation
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE calc_microbial_growth(options,nsoil,pool_dom,pool_microbial, & 
                                   rtm_mic_uptake,rmm_mic_uptake,rtm_hsc,rmm_hsc, &
                                   loss_dom,formation_microbial, &
                                   het_respiration,het_respiration_c13,het_respiration_c14, &
                                   net_mineralisation_n,net_mineralisation_n15,net_mineralisation_p)

    USE mo_sb_constants, ONLY: vmax_mic_uptake,km_mic_uptake,microbial_cue,microbial_cn,microbial_np

    IMPLICIT NONE

    ! 0.1 InOut
    TYPE(t_jsb_task_options),              INTENT(in)    :: options             !< model options
    INTEGER,                               INTENT(in)    :: nsoil               !< number of soil layers 
    TYPE(t_slm_struct_sb_element),         INTENT(in)    :: pool_dom            !< DOM [mol/m3]
    TYPE(t_slm_struct_sb_element),         INTENT(in)    :: pool_microbial      !< microbial biomass [mol/m3]
    REAL(wp), DIMENSION(options%nc,nsoil), INTENT(in)    :: rtm_mic_uptake      !< temperature rate modifier [-]
    REAL(wp), DIMENSION(options%nc,nsoil), INTENT(in)    :: rmm_mic_uptake      !< moisture rate modifier [-]
    REAL(wp), DIMENSION(options%nc,nsoil), INTENT(in)    :: rtm_hsc             !< temperature rate modifier [-]
    REAL(wp), DIMENSION(options%nc,nsoil), INTENT(in)    :: rmm_hsc             !< moisture rate modifier [-]
    TYPE(t_slm_struct_sb_element),         INTENT(inout) :: loss_dom            !< microbial death [mol/m3/timestep]
    TYPE(t_slm_struct_sb_element),         INTENT(inout) :: formation_microbial !< DOM formation from microbial death [mol/m3/timestep]
    REAL(wp), DIMENSION(options%nc,nsoil), INTENT(inout) :: het_respiration     !< respiration rate [mol C/m3/timestep]
    REAL(wp), DIMENSION(options%nc,nsoil), INTENT(inout) :: het_respiration_c13 !< respiration rate [mol 13C/m3/timestep] 
    REAL(wp), DIMENSION(options%nc,nsoil), INTENT(inout) :: het_respiration_c14 !< respiration rate [mol 14C/m3/timestep]
    REAL(wp), DIMENSION(options%nc,nsoil), INTENT(inout) :: net_mineralisation_n !< net mineralisation rate [mol N/m3/timestep]
    REAL(wp), DIMENSION(options%nc,nsoil), INTENT(inout) :: net_mineralisation_n15 !< net mineralisation rate [mol 15N/m3/timestep]
    REAL(wp), DIMENSION(options%nc,nsoil), INTENT(inout) :: net_mineralisation_p !< net mineralisation rate [mol P/m3/timestep]
    ! 0.2 Local 
    INTEGER                               :: iblk, ics, ice, nc
    REAL(wp), DIMENSION(options%nc,nsoil) :: hlp1, uptake_rate, f_uptake, uptake_rate_org, uptake_rate_min
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_microbial_growth'

    iblk  = options%iblk
    ics   = options%ics
    ice   = options%ice
    nc    = options%nc

    !> X.X microbial growth
    !! growth on DOM
    hlp1 = 1._wp   
    hlp1 = pool_dom%carbon%ptr(ics:ice,:,iblk) / (km_mic_uptake * rtm_hsc * rmm_hsc + pool_dom%carbon%ptr(ics:ice,:,iblk) )
    uptake_rate(:,:) = vmax_mic_uptake * rtm_mic_uptake(:,:) * rmm_mic_uptake(:,:) * hlp1(:,:) * dtime * & 
                       pool_microbial%carbon%ptr(ics:ice,:,iblk) 

    !> X.X nutrient limitation stuff ...

    !> X.X final rate 
    WHERE(pool_dom%carbon%ptr(ics:ice,:,iblk) > eps8)
       f_uptake(:,:) = uptake_rate(:,:) / pool_dom%carbon%ptr(ics:ice,:,iblk)
    ELSEWHERE
       f_uptake(:,:) = 0.0_wp
    ENDWHERE

    !> X.X update fluxes
    !! carbon
    uptake_rate_org(:,:) = f_uptake(:,:) * pool_dom%carbon%ptr(ics:ice,:,iblk)
    loss_dom%carbon%ptr(ics:ice,:,iblk)            = loss_dom%carbon%ptr(ics:ice,:,iblk) + & 
                                                     f_uptake(:,:) * pool_dom%carbon%ptr(ics:ice,:,iblk)
    formation_microbial%carbon%ptr(ics:ice,:,iblk) = formation_microbial%carbon%ptr(ics:ice,:,iblk) + & 
                                                     microbial_cue * uptake_rate_org 
    het_respiration(:,:)                           = het_respiration(:,:) + (1._wp - microbial_cue ) * uptake_rate_org

    !! nitrogen
    uptake_rate_min(:,:) = microbial_cue * uptake_rate_org / microbial_cn - f_uptake(:,:) * pool_dom%nitrogen%ptr(ics:ice,:,iblk) 
    uptake_rate_org(:,:) = f_uptake(:,:) * pool_dom%nitrogen%ptr(ics:ice,:,iblk)
    loss_dom%nitrogen%ptr(ics:ice,:,iblk)            = loss_dom%nitrogen%ptr(ics:ice,:,iblk) + & 
                                                       f_uptake(:,:) * pool_dom%nitrogen%ptr(ics:ice,:,iblk)
    formation_microbial%nitrogen%ptr(ics:ice,:,iblk) = formation_microbial%nitrogen%ptr(ics:ice,:,iblk) + & 
                                                       uptake_rate_org(:,:) + uptake_rate_min(:,:)
    net_mineralisation_n(:,:)                        = net_mineralisation_n(:,:) - uptake_rate_min(:,:) 
   
    !! phosphorus
    uptake_rate_min(:,:) = uptake_rate_org / microbial_np - f_uptake(:,:) * pool_dom%nitrogen%ptr(ics:ice,:,iblk) 
    uptake_rate_org(:,:) = f_uptake(:,:) * pool_dom%phosphorus%ptr(ics:ice,:,iblk)
    loss_dom%phosphorus%ptr(ics:ice,:,iblk)            = loss_dom%phosphorus%ptr(ics:ice,:,iblk) + & 
                                                         f_uptake(:,:) * pool_dom%phosphorus%ptr(ics:ice,:,iblk)
    formation_microbial%phosphorus%ptr(ics:ice,:,iblk) = formation_microbial%phosphorus%ptr(ics:ice,:,iblk) + & 
                                                         uptake_rate_org(:,:) + uptake_rate_min(:,:)
    net_mineralisation_p(:,:)                          = net_mineralisation_p(:,:) - uptake_rate_min(:,:) 
   
    !! carbon13
    uptake_rate_org(:,:) = f_uptake(:,:) * pool_dom%carbon13%ptr(ics:ice,:,iblk)
    loss_dom%carbon13%ptr(ics:ice,:,iblk)            = loss_dom%carbon13%ptr(ics:ice,:,iblk) + & 
                                                      f_uptake(:,:) * pool_dom%carbon13%ptr(ics:ice,:,iblk)
    formation_microbial%carbon13%ptr(ics:ice,:,iblk) = formation_microbial%carbon13%ptr(ics:ice,:,iblk) + & 
                                                       microbial_cue * uptake_rate_org 
    het_respiration_c13(:,:)                         = het_respiration_c13(:,:) + (1._wp - microbial_cue ) * uptake_rate_org

    !! carbon14
    uptake_rate_org(:,:) = f_uptake(:,:) * pool_dom%carbon14%ptr(ics:ice,:,iblk)
    loss_dom%carbon14%ptr(ics:ice,:,iblk)            = loss_dom%carbon14%ptr(ics:ice,:,iblk) + & 
                                                       f_uptake(:,:) * pool_dom%carbon14%ptr(ics:ice,:,iblk)
    formation_microbial%carbon14%ptr(ics:ice,:,iblk) = formation_microbial%carbon14%ptr(ics:ice,:,iblk) + & 
                                                       microbial_cue * uptake_rate_org 
    het_respiration_c14(:,:)                         = het_respiration_c14(:,:) + (1._wp - microbial_cue ) * uptake_rate_org

    !! nitrogen15
    uptake_rate_org(:,:) = f_uptake(:,:) * pool_dom%nitrogen15%ptr(ics:ice,:,iblk)
    loss_dom%nitrogen15%ptr(ics:ice,:,iblk)            = loss_dom%nitrogen15%ptr(ics:ice,:,iblk) + & 
                                                         f_uptake(:,:) * pool_dom%nitrogen15%ptr(ics:ice,:,iblk)
    formation_microbial%nitrogen15%ptr(ics:ice,:,iblk) = formation_microbial%nitrogen15%ptr(ics:ice,:,iblk) + & 
                                                         uptake_rate_org 
    net_mineralisation_n15(:,:)                        = net_mineralisation_n15(:,:) 
   
    !! phosphorus

  END SUBROUTINE calc_microbial_growth

  !-----------------------------------------------------------------------------------------------------
  ! Sub Task called from update_sb_pools
  !
  !-----------------------------------------------------------------------------------------------------
  !> calculate depolymerisation rate of polymeric litter and microbial residue 
  !! corresponds Ahrens et al. 2015, eq. 1 
  !! updated to the more recent version of JSM, including microbial residue as alternative polymeric pool
  !! 
  !! Input: maximum depolymerisation rate (vmax), temperature and moisture modifiers for 
  !!        depolymerisation and half-saturation, polymeric/residue pool and microbial biomass 
  !!
  !! Output: transfer of residue/polymeric pool to DOM 
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE calc_depolymerisation(options,nsoil,vmax_depolymerisation,rtm_depolymerisation,rmm_depolymerisation, & 
                                   rtm_hsc,rmm_hsc,microbial_carbon,pool_polymeric,loss_polymeric,formation_dom)

    USE mo_sb_constants,          ONLY: km_depolymerisation

    IMPLICIT NONE

    ! 0.1 InOut
    TYPE(t_jsb_task_options),              INTENT(in)    :: options               !< model options
    INTEGER,                               INTENT(in)    :: nsoil                 !< number of soil layers 
    REAL(wp),                              INTENT(in)    :: vmax_depolymerisation !< maximum depolymerisation rate
    REAL(wp), DIMENSION(options%nc,nsoil), INTENT(in)    :: rtm_depolymerisation  !< temperature rate modifier for depolymerisation
    REAL(wp), DIMENSION(options%nc,nsoil), INTENT(in)    :: rmm_depolymerisation  !< moisture rate modifier for depolymerisation
    REAL(wp), DIMENSION(options%nc,nsoil), INTENT(in)    :: rtm_hsc               !< temperature rate modifier for half-saturation constants
    REAL(wp), DIMENSION(options%nc,nsoil), INTENT(in)    :: rmm_hsc               !< moisture rate modifier for half-saturation constants
    REAL(wp), DIMENSION(options%nc,nsoil), INTENT(in)    :: microbial_carbon      !< microbial biomass carbon [mol C/m3]
    TYPE(t_slm_struct_sb_element),         INTENT(in)    :: pool_polymeric        !< polymeric/residue pool [mol/m3]
    TYPE(t_slm_struct_sb_element),         INTENT(inout) :: loss_polymeric        !< loss of polymeric/residue pool [mol/m3/timestep]
    TYPE(t_slm_struct_sb_element),         INTENT(inout) :: formation_dom         !< formation of DOM [mol/m3/timestep]
    ! 0.2 Local 
    INTEGER                               :: iblk, ics, ice, nc
    REAL(wp), DIMENSION(options%nc,nsoil) :: hlp1,hlp2,depolymerisation_rate
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_depolymerisation'

    iblk  = options%iblk
    ics   = options%ics
    ice   = options%ice
    nc    = options%nc

    !> 1.0 determine potential depolymerisation rate [mol / m3 / timestep]
    !!
    hlp1 = 1._wp   
    hlp1 = microbial_carbon(:,:) / (km_depolymerisation * rtm_hsc(:,:) * rmm_hsc(:,:) + microbial_carbon(:,:)) 
    hlp2 = vmax_depolymerisation * rtm_depolymerisation * rmm_depolymerisation * hlp1 * dtime * & 
           pool_polymeric%carbon%ptr(ics:ice,:,iblk) 
    
    !> 2.0 convert to turnover rate [1 / timestep] 
    !!
    WHERE(pool_polymeric%carbon%ptr(ics:ice,:,iblk) > eps8)
        depolymerisation_rate = hlp2 / pool_polymeric%carbon%ptr(ics:ice,:,iblk)
    ELSEWHERE
        depolymerisation_rate = 0.0_wp
    ENDWHERE

    !> 3.0 apply to all elements in this pool (should be a simple operator overload, but does currently not work safely)
    !!
    hlp2 = depolymerisation_rate * pool_polymeric%carbon%ptr(ics:ice,:,iblk)
    loss_polymeric%carbon%ptr(ics:ice,:,iblk) = loss_polymeric%carbon%ptr(ics:ice,:,iblk) + hlp2 
    formation_dom%carbon%ptr(ics:ice,:,iblk) = formation_dom%carbon%ptr(ics:ice,:,iblk) + hlp2 

    hlp2 = depolymerisation_rate * pool_polymeric%nitrogen%ptr(ics:ice,:,iblk)
    loss_polymeric%nitrogen%ptr(ics:ice,:,iblk) = loss_polymeric%nitrogen%ptr(ics:ice,:,iblk) + hlp2 
    formation_dom%nitrogen%ptr(ics:ice,:,iblk) = formation_dom%nitrogen%ptr(ics:ice,:,iblk) + hlp2 

    hlp2 = depolymerisation_rate * pool_polymeric%phosphorus%ptr(ics:ice,:,iblk)
    loss_polymeric%phosphorus%ptr(ics:ice,:,iblk) = loss_polymeric%phosphorus%ptr(ics:ice,:,iblk) + hlp2 
    formation_dom%phosphorus%ptr(ics:ice,:,iblk) = formation_dom%phosphorus%ptr(ics:ice,:,iblk) + hlp2 

    hlp2 = depolymerisation_rate * pool_polymeric%carbon13%ptr(ics:ice,:,iblk)
    loss_polymeric%carbon13%ptr(ics:ice,:,iblk) = loss_polymeric%carbon13%ptr(ics:ice,:,iblk) + hlp2 
    formation_dom%carbon13%ptr(ics:ice,:,iblk) = formation_dom%carbon13%ptr(ics:ice,:,iblk) + hlp2 

    hlp2 = depolymerisation_rate * pool_polymeric%carbon14%ptr(ics:ice,:,iblk)
    loss_polymeric%carbon14%ptr(ics:ice,:,iblk) = loss_polymeric%carbon14%ptr(ics:ice,:,iblk) + hlp2 
    formation_dom%carbon14%ptr(ics:ice,:,iblk) = formation_dom%carbon14%ptr(ics:ice,:,iblk) + hlp2 

    hlp2 = depolymerisation_rate * pool_polymeric%nitrogen15%ptr(ics:ice,:,iblk)
    loss_polymeric%nitrogen15%ptr(ics:ice,:,iblk) = loss_polymeric%nitrogen15%ptr(ics:ice,:,iblk) + hlp2 
    formation_dom%nitrogen15%ptr(ics:ice,:,iblk) = formation_dom%nitrogen15%ptr(ics:ice,:,iblk) + hlp2 

  END SUBROUTINE calc_depolymerisation

  !-----------------------------------------------------------------------------------------------------
  ! Sub Task called from update_sb_pools 
  !
  !-----------------------------------------------------------------------------------------------------
  !> calculates the bulk density of the soil corrected by its organic component 
  !! derived from Federer et al. 1993, Can J. For. Res. 23 1026-1032 
  !!
  !! Input: carbon content by layer 
  !!
  !! Output: corrected bulk density of the soil 
  !-----------------------------------------------------------------------------------------------------
  ELEMENTAL FUNCTION calc_bulk_density_correction(bulk_soil_carbon_sl,bulk_dens_l) RESULT(bulk_dens_corr_sl)

    USE mo_sb_constants,          ONLY: rho_bulk_org 
    USE mo_shared_constants,      ONLY: molar_mass_C 

    IMPLICIT NONE

    ! 0.1 InOut
    REAL(wp), INTENT(in) :: bulk_soil_carbon_sl     !< soil carbon contributing to density [mol C / m3] 
    REAL(wp), INTENT(in) :: bulk_dens_l             !< bulk density of the mineral soil [kg / m3]
    REAL(wp)             :: bulk_dens_corr_sl       !< combined bulk density of the mineral and organic soil [kg / m3]
    ! 0.2 Local
    REAL(wp) :: fract_org
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_bulk_density_correction'
 
    !> 1.0 calculate organic mass fraction per soil layer 
    !!
    IF(bulk_soil_carbon_sl * molar_mass_C + bulk_dens_l > eps8)THEN
       fract_org = bulk_soil_carbon_sl * molar_mass_C/1000._wp / & 
                   (bulk_soil_carbon_sl * molar_mass_C/1000._wp + bulk_dens_l)
    ELSE 
       fract_org = 0.0_wp
    ENDIF

    !> 2.0 actual bulk density of the soil 
    !!
    bulk_dens_corr_sl = bulk_dens_l * rho_bulk_org / & 
                          (fract_org * bulk_dens_l + (1._wp - fract_org ) * rho_bulk_org)

  END FUNCTION calc_bulk_density_correction 

  !-----------------------------------------------------------------------------------------------------
  ! Sub Task called from update_sb_pools 
  !
  !-----------------------------------------------------------------------------------------------------
  !> calculates the maximum sorption capacity for organic material of a soil layer
  !! derived from COMISSION, the implict assumption here is that no association happens on organic material
  !!
  !! Input: bulk density of the soil (mineral and organic combined) 
  !!
  !! Output: sorption capacity for organic material 
  !-----------------------------------------------------------------------------------------------------
  ELEMENTAL FUNCTION calc_qmax_org(bulk_dens_corr_sl) RESULT(qmax_org)

    USE mo_sb_constants,          ONLY: qmax_org_mineral_default 

    IMPLICIT NONE

    ! 0.1 InOut
    REAL(wp), INTENT(in) :: bulk_dens_corr_sl     !< combined bulk density of the mineral and organic soil [kg / m3] 
    REAL(wp)             :: qmax_org      !< maximum sorption capacity of the soil [kg C / m3]
    ! 0.2 Local
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_qmax_org'
   
    qmax_org = qmax_org_mineral_default * bulk_dens_corr_sl

  END FUNCTION calc_qmax_org

  !-----------------------------------------------------------------------------------------------------
  ! Sub Task called from update_sb_pools 
  !
  !-----------------------------------------------------------------------------------------------------
  !> calculates bulk soil carbon used for the calculation of soil bulk density and sorption capacity 
  !!
  !! Input: soil pools 
  !!
  !! Output: total_bulk_soil_carbon 
  !-----------------------------------------------------------------------------------------------------
  FUNCTION calc_bulk_soil_carbon(options,nsoil,sb_pool) RESULT (bulk_soil_carbon_sl)

    IMPLICIT NONE

    ! 0.1 InOut
    TYPE(t_jsb_task_options),              INTENT(in)    :: options
    INTEGER,                               INTENT(in)    :: nsoil
    TYPE(t_slm_struct_sb_pool),            INTENT(in)    :: sb_pool                 !< pools [mol/m3]
    REAL(wp), DIMENSION(options%nc,nsoil)                :: bulk_soil_carbon_sl
    ! 0.2 Local
    INTEGER                     :: iblk, ics, ice, nc
    INTEGER                     :: isoil
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_bulk_soil_carbon'

    iblk  = options%iblk
    ics   = options%ics
    ice   = options%ice

    DO isoil = 1,nsoil
       bulk_soil_carbon_sl(:,isoil) = &
            sb_pool%microbial%carbon%ptr(ics:ice,isoil,iblk) + &
            sb_pool%residue%carbon%ptr(ics:ice,isoil,iblk) + &
            sb_pool%residue_assoc%carbon%ptr(ics:ice,isoil,iblk) + &
            sb_pool%dom%carbon%ptr(ics:ice,isoil,iblk) + &
            sb_pool%dom_assoc%carbon%ptr(ics:ice,isoil,iblk) + &
            sb_pool%soluable_litter%carbon%ptr(ics:ice,isoil,iblk) + &
            sb_pool%polymeric_litter%carbon%ptr(ics:ice,isoil,iblk)
    !   IF(isoil > 1 ) THEN       
          bulk_soil_carbon_sl(:,isoil) = bulk_soil_carbon_sl(:,isoil) + &
               sb_pool%woody_litter%carbon%ptr(ics:ice,isoil,iblk) 
    !   ENDIF
    ENDDO

  END FUNCTION calc_bulk_soil_carbon

  !-----------------------------------------------------------------------------------------------------
  ! Sub Task called from update_sb_pools 
  !
  !-----------------------------------------------------------------------------------------------------
  !> calculates sorption and desorption of organic material 
  !! corresponds to the sorption and desorption terms in Ahrens et al. 2015, eq. 5 & 6 
  !!
  !! Input: sorption temperature and moisture modifiers, sorption pools
  !!
  !! Output: sorption and desorption rate (mol/m2/timestep) 
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE calc_om_sorption(options,nsoil,qmax_org,rtm_sorption,rmm_sorption,rtm_desorption, &
                              rmm_desorption,pool_dom,pool_dom_assoc,pool_residue, & 
                              pool_residue_assoc,formation_dom,formation_dom_assoc,formation_residue, &
                              formation_residue_assoc,loss_dom,loss_dom_assoc,loss_residue,loss_residue_assoc)

    USE mo_sb_constants,          ONLY: k_adsorpt_dom,k_desorpt_dom,k_adsorpt_det,k_desorpt_det 

    IMPLICIT NONE

    ! 0.1 InOut
    TYPE(t_jsb_task_options),              INTENT(in)    :: options
    INTEGER,                               INTENT(in)    :: nsoil
    REAL(wp), DIMENSION(options%nc,nsoil), INTENT(in)    :: qmax_org 
    REAL(wp), DIMENSION(options%nc,nsoil), INTENT(in)    :: rtm_sorption            !< temperature rate modifier [-]
    REAL(wp), DIMENSION(options%nc,nsoil), INTENT(in)    :: rmm_sorption            !< moisture rate modifier [-]
    REAL(wp), DIMENSION(options%nc,nsoil), INTENT(in)    :: rtm_desorption          !< temperature rate modifier [-]
    REAL(wp), DIMENSION(options%nc,nsoil), INTENT(in)    :: rmm_desorption          !< moisture rate modifier [-]
    TYPE(t_slm_struct_sb_element),         INTENT(in)    :: pool_dom                !< DOM pool [mol/m3]
    TYPE(t_slm_struct_sb_element),         INTENT(in)    :: pool_residue            !< DOM pool [mol/m3]
    TYPE(t_slm_struct_sb_element),         INTENT(in)    :: pool_dom_assoc          !< DOM pool [mol/m3]
    TYPE(t_slm_struct_sb_element),         INTENT(in)    :: pool_residue_assoc      !< DOM pool [mol/m3]
    TYPE(t_slm_struct_sb_element),         INTENT(inout) :: formation_dom           !< rate of DOM transport [mol/m3/timestep]
    TYPE(t_slm_struct_sb_element),         INTENT(inout) :: formation_dom_assoc     !< rate of DOM transport [mol/m3/timestep]
    TYPE(t_slm_struct_sb_element),         INTENT(inout) :: formation_residue       !< rate of DOM transport [mol/m3/timestep]
    TYPE(t_slm_struct_sb_element),         INTENT(inout) :: formation_residue_assoc !< rate of DOM transport [mol/m3/timestep]
    TYPE(t_slm_struct_sb_element),         INTENT(inout) :: loss_dom                !< rate of DOM transport [mol/m3/timestep]
    TYPE(t_slm_struct_sb_element),         INTENT(inout) :: loss_dom_assoc          !< rate of DOM transport [mol/m3/timestep]
    TYPE(t_slm_struct_sb_element),         INTENT(inout) :: loss_residue            !< rate of DOM transport [mol/m3/timestep]
    TYPE(t_slm_struct_sb_element),         INTENT(inout) :: loss_residue_assoc      !< rate of DOM transport [mol/m3/timestep]
    ! 0.2 Local
    REAL(wp), DIMENSION(options%nc,nsoil) :: hlp1,avail_sorption_sites,sorption_rate,desorption_rate
    INTEGER                               :: iblk, ics, ice, nc
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_om_sorption'

    iblk  = options%iblk
    ics   = options%ics
    ice   = options%ice
    nc    = options%nc

    !> 1.0 Determine available sorption sites [mol C m-3]
    !!
    avail_sorption_sites(:,:) = qmax_org(:,:) - & 
                                (pool_dom_assoc%carbon%ptr(ics:ice,:,iblk) + pool_residue_assoc%carbon%ptr(ics:ice,:,iblk))

    !> 2.0 mineral association of DOM
    sorption_rate(:,:)   = k_adsorpt_dom * rtm_sorption(:,:)   * rmm_sorption(:,:) * avail_sorption_sites(:,:) 
    desorption_rate(:,:) = k_desorpt_dom * rtm_desorption(:,:) * rmm_desorption(:,:) 

    hlp1 = sorption_rate(:,:) * pool_dom%carbon%ptr(ics:ice,:,iblk)
    formation_dom_assoc%carbon%ptr(ics:ice,:,iblk) = formation_dom_assoc%carbon%ptr(ics:ice,:,iblk) + hlp1 
    loss_dom%carbon%ptr(ics:ice,:,iblk)            = loss_dom%carbon%ptr(ics:ice,:,iblk) + hlp1 
    hlp1 = desorption_rate(:,:) * pool_dom_assoc%carbon%ptr(ics:ice,:,iblk)
    formation_dom%carbon%ptr(ics:ice,:,iblk)       = formation_dom%carbon%ptr(ics:ice,:,iblk) + hlp1 
    loss_dom_assoc%carbon%ptr(ics:ice,:,iblk)      = loss_dom_assoc%carbon%ptr(ics:ice,:,iblk) + hlp1 

    hlp1 = sorption_rate(:,:) * pool_dom%nitrogen%ptr(ics:ice,:,iblk)
    formation_dom_assoc%nitrogen%ptr(ics:ice,:,iblk) = formation_dom_assoc%nitrogen%ptr(ics:ice,:,iblk) + hlp1 
    loss_dom%nitrogen%ptr(ics:ice,:,iblk)            = loss_dom%nitrogen%ptr(ics:ice,:,iblk) + hlp1 
    hlp1 = desorption_rate(:,:) * pool_dom_assoc%nitrogen%ptr(ics:ice,:,iblk)
    formation_dom%nitrogen%ptr(ics:ice,:,iblk)       = formation_dom%nitrogen%ptr(ics:ice,:,iblk) + hlp1 
    loss_dom_assoc%nitrogen%ptr(ics:ice,:,iblk)      = loss_dom_assoc%nitrogen%ptr(ics:ice,:,iblk) + hlp1 

    hlp1 = sorption_rate(:,:) * pool_dom%phosphorus%ptr(ics:ice,:,iblk)
    formation_dom_assoc%phosphorus%ptr(ics:ice,:,iblk) = formation_dom_assoc%phosphorus%ptr(ics:ice,:,iblk) + hlp1 
    loss_dom%phosphorus%ptr(ics:ice,:,iblk)            = loss_dom%phosphorus%ptr(ics:ice,:,iblk) + hlp1 
    hlp1 = desorption_rate(:,:) * pool_dom_assoc%phosphorus%ptr(ics:ice,:,iblk)
    formation_dom%phosphorus%ptr(ics:ice,:,iblk)       = formation_dom%phosphorus%ptr(ics:ice,:,iblk) + hlp1 
    loss_dom_assoc%phosphorus%ptr(ics:ice,:,iblk)      = loss_dom_assoc%phosphorus%ptr(ics:ice,:,iblk) + hlp1 

    hlp1 = sorption_rate(:,:) * pool_dom%carbon13%ptr(ics:ice,:,iblk)
    formation_dom_assoc%carbon13%ptr(ics:ice,:,iblk) = formation_dom_assoc%carbon13%ptr(ics:ice,:,iblk) + hlp1 
    loss_dom%carbon13%ptr(ics:ice,:,iblk)            = loss_dom%carbon13%ptr(ics:ice,:,iblk) + hlp1 
    hlp1 = desorption_rate(:,:) * pool_dom_assoc%carbon13%ptr(ics:ice,:,iblk)
    formation_dom%carbon13%ptr(ics:ice,:,iblk)       = formation_dom%carbon13%ptr(ics:ice,:,iblk) + hlp1 
    loss_dom_assoc%carbon13%ptr(ics:ice,:,iblk)      = loss_dom_assoc%carbon13%ptr(ics:ice,:,iblk) + hlp1 

    hlp1 = sorption_rate(:,:) * pool_dom%carbon14%ptr(ics:ice,:,iblk)
    formation_dom_assoc%carbon14%ptr(ics:ice,:,iblk) = formation_dom_assoc%carbon14%ptr(ics:ice,:,iblk) + hlp1 
    loss_dom%carbon14%ptr(ics:ice,:,iblk)            = loss_dom%carbon14%ptr(ics:ice,:,iblk) + hlp1 
    hlp1 = desorption_rate(:,:) * pool_dom_assoc%carbon14%ptr(ics:ice,:,iblk)
    formation_dom%carbon14%ptr(ics:ice,:,iblk)       = formation_dom%carbon14%ptr(ics:ice,:,iblk) + hlp1 
    loss_dom_assoc%carbon14%ptr(ics:ice,:,iblk)      = loss_dom_assoc%carbon14%ptr(ics:ice,:,iblk) + hlp1 

    hlp1 = sorption_rate(:,:) * pool_dom%nitrogen15%ptr(ics:ice,:,iblk)
    formation_dom_assoc%nitrogen15%ptr(ics:ice,:,iblk) = formation_dom_assoc%nitrogen15%ptr(ics:ice,:,iblk) + hlp1 
    loss_dom%nitrogen15%ptr(ics:ice,:,iblk)            = loss_dom%nitrogen15%ptr(ics:ice,:,iblk) + hlp1 
    hlp1 = desorption_rate(:,:) * pool_dom_assoc%nitrogen15%ptr(ics:ice,:,iblk)
    formation_dom%nitrogen15%ptr(ics:ice,:,iblk)       = formation_dom%nitrogen15%ptr(ics:ice,:,iblk) + hlp1 
    loss_dom_assoc%nitrogen15%ptr(ics:ice,:,iblk)      = loss_dom_assoc%nitrogen15%ptr(ics:ice,:,iblk) + hlp1 

    !> 3.0 mineral association of microbial residue 
    sorption_rate(:,:)   = k_adsorpt_det * rtm_sorption(:,:)   * rmm_sorption(:,:) * avail_sorption_sites(:,:) 
    desorption_rate(:,:) = k_desorpt_det * rtm_desorption(:,:) * rmm_desorption(:,:) 

    hlp1 = sorption_rate(:,:) * pool_residue%carbon%ptr(ics:ice,:,iblk)
    formation_residue_assoc%carbon%ptr(ics:ice,:,iblk) = formation_residue_assoc%carbon%ptr(ics:ice,:,iblk) + hlp1 
    loss_residue%carbon%ptr(ics:ice,:,iblk)            = loss_residue%carbon%ptr(ics:ice,:,iblk) + hlp1 
    hlp1 = desorption_rate(:,:) * pool_residue_assoc%carbon%ptr(ics:ice,:,iblk)
    formation_residue%carbon%ptr(ics:ice,:,iblk)       = formation_residue%carbon%ptr(ics:ice,:,iblk) + hlp1 
    loss_residue_assoc%carbon%ptr(ics:ice,:,iblk)      = loss_residue_assoc%carbon%ptr(ics:ice,:,iblk) + hlp1 

    hlp1 = sorption_rate(:,:) * pool_residue%nitrogen%ptr(ics:ice,:,iblk)
    formation_residue_assoc%nitrogen%ptr(ics:ice,:,iblk) = formation_residue_assoc%nitrogen%ptr(ics:ice,:,iblk) + hlp1 
    loss_residue%nitrogen%ptr(ics:ice,:,iblk)            = loss_residue%nitrogen%ptr(ics:ice,:,iblk) + hlp1 
    hlp1 = desorption_rate(:,:) * pool_residue_assoc%nitrogen%ptr(ics:ice,:,iblk)
    formation_residue%nitrogen%ptr(ics:ice,:,iblk)       = formation_residue%nitrogen%ptr(ics:ice,:,iblk) + hlp1 
    loss_residue_assoc%nitrogen%ptr(ics:ice,:,iblk)      = loss_residue_assoc%nitrogen%ptr(ics:ice,:,iblk) + hlp1 

    hlp1 = sorption_rate(:,:) * pool_residue%phosphorus%ptr(ics:ice,:,iblk)
    formation_residue_assoc%phosphorus%ptr(ics:ice,:,iblk) = formation_residue_assoc%phosphorus%ptr(ics:ice,:,iblk) + hlp1 
    loss_residue%phosphorus%ptr(ics:ice,:,iblk)            = loss_residue%phosphorus%ptr(ics:ice,:,iblk) + hlp1 
    hlp1 = desorption_rate(:,:) * pool_residue_assoc%phosphorus%ptr(ics:ice,:,iblk)
    formation_residue%phosphorus%ptr(ics:ice,:,iblk)       = formation_residue%phosphorus%ptr(ics:ice,:,iblk) + hlp1 
    loss_residue_assoc%phosphorus%ptr(ics:ice,:,iblk)      = loss_residue_assoc%phosphorus%ptr(ics:ice,:,iblk) + hlp1 

    hlp1 = sorption_rate(:,:) * pool_residue%carbon13%ptr(ics:ice,:,iblk)
    formation_residue_assoc%carbon13%ptr(ics:ice,:,iblk) = formation_residue_assoc%carbon13%ptr(ics:ice,:,iblk) + hlp1 
    loss_residue%carbon13%ptr(ics:ice,:,iblk)            = loss_residue%carbon13%ptr(ics:ice,:,iblk) + hlp1 
    hlp1 = desorption_rate(:,:) * pool_residue_assoc%carbon13%ptr(ics:ice,:,iblk)
    formation_residue%carbon13%ptr(ics:ice,:,iblk)       = formation_residue%carbon13%ptr(ics:ice,:,iblk) + hlp1 
    loss_residue_assoc%carbon13%ptr(ics:ice,:,iblk)      = loss_residue_assoc%carbon13%ptr(ics:ice,:,iblk) + hlp1 

    hlp1 = sorption_rate(:,:) * pool_residue%carbon14%ptr(ics:ice,:,iblk)
    formation_residue_assoc%carbon14%ptr(ics:ice,:,iblk) = formation_residue_assoc%carbon14%ptr(ics:ice,:,iblk) + hlp1 
    loss_residue%carbon14%ptr(ics:ice,:,iblk)            = loss_residue%carbon14%ptr(ics:ice,:,iblk) + hlp1 
    hlp1 = desorption_rate(:,:) * pool_residue_assoc%carbon14%ptr(ics:ice,:,iblk)
    formation_residue%carbon14%ptr(ics:ice,:,iblk)       = formation_residue%carbon14%ptr(ics:ice,:,iblk) + hlp1 
    loss_residue_assoc%carbon14%ptr(ics:ice,:,iblk)      = loss_residue_assoc%carbon14%ptr(ics:ice,:,iblk) + hlp1 

    hlp1 = sorption_rate(:,:) * pool_residue%nitrogen15%ptr(ics:ice,:,iblk)
    formation_residue_assoc%nitrogen15%ptr(ics:ice,:,iblk) = formation_residue_assoc%nitrogen15%ptr(ics:ice,:,iblk) + hlp1 
    loss_residue%nitrogen15%ptr(ics:ice,:,iblk)            = loss_residue%nitrogen15%ptr(ics:ice,:,iblk) + hlp1 
    hlp1 = desorption_rate(:,:) * pool_residue_assoc%nitrogen15%ptr(ics:ice,:,iblk)
    formation_residue%nitrogen15%ptr(ics:ice,:,iblk)       = formation_residue%nitrogen15%ptr(ics:ice,:,iblk) + hlp1 
    loss_residue_assoc%nitrogen15%ptr(ics:ice,:,iblk)      = loss_residue_assoc%nitrogen15%ptr(ics:ice,:,iblk) + hlp1 

  END SUBROUTINE calc_om_sorption 

  !-----------------------------------------------------------------------------------------------------
  ! Sub Task called from update_sb_pools
  !
  !-----------------------------------------------------------------------------------------------------
  !> updates all SB pools given the JSM equations for formation, loss and transport 
  !! should be a simple pool operator overload function, but that doesn't work yet
  !! 
  !! Input: formation, loss and transport rate, pool sizes 
  !!
  !! Output: updated pool sizes 
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE calc_jsm_pool_change(options,pool,formation,loss,transport)

    IMPLICIT NONE

    ! 0.1 InOut
    TYPE(t_jsb_task_options), INTENT(in)          :: options      !< model options
    TYPE(t_slm_struct_sb_pool),  INTENT(in)       :: formation    !< rate of formation of a pool [mol/m3/timestep]
    TYPE(t_slm_struct_sb_pool),  INTENT(in)       :: loss         !< rate of loss of a pool [mol/m3/timestep]
    TYPE(t_slm_struct_sb_pool),  INTENT(in)       :: transport    !< rate of transport of a pool [mol/m3/timestep]
    TYPE(t_slm_struct_sb_pool),  INTENT(inout)    :: pool         !< the pool [mol/m3]
    ! 0.2 Local 
    INTEGER                               :: iblk, ics, ice, nc
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_jsm_pool_change'

    iblk  = options%iblk
    ics   = options%ics
    ice   = options%ice
    nc    = options%nc

    ! carbon
    pool%dom%carbon%ptr(ics:ice,:,iblk) = pool%dom%carbon%ptr(ics:ice,:,iblk) + & 
                                                formation%dom%carbon%ptr(ics:ice,:,iblk) - & 
                                                loss%dom%carbon%ptr(ics:ice,:,iblk) + & 
                                                transport%dom%carbon%ptr(ics:ice,:,iblk)  
    pool%dom_assoc%carbon%ptr(ics:ice,:,iblk) = pool%dom_assoc%carbon%ptr(ics:ice,:,iblk) + & 
                                                formation%dom_assoc%carbon%ptr(ics:ice,:,iblk) - & 
                                                loss%dom_assoc%carbon%ptr(ics:ice,:,iblk) + & 
                                                transport%dom_assoc%carbon%ptr(ics:ice,:,iblk)  
    pool%soluable_litter%carbon%ptr(ics:ice,:,iblk) = pool%soluable_litter%carbon%ptr(ics:ice,:,iblk) + & 
                                                formation%soluable_litter%carbon%ptr(ics:ice,:,iblk) - & 
                                                loss%soluable_litter%carbon%ptr(ics:ice,:,iblk) + & 
                                                transport%soluable_litter%carbon%ptr(ics:ice,:,iblk)  
    pool%polymeric_litter%carbon%ptr(ics:ice,:,iblk) = pool%polymeric_litter%carbon%ptr(ics:ice,:,iblk) + & 
                                                formation%polymeric_litter%carbon%ptr(ics:ice,:,iblk) - & 
                                                loss%polymeric_litter%carbon%ptr(ics:ice,:,iblk) + & 
                                                transport%polymeric_litter%carbon%ptr(ics:ice,:,iblk)  
    pool%woody_litter%carbon%ptr(ics:ice,:,iblk) = pool%woody_litter%carbon%ptr(ics:ice,:,iblk) + & 
                                                formation%woody_litter%carbon%ptr(ics:ice,:,iblk) - & 
                                                loss%woody_litter%carbon%ptr(ics:ice,:,iblk) + & 
                                                transport%woody_litter%carbon%ptr(ics:ice,:,iblk)  
    pool%microbial%carbon%ptr(ics:ice,:,iblk) = pool%microbial%carbon%ptr(ics:ice,:,iblk) + & 
                                                formation%microbial%carbon%ptr(ics:ice,:,iblk) - & 
                                                loss%microbial%carbon%ptr(ics:ice,:,iblk) + & 
                                                transport%microbial%carbon%ptr(ics:ice,:,iblk)  
    pool%residue%carbon%ptr(ics:ice,:,iblk) = pool%residue%carbon%ptr(ics:ice,:,iblk) + & 
                                                formation%residue%carbon%ptr(ics:ice,:,iblk) - & 
                                                loss%residue%carbon%ptr(ics:ice,:,iblk) + & 
                                                transport%residue%carbon%ptr(ics:ice,:,iblk)  
    pool%residue_assoc%carbon%ptr(ics:ice,:,iblk) = pool%residue_assoc%carbon%ptr(ics:ice,:,iblk) + & 
                                                formation%residue_assoc%carbon%ptr(ics:ice,:,iblk) - & 
                                                loss%residue_assoc%carbon%ptr(ics:ice,:,iblk) + & 
                                                transport%residue_assoc%carbon%ptr(ics:ice,:,iblk)  

    ! nitrogen
    pool%dom%nitrogen%ptr(ics:ice,:,iblk) = pool%dom%nitrogen%ptr(ics:ice,:,iblk) + & 
                                                formation%dom%nitrogen%ptr(ics:ice,:,iblk) - & 
                                                loss%dom%nitrogen%ptr(ics:ice,:,iblk) + & 
                                                transport%dom%nitrogen%ptr(ics:ice,:,iblk)  
    pool%dom_assoc%nitrogen%ptr(ics:ice,:,iblk) = pool%dom_assoc%nitrogen%ptr(ics:ice,:,iblk) + & 
                                                formation%dom_assoc%nitrogen%ptr(ics:ice,:,iblk) - & 
                                                loss%dom_assoc%nitrogen%ptr(ics:ice,:,iblk) + & 
                                                transport%dom_assoc%nitrogen%ptr(ics:ice,:,iblk)  
    pool%soluable_litter%nitrogen%ptr(ics:ice,:,iblk) = pool%soluable_litter%nitrogen%ptr(ics:ice,:,iblk) + & 
                                                formation%soluable_litter%nitrogen%ptr(ics:ice,:,iblk) - & 
                                                loss%soluable_litter%nitrogen%ptr(ics:ice,:,iblk) + & 
                                                transport%soluable_litter%nitrogen%ptr(ics:ice,:,iblk)  
    pool%polymeric_litter%nitrogen%ptr(ics:ice,:,iblk) = pool%polymeric_litter%nitrogen%ptr(ics:ice,:,iblk) + & 
                                                formation%polymeric_litter%nitrogen%ptr(ics:ice,:,iblk) - & 
                                                loss%polymeric_litter%nitrogen%ptr(ics:ice,:,iblk) + & 
                                                transport%polymeric_litter%nitrogen%ptr(ics:ice,:,iblk)  
    pool%woody_litter%nitrogen%ptr(ics:ice,:,iblk) = pool%woody_litter%nitrogen%ptr(ics:ice,:,iblk) + & 
                                                formation%woody_litter%nitrogen%ptr(ics:ice,:,iblk) - & 
                                                loss%woody_litter%nitrogen%ptr(ics:ice,:,iblk) + & 
                                                transport%woody_litter%nitrogen%ptr(ics:ice,:,iblk)  
    pool%microbial%nitrogen%ptr(ics:ice,:,iblk) = pool%microbial%nitrogen%ptr(ics:ice,:,iblk) + & 
                                                formation%microbial%nitrogen%ptr(ics:ice,:,iblk) - & 
                                                loss%microbial%nitrogen%ptr(ics:ice,:,iblk) + & 
                                                transport%microbial%nitrogen%ptr(ics:ice,:,iblk)  
    pool%residue%nitrogen%ptr(ics:ice,:,iblk) = pool%residue%nitrogen%ptr(ics:ice,:,iblk) + & 
                                                formation%residue%nitrogen%ptr(ics:ice,:,iblk) - & 
                                                loss%residue%nitrogen%ptr(ics:ice,:,iblk) + & 
                                                transport%residue%nitrogen%ptr(ics:ice,:,iblk)  
    pool%residue_assoc%nitrogen%ptr(ics:ice,:,iblk) = pool%residue_assoc%nitrogen%ptr(ics:ice,:,iblk) + & 
                                                formation%residue_assoc%nitrogen%ptr(ics:ice,:,iblk) - & 
                                                loss%residue_assoc%nitrogen%ptr(ics:ice,:,iblk) + & 
                                                transport%residue_assoc%nitrogen%ptr(ics:ice,:,iblk)  

    ! phosphorus
    pool%dom%phosphorus%ptr(ics:ice,:,iblk) = pool%dom%phosphorus%ptr(ics:ice,:,iblk) + & 
                                                formation%dom%phosphorus%ptr(ics:ice,:,iblk) - & 
                                                loss%dom%phosphorus%ptr(ics:ice,:,iblk) + & 
                                                transport%dom%phosphorus%ptr(ics:ice,:,iblk)  
    pool%dom_assoc%phosphorus%ptr(ics:ice,:,iblk) = pool%dom_assoc%phosphorus%ptr(ics:ice,:,iblk) + & 
                                                formation%dom_assoc%phosphorus%ptr(ics:ice,:,iblk) - & 
                                                loss%dom_assoc%phosphorus%ptr(ics:ice,:,iblk) + & 
                                                transport%dom_assoc%phosphorus%ptr(ics:ice,:,iblk)  
    pool%soluable_litter%phosphorus%ptr(ics:ice,:,iblk) = pool%soluable_litter%phosphorus%ptr(ics:ice,:,iblk) + & 
                                                formation%soluable_litter%phosphorus%ptr(ics:ice,:,iblk) - & 
                                                loss%soluable_litter%phosphorus%ptr(ics:ice,:,iblk) + & 
                                                transport%soluable_litter%phosphorus%ptr(ics:ice,:,iblk)  
    pool%polymeric_litter%phosphorus%ptr(ics:ice,:,iblk) = pool%polymeric_litter%phosphorus%ptr(ics:ice,:,iblk) + & 
                                                formation%polymeric_litter%phosphorus%ptr(ics:ice,:,iblk) - & 
                                                loss%polymeric_litter%phosphorus%ptr(ics:ice,:,iblk) + & 
                                                transport%polymeric_litter%phosphorus%ptr(ics:ice,:,iblk)  
    pool%woody_litter%phosphorus%ptr(ics:ice,:,iblk) = pool%woody_litter%phosphorus%ptr(ics:ice,:,iblk) + & 
                                                formation%woody_litter%phosphorus%ptr(ics:ice,:,iblk) - & 
                                                loss%woody_litter%phosphorus%ptr(ics:ice,:,iblk) + & 
                                                transport%woody_litter%phosphorus%ptr(ics:ice,:,iblk)  
    pool%microbial%phosphorus%ptr(ics:ice,:,iblk) = pool%microbial%phosphorus%ptr(ics:ice,:,iblk) + & 
                                                formation%microbial%phosphorus%ptr(ics:ice,:,iblk) - & 
                                                loss%microbial%phosphorus%ptr(ics:ice,:,iblk) + & 
                                                transport%microbial%phosphorus%ptr(ics:ice,:,iblk)  
    pool%residue%phosphorus%ptr(ics:ice,:,iblk) = pool%residue%phosphorus%ptr(ics:ice,:,iblk) + & 
                                                formation%residue%phosphorus%ptr(ics:ice,:,iblk) - & 
                                                loss%residue%phosphorus%ptr(ics:ice,:,iblk) + & 
                                                transport%residue%phosphorus%ptr(ics:ice,:,iblk)  
    pool%residue_assoc%phosphorus%ptr(ics:ice,:,iblk) = pool%residue_assoc%phosphorus%ptr(ics:ice,:,iblk) + & 
                                                formation%residue_assoc%phosphorus%ptr(ics:ice,:,iblk) - & 
                                                loss%residue_assoc%phosphorus%ptr(ics:ice,:,iblk) + & 
                                                transport%residue_assoc%phosphorus%ptr(ics:ice,:,iblk)  

    ! C13
    pool%dom%carbon13%ptr(ics:ice,:,iblk) = pool%dom%carbon13%ptr(ics:ice,:,iblk) + & 
                                                formation%dom%carbon13%ptr(ics:ice,:,iblk) - & 
                                                loss%dom%carbon13%ptr(ics:ice,:,iblk) + & 
                                                transport%dom%carbon13%ptr(ics:ice,:,iblk)  
    pool%dom_assoc%carbon13%ptr(ics:ice,:,iblk) = pool%dom_assoc%carbon13%ptr(ics:ice,:,iblk) + & 
                                                formation%dom_assoc%carbon13%ptr(ics:ice,:,iblk) - & 
                                                loss%dom_assoc%carbon13%ptr(ics:ice,:,iblk) + & 
                                                transport%dom_assoc%carbon13%ptr(ics:ice,:,iblk)  
    pool%soluable_litter%carbon13%ptr(ics:ice,:,iblk) = pool%soluable_litter%carbon13%ptr(ics:ice,:,iblk) + & 
                                                formation%soluable_litter%carbon13%ptr(ics:ice,:,iblk) - & 
                                                loss%soluable_litter%carbon13%ptr(ics:ice,:,iblk) + & 
                                                transport%soluable_litter%carbon13%ptr(ics:ice,:,iblk)  
    pool%polymeric_litter%carbon13%ptr(ics:ice,:,iblk) = pool%polymeric_litter%carbon13%ptr(ics:ice,:,iblk) + & 
                                                formation%polymeric_litter%carbon13%ptr(ics:ice,:,iblk) - & 
                                                loss%polymeric_litter%carbon13%ptr(ics:ice,:,iblk) + & 
                                                transport%polymeric_litter%carbon13%ptr(ics:ice,:,iblk)  
    pool%woody_litter%carbon13%ptr(ics:ice,:,iblk) = pool%woody_litter%carbon13%ptr(ics:ice,:,iblk) + & 
                                                formation%woody_litter%carbon13%ptr(ics:ice,:,iblk) - & 
                                                loss%woody_litter%carbon13%ptr(ics:ice,:,iblk) + & 
                                                transport%woody_litter%carbon13%ptr(ics:ice,:,iblk)  
    pool%microbial%carbon13%ptr(ics:ice,:,iblk) = pool%microbial%carbon13%ptr(ics:ice,:,iblk) + & 
                                                formation%microbial%carbon13%ptr(ics:ice,:,iblk) - & 
                                                loss%microbial%carbon13%ptr(ics:ice,:,iblk) + & 
                                                transport%microbial%carbon13%ptr(ics:ice,:,iblk)  
    pool%residue%carbon13%ptr(ics:ice,:,iblk) = pool%residue%carbon13%ptr(ics:ice,:,iblk) + & 
                                                formation%residue%carbon13%ptr(ics:ice,:,iblk) - & 
                                                loss%residue%carbon13%ptr(ics:ice,:,iblk) + & 
                                                transport%residue%carbon13%ptr(ics:ice,:,iblk)  
    pool%residue_assoc%carbon13%ptr(ics:ice,:,iblk) = pool%residue_assoc%carbon13%ptr(ics:ice,:,iblk) + & 
                                                formation%residue_assoc%carbon13%ptr(ics:ice,:,iblk) - & 
                                                loss%residue_assoc%carbon13%ptr(ics:ice,:,iblk) + & 
                                                transport%residue_assoc%carbon13%ptr(ics:ice,:,iblk)  
    ! C14
    pool%dom%carbon14%ptr(ics:ice,:,iblk) = pool%dom%carbon14%ptr(ics:ice,:,iblk) + & 
                                                formation%dom%carbon14%ptr(ics:ice,:,iblk) - & 
                                                loss%dom%carbon14%ptr(ics:ice,:,iblk) + & 
                                                transport%dom%carbon14%ptr(ics:ice,:,iblk)  
    pool%dom_assoc%carbon14%ptr(ics:ice,:,iblk) = pool%dom_assoc%carbon14%ptr(ics:ice,:,iblk) + & 
                                                formation%dom_assoc%carbon14%ptr(ics:ice,:,iblk) - & 
                                                loss%dom_assoc%carbon14%ptr(ics:ice,:,iblk) + & 
                                                transport%dom_assoc%carbon14%ptr(ics:ice,:,iblk)  
    pool%soluable_litter%carbon14%ptr(ics:ice,:,iblk) = pool%soluable_litter%carbon14%ptr(ics:ice,:,iblk) + & 
                                                formation%soluable_litter%carbon14%ptr(ics:ice,:,iblk) - & 
                                                loss%soluable_litter%carbon14%ptr(ics:ice,:,iblk) + & 
                                                transport%soluable_litter%carbon14%ptr(ics:ice,:,iblk)  
    pool%polymeric_litter%carbon14%ptr(ics:ice,:,iblk) = pool%polymeric_litter%carbon14%ptr(ics:ice,:,iblk) + & 
                                                formation%polymeric_litter%carbon14%ptr(ics:ice,:,iblk) - & 
                                                loss%polymeric_litter%carbon14%ptr(ics:ice,:,iblk) + & 
                                                transport%polymeric_litter%carbon14%ptr(ics:ice,:,iblk)  
    pool%woody_litter%carbon14%ptr(ics:ice,:,iblk) = pool%woody_litter%carbon14%ptr(ics:ice,:,iblk) + & 
                                                formation%woody_litter%carbon14%ptr(ics:ice,:,iblk) - & 
                                                loss%woody_litter%carbon14%ptr(ics:ice,:,iblk) + & 
                                                transport%woody_litter%carbon14%ptr(ics:ice,:,iblk)  
    pool%microbial%carbon14%ptr(ics:ice,:,iblk) = pool%microbial%carbon14%ptr(ics:ice,:,iblk) + & 
                                                formation%microbial%carbon14%ptr(ics:ice,:,iblk) - & 
                                                loss%microbial%carbon14%ptr(ics:ice,:,iblk) + & 
                                                transport%microbial%carbon14%ptr(ics:ice,:,iblk)  
    pool%residue%carbon14%ptr(ics:ice,:,iblk) = pool%residue%carbon14%ptr(ics:ice,:,iblk) + & 
                                                formation%residue%carbon14%ptr(ics:ice,:,iblk) - & 
                                                loss%residue%carbon14%ptr(ics:ice,:,iblk) + & 
                                                transport%residue%carbon14%ptr(ics:ice,:,iblk)  
    pool%residue_assoc%carbon14%ptr(ics:ice,:,iblk) = pool%residue_assoc%carbon14%ptr(ics:ice,:,iblk) + & 
                                                formation%residue_assoc%carbon14%ptr(ics:ice,:,iblk) - & 
                                                loss%residue_assoc%carbon14%ptr(ics:ice,:,iblk) + & 
                                                transport%residue_assoc%carbon14%ptr(ics:ice,:,iblk)  
    ! N15
    pool%dom%nitrogen15%ptr(ics:ice,:,iblk) = pool%dom%nitrogen15%ptr(ics:ice,:,iblk) + & 
                                                formation%dom%nitrogen15%ptr(ics:ice,:,iblk) - & 
                                                loss%dom%nitrogen15%ptr(ics:ice,:,iblk) + & 
                                                transport%dom%nitrogen15%ptr(ics:ice,:,iblk)  
    pool%dom_assoc%nitrogen15%ptr(ics:ice,:,iblk) = pool%dom_assoc%nitrogen15%ptr(ics:ice,:,iblk) + & 
                                                formation%dom_assoc%nitrogen15%ptr(ics:ice,:,iblk) - & 
                                                loss%dom_assoc%nitrogen15%ptr(ics:ice,:,iblk) + & 
                                                transport%dom_assoc%nitrogen15%ptr(ics:ice,:,iblk)  
    pool%soluable_litter%nitrogen15%ptr(ics:ice,:,iblk) = pool%soluable_litter%nitrogen15%ptr(ics:ice,:,iblk) + & 
                                                formation%soluable_litter%nitrogen15%ptr(ics:ice,:,iblk) - & 
                                                loss%soluable_litter%nitrogen15%ptr(ics:ice,:,iblk) + & 
                                                transport%soluable_litter%nitrogen15%ptr(ics:ice,:,iblk)  
    pool%polymeric_litter%nitrogen15%ptr(ics:ice,:,iblk) = pool%polymeric_litter%nitrogen15%ptr(ics:ice,:,iblk) + & 
                                                formation%polymeric_litter%nitrogen15%ptr(ics:ice,:,iblk) - & 
                                                loss%polymeric_litter%nitrogen15%ptr(ics:ice,:,iblk) + & 
                                                transport%polymeric_litter%nitrogen15%ptr(ics:ice,:,iblk)  
    pool%woody_litter%nitrogen15%ptr(ics:ice,:,iblk) = pool%woody_litter%nitrogen15%ptr(ics:ice,:,iblk) + & 
                                                formation%woody_litter%nitrogen15%ptr(ics:ice,:,iblk) - & 
                                                loss%woody_litter%nitrogen15%ptr(ics:ice,:,iblk) + & 
                                                transport%woody_litter%nitrogen15%ptr(ics:ice,:,iblk)  
    pool%microbial%nitrogen15%ptr(ics:ice,:,iblk) = pool%microbial%nitrogen15%ptr(ics:ice,:,iblk) + & 
                                                formation%microbial%nitrogen15%ptr(ics:ice,:,iblk) - & 
                                                loss%microbial%nitrogen15%ptr(ics:ice,:,iblk) + & 
                                                transport%microbial%nitrogen15%ptr(ics:ice,:,iblk)  
    pool%residue%nitrogen15%ptr(ics:ice,:,iblk) = pool%residue%nitrogen15%ptr(ics:ice,:,iblk) + & 
                                                formation%residue%nitrogen15%ptr(ics:ice,:,iblk) - & 
                                                loss%residue%nitrogen15%ptr(ics:ice,:,iblk) + & 
                                                transport%residue%nitrogen15%ptr(ics:ice,:,iblk)  
    pool%residue_assoc%nitrogen15%ptr(ics:ice,:,iblk) = pool%residue_assoc%nitrogen15%ptr(ics:ice,:,iblk) + & 
                                                formation%residue_assoc%nitrogen15%ptr(ics:ice,:,iblk) - & 
                                                loss%residue_assoc%nitrogen15%ptr(ics:ice,:,iblk) + & 
                                                transport%residue_assoc%nitrogen15%ptr(ics:ice,:,iblk)  


  END SUBROUTINE calc_jsm_pool_change

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
  !! Input: soil temperature, activiation energy 
  !!
  !! Output: rate modifier (unitless) 
  !-----------------------------------------------------------------------------------------------------
  ELEMENTAL FUNCTION calc_temperature_rate_modifier(t_soil,Ea) RESULT(rate_modifier)

    USE mo_shared_constants,    ONLY: r_gas
    USE mo_sb_constants,        ONLY: t_ref_tresponse 

    IMPLICIT NONE

    !0.1 InOut
    REAL(wp), INTENT(in)  :: t_soil        !< soil temperature [K]
    REAL(wp), INTENT(in)  :: Ea            !< activation energy [J mol-1]
    REAL(wp)              :: rate_modifier !< temperature rate modifier [unitless] 
    ! 0.2 Local
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_temperature_rate_modifier'

    IF(Ea>eps8)THEN
       rate_modifier = exp(-Ea/r_gas*(1._wp/t_soil - 1._wp/t_ref_tresponse))
    ELSE
       rate_modifier = 1._wp
    ENDIF   

  END FUNCTION calc_temperature_rate_modifier

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

    !0.1 InOut
    REAL(wp), INTENT(in)  :: w_soil_sl     !< soil moisture [m]
    REAL(wp), INTENT(in)  :: w_soil_fc_sl  !< soil moisture at field capacity [m]
    REAL(wp), INTENT(in)  :: param         !< 
    REAL(wp)              :: rate_modifier !< oxygen availability rate modifier [unitless] 
    ! 0.2 Local
    REAL(wp)              :: air_filled_pore_space
    REAL(wp)              :: hlp
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_oxygen_rate_modifier'

    IF(param>eps8) THEN
       IF(w_soil_fc_sl > eps8) THEN
          hlp = w_soil_fc_sl - w_soil_sl
          IF (hlp < 0.0_wp) THEN
            hlp = 0.0_wp
          END IF
          air_filled_pore_space = hlp / w_soil_fc_sl
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
  ELEMENTAL FUNCTION calc_convert_bulk_to_liquid(w_soil_sl,param) RESULT(factor)

    IMPLICIT NONE

    !0.1 InOut
    REAL(wp), INTENT(in)  :: w_soil_sl           !< soil moisture [m]
    REAL(wp), INTENT(in)  :: param               !< 
    REAL(wp)              :: factor              !< conversion factor [unitless] 
    ! 0.2 Local
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_convert_bulk_to_liquid'

    IF(param>eps8) THEN
       IF(w_soil_sl > eps4) THEN
         factor = 1._wp / w_soil_sl 
       ELSE
          factor = 1._wp / eps4 
       ENDIF
    ELSE 
       factor = 1._wp
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

    !0.1 InOut
    REAL(wp), INTENT(in)  :: w_soil_sl           !< soil moisture [m]
    REAL(wp), INTENT(in)  :: w_soil_fc_sl        !< soil moisture at field capacity [m]
    REAL(wp), INTENT(in)  :: param               !< 
    REAL(wp)              :: rate_modifier !< oxygen availability rate modifier [unitless] 
    ! 0.2 Local
    REAL(wp)              :: air_filled_pore_space
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_substrate_rate_modifier'
    
    IF(param>eps8.AND.w_soil_fc_sl > eps8) THEN
       air_filled_pore_space = (w_soil_sl - w_soil_fc_sl ) / w_soil_fc_sl
       IF(air_filled_pore_space > eps4) THEN
          rate_modifier = 1._wp / air_filled_pore_space
       ELSE 
          rate_modifier = 1._wp / eps4
       ENDIF
    ELSE
       rate_modifier = 1._wp
    ENDIF

  END FUNCTION calc_substrate_rate_modifier

END MODULE mo_sb_jsm_processes
