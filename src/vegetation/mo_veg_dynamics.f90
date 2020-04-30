!>
!!  @par Copyright
!!    QUINCY terrestrial biosphere model with coupled carbon, nitrogen, and phosphorus cycles
!!    Copyright (C) 2019 Max Planck Institute for Biogeochemistry (MPI-BGC), Jena, Germany 
!!    Contributors: Caldararu S., Engel J., Kern M., Pallandt M., Schnur R., Thum T., Yu L., Zaehle S.
!!
!!    This file is part of the QUINCY model.
!!    The QUINCY model is free software: you can redistribute it and/or modify
!!    it under the terms of the GNU General Public License version 3,
!!    as published by the Free Software Foundation.
!!
!!    The QUINCY model is distributed in the hope that it will be useful,
!!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!!    GNU General Public License for more details.
!!
!!    You should have received a copy of the GNU General Public License
!!    along with the QUINCY model.  If not, see <https://www.gnu.org/licenses/>.
!! 
!! @brief module to calculate vegetation dynamics 
!!
!! @author
!!  Sönke Zaehle, MPI-BGC Jena
!!
!!
MODULE mo_veg_dynamics

  USE mo_kind,                 ONLY: wp
  ! USE mo_jsb_impl_constants,   ONLY: true, false, test_false_true
  USE mo_jsb_control,          ONLY: debug_on
  USE mo_exception,            ONLY: message

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: update_veg_dynamics

  CHARACTER(len=*), PARAMETER :: modname = 'mo_veg_dynamics'

CONTAINS

  !-----------------------------------------------------------------------------------------------------
  ! Main Task
  ! 
  !-----------------------------------------------------------------------------------------------------
  !> Calculate establishment and mortality of a population (currently static mortality rate only) 
  !!
  !! a background_mort_rate exists for each trees and grasses, defined as mortality rate per year
  !!
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE update_veg_dynamics(tile, options)

    USE mo_jsb_lctlib_class,      ONLY: t_lctlib_element
    USE mo_jsb_class,             ONLY: Get_model
    USE mo_jsb_tile_class,        ONLY: t_jsb_tile_abstract
    USE mo_jsb_task_class,        ONLY: t_jsb_task_options
    ! USE mo_jsb_math_constants,    ONLY: one_year, one_day, eps4
    USE mo_jsb_model_class,       ONLY: t_jsb_model
    ! USE mo_jsb_process_class,     ONLY: A2L_, L2A_, SEB_, TURB_, SSE_, HYDRO_, VEG_, HD_, RAD_, ASSIMI_, PHENO_, SB_
    USE mo_jsb_process_class,     ONLY: ASSIMI_, VEG_
    ! USE mo_jsb4_forcing,          ONLY: forcing_options
    ! USE mo_veg_constants,         ONLY: itree

    ! Use of pools
    dsl4jsb_Use_pools(VEG_) t_veg_main_pool, t_veg_tissues_pool, t_veg_elements_pool
    dsl4jsb_Use_elements(VEG_) carbon_, nitrogen_, phosphorus_, carbon13_, carbon14_, nitrogen15_

    ! Use of process configurations (t_PROC_config)
    dsl4jsb_Use_config(ASSIMI_)
    dsl4jsb_Use_config(VEG_)

    ! Use of process memories
    dsl4jsb_Use_memory(VEG_)
    ! dsl4slm_Use_memory_a_veg
    ! USE mo_veg_memory_class,      ONLY: veg_pool, veg_litterfall, veg_establishment
    ! dsl4jsb_Use_memory(RAD_)
    ! dsl4jsb_Use_memory(PHENO_)
    ! dsl4jsb_Use_memory(SSE_)

    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    CLASS(t_jsb_tile_abstract), INTENT(inout)     :: tile         !< one tile with data structure for one lct
    TYPE(t_jsb_task_options),   INTENT(in)        :: options      !< model options
    ! ---------------------------
    ! 0.2 Local 
    TYPE(t_jsb_model),      POINTER         :: model
    TYPE(t_lctlib_element), POINTER         :: lctlib           !< land-cover-type library - parameter across pft's
    ! REAL(wp), DIMENSION(options%nc)         :: fpc              !< current foliage projective cover
    INTEGER                                 :: iblk, ics, ice, nc
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':update_veg_dynamics'
    ! ---------------------------
    ! 0.3 Declare Memory
    ! Declare process configuration and memory Pointers
    dsl4jsb_Def_config(ASSIMI_)
    dsl4jsb_Def_config(VEG_)
    dsl4jsb_Def_memory(VEG_)                              ! TYPE(t_veg_memory),      POINTER :: veg__mem
    ! dsl4jsb_Def_memory(RAD_)                              ! TYPE(t_rad_memory),      POINTER :: rad__mem
    ! dsl4jsb_Def_memory(PHENO_)                            ! TYPE(t_pheno_memory),    POINTER :: pheno__mem
    ! dsl4jsb_Def_memory(SSE_)                             ! TYPE(t_sse_memory),     POINTER :: sse__mem
    dsl4jsb_Def_pool(t_veg_tissues_pool) :: veg_pool
    dsl4jsb_Def_pool(t_veg_elements_pool) :: veg_pool_leaf
    ! ! Declare pointers to variables in memory
    ! ! RAD_
    ! dsl4jsb_Real2D_onChunk      :: rfr_ratio_boc_tvegdyn_mavg
    ! ! PHENO_
    ! dsl4jsb_Real2D_onChunk      :: growing_season
    ! ! SSE_
    ! dsl4jsb_Real2D_onChunk      :: w_soil_root_theta
    ! ! VEG_
    ! dsl4jsb_Real2D_onChunk      :: do_cohort_harvest
    ! dsl4jsb_Real2D_onChunk      :: dens_ind
    ! dsl4jsb_Real2D_onChunk      :: diameter 
    ! dsl4jsb_Real2D_onChunk      :: height 
    ! dsl4jsb_Real2D_onChunk      :: cohort_age
    ! dsl4jsb_Real2D_onChunk      :: mortality_rate
    ! dsl4jsb_Real2D_onChunk      :: delta_dens_ind
    ! dsl4jsb_Real2D_onChunk      :: t_air_week_mavg
    ! dsl4jsb_Real2D_onChunk      :: an_boc_tvegdyn_mavg
    ! dsl4jsb_Real2D_onChunk      :: net_growth_tvegdyn_mavg
    ! dsl4jsb_Real2D_onChunk      :: lai_tvegdyn_mavg
    ! Get local variables from options argument
    iblk    = options%iblk
    ics     = options%ics
    ice     = options%ice
    nc      = options%nc
    ! ---------------------------
    ! 0.4 Process Activity, Debug Option
    IF (.NOT. tile%Is_process_active(VEG_)) RETURN
    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')
    ! ---------------------------
    ! 0.5 Get Memory
    model  => Get_model(tile%owner_model_id)
    lctlib => model%lctlib(tile%lcts(1)%lib_id)
    ! Get process config
    dsl4jsb_Get_config(ASSIMI_)
    dsl4jsb_Get_config(VEG_)
    ! Get process memories
    dsl4jsb_Get_memory(VEG_)
!     dsl4jsb_Get_memory(RAD_)
!     dsl4jsb_Get_memory(PHENO_)
!     dsl4jsb_Get_memory(SSE_)
    dsl4jsb_Get_pool(VEG_, veg_main_pool%veg_pool,      veg_pool)
    dsl4jsb_Get_pool(VEG_, veg_main_pool%veg_pool%leaf, veg_pool_leaf)
!     ! Get process variables (Set pointers to variables in memory)
!     dsl4jsb_Get_var2D_onChunk(RAD_, rfr_ratio_boc_tvegdyn_mavg) ! in
!     ! ---------------------------
!     dsl4jsb_Get_var2D_onChunk(PHENO_, growing_season)       ! in
!     ! ---------------------------
!     dsl4jsb_Get_var2D_onChunk(SSE_, w_soil_root_theta)      ! in
!     ! ---------------------------
!     dsl4jsb_Get_var2D_onChunk(VEG_, do_cohort_harvest)      ! out
!     dsl4jsb_Get_var2D_onChunk(VEG_, dens_ind)               ! in
!     dsl4jsb_Get_var2D_onChunk(VEG_, diameter)               ! in
!     dsl4jsb_Get_var2D_onChunk(VEG_, height)                 ! in
!     dsl4jsb_Get_var2D_onChunk(VEG_, cohort_age)             ! inout
!     dsl4jsb_Get_var2D_onChunk(VEG_, mortality_rate)         ! out
!     dsl4jsb_Get_var2D_onChunk(VEG_, delta_dens_ind)         ! out
!     dsl4jsb_Get_var2D_onChunk(VEG_, t_air_week_mavg)        ! in
!     dsl4jsb_Get_var2D_onChunk(VEG_, an_boc_tvegdyn_mavg)    ! in
!     dsl4jsb_Get_var2D_onChunk(VEG_, net_growth_tvegdyn_mavg)! in
!     dsl4jsb_Get_var2D_onChunk(VEG_, lai_tvegdyn_mavg)       ! in
!     ! Get process memory pool structures
! #define Get_veg_pool
! #define Get_veg_litterfall
! #define Get_veg_establishment
! #include "slm_pools_a_incl.f90"
! #undef Get_veg_pool
! #undef Get_veg_litterfall
! #undef Get_veg_establishment



    ! ------------------------------------------------------------------------------------------------------------
    ! Go science
    ! ------------------------------------------------------------------------------------------------------------

    !! tested for VEG_ pools & fluxes !!
    veg_pool%leaf%carbon = 10.0_wp
    veg_pool_leaf%carbon = 20.0_wp

  END SUBROUTINE update_veg_dynamics


END MODULE mo_veg_dynamics
