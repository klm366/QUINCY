!> QUINCY simple harvest variables init
!>
!> ICON-Land
!>
!> ---------------------------------------
!> Copyright (C) 2013-2024, MPI-M, MPI-BGC
!>
!> Contact: icon-model.org
!> Authors: AUTHORS.md
!> See LICENSES/ for license information
!> SPDX-License-Identifier: BSD-3-Clause
!> ---------------------------------------
!>
!> For more information on the QUINCY model see: <https://doi.org/10.17871/quincy-model-2019>
!>
!>#### initialization of simple harvest memory variables using, e.g., ic & bc input files
!>
MODULE mo_q_sh_init
#ifndef __NO_QUINCY__

  USE mo_kind,                ONLY: wp
  USE mo_exception,           ONLY: message
  USE mo_jsb_control,         ONLY: debug_on

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: q_sh_init

  CHARACTER(len=*), PARAMETER :: modname = 'mo_q_sh_init'

CONTAINS

  ! ======================================================================================================= !
  !> Run simple harvest init
  !>
  SUBROUTINE q_sh_init(tile)
    USE mo_jsb_class,             ONLY: Get_model
    USE mo_jsb_tile_class,        ONLY: t_jsb_tile_abstract
    USE mo_jsb_model_class,       ONLY: t_jsb_model
    USE mo_jsb_tile_class,        ONLY: t_jsb_tile_abstract
    USE mo_jsb_process_class,     ONLY: Q_SH_, VEG_
    dsl4jsb_Use_memory(Q_SH_)
    dsl4jsb_Use_config(VEG_)
    ! ----------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    ! ----------------------------------------------------------------------------------------------------- !
    TYPE(t_jsb_model), POINTER :: model
    CHARACTER(len=*), PARAMETER :: routine = modname//':q_sh_init'
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Def_memory(Q_SH_)
    dsl4jsb_Def_config(VEG_)
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Real2D_onDomain :: fract_wood_to_slash
    dsl4jsb_Real2D_onDomain :: fract_wood_to_pp_fuel
    dsl4jsb_Real2D_onDomain :: fract_wood_to_pp_paper
    dsl4jsb_Real2D_onDomain :: fract_wood_to_pp_fiberboard
    dsl4jsb_Real2D_onDomain :: fract_wood_to_pp_oirw
    dsl4jsb_Real2D_onDomain :: fract_wood_to_pp_pv
    dsl4jsb_Real2D_onDomain :: fract_wood_to_pp_sawnwood
    ! ----------------------------------------------------------------------------------------------------- !
    model => Get_model(tile%owner_model_id)
    IF (.NOT. tile%Is_process_calculated(Q_SH_)) RETURN
    ! ----------------------------------------------------------------------------------------------------- !
    IF (debug_on()) CALL message(TRIM(routine), 'Setting initial conditions of sh memory (quincy) for tile '// &
      &                          TRIM(tile%name))
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Get_memory(Q_SH_)
    dsl4jsb_Get_config(VEG_)
    ! ----------------------------------------------------------------------------------------------------- !
    IF (dsl4jsb_Config(VEG_)%l_use_product_pools) THEN
      dsl4jsb_Get_var2D_onDomain(Q_SH_, fract_wood_to_slash)
      dsl4jsb_Get_var2D_onDomain(Q_SH_, fract_wood_to_pp_fuel)
      dsl4jsb_Get_var2D_onDomain(Q_SH_, fract_wood_to_pp_paper)
      dsl4jsb_Get_var2D_onDomain(Q_SH_, fract_wood_to_pp_fiberboard)
      dsl4jsb_Get_var2D_onDomain(Q_SH_, fract_wood_to_pp_oirw)
      dsl4jsb_Get_var2D_onDomain(Q_SH_, fract_wood_to_pp_pv)
      dsl4jsb_Get_var2D_onDomain(Q_SH_, fract_wood_to_pp_sawnwood)
    END IF
    ! ----------------------------------------------------------------------------------------------------- !

    ! Harvested wood allocation fractions vary strongly regionally and temporally and might be read from input maps
    IF (dsl4jsb_Config(VEG_)%l_use_product_pools) THEN
      ! QS used a wood_extraction_rate of 0.8, thus, we use an initial default slash of 0.2 (= 1.0 - 0.8) if no maps are provided
      fract_wood_to_slash(:,:) = 0.2_wp

      ! Note: allocation fractions to product pools need to add up to one (asserted in update_pools_on_harvest)
      ! Below some European longterm averages inspired by numbers from Nuetzel (2021, LMU master thesis)
      fract_wood_to_pp_fuel(:,:)       = 0.4_wp
      fract_wood_to_pp_paper(:,:)      = 0.25_wp
      fract_wood_to_pp_fiberboard(:,:) = 0.05_wp
      fract_wood_to_pp_oirw(:,:)       = 0.05_wp
      fract_wood_to_pp_pv(:,:)         = 0.05_wp
      fract_wood_to_pp_sawnwood(:,:)   = 0.2_wp
    END IF

  END SUBROUTINE q_sh_init

#endif
END MODULE mo_q_sh_init
