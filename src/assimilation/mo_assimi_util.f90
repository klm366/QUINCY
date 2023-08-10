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
!! @brief diverse subroutines and functions for the assimilation process
!!
!! @author
!!  Jan Engel, MPI-BGC Jena
!!
!! @par Revision History
!! First version by Jan Engel (2018-04-06)
!!
MODULE mo_assimi_util
#ifndef __NO_JSBACH__

  USE mo_kind,        ONLY: wp
  USE mo_jsb_control, ONLY: debug_on
  USE mo_exception,   ONLY: message

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: reset_assimi_fluxes, calculate_time_average_assimilation

  CHARACTER(len=*), PARAMETER :: modname = 'mo_assimi_util'


CONTAINS


  ! ====================================================================================================== !
  !
  !> init/reset assimilation fluxes (with/to zero) - called prior to other assimi tasks
  !
  SUBROUTINE reset_assimi_fluxes(tile, options)

    USE mo_jsb_class,             ONLY: Get_model
    USE mo_jsb_tile_class,        ONLY: t_jsb_tile_abstract
    USE mo_jsb_task_class,        ONLY: t_jsb_task_options
    USE mo_jsb_model_class,       ONLY: t_jsb_model
    USE mo_jsb_math_constants,    ONLY: zero

    dsl4jsb_Use_processes ASSIMI_
    dsl4jsb_Use_memory(ASSIMI_)
    ! -------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_tile_abstract), INTENT(inout)     :: tile         !< one tile with data structure for one lct
    TYPE(t_jsb_task_options),   INTENT(in)        :: options      !< model options
    ! -------------------------------------------------------------------------------------------------- !
    TYPE(t_jsb_model),      POINTER       :: model
    INTEGER                               :: iblk, ics, ice, nc
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':reset_assimi_fluxes'

    dsl4jsb_Def_memory(ASSIMI_)
    ! -------------------------------------------------------------------------------------------------- !
    ! Get local variables from options argument
    iblk    = options%iblk
    ics     = options%ics
    ice     = options%ice
    nc      = options%nc

    IF (.NOT. tile%Is_process_calculated(ASSIMI_)) RETURN

    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    ! Get process memories
    dsl4jsb_Get_memory(ASSIMI_)

    ! Get process variables (Set pointers to variables in memory)
    model  => Get_model(tile%owner_model_id)

    ! ------------------------------------------------------------------------------------------------------------
    ! Go science
    ! ------------------------------------------------------------------------------------------------------------


    ! set zero all flux variables of the assimilation process
    dsl4jsb_var2D_onChunk(ASSIMI_, gross_assimilation)           = zero
    dsl4jsb_var2D_onChunk(ASSIMI_, gross_assimilation_C13)       = zero
    dsl4jsb_var2D_onChunk(ASSIMI_, gross_assimilation_C14)       = zero
    dsl4jsb_var2D_onChunk(ASSIMI_, net_assimilation)             = zero
    dsl4jsb_var2D_onChunk(ASSIMI_, net_assimilation_boc)         = zero
    dsl4jsb_var2D_onChunk(ASSIMI_, maint_respiration_leaf)       = zero

    dsl4jsb_var3D_onChunk(ASSIMI_, gross_assimilation_cl)      = zero
    dsl4jsb_var3D_onChunk(ASSIMI_, net_assimilation_cl)        = zero
    dsl4jsb_var3D_onChunk(ASSIMI_, maint_respiration_leaf_cl)  = zero

    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE reset_assimi_fluxes


  ! ======================================================================================================= !
  !> calculate time moving averages and daytime averages for ASSIMI_
  !>
  SUBROUTINE calculate_time_average_assimilation(tile, options)

    USE mo_kind,                    ONLY: wp
    USE mo_jsb_math_constants,      ONLY: eps8
    USE mo_jsb_physical_constants,  ONLY: Tzero
    USE mo_jsb_impl_constants,      ONLY: test_false_true
    USE mo_jsb_control,             ONLY: debug_on
    USE mo_jsb_tile_class,          ONLY: t_jsb_tile_abstract
    USE mo_jsb_task_class,          ONLY: t_jsb_task_options
    USE mo_time_averages            ! e.g. calc_time_mavg, mavg_period_tphen, mavg_period_weekly
    dsl4jsb_Use_processes A2L_, ASSIMI_, PHENO_
    dsl4jsb_Use_memory(A2L_)
    dsl4jsb_Use_memory(ASSIMI_)
    dsl4jsb_Use_memory(PHENO_)
    ! ----------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options
    ! ----------------------------------------------------------------------------------------------------- !
    LOGICAL,  DIMENSION(options%nc)         :: l_growing_season     !< growing_season LOGICAL
    REAl(wp), DIMENSION(options%nc)         :: hlp                  !< helper
    INTEGER                                 :: ic                   !< looping
    INTEGER                                 :: iblk, ics, ice, nc   !< dimensions
    CHARACTER(len=*), PARAMETER             :: routine = modname//':calculate_time_average_assimilation'
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Def_memory(A2L_)
    dsl4jsb_Def_memory(ASSIMI_)
    dsl4jsb_Def_memory(PHENO_)
    ! ----------------------------------------------------------------------------------------------------- !
    ! A2L_
    dsl4jsb_Real2D_onChunk      :: t_air
    dsl4jsb_Real2D_onChunk      :: swpar_srf_down
    dsl4jsb_Real2D_onChunk      :: prev_one_cos_zenith_angle
    dsl4jsb_Real2D_onChunk      :: prev_two_cos_zenith_angle
    dsl4jsb_Real2D_onChunk      :: prev_three_cos_zenith_angle
    dsl4jsb_Real2D_onChunk      :: daytime_counter
    ! ASSIMI_
    dsl4jsb_Real2D_onChunk      :: beta_air
    dsl4jsb_Real2D_onChunk      :: beta_soil_ps
    dsl4jsb_Real2D_onChunk      :: beta_soil_gs
    dsl4jsb_Real2D_onChunk      :: beta_soa
    dsl4jsb_Real2D_onChunk      :: soa_tsoa_mavg
    dsl4jsb_Real2D_onChunk      :: beta_soa_tphen_mavg
    dsl4jsb_Real2D_onChunk      :: beta_soil_gs_tphen_mavg
    dsl4jsb_Real2D_onChunk      :: beta_soil_ps_tfrac_mavg
    dsl4jsb_Real2D_onChunk      :: beta_soil_gs_tfrac_mavg
    dsl4jsb_Real2D_onChunk      :: beta_soil_ps_tcnl_mavg
    dsl4jsb_Real2D_onChunk      :: beta_soil_gs_tcnl_mavg
    dsl4jsb_Real2D_onChunk      :: beta_air_tfrac_mavg
    dsl4jsb_Real2D_onChunk      :: beta_air_tcnl_mavg
    dsl4jsb_Real2D_onChunk      :: beta_air_daytime
    dsl4jsb_Real2D_onChunk      :: beta_soil_ps_daytime
    dsl4jsb_Real2D_onChunk      :: beta_soil_gs_daytime
    dsl4jsb_Real2D_onChunk      :: beta_air_daytime_dacc
    dsl4jsb_Real2D_onChunk      :: beta_soil_ps_daytime_dacc
    dsl4jsb_Real2D_onChunk      :: beta_soil_gs_daytime_dacc
    ! PHENO_
    dsl4jsb_Real2D_onChunk      :: growing_season
    ! ----------------------------------------------------------------------------------------------------- !
    iblk    = options%iblk
    ics     = options%ics
    ice     = options%ice
    nc      = options%nc
    ! ----------------------------------------------------------------------------------------------------- !
    IF (.NOT. tile%Is_process_calculated(ASSIMI_)) RETURN
    IF (tile%lcts(1)%lib_id == 0)              RETURN  ! only if the present tile is a pft
    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Get_memory(A2L_)
    dsl4jsb_Get_memory(ASSIMI_)
    dsl4jsb_Get_memory(PHENO_)
    ! ----------------------------------------------------------------------------------------------------- !
    ! A2L_
    dsl4jsb_Get_var2D_onChunk(A2L_, t_air)                            ! in
    dsl4jsb_Get_var2D_onChunk(A2L_, swpar_srf_down)                   ! in
    dsl4jsb_Get_var2D_onChunk(A2L_, prev_one_cos_zenith_angle)        ! in
    dsl4jsb_Get_var2D_onChunk(A2L_, prev_three_cos_zenith_angle)      ! in
    dsl4jsb_Get_var2D_onChunk(A2L_, prev_two_cos_zenith_angle)        ! in
    dsl4jsb_Get_var2D_onChunk(A2L_, daytime_counter)                  ! in
    ! ASSIMI_
    dsl4jsb_Get_var2D_onChunk(ASSIMI_, beta_air)                      ! in
    dsl4jsb_Get_var2D_onChunk(ASSIMI_, beta_soil_ps)                  ! in
    dsl4jsb_Get_var2D_onChunk(ASSIMI_, beta_soil_gs)                  ! in
    dsl4jsb_Get_var2D_onChunk(ASSIMI_, beta_soa)                      ! in
    dsl4jsb_Get_var2D_onChunk(ASSIMI_, soa_tsoa_mavg)                 ! out
    dsl4jsb_Get_var2D_onChunk(ASSIMI_, beta_soa_tphen_mavg)           ! out
    dsl4jsb_Get_var2D_onChunk(ASSIMI_, beta_soil_gs_tphen_mavg)       ! out
    dsl4jsb_Get_var2D_onChunk(ASSIMI_, beta_soil_ps_tfrac_mavg)       ! out
    dsl4jsb_Get_var2D_onChunk(ASSIMI_, beta_soil_gs_tfrac_mavg)       ! out
    dsl4jsb_Get_var2D_onChunk(ASSIMI_, beta_soil_ps_tcnl_mavg)        ! out
    dsl4jsb_Get_var2D_onChunk(ASSIMI_, beta_soil_gs_tcnl_mavg)        ! out
    dsl4jsb_Get_var2D_onChunk(ASSIMI_, beta_air_tfrac_mavg)           ! out
    dsl4jsb_Get_var2D_onChunk(ASSIMI_, beta_air_tcnl_mavg)            ! out
    dsl4jsb_Get_var2D_onChunk(ASSIMI_, beta_air_daytime)              ! out
    dsl4jsb_Get_var2D_onChunk(ASSIMI_, beta_soil_ps_daytime)          ! out
    dsl4jsb_Get_var2D_onChunk(ASSIMI_, beta_soil_gs_daytime)          ! out
    dsl4jsb_Get_var2D_onChunk(ASSIMI_, beta_air_daytime_dacc)         ! out
    dsl4jsb_Get_var2D_onChunk(ASSIMI_, beta_soil_ps_daytime_dacc)     ! out
    dsl4jsb_Get_var2D_onChunk(ASSIMI_, beta_soil_gs_daytime_dacc)     ! out
    ! PHENO_
    dsl4jsb_Get_var2D_onChunk(PHENO_, growing_season)                 ! in


    !>0.9 transform REAL growing_season(:) into LOGICAL l_growing_season(:)/_cl(:,:)
    !>
    DO ic = 1, nc
      IF (growing_season(ic) > test_false_true) THEN
        l_growing_season(ic) = .TRUE.
      ELSE
        l_growing_season(ic) = .FALSE.
      END IF
    END DO

    ! ----------------------------------------------------------------------------------------------------- !
    !>1.0 daytime averages
    !>
    !>
    
    !>  1.1 accumulate values - if light is available (i.e. daytime)
    !>
    DO ic = 1, nc
      IF (swpar_srf_down(ic) > eps8) THEN
        ! add values
        beta_air_daytime_dacc(ic)       = beta_air_daytime_dacc(ic)     + beta_air(ic)
        beta_soil_ps_daytime_dacc(ic)   = beta_soil_ps_daytime_dacc(ic) + beta_soil_ps(ic)
        beta_soil_gs_daytime_dacc(ic)   = beta_soil_gs_daytime_dacc(ic) + beta_soil_gs(ic)
      END IF
    END DO

    !>  1.2 calc daytime averages - at the timestep after local midnight (1st timestep of the new day)
    !>
    !>    calculate the average of the previous day and pass this value to the according variable
    !>
    !>    the timestep right after midnight is identified by the changes in cos_zenith_angle (turning from negative to positive)
    !>    these are traced (and shifted) in update_a2l/update_atm2land() or directly in the ICON-Land standalone driver
    !>    Thus the current cos_zenith_angle is actually stored in prev_one_cos_zenith_angle
    DO ic = 1, nc
      IF ((prev_one_cos_zenith_angle(ic) - prev_two_cos_zenith_angle(ic)) > eps8 .AND. &
          & (prev_two_cos_zenith_angle(ic) - prev_three_cos_zenith_angle(ic)) <= eps8) THEN
        ! check if at least one value had been assigned at the current day, avoiding division by zero
        IF (daytime_counter(ic) > eps8) THEN
          beta_air_daytime(ic)     = beta_air_daytime_dacc(ic)     / daytime_counter(ic)
          beta_soil_ps_daytime(ic) = beta_soil_ps_daytime_dacc(ic) / daytime_counter(ic)
          beta_soil_gs_daytime(ic) = beta_soil_gs_daytime_dacc(ic) / daytime_counter(ic)
        END IF
        ! zero the accumulation variables after daily average has been calculated
        beta_air_daytime_dacc(ic)     = 0.0_wp
        beta_soil_ps_daytime_dacc(ic) = 0.0_wp
        beta_soil_gs_daytime_dacc(ic) = 0.0_wp
      END IF
    END DO


    ! ----------------------------------------------------------------------------------------------------- !
    !>2.0 moving averages
    !>
    !>

    ! docu:
    ! calc_time_mavg(current average, new value, length of avg_period,  !
    !                do_calc=LOGICAL, avg_period_unit='day')            ! OPTIONAL
    !                RETURN(new current average)
    ! the unit of the averaging period is 'day' by default, but can also be 'week' or 'year'
    
    !>  2.1 tphen (averages at the time-scale of phenology)
    !>
    beta_soil_gs_tphen_mavg(:) = calc_time_mavg(beta_soil_gs_tphen_mavg(:), beta_soil_gs(:), mavg_period_tphen)
    beta_soa_tphen_mavg(:)     = calc_time_mavg(beta_soa_tphen_mavg(:), beta_soa(:), mavg_period_tphen)

    !>  2.2 soa time-scale
    !>
    hlp(:)                     = t_air(:) - Tzero     ! convert from Kelvin to deg Celsius
    soa_tsoa_mavg(:)           = calc_time_mavg(soa_tsoa_mavg(:), hlp(:), mavg_period_tsoa)

    !>  2.3 tfrac (averages at the time-scale of within-leaf N allocation fractions)
    !>
    beta_soil_ps_tfrac_mavg(:) = calc_time_mavg(beta_soil_ps_tfrac_mavg(:), beta_soil_ps_daytime(:), &
                                                mavg_period_tfrac)
    beta_soil_gs_tfrac_mavg(:) = calc_time_mavg(beta_soil_gs_tfrac_mavg(:), beta_soil_gs_daytime(:), &
                                                mavg_period_tfrac)
    beta_air_tfrac_mavg(:)     = calc_time_mavg(beta_air_tfrac_mavg(:), beta_air_daytime(:), &
                                                mavg_period_tfrac)

    !>  2.4 tcnl (averages at the time-scale of leaf N allocation fractions)
    !>
    beta_soil_ps_tcnl_mavg(:)  = calc_time_mavg(beta_soil_ps_tcnl_mavg(:), beta_soil_ps(:), &
                                                mavg_period_tcnl, do_calc=l_growing_season(:))
    beta_soil_gs_tcnl_mavg(:)  = calc_time_mavg(beta_soil_gs_tcnl_mavg(:), beta_soil_gs(:), &
                                                mavg_period_tcnl, do_calc=l_growing_season(:))
    beta_air_tcnl_mavg(:)      = calc_time_mavg(beta_air_tcnl_mavg(:), beta_air(:), &
                                                mavg_period_tcnl, do_calc=l_growing_season(:))

  END SUBROUTINE calculate_time_average_assimilation

#endif
END MODULE mo_assimi_util
