!>
!> Contains the interfaces to the megan scheme
!>
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
!>
!> author
!>  Stiig Wilkenskjeld, MPI-M Hamburg
!>
!> Revision History
!> Adopted from JSBACH3 code by Thomas Kleinen, 2024-01-24
!>

MODULE mo_megan_interface
#ifndef __NO_JSBACH__

  USE mo_kind,                   ONLY: lp
  USE mo_exception,              ONLY: message, finish, message_text
  USE mo_jsb_control,            ONLY: debug_on
  USE mo_jsb_model_class,        ONLY: t_jsb_model
  USE mo_jsb_class,              ONLY: Get_model
  USE mo_jsb_grid_class,         ONLY: t_jsb_grid, t_jsb_vgrid
  USE mo_jsb_grid,               ONLY: Get_vgrid, Get_grid
  USE mo_jsb_tile_class,         ONLY: t_jsb_tile_abstract, t_jsb_aggregator
  USE mo_jsb_process_class,      ONLY: t_jsb_process, TOP_, A2L_, RAD_, SSE_, HYDRO_, PHENO_, MEGAN_
  USE mo_jsb_task_class,         ONLY: t_jsb_process_task, t_jsb_task_options
  USE mo_jsb_math_constants,     ONLY: eps8
  USE mo_jsb_physical_constants, ONLY: tmelt
  USE mo_megan_process,          ONLY: calc_gamma_lai, calc_gamma_sm, calc_gamma_temp_iso, calc_gamma_temp_other, &
    &                                  calc_gamma_pfd, calc_gamma_co2, calc_gamma_age, update_megan_runmean

  dsl4jsb_Use_processes MEGAN_
  dsl4jsb_Use_config(MEGAN_)
  dsl4jsb_Use_memory(MEGAN_)
  dsl4jsb_Use_memory(A2L_)
  dsl4jsb_Use_memory(HYDRO_)
  dsl4jsb_Use_memory(PHENO_)

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: Register_megan_tasks

  !> Type definition for topmodel task
  TYPE, EXTENDS(t_jsb_process_task) :: tsk_megan
  CONTAINS
    PROCEDURE, NOPASS :: Integrate => update_megan
    PROCEDURE, NOPASS :: Aggregate => aggregate_megan
  END TYPE tsk_megan

  !> Constructor interface for topmodel task
  INTERFACE tsk_megan
    PROCEDURE Create_task_megan
  END INTERFACE tsk_megan

  CHARACTER(len=*), PARAMETER :: modname = 'mo_megan_interface'

CONTAINS

  ! ================================================================================================================================
  !> Constructors for tasks

  ! -------------------------------------------------------------------------------------------------------
  !> Constructor for megan
  !>
  !> param[in]     model_id     Model id
  !> return        return_ptr   Instance of process task "topmodel"
  !>
  FUNCTION Create_task_megan(model_id) RESULT(return_ptr)

    INTEGER,                   INTENT(in) :: model_id
    CLASS(t_jsb_process_task), POINTER    :: return_ptr

    ALLOCATE(tsk_megan::return_ptr)
    CALL return_ptr%Construct(name='megan_processes', process_id=MEGAN_, owner_model_id=model_id)

  END FUNCTION Create_task_megan

  !-----------------------------------------------------------------------------------------------------
  !> Register tasks: MEGAN_ process
  !>
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE Register_megan_tasks(process, model_id)

    USE mo_jsb_model_class,   ONLY: t_jsb_model
    USE mo_jsb_class,         ONLY: Get_model

    CLASS(t_jsb_process), INTENT(inout) :: process
    INTEGER,              INTENT(in)    :: model_id

!    TYPE(t_jsb_model), POINTER :: model

!    model   => Get_model(model_id)

    CALL process%Register_task(tsk_megan(model_id))

  END SUBROUTINE Register_megan_tasks

  !-----------------------------------------------------------------------------------------------------
  !> update_megan
  !>
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE update_megan(tile, options)

    USE mo_jsb_time,           ONLY: is_time_experiment_start, is_time_restart, is_newyear, is_newmonth, is_newday, &
      &                              timesteps_per_day, get_date_components, get_year_length, get_month_length
    USE mo_megan_memory_class, ONLY: t_megan_specie, t_megan_ptrs
    USE mo_jsb_control,        ONLY: jsbach_runs_standalone
    ! TODO: These are likely specific to the DWD-NWP physics and do likely not work with ECHAM physics
    USE mo_nwp_phy_types,      ONLY: t_nwp_phy_diag
    USE mo_nwp_phy_state,      ONLY: prm_diag
    ! And these likely only for ECHAM physics
!    USE mo_echam_phy_memory,   ONLY: prm_field
    USE mo_run_config,         ONLY: iforcing
    USE mo_impl_constants,     ONLY: inwp

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    ! Local variables
    !
    TYPE(t_jsb_model), POINTER :: model
    dsl4jsb_Def_config(MEGAN_)
    dsl4jsb_Def_memory(MEGAN_)
    dsl4jsb_Def_memory(A2L_)
    dsl4jsb_Def_memory(HYDRO_)
    dsl4jsb_Def_memory(PHENO_)

    ! Pointers to variables in memory
    !
    ! Inputs
    dsl4jsb_Real2D_onChunk :: t_air
    dsl4jsb_Real2D_onChunk :: CO2_air
    dsl4jsb_Real2D_onChunk :: cos_zenith_angle
    dsl4jsb_Real2D_onChunk :: swvis_srf_down
    dsl4jsb_Real2D_onChunk :: swnir_srf_down
    dsl4jsb_Real2D_onChunk :: sw_toa_down
    dsl4jsb_Real2D_onChunk :: lai_ta
    dsl4jsb_Real2D_onChunk :: wtr_rootzone_rel
    dsl4jsb_Real2D_onChunk :: temp_relax_short
    dsl4jsb_Real2D_onChunk :: lai_relax_short
    dsl4jsb_Real2D_onChunk :: pfd_relax_short
    dsl4jsb_Real2D_onChunk :: lai_relax_long
    dsl4jsb_Real2D_onChunk :: pfd_relax_long
    dsl4jsb_Real2D_onChunk :: gamma_lai
    dsl4jsb_Real2D_onChunk :: gamma_pfd
    dsl4jsb_Real2D_onChunk :: gamma_CO2
    dsl4jsb_Real2D_onChunk :: gamma_sm

    TYPE (t_megan_ptrs) :: &
      & iso,               &
      & CO,                &
      & NO

    TYPE (t_megan_specie), POINTER :: &
      & isos,                         &
      & COs,                          &
      & NOs

    !--- local declarations
    TYPE (t_nwp_phy_diag), POINTER :: atm_diag
    CLASS(t_jsb_tile_abstract), POINTER :: ptr_tile
    TYPE(t_jsb_grid), POINTER :: grid
    TYPE(t_jsb_vgrid), POINTER :: soil_e, soil_w
    LOGICAL  :: lstart, lrestart, lnew_year, lnew_month, lnew_day, lstop
    INTEGER  :: tspd, iblk, ics, ice, nc, nsoil, ic
    INTEGER  :: year, month, summer_state, ylen, mlen
    REAL(lp) :: dtime

    !    Parameters
    REAL(lp), PARAMETER :: Wm2_pfd = 4.766_lp          ! factor to convert W to umol s-1
    REAL(lp), PARAMETER :: rho_iso = 0.96_lp           ! MEGAN model parameter
    REAL(lp), PARAMETER :: cce = 1.00_lp               ! MEGAN canopy model factor - not yet used
    REAL(lp), PARAMETER :: unit_scale = 1.e-6_lp / 3600._lp ! Scaling from [ug h-1] to [g s-1]

    !--- local variables ---

    CHARACTER(len=*), PARAMETER :: routine = modname//':update_megan'

    IF (.NOT. tile%Is_process_calculated(MEGAN_)) RETURN

    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')
    lstart     = is_time_experiment_start(options%current_datetime)
    lrestart   = is_time_restart(options%current_datetime)

    dtime = options%dtime

    iblk  = options%iblk
    ics   = options%ics
    ice   = options%ice
    nc    = options%nc
    tspd  = timesteps_per_day(dtime)

    model  => Get_model(tile%owner_model_id)
    grid   => Get_grid(model%grid_id)
    lStop  = .FALSE.

    dsl4jsb_Get_Config(MEGAN_)

    ! Get reference to variables for current block
    !
    dsl4jsb_Get_memory(MEGAN_)
    dsl4jsb_Get_memory(A2L_)
    dsl4jsb_Get_memory(HYDRO_)
    dsl4jsb_Get_memory(PHENO_)

    ! A2L_ variables
    dsl4jsb_Get_var2D_onChunk(A2L_,   t_air)
    dsl4jsb_Get_var2D_onChunk(A2L_,   swvis_srf_down)
    dsl4jsb_Get_var2D_onChunk(A2L_,   swnir_srf_down)
    ! PHENO_ variables
    dsl4jsb_Get_var2D_onChunk(PHENO_, lai_ta)
    ! MEGAN_ variables for running means
    dsl4jsb_Get_var2D_onChunk(MEGAN_, temp_relax_short)
    dsl4jsb_Get_var2D_onChunk(MEGAN_,  lai_relax_short)
    dsl4jsb_Get_var2D_onChunk(MEGAN_,  pfd_relax_short)
    dsl4jsb_Get_var2D_onChunk(MEGAN_,  lai_relax_long)
    dsl4jsb_Get_var2D_onChunk(MEGAN_,  pfd_relax_long)

    IF (lstart .AND. (dsl4jsb_Config(MEGAN_)%dbg_process)) CALL message(routine,'First instance of update_megan')

    ! Get needed time state
    lnew_year  = is_newyear (options%current_datetime,dtime)
    lnew_month = is_newmonth(options%current_datetime,dtime)
    lnew_day   = is_newday  (options%current_datetime,dtime)

    ! Running means must be updated every timestep
    DO ic = 1, nc
      CALL update_megan_runmean(lstart,dtime,                                                                    &
        &                       dsl4jsb_Config(MEGAN_)%relax_short,                                              &
        &                       dsl4jsb_Config(MEGAN_)%relax_long,                                               &
        &                                              t_air(ic), temp_relax_short(ic),                          &
        &                                              lai_ta(ic),  lai_relax_short(ic), lai_relax_long(ic),     &
        &                       swvis_srf_down(ic)+swnir_srf_down(ic),  pfd_relax_short(ic), pfd_relax_long(ic))
    ENDDO

    ! A2L_ variables
    dsl4jsb_Get_var2D_onChunk(A2L_,   cos_zenith_angle)
    dsl4jsb_Get_var2D_onChunk(A2L_,   CO2_air)
    ! HYDRO_ variables
    dsl4jsb_Get_var2D_onChunk(HYDRO_, wtr_rootzone_rel)
    ! MEGAN_ variables
    dsl4jsb_Get_var2D_onChunk(MEGAN_, iso%total_out)
    dsl4jsb_Get_var2D_onChunk(MEGAN_, iso%gamma_age)
    dsl4jsb_Get_var2D_onChunk(MEGAN_, iso%gamma_temp)
    dsl4jsb_Get_var2D_onChunk(MEGAN_,  CO%total_out)
    dsl4jsb_Get_var2D_onChunk(MEGAN_,  CO%gamma_age)
    dsl4jsb_Get_var2D_onChunk(MEGAN_,  CO%gamma_temp)
    dsl4jsb_Get_var2D_onChunk(MEGAN_,  NO%total_out)
    dsl4jsb_Get_var2D_onChunk(MEGAN_, gamma_lai)
    dsl4jsb_Get_var2D_onChunk(MEGAN_, gamma_pfd)
    dsl4jsb_Get_var2D_onChunk(MEGAN_, gamma_CO2)
    dsl4jsb_Get_var2D_onChunk(MEGAN_, gamma_sm)
    IF (jsbach_runs_standalone()) THEN
      ! @TODO: Wrong! needs to be solved differently for standalone runs, but this is a quick fix for now
      sw_toa_down => dsl4jsb_var2D_onChunk(A2L_, swvis_srf_down)
    ELSE IF (iforcing == inwp) THEN
      sw_toa_down => prm_diag(model%id)%flxdwswtoa(ics:ice,iblk)
!    ELSE ! ECHAM physics no longer supported
!      sw_toa_down => prm_field(model%id)%rsdt(ics:ice,iblk)
    ENDIF
    isos => dsl4jsb_memory(MEGAN_)%iso
    COs  => dsl4jsb_memory(MEGAN_)%CO
    NOs  => dsl4jsb_memory(MEGAN_)%NO

    !--- calculate vegetation emissions
    !    Isoprene:
    !         Emission = [EF][GAMMA_LAI][GAMMA_P][GAMMA_T][GAMMA_age][GAMMA_SM]
    !
    !    non Isoprene:
    !         Emission = [EF][GAMMA_LAI][GAMMA_T][GAMMA_age][GAMMA_SM]*
    !                    { (1-LDF) + [LDF][GAMMA_P] }
    !--- first compute gamma factors independent of compound
    DO ic = 1, nc
      CALL calc_gamma_lai (lai_ta(ic), gamma_lai(ic))
      CALL calc_gamma_pfd(cos_zenith_angle(ic), swvis_srf_down(ic)+swnir_srf_down(ic), sw_toa_down(ic), &
                          pfd_relax_short(ic), gamma_pfd(ic))
      CALL calc_gamma_co2(CO2_air(ic), gamma_CO2(ic))
      CALL calc_gamma_sm(wtr_rootzone_rel(ic), gamma_sm(ic))

      ! Isoprene (only one using gamma_CO2 and gamma_sm)
      CALL calc_gamma_age(isos%cafac, lai_relax_short(ic), lai_relax_long(ic), temp_relax_short(ic), iso%gamma_age(ic))
      CALL calc_gamma_temp_iso(t_air(ic), temp_relax_short(ic), iso%gamma_temp(ic))
      iso%total_out(ic) = isos%ef * gamma_lai(ic) * iso%gamma_temp(ic) * iso%gamma_age(ic) * gamma_CO2(ic) * gamma_sm(ic) * &
                          (1._lp - isos%ldf + isos%ldf * gamma_pfd(ic)) * unit_scale
      ! CO
      CALL calc_gamma_age(COs%cafac, lai_relax_short(ic), lai_relax_long(ic), temp_relax_short(ic), CO%gamma_age(ic))
      CALL calc_gamma_temp_other(t_air(ic), COs%beta, CO%gamma_temp(ic))
      CO%total_out(ic) = COs%ef * gamma_lai(ic) * CO%gamma_temp(ic) * CO%gamma_age(ic) * &
                         (1._lp - COs%ldf + COs%ldf * gamma_pfd(ic)) * unit_scale
      ! NO
      NO%total_out(ic) = NOs%ef * gamma_lai(ic) * (1._lp - NOs%ldf + NOs%ldf * gamma_pfd(ic)) * unit_scale
    ENDDO

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE update_megan

  SUBROUTINE aggregate_megan(tile, options)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    TYPE(t_jsb_model), POINTER :: model
    dsl4jsb_Def_config(MEGAN_)
    dsl4jsb_Def_memory(MEGAN_)

    CLASS(t_jsb_aggregator), POINTER :: weighted_by_fract

    CHARACTER(len=*), PARAMETER :: routine = modname//':aggregate_megan'

    INTEGER :: iblk, ics, ice

    iblk = options%iblk
    ics  = options%ics
    ice  = options%ice

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    model => Get_model(tile%owner_model_id)

    dsl4jsb_Get_config(MEGAN_)
    dsl4jsb_Get_memory(MEGAN_)

    weighted_by_fract => tile%Get_aggregator("weighted_by_fract")

    dsl4jsb_Aggregate_onChunk(MEGAN_, iso%total_out,    weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(MEGAN_, iso%gamma_age,    weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(MEGAN_, iso%gamma_temp,   weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(MEGAN_,  CO%total_out,    weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(MEGAN_,  CO%gamma_age,    weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(MEGAN_,  CO%gamma_temp,   weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(MEGAN_,  NO%total_out,    weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(MEGAN_, gamma_lai,        weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(MEGAN_, gamma_pfd,        weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(MEGAN_, gamma_CO2,        weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(MEGAN_, gamma_sm,         weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(MEGAN_, temp_relax_short, weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(MEGAN_,  lai_relax_short, weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(MEGAN_,  pfd_relax_short, weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(MEGAN_,  lai_relax_long,  weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(MEGAN_,  pfd_relax_long,  weighted_by_fract)

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE aggregate_megan

#endif
END MODULE mo_megan_interface
