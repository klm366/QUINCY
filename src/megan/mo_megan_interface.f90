!>
!! @brief Contains the interfaces to the megan scheme
!!
!! @author
!!  Stiig Wilkenskjeld, MPI-M Hamburg
!!
!! @par Revision History
!! Adopted from JSBACH3 code by Thomas Kleinen, 2024-01-24
!!
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/license .
!!

MODULE mo_megan_interface
#ifndef __NO_JSBACH__

  USE mo_kind,                   ONLY: wp
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
  !! Constructors for tasks

  ! -------------------------------------------------------------------------------------------------------
  !> Constructor for megan
  !!
  !! @param[in]     model_id     Model id
  !! @return        return_ptr   Instance of process task "topmodel"
  !!
  FUNCTION Create_task_megan(model_id) RESULT(return_ptr)

    INTEGER,                   INTENT(in) :: model_id
    CLASS(t_jsb_process_task), POINTER    :: return_ptr

    ALLOCATE(tsk_megan::return_ptr)
    CALL return_ptr%Construct(name='megan_processes', process_id=MEGAN_, owner_model_id=model_id)

  END FUNCTION Create_task_megan

  !-----------------------------------------------------------------------------------------------------
  !> Register tasks: MEGAN_ process
  !!
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
  !!
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE update_megan(tile, options)

    USE mo_jsb_time,           ONLY: is_time_experiment_start, is_time_restart, is_newyear, is_newmonth, is_newday, &
      &                              timesteps_per_day, get_date_components, get_year_length, get_month_length
    USE mo_megan_memory_class, ONLY: t_megan_specie, t_megan_ptrs
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
    INTEGER  :: tspd, iblk, ics, ice, nc, nsoil, i
    INTEGER  :: year, month, summer_state, ylen, mlen
    REAL(wp) :: dtime

    !    Parameters
    REAL(wp), PARAMETER :: Wm2_pfd = 4.766_wp          ! factor to convert W to umol s-1
    REAL(wp), PARAMETER :: rho_iso = 0.96_wp           ! MEGAN model parameter
    REAL(wp), PARAMETER :: cce = 1.00_wp               ! MEGAN canopy model factor - not yet used
    REAL(wp), PARAMETER :: unit_scale = 1.e-6_wp / 3600._wp ! Scaling from [ug h-1] to [g s-1]

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
    CALL update_megan_runmean(lstart,dtime,                                                    &
      &                       dsl4jsb_Config(MEGAN_)%relax_short,                              &
      &                       dsl4jsb_Config(MEGAN_)%relax_long,                               &
      &                                               t_air, temp_relax_short,                 &
      &                                              lai_ta,  lai_relax_short, lai_relax_long, &
      &                       swvis_srf_down+swnir_srf_down,  pfd_relax_short, pfd_relax_long)

    IF (.NOT. lnew_day) RETURN ! Megan process is updated only once per day TODO: Check this!

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
!    IF (iforcing == inwp) THEN
      sw_toa_down => prm_diag(model%id)%flxdwswtoa(ics:ice,iblk)
!    ELSE ! ECHAM physics no longer supported
!      sw_toa_down => prm_field(model%id)%rsdt(ics:ice,iblk)
!    ENDIF
    isos => dsl4jsb_memory(MEGAN_)%iso
     COs => dsl4jsb_memory(MEGAN_)%CO
     NOs => dsl4jsb_memory(MEGAN_)%NO

    !--- calculate vegetation emissions
    !    Isoprene:
    !         Emission = [EF][GAMMA_LAI][GAMMA_P][GAMMA_T][GAMMA_age][GAMMA_SM]
    !
    !    non Isoprene:
    !         Emission = [EF][GAMMA_LAI][GAMMA_T][GAMMA_age][GAMMA_SM]*
    !                    { (1-LDF) + [LDF][GAMMA_P] }
    !--- first compute gamma factors independent of compound
    CALL calc_gamma_lai (lai_ta, gamma_lai)
    CALL calc_gamma_pfd(cos_zenith_angle, swvis_srf_down+swnir_srf_down, sw_toa_down, pfd_relax_short, gamma_pfd)
    CALL calc_gamma_co2(CO2_air, gamma_CO2)
    CALL calc_gamma_sm(wtr_rootzone_rel, gamma_sm)

    ! Isoprene (only one using gamma_CO2 and gamma_sm)
    CALL calc_gamma_age(isos%cafac, lai_relax_short, lai_relax_long, temp_relax_short, iso%gamma_age)
    CALL calc_gamma_temp_iso(t_air, temp_relax_short, iso%gamma_temp)
    iso%total_out = isos%ef * gamma_lai * iso%gamma_temp * iso%gamma_age * gamma_CO2 * gamma_sm * &
                     (1._wp - isos%ldf + isos%ldf * gamma_pfd) * unit_scale
    ! CO
    CALL calc_gamma_age(COs%cafac, lai_relax_short, lai_relax_long, temp_relax_short, CO%gamma_age)
    CALL calc_gamma_temp_other(t_air(:), COs%beta, CO%gamma_temp)
    CO%total_out = COs%ef * gamma_lai * CO%gamma_temp * CO%gamma_age * &
                    (1._wp - COs%ldf + COs%ldf * gamma_pfd) * unit_scale
    ! NO
    NO%total_out = NOs%ef * gamma_lai * (1._wp - NOs%ldf + NOs%ldf * gamma_pfd) * unit_scale

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
