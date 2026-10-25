!> Contains the interfaces to methane production and transport
!>
!> ICON-Land
!>
!> ---------------------------------------
!> Copyright (C) 2013-2026, MPI-M, MPI-BGC
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
!> Adopted from JSBACH3 code by Thomas Kleinen, 2023-03-17
!>

MODULE mo_methane_interface
#ifndef __NO_JSBACH__

  USE mo_kind,                   ONLY: lp, i4
  USE mo_exception,              ONLY: message, finish, message_text
  USE mo_jsb_control,            ONLY: debug_on
  USE mo_jsb_model_class,        ONLY: t_jsb_model
  USE mo_jsb_class,              ONLY: Get_model
  USE mo_jsb_grid_class,         ONLY: t_jsb_grid, t_jsb_vgrid
  USE mo_jsb_grid,               ONLY: Get_vgrid, Get_grid
  USE mo_jsb_tile_class,         ONLY: t_jsb_tile_abstract, t_jsb_aggregator
  USE mo_jsb_process_class,      ONLY: t_jsb_process, PHENO_, CH4_, TOP_, CARBON_, A2L_, CLIMBUF_, HYDRO_, SSE_, NLCC_, MEGAN_
  USE mo_jsb_task_class,         ONLY: t_jsb_process_task, t_jsb_task_options
  USE mo_jsb_math_constants,     ONLY: eps8
  USE mo_jsb_physical_constants, ONLY: tmelt
  USE mo_jsb_time,               ONLY: t_datetime, is_newmonth, get_time_next, is_newyear, timesteps_per_day, &
                                       is_time_experiment_start, get_time_start, is_time_restart, get_date_components
  USE mo_jsb_io_netcdf,          ONLY: jsb_netcdf_open_input
  USE mo_mpi,                    ONLY: p_bcast, p_io, my_process_is_stdio
  USE mo_ch4_memory_class,       ONLY: t_chem_ptrs, t_timeseries
  USE mo_methane_constants,      ONLY: t_specie_const, CH4_consts, CO2_consts, O2_consts, emission_factors_for_methane,  &
  &                                    LT_TROPICAL_FOREST, LT_EXTRATROPICAL_FOREST, LT_SAVANNAH, LT_AGRICULTURE, LT_PEAT
  USE mo_methane_util,           ONLY: concentration_mmr_2_mol_m3, concentration_vmr_2_mol_m3, concentration_mol_m3_2_vmr, &
                                       determine_Henry_cc
  USE mo_methane_process,        ONLY: plant_setup, calc_methane_termites, calc_methane_herbivores, calc_methane_fire, &
                                       soil_column_setup, diffusion_setup, production, plant_transport, oxidation_mm,  &
                                       ebullition, init_gas_pools
! TODO: Emergency "solution" to allow compilation for GPU
#ifndef _OPENACC
  USE mo_methane_process,        ONLY: reaction_diffusion
#endif

  USE mo_ch4_config_class,       ONLY: t_ch4_source_config, t_multi_file
  USE mo_util,                   ONLY: real2string
  USE mo_megan_memory_class,     ONLY: t_megan_ptrs
  USE mo_debug_stw

  dsl4jsb_Use_config(CH4_)
  dsl4jsb_Use_config(MEGAN_)
  dsl4jsb_Use_config(TOP_)
  dsl4jsb_Use_config(SSE_)
  dsl4jsb_Use_config(NLCC_)
  dsl4jsb_Use_memory(CH4_)
  dsl4jsb_Use_memory(MEGAN_)
  dsl4jsb_Use_memory(TOP_)
  dsl4jsb_Use_memory(SSE_)
  dsl4jsb_Use_memory(CARBON_)
  dsl4jsb_Use_memory(A2L_)
  dsl4jsb_Use_memory(CLIMBUF_)
  dsl4jsb_Use_memory(PHENO_)
  dsl4jsb_Use_memory(HYDRO_)

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: Register_methane_tasks

  !> Type definition for soil_temp_and_hydro task
  TYPE, EXTENDS(t_jsb_process_task) :: tsk_methane
  CONTAINS
    PROCEDURE, NOPASS :: Integrate => update_methane
    PROCEDURE, NOPASS :: Aggregate => aggregate_methane
  END TYPE tsk_methane

  !> Constructor interface for methane task
  INTERFACE tsk_methane
    PROCEDURE Create_task_methane
  END INTERFACE tsk_methane

  CHARACTER(len=*), PARAMETER :: modname = 'mo_methane_interface'

  ! StW: These variables are shared between tiles - not very nice! Reconsider!
  REAL(lp), POINTER, SAVE :: frac_non_crop(:,:) => NULL()
  REAL(lp), POINTER, SAVE :: sum_non_crop (:,:) => NULL()

CONTAINS

  ! ================================================================================================================================
  !> Constructors for tasks

  ! -------------------------------------------------------------------------------------------------------
  !> Constructor for methane
  !>
  !> param[in]     model_id     Model id
  !> return        return_ptr   Instance of process task "methane"
  !>
  FUNCTION Create_task_methane(model_id) RESULT(return_ptr)

    INTEGER,                   INTENT(in) :: model_id
    CLASS(t_jsb_process_task), POINTER    :: return_ptr

    ALLOCATE(tsk_methane::return_ptr)
    CALL return_ptr%Construct(name='methane_processes', process_id=CH4_, owner_model_id=model_id)

  END FUNCTION Create_task_methane

  !-----------------------------------------------------------------------------------------------------
  !> Register tasks: CH4_ process
  !>
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE Register_methane_tasks(process, model_id)

    USE mo_jsb_model_class,   ONLY: t_jsb_model
    USE mo_jsb_class,         ONLY: Get_model

    CLASS(t_jsb_process), INTENT(inout) :: process
    INTEGER,              INTENT(in)    :: model_id

    CALL process%Register_task(tsk_methane(model_id))

  END SUBROUTINE Register_methane_tasks

  SUBROUTINE reset_chem_specie(specie,nc,nsoil,dbg)

    TYPE (t_chem_ptrs), INTENT(inout) :: specie
    INTEGER,            INTENT(in)    :: nc, nsoil
    LOGICAL, OPTIONAL,  INTENT(in)    :: dbg

    INTEGER :: ic, is

    !$ACC PARALLEL LOOP DEFAULT(PRESENT) ASYNC(1) GANG VECTOR
    DO ic = 1, nc
      specie%diff_out   (ic) = 0._lp
      specie%prod_out   (ic) = 0._lp
      specie%plant_out  (ic) = 0._lp
      specie%oxid_out   (ic) = 0._lp
      specie%in_out     (ic) = 0._lp
      specie%total_out  (ic) = 0._lp
      specie%total_final(ic) = 0._lp
!      IF (ASSOCIATED(specie%other_out)) THEN
!        specie%other_out(ic)   = 0._lp
!      ENDIF
    ENDDO
    !$ACC END PARALLEL LOOP
    IF (PRESENT(dbg)) THEN
      IF (dbg) THEN
        !$ACC PARALLEL LOOP DEFAULT(PRESENT) ASYNC(1) GANG VECTOR
        DO ic = 1, nc
          specie%test2d(ic)       = 0._lp
        ENDDO
        !$ACC END PARALLEL LOOP
        !$ACC PARALLEL LOOP DEFAULT(PRESENT) ASYNC(2) COLLAPSE(2) GANG VECTOR
        DO is = 1, nsoil
          DO ic = 1, nc
            specie%prod_out_sl (ic,is) = 0._lp
            specie%plant_out_sl(ic,is) = 0._lp
            specie%oxid_out_sl (ic,is) = 0._lp
!            specie%test3d(ic,is)       = 0._lp
          ENDDO
        ENDDO
        !$ACC END PARALLEL LOOP
      ENDIF
    ENDIF
    !$ACC WAIT

  END SUBROUTINE reset_chem_specie

  SUBROUTINE update_chem_specie(dt,nly,i,frac,specie,diff,act_sl,prod_sl,plant_sl,oxid_sl,dbg,other_sl,to_other_i)

    REAL(lp),           INTENT(in)    :: dt
    INTEGER,            INTENT(in)    :: nly
    INTEGER,            INTENT(in)    :: i
    REAL(lp),           INTENT(in)    :: frac
    TYPE (t_chem_ptrs), INTENT(inout) :: specie
    REAL(lp),           INTENT(in)    :: diff
    REAL(lp),           INTENT(in)    :: act_sl(:)
    REAL(lp),           INTENT(in)    :: prod_sl(:)
    REAL(lp),           INTENT(in)    :: plant_sl(:)
    REAL(lp),           INTENT(in)    :: oxid_sl(:)
    LOGICAL,  OPTIONAL, INTENT(in)    :: dbg
    REAL(lp), OPTIONAL, INTENT(in)    :: other_sl(:)
    INTEGER,  OPTIONAL, INTENT(in)    :: to_other_i

    INTEGER  :: j
    REAL(lp) :: mul

    mul = frac / dt

    specie%dens_sl    (i,:) =     act_sl  ( :   )  * frac
    specie%total_final(i)   = SUM(act_sl  (1:nly)) * frac
    specie%prod_out   (i)   = SUM(prod_sl (1:nly)) * mul
    specie%plant_out  (i)   = SUM(plant_sl(1:nly)) * mul
    specie%oxid_out   (i)   = SUM(oxid_sl (1:nly)) * mul
    specie%diff_out   (i)   =     diff             * mul
    specie%total_out  (i)   = specie%diff_out(i) + specie%plant_out(i)
    IF (PRESENT(other_sl)) THEN
      specie%other_out(i)   = SUM(other_sl(1:nly)) * mul
      IF (PRESENT(to_other_i)) THEN
        IF (to_other_i /= 0) specie%total_out(i) = specie%total_out(i) + specie%other_out(i)
      ENDIF
    ENDIF
    specie%in_out   (i) = MAX(0._lp,specie%total_out(i))
    specie%total_out(i) = MIN(0._lp,specie%total_out(i))
    IF (PRESENT(dbg)) THEN
      IF (dbg) THEN
        DO j = 1, nly
          specie%prod_out_sl (i,j) = prod_sl (j) * mul
          specie%plant_out_sl(i,j) = plant_sl(j) * mul
          specie%oxid_out_sl (i,j) = oxid_sl (j) * mul
          IF (PRESENT(other_sl)) THEN
            specie%other_out_sl(i,j) = other_sl(j) * mul
          ENDIF
        ENDDO
      ENDIF
    ENDIF

  END SUBROUTINE update_chem_specie

  SUBROUTINE methane_update_atm_GHG_conc(options,GAS,source,atm_conc,nml_conc,GHG_atm,dst,dbg)

    TYPE(t_jsb_task_options),   INTENT(in)    :: options
    TYPE(t_specie_const),       INTENT(in)    :: GAS
    INTEGER,                    INTENT(in)    :: source
    REAL(lp),                   INTENT(in)    :: atm_conc(:)
    REAL(lp),                   INTENT(in)    :: nml_conc
    TYPE(t_timeseries),         INTENT(in)    :: GHG_atm
    REAL(lp),                   INTENT(inout) :: dst(:)
    LOGICAL,                    INTENT(in)    :: dbg

    LOGICAL :: lstart, lrestart, lyear, lupdate
    INTEGER :: src, ic, nc, iy, year
    LOGICAL, SAVE :: First = .TRUE.

    CHARACTER(len=*), PARAMETER :: routine = modname//':methane_update_atm_GHG_conc'

    src = MIN(source,3) ! Different file specifications are already converted
    nc  = SIZE(dst)

    lupdate = .FALSE.
    SELECT CASE (src)
      CASE (1) ! Constant from namelist - just copy if not already done
        IF (ABS(dst(1)-nml_conc) > 1.e-10_lp) THEN ! Only need to do this once
          !$ACC PARALLEL LOOP DEFAULT(PRESENT) ASYNC(1) GANG VECTOR
          DO ic = 1, nc
            dst(ic) = nml_conc
          ENDDO
          !$ACC END PARALLEL LOOP
          lupdate = .TRUE.
        ENDIF
      CASE (2) ! From atmospheric model - update every timestep
        !$ACC PARALLEL LOOP DEFAULT(PRESENT) ASYNC(1) GANG VECTOR
        DO ic = 1, nc
          dst(ic) = concentration_mmr_2_mol_m3(GAS,atm_conc(ic))
        ENDDO
        !$ACC END PARALLEL LOOP
      CASE (3) ! From timeseries file - update at year start
        lstart   = is_time_experiment_start(options%current_datetime)
        lrestart = is_time_restart(options%current_datetime)
        lyear    = is_newyear(options%current_datetime,options%dtime)
        IF (lrestart .OR. lstart .OR. lyear) THEN
          CALL get_date_components(options%current_datetime, year=year)
          iy = MAX(1,MIN(SIZE(GHG_atm%atm_vmr),year - GHG_atm%atm_year_min + 1))
          dst(1) = concentration_vmr_2_mol_m3(GAS,GHG_atm%atm_vmr(iy))
          !$ACC PARALLEL LOOP DEFAULT(PRESENT) ASYNC(1) GANG VECTOR
          DO ic = 2, nc
            dst(ic) = dst(1)
          ENDDO
          !$ACC END PARALLEL LOOP
          lupdate = .TRUE.
        ENDIF
    END SELECT
    !$ACC WAIT(1)
    IF ((lupdate .OR. First) .AND. dbg) &
      CALL message(routine,'Updated atmospheric '//TRIM(GAS%specie_name)//' vmr: '// &
                           real2string(concentration_mol_m3_2_vmr(GAS,dst(1))))
    First = .FALSE.

  END SUBROUTINE methane_update_atm_GHG_conc

  !-----------------------------------------------------------------------------------------------------
  !> update_methane
  !>
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE update_methane(tile, options)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    ! Local variables
    !
    CLASS(t_jsb_tile_abstract), POINTER :: curr_tile
    TYPE(t_ch4_source_config), POINTER :: curr
    TYPE(t_jsb_model), POINTER :: model
    LOGICAL, SAVE :: lfirstcall = .TRUE.
    dsl4jsb_Def_config(CH4_)
    dsl4jsb_Def_config(MEGAN_)
    dsl4jsb_Def_config(TOP_)
    dsl4jsb_Def_config(SSE_)
    dsl4jsb_Def_config(NLCC_)
    dsl4jsb_Def_memory(CARBON_)
    dsl4jsb_Def_memory(CH4_)
    dsl4jsb_Def_memory(MEGAN_)
    dsl4jsb_Def_memory(TOP_)
    dsl4jsb_Def_memory(SSE_)
    dsl4jsb_Def_memory(A2L_)
    dsl4jsb_Def_memory(CLIMBUF_)
    dsl4jsb_Def_memory(PHENO_)
    dsl4jsb_Def_memory(HYDRO_)

    ! Multivariables defined for chemical species
    TYPE (t_chem_ptrs) :: &
      & CH4_dry,          &
      & CH4_wet,          &
!      & CH4_tot,          &
      & CO2_dry,          &
      & CO2_wet,          &
      & O2_dry,           &
      & O2_wet

    TYPE (t_megan_ptrs) :: &
      & co,                &
      & no,                &
      & iso

    ! Pointers to variables in memory
    ! Inputs
    dsl4jsb_Real2D_onChunk :: CH4_air     ! Should be from interface but not yet available
    dsl4jsb_Real2D_onChunk :: CO2_air     ! Values from atmosphere
    dsl4jsb_Real2D_onChunk :: O2_air      ! Should be from interface but not yet available
    dsl4jsb_Real2D_onChunk :: CO2_conc
    dsl4jsb_Real2D_onChunk :: CH4_conc
    dsl4jsb_Real2D_onChunk :: O2_conc
    dsl4jsb_Real2D_onChunk :: press_srf
    dsl4jsb_Real2D_onChunk :: soil_respiration
    dsl4jsb_Real2D_onChunk :: soil_respiration_oxic_ta
    dsl4jsb_Real2D_onChunk :: soil_respiration_anoxic_ta
    dsl4jsb_Real2D_onChunk :: acrotelm_height_ag
    dsl4jsb_Real2D_onChunk :: acrotelm_height_bg
    dsl4jsb_Real2D_onChunk :: mxlay_above_rock
    dsl4jsb_Real2D_onChunk :: lai
    dsl4jsb_Real2D_onChunk :: frac_inun
    dsl4jsb_Real2D_onChunk :: frac_inun_pre
    dsl4jsb_Real2D_onChunk :: peat_water_table
    dsl4jsb_Real2D_onChunk :: wt_inun
    dsl4jsb_Real2D_onChunk :: soil_depth
    dsl4jsb_Real2D_onChunk :: vol_porosity
    dsl4jsb_Real2D_onChunk :: bclapp
    dsl4jsb_Real3D_onChunk :: wtr_soil_sl
    dsl4jsb_Real3D_onChunk :: ice_soil_sl
    dsl4jsb_Real3D_onChunk :: t_soil_sl
    dsl4jsb_Real3D_onChunk :: root_dist_sl
    dsl4jsb_Real2D_onChunk :: c_green
    dsl4jsb_Real2D_onChunk :: CH4_out_net
    dsl4jsb_Real2D_onChunk :: CO2_out_net
    dsl4jsb_Real2D_onChunk :: CH4_prod_ratio_out_wet
    dsl4jsb_Real2D_onChunk :: co2flux_herb_2_atm_ta
    dsl4jsb_Real2D_onChunk :: co2flux_fire_all_2_atm_ta
    dsl4jsb_Real2D_onChunk :: CH4_out_fire
    dsl4jsb_Real2D_onChunk :: CH4_out_termite
    dsl4jsb_Real2D_onChunk :: CH4_out_herbivore
    dsl4jsb_Real2D_onChunk :: fire_em_NO
    dsl4jsb_Real2D_onChunk :: fire_em_CO
    dsl4jsb_Real2D_onChunk :: fire_em_iso
    dsl4jsb_Real2D_onChunk :: termite_biomass
    dsl4jsb_Real2D_onChunk :: GPP_prev_year_ta
    dsl4jsb_Real2D_onChunk :: GDD_prev_year
    dsl4jsb_Real2D_onChunk :: GDD_sum
    dsl4jsb_Real2D_onChunk :: min_mmtemp_20yr
    dsl4jsb_Real2D_onChunk :: ch4_2_atm
    dsl4jsb_Real2D_onChunk :: co_2_atm
    dsl4jsb_Real2D_onChunk :: no_2_atm
    dsl4jsb_Real2D_onChunk :: iso_2_atm

    ! Pointers to allocated non-variables
    !
    REAL(lp), CONTIGUOUS, POINTER :: OxicDecompCarbon(:)
    REAL(lp), CONTIGUOUS, POINTER :: AnoxDecompCarbon(:)
    REAL(lp), CONTIGUOUS, POINTER :: water_table(:)
    REAL(lp), CONTIGUOUS, POINTER :: tiller_area_sl(:)
    REAL(lp), CONTIGUOUS, POINTER :: frac_tile(:), frac_tile_abs(:), frac1(:), tile_fac(:)
    REAL(lp), CONTIGUOUS, POINTER :: D_CH4(:), D_CO2(:), D_O2(:), D_used(:)
    REAL(lp), CONTIGUOUS, POINTER :: CH4_dry_sl(:), CO2_dry_sl(:), O2_dry_sl(:)
    REAL(lp), CONTIGUOUS, POINTER :: CH4_wet_sl(:), CO2_wet_sl(:), O2_wet_sl(:)
    REAL(lp), CONTIGUOUS, POINTER :: CH4_plant_sl(:), CO2_plant_sl(:), O2_plant_sl(:)
    REAL(lp), CONTIGUOUS, POINTER :: CH4_prod_sl(:), CO2_prod_sl(:), O2_prod_sl(:)
    REAL(lp), CONTIGUOUS, POINTER :: CH4_oxid_sl(:), CO2_oxid_sl(:), O2_oxid_sl(:)
!    REAL(lp), CONTIGUOUS, POINTER :: CH4_excess_sl(:), CO2_excess_sl(:), O2_excess_sl(:)
    REAL(lp), CONTIGUOUS, POINTER :: CH4_ebul_sl(:), O2_other_sl(:)
    REAL(lp), CONTIGUOUS, POINTER :: relative_moist_sl(:), relative_ice_sl(:)
    REAL(lp), CONTIGUOUS, POINTER :: relative_porevol_sl(:), relative_pore_air_sl(:), relative_pore_moist_sl(:)
    REAL(lp), CONTIGUOUS, POINTER :: dens2conc_sl(:), conc2dens_sl(:)
    REAL(lp), CONTIGUOUS, POINTER :: CH4_act(:), CO2_act(:), O2_act(:)
    ! Pointers for assignment
    REAL(lp), CONTIGUOUS, POINTER :: GDD(:)

    TYPE(t_datetime), POINTER :: time_next
    TYPE(t_jsb_grid), POINTER :: grid
    TYPE(t_jsb_vgrid), POINTER :: soil_e, soil_w
    REAL(lp) :: stored_C_sum
    LOGICAL  :: lnew_day, lstart, lrestart, labs_tile_calc
    INTEGER  :: year, hour, minute, sec
    LOGICAL  :: is_veg, is_wet, is_peat, compute_flag, l_stop, l_GPP_dep
    LOGICAL  :: lissue, lissue_tst, dbg, lmsg
    INTEGER  :: tspd, iblk, ics, ice, nc, nsoil, lt, ic, is, iy
    INTEGER  :: iWetness, iWetness_st, iWetness_en
    REAL(lp) :: dt
    REAL(lp) :: acro_height, frac_inun_loc, frac
    REAL(lp) :: D_CH4_surface, D_CO2_surface, D_O2_surface
    REAL(lp) :: CH4_Prod, CO2_Prod, O2_Prod, CH4_prod_ratio
    REAL(lp) :: CH4_diff, CO2_diff, O2_diff
    INTEGER  :: mxlay_above_wt, mxlay_acrotelm, nly

    REAL(lp), TARGET :: empty(options%nc)

    CHARACTER(len=*), PARAMETER :: routine = modname//':update_methane'

    IF (.NOT. tile%Is_process_calculated(CH4_)) RETURN
! TODO: Emergency exit, since diffusion doesn't work with GPUs (using non-GPU LAPACK routine dgtvs)
#ifdef _OPENACC
    RETURN
#endif

    model => Get_model(tile%owner_model_id)
    dsl4jsb_Get_Config(CH4_)
    dsl4jsb_Get_Config(TOP_)
    dsl4jsb_Get_Config(SSE_)
    dsl4jsb_Get_Config(NLCC_)
    dsl4jsb_Get_Config(MEGAN_)
    dsl4jsb_Get_memory(CH4_)
    dsl4jsb_Get_memory(TOP_)
    dsl4jsb_Get_memory(SSE_)
    dsl4jsb_Get_memory(HYDRO_)
    dsl4jsb_Get_memory(CARBON_)
    dsl4jsb_Get_memory(A2L_)
    dsl4jsb_Get_memory(CLIMBUF_)
    dsl4jsb_Get_memory(PHENO_)
    dsl4jsb_Get_memory(MEGAN_)

    lstart = is_time_experiment_start(options%current_datetime)
    lrestart = is_time_restart(options%current_datetime)
    iblk  = options%iblk
    ics   = options%ics
    ice   = options%ice
    nc    = options%nc
    dt    = 86400._lp / REAL(dsl4jsb_Config(CH4_)%nstep_per_day,lp)
    tspd  = timesteps_per_day(dt)
    dbg   = dsl4jsb_Config(CH4_)%dbg_process .AND. .NOT. ASSOCIATED(tile%get_prev_sibling_tile()) .AND. (iblk==1)

    ! Atmospheric trace gase concentration needs to be updated here, since restarts doesn't start at a methane time step

    ! A2L_ variables
    dsl4jsb_Get_var2D_onChunk(A2L_, CO2_air)
    CH4_air => empty(:)
    O2_air => empty(:)
    IF (dsl4jsb_Config(CH4_)%CH4_source > 1 .OR. dsl4jsb_Config(CH4_)%O2_source > 1) THEN
      CALL finish(routine,'CH4 and O2 currently only implemented as constant value. Please use *_source=1')
    END IF
    ! CH4_ variables
    dsl4jsb_Get_var2D_onChunk(CH4_, CH4_conc)
    dsl4jsb_Get_var2D_onChunk(CH4_, CO2_conc)
    dsl4jsb_Get_var2D_onChunk(CH4_, O2_conc)
    CALL methane_update_atm_GHG_conc(options,CH4_consts,dsl4jsb_Config(CH4_)%CH4_source,CH4_air, &
                                     dsl4jsb_Config(CH4_)%CH4_conc,dsl4jsb_memory(CH4_)%GHG_atm_CH4,CH4_conc,dbg)
    CALL methane_update_atm_GHG_conc(options,CO2_consts,dsl4jsb_Config(CH4_)%CO2_source,CO2_air, &
                                     dsl4jsb_Config(CH4_)%CO2_conc,dsl4jsb_memory(CH4_)%GHG_atm_CO2,CO2_conc,dbg)
    CALL methane_update_atm_GHG_conc(options, O2_consts,dsl4jsb_Config(CH4_)%O2_source, O2_air, &
                                     dsl4jsb_Config(CH4_)%O2_conc, dsl4jsb_memory(CH4_)%GHG_atm_O2,O2_conc,dbg)

    CALL get_date_components(options%current_datetime, year=year, hour=hour, minute=minute, second=sec)
    sec = sec + minute*60 + hour*3600
    IF (MOD(sec,INT(dt)) /= 0) RETURN ! Not a methane time step

    ! TODO: Could this test be safely done somewhere in the init module?
    IF (dsl4jsb_Config(CH4_)%do_termite .AND. .NOT. dsl4jsb_Config(NLCC_)%active) THEN
      IF (iblk==1) &
        CALL message(routine,'Termite methane can only be calculated when the NLCC process is active. Termites are switched off')
      dsl4jsb_Config(CH4_)%do_termite = .FALSE.
    ENDIF

    grid   => Get_grid(model%grid_id)
    soil_e => Get_vgrid('soil_depth_energy')
    soil_w => Get_vgrid('soil_depth_water')
    nsoil  = soil_e%n_levels

    ! TODO: This test should be moved outside the time/blk loop
    l_stop  = .FALSE.
    IF (nsoil /= soil_w%n_levels) l_stop = .TRUE.
    IF (.NOT. l_stop) l_stop = ANY(ABS(soil_e%dz - soil_w%dz) > eps8)
    IF (l_stop) CALL finish(routine,'Sorry, present implementation of methane requires water and energy soil levels to coincide')

    ! Allocate global arrays shared between tiles. NOTICE: TODO: These are NEVER deallocated and shared among all tiles
    IF (.NOT. ASSOCIATED(frac_non_crop) .AND. dsl4jsb_Config(CH4_)%do_termite) THEN
      ALLOCATE(frac_non_crop(grid%nproma,grid%nblks))
      ALLOCATE(sum_non_crop (grid%nproma,grid%nblks))
      !$ACC ENTER DATA PCREATE(frac_non_crop,sum_non_crop)
      !$ACC ENTER DATA ATTACH(frac_non_crop,sum_non_crop)
      !$ACC PARALLEL LOOP DEFAULT(PRESENT) ASYNC(1) GANG VECTOR
      DO ic = 1, nc
        frac_non_crop(ic,:) = 0._lp
        sum_non_crop (ic,:) = 0._lp
      ENDDO
      !$ACC END PARALLEL LOOP
      !$ACC WAIT(1)
    ENDIF

    ! Allocate local work arrays
    !
    ALLOCATE(OxicDecompCarbon      (nc))
    ALLOCATE(AnoxDecompCarbon      (nc))
    ALLOCATE(water_table           (nc))
    ALLOCATE(frac_tile             (nc))
    ALLOCATE(frac_tile_abs         (nc))
    ALLOCATE(frac1                 (nc))
    ALLOCATE(tiller_area_sl        (nsoil))
    ALLOCATE(CH4_dry_sl            (nsoil))
    ALLOCATE(CO2_dry_sl            (nsoil))
    ALLOCATE( O2_dry_sl            (nsoil))
    ALLOCATE(CH4_wet_sl            (nsoil))
    ALLOCATE(CO2_wet_sl            (nsoil))
    ALLOCATE( O2_wet_sl            (nsoil))
    ALLOCATE(CH4_plant_sl          (nsoil))
    ALLOCATE(CO2_plant_sl          (nsoil))
    ALLOCATE( O2_plant_sl          (nsoil))
    ALLOCATE(CH4_prod_sl           (nsoil))
    ALLOCATE(CO2_prod_sl           (nsoil))
    ALLOCATE( O2_prod_sl           (nsoil))
    ALLOCATE(CH4_oxid_sl           (nsoil))
    ALLOCATE(CO2_oxid_sl           (nsoil))
    ALLOCATE( O2_oxid_sl           (nsoil))
    ALLOCATE(CH4_ebul_sl           (nsoil))
    ALLOCATE( O2_other_sl          (nsoil))
!    ALLOCATE(CH4_excess_sl         (nsoil))
!    ALLOCATE(CO2_excess_sl         (nsoil))
!    ALLOCATE( O2_excess_sl         (nsoil))
    ALLOCATE(relative_moist_sl     (nsoil))
    ALLOCATE(relative_ice_sl       (nsoil))
    ALLOCATE(relative_porevol_sl   (nsoil))
    ALLOCATE(relative_pore_air_sl  (nsoil))
    ALLOCATE(relative_pore_moist_sl(nsoil))
    ALLOCATE(dens2conc_sl          (nsoil))
    ALLOCATE(conc2dens_sl          (nsoil))
    ALLOCATE(D_CH4                 (nsoil))
    ALLOCATE(D_CO2                 (nsoil))
    ALLOCATE(D_O2                  (nsoil))
    ALLOCATE(D_used                (nsoil))
    ALLOCATE(CH4_act               (nsoil))
    ALLOCATE(CO2_act               (nsoil))
    ALLOCATE( O2_act               (nsoil))
    !$ACC DATA CREATE (OxicDecompCarbon,AnoxDecompCarbon,water_table,frac_tile,frac_tile_abs,frac1,tiller_area_sl) &
    !$ACC      CREATE (CH4_dry_sl,CO2_dry_sl,O2_dry_sl,CH4_wet_sl,CO2_wet_sl,O2_wet_sl,CH4_plant_sl,CO2_plant_sl,O2_plant_sl) &
    !$ACC      CREATE (CH4_prod_sl,CO2_prod_sl,O2_prod_sl,CH4_oxid_sl,CO2_oxid_sl,O2_oxid_sl,CH4_ebul_sl,O2_other_sl) &
    !$ACC      CREATE (relative_moist_sl,relative_ice_sl,relative_porevol_sl,relative_pore_air_sl,relative_pore_moist_sl) &
    !$ACC      CREATE (dens2conc_sl,conc2dens_sl,D_CH4,D_CO2,D_O2,D_used,CH4_act,CO2_act,O2_act) ASYNC(1)

    ! Get reference to variables for current block
    !
    ! A2L_ variables
    dsl4jsb_Get_var2D_onChunk(A2L_, CO2_air)
    dsl4jsb_Get_var2D_onChunk(A2L_, press_srf)
    ! SSE_ variables
    dsl4jsb_Get_var3D_onChunk(SSE_, t_soil_sl)
    ! HYDRO_ variables
    dsl4jsb_Get_var2D_onChunk(HYDRO_, bclapp)
    dsl4jsb_Get_var2D_onChunk(HYDRO_, vol_porosity)
    dsl4jsb_Get_var2D_onChunk(HYDRO_, soil_depth)
    dsl4jsb_Get_var3D_onChunk(HYDRO_, wtr_soil_sl)
    dsl4jsb_Get_var3D_onChunk(HYDRO_, ice_soil_sl)
    ! TOP_ variables
    dsl4jsb_Get_var2D_onChunk(TOP_, frac_inun)
    dsl4jsb_Get_var2D_onChunk(TOP_, frac_inun_pre)
    dsl4jsb_Get_var2D_onChunk(TOP_, wt_inun)
    dsl4jsb_Get_var2D_onChunk(TOP_, peat_water_table)
    ! dsl4jsb_Get_var2D_onChunk(A2L_, CH4)
    ! CARBON_ variables
    dsl4jsb_Get_var2D_onChunk(CARBON_, soil_respiration_oxic_ta)
    dsl4jsb_Get_var2D_onChunk(CARBON_, soil_respiration_anoxic_ta)
    dsl4jsb_Get_var2D_onChunk(CARBON_, c_green)
    dsl4jsb_Get_var2D_onChunk(CARBON_, acrotelm_height_ag)
    dsl4jsb_Get_var2D_onChunk(CARBON_, acrotelm_height_bg)
    dsl4jsb_Get_var2D_onChunk(CARBON_, co2flux_fire_all_2_atm_ta)
    dsl4jsb_Get_var2D_onChunk(CARBON_, co2flux_herb_2_atm_ta)
    dsl4jsb_Get_var2D_onChunk(CARBON_, GPP_prev_year_ta)
    ! PHENO_ variables
    dsl4jsb_Get_var2D_onChunk(PHENO_, lai)
    dsl4jsb_Get_var2D_onChunk(PHENO_, GDD_prev_year)
    dsl4jsb_Get_var2D_onChunk(PHENO_, GDD_sum)
    ! CH4_ variables
    dsl4jsb_Get_var2D_onChunk(CH4_, CH4_conc)
    dsl4jsb_Get_var2D_onChunk(CH4_, CO2_conc)
    dsl4jsb_Get_var2D_onChunk(CH4_, O2_conc)
    dsl4jsb_Get_var2D_onChunk(CH4_, mxlay_above_rock)
    dsl4jsb_Get_var3D_onChunk(CH4_, root_dist_sl)
    dsl4jsb_Get_var2D_onChunk(CH4_, CH4_out_net)
    dsl4jsb_Get_var2D_onChunk(CH4_, CO2_out_net)
    dsl4jsb_Get_var2D_onChunk(CH4_, CH4_prod_ratio_out_wet)
    dsl4jsb_Get_var2D_onChunk(CH4_, CH4_dry%prod_out)
    dsl4jsb_Get_var2D_onChunk(CH4_, CH4_dry%oxid_out)
    dsl4jsb_Get_var2D_onChunk(CH4_, CH4_dry%diff_out)
    dsl4jsb_Get_var2D_onChunk(CH4_, CH4_dry%plant_out)
    dsl4jsb_Get_var2D_onChunk(CH4_, CH4_dry%in_out)
!    dsl4jsb_Get_var2D_onChunk(CH4_, CH4_dry%excess)
    dsl4jsb_Get_var2D_onChunk(CH4_, CH4_dry%total_out)
    dsl4jsb_Get_var2D_onChunk(CH4_, CH4_dry%total_initial)
    dsl4jsb_Get_var2D_onChunk(CH4_, CH4_dry%total_final)
    dsl4jsb_Get_var3D_onChunk(CH4_, CH4_dry%dens_sl)
    dsl4jsb_Get_var2D_onChunk(CH4_, CH4_wet%prod_out)
    dsl4jsb_Get_var2D_onChunk(CH4_, CH4_wet%oxid_out)
    dsl4jsb_Get_var2D_onChunk(CH4_, CH4_wet%diff_out)
    dsl4jsb_Get_var2D_onChunk(CH4_, CH4_wet%other_out) ! Ebullition
    dsl4jsb_Get_var2D_onChunk(CH4_, CH4_wet%plant_out)
    dsl4jsb_Get_var2D_onChunk(CH4_, CH4_wet%in_out)
!    dsl4jsb_Get_var2D_onChunk(CH4_, CH4_wet%excess)
    dsl4jsb_Get_var2D_onChunk(CH4_, CH4_wet%total_out)
    dsl4jsb_Get_var2D_onChunk(CH4_, CH4_wet%total_initial)
    dsl4jsb_Get_var2D_onChunk(CH4_, CH4_wet%total_final)
    dsl4jsb_Get_var3D_onChunk(CH4_, CH4_wet%dens_sl)
!    dsl4jsb_Get_var2D_onChunk(CH4_, CH4_tot%prod_out)
!    dsl4jsb_Get_var2D_onChunk(CH4_, CH4_tot%oxid_out)
!    dsl4jsb_Get_var2D_onChunk(CH4_, CH4_tot%diff_out)
!    dsl4jsb_Get_var2D_onChunk(CH4_, CH4_tot%other_out) ! Ebullition
!    dsl4jsb_Get_var2D_onChunk(CH4_, CH4_tot%plant_out)
!    dsl4jsb_Get_var2D_onChunk(CH4_, CH4_tot%total_out)
!    dsl4jsb_Get_var2D_onChunk(CH4_, CH4_tot%total_initial)
!    dsl4jsb_Get_var2D_onChunk(CH4_, CH4_tot%total_final)
!    dsl4jsb_Get_var3D_onChunk(CH4_, CH4_tot%dens_sl)
    dsl4jsb_Get_var2D_onChunk(CH4_, CO2_dry%prod_out)
    dsl4jsb_Get_var2D_onChunk(CH4_, CO2_dry%oxid_out)
    dsl4jsb_Get_var2D_onChunk(CH4_, CO2_dry%diff_out)
    dsl4jsb_Get_var2D_onChunk(CH4_, CO2_dry%plant_out)
    dsl4jsb_Get_var2D_onChunk(CH4_, CO2_dry%in_out)
!    dsl4jsb_Get_var2D_onChunk(CH4_, CO2_dry%excess)
    dsl4jsb_Get_var2D_onChunk(CH4_, CO2_dry%total_out)
    dsl4jsb_Get_var2D_onChunk(CH4_, CO2_dry%total_initial)
    dsl4jsb_Get_var2D_onChunk(CH4_, CO2_dry%total_final)
    dsl4jsb_Get_var3D_onChunk(CH4_, CO2_dry%dens_sl)
    dsl4jsb_Get_var2D_onChunk(CH4_, CO2_wet%prod_out)
    dsl4jsb_Get_var2D_onChunk(CH4_, CO2_wet%oxid_out)
    dsl4jsb_Get_var2D_onChunk(CH4_, CO2_wet%diff_out)
    dsl4jsb_Get_var2D_onChunk(CH4_, CO2_wet%plant_out)
    dsl4jsb_Get_var2D_onChunk(CH4_, CO2_wet%in_out)
!    dsl4jsb_Get_var2D_onChunk(CH4_, CO2_wet%excess)
    dsl4jsb_Get_var2D_onChunk(CH4_, CO2_wet%total_out)
    dsl4jsb_Get_var2D_onChunk(CH4_, CO2_wet%total_initial)
    dsl4jsb_Get_var2D_onChunk(CH4_, CO2_wet%total_final)
    dsl4jsb_Get_var3D_onChunk(CH4_, CO2_wet%dens_sl)
    dsl4jsb_Get_var2D_onChunk(CH4_, O2_dry%prod_out)
    dsl4jsb_Get_var2D_onChunk(CH4_, O2_dry%oxid_out)
    dsl4jsb_Get_var2D_onChunk(CH4_, O2_dry%diff_out)
    dsl4jsb_Get_var2D_onChunk(CH4_, O2_dry%other_out)
    dsl4jsb_Get_var2D_onChunk(CH4_, O2_dry%plant_out)
    dsl4jsb_Get_var2D_onChunk(CH4_, O2_dry%in_out)
!    dsl4jsb_Get_var2D_onChunk(CH4_, O2_dry%excess)
    dsl4jsb_Get_var2D_onChunk(CH4_, O2_dry%total_out)
    dsl4jsb_Get_var2D_onChunk(CH4_, O2_dry%total_initial)
    dsl4jsb_Get_var2D_onChunk(CH4_, O2_dry%total_final)
    dsl4jsb_Get_var3D_onChunk(CH4_, O2_dry%dens_sl)
    dsl4jsb_Get_var2D_onChunk(CH4_, O2_wet%prod_out)
    dsl4jsb_Get_var2D_onChunk(CH4_, O2_wet%oxid_out)
    dsl4jsb_Get_var2D_onChunk(CH4_, O2_wet%diff_out)
    dsl4jsb_Get_var2D_onChunk(CH4_, O2_wet%other_out)
    dsl4jsb_Get_var2D_onChunk(CH4_, O2_wet%plant_out)
    dsl4jsb_Get_var2D_onChunk(CH4_, O2_wet%in_out)
!    dsl4jsb_Get_var2D_onChunk(CH4_, O2_wet%excess)
    dsl4jsb_Get_var2D_onChunk(CH4_, O2_wet%total_out)
    dsl4jsb_Get_var2D_onChunk(CH4_, O2_wet%total_initial)
    dsl4jsb_Get_var2D_onChunk(CH4_, O2_wet%total_final)
    dsl4jsb_Get_var3D_onChunk(CH4_, O2_wet%dens_sl)
    dsl4jsb_Get_var2D_onChunk(CLIMBUF_, min_mmtemp_20yr)
    IF (dsl4jsb_Config(CH4_)%do_fire) THEN
      dsl4jsb_Get_var2D_onChunk(CH4_, CH4_out_fire)
      dsl4jsb_Get_var2D_onChunk(CH4_, fire_em_NO)
      dsl4jsb_Get_var2D_onChunk(CH4_, fire_em_CO)
      dsl4jsb_Get_var2D_onChunk(CH4_, fire_em_iso)
    ENDIF
    IF (dsl4jsb_Config(CH4_)%do_termite) THEN
      dsl4jsb_Get_var2D_onChunk(CH4_, termite_biomass)
      dsl4jsb_Get_var2D_onChunk(CH4_, CH4_out_termite)
    ENDIF
    IF (dsl4jsb_Config(CH4_)%do_herbivore) THEN
      dsl4jsb_Get_var2D_onChunk(CH4_, CH4_out_herbivore)
    ENDIF
    IF (dbg) THEN
      dsl4jsb_Get_var3D_onChunk(CH4_, CH4_wet%prod_out_sl)
      dsl4jsb_Get_var3D_onChunk(CH4_, CH4_wet%oxid_out_sl)
      dsl4jsb_Get_var3D_onChunk(CH4_, CH4_wet%other_out_sl) ! Ebullition
      dsl4jsb_Get_var3D_onChunk(CH4_, CH4_wet%plant_out_sl)
      dsl4jsb_Get_var3D_onChunk(CH4_, CH4_dry%prod_out_sl)
      dsl4jsb_Get_var3D_onChunk(CH4_, CH4_dry%oxid_out_sl)
      dsl4jsb_Get_var3D_onChunk(CH4_, CH4_dry%plant_out_sl)
      dsl4jsb_Get_var3D_onChunk(CH4_, CO2_wet%prod_out_sl)
      dsl4jsb_Get_var3D_onChunk(CH4_, CO2_wet%oxid_out_sl)
      dsl4jsb_Get_var3D_onChunk(CH4_, CO2_wet%plant_out_sl)
      dsl4jsb_Get_var3D_onChunk(CH4_, CO2_dry%prod_out_sl)
      dsl4jsb_Get_var3D_onChunk(CH4_, CO2_dry%oxid_out_sl)
      dsl4jsb_Get_var3D_onChunk(CH4_, CO2_dry%plant_out_sl)
      dsl4jsb_Get_var3D_onChunk(CH4_,  O2_wet%prod_out_sl)
      dsl4jsb_Get_var3D_onChunk(CH4_,  O2_wet%oxid_out_sl)
      dsl4jsb_Get_var3D_onChunk(CH4_,  O2_wet%other_out_sl)
      dsl4jsb_Get_var3D_onChunk(CH4_,  O2_wet%plant_out_sl)
      dsl4jsb_Get_var3D_onChunk(CH4_,  O2_dry%prod_out_sl)
      dsl4jsb_Get_var3D_onChunk(CH4_,  O2_dry%oxid_out_sl)
      dsl4jsb_Get_var3D_onChunk(CH4_,  O2_dry%other_out_sl)
      dsl4jsb_Get_var3D_onChunk(CH4_,  O2_dry%plant_out_sl)
      dsl4jsb_Get_var2D_onChunk(CH4_, CH4_dry%test2d)
      dsl4jsb_Get_var3D_onChunk(CH4_, CH4_dry%test3d)
      dsl4jsb_Get_var2D_onChunk(CH4_, CH4_wet%test2d)
      dsl4jsb_Get_var3D_onChunk(CH4_, CH4_wet%test3d)
      dsl4jsb_Get_var2D_onChunk(CH4_, CO2_dry%test2d)
      dsl4jsb_Get_var3D_onChunk(CH4_, CO2_dry%test3d)
      dsl4jsb_Get_var2D_onChunk(CH4_, CO2_wet%test2d)
      dsl4jsb_Get_var3D_onChunk(CH4_, CO2_wet%test3d)
      dsl4jsb_Get_var2D_onChunk(CH4_,  O2_dry%test2d)
      dsl4jsb_Get_var3D_onChunk(CH4_,  O2_dry%test3d)
      dsl4jsb_Get_var2D_onChunk(CH4_,  O2_wet%test2d)
      dsl4jsb_Get_var3D_onChunk(CH4_,  O2_wet%test3d)
    ENDIF

    IF (dbg) &
      CH4_wet%test3d = CH4_wet%oxid_out_sl

    CALL reset_chem_specie(CH4_wet,nc,nsoil,dbg)
    CALL reset_chem_specie(CO2_wet,nc,nsoil,dbg)
    CALL reset_chem_specie( O2_wet,nc,nsoil,dbg)
    CALL reset_chem_specie(CH4_dry,nc,nsoil,dbg)
    CALL reset_chem_specie(CO2_dry,nc,nsoil,dbg)
    CALL reset_chem_specie( O2_dry,nc,nsoil,dbg)

    IF (lstart .AND. dbg .AND. (iblk==1)) CALL message(routine,'First instance of update_methane')

    ! determine carbon input during time step for oxic and anoxic conditions
    ! convert respiration flux, respiration is defined NEGATIVE
    !$ACC PARALLEL LOOP DEFAULT(PRESENT) ASYNC(1) GANG VECTOR
    DO ic=1,nc
      OxicDecompCarbon(ic) = -soil_respiration_oxic_ta(ic)   * dt
      AnoxDecompCarbon(ic) = -soil_respiration_anoxic_ta(ic) * dt
    ENDDO
    !$ACC END PARALLEL LOOP

    is_peat = .FALSE. ! StW, TODO: Previous no peatland tiles are defined - to be replaced with some lctlib-constant

    IF (is_peat) THEN ! StW: Where to get is_peat?
      !$ACC PARALLEL LOOP DEFAULT(PRESENT) ASYNC(2) GANG VECTOR
      DO ic = 1, nc
        water_table(ic) = peat_water_table(ic)
      ENDDO
      !$ACC END PARALLEL LOOP
    ELSE
      !$ACC PARALLEL LOOP DEFAULT(PRESENT) ASYNC(2) GANG VECTOR
      DO ic = 1, nc
        water_table(ic) = wt_inun(ic)
      ENDDO
      !$ACC END PARALLEL LOOP
    ENDIF
    !$ACC WAIT

    ! Is tile any kind of vegentation? Should always be the case in JSB4 since this process is only active on veg-tiles
    is_veg = dsl4jsb_Lctlib_param(NaturalVegFlag) .OR. &
             dsl4jsb_Lctlib_param(PastureFlag)    .OR. &
             dsl4jsb_Lctlib_param(CropFlag)

    ! In case of dynamic roots, the root distribution should be re-calculated here
    IF (is_veg) THEN

      CALL tile%Get_fraction(ics,ice,iblk,frac_tile)

      IF (is_peat) THEN
        lt = LT_PEAT
      ELSEIF (dsl4jsb_Lctlib_param(PastureFlag) .OR. dsl4jsb_Lctlib_param(CropFlag)) THEN
        lt = LT_AGRICULTURE
      ELSEIF (dsl4jsb_Lctlib_param(GrassFlag)) THEN
        lt = LT_SAVANNAH
      ELSEIF (dsl4jsb_Lctlib_param(PhenologyType) == 3) THEN
        lt = LT_TROPICAL_FOREST
      ELSE
        lt = LT_EXTRATROPICAL_FOREST
      ENDIF

      ! Debug information
      lissue = .FALSE.

      ! TODO: Is it possible to parallelize such a loop on GPU?
      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR &
      !$ACC   PRIVATE( is_veg, is_wet, is_peat, compute_flag, l_stop, l_GPP_dep, lissue, lissue_tst, dbg, lmsg) &
      !$ACC   PRIVATE( tspd, iblk, ics, ice, nc, nsoil, lt, ic, is, iy, iWetness, iWetness_st, iWetness_en)     &
      !$ACC   PRIVATE( acro_height, frac_inun_loc, frac, D_CH4_surface, D_CO2_surface, D_O2_surface)            &
      !$ACC   PRIVATE( CH4_prod_ratio, CH4_diff, CO2_diff, O2_diff)                &
      !$ACC   PRIVATE( mxlay_above_wt, mxlay_acrotelm, nly)
      DO ic = 1, nc

        is_wet = .TRUE. ! StW, TODO: This is a workaround for is_wet not yet properly defined
        nly = INT(mxlay_above_rock(ic))
        CALL plant_setup(is_wet, dsl4jsb_Lctlib_param(woody_PFT), &
                                 dsl4jsb_Lctlib_param(MaxLAI),    &
                         nly, lai(ic), c_green(ic),               &
                         root_dist_sl(ic,:), tiller_area_sl(:))

        compute_flag = .TRUE.
        IF (nly < dsl4jsb_Config(CH4_)%minimum_no_layers) THEN
          compute_flag = .FALSE.
          CH4_out_net(ic) = 0._lp
          CO2_out_net(ic) = -(OxicDecompCarbon(ic) + AnoxDecompCarbon(ic)) / dt
        ENDIF

        IF (frac_tile(ic) < dsl4jsb_Config(CH4_)%tile_threshold) THEN
          compute_flag = .FALSE.
          CH4_out_net(ic) = 0._lp
          CO2_out_net(ic) = -(OxicDecompCarbon(ic) + AnoxDecompCarbon(ic)) / dt

          ! flush stored C out and reset to ensure conservation
          stored_C_sum = 0._lp
          IF (frac_inun_pre(ic) < 1._lp) THEN
            stored_C_sum = stored_C_sum + &
                           (SUM(CO2_dry%dens_sl(ic,1:nly)) +   &
                            SUM(CH4_dry%dens_sl(ic,1:nly))) *  &
                            (1._lp - frac_inun_pre(ic))
          END IF
          IF (frac_inun_pre(ic) > 0._lp) THEN
            stored_C_sum = stored_C_sum + &
                           (SUM(CO2_wet%dens_sl(ic,1:nly)) +   &
                            SUM(CH4_wet%dens_sl(ic,1:nly))) *  &
                            frac_inun_pre(ic)
          END IF

          CO2_out_net(ic) = CO2_out_net(ic) + stored_C_sum / dt

          ! TODO: Some smart acc-directives here?
          DO is=1,nly
            CO2_dry%dens_sl(ic,1:nly) = 0._lp
            CO2_wet%dens_sl(ic,1:nly) = 0._lp
            CH4_dry%dens_sl(ic,1:nly) = 0._lp
            CH4_wet%dens_sl(ic,1:nly) = 0._lp
            O2_dry%dens_sl (ic,1:nly) = 0._lp
            O2_wet%dens_sl (ic,1:nly) = 0._lp
          ENDDO
        ENDIF

        IF (compute_flag) THEN

          frac_inun_loc = frac_inun(ic)
          CALL init_gas_pools(is_peat, nly,                                                         &
                              dsl4jsb_Config(CH4_)%frac_inun_min, frac_inun_loc, frac_inun_pre(ic), &
                              CH4_dry%dens_sl(ic,:), CO2_dry%dens_sl(ic,:), O2_dry%dens_sl(ic,:),   &
                              CH4_wet%dens_sl(ic,:), CO2_wet%dens_sl(ic,:), O2_wet%dens_sl(ic,:),   &
                              CH4_dry_sl     (   :), CO2_dry_sl     (   :), O2_dry_sl     (   :),   &
                              CH4_wet_sl     (   :), CO2_wet_sl     (   :), O2_wet_sl     (   :))

          ! setup for wetness loop:
          ! For mineral soil ignore wet fraction, if inundated fraction too small (inun_frac_loc = 0)
          ! for organic soil or completely flooded grid cell ignore dry fraction
          iWetness_en = 1; IF (frac_inun_loc > 0._lp) iWetness_en = 2
          iWetness_st = 2; IF (frac_inun_loc < 1._lp) iWetness_st = 1

          ! wetness loop: iWetness = 1 -> dry part of tile, iWetness = 2 -> wet part
          DO iWetness = iWetness_st, iWetness_en

            is_Wet = iWetness == 2

            DO is = 1, nsoil
              O2_other_sl  (is) = 0._lp
              CH4_ebul_sl  (is) = 0._lp
              CH4_prod_sl  (is) = 0._lp
              CH4_plant_sl (is) = 0._lp
              CH4_oxid_sl  (is) = 0._lp
              CO2_prod_sl  (is) = 0._lp
              CO2_plant_sl (is) = 0._lp
              CO2_oxid_sl  (is) = 0._lp
              O2_prod_sl   (is) = 0._lp
              O2_plant_sl  (is) = 0._lp
              O2_oxid_sl   (is) = 0._lp
!              CH4_excess_sl(is) = 0._lp
!              CO2_excess_sl(is) = 0._lp
!              O2_excess_sl (is) = 0._lp
            ENDDO

            ! Select correct state variables and save initial state (StW: latter only for debug?)
            IF (is_Wet) THEN
              DO is = 1, nsoil
                CH4_act(is)   = CH4_wet_sl(is)
                CO2_act(is)   = CO2_wet_sl(is)
                 O2_act(is)   =  O2_wet_sl(is)
              ENDDO
              CH4_wet%total_initial(ic) = SUM(CH4_wet%dens_sl(ic,1:nly))
              CO2_wet%total_initial(ic) = SUM(CO2_wet%dens_sl(ic,1:nly))
               O2_wet%total_initial(ic) = SUM( O2_wet%dens_sl(ic,1:nly))
            ELSE
              DO is = 1, nsoil
                CH4_act(is)   = CH4_dry_sl(is)
                CO2_act(is)   = CO2_dry_sl(is)
                 O2_act(is)   =  O2_dry_sl(is)
              ENDDO
              CH4_dry%total_initial(ic) = SUM(CH4_dry%dens_sl(ic,1:nly))
              CO2_dry%total_initial(ic) = SUM(CO2_dry%dens_sl(ic,1:nly))
               O2_dry%total_initial(ic) = SUM( O2_dry%dens_sl(ic,1:nly))
            ENDIF

            ! set up soil column

            ! set up acrotelm height (= 0 for mineral soil)
            acro_height = 0._lp
            IF (is_peat) THEN
              acro_height = acrotelm_height_ag(ic)+acrotelm_height_bg(ic)
            END IF

            CALL soil_column_setup(dsl4jsb_Config(SSE_)%l_freeze, is_peat, is_Wet, nly,                &
                                   soil_w%dz(:), (soil_w%lbounds(:)+soil_w%ubounds(:))*0.5_lp,         &
                                   soil_w%ubounds(:), soil_depth(ic), vol_porosity(ic),                &
                                   water_table(ic), acro_height, wtr_soil_sl(ic,:), ice_soil_sl(ic,:), &
                                   relative_moist_sl(:), relative_ice_sl(:), relative_porevol_sl(:),   &
                                   relative_pore_air_sl(:), relative_pore_moist_sl(:),                 &
                                   dens2conc_sl(:), conc2dens_sl(:), mxlay_above_wt, mxlay_acrotelm)

            ! Calculate the diffusivities for each gas and for water
            ! and air once and then allocate it as necessary
            CALL diffusion_setup(is_peat, nly, mxlay_above_wt, bclapp(ic),                         &
                                 relative_porevol_sl(:), relative_pore_air_sl(:), t_soil_sl(ic,:), &
                                 D_CH4, D_CO2, D_O2, D_CH4_surface, D_CO2_surface, D_O2_surface)

            ! CH4 / CO2 production
            CALL production(is_peat, is_Wet, nly, mxlay_above_wt,       &
                            dsl4jsb_Config(CH4_)%frac_CH4_anox_mineral, &
                            dsl4jsb_Config(CH4_)%frac_CH4_anox_organic, &
                            dsl4jsb_Config(CH4_)%meth_q10,              &
                            root_dist_sl(ic,:),                         &
                            frac_inun_loc, t_soil_sl(ic,:),             &
                            OxicDecompCarbon(ic), AnoxDecompCarbon(ic), &
                            CH4_prod_sl(:), CO2_prod_sl(:),             &
                            O2_prod_sl(:), CH4_prod_ratio)

            ! -----------------------------------------------------------------------------------------------------------
            !              GAS TRANSPORT
            ! -----------------------------------------------------------------------------------------------------------
            ! PLANT TRANSPORT
            ! PLANT TRANSPORT OF CH4
            CALL plant_transport(CH4_consts, soil_w%dz(:), soil_w%ubounds(:), nly, mxlay_above_wt,   &
                                 dt, t_soil_sl(ic,:), relative_porevol_sl(:),                        &
                                 dens2conc_sl(:), conc2dens_sl(:), tiller_area_sl(:), D_CH4_surface, &
                                 CH4_conc(ic), CH4_act(:), CH4_Plant_sl(:)) !, CH4_excess_sl(:))

            ! PLANT TRANSPORT OF CO2
            CALL plant_transport(CO2_consts, soil_w%dz(:), soil_w%ubounds(:), nly, mxlay_above_wt,   &
                                 dt, t_soil_sl(ic,:), relative_porevol_sl(:),                        &
                                 dens2conc_sl(:), conc2dens_sl(:), tiller_area_sl(:), D_CO2_surface, &
                                 CO2_conc(ic), CO2_act(:), CO2_Plant_sl(:)) !, CO2_excess_sl(:))

            ! PLANT TRANSPORT OF O2
            CALL plant_transport(O2_consts, soil_w%dz(:), soil_w%ubounds(:), nly, mxlay_above_wt,    &
                                 dt, t_soil_sl(ic,:), relative_porevol_sl(:),                        &
                                 dens2conc_sl(:), conc2dens_sl(:), tiller_area_sl(:), D_O2_surface,  &
                                 O2_conc(ic), O2_act(:), O2_Plant_sl(:)) !, O2_excess_sl(:))

            ! EBULLITION OF CH4

            ! Ebullition only happens in wet part of cell (in CH4_wet field "other" is actually "ebullition")
            IF (is_Wet) THEN
              CALL ebullition(is_peat, nly, mxlay_above_wt, t_soil_sl(ic,:), press_srf(ic), &
                              dsl4jsb_memory(CH4_)%layer_midpoint_dz(:),                    &
                              relative_porevol_sl(:), water_table(ic),                      &
                              dens2conc_sl(:), conc2dens_sl(:), CH4_act(:), CH4_ebul_sl(:))
            END IF

            ! -------------------------------------------------------------------------------------------
            ! OXIDATION OF CH4
            ! TODO: Pass relevant constants for only wet or dry conditions. No distinguishing between cases in subroutine
            DO is = 1, nly
              CALL oxidation_mm(is_Wet,  dt,                                       &
                                dsl4jsb_Config(CH4_)%scale_oxid_dry,               &
                                dsl4jsb_Config(CH4_)%scale_oxid_wet,               &
                                dsl4jsb_Config(CH4_)%oxid_q10,                     &
                                t_soil_sl(ic,is), relative_moist_sl          (is), &
                                dens2conc_sl(is), conc2dens_sl               (is), &
                                CH4_act     (is), CO2_act    (is), O2_act    (is), &
                                O2_prod_sl  (is),                                  &
                                CH4_oxid_sl (is), CO2_oxid_sl(is), O2_oxid_sl(is), &
                                O2_other_sl (is))
            ENDDO

            ! -------------------------------------------------------------------------------------------
            ! REACTION-DIFFUSION OF ALL GASES
            ! Also: competition for CH4/O2
! TODO: Emergency "solution" to allow compilation for GPU
#ifndef _OPENACC
            CALL reaction_diffusion(dsl4jsb_Config(CH4_)%nstep_diffusion,             &
                                    dsl4jsb_Config(CH4_)%implicit_contrib,            &
                                    dt, nly, mxlay_above_wt,                          &
                                    soil_w%dz(:), t_soil_sl(ic,:),                    &
                                    relative_porevol_sl(:),                           &
                                    dens2conc_sl(:), conc2dens_sl(:),                 &
                                    CH4_prod_sl (:), CO2_prod_sl (:),                 &
                                    CH4_plant_sl(:), CO2_plant_sl(:), O2_plant_sl(:), &
                                    D_CH4,           D_CO2,           D_O2,           &
                                    D_CH4_surface,   D_CO2_surface,   D_O2_surface,   &
                                    CH4_conc   (ic), CO2_conc   (ic), O2_conc   (ic), &
                                    CH4_act     (:), CO2_act     (:), O2_act     (:), &
                                    CH4_oxid_sl (:), CO2_oxid_sl (:), O2_oxid_sl (:), &
                                    CH4_ebul_sl (:),                  O2_other_sl(:), &
                                    CH4_diff,        CO2_diff,        O2_diff,        &
                                    D_used, lissue_tst)
            lissue = lissue .OR. lissue_tst
#endif

            ! Fill global memory for output (fluxes in s-1)
            IF (is_Wet) THEN
              IF (dbg) THEN
                CH4_wet%test2d(ic) = REAL(mxlay_above_wt,lp)
                DO is = 1, nsoil
!                  CH4_wet%test3d(ic,is) = D_CH4(is)
                  CO2_wet%test3d(ic,is) = D_CO2(is)
                   O2_wet%test3d(ic,is) = D_O2(is)
                ENDDO
              ENDIF
              frac = frac_inun_loc
              CALL update_chem_specie(dt,nly,ic,frac,CH4_wet,CH4_diff,CH4_act,CH4_prod_sl,CH4_plant_sl,CH4_oxid_sl, &
                                      dbg,CH4_ebul_sl,1)
              CALL update_chem_specie(dt,nly,ic,frac,CO2_wet,CO2_diff,CO2_act,CO2_prod_sl,CO2_plant_sl,CO2_oxid_sl,dbg)
              CALL update_chem_specie(dt,nly,ic,frac, O2_wet, O2_diff, O2_act, O2_prod_sl, O2_plant_sl, O2_oxid_sl,dbg,O2_other_sl)
              CH4_prod_ratio_out_wet(ic) = CH4_prod_ratio
            ELSE
              IF (dbg) THEN
                CH4_dry%test2d(ic  ) = REAL(mxlay_above_wt,lp)
                DO is = 1, nsoil
                  CH4_dry%test3d(ic,is) = D_CH4(is)
                  CO2_dry%test3d(ic,is) = D_CO2(is)
                   O2_dry%test3d(ic,is) = D_O2(is)
                ENDDO
              ENDIF
              frac = 1._lp-frac_inun_loc
              CALL update_chem_specie(dt,nly,ic,frac,CH4_dry,CH4_diff,CH4_act,CH4_prod_sl,CH4_plant_sl,CH4_oxid_sl,dbg)
              CALL update_chem_specie(dt,nly,ic,frac,CO2_dry,CO2_diff,CO2_act,CO2_prod_sl,CO2_plant_sl,CO2_oxid_sl,dbg)
              CALL update_chem_specie(dt,nly,ic,frac, O2_dry, O2_diff, O2_act, O2_prod_sl, O2_plant_sl, O2_oxid_sl,dbg,O2_other_sl)
            ENDIF

            ! Peatlands updates densities differently
            IF (is_peat) THEN
              frac_inun_pre  (ic)   = 1._lp
              DO is = 1, nsoil
                CH4_dry%dens_sl(ic,is) = 0._lp
                CO2_dry%dens_sl(ic,is) = 0._lp
                 O2_dry%dens_sl(ic,is) = 0._lp
                CH4_wet%dens_sl(ic,is) = CH4_act(is)
                CO2_wet%dens_sl(ic,is) = CO2_act(is)
                 O2_wet%dens_sl(ic,is) =  O2_act(is)
              ENDDO
            ELSE
              frac_inun_pre(ic) = frac_inun_loc
            ENDIF

          ENDDO ! iWetness

        ENDIF ! Compute flag

!        CH4_tot%out          (ic) = CH4_dry%out          (ic) * (1._lp-frac_inun_loc) + CH4_wet%out          (ic) * frac_inun_loc
!        CH4_tot%prod_out     (ic) = CH4_dry%prod_out     (ic) * (1._lp-frac_inun_loc) + CH4_wet%prod_out     (ic) * frac_inun_loc
!        CH4_tot%oxid_out     (ic) = CH4_dry%oxid_out     (ic) * (1._lp-frac_inun_loc) + CH4_wet%oxid_out     (ic) * frac_inun_loc
!        CH4_tot%diff_out     (ic) = CH4_dry%diff_out     (ic) * (1._lp-frac_inun_loc) + CH4_wet%diff_out     (ic) * frac_inun_loc
!        CH4_tot%other_out    (ic) =                                                    CH4_wet%other_out    (ic) * frac_inun_loc
!        CH4_tot%plant_out    (ic) = CH4_dry%plant_out    (ic) * (1._lp-frac_inun_loc) + CH4_wet%plant_out    (ic) * frac_inun_loc
!        CH4_tot%total_out    (ic) = CH4_dry%total_out    (ic) * (1._lp-frac_inun_loc) + CH4_wet%total_out    (ic) * frac_inun_loc
!        CH4_tot%total_initial(ic) = CH4_dry%total_initial(ic) * (1._lp-frac_inun_loc) + CH4_wet%total_initial(ic) * frac_inun_loc
!        CH4_tot%total_final  (ic) = CH4_dry%total_final  (ic) * (1._lp-frac_inun_loc) + CH4_wet%total_final  (ic) * frac_inun_loc

      ENDDO ! ic = landpoints
      !$ACC END PARALLEL LOOP

      ! Other methane processes

      !   Fire
      IF (dsl4jsb_Config(CH4_)%do_fire) THEN
        !$ACC PARALLEL LOOP DEFAULT(PRESENT) ASYNC(1) GANG VECTOR
        DO ic = 1, nc
          CALL calc_methane_fire(co2flux_fire_all_2_atm_ta(ic),lt,CH4_out_fire(ic),fire_em_NO(ic),fire_em_CO(ic),fire_em_iso(ic))
          CH4_out_fire(ic) = -CH4_out_fire(ic) ! Sign convention: +ve downward
        ENDDO
        !$ACC END PARALLEL LOOP
      ENDIF

      !   Termites
      ! TODO, StW: The way the termite CH4 emissions are calculated for crops depends on the processing order of the tiles.
      ! Currently it only gives expected results if crop tiles are processed last.
      ! This is unsatisfactory and should be somehow replaced!>!
      IF ((dsl4jsb_Config(CH4_)%do_termite) .AND. (min_mmtemp_20yr(1) < 1000._lp)) THEN

        l_GPP_dep = dsl4jsb_Lctlib_param(PhenologyType)==3 .OR. dsl4jsb_Lctlib_param(C4flag) ! Crops are separate anyway
        IF (GDD_prev_year(1) < 0._lp) THEN
          GDD => GDD_sum      ! 1st year: Use sum during current year up to now
        ELSE
          GDD => GDD_prev_year ! Any other year: Use sum of previous year
        ENDIF
        !$ACC PARALLEL LOOP DEFAULT(PRESENT) ASYNC(2) GANG VECTOR
        DO ic = 1, nc
          CALL calc_methane_termites(dsl4jsb_Lctlib_param(CropFlag),                                            &
                                     frac_tile(ic), GPP_prev_year_ta(ic), GDD(ic), min_mmtemp_20yr(ic), l_GPP_dep, &
                                     termite_biomass(ic), CH4_out_termite(ic),                                  &
                                     frac_non_crop(ics+ic-1,iblk), sum_non_crop(ics+ic-1,iblk))
        ENDDO
        !$ACC END PARALLEL LOOP
        ! TODO: This assumes that crops occupy only the last tile processed
        IF (.NOT. ASSOCIATED(tile%next_sibling_tile)) THEN
          !$ACC PARALLEL LOOP DEFAULT(PRESENT) ASYNC(2) GANG VECTOR
          DO ic = ics, ice
            frac_non_crop(ic,iblk) = 0._lp
            sum_non_crop (ic,iblk) = 0._lp
          ENDDO
          !$ACC END PARALLEL LOOP
        ENDIF
        !$ACC PARALLEL LOOP DEFAULT(PRESENT) ASYNC(2) GANG VECTOR
        DO ic = 1, nc
          CH4_out_termite(ic) = CH4_out_termite(ic) * dsl4jsb_Config(CH4_)%scale_termite
        ENDDO
        !$ACC END PARALLEL LOOP
      ENDIF

      !   Herbivores
      IF (dsl4jsb_Config(CH4_)%do_herbivore) THEN
        !$ACC PARALLEL LOOP DEFAULT(PRESENT) ASYNC(3) GANG VECTOR
        DO ic = 1, nc
          CALL calc_methane_herbivores(co2flux_herb_2_atm_ta(ic)/CO2_consts%molarMass_kg,dsl4jsb_Lctlib_param(woody_pft), &
                                       CH4_out_herbivore(ic))
          CH4_out_herbivore(ic) = CH4_out_herbivore(ic) * dsl4jsb_Config(CH4_)%scale_herbivore
        ENDDO
        !$ACC END PARALLEL LOOP
      ENDIF
      !$ACC WAIT
    ELSE ! Vegetation tile
      ! Should never happen, since methane process only runs on vegetation tiles
    ENDIF

    IF (lissue) THEN
      CALL message(routine,TRIM(tile%name)//': Tridiagonal matrix solving failed for one or more points',all_print = .TRUE.)
    ENDIF

    IF (dsl4jsb_Config(CH4_)%atm_couple) THEN
      ! Sum total methane emissions for the atmosphere
      dsl4jsb_Get_var2D_onChunk(CH4_, ch4_2_atm)
      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR
      DO ic=1,nc
        ch4_2_atm(ic) = 0._lp
        frac1(ic) = 1._lp
      ENDDO
      !$ACC END PARALLEL LOOP
      labs_tile_calc = .FALSE.
      curr => dsl4jsb_Config(CH4_)%ch4_sources
      DO WHILE (ASSOCIATED(curr))
        IF (curr%isource > 0) THEN
          tile_fac => frac1
          IF (.NOT. ASSOCIATED(curr%dta)) CALL allocate_ch4_source_data(curr) ! Needed by isource in {1,2,5,6,7}
          IF (curr%isource >= 3) THEN
            SELECT CASE (curr%isource)
              CASE (3) ! Wetlands
                !$ACC PARALLEL LOOP DEFAULT(PRESENT) ASYNC(1) GANG VECTOR
                DO ic=1,nc
                  curr%dta(ics+ic-1,iblk) = CH4_wet%total_out(ic) + CH4_dry%total_out(ic)
                ENDDO
                !$ACC END PARALLEL LOOP
              CASE (4) ! Deposition
                !$ACC PARALLEL LOOP DEFAULT(PRESENT) ASYNC(2) GANG VECTOR
                DO ic=1,nc
                  curr%dta(ics+ic-1,iblk) = CH4_wet%in_out(ic) + CH4_dry%in_out(ic)
                ENDDO
                !$ACC END PARALLEL LOOP
              CASE (5) ! Termites
                curr%dta => dsl4jsb_var2d_onDomain(CH4_,CH4_out_termite)
              CASE (6) ! Herbivores
                curr%dta => dsl4jsb_var2d_onDomain(CH4_,CH4_out_herbivore)
              CASE (7) ! Fires
                curr%dta => dsl4jsb_var2d_onDomain(CH4_,CH4_out_fire)
            END SELECT
          ELSE ! Something is to be read from file (climatology (isource==1) or time series (isource==2))
            IF ((iblk==1) .AND. (lfirstcall .OR. is_newmonth(options%current_datetime,options%dtime)) .AND. &
               .NOT. ASSOCIATED(tile%get_prev_sibling_tile())) THEN
              IF (.NOT. curr%ifile%is_readable) CALL open_ch4_source_file(curr%ifile,curr%name_file,year) ! Open file if needed
              CALL read_new_ch4_source_field(curr) ! Read next field
              ! Prepare scaling for correct aggregation
              IF (.NOT. labs_tile_calc) THEN
                CALL model%Get_top_tile(curr_tile)
                curr_tile => curr_tile%first_child_tile
                DO WHILE (ASSOCIATED(curr_tile))
                  IF (trim(curr_tile%name) == 'lake') EXIT
                  curr_tile => curr_tile%next_sibling_tile
                ENDDO
                IF (ASSOCIATED(curr_tile)) THEN
                  CALL curr_tile%Get_fraction(ics,ice,iblk,frac_tile)
                  !$ACC PARALLEL LOOP DEFAULT(PRESENT) ASYNC(3) GANG VECTOR
                  DO ic=1,nc
                    frac_tile_abs(ic) = 1._lp - frac_tile(ic)
                  ENDDO
                  !$ACC END PARALLEL LOOP
                  !$ACC WAIT(3)
                  labs_tile_calc = .TRUE.
                ENDIF
              ENDIF
              tile_fac => frac_tile_abs
            ENDIF
          ENDIF
          !$ACC WAIT
          ! Update the variable for the atmosphere
          !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR
          DO ic=1,nc
            ch4_2_atm(ic) = ch4_2_atm(ic) + curr%dta(ics+ic-1,iblk) / tile_fac(ic) * curr%fac
          ENDDO
          !$ACC END PARALLEL LOOP
        ENDIF
        curr => curr%next
      ENDDO

      ! Sum of sources of emissions of species which have influence on atmospheric methane lifetime
      dsl4jsb_Get_var2D_onChunk(CH4_, co_2_atm)
      dsl4jsb_Get_var2D_onChunk(CH4_, no_2_atm)
      dsl4jsb_Get_var2D_onChunk(CH4_, iso_2_atm)
      !$ACC PARALLEL LOOP DEFAULT(PRESENT) ASYNC(1) GANG VECTOR
      DO ic = 1, nc
        co_2_atm (ic) = 0_lp
        no_2_atm (ic) = 0_lp
        iso_2_atm(ic) = 0_lp
      ENDDO
      !$ACC END PARALLEL LOOP

      ! Fire emissions
      IF (dsl4jsb_Config(CH4_)%do_fire) THEN
        !$ACC PARALLEL LOOP DEFAULT(PRESENT) ASYNC(1) GANG VECTOR
        DO ic = 1, nc
          co_2_atm(ic)  = fire_em_CO(ic)
          no_2_atm(ic)  = fire_em_NO(ic)
          iso_2_atm(ic) = fire_em_iso(ic)
        ENDDO
        !$ACC END PARALLEL LOOP
      ENDIF

      ! Biogenic emissions (MEGAN should be updated before methane)
      IF (dsl4jsb_Config(MEGAN_)%active) THEN
        dsl4jsb_Get_var2D_onChunk(MEGAN_, co%total_out)
        dsl4jsb_Get_var2D_onChunk(MEGAN_, no%total_out)
        dsl4jsb_Get_var2D_onChunk(MEGAN_, iso%total_out)

        !$ACC PARALLEL LOOP DEFAULT(PRESENT) ASYNC(1) GANG VECTOR
        DO ic = 1, nc
          co_2_atm(ic)  = (co_2_atm(ic)  + co%total_out(ic) *1e-3) / 28.01_lp * 12.01_lp    ![(k)g(CO) m-2 s-1] to [kg(C) m-2 s-1]
          no_2_atm(ic)  = (no_2_atm(ic)  + no%total_out(ic) *1e-3) / 30.01_lp * 14.01_lp    ![(k)g(NO) m-2 s-1] to [kg(N) m-2 s-1]
          iso_2_atm(ic) = (iso_2_atm(ic) + iso%total_out(ic)*1e-3) / 68.12_lp*12.01_lp*5._lp![(k)g(iso) m-2 s-1]-> [kg(C) m-2 s-1]
        ENDDO
        !$ACC END PARALLEL LOOP
      ENDIF
      !$ACC WAIT(1)
      lfirstcall = .FALSE.
    ENDIF

    ! Deallocate local work arrays
    !
    !$ACC END DATA
    DEALLOCATE(OxicDecompCarbon)
    DEALLOCATE(AnoxDecompCarbon)
    DEALLOCATE(water_table)
    DEALLOCATE(tiller_area_sl)
    DEALLOCATE(frac_tile)
    DEALLOCATE(frac_tile_abs)
    DEALLOCATE(frac1)
    DEALLOCATE(CH4_dry_sl)
    DEALLOCATE(CO2_dry_sl)
    DEALLOCATE( O2_dry_sl)
    DEALLOCATE(CH4_wet_sl)
    DEALLOCATE(CO2_wet_sl)
    DEALLOCATE( O2_wet_sl)
    DEALLOCATE(CH4_plant_sl)
    DEALLOCATE(CO2_plant_sl)
    DEALLOCATE( O2_plant_sl)
    DEALLOCATE(CH4_prod_sl)
    DEALLOCATE(CO2_prod_sl)
    DEALLOCATE( O2_prod_sl)
    DEALLOCATE(CH4_oxid_sl)
    DEALLOCATE(CO2_oxid_sl)
    DEALLOCATE( O2_oxid_sl)
!    DEALLOCATE(CH4_excess_sl)
!    DEALLOCATE(CO2_excess_sl)
!    DEALLOCATE( O2_excess_sl)
    DEALLOCATE(CH4_ebul_sl)
    DEALLOCATE( O2_other_sl)
    DEALLOCATE(relative_moist_sl)
    DEALLOCATE(relative_ice_sl)
    DEALLOCATE(relative_porevol_sl)
    DEALLOCATE(relative_pore_air_sl)
    DEALLOCATE(relative_pore_moist_sl)
    DEALLOCATE(dens2conc_sl)
    DEALLOCATE(conc2dens_sl)
    DEALLOCATE(D_CH4)
    DEALLOCATE(D_CO2)
    DEALLOCATE(D_O2)
    DEALLOCATE(D_used)
    DEALLOCATE(CH4_act)
    DEALLOCATE(CO2_act)
    DEALLOCATE( O2_act)

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  CONTAINS

    FUNCTION tile_as_ptr(tile) RESULT(ret)

      CLASS(t_jsb_tile_abstract), POINTER :: ret
      CLASS(t_jsb_tile_abstract), TARGET, INTENT(in) :: tile

      ret => tile

    END FUNCTION tile_as_ptr

    SUBROUTINE allocate_ch4_source_data(src)
      TYPE (t_ch4_source_config), POINTER :: src

      ALLOCATE(src%buf(grid%nproma,grid%nblks,1)) ! Only a function to read multiple available - therefore 1 extra dummy dim.
      src%dta => src%buf(:,:,1)                   ! This is how we want the data
      src%dta_alloc = .TRUE.

    END SUBROUTINE allocate_ch4_source_data

    SUBROUTINE read_new_ch4_source_field(src)
      TYPE (t_ch4_source_config), POINTER :: src

      REAL(lp), POINTER :: tmp(:,:,:)

      CALL message(routine,'Reading '//trim(src%name_source)//' from rec. '//TRIM(i2s(ABS(src%crec)))// &
              ' of '//trim(src%name_var)//''//trim(src%name_file))
      tmp => src%ifile%Read_2d_time(src%name_var,src%buf,src%crec,src%crec)
      src%crec = src%crec + 1
      IF (src%crec > src%ifile%nrec) THEN
        src%crec = 1
        IF (src%isource == 2) THEN
          CALL src%ifile%close() ! New file needed if not a climatology is used
          src%ifile%is_readable = .FALSE.
        ENDIF
      ENDIF

    END SUBROUTINE read_new_ch4_source_field

    SUBROUTINE open_ch4_source_file(ifile,fname,year)
      TYPE(t_multi_file), INTENT(inout) :: ifile
      CHARACTER(len=*),   INTENT(in)    :: fname
      INTEGER,            INTENT(in)    :: year

      CHARACTER(:), ALLOCATABLE :: rname
      CHARACTER(4) :: syear
      INTEGER  :: i, ind, nc_id

      CHARACTER(len=*), PARAMETER :: routine = modname//':open_ch4_source_file'

      rname = TRIM(fname)
      ind = INDEX(fname,'<yyyy>')
      IF (ind > 0) THEN
        WRITE(syear,'(i4)') year
        DO i=1,4
          IF (syear(i:i) == ' ') syear(i:i) = '0'
        ENDDO
        rname = fname(1:ind-1)//syear//TRIM(fname(ind+6:))
      ENDIF
      ! Temporarily disable until it is clear whether or not this will be used
      ! IF (ASSOCIATED(jsb_netcdf_open_input(rname,ifile,model%grid_id))) THEN
      !   IF (my_process_is_stdio()) THEN
      !     ifile%is_readable = ifile%is_open
      !     IF (ifile%is_open) THEN
      !       ifile%nrec = ifile%Get_dimlen('time')
      !     ELSE
      !       CALL finish(routine,'Source file '//TRIM(rname)//' could not be opened')
      !     ENDIF
      !   ENDIF
      ! ENDIF
      CALL p_bcast(ifile%nrec,p_io)
      CALL p_bcast(ifile%is_readable,p_io)

    END SUBROUTINE open_ch4_source_file

  END SUBROUTINE update_methane

  SUBROUTINE aggregate_methane(tile, options)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    TYPE(t_jsb_model), POINTER :: model
    dsl4jsb_Def_config(CH4_)
    dsl4jsb_Def_memory(CH4_)

    LOGICAL :: dbg
    CLASS(t_jsb_aggregator), POINTER :: weighted_by_fract

    CHARACTER(len=*), PARAMETER :: routine = modname//':aggregate_methane'

    INTEGER :: iblk, ics, ice

    iblk = options%iblk
    ics  = options%ics
    ice  = options%ice

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    model => Get_model(tile%owner_model_id)

    dsl4jsb_Get_config(CH4_)
    dsl4jsb_Get_memory(CH4_)

    weighted_by_fract => tile%Get_aggregator("weighted_by_fract")
    dbg = dsl4jsb_Config(CH4_)%dbg_process

!    dsl4jsb_Aggregate_onChunk(CH4_, CH4_tot%out          ,  weighted_by_fract)
!    dsl4jsb_Aggregate_onChunk(CH4_, CH4_tot%prod_out     ,  weighted_by_fract)
!    dsl4jsb_Aggregate_onChunk(CH4_, CH4_tot%oxid_out     ,  weighted_by_fract)
!    dsl4jsb_Aggregate_onChunk(CH4_, CH4_tot%diff_out     ,  weighted_by_fract)
!    dsl4jsb_Aggregate_onChunk(CH4_, CH4_tot%other_out    ,  weighted_by_fract)
!    dsl4jsb_Aggregate_onChunk(CH4_, CH4_tot%plant_out    ,  weighted_by_fract)
!    dsl4jsb_Aggregate_onChunk(CH4_, CH4_tot%total_out    ,  weighted_by_fract)
!    dsl4jsb_Aggregate_onChunk(CH4_, CH4_tot%total_initial,  weighted_by_fract)
!    dsl4jsb_Aggregate_onChunk(CH4_, CH4_tot%total_final  ,  weighted_by_fract)
    ! TODO: Do we need the next three for production runs?
    dsl4jsb_Aggregate_onChunk(CH4_, CH4_conc,               weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CH4_, CO2_conc,               weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CH4_, O2_conc,                weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CH4_, mxlay_above_rock,       weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CH4_, root_dist_sl,           weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CH4_, CH4_out_net,            weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CH4_, CO2_out_net,            weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CH4_, CH4_prod_ratio_out_wet, weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CH4_, CH4_dry%prod_out,       weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CH4_, CH4_dry%oxid_out,       weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CH4_, CH4_dry%diff_out,       weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CH4_, CH4_dry%plant_out,      weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CH4_, CH4_dry%in_out,         weighted_by_fract)
!    dsl4jsb_Aggregate_onChunk(CH4_, CH4_dry%excess,         weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CH4_, CH4_dry%total_out,      weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CH4_, CH4_dry%total_initial,  weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CH4_, CH4_dry%total_final,    weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CH4_, CH4_dry%dens_sl,        weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CH4_, CH4_wet%prod_out,       weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CH4_, CH4_wet%oxid_out,       weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CH4_, CH4_wet%diff_out,       weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CH4_, CH4_wet%in_out,         weighted_by_fract)
!    dsl4jsb_Aggregate_onChunk(CH4_, CH4_wet%excess,         weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CH4_, CH4_wet%other_out,      weighted_by_fract) ! Ebullition
    dsl4jsb_Aggregate_onChunk(CH4_, CH4_wet%plant_out,      weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CH4_, CH4_wet%total_out,      weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CH4_, CH4_wet%total_initial,  weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CH4_, CH4_wet%total_final,    weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CH4_, CH4_wet%dens_sl,        weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CH4_, CO2_dry%prod_out,       weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CH4_, CO2_dry%oxid_out,       weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CH4_, CO2_dry%diff_out,       weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CH4_, CO2_dry%plant_out,      weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CH4_, CO2_dry%in_out,         weighted_by_fract)
!    dsl4jsb_Aggregate_onChunk(CH4_, CO2_dry%excess,         weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CH4_, CO2_dry%total_out,      weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CH4_, CO2_dry%total_initial,  weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CH4_, CO2_dry%total_final,    weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CH4_, CO2_dry%dens_sl,        weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CH4_, CO2_wet%prod_out,       weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CH4_, CO2_wet%oxid_out,       weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CH4_, CO2_wet%diff_out,       weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CH4_, CO2_wet%plant_out,      weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CH4_, CO2_wet%in_out,         weighted_by_fract)
!    dsl4jsb_Aggregate_onChunk(CH4_, CO2_wet%excess,         weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CH4_, CO2_wet%total_out,      weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CH4_, CO2_wet%total_initial,  weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CH4_, CO2_wet%total_final,    weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CH4_, CO2_wet%dens_sl,        weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CH4_, O2_dry%prod_out,        weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CH4_, O2_dry%oxid_out,        weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CH4_, O2_dry%diff_out,        weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CH4_, O2_dry%other_out,       weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CH4_, O2_dry%plant_out,       weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CH4_, O2_dry%in_out,          weighted_by_fract)
!    dsl4jsb_Aggregate_onChunk(CH4_, O2_dry%excess,          weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CH4_, O2_dry%total_out,       weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CH4_, O2_dry%total_initial,   weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CH4_, O2_dry%total_final,     weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CH4_, O2_dry%dens_sl,         weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CH4_, O2_wet%prod_out,        weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CH4_, O2_wet%oxid_out,        weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CH4_, O2_wet%diff_out,        weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CH4_, O2_wet%other_out,       weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CH4_, O2_wet%plant_out,       weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CH4_, O2_wet%in_out,          weighted_by_fract)
!    dsl4jsb_Aggregate_onChunk(CH4_, O2_wet%excess,          weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CH4_, O2_wet%total_out,       weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CH4_, O2_wet%total_initial,   weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CH4_, O2_wet%total_final,     weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CH4_, O2_wet%dens_sl,         weighted_by_fract)
    IF (dsl4jsb_Config(CH4_)%do_fire) THEN
      dsl4jsb_Aggregate_onChunk(CH4_, CH4_out_fire, weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(CH4_, fire_em_NO,   weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(CH4_, fire_em_CO,   weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(CH4_, fire_em_iso,  weighted_by_fract)
    ENDIF
    IF (dsl4jsb_Config(CH4_)%do_termite) THEN
      dsl4jsb_Aggregate_onChunk(CH4_, CH4_out_termite, weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(CH4_, termite_biomass, weighted_by_fract)
    ENDIF
    IF (dsl4jsb_Config(CH4_)%do_herbivore) THEN
      dsl4jsb_Aggregate_onChunk(CH4_, CH4_out_herbivore, weighted_by_fract)
    ENDIF
    IF (dsl4jsb_Config(CH4_)%atm_couple) THEN
      dsl4jsb_Aggregate_onChunk(CH4_, ch4_2_atm, weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(CH4_, co_2_atm,  weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(CH4_, no_2_atm,  weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(CH4_, iso_2_atm, weighted_by_fract)
    ENDIF
    IF (dbg) THEN
      dsl4jsb_Aggregate_onChunk(CH4_,CH4_wet%prod_out_sl,  weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(CH4_,CH4_wet%oxid_out_sl,  weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(CH4_,CH4_wet%other_out_sl, weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(CH4_,CH4_wet%plant_out_sl, weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(CH4_,CO2_wet%prod_out_sl,  weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(CH4_,CO2_wet%oxid_out_sl,  weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(CH4_,CO2_wet%plant_out_sl, weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(CH4_, O2_wet%prod_out_sl,  weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(CH4_, O2_wet%oxid_out_sl,  weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(CH4_, O2_wet%other_out_sl, weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(CH4_, O2_wet%plant_out_sl, weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(CH4_,CH4_dry%prod_out_sl,  weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(CH4_,CH4_dry%oxid_out_sl,  weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(CH4_,CH4_dry%other_out_sl, weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(CH4_,CH4_dry%plant_out_sl, weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(CH4_,CO2_dry%prod_out_sl,  weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(CH4_,CO2_dry%oxid_out_sl,  weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(CH4_,CO2_dry%plant_out_sl, weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(CH4_, O2_dry%prod_out_sl,  weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(CH4_, O2_dry%oxid_out_sl,  weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(CH4_, O2_dry%other_out_sl, weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(CH4_, O2_dry%plant_out_sl, weighted_by_fract)
    ENDIF

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE aggregate_methane

#endif
END MODULE mo_methane_interface
