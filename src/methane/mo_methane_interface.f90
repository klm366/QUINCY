!>
!! @brief Contains the interfaces to methane production and transport
!!
!! @author
!!  Stiig Wilkenskjeld, MPI-M Hamburg
!!
!! @par Revision History
!! Adopted from JSBACH3 code by Thomas Kleinen, 2023-03-17
!!
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/license .
!!

MODULE mo_methane_interface
#ifndef __NO_JSBACH__

  USE mo_kind,                   ONLY: wp, i4
  USE mo_exception,              ONLY: message, finish, message_text
  USE mo_jsb_control,            ONLY: debug_on
  USE mo_jsb_model_class,        ONLY: t_jsb_model
  USE mo_jsb_class,              ONLY: Get_model
  USE mo_jsb_grid_class,         ONLY: t_jsb_grid, t_jsb_vgrid
  USE mo_jsb_grid,               ONLY: Get_vgrid, Get_grid
  USE mo_jsb_tile_class,         ONLY: t_jsb_tile_abstract, t_jsb_aggregator
  USE mo_jsb_process_class,      ONLY: t_jsb_process, PHENO_, CH4_, TOP_, CARBON_, A2L_, HYDRO_, SSE_, NLCC_
  USE mo_jsb_task_class,         ONLY: t_jsb_process_task, t_jsb_task_options
  USE mo_jsb_math_constants,     ONLY: eps8
  USE mo_jsb_physical_constants, ONLY: tmelt

  dsl4jsb_Use_config(CH4_)
  dsl4jsb_Use_config(TOP_)
  dsl4jsb_Use_config(SSE_)
  dsl4jsb_Use_memory(CH4_)
  dsl4jsb_Use_memory(TOP_)
  dsl4jsb_Use_memory(SSE_)
  dsl4jsb_Use_memory(CARBON_)
  dsl4jsb_Use_memory(A2L_)
  dsl4jsb_Use_memory(PHENO_)
  dsl4jsb_Use_memory(HYDRO_)
  dsl4jsb_Use_memory(NLCC_)

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
  REAL(wp), POINTER, SAVE :: frac_non_crop(:,:) => NULL()
  REAL(wp), POINTER, SAVE :: sum_non_crop (:,:) => NULL()

CONTAINS

  ! ================================================================================================================================
  !! Constructors for tasks

  ! -------------------------------------------------------------------------------------------------------
  !> Constructor for methane
  !!
  !! @param[in]     model_id     Model id
  !! @return        return_ptr   Instance of process task "methane"
  !!
  FUNCTION Create_task_methane(model_id) RESULT(return_ptr)

    INTEGER,                   INTENT(in) :: model_id
    CLASS(t_jsb_process_task), POINTER    :: return_ptr

    ALLOCATE(tsk_methane::return_ptr)
    CALL return_ptr%Construct(name='methane_processes', process_id=CH4_, owner_model_id=model_id)

  END FUNCTION Create_task_methane

  !-----------------------------------------------------------------------------------------------------
  !> Register tasks: TOP_ process
  !!
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE Register_methane_tasks(process, model_id)

    USE mo_jsb_model_class,   ONLY: t_jsb_model
    USE mo_jsb_class,         ONLY: Get_model

    CLASS(t_jsb_process), INTENT(inout) :: process
    INTEGER,              INTENT(in)    :: model_id

!    TYPE(t_jsb_model), POINTER :: model

!    model   => Get_model(model_id)

    CALL process%Register_task(tsk_methane(model_id))

  END SUBROUTINE Register_methane_tasks

  SUBROUTINE reset_chem_specie(specie)

    USE mo_ch4_memory_class,  ONLY: t_chem_ptrs

    TYPE (t_chem_ptrs), INTENT(inout) :: specie

    specie%diff_out    = 0._wp
    specie%prod_out    = 0._wp
    specie%plant_out   = 0._wp
    specie%oxid_out    = 0._wp
    specie%total_out   = 0._wp
    specie%total_final = 0._wp
!    IF (ASSOCIATED(specie%other_out)) THEN
!      specie%other_out   = 0._wp
!    ENDIF
! >>> Debug
    specie%prod_out_sl = 0._wp
    specie%plant_out_sl= 0._wp
    specie%oxid_out_sl = 0._wp
! <<<

  END SUBROUTINE reset_chem_specie

  SUBROUTINE update_chem_specie(dt,nly,i,frac,specie,diff,act_sl,prod_sl,plant_sl,oxid_sl,other_sl,to_other_i)

    USE mo_ch4_memory_class, ONLY: t_chem_ptrs

    REAL(wp),           INTENT(in)    :: dt
    INTEGER,            INTENT(in)    :: nly
    INTEGER,            INTENT(in)    :: i
    REAL(wp),           INTENT(in)    :: frac
    TYPE (t_chem_ptrs), INTENT(inout) :: specie
    REAL(wp),           INTENT(in)    :: diff
    REAL(wp),           INTENT(in)    :: act_sl(:)
    REAL(wp),           INTENT(in)    :: prod_sl(:)
    REAL(wp),           INTENT(in)    :: plant_sl(:)
    REAL(wp),           INTENT(in)    :: oxid_sl(:)
    REAL(wp), OPTIONAL, INTENT(in)    :: other_sl(:)
    INTEGER,  OPTIONAL, INTENT(in)    :: to_other_i

    INTEGER  :: j
    REAL(wp) :: mul

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
! >>> Debug
    DO j = 1, nly
      specie%prod_out_sl (i,j) = prod_sl (j) * mul
      specie%plant_out_sl(i,j) = plant_sl(j) * mul
      specie%oxid_out_sl (i,j) = oxid_sl (j) * mul
      IF (PRESENT(other_sl)) THEN
        specie%other_out_sl(i,j) = other_sl(j) * mul
      ENDIF
    ENDDO
! <<<

  END SUBROUTINE update_chem_specie

  !-----------------------------------------------------------------------------------------------------
  !> update_methane
  !!
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE update_methane(tile, options)

    USE mo_jsb_time,          ONLY: is_newday, timesteps_per_day, is_time_experiment_start, is_time_restart, get_date_components
    USE mo_ch4_memory_class,  ONLY: t_chem_ptrs
    USE mo_methane_constants, ONLY: CH4_consts, CO2_consts, O2_consts, emission_factors_for_methane,                   &
    &                               LT_TROPICAL_FOREST, LT_EXTRATROPICAL_FOREST, LT_SAVANNAH, LT_AGRICULTURE, LT_PEAT
    USE mo_methane_util,      ONLY: concentration_mmr_2_mol_m3, determine_Henry_cc
    USE mo_methane_process,   ONLY: plant_setup, calc_methane_termites, calc_methane_herbivores, calc_methane_fire
    USE mo_methane_process,   ONLY: soil_column_setup, diffusion_setup, production, plant_transport, oxidation_mm
    USE mo_methane_process,   ONLY: reaction_diffusion, ebullition, init_gas_pools
    USE mo_mpi,               ONLY: p_sum
use mo_util, only: i2s => int2string, r2s => real2string, l2s => logical2string
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    ! Local variables
    !
    TYPE(t_jsb_model), POINTER :: model
    dsl4jsb_Def_config(CH4_)
    dsl4jsb_Def_config(TOP_)
    dsl4jsb_Def_config(SSE_)
    dsl4jsb_Def_memory(CARBON_)
    dsl4jsb_Def_memory(CH4_)
    dsl4jsb_Def_memory(TOP_)
    dsl4jsb_Def_memory(SSE_)
    dsl4jsb_Def_memory(A2L_)
    dsl4jsb_Def_memory(PHENO_)
    dsl4jsb_Def_memory(HYDRO_)
    dsl4jsb_Def_memory(NLCC_)

    ! Multivariables defined for chemical species
    TYPE (t_chem_ptrs) :: &
      & CH4_dry,          &
      & CH4_wet,          &
!      & CH4_tot,          &
      & CO2_dry,          &
      & CO2_wet,          &
      & O2_dry,           &
      & O2_wet

    ! Pointers to variables in memory
    ! Inputs
    dsl4jsb_Real2D_onChunk :: CO2_air_ptr
    dsl4jsb_Real2D_onChunk :: CH4_air_ptr ! StW: TODO: Presently not available from interface - always using namelist value
    dsl4jsb_Real2D_onChunk :: CO2_air
    dsl4jsb_Real2D_onChunk :: CH4_air
    dsl4jsb_Real2D_onChunk :: O2_air
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
    dsl4jsb_Real2D_onChunk :: cflux_herb_2_atm
    dsl4jsb_Real2D_onChunk :: co2flux_fire_all_2_atm
    dsl4jsb_Real2D_onChunk :: CH4_out_fire
    dsl4jsb_Real2D_onChunk :: CH4_out_termite
    dsl4jsb_Real2D_onChunk :: CH4_out_herbivore
    dsl4jsb_Real2D_onChunk :: fire_em_NO
    dsl4jsb_Real2D_onChunk :: fire_em_CO
    dsl4jsb_Real2D_onChunk :: fire_em_iso
    dsl4jsb_Real2D_onChunk :: termite_biomass
    dsl4jsb_Real2D_onChunk :: GPP_prev_year
    dsl4jsb_Real2D_onChunk :: GDD_prev_year
    dsl4jsb_Real2D_onChunk :: GDD_sum
    dsl4jsb_Real2D_onChunk :: min_mmtemp20
    dsl4jsb_Real2D_onChunk :: layer_wt
    dsl4jsb_Real2D_onChunk :: test

    ! Pointers to allocated non-variables
    !
    REAL(wp), POINTER :: OxicDecompCarbon(:)
    REAL(wp), POINTER :: AnoxDecompCarbon(:)
    REAL(wp), POINTER :: water_table(:)
    REAL(wp), POINTER :: tiller_area_sl(:)
    REAL(wp), POINTER :: frac_tile(:)
    REAL(wp), POINTER :: D_CH4(:), D_CO2(:), D_O2(:), D_used(:)
    REAL(wp), POINTER :: CH4_dry_sl(:), CO2_dry_sl(:), O2_dry_sl(:)
    REAL(wp), POINTER :: CH4_wet_sl(:), CO2_wet_sl(:), O2_wet_sl(:)
    REAL(wp), POINTER :: CH4_plant_sl(:), CO2_plant_sl(:), O2_plant_sl(:)
    REAL(wp), POINTER :: CH4_prod_sl(:), CO2_prod_sl(:), O2_prod_sl(:)
    REAL(wp), POINTER :: CH4_oxid_sl(:), CO2_oxid_sl(:), O2_oxid_sl(:)
!    REAL(wp), POINTER :: CH4_excess_sl(:), CO2_excess_sl(:), O2_excess_sl(:)
    REAL(wp), POINTER :: CH4_ebul_sl(:), O2_other_sl(:)
    REAL(wp), POINTER :: relative_moist_sl(:), relative_ice_sl(:)
    REAL(wp), POINTER :: relative_porevol_sl(:), relative_pore_air_sl(:), relative_pore_moist_sl(:)
    REAL(wp), POINTER :: dens2conc_sl(:), conc2dens_sl(:)
    REAL(wp), POINTER :: CH4_act(:), CO2_act(:), O2_act(:)
    ! Pointers for assignment
    REAL(wp), POINTER :: GDD(:)

    TYPE(t_jsb_grid), POINTER :: grid
    TYPE(t_jsb_vgrid), POINTER :: soil_e, soil_w
    REAL(wp) :: stored_C_sum
    LOGICAL  :: lnew_day, lstart, lrestart
    INTEGER  :: hour, minute, sec
    LOGICAL  :: is_veg, is_wet, is_peat, compute_flag, l_stop, l_GPP_dep
    LOGICAL  :: lissue, lissue_tst
    INTEGER  :: tspd, iblk, ics, ice, nc, nsoil, i, lt
    INTEGER  :: iWetness, iWetness_st, iWetness_en
    REAL(wp) :: dt
    REAL(wp) :: acro_height, frac_inun_loc, frac
    REAL(wp) :: D_CH4_surface, D_CO2_surface, D_O2_surface
    REAL(wp) :: CH4_Prod, CO2_Prod, O2_Prod, CH4_prod_ratio
    REAL(wp) :: CH4_diff, CO2_diff, O2_diff
    INTEGER  :: mxlay_above_wt, mxlay_acrotelm, nly

    CHARACTER(len=*), PARAMETER :: routine = modname//':update_methane'

    IF (.NOT. tile%Is_process_calculated(CH4_)) RETURN


    model => Get_model(tile%owner_model_id)
    dsl4jsb_Get_Config(CH4_)
    dsl4jsb_Get_Config(TOP_)
    dsl4jsb_Get_Config(SSE_)

    dt    = 86400._wp / REAL(dsl4jsb_Config(CH4_)%nstep_per_day,wp)

    CALL get_date_components(options%current_datetime, hour=hour, minute=minute, second=sec)
    sec = sec + minute*60 + hour*3600
    IF (MOD(sec,INT(dt)) /= 0) RETURN ! Not a methane time step

    iblk  = options%iblk
    ics   = options%ics
    ice   = options%ice
    nc    = options%nc
    tspd  = timesteps_per_day(dt)

    lstart = is_time_experiment_start(options%current_datetime)
    lrestart = is_time_restart(options%current_datetime)

    IF ((debug_on() .OR. (dsl4jsb_Config(CH4_)%dbg_process .AND. lstart)) .AND. iblk == 1) &
      CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    model  => Get_model(tile%owner_model_id)
    grid   => Get_grid(model%grid_id)
    soil_e => Get_vgrid('soil_depth_energy')
    soil_w => Get_vgrid('soil_depth_water')
    nsoil  = soil_e%n_levels
    ! TODO: This test should be moved outside the time/blk loop
    l_stop  = .FALSE.
    IF (nsoil /= soil_w%n_levels) l_stop = .TRUE.
    IF (.NOT. l_stop) l_stop = ANY(ABS(soil_e%dz - soil_w%dz) > eps8)
    IF (l_stop) CALL finish(routine,'Sorry, present implementation of methane requires water and energy soil levels to coincide')

    ! Allocate global arrays shared between tiles. NOTICE: TODO: These are NEVER deallocated
    IF (.NOT. ASSOCIATED(frac_non_crop) .AND. dsl4jsb_Config(CH4_)%do_termite) THEN
      ALLOCATE(frac_non_crop(grid%nproma,grid%nblks))
      ALLOCATE(sum_non_crop (grid%nproma,grid%nblks))
      frac_non_crop(:,:) = 0._wp
      sum_non_crop (:,:) = 0._wp
    ENDIF

    ! Allocate local work arrays
    !
    ALLOCATE(OxicDecompCarbon      (nc))
    ALLOCATE(AnoxDecompCarbon      (nc))
    ALLOCATE(water_table           (nc))
    ALLOCATE(frac_tile             (nc))
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
    ! TODO: $acc stuff
    ! Get reference to variables for current block
    !
    dsl4jsb_Get_memory(CH4_)
    dsl4jsb_Get_memory(TOP_)
    dsl4jsb_Get_memory(SSE_)
    dsl4jsb_Get_memory(HYDRO_)
    dsl4jsb_Get_memory(CARBON_)
    dsl4jsb_Get_memory(A2L_)
    dsl4jsb_Get_memory(PHENO_)
    dsl4jsb_Get_memory(NLCC_)

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
    ! A2L_ variables
    dsl4jsb_Get_var2D_onChunk(A2L_, CO2_air)
    dsl4jsb_Get_var2D_onChunk(A2L_, press_srf)
    ! dsl4jsb_Get_var2D_onChunk(A2L_, CH4)
    ! CARBON_ variables
    dsl4jsb_Get_var2D_onChunk(CARBON_, soil_respiration_oxic_ta)
    dsl4jsb_Get_var2D_onChunk(CARBON_, soil_respiration_anoxic_ta)
    dsl4jsb_Get_var2D_onChunk(CARBON_, c_green)
    dsl4jsb_Get_var2D_onChunk(CARBON_, acrotelm_height_ag)
    dsl4jsb_Get_var2D_onChunk(CARBON_, acrotelm_height_bg)
    dsl4jsb_Get_var2D_onChunk(CARBON_, co2flux_fire_all_2_atm)
    dsl4jsb_Get_var2D_onChunk(CARBON_, cflux_herb_2_atm)
    dsl4jsb_Get_var2D_onChunk(CARBON_, GPP_prev_year)
    ! PHENO_ variables
    dsl4jsb_Get_var2D_onChunk(PHENO_, lai)
    dsl4jsb_Get_var2D_onChunk(PHENO_, GDD_prev_year)
    dsl4jsb_Get_var2D_onChunk(PHENO_, GDD_sum)
    ! CH4_ variables
    dsl4jsb_Get_var2D_onChunk(CH4_, CH4_air)
    dsl4jsb_Get_var2D_onChunk(CH4_, CO2_air)
    dsl4jsb_Get_var2D_onChunk(CH4_, O2_air)
    dsl4jsb_Get_var2D_onChunk(CH4_, mxlay_above_rock)
    dsl4jsb_Get_var3D_onChunk(CH4_, root_dist_sl)
    dsl4jsb_Get_var2D_onChunk(CH4_, CH4_out_net)
    dsl4jsb_Get_var2D_onChunk(CH4_, CO2_out_net)
    dsl4jsb_Get_var2D_onChunk(CH4_, CH4_prod_ratio_out_wet)
    dsl4jsb_Get_var2D_onChunk(CH4_, CH4_dry%out)
    dsl4jsb_Get_var2D_onChunk(CH4_, CH4_dry%prod_out)
    dsl4jsb_Get_var2D_onChunk(CH4_, CH4_dry%oxid_out)
    dsl4jsb_Get_var2D_onChunk(CH4_, CH4_dry%diff_out)
    dsl4jsb_Get_var2D_onChunk(CH4_, CH4_dry%plant_out)
!    dsl4jsb_Get_var2D_onChunk(CH4_, CH4_dry%excess)
    dsl4jsb_Get_var2D_onChunk(CH4_, CH4_dry%total_out)
    dsl4jsb_Get_var2D_onChunk(CH4_, CH4_dry%total_initial)
    dsl4jsb_Get_var2D_onChunk(CH4_, CH4_dry%total_final)
    dsl4jsb_Get_var3D_onChunk(CH4_, CH4_dry%dens_sl)
    dsl4jsb_Get_var2D_onChunk(CH4_, CH4_wet%out)
    dsl4jsb_Get_var2D_onChunk(CH4_, CH4_wet%prod_out)
    dsl4jsb_Get_var2D_onChunk(CH4_, CH4_wet%oxid_out)
    dsl4jsb_Get_var2D_onChunk(CH4_, CH4_wet%diff_out)
    dsl4jsb_Get_var2D_onChunk(CH4_, CH4_wet%other_out) ! Ebullition
    dsl4jsb_Get_var2D_onChunk(CH4_, CH4_wet%plant_out)
!    dsl4jsb_Get_var2D_onChunk(CH4_, CH4_wet%excess)
    dsl4jsb_Get_var2D_onChunk(CH4_, CH4_wet%total_out)
    dsl4jsb_Get_var2D_onChunk(CH4_, CH4_wet%total_initial)
    dsl4jsb_Get_var2D_onChunk(CH4_, CH4_wet%total_final)
    dsl4jsb_Get_var3D_onChunk(CH4_, CH4_wet%dens_sl)
!    dsl4jsb_Get_var2D_onChunk(CH4_, CH4_tot%out)
!    dsl4jsb_Get_var2D_onChunk(CH4_, CH4_tot%prod_out)
!    dsl4jsb_Get_var2D_onChunk(CH4_, CH4_tot%oxid_out)
!    dsl4jsb_Get_var2D_onChunk(CH4_, CH4_tot%diff_out)
!    dsl4jsb_Get_var2D_onChunk(CH4_, CH4_tot%other_out) ! Ebullition
!    dsl4jsb_Get_var2D_onChunk(CH4_, CH4_tot%plant_out)
!    dsl4jsb_Get_var2D_onChunk(CH4_, CH4_tot%total_out)
!    dsl4jsb_Get_var2D_onChunk(CH4_, CH4_tot%total_initial)
!    dsl4jsb_Get_var2D_onChunk(CH4_, CH4_tot%total_final)
!    dsl4jsb_Get_var3D_onChunk(CH4_, CH4_tot%dens_sl)
    dsl4jsb_Get_var2D_onChunk(CH4_, CO2_dry%out)
    dsl4jsb_Get_var2D_onChunk(CH4_, CO2_dry%prod_out)
    dsl4jsb_Get_var2D_onChunk(CH4_, CO2_dry%oxid_out)
    dsl4jsb_Get_var2D_onChunk(CH4_, CO2_dry%diff_out)
    dsl4jsb_Get_var2D_onChunk(CH4_, CO2_dry%plant_out)
!    dsl4jsb_Get_var2D_onChunk(CH4_, CO2_dry%excess)
    dsl4jsb_Get_var2D_onChunk(CH4_, CO2_dry%total_out)
    dsl4jsb_Get_var2D_onChunk(CH4_, CO2_dry%total_initial)
    dsl4jsb_Get_var2D_onChunk(CH4_, CO2_dry%total_final)
    dsl4jsb_Get_var3D_onChunk(CH4_, CO2_dry%dens_sl)
    dsl4jsb_Get_var2D_onChunk(CH4_, CO2_wet%out)
    dsl4jsb_Get_var2D_onChunk(CH4_, CO2_wet%prod_out)
    dsl4jsb_Get_var2D_onChunk(CH4_, CO2_wet%oxid_out)
    dsl4jsb_Get_var2D_onChunk(CH4_, CO2_wet%diff_out)
    dsl4jsb_Get_var2D_onChunk(CH4_, CO2_wet%plant_out)
!    dsl4jsb_Get_var2D_onChunk(CH4_, CO2_wet%excess)
    dsl4jsb_Get_var2D_onChunk(CH4_, CO2_wet%total_out)
    dsl4jsb_Get_var2D_onChunk(CH4_, CO2_wet%total_initial)
    dsl4jsb_Get_var2D_onChunk(CH4_, CO2_wet%total_final)
    dsl4jsb_Get_var3D_onChunk(CH4_, CO2_wet%dens_sl)
    dsl4jsb_Get_var2D_onChunk(CH4_, O2_dry%out)
    dsl4jsb_Get_var2D_onChunk(CH4_, O2_dry%prod_out)
    dsl4jsb_Get_var2D_onChunk(CH4_, O2_dry%oxid_out)
    dsl4jsb_Get_var2D_onChunk(CH4_, O2_dry%diff_out)
    dsl4jsb_Get_var2D_onChunk(CH4_, O2_dry%other_out)
    dsl4jsb_Get_var2D_onChunk(CH4_, O2_dry%plant_out)
!    dsl4jsb_Get_var2D_onChunk(CH4_, O2_dry%excess)
    dsl4jsb_Get_var2D_onChunk(CH4_, O2_dry%total_out)
    dsl4jsb_Get_var2D_onChunk(CH4_, O2_dry%total_initial)
    dsl4jsb_Get_var2D_onChunk(CH4_, O2_dry%total_final)
    dsl4jsb_Get_var3D_onChunk(CH4_, O2_dry%dens_sl)
    dsl4jsb_Get_var2D_onChunk(CH4_, O2_wet%out)
    dsl4jsb_Get_var2D_onChunk(CH4_, O2_wet%prod_out)
    dsl4jsb_Get_var2D_onChunk(CH4_, O2_wet%oxid_out)
    dsl4jsb_Get_var2D_onChunk(CH4_, O2_wet%diff_out)
    dsl4jsb_Get_var2D_onChunk(CH4_, O2_wet%other_out)
    dsl4jsb_Get_var2D_onChunk(CH4_, O2_wet%plant_out)
!    dsl4jsb_Get_var2D_onChunk(CH4_, O2_wet%excess)
    dsl4jsb_Get_var2D_onChunk(CH4_, O2_wet%total_out)
    dsl4jsb_Get_var2D_onChunk(CH4_, O2_wet%total_initial)
    dsl4jsb_Get_var2D_onChunk(CH4_, O2_wet%total_final)
    dsl4jsb_Get_var3D_onChunk(CH4_, O2_wet%dens_sl)
    dsl4jsb_Get_var2D_onChunk(CH4_, CH4_out_termite)
    dsl4jsb_Get_var2D_onChunk(CH4_, CH4_out_herbivore)
    dsl4jsb_Get_var2D_onChunk(CH4_, CH4_out_fire)
    dsl4jsb_Get_var2D_onChunk(CH4_, fire_em_NO)
    dsl4jsb_Get_var2D_onChunk(CH4_, fire_em_CO)
    dsl4jsb_Get_var2D_onChunk(CH4_, fire_em_iso)
    dsl4jsb_Get_var2D_onChunk(CH4_, termite_biomass)
    dsl4jsb_Get_var2D_onChunk(CH4_, layer_wt)
    dsl4jsb_Get_var2D_onChunk(CH4_, test)
! >>> Debug
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
! <<<

    CALL reset_chem_specie(CH4_wet)
    CALL reset_chem_specie(CO2_wet)
    CALL reset_chem_specie( O2_wet)
    CALL reset_chem_specie(CH4_dry)
    CALL reset_chem_specie(CO2_dry)
    CALL reset_chem_specie( O2_dry)

    IF (lstart .AND. (dsl4jsb_Config(CH4_)%dbg_process) .AND. (iblk==1)) CALL message(routine,'First instance of update_methane')

    IF (dsl4jsb_Config(CH4_)%GHG_from_atm) THEN
      IF (.FALSE.) THEN ! TODO: Some test for atmospheric CH4 concentration available
        CH4_air_ptr => dsl4jsb_var2D_onChunk(CH4_, CH4_air)
        CH4_air(:) = concentration_mmr_2_mol_m3(CH4_consts,CH4_air_ptr)
      ELSE
        CH4_air(:) = dsl4jsb_Config(CH4_)%CH4_conc ! Equal to below
      ENDIF
      CO2_air_ptr => dsl4jsb_var2D_onChunk(A2L_, CO2_air)
      CO2_air(:) = concentration_mmr_2_mol_m3(CO2_consts,CO2_air_ptr)
    ELSE
      CH4_air(:) = dsl4jsb_Config(CH4_)%CH4_conc
      CO2_air(:) = dsl4jsb_Config(CH4_)%CO2_conc
    ENDIF
    O2_air(:) = dsl4jsb_Config(CH4_)%O2_conc

    ! determine carbon input during time step for oxic and anoxic conditions
    ! convert respiration flux, respiration is defined NEGATIVE
    OxicDecompCarbon = -soil_respiration_oxic_ta   * dt
    AnoxDecompCarbon = -soil_respiration_anoxic_ta * dt

    is_peat = .FALSE. ! StW, TODO: Previous no peatland tiles are defined - to be replaced with some lctlib-constant

    IF (is_peat) THEN ! StW: Where to get is_peat?
      water_table(:) = peat_water_table(:)
    ELSE
      water_table(:) = wt_inun(:)
    ENDIF

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

      DO i=1,nc

        is_wet = .TRUE. ! StW, TODO: This is a workaround for is_wet not yet properly defined
        nly = INT(mxlay_above_rock(i))
        CALL plant_setup(is_wet, dsl4jsb_Lctlib_param(woody_PFT), &
                                 dsl4jsb_Lctlib_param(MaxLAI),    &
                         nly, lai(i), c_green(i),                 &
                         root_dist_sl(i,:), tiller_area_sl(:))

        compute_flag = .TRUE.
        IF (nly < dsl4jsb_Config(CH4_)%minimum_no_layers) THEN
          compute_flag = .FALSE.
          CH4_out_net(i) = 0._wp
          CH4_wet%out(i) = 0._wp
          CH4_dry%out(i) = 0._wp
          CO2_out_net(i) = -(OxicDecompCarbon(i) + AnoxDecompCarbon(i)) / dt
          CO2_wet%out(i) = CO2_out_net(i)
          CO2_dry%out(i) = 0._wp
        ENDIF

        IF (frac_tile(i) < dsl4jsb_Config(CH4_)%tile_threshold) THEN
          compute_flag = .FALSE.
          CH4_out_net(i) = 0._wp
          CH4_wet%out(i) = 0._wp
          CH4_dry%out(i) = 0._wp
          CO2_out_net(i) = -(OxicDecompCarbon(i) + AnoxDecompCarbon(i)) / dt

          ! flush stored C out and reset to ensure conservation
          stored_C_sum = 0._wp
          IF (frac_inun_pre(i) < 1._wp) THEN
            stored_C_sum = stored_C_sum + &
                           (SUM(CO2_dry%dens_sl(i,1:nly)) +   &
                            SUM(CH4_dry%dens_sl(i,1:nly))) *  &
                            (1._wp - frac_inun_pre(i))
          END IF
          IF (frac_inun_pre(i) > 0._wp) THEN
            stored_C_sum = stored_C_sum + &
                           (SUM(CO2_wet%dens_sl(i,1:nly)) +   &
                            SUM(CH4_wet%dens_sl(i,1:nly))) *  &
                            frac_inun_pre(i)
          END IF

          CO2_out_net(i) = CO2_out_net(i) + stored_C_sum / dt
          CO2_dry%out(i) = CO2_out_net(i)
          CO2_wet%out(i) = 0._wp

          CO2_dry%dens_sl(i,1:nly) = 0._wp
          CO2_wet%dens_sl(i,1:nly) = 0._wp
          CH4_dry%dens_sl(i,1:nly) = 0._wp
          CH4_wet%dens_sl(i,1:nly) = 0._wp
          O2_dry%dens_sl (i,1:nly) = 0._wp
          O2_wet%dens_sl (i,1:nly) = 0._wp
        ENDIF

        IF (compute_flag) THEN

          frac_inun_loc = frac_inun(i)
          CALL init_gas_pools(is_peat, nly,                                                        &
                              dsl4jsb_Config(CH4_)%frac_inun_min, frac_inun_loc, frac_inun_pre(i), &
                              CH4_dry%dens_sl(i,:), CO2_dry%dens_sl(i,:), O2_dry%dens_sl(i,:),     &
                              CH4_wet%dens_sl(i,:), CO2_wet%dens_sl(i,:), O2_wet%dens_sl(i,:),     &
                              CH4_dry_sl     (  :), CO2_dry_sl     (  :), O2_dry_sl     (  :),     &
                              CH4_wet_sl     (  :), CO2_wet_sl     (  :), O2_wet_sl     (  :))

          ! setup for wetness loop:
          ! For mineral soil ignore wet fraction, if inundated fraction too small (inun_frac_loc = 0)
          ! for organic soil or completely flooded grid cell ignore dry fraction
          iWetness_en = 1; IF (frac_inun_loc > 0._wp) iWetness_en = 2
          iWetness_st = 2; IF (frac_inun_loc < 1._wp) iWetness_st = 1

          ! wetness loop: iWetness = 1 -> dry part of tile, iWetness = 2 -> wet part
          DO iWetness = iWetness_st, iWetness_en

            is_Wet = iWetness == 2

            O2_other_sl  (:) = 0._wp
            CH4_ebul_sl  (:) = 0._wp
            CH4_prod_sl  (:) = 0._wp
            CH4_plant_sl (:) = 0._wp
            CH4_oxid_sl  (:) = 0._wp
            CO2_prod_sl  (:) = 0._wp
            CO2_plant_sl (:) = 0._wp
            CO2_oxid_sl  (:) = 0._wp
            O2_prod_sl   (:) = 0._wp
            O2_plant_sl  (:) = 0._wp
            O2_oxid_sl   (:) = 0._wp
!            CH4_excess_sl(:) = 0._wp
!            CO2_excess_sl(:) = 0._wp
!            O2_excess_sl (:) = 0._wp

            ! Select correct state variables and save initial state (StW: latter only for debug?)
            IF (is_Wet) THEN
              CH4_act   = CH4_wet_sl
              CO2_act   = CO2_wet_sl
               O2_act   =  O2_wet_sl
              CH4_wet%total_initial(i) = SUM(CH4_wet%dens_sl(i,1:nly))
              CO2_wet%total_initial(i) = SUM(CO2_wet%dens_sl(i,1:nly))
               O2_wet%total_initial(i) = SUM( O2_wet%dens_sl(i,1:nly))
            ELSE
              CH4_act   = CH4_dry_sl
              CO2_act   = CO2_dry_sl
               O2_act   =  O2_dry_sl
              CH4_dry%total_initial(i) = SUM(CH4_dry%dens_sl(i,1:nly))
              CO2_dry%total_initial(i) = SUM(CO2_dry%dens_sl(i,1:nly))
               O2_dry%total_initial(i) = SUM( O2_dry%dens_sl(i,1:nly))
            ENDIF

            ! set up soil column 

            ! set up acrotelm height (= 0 for mineral soil)
            acro_height = 0._wp
            IF (is_peat) THEN
              acro_height = acrotelm_height_ag(i)+acrotelm_height_bg(i)
            END IF

            CALL soil_column_setup(dsl4jsb_Config(SSE_)%l_freeze, is_peat, is_Wet, nly,              &
                                   soil_w%dz(:), (soil_w%lbounds(:)+soil_w%ubounds(:))*0.5_wp,       &
                                   soil_w%ubounds(:), soil_depth(i), vol_porosity(i),                &
                                   water_table(i), acro_height, wtr_soil_sl(i,:), ice_soil_sl(i,:),  &
                                   relative_moist_sl(:), relative_ice_sl(:), relative_porevol_sl(:), &
                                   relative_pore_air_sl(:), relative_pore_moist_sl(:),               &
                                   dens2conc_sl(:), conc2dens_sl(:), mxlay_above_wt, mxlay_acrotelm)

            ! Calculate the diffusivities for each gas and for water
            ! and air once and then allocate it as necessary
            CALL diffusion_setup(is_peat, nly, mxlay_above_wt, bclapp(i),                         &
                                 relative_porevol_sl(:), relative_pore_air_sl(:), t_soil_sl(i,:), &
                                 D_CH4, D_CO2, D_O2, D_CH4_surface, D_CO2_surface, D_O2_surface)

            ! CH4 / CO2 production
            CALL production(is_peat, is_Wet, nly, mxlay_above_wt,       &
                            dsl4jsb_Config(CH4_)%frac_CH4_anox_mineral, &
                            dsl4jsb_Config(CH4_)%frac_CH4_anox_organic, &
                            dsl4jsb_Config(CH4_)%meth_q10,              &
                            root_dist_sl(i,:),                          &
                            frac_inun_loc, t_soil_sl(i,:),              &
                            OxicDecompCarbon(i), AnoxDecompCarbon(i),   &
                            CH4_prod_sl(:), CO2_prod_sl(:),             &
                            O2_prod_sl(:), CH4_prod_ratio)

            ! -----------------------------------------------------------------------------------------------------------
            !              GAS TRANSPORT
            ! -----------------------------------------------------------------------------------------------------------
            ! PLANT TRANSPORT
            ! PLANT TRANSPORT OF CH4
            CALL plant_transport(CH4_consts, soil_w%dz(:), soil_w%ubounds(:), nly, mxlay_above_wt,   &
                                 dt, t_soil_sl(i,:), relative_porevol_sl(:),                         &
                                 dens2conc_sl(:), conc2dens_sl(:), tiller_area_sl(:), D_CH4_surface, &
                                 CH4_air(i), CH4_act(:), CH4_Plant_sl(:)) !, CH4_excess_sl(:))

            ! PLANT TRANSPORT OF CO2
            CALL plant_transport(CO2_consts, soil_w%dz(:), soil_w%ubounds(:), nly, mxlay_above_wt,   &
                                 dt, t_soil_sl(i,:), relative_porevol_sl(:),                         &
                                 dens2conc_sl(:), conc2dens_sl(:), tiller_area_sl(:), D_CO2_surface, &
                                 CO2_air(i), CO2_act(:), CO2_Plant_sl(:)) !, CO2_excess_sl(:))

            ! PLANT TRANSPORT OF O2
            CALL plant_transport(O2_consts, soil_w%dz(:), soil_w%ubounds(:), nly, mxlay_above_wt,    &
                                 dt, t_soil_sl(i,:), relative_porevol_sl(:),                         &
                                 dens2conc_sl(:), conc2dens_sl(:), tiller_area_sl(:), D_O2_surface,  &
                                 O2_air(i), O2_act(:), O2_Plant_sl(:)) !, O2_excess_sl(:))

            ! EBULLITION OF CH4

            ! Ebullition only happens in wet part of cell (in CH4_wet field "other" is actually "ebullition")
            IF (is_Wet) THEN
              CALL ebullition(is_peat, nly, mxlay_above_wt, t_soil_sl(i,:), press_srf(i),   &
                              dsl4jsb_memory(CH4_)%layer_midpoint_dz(:),                    &
                              relative_porevol_sl(:), water_table(i),                       &
                              dens2conc_sl(:), conc2dens_sl(:), CH4_act(:), CH4_ebul_sl(:))
            ELSE
              CH4_ebul_sl(:) = 0._wp
            END IF

            ! -------------------------------------------------------------------------------------------
            ! OXIDATION OF CH4
            ! TODO: Pass relevant constants for only wet or dry conditions. No distinguishing between cases in subroutine
            CALL oxidation_mm(is_Wet,  dt,                                                &
                              dsl4jsb_Config(CH4_)%scale_oxid_dry,                        &
                              dsl4jsb_Config(CH4_)%scale_oxid_wet,                        &
                              t_soil_sl (i,1:nly), relative_moist_sl             (1:nly), &
                              dens2conc_sl(1:nly), conc2dens_sl                  (1:nly), &
                              CH4_act     (1:nly), CO2_act    (1:nly), O2_act    (1:nly), &
                              O2_prod_sl  (1:nly),                                        &
                              CH4_oxid_sl (1:nly), CO2_oxid_sl(1:nly), O2_oxid_sl(1:nly), &
                              O2_other_sl (1:nly))

            ! -------------------------------------------------------------------------------------------
            ! REACTION-DIFFUSION OF ALL GASES
            ! Also: competition for CH4/O2
            CALL reaction_diffusion(dsl4jsb_Config(CH4_)%nstep_diffusion,             &
                                    dsl4jsb_Config(CH4_)%implicit_contrib,            &
                                    dt, nly, mxlay_above_wt,                          &
                                    soil_w%dz(:), t_soil_sl(i,:),                     &
                                    relative_porevol_sl(:),                           &
                                    dens2conc_sl(:), conc2dens_sl(:),                 &
                                    CH4_prod_sl (:), CO2_prod_sl (:),                 &
                                    CH4_plant_sl(:), CO2_plant_sl(:), O2_plant_sl(:), &
                                    D_CH4,           D_CO2,           D_O2,           &
                                    D_CH4_surface,   D_CO2_surface,   D_O2_surface,   &
                                    CH4_air     (i), CO2_air     (i), O2_air     (i), &
                                    CH4_act     (:), CO2_act     (:), O2_act     (:), &
                                    CH4_oxid_sl (:), CO2_oxid_sl (:), O2_oxid_sl (:), &
                                    CH4_ebul_sl (:),                  O2_other_sl(:), &
                                    CH4_diff,        CO2_diff,        O2_diff,        &
                                    D_used, lissue_tst)
            lissue = lissue .OR. lissue_tst

            ! Fill global memory for output (fluxes in s-1)
            IF (is_Wet) THEN
              frac = frac_inun_loc
              CALL update_chem_specie(dt,nly,i,frac,CH4_wet,CH4_diff,CH4_act,CH4_prod_sl,CH4_plant_sl,CH4_oxid_sl,CH4_ebul_sl,1)
              CALL update_chem_specie(dt,nly,i,frac,CO2_wet,CO2_diff,CO2_act,CO2_prod_sl,CO2_plant_sl,CO2_oxid_sl)
              CALL update_chem_specie(dt,nly,i,frac, O2_wet, O2_diff, O2_act, O2_prod_sl, O2_plant_sl, O2_oxid_sl,O2_other_sl)
              CH4_prod_ratio_out_wet(i) = CH4_prod_ratio
            ELSE
              frac = 1._wp-frac_inun_loc
              CALL update_chem_specie(dt,nly,i,frac,CH4_dry,CH4_diff,CH4_act,CH4_prod_sl,CH4_plant_sl,CH4_oxid_sl)
              CALL update_chem_specie(dt,nly,i,frac,CO2_dry,CO2_diff,CO2_act,CO2_prod_sl,CO2_plant_sl,CO2_oxid_sl)
              CALL update_chem_specie(dt,nly,i,frac, O2_dry, O2_diff, O2_act, O2_prod_sl, O2_plant_sl, O2_oxid_sl,O2_other_sl)
            ENDIF

            ! Peatlands updates densities differently
            IF (is_peat) THEN
              frac_inun_pre  (i)   = 1._wp
              CH4_dry%dens_sl(i,:) = 0._wp
              CO2_dry%dens_sl(i,:) = 0._wp
               O2_dry%dens_sl(i,:) = 0._wp
              CH4_wet%dens_sl(i,:) = CH4_act(:)
              CO2_wet%dens_sl(i,:) = CO2_act(:)
               O2_wet%dens_sl(i,:) =  O2_act(:)
            ELSE
              frac_inun_pre(i) = frac_inun_loc
            ENDIF

          ENDDO ! iWetness

        ENDIF ! Compute flag

        CH4_dry%out(i) = CH4_dry%total_out(i)! * (1._wp - frac_inun_loc)
        CO2_dry%out(i) = CO2_dry%total_out(i)! * (1._wp - frac_inun_loc)
         O2_dry%out(i) =  O2_dry%total_out(i)! * (1._wp - frac_inun_loc)
        CH4_wet%out(i) = CH4_wet%total_out(i)! *          frac_inun_loc
        CO2_wet%out(i) = CO2_wet%total_out(i)! *          frac_inun_loc
         O2_wet%out(i) =  O2_wet%total_out(i)! *          frac_inun_loc
        layer_wt(i) = REAL(mxlay_above_wt,wp)

!        CH4_tot%out          (i) = CH4_dry%out          (i) * (1._wp-frac_inun_loc) + CH4_wet%out          (i) * frac_inun_loc
!        CH4_tot%prod_out     (i) = CH4_dry%prod_out     (i) * (1._wp-frac_inun_loc) + CH4_wet%prod_out     (i) * frac_inun_loc
!        CH4_tot%oxid_out     (i) = CH4_dry%oxid_out     (i) * (1._wp-frac_inun_loc) + CH4_wet%oxid_out     (i) * frac_inun_loc
!        CH4_tot%diff_out     (i) = CH4_dry%diff_out     (i) * (1._wp-frac_inun_loc) + CH4_wet%diff_out     (i) * frac_inun_loc
!        CH4_tot%other_out    (i) =                                                    CH4_wet%other_out    (i) * frac_inun_loc
!        CH4_tot%plant_out    (i) = CH4_dry%plant_out    (i) * (1._wp-frac_inun_loc) + CH4_wet%plant_out    (i) * frac_inun_loc
!        CH4_tot%total_out    (i) = CH4_dry%total_out    (i) * (1._wp-frac_inun_loc) + CH4_wet%total_out    (i) * frac_inun_loc
!        CH4_tot%total_initial(i) = CH4_dry%total_initial(i) * (1._wp-frac_inun_loc) + CH4_wet%total_initial(i) * frac_inun_loc
!        CH4_tot%total_final  (i) = CH4_dry%total_final  (i) * (1._wp-frac_inun_loc) + CH4_wet%total_final  (i) * frac_inun_loc

      ENDDO ! i = landpoints

      ! Other methane processes

      !   Fire
      IF (dsl4jsb_Config(CH4_)%do_fire) THEN
        CALL calc_methane_fire(co2flux_fire_all_2_atm(:),lt,CH4_out_fire(:),fire_em_NO(:),fire_em_CO(:),fire_em_iso(:))
      ENDIF

      !   Termites
      ! TODO, StW: The way the termite CH4 emissions are calculated for crops depends on the processing order of the tiles.
      ! Currently it only gives expected results if crop tiles are processed last.
      ! This is unsatisfactory and should be somehow replaced!!!
      IF ((dsl4jsb_Config(CH4_)%do_termite) .AND. (min_mmtemp20(1) < 1000._wp)) THEN

        ! NLCC_ variables (TODO: What if dynveg/NLCC is not active?)
        dsl4jsb_Get_var2D_onChunk(NLCC_, min_mmtemp20)

        l_GPP_dep = dsl4jsb_Lctlib_param(PhenologyType)==3 .OR. dsl4jsb_Lctlib_param(C4flag) ! Crops are separate anyway
        IF (GDD_prev_year(1) < 0._wp) THEN
          GDD => GDD_sum      ! 1st year: Use sum during current year up to now
        ELSE
          GDD => GDD_prev_year ! Any other year: Use sum of previous year
        ENDIF
        CALL calc_methane_termites(dsl4jsb_Lctlib_param(CropFlag),                                     &
                                   frac_tile(:), GPP_prev_year(:), GDD(:), min_mmtemp20(:), l_GPP_dep, &
                                   termite_biomass(:), CH4_out_termite(:),                             &
                                   frac_non_crop(ics:ice,iblk), sum_non_crop(ics:ice,iblk))
        IF (.NOT. ASSOCIATED(tile%next_sibling_tile)) THEN
          frac_non_crop(ics:ice,iblk) = 0._wp
          sum_non_crop (ics:ice,iblk) = 0._wp
        ENDIF
        CH4_out_termite(:) = CH4_out_termite(:) * dsl4jsb_Config(CH4_)%scale_termite
      ENDIF
 
      !   Herbivores
      IF (dsl4jsb_Config(CH4_)%do_herbivore) THEN
        CALL calc_methane_herbivores(cflux_herb_2_atm(:),dsl4jsb_Lctlib_param(woody_pft),CH4_out_herbivore(:))
        CH4_out_herbivore(:) = CH4_out_herbivore(:) * dsl4jsb_Config(CH4_)%scale_herbivore
      ENDIF

    ELSE ! Vegetation tile
      ! Should never happen, since methane process only runs on vegetation tiles
    ENDIF

    IF (lissue) THEN
      CALL message(routine,TRIM(tile%name)//': Tridiagonal matrix solving failed for one or more points',all_print = .TRUE.)
    ENDIF

    ! Deallocate local work arrays
    !
    DEALLOCATE(OxicDecompCarbon)
    DEALLOCATE(AnoxDecompCarbon)
    DEALLOCATE(water_table)
    DEALLOCATE(tiller_area_sl)
    DEALLOCATE(frac_tile)
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

  END SUBROUTINE update_methane

  SUBROUTINE aggregate_methane(tile, options)
    USE mo_jsb_time,          ONLY: is_newday, timesteps_per_day, is_time_experiment_start
    USE mo_ch4_memory_class,  ONLY: t_chem_ptrs

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    TYPE(t_jsb_model), POINTER :: model
    dsl4jsb_Def_config(CH4_)
    dsl4jsb_Def_memory(CH4_)

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

!    dsl4jsb_Aggregate_onChunk(CH4_, CH4_tot%out          , weighted_by_fract)
!    dsl4jsb_Aggregate_onChunk(CH4_, CH4_tot%prod_out     , weighted_by_fract)
!    dsl4jsb_Aggregate_onChunk(CH4_, CH4_tot%oxid_out     , weighted_by_fract)
!    dsl4jsb_Aggregate_onChunk(CH4_, CH4_tot%diff_out     , weighted_by_fract)
!    dsl4jsb_Aggregate_onChunk(CH4_, CH4_tot%other_out    , weighted_by_fract)
!    dsl4jsb_Aggregate_onChunk(CH4_, CH4_tot%plant_out    , weighted_by_fract)
!    dsl4jsb_Aggregate_onChunk(CH4_, CH4_tot%total_out    , weighted_by_fract)
!    dsl4jsb_Aggregate_onChunk(CH4_, CH4_tot%total_initial, weighted_by_fract)
!    dsl4jsb_Aggregate_onChunk(CH4_, CH4_tot%total_final  , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CH4_, CH4_air, weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CH4_, CO2_air, weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CH4_, O2_air, weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CH4_, mxlay_above_rock, weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CH4_, layer_wt, weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CH4_, root_dist_sl, weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CH4_, CH4_out_net, weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CH4_, CO2_out_net, weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CH4_, CH4_prod_ratio_out_wet, weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CH4_, CH4_dry%out, weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CH4_, CH4_dry%prod_out, weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CH4_, CH4_dry%oxid_out, weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CH4_, CH4_dry%diff_out, weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CH4_, CH4_dry%plant_out, weighted_by_fract)
!    dsl4jsb_Aggregate_onChunk(CH4_, CH4_dry%excess, weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CH4_, CH4_dry%total_out, weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CH4_, CH4_dry%total_initial, weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CH4_, CH4_dry%total_final, weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CH4_, CH4_dry%dens_sl, weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CH4_, CH4_wet%out, weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CH4_, CH4_wet%prod_out, weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CH4_, CH4_wet%oxid_out, weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CH4_, CH4_wet%diff_out, weighted_by_fract)
!    dsl4jsb_Aggregate_onChunk(CH4_, CH4_wet%excess, weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CH4_, CH4_wet%other_out, weighted_by_fract) ! Ebullition
    dsl4jsb_Aggregate_onChunk(CH4_, CH4_wet%plant_out, weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CH4_, CH4_wet%total_out, weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CH4_, CH4_wet%total_initial, weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CH4_, CH4_wet%total_final, weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CH4_, CH4_wet%dens_sl, weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CH4_, CO2_dry%out, weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CH4_, CO2_dry%prod_out, weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CH4_, CO2_dry%oxid_out, weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CH4_, CO2_dry%diff_out, weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CH4_, CO2_dry%plant_out, weighted_by_fract)
!    dsl4jsb_Aggregate_onChunk(CH4_, CO2_dry%excess, weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CH4_, CO2_dry%total_out, weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CH4_, CO2_dry%total_initial, weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CH4_, CO2_dry%total_final, weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CH4_, CO2_dry%dens_sl, weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CH4_, CO2_wet%out, weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CH4_, CO2_wet%prod_out, weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CH4_, CO2_wet%oxid_out, weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CH4_, CO2_wet%diff_out, weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CH4_, CO2_wet%plant_out, weighted_by_fract)
!    dsl4jsb_Aggregate_onChunk(CH4_, CO2_wet%excess, weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CH4_, CO2_wet%total_out, weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CH4_, CO2_wet%total_initial, weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CH4_, CO2_wet%total_final, weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CH4_, CO2_wet%dens_sl, weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CH4_, O2_dry%out, weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CH4_, O2_dry%prod_out, weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CH4_, O2_dry%oxid_out, weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CH4_, O2_dry%diff_out, weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CH4_, O2_dry%other_out, weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CH4_, O2_dry%plant_out, weighted_by_fract)
!    dsl4jsb_Aggregate_onChunk(CH4_, O2_dry%excess, weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CH4_, O2_dry%total_out, weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CH4_, O2_dry%total_initial, weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CH4_, O2_dry%total_final, weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CH4_, O2_dry%dens_sl, weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CH4_, O2_wet%out, weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CH4_, O2_wet%prod_out, weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CH4_, O2_wet%oxid_out, weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CH4_, O2_wet%diff_out, weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CH4_, O2_wet%other_out, weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CH4_, O2_wet%plant_out, weighted_by_fract)
!    dsl4jsb_Aggregate_onChunk(CH4_, O2_wet%excess, weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CH4_, O2_wet%total_out, weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CH4_, O2_wet%total_initial, weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CH4_, O2_wet%total_final, weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CH4_, O2_wet%dens_sl, weighted_by_fract)
    IF (dsl4jsb_Config(CH4_)%do_fire) THEN
      dsl4jsb_Aggregate_onChunk(CH4_, CH4_out_fire, weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(CH4_, fire_em_NO, weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(CH4_, fire_em_CO, weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(CH4_, fire_em_iso, weighted_by_fract)
    ENDIF
    IF (dsl4jsb_Config(CH4_)%do_termite) THEN
      dsl4jsb_Aggregate_onChunk(CH4_, CH4_out_termite, weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(CH4_, termite_biomass, weighted_by_fract)
    ENDIF
    IF (dsl4jsb_Config(CH4_)%do_herbivore) THEN
      dsl4jsb_Aggregate_onChunk(CH4_, CH4_out_herbivore, weighted_by_fract)
    ENDIF
! >>> Debug
    dsl4jsb_Aggregate_onChunk(CH4_,CH4_wet%prod_out_sl, weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CH4_,CH4_wet%oxid_out_sl, weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CH4_,CH4_wet%other_out_sl, weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CH4_,CH4_wet%plant_out_sl, weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CH4_,CO2_wet%prod_out_sl, weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CH4_,CO2_wet%oxid_out_sl, weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CH4_,CO2_wet%plant_out_sl, weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CH4_, O2_wet%prod_out_sl, weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CH4_, O2_wet%oxid_out_sl, weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CH4_, O2_wet%other_out_sl, weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CH4_, O2_wet%plant_out_sl, weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CH4_,CH4_dry%prod_out_sl, weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CH4_,CH4_dry%oxid_out_sl, weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CH4_,CH4_dry%other_out_sl, weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CH4_,CH4_dry%plant_out_sl, weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CH4_,CO2_dry%prod_out_sl, weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CH4_,CO2_dry%oxid_out_sl, weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CH4_,CO2_dry%plant_out_sl, weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CH4_, O2_dry%prod_out_sl, weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CH4_, O2_dry%oxid_out_sl, weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CH4_, O2_dry%other_out_sl, weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CH4_, O2_dry%plant_out_sl, weighted_by_fract)
! <<<
    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE aggregate_methane

#endif
END MODULE mo_methane_interface
