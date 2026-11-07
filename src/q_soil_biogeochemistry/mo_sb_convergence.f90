!> Restart-persistent, report-only coupled SOC convergence on complete forcing cycles
!>
!> ICON-Land
!>
!> ---------------------------------------
!> Copyright (C) 2013-2026, MPI-M, MPI-BGC
!> Contact: icon-model.org
!> Authors: AUTHORS.md
!> See LICENSES/ for license information
!> SPDX-License-Identifier: BSD-3-Clause
!> ---------------------------------------
MODULE mo_sb_convergence
#ifndef __NO_JSBACH__
  USE mo_kind, ONLY: lp
  USE mo_exception, ONLY: message, finish
  USE mo_jsb_class, ONLY: Get_model
  USE mo_jsb_model_class, ONLY: t_jsb_model, MODEL_QUINCY
  USE mo_jsb_tile_class, ONLY: t_jsb_tile_abstract
  USE mo_jsb_process_class, ONLY: SB_, HYDRO_
  USE mo_sb_memory_class, ONLY: t_sb_memory
  USE mo_sb_convergence_kernel, ONLY: stock_relative_change, soc_window_pass, &
    & soc_relative_tolerance, soc_passing_fraction
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: init_sb_convergence, report_sb_convergence
  ! History stores five g/m2 inventories, year, event serial, active fraction,
  ! first qualifying year, forcing origin, calendar signature and inventory version.
  ! Metadata is exact integer-valued real where relevant.
  INTEGER, PARAMETER :: stride = 12 !< Fields per annual grid-cell sample [-]
CONTAINS
  !> Called by the standalone driver after reading forcing, BEFORE restart input.
  !! No additional forcing cache: only 12*(N+1) scalars per grid cell on the root.
  !! Changing N changes restart shape and requires a fresh compatible history.
  SUBROUTINE init_sb_convergence(model_id, cycle_years, forcing_start_year)
    USE mo_jsb_grid, ONLY: Get_grid, Register_vgrid
    USE mo_jsb_grid_class, ONLY: t_jsb_grid, t_jsb_vgrid, new_vgrid
    USE mo_jsb_io, ONLY: ZAXIS_GENERIC, t_cf, t_grib1, t_grib2, tables, grib_bits
    USE mo_jsb_varlist, ONLY: NONE
    USE mo_sb_config_class, ONLY: t_sb_config
    INTEGER, INTENT(in) :: model_id, cycle_years !< Model index and configured climate cycle [years]
    INTEGER, INTENT(in) :: forcing_start_year !< First year of cyclic climate input [calendar year]
    TYPE(t_jsb_model), POINTER :: model
    CLASS(t_jsb_tile_abstract), POINTER :: root
    TYPE(t_jsb_grid), POINTER :: grid
    TYPE(t_jsb_vgrid), POINTER :: history_grid
    TYPE(t_sb_memory), POINTER :: mem
    CHARACTER(len=40) :: grid_name
    CHARACTER(len=256) :: text
    model => Get_model(model_id)
    IF (model%config%model_scheme /= MODEL_QUINCY .OR. .NOT. model%config%l_spinup) RETURN
    IF (.NOT. model%processes(SB_)%p%config%active) RETURN
    SELECT TYPE (conf => model%processes(SB_)%p%config)
    TYPE IS (t_sb_config)
      conf%sb_convergence_cycle_years = cycle_years
      IF (conf%flag_sb_model_spinup_accelerator .AND. cycle_years <= 0) &
        & CALL finish('init_sb_convergence', 'SOC replay convergence requires positive cyclic_nyears.')
    END SELECT
    IF (cycle_years <= 0) THEN
      CALL message('init_sb_convergence', 'No cyclic climate forcing: coupled SOC monitor disabled.')
      RETURN
    END IF
    CALL model%Get_top_tile(root)
    SELECT TYPE (m => root%mem(SB_)%p)
    TYPE IS (t_sb_memory)
      mem => m
    CLASS DEFAULT
      CALL finish('init_sb_convergence', 'Missing root SB memory.')
    END SELECT
    grid => Get_grid(model%grid_id)
    WRITE(grid_name,'(A,I0)') 'sb_soc_history_',model_id
    history_grid => new_vgrid(TRIM(grid_name), ZAXIS_GENERIC, stride*(cycle_years+1), &
      & longname='annual SOC convergence stocks and metadata', units='1')
    CALL Register_vgrid(history_grid)
    CALL mem%Add_var('convergence_history', mem%convergence_history, grid, history_grid, &
      & t_cf('convergence_history', '1', 'packed annual stocks in g/m2 and convergence metadata'), &
      & t_grib1(tables(1),255,grib_bits), t_grib2(255,255,255,grib_bits), 'sb', TRIM(root%name), &
      & output_level=NONE, loutput=.FALSE., lrestart=.TRUE., lrestart_cont=.TRUE., initval_r=-1._lp)
    mem%convergence_cycle_years = cycle_years
    mem%convergence_forcing_start_year = forcing_start_year
    mem%convergence_history%ptr = -1._lp
    WRITE(text,*) 'Report-only SOC monitoring: model, cycle years, tolerance, passing fraction = ', &
      & model_id, cycle_years, soc_relative_tolerance, soc_passing_fraction
    CALL message('init_sb_convergence', TRIM(text))
    CALL message('init_sb_convergence', &
      & 'Missing restart convergence history starts a new pending window; scientific replay state is unchanged.')
    CALL message('init_sb_convergence', &
      & 'Replay convergence uses cyclic_nyears; legacy sb_model_spinup_accelerator_convergence_years is ignored.')
  END SUBROUTINE init_sb_convergence

  !> Sample after all tile tasks and artificial events, at the same annual phase.
  !! Depth integration precedes tile weighting. Use only calculated SB tiles, not
  !! both their parents and children. Root/child fractions are absolute grid fractions.
  !! Event-containing windows retain their measured drift but cannot claim clean
  !! coupled equilibrium. Nutrients never veto SOC. Nothing here stops a run.
  SUBROUTINE report_sb_convergence(model_id, current_datetime, dtime)
    USE mo_jsb_time, ONLY: t_datetime, is_newyear, get_year, get_year_length
    USE mo_jsb_grid, ONLY: Get_grid
    USE mo_jsb_grid_class, ONLY: t_jsb_grid
    USE mo_hydro_memory_class, ONLY: t_hydro_memory
    USE mo_lnd_bgcm_store_class, ONLY: SB_BGCM_POOL_ID
    USE mo_lnd_bgcm_idx, ONLY: ixC, ixN, ixP
    USE mo_q_sb_update_pools, ONLY: calculate_sb_pool_bg_soil_for_element
    USE mo_jsb_physical_constants, ONLY: molar_mass_C, molar_mass_N, molar_mass_P
    USE mo_jsb_parallel, ONLY: sum_land_diagnostic
    INTEGER, INTENT(in) :: model_id !< Independent land-model/domain identifier [-]
    TYPE(t_datetime), POINTER, INTENT(in) :: current_datetime !< Unmodified global model time
    REAL(lp), INTENT(in) :: dtime !< Coupled timestep [s], never artificial replay duration
    TYPE(t_jsb_model), POINTER :: model
    CLASS(t_jsb_tile_abstract), POINTER :: root, tile
    TYPE(t_jsb_grid), POINTER :: grid
    TYPE(t_sb_memory), POINTER :: mem, soil
    TYPE(t_hydro_memory), POINTER :: hydro
    REAL(lp), POINTER :: history(:,:,:), pool(:,:,:,:)
    REAL(lp), ALLOCATABLE :: stocks(:,:,:), fraction(:,:), events(:,:), layer(:,:)
    REAL(lp), ALLOCATABLE :: active(:,:), passed(:,:), eligible(:,:), pending(:,:), duplicate(:,:)
    REAL(lp), ALLOCATABLE :: affected(:,:), threshold_pass(:,:)
    REAL(lp), ALLOCATABLE :: drift(:,:,:) !< Signed endpoint changes for valid windows [g m-2]
    REAL(lp) :: old(5), change(5), relative(5), weights(5), area_total, nactive, npass, first_pass
    REAL(lp) :: pass_area, npending, naffected, raw_pass, area_fraction
    REAL(lp) :: global_drift(5), valid_area !< Valid-window inventory changes [g] and area [m2]
    INTEGER :: year, n, slot, old_slot, latest_slot, it, b, c, e, nz, nc, j
    LOGICAL :: has_history, event_affected, stock_pass
    CHARACTER(len=1024) :: text
    model => Get_model(model_id)
    IF (model%config%model_scheme /= MODEL_QUINCY .OR. .NOT. model%config%l_spinup) RETURN
    IF (.NOT. model%processes(SB_)%p%config%active) RETURN
    IF (.NOT. is_newyear(current_datetime,dtime)) RETURN
    CALL model%Get_top_tile(root)
    SELECT TYPE (m => root%mem(SB_)%p)
    TYPE IS (t_sb_memory)
      mem => m
    CLASS DEFAULT
      RETURN
    END SELECT
    n = mem%convergence_cycle_years
    IF (n <= 0) RETURN
    grid => Get_grid(model%grid_id)
    history => mem%convergence_history%ptr
    year = get_year(current_datetime)
    slot = stride*MODULO(year,n+1)
    old_slot = stride*MODULO(year-n,n+1)
    ALLOCATE(stocks(grid%nproma,5,grid%nblks),fraction(grid%nproma,grid%nblks), &
      & events(grid%nproma,grid%nblks),active(grid%nproma,grid%nblks), &
      & passed(grid%nproma,grid%nblks),eligible(grid%nproma,grid%nblks), &
      & pending(grid%nproma,grid%nblks),duplicate(grid%nproma,grid%nblks), &
      & affected(grid%nproma,grid%nblks),threshold_pass(grid%nproma,grid%nblks), &
      & drift(grid%nproma,5,grid%nblks))
    stocks=0._lp; fraction=0._lp; events=0._lp
    active=0._lp; passed=0._lp; eligible=0._lp; pending=0._lp; duplicate=0._lp
    affected=0._lp; threshold_pass=0._lp
    drift=0._lp
    weights = [molar_mass_C,molar_mass_N,molar_mass_P,molar_mass_N,molar_mass_P]
    DO it=1,SIZE(model%tiles)
      tile => model%tiles(it)%p
      IF (.NOT. tile%Is_process_calculated(SB_)) CYCLE
      IF (tile%lcts(1)%lib_id == 0) CYCLE
      SELECT TYPE (m => tile%mem(SB_)%p)
      TYPE IS (t_sb_memory)
        soil => m
      END SELECT
      SELECT TYPE (h => tile%mem(HYDRO_)%p)
      TYPE IS (t_hydro_memory)
        hydro => h
      END SELECT
      DO b=1,grid%nblks
        nc=grid%nproma
        IF (b == grid%nblks) nc=grid%npromz
        pool => tile%bgcm_store%Get_matrix_from_store_2l_3d(SB_BGCM_POOL_ID,1,nc,b,tile%name)
        nz=SIZE(pool,4)
        ALLOCATE(layer(nc,nz))
        DO e=1,3
          SELECT CASE(e)
          CASE(1)
            CALL calculate_sb_pool_bg_soil_for_element(pool,nz,ixC,layer)
          CASE(2)
            CALL calculate_sb_pool_bg_soil_for_element(pool,nz,ixN,layer)
          CASE(3)
            CALL calculate_sb_pool_bg_soil_for_element(pool,nz,ixP,layer)
          END SELECT
          DO c=1,nc
            IF (.NOT. grid%patch%cells%decomp_info%owner_mask(c,b) .OR. tile%fract(c,b) <= 0._lp) CYCLE
            stocks(c,e,b)=stocks(c,e,b)+tile%fract(c,b)*SUM(layer(c,:)*hydro%soil_depth_sl%ptr(c,:,b))
          END DO
        END DO
        DO c=1,nc
          IF (.NOT. grid%patch%cells%decomp_info%owner_mask(c,b) .OR. tile%fract(c,b) <= 0._lp) CYCLE
          fraction(c,b)=fraction(c,b)+tile%fract(c,b)
          events(c,b)=events(c,b)+soil%convergence_event_count%ptr(c,b)
          stocks(c,4,b)=stocks(c,4,b)+tile%fract(c,b)*SUM(hydro%soil_depth_sl%ptr(c,:,b)* &
            & (soil%nh4_solute%ptr(c,:,b)+soil%nh4_assoc%ptr(c,:,b)+soil%no3_solute%ptr(c,:,b)))
          stocks(c,5,b)=stocks(c,5,b)+tile%fract(c,b)*SUM(hydro%soil_depth_sl%ptr(c,:,b)* &
            & (soil%po4_solute%ptr(c,:,b)+soil%po4_assoc_fast%ptr(c,:,b)+soil%po4_assoc_slow%ptr(c,:,b)))
        END DO
        DEALLOCATE(layer)
      END DO
    END DO
    DO b=1,grid%nblks
      DO c=1,grid%nproma
        IF (fraction(c,b) <= 0._lp) CYCLE
        active(c,b)=1._lp
        stocks(c,:,b)=stocks(c,:,b)*weights/fraction(c,b)
        ! Calendar/forcing/inventory changes never silently reuse old history.
        DO j=0,n
          IF (history(c,stride*j+8,b) < 0._lp) CYCLE
          IF (history(c,stride*j+10,b) /= REAL(mem%convergence_forcing_start_year,lp) .OR. &
            & history(c,stride*j+11,b) /= REAL(get_year_length(2000),lp) .OR. &
            & history(c,stride*j+12,b) /= 1._lp) THEN
            history(c,:,b)=-1._lp
            EXIT
          END IF
        END DO
        IF (history(c,slot+8,b) >= 0._lp .AND. history(c,slot+6,b) == REAL(year,lp)) THEN
          duplicate(c,b)=1._lp
        END IF
      END DO
    END DO
    nactive=sum_land_diagnostic(active)
    IF (nactive == 0._lp) RETURN
    IF (sum_land_diagnostic(duplicate) == nactive) RETURN
    DO b=1,grid%nblks
      DO c=1,grid%nproma
        IF (active(c,b) == 0._lp) CYCLE
        first_pass=-1._lp
        latest_slot=-1
        DO j=0,n
          IF (history(c,stride*j+8,b) < 0._lp) CYCLE
          IF (latest_slot < 0) THEN
            latest_slot=stride*j
          ELSE IF (history(c,stride*j+6,b) > history(c,latest_slot+6,b)) THEN
            latest_slot=stride*j
          END IF
        END DO
        IF (latest_slot >= 0) first_pass=history(c,latest_slot+9,b)
        has_history=history(c,old_slot+8,b) >= 0._lp .AND. &
          & history(c,old_slot+6,b) == REAL(year-n,lp)
        event_affected=.FALSE.; stock_pass=.FALSE.
        relative=HUGE(1._lp); change=0._lp
        IF (has_history) THEN
          old=history(c,old_slot+1:old_slot+5,b)
          change=stocks(c,:,b)-old
          drift(c,:,b)=change
          eligible(c,b)=1._lp
          relative=stock_relative_change(old,stocks(c,:,b))
          event_affected=events(c,b) /= history(c,old_slot+7,b) .OR. &
            & fraction(c,b) /= history(c,old_slot+8,b)
          stock_pass=soc_window_pass(old(1),stocks(c,1,b),soc_relative_tolerance)
          IF (stock_pass) threshold_pass(c,b)=1._lp
          IF (event_affected) THEN
            affected(c,b)=1._lp
          ELSE
            IF (stock_pass) THEN
              passed(c,b)=1._lp
              IF (first_pass < 0._lp) first_pass=REAL(year,lp)
            END IF
          END IF
        ELSE
          pending(c,b)=1._lp
        END IF
        ! Global counts remain unconditional below. Per-cell inventories and
        ! nonpassing coordinates are detailed spinup diagnostics, not fatal errors.
        IF (model%config%flag_spinup_debug .AND. passed(c,b) == 0._lp) THEN
          WRITE(text,*) 'SOC monitor cell: model/year/ic/iblk/lon/lat=',model_id,year,c,b, &
            & grid%lon(c,b),grid%lat(c,b),'; history/event-affected/threshold-pass=', &
            & has_history,event_affected,stock_pass,'; stocks g/m2 C,N,P,SIN,SIP=',stocks(c,:,b), &
            & '; signed changes=',change,'; relative changes=',relative
          CALL message('report_sb_convergence',TRIM(text),all_print=.TRUE.)
        END IF
        history(c,slot+1:slot+5,b)=stocks(c,:,b)
        history(c,slot+6,b)=REAL(year,lp)
        history(c,slot+7,b)=events(c,b)
        history(c,slot+8,b)=fraction(c,b)
        history(c,slot+9,b)=first_pass
        history(c,slot+10,b)=REAL(mem%convergence_forcing_start_year,lp)
        history(c,slot+11,b)=REAL(get_year_length(2000),lp)
        history(c,slot+12,b)=1._lp
      END DO
    END DO
    npass=sum_land_diagnostic(passed)
    npending=sum_land_diagnostic(pending)
    naffected=sum_land_diagnostic(affected)
    raw_pass=sum_land_diagnostic(threshold_pass)
    area_total=sum_land_diagnostic(fraction*grid%area)
    pass_area=sum_land_diagnostic(passed*fraction*grid%area)
    area_fraction=0._lp
    IF (area_total > 0._lp) area_fraction=pass_area/area_total
    WRITE(text,*) 'SOC report-only: model/year/N=',model_id,year,n, &
      & '; active/clean-passing/pending/event-affected/threshold-only=', &
      & nactive,npass,npending,naffected,raw_pass,'; passing count/area fraction=', &
      & npass/nactive,area_fraction,'; 98-percent SOC status=',npass/nactive >= soc_passing_fraction, &
      & '; nutrient diagnostics only; no stop requested.'
    CALL message('report_sb_convergence',TRIM(text))
    valid_area=sum_land_diagnostic(eligible*fraction*grid%area)
    DO e=1,5
      global_drift(e)=sum_land_diagnostic(drift(:,e,:)*fraction*grid%area)
    END DO
    IF (valid_area > 0._lp) THEN
      WRITE(text,*) 'SOC monitor valid-window area-mean signed changes [g/m2] C,N,P,SIN,SIP=', &
        & global_drift/valid_area,'; includes flagged event windows; nutrients diagnostic only.'
      CALL message('report_sb_convergence',TRIM(text))
    END IF
  END SUBROUTINE report_sb_convergence
#endif
END MODULE mo_sb_convergence
