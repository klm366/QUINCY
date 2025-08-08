!> Initialization of the lateral (intertile) exchange memory.
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
!>#### Initialisation and reading of lex (lateral [intertile] exchange) variables.
!>
MODULE mo_lex_init
#ifndef __NO_JSBACH__

  USE mo_kind,               ONLY: wp
  USE mo_exception,          ONLY: message, finish
  USE mo_jsb_control,        ONLY: debug_on
  USE mo_util,               ONLY: int2string

  USE mo_jsb_model_class,    ONLY: t_jsb_model
  USE mo_jsb_grid_class,     ONLY: t_jsb_grid, t_jsb_vgrid
  USE mo_jsb_grid,           ONLY: Get_grid, Get_vgrid
  USE mo_jsb_tile_class,     ONLY: t_jsb_tile_abstract
  USE mo_jsb_class,          ONLY: get_model

  ! Use of processes in this module
  dsl4jsb_Use_processes LEX_, HYDRO_

  ! Use of process configurations
  dsl4jsb_Use_config(LEX_)
  dsl4jsb_Use_config(HYDRO_)

  ! Use of process memories
  dsl4jsb_Use_memory(LEX_)
  dsl4jsb_Use_memory(HYDRO_)

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: lex_init

  CHARACTER(len=*), PARAMETER :: modname = 'mo_lex_init'

CONTAINS

  !
  !> Intialize lateral exchange process (after memory has been set up)
  !
  SUBROUTINE lex_init(tile)

    USE mo_lex_process,        ONLY: check_flow_matrix, get_geometry, &
      &                              get_flow_lags
    USE mo_hydro_constants,    ONLY: vel_lat_srf, vel_lat_blg

    USE mo_lex_connectivities, ONLY: &
      & cmx_micro_ctc_2tilelex, cmx_micro_hfp_srf_2tilelex, cmx_micro_hfp_blg_2tilelex, &
      & cmx_macro_hfp_srf_2tilelex, cmx_macro_hfp_blg_2tilelex,                         &
      & cmx_micro_ctc_3tilelex, cmx_micro_hfp_srf_3tilelex, cmx_micro_hfp_blg_3tilelex, &
      & cmx_macro_hfp_srf_3tilelex, cmx_macro_hfp_blg_3tilelex,                         &
      & cmx_micro_ctc_topolex,  cmx_micro_hfp_srf_topolex, cmx_micro_hfp_blg_topolex,   &
      & cmx_macro_hfp_srf_topolex, cmx_macro_hfp_blg_topolex

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile ! Tile with LEX_ process running


    dsl4jsb_Def_config(LEX_)
    dsl4jsb_Def_config(HYDRO_)

    dsl4jsb_Def_memory(LEX_)
    dsl4jsb_Def_memory(HYDRO_)
    dsl4jsb_Def_memory_tile(HYDRO_, sibl)

    CLASS(t_jsb_tile_abstract), POINTER :: sibl,       & ! Sibling tile with LEX_ process running
      &                                    dummy         ! Dummy tile
    TYPE(t_jsb_model), POINTER :: model
    TYPE(t_jsb_grid),  POINTER :: grid
    TYPE(t_jsb_vgrid), POINTER :: connecttiles

    LOGICAL :: has_loop     !< .TRUE. Flow matrix contains closed loop

    INTEGER :: nproma,   &  !< nproma
      &        nblks,    &  !< nblocks
      &        nsibl,    &  !< number of siblings on current tile level
      &        nconti,   &  !< number of connected tiles on current tile level
      &        it,       &  !< counter for loop over tiles
      &        is,       &  !< counter for loop over tiles when within loop over tiles
      &        src_id,   &  !< source lex_id
      &        trgt_id,  &  !< target lex_id
      &        lex_idt,  &  !< lex id of of tile
      &        lex_ids,  &  !< lex id of sibling
      &        lex_idd      !< lex id of dummy

    REAL(wp), POINTER :: &
      &        test_matrix(:,:),       & !< matrix used to check for inconsitances in base matrices
      &        cmx_micro_ctc(:,:),     & !< base matrix for micro-scale connectivity via direct contact
      &        cmx_micro_hfp_srf(:,:), & !< base matrix for micro-scale connectivity via surface hydrological flow paths
      &        cmx_micro_hfp_blg(:,:), & !< base matrix for micro-scale connectivity via below-ground hydrological flow paths
      &        cmx_macro_hfp_srf(:,:), & !< base matrix for macro-scale connectivity via surface hydrological flow paths
      &        cmx_macro_hfp_blg(:,:)    !< base matrix for macro-scale connectivity via below-ground hydrological flow paths

    CHARACTER(len=*), PARAMETER :: routine = modname//':lex_init'

    ! Pointers to variables in memory for tile
    dsl4jsb_Real2D_onDomain :: &
      & h_surf,                &
      & area,                  &
      & macro_slope,           &
      & slope_eff,             &
      & fact_lag_srf,          &
      & fact_lag_blg

    dsl4jsb_Real3D_onDomain :: &
      & h_offset_ct,           &
      & length_scale_ct,       &
      & length_iface_ct,       &
      & micro_slope_ct,        &
      & fract_latflow_srf_ct,  &
      & fract_latflow_blg_ct,  &
      & microcon_ctc_ct,       &
      & microcon_hfp_srf_ct,   &
      & microcon_hfp_blg_ct,   &
      & macrocon_hfp_srf_ct,   &
      & macrocon_hfp_blg_ct

    ! Pointers to variables in memory for sibling
    dsl4jsb_Real2D_onDomain :: &
      & h_surf_sibl,           &
      & area_sibl

    ! Locally allocated vectors
    dsl4jsb_Real2D_onDomain :: &
      & fract,                 &
      & fract_sibl,            &
      & hlp1,                  &
      & hlp2,                  &
      & hlp3,                  &
      & hlp4,                  &
      & fract_micro_hfp_srf,   &
      & fract_micro_hfp_blg,   &
      & fract_macro_hfp_srf,   &
      & fract_macro_hfp_blg,   &
      & fact_lag_micro_srf,    &
      & fact_lag_micro_blg,    &
      & fact_lag_macro_srf,    &
      & fact_lag_macro_blg

    dsl4jsb_Real3D_onDomain :: &
      & fracts,                &
      & area_nclsd,            &
      & fract_nclsd,           &
      & fact_lag_srf_ct,       &
      & fact_lag_blg_ct

    IF (.NOT. tile%Is_process_active(LEX_)) RETURN

    IF (debug_on()) CALL message(routine, 'Initializing variables for lateral exchange on tile ' // &
      & TRIM(tile%name) // ' (lex_id=' // TRIM(int2string(tile%lcts(1)%lex_id)) // ')')

    lex_idt = tile%lcts(1)%lex_id
    IF (lex_idt > 0) THEN

      dsl4jsb_Get_memory(LEX_)
      dsl4jsb_Get_memory(HYDRO_)

      dsl4jsb_Get_var2D_onDomain(HYDRO_, h_surf)               ! in
      dsl4jsb_Get_var2D_onDomain(HYDRO_, area)                 ! in
      dsl4jsb_Get_var2D_onDomain(HYDRO_, macro_slope)          ! in

      dsl4jsb_Get_var2D_onDomain(LEX_,   fact_lag_srf)         ! out
      dsl4jsb_Get_var2D_onDomain(LEX_,   fact_lag_blg)         ! out
      dsl4jsb_Get_var2D_onDomain(LEX_,   slope_eff)            ! out

      dsl4jsb_Get_var3D_onDomain(LEX_,   h_offset_ct)          ! out
      dsl4jsb_Get_var3D_onDomain(LEX_,   length_scale_ct)      ! out
      dsl4jsb_Get_var3D_onDomain(LEX_,   length_iface_ct)      ! out
      dsl4jsb_Get_var3D_onDomain(LEX_,   micro_slope_ct)       ! out
      dsl4jsb_Get_var3D_onDomain(LEX_,   microcon_ctc_ct)      ! out
      dsl4jsb_Get_var3D_onDomain(LEX_,   microcon_hfp_srf_ct)  ! out
      dsl4jsb_Get_var3D_onDomain(LEX_,   microcon_hfp_blg_ct)  ! out
      dsl4jsb_Get_var3D_onDomain(LEX_,   macrocon_hfp_srf_ct)  ! out
      dsl4jsb_Get_var3D_onDomain(LEX_,   macrocon_hfp_blg_ct)  ! out
      dsl4jsb_Get_var3D_onDomain(LEX_,   fract_latflow_srf_ct) ! out
      dsl4jsb_Get_var3D_onDomain(LEX_,   fract_latflow_blg_ct) ! out


      model  => get_model(tile%owner_model_id)
      grid   => Get_grid(model%grid_id)
      nproma = grid%Get_nproma()
      nblks  = grid%Get_nblks()
      nsibl  = tile%parent_tile%Get_no_of_children()
      connecttiles => Get_vgrid('connecttiles')
      nconti = connecttiles%n_levels

      ALLOCATE(fract(nproma, nblks))
      ALLOCATE(fract_sibl(nproma, nblks))
      ALLOCATE(test_matrix(nconti,nconti))
      ALLOCATE(cmx_micro_ctc(nconti,nconti))
      ALLOCATE(cmx_micro_hfp_srf(nconti,nconti))
      ALLOCATE(cmx_micro_hfp_blg(nconti,nconti))
      ALLOCATE(cmx_macro_hfp_srf(nconti,nconti))
      ALLOCATE(cmx_macro_hfp_blg(nconti,nconti))
      ALLOCATE(fract_micro_hfp_srf(nproma, nblks))
      ALLOCATE(fract_micro_hfp_blg(nproma, nblks))
      ALLOCATE(fract_macro_hfp_srf(nproma, nblks))
      ALLOCATE(fract_macro_hfp_blg(nproma, nblks))

      ! Get config
      dsl4jsb_Get_config(LEX_)
      dsl4jsb_Get_config(HYDRO_)

      ! Get the general connection matrixes
      IF  (TRIM(model%config%usecase) == 'jsbach_topolex') THEN

        cmx_micro_ctc = reshape(cmx_micro_ctc_topolex, shape(cmx_micro_ctc), order=(/2,1/))
        cmx_micro_hfp_srf = reshape(cmx_micro_hfp_srf_topolex, shape(cmx_micro_hfp_srf), order=(/2,1/))
        cmx_macro_hfp_srf = reshape(cmx_macro_hfp_srf_topolex, shape(cmx_macro_hfp_srf), order=(/2,1/))

        cmx_micro_hfp_blg = reshape(cmx_micro_hfp_blg_topolex, shape(cmx_micro_hfp_blg), order=(/2,1/))
        cmx_macro_hfp_blg = reshape(cmx_macro_hfp_blg_topolex, shape(cmx_macro_hfp_blg), order=(/2,1/))

      ELSE IF  (TRIM(model%config%usecase) == 'jsbach_topolex_sebt') THEN

        cmx_micro_ctc = reshape(cmx_micro_ctc_topolex, shape(cmx_micro_ctc), order=(/2,1/))
        cmx_micro_hfp_srf = reshape(cmx_micro_hfp_srf_topolex, shape(cmx_micro_hfp_srf), order=(/2,1/))
        cmx_macro_hfp_srf = reshape(cmx_macro_hfp_srf_topolex, shape(cmx_macro_hfp_srf), order=(/2,1/))

        cmx_micro_hfp_blg = reshape(cmx_micro_hfp_blg_topolex, shape(cmx_micro_hfp_blg), order=(/2,1/))
        cmx_macro_hfp_blg = reshape(cmx_macro_hfp_blg_topolex, shape(cmx_macro_hfp_blg), order=(/2,1/))

      ELSE IF (TRIM(model%config%usecase) == 'jsbach_2tilelex') THEN

        cmx_micro_ctc = reshape(cmx_micro_ctc_2tilelex, shape(cmx_micro_ctc), order=(/2,1/))
        cmx_micro_hfp_srf = reshape(cmx_micro_hfp_srf_2tilelex, shape(cmx_micro_hfp_srf), order=(/2,1/))
        cmx_macro_hfp_srf = reshape(cmx_macro_hfp_srf_2tilelex, shape(cmx_macro_hfp_srf), order=(/2,1/))

        cmx_micro_hfp_blg = reshape(cmx_micro_hfp_blg_2tilelex, shape(cmx_micro_hfp_blg), order=(/2,1/))
        cmx_macro_hfp_blg = reshape(cmx_macro_hfp_blg_2tilelex, shape(cmx_macro_hfp_blg), order=(/2,1/))

      ELSE IF (TRIM(model%config%usecase) == 'jsbach_2tilelex_sebt') THEN

        cmx_micro_ctc = reshape(cmx_micro_ctc_2tilelex, shape(cmx_micro_ctc), order=(/2,1/))
        cmx_micro_hfp_srf = reshape(cmx_micro_hfp_srf_2tilelex, shape(cmx_micro_hfp_srf), order=(/2,1/))
        cmx_macro_hfp_srf = reshape(cmx_macro_hfp_srf_2tilelex, shape(cmx_macro_hfp_srf), order=(/2,1/))

        cmx_micro_hfp_blg = reshape(cmx_micro_hfp_blg_2tilelex, shape(cmx_micro_hfp_blg), order=(/2,1/))
        cmx_macro_hfp_blg = reshape(cmx_macro_hfp_blg_2tilelex, shape(cmx_macro_hfp_blg), order=(/2,1/))

      ELSE IF (TRIM(model%config%usecase) == 'jsbach_3tilelex') THEN

        cmx_micro_ctc = reshape(cmx_micro_ctc_3tilelex, shape(cmx_micro_ctc), order=(/2,1/))
        cmx_micro_hfp_srf = reshape(cmx_micro_hfp_srf_3tilelex, shape(cmx_micro_hfp_srf), order=(/2,1/))
        cmx_macro_hfp_srf = reshape(cmx_macro_hfp_srf_3tilelex, shape(cmx_macro_hfp_srf), order=(/2,1/))

        cmx_micro_hfp_blg = reshape(cmx_micro_hfp_blg_3tilelex, shape(cmx_micro_hfp_blg), order=(/2,1/))
        cmx_macro_hfp_blg = reshape(cmx_macro_hfp_blg_3tilelex, shape(cmx_macro_hfp_blg), order=(/2,1/))

      ELSE IF (TRIM(model%config%usecase) == 'jsbach_3tilelex_sebt') THEN

        cmx_micro_ctc = reshape(cmx_micro_ctc_3tilelex, shape(cmx_micro_ctc), order=(/2,1/))
        cmx_micro_hfp_srf = reshape(cmx_micro_hfp_srf_3tilelex, shape(cmx_micro_hfp_srf), order=(/2,1/))
        cmx_macro_hfp_srf = reshape(cmx_macro_hfp_srf_3tilelex, shape(cmx_macro_hfp_srf), order=(/2,1/))

        cmx_micro_hfp_blg = reshape(cmx_micro_hfp_blg_3tilelex, shape(cmx_micro_hfp_blg), order=(/2,1/))
        cmx_macro_hfp_blg = reshape(cmx_macro_hfp_blg_3tilelex, shape(cmx_macro_hfp_blg), order=(/2,1/))

      ELSE
        CALL finish(TRIM(routine), '*** Lateral connection matrix not defined for usecase')
      END IF

      ! Check the base matrices for inconsistancies
      ! First check wether any matrix contains a closed loop
      test_matrix(:,:) = cmx_micro_hfp_srf(:,:)
      CALL check_flow_matrix(               &
        & nconti,                           &
        & test_matrix(:,:),                 &
        & has_loop                          &
          )
      IF (has_loop) CALL finish(TRIM(routine), '*** Circular connections detected: cmx_micro_hfp_srf')

      test_matrix(:,:) = cmx_micro_hfp_blg(:,:)
      CALL check_flow_matrix(               &
        & nconti,                           &
        & test_matrix(:,:),                 &
        & has_loop                          &
          )
      IF (has_loop) CALL finish(TRIM(routine), '*** Circular connections detected: cmx_micro_hfp_blg')

      test_matrix(:,:) = cmx_macro_hfp_srf(:,:)
      CALL check_flow_matrix(               &
        & nconti,                           &
        & test_matrix(:,:),                 &
        & has_loop                          &
          )
      IF (has_loop) CALL finish(TRIM(routine), '*** Circular connections detected: cmx_macro_hfp_srf')

      test_matrix(:,:) = cmx_macro_hfp_blg(:,:)
      CALL check_flow_matrix(               &
        & nconti,                           &
        & test_matrix(:,:),                 &
        & has_loop                          &
          )
      IF (has_loop) CALL finish(TRIM(routine), '*** Circular connections detected: cmx_macro_hfp_blg')

      ! Second check wether cmx_micro_ctc is consistent with current assumptions on geometric relations
      test_matrix(:,:) = 0._wp
      WHERE(cmx_micro_ctc(:,:) > 0.0)
        test_matrix(:,:) = 1.0_wp
      END WHERE
      DO it = 1, nconti
        IF (SUM(cmx_micro_ctc(it,:)) > 1._wp) THEN
          CALL finish(TRIM(routine), '*** Connected fraction exceeds 1.')
        ELSE IF (sum(cmx_micro_ctc(it,:)) > 0._wp .AND. dsl4jsb_Config(LEX_)%l_nested) THEN
          IF (SUM(test_matrix(:,it)) > 1._wp) THEN
            CALL finish(TRIM(routine), '*** Contained tile encloses more than one tile')
          ELSE IF (sum(cmx_micro_ctc(:,it)) > 0._wp .AND. SUM(cmx_micro_ctc(:,it)) < 1._wp) THEN
            CALL finish(TRIM(routine), '*** Contained tile encloses another tile only partially')
          ELSE IF (SUM(cmx_micro_ctc(:,it)) == 1._wp) THEN
            src_id=MAXLOC(cmx_micro_ctc(:,it), 1)
            IF (SUM(cmx_micro_ctc(:,src_id)) > 0._wp) THEN
              CALL finish(TRIM(routine), '*** 3-level tile-nesting not implemented')
            END IF
          END IF
        END IF
      END DO

      CALL tile%Get_fraction(fract=fract(:,:))

      fract_micro_hfp_srf(:,:) = 0._wp
      fract_micro_hfp_blg(:,:) = 0._wp
      fract_macro_hfp_srf(:,:) = 0._wp
      fract_macro_hfp_blg(:,:) = 0._wp

      ! First Loop over tiles to set up field containing all cover fractions and
      ! to determine the areas and cover fractions of any tile enclosed by annother tile.
      ! The latter are needed for nested constellations.
      ! Note that at the moment this only works for a setup in which a single surface
      ! feature SF1 is fully enclosed by surface feature SF2, which in turn is enclosed by
      ! another (other) surface feature (s).
      ALLOCATE(fracts(nproma, nconti, nblks))
      ALLOCATE(area_nclsd(nproma, nconti, nblks))
      ALLOCATE(fract_nclsd(nproma, nconti, nblks))
      area_nclsd(:,:,:)   = 0._wp
      fract_nclsd(:,:,:)  = 0._wp
      DO it = 1, nsibl
        IF (it == 1) THEN
          sibl => tile%parent_tile%Get_first_child_tile()
        ELSE
          sibl => sibl%Get_next_sibling_tile()
        END IF
        lex_ids = sibl%lcts(1)%lex_id

        IF (lex_ids == 0 .OR. lex_ids /= lex_ids) CYCLE

        CALL sibl%Get_fraction(fract=fract_sibl(:,:))
        fracts(:,lex_ids,:) = fract_sibl(:,:)

        IF (dsl4jsb_Config(LEX_)%l_nested) THEN

          dsl4jsb_Get_memory_tile(HYDRO_, sibl)
          dsl4jsb_Get_var2D_onDomain_tile_name(HYDRO_, area, sibl)    ! in


          trgt_id = INT(FINDLOC(cmx_micro_ctc(lex_ids,:), 1._wp, 1))
          IF (trgt_id /= 0) THEN
            area_nclsd(:,trgt_id,:)  = area_sibl(:,:)
            fract_nclsd(:,trgt_id,:) = fract_sibl(:,:)
          END IF
        END IF
      END DO


      ! Second Loop over neighbour tiles to get the
      ! the geometry and the effective connectedness of a tile.
      ALLOCATE(hlp1(nproma, nblks))
      ALLOCATE(hlp2(nproma, nblks))
      DO it = 1, nsibl

        IF (it == 1) THEN
          sibl => tile%parent_tile%Get_first_child_tile()
        ELSE
          sibl => sibl%Get_next_sibling_tile()
        END IF
        lex_ids = sibl%lcts(1)%lex_id

        IF (lex_ids == 0 .OR. lex_ids /= lex_ids) CYCLE

        dsl4jsb_Get_memory_tile(HYDRO_, sibl)
        dsl4jsb_Get_var2D_onDomain_tile_name(HYDRO_, h_surf, sibl)  ! in
        dsl4jsb_Get_var2D_onDomain_tile_name(HYDRO_, area, sibl)    ! in

        CALL sibl%Get_fraction(fract=fract_sibl(:,:))

        ! Setting up fields of micro scale connectivities
        ! corresponding to a direct contact of tiles:
        ! One possibility is to use cmx_micro_ctc directly, the other is
        ! to assume a well-mixed distirbution of tiles, with cmx_micro_ctc
        ! being weighted by the cover fraction of the connected sibing
        ! relative to the cover fractions of all connected siblings.

        hlp1(:,:) = 0._wp
        hlp2(:,:) = 0._wp

        IF (dsl4jsb_Config(LEX_)%l_wellmixed) THEN
          DO is = 1, nsibl
            IF (is == 1) THEN
              dummy => tile%parent_tile%Get_first_child_tile()
            ELSE
              dummy => dummy%Get_next_sibling_tile()
            END IF

            lex_idd = dummy%lcts(1)%lex_id

            IF (lex_idd == 0 .OR. lex_idd /= lex_idd) CYCLE

            IF (cmx_micro_ctc(lex_idt,lex_idd) > 0._wp) THEN
              hlp1(:,:) = hlp1(:,:) + fracts(:,lex_idd,:) * cmx_micro_ctc(lex_idt,lex_idd)
            END IF

            IF (cmx_micro_ctc(lex_ids,lex_idd) > 0._wp) THEN
              hlp2(:,:) = hlp2(:,:) + fracts(:,lex_idd,:) * cmx_micro_ctc(lex_ids,lex_idd)
            END IF

          END DO
          WHERE (hlp1(:,:) > 0._wp)
            hlp1(:,:) = cmx_micro_ctc(lex_idt,lex_ids) * fract_sibl(:,:) &
              &       / hlp1(:,:) * SUM(cmx_micro_ctc(lex_idt,:))
          END WHERE

          WHERE (hlp2(:,:) > 0._wp)
            hlp2(:,:) = cmx_micro_ctc(lex_ids,lex_idt) * fract(:,:) &
              &       / hlp2(:,:) * SUM(cmx_micro_ctc(lex_ids,:))
          END WHERE
        ELSE
          hlp1(:,:) = cmx_micro_ctc(lex_idt,lex_ids)
          hlp2(:,:) = cmx_micro_ctc(lex_ids,lex_idt)
        END IF

        ! Determine the geometric relationship between two connected tiles, i.e.
        ! the offset in surface height (h_offset_ct), the characteric distance between the tile centres (length_scale_ct),
        ! the length of the interface between the two tiles (length_iface_ct) and the slope between the tiles (micro_slope_ct)
        CALL get_geometry(                                      &
          & hlp1(:,:), hlp2(:,:),                               &
          & fract(:,:), fract_sibl(:,:),                        &
          & fract_nclsd(:,lex_idt,:), fract_nclsd(:,lex_ids,:), &
          & area(:,:), area_sibl(:,:),                          &
          & area_nclsd(:,lex_idt,:), area_nclsd(:,lex_ids,:),   &
          & h_surf(:,:), h_surf_sibl(:,:),                      &
          & h_offset_ct(:,lex_ids,:),                           &
          & length_scale_ct(:,lex_ids,:),                       &
          & length_iface_ct(:,lex_ids,:),                       &
          & micro_slope_ct(:,lex_ids,:)                         &
          & )

        ! Check if h_offset_ct is consistant with assumed micro-scale connection.
        ! If not adjust the direction of flow:
        ! For hydrological flow paths, connections are only assumed to be present
        ! for positive differences in surface height (for negative offsets, cmx_micro_hfp_x is inverted
        ! and for a zero-offset, cmx_micro_x is set to 0).
        ! For direct conections, connectivities are reversed in case of negative differences in surface height,
        ! to simplify the call to subroutines in update_lex and the determination of the effective slopes

        microcon_ctc_ct(:,lex_ids,:)     = hlp1(:,:)
        microcon_hfp_srf_ct(:,lex_ids,:) = cmx_micro_hfp_srf(lex_idt,lex_ids)
        microcon_hfp_blg_ct(:,lex_ids,:) = cmx_micro_hfp_blg(lex_idt,lex_ids)
        macrocon_hfp_srf_ct(:,lex_ids,:) = cmx_macro_hfp_srf(lex_idt,lex_ids)
        macrocon_hfp_blg_ct(:,lex_ids,:) = cmx_macro_hfp_blg(lex_idt,lex_ids)

        ! Setting zero connectivity where height offsets negative (and zero in case of microcon_hfp_x_ct)
        IF (cmx_micro_hfp_srf(lex_idt,lex_ids) > 0.0_wp) THEN
          WHERE (h_offset_ct(:,lex_ids,:) <= 0.0_wp)
            microcon_hfp_srf_ct(:,lex_ids,:) = 0._wp
          END WHERE
        END IF
        IF (cmx_micro_hfp_blg(lex_idt,lex_ids) > 0.0_wp) THEN
          WHERE (h_offset_ct(:,lex_ids,:) <= 0.0_wp)
            microcon_hfp_blg_ct(:,lex_ids,:) = 0._wp
          END WHERE
        END IF
        IF (cmx_micro_ctc(lex_idt,lex_ids) > 0.0_wp) THEN
          WHERE (h_offset_ct(:,lex_ids,:) < 0.0_wp)
            microcon_ctc_ct(:,lex_ids,:) = 0._wp
          END WHERE
        END IF

        ! Setting the micro connectivity where initially the connection is asumed to be in other direction
        ! but is changed due to sign of height offsets
        IF (cmx_micro_hfp_srf(lex_ids,lex_idt) > 0.0_wp) THEN ! here indexing is sibl/tile
          WHERE (h_offset_ct(:,lex_ids,:) > 0.0_wp) ! h_offset tile/sibl is + so height offset sibl/tile is -
            microcon_hfp_srf_ct(:,lex_ids,:) = cmx_micro_hfp_srf(lex_ids,lex_idt)
          END WHERE
        END IF

        IF (cmx_micro_hfp_blg(lex_ids,lex_idt) > 0.0_wp) THEN ! here indexing is sibl/tile
          WHERE (h_offset_ct(:,lex_ids,:) > 0.0_wp) ! h_offset tile/sibl is + so height offset sibl/tile is -
            microcon_hfp_blg_ct(:,lex_ids,:) = cmx_micro_hfp_blg(lex_ids,lex_idt)
          END WHERE
        END IF

        IF (cmx_micro_ctc(lex_ids,lex_idt) > 0.0_wp) THEN ! here indexing is sibl/tile
          WHERE (h_offset_ct(:,lex_ids,:) > 0.0_wp) ! h_offset tile/sibl is + so height offset sibl/tile  is -
            microcon_ctc_ct(:,lex_ids,:) = hlp2(:,:)
          END WHERE
        END IF

       ! Summing up the fractions of a tile that are connected to siblings via hydrological flow paths
        fract_micro_hfp_srf(:,:) = fract_micro_hfp_srf(:,:) &
          &                      + fract_sibl(:,:) * microcon_hfp_srf_ct(:,lex_ids,:)
        fract_micro_hfp_blg(:,:) = fract_micro_hfp_blg(:,:) &
          &                      + fract_sibl(:,:) * microcon_hfp_blg_ct(:,lex_ids,:)
        fract_macro_hfp_srf(:,:) = fract_macro_hfp_srf(:,:) &
          &                      + fract_sibl(:,:) * macrocon_hfp_srf_ct(:,lex_ids,:)
        fract_macro_hfp_blg(:,:) = fract_macro_hfp_blg(:,:) &
          &                      + fract_sibl(:,:) * macrocon_hfp_blg_ct(:,lex_ids,:)
      END DO

      ! Ensuring that sum of microcon_ctc_ct of a given tile still <= 1. after connections have been adjusted
      ! Here microcon_ctc is recalculated based on length_iface_ct calculated above

      ALLOCATE(hlp3(nproma, nblks))
      ALLOCATE(hlp4(nproma, nblks))
      hlp3(:,:) = SUM(length_iface_ct(:,:,:), DIM=2, MASK=microcon_ctc_ct(:,:,:) > 0._wp)
      hlp4(:,:) = MIN(1._wp, SUM(microcon_ctc_ct(:,:,:), DIM=2))
      microcon_ctc_ct(:,:,:) = 0._wp
      DO it = 1, nconti
        WHERE (hlp3(:,:) > 0._wp)
          microcon_ctc_ct(:,it,:) = length_iface_ct(:,it,:) / hlp3(:,:) * hlp4(:,:)
        END WHERE
      END DO

      !$ACC UPDATE DEVICE(h_offset_ct, length_scale_ct, length_iface_ct, micro_slope_ct) ASYNC(1)
      !$ACC UPDATE DEVICE(microcon_ctc_ct, microcon_hfp_srf_ct, microcon_hfp_blg_ct, macrocon_hfp_srf_ct, macrocon_hfp_blg_ct) ASYNC(1)

      ALLOCATE(fact_lag_micro_srf(nproma, nblks))
      ALLOCATE(fact_lag_micro_blg(nproma, nblks))
      ALLOCATE(fact_lag_macro_srf(nproma, nblks))
      ALLOCATE(fact_lag_macro_blg(nproma, nblks))

      ALLOCATE(fact_lag_srf_ct(nproma, nconti, nblks))
      ALLOCATE(fact_lag_blg_ct(nproma, nconti, nblks))

      ! Third loop over neighbour tiles:
      ! Get effective flow lags of tile and get effective slope
      slope_eff(:,:) = macro_slope(:,:)
      DO it = 1, nsibl
        IF (it == 1) THEN
          sibl => tile%parent_tile%Get_first_child_tile()
        ELSE
          sibl => sibl%Get_next_sibling_tile()
        END IF
        lex_ids = sibl%lcts(1)%lex_id

        IF (lex_ids == 0 .OR. lex_ids /= lex_ids) CYCLE

        CALL sibl%Get_fraction(fract=fract_sibl(:,:))

        ! Get effective flow lag between tile and sibling
        CALL get_flow_lags(                                     &
          & lex_idt, lex_ids,                                   &
          & vel_lat_srf, vel_lat_blg,                           &
          & dsl4jsb_Config(HYDRO_)%ret_macro_srf,               &
          & dsl4jsb_Config(HYDRO_)%ret_macro_blg,               &
          & length_scale_ct(:,lex_ids,:),                       &
          & fract(:,:), fract_sibl(:,:),                        &
          & microcon_hfp_srf_ct(:,lex_ids,:),                   &
          & microcon_hfp_blg_ct(:,lex_ids,:),                   &
          & macrocon_hfp_srf_ct(:,lex_ids,:),                   &
          & macrocon_hfp_blg_ct(:,lex_ids,:),                   &
          & fract_micro_hfp_srf(:,:), fract_micro_hfp_blg(:,:), &
          & fract_macro_hfp_srf(:,:), fract_macro_hfp_blg(:,:), &
          & fact_lag_micro_srf(:,:), fact_lag_micro_blg(:,:),   &
          & fact_lag_macro_srf(:,:), fact_lag_macro_blg(:,:)    &
          )

        fact_lag_srf_ct(:,lex_ids,:) = fact_lag_micro_srf(:,:) + fact_lag_macro_srf(:,:)
        fact_lag_blg_ct(:,lex_ids,:) = fact_lag_micro_blg(:,:) + fact_lag_macro_blg(:,:)

        ! Switching of individual flow components
        IF (.NOT. dsl4jsb_Config(LEX_)%l_macrtrns) THEN
          fact_lag_srf_ct(:,lex_ids,:) = fact_lag_srf_ct(:,lex_ids,:) - fact_lag_macro_srf(:,:)
          fact_lag_blg_ct(:,lex_ids,:) = fact_lag_blg_ct(:,lex_ids,:) - fact_lag_macro_blg(:,:)
          fact_lag_srf_ct(:,lex_ids,:) = MAX(0._wp, fact_lag_srf_ct(:,lex_ids,:))
          fact_lag_blg_ct(:,lex_ids,:) = MAX(0._wp, fact_lag_blg_ct(:,lex_ids,:))
        END IF
        IF (.NOT. dsl4jsb_Config(LEX_)%l_micrtrns) THEN
          fact_lag_srf_ct(:,lex_ids,:) = fact_lag_srf_ct(:,lex_ids,:) - fact_lag_micro_srf(:,:)
          fact_lag_blg_ct(:,lex_ids,:) = fact_lag_blg_ct(:,lex_ids,:) - fact_lag_micro_blg(:,:)
          fact_lag_srf_ct(:,lex_ids,:) = MAX(0._wp, fact_lag_srf_ct(:,lex_ids,:))
          fact_lag_blg_ct(:,lex_ids,:) = MAX(0._wp, fact_lag_blg_ct(:,lex_ids,:))
        END IF


        ! Determine effective lag factor as the sum of all (weighted) flow lag factors for all siblings
        fact_lag_srf(:,:) = fact_lag_srf(:,:) + fact_lag_srf_ct(:,lex_ids,:)
        fact_lag_blg(:,:) = fact_lag_blg(:,:) + fact_lag_blg_ct(:,lex_ids,:)


        ! Determine effective slope based on macro slope and weighted sum of micro slopes
        WHERE (micro_slope_ct(:,lex_ids,:) > 0.0_wp)
          slope_eff(:,:) = slope_eff(:,:)  &
            &            + micro_slope_ct(:,lex_ids,:) * microcon_ctc_ct(:,lex_ids,:)
        END WHERE
      END DO

      !$ACC UPDATE DEVICE(slope_eff, fact_lag_srf, fact_lag_blg) ASYNC(1)

      ! Fourth loop over neighbour tiles:
      ! Determines the fraction of flow from tile to sibling
      DO it = 1, nsibl
        IF (it == 1) THEN
          sibl => tile%parent_tile%Get_first_child_tile()
        ELSE
          sibl => sibl%Get_next_sibling_tile()
        END IF
        lex_ids = sibl%lcts(1)%lex_id

        IF (lex_ids == 0 .OR. lex_ids /= lex_ids) CYCLE

        ! Determine the fraction of flow from tile to sibling as:
        ! Ratio of flow lag between tile and sibling and effective lag factor for tile
        WHERE((fract_macro_hfp_srf(:,:) == 0._wp) .AND. (fract_micro_hfp_srf(:,:) == 0._wp) .OR. &
          & ABS(fact_lag_srf(:,:)) < EPSILON(1._wp))
          fract_latflow_srf_ct(:,lex_ids,:) = 0.0
        ELSE WHERE
          fract_latflow_srf_ct(:,lex_ids,:) = fact_lag_srf_ct(:,lex_ids,:) / fact_lag_srf(:,:)
        END WHERE
        WHERE((fract_macro_hfp_blg(:,:) == 0._wp) .AND. (fract_micro_hfp_blg(:,:) == 0._wp) .OR. &
          & ABS(fact_lag_blg(:,:)) < EPSILON(1._wp))
          fract_latflow_blg_ct(:,lex_ids,:) = 0.0
        ELSE WHERE
          fract_latflow_blg_ct(:,lex_ids,:) = fact_lag_blg_ct(:,lex_ids,:) / fact_lag_blg(:,:)
        END WHERE

        ! If a tile only partially connected to downstream tiles (cmx_micro_hfp + cmx_macro_hfp < 1.0)
        ! A fraction of the generated outflow is given to streamflow this is done by scaling
        ! the fractions lateral flow directed to downstream a downstream tile
        IF (dsl4jsb_Config(LEX_)%l_prtlcon) THEN
          WHERE   (SUM(microcon_hfp_srf_ct(:,:,:), DIM=2) + SUM(macrocon_hfp_srf_ct(:,:,:), DIM=2) < 1._wp)
            fract_latflow_srf_ct(:,lex_ids,:) = fract_latflow_srf_ct(:,lex_ids,:)  &
              & * (SUM(microcon_hfp_srf_ct(:,:,:), DIM=2) + SUM(macrocon_hfp_srf_ct(:,:,:), DIM=2))
          END WHERE
          WHERE   (SUM(microcon_hfp_blg_ct(:,:,:), DIM=2) + SUM(macrocon_hfp_blg_ct(:,:,:), DIM=2) < 1._wp)
            fract_latflow_blg_ct(:,lex_ids,:) =  fract_latflow_blg_ct(:,lex_ids,:) &
              & * (SUM(microcon_hfp_blg_ct(:,:,:), DIM=2) + SUM(macrocon_hfp_blg_ct(:,:,:), DIM=2))
          END WHERE
        END IF
      END DO

      !$ACC UPDATE DEVICE(fract_latflow_srf_ct, fract_latflow_blg_ct) ASYNC(1)

      DEALLOCATE(fract)
      DEALLOCATE(fract_sibl)
      DEALLOCATE(fracts)
      DEALLOCATE(hlp1)
      DEALLOCATE(hlp2)
      DEALLOCATE(hlp3)
      DEALLOCATE(hlp4)
      DEALLOCATE(test_matrix)
      DEALLOCATE(cmx_micro_ctc)
      DEALLOCATE(cmx_micro_hfp_srf)
      DEALLOCATE(cmx_micro_hfp_blg)
      DEALLOCATE(cmx_macro_hfp_srf)
      DEALLOCATE(cmx_macro_hfp_blg)
      DEALLOCATE(fract_micro_hfp_srf)
      DEALLOCATE(fract_micro_hfp_blg)
      DEALLOCATE(fract_macro_hfp_srf)
      DEALLOCATE(fract_macro_hfp_blg)
      DEALLOCATE(area_nclsd)
      DEALLOCATE(fract_nclsd)
      DEALLOCATE(fact_lag_micro_srf)
      DEALLOCATE(fact_lag_micro_blg)
      DEALLOCATE(fact_lag_macro_srf)
      DEALLOCATE(fact_lag_macro_blg)
      DEALLOCATE(fact_lag_srf_ct)
      DEALLOCATE(fact_lag_blg_ct)
    END IF

  END SUBROUTINE lex_init

#endif
END MODULE mo_lex_init
