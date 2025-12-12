!> Initialization of the lateral exchange process.
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
!>#### Initialisation of connections between tiles for lateral exchange process
!>
!> This module contains the subroutine lex_init which initializes the lateral exchange
!> with respect to geometric properties and relationships between connected tiles.
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
  !>#### Intialize lateral exchange process
  !>
  !> The subroutine loops over all tiles where the LEX process is active and determines
  !> the geometric properties and relationships between connected tiles. Furthermore,
  !> it computes the effective flow lags between connected tiles used later for the lateral
  !> flow calculations.
  !>
  SUBROUTINE lex_init(tile)

    USE mo_lex_process,        ONLY: check_flow_matrix, get_geometry, &
      &                              get_flow_lags
    USE mo_hydro_constants,    ONLY: vel_lat_srf, vel_lat_blg

    USE mo_lex_constants,      ONLY: &
      & cmx_micro_ctc_2tilelex, cmx_micro_hfp_srf_2tilelex, cmx_micro_hfp_blg_2tilelex, &
      & cmx_macro_hfp_srf_2tilelex, cmx_macro_hfp_blg_2tilelex,                         &
      & cmx_micro_ctc_3tilelex, cmx_micro_hfp_srf_3tilelex, cmx_micro_hfp_blg_3tilelex, &
      & cmx_macro_hfp_srf_3tilelex, cmx_macro_hfp_blg_3tilelex,                         &
      & cmx_micro_ctc_topolex,  cmx_micro_hfp_srf_topolex, cmx_micro_hfp_blg_topolex,   &
      & cmx_macro_hfp_srf_topolex, cmx_macro_hfp_blg_topolex

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile    !< current tile

    dsl4jsb_Def_config(LEX_)
    dsl4jsb_Def_config(HYDRO_)

    dsl4jsb_Def_memory(LEX_)
    dsl4jsb_Def_memory(HYDRO_)
    dsl4jsb_Def_memory_tile(HYDRO_, sibl)

    CLASS(t_jsb_tile_abstract), POINTER :: sibl, &  !< Sibling tile
      &                                    first    !< First tile
    TYPE(t_jsb_model), POINTER :: model             !< model instance
    TYPE(t_jsb_grid),  POINTER :: grid              !< model grid
    TYPE(t_jsb_vgrid), POINTER :: connecttiles      !< grid for connected tiles

    LOGICAL :: has_loop     !< T: Flow matrix contains closed loop

    INTEGER ::                   &
      &  nproma,                 &  !< nproma
      &  nblks,                  &  !< nblocks
      &  nsibl,                  &  !< number of siblings on current tile level
      &  nconti,                 &  !< number of connected tiles on current tile level
      &  it,                     &  !< counter for loop over tiles
      &  is,                     &  !< counter for loop over tiles when within loop over tiles
      &  ic,                     &  !< counter for loop over cells
      &  ib,                     &  !< counter for loop over blocks
      &  src_id,                 &  !< source lex_id
      &  trgt_id,                &  !< target lex_id
      &  lex_idt,                &  !< lex id of of tile
      &  lex_ids,                &  !< lex id of sibling
      &  lex_idf                    !< lex id of the first tile on current tile level

    REAL(wp), POINTER ::         &
      &  test_matrix(:,:),       &  !< matrix used to check for inconsistances in base matrices
      &  cmx_micro_ctc(:,:),     &  !< base matrix for micro-scale connectivity via direct contact
      &  cmx_micro_hfp_srf(:,:), &  !< base matrix for micro-scale connectivity via surface hydrological flow paths
      &  cmx_micro_hfp_blg(:,:), &  !< base matrix for micro-scale connectivity via below-ground hydrological flow paths
      &  cmx_macro_hfp_srf(:,:), &  !< base matrix for macro-scale connectivity via surface hydrological flow paths
      &  cmx_macro_hfp_blg(:,:)     !< base matrix for macro-scale connectivity via below-ground hydrological flow paths

    CHARACTER(len=*), PARAMETER :: routine = modname//':lex_init' !< name of this routine

    dsl4jsb_Real2D_onDomain :: &
      & fract,                 &    !< fraction of tile []
      & char_h_surf,           &    !< vertical offset of tile surface [m]
      & char_area,             &    !< characteristic average size of tile feature (e.g. depressions) [m2]
      & macro_slope,           &    !< macro scale slope of tile []
      & slope_eff,             &    !< effective slope (micro- and macro-scale) used for lateral flow calculations []
      & fact_lag_srf,          &    !< normalized lag for lateral surface fluxes [s-1]
      & fact_lag_blg                !< normalized lag for lateral below-ground fluxes [s-1]

    dsl4jsb_Real3D_onDomain :: &
      & h_offset_ct,           &    !< difference in surface elevation between adjacent tiles [m]
      & length_scale_ct,       &    !< char. length scale of surface heterogeneity [m]
      & length_iface_ct,       &    !< char. interface length between adjacent tiles [m m-2]
      & micro_slope_ct,        &    !< micro-scale slope between connected tiles []
      & fract_latflow_srf_ct,  &    !< fraction of lateral surface flow to adjacent tile []
      & fract_latflow_blg_ct,  &    !< fraction of lateral below-ground flow to adjacent tile []
      & microcon_ctc_ct,       &    !< micro-scale connectivity factors for direct contact []
      & microcon_hfp_srf_ct,   &    !< micro-scale connectivity factors for surface hydrological flow paths []
      & microcon_hfp_blg_ct,   &    !< micro-scale connectivity factors for below-ground hydrological flow paths []
      & macrocon_hfp_srf_ct,   &    !< macro-scale connectivity factors for surface hydrological flow paths []
      & macrocon_hfp_blg_ct         !< macro-scale connectivity factors for below-ground hydrological flow paths []

    dsl4jsb_Real2D_onDomain :: &
      & fract_sibl,            &    !< fraction of sibling tile []
      & char_h_surf_sibl,      &    !< vertical offset of sibling tile surface [m]
      & char_area_sibl              !< characteristic average size of sibling tile feature (e.g. depressions) [m2]

    dsl4jsb_Real2D_onDomain :: &
      & hlp1,                  &    !< helper array 1
      & hlp2,                  &    !< helper array 2
      & hlp3,                  &    !< helper array 3
      & hlp4,                  &    !< helper array 4
      & fract_micro_hfp_srf,   &    !< weighted sum of micro-scale connect. fractions (surface) on tile []
      & fract_micro_hfp_blg,   &    !< weighted sum of micro-scale connect. fractions (below-ground) on tile []
      & fract_macro_hfp_srf,   &    !< weighted sum of macro-scale connect. fractions (surface) on tile []
      & fract_macro_hfp_blg,   &    !< weighted sum of macro-scale connect. fractions (below-ground) on tile []
      & fact_lag_micro_srf,    &    !< lag factor for micro-scale surface fluxes between tile and sibling [s-1]
      & fact_lag_micro_blg,    &    !< lag factor for micro-scale below-ground fluxes between tile and sibling [s-1]
      & fact_lag_macro_srf,    &    !< lag factor for macro-scale surface fluxes between tile and sibling [s-1]
      & fact_lag_macro_blg          !< lag factor for macro-scale below-ground fluxes between tile and sibling [s-1]

    dsl4jsb_Real3D_onDomain :: &
      & fracts,                &    !< cover fractions of all tiles []
      & area_nclsd,            &    !< characteristic area enclosed by tile [m2]
      & fract_nclsd,           &    !< cover fraction enclosed by tile [-]
      & fact_lag_srf_ct,       &    !< lag fact for micro+macro scale surface fluxes between tile and sibling [s-1]
      & fact_lag_blg_ct             !< lag fact for micro+macro scale below-ground fluxes between tile and sibling [s-1]

    IF (.NOT. tile%Is_process_active(LEX_)) RETURN

    IF (debug_on()) CALL message(routine, 'Initializing variables for lateral exchange on tile ' // &
      & TRIM(tile%name) // ' (lex_id=' // TRIM(int2string(tile%lcts(1)%lex_id)) // ')')

    ! Get LEX ID of current tile
    lex_idt = tile%lcts(1)%lex_id
    IF (lex_idt > 0) THEN

      dsl4jsb_Get_memory(LEX_)
      dsl4jsb_Get_memory(HYDRO_)

      dsl4jsb_Get_var2D_onDomain(HYDRO_, char_h_surf)          ! in
      dsl4jsb_Get_var2D_onDomain(HYDRO_, char_area)            ! in
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

      !> Retrieve the general connection matrixes for the actual usecase
      IF  (TRIM(model%config%usecase) == 'jsbach_topolex') THEN
        cmx_micro_ctc     = reshape(cmx_micro_ctc_topolex,     shape(cmx_micro_ctc),     order=(/2,1/))
        cmx_micro_hfp_srf = reshape(cmx_micro_hfp_srf_topolex, shape(cmx_micro_hfp_srf), order=(/2,1/))
        cmx_macro_hfp_srf = reshape(cmx_macro_hfp_srf_topolex, shape(cmx_macro_hfp_srf), order=(/2,1/))
        cmx_micro_hfp_blg = reshape(cmx_micro_hfp_blg_topolex, shape(cmx_micro_hfp_blg), order=(/2,1/))
        cmx_macro_hfp_blg = reshape(cmx_macro_hfp_blg_topolex, shape(cmx_macro_hfp_blg), order=(/2,1/))
      ELSE
        CALL finish(TRIM(routine), '*** Lateral connection matrix not defined for usecase')
      END IF

      !> Check the base matrices for inconsistancies
      ! First check whether any matrix contains a closed loop
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

      ! Second check whether cmx_micro_ctc is consistent with current assumptions on geometric relations
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

      !> First loop over tiles is used to set up arrays containing all cover fractions and
      !> to determine the areas and cover fractions of any tile enclosed by another tile.
      !> The latter are needed for nested constellations.
      !>
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
          dsl4jsb_Get_var2D_onDomain_tile_name(HYDRO_, char_area, sibl)    ! in


          trgt_id = INT(FINDLOC(cmx_micro_ctc(lex_ids,:), 1._wp, 1))
          IF (trgt_id /= 0) THEN
            area_nclsd(:,trgt_id,:)  = char_area_sibl(:,:)
            fract_nclsd(:,trgt_id,:) = fract_sibl(:,:)
          END IF
        END IF
      END DO


      !> The second loop over tiles derives the
      !> geometry and the effective connectedness of a tile.
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
        dsl4jsb_Get_var2D_onDomain_tile_name(HYDRO_, char_h_surf, sibl)  ! in
        dsl4jsb_Get_var2D_onDomain_tile_name(HYDRO_, char_area, sibl)    ! in

        CALL sibl%Get_fraction(fract=fract_sibl(:,:))

        ! Setting up arrays of micro-scale connectivities
        ! corresponding to a direct contact of tiles:
        ! One possibility is to use cmx_micro_ctc directly, the other is
        ! to assume a well-mixed distribution of tiles, with cmx_micro_ctc
        ! being weighted by the cover fraction of the connected sibling
        ! relative to the cover fractions of all connected siblings.

        hlp1(:,:) = 0._wp
        hlp2(:,:) = 0._wp

        IF (dsl4jsb_Config(LEX_)%l_wellmixed) THEN
          DO is = 1, nsibl
            IF (is == 1) THEN
              first => tile%parent_tile%Get_first_child_tile()
            ELSE
              first => first%Get_next_sibling_tile()
            END IF

            lex_idf = first%lcts(1)%lex_id

            IF (lex_idf == 0 .OR. lex_idf /= lex_idf) CYCLE

            IF (cmx_micro_ctc(lex_idt,lex_idf) > 0._wp) THEN
              hlp1(:,:) = hlp1(:,:) + fracts(:,lex_idf,:) * cmx_micro_ctc(lex_idt,lex_idf)
            END IF

            IF (cmx_micro_ctc(lex_ids,lex_idf) > 0._wp) THEN
              hlp2(:,:) = hlp2(:,:) + fracts(:,lex_idf,:) * cmx_micro_ctc(lex_ids,lex_idf)
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

        ! Determine the geometric relationship between two connected tiles, i.e. the offset
        ! in surface height (h_offset_ct), the characteric distance between the tile
        ! centres (length_scale_ct), the length of the interface between the two
        ! tiles (length_iface_ct), and the slope between the tiles (micro_slope_ct)
        DO ib = 1, nblks
          DO ic = 1, nproma
            CALL get_geometry(                                          &
              & hlp1(ic,ib), hlp2(ic,ib),                               &
              & fract(ic,ib), fract_sibl(ic,ib),                        &
              & fract_nclsd(ic,lex_idt,ib), fract_nclsd(ic,lex_ids,ib), &
              & char_area(ic,ib), char_area_sibl(ic,ib),                &
              & area_nclsd(ic,lex_idt,ib), area_nclsd(ic,lex_ids,ib),   &
              & char_h_surf(ic,ib), char_h_surf_sibl(ic,ib),            &
              & h_offset_ct(ic,lex_ids,ib),                             &
              & length_scale_ct(ic,lex_ids,ib),                         &
              & length_iface_ct(ic,lex_ids,ib),                         &
              & micro_slope_ct(ic,lex_ids,ib)                           &
              & )
          END DO
        END DO
        !$ACC UPDATE ASYNC(1) DEVICE(h_offset_ct, length_scale_ct, length_iface_ct, micro_slope_ct)

        ! Check if h_offset_ct is consistant with assumed micro-scale connection.
        ! If not adjust the direction of flow:
        ! For hydrological flow paths, connections are only assumed to be present
        ! for positive differences in surface height (for negative offsets, cmx_micro_hfp_x is inverted
        ! and for a zero-offset, cmx_micro_x is set to 0).
        ! For direct connections, connectivities are reversed in case of negative differences in surface height,
        ! to simplify the call to subroutines in update_lex and the determination of the effective slopes.

        microcon_ctc_ct(:,lex_ids,:)     = hlp1(:,:)
        microcon_hfp_srf_ct(:,lex_ids,:) = cmx_micro_hfp_srf(lex_idt,lex_ids)
        microcon_hfp_blg_ct(:,lex_ids,:) = cmx_micro_hfp_blg(lex_idt,lex_ids)
        macrocon_hfp_srf_ct(:,lex_ids,:) = cmx_macro_hfp_srf(lex_idt,lex_ids)
        macrocon_hfp_blg_ct(:,lex_ids,:) = cmx_macro_hfp_blg(lex_idt,lex_ids)

        ! Setting zero connectivity where height offsets are negative (and zero in case of microcon_hfp_x_ct)
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
          WHERE (h_offset_ct(:,lex_ids,:) > 0.0_wp) ! h_offset tile/sibl is + so height offset sibl/tile is -
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

      ! Ensuring that sum of microcon_ctc_ct of a given tile is still <= 1. after connections have been adjusted.
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

      !$ACC UPDATE ASYNC(1) &
      !$ACC   DEVICE(microcon_ctc_ct, microcon_hfp_srf_ct, microcon_hfp_blg_ct) &
      !$ACC   DEVICE(macrocon_hfp_srf_ct, macrocon_hfp_blg_ct)

      ALLOCATE(fact_lag_micro_srf(nproma, nblks))
      ALLOCATE(fact_lag_micro_blg(nproma, nblks))
      ALLOCATE(fact_lag_macro_srf(nproma, nblks))
      ALLOCATE(fact_lag_macro_blg(nproma, nblks))

      ALLOCATE(fact_lag_srf_ct(nproma, nconti, nblks))
      ALLOCATE(fact_lag_blg_ct(nproma, nconti, nblks))

      !> In the third loop over tiles their effective flow lags and effective slope is computed
      !>
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
        DO ib = 1, nblks
          DO ic = 1, nproma
            CALL get_flow_lags(                                         &
              & lex_idt, lex_ids,                                       &
              & vel_lat_srf, vel_lat_blg,                               &
              & dsl4jsb_Config(HYDRO_)%ret_macro_srf,                   &
              & dsl4jsb_Config(HYDRO_)%ret_macro_blg,                   &
              & length_scale_ct(ic,lex_ids,ib),                         &
              & fract(ic,ib), fract_sibl(ic,ib),                        &
              & microcon_hfp_srf_ct(ic,lex_ids,ib),                     &
              & microcon_hfp_blg_ct(ic,lex_ids,ib),                     &
              & macrocon_hfp_srf_ct(ic,lex_ids,ib),                     &
              & macrocon_hfp_blg_ct(ic,lex_ids,ib),                     &
              & fract_micro_hfp_srf(ic,ib), fract_micro_hfp_blg(ic,ib), &
              & fract_macro_hfp_srf(ic,ib), fract_macro_hfp_blg(ic,ib), &
              & fact_lag_micro_srf(ic,ib), fact_lag_micro_blg(ic,ib),   &
              & fact_lag_macro_srf(ic,ib), fact_lag_macro_blg(ic,ib)    &
              )
          END DO
        END DO

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

      !$ACC UPDATE ASYNC(1) DEVICE(slope_eff, fact_lag_srf, fact_lag_blg)

      !> In the fourth loop over the tiles, the ratio of flow from tile to sibling compared
      !> to the total lateral outflow from the tile is computed.
      !>
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

        ! If a tile is only partially connected to downstream tiles (cmx_micro_hfp + cmx_macro_hfp < 1.0),
        ! only a fraction of the generated outflow is given to streamflow.
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

      !$ACC UPDATE ASYNC(1) DEVICE(fract_latflow_srf_ct, fract_latflow_blg_ct)

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
