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
!>
!>#### Initialisation and reading of lex (lateral [intertile] exchange) variables.
!>
MODULE mo_lex_init
#ifndef __NO_JSBACH__

  USE mo_kind,               ONLY: wp
  USE mo_exception,          ONLY: message, finish
  USE mo_jsb_control,        ONLY: debug_on

  USE mo_jsb_model_class,    ONLY: t_jsb_model
  USE mo_jsb_grid_class,     ONLY: t_jsb_grid, t_jsb_vgrid
  USE mo_jsb_grid,           ONLY: Get_grid, Get_vgrid
  USE mo_jsb_tile_class,     ONLY: t_jsb_tile_abstract
  USE mo_jsb_class,          ONLY: get_model
  USE mo_jsb_io_netcdf,      ONLY: t_input_file, jsb_netcdf_open_input
  USE mo_jsb_io,             ONLY: missval
  USE mo_jsb_impl_constants, ONLY: ifs_nsoil, ifs_soil_depth
  USE mo_util,               ONLY: int2string, soil_init_from_texture

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
    USE mo_hydro_constants,    ONLY: vel_lat_srf, vel_lat_blg,        &
      &                              ret_micro_srf, ret_macro_srf,    &
      &                              ret_micro_blg, ret_macro_blg

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile ! Tile with LEX_ process running


    dsl4jsb_Def_config(LEX_)

    dsl4jsb_Def_memory(LEX_)
    dsl4jsb_Def_memory(HYDRO_)
    dsl4jsb_Def_memory_tile(HYDRO_, sibl)

    CHARACTER(len=*), PARAMETER :: routine = modname//':lex_init'
    CHARACTER(len=5000)         :: message_text_long

    CLASS(t_jsb_tile_abstract), POINTER :: sibl ! Sibling tile with LEX_ process running

    TYPE(t_jsb_model), POINTER :: model
    TYPE(t_jsb_grid),  POINTER :: grid
    TYPE(t_jsb_vgrid), POINTER :: connecttiles

    LOGICAL :: has_loop     !< .TRUE. Flow matrix contains closed loop

    INTEGER :: nproma,   &  !< nproma
      &        nblks,    &  !< nblocks
      &        nsibl,    &  !< number of siblings on current tile level
      &        nconti,   &  !< number of connected tiles on current tile level
      &        it,       &  !< counter for loop over tiles
      &        is,       &  !< counter for loop over siblings (when on loop over tiles)
      &        cntr,     &  !< counter used for checkin for circular tile conections
      &        src_id,   &  !< source lex_id
      &        trgt_id,  &  !< target lex_id
      &        lex_idt,  &  !< lex id of of tile
      &        lex_ids      !< lex id of sibling

    REAL(wp), POINTER :: test_matrix(:,:),     & !< matrix used to check for inconsitances in base matrices
      &                  micro_cnct_grvt(:,:), & !< base matrix for micro-scale connectivity between connected tiles (slopes; used for runoff and drainage)
      &                  micro_cnct_dfsn(:,:), & !< base matrix for micro-scale connectivity between connected tiles (diffusion)
      &                  macro_cnct_grvt(:,:)    !< base matrix for macro-scale connectivity between connected tiles
   
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
      & microcon_grvt_ct,      &
      & microcon_dfsn_ct,      &
      & macrocon_grvt_ct

    ! Pointers to variables in memory for sibling
    dsl4jsb_Real2D_onDomain :: &
      & h_surf_sibl,           &
      & area_sibl

    ! Locally allocated vectors
    dsl4jsb_Real2D_onDomain :: &
      & fract,                 &
      & fract_sibl,            &
      & fract_microcon,        &
      & fract_macrocon,        &
      & fact_lag_micro_srf,    &
      & fact_lag_micro_blg,    &
      & fact_lag_macro_srf,    &
      & fact_lag_macro_blg

    dsl4jsb_Real3D_onDomain :: &
      & area_nclsd,            &
      & fract_nclsd,           &
      & fact_lag_srf_ct,       &
      & fact_lag_blg_ct


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
      dsl4jsb_Get_var3D_onDomain(LEX_,   microcon_grvt_ct)     ! out
      dsl4jsb_Get_var3D_onDomain(LEX_,   microcon_dfsn_ct)     ! out
      dsl4jsb_Get_var3D_onDomain(LEX_,   macrocon_grvt_ct)     ! out
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
      ALLOCATE(fract_microcon(nproma, nblks))
      ALLOCATE(fract_macrocon(nproma, nblks))
      ALLOCATE(test_matrix(nconti,nconti))
      ALLOCATE(micro_cnct_grvt(nconti,nconti))
      ALLOCATE(micro_cnct_dfsn(nconti,nconti))
      ALLOCATE(macro_cnct_grvt(nconti,nconti))

      ! Get lateral exchange config
      dsl4jsb_Get_config(LEX_)


      ! Base connection matrixes: 
      ! In case of micro_cnct_grvt and macro_cnct_grvt, they give the fraction of fluxes (runoff and drainage) 
      ! from a given tile to a specific, connected (downstream) tile. In case of micro_cnct_dfsn they indicate
      ! the fraction of the outer border of a tile that is assumed to be connected to a specific tile. 
      ! Note that connectivties are not directly used in the model, but adjusted based on their consitance
      ! with the topographic information provided. 
      !
      ! micro_cnct_grvt: Connectivity for gravity driven fluxes (surface runoff & lateral drainage)
      !                  between two interfacing tiles. Here, micro scale connections reflect landscapes where
      !                  the tiles are well mixed and have comparativly small characteristic surface areas. 
      !                  such as ponds and small scale depressions. 
      !                  The matrix elements need to be provided to indicate the fractions of fluxes originating
      !                  on the tile given by the row to the tile indicated by the column.
      !
      ! micro_cnct_dfsn: Connectivity for diffusive fluxes.
      !                  between two interfacing tiles. Here, micro scale connections reflect landscapes where
      !                  the tiles are well mixed and have comparativly small characteristic surface areas. 
      !                  such as ponds and small scale depressions. 
      !                  The matrix elements need to be provided to so that entries are given for the tile with the 
      !                  relativly small characteristic surface area (given by the row)  
      !                  which is assumed to be partially encompassed by a larger tile (indicated by the column)
      !
      ! macro_cnct_grvt: Connectivity for gravity driven fluxes (surface runoff & lateral drainage)
      !                  between two tiles that are spatially seperated across larger distances, such as uplands, slopes and lowlands.
      !                  The matrix elements need to be provided to so that entries are given for the tile with the 
      !                  The matrix elements need to be provided to indicate the fractions of fluxes originating
      !                  on the tile given by the row to the tile indicated by the column.
      !
      ! For micro_cnct_grvt and macro_cnct_grvt, indiviual matrix entries should not exceed 1. 
      ! However, it is intended that the sum of a row can be larger than 1. The flow lags assumed for surface runoof and drainage
      ! are linearly dependend on the matrix elements of micro_cnct_grvt and macro_cnct_grvt. Hence, the respective matrix entries 
      ! determine the overall retention time hence the outgoing fluxes. But the distribution amongst all receiving tiles
      ! is done relative to the weighted sum of all entries on a row, so the latter only depends on the values relatve to each other.
      ! If, however, the sum of all entries on a row is smaller than 1, a fraction (1 - SUM(row elements)) of the outgoing fluxes
      ! is directily given to the streamflow.
      ! For micro_cnct_dfsn the sum of entries on a row should not exeed 1. Furthermore, the model does not (yet) properly considder nested 
      ! configurations. These options for such a setup are limited to a case in which a single surface feature SF1 is fully encompassed by a second surface feature SF2, which in turn is encompassed by another (other) surface feature(s).

      IF  (TRIM(model%config%usecase) == 'jsbach_topolex') THEN

        micro_cnct_grvt=reshape( (/ 0.0, 1.0, 0.0 ,0.0,    &
                                    0.0, 0.0, 0.0 ,0.0,    &
                                    0.0, 0.0, 0.0 ,1.0,    &
                                    0.0, 0.0, 0.0 ,0.0 /), &
                       shape(micro_cnct_grvt), order=(/2,1/) )

        micro_cnct_dfsn=reshape( (/ 0.0, 0.0, 0.0 ,0.0,    &
                                    1.0, 0.0, 0.0 ,0.0,    &
                                    0.0, 0.0, 0.0 ,0.0,    &
                                    0.0, 0.0, 1.0 ,0.0 /), &
                       shape(micro_cnct_grvt), order=(/2,1/) )

        macro_cnct_grvt=reshape( (/ 0.0, 0.0, 0.0 ,1.0,    &
                                    0.0, 0.0, 0.0 ,1.0,    &
                                    0.0, 0.0, 0.0 ,0.0,    &
                                    0.0, 0.0, 0.0 ,0.0 /), &
                       shape(micro_cnct_grvt), order=(/2,1/) )
      ELSE IF (TRIM(model%config%usecase) == 'jsbach_2tilelex') THEN

        micro_cnct_grvt=reshape( (/ 0.0, 1.0,    &
                                    0.0, 0.0 /), &
                       shape(micro_cnct_grvt), order=(/2,1/) )

        micro_cnct_dfsn=reshape( (/ 0.0, 0.0,    &
                                    1.0, 0.0 /), &
                       shape(micro_cnct_grvt), order=(/2,1/) )

        macro_cnct_grvt=reshape( (/ 0.0, 0.0,    &
                                    0.0 ,0.0 /), &
                       shape(micro_cnct_grvt), order=(/2,1/) )

      ELSE
        CALL finish(TRIM(routine), '*** Lateral connection matrix not defined for usecase')
      END IF

      ! Check the base matrices for inconsistancies
      ! First check wether any matrix contains a closed loop
      test_matrix(:,:) = micro_cnct_grvt(:,:)
      CALL check_flow_matrix(               &    
        & nconti,                           &
        & test_matrix(:,:),                 &
        & has_loop                          &
          )
      IF (has_loop) CALL finish(TRIM(routine), '*** Circular connections detected: micro_cnct_grvt')

      test_matrix(:,:) = macro_cnct_grvt(:,:)
      CALL check_flow_matrix(               &    
        & nconti,                           &
        & test_matrix(:,:),                 &
        & has_loop                          &
          )
      IF (has_loop) CALL finish(TRIM(routine), '*** Circular connections detected: macro_cnct_grvt')


      ! Second check wether micro_cnct_dfsn is consistent with current assumptions on geometric relations
      test_matrix(:,:) = 0._wp 
      WHERE(micro_cnct_dfsn(:,:) > 0.0)
        test_matrix(:,:) = 1.0_wp
      END WHERE
      DO it = 1, nconti
        IF (SUM(micro_cnct_dfsn(it,:)) > 1._wp) THEN
          CALL finish(TRIM(routine), '*** Connected fraction exceeds 1.')             
        ELSE IF (sum(micro_cnct_dfsn(it,:)) > 0._wp .AND. dsl4jsb_Config(LEX_)%l_nested) THEN
          IF (SUM(test_matrix(:,it)) > 1._wp) THEN
            CALL finish(TRIM(routine), '*** Contained tile encloses more than one tile') 
          ELSE IF (sum(micro_cnct_dfsn(:,it)) > 0._wp .AND. SUM(micro_cnct_dfsn(:,it)) < 1._wp) THEN
            CALL finish(TRIM(routine), '*** Contained tile encloses another tile only partially')
          ELSE IF (SUM(micro_cnct_dfsn(:,it)) == 1._wp) THEN
            src_id=MAXLOC(micro_cnct_dfsn(:,it), 1)
            IF (SUM(micro_cnct_dfsn(:,src_id)) > 0._wp) THEN
              CALL finish(TRIM(routine), '*** 3-level tile-nesting not implemented')            
            END IF
          END IF
        END IF     
      END DO

      CALL tile%Get_fraction(fract=fract(:,:))
      
      fract_microcon(:,:) = 0._wp
      fract_macrocon(:,:) = 0._wp

      ! First Loop over tiles to determine the areas and cover fractions 
      ! of any tile enclosed by annother tile. Needed for nested constellations
      ! note that at the moment this only works for a setup in which a single surface 
      ! feature SF1 is fully enclosed by surface feature SF2, which in turn is enclosed by 
      ! another surface features.
      ALLOCATE(area_nclsd(nproma, nconti, nblks))
      ALLOCATE(fract_nclsd(nproma, nconti, nblks))
      area_nclsd(:,:,:)   = 0._wp
      fract_nclsd(:,:,:)  = 0._wp
      IF (dsl4jsb_Config(LEX_)%l_nested) THEN
        DO it = 1, nsibl
 
          IF (it == 1) THEN
            sibl => tile%parent_tile%Get_first_child_tile()
          ELSE
            sibl => sibl%Get_next_sibling_tile()
          END IF
          lex_ids = sibl%lcts(1)%lex_id

          IF (lex_ids == 0 .OR. lex_ids /= lex_ids) CYCLE  

          dsl4jsb_Get_memory_tile(HYDRO_, sibl)
          dsl4jsb_Get_var2D_onDomain_tile_name(HYDRO_, area, sibl)    ! in

          ALLOCATE(fract_sibl(nproma, nblks))
          CALL sibl%Get_fraction(fract=fract_sibl(:,:))


          trgt_id = INT(FINDLOC(micro_cnct_dfsn(lex_ids,:), 1._wp, 1))
          IF (trgt_id /= 0) THEN
            area_nclsd(:,trgt_id,:)  = area_sibl(:,:)
            fract_nclsd(:,trgt_id,:) = fract_sibl(:,:)
          END IF
          DEALLOCATE(fract_sibl)
        END DO
      END IF


      ! First Loop over neighbour tiles
      ! Get geometry and the effective connectedness of tile
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

        ALLOCATE(fract_sibl(nproma, nblks))
        CALL sibl%Get_fraction(fract=fract_sibl(:,:))

        ! Determine the geometric relationship between two connected tiles, i.e.
        ! the offset in surface height (h_offset_ct), the characteric distance between the tile centres (length_scale_ct),
        ! the length of the interface between the two tiles (length_iface_ct) and the slope between the tiles (micro_slope_ct)
        CALL get_geometry(                                      &
          & micro_cnct_grvt(lex_idt,lex_ids),                   &
          & micro_cnct_dfsn(lex_idt,lex_ids),                   &
          & micro_cnct_dfsn(lex_ids,lex_idt),                   &
          & macro_cnct_grvt(lex_idt,lex_ids),                   &
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
        ! For gravity driven fluxes, connections are only assumed to be present
        ! for positive differences in surface height (for negative offsets, micro_cnct_grvt is inverted
        ! and for a zero-offset, micro_cnct_grvt is set to 0). 
        ! For diffusive fluxes, connections are reversed in case of negative differences in surface height,
        ! to simplify the call to subroutines in update_lex and the determination of the effective slopes

        microcon_grvt_ct(:,lex_ids,:) = micro_cnct_grvt(lex_idt,lex_ids)
        microcon_dfsn_ct(:,lex_ids,:) = micro_cnct_dfsn(lex_idt,lex_ids)
        macrocon_grvt_ct(:,lex_ids,:) = macro_cnct_grvt(lex_idt,lex_ids)

        ! Setting zero connectivity where height offsets negative (and zero in case of microcon_grvt_ct)
        IF (micro_cnct_grvt(lex_idt,lex_ids) > 0.0_wp) THEN
          WHERE (h_offset_ct(:,lex_ids,:) <= 0.0_wp)
            microcon_grvt_ct(:,lex_ids,:) = 0._wp
          END WHERE
          WHERE (h_offset_ct(:,lex_ids,:) < 0.0_wp)
            microcon_dfsn_ct(:,lex_ids,:) = 0._wp
          END WHERE
        END IF
        ! Setting the micro connectivity where initially the connection is asumed to be in other direction
        ! but is changed due to sign of height offsets 
        IF (micro_cnct_grvt(lex_ids,lex_idt) > 0.0_wp) THEN ! here indexing is sibl/tile
          WHERE (h_offset_ct(:,lex_ids,:) > 0.0_wp) ! h_offset tile/sibl is + so height offset sibl/tile -
            microcon_grvt_ct(:,lex_ids,:) = micro_cnct_grvt(lex_ids,lex_idt)
          END WHERE
          WHERE (h_offset_ct(:,lex_ids,:) > 0.0_wp)
            microcon_dfsn_ct(:,lex_ids,:) = micro_cnct_dfsn(lex_ids,lex_idt)
          END WHERE
        END IF

       ! Summing up the fractions of a tile that are connected to siblings via slopes
        fract_microcon(:,:) = fract_microcon(:,:) &
                          & + fract_sibl(:,:) * microcon_grvt_ct(:,lex_ids,:)
        fract_macrocon(:,:) = fract_macrocon(:,:) &
                          & + fract_sibl(:,:) * macrocon_grvt_ct(:,lex_ids,:)

      END DO
      
      ! Ensuring that sum of microcon_dfsn_ct of a given tile still <= 1. after connections have been adjusted
      DO it = 1, nconti
        WHERE (SUM(microcon_dfsn_ct(:,:,:), DIM=2) > 1._wp)
          microcon_dfsn_ct(:,it,:) = microcon_dfsn_ct(:,it,:) / SUM(microcon_dfsn_ct(:,:,:), DIM=2)
        END WHERE
      END DO

      !$ACC UPDATE DEVICE(h_offset_ct, length_scale_ct, length_iface_ct, micro_slope_ct) ASYNC(1)
      !$ACC UPDATE DEVICE(microcon_grvt_ct, microcon_dfsn_ct, macrocon_grvt_ct) ASYNC(1)

      ! Second loop over neighbour tiles
      ! Get effective flow lags of tile and get effective slope

      ALLOCATE(fact_lag_micro_srf(nproma, nblks))
      ALLOCATE(fact_lag_micro_blg(nproma, nblks))
      ALLOCATE(fact_lag_macro_srf(nproma, nblks))
      ALLOCATE(fact_lag_macro_blg(nproma, nblks))

      ALLOCATE(fact_lag_srf_ct(nproma, nconti, nblks))
      ALLOCATE(fact_lag_blg_ct(nproma, nconti, nblks))
      
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
        CALL get_flow_lags(                                    &
          & lex_idt, lex_ids,                                  &
          & vel_lat_srf, vel_lat_blg,                          &
          & ret_micro_srf, ret_micro_blg,                      &
          & ret_macro_srf, ret_macro_blg,                      &
          & length_scale_ct(:,lex_ids,:),                      &
          & microcon_grvt_ct(:,lex_ids,:),                     &
          & macrocon_grvt_ct(:,lex_ids,:),                     &
          & fract(:,:), fract_sibl(:,:),                       &
          & fract_microcon(:,:), fract_macrocon(:,:),          &
          & fact_lag_micro_srf(:,:), fact_lag_micro_blg(:,:),  &
          & fact_lag_macro_srf(:,:), fact_lag_macro_blg(:,:)   &
          )

        fact_lag_srf_ct(:,lex_ids,:) = fact_lag_micro_srf(:,:) + fact_lag_macro_srf(:,:)
        fact_lag_blg_ct(:,lex_ids,:) = fact_lag_micro_blg(:,:) + fact_lag_macro_blg(:,:)


        ! Determine effective lag factor as the sum of all (weighted) flow lag factors for all siblings
        fact_lag_srf(:,:) = fact_lag_srf(:,:) + fact_lag_srf_ct(:,lex_ids,:)
        fact_lag_blg(:,:) = fact_lag_blg(:,:) + fact_lag_blg_ct(:,lex_ids,:)

        ! Determine effective slope based on macro slope and weighted sum of micro slopes
        WHERE (micro_slope_ct(:,lex_ids,:) > 0.0_wp)
          slope_eff(:,:) = slope_eff(:,:)  &
            &            + micro_slope_ct(:,lex_ids,:) * microcon_dfsn_ct(:,lex_ids,:)
        END WHERE 

      END DO

      !$ACC UPDATE DEVICE(slope_eff, fact_lag_srf, fact_lag_blg) ASYNC(1)

      ! Third loop over neighbour tiles
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
        WHERE((fract_macrocon(:,:) == 0._wp) .AND. (fract_microcon(:,:) == 0._wp))
          fract_latflow_srf_ct(:,lex_ids,:) = 0.0
          fract_latflow_blg_ct(:,lex_ids,:) = 0.0  
        ELSE WHERE
          fract_latflow_srf_ct(:,lex_ids,:) = fact_lag_srf_ct(:,lex_ids,:) / fact_lag_srf(:,:)
          fract_latflow_blg_ct(:,lex_ids,:) = fact_lag_blg_ct(:,lex_ids,:) / fact_lag_blg(:,:)
        END WHERE

        ! If a tile is not fully connected to downstream tiles (micro_cnct_grvt + macro_cnct_grvt > 1.0)
        ! A fraction of the generated outflow is given to streamflow this is done by scaling 
        ! the fractions lateral flow directed to downstream a downstream tile        
        WHERE   (SUM(microcon_grvt_ct(:,:,:), DIM=2) + SUM(macrocon_grvt_ct(:,:,:), DIM=2) < 1._wp)
          fract_latflow_srf_ct(:,lex_ids,:) = fract_latflow_srf_ct(:,lex_ids,:)  &
            & * (SUM(microcon_grvt_ct(:,:,:), DIM=2) + SUM(macrocon_grvt_ct(:,:,:), DIM=2))
          fract_latflow_blg_ct(:,lex_ids,:) =  fract_latflow_blg_ct(:,lex_ids,:) &
            & * (SUM(microcon_grvt_ct(:,:,:), DIM=2) + SUM(macrocon_grvt_ct(:,:,:), DIM=2))
        END WHERE
       
      END DO

      !$ACC UPDATE DEVICE(fract_latflow_srf_ct, fract_latflow_blg_ct) ASYNC(1)

      DEALLOCATE(fract)
      DEALLOCATE(fract_sibl)
      DEALLOCATE(fract_microcon)
      DEALLOCATE(fract_macrocon)
      DEALLOCATE(test_matrix)
      DEALLOCATE(micro_cnct_grvt)
      DEALLOCATE(macro_cnct_grvt)
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
