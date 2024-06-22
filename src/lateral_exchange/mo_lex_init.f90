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

  dsl4jsb_Use_processes LEX_, HYDRO_

  dsl4jsb_Use_config(LEX_)
  dsl4jsb_Use_memory(LEX_)

  dsl4jsb_Use_config(HYDRO_)
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

    USE mo_lex_process,        ONLY: get_geometry, get_flow_lags
    USE mo_hydro_constants,    ONLY: ret_micro_srf, ret_macro_srf, &
                                   & ret_micro_blg, ret_macro_blg

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile ! Tile with LEX_ process running

    ! Local variables
    !
    dsl4jsb_Def_memory(LEX_)
    dsl4jsb_Def_memory(HYDRO_)
    dsl4jsb_Def_memory_tile(HYDRO_, sibl)

    CHARACTER(len=*), PARAMETER :: routine = modname//':lex_init'

    CLASS(t_jsb_tile_abstract), POINTER :: sibl ! Sibling tile with LEX_ process running

    TYPE(t_jsb_model), POINTER :: model
    TYPE(t_jsb_grid),  POINTER :: grid
    TYPE(t_jsb_vgrid), POINTER :: connecttiles

    INTEGER :: nproma,   &  !< nproma
      &        nblks,    &  !< nblocks
      &        nsibl,    &  !< number of siblings on current tile level
      &        nconti,   &  !< number of connected tiles on current tile level
      &        it,       &  !< counter for loop over tiles
      &        lex_idt,  &  !< lex id of of tile
      &        lex_ids      !< lex id of sibling

    REAL(wp), POINTER :: micro_connect(:,:), & !< base matrix for micro-scale connectivity between connected tiles
      &                  macro_connect(:,:)    !< base matrix for macro-scale connectivity between connected tiles
   
    ! Pointers to variables in memory for tile
    dsl4jsb_Real2D_onDomain :: &
      & h_surf,                &
      & area,                  &
      & macro_slope,           &
      & slope_eff,             &
      & lag_fact_srf,          &
      & lag_fact_blg

    dsl4jsb_Real3D_onDomain :: &
      & h_offset_ct,           &
      & length_scale_ct,       &
      & length_iface_ct,       &
      & micro_slope_ct,        &
      & fract_latflow_srf_ct,  &
      & fract_latflow_blg_ct,  &
      & microcon_slope_ct,     &
      & microcon_dffsn_ct,     &
      & macrocon_used_ct

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
      & lag_micro_srf,         &
      & lag_micro_blg,         &
      & lag_macro_srf,         &
      & lag_macro_blg

    dsl4jsb_Real3D_onDomain :: &
      & lag_srf_ct,            &
      & lag_blg_ct

    lex_idt = tile%lcts(1)%lex_id
    IF (lex_idt > 0) THEN

      dsl4jsb_Get_memory(LEX_)
      dsl4jsb_Get_memory(HYDRO_)

      dsl4jsb_Get_var2D_onDomain(HYDRO_, h_surf)               ! in
      dsl4jsb_Get_var2D_onDomain(HYDRO_, area)                 ! in
      dsl4jsb_Get_var2D_onDomain(HYDRO_, macro_slope)          ! in

      dsl4jsb_Get_var2D_onDomain(LEX_,   lag_fact_srf)         ! out
      dsl4jsb_Get_var2D_onDomain(LEX_,   lag_fact_blg)         ! out
      dsl4jsb_Get_var2D_onDomain(LEX_,   slope_eff)            ! out

      dsl4jsb_Get_var3D_onDomain(LEX_,   h_offset_ct)          ! out
      dsl4jsb_Get_var3D_onDomain(LEX_,   length_scale_ct)      ! out
      dsl4jsb_Get_var3D_onDomain(LEX_,   length_iface_ct)      ! out
      dsl4jsb_Get_var3D_onDomain(LEX_,   micro_slope_ct)       ! out
      dsl4jsb_Get_var3D_onDomain(LEX_,   microcon_slope_ct)    ! out
      dsl4jsb_Get_var3D_onDomain(LEX_,   microcon_dffsn_ct)    ! out
      dsl4jsb_Get_var3D_onDomain(LEX_,   macrocon_used_ct)     ! out
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
      ALLOCATE(micro_connect(nconti,nconti))
      ALLOCATE(macro_connect(nconti,nconti))


      ! Base connection matrixes: Give the fraction of fluxes to from tile to connected downstream tiles.
      ! These connectivties are not directly used in the model, but adjusted based on their
      ! consitance with the topographic information provided. 
     
      IF  (TRIM(model%config%usecase) == 'jsbach_topolex') THEN
        micro_connect=reshape( (/ 0.0, 1.0, 0.0 ,0.0,    &
                                  0.0, 0.0, 0.0 ,0.0,    &
                                  0.0, 0.0, 0.0 ,1.0,    &
                                  0.0, 0.0, 0.0 ,0.0 /), &
                      shape(micro_connect), order=(/2,1/) )


        macro_connect=reshape( (/ 0.0, 0.0, 0.0 ,1.0,    &
                                  0.0, 0.0, 0.0 ,1.0,    &
                                  0.0, 0.0, 0.0 ,0.0,    &
                                  0.0, 0.0, 0.0 ,0.0 /), &
                      shape(micro_connect), order=(/2,1/) )
      ELSE IF (TRIM(model%config%usecase) == 'jsbach_2tilelex') THEN
        micro_connect=reshape( (/ 0.0, 1.0,    &
                                  0.0, 0.0 /), &
                      shape(micro_connect), order=(/2,1/) )


        macro_connect=reshape( (/ 0.0, 1.0,    &
                                  0.0 ,0.0 /), &
                      shape(micro_connect), order=(/2,1/) )
      ELSE
        CALL finish(TRIM(routine), '*** Lateral connection matrix not defined for usecase')
      END IF



      CALL tile%Get_fraction(fract=fract(:,:))
      
      fract_microcon(:,:) = 0._wp
      fract_macrocon(:,:) = 0._wp

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
        CALL get_geometry(                                  &
          & micro_connect(lex_idt,lex_ids),                 &
          & macro_connect(lex_idt,lex_ids),                 &
          & fract(:,:), fract_sibl(:,:),                    &
          & h_surf(:,:), h_surf_sibl(:,:),                  &
          & area(:,:), area_sibl(:,:),                      &
          & h_offset_ct(:,lex_ids,:),                       &
          & length_scale_ct(:,lex_ids,:),                   &
          & length_iface_ct(:,lex_ids,:),                   &
          & micro_slope_ct(:,lex_ids,:)                     &
          & )        
        
        ! Check if h_offset_ct is consistant with assumed micro-scale connection.
        ! If not adjust the direction of flow:
        ! For gravitationally driven fluxes, connections are only assumed to be present
        ! for positive differences in surface height (for negative offsets, micro_connect is inverted and
        ! for a zero-offset, micro_connect is set to 0). 
        ! For diffusive fluxes, connections are reversed in case of negative differences in surface height,
        ! to simplify the call to subroutines in update_lex

        microcon_slope_ct(:,lex_ids,:) = micro_connect(lex_idt,lex_ids)
        microcon_dffsn_ct(:,lex_ids,:) = micro_connect(lex_idt,lex_ids)
        macrocon_used_ct(:,lex_ids,:)  = macro_connect(lex_idt,lex_ids)

        IF (micro_connect(lex_idt,lex_ids) > 0.0_wp) THEN
          WHERE (h_offset_ct(:,lex_ids,:) <= 0.0_wp)
            microcon_slope_ct(:,lex_ids,:) = 0._wp
          END WHERE
          WHERE (h_offset_ct(:,lex_ids,:) < 0.0_wp)
            microcon_dffsn_ct(:,lex_ids,:) = 0._wp
          END WHERE
        END IF
        IF (micro_connect(lex_ids,lex_idt) > 0.0_wp) THEN
          WHERE (h_offset_ct(:,lex_ids,:) > 0.0_wp)
            microcon_slope_ct(:,lex_ids,:) = micro_connect(lex_ids,lex_idt)
          END WHERE
          WHERE (h_offset_ct(:,lex_ids,:) > 0.0_wp)
            microcon_dffsn_ct(:,lex_ids,:) = micro_connect(lex_ids,lex_idt)
          END WHERE
        END IF

        ! Summing up the fractions of a tile that are connected to siblings 
        fract_microcon(:,:) = fract_microcon(:,:) &
                          & + fract_sibl(:,:) * microcon_slope_ct(:,lex_ids,:)
        fract_macrocon(:,:) = fract_macrocon(:,:) &
                          & + fract_sibl(:,:) * macrocon_used_ct(:,lex_ids,:)

      END DO

      !$ACC UPDATE DEVICE(h_offset_ct, length_scale_ct, length_iface_ct, micro_slope_ct) ASYNC(1)
      !$ACC UPDATE DEVICE(microcon_slope_ct, microcon_dffsn_ct, macrocon_used_ct) ASYNC(1)

      ! Second loop over neighbour tiles
      ! Get effective flow lags of tile and get effective slope

      ALLOCATE(lag_micro_srf(nproma, nblks))
      ALLOCATE(lag_micro_blg(nproma, nblks))
      ALLOCATE(lag_macro_srf(nproma, nblks))
      ALLOCATE(lag_macro_blg(nproma, nblks))

      ALLOCATE(lag_srf_ct(nproma, nconti, nblks))
      ALLOCATE(lag_blg_ct(nproma, nconti, nblks))
      
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
        CALL get_flow_lags(                                 &
          & lex_idt, lex_ids,                               &
          & ret_micro_srf, ret_micro_blg,                   &
          & ret_macro_srf, ret_macro_blg,                   &
          & microcon_slope_ct(:,lex_ids,:),                 &
          & macrocon_used_ct(:,lex_ids,:),                  &
          & fract(:,:), fract_sibl(:,:),                    &
          & fract_microcon(:,:), fract_macrocon(:,:),       &
          & lag_micro_srf(:,:), lag_micro_blg(:,:),         &
          & lag_macro_srf(:,:), lag_macro_blg(:,:)          &
          )

        lag_srf_ct(:,lex_ids,:) = lag_micro_srf(:,:) + lag_macro_srf(:,:)
        lag_blg_ct(:,lex_ids,:) = lag_micro_blg(:,:) + lag_macro_blg(:,:)


        ! Determine effective lag factor as the sum of all (weighted) flow lags for tile
        lag_fact_srf(:,:) = lag_fact_srf(:,:) + lag_srf_ct(:,lex_ids,:)
        lag_fact_blg(:,:) = lag_fact_blg(:,:) + lag_blg_ct(:,lex_ids,:)

        ! Determine effective slope based in macro slope and weighted sum of micro slopes        
        WHERE ((fract(:,:)+fract_microcon(:,:) > 0.0_wp) .AND.              &
            &  (micro_slope_ct(:,lex_ids,:)    > 0.0_wp))
          slope_eff(:,:) = slope_eff(:,:)                                   &
            &  + microcon_slope_ct(:,lex_ids,:)*micro_slope_ct(:,lex_ids,:) &
            &  * fract_sibl(:,:)/(fract(:,:)+fract_microcon(:,:))
        END WHERE 

      END DO

      !$ACC UPDATE DEVICE(slope_eff, lag_fact_srf, lag_fact_blg) ASYNC(1)

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
          fract_latflow_srf_ct(:,lex_ids,:) = lag_srf_ct(:,lex_ids,:) / lag_fact_srf(:,:)
          fract_latflow_blg_ct(:,lex_ids,:) = lag_blg_ct(:,lex_ids,:) / lag_fact_blg(:,:)
        END WHERE

      END DO

      !$ACC UPDATE DEVICE(fract_latflow_srf_ct, fract_latflow_blg_ct) ASYNC(1)

      DEALLOCATE(fract)
      DEALLOCATE(fract_sibl)
      DEALLOCATE(fract_microcon)
      DEALLOCATE(fract_macrocon)
      DEALLOCATE(lag_micro_srf)
      DEALLOCATE(lag_micro_blg)
      DEALLOCATE(lag_macro_srf)
      DEALLOCATE(lag_macro_blg)
      DEALLOCATE(micro_connect)
      DEALLOCATE(macro_connect)

      DEALLOCATE(lag_srf_ct)
      DEALLOCATE(lag_blg_ct)
    END IF

  END SUBROUTINE lex_init

#endif
END MODULE mo_lex_init
