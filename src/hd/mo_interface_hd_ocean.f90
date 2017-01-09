!>
!! @brief Interface between HD and the ocean, through a coupler
!!
!! @author Veronika Gayler (MPI-M)
!!
!! @par Revision History
!! First version by Veronika Gayler (2013-06-29)

MODULE mo_interface_hd_ocean
#if !defined(__NO_JSBACH__) && !defined(__NO_JSBACH_HD__)

  USE mo_kind,                ONLY: wp
  USE mo_exception,           ONLY: warning

  USE mo_jsb_model_class,     ONLY: t_jsb_model
  USE mo_jsb_class,           ONLY: get_model
  USE mo_jsb_tile_class,      ONLY: t_jsb_tile
  USE mo_jsb_grid,            ONLY: Get_grid
  USE mo_jsb_grid_class,      ONLY: t_jsb_grid
  USE mo_coupling_config,     ONLY: is_coupled_run
#ifdef YAC_coupling
  USE mo_yac_finterface,      ONLY: yac_fdef_field, yac_fput, yac_fget_field_id_from_name
#endif

  dsl4jsb_Use_processes SOIL
  dsl4jsb_Use_memory(SOIL)

  IMPLICIT NONE
  PRIVATE

#ifdef YAC_coupling
  PUBLIC :: jsb_fdef_hd_fields
#endif
  PUBLIC :: interface_hd_ocean

  CHARACTER(len=*), PARAMETER :: modname = 'mo_interface_hd_ocean'

CONTAINS

#ifdef YAC_coupling
  SUBROUTINE jsb_fdef_hd_fields(comp_id, domain_id, cell_point_ids, cell_mask_ids)

    INTEGER, INTENT(in) :: comp_id, domain_id
    INTEGER, INTENT(in) :: cell_point_ids(:), cell_mask_ids(:)

    INTEGER :: field_id

    CALL yac_fdef_field (      &
      & 'river_runoff',        &
      & comp_id,               &
      & domain_id,             &
      & cell_point_ids,        &
      & cell_mask_ids,         &
      & 1,                     &
      & field_id )

  END SUBROUTINE jsb_fdef_hd_fields
#endif

  SUBROUTINE interface_hd_ocean(model_id)

    INTEGER, INTENT(in) :: model_id

    ! Local variables
    !
    TYPE(t_jsb_model), POINTER    :: model
    TYPE(t_jsb_tile),  POINTER    :: tile
    TYPE(t_jsb_grid),  POINTER    :: grid

    dsl4jsb_Def_memory(SOIL)

    CHARACTER(len=*),  PARAMETER  :: routine = TRIM(modname)//':interface_hd_ocean'

    ! Pointers to variables in memory
    dsl4jsb_Real2D_onDomain ::discharge_ocean

    ! Local variables

    LOGICAL               :: write_coupler_restart
    INTEGER               :: field_id
    INTEGER               :: nbr_hor_cells
    REAL(wp), ALLOCATABLE :: buffer(:,:)
    INTEGER               :: collection_size, nbr_pointsets, nbr_subdomains
    INTEGER               :: info, ierror !< return values form cpl_put/get calls
    INTEGER               :: iblk, n, nn, nlen


    IF ( .NOT. is_coupled_run() ) RETURN

    model => get_model(model_id)
    tile => model%Get_top_tile()
    grid => Get_grid(model%grid_id)

    ! Get reference to variables in memory
    dsl4jsb_Get_memory(SOIL)

    dsl4jsb_Get_var2D_onDomain(SOIL, discharge_ocean, IN)

    ! Send river discharge
    ! --------------------

    write_coupler_restart = .FALSE.

    ALLOCATE(buffer(grid%nproma * grid%nblks, 1))
    buffer(:,:) = 0.0_wp
    DO iblk = 1, grid%nblks
      nn = (iblk - 1) * grid%nproma
      IF (iblk /= grid%nblks) THEN
        nlen = grid%nproma
      ELSE
        nlen = grid%npromz
      END IF
      DO n = 1, grid%nproma
        buffer(nn+n,1) = discharge_ocean(n,iblk)
      END DO
    END DO

#ifdef YAC_coupling
    CALL yac_fget_field_id_from_name('river_runoff', field_id)
#endif
    nbr_hor_cells = grid%ntotal
    collection_size = 1
    nbr_pointsets   = 1
    nbr_subdomains  = 1
#ifdef YAC_coupling
    CALL yac_fput(field_id, nbr_hor_cells, collection_size, nbr_pointsets, nbr_subdomains, buffer, info, ierror)
#endif
    IF ( info > 1 .AND. info < 7 ) write_coupler_restart = .TRUE.

    IF ( write_coupler_restart ) THEN
       CALL warning(TRIM(routine), 'YAC says it is put for restart')
    ENDIF

    DEALLOCATE(buffer)

  END SUBROUTINE interface_hd_ocean

#endif
END MODULE mo_interface_hd_ocean
