!>
!! @brief Interface between HD and the ocean, through a coupler
!!
!! @author Veronika Gayler (MPI-M)
!!
!! @par Revision History
!! First version by Veronika Gayler (2013-06-29)

#ifndef __NO_JSBACH__
!----------------------------
#include "omp_definitions.inc"
!----------------------------

MODULE mo_interface_hd_ocean

  USE mo_jsb_model_class,     ONLY: t_jsb_model
  USE mo_jsb_tile_class,      ONLY: t_jsb_tile
  USE mo_jsb_grid_class,      ONLY: t_jsb_grid
  USE mo_jsb_base,            ONLY: get_model
  USE mo_kind,                ONLY: wp
  USE mo_coupling_config,     ONLY: is_coupled_run
  USE mo_exception,           ONLY: warning
#ifdef YAC_coupling
  USE mo_yac_finterface,      ONLY: yac_fput, yac_fget_field_id_from_name
#endif

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: interface_hd_ocean
  CHARACTER(len=*), PARAMETER :: modname = 'mo_interface_hd_ocean'

CONTAINS

  SUBROUTINE interface_hd_ocean(model_id)

    INTEGER, INTENT(in) :: model_id

    ! Local variables
    !
    TYPE(t_jsb_model), POINTER    :: model
    TYPE(t_jsb_tile),  POINTER    :: tile
    TYPE(t_jsb_grid),  POINTER    :: grid
    CHARACTER(len=*),  PARAMETER  :: routine = TRIM(modname)//':interface_hd_ocean'

    ! Pointers to variables in memory
    REAL(wp), POINTER :: &
      & discharge_ocean   (:,:)

    ! Local variables

    LOGICAL               :: write_coupler_restart
    INTEGER               :: field_id
    INTEGER               :: nbr_hor_cells
    REAL(wp), ALLOCATABLE :: buffer(:,:)
    INTEGER               :: collection_size, nbr_pointsets, nbr_subdomains
    INTEGER               :: info, ierror !< return values form cpl_put/get calls


    IF ( .NOT. is_coupled_run() ) RETURN

    model => get_model(model_id)
    tile => model%land
    grid => model%grid

    ! Get reference to variables in memory
    discharge_ocean      => tile%mem%soil %water_flux %discharge_ocean   (:,:)

    ! Send river discharge
    ! --------------------

    write_coupler_restart = .FALSE.

    ALLOCATE(buffer(grid%nproma * grid%nblks, 1))
    buffer(:,:) = 0.0_wp
!vg     DO iblk = 1, grid%nblks
!vg       nn = (iblk - 1) * grid%nproma
!vg       DO n = 1, nproma
!vg         buffer(nn+n,1) = discharge_ocean(n,iblk)
!vg       END DO
!vg     END DO
    nbr_hor_cells = model%grid%ntotal
    buffer(:,:) = RESHAPE(discharge_ocean(:,:), (/ nbr_hor_cells, 1 /))

#ifdef YAC_coupling
    CALL yac_fget_field_id_from_name('river_runoff', field_id)
#endif
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
