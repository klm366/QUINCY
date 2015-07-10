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

  USE mo_jsb_model_class, ONLY: t_jsb_model
  USE mo_jsb_tile_class,  ONLY: t_jsb_tile
  USE mo_jsb_base,        ONLY: get_model
  USE mo_run_config,      ONLY: nlev, ltimer
  USE mo_kind                ,ONLY: wp
  USE mo_model_domain        ,ONLY: t_patch
  USE mo_echam_phy_memory    ,ONLY: prm_field
                                
  USE mo_parallel_config     ,ONLY: nproma
  USE mo_echam_sfc_indices   ,ONLY: iwtr, iice

  USE mo_sync                ,ONLY: SYNC_C, sync_patch_array
  USE mo_coupling_config,     ONLY: is_coupled_run
  USE mo_exception           ,ONLY: finish
#ifdef YAC_coupling
  USE mo_yac_finterface      ,ONLY: yac_fput, yac_fget, yac_fget_nbr_fields, yac_fget_field_ids
#else
  USE mo_icon_cpl_exchg      ,ONLY: ICON_cpl_put, ICON_cpl_get
  USE mo_icon_cpl_def_field  ,ONLY: ICON_cpl_get_nbr_fields, ICON_cpl_get_field_ids
  USE mo_icon_cpl_restart    ,ONLY: icon_cpl_write_restart
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
    CHARACTER(len=*),  PARAMETER  :: routine = TRIM(modname)//':interface_hd_ocean'

    ! Pointers to variables in memory
    REAL(wp), POINTER :: &
      & discharge_ocean   (:,:)

    ! Local variables

    LOGICAL               :: write_coupler_restart
    INTEGER               :: nbr_hor_points ! = inner and halo points
    INTEGER               :: field_shape(3)
    REAL(wp), ALLOCATABLE :: buffer(:,:)
    INTEGER               :: info, ierror !< return values form cpl_put/get calls


    IF ( .NOT. is_coupled_run() ) RETURN

    model => get_model(model_id)
    tile => model%land

    ! Get reference to variables in memory
    discharge_ocean      => tile%mem%soil %water_flux %discharge_ocean   (:,:)

    ! Send river discharge
    ! --------------------

    write_coupler_restart = .FALSE.

!vg     DO i_blk = 1, p_patch%nblk_c
!vg       nn = (i_blk - 1) * nproma
!vg       DO n = 1, nproma
!vg         buffer(nn+n,1) = discharge_ocean(n,i_blk)
!vg       END DO
!vg     END DO

#ifdef YAC_coupling

    !vg: field_id, nbr_hor_points, collection_size, nbr_pointsets, nbr_subdomains, send_field, info, ierror)
    CALL yac_fput ('DISCHARGE', nbr_hor_points, 1, 1, 1, buffer, info, ierror )
    IF ( info > 1 .AND. info < 7 ) THEN
       write_coupler_restart = .TRUE.
       WRITE ( 6 , * ) 'DISCHARGE', minval(buffer(1:nbr_hor_points,1)), maxval(buffer(1:nbr_hor_points,1))
    ELSE IF ( info == 7 ) THEN
       CALL warning "..."
    ENDIF
#else
    CALL finish(TRIM(routine), 'ERROR: the old ICON coupler is no longer supported. Use YAC.')
#endif

    IF ( write_coupler_restart ) THEN
#ifdef YAC_coupling
       CALL finish(TRIM(routine), 'ERROR: restart not supported yet for YAC')
#endif
    ENDIF

    DEALLOCATE(buffer)

  END SUBROUTINE interface_hd_ocean

#endif
END MODULE mo_interface_hd_ocean
