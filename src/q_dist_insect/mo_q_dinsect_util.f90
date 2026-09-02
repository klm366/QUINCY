!> Utility functions for insect disturbance module
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
MODULE mo_q_dinsect_util
#ifndef __NO_JSBACH__

  USE mo_kind,                    ONLY: wp
  USE mo_jsb_control,             ONLY: debug_on
  USE mo_exception,               ONLY: message
  USE mo_jsb_math_constants,      ONLY: zero, eps8
  USE mo_lnd_bgcm_idx
  USE mo_lnd_bgcm_store,          ONLY: t_lnd_bgcm_store
  USE mo_lnd_bgcm_store_class,    ONLY: DINSECT_BGCM_VEG_EXTRACT_ID, &
  & DINSECT_BGCM_FAST_LITTER_ID, &
  & DINSECT_BGCM_SEEDB_FASTL_ID, &
  & DINSECT_BGCM_SLOW_LITTER_ID, &
  & DINSECT_BGCM_DROPPING_ID

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: reset_dinsect_fluxes, Write2DMatrixToFile
  PUBLIC :: print_dinsect_info

  CHARACTER(len=*), PARAMETER :: modname = 'mo_q_dinsect_util'

CONTAINS

  ! ======================================================================================================= !
  !>reset dinsect fluxes (to zero)
  ! ----------------------------------------------------------------------------------------------------- !
  SUBROUTINE reset_dinsect_fluxes(tile, options)

    USE mo_jsb_class,             ONLY: Get_model
    USE mo_jsb_tile_class,        ONLY: t_jsb_tile_abstract
    USE mo_jsb_task_class,        ONLY: t_jsb_task_options
    USE mo_jsb_model_class,       ONLY: t_jsb_model
    USE mo_jsb_process_class,     ONLY: Q_DINSECT_
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Use_memory(Q_DINSECT_)
    ! ----------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile         !< one tile with data structure for one lct
    TYPE(t_jsb_task_options),   INTENT(in)    :: options      !< model options
    ! ----------------------------------------------------------------------------------------------------- !
    TYPE(t_jsb_model),      POINTER   :: model                !< the model
    TYPE(t_lnd_bgcm_store), POINTER   :: bgcm_store           !< the bgcm store of this tile
    INTEGER                           :: iblk, ics, ice, nc   !< dimensions

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':reset_dinsect_fluxes'
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Def_mt1L2D :: dinsect_dropping_mt
    dsl4jsb_Def_mt2L2D :: dinsect_veg_extract_mt
    dsl4jsb_Def_mt2L2D :: dinsect_fast_litter_mt
    dsl4jsb_Def_mt2L2D :: dinsect_slow_litter_mt
    dsl4jsb_Real2D_onChunk :: dinsect_pool_total_c          ! total C from insect pools [mol C/m2]
    dsl4jsb_Real2D_onChunk :: dinsect_total_mort_dt         ! fraction of mortality rate caused by insect disturbances [fraction/timestep]
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Def_memory(Q_DINSECT_)
    ! ----------------------------------------------------------------------------------------------------- !
    iblk    = options%iblk
    ics     = options%ics
    ice     = options%ice
    nc      = options%nc
    ! ----------------------------------------------------------------------------------------------------- !
    !@TODO: comment in when new tile module version is merged into QS?!
    IF (.NOT. tile%Is_process_calculated(Q_DINSECT_)) RETURN
    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')
    ! ----------------------------------------------------------------------------------------------------- !
    model  => Get_model(tile%owner_model_id)
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Get_memory(Q_DINSECT_)
    ! ----------------------------------------------------------------------------------------------------- !

    bgcm_store => tile%bgcm_store
    dsl4jsb_Get_mt2L2D(DINSECT_BGCM_VEG_EXTRACT_ID, dinsect_veg_extract_mt)
    dsl4jsb_Get_mt2L2D(DINSECT_BGCM_FAST_LITTER_ID, dinsect_fast_litter_mt)
    dsl4jsb_Get_mt2L2D(DINSECT_BGCM_SLOW_LITTER_ID, dinsect_slow_litter_mt)
    dsl4jsb_Get_mt1L2D(DINSECT_BGCM_DROPPING_ID, dinsect_dropping_mt)
    dsl4jsb_Get_var2D_onChunk(Q_DINSECT_, dinsect_pool_total_c) ! out
    dsl4jsb_Get_var2D_onChunk(Q_DINSECT_, dinsect_total_mort_dt) ! out

    ! ----------------------------------------------------------------------------------------------------- !
    !>1.0 bgcm fluxes
    !>
    ! set the bgcm matrices to zero, and write these values back to the bgcm below
    ! no need to get the values of the bgcm (`CALL bgcm_store%Get_matrix_from_bgcm_*`)
    !   as these values are zeroed anyway
    dinsect_veg_extract_mt(:,:,:) = zero
    dinsect_fast_litter_mt(:,:,:) = zero
    dinsect_slow_litter_mt(:,:,:) = zero
    dinsect_dropping_mt(:,:) = zero
    dinsect_pool_total_c = zero
    dinsect_total_mort_dt = zero

  END SUBROUTINE reset_dinsect_fluxes


  ! ======================================================================================================= !
  !> Write2DMatrixToFile
  !> Write the 2D matrix to a file
  !> Used for printing the damageMatrix for debug
  ! ----------------------------------------------------------------------------------------------------- !
  SUBROUTINE Write2DMatrixToFile(filename, matrix)
      IMPLICIT NONE
      CHARACTER(LEN=*), INTENT(IN) :: filename
      REAL(wp), DIMENSION(:,:), INTENT(IN) :: matrix
      INTEGER :: i, j, nrows, ncols, ios
      LOGICAL :: file_exists

      !Check if the file already exists
      INQUIRE(FILE=filename, EXIST=file_exists)
      IF (file_exists) THEN
           !Delete the file if it already exists using a system command
          !CALL SYSTEM('rm -f ' // TRIM(filename))
          ! So I need to open the file for appending
          OPEN(UNIT=10, FILE=filename, STATUS='OLD', ACTION='WRITE', POSITION='APPEND', IOSTAT=ios)
      ELSE
          OPEN(UNIT=10, FILE=filename, STATUS='NEW', ACTION='WRITE', IOSTAT=ios)
      END IF

      IF (ios /= 0) THEN
          PRINT *, 'Error opening file.'
          RETURN
      END IF



      ! Get the dimensions of the matrix
      nrows = SIZE(matrix, 1)
      ncols = SIZE(matrix, 2)

      ! Write the matrix to the file
      DO i = 1, ncols
          WRITE(10, '(I4, 1X, I3, 1X, F5.1, 1X, F12.8)') int(matrix(1, i)), int(matrix(2, i)), matrix(3, i), matrix(4, i)
      END DO

      ! Close the file
      CLOSE(UNIT=10)

  END SUBROUTINE Write2DMatrixToFile
! ======================================================================================================= !
! ======================================================================================================= !
  SUBROUTINE print_dinsect_info(type, prescribe, from_namelist, year, year_span)
    USE mo_exception,   ONLY: message
 ! ----------------------------------------------------------------------------------------------------- !
    INTEGER, INTENT(IN) :: type, year, year_span
    LOGICAL, INTENT(IN) :: prescribe, from_namelist
    CHARACTER(LEN=10) :: year_str, year_span_str
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':print_dinsect_info'
 ! ----------------------------------------------------------------------------------------------------- !

    ! Convert integers to strings
    WRITE(year_str, '(I10)') year
    WRITE(year_span_str, '(I10)') year_span

    ! For testbed screen output
    PRINT *, "============================================="
    PRINT *, "QUINCY Insect Disturbance Module (DINSECT)"
    PRINT *, "============================================="
    PRINT *, "Module Configuration:"
    PRINT *, "  - Insect Type: ", MERGE("Bark Beetle", "Defoliator ", type == 1)
    PRINT *, "  - Damage Mode: ", MERGE("Prescribed ", "Dynamic    ", prescribe)

    IF (from_namelist) THEN
        PRINT *, "  - Data Source: Namelist parameters"
    ELSE
        PRINT *, "  - Data Source: Default values"
    END IF

    PRINT *, "  - Start Year: ", year
    PRINT *, "  - Duration (years): ", year_span
    PRINT *, "============================================="
    PRINT *, "DINSECT module initialized"

    ! For normal log
    CALL message(routine, '=============================================')
    CALL message(routine, 'QUINCY Insect Disturbance Module (DINSECT)')
    CALL message(routine, '=============================================')
    CALL message(routine, 'Module Configuration:')
    CALL message(routine, '  - Insect Type: ' // TRIM(MERGE("Bark Beetle", "Defoliator ", type == 1)))
    CALL message(routine, '  - Damage Mode: ' // TRIM(MERGE("Prescribed ", "Dynamic    ", prescribe)))

    IF (from_namelist) THEN
        CALL message(routine, '  - Data Source: Namelist parameters')
    ELSE
        CALL message(routine, '  - Data Source: Default values')
    END IF

    CALL message(routine, '  - Start Year: ' // TRIM(ADJUSTL(year_str)))
    CALL message(routine, '  - Duration (years): ' // TRIM(ADJUSTL(year_span_str)))
    CALL message(routine, '=============================================')
    CALL message(routine, 'DINSECT module initialized')


  END SUBROUTINE print_dinsect_info
! ======================================================================================================= !

#endif
END MODULE mo_q_dinsect_util
