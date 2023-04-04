!> Summary: Contains methods for collecting, storing and providing bgc materials
!>
!>### Contains methods (used on tiles) to ease working with bgc materials
!>
!> Copyright (c) 2023 Julia Nabel
!>
!>#### Revision History
!> First version and ongoing changes   by Julia Nabel (2023-01-12 - present)
!>
!>#### License
!>
!> This code is subject to the MPI-M-ICON - License - Agreement in it's most recent form.
!> Please see [this link](https://code.mpimet.mpg.de/attachments/download/20888/MPI-M-ICONLizenzvertragV2.6.pdf) for details.
!> Where software is supplied by third parties, it is indicated in the headers of the routines.
!>
MODULE mo_lnd_bgcm_store
#ifndef __NO_JSBACH__

  USE mo_exception,             ONLY: finish, message
  USE mo_kind,                  ONLY: wp
  USE mo_jsb_control,           ONLY: debug_on
  USE mo_jsb_parallel,          ONLY: Get_omp_thread
  USE mo_util,                  ONLY: int2string
  USE mo_lnd_bgcm_store_class,  ONLY: t_lnd_bgcm_store_abstract

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: t_lnd_bgcm_store

  
  ! ======================================================================================================= !
  !>
  !> bgc material store - stores pointers to the element variables of all bgcms across all processes
  !>
  TYPE, EXTENDS(t_lnd_bgcm_store_abstract) :: t_lnd_bgcm_store
  CONTAINS
    !> These deferred routines must match the interfaces provided in mo_lnd_bgcm_store_class
    PROCEDURE :: Init           => Init_bgcm_store
    PROCEDURE :: Append_bgcms   => Append_process_bgcms_to_store
    PROCEDURE :: Get_matrix_from_bgcm_2l_2d
    PROCEDURE :: Get_matrix_from_bgcm_2l_3d
    GENERIC :: Get_matrix_from_bgcm_2l => Get_matrix_from_bgcm_2l_2d, Get_matrix_from_bgcm_2l_3d
    PROCEDURE :: Get_matrix_from_bgcm_1l_2d
    PROCEDURE :: Get_matrix_from_bgcm_1l_3d
    GENERIC :: Get_matrix_from_bgcm_1l => Get_matrix_from_bgcm_1l_2d, Get_matrix_from_bgcm_1l_3d
    PROCEDURE :: Write_matrix_to_bgcm_2l_2d
    PROCEDURE :: Write_matrix_to_bgcm_2l_3d
    GENERIC :: Write_matrix_to_bgcm_2l => Write_matrix_to_bgcm_2l_2d, Write_matrix_to_bgcm_2l_3d
    PROCEDURE :: Write_matrix_to_bgcm_1l_2d
    PROCEDURE :: Write_matrix_to_bgcm_1l_3d
    GENERIC :: Write_matrix_to_bgcm_1l => Write_matrix_to_bgcm_1l_2d, Write_matrix_to_bgcm_1l_3d

    FINAL   :: Finalize_lnd_bgcm_store
  END TYPE

  CHARACTER(len=*), PARAMETER :: modname = 'mo_lnd_bgcm_store'

CONTAINS

  ! ======================================================================================================= !
  !>
  !> Init this instance of a bgc material store
  !>
  SUBROUTINE Init_bgcm_store(this)

    USE mo_jsb_parallel,         ONLY: Get_omp_no_of_threads
    ! -------------------------------------------------------------------------------------------------- !
    CLASS(t_lnd_bgcm_store), INTENT(INOUT) :: this !< to be initialised bgcm store
    ! -------------------------------------------------------------------------------------------------- !

    ! After counting the materials the store can be prepared to 
    ! 1. collect bookkeeping info about them
    ALLOCATE(this%bgc_material_name(this%nr_of_bgc_materials))
    ALLOCATE(this%bgc_material_ID(this%nr_of_bgc_materials))
    ALLOCATE(this%nr_of_levels(this%nr_of_bgc_materials))
    ALLOCATE(this%nr_of_dimensions(this%nr_of_bgc_materials))
    ALLOCATE(this%idx_in_store(this%nr_of_bgc_materials))

    ! 2. collect a flag used in assertions in the getter
    ALLOCATE(this%wrote_back_to_store(this%nr_of_bgc_materials,Get_omp_no_of_threads()))
    ! (initially no bgcms has been in use)
    this%wrote_back_to_store = .TRUE.

    ! 3. collect pointers to the bgc materials
    IF(this%nr_of_1l_2d_bgc_materials > 0) &
      & ALLOCATE(this%store_1l_2d_bgcms(this%nr_of_1l_2d_bgc_materials))
    IF(this%nr_of_1l_3d_bgc_materials > 0) &
      & ALLOCATE(this%store_1l_3d_bgcms(this%nr_of_1l_3d_bgc_materials))
    IF(this%nr_of_2l_2d_bgc_materials > 0) &
      & ALLOCATE(this%store_2l_2d_bgcms(this%nr_of_2l_2d_bgc_materials))
    IF(this%nr_of_2l_3d_bgc_materials > 0) &
      & ALLOCATE(this%store_2l_3d_bgcms(this%nr_of_2l_3d_bgc_materials))

  END SUBROUTINE Init_bgcm_store

  ! ======================================================================================================= !
  !>
  !> Append the bgcms of a process to this instance of a bgc material store
  !>
  SUBROUTINE Append_process_bgcms_to_store(this, main_bgcm, tile_name, idx_bgcm, idx_1l_2d_bgcm, &
      &                                    idx_2l_2d_bgcm, idx_1l_3d_bgcm, idx_2l_3d_bgcm)

    USE mo_jsb_pool_class,      ONLY: t_jsb_pool
    USE mo_lnd_bgcm_class,      ONLY: get_dimension_of_element_variables
    ! -------------------------------------------------------------------------------------------------- !
    CLASS(t_lnd_bgcm_store), INTENT(INOUT) :: this           !< this bgcm store instance
    CLASS(t_jsb_pool),       INTENT(IN)    :: main_bgcm      !< 'main' bgcm of a process
    CHARACTER(len=*),        INTENT(IN)    :: tile_name      !< name of the tile with this bgcm store
    INTEGER,                 INTENT(INOUT) :: idx_bgcm       !< idx of last added bgcm in the bgcm store
    INTEGER,                 INTENT(INOUT) :: idx_1l_2d_bgcm !< idx of last 1l 2d bgcm
    INTEGER,                 INTENT(INOUT) :: idx_2l_2d_bgcm !< idx of last 2l 2d bgcm
    INTEGER,                 INTENT(INOUT) :: idx_1l_3d_bgcm !< idx of last 1l 3d bgcm
    INTEGER,                 INTENT(INOUT) :: idx_2l_3d_bgcm !< idx of last 2l 3d bgcm
    ! -------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_pool),   POINTER  :: bgcm
    INTEGER                       :: i_bgcm, i_compartment, i_element
    INTEGER                       :: nr_of_bgcms_this_proc, nr_of_elements, nr_of_compartments

    CHARACTER(len=*), PARAMETER   :: routine = modname//':Append_process_bgcms_to_store'
    ! -------------------------------------------------------------------------------------------------- !

    IF (debug_on()) & 
      & CALL message(TRIM(routine), 'Append bgcms: '//TRIM(main_bgcm%name)//' - tile '//TRIM(tile_name))

    nr_of_bgcms_this_proc = SIZE(main_bgcm%pool_list)
    DO i_bgcm = 1, nr_of_bgcms_this_proc
      idx_bgcm = idx_bgcm + 1
      bgcm => main_bgcm%pool_list(i_bgcm)%p

      this%bgc_material_name(idx_bgcm) = bgcm%name
      IF(bgcm%id == -1) THEN
        CALL finish(TRIM(routine), 'Found bgcm with no specified id '//TRIM(bgcm%name)// &
          & ' in tile '//TRIM(tile_name)//', please check!')
      ENDIF
      this%bgc_material_ID(idx_bgcm) = bgcm%id

      ! Currently only two cases are expected
      ! - either directly elements and no compartments (1l)
      ! - or compartments with elements but no direct elements (2l)
      IF(bgcm%contains_elements) THEN
        ! i.e. directly elements and no compartments (asserted in the counting routine)
        this%nr_of_levels(idx_bgcm) = 1
        this%nr_of_dimensions(idx_bgcm) = get_dimension_of_element_variables(bgcm, tile_name, routine)
        nr_of_elements = SIZE(bgcm%element_list)

        ! Currently two cases are expected: 2D or 3D variables
        SELECT CASE(this%nr_of_dimensions(idx_bgcm))
        CASE(2)
          idx_1l_2d_bgcm = idx_1l_2d_bgcm + 1
          ALLOCATE(this%store_1l_2d_bgcms(idx_1l_2d_bgcm)%collection_1l_p(nr_of_elements))
          this%store_1l_2d_bgcms(idx_1l_2d_bgcm)%nr_of_elems = nr_of_elements
          this%idx_in_store(idx_bgcm) = idx_1l_2d_bgcm
          DO i_element = 1, nr_of_elements
            this%store_1l_2d_bgcms(idx_1l_2d_bgcm)%collection_1l_p(i_element)%p &
              & => bgcm%element_list(i_element)%p
          ENDDO
        CASE(3)
          idx_1l_3d_bgcm = idx_1l_3d_bgcm + 1
          ALLOCATE(this%store_1l_3d_bgcms(idx_1l_3d_bgcm)%collection_1l_p(nr_of_elements))
          this%store_1l_3d_bgcms(idx_1l_3d_bgcm)%nr_of_elems = nr_of_elements
          this%store_1l_3d_bgcms(idx_1l_3d_bgcm)%dim_len = bgcm%element_list(1)%p%len_dim3
          this%idx_in_store(idx_bgcm) = idx_1l_3d_bgcm
          DO i_element = 1, nr_of_elements
            this%store_1l_3d_bgcms(idx_1l_3d_bgcm)%collection_1l_p(i_element)%p &
              & => bgcm%element_list(i_element)%p
          ENDDO
        ENDSELECT

      ELSE
        ! i.e. no elements but one layer of compartments which are then expected to have elements
        ! but no further compartments (asserted in the counting routine) 
        this%nr_of_levels(idx_bgcm) = 2
        ! current assumption (asserted in the counting routine): all elements of all compartments have the same dimension
        this%nr_of_dimensions(idx_bgcm) &
          & = get_dimension_of_element_variables(bgcm%pool_list(1)%p, tile_name, routine)

        nr_of_compartments = SIZE(bgcm%pool_list)
        nr_of_elements = SIZE(bgcm%pool_list(1)%p%element_list)

        SELECT CASE(this%nr_of_dimensions(idx_bgcm))
        CASE(2)
          idx_2l_2d_bgcm = idx_2l_2d_bgcm + 1
          ALLOCATE(this%store_2l_2d_bgcms(idx_2l_2d_bgcm)%collection_2l_p(nr_of_compartments, nr_of_elements))
          this%store_2l_2d_bgcms(idx_2l_2d_bgcm)%nr_of_elems = nr_of_elements
          this%store_2l_2d_bgcms(idx_2l_2d_bgcm)%nr_of_parts = nr_of_compartments
          this%idx_in_store(idx_bgcm) = idx_2l_2d_bgcm
        CASE(3)
          idx_2l_3d_bgcm = idx_2l_3d_bgcm + 1
          ALLOCATE(this%store_2l_3d_bgcms(idx_2l_3d_bgcm)%collection_2l_p(nr_of_compartments, nr_of_elements))
          this%store_2l_3d_bgcms(idx_2l_3d_bgcm)%nr_of_elems = nr_of_elements
          this%store_2l_3d_bgcms(idx_2l_3d_bgcm)%nr_of_parts = nr_of_compartments
          this%idx_in_store(idx_bgcm) = idx_2l_3d_bgcm
        ENDSELECT

        DO i_compartment = 1, nr_of_compartments
          SELECT CASE(this%nr_of_dimensions(idx_bgcm))
          CASE(2)
            DO i_element = 1, nr_of_elements
              this%store_2l_2d_bgcms(idx_2l_2d_bgcm)%collection_2l_p(i_compartment, i_element)%p &
                & => bgcm%pool_list(i_compartment)%p%element_list(i_element)%p
            ENDDO
          CASE(3)
            this%store_2l_3d_bgcms(idx_2l_3d_bgcm)%dim_len &
              &  = bgcm%pool_list(1)%p%element_list(1)%p%len_dim3
            DO i_element = 1, nr_of_elements
              this%store_2l_3d_bgcms(idx_2l_3d_bgcm)%collection_2l_p(i_compartment, i_element)%p &
                & => bgcm%pool_list(i_compartment)%p%element_list(i_element)%p
            ENDDO
          ENDSELECT
        ENDDO ! compartments
      ENDIF ! bgcm type
    ENDDO ! bgcms

  END SUBROUTINE Append_process_bgcms_to_store


  ! ======================================================================================================= !
  !>
  !> Get a matrix for the chunk (ics:ice, iblk) of this 2l 2d bgcm
  !>
  SUBROUTINE Get_matrix_from_bgcm_2l_2d(this, bgcm_id, ics, ice, iblk, tile_name, bgcm_matrix)
    ! ----------------------------------------------------------------------------------------------------- !
    CLASS(t_lnd_bgcm_store),    INTENT(INOUT) :: this               !< this bgcm store instance
    INTEGER,                    INTENT(IN)    :: bgcm_id            !< id of bgcm queried from the store
    INTEGER,                    INTENT(IN)    :: ics                !< start index
    INTEGER,                    INTENT(IN)    :: ice                !< end index
    INTEGER,                    INTENT(IN)    :: iblk               !< current block
    CHARACTER (len=*),          INTENT(IN)    :: tile_name          !< name of the current tile
    REAL(wp),                   INTENT(INOUT) :: bgcm_matrix(:,:,:) !< dims: (part,elem,nc)
    ! ----------------------------------------------------------------------------------------------------- !
    INTEGER                                   :: bgcm_idx           !< index of bgcm queried from the store
    INTEGER                                   :: i_part, i_elem, nc

    CHARACTER(len=*), PARAMETER :: routine = modname//':Get_matrix_from_bgcm_2l_2d'
    ! ----------------------------------------------------------------------------------------------------- !
    bgcm_idx = get_bgcm_idx(this, bgcm_id, tile_name, routine)

    ! Check that this bgc material complies with the expectations
    nc = ice - ics + 1
    CALL assert_matrix_shape(this, bgcm_idx, nc, SHAPE(bgcm_matrix), tile_name, routine)
    CALL assert_wrote_back_to_store(this, bgcm_idx, Get_omp_thread(), tile_name, routine)

    ! 2levels 2D -> fill matrix (parts, elems, nc)
    DO i_part = 1,this%store_2l_2d_bgcms(this%idx_in_store(bgcm_idx))%nr_of_parts
      DO i_elem = 1,this%store_2l_2d_bgcms(this%idx_in_store(bgcm_idx))%nr_of_elems
        bgcm_matrix(i_part,i_elem,:) = &
          & this%store_2l_2d_bgcms(this%idx_in_store(bgcm_idx))%collection_2l_p(i_part, i_elem)%p%ptr2d(ics:ice,iblk)
      ENDDO
    ENDDO

    ! After getting the matrix the boolean is set to false to indicate that the bgcm is in use
    this%wrote_back_to_store(bgcm_idx, Get_omp_thread()) = .FALSE.

  END SUBROUTINE Get_matrix_from_bgcm_2l_2d


  ! ======================================================================================================= !
  !>
  !> Get a matrix for the chunk (ics:ice, iblk) of this 2l 3d bgcm
  !>
  SUBROUTINE Get_matrix_from_bgcm_2l_3d(this, bgcm_id, ics, ice, iblk, tile_name, bgcm_matrix)
    ! ----------------------------------------------------------------------------------------------------- !
    CLASS(t_lnd_bgcm_store),    INTENT(INOUT) :: this                 !< this bgcm store instance
    INTEGER,                    INTENT(IN)    :: bgcm_id              !< id of bgcm queried from the store
    INTEGER,                    INTENT(IN)    :: ics                  !< start index
    INTEGER,                    INTENT(IN)    :: ice                  !< end index
    INTEGER,                    INTENT(IN)    :: iblk                 !< current block
    CHARACTER (len=*),          INTENT(IN)    :: tile_name            !< name of the current tile
    REAL(wp),                   INTENT(INOUT) :: bgcm_matrix(:,:,:,:) !< dims: (part,elem,nc,layer)
    ! ----------------------------------------------------------------------------------------------------- !
    INTEGER                                   :: bgcm_idx           !< index of bgcm queried from the store
    INTEGER                                   :: i_part, i_elem, nc

    CHARACTER(len=*), PARAMETER :: routine = modname//':Get_matrix_from_bgcm_2l_2d'
    ! ----------------------------------------------------------------------------------------------------- !
    bgcm_idx = get_bgcm_idx(this, bgcm_id, tile_name, routine)

    ! Check that this bgc material complies with the expectations
    nc = ice - ics + 1
    CALL assert_matrix_shape(this, bgcm_idx, nc, SHAPE(bgcm_matrix), tile_name, routine)
    CALL assert_wrote_back_to_store(this, bgcm_idx, Get_omp_thread(), tile_name, routine)

    ! 2levels 3D -> fill matrix (parts, elems, nc, var_dim3)
    DO i_part = 1,this%store_2l_3d_bgcms(this%idx_in_store(bgcm_idx))%nr_of_parts
      DO i_elem = 1,this%store_2l_3d_bgcms(this%idx_in_store(bgcm_idx))%nr_of_elems
        bgcm_matrix(i_part,i_elem,:,:) &
          & = this%store_2l_3d_bgcms(this%idx_in_store(bgcm_idx))%collection_2l_p(i_part, i_elem)%p%ptr3d(ics:ice,:,iblk)
      ENDDO
    ENDDO

    ! After getting the matrix the boolean is set to false to indicate that the bgcm is in use
    this%wrote_back_to_store(bgcm_idx, Get_omp_thread()) = .FALSE.

  END SUBROUTINE Get_matrix_from_bgcm_2l_3d


  ! ======================================================================================================= !
  !>
  !> Get a matrix for the chunk (ics:ice, iblk) of this 1l 2d bgcm
  !>
  SUBROUTINE Get_matrix_from_bgcm_1l_2d(this, bgcm_id, ics, ice, iblk, tile_name, bgcm_matrix)
    ! ----------------------------------------------------------------------------------------------------- !
    CLASS(t_lnd_bgcm_store),    INTENT(INOUT) :: this             !< this bgcm store instance
    INTEGER,                    INTENT(IN)    :: bgcm_id          !< id of bgcm queried from the store
    INTEGER,                    INTENT(IN)    :: ics              !< start index
    INTEGER,                    INTENT(IN)    :: ice              !< end index
    INTEGER,                    INTENT(IN)    :: iblk             !< current block
    CHARACTER (len=*),          INTENT(IN)    :: tile_name        !< name of the current tile
    REAL(wp),                   INTENT(INOUT) :: bgcm_matrix(:,:) !< dims: (elem,nc)
    ! ----------------------------------------------------------------------------------------------------- !
    INTEGER                                   :: bgcm_idx           !< index of bgcm queried from the store
    INTEGER                                   :: i_elem, nc

    CHARACTER(len=*), PARAMETER :: routine = modname//':Get_matrix_from_bgcm_1l_2d'
    ! ----------------------------------------------------------------------------------------------------- !
    bgcm_idx = get_bgcm_idx(this, bgcm_id, tile_name, routine)

    ! Check that this bgc material complies with the expectations
    nc = ice - ics + 1
    CALL assert_matrix_shape(this, bgcm_idx, nc, SHAPE(bgcm_matrix), tile_name, routine)
    CALL assert_wrote_back_to_store(this, bgcm_idx, Get_omp_thread(), tile_name, routine)

    ! 1level 2D -> fill matrix (elems, nc, var_dim2)
    DO i_elem = 1,this%store_1l_2d_bgcms(this%idx_in_store(bgcm_idx))%nr_of_elems
      bgcm_matrix(i_elem,:) = &
        & this%store_1l_2d_bgcms(this%idx_in_store(bgcm_idx))%collection_1l_p(i_elem)%p%ptr2d(ics:ice,iblk)
    ENDDO

    ! After getting the matrix the boolean is set to false to indicate that the bgcm is in use
    this%wrote_back_to_store(bgcm_idx, Get_omp_thread()) = .FALSE.

  END SUBROUTINE Get_matrix_from_bgcm_1l_2d


  ! ======================================================================================================= !
  !>
  !> Get a matrix for the chunk (ics:ice, iblk) of this 1l 3d bgcm
  !>
  SUBROUTINE Get_matrix_from_bgcm_1l_3d(this, bgcm_id, ics, ice, iblk, tile_name, bgcm_matrix)
    ! ----------------------------------------------------------------------------------------------------- !
    CLASS(t_lnd_bgcm_store),    INTENT(INOUT) :: this               !< this bgcm store instance
    INTEGER,                    INTENT(IN)    :: bgcm_id            !< id of bgcm queried from the store
    INTEGER,                    INTENT(IN)    :: ics                !< start index
    INTEGER,                    INTENT(IN)    :: ice                !< end index
    INTEGER,                    INTENT(IN)    :: iblk               !< current block
    CHARACTER (len=*),          INTENT(IN)    :: tile_name          !< name of the current tile
    REAL(wp),                   INTENT(INOUT) :: bgcm_matrix(:,:,:) !< dims: (elem,nc,layer)
    ! ----------------------------------------------------------------------------------------------------- !
    INTEGER                                   :: bgcm_idx           !< index of bgcm queried from the store
    INTEGER                                   :: i_elem, nc

    CHARACTER(len=*), PARAMETER :: routine = modname//':Get_matrix_from_bgcm_1l_3d'
    ! ----------------------------------------------------------------------------------------------------- !
    bgcm_idx = get_bgcm_idx(this, bgcm_id, tile_name, routine)

    ! Check that this bgc material complies with the expectations
    nc = ice - ics + 1
    CALL assert_matrix_shape(this, bgcm_idx, nc, SHAPE(bgcm_matrix), tile_name, routine)
    CALL assert_wrote_back_to_store(this, bgcm_idx, Get_omp_thread(), tile_name, routine)

    ! 1level 3D -> fill matrix (elems, nc, var_dim3)
    DO i_elem = 1,this%store_1l_3d_bgcms(this%idx_in_store(bgcm_idx))%nr_of_elems
      bgcm_matrix(i_elem,:,:) = &
        & this%store_1l_3d_bgcms(this%idx_in_store(bgcm_idx))%collection_1l_p(i_elem)%p%ptr3d(ics:ice,:,iblk)
    ENDDO

    ! After getting the matrix the boolean is set to false to indicate that the bgcm is in use
    this%wrote_back_to_store(bgcm_idx, Get_omp_thread()) = .FALSE.

  END SUBROUTINE Get_matrix_from_bgcm_1l_3d


  ! ======================================================================================================= !
  !>
  !> Write matrix for the chunk (ics:ice, iblk) of this 2l 2d bgcm back
  !>
  SUBROUTINE Write_matrix_to_bgcm_2l_2d(this, bgcm_id, ics, ice, iblk, tile_name, bgcm_matrix)
    ! ----------------------------------------------------------------------------------------------------- !
    CLASS(t_lnd_bgcm_store),    INTENT(INOUT) :: this               !< this bgcm store instance
    INTEGER,                    INTENT(IN)    :: bgcm_id            !< id of bgcm queried from the store
    INTEGER,                    INTENT(IN)    :: ics                !< start index
    INTEGER,                    INTENT(IN)    :: ice                !< end index
    INTEGER,                    INTENT(IN)    :: iblk               !< current block
    CHARACTER (len=*),          INTENT(IN)    :: tile_name          !< name of the current tile
    REAL(wp),                   INTENT(IN)    :: bgcm_matrix(:,:,:) !< dims: (part,elem,nc)
    ! ----------------------------------------------------------------------------------------------------- !
    INTEGER                                   :: bgcm_idx           !< index of bgcm queried from the store
    INTEGER                                   :: i_part, i_elem, nc

    CHARACTER(len=*), PARAMETER :: routine = modname//':Write_matrix_to_bgcm_2l_2d'
    ! ----------------------------------------------------------------------------------------------------- !
    bgcm_idx = get_bgcm_idx(this, bgcm_id, tile_name, routine)

    ! Check that this bgc material complies with the expectations
    nc = ice - ics + 1
    CALL assert_matrix_shape(this, bgcm_idx, nc, SHAPE(bgcm_matrix), tile_name, routine)

    ! 2levels 2D -> write matrix back (parts, elems, nc)
    DO i_part = 1,this%store_2l_2d_bgcms(this%idx_in_store(bgcm_idx))%nr_of_parts
      DO i_elem = 1,this%store_2l_2d_bgcms(this%idx_in_store(bgcm_idx))%nr_of_elems
        this%store_2l_2d_bgcms(this%idx_in_store(bgcm_idx))%collection_2l_p(i_part, i_elem)%p%ptr2d(ics:ice,iblk) &
          & = bgcm_matrix(i_part,i_elem,:)
      ENDDO
    ENDDO

    ! After writing back the matrix the boolean is reset to true
    ! (This is no guarantee against an erroneous further use in the calling routine after writing it back!)
    this%wrote_back_to_store(bgcm_idx, Get_omp_thread()) = .TRUE.

  END SUBROUTINE Write_matrix_to_bgcm_2l_2d


  ! ======================================================================================================= !
  !>
  !> Write matrix for the chunk (ics:ice, iblk) of this 2l 3d bgcm back
  !>
  SUBROUTINE Write_matrix_to_bgcm_2l_3d(this, bgcm_id, ics, ice, iblk, tile_name, bgcm_matrix)
    ! ----------------------------------------------------------------------------------------------------- !
    CLASS(t_lnd_bgcm_store),    INTENT(INOUT) :: this                 !< this bgcm store instance
    INTEGER,                    INTENT(IN)    :: bgcm_id              !< id of bgcm queried from the store
    INTEGER,                    INTENT(IN)    :: ics                  !< start index
    INTEGER,                    INTENT(IN)    :: ice                  !< end index
    INTEGER,                    INTENT(IN)    :: iblk                 !< current block
    CHARACTER (len=*),          INTENT(IN)    :: tile_name            !< name of the current tile
    REAL(wp),                   INTENT(IN)    :: bgcm_matrix(:,:,:,:) !< dims: (part,elem,nc,layer)
    ! ----------------------------------------------------------------------------------------------------- !
    INTEGER                                   :: bgcm_idx             !< index of bgcm queried from the store
    INTEGER                                   :: i_part, i_elem, nc

    CHARACTER(len=*), PARAMETER :: routine = modname//':Write_matrix_to_bgcm_2l_2d'
    ! ----------------------------------------------------------------------------------------------------- !
    bgcm_idx = get_bgcm_idx(this, bgcm_id, tile_name, routine)

    ! Check that this bgc material complies with the expectations
    nc = ice - ics + 1
    CALL assert_matrix_shape(this, bgcm_idx, nc, SHAPE(bgcm_matrix), tile_name, routine)

    ! 2levels 3D -> write matrix back (parts, elems, nc, var_dim3)
    DO i_part = 1,this%store_2l_3d_bgcms(this%idx_in_store(bgcm_idx))%nr_of_parts
      DO i_elem = 1,this%store_2l_3d_bgcms(this%idx_in_store(bgcm_idx))%nr_of_elems
        this%store_2l_3d_bgcms(this%idx_in_store(bgcm_idx))%collection_2l_p(i_part, i_elem)%p%ptr3d(ics:ice,:,iblk) &
          & = bgcm_matrix(i_part,i_elem,:,:)
      ENDDO
    ENDDO

    ! After writing back the matrix the boolean is reset to true
    ! (This is no guarantee against an erroneous further use in the calling routine after writing it back!)
    this%wrote_back_to_store(bgcm_idx, Get_omp_thread()) = .TRUE.

  END SUBROUTINE Write_matrix_to_bgcm_2l_3d


  ! ======================================================================================================= !
  !>
  !> Write matrix for the chunk (ics:ice, iblk) of this 1l 2d bgcm back
  !>
  SUBROUTINE Write_matrix_to_bgcm_1l_2d(this, bgcm_id, ics, ice, iblk, tile_name, bgcm_matrix)
    ! ----------------------------------------------------------------------------------------------------- !
    CLASS(t_lnd_bgcm_store),    INTENT(INOUT) :: this             !< this bgcm store instance
    INTEGER,                    INTENT(IN)    :: bgcm_id          !< id of bgcm queried from the store
    INTEGER,                    INTENT(IN)    :: ics              !< start index
    INTEGER,                    INTENT(IN)    :: ice              !< end index
    INTEGER,                    INTENT(IN)    :: iblk             !< current block
    CHARACTER (len=*),          INTENT(IN)    :: tile_name        !< name of the current tile
    REAL(wp),                   INTENT(IN)    :: bgcm_matrix(:,:) !< dims: (elem,nc)
    ! ----------------------------------------------------------------------------------------------------- !
    INTEGER                                   :: bgcm_idx         !< index of bgcm queried from the store
    INTEGER                                   :: i_elem, nc

    CHARACTER(len=*), PARAMETER :: routine = modname//':Write_matrix_to_bgcm_1l_2d'
    ! ----------------------------------------------------------------------------------------------------- !
    bgcm_idx = get_bgcm_idx(this, bgcm_id, tile_name, routine)

    ! Check that this bgc material complies with the expectations
    nc = ice - ics + 1
    CALL assert_matrix_shape(this, bgcm_idx, nc, SHAPE(bgcm_matrix), tile_name, routine)

    ! 1level 2D -> write matrix back (elems, nc, var_dim2)
    DO i_elem = 1,this%store_1l_2d_bgcms(this%idx_in_store(bgcm_idx))%nr_of_elems
      this%store_1l_2d_bgcms(this%idx_in_store(bgcm_idx))%collection_1l_p(i_elem)%p%ptr2d(ics:ice,iblk) &
        & = bgcm_matrix(i_elem,:)
    ENDDO

    ! After writing back the matrix the boolean is reset to true
    ! (This is no guarantee against an erroneous further use in the calling routine after writing it back!)
    this%wrote_back_to_store(bgcm_idx, Get_omp_thread()) = .TRUE.

  END SUBROUTINE Write_matrix_to_bgcm_1l_2d


  ! ======================================================================================================= !
  !>
  !> Write matrix for the chunk (ics:ice, iblk) of this 1l 3d bgcm back
  !>
  SUBROUTINE Write_matrix_to_bgcm_1l_3d(this, bgcm_id, ics, ice, iblk, tile_name, bgcm_matrix)
    ! ----------------------------------------------------------------------------------------------------- !
    CLASS(t_lnd_bgcm_store),    INTENT(INOUT) :: this               !< this bgcm store instance
    INTEGER,                    INTENT(IN)    :: bgcm_id            !< id of bgcm queried from the store
    INTEGER,                    INTENT(IN)    :: ics                !< start index
    INTEGER,                    INTENT(IN)    :: ice                !< end index
    INTEGER,                    INTENT(IN)    :: iblk               !< current block
    CHARACTER (len=*),          INTENT(IN)    :: tile_name          !< name of the current tile
    REAL(wp),                   INTENT(IN)    :: bgcm_matrix(:,:,:) !< dims: (elem,nc,layer)
    ! ----------------------------------------------------------------------------------------------------- !
    INTEGER                                   :: bgcm_idx           !< index of bgcm queried from the store
    INTEGER                                   :: i_elem, nc

    CHARACTER(len=*), PARAMETER :: routine = modname//':Write_matrix_to_bgcm_1l_3d'
    ! ----------------------------------------------------------------------------------------------------- !
    bgcm_idx = get_bgcm_idx(this, bgcm_id, tile_name, routine)

    ! Check that this bgc material complies with the expectations
    nc = ice - ics + 1
    CALL assert_matrix_shape(this, bgcm_idx, nc, SHAPE(bgcm_matrix), tile_name, routine)

    ! 1level 3D -> write matrix back (elems, nc, var_dim3)
    DO i_elem = 1,this%store_1l_3d_bgcms(this%idx_in_store(bgcm_idx))%nr_of_elems
      this%store_1l_3d_bgcms(this%idx_in_store(bgcm_idx))%collection_1l_p(i_elem)%p%ptr3d(ics:ice,:,iblk) &
        & = bgcm_matrix(i_elem,:,:)
    ENDDO

    ! After writing back the matrix the boolean is reset to true
    ! (This is no guarantee against an erroneous further use in the calling routine after writing it back!)
    this%wrote_back_to_store(bgcm_idx, Get_omp_thread()) = .TRUE.

  END SUBROUTINE Write_matrix_to_bgcm_1l_3d


  ! ======================================================================================================= !
  !>
  !> Get the index of this bgcm (id) from the store
  !>
  FUNCTION get_bgcm_idx(this, bgcm_id, tile_name, caller) RESULT(bgcm_idx)
    ! ----------------------------------------------------------------------------------------------------- !
    CLASS(t_lnd_bgcm_store),    INTENT(IN)    :: this        !< this bgcm store instance
    INTEGER,                    INTENT(IN)    :: bgcm_id     !< id of bgcm queried from the store
    CHARACTER (len=*),          INTENT(IN)    :: tile_name   !< name of the current tile
    CHARACTER(len=*),           INTENT(in)    :: caller      !< calling routine
    ! ----------------------------------------------------------------------------------------------------- !
    INTEGER  :: bgcm_idx                         !< index of the queried bgcm in the store
    INTEGER  :: bgcm_min_idx(1), bgcm_max_idx(1) !< min and max index to check for uniqueness
    ! ----------------------------------------------------------------------------------------------------- !
    bgcm_min_idx = FINDLOC(this%bgc_material_ID, bgcm_id)
    bgcm_idx = bgcm_min_idx(1)

    bgcm_max_idx = FINDLOC(this%bgc_material_ID, bgcm_id, BACK=.TRUE.)
    IF (bgcm_idx == 0) THEN
      CALL finish(TRIM(caller), 'Did not find bgcm with given id '//int2string(bgcm_id)// &
        & ' in bgcm store of tile '//TRIM(tile_name)//', please check!')
    ELSEIF(bgcm_min_idx(1) /= bgcm_max_idx(1)) THEN
      CALL finish(TRIM(caller), 'Found not unique id '//int2string(bgcm_id)//                  &
        & ' used for '//TRIM(this%bgc_material_name(bgcm_idx))//' and '//                      &
        & TRIM(this%bgc_material_name(bgcm_max_idx(1)))//                                      &
        & ' in bgcm store of tile '//TRIM(tile_name)//', please check!')
    ENDIF

  END FUNCTION get_bgcm_idx


  ! ======================================================================================================= !
  !>
  !> Assert matching size and dimensions of the matrix to be filled from the store
  !>
  SUBROUTINE assert_matrix_shape(this, bgcm_idx, nc, dims, tile_name, caller)
    ! ----------------------------------------------------------------------------------------------------- !
    CLASS(t_lnd_bgcm_store),    INTENT(IN)    :: this       !< this bgcm store instance
    INTEGER,                    INTENT(IN)    :: bgcm_idx   !< index of the demanded bgcm in the store
    INTEGER,                    INTENT(IN)    :: nc         !< vector length
    INTEGER,                    INTENT(IN)    :: dims(:)    !< dimensions of the 
    CHARACTER (len=*),          INTENT(IN)    :: tile_name  !< name of the tile with this bgcm store
    CHARACTER (len=*),          INTENT(IN)    :: caller     !< calling routine
    ! ----------------------------------------------------------------------------------------------------- !
    INTEGER  :: i, nr_of_levels, nr_of_dims_bgcm, nr_dims_mt, nr_dims_expected
    INTEGER, ALLOCATABLE :: dims_expected(:)
    ! ----------------------------------------------------------------------------------------------------- !
    nr_of_levels = this%nr_of_levels(bgcm_idx)
    nr_of_dims_bgcm = this%nr_of_dimensions(bgcm_idx)
    nr_dims_mt = SIZE(dims)

    ! Currently there are 4 possible expected combinations
    IF(nr_of_levels == 2 .AND. nr_of_dims_bgcm == 3) THEN
      ! 1. 2levels 3D -> matrix with size 4
      nr_dims_expected = 4
      ALLOCATE(dims_expected(4))
      dims_expected(1) = this%store_2l_3d_bgcms(this%idx_in_store(bgcm_idx))%nr_of_parts
      dims_expected(2) = this%store_2l_3d_bgcms(this%idx_in_store(bgcm_idx))%nr_of_elems
      dims_expected(3) = nc
      dims_expected(4) = this%store_2l_3d_bgcms(this%idx_in_store(bgcm_idx))%dim_len
    ELSEIF(nr_of_levels == 2 .AND. nr_of_dims_bgcm == 2) THEN
      ! 2. 2levels 2D -> matrix with size 3
      nr_dims_expected = 3
      ALLOCATE(dims_expected(3))
      dims_expected(1) = this%store_2l_2d_bgcms(this%idx_in_store(bgcm_idx))%nr_of_parts
      dims_expected(2) = this%store_2l_2d_bgcms(this%idx_in_store(bgcm_idx))%nr_of_elems
      dims_expected(3) = nc
    ELSEIF(nr_of_levels == 1 .AND. nr_of_dims_bgcm == 3) THEN
      ! 3. 1level  3D -> matrix with size 3
      nr_dims_expected = 3
      ALLOCATE(dims_expected(3))
      dims_expected(1) = this%store_1l_3d_bgcms(this%idx_in_store(bgcm_idx))%nr_of_elems
      dims_expected(2) = nc
      dims_expected(3) = this%store_1l_3d_bgcms(this%idx_in_store(bgcm_idx))%dim_len
    ELSEIF(nr_of_levels == 1 .AND. nr_of_dims_bgcm == 2) THEN
      ! 4. 1level  2D -> matrix with size 2
      nr_dims_expected = 2
      ALLOCATE(dims_expected(2))
      dims_expected(1) = this%store_1l_2d_bgcms(this%idx_in_store(bgcm_idx))%nr_of_elems
      dims_expected(2) = nc
    ELSE
      CALL finish(TRIM(caller), 'Found unexpected combination of bgcm levels and dimensions ('// &
        & int2string(nr_of_levels)//' ,'//int2string(nr_of_dims_bgcm)//')'//                     &
        & ' for '//TRIM(this%bgc_material_name(bgcm_idx))//                                      &
        & ' in bgcm store of tile '//TRIM(tile_name)//', please check!')
    ENDIF

    IF(nr_dims_expected /= nr_dims_mt) THEN
      CALL finish(TRIM(caller), 'The number of dimensions of the passed bgcm matrix ('//  &
        & int2string(nr_dims_mt)//' does not comply with the expected dimension '//       &
        & int2string(nr_dims_expected)//' for '//TRIM(this%bgc_material_name(bgcm_idx))// &
        & ' in bgcm store of tile '//TRIM(tile_name)//', please check!')
    ENDIF

    DO i = 1, nr_dims_expected
      IF(dims_expected(i) /= dims(i)) THEN
        CALL finish(TRIM(caller), 'Mismatching bgcm matrix passed! For bgcm '//               &
          & TRIM(this%bgc_material_name(bgcm_idx))//                                          &
          & ' in bgcm store of tile '//TRIM(tile_name)// 'The size of dim '//                 &
          & int2string(i)//' does not comply with the expected dimension. Expected: '//       &
          & int2string(dims_expected(i))//', but found '//int2string(dims(i))//', please check!')
      ENDIF
    ENDDO

    DEALLOCATE(dims_expected)

  END SUBROUTINE assert_matrix_shape

  ! ======================================================================================================= !
  !>
  !> Assert wrote back (used to assert that a getter is used together with a write_back)
  !>
  SUBROUTINE assert_wrote_back_to_store(this, bgcm_idx, omp_thread, tile_name, caller)
    ! ----------------------------------------------------------------------------------------------------- !
    CLASS(t_lnd_bgcm_store),    INTENT(IN)    :: this       !< this bgcm store instance
    INTEGER,                    INTENT(IN)    :: bgcm_idx   !< index of the demanded bgcm in the store
    INTEGER,                    INTENT(IN)    :: omp_thread !< nr of this omp thread
    CHARACTER (len=*),          INTENT(IN)    :: tile_name  !< name of the tile with this bgcm
    CHARACTER (len=*),          INTENT(IN)    :: caller     !< calling routine
    ! ----------------------------------------------------------------------------------------------------- !
    IF(.NOT. this%wrote_back_to_store(bgcm_idx, omp_thread)) THEN
      CALL finish(TRIM(caller), 'Getter called for bgcm '//TRIM(this%bgc_material_name(bgcm_idx))// &
        & ' on thread '//int2string(omp_thread)//' on tile '//TRIM(tile_name)//                     &
        & ' with flag indicating missing wrote back from previous call, please check: '//           &
        & ' a getter should be followed by a write_back!')
    ENDIF

  END SUBROUTINE assert_wrote_back_to_store

  
  ! ======================================================================================================= !
  !>
  !> deallocate bgcm store components
  !>
  SUBROUTINE Finalize_lnd_bgcm_store(bgcm_store)
    ! ----------------------------------------------------------------------------------------------------- !
    TYPE(t_lnd_bgcm_store) :: bgcm_store
    ! ----------------------------------------------------------------------------------------------------- !
    INTEGER :: i
    ! ----------------------------------------------------------------------------------------------------- !

    IF(ALLOCATED(bgcm_store%bgc_material_ID)) DEALLOCATE(bgcm_store%bgc_material_ID)
    IF(ALLOCATED(bgcm_store%nr_of_levels)) DEALLOCATE(bgcm_store%nr_of_levels)
    IF(ALLOCATED(bgcm_store%nr_of_dimensions)) DEALLOCATE(bgcm_store%nr_of_dimensions)
    IF(ALLOCATED(bgcm_store%idx_in_store)) DEALLOCATE(bgcm_store%idx_in_store)
    IF(ALLOCATED(bgcm_store%wrote_back_to_store)) DEALLOCATE(bgcm_store%wrote_back_to_store)

    IF(ALLOCATED(bgcm_store%store_1l_2d_bgcms)) THEN
      DO i = 1,SIZE(bgcm_store%store_1l_2d_bgcms)
        IF(ALLOCATED(bgcm_store%store_1l_2d_bgcms(i)%collection_1l_p)) &
          & DEALLOCATE(bgcm_store%store_1l_2d_bgcms(i)%collection_1l_p)
      ENDDO
      DEALLOCATE(bgcm_store%store_1l_2d_bgcms)
    ENDIF

    IF(ALLOCATED(bgcm_store%store_1l_3d_bgcms)) THEN
      DO i = 1,SIZE(bgcm_store%store_1l_3d_bgcms)
        IF(ALLOCATED(bgcm_store%store_1l_3d_bgcms(i)%collection_1l_p)) &
          & DEALLOCATE(bgcm_store%store_1l_3d_bgcms(i)%collection_1l_p)
      ENDDO
      DEALLOCATE(bgcm_store%store_1l_3d_bgcms)
    ENDIF

    IF(ALLOCATED(bgcm_store%store_2l_2d_bgcms)) THEN
      DO i = 1,SIZE(bgcm_store%store_2l_2d_bgcms)
        IF(ALLOCATED(bgcm_store%store_2l_2d_bgcms(i)%collection_2l_p)) &
          & DEALLOCATE(bgcm_store%store_2l_2d_bgcms(i)%collection_2l_p)
      ENDDO
      DEALLOCATE(bgcm_store%store_2l_2d_bgcms)
    ENDIF

    IF(ALLOCATED(bgcm_store%store_2l_3d_bgcms)) THEN
      DO i = 1,SIZE(bgcm_store%store_2l_3d_bgcms)
        IF(ALLOCATED(bgcm_store%store_2l_3d_bgcms(i)%collection_2l_p)) &
          & DEALLOCATE(bgcm_store%store_2l_3d_bgcms(i)%collection_2l_p)
      ENDDO
      DEALLOCATE(bgcm_store%store_2l_3d_bgcms)
    ENDIF

  END SUBROUTINE Finalize_lnd_bgcm_store
  

#endif
END MODULE mo_lnd_bgcm_store
