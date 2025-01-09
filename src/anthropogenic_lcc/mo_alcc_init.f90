!> alcc (anthropogenic land cover change) memory initialisation
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
!>#### Initialisation and reading of alcc (anthropogenic land cover change) variables
!>
MODULE mo_alcc_init
#ifndef __NO_JSBACH__

  USE mo_kind,             ONLY: wp, dp
  USE mo_exception,         ONLY: message, finish

  USE mo_jsb_model_class,   ONLY: t_jsb_model
  USE mo_jsb_class,         ONLY: get_model
  USE mo_jsb_tile_class,    ONLY: t_jsb_tile_abstract
  USE mo_jsb_io_netcdf,     ONLY: t_input_file, jsb_netcdf_open_input
  USE mo_jsb_time,          ONLY: get_year
  USE mo_jsb_control,       ONLY: debug_on


  dsl4jsb_Use_processes ALCC_
  dsl4jsb_Use_memory(ALCC_)
  dsl4jsb_Use_config(ALCC_)

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: alcc_init, read_land_use_data

  CHARACTER(len=*), PARAMETER :: modname = 'mo_alcc_init'

CONTAINS

  ! ====================================================================================================== !
  !
  !> Initialize alcc process (after memory has been set up)
  !
  SUBROUTINE alcc_init(tile)

    ! -------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    ! -------------------------------------------------------------------------------------------------- !
    CHARACTER(len=*), PARAMETER :: routine = modname//':alcc_init'

    dsl4jsb_Def_memory(ALCC_)

!    dsl4jsb_Real2D_onDomain ::
    ! -------------------------------------------------------------------------------------------------- !

    dsl4jsb_Get_memory(ALCC_)

!    dsl4jsb_Get_var2D_onDomain(ALCC_, )

    ! currently nothing to do here

  END SUBROUTINE alcc_init


  ! ====================================================================================================== !
  !
  !> Reading of land use data
  !
  SUBROUTINE read_land_use_data(model_id, current_datetime, options)   ! Mingyue Zhang added "options" as input argument
  
    USE mo_kind,               ONLY: wp
    USE mo_jsb_grid_class,      ONLY: t_jsb_grid
    USE mo_jsb_grid,            ONLY: Get_grid
    USE mo_jsb_time_iface,      ONLY: t_datetime
    USE mo_jsb_varlist,         ONLY: VARNAME_LEN
    USE mo_io_units,            ONLY: filename_max
    USE mo_jsb_parallel,        ONLY: Get_omp_thread
    USE mo_jsb_lcc_class,       ONLY: min_tolerated_cf_mismatch
    USE mo_jsb_task_class,      ONLY: t_jsb_process_task, t_jsb_task_options

    ! -------------------------------------------------------------------------------------------------- !
    INTEGER,                   INTENT(in) :: model_id
    TYPE(t_datetime), POINTER, INTENT(in) :: current_datetime
    ! -------------------------------------------------------------------------------------------------- !
    TYPE(t_jsb_model),           POINTER       :: model
    TYPE(t_jsb_grid),            POINTER       :: hgrid
    CLASS(t_jsb_tile_abstract),  POINTER       :: tile, current_tile
    TYPE(t_jsb_task_options),    INTENT(in)    :: options  !< Additional run-time parameters
    ! -------------------------------------------------------------------------------------------------- !

    dsl4jsb_Def_config(ALCC_)
    dsl4jsb_Def_memory(ALCC_)
    dsl4jsb_Real3D_onDomain :: cf_current_year
    dsl4jsb_Real3D_onDomain :: cf_transition            ! Mingyue add dsl4jsb for transition vars
    dsl4jsb_Real3D_onDomain :: cf_transition_pfts       ! Mingyue add dsl4jsb for transition vars
    dsl4jsb_Real3D_onDomain :: cf_transition_gain       ! Mingyue add dsl4jsb for transition vars
    dsl4jsb_Real3D_onDomain :: cf_transition_loss       ! Mingyue add dsl4jsb for transition vars
    REAL(wp),                    POINTER :: ptr_2D(:,:)  ! tmp pointer
    REAL(wp), ALLOCATABLE   :: cf_correction(:,:) ! in case of relative fracts in map files: cf of the veg tile
    REAL(wp), ALLOCATABLE   :: transition_prob(:,:)  ! Mingyue add for transition probability to extend transition matrix

    INTEGER :: ic, nc, current_year, i_tile, npft, no_omp_thread, nproma, nblks, ics, ice
    INTEGER :: idx, i, j, row, col, trans_idx, pft_idx    ! Mingyue add for transition matrix index
    TYPE(t_input_file)          :: input_file
    TYPE(t_input_file)          :: input_file_trans                    ! Mingyue Zhang
    CHARACTER(len=filename_max) :: filename_land_use_data
    CHARACTER(len=filename_max) :: filename_land_transition_data       ! Mingyue Zhang
    CHARACTER(len=15), DIMENSION(9)  :: transitions                    ! Mingyue Zhang
    INTEGER,           DIMENSION(11) :: LC_type2LC_class_mapping       ! Mingyue Zhang
    CHARACTER(len=VARNAME_LEN)  :: varname_within_file
    CHARACTER(len=*), PARAMETER :: routine = modname//':read_land_use_data'
    ! -------------------------------------------------------------------------------------------------- !
    ! Get local variables from options argument
    
    model => Get_model(model_id)
    no_omp_thread = Get_omp_thread()
    hgrid => Get_grid(model%grid_id)

    nproma = hgrid%nproma
    nblks  = hgrid%nblks
    ics     = options%ics
    ice     = options%ice

    ! nc     = options%nc  nc = ice- ics + 1
    IF (debug_on()) CALL message( TRIM(routine), 'Starting routine')
    current_year  = get_year(current_datetime)

    ! alcc is expected to run on the veg tile
    ! -> search for the veg tile
    CALL model%Get_top_tile(tile)
    DO WHILE (ASSOCIATED(tile))
      IF (tile%visited(no_omp_thread)) THEN  ! Have been here before
        CALL model%Goto_next_tile(tile)
        CYCLE
      END IF

      ! JN-TODO: should such often required tile names be constants somewhere?
      IF (tile%name .EQ. 'veg') EXIT

      CALL model%Goto_next_tile(tile)
    ENDDO

    !>
    !> Assertion: alcc is expected to run on the veg tile
    !>
    IF (.NOT. ASSOCIATED(tile)) THEN
      CALL finish(TRIM(routine), &
        & 'Violation of assertion: veg tile required but not found in tile tree.')
    ENDIF

    ALLOCATE(cf_correction(nproma, nblks))

    dsl4jsb_Get_config(ALCC_)

    npft = tile%Get_no_of_children()
    IF (.NOT. npft .EQ. dsl4jsb_Config(ALCC_)%nr_of_pfts) THEN
      CALL finish(TRIM(routine), &
        & 'Violation of assertion: Expected nr of pfts for alcc operations '&
        & //'did not match number of veg tile children, please check.')
    ENDIF

    ! check if model runs with relative or absolute fractions
    IF(model%config%relative_fractions_in_file) THEN
      ! In this case read fractions have to be corrected with veg fract share
      CALL tile%Get_fraction(fract=cf_correction(:,:))
    ELSE
      ! If anyway already absolute fracts, cf_correction is not required
      cf_correction(:,:) = 1.0_wp
    ENDIF

    ! Read land use data
    dsl4jsb_Get_memory(ALCC_)
    dsl4jsb_Get_var3D_onDomain(ALCC_, cf_current_year)  
    dsl4jsb_Get_var3D_onDomain(ALCC_, cf_transition)        ! Mingyue add dsl4jsb for transition vars
    dsl4jsb_Get_var3D_onDomain(ALCC_, cf_transition_pfts)   ! Mingyue add dsl4jsb for transition vars
    dsl4jsb_Get_var3D_onDomain(ALCC_, cf_transition_gain)   ! Mingyue add dsl4jsb for transition vars
    dsl4jsb_Get_var3D_onDomain(ALCC_, cf_transition_loss)   ! Mingyue add dsl4jsb for transition vars
    !>
    !> Assertion: alcc process currently expects filenames with 4 digits
    !>
    IF (( current_year > 9999) .OR. (current_year < 1000)) THEN
      CALL finish(TRIM(routine), &
        & 'Violation of assertion: alcc process currently expects filenames with 4 digits.')
    ENDIF
    WRITE(filename_land_use_data,'(a,I4.4,a)') &
      & TRIM(dsl4jsb_Config(ALCC_)%alcc_filename_prefix), current_year, ".nc"

    input_file = jsb_netcdf_open_input(TRIM(filename_land_use_data), model%grid_id)
    
    !Mingyue Zhang define transition data name start
    WRITE(filename_land_transition_data,'(a,I4.4,a)') &
      & TRIM(dsl4jsb_Config(ALCC_)%alcc_transition_filename_prefix), current_year, ".nc"

    input_file_trans = jsb_netcdf_open_input(TRIM(filename_land_transition_data), model%grid_id)
    !Mingyue Zhang define transition data name end


    SELECT CASE (TRIM(dsl4jsb_Config(ALCC_)%scheme))
    CASE ('maps')
      !-- start: parts from Rainer Schnecks alcc implementation
        ! Read this years pft-fraction from map file
      IF (debug_on()) CALL message(TRIM(routine), '...reading annual cover_fract maps for '//TRIM(current_tile%name))
      WRITE(varname_within_file,'(a,a)') "fract_", TRIM(current_tile%name)

      current_tile => tile%Get_first_child_tile()
      i_tile = 0
      DO WHILE (ASSOCIATED(current_tile))
        i_tile = i_tile + 1
        
        ptr_2D => input_file%Read_2d(variable_name=TRIM(varname_within_file), &
        & fill_array = cf_current_year(:,i_tile,:) )
        !-- end: rewritten remainder from Rainer Schnecks alcc implementation

        ! potentially correct read values (if relative fracts in map, else cf_correction has been set to 1.0_wp)
        cf_current_year(:,i_tile,:) = cf_current_year(:,i_tile,:) * cf_correction
         
        ! Assert: fractions >= 0 and < 1
        IF (ANY(cf_current_year(:,i_tile,:) < 0.0_wp)) THEN
          CALL finish(TRIM(routine), &
          & 'Violation of assertion: Alcc land use data '//TRIM(filename_land_use_data)//' leads to fractions <0 for '&
          & // TRIM(current_tile%name)//'. Please check.') 
        ENDIF
        IF (ANY(cf_current_year(:,i_tile,:) > 1.0_wp + EPSILON(1._wp))) THEN
          CALL finish(TRIM(routine), &
          & 'Violation of assertion: Alcc land use data '//TRIM(filename_land_use_data)//' leads to fractions >1 for '&
          & // TRIM(current_tile%name)//'. Please check.') 
        ENDIF
    
        current_tile => current_tile%Get_next_sibling_tile()
      ENDDO
        
    CASE ('transitions')
      !Mingyue Zhang added for reading in transition data start 
      !here the loop is over each tile, meaning each pfts, but the matrix should be based on each grid box instead of tiles, 
      !so possibly here the CASE ('transitions') should not be inside the do while loop
      IF (debug_on()) CALL message(TRIM(routine), '...reading annual transitions for '//TRIM(current_tile%name))
      WRITE(varname_within_file,'(a,a)') "transition_", TRIM(current_tile%name)
      

      ! Loop over gridcell in the chunk to read in transition matrix data
      ! Define a 1D vector `TransMtrxVec` for the 3x3 matrix, dimensioned as 9 elements
      ! Read and set transition matrix
      DO idx = 1, 9 
        ! Calculate the original row (`i`) and column (`j`) positions from the 1D index `idx
        i = (idx - 1) / 3 + 1
        j = MOD(idx - 1, 3) + 1
      
        ! Set diagonal elements to 0.0; otherwise, read from input
        IF (i .EQ. j) THEN
          cf_transition(:, idx, :) = 0.0
        ELSE
          ptr_2D => input_file_trans%Read_2d(variable_name = TRIM(transitions(i*j)), &
          & fill_array = cf_transition(:, idx, :) )
        END IF
      END DO
      
      ! Loop over all pairs of PFTs and extend the transition matrix for PFTs
      DO idx = 1, npft * npft
        ! Calculate the row (`i`) and column (`j`) indices in the `npft x npft` matrix
        i = (idx - 1) / npft + 1
        j = MOD(idx - 1, npft) + 1
      
        ! Use the mappings to access the corresponding element in `TransMtrxVec`
        row = LC_type2LC_class_mapping(i)
        col = LC_type2LC_class_mapping(j)
        trans_idx = (row - 1) * 3 + col  ! Mapping to 3x3 transition matrix vector index
      
        transition_prob = cf_transition(:, trans_idx, :)
      
        ! Calculate transition probabilities for PFTs
        cf_transition_pfts(:, idx, :) = cf_current_year(:, i, :) * transition_prob * cf_current_year(:, j, :)
      END DO
      
      ! Mingyue Zhang added for reading in transition data end
      
      ! Initialize transition gain and loss
      cf_transition_gain = 0.0
      cf_transition_loss = 0.0
      
      ! Calculate transition gain and loss of each PFT in one loop
      DO pft_idx = 1, npft
        ! Calculate gain(pft_idx)
        cf_transition_gain(:, pft_idx, :) = 0.0
        DO j = 1, npft
          cf_transition_gain(:, pft_idx, :) = cf_transition_gain(:, pft_idx, :) + cf_current_year(:, j, :) * cf_transition_pfts(:, ((j - 1) * npft + pft_idx), :)
        END DO
      
        ! Calculate loss(pft_idx) within the inner loop
        cf_transition_loss(:, pft_idx, :) = 0.0
        DO j = 1, npft
          cf_transition_loss(:, pft_idx, :) = cf_transition_loss(:, pft_idx, :) + cf_current_year(:, pft_idx, :) * cf_transition_pfts(:, ((pft_idx - 1) * npft + j), :)
        END DO
      END DO

      ! Loop over the vegetation grid to update the current year's cover fractions
      ! by subtracting transition losses and adding transition gains for each tile.
      current_tile => tile%Get_first_child_tile()
      i_tile = 0
      
      DO WHILE (ASSOCIATED(current_tile))   
        i_tile = i_tile + 1
   
        cf_current_year(:, i_tile, :) = MAX(0.0_wp, cf_current_year(:, i_tile, :) + cf_transition_gain(:, i_tile, :) - cf_transition_loss(:, i_tile, :))
        !cf_current_year(:, i_tile, :) = cf_current_year(:, i_tile, :) + cf_transition_gain(:, i_tile, :) - cf_transition_loss(:, i_tile, :)
     
        ! Get the next sibling tile
        current_tile => current_tile%Get_next_sibling_tile()
      END DO
      
      
    CASE DEFAULT 
      CALL finish(TRIM(routine), 'Unimplemented alcc scheme specified: '// TRIM(dsl4jsb_Config(ALCC_)%scheme) )   
      
    
    END SELECT

    CALL input_file%Close()

    DEALLOCATE(cf_correction)

    ! Assert: sum of absolute fractions needs to be <=1
    IF (ANY(SUM(cf_current_year(:,:,:),DIM=2) > 1.0_wp + min_tolerated_cf_mismatch)) THEN
      CALL finish(TRIM(routine), &
      & 'Violation of assertion: PFT fractions derived from '//TRIM(filename_land_use_data) &
      & //' sum to > 1.0_wp. Please check.')
    ENDIF

    IF (debug_on()) CALL message( TRIM(routine), 'Finishing routine')

  
  END SUBROUTINE read_land_use_data
#endif
END MODULE mo_alcc_init
