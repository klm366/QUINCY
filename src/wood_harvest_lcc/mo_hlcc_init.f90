!>MR: work in progress: 
!> Status: 
!> 1st version: reads in annual harvesting area fractions (fraction of veg), for each age class separately -> "disturbance" on interface on AC level
!> this means, that needs to be read not just annually, but for every pft...should I make separate variables for each pfts? -> maybe easiest, and the name of pft could be read
!>  -issues: by default this module queued once a year, could theoretically probs have more dimensions for the input matrix, but current dsls for 3D I think
!>
!>
!>
!>To DO: 
!> check that number of age classes (nacs) match to children
!> IMPORTANT: need to somewhere check that the input tiles match the ac tiles - in practice, only number seems practical...distribution trickier - maybe through name? eg. IAS11/EAS11 vs FAGE?
!> read in "real area?" -> convert to fract of veg tile (or absolute grid)
!> read in biomass -> same, convert
!> make read in both AC level harv. and just pft harv. -> actually the latter could also work with FAGE on, so that there's just a rule for where to harvest first (this would happen in interface probs. maybe with process)
!> read in daily input
!>
!> 
!>
!>hlcc (wood harvesting land cover change) memory initialisation
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
!>#### Initialisation and reading of hlcc (wood harvesting) variables
!>
MODULE mo_hlcc_init
#ifndef __NO_JSBACH__

  USE mo_kind,              ONLY: wp
  USE mo_exception,         ONLY: message, finish

  USE mo_jsb_model_class,   ONLY: t_jsb_model
  USE mo_jsb_class,         ONLY: get_model
  USE mo_jsb_tile_class,    ONLY: t_jsb_tile_abstract
  USE mo_jsb_io_netcdf,     ONLY: t_input_file, jsb_netcdf_open_input
  USE mo_jsb_time,          ONLY: get_year
  USE mo_jsb_control,       ONLY: debug_on

  dsl4jsb_Use_processes HLCC_
  dsl4jsb_Use_memory(HLCC_)
  dsl4jsb_Use_config(HLCC_)
  dsl4jsb_Use_config(FAGE_)


  IMPLICIT NONE
  PRIVATE

  PUBLIC :: hlcc_init, read_harvest_data

  CHARACTER(len=*), PARAMETER :: modname = 'mo_hlcc_init'

CONTAINS

  ! ====================================================================================================== !
  !
  !> Initialize hlcc process (after memory has been set up)
  !
  SUBROUTINE hlcc_init(tile)

    ! -------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    ! -------------------------------------------------------------------------------------------------- !
    CHARACTER(len=*), PARAMETER :: routine = modname//':hlcc_init'

    dsl4jsb_Def_memory(HLCC_)

!    dsl4jsb_Real2D_onDomain ::
    ! -------------------------------------------------------------------------------------------------- !

    dsl4jsb_Get_memory(HLCC_)

!    dsl4jsb_Get_var2D_onDomain(HLCC_, )

    ! currently nothing to do here

  END SUBROUTINE hlcc_init


  ! ====================================================================================================== !
  !
  !> Reading of land use data
  !
  SUBROUTINE read_harvest_data(model_id, current_datetime)

    USE mo_jsb_grid_class,      ONLY: t_jsb_grid
    USE mo_jsb_grid,            ONLY: Get_grid
    USE mo_jsb_time_iface,      ONLY: t_datetime
    USE mo_jsb_varlist,         ONLY: VARNAME_LEN
    USE mo_io_units,            ONLY: filename_max
    USE mo_jsb_parallel,        ONLY: Get_omp_thread
!    USE mo_jsb_lcc_class,       ONLY: min_tolerated_cf_mismatch !-----MR-TODO: this is from alcc - should there be some similar test? It's removal (harvest) but I guess the absolute frac still def. can't be >1

    ! -------------------------------------------------------------------------------------------------- !
    INTEGER,                   INTENT(in) :: model_id
    TYPE(t_datetime), POINTER, INTENT(in) :: current_datetime
    ! -------------------------------------------------------------------------------------------------- !
    TYPE(t_jsb_model),           POINTER :: model
    TYPE(t_jsb_grid),            POINTER :: hgrid
    CLASS(t_jsb_tile_abstract),  POINTER :: tile, current_tile

    dsl4jsb_Def_config(HLCC_)
    dsl4jsb_Def_memory(HLCC_)
    dsl4jsb_Def_config(FAGE_)

    dsl4jsb_Real3D_onDomain :: harvest_fract_year
    REAL(wp),       POINTER :: ptr_2D(:,:)  ! tmp pointer
    REAL(wp), ALLOCATABLE   :: cf_correction(:,:) ! in case of relative fracts in harvest files: cf of the veg tile 

    INTEGER :: current_year, i_tile, npft, no_omp_thread, nproma, nblks !MR-TODO: check if all used
    CHARACTER(LEN=4) :: ac_name !MR
    TYPE(t_input_file)          :: input_file
    CHARACTER(len=filename_max) :: filename_harvesting_data
    CHARACTER(len=VARNAME_LEN)  :: varname_within_file
    CHARACTER(len=*), PARAMETER :: routine = modname//':read_harvest_data'
    ! -------------------------------------------------------------------------------------------------- !
    model => Get_model(model_id)
    no_omp_thread = Get_omp_thread()
    hgrid => Get_grid(model%grid_id)

    nproma = hgrid%nproma
    nblks  = hgrid%nblks

    IF (debug_on()) CALL message( TRIM(routine), 'Starting routine')
    current_year  = get_year(current_datetime)
 
    DO WHILE (ASSOCIATED(tile))  !MR_TODO: why is this necessary? what is the structure like, if in the init, need to check the tile is right? does it just shove all tiles through here?
      IF (tile%visited(no_omp_thread)) THEN  ! Have been here before
        CALL model%Goto_next_tile(tile)
        CYCLE
      END IF

      ! Check that the tile is a pft tile, and that it's a forest tile
      IF (INDEX(tile%name, 'pft') > 0 .AND. (dsl4jsb_Lctlib_param(woody_PFT)) .AND. tile%Has_children()) EXIT ! MR-TODO: one of these is probably redundant? (earlier: pft_tile%Is_process_calculated(FAGE_)). Does woody include shrubs?
 
      CALL model%Goto_next_tile(tile)
    END DO

    !>
    !> Assertion: hlcc is expected to run on forest pft tile 
    !>
    IF (.NOT. ASSOCIATED(tile)) THEN
      CALL finish(TRIM(routine), &
        & 'Violation of assertion: forest tile (with fage on*) required but not found in tile tree.') ! (MR-TODO: update printout later *)
    ENDIF

    ALLOCATE(cf_correction(nproma, nblks))

    dsl4jsb_Get_config(HLCC_)
    dsl4jsb_Get_config(FAGE_)

    
    ! Check that number of tile children is what's expected based on the number of ACs ! MR-TODO: need to check that input n:o of ACs matches fage config ACs (more important check - this prob in fage plus could be just crit)
    IF (.NOT. tile%Get_no_of_children() .EQ. dsl4jsb_Config(FAGE_)%nacs) THEN
      CALL finish(TRIM(routine), &
              & 'Violation of assertion: Expected nr of ACs for hlcc operations '&
              & //'did not match number of pft tile children, please check.')
    END IF

    ! check if model runs with relative or absolute fractions 
    IF(model%config%relative_fractions_in_file) THEN
      ! In this case read fractions have to be corrected with veg fract share
      CALL tile%Get_fraction(fract=cf_correction(:,:))
    ELSE
      ! If anyway already absolute fracts, cf_correction is not required
      cf_correction(:,:) = 1.0_wp
    ENDIF

    ! Read land use data
    dsl4jsb_Get_memory(HLCC_)
    dsl4jsb_Get_var3D_onDomain(HLCC_, harvest_fract_year) !-------------------------------------------MR-TODO: now vgrid should be AC's, but worth double checking later that it rly works
    !>
    !> Assertion: hlcc process currently expects filenames with 4 digits
    !>
    IF (( current_year > 9999) .OR. (current_year < 1000)) THEN 
      CALL finish(TRIM(routine), &
        & 'Violation of assertion: hlcc process currently expects filenames with 4 digits.')
    ENDIF
    WRITE(filename_harvesting_data,'(a,I4.4,a)') &
      & TRIM(dsl4jsb_Config(HLCC_)%hlcc_filename_prefix), current_year, ".nc" !---------------------MR-TODO: make sure file names fit this

    input_file = jsb_netcdf_open_input(TRIM(filename_harvesting_data), model%grid_id)

    current_tile => tile%Get_first_child_tile() !--------------------------------------------------MR-TODO: since tiles so far should've been pft, this should then go to AC level
    i_tile = 0 ! this should point to ACs
    
    DO WHILE (ASSOCIATED(current_tile))
      i_tile = i_tile + 1
      
      SELECT CASE (TRIM(dsl4jsb_Config(HLCC_)%scheme)) ! ------------------------------------------MR-TODO: needs to exist ("scheme") - not yet? does it matter here?? probs not, staying within pft tile
      CASE ('maps')
        !-- start: parts from Rainer Schnecks hlcc implementation
        ! Read this years pft-fraction from map file
        IF (debug_on()) CALL message(TRIM(routine), '...reading annual cover_fract maps for '//TRIM(current_tile%name))  !--MR-TODO: current_tile%name

        WRITE(ac_name, '(A2,I2.2)') '_ac', i_tile                          ! MR-TODO: added AC! - this is how the input should be formatted now
        WRITE(varname_within_file,'(a,a,a)') "harv_fract_", TRIM(current_tile%name), ac_name ! harv_fract_pftXX_acXX

        ptr_2D => input_file%Read_2d(variable_name=TRIM(varname_within_file), & !-----------------MR:TODO: fill array here!
          & fill_array = harvest_fract_year(:,i_tile,:) )   !MR-TODO: cf_current_year -> harvest_fract_year here
        !-- end: rewritten remainder from Rainer Schnecks hlcc implementation

        ! potentially correct read values (if relative fracts in map, else cf_correction has been set to 1.0_wp) !MR-TODO: confusingly phrased, which one's "correct"
          harvest_fract_year(:,i_tile,:) = harvest_fract_year(:,i_tile,:) * cf_correction !
      CASE DEFAULT
        CALL finish(TRIM(routine), 'Unimplemented hlcc scheme specified: '// TRIM(dsl4jsb_Config(HLCC_)%scheme) )
      END SELECT
      ! Assert: fractions >= 0 
      IF (ANY(harvest_fract_year(:,i_tile,:) < 0.0_wp)) THEN 
        CALL finish(TRIM(routine), &
        & 'Violation of assertion: hlcc land use data '//TRIM(filename_harvesting_data)//' attempts to harvest fractions <0 for '&
        & // TRIM(current_tile%name)//'. Please check.')
      ENDIF

      current_tile => current_tile%Get_next_sibling_tile()
    ENDDO

    CALL input_file%Close()

!    DEALLOCATE(cf_correction) !MR-TODO: did I comment this out? (probs, if no cf cor.))

!    ! Assert: sum of absolute fractions needs to be <=1
!    IF (ANY(SUM(cf_current_year(:,:,:),DIM=2) > 1.0_wp + min_tolerated_cf_mismatch)) THEN  !MR-TODO: -----edit this for some own test..? harvested frac is removal not the cf... but I guess if produces abs fractions>1, something def wrong. or harv>whole pft? (fixed/rule elsewhere?) !MR-TODO: cf_current_year -> harvest_fract_year here
!      CALL finish(TRIM(routine), &
!      & 'Violation of assertion: PFT fractions derived from '//TRIM(filename_harvesting_data) &
!      & //' sum to > 1.0_wp. Please check.')
!    ENDIF

    IF (debug_on()) CALL message( TRIM(routine), 'Finishing routine')

  END SUBROUTINE read_harvest_data

#endif
END MODULE mo_hlcc_init
