!> QUINCY agriculture variables init
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
!> For more information on the QUINCY model see: <https://doi.org/10.17871/quincy-model-2019>
!>
!>#### initialization of agriculture memory variables using, e.g., ic & bc input files
!>
MODULE mo_q_agr_init
#ifndef __NO_JSBACH__

  USE mo_kind,                ONLY: lp
  USE mo_exception,           ONLY: message, finish
  USE mo_jsb_control,         ONLY: debug_on
  USE mo_jsb_class,           ONLY: Get_model
  USE mo_jsb_model_class,     ONLY: t_jsb_model
  USE mo_jsb_grid_class,      ONLY: t_jsb_grid
  USE mo_jsb_grid,            ONLY: Get_grid
  USE mo_jsb_tile_class,      ONLY: t_jsb_tile_abstract
  USE mo_jsb_process_class,   ONLY: Q_AGR_

  dsl4jsb_Use_config(Q_AGR_)

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: q_agr_init
#ifndef __QUINCY_STANDALONE__
  PUBLIC read_fertiliser_input
#endif

  ! -------------------------------------------------------------------------------------------------- !
  !
  !> Type to hold variables read from input files
  TYPE t_q_agr_init_vars
    REAL(lp), POINTER ::  &
      & c3_crop_type(:,:), &    !< major c3 crop
      & c4_crop_type(:,:)       !< major c4 crop
  END TYPE t_q_agr_init_vars

  TYPE(t_q_agr_init_vars) :: q_agr_init_vars

  CHARACTER(len=*), PARAMETER :: modname = 'mo_q_agr_init'

CONTAINS

  ! ======================================================================================================= !
  !> Run agriculture init
  !>
  SUBROUTINE q_agr_init(tile)
    USE mo_q_agr_constants,       ONLY: def_n_fertiliser_rate, def_n_manure_rate, ix_wheat, ix_tecorn, ix_trcorn, ix_trsoy, &
        &                               frac_fertil_nh4
    USE mo_q_assimi_constants,    ONLY: ic3phot
    dsl4jsb_Use_memory(Q_AGR_)
    ! ----------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    ! ----------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_tile_abstract), POINTER :: box_tile
    TYPE(t_jsb_model), POINTER          :: model
    TYPE(t_jsb_grid), POINTER           :: grid
    INTEGER                             :: iblk, ic, is         !< loop dimensions
    INTEGER                             :: nproma, nblks        !< dimensions
    INTEGER                             :: lctlib_ps_pathway
    CHARACTER(len=*), PARAMETER :: routine = modname//':q_agr_init'
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Def_memory_tile(Q_AGR_, box_tile)
    dsl4jsb_Real2D_onDomain :: nh4_fertiliser_c3crops
    dsl4jsb_Real2D_onDomain :: no3_fertiliser_c3crops
    dsl4jsb_Real2D_onDomain :: nh4_fertiliser_c4crops
    dsl4jsb_Real2D_onDomain :: no3_fertiliser_c4crops
    dsl4jsb_Real2D_onDomain :: n_manure_crop
    dsl4jsb_Def_memory(Q_AGR_)
    dsl4jsb_Real2D_onDomain :: crop_type_index
    dsl4jsb_Def_config(Q_AGR_)
    ! ----------------------------------------------------------------------------------------------------- !
    model  => Get_model(tile%owner_model_id)
    ! ----------------------------------------------------------------------------------------------------- !
    IF (debug_on()) CALL message(TRIM(routine), 'Setting initial conditions of agr memory (quincy) for tile '// &
      &                          TRIM(tile%name))
    ! ----------------------------------------------------------------------------------------------------- !
    grid   => Get_grid(model%grid_id)
    nblks  =  grid%nblks
    nproma =  grid%nproma
    CALL model%Get_top_tile(box_tile)
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Get_memory_tile(Q_AGR_, box_tile)
    dsl4jsb_Get_var2D_onDomain_tile(Q_AGR_, nh4_fertiliser_c3crops, box_tile) ! out
    dsl4jsb_Get_var2D_onDomain_tile(Q_AGR_, no3_fertiliser_c3crops, box_tile) ! out
    dsl4jsb_Get_var2D_onDomain_tile(Q_AGR_, nh4_fertiliser_c4crops, box_tile) ! out
    dsl4jsb_Get_var2D_onDomain_tile(Q_AGR_, no3_fertiliser_c4crops, box_tile) ! out
    dsl4jsb_Get_var2D_onDomain_tile(Q_AGR_, n_manure_crop, box_tile)          ! out
    dsl4jsb_Get_memory(Q_AGR_)
    dsl4jsb_Get_var2D_onDomain(Q_AGR_, crop_type_index)                ! out
    dsl4jsb_Get_config(Q_AGR_)
    ! ----------------------------------------------------------------------------------------------------- !

#ifdef __QUINCY_STANDALONE__
    lctlib_ps_pathway = dsl4jsb_Lctlib_param(ps_pathway)

    !---- For QS: Determine crop type index from latitudes
    ! set defaults based on photosynthesis type and latitude
    IF (lctlib_ps_pathway == ic3phot) THEN ! C3 Crops
      DO ic = 1,nproma
        DO iblk = 1,nblks
          IF (grid%lat(ic,iblk) < 30.0_lp .AND. grid%lat(ic,iblk) > -30.0_lp ) THEN
            crop_type_index(ic,iblk) = ix_trsoy ! tropical C3 crop = soybean
          ELSE
            crop_type_index(ic,iblk) = ix_wheat ! temperate C3 crop = wheat
          END IF
        END DO
      END DO
    ELSE ! C4 crops
      DO ic = 1,nproma
        DO iblk = 1,nblks
          IF (grid%lat(ic,iblk) < 30.0_lp .AND. grid%lat(ic,iblk) > -30.0_lp ) THEN
            crop_type_index(ic,iblk) = ix_trcorn ! tropical C4 crop = corn
          ELSE
            crop_type_index(ic,iblk) = ix_tecorn ! temperate C4 crop = corn
          END IF
        END DO
      END DO
    END IF

    ! default fertiliser used in QS (overwritten in IQ with constant or transient data in read_fertiliser_input)
    nh4_fertiliser_c3crops(:,:) = frac_fertil_nh4 * def_n_fertiliser_rate
    no3_fertiliser_c3crops(:,:) = (1._lp - frac_fertil_nh4) * def_n_fertiliser_rate
    nh4_fertiliser_c4crops(:,:) = frac_fertil_nh4 * def_n_fertiliser_rate
    no3_fertiliser_c4crops(:,:) = (1._lp - frac_fertil_nh4) * def_n_fertiliser_rate
    n_manure_crops(:,:)         = def_n_manure_rate
    ! ----------------------------------------------------------------------------------------------------- !
#else
    !---- For IQ: Read crop type from bc file (fertiliser is also read but in the time-loop not in the init)
    IF (.NOT. ASSOCIATED(tile%parent_tile)) THEN
      ! Read bc file on the box tile
      CALL q_agr_read_init_vars(tile)
    ELSE IF (tile%lcts(1)%lib_id /= 0) THEN
      IF (dsl4jsb_Lctlib_param(CropFlag)) THEN
        ! Determine crop_type_index for crops
        lctlib_ps_pathway = dsl4jsb_Lctlib_param(ps_pathway)

        IF (lctlib_ps_pathway == ic3phot) THEN
          ! C3 Crop
          crop_type_index(:,:) = q_agr_init_vars%c3_crop_type(:,:)
        ELSE
          ! C4 crop
          crop_type_index(:,:) = q_agr_init_vars%c4_crop_type(:,:)
        END IF
      END IF
    END IF

    IF (tile%Is_last_process_tile(Q_AGR_)) THEN
      CALL q_agr_finalise_init_vars()
    END IF

    !$ACC UPDATE DEVICE(crop_type_index) ASYNC(1)
    !$ACC WAIT(1)

#endif
  END SUBROUTINE q_agr_init

#ifndef __QUINCY_STANDALONE__
  ! ====================================================================================================== !
  !
  !> Reading of this years fertiliser input from input file
  !
  SUBROUTINE read_fertiliser_input(model_id, current_datetime)
    USE mo_io_units,                  ONLY: filename_max
    USE mo_jsb_time_iface,            ONLY: t_datetime
    USE mo_jsb_time,                  ONLY: is_time_experiment_start, get_year
    USE mo_jsb_io_netcdf,             ONLY: t_input_file, jsb_netcdf_open_input
    USE mo_q_agr_constants,           ONLY: background_n_fertiliser_rate, frac_fertil_nh4
    USE mo_jsb_physical_constants,    ONLY: molar_mass_N
    USE mo_jsb_math_constants,        ONLY: eps8

    dsl4jsb_Use_memory(Q_AGR_)
    ! -------------------------------------------------------------------------------------------------- !
    INTEGER,                   INTENT(in) :: model_id
    TYPE(t_datetime), POINTER, INTENT(in) :: current_datetime
    ! -------------------------------------------------------------------------------------------------- !
    TYPE(t_jsb_model),            POINTER :: model
    TYPE(t_jsb_grid),             POINTER :: hgrid
    CLASS(t_jsb_tile_abstract),   POINTER :: tile
    INTEGER                               :: nproma
    INTEGER                               :: nblks
    INTEGER                               :: iblk, ic         !< loop dimensions

    dsl4jsb_Def_config(Q_AGR_)
    dsl4jsb_Def_memory(Q_AGR_)
    dsl4jsb_Real2D_onDomain :: nh4_fertiliser_c3crops
    dsl4jsb_Real2D_onDomain :: no3_fertiliser_c3crops
    dsl4jsb_Real2D_onDomain :: nh4_fertiliser_c4crops
    dsl4jsb_Real2D_onDomain :: no3_fertiliser_c4crops
    dsl4jsb_Real2D_onDomain :: nh4_fertiliser_pasture
    dsl4jsb_Real2D_onDomain :: no3_fertiliser_pasture
    dsl4jsb_Real2D_onDomain :: n_manure_crop
    dsl4jsb_Real2D_onDomain :: n_manure_pasture

    ! -------------------------------------------------------------------------------------------------- !
    REAL(lp),           POINTER :: ptr_2D(:,:)  ! tmp pointer
    REAL(lp), ALLOCATABLE       :: hlp1(:,:)
    INTEGER                     :: current_year
    TYPE(t_input_file)          :: input_file
    CHARACTER(len=filename_max) :: filename_fertiliser_data
    CHARACTER(len=*), PARAMETER :: routine = modname//':read_fertiliser_input'
    ! -------------------------------------------------------------------------------------------------- !
    model => Get_model(model_id)
    hgrid  => Get_grid(model%grid_id)
    nproma = hgrid%Get_nproma()
    nblks  = hgrid%Get_nblks()
    CALL model%Get_top_tile(tile)

    IF (debug_on()) CALL message( TRIM(routine), 'Starting routine')

    dsl4jsb_Get_config(Q_AGR_)
    dsl4jsb_Get_memory(Q_AGR_)
    dsl4jsb_Get_var2D_onDomain(Q_AGR_, nh4_fertiliser_c3crops)
    dsl4jsb_Get_var2D_onDomain(Q_AGR_, no3_fertiliser_c3crops)
    dsl4jsb_Get_var2D_onDomain(Q_AGR_, nh4_fertiliser_c4crops)
    dsl4jsb_Get_var2D_onDomain(Q_AGR_, no3_fertiliser_c4crops)
    dsl4jsb_Get_var2D_onDomain(Q_AGR_, nh4_fertiliser_pasture)
    dsl4jsb_Get_var2D_onDomain(Q_AGR_, no3_fertiliser_pasture)
    dsl4jsb_Get_var2D_onDomain(Q_AGR_, n_manure_crop)
    dsl4jsb_Get_var2D_onDomain(Q_AGR_, n_manure_pasture)

    IF (.NOT. dsl4jsb_Config(Q_AGR_)%l_read_fertiliser) THEN
      ! ... use background rate
      nh4_fertiliser_c3crops(:,:) = frac_fertil_nh4 * background_n_fertiliser_rate
      no3_fertiliser_c3crops(:,:) = (1._lp - frac_fertil_nh4) * background_n_fertiliser_rate
      nh4_fertiliser_c4crops(:,:) = frac_fertil_nh4 * background_n_fertiliser_rate
      no3_fertiliser_c4crops(:,:) = (1._lp - frac_fertil_nh4) * background_n_fertiliser_rate
      ! no backround rate fertiliser on pasture
      nh4_fertiliser_pasture(:,:) = 0.0_lp
      no3_fertiliser_pasture(:,:) = 0.0_lp
    ELSE
      ! search a file ending on the current year
      current_year  = get_year(current_datetime)

      !>
      !> Assertion: routine currently expects filenames with 4 digits
      !>
      IF (( current_year > 9999) .OR. (current_year < 1000)) THEN
        CALL finish(TRIM(routine), 'Violation of assertion: this routine currently expects filenames with 4 digits.')
      END IF
      WRITE (filename_fertiliser_data,'(a,a,I4.4,a)') &
        & TRIM(dsl4jsb_Config(Q_AGR_)%fertiliser_filename_prefix), '_', current_year, ".nc"

      !> The 'LUH' fertiliser data is from LUH3, where the fertiliser is categorized as C3 and C4 fertiliser. No manure data in LUH.
      !> In constrast, the 'NMIP' data is from HaNi, where the fertiliser is categorized as NH4 and NO3, but no differences for C3 and C4 crops
      !> The 'NMIP' data also includes manure data, which further categorized as crop manure and pasture manure
      IF (TRIM(dsl4jsb_Config(Q_AGR_)%fertiliser_data_source) == 'LUH') THEN

        input_file = jsb_netcdf_open_input(TRIM(filename_fertiliser_data), model%grid_id)
        ptr_2D => input_file%Read_2d(variable_name='fertl_c3crops', fill_array = nh4_fertiliser_c3crops(:,:))
        ptr_2D => input_file%Read_2d(variable_name='fertl_c4crops', fill_array = nh4_fertiliser_c4crops(:,:))
        CALL input_file%Close()

        ! The actually applied fertiliser is the maximum of the default and the (converted) read fertiliser
        nh4_fertiliser_c3crops(:,:) = MAX(background_n_fertiliser_rate, &
          &                               nh4_fertiliser_c3crops(:,:) / 10._lp / molar_mass_N)
        nh4_fertiliser_c4crops(:,:) = MAX(background_n_fertiliser_rate, &
          &                               nh4_fertiliser_c4crops(:,:) / 10._lp / molar_mass_N)

        ! split LUH fertiliser into NH4 and NO3, since it's not known from the data source
        no3_fertiliser_c3crops(:,:) = (1._lp - frac_fertil_nh4) * nh4_fertiliser_c3crops(:,:)
        no3_fertiliser_c4crops(:,:) = (1._lp - frac_fertil_nh4) * nh4_fertiliser_c4crops(:,:)
        nh4_fertiliser_c3crops(:,:) = frac_fertil_nh4 * nh4_fertiliser_c3crops(:,:)
        no3_fertiliser_c4crops(:,:) = frac_fertil_nh4 * nh4_fertiliser_c4crops(:,:)

        ! no fertiliser on pasture in this data set
        nh4_fertiliser_pasture(:,:) = 0.0_lp
        no3_fertiliser_pasture(:,:) = 0.0_lp

      ELSE IF (TRIM(dsl4jsb_Config(Q_AGR_)%fertiliser_data_source) == 'NMIP') THEN

        ALLOCATE(hlp1(nproma, nblks))

        input_file = jsb_netcdf_open_input(TRIM(filename_fertiliser_data), model%grid_id)
        ptr_2D => input_file%Read_2d(variable_name='fertl_crops_nh4', fill_array = nh4_fertiliser_c3crops(:,:))
        ptr_2D => input_file%Read_2d(variable_name='fertl_crops_no3', fill_array = no3_fertiliser_c3crops(:,:))
        ptr_2D => input_file%Read_2d(variable_name='fertl_pasture_nh4', fill_array = nh4_fertiliser_pasture(:,:))
        ptr_2D => input_file%Read_2d(variable_name='fertl_pasture_no3', fill_array = no3_fertiliser_pasture(:,:))
        CALL input_file%Close()

        ! The actually applied fertiliser is the maximum of the default and the read fertiliser
        hlp1(:,:) = nh4_fertiliser_c3crops(:,:) + no3_fertiliser_c3crops(:,:)
        ! If fertiliser is zero in the input, applying the minimum fertiliser rate and default NH4/NO3 ratio
        DO ic = 1,nproma
          DO iblk = 1,nblks
            IF (hlp1(ic,iblk) <= eps8) THEN
              nh4_fertiliser_c3crops(ic,iblk) = frac_fertil_nh4 * background_n_fertiliser_rate
              no3_fertiliser_c3crops(ic,iblk) = (1._lp - frac_fertil_nh4) * background_n_fertiliser_rate
            ! IF fertiliser is non-zero but low, applying the minimum fertiliser rate but original NH4/NO3 ratio
            ELSE IF (hlp1(ic,iblk) > eps8 .AND. hlp1(ic,iblk) <= background_n_fertiliser_rate * 10._lp * molar_mass_N) THEN
              nh4_fertiliser_c3crops(ic,iblk) = nh4_fertiliser_c3crops(ic,iblk) / hlp1(ic,iblk) &
                &                           * background_n_fertiliser_rate
              no3_fertiliser_c3crops(ic,iblk) = no3_fertiliser_c3crops(ic,iblk) / hlp1(ic,iblk) &
                &                           * background_n_fertiliser_rate
            ELSE
          ! Convert unit from kg N ha-1 yr-1 to mol N m-2 yr-1
              nh4_fertiliser_c3crops(ic,iblk) = nh4_fertiliser_c3crops(ic,iblk) / 10._lp / molar_mass_N
              no3_fertiliser_c3crops(ic,iblk) = no3_fertiliser_c3crops(ic,iblk) / 10._lp / molar_mass_N
            END IF
          END DO
        END DO
        ! NMIP does not provide information for C3 vs. C4 crops, assert rates are the same
        ! Clean the temporary nh4_fertiliser_c4crops and use it to represent C4 crop nh4 fertiliser again
        nh4_fertiliser_c4crops(:,:) = nh4_fertiliser_c3crops(:,:)
        no3_fertiliser_c4crops(:,:) = no3_fertiliser_c3crops(:,:)

        ! The actually applied fertiliser is the (converted) read fertiliser, no background fertilisation assumed on pasture
        nh4_fertiliser_pasture(:,:) = nh4_fertiliser_pasture(:,:) / 10._lp / molar_mass_N
        no3_fertiliser_pasture(:,:) = no3_fertiliser_pasture(:,:) / 10._lp / molar_mass_N

        DEALLOCATE(hlp1)
      ELSE
        CALL finish(TRIM(routine), 'Invalid fertiliser data source. Please check your settings.')
      END IF

    END IF

    !$ACC UPDATE DEVICE(nh4_fertiliser_c3crops, nh4_fertiliser_c4crops, nh4_fertiliser_pasture) ASYNC(1)
    !$ACC UPDATE DEVICE(no3_fertiliser_c3crops, no3_fertiliser_c4crops, no3_fertiliser_pasture) ASYNC(1)
    !$ACC WAIT(1)

    IF (TRIM(dsl4jsb_Config(Q_AGR_)%fertiliser_data_source)=='NMIP' .AND. &
        &                             dsl4jsb_Config(Q_AGR_)%l_read_manure) THEN
      ! search a file ending on the current year
      current_year  = get_year(current_datetime)

      !
      ! Assertion: routine currently expects filenames with 4 digits, and that manure is provided in the same input
      !            file as N fertiliser
      !
      IF (( current_year > 9999) .OR. (current_year < 1000)) THEN
        CALL finish(TRIM(routine), 'Violation of assertion: this routine currently expects filenames with 4 digits.')
      END IF
      WRITE (filename_fertiliser_data,'(a,a,I4.4,a)') &
        & TRIM(dsl4jsb_Config(Q_AGR_)%fertiliser_filename_prefix), '_', current_year, ".nc"

      input_file = jsb_netcdf_open_input(TRIM(filename_fertiliser_data), model%grid_id)
      ptr_2D => input_file%Read_2d(variable_name='manure_crops', fill_array = n_manure_crop(:,:))
      ptr_2D => input_file%Read_2d(variable_name='manure_pasture', fill_array = n_manure_pasture(:,:))
      CALL input_file%Close()

      ! Convert unit of manure to model unit (kg / ha to mol/m2)
      n_manure_crop(:,:)    = n_manure_crop(:,:) / 10._lp / molar_mass_N
      n_manure_pasture(:,:) = n_manure_pasture(:,:) / 10._lp / molar_mass_N
    ELSE
      ! no manure is applied if not provided as forcing
      n_manure_crop(:,:)    = 0.0_lp
      n_manure_pasture(:,:) = 0.0_lp

    END IF

    !$ACC UPDATE DEVICE(n_manure_crop, n_manure_pasture) ASYNC(1)
    !$ACC WAIT(1)

    IF (debug_on()) CALL message(TRIM(routine), 'Finishing routine')

  END SUBROUTINE read_fertiliser_input

  ! ====================================================================================================== !
  !>
  !> Finalise initialisation variables
  !>
  SUBROUTINE q_agr_finalise_init_vars

    !> Deallocation of variables that had only been needed during the initialisation phase for reading from input
    DEALLOCATE( &
      & q_agr_init_vars%c3_crop_type, &
      & q_agr_init_vars%c4_crop_type  &
      & )

  END SUBROUTINE q_agr_finalise_init_vars


  ! ====================================================================================================== !
  !
  !> Read variables from bc file to module type variable used to initialise fertiliser input
  !>
  SUBROUTINE q_agr_read_init_vars(tile)
    ! ----------------------------------------------------------------------------------------------------- !
    USE mo_jsb_io_netcdf,       ONLY: t_input_file, jsb_netcdf_open_input
    USE mo_jsb_io,              ONLY: missval
    ! ----------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    ! ----------------------------------------------------------------------------------------------------- !
    TYPE(t_jsb_model), POINTER :: model
    TYPE(t_jsb_grid),  POINTER :: grid

    REAL(lp), POINTER  :: ptr_2D(:,:)
    TYPE(t_input_file) :: input_file
    INTEGER            :: nproma, nblks

    CHARACTER(len=*), PARAMETER :: routine = modname//':q_agr_read_init_vars'
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Def_config(Q_AGR_)
    ! ----------------------------------------------------------------------------------------------------- !
    model => get_model(tile%owner_model_id)
    grid   => Get_grid(model%grid_id)
    nproma = grid%Get_nproma()
    nblks  = grid%Get_nblks()

    ALLOCATE( &
      & q_agr_init_vars%c3_crop_type(nproma, nblks), &
      & q_agr_init_vars%c4_crop_type(nproma, nblks)  &
      & )

    q_agr_init_vars%c3_crop_type(:,:) = missval
    q_agr_init_vars%c4_crop_type(:,:) = missval

    dsl4jsb_Get_config(Q_AGR_)

    IF (debug_on()) CALL message(TRIM(routine), 'Reading agriculture init vars from ' &
      &                          //TRIM(dsl4jsb_Config(Q_AGR_)%bc_filename))

    input_file = jsb_netcdf_open_input(TRIM(dsl4jsb_Config(Q_AGR_)%bc_filename), model%grid_id)

    IF (debug_on()) CALL message(TRIM(routine), 'reading c3 crop type ...')
    ptr_2D => input_file%Read_2d(     &
      & variable_name='c3_crop_type', &
      & fill_array = q_agr_init_vars%c3_crop_type)

    IF (debug_on()) CALL message(TRIM(routine), 'reading c4 crop type ...')
    ptr_2D => input_file%Read_2d(     &
      & variable_name='c4_crop_type', &
      & fill_array = q_agr_init_vars%c4_crop_type)

  END SUBROUTINE q_agr_read_init_vars

#endif
#endif
END MODULE mo_q_agr_init
