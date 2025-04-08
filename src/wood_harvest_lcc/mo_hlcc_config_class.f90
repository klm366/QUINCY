! Meri's to-be-deleted notes:
! Done so far: 
! - copied the flcc process and swapped all flcc -> hlcc
!- updated some USE statements for HLCC
!
!> hlcc (lcc due to harvest) config
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
MODULE mo_hlcc_config_class
  #ifndef __NO_JSBACH__
  
    USE mo_exception,         ONLY: message, finish !MR-TODO: check if used (finish)
    USE mo_jsb_control,       ONLY: debug_on  !--check if used
    USE mo_io_units,          ONLY: filename_max
    USE mo_jsb_impl_constants, ONLY: SHORT_NAME_LEN
    !USE mo_kind,               ONLY: wp
    USE mo_jsb_config_class,  ONLY: t_jsb_config
    USE mo_jsb_model_class,    ONLY: t_jsb_model !MR: for reading the FAGE variable (->model) can be removed, if just having one fage 'scheme' as default (end checking elsewhere)!

    dsl4jsb_Use_processes FAGE_ !MR-TODO: not sure if this is good - "using other processes" in config
    dsl4jsb_Use_config(FAGE_)

  
    IMPLICIT NONE
    PRIVATE
    PUBLIC :: t_hlcc_config
  
    TYPE, EXTENDS(t_jsb_config) :: t_hlcc_config
        !LOGICAL :: separate  !MR-TODO: right now this is just kinda floating around, and should be sorted
        !INTEGER                       :: nacs ! MR-TODO: I don't think this is needed, cos number of ACs can be read from FAGE config (nacs) (I guess npft could be read elsewhere too...)
        CHARACTER(len=filename_max-7) :: hlcc_filename_prefix
          !< Prefix of hlcc (harvest input) filename prefix with len filename_max-7 (to add yyyy.nc)
        CHARACTER(len=SHORT_NAME_LEN) :: scheme
        LOGICAL                       :: l_daily_hlcc
          !< determines if land use change is to be applied on a daily (TRUE) or annual (FALSE) basis    
      CONTAINS
        PROCEDURE :: Init => Init_hlcc_config
    END type t_hlcc_config
  
    CHARACTER(len=*), PARAMETER :: modname = 'mo_hlcc_config_class'
  
  CONTAINS
  
    ! ====================================================================================================== !
    !
    !> Initialize hlcc process
    !
    ! -------------------------------------------------------------------------------------------------------
    SUBROUTINE Init_hlcc_config(config)
  
      USE mo_jsb_namelist_iface, ONLY: open_nml, POSITIONED, position_nml, close_nml
      USE mo_jsb_grid_class,     ONLY: t_jsb_vgrid!, new_vgrid !MR-TODO: some stuff not needed (commented) unless creating axis
      !USE mo_jsb_grid,           ONLY: Register_vgrid
      !USE mo_jsb_io,             ONLY: ZAXIS_GENERIC

      ! -------------------------------------------------------------------------------------------------- !
      CLASS(t_hlcc_config), INTENT(inout) :: config !<Configuration type of process (t_hlcc_config)
      ! -------------------------------------------------------------------------------------------------- !
      LOGICAL                       :: active, lrestart_cont
      CHARACTER(len=filename_max)   :: ic_filename, bc_filename
      !LOGICAL                      :: separate !< whether burned tiles should be treated separately or at once
      CHARACTER(len=filename_max-7) :: hlcc_filename_prefix
      CHARACTER(len=SHORT_NAME_LEN) :: scheme
      LOGICAL                       :: l_daily_hlcc
      !TYPE(t_jsb_vgrid), POINTER   :: ac_layer
      !INTEGER                      :: i, nacs - already dsl4jsb_Config(FAGE_)%nacs
      TYPE(t_jsb_model), POINTER    :: model ! MR: Get_config for FAGE_ uses the pointer "model" 


      dsl4jsb_Def_config(FAGE_)

  
      NAMELIST /jsb_hlcc_nml/      &
           active,                 &
           lrestart_cont,          & 
           ic_filename,            &
           bc_filename,            &
           !separate,               &
           hlcc_filename_prefix,   &
           scheme,                 &
           l_daily_hlcc           
           !nacs
  
      INTEGER :: nml_handler, nml_unit, istat
  
      CHARACTER(len=*), PARAMETER :: routine = modname//':Init_hlcc_config'
      ! -------------------------------------------------------------------------------------------------- !
  
      IF (debug_on()) CALL message(TRIM(routine), 'Starting hlcc configuration')
  
      ! Set defaults
      !separate            = .TRUE.
      active              = .FALSE.
      lrestart_cont        = .FALSE.       ! TRUE: Continue although HLCC variables are missing in restart file
      bc_filename         = 'bc_land_hlcc.nc'
      ic_filename         = 'ic_land_hlcc.nc'
      !nacs           = 11                     !MR-TODO: this is already set for fage, and for AC harvesting it should match
      scheme               = 'maps'  !MR-TODO: could this just be read from some other config, since I assume the scheme should be the same between processes..?
      l_daily_hlcc         = .FALSE. ! FALSE: annual instead as a default

      ! Age class spacing should match that of FAGE !MR: could also define in namelist though (wouldn't need FAGE_ here), and just have a check (in input) that they match

      dsl4jsb_Get_config(FAGE_)
      IF (dsl4jsb_Config(FAGE_)%age_class_scheme .EQ. "shorterYoung") THEN
        hlcc_filename_prefix = 'harv_frac_acs_IAS11_' ! Increasing age-spacing scheme (Nabel et al., 2020)
      ELSE IF (dsl4jsb_Config(FAGE_)%age_class_scheme .EQ. "equal") THEN
        hlcc_filename_prefix = 'harv_frac_acs_EAS11_' ! Filename for equal age spacing scheme
      END IF

      nml_handler = open_nml(TRIM(config%namelist_filename))
  
      nml_unit = position_nml('jsb_hlcc_nml', nml_handler, STATUS=istat)
      IF (istat == POSITIONED) READ(nml_unit, jsb_hlcc_nml)
  
      CALL close_nml(nml_handler)
  
      !config%separate             = separate
      config%active               = active
      config%ic_filename          = ic_filename
      config%bc_filename          = bc_filename
      config%lrestart_cont        = lrestart_cont
      config%hlcc_filename_prefix = hlcc_filename_prefix
      config%scheme               = scheme
      config%l_daily_hlcc         = l_daily_hlcc
      ! config%nacs            = nacs 
  
      SELECT CASE (TRIM(config%scheme))
      CASE ('maps')
        CALL message(TRIM(routine), 'Reading harvesting data from maps.')  !MR-TODO: is this even relevant for transitions??
      CASE DEFAULT
        CALL finish(TRIM(routine), 'Not implemented hlcc scheme specified: '// TRIM(config%scheme))
      END SELECT
  
      IF(l_daily_hlcc) THEN
        CALL message(TRIM(routine), '... harvesting will be applied on a daily basis.')
      ELSE
        CALL message(TRIM(routine), '... harvesting will be applied on an annual basis.')
      END IF
  
      ! MR-TODO: compare to alcc: here, vertical axis (vgrid) was created - 
      !...but, age_classes_grid is created by FAGE config, and as long as it's ran first, I think that's enough?
      !ac_layer  => new_vgrid('age_classes_grid', ZAXIS_GENERIC, nacs, &
      !& levels=(/ (REAL(i,kind=wp),i=1,nacs) /), units='')
      !CALL register_vgrid(ac_layer)

    END SUBROUTINE Init_hlcc_config
  
  #endif
  END MODULE mo_hlcc_config_class