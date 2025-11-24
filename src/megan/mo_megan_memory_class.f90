!>
!! @brief Contains the memory class for the topmodel wetland diagnostics
!!
!! @author
!!  Stiig Wilkenskjeld, MPI-M Hamburg
!!
!! @par Revision History
!! Adopted from JSBACH3 code by Thomas Kleinen, 2024-01-24
!!
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/license .
!!

MODULE mo_megan_memory_class
#ifndef __NO_JSBACH__

  USE mo_kind,                   ONLY: wp
  USE mo_exception,              ONLY: finish
  USE mo_util,                   ONLY: One_of
  USE mo_jsb_memory_class,       ONLY: t_jsb_memory
  USE mo_jsb_lct_class,          ONLY: VEG_TYPE, LAND_TYPE, BARE_TYPE
  USE mo_jsb_var_class,          ONLY: t_jsb_var_real2d
  USE mo_io_units,               ONLY: ll => filename_max


  ! Use of prcesses in this module
  dsl4jsb_Use_processes MEGAN_

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: t_megan_memory, max_no_of_vars
  PUBLIC :: t_megan_specie, t_megan_ptrs

  INTEGER, PARAMETER :: max_no_of_vars = 19

  ! ======================================================================================================= !
  !>Type definition for megan memory
  !>
  TYPE t_megan_specie
    ! Constants
    CHARACTER(len=ll)      :: longname        ! name of species
    CHARACTER(len=8)       :: shortname       ! chemical formula or other short name
    INTEGER                :: class           ! MEGAN compound class (1..20)
    INTEGER                :: id_tr           ! index of compound tracer in tracer list
    INTEGER                :: id_sp           ! index of compound species in species list
    REAL(wp)               :: mweight         ! molecular weight
    INTEGER                :: cnumber         ! carbon number
    REAL(wp)               :: beta            ! temperature adjustment factor (2012 paper, Table 4)
    REAL(wp)               :: ldf             ! light-dependent fraction (2012 paper, Table 4)
    INTEGER                :: cafac           ! relative emission rates for leaf age (2012 paper, Table 4)
    REAL(wp)               :: ef              ! PFT specific emission factor @ standard conditions
    REAL(wp)               :: scale           ! Emission factor scaling

    ! 2D Fields
    TYPE(t_jsb_var_real2d) ::              &
      & total_out,                         &  ! Emission [g m-2 s-1]
      & gamma_age,                         &
      & gamma_temp
  END TYPE t_megan_specie

  TYPE t_megan_ptrs
    REAL(wp), DIMENSION(:), POINTER ::     &
      & total_out,                         &  ! Emission [g m-2 s-1]
      & gamma_age,                         &
      & gamma_temp
  END TYPE t_megan_ptrs

  TYPE, EXTENDS(t_jsb_memory) :: t_megan_memory

    TYPE(t_megan_specie) ::                &
      & iso,                               &  ! Values and fields for isoprene
      & CO,                                &  ! Values and fields for CO
      & NO                                    ! Values and fields for NO

    TYPE(t_jsb_var_real2d) ::              &
      temp_relax_short,                    &  ! air temperature, 24h pseudo-mean [K]
      lai_relax_short,                     &  ! leaf area index, 24h pseudo-mean
      lai_relax_long,                      &  ! leaf area index, 30d pseudo-mean
      pfd_relax_short,                     &  ! photon flux density, 24h pseudo-mean
      pfd_relax_long,                      &  ! photon flux density, 30d pseudo-mean
      gamma_lai,                           &
      gamma_pfd,                           &
      gamma_CO2,                           &
      gamma_sm

  CONTAINS

    PROCEDURE :: Init => Init_megan_memory

  END TYPE t_megan_memory

  CHARACTER(len=*), PARAMETER :: modname = 'mo_megan_memory_class'

CONTAINS

  ! ======================================================================================================= !
  !>initialize memory (variables) for the process: megan
  !>
  SUBROUTINE Init_megan_memory(mem, prefix, suffix, lct_ids, lib_id, model_id)

    USE mo_jsb_varlist,       ONLY: BASIC, MEDIUM, FULL
    USE mo_jsb_io,            ONLY: grib_bits, t_cf, t_grib1, t_grib2, TSTEP_CONSTANT, TSTEP_INSTANT, tables
    USE mo_jsb_grid_class,    ONLY: t_jsb_grid, t_jsb_vgrid
    USE mo_jsb_grid,          ONLY: Get_grid, Get_vgrid
    USE mo_jsb_model_class,   ONLY: t_jsb_model
    USE mo_jsb_class,         ONLY: Get_model
    dsl4jsb_Use_config(MEGAN_)
    ! ----------------------------------------------------------------------------------------------------- !
    CLASS(t_megan_memory), INTENT(inout), TARGET :: mem
    CHARACTER(len=*),      INTENT(in)            :: prefix
    CHARACTER(len=*),      INTENT(in)            :: suffix
    INTEGER,               INTENT(in)            :: lct_ids(:)
    INTEGER,               INTENT(in)            :: lib_id        !< ID of the tile's primary land cover type (e.g. LAND_TYPE)
    INTEGER,               INTENT(in)            :: model_id
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Def_config(MEGAN_)
    TYPE(t_jsb_model), POINTER  :: model
    TYPE(t_jsb_grid),  POINTER  :: hgrid                        ! Horizontal grid
    TYPE(t_jsb_vgrid), POINTER  :: surface !, soil_w              ! Vertical grids
    INTEGER                     :: table
    CHARACTER(len=*), PARAMETER :: routine = modname//':Init_megang_memory'
    INTEGER :: default_output_level = BASIC
    ! ----------------------------------------------------------------------------------------------------- !
    IF (SIZE(lct_ids) < 1) CALL finish(routine, 'lcd_ids not specified')
    ! ----------------------------------------------------------------------------------------------------- !
    model         => Get_model(model_id)
    dsl4jsb_Get_Config(MEGAN_)
    table         = tables(1)
    hgrid         => Get_grid(mem%grid_id)
    surface       => Get_vgrid('surface')

    !! !$acc enter data pcreate(mem)

    ! ----------------------------------------------------------------------------------------------------- !
    ! jsbach4 memory
    !

    CALL Init_megan_specie('iso',mem%iso,.TRUE.)
    CALL Init_megan_specie('CO', mem%CO, .TRUE.)
    CALL Init_megan_specie('NO', mem%NO, .FALSE.)

    CALL mem%Add_var( 'gamma_lai', mem%gamma_lai,                            &
      & hgrid, surface,                                                      &
      & t_cf('gamma_lai','', '?'),                                           &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
      & prefix, suffix,                                                      &
      & lrestart=.TRUE.,                                                     &
      & output_level = FULL, l_aggregate_all = .FALSE.,                      &
      & initval_r=0.0_wp, isteptype=TSTEP_INSTANT )

    CALL mem%Add_var( 'gamma_pfd', mem%gamma_pfd,                            &
      & hgrid, surface,                                                      &
      & t_cf('gamma_pfd','', '?'),                                           &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
      & prefix, suffix,                                                      &
      & lrestart=.TRUE.,                                                     &
      & output_level = FULL, l_aggregate_all = .FALSE.,                      &
      & initval_r=0.0_wp, isteptype=TSTEP_INSTANT )

    CALL mem%Add_var( 'gamma_co2', mem%gamma_co2,                            &
      & hgrid, surface,                                                      &
      & t_cf('gamma_co2','', '?'),                                           &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
      & prefix, suffix,                                                      &
      & lrestart=.TRUE.,                                                     &
      & output_level = FULL, l_aggregate_all = .FALSE.,                      &
      & initval_r=0.0_wp, isteptype=TSTEP_INSTANT )

    CALL mem%Add_var( 'gamma_sm', mem%gamma_sm,                              &
      & hgrid, surface,                                                      &
      & t_cf('gamma_sm','', '?'),                                            &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
      & prefix, suffix,                                                      &
      & lrestart=.TRUE.,                                                     &
      & output_level = FULL, l_aggregate_all = .FALSE.,                      &
      & initval_r=0.0_wp, isteptype=TSTEP_INSTANT )

    CALL mem%Add_var( 'temp_relax_short', mem%temp_relax_short,              &
      & hgrid, surface,                                                      &
      & t_cf('temp_relax_short','K', 'Short term mean of surface temperature'), &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
      & prefix, suffix,                                                      &
      & lrestart=.TRUE.,                                                     &
      & output_level = FULL, l_aggregate_all = .FALSE.,                      &
      & initval_r=0.0_wp, isteptype=TSTEP_INSTANT )

    CALL mem%Add_var( 'lai_relax_short', mem%lai_relax_short,                &
      & hgrid, surface,                                                      &
      & t_cf('lai_relax_short','K', 'Short term mean of LAI'),               &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
      & prefix, suffix,                                                      &
      & lrestart=.TRUE.,                                                     &
      & output_level = FULL, l_aggregate_all = .FALSE.,                      &
      & initval_r=0.0_wp, isteptype=TSTEP_INSTANT )

    CALL mem%Add_var( 'lai_relax_long', mem%lai_relax_long,                  &
      & hgrid, surface,                                                      &
      & t_cf('lai_relax_long','K', 'Long term mean of LAI'),                 &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
      & prefix, suffix,                                                      &
      & lrestart=.TRUE.,                                                     &
      & output_level = FULL, l_aggregate_all = .FALSE.,                      &
      & initval_r=0.0_wp, isteptype=TSTEP_INSTANT )

    CALL mem%Add_var( 'pfd_relax_short', mem%pfd_relax_short,                &
      & hgrid, surface,                                                      &
      & t_cf('pfd_relax_short','K', 'Short term mean of pfd'),               &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
      & prefix, suffix,                                                      &
      & lrestart=.TRUE.,                                                     &
      & output_level = FULL, l_aggregate_all = .FALSE.,                      &
      & initval_r=0.0_wp, isteptype=TSTEP_INSTANT )

    CALL mem%Add_var( 'pfd_relax_long', mem%pfd_relax_long,                  &
      & hgrid, surface,                                                      &
      & t_cf('pfd_relax_long','K', 'Long term mean of pfd'),                 &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
      & prefix, suffix,                                                      &
      & lrestart=.TRUE.,                                                     &
      & output_level = FULL, l_aggregate_all = .FALSE.,                      &
      & initval_r=0.0_wp, isteptype=TSTEP_INSTANT )

!    IF (dsl4jsb_Config(TOP_)%lpeat) THEN

!    ENDIF
  CONTAINS

    SUBROUTINE Init_megan_specie(specie_name, specie, use_other, other_name)

      USE mo_jsb_grid_class,    ONLY: t_jsb_grid, t_jsb_vgrid
      USE mo_jsb_grid,          ONLY: Get_vgrid
      USE mo_exception,         ONLY: message_text, message
      USE mo_jsb_io,            ONLY: TSTEP_CONSTANT, TSTEP_INSTANT
      dsl4jsb_Use_config(MEGAN_)

      CHARACTER(len=*),            INTENT(in)      :: specie_name
      TYPE (t_megan_specie),       INTENT(inout)   :: specie
      LOGICAL,          OPTIONAL,  INTENT(in)      :: use_other
      CHARACTER(len=*), OPTIONAL,  INTENT(in)      :: other_name

      TYPE(t_jsb_vgrid), POINTER  :: surface                ! Vertical grids
      CHARACTER(len=256) :: str
      LOGICAL :: uother

      dsl4jsb_Def_config(MEGAN_)
      dsl4jsb_Get_config(MEGAN_)

      uother = .TRUE.; IF (PRESENT(use_other)) uother = use_other

      surface       => Get_vgrid('surface')

      specie%shortname = TRIM(specie_name)
      specie%scale = 0._wp
      specie%beta  = 0._wp

      CALL mem%Add_var( specie_name//'_total_out', specie%total_out,           &
        & hgrid, surface,                                                      &
        & t_cf(specie_name//'_total_out', 'kg m-2 s-1', ''),                   &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
        & prefix, suffix, lrestart = .FALSE.,                                  &
        & output_level = default_output_level, l_aggregate_all=.TRUE.,         &
        & initval_r=0.0_wp, isteptype=TSTEP_INSTANT )

      IF (uother) THEN
        CALL mem%Add_var( specie_name//'_gamma_age', specie%gamma_age,           &
          & hgrid, surface,                                                      &
          & t_cf(specie_name//'_gamma_age', '', ''),                             &
          & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
          & prefix, suffix, lrestart = .FALSE.,                                  &
          & output_level = FULL, l_aggregate_all=.FALSE.,                        &
          & initval_r=0.0_wp, isteptype=TSTEP_INSTANT )

        CALL mem%Add_var( specie_name//'_gamma_temp', specie%gamma_temp,         &
          & hgrid, surface,                                                      &
          & t_cf(specie_name//'_gamma_temp', '', ''),                            &
          & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
          & prefix, suffix, lrestart = .FALSE.,                                  &
          & output_level = FULL, l_aggregate_all=.FALSE.,                        &
          & initval_r=0.0_wp, isteptype=TSTEP_INSTANT )
      ENDIF

    END SUBROUTINE Init_megan_specie

  END SUBROUTINE Init_megan_memory

#endif
END MODULE mo_megan_memory_class
