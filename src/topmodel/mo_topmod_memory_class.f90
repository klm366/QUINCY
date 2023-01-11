!>
!! @brief Contains the memory class for the soil and snow energy process.
!!
!! @author
!!  Stiig Wilkenskjeld, MPI-M Hamburg
!!
!! @par Revision History
!! Adopted from JSBACH3 code by Thomas Kleinen, 2022-12-06
!!
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/license .
!!

MODULE mo_top_memory_class
#ifndef __NO_JSBACH__

  USE mo_kind,                   ONLY: wp
  USE mo_exception,              ONLY: finish
  USE mo_util,                   ONLY: One_of
  USE mo_jsb_memory_class,       ONLY: t_jsb_memory
  USE mo_jsb_lct_class,          ONLY: VEG_TYPE, LAND_TYPE, BARE_TYPE
  USE mo_jsb_var_class,          ONLY: t_jsb_var_real2d, t_jsb_var_real3d

  ! Use of prcesses in this module
  dsl4jsb_Use_processes TOP_

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: t_top_memory, max_no_of_vars

  INTEGER, PARAMETER :: max_no_of_vars = 80

  ! ======================================================================================================= !
  !>Type definition for topmodel  memory
  !>
  TYPE, EXTENDS(t_jsb_memory) :: t_top_memory

    ! Parameters
    TYPE(t_jsb_var_real2d) :: &
      & frac_inun,            &                ! instant inundated fraction []
      & frac_inun_uncorr,     & ! StW: State?  ! instant inundated fraction, not corrected for liquid water fraction []
      & frac_inun_min_year,   &                ! minimum frac_inun in a year []
      & frac_inun_summ_year,  &                ! sum/avg of summer month mean of frac_inun []
      & frac_inun_min_lt,     &                ! multiyear avg of mimimum monthly frac_inun []
      & frac_inun_sum_month,  &                ! intramonthly sum of frac_inun []
      & frac_inun_mean_month, &                ! monthly mean of frac_inun []
      & frac_inun_summ_lt,    &                ! multiyear avg of mean monthly frac_inun []
      & water_table,          &                ! instant tile mean water table depth [m]
      & wt_inun,              &                ! instant water table depth [m] in inundated fraction [m]
      & wt_inun_min_year,     &                ! minimum monthly mean of wt_inun in current year [m]
      & wt_inun_sum_month,    &                ! intramonthly sum of wt_inun [m]
      & wt_inun_mean_month,   &                ! monthly mean of wt_inun [m]
      & wt_inun_min_lt,       &                ! multiyear avg of minimum monthly wt_inun [m] (StW: Long term mean)
      & aet,                  &                ! equilibrium evapotranspiration  (StW: Temporary storage)
!      & evaporation_pot,      &                ! potential evaporation (StW: avg_tiles output, diag)
      & poaet_lt,             &                ! (StW: Long term mean)
      & aet_sum_day,          &                ! equilibrium evapotranspiration sum over day (StW: To calc aet_mean_day)
      & aet_mean_day,         &                ! mean equilibrium evapotranspiration during day (StW: To calc poaet_mean_day)
      & poaet_mean_day,       &                ! (StW: to calc poaet_sum_year)
      & poaet_mean_year,      &                ! (StW: For long term means)
      & poaet_sum_year,       &                ! (StW: to calc poaet_mean_year)
      & precip_sum_day,       &                ! precipitation sum over day (StW: To calc precip_mean_day)
      & precip_mean_day,      &                ! mean precipitation during day (StW: To calc aet_mean_day)

      & frac_liquid_top,    &                  ! fraction of water in liquid form in uppermost soil layers (StW: avg_tiles output)
      & frac_liquid_top_uncorr,  &             ! fraction of water in liquid form in uppermost soil layers (uncorrected) (StW: diag)
      & peat_water_balance,   &                ! water balance of peatland tile (StW: to calc water_volume_peat)
      & peat_water_volume,    &                ! water volume in peatland tile (StW: state?)

      & peat_water_table,     &                ! peatland water table depth (StW: state?)
      & frac_permanent_wetland                 ! fraction of tile with permanent wetland (StW: state)

    TYPE(t_jsb_var_real3d) :: &
!      &  water_table_tile      &
      &  frac_liquid_sl

    REAL(wp), POINTER ::                &
      & cti_mean     (:,:)   => NULL(), &
      & cti_mu       (:,:)   => NULL(), &
      & cti_phi      (:,:)   => NULL(), &
      & cti_chi      (:,:)   => NULL(), &
      & cti_basin    (:,:)   => NULL(), &
      & cti_openwater(:,:)   => NULL(), &
      & inun_mon     (:,:,:) => NULL()

  CONTAINS
    PROCEDURE :: Init => Init_topmodel_memory
  END TYPE t_top_memory

  CHARACTER(len=*), PARAMETER :: modname = 'mo_top_memory_class'

CONTAINS

  ! ======================================================================================================= !
  !>initialize memory (variables) for the process: topmodel
  !>
  SUBROUTINE Init_topmodel_memory(mem, prefix, suffix, lct_ids, model_id)

    USE mo_jsb_varlist,       ONLY: BASIC, MEDIUM, FULL
    USE mo_jsb_io,            ONLY: grib_bits, t_cf, t_grib1, t_grib2, TSTEP_CONSTANT, TSTEP_INSTANT, tables
    USE mo_jsb_grid_class,    ONLY: t_jsb_grid, t_jsb_vgrid
    USE mo_jsb_grid,          ONLY: Get_grid, Get_vgrid
    USE mo_jsb_model_class,   ONLY: t_jsb_model
    USE mo_jsb_class,         ONLY: Get_model
    dsl4jsb_Use_config(TOP_)
    ! ----------------------------------------------------------------------------------------------------- !
    CLASS(t_top_memory), INTENT(inout), TARGET :: mem
    CHARACTER(len=*),    INTENT(in)            :: prefix
    CHARACTER(len=*),    INTENT(in)            :: suffix
    INTEGER,             INTENT(in)            :: lct_ids(:)
    INTEGER,             INTENT(in)            :: model_id
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Def_config(TOP_)
    TYPE(t_jsb_model), POINTER  :: model
    TYPE(t_jsb_grid),  POINTER  :: hgrid                        ! Horizontal grid
    TYPE(t_jsb_vgrid), POINTER  :: surface, soil_w              ! Vertical grids
    INTEGER                     :: table
    CHARACTER(len=*), PARAMETER :: routine = modname//':Init_topmodel_memory'
    INTEGER :: default_output_level = BASIC
    ! ----------------------------------------------------------------------------------------------------- !
    IF (SIZE(lct_ids) < 1) CALL finish(routine, 'lcd_ids not specified')
    ! ----------------------------------------------------------------------------------------------------- !
    model         => Get_model(model_id)
    dsl4jsb_Get_Config(TOP_)
    table         = tables(1)
    hgrid         => Get_grid(mem%grid_id)
    surface       => Get_vgrid('surface')
    soil_w        => Get_vgrid('soil_depth_water')

    !! !$acc enter data pcreate(mem)

    ! ----------------------------------------------------------------------------------------------------- !
    ! jsbach4 memory
    !
    IF (.NOT. dsl4jsb_Config(TOP_)%inun_prescribed) THEN
      ALLOCATE(mem%cti_mean     (hgrid%nproma, hgrid%nblks))
      ALLOCATE(mem%cti_basin    (hgrid%nproma, hgrid%nblks))
      ALLOCATE(mem%cti_openwater(hgrid%nproma, hgrid%nblks))
      ALLOCATE(mem%cti_mu       (hgrid%nproma, hgrid%nblks))
      ALLOCATE(mem%cti_phi      (hgrid%nproma, hgrid%nblks))
      ALLOCATE(mem%cti_chi      (hgrid%nproma, hgrid%nblks))
      !$acc enter data pcreate(mem%cti_mean,mem%cti_basin,mem%cti_openwater,mem%cti_mu,mem%cti_phi,mem%cti_chi)
      !$acc enter data attach(mem%cti_mean,mem%cti_basin,mem%cti_openwater,mem%cti_mu,mem%cti_phi,mem%cti_chi)
    ENDIF

    CALL mem%Add_var( 'water_table', mem%water_table,                        &
      & hgrid, surface,                                                      &
      & t_cf('water_table', 'm', ''),                                        & ! StW: To be refined
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
      & prefix, suffix,                                                      &
      & lrestart=.FALSE.,                                                    & ! StW: Right?
      & output_level = default_output_level,                                 &
      & initval_r=0.0_wp, isteptype=TSTEP_INSTANT )

    CALL mem%Add_var( 'frac_inun', mem%frac_inun,                            &
      & hgrid, surface,                                                      &
      & t_cf('frac_inun', '?', ''),                                          & ! StW: To be refined
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
      & prefix, suffix,                                                      &
      & lrestart=.TRUE.,                                                     & ! StW: Right?
      & output_level = default_output_level,                                 &
      & initval_r=0.0_wp, isteptype=TSTEP_INSTANT )

    CALL mem%Add_var( 'wt_inun', mem%wt_inun,                                &
      & hgrid, surface,                                                      &
      & t_cf('wt_inun', 'm', ''),                                            & ! StW: To be refined
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
      & prefix, suffix,                                                      &
      & lrestart=.FALSE.,                                                    & ! StW: Right?
      & output_level = default_output_level,                                 &
      & initval_r=0.0_wp, isteptype=TSTEP_INSTANT )

    CALL mem%Add_var( 'frac_inun_uncorr', mem%frac_inun_uncorr,              &
      & hgrid, surface,                                                      &
      & t_cf('frac_inun_uncorr', '?', ''),                                   & ! StW: To be refined
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
      & prefix, suffix,                                                      &
      & lrestart=.FALSE.,                                                    & ! StW: Right?
      & output_level = default_output_level,                                 &
      & initval_r=0.0_wp, isteptype=TSTEP_INSTANT )

    CALL mem%Add_var( 'frac_liquid_top', mem%frac_liquid_top,                &
      & hgrid, surface,                                                      &
      & t_cf('frac_liquid_top', '[]', ''),                                   & ! StW: To be refined
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
      & prefix, suffix,                                                      &
      & lrestart=.FALSE.,                                                    & ! StW: Right?
      & output_level = default_output_level,                                 &
      & initval_r=0.0_wp, isteptype=TSTEP_INSTANT )

    CALL mem%Add_var( 'frac_liquid_sl', mem%frac_liquid_sl,                  &
      & hgrid, soil_w,                                                       &
      & t_cf('frac_liquid_sl', '[]', ''),                                    & ! StW: To be refined
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
      & prefix, suffix,                                                      &
      & lrestart=.FALSE.,                                                    & ! StW: Right?
      & output_level = default_output_level,                                 &
      & initval_r=0.0_wp, isteptype=TSTEP_INSTANT )

    CALL mem%Add_var( 'frac_liquid_top_uncorr', mem%frac_liquid_top_uncorr,  &
      & hgrid, surface,                                                      &
      & t_cf('frac_liquid_top_uncorr', '[]', ''),                            & ! StW: To be refined
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
      & prefix, suffix,                                                      &
      & lrestart=.FALSE.,                                                    & ! StW: Right?
      & output_level = default_output_level,                                 &
      & initval_r=0.0_wp, isteptype=TSTEP_INSTANT )

    IF (dsl4jsb_Config(TOP_)%ltime_means) THEN
      CALL mem%Add_var( 'frac_inun_min_year', mem%frac_inun_min_year,          &
        & hgrid, surface,                                                      &
        & t_cf('frac_inun_min_year', '?', ''),                                 & ! StW: To be refined
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
        & prefix, suffix,                                                      &
        & lrestart=.FALSE.,                                                    & ! StW: Right?
        & output_level = default_output_level,                                 &
        & initval_r=1000.0_wp, isteptype=TSTEP_INSTANT )

      CALL mem%Add_var( 'frac_inun_summ_year', mem%frac_inun_summ_year,        &
        & hgrid, surface,                                                      &
        & t_cf('frac_inun_summ_year', '?', ''),                                & ! StW: To be refined
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
        & prefix, suffix,                                                      &
        & lrestart=.FALSE.,                                                    & ! StW: Right?
        & output_level = default_output_level,                                 &
        & initval_r=0.0_wp, isteptype=TSTEP_INSTANT )

      CALL mem%Add_var( 'frac_inun_sum_month', mem%frac_inun_sum_month,        &
        & hgrid, surface,                                                      &
        & t_cf('frac_inun_sum_month', '?', ''),                                & ! StW: To be refined
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
        & prefix, suffix,                                                      &
        & lrestart=.FALSE.,                                                    & ! StW: Right?
        & output_level = default_output_level,                                 &
        & initval_r=0.0_wp, isteptype=TSTEP_INSTANT )

      CALL mem%Add_var( 'frac_inun_mean_month', mem%frac_inun_mean_month,      &
        & hgrid, surface,                                                      &
        & t_cf('frac_inun_mean_month', '?', ''),                               & ! StW: To be refined
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
        & prefix, suffix,                                                      &
        & lrestart=.FALSE.,                                                    & ! StW: Right?
        & output_level = default_output_level,                                 &
        & initval_r=0.0_wp, isteptype=TSTEP_INSTANT )

      CALL mem%Add_var( 'frac_inun_min_lt', mem%frac_inun_min_lt,              &
        & hgrid, surface,                                                      &
        & t_cf('frac_inun_min_lt', '?', ''),                                   & ! StW: To be refined
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
        & prefix, suffix,                                                      &
        & lrestart=.FALSE.,                                                    & ! StW: Right?
        & output_level = default_output_level,                                 &
        & initval_r=0.0_wp, isteptype=TSTEP_INSTANT )

      CALL mem%Add_var( 'frac_inun_summ_lt', mem%frac_inun_summ_lt,            &
        & hgrid, surface,                                                      &
        & t_cf('frac_inun_summ_lt', '?', ''),                                  & ! StW: To be refined
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
        & prefix, suffix,                                                      &
        & lrestart=.FALSE.,                                                    & ! StW: Right?
        & output_level = default_output_level,                                 &
        & initval_r=0.0_wp, isteptype=TSTEP_INSTANT )

      CALL mem%Add_var( 'wt_inun_min_year', mem%wt_inun_min_year,              &
        & hgrid, surface,                                                      &
        & t_cf('wt_inun_min_year', 'm', ''),                                   & ! StW: To be refined
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
        & prefix, suffix,                                                      &
        & lrestart=.FALSE.,                                                    & ! StW: Right?
        & output_level = default_output_level,                                 &
        & initval_r=1000.0_wp, isteptype=TSTEP_INSTANT )

      CALL mem%Add_var( 'wt_inun_sum_month', mem%wt_inun_sum_month,            &
        & hgrid, surface,                                                      &
        & t_cf('wt_inun_sum_month', 'm', ''),                                  & ! StW: To be refined
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
        & prefix, suffix,                                                      &
        & lrestart=.FALSE.,                                                    & ! StW: Right?
        & output_level = default_output_level,                                 &
        & initval_r=0.0_wp, isteptype=TSTEP_INSTANT )

      CALL mem%Add_var( 'wt_inun_mean_month', mem%wt_inun_mean_month,          &
        & hgrid, surface,                                                      &
        & t_cf('wt_inun_mean_month', 'm', ''),                                 & ! StW: To be refined
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
        & prefix, suffix,                                                      &
        & lrestart=.FALSE.,                                                    & ! StW: Right?
        & output_level = default_output_level,                                 &
        & initval_r=0.0_wp, isteptype=TSTEP_INSTANT )

      CALL mem%Add_var( 'wt_inun_min_lt', mem%wt_inun_min_lt,                  &
        & hgrid, surface,                                                      &
        & t_cf('wt_inun_min_lt', 'm', ''),                                     & ! StW: To be refined
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
        & prefix, suffix,                                                      &
        & lrestart=.FALSE.,                                                    & ! StW: Right?
        & output_level = default_output_level,                                 &
        & initval_r=0.0_wp, isteptype=TSTEP_INSTANT )


      CALL mem%Add_var( 'frac_permanent_wetland', mem%frac_permanent_wetland,  &
        & hgrid, surface,                                                      &
        & t_cf('frac_permanent_wetland', '?', ''),                             & ! StW: To be refined
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
        & prefix, suffix,                                                      &
        & lrestart=.TRUE.,                                                     &
        & output_level = default_output_level,                                 &
        & initval_r=0.0_wp, isteptype=TSTEP_INSTANT )

    ENDIF
!      & evaporation_pot,      &                ! potential evaporation (StW: avg_tiles output, diag)

    IF (dsl4jsb_Config(TOP_)%lpeat) THEN

      CALL mem%Add_var( 'precip_sum_day', mem%precip_sum_day,                  &
        & hgrid, surface,                                                      &
        & t_cf('precip_sum_day', 'kg m-2 s-1', ''),                            &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
        & prefix, suffix,                                                      &
        & lrestart=.FALSE.,l_aggregate_all=.FALSE.,                            & ! Only works if restart is done at day start
        & output_level = default_output_level,                                 &
        & initval_r=0.0_wp, isteptype=TSTEP_INSTANT )

      CALL mem%Add_var( 'precip_mean_day', mem%precip_mean_day,                &
        & hgrid, surface,                                                      &
        & t_cf('precip_mean_day', 'kg m-2 s-1', ''),                           &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
        & prefix, suffix,                                                      &
        & lrestart=.TRUE.,l_aggregate_all=.TRUE.,                              &
        & output_level = default_output_level,                                 &
        & initval_r=0.0_wp, isteptype=TSTEP_INSTANT )

      CALL mem%Add_var( 'aet', mem%aet,                                        &
        & hgrid, surface,                                                      &
        & t_cf('aet', '?', ''),                                                & ! StW: To be refined
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
        & prefix, suffix,                                                      &
        & lrestart=.TRUE.,l_aggregate_all=.TRUE.,                              &
        & output_level = default_output_level,                                 &
        & initval_r=0.0_wp, isteptype=TSTEP_INSTANT )

      CALL mem%Add_var( 'peat_water_balance', mem%peat_water_balance,          &
        & hgrid, surface,                                                      &
        & t_cf('peat_water_balance', '?', ''),                                 & ! StW: To be refined
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
        & prefix, suffix,                                                      &
        & lrestart=.FALSE.,                                                    & ! StW: Right?
        & output_level = default_output_level,                                 &
        & initval_r=0.0_wp, isteptype=TSTEP_INSTANT )

      CALL mem%Add_var( 'peat_water_volume', mem%peat_water_volume,            &
        & hgrid, surface,                                                      &
        & t_cf('peat_water_volume', '?', ''),                                  & ! StW: To be refined
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
        & prefix, suffix,                                                      &
        & lrestart=.FALSE.,                                                    & ! StW: Right?
        & output_level = default_output_level,                                 &
        & initval_r=0.246_wp, isteptype=TSTEP_INSTANT )

      CALL mem%Add_var( 'peat_water_table', mem%peat_water_table,              &
        & hgrid, surface,                                                      &
        & t_cf('peat_water_table', '?', ''),                                   & ! StW: To be refined
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
        & prefix, suffix,                                                      &
        & lrestart=.FALSE.,                                                    & ! StW: Right?
        & output_level = default_output_level,                                 &
        & initval_r=0.0_wp, isteptype=TSTEP_INSTANT )

      IF (dsl4jsb_Config(TOP_)%ltime_means) THEN

        CALL mem%Add_var( 'aet_sum_day', mem%aet_sum_day,                        &
          & hgrid, surface,                                                      &
          & t_cf('aet_sum_day', '?', ''),                                        & ! StW: To be refined
          & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
          & prefix, suffix,                                                      &
          & lrestart=.TRUE.,l_aggregate_all=.TRUE.,                              &
          & output_level = default_output_level,                                 &
          & initval_r=0.0_wp, isteptype=TSTEP_INSTANT )

        CALL mem%Add_var( 'aet_mean_day', mem%aet_mean_day,                      &
          & hgrid, surface,                                                      &
          & t_cf('aet_mean_day', '?', ''),                                    & ! StW: To be refined
          & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
          & prefix, suffix,                                                      &
          & lrestart=.TRUE.,l_aggregate_all=.TRUE.,                              &
          & output_level = default_output_level,                                 &
          & initval_r=0.0_wp, isteptype=TSTEP_INSTANT )

        CALL mem%Add_var( 'poaet_mean_day', mem%poaet_mean_day,                  &
          & hgrid, surface,                                                      &
          & t_cf('poaet_mean_day', '?', ''),                                     & ! StW: To be refined
          & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
          & prefix, suffix,                                                      &
          & lrestart=.TRUE.,l_aggregate_all=.TRUE.,                              &
          & output_level = default_output_level,                                 &
          & initval_r=0.0_wp, isteptype=TSTEP_INSTANT )

        CALL mem%Add_var( 'poaet_sum_year', mem%poaet_sum_year,                  &
          & hgrid, surface,                                                      &
          & t_cf('poaet_sum_year', '?', ''),                                     & ! StW: To be refined
          & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
          & prefix, suffix,                                                      &
          & lrestart=.TRUE.,l_aggregate_all=.TRUE.,                              &
          & output_level = default_output_level,                                 &
          & initval_r=0.0_wp, isteptype=TSTEP_INSTANT )

        CALL mem%Add_var( 'poaet_mean_year', mem%poaet_mean_year,                &
          & hgrid, surface,                                                      &
          & t_cf('poaet_mean_year', '?', ''),                                    & ! StW: To be refined
          & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
          & prefix, suffix,                                                      &
          & lrestart=.TRUE.,l_aggregate_all=.TRUE.,                              &
          & output_level = default_output_level,                                 &
          & initval_r=0.0_wp, isteptype=TSTEP_INSTANT )

        CALL mem%Add_var( 'poaet_lt', mem%poaet_lt,                              &
          & hgrid, surface,                                                      &
          & t_cf('poaet_lt', '?', ''),                                           & ! StW: To be refined
          & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
          & prefix, suffix,                                                      &
          & lrestart=.TRUE.,l_aggregate_all=.TRUE.,                              &
          & output_level = default_output_level,                                 &
          & initval_r=0.0_wp, isteptype=TSTEP_INSTANT )

      ENDIF

    ENDIF

   END SUBROUTINE Init_topmodel_memory

#endif
END MODULE mo_top_memory_class
