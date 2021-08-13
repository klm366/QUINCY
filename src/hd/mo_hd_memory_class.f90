!>
!! @brief Contains the memory class for the surface.
!!
!! @author
!!  Veronika Gayler, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Veronika Gayler (2015-05-06)
!!
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/license .
!!

MODULE mo_hd_memory_class
!#ifndef __NO_JSBACH__
#if !defined(__NO_JSBACH__) && !defined(__NO_JSBACH_HD__)

  USE mo_kind, ONLY: wp

  !USE mo_jsb_model_class,  ONLY: t_jsb_model
  !USE mo_jsb_class,        ONLY: Get_model
  USE mo_jsb_memory_class, ONLY: t_jsb_memory
  USE mo_jsb_var_class,    ONLY: t_jsb_var_real2d, t_jsb_var_real3d
  USE mo_jsb_varlist,      ONLY: BASIC !, MEDIUM, FULL

  ! Use of prcesses in this module
  !dsl4jsb_Use_processes HD_

  ! Use process configurations
  !dsl4jsb_Use_config(HD_)

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: t_hd_memory, max_no_of_vars

  INTEGER, PARAMETER :: max_no_of_vars = 40

  TYPE, EXTENDS(t_jsb_memory) :: t_hd_memory

    TYPE(t_jsb_var_real2d) :: &
      hd_mask         , &  !< HD mask
      coast_ocean  ,    &  !< mask for coastal cells of the ocean
      ret_overlflow   , &  !< retention constant of overland flow
      ret_baseflow    , &  !< retention constant of baseflow
      ret_riverflow   , &  !< retention constant of riverflow
      nres_overlflow  , &  !< number of overland flow reservoirs
      nres_baseflow   , &  !< number of baseflow reservoirs
      nres_riverflow  , &  !< number of riverflow reservoirs
      outflow_runoff  , &  !< grid cell outflow due to runoff
      outflow_drainage, &  !< grid cell outflow due to drainage
      outflow_rivers  , &  !< grid cell outflow due to rivers
      volume_lake    , &  !< lake water volume storage of endorheic lakes [m^3]
      lake_outflow    , &  !< lake water outflow to neighbor grid cell [m^3/s]
      height_lake     , &  !< actual lake water height of endorheic lakes [m]
      area_lake       , &  !< actual lake water area of endorheic lakes [m^2]
      height_lake_max , &  !< maximum lake water height of endorheic lakes [m]
      area_lake_max   , &  !< maximum lake water area of endorheic lakes [m^2]
      volume_lake_max , &  !< maximum lake water volume of endorheic lakes [m^3]
      height_oro_min        , &  !< minimum elevation within ICON grid cell [m]
      height_dyn      , &  !< minimum dynamic height within ICON grid cell [m]
      alpha           , &  !< parameter to describe lake volume-height relation
      beta            , &  !< parameter to describe lake volume-area relation
      basin_mask      , &  !< lake basin grid points flagged as 1 else 0 
      lake_inflow     , &  !< lake inflow from runoff (precipitation minus evaporation) [m^3/s]
      lake_drain      , &  !< lake inflow from runoff and discharge into lake storage [m^3/s]
      water_budget    , &
      water_budget_old, &
      water_flux      , &
      water_error     , &
      outflow_resid        !< residual outflow that's not routed to discharge - only used for debugging (see mo_hd_interface.f90)

    TYPE(t_jsb_var_real3d) :: &
      nidx_upstream,          &  !< 1st dim. of local index of grid cell upstream (nproma)
      bidx_upstream,          &  !< 2nd dim. of local index of grid cell upstream (nblks)
      overlflow_res,          &  !< content of overland flow reservoir cascade
      baseflow_res,           &  !< content of baseflow reservoir cascade
      riverflow_res              !< content of groundflow reservoir cascade

    INTEGER ::               &
      nneigh,                &  !< maximum number of neighboring (inflow) cells
      nres_o_max,            &  !< maximum number of overland flow reservoirs
      nres_b_max,            &  !< maximum number of baseflow reservoirs
      nres_r_max                !< maximum number of riverflow reservoirs
  CONTAINS
    PROCEDURE :: Init => Init_hd_memory
  END TYPE t_hd_memory

  CHARACTER(len=*), PARAMETER :: modname = 'mo_hd_memory_class'

CONTAINS

  SUBROUTINE Init_hd_memory(mem, prefix, suffix, lct_ids, model_id)

    USE mo_jsb_io,            ONLY: grib_bits, t_cf, t_grib1, t_grib2, &
                                    TSTEP_CONSTANT, tables
    USE mo_jsb_grid_class,    ONLY: t_jsb_grid, t_jsb_vgrid
    USE mo_jsb_grid,          ONLY: Get_grid, Get_vgrid

    CLASS(t_hd_memory), INTENT(inout), TARGET :: mem
    CHARACTER(len=*),   INTENT(in)    :: prefix
    CHARACTER(len=*),   INTENT(in)    :: suffix
    INTEGER,            INTENT(in)    :: lct_ids(:)
    INTEGER,            INTENT(in)    :: model_id

    !dsl4jsb_Def_config(HD_)

    TYPE(t_jsb_grid),  POINTER :: hgrid                      ! Horizontal grid
    TYPE(t_jsb_vgrid), POINTER :: surface, hd_o, hd_b, hd_r  ! Vertical grids
    INTEGER :: table

    CHARACTER(len=*), PARAMETER :: routine = modname//':Init_hd_memory'

    IF (lct_ids(1) > 0)  CONTINUE ! avoid compiler warning about dummy argument not being used
    IF (model_id > 0)    CONTINUE ! avoid compiler warning about dummy argument not being used

    table = tables(1)

    hgrid   => Get_grid(mem%grid_id)
    surface => get_vgrid('surface')

    !vg   wie kann ich die level dynamisch definieren??
    hd_o => Get_vgrid('hd_nres_overlflow')
    hd_b => get_vgrid('hd_nres_baseflow')
    hd_r => get_vgrid('hd_nres_riverflow')

    !dsl4jsb_Get_config(HD_)

    CALL mem%Add_var( 'hd_mask', mem%hd_mask,     &
      & hgrid, surface,                                                                    &
      & t_cf('hd_mask', '', 'HD_ mask: -1 ocean, 0 ocean inflow, 1 land, 2 internal drainage'),                  &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                 &
      & prefix, suffix,                                                      &
      & lrestart=.FALSE.,                                                    &
      & initval_r=0.0_wp, isteptype=TSTEP_CONSTANT )

    CALL mem%Add_var( 'coast_ocean', mem%coast_ocean,   &
      & hgrid, surface,                                                               &
      & t_cf('coast_ocean', '', 'coastal grid cells of the ocean '),                  &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),            &
      & prefix, suffix,                                                               &
      & lrestart=.FALSE.,                                                    &
      & initval_r=0.0_wp, isteptype=TSTEP_CONSTANT )

    CALL mem%Add_var( prefix//'ret_overlflow'//suffix, mem%ret_overlflow,     &
      & hgrid, surface,                                                                    &
      & t_cf('ret_overlflow', 'day', 'overland flow retention constant'),                  &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                 &
      & prefix, suffix,                                                      &
      & lrestart=.FALSE.,                                                    &
      & initval_r=0.0_wp, isteptype=TSTEP_CONSTANT )

    CALL mem%Add_var( 'ret_baseflow', mem%ret_baseflow,       &
      & hgrid, surface,                                                                    &
      & t_cf('ret_baseflow', 'day', 'baseflow retention constant'),                        &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                 &
      & prefix, suffix,                                                      &
      & lrestart=.FALSE.,                                                    &
      & initval_r=0.0_wp, isteptype=TSTEP_CONSTANT )

    CALL mem%Add_var( 'ret_riverflow', mem%ret_riverflow,     &
      & hgrid, surface,                                                                    &
      & t_cf('ret_riverflow', 'day', 'riverflow retention constant'),                      &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                 &
      & prefix, suffix,                                                      &
      & lrestart=.FALSE.,                                                    &
      & initval_r=0.0_wp, isteptype=TSTEP_CONSTANT )

    CALL mem%Add_var( 'nres_overlflow', mem%nres_overlflow,   &
      & hgrid, surface,                                                                    &
      & t_cf('nres_overlflow', '', 'number of overland flow reservoirs'),                  &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                 &
      & prefix, suffix,                                                      &
      & lrestart=.FALSE.,                                                    &
      & initval_r=0.0_wp, isteptype=TSTEP_CONSTANT )

    CALL mem%Add_var( 'nres_baseflow', mem%nres_baseflow,     &
      & hgrid, surface,                                                                    &
      & t_cf('nres_baseflow', '', 'number of baseflow reservoirs'),                        &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                 &
      & prefix, suffix,                                                      &
      & lrestart=.FALSE.,                                                    &
      & initval_r=0.0_wp, isteptype=TSTEP_CONSTANT )

    CALL mem%Add_var( 'nres_riverflow', mem%nres_riverflow,           &
      & hgrid, surface,                                                                    &
      & t_cf('nres_riverflow', '', 'number of riverflow reservoirs'),                      &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                 &
      & prefix, suffix,                                                      &
      & lrestart=.FALSE.,                                                    &
      & initval_r=0.0_wp, isteptype=TSTEP_CONSTANT )

    CALL mem%Add_var( 'outflow_runoff', mem%outflow_runoff,           &
      & hgrid, surface,                                                                    &
      & t_cf('outflow_runoff', 'm3/s', 'grid cell outflow due to runoff'),                 &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                 &
      & prefix, suffix,                                                      &
      & initval_r=0.0_wp )

    CALL mem%Add_var( 'outflow_drainage', mem%outflow_drainage,       &
      & hgrid, surface,                                                                    &
      & t_cf('outflow_drainage', 'm3/s', 'grid cell outflow due to drainage'),             &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                 &
      & prefix, suffix,                                                      &
      & initval_r=0.0_wp )

    CALL mem%Add_var( 'outflow_rivers', mem%outflow_rivers,              &
      & hgrid, surface,                                                                       &
      & t_cf('outflow_rivers', 'm3/s', 'grid cell outflow due to rivers'),                    &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                    &
      & prefix, suffix,                                                      &
      & initval_r=0.0_wp )

     CALL mem%Add_var( 'volume_lake', mem%volume_lake,                     &
      & hgrid, surface,                                                      &
      & t_cf('volume_lake', 'm3', 'water content of the endorheic lake'),   &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
      & prefix, suffix,                                                      &
      & output_level=BASIC,                                                  &
      & initval_r=0.0_wp )

     CALL mem%Add_var( 'lake_outflow', mem%lake_outflow,                     &
      & hgrid, surface,                                                      &
      & t_cf('lake_outflow', 'm3/s', 'lake outflow to neighbor grid cell'),  &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
      & prefix, suffix,                                                      &
      & output_level=BASIC,                                                  &
      & initval_r=0.0_wp )

     CALL mem%Add_var( 'height_lake', mem%height_lake,                       &
      & hgrid, surface,                                                      &
      & t_cf('height_lake', 'm', 'grid cell lake height'),                   &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
      & prefix, suffix,                                                      & 
      & output_level=BASIC,                                                  &
      & initval_r=0.0_wp )

     CALL mem%Add_var( 'area_lake', mem%area_lake,                           &
      & hgrid, surface,                                                      &
      & t_cf('area_lake', 'm2', 'grid cell lake area'),                      &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
      & prefix, suffix,                                                      &
      & output_level=BASIC,                                                  &
      & initval_r=0.0_wp )

    CALL mem%Add_var( 'height_lake_max', mem%height_lake_max,                &
      & hgrid, surface,                                                      &
      & t_cf('height_lake_max', 'm', 'grid cell maxium lake height'),        &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
      & prefix, suffix,                                                      &
      & lrestart=.FALSE.,                                                    &
      & initval_r=0.0_wp, isteptype=TSTEP_CONSTANT )

    CALL mem%Add_var( 'area_lake_max', mem%area_lake_max,                    &
      & hgrid, surface,                                                      &
      & t_cf('area_lake_max', 'm2', 'grid cell maxium lake area'),           & 
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
      & prefix, suffix,                                                      &
      & lrestart=.FALSE.,                                                    &
      & initval_r=0.0_wp, isteptype=TSTEP_CONSTANT )

    CALL mem%Add_var( 'volume_lake_max', mem%volume_lake_max,                &
      & hgrid, surface,                                                      &
      & t_cf('volume_lake_max', 'm3', 'grid cell maxium lake volume'),       &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
      & prefix, suffix,                                                      &
      & lrestart=.FALSE.,                                                    &
      & initval_r=0.0_wp, isteptype=TSTEP_CONSTANT )

    CALL mem%Add_var( 'height_oro_min', mem%height_oro_min,                              &
      & hgrid, surface,                                                      &
      & t_cf('height_oro_min', 'm', 'minimum elevation within grid cell'),         &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
      & prefix, suffix,                                                      &
      & lrestart=.FALSE.,                                                    &
      & initval_r=0.0_wp, isteptype=TSTEP_CONSTANT )

    CALL mem%Add_var( 'height_dyn', mem%height_dyn,                          &
      & hgrid, surface,                                                      &
      & t_cf('height_dyn', 'm', 'dynamic minimum height in each grid cell'), &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
      & prefix, suffix,                                                      &
      & output_level=BASIC,                                                  &
      & initval_r=0.0_wp )

    CALL mem%Add_var( 'alpha', mem%alpha,                                    &
      & hgrid, surface,                                                      &
      & t_cf('alpha', 'none', 'parameter of height-volume relations'),       &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
      & prefix, suffix,                                                      &
      & lrestart=.FALSE.,                                                    &
      & initval_r=0.0_wp, isteptype=TSTEP_CONSTANT )

    CALL mem%Add_var( 'beta', mem%beta,                                      &
      & hgrid, surface,                                                      &
      & t_cf('beta', 'none', 'parameter of area-volume relations'),          &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
      & prefix, suffix,                                                      &
      & lrestart=.FALSE.,                                                    &
      & initval_r=0.0_wp, isteptype=TSTEP_CONSTANT )

    CALL mem%Add_var( 'basin_mask', mem%basin_mask,                          &
      & hgrid, surface,                                                      &
      & t_cf('basin_mask', 'none', 'lake basin grid point flagged'),         &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
      & prefix, suffix,                                                      &
      & lrestart=.FALSE.,                                                    &
      & initval_r=0.0_wp, isteptype=TSTEP_CONSTANT )

     CALL mem%Add_var( 'lake_inflow', mem%lake_inflow,                       &
      & hgrid, surface,                                                      &
      & t_cf('lake_inflow', 'm3/s', 'lake inflow from runoff'),              &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
      & prefix, suffix,                                                      &
      & output_level=BASIC,                                                  &
      & initval_r=0.0_wp )

     CALL mem%Add_var( 'lake_drain', mem%lake_drain,                         &
      & hgrid, surface,                                                      &
      & t_cf('lake_drain', 'm3/s', 'lake inflow from runoff and discharge'), &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
      & prefix, suffix,                                                      &
      & output_level=BASIC,                                                  &
      & initval_r=0.0_wp )

    CALL mem%Add_var( 'outflow_resid', mem%outflow_resid,                          &
      & hgrid, surface,                                                                       &
      & t_cf('outflow_resid', 'm3/s', 'grid cell outflow residual'),                    &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                    &
      & prefix, suffix,                                                                 &
      & initval_r=0.0_wp )

    CALL mem%Add_var( 'overlflow_res', mem%overlflow_res,                           &
      & hgrid, hd_o,                                                                     &
      & t_cf('overlflow_res', '', 'content of the overland flow reservoir'),             &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),               &
      & prefix, suffix,                                                      &
      & initval_r=0.0_wp )

    CALL mem%Add_var( 'baseflow_res', mem%baseflow_res,     &
      & hgrid, hd_b,                                                                     &
      & t_cf('baseflow_res', '', 'content of the baseflow reservoir'),                   &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),               &
      & prefix, suffix,                                                      &
      & initval_r=0.0_wp )

    CALL mem%Add_var( 'riverflow_res', mem%riverflow_res,   &
      & hgrid, hd_r,                                                                     &
      & t_cf('riverflow_res', '', 'content of the riverflow reservoir'),                 &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),               &
      & prefix, suffix,                                                      &
      & initval_r=0.0_wp )

    CALL mem%Add_var( 'water_budget', mem%water_budget,                             &
      & hgrid, surface,                                                                  &
      & t_cf('hd_water_budget', 'm3', 'HD water budget'),                                &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),               &
      & prefix, suffix,                                                                  &
      & lrestart=.TRUE.,                                                                 &
      & initval_r=0.0_wp )

    CALL mem%Add_var( 'water_budget_old', mem%water_budget_old,                     &
      & hgrid, surface,                                                                  &
      & t_cf('hd_water_budget_old', 'm3', 'HD water budget (prev. time step)'),          &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),               &
      & prefix, suffix,                                                                  &
      & lrestart=.FALSE.,                                                                &
      & initval_r=0.0_wp )

    CALL mem%Add_var( 'water_flux', mem%water_flux,                                 &
      & hgrid, surface,                                                                  &
      & t_cf('hd_water_flux', 'm3/s', 'HD water flux'),                                  &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),               &
      & prefix, suffix,                                                                  &
      & lrestart=.FALSE.,                                                                &
      & initval_r=0.0_wp )

    CALL mem%Add_var( 'water_error', mem%water_error,                               &
      & hgrid, surface,                                                                  &
      & t_cf('hd_water_error', 'm3', 'HD water imbalance during time step'),             &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),               &
      & prefix, suffix,                                                                  &
      & lrestart=.FALSE.,                                                                &
      & initval_r=0.0_wp )

  END SUBROUTINE Init_hd_memory

#endif
END MODULE mo_hd_memory_class
