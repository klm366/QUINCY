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
  USE mo_jsb_var_class,    ONLY: t_jsb_var_real1d, t_jsb_var_real2d, t_jsb_var_real3d
  USE mo_jsb_varlist,      ONLY: BASIC !, MEDIUM, FULL

  ! Use of prcesses in this module
  !dsl4jsb_Use_processes HD_

  ! Use process configurations
  !dsl4jsb_Use_config(HD_)

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: t_hd_memory, max_no_of_vars

  INTEGER, PARAMETER :: max_no_of_vars = 60

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
      nsplit          , &  !< number of channels river splits into 
      outflow_runoff  , &  !< grid cell outflow due to runoff
      outflow_drainage, &  !< grid cell outflow due to drainage
      outflow_rivers  , &  !< grid cell outflow due to rivers
      outflow_outlet  , &  !< outlet outflow from lake [m3/s]
      inflow_outlet   , &  !< lake outlet outflow to HD [m3/s]
      discharge_lake  , &  !< lake discharge (storage changes) [m3/s]
      runoff_lake     , &  !< lake inflow from over lake runoff (P-E) [m3/s]
      riverflow_lake  , &  !< river discharge into lake [m3/s]
      overlflow_lake  , &  !< overland flow into lake [m3/s]
      baseflow_lake   , &  !< base flow discharge into lake [m3/s]
      inflow_lake     , &  !< lake inflow from neighbor grid cell [m3/s]
      outflow_lake    , &  !< lake outflow to neighbor grid cell [m3/s]
      overflow_lake   , &  !< lake overflow outflow [m3/s]
      volume_lake_old , &  !< lake water volume from previous time step [m3]
      volume_lake     , &  !< lake water volume storage [m3]
      area_lake       , &  !< lake area of endorheic lakes [m^2]
      height_lake     , &  !< lake height of endorheic lakes [m]
      fract_lake      , &  !< lake fraction of endorheic lakes [-]
      volume_lake_max , &  !< maximum lake water volume of endorheic lakes [m^3]
      height_oro_min  , &  !< minimum elevation within ICON grid cell [m]
      height_dyn      , &  !< minimum dynamic height within ICON grid cell [m]
      basin_mask      , &  !< endorheic lake basin grid points flagged with [0;>0] 
      outlets         , &  !< endorheic lake outlet grid points flagged with [0;>0] 
      w_river_res     , &  !< river reservoir water content [m] 
      w_overl_res     , &  !< overland reservoir water content [m] 
      w_river_res_old , &  !< river reservoir water content [m]
      w_overl_res_old , &  !< overland reservoir water content [m]
      fract_river     , &  !< river reservoir wetskin fraction [-]
      fract_overl     , &  !< overland reservoir wetskin fraction [-]
      volume_river_seas_mean, &  !< month mean river water volume storage [m]
      volume_overl_seas_mean, &  !< month mean overland water volume storage [m]
      local_budget    , &  !< total HD reservoir content
      local_fluxes    , &  !< grid cell in- and outflow
      local_wbal_error, &  !< grid cell water balance error
      water_budget    , &
      water_budget_change, &
      outflow_resid,    &  !< residual outflow that's not routed to discharge - only used for debugging (see mo_hd_interface.f90)
      outflow_count        !< count of outflows from a given cell - only used for debugging (see mo_hd_interface.f90)

    TYPE(t_jsb_var_real3d) :: &
      nidx_upstream,          &  !< 1st dim. of local index of grid cell upstream (nproma)
      bidx_upstream,          &  !< 2nd dim. of local index of grid cell upstream (nblks)
      overlflow_res,          &  !< content of overland flow reservoir cascade
      baseflow_res,           &  !< content of baseflow reservoir cascade
      riverflow_res,          &  !< content of groundflow reservoir cascade
      volumes_lake,           &  !< lake volumes from lookup table [m3]
      areas_lake,             &  !< lake areas from looup table [m2]
      heights_lake,           &  !< lake heights from loopup table [m]
      volume_river_seas,      &  !< lake water volume storage [m]
      volume_overl_seas,      &  !< lake water volume storage [m]      
      outflow_min                !< outflow minimum height [m]

    ! Diagnostic global land means/sums e.g. for monitoring (only available with ICON)
    TYPE(t_jsb_var_real1d) :: &
      water_error_gsum           !< Total HD water imbalance during time step       [m3]

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
    TYPE(t_jsb_vgrid), POINTER :: surface, hd_o, hd_b, hd_r, dl_l, dg_n, day_1y  ! Vertical grids
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
    dl_l => get_vgrid('dynamic_lake_level')
    dg_n => get_vgrid('edge_neighbors')
    day_1y => get_vgrid('days_in_year')

    CALL mem%Add_var( 'hd_mask', mem%hd_mask,                                                   &
      & hgrid, surface,                                                                         &
      & t_cf('hd_mask', '', 'HD_ mask: -1 ocean, 0 ocean inflow, 1 land, 2 internal drainage'), &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                      & 
      & prefix, suffix,                                                                         &
      & lrestart=.FALSE.,                                                                       &
      & initval_r=0.0_wp, isteptype=TSTEP_CONSTANT )

    CALL mem%Add_var( 'coast_ocean', mem%coast_ocean,                                    &
      & hgrid, surface,                                                                  &
      & t_cf('coast_ocean', '', 'coastal grid cells of the ocean '),                     &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),               &
      & prefix, suffix,                                                                  &
      & lrestart=.FALSE.,                                                                &
      & initval_r=0.0_wp, isteptype=TSTEP_CONSTANT )

    CALL mem%Add_var( prefix//'ret_overlflow'//suffix, mem%ret_overlflow,                &
      & hgrid, surface,                                                                  &
      & t_cf('ret_overlflow', 'day', 'overland flow retention constant'),                &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),               &
      & prefix, suffix,                                                                  &
      & lrestart=.FALSE.,                                                                &
      & initval_r=0.0_wp, isteptype=TSTEP_CONSTANT )

    CALL mem%Add_var( 'ret_baseflow', mem%ret_baseflow,                                  &
      & hgrid, surface,                                                                  &
      & t_cf('ret_baseflow', 'day', 'baseflow retention constant'),                      &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),               &
      & prefix, suffix,                                                                  &
      & lrestart=.FALSE.,                                                                &
      & initval_r=0.0_wp, isteptype=TSTEP_CONSTANT )

    CALL mem%Add_var( 'ret_riverflow', mem%ret_riverflow,                                &
      & hgrid, surface,                                                                  &
      & t_cf('ret_riverflow', 'day', 'riverflow retention constant'),                    &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),               &
      & prefix, suffix,                                                                  &
      & lrestart=.FALSE.,                                                                &
      & initval_r=0.0_wp, isteptype=TSTEP_CONSTANT )

    CALL mem%Add_var( 'nres_overlflow', mem%nres_overlflow,                              &
      & hgrid, surface,                                                                  &
      & t_cf('nres_overlflow', '', 'number of overland flow reservoirs'),                &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),               &
      & prefix, suffix,                                                                  &
      & lrestart=.FALSE.,                                                                &
      & initval_r=0.0_wp, isteptype=TSTEP_CONSTANT )

    CALL mem%Add_var( 'nres_baseflow', mem%nres_baseflow,                                &
      & hgrid, surface,                                                                  &
      & t_cf('nres_baseflow', '', 'number of baseflow reservoirs'),                      &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),               &
      & prefix, suffix,                                                                  &
      & lrestart=.FALSE.,                                                                &
      & initval_r=0.0_wp, isteptype=TSTEP_CONSTANT )

    CALL mem%Add_var( 'nres_riverflow', mem%nres_riverflow,                              &
      & hgrid, surface,                                                                  &
      & t_cf('nres_riverflow', '', 'number of riverflow reservoirs'),                    &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),               &
      & prefix, suffix,                                                                  &
      & lrestart=.FALSE.,                                                                &
      & initval_r=0.0_wp, isteptype=TSTEP_CONSTANT )

    CALL mem%Add_var( 'nsplit', mem%nsplit,                                              &
      & hgrid, surface,                                                                  &
      & t_cf('nsplit', '', 'number of channels river splits into'),                      &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),               &
      & prefix, suffix,                                                                  &
      & lrestart=.FALSE.,                                                                &
      & initval_r=0.0_wp, isteptype=TSTEP_CONSTANT )

    CALL mem%Add_var( 'outflow_runoff', mem%outflow_runoff,                              &
      & hgrid, surface,                                                                  &
      & t_cf('outflow_runoff', 'm3/s', 'grid cell outflow due to runoff'),               &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),               &
      & prefix, suffix,                                                                  &
      & output_level=BASIC,                                                              &
      & initval_r=0.0_wp )

    CALL mem%Add_var( 'outflow_drainage', mem%outflow_drainage,                          &
      & hgrid, surface,                                                                  &
      & t_cf('outflow_drainage', 'm3/s', 'grid cell outflow due to drainage'),           &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),               &
      & prefix, suffix,                                                                  &
      & output_level=BASIC,                                                              &
      & initval_r=0.0_wp )

    CALL mem%Add_var( 'outflow_rivers', mem%outflow_rivers,                              &
      & hgrid, surface,                                                                  &
      & t_cf('outflow_rivers', 'm3/s', 'grid cell outflow due to rivers'),               &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),               &
      & prefix, suffix,                                                                  &
      & output_level=BASIC,                                                              &
      & initval_r=0.0_wp )

     CALL mem%Add_var( 'inflow_outlet', mem%inflow_outlet,                   &
      & hgrid, surface,                                                      &
      & t_cf('inflow_outlet', 'm3/s', 'outlet outflow to HD'),                 &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
      & prefix, suffix,                                                      &
      & output_level=BASIC,                                                  &
      & initval_r=0.0_wp )

     CALL mem%Add_var( 'outflow_outlet', mem%outflow_outlet,                 &
      & hgrid, surface,                                                      &
      & t_cf('outflow_outlet', 'm3/s', 'lake outflow to outlet'),              &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
      & prefix, suffix,                                                      &
      & output_level=BASIC,                                                  &
      & initval_r=0.0_wp )

     CALL mem%Add_var( 'discharge_lake', mem%discharge_lake,                 &
      & hgrid, surface,                                                      &
      & t_cf('discharge_lake', 'm3/s', 'lake discharge'),                    &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
      & prefix, suffix,                                                      &
      & initval_r=0.0_wp )

     CALL mem%Add_var( 'runoff_lake', mem%runoff_lake,                       &
      & hgrid, surface,                                                      &
      & t_cf('runoff_lake', 'm3/s', 'lake inflow from runoff'),              &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
      & prefix, suffix,                                                      &
      & output_level=BASIC,                                                  &
      & initval_r=0.0_wp )

     CALL mem%Add_var( 'riverflow_lake', mem%riverflow_lake,                 &
      & hgrid, surface,                                                      &
      & t_cf('riverflow_lake', 'm3/s', 'river discharge into lake'),         &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
      & prefix, suffix,                                                      &
      & output_level=BASIC,                                                  &
      & initval_r=0.0_wp )

     CALL mem%Add_var( 'overlflow_lake', mem%overlflow_lake,                 &
      & hgrid, surface,                                                      &
      & t_cf('overlflow_lake', 'm3/s', 'overland flow into lake'),           &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
      & prefix, suffix,                                                      &
      & output_level=BASIC,                                                  &
      & initval_r=0.0_wp )

     CALL mem%Add_var( 'baseflow_lake', mem%baseflow_lake,                   &
      & hgrid, surface,                                                      &
      & t_cf('baseflow_lake', 'm3/s', 'baseflow into lake'),                 &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
      & prefix, suffix,                                                      &
      & output_level=BASIC,                                                  &
      & initval_r=0.0_wp )

     CALL mem%Add_var( 'inflow_lake', mem%inflow_lake,                       &
      & hgrid, surface,                                                      &
      & t_cf('inflow_lake', 'm3/s', 'lake inflow to neighbor grid cell'),    &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
      & prefix, suffix,                                                      &
      & output_level=BASIC,                                                  &
      & initval_r=0.0_wp )

     CALL mem%Add_var( 'outflow_lake', mem%outflow_lake,                     &
      & hgrid, surface,                                                      &
      & t_cf('outflow_lake', 'm3/s', 'lake outflow to neighbor grid cell'),  &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
      & prefix, suffix,                                                      &
      & output_level=BASIC,                                                  &
      & initval_r=0.0_wp )

     CALL mem%Add_var( 'overflow_lake', mem%overflow_lake,                   &
      & hgrid, surface,                                                      &
      & t_cf('overflow_lake', 'm3/s', 'lake outflow to neighbor grid cell'), &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
      & prefix, suffix,                                                      &
      & output_level=BASIC,                                                  &
      & initval_r=0.0_wp )

    CALL mem%Add_var( 'volume_lake_old', mem%volume_lake_old,                &
      & hgrid, surface,                                                      &
      & t_cf('volume_lake_old', 'm3', 'lake volume from previous time step'),&
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
      & prefix, suffix,                                                      &
      & output_level=BASIC,                                                  &
      & initval_r=0.0_wp )

    CALL mem%Add_var( 'volume_lake', mem%volume_lake,                        &
      & hgrid, surface,                                                      &
      & t_cf('volume_lake', 'm3', 'grid cell lake volume'),                  &
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

      CALL mem%Add_var( 'height_lake', mem%height_lake,                      &
      & hgrid, surface,                                                      &
      & t_cf('height_lake', 'm', 'grid cell lake height'),                   &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
      & prefix, suffix,                                                      & 
      & output_level=BASIC,                                                  &
      & initval_r=0.0_wp )

     CALL mem%Add_var( 'fract_lake', mem%fract_lake,                         &
      & hgrid, surface,                                                      &
      & t_cf('fract_lake', 'none', 'grid cell lake fraction'),               &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
      & prefix, suffix,                                                      &
      & output_level=BASIC,                                                  &
      & initval_r=0.0_wp )

    CALL mem%Add_var( 'volume_lake_max', mem%volume_lake_max,                &
      & hgrid, surface,                                                      &
      & t_cf('volume_lake_max', 'm3', 'grid cell maxium lake volume'),       &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
      & prefix, suffix,                                                      &
      & lrestart=.FALSE.,                                                    &
      & initval_r=0.0_wp, isteptype=TSTEP_CONSTANT )

    CALL mem%Add_var( 'height_oro_min', mem%height_oro_min,                  &
      & hgrid, surface,                                                      &
      & t_cf('height_oro_min', 'm', 'minimum elevation within grid cell'),   &
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

    CALL mem%Add_var( 'basin_mask', mem%basin_mask,                          &
      & hgrid, surface,                                                      &
      & t_cf('basin_mask', 'none', 'lake basin grid point flagged'),         &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
      & prefix, suffix,                                                      &
      & lrestart=.FALSE.,                                                    &
      & initval_r=0.0_wp, isteptype=TSTEP_CONSTANT )

    CALL mem%Add_var( 'outlets', mem%outlets,                                &
      & hgrid, surface,                                                      &
      & t_cf('outlets', 'none', 'lake outlet grid point flagged'),           &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
      & prefix, suffix,                                                      &
      & lrestart=.FALSE.,                                                    &
      & initval_r=0.0_wp, isteptype=TSTEP_CONSTANT )

    CALL mem%Add_var( 'w_river_res', mem%w_river_res,                        &
      & hgrid, surface,                                                      &
      & t_cf('w_river_res', 'm', 'river reservoir water content'),           &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
      & prefix, suffix,                                                      &
      & lrestart=.TRUE.,                                                     &
      & output_level=BASIC,                                                  &
      & initval_r=0.0_wp)

    CALL mem%Add_var( 'w_overl_res', mem%w_overl_res,                        &
      & hgrid, surface,                                                      &
      & t_cf('w_overl_res', 'm', 'overland reservoir water content'),        &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
      & prefix, suffix,                                                      &
      & lrestart=.TRUE.,                                                     &
      & output_level=BASIC,                                                  &
      & initval_r=0.0_wp)

    CALL mem%Add_var( 'w_river_res_old', mem%w_river_res_old,                &
      & hgrid, surface,                                                      &
      & t_cf('w_river_res_old', 'm', 'river reservoir water content'),       &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
      & prefix, suffix,                                                      &
      & lrestart=.TRUE.,                                                     &
      & initval_r=0.0_wp)

    CALL mem%Add_var( 'w_overl_res_old', mem%w_overl_res_old,                &
      & hgrid, surface,                                                      &
      & t_cf('w_overl_res_old', 'm', 'overland reservoir water content'),    &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
      & prefix, suffix,                                                      &
      & lrestart=.TRUE.,                                                     &
      & initval_r=0.0_wp)

    CALL mem%Add_var( 'fract_river', mem%fract_river,                        &
      & hgrid, surface,                                                      &
      & t_cf('fract_river', 'none', 'river reservoir wetskin fraction'),     &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
      & prefix, suffix,                                                      &
      & lrestart=.TRUE.,                                                     &
      & output_level=BASIC,                                                  &
      & initval_r=0.0_wp)

    CALL mem%Add_var( 'fract_overl', mem%fract_overl,                        &
      & hgrid, surface,                                                      &
      & t_cf('fract_overl', 'none', 'overland reservoir wetskin fraction'),  &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
      & prefix, suffix,                                                      &
      & lrestart=.TRUE.,                                                     &
      & output_level=BASIC,                                                  &
      & initval_r=0.0_wp)

    CALL mem%Add_var( 'volume_river_seas_mean', mem%volume_river_seas_mean,  &
      & hgrid, surface,                                                      &
      & t_cf('volume_river_seas_mean', 'm', 'seas mean river volume'),      &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
      & prefix, suffix,                                                      &
      & output_level=BASIC,                                                  &
      & initval_r=0.0_wp)

    CALL mem%Add_var( 'volume_overl_seas_mean', mem%volume_overl_seas_mean,  &
      & hgrid, surface,                                                      &
      & t_cf('volume_overl_seas_mean', 'm', 'seas mean overland volume'),   &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
      & prefix, suffix,                                                      &
      & output_level=BASIC,                                                  &
      & initval_r=0.0_wp)

    CALL mem%Add_var( 'local_budget', mem%local_budget,                                  &
      & hgrid, surface,                                                                  &
      & t_cf('local_budget', 'm3/s', 'Total HD reservoir content'),                      &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),               &
      & prefix, suffix,                                                                  &
      & initval_r=0.0_wp )

    CALL mem%Add_var( 'local_fluxes', mem%local_fluxes,                                  &
      & hgrid, surface,                                                                  &
      & t_cf('local_fluxes', 'm3/s', 'Fluxes for local water balance check'),            &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),               &
      & prefix, suffix,                                                                  &
      & lrestart=.FALSE.,                                                                &
      & initval_r=0.0_wp )

    CALL mem%Add_var( 'local_wbal_error', mem%local_wbal_error,                          &
      & hgrid, surface,                                                                  &
      & t_cf('local_wbal_error', 'm3', 'Local water balance error'),                     &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),               &
      & prefix, suffix,                                                                  &
      & loutput=.TRUE., output_level=BASIC,                                              &
      & lrestart=.FALSE.,                                                                &
      & initval_r=0.0_wp )

    CALL mem%Add_var( 'water_budget', mem%water_budget,                                  &
      & hgrid, surface,                                                                  &
      & t_cf('hd_water_budget', 'm3', 'HD water budget'),                                &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),               &
      & prefix, suffix,                                                                  &
      & lrestart=.TRUE.,                                                                 &
      & initval_r=0.0_wp )

    CALL mem%Add_var( 'water_budget_change', mem%water_budget_change,                    &
      & hgrid, surface,                                                                  &
      & t_cf('hd_water_budget_change', 'm3/timestep', 'HD water budget change within time step'), &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),               &
      & prefix, suffix,                                                                  &
      & lrestart=.FALSE.,                                                                &
      & initval_r=0.0_wp )

    CALL mem%Add_var( 'outflow_resid', mem%outflow_resid,                                &
      & hgrid, surface,                                                                  &
      & t_cf('outflow_resid', 'm3/s', 'grid cell outflow residual: should be zero.'),    &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),               &
      & prefix, suffix,                                                                  &
      & initval_r=0.0_wp )

    CALL mem%Add_var( 'outflow_count', mem%outflow_count,                                &
      & hgrid, surface,                                                                  &
      & t_cf('outflow_count', 'm3/s', 'count of outflows from grid cell'),               &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),               &
      & prefix, suffix,                                                                  &
      & initval_r=0.0_wp )

    CALL mem%Add_var( 'overlflow_res', mem%overlflow_res,                                &
      & hgrid, hd_o,                                                                     &
      & t_cf('overlflow_res', '', 'content of the overland flow reservoir'),             &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),               &
      & prefix, suffix,                                                                  &
      & output_level=BASIC,                                                              &
      & initval_r=0.0_wp )

    CALL mem%Add_var( 'baseflow_res', mem%baseflow_res,                                  &
      & hgrid, hd_b,                                                                     &
      & t_cf('baseflow_res', '', 'content of the baseflow reservoir'),                   &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),               &
      & prefix, suffix,                                                                  &
      & output_level=BASIC,                                                              &
      & initval_r=0.0_wp )

    CALL mem%Add_var( 'riverflow_res', mem%riverflow_res,                                &
      & hgrid, hd_r,                                                                     &
      & t_cf('riverflow_res', '', 'content of the riverflow reservoir'),                 &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),               &
      & prefix, suffix,                                                                  &
      & output_level=BASIC,                                                              &
      & initval_r=0.0_wp )

    CALL mem%Add_var( 'volumes_lake', mem%volumes_lake,                      &
      & hgrid, dl_l,                                                         &
      & t_cf('volumes_lake', 'm3', 'lookup lake volumes'),                   &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
      & prefix, suffix,                                                      &
      & initval_r=0.0_wp, isteptype=TSTEP_CONSTANT )

    CALL mem%Add_var( 'areas_lake', mem%areas_lake,                          &
      & hgrid, dl_l,                                                         &
      & t_cf('areas_lake', 'm2', 'lookup lake areas'),                       &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &           
      & prefix, suffix,                                                      &           
      & initval_r=0.0_wp, isteptype=TSTEP_CONSTANT )                                     

    CALL mem%Add_var( 'heights_lake', mem%heights_lake,                      &
      & hgrid, dl_l,                                                         &
      & t_cf('heights_lake', 'm', 'lookup lake heights'),                   &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &           
      & prefix, suffix,                                                      &           
      & initval_r=0.0_wp, isteptype=TSTEP_CONSTANT )  

    CALL mem%Add_var( 'volume_river_seas', mem%volume_river_seas,            &
      & hgrid, day_1y,                                                       &
      & t_cf('volume_river_seas', 'm', 'lake minum volume'),                &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
      & prefix, suffix,                                                      &
      & lrestart=.TRUE.,                                                     &
      & initval_r=0.0_wp)

    CALL mem%Add_var( 'volume_overl_seas', mem%volume_overl_seas,            &
      & hgrid, day_1y,                                                       &
      & t_cf('volume_overl_seas', 'm', 'lake minum volume'),                &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
      & prefix, suffix,                                                      &
      & lrestart=.TRUE.,                                                     &
      & initval_r=0.0_wp)                      

    CALL mem%Add_var( 'outflow_min', mem%outflow_min,                        &
      & hgrid, dg_n,                                                         &
      & t_cf('outflow_min', 'm', 'outflow minimum height'),                  &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
      & prefix, suffix,                                                      &
      & initval_r=0.0_wp, isteptype=TSTEP_CONSTANT )

#ifdef __ICON__
    ! Diagnostic 1d global land variables for experiment monitoring
    !       (1d stream variables are not supported with echam)
    !
    IF ( TRIM(suffix) == 'box' ) THEN
      CALL mem%Add_var( 'water_error_gsum', mem%water_error_gsum,                          &
        & hgrid, surface,                                                                  &
        & t_cf('hd_water_error_gsum', 'm3', 'Total HD water imbalance during time step'),  &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),               &
        & prefix, suffix,                                                                  &
        & lrestart=.FALSE.,                                                                &
        & initval_r=0.0_wp )
    ENDIF
#endif  

  END SUBROUTINE Init_hd_memory

#endif
END MODULE mo_hd_memory_class
