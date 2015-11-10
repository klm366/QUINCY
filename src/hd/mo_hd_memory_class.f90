!>
!! @brief Contains the memory class for the surface.
!!
!! @author
!!  Veronika Gayler, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Veronika Gayler (2015-05-06)
!!
MODULE mo_hd_memory_class
!#ifndef __NO_JSBACH__
#if !defined(__NO_JSBACH__) && !defined(__NO_JSBACH_HD__)

  USE mo_kind, ONLY: wp
  USE mo_jsb_memory_class, ONLY: t_jsb_memory
  USE mo_jsb_var_class,    ONLY: t_jsb_var_real2d, t_jsb_var_real3d

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: t_hd_memory, max_no_of_vars

  INTEGER, PARAMETER :: max_no_of_vars = 40

  TYPE, EXTENDS(t_jsb_memory) :: t_hd_memory

    TYPE(t_jsb_var_real2d) :: &
      hd_mask         , &  !< HD mask
      ret_overlflow   , &  !< retention constant of overland flow
      ret_baseflow    , &  !< retention constant of baseflow
      ret_riverflow   , &  !< retention constant of riverflow
      nres_overlflow  , &  !< number of overland flow reservoirs
      nres_baseflow   , &  !< number of baseflow reservoirs
      nres_riverflow  , &  !< number of riverflow reservoirs
      outflow_runoff  , &  !< grid cell outflow due to runoff
      outflow_drainage, &  !< grid cell outflow due to drainage
      outflow_rivers  , &  !< grid cell outflow due to rivers
      water_budget

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

  SUBROUTINE Init_hd_memory(mem, prefix, suffix)

    USE mo_jsb_varlist,       ONLY: jsb_add_var
    USE mo_jsb_io,            ONLY: grib_bits, nc_bits, t_cf, t_grib1, t_grib2, &
                                    TSTEP_CONSTANT, tables
    USE mo_jsb_grid_class,    ONLY: t_jsb_grid, t_jsb_vgrid
    USE mo_jsb_grid,          ONLY: Get_grid, Get_vgrid

    CLASS(t_hd_memory), INTENT(inout) :: mem
    CHARACTER(len=*),   INTENT(in)    :: prefix
    CHARACTER(len=*),   INTENT(in)    :: suffix

    TYPE(t_jsb_grid),  POINTER :: hgrid                      ! Horizontal grid
    TYPE(t_jsb_vgrid), POINTER :: surface, hd_o, hd_b, hd_r  ! Vertical grids
    INTEGER :: table

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':Init_hd_memory'

    table = tables(1)

    hgrid   => Get_grid(mem%grid_id)
    surface => get_vgrid('surface')

    !vg   wie kann ich die level dynamisch definieren??
    hd_o => Get_vgrid('hd_nres_overlflow')
    hd_b => get_vgrid('hd_nres_baseflow')
    hd_r => get_vgrid('hd_nres_riverflow')

    CALL jsb_add_var( mem, prefix//'hd_mask'//suffix, mem%hd_mask,     &
      & hgrid, surface,                                                                    &
      & t_cf('hd_mask', '', 'HD mask: -1 ocean, 0 ocean inflow, 1 land, 2 internal drainage', nc_bits),         &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                 &
      & initval_r=0.0_wp, isteptype=TSTEP_CONSTANT )

    CALL jsb_add_var( mem, prefix//'ret_overlflow'//suffix, mem%ret_overlflow,     &
      & hgrid, surface,                                                                    &
      & t_cf('ret_overlflow', 'day', 'overland flow retention constant', nc_bits),         &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                 &
      & initval_r=0.0_wp, isteptype=TSTEP_CONSTANT )

    CALL jsb_add_var( mem, prefix//'ret_baseflow'//suffix, mem%ret_baseflow,       &
      & hgrid, surface,                                                                    &
      & t_cf('ret_baseflow', 'day', 'baseflow retention constant', nc_bits),               &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                 &
      & initval_r=0.0_wp, isteptype=TSTEP_CONSTANT )

    CALL jsb_add_var( mem, prefix//'ret_riverflow'//suffix, mem%ret_riverflow,     &
      & hgrid, surface,                                                                    &
      & t_cf('ret_riverflow', 'day', 'riverflow retention constant', nc_bits),             &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                 &
      & initval_r=0.0_wp, isteptype=TSTEP_CONSTANT )

    CALL jsb_add_var( mem, prefix//'nres_overlflow'//suffix, mem%nres_overlflow,   &
      & hgrid, surface,                                                                    &
      & t_cf('nres_overlflow', '', 'number of overland flow reservoirs', nc_bits),         &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                 &
      & initval_r=0.0_wp, isteptype=TSTEP_CONSTANT )

    CALL jsb_add_var( mem, prefix//'nres_baseflow'//suffix, mem%nres_baseflow,     &
      & hgrid, surface,                                                                    &
      & t_cf('nres_baseflow', '', 'number of baseflow reservoirs', nc_bits),               &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                 &
      & initval_r=0.0_wp, isteptype=TSTEP_CONSTANT )

    CALL jsb_add_var( mem, prefix//'nres_riverflow'//suffix, mem%nres_riverflow,           &
      & hgrid, surface,                                                                    &
      & t_cf('nres_riverflow', '', 'number of riverflow reservoirs', nc_bits),             &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                 &
      & initval_r=0.0_wp, isteptype=TSTEP_CONSTANT )

    CALL jsb_add_var( mem, prefix//'outflow_runoff'//suffix, mem%outflow_runoff,           &
      & hgrid, surface,                                                                    &
      & t_cf('outflow_runoff', 'm3/s', 'grid cell outflow due to runoff', nc_bits),        &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                 &
      & initval_r=0.0_wp )

    CALL jsb_add_var( mem, prefix//'outflow_drainage'//suffix, mem%outflow_drainage,       &
      & hgrid, surface,                                                                    &
      & t_cf('outflow_drainage', 'm3/s', 'grid cell outflow due to drainage', nc_bits),    &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                 &
      & initval_r=0.0_wp )

    CALL jsb_add_var( mem, prefix//'outflow_rivers'//suffix, mem%outflow_rivers,              &
      & hgrid, surface,                                                                       &
      & t_cf('outflow_rivers', 'm3/s', 'grid cell outflow due to rivers', nc_bits),           &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                    &
      & initval_r=0.0_wp )

    CALL jsb_add_var( mem, prefix//'overlflow_res'//suffix, mem%overlflow_res,   &
      & hgrid, hd_o,                                                                     &
      & t_cf('overlflow_res', '', 'content of the overland flow reservoir', nc_bits),    &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),               &
      & initval_r=0.0_wp )

    CALL jsb_add_var( mem, prefix//'baseflow_res'//suffix, mem%baseflow_res,     &
      & hgrid, hd_b,                                                                     &
      & t_cf('baseflow_res', '', 'content of the baseflow reservoir', nc_bits),          &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),               &
      & initval_r=0.0_wp )

    CALL jsb_add_var( mem, prefix//'riverflow_res'//suffix, mem%riverflow_res,   &
      & hgrid, hd_r,                                                                     &
      & t_cf('riverflow_res', '', 'content of the riverflow reservoir', nc_bits),        &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),               &
      & initval_r=0.0_wp )

    ! It would be sufficient to store just one grid point
    CALL jsb_add_var( mem, prefix//'water_budget'//suffix, mem%water_budget,             &
      & hgrid, surface,                                                                  &
      & t_cf('water_budget', 'm3', 'global water budget', nc_bits),                      &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),               &
      & initval_r=0.0_wp )

  END SUBROUTINE Init_hd_memory

#endif
END MODULE mo_hd_memory_class
