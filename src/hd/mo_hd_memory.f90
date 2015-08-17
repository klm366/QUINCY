!>
!! @brief Contains methods for the soil memory.
!!
!! @author
!!  Veronika Gayler, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Veronika Gayler (2015-07-17)
!!
MODULE mo_hd_memory
!#ifndef __NO_JSBACH__
#if !defined(__NO_JSBACH__) && !defined(__NO_JSBACH_HD__)

  USE mo_hd_memory_class,   ONLY: t_hd_param, t_hd_state, t_hd_memory
  USE mo_kind,              ONLY: wp
  USE mo_exception,         ONLY: message

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: new_hd_memory

  CHARACTER(len=*), PARAMETER :: modname = 'mo_hd_memory'

CONTAINS

  FUNCTION new_hd_memory(varlist_name, prefix, suffix, tile) RESULT(mem)

    USE mo_jsb_varlist,       ONLY: t_var_list, jsb_get_var_list, jsb_new_var_list, jsb_add_var
    USE mo_jsb_io,            ONLY: grib_bits, nc_bits, t_cf, t_grib1, t_grib2, &
                                    FILETYPE_NC2, TSTEP_CONSTANT, ZAXIS_SURFACE, ZAXIS_DEPTH_BELOW_LAND, tables
    USE mo_jsb_grid_class,    ONLY: t_jsb_grid, t_jsb_vgrid, new_vgrid
    USE mo_jsb_tile_class,    ONLY: t_jsb_tile
    USE mo_jsb_model_class,   ONLY: t_jsb_model, get_vgrid, register_vgrid
    USE mo_jsb_base,          ONLY: get_model

    CHARACTER(len=*),           INTENT(in) :: varlist_name
    CHARACTER(len=*),           INTENT(in) :: prefix
    CHARACTER(len=*),           INTENT(in) :: suffix
#ifdef HAVE_F2003
    TYPE(t_jsb_tile), POINTER , INTENT(in) :: tile
#else
    TYPE(t_jsb_tile), POINTER              :: tile
#endif
    TYPE(t_hd_memory), POINTER             :: mem

    TYPE(t_hd_memory), POINTER  :: new
    TYPE(t_var_list),  POINTER  :: this_list
    TYPE(t_jsb_model), POINTER  :: model

    INTEGER :: nsoil_energy, nsoil_water
    TYPE(t_jsb_grid),  POINTER :: hgrid                      ! Horizontal grid
    TYPE(t_jsb_vgrid), POINTER :: surface, hd_o, hd_b, hd_r  ! Vertical grids
    INTEGER :: table

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':new_hd_memory'

    model => get_model(tile%owner_model_id)

    ALLOCATE(new)

    table = tables(1)

    !this_list => get_var_list(TRIM(varlist_name))
    CALL jsb_get_var_list(this_list, TRIM(varlist_name))
    IF (.NOT. ASSOCIATED(this_list)) THEN
      !this_list => new_var_list(TRIM(varlist_name), patch_id=patch_id)
      ALLOCATE(this_list)
      CALL jsb_new_var_list(this_list, TRIM(varlist_name),              &
                            lrestart=.TRUE., restart_type=FILETYPE_NC2, &
                            patch_id=model%grid%patch%id, table=table)
    END IF

    hgrid        => model%grid

    surface => get_vgrid(model, 'surface')
    IF (.NOT. ASSOCIATED(surface)) THEN
      surface => new_vgrid('surface', ZAXIS_SURFACE, 1)
      CALL register_vgrid(model, surface)
    END IF

    !vg   wie kann ich die level dynamisch definieren??
    hd_o => get_vgrid(model, 'hd_nres_overlflow')
    IF (.NOT. ASSOCIATED(hd_o)) THEN
      hd_o  => new_vgrid('hd_nres_overlflow', ZAXIS_DEPTH_BELOW_LAND, 1, levels=(/1._wp/), units='')
      CALL register_vgrid(model, hd_o)
    END IF
    hd_b => get_vgrid(model, 'hd_nres_baseflow')
    IF (.NOT. ASSOCIATED(hd_b)) THEN
      hd_b  => new_vgrid('hd_nres_baseflow', ZAXIS_DEPTH_BELOW_LAND, 1, levels=(/1._wp/), units='')
      CALL register_vgrid(model, hd_b)
    END IF
    hd_r => get_vgrid(model, 'hd_nres_riverflow')
    IF (.NOT. ASSOCIATED(hd_r)) THEN
      hd_r  => new_vgrid('hd_nres_riverflow', ZAXIS_DEPTH_BELOW_LAND, 5, levels=(/1._wp,2._wp,3._wp,4._wp,5._wp/), units='')
      CALL register_vgrid(model, hd_r)
    END IF

    IF (ASSOCIATED(tile%parent_tile)) THEN
       new%hd_param => tile%parent_tile%mem%hd%hd_param
    ELSE
       new%hd_param => new_param()
    END IF

    new%hd_state => new_state()

    mem => new

  CONTAINS

    FUNCTION new_param() RESULT(param)

      TYPE(t_hd_param), POINTER :: param

      ALLOCATE(param)

      CALL jsb_add_var( this_list, prefix//'ret_overlflow'//suffix, param%ret_overlflow,     &
        & hgrid, surface,                                                                    &
        & t_cf('ret_overlflow', 'day', 'overland flow retention constant', nc_bits),         &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                 &
        & initval_r=0.0_wp, isteptype=TSTEP_CONSTANT )

      CALL jsb_add_var( this_list, prefix//'ret_baseflow'//suffix, param%ret_baseflow,       &
        & hgrid, surface,                                                                    &
        & t_cf('ret_baseflow', 'day', 'baseflow retention constant', nc_bits),               &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                 &
        & initval_r=0.0_wp, isteptype=TSTEP_CONSTANT )

      CALL jsb_add_var( this_list, prefix//'ret_riverflow'//suffix, param%ret_riverflow,     &
        & hgrid, surface,                                                                    &
        & t_cf('ret_riverflow', 'day', 'riverflow retention constant', nc_bits),             &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                 &
        & initval_r=0.0_wp, isteptype=TSTEP_CONSTANT )

      CALL jsb_add_var( this_list, prefix//'nres_overlflow'//suffix, param%nres_overlflow,   &
        & hgrid, surface,                                                                    &
        & t_cf('nres_overlflow', '', 'number of overland flow reservoirs', nc_bits),         &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                 &
        & initval_r=0.0_wp, isteptype=TSTEP_CONSTANT )

      CALL jsb_add_var( this_list, prefix//'nres_baseflow'//suffix, param%nres_baseflow,     &
        & hgrid, surface,                                                                    &
        & t_cf('nres_baseflow', '', 'number of baseflow reservoirs', nc_bits),               &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                 &
        & initval_r=0.0_wp, isteptype=TSTEP_CONSTANT )

      CALL jsb_add_var( this_list, prefix//'nres_riverflow'//suffix, param%nres_riverflow,   &
        & hgrid, surface,                                                                    &
        & t_cf('nres_riverflow', '', 'number of riverflow reservoirs', nc_bits),             &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                 &
        & initval_r=0.0_wp, isteptype=TSTEP_CONSTANT )

    END FUNCTION new_param

    FUNCTION new_state() RESULT(state)

      TYPE(t_hd_state), POINTER :: state

      ALLOCATE(state)

      CALL jsb_add_var( this_list, prefix//'overlflow_res'//suffix, state%overlflow_res,   &
        & hgrid, hd_o,                                                                     &
        & t_cf('overlflow_res', '', 'content of the overland flow reservoir', nc_bits),    &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),               &
        & initval_r=0.0_wp )

      CALL jsb_add_var( this_list, prefix//'baseflow_res'//suffix, state%baseflow_res,     &
        & hgrid, hd_b,                                                                     &
        & t_cf('baseflow_res', '', 'content of the baseflow reservoir', nc_bits),          &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),               &
        & initval_r=0.0_wp )

      CALL jsb_add_var( this_list, prefix//'riverflow_res'//suffix, state%riverflow_res,   &
        & hgrid, hd_r,                                                                     &
        & t_cf('riverflow_res', '', 'content of the riverflow reservoir', nc_bits),        &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),               &
        & initval_r=0.0_wp )

    END FUNCTION new_state

  END FUNCTION new_hd_memory

#endif
END MODULE mo_hd_memory
