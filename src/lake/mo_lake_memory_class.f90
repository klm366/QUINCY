!>
!! @brief Contains the memory class for the lake process.
!!
!! @author
!!  Reiner Schnur, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Reiner Schnur (2013-04-04)
!!
MODULE mo_lake_memory_class
#ifndef __NO_JSBACH__

  USE mo_kind, ONLY: wp
  USE mo_jsb_memory_class, ONLY: t_jsb_memory
  USE mo_jsb_var_class,    ONLY: t_jsb_var_real2d, t_jsb_var_real3d

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: t_lake_memory, max_no_of_vars
!!$            copy_srf_mem, add_srf_mem

  INTEGER, PARAMETER :: max_no_of_vars = 10

  TYPE, EXTENDS(t_jsb_memory) :: t_lake_memory
    !
    ! Parameters
    TYPE(t_jsb_var_real2d) :: &
      & land_fract,       & !< Land fraction: 1: land/lake/glacier, 0: sea
      & lake_fract          !< Lake fraction
    !
    ! Fractions
    TYPE(t_jsb_var_real2d) :: &
      & fract,    &  !< Cover fraction of tile (relative to parent tile)
      & fract_max    !< Maximum potential cover fraction for tile (rel. to parent)
    !
    TYPE(t_jsb_var_real2d) :: &
      & t_srf,          & !< Temperature of lake surface
      & hflx_res,       & !! Residual heat flux
      & ice_fract         !! Fraction of ice cover on lake
  CONTAINS
    PROCEDURE :: Init => Init_lake_memory
  END TYPE t_lake_memory

  ! Note: Would be nicer to be able to overload the "=" operator, but that would be more
  ! complicated because we only want to copy parts of the arrays (ics:ice,iblk)
!!$  INTERFACE copy_srf_mem
!!$    MODULE PROCEDURE copy_t_srf_energy_state
!!$  END INTERFACE copy_srf_mem
!!$
!!$  INTERFACE add_srf_mem
!!$    MODULE PROCEDURE add_t_srf_energy_state
!!$  END INTERFACE add_srf_mem

  CHARACTER(len=*), PARAMETER :: modname = 'mo_lake_memory_class'

CONTAINS

  SUBROUTINE Init_lake_memory(mem, prefix, suffix)

    USE mo_jsb_varlist,       ONLY: t_jsb_varlist, jsb_add_var
    USE mo_jsb_io,            ONLY: grib_bits, t_cf, t_grib1, t_grib2, &
                                    TSTEP_CONSTANT, tables
    USE mo_jsb_grid_class,    ONLY: t_jsb_grid, t_jsb_vgrid
    USE mo_jsb_grid,          ONLY: Get_grid, Get_vgrid

    CLASS(t_lake_memory), INTENT(inout) :: mem
    CHARACTER(len=*),     INTENT(in)    :: prefix
    CHARACTER(len=*),     INTENT(in)    :: suffix

    TYPE(t_jsb_grid),  POINTER :: hgrid                        ! Horizontal grid
    TYPE(t_jsb_vgrid), POINTER :: surface                      ! Vertical grid
    INTEGER :: table

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':Init_lake_memory'

    table = tables(1)

    hgrid   => Get_grid(mem%grid_id)
    surface => Get_vgrid('surface')

    CALL jsb_add_var( mem, 'land_frac', mem%land_fract,                     &
      & hgrid, surface,                                                               &
      & t_cf('land_mask', '', 'land fraction'),                                       &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),            &
      & prefix, suffix,                                                      &
      & initval_r=0.0_wp, isteptype=TSTEP_CONSTANT )

    CALL jsb_add_var( mem, 'lake', mem%lake_fract,                          &
      & hgrid, surface,                                                               &
      & t_cf('lake_mask', '', 'lake fraction'),                                       &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),            &
      & prefix, suffix,                                                      &
      & initval_r=0.0_wp, isteptype=TSTEP_CONSTANT )

    CALL jsb_add_var( mem, 't_srf', mem%t_srf,                                   &
      & hgrid, surface,                                                           &
      & t_cf('surface_temperature', 'K', ''),                                              &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                       &
      & prefix, suffix,                                                      &
      & initval_r=280.0_wp )

    CALL jsb_add_var( mem, 'ice_frac', mem%ice_fract,                                   &
      & hgrid, surface,                                                           &
      & t_cf('ice fraction on lake', 'K', ''),                                              &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                       &
      & prefix, suffix,                                                      &
      & initval_r=0.0_wp )

  END SUBROUTINE Init_lake_memory

#endif
END MODULE mo_lake_memory_class
