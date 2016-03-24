!>
!! @brief Contains the memory class for the surface energy balance process.
!!
!! @author
!!  Reiner Schnur, MPI-M Hamburg
!!
!! @par Revision History
!! First version in mo_srf_memory_class.f90 by Reiner Schnur (2013-04-04)
!! Separated out the surface energy balance part by Reiner Schnur (2016-03-21)
!!
MODULE mo_seb_memory_class
#ifndef __NO_JSBACH__

  USE mo_kind, ONLY: wp

  USE mo_jsb_memory_class, ONLY: t_jsb_memory
  USE mo_jsb_var_class,    ONLY: t_jsb_var_real2d, t_jsb_var_real3d

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: t_seb_memory, max_no_of_vars

  INTEGER, PARAMETER :: max_no_of_vars = 40

  !> Type definition for memory
  TYPE, EXTENDS(t_jsb_memory) :: t_seb_memory
    TYPE(t_jsb_var_real2d) :: &
      & t,                    & !< Surface temperature                       [K]
      & t_old,                & !< Surface temperature at previous timestep  [K]
      & s,                    & !< Surface dry static energy                 [m2 s-2]
      & heat_cap,             & !< Heat capacity of surface layer            [J m-2 K-1]
      & latent_hflx,          & !< Latent heat flux                          [W m-2]
      & sensible_hflx           !< Sensible heat flux                        [W m-2]
  CONTAINS
    PROCEDURE :: Init => Init_seb_memory
  END TYPE t_seb_memory

  CHARACTER(len=*), PARAMETER :: modname = 'mo_seb_memory_class'

CONTAINS

  SUBROUTINE Init_seb_memory(mem, prefix, suffix)

    USE mo_jsb_varlist,       ONLY: t_jsb_varlist, jsb_add_var
    USE mo_jsb_io,            ONLY: grib_bits, t_cf, t_grib1, t_grib2, TSTEP_CONSTANT, tables
    USE mo_jsb_grid_class,    ONLY: t_jsb_grid, t_jsb_vgrid
    USE mo_jsb_grid,          ONLY: Get_grid, Get_vgrid

    CLASS(t_seb_memory), INTENT(inout) :: mem
    CHARACTER(len=*),    INTENT(in)    :: prefix
    CHARACTER(len=*),    INTENT(in)    :: suffix

    TYPE(t_jsb_grid),  POINTER :: hgrid                        ! Horizontal grid
    TYPE(t_jsb_vgrid), POINTER :: surface                      ! Vertical grid
    INTEGER :: table

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':Init_seb_memory'

    table = tables(1)

    hgrid   => Get_grid(mem%grid_id)
    surface => Get_vgrid('surface')

    CALL jsb_add_var( mem, 't', mem%t,                                               &
      & hgrid, surface,                                                              &
      & t_cf('sfc_temp', 'K', 'surface temperature'),                                &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),           &
      & prefix, suffix,                                                              &
      & initval_r=280.0_wp )

    CALL jsb_add_var( mem, 't_old', mem%t_old,                                       &
      & hgrid, surface,                                                              &
      & t_cf('sfc_temp_old', 'K', 'surface temperature at previous timestep'),       &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),           &
      & prefix, suffix,                                                              &
      & initval_r=280.0_wp )

    CALL jsb_add_var( mem, 's', mem%s,                                               &
      & hgrid, surface,                                                              &
      & t_cf('sfc_dry_static_energy', 'm2 s-2', 'surface dry static energy'),        &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),           &
      & prefix, suffix,                                                              &
      & initval_r=2.9E5_wp )

    CALL jsb_add_var( mem, 'heat_cap', mem%heat_cap,                                 &
      & hgrid, surface,                                                              &
      & t_cf('sfc_heat_cap', 'J m-2 K-1)', 'surface layer heat capacity'),           &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),           &
      & prefix, suffix,                                                              &
      & initval_r=0.0_wp )

    CALL jsb_add_var( mem, 'latent_hflx', mem%latent_hflx,                           &
      & hgrid, surface,                                                              &
      & t_cf('latent_hflx', 'W m-2', 'latent heat flux at surface'),                 &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),           &
      & prefix, suffix,                                                              &
      & initval_r=0.0_wp )

    CALL jsb_add_var( mem, 'sensible_hflx', mem%sensible_hflx,                       &
      & hgrid, surface,                                                              &
      & t_cf('sensible_hflx', 'W m-2', 'sensible heat flux at surface'),             &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),           &
      & prefix, suffix,                                                              &
      & initval_r=280.0_wp )

  END SUBROUTINE Init_seb_memory

#endif
END MODULE mo_seb_memory_class
