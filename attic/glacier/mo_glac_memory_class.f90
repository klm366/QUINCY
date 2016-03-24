!>
!! @brief Contains the memory class for the glacier process.
!!
!! @author
!!  Reiner Schnur, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Reiner Schnur (2013-04-04)
!!
MODULE mo_glacier_memory_class
#ifndef __NO_JSBACH__

  USE mo_kind, ONLY: wp
  USE mo_jsb_memory_class, ONLY: t_jsb_memory
  USE mo_jsb_var_class,    ONLY: t_jsb_var_real2d, t_jsb_var_real3d

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: t_glacier_memory, max_no_of_vars
!!$            copy_srf_mem, add_srf_mem

  INTEGER, PARAMETER :: max_no_of_vars = 10

  TYPE, EXTENDS(t_jsb_memory) :: t_glacier_memory
    !
    ! Parameters
    TYPE(t_jsb_var_real2d) :: &
      land         ,    & !< Land fraction: 1: land/lake/glacier, 0: sea
      glacier             !< Glacier fraction
    !
    ! Fractions
    TYPE(t_jsb_var_real2d) :: &
      fract,    &  !< Cover fraction of tile (relative to parent tile)
      fract_max    !< Maximum potential cover fraction for tile (rel. to parent)
  CONTAINS
    PROCEDURE :: Init => Init_glacier_memory
  END TYPE t_glacier_memory

  ! Note: Would be nicer to be able to overload the "=" operator, but that would be more
  ! complicated because we only want to copy parts of the arrays (ics:ice,iblk)
!!$  INTERFACE copy_srf_mem
!!$    MODULE PROCEDURE copy_t_srf_energy_state
!!$  END INTERFACE copy_srf_mem
!!$
!!$  INTERFACE add_srf_mem
!!$    MODULE PROCEDURE add_t_srf_energy_state
!!$  END INTERFACE add_srf_mem

  CHARACTER(len=*), PARAMETER :: modname = 'mo_glacier_memory_class'

CONTAINS

  SUBROUTINE Init_glacier_memory(mem, prefix, suffix, lct_type)

    USE mo_jsb_varlist,       ONLY: t_jsb_varlist, jsb_add_var
    USE mo_jsb_io,            ONLY: grib_bits, t_cf, t_grib1, t_grib2, &
                                    TSTEP_CONSTANT, tables
    USE mo_jsb_grid_class,    ONLY: t_jsb_grid, t_jsb_vgrid
    USE mo_jsb_grid,          ONLY: Get_grid, Get_vgrid

    CLASS(t_glacier_memory), INTENT(inout) :: mem
    CHARACTER(len=*),     INTENT(in)    :: prefix
    CHARACTER(len=*),     INTENT(in)    :: suffix
    INTEGER, OPTIONAL, INTENT(in) :: lct_type

    TYPE(t_jsb_grid),  POINTER :: hgrid                        ! Horizontal grid
    TYPE(t_jsb_vgrid), POINTER :: surface                      ! Vertical grid
    INTEGER :: table

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':Init_glacier_memory'

    table = tables(1)

    hgrid   => Get_grid(mem%grid_id)
    surface => Get_vgrid('surface')

    CALL jsb_add_var( mem, 'land_frac', mem%land,                     &
      & hgrid, surface,                                                               &
      & t_cf('land_mask', '', 'land fraction'),                                       &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),            &
      & prefix, suffix,                                                      &
      & initval_r=0.0_wp, isteptype=TSTEP_CONSTANT )

    CALL jsb_add_var( mem, 'glacier', mem%glacier,           &
      & hgrid, surface,                                                               &
      & t_cf('glacier_mask', '', 'glacier fraction'),                                 &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),            &
      & prefix, suffix,                                                      &
      & initval_r=0.0_wp, isteptype=TSTEP_CONSTANT )

  END SUBROUTINE Init_glacier_memory

!!$  SUBROUTINE copy_t_srf_energy_state(in, out, iblk, ics, ice)
!!$
!!$    TYPE(t_srf_energy_state), INTENT(in)    :: in
!!$    TYPE(t_srf_energy_state), INTENT(inout) :: out
!!$    INTEGER,                  INTENT(in)    :: iblk, ics, ice
!!$
!!$    out%t_srf         (ics:ice,    iblk) = in%t_srf         (ics:ice,    iblk)
!!$    out%t_srf_old     (ics:ice,    iblk) = in%t_srf_old     (ics:ice,    iblk)
!!$    out%t_srf_unfilt  (ics:ice,    iblk) = in%t_srf_unfilt  (ics:ice,    iblk)
!!$    out%t_eff_srf     (ics:ice,    iblk) = in%t_eff_srf     (ics:ice,    iblk)
!!$    out%s_srf         (ics:ice,    iblk) = in%s_srf         (ics:ice,    iblk)
!!$    !out%s_srf_old     (ics:ice,    iblk) = in%s_srf_old     (ics:ice,    iblk)
!!$    out%qsat_srf      (ics:ice,    iblk) = in%qsat_srf      (ics:ice,    iblk)
!!$    !out%qsat_srf_old  (ics:ice,    iblk) = in%qsat_srf_old  (ics:ice,    iblk)
!!$    out%tte_corr      (ics:ice,    iblk) = in%tte_corr      (ics:ice,    iblk)
!!$    out%fact_q_air    (ics:ice,    iblk) = in%fact_q_air    (ics:ice,    iblk)
!!$    out%fact_qsat_srf (ics:ice,    iblk) = in%fact_qsat_srf (ics:ice,    iblk)
!!$    out%fact_qsat_trans_srf (ics:ice,    iblk) = in%fact_qsat_trans_srf (ics:ice,    iblk)
!!$
!!$  END SUBROUTINE copy_t_srf_energy_state
!!$
!!$  SUBROUTINE add_t_srf_energy_state(in1, in2, out, iblk, ics, ice)
!!$
!!$    TYPE(t_srf_energy_state), INTENT(in)    :: in1, in2
!!$    TYPE(t_srf_energy_state), INTENT(inout) :: out
!!$    INTEGER,                  INTENT(in)    :: iblk, ics, ice
!!$
!!$    out%t_srf         (ics:ice,    iblk) = in1%t_srf         (ics:ice,    iblk) + in2%t_srf         (ics:ice,    iblk)
!!$    out%t_srf_old     (ics:ice,    iblk) = in1%t_srf_old     (ics:ice,    iblk) + in2%t_srf_old     (ics:ice,    iblk)
!!$    out%t_srf_unfilt  (ics:ice,    iblk) = in1%t_srf_unfilt  (ics:ice,    iblk) + in2%t_srf_unfilt  (ics:ice,    iblk)
!!$    out%t_eff_srf     (ics:ice,    iblk) = in1%t_eff_srf     (ics:ice,    iblk) + in2%t_eff_srf     (ics:ice,    iblk)
!!$    out%s_srf         (ics:ice,    iblk) = in1%s_srf         (ics:ice,    iblk) + in2%s_srf         (ics:ice,    iblk)
!!$    !out%s_srf_old     (ics:ice,    iblk) = in1%s_srf_old     (ics:ice,    iblk) + in2%s_srf_old     (ics:ice,    iblk)
!!$    out%qsat_srf      (ics:ice,    iblk) = in1%qsat_srf      (ics:ice,    iblk) + in2%qsat_srf      (ics:ice,    iblk)
!!$    !out%qsat_srf_old  (ics:ice,    iblk) = in1%qsat_srf_old  (ics:ice,    iblk) + in2%qsat_srf_old  (ics:ice,    iblk)
!!$    out%tte_corr      (ics:ice,    iblk) = in1%tte_corr      (ics:ice,    iblk) + in2%tte_corr      (ics:ice,    iblk)
!!$    out%fact_q_air    (ics:ice,    iblk) = in1%fact_q_air    (ics:ice,    iblk) + in2%fact_q_air    (ics:ice,    iblk)
!!$    out%fact_qsat_srf (ics:ice,    iblk) = in1%fact_qsat_srf (ics:ice,    iblk) + in2%fact_qsat_srf (ics:ice,    iblk)
!!$    out%fact_qsat_trans_srf (ics:ice,  iblk) = in1%fact_qsat_trans_srf (ics:ice,  iblk) + in2%fact_qsat_trans_srf (ics:ice,  iblk)
!!$
!!$  END SUBROUTINE add_t_srf_energy_state

#endif
END MODULE mo_glacier_memory_class
