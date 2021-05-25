!>
!! @brief Contains the memory class for anthropogenic land cover change (ALCC).
!!
!! @author
!!  Rainer Schneck, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Rainer Schneck (2018-12-03)


MODULE mo_alcc_memory_class
#ifndef __NO_JSBACH__

  USE mo_kind, ONLY: wp
  USE mo_util, ONLY: One_of

  USE mo_jsb_model_class,        ONLY: t_jsb_model
  USE mo_jsb_class,              ONLY: Get_model
  USE mo_jsb_memory_class,       ONLY: t_jsb_memory
  USE mo_jsb_lct_class,          ONLY: VEG_TYPE, LAND_TYPE
  USE mo_jsb_var_class,          ONLY: t_jsb_var_real2d, t_jsb_var_real3d
  USE mo_jsb_varlist,            ONLY: BASIC, FULL !, MEDIUM
  !USE mo_jsb_physical_constants, ONLY: tmelt

  ! Use of processes in this module
  dsl4jsb_Use_processes ALCC_

  ! Use process configurations
  dsl4jsb_Use_config(ALCC_)

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: t_alcc_memory, max_no_of_vars

  INTEGER, PARAMETER :: max_no_of_vars = 40

  TYPE, EXTENDS(t_jsb_memory) :: t_alcc_memory

    TYPE(t_jsb_var_real2d) :: C_pools_sum_check
    TYPE(t_jsb_var_real3d)           :: &
      & fract_pfts_currentyear_maps,    &   !< cover fractions for the pfts. These fractions are read in yearly from a map file. [-]
      & fract_pfts_before_alccstep          !< cover fractions for the pfts. These fractions are those before alcc changed them within
                                            !< this timestep. [-]

  CONTAINS
    PROCEDURE :: Init => Init_alcc_memory
  END TYPE t_alcc_memory

  CHARACTER(len=*), PARAMETER :: modname = 'mo_alcc_memory_class'

CONTAINS


  SUBROUTINE Init_alcc_memory(mem, prefix, suffix, lct_ids, model_id, l_fixed_fractions)

    ! Declarations
    USE mo_jsb_io,             ONLY: grib_bits, t_cf, t_grib1, t_grib2, tables
    USE mo_jsb_grid_class,     ONLY: t_jsb_grid, t_jsb_vgrid
    USE mo_jsb_grid,           ONLY: Get_grid, Get_vgrid
    !USE mo_jsb_model_usecases, ONLY: npft ! R: Note, this leads to circular dependency

    ! Arguments
    CLASS(t_alcc_memory), INTENT(inout), TARGET :: mem
    CHARACTER(len=*),      INTENT(in)             :: prefix
    CHARACTER(len=*),      INTENT(in)             :: suffix
    INTEGER,               INTENT(in)             :: lct_ids(:)
    INTEGER,               INTENT(in)             :: model_id
    LOGICAL, OPTIONAL,     INTENT(in)             :: l_fixed_fractions

    ! Local variables
    dsl4jsb_Def_config(ALCC_)

    TYPE(t_jsb_model), POINTER  :: model
    TYPE(t_jsb_grid),  POINTER  :: hgrid                        ! Horizontal grid
    TYPE(t_jsb_vgrid), POINTER  :: surface, pft_vgrid           ! Vertical grid
    INTEGER                     :: table, npft

    CHARACTER(len=*), PARAMETER :: routine = modname//':Init_alcc_memory'

    ! ---------------------------
    ! Go

    model => Get_model(model_id)

    table = tables(1)

    ! Set pointers to variables in the "grid"
    hgrid   => Get_grid(mem%grid_id)
    surface => Get_vgrid('surface')
    pft_vgrid    => Get_vgrid('pfts')

!    dsl4jsb_Get_config(ALCC_)

    CALL mem%Add_var('C_pools_sum_check', mem%C_pools_sum_check,             &
      & hgrid, surface,                                                      &
      & t_cf('C_pools_sum_check in alcc', '[-]', ''),                        &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
      & prefix, suffix,                                                      &
      & output_level=BASIC,                                                  &
      & loutput=.TRUE., lrestart=.FALSE., initval_r=0._wp)
      
    CALL mem%Add_var('fract_pfts_currentyear_maps', mem%fract_pfts_currentyear_maps,                 &
      & hgrid, pft_vgrid,                                                                            &
      & t_cf('Cover fractions for the pfts, which are read in yearly from a map file.', '[-]', ''),  &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                           &
      & prefix, suffix,                                                                              &
      & output_level=BASIC,                                                                          &
      & loutput=.TRUE., lrestart=.FALSE., initval_r=-0.1_wp)
      
    CALL mem%Add_var('fract_pfts_before_alccstep', mem%fract_pfts_before_alccstep,                      &
      & hgrid, pft_vgrid,                                                                               &
      & t_cf('Former cover fractions for the pfts, before read in yearly from a map file.', '[-]', ''), &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                              &
      & prefix, suffix,                                                                                 &
      & output_level=BASIC,                                                                             &
      & loutput=.TRUE., lrestart=.FALSE., initval_r=-0.1_wp)
      
  END SUBROUTINE Init_alcc_memory


#endif
END MODULE mo_alcc_memory_class
