!>
!! @brief Contains the memory class for the forest-age process.
!!
!! @author
!!  Julia Nabel, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Julia Nabel (2017-11-21)
!!

MODULE mo_fage_memory_class
#ifndef __NO_JSBACH__

  USE mo_kind, ONLY: wp
  USE mo_util, ONLY: One_of

  USE mo_jsb_model_class,        ONLY: t_jsb_model
  USE mo_jsb_class,              ONLY: Get_model
  USE mo_jsb_memory_class,       ONLY: t_jsb_memory
  USE mo_jsb_lct_class,          ONLY: VEG_TYPE, LAND_TYPE
  USE mo_jsb_var_class,          ONLY: t_jsb_var_real2d, t_jsb_var_real3d

  ! Use of processes in this module
  dsl4jsb_Use_processes FAGE_, DISTURB_

  ! Use process configurations
  dsl4jsb_Use_config(FAGE_)
  dsl4jsb_Use_config(DISTURB_)
  

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: t_fage_memory, max_no_of_vars

  INTEGER, PARAMETER :: max_no_of_vars = 40

  TYPE, EXTENDS(t_jsb_memory) :: t_fage_memory
    TYPE(t_jsb_var_real2d) :: &
      & mean_age                    ! Mean age of child cohorts
    TYPE(t_jsb_var_real3d) :: &
      & fract_per_age,        &     ! Area fraction having a certain age - age = index third dim
      & cohort_max_age,       &     ! maximum age of the child cohorts 
      & disturbed_area              ! If run with disturbances: disturbed area 
                                    ! -> TODO-JN: this only makes sense with fixed cohort boundaries!? 
  CONTAINS
    PROCEDURE :: Init => Init_fage_memory
  END TYPE t_fage_memory

  CHARACTER(len=*), PARAMETER :: modname = 'mo_fage_memory_class'

CONTAINS


  SUBROUTINE Init_fage_memory(mem, prefix, suffix, lct_ids, model_id)

    ! Declarations
    USE mo_jsb_varlist,       ONLY: jsb_add_var
    USE mo_jsb_io,            ONLY: grib_bits, t_cf, t_grib1, t_grib2, tables
    USE mo_jsb_grid_class,    ONLY: t_jsb_grid, t_jsb_vgrid
    USE mo_jsb_grid,          ONLY: Get_grid, Get_vgrid

    ! Arguments
    CLASS(t_fage_memory),  INTENT(inout), TARGET :: mem
    CHARACTER(len=*),      INTENT(in)            :: prefix
    CHARACTER(len=*),      INTENT(in)            :: suffix
    INTEGER,               INTENT(in)            :: lct_ids(:)
    INTEGER,               INTENT(in)            :: model_id

    ! Local variables
    dsl4jsb_Def_config(FAGE_)

    TYPE(t_jsb_model), POINTER :: model
    TYPE(t_jsb_grid),  POINTER :: hgrid                            ! Horizontal grid
    TYPE(t_jsb_vgrid), POINTER :: surface, forest_age, cohort_grid ! Vertical grids
    INTEGER :: table

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':Init_fage_memory'

    ! ---------------------------
    ! Go

    model => Get_model(model_id)

    table = tables(1)

    ! Set pointers to variables in the "grid"
    hgrid       => Get_grid(mem%grid_id)
    surface     => Get_vgrid('surface')
    forest_age   => Get_vgrid('forest_age')
    cohort_grid   => Get_vgrid('cohort_grid')

    dsl4jsb_Get_config(FAGE_)

    IF ( (     One_of(LAND_TYPE, lct_ids(:)) > 0 &
        &   .OR. One_of(VEG_TYPE,  lct_ids(:)) > 0 &
        &  )                                       &
        & .AND. .NOT. ASSOCIATED(mem%mean_age%ptr)) THEN

      CALL jsb_add_var( mem, 'mean_age', mem%mean_age,                         &
        & hgrid, surface,                                                      &
        & t_cf('mean_age', '', ''),                                            &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
        & prefix, suffix,                                                      &
        & loutput = .TRUE.,initval_r=0.0_wp )

      CALL jsb_add_var( mem, 'fract_per_age', mem%fract_per_age,               &
        & hgrid, forest_age,                                                   &
        & t_cf('fract_per_age', 'y', ''),                                      &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
        & prefix, suffix,                                                      &
        & loutput = .TRUE.,initval_r=0.0_wp )

      CALL jsb_add_var( mem, 'cohort_max_age', mem%cohort_max_age,             &
        & hgrid, cohort_grid,                                                  &
        & t_cf('cohort_max_age', 'y', ''),                                     &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
        & prefix, suffix,                                                      &
        & initval_r=0.0_wp )

      IF (model%processes(DISTURB_)%p%config%active) THEN

        CALL jsb_add_var( mem, 'disturbed_area', mem%disturbed_area,             &
          & hgrid, cohort_grid,                                                  &
          & t_cf('disturbed_area', 'y', ''),                                     &
          & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
          & prefix, suffix,                                                      &
          & loutput = .TRUE., initval_r=0.0_wp )

      ENDIF


    END IF

  END SUBROUTINE Init_fage_memory


#endif
END MODULE mo_fage_memory_class