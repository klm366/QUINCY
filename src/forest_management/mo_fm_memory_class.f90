!>
!! @brief Contains the memory class for the fm process.
!!
!! @author
!!  Julia Nabel, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Julia Nabel (2018-03-19)
!!
MODULE mo_fm_memory_class
#ifndef __NO_JSBACH__

  USE mo_kind, ONLY: wp
  USE mo_util, ONLY: One_of

  USE mo_jsb_model_class,        ONLY: t_jsb_model
  USE mo_jsb_class,              ONLY: Get_model
  USE mo_jsb_memory_class,       ONLY: t_jsb_memory
  USE mo_jsb_lct_class,          ONLY: VEG_TYPE
  USE mo_jsb_var_class,          ONLY: t_jsb_var_real2d, t_jsb_var_real3d

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: t_fm_memory, max_no_of_vars

  INTEGER, PARAMETER :: max_no_of_vars = 5

  !> Type definition for memory
  TYPE, EXTENDS(t_jsb_memory) :: t_fm_memory
    TYPE(t_jsb_var_real3d) ::   &
          to_harvest                  !Array to keep harvest per pft on the grid cell for now 
  CONTAINS
    PROCEDURE :: Init => Init_fm_memory
  END TYPE t_fm_memory


  CHARACTER(len=*), PARAMETER :: modname = 'mo_fm_memory_class'

CONTAINS

  SUBROUTINE Init_fm_memory(mem, prefix, suffix, lct_ids, model_id)

    USE mo_jsb_varlist,       ONLY: jsb_add_var, BASIC, MEDIUM, FULL
    USE mo_jsb_io,            ONLY: grib_bits, t_cf, t_grib1, t_grib2, TSTEP_CONSTANT, tables
    USE mo_jsb_grid_class,    ONLY: t_jsb_grid, t_jsb_vgrid
    USE mo_jsb_grid,          ONLY: Get_grid, Get_vgrid

    CLASS(t_fm_memory), INTENT(inout), TARGET :: mem
    CHARACTER(len=*),     INTENT(in)    :: prefix
    CHARACTER(len=*),     INTENT(in)    :: suffix
    INTEGER,              INTENT(in)    :: lct_ids(:)
    INTEGER,              INTENT(in)    :: model_id

    TYPE(t_jsb_model), POINTER :: model
    TYPE(t_jsb_grid),  POINTER :: hgrid          ! Horizontal grid
    TYPE(t_jsb_vgrid), POINTER :: soil_depth_water    ! Vertical grid
!    TYPE(t_jsb_vgrid), POINTER :: cohort_grid    ! Vertical grid
    INTEGER :: table

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':Init_fm_memory'

    model => Get_model(model_id)

    table = tables(1)

    hgrid   => Get_grid(mem%grid_id)

    !TODO-JN: this is a temporary solution - when replacing reading with e.g. Stiigs solution, here an own grid should be defined if still required 
    soil_depth_water   => Get_vgrid('soil_depth_water')

    IF (.NOT. ASSOCIATED(mem%to_harvest%ptr)) THEN

      CALL jsb_add_var( mem, 'to_harvest', mem%to_harvest,                   &
!        & hgrid, cohort_grid,                                                                    &
        & hgrid, soil_depth_water, &
        & t_cf('to_harvest', '', 'Whatever to Harvest.'), &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                 &
        & prefix, suffix,                                                                    &
        & loutput = .TRUE.,lrestart=.FALSE., initval_r=0.0_wp ) 

    END IF


  END SUBROUTINE Init_fm_memory

#endif
END MODULE mo_fm_memory_class
