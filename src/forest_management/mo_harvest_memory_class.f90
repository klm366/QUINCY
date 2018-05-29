!>
!! @brief Contains the memory class for the forest harvest process.
!!
!! @author
!!  Julia Nabel, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Julia Nabel (2018-03-19)
!!
MODULE mo_harv_memory_class
#ifndef __NO_JSBACH__

  USE mo_kind, ONLY: wp
  USE mo_util, ONLY: One_of

  USE mo_jsb_model_class,        ONLY: t_jsb_model
  USE mo_jsb_class,              ONLY: Get_model
  USE mo_jsb_memory_class,       ONLY: t_jsb_memory
  USE mo_jsb_lct_class,          ONLY: VEG_TYPE
  USE mo_jsb_var_class,          ONLY: t_jsb_var_real2d, t_jsb_var_real3d

  dsl4jsb_Use_processes FAGE_

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: t_harv_memory, max_no_of_vars

  INTEGER, PARAMETER :: max_no_of_vars = 5

  !> Type definition for memory
  TYPE, EXTENDS(t_jsb_memory) :: t_harv_memory
   TYPE(t_jsb_var_real2d) ::   &
            area_to_harvest,   &   ! are that should be harvested
            harvested_area,    &   ! are that was harvested
            cconservation_harvest  ! field for c conservation check for this process
  CONTAINS
    PROCEDURE :: Init => Init_harv_memory
  END TYPE t_harv_memory


  CHARACTER(len=*), PARAMETER :: modname = 'mo_harv_memory_class'

CONTAINS

  SUBROUTINE Init_harv_memory(mem, prefix, suffix, lct_ids, model_id)

    USE mo_jsb_varlist,       ONLY: jsb_add_var, BASIC, MEDIUM, FULL
    USE mo_jsb_io,            ONLY: grib_bits, t_cf, t_grib1, t_grib2, TSTEP_CONSTANT, tables
    USE mo_jsb_grid_class,    ONLY: t_jsb_grid, t_jsb_vgrid
    USE mo_jsb_grid,          ONLY: Get_grid, Get_vgrid

    CLASS(t_harv_memory), INTENT(inout), TARGET :: mem
    CHARACTER(len=*),     INTENT(in)    :: prefix
    CHARACTER(len=*),     INTENT(in)    :: suffix
    INTEGER,              INTENT(in)    :: lct_ids(:)
    INTEGER,              INTENT(in)    :: model_id

    TYPE(t_jsb_model), POINTER :: model
    TYPE(t_jsb_grid),  POINTER :: hgrid          ! Horizontal grid
    TYPE(t_jsb_vgrid), POINTER :: surface        ! Vertical grid
    INTEGER :: table

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':Init_harv_memory'

    model => Get_model(model_id)

    table = tables(1)

    hgrid   => Get_grid(mem%grid_id)
    surface => Get_vgrid('surface')

    IF ( (     One_of(VEG_TYPE,  lct_ids(:)) > 0 &
        &  )                                       &
        & .AND. .NOT. ASSOCIATED(mem%area_to_harvest%ptr)) THEN

      CALL jsb_add_var( mem, 'area_to_harvest', mem%area_to_harvest,                   &
        & hgrid, surface,                                                                    &
        & t_cf('area_to_harvest', '', 'Area to be harvested.'), &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                 &
        & prefix, suffix,                                                                    &
        & loutput = .TRUE.,lrestart=.FALSE., initval_r=0.0_wp ) 

      CALL jsb_add_var( mem, 'harvested_area', mem%harvested_area,                   &
        & hgrid, surface,                                                                    &
        & t_cf('harvested_area', '', 'Actually harvested area.'), &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                 &
        & prefix, suffix,                                                                    &
        & loutput = .TRUE.,lrestart=.FALSE., initval_r=0.0_wp ) 

      IF (.NOT. model%processes(FAGE_)%p%config%active) THEN
        ! C conservation on the tile can only be tested for harvest, 
        !    if the harvested area stays on the tile, i.e. if fage is not active
        CALL jsb_add_var( mem, 'cconservation_harvest', mem%cconservation_harvest,    &
          & hgrid, surface,                                                           &
          & t_cf('cconservation_harvest', 'mol(C) m-2(canopy)', &
            & 'Test for c conservation on relocation for harvest.'), &
          & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),        &
          & prefix, suffix,                                                           &
          & loutput = .TRUE.,lrestart=.FALSE., initval_r=0.0_wp ) 
      END IF

    END IF

  END SUBROUTINE Init_harv_memory

#endif
END MODULE mo_harv_memory_class
