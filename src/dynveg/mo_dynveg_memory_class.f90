!>
!! @brief Contains the memory class for dynveglogy.
!!
!! @author
!!  Rainer Schneck, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Rainer Schneck (2018-12-03)


MODULE mo_dynveg_memory_class
#ifndef __NO_JSBACH__

  USE mo_kind, ONLY: wp
  USE mo_util, ONLY: One_of

  USE mo_jsb_model_class,        ONLY: t_jsb_model
  USE mo_jsb_class,              ONLY: Get_model
  USE mo_jsb_memory_class,       ONLY: t_jsb_memory
  USE mo_jsb_lct_class,          ONLY: VEG_TYPE, LAND_TYPE
  USE mo_jsb_var_class,          ONLY: t_jsb_var_real2d !, t_jsb_var_real3d
  USE mo_jsb_varlist,            ONLY: BASIC, FULL !, MEDIUM
  !USE mo_jsb_physical_constants, ONLY: tmelt

  ! Use of processes in this module
  dsl4jsb_Use_processes DYNVEG_

  ! Use process configurations
  dsl4jsb_Use_config(DYNVEG_)

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: t_dynveg_memory, max_no_of_vars

  INTEGER, PARAMETER :: max_no_of_vars = 40

  TYPE, EXTENDS(t_jsb_memory) :: t_dynveg_memory

    TYPE(t_jsb_var_real2d) :: C_pools_sum_check
!   TYPE(t_jsb_var_real3d)           :: &
!     & fract_pfts_currentyear_maps,    &   !< cover fractions for the pfts. These fractions are read in yearly from a map file. [-]
!     & fract_pfts_before_dynvegstep        !< cover fractions for the pfts. These fractions are those before dynveg changed them within
!                                           !< this timestep. [-]

     ! If we do not want this vars in the output or restart file
     REAL(wp), ALLOCATABLE :: fract_pfts_currentyear_maps(:,:,:)  !< cover fractions for the pfts. These fractions are read in yearly from a map file. [-]
     REAL(wp), ALLOCATABLE :: fract_pfts_before_dynvegstep(:,:,:) !< cover fractions for the pfts. These fractions are those before dynveg changed them within
                                                                  !< this timestep. [-]

     ! Additional variable to avoid stream output error
!     TYPE(t_jsb_var_real2d) :: testvar1 


  CONTAINS
    PROCEDURE :: Init => Init_dynveg_memory
  END TYPE t_dynveg_memory

  CHARACTER(len=*), PARAMETER :: modname = 'mo_dynveg_memory_class'

CONTAINS


  SUBROUTINE Init_dynveg_memory(mem, prefix, suffix, lct_ids, model_id)

    ! Declarations
    USE mo_jsb_varlist,        ONLY: jsb_add_var
    USE mo_jsb_io,             ONLY: grib_bits, t_cf, t_grib1, t_grib2, tables
    USE mo_jsb_grid_class,     ONLY: t_jsb_grid, t_jsb_vgrid
    USE mo_jsb_grid,           ONLY: Get_grid, Get_vgrid
    !USE mo_jsb_model_usecases, ONLY: npft ! R: Note, this leads to circular dependency

    ! Arguments
    CLASS(t_dynveg_memory), INTENT(inout), TARGET :: mem
    CHARACTER(len=*),      INTENT(in)             :: prefix
    CHARACTER(len=*),      INTENT(in)             :: suffix
    INTEGER,               INTENT(in)             :: lct_ids(:)
    INTEGER,               INTENT(in)             :: model_id

    ! Local variables
    dsl4jsb_Def_config(DYNVEG_)

    TYPE(t_jsb_model), POINTER  :: model
    TYPE(t_jsb_grid),  POINTER  :: hgrid                        ! Horizontal grid
    TYPE(t_jsb_vgrid), POINTER  :: surface !, pft_vgrid         ! Vertical grid
    INTEGER                     :: table, npft

    CHARACTER(len=*), PARAMETER :: routine = modname//':Init_dynveg_memory'

    ! ---------------------------
    ! Go

    model => Get_model(model_id)

    table = tables(1)

    ! Set pointers to variables in the "grid"
    hgrid   => Get_grid(mem%grid_id)
    surface => Get_vgrid('surface')
!    pft_vgrid    => Get_vgrid('pfts')

!    dsl4jsb_Get_config(DYNVEG_)

    CALL jsb_add_var( mem, 'C_pools_sum_check', mem%C_pools_sum_check,       &
      & hgrid, surface,                                                      &
      & t_cf('C_pools_sum_check in dynveg', '[-]', ''),                      &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),   &
      & prefix, suffix,                                                      &
      & output_level=BASIC,                                                  &
      & loutput=.TRUE., lrestart=.FALSE., initval_r=0._wp)
      
!     CALL jsb_add_var( mem, 'fract_pfts_currentyear_maps', mem%fract_pfts_currentyear_maps,           &
!       & hgrid, pft_vgrid,                                                                            &
!       & t_cf('Cover fractions for the pfts, which are read in yearly from a map file.', '[-]', ''),  &
!       & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                           &
!       & prefix, suffix,                                                                              &
!       & output_level=BASIC,                                                                          &
!       & loutput=.FALSE., lrestart=.TRUE., initval_r=7._wp)
!       
!     CALL jsb_add_var( mem, 'fract_pfts_before_dynvegstep', mem%fract_pfts_before_dynvegstep,            &
!       & hgrid, pft_vgrid,                                                                               &
!       & t_cf('Former cover fractions for the pfts, before read in yearly from a map file.', '[-]', ''), &
!       & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                              &
!       & prefix, suffix,                                                                                 &
!       & output_level=BASIC,                                                                             &
!       & loutput=.FALSE., lrestart=.TRUE., initval_r=7._wp)
!       
    ! If we do not want this vars in the output or restart file
    !npft =tile%no_of_children
    npft=mem%no_of_children
    ALLOCATE(mem%fract_pfts_currentyear_maps  (hgrid%nproma, hgrid%nblks, npft))
    ALLOCATE(mem%fract_pfts_before_dynvegstep (hgrid%nproma, hgrid%nblks, npft))

!     IF (ASSOCIATED(mem%testvar1%ptr)) RETURN
!     CALL jsb_add_var( mem, 'testvar1', mem%testvar1,                           &
!          & hgrid, surface,                                                     &
!          & t_cf('testvar1', 'depends', 'Testvar to avoid stream error.'),      &
!          & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),  &
!          & prefix, suffix,                                                     &
!          & output_level=FULL,                                                  &
!          & lrestart=.FALSE.,                                                   &
!          & initval_r=0.0_wp, l_aggregate_all=.FALSE.  )


  END SUBROUTINE Init_dynveg_memory


#endif
END MODULE mo_dynveg_memory_class
