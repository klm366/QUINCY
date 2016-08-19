!>
!! @brief Contains the memory class for the assimilation process.
!!
!! @author
!!  Jan Engel, MPI-BGC Jena (based on the JSBACH4 version from Rainer Schneck, MPI-M Hamburg)
!!
!! @par Revision History
!! First version by Jan Engel (20160815)
!!
MODULE mo_assimi_memory_class

  USE mo_kind, ONLY: wp

  USE mo_jsb_memory_class,       ONLY: t_jsb_memory
  USE mo_jsb_var_class,          ONLY: t_jsb_var_real2d, t_jsb_var_real3d, t_jsb_var_real1d

  USE mo_slm_varlist,            ONLY: slm_add_var   ! replace by mo_jsb_varlist later



  IMPLICIT NONE
  PRIVATE
  PUBLIC :: t_assimi_memory, max_no_of_vars

  INTEGER, PARAMETER :: max_no_of_vars = 40   ! default JSBACH4

  !-----------------------------------------------------------------------------------------------------
  !> Type definition for assimi memory
  !! 
  !-----------------------------------------------------------------------------------------------------
  TYPE, EXTENDS(t_jsb_memory) :: t_assimi_memory

!    TYPE(t_jsb_var_real2d)  :: &
!                                gpp, &      ! plant%flux%gpp
!                                gpp_c13     ! plant%flux%gpp_c13

    TYPE(t_jsb_var_real1d)    :: &
                                gpp, &
                                gpp_c13









  CONTAINS
    PROCEDURE :: Init => Init_assimi_memory
  END TYPE t_assimi_memory


  CHARACTER(len=*), PARAMETER :: modname = 'mo_assimi_memory_class'




CONTAINS


  !-----------------------------------------------------------------------------------------------------
  !> initialize memory (variables) for the process: assimilation
  !!
  !!  
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE Init_assimi_memory(mem, prefix, suffix, lct_ids, model_id, ics, ice, iblk)

!    USE mo_jsb_varlist,       ONLY: t_jsb_varlist, jsb_add_var
!    USE mo_jsb_io,            ONLY: grib_bits, t_cf, t_grib1, t_grib2, TSTEP_CONSTANT, tables
!    USE mo_jsb_grid_class,    ONLY: t_jsb_grid, t_jsb_vgrid
!    USE mo_jsb_grid,          ONLY: Get_grid, Get_vgrid

    CLASS(t_assimi_memory), INTENT(inout), TARGET :: mem
    CHARACTER(len=*),    INTENT(in)    :: prefix
    CHARACTER(len=*),    INTENT(in)    :: suffix
    INTEGER,             INTENT(in)    :: lct_ids !(:)
    INTEGER,             INTENT(in)    :: model_id
    INTEGER,             INTENT(in)    :: ics   ! model%options -> for SLM compatibility
    INTEGER,             INTENT(in)    :: ice   ! model%options -> for SLM compatibility
    INTEGER,             INTENT(in)    :: iblk  ! model%options -> for SLM compatibility

!    TYPE(t_jsb_model), POINTER :: model
!    TYPE(t_jsb_grid),  POINTER :: hgrid                        ! Horizontal grid
!    TYPE(t_jsb_vgrid), POINTER :: surface, canopy_layer        ! Vertical grid
!    INTEGER :: table

!    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':Init_assimi_memory'

!    model => Get_model(model_id)

!    table = tables(1)

!    hgrid   => Get_grid(mem%grid_id)
!    surface => Get_vgrid('surface')
!    canopy_layer => Get_vgrid('canopy_layer')

    IF (.NOT. ASSOCIATED(mem%gpp%ptr)) THEN

      CALL slm_add_var(mem, 'gross_assimilation', mem%gpp, &   ! i.e. gross photosynthesis
            ics, ice, iblk, &    ! model%options -> for SLM compatibility
!        & hgrid, surface,                                                                              &
!        & t_cf('gross_assimilation', 'µmol CO2 m-2 s-1', 'Gross photosynthesis'),    &
!        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                           &
!        & prefix, suffix,                                                                                &
!        & lrestart=.FALSE., 
          initval_r=0.0_wp) ! initval taken from JSBACH3

      CALL slm_add_var(mem, 'gross_assimilation_c13', mem%gpp_c13, &   ! i.e. gross photosynthesis
            ics, ice, iblk, &    ! model%options -> for SLM compatibility
!        & hgrid, surface,                                                                              &
!        & t_cf('gross_assimilation', 'µmol CO2 m-2 s-1', 'Gross photosynthesis'),    &
!        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                           &
!        & prefix, suffix,                                                                                &
!        & lrestart=.FALSE., 
          initval_r=0.0_wp) ! initval taken from JSBACH3

    END IF

  END SUBROUTINE Init_assimi_memory

END MODULE mo_assimi_memory_class
