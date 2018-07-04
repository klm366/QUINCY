!>
!! @brief Contains the memory class for the phenology process.
!!
!! @author
!!  Jan Engel, MPI-BGC Jena (based on the JSBACH4 version from Rainer Schneck, MPI-M Hamburg)
!!
!! @par Revision History
!! First version by Jan Engel (20161214)
!!
MODULE mo_pheno_memory_class
#ifndef __NO_JSBACH__

  USE mo_kind,                   ONLY: wp, l_false
  USE mo_exception,              ONLY: message
  USE mo_jsb_class,              ONLY: Get_model
  USE mo_jsb_memory_class,       ONLY: t_jsb_memory
  USE mo_jsb_var_class,          ONLY: t_jsb_var_real2d, t_jsb_var_real3d              ! *slmjsb4
  USE mo_slm_varlist,            ONLY: jsb_add_var   ! replace by mo_jsb_varlist later  *slmtemp



  IMPLICIT NONE
  PRIVATE
  PUBLIC :: t_pheno_memory, max_no_of_vars

  INTEGER, PARAMETER :: max_no_of_vars = 20   ! default JSBACH4


  !-----------------------------------------------------------------------------------------------------
  !> Type definition for pheno memory
  !! 
  !-----------------------------------------------------------------------------------------------------
  TYPE, EXTENDS(t_jsb_memory) :: t_pheno_memory


    ! used as logical: 0 = FALSE / 1 = TRUE
    TYPE(t_jsb_var_real2d)          :: &
                                    l_growing_season        !< logical to identify whether the plant is in the growing season or not


    TYPE(t_jsb_var_real2d)          :: &
                                    gdd                 , & !< growing degree days above 5°C
                                    nd_dormance         , & !< number of days since last growing season 
                                    lai_max                 !< maximum lai - from site specific lai data


    !! test variable introduced by Steffen to have a test for the pheno_interface aggregate & integrate routines !!
    TYPE(t_jsb_var_real2d)          :: testvar

  CONTAINS
    PROCEDURE :: Init => Init_pheno_memory
  END TYPE t_pheno_memory


  CHARACTER(len=*), PARAMETER :: modname = 'mo_pheno_memory_class'


CONTAINS


  !-----------------------------------------------------------------------------------------------------
  !> initialize memory (variables) for the process: phenology
  !!
  !!  
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE Init_pheno_memory(mem, prefix, suffix, lct_ids, model_id)

!    USE mo_jsb_varlist,       ONLY: t_jsb_varlist, jsb_add_var
    USE mo_jsb_io,            ONLY: grib_bits, t_cf, t_grib1, t_grib2   ! , TSTEP_CONSTANT, tables
    USE mo_jsb_grid_class,    ONLY: t_jsb_grid, t_jsb_vgrid
!    USE mo_jsb_grid,          ONLY: Get_grid, Get_vgrid
    USE mo_jsb_model_class,   ONLY: t_jsb_model

    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    CLASS(t_pheno_memory), INTENT(inout), TARGET :: mem             !> pheno memory
    CHARACTER(len=*),      INTENT(in)            :: prefix          !> ?
    CHARACTER(len=*),      INTENT(in)            :: suffix          !> ?
    INTEGER,               INTENT(in)            :: lct_ids !(:)    !> in jsbach4: vector of land cover type IDs
    INTEGER,               INTENT(in)            :: model_id        !> model ID model\%id
    ! ---------------------------
    ! 0.2 Local
    TYPE(t_jsb_model), POINTER   :: model
    TYPE(t_jsb_grid),  POINTER   :: hgrid                        ! Horizontal grid
!    TYPE(t_jsb_vgrid), POINTER   :: canopy_layer                 ! Vertical grid
    TYPE(t_jsb_vgrid), POINTER   :: surface                      ! Vertical grid
    INTEGER                      :: table
    CHARACTER(len=*),  PARAMETER :: routine = TRIM(modname)//':Init_pheno_memory'

    model => get_model(model_id)

!    table = tables(1)
    table = 1

!    hgrid   => Get_grid(mem%grid_id)
!    surface => Get_vgrid('surface')
!    canopy_layer => Get_vgrid('canopy_layer')

    IF (.NOT. ASSOCIATED(mem%l_growing_season%ptr)) THEN

    !! test variable introduced by Steffen to have a test for the pheno_interface aggregate & integrate routines !!
    CALL jsb_add_var( mem, 'testvar', mem%testvar,      &
          hgrid, surface,                               &
          t_cf('testvariable', '', ''),                 &
          t_grib1(table, 255, grib_bits),               &
          t_grib2(255, 255, 255, grib_bits),            &
          prefix, suffix,                               &
          initval_r=0.0_wp )

      !> add variables to pheno memory pheno\%VAR
      !! 
      !! jsb_add_var -> jsb_add_var_list_element_r1d    *slmtemp 1D -> 2D
      !!  
      CALL jsb_add_var(mem, 'l_growing_season', mem%l_growing_season, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=l_false) ! initvalue

      CALL jsb_add_var(mem, 'gdd', mem%gdd, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'nd_dormance', mem%nd_dormance, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue

      CALL jsb_add_var(mem, 'lai_max', mem%lai_max, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue

    END IF

  END SUBROUTINE Init_pheno_memory

#endif
END MODULE mo_pheno_memory_class
