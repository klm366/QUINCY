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

  USE mo_kind,                   ONLY: wp
  USE mo_exception,              ONLY: message
  USE mo_jsb_class,              ONLY: Get_model

  USE mo_jsb_memory_class,       ONLY: t_jsb_memory
  USE mo_jsb_var_class,          ONLY: t_jsb_var_real2d, t_jsb_var_real3d, t_jsb_var_real1d

  USE mo_slm_varlist,            ONLY: slm_add_var   ! replace by mo_jsb_varlist later  *slmtemp



  IMPLICIT NONE
  PRIVATE
  PUBLIC :: t_pheno_memory, max_no_of_vars

  INTEGER, PARAMETER :: max_no_of_vars = 20   ! default JSBACH4


  !-----------------------------------------------------------------------------------------------------
  !> Type definition for pheno memory
  !! 
  !-----------------------------------------------------------------------------------------------------
  TYPE, EXTENDS(t_jsb_memory) :: t_pheno_memory

    TYPE(t_jsb_var_real1d)    :: &
                              test_pheno  !< test !


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
!    USE mo_jsb_io,            ONLY: grib_bits, t_cf, t_grib1, t_grib2, TSTEP_CONSTANT, tables
!    USE mo_jsb_grid_class,    ONLY: t_jsb_grid, t_jsb_vgrid
!    USE mo_jsb_grid,          ONLY: Get_grid, Get_vgrid
    USE mo_base_memory,         ONLY: t_model

    ! 0.1.InOut
    CLASS(t_pheno_memory), INTENT(inout), TARGET :: mem   !> pheno memory
    CHARACTER(len=*),    INTENT(in)    :: prefix          !> ?
    CHARACTER(len=*),    INTENT(in)    :: suffix          !> ?
    INTEGER,             INTENT(in)    :: lct_ids !(:)    !> in jsbach4: vector of land cover type IDs
    INTEGER,             INTENT(in)    :: model_id        !> model ID model\%id
    ! 0.2 Local
    TYPE(t_model), POINTER :: model
!    TYPE(t_jsb_model), POINTER :: model
!    TYPE(t_jsb_grid),  POINTER :: hgrid                        ! Horizontal grid
!    TYPE(t_jsb_vgrid), POINTER :: surface, canopy_layer        ! Vertical grid
!    INTEGER :: table

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':Init_pheno_memory'

    model => Get_model(model_id)

!    table = tables(1)

!    hgrid   => Get_grid(mem%grid_id)
!    surface => Get_vgrid('surface')
!    canopy_layer => Get_vgrid('canopy_layer')

    IF (.NOT. ALLOCATED(mem%test_pheno%ptr)) THEN

      !> add variables to pheno memory pheno\%VAR
      !! 
      !! slm_add_var -> slm_add_var_list_element_r1d    *slmtemp 1D -> 2D
      !!  
      CALL slm_add_var(mem, 'test_pheno', mem%test_pheno, &
            initval_r=0.0_wp)     ! initvalue 


    END IF

  END SUBROUTINE Init_pheno_memory

END MODULE mo_pheno_memory_class
