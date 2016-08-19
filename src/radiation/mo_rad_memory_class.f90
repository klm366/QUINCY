!>
!! @brief Contains the memory class for the radiation process.
!!
!! @author
!!  Jan Engel, MPI-BGC Jena (based on the JSBACH4 version from Rainer Schneck, MPI-M Hamburg)
!!
!! @par Revision History
!! First version by Jan Engel (20160819)
!!
MODULE mo_rad_memory_class

  USE mo_kind, ONLY: wp

  USE mo_jsb_memory_class,       ONLY: t_jsb_memory
  USE mo_jsb_var_class,          ONLY: t_jsb_var_real2d, t_jsb_var_real3d, t_jsb_var_real1d

  USE mo_slm_varlist,            ONLY: slm_add_var   ! replace by mo_jsb_varlist later



  IMPLICIT NONE
  PRIVATE
  PUBLIC :: t_rad_memory, max_no_of_vars

  INTEGER, PARAMETER :: max_no_of_vars = 40   ! default JSBACH4

  !-----------------------------------------------------------------------------------------------------
  !> Type definition for rad memory
  !! 
  !-----------------------------------------------------------------------------------------------------
  TYPE, EXTENDS(t_jsb_memory) :: t_rad_memory

    TYPE(t_jsb_var_real1d)    :: &        ! REAL(wp) 1D ; implemented for SLM use only, see mo_jsb_var_class 
                                alb_nir_dif, &
                                alb_nir_dir, &
                                alb_vis_dif, &
                                alb_vis_dir


  CONTAINS
    PROCEDURE :: Init => Init_rad_memory
  END TYPE t_rad_memory


  CHARACTER(len=*), PARAMETER :: modname = 'mo_rad_memory_class'




CONTAINS


  !-----------------------------------------------------------------------------------------------------
  !> initialize memory (variables) for the process: radiation
  !!
  !! this is a radically simplified version of the init routine when compared to JSBACH4 (20160819) \n
  !! please see Init_assimi_memory where the missign code is just commeted out (i.e. not deleted)
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE Init_rad_memory(mem, prefix, suffix, lct_ids, model_id, ics, ice, iblk)

    CLASS(t_rad_memory), INTENT(inout), TARGET :: mem
    CHARACTER(len=*),    INTENT(in)    :: prefix
    CHARACTER(len=*),    INTENT(in)    :: suffix
    INTEGER,             INTENT(in)    :: lct_ids !(:)
    INTEGER,             INTENT(in)    :: model_id
    INTEGER,             INTENT(in)    :: ics   ! model%options -> for SLM compatibility
    INTEGER,             INTENT(in)    :: ice   ! model%options -> for SLM compatibility
    INTEGER,             INTENT(in)    :: iblk  ! model%options -> for SLM compatibility


    IF (.NOT. ASSOCIATED(mem%alb_nir_dif%ptr)) THEN

      CALL slm_add_var(mem, 'descr_alb_nir_dif', mem%alb_nir_dif, &   
            ics, ice, iblk, &   ! model%options -> for SLM compatibility
            initval_r=0.0_wp)   ! initval taken from JSBACH3

      CALL slm_add_var(mem, 'descr_alb_nir_dir', mem%alb_nir_dir, &   
            ics, ice, iblk, &   ! model%options -> for SLM compatibility
            initval_r=0.0_wp)   ! initval taken from JSBACH3

      CALL slm_add_var(mem, 'descr_alb_vis_dif', mem%alb_vis_dif, &   
            ics, ice, iblk, &   ! model%options -> for SLM compatibility
            initval_r=0.0_wp)   ! initval taken from JSBACH3

      CALL slm_add_var(mem, 'descr_alb_vis_dir', mem%alb_vis_dir, &   
            ics, ice, iblk, &   ! model%options -> for SLM compatibility
            initval_r=0.0_wp)   ! initval taken from JSBACH3


    END IF

  END SUBROUTINE Init_rad_memory

END MODULE mo_rad_memory_class
