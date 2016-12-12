!>
!! @brief Contains the interface to the soil processes (the Commission model)
!!
!! @author
!!  Jan Engel, MPI-BGC Jena
!!
!! @par Revision History
!! First version by Jan Engel (2016-07-19)
!!
MODULE mo_soil_interface

#include "dsl_jsbach.inc"

  USE mo_kind,                ONLY: wp
  USE mo_exception,           ONLY: message, finish

!  USE 

  IMPLICIT NONE

!  PRIVATE
!  PUBLIC :: Register_soil_tasks


CONTAINS

  subroutine soil_1st_test_routine(tile)

    USE mo_jsb_class,               ONLY: Get_model
    USE mo_base_memory,             ONLY: t_slm_tile, t_model
    USE mo_jsb_lctlib_class,        ONLY: t_lctlib_element
    USE mo_jsb_process_class,       ONLY: ATM2LAND_, LAND2ATM_, SEB_, TURB_, SOIL_, HYDRO_, VEG_, HD_, RAD_, ASSIMI_, PHENO_

    dsl4jsb_Use_memory(SOIL_)        !  USE mo_soil_memory_class,      ONLY: t_soil_memory

  implicit none

    ! 0.1 InOut
    TYPE(t_slm_tile),INTENT(inout)         :: tile         !< one tile
    ! 0.2 Local
    TYPE(t_model),          POINTER        :: model
    TYPE(t_lctlib_element), POINTER        :: lctlib     !< land-cover-type library - parameter across pft's
    integer         :: test_count = 1
    ! 0.3 Memory allocation
    ! Declare process configuration and memory Pointers
    dsl4jsb_Def_memory(SOIL_)                       ! TYPE(t_soil_memory),      POINTER :: soil__mem
    TYPE(t_soil_memory),     POINTER :: soil        ! local pointer to veg__mem for SLM usability
    dsl4jsb_Get_memory(SOIL_)                       ! soil__mem => tile%mem(SOIL_)%p

    soil    => soil__mem 
    model   => Get_model(tile%owner_model_id)
    lctlib  => model%lctlib(tile%lcts%lib_id)


    soil%pool = soil%pool + 0.001_wp
    soil%pool = soil%pool + soil%pool / 2.0_wp

    soil%soil_depth = soil%soil_depth + 0.001_wp

    write (999,*) " Hi, I am SOIL", soil%soil_depth%ptr(1,1), test_count
    test_count = test_count + 1

  end subroutine soil_1st_test_routine


END MODULE mo_soil_interface
