!>
!! @brief Contains the routines for the scope model to calc radiative transfer
!!
!! @author
!!  Jan Engel, MPI-BGC Jena - first bare version of this MODULE
!!  Javier Pacheco-Labrador, MPI-BGC Jena - adds actual science
!!
!! @par Revision History
!!  First version by Javier Pacheco-Labrador (dezember 2018)
!!
!!
MODULE mo_rad_scope

#include "dsl_jsbach.inc"

  USE mo_kind,                  ONLY: wp

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: calc_radiative_transfer

  CHARACTER(len=*), PARAMETER :: modname = 'mo_rad_scope'


CONTAINS


  !-----------------------------------------------------------------------------------------------------
  !> calc_radiative_transfer
  !! 
  !!  ... describe what this routine does (major in, calc, major out)
  !! 
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE calc_radiative_transfer( options                   , &   ! in
                                      ncanopy                   , &
                                      nspec                     , &
                                      t_air                     , &
                                      canopy_reflectance_sp       )   ! inout

    USE mo_jsb_task_class,            ONLY: t_jsb_task_options

    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    TYPE(t_jsb_task_options)               , INTENT(in)     :: options                !< model options
    INTEGER                                , INTENT(in)     :: ncanopy                !< number of canopy layers
    INTEGER                                , INTENT(in)     :: nspec                  !< number of spectral bands
    REAL(wp), DIMENSION(options%nc)        , INTENT(in)     :: t_air                  !< air temperature
    REAL(wp), DIMENSION(options%nc,nspec)  , INTENT(inout)  :: canopy_reflectance_sp  !< top of canopy reflectance factor
    ! ---------------------------
    ! 0.2 Local
    ! ...
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_radiative_transfer'


    ! ------------------------------------------------------------------------------------------------------------
    ! Go science
    ! ------------------------------------------------------------------------------------------------------------


    !> 1.0 ...
    !! 
    print*," "
    print*,"Hello World!"
    print*,"The air temp is: ", t_air
    print*,"Number of spectro bands is: ", nspec
    print*,"Init value of the canopy reflectance of the first 10 spectral bands is: ", canopy_reflectance_sp(1,1:10)
    

  END SUBROUTINE calc_radiative_transfer

END MODULE mo_rad_scope
