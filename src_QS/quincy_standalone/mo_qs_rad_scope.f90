!> Contains the routines for the scope model to calc radiative transfer
!>
!> QUINCY standalone
!>
!> ---------------------------------------
!> Copyright (C) 2013-2026, MPI-M, MPI-BGC
!>
!> Contact: https://doi.org/10.17871/quincy-model-2019
!> Authors: AUTHORS.md
!> See LICENSES/ for license information
!> SPDX-License-Identifier: BSD-3-Clause
!> ---------------------------------------
!>
MODULE mo_qs_rad_scope

  USE mo_kind,                  ONLY: wp

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: calc_radiative_transfer

  CHARACTER(len=*), PARAMETER :: modname = 'mo_qs_rad_scope'


CONTAINS


  !-----------------------------------------------------------------------------------------------------
  !> calc_radiative_transfer
  !! 
  !!  ... describe what this routine does (major in, calc, major out)
  !! 
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE calc_radiative_transfer( nc                        , &   ! in
                                      ncanopy                   , &
                                      nspec                     , &
                                      t_air                     , &
                                      canopy_reflectance_sp       )   ! inout

    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    INTEGER                                , INTENT(in)     :: nc             !< number of points in the chunk
    INTEGER                                , INTENT(in)     :: ncanopy        !< number of canopy layers
    INTEGER                                , INTENT(in)     :: nspec          !< number of spectral bands
    REAL(wp), DIMENSION(nc)        , INTENT(in)     :: t_air                  !< air temperature
    REAL(wp), DIMENSION(nc,nspec)  , INTENT(inout)  :: canopy_reflectance_sp  !< top of canopy reflectance factor
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

END MODULE mo_qs_rad_scope
