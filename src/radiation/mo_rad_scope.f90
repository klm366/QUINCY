!>
!!  @par Copyright
!!    QUINCY terrestrial biosphere model with coupled carbon, nitrogen, and phosphorus cycles
!!    Copyright (C) 2019 Caldararu S., Engel J., Kern M., Pallandt M., Schnur R., Thum T., Yu L., and Zaehle S.
!!
!!    This file is part of the QUINCY model.
!!    The QUINCY model is free software: you can redistribute it and/or modify
!!    it under the terms of the GNU General Public License as published by
!!    the Free Software Foundation, either version 3 of the License, or
!!    (at your option) any later version.
!!
!!    The QUINCY model is distributed in the hope that it will be useful,
!!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!!    GNU General Public License for more details.
!!
!!    You should have received a copy of the GNU General Public License
!!    along with the QUINCY model.  If not, see <https://www.gnu.org/licenses/>.
!! 
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

END MODULE mo_rad_scope
