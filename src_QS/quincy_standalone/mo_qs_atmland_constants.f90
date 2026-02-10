!> Atmosphere-land constants, orbital parameters, conversion factors, plus init routine 
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
MODULE mo_qs_atmland_constants

  USE mo_kind,            ONLY: wp

  IMPLICIT NONE

  ! orbital parameters (not calculated - do I need to include this dynamically later?)
  REAL(wp),SAVE :: obliquity                   !< obliquity of the Earth's orbit
  REAL(wp),SAVE :: perihel                     !< perihel of the Earth's orbit
  REAL(wp),SAVE :: eccentricity                !< eccentricity of the Earth's orbit
  REAL(wp),SAVE :: solstice                    !< day of the year at which solstice occurs 
  REAL(wp),SAVE :: declination                 !< solar declination (degrees, angle sun at equator)
                                               !! not declared here, but calculated in mo_qs_atmland_util:update_orbital_parameters
                                               !! which is called in mo_qs_atmland_forcing:update_forcing
  
  ! solar angle calculation
  LOGICAL,SAVE  :: angle_at_GMT                !< logical to simulate solar angle according to local time or GMT 
  REAL(wp),SAVE :: reference_longitude_GMT     !< reference longitude for GMT correction 
  REAL(wp),SAVE :: min_cos_zenith_angle        !< minimum zenith angle for which PAR is calculated, see cos_zenith_angle
  
  ! solar constant
  REAL(wp),SAVE :: solar_constant              !< solar shortwave radiation (W m-2) 
  
  ! diffuse fraction calculation
  REAL(wp),SAVE :: k1_fdiff                    !< constants for calculation of diffuse radiation
  REAL(wp),SAVE :: k2_fdiff                    !< calculation from transmission as in 
  REAL(wp),SAVE :: k3_fdiff                    !< calc_diffuse_radiation_fraction
  REAL(wp),SAVE :: k4_fdiff                    !<  
  REAL(wp),SAVE :: min_fract_par_diffuse
    !< minimum diffuse fraction in shortwave downward flux. Note: sw diffuse fraction assumed to be identical to PAR diffuse fract
  REAL(wp),SAVE :: max_transmission_fdiff      !< atmospheric transmission rate above which the diffuse fraction becomes minimal
  
  ! diurnal temperature calculation
  REAL(wp),SAVE :: min_daylength               !< ?
  REAL(wp),SAVE :: dusk_time                   !< ?
  REAL(wp),SAVE :: tmp_peak_time               !< ?

  PUBLIC

  CHARACTER(len=*), PARAMETER, PRIVATE :: modname = 'mo_qs_atmland_constants'
 
CONTAINS

  !-----------------------------------------------------------------------------------------------------
  !> initialise atmland constants
  !! 
  !! initialise Earth's orbit; solar angle calculation, solar constant, diffuse fraction calculation
  !! diurnal temperature calculation, conversion factors, default values
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE init_qs_atmland_constants

    USE mo_jsb_math_constants,    ONLY: hourday

    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    !
    ! ---------------------------
    ! 0.2 Local
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':init_qs_atmland_constants'


    !> initialise Earth's orbit
    !!  
    obliquity                 = 23.446_wp
    perihel                   = 102.04_wp
    eccentricity              = 0.016724_wp
    solstice                  = 174.0_wp   ! 23rd of June

    !> initialise solar angle calculation
    !! 
    angle_at_GMT              = .FALSE.
    reference_longitude_GMT   = 0.0_wp
    min_cos_zenith_angle      = 0.02_wp   ! in JSBACH they use: mo_bethy_constants.f90:  REAL(dp), PARAMETER :: ZenithMinPar = 1.E-3_dp

    !> initialise solar constant
    !! 
    solar_constant            = 1370.0_wp

    !> initalise diffuse fraction calculation
    !! 
    k1_fdiff                  = 1.0045_wp
    k2_fdiff                  = 0.0435_wp
    k3_fdiff                  = -3.5227_wp
    k4_fdiff                  = 2.6313_wp
    min_fract_par_diffuse     = 0.1667_wp 
    max_transmission_fdiff    = 0.75_wp

    !> initialise diurnal temperature calculation
    !! 
    min_daylength             = 4._wp / hourday
    tmp_peak_time             = 14.0_wp / hourday
    dusk_time                 = 1.0_wp / hourday

  END SUBROUTINE init_qs_atmland_constants

END MODULE mo_qs_atmland_constants
