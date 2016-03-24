!!
!! @brief memory structure for atm->land and land->atm processes (not used yet) 
!!
!! @author
!!  Sönke Zaehle, MPI-BGC Jena 
!!
!! @par Revision History
!! First version by Sönke Zaehle (2015-12-04)
!!
MODULE mo_atmland_constants

  USE mo_kind, ONLY: wp
  USE mo_base_constants, ONLY: hourday

  ! orbital parameters (not calculated - do I need to include this dynamically later?)
  REAL(wp),SAVE :: eccentricity                ! eccentricity of the Earth's orbit
  REAL(wp),SAVE :: perihel                     ! perihel of the Earth's orbit
  REAL(wp),SAVE :: solstice                    ! day of the year at which solstice occurs 
  REAL(wp),SAVE :: obliquity                   ! obliquity of the Earth's orbit
  REAL(wp),SAVE :: declination                 ! solar declination (degrees, angle sun at equator)
  LOGICAL,SAVE  :: angle_at_GMT                ! logical to simulate solar angle according to local time or GMT 
  REAL(wp),SAVE :: reference_longitude_GMT     ! reference longitude for GMT correction 
  REAL(wp),SAVE :: solar_constant              ! solar shortwave radiation (W m-2) 
  REAL(wp),SAVE :: k1_fdiff                    ! constants for calculation of diffuse radiation
  REAL(wp),SAVE :: k2_fdiff                    ! calculation from transmission as in 
  REAL(wp),SAVE :: k3_fdiff                    ! calc_diffuse_radiation_fraction
  REAL(wp),SAVE :: k4_fdiff                    ! 
  REAL(wp),SAVE :: min_sw_fdiffuse             ! minimum diffuse fraction in shortwave downward flux
  REAL(wp),SAVE :: max_transmission_fdiff      ! atmospheric transmission rate above which the 
                                               ! diffuse fraction becomes minimal
  REAL(wp),SAVE :: min_daylength
  REAL(wp),SAVE :: dusk_time 
  REAL(wp),SAVE :: tmp_peak_time 

  ! conversion factors
  REAL(wp),SAVE :: sw2par                      ! W / µmol / s
  REAL(wp),SAVE :: Tzero                       ! melting point in K
  REAL(wp),SAVE :: standard_press_srf          ! standard pressure in Pa
  REAL(wp),SAVE :: eps_vpd                     ! eps_vpd used for unit conversion Pa -> g g-1

  REAL(wp),SAVE :: def_t_air                   ! temperature in deg K
  REAL(wp),SAVE :: def_vpd                     ! vpd in Pa
  REAL(wp),SAVE :: def_rh                      ! relative humidity in fractions
  REAL(wp),SAVE :: def_swdown                  ! PAR in W m-2
  REAL(wp),SAVE :: def_fdiffuse                ! diffuse fraction
  REAL(wp),SAVE :: def_cos_angle               ! cosine of zenith angle
  REAL(wp),SAVE :: def_co2                     ! co2 in ppm
  REAL(wp),SAVE :: def_wind                    ! wind in m/s
  REAL(wp),SAVE :: def_latitude                ! latitude in degrees
  REAL(wp),SAVE :: def_longitude               ! longitude in degrees

  PUBLIC :: init_atmland_constants
 
CONTAINS

  SUBROUTINE init_atmland_constants
  ! 
  ! subroutine to initialise atmland constants
  !
  IMPLICIT NONE
  ! 
  ! initialise Earth's orbit
  !  
  obliquity = 23.446_wp
  perihel = 102.04_wp
  eccentricity = 0.016724_wp
  solstice = 174.0_wp   ! 23rd June
  !
  ! initialise solar angle calculation
  !
  angle_at_GMT = .FALSE.
  reference_longitude_GMT = 0.0_wp
  ! 
  ! initialise solar constant
  ! 
  solar_constant = 1370.0_wp
  ! 
  ! initalise diffuse fraction calculation
  k1_fdiff = 1.0045_wp
  k2_fdiff = 0.0435_wp
  k3_fdiff = -3.5227_wp
  k4_fdiff = 2.6313_wp
  min_sw_fdiffuse = 0.1667_wp 
  max_transmission_fdiff = 0.75_wp
  ! 
  ! initialise diurnal temperature calculation
  min_daylength = 4._wp / hourday
  tmp_peak_time = 14.0_wp / hourday
  dusk_time = 1.0_wp / hourday
  ! 
  ! conversion factors
  sw2par              = 2.3_wp  
  Tzero               = 273.15_wp  
  standard_press_srf  = 101300.0_wp  
  eps_vpd             = 0.622_wp  
  ! 
  ! default values
  def_t_air = 298.15_wp
  def_vpd = 1000._wp   
  def_rh = 0.8_wp      
  def_swdown = 652.173913_wp
  def_fdiffuse = 0.2_wp     
  def_cos_angle = 1.0_wp    
  def_co2 = 380.0_wp        
  def_wind = 1.0_wp         
  def_latitude = 50.0_wp    
  def_longitude = 0.0_wp    


  END SUBROUTINE init_atmland_constants

END MODULE mo_atmland_constants
