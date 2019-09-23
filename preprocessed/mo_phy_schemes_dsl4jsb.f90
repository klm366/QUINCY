MODULE mo_phy_schemes
 USE mo_kind , ONLY: wp
 PUBLIC :: qsat
 PUBLIC :: q_effective
 PUBLIC :: surface_dry_static_energy
 PUBLIC :: heat_transfer_coef
 PUBLIC :: thermal_radiation
 PUBLIC :: lwnet_from_lwdown
 CHARACTER ( LEN= * ) , PRIVATE , PARAMETER :: modname ="mo_phy_schemes"

CONTAINS
 FUNCTION qsat ( temperature , pressure )
  USE mo_echam_convect_tables , ONLY: tlucua , jptlucu1 , jptlucu2
  USE mo_jsb_physical_constants , ONLY: rvd1
  REAL ( KIND= wp ) :: qsat
  REAL ( KIND= wp ) , INTENT(IN) :: temperature
  REAL ( KIND= wp ) , INTENT(IN) :: pressure
  INTEGER :: it
  REAL ( KIND= wp ) :: tluc

!$acc routine seq
  it = nint ( temperature * 1000._wp )
  tluc = tlucua ( it )
  IF ( it >= jptlucu1 .AND. it <= jptlucu2 ) THEN
   qsat = tluc / ( pressure - rvd1 * tluc )
  ELSE
   qsat = huge ( 1._wp )
  END IF
 END FUNCTION qsat

 FUNCTION q_effective ( qsat , qair , fsat , fair )
  REAL ( KIND= wp ) :: q_effective
  REAL ( KIND= wp ) , INTENT(IN) :: qsat
  REAL ( KIND= wp ) , INTENT(IN) :: qair
  REAL ( KIND= wp ) , INTENT(IN) :: fsat
  REAL ( KIND= wp ) , INTENT(IN) :: fair

!$acc routine seq
  q_effective = fsat * qsat + ( 1._wp - fair ) * qair
 END FUNCTION q_effective

 FUNCTION surface_dry_static_energy ( t_srf , qsat_srf )
  USE mo_jsb_physical_constants , ONLY: cpd , cpvd1
  REAL ( KIND= wp ) :: surface_dry_static_energy
  REAL ( KIND= wp ) , INTENT(IN) :: t_srf
  REAL ( KIND= wp ) , INTENT(IN) :: qsat_srf

!$acc routine seq
  surface_dry_static_energy = t_srf * cpd * ( 1._wp + cpvd1 * qsat_srf )
 END FUNCTION surface_dry_static_energy

 ELEMENTAL FUNCTION thermal_radiation ( t_srf )
  USE mo_jsb_physical_constants , ONLY: stbo , zemiss_def
  REAL ( KIND= wp ) :: thermal_radiation
  REAL ( KIND= wp ) , INTENT(IN) :: t_srf

  thermal_radiation = stbo * zemiss_def * t_srf ** ( 4 )
 END FUNCTION thermal_radiation

 FUNCTION lwnet_from_lwdown ( lwdown , t_srf )
  USE mo_physical_constants , ONLY: stbo , zemiss_def
  REAL ( KIND= wp ) :: lwnet_from_lwdown
  REAL ( KIND= wp ) , INTENT(IN) :: lwdown
  REAL ( KIND= wp ) , INTENT(IN) :: t_srf

!$acc routine seq
  lwnet_from_lwdown = zemiss_def * ( lwdown - stbo * t_srf ** ( 4 ) )
 END FUNCTION lwnet_from_lwdown

 FUNCTION heat_transfer_coef ( drag , steplen )
  USE mo_jsb_physical_constants , ONLY: grav , cvdifts
  REAL ( KIND= wp ) :: heat_transfer_coef
  REAL ( KIND= wp ) , INTENT(IN) :: drag
  REAL ( KIND= wp ) , INTENT(IN) :: steplen

!$acc routine seq
  heat_transfer_coef = drag / ( cvdifts * grav * steplen )
 END FUNCTION heat_transfer_coef

END MODULE mo_phy_schemes

