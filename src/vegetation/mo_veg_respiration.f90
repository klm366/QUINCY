!>
!! @brief Contains the routines for the plant respiration 
!!
!! @author
!!  Sönke Zaehle, MPI-BGC Jena 
!!
!! @par Revision History
!! First version by Sönke Zaehle (2016-02-17)
!!
MODULE mo_veg_respiration

  USE mo_kind

  PRIVATE

  PUBLIC :: calc_maintenance_respiration,temperature_response_respiration

CONTAINS

  ELEMENTAL FUNCTION temperature_response_respiration(temp) RESULT(temp_response)
    ! 
    ! function to calculate the instantaneous temperature response of maintenance respiration
    ! following the temperature function of Lloyd and Taylor, 1996, Functional Ecology
    !  
    USE mo_veg_parameters, ONLY: t_maint_k1, t_maint_k2, t_maint_k3

    IMPLICIT NONE

    ! 0.1 Input
    REAL(wp), INTENT(in) :: temp          ! instantaneous temperature (K)
    ! 0.3 Output
    REAL(wp) :: temp_response             ! temperature response of respiration (1 at 10°C)

    temp_response = exp(t_maint_k1 * (1._wp/t_maint_k2 - 1._wp / (temp - t_maint_k3)))

  END FUNCTION temperature_response_respiration


  SUBROUTINE calc_maintenance_respiration(t_stem,t_root,sap_wood_nitrogen,fine_root_nitrogen,coarse_root_nitrogen,&
                  labile_carbon,labile_carbon13,labile_carbon14,maint_respiration_leaf,& 
                  maint_respiration,maint_respiration_c13,maint_respiration_c14)
    ! 
    ! Routine to calculate plant maintenance respiration based on the tissue N contents
    ! 
    USE mo_veg_parameters, ONLY: fmaint_rate_base

    IMPLICIT NONE

    ! 0.1 Input
    REAL(wp), INTENT(in) :: t_stem                   ! temperature of the stem (K) 
    REAL(wp), INTENT(in) :: t_root                   ! temperature of the roots (K)
    REAL(wp), INTENT(in) :: sap_wood_nitrogen        ! nitrogen content of sapwood (mol N / m2)
    REAL(wp), INTENT(in) :: coarse_root_nitrogen     ! nitrogen content of sapwood (mol N / m2)
    REAL(wp), INTENT(in) :: fine_root_nitrogen       ! nitrogen content of sapwood (mol N / m2)
    REAL(wp), INTENT(in) :: labile_carbon            ! labile carbon (mol C / m2)
    REAL(wp), INTENT(in) :: labile_carbon13          ! labile carbon13 (mol 13C / m2)
    REAL(wp), INTENT(in) :: labile_carbon14          ! labile carbon14 (mod 14C / m2)
    REAL(wp), INTENT(in) :: maint_respiration_leaf   ! maintenance respiration of leaves 
                                                     ! (calculated in photosynthesis; µmol C / m2 / s)
    ! 0.3 Output
    REAL(wp), INTENT(out) :: maint_respiration       ! maintenance respiration (µmol C / m2 / s)
    REAL(wp), INTENT(out) :: maint_respiration_c13   ! maintenance respiration (µmol 13C / m2 / s)
    REAL(wp), INTENT(out) :: maint_respiration_c14   ! maintenance respiration (µmol 14C / m2 / s)
    ! 
    ! 1.0 maintenance respiration for the stem
    ! 
    maint_respiration = temperature_response_respiration(t_stem) * & 
                        fmaint_rate_base * sap_wood_nitrogen
    ! 
    ! 2.0 maintenance respiration for the roots 
    ! 
    maint_respiration = maint_respiration + &
        temperature_response_respiration(t_root) * fmaint_rate_base * &  
        ( coarse_root_nitrogen + fine_root_nitrogen ) 
    ! 
    ! 3.0 maintenance respiration for the leaves (already calculated during photosynthesis, parsed in here) 
    ! 
    maint_respiration = maint_respiration + maint_respiration_leaf 
    !
    ! 4.0 Isotopic signature of maintenance respiration is that of the labile pool 
    !
    IF(labile_carbon > eps8) THEN 
      maint_respiration_c13 = maint_respiration * labile_carbon13/labile_carbon
      maint_respiration_c14 = maint_respiration * labile_carbon14/labile_carbon
    ELSE 
      maint_respiration_c13 = 0.0_wp
      maint_respiration_c14 = 0.0_wp
    ENDIF

  END SUBROUTINE calc_maintenance_respiration

END MODULE mo_veg_respiration
