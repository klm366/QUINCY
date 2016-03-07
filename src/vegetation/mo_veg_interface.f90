!!
!! @brief Contains the interface for plant growth and turnover processes 
!!
!! @author
!!  Sönke Zaehle, MPI-BGC Jena 
!!
!! @par Revision History
!! First version by Sönke Zaehle (2015-12-04)
!!
MODULE mo_veg_interface

  USE mo_kind
  USE mo_base_memory

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: update_vegetation 

CONTAINS

  SUBROUTINE update_vegetation(tile,t_air,veg_param)
    !
    ! calls the vegetation interface
    ! TBD: Merge with interface? Somewhat a strategic decision for all header routines
    ! where will the loop over tiles and gridcells happen. For the case of vegeation it may
    ! be beneficial to keep the loop outside the actual growth routine
    ! 
    USE mo_land_constants
    USE mo_veg_parameters, ONLY: VegetationParameterClass

    IMPLICIT NONE

    ! 0.1 Input
    TYPE(VegetationParameterClass),INTENT(in) :: veg_param 
    REAL(wp),INTENT(in) :: t_air 
    ! 0.2 In- and Output
    TYPE(TileClass), DIMENSION(ntile), INTENT(inout) :: tile
    ! 0.4 Local
    INTEGER :: t

    DO t = 1,ntile
 
       CALL vegetation_interface(tile(t)%plant,t_air,veg_param)

    END DO
 
  END SUBROUTINE update_vegetation

  SUBROUTINE vegetation_interface(plant,t_air,veg_param)
    ! 
    ! update vegetation pools given pre-calculated rates of 
    ! net C, N and P uptake 
    ! !! OBS: incomplete!!
    ! 
    USE mo_veg_memory, ONLY: VegPoolAdd,VegPoolSub
    USE mo_veg_parameters, ONLY: VegetationParameterClass
    USE mo_veg_turnover
    USE mo_veg_growth
    USE mo_veg_respiration, ONLY: calc_maintenance_respiration
    IMPLICIT NONE
    ! 0.1 Input
    TYPE(VegetationParameterClass),INTENT(in) :: veg_param 
    REAL(wp),INTENT(in) :: t_air 
    ! 0.2 In- and Output
    TYPE(PlantClass),INTENT(inout) :: plant
    ! 
    ! 1.0 Calculate growth
    ! 
    CALL  calc_growth(plant,veg_param,plant%flux%growth, & 
            plant%flux%growth_respiration,plant%flux%growth_respiration_c13,plant%flux%growth_respiration_c14)
    ! 
    ! 1.5 Calculate maintenance respiration
    ! 
    CALL  calc_maintenance_respiration(t_air,t_air,plant%pool%sap_wood%nitrogen,plant%pool%fine_root%nitrogen, &
            plant%pool%coarse_root%nitrogen,plant%pool%labile%carbon,plant%pool%labile%carbon13, & 
            plant%pool%labile%carbon14,plant%flux%maint_respiration_leaf,&
            plant%flux%maint_respiration,plant%flux%maint_respiration_c13,plant%flux%maint_respiration_c14)
    ! 
    ! 2.0 Calculate turnover from tissue aging and senescence
    !
    CALL calc_turnover(plant,veg_param)
    !
    ! 3.0 Add carbon from photosynthesis to the labile pool
    !
    plant%flux%npp = plant%flux%gpp - plant%flux%growth_respiration - plant%flux%maint_respiration 
    plant%flux%npp_c13 = plant%flux%gpp_c13 - plant%flux%growth_respiration_c13 - plant%flux%maint_respiration_c13 
    plant%flux%npp_c14 = plant%flux%gpp_c14 - plant%flux%growth_respiration_c14 - plant%flux%maint_respiration_c14 
    plant%flux%growth%labile%carbon = plant%flux%growth%labile%carbon + plant%flux%npp * dt / 1.e6
    plant%flux%growth%labile%carbon13 = plant%flux%growth%labile%carbon13 + & 
       plant%flux%npp_c13 * dt / 1.e6
    plant%flux%growth%labile%carbon14 = plant%flux%growth%labile%carbon14 + & 
       plant%flux%npp_c14 * dt / 1.e6
    !
    ! 4.0 Add nutrients from uptake calculations
    !
    plant%flux%uptake_n = MAX( ( plant%pool%labile%carbon + plant%flux%growth%labile%carbon) /25.0_wp - & 
                plant%pool%labile%nitrogen, 0.0_wp ) 
    plant%flux%uptake_n15 = plant%flux%uptake_n * 0.0002_wp 
    plant%flux%uptake_p = MAX( (plant%pool%labile%nitrogen+plant%flux%uptake_n)/14.0_wp - &
                plant%pool%labile%phosphorus, 0.0_wp ) 
    plant%flux%growth%labile%nitrogen = plant%flux%growth%labile%nitrogen + plant%flux%uptake_n
    plant%flux%growth%labile%nitrogen15 = plant%flux%growth%labile%nitrogen15 + plant%flux%uptake_n15
    plant%flux%growth%labile%phosphorus = plant%flux%growth%labile%phosphorus + plant%flux%uptake_p
    !
    ! 5.0 Update plant biogeochemical pools with growth and litterfall
    ! 
    plant%pool = VegPoolSub(plant%pool,plant%flux%litterfall)
    plant%pool = VegPoolAdd(plant%pool,plant%flux%growth)
    ! 
    ! 6.0 Update plant diagnostics
    !
    plant%lai = plant%pool%leaf%carbon * veg_param%sla(plant%pftid)

  END SUBROUTINE vegetation_interface

END MODULE mo_veg_interface
