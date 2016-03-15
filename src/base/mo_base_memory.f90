!!
!! @brief land model memory structure
!!
!! @author
!!  Sönke Zaehle, MPI-BGC Jena 
!!
!! @par Revision History
!! First version by Sönke Zaehle (2015-12-04)
!!
MODULE mo_base_memory

  USE mo_kind
  USE mo_base_constants

  TYPE Element
     REAL(wp) :: carbon = 0.0_wp     ! total carbon concentration of pool in mol / m2
     REAL(wp) :: nitrogen = 0.0_wp   ! totol nitrogen concentration of pool in mol / m2
     REAL(wp) :: phosphorus = 0.0_wp ! phosphorus concentration of pool in mol / m2
     REAL(wp) :: carbon13 = 0.0_wp   ! 13carbon concentration of pool in mol / m2
     REAL(wp) :: carbon14 = 0.0_wp   ! 14carbon concentration of pool in mol / m2
     REAL(wp) :: nitrogen15 = 0.0_wp ! 15nitrogen concentration of pool in mol / m2
  END TYPE Element

  TYPE CanopyPartitioningClass       ! partitioning of canopy N into photosynthetic fractions
     REAL(wp) :: nl_tot(ncanop) = 0.0_wp     ! N in canopy layer (mmol / m2 (LAI))
     REAL(wp) :: fnl_chl(ncanop) = 0.0_wp    ! fraction of layer N in chlorophyll (pigments only!)
     REAL(wp) :: fnl_et(ncanop) = 0.0_wp     ! fraction of layer N in electron transport 
     REAL(wp) :: fnl_rub(ncanop) = 0.0_wp    ! fraction of layer N in Rubisco 
  END TYPE CanopyPartitioningClass

  TYPE VegPool                       ! vegetation pools 
     TYPE(Element) :: leaf           ! foliage
     TYPE(Element) :: fine_root      ! fine roots (< 2mm)
     TYPE(Element) :: coarse_root    ! coarse roots (> 2mm)
     TYPE(Element) :: sap_wood       ! living wood for water and phloem transport
     TYPE(Element) :: heart_wood     ! dead standing wood
     TYPE(Element) :: labile         ! meristem carbon
     TYPE(Element) :: reserve        ! storage
     TYPE(Element) :: fruit          ! live propagation pool
     TYPE(Element) :: seed_bed       ! fruit mass on ground
     TYPE(CanopyPartitioningClass) :: foliarN ! probably undesirable feature to have this here
  END TYPE VegPool

  TYPE VegFlux
     TYPE(VegPool) :: growth
     TYPE(VegPool) :: litterfall
     REAL(wp) :: gpp                      ! gross photosynthesis (µmol CO2 m-2 s-1)
     REAL(wp) :: gpp_c13                  ! gross photosynthesis (µmol 13CO2 m-2 s-1)
     REAL(wp) :: gpp_c14                  ! gross photosynthesis (µmol 14CO2 m-2 s-1)
     REAL(wp) :: maint_respiration_leaf   ! foliar maintenance respiration (µmol CO2 m-2 s-1)
     REAL(wp) :: maint_respiration        ! maintenance respiration (µmol CO2 m-2 s-1)
     REAL(wp) :: maint_respiration_c13    ! maintenance respiration (µmol 13CO2 m-2 s-1)
     REAL(wp) :: maint_respiration_c14    ! maintenance respiration (µmol 14CO2 m-2 s-1)
     REAL(wp) :: growth_respiration       ! growth respiration (µmol CO2 m-2 s-1)
     REAL(wp) :: growth_respiration_c13   ! growth respiration (µmol 13CO2 m-2 s-1)
     REAL(wp) :: growth_respiration_c14   ! growth respiration (µmol 14CO2 m-2 s-1)
     REAL(wp) :: npp                      ! net primary production (µmol CO2 m-2 s-1)
     REAL(wp) :: npp_c13                  ! net primary production (µmol 13CO2 m-2 s-1)
     REAL(wp) :: npp_c14                  ! net primary production (µmol 14CO2 m-2 s-1)
     REAL(wp) :: uptake_n
     REAL(wp) :: uptake_n15
     REAL(wp) :: uptake_p
     REAL(wp) :: gc
     REAL(wp) :: ci
  END TYPE VegFlux

  ! 
  ! Plant
  ! 
  TYPE AccumulationPlantClass
    ! definition of averages of variable across a certain period (model time steps)
    ! Light in the canopy
    REAL(wp) :: par_sunlit = 1._wp
    REAL(wp) :: par_shaded = 1._wp
    REAL(wp) :: fleaf_sunlit = 1._wp
    ! Soil moisture (currently in mo_veg_photosynthesis)
    REAL(wp) :: beta_air = 1._wp
    REAL(wp) :: beta_soil_ps = 1._wp
    REAL(wp) :: beta_soil_gs = 1._wp
  END TYPE AccumulationPlantClass

  TYPE PlantClass
     TYPE(VegPool) :: pool
     TYPE(VegFlux) :: flux
     TYPE(AccumulationPlantClass) :: accumulated
     INTEGER :: pftid = 1
     REAL(wp) :: lai
     REAL(wp) :: par_sunlit(ncanop)
     REAL(wp) :: par_shaded(ncanop)
     REAL(wp) :: fleaf_sunlit(ncanop)
  END TYPE PlantClass

  ! 
  ! Soil
  ! 
! activate when soil model needs it ! no init subroutine implemented yet
!  TYPE AccumulationSoilClass
!    ! definition of averages of variable across a certain period (model time steps)
!    REAL(wp) :: x_something = 1.
!  END TYPE AccumulationSoilClass

  TYPE SoilClass
    !TYPE(AccumulationSoilClass) :: accumulated    ! activate when soil model needs it
    REAL(wp) :: theta
  END TYPE SoilClass

  ! 
  ! Tile (integrates plant & soil)
  ! 
  TYPE AccumulationTileClass
    ! definition of averages of variable across a certain period (model time steps)
    ! Environmental forcing
    REAL(wp) :: t_air = 1._wp
    REAL(wp) :: q_air = 1._wp
    REAL(wp) :: press_surf = 1._wp
  END TYPE AccumulationTileClass

  TYPE TileClass
     TYPE(PlantClass) :: plant
     TYPE(SoilClass) :: soil
     TYPE(AccumulationTileClass) :: accumulated
  END TYPE TileClass  
 
! PUBLIC :: ElementAdd,ElementSub,ElementMul,ElementDiv
 
CONTAINS

  FUNCTION ElementAdd(e1,e2) RESULT(e3)
     IMPLICIT NONE
     TYPE(Element), INTENT(in) :: e1,e2
     TYPE(Element) :: e3
     e3%carbon = e1%carbon + e2%carbon 
     e3%nitrogen = e1%nitrogen + e2%nitrogen 
     e3%phosphorus = e1%phosphorus + e2%phosphorus 
     e3%carbon13 = e1%carbon13 + e2%carbon13 
     e3%carbon14 = e1%carbon14 + e2%carbon14
     e3%nitrogen15 = e1%nitrogen15 + e2%nitrogen15
  END FUNCTION ElementAdd  

  FUNCTION ElementSub(e1,e2) RESULT(e3)
     IMPLICIT NONE
     TYPE(Element), INTENT(in) :: e1,e2
     TYPE(Element) :: e3
     e3%carbon = e1%carbon - e2%carbon 
     e3%nitrogen = e1%nitrogen - e2%nitrogen 
     e3%phosphorus = e1%phosphorus - e2%phosphorus 
     e3%carbon13 = e1%carbon13 - e2%carbon13 
     e3%carbon14 = e1%carbon14 - e2%carbon14
     e3%nitrogen15 = e1%nitrogen15 - e2%nitrogen15
  END FUNCTION ElementSub  

  FUNCTION ElementMul(e1,e2) RESULT(e3)
     IMPLICIT NONE
     TYPE(Element), INTENT(in) :: e1,e2
     TYPE(Element) :: e3
     e3%carbon = e1%carbon * e2%carbon 
     e3%nitrogen = e1%nitrogen * e2%nitrogen 
     e3%phosphorus = e1%phosphorus * e2%phosphorus 
     e3%carbon13 = e1%carbon13 * e2%carbon13 
     e3%carbon14 = e1%carbon14 * e2%carbon14
     e3%nitrogen15 = e1%nitrogen15 * e2%nitrogen15
  END FUNCTION ElementMul

  FUNCTION ElementDiv(e1,e2) RESULT(e3)
     IMPLICIT NONE
     TYPE(Element), INTENT(in) :: e1,e2
     TYPE(Element) :: e3
     e3%carbon = e1%carbon / e2%carbon 
     e3%nitrogen = e1%nitrogen / e2%nitrogen 
     e3%phosphorus = e1%phosphorus / e2%phosphorus 
     e3%carbon13 = e1%carbon13 / e2%carbon13 
     e3%carbon14 = e1%carbon14 / e2%carbon14
     e3%nitrogen15 = e1%nitrogen15 / e2%nitrogen15
  END FUNCTION ElementDiv  

END MODULE
