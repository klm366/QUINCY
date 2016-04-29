!>
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
     REAL(wp) :: carbon      ! total carbon concentration of pool in mol / m2
     REAL(wp) :: nitrogen    ! totol nitrogen concentration of pool in mol / m2
     REAL(wp) :: phosphorus  ! phosphorus concentration of pool in mol / m2
     REAL(wp) :: carbon13    ! 13carbon concentration of pool in mol / m2
     REAL(wp) :: carbon14    ! 14carbon concentration of pool in mol / m2
     REAL(wp) :: nitrogen15  ! 15nitrogen concentration of pool in mol / m2
  END TYPE Element

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
  END TYPE VegPool

  TYPE VegFlux
     TYPE(VegPool) :: growth         !< change in plant%pool due to growth per timestep growth of tissues (mol m-2 timestep-1)
     TYPE(VegPool) :: litterfall     !< change in plant%pool due to litterfall (incl. mortality) per timestep (mol m-2 timestep-1)
     REAL(wp) :: gpp                      !< gross photosynthesis (µmol CO2 m-2 s-1)
     REAL(wp) :: an                       !< net photosynthesis (µmol CO2 m-2 s-1)
     REAL(wp) :: gpp_c13                  !< gross photosynthesis (µmol 13CO2 m-2 s-1)
     REAL(wp) :: gpp_c14                  !< gross photosynthesis (µmol 14CO2 m-2 s-1)
     REAL(wp) :: maint_respiration_leaf   !< foliar maintenance respiration (µmol CO2 m-2 s-1)
     REAL(wp) :: maint_respiration        !< maintenance respiration (µmol CO2 m-2 s-1)
     REAL(wp) :: maint_respiration_c13    !< maintenance respiration (µmol 13CO2 m-2 s-1)
     REAL(wp) :: maint_respiration_c14    !< maintenance respiration (µmol 14CO2 m-2 s-1)
     REAL(wp) :: growth_respiration       !< growth respiration (µmol CO2 m-2 s-1)
     REAL(wp) :: growth_respiration_c13   !< growth respiration (µmol 13CO2 m-2 s-1)
     REAL(wp) :: growth_respiration_c14   !< growth respiration (µmol 14CO2 m-2 s-1)
     REAL(wp) :: npp                      !< net primary production (µmol CO2 m-2 s-1)
     REAL(wp) :: npp_c13                  !< net primary production (µmol 13CO2 m-2 s-1)
     REAL(wp) :: npp_c14                  !< net primary production (µmol 14CO2 m-2 s-1)
     REAL(wp) :: uptake_n
     REAL(wp) :: uptake_n15
     REAL(wp) :: uptake_p
     REAL(wp) :: leaf_resorption_n        !< resorption of N upon leaf senescence (mol m-2 timestep-1)
     REAL(wp) :: leaf_resorption_n15      !< resorption of 15N upon leaf senescence (mol m-2 timestep-1)
     REAL(wp) :: leaf_resorption_p        !< resorption of P upon leaf senescence (mol m-2 timestep-1)
     REAL(wp) :: transpiration            !< transpiration in kg m-2 s-1
     REAL(wp) :: ci
     TYPE(Element) :: reserve_use         !< net flux between reserve and labile pool (mol m-2 timestep-1)
                                          !!  (+ve: flux to labile; -ve flux from labile) 
  END TYPE VegFlux

  !-----------------------------------------------------------------------------------------------------
  !> vegetation canopy structure (for each layer)
  !! 
  !! number of layers is: ncanop
  !-----------------------------------------------------------------------------------------------------
  TYPE VegCanopy
     REAL(wp) :: nl_tot                !< N in canopy layer (mmol / m2 (LAI))
     REAL(wp) :: fnl_chl               !< fraction of N in chlorophyll (pigments only!)
     REAL(wp) :: fnl_et                !< fraction of N in electron transport 
     REAL(wp) :: fnl_rub               !< fraction of N in Rubisco 
     REAL(wp) :: par_sunlit            !< PAR on sunlit leaves (µmol m-2 s-1)
     REAL(wp) :: par_shaded            !< PAR on shaded leaves (µmol m-2 s-1)
     REAL(wp) :: fleaf_sunlit          !< fraction of layer which is sunlit
  END TYPE VegCanopy

 TYPE AccumulationPlantClass
    ! definition of averages of variable across a certain period (model time steps)
    ! init subroutine in mo_veg_memory
    ! Light in the canopy
    REAL(wp),DIMENSION(ncanop) :: par_sunlit = 1._wp
    REAL(wp),DIMENSION(ncanop) :: par_shaded = 1._wp
    REAL(wp),DIMENSION(ncanop) :: fleaf_sunlit = 1._wp
    ! Soil moisture (currently in mo_veg_photosynthesis)
    REAL(wp) :: beta_air = 0._wp
    REAL(wp) :: beta_soil_ps = 0._wp
    REAL(wp) :: beta_soil_gs = 0._wp
    ! memory of fluxes
    REAL(wp) :: gpp_week = 0._wp                   ! gpp of the last week (µmol m-2 s-1)
    REAL(wp) :: maint_respiration_week = 0._wp     ! maintenance respiration of the last week (µmol m-2 s-1)
    ! nutrient demand for plant uptake
    REAL(wp) :: demand_uptake_n = 0.0_wp           ! scalar of Vmax for root N uptake 
    REAL(wp) :: demand_uptake_p = 0.0_wp           ! scalar of Vmax for root P uptake
    ! phenological variables
    REAL(wp) :: gdd                                ! growing degree days above 5°C
    REAL(wp) :: nd_dormance                        ! number of days since last growing season 
  END TYPE AccumulationPlantClass

  TYPE PlantClass
     TYPE(VegPool) :: pool               !< the actual vegetation pools
     TYPE(VegFlux) :: flux               !< contains growth and litter fall fluxes per timestep for each pool
     TYPE(VegCanopy),DIMENSION(ncanop) :: canopy_layer
     TYPE(AccumulationPlantClass)      :: accumulated
     INTEGER  :: pftid                             !< ID of the plant functional type; namelist parameter
     REAL(wp) :: height                            !< plant height in m
     REAL(wp) :: diameter                          !< plant diameter in m
     REAL(wp) :: lai                               !< actual leaf area index
     REAL(wp) :: target_cn_leaf                    !< target foliar C:N ratio
     REAL(wp) :: target_np_leaf                    !< target foliar N:P ratio 
     REAL(wp) :: mean_leaf_age                     !< average foliage age in days 
     REAL(wp) :: gc                                !< ??
     LOGICAL  :: growing_season                    !< logical to identify whether the plant is in the growing season or not
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
    REAL(wp) :: theta                 ! fractional water holding capacity
    REAL(wp) :: moisture              ! actual plant available water content (mm)
    REAL(wp) :: ammonium              ! soil ammonium concentration (mmol/L)
    REAL(wp) :: phosphate             ! soil P concentration (mmol/L)
  END TYPE SoilClass

  ! 
  ! Tile ( TileClass: integrates tile, plant & soil )
  ! 
  TYPE AccumulationTileClass
    ! definition of averages of variable across a certain period (model time steps)
    ! see mo_land_memory for initialization of default values
    ! Environmental forcing
    REAL(wp) :: t_air = 1._wp
    REAL(wp) :: q_air = 1._wp
    REAL(wp) :: press_surf = 1._wp
    ! memory of certain variables
    REAL(wp) :: t_air_week     = 273.15_wp   ! 7 days
    REAL(wp) :: t_air_month    = 273.15_wp   ! 30 days
    REAL(wp) :: t_air_longterm = 273.15_wp   ! 30 years == 10950 days
  END TYPE AccumulationTileClass

  TYPE TileClass
     TYPE(PlantClass)   :: plant
     TYPE(SoilClass)    :: soil
     TYPE(AccumulationTileClass) :: accumulated   !! <-  is going to be an own type with procedures soon -> Jan
  END TYPE TileClass  

  TYPE t_slm_grid           ! temporary ! -> will be stored in a /base/mo_grid_class later on (see JSBACH4)
     REAL(wp) :: lon
     REAL(wp) :: lat
  END TYPE t_slm_grid
  TYPE(t_slm_grid) :: grid  !< currently initialized with namelist parameter values

  TYPE(Element) :: element_zero   !< all variables are zero, is used for setting other elements zero frequently
! PUBLIC :: ElementAdd,ElementSub,ElementMul,ElementDiv
 
CONTAINS

  !-----------------------------------------------------------------------------------------------------
  !> initialize the derived TYPE Element element_zero with zero
  !! 
  !! this is done only one time, this element_zero is used to set other derived TYPE's of Element to zero
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE init_element_zero
     ! TYPE(Element) :: element_zero | defined above
     element_zero%carbon     = 0.0_wp
     element_zero%nitrogen   = 0.0_wp
     element_zero%phosphorus = 0.0_wp
     element_zero%carbon13   = 0.0_wp
     element_zero%carbon14   = 0.0_wp
     element_zero%nitrogen15 = 0.0_wp
  END SUBROUTINE init_element_zero

  !-----------------------------------------------------------------------------------------------------
  !> initialize all plant variables for pools and so, not fluxes
  !! 
  !! called in the driver; for fluxes see mo_veg_memory:set_zero_fluxes_of_plant_pools
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE init_plant(plant)
    IMPLICIT NONE
    ! 0.1 InOut
    TYPE(PlantClass), INTENT(inout) :: plant

    ! pools
    plant%pool%leaf         = element_zero   ! TYPE(Element) :: element_zero
    plant%pool%fine_root    = element_zero
    plant%pool%coarse_root  = element_zero
    plant%pool%sap_wood     = element_zero
    plant%pool%heart_wood   = element_zero
    plant%pool%labile       = element_zero
    plant%pool%reserve      = element_zero
    plant%pool%fruit        = element_zero
    plant%pool%seed_bed     = element_zero
    ! canopy layers
    plant%canopy_layer(:)%nl_tot       = 0.0_wp
    plant%canopy_layer(:)%fnl_chl      = 0.0_wp
    plant%canopy_layer(:)%fnl_et       = 0.0_wp
    plant%canopy_layer(:)%fnl_rub      = 0.0_wp
    plant%canopy_layer(:)%par_sunlit   = 0.0_wp
    plant%canopy_layer(:)%par_shaded   = 0.0_wp
    plant%canopy_layer(:)%fleaf_sunlit = 0.0_wp
    ! plant variables
    plant%pftid            = 1          ! INTEGER
    plant%height           = 0.0_wp
    plant%diameter         = 0.0_wp
    plant%lai              = 0.0_wp
    plant%target_cn_leaf   = 0.0_wp
    plant%target_np_leaf   = 0.0_wp
    plant%mean_leaf_age    = 0.0_wp
    plant%gc               = 0.0_wp
    plant%growing_season   = .TRUE.     ! LOGICAL
  END SUBROUTINE init_plant

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
