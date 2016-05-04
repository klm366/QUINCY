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

  !-----------------------------------------------------------------------------------------------------
  !> elements used across veg & soil of the model (C N P 13C 14C 15N)
  !! 
  !! this type is used for all tissues, pools, fluxes in vegetation and soil;
  !! changes made here affect the whole model
  !-----------------------------------------------------------------------------------------------------
  TYPE t_element 
     REAL(wp) :: carbon      !< total carbon concentration of pool in mol / m2
     REAL(wp) :: nitrogen    !< totol nitrogen concentration of pool in mol / m2
     REAL(wp) :: phosphorus  !< phosphorus concentration of pool in mol / m2
     REAL(wp) :: carbon13    !< 13carbon concentration of pool in mol / m2
     REAL(wp) :: carbon14    !< 14carbon concentration of pool in mol / m2
     REAL(wp) :: nitrogen15  !< 15nitrogen concentration of pool in mol / m2
  END TYPE t_element

  !-----------------------------------------------------------------------------------------------------
  !> vegetation pools storing elements in tissues (e.g. leaf, labile, reserve)
  !! 
  !! this type is used for vegetation pools, and fluxes in & out of pools;
  !! pools: there is one pool per plant and tissue
  !! flux: there are two fluxes per pool of a tissue, growth & litterfall (in & out of pool)
  !-----------------------------------------------------------------------------------------------------
  TYPE t_slm_veg_pool 
     TYPE(t_element) :: leaf           !< foliage
     TYPE(t_element) :: fine_root      !< fine roots (< 2mm diameter)
     TYPE(t_element) :: coarse_root    !< coarse roots (> 2mm diameter)
     TYPE(t_element) :: sap_wood       !< living wood for water and phloem transport
     TYPE(t_element) :: heart_wood     !< dead standing wood
     TYPE(t_element) :: labile         !< meristem carbon
     TYPE(t_element) :: reserve        !< storage
     TYPE(t_element) :: fruit          !< live propagation pool
     TYPE(t_element) :: seed_bed       !< fruit mass on ground
  END TYPE t_slm_veg_pool

  !-----------------------------------------------------------------------------------------------------
  !> all fluxes of vegetation (in & out of pools)
  !! 
  !! uses the TYPE t_slm_veg_pool to define fluxes: growth & litterfall \n
  !! plant\%flux\%growth & plant\%flux\%litterfall
  !! 
  !! includes additionally: TYPE(t_element) :: reserve_use \n
  !! plant\%flux\%reserve_use
  !-----------------------------------------------------------------------------------------------------
  TYPE t_slm_veg_flux 
     TYPE(t_element)      :: reserve_use  !< net flux between reserve and labile pool (mol m-2 timestep-1)
                                          !!  (+ve: flux to labile; -ve flux from labile) 
     TYPE(t_slm_veg_pool) :: growth       !< add to      plant%pool due to growth     per timestep (mol m-2 timestep-1)
     TYPE(t_slm_veg_pool) :: litterfall   !< remove from plant%pool due to litterfall per timestep (mol m-2 timestep-1)
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
     REAL(wp) :: uptake_n                 !< N uptake from soil ( is that correct? )
     REAL(wp) :: uptake_n15               !< 15N uptake from soil ( is that correct? ) 
     REAL(wp) :: uptake_p                 !< P uptake from soil ( is that correct? )
     REAL(wp) :: leaf_resorption_n        !< resorption of N upon leaf senescence (mol m-2 timestep-1)
     REAL(wp) :: leaf_resorption_n15      !< resorption of 15N upon leaf senescence (mol m-2 timestep-1)
     REAL(wp) :: leaf_resorption_p        !< resorption of P upon leaf senescence (mol m-2 timestep-1)
     REAL(wp) :: transpiration            !< transpiration in kg m-2 s-1
     REAL(wp) :: ci                       !< ??
  END TYPE t_slm_veg_flux

  !-----------------------------------------------------------------------------------------------------
  !> vegetation canopy structure (for each layer)
  !! 
  !! number of layers is: ncanop
  !-----------------------------------------------------------------------------------------------------
  TYPE t_slm_canopy_layer
     REAL(wp) :: nl_tot                !< N in canopy layer (mmol / m2 (LAI))
     REAL(wp) :: fnl_chl               !< fraction of N in chlorophyll (pigments only!)
     REAL(wp) :: fnl_et                !< fraction of N in electron transport 
     REAL(wp) :: fnl_rub               !< fraction of N in Rubisco 
     REAL(wp) :: par_sunlit            !< PAR on sunlit leaves (µmol m-2 s-1)
     REAL(wp) :: par_shaded            !< PAR on shaded leaves (µmol m-2 s-1)
     REAL(wp) :: fleaf_sunlit          !< fraction of layer which is sunlit
  END TYPE t_slm_canopy_layer

  !-----------------------------------------------------------------------------------------------------
  !> all variables of vegetation per tile
  !! 
  !! including all pools and all fluxes
  !-----------------------------------------------------------------------------------------------------
  TYPE t_slm_plant 
     TYPE(t_slm_veg_pool) :: pool                  !< the actual vegetation pools
     TYPE(t_slm_veg_flux) :: flux                  !< contains growth and litter fall fluxes per timestep for each pool
     TYPE(t_slm_canopy_layer),DIMENSION(ncanop) :: canopy_layer  !< vegetation variables per canopy layer
     INTEGER  :: pftid                             !< ID of the plant functional type; namelist parameter
     REAL(wp) :: height                            !< plant height in m
     REAL(wp) :: diameter                          !< plant diameter in m
     REAL(wp) :: lai                               !< actual leaf area index
     REAL(wp) :: target_cn_leaf                    !< target foliar C:N ratio
     REAL(wp) :: target_np_leaf                    !< target foliar N:P ratio 
     REAL(wp) :: mean_leaf_age                     !< average foliage age in days 
     REAL(wp) :: gc                                !< ??
     REAL(wp) :: f_n_demand                        !< relative plant demand for N; see update_vegsoil_bgc
     REAL(wp) :: f_p_demand                        !< relative plant demand for P; see update_vegsoil_bgc
     REAL(wp) :: beta_air                          !< scaling factor to account for air humidity effects on conductance
     REAL(wp) :: beta_soil_ps                      !< scaling factor to account for soil moisture constraints on conductance
     REAL(wp) :: beta_soil_gs                      !< scaling factor to account for soil moisture constraints on conductance
     LOGICAL  :: growing_season                    !< logical to identify whether the plant is in the growing season or not
   END TYPE t_slm_plant

  !-----------------------------------------------------------------------------------------------------
  !> all (or most ?) variables of soil per tile
  !! 
  !! including all pools and all fluxes \n
  !! @todo this will be extended by an extra soil model
  !-----------------------------------------------------------------------------------------------------
  TYPE t_slm_soil
    REAL(wp) :: theta                 !< fractional water holding capacity
    REAL(wp) :: moisture              !< actual plant available water content (mm)
    REAL(wp) :: ammonium              !< soil ammonium concentration (mmol/L)
    REAL(wp) :: phosphate             !< soil P concentration (mmol/L)
  END TYPE t_slm_soil

  !-----------------------------------------------------------------------------------------------------
  !> integrates all variables of one tile
  !!
  !! including plant and soil \n
  !! TYPE(t_slm_tile) :: tile --> declaration of tile in the driver \n
  !! note: during model simulations tiles may run parallel !
  !-----------------------------------------------------------------------------------------------------
  TYPE t_slm_tile
     TYPE(t_slm_plant)   :: plant
     TYPE(t_slm_soil)    :: soil
  END TYPE t_slm_tile  

  !-----------------------------------------------------------------------------------------------------
  !> longitude and latitude of a tile
  !!
  !! note: this TYPE is just a temporary solution \n
  !! migration to /base/mo_grid_class may follow, depending on structure in JSBACH4
  !-----------------------------------------------------------------------------------------------------
  TYPE t_slm_grid 
     REAL(wp) :: lon
     REAL(wp) :: lat
  END TYPE t_slm_grid

  TYPE(t_slm_grid)  :: grid           !< currently initialized with namelist parameter values
  TYPE(t_element)   :: element_zero   !< all variables are zero, is used for setting other elements zero frequently

! PUBLIC :: ElementAdd,ElementSub,ElementMul,ElementDiv
 
CONTAINS

  !-----------------------------------------------------------------------------------------------------
  !> initialize the derived TYPE t_element element_zero with zero
  !! 
  !! this is done only one time, this element_zero is used to set other derived TYPE's of t_element to zero
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE init_element_zero
     ! TYPE(t_element) :: element_zero | defined above
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
    TYPE(t_slm_plant), INTENT(inout) :: plant

    ! pools
    plant%pool%leaf         = element_zero   ! TYPE(t_element) :: element_zero
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
    plant%f_n_demand       = 0.0_wp
    plant%f_p_demand       = 0.0_wp
    plant%beta_air         = 0.0_wp
    plant%beta_soil_ps     = 0.0_wp
    plant%beta_soil_gs     = 0.0_wp
    plant%growing_season   = .TRUE.     ! LOGICAL
  END SUBROUTINE init_plant

  FUNCTION ElementAdd(e1,e2) RESULT(e3)
     IMPLICIT NONE
     TYPE(t_element), INTENT(in) :: e1,e2
     TYPE(t_element) :: e3
     e3%carbon = e1%carbon + e2%carbon 
     e3%nitrogen = e1%nitrogen + e2%nitrogen 
     e3%phosphorus = e1%phosphorus + e2%phosphorus 
     e3%carbon13 = e1%carbon13 + e2%carbon13 
     e3%carbon14 = e1%carbon14 + e2%carbon14
     e3%nitrogen15 = e1%nitrogen15 + e2%nitrogen15
  END FUNCTION ElementAdd  

  FUNCTION ElementSub(e1,e2) RESULT(e3)
     IMPLICIT NONE
     TYPE(t_element), INTENT(in) :: e1,e2
     TYPE(t_element) :: e3
     e3%carbon = e1%carbon - e2%carbon 
     e3%nitrogen = e1%nitrogen - e2%nitrogen 
     e3%phosphorus = e1%phosphorus - e2%phosphorus 
     e3%carbon13 = e1%carbon13 - e2%carbon13 
     e3%carbon14 = e1%carbon14 - e2%carbon14
     e3%nitrogen15 = e1%nitrogen15 - e2%nitrogen15
  END FUNCTION ElementSub  

  FUNCTION ElementMul(e1,e2) RESULT(e3)
     IMPLICIT NONE
     TYPE(t_element), INTENT(in) :: e1,e2
     TYPE(t_element) :: e3
     e3%carbon = e1%carbon * e2%carbon 
     e3%nitrogen = e1%nitrogen * e2%nitrogen 
     e3%phosphorus = e1%phosphorus * e2%phosphorus 
     e3%carbon13 = e1%carbon13 * e2%carbon13 
     e3%carbon14 = e1%carbon14 * e2%carbon14
     e3%nitrogen15 = e1%nitrogen15 * e2%nitrogen15
  END FUNCTION ElementMul

  FUNCTION ElementDiv(e1,e2) RESULT(e3)
     IMPLICIT NONE
     TYPE(t_element), INTENT(in) :: e1,e2
     TYPE(t_element) :: e3
     e3%carbon = e1%carbon / e2%carbon 
     e3%nitrogen = e1%nitrogen / e2%nitrogen 
     e3%phosphorus = e1%phosphorus / e2%phosphorus 
     e3%carbon13 = e1%carbon13 / e2%carbon13 
     e3%carbon14 = e1%carbon14 / e2%carbon14
     e3%nitrogen15 = e1%nitrogen15 / e2%nitrogen15
  END FUNCTION ElementDiv  

END MODULE
