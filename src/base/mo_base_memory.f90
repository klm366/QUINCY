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

  USE mo_kind,           ONLY: wp
  USE mo_veg_memory,     ONLY: t_slm_plant
  USE mo_shared_memory,    ONLY: t_element

  IMPLICIT NONE

  !-----------------------------------------------------------------------------------------------------
  !> soil pool
  !! 
  !-----------------------------------------------------------------------------------------------------
  TYPE t_slm_soil_pool
     TYPE(t_element)  :: litter         !< soil litter
     TYPE(t_element)  :: som            !< soil organic mater
  END TYPE t_slm_soil_pool

  !-----------------------------------------------------------------------------------------------------
  !> soil flux
  !! 
  !-----------------------------------------------------------------------------------------------------
  TYPE t_slm_soil_flux
     TYPE(t_slm_soil_pool) :: loss      !< decomposition flux from each pool
  END TYPE t_slm_soil_flux

  !-----------------------------------------------------------------------------------------------------
  !> all (or most ?) variables of soil per tile
  !! 
  !! including all pools and all fluxes \n
  !! @todo this will be extended by an extra soil model
  !-----------------------------------------------------------------------------------------------------
  TYPE t_slm_soil
    TYPE(t_slm_soil_pool) :: pool     !< soil pool
    TYPE(t_slm_soil_flux) :: flux     !< soil flux
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

 
!CONTAINS



END MODULE mo_base_memory
