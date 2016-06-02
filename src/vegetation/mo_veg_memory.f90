!>
!! @brief Contains the routines for plant memory operations 
!!
!! @author
!!  Sönke Zaehle, MPI-BGC Jena 
!!
!! @par Revision History
!! First version by Sönke Zaehle (2015-12-04)
!!
MODULE mo_veg_memory

  USE mo_kind,            ONLY: wp, eps8
  USE mo_base_constants,  ONLY: ncanop, dt, one_day
  USE mo_shared_memory,     ONLY: t_element

  IMPLICIT NONE

  !-----------------------------------------------------------------------------------------------------
  !> Type for storing and calculating averages per day for multiple variables
  !!
  !! is part of tile\%plant now and declared there
  !!
  !! includes the PROCEDURE's necessary to do so \n 
  !! makes the average of the previous day globally available
  !-----------------------------------------------------------------------------------------------------
  TYPE :: t_daily_average
    ! the 1st index of an array can be addressed with 1 (i.e. not 0 !)
    ! averages
    REAL(wp), DIMENSION(ncanop)     :: avrg_prev_par_sunlit_day      !<  average of the previous day
    REAL(wp), DIMENSION(ncanop)     :: avrg_prev_par_shaded_day      !<  average of the previous day
    REAL(wp), DIMENSION(ncanop)     :: avrg_prev_fleaf_sunlit_day    !<  average of the previous day
    ! accumulated value of the current day & number of values accumulated
    REAL(wp), DIMENSION(ncanop)     :: curr_val_par_sunlit_day        !<  (1) value accumulated over current day
    REAL(wp), DIMENSION(ncanop)     :: curr_val_par_shaded_day        !<  (1) value accumulated over current day
    REAL(wp), DIMENSION(ncanop)     :: curr_val_fleaf_sunlit_day      !<  (1) value accumulated over current day
    REAL(wp), DIMENSION(ncanop)     :: curr_count_par_sunlit_day      !<  (2) counter of acc values over current day
    REAL(wp), DIMENSION(ncanop)     :: curr_count_par_shaded_day      !<  (2) counter of acc values over current day
    REAL(wp), DIMENSION(ncanop)     :: curr_count_fleaf_sunlit_day    !<  (2) counter of acc values over current day
    ! condition per variable and timestep
    LOGICAL                         :: l_par_sunlit_day              !< whether sunlight is available at canopy top
    LOGICAL                         :: l_par_shaded_day              !< whether sunlight is available at canopy top
    LOGICAL                         :: l_fleaf_sunlit_day            !< whether sunlight is available at canopy top
  CONTAINS
    PROCEDURE :: init              => init_daily_avrg_variables
    PROCEDURE :: update_conditions => update_conditions_current_timestep
    PROCEDURE :: accumulate        => accumulate_values_over_current_day
    PROCEDURE :: calc_avrg         => calculate_avrg_previous_day  ! at the end of calculations of the last timestep per day
  END TYPE t_daily_average

  !-----------------------------------------------------------------------------------------------------
  !> Type for storing and calculating the moving averages for multiple variables across the model
  !!
  !! is part of tile\%plant now and declared there
  !!
  !! initial default values and accumulation period per variable are set
  !!
  !! includes the PROCEDURE's necessary to do so \n
  !! makes the moving averages globally available
  !-----------------------------------------------------------------------------------------------------
  TYPE t_moving_average
    !< 1.0 variables that store averages
    !! 
    ! 1.1 light in the canopy
    REAL(wp),DIMENSION(ncanop) :: par_sunlit
    REAL(wp),DIMENSION(ncanop) :: par_shaded
    REAL(wp),DIMENSION(ncanop) :: fleaf_sunlit
    ! 1.2 soil moisture (currently in mo_veg_photosynthesis)
    REAL(wp) :: beta_air                     !< not used (scaling factor to account for air humidity effects on conductance)
    REAL(wp) :: beta_soil_ps                 !< not used (scaling factor to account for soil moisture constraints on conductance)
    REAL(wp) :: beta_soil_gs                 !< scaling factor to account for soil moisture constraints on conductance
    ! 1.3 memory of fluxes
    REAL(wp) :: gpp_week                                   !< gpp of the last week (µmol m-2 s-1)
    REAL(wp) :: maint_respiration_week                     !< maintenance respiration of the last week (µmol m-2 s-1)
    ! 1.4 nutrient demand for plant uptake (average demand over a couple of days to get rid of diurnal C cycle)
    REAL(wp) :: demand_uptake_n                            !< scalar of Vmax for root N uptake 
    REAL(wp) :: demand_uptake_p                            !< scalar of Vmax for root P uptake
    REAL(wp) :: growth_req_n_week         !< moles N required for a unit of C growth under current allocation fractions
    REAL(wp) :: growth_req_p_week         !< moles P required for a unit of N growth under current allocation fractions
    ! 1.5 environmental forcing
    REAL(wp) :: t_air             ! = 1._wp
    REAL(wp) :: q_air             ! = 1._wp
    REAL(wp) :: press_surf        ! = 1._wp
    ! 1.6 memory of certain variables, temperature, ...
    REAL(wp) :: t_air_week         ! = 273.15_wp   ! 7 days
    REAL(wp) :: t_air_month        ! = 273.15_wp   ! 30 days
    REAL(wp) :: t_air_longterm     ! = 273.15_wp   ! 30 years == 10950 days

    !< 2.0 accumulation period per variable in days !!
    !! 
    ! 2.0a some default values
    REAL(wp) :: p_weekly                 !< 7 days; constant across the model
    REAL(wp) :: p_monthly                !< 30 days; constant across the model
    REAL(wp) :: p_longterm               !< 10950 days (30 years); constant across the model
    ! 2.1 light in the canopy (for all canopy_layer the same period)
    REAL(wp) :: p_par_sunlit 
    REAL(wp) :: p_par_shaded 
    REAL(wp) :: p_fleaf_sunlit 
    ! 2.2 soil moisture
    REAL(wp) :: p_beta_air 
    REAL(wp) :: p_beta_soil_ps 
    REAL(wp) :: p_beta_soil_gs 
    ! 2.3 memory of fluxes are averaged weekly, thus, p_weekly is used

    ! 2.4 nutrient demand for plant uptake
    REAL(wp) :: p_demand_uptake_n 
    REAL(wp) :: p_demand_uptake_p 
    ! 2.5 environmental forcing
    REAL(wp) :: p_t_air 
    REAL(wp) :: p_q_air 
    REAL(wp) :: p_press_surf 
    ! 2.6 memory of certain variables: default averaging period values are used, according to variable names (week, month, ..)

  CONTAINS
    !< 4.0 procedures
    !! 
    PROCEDURE :: init                => init_moving_avrg_variables
    PROCEDURE :: calc_moving_avrg    => calculate_moving_average_of_variable
    PROCEDURE :: update_veg          => update_veg_moving_avrg_variables
    PROCEDURE :: update_soil         => update_soil_moving_avrg_variables
    PROCEDURE :: update_tile         => update_tile_moving_avrg_variables
  END TYPE t_moving_average

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
  CONTAINS
    PRIVATE
    PROCEDURE :: add_veg_pools_type_type
    PROCEDURE :: add_veg_pools_type_scalar
    PROCEDURE, PASS(this2) :: add_veg_pools_scalar_type
    GENERIC, PUBLIC :: OPERATOR(+) => add_veg_pools_type_type, add_veg_pools_type_scalar, add_veg_pools_scalar_type
    PROCEDURE :: subtract_veg_pools_type_type
    PROCEDURE :: subtract_veg_pools_type_scalar
    PROCEDURE, PASS(this2) :: subtract_veg_pools_scalar_type
    GENERIC, PUBLIC :: OPERATOR(-) => subtract_veg_pools_type_type, subtract_veg_pools_type_scalar, subtract_veg_pools_scalar_type
    PROCEDURE :: multiply_veg_pools_type_type
    PROCEDURE :: multiply_veg_pools_type_scalar
    PROCEDURE, PASS(this2) :: multiply_veg_pools_scalar_type
    GENERIC, PUBLIC :: OPERATOR(*) => multiply_veg_pools_type_type, multiply_veg_pools_type_scalar, multiply_veg_pools_scalar_type
    PROCEDURE :: divide_veg_pools_type_type
    PROCEDURE :: divide_veg_pools_type_scalar
    PROCEDURE, PASS(this2) :: divide_veg_pools_scalar_type
    GENERIC, PUBLIC :: OPERATOR(/) => divide_veg_pools_type_type, divide_veg_pools_type_scalar, divide_veg_pools_scalar_type
    PROCEDURE, PUBLIC, PASS(this2) :: init => init_veg_pool
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
     REAL(wp) :: maint_respiration_pot    !< potential maintenance respiration (µmol CO2 m-2 s-1)
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
     REAL(wp) :: resorption_n             !< resorption of N upon tissue senescence (mol m-2 timestep-1)
     REAL(wp) :: resorption_n15           !< resorption of 15N upon tissue senescence (mol m-2 timestep-1)
     REAL(wp) :: resorption_p             !< resorption of P upon tissue senescence (mol m-2 timestep-1)
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
     REAL(wp) :: fnl_pepc              !< fraction of N in PEP carboylase (C4 plants only)
     REAL(wp) :: fnl_oth               !< fraction of N not associated with photosynthesis 
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
    TYPE(t_slm_veg_pool)   :: pool                !< the actual vegetation pools
    TYPE(t_slm_veg_flux)   :: flux                !< contains growth and litter fall fluxes per timestep for each pool
    TYPE(t_slm_canopy_layer),DIMENSION(ncanop) :: canopy_layer  !< vegetation variables per canopy layer
    TYPE(t_daily_average)  :: daily_average       !< structure to calculate and store all moving average values
    TYPE(t_moving_average) :: moving_average      !< structure to calculate and store all moving average values
    INTEGER  :: pftid                             !< ID of the plant functional type; namelist parameter
    REAL(wp) :: height                            !< plant height in m
    REAL(wp) :: diameter                          !< plant diameter in m
    REAL(wp) :: lai                               !< actual leaf area index
    REAL(wp) :: target_cn_leaf                    !< target foliar C:N ratio
    REAL(wp) :: target_np_leaf                    !< target foliar N:P ratio 
    REAL(wp) :: mean_leaf_age                     !< average foliage age in days 
    REAL(wp) :: gc                                !< some kind of stomatal conductance ?
    REAL(wp) :: f_n_demand                        !< relative plant demand for N; see update_vegsoil_bgc
    REAL(wp) :: f_p_demand                        !< relative plant demand for P; see update_vegsoil_bgc
    REAL(wp) :: beta_air                          !< scaling factor to account for air humidity effects on conductance
    REAL(wp) :: beta_soil_ps                      !< scaling factor to account for soil moisture constraints on conductance
    REAL(wp) :: beta_soil_gs                      !< scaling factor to account for soil moisture constraints on conductance
    LOGICAL  :: growing_season                    !< logical to identify whether the plant is in the growing season or not
    REAL(wp) :: growth_req_n                      !< moles N required for a unit of C growth under current allocation fractions 
    REAL(wp) :: growth_req_p                      !< moles P required for a unit of N growth under current allocation fractions 
    REAL(wp) :: gdd                               !< growing degree days above 5°C
    REAL(wp) :: nd_dormance                       !< number of days since last growing season 
   END TYPE t_slm_plant


CONTAINS

  !-----------------------------------------------------------------------------------------------------
  !> 1.0 addition of two vegetation pools: type + type
  !! 
  !-----------------------------------------------------------------------------------------------------
  PURE ELEMENTAL FUNCTION add_veg_pools_type_type(this,t2) RESULT(t3) 

    USE mo_shared_memory,   ONLY: t_element

     IMPLICIT NONE
     CLASS(t_slm_veg_pool), INTENT(in) :: this, t2         !< any t_slm_veg_pool
     TYPE(t_slm_veg_pool)              :: t3               !< any t_slm_veg_pool
     t3%leaf          = this%leaf        + t2%leaf         ! operator overload, see t_element
     t3%fine_root     = this%fine_root   + t2%fine_root
     t3%coarse_root   = this%coarse_root + t2%coarse_root
     t3%sap_wood      = this%sap_wood    + t2%sap_wood
     t3%heart_wood    = this%heart_wood  + t2%heart_wood
     t3%labile        = this%labile      + t2%labile
     t3%reserve       = this%reserve     + t2%reserve
     t3%fruit         = this%fruit       + t2%fruit
     t3%seed_bed      = this%seed_bed    + t2%seed_bed
  END FUNCTION add_veg_pools_type_type  

  !-----------------------------------------------------------------------------------------------------
  !> 1.1 addition of two vegetation pools: type + scalar
  !! 
  !-----------------------------------------------------------------------------------------------------
  PURE ELEMENTAL FUNCTION add_veg_pools_type_scalar(this,s2) RESULT(t3) 

    USE mo_shared_memory,   ONLY: t_element

     IMPLICIT NONE
     CLASS(t_slm_veg_pool), INTENT(in) :: this        !< any t_slm_veg_pool
     REAL(wp), INTENT(in)              :: s2          !< any scalar
     TYPE(t_slm_veg_pool)              :: t3
     t3%leaf          = this%leaf        + s2         ! operator overload, see t_element
     t3%fine_root     = this%fine_root   + s2
     t3%coarse_root   = this%coarse_root + s2
     t3%sap_wood      = this%sap_wood    + s2
     t3%heart_wood    = this%heart_wood  + s2
     t3%labile        = this%labile      + s2
     t3%reserve       = this%reserve     + s2
     t3%fruit         = this%fruit       + s2
     t3%seed_bed      = this%seed_bed    + s2
  END FUNCTION add_veg_pools_type_scalar  

  !-----------------------------------------------------------------------------------------------------
  !> 1.2 addition of two vegetation pools: scalar + type
  !! 
  !-----------------------------------------------------------------------------------------------------
  PURE ELEMENTAL FUNCTION add_veg_pools_scalar_type(s1,this2) RESULT(t3) 

    USE mo_shared_memory,   ONLY: t_element

     IMPLICIT NONE
     REAL(wp), INTENT(in)              :: s1          !< any scalar
     CLASS(t_slm_veg_pool), INTENT(in) :: this2       !< any t_slm_veg_pool
     TYPE(t_slm_veg_pool)              :: t3
     t3%leaf          = s1 + this2%leaf               ! operator overload, see t_element
     t3%fine_root     = s1 + this2%fine_root 
     t3%coarse_root   = s1 + this2%coarse_root
     t3%sap_wood      = s1 + this2%sap_wood  
     t3%heart_wood    = s1 + this2%heart_wood
     t3%labile        = s1 + this2%labile 
     t3%reserve       = s1 + this2%reserve 
     t3%fruit         = s1 + this2%fruit 
     t3%seed_bed      = s1 + this2%seed_bed 
  END FUNCTION add_veg_pools_scalar_type  

  !-----------------------------------------------------------------------------------------------------
  !> 2.0 subtraction of two vegetation pools: type - type
  !! 
  !-----------------------------------------------------------------------------------------------------
  PURE ELEMENTAL FUNCTION subtract_veg_pools_type_type(this,t2) RESULT(t3) 

    USE mo_shared_memory,   ONLY: t_element

     IMPLICIT NONE
     CLASS(t_slm_veg_pool), INTENT(in) :: this, t2         !< any t_slm_veg_pool
     TYPE(t_slm_veg_pool)              :: t3               !< any t_slm_veg_pool
     t3%leaf          = this%leaf        - t2%leaf         ! operator overload, see t_element
     t3%fine_root     = this%fine_root   - t2%fine_root
     t3%coarse_root   = this%coarse_root - t2%coarse_root
     t3%sap_wood      = this%sap_wood    - t2%sap_wood
     t3%heart_wood    = this%heart_wood  - t2%heart_wood
     t3%labile        = this%labile      - t2%labile
     t3%reserve       = this%reserve     - t2%reserve
     t3%fruit         = this%fruit       - t2%fruit
     t3%seed_bed      = this%seed_bed    - t2%seed_bed
  END FUNCTION subtract_veg_pools_type_type  

  !-----------------------------------------------------------------------------------------------------
  !> 2.1 subtraction of two vegetation pools: type - scalar
  !! 
  !-----------------------------------------------------------------------------------------------------
  PURE ELEMENTAL FUNCTION subtract_veg_pools_type_scalar(this,s2) RESULT(t3) 

    USE mo_shared_memory,   ONLY: t_element

     IMPLICIT NONE
     CLASS(t_slm_veg_pool), INTENT(in) :: this        !< any t_slm_veg_pool
     REAL(wp), INTENT(in)              :: s2          !< any scalar
     TYPE(t_slm_veg_pool)              :: t3
     t3%leaf          = this%leaf        - s2         ! operator overload, see t_element
     t3%fine_root     = this%fine_root   - s2
     t3%coarse_root   = this%coarse_root - s2
     t3%sap_wood      = this%sap_wood    - s2
     t3%heart_wood    = this%heart_wood  - s2
     t3%labile        = this%labile      - s2
     t3%reserve       = this%reserve     - s2
     t3%fruit         = this%fruit       - s2
     t3%seed_bed      = this%seed_bed    - s2
  END FUNCTION subtract_veg_pools_type_scalar  

  !-----------------------------------------------------------------------------------------------------
  !> 2.2 subtraction of two vegetation pools: scalar - type
  !! 
  !-----------------------------------------------------------------------------------------------------
  PURE ELEMENTAL FUNCTION subtract_veg_pools_scalar_type(s1,this2) RESULT(t3) 

    USE mo_shared_memory,   ONLY: t_element

     IMPLICIT NONE
     REAL(wp), INTENT(in)              :: s1          !< any scalar
     CLASS(t_slm_veg_pool), INTENT(in) :: this2       !< any t_slm_veg_pool
     TYPE(t_slm_veg_pool)              :: t3
     t3%leaf          = s1 - this2%leaf               ! operator overload, see t_element
     t3%fine_root     = s1 - this2%fine_root 
     t3%coarse_root   = s1 - this2%coarse_root
     t3%sap_wood      = s1 - this2%sap_wood  
     t3%heart_wood    = s1 - this2%heart_wood
     t3%labile        = s1 - this2%labile 
     t3%reserve       = s1 - this2%reserve 
     t3%fruit         = s1 - this2%fruit 
     t3%seed_bed      = s1 - this2%seed_bed 
  END FUNCTION subtract_veg_pools_scalar_type  

  !-----------------------------------------------------------------------------------------------------
  !> 3.0 multiplication of two vegetation pools: type * type
  !! 
  !-----------------------------------------------------------------------------------------------------
  PURE ELEMENTAL FUNCTION multiply_veg_pools_type_type(this,t2) RESULT(t3) 

    USE mo_shared_memory,   ONLY: t_element

     IMPLICIT NONE
     CLASS(t_slm_veg_pool), INTENT(in) :: this, t2         !< any t_slm_veg_pool
     TYPE(t_slm_veg_pool)              :: t3               !< any t_slm_veg_pool
     t3%leaf          = this%leaf        * t2%leaf         ! operator overload, see t_element
     t3%fine_root     = this%fine_root   * t2%fine_root
     t3%coarse_root   = this%coarse_root * t2%coarse_root
     t3%sap_wood      = this%sap_wood    * t2%sap_wood
     t3%heart_wood    = this%heart_wood  * t2%heart_wood
     t3%labile        = this%labile      * t2%labile
     t3%reserve       = this%reserve     * t2%reserve
     t3%fruit         = this%fruit       * t2%fruit
     t3%seed_bed      = this%seed_bed    * t2%seed_bed
  END FUNCTION multiply_veg_pools_type_type  

  !-----------------------------------------------------------------------------------------------------
  !> 3.1 multiplication of two vegetation pools: type * scalar
  !! 
  !-----------------------------------------------------------------------------------------------------
  PURE ELEMENTAL FUNCTION multiply_veg_pools_type_scalar(this,s2) RESULT(t3) 

    USE mo_shared_memory,   ONLY: t_element

     IMPLICIT NONE
     CLASS(t_slm_veg_pool), INTENT(in) :: this        !< any t_slm_veg_pool
     REAL(wp), INTENT(in)              :: s2          !< any scalar
     TYPE(t_slm_veg_pool)              :: t3
     t3%leaf          = this%leaf        * s2         ! operator overload, see t_element
     t3%fine_root     = this%fine_root   * s2
     t3%coarse_root   = this%coarse_root * s2
     t3%sap_wood      = this%sap_wood    * s2
     t3%heart_wood    = this%heart_wood  * s2
     t3%labile        = this%labile      * s2
     t3%reserve       = this%reserve     * s2
     t3%fruit         = this%fruit       * s2
     t3%seed_bed      = this%seed_bed    * s2
  END FUNCTION multiply_veg_pools_type_scalar  

  !-----------------------------------------------------------------------------------------------------
  !> 3.2 multiplication of two vegetation pools: scalar * type
  !! 
  !-----------------------------------------------------------------------------------------------------
  PURE ELEMENTAL FUNCTION multiply_veg_pools_scalar_type(s1,this2) RESULT(t3) 

    USE mo_shared_memory,   ONLY: t_element

     IMPLICIT NONE
     REAL(wp), INTENT(in)              :: s1          !< any scalar
     CLASS(t_slm_veg_pool), INTENT(in) :: this2       !< any t_slm_veg_pool
     TYPE(t_slm_veg_pool)              :: t3
     t3%leaf          = s1 * this2%leaf               ! operator overload, see t_element
     t3%fine_root     = s1 * this2%fine_root 
     t3%coarse_root   = s1 * this2%coarse_root
     t3%sap_wood      = s1 * this2%sap_wood  
     t3%heart_wood    = s1 * this2%heart_wood
     t3%labile        = s1 * this2%labile 
     t3%reserve       = s1 * this2%reserve 
     t3%fruit         = s1 * this2%fruit 
     t3%seed_bed      = s1 * this2%seed_bed 
  END FUNCTION multiply_veg_pools_scalar_type  

  !-----------------------------------------------------------------------------------------------------
  !> 4.0 division of two vegetation pools: type / type
  !! 
  !-----------------------------------------------------------------------------------------------------
  PURE ELEMENTAL FUNCTION divide_veg_pools_type_type(this,t2) RESULT(t3) 

    USE mo_shared_memory,   ONLY: t_element

     IMPLICIT NONE
     CLASS(t_slm_veg_pool), INTENT(in) :: this, t2         !< any t_slm_veg_pool
     TYPE(t_slm_veg_pool)              :: t3               !< any t_slm_veg_pool
     t3%leaf          = this%leaf        / t2%leaf         ! operator overload, see t_element
     t3%fine_root     = this%fine_root   / t2%fine_root
     t3%coarse_root   = this%coarse_root / t2%coarse_root
     t3%sap_wood      = this%sap_wood    / t2%sap_wood
     t3%heart_wood    = this%heart_wood  / t2%heart_wood
     t3%labile        = this%labile      / t2%labile
     t3%reserve       = this%reserve     / t2%reserve
     t3%fruit         = this%fruit       / t2%fruit
     t3%seed_bed      = this%seed_bed    / t2%seed_bed
  END FUNCTION divide_veg_pools_type_type  

  !-----------------------------------------------------------------------------------------------------
  !> 4.1 division of two vegetation pools: type / scalar
  !! 
  !-----------------------------------------------------------------------------------------------------
  PURE ELEMENTAL FUNCTION divide_veg_pools_type_scalar(this,s2) RESULT(t3) 

    USE mo_shared_memory,   ONLY: t_element

     IMPLICIT NONE
     CLASS(t_slm_veg_pool), INTENT(in) :: this        !< any t_slm_veg_pool
     REAL(wp), INTENT(in)              :: s2          !< any scalar
     TYPE(t_slm_veg_pool)              :: t3
     t3%leaf          = this%leaf        / s2         ! operator overload, see t_element
     t3%fine_root     = this%fine_root   / s2
     t3%coarse_root   = this%coarse_root / s2
     t3%sap_wood      = this%sap_wood    / s2
     t3%heart_wood    = this%heart_wood  / s2
     t3%labile        = this%labile      / s2
     t3%reserve       = this%reserve     / s2
     t3%fruit         = this%fruit       / s2
     t3%seed_bed      = this%seed_bed    / s2
  END FUNCTION divide_veg_pools_type_scalar  

  !-----------------------------------------------------------------------------------------------------
  !> 4.2 division of two vegetation pools: scalar / type
  !! 
  !-----------------------------------------------------------------------------------------------------
  PURE ELEMENTAL FUNCTION divide_veg_pools_scalar_type(s1,this2) RESULT(t3) 

    USE mo_shared_memory,   ONLY: t_element

     IMPLICIT NONE
     REAL(wp), INTENT(in)              :: s1          !< any scalar
     CLASS(t_slm_veg_pool), INTENT(in) :: this2       !< any t_slm_veg_pool
     TYPE(t_slm_veg_pool)              :: t3
     t3%leaf          = s1 / this2%leaf               ! operator overload, see t_element
     t3%fine_root     = s1 / this2%fine_root 
     t3%coarse_root   = s1 / this2%coarse_root
     t3%sap_wood      = s1 / this2%sap_wood  
     t3%heart_wood    = s1 / this2%heart_wood
     t3%labile        = s1 / this2%labile 
     t3%reserve       = s1 / this2%reserve 
     t3%fruit         = s1 / this2%fruit 
     t3%seed_bed      = s1 / this2%seed_bed 
  END FUNCTION divide_veg_pools_scalar_type  

  !-----------------------------------------------------------------------------------------------------
  !> 5.0 initialise veg pool with scalar 
  !! 
  !-----------------------------------------------------------------------------------------------------
  PURE ELEMENTAL SUBROUTINE init_veg_pool(this2,s1)

    USE mo_shared_memory,   ONLY: t_element

    IMPLICIT NONE
    REAL(wp), INTENT(in)         :: s1
    CLASS(t_slm_veg_pool), INTENT(inout) :: this2 
     CALL this2%leaf%init(s1)
     CALL this2%fine_root%init(s1)
     CALL this2%coarse_root%init(s1)
     CALL this2%sap_wood%init(s1)
     CALL this2%heart_wood%init(s1)
     CALL this2%labile%init(s1)
     CALL this2%reserve%init(s1)
     CALL this2%fruit%init(s1)
     CALL this2%seed_bed%init(s1)
  END SUBROUTINE init_veg_pool 

  !-----------------------------------------------------------------------------------------------------
  !> intialize all variables (TYPE t_daily_average)
  !! 
  !! just set to zero or FALSE, called in the driver
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE init_daily_avrg_variables (this)
    IMPLICIT NONE 
    CLASS(t_daily_average), INTENT(inout) :: this

    this%avrg_prev_par_sunlit_day    = 0.0_wp    ! this sets the complete array to 0.0_wp
    this%avrg_prev_par_shaded_day    = 0.0_wp
    this%avrg_prev_fleaf_sunlit_day  = 0.0_wp

    this%curr_val_par_sunlit_day      = 0.0_wp
    this%curr_val_par_shaded_day      = 0.0_wp
    this%curr_val_fleaf_sunlit_day    = 0.0_wp
    this%curr_count_par_sunlit_day    = 0.0_wp
    this%curr_count_par_shaded_day    = 0.0_wp
    this%curr_count_fleaf_sunlit_day  = 0.0_wp

    this%l_par_sunlit_day    = .FALSE.
    this%l_par_shaded_day    = .FALSE.
    this%l_fleaf_sunlit_day  = .FALSE.
  END SUBROUTINE init_daily_avrg_variables

  !-----------------------------------------------------------------------------------------------------
  !> definition and check of conditions for each variable (TYPE t_daily_average)
  !! 
  !! for each variable a condition is defined whether the value of the current time step is used for
  !! average calculation or not;
  !! during model-runtime this subroutine is used to check the condition and set the particular
  !! variable TRUE or FALSE, valid for one timestep
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE update_conditions_current_timestep (this, plant)
    USE mo_atmland_memory, ONLY: atm2land

    IMPLICIT NONE 
    TYPE(t_slm_plant), INTENT(in) :: plant          !< mo_veg_memory: TYPE(t_slm_plant) :: plant; initialized in TYPE t_slm_tile
    CLASS(t_daily_average), INTENT(inout) :: this   !< TYPE t_daily_average

    ! condition type 1: light > zero
    IF(atm2land%swpar_srf_down > eps8) THEN
      this%l_par_sunlit_day    = .TRUE.
      this%l_par_shaded_day    = .TRUE.
      this%l_fleaf_sunlit_day  = .TRUE.
    ELSE
      this%l_par_sunlit_day    = .FALSE.
      this%l_par_shaded_day    = .FALSE.
      this%l_fleaf_sunlit_day  = .FALSE.      
    END IF
    ! condition type 2:
        ! ... currently none
  END SUBROUTINE update_conditions_current_timestep

  !-----------------------------------------------------------------------------------------------------
  !> accumulation of values over one day (TYPE t_daily_average)
  !! 
  !! per variable the values are summed up and the number of summed values is stored
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE accumulate_values_over_current_day (this, plant)
    IMPLICIT NONE 
    ! 0.1 InOut
    TYPE(t_slm_plant), INTENT(in) :: plant         !< mo_veg_memory: TYPE(t_slm_plant) :: plant; initialized in TYPE t_slm_tile
    CLASS(t_daily_average), INTENT(inout) :: this  !< TYPE t_daily_average
    ! 0.2 Local
    INTEGER  :: i                                  ! loop over canopy layers

    CALL this%update_conditions(plant)

    IF(this%l_par_sunlit_day) THEN
      DO i=1,ncanop
        this%curr_val_par_sunlit_day(i) = this%curr_val_par_sunlit_day(i) + plant%canopy_layer(i)%par_sunlit
        this%curr_count_par_sunlit_day(i) = this%curr_count_par_sunlit_day(i) + 1.0_wp
      END DO
    END IF

    IF(this%l_par_shaded_day) THEN
      DO i=1,ncanop
        this%curr_val_par_shaded_day(i) = this%curr_val_par_shaded_day(i) + plant%canopy_layer(i)%par_shaded
        this%curr_count_par_shaded_day(i) = this%curr_count_par_shaded_day(i) + 1.0_wp
      END DO
    END IF

    IF(this%l_fleaf_sunlit_day) THEN
      DO i=1,ncanop
        this%curr_val_fleaf_sunlit_day(i) = this%curr_val_fleaf_sunlit_day(i) + plant%canopy_layer(i)%fleaf_sunlit
        this%curr_count_fleaf_sunlit_day(i) = this%curr_count_fleaf_sunlit_day(i) + 1.0_wp
      END DO
    END IF

  END SUBROUTINE accumulate_values_over_current_day

  !-----------------------------------------------------------------------------------------------------
  !> calc the average of variables for the previous day (TYPE t_daily_average)
  !! 
  !! calculates the average at the end of the day (end of last time step of day), and sets all other
  !! variables to zero then
  !! 
  !! is called in the driver with: CALL tile(itile)\%plant\%daily_average\%calc_avrg
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE calculate_avrg_previous_day (this, plant, time_day, timestep_prop_day)
    IMPLICIT NONE 
    ! 0.1 InOut
    TYPE(t_slm_plant), INTENT(in) :: plant        !< mo_veg_memory: TYPE(t_slm_plant) :: plant; initialized in TYPE t_slm_tile
    REAL(wp), INTENT(in)          :: time_day           !< time in days
    REAL(wp), INTENT(in)          :: timestep_prop_day  !< length of one time step (proportion of one day)
    CLASS(t_daily_average), INTENT(inout) :: this       !< TYPE t_daily_average
    ! 0.2 Local
    INTEGER  :: i                                   ! loop over canopy layers

    ! the accumulation over the current day
    CALL this%accumulate(plant)   ! calls this%update_conditions(plant)

    ! if this is the last timestep of the current day (i.e. next timestep belongs to the next day):
    ! actually it is tested if the next timestep would be midnight: next_step(time_day) - round_to_int(time_day) < eps8
    IF( ABS(((time_day + timestep_prop_day) - REAL(NINT(time_day)))) < eps8 ) THEN
      ! calculate the average of the current day and pass it to the according avrg_prev_* variable; 
      ! check if at least one value had been assigned at the current day, avoiding division by zero 
      DO i=1,ncanop
        ! avrg_prev_par_sunlit_day
        IF( this%curr_count_par_sunlit_day(i) < eps8 ) THEN
          this%avrg_prev_par_sunlit_day(i)    = 0.0_wp
          ELSE
          this%avrg_prev_par_sunlit_day(i)    = this%curr_val_par_sunlit_day(i) / this%curr_count_par_sunlit_day(i)
        END IF
        ! avrg_prev_par_shaded_day
        IF( this%curr_count_par_shaded_day(i) < eps8 ) THEN
          this%avrg_prev_par_shaded_day(i)    = 0.0_wp
          ELSE
          this%avrg_prev_par_shaded_day(i)    = this%curr_val_par_shaded_day(i)   / this%curr_count_par_shaded_day(i)
        END IF
        ! avrg_prev_fleaf_sunlit_day
        IF( this%curr_count_fleaf_sunlit_day(i) < eps8 ) THEN
          this%avrg_prev_fleaf_sunlit_day(i)  = 0.0_wp
          ELSE
          this%avrg_prev_fleaf_sunlit_day(i)  = this%curr_val_fleaf_sunlit_day(i) / this%curr_count_fleaf_sunlit_day(i)
        END IF
      END DO
      
      !write(*,*) time_day, &   ! may enable for testing
      !           NINT(this%curr_count_par_sunlit_day(1)), this%avrg_prev_par_sunlit_day(1), &
      !           NINT(this%curr_count_par_shaded_day(1)), this%avrg_prev_par_shaded_day(1), &
      !           NINT(this%curr_count_fleaf_sunlit_day(1)), this%avrg_prev_fleaf_sunlit_day(1)

      ! set all values of the current day to zero (code can be copied from initialization_of_variables above)
      this%curr_val_par_sunlit_day      = 0.0_wp
      this%curr_val_par_shaded_day      = 0.0_wp
      this%curr_val_fleaf_sunlit_day    = 0.0_wp
      this%curr_count_par_sunlit_day    = 0.0_wp
      this%curr_count_par_shaded_day    = 0.0_wp
      this%curr_count_fleaf_sunlit_day  = 0.0_wp
    END IF

  END SUBROUTINE calculate_avrg_previous_day

  !-----------------------------------------------------------------------------------------------------
  !> intialize all variables (TYPE t_daily_average)
  !! 
  !! just set to zero or FALSE, called in the driver
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE init_moving_avrg_variables (this)
    IMPLICIT NONE 
    CLASS(t_moving_average), INTENT(inout) :: this        !< calculates and stores all moving average values

    !< 1.0 init/set reasonable default values
    !!
    ! 1.1 light in the canopy
    this%par_sunlit                    = 1.0_wp
    this%par_shaded                    = 1.0_wp
    this%fleaf_sunlit                  = 1.0_wp
    ! 1.2 soil moisture
    this%beta_air                      = 1.0_wp
    this%beta_soil_ps                  = 1.0_wp
    this%beta_soil_gs                  = 1.0_wp
    ! 1.3 memory of fluxes
    this%gpp_week                      = 1.0_wp
    this%maint_respiration_week        = 1.0_wp
    ! 1.4 nutrient demand for plant uptake
    this%demand_uptake_n               = 0.1_wp
    this%demand_uptake_p               = 0.1_wp
    this%growth_req_n_week             = 1._wp / 25._wp 
    this%growth_req_p_week             = 1._wp / 14._wp
    ! 1.5 environmental forcing
    this%t_air                         = 283.15_wp   ! 10°C; temperature in K
    this%q_air                         = 283.15_wp   ! 10°C; temperature in K
    this%press_surf                    = 1._wp
    ! 1.6 memory of certain variables, temperature [K], ...
    this%t_air_week                    = 283.15_wp   ! 10°C; 7  days
    this%t_air_month                   = 283.15_wp   ! 10°C; 30 days
    this%t_air_longterm                = 283.15_wp   ! 10°C; 30 years == 10950 days

    !< 2.0 init/set accumulation period per variable in days !!
    !!
    ! 2.0a some default values
    this%p_weekly                      = 7._wp     ! 7 days; constant across the model
    this%p_monthly                     = 30._wp    ! 30 days; constant across the model
    this%p_longterm                    = 10950._wp ! 10950 days (30 years); constant across the model
    ! 2.1 light in the canopy (for all canopy_layer the same period)
    this%p_par_sunlit                  = 1.0_wp 
    this%p_par_shaded                  = 1.0_wp 
    this%p_fleaf_sunlit                = 1.0_wp  
    ! 2.2 soil moisture
    this%p_beta_air                    = this%p_weekly  
    this%p_beta_soil_ps                = this%p_weekly  
    this%p_beta_soil_gs                = this%p_weekly  

    ! 2.4 nutrient demand for plant uptake (average demand over a couple of days to get rid of diurnal C cycle)
    this%p_demand_uptake_n             = 3.0_wp 
    this%p_demand_uptake_p             = 3.0_wp 

    ! 2.5 environmental forcing
    this%p_t_air                       = 1.0_wp
    this%p_q_air                       = 1.0_wp
    this%p_press_surf                  = 1.0_wp
  END SUBROUTINE init_moving_avrg_variables

  !-----------------------------------------------------------------------------------------------------
  !> calculate and return the moving average of a single variable or an array
  !!
  !! calculates the average across a certain period preceding to the current value; 
  !! actually, the average gets rather estimated than exactly calculated using the weighted arithmetic mean 
  !-----------------------------------------------------------------------------------------------------
  PURE ELEMENTAL FUNCTION calculate_moving_average_of_variable(this, cav,nv,w,mask) RESULT(nav)

    IMPLICIT NONE
    ! 0.1 InOut
    CLASS(t_moving_average), INTENT(in) :: this  !< calculates and stores all moving average values
    REAL(wp), INTENT(in)                :: cav   !< current average value
    REAL(wp), INTENT(in)                :: nv    !< new value
    REAL(wp), INTENT(in)                :: w     !< weight of the new value
    LOGICAL, OPTIONAL, INTENT(in)       :: mask  !< bool, whether the new average is calculated from the nv, or av is used as nav
    REAL(wp)                            :: nav   !< function output: the new average
    ! 0.2 Local
    LOGICAL                             :: l_mask ! True/ False of the Optional variable mask present

    ! check whether the optional statement is provided
    l_mask = .FALSE.
    IF(PRESENT(mask))  l_mask = mask

    ! calc average as weighted mean
    IF(l_mask) THEN
      nav = cav
    ELSE
      nav = cav * (1.0_wp - w) + nv * w
    END IF
  END FUNCTION calculate_moving_average_of_variable

  !-----------------------------------------------------------------------------------------------------
  !> calculate moving average of vegetation variables
  !!
  !! called in the driver at the end of each timestep
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE update_veg_moving_avrg_variables (this, plant)

    USE mo_base_constants,    ONLY: dt, one_day

    IMPLICIT NONE
    ! 0.1 InOut
    CLASS(t_moving_average), INTENT(inout) :: this        !< calculates and stores all moving average values
    TYPE(t_slm_plant), INTENT(inout)       :: plant       !< plant variables and data

    ! moving_avrg( 
    !             current average, new value, weight of new value, 
    !                 use function TRUE FALSE [OPTIONAL argument] )
    ! returns: new current average

    this%par_sunlit(:)    = this%calc_moving_avrg(this%par_sunlit(:), plant%daily_average%avrg_prev_par_sunlit_day(:), &
                                               1.0_wp / (this%p_par_sunlit * one_day/ dt))
    this%par_shaded(:)    = this%calc_moving_avrg(this%par_shaded(:), plant%daily_average%avrg_prev_par_shaded_day(:), &
                                               1.0_wp / (this%p_par_shaded * one_day/ dt))
    this%fleaf_sunlit(:)  = this%calc_moving_avrg(this%fleaf_sunlit(:), plant%daily_average%avrg_prev_fleaf_sunlit_day(:), &
                                               1.0_wp / (this%p_fleaf_sunlit * one_day/ dt))

    this%beta_soil_gs     = this%calc_moving_avrg(this%beta_soil_gs, plant%beta_soil_gs, &
                                               1._wp / ( this%p_beta_soil_gs * one_day/dt))

    this%gpp_week                = this%calc_moving_avrg(this%gpp_week, plant%flux%gpp, &
                                               1.0_wp / (this%p_weekly * one_day/ dt))
    this%maint_respiration_week  = this%calc_moving_avrg(this%maint_respiration_week, plant%flux%maint_respiration_pot, &
                                               1.0_wp / (this%p_weekly * one_day/ dt))
    this%growth_req_n_week       = this%calc_moving_avrg(this%growth_req_n_week, plant%growth_req_n, &
                                               1.0_wp / (this%p_weekly * one_day/ dt), &
                                               plant%growth_req_n < eps8)
    this%growth_req_p_week       = this%calc_moving_avrg(this%growth_req_p_week, plant%growth_req_p, &
                                               1.0_wp / (this%p_weekly * one_day/ dt), & 
                                               plant%growth_req_p < eps8)

  END SUBROUTINE update_veg_moving_avrg_variables

  !-----------------------------------------------------------------------------------------------------
  !> calculate moving average of soil variables
  !!
  !! called in the driver at the end of each timestep
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE update_soil_moving_avrg_variables (this, plant)

    USE mo_base_constants,    ONLY: dt, one_day

    IMPLICIT NONE
    ! 0.1 InOut
    CLASS(t_moving_average), INTENT(inout) :: this        !< calculates and stores all moving average values
    TYPE(t_slm_plant), INTENT(inout)       :: plant       !< plant variables and data

    ! moving_avrg( 
    !             current average, new value, weight of new value, 
    !                 use function TRUE FALSE [OPTIONAL argument] )
    ! returns: new current average

    this%demand_uptake_n = this%calc_moving_avrg(this%demand_uptake_n, plant%f_n_demand, &
                                                 1.0_wp / (this%p_demand_uptake_n * one_day/ dt), &
                                                 plant%pool%fine_root%carbon < eps8)
    this%demand_uptake_p = this%calc_moving_avrg(this%demand_uptake_p, plant%f_p_demand, &
                                                 1.0_wp / (this%p_demand_uptake_p * one_day/ dt), &
                                                 plant%pool%fine_root%carbon < eps8)
  END SUBROUTINE update_soil_moving_avrg_variables

  !-----------------------------------------------------------------------------------------------------
  !> calculate moving average of tile variables
  !!
  !! called in the driver at the end of each timestep
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE update_tile_moving_avrg_variables (this, t_air, press_surf)

    USE mo_base_constants,    ONLY: dt, one_day

    IMPLICIT NONE
    ! 0.1 InOut
    REAL(wp), INTENT(in)                   :: t_air       !< current air temperature from atm2land:mo_atmland_memory
    REAL(wp), INTENT(in)                   :: press_surf  !< variable from atm2land:mo_atmland_memory
    CLASS(t_moving_average), INTENT(inout) :: this        !< calculates and stores all moving average values

    ! moving_avrg( 
    !             current average, new value, weight of new value, 
    !                 use function TRUE FALSE [OPTIONAL argument] )
    ! returns: new current average

    this%t_air            = this%calc_moving_avrg(this%t_air, t_air, &
                                                  1.0_wp / (this%p_t_air * one_day/ dt))
    this%press_surf       = this%calc_moving_avrg(this%press_surf, press_surf, &
                                                  1.0_wp / (this%p_press_surf * one_day/ dt))

    this%t_air_week       = this%calc_moving_avrg(this%t_air_week, t_air, &
                                                  1.0_wp / (this%p_weekly * one_day/ dt))
    this%t_air_month      = this%calc_moving_avrg(this%t_air_month, t_air, &
                                                  1.0_wp / (this%p_monthly * one_day/ dt))
    this%t_air_longterm   = this%calc_moving_avrg(this%t_air_longterm, t_air, &
                                                  1.0_wp / (this%p_longterm * one_day/ dt))
  END SUBROUTINE update_tile_moving_avrg_variables


END MODULE mo_veg_memory
