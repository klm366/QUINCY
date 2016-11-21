!>
!! @brief module to calculate vegetation dynamics 
!!
!! @author
!!  Sönke Zaehle, MPI-BGC Jena
!!
!! @par Revision History
!! First version by Sönke Zaehle (2016-04-14)
!!
MODULE mo_veg_dynamics

#include "dsl_jsbach.inc"

  USE mo_kind,                 ONLY: wp, eps4


  IMPLICIT NONE

  PRIVATE
  PUBLIC :: update_veg_dynamics

CONTAINS

  !-----------------------------------------------------------------------------------------------------
  ! Main Task
  ! 
  !-----------------------------------------------------------------------------------------------------
  !> Calculate establishment and mortality of a population (currently static mortality rate only) 
  !!
  !! a background_mort_rate exists for each trees and grasses, defined as mortality rate per year
  !!
  !! fluxes are set to zero at each time step: set_zero_vegetation_fluxes:mo_assimi_update_canopy_fluxes
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE update_veg_dynamics(tile)

    USE mo_jsb_lctlib_class,      ONLY: t_lctlib_element
    USE mo_jsb_class,             ONLY: Get_model
    USE mo_base_memory,           ONLY: t_slm_tile, t_model
    USE mo_jsb_process_class,     ONLY: ATM2LAND_, LAND2ATM_, SEB_, TURB_, SOIL_, HYDRO_, VEG_, HD_, RAD_, ASSIMI_, PHENO_

    ! Use of process memories
    dsl4jsb_Use_memory(VEG_)        !  USE mo_veg_memory_class,      ONLY: t_veg_memory


    IMPLICIT NONE
    ! 0.1 InOut
    TYPE(t_slm_tile),       INTENT(inout)     :: tile           !< current tile 
    ! 0.2 Local
    TYPE(t_model),          POINTER           :: model
    TYPE(t_lctlib_element), POINTER           :: lctlib     !< land-cover-type library - parameter across pft's
    ! 0.3 Memory allocation
    ! Declare process configuration and memory Pointers
    dsl4jsb_Def_memory(VEG_)                          ! TYPE(t_veg_memory),      POINTER :: veg__mem
    TYPE(t_veg_memory), POINTER         :: veg        ! local pointer to veg__mem for SLM usability
    dsl4jsb_Get_memory(VEG_)                          ! veg__mem => veg_memory

    veg => veg__mem 
    model  => Get_model(tile%owner_model_id)
    lctlib => model%lctlib(tile%lcts%lib_id)


    !> Sub-Task 1
    !! 
    CALL calc_veg_mortality(lctlib                      , &
                            veg%an_boc_tvegdyn_mavg%ptr , &
                            veg%pool                    , &
                            veg%litterfall              , &
                            veg%mortality_rate%ptr )
    !> Sub-Task 2
    !! 
    CALL calc_veg_establishment(lctlib                        , &
                                 tile%soil%theta              , &
                                 veg%growing_season%ptr       , &
                                 veg%t_air_week_mavg%ptr      , &
                                 veg%nl_tot_cl%ptr(:,1)       , &
                                 veg%fnl_chl_cl%ptr(:,1)      , &
                                 veg%lai_cl%ptr(:,1)          , &
                                 veg%pool%seed_bed            , &
                                 veg%establishment )
    !> Sub-Task 3
    !! 
    CALL calc_delta_density_individuals(lctlib%seed_size              , &
                                        lctlib%growthform             , &
                                        veg%dens_ind%ptr              , &
                                        veg%mortality_rate%ptr        , &
                                        veg%establishment%carbon%ptr  , &
                                        veg%delta_dens_ind%ptr )


  END SUBROUTINE update_veg_dynamics


  !-----------------------------------------------------------------------------------------------------
  ! Sub Task of update_veg_dynamics
  ! 
  !-----------------------------------------------------------------------------------------------------
  !> calculate vegetation mortality
  !!
  !! veg\%litterfall
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE calc_veg_mortality(lctlib                    , &
                                an_boc_tvegdyn_mavg       , &
                                veg_pool                  , &
                                veg_litterfall            , &
                                mortality_rate )

    USE mo_base_constants
    USE mo_veg_constants
    USE mo_jsb_lctlib_class,      ONLY: t_lctlib_element
    USE mo_veg_config_class,      ONLY: veg_config
    ! veg memory
    USE mo_veg_memory_class,      ONLY: t_slm_struct_veg_pool

    IMPLICIT NONE

    ! 0.1 InOut
    TYPE(t_lctlib_element),       INTENT(in)    :: lctlib                   !< land-cover-type library - parameter across pft's
    REAL(wp),                     INTENT(in)    :: an_boc_tvegdyn_mavg      !< moving average ...
    TYPE(t_slm_struct_veg_pool),  INTENT(in)    :: veg_pool                 !< vegetation pool
    TYPE(t_slm_struct_veg_pool),  INTENT(inout) :: veg_litterfall           !< vegetation litterfall
    REAL(wp),                     INTENT(out)   :: mortality_rate           !< background_mort_rate_tree/grass converted to 1/timestep
    ! 0.2 Local
    REAL(wp)                                    :: boc_cbalance


    !> 1.0 background mortality 
    !! 
    !! select rate and convert from 1/year to 1/timestep
    SELECT CASE (lctlib%growthform)
      CASE (itree)
         mortality_rate = ((background_mort_rate_tree /  one_year) / one_day) * dt
      CASE (igrass)
         mortality_rate = ((background_mort_rate_grass / one_year) / one_day) * dt
    ENDSELECT 

    IF(veg_config%flag_tile_veg_dynamics) THEN 

       !> 2.0 self-thinning related mortality
       !! 
       !! depends on the carbon balance of lowest canopy layer, i.e. net photosynthesis minus average construction cost \n
       !! in mol / m2 LAI / year \n
       !! transformed to a mortality estimate using a Weibull function
       boc_cbalance = an_boc_tvegdyn_mavg * one_day * one_year / 1.e6_wp &  
                    - ( 1._wp + fresp_growth ) * lctlib%tau_leaf/lctlib%sla
       IF(boc_cbalance < -eps4 )THEN
          mortality_rate = mortality_rate & 
                           + ( 1._wp - exp ( - ( - lambda_mort_light * boc_cbalance ) ** k_mort_light)) & 
                           / one_year / one_day * dt
       ENDIF

    ENDIF

    !> 4.0 calculate litter fall from mortality
    !!
    !! for all tissues but seed_bed
    veg_litterfall%leaf         = veg_litterfall%leaf         + veg_pool%leaf         * mortality_rate
    veg_litterfall%fine_root    = veg_litterfall%fine_root    + veg_pool%fine_root    * mortality_rate
    veg_litterfall%coarse_root  = veg_litterfall%coarse_root  + veg_pool%coarse_root  * mortality_rate
    veg_litterfall%sap_wood     = veg_litterfall%sap_wood     + veg_pool%sap_wood     * mortality_rate
    veg_litterfall%heart_wood   = veg_litterfall%heart_wood   + veg_pool%heart_wood   * mortality_rate
    veg_litterfall%labile       = veg_litterfall%labile       + veg_pool%labile       * mortality_rate
    veg_litterfall%reserve      = veg_litterfall%reserve      + veg_pool%reserve      * mortality_rate
    veg_litterfall%fruit        = veg_litterfall%fruit        + veg_pool%fruit        * mortality_rate

  END SUBROUTINE calc_veg_mortality


  !-----------------------------------------------------------------------------------------------------
  ! Sub Task of update_veg_dynamics
  ! 
  !-----------------------------------------------------------------------------------------------------
  !> calculate vegetation establishment 
  !!
  !!
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE calc_veg_establishment(lctlib              , &
                                    soil_theta          , &
                                    growing_season      , &
                                    t_air_week_mavg     , &
                                    nl_tot_cl           , &
                                    fnl_chl_cl          , &
                                    lai_cl              , &
                                    veg_pool_seed_bed   , &
                                    veg_establishment )


    USE mo_base_constants,        ONLY: dt, one_day, one_year, ncanop
    USE mo_veg_constants
    USE mo_jsb_lctlib_class,      ONLY: t_lctlib_element
    USE mo_veg_config_class,      ONLY: veg_config
    USE mo_shared_constants,      ONLY: Tzero
    ! veg memory
    USE mo_veg_memory_class,      ONLY: t_slm_struct_veg_element

    IMPLICIT NONE

    ! 0.1 InOut
    TYPE(t_lctlib_element),         INTENT(in)    :: lctlib             !< land-cover-type library - parameter across pft's
    REAL(wp),                       INTENT(in)    :: soil_theta         !< soil moisture 
    LOGICAL,                        INTENT(in)    :: growing_season     !< growing season Yes/No
    REAL(wp),                       INTENT(in)    :: t_air_week_mavg    !< moving average ...
    REAL(wp), DIMENSION(:),         INTENT(in)    :: nl_tot_cl      , & !< ..
                                                     fnl_chl_cl     , & !< ..
                                                     lai_cl             !< ..
    TYPE(t_slm_struct_veg_element), INTENT(in)    :: veg_pool_seed_bed    !< seed_bed in vegetation pool
    TYPE(t_slm_struct_veg_element), INTENT(out)   :: veg_establishment    !< vegetation establishment
    ! 0.2 Local
    REAL(wp)                                      :: establishment_rate
    REAL(wp)                                      :: tc
    REAL(wp)                                      :: rfr_ratio
    REAL(wp)                                      :: flim_light   , &
                                                     flim_moist   , &
                                                     flim_temp
    INTEGER                                       :: icanop


    IF(veg_config%flag_tile_veg_dynamics) THEN 

       !> 1.0 establishment rate
       !!

       !> 1.1 light limitation
       !! depending on the red to far-red ratio of the light at the bottom of the canopy
       tc = 0.0_wp
       DO icanop=1,ncanop
         tc = tc + nl_tot_cl(icanop) * fnl_chl_cl(icanop) * lai_cl(icanop)
       ENDDO
       rfr_ratio = toc_rfr_ratio * exp( k_r2fr_chl * tc) 
       flim_light = 1._wp-exp(-(lctlib%lambda_est_light*rfr_ratio)**lctlib%k_est_light)

       !> 1.2 temperature and soil moisture limitation 
       !! 
       !! OBS: should be top soil moisture!
       IF(t_air_week_mavg - Tzero < 0.0_wp .OR..NOT. growing_season) THEN 
          flim_temp = 0.0_wp
       ELSE
          flim_temp = 1.0_wp - exp(-(lambda_est_temp * (t_air_week_mavg - Tzero))**k_est_temp)
       ENDIF
       flim_moist   = 1.0_wp - exp(-(lambda_est_moist * soil_theta)**k_est_moist) 
   
       !> 1.3 actual establishment rate and associated matter flux from seed bed to reserve pool
       !!
       establishment_rate = flim_light * flim_temp * flim_moist * & 
                            lctlib%tau_seed_est / one_year / one_day * dt
       veg_establishment  = veg_pool_seed_bed * establishment_rate

    ENDIF

  END SUBROUTINE calc_veg_establishment


  !-----------------------------------------------------------------------------------------------------
  ! Sub Task of update_veg_dynamics
  ! 
  !-----------------------------------------------------------------------------------------------------
  !> calculate delta of individual density
  !!
  !!
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE calc_delta_density_individuals(lctlib_seed_size          , &
                                            lctlib_growthform         , &
                                            dens_ind                  , &
                                            mortality_rate            , &
                                            veg_establishment_carbon  , &
                                            delta_dens_ind )

    USE mo_veg_constants,           ONLY: itree

    IMPLICIT NONE

    ! 0.1 InOut
    REAL(wp),                       INTENT(in)    :: lctlib_seed_size         , & !< seed size
                                                     dens_ind                 , & !< soil moisture 
                                                     mortality_rate           , & !< mortality per timestep
                                                     veg_establishment_carbon     !< egetation establishment measured in C
    INTEGER,                        INTENT(in)    :: lctlib_growthform            !< tree / grass / .. 
    REAL(wp),                       INTENT(out)   :: delta_dens_ind               !< diff in individual density


    !> 1.0 calc delta of individual density (#/m2)
    !! 

    !> 1.1 trees
    !!
    IF(lctlib_growthform == itree) THEN
	    delta_dens_ind = (veg_establishment_carbon / lctlib_seed_size) - (mortality_rate * dens_ind)
    END IF

  END SUBROUTINE calc_delta_density_individuals

END MODULE mo_veg_dynamics
