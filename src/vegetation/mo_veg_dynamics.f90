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

  USE mo_kind,                 ONLY: wp, eps4, eps8
  USE mo_jsb_control,          ONLY: debug_on
  USE mo_exception,            ONLY: message

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: update_veg_dynamics

  CHARACTER(len=*), PARAMETER :: modname = 'mo_veg_dynamics'

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
  SUBROUTINE update_veg_dynamics(tile, options)

    USE mo_jsb_lctlib_class,      ONLY: t_lctlib_element
    USE mo_jsb_class,             ONLY: Get_model
    USE mo_base_memory,           ONLY: t_slm_tile, t_jsb_task_options
    USE mo_base_config_class,     ONLY: base_config
    USE mo_jsb_model_class,       ONLY: t_jsb_model
    USE mo_jsb_process_class,     ONLY: A2L_, L2A_, SEB_, TURB_, SOIL_, HYDRO_, VEG_, HD_, RAD_, ASSIMI_, PHENO_, SB_

    ! Use of process configurations (t_PROC_config)
    dsl4jsb_Use_config(ASSIMI_)
    dsl4jsb_Use_config(VEG_)

    ! Use of process memories
    dsl4jsb_Use_memory(VEG_)        !  USE mo_veg_memory_class,      ONLY: t_veg_memory
    dsl4slm_Use_memory_a_veg        !  USE mo_veg_memory_class,      ONLY: t_veg_element_a, t_veg_pool_a
    USE mo_veg_memory_class,      ONLY: veg_pool, veg_litterfall, veg_establishment
    dsl4jsb_Use_memory(RAD_)        !  USE mo_rad_memory_class,      ONLY: t_rad_memory
    dsl4jsb_Use_memory(PHENO_)      !  USE mo_pheno_memory_class,    ONLY: t_pheno_memory
    dsl4jsb_Use_memory(SOIL_)       !  USE mo_soil_memory_class,     ONLY: t_soil_memory

    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    TYPE(t_slm_tile),         INTENT(inout) :: tile             !< one tile
    TYPE(t_jsb_task_options), INTENT(in)    :: options          !< model options
    ! ---------------------------
    ! 0.2 Local 
    TYPE(t_jsb_model),      POINTER       :: model
    TYPE(t_lctlib_element), POINTER       :: lctlib     !< land-cover-type library - parameter across pft's
    INTEGER                               :: iblk, ics, ice, nc
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':update_veg_dynamics'
    ! ---------------------------
    ! 0.3 Declare Memory
    ! Declare process configuration and memory Pointers
    dsl4jsb_Def_config(ASSIMI_)
    dsl4jsb_Def_config(VEG_)
    dsl4jsb_Def_memory(VEG_)                              ! TYPE(t_veg_memory),      POINTER :: veg__mem
    dsl4jsb_Def_memory(RAD_)                              ! TYPE(t_rad_memory),      POINTER :: rad__mem
    dsl4jsb_Def_memory(PHENO_)                            ! TYPE(t_pheno_memory),    POINTER :: pheno__mem
    dsl4jsb_Def_memory(SOIL_)                             ! TYPE(t_soil_memory),     POINTER :: soil__mem
    ! Declare pointers to variables in memory
    ! RAD_
    dsl4jsb_Real2D_onChunk      :: rfr_ratio_boc_tvegdyn_mavg
    ! PHENO_
    dsl4jsb_Logical2D_onChunk   :: growing_season
    ! SOIL_
    dsl4jsb_Real2D_onChunk      :: w_soil_root_theta
    ! VEG_
    dsl4jsb_Real2D_onChunk      :: dens_ind
    dsl4jsb_Real2D_onChunk      :: mortality_rate
    dsl4jsb_Real2D_onChunk      :: delta_dens_ind
    dsl4jsb_Real2D_onChunk      :: t_air_week_mavg
    dsl4jsb_Real2D_onChunk      :: an_boc_tvegdyn_mavg
    ! Get local variables from options argument
    iblk    = options%iblk
    ics     = options%ics
    ice     = options%ice
    nc      = options%nc
    ! ---------------------------
    ! 0.4 Process Activity, Debug Option
    IF (.NOT. tile%Is_process_active(VEG_)) RETURN
    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')
    ! ---------------------------
    ! 0.5 Get Memory
    model  => Get_model(tile%owner_model_id)
    lctlib => model%lctlib(tile%lcts%lib_id)
    ! Get process config
    dsl4jsb_Get_config(ASSIMI_)
    dsl4jsb_Get_config(VEG_)
    ! Get process memories
    dsl4jsb_Get_memory(VEG_)                              ! veg__mem => veg_memory
    dsl4jsb_Get_memory(RAD_)                              ! rad__mem => rad_memory
    dsl4jsb_Get_memory(PHENO_)                            ! pheno__mem => pheno_memory
    dsl4jsb_Get_memory(SOIL_)                             ! soil__mem => soil_memory
    ! Get process variables (Set pointers to variables in memory)
    dsl4jsb_Get_var2D_onChunk(RAD_, rfr_ratio_boc_tvegdyn_mavg) ! in
    ! ---------------------------
    dsl4jsb_Get_var2D_onChunk(PHENO_, growing_season)       ! in
    ! ---------------------------
    dsl4jsb_Get_var2D_onChunk(SOIL_, w_soil_root_theta)     ! in
    ! ---------------------------
    dsl4jsb_Get_var2D_onChunk(VEG_, dens_ind)               ! in
    dsl4jsb_Get_var2D_onChunk(VEG_, mortality_rate)         ! out
    dsl4jsb_Get_var2D_onChunk(VEG_, delta_dens_ind)         ! out
    dsl4jsb_Get_var2D_onChunk(VEG_, t_air_week_mavg)        ! in
    dsl4jsb_Get_var2D_onChunk(VEG_, an_boc_tvegdyn_mavg)    ! in
    ! Get process memory pool structures
#define Get_veg_pool
#define Get_veg_litterfall
#define Get_veg_establishment
#include "slm_pools_a_incl.f90"
#undef Get_veg_pool
#undef Get_veg_litterfall
#undef Get_veg_establishment



    ! ------------------------------------------------------------------------------------------------------------
    ! Go science
    ! ------------------------------------------------------------------------------------------------------------


    !> 1.0 Calculate mortality rate and veg\%litterfall
    !! 
    !! retruns constant mortality_rate if flag_tile_veg_dynamics == false \n
    !! but calcs mortality_rate dynamically if flag_tile_veg_dynamics == true
    CALL calc_veg_mortality(options                                     , & ! in
                            lctlib                                      , & 
                            dsl4jsb_Config(VEG_)%flag_tile_veg_dynamics , &
                            an_boc_tvegdyn_mavg(:)                      , &
                            veg_pool                                    , &
                            veg_litterfall                              , & ! inout
                            mortality_rate(:)                           )   ! out

    !> 1.1 Calculate stand harvest_rate and veg\%litterfall
    !! 
    IF(base_config%l_do_stand_replacing_harvest)THEN 
       CALL calc_veg_harvest(options                                     , & ! in
                             veg_pool                                    , &
                             veg_litterfall                              , & ! inout
                             mortality_rate(:)                           )   ! out
    ENDIF


    !> 2.0 with vegetation dynamics turned on (flag_tile_veg_dynamics == true)
    !! 
    IF(dsl4jsb_Config(VEG_)%flag_tile_veg_dynamics.AND..NOT.base_config%l_do_stand_replacing_harvest) THEN 

      !> 2.1 Calculate establishment_rate and veg\%establishment
      !! 
      CALL calc_veg_establishment(options                         , & ! in
                                  lctlib                          , & 
                                  dsl4jsb_Config(ASSIMI_)%ncanopy , &
                                  w_soil_root_theta(:)            , &
                                  growing_season(:)               , &
                                  t_air_week_mavg(:)              , &
                                  rfr_ratio_boc_tvegdyn_mavg(:)   , &
                                  veg_pool%seed_bed               , &
                                  veg_establishment               )   ! out

      !> 2.2 Calculate change in individual density from mortality and establishment
      !! 
      CALL calc_delta_density_individuals(lctlib                        , & ! in
                                          dens_ind(:)                   , &
                                          mortality_rate(:)             , &
                                          veg_establishment%carbon(:)   , &
                                          delta_dens_ind(:)             )   ! out

    END IF  ! dsl4jsb_Config(VEG_)%flag_tile_veg_dynamics

! Pass back values to process memory
#define Passback_veg_pool
#define Passback_veg_litterfall
#define Passback_veg_establishment
#include "slm_pools_a_incl.f90"
#undef Passback_veg_pool
#undef Passback_veg_litterfall
#undef Passback_veg_establishment

  END SUBROUTINE update_veg_dynamics


  !-----------------------------------------------------------------------------------------------------
  ! Sub Task of update_veg_dynamics
  ! 
  !-----------------------------------------------------------------------------------------------------
  !> calculate vegetation mortality
  !!
  !! veg\%litterfall
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE calc_veg_mortality(options                   , &
                                lctlib                    , &
                                flag_tile_veg_dynamics    , &
                                an_boc_tvegdyn_mavg       , &
                                veg_pool                  , &
                                veg_litterfall            , &
                                mortality_rate )

    USE mo_base_constants
    USE mo_base_memory,           ONLY: t_jsb_task_options
    USE mo_veg_constants
    USE mo_jsb_lctlib_class,      ONLY: t_lctlib_element
    dsl4slm_Use_memory_a_veg

    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    TYPE(t_jsb_task_options),     INTENT(in)    :: options                  !< model options
    TYPE(t_lctlib_element),       INTENT(in)    :: lctlib                   !< land-cover-type library - parameter across pft's
    LOGICAL,                      INTENT(in)    :: flag_tile_veg_dynamics
    REAL(wp), DIMENSION(:),       INTENT(in)    :: an_boc_tvegdyn_mavg      !< moving average ...
    TYPE(t_veg_pool_a),           INTENT(in)    :: veg_pool                 !< vegetation pool
    TYPE(t_veg_pool_a),           INTENT(inout) :: veg_litterfall           !< vegetation litterfall
    REAL(wp), DIMENSION(:),       INTENT(out)   :: mortality_rate           !< background_mort_rate_tree/grass converted to 1/timestep
    ! ---------------------------
    ! 0.2 Local
    REAL(wp), DIMENSION(options%nc)              :: boc_cbalance
    REAL(wp), DIMENSION(options%nc)              :: fscal 
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_veg_mortality'



    !> 1.0 background mortality 
    !! 
    !! select rate and convert from 1/year to 1/timestep
    SELECT CASE (lctlib%growthform)
      CASE (itree)
         mortality_rate(:) = ((background_mort_rate_tree /  one_year) / one_day) * dtime
      CASE (igrass)
         mortality_rate(:) = ((background_mort_rate_grass / one_year) / one_day) * dtime
    ENDSELECT 

    !> 2.0 self-thinning related mortality
    !! 
    !! depends on the carbon balance of lowest canopy layer, i.e. net photosynthesis minus average construction cost \n
    !! in mol / m2 LAI / year \n
    !! transformed to a mortality estimate using a Weibull function
    IF(flag_tile_veg_dynamics) THEN 
       boc_cbalance(:) = an_boc_tvegdyn_mavg(:) * one_day * one_year / 1.e6_wp &  
                         - ( 1._wp + fresp_growth ) * 1.0_wp/lctlib%tau_leaf/lctlib%sla
       WHERE(boc_cbalance(:) < -eps4 )
          mortality_rate(:) = mortality_rate(:) & 
                              + ( 1._wp - exp ( - ( - lambda_mort_light * boc_cbalance(:) ) ** k_mort_light)) & 
                              / one_year / one_day * dtime
       END WHERE
    ENDIF

    !> 3.0 calculate litter fall from mortality
    !!
    !! for all tissues but seed_bed
    WHERE(veg_pool%leaf%carbon(:) > eps8)
       fscal(:) = (veg_pool%leaf%carbon(:) - veg_litterfall%leaf%carbon(:)) / veg_pool%leaf%carbon(:)
    ELSEWHERE
       fscal(:) = 0.0_wp 
    ENDWHERE
    WHERE(mortality_rate(:)<fscal(:))
       fscal(:) = mortality_rate
    ENDWHERE 
    veg_litterfall%leaf         = veg_litterfall%leaf         + veg_pool%leaf         * fscal(:)
    veg_litterfall%fine_root    = veg_litterfall%fine_root    + veg_pool%fine_root    * mortality_rate(:)
    veg_litterfall%coarse_root  = veg_litterfall%coarse_root  + veg_pool%coarse_root  * mortality_rate(:)
    veg_litterfall%sap_wood     = veg_litterfall%sap_wood     + veg_pool%sap_wood     * mortality_rate(:)
    veg_litterfall%heart_wood   = veg_litterfall%heart_wood   + veg_pool%heart_wood   * mortality_rate(:)
    veg_litterfall%labile       = veg_litterfall%labile       + veg_pool%labile       * mortality_rate(:)
    veg_litterfall%reserve      = veg_litterfall%reserve      + veg_pool%reserve      * mortality_rate(:)
    veg_litterfall%fruit        = veg_litterfall%fruit        + veg_pool%fruit        * mortality_rate(:)

  END SUBROUTINE calc_veg_mortality

  !-----------------------------------------------------------------------------------------------------
  ! Sub Task of update_veg_dynamics
  ! 
  !-----------------------------------------------------------------------------------------------------
  !> calculate vegetation harvest 
  !!
  !! veg\%litterfall
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE calc_veg_harvest(options                   , &
                              veg_pool                  , &
                              veg_litterfall            , &
                              mortality_rate)

    USE mo_base_constants
    USE mo_base_memory,           ONLY: t_jsb_task_options
    USE mo_veg_constants,         ONLY: wood_extraction_rate
    dsl4slm_Use_memory_a_veg

    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    TYPE(t_jsb_task_options),     INTENT(in)    :: options                  !< model options
    TYPE(t_veg_pool_a),           INTENT(in)    :: veg_pool                 !< vegetation pool
    TYPE(t_veg_pool_a),           INTENT(inout) :: veg_litterfall           !< vegetation litterfall
    REAL(wp), DIMENSION(:),       INTENT(inout) :: mortality_rate           !< background_mort_rate_tree/grass converted to 1/timestep
    ! ---------------------------
    ! 0.2 Local
    REAL(wp), DIMENSION(options%nc)             :: rest_mortality_rate 
    REAL(wp)                                    :: export_rate
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_veg_harvest'

    !> 1.0 calculate rest mortality to achieve stand-replacing harvest
    !!
    rest_mortality_rate(:) = 1._wp - mortality_rate(:)

    veg_litterfall%leaf         = veg_litterfall%leaf         + veg_pool%leaf         * rest_mortality_rate(:)
    veg_litterfall%fine_root    = veg_litterfall%fine_root    + veg_pool%fine_root    * rest_mortality_rate(:)
    veg_litterfall%coarse_root  = veg_litterfall%coarse_root  + veg_pool%coarse_root  * rest_mortality_rate(:)
    veg_litterfall%sap_wood     = veg_litterfall%sap_wood     + veg_pool%sap_wood     * (1._wp - wood_extraction_rate) * rest_mortality_rate(:)
    veg_litterfall%heart_wood   = veg_litterfall%heart_wood   + veg_pool%heart_wood   * (1._wp - wood_extraction_rate) * rest_mortality_rate(:)
    veg_litterfall%labile       = veg_litterfall%labile       + veg_pool%labile       * (1._wp - wood_extraction_rate) * rest_mortality_rate(:)
    veg_litterfall%reserve      = veg_litterfall%reserve      + veg_pool%reserve      * (1._wp - wood_extraction_rate) * rest_mortality_rate(:)
    veg_litterfall%fruit        = veg_litterfall%fruit        + veg_pool%fruit        * rest_mortality_rate(:)

    mortality_rate(:) = 1.0_wp

  END SUBROUTINE calc_veg_harvest

  !-----------------------------------------------------------------------------------------------------
  ! Sub Task of update_veg_dynamics
  ! 
  !-----------------------------------------------------------------------------------------------------
  !> calculate vegetation establishment 
  !!
  !!
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE calc_veg_establishment(options             , &
                                    lctlib              , &
                                    ncanopy             , &
                                    w_soil_root_theta   , &
                                    growing_season      , &
                                    t_air_week_mavg     , &
                                    rfr_ratio_boc_tvegdyn_mavg , &
                                    veg_pool_seed_bed   , &
                                    veg_establishment )


    USE mo_base_constants,        ONLY: dtime, one_day, one_year
    USE mo_base_memory,           ONLY: t_jsb_task_options
    USE mo_veg_constants
    USE mo_jsb_lctlib_class,      ONLY: t_lctlib_element
    USE mo_shared_constants,      ONLY: Tzero
    dsl4slm_Use_memory_a_veg

    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    TYPE(t_jsb_task_options),       INTENT(in)    :: options            !< model options
    TYPE(t_lctlib_element),         INTENT(in)    :: lctlib             !< land-cover-type library - parameter across pft's
    INTEGER,                        INTENT(in)    :: ncanopy            !< number of canopy layer
    REAL(wp), DIMENSION(:),         INTENT(in)    :: w_soil_root_theta  !< plant available water [fraction of maximum] 
    LOGICAL , DIMENSION(:),         INTENT(in)    :: growing_season     !< growing season 
    REAL(wp), DIMENSION(:),         INTENT(in)    :: t_air_week_mavg    !< weekly air temperature [K] 
    REAL(wp), DIMENSION(:),         INTENT(in)    :: rfr_ratio_boc_tvegdyn_mavg !< red-farred ratio at the bottom of the canopy 
    TYPE(t_veg_element_a),          INTENT(in)    :: veg_pool_seed_bed  !< seed_bed in vegetation pool
    TYPE(t_veg_element_a),          INTENT(inout) :: veg_establishment  !< vegetation establishment
    ! ---------------------------
    ! 0.2 Local
    REAL(wp), DIMENSION(options%nc)               :: establishment_rate , &
                                                     tc                 , &
                                                     rfr_ratio          , &
                                                     flim_light         , &
                                                     flim_moist         , &
                                                     flim_temp
    INTEGER                                       :: icanopy
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_veg_establishment'


    !> 1.0 establishment rate
    !!

    !> 1.1 light limitation
    !! depending on the red to far-red ratio of the light at the bottom of the canopy
    flim_light(:) = 1._wp-exp(-(lctlib%lambda_est_light*rfr_ratio_boc_tvegdyn_mavg(:))**lctlib%k_est_light)

    !> 1.2 temperature and soil moisture limitation 
    !! 
    !! OBS: should be top soil moisture!
    WHERE(t_air_week_mavg(:) - Tzero < 0.0_wp .OR..NOT. growing_season(:))
      flim_temp(:) = 0.0_wp
    ELSE WHERE
      flim_temp(:) = 1.0_wp - exp(-(lambda_est_temp * (t_air_week_mavg(:) - Tzero))**k_est_temp)
    END WHERE
    flim_moist(:)  = 1.0_wp - exp(-(lambda_est_moist * w_soil_root_theta(:))**k_est_moist) 
   
    !> 1.3 actual establishment rate and associated matter flux from seed bed to reserve pool
    !!
    establishment_rate(:) = flim_light(:) * flim_temp(:) * flim_moist(:) * & 
                            1.0_wp/lctlib%tau_seed_est/one_day/one_year * dtime
    veg_establishment     = veg_pool_seed_bed * establishment_rate(:)

  END SUBROUTINE calc_veg_establishment


  !-----------------------------------------------------------------------------------------------------
  ! Sub Task of update_veg_dynamics
  ! 
  !-----------------------------------------------------------------------------------------------------
  !> calculate delta of individual density
  !!
  !!
  !-----------------------------------------------------------------------------------------------------
  ELEMENTAL SUBROUTINE calc_delta_density_individuals(lctlib          , &
                                            dens_ind                  , &
                                            mortality_rate            , &
                                            veg_establishment_carbon  , &
                                            delta_dens_ind )

    USE mo_jsb_lctlib_class,        ONLY: t_lctlib_element
    USE mo_veg_constants,           ONLY: itree

    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    TYPE(t_lctlib_element),         INTENT(in)    :: lctlib                       !< land-cover-type library - parameter across pft's
    REAL(wp),                       INTENT(in)    :: dens_ind                 , & !< soil moisture 
                                                     mortality_rate           , & !< mortality per timestep
                                                     veg_establishment_carbon     !< egetation establishment measured in C
    REAL(wp),                       INTENT(out)   :: delta_dens_ind               !< diff in individual density
    ! ---------------------------
    ! 0.2 Local
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_delta_density_individuals'



    !> 1.0 calc delta of individual density (#/m2)
    !! 

    !> 1.1 trees
    !!
    IF(lctlib%growthform == itree) THEN
       delta_dens_ind = (veg_establishment_carbon / lctlib%seed_size) - (mortality_rate * dens_ind)
    END IF

  END SUBROUTINE calc_delta_density_individuals

END MODULE mo_veg_dynamics
