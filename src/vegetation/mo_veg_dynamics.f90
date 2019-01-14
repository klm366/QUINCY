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

  USE mo_kind,                 ONLY: wp, eps4, eps8, false, true, test_false_true
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
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE update_veg_dynamics(tile, options)

    USE mo_jsb_lctlib_class,      ONLY: t_lctlib_element
    USE mo_jsb_class,             ONLY: Get_model
    USE mo_jsb_tile_class,        ONLY: t_jsb_tile_abstract
    USE mo_jsb_task_class,        ONLY: t_jsb_task_options
    USE mo_base_constants,        ONLY: one_year, one_day
    USE mo_jsb_model_class,       ONLY: t_jsb_model
    USE mo_jsb_process_class,     ONLY: A2L_, L2A_, SEB_, TURB_, SSE_, HYDRO_, VEG_, HD_, RAD_, ASSIMI_, PHENO_, SB_
    USE mo_jsb4_forcing,          ONLY: forcing_options
    USE mo_veg_constants,         ONLY: itree

    ! Use of process configurations (t_PROC_config)
    dsl4jsb_Use_config(ASSIMI_)
    dsl4jsb_Use_config(VEG_)

    ! Use of process memories
    dsl4jsb_Use_memory(VEG_)        !  USE mo_veg_memory_class,      ONLY: t_veg_memory
    dsl4slm_Use_memory_a_veg        !  USE mo_veg_memory_class,      ONLY: t_veg_element_a, t_veg_pool_a
    USE mo_veg_memory_class,      ONLY: veg_pool, veg_litterfall, veg_establishment
    dsl4jsb_Use_memory(RAD_)        !  USE mo_rad_memory_class,      ONLY: t_rad_memory
    dsl4jsb_Use_memory(PHENO_)      !  USE mo_pheno_memory_class,    ONLY: t_pheno_memory
    dsl4jsb_Use_memory(SSE_)       !  USE mo_sse_memory_class,     ONLY: t_sse_memory

    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    CLASS(t_jsb_tile_abstract), INTENT(inout)     :: tile         !< one tile with data structure for one lct
    TYPE(t_jsb_task_options),   INTENT(in)        :: options      !< model options
    ! ---------------------------
    ! 0.2 Local 
    TYPE(t_jsb_model),      POINTER         :: model
    TYPE(t_lctlib_element), POINTER         :: lctlib           !< land-cover-type library - parameter across pft's
    REAL(wp), DIMENSION(options%nc)         :: fpc              !< current foliage projective cover
    INTEGER                                 :: iblk, ics, ice, nc
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':update_veg_dynamics'
    ! ---------------------------
    ! 0.3 Declare Memory
    ! Declare process configuration and memory Pointers
    dsl4jsb_Def_config(ASSIMI_)
    dsl4jsb_Def_config(VEG_)
    dsl4jsb_Def_memory(VEG_)                              ! TYPE(t_veg_memory),      POINTER :: veg__mem
    dsl4jsb_Def_memory(RAD_)                              ! TYPE(t_rad_memory),      POINTER :: rad__mem
    dsl4jsb_Def_memory(PHENO_)                            ! TYPE(t_pheno_memory),    POINTER :: pheno__mem
    dsl4jsb_Def_memory(SSE_)                             ! TYPE(t_sse_memory),     POINTER :: sse__mem
    ! Declare pointers to variables in memory
    ! RAD_
    dsl4jsb_Real2D_onChunk      :: rfr_ratio_boc_tvegdyn_mavg
    ! PHENO_
    dsl4jsb_Real2D_onChunk      :: growing_season
    ! SSE_
    dsl4jsb_Real2D_onChunk      :: w_soil_root_theta
    ! VEG_
    dsl4jsb_Real2D_onChunk      :: do_cohort_harvest
    dsl4jsb_Real2D_onChunk      :: dens_ind
    dsl4jsb_Real2D_onChunk      :: diameter 
    dsl4jsb_Real2D_onChunk      :: height 
    dsl4jsb_Real2D_onChunk      :: cohort_age
    dsl4jsb_Real2D_onChunk      :: mortality_rate
    dsl4jsb_Real2D_onChunk      :: delta_dens_ind
    dsl4jsb_Real2D_onChunk      :: t_air_week_mavg
    dsl4jsb_Real2D_onChunk      :: an_boc_tvegdyn_mavg
    dsl4jsb_Real2D_onChunk      :: net_growth_tvegdyn_mavg
    dsl4jsb_Real2D_onChunk      :: lai_tvegdyn_mavg
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
    dsl4jsb_Get_memory(SSE_)                             ! sse__mem => sse_memory
    ! Get process variables (Set pointers to variables in memory)
    dsl4jsb_Get_var2D_onChunk(RAD_, rfr_ratio_boc_tvegdyn_mavg) ! in
    ! ---------------------------
    dsl4jsb_Get_var2D_onChunk(PHENO_, growing_season)       ! in
    ! ---------------------------
    dsl4jsb_Get_var2D_onChunk(SSE_, w_soil_root_theta)      ! in
    ! ---------------------------
    dsl4jsb_Get_var2D_onChunk(VEG_, do_cohort_harvest)      ! out
    dsl4jsb_Get_var2D_onChunk(VEG_, dens_ind)               ! in
    dsl4jsb_Get_var2D_onChunk(VEG_, diameter)               ! in
    dsl4jsb_Get_var2D_onChunk(VEG_, height)                 ! in
    dsl4jsb_Get_var2D_onChunk(VEG_, cohort_age)             ! inout
    dsl4jsb_Get_var2D_onChunk(VEG_, mortality_rate)         ! out
    dsl4jsb_Get_var2D_onChunk(VEG_, delta_dens_ind)         ! out
    dsl4jsb_Get_var2D_onChunk(VEG_, t_air_week_mavg)        ! in
    dsl4jsb_Get_var2D_onChunk(VEG_, an_boc_tvegdyn_mavg)    ! in
    dsl4jsb_Get_var2D_onChunk(VEG_, net_growth_tvegdyn_mavg)! in
    dsl4jsb_Get_var2D_onChunk(VEG_, lai_tvegdyn_mavg)       ! in
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


    !> 1.0 Calculate foliage projective cover if vegetation dynamics are simulated, i.e., not needed for constant mortality
    !! 
    SELECT CASE (TRIM(dsl4jsb_Config(VEG_)%veg_dynamics_scheme))
      CASE ("population","cohort")
        fpc(:) = calc_foliage_projective_cover_tile(options, lctlib, veg_pool%sap_wood%carbon(:), & 
                                                    dens_ind(:), height(:), diameter(:))
    ENDSELECT


    !> 1.1 check stand harvest year (if this functionality is enabled: model%config\%flag_stand_harvest)
    !! 
    !! but enable harvesting only after transient spinup, i.e., when l_transient_spinup=FALSE
    IF (model%config%flag_stand_harvest .AND. .NOT. model%config%l_transient_spinup) THEN
      model%config%l_do_stand_replacing_harvest = .FALSE.                                      ! per default this flag is FALSE
      ! set stand_harvest to TRUE for the 1st timestep of the stand_replacing_year only
      IF (forcing_options%forcing_year .EQ. model%config%stand_replacing_year .AND. &
          forcing_options%forcing_doy  .EQ. 1                                .AND. &
          forcing_options%forcing_hour .LT. eps4                                    ) THEN
        model%config%l_do_stand_replacing_harvest = .TRUE. 
      END IF
    END IF


    !> 1.2 Increase age of the cohort [years] & check if cohort harvest is reached
    !! 
    !! apply this only after spinup of transient mode has finished
    IF (.NOT. model%config%l_transient_spinup) THEN
      SELECT CASE (TRIM(dsl4jsb_Config(VEG_)%veg_dynamics_scheme))
        CASE ("cohort")
          cohort_age(:)        = cohort_age(:) + model%config%dtime / (one_year * one_day)
          do_cohort_harvest(:) = false                                               ! per default it is set to FALSE
          WHERE (cohort_age(:) > dsl4jsb_Config(VEG_)%cohort_harvest_interval)
            do_cohort_harvest(:)   = true
            cohort_age(:)          = 0.0_wp
          ENDWHERE
      ENDSELECT
    ENDIF


    !> 2.0 Calculate mortality rate and veg\%litterfall
    !! 
    !! retruns constant mortality_rate if veg_dynamics_scheme == none \n
    !! but calcs mortality_rate dynamically if veg_dynamics_scheme == population or cohort
    CALL calc_veg_mortality(options                                         , & ! in
                            lctlib                                          , & 
                            TRIM(dsl4jsb_Config(VEG_)%veg_dynamics_scheme)  , &
                            an_boc_tvegdyn_mavg(:)                          , &
                            net_growth_tvegdyn_mavg(:)                      , &
                            lai_tvegdyn_mavg(:)                             , &
                            fpc(:)                                          , &
                            veg_pool                                        , &
                            veg_litterfall                                  , & ! inout
                            mortality_rate(:)                               )   ! out

    !> 2.1 Calculate stand harvest_rate and veg\%litterfall
    !! 
    !! if either logical is true
    !! @todo update_veg_dynamics() do_cohort_harvest dimension should be used as nc
    IF (model%config%l_do_stand_replacing_harvest .OR. do_cohort_harvest(options%ics) > test_false_true) THEN
      CALL calc_veg_harvest(options                                     , & ! in
                            veg_pool                                    , &
                            veg_litterfall                              , & ! inout
                            mortality_rate(:)                           )   ! inout
    ENDIF


    !> 3.0 Calculate establishment_rate and veg\%establishment
    !! 
    SELECT CASE (TRIM(dsl4jsb_Config(VEG_)%veg_dynamics_scheme))
      !CASE ("cohort")
        !> no establishment in cohort mode, but after cohort harvesting \n 
        !!  the dens_ind is set to 1.0 and veg_pool%reserve gets some amount of CNP same as during the veg_init \n 
        !!  in the routine update_veg_pools section 3.0
      CASE ("population")
        IF(.NOT. model%config%l_do_stand_replacing_harvest) THEN 
          CALL calc_veg_establishment(options                         , & ! in
                                      lctlib                          , & 
                                      dsl4jsb_Config(ASSIMI_)%ncanopy , &
                                      w_soil_root_theta(:)            , &
                                      growing_season(:)               , &
                                      t_air_week_mavg(:)              , &
                                      fpc(:)                          , &
                                      rfr_ratio_boc_tvegdyn_mavg(:)   , &
                                      veg_pool%seed_bed               , &
                                      veg_establishment               )   ! out
        ENDIF
    ENDSELECT


    !> 4.0 Calculate change in individual density from mortality and establishment
    !! 
    !! only if no stand/cohort harvesting is done
    SELECT CASE (TRIM(dsl4jsb_Config(VEG_)%veg_dynamics_scheme))
    CASE ("population","cohort")
      !! @todo update_veg_dynamics() do_cohort_harvest dimension should be used as nc
      IF(.NOT. model%config%l_do_stand_replacing_harvest .AND. do_cohort_harvest(options%ics) < test_false_true) THEN 
        IF(lctlib%growthform == itree) THEN
           delta_dens_ind(:) = calc_delta_dens_ind(lctlib%seed_size           , & 
                                                   dens_ind(:)                , &
                                                   mortality_rate(:)          , &
                                                   veg_establishment%carbon(:)) 
        ELSE 
           delta_dens_ind(:) = 0.0_wp
        ENDIF
      ENDIF
    ENDSELECT

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
  !> calculates current foliage projective cover for a tile 
  !!
  !! Follows Sitch et al. 2003, eq 8)
  !-----------------------------------------------------------------------------------------------------
  FUNCTION calc_foliage_projective_cover_tile(options                 , &
                                              lctlib                  , &
                                              veg_pool_sap_wood_carbon, &
                                              dens_ind                , &
                                              height                  , &
                                              diameter) RESULT(fpc) 

    USE mo_jsb_task_class,        ONLY: t_jsb_task_options
    USE mo_jsb_lctlib_class,      ONLY: t_lctlib_element
    USE mo_veg_constants,         ONLY: max_crown_area, min_diameter, k_crown_area, k_rp, k_fpc, itree

    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    TYPE(t_jsb_task_options),        INTENT(in) :: options                     !< model options
    TYPE(t_lctlib_element),          INTENT(in) :: lctlib                      !< land-cover-type library - parameter across pft's
    REAL(wp), DIMENSION(options%nc), INTENT(in) :: veg_pool_sap_wood_carbon, & !< sapwood mass (mol C / m2)
                                                   dens_ind                , & !< tree density (#/m2)
                                                   height                  , & !< tree height (m)
                                                   diameter                    !< tree diameter (m)
    REAL(wp), DIMENSION(options%nc)             :: fpc                         !< foliage projective cover (fraction)
    ! ---------------------------
    ! 0.2 Local
    REAL(wp), DIMENSION(options%nc)             :: crown_area, &
                                                   lai_ind 
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_foliage_projective_cover_tile'


    WHERE(lctlib%growthform == itree .AND. dens_ind(:) .GT. eps4 .AND. height(:) .GT. eps4)
      ! crown area derived from diameter following Sitch et al. 2003; eq 4
      crown_area(:) = MIN(max_crown_area, k_crown_area * MAX(min_diameter,diameter(:)) ** k_rp)
      ! maximum annual, individual LAI, derived from the sapwood area of an individual tree, devided by 
      ! its crown area 
      lai_ind(:)    = veg_pool_sap_wood_carbon(:) / dens_ind(:) / lctlib%wood_density / height(:) * lctlib%k_latosa / & 
                      crown_area(:)
      ! foliage projective cover is then the product of density and individual crown area
      fpc(:)        = crown_area(:) * dens_ind(:) * (1._wp - exp(-k_fpc * lai_ind(:)))
    ELSEWHERE
      lai_ind(:) = 0.0_wp 
      fpc(:)     = 0.0_wp
    ENDWHERE
   
  END FUNCTION calc_foliage_projective_cover_tile


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
                                veg_dynamics_scheme       , &
                                an_boc_tvegdyn_mavg       , &
                                net_growth_tvegdyn_mavg   , &
                                lai_tvegdyn_mavg          , &
                                fpc                       , &
                                veg_pool                  , &
                                veg_litterfall            , &
                                mortality_rate )

    USE mo_base_constants,        ONLY: one_day, one_year, dtime
    USE mo_jsb_task_class,        ONLY: t_jsb_task_options
    USE mo_veg_constants,         ONLY: min_greff, k1_mort_greff, k2_mort_greff, k3_mort_greff, fpc_max, &
                                        itree, igrass, background_mort_rate_tree, background_mort_rate_grass
    USE mo_jsb_lctlib_class,      ONLY: t_lctlib_element
    dsl4slm_Use_memory_a_veg

    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    TYPE(t_jsb_task_options),        INTENT(in)    :: options                    !< model options
    TYPE(t_lctlib_element),          INTENT(in)    :: lctlib                     !< land-cover-type library - parameter across pft's
    CHARACTER(len=*),                INTENT(in)    :: veg_dynamics_scheme        !< vegetation dynamics: none population cohort
    REAL(wp), DIMENSION(options%nc), INTENT(in)    :: an_boc_tvegdyn_mavg    , & !< long-term net C balance at the bottom of the canopy (µmol / m2 / s)
                                                      net_growth_tvegdyn_mavg, & !< long-term net growth of the plant (µmol / m2 / s) 
                                                      lai_tvegdyn_mavg,        & !< long-term LAI of the plant
                                                      fpc                        !< foliage projective cover
    TYPE(t_veg_pool_a),              INTENT(in)    :: veg_pool                   !< vegetation pool
    TYPE(t_veg_pool_a),              INTENT(inout) :: veg_litterfall             !< vegetation litterfall
    REAL(wp), DIMENSION(options%nc), INTENT(out)   :: mortality_rate             !< current mortality rate (1/timestep)
    ! ---------------------------
    ! 0.2 Local
    REAL(wp), DIMENSION(options%nc)                :: boc_cbalance               !< C balance at bottom of canopy
    REAL(wp), DIMENSION(options%nc)                :: greff                      !< long-term growth efficiency
    REAL(wp), DIMENSION(options%nc)                :: fscal
    REAL(wp), DIMENSION(options%nc)                :: greff_mortality_rate, &
                                                      self_thinning_mortality_rate, &
                                                      bg_mortality_rate 
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_veg_mortality'
    

    !> differ between veg_dynamics_schemes
    !! population/cohort & none
    !! 
    SELECT CASE (veg_dynamics_scheme)
      CASE ("population","cohort")
        !> 1.0 mortality related to growth efficiency mortality 
        !! 
        WHERE(lai_tvegdyn_mavg(:) > eps8)
          !> growth efficiency per unit leaf area in mol C / m2 LAI / yr
          !!  the application of 'MAX()' avoids negative values of greff_mortality_rate(:) and mortality_rate(:)
          greff(:) = MAX(net_growth_tvegdyn_mavg(:) / lai_tvegdyn_mavg(:) * one_day * one_year / 1.e6_wp, 0.0_wp)
          !> 1.1 deduce mortality from asymptoptic mortality rate given growth efficiency, 
          !!     below a minimum threshold mortality rises to 100%
          !!
          WHERE(greff(:) > min_greff)
            greff_mortality_rate(:) = k1_mort_greff / (k2_mort_greff * greff(:) + 1.0_wp)
          ELSEWHERE
            greff_mortality_rate(:) = 1.0_wp / (k3_mort_greff * greff(:) + 1.0_wp)
          ENDWHERE
        ENDWHERE 

        !> 1.2 self-thinning related mortality
        !! 
        self_thinning_mortality_rate(:) = MAX(fpc(:) - fpc_max,0.0_wp)

!!       !! depends on the carbon balance of lowest canopy layer, i.e. net photosynthesis minus average construction cost \n
!!       !! in mol / m2 LAI / year \n
!!       !! transformed to a mortality estimate using a Weibull function
!!       boc_cbalance(:) = an_boc_tvegdyn_mavg(:) * one_day * one_year / 1.e6_wp &  
!!                         - ( 1._wp + fresp_growth ) * 1.0_wp/lctlib%tau_leaf/lctlib%sla
!!       WHERE(boc_cbalance(:) < (-eps4))
!!          mortality_rate(:) = mortality_rate(:) & 
!!                              + ( 1._wp - exp ( - ( - lambda_mort_light * boc_cbalance(:) ) ** k_mort_light))
!!       END WHERE

        !> 1.3 quasi-stochastic background mortality
        !!
        SELECT CASE (lctlib%growthform)
          CASE (itree)
            bg_mortality_rate(:) = background_mort_rate_tree
          CASE (igrass)
            bg_mortality_rate(:) = background_mort_rate_grass
        ENDSELECT
 
        !> 1.4 total mortality is the sum of all mortality terms, and limited to one
        !! 
        mortality_rate(:) = MIN(greff_mortality_rate(:)+self_thinning_mortality_rate(:)+bg_mortality_rate(:),1.0_wp) & 
                           / one_year / one_day * dtime 

      CASE ("none")
        SELECT CASE (lctlib%growthform)
          CASE (itree)
            mortality_rate(:) = ((background_mort_rate_tree /  one_year) / one_day) * dtime
          CASE (igrass)
            mortality_rate(:) = ((background_mort_rate_grass / one_year) / one_day) * dtime
        ENDSELECT
    ENDSELECT ! veg_dynamics_scheme


    !> 2.0 calculate litter fall from mortality
    !!
    !! for all tissues but seed_bed
    WHERE(veg_pool%leaf%carbon(:) > eps8)
       fscal(:) = (veg_pool%leaf%carbon(:) - veg_litterfall%leaf%carbon(:)) / veg_pool%leaf%carbon(:)
    ELSEWHERE
       fscal(:) = 0.0_wp 
    ENDWHERE
    WHERE(mortality_rate(:) < fscal(:))
       fscal(:) = mortality_rate(:)
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

    USE mo_jsb_task_class,        ONLY: t_jsb_task_options
    USE mo_veg_constants,         ONLY: wood_extraction_rate
    dsl4slm_Use_memory_a_veg

    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    TYPE(t_jsb_task_options),         INTENT(in)    :: options                  !< model options
    TYPE(t_veg_pool_a),               INTENT(in)    :: veg_pool                 !< vegetation pool
    TYPE(t_veg_pool_a),               INTENT(inout) :: veg_litterfall           !< vegetation litterfall
    REAL(wp), DIMENSION(options%nc),  INTENT(inout) :: mortality_rate           !< background_mort_rate_tree/grass converted to 1/timestep
    ! ---------------------------
    ! 0.2 Local
    REAL(wp), DIMENSION(options%nc)                 :: rest_mortality_rate 
    REAL(wp)                                        :: export_rate
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_veg_harvest'


    !> 1.0 calculate rest mortality to achieve stand-replacing harvest
    !!
    rest_mortality_rate(:)      = 1._wp - mortality_rate(:)

    veg_litterfall%leaf         = veg_litterfall%leaf         + veg_pool%leaf         * rest_mortality_rate(:)
    veg_litterfall%fine_root    = veg_litterfall%fine_root    + veg_pool%fine_root    * rest_mortality_rate(:)
    veg_litterfall%coarse_root  = veg_litterfall%coarse_root  + veg_pool%coarse_root  * rest_mortality_rate(:)
    veg_litterfall%sap_wood     = veg_litterfall%sap_wood     + veg_pool%sap_wood     * (1._wp - wood_extraction_rate) * rest_mortality_rate(:)
    veg_litterfall%heart_wood   = veg_litterfall%heart_wood   + veg_pool%heart_wood   * (1._wp - wood_extraction_rate) * rest_mortality_rate(:)
    veg_litterfall%labile       = veg_litterfall%labile       + veg_pool%labile       * (1._wp - wood_extraction_rate) * rest_mortality_rate(:)
    veg_litterfall%reserve      = veg_litterfall%reserve      + veg_pool%reserve      * (1._wp - wood_extraction_rate) * rest_mortality_rate(:)
    veg_litterfall%fruit        = veg_litterfall%fruit        + veg_pool%fruit        * rest_mortality_rate(:)

    mortality_rate(:)           = 1.0_wp

  END SUBROUTINE calc_veg_harvest


  !-----------------------------------------------------------------------------------------------------
  ! Sub Task of update_veg_dynamics
  ! 
  !-----------------------------------------------------------------------------------------------------
  !> calculate vegetation establishment 
  !!
  !!
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE calc_veg_establishment(options                     , &
                                    lctlib                      , &
                                    ncanopy                     , &
                                    w_soil_root_theta           , &
                                    growing_season              , &
                                    t_air_week_mavg             , &
                                    fpc                         , & 
                                    rfr_ratio_boc_tvegdyn_mavg  , &
                                    veg_pool_seed_bed           , &
                                    veg_establishment           )

    USE mo_jsb_task_class,        ONLY: t_jsb_task_options
    USE mo_jsb_lctlib_class,      ONLY: t_lctlib_element
    USE mo_base_constants,        ONLY: dtime, one_day, one_year
    USE mo_rad_constants,         ONLY: rfr_ratio_toc 
    USE mo_veg_constants,         ONLY: fpc_max, lambda_est_temp, k_est_temp, lambda_est_moist, k_est_moist
    USE mo_jsb_phys_constants,    ONLY: Tzero
    dsl4slm_Use_memory_a_veg

    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    TYPE(t_jsb_task_options),         INTENT(in)    :: options                      !< model options
    TYPE(t_lctlib_element),           INTENT(in)    :: lctlib                       !< land-cover-type library - parameter across pft's
    INTEGER,                          INTENT(in)    :: ncanopy                      !< number of canopy layer
    REAL(wp), DIMENSION(options%nc),  INTENT(in)    :: w_soil_root_theta            !< plant available water [fraction of maximum] 
    REAL(wp), DIMENSION(options%nc),  INTENT(in)    :: growing_season               !< growing season 
    REAL(wp), DIMENSION(options%nc),  INTENT(in)    :: t_air_week_mavg              !< weekly air temperature [K] 
    REAL(wp), DIMENSION(options%nc),  INTENT(in)    :: fpc                          !<  ...
    REAL(wp), DIMENSION(options%nc),  INTENT(in)    :: rfr_ratio_boc_tvegdyn_mavg   !< red-farred ratio at the bottom of the canopy 
    TYPE(t_veg_element_a),            INTENT(in)    :: veg_pool_seed_bed            !< seed_bed in vegetation pool
    TYPE(t_veg_element_a),            INTENT(inout) :: veg_establishment            !< vegetation establishment
    ! ---------------------------
    ! 0.2 Local
    REAL(wp), DIMENSION(options%nc)                 :: establishment_rate , &
                                                       tc                 , &
                                                       rfr_ratio          , &
                                                       flim_light         , &
                                                       flim_moist         , &
                                                       flim_temp
    INTEGER                                         :: icanopy
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_veg_establishment'


    !> 1.0 establishment rate
    !!

    !> 1.1 light limitation, 
    !!
    flim_light(:) = MAX(fpc_max - fpc(:),0.0_wp)

!!    !! depending on the red to far-red ratio of the light at the bottom of the canopy
!!    flim_light(:) = exp(-(lctlib%lambda_est_light * & 
!!                          (rfr_ratio_toc-rfr_ratio_boc_tvegdyn_mavg(:)))**lctlib%k_est_light)

    !> 1.2 temperature and soil moisture limitation 
    !! 
    !! OBS: should be top soil moisture!
    WHERE((t_air_week_mavg(:) - Tzero) < 0.0_wp .OR. growing_season(:) < test_false_true)
      flim_temp(:) = 0.0_wp
    ELSEWHERE
      flim_temp(:) = 1.0_wp - exp(-(lambda_est_temp * (t_air_week_mavg(:) - Tzero))**k_est_temp)
    END WHERE
    flim_moist(:)  = 1.0_wp - exp(-(lambda_est_moist * w_soil_root_theta(:))**k_est_moist) 
   
    !> 1.3 actual establishment rate and associated matter flux from seed bed to reserve pool
    !!
    establishment_rate(:) = flim_light(:) * flim_temp(:) * flim_moist(:) * & 
                            1.0_wp/lctlib%tau_seed_est/one_day/one_year * dtime
    WHERE(establishment_rate(:) < eps8)
       establishment_rate(:) = 0.0_wp
    ENDWHERE

    veg_establishment     = veg_pool_seed_bed * establishment_rate(:)

  END SUBROUTINE calc_veg_establishment


  !-----------------------------------------------------------------------------------------------------
  ! Sub Task of update_veg_dynamics
  ! 
  !-----------------------------------------------------------------------------------------------------
  !> calculate change in tree density 
  !!
  !-----------------------------------------------------------------------------------------------------
  PURE ELEMENTAL FUNCTION calc_delta_dens_ind(seed_size                , & ! lctlib
                                              dens_ind                 , &
                                              mortality_rate           , &
                                              veg_establishment_carbon)  &
                                              RESULT (delta_dens_ind)

    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    REAL(wp), INTENT(in) :: seed_size                , & !< seed size to convert C into individuals, lctlib parameter
                            dens_ind                 , & !< current individuum density (#/m2) 
                            mortality_rate           , & !< mortality per timestep
                            veg_establishment_carbon     !< vegetation establishment measured in C
    REAL(wp)             :: delta_dens_ind               !< change in individual density
    ! ---------------------------
    ! 0.2 Local
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_delta_density_individuals'


    delta_dens_ind = (veg_establishment_carbon / seed_size) - (mortality_rate * dens_ind)

  END FUNCTION calc_delta_dens_ind

END MODULE mo_veg_dynamics
