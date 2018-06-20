!>
!! @brief module to determine the phenological state of a plant 
!!
!! @author
!!  Sönke Zaehle, MPI-BGC Jena
!!
!! @par Revision History
!! First version by Sönke Zaehle (2016-04-11)
!! Second version by Jan Engel (20161214): just moved from vegetation to phenology making it an individual process
!!
MODULE mo_pheno_update_phenology

#include "dsl_jsbach.inc"

  USE mo_kind,          ONLY: wp, eps4
  USE mo_jsb_control,   ONLY: debug_on
  USE mo_exception,     ONLY: message

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: update_phenology

  CHARACTER(len=*), PARAMETER :: modname = 'mo_pheno_update_phenology'

CONTAINS

  !-----------------------------------------------------------------------------------------------------
  ! Main Task
  ! 
  !-----------------------------------------------------------------------------------------------------
  !> determine onset and end of the growing season for some PFT
  !! 
  !! PFT: evergreen, raingreen, cold-deciduous trees and perennial grasses \n
  !! 
  !! @todo OBS: currently uncalibrated, in particular rain green phenology \n
  !!  Dev: Have more cost-benefit in decision to start/stop growing season \n
  !!      what about seasonality beyond off/on behaviour? -> implies some collaboration between
  !!      turnover rates and phenology -> should the two be linked? 
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE update_phenology(tile, options)

    USE mo_jsb_class,               ONLY: Get_model
    USE mo_jsb_tile_class,          ONLY: t_jsb_tile_abstract
    USE mo_jsb_task_class,          ONLY: t_jsb_task_options
    USE mo_jsb_model_class,         ONLY: t_jsb_model
    USE mo_jsb_lctlib_class,        ONLY: t_lctlib_element
    USE mo_jsb_process_class,       ONLY: A2L_, L2A_, SEB_, TURB_, SOIL_, HYDRO_, VEG_, HD_, RAD_, ASSIMI_, PHENO_, SB_

    ! Use of process configurations (t_PROC_config)
    dsl4jsb_Use_config(PHENO_)    

    ! Use of process memories
    dsl4jsb_Use_memory(A2L_)        !  USE mo_a2l_memory_class, ONLY: t_a2l_memory
    dsl4jsb_Use_memory(VEG_)        !  USE mo_veg_memory_class,      ONLY: t_veg_memory
    dsl4jsb_Use_memory(PHENO_)      !  USE mo_pheno_memory_class,    ONLY: t_pheno_memory
    dsl4jsb_Use_memory(ASSIMI_)     !  USE mo_assimi_memory_class,   ONLY: t_assimi_memory

    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    CLASS(t_jsb_tile_abstract), INTENT(inout)     :: tile         !< one tile with data structure for one lct
    TYPE(t_jsb_task_options),   INTENT(in)        :: options      !< model options
    ! ---------------------------
    ! 0.2 Local
    TYPE(t_jsb_model),      POINTER       :: model
    TYPE(t_lctlib_element), POINTER       :: lctlib     !< land-cover-type library - parameter across pft's
    INTEGER                               :: iblk, ics, ice, nc
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':update_phenology'
    ! ---------------------------
    ! 0.3 Declare Memory
    ! Declare process configuration and memory Pointers
    dsl4jsb_Def_config(PHENO_)
    dsl4jsb_Def_memory(A2L_)                        !     TYPE(t_a2l_memory),      POINTER :: a2l__mem
    dsl4jsb_Def_memory(VEG_)                        !     TYPE(t_veg_memory),      POINTER :: veg__mem
    dsl4jsb_Def_memory(PHENO_)                      !     TYPE(t_pheno_memory),    POINTER :: pheno__mem
    dsl4jsb_Def_memory(ASSIMI_)                     !     TYPE(t_assimi_memory),   POINTER :: assimi__mem
    ! Declare pointers to variables in memory
    ! A2L_
    dsl4jsb_Real2D_onChunk      :: t_air
    ! PHENO_
    dsl4jsb_Logical2D_onChunk   :: growing_season
    dsl4jsb_Real2D_onChunk      :: gdd
    dsl4jsb_Real2D_onChunk      :: nd_dormance
    ! ASSIMI_
    dsl4jsb_Real2D_onChunk      :: beta_soil_gs_tphen_mavg
    ! VEG_
    dsl4jsb_Real2D_onChunk      :: t_air_tphen_mavg
    dsl4jsb_Real2D_onChunk      :: t_air_week_mavg
    dsl4jsb_Real2D_onChunk      :: t_air_month_mavg
    dsl4jsb_Real2D_onChunk      :: mean_leaf_age
    dsl4jsb_Real2D_onChunk      :: gpp_tlabile_mavg
    dsl4jsb_Real2D_onChunk      :: maint_respiration_tlabile_mavg
    dsl4jsb_Real2D_onChunk      :: growth_req_n_tlabile_mavg
    dsl4jsb_Real2D_onChunk      :: growth_req_p_tlabile_mavg
    ! Get local variables from options argument
    iblk    = options%iblk
    ics     = options%ics
    ice     = options%ice
    nc      = options%nc
    ! ---------------------------
    ! 0.4 Process Activity, Debug Option
    IF (.NOT. tile%Is_process_active(PHENO_)) RETURN
    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')
    ! ---------------------------
    ! 0.5 Get Memory
    ! Get process variables (Set pointers to variables in memory)
    model  => Get_model(tile%owner_model_id)
    lctlib => model%lctlib(tile%lcts%lib_id) 
    ! Get process config
    dsl4jsb_Get_config(PHENO_)          ! PHENO_conf => model%processes(PROC)%p%config | use: dsl4jsb_Config(PHENO_)%var or PHENO__conf%var
    ! Get process memories
    dsl4jsb_Get_memory(A2L_)            ! a2l__mem => a2l_memory
    dsl4jsb_Get_memory(VEG_)            ! veg__mem => veg_memory
    dsl4jsb_Get_memory(PHENO_)          ! pheno__mem => pheno_memory
    dsl4jsb_Get_memory(ASSIMI_)
    ! ---------------------------
    dsl4jsb_Get_var2D_onChunk(PHENO_, growing_season)                 ! inout
    dsl4jsb_Get_var2D_onChunk(PHENO_, gdd)                            ! inout
    dsl4jsb_Get_var2D_onChunk(PHENO_, nd_dormance)                    ! inout
    ! ---------------------------
    dsl4jsb_Get_var2D_onChunk(A2L_, t_air)                            ! in
    ! ---------------------------
    dsl4jsb_Get_var2D_onChunk(ASSIMI_, beta_soil_gs_tphen_mavg)       ! in
    ! ---------------------------
    dsl4jsb_Get_var2D_onChunk(VEG_, t_air_tphen_mavg)                 ! in
    dsl4jsb_Get_var2D_onChunk(VEG_, t_air_week_mavg)                  ! in
    dsl4jsb_Get_var2D_onChunk(VEG_, t_air_month_mavg)                 ! in
    dsl4jsb_Get_var2D_onChunk(VEG_, mean_leaf_age)                    ! in
    dsl4jsb_Get_var2D_onChunk(VEG_, gpp_tlabile_mavg)                 ! in
    dsl4jsb_Get_var2D_onChunk(VEG_, maint_respiration_tlabile_mavg)   ! in
    dsl4jsb_Get_var2D_onChunk(VEG_, growth_req_n_tlabile_mavg)        ! inout
    dsl4jsb_Get_var2D_onChunk(VEG_, growth_req_p_tlabile_mavg)        ! inout



    ! ------------------------------------------------------------------------------------------------------------
    ! Go science
    ! ------------------------------------------------------------------------------------------------------------


    !> 1.0 Check for start and end of growing season, and update both gdd and nd_dormance counter
    !! 
    CALL calc_phenology(lctlib%phenology_type              , &  ! lctlib
                        lctlib%beta_soil_flush             , &
                        lctlib%cn_leaf                     , &
                        lctlib%np_leaf                     , &
                        lctlib%gdd_req_max                 , &
                        lctlib%gdd_t_air_threshold         , &
                        lctlib%k_gdd_dormance              , &
                        lctlib%beta_soil_senescence        , &
                        lctlib%t_air_senescence            , &
                        lctlib%min_leaf_age                , &  ! lctlib
                        t_air(:)                           , &  ! in
                        beta_soil_gs_tphen_mavg(:)         , &
                        t_air_tphen_mavg(:)                , &
                        t_air_week_mavg(:)                 , &
                        t_air_month_mavg(:)                , &
                        mean_leaf_age(:)                   , &
                        gpp_tlabile_mavg(:)                , &
                        maint_respiration_tlabile_mavg(:)  , &
                        growing_season(:)                  , &  ! inout
                        growth_req_n_tlabile_mavg(:)       , &
                        growth_req_p_tlabile_mavg(:)       , &
                        nd_dormance(:)                     , &
                        gdd(:)                             )

  END SUBROUTINE update_phenology 


  !-----------------------------------------------------------------------------------------------------
  ! Sub Task of update_phenology
  ! 
  !-----------------------------------------------------------------------------------------------------
  !> calculate phenology 
  !!
  !!
  !-----------------------------------------------------------------------------------------------------
  ELEMENTAL SUBROUTINE calc_phenology(phenology_type                  , &  ! lctlib
                                      beta_soil_flush                 , &
                                      cn_leaf                         , &
                                      np_leaf                         , &
                                      gdd_req_max                     , &
                                      gdd_t_air_threshold             , &
                                      k_gdd_dormance                  , &
                                      beta_soil_senescence            , &
                                      t_air_senescence                , &
                                      min_leaf_age                    , &  ! lctlib
                                      t_air                           , &
                                      beta_soil_gs_tphen_mavg         , &
                                      t_air_tphen_mavg                , &
                                      t_air_week_mavg                 , &
                                      t_air_month_mavg                , &
                                      mean_leaf_age                   , &
                                      gpp_tlabile_mavg                , &
                                      maint_respiration_tlabile_mavg  , &
                                      growing_season                  , &
                                      growth_req_n_tlabile_mavg       , &
                                      growth_req_p_tlabile_mavg       , &
                                      nd_dormance                     , &
                                      gdd )

    USE mo_base_constants,          ONLY: dtime, one_day
    USE mo_shared_constants,        ONLY: Tzero
    USE mo_pheno_constants,         ONLY: ievergreen, iraingreen, isummergreen, iperennial

    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    INTEGER,                      INTENT(in)    :: phenology_type                      !< land-cover-type library parameter
    REAL(wp),                     INTENT(in)    :: beta_soil_flush                 , & !< land-cover-type library parameter
                                                   cn_leaf                         , & !< land-cover-type library parameter
                                                   np_leaf                         , & !< land-cover-type library parameter
                                                   gdd_req_max                     , & !< land-cover-type library parameter
                                                   gdd_t_air_threshold             , & !< land-cover-type library parameter
                                                   k_gdd_dormance                  , & !< land-cover-type library parameter
                                                   beta_soil_senescence            , & !< land-cover-type library parameter
                                                   t_air_senescence                , & !< land-cover-type library parameter
                                                   min_leaf_age                        !< land-cover-type library parameter
    REAL(wp),                     INTENT(in)    :: t_air                           , & !< air temp
                                                   beta_soil_gs_tphen_mavg         , & !< soil ..
                                                   t_air_tphen_mavg                , & !< moving average air temperature over ...
                                                   t_air_week_mavg                 , & !< moving average air temperature over a week
                                                   t_air_month_mavg                , & !< moving average air temperature over a month
                                                   mean_leaf_age                   , & !< mean leaf age
                                                   gpp_tlabile_mavg                , & !< moving average gpp over tlabile
                                                   maint_respiration_tlabile_mavg      !< moving average maintenance respiration over tlabile
    LOGICAL,                      INTENT(inout) :: growing_season                      !< growing season Yes/No
    REAL(wp),                     INTENT(inout) :: growth_req_n_tlabile_mavg       , & !< moles N required for a unit of C growth to determine labile pool size
                                                   growth_req_p_tlabile_mavg       , & !< moles P required for a unit of N growth to determine labile pool size
                                                   nd_dormance                     , & !< number of days in dormancy since last growing season
                                                   gdd                                 !< growing degree days 
    ! ---------------------------
    ! 0.2 Local
    REAL(wp) :: hlp1                                      ! helper variable for tmp values
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_phenology'



    ! ------------------------------------------------------------------------------------------------------------
    ! Go science
    ! ------------------------------------------------------------------------------------------------------------


    !> 1.0 check growing season has started or ended
    !! 
    IF(.NOT.growing_season) THEN ! If not in growing season, check whether growing season has started
       SELECT CASE (phenology_type)
         CASE (ievergreen)       ! nothing to do here
         CASE (iraingreen)       ! if soil moisture is above specified threshold
            IF(beta_soil_gs_tphen_mavg > beta_soil_flush) THEN 
                  growing_season = .TRUE.
                  growth_req_n_tlabile_mavg = 1._wp / cn_leaf
                  growth_req_p_tlabile_mavg = 1._wp / np_leaf
            ENDIF
         CASE (isummergreen)     ! if GDD requirement is exceeded. The latter depends on the length of the 
                                 ! dormancy period (aka "chilling requirement"). Parameterisation follows ORCHIDEE
            hlp1 = gdd_req_max * exp( - k_gdd_dormance * nd_dormance ) 
            IF(gdd > hlp1) THEN 
                  growing_season = .TRUE.
                  growth_req_n_tlabile_mavg = 1._wp / cn_leaf
                  growth_req_p_tlabile_mavg = 1._wp / np_leaf
            ENDIF
         CASE (iperennial)       ! if temperature and soil moisture are not limiting growth
            hlp1 = gdd_req_max * exp( - k_gdd_dormance * nd_dormance ) 
            IF(gdd > hlp1 .AND. &
               beta_soil_gs_tphen_mavg > beta_soil_flush)THEN 
                  growing_season = .TRUE.
                  growth_req_n_tlabile_mavg = 1._wp / cn_leaf
                  growth_req_p_tlabile_mavg = 1._wp / np_leaf
            ENDIF
       END SELECT
    ELSE                         ! If in the growing season, check whether growing season has ended
       SELECT CASE (phenology_type)
         CASE (ievergreen)       ! nothing to do here
         CASE (iraingreen)       ! if soil moisture is below given threshold AND leaves have been out a given number of days
            IF(beta_soil_gs_tphen_mavg < beta_soil_senescence .AND. & 
               mean_leaf_age > min_leaf_age ) & 
                  growing_season = .FALSE.
         CASE (isummergreen)     ! if temperature are declining below a certain threshold AND 
            ! leaves have been out a given number of days
            IF(t_air_week_mavg < t_air_month_mavg .AND. &
               t_air_tphen_mavg < t_air_senescence .AND. & 
               mean_leaf_age > min_leaf_age) & 
                  growing_season = .FALSE.
         CASE (iperennial)       ! if carbon balance of the plant becomes negative OR leaves have been damaged by frost/drought AND 
                                 ! leaves have been out a given number of days
            IF((t_air_tphen_mavg < t_air_senescence .OR. & 
                beta_soil_gs_tphen_mavg < beta_soil_senescence .OR. &
                gpp_tlabile_mavg < maint_respiration_tlabile_mavg).AND. &
                mean_leaf_age > min_leaf_age) &
                  growing_season = .FALSE.
       END SELECT
    ENDIF   

    !> 2.1 gdd counter
    !! 
    IF(t_air > gdd_t_air_threshold) THEN
       gdd = gdd + ((t_air - gdd_t_air_threshold) * dtime/one_day)
    ENDIF
    IF(t_air_tphen_mavg < Tzero) gdd = 0.0_wp

    !> 2.2 nd_dormance counter
    !! 
    IF (growing_season) THEN 
      gdd         = 0.0_wp
      nd_dormance = 0.0_wp
    ELSE
      nd_dormance = nd_dormance + dtime/one_day 
    ENDIF

  END SUBROUTINE calc_phenology 

END MODULE mo_pheno_update_phenology
