!>
!! @brief Contains the routines to calculate the canopy layer C and N content, and routines using this functionality 
!!
!! @author
!!  Sönke Zaehle, Silvia Caldararu, MPI-BGC Jena 
!!
!! @par Revision History
!! First version by Sönke Zaehle (2016-29-08)
!!
!!
MODULE mo_veg_canopy

#include "dsl_jsbach.inc"

  USE mo_kind
  USE mo_jsb_lctlib_class,            ONLY: t_lctlib_element
  USE mo_exception,                   ONLY: message_text, finish
  USE mo_assimi_update_canopy_fluxes, ONLY: calc_photosynthesis

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: calc_canopy_layers, calc_marginal_canopy_flux_increment, calc_dir_optimal_cn_leaf

CONTAINS

  !-----------------------------------------------------------------------------------------------------
  ! Sub Task (used by update vegetation pools, and some routines inquiring the effect of 
  !           changes in leaf C or leaf N)
  ! 
  !-----------------------------------------------------------------------------------------------------
  !> calc canopy layers 
  !! updates the canopy layer content (LAI, cummulative LAI above layer-mid-point, foliar N, fractional
  !!                                   allocation to photosynthetic functioning)
  !! 
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE calc_canopy_layers(lctlib            , &
                                leaf_nitrogen     , &
                                lai               , &
                                ppfd_sunlit_cl     , &
                                ppfd_shaded_cl     , & 
                                fleaf_sunlit_cl   , &
                                t_air             , &
                                t_acclim          , &
                                press_srf         , &
                                co2_mixing_ratio  , &
                                ga                , & 
                                beta_air          , &
                                beta_soil_ps      , &
                                beta_sinklim_ps   , &
                                beta_soil_gs      , &
                                t_jmax_opt        , & 
                                lai_cl            , &
                                cumm_lai_cl       , &
                                leaf_nitrogen_cl          , &
                                fn_rub_cl         , &
                                fn_et_cl          , &
                                fn_pepc_cl        , &
                                fn_chl_cl         , &
                                fn_oth_cl         , &
                                inquire_n_fractions)

    USE mo_assimi_constants
    USE mo_veg_constants
    USE mo_assimi_config_class,     ONLY: assimi_config     ! e.g. ncanopy, canopy_layer_depth, flag_optimal_Nfraction


    IMPLICIT NONE
    ! 0.1 InOut
    TYPE(t_lctlib_element),INTENT(in)            :: lctlib                !< land-cover-type library
    REAL(wp),INTENT(in)                          :: leaf_nitrogen         !< total canopy N (mol/m2)
    REAL(wp),INTENT(in)                          :: lai                   !< total leaf area index
    REAL(wp),INTENT(in),OPTIONAL                 :: t_air                 !< air temperature (avg)
    REAL(wp),INTENT(in),OPTIONAL                 :: t_acclim              !< acclimation temperature (avg)
    REAL(wp),INTENT(in),OPTIONAL                 :: press_srf             !< air pressure (avg) 
    REAL(wp),INTENT(in),OPTIONAL                 :: co2_mixing_ratio      !< CO2 concentration (avg)
    REAL(wp),INTENT(in),OPTIONAL                 :: ga                    !< aerodynamic conductance (avg)
    REAL(wp),INTENT(in),OPTIONAL                 :: beta_air              !< beta air (avg)
    REAL(wp),INTENT(in),OPTIONAL                 :: beta_soil_ps          !< beta soil ps (avg)
    REAL(wp),INTENT(in),OPTIONAL                 :: beta_sinklim_ps       !< beta sinklim ps (avg)
    REAL(wp),INTENT(in),OPTIONAL                 :: beta_soil_gs          !< beta soil gs (avg)
    REAL(wp),INTENT(in),OPTIONAL                 :: t_jmax_opt            !< Jmax optimum temperature (avg)
    REAL(wp),DIMENSION(assimi_config%ncanopy),INTENT(inout),OPTIONAL :: ppfd_sunlit_cl !< absorbed Photosynthetically Active Photon Flux Density
                                                                          !!  on sunlit leafes (average over time)
    REAL(wp),DIMENSION(assimi_config%ncanopy),INTENT(inout),OPTIONAL :: ppfd_shaded_cl !< absorbed Photosynthetically Active Photon Flux Density
                                                                          !!  on shaded leafes (average over time)
    REAL(wp),DIMENSION(:),INTENT(inout),OPTIONAL :: fleaf_sunlit_cl       !< fraction of sunlit leaves
                                                                          !  (average over time)
    REAL(wp),DIMENSION(:),INTENT(OUT)            :: lai_cl                !< LAI of canopy layer
    REAL(wp),DIMENSION(:),INTENT(OUT)            :: cumm_lai_cl           !< cummlative LAI above canopy layer
    REAL(wp),DIMENSION(:),INTENT(inout)          :: leaf_nitrogen_cl              !< total N in canopy layer (mmol N/m2)
    REAL(wp),DIMENSION(:),INTENT(inout)          :: fn_rub_cl             !< fraction of N in Rubisco
    REAL(wp),DIMENSION(:),INTENT(inout)          :: fn_et_cl              !< fraction of N in Electron Transport
    REAL(wp),DIMENSION(:),INTENT(inout)          :: fn_pepc_cl            !< fraction of N in PepC
    REAL(wp),DIMENSION(:),INTENT(inout)          :: fn_chl_cl             !< fraction of N in Chloroplast
    REAL(wp),DIMENSION(:),INTENT(inout)          :: fn_oth_cl             !< fraction of N not associated with PS
    LOGICAL,OPTIONAL                             :: inquire_n_fractions   !< whether to return fractions or 
                                                                          !  also run optimisation 
    ! 0.2 Local
    REAL(wp)  :: lc                              ! LAI at top of the canopy layer
    REAL(wp)  :: leaf_nitrogen_tc                        ! top-of-the-canopy N content
    REAL(wp)  :: k0_fn_chl,k1_fn_chl           ! factors in chlorophyll profile calculation
    REAL(wp)  :: jmax2vcmax                      ! jmax to vcmax ratio
    REAL(wp)  :: f_pepc                          ! fraction of N in PEP Carboxylase
    LOGICAL   :: optimise_n_fractions
    INTEGER   :: icanopy
    ! 0.3 Memory allocation
    ! ...

    !> 0.4 some initialisation and input checks
    !!

    !! @todo use lctlib for pathway specific constants in calc_canopy_layers
    ! set pathway specific constants (BLARRP!! this should not be here -> lctlib)
    IF(lctlib%ps_pathway == ic4phot) THEN
        jmax2vcmax = jmax2vcmax_C4
        k0_fn_chl = k0_fn_chl_C4
        k1_fn_chl = k1_fn_chl_C4
        f_pepc     = k0_fn_pepc
    ELSE
        jmax2vcmax = jmax2vcmax_C3
        k0_fn_chl = k0_fn_chl_C3
        k1_fn_chl = k1_fn_chl_C3
        f_pepc     = 0.0_wp
    ENDIF

    ! configure options for fractions optimisation 
    optimise_n_fractions = assimi_config%flag_optimal_Nfraction
    IF (PRESENT(inquire_n_fractions)) optimise_n_fractions = .NOT.inquire_n_fractions

    ! check if all input is available to run optimal fractions
    IF(assimi_config%flag_optimal_Nfraction)THEN
       IF(.NOT.PRESENT(ppfd_sunlit_cl).OR..NOT.PRESENT(ppfd_shaded_cl).OR..NOT.PRESENT(fleaf_sunlit_cl)) THEN
          WRITE(message_text,'(a)') 'called with fractional optimisation flag, but without par[sunlit|shaded]_cl'
          CALL finish("mo_veg_canopy:calc_canopy_layer",message_text)
       ENDIF
    ENDIF
    IF(optimise_n_fractions)THEN
       IF(.NOT.PRESENT(t_air).OR..NOT.PRESENT(t_acclim).OR..NOT.PRESENT(press_srf).OR.& 
          .NOT.PRESENT(co2_mixing_ratio).OR..NOT.PRESENT(ga).OR.&
          .NOT.PRESENT(beta_soil_ps).OR..NOT.PRESENT(beta_sinklim_ps).OR..NOT.PRESENT(beta_soil_gs).OR.&
          .NOT.PRESENT(beta_air).OR..NOT.PRESENT(t_jmax_opt)) THEN

          WRITE(message_text,'(a)') 'called to optimise N fractions, but without [climate] variables'
          CALL finish("mo_veg_canopy:calc_canopy_layer",message_text)
       ENDIF
    ENDIF


    ! initialise output
    lai_cl(:) = 0._wp
    cumm_lai_cl(:) = 0._wp
    leaf_nitrogen_cl(:) = 0._wp

    !> 1.0 calculate top of the canopy N
    !!
    IF(lai > assimi_config%canopy_layer_depth(1)) THEN
       leaf_nitrogen_tc = leaf_nitrogen * 1000._wp * kn / (1._wp - exp(-kn * lai))
    ELSEIF(lai > eps8) THEN ! if LAI does not fill the first canopy layer, top canopy N = leaf N
       leaf_nitrogen_tc = leaf_nitrogen * 1000._wp /lai 
    ENDIF

    !> 2. calculate LAI and N content of canopy layers 
    !!
    IF(lai > eps8) THEN ! avoid any calculation if there are no leaves
       lc = 0._wp
       DO icanopy = 1,assimi_config%ncanopy
          !> 2.1 LAI of canopy layer 
          !!
          IF(lc+assimi_config%canopy_layer_depth(icanopy) < lai) THEN 
             ! if canopy layer is fully developed 
             lai_cl(icanopy) = assimi_config%canopy_layer_depth(icanopy)
          ELSE IF(lc < lai.AND.lc+assimi_config%canopy_layer_depth(icanopy) >= lai) THEN 
             ! if canopy layer is only partially filled assign difference between full LAI and lc to canopy layer
             lai_cl(icanopy) = lai - lc 
          ENDIF
          lc = lc + lai_cl(icanopy)
          ! make sure all LAI is represented by letting last layer go to the end of the canopy
          IF(icanopy==assimi_config%ncanopy .AND. lc < lai ) THEN
            lai_cl(assimi_config%ncanopy) = lai_cl(assimi_config%ncanopy) + lai - lc 
          END IF

          !> 2.2 cummulative LAI above mid-point of canopy layer
          !!
          IF(lai_cl(icanopy) > eps8)THEN
             IF(icanopy==1) THEN
                cumm_lai_cl(icanopy) = lai_cl(icanopy) / 2.0_wp
             ELSE
                cumm_lai_cl(icanopy) = cumm_lai_cl(icanopy-1) + lai_cl(icanopy-1) / 2.0_wp + lai_cl(icanopy) / 2.0_wp
             ENDIF
          ENDIF

          !> 2.3 N content of canopy layer
          IF(lai > assimi_config%canopy_layer_depth(1)) THEN
             leaf_nitrogen_cl(icanopy) = leaf_nitrogen_tc * exp(-kn * cumm_lai_cl(icanopy))
          ELSE ! if LAI does not fill the first canopy layer, top canopy N = leaf N
             leaf_nitrogen_cl(icanopy) = leaf_nitrogen_tc
          ENDIF
       ENDDO
    ENDIF


    !> 3.0 fractionation of N within canopy laer
    !!
    !> 3.1 fraction not associated with photosynthesis 
    WHERE(leaf_nitrogen_cl > 0.0_wp) 
       fn_oth_cl = lctlib%k0_fn_struc - k1_fn_struc * leaf_nitrogen_cl   
    ELSEWHERE
       fn_oth_cl = 1.0_wp
    END WHERE

      ! structural limitation to the indefinite increase in phot. N described by above equation 
      WHERE (fn_oth_cl < lctlib%fn_oth_min) fn_oth_cl = lctlib%fn_oth_min


    !> 3.2 photosynthetic fractions
    !! prescribed, or dynamically calculated if assimi_config%flag_optimal_Nfraction = .TRUE.
    IF (.NOT.assimi_config%flag_optimal_Nfraction) THEN
        !> 3.2.1 case of presribed fractions
        WHERE(lai_cl > eps8) ! canopy layers with leaves
           ! fraction in chlorophyll 
           fn_chl_cl = (k0_fn_chl - k1_fn_chl * exp(-kfn_chl * cumm_lai_cl)) / chl2n
           ! fraction in PepC
           fn_pepc_cl = f_pepc
           ! fraction of electron-transport is dependent on N in electron transport versus N in light harvesting
           fn_et_cl = (1._wp - (fn_oth_cl + fn_chl_cl + fn_pepc_cl )) / &
                      (1._wp + jmax2n/(jmax2vcmax * vcmax2n))
           ! fraction of Rubisco is determined as the remainder
           fn_rub_cl = 1._wp - fn_oth_cl - fn_chl_cl - fn_et_cl - fn_pepc_cl
        ELSEWHERE ! canopy layers without leaves
           leaf_nitrogen_cl   = 0.0_wp
           fn_oth_cl  = 1.0_wp
           fn_chl_cl  = 0.0_wp
           fn_et_cl   = 0.0_wp
           fn_rub_cl  = 0.0_wp
           fn_pepc_cl = 0.0_wp
        ENDWHERE
     ELSE
        !> 3.2.1 optimal case 
        DO icanopy = 1,assimi_config%ncanopy
           IF(lai_cl(icanopy) > eps8)THEN
              IF(fn_chl_cl(icanopy) < eps8)THEN ! new layer with new leaves
                 IF(icanopy == 1)THEN ! first layers: use static allocation as first guess
                    fn_chl_cl(icanopy) = (k0_fn_chl - k1_fn_chl * exp(-kfn_chl * lc)) / chl2n
                    ! fraction of electron-transport is dependent on N in electron transport versus N in light harvesting
                    fn_et_cl(icanopy) = (1._wp - (fn_oth_cl(icanopy) + fn_chl_cl(icanopy) + fn_pepc_cl(icanopy) )) / &
                                      (1._wp + jmax2n/(jmax2vcmax * vcmax2n))
                    ! fraction in PepC
                    fn_pepc_cl(icanopy)=f_pepc
                    ! fraction of Rubisco is determined as the remainder
                    fn_rub_cl(icanopy) = 1._wp - fn_oth_cl(icanopy) - fn_chl_cl(icanopy) - fn_et_cl(icanopy) - fn_pepc_cl(icanopy)

                 ELSE ! other layers: copy from above layer
                    fn_chl_cl(icanopy)        = fn_chl_cl(icanopy-1)
                    fn_et_cl(icanopy)         = fn_et_cl(icanopy-1)
                    fn_rub_cl(icanopy)        = fn_rub_cl(icanopy-1)
                    fn_pepc_cl(icanopy)       = fn_pepc_cl(icanopy-1)
                    ppfd_sunlit_cl(icanopy)    = ppfd_sunlit_cl(icanopy-1)
                    ppfd_shaded_cl(icanopy)    = ppfd_shaded_cl(icanopy-1) 
                    fleaf_sunlit_cl(icanopy)  = fleaf_sunlit_cl(icanopy-1)
                 ENDIF
              ENDIF

              IF(optimise_n_fractions)THEN

                 CALL optimise_canopy_n_fractions(lctlib                  , & ! in
                                                  leaf_nitrogen_cl(icanopy)       , &
                                                  ppfd_sunlit_cl(icanopy)  , &
                                                  ppfd_shaded_cl(icanopy)  , &
                                                  fleaf_sunlit_cl(icanopy), &
                                                  t_air                   , &
                                                  t_acclim                , &
                                                  press_srf               , &
                                                  co2_mixing_ratio        , &
                                                  ga                      , &
                                                  beta_air                , &
                                                  beta_soil_ps            , &
                                                  beta_sinklim_ps         , &
                                                  beta_soil_gs            , &
                                                  t_jmax_opt              , &
                                                  fn_oth_cl(icanopy)      , &
                                                  fn_chl_cl(icanopy)      , & ! inout
                                                  fn_rub_cl(icanopy)      , &
                                                  fn_et_cl(icanopy)       , &
                                                  fn_pepc_cl(icanopy) )           
              ENDIF
           ELSE ! no leaves in canopy layer
              leaf_nitrogen_cl(icanopy)   = 0.0_wp
              fn_oth_cl(icanopy)  = 1.0_wp
              fn_chl_cl(icanopy)  = 0.0_wp
              fn_et_cl(icanopy)   = 0.0_wp
              fn_rub_cl(icanopy)  = 0.0_wp
              fn_pepc_cl(icanopy) = 0.0_wp
           ENDIF
        ENDDO
     ENDIF

  END SUBROUTINE calc_canopy_layers

  !-----------------------------------------------------------------------------------------------------
  ! Sub Task to calc_canopy_layers 
  ! 
  ! --------------------------------------------------------------------------------------------------
  !> Subroutine to opimise fraction of leaf N allocated to each photosynthetic component 
  !! Calculates carbon assimilation (an_cl) by varying fractions in either directions and determines
  !! fractions that give maximum assimilation
  !-------------------------------------------------------------------------------------------------------
  ELEMENTAL SUBROUTINE optimise_canopy_n_fractions(lctlib     , &
                                    leaf_nitrogen_cl                  , &
                                    ppfd_sunlit_cl             , &
                                    ppfd_shaded_cl             , &
                                    fleaf_sunlit_cl           , &
                                    t_air                     , &
                                    t_acclim                  , &
                                    press_srf                 , &
                                    co2_mixing_ratio          , &
                                    ga                        , & 
                                    beta_air                  , &
                                    beta_soil_ps              , &
                                    beta_sinklim_ps           , &
                                    beta_soil_gs              , &
                                    t_jmax_opt                , &
                                    fn_oth_cl                 , &
                                    fn_chl_cl                 , &
                                    fn_rub_cl                 , &
                                    fn_et_cl                  , &
                                    fn_pepc_cl)                       

    USE mo_base_constants,      ONLY: dtime,one_day
    USE mo_veg_constants,       ONLY: delta_n_fraction

    IMPLICIT NONE
    ! 0.1 InOut
    TYPE(t_lctlib_element),INTENT(in) :: lctlib         !< land-cover-type library - parameters for one pft
    REAL(wp),              INTENT(in) :: leaf_nitrogen_cl, &    !< total canopy N content
                                       & ppfd_sunlit_cl,   &    !< Photosynthetically Active Photon Flux Density 
                                                                !!  on sunlit leaves (day-time running mean, muemol/ms/s) 
                                       & ppfd_shaded_cl,   &    !< Photosynthetically Active Photon Flux Density
                                                                !!  on shaded leaves (day-time running mean, muemol/ms/s) 
                                       & fleaf_sunlit_cl,  &    !< fraction of sunlit leaves (day-time running mean, muemol/ms/s) 
                                       & t_air,            &    !< temperature (day-time running mean, K)
                                       & t_acclim,         &    !< ..
                                       & press_srf,        &    !< air pressure (day-time running mean, Pa)
                                       & co2_mixing_ratio, &    !< CO2 mixing ratio (day-time running mean, ppm)
                                       & ga,               &    !< aerodynamic conductance (day-time running mean, m/s)
                                       & beta_air,         &    !< beta air (day-time running mean, --)
                                       & beta_soil_ps,     &    !< beta soil for PS (day-time running mean, --)
                                       & beta_sinklim_ps,  &    !< beta sinklim for PS (day-time running mean, --)
                                       & beta_soil_gs,     &    !< beta soil for GC (day-time running mean, --) 
                                       & t_jmax_opt             !< optimum temperature for Jmax (day-time running mean, degreeC)
    REAL(wp),              INTENT(in) :: fn_oth_cl              !< fraction of leaf N not associated with PS
    REAL(wp),              INTENT(inout) :: & 
                                       & fn_chl_cl,        &    !< fraction of leaf N in chlorophyll
                                       & fn_rub_cl,        &    !< fraction of leaf N in rubisco
                                       & fn_et_cl,         &    !< fraction of leaf N in electron transport
                                       & fn_pepc_cl             !< fraction of leaf N in electron transport


    ! 0.2 Local
    INTEGER :: i
    REAL(wp):: hlp
    REAl(wp) :: ag0_cl, an0_cl, r0_cl, ci0_cl,gs0_cl           ! photosynthesis outputs of initial guess
    REAL(wp) :: fn_chl0_cl, fn_rub0_cl, fn_et0_cl, fn_pepc0_cl ! N fractions of initial guess
    REAL(wp) :: ag_cl, an_cl, r_cl, ci_cl,  gs_cl              ! photosynthesis outputs of optimisation runs
    REAL(wp) :: fn_rub_opt_cl, fn_et_opt_cl, fn_pepc_opt_cl    ! optimal N fractions
    REAL(wp) :: m_rub, m_et, m_pepc                            ! specific PS rates

    !> 1.0 Set N fractions to previous values and calculate photosynthesis for avererage conditions 
    ! set fractions to previous values. In current code chl + et+ rub=1
    ! all other routines chl +et +rub +oth = 1

    IF(fn_oth_cl > eps8) THEN
       fn_chl0_cl   = fn_chl_cl/(1._wp-fn_oth_cl)
       fn_rub0_cl   = fn_rub_cl/(1._wp-fn_oth_cl)
       fn_et0_cl    = fn_et_cl/(1._wp-fn_oth_cl)
       fn_pepc0_cl  = fn_pepc_cl/(1._wp-fn_oth_cl)
    ENDIF


    fn_chl_cl   = fn_chl0_cl
    fn_et_cl    = fn_et0_cl
    fn_rub_cl   = fn_rub0_cl
    fn_pepc_cl  = fn_pepc0_cl

    ! calculate photosynthesis with previous fractions
    CALL calc_photosynthesis( t_air                         , &
                              press_srf                     , &
                              co2_mixing_ratio              , &
                              ga                            , & 
                              t_acclim                      , &
                              ppfd_sunlit_cl                 , &
                              ppfd_shaded_cl                 , &
                              fleaf_sunlit_cl               , & 
                              beta_air                      , &
                              beta_soil_ps                  , &
                              beta_sinklim_ps               , &
                              beta_soil_gs                  , &
                              fn_chl0_cl*(1._wp-fn_oth_cl)  , &
                              fn_et0_cl*(1._wp-fn_oth_cl)   , &
                              fn_rub0_cl*(1._wp-fn_oth_cl)  , & 
                              fn_pepc0_cl*(1._wp-fn_oth_cl) , &
                              leaf_nitrogen_cl                      , &
                              lctlib%ps_pathway             , & 
                              t_jmax_opt                    , &
                              lctlib%t_jmax_omega           , &
                              lctlib%g0                     , &
                              lctlib%g1                     , &
                              ag0_cl                        , &
                              an0_cl                        , &
                              r0_cl                         , &
                              gs0_cl                        , &
                              ci0_cl                        , &
                              m_rub                         , &
                              m_et                          , &
                              m_pepc )

    !> 2.0 Attempt to increase PS by increasing the N fraction in chlorophyll 
    !! 
    fn_chl_cl = MIN(fn_chl0_cl * (1.0_wp + delta_n_fraction*dtime/one_day), 1.0_wp)

    ! optimal N fraction assuming Aj=Av=Apepc
    CALL calc_coordination_n_fractions(fn_chl_cl,m_rub,m_et,m_pepc, & 
                                       fn_rub_opt_cl,fn_et_opt_cl,fn_pepc_opt_cl)

    ! check if changing fraction to optimal values is possible in one timestep
    CALL calc_update_n_fractions(delta_n_fraction*dtime/one_day,fn_chl_cl,& 
                                  fn_rub_opt_cl,fn_et_opt_cl,fn_pepc_opt_cl, & 
                                  fn_rub_cl,fn_et_cl,fn_pepc_cl) 
                   
    ! re-calculate photosynthesis
    CALL calc_photosynthesis(t_air,press_srf,co2_mixing_ratio,ga,t_acclim, & 
              ppfd_sunlit_cl,ppfd_shaded_cl,fleaf_sunlit_cl, & 
              beta_air,beta_soil_ps,beta_sinklim_ps,beta_soil_gs, &
              fn_chl_cl*(1._wp-fn_oth_cl),fn_et_cl*(1._wp-fn_oth_cl),fn_rub_cl*(1._wp-fn_oth_cl), & 
              fn_pepc_cl*(1._wp-fn_oth_cl),leaf_nitrogen_cl,lctlib%ps_pathway, & 
              t_jmax_opt,lctlib%t_jmax_omega,lctlib%g0,lctlib%g1, &
              ag_cl,an_cl,r_cl,gs_cl,ci_cl,hlp,hlp,hlp)
    
    !> 3.0 check if increased chlorophyll has increased photosynthesis. If not, try decreasing it
    IF (an_cl < an0_cl) THEN
       
       ! decrease N fraction in chlorophyll
       fn_chl_cl = MAX(fn_chl0_cl * (1._wp - delta_n_fraction*dtime/one_day), eps8)
       fn_et_cl = fn_et0_cl
       fn_rub_cl = fn_rub0_cl
       fn_pepc_cl = fn_pepc0_cl

       ! optimal N fraction assuming Aj=Av=Apepc
       CALL calc_coordination_n_fractions(fn_chl_cl,m_rub,m_et,m_pepc, & 
                                          fn_rub_opt_cl,fn_et_opt_cl,fn_pepc_opt_cl)
       ! check if changing fraction to optimal values is possible in one timestep
       CALL calc_update_n_fractions(delta_n_fraction*dtime/one_day,fn_chl_cl,& 
                                     fn_rub_opt_cl,fn_et_opt_cl,fn_pepc_opt_cl, & 
                                     fn_rub_cl,fn_et_cl,fn_pepc_cl) 
     
       ! re-calculate photosynthesis
       CALL calc_photosynthesis(t_air,press_srf,co2_mixing_ratio,ga,t_acclim, & 
            ppfd_sunlit_cl,ppfd_shaded_cl,fleaf_sunlit_cl,& 
            beta_air,beta_soil_ps,beta_sinklim_ps,beta_soil_gs, &
            fn_chl_cl*(1._wp-fn_oth_cl), fn_et_cl*(1._wp-fn_oth_cl), fn_rub_cl*(1._wp-fn_oth_cl), & 
            fn_pepc_cl*(1._wp-fn_oth_cl),leaf_nitrogen_cl, lctlib%ps_pathway, & 
            t_jmax_opt,lctlib%t_jmax_omega,lctlib%g0,lctlib%g1, &
            ag_cl,an_cl,r_cl,gs_cl,ci_cl,hlp, hlp,hlp)

       ! if photosyntheis has still not increased, do not change fn_chl_cl
       IF (an_cl < an0_cl) THEN
          fn_chl_cl = fn_chl0_cl
          fn_et_cl = fn_et0_cl
          fn_rub_cl = fn_rub0_cl
          fn_pepc_cl = fn_pepc0_cl

          ! optimal N fraction assuming Aj=Av=Apepc
          CALL calc_coordination_n_fractions(fn_chl_cl,m_rub,m_et,m_pepc, & 
                                             fn_rub_opt_cl,fn_et_opt_cl,fn_pepc_opt_cl)
          ! check if changing fraction to optimal values is possible in one timestep
          CALL calc_update_n_fractions(delta_n_fraction*dtime/one_day,fn_chl_cl,& 
                                        fn_rub_opt_cl,fn_et_opt_cl,fn_pepc_opt_cl, & 
                                        fn_rub_cl,fn_et_cl,fn_pepc_cl) 

       ENDIF
    ENDIF
    !> 4.0 Transform optimal fractions to fractions of total canopy N 
    fn_chl_cl = fn_chl_cl*(1._wp-fn_oth_cl)
    fn_rub_cl = fn_rub_cl*(1._wp-fn_oth_cl)
    fn_et_cl = fn_et_cl*(1._wp-fn_oth_cl)
    fn_pepc_cl = fn_pepc_cl*(1._wp-fn_oth_cl)

  END SUBROUTINE optimise_canopy_n_fractions

  !-----------------------------------------------------------------------------------------------------
  ! Sub Task to optimise_canopy_n_fractions
  ! 
  ! --------------------------------------------------------------------------------------------------
  !> Subroutine to calculate fractions of rubisco, electron transport and PepC assuming them to
  !! be coordinated (Aj=Av=Apepc) 
  !@todo: veg_canopy:calc_coordination_n_fractions: unsafe division by m_et/m_pepc
  !-------------------------------------------------------------------------------------------------------
  ELEMENTAL SUBROUTINE calc_coordination_n_fractions(fn_chl_cl,m_rub,m_et,m_pepc, &
                                                     fn_rub_opt_cl,fn_et_opt_cl,fn_pepc_opt_cl)
  
    IMPLICIT NONE
  
    ! 0.1 InOut
    REAL(wp), INTENT(in)    :: fn_chl_cl      , &  !< fraction of N in chlorophyll
                               m_rub          , &  !< N-specific rate constant for Rubisco
                               m_pepc              !< N-speficic rate constant for PepC (C4-only)
    REAL(wp), INTENT(inout) :: m_et                !< N-specific rate constant for ET
    REAL(wp), INTENT(out)   :: fn_rub_opt_cl  , &  !< optimal fraction of N in Rubisco 
                               fn_et_opt_cl   , &  !< optimal fraction of N in ET 
                               fn_pepc_opt_cl      !< optimal fraction of N in PepC (C4-only)

    ! optimal N fraction assuming Aj=Av=Apepc
!@todo: veg_canopy:optimise_n_fractions: should this not replaced by lctlib%ps_pathway==ic4phot?
    ! ensure m_et is > zero
    IF (m_et < eps8) THEN 
      m_et = eps8
    END IF
    IF (m_pepc > 0._wp) THEN
       ! fraction in Rubisco
       fn_rub_opt_cl = (1._wp - fn_chl_cl)/(1._wp + m_rub/m_et + m_rub/m_pepc)
       ! fraction in electron transport
       fn_et_opt_cl = fn_rub_opt_cl * m_rub/m_et
       ! fraction in PEPc
       fn_pepc_opt_cl = fn_rub_opt_cl * m_rub/m_pepc
    ELSE
       ! fraction in Rubisco
       fn_rub_opt_cl = (1._wp - fn_chl_cl)/(1._wp + (m_rub)/(m_et))
       ! fraction in electron transport
       fn_et_opt_cl = fn_rub_opt_cl * (m_rub)/(m_et)
       ! fraction in PEPc
       fn_pepc_opt_cl = 0._wp
    ENDIF

  END SUBROUTINE calc_coordination_n_fractions

  !-----------------------------------------------------------------------------------------------------
  ! Sub Task to optimise_canopy_n_fractions
  ! 
  ! --------------------------------------------------------------------------------------------------
  !> Subroutine to update N fractions towards their optimal values in the limit of the 
  !! impose maximum rate of change  
  !-------------------------------------------------------------------------------------------------------
  ELEMENTAL SUBROUTINE calc_update_n_fractions(delta_n_fraction,fn_chl_cl, &  
                                                fn_rub_opt_cl,fn_et_opt_cl,fn_pepc_opt_cl, & 
                                                fn_rub_cl,fn_et_cl,fn_pepc_cl)

    IMPLICIT NONE

    ! 0.1 InOut
    REAL(wp), INTENT(in) :: delta_n_fraction, &  !< maximum fractional change in this time-step
     &     fn_chl_cl,      &                     !< fraction of N in chlorophyll
     &     fn_rub_opt_cl,  &                     !< optimal fraction of N in Rubisco
     &     fn_et_opt_cl,   &                     !< optimal fraction of N in ET
     &     fn_pepc_opt_cl                        !< optimal fraction of N in PepC (C4-only)
    REAL(wp), INTENT(inout) :: fn_rub_cl, &      !< updated fraction of N in Rubisco 
     &     fn_et_cl,       &                     !< updated fraction of N in ET
     &     fn_pepc_cl                            !< updated fraction of N in PepC (C4-only)

    ! check if changing fraction to optimal values is possible in one timestep
    IF (ABS(fn_rub_opt_cl-fn_rub_cl) / fn_rub_cl <= delta_n_fraction .AND. &
        ABS(fn_et_opt_cl-fn_et_cl) / fn_et_cl <= delta_n_fraction .AND. &
        ABS(fn_pepc_opt_cl-fn_pepc_cl) / fn_pepc_cl <= delta_n_fraction) THEN
       fn_rub_cl = fn_rub_opt_cl
       fn_et_cl = fn_et_opt_cl
       fn_pepc_cl = fn_pepc_opt_cl
    ! if not, change by maximum amount into right direction
    ELSEIF (fn_rub_cl < fn_rub_opt_cl) THEN
       fn_rub_cl = MIN(fn_rub_cl * (1.0_wp+delta_n_fraction), fn_rub_opt_cl)
       IF (fn_pepc_cl < fn_pepc_opt_cl) THEN
          fn_pepc_cl = MAX(fn_pepc_cl * (1.0_wp+delta_n_fraction), fn_pepc_opt_cl)
          fn_et_cl = MAX(1.0_wp -fn_rub_cl - fn_pepc_cl - fn_chl_cl, eps8)
       ELSE
          fn_pepc_cl = MIN(fn_pepc_cl * (1.0_wp-delta_n_fraction), fn_pepc_opt_cl)
          fn_et_cl = MAX(1.0_wp -fn_rub_cl - fn_pepc_cl - fn_chl_cl, eps8)
       ENDIF
    ELSE
       fn_rub_cl = MIN(fn_rub_cl * (1.0_wp-delta_n_fraction), fn_rub_opt_cl)
       IF (fn_pepc_cl < fn_pepc_opt_cl) THEN
          fn_pepc_cl = MAX(fn_pepc_cl * (1.0_wp+delta_n_fraction), fn_pepc_opt_cl)
          fn_et_cl = MAX(1.0_wp -fn_rub_cl - fn_pepc_cl - fn_chl_cl, eps8)
       ELSE
          fn_pepc_cl = MIN(fn_pepc_cl * (1.0_wp-delta_n_fraction), fn_pepc_opt_cl)
          fn_et_cl = MAX(1.0_wp -fn_rub_cl - fn_pepc_cl - fn_chl_cl, eps8)
       ENDIF
    ENDIF

  END SUBROUTINE calc_update_n_fractions                   

  !-----------------------------------------------------------------------------------------------------
  ! Sub Task to calc_growth (currently in mo_veg_canopy to get rid of circular dependencies)
  ! 
  !------------------------------------------------------------------------------------------------------
  !> Subroutine to calculate change in net photosynthesis and transpiration with small change in LAI
  !! Needed for calculation of optimal leaf and root allocation fraction. See mo_veg_growth
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE calc_marginal_canopy_flux_increment(lctlib             , &
                                                  t_air             , &
                                                  press_srf         , & 
                                                  q_air             , & 
                                                  co2_mixing_ratio  , & 
                                                  ga                , &
                                                  net_assimilation  , &
                                                  canopy_cond       , &
                                                  beta_air          , &
                                                  beta_soil_ps      , &
                                                  beta_soil_gs      , &
                                                  rad_ppfd_sunlit_cl, &
                                                  rad_ppfd_shaded_cl, &
                                                  veg               , &
                                                  unit_npp          , &
                                                  unit_transpiration )

    USE mo_assimi_config_class,   ONLY: assimi_config
    USE mo_atmland_util,          ONLY: calc_spec_humidity_sat
    USE mo_shared_constants,      ONLY: r_gas_dryair

    ! Use of process memories
    dsl4jsb_Use_memory(VEG_)        !  USE mo_veg_memory_class,      ONLY: t_veg_memory

    IMPLICIT NONE
    ! 0.1 InOut
    TYPE(t_lctlib_element),INTENT(in) :: lctlib              !< land-cover-type library - parameters for one pft
    TYPE(t_veg_memory),    INTENT(in) :: veg                 !< veg process memory
    REAL(wp),              INTENT(in) :: t_air               !< air temperature (K) 
    REAL(wp),              INTENT(in) :: press_srf           !< air pressure (Pa) 
    REAL(wp),              INTENT(in) :: q_air               !< air specific humidity (--) 
    REAL(wp),              INTENT(in) :: co2_mixing_ratio    !< co2 mixing ratio (ppm) 
    REAL(wp),              INTENT(in) :: ga                  !< air temperature (K) 
    REAL(wp),              INTENT(in) :: net_assimilation    !< net assimilation rate (µmol CO2 / m2 / s) 
    REAL(wp),              INTENT(in) :: canopy_cond         !< canopy conductance (m/s) 
    REAL(wp),              INTENT(in) :: beta_air            !< air humidity limiting factor for photosynthesis
    REAL(wp),              INTENT(in) :: beta_soil_ps        !< soil moisture limiting factor for photosynthesis
    REAL(wp),              INTENT(in) :: beta_soil_gs        !< soil moisture limiting factor for canopy conductance 
    REAL(wp), DIMENSION(assimi_config%ncanopy), INTENT(in)    :: rad_ppfd_sunlit_cl , &
                                                                 rad_ppfd_shaded_cl
    REAL(wp),              INTENT(out) :: unit_npp           !<- unit change of NPP per ??
    REAL(wp),              INTENT(out) :: unit_transpiration !<- unit change of NPP per ??
    ! 0.2 Local
    REAL(wp) :: lai_new, leaf_nitrogen_new
    REAL(wp) :: net_assimilation_new 
    REAL(wp) :: canopy_cond_new
    REAL(wp),DIMENSION(assimi_config%ncanopy) :: lai_new_cl
    REAL(wp),DIMENSION(assimi_config%ncanopy) :: cumm_lai_new_cl 
    REAL(wp),DIMENSION(assimi_config%ncanopy) :: ppfd_sunlit_cl,ppfd_shaded_cl,fleaf_sunlit_cl 
    REAl(wp),DIMENSION(assimi_config%ncanopy) :: fn_oth_new_cl, leaf_nitrogen_new_cl, fn_chl_new_cl, fn_et_new_cl, fn_rub_new_cl, fn_pepc_new_cl
    REAL(wp) :: ag_cl, an_cl, r_cl, gs_cl
    REAL(wp) :: q_sat 
    REAL(wp) :: hlp1 
    INTEGER :: icanopy
 
    !> 1.0 copy current plant state into dummy variables
    !!
    !fn_chl_cl       = veg%fn_chl_cl%ptr(:,1)
    !fn_et_cl        = veg%fn_et_cl%ptr(:,1)
    !fn_rub_cl       = veg%fn_rub_cl%ptr(:,1)
    !fn_pepc_cl      = veg%fn_pepc_cl%ptr(:,1)
    ppfd_sunlit_cl   = rad_ppfd_sunlit_cl(:)
    ppfd_shaded_cl   = rad_ppfd_shaded_cl(:)
    fleaf_sunlit_cl = veg%fleaf_sunlit_cl%ptr(:,1)
 
 !> @todo: mo_veg_canopy calc_npp_increment: I think this is wrong currently, as increaseing leaf area require an
 !! increase in roots and sapwood mass and thus respiration... Also, should be based on current mass
 !! ratios, not a fixed constant. ie
 !! f_increment = f_alloc_leaf/f_alloc_root+sapwood
 !! maint_resp = maint_resp * (1 + f_alloc_root+sapwood*1/totalmassofpools)
 !! probably only relevant for BNF and uptake calculation?

    !> 2.0 Estimate effect of a small increase in leaf C on LAI and leaf N
    !!
    ! increase LAI by 1 mol of C 
    lai_new = veg%lai%ptr + lctlib%sla
    ! equivalent leaf N assuming any new leaf growth will have target CN
    IF(veg%pool%leaf%carbon%ptr > eps8)THEN
       leaf_nitrogen_new = veg%pool%leaf%nitrogen%ptr + veg%pool%leaf%nitrogen%ptr / veg%pool%leaf%carbon%ptr
    ELSE
       leaf_nitrogen_new = 1._wp / lctlib%cn_leaf
    ENDIF

    ! recompute canopy layer LAI and foliar N with new numbers
    CALL calc_canopy_layers(lctlib = lctlib, & 
                            leaf_nitrogen = leaf_nitrogen_new, &
                            lai = lai_new,&
                            ppfd_sunlit_cl = ppfd_sunlit_cl, &
                            ppfd_shaded_cl = ppfd_shaded_cl, &
                            fleaf_sunlit_cl = fleaf_sunlit_cl, &
                            lai_cl = lai_new_cl, &
                            cumm_lai_cl = cumm_lai_new_cl, &
                            leaf_nitrogen_cl = leaf_nitrogen_new_cl, &
                            fn_rub_cl = fn_rub_new_cl, &
                            fn_et_cl = fn_et_new_cl, &
                            fn_pepc_cl = fn_pepc_new_cl, &
                            fn_chl_cl = fn_chl_new_cl, &
                            fn_oth_cl = fn_oth_new_cl, &
                            inquire_n_fractions=.TRUE.)
 
    ! calculate photosynthesis by integration over canopy layers containing leaves 
    net_assimilation_new = 0._wp
    canopy_cond_new = 0.0_wp
    DO icanopy = 1,assimi_config%ncanopy
       IF(lai_new_cl(icanopy).GT.eps8)THEN   
          CALL calc_photosynthesis( t_air                               , & ! in
                                    press_srf                           , &
                                    co2_mixing_ratio                    , &
                                    ga                                  , &
                                    veg%t_air_tacclim_mavg%ptr          , &
                                    rad_ppfd_sunlit_cl(icanopy)    , &
                                    rad_ppfd_shaded_cl(icanopy)    , &
                                    veg%fleaf_sunlit_cl%ptr(icanopy,1)  , &
                                    beta_air                    , &
                                    beta_soil_ps                , &
                                    veg%beta_sinklim_ps%ptr             , &
                                    beta_soil_gs                , &
                                    veg%fn_chl_cl%ptr(icanopy,1)       , &
                                    veg%fn_et_cl%ptr(icanopy,1)        , &
                                    veg%fn_rub_cl%ptr(icanopy,1)       , &
                                    veg%fn_pepc_cl%ptr(icanopy,1)      , &
                                    leaf_nitrogen_new_cl(icanopy)                  , &
                                    lctlib%ps_pathway                   , &
                                    veg%t_jmax_opt_mavg%ptr             , &
                                    lctlib%t_jmax_omega                 , &
                                    lctlib%g0                           , &
                                    lctlib%g1                           , &
                                    ag_cl                               , & ! out
                                    an_cl                               , &
                                    r_cl                                , &
                                    gs_cl                               , &
                                    hlp1,hlp1,hlp1,hlp1)             
    
            net_assimilation_new = net_assimilation_new + an_cl * lai_new_cl(icanopy) 
            canopy_cond_new = canopy_cond_new + gs_cl * lai_new_cl(icanopy) 

       ENDIF
    ENDDO
    
    !> 3.0 Calculate change in canopy fluxes from altered LAI 

    !> 3.1 net carbon gain from increased LAI
    !!
    unit_npp = net_assimilation_new - net_assimilation 

    !> 3.2 Calculate increment in net water loss from increased LAI
    !!
    q_sat = calc_spec_humidity_sat(t_air,press_srf) 
    IF(canopy_cond > eps8 .AND. canopy_cond_new > eps8)THEN 
       unit_transpiration = press_srf / ( r_gas_dryair * t_air ) * &  
                        (q_sat - q_air) / (1._wp / canopy_cond_new + 1._wp / ga ) - & 
                        press_srf / ( r_gas_dryair * t_air ) * &  
                        (q_sat - q_air) / (1._wp / canopy_cond + 1._wp / ga )  
    ELSEIF (canopy_cond_new > eps8) THEN 
       unit_transpiration =  press_srf / ( r_gas_dryair * t_air ) * &  
                         (q_sat - q_air) / (1._wp / canopy_cond_new + 1._wp / ga ) 
    ELSE 
      unit_transpiration = 0.0_wp  
    ENDIF 

   END SUBROUTINE calc_marginal_canopy_flux_increment

  !-----------------------------------------------------------------------------------------------------
  ! Sub Task to calc_growth (currently in mo_veg_canopy to get rid of circular dependencies)
  ! 
  !------------------------------------------------------------------------------------------------------
  !> Subroutine to determine direction in which total leaf N should change in order to increase NPP
  !! 
  !! Parameter leaf_cn_direction
  !----------------------------------------------------------------------------------------------------
  SUBROUTINE calc_dir_optimal_cn_leaf(lctlib                          , &
                                        veg                           , &
                                        beta_air_tcnl_mavg            , &
                                        beta_soil_ps_tcnl_mavg        , &
                                        beta_soil_gs_tcnl_mavg        , &
                                        rad_ppfd_sunlit_tcnl_mavg_cl  , &
                                        rad_ppfd_shaded_tcnl_mavg_cl )

    USE mo_base_constants,        ONLY: dtime, one_day 
    USE mo_veg_constants,         ONLY: delta_n_leaf,root2leaf_cn,wood2leaf_cn 
    USE mo_veg_respiration,       ONLY: calc_maintenance_respiration
    USE mo_assimi_config_class,   ONLY: assimi_config

    ! Use of process memories
    dsl4jsb_Use_memory(VEG_)        !  USE mo_veg_memory_class,      ONLY: t_veg_memory

    IMPLICIT NONE
    ! 0.1 InOut
    TYPE(t_lctlib_element),                     INTENT(in)    :: lctlib       !< land-cover-type library - parameters for one pft
    TYPE(t_veg_memory),                         INTENT(inout) :: veg          !< veg process memory 
    REAL(wp),                                   INTENT(in)    :: beta_air_tcnl_mavg
    REAL(wp),                                   INTENT(in)    :: beta_soil_ps_tcnl_mavg
    REAL(wp),                                   INTENT(in)    :: beta_soil_gs_tcnl_mavg
    REAL(wp), DIMENSION(assimi_config%ncanopy), INTENT(in)    :: rad_ppfd_sunlit_tcnl_mavg_cl , &
                                                                 rad_ppfd_shaded_tcnl_mavg_cl
    ! 0.2 Local
    REAL(wp)                                   :: leaf_nitrogen,fine_root_nitrogen,coarse_root_nitrogen,sap_wood_nitrogen
    REAL(wp), DIMENSION(assimi_config%ncanopy) :: ppfd_sunlit_tcnl_mavg_cl,ppfd_shaded_tcnl_mavg_cl,fleaf_sunlit_tcnl_mavg_cl 
    REAL(wp), DIMENSION(assimi_config%ncanopy) :: leaf_nitrogen_cl,fn_oth_cl,fn_chl_cl,fn_et_cl, fn_rub_cl, fn_pepc_cl,lai_cl,cumm_lai_cl
    REAL(wp), DIMENSION(3)                     :: npp_direction, maint_respiration
    REAL(wp) :: anl, agl, rl
    REAL(wp)  :: hlp1
    REAL(wp)  :: direction
    INTEGER   :: icanopy, dd


    IF(veg%pool%leaf%nitrogen%ptr > eps4)THEN

    fn_chl_cl                 = veg%fn_chl_cl%ptr(:,1)
    fn_et_cl                  = veg%fn_et_cl%ptr(:,1)
    fn_rub_cl                 = veg%fn_rub_cl%ptr(:,1)
    ppfd_sunlit_tcnl_mavg_cl   = rad_ppfd_sunlit_tcnl_mavg_cl(:)
    ppfd_shaded_tcnl_mavg_cl   = rad_ppfd_shaded_tcnl_mavg_cl(:)
    fleaf_sunlit_tcnl_mavg_cl = veg%fleaf_sunlit_tcnl_mavg_cl%ptr(:,1)

    ! Calculate photosynthesis with lower, current and higher target leaf N to determine optimal canopy N
    npp_direction(:) = 0.0_wp
    DO dd = 1,3 
        
       direction = REAL(dd,KIND=wp) - 2.0_wp

       leaf_nitrogen = veg%pool%leaf%carbon%ptr / veg%target_cn_leaf%ptr * (1.0_wp + direction * delta_n_leaf * dtime/one_day)
       
       ! update N distribution in the canopy
       CALL calc_canopy_layers(lctlib = lctlib, & 
                          leaf_nitrogen = leaf_nitrogen, &
                          lai = veg%lai%ptr,&
                          ppfd_sunlit_cl = ppfd_sunlit_tcnl_mavg_cl, &
                          ppfd_shaded_cl = ppfd_shaded_tcnl_mavg_cl, &
                          fleaf_sunlit_cl = fleaf_sunlit_tcnl_mavg_cl, &
                          lai_cl = lai_cl, &
                          cumm_lai_cl = cumm_lai_cl, &
                          leaf_nitrogen_cl = leaf_nitrogen_cl, &
                          fn_rub_cl = fn_rub_cl, &
                          fn_et_cl = fn_et_cl, &
                          fn_pepc_cl = fn_pepc_cl, &
                          fn_chl_cl = fn_chl_cl, &
                          fn_oth_cl = fn_oth_cl, &
                          inquire_n_fractions=.TRUE.)

       ! calculate photosynthesis for each layer
       DO icanopy =1, assimi_config%ncanopy
          IF(lai_cl(icanopy) > eps8)THEN

             CALL calc_photosynthesis(veg%t_air_tcnl_mavg%ptr,& 
                       veg%press_srf_tcnl_mavg%ptr, &
                       veg%co2_mixing_ratio_tcnl_mavg%ptr,&
                       veg%ga_tcnl_mavg%ptr, & 
                       veg%t_air_tacclim_mavg%ptr, &
                       rad_ppfd_sunlit_tcnl_mavg_cl(icanopy),&
                       rad_ppfd_shaded_tcnl_mavg_cl(icanopy),& 
                       veg%fleaf_sunlit_tcnl_mavg_cl%ptr(icanopy,1), &
                       beta_air_tcnl_mavg,&
                       beta_soil_ps_tcnl_mavg,&
                       veg%beta_sinklim_ps_tcnl_mavg%ptr,&
                       beta_soil_gs_tcnl_mavg, &
                       veg%fn_chl_cl%ptr(icanopy,1)             , &
                       veg%fn_et_cl%ptr(icanopy,1)              , &
                       veg%fn_rub_cl%ptr(icanopy,1)             , &
                       veg%fn_pepc_cl%ptr(icanopy,1)            , &
                       leaf_nitrogen_cl(icanopy),&
                       lctlib%ps_pathway, & 
                       veg%t_jmax_opt_mavg%ptr, &                  
                       lctlib%t_jmax_omega,lctlib%g0,lctlib%g1, &
                       agl,anl,rl,hlp1,hlp1,&
                       hlp1,hlp1,hlp1)
               
             npp_direction(dd) = npp_direction(dd) + anl * lai_cl(icanopy)
          
          ENDIF
       ENDDO
     
     ENDDO  ! dd = 1,3 
        
     ! Direction of change in leaf N. Actual leaf N is modified in interfacte routine
!@todo:veg_canopy:calc_dir_optimal_cn_leaf:potentially dangerous use of MAXLOC - ascertain functiong in 4D
!@todo:veg%leaf_cn_direction%ptr rename to veg%opt_leaf_cn_direction%ptr

     !IF (veg%labile_carbon_talloc_mavg%ptr/veg%labile_nitrogen_talloc_mavg%ptr > veg%growth_cn_talloc_mavg%ptr) THEN
     !   veg%leaf_cn_direction%ptr = -1._wp
     !ELSE
        veg%leaf_cn_direction%ptr = REAL(MAXLOC(npp_direction, 1),KIND=wp) - 2.0_wp
     !ENDIF

  ELSE
     veg%leaf_cn_direction%ptr = 0.0_wp
  END IF  ! IF(veg%pool%leaf%nitrogen%ptr > eps4)

  END SUBROUTINE calc_dir_optimal_cn_leaf
END MODULE mo_veg_canopy
