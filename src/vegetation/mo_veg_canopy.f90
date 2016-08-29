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

  USE mo_kind
  USE mo_jsb_lctlib_class,            ONLY: t_lctlib_element
  USE mo_exception,                   ONLY: message_text, finish
  USE mo_assimi_update_canopy_fluxes, ONLY: calc_photosynthesis

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: calc_canopy_layers,calc_npp_increment,optimize_Nleaf

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
  SUBROUTINE calc_canopy_layers(lctlib,leaf_nitrogen,lai,growing_season,par_sunlit_mavg_cl,par_shaded_mavg_cl, & 
                         fleaf_sunlit_mavg_cl,t_air_mavg,press_srf_mavg,co2_mixing_ratio_mavg,ga_mavg,& 
                         beta_air_mavg,beta_soil_ps_mavg,beta_soil_gs_mavg,t_jmax_opt_mavg,& 
                         lai_cl,cumm_lai_cl,n_tot_cl,fn_rub_cl,fn_et_cl,fn_pepc_cl,fn_chl_cl,fn_oth_cl,&
                         inquire_n_fractions)

    USE mo_base_constants,      ONLY: ncanop, canopy_layer_depth
    USE mo_assimi_constants
    USE mo_veg_constants
    USE mo_assimi_config_class

    IMPLICIT NONE
    ! 0.1 InOut
    TYPE(t_lctlib_element),INTENT(in)            :: lctlib                !< land-cover-type library
    REAL(wp),INTENT(in)                          :: leaf_nitrogen         !< total canopy N (mol/m2)
    REAL(wp),INTENT(in)                          :: lai                   !< total leaf area index
    LOGICAL,OPTIONAL                             :: growing_season        !< is plant in growing season?
    REAL(wp),INTENT(in),OPTIONAL                 :: t_air_mavg            !< air temperature (avg)
    REAL(wp),INTENT(in),OPTIONAL                 :: press_srf_mavg        !< air pressure (avg) 
    REAL(wp),INTENT(in),OPTIONAL                 :: co2_mixing_ratio_mavg !< CO2 concentration (avg)
    REAL(wp),INTENT(in),OPTIONAL                 :: ga_mavg               !< aerodynamic conductance (avg)
    REAL(wp),INTENT(in),OPTIONAL                 :: beta_air_mavg         !< beta air (avg)
    REAL(wp),INTENT(in),OPTIONAL                 :: beta_soil_ps_mavg     !< beta soil ps (avg)
    REAL(wp),INTENT(in),OPTIONAL                 :: beta_soil_gs_mavg     !< beta soil gs (avg)
    REAL(wp),INTENT(in),OPTIONAL                 :: t_jmax_opt_mavg       !< Jmax optimum temperature (avg)
    REAL(wp),DIMENSION(:),INTENT(inout),OPTIONAL :: par_sunlit_mavg_cl    !< absorbed PAR on sunlit
                                                                          !  leafes (average over time)
    REAL(wp),DIMENSION(:),INTENT(inout),OPTIONAL :: par_shaded_mavg_cl    !< absorbed PAR on shaded
                                                                          !  leafes (average over time)
    REAL(wp),DIMENSION(:),INTENT(inout),OPTIONAL :: fleaf_sunlit_mavg_cl  !< fraction of sunlit leaves
                                                                          !  (average over time)
    REAL(wp),DIMENSION(:),INTENT(OUT)            :: lai_cl                !< LAI of canopy layer
    REAL(wp),DIMENSION(:),INTENT(OUT)            :: cumm_lai_cl           !< cummlative LAI above canopy layer
    REAL(wp),DIMENSION(:),INTENT(inout)          :: n_tot_cl              !< total N in canopy layer (mmol N/m2)
    REAL(wp),DIMENSION(:),INTENT(inout)          :: fn_rub_cl             !< fraction of N in Rubisco
    REAL(wp),DIMENSION(:),INTENT(inout)          :: fn_et_cl              !< fraction of N in Electron Transport
    REAL(wp),DIMENSION(:),INTENT(inout)          :: fn_pepc_cl            !< fraction of N in PepC
    REAL(wp),DIMENSION(:),INTENT(inout)          :: fn_chl_cl             !< fraction of N in Chloroplast
    REAL(wp),DIMENSION(:),INTENT(inout)          :: fn_oth_cl             !< fraction of N not associated with PS
    LOGICAL,OPTIONAL                             :: inquire_n_fractions   !< whether to return fractions or 
                                                                          !  also run optimisation 
    ! 0.2 Local
    REAL(wp)  :: lc                             ! LAI at top of the canopy layer
    REAL(wp)  :: n_tot_tc                       ! top-of-the-canopy N content
    REAL(wp) :: k0_fnl_chl,k1_fnl_chl           ! factors in chlorophyll profile calculation
    REAL(wp) :: jmax2vcmax                      ! jmax to vcmax ratio
    REAL(wp) :: f_pepc                          ! fraction of N in PEP Carboxylase
    LOGICAL:: optimise_n_fractions
    INTEGER   :: icanop

    !> 0.3 some initialisation and input checks
    !!
    ! set pathway specific constants (BLARRP!! this should not be here -> lctlib)
    IF(lctlib%ps_pathway == ic4phot) THEN
        jmax2vcmax = jmax2vcmax_C4
        k0_fnl_chl = k0_fnl_chl_C4
        k1_fnl_chl = k1_fnl_chl_C4
        f_pepc     = k0_fnl_pepc
    ELSE
        jmax2vcmax = jmax2vcmax_C3
        k0_fnl_chl = k0_fnl_chl_C3
        k1_fnl_chl = k1_fnl_chl_C3
        f_pepc     = 0.0_wp
    ENDIF

    ! configure options for fractions optimisation 
    optimise_n_fractions = assimi_config%flag_optimal_Nfraction
    IF (PRESENT(inquire_n_fractions)) optimise_n_fractions = .NOT.inquire_n_fractions
    IF (PRESENT(growing_season)) optimise_n_fractions = optimise_n_fractions.AND.growing_season

    ! check if all input is available to run optimal fractions
    IF(assimi_config%flag_optimal_Nfraction)THEN
       IF(.NOT.PRESENT(par_sunlit_mavg_cl).OR..NOT.PRESENT(par_shaded_mavg_cl).OR..NOT.PRESENT(fleaf_sunlit_mavg_cl)) THEN
          WRITE(message_text,'(a)') 'called with fractional optimisation flag, but without par[sunlit|shaded]_mavg_cl'
          CALL finish("mo_veg_canopy:calc_canopy_layer",message_text)
       ENDIF
    ENDIF
    IF(optimise_n_fractions)THEN
       IF(.NOT.PRESENT(t_air_mavg).OR..NOT.PRESENT(press_srf_mavg).OR.& 
          .NOT.PRESENT(co2_mixing_ratio_mavg).OR..NOT.PRESENT(ga_mavg).OR.&
          .NOT.PRESENT(beta_soil_ps_mavg).OR..NOT.PRESENT(beta_soil_gs_mavg).OR.&
          .NOT.PRESENT(beta_air_mavg).OR..NOT.PRESENT(t_jmax_opt_mavg)) THEN

          WRITE(message_text,'(a)') 'called to optimise N fractions, but without [climate]_mavg'
          CALL finish("mo_veg_canopy:calc_canopy_layer",message_text)
       ENDIF
    ENDIF

    ! initialise output
    lai_cl(:) = 0._wp
    cumm_lai_cl(:) = 0._wp
    n_tot_cl(:) = 0._wp

    !> 1.0 calculate top of the canopy N
    !!
    IF(lai > canopy_layer_depth(1)) THEN
       n_tot_tc = leaf_nitrogen * 1000._wp * kn / (1._wp - exp(-kn * lai))
    ELSE ! if LAI does not fill the first canopy layer, top canopy N = leaf N
       n_tot_tc = leaf_nitrogen * 1000._wp
    ENDIF

    !> 2. calculate LAI and N content of canopy layers 
    !!
    IF(lai > eps8) THEN ! avoid any calculation if there are no leaves 
       lc = 0._wp
       DO icanop = 1,ncanop
          !> 2.1 LAI of canopy layer 
          !!
          IF(lc+canopy_layer_depth(icanop) < lai) THEN 
             ! if canopy layer is fully developed 
             lai_cl(icanop) = canopy_layer_depth(icanop)
          ELSE IF(lc < lai.AND.lc+canopy_layer_depth(icanop) >= lai) THEN 
             ! if canopy layer is only partially filled assign difference between full LAI and lc to canopy layer
             lai_cl(icanop) = lai - lc 
          ENDIF
          lc = lc + lai_cl(icanop)
          ! make sure all LAI is represented by letting last layer go to the end of the canopy
          IF(icanop==ncanop .AND. lc < lai ) lai_cl(ncanop) = lai_cl(ncanop) + lai - lc 

          !> 2.2 cummulative LAI above mid-point of canopy layer
          !!
          IF(lai_cl(icanop) > eps8)THEN
             IF(icanop==1) THEN
                cumm_lai_cl(icanop) = lai_cl(icanop) / 2.0_wp
             ELSE
                cumm_lai_cl(icanop) = cumm_lai_cl(icanop-1) + lai_cl(icanop-1) / 2.0_wp + lai_cl(icanop) / 2.0_wp
             ENDIF
          ENDIF

          !> 2.3 N content of canopy layer
          IF(cumm_lai_cl(icanop) > eps8) THEN
             n_tot_cl(icanop) = n_tot_tc * exp(-kn * cumm_lai_cl(icanop))
         ENDIF
       ENDDO
    ENDIF

    !> 3.0 fractionation of N within canopy laer
    !!
    !> 3.1 fraction not associated with photosynthesis 
    WHERE(n_tot_cl > 0.0_wp) 
       fn_oth_cl = lctlib%k0_fnl_struc - k1_fnl_struc * n_tot_cl   
    ELSEWHERE
       fn_oth_cl = 1.0_wp
    ENDWHERE


    !> 3.2 photosynthetic fractions
    !! prescribed, or dynamically calculated if assimi_config%flag_optimal_Nfraction = .TRUE.
    IF (.NOT.assimi_config%flag_optimal_Nfraction) THEN
        !> 3.2.1 case of presribed fractions
        WHERE(n_tot_cl > 0.0_wp) ! canopy layers with leaves
           ! fraction in chlorophyll 
           fn_chl_cl = (k0_fnl_chl - k1_fnl_chl * exp(-kfnl_chl * cumm_lai_cl)) / chl2n
           ! fraction in PepC
           fn_pepc_cl = f_pepc
           ! fraction of electron-transport is dependent on N in electron transport versus N in light harvesting
           fn_et_cl = (1._wp - (fn_oth_cl + fn_chl_cl + fn_pepc_cl )) / &
                      (1._wp + jmax2n/(jmax2vcmax * vcmax2n))
           ! fraction of Rubisco is determined as the remainder
           fn_rub_cl = 1._wp - fn_oth_cl - fn_chl_cl - fn_et_cl - fn_pepc_cl
        ELSEWHERE ! canopy layers without leaves
           n_tot_cl = 0.0_wp
           fn_oth_cl = 1.0_wp
           fn_chl_cl = 0.0_wp
           fn_et_cl =0.0_wp
           fn_rub_cl = 0.0_wp
           fn_pepc_cl = 0.0_wp
        ENDWHERE
     ELSE
        !> 3.2.1 optimal case 
        DO icanop = 1,ncanop
           IF(lai_cl(icanop) > 0.0_wp)THEN
              IF(fn_chl_cl(icanop) < eps8)THEN ! new layer with new leaves
                 IF(icanop == 1)THEN ! first layers: use static allocation as first guess
                    fn_chl_cl(icanop) = (k0_fnl_chl - k1_fnl_chl * exp(-kfnl_chl * lc)) / chl2n
                    ! fraction of electron-transport is dependent on N in electron transport versus N in light harvesting
                    fn_et_cl(icanop) = (1._wp - (fn_oth_cl(icanop) + fn_chl_cl(icanop) + fn_pepc_cl(icanop) )) / &
                                      (1._wp + jmax2n/(jmax2vcmax * vcmax2n))
                    ! fraction in PepC
                    fn_pepc_cl(icanop)=f_pepc
                    ! fraction of Rubisco is determined as the remainder
                    fn_rub_cl(icanop) = 1._wp - fn_oth_cl(icanop) - fn_chl_cl(icanop) - fn_et_cl(icanop) - fn_pepc_cl(icanop)
                 ELSE ! other layers: copy from above layer
                    fn_chl_cl(icanop)  = fn_chl_cl(icanop-1)
                    fn_et_cl(icanop)   = fn_et_cl(icanop-1)
                    fn_rub_cl(icanop)  = fn_rub_cl(icanop-1)
                    fn_pepc_cl(icanop) = fn_pepc_cl(icanop-1)
                    par_sunlit_mavg_cl(icanop) = par_sunlit_mavg_cl(icanop-1)
                    par_shaded_mavg_cl(icanop) = par_shaded_mavg_cl(icanop-1) 
                    fleaf_sunlit_mavg_cl(icanop) = fleaf_sunlit_mavg_cl(icanop-1)
                 ENDIF
              ENDIF
              IF(optimise_n_fractions)THEN

                 CALL optimise_canopy_n_fractions(lctlib, &
                            n_tot_cl(icanop), &
                            par_sunlit_mavg_cl(icanop),&
                            par_shaded_mavg_cl(icanop), &
                            fleaf_sunlit_mavg_cl(icanop),&
                            t_air_mavg, &
                            press_srf_mavg,&
                            co2_mixing_ratio_mavg,&
                            ga_mavg,&
                            beta_air_mavg,&
                            beta_soil_ps_mavg,&
                            beta_soil_gs_mavg,&
                            t_jmax_opt_mavg,&
                            fn_oth_cl(icanop),&
                            fn_chl_cl(icanop),&
                            fn_rub_cl(icanop),&
                            fn_et_cl(icanop),&
                            fn_pepc_cl(icanop))
   
              ENDIF
           ELSE ! no leaves in canopy layer
              n_tot_cl(icanop)   = 0.0_wp
              fn_oth_cl(icanop)  = 1.0_wp
              fn_chl_cl(icanop)  = 0.0_wp
              fn_et_cl(icanop)   = 0.0_wp
              fn_rub_cl(icanop)  = 0.0_wp
              fn_pepc_cl(icanop) = 0.0_wp
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
  ELEMENTAL SUBROUTINE optimise_canopy_n_fractions(lctlib,n_tot_cl,par_sunlit_mavg_cl,par_shaded_mavg_cl,&
                                 fleaf_sunlit_mavg_cl,t_air_mavg,press_srf_mavg,co2_mixing_ratio_mavg,ga_mavg,& 
                                 beta_air_mavg,beta_soil_ps_mavg,beta_soil_gs_mavg,t_jmax_opt_mavg,&
                                 fn_oth_cl,fn_chl_cl,fn_rub_cl,fn_et_cl,fn_pepc_cl)

    USE mo_base_constants, ONLY: dt,one_day
    USE mo_veg_constants,  ONLY: delta_n_fraction

    IMPLICIT NONE
    ! 0.1 InOut
    TYPE(t_lctlib_element),INTENT(in) :: lctlib         !< land-cover-type library - parameters for one pft
    REAL(wp),              INTENT(in) :: n_tot_cl, &    !< total canopy N content
      & par_sunlit_mavg_cl,    &    !< PAR on sunlit leaves (day-time running mean, µmol/ms/s) 
      & par_shaded_mavg_cl,    &    !< PAR on shaded leaves (day-time running mean, µmol/ms/s) 
      & fleaf_sunlit_mavg_cl,  &    !< fraction of sunlit leaves (day-time running mean, µmol/ms/s)
      & t_air_mavg,            &    !< temperature (day-time running mean, K)
      & press_srf_mavg,        &    !< air pressure (day-time running mean, Pa)
      & co2_mixing_ratio_mavg, &    !< CO2 mixing ratio (day-time running mean, ppm)
      & ga_mavg,               &    !< aerodynamic conductance (day-time running mean, m/s)
      & beta_air_mavg,         &    !< beta air (day-time running mean, --)
      & beta_soil_ps_mavg,     &    !< beta soil for PS (day-time running mean, --)
      & beta_soil_gs_mavg,     &    !< beta soil for GC (day-time running mean, --) 
      & t_jmax_opt_mavg             !< optimum temperature for Jmax (day-time running mean, °C)
    REAL(wp),              INTENT(in) :: fn_oth_cl      !< fraction of leaf N not associated with PS
    REAL(wp),              INTENT(inout) :: & 
       & fn_chl_cl,            &    !< fraction of leaf N in chlorophyll
       & fn_rub_cl,            &    !< fraction of leaf N in rubisco
       & fn_et_cl,             &    !< fraction of leaf N in electron transport
       & fn_pepc_cl                 !< fraction of leaf N in electron transport
    ! 0.2 Local
    REAl(wp) :: ag0_cl, an0_cl, r0_cl, ci0_cl,gs0_cl           ! photosynthesis outputs of initial guess
    REAL(wp) :: fn_chl0_cl, fn_rub0_cl, fn_et0_cl, fn_pepc0_cl ! N fractions of initial guess
    REAL(wp) :: ag_cl, an_cl, r_cl, ci_cl,  gs_cl              ! photosynthesis outputs of optimisation runs
    REAL(wp) :: fn_rub_opt_cl, fn_et_opt_cl, fn_pepc_opt_cl    ! optimal N fractions
    REAL(wp) :: m_rub, m_et, m_pepc                            ! specific PS rates
    INTEGER :: i
    REAL(wp):: hlp
 
    !> 1.0 Set N fractions to previous values and calculate photosynthesis for avererage conditions 
    ! set fractions to previous values. In current code chl + et+ rub=1
    ! all other routines chl +et +rub +oth = 1
    IF(fn_oth_cl > eps8) THEN
       fn_chl0_cl = fn_chl_cl/(1._wp-fn_oth_cl)
       fn_rub0_cl = fn_rub_cl/(1._wp-fn_oth_cl)
       fn_et0_cl = fn_et_cl/(1._wp-fn_oth_cl)
       fn_pepc0_cl = fn_pepc_cl/(1._wp-fn_oth_cl)
    ENDIF
    fn_chl_cl = fn_chl0_cl
    fn_et_cl = fn_et0_cl
    fn_rub_cl = fn_rub0_cl
    fn_pepc_cl = fn_pepc0_cl

    ! calculate photosynthesis with previous fractions
    CALL calc_photosynthesis(t_air_mavg,press_srf_mavg,co2_mixing_ratio_mavg,ga_mavg, & 
              par_sunlit_mavg_cl,par_shaded_mavg_cl,fleaf_sunlit_mavg_cl, & 
              beta_air_mavg,beta_soil_ps_mavg,beta_soil_gs_mavg, &
              fn_chl0_cl*(1._wp-fn_oth_cl),fn_et0_cl*(1._wp-fn_oth_cl),fn_rub0_cl*(1._wp-fn_oth_cl), & 
              fn_pepc0_cl*(1._wp-fn_oth_cl),n_tot_cl,lctlib%ps_pathway, & 
              t_jmax_opt_mavg,lctlib%t_jmax_omega,lctlib%g0,lctlib%g1, &
              ag0_cl,an0_cl,r0_cl,gs0_cl,ci0_cl,m_rub,m_et,m_pepc)
   
    !> 2.0 Attempt to increase PS by increasing the N fraction in chlorophyll 
    !! 
    fn_chl_cl = MIN(fn_chl0_cl * (1.0_wp + delta_n_fraction*dt/one_day), 1.0_wp)
    
    ! optimal N fraction assuming Aj=Av=Apepc
    CALL calc_coordination_n_fractions(fn_chl_cl,m_rub,m_et,m_pepc, & 
                                       fn_rub_opt_cl,fn_et_opt_cl,fn_pepc_opt_cl)
    ! check if changing fraction to optimal values is possible in one timestep
    CALL calc_update_n_fractions(delta_n_fraction*dt/one_day,fn_chl_cl,& 
                                  fn_rub_opt_cl,fn_et_opt_cl,fn_pepc_opt_cl, & 
                                  fn_rub_cl,fn_et_cl,fn_pepc_cl) 
                   
    ! re-calculate photosynthesis
    CALL calc_photosynthesis(t_air_mavg,press_srf_mavg,co2_mixing_ratio_mavg,ga_mavg, & 
              par_sunlit_mavg_cl,par_shaded_mavg_cl,fleaf_sunlit_mavg_cl, & 
              beta_air_mavg,beta_soil_ps_mavg,beta_soil_gs_mavg, &
              fn_chl_cl*(1._wp-fn_oth_cl),fn_et_cl*(1._wp-fn_oth_cl),fn_rub_cl*(1._wp-fn_oth_cl), & 
              fn_pepc_cl*(1._wp-fn_oth_cl),n_tot_cl,lctlib%ps_pathway, & 
              t_jmax_opt_mavg,lctlib%t_jmax_omega,lctlib%g0,lctlib%g1, &
              ag_cl,an_cl,r_cl,gs_cl,ci_cl,hlp,hlp,hlp)
    
    !> 3.0 check if increased chlorophyll has increased photosynthesis. If not, try decreasing it
    IF (an_cl < an0_cl) THEN
       
       ! decrease N fraction in chlorophyll
       fn_chl_cl = MAX(fn_chl0_cl * (1._wp - delta_n_fraction*dt/one_day), eps8)
       fn_et_cl = fn_et0_cl
       fn_rub_cl = fn_rub0_cl
       fn_pepc_cl = fn_pepc0_cl

       ! optimal N fraction assuming Aj=Av=Apepc
       CALL calc_coordination_n_fractions(fn_chl_cl,m_rub,m_et,m_pepc, & 
                                          fn_rub_opt_cl,fn_et_opt_cl,fn_pepc_opt_cl)
       ! check if changing fraction to optimal values is possible in one timestep
       CALL calc_update_n_fractions(delta_n_fraction*dt/one_day,fn_chl_cl,& 
                                     fn_rub_opt_cl,fn_et_opt_cl,fn_pepc_opt_cl, & 
                                     fn_rub_cl,fn_et_cl,fn_pepc_cl) 
     
       ! re-calculate photosynthesis
       CALL calc_photosynthesis(t_air_mavg,press_srf_mavg,co2_mixing_ratio_mavg,ga_mavg, & 
            par_sunlit_mavg_cl,par_shaded_mavg_cl,fleaf_sunlit_mavg_cl,& 
            beta_air_mavg,beta_soil_ps_mavg,beta_soil_gs_mavg, &
            fn_chl_cl*(1._wp-fn_oth_cl), fn_et_cl*(1._wp-fn_oth_cl), fn_rub_cl*(1._wp-fn_oth_cl), & 
            fn_pepc_cl*(1._wp-fn_oth_cl),n_tot_cl, lctlib%ps_pathway, & 
            t_jmax_opt_mavg,lctlib%t_jmax_omega,lctlib%g0,lctlib%g1, &
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
          CALL calc_update_n_fractions(delta_n_fraction*dt/one_day,fn_chl_cl,& 
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
    REAL(wp), INTENT(in) :: fn_chl_cl, &       !< fraction of N in chlorophyll
     &     m_rub, &                            !< N-specific rate constant for Rubisco
     &     m_et,  &                            !< N-specific rate constant for ET
     &     m_pepc                              !< N-speficic rate constant for PepC (C4-only)
    REAL(wp), INTENT(out) :: fn_rub_opt_cl, &  !< optimal fraction of N in Rubisco 
     &     fn_et_opt_cl, &                     !< optimal fraction of N in ET 
     &     fn_pepc_opt_cl                      !< optimal fraction of N in PepC (C4-only)

    ! optimal N fraction assuming Aj=Av=Apepc
!@todo: veg_canopy:optimise_n_fractions: should this not replaced by lctlib%ps_pathway==ic4phot?
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
  ! Sub Task to update_canopy_fluxes
  ! 
  !------------------------------------------------------------------------------------------------------
  !> Subroutine to calculate change in NPP with very small change in LAI
  !! Needed for calculation of optimal leaf and root allocation fraction. See mo_veg_growth
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE calc_npp_increment(plant, lctlib, dNPP)

    USE mo_base_memory
    USE mo_base_config_class,       ONLY: base_config
    USE mo_base_constants
    USE mo_shared_constants        ! only lambda_C14 ?
    USE mo_veg_memory              ! only t_slm_plant ?
   ! USE mo_veg_growth,          ONLY: calc_allometry
   !USE mo_rad_update_radiation, ONLY: calc_canopy_layer_depth

   IMPLICIT NONE
   ! 0.1 InOut
   TYPE(t_slm_plant), INTENT(in)          :: plant
   TYPE(t_lctlib_element),INTENT(in)      :: lctlib         !< land-cover-type library - parameters for one pft
   REAL(wp), INTENT(out)                  :: dNPP
   ! 0.2 Local
   REAL(wp) :: lai, nitrogen
   REAL(wp), DIMENSION(ncanop) :: canopy_layer_depth_act
   REAL(wp), DIMENSION(ncanop) :: lai_cl
   REAL(wp), DIMENSION(ncanop) :: cumm_lai_cl 
   REAL(wp), DIMENSION(ncanop) :: par_sunlit_mavg_cl,par_shaded_mavg_cl,fleaf_sunlit_mavg_cl 
   REAl(wp), DIMENSION(ncanop) :: fnl_oth_cl, nl_tot_cl, fnl_chl_cl, fnl_et_cl, fnl_rub_cl, fnl_pepc_cl
   REAL(wp) :: gpp, maint_respiration_leaf
   REAL(wp), DIMENSION(ncanop) :: agl, anl, rl, gsl, cil
   REAL(wp) :: hlp1 
   INTEGER :: icanop

   fnl_chl_cl = plant%fnl_chl_cl(:)
   fnl_et_cl = plant%fnl_et_cl(:)
   fnl_rub_cl = plant%fnl_rub_cl(:)
   fnl_pepc_cl = plant%fnl_pepc_cl(:)
   par_sunlit_mavg_cl = plant%moving_average%par_sunlit_cl
   par_shaded_mavg_cl = plant%moving_average%par_shaded_cl
   fleaf_sunlit_mavg_cl = plant%moving_average%fleaf_sunlit_cl

   ! increase LAI by small increment
   lai = plant%lai + 0.3_wp * lctlib%sla
   ! equivalent leaf N assuming any new leaf growth will have target CN
   nitrogen = plant%pool%leaf%nitrogen + 0.3_wp * plant%target_cn_leaf
   ! recompute LAI and foliar N with new numbers
   CALL calc_canopy_layers(lctlib = lctlib, & 
                          leaf_nitrogen = nitrogen, &
                          lai = lai,&
                          par_sunlit_mavg_cl = par_sunlit_mavg_cl, &
                          par_shaded_mavg_cl = par_shaded_mavg_cl, &
                          fleaf_sunlit_mavg_cl = fleaf_sunlit_mavg_cl, &
                          lai_cl = lai_cl, &
                          cumm_lai_cl = cumm_lai_cl, &
                          n_tot_cl = nl_tot_cl, &
                          fn_rub_cl = fnl_rub_cl, &
                          fn_et_cl = fnl_et_cl, &
                          fn_pepc_cl = fnl_pepc_cl, &
                          fn_chl_cl = fnl_chl_cl, &
                          fn_oth_cl = fnl_oth_cl, &
                          inquire_n_fractions=.TRUE.)
!   ! recalculate number of canopy layers with new LAI
!   canopy_layer_depth_act = calc_canopy_layer_depth(lai)
!
!   ! re-calculate N distribution in the canopy, needed because of new LAI and leaf nitrogen
!   ! photosynthetic N fractions should NOT be recalculated to allow interaction with optimal fractions
!    CALL calc_canopy_n_distribution(nitrogen,lai, & 
!            canopy_layer_depth_act, lctlib%ps_pathway, & 
!            lctlib%k0_fnl_struc, nl_tot_cl, fnl_chl_cl, fnl_et_cl, fnl_rub_cl, fnl_pepc_cl, fnl_oth_cl)    


   ! Loop over canopy layers
   ! calculate photosynthesis per layer
   DO icanop = 1,ncanop
     
        ! only calculate phosynthesis if leaves present in layer
        !IF (canopy_layer_depth_act(icanop).GT.eps8)THEN   
        IF (plant%lai_cl(icanop).GT.eps8)THEN   
           
           CALL calc_photosynthesis(plant%moving_average%t_air,plant%moving_average%press_srf,& 
                plant%moving_average%co2_mixing_ratio,plant%moving_average%ga, & 
                plant%par_sunlit_cl(icanop), plant%par_shaded_cl(icanop),& 
                plant%fleaf_sunlit_cl(icanop), &
                plant%beta_air,plant%beta_soil_ps,plant%beta_soil_gs, &
                fnl_chl_cl(icanop),fnl_et_cl(icanop),fnl_rub_cl(icanop),fnl_pepc_cl(icanop), nl_tot_cl(icanop),& 
                lctlib%ps_pathway, &
                plant%moving_average%t_jmax_opt, &                  
                lctlib%t_jmax_omega,lctlib%g0,lctlib%g1, &
                agl(icanop),anl(icanop),rl(icanop),gsl(icanop),cil(icanop),&
                hlp1,hlp1,hlp1)
        END IF
        
     END DO ! icanop = 1,ncanop
     
     ! integrate fluxes over the canopy
     gpp = 0._wp
     maint_respiration_leaf = 0._wp
     DO icanop = 1,ncanop
        !gpp = gpp + agl(icanop) * canopy_layer_depth_act(icanop)  
        gpp = gpp + agl(icanop) * plant%lai_cl(icanop)  
        maint_respiration_leaf = maint_respiration_leaf + & 
            !rl(icanop) * canopy_layer_depth_act(icanop) 
            rl(icanop) * plant%lai_cl(icanop) 
     ENDDO

     ! increment in NPP
     dNPP = (gpp - plant%flux%gpp) - (maint_respiration_leaf - plant%flux%maint_respiration_leaf)

  END SUBROUTINE calc_npp_increment
  !> Subroutine to determine direction in which total leaf N should change in order to increase NPP
  !! 
  !! Parameter leaf_cn_direction
  !----------------------------------------------------------------------------------------------------
  SUBROUTINE optimize_Nleaf(plant, canopy_layer_depth_act,lctlib, t_air, press_srf,co2_mixing_ratio,ga,beta_air,beta_soil_ps, & 
                          beta_soil_gs)

    USE mo_base_constants    
    USE mo_veg_respiration
    USE mo_veg_memory, ONLY: t_slm_plant, t_slm_veg_pool
    USE mo_veg_constants, ONLY: delta_n_leaf
    
    IMPLICIT NONE
    ! 0.1 InOut
    TYPE(t_slm_plant), INTENT(inout)        :: plant                  ! plant variables per tile
    TYPE(t_lctlib_element),INTENT(in)        :: lctlib         !< land-cover-type library - parameters for one pft
    REAL(wp), DIMENSION(ncanop), INTENT(in) :: canopy_layer_depth_act
    REAL(wp) :: t_air, press_Srf,co2_mixing_ratio,ga
    REAL(wp) ::  beta_air,beta_soil_ps,beta_soil_gs
    ! 0.2 Local
    INTEGER :: icanop, dd, t
    TYPE(t_slm_veg_pool) :: new
    REAL(wp),DIMENSION(ncanop, 3) :: agl,anl,rl,gsl,cil
    REAL(wp) :: m_rub, m_et, m_pepc
   REAL(wp), DIMENSION(ncanop) :: par_sunlit_mavg_cl,par_shaded_mavg_cl,fleaf_sunlit_mavg_cl 
    REAL(wp), DIMENSION(3) :: npp_direction, gpp_direction, rl_direction
    REAL(wp), DIMENSION(ncanop) :: nl_tot_cl, fnl_oth_cl, fnl_chl_cl, fnl_et_cl, fnl_rub_cl, fnl_pepc_cl,lai_cl,cumm_lai_cl
    REAL(wp), DIMENSION(3) :: maint_respiration
    REAL(Wp) :: maint_respiration_c13,maint_respiration_c14

    fnl_chl_cl = plant%fnl_chl_cl(:)
    fnl_et_cl = plant%fnl_et_cl(:)
    fnl_rub_cl = plant%fnl_rub_cl(:)
   par_sunlit_mavg_cl = plant%moving_average%par_sunlit_cl
   par_shaded_mavg_cl = plant%moving_average%par_shaded_cl
   fleaf_sunlit_mavg_cl = plant%moving_average%fleaf_sunlit_cl
    ! initialse output
    agl(:, :) = 0.0_wp 
    anl(:, :) = 0.0_wp 
    rl(:, :) = 0.0_wp
    gsl(:, :) = 0.0_wp
    cil(:, :) = 0.0_wp 
    ! Calculate photosynthesis with lower, current and higher leaf N to determine optimal canopy N
    DO dd = -1, 1
        
       new%leaf%nitrogen = plant%pool%leaf%nitrogen * (1.0_wp + REAL(dd) * delta_n_leaf)
       
       ! update N distribution in the canopy
       CALL calc_canopy_layers(lctlib = lctlib, & 
                          leaf_nitrogen = new%leaf%nitrogen, &
                          lai = plant%lai,&
                          par_sunlit_mavg_cl = par_sunlit_mavg_cl, &
                          par_shaded_mavg_cl = par_shaded_mavg_cl, &
                          fleaf_sunlit_mavg_cl = fleaf_sunlit_mavg_cl, &
                          lai_cl = lai_cl, &
                          cumm_lai_cl = cumm_lai_cl, &
                          n_tot_cl = nl_tot_cl, &
                          fn_rub_cl = fnl_rub_cl, &
                          fn_et_cl = fnl_et_cl, &
                          fn_pepc_cl = fnl_pepc_cl, &
                          fn_chl_cl = fnl_chl_cl, &
                          fn_oth_cl = fnl_oth_cl, &
                          inquire_n_fractions=.TRUE.)

        ! calculate photosynthesis for each layer

       DO icanop =1, ncanop
          IF(canopy_layer_depth_act(icanop).GT.0.0)THEN
             CALL calc_photosynthesis(t_air,press_srf,co2_mixing_ratio,ga, & 
                  plant%moving_average%par_sunlit_cl(icanop),plant%moving_average%par_shaded_cl(icanop),& 
                  plant%moving_average%fleaf_sunlit_cl(icanop), &
                  beta_air,beta_soil_ps,beta_soil_gs, &
                  plant%fnl_chl_cl(icanop),&
                  plant%fnl_et_cl(icanop),& 
                  plant%fnl_rub_cl(icanop),&
                  plant%fnl_pepc_cl(icanop),&
                  nl_tot_cl(icanop),&
                  lctlib%ps_pathway, & 
                  plant%moving_average%t_jmax_opt, &                  
                  lctlib%t_jmax_omega,lctlib%g0,lctlib%g1, &
                  agl(icanop, dd+2),anl(icanop, dd+2),rl(icanop, dd+2),gsl(icanop, dd+2),cil(icanop, dd+2),&
                  m_rub, m_et, m_pepc)
          
                 
          ENDIF
       ENDDO
    ENDDO

       ! calculate NPP for each leaf N value (lower, equal, higher)
       gpp_direction(:) = 0.0_wp
       rl_direction(:)=0.0_wp
       DO t=1, 3
          DO icanop=1, ncanop
             gpp_direction(t) = gpp_direction(t) + agl(icanop, t)
             rl_direction(t) = rl_direction(t) + rl(icanop,t)
          ENDDO
       ENDDO


       DO dd=-1, 1
       ! Calculate nitrogen in all other pools assuming lower leaf nitrogen
       ! OBS: currently assuming N is atken from fine root pool since respiration 
       ! fucntion is linear with sum of pools. Should be changed if adding stem and root temperature
       new%fine_root%nitrogen = plant%pool%fine_root%nitrogen * (1 - REAL(dd) * delta_n_leaf)
       new%coarse_root%nitrogen = plant%pool%coarse_root%nitrogen
       new%sap_wood%nitrogen = plant%pool%sap_wood%nitrogen
        !Calculate whole plant respiration
 
        maint_respiration(dd+2) = calc_maintenance_respiration(t_air,t_air,new%sap_wood%nitrogen,&
             new%fine_root%nitrogen, new%coarse_root%nitrogen) + rl_direction(dd+2)
        ! calculate NPP
          npp_direction(dd+2) = gpp_direction(dd+2)-maint_respiration(dd+2)
     ENDDO
        
       ! Direction of change in leaf N. Actual leaf N is modified in interfacte routine
       plant%leaf_cn_direction = REAL(MAXLOC(npp_direction, 1)) - 2.0_wp

  END SUBROUTINE optimize_Nleaf

END MODULE mo_veg_canopy
