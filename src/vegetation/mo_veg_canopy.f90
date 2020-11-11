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
#ifndef __NO_JSBACH__

  USE mo_kind,                        ONLY: wp
  USE mo_jsb_impl_constants,          ONLY: test_false_true
  USE mo_exception,                   ONLY: message_text, finish
  !USE mo_assimi_update_canopy_fluxes, ONLY: calc_photosynthesis

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: calc_canopy_layers !, calc_marginal_canopy_flux_increment, calc_dir_optimal_cn_leaf

  CHARACTER(len=*), PARAMETER :: modname = 'mo_veg_canopy'

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
  SUBROUTINE calc_canopy_layers(nc                    , &
                                ncanopy               , &
                                dz                    , &     
                                lbounds               , &
                                ubounds               , &
                                lctlib_ps_pathway     , &
                                lctlib_k0_fn_struc    , &
                                lctlib_fn_oth_min     , &
                                lctlib_sla            , &
                                lctlib_np_leaf        , &
                                lctlib_gmin           , &
                                lctlib_g0             , &
                                lctlib_g1             , &
                                lctlib_t_jmax_omega   , &
                                flag_optimal_Nfraction, &
                                leaf_nitrogen     , &
                                lai               , &
                                ppfd_sunlit_cl    , &
                                ppfd_shaded_cl    , &
                                fleaf_sunlit_cl   , &
                                t_air             , &
                                t_acclim          , &
                                press_srf         , &
                                co2_mixing_ratio  , &
                                ga                , &
                                beta_air          , &
                                beta_soa          , &
                                beta_soil_ps      , &
                                beta_sinklim_ps   , &
                                beta_soil_gs      , &
                                t_jmax_opt        , &
                                lai_cl            , &
                                cumm_lai_cl       , &
                                leaf_nitrogen_cl  , &
                                fn_rub_cl         , &
                                fn_et_cl          , &
                                fn_pepc_cl        , &
                                fn_chl_cl         , &
                                fn_oth_cl         , &
                                inquire_n_fractions)

    USE mo_veg_constants,           ONLY: carbon_per_dryweight_leaf, &
                                          jmax2vcmax_C4, k0_fn_chl_C4, k1_fn_chl_C4, k0_fn_pepc_C4, &
                                          jmax2vcmax_C3, k0_fn_chl_C3, k1_fn_chl_C3, k0_fn_pepc_C3, kfn_chl
    USE mo_assimi_constants,        ONLY: kn, chl2n, ic4phot, jmax2n, k1_fn_struc, vcmax2n
    USE mo_jsb_physical_constants,  ONLY: molar_mass_C, molar_mass_N, molar_mass_P
    USE mo_jsb_math_constants,      ONLY: eps8, eps12

    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    INTEGER,                                INTENT(in)              :: nc                    !< dimensions
    INTEGER,                                INTENT(in)              :: ncanopy               !< number of canopy layers
    REAL(wp),DIMENSION(ncanopy),            INTENT(in)              :: dz, &                 !< canopy layer thickness  [m2 m-2 LAI]
                                                                       lbounds, &            !< lower bound of canopy layer (lbounds < ubounds; lbounds is at top of canopy layer)
                                                                       ubounds               !< upper bound of canopy layer (lbounds < ubounds; ubounds is at bottom of canopy layer)
    INTEGER,                                INTENT(in)              :: lctlib_ps_pathway     !< dimensions
    REAL(wp),                               INTENT(in)              :: lctlib_k0_fn_struc    !< lctlib parameter
    REAL(wp),                               INTENT(in)              :: lctlib_fn_oth_min     !< lctlib parameter
    REAL(wp),                               INTENT(in)              :: lctlib_sla            !< lctlib parameter
    REAL(wp),                               INTENT(in)              :: lctlib_np_leaf        !< lctlib parameter
    REAL(wp),                               INTENT(in)              :: lctlib_gmin           !< lctlib parameter
    REAL(wp),                               INTENT(in)              :: lctlib_g0             !< lctlib parameter
    REAL(wp),                               INTENT(in)              :: lctlib_g1             !< lctlib parameter
    REAL(wp),                               INTENT(in)              :: lctlib_t_jmax_omega   !< lctlib parameter  
    LOGICAL,                                INTENT(in)              :: flag_optimal_Nfraction!< on/off optimise leaf internal N allocation
    REAL(wp),DIMENSION(nc),                 INTENT(in)              :: leaf_nitrogen         !< total canopy N (mol/m2)
    REAL(wp),DIMENSION(nc),                 INTENT(in)              :: lai                   !< total leaf area index
    REAL(wp),DIMENSION(nc),                 INTENT(in),   OPTIONAL  :: t_air                 !< air temperature (avg)
    REAL(wp),DIMENSION(nc),                 INTENT(in),   OPTIONAL  :: t_acclim              !< acclimation temperature (avg)
    REAL(wp),DIMENSION(nc),                 INTENT(in),   OPTIONAL  :: press_srf             !< air pressure (avg)
    REAL(wp),DIMENSION(nc),                 INTENT(in),   OPTIONAL  :: co2_mixing_ratio      !< CO2 concentration (avg)
    REAL(wp),DIMENSION(nc),                 INTENT(in),   OPTIONAL  :: ga                    !< aerodynamic conductance (avg)
    REAL(wp),DIMENSION(nc),                 INTENT(in),   OPTIONAL  :: beta_air              !< beta air (avg)
    REAL(wp),DIMENSION(nc),                 INTENT(in),   OPTIONAL  :: beta_soa              !< beta soa (avg)
    REAL(wp),DIMENSION(nc),                 INTENT(in),   OPTIONAL  :: beta_soil_ps          !< beta soil ps (avg)
    REAL(wp),DIMENSION(nc),                 INTENT(in),   OPTIONAL  :: beta_sinklim_ps       !< beta sinklim ps (avg)
    REAL(wp),DIMENSION(nc),                 INTENT(in),   OPTIONAL  :: beta_soil_gs          !< beta soil gs (avg)
    REAL(wp),DIMENSION(nc),                 INTENT(in),   OPTIONAL  :: t_jmax_opt            !< Jmax optimum temperature (avg)
    REAL(wp),DIMENSION(nc,ncanopy),         INTENT(inout)           :: ppfd_sunlit_cl        !< absorbed Photosynthetically Active Photon Flux Density on sunlit leafes (average over time)
    REAL(wp),DIMENSION(nc,ncanopy),         INTENT(inout)           :: ppfd_shaded_cl        !< absorbed Photosynthetically Active Photon Flux Density on shaded leafes (average over time)
    REAL(wp),DIMENSION(nc,ncanopy),         INTENT(inout)           :: fleaf_sunlit_cl       !< fraction of sunlit leaves (average over time)
    REAL(wp),DIMENSION(nc,ncanopy),         INTENT(out)             :: lai_cl                !< LAI of canopy layer
    REAL(wp),DIMENSION(nc,ncanopy),         INTENT(out)             :: cumm_lai_cl           !< cummlative LAI above canopy layer
    REAL(wp),DIMENSION(nc,ncanopy),         INTENT(out)             :: leaf_nitrogen_cl      !< total N in canopy layer (mmol N/m2)
    REAL(wp),DIMENSION(nc,ncanopy),         INTENT(inout)           :: fn_rub_cl             !< fraction of N in Rubisco
    REAL(wp),DIMENSION(nc,ncanopy),         INTENT(inout)           :: fn_et_cl              !< fraction of N in Electron Transport
    REAL(wp),DIMENSION(nc,ncanopy),         INTENT(inout)           :: fn_pepc_cl            !< fraction of N in PepC
    REAL(wp),DIMENSION(nc,ncanopy),         INTENT(inout)           :: fn_chl_cl             !< fraction of N in Chloroplast
    REAL(wp),DIMENSION(nc,ncanopy),         INTENT(inout)           :: fn_oth_cl             !< fraction of N not associated with PS
    LOGICAL,                                INTENT(in),   OPTIONAL  :: inquire_n_fractions   !< whether to return fractions or also run optimisation
    ! ---------------------------
    ! 0.2 Local
    REAL(wp),DIMENSION(nc)          :: leaf_nitrogen_tc                              ! top-of-the-canopy N content
    REAL(wp)                        :: k0_fn_chl, k1_fn_chl                          ! factors in chlorophyll profile calculation
    REAL(wp)                        :: jmax2vcmax                                    ! jmax to vcmax ratio
    REAL(wp)                        :: f_pepc                                        ! fraction of N in PEP Carboxylase
    LOGICAL                         :: optimise_n_fractions
    INTEGER                         :: icanopy, ipoint                               ! looping over ncanopy and nc
    REAL(wp),DIMENSION(nc,ncanopy)  :: hlp1,hlp2,hlp3,hlp4
    REAL(wp),DIMENSION(nc,ncanopy)  :: mid_point_lai_cl
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_canopy_layers'
    ! ---------------------------
    ! 0.3 Memory allocation
    ! ...
    ! ---------------------------
    !> 0.4 some initialisation and input checks
    !!



    ! ------------------------------------------------------------------------------------------------------------
    ! Go science
    ! ------------------------------------------------------------------------------------------------------------
    

    !> 0.9 initialise output (to ensure that 0.0 is returned in case of the below statements is TRUE)
    !! 
    lai_cl(:,:)           = 0._wp
    cumm_lai_cl(:,:)      = 0._wp
    leaf_nitrogen_cl(:,:) = 0._wp
    mid_point_lai_cl(:,:) = 0._wp

    !> 0.9.1  set pathway specific constants 
    !! 
    IF(lctlib_ps_pathway == ic4phot) THEN
        jmax2vcmax = jmax2vcmax_C4
        k0_fn_chl  = k0_fn_chl_C4
        k1_fn_chl  = k1_fn_chl_C4
        f_pepc     = k0_fn_pepc_C4
    ELSE
        jmax2vcmax = jmax2vcmax_C3
        k0_fn_chl  = k0_fn_chl_C3
        k1_fn_chl  = k1_fn_chl_C3
        f_pepc     = k0_fn_pepc_C3
    ENDIF

    ! !> 0.9.2 configure options for fractions optimisation
    ! !! 
    ! optimise_n_fractions = flag_optimal_Nfraction
    ! IF (PRESENT(inquire_n_fractions)) optimise_n_fractions = .NOT.inquire_n_fractions
    ! 
    ! !> 0.9.3  check presence of all arguments that are OPTIONAL but needed for 'CALL optimise_canopy_n_fractions()' \n
    ! !!         only in case optimise_n_fractions = TRUE \n
    ! !!         note, all variables are of DIMENSION(nc)         but compiler gives an error if '(:)' is added here
    ! IF(optimise_n_fractions)THEN
    !    IF(.NOT.PRESENT(t_air)            .OR..NOT.PRESENT(t_acclim)         .OR..NOT.PRESENT(press_srf)    .OR. &
    !       .NOT.PRESENT(co2_mixing_ratio) .OR..NOT.PRESENT(ga)               .OR. &
    !       .NOT.PRESENT(beta_air)         .OR..NOT.PRESENT(beta_soa)         .OR. &
    !       .NOT.PRESENT(beta_soil_ps)     .OR..NOT.PRESENT(beta_sinklim_ps)  .OR..NOT.PRESENT(beta_soil_gs) .OR. &
    !       .NOT.PRESENT(t_jmax_opt)) THEN
    ! 
    !       WRITE(message_text,'(a)') 'called to optimise N fractions, but without [climate] variables'
    !       CALL finish("mo_veg_canopy:calc_canopy_layer",message_text)
    !    ENDIF
    ! ENDIF


    !> 1.0 calculate top of the canopy N
    !!
    !! dz(1) == ubounds(1) -> thickness of uppermost canopy layer
    WHERE(lai(:) > dz(1))
       leaf_nitrogen_tc(:) = leaf_nitrogen(:) * 1000._wp * kn / (1._wp - exp(-kn * lai(:)))
    ELSEWHERE(lai(:) > eps8)  ! if LAI does not fill the first canopy layer, top canopy N = leaf N
       leaf_nitrogen_tc(:) = leaf_nitrogen(:) * 1000._wp /lai(:)
    ENDWHERE


    !> 2.0 calculate LAI and N content of canopy layers
    !!
    !> 2.1 LAI of canopy layer
    !!
    DO icanopy = 1,ncanopy
      WHERE(lai(:) > eps8) ! avoid any calculation if there are no leaves
        ! LAI calculation
        WHERE(ubounds(icanopy) < lai(:) .OR. ABS(ubounds(icanopy) - lai(:)) < eps12 )
          ! ubounds contains cummulative LAI at the bottom of this layer 
          ! if canopy layer is fully developed
          lai_cl(:,icanopy) = dz(icanopy)
        ELSEWHERE(lbounds(icanopy) < lai(:) .AND. ubounds(icanopy) > lai(:))
          ! if canopy layer is only partially filled assign difference between full LAI and lc to canopy layer
          lai_cl(:,icanopy) = lai(:) - lbounds(icanopy)
        ENDWHERE
        
        ! make sure all LAI is represented by letting last layer go to the end of the canopy
        WHERE(icanopy==ncanopy .AND. ubounds(ncanopy) < lai(:) )
          lai_cl(:,ncanopy) = lai(:)
        ENDWHERE

        ! LAI mid-point of canopy layer
        WHERE(lai_cl(:,icanopy) > 0.0_wp)
          mid_point_lai_cl(:,icanopy) = (SUM(lai_cl(:,1:icanopy), DIM=2) - lbounds(icanopy)) / 2.0_wp
        ENDWHERE
      ENDWHERE ! lai(:) > eps8
    END DO   ! loop ncanopy
    
    !> 2.2 cummulative LAI above mid-point of canopy layer
    !!
    ! 1st layer
    WHERE(lai(:) > eps8) ! avoid any calculation if there are no leaves    
      cumm_lai_cl(:,1) = mid_point_lai_cl(:,1)
    ENDWHERE
    ! 2-x layers
    DO icanopy = 2,ncanopy
      WHERE(lai(:) > eps8) ! avoid any calculation if there are no leaves
        WHERE(lai_cl(:,icanopy) > eps8)
          cumm_lai_cl(:,icanopy) = lbounds(icanopy) + mid_point_lai_cl(:,icanopy)
        ENDWHERE
      ENDWHERE
    ENDDO 
    
    !> 2.3 N content of the canopy layer
    !!
    DO icanopy = 1,ncanopy
      WHERE(lai(:) > eps8)                        ! avoid any calculation if there are no leaves
        WHERE(lai(:) > ubounds(1))                ! LAI is present and LAI fills more than first layer 
          WHERE(cumm_lai_cl(:,icanopy) > eps12)   ! this layer has some leafs
            leaf_nitrogen_cl(:,icanopy) = leaf_nitrogen_tc(:) * exp(-kn * cumm_lai_cl(:,icanopy))
          ELSEWHERE
            leaf_nitrogen_cl(:,icanopy) = 0.0_wp  ! set N to zero if layer has no leafs 
          ENDWHERE
        ELSEWHERE
          leaf_nitrogen_cl(:,icanopy)   = 0.0_wp  ! set N to zero if layer has no leafs 
        ENDWHERE
      ENDWHERE
    ENDDO
    WHERE(lai(:) > eps8 .AND..NOT. lai(:) > ubounds(1))   ! LAI is present but 1st layer is not filled
      leaf_nitrogen_cl(:,1) = leaf_nitrogen_tc(:)         ! set 1st layer N to leaf_nitrogen_tc, but leave all other layers with zero N
    ENDWHERE


    !> 3.0 fractionation of N within canopy laer
    !!
    !> 3.1 fraction not associated with photosynthesis
    !! 
    WHERE(leaf_nitrogen_cl(:,:) > 0.0_wp)
       fn_oth_cl(:,:) = lctlib_k0_fn_struc - k1_fn_struc * leaf_nitrogen_cl(:,:)
    ELSEWHERE
       fn_oth_cl(:,:) = 1.0_wp
    END WHERE
    ! structural limitation to the indefinite increase in phot. N described by above equation
    WHERE (fn_oth_cl(:,:) < lctlib_fn_oth_min) fn_oth_cl(:,:) = lctlib_fn_oth_min


    !> 3.2 photosynthetic fractions
    !! prescribed, or dynamically calculated if flag_optimal_Nfraction = .TRUE.
    ! IF (.NOT.flag_optimal_Nfraction) THEN
        !> 3.2.1 case of presribed fractions
        WHERE(lai_cl(:,:) > eps8) ! canopy layers with leaves
           ! fraction in chlorophyll
           fn_chl_cl(:,:) = MIN(0.99_wp - fn_oth_cl(:,:), & 
                                (k0_fn_chl - k1_fn_chl * exp(-kfn_chl * cumm_lai_cl(:,:))) / chl2n)
           ! fraction in PepC
           fn_pepc_cl(:,:) = f_pepc
           ! fraction of electron-transport is dependent on N in electron transport versus N in light harvesting
           fn_et_cl(:,:) = (1._wp - (fn_oth_cl(:,:) + fn_chl_cl(:,:) + fn_pepc_cl(:,:) )) / &
                           (1._wp + jmax2n / (jmax2vcmax * vcmax2n))
           ! fraction of Rubisco is determined as the remainder
           fn_rub_cl(:,:) = 1._wp - fn_oth_cl(:,:) - fn_chl_cl(:,:) - fn_et_cl(:,:) - fn_pepc_cl(:,:)

           ! convert mmol N / m2 to mg N / g DW
           hlp3(:,:) = leaf_nitrogen_cl(:,:) * molar_mass_N / &
                         ((molar_mass_C / lctlib_sla) / carbon_per_dryweight_leaf)
           hlp4(:,:) = leaf_nitrogen_cl(:,:) / lctlib_np_leaf * molar_mass_P / &
                         ((molar_mass_C / lctlib_sla) / carbon_per_dryweight_leaf)
           !DE stuff Vmax = exp(k1 + k2 * log(P) + k3 * log(N))
           ! implied Vmax and Jmax in nmol / g DW / s
           !hlp1(:,:) = exp(4.44904_wp + 0.3472_wp*log10(hlp4(:,:)) + 0.49078*log10(hlp3(:,:)))
           !hlp2(:,:) = exp(5.49435_wp + 0.37345_wp*log10(hlp4(:,:)) + 0.41435_wp*log10(hlp3(:,:)))
           ! convert to f_rub_cl: Vmax = vmax2n * frub * leaf_n | µmol / m2 / s
           !                      Vmax1 = vmax2n * frub * leaf_n / DW * 1000 | nmol / gDW / s
           !                      exp(X) = vmax2n * frub * leaf_n
           !                      frub = exp(X) / leaf_n / vmax2n         
           !fn_rub_cl(:,:) = exp(4.44904_wp + 0.3472_wp*log(hlp4(:,:)) + 0.49078*log(hlp3(:,:))) / & 
           !                 (leaf_nitrogen_cl(:,:) * vcmax2n * 1000._wp / (molar_mass_C / lctlib_sla / 0.48_wp))
           !fn_et_cl(:,:) = exp(5.49435_wp + 0.37345_wp*log(hlp4(:,:)) + 0.41435_wp*log(hlp3(:,:))) / 4.0_wp / &
           !                 (leaf_nitrogen_cl(:,:) * jmax2n * 1000._wp / (molar_mass_C / lctlib_sla / 0.48_wp))
           !fn_oth_cl(:,:) = 1.0_wp - fn_chl_cl(:,:) - fn_et_cl(:,:) - fn_pepc_cl(:,:) - fn_rub_cl(:,:)

        ELSEWHERE ! canopy layers without leaves
           leaf_nitrogen_cl(:,:)  = 0.0_wp
           fn_oth_cl(:,:)         = 1.0_wp
           fn_chl_cl(:,:)         = 0.0_wp
           fn_et_cl(:,:)          = 0.0_wp
           fn_rub_cl(:,:)         = 0.0_wp
           fn_pepc_cl(:,:)        = 0.0_wp
        ENDWHERE

        ! for canopy layer > 1
        DO icanopy = 2,ncanopy
           WHERE(lai_cl(:,icanopy) > eps8)
              WHERE(ppfd_sunlit_cl(:,icanopy) < eps8) ! no radiation, take from above
                 ppfd_sunlit_cl(:,icanopy)   = ppfd_sunlit_cl(:,icanopy-1)
                 ppfd_shaded_cl(:,icanopy)   = ppfd_shaded_cl(:,icanopy-1)
                 fleaf_sunlit_cl(:,icanopy)  = fleaf_sunlit_cl(:,icanopy-1)
              ENDWHERE
           ENDWHERE
        ENDDO
                
    ! ELSE ! IF (.NOT.flag_optimal_Nfraction)
    !     !> 3.2.1 optimal N fraction case
    !     ! for canopy layer = 1
    !     icanopy = 1
    !     WHERE(lai_cl(:,icanopy) > eps8)
    !       WHERE(fn_chl_cl(:,icanopy) < eps8) ! new layer with new leaves
    !         ! first layers: use static allocation as first guess
    !         fn_chl_cl(:,icanopy)  = (k0_fn_chl - k1_fn_chl * exp(-kfn_chl * cumm_lai_cl(:,icanopy))) / chl2n
    !         ! fraction in PepC
    !         fn_pepc_cl(:,icanopy) = f_pepc
    !         ! fraction of electron-transport is dependent on N in electron transport versus N in light harvesting
    !         fn_et_cl(:,icanopy)   = (1._wp - (fn_oth_cl(:,icanopy) + fn_chl_cl(:,icanopy) + fn_pepc_cl(:,icanopy) )) / &
    !                                 (1._wp + jmax2n/(jmax2vcmax * vcmax2n))
    !         ! fraction of Rubisco is determined as the remainder
    !         fn_rub_cl(:,icanopy)  = 1._wp - fn_oth_cl(:,icanopy) - fn_chl_cl(:,icanopy) - fn_et_cl(:,icanopy) - fn_pepc_cl(:,icanopy)
    !       ENDWHERE
    !     ELSEWHERE
    !           leaf_nitrogen_cl(:,icanopy) = 0.0_wp
    !           fn_oth_cl(:,icanopy)        = 1.0_wp
    !           fn_chl_cl(:,icanopy)        = 0.0_wp
    !           fn_et_cl(:,icanopy)         = 0.0_wp
    !           fn_rub_cl(:,icanopy)        = 0.0_wp
    !           fn_pepc_cl(:,icanopy)       = 0.0_wp
    !     ENDWHERE
    ! 
    !     ! for canopy layer > 1
    !     DO icanopy = 2,ncanopy
    !        WHERE(lai_cl(:,icanopy) > eps8)
    !           WHERE(fn_chl_cl(:,icanopy) < eps8) ! other layers: copy from above layer
    !              fn_chl_cl(:,icanopy)        = fn_chl_cl(:,icanopy-1)
    !              fn_et_cl(:,icanopy)         = fn_et_cl(:,icanopy-1)
    !              fn_rub_cl(:,icanopy)        = fn_rub_cl(:,icanopy-1)
    !              fn_pepc_cl(:,icanopy)       = fn_pepc_cl(:,icanopy-1)
    !              ppfd_sunlit_cl(:,icanopy)   = ppfd_sunlit_cl(:,icanopy-1)
    !              ppfd_shaded_cl(:,icanopy)   = ppfd_shaded_cl(:,icanopy-1)
    !              fleaf_sunlit_cl(:,icanopy)  = fleaf_sunlit_cl(:,icanopy-1)
    !           ENDWHERE
    !        ELSEWHERE
    !           leaf_nitrogen_cl(:,icanopy) = 0.0_wp
    !           fn_oth_cl(:,icanopy)        = 1.0_wp
    !           fn_chl_cl(:,icanopy)        = 0.0_wp
    !           fn_et_cl(:,icanopy)         = 0.0_wp
    !           fn_rub_cl(:,icanopy)        = 0.0_wp
    !           fn_pepc_cl(:,icanopy)       = 0.0_wp
    !        ENDWHERE
    !     ENDDO ! icanopy
    ! 
    !     IF(optimise_n_fractions) THEN
    !       DO icanopy = 1,ncanopy
    !         CALL optimise_canopy_n_fractions(  lctlib_gmin                , & ! lctlib
    !                                            lctlib_g0                  , &
    !                                            lctlib_g1                  , &
    !                                            lctlib_t_jmax_omega        , &
    !                                            lctlib_ps_pathway          , & ! lctlib
    !                                            lai_cl(:,icanopy)          , & ! in
    !                                            leaf_nitrogen_cl(:,icanopy), &
    !                                            ppfd_sunlit_cl(:,icanopy)  , &
    !                                            ppfd_shaded_cl(:,icanopy)  , &
    !                                            fleaf_sunlit_cl(:,icanopy) , &
    !                                            t_air(:)                   , &
    !                                            t_acclim(:)                , &
    !                                            press_srf(:)               , &
    !                                            co2_mixing_ratio(:)        , &
    !                                            ga(:)                      , &
    !                                            beta_air(:)                , &
    !                                            beta_soa(:)                , &
    !                                            beta_soil_ps(:)            , &
    !                                            beta_sinklim_ps(:)         , &
    !                                            beta_soil_gs(:)            , &
    !                                            t_jmax_opt(:)              , &
    !                                            fn_oth_cl(:,icanopy)       , &
    !                                            fn_chl_cl(:,icanopy)       , & ! inout
    !                                            fn_rub_cl(:,icanopy)       , &
    !                                            fn_et_cl(:,icanopy)        , &
    !                                            fn_pepc_cl(:,icanopy)      )
    !       ENDDO
    !     ENDIF ! optimise_n_fractions
    ! ENDIF     ! end ELSE of IF (.NOT.flag_optimal_Nfraction)

  END SUBROUTINE calc_canopy_layers


#endif  
END MODULE mo_veg_canopy
