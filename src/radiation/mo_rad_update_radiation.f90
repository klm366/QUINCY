!>
!! @brief Contains the routines for the shortwave surface radiation budget (fapar & albedo) 
!!
!! @author
!!  Sönke Zaehle, MPI-BGC Jena 
!!
!! @par Revision History
!! First version by Sönke Zaehle (2015-12-04)
!!
!! 18-08-2016: update to include VIS And NIR albedo calculation for the canopy
!! 24-08-2016: update to include two-stream approximation and soil albedo 
!! 02-09-2016: update to include clumping 
!!
MODULE mo_rad_update_radiation

#include "dsl_jsbach.inc"

  USE mo_kind,                  ONLY: wp, eps8
  USE mo_rad_constants


  IMPLICIT NONE

  PRIVATE
  PUBLIC :: update_radiation

  CHARACTER(len=*), PARAMETER :: modname = 'mo_rad_update_radiation'

CONTAINS

  !-----------------------------------------------------------------------------------------------------
  ! Main Task
  ! 
  !-----------------------------------------------------------------------------------------------------
  !> calculates surface radiation reflection and absorption
  !! 
  !! @todo currently no soil and snow albedo calculation in the update_radiation routine
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE update_radiation(tile, options)

    USE mo_jsb_class,             ONLY: Get_model
    USE mo_base_memory,           ONLY: t_slm_tile, t_model, t_jsb_task_options
    USE mo_assimi_config_class,   ONLY: assimi_config   ! e.g. ncanopy
    USE mo_jsb_process_class,     ONLY: ATM2LAND_, LAND2ATM_, SEB_, TURB_, SOIL_, HYDRO_, VEG_, HD_, RAD_, ASSIMI_, PHENO_, SB_
    USE mo_jsb_lctlib_class,      ONLY: t_lctlib_element
    USE mo_rad_constants,         ONLY: rad2ppfd, def_alb_vis_soil, def_alb_nir_soil, kbl0_vis, kdf0_vis, kbl0_nir, kdf0_nir, &
                                        rho2sbeta, k_csf, min_cos_zenith_angle_rad

    ! Use of process memories
    dsl4jsb_Use_memory(ATM2LAND_)   !  USE mo_atm2land_memory_class, ONLY: t_atm2land_memory
    dsl4jsb_Use_memory(LAND2ATM_)   !  USE mo_land2atm_memory_class, ONLY: t_land2atm_memory
    dsl4jsb_Use_memory(RAD_)        !  USE mo_rad_memory_class,      ONLY: t_rad_memory
    dsl4jsb_Use_memory(VEG_)        !  USE mo_veg_memory_class,      ONLY: t_veg_memory
    

    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    TYPE(t_slm_tile),         INTENT(inout) :: tile             !< one tile
    TYPE(t_jsb_task_options), INTENT(in)  :: options            !< model options
    ! ---------------------------
    ! 0.2 Local
    TYPE(t_model),          POINTER        :: model
    TYPE(t_lctlib_element), POINTER        :: lctlib     !< land-cover-type library - parameter across pft's
    INTEGER                               :: iblk, ics, ice, nc
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':update_radiation'
    ! ---------------------------
    ! 0.3 Declare Memory
    ! Declare process configuration and memory Pointers
    dsl4jsb_Def_memory(ATM2LAND_)                   !     TYPE(t_atm2land_memory), POINTER :: atm2land__mem
    dsl4jsb_Def_memory(LAND2ATM_)                   !     TYPE(t_land2atm_memory), POINTER :: land2atm__mem
    dsl4jsb_Def_memory(RAD_)                        !     TYPE(t_rad_memory),      POINTER :: rad__mem
    dsl4jsb_Def_memory(VEG_)                        !     TYPE(t_veg_memory),      POINTER :: veg__mem
    ! Declare pointers to variables in memory
    dsl4jsb_Real2D_onChunk :: sw_srf_net              ! RAD_
    dsl4jsb_Real2D_onChunk :: swvis_srf_net
    dsl4jsb_Real2D_onChunk :: swnir_srf_net
    dsl4jsb_Real2D_onChunk :: rad_srf_net
    dsl4jsb_Real2D_onChunk :: alb_vis
    dsl4jsb_Real2D_onChunk :: alb_nir
    dsl4jsb_Real2D_onChunk :: alb_vis_soil
    dsl4jsb_Real2D_onChunk :: alb_nir_soil
    dsl4jsb_Real2D_onChunk :: alb_vis_can
    dsl4jsb_Real2D_onChunk :: alb_nir_can
    dsl4jsb_Real2D_onChunk :: arad_vis_soil
    dsl4jsb_Real2D_onChunk :: arad_nir_soil
    dsl4jsb_Real2D_onChunk :: arad_vis_can
    dsl4jsb_Real2D_onChunk :: arad_nir_can
    dsl4jsb_Real2D_onChunk :: appfd
    dsl4jsb_Real3D_onChunk :: ppfd_sunlit_cl
    dsl4jsb_Real3D_onChunk :: ppfd_shaded_cl
    dsl4jsb_Real3D_onChunk :: arad_sunlit_vis_cl
    dsl4jsb_Real3D_onChunk :: arad_shaded_vis_cl
    dsl4jsb_Real3D_onChunk :: fleaf_sunlit_vis_cl
    dsl4jsb_Real3D_onChunk :: arad_sunlit_nir_cl
    dsl4jsb_Real3D_onChunk :: arad_shaded_nir_cl
    dsl4jsb_Real3D_onChunk :: fleaf_sunlit_nir_cl
    dsl4jsb_Real2D_onChunk :: swpar_srf_down          ! ATM2LAND_
    dsl4jsb_Real2D_onChunk :: swpar_fdiffuse
    dsl4jsb_Real2D_onChunk :: cos_zenith_angle
    dsl4jsb_Real2D_onChunk :: swnir_srf_down
    dsl4jsb_Real2D_onChunk :: swvis_srf_down
    dsl4jsb_Real2D_onChunk :: albedo                  ! LAND2ATM_
    dsl4jsb_Real2D_onChunk :: net_radiation
!! *slmtemp
    TYPE(t_veg_memory),      POINTER :: veg         ! local pointer to veg__mem for SLM usability
    ! Get local variables from options argument
    iblk    = options%iblk
    ics     = options%ics
    ice     = options%ice
    nc      = options%nc
    ! ---------------------------
    ! 0.4 Process Activity, Debug Option
!    IF (.NOT. tile%Is_process_active(VEG_)) RETURN   ! ask for all processes considered here ?
!    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')
    ! ---------------------------
    ! 0.5 Get Memory
    ! Get process config
!    dsl4jsb_Get_config(xx_)
    ! Get process memories
    dsl4jsb_Get_memory(ATM2LAND_)
    dsl4jsb_Get_memory(LAND2ATM_)
    dsl4jsb_Get_memory(RAD_)
    dsl4jsb_Get_memory(VEG_)
    ! Get process variables (Set pointers to variables in memory)
    model  => Get_model(tile%owner_model_id)
    lctlib => model%lctlib(tile%lcts%lib_id)
    ! ---------------------------
    dsl4jsb_Get_var2D_onChunk(ATM2LAND_, swpar_srf_down)
    dsl4jsb_Get_var2D_onChunk(ATM2LAND_, swpar_fdiffuse)
    dsl4jsb_Get_var2D_onChunk(ATM2LAND_, cos_zenith_angle)
    dsl4jsb_Get_var2D_onChunk(ATM2LAND_, swnir_srf_down)
    dsl4jsb_Get_var2D_onChunk(ATM2LAND_, swvis_srf_down)
    ! ---------------------------
    dsl4jsb_Get_var2D_onChunk(LAND2ATM_, albedo)
    dsl4jsb_Get_var2D_onChunk(LAND2ATM_, net_radiation)
    ! ---------------------------
    dsl4jsb_Get_var2D_onChunk(RAD_, sw_srf_net)
    dsl4jsb_Get_var2D_onChunk(RAD_, swvis_srf_net)
    dsl4jsb_Get_var2D_onChunk(RAD_, swnir_srf_net)
    dsl4jsb_Get_var2D_onChunk(RAD_, rad_srf_net)
    dsl4jsb_Get_var2D_onChunk(RAD_, alb_vis)
    dsl4jsb_Get_var2D_onChunk(RAD_, alb_nir)
    dsl4jsb_Get_var2D_onChunk(RAD_, alb_vis_soil)
    dsl4jsb_Get_var2D_onChunk(RAD_, alb_nir_soil)
    dsl4jsb_Get_var2D_onChunk(RAD_, alb_vis_can)
    dsl4jsb_Get_var2D_onChunk(RAD_, alb_nir_can)
    dsl4jsb_Get_var2D_onChunk(RAD_, arad_vis_soil)
    dsl4jsb_Get_var2D_onChunk(RAD_, arad_nir_soil)
    dsl4jsb_Get_var2D_onChunk(RAD_, arad_vis_can)
    dsl4jsb_Get_var2D_onChunk(RAD_, arad_nir_can)
    dsl4jsb_Get_var2D_onChunk(RAD_, appfd)
    dsl4jsb_Get_var3D_onChunk(RAD_, ppfd_sunlit_cl)
    dsl4jsb_Get_var3D_onChunk(RAD_, ppfd_shaded_cl)
    dsl4jsb_Get_var3D_onChunk(RAD_, arad_sunlit_vis_cl)
    dsl4jsb_Get_var3D_onChunk(RAD_, arad_shaded_vis_cl)
    dsl4jsb_Get_var3D_onChunk(RAD_, fleaf_sunlit_vis_cl)
    dsl4jsb_Get_var3D_onChunk(RAD_, arad_sunlit_nir_cl)
    dsl4jsb_Get_var3D_onChunk(RAD_, arad_shaded_nir_cl)
    dsl4jsb_Get_var3D_onChunk(RAD_, fleaf_sunlit_nir_cl)


!! *slmtemp
    veg       => veg__mem 


    ! ------------------------------------------------------------------------------------------------------------
    ! Go science
    ! ------------------------------------------------------------------------------------------------------------



    !> 1.0 Snow and soil albedo 
    !! 
    !! for now, prescribed, should be based on JSBACH 4 routines, and include the effect of snow on the
    !! soil's albedo
    alb_vis_soil(:) = def_alb_vis_soil
    alb_nir_soil(:) = def_alb_nir_soil

    !> 2.0 calculate VIS radiation budget (absorbed radiation on canopy layers and soil + albedo)
    !! 
    CALL calc_sw_srf_net(swvis_srf_down(1), & 
                     swpar_fdiffuse(1), cos_zenith_angle(1), & 
                     veg%lai%ptr, veg%lai_cl%ptr(:,1), veg%cumm_lai_cl%ptr(:,1),& 
                     lctlib%sigma_vis, lctlib%omega_clumping, lctlib%crown_shape_factor, & 
                     kbl0_vis, kdf0_vis, alb_vis_soil(1), &
                     arad_sunlit_vis_cl(1,:),arad_shaded_vis_cl(1,:),fleaf_sunlit_vis_cl(1,:), &
                     arad_vis_can(1),arad_vis_soil(1), &
                     alb_vis_can(1),alb_vis(1),swvis_srf_net(1)) 

    !> 3.0 calculate NIR radiation budget (absorbed radiation on canopy layers and soil + albedo)
    !! 
    !! Following Weiss & Norman, Agri. For. Met. 1985, the fraction of VIS and NIR 
    !! diffuse radiation are very similar 
    CALL calc_sw_srf_net(swnir_srf_down(1), & 
                     swpar_fdiffuse(1), cos_zenith_angle(1), & 
                     veg%lai%ptr, veg%lai_cl%ptr(:,1), veg%cumm_lai_cl%ptr(:,1), &
                     lctlib%sigma_nir, lctlib%omega_clumping, lctlib%crown_shape_factor, & 
                     kbl0_nir, kdf0_nir, alb_nir_soil(1), &
                     arad_sunlit_nir_cl(1,:),arad_shaded_nir_cl(1,:),fleaf_sunlit_nir_cl(1,:), &
                     arad_nir_can(1),arad_nir_soil(1), &
                     alb_nir_can(1),alb_nir(1),swnir_srf_net(1)) 

    !> 4.0 calculate albedo and net radiation for total shortwave radiation
    !!
    WHERE(swvis_srf_down(:) + swnir_srf_down(:) > eps8)
       albedo(:) = ( alb_vis(:) * swvis_srf_down(:)   + & 
                    alb_nir(:) * swnir_srf_down(:) ) / & 
                    (swvis_srf_down(:) + swnir_srf_down(:)) 
    END WHERE 
    sw_srf_net(:) = swvis_srf_net(1) + swnir_srf_net(1)

    !> 5.0 diagnose net radiation (temporarily here for diagnostic purposes)
    !! 
    !! Only for now longwave balance is zero
    rad_srf_net(:) = sw_srf_net(:)
    net_radiation(:) = rad_srf_net(:)

    !> 6.0 convert radiation to Photosynthetically Active Photon Flux Density (ppfd)
    !!
    veg%fleaf_sunlit_cl%ptr(:,1) = fleaf_sunlit_vis_cl(1,:) 
    IF(swvis_srf_down(1) > eps8) THEN
       ppfd_sunlit_cl(1,:) = arad_sunlit_vis_cl(1,:) * swpar_srf_down(1)/swvis_srf_down(1) * rad2ppfd
       ppfd_shaded_cl(1,:) = arad_shaded_vis_cl(1,:) * swpar_srf_down(1)/swvis_srf_down(1) * rad2ppfd
       appfd(:)            = arad_vis_can(:)         * swpar_srf_down(:)/swvis_srf_down(:) * rad2ppfd 
    ELSE
       ppfd_sunlit_cl(:,:) = 0.0_wp
       ppfd_shaded_cl(:,:) = 0.0_wp
       appfd(:)            = 0.0_wp
    END IF

  END SUBROUTINE update_radiation


  !-----------------------------------------------------------------------------------------------------
  !> calculate shortwave radiation budget of the tile
  !! 
  !! @par info on routine
  !!  This routine has been adapted from ADFs implementation of the radiation calculations
  !!  by Spitters, C J T. 1986. \n “Separating the Diffuse and Direct Component of Global Radiation
  !!  and Its Implications for Modeling Canopy Photosynthesis .2. Calculation of Canopy Photosynthesis.”
  !!  Agricultural and Forest Meteorology 38 (1-3): 231–42. \n
  !!  Equations numbers refer to the Spitters paper \n 
  !!  - The routine has been extended to diagnose albedo in the visible and NIR range. The implicit assumption 
  !!  is made that albedo in the PAR and VIS range are identical \n
  !!  - The routine has been adjusted to a two-stream approximation scheme by adding the upward flux 
  !!  from soil reflection \n
  !!  - The routine has been adjusted to account for canopy clumping (Campbell & Norman 1998) \n
  !! @todo: calc_sw_srf_net does not account non-leaf surfaces
  !! 
  !! @par Spitters Eq. 12 is wrong
  !!  accordingly to Wang2003, "radf4" was replaced by "-kbl" when calculating "idrdra"
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE calc_sw_srf_net( sw_srf_down         , &   ! in
                              sw_fdiffuse         , &   ! in
                              cos_zenith_angle    , &   ! in 
                              lai                 , &   ! in
                              lai_cl              , &   ! in
                              cumm_lai_cl         , &   ! in 
                              sigma               , &   ! in
                              omega_clumping      , &   ! in
                              crown_shape_factor  , &   ! in
                              kbl0                , &   ! in
                              kdf0                , &   ! in
                              alb_soil            , &   ! in
                              arad_sunlit_cl      , &   ! out
                              arad_shaded_cl      , &   ! out
                              fleaf_sunlit_cl     , &   ! out 
                              arad_can            , &   ! out
                              arad_soil           , &   ! out
                              alb_can             , &   ! out
                              alb                 , &   ! out
                              sw_srf_net)

    USE mo_assimi_config_class,   ONLY: assimi_config   ! e.g. ncanopy

    IMPLICIT NONE
    ! 0.1 InOut
    REAL(wp)              , INTENT(in) :: sw_srf_down             !< shortwave downward radiation flux (W m-2)
    REAL(wp)              , INTENT(in) :: sw_fdiffuse             !< fraction of diffuse radiation (--)
    REAL(wp)              , INTENT(in) :: cos_zenith_angle        !< solar angle (radians)
    REAL(wp)              , INTENT(in) :: lai                     !< leaf area index (m-2 m-2)
    REAL(wp), DIMENSION(:), INTENT(in) :: lai_cl                  !< leaf area index of canopy layer (m-2 m-2) 
    REAL(wp), DIMENSION(:), INTENT(in) :: cumm_lai_cl             !< cummulative leaf area index above centre of canopy layer
    REAL(wp)              , INTENT(in) :: sigma                   !< single leaf scattering coefficient 
    REAL(wp)              , INTENT(in) :: omega_clumping          !< clumping factor 
    REAL(wp)              , INTENT(in) :: crown_shape_factor      !< crown-shape factor 
    REAL(wp)              , INTENT(in) :: kbl0                    !< extinction coefficient over black leaves 
    REAL(wp)              , INTENT(in) :: kdf0                    !< extinction coefficient for diffuse radiation 
    REAL(wp)              , INTENT(in) :: alb_soil                !< soil albedo 
    REAL(wp), DIMENSION(assimi_config%ncanopy), INTENT(out) :: arad_sunlit_cl    !< absorbed radiation on sunlit leaves (W m-2) 
    REAL(wp), DIMENSION(assimi_config%ncanopy), INTENT(out) :: arad_shaded_cl    !< absorbed radiation of shaded leaves (W m-2) 
    REAL(wp), DIMENSION(assimi_config%ncanopy), INTENT(out) :: fleaf_sunlit_cl   !< fraction of leaves that is sunlit (--)
    REAL(wp)                                  , INTENT(out) :: arad_can          !< absorbed radiation by the canopy (W m-2) 
    REAL(wp)                                  , INTENT(out) :: arad_soil         !< absorbed radiation by the soil (W m-2)
    REAL(wp)                                  , INTENT(out) :: alb_can           !< albedo of the canopy 
    REAL(wp)                                  , INTENT(out) :: alb               !< albedo of the tile 
    REAL(wp)                                  , INTENT(out) :: sw_srf_net        !< net shortwave radiation of the tile 
    ! 0.2 Local
    REAL(wp) :: irad_dir_down,irad_dif_down            ! downward direct and diffuse radiation 
    REAL(wp) :: irad_dir_down_bc,irad_dif_down_bc      ! light intensity at the bottom of the canopy
    REAL(wp) :: irad_dif_up                            ! upward diffuse shortwave radiation
    REAL(wp) :: rho                                    ! fraction reflected by green canopy (--)
    REAL(wp) :: kdf,kbl                                ! extinction coefficients for radiation over black canopy
    REAL(wp) :: lc                                     ! cummulative LAI over mid point of a canopy layer
    REAL(wp) :: idfa                                   ! diffuse radiation penetrating foliage at lc 
    REAL(wp) :: idra                                   ! direct radiation penetrating sunlit foliage at lc 
    REAL(wp) :: idrdra                                 ! direct radiation remaining direct at lc 
    REAL(wp) :: arad_can_up                            ! absorbed radiation from the upward stream
    REAL(wp), DIMENSION(assimi_config%ncanopy) :: arad_shaded_up_cl   ! idem for canopy layers
    REAL(wp) :: tmp_cos_zenith_angle                   ! local variable of solar angle (radians)
    REAL(wp) :: k4_csf                                 ! ..
    REAL(wp) :: omega                                  ! .. 
    INTEGER  :: icanopy                                 ! iteration over canopy layers
    REAL(wp) :: temp1,temp2,radf1,radf2,radf3,radf3b,radf4,radf5 ! dummy variables
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_sw_srf_net'

    ! Initialise output 
    alb                   = 0.0_wp
    alb_can               = 0.0_wp
    arad_can              = 0.0_wp
    arad_can_up           = 0.0_wp
    arad_soil             = 0.0_wp
    arad_sunlit_cl(:)     = 0.0_wp 
    arad_shaded_cl(:)     = 0.0_wp
    arad_shaded_up_cl(:)  = 0.0_wp
    fleaf_sunlit_cl(:)    = 0.0_wp

    !> enforce tmp_cos_zenith_angle to be larger than 'min_cos_zenith_angle_rad'
    tmp_cos_zenith_angle = cos_zenith_angle
    IF(cos_zenith_angle < min_cos_zenith_angle_rad) tmp_cos_zenith_angle = min_cos_zenith_angle_rad  

    !> @par IF only do calculations if daylight is available, either direct or diffuse (in case of negative cos_zenith_angle)
    IF(sw_srf_down > eps8) THEN
    
       !> Downwelling direct radiation 
       irad_dir_down = (1.0_wp - sw_fdiffuse) * sw_srf_down
       !> Downwelling diffuse radiation 
       irad_dif_down = sw_fdiffuse * sw_srf_down

       !> @par IF lai > eps8
       IF( lai > eps8) THEN  

          !> Eq.1 constants for albedo calculations
          radf3b = (1.0_wp - sigma)
          temp1 = SQRT(radf3b)
          temp2 = 2.0_wp * (1.0_wp - temp1)/(1.0_wp + temp1)
   
          !> Eq.1  Fraction of radiation reflected by green closed canopy (fraction)
          rho = temp2/(1.0_wp + rho2sbeta * tmp_cos_zenith_angle)
  
          !> clumping factor (eq 15.35) Campbell & Norman (1998) as a function of sun-angle
          omega = omega_clumping / (omega_clumping + (1._wp - omega_clumping) * & 
                  exp(-k_csf * ACOS(tmp_cos_zenith_angle)**crown_shape_factor)) 
 
          !> Extinction coefficient for radiation over black leaves in canopy
          kdf = kdf0 * temp1 * omega ! temp1 == SQRT(1.0_wp - sigma)
          kbl = kbl0 * kdf/(0.8_wp * temp1 * tmp_cos_zenith_angle) 

          !> Eq.10 & 11 Some factors for radiation calculations that change with time
          radf1 = (1.0_wp - rho) * irad_dif_down * kdf            ! eq 10
          radf2 = (1.0_wp - rho) * temp1 * irad_dir_down * kbl    ! eq 11
          !> Eq.12 absorbed direct radiation
          radf3 = radf3b * irad_dir_down * kbl                    ! eq 12
          radf4 = -temp1 * kbl
   
          !> light intensity at bottom of canopy
          irad_dif_down_bc = (1._wp - rho) * irad_dif_down * exp(-kdf * lai)
          irad_dir_down_bc = (1._wp - rho) * irad_dir_down * exp(radf4 * lai)
   
          !> diffuse upward shortwave flux from soil reflection and absorbed radiation of the soil
          irad_dif_up = alb_soil * (irad_dif_down_bc + irad_dir_down_bc)
          arad_soil = (1._wp - alb_soil) * (irad_dif_down_bc + irad_dir_down_bc)

          !> factor for the absorption calculation of the upward stream
          radf5 = (1.0_wp - rho) * irad_dif_up * kdf            ! eq 10

          !> upward shortwave diffuse flux reaching the top of the canopy
          irad_dif_up = (1.0_wp - rho) * irad_dif_up * exp(-kdf*lai)

          !> @par loop over canopy layers. All fluxes are evaluated at the centre of the canopy layer
          DO icanopy = 1,assimi_config%ncanopy
             !> @par IF leaf area index of canopy layer > eps8: Integrate absorbed radiation over actually existing canopy layers
             IF(lai_cl(icanopy)>eps8)THEN
     
                !> @par IF radf2 <= eps8
                IF (radf2 <= eps8) THEN ! case of only diffuse radiation
                   !> Eq.10 Mean radiation penetrating foliage at lc 
                   arad_shaded_cl(icanopy) = radf1 * exp(-kdf * cumm_lai_cl(icanopy)) + &    ! eq 10
                          radf5 * exp(-kdf * (max(lai-cumm_lai_cl(icanopy),0.0_wp)))
                   !> radiation penetrating foliage from below for albedo cacumm_lai_cl(icanopy)ulation
                   arad_shaded_up_cl(icanopy) = radf5 * exp(-kdf * (max(lai-cumm_lai_cl(icanopy),0.0_wp)))
                !> @par ELSE Case of direct and diffuse radiation
                ELSE 
                   !> Sunlit foliage at cumm_lai_cl(icanopy) (fraction)
                   !! kbl: Extinction coefficient for radiation over black leaves in canopy
                   fleaf_sunlit_cl(icanopy) = exp(-kbl * cumm_lai_cl(icanopy))
                   !> Eq.10 Diffuse radiation penetrating foliage at cumm_lai_cl(icanopy) 
                   idfa = radf1 * exp(-kdf * cumm_lai_cl(icanopy)) + & 
                          radf5 * exp(-kdf * (max(lai-cumm_lai_cl(icanopy),0.0_wp)))
                   !> radiation penetrating foliage from below for albedo accumulation ?
                   arad_shaded_up_cl(icanopy) = radf5 * exp(-kdf * (max(lai-cumm_lai_cl(icanopy),0.0_wp)))
                   !> Direct radiation penetrating sunlit foliage at cumm_lai_cl(icanopy) 
                   idra = radf2 * exp(radf4 * cumm_lai_cl(icanopy))
                   !> Eq.12 Mean direct radiation remaining direct at cumm_lai_cl(icanopy) 
                   !! radf4 replaced by -kbl according to Wang2003 (error in Spitters 1986 Eqn.12)
                   idrdra = radf3 * exp(-kbl * cumm_lai_cl(icanopy))  
                   !> Radiation penetrating shaded foliage at cumm_lai_cl(icanopy) (umol/m2/s)
                   !! Eq.13 Sum of diffuse radiation and diffuse component of direct flux
                   arad_shaded_cl(icanopy) = idfa + max((idra - idrdra), 0.0_wp)
                   !> Radiation penetrating sunlit foliage at lc (umol/m2/s)
                   !! sum of diffuse and direct radiation
                   arad_sunlit_cl(icanopy) = arad_shaded_cl(icanopy) + radf3
                ENDIF
                !> ENDIF

                !> total absorbed radiation
                arad_can = arad_can + ( arad_shaded_cl(icanopy) * (1._wp-fleaf_sunlit_cl(icanopy)) + & 
                                        arad_sunlit_cl(icanopy) * fleaf_sunlit_cl(icanopy) )       * & 
                                        lai_cl(icanopy)
                !> for albedo calculation: radiation absorbed from the upward stream
                arad_can_up = arad_can_up + arad_shaded_up_cl(icanopy) * lai_cl(icanopy)
             ENDIF
             !> ENDIF
          ENDDO
          !> END LOOP

          !> albedo of the canopy, as seen from the top of the canopy
          !! inferred from the difference between top and bottom of the 
          !! canopy light intensities minus the amount of absorbed downward radiation
          alb_can = 1.0_wp - (irad_dif_down_bc + irad_dir_down_bc + (arad_can-arad_can_up)) / & 
                             (irad_dir_down + irad_dif_down) 

          !> total upward shortwave radiation flux at the top of the canopy:
          !! soil-based flux + reflectance of the canopy 
          irad_dif_up = irad_dif_up + & 
                        alb_can * (irad_dir_down + irad_dif_down - (irad_dif_down_bc + irad_dir_down_bc))
       ELSE
          !> case of no LAI: absorped radiation and upward flux from soil only 
          irad_dif_up = alb_soil            * (irad_dif_down + irad_dir_down) 
          arad_soil   = (1._wp - alb_soil ) * (irad_dif_down + irad_dir_down)
       ENDIF

       !> albedo and net shortwave radiation of the tile
       alb        = irad_dif_up / (irad_dir_down + irad_dif_down)
       sw_srf_net = (1.0_wp - alb) * sw_srf_down
 
    ELSE 

       !> albedo and net shortwave radiation of the tile
       alb        = 0.0_wp 
       sw_srf_net = 0.0_wp 

    ENDIF

  END SUBROUTINE calc_sw_srf_net

END MODULE mo_rad_update_radiation
