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
!!
MODULE mo_rad_update_radiation

  USE mo_kind,                  ONLY: wp, eps8
  USE mo_rad_constants


  IMPLICIT NONE

  PRIVATE
  PUBLIC :: update_radiation

CONTAINS

  !-----------------------------------------------------------------------------------------------------
  ! Main Task
  ! 
  !-----------------------------------------------------------------------------------------------------
  !> calculates surface radiation reflection and absorption
  !! 
  !! @todo currently no soil and snow albedo calculation in the update_radiation routine
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE update_radiation(tile, lctlib)

    USE mo_base_memory,           ONLY: t_slm_tile
    USE mo_base_constants,        ONLY: ncanop 
    USE mo_jsb_process_class,     ONLY: ATM2LAND, LAND2ATM, SEB, TURB, SOIL, HYDRO, VEG, HD, RAD, ASSIMI, PHENO

    USE mo_jsb_lctlib_class,      ONLY: t_lctlib_element
    USE mo_jsb_var_class,         ONLY: t_jsb_var_real2d, t_jsb_var_real3d, t_jsb_var_real1d

    USE mo_veg_memory,            ONLY: t_slm_plant
    USE mo_rad_memory_class,      ONLY: t_rad_memory
    USE mo_atm2land_memory_class, ONLY: t_atm2land_memory
    USE mo_land2atm_memory_class, ONLY: t_land2atm_memory

    USE mo_rad_constants 
    USE mo_atmland_constants,     ONLY: rad2par
   

    IMPLICIT NONE
    ! 0.1 InOut
    TYPE(t_slm_tile),INTENT(inout), TARGET :: tile    !< one tile
    TYPE(t_lctlib_element),INTENT(in)      :: lctlib  !< land-cover-type library - parameters for one pft
    ! 0.2 Local
    TYPE(t_slm_plant),POINTER              :: plant   ! plant of one tile
    REAL(wp),DIMENSION(ncanop)             :: arad_sunlit_vis_cl, & ! these variables need to be added to radiation memeory once 3D variables are allowed for
                                              arad_shaded_vis_cl, &
                                              fleaf_sunlit_vis_cl, & 
                                              arad_sunlit_nir_cl, &
                                              arad_shaded_nir_cl, &
                                              fleaf_sunlit_nir_cl

    ! 0.3 Memory allocation
    TYPE(t_atm2land_memory), POINTER :: atm2land__mem   ! dsl4jsb_Def_memory(ATM2LAND)
    TYPE(t_land2atm_memory), POINTER :: land2atm__mem
    TYPE(t_rad_memory), POINTER :: rad__mem
    REAL(wp), POINTER           :: swpar_srf_down     , &
                                   swpar_fdiffuse     , &
                                   cos_zenith_angle   , &
                                   swnir_srf_down     , &
                                   swvis_srf_down
    REAL(wp), POINTER           :: albedo             , &
                                   net_radiation      , &
                                   alb_vis_dir        , &
                                   alb_vis_dif        , &
                                   alb_nir_dir        , &
                                   alb_nir_dif
    REAL(wp), POINTER           :: alb_vis, & 
                                   alb_nir, &
                                   alb_vis_can, & 
                                   alb_nir_can, & 
                                   alb_vis_soil, & 
                                   alb_nir_soil
    REAL(wp), POINTER           :: arad_vis_can, & 
                                   arad_nir_can, & 
                                   arad_vis_soil, & 
                                   arad_nir_soil
    REAL(wp), POINTER           :: sw_srf_net, & 
                                   swvis_srf_net, & 
                                   swnir_srf_net, & 
                                   rad_srf_net

    SELECT TYPE (mm => tile%mem(ATM2LAND)%p)            ! dsl4jsb_Get_memory(ATM2LAND)
      TYPE IS(t_atm2land_memory)
        atm2land__mem => mm
    END SELECT
    SELECT TYPE (mm => tile%mem(LAND2ATM)%p)
      TYPE IS(t_land2atm_memory)
        land2atm__mem => mm
    END SELECT
    SELECT TYPE (mm => tile%mem(RAD)%p)
      TYPE IS(t_rad_memory)
        rad__mem => mm
    END SELECT

    swpar_srf_down   => atm2land__mem%swpar_srf_down%ptr
    swpar_fdiffuse   => atm2land__mem%swpar_fdiffuse%ptr
    cos_zenith_angle => atm2land__mem%cos_zenith_angle%ptr
    swnir_srf_down   => atm2land__mem%swnir_srf_down%ptr
    swvis_srf_down   => atm2land__mem%swvis_srf_down%ptr
    albedo            => land2atm__mem%albedo%ptr
    net_radiation     => land2atm__mem%net_radiation%ptr
    alb_vis_dir       => land2atm__mem%alb_vis_dir%ptr
    alb_vis_dif       => land2atm__mem%alb_vis_dif%ptr
    alb_nir_dir       => land2atm__mem%alb_nir_dir%ptr
    alb_nir_dif       => land2atm__mem%alb_nir_dif%ptr
    sw_srf_net      => rad__mem%sw_srf_net%ptr
    swvis_srf_net   => rad__mem%swvis_srf_net%ptr
    swnir_srf_net   => rad__mem%swnir_srf_net%ptr
    rad_srf_net     => rad__mem%rad_srf_net%ptr
    alb_vis         => rad__mem%alb_vis%ptr
    alb_nir         => rad__mem%alb_nir%ptr
    alb_vis_soil    => rad__mem%alb_vis_soil%ptr
    alb_nir_soil    => rad__mem%alb_nir_soil%ptr
    alb_vis_can     => rad__mem%alb_vis_can%ptr
    alb_nir_can     => rad__mem%alb_nir_can%ptr
    arad_vis_soil   => rad__mem%arad_vis_soil%ptr
    arad_nir_soil   => rad__mem%arad_nir_soil%ptr
    arad_vis_can    => rad__mem%arad_vis_can%ptr
    arad_nir_can    => rad__mem%arad_nir_can%ptr

    !  bi-directional pointer
    plant => tile%plant


    !> 2.0 Snow and soil albedo 
    !! for now, prescribed, should be based on JSBACH 4 routines, and include the effect of snow on the
    !! soil's albedo
    alb_vis_soil = def_alb_vis_soil 
    alb_nir_soil = def_alb_nir_soil

    !> 3.0 calculate VIS radiation budget (absorbed radiation on canopy layers and soil + albedo)
    !! 
    CALL calc_sw_srf_net(swvis_srf_down, & 
                     swpar_fdiffuse, cos_zenith_angle, & 
                     plant%lai, lctlib%sigma_vis, kbl0_vis, kdf0_vis, alb_vis_soil, &
                     arad_sunlit_vis_cl,arad_shaded_vis_cl,fleaf_sunlit_vis_cl, &
                     arad_vis_can,arad_vis_soil, &
                     alb_vis_can,alb_vis,swvis_srf_net,plant%lai_cl,plant%cumm_lai_cl) 

    !> 4.0 calculate NIR radiation budget (absorbed radiation on canopy layers and soil + albedo).
    !! Following Weiss & Norman, Agri. For. Met. 1985, the fraction of VIS and NIR 
    !! diffuse radiation are very similar 
    CALL calc_sw_srf_net(swnir_srf_down, & 
                     swpar_fdiffuse, cos_zenith_angle, & 
                     plant%lai, lctlib%sigma_nir, kbl0_nir, kdf0_nir, alb_nir_soil, &
                     arad_sunlit_nir_cl,arad_shaded_nir_cl,fleaf_sunlit_nir_cl, &
                     arad_nir_can,arad_nir_soil, &
                     alb_nir_can,alb_nir,swnir_srf_net,plant%lai_cl,plant%cumm_lai_cl) 

    !> 5.0 calculate albedo and net radiation for total shortwave radiation
    !!
    IF(swvis_srf_down + swnir_srf_down>eps8)THEN 
       albedo = ( alb_vis * swvis_srf_down + & 
                           alb_nir * swnir_srf_down ) / & 
                         (swvis_srf_down + swnir_srf_down) 
    ENDIF 
    sw_srf_net = swvis_srf_net + swnir_srf_net

    !> 6.0 diagnose net radiation (temporarily here for diagnostic purposes)
    !!
    !! Only for now longwave balance is zero
    rad_srf_net = sw_srf_net
    net_radiation = rad_srf_net

    !> 7.0 communitcate APAR to plant type
    !!
    plant%fleaf_sunlit_cl(:) = fleaf_sunlit_vis_cl(:) 
    IF(swvis_srf_down > eps8) THEN 
       plant%par_sunlit_cl(:) = arad_sunlit_vis_cl(:) * swpar_srf_down/swvis_srf_down * rad2par
       plant%par_shaded_cl(:) = arad_shaded_vis_cl(:) * swpar_srf_down/swvis_srf_down * rad2par
       plant%apar = arad_vis_can * swpar_srf_down/swvis_srf_down * rad2par 
    ELSE 
       plant%par_sunlit_cl(:) = 0.0_wp
       plant%par_shaded_cl(:) = 0.0_wp
       plant%apar = 0.0_wp
    ENDIF

  END SUBROUTINE update_radiation

  !-----------------------------------------------------------------------------------------------------
  !> calculate shortwave radiation budget of the tile
  !! 
  !! @par info on routine
  !!  This routine has been adapted from ADFs implementation of the radiation calculations
  !!  by Spitters, C J T. 1986. \n “Separating the Diffuse and Direct Component of Global Radiation
  !!  and Its Implications for Modeling Canopy Photosynthesis .2. Calculation of Canopy Photosynthesis.”
  !!  Agricultural and Forest Meteorology 38 (1-3): 231–42. \n
  !!  Equations numbers refer to the Spitters paper\n 
  !!  - The routine has been extended to diagnose albedo in the visible and NIR range. The implicit assumption 
  !!  is made that albedo in the PAR and VIS range are identical\n
  !!  - The routine has been adjusted to a two-stream approximation scheme by adding the upward flux 
  !!  from soil reflection\n
  !! @todo: calc_sw_srf_net does not account for clumping and non-leaf surfaces
  !! 
  !! @par Spitters Eq. 12 is wrong
  !!  accordingly to Wang2003, "radf4" was replaced by "-kbl" when calculating "idrdra"
  !! @par: APAR is calculated for all canopy layers irrespective of whether or not there are leaves 
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE calc_sw_srf_net(sw_srf_down,sw_fdiffuse,cos_zenith_angle,& 
                        lai,sigma,kbl0,kdf0,alb_soil,arad_sunlit_cl,arad_shaded_cl,fleaf_sunlit_cl, & 
                        arad_can,arad_soil,& 
                        alb_can,alb,sw_srf_net,lai_cl,cumm_lai_cl)

    USE mo_base_constants,     ONLY: ncanop,canopy_layer_depth

    IMPLICIT NONE
    ! 0.1 InOut
    REAL(wp), INTENT(in) ::    sw_srf_down             !< shortwave downward radiation flux (W m-2)
    REAL(wp), INTENT(in) ::    sw_fdiffuse             !< fraction of diffuse radiation (--)
    REAL(wp), INTENT(in) ::    cos_zenith_angle        !< solar angle (radians)
    REAL(wp), INTENT(in) ::    lai                     !< leaf area index (m-2 m-2)
    REAL(wp), INTENT(in) ::    sigma                   !< single leaf scattering coefficient 
    REAL(wp), INTENT(in) ::    kbl0                    !< extinction coefficient over black leaves 
    REAL(wp), INTENT(in) ::    kdf0                    !< extinction coefficient for diffuse radiation 
    REAL(wp), INTENT(in) ::    alb_soil                !< soil albedo 
    REAL(wp), DIMENSION(ncanop), INTENT(out) :: arad_sunlit_cl     !< absorbed radiation on sunlit leaves (W m-2) 
    REAL(wp), DIMENSION(ncanop), INTENT(out) :: arad_shaded_cl     !< absorbed radiation of shaded leaves (W m-2) 
    REAL(wp), DIMENSION(ncanop), INTENT(out) :: fleaf_sunlit_cl    !< fraction of leaves that is sunlit (--)
    REAL(wp), INTENT(out) ::   arad_can                !< absorbed radiation by the canopy (W m-2) 
    REAL(wp), INTENT(out) ::   arad_soil               !< absorbed radiation by the soil (W m-2)
    REAL(wp), INTENT(out) ::   alb_can                 !< albedo of the canopy 
    REAL(wp), INTENT(out) ::   alb                     !< albedo of the tile 
    REAL(wp), INTENT(out) ::   sw_srf_net              !< net shortwave radiation of the tile 
    REAL(wp), DIMENSION(:),INTENT(in) ::    lai_cl                !< soil albedo 
    REAL(wp), DIMENSION(:),INTENT(in) ::    cumm_lai_cl                !< soil albedo 
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
    REAL(wp), DIMENSION(ncanop) :: arad_shaded_up_cl   ! idem for canopy layers
    REAL(wp) :: temp1,temp2,radf1,radf2,radf3,radf3b,radf4,radf5 ! dummy variables
    REAL(wp), DIMENSION(ncanop) :: canopy_layer_depth_act  ! actual depth of canopy layers (account for incomplete lower bound) 
    REAL(wp) :: tmp_cos_zenith_angle                   ! local variable of solar angle (radians)
    INTEGER  :: icanop

    ! Initialise output 
    alb = 0.0_wp
    alb_can = 0.0_wp
    arad_can = 0.0_wp
    arad_can_up = 0.0_wp
    arad_soil = 0.0_wp
    arad_sunlit_cl(:) = 0.0_wp 
    arad_shaded_cl(:) = 0.0_wp
    arad_shaded_up_cl(:) = 0.0_wp
    fleaf_sunlit_cl(:) = 0.0_wp

    tmp_cos_zenith_angle = cos_zenith_angle
    IF(cos_zenith_angle < min_cos_zenith_angle_rad) & 
       tmp_cos_zenith_angle = min_cos_zenith_angle_rad  

    ! only do calculations if daylight is available, either direct or diffuse (in case of negative cos_zenith_angle)
    IF(sw_srf_down > eps8) THEN
    
       ! Downwelling direct radiation 
       irad_dir_down = (1.0_wp - sw_fdiffuse) * sw_srf_down
       ! Downwelling diffuse radiation 
       irad_dif_down = sw_fdiffuse * sw_srf_down

       IF( lai > eps8) THEN  

          !> Eq.1 constants for albedo calculations
          radf3b = (1.0_wp - sigma)
          temp1 = SQRT(radf3b)
          temp2 = 2.0_wp * (1.0_wp - temp1)/(1.0_wp + temp1)
   
          !> Eq.1  Fraction of radiation reflected by green closed canopy (fraction)
          rho = temp2/(1.0_wp + rho2sbeta * tmp_cos_zenith_angle)
   
          ! Extinction coefficient for radiation over black leaves in canopy
          kdf = kdf0 * temp1  ! temp1 == SQRT(1.0_wp - sigma)
          kbl = kbl0 * kdf/(0.8_wp * temp1 * tmp_cos_zenith_angle)
   
          !> Eq.10 & 11 Some factors for radiation calculations that change with time
          radf1 = (1.0_wp - rho) * irad_dif_down * kdf            ! eq 10
          radf2 = (1.0_wp - rho) * temp1 * irad_dir_down * kbl    ! eq 11
          !> Eq.12 absorbed direct radiation
          radf3 = radf3b * irad_dir_down * kbl                    ! eq 12
          radf4 = -temp1 * kbl
   
          ! light intensity at bottom of canopy
          irad_dif_down_bc = (1._wp - rho) * irad_dif_down * exp(-kdf * lai)
          irad_dir_down_bc = (1._wp - rho) * irad_dir_down * exp(radf4 * lai)
   
          ! diffuse upward shortwave flux from soil reflection and absorbed radiation of the soil
          irad_dif_up = alb_soil * (irad_dif_down_bc + irad_dir_down_bc)
          arad_soil = (1._wp - alb_soil) * (irad_dif_down_bc + irad_dir_down_bc)

          ! factor for the absorption calculation of the upward stream
          radf5 = (1.0_wp - rho) * irad_dif_up * kdf            ! eq 10

          ! upward shortwave diffuse flux reaching the top of the canopy
          irad_dif_up = (1.0_wp - rho) * irad_dif_up * exp(-kdf*lai)

          ! loop over canopy layers. All fluxes are evaluated at the centre of the canopy layer
          DO icanop = 1,ncanop
     
             IF (radf2 <= eps8) THEN ! case of only diffuse radiation
                ! Eq.10 Mean radiation penetrating foliage at lc 
                arad_shaded_cl(icanop) = radf1 * exp(-kdf * cumm_lai_cl(icanop)) + &    ! eq 10
                       radf5 * exp(-kdf * (max(lai-cumm_lai_cl(icanop),0.0_wp)))
                ! radiation penetrating foliage from below for albedo cacumm_lai_cl(icanop)ulation
                arad_shaded_up_cl(icanop) = radf5 * exp(-kdf * (max(lai-cumm_lai_cl(icanop),0.0_wp)))
             ELSE ! Case of direct and diffuse radiation
                ! Sunlit foliage at cumm_lai_cl(icanop) (fraction)
                ! kbl: Extinction coefficient for radiation over black leaves in canopy
                fleaf_sunlit_cl(icanop) = exp(-kbl * cumm_lai_cl(icanop))
                ! Eq.10 Diffuse radiation penetrating foliage at cumm_lai_cl(icanop) 
                idfa = radf1 * exp(-kdf * cumm_lai_cl(icanop)) + & 
                       radf5 * exp(-kdf * (max(lai-cumm_lai_cl(icanop),0.0_wp)))
                ! radiation penetrating foliage from below for albedo cacumm_lai_cl(icanop)ulation
                arad_shaded_up_cl(icanop) = radf5 * exp(-kdf * (max(lai-cumm_lai_cl(icanop),0.0_wp)))
                ! Direct radiation penetrating sunlit foliage at cumm_lai_cl(icanop) 
                idra = radf2 * exp(radf4 * cumm_lai_cl(icanop))
                ! Eq.12 Mean direct radiation remaining direct at cumm_lai_cl(icanop) 
                ! radf4 replaced by -kbl according to Wang2003 (error in Spitters 1986 Eqn.12)
                idrdra = radf3 * exp(-kbl * cumm_lai_cl(icanop))  
                ! Radiation penetrating shaded foliage at cumm_lai_cl(icanop) (umol/m2/s)
                ! Eq.13 Sum of diffuse radiation and diffuse component of direct flux
                arad_shaded_cl(icanop) = idfa + max((idra - idrdra), 0.0_wp)
                ! Radiation penetrating sunlit foliage at lc (umol/m2/s)
                ! sum of diffuse and direct radiation
                arad_sunlit_cl(icanop) = arad_shaded_cl(icanop) + radf3
             ENDIF 
   
             ! Integrate absorbed radiation over actually existing canopy layers
             IF(lai_cl(icanop)>eps8)THEN
                ! total absorbed radiation
                arad_can = arad_can + ( arad_shaded_cl(icanop) * (1._wp-fleaf_sunlit_cl(icanop)) +  & 
                                       arad_sunlit_cl(icanop) * fleaf_sunlit_cl(icanop) ) * & 
                                      canopy_layer_depth_act(icanop)
                ! for albedo calculation: radiation absorbed from the upward stream
                arad_can_up = arad_can_up + arad_shaded_up_cl(icanop) * canopy_layer_depth_act(icanop)
             ENDIF
  
          ENDDO
   
          ! albedo of the canopy, as seen from the top of the canopy
          ! inferred from the difference between top and bottom of the 
          ! canopy light intensities minus the amount of absorbed downward radiation
          alb_can = 1.0_wp - (irad_dif_down_bc + irad_dir_down_bc + (arad_can-arad_can_up)) / & 
                             (irad_dir_down + irad_dif_down) 

          ! total upward shortwave radiation flux at the top of the canopy:
          ! soil-based flux + reflectance of the canopy 
          irad_dif_up = irad_dif_up + & 
                        alb_can * (irad_dir_down + irad_dif_down - (irad_dif_down_bc + irad_dir_down_bc))
       ELSE
          ! case of no LAI: absorped radiation and upward flux from soil only 
          irad_dif_up = alb_soil * (irad_dif_down + irad_dir_down) 
          arad_soil = (1._wp - alb_soil ) * (irad_dif_down + irad_dir_down)
       ENDIF

       ! albedo and net shortwave radiation of the tile
       alb = irad_dif_up / (irad_dir_down + irad_dif_down)
       sw_srf_net = (1.0_wp - alb) * sw_srf_down

    ENDIF

  END SUBROUTINE calc_sw_srf_net

END MODULE mo_rad_update_radiation
