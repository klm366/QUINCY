!>
!! @brief diverse subroutines and functions to mimick soil for the stand-alone plant-model
!!
!! @author
!!  Sönke Zaehle, MPI-BGC Jena
!!
!! @par Revision History
!! First version by Sönke Zaehle (2016-03-31)
!!
MODULE mo_veg_soil

  USE mo_kind,               ONLY: wp, eps8
  USE mo_jsb_lctlib_class,   ONLY: t_lctlib_element
  USE mo_veg_config_class,   ONLY: veg_config
  USE mo_veg_constants


  IMPLICIT NONE

  PRIVATE
  PUBLIC :: update_vegsoil_bgc, update_vegsoil_physics, update_vegsoil_pool


CONTAINS

  !-----------------------------------------------------------------------------------------------------
  ! Main Task
  ! 
  !-----------------------------------------------------------------------------------------------------
  !> routine to simulate N and P uptake from soil in the absence of a proper soil model
  !! 
  !! @par alternative use
  !! supply N and P to maintain a certain target C:N:P of the labile pool
  !! 
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE update_vegsoil_bgc(tile, lctlib)

    USE mo_base_memory,      ONLY: t_slm_soil, t_slm_tile
    USE mo_veg_memory,       ONLY: t_slm_plant
    USE mo_base_constants,   ONLY: dt

    IMPLICIT NONE
    ! 0.1 InOut
    TYPE(t_lctlib_element),INTENT(in)        :: lctlib        !< land-cover-type library - parameters for one pft
    TYPE(t_slm_tile),INTENT(inout), TARGET   :: tile          !< one tile
    ! 0.2 Local
    TYPE(t_slm_soil),POINTER      :: soil               ! soil variables per tile
    TYPE(t_slm_plant),POINTER     :: veg                ! plant variables per tile
    REAL(wp) :: dcost                                   ! cost_n_uptake_root - cost_n_fixation (mol C / mol N)
    REAL(wp) :: cost_nacq                               ! total C cost of N acquisition (mol C / mol N)
    REAL(wp) :: f_nacq_bnf,f_nacq_nup                   ! fraction of N uptake as BFN or Nuptake

    !  bi-directional pointer
    veg   => tile%veg
    soil  => tile%soil

    !> 1.0 Calculate an instantaneous scalar of nutrient demand
    !!     taken as difference of the the weekly growth requirement of the nutrient to the stoichiometry of the labile pool
    !!
    !! this probably needs some thought - f_n_demand should not grow too large, but go to zero, if the plant has too much N/P 
    IF(veg%pool%labile%carbon > eps8) THEN
       veg%f_n_demand = MAX(veg%growth_req_n_tuptake_mavg - veg%pool%labile%nitrogen / veg%pool%labile%carbon, 0.0_wp) 
    ELSE
       veg%f_n_demand = 1.0_wp
    ENDIF
    IF(veg%pool%labile%nitrogen > eps8) THEN
       veg%f_p_demand = MAX(veg%growth_req_p_tuptake_mavg - veg%pool%labile%phosphorus / veg%pool%labile%nitrogen, 0.0_wp)
    ELSE
       veg%f_p_demand = 1.0_wp
    ENDIF

    !> 2.0 inorganic N uptake (first guess)
    !!
    !!
    SELECT CASE (veg_config%n_uptake_scheme)
      !> 2.1 Case 1  maintain prescribed C:N:P irrespective of nutrient uptake capacity == no nutrient limitation case 
      !!
      CASE ("prescribed")
         ! uptake required to maintain a presribed C:N ratio of the labile pool in µmol N / m2 / s
         veg%uptake_n = MAX(veg%pool%labile%carbon / lctlib%cn_leaf - & 
                                           veg%pool%labile%nitrogen, 0.0_wp) * 1.e6_wp / dt
         ! uptake required to maintain a presribed N:P ratio of the labile pool in µmol P / m2 / s
         veg%uptake_p = MAX(veg%pool%labile%nitrogen / lctlib%np_leaf + & 
                                              veg%uptake_n * dt / 1.e6 - &
                                              veg%pool%labile%phosphorus, 0.0_wp) * 1.e6_wp / dt

      !> 2.2 Case 2. Calculate inorganic N uptake according to Michaelis-Mnten Kinetics  
      !!
      CASE ("standard","optimal","resistance")
         ! OCN's MM Kinetics
         veg%uptake_n = lctlib%vmax_uptake_n * &
                 veg%demand_uptake_n_tuptake_mavg * & 
                 soil%ammonium * ( lctlib%soil_n_aff1 + & 
                   1._wp / ( lctlib%soil_n_aff2 + soil%ammonium )) * & 
                 veg%pool%fine_root%carbon
         veg%uptake_p = lctlib%vmax_uptake_p * & 
                 veg%demand_uptake_p_tuptake_mavg * & 
                 soil%phosphate * ( lctlib%soil_p_aff1 + & 
                   1._wp / ( lctlib%soil_p_aff2 + soil%phosphate )) * & 
                 veg%pool%fine_root%carbon

    ENDSELECT

    !> 3.0 N fixation
    !!
    !!
    SELECT CASE (veg_config%n_uptake_scheme)
      !> 3.1 Case 1: no N fixation
      !!
      CASE ("prescribed") 
        veg%n_fixation = 0.0_wp
      CASE("standard")
        IF(soil%ammonium < 0.5_wp) veg%n_fixation = lctlib%bnf_base 
      !> 3.2 Case 2: optimal N fixation
      !! BNF calculation with Michaelis Menten function (Rastetter, 2001) (MM)
      CASE ("optimal")
        IF (veg%cost_n_fixation_tuptake_mavg <= veg%cost_n_uptake_root) THEN 
           dcost  = veg%cost_n_uptake_root - veg%cost_n_fixation_tuptake_mavg
           !@BLARPP!! Don'T understand the term veg%growth_req_n_tuptake_mavg*veg%npp_tuptake_mavg here
           !veg%n_fixation = c_root * vmax_bnf * dcost / (k_bnf + dcost) * veg%growth_req_n_tuptake_mavg*veg%npp_tuptake_mavg
           veg%n_fixation = veg%pool%fine_root%carbon * & 
                              vmax_bnf * dcost / (k_bnf + dcost)
        ELSE
           veg%n_fixation = 0.0_wp
        ENDIF
      !> 3.3 Case 3: Ohm's law analogy 
      !! BNF calculation after FUN2.0 (Brzostek, 2014) (RA)
      CASE ("resistance")
         ! total cost according to Ohm's law: parallel resistors
         cost_nacq = 1._wp / (1._wp/veg%cost_n_uptake_root + 1._wp/veg%cost_n_fixation_tuptake_mavg) 
         IF (veg%cost_n_fixation_tuptake_mavg > 5._wp * veg%cost_n_uptake_root) THEN
            f_nacq_bnf = 0.0_wp
            f_nacq_nup = 1.0_wp
         ELSE 
            f_nacq_bnf = veg%cost_n_uptake_root/(veg%cost_n_uptake_root+veg%cost_n_fixation_tuptake_mavg)
            f_nacq_nup = veg%cost_n_fixation_tuptake_mavg/(veg%cost_n_uptake_root+veg%cost_n_fixation_tuptake_mavg)
         ENDIF
         veg%n_fixation = f_nacq_bnf * MAX(veg%growth_req_n_tuptake_mavg*veg%npp_tuptake_mavg,0.0_wp) 
         veg%uptake_n   = MIN(f_nacq_nup * MAX(veg%growth_req_n_tuptake_mavg*veg%npp_tuptake_mavg,0.0_wp), & 
                                veg%uptake_n)
    ENDSELECT 

    !> 3.0 Calculate potential N uptake per unit root
    !!
    veg%unit_uptake_n = lctlib%vmax_uptake_n * &
                               soil%ammonium * ( lctlib%soil_n_aff1 + & 
                               1._wp / ( lctlib%soil_n_aff2 + soil%ammonium ))
    veg%unit_uptake_p = lctlib%vmax_uptake_p * &
                               soil%phosphate * ( lctlib%soil_p_aff1 + & 
                               1._wp / ( lctlib%soil_p_aff2 + soil%phosphate ))

  END SUBROUTINE update_vegsoil_bgc

  !-----------------------------------------------------------------------------------------------------
  ! Main Task
  ! 
  !-----------------------------------------------------------------------------------------------------
  !> routine to update soil litter and mineral pools
  !! 
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE update_vegsoil_pool(tile, lctlib)

    USE mo_base_memory,      ONLY: t_slm_soil, t_slm_tile
    USE mo_base_constants,   ONLY: dt
    USE mo_veg_memory,       ONLY: t_slm_plant

    IMPLICIT NONE
    ! 0.1 InOut
    TYPE(t_lctlib_element),INTENT(in)         :: lctlib        !< land-cover-type library - parameters for one pft
    TYPE(t_slm_tile),INTENT(inout), TARGET    :: tile          !< one tile
    ! 0.2 Local
    TYPE(t_slm_soil),POINTER      :: soil               ! soil variables per tile
    TYPE(t_slm_plant),POINTER     :: veg              ! plant variables per tile
    REAl(wp) :: soil_depth = 1.0 ! soil depth, m

    !  bi-directional pointer
    veg   => tile%veg
    soil  => tile%soil


    !> 1 Calculate soil litter pool decomposition. Currently only for N
    !!    
    !! Needs temperature and moisture correction
    soil%loss%litter = 1.0_wp/tau_litter * soil%pool%litter
    
    !> 2 Update soil litter pool. Probably needs splitting into structural and non-structural litter
    !!
    soil%pool%litter = soil%pool%litter + veg%litterfall%leaf
    soil%pool%litter = soil%pool%litter + veg%litterfall%fine_root
    soil%pool%litter = soil%pool%litter + veg%litterfall%coarse_root
    soil%pool%litter = soil%pool%litter + veg%litterfall%heart_wood

    soil%pool%litter = soil%pool%litter - soil%loss%litter * dt
    
    !> 3 Update mineral pool
    !! 
!SZ BLARPP. check units
    IF(veg_config%n_uptake_scheme=='prescribed') THEN 
       soil%ammonium = 0.2_wp
    ELSE
       soil%ammonium = soil%ammonium   + soil%loss%litter%nitrogen*dt*  1.e3 / soil_depth * 1.e-3
       soil%ammonium = soil%ammonium   - veg%uptake_n * dt/1.e6 * 1.e3 / soil_depth * 1.e-3 
    ENDIF

    !soil%phosphate = soil%phosphate + soil%loss%litter%phosphorus*dt*  1.e3 / soil_depth * 1.e-3
    !soil%phosphate = soil%phosphate - veg%uptake_p * dt/1.e6 * 1.e3 / soil_depth * 1.e-3 

    END SUBROUTINE update_vegsoil_pool

  !-----------------------------------------------------------------------------------------------------
  ! Main Task
  ! 
  !-----------------------------------------------------------------------------------------------------
  !> update soil moisture and theta for the plant only model
  !! 
  !! two cases considered: \n a) dynamic water uptake calculation, \n b) prescribed values from namelist
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE update_vegsoil_physics(tile, lctlib)

    USE mo_jsb_process_class,     ONLY: ATM2LAND_, LAND2ATM_, SEB_, TURB_, HYDRO_, VEG_, HD_, RAD_, ASSIMI_, PHENO_
!    USE mo_jsb_process_class,     ONLY: ATM2LAND_, LAND2ATM_, SEB_, TURB_, SOIL_, HYDRO_, VEG_, HD_, RAD_, ASSIMI_, PHENO_ --> soil creates conflict

    USE mo_base_memory,           ONLY: t_slm_soil, t_slm_tile
    USE mo_veg_memory,            ONLY: t_slm_plant
    USE mo_atm2land_memory_class, ONLY: t_atm2land_memory
    USE mo_land2atm_memory_class, ONLY: t_land2atm_memory

    USE mo_base_constants,        ONLY: dt
    USE mo_atmland_util,          ONLY: calc_spec_humidity_sat
    USE mo_shared_constants,      ONLY: r_gas_dryair,LatentHeatVaporization,rhoh2o,grav

    IMPLICIT NONE
    ! 0.1 InOut
    TYPE(t_slm_tile),INTENT(inout), TARGET  :: tile     !< one tile
    TYPE(t_lctlib_element),INTENT(in)       :: lctlib   !< land-cover-type library - parameters for one pft
    ! 0.2 Local
    TYPE(t_slm_soil),POINTER                :: soil     ! soil variables per tile
    TYPE(t_slm_plant),POINTER               :: veg    ! plant variables per tile
    INTEGER  :: t
    REAL(wp) :: q_sat
    REAL(wp) :: evaporation
    ! 0.3 Memory allocation
    TYPE(t_atm2land_memory), POINTER :: atm2land__mem   ! dsl4jsb_Def_memory(ATM2LAND_)
    TYPE(t_land2atm_memory), POINTER :: land2atm__mem
    REAL(wp), POINTER                :: t_air
    REAL(wp), POINTER                :: press_srf
    REAL(wp), POINTER                :: q_air
    REAL(wp), POINTER                :: ga
    REAL(wp), POINTER                :: rain
    REAL(wp), POINTER                :: snow
    REAL(wp), POINTER                :: net_radiation
    REAL(wp), POINTER                :: latent_heat_flx
    REAL(wp), POINTER                :: sensible_heat_flx
    SELECT TYPE (mm => tile%mem(ATM2LAND_)%p)            ! dsl4jsb_Get_memory(ATM2LAND_)
      TYPE IS(t_atm2land_memory)
        atm2land__mem => mm
    END SELECT
    SELECT TYPE (mm => tile%mem(LAND2ATM_)%p)
      TYPE IS(t_land2atm_memory)
        land2atm__mem => mm
    END SELECT
    t_air       => atm2land__mem%t_air%ptr
    press_srf   => atm2land__mem%press_srf%ptr
    q_air       => atm2land__mem%q_air%ptr
    ga          => atm2land__mem%ga%ptr
    rain        => atm2land__mem%rain%ptr
    snow        => atm2land__mem%snow%ptr
    net_radiation      => land2atm__mem%net_radiation%ptr
    latent_heat_flx    => land2atm__mem%latent_heat_flx%ptr
    sensible_heat_flx  => land2atm__mem%sensible_heat_flx%ptr


    !  bi-directional pointer
    veg   => tile%veg
    soil  => tile%soil


    ! saturating specific humidity 
    q_sat = calc_spec_humidity_sat(t_air,press_srf)

    IF(veg%gc > eps8)THEN
       ! T = rho * (qsat - qair ) * Gtot in mm/s 
       veg%transpiration = press_srf / ( r_gas_dryair * t_air ) * & 
                  (q_sat - q_air) / (1._wp / veg%gc + 1._wp / ga )  
    ELSE 
       veg%transpiration = 0.0_wp 
    ENDIF
    ! E = rho * (qsat - qair) * theta * 1/ga_soil (taken as ga_soil = ga / 10)
    evaporation = press_srf / ( r_gas_dryair * t_air ) * & 
                  (q_sat - q_air) * soil%theta / 100.0_wp / ( 100._wp / ga)
    IF(net_radiation < eps8) evaporation = 0.0_wp

    ! calculate latent heat flux (W m-2) 
    latent_heat_flx = LatentHeatVaporization * (evaporation + veg%transpiration)
    ! diagnose sensible heat flux (shortwave compontent only)
    sensible_heat_flx = net_radiation - latent_heat_flx 

    IF(veg_config%flag_calc_soil_moisture) THEN ! case of dynamic water uptake calculation 

       ! update soil moisture dS = rain - transpiration - evaporation 
       soil%moisture = soil%moisture + (rain + snow & 
            - veg%transpiration - evaporation) * dt
       ! constrain to possible solutions
       soil%moisture = MAX(MIN(soil%moisture,veg_config%soil_awc_prescribe),0.0_wp)
       ! scale to relative water availability (%)
       soil%theta = soil%moisture / veg_config%soil_awc_prescribe * 100._wp

   ELSE 
       ! use prescribed value of theta from namelist
       soil%theta = veg_config%soil_theta_prescribe 
 
       ! calculate corresponding soil moisture value (mm)     
       soil%moisture = veg_config%soil_awc_prescribe * soil%theta / 100.0_wp

   ENDIF 

   ! calculate soil wtaer potential according to Campbell et al 1974 
   soil%water_potential = k_soil_phi1 * & 
         (( theta_wilt + 0.01_wp * soil%theta * (theta_sat - theta_wilt))/theta_sat)**k_soil_phi2
   veg%dphi = soil%water_potential - lctlib%phi_leaf_min - grav * rhoh2o * veg%height*1.e-6 

  END SUBROUTINE update_vegsoil_physics

END MODULE mo_veg_soil
