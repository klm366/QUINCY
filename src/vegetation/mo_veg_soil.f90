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

#include "dsl_jsbach.inc"

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
  SUBROUTINE update_vegsoil_bgc(tile)

    USE mo_jsb_class,             ONLY: Get_model
    USE mo_base_memory,           ONLY: t_slm_soil, t_slm_tile, t_model
    USE mo_base_constants,        ONLY: dtime
    USE mo_jsb_process_class,     ONLY: ATM2LAND_, LAND2ATM_, SEB_, TURB_, SOIL_, HYDRO_, VEG_, HD_, RAD_, ASSIMI_, PHENO_

    ! Use of process memories
    dsl4jsb_Use_memory(VEG_)        !  USE mo_veg_memory_class,      ONLY: t_veg_memory

    IMPLICIT NONE
    ! 0.1 InOut
    TYPE(t_slm_tile),INTENT(inout), TARGET   :: tile      !< one tile
    ! 0.2 Local
    TYPE(t_model),          POINTER        :: model
    TYPE(t_lctlib_element), POINTER        :: lctlib      !< land-cover-type library - parameter across pft's
    TYPE(t_slm_soil),       POINTER        :: soil        ! soil variables per tile
    REAL(wp) :: dcost                                     ! cost_n_uptake_root - cost_n_fixation (mol C / mol N)
    REAL(wp) :: cost_nacq                                 ! total C cost of N acquisition (mol C / mol N)
    REAL(wp) :: f_nacq_bnf,f_nacq_nup                     ! fraction of N uptake as BFN or Nuptake
    ! 0.3 Memory allocation
    ! Declare process configuration and memory Pointers
    dsl4jsb_Def_memory(VEG_)                          ! TYPE(t_veg_memory),      POINTER :: veg__mem
    TYPE(t_veg_memory), POINTER             :: veg        ! local pointer to veg__mem for SLM usability
    dsl4jsb_Get_memory(VEG_)
!    SELECT TYPE (mm => tile%mem(VEG_)%p) 
!      TYPE IS(t_veg_memory)
!        veg__mem => mm
!    END SELECT
    veg => veg__mem 

    model  => Get_model(tile%owner_model_id)
    lctlib => model%lctlib(tile%lcts%lib_id)

    !  bi-directional pointer
    soil  => tile%soil

    !> 1.0 Calculate an instantaneous scalar of nutrient demand
    !!     taken as difference of the the weekly growth requirement of the nutrient to the stoichiometry of the labile pool
    !!
    !! this probably needs some thought - f_n_demand should not grow too large, but go to zero, if the plant has too much N/P 
    IF(veg%pool%labile%carbon%ptr > eps8) THEN
       veg%f_n_demand%ptr = MAX(veg%growth_req_n_tuptake_mavg%ptr - veg%pool%labile%nitrogen%ptr / veg%pool%labile%carbon%ptr, 0.0_wp) 
    ELSE
       veg%f_n_demand%ptr = 1.0_wp
    ENDIF
    IF(veg%pool%labile%nitrogen%ptr > eps8) THEN
       veg%f_p_demand%ptr = MAX(veg%growth_req_p_tuptake_mavg%ptr - veg%pool%labile%phosphorus%ptr / veg%pool%labile%nitrogen%ptr, 0.0_wp)
    ELSE
       veg%f_p_demand%ptr = 1.0_wp
    ENDIF

    !> 2.0 inorganic N uptake (first guess)
    !!
    !!
    SELECT CASE (veg_config%n_uptake_scheme)
      !> 2.1 Case 1  maintain prescribed C:N:P irrespective of nutrient uptake capacity == no nutrient limitation case 
      !!
      CASE ("prescribed")
         ! uptake required to maintain a presribed C:N ratio of the labile pool in µmol N / m2 / s
         veg%uptake_n%ptr = MAX(veg%pool%labile%carbon%ptr / lctlib%cn_leaf - & 
                                           veg%pool%labile%nitrogen%ptr, 0.0_wp) * 1.e6_wp / dtime
         ! uptake required to maintain a presribed N:P ratio of the labile pool in µmol P / m2 / s
         veg%uptake_p%ptr = MAX(veg%pool%labile%nitrogen%ptr / lctlib%np_leaf + & 
                                              veg%uptake_n%ptr * dtime / 1.e6 - &
                                              veg%pool%labile%phosphorus%ptr, 0.0_wp) * 1.e6_wp / dtime

      !> 2.2 Case 2. Calculate inorganic N uptake according to Michaelis-Mnten Kinetics  
      !!
      CASE ("standard","optimal","resistance")
         ! OCN's MM Kinetics
         veg%uptake_n%ptr = lctlib%vmax_uptake_n * &
                 veg%demand_uptake_n_tuptake_mavg%ptr * & 
                 soil%ammonium * ( lctlib%soil_n_aff1 + & 
                   1._wp / ( lctlib%soil_n_aff2 + soil%ammonium )) * & 
                 veg%pool%fine_root%carbon%ptr
         veg%uptake_p%ptr = lctlib%vmax_uptake_p * & 
                 veg%demand_uptake_p_tuptake_mavg%ptr * & 
                 soil%phosphate * ( lctlib%soil_p_aff1 + & 
                   1._wp / ( lctlib%soil_p_aff2 + soil%phosphate )) * & 
                 veg%pool%fine_root%carbon%ptr

    ENDSELECT

    !> 3.0 N fixation
    !!
    !!
    SELECT CASE (veg_config%n_uptake_scheme)
      !> 3.1 Case 1: no N fixation
      !!
      CASE ("prescribed") 
        veg%n_fixation%ptr = 0.0_wp
      CASE("standard")
        IF(soil%ammonium < 0.5_wp) veg%n_fixation%ptr = lctlib%bnf_base 
      !> 3.2 Case 2: optimal N fixation
      !! BNF calculation with Michaelis Menten function (Rastetter, 2001) (MM)
      CASE ("optimal")
        IF (veg%cost_n_fixation_tuptake_mavg%ptr <= veg%cost_n_uptake_root%ptr) THEN 
           dcost  = veg%cost_n_uptake_root%ptr - veg%cost_n_fixation_tuptake_mavg%ptr
           !@BLARPP!! Don'T understand the term veg%growth_req_n_tuptake_mavg%ptr*veg%npp_tuptake_mavg%ptr here
           !veg%n_fixation%ptr = c_root * vmax_bnf * dcost / (k_bnf + dcost) * veg%growth_req_n_tuptake_mavg%ptr*veg%npp_tuptake_mavg%ptr
           veg%n_fixation%ptr = veg%pool%fine_root%carbon%ptr * & 
                              vmax_bnf * dcost / (k_bnf + dcost)
        ELSE
           veg%n_fixation%ptr = 0.0_wp
        ENDIF
      !> 3.3 Case 3: Ohm's law analogy 
      !! BNF calculation after FUN2.0 (Brzostek, 2014) (RA)
      CASE ("resistance")
         ! total cost according to Ohm's law: parallel resistors
         cost_nacq = 1._wp / (1._wp/veg%cost_n_uptake_root%ptr + 1._wp/veg%cost_n_fixation_tuptake_mavg%ptr) 
         IF (veg%cost_n_fixation_tuptake_mavg%ptr > 5._wp * veg%cost_n_uptake_root%ptr) THEN
            f_nacq_bnf = 0.0_wp
            f_nacq_nup = 1.0_wp
         ELSE 
            f_nacq_bnf = veg%cost_n_uptake_root%ptr/(veg%cost_n_uptake_root%ptr+veg%cost_n_fixation_tuptake_mavg%ptr)
            f_nacq_nup = veg%cost_n_fixation_tuptake_mavg%ptr/(veg%cost_n_uptake_root%ptr+veg%cost_n_fixation_tuptake_mavg%ptr)
         ENDIF
         veg%n_fixation%ptr = f_nacq_bnf * MAX(veg%growth_req_n_tuptake_mavg%ptr*veg%npp_tuptake_mavg%ptr,0.0_wp) 
         veg%uptake_n%ptr   = MIN(f_nacq_nup * MAX(veg%growth_req_n_tuptake_mavg%ptr*veg%npp_tuptake_mavg%ptr,0.0_wp), & 
                                veg%uptake_n%ptr)
    ENDSELECT 

    !> 3.0 Calculate potential N uptake per unit root
    !!
    veg%unit_uptake_n%ptr = lctlib%vmax_uptake_n * &
                               soil%ammonium * ( lctlib%soil_n_aff1 + & 
                               1._wp / ( lctlib%soil_n_aff2 + soil%ammonium ))
    veg%unit_uptake_p%ptr = lctlib%vmax_uptake_p * &
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
  SUBROUTINE update_vegsoil_pool(tile)

    USE mo_jsb_class,             ONLY: Get_model
    USE mo_base_memory,           ONLY: t_slm_soil, t_slm_tile, t_model
    USE mo_base_constants,        ONLY: dtime,one_day,one_year
    USE mo_jsb_process_class,     ONLY: ATM2LAND_, LAND2ATM_, SEB_, TURB_, SOIL_, HYDRO_, VEG_, HD_, RAD_, ASSIMI_, PHENO_

    ! Use of process memories
    dsl4jsb_Use_memory(VEG_)        !  USE mo_veg_memory_class,      ONLY: t_veg_memory

    IMPLICIT NONE
    ! 0.1 InOut
    TYPE(t_slm_tile),INTENT(inout), TARGET    :: tile          !< one tile
    ! 0.2 Local
    TYPE(t_model),          POINTER       :: model
    TYPE(t_lctlib_element), POINTER       :: lctlib     !< land-cover-type library - parameter across pft's
    TYPE(t_slm_soil),       POINTER       :: soil       ! soil variables per tile
    REAl(wp) :: soil_depth = 1.0 ! soil depth, m
    ! 0.3 Memory allocation
    ! Declare process configuration and memory Pointers
    dsl4jsb_Def_memory(VEG_)                          ! TYPE(t_veg_memory),      POINTER :: veg__mem
    TYPE(t_veg_memory), POINTER             :: veg        ! local pointer to veg__mem for SLM usability
    dsl4jsb_Get_memory(VEG_)
!    SELECT TYPE (mm => tile%mem(VEG_)%p) 
!      TYPE IS(t_veg_memory)
!        veg__mem => mm
!    END SELECT
    veg => veg__mem 

    model  => Get_model(tile%owner_model_id)
    lctlib => model%lctlib(tile%lcts%lib_id)

    !  bi-directional pointer
    soil  => tile%soil


    !> 1 Calculate soil litter pool decomposition. Currently only for N
    !!    
    !! Needs temperature and moisture correction
    soil%loss%litter = 1.0_wp/tau_litter/one_day/one_year * soil%pool%litter * dtime
    
    !> 2 Update soil litter pool. Probably needs splitting into structural and non-structural litter
    !!
!! not working because TYPE(t_element)  :: litter has not implemented any operator overload for t_jsb_var_real1d
! *slmdev
!    soil%pool%litter = soil%pool%litter + veg%litterfall%leaf
!    soil%pool%litter = soil%pool%litter + veg%litterfall%fine_root
!    soil%pool%litter = soil%pool%litter + veg%litterfall%coarse_root
!    soil%pool%litter = soil%pool%litter + veg%litterfall%heart_wood
!! writing the code explicit instead and making it a operator overload with a scalar by using ...%ptr
    soil%pool%litter%carbon = soil%pool%litter%carbon + veg%litterfall%leaf%carbon%ptr
    soil%pool%litter%nitrogen = soil%pool%litter%nitrogen + veg%litterfall%leaf%nitrogen%ptr
    soil%pool%litter%phosphorus = soil%pool%litter%phosphorus + veg%litterfall%leaf%phosphorus%ptr
    soil%pool%litter%carbon13 = soil%pool%litter%carbon13 + veg%litterfall%leaf%carbon13%ptr
    soil%pool%litter%carbon14 = soil%pool%litter%carbon14 + veg%litterfall%leaf%carbon14%ptr
    soil%pool%litter%nitrogen15 = soil%pool%litter%nitrogen15 + veg%litterfall%leaf%nitrogen15%ptr
 
    soil%pool%litter%carbon = soil%pool%litter%carbon + veg%litterfall%fine_root%carbon%ptr
    soil%pool%litter%nitrogen = soil%pool%litter%nitrogen + veg%litterfall%fine_root%nitrogen%ptr
    soil%pool%litter%phosphorus = soil%pool%litter%phosphorus + veg%litterfall%fine_root%phosphorus%ptr
    soil%pool%litter%carbon13 = soil%pool%litter%carbon13 + veg%litterfall%fine_root%carbon13%ptr
    soil%pool%litter%carbon14 = soil%pool%litter%carbon14 + veg%litterfall%fine_root%carbon14%ptr
    soil%pool%litter%nitrogen15 = soil%pool%litter%nitrogen15 + veg%litterfall%fine_root%nitrogen15%ptr
 
    soil%pool%litter%carbon = soil%pool%litter%carbon + veg%litterfall%coarse_root%carbon%ptr
    soil%pool%litter%nitrogen = soil%pool%litter%nitrogen + veg%litterfall%coarse_root%nitrogen%ptr
    soil%pool%litter%phosphorus = soil%pool%litter%phosphorus + veg%litterfall%coarse_root%phosphorus%ptr
    soil%pool%litter%carbon13 = soil%pool%litter%carbon13 + veg%litterfall%coarse_root%carbon13%ptr
    soil%pool%litter%carbon14 = soil%pool%litter%carbon14 + veg%litterfall%coarse_root%carbon14%ptr
    soil%pool%litter%nitrogen15 = soil%pool%litter%nitrogen15 + veg%litterfall%coarse_root%nitrogen15%ptr
 
    soil%pool%litter%carbon = soil%pool%litter%carbon + veg%litterfall%sap_wood%carbon%ptr
    soil%pool%litter%nitrogen = soil%pool%litter%nitrogen + veg%litterfall%sap_wood%nitrogen%ptr
    soil%pool%litter%phosphorus = soil%pool%litter%phosphorus + veg%litterfall%sap_wood%phosphorus%ptr
    soil%pool%litter%carbon13 = soil%pool%litter%carbon13 + veg%litterfall%sap_wood%carbon13%ptr
    soil%pool%litter%carbon14 = soil%pool%litter%carbon14 + veg%litterfall%sap_wood%carbon14%ptr
    soil%pool%litter%nitrogen15 = soil%pool%litter%nitrogen15 + veg%litterfall%sap_wood%nitrogen15%ptr

    soil%pool%litter%carbon = soil%pool%litter%carbon + veg%litterfall%heart_wood%carbon%ptr
    soil%pool%litter%nitrogen = soil%pool%litter%nitrogen + veg%litterfall%heart_wood%nitrogen%ptr
    soil%pool%litter%phosphorus = soil%pool%litter%phosphorus + veg%litterfall%heart_wood%phosphorus%ptr
    soil%pool%litter%carbon13 = soil%pool%litter%carbon13 + veg%litterfall%heart_wood%carbon13%ptr
    soil%pool%litter%carbon14 = soil%pool%litter%carbon14 + veg%litterfall%heart_wood%carbon14%ptr
    soil%pool%litter%nitrogen15 = soil%pool%litter%nitrogen15 + veg%litterfall%heart_wood%nitrogen15%ptr

    soil%pool%litter%carbon = soil%pool%litter%carbon + veg%litterfall%fruit%carbon%ptr
    soil%pool%litter%nitrogen = soil%pool%litter%nitrogen + veg%litterfall%fruit%nitrogen%ptr
    soil%pool%litter%phosphorus = soil%pool%litter%phosphorus + veg%litterfall%fruit%phosphorus%ptr
    soil%pool%litter%carbon13 = soil%pool%litter%carbon13 + veg%litterfall%fruit%carbon13%ptr
    soil%pool%litter%carbon14 = soil%pool%litter%carbon14 + veg%litterfall%fruit%carbon14%ptr
    soil%pool%litter%nitrogen15 = soil%pool%litter%nitrogen15 + veg%litterfall%fruit%nitrogen15%ptr

    soil%pool%litter%carbon = soil%pool%litter%carbon + veg%litterfall%seed_bed%carbon%ptr
    soil%pool%litter%nitrogen = soil%pool%litter%nitrogen + veg%litterfall%seed_bed%nitrogen%ptr
    soil%pool%litter%phosphorus = soil%pool%litter%phosphorus + veg%litterfall%seed_bed%phosphorus%ptr
    soil%pool%litter%carbon13 = soil%pool%litter%carbon13 + veg%litterfall%seed_bed%carbon13%ptr
    soil%pool%litter%carbon14 = soil%pool%litter%carbon14 + veg%litterfall%seed_bed%carbon14%ptr
    soil%pool%litter%nitrogen15 = soil%pool%litter%nitrogen15 + veg%litterfall%seed_bed%nitrogen15%ptr

    soil%pool%litter%carbon = soil%pool%litter%carbon + veg%litterfall%labile%carbon%ptr
    soil%pool%litter%nitrogen = soil%pool%litter%nitrogen + veg%litterfall%labile%nitrogen%ptr
    soil%pool%litter%phosphorus = soil%pool%litter%phosphorus + veg%litterfall%labile%phosphorus%ptr
    soil%pool%litter%carbon13 = soil%pool%litter%carbon13 + veg%litterfall%labile%carbon13%ptr
    soil%pool%litter%carbon14 = soil%pool%litter%carbon14 + veg%litterfall%labile%carbon14%ptr
    soil%pool%litter%nitrogen15 = soil%pool%litter%nitrogen15 + veg%litterfall%labile%nitrogen15%ptr

    soil%pool%litter%carbon = soil%pool%litter%carbon + veg%litterfall%reserve%carbon%ptr
    soil%pool%litter%nitrogen = soil%pool%litter%nitrogen + veg%litterfall%reserve%nitrogen%ptr
    soil%pool%litter%phosphorus = soil%pool%litter%phosphorus + veg%litterfall%reserve%phosphorus%ptr
    soil%pool%litter%carbon13 = soil%pool%litter%carbon13 + veg%litterfall%reserve%carbon13%ptr
    soil%pool%litter%carbon14 = soil%pool%litter%carbon14 + veg%litterfall%reserve%carbon14%ptr
    soil%pool%litter%nitrogen15 = soil%pool%litter%nitrogen15 + veg%litterfall%reserve%nitrogen15%ptr

    soil%pool%litter = soil%pool%litter - soil%loss%litter 
    
    !> 3 Update mineral pool
    !! 
!SZ BLARPP. check units
    IF(veg_config%n_uptake_scheme=='prescribed') THEN 
       soil%ammonium = 0.2_wp
    ELSE
       soil%ammonium = soil%ammonium   + soil%loss%litter%nitrogen *  1.e3_wp / soil_depth * 1.e-3_wp
       soil%ammonium = soil%ammonium   - veg%uptake_n%ptr * dtime/1.e6_wp * 1.e3_wp / soil_depth * 1.e-3_wp 
    ENDIF

    !soil%phosphate = soil%phosphate + soil%loss%litter%phosphorus*dtime*  1.e3_wp / soil_depth * 1.e-3_wp
    !soil%phosphate = soil%phosphate - veg%uptake_p%ptr * dtime/1.e6_wp * 1.e3_wp / soil_depth * 1.e-3_wp 

    END SUBROUTINE update_vegsoil_pool

  !-----------------------------------------------------------------------------------------------------
  ! Main Task
  ! 
  !-----------------------------------------------------------------------------------------------------
  !> update soil moisture and theta for the plant only model
  !! 
  !! two cases considered: \n a) dynamic water uptake calculation, \n b) prescribed values from namelist
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE update_vegsoil_physics(tile)

    USE mo_jsb_class,             ONLY: Get_model
    USE mo_jsb_process_class,     ONLY: ATM2LAND_, LAND2ATM_, SEB_, TURB_, SOIL_, HYDRO_, VEG_, HD_, RAD_, ASSIMI_, PHENO_
    USE mo_base_memory,           ONLY: t_slm_soil, t_slm_tile, t_model
    USE mo_base_constants,        ONLY: dtime
    USE mo_atmland_util,          ONLY: calc_spec_humidity_sat
    USE mo_shared_constants,      ONLY: r_gas_dryair, LatentHeatVaporization, rhoh2o, grav

    ! Use of process memories
    dsl4jsb_Use_memory(ATM2LAND_)   !  USE mo_atm2land_memory_class, ONLY: t_atm2land_memory
    dsl4jsb_Use_memory(LAND2ATM_)   !  USE mo_land2atm_memory_class, ONLY: t_land2atm_memory
    dsl4jsb_Use_memory(VEG_)        !  USE mo_veg_memory_class,      ONLY: t_veg_memory

    IMPLICIT NONE
    ! 0.1 InOut
    TYPE(t_slm_tile),INTENT(inout), TARGET  :: tile     !< one tile
    ! 0.2 Local
    TYPE(t_model),          POINTER       :: model
    TYPE(t_lctlib_element), POINTER       :: lctlib     !< land-cover-type library - parameter across pft's
    TYPE(t_slm_soil),       POINTER       :: soil       ! soil variables per tile
    INTEGER  :: t
    REAL(wp) :: q_sat
    REAL(wp) :: evaporation
    ! 0.3 Memory allocation
    ! Declare process configuration and memory Pointers
    dsl4jsb_Def_memory(ATM2LAND_)                   !     TYPE(t_atm2land_memory), POINTER :: atm2land__mem
    dsl4jsb_Def_memory(LAND2ATM_)                   !     TYPE(t_land2atm_memory), POINTER :: land2atm__mem
    dsl4jsb_Def_memory(VEG_)                        !     TYPE(t_veg_memory),      POINTER :: veg__mem
    TYPE(t_veg_memory),      POINTER :: veg             ! local pointer to veg__mem for SLM usability
    REAL(wp), POINTER                :: t_air
    REAL(wp), POINTER                :: press_srf
    REAL(wp), POINTER                :: q_air
    REAL(wp), POINTER                :: ga
    REAL(wp), POINTER                :: rain
    REAL(wp), POINTER                :: snow
    REAL(wp), POINTER                :: net_radiation
    REAL(wp), POINTER                :: latent_heat_flx
    REAL(wp), POINTER                :: sensible_heat_flx
    dsl4jsb_Get_memory(ATM2LAND_)
    dsl4jsb_Get_memory(LAND2ATM_)
    dsl4jsb_Get_memory(VEG_)
!    SELECT TYPE (mm => tile%mem(ATM2LAND_)%p)
!      TYPE IS(t_atm2land_memory)
!        atm2land__mem => mm
!    END SELECT
!    SELECT TYPE (mm => tile%mem(LAND2ATM_)%p)
!      TYPE IS(t_land2atm_memory)
!        land2atm__mem => mm
!    END SELECT
!    SELECT TYPE (mm => tile%mem(VEG_)%p)
!      TYPE IS(t_veg_memory)
!        veg__mem => mm
!    END SELECT
    veg => veg__mem 
    t_air       => atm2land__mem%t_air%ptr
    press_srf   => atm2land__mem%press_srf%ptr
    q_air       => atm2land__mem%q_air%ptr
    ga          => atm2land__mem%ga%ptr
    rain        => atm2land__mem%rain%ptr
    snow        => atm2land__mem%snow%ptr
    net_radiation      => land2atm__mem%net_radiation%ptr
    latent_heat_flx    => land2atm__mem%latent_heat_flx%ptr
    sensible_heat_flx  => land2atm__mem%sensible_heat_flx%ptr

    model  => Get_model(tile%owner_model_id)
    lctlib => model%lctlib(tile%lcts%lib_id)

    !  bi-directional pointer
    soil  => tile%soil


    ! saturating specific humidity 
    q_sat = calc_spec_humidity_sat(t_air,press_srf)

    IF(veg%gc%ptr > eps8)THEN
       ! T = rho * (qsat - qair ) * Gtot in mm/s 
       veg%transpiration%ptr = press_srf / ( r_gas_dryair * t_air ) * & 
                  (q_sat - q_air) / (1._wp / veg%gc%ptr + 1._wp / ga )  
    ELSE 
       veg%transpiration%ptr = 0.0_wp 
    ENDIF
    ! E = rho * (qsat - qair) * theta * 1/ga_soil (taken as ga_soil = ga / 10)
    evaporation = press_srf / ( r_gas_dryair * t_air ) * & 
                  (q_sat - q_air) * soil%theta / 100.0_wp / ( 100._wp / ga)
    IF(net_radiation < eps8) evaporation = 0.0_wp

    ! calculate latent heat flux (W m-2) 
    latent_heat_flx = LatentHeatVaporization * (evaporation + veg%transpiration%ptr)
    ! diagnose sensible heat flux (shortwave compontent only)
    sensible_heat_flx = net_radiation - latent_heat_flx 

    IF(veg_config%flag_calc_soil_moisture) THEN ! case of dynamic water uptake calculation 

       ! update soil moisture dS = rain - transpiration - evaporation 
       soil%moisture = soil%moisture + (rain + snow & 
            - veg%transpiration%ptr - evaporation) * dtime
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
   veg%dphi%ptr = soil%water_potential - lctlib%phi_leaf_min - grav * rhoh2o * veg%height%ptr * 1.e-6_wp

  END SUBROUTINE update_vegsoil_physics

END MODULE mo_veg_soil
