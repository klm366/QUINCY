!>
!!  @par Copyright
!!    QUINCY terrestrial biosphere model with coupled carbon, nitrogen, and phosphorus cycles
!!    Copyright (C) 2019 Max Planck Institute for Biogeochemistry (MPI-BGC), Jena, Germany 
!!    Contributors: Caldararu S., Engel J., Kern M., Pallandt M., Schnur R., Thum T., Yu L., Zaehle S.
!!
!!    This file is part of the QUINCY model.
!!    The QUINCY model is free software: you can redistribute it and/or modify
!!    it under the terms of the GNU General Public License version 3,
!!    as published by the Free Software Foundation.
!!
!!    The QUINCY model is distributed in the hope that it will be useful,
!!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!!    GNU General Public License for more details.
!!
!!    You should have received a copy of the GNU General Public License
!!    along with the QUINCY model.  If not, see <https://www.gnu.org/licenses/>.
!! 
!! @brief Contains the memory class for the vegetation processes
!!
!! @author
!!  Jan Engel, MPI-BGC Jena (based on the JSBACH4 version from Rainer Schneck, MPI-M Hamburg)
!!
!! @par Revision History
!! First version by Jan Engel (20160830)
!!
MODULE mo_veg_memory_class

  USE mo_kind,                   ONLY: wp
  USE mo_exception,              ONLY: finish, message
  USE mo_util,                   ONLY: One_of
  USE mo_jsb_utils_iface,        ONLY: assign_if_present_allocatable
  USE mo_jsb_impl_constants,     ONLY: false
  
  USE mo_jsb_model_class,        ONLY: t_jsb_model
  USE mo_jsb_class,              ONLY: Get_model
  USE mo_jsb_memory_class,       ONLY: t_jsb_memory
  USE mo_jsb_var_class,          ONLY: t_jsb_var, t_jsb_var_real2d, t_jsb_var_real3d
  USE mo_jsb_pool_class,         ONLY: t_jsb_pool, t_jsb_pool_elem2d, ELEM_C, ELEM_N, ELEM_P, ELEM_C13, ELEM_C14, ELEM_N15
  USE mo_jsb_lct_class,          ONLY: LAND_TYPE, VEG_TYPE
  USE mo_jsb_subset,             ONLY: t_subset, ON_DOMAIN, ON_CHUNK

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: t_veg_memory, max_no_of_vars
  PUBLIC :: t_veg_main_pool, t_veg_tissues_pool, t_veg_elements_pool
 
  INTEGER, PARAMETER :: max_no_of_vars = 530    !> TYPE t_jsb_var | 'normal' var: ~ 220 | 5 pools: 270 var | 4 tissues: 24

  !-----------------------------------------------------------------------------------------------------
  !> Type definition for veg_elements_pool
  !! pool with only variables
  !! last pool in line
  !-----------------------------------------------------------------------------------------------------
  TYPE, EXTENDS(t_jsb_pool_elem2d) :: t_veg_elements_pool
    !
    ! element variable pointers are declared in mo_jsb_pool_class:t_jsb_pool_elem2d (see also mo_veg_memory_class:all_elements(:))
    !
    !! var for operator overload
    ! 1d ON_CHUNK
    REAL(wp), ALLOCATABLE   :: op_ov_hlp_1d_c(:)
    REAL(wp), ALLOCATABLE   :: op_ov_hlp_1d_n(:)
    REAL(wp), ALLOCATABLE   :: op_ov_hlp_1d_p(:)
    REAL(wp), ALLOCATABLE   :: op_ov_hlp_1d_c13(:)
    REAL(wp), ALLOCATABLE   :: op_ov_hlp_1d_c14(:)
    REAL(wp), ALLOCATABLE   :: op_ov_hlp_1d_n15(:)
    ! 2d ON_DOMAIN
    REAL(wp), ALLOCATABLE   :: op_ov_hlp_2d_c(:,:)
    REAL(wp), ALLOCATABLE   :: op_ov_hlp_2d_n(:,:)
    REAL(wp), ALLOCATABLE   :: op_ov_hlp_2d_p(:,:)
    REAL(wp), ALLOCATABLE   :: op_ov_hlp_2d_c13(:,:)
    REAL(wp), ALLOCATABLE   :: op_ov_hlp_2d_c14(:,:)
    REAL(wp), ALLOCATABLE   :: op_ov_hlp_2d_n15(:,:)
  CONTAINS
    !! infrastructure
    PROCEDURE :: Get_element_name_by_id => t_veg_elements_pool_get_element_name_by_id
    PROCEDURE :: Init                   => t_veg_elements_pool_init
    PROCEDURE :: Init_op_ov_hlp         => t_veg_elements_pool_op_ov_hlp_init             ! init hlp object for operator overload
    !! operator overload  
    !! [def: field1d(nc) | field2d(nc,nblks)]
    ! +
    PROCEDURE                 :: Add_veg_elem_p_and_scalar            ! elements_pool + scalar        | ON_CHUNK & ON_DOMAIN
    PROCEDURE, PASS(struct1)  :: Add_scalar_and_veg_elem_p            ! scalar        + elements_pool | ON_CHUNK & ON_DOMAIN
    PROCEDURE                 :: Add_veg_elem_p_and_field1d           ! elements_pool + field1d       | ON_CHUNK
    PROCEDURE, PASS(struct1)  :: Add_field1d_and_veg_elem_p           ! field1d       + elements_pool | ON_CHUNK
    PROCEDURE                 :: Add_veg_elem_p_and_field2d           ! elements_pool + field2d       | ON_CHUNK & ON_DOMAIN | field2d(:,1) & field2d(:,:)
    PROCEDURE                 :: Add_veg_elem_p_and_veg_elem_p        ! elements_pool + elements_pool | ON_CHUNK & ON_DOMAIN
    GENERIC                   :: OPERATOR(+) => Add_veg_elem_p_and_scalar,  Add_scalar_and_veg_elem_p, &
      &                                         Add_veg_elem_p_and_field1d, Add_field1d_and_veg_elem_p, &
      &                                         Add_veg_elem_p_and_field2d, &
      &                                         Add_veg_elem_p_and_veg_elem_p
    ! -
    PROCEDURE                 :: Subtract_veg_elem_p_minus_scalar     ! elements_pool - scalar        | ON_CHUNK & ON_DOMAIN
    PROCEDURE                 :: Subtract_veg_elem_p_minus_field1d    ! elements_pool - field1d       | ON_CHUNK
    PROCEDURE                 :: Subtract_veg_elem_p_minus_field2d    ! elements_pool - field2d       | ON_CHUNK & ON_DOMAIN | field2d(:,1) & field2d(:,:)
    PROCEDURE                 :: Subtract_veg_elem_p_minus_veg_elem_p ! elements_pool - elements_pool | ON_CHUNK & ON_DOMAIN
    GENERIC                   :: OPERATOR(-) => Subtract_veg_elem_p_minus_scalar, &
      &                                         Subtract_veg_elem_p_minus_field1d, &
      &                                         Subtract_veg_elem_p_minus_field2d, &
      &                                         Subtract_veg_elem_p_minus_veg_elem_p
    ! *
    PROCEDURE                 :: Multiply_veg_elem_p_with_scalar      ! elements_pool * scalar        | ON_CHUNK & ON_DOMAIN
    PROCEDURE, PASS(struct1)  :: Multiply_scalar_with_veg_elem_p      ! scalar        * elements_pool | ON_CHUNK & ON_DOMAIN
    PROCEDURE                 :: Multiply_veg_elem_p_with_field1d     ! elements_pool * field1d       | ON_CHUNK
    PROCEDURE, PASS(struct1)  :: Multiply_field1d_with_veg_elem_p     ! field1d       * elements_pool | ON_CHUNK
    PROCEDURE                 :: Multiply_veg_elem_p_with_field2d     ! elements_pool * field2d       | ON_CHUNK & ON_DOMAIN | field2d(:,1) & field2d(:,:)
    PROCEDURE, PASS(struct1)  :: Multiply_field2d_with_veg_elem_p     ! field2d       * elements_pool | ON_CHUNK & ON_DOMAIN | field2d(:,1) & field2d(:,:)
    PROCEDURE                 :: Multiply_veg_elem_p_with_veg_elem_p  ! elements_pool * elements_pool | ON_CHUNK & ON_DOMAIN
    GENERIC                   :: OPERATOR(*) => Multiply_veg_elem_p_with_scalar,  Multiply_scalar_with_veg_elem_p, &
      &                                         Multiply_veg_elem_p_with_field1d, Multiply_field1d_with_veg_elem_p, &
      &                                         Multiply_veg_elem_p_with_field2d, Multiply_field2d_with_veg_elem_p, &
      &                                         Multiply_veg_elem_p_with_veg_elem_p
    ! /
    PROCEDURE                 :: Div_veg_elem_p_by_scalar             ! elements_pool / scalar        | ON_CHUNK & ON_DOMAIN
    PROCEDURE                 :: Div_veg_elem_p_by_field1d            ! elements_pool / field1d       | ON_CHUNK
    PROCEDURE                 :: Div_veg_elem_p_by_field2d            ! elements_pool / field2d       | ON_CHUNK & ON_DOMAIN | field2d(:,1) & field2d(:,:)
    PROCEDURE                 :: Div_veg_elem_p_by_veg_elem_p         ! elements_pool / elements_pool | ON_CHUNK & ON_DOMAIN
    GENERIC                   :: OPERATOR(/) => Div_veg_elem_p_by_scalar, &
      &                                         Div_veg_elem_p_by_field1d, &
      &                                         Div_veg_elem_p_by_field2d, &
      &                                         Div_veg_elem_p_by_veg_elem_p
    ! =
    PROCEDURE                 :: Assign_scalar_to_veg_elem_p          ! elements_pool = scalar        | ON_CHUNK & ON_DOMAIN
    PROCEDURE                 :: Assign_field1d_to_veg_elem_p         ! elements_pool = field1d(:)    | ON_CHUNK
    PROCEDURE                 :: Assign_field2d_to_veg_elem_p         ! elements_pool = field2d(:,:)  | ON_CHUNK & ON_DOMAIN | field2d(:,1) & field2d(:,:)
    PROCEDURE                 :: Assign_veg_elem_p_to_veg_elem_p      ! elements_pool = elements_pool | ON_CHUNK & ON_DOMAIN | only numeric values of elements are assigned
    GENERIC                   :: ASSIGNMENT(=) => Assign_scalar_to_veg_elem_p, &
      &                                           Assign_field1d_to_veg_elem_p, &
      &                                           Assign_field2d_to_veg_elem_p, &
      &                                           Assign_veg_elem_p_to_veg_elem_p
    !! finalize
    FINAL     :: finalize_veg_elements_pool
  END TYPE t_veg_elements_pool

  !-----------------------------------------------------------------------------------------------------
  !> Type definition for t_veg_tissues_pool
  !! vegetation pool of tissues
  !-----------------------------------------------------------------------------------------------------
  TYPE, EXTENDS(t_jsb_pool_elem2d) :: t_veg_tissues_pool
    !
    ! element variable pointers are declared in mo_jsb_pool_class:t_jsb_pool_elem2d (see also mo_veg_memory_class:all_elements(:))
    !
    !! var for operator overload
    ! 1d ON_CHUNK
    REAL(wp), ALLOCATABLE   :: op_ov_hlp_1d_c(:)
    REAL(wp), ALLOCATABLE   :: op_ov_hlp_1d_n(:)
    REAL(wp), ALLOCATABLE   :: op_ov_hlp_1d_p(:)
    REAL(wp), ALLOCATABLE   :: op_ov_hlp_1d_c13(:)
    REAL(wp), ALLOCATABLE   :: op_ov_hlp_1d_c14(:)
    REAL(wp), ALLOCATABLE   :: op_ov_hlp_1d_n15(:)
    ! 2d ON_DOMAIN
    REAL(wp), ALLOCATABLE   :: op_ov_hlp_2d_c(:,:)
    REAL(wp), ALLOCATABLE   :: op_ov_hlp_2d_n(:,:)
    REAL(wp), ALLOCATABLE   :: op_ov_hlp_2d_p(:,:)
    REAL(wp), ALLOCATABLE   :: op_ov_hlp_2d_c13(:,:)
    REAL(wp), ALLOCATABLE   :: op_ov_hlp_2d_c14(:,:)
    REAL(wp), ALLOCATABLE   :: op_ov_hlp_2d_n15(:,:)
    !
    ! element pools
    !
    TYPE(t_veg_elements_pool), POINTER :: leaf       
    TYPE(t_veg_elements_pool), POINTER :: fine_root  
    TYPE(t_veg_elements_pool), POINTER :: coarse_root
    TYPE(t_veg_elements_pool), POINTER :: sap_wood   
    TYPE(t_veg_elements_pool), POINTER :: heart_wood 
    TYPE(t_veg_elements_pool), POINTER :: labile     
    TYPE(t_veg_elements_pool), POINTER :: reserve    
    TYPE(t_veg_elements_pool), POINTER :: fruit      
    TYPE(t_veg_elements_pool), POINTER :: seed_bed   
  CONTAINS
    !! infrastructure
    PROCEDURE :: Get_element_name_by_id   => t_veg_tissues_pool_get_element_name_by_id
    PROCEDURE :: Init                     => t_veg_tissues_pool_init                      ! pool_object created with FUNCTION Create_veg_tissues_pool()
    PROCEDURE :: Init_op_ov_hlp           => t_veg_tissues_pool_op_ov_hlp_create_init     ! create and init hlp object for operator overload
    !! operator overload  
    !! [def: field1d(nc) | field2d(nc,nblks)]
    ! +
    PROCEDURE                 :: Add_veg_tiss_p_and_scalar            ! tissues_pool + scalar       | ON_CHUNK & ON_DOMAIN
    PROCEDURE, PASS(struct1)  :: Add_scalar_and_veg_tiss_p            ! scalar       + tissues_pool | ON_CHUNK & ON_DOMAIN
    PROCEDURE                 :: Add_veg_tiss_p_and_field1d           ! tissues_pool + field1d      | ON_CHUNK
    PROCEDURE, PASS(struct1)  :: Add_field1d_and_veg_tiss_p           ! field1d      + tissues_pool | ON_CHUNK
    PROCEDURE                 :: Add_veg_tiss_p_and_field2d           ! tissues_pool + field2d      | ON_CHUNK & ON_DOMAIN | field2d(:,1) & field2d(:,:)
    PROCEDURE                 :: Add_veg_tiss_p_and_veg_tiss_p        ! tissues_pool + tissues_pool | ON_CHUNK & ON_DOMAIN
    GENERIC                   :: OPERATOR(+) => Add_veg_tiss_p_and_scalar,  Add_scalar_and_veg_tiss_p, &
      &                                         Add_veg_tiss_p_and_field1d, Add_field1d_and_veg_tiss_p, &
      &                                         Add_veg_tiss_p_and_field2d, &
      &                                         Add_veg_tiss_p_and_veg_tiss_p
    ! -
    PROCEDURE                 :: Subtract_veg_tiss_p_minus_scalar     ! tissues_pool - scalar       | ON_CHUNK & ON_DOMAIN
    PROCEDURE                 :: Subtract_veg_tiss_p_minus_field1d    ! tissues_pool - field1d      | ON_CHUNK
    PROCEDURE                 :: Subtract_veg_tiss_p_minus_field2d    ! tissues_pool - field2d      | ON_CHUNK & ON_DOMAIN | field2d(:,1) & field2d(:,:)
    PROCEDURE                 :: Subtract_veg_tiss_p_minus_veg_tiss_p ! tissues_pool - tissues_pool | ON_CHUNK & ON_DOMAIN
    GENERIC                   :: OPERATOR(-) => Subtract_veg_tiss_p_minus_scalar, &
      &                                         Subtract_veg_tiss_p_minus_field1d, &
      &                                         Subtract_veg_tiss_p_minus_field2d, &
      &                                         Subtract_veg_tiss_p_minus_veg_tiss_p
    ! *
    PROCEDURE                 :: Multiply_veg_tiss_p_with_scalar      ! tissues_pool * scalar       | ON_CHUNK & ON_DOMAIN
    PROCEDURE, PASS(struct1)  :: Multiply_scalar_with_veg_tiss_p      ! scalar       * tissues_pool | ON_CHUNK & ON_DOMAIN
    PROCEDURE                 :: Multiply_veg_tiss_p_with_field1d     ! tissues_pool * field1d      | ON_CHUNK
    PROCEDURE, PASS(struct1)  :: Multiply_field1d_with_veg_tiss_p     ! field1d      * tissues_pool | ON_CHUNK
    PROCEDURE                 :: Multiply_veg_tiss_p_with_field2d     ! tissues_pool * field2d      | ON_CHUNK & ON_DOMAIN | field2d(:,1) & field2d(:,:)
    PROCEDURE, PASS(struct1)  :: Multiply_field2d_with_veg_tiss_p     ! field2d      * tissues_pool | ON_CHUNK & ON_DOMAIN | field2d(:,1) & field2d(:,:)
    PROCEDURE                 :: Multiply_veg_tiss_p_with_veg_tiss_p  ! tissues_pool * tissues_pool | ON_CHUNK & ON_DOMAIN
    GENERIC                   :: OPERATOR(*) => Multiply_veg_tiss_p_with_scalar,  Multiply_scalar_with_veg_tiss_p, &
      &                                         Multiply_veg_tiss_p_with_field1d, Multiply_field1d_with_veg_tiss_p, &
      &                                         Multiply_veg_tiss_p_with_field2d, Multiply_field2d_with_veg_tiss_p, &
      &                                         Multiply_veg_tiss_p_with_veg_tiss_p
    ! /
    PROCEDURE                 :: Div_veg_tiss_p_by_scalar             ! tissues_pool / scalar       | ON_CHUNK & ON_DOMAIN
    PROCEDURE                 :: Div_veg_tiss_p_by_field1d            ! tissues_pool / field1d      | ON_CHUNK
    PROCEDURE                 :: Div_veg_tiss_p_by_field2d            ! tissues_pool / field2d      | ON_CHUNK & ON_DOMAIN | field2d(:,1) & field2d(:,:)
    PROCEDURE                 :: Div_veg_tiss_p_by_veg_tiss_p         ! tissues_pool / tissues_pool | ON_CHUNK & ON_DOMAIN
    GENERIC                   :: OPERATOR(/) => Div_veg_tiss_p_by_scalar, &
      &                                         Div_veg_tiss_p_by_field1d, &
      &                                         Div_veg_tiss_p_by_field2d, &
      &                                         Div_veg_tiss_p_by_veg_tiss_p
    ! =
    PROCEDURE                 :: Assign_scalar_to_veg_tiss_p          ! tissues_pool = scalar       | ON_CHUNK & ON_DOMAIN
    PROCEDURE                 :: Assign_field1d_to_veg_tiss_p         ! tissues_pool = field1d(:)   | ON_CHUNK
    PROCEDURE                 :: Assign_field2d_to_veg_tiss_p         ! tissues_pool = field2d(:,:) | ON_CHUNK & ON_DOMAIN | field2d(:,1) & field2d(:,:)
    PROCEDURE                 :: Assign_veg_tiss_p_to_veg_tiss_p      ! tissues_pool = tissues_pool | ON_CHUNK & ON_DOMAIN | only numeric values of elements are assigned
    GENERIC                   :: ASSIGNMENT(=) => Assign_scalar_to_veg_tiss_p, &
      &                                           Assign_field1d_to_veg_tiss_p, &
      &                                           Assign_field2d_to_veg_tiss_p, &
      &                                           Assign_veg_tiss_p_to_veg_tiss_p
    !! finalize
    FINAL     :: finalize_veg_tissues_pool
  END TYPE t_veg_tissues_pool

  !-----------------------------------------------------------------------------------------------------
  !> Type definition for t_veg_main_pool
  !! pool-structure containing pools and fluxes
  !! 1st pool in line
  !! added directly into the t_veg_memory as "CLASS(t_jsb_pool), POINTER :: pools"
  !-----------------------------------------------------------------------------------------------------
  TYPE, EXTENDS(t_veg_elements_pool) :: t_veg_main_pool
    !> 1.0 element variables
    !! inherited from t_veg_elements_pool
    !! 
    !> 2.0 tissues_pools
    !! pools
    TYPE(t_veg_tissues_pool),   POINTER :: veg_pool
    !! fluxes
    TYPE(t_veg_tissues_pool),   POINTER :: veg_growth    
    TYPE(t_veg_tissues_pool),   POINTER :: veg_litterfall
    TYPE(t_veg_tissues_pool),   POINTER :: veg_falloc    
    TYPE(t_veg_tissues_pool),   POINTER :: veg_total_flux
    !> 3.0 element_pools
    !! pools
    TYPE(t_veg_elements_pool),  POINTER :: veg_total_biomass
    !! fluxes
    TYPE(t_veg_elements_pool),  POINTER :: veg_exudation    
    TYPE(t_veg_elements_pool),  POINTER :: veg_establishment
    TYPE(t_veg_elements_pool),  POINTER :: veg_reserve_use  
  END TYPE t_veg_main_pool

  !> list with the element_id of all elements used in VEG_ as defined above
  !! is allocated in Init_PROC_memory
  INTEGER, ALLOCATABLE, SAVE :: all_elements(:)

  !-----------------------------------------------------------------------------------------------------
  !> Type definition for veg memory
  !!
  !! 
  !-----------------------------------------------------------------------------------------------------
  TYPE, EXTENDS(t_jsb_memory) :: t_veg_memory

    !> pools & fluxes
    !! 
    ! the veg_main_pool does contain all other pools & fluxes 
    TYPE(t_veg_main_pool), POINTER :: veg_main_pool

    !> 2D & 3D variables
    !! 
      ! basic plant 2D
      TYPE(t_jsb_var_real2d)          :: &   
                                    height            , &            !< plant height [m]
                                    ! diameter          , &            !< plant diameter [m]
                                    ! dens_ind          , &            !< density of individuals [# m-2]
                                    ! lai               , &            !< actual leaf area index [m2 m-2]
                                    ! target_cn_leaf    , &            !< target foliar C:N ratio [mol C mol-1 N]
                                    ! target_np_leaf    , &            !< target foliar N:P ratio [mol N mol-1 P]
                                    ! mean_leaf_age     , &            !< average foliage age in days [days]
                                    ! cohort_age        , &            !< age of the cohort in years, used with veg_dynamics_scheme cohort [years]
                                    ! f_n_demand        , &            !< relative plant demand for N
                                    ! f_p_demand        , &            !< relative plant demand for P
                                    ! beta_sinklim_ps   , &            !< scaling factor to account for a direct sink limitation contraint on photosynthesis [unitless]
                                    ! dphi              , &            !< plant-soil water potential gradient [MPa]
                                    growth_req_n      , &            !< moles N required for a unit of C growth under current allocation fractions [mol N m-2 -1]
                                    growth_req_p      !, &            !< moles P required for a unit of N growth under current allocation fractions [mol P m-2 -1]
                                    ! leaf_cn_direction , &            !< indicates direction of change in leaf CN ratio to maximise NPP [unitless]
                                    ! target_lai        , &            !< LAI implied from sap wood area [m2 m-2]
                                    ! kstar_labile                     !< labile pool turnover rate assuming only temperature and moisture constraints, no stoichiometry [day-1]
      ! 1.1 memory of conditions for the labile pool dynamics
      TYPE(t_jsb_var_real2d)          :: &   
                                    gpp_tlabile_mavg               , &  !< gpp  [µmol m-2 s-1]
                                    maint_respiration_tlabile_mavg , &  !< maintenance respiration [µmol m-2 s-1]
                                    growth_req_n_tlabile_mavg      , &  !< moles N required for a unit of C growth under current allocation fractions time-averaged [mol N m-2 s-1]
                                    growth_req_p_tlabile_mavg           !< moles P required for a unit of N growth under current allocation fractions time-averaged [mol P m-2 s-1]
      TYPE(t_jsb_var_real2d)          :: &   
                                    t_air_tphen_mavg               , &  !< air temperature time-averaged [K]
                                    t_air_week_mavg                , &  !< air temperature time-averaged [K]
                                    t_air_month_mavg               , &  !< air temperature time-averaged [K]
                                    mean_leaf_age                       !< average foliage age in days [days]
      ! various var needed in update_canopy_fluxes()
      TYPE(t_jsb_var_real2d)          :: &
                                    beta_sinklim_ps, &
                                    t_jmax_opt_mavg, &
                                    t_air_tacclim_mavg
      ! fluxes
      TYPE(t_jsb_var_real2d)          :: &   
                                    maint_respiration_pot         !, & !< potential maintenance respiration [µmol CO2 m-2 s-1]
                                    ! maint_respiration             , & !< maintenance respiration [µmol CO2 m-2 s-1]
                                    ! maint_respiration_c13         , & !< maintenance respiration [µmol 13CO2 m-2 s-1]
                                    ! maint_respiration_c14         , & !< maintenance respiration [µmol 14CO2 m-2 s-1]
                                    ! growth_respiration            , & !< growth respiration [µmol CO2 m-2 s-1]
                                    ! growth_respiration_c13        , & !< growth respiration [µmol 13CO2 m-2 s-1]
                                    ! growth_respiration_c14        , & !< growth respiration [µmol 14CO2 m-2 s-1]
                                    ! n_transform_respiration       , & !< respiration associated with N transformation [µmol CO2 m-2 s-1]
                                    ! n_fixation_respiration        , & !< respiration associated with N fixation [µmol CO2 m-2 s-1]
                                    ! n_processing_respiration      , & !< respiration associated with N processing [µmol CO2 m-2 s-1]
                                    ! n_processing_respiration_c13  , & !< respiration associated with N processing [µmol 13CO2 m-2 s-1]
                                    ! n_processing_respiration_c14  , & !< respiration associated with N processing [µmol 14CO2 m-2 s-1]
                                    ! npp                           , & !< net primary production [µmol CO2 m-2 s-1]
                                    ! npp_c13                       , & !< net primary production [µmol 13CO2 m-2 s-1]
                                    ! npp_c14                       , & !< net primary production [µmol 14CO2 m-2 s-1]
                                    ! uptake_nh4                    , & !< NH4 uptake from soil [µmol N m-2 s-1]
                                    ! uptake_nh4_n15                , & !< 15NH4 uptake from soil [µmol N m-2 s-1]
                                    ! uptake_no3                    , & !< NO3 uptake from soil [µmol N m-2 s-1]
                                    ! uptake_no3_n15                , & !< 15NO3 uptake from soil [µmol N m-2 s-1]
                                    ! uptake_po4                    , & !< PO4 uptake from soil [µmol P m-2 s-1]
                                    ! n_fixation                    , & !< symbiontic N fixation [µmol N m-2 s-1]
                                    ! n_fixation_n15                , & !< symbiontic N fixation [µmol 15N m-2 s-1]
                                    ! recycling_leaf_n              , & !< net flux of N from leaf to labile (senescence and maintenance) [mol m-2 timestep-1]
                                    ! recycling_leaf_n15            , & !< net flux of 15N from leaf to labile (senescence and maintenance) [mol m-2 timestep-1]
                                    ! recycling_leaf_p              , & !< net flux of P from leaf to labile (senescence and maintenance) [mol m-2 timestep-1]
                                    ! recycling_fine_root_n         , & !< net flux of N from fine roots to labile (senescence and maintenance) [mol m-2 timestep-1]
                                    ! recycling_fine_root_n15       , & !< net flux of 15N from fine roots to labile (senescence and maintenance) [mol m-2 timestep-1]
                                    ! recycling_fine_root_p         , & !< net flux of P from fine roots to labile (senescence and maintenance) [mol m-2 timestep-1]
                                    ! recycling_heart_wood_n        , & !< net flux of N from sapwood to labile (senescence and maintenance) [mol m-2 timestep-1]
                                    ! recycling_heart_wood_n15      , & !< net flux of 15N from sapwood to labile (senescence and maintenance) [mol m-2 timestep-1]
                                    ! recycling_heart_wood_p        , & !< net flux of P from sapwood to labile (senescence and maintenance) [mol m-2 timestep-1]
                                    ! net_growth                        !< npp minus turnover [µmol CO2 m-2 s-1]
      TYPE(t_jsb_var_real3d)          :: &
                                    fleaf_sunlit_tfrac_mavg_cl
      ! canopy layers
      TYPE(t_jsb_var_real3d)          :: &   ! DIMENSION(ncanopy)
                                    leaf_nitrogen_cl  , &           !< N in canopy layer [mmol m-2 LAI]
                                    fn_chl_cl         , &           !< fraction of N in chlorophyll, pigments only [unitless]
                                    fn_et_cl          , &           !< fraction of N in electron transport [unitless]
                                    fn_rub_cl         , &           !< fraction of N in Rubisco [unitless]
                                    fn_pepc_cl        , &           !< fraction of N in PEP carboylase, C4 plants only [unitless]
                                    fn_oth_cl         , &           !< fraction of N not associated with photosynthesis [unitless]
                                    fleaf_sunlit_cl   , &           !< fraction of layer which is sunlit [unitless]
                                    lai_cl            , &           !< leaf area in the canopy layer [m2 m-2]
                                    cumm_lai_cl                     !< cummulative leaf area above the mid-point of the canopy layer [m2 m-2]
      ! 2D var used for vegetation Pool/Flux replacement
      TYPE(t_jsb_var_real2d)          :: &   
                                    veg_pool_leaf_carbon          , &
                                    veg_pool_leaf_nitrogen        , &
                                    veg_pool_leaf_phosphorus      , &
                                    veg_growth_leaf_carbon        , &
                                    veg_litterfall_leaf_carbon
  CONTAINS
    PROCEDURE :: Init => Init_veg_memory
  END TYPE t_veg_memory

  CHARACTER(len=*), PARAMETER :: modname = 'mo_veg_memory_class'

CONTAINS

  !-----------------------------------------------------------------------------------------------------
  !> helper function for t_veg_elements_pool
  !! PROCEDURE :: Get_element_name_by_id => t_veg_elements_pool_get_element_name_by_id
  !-----------------------------------------------------------------------------------------------------
  FUNCTION t_veg_elements_pool_get_element_name_by_id(this, id) RESULT(name)

    CLASS(t_veg_elements_pool), INTENT(in)  :: this
    INTEGER,                    INTENT(in)  :: id
    CHARACTER(LEN=:),           ALLOCATABLE :: name

    ! name may be consistent with the "name" provided to "CALL this%Add_element()"
    SELECT CASE (id)
    CASE (ELEM_C)
      name='carbon'
    CASE (ELEM_N)
      name='nitrogen'
    CASE (ELEM_P)
      name='phosphorus'
    CASE (ELEM_C13)
      name='carbon13'
    CASE (ELEM_C14)
      name='carbon14'
    CASE (ELEM_N15)
      name='nitrogen15'
    END SELECT

  END FUNCTION t_veg_elements_pool_get_element_name_by_id

  !-----------------------------------------------------------------------------------------------------
  !> initializing an object of the VEG_ elements pool "t_veg_elements_pool"
  !! "pool_object%Init()"
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE t_veg_elements_pool_init(this, name, shortname, l_elements, element_list, element_unit)

    CLASS(t_veg_elements_pool), INTENT(inout) :: this               ! object of TYPE t_veg_elements_pool
    CHARACTER(LEN=*),           INTENT(in)    :: name               ! pool name       (No spaces !)
    CHARACTER(LEN=*),           INTENT(in)    :: shortname          ! pool short name (No spaces !)
    LOGICAL,          OPTIONAL, INTENT(in)    :: l_elements         ! IF pool has elements 
    INTEGER,          OPTIONAL, INTENT(in)    :: element_list(:)    ! list of elements; use if only selected elements should be added
    CHARACTER(LEN=*), OPTIONAL, INTENT(in)    :: element_unit       ! unit of elements

    LOGICAL                     :: l_add_elements
    INTEGER, ALLOCATABLE        :: element_list_loc(:)
    INTEGER                     :: i
    CHARACTER(len=*), PARAMETER :: routine = modname//':t_veg_elements_pool_init'

    !--------------------------------------------
    this%name = name
    this%shortname = shortname

    l_add_elements = .TRUE.
    IF (PRESENT(l_elements)) l_add_elements = l_elements
    IF (.NOT. l_add_elements) RETURN

    IF (.NOT. PRESENT(element_unit)) CALL finish(routine, 'Missing unit for elements')

    ALLOCATE(element_list_loc(SIZE(all_elements)))
    element_list_loc(:) = all_elements(:)
    CALL assign_if_present_allocatable(element_list_loc, element_list)  ! selects element_list if available
  
    this%element_unit                  = element_unit
    this%contains_elem_1d_ptr_to_chunk = .TRUE.
  
    !> carbon
    IF (One_of(ELEM_C, element_list_loc) > 0) THEN
      CALL this%Add_element('carbon', 'C', ELEM_C, dim=2)
      i = this%Find_element_position(ELEM_C)
      IF (i > 0) THEN
        SELECT TYPE (selector => this%element_list(i)%p)
        CLASS IS (t_jsb_var_real2d)
          this%carbon     => selector%ptrm1
          this%domain_c   => selector%ptr2d
        END SELECT
      END IF
    END IF
    !> nitrogen
    IF (One_of(ELEM_N, element_list_loc) > 0) THEN
      CALL this%Add_element('nitrogen', 'N', ELEM_N, dim=2)
      i = this%Find_element_position(ELEM_N)
      IF (i > 0) THEN
        SELECT TYPE (selector => this%element_list(i)%p)
        CLASS IS (t_jsb_var_real2d)
          this%nitrogen   => selector%ptrm1
          this%domain_n   => selector%ptr2d
        END SELECT
      END IF
    END IF
    !> phosphorus
    IF (One_of(ELEM_P, element_list_loc) > 0) THEN
      CALL this%Add_element('phosphorus', 'P', ELEM_P, dim=2)
      i = this%Find_element_position(ELEM_P)
      IF (i > 0) THEN
        SELECT TYPE (selector => this%element_list(i)%p)
        CLASS IS (t_jsb_var_real2d)
          this%phosphorus => selector%ptrm1
          this%domain_p   => selector%ptr2d
        END SELECT
      END IF
    END IF
    !> carbon13
    IF (One_of(ELEM_C13, element_list_loc) > 0) THEN
      CALL this%Add_element('carbon13', 'C13', ELEM_C13, dim=2)
      i = this%Find_element_position(ELEM_C13)
      IF (i > 0) THEN
        SELECT TYPE (selector => this%element_list(i)%p)
        CLASS IS (t_jsb_var_real2d)
          this%carbon13   => selector%ptrm1
          this%domain_c13 => selector%ptr2d
        END SELECT
      END IF
    END IF
    !> carbon14
    IF (One_of(ELEM_C14, element_list_loc) > 0) THEN
      CALL this%Add_element('carbon14', 'C14', ELEM_C14, dim=2)
      i = this%Find_element_position(ELEM_C14)
      IF (i > 0) THEN
        SELECT TYPE (selector => this%element_list(i)%p)
        CLASS IS (t_jsb_var_real2d)
          this%carbon14   => selector%ptrm1
          this%domain_c14 => selector%ptr2d
        END SELECT
      END IF
    END IF
    !> nitrogen15
    IF (One_of(ELEM_N15, element_list_loc) > 0) THEN
      CALL this%Add_element('nitrogen15', 'N15', ELEM_N15, dim=2)
      i = this%Find_element_position(ELEM_N15)
      IF (i > 0) THEN
        SELECT TYPE (selector => this%element_list(i)%p)
        CLASS IS (t_jsb_var_real2d)
          this%nitrogen15 => selector%ptrm1
          this%domain_n15 => selector%ptr2d
        END SELECT
      END IF
    END IF

  END SUBROUTINE t_veg_elements_pool_init

  !-----------------------------------------------------------------------------------------------------
  !> initializing a hlp object of the VEG_ elements pool "t_veg_elements_pool" for operator overload
  !! "pool_object%Init_op_ov_hlp()"
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE t_veg_elements_pool_op_ov_hlp_init(this, name, shortname, owner_model_id, subset_type, &
    &                                           nc, nb, l_elements, element_unit)

    CLASS(t_veg_elements_pool), INTENT(inout), TARGET :: this               ! object of TYPE t_veg_elements_pool
    CHARACTER(LEN=*),           INTENT(in)            :: name               ! pool name       (No spaces !)
    CHARACTER(LEN=*),           INTENT(in)            :: shortname          ! pool short name (No spaces !)
    INTEGER,                    INTENT(in)            :: owner_model_id     ! model ID
    INTEGER,                    INTENT(in)            :: subset_type        ! ON_DOMAIN=1 / ON_CHUNK=2
    INTEGER,                    INTENT(in)            :: nc                 ! number of points in chunk
    INTEGER,                    INTENT(in)            :: nb                 ! number of blocks
    LOGICAL,                    INTENT(in)            :: l_elements         ! IF pool has elements 
    CHARACTER(LEN=*),           INTENT(in)            :: element_unit       ! unit of elements

    CHARACTER(len=*), PARAMETER :: routine = modname//':t_veg_elements_pool_op_ov_hlp_init'

    IF (.NOT. l_elements) THEN
      CALL message(routine, 'Not adding elements to the pool structure. Might be an error!')
      RETURN
    ENDIF
    this%name           = name            ! needed just for debugging
    this%shortname      = shortname       ! needed just for debugging
    this%owner_model_id = owner_model_id  ! needed for 'sub = pool%Get_subset()'
    this%element_unit   = element_unit    ! not needed but setting it anyway

    IF (subset_type == ON_CHUNK) THEN
      ! all_elements(:) defined in PROC_memory_class: contains the elements used by this PROC
      ! NOT checking elements; simply assuming all 'all_elements(:)' are used 
      ! [simply to improve runtime by avoiding: 'IF (One_of(ELEM_C, all_elements(:)) > 0) THEN']
      ALLOCATE(this%op_ov_hlp_1d_c(nc))
      this%carbon     => this%op_ov_hlp_1d_c
      ALLOCATE(this%op_ov_hlp_1d_n(nc))
      this%nitrogen   => this%op_ov_hlp_1d_n
      ALLOCATE(this%op_ov_hlp_1d_p(nc))
      this%phosphorus => this%op_ov_hlp_1d_p
      ALLOCATE(this%op_ov_hlp_1d_c13(nc))
      this%carbon13   => this%op_ov_hlp_1d_c13
      ALLOCATE(this%op_ov_hlp_1d_c14(nc))
      this%carbon14   => this%op_ov_hlp_1d_c14
      ALLOCATE(this%op_ov_hlp_1d_n15(nc))
      this%nitrogen15 => this%op_ov_hlp_1d_n15
    ELSE
      ALLOCATE(this%op_ov_hlp_2d_c(nc,nb))
      this%domain_c   => this%op_ov_hlp_2d_c
      ALLOCATE(this%op_ov_hlp_2d_n(nc,nb))
      this%domain_n   => this%op_ov_hlp_2d_n
      ALLOCATE(this%op_ov_hlp_2d_p(nc,nb))
      this%domain_p   => this%op_ov_hlp_2d_p
      ALLOCATE(this%op_ov_hlp_2d_c13(nc,nb))
      this%domain_c13 => this%op_ov_hlp_2d_c13
      ALLOCATE(this%op_ov_hlp_2d_c14(nc,nb))
      this%domain_c14 => this%op_ov_hlp_2d_c14
      ALLOCATE(this%op_ov_hlp_2d_n15(nc,nb))
      this%domain_n15 => this%op_ov_hlp_2d_n15
    ENDIF

  END SUBROUTINE t_veg_elements_pool_op_ov_hlp_init

  !-----------------------------------------------------------------------------------------------------
  !> helper function for t_veg_tissues_pool
  !! PROCEDURE :: Get_element_name_by_id => t_veg_tissues_pool_get_element_name_by_id
  !-----------------------------------------------------------------------------------------------------
  FUNCTION t_veg_tissues_pool_get_element_name_by_id(this, id) RESULT(name)

    CLASS(t_veg_tissues_pool),  INTENT(in)  :: this
    INTEGER,                    INTENT(in)  :: id
    CHARACTER(LEN=:),           ALLOCATABLE :: name

    ! name may be consistent with the "name" provided to "CALL this%Add_element()"
    SELECT CASE (id)
    CASE (ELEM_C)
      name='carbon'
    CASE (ELEM_N)
      name='nitrogen'
    CASE (ELEM_P)
      name='phosphorus'
    CASE (ELEM_C13)
      name='carbon13'
    CASE (ELEM_C14)
      name='carbon14'
    CASE (ELEM_N15)
      name='nitrogen15'
    END SELECT

  END FUNCTION t_veg_tissues_pool_get_element_name_by_id

  !-----------------------------------------------------------------------------------------------------
  !> initializing the element-variables of an object of the VEG_ tissues pool "t_veg_tissues_pool"
  !! "pool_object%Init()"
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE t_veg_tissues_pool_init(this, name, shortname, l_elements, element_list, element_unit)

    CLASS(t_veg_tissues_pool),  INTENT(inout) :: this               ! object of TYPE t_veg_tissues_pool
    CHARACTER(LEN=*),           INTENT(in)    :: name               ! pool name       (No spaces !)
    CHARACTER(LEN=*),           INTENT(in)    :: shortname          ! pool short name (No spaces !)
    LOGICAL,          OPTIONAL, INTENT(in)    :: l_elements         ! IF pool has elements 
    INTEGER,          OPTIONAL, INTENT(in)    :: element_list(:)    ! list of elements; use if only selected elements should be added
    CHARACTER(LEN=*), OPTIONAL, INTENT(in)    :: element_unit       ! unit of elements

    LOGICAL                     :: l_add_elements
    INTEGER, ALLOCATABLE        :: element_list_loc(:)
    INTEGER                     :: i
    CHARACTER(len=*), PARAMETER :: routine = modname//':t_veg_tissues_pool_init'

    !--------------------------------------------
    this%name = name
    this%shortname = shortname

    l_add_elements = .TRUE.
    IF (PRESENT(l_elements)) l_add_elements = l_elements
    IF (.NOT. l_add_elements) RETURN

    IF (.NOT. PRESENT(element_unit)) CALL finish(routine, 'Missing unit for elements')

    ALLOCATE(element_list_loc(SIZE(all_elements)))
    element_list_loc(:) = all_elements(:)
    CALL assign_if_present_allocatable(element_list_loc, element_list)  ! selects element_list if available
  
    this%element_unit                  = element_unit
    this%contains_elem_1d_ptr_to_chunk = .TRUE.
  
    !> carbon
    IF (One_of(ELEM_C, element_list_loc) > 0) THEN
      CALL this%Add_element('carbon', 'C', ELEM_C, dim=2)
      i = this%Find_element_position(ELEM_C)
      IF (i > 0) THEN
        SELECT TYPE (selector => this%element_list(i)%p)
        CLASS IS (t_jsb_var_real2d)
          this%carbon     => selector%ptrm1
          this%domain_c   => selector%ptr2d
        END SELECT
      END IF
    END IF
    !> nitrogen
    IF (One_of(ELEM_N, element_list_loc) > 0) THEN
      CALL this%Add_element('nitrogen', 'N', ELEM_N, dim=2)
      i = this%Find_element_position(ELEM_N)
      IF (i > 0) THEN
        SELECT TYPE (selector => this%element_list(i)%p)
        CLASS IS (t_jsb_var_real2d)
          this%nitrogen   => selector%ptrm1
          this%domain_n   => selector%ptr2d
        END SELECT
      END IF
    END IF
    !> phosphorus
    IF (One_of(ELEM_P, element_list_loc) > 0) THEN
      CALL this%Add_element('phosphorus', 'P', ELEM_P, dim=2)
      i = this%Find_element_position(ELEM_P)
      IF (i > 0) THEN
        SELECT TYPE (selector => this%element_list(i)%p)
        CLASS IS (t_jsb_var_real2d)
          this%phosphorus => selector%ptrm1
          this%domain_p   => selector%ptr2d
        END SELECT
      END IF
    END IF
    !> carbon13
    IF (One_of(ELEM_C13, element_list_loc) > 0) THEN
      CALL this%Add_element('carbon13', 'C13', ELEM_C13, dim=2)
      i = this%Find_element_position(ELEM_C13)
      IF (i > 0) THEN
        SELECT TYPE (selector => this%element_list(i)%p)
        CLASS IS (t_jsb_var_real2d)
          this%carbon13   => selector%ptrm1
          this%domain_c13 => selector%ptr2d
        END SELECT
      END IF
    END IF
    !> carbon14
    IF (One_of(ELEM_C14, element_list_loc) > 0) THEN
      CALL this%Add_element('carbon14', 'C14', ELEM_C14, dim=2)
      i = this%Find_element_position(ELEM_C14)
      IF (i > 0) THEN
        SELECT TYPE (selector => this%element_list(i)%p)
        CLASS IS (t_jsb_var_real2d)
          this%carbon14   => selector%ptrm1
          this%domain_c14 => selector%ptr2d
        END SELECT
      END IF
    END IF
    !> nitrogen15
    IF (One_of(ELEM_N15, element_list_loc) > 0) THEN
      CALL this%Add_element('nitrogen15', 'N15', ELEM_N15, dim=2)
      i = this%Find_element_position(ELEM_N15)
      IF (i > 0) THEN
        SELECT TYPE (selector => this%element_list(i)%p)
        CLASS IS (t_jsb_var_real2d)
          this%nitrogen15 => selector%ptrm1
          this%domain_n15 => selector%ptr2d
        END SELECT
      END IF
    END IF

  END SUBROUTINE t_veg_tissues_pool_init

  !-----------------------------------------------------------------------------------------------------
  !> create and init a hlp object of the VEG_ tissues pool "t_veg_tissues_pool" for operator overload
  !! "pool_object%Init_op_ov_hlp()"
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE t_veg_tissues_pool_op_ov_hlp_create_init(this, name, shortname, owner_model_id, subset_type, &
    &                                                 nc, nb, l_elements, element_unit, pool_suffix)

    CLASS(t_veg_tissues_pool),  INTENT(inout), TARGET :: this
    CHARACTER(LEN=*),           INTENT(in)            :: name               ! pool name       (No spaces !)
    CHARACTER(LEN=*),           INTENT(in)            :: shortname          ! pool short name (No spaces !)
    INTEGER,                    INTENT(in)            :: owner_model_id     ! model ID
    INTEGER,                    INTENT(in)            :: subset_type        ! ON_DOMAIN=1 / ON_CHUNK=2
    INTEGER,                    INTENT(in)            :: nc                 ! number of points in chunk
    INTEGER,                    INTENT(in)            :: nb                 ! number of blocks
    LOGICAL,                    INTENT(in)            :: l_elements         ! IF pool has elements 
    CHARACTER(LEN=*),           INTENT(in)            :: element_unit       ! unit of elements
    CHARACTER(LEN=*),           INTENT(in)            :: pool_suffix        ! add to the name of each elements_pool created here

    TYPE(t_veg_elements_pool)   :: veg_elements_pool
    CHARACTER(len=*), PARAMETER :: routine = modname//':t_veg_tissues_pool_op_ov_hlp_create_init'

    IF (.NOT. l_elements) THEN
      CALL message(routine, 'Not adding elements to the pool structure. Might be an error!')
      RETURN
    ENDIF
    this%name           = name            ! needed just for debugging
    this%shortname      = shortname       ! needed just for debugging
    this%owner_model_id = owner_model_id  ! needed for 'sub = pool%Get_subset()'
    this%element_unit   = element_unit    ! not needed but setting it anyway

    !> 1.0 tissue pool: init element var
    !! 
    IF (subset_type == ON_CHUNK) THEN
      ! all_elements(:) defined in PROC_memory_class: contains the elements used by this PROC
      ! NOT checking elements; simply assuming all 'all_elements(:)' are used 
      ! [simply to improve runtime by avoiding: 'IF (One_of(ELEM_C, all_elements(:)) > 0) THEN']
      ALLOCATE(this%op_ov_hlp_1d_c(nc))
      this%carbon     => this%op_ov_hlp_1d_c
      ALLOCATE(this%op_ov_hlp_1d_n(nc))
      this%nitrogen   => this%op_ov_hlp_1d_n
      ALLOCATE(this%op_ov_hlp_1d_p(nc))
      this%phosphorus => this%op_ov_hlp_1d_p
      ALLOCATE(this%op_ov_hlp_1d_c13(nc))
      this%carbon13   => this%op_ov_hlp_1d_c13
      ALLOCATE(this%op_ov_hlp_1d_c14(nc))
      this%carbon14   => this%op_ov_hlp_1d_c14
      ALLOCATE(this%op_ov_hlp_1d_n15(nc))
      this%nitrogen15 => this%op_ov_hlp_1d_n15
    ELSE
      ALLOCATE(this%op_ov_hlp_2d_c(nc,nb))
      this%domain_c   => this%op_ov_hlp_2d_c
      ALLOCATE(this%op_ov_hlp_2d_n(nc,nb))
      this%domain_n   => this%op_ov_hlp_2d_n
      ALLOCATE(this%op_ov_hlp_2d_p(nc,nb))
      this%domain_p   => this%op_ov_hlp_2d_p
      ALLOCATE(this%op_ov_hlp_2d_c13(nc,nb))
      this%domain_c13 => this%op_ov_hlp_2d_c13
      ALLOCATE(this%op_ov_hlp_2d_c14(nc,nb))
      this%domain_c14 => this%op_ov_hlp_2d_c14
      ALLOCATE(this%op_ov_hlp_2d_n15(nc,nb))
      this%domain_n15 => this%op_ov_hlp_2d_n15
    ENDIF

    !> 2.0 element pools: allocate & init & add to tissue pool
    !!     
    ALLOCATE(this%leaf)
    CALL this%leaf%Init_op_ov_hlp(        'op_ov_hlp_leaf'//pool_suffix, 'op_ov_hlp_leaf', &
      &                                   owner_model_id, subset_type, nc, nb, .TRUE., element_unit)
    ALLOCATE(this%fine_root)
    CALL this%fine_root%Init_op_ov_hlp(   'op_ov_hlp_froot'//pool_suffix, 'op_ov_hlp_froot', &
      &                                   owner_model_id, subset_type, nc, nb, .TRUE., element_unit)
    ALLOCATE(this%coarse_root)
    CALL this%coarse_root%Init_op_ov_hlp( 'op_ov_hlp_croot'//pool_suffix, 'op_ov_hlp_croot', &
      &                                   owner_model_id, subset_type, nc, nb, .TRUE., element_unit)
    ALLOCATE(this%sap_wood)
    CALL this%sap_wood%Init_op_ov_hlp(    'op_ov_hlp_swood'//pool_suffix, 'op_ov_hlp_swood', &
      &                                   owner_model_id, subset_type, nc, nb, .TRUE., element_unit)
    ALLOCATE(this%heart_wood)
    CALL this%heart_wood%Init_op_ov_hlp(  'op_ov_hlp_hwood'//pool_suffix, 'op_ov_hlp_hwood', &
      &                                   owner_model_id, subset_type, nc, nb, .TRUE., element_unit)
    ALLOCATE(this%labile)
    CALL this%labile%Init_op_ov_hlp(      'op_ov_hlp_labile'//pool_suffix, 'op_ov_hlp_labile', &
      &                                   owner_model_id, subset_type, nc, nb, .TRUE., element_unit)
    ALLOCATE(this%reserve)
    CALL this%reserve%Init_op_ov_hlp(     'op_ov_hlp_reserve'//pool_suffix, 'op_ov_hlp_reserve', &
      &                                   owner_model_id, subset_type, nc, nb, .TRUE., element_unit)
    ALLOCATE(this%fruit)
    CALL this%fruit%Init_op_ov_hlp(       'op_ov_hlp_fruit'//pool_suffix, 'op_ov_hlp_fruit', &
      &                                   owner_model_id, subset_type, nc, nb, .TRUE., element_unit)
    ALLOCATE(this%seed_bed)
    CALL this%seed_bed%Init_op_ov_hlp(    'op_ov_hlp_sbed'//pool_suffix, 'op_ov_hlp_sbed', &
      &                                   owner_model_id, subset_type, nc, nb, .TRUE., element_unit)

  END SUBROUTINE t_veg_tissues_pool_op_ov_hlp_create_init

  !-----------------------------------------------------------------------------------------------------
  !> create an object of the VEG_ tissues pool "t_veg_tissues_pool"
  !! used for VEG_ memory init
  !!   add and init all its objects of the "t_veg_elements_pool" with "pool_object%Init()"
  !!   "object => Create_veg_tissues_pool()"
  !-----------------------------------------------------------------------------------------------------
  FUNCTION Create_veg_tissues_pool(name, shortname, unit, lct_ids, pool_prefix, pool_suffix) RESULT(veg_tissues_pool)

    CHARACTER(LEN=*),           INTENT(in)  :: name
    CHARACTER(LEN=*),           INTENT(in)  :: shortname
    CHARACTER(LEN=*),           INTENT(in)  :: unit
    INTEGER,                    INTENT(in)  :: lct_ids(:)
    CHARACTER(LEN=*), OPTIONAL, INTENT(in)  :: pool_prefix
    CHARACTER(LEN=*), OPTIONAL, INTENT(in)  :: pool_suffix
    TYPE(t_veg_tissues_pool),   POINTER     :: veg_tissues_pool

    TYPE(t_veg_elements_pool),  POINTER     :: veg_elements_pool
    CHARACTER(LEN=:),           ALLOCATABLE :: prefix_loc, suffix_loc

    prefix_loc = ''
    suffix_loc = ''
    CALL assign_if_present_allocatable(prefix_loc, pool_prefix)
    CALL assign_if_present_allocatable(suffix_loc, pool_suffix)
    IF (prefix_loc /= '') prefix_loc = prefix_loc // '_'
    IF (suffix_loc /= '') suffix_loc = '_' // suffix_loc

    ! tissue pool init
    ALLOCATE(veg_tissues_pool)
    CALL veg_tissues_pool%Init(name, shortname, l_elements=.TRUE., element_unit=unit)
    
    ! element pools allocate & init & add to tissue pool
    ! "shortname" is used for model output
    ! "element_list" is optional and can be used to select specific elements (not needed when all elements may be added)
    ! leaf
    ALLOCATE(veg_elements_pool)
    CALL veg_elements_pool%Init(name='leaf_pool'//suffix_loc, shortname='leaf', l_elements=.TRUE., element_unit=unit)
    CALL veg_tissues_pool%Add_pool(veg_elements_pool)
    veg_tissues_pool%leaf => veg_elements_pool
    NULLIFY(veg_elements_pool)
    ! fine_root
    ALLOCATE(veg_elements_pool)
    CALL veg_elements_pool%Init(name='fine_root_pool'//suffix_loc, shortname='fine_root', l_elements=.TRUE., element_unit=unit)
    CALL veg_tissues_pool%Add_pool(veg_elements_pool)
    veg_tissues_pool%fine_root => veg_elements_pool
    NULLIFY(veg_elements_pool)
    ! coarse_root
    ALLOCATE(veg_elements_pool)
    CALL veg_elements_pool%Init(name='coarse_root_pool'//suffix_loc, shortname='coarse_root', l_elements=.TRUE., element_unit=unit)
    CALL veg_tissues_pool%Add_pool(veg_elements_pool)
    veg_tissues_pool%coarse_root => veg_elements_pool
    NULLIFY(veg_elements_pool)
    ! sap_wood
    ALLOCATE(veg_elements_pool)
    CALL veg_elements_pool%Init(name='sap_wood_pool'//suffix_loc, shortname='sap_wood', l_elements=.TRUE., element_unit=unit)
    CALL veg_tissues_pool%Add_pool(veg_elements_pool)
    veg_tissues_pool%sap_wood => veg_elements_pool
    NULLIFY(veg_elements_pool)
    ! heart_wood
    ALLOCATE(veg_elements_pool)
    CALL veg_elements_pool%Init(name='heart_wood_pool'//suffix_loc, shortname='heart_wood', l_elements=.TRUE., element_unit=unit)
    CALL veg_tissues_pool%Add_pool(veg_elements_pool)
    veg_tissues_pool%heart_wood => veg_elements_pool
    NULLIFY(veg_elements_pool)
    ! labile
    ALLOCATE(veg_elements_pool)
    CALL veg_elements_pool%Init(name='labile_pool'//suffix_loc, shortname='labile', l_elements=.TRUE., element_unit=unit)
    CALL veg_tissues_pool%Add_pool(veg_elements_pool)
    veg_tissues_pool%labile => veg_elements_pool
    NULLIFY(veg_elements_pool)
    ! reserve
    ALLOCATE(veg_elements_pool)
    CALL veg_elements_pool%Init(name='reserve_pool'//suffix_loc, shortname='reserve', l_elements=.TRUE., element_unit=unit)
    CALL veg_tissues_pool%Add_pool(veg_elements_pool)
    veg_tissues_pool%reserve => veg_elements_pool
    NULLIFY(veg_elements_pool)
    ! fruit
    ALLOCATE(veg_elements_pool)
    CALL veg_elements_pool%Init(name='fruit_pool'//suffix_loc, shortname='fruit', l_elements=.TRUE., element_unit=unit)
    CALL veg_tissues_pool%Add_pool(veg_elements_pool)
    veg_tissues_pool%fruit => veg_elements_pool
    NULLIFY(veg_elements_pool)
    ! seed_bed
    ALLOCATE(veg_elements_pool)
    CALL veg_elements_pool%Init(name='seed_bed_pool'//suffix_loc, shortname='seed_bed', l_elements=.TRUE., element_unit=unit)
    CALL veg_tissues_pool%Add_pool(veg_elements_pool)
    veg_tissues_pool%seed_bed => veg_elements_pool
    NULLIFY(veg_elements_pool)

  END FUNCTION Create_veg_tissues_pool

  !-----------------------------------------------------------------------------------------------------
  !> initialize memory for the VEG_ process
  !! pools and 
  !! 2d/3d variables
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE Init_veg_memory(mem, prefix, suffix, lct_ids, model_id)

    USE mo_jsb_varlist,       ONLY: BASIC, MEDIUM, FULL
    USE mo_jsb_io,            ONLY: grib_bits, t_cf, t_grib1, t_grib2, tables
    USE mo_jsb_grid_class,    ONLY: t_jsb_grid, t_jsb_vgrid
    USE mo_jsb_grid,          ONLY: Get_grid, Get_vgrid
    USE mo_jsb_model_class,   ONLY: t_jsb_model
    USE mo_jsb_class,         ONLY: Get_model

    CLASS(t_veg_memory),  INTENT(inout), TARGET  :: mem             !> veg memory
    CHARACTER(len=*),     INTENT(in)             :: prefix          !> process name
    CHARACTER(len=*),     INTENT(in)             :: suffix          !> tile name
    INTEGER,              INTENT(in)             :: lct_ids(:)      !> vector of land cover type IDs
    INTEGER,              INTENT(in)             :: model_id        !> model ID model\%id
    ! ---------------------------
    TYPE(t_jsb_model),          POINTER :: model
    TYPE(t_jsb_grid),           POINTER :: hgrid                        ! Horizontal grid
    TYPE(t_jsb_vgrid),          POINTER :: vgrid_surface                ! Vertical grid
    TYPE(t_jsb_vgrid),          POINTER :: vgrid_canopy                 ! Vertical grid
    TYPE(t_veg_elements_pool),  POINTER :: veg_elements_pool            ! pools /fluxes like veg_exudation & veg_total_biomass
    TYPE(t_veg_tissues_pool),   POINTER :: veg_tissues_pool             ! pools /fluxes like veg_pool & veg_growth
    INTEGER                             :: table
    CHARACTER(len=*), PARAMETER         :: routine = TRIM(modname)//':Init_veg_memory'
    ! ---------------------------
    model         => Get_model(model_id)
    table         =  tables(1)
    hgrid         => Get_grid(model%grid_id)
    vgrid_surface => Get_vgrid('surface')
    vgrid_canopy  => Get_vgrid('canopy_layer_q_')

    !> add memory only for LAND & PFT tiles
    !!
    IF ( One_of(LAND_TYPE, lct_ids(:)) > 0 .OR. &
       & One_of(VEG_TYPE,  lct_ids(:)) > 0) THEN
      
      !=============================================================================================================================
      !> pools & fluxes
      !=============================================================================================================================

      !-----------------------------------------------------------------------------------------------------------------------------
      !> populate all_elements(:) with the element_id of all elements
      !! QUINCY is using always all of the below six elements
      IF (.NOT. ALLOCATED(all_elements)) THEN
        ALLOCATE(all_elements(6))
        all_elements(:) = [ELEM_C, ELEM_N, ELEM_P, ELEM_C13, ELEM_C14, ELEM_N15]
      END IF

      !-----------------------------------------------------------------------------------------------------------------------------
      !> create and populate the main pool structure mem%pools that contains all other pools
      !! note: no spaces in 'name' and/or 'short name'
      ALLOCATE(t_veg_main_pool :: mem%pools)
      CALL mem%pools%Init(name='vegetation_main_pool', shortname='veg_main_pool', l_elements=.FALSE.)

      SELECT TYPE (selected => mem%pools)
      CLASS IS (t_veg_main_pool)
        mem%veg_main_pool => selected
      END SELECT

      !-----------------------------------------------------------------------------------------------------------------------------
      !> create the vegetation tissue pools & fluxes (veg_tissues_pool) 
      !! veg_pool / veg_growth / ...
      ! create pools & fluxes
      ! "shortname" is used for model output names
      ! veg_pool
      veg_tissues_pool => Create_veg_tissues_pool(name='vegetation_pool', shortname='veg_pool', unit='mol m-2', &
        & lct_ids=lct_ids, pool_suffix='vegpool') ! l_elements=.TRUE. by default
      CALL mem%veg_main_pool%Add_pool(veg_tissues_pool)
      mem%veg_main_pool%veg_pool => veg_tissues_pool
      ! veg_growth
      veg_tissues_pool => Create_veg_tissues_pool(name='veg_growth_flux', shortname='veg_growth', unit='mol m-2', &
        & lct_ids=lct_ids, pool_suffix='groflux') ! l_elements=.TRUE. by default
      CALL mem%veg_main_pool%Add_pool(veg_tissues_pool)
      mem%veg_main_pool%veg_growth => veg_tissues_pool
      ! veg_litterfall
      veg_tissues_pool => Create_veg_tissues_pool(name='veg_litterfall_flux', shortname='veg_litterfall', unit='mol m-2', &
        & lct_ids=lct_ids, pool_suffix='litflux') ! l_elements=.TRUE. by default
      CALL mem%veg_main_pool%Add_pool(veg_tissues_pool)
      mem%veg_main_pool%veg_litterfall => veg_tissues_pool
      ! veg_falloc
      veg_tissues_pool => Create_veg_tissues_pool(name='veg_allocation_fraction_flux', shortname='veg_falloc', unit='mol m-2', &
        & lct_ids=lct_ids, pool_suffix='falflux') ! l_elements=.TRUE. by default
      CALL mem%veg_main_pool%Add_pool(veg_tissues_pool)
      mem%veg_main_pool%veg_falloc => veg_tissues_pool
      ! veg_total_flux
      veg_tissues_pool => Create_veg_tissues_pool(name='vegetation_total_flux', shortname='veg_total_flux', unit='mol m-2', &
        & lct_ids=lct_ids, pool_suffix='totflux') ! l_elements=.TRUE. by default
      CALL mem%veg_main_pool%Add_pool(veg_tissues_pool)
      mem%veg_main_pool%veg_total_flux => veg_tissues_pool
      NULLIFY(veg_tissues_pool)
      
      ! init and add element pools to the veg_main_pool
      ! (i.e., add them at the same level as the above veg_tissues_pool veg_pool / veg_growth / ...)
      ! "shortname" is used for model output names
      ! veg_total_biomass
      ALLOCATE(veg_elements_pool)
      CALL veg_elements_pool%Init(name='vegetation_total_biomass', shortname='veg_total_biomass', &
        &                         l_elements=.TRUE., element_unit='mol m-2')
      CALL mem%veg_main_pool%Add_pool(veg_elements_pool)
      mem%veg_main_pool%veg_total_biomass => veg_elements_pool
      NULLIFY(veg_elements_pool)
      ! veg_exudation
      ALLOCATE(veg_elements_pool)
      CALL veg_elements_pool%Init(name='veg_exudation', shortname='veg_exudation', &
        &                         l_elements=.TRUE., element_unit='mol m-2')
      CALL mem%veg_main_pool%Add_pool(veg_elements_pool)
      mem%veg_main_pool%veg_exudation => veg_elements_pool
      NULLIFY(veg_elements_pool)
      ! veg_establishment
      ALLOCATE(veg_elements_pool)
      CALL veg_elements_pool%Init(name='veg_establishment', shortname='veg_establishment', &
        &                         l_elements=.TRUE., element_unit='mol m-2')
      CALL mem%veg_main_pool%Add_pool(veg_elements_pool)
      mem%veg_main_pool%veg_establishment => veg_elements_pool
      NULLIFY(veg_elements_pool)
      ! veg_reserve_use
      ALLOCATE(veg_elements_pool)
      CALL veg_elements_pool%Init(name='veg_reserve_use', shortname='veg_reserve_use', &
        &                         l_elements=.TRUE., element_unit='mol m-2')
      CALL mem%veg_main_pool%Add_pool(veg_elements_pool)
      mem%veg_main_pool%veg_reserve_use => veg_elements_pool
      NULLIFY(veg_elements_pool)

      ! set paths of all above initialized & added pools (type-bound procedure of t_jsb_pool)
      CALL mem%pools%Set_paths()

      ! add all variables in the veg_tissues_pool to the VEG_ memory
      CALL mem%Add_var(mem%veg_main_pool%veg_pool, &
        & hgrid, vgrid_surface, vgrid_canopy,                                                          &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                           &
        & prefix, suffix,                                                                              &
        & loutput=.TRUE.,                                                                              &
        & lrestart=.TRUE., initval_r=0.0_wp )
      CALL mem%Add_var(mem%veg_main_pool%veg_growth, &
        & hgrid, vgrid_surface, vgrid_canopy,                                                          &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                           &
        & prefix, suffix,                                                                              &
        & loutput=.TRUE.,                                                                              &
        & lrestart=.TRUE., initval_r=0.0_wp )
      CALL mem%Add_var(mem%veg_main_pool%veg_litterfall, &
        & hgrid, vgrid_surface, vgrid_canopy,                                                          &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                           &
        & prefix, suffix,                                                                              &
        & loutput=.TRUE.,                                                                              &
        & lrestart=.TRUE., initval_r=0.0_wp )
      CALL mem%Add_var(mem%veg_main_pool%veg_falloc, &
        & hgrid, vgrid_surface, vgrid_canopy,                                                          &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                           &
        & prefix, suffix,                                                                              &
        & loutput=.TRUE.,                                                                              &
        & lrestart=.TRUE., initval_r=0.0_wp )
      CALL mem%Add_var(mem%veg_main_pool%veg_total_flux, &
        & hgrid, vgrid_surface, vgrid_canopy,                                                          &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                           &
        & prefix, suffix,                                                                              &
        & loutput=.TRUE.,                                                                              &
        & lrestart=.TRUE., initval_r=0.0_wp )

      ! add all variables in the veg_elements_pool to the VEG_ memory
      CALL mem%Add_var(mem%veg_main_pool%veg_total_biomass, &
        & hgrid, vgrid_surface, vgrid_canopy,                                                          &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                           &
        & prefix, suffix,                                                                              &
        & loutput=.TRUE.,                                                                              &
        & lrestart=.TRUE., initval_r=0.0_wp )
      CALL mem%Add_var(mem%veg_main_pool%veg_exudation, &
        & hgrid, vgrid_surface, vgrid_canopy,                                                          &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                           &
        & prefix, suffix,                                                                              &
        & loutput=.TRUE.,                                                                              &
        & lrestart=.TRUE., initval_r=0.0_wp )
      CALL mem%Add_var(mem%veg_main_pool%veg_establishment, &
        & hgrid, vgrid_surface, vgrid_canopy,                                                          &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                           &
        & prefix, suffix,                                                                              &
        & loutput=.TRUE.,                                                                              &
        & lrestart=.TRUE., initval_r=0.0_wp )
      CALL mem%Add_var(mem%veg_main_pool%veg_reserve_use, &
        & hgrid, vgrid_surface, vgrid_canopy,                                                          &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                           &
        & prefix, suffix,                                                                              &
        & loutput=.TRUE.,                                                                              &
        & lrestart=.TRUE., initval_r=0.0_wp )

      ! output to logfile all pools & their variables at model-init 
      ! (maybe use for debugging; its quite a long list)
      !IF (.NOT. ASSOCIATED(mem%parent)) THEN
      !  CALL mem%veg_main_pool%Print()
      !  ! CALL mem%veg_main_pool%veg_pool%Print()
      !END IF

      !=============================================================================================================================
      !> 2D & 3D variables
      !=============================================================================================================================

      CALL mem%Add_var('height', mem%height, &
            hgrid, vgrid_surface, &
            t_cf('height', 'm', 'plant height'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            loutput=.TRUE., lrestart=.TRUE., &
            initval_r=0.0_wp) ! initvalue
            
      CALL mem%Add_var('growth_req_n', mem%growth_req_n, &
            hgrid, vgrid_surface, &
            t_cf('growth_req_n', 'mol N m-2 -1', 'moles N required for a unit of C growth under current allocation fractions '), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            loutput=.TRUE., lrestart=.TRUE., &
            initval_r=0.0_wp) ! initvalue
  
      CALL mem%Add_var('growth_req_p', mem%growth_req_p, &
            hgrid, vgrid_surface, &
            t_cf('growth_req_p', 'mol P m-2 -1', 'moles P required for a unit of N growth under current allocation fractions '), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            loutput=.TRUE., lrestart=.TRUE., &
            initval_r=0.0_wp) ! initvalue
            
      CALL mem%Add_var('gpp_tlabile_mavg', mem%gpp_tlabile_mavg, &
            hgrid, vgrid_surface, &
            t_cf('gpp_tlabile_mavg', 'µmol m-2 s-1', 'gpp '), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            loutput=.TRUE., lrestart=.TRUE., &
            initval_r=0.0_wp) ! initvalue
  
      CALL mem%Add_var('maint_respiration_tlabile_mavg', mem%maint_respiration_tlabile_mavg, &
            hgrid, vgrid_surface, &
            t_cf('maint_respiration_tlabile_mavg', 'µmol m-2 s-1', 'maintenance respiration'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            loutput=.TRUE., lrestart=.TRUE., &
            initval_r=0.0_wp) ! initvalue
  
      CALL mem%Add_var('growth_req_n_tlabile_mavg', mem%growth_req_n_tlabile_mavg, &
            hgrid, vgrid_surface, &
            t_cf('growth_req_n_tlabile_mavg', '', 'moles N required for a unit of C growth under current allocation fractions'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            loutput=.TRUE., lrestart=.TRUE., &
            initval_r=0.0_wp) ! initvalue
  
      CALL mem%Add_var('growth_req_p_tlabile_mavg', mem%growth_req_p_tlabile_mavg, &
            hgrid, vgrid_surface, &
            t_cf('growth_req_p_tlabile_mavg', '', 'moles P required for a unit of N growth under current allocation fractions'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            loutput=.TRUE., lrestart=.TRUE., &
            initval_r=0.0_wp) ! initvalue
  
      CALL mem%Add_var('t_air_tphen_mavg', mem%t_air_tphen_mavg, &
            hgrid, vgrid_surface, &
            t_cf('t_air_tphen_mavg', '', 'air temp mavg for phenology calc'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            loutput=.TRUE., lrestart=.TRUE., &
            initval_r=0.0_wp) ! initvalue

      CALL mem%Add_var('t_air_week_mavg', mem%t_air_week_mavg, &
            hgrid, vgrid_surface, &
            t_cf('t_air_week_mavg', '', 'avrg air temperature across 7 days'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            loutput=.TRUE., lrestart=.TRUE., &
            initval_r=0.0_wp) ! initvalue
  
      CALL mem%Add_var('t_air_month_mavg', mem%t_air_month_mavg, &
            hgrid, vgrid_surface, &
            t_cf('t_air_month_mavg', '', 'avrg air temperature across 30 days'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            loutput=.TRUE., lrestart=.TRUE., &
            initval_r=0.0_wp) ! initvalue

      CALL mem%Add_var('mean_leaf_age', mem%mean_leaf_age, &
            hgrid, vgrid_surface, &
            t_cf('mean_leaf_age', '', 'average foliage age in days '), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            loutput=.TRUE., lrestart=.TRUE., &
            initval_r=0.0_wp) ! initvalue
            

      CALL mem%Add_var('beta_sinklim_ps', mem%beta_sinklim_ps, &
            hgrid, vgrid_surface, &
            t_cf('beta_sinklim_ps', '', 'scaling factor to account for a direct sink limitation contraint on photosynthesis'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            loutput=.TRUE., lrestart=.TRUE., &
            initval_r=0.0_wp) ! initvalue

      CALL mem%Add_var('t_jmax_opt_mavg', mem%t_jmax_opt_mavg, &
            hgrid, vgrid_surface, &
            t_cf('t_jmax_opt_mavg', '', 'optimim temperature for electron transport of photosynthesis'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            loutput=.TRUE., lrestart=.TRUE., &
            initval_r=0.0_wp) ! initvalue
            
      CALL mem%Add_var('t_air_tacclim_mavg', mem%t_air_tacclim_mavg, &
            hgrid, vgrid_surface, &
            t_cf('t_air_tacclim_mavg', '', 'average air temperature for respiration acclimation'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            loutput=.TRUE., lrestart=.TRUE., &
            initval_r=0.0_wp) ! initvalue

      CALL mem%Add_var('maint_respiration_pot', mem%maint_respiration_pot, &
            hgrid, vgrid_surface, &
            t_cf('maint_respiration_pot', 'µmol CO2 m-2 s-1', 'potential maintenance respiration'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            loutput=.TRUE., lrestart=.TRUE., &
            initval_r=0.0_wp) ! initvalue
            
      CALL mem%Add_var('fleaf_sunlit_tfrac_mavg_cl', mem%fleaf_sunlit_tfrac_mavg_cl, &
            hgrid, vgrid_canopy, &
            t_cf('fleaf_sunlit_tfrac_mavg_cl', '', 'fraction of sunlit leaves'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            loutput=.TRUE., lrestart=.TRUE., &
            initval_r=0.0_wp) ! initvalue
            
      CALL mem%Add_var('leaf_nitrogen_cl', mem%leaf_nitrogen_cl, &
            hgrid, vgrid_canopy, &
            t_cf('leaf_nitrogen_cl', 'mmol m-2 LAI', 'N in canopy layer'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            loutput=.TRUE., lrestart=.TRUE., &
            initval_r=0.0_wp) ! initvalue
  
      CALL mem%Add_var('fn_chl_cl', mem%fn_chl_cl, &
            hgrid, vgrid_canopy, &
            t_cf('fn_chl_cl', '', 'fraction of N in chlorophyll, pigments only'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            loutput=.TRUE., lrestart=.TRUE., &
            initval_r=0.0_wp) ! initvalue
  
      CALL mem%Add_var('fn_et_cl', mem%fn_et_cl, &
            hgrid, vgrid_canopy, &
            t_cf('fn_et_cl', '', 'fraction of N in electron transport '), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            loutput=.TRUE., lrestart=.TRUE., &
            initval_r=0.0_wp) ! initvalue
  
      CALL mem%Add_var('fn_rub_cl', mem%fn_rub_cl, &
            hgrid, vgrid_canopy, &
            t_cf('fn_rub_cl', '', 'fraction of N in Rubisco '), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            loutput=.TRUE., lrestart=.TRUE., &
            initval_r=0.0_wp) ! initvalue
  
      CALL mem%Add_var('fn_pepc_cl', mem%fn_pepc_cl, &
            hgrid, vgrid_canopy, &
            t_cf('fn_pepc_cl', '', 'fraction of N in PEP carboylase, C4 plants only'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            loutput=.TRUE., lrestart=.TRUE., &
            initval_r=0.0_wp) ! initvalue
  
      CALL mem%Add_var('fn_oth_cl', mem%fn_oth_cl, &
            hgrid, vgrid_canopy, &
            t_cf('fn_oth_cl', '', 'fraction of N not associated with photosynthesis '), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            loutput=.TRUE., lrestart=.TRUE., &
            initval_r=0.0_wp) ! initvalue
  
      CALL mem%Add_var('fleaf_sunlit_cl', mem%fleaf_sunlit_cl, &
            hgrid, vgrid_canopy, &
            t_cf('fleaf_sunlit_cl', '', 'fraction of layer which is sunlit'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            loutput=.TRUE., lrestart=.TRUE., &
            initval_r=0.0_wp) ! initvalue
  
      CALL mem%Add_var('lai_cl', mem%lai_cl, &
            hgrid, vgrid_canopy, &
            t_cf('lai_cl', 'm2 m-2', 'leaf area in the canopy layer '), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            loutput=.TRUE., lrestart=.TRUE., &
            initval_r=0.0_wp) ! initvalue
  
      CALL mem%Add_var('cumm_lai_cl', mem%cumm_lai_cl, &
            hgrid, vgrid_canopy, &
            t_cf('cumm_lai_cl', 'm2 m-2', 'cummulative leaf area above the mid-point of the canopy layer '), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            loutput=.TRUE., lrestart=.TRUE., &
            initval_r=0.0_wp) ! initvalue

      CALL mem%Add_var('veg_pool_leaf_carbon', mem%veg_pool_leaf_carbon, &
            hgrid, vgrid_surface, &
            t_cf('veg_pool_leaf_carbon', '[]', 'leaf carbon veg pool'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            loutput=.TRUE., lrestart=.TRUE., &
            initval_r=0.0_wp) ! initvalue

      CALL mem%Add_var('veg_pool_leaf_nitrogen', mem%veg_pool_leaf_nitrogen, &
            hgrid, vgrid_surface, &
            t_cf('veg_pool_leaf_nitrogen', '[]', 'leaf nitrogen veg pool'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            loutput=.TRUE., lrestart=.TRUE., &
            initval_r=0.0_wp) ! initvalue
            
      CALL mem%Add_var('veg_pool_leaf_phosphorus', mem%veg_pool_leaf_phosphorus, &
            hgrid, vgrid_surface, &
            t_cf('veg_pool_leaf_phosphorus', '[]', 'leaf phosphorus veg pool'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            loutput=.TRUE., lrestart=.TRUE., &
            initval_r=0.0_wp) ! initvalue
            
      CALL mem%Add_var('veg_growth_leaf_carbon', mem%veg_growth_leaf_carbon, &
            hgrid, vgrid_surface, &
            t_cf('veg_growth_leaf_carbon', '[]', 'leaf carbon growth flux'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            loutput=.TRUE., lrestart=.TRUE., &
            initval_r=0.0_wp) ! initvalue
            
      CALL mem%Add_var('veg_litterfall_leaf_carbon', mem%veg_litterfall_leaf_carbon, &
            hgrid, vgrid_surface, &
            t_cf('veg_litterfall_leaf_carbon', '[]', 'leaf crabon litterfall flux'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            loutput=.TRUE., lrestart=.TRUE., &
            initval_r=0.0_wp) ! initvalue

    ENDIF  ! IF One_of(LAND_TYPE  OR VEG_TYPE, lct_ids)

  END SUBROUTINE Init_veg_memory

  !--------------------------------------------------------------------------------------------------------------------------------
  !> operator overload and assignment routines for element pools
  !!  TYPE t_veg_elements_pool
  !! 
  !--------------------------------------------------------------------------------------------------------------------------------

  !--------------------------------------------------------------------------------------------------------------------------------
  ! finalize after operator overload
  ! t_veg_elements_pool
  !--------------------------------------------------------------------------------------------------------------------------------
  SUBROUTINE finalize_veg_elements_pool(this)

    TYPE(t_veg_elements_pool) :: this

    CHARACTER(len=*), PARAMETER :: routine = modname//':finalize_veg_elements_pool'

    ! CALL message(routine, 'Finalizing '//TRIM(this%name))
    IF (ASSOCIATED(this%carbon))            NULLIFY(this%carbon)
    IF (ASSOCIATED(this%nitrogen))          NULLIFY(this%nitrogen)
    IF (ASSOCIATED(this%phosphorus))        NULLIFY(this%phosphorus)
    IF (ASSOCIATED(this%carbon13))          NULLIFY(this%carbon13)
    IF (ASSOCIATED(this%carbon14))          NULLIFY(this%carbon14)
    IF (ASSOCIATED(this%nitrogen15))        NULLIFY(this%nitrogen15)
    IF (ASSOCIATED(this%domain_c))          NULLIFY(this%domain_c)
    IF (ASSOCIATED(this%domain_n))          NULLIFY(this%domain_n)
    IF (ASSOCIATED(this%domain_p))          NULLIFY(this%domain_p)
    IF (ASSOCIATED(this%domain_c13))        NULLIFY(this%domain_c13)
    IF (ASSOCIATED(this%domain_c14))        NULLIFY(this%domain_c14)
    IF (ASSOCIATED(this%domain_n15))        NULLIFY(this%domain_n15)

    IF (ALLOCATED(this%op_ov_hlp_1d_c))     DEALLOCATE(this%op_ov_hlp_1d_c)
    IF (ALLOCATED(this%op_ov_hlp_1d_n))     DEALLOCATE(this%op_ov_hlp_1d_n)
    IF (ALLOCATED(this%op_ov_hlp_1d_p))     DEALLOCATE(this%op_ov_hlp_1d_p)
    IF (ALLOCATED(this%op_ov_hlp_1d_c13))   DEALLOCATE(this%op_ov_hlp_1d_c13)
    IF (ALLOCATED(this%op_ov_hlp_1d_c14))   DEALLOCATE(this%op_ov_hlp_1d_c14)
    IF (ALLOCATED(this%op_ov_hlp_1d_n15))   DEALLOCATE(this%op_ov_hlp_1d_n15)
    IF (ALLOCATED(this%op_ov_hlp_2d_c))     DEALLOCATE(this%op_ov_hlp_2d_c)
    IF (ALLOCATED(this%op_ov_hlp_2d_n))     DEALLOCATE(this%op_ov_hlp_2d_n)
    IF (ALLOCATED(this%op_ov_hlp_2d_p))     DEALLOCATE(this%op_ov_hlp_2d_p)
    IF (ALLOCATED(this%op_ov_hlp_2d_c13))   DEALLOCATE(this%op_ov_hlp_2d_c13)
    IF (ALLOCATED(this%op_ov_hlp_2d_c14))   DEALLOCATE(this%op_ov_hlp_2d_c14)
    IF (ALLOCATED(this%op_ov_hlp_2d_n15))   DEALLOCATE(this%op_ov_hlp_2d_n15)
    
    IF (ALLOCATED(this%shortname))          DEALLOCATE(this%shortname)
    IF (ALLOCATED(this%name))               DEALLOCATE(this%name)
    ! not allocated: this%path
    IF (ALLOCATED(this%element_unit))       DEALLOCATE(this%element_unit)
    IF (ASSOCIATED(this%pool_list))         DEALLOCATE(this%pool_list)
    IF (ASSOCIATED(this%leaf_list))         DEALLOCATE(this%leaf_list)
    IF (ALLOCATED(this%element_list))       DEALLOCATE(this%element_list)
    IF (ALLOCATED(this%var_list))           DEALLOCATE(this%var_list)    
    ! this%temporary of t_hermetic is PRIVATE

  END SUBROUTINE finalize_veg_elements_pool
  
  !--------------------------------------------------------------------------------------------------------------------------------
  ! = assign
  !--------------------------------------------------------------------------------------------------------------------------------
  !--------------------------------------------------------------------------------------------------------------------------------
  ! t_veg_elements_pool = scalar | ON_CHUNK & ON_DOMAIN
  !--------------------------------------------------------------------------------------------------------------------------------
  SUBROUTINE Assign_scalar_to_veg_elem_p(lhs, rhs)

    CLASS(t_veg_elements_pool), INTENT(inout) :: lhs
    REAL(wp),                   INTENT(in)    :: rhs

    TYPE(t_subset)              :: sub
    CHARACTER(len=*), PARAMETER :: routine = modname//':Assign_scalar_to_veg_elem_p'

    sub = lhs%Get_subset()
    
    IF (sub%type == ON_CHUNK) THEN
      lhs%carbon     = rhs
      lhs%nitrogen   = rhs
      lhs%phosphorus = rhs
      lhs%carbon13   = rhs
      lhs%carbon14   = rhs
      lhs%nitrogen15 = rhs
    ELSE
      lhs%domain_c   = rhs
      lhs%domain_n   = rhs
      lhs%domain_p   = rhs
      lhs%domain_c13 = rhs
      lhs%domain_c14 = rhs
      lhs%domain_n15 = rhs
    END IF

  END SUBROUTINE Assign_scalar_to_veg_elem_p
  !--------------------------------------------------------------------------------------------------------------------------------
  ! t_veg_elements_pool = field1d(:) | ON_CHUNK
  !--------------------------------------------------------------------------------------------------------------------------------
  SUBROUTINE Assign_field1d_to_veg_elem_p(lhs, rhs)
  
    CLASS(t_veg_elements_pool), INTENT(inout) :: lhs
    REAL(wp),                   INTENT(in)    :: rhs(:)

    TYPE(t_subset)              :: sub
    CHARACTER(len=*), PARAMETER :: routine = modname//':Assign_field1d_to_veg_elem_p'
    
    sub = lhs%Get_subset()
    
    IF (sub%type == ON_CHUNK) THEN
      lhs%carbon     = rhs(:)
      lhs%nitrogen   = rhs(:)
      lhs%phosphorus = rhs(:)
      lhs%carbon13   = rhs(:)
      lhs%carbon14   = rhs(:)
      lhs%nitrogen15 = rhs(:)
    ELSE
      CALL finish(routine, lhs%shortname//' Dimension mismatch, cannot assign 1d-field to jsb_var real2d on domain')
    END IF
  
  END SUBROUTINE Assign_field1d_to_veg_elem_p
  !--------------------------------------------------------------------------------------------------------------------------------
  ! t_veg_elements_pool = field2d(:,:) | ON_CHUNK & ON_DOMAIN | field2d(:,1) & field2d(:,:)
  !--------------------------------------------------------------------------------------------------------------------------------
  SUBROUTINE Assign_field2d_to_veg_elem_p(lhs, rhs)
  
    CLASS(t_veg_elements_pool), INTENT(inout) :: lhs
    REAL(wp),                   INTENT(in)    :: rhs(:,:)

    TYPE(t_subset)              :: sub
    CHARACTER(len=*), PARAMETER :: routine = modname//':Assign_field2d_to_veg_elem_p'

    sub = lhs%Get_subset()

    IF (sub%type == ON_CHUNK) THEN
      IF (SIZE(rhs,2) > 1) CALL finish(routine, 'Dimension mismatch, field2d with nblks > 1; assignment to REAL pool%element(:)')
      lhs%carbon     = rhs(:,1)
      lhs%nitrogen   = rhs(:,1)
      lhs%phosphorus = rhs(:,1)
      lhs%carbon13   = rhs(:,1)
      lhs%carbon14   = rhs(:,1)
      lhs%nitrogen15 = rhs(:,1)
    ELSE
      lhs%domain_c   = rhs(:,:)
      lhs%domain_n   = rhs(:,:)
      lhs%domain_p   = rhs(:,:)
      lhs%domain_c13 = rhs(:,:)
      lhs%domain_c14 = rhs(:,:)
      lhs%domain_n15 = rhs(:,:)
    END IF
  
  END SUBROUTINE Assign_field2d_to_veg_elem_p
  !--------------------------------------------------------------------------------------------------------------------------------
  ! t_veg_elements_pool = t_veg_elements_pool | ON_CHUNK & ON_DOMAIN
  !--------------------------------------------------------------------------------------------------------------------------------
  SUBROUTINE Assign_veg_elem_p_to_veg_elem_p(lhs, rhs)

    CLASS(t_veg_elements_pool), INTENT(inout) :: lhs
    CLASS(t_veg_elements_pool), INTENT(in)    :: rhs

    TYPE(t_subset)              :: sub
    CHARACTER(len=*), PARAMETER :: routine = modname//':Assign_veg_elem_p_to_veg_elem_p'
    
    sub = lhs%Get_subset()

    IF (sub%type == ON_CHUNK) THEN
      lhs%carbon     = rhs%carbon    
      lhs%nitrogen   = rhs%nitrogen  
      lhs%phosphorus = rhs%phosphorus
      lhs%carbon13   = rhs%carbon13  
      lhs%carbon14   = rhs%carbon14  
      lhs%nitrogen15 = rhs%nitrogen15
    ELSE
      lhs%domain_c   = rhs%domain_c  
      lhs%domain_n   = rhs%domain_n  
      lhs%domain_p   = rhs%domain_p  
      lhs%domain_c13 = rhs%domain_c13
      lhs%domain_c14 = rhs%domain_c14
      lhs%domain_n15 = rhs%domain_n15
    END IF

  END SUBROUTINE Assign_veg_elem_p_to_veg_elem_p

  !--------------------------------------------------------------------------------------------------------------------------------
  ! + addition
  !--------------------------------------------------------------------------------------------------------------------------------
  !--------------------------------------------------------------------------------------------------------------------------------
  ! t_veg_elements_pool + scalar        | ON_CHUNK & ON_DOMAIN
  !--------------------------------------------------------------------------------------------------------------------------------
  FUNCTION Add_veg_elem_p_and_scalar(struct1, scalar) RESULT(structOut)
  
    CLASS(t_veg_elements_pool),     INTENT(in)  :: struct1
    REAL(wp),                       INTENT(in)  :: scalar
    TYPE(t_veg_elements_pool)                   :: structOut

    TYPE(t_subset)              :: sub  
    CHARACTER(len=*), PARAMETER :: routine = modname//':Add_veg_elem_p_and_scalar'

    sub = struct1%Get_subset()
    CALL structOut%Init_op_ov_hlp('op_ov_hlp', 'op_ov_hlp', struct1%owner_model_id, &
      &                           sub%type, sub%nc, sub%nb, .TRUE., TRIM(struct1%element_unit))

    ! run calculation
    IF (sub%type == ON_CHUNK) THEN
      structOut%carbon          = struct1%carbon      + scalar
      structOut%nitrogen        = struct1%nitrogen    + scalar
      structOut%phosphorus      = struct1%phosphorus  + scalar
      structOut%carbon13        = struct1%carbon13    + scalar
      structOut%carbon14        = struct1%carbon14    + scalar
      structOut%nitrogen15      = struct1%nitrogen15  + scalar
    ELSE
      structOut%domain_c        = struct1%domain_c    + scalar
      structOut%domain_n        = struct1%domain_n    + scalar
      structOut%domain_p        = struct1%domain_p    + scalar
      structOut%domain_c13      = struct1%domain_c13  + scalar
      structOut%domain_c14      = struct1%domain_c14  + scalar
      structOut%domain_n15      = struct1%domain_n15  + scalar
    END IF

  END FUNCTION Add_veg_elem_p_and_scalar
  !--------------------------------------------------------------------------------------------------------------------------------
  ! scalar + t_veg_elements_pool        | ON_CHUNK & ON_DOMAIN
  !--------------------------------------------------------------------------------------------------------------------------------
  FUNCTION Add_scalar_and_veg_elem_p(scalar, struct1) RESULT(structOut)
  
    REAL(wp),                       INTENT(in)  :: scalar
    CLASS(t_veg_elements_pool),     INTENT(in)  :: struct1
    TYPE(t_veg_elements_pool)                   :: structOut

    TYPE(t_subset)              :: sub  
    CHARACTER(len=*), PARAMETER :: routine = modname//':Add_scalar_and_veg_elem_p'

    sub = struct1%Get_subset()
    CALL structOut%Init_op_ov_hlp('op_ov_hlp', 'op_ov_hlp', struct1%owner_model_id, &
      &                           sub%type, sub%nc, sub%nb, .TRUE., TRIM(struct1%element_unit))

    ! run calculation
    IF (sub%type == ON_CHUNK) THEN
      structOut%carbon          = scalar + struct1%carbon    
      structOut%nitrogen        = scalar + struct1%nitrogen  
      structOut%phosphorus      = scalar + struct1%phosphorus
      structOut%carbon13        = scalar + struct1%carbon13  
      structOut%carbon14        = scalar + struct1%carbon14  
      structOut%nitrogen15      = scalar + struct1%nitrogen15
    ELSE
      structOut%domain_c        = scalar + struct1%domain_c  
      structOut%domain_n        = scalar + struct1%domain_n  
      structOut%domain_p        = scalar + struct1%domain_p  
      structOut%domain_c13      = scalar + struct1%domain_c13
      structOut%domain_c14      = scalar + struct1%domain_c14
      structOut%domain_n15      = scalar + struct1%domain_n15
    END IF

  END FUNCTION Add_scalar_and_veg_elem_p
  !--------------------------------------------------------------------------------------------------------------------------------
  ! t_veg_elements_pool + field1d(:)        | ON_CHUNK
  !--------------------------------------------------------------------------------------------------------------------------------
  FUNCTION Add_veg_elem_p_and_field1d(struct1, field1d) RESULT(structOut)
  
    CLASS(t_veg_elements_pool),     INTENT(in)  :: struct1
    REAL(wp),                       INTENT(in)  :: field1d(:)
    TYPE(t_veg_elements_pool)                   :: structOut

    TYPE(t_subset)              :: sub  
    CHARACTER(len=*), PARAMETER :: routine = modname//':Add_veg_elem_p_and_field1d'

    sub = struct1%Get_subset()
    CALL structOut%Init_op_ov_hlp('op_ov_hlp', 'op_ov_hlp', struct1%owner_model_id, &
      &                           sub%type, sub%nc, sub%nb, .TRUE., TRIM(struct1%element_unit))

    ! run calculation
    IF (sub%type == ON_CHUNK) THEN
      structOut%carbon          = struct1%carbon      + field1d(:)
      structOut%nitrogen        = struct1%nitrogen    + field1d(:)
      structOut%phosphorus      = struct1%phosphorus  + field1d(:)
      structOut%carbon13        = struct1%carbon13    + field1d(:)
      structOut%carbon14        = struct1%carbon14    + field1d(:)
      structOut%nitrogen15      = struct1%nitrogen15  + field1d(:)
    ELSE
      CALL finish(routine, 'Dimension mismatch, cannot calc with 1d-field and jsb_var real2d on domain')
    END IF

  END FUNCTION Add_veg_elem_p_and_field1d
  !--------------------------------------------------------------------------------------------------------------------------------
  ! field1d(:) + t_veg_elements_pool        | ON_CHUNK
  !--------------------------------------------------------------------------------------------------------------------------------
  FUNCTION Add_field1d_and_veg_elem_p(field1d, struct1) RESULT(structOut)
  
    REAL(wp),                       INTENT(in)  :: field1d(:)
    CLASS(t_veg_elements_pool),     INTENT(in)  :: struct1
    TYPE(t_veg_elements_pool)                   :: structOut

    TYPE(t_subset)              :: sub  
    CHARACTER(len=*), PARAMETER :: routine = modname//':Add_field1d_and_veg_elem_p'

    sub = struct1%Get_subset()
    CALL structOut%Init_op_ov_hlp('op_ov_hlp', 'op_ov_hlp', struct1%owner_model_id, &
      &                           sub%type, sub%nc, sub%nb, .TRUE., TRIM(struct1%element_unit))

    ! run calculation
    IF (sub%type == ON_CHUNK) THEN
      structOut%carbon          = field1d(:) + struct1%carbon    
      structOut%nitrogen        = field1d(:) + struct1%nitrogen  
      structOut%phosphorus      = field1d(:) + struct1%phosphorus
      structOut%carbon13        = field1d(:) + struct1%carbon13  
      structOut%carbon14        = field1d(:) + struct1%carbon14  
      structOut%nitrogen15      = field1d(:) + struct1%nitrogen15
    ELSE
      CALL finish(routine, 'Dimension mismatch, cannot calc with 1d-field and jsb_var real2d on domain')
    END IF

  END FUNCTION Add_field1d_and_veg_elem_p
  !--------------------------------------------------------------------------------------------------------------------------------
  ! t_veg_elements_pool + field2d(:,:)        | ON_CHUNK & ON_DOMAIN | field2d(:,1) & field2d(:,:)
  !--------------------------------------------------------------------------------------------------------------------------------
  FUNCTION Add_veg_elem_p_and_field2d(struct1, field2d) RESULT(structOut)
  
    CLASS(t_veg_elements_pool),     INTENT(in)  :: struct1
    REAL(wp),                       INTENT(in)  :: field2d(:,:)
    TYPE(t_veg_elements_pool)                   :: structOut

    TYPE(t_subset)              :: sub  
    CHARACTER(len=*), PARAMETER :: routine = modname//':Add_veg_elem_p_and_field2d'

    sub = struct1%Get_subset()
    CALL structOut%Init_op_ov_hlp('op_ov_hlp', 'op_ov_hlp', struct1%owner_model_id, &
      &                           sub%type, sub%nc, sub%nb, .TRUE., TRIM(struct1%element_unit))

    ! run calculation
    IF (sub%type == ON_CHUNK) THEN
      IF (SIZE(field2d,2) > 1) CALL finish(routine, 'Dimension mismatch, field2d with nblks > 1; calc with REAL pool%element(:)')
      structOut%carbon          = struct1%carbon      + field2d(:,1)
      structOut%nitrogen        = struct1%nitrogen    + field2d(:,1)
      structOut%phosphorus      = struct1%phosphorus  + field2d(:,1)
      structOut%carbon13        = struct1%carbon13    + field2d(:,1)
      structOut%carbon14        = struct1%carbon14    + field2d(:,1)
      structOut%nitrogen15      = struct1%nitrogen15  + field2d(:,1)
    ELSE
      structOut%domain_c        = struct1%domain_c    + field2d(:,:)
      structOut%domain_n        = struct1%domain_n    + field2d(:,:)
      structOut%domain_p        = struct1%domain_p    + field2d(:,:)
      structOut%domain_c13      = struct1%domain_c13  + field2d(:,:)
      structOut%domain_c14      = struct1%domain_c14  + field2d(:,:)
      structOut%domain_n15      = struct1%domain_n15  + field2d(:,:)
    END IF

  END FUNCTION Add_veg_elem_p_and_field2d
  !--------------------------------------------------------------------------------------------------------------------------------
  ! t_veg_elements_pool + t_veg_elements_pool | ON_CHUNK & ON_DOMAIN
  !--------------------------------------------------------------------------------------------------------------------------------
  FUNCTION Add_veg_elem_p_and_veg_elem_p(struct1, struct2) RESULT(structOut)
  
    CLASS(t_veg_elements_pool),     INTENT(in)  :: struct1, struct2
    TYPE(t_veg_elements_pool)                   :: structOut

    TYPE(t_subset)              :: sub  
    CHARACTER(len=*), PARAMETER :: routine = modname//':Add_veg_elem_p_and_veg_elem_p'

    sub = struct1%Get_subset()
    CALL structOut%Init_op_ov_hlp('op_ov_hlp', 'op_ov_hlp', struct1%owner_model_id, &
      &                           sub%type, sub%nc, sub%nb, .TRUE., TRIM(struct1%element_unit))

    ! run calculation
    IF (sub%type == ON_CHUNK) THEN
      structOut%carbon          = struct1%carbon      + struct2%carbon
      structOut%nitrogen        = struct1%nitrogen    + struct2%nitrogen  
      structOut%phosphorus      = struct1%phosphorus  + struct2%phosphorus
      structOut%carbon13        = struct1%carbon13    + struct2%carbon13  
      structOut%carbon14        = struct1%carbon14    + struct2%carbon14  
      structOut%nitrogen15      = struct1%nitrogen15  + struct2%nitrogen15
    ELSE
      structOut%domain_c        = struct1%domain_c    + struct2%domain_c  
      structOut%domain_n        = struct1%domain_n    + struct2%domain_n  
      structOut%domain_p        = struct1%domain_p    + struct2%domain_p  
      structOut%domain_c13      = struct1%domain_c13  + struct2%domain_c13
      structOut%domain_c14      = struct1%domain_c14  + struct2%domain_c14
      structOut%domain_n15      = struct1%domain_n15  + struct2%domain_n15
    END IF

  END FUNCTION Add_veg_elem_p_and_veg_elem_p

  !--------------------------------------------------------------------------------------------------------------------------------
  ! - minus
  !--------------------------------------------------------------------------------------------------------------------------------
  !--------------------------------------------------------------------------------------------------------------------------------
  ! t_veg_elements_pool - scalar        | ON_CHUNK & ON_DOMAIN
  !--------------------------------------------------------------------------------------------------------------------------------
  FUNCTION Subtract_veg_elem_p_minus_scalar(struct1, scalar) RESULT(structOut)
  
    CLASS(t_veg_elements_pool),     INTENT(in)  :: struct1
    REAL(wp),                       INTENT(in)  :: scalar
    TYPE(t_veg_elements_pool)                   :: structOut

    TYPE(t_subset)              :: sub  
    CHARACTER(len=*), PARAMETER :: routine = modname//':Subtract_veg_elem_p_minus_scalar'

    sub = struct1%Get_subset()
    CALL structOut%Init_op_ov_hlp('op_ov_hlp', 'op_ov_hlp', struct1%owner_model_id, &
      &                           sub%type, sub%nc, sub%nb, .TRUE., TRIM(struct1%element_unit))

    ! run calculation
    IF (sub%type == ON_CHUNK) THEN
      structOut%carbon          = struct1%carbon      - scalar
      structOut%nitrogen        = struct1%nitrogen    - scalar
      structOut%phosphorus      = struct1%phosphorus  - scalar
      structOut%carbon13        = struct1%carbon13    - scalar
      structOut%carbon14        = struct1%carbon14    - scalar
      structOut%nitrogen15      = struct1%nitrogen15  - scalar
    ELSE
      structOut%domain_c        = struct1%domain_c    - scalar
      structOut%domain_n        = struct1%domain_n    - scalar
      structOut%domain_p        = struct1%domain_p    - scalar
      structOut%domain_c13      = struct1%domain_c13  - scalar
      structOut%domain_c14      = struct1%domain_c14  - scalar
      structOut%domain_n15      = struct1%domain_n15  - scalar
    END IF

  END FUNCTION Subtract_veg_elem_p_minus_scalar
  !--------------------------------------------------------------------------------------------------------------------------------
  ! t_veg_elements_pool - field1d(:)        | ON_CHUNK
  !--------------------------------------------------------------------------------------------------------------------------------
  FUNCTION Subtract_veg_elem_p_minus_field1d(struct1, field1d) RESULT(structOut)
  
    CLASS(t_veg_elements_pool),     INTENT(in)  :: struct1
    REAL(wp),                       INTENT(in)  :: field1d(:)
    TYPE(t_veg_elements_pool)                   :: structOut

    TYPE(t_subset)              :: sub  
    CHARACTER(len=*), PARAMETER :: routine = modname//':Subtract_veg_elem_p_minus_field1d'

    sub = struct1%Get_subset()
    CALL structOut%Init_op_ov_hlp('op_ov_hlp', 'op_ov_hlp', struct1%owner_model_id, &
      &                           sub%type, sub%nc, sub%nb, .TRUE., TRIM(struct1%element_unit))

    ! run calculation
    IF (sub%type == ON_CHUNK) THEN
      structOut%carbon          = struct1%carbon      - field1d(:)
      structOut%nitrogen        = struct1%nitrogen    - field1d(:)
      structOut%phosphorus      = struct1%phosphorus  - field1d(:)
      structOut%carbon13        = struct1%carbon13    - field1d(:)
      structOut%carbon14        = struct1%carbon14    - field1d(:)
      structOut%nitrogen15      = struct1%nitrogen15  - field1d(:)
    ELSE
      CALL finish(routine, 'Dimension mismatch, cannot calc with 1d-field and jsb_var real2d on domain')
    END IF

  END FUNCTION Subtract_veg_elem_p_minus_field1d
  !--------------------------------------------------------------------------------------------------------------------------------
  ! t_veg_elements_pool - field2d(:,:)        | ON_CHUNK & ON_DOMAIN | field2d(:,1) & field2d(:,:)
  !--------------------------------------------------------------------------------------------------------------------------------
  FUNCTION Subtract_veg_elem_p_minus_field2d(struct1, field2d) RESULT(structOut)
  
    CLASS(t_veg_elements_pool),     INTENT(in)  :: struct1
    REAL(wp),                       INTENT(in)  :: field2d(:,:)
    TYPE(t_veg_elements_pool)                   :: structOut

    TYPE(t_subset)              :: sub  
    CHARACTER(len=*), PARAMETER :: routine = modname//':Subtract_veg_elem_p_minus_field2d'

    sub = struct1%Get_subset()
    CALL structOut%Init_op_ov_hlp('op_ov_hlp', 'op_ov_hlp', struct1%owner_model_id, &
      &                           sub%type, sub%nc, sub%nb, .TRUE., TRIM(struct1%element_unit))

    ! run calculation
    IF (sub%type == ON_CHUNK) THEN
      IF (SIZE(field2d,2) > 1) CALL finish(routine, 'Dimension mismatch, field2d with nblks > 1; calc with REAL pool%element(:)')
      structOut%carbon          = struct1%carbon      - field2d(:,1)
      structOut%nitrogen        = struct1%nitrogen    - field2d(:,1)
      structOut%phosphorus      = struct1%phosphorus  - field2d(:,1)
      structOut%carbon13        = struct1%carbon13    - field2d(:,1)
      structOut%carbon14        = struct1%carbon14    - field2d(:,1)
      structOut%nitrogen15      = struct1%nitrogen15  - field2d(:,1)
    ELSE
      structOut%domain_c        = struct1%domain_c    - field2d(:,:)
      structOut%domain_n        = struct1%domain_n    - field2d(:,:)
      structOut%domain_p        = struct1%domain_p    - field2d(:,:)
      structOut%domain_c13      = struct1%domain_c13  - field2d(:,:)
      structOut%domain_c14      = struct1%domain_c14  - field2d(:,:)
      structOut%domain_n15      = struct1%domain_n15  - field2d(:,:)
    END IF

  END FUNCTION Subtract_veg_elem_p_minus_field2d
  !--------------------------------------------------------------------------------------------------------------------------------
  ! t_veg_elements_pool - t_veg_elements_pool | ON_CHUNK & ON_DOMAIN
  !--------------------------------------------------------------------------------------------------------------------------------
  FUNCTION Subtract_veg_elem_p_minus_veg_elem_p(struct1, struct2) RESULT(structOut)
  
    CLASS(t_veg_elements_pool),     INTENT(in)  :: struct1, struct2
    TYPE(t_veg_elements_pool)                   :: structOut

    TYPE(t_subset)              :: sub  
    CHARACTER(len=*), PARAMETER :: routine = modname//':Subtract_veg_elem_p_minus_veg_elem_p'

    sub = struct1%Get_subset()
    CALL structOut%Init_op_ov_hlp('op_ov_hlp', 'op_ov_hlp', struct1%owner_model_id, &
      &                           sub%type, sub%nc, sub%nb, .TRUE., TRIM(struct1%element_unit))

    ! run calculation
    IF (sub%type == ON_CHUNK) THEN
      structOut%carbon          = struct1%carbon      - struct2%carbon
      structOut%nitrogen        = struct1%nitrogen    - struct2%nitrogen  
      structOut%phosphorus      = struct1%phosphorus  - struct2%phosphorus
      structOut%carbon13        = struct1%carbon13    - struct2%carbon13  
      structOut%carbon14        = struct1%carbon14    - struct2%carbon14  
      structOut%nitrogen15      = struct1%nitrogen15  - struct2%nitrogen15
    ELSE
      structOut%domain_c        = struct1%domain_c    - struct2%domain_c  
      structOut%domain_n        = struct1%domain_n    - struct2%domain_n  
      structOut%domain_p        = struct1%domain_p    - struct2%domain_p  
      structOut%domain_c13      = struct1%domain_c13  - struct2%domain_c13
      structOut%domain_c14      = struct1%domain_c14  - struct2%domain_c14
      structOut%domain_n15      = struct1%domain_n15  - struct2%domain_n15
    END IF

  END FUNCTION Subtract_veg_elem_p_minus_veg_elem_p

  !--------------------------------------------------------------------------------------------------------------------------------
  ! * multiply
  !--------------------------------------------------------------------------------------------------------------------------------
  !--------------------------------------------------------------------------------------------------------------------------------
  ! t_veg_elements_pool * scalar        | ON_CHUNK & ON_DOMAIN
  !--------------------------------------------------------------------------------------------------------------------------------
  FUNCTION Multiply_veg_elem_p_with_scalar(struct1, scalar) RESULT(structOut)
  
    CLASS(t_veg_elements_pool),     INTENT(in)  :: struct1
    REAL(wp),                       INTENT(in)  :: scalar
    TYPE(t_veg_elements_pool)                   :: structOut

    TYPE(t_subset)              :: sub  
    CHARACTER(len=*), PARAMETER :: routine = modname//':Multiply_veg_elem_p_with_scalar'

    sub = struct1%Get_subset()
    CALL structOut%Init_op_ov_hlp('op_ov_hlp', 'op_ov_hlp', struct1%owner_model_id, &
      &                           sub%type, sub%nc, sub%nb, .TRUE., TRIM(struct1%element_unit))

    ! run calculation
    IF (sub%type == ON_CHUNK) THEN
      structOut%carbon          = struct1%carbon      * scalar
      structOut%nitrogen        = struct1%nitrogen    * scalar
      structOut%phosphorus      = struct1%phosphorus  * scalar
      structOut%carbon13        = struct1%carbon13    * scalar
      structOut%carbon14        = struct1%carbon14    * scalar
      structOut%nitrogen15      = struct1%nitrogen15  * scalar
    ELSE
      structOut%domain_c        = struct1%domain_c    * scalar
      structOut%domain_n        = struct1%domain_n    * scalar
      structOut%domain_p        = struct1%domain_p    * scalar
      structOut%domain_c13      = struct1%domain_c13  * scalar
      structOut%domain_c14      = struct1%domain_c14  * scalar
      structOut%domain_n15      = struct1%domain_n15  * scalar
    END IF

  END FUNCTION Multiply_veg_elem_p_with_scalar
  !--------------------------------------------------------------------------------------------------------------------------------
  ! scalar * t_veg_elements_pool        | ON_CHUNK & ON_DOMAIN
  !--------------------------------------------------------------------------------------------------------------------------------
  FUNCTION Multiply_scalar_with_veg_elem_p(scalar, struct1) RESULT(structOut)
  
    REAL(wp),                       INTENT(in)  :: scalar
    CLASS(t_veg_elements_pool),     INTENT(in)  :: struct1
    TYPE(t_veg_elements_pool)                   :: structOut

    TYPE(t_subset)              :: sub  
    CHARACTER(len=*), PARAMETER :: routine = modname//':Multiply_scalar_with_veg_elem_p'

    sub = struct1%Get_subset()
    CALL structOut%Init_op_ov_hlp('op_ov_hlp', 'op_ov_hlp', struct1%owner_model_id, &
      &                           sub%type, sub%nc, sub%nb, .TRUE., TRIM(struct1%element_unit))

    ! run calculation
    IF (sub%type == ON_CHUNK) THEN
      structOut%carbon          = scalar * struct1%carbon    
      structOut%nitrogen        = scalar * struct1%nitrogen  
      structOut%phosphorus      = scalar * struct1%phosphorus
      structOut%carbon13        = scalar * struct1%carbon13  
      structOut%carbon14        = scalar * struct1%carbon14  
      structOut%nitrogen15      = scalar * struct1%nitrogen15
    ELSE
      structOut%domain_c        = scalar * struct1%domain_c  
      structOut%domain_n        = scalar * struct1%domain_n  
      structOut%domain_p        = scalar * struct1%domain_p  
      structOut%domain_c13      = scalar * struct1%domain_c13
      structOut%domain_c14      = scalar * struct1%domain_c14
      structOut%domain_n15      = scalar * struct1%domain_n15
    END IF

  END FUNCTION Multiply_scalar_with_veg_elem_p
  !--------------------------------------------------------------------------------------------------------------------------------
  ! t_veg_elements_pool * field1d(:)        | ON_CHUNK
  !--------------------------------------------------------------------------------------------------------------------------------
  FUNCTION Multiply_veg_elem_p_with_field1d(struct1, field1d) RESULT(structOut)
  
    CLASS(t_veg_elements_pool),     INTENT(in)  :: struct1
    REAL(wp),                       INTENT(in)  :: field1d(:)
    TYPE(t_veg_elements_pool)                   :: structOut

    TYPE(t_subset)              :: sub  
    CHARACTER(len=*), PARAMETER :: routine = modname//':Multiply_veg_elem_p_with_field1d'

    sub = struct1%Get_subset()
    CALL structOut%Init_op_ov_hlp('op_ov_hlp', 'op_ov_hlp', struct1%owner_model_id, &
      &                           sub%type, sub%nc, sub%nb, .TRUE., TRIM(struct1%element_unit))

    ! run calculation
    IF (sub%type == ON_CHUNK) THEN
      structOut%carbon          = struct1%carbon      * field1d(:)
      structOut%nitrogen        = struct1%nitrogen    * field1d(:)
      structOut%phosphorus      = struct1%phosphorus  * field1d(:)
      structOut%carbon13        = struct1%carbon13    * field1d(:)
      structOut%carbon14        = struct1%carbon14    * field1d(:)
      structOut%nitrogen15      = struct1%nitrogen15  * field1d(:)
    ELSE
      CALL finish(routine, 'Dimension mismatch, cannot calc with 1d-field and jsb_var real2d on domain')
    END IF

  END FUNCTION Multiply_veg_elem_p_with_field1d
  !--------------------------------------------------------------------------------------------------------------------------------
  ! field1d(:) * t_veg_elements_pool        | ON_CHUNK
  !--------------------------------------------------------------------------------------------------------------------------------
  FUNCTION Multiply_field1d_with_veg_elem_p(field1d, struct1) RESULT(structOut)
  
    REAL(wp),                       INTENT(in)  :: field1d(:)
    CLASS(t_veg_elements_pool),     INTENT(in)  :: struct1
    TYPE(t_veg_elements_pool)                   :: structOut

    TYPE(t_subset)              :: sub  
    CHARACTER(len=*), PARAMETER :: routine = modname//':Multiply_field1d_with_veg_elem_p'

    sub = struct1%Get_subset()
    CALL structOut%Init_op_ov_hlp('op_ov_hlp', 'op_ov_hlp', struct1%owner_model_id, &
      &                           sub%type, sub%nc, sub%nb, .TRUE., TRIM(struct1%element_unit))

    ! run calculation
    IF (sub%type == ON_CHUNK) THEN
      structOut%carbon          = field1d(:) * struct1%carbon    
      structOut%nitrogen        = field1d(:) * struct1%nitrogen  
      structOut%phosphorus      = field1d(:) * struct1%phosphorus
      structOut%carbon13        = field1d(:) * struct1%carbon13  
      structOut%carbon14        = field1d(:) * struct1%carbon14  
      structOut%nitrogen15      = field1d(:) * struct1%nitrogen15
    ELSE
      CALL finish(routine, 'Dimension mismatch, cannot calc with 1d-field and jsb_var real2d on domain')
    END IF

  END FUNCTION Multiply_field1d_with_veg_elem_p
  !--------------------------------------------------------------------------------------------------------------------------------
  ! t_veg_elements_pool * field2d(:,:)        | ON_CHUNK & ON_DOMAIN | field2d(:,1) & field2d(:,:)
  !--------------------------------------------------------------------------------------------------------------------------------
  FUNCTION Multiply_veg_elem_p_with_field2d(struct1, field2d) RESULT(structOut)
  
    CLASS(t_veg_elements_pool),     INTENT(in)  :: struct1
    REAL(wp),                       INTENT(in)  :: field2d(:,:)
    TYPE(t_veg_elements_pool)                   :: structOut

    TYPE(t_subset)              :: sub  
    CHARACTER(len=*), PARAMETER :: routine = modname//':Multiply_veg_elem_p_with_field2d'

    sub = struct1%Get_subset()
    CALL structOut%Init_op_ov_hlp('op_ov_hlp', 'op_ov_hlp', struct1%owner_model_id, &
      &                           sub%type, sub%nc, sub%nb, .TRUE., TRIM(struct1%element_unit))

    ! run calculation
    IF (sub%type == ON_CHUNK) THEN
      IF (SIZE(field2d,2) > 1) CALL finish(routine, 'Dimension mismatch, field2d with nblks > 1; calc with REAL pool%element(:)')
      structOut%carbon          = struct1%carbon      * field2d(:,1)
      structOut%nitrogen        = struct1%nitrogen    * field2d(:,1)
      structOut%phosphorus      = struct1%phosphorus  * field2d(:,1)
      structOut%carbon13        = struct1%carbon13    * field2d(:,1)
      structOut%carbon14        = struct1%carbon14    * field2d(:,1)
      structOut%nitrogen15      = struct1%nitrogen15  * field2d(:,1)
    ELSE
      structOut%domain_c        = struct1%domain_c    * field2d(:,:)
      structOut%domain_n        = struct1%domain_n    * field2d(:,:)
      structOut%domain_p        = struct1%domain_p    * field2d(:,:)
      structOut%domain_c13      = struct1%domain_c13  * field2d(:,:)
      structOut%domain_c14      = struct1%domain_c14  * field2d(:,:)
      structOut%domain_n15      = struct1%domain_n15  * field2d(:,:)
    END IF

  END FUNCTION Multiply_veg_elem_p_with_field2d
  !--------------------------------------------------------------------------------------------------------------------------------
  ! field2d(:,:) * t_veg_elements_pool        | ON_CHUNK & ON_DOMAIN | field2d(:,1) & field2d(:,:)
  !--------------------------------------------------------------------------------------------------------------------------------
  FUNCTION Multiply_field2d_with_veg_elem_p(field2d, struct1) RESULT(structOut)
  
    REAL(wp),                       INTENT(in)  :: field2d(:,:)
    CLASS(t_veg_elements_pool),     INTENT(in)  :: struct1
    TYPE(t_veg_elements_pool)                   :: structOut

    TYPE(t_subset)              :: sub  
    CHARACTER(len=*), PARAMETER :: routine = modname//':Multiply_field2d_with_veg_elem_p'

    sub = struct1%Get_subset()
    CALL structOut%Init_op_ov_hlp('op_ov_hlp', 'op_ov_hlp', struct1%owner_model_id, &
      &                           sub%type, sub%nc, sub%nb, .TRUE., TRIM(struct1%element_unit))

    ! run calculation
    IF (sub%type == ON_CHUNK) THEN
      IF (SIZE(field2d,2) > 1) CALL finish(routine, 'Dimension mismatch, field2d with nblks > 1; calc with REAL pool%element(:)')
      structOut%carbon          = field2d(:,1) * struct1%carbon    
      structOut%nitrogen        = field2d(:,1) * struct1%nitrogen  
      structOut%phosphorus      = field2d(:,1) * struct1%phosphorus
      structOut%carbon13        = field2d(:,1) * struct1%carbon13  
      structOut%carbon14        = field2d(:,1) * struct1%carbon14  
      structOut%nitrogen15      = field2d(:,1) * struct1%nitrogen15
    ELSE
      structOut%domain_c        = field2d(:,:) * struct1%domain_c  
      structOut%domain_n        = field2d(:,:) * struct1%domain_n  
      structOut%domain_p        = field2d(:,:) * struct1%domain_p  
      structOut%domain_c13      = field2d(:,:) * struct1%domain_c13
      structOut%domain_c14      = field2d(:,:) * struct1%domain_c14
      structOut%domain_n15      = field2d(:,:) * struct1%domain_n15
    END IF

  END FUNCTION Multiply_field2d_with_veg_elem_p
  !--------------------------------------------------------------------------------------------------------------------------------
  ! t_veg_elements_pool * t_veg_elements_pool | ON_CHUNK & ON_DOMAIN
  !--------------------------------------------------------------------------------------------------------------------------------
  FUNCTION Multiply_veg_elem_p_with_veg_elem_p(struct1, struct2) RESULT(structOut)
  
    CLASS(t_veg_elements_pool),     INTENT(in)  :: struct1, struct2
    TYPE(t_veg_elements_pool)                   :: structOut

    TYPE(t_subset)              :: sub  
    CHARACTER(len=*), PARAMETER :: routine = modname//':Multiply_veg_elem_p_with_veg_elem_p'

    sub = struct1%Get_subset()
    CALL structOut%Init_op_ov_hlp('op_ov_hlp', 'op_ov_hlp', struct1%owner_model_id, &
      &                           sub%type, sub%nc, sub%nb, .TRUE., TRIM(struct1%element_unit))

    ! run calculation
    IF (sub%type == ON_CHUNK) THEN
      structOut%carbon          = struct1%carbon      * struct2%carbon
      structOut%nitrogen        = struct1%nitrogen    * struct2%nitrogen  
      structOut%phosphorus      = struct1%phosphorus  * struct2%phosphorus
      structOut%carbon13        = struct1%carbon13    * struct2%carbon13  
      structOut%carbon14        = struct1%carbon14    * struct2%carbon14  
      structOut%nitrogen15      = struct1%nitrogen15  * struct2%nitrogen15
    ELSE
      structOut%domain_c        = struct1%domain_c    * struct2%domain_c  
      structOut%domain_n        = struct1%domain_n    * struct2%domain_n  
      structOut%domain_p        = struct1%domain_p    * struct2%domain_p  
      structOut%domain_c13      = struct1%domain_c13  * struct2%domain_c13
      structOut%domain_c14      = struct1%domain_c14  * struct2%domain_c14
      structOut%domain_n15      = struct1%domain_n15  * struct2%domain_n15
    END IF

  END FUNCTION Multiply_veg_elem_p_with_veg_elem_p

  !--------------------------------------------------------------------------------------------------------------------------------
  ! / divide
  !--------------------------------------------------------------------------------------------------------------------------------
  !--------------------------------------------------------------------------------------------------------------------------------
  ! t_veg_elements_pool / scalar        | ON_CHUNK & ON_DOMAIN
  !--------------------------------------------------------------------------------------------------------------------------------
  FUNCTION Div_veg_elem_p_by_scalar(struct1, scalar) RESULT(structOut)
  
    CLASS(t_veg_elements_pool),     INTENT(in)  :: struct1
    REAL(wp),                       INTENT(in)  :: scalar
    TYPE(t_veg_elements_pool)                   :: structOut

    TYPE(t_subset)              :: sub  
    CHARACTER(len=*), PARAMETER :: routine = modname//':Div_veg_elem_p_by_scalar'

    sub = struct1%Get_subset()
    CALL structOut%Init_op_ov_hlp('op_ov_hlp', 'op_ov_hlp', struct1%owner_model_id, &
      &                           sub%type, sub%nc, sub%nb, .TRUE., TRIM(struct1%element_unit))

    ! run calculation
    IF (sub%type == ON_CHUNK) THEN
      structOut%carbon          = struct1%carbon      / scalar
      structOut%nitrogen        = struct1%nitrogen    / scalar
      structOut%phosphorus      = struct1%phosphorus  / scalar
      structOut%carbon13        = struct1%carbon13    / scalar
      structOut%carbon14        = struct1%carbon14    / scalar
      structOut%nitrogen15      = struct1%nitrogen15  / scalar
    ELSE
      structOut%domain_c        = struct1%domain_c    / scalar
      structOut%domain_n        = struct1%domain_n    / scalar
      structOut%domain_p        = struct1%domain_p    / scalar
      structOut%domain_c13      = struct1%domain_c13  / scalar
      structOut%domain_c14      = struct1%domain_c14  / scalar
      structOut%domain_n15      = struct1%domain_n15  / scalar
    END IF

  END FUNCTION Div_veg_elem_p_by_scalar
  !--------------------------------------------------------------------------------------------------------------------------------
  ! t_veg_elements_pool / field1d(:)        | ON_CHUNK
  !--------------------------------------------------------------------------------------------------------------------------------
  FUNCTION Div_veg_elem_p_by_field1d(struct1, field1d) RESULT(structOut)
  
    CLASS(t_veg_elements_pool),     INTENT(in)  :: struct1
    REAL(wp),                       INTENT(in)  :: field1d(:)
    TYPE(t_veg_elements_pool)                   :: structOut

    TYPE(t_subset)              :: sub  
    CHARACTER(len=*), PARAMETER :: routine = modname//':Div_veg_elem_p_by_field1d'

    sub = struct1%Get_subset()
    CALL structOut%Init_op_ov_hlp('op_ov_hlp', 'op_ov_hlp', struct1%owner_model_id, &
      &                           sub%type, sub%nc, sub%nb, .TRUE., TRIM(struct1%element_unit))

    ! run calculation
    IF (sub%type == ON_CHUNK) THEN
      structOut%carbon          = struct1%carbon      / field1d(:)
      structOut%nitrogen        = struct1%nitrogen    / field1d(:)
      structOut%phosphorus      = struct1%phosphorus  / field1d(:)
      structOut%carbon13        = struct1%carbon13    / field1d(:)
      structOut%carbon14        = struct1%carbon14    / field1d(:)
      structOut%nitrogen15      = struct1%nitrogen15  / field1d(:)
    ELSE
      CALL finish(routine, 'Dimension mismatch, cannot calc with 1d-field and jsb_var real2d on domain')
    END IF

  END FUNCTION Div_veg_elem_p_by_field1d
  !--------------------------------------------------------------------------------------------------------------------------------
  ! t_veg_elements_pool / field2d(:,:)        | ON_CHUNK & ON_DOMAIN | field2d(:,1) & field2d(:,:)
  !--------------------------------------------------------------------------------------------------------------------------------
  FUNCTION Div_veg_elem_p_by_field2d(struct1, field2d) RESULT(structOut)
  
    CLASS(t_veg_elements_pool),     INTENT(in)  :: struct1
    REAL(wp),                       INTENT(in)  :: field2d(:,:)
    TYPE(t_veg_elements_pool)                   :: structOut

    TYPE(t_subset)              :: sub  
    CHARACTER(len=*), PARAMETER :: routine = modname//':Div_veg_elem_p_by_field2d'

    sub = struct1%Get_subset()
    CALL structOut%Init_op_ov_hlp('op_ov_hlp', 'op_ov_hlp', struct1%owner_model_id, &
      &                           sub%type, sub%nc, sub%nb, .TRUE., TRIM(struct1%element_unit))

    ! run calculation
    IF (sub%type == ON_CHUNK) THEN
      IF (SIZE(field2d,2) > 1) CALL finish(routine, 'Dimension mismatch, field2d with nblks > 1; calc with REAL pool%element(:)')
      structOut%carbon          = struct1%carbon      / field2d(:,1)
      structOut%nitrogen        = struct1%nitrogen    / field2d(:,1)
      structOut%phosphorus      = struct1%phosphorus  / field2d(:,1)
      structOut%carbon13        = struct1%carbon13    / field2d(:,1)
      structOut%carbon14        = struct1%carbon14    / field2d(:,1)
      structOut%nitrogen15      = struct1%nitrogen15  / field2d(:,1)
    ELSE
      structOut%domain_c        = struct1%domain_c    / field2d(:,:)
      structOut%domain_n        = struct1%domain_n    / field2d(:,:)
      structOut%domain_p        = struct1%domain_p    / field2d(:,:)
      structOut%domain_c13      = struct1%domain_c13  / field2d(:,:)
      structOut%domain_c14      = struct1%domain_c14  / field2d(:,:)
      structOut%domain_n15      = struct1%domain_n15  / field2d(:,:)
    END IF

  END FUNCTION Div_veg_elem_p_by_field2d
  !--------------------------------------------------------------------------------------------------------------------------------
  ! t_veg_elements_pool / t_veg_elements_pool | ON_CHUNK & ON_DOMAIN
  !--------------------------------------------------------------------------------------------------------------------------------
  FUNCTION Div_veg_elem_p_by_veg_elem_p(struct1, struct2) RESULT(structOut)
  
    CLASS(t_veg_elements_pool),     INTENT(in)  :: struct1, struct2
    TYPE(t_veg_elements_pool)                   :: structOut

    TYPE(t_subset)              :: sub  
    CHARACTER(len=*), PARAMETER :: routine = modname//':Add_veg_elem_p_and_veg_elem_p'

    sub = struct1%Get_subset()
    CALL structOut%Init_op_ov_hlp('op_ov_hlp', 'op_ov_hlp', struct1%owner_model_id, &
      &                           sub%type, sub%nc, sub%nb, .TRUE., TRIM(struct1%element_unit))

    ! run calculation
    IF (sub%type == ON_CHUNK) THEN
      structOut%carbon          = struct1%carbon      / struct2%carbon
      structOut%nitrogen        = struct1%nitrogen    / struct2%nitrogen  
      structOut%phosphorus      = struct1%phosphorus  / struct2%phosphorus
      structOut%carbon13        = struct1%carbon13    / struct2%carbon13  
      structOut%carbon14        = struct1%carbon14    / struct2%carbon14  
      structOut%nitrogen15      = struct1%nitrogen15  / struct2%nitrogen15
    ELSE
      structOut%domain_c        = struct1%domain_c    / struct2%domain_c  
      structOut%domain_n        = struct1%domain_n    / struct2%domain_n  
      structOut%domain_p        = struct1%domain_p    / struct2%domain_p  
      structOut%domain_c13      = struct1%domain_c13  / struct2%domain_c13
      structOut%domain_c14      = struct1%domain_c14  / struct2%domain_c14
      structOut%domain_n15      = struct1%domain_n15  / struct2%domain_n15
    END IF

  END FUNCTION Div_veg_elem_p_by_veg_elem_p

  !--------------------------------------------------------------------------------------------------------------------------------
  !> operator overload and assignment routines for tissue pools
  !!  TYPE t_veg_tissues_pool
  !! 
  !--------------------------------------------------------------------------------------------------------------------------------

  !--------------------------------------------------------------------------------------------------------------------------------
  ! finalize after operator overload
  ! t_veg_tissues_pool
  !--------------------------------------------------------------------------------------------------------------------------------
  SUBROUTINE finalize_veg_tissues_pool(this)

    TYPE(t_veg_tissues_pool) :: this

    CHARACTER(len=*), PARAMETER :: routine = modname//':finalize_veg_tissues_pool'

    ! CALL message(routine, 'Finalizing '//TRIM(this%name))
    IF (ASSOCIATED(this%carbon))            NULLIFY(this%carbon)
    IF (ASSOCIATED(this%nitrogen))          NULLIFY(this%nitrogen)
    IF (ASSOCIATED(this%phosphorus))        NULLIFY(this%phosphorus)
    IF (ASSOCIATED(this%carbon13))          NULLIFY(this%carbon13)
    IF (ASSOCIATED(this%carbon14))          NULLIFY(this%carbon14)
    IF (ASSOCIATED(this%nitrogen15))        NULLIFY(this%nitrogen15)
    IF (ASSOCIATED(this%domain_c))          NULLIFY(this%domain_c)
    IF (ASSOCIATED(this%domain_n))          NULLIFY(this%domain_n)
    IF (ASSOCIATED(this%domain_p))          NULLIFY(this%domain_p)
    IF (ASSOCIATED(this%domain_c13))        NULLIFY(this%domain_c13)
    IF (ASSOCIATED(this%domain_c14))        NULLIFY(this%domain_c14)
    IF (ASSOCIATED(this%domain_n15))        NULLIFY(this%domain_n15)

    IF (ALLOCATED(this%op_ov_hlp_1d_c))     DEALLOCATE(this%op_ov_hlp_1d_c)
    IF (ALLOCATED(this%op_ov_hlp_1d_n))     DEALLOCATE(this%op_ov_hlp_1d_n)
    IF (ALLOCATED(this%op_ov_hlp_1d_p))     DEALLOCATE(this%op_ov_hlp_1d_p)
    IF (ALLOCATED(this%op_ov_hlp_1d_c13))   DEALLOCATE(this%op_ov_hlp_1d_c13)
    IF (ALLOCATED(this%op_ov_hlp_1d_c14))   DEALLOCATE(this%op_ov_hlp_1d_c14)
    IF (ALLOCATED(this%op_ov_hlp_1d_n15))   DEALLOCATE(this%op_ov_hlp_1d_n15)
    IF (ALLOCATED(this%op_ov_hlp_2d_c))     DEALLOCATE(this%op_ov_hlp_2d_c)
    IF (ALLOCATED(this%op_ov_hlp_2d_n))     DEALLOCATE(this%op_ov_hlp_2d_n)
    IF (ALLOCATED(this%op_ov_hlp_2d_p))     DEALLOCATE(this%op_ov_hlp_2d_p)
    IF (ALLOCATED(this%op_ov_hlp_2d_c13))   DEALLOCATE(this%op_ov_hlp_2d_c13)
    IF (ALLOCATED(this%op_ov_hlp_2d_c14))   DEALLOCATE(this%op_ov_hlp_2d_c14)
    IF (ALLOCATED(this%op_ov_hlp_2d_n15))   DEALLOCATE(this%op_ov_hlp_2d_n15)
    
    IF (ASSOCIATED(this%leaf))              DEALLOCATE(this%leaf)
    IF (ASSOCIATED(this%fine_root))         DEALLOCATE(this%fine_root)
    IF (ASSOCIATED(this%coarse_root))       DEALLOCATE(this%coarse_root)
    IF (ASSOCIATED(this%sap_wood))          DEALLOCATE(this%sap_wood)
    IF (ASSOCIATED(this%heart_wood))        DEALLOCATE(this%heart_wood)
    IF (ASSOCIATED(this%labile))            DEALLOCATE(this%labile)
    IF (ASSOCIATED(this%reserve))           DEALLOCATE(this%reserve)
    IF (ASSOCIATED(this%fruit))             DEALLOCATE(this%fruit)
    IF (ASSOCIATED(this%seed_bed))          DEALLOCATE(this%seed_bed)
    
    IF (ALLOCATED(this%shortname))          DEALLOCATE(this%shortname)
    IF (ALLOCATED(this%name))               DEALLOCATE(this%name)
    ! not allocated: this%path
    IF (ALLOCATED(this%element_unit))       DEALLOCATE(this%element_unit)
    IF (ASSOCIATED(this%pool_list))         DEALLOCATE(this%pool_list)
    IF (ASSOCIATED(this%leaf_list))         DEALLOCATE(this%leaf_list)
    IF (ALLOCATED(this%element_list))       DEALLOCATE(this%element_list)
    IF (ALLOCATED(this%var_list))           DEALLOCATE(this%var_list)    
    ! this%temporary of t_hermetic is PRIVATE
   
  END SUBROUTINE finalize_veg_tissues_pool

  !--------------------------------------------------------------------------------------------------------------------------------
  ! = assign
  !--------------------------------------------------------------------------------------------------------------------------------
  !--------------------------------------------------------------------------------------------------------------------------------
  ! veg_tissues_pool = scalar | ON_CHUNK & ON_DOMAIN
  !--------------------------------------------------------------------------------------------------------------------------------
  SUBROUTINE Assign_scalar_to_veg_tiss_p(lhs, rhs)

    CLASS(t_veg_tissues_pool),  INTENT(inout) :: lhs
    REAL(wp),                   INTENT(in)    :: rhs

    TYPE(t_subset)              :: sub
    CHARACTER(len=*), PARAMETER :: routine = modname//':Assign_scalar_to_veg_tiss_p'
    
    sub = lhs%Get_subset()

    IF (sub%type == ON_CHUNK) THEN
      lhs%carbon     = rhs
      lhs%nitrogen   = rhs
      lhs%phosphorus = rhs
      lhs%carbon13   = rhs
      lhs%carbon14   = rhs
      lhs%nitrogen15 = rhs
    ELSE
      lhs%domain_c   = rhs
      lhs%domain_n   = rhs
      lhs%domain_p   = rhs
      lhs%domain_c13 = rhs
      lhs%domain_c14 = rhs
      lhs%domain_n15 = rhs
    END IF
    
    lhs%leaf        = rhs
    lhs%fine_root   = rhs
    lhs%coarse_root = rhs
    lhs%sap_wood    = rhs
    lhs%heart_wood  = rhs
    lhs%labile      = rhs
    lhs%reserve     = rhs
    lhs%fruit       = rhs
    lhs%seed_bed    = rhs
    
  END SUBROUTINE Assign_scalar_to_veg_tiss_p
  !--------------------------------------------------------------------------------------------------------------------------------
  ! veg_tissues_pool = field1d(:) | ON_CHUNK
  !--------------------------------------------------------------------------------------------------------------------------------
  SUBROUTINE Assign_field1d_to_veg_tiss_p(lhs, rhs)
  
    CLASS(t_veg_tissues_pool),  INTENT(inout) :: lhs
    REAL(wp),                   INTENT(in)    :: rhs(:)

    TYPE(t_subset)              :: sub
    CHARACTER(len=*), PARAMETER :: routine = modname//':Assign_field1d_to_veg_tiss_p'
    
    sub = lhs%Get_subset()

    IF (sub%type == ON_CHUNK) THEN
      lhs%carbon     = rhs(:)
      lhs%nitrogen   = rhs(:)
      lhs%phosphorus = rhs(:)
      lhs%carbon13   = rhs(:)
      lhs%carbon14   = rhs(:)
      lhs%nitrogen15 = rhs(:)
    ELSE
      CALL finish(routine, lhs%shortname//' Dimension mismatch, cannot assign 1d-field to jsb_var real2d on domain')
    END IF

    lhs%leaf        = rhs(:)
    lhs%fine_root   = rhs(:)
    lhs%coarse_root = rhs(:)
    lhs%sap_wood    = rhs(:)
    lhs%heart_wood  = rhs(:)
    lhs%labile      = rhs(:)
    lhs%reserve     = rhs(:)
    lhs%fruit       = rhs(:)
    lhs%seed_bed    = rhs(:)
    
  END SUBROUTINE Assign_field1d_to_veg_tiss_p
  !--------------------------------------------------------------------------------------------------------------------------------
  ! veg_tissues_pool = field2d(:,:) | ON_CHUNK & ON_DOMAIN | field2d(:,1) & field2d(:,:)
  !--------------------------------------------------------------------------------------------------------------------------------
  SUBROUTINE Assign_field2d_to_veg_tiss_p(lhs, rhs)
  
    CLASS(t_veg_tissues_pool),  INTENT(inout) :: lhs
    REAL(wp),                   INTENT(in)    :: rhs(:,:)

    TYPE(t_subset)              :: sub
    CHARACTER(len=*), PARAMETER :: routine = modname//':Assign_field2d_to_veg_tiss_p'

    sub = lhs%Get_subset()

    IF (sub%type == ON_CHUNK) THEN
      IF (SIZE(rhs,2) > 1) CALL finish(routine, 'Dimension mismatch, field2d with nblks > 1; assignment to REAL pool%element(:)')
      lhs%carbon     = rhs(:,1)
      lhs%nitrogen   = rhs(:,1)
      lhs%phosphorus = rhs(:,1)
      lhs%carbon13   = rhs(:,1)
      lhs%carbon14   = rhs(:,1)
      lhs%nitrogen15 = rhs(:,1)
    ELSE
      lhs%domain_c   = rhs(:,:)
      lhs%domain_n   = rhs(:,:)
      lhs%domain_p   = rhs(:,:)
      lhs%domain_c13 = rhs(:,:)
      lhs%domain_c14 = rhs(:,:)
      lhs%domain_n15 = rhs(:,:)
    END IF

    lhs%leaf        = rhs(:,:)
    lhs%fine_root   = rhs(:,:)
    lhs%coarse_root = rhs(:,:)
    lhs%sap_wood    = rhs(:,:)
    lhs%heart_wood  = rhs(:,:)
    lhs%labile      = rhs(:,:)
    lhs%reserve     = rhs(:,:)
    lhs%fruit       = rhs(:,:)
    lhs%seed_bed    = rhs(:,:)

  END SUBROUTINE Assign_field2d_to_veg_tiss_p
  !--------------------------------------------------------------------------------------------------------------------------------
  ! veg_tissues_pool = veg_tissues_pool | ON_CHUNK & ON_DOMAIN
  !--------------------------------------------------------------------------------------------------------------------------------
  SUBROUTINE Assign_veg_tiss_p_to_veg_tiss_p(lhs, rhs)

    CLASS(t_veg_tissues_pool), INTENT(inout) :: lhs
    CLASS(t_veg_tissues_pool), INTENT(in)    :: rhs

    TYPE(t_subset)              :: sub
    CHARACTER(len=*), PARAMETER :: routine = modname//':Assign_veg_tiss_p_to_veg_tiss_p'
    
    sub = lhs%Get_subset()

    IF (sub%type == ON_CHUNK) THEN
      lhs%carbon     = rhs%carbon    
      lhs%nitrogen   = rhs%nitrogen  
      lhs%phosphorus = rhs%phosphorus
      lhs%carbon13   = rhs%carbon13  
      lhs%carbon14   = rhs%carbon14  
      lhs%nitrogen15 = rhs%nitrogen15
    ELSE
      lhs%domain_c   = rhs%domain_c  
      lhs%domain_n   = rhs%domain_n  
      lhs%domain_p   = rhs%domain_p  
      lhs%domain_c13 = rhs%domain_c13
      lhs%domain_c14 = rhs%domain_c14
      lhs%domain_n15 = rhs%domain_n15
    END IF
    
    lhs%leaf        = rhs%leaf       
    lhs%fine_root   = rhs%fine_root  
    lhs%coarse_root = rhs%coarse_root
    lhs%sap_wood    = rhs%sap_wood   
    lhs%heart_wood  = rhs%heart_wood 
    lhs%labile      = rhs%labile     
    lhs%reserve     = rhs%reserve    
    lhs%fruit       = rhs%fruit      
    lhs%seed_bed    = rhs%seed_bed   

  END SUBROUTINE Assign_veg_tiss_p_to_veg_tiss_p

  !--------------------------------------------------------------------------------------------------------------------------------
  ! + addition
  !--------------------------------------------------------------------------------------------------------------------------------
  !--------------------------------------------------------------------------------------------------------------------------------
  ! veg_tissues_pool + scalar        | ON_CHUNK & ON_DOMAIN
  !--------------------------------------------------------------------------------------------------------------------------------
  FUNCTION Add_veg_tiss_p_and_scalar(struct1, scalar) RESULT(structOut)

    CLASS(t_veg_tissues_pool),      INTENT(in)  :: struct1
    REAL(wp),                       INTENT(in)  :: scalar
    TYPE(t_veg_tissues_pool)                    :: structOut

    TYPE(t_subset)              :: sub  
    CHARACTER(len=*), PARAMETER :: routine = modname//':Add_veg_tiss_p_and_scalar'

    sub = struct1%Get_subset()
    CALL structOut%Init_op_ov_hlp('op_ov_hlp', 'op_ov_hlp', struct1%owner_model_id, &
      &                           sub%type, sub%nc, sub%nb, .TRUE., TRIM(struct1%element_unit), 'veg_p')

    ! run calculation
    IF (sub%type == ON_CHUNK) THEN
      structOut%carbon          = struct1%carbon      + scalar
      structOut%nitrogen        = struct1%nitrogen    + scalar
      structOut%phosphorus      = struct1%phosphorus  + scalar
      structOut%carbon13        = struct1%carbon13    + scalar
      structOut%carbon14        = struct1%carbon14    + scalar
      structOut%nitrogen15      = struct1%nitrogen15  + scalar
    ELSE
      structOut%domain_c        = struct1%domain_c    + scalar
      structOut%domain_n        = struct1%domain_n    + scalar
      structOut%domain_p        = struct1%domain_p    + scalar
      structOut%domain_c13      = struct1%domain_c13  + scalar
      structOut%domain_c14      = struct1%domain_c14  + scalar
      structOut%domain_n15      = struct1%domain_n15  + scalar
    END IF

    structOut%leaf            = struct1%leaf        + scalar
    structOut%fine_root       = struct1%fine_root   + scalar
    structOut%coarse_root     = struct1%coarse_root + scalar
    structOut%sap_wood        = struct1%sap_wood    + scalar
    structOut%heart_wood      = struct1%heart_wood  + scalar
    structOut%labile          = struct1%labile      + scalar
    structOut%reserve         = struct1%reserve     + scalar
    structOut%fruit           = struct1%fruit       + scalar
    structOut%seed_bed        = struct1%seed_bed    + scalar

  END FUNCTION Add_veg_tiss_p_and_scalar
  !--------------------------------------------------------------------------------------------------------------------------------
  ! scalar + veg_tissues_pool        | ON_CHUNK & ON_DOMAIN
  !--------------------------------------------------------------------------------------------------------------------------------
  FUNCTION Add_scalar_and_veg_tiss_p(scalar, struct1) RESULT(structOut)
  
    REAL(wp),                       INTENT(in)  :: scalar
    CLASS(t_veg_tissues_pool),      INTENT(in)  :: struct1
    TYPE(t_veg_tissues_pool)                    :: structOut
  
    TYPE(t_subset)              :: sub  
    CHARACTER(len=*), PARAMETER :: routine = modname//':Add_scalar_and_veg_tiss_p'
  
    sub = struct1%Get_subset()
    CALL structOut%Init_op_ov_hlp('op_ov_hlp', 'op_ov_hlp', struct1%owner_model_id, &
      &                           sub%type, sub%nc, sub%nb, .TRUE., TRIM(struct1%element_unit), 'veg_p')

    ! run calculation
    IF (sub%type == ON_CHUNK) THEN
      structOut%carbon          = scalar + struct1%carbon    
      structOut%nitrogen        = scalar + struct1%nitrogen  
      structOut%phosphorus      = scalar + struct1%phosphorus
      structOut%carbon13        = scalar + struct1%carbon13  
      structOut%carbon14        = scalar + struct1%carbon14  
      structOut%nitrogen15      = scalar + struct1%nitrogen15
    ELSE
      structOut%domain_c        = scalar + struct1%domain_c  
      structOut%domain_n        = scalar + struct1%domain_n  
      structOut%domain_p        = scalar + struct1%domain_p  
      structOut%domain_c13      = scalar + struct1%domain_c13
      structOut%domain_c14      = scalar + struct1%domain_c14
      structOut%domain_n15      = scalar + struct1%domain_n15
    END IF

    structOut%leaf            = scalar + struct1%leaf       
    structOut%fine_root       = scalar + struct1%fine_root  
    structOut%coarse_root     = scalar + struct1%coarse_root
    structOut%sap_wood        = scalar + struct1%sap_wood   
    structOut%heart_wood      = scalar + struct1%heart_wood 
    structOut%labile          = scalar + struct1%labile     
    structOut%reserve         = scalar + struct1%reserve    
    structOut%fruit           = scalar + struct1%fruit      
    structOut%seed_bed        = scalar + struct1%seed_bed   
  
  END FUNCTION Add_scalar_and_veg_tiss_p
  !--------------------------------------------------------------------------------------------------------------------------------
  ! veg_tissues_pool + field1d(:)        | ON_CHUNK
  !--------------------------------------------------------------------------------------------------------------------------------
  FUNCTION Add_veg_tiss_p_and_field1d(struct1, field1d) RESULT(structOut)
  
    CLASS(t_veg_tissues_pool),      INTENT(in)  :: struct1
    REAL(wp),                       INTENT(in)  :: field1d(:)
    TYPE(t_veg_tissues_pool)                    :: structOut
  
    TYPE(t_subset)              :: sub  
    CHARACTER(len=*), PARAMETER :: routine = modname//':Add_veg_tiss_p_and_field1d'
  
    sub = struct1%Get_subset()
    CALL structOut%Init_op_ov_hlp('op_ov_hlp', 'op_ov_hlp', struct1%owner_model_id, &
      &                           sub%type, sub%nc, sub%nb, .TRUE., TRIM(struct1%element_unit), 'veg_p')
  
    ! run calculation
    IF (sub%type == ON_CHUNK) THEN
      structOut%carbon          = struct1%carbon      + field1d(:)
      structOut%nitrogen        = struct1%nitrogen    + field1d(:)
      structOut%phosphorus      = struct1%phosphorus  + field1d(:)
      structOut%carbon13        = struct1%carbon13    + field1d(:)
      structOut%carbon14        = struct1%carbon14    + field1d(:)
      structOut%nitrogen15      = struct1%nitrogen15  + field1d(:)
    ELSE
      CALL finish(routine, 'Dimension mismatch, cannot calc with 1d-field and jsb_var real2d on domain')
    END IF

    structOut%leaf            = struct1%leaf        + field1d(:)
    structOut%fine_root       = struct1%fine_root   + field1d(:)
    structOut%coarse_root     = struct1%coarse_root + field1d(:)
    structOut%sap_wood        = struct1%sap_wood    + field1d(:)
    structOut%heart_wood      = struct1%heart_wood  + field1d(:)
    structOut%labile          = struct1%labile      + field1d(:)
    structOut%reserve         = struct1%reserve     + field1d(:)
    structOut%fruit           = struct1%fruit       + field1d(:)
    structOut%seed_bed        = struct1%seed_bed    + field1d(:)
  
  END FUNCTION Add_veg_tiss_p_and_field1d
  !--------------------------------------------------------------------------------------------------------------------------------
  ! field1d(:) + veg_tissues_pool        | ON_CHUNK
  !--------------------------------------------------------------------------------------------------------------------------------
  FUNCTION Add_field1d_and_veg_tiss_p(field1d, struct1) RESULT(structOut)
  
    REAL(wp),                       INTENT(in)  :: field1d(:)
    CLASS(t_veg_tissues_pool),      INTENT(in)  :: struct1
    TYPE(t_veg_tissues_pool)                    :: structOut
  
    TYPE(t_subset)              :: sub  
    CHARACTER(len=*), PARAMETER :: routine = modname//':Add_field1d_and_veg_tiss_p'
  
    sub = struct1%Get_subset()
    CALL structOut%Init_op_ov_hlp('op_ov_hlp', 'op_ov_hlp', struct1%owner_model_id, &
      &                           sub%type, sub%nc, sub%nb, .TRUE., TRIM(struct1%element_unit), 'veg_p')
  
    ! run calculation
    IF (sub%type == ON_CHUNK) THEN
      structOut%carbon          = field1d(:) + struct1%carbon    
      structOut%nitrogen        = field1d(:) + struct1%nitrogen  
      structOut%phosphorus      = field1d(:) + struct1%phosphorus
      structOut%carbon13        = field1d(:) + struct1%carbon13  
      structOut%carbon14        = field1d(:) + struct1%carbon14  
      structOut%nitrogen15      = field1d(:) + struct1%nitrogen15
    ELSE
      CALL finish(routine, 'Dimension mismatch, cannot calc with 1d-field and jsb_var real2d on domain')
    END IF

    structOut%leaf            = field1d(:) + struct1%leaf       
    structOut%fine_root       = field1d(:) + struct1%fine_root  
    structOut%coarse_root     = field1d(:) + struct1%coarse_root
    structOut%sap_wood        = field1d(:) + struct1%sap_wood   
    structOut%heart_wood      = field1d(:) + struct1%heart_wood 
    structOut%labile          = field1d(:) + struct1%labile     
    structOut%reserve         = field1d(:) + struct1%reserve    
    structOut%fruit           = field1d(:) + struct1%fruit      
    structOut%seed_bed        = field1d(:) + struct1%seed_bed   
  
  END FUNCTION Add_field1d_and_veg_tiss_p
  !--------------------------------------------------------------------------------------------------------------------------------
  ! veg_tissues_pool + field2d(:,:)        | ON_CHUNK & ON_DOMAIN | field2d(:,1) & field2d(:,:)
  !--------------------------------------------------------------------------------------------------------------------------------
  FUNCTION Add_veg_tiss_p_and_field2d(struct1, field2d) RESULT(structOut)
  
    CLASS(t_veg_tissues_pool),      INTENT(in)  :: struct1
    REAL(wp),                       INTENT(in)  :: field2d(:,:)
    TYPE(t_veg_tissues_pool)                    :: structOut
  
    TYPE(t_subset)              :: sub  
    CHARACTER(len=*), PARAMETER :: routine = modname//':Add_veg_tiss_p_and_field2d'
  
    sub = struct1%Get_subset()
    CALL structOut%Init_op_ov_hlp('op_ov_hlp', 'op_ov_hlp', struct1%owner_model_id, &
      &                           sub%type, sub%nc, sub%nb, .TRUE., TRIM(struct1%element_unit), 'veg_p')
  
    ! run calculation
    IF (sub%type == ON_CHUNK) THEN
      IF (SIZE(field2d,2) > 1) CALL finish(routine, 'Dimension mismatch, field2d with nblks > 1; calc with REAL pool%element(:)')
      structOut%carbon          = struct1%carbon      + field2d(:,1)
      structOut%nitrogen        = struct1%nitrogen    + field2d(:,1)
      structOut%phosphorus      = struct1%phosphorus  + field2d(:,1)
      structOut%carbon13        = struct1%carbon13    + field2d(:,1)
      structOut%carbon14        = struct1%carbon14    + field2d(:,1)
      structOut%nitrogen15      = struct1%nitrogen15  + field2d(:,1)
    ELSE
      structOut%domain_c        = struct1%domain_c    + field2d(:,:)
      structOut%domain_n        = struct1%domain_n    + field2d(:,:)
      structOut%domain_p        = struct1%domain_p    + field2d(:,:)
      structOut%domain_c13      = struct1%domain_c13  + field2d(:,:)
      structOut%domain_c14      = struct1%domain_c14  + field2d(:,:)
      structOut%domain_n15      = struct1%domain_n15  + field2d(:,:)
    END IF

    structOut%leaf            = struct1%leaf        + field2d(:,:)
    structOut%fine_root       = struct1%fine_root   + field2d(:,:)
    structOut%coarse_root     = struct1%coarse_root + field2d(:,:)
    structOut%sap_wood        = struct1%sap_wood    + field2d(:,:)
    structOut%heart_wood      = struct1%heart_wood  + field2d(:,:)
    structOut%labile          = struct1%labile      + field2d(:,:)
    structOut%reserve         = struct1%reserve     + field2d(:,:)
    structOut%fruit           = struct1%fruit       + field2d(:,:)
    structOut%seed_bed        = struct1%seed_bed    + field2d(:,:)
  
  END FUNCTION Add_veg_tiss_p_and_field2d
  !--------------------------------------------------------------------------------------------------------------------------------
  ! veg_tissues_pool + veg_tissues_pool        | ON_CHUNK & ON_DOMAIN
  !--------------------------------------------------------------------------------------------------------------------------------
  FUNCTION Add_veg_tiss_p_and_veg_tiss_p(struct1, struct2) RESULT(structOut)
  
    CLASS(t_veg_tissues_pool),      INTENT(in)  :: struct1, struct2
    TYPE(t_veg_tissues_pool)                    :: structOut
  
    TYPE(t_subset)              :: sub  
    CHARACTER(len=*), PARAMETER :: routine = modname//':Add_veg_tiss_p_and_veg_tiss_p'
  
    sub = struct1%Get_subset()
    CALL structOut%Init_op_ov_hlp('op_ov_hlp', 'op_ov_hlp', struct1%owner_model_id, &
      &                           sub%type, sub%nc, sub%nb, .TRUE., TRIM(struct1%element_unit), 'veg_p')
  
    ! run calculation
    IF (sub%type == ON_CHUNK) THEN
      structOut%carbon          = struct1%carbon      + struct2%carbon
      structOut%nitrogen        = struct1%nitrogen    + struct2%nitrogen  
      structOut%phosphorus      = struct1%phosphorus  + struct2%phosphorus
      structOut%carbon13        = struct1%carbon13    + struct2%carbon13  
      structOut%carbon14        = struct1%carbon14    + struct2%carbon14  
      structOut%nitrogen15      = struct1%nitrogen15  + struct2%nitrogen15
    ELSE
      structOut%domain_c        = struct1%domain_c    + struct2%domain_c  
      structOut%domain_n        = struct1%domain_n    + struct2%domain_n  
      structOut%domain_p        = struct1%domain_p    + struct2%domain_p  
      structOut%domain_c13      = struct1%domain_c13  + struct2%domain_c13
      structOut%domain_c14      = struct1%domain_c14  + struct2%domain_c14
      structOut%domain_n15      = struct1%domain_n15  + struct2%domain_n15
    END IF

    structOut%leaf            = struct1%leaf        + struct2%leaf       
    structOut%fine_root       = struct1%fine_root   + struct2%fine_root  
    structOut%coarse_root     = struct1%coarse_root + struct2%coarse_root
    structOut%sap_wood        = struct1%sap_wood    + struct2%sap_wood   
    structOut%heart_wood      = struct1%heart_wood  + struct2%heart_wood 
    structOut%labile          = struct1%labile      + struct2%labile     
    structOut%reserve         = struct1%reserve     + struct2%reserve    
    structOut%fruit           = struct1%fruit       + struct2%fruit      
    structOut%seed_bed        = struct1%seed_bed    + struct2%seed_bed   
  
  END FUNCTION Add_veg_tiss_p_and_veg_tiss_p
  
  !--------------------------------------------------------------------------------------------------------------------------------
  ! - minus
  !--------------------------------------------------------------------------------------------------------------------------------
  !--------------------------------------------------------------------------------------------------------------------------------
  ! veg_tissues_pool - scalar        | ON_CHUNK & ON_DOMAIN
  !--------------------------------------------------------------------------------------------------------------------------------
  FUNCTION Subtract_veg_tiss_p_minus_scalar(struct1, scalar) RESULT(structOut)
  
    CLASS(t_veg_tissues_pool),      INTENT(in)  :: struct1
    REAL(wp),                       INTENT(in)  :: scalar
    TYPE(t_veg_tissues_pool)                    :: structOut

    TYPE(t_subset)              :: sub  
    CHARACTER(len=*), PARAMETER :: routine = modname//':Subtract_veg_tiss_p_minus_scalar'

    sub = struct1%Get_subset()
    CALL structOut%Init_op_ov_hlp('op_ov_hlp', 'op_ov_hlp', struct1%owner_model_id, &
      &                           sub%type, sub%nc, sub%nb, .TRUE., TRIM(struct1%element_unit), 'veg_p')

    ! run calculation
    IF (sub%type == ON_CHUNK) THEN
      structOut%carbon          = struct1%carbon      - scalar
      structOut%nitrogen        = struct1%nitrogen    - scalar
      structOut%phosphorus      = struct1%phosphorus  - scalar
      structOut%carbon13        = struct1%carbon13    - scalar
      structOut%carbon14        = struct1%carbon14    - scalar
      structOut%nitrogen15      = struct1%nitrogen15  - scalar
    ELSE
      structOut%domain_c        = struct1%domain_c    - scalar
      structOut%domain_n        = struct1%domain_n    - scalar
      structOut%domain_p        = struct1%domain_p    - scalar
      structOut%domain_c13      = struct1%domain_c13  - scalar
      structOut%domain_c14      = struct1%domain_c14  - scalar
      structOut%domain_n15      = struct1%domain_n15  - scalar
    END IF

    structOut%leaf            = struct1%leaf        - scalar
    structOut%fine_root       = struct1%fine_root   - scalar
    structOut%coarse_root     = struct1%coarse_root - scalar
    structOut%sap_wood        = struct1%sap_wood    - scalar
    structOut%heart_wood      = struct1%heart_wood  - scalar
    structOut%labile          = struct1%labile      - scalar
    structOut%reserve         = struct1%reserve     - scalar
    structOut%fruit           = struct1%fruit       - scalar
    structOut%seed_bed        = struct1%seed_bed    - scalar
  
  END FUNCTION Subtract_veg_tiss_p_minus_scalar
  !--------------------------------------------------------------------------------------------------------------------------------
  ! veg_tissues_pool - field1d(:)        | ON_CHUNK
  !--------------------------------------------------------------------------------------------------------------------------------
  FUNCTION Subtract_veg_tiss_p_minus_field1d(struct1, field1d) RESULT(structOut)
  
    CLASS(t_veg_tissues_pool),      INTENT(in)  :: struct1
    REAL(wp),                       INTENT(in)  :: field1d(:)
    TYPE(t_veg_tissues_pool)                    :: structOut
  
    TYPE(t_subset)              :: sub  
    CHARACTER(len=*), PARAMETER :: routine = modname//':Subtract_veg_tiss_p_minus_field1d'
  
    sub = struct1%Get_subset()
    CALL structOut%Init_op_ov_hlp('op_ov_hlp', 'op_ov_hlp', struct1%owner_model_id, &
      &                           sub%type, sub%nc, sub%nb, .TRUE., TRIM(struct1%element_unit), 'veg_p')
  
    ! run calculation
    IF (sub%type == ON_CHUNK) THEN
      structOut%carbon          = struct1%carbon      - field1d(:)
      structOut%nitrogen        = struct1%nitrogen    - field1d(:)
      structOut%phosphorus      = struct1%phosphorus  - field1d(:)
      structOut%carbon13        = struct1%carbon13    - field1d(:)
      structOut%carbon14        = struct1%carbon14    - field1d(:)
      structOut%nitrogen15      = struct1%nitrogen15  - field1d(:)
    ELSE
      CALL finish(routine, 'Dimension mismatch, cannot calc with 1d-field and jsb_var real2d on domain')
    END IF

    structOut%leaf            = struct1%leaf        - field1d(:)
    structOut%fine_root       = struct1%fine_root   - field1d(:)
    structOut%coarse_root     = struct1%coarse_root - field1d(:)
    structOut%sap_wood        = struct1%sap_wood    - field1d(:)
    structOut%heart_wood      = struct1%heart_wood  - field1d(:)
    structOut%labile          = struct1%labile      - field1d(:)
    structOut%reserve         = struct1%reserve     - field1d(:)
    structOut%fruit           = struct1%fruit       - field1d(:)
    structOut%seed_bed        = struct1%seed_bed    - field1d(:)
  
  END FUNCTION Subtract_veg_tiss_p_minus_field1d
  !--------------------------------------------------------------------------------------------------------------------------------
  ! veg_tissues_pool - field2d(:,:)        | ON_CHUNK & ON_DOMAIN | field2d(:,1) & field2d(:,:)
  !--------------------------------------------------------------------------------------------------------------------------------
  FUNCTION Subtract_veg_tiss_p_minus_field2d(struct1, field2d) RESULT(structOut)
  
    CLASS(t_veg_tissues_pool),      INTENT(in)  :: struct1
    REAL(wp),                       INTENT(in)  :: field2d(:,:)
    TYPE(t_veg_tissues_pool)                    :: structOut
  
    TYPE(t_subset)              :: sub  
    CHARACTER(len=*), PARAMETER :: routine = modname//':Subtract_veg_tiss_p_minus_field2d'
  
    sub = struct1%Get_subset()
    CALL structOut%Init_op_ov_hlp('op_ov_hlp', 'op_ov_hlp', struct1%owner_model_id, &
      &                           sub%type, sub%nc, sub%nb, .TRUE., TRIM(struct1%element_unit), 'veg_p')
  
    ! run calculation
    IF (sub%type == ON_CHUNK) THEN
      IF (SIZE(field2d,2) > 1) CALL finish(routine, 'Dimension mismatch, field2d with nblks > 1; calc with REAL pool%element(:)')
      structOut%carbon          = struct1%carbon      - field2d(:,1)
      structOut%nitrogen        = struct1%nitrogen    - field2d(:,1)
      structOut%phosphorus      = struct1%phosphorus  - field2d(:,1)
      structOut%carbon13        = struct1%carbon13    - field2d(:,1)
      structOut%carbon14        = struct1%carbon14    - field2d(:,1)
      structOut%nitrogen15      = struct1%nitrogen15  - field2d(:,1)
    ELSE
      structOut%domain_c        = struct1%domain_c    - field2d(:,:)
      structOut%domain_n        = struct1%domain_n    - field2d(:,:)
      structOut%domain_p        = struct1%domain_p    - field2d(:,:)
      structOut%domain_c13      = struct1%domain_c13  - field2d(:,:)
      structOut%domain_c14      = struct1%domain_c14  - field2d(:,:)
      structOut%domain_n15      = struct1%domain_n15  - field2d(:,:)
    END IF

    structOut%leaf            = struct1%leaf        - field2d(:,:)
    structOut%fine_root       = struct1%fine_root   - field2d(:,:)
    structOut%coarse_root     = struct1%coarse_root - field2d(:,:)
    structOut%sap_wood        = struct1%sap_wood    - field2d(:,:)
    structOut%heart_wood      = struct1%heart_wood  - field2d(:,:)
    structOut%labile          = struct1%labile      - field2d(:,:)
    structOut%reserve         = struct1%reserve     - field2d(:,:)
    structOut%fruit           = struct1%fruit       - field2d(:,:)
    structOut%seed_bed        = struct1%seed_bed    - field2d(:,:)
  
  END FUNCTION Subtract_veg_tiss_p_minus_field2d
  !--------------------------------------------------------------------------------------------------------------------------------
  ! veg_tissues_pool - veg_tissues_pool        | ON_CHUNK & ON_DOMAIN
  !--------------------------------------------------------------------------------------------------------------------------------
  FUNCTION Subtract_veg_tiss_p_minus_veg_tiss_p(struct1, struct2) RESULT(structOut)
  
    CLASS(t_veg_tissues_pool),      INTENT(in)  :: struct1, struct2
    TYPE(t_veg_tissues_pool)                    :: structOut
  
    TYPE(t_subset)              :: sub  
    CHARACTER(len=*), PARAMETER :: routine = modname//':Subtract_veg_tiss_p_minus_veg_tiss_p'
  
    sub = struct1%Get_subset()
    CALL structOut%Init_op_ov_hlp('op_ov_hlp', 'op_ov_hlp', struct1%owner_model_id, &
      &                           sub%type, sub%nc, sub%nb, .TRUE., TRIM(struct1%element_unit), 'veg_p')
  
    ! run calculation
    IF (sub%type == ON_CHUNK) THEN
      structOut%carbon          = struct1%carbon      - struct2%carbon
      structOut%nitrogen        = struct1%nitrogen    - struct2%nitrogen  
      structOut%phosphorus      = struct1%phosphorus  - struct2%phosphorus
      structOut%carbon13        = struct1%carbon13    - struct2%carbon13  
      structOut%carbon14        = struct1%carbon14    - struct2%carbon14  
      structOut%nitrogen15      = struct1%nitrogen15  - struct2%nitrogen15
    ELSE
      structOut%domain_c        = struct1%domain_c    - struct2%domain_c  
      structOut%domain_n        = struct1%domain_n    - struct2%domain_n  
      structOut%domain_p        = struct1%domain_p    - struct2%domain_p  
      structOut%domain_c13      = struct1%domain_c13  - struct2%domain_c13
      structOut%domain_c14      = struct1%domain_c14  - struct2%domain_c14
      structOut%domain_n15      = struct1%domain_n15  - struct2%domain_n15
    END IF

    structOut%leaf            = struct1%leaf        - struct2%leaf       
    structOut%fine_root       = struct1%fine_root   - struct2%fine_root  
    structOut%coarse_root     = struct1%coarse_root - struct2%coarse_root
    structOut%sap_wood        = struct1%sap_wood    - struct2%sap_wood   
    structOut%heart_wood      = struct1%heart_wood  - struct2%heart_wood 
    structOut%labile          = struct1%labile      - struct2%labile     
    structOut%reserve         = struct1%reserve     - struct2%reserve    
    structOut%fruit           = struct1%fruit       - struct2%fruit      
    structOut%seed_bed        = struct1%seed_bed    - struct2%seed_bed   
  
  END FUNCTION Subtract_veg_tiss_p_minus_veg_tiss_p
  
  !--------------------------------------------------------------------------------------------------------------------------------
  ! * multiply
  !--------------------------------------------------------------------------------------------------------------------------------
  !--------------------------------------------------------------------------------------------------------------------------------
  ! veg_tissues_pool * scalar        | ON_CHUNK & ON_DOMAIN
  !--------------------------------------------------------------------------------------------------------------------------------
  FUNCTION Multiply_veg_tiss_p_with_scalar(struct1, scalar) RESULT(structOut)
  
    CLASS(t_veg_tissues_pool),      INTENT(in)  :: struct1
    REAL(wp),                       INTENT(in)  :: scalar
    TYPE(t_veg_tissues_pool)                    :: structOut

    TYPE(t_subset)              :: sub  
    CHARACTER(len=*), PARAMETER :: routine = modname//':Multiply_veg_tiss_p_with_scalar'

    sub = struct1%Get_subset()
    CALL structOut%Init_op_ov_hlp('op_ov_hlp', 'op_ov_hlp', struct1%owner_model_id, &
      &                           sub%type, sub%nc, sub%nb, .TRUE., TRIM(struct1%element_unit), 'veg_p')

    ! run calculation
    IF (sub%type == ON_CHUNK) THEN
      structOut%carbon          = struct1%carbon      * scalar
      structOut%nitrogen        = struct1%nitrogen    * scalar
      structOut%phosphorus      = struct1%phosphorus  * scalar
      structOut%carbon13        = struct1%carbon13    * scalar
      structOut%carbon14        = struct1%carbon14    * scalar
      structOut%nitrogen15      = struct1%nitrogen15  * scalar
    ELSE
      structOut%domain_c        = struct1%domain_c    * scalar
      structOut%domain_n        = struct1%domain_n    * scalar
      structOut%domain_p        = struct1%domain_p    * scalar
      structOut%domain_c13      = struct1%domain_c13  * scalar
      structOut%domain_c14      = struct1%domain_c14  * scalar
      structOut%domain_n15      = struct1%domain_n15  * scalar
    END IF

    structOut%leaf            = struct1%leaf        * scalar
    structOut%fine_root       = struct1%fine_root   * scalar
    structOut%coarse_root     = struct1%coarse_root * scalar
    structOut%sap_wood        = struct1%sap_wood    * scalar
    structOut%heart_wood      = struct1%heart_wood  * scalar
    structOut%labile          = struct1%labile      * scalar
    structOut%reserve         = struct1%reserve     * scalar
    structOut%fruit           = struct1%fruit       * scalar
    structOut%seed_bed        = struct1%seed_bed    * scalar
  
  END FUNCTION Multiply_veg_tiss_p_with_scalar
  !--------------------------------------------------------------------------------------------------------------------------------
  ! scalar * veg_tissues_pool        | ON_CHUNK & ON_DOMAIN
  !--------------------------------------------------------------------------------------------------------------------------------
  FUNCTION Multiply_scalar_with_veg_tiss_p(scalar, struct1) RESULT(structOut)
  
    REAL(wp),                       INTENT(in)  :: scalar
    CLASS(t_veg_tissues_pool),      INTENT(in)  :: struct1
    TYPE(t_veg_tissues_pool)                    :: structOut
  
    TYPE(t_subset)              :: sub  
    CHARACTER(len=*), PARAMETER :: routine = modname//':Multiply_scalar_with_veg_tiss_p'
  
    sub = struct1%Get_subset()
    CALL structOut%Init_op_ov_hlp('op_ov_hlp', 'op_ov_hlp', struct1%owner_model_id, &
      &                           sub%type, sub%nc, sub%nb, .TRUE., TRIM(struct1%element_unit), 'veg_p')

    ! run calculation
    IF (sub%type == ON_CHUNK) THEN
      structOut%carbon          = scalar * struct1%carbon    
      structOut%nitrogen        = scalar * struct1%nitrogen  
      structOut%phosphorus      = scalar * struct1%phosphorus
      structOut%carbon13        = scalar * struct1%carbon13  
      structOut%carbon14        = scalar * struct1%carbon14  
      structOut%nitrogen15      = scalar * struct1%nitrogen15
    ELSE
      structOut%domain_c        = scalar * struct1%domain_c  
      structOut%domain_n        = scalar * struct1%domain_n  
      structOut%domain_p        = scalar * struct1%domain_p  
      structOut%domain_c13      = scalar * struct1%domain_c13
      structOut%domain_c14      = scalar * struct1%domain_c14
      structOut%domain_n15      = scalar * struct1%domain_n15
    END IF

    structOut%leaf            = scalar * struct1%leaf       
    structOut%fine_root       = scalar * struct1%fine_root  
    structOut%coarse_root     = scalar * struct1%coarse_root
    structOut%sap_wood        = scalar * struct1%sap_wood   
    structOut%heart_wood      = scalar * struct1%heart_wood 
    structOut%labile          = scalar * struct1%labile     
    structOut%reserve         = scalar * struct1%reserve    
    structOut%fruit           = scalar * struct1%fruit      
    structOut%seed_bed        = scalar * struct1%seed_bed   
  
  END FUNCTION Multiply_scalar_with_veg_tiss_p
  !--------------------------------------------------------------------------------------------------------------------------------
  ! veg_tissues_pool * field1d(:)        | ON_CHUNK
  !--------------------------------------------------------------------------------------------------------------------------------
  FUNCTION Multiply_veg_tiss_p_with_field1d(struct1, field1d) RESULT(structOut)
  
    CLASS(t_veg_tissues_pool),      INTENT(in)  :: struct1
    REAL(wp),                       INTENT(in)  :: field1d(:)
    TYPE(t_veg_tissues_pool)                    :: structOut
  
    TYPE(t_subset)              :: sub  
    CHARACTER(len=*), PARAMETER :: routine = modname//':Multiply_veg_tiss_p_with_field1d'
  
    sub = struct1%Get_subset()
    CALL structOut%Init_op_ov_hlp('op_ov_hlp', 'op_ov_hlp', struct1%owner_model_id, &
      &                           sub%type, sub%nc, sub%nb, .TRUE., TRIM(struct1%element_unit), 'veg_p')
  
    ! run calculation
    IF (sub%type == ON_CHUNK) THEN
      structOut%carbon          = struct1%carbon      * field1d(:)
      structOut%nitrogen        = struct1%nitrogen    * field1d(:)
      structOut%phosphorus      = struct1%phosphorus  * field1d(:)
      structOut%carbon13        = struct1%carbon13    * field1d(:)
      structOut%carbon14        = struct1%carbon14    * field1d(:)
      structOut%nitrogen15      = struct1%nitrogen15  * field1d(:)
    ELSE
      CALL finish(routine, 'Dimension mismatch, cannot calc with 1d-field and jsb_var real2d on domain')
    END IF

    structOut%leaf            = struct1%leaf        * field1d(:)
    structOut%fine_root       = struct1%fine_root   * field1d(:)
    structOut%coarse_root     = struct1%coarse_root * field1d(:)
    structOut%sap_wood        = struct1%sap_wood    * field1d(:)
    structOut%heart_wood      = struct1%heart_wood  * field1d(:)
    structOut%labile          = struct1%labile      * field1d(:)
    structOut%reserve         = struct1%reserve     * field1d(:)
    structOut%fruit           = struct1%fruit       * field1d(:)
    structOut%seed_bed        = struct1%seed_bed    * field1d(:)
  
  END FUNCTION Multiply_veg_tiss_p_with_field1d
  !--------------------------------------------------------------------------------------------------------------------------------
  ! field1d(:) * veg_tissues_pool        | ON_CHUNK
  !--------------------------------------------------------------------------------------------------------------------------------
  FUNCTION Multiply_field1d_with_veg_tiss_p(field1d, struct1) RESULT(structOut)
  
    REAL(wp),                       INTENT(in)  :: field1d(:)
    CLASS(t_veg_tissues_pool),      INTENT(in)  :: struct1
    TYPE(t_veg_tissues_pool)                    :: structOut
  
    TYPE(t_subset)              :: sub  
    CHARACTER(len=*), PARAMETER :: routine = modname//':Multiply_field1d_with_veg_tiss_p'
  
    sub = struct1%Get_subset()
    CALL structOut%Init_op_ov_hlp('op_ov_hlp', 'op_ov_hlp', struct1%owner_model_id, &
      &                           sub%type, sub%nc, sub%nb, .TRUE., TRIM(struct1%element_unit), 'veg_p')
  
    ! run calculation
    IF (sub%type == ON_CHUNK) THEN
      structOut%carbon          = field1d(:) * struct1%carbon    
      structOut%nitrogen        = field1d(:) * struct1%nitrogen  
      structOut%phosphorus      = field1d(:) * struct1%phosphorus
      structOut%carbon13        = field1d(:) * struct1%carbon13  
      structOut%carbon14        = field1d(:) * struct1%carbon14  
      structOut%nitrogen15      = field1d(:) * struct1%nitrogen15
    ELSE
      CALL finish(routine, 'Dimension mismatch, cannot calc with 1d-field and jsb_var real2d on domain')
    END IF

    structOut%leaf            = field1d(:) * struct1%leaf       
    structOut%fine_root       = field1d(:) * struct1%fine_root  
    structOut%coarse_root     = field1d(:) * struct1%coarse_root
    structOut%sap_wood        = field1d(:) * struct1%sap_wood   
    structOut%heart_wood      = field1d(:) * struct1%heart_wood 
    structOut%labile          = field1d(:) * struct1%labile     
    structOut%reserve         = field1d(:) * struct1%reserve    
    structOut%fruit           = field1d(:) * struct1%fruit      
    structOut%seed_bed        = field1d(:) * struct1%seed_bed   
  
  END FUNCTION Multiply_field1d_with_veg_tiss_p
  !--------------------------------------------------------------------------------------------------------------------------------
  ! veg_tissues_pool * field2d(:,:)        | ON_CHUNK & ON_DOMAIN | field2d(:,1) & field2d(:,:)
  !--------------------------------------------------------------------------------------------------------------------------------
  FUNCTION Multiply_veg_tiss_p_with_field2d(struct1, field2d) RESULT(structOut)
  
    CLASS(t_veg_tissues_pool),      INTENT(in)  :: struct1
    REAL(wp),                       INTENT(in)  :: field2d(:,:)
    TYPE(t_veg_tissues_pool)                    :: structOut
  
    TYPE(t_subset)              :: sub  
    CHARACTER(len=*), PARAMETER :: routine = modname//':Multiply_veg_tiss_p_with_field2d'
  
    sub = struct1%Get_subset()
    CALL structOut%Init_op_ov_hlp('op_ov_hlp', 'op_ov_hlp', struct1%owner_model_id, &
      &                           sub%type, sub%nc, sub%nb, .TRUE., TRIM(struct1%element_unit), 'veg_p')
  
    ! run calculation
    IF (sub%type == ON_CHUNK) THEN
      IF (SIZE(field2d,2) > 1) CALL finish(routine, 'Dimension mismatch, field2d with nblks > 1; calc with REAL pool%element(:)')
      structOut%carbon          = struct1%carbon      * field2d(:,1)
      structOut%nitrogen        = struct1%nitrogen    * field2d(:,1)
      structOut%phosphorus      = struct1%phosphorus  * field2d(:,1)
      structOut%carbon13        = struct1%carbon13    * field2d(:,1)
      structOut%carbon14        = struct1%carbon14    * field2d(:,1)
      structOut%nitrogen15      = struct1%nitrogen15  * field2d(:,1)
    ELSE
      structOut%domain_c        = struct1%domain_c    * field2d(:,:)
      structOut%domain_n        = struct1%domain_n    * field2d(:,:)
      structOut%domain_p        = struct1%domain_p    * field2d(:,:)
      structOut%domain_c13      = struct1%domain_c13  * field2d(:,:)
      structOut%domain_c14      = struct1%domain_c14  * field2d(:,:)
      structOut%domain_n15      = struct1%domain_n15  * field2d(:,:)
    END IF

    structOut%leaf            = struct1%leaf        * field2d(:,:)
    structOut%fine_root       = struct1%fine_root   * field2d(:,:)
    structOut%coarse_root     = struct1%coarse_root * field2d(:,:)
    structOut%sap_wood        = struct1%sap_wood    * field2d(:,:)
    structOut%heart_wood      = struct1%heart_wood  * field2d(:,:)
    structOut%labile          = struct1%labile      * field2d(:,:)
    structOut%reserve         = struct1%reserve     * field2d(:,:)
    structOut%fruit           = struct1%fruit       * field2d(:,:)
    structOut%seed_bed        = struct1%seed_bed    * field2d(:,:)
  
  END FUNCTION Multiply_veg_tiss_p_with_field2d
  !--------------------------------------------------------------------------------------------------------------------------------
  ! field2d(:,:) * veg_tissues_pool        | ON_CHUNK & ON_DOMAIN | field2d(:,1) & field2d(:,:)
  !--------------------------------------------------------------------------------------------------------------------------------
  FUNCTION Multiply_field2d_with_veg_tiss_p(field2d, struct1) RESULT(structOut)
  
    REAL(wp),                       INTENT(in)  :: field2d(:,:)
    CLASS(t_veg_tissues_pool),      INTENT(in)  :: struct1
    TYPE(t_veg_tissues_pool)                    :: structOut
  
    TYPE(t_subset)              :: sub
    CHARACTER(len=*), PARAMETER :: routine = modname//':Multiply_field2d_with_veg_tiss_p'
  
    sub = struct1%Get_subset()
    CALL structOut%Init_op_ov_hlp('op_ov_hlp', 'op_ov_hlp', struct1%owner_model_id, &
      &                           sub%type, sub%nc, sub%nb, .TRUE., TRIM(struct1%element_unit), 'veg_p')
  
    ! run calculation
    IF (sub%type == ON_CHUNK) THEN
      IF (SIZE(field2d,2) > 1) CALL finish(routine, 'Dimension mismatch, field2d with nblks > 1; calc with REAL pool%element(:)')
      structOut%carbon          = field2d(:,1) * struct1%carbon    
      structOut%nitrogen        = field2d(:,1) * struct1%nitrogen  
      structOut%phosphorus      = field2d(:,1) * struct1%phosphorus
      structOut%carbon13        = field2d(:,1) * struct1%carbon13  
      structOut%carbon14        = field2d(:,1) * struct1%carbon14  
      structOut%nitrogen15      = field2d(:,1) * struct1%nitrogen15
    ELSE
      structOut%domain_c        = field2d(:,:) * struct1%domain_c  
      structOut%domain_n        = field2d(:,:) * struct1%domain_n  
      structOut%domain_p        = field2d(:,:) * struct1%domain_p  
      structOut%domain_c13      = field2d(:,:) * struct1%domain_c13
      structOut%domain_c14      = field2d(:,:) * struct1%domain_c14
      structOut%domain_n15      = field2d(:,:) * struct1%domain_n15
    END IF

    structOut%leaf            = field2d(:,:) * struct1%leaf       
    structOut%fine_root       = field2d(:,:) * struct1%fine_root  
    structOut%coarse_root     = field2d(:,:) * struct1%coarse_root
    structOut%sap_wood        = field2d(:,:) * struct1%sap_wood   
    structOut%heart_wood      = field2d(:,:) * struct1%heart_wood 
    structOut%labile          = field2d(:,:) * struct1%labile     
    structOut%reserve         = field2d(:,:) * struct1%reserve    
    structOut%fruit           = field2d(:,:) * struct1%fruit      
    structOut%seed_bed        = field2d(:,:) * struct1%seed_bed   
  
  END FUNCTION Multiply_field2d_with_veg_tiss_p
  !--------------------------------------------------------------------------------------------------------------------------------
  ! veg_tissues_pool * veg_tissues_pool        | ON_CHUNK & ON_DOMAIN
  !--------------------------------------------------------------------------------------------------------------------------------
  FUNCTION Multiply_veg_tiss_p_with_veg_tiss_p(struct1, struct2) RESULT(structOut)
  
    CLASS(t_veg_tissues_pool),      INTENT(in)  :: struct1, struct2
    TYPE(t_veg_tissues_pool)                    :: structOut
  
    TYPE(t_subset)              :: sub  
    CHARACTER(len=*), PARAMETER :: routine = modname//':Multiply_veg_tiss_p_with_veg_tiss_p'
  
    sub = struct1%Get_subset()
    CALL structOut%Init_op_ov_hlp('op_ov_hlp', 'op_ov_hlp', struct1%owner_model_id, &
      &                           sub%type, sub%nc, sub%nb, .TRUE., TRIM(struct1%element_unit), 'veg_p')
  
    ! run calculation
    IF (sub%type == ON_CHUNK) THEN
      structOut%carbon          = struct1%carbon      * struct2%carbon
      structOut%nitrogen        = struct1%nitrogen    * struct2%nitrogen  
      structOut%phosphorus      = struct1%phosphorus  * struct2%phosphorus
      structOut%carbon13        = struct1%carbon13    * struct2%carbon13  
      structOut%carbon14        = struct1%carbon14    * struct2%carbon14  
      structOut%nitrogen15      = struct1%nitrogen15  * struct2%nitrogen15
    ELSE
      structOut%domain_c        = struct1%domain_c    * struct2%domain_c  
      structOut%domain_n        = struct1%domain_n    * struct2%domain_n  
      structOut%domain_p        = struct1%domain_p    * struct2%domain_p  
      structOut%domain_c13      = struct1%domain_c13  * struct2%domain_c13
      structOut%domain_c14      = struct1%domain_c14  * struct2%domain_c14
      structOut%domain_n15      = struct1%domain_n15  * struct2%domain_n15
    END IF

    structOut%leaf            = struct1%leaf        * struct2%leaf       
    structOut%fine_root       = struct1%fine_root   * struct2%fine_root  
    structOut%coarse_root     = struct1%coarse_root * struct2%coarse_root
    structOut%sap_wood        = struct1%sap_wood    * struct2%sap_wood   
    structOut%heart_wood      = struct1%heart_wood  * struct2%heart_wood 
    structOut%labile          = struct1%labile      * struct2%labile     
    structOut%reserve         = struct1%reserve     * struct2%reserve    
    structOut%fruit           = struct1%fruit       * struct2%fruit      
    structOut%seed_bed        = struct1%seed_bed    * struct2%seed_bed   
  
  END FUNCTION Multiply_veg_tiss_p_with_veg_tiss_p
  
  !--------------------------------------------------------------------------------------------------------------------------------
  ! / divide
  !--------------------------------------------------------------------------------------------------------------------------------
  !--------------------------------------------------------------------------------------------------------------------------------
  ! veg_tissues_pool / scalar        | ON_CHUNK & ON_DOMAIN
  !--------------------------------------------------------------------------------------------------------------------------------
  FUNCTION Div_veg_tiss_p_by_scalar(struct1, scalar) RESULT(structOut)
  
    CLASS(t_veg_tissues_pool),      INTENT(in)  :: struct1
    REAL(wp),                       INTENT(in)  :: scalar
    TYPE(t_veg_tissues_pool)                    :: structOut

    TYPE(t_subset)              :: sub  
    CHARACTER(len=*), PARAMETER :: routine = modname//':Div_veg_tiss_p_by_scalar'

    sub = struct1%Get_subset()
    CALL structOut%Init_op_ov_hlp('op_ov_hlp', 'op_ov_hlp', struct1%owner_model_id, &
      &                           sub%type, sub%nc, sub%nb, .TRUE., TRIM(struct1%element_unit), 'veg_p')

    ! run calculation
    IF (sub%type == ON_CHUNK) THEN
      structOut%carbon          = struct1%carbon      / scalar
      structOut%nitrogen        = struct1%nitrogen    / scalar
      structOut%phosphorus      = struct1%phosphorus  / scalar
      structOut%carbon13        = struct1%carbon13    / scalar
      structOut%carbon14        = struct1%carbon14    / scalar
      structOut%nitrogen15      = struct1%nitrogen15  / scalar
    ELSE
      structOut%domain_c        = struct1%domain_c    / scalar
      structOut%domain_n        = struct1%domain_n    / scalar
      structOut%domain_p        = struct1%domain_p    / scalar
      structOut%domain_c13      = struct1%domain_c13  / scalar
      structOut%domain_c14      = struct1%domain_c14  / scalar
      structOut%domain_n15      = struct1%domain_n15  / scalar
    END IF

    structOut%leaf            = struct1%leaf        / scalar
    structOut%fine_root       = struct1%fine_root   / scalar
    structOut%coarse_root     = struct1%coarse_root / scalar
    structOut%sap_wood        = struct1%sap_wood    / scalar
    structOut%heart_wood      = struct1%heart_wood  / scalar
    structOut%labile          = struct1%labile      / scalar
    structOut%reserve         = struct1%reserve     / scalar
    structOut%fruit           = struct1%fruit       / scalar
    structOut%seed_bed        = struct1%seed_bed    / scalar
  
  END FUNCTION Div_veg_tiss_p_by_scalar
  !--------------------------------------------------------------------------------------------------------------------------------
  ! veg_tissues_pool / field1d(:)        | ON_CHUNK
  !--------------------------------------------------------------------------------------------------------------------------------
  FUNCTION Div_veg_tiss_p_by_field1d(struct1, field1d) RESULT(structOut)
  
    CLASS(t_veg_tissues_pool),      INTENT(in)  :: struct1
    REAL(wp),                       INTENT(in)  :: field1d(:)
    TYPE(t_veg_tissues_pool)                    :: structOut
  
    TYPE(t_subset)              :: sub  
    CHARACTER(len=*), PARAMETER :: routine = modname//':Div_veg_tiss_p_by_field1d'
  
    sub = struct1%Get_subset()
    CALL structOut%Init_op_ov_hlp('op_ov_hlp', 'op_ov_hlp', struct1%owner_model_id, &
      &                           sub%type, sub%nc, sub%nb, .TRUE., TRIM(struct1%element_unit), 'veg_p')
  
    ! run calculation
    IF (sub%type == ON_CHUNK) THEN
      structOut%carbon          = struct1%carbon      / field1d(:)
      structOut%nitrogen        = struct1%nitrogen    / field1d(:)
      structOut%phosphorus      = struct1%phosphorus  / field1d(:)
      structOut%carbon13        = struct1%carbon13    / field1d(:)
      structOut%carbon14        = struct1%carbon14    / field1d(:)
      structOut%nitrogen15      = struct1%nitrogen15  / field1d(:)
    ELSE
      CALL finish(routine, 'Dimension mismatch, cannot calc with 1d-field and jsb_var real2d on domain')
    END IF

    structOut%leaf            = struct1%leaf        / field1d(:)
    structOut%fine_root       = struct1%fine_root   / field1d(:)
    structOut%coarse_root     = struct1%coarse_root / field1d(:)
    structOut%sap_wood        = struct1%sap_wood    / field1d(:)
    structOut%heart_wood      = struct1%heart_wood  / field1d(:)
    structOut%labile          = struct1%labile      / field1d(:)
    structOut%reserve         = struct1%reserve     / field1d(:)
    structOut%fruit           = struct1%fruit       / field1d(:)
    structOut%seed_bed        = struct1%seed_bed    / field1d(:)
  
  END FUNCTION Div_veg_tiss_p_by_field1d
  !--------------------------------------------------------------------------------------------------------------------------------
  ! veg_tissues_pool / field2d(:,:)        | ON_CHUNK & ON_DOMAIN | field2d(:,1) & field2d(:,:)
  !--------------------------------------------------------------------------------------------------------------------------------
  FUNCTION Div_veg_tiss_p_by_field2d(struct1, field2d) RESULT(structOut)
  
    CLASS(t_veg_tissues_pool),      INTENT(in)  :: struct1
    REAL(wp),                       INTENT(in)  :: field2d(:,:)
    TYPE(t_veg_tissues_pool)                    :: structOut
  
    TYPE(t_subset)              :: sub  
    CHARACTER(len=*), PARAMETER :: routine = modname//':Div_veg_tiss_p_by_field2d'
  
    sub = struct1%Get_subset()
    CALL structOut%Init_op_ov_hlp('op_ov_hlp', 'op_ov_hlp', struct1%owner_model_id, &
      &                           sub%type, sub%nc, sub%nb, .TRUE., TRIM(struct1%element_unit), 'veg_p')
  
    ! run calculation
    IF (sub%type == ON_CHUNK) THEN
      IF (SIZE(field2d,2) > 1) CALL finish(routine, 'Dimension mismatch, field2d with nblks > 1; calc with REAL pool%element(:)')
      structOut%carbon          = struct1%carbon      / field2d(:,1)
      structOut%nitrogen        = struct1%nitrogen    / field2d(:,1)
      structOut%phosphorus      = struct1%phosphorus  / field2d(:,1)
      structOut%carbon13        = struct1%carbon13    / field2d(:,1)
      structOut%carbon14        = struct1%carbon14    / field2d(:,1)
      structOut%nitrogen15      = struct1%nitrogen15  / field2d(:,1)
    ELSE
      structOut%domain_c        = struct1%domain_c    / field2d(:,:)
      structOut%domain_n        = struct1%domain_n    / field2d(:,:)
      structOut%domain_p        = struct1%domain_p    / field2d(:,:)
      structOut%domain_c13      = struct1%domain_c13  / field2d(:,:)
      structOut%domain_c14      = struct1%domain_c14  / field2d(:,:)
      structOut%domain_n15      = struct1%domain_n15  / field2d(:,:)
    END IF

    structOut%leaf            = struct1%leaf        / field2d(:,:)
    structOut%fine_root       = struct1%fine_root   / field2d(:,:)
    structOut%coarse_root     = struct1%coarse_root / field2d(:,:)
    structOut%sap_wood        = struct1%sap_wood    / field2d(:,:)
    structOut%heart_wood      = struct1%heart_wood  / field2d(:,:)
    structOut%labile          = struct1%labile      / field2d(:,:)
    structOut%reserve         = struct1%reserve     / field2d(:,:)
    structOut%fruit           = struct1%fruit       / field2d(:,:)
    structOut%seed_bed        = struct1%seed_bed    / field2d(:,:)
  
  END FUNCTION Div_veg_tiss_p_by_field2d
  !--------------------------------------------------------------------------------------------------------------------------------
  ! veg_tissues_pool / veg_tissues_pool        | ON_CHUNK & ON_DOMAIN
  !--------------------------------------------------------------------------------------------------------------------------------
  FUNCTION Div_veg_tiss_p_by_veg_tiss_p(struct1, struct2) RESULT(structOut)
  
    CLASS(t_veg_tissues_pool),      INTENT(in)  :: struct1, struct2
    TYPE(t_veg_tissues_pool)                    :: structOut
  
    TYPE(t_subset)              :: sub  
    CHARACTER(len=*), PARAMETER :: routine = modname//':Div_veg_tiss_p_by_veg_tiss_p'
  
    sub = struct1%Get_subset()
    CALL structOut%Init_op_ov_hlp('op_ov_hlp', 'op_ov_hlp', struct1%owner_model_id, &
      &                           sub%type, sub%nc, sub%nb, .TRUE., TRIM(struct1%element_unit), 'veg_p')
  
    ! run calculation
    IF (sub%type == ON_CHUNK) THEN
      structOut%carbon          = struct1%carbon      / struct2%carbon
      structOut%nitrogen        = struct1%nitrogen    / struct2%nitrogen  
      structOut%phosphorus      = struct1%phosphorus  / struct2%phosphorus
      structOut%carbon13        = struct1%carbon13    / struct2%carbon13  
      structOut%carbon14        = struct1%carbon14    / struct2%carbon14  
      structOut%nitrogen15      = struct1%nitrogen15  / struct2%nitrogen15
    ELSE
      structOut%domain_c        = struct1%domain_c    / struct2%domain_c  
      structOut%domain_n        = struct1%domain_n    / struct2%domain_n  
      structOut%domain_p        = struct1%domain_p    / struct2%domain_p  
      structOut%domain_c13      = struct1%domain_c13  / struct2%domain_c13
      structOut%domain_c14      = struct1%domain_c14  / struct2%domain_c14
      structOut%domain_n15      = struct1%domain_n15  / struct2%domain_n15
    END IF

    structOut%leaf            = struct1%leaf        / struct2%leaf       
    structOut%fine_root       = struct1%fine_root   / struct2%fine_root  
    structOut%coarse_root     = struct1%coarse_root / struct2%coarse_root
    structOut%sap_wood        = struct1%sap_wood    / struct2%sap_wood   
    structOut%heart_wood      = struct1%heart_wood  / struct2%heart_wood 
    structOut%labile          = struct1%labile      / struct2%labile     
    structOut%reserve         = struct1%reserve     / struct2%reserve    
    structOut%fruit           = struct1%fruit       / struct2%fruit      
    structOut%seed_bed        = struct1%seed_bed    / struct2%seed_bed   
  
  END FUNCTION Div_veg_tiss_p_by_veg_tiss_p

END MODULE mo_veg_memory_class
