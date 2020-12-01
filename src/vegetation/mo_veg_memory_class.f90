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
  USE mo_jsb_pool_class,         ONLY: t_jsb_pool, ELEM_C, ELEM_N, ELEM_P, ELEM_C13, ELEM_C14, ELEM_N15
  USE mo_jsb_lct_class,          ONLY: LAND_TYPE, VEG_TYPE


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
  TYPE, EXTENDS(t_jsb_pool) :: t_veg_elements_pool
    ! element var
    TYPE(t_jsb_var_real2d), POINTER :: carbon         ! element
    TYPE(t_jsb_var_real2d), POINTER :: nitrogen       ! element
    TYPE(t_jsb_var_real2d), POINTER :: phosphorus     ! element
    TYPE(t_jsb_var_real2d), POINTER :: carbon13       ! element
    TYPE(t_jsb_var_real2d), POINTER :: carbon14       ! element
    TYPE(t_jsb_var_real2d), POINTER :: nitrogen15     ! element
  CONTAINS
    !! infrastructure
    PROCEDURE :: Init                   => t_veg_elements_pool_init
    PROCEDURE :: Get_element_name_by_id => t_veg_elements_pool_get_element_name_by_id
    PROCEDURE :: Init_elem_var_elem_p   => Init_element_var_veg_elements_pool           ! used in t_veg_tissues_pool operator overload routines to init temporary objects
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
  TYPE, EXTENDS(t_jsb_pool) :: t_veg_tissues_pool
    ! element variables
    TYPE(t_jsb_var_real2d), POINTER :: carbon
    TYPE(t_jsb_var_real2d), POINTER :: nitrogen
    TYPE(t_jsb_var_real2d), POINTER :: phosphorus
    TYPE(t_jsb_var_real2d), POINTER :: carbon13
    TYPE(t_jsb_var_real2d), POINTER :: carbon14
    TYPE(t_jsb_var_real2d), POINTER :: nitrogen15
    ! element pools
    TYPE(t_veg_elements_pool), POINTER :: leaf             !< pool with elements of leaf
    TYPE(t_veg_elements_pool), POINTER :: fine_root        !< pool with elements of fine_root
    TYPE(t_veg_elements_pool), POINTER :: coarse_root      !< pool with elements of coarse_root
    TYPE(t_veg_elements_pool), POINTER :: sap_wood         !< pool with elements of sap_wood
    TYPE(t_veg_elements_pool), POINTER :: heart_wood       !< pool with elements of heart_wood
    TYPE(t_veg_elements_pool), POINTER :: labile           !< pool with elements of labile
    TYPE(t_veg_elements_pool), POINTER :: reserve          !< pool with elements of reserve
    TYPE(t_veg_elements_pool), POINTER :: fruit            !< pool with elements of fruit
    TYPE(t_veg_elements_pool), POINTER :: seed_bed         !< pool with elements of seed_bed
  CONTAINS
    !! infrastructure
    PROCEDURE :: Init                     => t_veg_tissues_pool_init
    PROCEDURE :: Get_element_name_by_id   => t_veg_tissues_pool_get_element_name_by_id
    PROCEDURE :: Init_tiss_pool_for_op_ov => Create_veg_tissues_pool_temporary_object_for_op_ov    ! create 
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
    ! TYPE(t_jsb_var_real2d), POINTER :: carbon, nitrogen, phosphorus, carbon13, carbon14, nitrogen15
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

  !> list with the element_id of all elements as defined above
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

    ! basic plant 2D
    TYPE(t_jsb_var_real2d)          :: &   
                                    height            , &            !< plant height [m]
                                    diameter          , &            !< plant diameter [m]
                                    dens_ind          , &            !< density of individuals [# m-2]
                                    lai               , &            !< actual leaf area index [m2 m-2]
                                    target_cn_leaf    , &            !< target foliar C:N ratio [mol C mol-1 N]
                                    target_np_leaf    , &            !< target foliar N:P ratio [mol N mol-1 P]
                                    mean_leaf_age     , &            !< average foliage age in days [days]
                                    cohort_age        , &            !< age of the cohort in years, used with veg_dynamics_scheme cohort [years]
                                    f_n_demand        , &            !< relative plant demand for N
                                    f_p_demand        , &            !< relative plant demand for P
                                    beta_sinklim_ps   , &            !< scaling factor to account for a direct sink limitation contraint on photosynthesis [unitless]
                                    dphi              , &            !< plant-soil water potential gradient [MPa]
                                    growth_req_n      , &            !< moles N required for a unit of C growth under current allocation fractions [mol N m-2 -1]
                                    growth_req_p      , &            !< moles P required for a unit of N growth under current allocation fractions [mol P m-2 -1]
                                    leaf_cn_direction , &            !< indicates direction of change in leaf CN ratio to maximise NPP [unitless]
                                    target_lai        , &            !< LAI implied from sap wood area [m2 m-2]
                                    kstar_labile                     !< labile pool turnover rate assuming only temperature and moisture constraints, no stoichiometry [day-1]

    ! used as logical: 0 = FALSE / 1 = TRUE
    TYPE(t_jsb_var_real2d)          :: &
                                    do_cohort_harvest                !< logical to identify whether the cohort may be harvested, in veg_dynamics_scheme cohort [logical]

    ! basic plant 3D 
    TYPE(t_jsb_var_real3d)          :: &
                                    root_fraction_sl                 !< fractions of root per layer [unitless]  DIMENSION(soil_layer_sb)

    ! target C:N and N:P ratios for various non-leaf tissues (mol/mol)
    TYPE(t_jsb_var_real2d)          :: &
                                    target_cn_fine_root       , &   !< 
                                    target_cn_coarse_root     , &   !< 
                                    target_cn_sap_wood        , &   !< 
                                    target_cn_fruit           , &   !< 
                                    target_np_fine_root       , &   !< 
                                    target_np_coarse_root     , &   !< 
                                    target_np_sap_wood        , &   !< 
                                    target_np_fruit                 !< 

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

    ! fluxes
    TYPE(t_jsb_var_real2d)          :: &   
                                    maint_respiration_pot         , & !< potential maintenance respiration [µmol CO2 m-2 s-1]
                                    maint_respiration             , & !< maintenance respiration [µmol CO2 m-2 s-1]
                                    maint_respiration_c13         , & !< maintenance respiration [µmol 13CO2 m-2 s-1]
                                    maint_respiration_c14         , & !< maintenance respiration [µmol 14CO2 m-2 s-1]
                                    growth_respiration            , & !< growth respiration [µmol CO2 m-2 s-1]
                                    growth_respiration_c13        , & !< growth respiration [µmol 13CO2 m-2 s-1]
                                    growth_respiration_c14        , & !< growth respiration [µmol 14CO2 m-2 s-1]
                                    n_transform_respiration       , & !< respiration associated with N transformation [µmol CO2 m-2 s-1]
                                    n_fixation_respiration        , & !< respiration associated with N fixation [µmol CO2 m-2 s-1]
                                    n_processing_respiration      , & !< respiration associated with N processing [µmol CO2 m-2 s-1]
                                    n_processing_respiration_c13  , & !< respiration associated with N processing [µmol 13CO2 m-2 s-1]
                                    n_processing_respiration_c14  , & !< respiration associated with N processing [µmol 14CO2 m-2 s-1]
                                    npp                           , & !< net primary production [µmol CO2 m-2 s-1]
                                    npp_c13                       , & !< net primary production [µmol 13CO2 m-2 s-1]
                                    npp_c14                       , & !< net primary production [µmol 14CO2 m-2 s-1]
                                    uptake_nh4                    , & !< NH4 uptake from soil [µmol N m-2 s-1]
                                    uptake_nh4_n15                , & !< 15NH4 uptake from soil [µmol N m-2 s-1]
                                    uptake_no3                    , & !< NO3 uptake from soil [µmol N m-2 s-1]
                                    uptake_no3_n15                , & !< 15NO3 uptake from soil [µmol N m-2 s-1]
                                    uptake_po4                    , & !< PO4 uptake from soil [µmol P m-2 s-1]
                                    n_fixation                    , & !< symbiontic N fixation [µmol N m-2 s-1]
                                    n_fixation_n15                , & !< symbiontic N fixation [µmol 15N m-2 s-1]
                                    recycling_leaf_n              , & !< net flux of N from leaf to labile (senescence and maintenance) [mol m-2 timestep-1]
                                    recycling_leaf_n15            , & !< net flux of 15N from leaf to labile (senescence and maintenance) [mol m-2 timestep-1]
                                    recycling_leaf_p              , & !< net flux of P from leaf to labile (senescence and maintenance) [mol m-2 timestep-1]
                                    recycling_fine_root_n         , & !< net flux of N from fine roots to labile (senescence and maintenance) [mol m-2 timestep-1]
                                    recycling_fine_root_n15       , & !< net flux of 15N from fine roots to labile (senescence and maintenance) [mol m-2 timestep-1]
                                    recycling_fine_root_p         , & !< net flux of P from fine roots to labile (senescence and maintenance) [mol m-2 timestep-1]
                                    recycling_heart_wood_n        , & !< net flux of N from sapwood to labile (senescence and maintenance) [mol m-2 timestep-1]
                                    recycling_heart_wood_n15      , & !< net flux of 15N from sapwood to labile (senescence and maintenance) [mol m-2 timestep-1]
                                    recycling_heart_wood_p        , & !< net flux of P from sapwood to labile (senescence and maintenance) [mol m-2 timestep-1]
                                    net_growth                        !< npp minus turnover [µmol CO2 m-2 s-1]

    ! everything else 2D
    TYPE(t_jsb_var_real2d)          :: &   
                                    fmaint_rate_root          , & !< current N-specific maintentance respiration [~µmol C mol-1 N s-1]
                                    unit_npp                  , & !< incremental change in NPP [µmol CO2 m-2 s-1]
                                    unit_transpiration        , & !< incremental change in transpiration [kg m-2 s-1]
                                    unit_uptake_n_pot         , & !< incremental change in potential N uptake per unit root carbon [µmol N m-2 s-1]
                                    unit_uptake_p_pot         , & !< incremental change in potential P uptake per unit root carbon [µmol P m-2 s-1]
                                    unit_uptake_n_act         , & !< incremental change in actual N uptake per unit root carbon [µmol N m-2 s-1]
                                    unit_uptake_p_act         , & !< incremental change in actual P uptake per unit root carbon [µmol P m-2 s-1]
                                    delta_dens_ind            , & !< change in density of individuals [# m-2 timestep-1]
                                    cost_n_uptake_root        , & !< actual cost for N uptake by roots [mol C mol-1 N]
                                    cost_n_fixation           , & !< actual cost for N uptake through fixation [mol C mol-1 N]
                                    mortality_rate            , & !< background_mort_rate_tree / grass [timestep-1]
                                    k1_opt                    , & !< coefficient for optimal leaf:root ratio
                                    k2_opt                    , & !< coefficient for optimal leaf:root ratio
                                    leaf2sapwood_mass_ratio   , & !< leaf to shoot mass ratio [unitless]
                                    leaf2root_mass_ratio          !< leaf to root mass ratio [unitless]

    !> daytime averages
    !! 
    TYPE(t_jsb_var_real3d)          :: &   ! DIMENSION(ncanopy)
                                    fleaf_sunlit_daytime_cl          !< fraction of layer which is sunlit average of the previous day [unitless]
    TYPE(t_jsb_var_real2d)          :: &   
                                    t_air_daytime             , &    !< air temperature average of the previous day [K]
                                    press_srf_daytime         , &    !< air pressure at surface average of the previous day [Pa]
                                    co2_mixing_ratio_daytime  , &    !< co2 mixing ratio of the air average of the previous day [co2 in ppm]
                                    ga_daytime                , &    !< aerodynamic conductance average of the previous day [m s-1]
                                    beta_sinklim_ps_daytime   , &    !< scaling factor to account for soil moisture constraints on photosynthesis average of the previous day [unitless]
                                    t_jmax_opt_daytime               !< optimal temperature for electron transport of photosynthesis average of the previous day [deg C]

    !> accumlated over current daytime
    !! 
    TYPE(t_jsb_var_real3d)          :: &   ! DIMENSION(ncanopy)
                                    fleaf_sunlit_daytime_dacc_cl     !< fraction of layer which is sunlit value accumulated over current day [unitless]
    TYPE(t_jsb_var_real2d)          :: &   
                                    t_air_daytime_dacc            , & !< air temperature average of the previous day [K]
                                    press_srf_daytime_dacc        , & !< air pressure at surface average of the previous day [Pa]
                                    co2_mixing_ratio_daytime_dacc , & !< co2 mixing ratio of the air average of the previous day [co2 in ppm]
                                    ga_daytime_dacc               , & !< aerodynamic conductance average of the previous day [m s-1]
                                    beta_sinklim_ps_daytime_dacc  , & !< scaling factor to account for soil moisture constraints on photosynthesis average of the previous day [unitless]
                                    t_jmax_opt_daytime_dacc           !< optimal temperature for electron transport of photosynthesis value accumulated over current day [deg C]
    TYPE(t_jsb_var_real2d)          :: &   
                                    daytime_counter                   !< number of timesteps of current day with light available [#]

    !> moving averages
    !! 
    ! 1.1 memory of conditions for the labile pool dynamics
    TYPE(t_jsb_var_real2d)          :: &   
                                    gpp_tlabile_mavg               , &  !< gpp  [µmol m-2 s-1]
                                    maint_respiration_tlabile_mavg , &  !< maintenance respiration [µmol m-2 s-1]
                                    growth_req_n_tlabile_mavg      , &  !< moles N required for a unit of C growth under current allocation fractions time-averaged [mol N m-2 s-1]
                                    growth_req_p_tlabile_mavg           !< moles P required for a unit of N growth under current allocation fractions time-averaged [mol P m-2 s-1]
    ! 1.2 memory of conditions for phenology calculations
    TYPE(t_jsb_var_real2d)          :: &   
                                    t_air_tphen_mavg                    !< air temperature time-averaged [K]
    ! 1.3 memory for plant uptake (average demand over a couple of days to get rid of diurnal C cycle)
    TYPE(t_jsb_var_real2d)          :: &   
                                    npp_tuptake_mavg               , &  !< npp of the last week [µmol m-2 s-1]
                                    demand_uptake_n_tuptake_mavg   , &  !< scalar of Vmax for root N uptake time-averaged [unitless]
                                    demand_uptake_p_tuptake_mavg   , &  !< scalar of Vmax for root P uptake time-averaged [unitless]
                                    growth_req_n_tuptake_mavg      , &  !< moles N required for a unit of C growth under current allocation fractions time-averaged [mol N m-2 s-1]
                                    growth_req_p_tuptake_mavg      , &  !< moles P required for a unit of N growth under current allocation fractions time-averaged [mol P m-2 s-1]
                                    cost_n_fixation_talloc_mavg    , &  !< moving average of N fixation cost [mol C mol-1 N]
                                    exudation_c_tmyc_mavg               !< moving average of C allocation to mycorrhizal fungi [mol C]
    ! 1.4 memory of conditions for photosynthesis for the optimisation of foliar N fractions 
    TYPE(t_jsb_var_real3d)          :: &   ! DIMENSION(ncanopy)
                                    fleaf_sunlit_tfrac_mavg_cl          !< fraction of layer which is sunlit time-averaged [unitless]
    TYPE(t_jsb_var_real2d)          :: &   
                                    t_air_tfrac_mavg               , &  !< air temperature time-averaged [K]
                                    press_srf_tfrac_mavg           , &  !< air pressure at surface time-averaged [Pa]
                                    co2_mixing_ratio_tfrac_mavg    , &  !< co2 mixing ratio of the air time-averaged [co2 in ppm]
                                    ga_tfrac_mavg                  , &  !< aerodynamic conductance time-averaged [m s-1]
                                    beta_sinklim_ps_tfrac_mavg          !< scaling factor to account for soil moisture constraints on photosynthesis time-averaged [unitless]
    ! 1.5 memory of conditions for photosynthesis for the optimisation of foliar N  
    TYPE(t_jsb_var_real3d)          :: &   ! DIMENSION(ncanopy)
                                    fleaf_sunlit_tcnl_mavg_cl           !< fraction of layer which is sunlit time-averaged [unitless]
    TYPE(t_jsb_var_real2d)          :: &   
                                    t_air_tcnl_mavg                , &  !< air temperature time-averaged [K]
                                    press_srf_tcnl_mavg            , &  !< air pressure at surface time-averaged [Pa]
                                    co2_mixing_ratio_tcnl_mavg     , &  !< co2 mixing ratio of the air time-averaged [co2 in ppm]
                                    ga_tcnl_mavg                   , &  !< aerodynamic conductance time-averaged [m s-1]
                                    beta_sinklim_ps_tcnl_mavg      , &  !< scaling factor to account for soil moisture constraints on conductance time-averaged [unitless]
                                    uptake_n_tcnl_mavg             , &  !< root N uptake time-averaged [µmol N m-2 s-1]
                                    growth_cn_tcnl_mavg            , &  !< new growth CN ratio time-averaged [mol C mol-1 N]
                                    growth_np_tcnl_mavg            , &  !< new growth NP ratio time-averaged [mol N mol-1 P]
                                    npp_tcnl_mavg                  , &  !< NPP time-averaged [µmol CO2 m-2 s-1]
                                    fmaint_rate_tcnl_mavg          , &  !< maintenance respiration rate per unit N for roots
                                    labile_nitrogen_tcnl_mavg      , &  !< labile nitrogen vegetation pool time-averaged [mol N m-2]
                                    labile_carbon_tcnl_mavg        , &  !< labile carbon vegetation pool time-averaged [mol C m-2]
                                    labile_phosphorus_tcnl_mavg         !< labile phosphorus vegetation pool time-averaged [mol P m-2]

    ! 1.6 memory of conditions for photosynthesis for the optimisation of biomass allocation fractions
    TYPE(t_jsb_var_real2d)          :: &   
                                    npp_talloc_mavg                , &  !< NPP [µmol m-2 s-1]
                                    n_fixation_talloc_mavg         , &  !< N fixation time-averaged [µmol N m-2 s-1]
                                    transpiration_talloc_mavg      , &  !< transpiration time-averaged [kg m-2 s-1]
                                    unit_npp_talloc_mavg           , &  !< rate of change in NPP with  increment in leaf mass [µmol C mol-1]
                                    unit_transpiration_talloc_mavg , &  !< incremental change in transpiration time-averaged [kg m-2 s-1]
                                    unit_uptake_n_talloc_mavg      , &  !< actual rate of N uptkae per unit root [µmol N mol-1]
                                    unit_uptake_p_talloc_mavg      , &  !< actual rate of P uptake per unit root [µmol P mol-1]
                                    dphi_talloc_mavg               , &  !< leaf-soil water potential gradient [MPa]
                                    growth_cn_talloc_mavg          , &  !< plant growth cn ratio, equivalent to 1/growth_req_n time-averaged [mol C mol-1 N]
                                    growth_cp_talloc_mavg          , &  !< plant growth cp ratio, equivalent to 1/(growthreq_n*growth_req_p) time-averaged [mol C mol-1 P]
                                    growth_np_talloc_mavg          , &  !< plant growth np ratio, equivalent to 1/(growth_req_p) time-averaged [mol N mol-1 P]
                                    labile_carbon_talloc_mavg      , &  !< labile carbon pool needed for allometry calculations time-averaged [mol C m-2]
                                    labile_nitrogen_talloc_mavg    , &  !< labile nitrogen pool needed for allometry calculations time-averaged [mol N m-2]
                                    labile_phosphorus_talloc_mavg  , &  !< labile phosphorus pool needed for allometry calculations time-averaged [mol P m-2]
                                    reserve_carbon_talloc_mavg     , &  !< reserve carbon needed for optimal allocation time-averaged [mol C m-2]
                                    reserve_nitrogen_talloc_mavg   , &  !< reserve nitrogen needed for optimal allocation time-averaged [mol N m-2]
                                    reserve_phosphorus_talloc_mavg , &  !< reserve phosphorus needed for optimal allocation time-averaged [mol P m-2]
                                    w_root_lim_talloc_mavg              !< root water limitation below which root allocation starts responding time-averaged
                                    
                                    
    ! 1.7 memory over root life time 
    TYPE(t_jsb_var_real2d)          :: &   
                                    leaf2root_troot_mavg           , &  !< leaf to root mass ratio time-averaged [unitless]
                                    unit_npp_troot_mavg            , &  !< rate of change in NPP with increment in leaf mass [µmol C mol-1]
                                    unit_uptake_n_pot_troot_mavg   , &  !< rate of N uptkae per unit root [µmol N mol-1]
                                    unit_uptake_p_pot_troot_mavg   , &  !< rate of P uptake per unit root [µmol P mol-1]
                                    fmaint_rate_troot_mavg              !< moving average of N-specific maintentance respiration [~µmol C mol-1 N s-1]
    ! 1.8 memory over timescale of vegetation dynamics 
    TYPE(t_jsb_var_real2d)          :: &   
                                    an_boc_tvegdyn_mavg            , &  !< net photosynthesis at bottom of canopy, lowest canopy layer existing time-averaged [µmol m-2 s-1]
                                    net_growth_tvegdyn_mavg        , &  !< net growth, npp minus turnover [µmol C m-2 s-1]
                                    lai_tvegdyn_mavg                    !< lai time-averaged [m2 m-2]
    ! 1.9 memory of variables for respiration accclimation 
    TYPE(t_jsb_var_real2d)          :: &   
                                    t_air_tacclim_mavg             , &  !< air temperature for respiration acclimation time-averaged [K]
                                    t_soil_root_tacclim_mavg            !< root temperature for respiration acclimation time-averaged [deg K]
    ! 1.10 memory of certain variables, temperature, ...
    TYPE(t_jsb_var_real2d)          :: &   
                                    t_air_week_mavg                , &  !< air temperature time-averaged [K]
                                    t_air_month_mavg               , &  !< air temperature time-averaged [K]
                                    t_jmax_opt_mavg                     !< optimal temperature for electron transport of photosynthesis time-averaged [deg C]

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
  
    this%element_unit = element_unit
  
    !> carbon
    IF (One_of(ELEM_C, element_list_loc) > 0) THEN
      CALL this%Add_element('carbon', 'C', ELEM_C, dim=2)
      i = this%Find_element_position(ELEM_C)
      IF (i > 0) THEN
        SELECT TYPE (selector => this%element_list(i)%p)
        CLASS IS (t_jsb_var_real2d)
          this%carbon => selector
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
          this%nitrogen => selector
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
          this%phosphorus => selector
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
          this%carbon13 => selector
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
          this%carbon14 => selector
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
          this%nitrogen15 => selector
        END SELECT
      END IF
    END IF

  END SUBROUTINE t_veg_elements_pool_init

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
  
    this%element_unit = element_unit
  
    !> carbon
    IF (One_of(ELEM_C, element_list_loc) > 0) THEN
      CALL this%Add_element('carbon', 'C', ELEM_C, dim=2)
      i = this%Find_element_position(ELEM_C)
      IF (i > 0) THEN
        SELECT TYPE (selector => this%element_list(i)%p)
        CLASS IS (t_jsb_var_real2d)
          this%carbon => selector
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
          this%nitrogen => selector
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
          this%phosphorus => selector
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
          this%carbon13 => selector
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
          this%carbon14 => selector
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
          this%nitrogen15 => selector
        END SELECT
      END IF
    END IF

  END SUBROUTINE t_veg_tissues_pool_init

  !-----------------------------------------------------------------------------------------------------
  !> create and init an object of the VEG_ tissues pool "t_veg_tissues_pool"
  !! init and add all its objects of the "t_veg_elements_pool" with "pool_object%Init()"
  !! "object => Create_veg_tissues_pool()"
  !! used for VEG_ memory init
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
  !> create an elements_pool that is part of a tissue pool
  !! used for operator overload routines of tissue pools
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE Create_veg_tissues_pool_temporary_object_for_op_ov(this, name, shortname, unit, pool_suffix)

    CLASS(t_veg_tissues_pool),  INTENT(inout) :: this
    CHARACTER(LEN=*),           INTENT(in)    :: name
    CHARACTER(LEN=*),           INTENT(in)    :: shortname
    CHARACTER(LEN=*),           INTENT(in)    :: unit
    CHARACTER(LEN=*), OPTIONAL, INTENT(in)    :: pool_suffix

    TYPE(t_veg_elements_pool),  POINTER       :: veg_elements_pool
    CHARACTER(LEN=:),           ALLOCATABLE   :: suffix_loc

    ! pool name suffix
    suffix_loc = ''
    CALL assign_if_present_allocatable(suffix_loc, pool_suffix)
    IF (suffix_loc /= '') suffix_loc = '_' // suffix_loc

    ! tissue pool init element var
    CALL this%Init(name, shortname, l_elements=.TRUE., element_unit=unit)
    
    ! element pools allocate & init & add to tissue pool
    ! "shortname" is used for model output
    ! "element_list" is optional and can be used to select specific elements (not needed when all elements may be added)
    ! leaf
    ALLOCATE(veg_elements_pool)
    CALL veg_elements_pool%Init(name='leaf_pool'//suffix_loc, shortname='leaf', l_elements=.TRUE., element_unit=unit)
    CALL this%Add_pool(veg_elements_pool)
    this%leaf => veg_elements_pool
    NULLIFY(veg_elements_pool)
    ! fine_root
    ALLOCATE(veg_elements_pool)
    CALL veg_elements_pool%Init(name='fine_root_pool'//suffix_loc, shortname='fine_root', l_elements=.TRUE., element_unit=unit)
    CALL this%Add_pool(veg_elements_pool)
    this%fine_root => veg_elements_pool
    NULLIFY(veg_elements_pool)
    ! coarse_root
    ALLOCATE(veg_elements_pool)
    CALL veg_elements_pool%Init(name='coarse_root_pool'//suffix_loc, shortname='coarse_root', l_elements=.TRUE., element_unit=unit)
    CALL this%Add_pool(veg_elements_pool)
    this%coarse_root => veg_elements_pool
    NULLIFY(veg_elements_pool)
    ! sap_wood
    ALLOCATE(veg_elements_pool)
    CALL veg_elements_pool%Init(name='sap_wood_pool'//suffix_loc, shortname='sap_wood', l_elements=.TRUE., element_unit=unit)
    CALL this%Add_pool(veg_elements_pool)
    this%sap_wood => veg_elements_pool
    NULLIFY(veg_elements_pool)
    ! heart_wood
    ALLOCATE(veg_elements_pool)
    CALL veg_elements_pool%Init(name='heart_wood_pool'//suffix_loc, shortname='heart_wood', l_elements=.TRUE., element_unit=unit)
    CALL this%Add_pool(veg_elements_pool)
    this%heart_wood => veg_elements_pool
    NULLIFY(veg_elements_pool)
    ! labile
    ALLOCATE(veg_elements_pool)
    CALL veg_elements_pool%Init(name='labile_pool'//suffix_loc, shortname='labile', l_elements=.TRUE., element_unit=unit)
    CALL this%Add_pool(veg_elements_pool)
    this%labile => veg_elements_pool
    NULLIFY(veg_elements_pool)
    ! reserve
    ALLOCATE(veg_elements_pool)
    CALL veg_elements_pool%Init(name='reserve_pool'//suffix_loc, shortname='reserve', l_elements=.TRUE., element_unit=unit)
    CALL this%Add_pool(veg_elements_pool)
    this%reserve => veg_elements_pool
    NULLIFY(veg_elements_pool)
    ! fruit
    ALLOCATE(veg_elements_pool)
    CALL veg_elements_pool%Init(name='fruit_pool'//suffix_loc, shortname='fruit', l_elements=.TRUE., element_unit=unit)
    CALL this%Add_pool(veg_elements_pool)
    this%fruit => veg_elements_pool
    NULLIFY(veg_elements_pool)
    ! seed_bed
    ALLOCATE(veg_elements_pool)
    CALL veg_elements_pool%Init(name='seed_bed_pool'//suffix_loc, shortname='seed_bed', l_elements=.TRUE., element_unit=unit)
    CALL this%Add_pool(veg_elements_pool)
    this%seed_bed => veg_elements_pool
    NULLIFY(veg_elements_pool)

  END SUBROUTINE Create_veg_tissues_pool_temporary_object_for_op_ov

  !-----------------------------------------------------------------------------------------------------
  !> initialize an elements_pool that is part of a tissue pool
  !! used for operator overload routines of tissue pools
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE Init_element_var_veg_elements_pool(this, struct1)

    CLASS(t_veg_elements_pool),  INTENT(inout)  :: this
    CLASS(t_veg_elements_pool),  INTENT(in)     :: struct1

    this%carbon%owner_model_id      = struct1%carbon%owner_model_id
    this%carbon%owner_proc_id       = struct1%carbon%owner_proc_id
    this%carbon%owner_tile_path     = struct1%carbon%owner_tile_path
    ALLOCATE(this%carbon%ptr2d,     source=struct1%carbon%ptr2d)
    this%nitrogen%owner_model_id    = struct1%nitrogen%owner_model_id
    this%nitrogen%owner_proc_id     = struct1%nitrogen%owner_proc_id
    this%nitrogen%owner_tile_path   = struct1%nitrogen%owner_tile_path
    ALLOCATE(this%nitrogen%ptr2d,   source=struct1%nitrogen%ptr2d)
    this%phosphorus%owner_model_id  = struct1%phosphorus%owner_model_id
    this%phosphorus%owner_proc_id   = struct1%phosphorus%owner_proc_id
    this%phosphorus%owner_tile_path = struct1%phosphorus%owner_tile_path
    ALLOCATE(this%phosphorus%ptr2d, source=struct1%phosphorus%ptr2d)
    this%carbon13%owner_model_id    = struct1%carbon13%owner_model_id
    this%carbon13%owner_proc_id     = struct1%carbon13%owner_proc_id
    this%carbon13%owner_tile_path   = struct1%carbon13%owner_tile_path
    ALLOCATE(this%carbon13%ptr2d,   source=struct1%carbon13%ptr2d)
    this%carbon14%owner_model_id    = struct1%carbon14%owner_model_id
    this%carbon14%owner_proc_id     = struct1%carbon14%owner_proc_id
    this%carbon14%owner_tile_path   = struct1%carbon14%owner_tile_path
    ALLOCATE(this%carbon14%ptr2d,   source=struct1%carbon14%ptr2d)
    this%nitrogen15%owner_model_id  = struct1%nitrogen15%owner_model_id
    this%nitrogen15%owner_proc_id   = struct1%nitrogen15%owner_proc_id
    this%nitrogen15%owner_tile_path = struct1%nitrogen15%owner_tile_path
    ALLOCATE(this%nitrogen15%ptr2d, source=struct1%nitrogen15%ptr2d)

  END SUBROUTINE Init_element_var_veg_elements_pool

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

      CALL mem%Add_var(mem, 'height', mem%height, &
            hgrid, vgrid_surface, &
            t_cf('height', 'm', 'plant height'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            loutput=.TRUE., lrestart=.TRUE., &
            initval_r=0.0_wp) ! initvalue
  
      CALL mem%Add_var(mem, 'diameter', mem%diameter, &
            hgrid, vgrid_surface, &
            t_cf('diameter', 'm', 'plant diameter'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            loutput=.TRUE., lrestart=.TRUE., &
            initval_r=0.0_wp) ! initvalue
  
      CALL mem%Add_var(mem, 'dens_ind', mem%dens_ind, &
            hgrid, vgrid_surface, &
            t_cf('dens_ind', '# m-2', 'density of individuals'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            loutput=.TRUE., lrestart=.TRUE., &
            initval_r=0.0_wp) ! initvalue
  
      CALL mem%Add_var(mem, 'lai', mem%lai, &
            hgrid, vgrid_surface, &
            t_cf('lai', 'm2 m-2', 'actual leaf area index'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            loutput=.TRUE., lrestart=.TRUE., &
            initval_r=0.0_wp) ! initvalue
  
      CALL mem%Add_var(mem, 'target_cn_leaf', mem%target_cn_leaf, &
            hgrid, vgrid_surface, &
            t_cf('target_cn_leaf', '', 'target foliar C:N ratio'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            loutput=.TRUE., lrestart=.TRUE., &
            initval_r=0.0_wp) ! initvalue
  
      CALL mem%Add_var(mem, 'target_np_leaf', mem%target_np_leaf, &
            hgrid, vgrid_surface, &
            t_cf('target_np_leaf', '', 'target foliar N:P ratio '), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            loutput=.TRUE., lrestart=.TRUE., &
            initval_r=0.0_wp) ! initvalue
  
      CALL mem%Add_var(mem, 'mean_leaf_age', mem%mean_leaf_age, &
            hgrid, vgrid_surface, &
            t_cf('mean_leaf_age', '', 'average foliage age in days '), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            loutput=.TRUE., lrestart=.TRUE., &
            initval_r=0.0_wp) ! initvalue
  
      CALL mem%Add_var(mem, 'cohort_age', mem%cohort_age, &
            hgrid, vgrid_surface, &
            t_cf('cohort_age', '', 'age of the cohort in years, used with veg_dynamics_scheme cohort'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            loutput=.TRUE., lrestart=.TRUE., &
            initval_r=0.0_wp) ! initvalue
  
      CALL mem%Add_var(mem, 'f_n_demand', mem%f_n_demand, &
            hgrid, vgrid_surface, &
            t_cf('f_n_demand', '', 'relative plant demand for N'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            loutput=.TRUE., lrestart=.TRUE., &
            initval_r=0.0_wp) ! initvalue
  
      CALL mem%Add_var(mem, 'f_p_demand', mem%f_p_demand, &
            hgrid, vgrid_surface, &
            t_cf('f_p_demand', '', 'relative plant demand for P'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            loutput=.TRUE., lrestart=.TRUE., &
            initval_r=0.0_wp) ! initvalue
  
      CALL mem%Add_var(mem, 'beta_sinklim_ps', mem%beta_sinklim_ps, &
            hgrid, vgrid_surface, &
            t_cf('beta_sinklim_ps', '', 'scaling factor to account for a direct sink limitation contraint on photosynthesis'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            loutput=.TRUE., lrestart=.TRUE., &
            initval_r=0.0_wp) ! initvalue
  
      CALL mem%Add_var(mem, 'dphi', mem%dphi, &
            hgrid, vgrid_surface, &
            t_cf('dphi', 'MPa', 'plant-soil water potential gradient'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            loutput=.TRUE., lrestart=.TRUE., &
            initval_r=0.0_wp) ! initvalue
  
      CALL mem%Add_var(mem, 'growth_req_n', mem%growth_req_n, &
            hgrid, vgrid_surface, &
            t_cf('growth_req_n', '', 'moles N required for a unit of C growth under current allocation fractions '), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            loutput=.TRUE., lrestart=.TRUE., &
            initval_r=0.0_wp) ! initvalue
  
      CALL mem%Add_var(mem, 'growth_req_p', mem%growth_req_p, &
            hgrid, vgrid_surface, &
            t_cf('growth_req_p', '', 'moles P required for a unit of N growth under current allocation fractions '), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            loutput=.TRUE., lrestart=.TRUE., &
            initval_r=0.0_wp) ! initvalue
  
      CALL mem%Add_var(mem, 'leaf_cn_direction', mem%leaf_cn_direction, &
            hgrid, vgrid_surface, &
            t_cf('leaf_cn_direction', '', 'indicates direction of change in leaf CN ratio to maximise NPP'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            loutput=.TRUE., lrestart=.TRUE., &
            initval_r=0.0_wp) ! initvalue
  
      CALL mem%Add_var(mem, 'target_lai', mem%target_lai, &
            hgrid, vgrid_surface, &
            t_cf('target_lai', 'm2 m-2', 'LAI implied from sap wood area'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            loutput=.TRUE., lrestart=.TRUE., &
            initval_r=0.0_wp) ! initvalue
  
      CALL mem%Add_var(mem, 'kstar_labile', mem%kstar_labile, &
            hgrid, vgrid_surface, &
            t_cf('kstar_labile', '', 'labile pool turnover rate assuming temperature and moisture constraints, no stoichiometry'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            loutput=.TRUE., lrestart=.TRUE., &
            initval_r=0.0_wp) ! initvalue

      CALL mem%Add_var(mem, 'do_cohort_harvest', mem%do_cohort_harvest, &
            hgrid, vgrid_surface, &
            t_cf('do_cohort_harvest', '', 'logical: whether the cohort may be harvested, in veg_dynamics_scheme cohort'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            loutput=.TRUE., lrestart=.TRUE., &
            initval_r=false) ! initvalue

      CALL mem%Add_var(mem, 'root_fraction_sl', mem%root_fraction_sl, &
            hgrid, vgrid_soil_sb, &
            t_cf('root_fraction_sl', '', 'fractions of root per layer'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            loutput=.TRUE., lrestart=.TRUE., &
            initval_r=0.0_wp) ! initvalue

      CALL mem%Add_var(mem, 'target_cn_fine_root', mem%target_cn_fine_root, &
            hgrid, vgrid_surface, &
            t_cf('target_cn_fine_root', '', ' '), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            loutput=.TRUE., lrestart=.TRUE., &
            initval_r=0.0_wp) ! initvalue
  
      CALL mem%Add_var(mem, 'target_cn_coarse_root', mem%target_cn_coarse_root, &
            hgrid, vgrid_surface, &
            t_cf('target_cn_coarse_root', '', ' '), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            loutput=.TRUE., lrestart=.TRUE., &
            initval_r=0.0_wp) ! initvalue
  
      CALL mem%Add_var(mem, 'target_cn_sap_wood', mem%target_cn_sap_wood, &
            hgrid, vgrid_surface, &
            t_cf('target_cn_sap_wood', '', ' '), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            loutput=.TRUE., lrestart=.TRUE., &
            initval_r=0.0_wp) ! initvalue
  
      CALL mem%Add_var(mem, 'target_cn_fruit', mem%target_cn_fruit, &
            hgrid, vgrid_surface, &
            t_cf('target_cn_fruit', '', ' '), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            loutput=.TRUE., lrestart=.TRUE., &
            initval_r=0.0_wp) ! initvalue
  
      CALL mem%Add_var(mem, 'target_np_fine_root', mem%target_np_fine_root, &
            hgrid, vgrid_surface, &
            t_cf('target_np_fine_root', '', ' '), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            loutput=.TRUE., lrestart=.TRUE., &
            initval_r=0.0_wp) ! initvalue
  
      CALL mem%Add_var(mem, 'target_np_coarse_root', mem%target_np_coarse_root, &
            hgrid, vgrid_surface, &
            t_cf('target_np_coarse_root', '', ' '), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            loutput=.TRUE., lrestart=.TRUE., &
            initval_r=0.0_wp) ! initvalue
  
      CALL mem%Add_var(mem, 'target_np_sap_wood', mem%target_np_sap_wood, &
            hgrid, vgrid_surface, &
            t_cf('target_np_sap_wood', '', ' '), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            loutput=.TRUE., lrestart=.TRUE., &
            initval_r=0.0_wp) ! initvalue
  
      CALL mem%Add_var(mem, 'target_np_fruit', mem%target_np_fruit, &
            hgrid, vgrid_surface, &
            t_cf('target_np_fruit', '', ' '), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            loutput=.TRUE., lrestart=.TRUE., &
            initval_r=0.0_wp) ! initvalue
  
      CALL mem%Add_var(mem, 'leaf_nitrogen_cl', mem%leaf_nitrogen_cl, &
            hgrid, vgrid_canopy, &
            t_cf('leaf_nitrogen_cl', 'mmol m-2 LAI', 'N in canopy layer'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            loutput=.TRUE., lrestart=.TRUE., &
            initval_r=0.0_wp) ! initvalue
  
      CALL mem%Add_var(mem, 'fn_chl_cl', mem%fn_chl_cl, &
            hgrid, vgrid_canopy, &
            t_cf('fn_chl_cl', '', 'fraction of N in chlorophyll, pigments only'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            loutput=.TRUE., lrestart=.TRUE., &
            initval_r=0.0_wp) ! initvalue
  
      CALL mem%Add_var(mem, 'fn_et_cl', mem%fn_et_cl, &
            hgrid, vgrid_canopy, &
            t_cf('fn_et_cl', '', 'fraction of N in electron transport '), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            loutput=.TRUE., lrestart=.TRUE., &
            initval_r=0.0_wp) ! initvalue
  
      CALL mem%Add_var(mem, 'fn_rub_cl', mem%fn_rub_cl, &
            hgrid, vgrid_canopy, &
            t_cf('fn_rub_cl', '', 'fraction of N in Rubisco '), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            loutput=.TRUE., lrestart=.TRUE., &
            initval_r=0.0_wp) ! initvalue
  
      CALL mem%Add_var(mem, 'fn_pepc_cl', mem%fn_pepc_cl, &
            hgrid, vgrid_canopy, &
            t_cf('fn_pepc_cl', '', 'fraction of N in PEP carboylase, C4 plants only'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            loutput=.TRUE., lrestart=.TRUE., &
            initval_r=0.0_wp) ! initvalue
  
      CALL mem%Add_var(mem, 'fn_oth_cl', mem%fn_oth_cl, &
            hgrid, vgrid_canopy, &
            t_cf('fn_oth_cl', '', 'fraction of N not associated with photosynthesis '), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            loutput=.TRUE., lrestart=.TRUE., &
            initval_r=0.0_wp) ! initvalue
  
      CALL mem%Add_var(mem, 'fleaf_sunlit_cl', mem%fleaf_sunlit_cl, &
            hgrid, vgrid_canopy, &
            t_cf('fleaf_sunlit_cl', '', 'fraction of layer which is sunlit'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            loutput=.TRUE., lrestart=.TRUE., &
            initval_r=0.0_wp) ! initvalue
  
      CALL mem%Add_var(mem, 'lai_cl', mem%lai_cl, &
            hgrid, vgrid_canopy, &
            t_cf('lai_cl', 'm2 m-2', 'leaf area in the canopy layer '), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            loutput=.TRUE., lrestart=.TRUE., &
            initval_r=0.0_wp) ! initvalue
  
      CALL mem%Add_var(mem, 'cumm_lai_cl', mem%cumm_lai_cl, &
            hgrid, vgrid_canopy, &
            t_cf('cumm_lai_cl', 'm2 m-2', 'cummulative leaf area above the mid-point of the canopy layer '), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            loutput=.TRUE., lrestart=.TRUE., &
            initval_r=0.0_wp) ! initvalue
  
      CALL mem%Add_var(mem, 'maint_respiration_pot', mem%maint_respiration_pot, &
            hgrid, vgrid_surface, &
            t_cf('maint_respiration_pot', 'µmol CO2 m-2 s-1', 'potential maintenance respiration'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            loutput=.TRUE., lrestart=.TRUE., &
            initval_r=0.0_wp) ! initvalue
  
      CALL mem%Add_var(mem, 'maint_respiration', mem%maint_respiration, &
            hgrid, vgrid_surface, &
            t_cf('maint_respiration', 'µmol CO2 m-2 s-1', 'maintenance respiration'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            loutput=.TRUE., lrestart=.TRUE., &
            initval_r=0.0_wp) ! initvalue
  
      CALL mem%Add_var(mem, 'maint_respiration_c13', mem%maint_respiration_c13, &
            hgrid, vgrid_surface, &
            t_cf('maint_respiration_c13', 'µmol 13CO2 m-2 s-1', 'maintenance respiration'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            loutput=.TRUE., lrestart=.TRUE., &
            initval_r=0.0_wp) ! initvalue
  
      CALL mem%Add_var(mem, 'maint_respiration_c14', mem%maint_respiration_c14, &
            hgrid, vgrid_surface, &
            t_cf('maint_respiration_c14', 'µmol 14CO2 m-2 s-1', 'maintenance respiration'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            loutput=.TRUE., lrestart=.TRUE., &
            initval_r=0.0_wp) ! initvalue
  
      CALL mem%Add_var(mem, 'growth_respiration', mem%growth_respiration, &
            hgrid, vgrid_surface, &
            t_cf('growth_respiration', 'µmol CO2 m-2 s-1', 'growth respiration'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            loutput=.TRUE., lrestart=.TRUE., &
            initval_r=0.0_wp) ! initvalue
  
      CALL mem%Add_var(mem, 'growth_respiration_c13', mem%growth_respiration_c13, &
            hgrid, vgrid_surface, &
            t_cf('growth_respiration_c13', 'µmol 13CO2 m-2 s-1', 'growth respiration'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            loutput=.TRUE., lrestart=.TRUE., &
            initval_r=0.0_wp) ! initvalue
  
      CALL mem%Add_var(mem, 'growth_respiration_c14', mem%growth_respiration_c14, &
            hgrid, vgrid_surface, &
            t_cf('growth_respiration_c14', 'µmol 14CO2 m-2 s-1', 'growth respiration'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            loutput=.TRUE., lrestart=.TRUE., &
            initval_r=0.0_wp) ! initvalue
  
      CALL mem%Add_var(mem, 'n_transform_respiration', mem%n_transform_respiration, &
            hgrid, vgrid_surface, &
            t_cf('n_transform_respiration', 'µmol CO2 m-2 s-1', 'respiration associated with N transformation'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            loutput=.TRUE., lrestart=.TRUE., &
            initval_r=0.0_wp) ! initvalue
  
      CALL mem%Add_var(mem, 'n_fixation_respiration', mem%n_fixation_respiration, &
            hgrid, vgrid_surface, &
            t_cf('n_fixation_respiration', 'µmol CO2 m-2 s-1', 'respiration associated with N fixation'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            loutput=.TRUE., lrestart=.TRUE., &
            initval_r=0.0_wp) ! initvalue
  
      CALL mem%Add_var(mem, 'n_processing_respiration', mem%n_processing_respiration, &
            hgrid, vgrid_surface, &
            t_cf('n_processing_respiration', 'µmol CO2 m-2 s-1', 'respiration associated with N processing'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            loutput=.TRUE., lrestart=.TRUE., &
            initval_r=0.0_wp) ! initvalue
  
      CALL mem%Add_var(mem, 'n_processing_respiration_c13', mem%n_processing_respiration_c13, &
            hgrid, vgrid_surface, &
            t_cf('n_processing_respiration_c13', 'µmol 13CO2 m-2 s-1', 'respiration associated with N processing'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            loutput=.TRUE., lrestart=.TRUE., &
            initval_r=0.0_wp) ! initvalue
  
      CALL mem%Add_var(mem, 'n_processing_respiration_c14', mem%n_processing_respiration_c14, &
            hgrid, vgrid_surface, &
            t_cf('n_processing_respiration_c14', 'µmol 14CO2 m-2 s-1', 'respiration associated with N processing'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            loutput=.TRUE., lrestart=.TRUE., &
            initval_r=0.0_wp) ! initvalue
  
      CALL mem%Add_var(mem, 'npp', mem%npp, &
            hgrid, vgrid_surface, &
            t_cf('npp', 'µmol CO2 m-2 s-1', 'net primary production'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            loutput=.TRUE., lrestart=.TRUE., &
            initval_r=0.0_wp) ! initvalue
  
      CALL mem%Add_var(mem, 'npp_c13', mem%npp_c13, &
            hgrid, vgrid_surface, &
            t_cf('npp_c13', 'µmol 13CO2 m-2 s-1', 'net primary production'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            loutput=.TRUE., lrestart=.TRUE., &
            initval_r=0.0_wp) ! initvalue
  
      CALL mem%Add_var(mem, 'npp_c14', mem%npp_c14, &
            hgrid, vgrid_surface, &
            t_cf('npp_c14', 'µmol 14CO2 m-2 s-1', 'net primary production'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            loutput=.TRUE., lrestart=.TRUE., &
            initval_r=0.0_wp) ! initvalue
  
      CALL mem%Add_var(mem, 'uptake_nh4', mem%uptake_nh4, &
            hgrid, vgrid_surface, &
            t_cf('uptake_nh4', 'µmol N m-2 s-1', 'NH4 uptake from soil'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            loutput=.TRUE., lrestart=.TRUE., &
            initval_r=0.0_wp) ! initvalue
  
      CALL mem%Add_var(mem, 'uptake_nh4_n15', mem%uptake_nh4_n15, &
            hgrid, vgrid_surface, &
            t_cf('uptake_nh4_n15', 'µmol N m-2 s-1', '15NH4 uptake from soil'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            loutput=.TRUE., lrestart=.TRUE., &
            initval_r=0.0_wp) ! initvalue
  
      CALL mem%Add_var(mem, 'uptake_no3', mem%uptake_no3, &
            hgrid, vgrid_surface, &
            t_cf('uptake_no3', 'µmol N m-2 s-1', 'NO3 uptake from soil'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            loutput=.TRUE., lrestart=.TRUE., &
            initval_r=0.0_wp) ! initvalue
  
      CALL mem%Add_var(mem, 'uptake_no3_n15', mem%uptake_no3_n15, &
            hgrid, vgrid_surface, &
            t_cf('uptake_no3_n15', 'µmol N m-2 s-1', '15NO3 uptake from soil'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            loutput=.TRUE., lrestart=.TRUE., &
            initval_r=0.0_wp) ! initvalue

      CALL mem%Add_var(mem, 'uptake_po4', mem%uptake_po4, &
            hgrid, vgrid_surface, &
            t_cf('uptake_po4', 'µmol P m-2 s-1', 'PO4 uptake from soil'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            loutput=.TRUE., lrestart=.TRUE., &
            initval_r=0.0_wp) ! initvalue
  
      CALL mem%Add_var(mem, 'n_fixation', mem%n_fixation, &
            hgrid, vgrid_surface, &
            t_cf('n_fixation', 'µmol N m-2 s-1', 'symbiontic N fixation'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            loutput=.TRUE., lrestart=.TRUE., &
            initval_r=0.0_wp) ! initvalue
  
      CALL mem%Add_var(mem, 'n_fixation_n15', mem%n_fixation_n15, &
            hgrid, vgrid_surface, &
            t_cf('n_fixation_n15', 'µmol 15N m-2 s-1', 'symbiontic N fixation'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            loutput=.TRUE., lrestart=.TRUE., &
            initval_r=0.0_wp) ! initvalue
  
      CALL mem%Add_var(mem, 'recycling_leaf_n', mem%recycling_leaf_n, &
            hgrid, vgrid_surface, &
            t_cf('recycling_leaf_n', 'mol m-2 timestep-1', 'N net flux from leaf to labile (senescence and maintenance)'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            loutput=.TRUE., lrestart=.TRUE., &
            initval_r=0.0_wp) ! initvalue
  
      CALL mem%Add_var(mem, 'recycling_leaf_n15', mem%recycling_leaf_n15, &
            hgrid, vgrid_surface, &
            t_cf('recycling_leaf_n15', 'mol m-2 timestep-1', '15N net flux from leaf to labile (senescence and maintenance)'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            loutput=.TRUE., lrestart=.TRUE., &
            initval_r=0.0_wp) ! initvalue
  
      CALL mem%Add_var(mem, 'recycling_leaf_p', mem%recycling_leaf_p, &
            hgrid, vgrid_surface, &
            t_cf('recycling_leaf_p', 'mol m-2 timestep-1', 'P net flux from leaf to labile (senescence and maintenance)'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            loutput=.TRUE., lrestart=.TRUE., &
            initval_r=0.0_wp) ! initvalue
  
      CALL mem%Add_var(mem, 'recycling_fine_root_n', mem%recycling_fine_root_n, &
            hgrid, vgrid_surface, &
            t_cf('recycling_fine_root_n', 'mol m-2 timestep-1', 'N net flux from fine roots to labile (senescence & maint)'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            loutput=.TRUE., lrestart=.TRUE., &
            initval_r=0.0_wp) ! initvalue
  
      CALL mem%Add_var(mem, 'recycling_fine_root_n15', mem%recycling_fine_root_n15, &
            hgrid, vgrid_surface, &
            t_cf('recycling_fine_root_n15', 'mol m-2 timestep-1', '15N net flux from fine roots to labile (senescence & maint)'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            loutput=.TRUE., lrestart=.TRUE., &
            initval_r=0.0_wp) ! initvalue
  
      CALL mem%Add_var(mem, 'recycling_fine_root_p', mem%recycling_fine_root_p, &
            hgrid, vgrid_surface, &
            t_cf('recycling_fine_root_p', 'mol m-2 timestep-1', 'P net flux from fine roots to labile (senescence & maint)'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            loutput=.TRUE., lrestart=.TRUE., &
            initval_r=0.0_wp) ! initvalue
  
      CALL mem%Add_var(mem, 'recycling_heart_wood_n', mem%recycling_heart_wood_n, &
            hgrid, vgrid_surface, &
            t_cf('recycling_heart_wood_n', 'mol m-2 timestep-1', 'N net flux from sapwood to labile (senescence & maint)'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            loutput=.TRUE., lrestart=.TRUE., &
            initval_r=0.0_wp) ! initvalue
  
      CALL mem%Add_var(mem, 'recycling_heart_wood_n15', mem%recycling_heart_wood_n15, &
            hgrid, vgrid_surface, &
            t_cf('recycling_heart_wood_n15', 'mol m-2 timestep-1', '15N net flux from sapwood to labile (senescence & maint)'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            loutput=.TRUE., lrestart=.TRUE., &
            initval_r=0.0_wp) ! initvalue
  
      CALL mem%Add_var(mem, 'recycling_heart_wood_p', mem%recycling_heart_wood_p, &
            hgrid, vgrid_surface, &
            t_cf('recycling_heart_wood_p', 'mol m-2 timestep-1', 'P net flux from sapwood to labile (senescence & maint)'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            loutput=.TRUE., lrestart=.TRUE., &
            initval_r=0.0_wp) ! initvalue
  
      CALL mem%Add_var(mem, 'net_growth', mem%net_growth, &
            hgrid, vgrid_surface, &
            t_cf('net_growth', 'µmol CO2 m-2 s-1', 'npp minus turnover'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            loutput=.TRUE., lrestart=.TRUE., &
            initval_r=0.0_wp) ! initvalue
  
      CALL mem%Add_var(mem, 'fmaint_rate_root', mem%fmaint_rate_root, &
            hgrid, vgrid_surface, &
            t_cf('fmaint_rate_root', '~µmol C mol-1 N s-1', 'current N-specific maintentance respiration'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            loutput=.TRUE., lrestart=.TRUE., &
            initval_r=0.0_wp) ! initvalue
  
      CALL mem%Add_var(mem, 'unit_npp', mem%unit_npp, &
            hgrid, vgrid_surface, &
            t_cf('unit_npp', '', 'incremental change in NPP. Needed for optimal allocation'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            loutput=.TRUE., lrestart=.TRUE., &
            initval_r=0.0_wp) ! initvalue
  
      CALL mem%Add_var(mem, 'unit_transpiration', mem%unit_transpiration, &
            hgrid, vgrid_surface, &
            t_cf('unit_transpiration', '', 'incremental change in transpiration. Needed for optimal allocation'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            loutput=.TRUE., lrestart=.TRUE., &
            initval_r=0.0_wp) ! initvalue
  
      CALL mem%Add_var(mem, 'unit_uptake_n_pot', mem%unit_uptake_n_pot, &
            hgrid, vgrid_surface, &
            t_cf('unit_uptake_n_pot', '', 'incremental change in potential N uptake / unit root carbon (N uptake calculation)'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            loutput=.TRUE., lrestart=.TRUE., &
            initval_r=0.0_wp) ! initvalue
  
      CALL mem%Add_var(mem, 'unit_uptake_p_pot', mem%unit_uptake_p_pot, &
            hgrid, vgrid_surface, &
            t_cf('unit_uptake_p_pot', '', 'incremental change in potential P uptake per unit root carbon (currently NOT used)'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            loutput=.TRUE., lrestart=.TRUE., &
            initval_r=0.0_wp) ! initvalue
  
      CALL mem%Add_var(mem, 'unit_uptake_n_act', mem%unit_uptake_n_act, &
            hgrid, vgrid_surface, &
            t_cf('unit_uptake_n_act', '', 'incremental change in actual N uptake per unit root carbon (optimal allocation)'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            loutput=.TRUE., lrestart=.TRUE., &
            initval_r=0.0_wp) ! initvalue
  
      CALL mem%Add_var(mem, 'unit_uptake_p_act', mem%unit_uptake_p_act, &
            hgrid, vgrid_surface, &
            t_cf('unit_uptake_p_act', '', 'incremental change in actual P uptake per unit root carbon (optimal allocation)'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            loutput=.TRUE., lrestart=.TRUE., &
            initval_r=0.0_wp) ! initvalue
  
      CALL mem%Add_var(mem, 'delta_dens_ind', mem%delta_dens_ind, &
            hgrid, vgrid_surface, &
            t_cf('delta_dens_ind', '# m-2 timestep-1', 'change in density of individuals'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            loutput=.TRUE., lrestart=.TRUE., &
            initval_r=0.0_wp) ! initvalue
  
      CALL mem%Add_var(mem, 'cost_n_uptake_root', mem%cost_n_uptake_root, &
            hgrid, vgrid_surface, &
            t_cf('cost_n_uptake_root', 'mol C mol-1 N', 'actual cost for N uptake by roots'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            loutput=.TRUE., lrestart=.TRUE., &
            initval_r=0.0_wp) ! initvalue
  
      CALL mem%Add_var(mem, 'cost_n_fixation', mem%cost_n_fixation, &
            hgrid, vgrid_surface, &
            t_cf('cost_n_fixation', 'mol C mol-1 N', 'actual cost for N uptake through fixation'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            loutput=.TRUE., lrestart=.TRUE., &
            initval_r=0.0_wp) ! initvalue
  
      CALL mem%Add_var(mem, 'mortality_rate', mem%mortality_rate, &
            hgrid, vgrid_surface, &
            t_cf('mortality_rate', 'timestep-1', 'background_mort_rate_tree / grass'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            loutput=.TRUE., lrestart=.TRUE., &
            initval_r=0.0_wp) ! initvalue
  
      CALL mem%Add_var(mem, 'k1_opt', mem%k1_opt, &
            hgrid, vgrid_surface, &
            t_cf('k1_opt', '', 'coefficient for optimal leaf:root ratio'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            loutput=.TRUE., lrestart=.TRUE., &
            initval_r=0.0_wp) ! initvalue
  
      CALL mem%Add_var(mem, 'k2_opt', mem%k2_opt, &
            hgrid, vgrid_surface, &
            t_cf('k2_opt', '', 'coefficient for optimal leaf:root ratio'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            loutput=.TRUE., lrestart=.TRUE., &
            initval_r=0.0_wp) ! initvalue
  
      CALL mem%Add_var(mem, 'leaf2sapwood_mass_ratio', mem%leaf2sapwood_mass_ratio, &
            hgrid, vgrid_surface, &
            t_cf('leaf2sapwood_mass_ratio', '', 'leaf to shoot mass ratio'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            loutput=.TRUE., lrestart=.TRUE., &
            initval_r=0.0_wp) ! initvalue
  
      CALL mem%Add_var(mem, 'leaf2root_mass_ratio', mem%leaf2root_mass_ratio, &
            hgrid, vgrid_surface, &
            t_cf('leaf2root_mass_ratio', '', 'leaf to root mass ratio'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            loutput=.TRUE., lrestart=.TRUE., &
            initval_r=0.0_wp) ! initvalue
  
      CALL mem%Add_var(mem, 'fleaf_sunlit_daytime_cl', mem%fleaf_sunlit_daytime_cl, &
            hgrid, vgrid_canopy, &
            t_cf('fleaf_sunlit_daytime_cl', '', 'average of the previous day'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            loutput=.TRUE., lrestart=.TRUE., &
            initval_r=0.0_wp) ! initvalue
  
      CALL mem%Add_var(mem, 't_air_daytime', mem%t_air_daytime, &
            hgrid, vgrid_surface, &
            t_cf('t_air_daytime', '', 'average of the previous day'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            loutput=.TRUE., lrestart=.TRUE., &
            initval_r=0.0_wp) ! initvalue
  
      CALL mem%Add_var(mem, 'press_srf_daytime', mem%press_srf_daytime, &
            hgrid, vgrid_surface, &
            t_cf('press_srf_daytime', '', 'average of the previous day'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            loutput=.TRUE., lrestart=.TRUE., &
            initval_r=0.0_wp) ! initvalue
  
      CALL mem%Add_var(mem, 'co2_mixing_ratio_daytime', mem%co2_mixing_ratio_daytime, &
            hgrid, vgrid_surface, &
            t_cf('co2_mixing_ratio_daytime', '', 'average of the previous day'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            loutput=.TRUE., lrestart=.TRUE., &
            initval_r=0.0_wp) ! initvalue
  
      CALL mem%Add_var(mem, 'ga_daytime', mem%ga_daytime, &
            hgrid, vgrid_surface, &
            t_cf('ga_daytime', '', 'average of the previous day'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            loutput=.TRUE., lrestart=.TRUE., &
            initval_r=0.0_wp) ! initvalue
  
      CALL mem%Add_var(mem, 'beta_sinklim_ps_daytime', mem%beta_sinklim_ps_daytime, &
            hgrid, vgrid_surface, &
            t_cf('beta_sinklim_ps_daytime', 'unitless', 'average of the previous day'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            loutput=.TRUE., lrestart=.TRUE., &
            initval_r=0.0_wp) ! initvalue
  
      CALL mem%Add_var(mem, 't_jmax_opt_daytime', mem%t_jmax_opt_daytime, &
            hgrid, vgrid_surface, &
            t_cf('t_jmax_opt_daytime', '', 'average of the previous day'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            loutput=.TRUE., lrestart=.TRUE., &
            initval_r=0.0_wp) ! initvalue
  
      CALL mem%Add_var(mem, 'fleaf_sunlit_daytime_dacc_cl', mem%fleaf_sunlit_daytime_dacc_cl, &
            hgrid, vgrid_canopy, &
            t_cf('fleaf_sunlit_daytime_dacc_cl', '', 'value accumulated over current day'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            loutput=.TRUE., lrestart=.TRUE., &
            initval_r=0.0_wp) ! initvalue
  
      CALL mem%Add_var(mem, 't_air_daytime_dacc', mem%t_air_daytime_dacc, &
            hgrid, vgrid_surface, &
            t_cf('t_air_daytime_dacc', '', 'average of the previous day'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            loutput=.TRUE., lrestart=.TRUE., &
            initval_r=0.0_wp) ! initvalue
  
      CALL mem%Add_var(mem, 'press_srf_daytime_dacc', mem%press_srf_daytime_dacc, &
            hgrid, vgrid_surface, &
            t_cf('press_srf_daytime_dacc', '', 'average of the previous day'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            loutput=.TRUE., lrestart=.TRUE., &
            initval_r=0.0_wp) ! initvalue
  
      CALL mem%Add_var(mem, 'co2_mixing_ratio_daytime_dacc', mem%co2_mixing_ratio_daytime_dacc, &
            hgrid, vgrid_surface, &
            t_cf('co2_mixing_ratio_daytime_dacc', '', 'average of the previous day'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            loutput=.TRUE., lrestart=.TRUE., &
            initval_r=0.0_wp) ! initvalue
  
      CALL mem%Add_var(mem, 'ga_daytime_dacc', mem%ga_daytime_dacc, &
            hgrid, vgrid_surface, &
            t_cf('ga_daytime_dacc', '', 'average of the previous day'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            loutput=.TRUE., lrestart=.TRUE., &
            initval_r=0.0_wp) ! initvalue
  
      CALL mem%Add_var(mem, 'beta_sinklim_ps_daytime_dacc', mem%beta_sinklim_ps_daytime_dacc, &
            hgrid, vgrid_surface, &
            t_cf('beta_sinklim_ps_daytime_dacc', '', 'average of the previous day'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            loutput=.TRUE., lrestart=.TRUE., &
            initval_r=0.0_wp) ! initvalue
  
      CALL mem%Add_var(mem, 't_jmax_opt_daytime_dacc', mem%t_jmax_opt_daytime_dacc, &
            hgrid, vgrid_surface, &
            t_cf('t_jmax_opt_daytime_dacc', '', 'value accumulated over current day '), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            loutput=.TRUE., lrestart=.TRUE., &
            initval_r=0.0_wp) ! initvalue
  
      CALL mem%Add_var(mem, 'daytime_counter', mem%daytime_counter, &
            hgrid, vgrid_surface, &
            t_cf('daytime_counter', '', 'number of timesteps of current day with light available'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            loutput=.TRUE., lrestart=.TRUE., &
            initval_r=0.0_wp) ! initvalue
  
      CALL mem%Add_var(mem, 'gpp_tlabile_mavg', mem%gpp_tlabile_mavg, &
            hgrid, vgrid_surface, &
            t_cf('gpp_tlabile_mavg', 'µmol m-2 s-1', 'gpp '), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            loutput=.TRUE., lrestart=.TRUE., &
            initval_r=0.0_wp) ! initvalue
  
      CALL mem%Add_var(mem, 'maint_respiration_tlabile_mavg', mem%maint_respiration_tlabile_mavg, &
            hgrid, vgrid_surface, &
            t_cf('maint_respiration_tlabile_mavg', 'µmol m-2 s-1', 'maintenance respiration'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            loutput=.TRUE., lrestart=.TRUE., &
            initval_r=0.0_wp) ! initvalue
  
      CALL mem%Add_var(mem, 'growth_req_n_tlabile_mavg', mem%growth_req_n_tlabile_mavg, &
            hgrid, vgrid_surface, &
            t_cf('growth_req_n_tlabile_mavg', '', 'moles N required for a unit of C growth under current allocation fractions'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            loutput=.TRUE., lrestart=.TRUE., &
            initval_r=0.0_wp) ! initvalue
  
      CALL mem%Add_var(mem, 'growth_req_p_tlabile_mavg', mem%growth_req_p_tlabile_mavg, &
            hgrid, vgrid_surface, &
            t_cf('growth_req_p_tlabile_mavg', '', 'moles P required for a unit of N growth under current allocation fractions'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            loutput=.TRUE., lrestart=.TRUE., &
            initval_r=0.0_wp) ! initvalue
  
      CALL mem%Add_var(mem, 't_air_tphen_mavg', mem%t_air_tphen_mavg, &
            hgrid, vgrid_surface, &
            t_cf('t_air_tphen_mavg', '', 'air temp mavg for phenology calc'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            loutput=.TRUE., lrestart=.TRUE., &
            initval_r=0.0_wp) ! initvalue
  
      CALL mem%Add_var(mem, 'npp_tuptake_mavg', mem%npp_tuptake_mavg, &
            hgrid, vgrid_surface, &
            t_cf('npp_tuptake_mavg', 'µmol m-2 s-1', 'npp of the last week'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            loutput=.TRUE., lrestart=.TRUE., &
            initval_r=0.0_wp) ! initvalue
  
      CALL mem%Add_var(mem, 'demand_uptake_n_tuptake_mavg', mem%demand_uptake_n_tuptake_mavg, &
            hgrid, vgrid_surface, &
            t_cf('demand_uptake_n_tuptake_mavg', '', 'scalar of Vmax for root N uptake '), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            loutput=.TRUE., lrestart=.TRUE., &
            initval_r=0.0_wp) ! initvalue
  
      CALL mem%Add_var(mem, 'demand_uptake_p_tuptake_mavg', mem%demand_uptake_p_tuptake_mavg, &
            hgrid, vgrid_surface, &
            t_cf('demand_uptake_p_tuptake_mavg', '', 'scalar of Vmax for root P uptake'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            loutput=.TRUE., lrestart=.TRUE., &
            initval_r=0.0_wp) ! initvalue
  
      CALL mem%Add_var(mem, 'growth_req_n_tuptake_mavg', mem%growth_req_n_tuptake_mavg, &
            hgrid, vgrid_surface, &
            t_cf('growth_req_n_tuptake_mavg', '', 'moles N required for a unit of C growth under current allocation fractions'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            loutput=.TRUE., lrestart=.TRUE., &
            initval_r=0.0_wp) ! initvalue
  
      CALL mem%Add_var(mem, 'growth_req_p_tuptake_mavg', mem%growth_req_p_tuptake_mavg, &
            hgrid, vgrid_surface, &
            t_cf('growth_req_p_tuptake_mavg', '', 'moles P required for a unit of N growth under current allocation fractions'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            loutput=.TRUE., lrestart=.TRUE., &
            initval_r=0.0_wp) ! initvalue
  
      CALL mem%Add_var(mem, 'cost_n_fixation_talloc_mavg', mem%cost_n_fixation_talloc_mavg, &
            hgrid, vgrid_surface, &
            t_cf('cost_n_fixation_talloc_mavg', 'mol C mol-1 N', 'moving average of N fixation cost'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            loutput=.TRUE., lrestart=.TRUE., &
            initval_r=0.0_wp) ! initvalue

      CALL mem%Add_var(mem, 'exudation_c_tmyc_mavg', mem%exudation_c_tmyc_mavg, &
            hgrid, vgrid_surface, &
            t_cf('exudation_c_tmyc_mavg', 'mol C', 'moving average of C allocation to mycorrhizae'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            loutput=.TRUE., lrestart=.TRUE., &
            initval_r=0.0_wp) ! initvalue

      CALL mem%Add_var(mem, 'fleaf_sunlit_tfrac_mavg_cl', mem%fleaf_sunlit_tfrac_mavg_cl, &
            hgrid, vgrid_canopy, &
            t_cf('fleaf_sunlit_tfrac_mavg_cl', '', 'fraction of sunlit leaves'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            loutput=.TRUE., lrestart=.TRUE., &
            initval_r=0.0_wp) ! initvalue
  
      CALL mem%Add_var(mem, 't_air_tfrac_mavg', mem%t_air_tfrac_mavg, &
            hgrid, vgrid_surface, &
            t_cf('t_air_tfrac_mavg', '', 'air temp'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            loutput=.TRUE., lrestart=.TRUE., &
            initval_r=0.0_wp) ! initvalue
  
      CALL mem%Add_var(mem, 'press_srf_tfrac_mavg', mem%press_srf_tfrac_mavg, &
            hgrid, vgrid_surface, &
            t_cf('press_srf_tfrac_mavg', '', 'air pressure at surface'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            loutput=.TRUE., lrestart=.TRUE., &
            initval_r=0.0_wp) ! initvalue
  
      CALL mem%Add_var(mem, 'co2_mixing_ratio_tfrac_mavg', mem%co2_mixing_ratio_tfrac_mavg, &
            hgrid, vgrid_surface, &
            t_cf('co2_mixing_ratio_tfrac_mavg', '', 'co2 mixing ratio of the air'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            loutput=.TRUE., lrestart=.TRUE., &
            initval_r=0.0_wp) ! initvalue
  
      CALL mem%Add_var(mem, 'ga_tfrac_mavg', mem%ga_tfrac_mavg, &
            hgrid, vgrid_surface, &
            t_cf('ga_tfrac_mavg', '', 'aerodynamic conductance'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            loutput=.TRUE., lrestart=.TRUE., &
            initval_r=0.0_wp) ! initvalue
  
      CALL mem%Add_var(mem, 'beta_sinklim_ps_tfrac_mavg', mem%beta_sinklim_ps_tfrac_mavg, &
            hgrid, vgrid_surface, &
            t_cf('beta_sinklim_ps_tfrac_mavg', '', 'scaling factor to account for soil moisture constraints on photosynthesis '), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            loutput=.TRUE., lrestart=.TRUE., &
            initval_r=0.0_wp) ! initvalue
  
      CALL mem%Add_var(mem, 'fleaf_sunlit_tcnl_mavg_cl', mem%fleaf_sunlit_tcnl_mavg_cl, &
            hgrid, vgrid_canopy, &
            t_cf('fleaf_sunlit_tcnl_mavg_cl', '', 'fraction of sunlit leaves'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            loutput=.TRUE., lrestart=.TRUE., &
            initval_r=0.0_wp) ! initvalue
  
      CALL mem%Add_var(mem, 't_air_tcnl_mavg', mem%t_air_tcnl_mavg, &
            hgrid, vgrid_surface, &
            t_cf('t_air_tcnl_mavg', '', 'air temp'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            loutput=.TRUE., lrestart=.TRUE., &
            initval_r=0.0_wp) ! initvalue
  
      CALL mem%Add_var(mem, 'press_srf_tcnl_mavg', mem%press_srf_tcnl_mavg, &
            hgrid, vgrid_surface, &
            t_cf('press_srf_tcnl_mavg', '', 'air pressure at surface'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            loutput=.TRUE., lrestart=.TRUE., &
            initval_r=0.0_wp) ! initvalue
  
      CALL mem%Add_var(mem, 'co2_mixing_ratio_tcnl_mavg', mem%co2_mixing_ratio_tcnl_mavg, &
            hgrid, vgrid_surface, &
            t_cf('co2_mixing_ratio_tcnl_mavg', '', 'co2 mixing ratio of the air'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            loutput=.TRUE., lrestart=.TRUE., &
            initval_r=0.0_wp) ! initvalue
  
      CALL mem%Add_var(mem, 'ga_tcnl_mavg', mem%ga_tcnl_mavg, &
            hgrid, vgrid_surface, &
            t_cf('ga_tcnl_mavg', '', 'aerodynamic conductance'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            loutput=.TRUE., lrestart=.TRUE., &
            initval_r=0.0_wp) ! initvalue
  
      CALL mem%Add_var(mem, 'beta_sinklim_ps_tcnl_mavg', mem%beta_sinklim_ps_tcnl_mavg, &
            hgrid, vgrid_surface, &
            t_cf('beta_sinklim_ps_tcnl_mavg', '', 'scaling factor to account for soil moisture constraints on conductance'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            loutput=.TRUE., lrestart=.TRUE., &
            initval_r=0.0_wp) ! initvalue

      CALL mem%Add_var(mem, 'uptake_n_tcnl_mavg', mem%uptake_n_tcnl_mavg, &
            hgrid, vgrid_surface, &
            t_cf('uptake_n_tcnl_mavg', '', 'root N uptake'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            loutput=.TRUE., lrestart=.TRUE., &
            initval_r=0.0_wp) ! initvalue

      CALL mem%Add_var(mem, 'growth_cn_tcnl_mavg', mem%growth_cn_tcnl_mavg, &
            hgrid, vgrid_surface, &
            t_cf('growth_cn_tcnl_mavg', '', 'new growth CN ratio'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            loutput=.TRUE., lrestart=.TRUE., &
            initval_r=0.0_wp) ! initvalue

      CALL mem%Add_var(mem, 'growth_np_tcnl_mavg', mem%growth_np_tcnl_mavg, &
            hgrid, vgrid_surface, &
            t_cf('growth_np_tcnl_mavg', '', ' '), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            loutput=.TRUE., lrestart=.TRUE., &
            initval_r=0.0_wp) ! initvalue

      CALL mem%Add_var(mem, 'npp_tcnl_mavg', mem%npp_tcnl_mavg, &
            hgrid, vgrid_surface, &
            t_cf('npp_tcnl_mavg', '', 'NPP'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            loutput=.TRUE., lrestart=.TRUE., &
            initval_r=0.0_wp) ! initvalue

      CALL mem%Add_var(mem, 'fmaint_rate_tcnl_mavg', mem%fmaint_rate_tcnl_mavg, &
            hgrid, vgrid_surface, &
            t_cf('fmaint_rate_tcnl_mavg', '~µmol C mol-1 N s-1', 'N-specific maintentance respiration for optmal leaf N'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            loutput=.TRUE., lrestart=.TRUE., &
            initval_r=0.0_wp) ! initvalue

      CALL mem%Add_var(mem, 'labile_nitrogen_tcnl_mavg', mem%labile_nitrogen_tcnl_mavg, &
            hgrid, vgrid_surface, &
            t_cf('labile_nitrogen_tcnl_mavg', '', ' '), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            loutput=.TRUE., lrestart=.TRUE., &
            initval_r=0.0_wp) ! initvalue

      CALL mem%Add_var(mem, 'labile_carbon_tcnl_mavg', mem%labile_carbon_tcnl_mavg, &
            hgrid, vgrid_surface, &
            t_cf('labile_carbon_tcnl_mavg', '', ' '), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            loutput=.TRUE., lrestart=.TRUE., &
            initval_r=0.0_wp) ! initvalue

      CALL mem%Add_var(mem, 'labile_phosphorus_tcnl_mavg', mem%labile_phosphorus_tcnl_mavg, &
            hgrid, vgrid_surface, &
            t_cf('labile_phosphorus_tcnl_mavg', '', ' '), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            loutput=.TRUE., lrestart=.TRUE., &
            initval_r=0.0_wp) ! initvalue
  
  
      CALL mem%Add_var(mem, 'npp_talloc_mavg', mem%npp_talloc_mavg, &
            hgrid, vgrid_surface, &
            t_cf('npp_talloc_mavg', 'µmol m-2 s-1', 'NPP'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            loutput=.TRUE., lrestart=.TRUE., &
            initval_r=0.0_wp) ! initvalue

      CALL mem%Add_var(mem, 'n_fixation_talloc_mavg', mem%n_fixation_talloc_mavg, &
            hgrid, vgrid_surface, &
            t_cf('n_fixation_talloc_mavg', '', 'N fixation averaged at the allocation timescale'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            loutput=.TRUE., lrestart=.TRUE., &
            initval_r=0.0_wp) ! initvalue
  
      CALL mem%Add_var(mem, 'transpiration_talloc_mavg', mem%transpiration_talloc_mavg, &
            hgrid, vgrid_surface, &
            t_cf('transpiration_talloc_mavg', '', 'average transpiration at allocation time-scale'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            loutput=.TRUE., lrestart=.TRUE., &
            initval_r=0.0_wp) ! initvalue
  
      CALL mem%Add_var(mem, 'unit_npp_talloc_mavg', mem%unit_npp_talloc_mavg, &
            hgrid, vgrid_surface, &
            t_cf('unit_npp_talloc_mavg', 'µmol C mol-1', 'rate of change in NPP with  increment in leaf mass'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            loutput=.TRUE., lrestart=.TRUE., &
            initval_r=0.0_wp) ! initvalue
  
      CALL mem%Add_var(mem, 'unit_transpiration_talloc_mavg', mem%unit_transpiration_talloc_mavg, &
            hgrid, vgrid_surface, &
            t_cf('unit_transpiration_talloc_mavg', '', 'rate of tranpiration change per unit LAI change?'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            loutput=.TRUE., lrestart=.TRUE., &
            initval_r=0.0_wp) ! initvalue
  
      CALL mem%Add_var(mem, 'unit_uptake_n_talloc_mavg', mem%unit_uptake_n_talloc_mavg, &
            hgrid, vgrid_surface, &
            t_cf('unit_uptake_n_talloc_mavg', 'µmol N mol-1', 'actual rate of N uptkae per unit root'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            loutput=.TRUE., lrestart=.TRUE., &
            initval_r=0.0_wp) ! initvalue
  
      CALL mem%Add_var(mem, 'unit_uptake_p_talloc_mavg', mem%unit_uptake_p_talloc_mavg, &
            hgrid, vgrid_surface, &
            t_cf('unit_uptake_p_talloc_mavg', 'µmol P mol-1', 'actual rate of P uptake per unit root'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            loutput=.TRUE., lrestart=.TRUE., &
            initval_r=0.0_wp) ! initvalue
  
      CALL mem%Add_var(mem, 'dphi_talloc_mavg', mem%dphi_talloc_mavg, &
            hgrid, vgrid_surface, &
            t_cf('dphi_talloc_mavg', 'MPa', 'leaf-soil water potential gradient'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            loutput=.TRUE., lrestart=.TRUE., &
            initval_r=0.0_wp) ! initvalue
  
      CALL mem%Add_var(mem, 'growth_cn_talloc_mavg', mem%growth_cn_talloc_mavg, &
            hgrid, vgrid_surface, &
            t_cf('growth_cn_talloc_mavg', '', 'plant growth cn ratio, equivalent to 1/growth_req_n'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            loutput=.TRUE., lrestart=.TRUE., &
            initval_r=0.0_wp) ! initvalue
  
      CALL mem%Add_var(mem, 'growth_cp_talloc_mavg', mem%growth_cp_talloc_mavg, &
            hgrid, vgrid_surface, &
            t_cf('growth_cp_talloc_mavg', '', 'plant growth cp ratio, equivalent to 1/(growthreq_n*growth_req_p)'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            loutput=.TRUE., lrestart=.TRUE., &
            initval_r=0.0_wp) ! initvalue

      CALL mem%Add_var(mem, 'growth_np_talloc_mavg', mem%growth_np_talloc_mavg, &
            hgrid, vgrid_surface, &
            t_cf('growth_np_talloc_mavg', '', 'plant growth hp ratio, equivalent to 1/(growth_req_p)'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            loutput=.TRUE., lrestart=.TRUE., &
            initval_r=0.0_wp) ! initvalue
  
      CALL mem%Add_var(mem, 'labile_carbon_talloc_mavg', mem%labile_carbon_talloc_mavg, &
            hgrid, vgrid_surface, &
            t_cf('labile_carbon_talloc_mavg', '', 'moving average of labile carbon pool (allometry calculations)'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            loutput=.TRUE., lrestart=.TRUE., &
            initval_r=0.0_wp) ! initvalue
  
      CALL mem%Add_var(mem, 'labile_nitrogen_talloc_mavg', mem%labile_nitrogen_talloc_mavg, &
            hgrid, vgrid_surface, &
            t_cf('labile_nitrogen_talloc_mavg', '', 'moving average of labile nitrogen pool (allometry calculations)'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            loutput=.TRUE., lrestart=.TRUE., &
            initval_r=0.0_wp) ! initvalue

      CALL mem%Add_var(mem, 'labile_phosphorus_talloc_mavg', mem%labile_phosphorus_talloc_mavg, &
            hgrid, vgrid_surface, &
            t_cf('labile_phosphorus_talloc_mavg', '', 'moving average of labile phosphorus pool (allometry calculations)'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            loutput=.TRUE., lrestart=.TRUE., &
            initval_r=0.0_wp) ! initvalue

      CALL mem%Add_var(mem, 'reserve_carbon_talloc_mavg', mem%reserve_carbon_talloc_mavg, &
            hgrid, vgrid_surface, &
            t_cf('reserve_carbon_talloc_mavg', '', 'moving average of reserve carbon (optimal allocation)'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            loutput=.TRUE., lrestart=.TRUE., &
            initval_r=0.0_wp) ! initvalue

      CALL mem%Add_var(mem, 'reserve_nitrogen_talloc_mavg', mem%reserve_nitrogen_talloc_mavg, &
            hgrid, vgrid_surface, &
            t_cf('reserve_nitrogen_talloc_mavg', '', 'moving average of reserve nitrogen (optimal allocation)'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            loutput=.TRUE., lrestart=.TRUE., &
            initval_r=0.0_wp) ! initvalue

      CALL mem%Add_var(mem, 'reserve_phosphorus_talloc_mavg', mem%reserve_phosphorus_talloc_mavg, &
            hgrid, vgrid_surface, &
            t_cf('reserve_phosphorus_talloc_mavg', '', 'moving average of reserve phosphorus (optimal allocation)'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            loutput=.TRUE., lrestart=.TRUE., &
            initval_r=0.0_wp) ! initvalue
  
      CALL mem%Add_var(mem, 'w_root_lim_talloc_mavg', mem%w_root_lim_talloc_mavg, &
            hgrid, vgrid_surface, &
            t_cf('w_root_lim_talloc_mavg', '', 'moving average of ...'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            loutput=.TRUE., lrestart=.TRUE., &
            initval_r=1.0_wp) ! initvalue
      
      CALL mem%Add_var(mem, 'leaf2root_troot_mavg', mem%leaf2root_troot_mavg, &
            hgrid, vgrid_surface, &
            t_cf('leaf2root_troot_mavg', '', 'leaf to root mass ratio'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            loutput=.TRUE., lrestart=.TRUE., &
            initval_r=0.0_wp) ! initvalue
  
      CALL mem%Add_var(mem, 'unit_npp_troot_mavg', mem%unit_npp_troot_mavg, &
            hgrid, vgrid_surface, &
            t_cf('unit_npp_troot_mavg', 'µmol C mol-1', 'rate of change in NPP with increment in leaf mass'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            loutput=.TRUE., lrestart=.TRUE., &
            initval_r=0.0_wp) ! initvalue
  
      CALL mem%Add_var(mem, 'unit_uptake_n_pot_troot_mavg', mem%unit_uptake_n_pot_troot_mavg, &
            hgrid, vgrid_surface, &
            t_cf('unit_uptake_n_pot_troot_mavg', 'µmol N mol-1', 'rate of N uptkae per unit root'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            loutput=.TRUE., lrestart=.TRUE., &
            initval_r=0.0_wp) ! initvalue
  
      CALL mem%Add_var(mem, 'unit_uptake_p_pot_troot_mavg', mem%unit_uptake_p_pot_troot_mavg, &
            hgrid, vgrid_surface, &
            t_cf('unit_uptake_p_pot_troot_mavg', 'µmol P mol-1', 'rate of P uptake per unit root'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            loutput=.TRUE., lrestart=.TRUE., &
            initval_r=0.0_wp) ! initvalue
  
      CALL mem%Add_var(mem, 'fmaint_rate_troot_mavg', mem%fmaint_rate_troot_mavg, &
            hgrid, vgrid_surface, &
            t_cf('fmaint_rate_troot_mavg', '~µmol C mol-1 N s-1', 'moving average of N-specific maintentance respiration'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            loutput=.TRUE., lrestart=.TRUE., &
            initval_r=0.0_wp) ! initvalue
  
      CALL mem%Add_var(mem, 'an_boc_tvegdyn_mavg', mem%an_boc_tvegdyn_mavg, &
            hgrid, vgrid_surface, &
            t_cf('an_boc_tvegdyn_mavg', '', 'net photosynthesis at bottom of canopy, lowest canopy layer existing'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            loutput=.TRUE., lrestart=.TRUE., &
            initval_r=0.0_wp) ! initvalue
  
      CALL mem%Add_var(mem, 'net_growth_tvegdyn_mavg', mem%net_growth_tvegdyn_mavg, &
            hgrid, vgrid_surface, &
            t_cf('net_growth_tvegdyn_mavg', 'µmol C m-2 s-1', 'net growth, npp minus turnover'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            loutput=.TRUE., lrestart=.TRUE., &
            initval_r=0.0_wp) ! initvalue
  
      CALL mem%Add_var(mem, 'lai_tvegdyn_mavg', mem%lai_tvegdyn_mavg, &
            hgrid, vgrid_surface, &
            t_cf('lai_tvegdyn_mavg', 'm2 m-2', 'lai at vegdynamics timescale'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            loutput=.TRUE., lrestart=.TRUE., &
            initval_r=0.0_wp) ! initvalue
  
      CALL mem%Add_var(mem, 't_air_tacclim_mavg', mem%t_air_tacclim_mavg, &
            hgrid, vgrid_surface, &
            t_cf('t_air_tacclim_mavg', '', 'average air temperature for respiration acclimation'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            loutput=.TRUE., lrestart=.TRUE., &
            initval_r=0.0_wp) ! initvalue
  
      CALL mem%Add_var(mem, 't_soil_root_tacclim_mavg', mem%t_soil_root_tacclim_mavg, &
            hgrid, vgrid_surface, &
            t_cf('t_soil_root_tacclim_mavg', '', 'average root temperature for respiration acclimation'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            loutput=.TRUE., lrestart=.TRUE., &
            initval_r=0.0_wp) ! initvalue
  
      CALL mem%Add_var(mem, 't_air_week_mavg', mem%t_air_week_mavg, &
            hgrid, vgrid_surface, &
            t_cf('t_air_week_mavg', '', 'avrg air temperature across 7 days'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            loutput=.TRUE., lrestart=.TRUE., &
            initval_r=0.0_wp) ! initvalue
  
      CALL mem%Add_var(mem, 't_air_month_mavg', mem%t_air_month_mavg, &
            hgrid, vgrid_surface, &
            t_cf('t_air_month_mavg', '', 'avrg air temperature across 30 days'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            loutput=.TRUE., lrestart=.TRUE., &
            initval_r=0.0_wp) ! initvalue
  
      CALL mem%Add_var(mem, 't_jmax_opt_mavg', mem%t_jmax_opt_mavg, &
            hgrid, vgrid_surface, &
            t_cf('t_jmax_opt_mavg', '', 'optimum temperature for electron transport of photosynthesis'), &
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

    ! not allocated: this%path
    ! no need to dealloc the variables 'this%ELEM%VAR' because these are dealloc by the FINAL routine of var_real2d & var_real3d
    IF (ASSOCIATED(this%carbon))        DEALLOCATE(this%carbon)
    IF (ASSOCIATED(this%nitrogen))      DEALLOCATE(this%nitrogen)
    IF (ASSOCIATED(this%phosphorus))    DEALLOCATE(this%phosphorus)
    IF (ASSOCIATED(this%carbon13))      DEALLOCATE(this%carbon13)
    IF (ASSOCIATED(this%carbon14))      DEALLOCATE(this%carbon14)
    IF (ASSOCIATED(this%nitrogen15))    DEALLOCATE(this%nitrogen15)
    IF (ALLOCATED(this%shortname))      DEALLOCATE(this%shortname)
    IF (ALLOCATED(this%name))           DEALLOCATE(this%name)
    IF (ALLOCATED(this%element_unit))   DEALLOCATE(this%element_unit)
    IF (ALLOCATED(this%element_list))   DEALLOCATE(this%element_list)
    IF (ALLOCATED(this%var_list))       DEALLOCATE(this%var_list)
    IF (ASSOCIATED(this%pool_list))     DEALLOCATE(this%pool_list)
    IF (ASSOCIATED(this%leaf_list))     DEALLOCATE(this%leaf_list)
   
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

    CHARACTER(len=*), PARAMETER :: routine = modname//':Assign_scalar_to_veg_elem_p'

    lhs%carbon     = rhs
    lhs%nitrogen   = rhs
    lhs%phosphorus = rhs
    lhs%carbon13   = rhs
    lhs%carbon14   = rhs
    lhs%nitrogen15 = rhs

  END SUBROUTINE Assign_scalar_to_veg_elem_p
  !--------------------------------------------------------------------------------------------------------------------------------
  ! t_veg_elements_pool = field1d(:) | ON_CHUNK
  !--------------------------------------------------------------------------------------------------------------------------------
  SUBROUTINE Assign_field1d_to_veg_elem_p(lhs, rhs)
  
    CLASS(t_veg_elements_pool), INTENT(inout) :: lhs
    REAL(wp),                   INTENT(in)    :: rhs(:)

    CHARACTER(len=*), PARAMETER :: routine = modname//':Assign_field1d_to_veg_elem_p'

    lhs%carbon     = rhs(:)
    lhs%nitrogen   = rhs(:)
    lhs%phosphorus = rhs(:)
    lhs%carbon13   = rhs(:)
    lhs%carbon14   = rhs(:)
    lhs%nitrogen15 = rhs(:)
  
  END SUBROUTINE Assign_field1d_to_veg_elem_p
  !--------------------------------------------------------------------------------------------------------------------------------
  ! t_veg_elements_pool = field2d(:,:) | ON_CHUNK & ON_DOMAIN | field2d(:,1) & field2d(:,:)
  !--------------------------------------------------------------------------------------------------------------------------------
  SUBROUTINE Assign_field2d_to_veg_elem_p(lhs, rhs)
  
    CLASS(t_veg_elements_pool), INTENT(inout) :: lhs
    REAL(wp),                   INTENT(in)    :: rhs(:,:)

    CHARACTER(len=*), PARAMETER :: routine = modname//':Assign_field2d_to_veg_elem_p'

    lhs%carbon     = rhs(:,:)
    lhs%nitrogen   = rhs(:,:)
    lhs%phosphorus = rhs(:,:)
    lhs%carbon13   = rhs(:,:)
    lhs%carbon14   = rhs(:,:)
    lhs%nitrogen15 = rhs(:,:)
  
  END SUBROUTINE Assign_field2d_to_veg_elem_p
  !--------------------------------------------------------------------------------------------------------------------------------
  ! t_veg_elements_pool = t_veg_elements_pool | ON_CHUNK & ON_DOMAIN
  !--------------------------------------------------------------------------------------------------------------------------------
  SUBROUTINE Assign_veg_elem_p_to_veg_elem_p(lhs, rhs)

    CLASS(t_veg_elements_pool), INTENT(inout) :: lhs
    CLASS(t_veg_elements_pool), INTENT(in)    :: rhs

    CHARACTER(len=*), PARAMETER :: routine = modname//':Assign_veg_elem_p_to_veg_elem_p'

    lhs%carbon     = rhs%carbon    
    lhs%nitrogen   = rhs%nitrogen  
    lhs%phosphorus = rhs%phosphorus
    lhs%carbon13   = rhs%carbon13  
    lhs%carbon14   = rhs%carbon14  
    lhs%nitrogen15 = rhs%nitrogen15

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
    TYPE(t_veg_elements_pool),      TARGET      :: structOut
  
    CHARACTER(len=*), PARAMETER :: routine = modname//':Add_veg_elem_p_and_scalar'
  
    ! init temporary pool structOut
    CALL structOut%Init(name='op_ov_help', shortname='op_ov_help', l_elements=.TRUE., element_unit=TRIM(struct1%element_unit))
    structOut%carbon%owner_model_id      = struct1%carbon%owner_model_id
    structOut%carbon%owner_proc_id       = struct1%carbon%owner_proc_id
    structOut%carbon%owner_tile_path     = struct1%carbon%owner_tile_path
    ALLOCATE(structOut%carbon%ptr2d,     source=struct1%carbon%ptr2d)
    structOut%nitrogen%owner_model_id    = struct1%nitrogen%owner_model_id
    structOut%nitrogen%owner_proc_id     = struct1%nitrogen%owner_proc_id
    structOut%nitrogen%owner_tile_path   = struct1%nitrogen%owner_tile_path
    ALLOCATE(structOut%nitrogen%ptr2d,   source=struct1%nitrogen%ptr2d)
    structOut%phosphorus%owner_model_id  = struct1%phosphorus%owner_model_id
    structOut%phosphorus%owner_proc_id   = struct1%phosphorus%owner_proc_id
    structOut%phosphorus%owner_tile_path = struct1%phosphorus%owner_tile_path
    ALLOCATE(structOut%phosphorus%ptr2d, source=struct1%phosphorus%ptr2d)
    structOut%carbon13%owner_model_id    = struct1%carbon13%owner_model_id
    structOut%carbon13%owner_proc_id     = struct1%carbon13%owner_proc_id
    structOut%carbon13%owner_tile_path   = struct1%carbon13%owner_tile_path
    ALLOCATE(structOut%carbon13%ptr2d,   source=struct1%carbon13%ptr2d)
    structOut%carbon14%owner_model_id    = struct1%carbon14%owner_model_id
    structOut%carbon14%owner_proc_id     = struct1%carbon14%owner_proc_id
    structOut%carbon14%owner_tile_path   = struct1%carbon14%owner_tile_path
    ALLOCATE(structOut%carbon14%ptr2d,   source=struct1%carbon14%ptr2d)
    structOut%nitrogen15%owner_model_id  = struct1%nitrogen15%owner_model_id
    structOut%nitrogen15%owner_proc_id   = struct1%nitrogen15%owner_proc_id
    structOut%nitrogen15%owner_tile_path = struct1%nitrogen15%owner_tile_path
    ALLOCATE(structOut%nitrogen15%ptr2d, source=struct1%nitrogen15%ptr2d)

    ! run calculation
    ! the differentiation between ON_CHUNK/ON_DOMAIN is done at the level of variables (i.e., the op_ov routines of the jsb_var)
    structOut%carbon          = struct1%carbon      + scalar
    structOut%nitrogen        = struct1%nitrogen    + scalar
    structOut%phosphorus      = struct1%phosphorus  + scalar
    structOut%carbon13        = struct1%carbon13    + scalar
    structOut%carbon14        = struct1%carbon14    + scalar
    structOut%nitrogen15      = struct1%nitrogen15  + scalar

    ! is this needed here too (not only with the jsb_var op_ov)
    structOut%carbon%ptr      => structOut%carbon%ptr2d(:,:)
    structOut%nitrogen%ptr    => structOut%nitrogen%ptr2d(:,:)
    structOut%phosphorus%ptr  => structOut%phosphorus%ptr2d(:,:)
    structOut%carbon13%ptr    => structOut%carbon13%ptr2d(:,:)
    structOut%carbon14%ptr    => structOut%carbon14%ptr2d(:,:)
    structOut%nitrogen15%ptr  => structOut%nitrogen15%ptr2d(:,:)
  
  END FUNCTION Add_veg_elem_p_and_scalar
  !--------------------------------------------------------------------------------------------------------------------------------
  ! scalar + t_veg_elements_pool        | ON_CHUNK & ON_DOMAIN
  !--------------------------------------------------------------------------------------------------------------------------------
  FUNCTION Add_scalar_and_veg_elem_p(scalar, struct1) RESULT(structOut)
  
    REAL(wp),                       INTENT(in)  :: scalar
    CLASS(t_veg_elements_pool),     INTENT(in)  :: struct1
    TYPE(t_veg_elements_pool),      TARGET      :: structOut
  
    CHARACTER(len=*), PARAMETER :: routine = modname//':Add_scalar_and_veg_elem_p'
  
    ! init temporary pool structOut
    CALL structOut%Init(name='op_ov_help', shortname='op_ov_help', l_elements=.TRUE., element_unit=TRIM(struct1%element_unit))
    structOut%carbon%owner_model_id      = struct1%carbon%owner_model_id
    structOut%carbon%owner_proc_id       = struct1%carbon%owner_proc_id
    structOut%carbon%owner_tile_path     = struct1%carbon%owner_tile_path
    ALLOCATE(structOut%carbon%ptr2d,     source=struct1%carbon%ptr2d)
    structOut%nitrogen%owner_model_id    = struct1%nitrogen%owner_model_id
    structOut%nitrogen%owner_proc_id     = struct1%nitrogen%owner_proc_id
    structOut%nitrogen%owner_tile_path   = struct1%nitrogen%owner_tile_path
    ALLOCATE(structOut%nitrogen%ptr2d,   source=struct1%nitrogen%ptr2d)
    structOut%phosphorus%owner_model_id  = struct1%phosphorus%owner_model_id
    structOut%phosphorus%owner_proc_id   = struct1%phosphorus%owner_proc_id
    structOut%phosphorus%owner_tile_path = struct1%phosphorus%owner_tile_path
    ALLOCATE(structOut%phosphorus%ptr2d, source=struct1%phosphorus%ptr2d)
    structOut%carbon13%owner_model_id    = struct1%carbon13%owner_model_id
    structOut%carbon13%owner_proc_id     = struct1%carbon13%owner_proc_id
    structOut%carbon13%owner_tile_path   = struct1%carbon13%owner_tile_path
    ALLOCATE(structOut%carbon13%ptr2d,   source=struct1%carbon13%ptr2d)
    structOut%carbon14%owner_model_id    = struct1%carbon14%owner_model_id
    structOut%carbon14%owner_proc_id     = struct1%carbon14%owner_proc_id
    structOut%carbon14%owner_tile_path   = struct1%carbon14%owner_tile_path
    ALLOCATE(structOut%carbon14%ptr2d,   source=struct1%carbon14%ptr2d)
    structOut%nitrogen15%owner_model_id  = struct1%nitrogen15%owner_model_id
    structOut%nitrogen15%owner_proc_id   = struct1%nitrogen15%owner_proc_id
    structOut%nitrogen15%owner_tile_path = struct1%nitrogen15%owner_tile_path
    ALLOCATE(structOut%nitrogen15%ptr2d, source=struct1%nitrogen15%ptr2d)

    ! run calculation
    ! the differentiation between ON_CHUNK/ON_DOMAIN is done at the level of variables (i.e., the op_ov routines of the jsb_var)
    structOut%carbon          = scalar + struct1%carbon    
    structOut%nitrogen        = scalar + struct1%nitrogen  
    structOut%phosphorus      = scalar + struct1%phosphorus
    structOut%carbon13        = scalar + struct1%carbon13  
    structOut%carbon14        = scalar + struct1%carbon14  
    structOut%nitrogen15      = scalar + struct1%nitrogen15

    ! is this needed here too (not only with the jsb_var op_ov)
    structOut%carbon%ptr      => structOut%carbon%ptr2d(:,:)
    structOut%nitrogen%ptr    => structOut%nitrogen%ptr2d(:,:)
    structOut%phosphorus%ptr  => structOut%phosphorus%ptr2d(:,:)
    structOut%carbon13%ptr    => structOut%carbon13%ptr2d(:,:)
    structOut%carbon14%ptr    => structOut%carbon14%ptr2d(:,:)
    structOut%nitrogen15%ptr  => structOut%nitrogen15%ptr2d(:,:)
  
  END FUNCTION Add_scalar_and_veg_elem_p
  !--------------------------------------------------------------------------------------------------------------------------------
  ! t_veg_elements_pool + field1d(:)        | ON_CHUNK
  !--------------------------------------------------------------------------------------------------------------------------------
  FUNCTION Add_veg_elem_p_and_field1d(struct1, field1d) RESULT(structOut)
  
    CLASS(t_veg_elements_pool),     INTENT(in)  :: struct1
    REAL(wp),                       INTENT(in)  :: field1d(:)
    TYPE(t_veg_elements_pool),      TARGET      :: structOut
  
    CHARACTER(len=*), PARAMETER :: routine = modname//':Add_veg_elem_p_and_field1d'
  
    ! init temporary pool structOut
    CALL structOut%Init(name='op_ov_help', shortname='op_ov_help', l_elements=.TRUE., element_unit=TRIM(struct1%element_unit))
    structOut%carbon%owner_model_id      = struct1%carbon%owner_model_id
    structOut%carbon%owner_proc_id       = struct1%carbon%owner_proc_id
    structOut%carbon%owner_tile_path     = struct1%carbon%owner_tile_path
    ALLOCATE(structOut%carbon%ptr2d,     source=struct1%carbon%ptr2d)
    structOut%nitrogen%owner_model_id    = struct1%nitrogen%owner_model_id
    structOut%nitrogen%owner_proc_id     = struct1%nitrogen%owner_proc_id
    structOut%nitrogen%owner_tile_path   = struct1%nitrogen%owner_tile_path
    ALLOCATE(structOut%nitrogen%ptr2d,   source=struct1%nitrogen%ptr2d)
    structOut%phosphorus%owner_model_id  = struct1%phosphorus%owner_model_id
    structOut%phosphorus%owner_proc_id   = struct1%phosphorus%owner_proc_id
    structOut%phosphorus%owner_tile_path = struct1%phosphorus%owner_tile_path
    ALLOCATE(structOut%phosphorus%ptr2d, source=struct1%phosphorus%ptr2d)
    structOut%carbon13%owner_model_id    = struct1%carbon13%owner_model_id
    structOut%carbon13%owner_proc_id     = struct1%carbon13%owner_proc_id
    structOut%carbon13%owner_tile_path   = struct1%carbon13%owner_tile_path
    ALLOCATE(structOut%carbon13%ptr2d,   source=struct1%carbon13%ptr2d)
    structOut%carbon14%owner_model_id    = struct1%carbon14%owner_model_id
    structOut%carbon14%owner_proc_id     = struct1%carbon14%owner_proc_id
    structOut%carbon14%owner_tile_path   = struct1%carbon14%owner_tile_path
    ALLOCATE(structOut%carbon14%ptr2d,   source=struct1%carbon14%ptr2d)
    structOut%nitrogen15%owner_model_id  = struct1%nitrogen15%owner_model_id
    structOut%nitrogen15%owner_proc_id   = struct1%nitrogen15%owner_proc_id
    structOut%nitrogen15%owner_tile_path = struct1%nitrogen15%owner_tile_path
    ALLOCATE(structOut%nitrogen15%ptr2d, source=struct1%nitrogen15%ptr2d)

    ! run calculation
    ! the differentiation between ON_CHUNK/ON_DOMAIN is done at the level of variables (i.e., the op_ov routines of the jsb_var)
    structOut%carbon          = struct1%carbon      + field1d(:)
    structOut%nitrogen        = struct1%nitrogen    + field1d(:)
    structOut%phosphorus      = struct1%phosphorus  + field1d(:)
    structOut%carbon13        = struct1%carbon13    + field1d(:)
    structOut%carbon14        = struct1%carbon14    + field1d(:)
    structOut%nitrogen15      = struct1%nitrogen15  + field1d(:)

    ! is this needed here too (not only with the jsb_var op_ov)
    structOut%carbon%ptr      => structOut%carbon%ptr2d(:,:)
    structOut%nitrogen%ptr    => structOut%nitrogen%ptr2d(:,:)
    structOut%phosphorus%ptr  => structOut%phosphorus%ptr2d(:,:)
    structOut%carbon13%ptr    => structOut%carbon13%ptr2d(:,:)
    structOut%carbon14%ptr    => structOut%carbon14%ptr2d(:,:)
    structOut%nitrogen15%ptr  => structOut%nitrogen15%ptr2d(:,:)
  
  END FUNCTION Add_veg_elem_p_and_field1d
  !--------------------------------------------------------------------------------------------------------------------------------
  ! field1d(:) + t_veg_elements_pool        | ON_CHUNK
  !--------------------------------------------------------------------------------------------------------------------------------
  FUNCTION Add_field1d_and_veg_elem_p(field1d, struct1) RESULT(structOut)
  
    REAL(wp),                       INTENT(in)  :: field1d(:)
    CLASS(t_veg_elements_pool),     INTENT(in)  :: struct1
    TYPE(t_veg_elements_pool),      TARGET      :: structOut
  
    CHARACTER(len=*), PARAMETER :: routine = modname//':Add_field1d_and_veg_elem_p'
  
    ! init temporary pool structOut
    CALL structOut%Init(name='op_ov_help', shortname='op_ov_help', l_elements=.TRUE., element_unit=TRIM(struct1%element_unit))
    structOut%carbon%owner_model_id      = struct1%carbon%owner_model_id
    structOut%carbon%owner_proc_id       = struct1%carbon%owner_proc_id
    structOut%carbon%owner_tile_path     = struct1%carbon%owner_tile_path
    ALLOCATE(structOut%carbon%ptr2d,     source=struct1%carbon%ptr2d)
    structOut%nitrogen%owner_model_id    = struct1%nitrogen%owner_model_id
    structOut%nitrogen%owner_proc_id     = struct1%nitrogen%owner_proc_id
    structOut%nitrogen%owner_tile_path   = struct1%nitrogen%owner_tile_path
    ALLOCATE(structOut%nitrogen%ptr2d,   source=struct1%nitrogen%ptr2d)
    structOut%phosphorus%owner_model_id  = struct1%phosphorus%owner_model_id
    structOut%phosphorus%owner_proc_id   = struct1%phosphorus%owner_proc_id
    structOut%phosphorus%owner_tile_path = struct1%phosphorus%owner_tile_path
    ALLOCATE(structOut%phosphorus%ptr2d, source=struct1%phosphorus%ptr2d)
    structOut%carbon13%owner_model_id    = struct1%carbon13%owner_model_id
    structOut%carbon13%owner_proc_id     = struct1%carbon13%owner_proc_id
    structOut%carbon13%owner_tile_path   = struct1%carbon13%owner_tile_path
    ALLOCATE(structOut%carbon13%ptr2d,   source=struct1%carbon13%ptr2d)
    structOut%carbon14%owner_model_id    = struct1%carbon14%owner_model_id
    structOut%carbon14%owner_proc_id     = struct1%carbon14%owner_proc_id
    structOut%carbon14%owner_tile_path   = struct1%carbon14%owner_tile_path
    ALLOCATE(structOut%carbon14%ptr2d,   source=struct1%carbon14%ptr2d)
    structOut%nitrogen15%owner_model_id  = struct1%nitrogen15%owner_model_id
    structOut%nitrogen15%owner_proc_id   = struct1%nitrogen15%owner_proc_id
    structOut%nitrogen15%owner_tile_path = struct1%nitrogen15%owner_tile_path
    ALLOCATE(structOut%nitrogen15%ptr2d, source=struct1%nitrogen15%ptr2d)

    ! run calculation
    ! the differentiation between ON_CHUNK/ON_DOMAIN is done at the level of variables (i.e., the op_ov routines of the jsb_var)
    structOut%carbon          = field1d(:) + struct1%carbon    
    structOut%nitrogen        = field1d(:) + struct1%nitrogen  
    structOut%phosphorus      = field1d(:) + struct1%phosphorus
    structOut%carbon13        = field1d(:) + struct1%carbon13  
    structOut%carbon14        = field1d(:) + struct1%carbon14  
    structOut%nitrogen15      = field1d(:) + struct1%nitrogen15

    ! is this needed here too (not only with the jsb_var op_ov)
    structOut%carbon%ptr      => structOut%carbon%ptr2d(:,:)
    structOut%nitrogen%ptr    => structOut%nitrogen%ptr2d(:,:)
    structOut%phosphorus%ptr  => structOut%phosphorus%ptr2d(:,:)
    structOut%carbon13%ptr    => structOut%carbon13%ptr2d(:,:)
    structOut%carbon14%ptr    => structOut%carbon14%ptr2d(:,:)
    structOut%nitrogen15%ptr  => structOut%nitrogen15%ptr2d(:,:)
  
  END FUNCTION Add_field1d_and_veg_elem_p
  !--------------------------------------------------------------------------------------------------------------------------------
  ! t_veg_elements_pool + field2d(:,:)        | ON_CHUNK & ON_DOMAIN | field2d(:,1) & field2d(:,:)
  !--------------------------------------------------------------------------------------------------------------------------------
  FUNCTION Add_veg_elem_p_and_field2d(struct1, field2d) RESULT(structOut)
  
    CLASS(t_veg_elements_pool),     INTENT(in)  :: struct1
    REAL(wp),                       INTENT(in)  :: field2d(:,:)
    TYPE(t_veg_elements_pool),      TARGET      :: structOut
  
    CHARACTER(len=*), PARAMETER :: routine = modname//':Add_veg_elem_p_and_field2d'
  
    ! init temporary pool structOut
    CALL structOut%Init(name='op_ov_help', shortname='op_ov_help', l_elements=.TRUE., element_unit=TRIM(struct1%element_unit))
    structOut%carbon%owner_model_id      = struct1%carbon%owner_model_id
    structOut%carbon%owner_proc_id       = struct1%carbon%owner_proc_id
    structOut%carbon%owner_tile_path     = struct1%carbon%owner_tile_path
    ALLOCATE(structOut%carbon%ptr2d,     source=struct1%carbon%ptr2d)
    structOut%nitrogen%owner_model_id    = struct1%nitrogen%owner_model_id
    structOut%nitrogen%owner_proc_id     = struct1%nitrogen%owner_proc_id
    structOut%nitrogen%owner_tile_path   = struct1%nitrogen%owner_tile_path
    ALLOCATE(structOut%nitrogen%ptr2d,   source=struct1%nitrogen%ptr2d)
    structOut%phosphorus%owner_model_id  = struct1%phosphorus%owner_model_id
    structOut%phosphorus%owner_proc_id   = struct1%phosphorus%owner_proc_id
    structOut%phosphorus%owner_tile_path = struct1%phosphorus%owner_tile_path
    ALLOCATE(structOut%phosphorus%ptr2d, source=struct1%phosphorus%ptr2d)
    structOut%carbon13%owner_model_id    = struct1%carbon13%owner_model_id
    structOut%carbon13%owner_proc_id     = struct1%carbon13%owner_proc_id
    structOut%carbon13%owner_tile_path   = struct1%carbon13%owner_tile_path
    ALLOCATE(structOut%carbon13%ptr2d,   source=struct1%carbon13%ptr2d)
    structOut%carbon14%owner_model_id    = struct1%carbon14%owner_model_id
    structOut%carbon14%owner_proc_id     = struct1%carbon14%owner_proc_id
    structOut%carbon14%owner_tile_path   = struct1%carbon14%owner_tile_path
    ALLOCATE(structOut%carbon14%ptr2d,   source=struct1%carbon14%ptr2d)
    structOut%nitrogen15%owner_model_id  = struct1%nitrogen15%owner_model_id
    structOut%nitrogen15%owner_proc_id   = struct1%nitrogen15%owner_proc_id
    structOut%nitrogen15%owner_tile_path = struct1%nitrogen15%owner_tile_path
    ALLOCATE(structOut%nitrogen15%ptr2d, source=struct1%nitrogen15%ptr2d)

    ! run calculation
    ! the differentiation between ON_CHUNK/ON_DOMAIN is done at the level of variables (i.e., the op_ov routines of the jsb_var)
    structOut%carbon          = struct1%carbon      + field2d(:,:)
    structOut%nitrogen        = struct1%nitrogen    + field2d(:,:)
    structOut%phosphorus      = struct1%phosphorus  + field2d(:,:)
    structOut%carbon13        = struct1%carbon13    + field2d(:,:)
    structOut%carbon14        = struct1%carbon14    + field2d(:,:)
    structOut%nitrogen15      = struct1%nitrogen15  + field2d(:,:)

    ! is this needed here too (not only with the jsb_var op_ov)
    structOut%carbon%ptr      => structOut%carbon%ptr2d(:,:)
    structOut%nitrogen%ptr    => structOut%nitrogen%ptr2d(:,:)
    structOut%phosphorus%ptr  => structOut%phosphorus%ptr2d(:,:)
    structOut%carbon13%ptr    => structOut%carbon13%ptr2d(:,:)
    structOut%carbon14%ptr    => structOut%carbon14%ptr2d(:,:)
    structOut%nitrogen15%ptr  => structOut%nitrogen15%ptr2d(:,:)
  
  END FUNCTION Add_veg_elem_p_and_field2d
  !--------------------------------------------------------------------------------------------------------------------------------
  ! t_veg_elements_pool + t_veg_elements_pool | ON_CHUNK & ON_DOMAIN
  !--------------------------------------------------------------------------------------------------------------------------------
  FUNCTION Add_veg_elem_p_and_veg_elem_p(struct1, struct2) RESULT(structOut)
  
    CLASS(t_veg_elements_pool),     INTENT(in)  :: struct1, struct2
    TYPE(t_veg_elements_pool),      TARGET      :: structOut
  
    CHARACTER(len=*), PARAMETER :: routine = modname//':Add_veg_elem_p_and_veg_elem_p'
  
    ! init temporary pool structOut
    CALL structOut%Init(name='op_ov_help', shortname='op_ov_help', l_elements=.TRUE., element_unit=TRIM(struct1%element_unit))
    structOut%carbon%owner_model_id      = struct1%carbon%owner_model_id
    structOut%carbon%owner_proc_id       = struct1%carbon%owner_proc_id
    structOut%carbon%owner_tile_path     = struct1%carbon%owner_tile_path
    ALLOCATE(structOut%carbon%ptr2d,     source=struct1%carbon%ptr2d)
    structOut%nitrogen%owner_model_id    = struct1%nitrogen%owner_model_id
    structOut%nitrogen%owner_proc_id     = struct1%nitrogen%owner_proc_id
    structOut%nitrogen%owner_tile_path   = struct1%nitrogen%owner_tile_path
    ALLOCATE(structOut%nitrogen%ptr2d,   source=struct1%nitrogen%ptr2d)
    structOut%phosphorus%owner_model_id  = struct1%phosphorus%owner_model_id
    structOut%phosphorus%owner_proc_id   = struct1%phosphorus%owner_proc_id
    structOut%phosphorus%owner_tile_path = struct1%phosphorus%owner_tile_path
    ALLOCATE(structOut%phosphorus%ptr2d, source=struct1%phosphorus%ptr2d)
    structOut%carbon13%owner_model_id    = struct1%carbon13%owner_model_id
    structOut%carbon13%owner_proc_id     = struct1%carbon13%owner_proc_id
    structOut%carbon13%owner_tile_path   = struct1%carbon13%owner_tile_path
    ALLOCATE(structOut%carbon13%ptr2d,   source=struct1%carbon13%ptr2d)
    structOut%carbon14%owner_model_id    = struct1%carbon14%owner_model_id
    structOut%carbon14%owner_proc_id     = struct1%carbon14%owner_proc_id
    structOut%carbon14%owner_tile_path   = struct1%carbon14%owner_tile_path
    ALLOCATE(structOut%carbon14%ptr2d,   source=struct1%carbon14%ptr2d)
    structOut%nitrogen15%owner_model_id  = struct1%nitrogen15%owner_model_id
    structOut%nitrogen15%owner_proc_id   = struct1%nitrogen15%owner_proc_id
    structOut%nitrogen15%owner_tile_path = struct1%nitrogen15%owner_tile_path
    ALLOCATE(structOut%nitrogen15%ptr2d, source=struct1%nitrogen15%ptr2d)

    ! run calculation
    ! the differentiation between ON_CHUNK/ON_DOMAIN is done at the level of variables (i.e., the op_ov routines of the jsb_var)
    structOut%carbon          = struct1%carbon      + struct2%carbon
    structOut%nitrogen        = struct1%nitrogen    + struct2%nitrogen  
    structOut%phosphorus      = struct1%phosphorus  + struct2%phosphorus
    structOut%carbon13        = struct1%carbon13    + struct2%carbon13  
    structOut%carbon14        = struct1%carbon14    + struct2%carbon14  
    structOut%nitrogen15      = struct1%nitrogen15  + struct2%nitrogen15

    ! is this needed here too (not only with the jsb_var op_ov)
    structOut%carbon%ptr      => structOut%carbon%ptr2d(:,:)
    structOut%nitrogen%ptr    => structOut%nitrogen%ptr2d(:,:)
    structOut%phosphorus%ptr  => structOut%phosphorus%ptr2d(:,:)
    structOut%carbon13%ptr    => structOut%carbon13%ptr2d(:,:)
    structOut%carbon14%ptr    => structOut%carbon14%ptr2d(:,:)
    structOut%nitrogen15%ptr  => structOut%nitrogen15%ptr2d(:,:)
  
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
    TYPE(t_veg_elements_pool),      TARGET      :: structOut
  
    CHARACTER(len=*), PARAMETER :: routine = modname//':Subtract_veg_elem_p_minus_scalar'
  
    ! init temporary pool structOut
    CALL structOut%Init(name='op_ov_help', shortname='op_ov_help', l_elements=.TRUE., element_unit=TRIM(struct1%element_unit))
    structOut%carbon%owner_model_id      = struct1%carbon%owner_model_id
    structOut%carbon%owner_proc_id       = struct1%carbon%owner_proc_id
    structOut%carbon%owner_tile_path     = struct1%carbon%owner_tile_path
    ALLOCATE(structOut%carbon%ptr2d,     source=struct1%carbon%ptr2d)
    structOut%nitrogen%owner_model_id    = struct1%nitrogen%owner_model_id
    structOut%nitrogen%owner_proc_id     = struct1%nitrogen%owner_proc_id
    structOut%nitrogen%owner_tile_path   = struct1%nitrogen%owner_tile_path
    ALLOCATE(structOut%nitrogen%ptr2d,   source=struct1%nitrogen%ptr2d)
    structOut%phosphorus%owner_model_id  = struct1%phosphorus%owner_model_id
    structOut%phosphorus%owner_proc_id   = struct1%phosphorus%owner_proc_id
    structOut%phosphorus%owner_tile_path = struct1%phosphorus%owner_tile_path
    ALLOCATE(structOut%phosphorus%ptr2d, source=struct1%phosphorus%ptr2d)
    structOut%carbon13%owner_model_id    = struct1%carbon13%owner_model_id
    structOut%carbon13%owner_proc_id     = struct1%carbon13%owner_proc_id
    structOut%carbon13%owner_tile_path   = struct1%carbon13%owner_tile_path
    ALLOCATE(structOut%carbon13%ptr2d,   source=struct1%carbon13%ptr2d)
    structOut%carbon14%owner_model_id    = struct1%carbon14%owner_model_id
    structOut%carbon14%owner_proc_id     = struct1%carbon14%owner_proc_id
    structOut%carbon14%owner_tile_path   = struct1%carbon14%owner_tile_path
    ALLOCATE(structOut%carbon14%ptr2d,   source=struct1%carbon14%ptr2d)
    structOut%nitrogen15%owner_model_id  = struct1%nitrogen15%owner_model_id
    structOut%nitrogen15%owner_proc_id   = struct1%nitrogen15%owner_proc_id
    structOut%nitrogen15%owner_tile_path = struct1%nitrogen15%owner_tile_path
    ALLOCATE(structOut%nitrogen15%ptr2d, source=struct1%nitrogen15%ptr2d)

    ! run calculation
    ! the differentiation between ON_CHUNK/ON_DOMAIN is done at the level of variables (i.e., the op_ov routines of the jsb_var)
    structOut%carbon          = struct1%carbon      - scalar
    structOut%nitrogen        = struct1%nitrogen    - scalar
    structOut%phosphorus      = struct1%phosphorus  - scalar
    structOut%carbon13        = struct1%carbon13    - scalar
    structOut%carbon14        = struct1%carbon14    - scalar
    structOut%nitrogen15      = struct1%nitrogen15  - scalar

    ! is this needed here too (not only with the jsb_var op_ov)
    structOut%carbon%ptr      => structOut%carbon%ptr2d(:,:)
    structOut%nitrogen%ptr    => structOut%nitrogen%ptr2d(:,:)
    structOut%phosphorus%ptr  => structOut%phosphorus%ptr2d(:,:)
    structOut%carbon13%ptr    => structOut%carbon13%ptr2d(:,:)
    structOut%carbon14%ptr    => structOut%carbon14%ptr2d(:,:)
    structOut%nitrogen15%ptr  => structOut%nitrogen15%ptr2d(:,:)
  
  END FUNCTION Subtract_veg_elem_p_minus_scalar
  !--------------------------------------------------------------------------------------------------------------------------------
  ! t_veg_elements_pool - field1d(:)        | ON_CHUNK
  !--------------------------------------------------------------------------------------------------------------------------------
  FUNCTION Subtract_veg_elem_p_minus_field1d(struct1, field1d) RESULT(structOut)
  
    CLASS(t_veg_elements_pool),     INTENT(in)  :: struct1
    REAL(wp),                       INTENT(in)  :: field1d(:)
    TYPE(t_veg_elements_pool),      TARGET      :: structOut
  
    CHARACTER(len=*), PARAMETER :: routine = modname//':Subtract_veg_elem_p_minus_field1d'
  
    ! init temporary pool structOut
    CALL structOut%Init(name='op_ov_help', shortname='op_ov_help', l_elements=.TRUE., element_unit=TRIM(struct1%element_unit))
    structOut%carbon%owner_model_id      = struct1%carbon%owner_model_id
    structOut%carbon%owner_proc_id       = struct1%carbon%owner_proc_id
    structOut%carbon%owner_tile_path     = struct1%carbon%owner_tile_path
    ALLOCATE(structOut%carbon%ptr2d,     source=struct1%carbon%ptr2d)
    structOut%nitrogen%owner_model_id    = struct1%nitrogen%owner_model_id
    structOut%nitrogen%owner_proc_id     = struct1%nitrogen%owner_proc_id
    structOut%nitrogen%owner_tile_path   = struct1%nitrogen%owner_tile_path
    ALLOCATE(structOut%nitrogen%ptr2d,   source=struct1%nitrogen%ptr2d)
    structOut%phosphorus%owner_model_id  = struct1%phosphorus%owner_model_id
    structOut%phosphorus%owner_proc_id   = struct1%phosphorus%owner_proc_id
    structOut%phosphorus%owner_tile_path = struct1%phosphorus%owner_tile_path
    ALLOCATE(structOut%phosphorus%ptr2d, source=struct1%phosphorus%ptr2d)
    structOut%carbon13%owner_model_id    = struct1%carbon13%owner_model_id
    structOut%carbon13%owner_proc_id     = struct1%carbon13%owner_proc_id
    structOut%carbon13%owner_tile_path   = struct1%carbon13%owner_tile_path
    ALLOCATE(structOut%carbon13%ptr2d,   source=struct1%carbon13%ptr2d)
    structOut%carbon14%owner_model_id    = struct1%carbon14%owner_model_id
    structOut%carbon14%owner_proc_id     = struct1%carbon14%owner_proc_id
    structOut%carbon14%owner_tile_path   = struct1%carbon14%owner_tile_path
    ALLOCATE(structOut%carbon14%ptr2d,   source=struct1%carbon14%ptr2d)
    structOut%nitrogen15%owner_model_id  = struct1%nitrogen15%owner_model_id
    structOut%nitrogen15%owner_proc_id   = struct1%nitrogen15%owner_proc_id
    structOut%nitrogen15%owner_tile_path = struct1%nitrogen15%owner_tile_path
    ALLOCATE(structOut%nitrogen15%ptr2d, source=struct1%nitrogen15%ptr2d)

    ! run calculation
    ! the differentiation between ON_CHUNK/ON_DOMAIN is done at the level of variables (i.e., the op_ov routines of the jsb_var)
    structOut%carbon          = struct1%carbon      - field1d(:)
    structOut%nitrogen        = struct1%nitrogen    - field1d(:)
    structOut%phosphorus      = struct1%phosphorus  - field1d(:)
    structOut%carbon13        = struct1%carbon13    - field1d(:)
    structOut%carbon14        = struct1%carbon14    - field1d(:)
    structOut%nitrogen15      = struct1%nitrogen15  - field1d(:)

    ! is this needed here too (not only with the jsb_var op_ov)
    structOut%carbon%ptr      => structOut%carbon%ptr2d(:,:)
    structOut%nitrogen%ptr    => structOut%nitrogen%ptr2d(:,:)
    structOut%phosphorus%ptr  => structOut%phosphorus%ptr2d(:,:)
    structOut%carbon13%ptr    => structOut%carbon13%ptr2d(:,:)
    structOut%carbon14%ptr    => structOut%carbon14%ptr2d(:,:)
    structOut%nitrogen15%ptr  => structOut%nitrogen15%ptr2d(:,:)
  
  END FUNCTION Subtract_veg_elem_p_minus_field1d
  !--------------------------------------------------------------------------------------------------------------------------------
  ! t_veg_elements_pool - field2d(:,:)        | ON_CHUNK & ON_DOMAIN | field2d(:,1) & field2d(:,:)
  !--------------------------------------------------------------------------------------------------------------------------------
  FUNCTION Subtract_veg_elem_p_minus_field2d(struct1, field2d) RESULT(structOut)
  
    CLASS(t_veg_elements_pool),     INTENT(in)  :: struct1
    REAL(wp),                       INTENT(in)  :: field2d(:,:)
    TYPE(t_veg_elements_pool),      TARGET      :: structOut
  
    CHARACTER(len=*), PARAMETER :: routine = modname//':Subtract_veg_elem_p_minus_field2d'
  
    ! init temporary pool structOut
    CALL structOut%Init(name='op_ov_help', shortname='op_ov_help', l_elements=.TRUE., element_unit=TRIM(struct1%element_unit))
    structOut%carbon%owner_model_id      = struct1%carbon%owner_model_id
    structOut%carbon%owner_proc_id       = struct1%carbon%owner_proc_id
    structOut%carbon%owner_tile_path     = struct1%carbon%owner_tile_path
    ALLOCATE(structOut%carbon%ptr2d,     source=struct1%carbon%ptr2d)
    structOut%nitrogen%owner_model_id    = struct1%nitrogen%owner_model_id
    structOut%nitrogen%owner_proc_id     = struct1%nitrogen%owner_proc_id
    structOut%nitrogen%owner_tile_path   = struct1%nitrogen%owner_tile_path
    ALLOCATE(structOut%nitrogen%ptr2d,   source=struct1%nitrogen%ptr2d)
    structOut%phosphorus%owner_model_id  = struct1%phosphorus%owner_model_id
    structOut%phosphorus%owner_proc_id   = struct1%phosphorus%owner_proc_id
    structOut%phosphorus%owner_tile_path = struct1%phosphorus%owner_tile_path
    ALLOCATE(structOut%phosphorus%ptr2d, source=struct1%phosphorus%ptr2d)
    structOut%carbon13%owner_model_id    = struct1%carbon13%owner_model_id
    structOut%carbon13%owner_proc_id     = struct1%carbon13%owner_proc_id
    structOut%carbon13%owner_tile_path   = struct1%carbon13%owner_tile_path
    ALLOCATE(structOut%carbon13%ptr2d,   source=struct1%carbon13%ptr2d)
    structOut%carbon14%owner_model_id    = struct1%carbon14%owner_model_id
    structOut%carbon14%owner_proc_id     = struct1%carbon14%owner_proc_id
    structOut%carbon14%owner_tile_path   = struct1%carbon14%owner_tile_path
    ALLOCATE(structOut%carbon14%ptr2d,   source=struct1%carbon14%ptr2d)
    structOut%nitrogen15%owner_model_id  = struct1%nitrogen15%owner_model_id
    structOut%nitrogen15%owner_proc_id   = struct1%nitrogen15%owner_proc_id
    structOut%nitrogen15%owner_tile_path = struct1%nitrogen15%owner_tile_path
    ALLOCATE(structOut%nitrogen15%ptr2d, source=struct1%nitrogen15%ptr2d)

    ! run calculation
    ! the differentiation between ON_CHUNK/ON_DOMAIN is done at the level of variables (i.e., the op_ov routines of the jsb_var)
    structOut%carbon          = struct1%carbon      - field2d(:,:)
    structOut%nitrogen        = struct1%nitrogen    - field2d(:,:)
    structOut%phosphorus      = struct1%phosphorus  - field2d(:,:)
    structOut%carbon13        = struct1%carbon13    - field2d(:,:)
    structOut%carbon14        = struct1%carbon14    - field2d(:,:)
    structOut%nitrogen15      = struct1%nitrogen15  - field2d(:,:)

    ! is this needed here too (not only with the jsb_var op_ov)
    structOut%carbon%ptr      => structOut%carbon%ptr2d(:,:)
    structOut%nitrogen%ptr    => structOut%nitrogen%ptr2d(:,:)
    structOut%phosphorus%ptr  => structOut%phosphorus%ptr2d(:,:)
    structOut%carbon13%ptr    => structOut%carbon13%ptr2d(:,:)
    structOut%carbon14%ptr    => structOut%carbon14%ptr2d(:,:)
    structOut%nitrogen15%ptr  => structOut%nitrogen15%ptr2d(:,:)
  
  END FUNCTION Subtract_veg_elem_p_minus_field2d
  !--------------------------------------------------------------------------------------------------------------------------------
  ! t_veg_elements_pool - t_veg_elements_pool | ON_CHUNK & ON_DOMAIN
  !--------------------------------------------------------------------------------------------------------------------------------
  FUNCTION Subtract_veg_elem_p_minus_veg_elem_p(struct1, struct2) RESULT(structOut)
  
    CLASS(t_veg_elements_pool),     INTENT(in)  :: struct1, struct2
    TYPE(t_veg_elements_pool),      TARGET      :: structOut
  
    CHARACTER(len=*), PARAMETER :: routine = modname//':Subtract_veg_elem_p_minus_veg_elem_p'
  
    ! init temporary pool structOut
    CALL structOut%Init(name='op_ov_help', shortname='op_ov_help', l_elements=.TRUE., element_unit=TRIM(struct1%element_unit))
    structOut%carbon%owner_model_id      = struct1%carbon%owner_model_id
    structOut%carbon%owner_proc_id       = struct1%carbon%owner_proc_id
    structOut%carbon%owner_tile_path     = struct1%carbon%owner_tile_path
    ALLOCATE(structOut%carbon%ptr2d,     source=struct1%carbon%ptr2d)
    structOut%nitrogen%owner_model_id    = struct1%nitrogen%owner_model_id
    structOut%nitrogen%owner_proc_id     = struct1%nitrogen%owner_proc_id
    structOut%nitrogen%owner_tile_path   = struct1%nitrogen%owner_tile_path
    ALLOCATE(structOut%nitrogen%ptr2d,   source=struct1%nitrogen%ptr2d)
    structOut%phosphorus%owner_model_id  = struct1%phosphorus%owner_model_id
    structOut%phosphorus%owner_proc_id   = struct1%phosphorus%owner_proc_id
    structOut%phosphorus%owner_tile_path = struct1%phosphorus%owner_tile_path
    ALLOCATE(structOut%phosphorus%ptr2d, source=struct1%phosphorus%ptr2d)
    structOut%carbon13%owner_model_id    = struct1%carbon13%owner_model_id
    structOut%carbon13%owner_proc_id     = struct1%carbon13%owner_proc_id
    structOut%carbon13%owner_tile_path   = struct1%carbon13%owner_tile_path
    ALLOCATE(structOut%carbon13%ptr2d,   source=struct1%carbon13%ptr2d)
    structOut%carbon14%owner_model_id    = struct1%carbon14%owner_model_id
    structOut%carbon14%owner_proc_id     = struct1%carbon14%owner_proc_id
    structOut%carbon14%owner_tile_path   = struct1%carbon14%owner_tile_path
    ALLOCATE(structOut%carbon14%ptr2d,   source=struct1%carbon14%ptr2d)
    structOut%nitrogen15%owner_model_id  = struct1%nitrogen15%owner_model_id
    structOut%nitrogen15%owner_proc_id   = struct1%nitrogen15%owner_proc_id
    structOut%nitrogen15%owner_tile_path = struct1%nitrogen15%owner_tile_path
    ALLOCATE(structOut%nitrogen15%ptr2d, source=struct1%nitrogen15%ptr2d)

    ! run calculation
    ! the differentiation between ON_CHUNK/ON_DOMAIN is done at the level of variables (i.e., the op_ov routines of the jsb_var)
    structOut%carbon          = struct1%carbon      - struct2%carbon
    structOut%nitrogen        = struct1%nitrogen    - struct2%nitrogen  
    structOut%phosphorus      = struct1%phosphorus  - struct2%phosphorus
    structOut%carbon13        = struct1%carbon13    - struct2%carbon13  
    structOut%carbon14        = struct1%carbon14    - struct2%carbon14  
    structOut%nitrogen15      = struct1%nitrogen15  - struct2%nitrogen15

    ! is this needed here too (not only with the jsb_var op_ov)
    structOut%carbon%ptr      => structOut%carbon%ptr2d(:,:)
    structOut%nitrogen%ptr    => structOut%nitrogen%ptr2d(:,:)
    structOut%phosphorus%ptr  => structOut%phosphorus%ptr2d(:,:)
    structOut%carbon13%ptr    => structOut%carbon13%ptr2d(:,:)
    structOut%carbon14%ptr    => structOut%carbon14%ptr2d(:,:)
    structOut%nitrogen15%ptr  => structOut%nitrogen15%ptr2d(:,:)
  
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
    TYPE(t_veg_elements_pool),      TARGET      :: structOut
  
    CHARACTER(len=*), PARAMETER :: routine = modname//':Multiply_veg_elem_p_with_scalar'
  
    ! init temporary pool structOut
    CALL structOut%Init(name='op_ov_help', shortname='op_ov_help', l_elements=.TRUE., element_unit=TRIM(struct1%element_unit))
    structOut%carbon%owner_model_id      = struct1%carbon%owner_model_id
    structOut%carbon%owner_proc_id       = struct1%carbon%owner_proc_id
    structOut%carbon%owner_tile_path     = struct1%carbon%owner_tile_path
    ALLOCATE(structOut%carbon%ptr2d,     source=struct1%carbon%ptr2d)
    structOut%nitrogen%owner_model_id    = struct1%nitrogen%owner_model_id
    structOut%nitrogen%owner_proc_id     = struct1%nitrogen%owner_proc_id
    structOut%nitrogen%owner_tile_path   = struct1%nitrogen%owner_tile_path
    ALLOCATE(structOut%nitrogen%ptr2d,   source=struct1%nitrogen%ptr2d)
    structOut%phosphorus%owner_model_id  = struct1%phosphorus%owner_model_id
    structOut%phosphorus%owner_proc_id   = struct1%phosphorus%owner_proc_id
    structOut%phosphorus%owner_tile_path = struct1%phosphorus%owner_tile_path
    ALLOCATE(structOut%phosphorus%ptr2d, source=struct1%phosphorus%ptr2d)
    structOut%carbon13%owner_model_id    = struct1%carbon13%owner_model_id
    structOut%carbon13%owner_proc_id     = struct1%carbon13%owner_proc_id
    structOut%carbon13%owner_tile_path   = struct1%carbon13%owner_tile_path
    ALLOCATE(structOut%carbon13%ptr2d,   source=struct1%carbon13%ptr2d)
    structOut%carbon14%owner_model_id    = struct1%carbon14%owner_model_id
    structOut%carbon14%owner_proc_id     = struct1%carbon14%owner_proc_id
    structOut%carbon14%owner_tile_path   = struct1%carbon14%owner_tile_path
    ALLOCATE(structOut%carbon14%ptr2d,   source=struct1%carbon14%ptr2d)
    structOut%nitrogen15%owner_model_id  = struct1%nitrogen15%owner_model_id
    structOut%nitrogen15%owner_proc_id   = struct1%nitrogen15%owner_proc_id
    structOut%nitrogen15%owner_tile_path = struct1%nitrogen15%owner_tile_path
    ALLOCATE(structOut%nitrogen15%ptr2d, source=struct1%nitrogen15%ptr2d)

    ! run calculation
    ! the differentiation between ON_CHUNK/ON_DOMAIN is done at the level of variables (i.e., the op_ov routines of the jsb_var)
    structOut%carbon          = struct1%carbon      * scalar
    structOut%nitrogen        = struct1%nitrogen    * scalar
    structOut%phosphorus      = struct1%phosphorus  * scalar
    structOut%carbon13        = struct1%carbon13    * scalar
    structOut%carbon14        = struct1%carbon14    * scalar
    structOut%nitrogen15      = struct1%nitrogen15  * scalar

    ! is this needed here too (not only with the jsb_var op_ov)
    structOut%carbon%ptr      => structOut%carbon%ptr2d(:,:)
    structOut%nitrogen%ptr    => structOut%nitrogen%ptr2d(:,:)
    structOut%phosphorus%ptr  => structOut%phosphorus%ptr2d(:,:)
    structOut%carbon13%ptr    => structOut%carbon13%ptr2d(:,:)
    structOut%carbon14%ptr    => structOut%carbon14%ptr2d(:,:)
    structOut%nitrogen15%ptr  => structOut%nitrogen15%ptr2d(:,:)
  
  END FUNCTION Multiply_veg_elem_p_with_scalar
  !--------------------------------------------------------------------------------------------------------------------------------
  ! scalar * t_veg_elements_pool        | ON_CHUNK & ON_DOMAIN
  !--------------------------------------------------------------------------------------------------------------------------------
  FUNCTION Multiply_scalar_with_veg_elem_p(scalar, struct1) RESULT(structOut)
  
    REAL(wp),                       INTENT(in)  :: scalar
    CLASS(t_veg_elements_pool),     INTENT(in)  :: struct1
    TYPE(t_veg_elements_pool),      TARGET      :: structOut
  
    CHARACTER(len=*), PARAMETER :: routine = modname//':Multiply_scalar_with_veg_elem_p'
  
    ! init temporary pool structOut
    CALL structOut%Init(name='op_ov_help', shortname='op_ov_help', l_elements=.TRUE., element_unit=TRIM(struct1%element_unit))
    structOut%carbon%owner_model_id      = struct1%carbon%owner_model_id
    structOut%carbon%owner_proc_id       = struct1%carbon%owner_proc_id
    structOut%carbon%owner_tile_path     = struct1%carbon%owner_tile_path
    ALLOCATE(structOut%carbon%ptr2d,     source=struct1%carbon%ptr2d)
    structOut%nitrogen%owner_model_id    = struct1%nitrogen%owner_model_id
    structOut%nitrogen%owner_proc_id     = struct1%nitrogen%owner_proc_id
    structOut%nitrogen%owner_tile_path   = struct1%nitrogen%owner_tile_path
    ALLOCATE(structOut%nitrogen%ptr2d,   source=struct1%nitrogen%ptr2d)
    structOut%phosphorus%owner_model_id  = struct1%phosphorus%owner_model_id
    structOut%phosphorus%owner_proc_id   = struct1%phosphorus%owner_proc_id
    structOut%phosphorus%owner_tile_path = struct1%phosphorus%owner_tile_path
    ALLOCATE(structOut%phosphorus%ptr2d, source=struct1%phosphorus%ptr2d)
    structOut%carbon13%owner_model_id    = struct1%carbon13%owner_model_id
    structOut%carbon13%owner_proc_id     = struct1%carbon13%owner_proc_id
    structOut%carbon13%owner_tile_path   = struct1%carbon13%owner_tile_path
    ALLOCATE(structOut%carbon13%ptr2d,   source=struct1%carbon13%ptr2d)
    structOut%carbon14%owner_model_id    = struct1%carbon14%owner_model_id
    structOut%carbon14%owner_proc_id     = struct1%carbon14%owner_proc_id
    structOut%carbon14%owner_tile_path   = struct1%carbon14%owner_tile_path
    ALLOCATE(structOut%carbon14%ptr2d,   source=struct1%carbon14%ptr2d)
    structOut%nitrogen15%owner_model_id  = struct1%nitrogen15%owner_model_id
    structOut%nitrogen15%owner_proc_id   = struct1%nitrogen15%owner_proc_id
    structOut%nitrogen15%owner_tile_path = struct1%nitrogen15%owner_tile_path
    ALLOCATE(structOut%nitrogen15%ptr2d, source=struct1%nitrogen15%ptr2d)

    ! run calculation
    ! the differentiation between ON_CHUNK/ON_DOMAIN is done at the level of variables (i.e., the op_ov routines of the jsb_var)
    structOut%carbon          = scalar * struct1%carbon    
    structOut%nitrogen        = scalar * struct1%nitrogen  
    structOut%phosphorus      = scalar * struct1%phosphorus
    structOut%carbon13        = scalar * struct1%carbon13  
    structOut%carbon14        = scalar * struct1%carbon14  
    structOut%nitrogen15      = scalar * struct1%nitrogen15

    ! is this needed here too (not only with the jsb_var op_ov)
    structOut%carbon%ptr      => structOut%carbon%ptr2d(:,:)
    structOut%nitrogen%ptr    => structOut%nitrogen%ptr2d(:,:)
    structOut%phosphorus%ptr  => structOut%phosphorus%ptr2d(:,:)
    structOut%carbon13%ptr    => structOut%carbon13%ptr2d(:,:)
    structOut%carbon14%ptr    => structOut%carbon14%ptr2d(:,:)
    structOut%nitrogen15%ptr  => structOut%nitrogen15%ptr2d(:,:)
  
  END FUNCTION Multiply_scalar_with_veg_elem_p
  !--------------------------------------------------------------------------------------------------------------------------------
  ! t_veg_elements_pool * field1d(:)        | ON_CHUNK
  !--------------------------------------------------------------------------------------------------------------------------------
  FUNCTION Multiply_veg_elem_p_with_field1d(struct1, field1d) RESULT(structOut)
  
    CLASS(t_veg_elements_pool),     INTENT(in)  :: struct1
    REAL(wp),                       INTENT(in)  :: field1d(:)
    TYPE(t_veg_elements_pool),      TARGET      :: structOut
  
    CHARACTER(len=*), PARAMETER :: routine = modname//':Multiply_veg_elem_p_with_field1d'
  
    ! init temporary pool structOut
    CALL structOut%Init(name='op_ov_help', shortname='op_ov_help', l_elements=.TRUE., element_unit=TRIM(struct1%element_unit))
    structOut%carbon%owner_model_id      = struct1%carbon%owner_model_id
    structOut%carbon%owner_proc_id       = struct1%carbon%owner_proc_id
    structOut%carbon%owner_tile_path     = struct1%carbon%owner_tile_path
    ALLOCATE(structOut%carbon%ptr2d,     source=struct1%carbon%ptr2d)
    structOut%nitrogen%owner_model_id    = struct1%nitrogen%owner_model_id
    structOut%nitrogen%owner_proc_id     = struct1%nitrogen%owner_proc_id
    structOut%nitrogen%owner_tile_path   = struct1%nitrogen%owner_tile_path
    ALLOCATE(structOut%nitrogen%ptr2d,   source=struct1%nitrogen%ptr2d)
    structOut%phosphorus%owner_model_id  = struct1%phosphorus%owner_model_id
    structOut%phosphorus%owner_proc_id   = struct1%phosphorus%owner_proc_id
    structOut%phosphorus%owner_tile_path = struct1%phosphorus%owner_tile_path
    ALLOCATE(structOut%phosphorus%ptr2d, source=struct1%phosphorus%ptr2d)
    structOut%carbon13%owner_model_id    = struct1%carbon13%owner_model_id
    structOut%carbon13%owner_proc_id     = struct1%carbon13%owner_proc_id
    structOut%carbon13%owner_tile_path   = struct1%carbon13%owner_tile_path
    ALLOCATE(structOut%carbon13%ptr2d,   source=struct1%carbon13%ptr2d)
    structOut%carbon14%owner_model_id    = struct1%carbon14%owner_model_id
    structOut%carbon14%owner_proc_id     = struct1%carbon14%owner_proc_id
    structOut%carbon14%owner_tile_path   = struct1%carbon14%owner_tile_path
    ALLOCATE(structOut%carbon14%ptr2d,   source=struct1%carbon14%ptr2d)
    structOut%nitrogen15%owner_model_id  = struct1%nitrogen15%owner_model_id
    structOut%nitrogen15%owner_proc_id   = struct1%nitrogen15%owner_proc_id
    structOut%nitrogen15%owner_tile_path = struct1%nitrogen15%owner_tile_path
    ALLOCATE(structOut%nitrogen15%ptr2d, source=struct1%nitrogen15%ptr2d)

    ! run calculation
    ! the differentiation between ON_CHUNK/ON_DOMAIN is done at the level of variables (i.e., the op_ov routines of the jsb_var)
    structOut%carbon          = struct1%carbon      * field1d(:)
    structOut%nitrogen        = struct1%nitrogen    * field1d(:)
    structOut%phosphorus      = struct1%phosphorus  * field1d(:)
    structOut%carbon13        = struct1%carbon13    * field1d(:)
    structOut%carbon14        = struct1%carbon14    * field1d(:)
    structOut%nitrogen15      = struct1%nitrogen15  * field1d(:)

    ! is this needed here too (not only with the jsb_var op_ov)
    structOut%carbon%ptr      => structOut%carbon%ptr2d(:,:)
    structOut%nitrogen%ptr    => structOut%nitrogen%ptr2d(:,:)
    structOut%phosphorus%ptr  => structOut%phosphorus%ptr2d(:,:)
    structOut%carbon13%ptr    => structOut%carbon13%ptr2d(:,:)
    structOut%carbon14%ptr    => structOut%carbon14%ptr2d(:,:)
    structOut%nitrogen15%ptr  => structOut%nitrogen15%ptr2d(:,:)
  
  END FUNCTION Multiply_veg_elem_p_with_field1d
  !--------------------------------------------------------------------------------------------------------------------------------
  ! field1d(:) * t_veg_elements_pool        | ON_CHUNK
  !--------------------------------------------------------------------------------------------------------------------------------
  FUNCTION Multiply_field1d_with_veg_elem_p(field1d, struct1) RESULT(structOut)
  
    REAL(wp),                       INTENT(in)  :: field1d(:)
    CLASS(t_veg_elements_pool),     INTENT(in)  :: struct1
    TYPE(t_veg_elements_pool),      TARGET      :: structOut
  
    CHARACTER(len=*), PARAMETER :: routine = modname//':Multiply_field1d_with_veg_elem_p'
  
    ! init temporary pool structOut
    CALL structOut%Init(name='op_ov_help', shortname='op_ov_help', l_elements=.TRUE., element_unit=TRIM(struct1%element_unit))
    structOut%carbon%owner_model_id      = struct1%carbon%owner_model_id
    structOut%carbon%owner_proc_id       = struct1%carbon%owner_proc_id
    structOut%carbon%owner_tile_path     = struct1%carbon%owner_tile_path
    ALLOCATE(structOut%carbon%ptr2d,     source=struct1%carbon%ptr2d)
    structOut%nitrogen%owner_model_id    = struct1%nitrogen%owner_model_id
    structOut%nitrogen%owner_proc_id     = struct1%nitrogen%owner_proc_id
    structOut%nitrogen%owner_tile_path   = struct1%nitrogen%owner_tile_path
    ALLOCATE(structOut%nitrogen%ptr2d,   source=struct1%nitrogen%ptr2d)
    structOut%phosphorus%owner_model_id  = struct1%phosphorus%owner_model_id
    structOut%phosphorus%owner_proc_id   = struct1%phosphorus%owner_proc_id
    structOut%phosphorus%owner_tile_path = struct1%phosphorus%owner_tile_path
    ALLOCATE(structOut%phosphorus%ptr2d, source=struct1%phosphorus%ptr2d)
    structOut%carbon13%owner_model_id    = struct1%carbon13%owner_model_id
    structOut%carbon13%owner_proc_id     = struct1%carbon13%owner_proc_id
    structOut%carbon13%owner_tile_path   = struct1%carbon13%owner_tile_path
    ALLOCATE(structOut%carbon13%ptr2d,   source=struct1%carbon13%ptr2d)
    structOut%carbon14%owner_model_id    = struct1%carbon14%owner_model_id
    structOut%carbon14%owner_proc_id     = struct1%carbon14%owner_proc_id
    structOut%carbon14%owner_tile_path   = struct1%carbon14%owner_tile_path
    ALLOCATE(structOut%carbon14%ptr2d,   source=struct1%carbon14%ptr2d)
    structOut%nitrogen15%owner_model_id  = struct1%nitrogen15%owner_model_id
    structOut%nitrogen15%owner_proc_id   = struct1%nitrogen15%owner_proc_id
    structOut%nitrogen15%owner_tile_path = struct1%nitrogen15%owner_tile_path
    ALLOCATE(structOut%nitrogen15%ptr2d, source=struct1%nitrogen15%ptr2d)

    ! run calculation
    ! the differentiation between ON_CHUNK/ON_DOMAIN is done at the level of variables (i.e., the op_ov routines of the jsb_var)
    structOut%carbon          = field1d(:) * struct1%carbon    
    structOut%nitrogen        = field1d(:) * struct1%nitrogen  
    structOut%phosphorus      = field1d(:) * struct1%phosphorus
    structOut%carbon13        = field1d(:) * struct1%carbon13  
    structOut%carbon14        = field1d(:) * struct1%carbon14  
    structOut%nitrogen15      = field1d(:) * struct1%nitrogen15

    ! is this needed here too (not only with the jsb_var op_ov)
    structOut%carbon%ptr      => structOut%carbon%ptr2d(:,:)
    structOut%nitrogen%ptr    => structOut%nitrogen%ptr2d(:,:)
    structOut%phosphorus%ptr  => structOut%phosphorus%ptr2d(:,:)
    structOut%carbon13%ptr    => structOut%carbon13%ptr2d(:,:)
    structOut%carbon14%ptr    => structOut%carbon14%ptr2d(:,:)
    structOut%nitrogen15%ptr  => structOut%nitrogen15%ptr2d(:,:)
  
  END FUNCTION Multiply_field1d_with_veg_elem_p
  !--------------------------------------------------------------------------------------------------------------------------------
  ! t_veg_elements_pool * field2d(:,:)        | ON_CHUNK & ON_DOMAIN | field2d(:,1) & field2d(:,:)
  !--------------------------------------------------------------------------------------------------------------------------------
  FUNCTION Multiply_veg_elem_p_with_field2d(struct1, field2d) RESULT(structOut)
  
    CLASS(t_veg_elements_pool),     INTENT(in)  :: struct1
    REAL(wp),                       INTENT(in)  :: field2d(:,:)
    TYPE(t_veg_elements_pool),      TARGET      :: structOut
  
    CHARACTER(len=*), PARAMETER :: routine = modname//':Multiply_veg_elem_p_with_field2d'
  
    ! init temporary pool structOut
    CALL structOut%Init(name='op_ov_help', shortname='op_ov_help', l_elements=.TRUE., element_unit=TRIM(struct1%element_unit))
    structOut%carbon%owner_model_id      = struct1%carbon%owner_model_id
    structOut%carbon%owner_proc_id       = struct1%carbon%owner_proc_id
    structOut%carbon%owner_tile_path     = struct1%carbon%owner_tile_path
    ALLOCATE(structOut%carbon%ptr2d,     source=struct1%carbon%ptr2d)
    structOut%nitrogen%owner_model_id    = struct1%nitrogen%owner_model_id
    structOut%nitrogen%owner_proc_id     = struct1%nitrogen%owner_proc_id
    structOut%nitrogen%owner_tile_path   = struct1%nitrogen%owner_tile_path
    ALLOCATE(structOut%nitrogen%ptr2d,   source=struct1%nitrogen%ptr2d)
    structOut%phosphorus%owner_model_id  = struct1%phosphorus%owner_model_id
    structOut%phosphorus%owner_proc_id   = struct1%phosphorus%owner_proc_id
    structOut%phosphorus%owner_tile_path = struct1%phosphorus%owner_tile_path
    ALLOCATE(structOut%phosphorus%ptr2d, source=struct1%phosphorus%ptr2d)
    structOut%carbon13%owner_model_id    = struct1%carbon13%owner_model_id
    structOut%carbon13%owner_proc_id     = struct1%carbon13%owner_proc_id
    structOut%carbon13%owner_tile_path   = struct1%carbon13%owner_tile_path
    ALLOCATE(structOut%carbon13%ptr2d,   source=struct1%carbon13%ptr2d)
    structOut%carbon14%owner_model_id    = struct1%carbon14%owner_model_id
    structOut%carbon14%owner_proc_id     = struct1%carbon14%owner_proc_id
    structOut%carbon14%owner_tile_path   = struct1%carbon14%owner_tile_path
    ALLOCATE(structOut%carbon14%ptr2d,   source=struct1%carbon14%ptr2d)
    structOut%nitrogen15%owner_model_id  = struct1%nitrogen15%owner_model_id
    structOut%nitrogen15%owner_proc_id   = struct1%nitrogen15%owner_proc_id
    structOut%nitrogen15%owner_tile_path = struct1%nitrogen15%owner_tile_path
    ALLOCATE(structOut%nitrogen15%ptr2d, source=struct1%nitrogen15%ptr2d)

    ! run calculation
    ! the differentiation between ON_CHUNK/ON_DOMAIN is done at the level of variables (i.e., the op_ov routines of the jsb_var)
    structOut%carbon          = struct1%carbon      * field2d(:,:)
    structOut%nitrogen        = struct1%nitrogen    * field2d(:,:)
    structOut%phosphorus      = struct1%phosphorus  * field2d(:,:)
    structOut%carbon13        = struct1%carbon13    * field2d(:,:)
    structOut%carbon14        = struct1%carbon14    * field2d(:,:)
    structOut%nitrogen15      = struct1%nitrogen15  * field2d(:,:)

    ! is this needed here too (not only with the jsb_var op_ov)
    structOut%carbon%ptr      => structOut%carbon%ptr2d(:,:)
    structOut%nitrogen%ptr    => structOut%nitrogen%ptr2d(:,:)
    structOut%phosphorus%ptr  => structOut%phosphorus%ptr2d(:,:)
    structOut%carbon13%ptr    => structOut%carbon13%ptr2d(:,:)
    structOut%carbon14%ptr    => structOut%carbon14%ptr2d(:,:)
    structOut%nitrogen15%ptr  => structOut%nitrogen15%ptr2d(:,:)
  
  END FUNCTION Multiply_veg_elem_p_with_field2d
  !--------------------------------------------------------------------------------------------------------------------------------
  ! field2d(:,:) * t_veg_elements_pool        | ON_CHUNK & ON_DOMAIN | field2d(:,1) & field2d(:,:)
  !--------------------------------------------------------------------------------------------------------------------------------
  FUNCTION Multiply_field2d_with_veg_elem_p(field2d, struct1) RESULT(structOut)
  
    REAL(wp),                       INTENT(in)  :: field2d(:,:)
    CLASS(t_veg_elements_pool),     INTENT(in)  :: struct1
    TYPE(t_veg_elements_pool),      TARGET      :: structOut
  
    CHARACTER(len=*), PARAMETER :: routine = modname//':Multiply_field2d_with_veg_elem_p'
  
    ! init temporary pool structOut
    CALL structOut%Init(name='op_ov_help', shortname='op_ov_help', l_elements=.TRUE., element_unit=TRIM(struct1%element_unit))
    structOut%carbon%owner_model_id      = struct1%carbon%owner_model_id
    structOut%carbon%owner_proc_id       = struct1%carbon%owner_proc_id
    structOut%carbon%owner_tile_path     = struct1%carbon%owner_tile_path
    ALLOCATE(structOut%carbon%ptr2d,     source=struct1%carbon%ptr2d)
    structOut%nitrogen%owner_model_id    = struct1%nitrogen%owner_model_id
    structOut%nitrogen%owner_proc_id     = struct1%nitrogen%owner_proc_id
    structOut%nitrogen%owner_tile_path   = struct1%nitrogen%owner_tile_path
    ALLOCATE(structOut%nitrogen%ptr2d,   source=struct1%nitrogen%ptr2d)
    structOut%phosphorus%owner_model_id  = struct1%phosphorus%owner_model_id
    structOut%phosphorus%owner_proc_id   = struct1%phosphorus%owner_proc_id
    structOut%phosphorus%owner_tile_path = struct1%phosphorus%owner_tile_path
    ALLOCATE(structOut%phosphorus%ptr2d, source=struct1%phosphorus%ptr2d)
    structOut%carbon13%owner_model_id    = struct1%carbon13%owner_model_id
    structOut%carbon13%owner_proc_id     = struct1%carbon13%owner_proc_id
    structOut%carbon13%owner_tile_path   = struct1%carbon13%owner_tile_path
    ALLOCATE(structOut%carbon13%ptr2d,   source=struct1%carbon13%ptr2d)
    structOut%carbon14%owner_model_id    = struct1%carbon14%owner_model_id
    structOut%carbon14%owner_proc_id     = struct1%carbon14%owner_proc_id
    structOut%carbon14%owner_tile_path   = struct1%carbon14%owner_tile_path
    ALLOCATE(structOut%carbon14%ptr2d,   source=struct1%carbon14%ptr2d)
    structOut%nitrogen15%owner_model_id  = struct1%nitrogen15%owner_model_id
    structOut%nitrogen15%owner_proc_id   = struct1%nitrogen15%owner_proc_id
    structOut%nitrogen15%owner_tile_path = struct1%nitrogen15%owner_tile_path
    ALLOCATE(structOut%nitrogen15%ptr2d, source=struct1%nitrogen15%ptr2d)

    ! run calculation
    ! the differentiation between ON_CHUNK/ON_DOMAIN is done at the level of variables (i.e., the op_ov routines of the jsb_var)
    structOut%carbon          = field2d(:,:) * struct1%carbon    
    structOut%nitrogen        = field2d(:,:) * struct1%nitrogen  
    structOut%phosphorus      = field2d(:,:) * struct1%phosphorus
    structOut%carbon13        = field2d(:,:) * struct1%carbon13  
    structOut%carbon14        = field2d(:,:) * struct1%carbon14  
    structOut%nitrogen15      = field2d(:,:) * struct1%nitrogen15

    ! is this needed here too (not only with the jsb_var op_ov)
    structOut%carbon%ptr      => structOut%carbon%ptr2d(:,:)
    structOut%nitrogen%ptr    => structOut%nitrogen%ptr2d(:,:)
    structOut%phosphorus%ptr  => structOut%phosphorus%ptr2d(:,:)
    structOut%carbon13%ptr    => structOut%carbon13%ptr2d(:,:)
    structOut%carbon14%ptr    => structOut%carbon14%ptr2d(:,:)
    structOut%nitrogen15%ptr  => structOut%nitrogen15%ptr2d(:,:)
  
  END FUNCTION Multiply_field2d_with_veg_elem_p
  !--------------------------------------------------------------------------------------------------------------------------------
  ! t_veg_elements_pool * t_veg_elements_pool | ON_CHUNK & ON_DOMAIN
  !--------------------------------------------------------------------------------------------------------------------------------
  FUNCTION Multiply_veg_elem_p_with_veg_elem_p(struct1, struct2) RESULT(structOut)
  
    CLASS(t_veg_elements_pool),     INTENT(in)  :: struct1, struct2
    TYPE(t_veg_elements_pool),      TARGET      :: structOut
  
    CHARACTER(len=*), PARAMETER :: routine = modname//':Multiply_veg_elem_p_with_veg_elem_p'
  
    ! init temporary pool structOut
    CALL structOut%Init(name='op_ov_help', shortname='op_ov_help', l_elements=.TRUE., element_unit=TRIM(struct1%element_unit))
    structOut%carbon%owner_model_id      = struct1%carbon%owner_model_id
    structOut%carbon%owner_proc_id       = struct1%carbon%owner_proc_id
    structOut%carbon%owner_tile_path     = struct1%carbon%owner_tile_path
    ALLOCATE(structOut%carbon%ptr2d,     source=struct1%carbon%ptr2d)
    structOut%nitrogen%owner_model_id    = struct1%nitrogen%owner_model_id
    structOut%nitrogen%owner_proc_id     = struct1%nitrogen%owner_proc_id
    structOut%nitrogen%owner_tile_path   = struct1%nitrogen%owner_tile_path
    ALLOCATE(structOut%nitrogen%ptr2d,   source=struct1%nitrogen%ptr2d)
    structOut%phosphorus%owner_model_id  = struct1%phosphorus%owner_model_id
    structOut%phosphorus%owner_proc_id   = struct1%phosphorus%owner_proc_id
    structOut%phosphorus%owner_tile_path = struct1%phosphorus%owner_tile_path
    ALLOCATE(structOut%phosphorus%ptr2d, source=struct1%phosphorus%ptr2d)
    structOut%carbon13%owner_model_id    = struct1%carbon13%owner_model_id
    structOut%carbon13%owner_proc_id     = struct1%carbon13%owner_proc_id
    structOut%carbon13%owner_tile_path   = struct1%carbon13%owner_tile_path
    ALLOCATE(structOut%carbon13%ptr2d,   source=struct1%carbon13%ptr2d)
    structOut%carbon14%owner_model_id    = struct1%carbon14%owner_model_id
    structOut%carbon14%owner_proc_id     = struct1%carbon14%owner_proc_id
    structOut%carbon14%owner_tile_path   = struct1%carbon14%owner_tile_path
    ALLOCATE(structOut%carbon14%ptr2d,   source=struct1%carbon14%ptr2d)
    structOut%nitrogen15%owner_model_id  = struct1%nitrogen15%owner_model_id
    structOut%nitrogen15%owner_proc_id   = struct1%nitrogen15%owner_proc_id
    structOut%nitrogen15%owner_tile_path = struct1%nitrogen15%owner_tile_path
    ALLOCATE(structOut%nitrogen15%ptr2d, source=struct1%nitrogen15%ptr2d)

    ! run calculation
    ! the differentiation between ON_CHUNK/ON_DOMAIN is done at the level of variables (i.e., the op_ov routines of the jsb_var)
    structOut%carbon          = struct1%carbon      * struct2%carbon
    structOut%nitrogen        = struct1%nitrogen    * struct2%nitrogen  
    structOut%phosphorus      = struct1%phosphorus  * struct2%phosphorus
    structOut%carbon13        = struct1%carbon13    * struct2%carbon13  
    structOut%carbon14        = struct1%carbon14    * struct2%carbon14  
    structOut%nitrogen15      = struct1%nitrogen15  * struct2%nitrogen15

    ! is this needed here too (not only with the jsb_var op_ov)
    structOut%carbon%ptr      => structOut%carbon%ptr2d(:,:)
    structOut%nitrogen%ptr    => structOut%nitrogen%ptr2d(:,:)
    structOut%phosphorus%ptr  => structOut%phosphorus%ptr2d(:,:)
    structOut%carbon13%ptr    => structOut%carbon13%ptr2d(:,:)
    structOut%carbon14%ptr    => structOut%carbon14%ptr2d(:,:)
    structOut%nitrogen15%ptr  => structOut%nitrogen15%ptr2d(:,:)
  
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
    TYPE(t_veg_elements_pool),      TARGET      :: structOut
  
    CHARACTER(len=*), PARAMETER :: routine = modname//':Div_veg_elem_p_by_scalar'
  
    ! init temporary pool structOut
    CALL structOut%Init(name='op_ov_help', shortname='op_ov_help', l_elements=.TRUE., element_unit=TRIM(struct1%element_unit))
    structOut%carbon%owner_model_id      = struct1%carbon%owner_model_id
    structOut%carbon%owner_proc_id       = struct1%carbon%owner_proc_id
    structOut%carbon%owner_tile_path     = struct1%carbon%owner_tile_path
    ALLOCATE(structOut%carbon%ptr2d,     source=struct1%carbon%ptr2d)
    structOut%nitrogen%owner_model_id    = struct1%nitrogen%owner_model_id
    structOut%nitrogen%owner_proc_id     = struct1%nitrogen%owner_proc_id
    structOut%nitrogen%owner_tile_path   = struct1%nitrogen%owner_tile_path
    ALLOCATE(structOut%nitrogen%ptr2d,   source=struct1%nitrogen%ptr2d)
    structOut%phosphorus%owner_model_id  = struct1%phosphorus%owner_model_id
    structOut%phosphorus%owner_proc_id   = struct1%phosphorus%owner_proc_id
    structOut%phosphorus%owner_tile_path = struct1%phosphorus%owner_tile_path
    ALLOCATE(structOut%phosphorus%ptr2d, source=struct1%phosphorus%ptr2d)
    structOut%carbon13%owner_model_id    = struct1%carbon13%owner_model_id
    structOut%carbon13%owner_proc_id     = struct1%carbon13%owner_proc_id
    structOut%carbon13%owner_tile_path   = struct1%carbon13%owner_tile_path
    ALLOCATE(structOut%carbon13%ptr2d,   source=struct1%carbon13%ptr2d)
    structOut%carbon14%owner_model_id    = struct1%carbon14%owner_model_id
    structOut%carbon14%owner_proc_id     = struct1%carbon14%owner_proc_id
    structOut%carbon14%owner_tile_path   = struct1%carbon14%owner_tile_path
    ALLOCATE(structOut%carbon14%ptr2d,   source=struct1%carbon14%ptr2d)
    structOut%nitrogen15%owner_model_id  = struct1%nitrogen15%owner_model_id
    structOut%nitrogen15%owner_proc_id   = struct1%nitrogen15%owner_proc_id
    structOut%nitrogen15%owner_tile_path = struct1%nitrogen15%owner_tile_path
    ALLOCATE(structOut%nitrogen15%ptr2d, source=struct1%nitrogen15%ptr2d)

    ! run calculation
    ! the differentiation between ON_CHUNK/ON_DOMAIN is done at the level of variables (i.e., the op_ov routines of the jsb_var)
    structOut%carbon          = struct1%carbon      / scalar
    structOut%nitrogen        = struct1%nitrogen    / scalar
    structOut%phosphorus      = struct1%phosphorus  / scalar
    structOut%carbon13        = struct1%carbon13    / scalar
    structOut%carbon14        = struct1%carbon14    / scalar
    structOut%nitrogen15      = struct1%nitrogen15  / scalar

    ! is this needed here too (not only with the jsb_var op_ov)
    structOut%carbon%ptr      => structOut%carbon%ptr2d(:,:)
    structOut%nitrogen%ptr    => structOut%nitrogen%ptr2d(:,:)
    structOut%phosphorus%ptr  => structOut%phosphorus%ptr2d(:,:)
    structOut%carbon13%ptr    => structOut%carbon13%ptr2d(:,:)
    structOut%carbon14%ptr    => structOut%carbon14%ptr2d(:,:)
    structOut%nitrogen15%ptr  => structOut%nitrogen15%ptr2d(:,:)
  
  END FUNCTION Div_veg_elem_p_by_scalar
  !--------------------------------------------------------------------------------------------------------------------------------
  ! t_veg_elements_pool / field1d(:)        | ON_CHUNK
  !--------------------------------------------------------------------------------------------------------------------------------
  FUNCTION Div_veg_elem_p_by_field1d(struct1, field1d) RESULT(structOut)
  
    CLASS(t_veg_elements_pool),     INTENT(in)  :: struct1
    REAL(wp),                       INTENT(in)  :: field1d(:)
    TYPE(t_veg_elements_pool),      TARGET      :: structOut
  
    CHARACTER(len=*), PARAMETER :: routine = modname//':Div_veg_elem_p_by_field1d'
  
    ! init temporary pool structOut
    CALL structOut%Init(name='op_ov_help', shortname='op_ov_help', l_elements=.TRUE., element_unit=TRIM(struct1%element_unit))
    structOut%carbon%owner_model_id      = struct1%carbon%owner_model_id
    structOut%carbon%owner_proc_id       = struct1%carbon%owner_proc_id
    structOut%carbon%owner_tile_path     = struct1%carbon%owner_tile_path
    ALLOCATE(structOut%carbon%ptr2d,     source=struct1%carbon%ptr2d)
    structOut%nitrogen%owner_model_id    = struct1%nitrogen%owner_model_id
    structOut%nitrogen%owner_proc_id     = struct1%nitrogen%owner_proc_id
    structOut%nitrogen%owner_tile_path   = struct1%nitrogen%owner_tile_path
    ALLOCATE(structOut%nitrogen%ptr2d,   source=struct1%nitrogen%ptr2d)
    structOut%phosphorus%owner_model_id  = struct1%phosphorus%owner_model_id
    structOut%phosphorus%owner_proc_id   = struct1%phosphorus%owner_proc_id
    structOut%phosphorus%owner_tile_path = struct1%phosphorus%owner_tile_path
    ALLOCATE(structOut%phosphorus%ptr2d, source=struct1%phosphorus%ptr2d)
    structOut%carbon13%owner_model_id    = struct1%carbon13%owner_model_id
    structOut%carbon13%owner_proc_id     = struct1%carbon13%owner_proc_id
    structOut%carbon13%owner_tile_path   = struct1%carbon13%owner_tile_path
    ALLOCATE(structOut%carbon13%ptr2d,   source=struct1%carbon13%ptr2d)
    structOut%carbon14%owner_model_id    = struct1%carbon14%owner_model_id
    structOut%carbon14%owner_proc_id     = struct1%carbon14%owner_proc_id
    structOut%carbon14%owner_tile_path   = struct1%carbon14%owner_tile_path
    ALLOCATE(structOut%carbon14%ptr2d,   source=struct1%carbon14%ptr2d)
    structOut%nitrogen15%owner_model_id  = struct1%nitrogen15%owner_model_id
    structOut%nitrogen15%owner_proc_id   = struct1%nitrogen15%owner_proc_id
    structOut%nitrogen15%owner_tile_path = struct1%nitrogen15%owner_tile_path
    ALLOCATE(structOut%nitrogen15%ptr2d, source=struct1%nitrogen15%ptr2d)

    ! run calculation
    ! the differentiation between ON_CHUNK/ON_DOMAIN is done at the level of variables (i.e., the op_ov routines of the jsb_var)
    structOut%carbon          = struct1%carbon      / field1d(:)
    structOut%nitrogen        = struct1%nitrogen    / field1d(:)
    structOut%phosphorus      = struct1%phosphorus  / field1d(:)
    structOut%carbon13        = struct1%carbon13    / field1d(:)
    structOut%carbon14        = struct1%carbon14    / field1d(:)
    structOut%nitrogen15      = struct1%nitrogen15  / field1d(:)

    ! is this needed here too (not only with the jsb_var op_ov)
    structOut%carbon%ptr      => structOut%carbon%ptr2d(:,:)
    structOut%nitrogen%ptr    => structOut%nitrogen%ptr2d(:,:)
    structOut%phosphorus%ptr  => structOut%phosphorus%ptr2d(:,:)
    structOut%carbon13%ptr    => structOut%carbon13%ptr2d(:,:)
    structOut%carbon14%ptr    => structOut%carbon14%ptr2d(:,:)
    structOut%nitrogen15%ptr  => structOut%nitrogen15%ptr2d(:,:)
  
  END FUNCTION Div_veg_elem_p_by_field1d
  !--------------------------------------------------------------------------------------------------------------------------------
  ! t_veg_elements_pool / field2d(:,:)        | ON_CHUNK & ON_DOMAIN | field2d(:,1) & field2d(:,:)
  !--------------------------------------------------------------------------------------------------------------------------------
  FUNCTION Div_veg_elem_p_by_field2d(struct1, field2d) RESULT(structOut)
  
    CLASS(t_veg_elements_pool),     INTENT(in)  :: struct1
    REAL(wp),                       INTENT(in)  :: field2d(:,:)
    TYPE(t_veg_elements_pool),      TARGET      :: structOut
  
    CHARACTER(len=*), PARAMETER :: routine = modname//':Div_veg_elem_p_by_field2d'
  
    ! init temporary pool structOut
    CALL structOut%Init(name='op_ov_help', shortname='op_ov_help', l_elements=.TRUE., element_unit=TRIM(struct1%element_unit))
    structOut%carbon%owner_model_id      = struct1%carbon%owner_model_id
    structOut%carbon%owner_proc_id       = struct1%carbon%owner_proc_id
    structOut%carbon%owner_tile_path     = struct1%carbon%owner_tile_path
    ALLOCATE(structOut%carbon%ptr2d,     source=struct1%carbon%ptr2d)
    structOut%nitrogen%owner_model_id    = struct1%nitrogen%owner_model_id
    structOut%nitrogen%owner_proc_id     = struct1%nitrogen%owner_proc_id
    structOut%nitrogen%owner_tile_path   = struct1%nitrogen%owner_tile_path
    ALLOCATE(structOut%nitrogen%ptr2d,   source=struct1%nitrogen%ptr2d)
    structOut%phosphorus%owner_model_id  = struct1%phosphorus%owner_model_id
    structOut%phosphorus%owner_proc_id   = struct1%phosphorus%owner_proc_id
    structOut%phosphorus%owner_tile_path = struct1%phosphorus%owner_tile_path
    ALLOCATE(structOut%phosphorus%ptr2d, source=struct1%phosphorus%ptr2d)
    structOut%carbon13%owner_model_id    = struct1%carbon13%owner_model_id
    structOut%carbon13%owner_proc_id     = struct1%carbon13%owner_proc_id
    structOut%carbon13%owner_tile_path   = struct1%carbon13%owner_tile_path
    ALLOCATE(structOut%carbon13%ptr2d,   source=struct1%carbon13%ptr2d)
    structOut%carbon14%owner_model_id    = struct1%carbon14%owner_model_id
    structOut%carbon14%owner_proc_id     = struct1%carbon14%owner_proc_id
    structOut%carbon14%owner_tile_path   = struct1%carbon14%owner_tile_path
    ALLOCATE(structOut%carbon14%ptr2d,   source=struct1%carbon14%ptr2d)
    structOut%nitrogen15%owner_model_id  = struct1%nitrogen15%owner_model_id
    structOut%nitrogen15%owner_proc_id   = struct1%nitrogen15%owner_proc_id
    structOut%nitrogen15%owner_tile_path = struct1%nitrogen15%owner_tile_path
    ALLOCATE(structOut%nitrogen15%ptr2d, source=struct1%nitrogen15%ptr2d)

    ! run calculation
    ! the differentiation between ON_CHUNK/ON_DOMAIN is done at the level of variables (i.e., the op_ov routines of the jsb_var)
    structOut%carbon          = struct1%carbon      / field2d(:,:)
    structOut%nitrogen        = struct1%nitrogen    / field2d(:,:)
    structOut%phosphorus      = struct1%phosphorus  / field2d(:,:)
    structOut%carbon13        = struct1%carbon13    / field2d(:,:)
    structOut%carbon14        = struct1%carbon14    / field2d(:,:)
    structOut%nitrogen15      = struct1%nitrogen15  / field2d(:,:)

    ! is this needed here too (not only with the jsb_var op_ov)
    structOut%carbon%ptr      => structOut%carbon%ptr2d(:,:)
    structOut%nitrogen%ptr    => structOut%nitrogen%ptr2d(:,:)
    structOut%phosphorus%ptr  => structOut%phosphorus%ptr2d(:,:)
    structOut%carbon13%ptr    => structOut%carbon13%ptr2d(:,:)
    structOut%carbon14%ptr    => structOut%carbon14%ptr2d(:,:)
    structOut%nitrogen15%ptr  => structOut%nitrogen15%ptr2d(:,:)
  
  END FUNCTION Div_veg_elem_p_by_field2d
  !--------------------------------------------------------------------------------------------------------------------------------
  ! t_veg_elements_pool / t_veg_elements_pool | ON_CHUNK & ON_DOMAIN
  !--------------------------------------------------------------------------------------------------------------------------------
  FUNCTION Div_veg_elem_p_by_veg_elem_p(struct1, struct2) RESULT(structOut)
  
    CLASS(t_veg_elements_pool),     INTENT(in)  :: struct1, struct2
    TYPE(t_veg_elements_pool),      TARGET      :: structOut
  
    CHARACTER(len=*), PARAMETER :: routine = modname//':Div_veg_elem_p_by_veg_elem_p'
  
    ! init temporary pool structOut
    CALL structOut%Init(name='op_ov_help', shortname='op_ov_help', l_elements=.TRUE., element_unit=TRIM(struct1%element_unit))
    structOut%carbon%owner_model_id      = struct1%carbon%owner_model_id
    structOut%carbon%owner_proc_id       = struct1%carbon%owner_proc_id
    structOut%carbon%owner_tile_path     = struct1%carbon%owner_tile_path
    ALLOCATE(structOut%carbon%ptr2d,     source=struct1%carbon%ptr2d)
    structOut%nitrogen%owner_model_id    = struct1%nitrogen%owner_model_id
    structOut%nitrogen%owner_proc_id     = struct1%nitrogen%owner_proc_id
    structOut%nitrogen%owner_tile_path   = struct1%nitrogen%owner_tile_path
    ALLOCATE(structOut%nitrogen%ptr2d,   source=struct1%nitrogen%ptr2d)
    structOut%phosphorus%owner_model_id  = struct1%phosphorus%owner_model_id
    structOut%phosphorus%owner_proc_id   = struct1%phosphorus%owner_proc_id
    structOut%phosphorus%owner_tile_path = struct1%phosphorus%owner_tile_path
    ALLOCATE(structOut%phosphorus%ptr2d, source=struct1%phosphorus%ptr2d)
    structOut%carbon13%owner_model_id    = struct1%carbon13%owner_model_id
    structOut%carbon13%owner_proc_id     = struct1%carbon13%owner_proc_id
    structOut%carbon13%owner_tile_path   = struct1%carbon13%owner_tile_path
    ALLOCATE(structOut%carbon13%ptr2d,   source=struct1%carbon13%ptr2d)
    structOut%carbon14%owner_model_id    = struct1%carbon14%owner_model_id
    structOut%carbon14%owner_proc_id     = struct1%carbon14%owner_proc_id
    structOut%carbon14%owner_tile_path   = struct1%carbon14%owner_tile_path
    ALLOCATE(structOut%carbon14%ptr2d,   source=struct1%carbon14%ptr2d)
    structOut%nitrogen15%owner_model_id  = struct1%nitrogen15%owner_model_id
    structOut%nitrogen15%owner_proc_id   = struct1%nitrogen15%owner_proc_id
    structOut%nitrogen15%owner_tile_path = struct1%nitrogen15%owner_tile_path
    ALLOCATE(structOut%nitrogen15%ptr2d, source=struct1%nitrogen15%ptr2d)

    ! run calculation
    ! the differentiation between ON_CHUNK/ON_DOMAIN is done at the level of variables (i.e., the op_ov routines of the jsb_var)
    structOut%carbon          = struct1%carbon      / struct2%carbon
    structOut%nitrogen        = struct1%nitrogen    / struct2%nitrogen  
    structOut%phosphorus      = struct1%phosphorus  / struct2%phosphorus
    structOut%carbon13        = struct1%carbon13    / struct2%carbon13  
    structOut%carbon14        = struct1%carbon14    / struct2%carbon14  
    structOut%nitrogen15      = struct1%nitrogen15  / struct2%nitrogen15

    ! is this needed here too (not only with the jsb_var op_ov)
    structOut%carbon%ptr      => structOut%carbon%ptr2d(:,:)
    structOut%nitrogen%ptr    => structOut%nitrogen%ptr2d(:,:)
    structOut%phosphorus%ptr  => structOut%phosphorus%ptr2d(:,:)
    structOut%carbon13%ptr    => structOut%carbon13%ptr2d(:,:)
    structOut%carbon14%ptr    => structOut%carbon14%ptr2d(:,:)
    structOut%nitrogen15%ptr  => structOut%nitrogen15%ptr2d(:,:)
  
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

    ! not allocated: this%path
    ! no need to dealloc the element_var & element_pools because these are dealloc by their own FINAL routine
    IF (ASSOCIATED(this%carbon))        DEALLOCATE(this%carbon)
    IF (ASSOCIATED(this%nitrogen))      DEALLOCATE(this%nitrogen)
    IF (ASSOCIATED(this%phosphorus))    DEALLOCATE(this%phosphorus)
    IF (ASSOCIATED(this%carbon13))      DEALLOCATE(this%carbon13)
    IF (ASSOCIATED(this%carbon14))      DEALLOCATE(this%carbon14)
    IF (ASSOCIATED(this%nitrogen15))    DEALLOCATE(this%nitrogen15)
    IF (ALLOCATED(this%shortname))      DEALLOCATE(this%shortname)
    IF (ALLOCATED(this%name))           DEALLOCATE(this%name)
    IF (ALLOCATED(this%element_unit))   DEALLOCATE(this%element_unit)
    IF (ALLOCATED(this%element_list))   DEALLOCATE(this%element_list)
    IF (ALLOCATED(this%var_list))       DEALLOCATE(this%var_list)
    IF (ASSOCIATED(this%pool_list))     DEALLOCATE(this%pool_list)
    IF (ASSOCIATED(this%leaf_list))     DEALLOCATE(this%leaf_list)
    
    IF (ASSOCIATED(this%leaf))          DEALLOCATE(this%leaf)
    IF (ASSOCIATED(this%fine_root))     DEALLOCATE(this%fine_root)
    IF (ASSOCIATED(this%coarse_root))   DEALLOCATE(this%coarse_root)
    IF (ASSOCIATED(this%sap_wood))      DEALLOCATE(this%sap_wood)
    IF (ASSOCIATED(this%heart_wood))    DEALLOCATE(this%heart_wood)
    IF (ASSOCIATED(this%labile))        DEALLOCATE(this%labile)
    IF (ASSOCIATED(this%reserve))       DEALLOCATE(this%reserve)
    IF (ASSOCIATED(this%fruit))         DEALLOCATE(this%fruit)
    IF (ASSOCIATED(this%seed_bed))      DEALLOCATE(this%seed_bed)
   
  END SUBROUTINE finalize_veg_tissues_pool

  !--------------------------------------------------------------------------------------------------------------------------------
  ! = assign
  !--------------------------------------------------------------------------------------------------------------------------------
  !--------------------------------------------------------------------------------------------------------------------------------
  ! veg_tissues_pool = scalar | ON_CHUNK & ON_DOMAIN
  !--------------------------------------------------------------------------------------------------------------------------------
  RECURSIVE SUBROUTINE Assign_scalar_to_veg_tiss_p(lhs, rhs)

    CLASS(t_veg_tissues_pool),  INTENT(inout) :: lhs
    REAL(wp),                   INTENT(in)    :: rhs

    CHARACTER(len=*), PARAMETER :: routine = modname//':Assign_scalar_to_veg_tiss_p'

    lhs%carbon      = rhs
    lhs%nitrogen    = rhs
    lhs%phosphorus  = rhs
    lhs%carbon13    = rhs
    lhs%carbon14    = rhs
    lhs%nitrogen15  = rhs
    
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

    CHARACTER(len=*), PARAMETER :: routine = modname//':Assign_field1d_to_veg_tiss_p'

    lhs%carbon      = rhs(:)
    lhs%nitrogen    = rhs(:)
    lhs%phosphorus  = rhs(:)
    lhs%carbon13    = rhs(:)
    lhs%carbon14    = rhs(:)
    lhs%nitrogen15  = rhs(:)

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

    CHARACTER(len=*), PARAMETER :: routine = modname//':Assign_field2d_to_veg_tiss_p'

    lhs%carbon      = rhs(:,:)
    lhs%nitrogen    = rhs(:,:)
    lhs%phosphorus  = rhs(:,:)
    lhs%carbon13    = rhs(:,:)
    lhs%carbon14    = rhs(:,:)
    lhs%nitrogen15  = rhs(:,:)

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

    CHARACTER(len=*), PARAMETER :: routine = modname//':Assign_veg_tiss_p_to_veg_tiss_p'

    lhs%carbon      = rhs%carbon    
    lhs%nitrogen    = rhs%nitrogen  
    lhs%phosphorus  = rhs%phosphorus
    lhs%carbon13    = rhs%carbon13  
    lhs%carbon14    = rhs%carbon14  
    lhs%nitrogen15  = rhs%nitrogen15
    
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
    TYPE(t_veg_tissues_pool),       TARGET      :: structOut
  
    CHARACTER(len=*), PARAMETER :: routine = modname//':Add_veg_tiss_p_and_scalar'
  
    !> 1.0 init temporary pool structOut
    !! init tissue pool with element variables
    CALL structOut%Init_tiss_pool_for_op_ov(name='op_ov_help', shortname='op_ov_help', &
      & unit=TRIM(struct1%element_unit), pool_suffix='veg_p') ! l_elements=.TRUE. by default
    structOut%carbon%owner_model_id      = struct1%carbon%owner_model_id
    structOut%carbon%owner_proc_id       = struct1%carbon%owner_proc_id
    structOut%carbon%owner_tile_path     = struct1%carbon%owner_tile_path
    ALLOCATE(structOut%carbon%ptr2d,     source=struct1%carbon%ptr2d)
    structOut%nitrogen%owner_model_id    = struct1%nitrogen%owner_model_id
    structOut%nitrogen%owner_proc_id     = struct1%nitrogen%owner_proc_id
    structOut%nitrogen%owner_tile_path   = struct1%nitrogen%owner_tile_path
    ALLOCATE(structOut%nitrogen%ptr2d,   source=struct1%nitrogen%ptr2d)
    structOut%phosphorus%owner_model_id  = struct1%phosphorus%owner_model_id
    structOut%phosphorus%owner_proc_id   = struct1%phosphorus%owner_proc_id
    structOut%phosphorus%owner_tile_path = struct1%phosphorus%owner_tile_path
    ALLOCATE(structOut%phosphorus%ptr2d, source=struct1%phosphorus%ptr2d)
    structOut%carbon13%owner_model_id    = struct1%carbon13%owner_model_id
    structOut%carbon13%owner_proc_id     = struct1%carbon13%owner_proc_id
    structOut%carbon13%owner_tile_path   = struct1%carbon13%owner_tile_path
    ALLOCATE(structOut%carbon13%ptr2d,   source=struct1%carbon13%ptr2d)
    structOut%carbon14%owner_model_id    = struct1%carbon14%owner_model_id
    structOut%carbon14%owner_proc_id     = struct1%carbon14%owner_proc_id
    structOut%carbon14%owner_tile_path   = struct1%carbon14%owner_tile_path
    ALLOCATE(structOut%carbon14%ptr2d,   source=struct1%carbon14%ptr2d)
    structOut%nitrogen15%owner_model_id  = struct1%nitrogen15%owner_model_id
    structOut%nitrogen15%owner_proc_id   = struct1%nitrogen15%owner_proc_id
    structOut%nitrogen15%owner_tile_path = struct1%nitrogen15%owner_tile_path
    ALLOCATE(structOut%nitrogen15%ptr2d, source=struct1%nitrogen15%ptr2d)
    !! init element pools
    CALL structOut%leaf%Init_elem_var_elem_p(struct1%leaf)
    CALL structOut%fine_root%Init_elem_var_elem_p(struct1%fine_root)
    CALL structOut%coarse_root%Init_elem_var_elem_p(struct1%coarse_root)
    CALL structOut%sap_wood%Init_elem_var_elem_p(struct1%sap_wood)
    CALL structOut%heart_wood%Init_elem_var_elem_p(struct1%heart_wood)
    CALL structOut%labile%Init_elem_var_elem_p(struct1%labile)
    CALL structOut%reserve%Init_elem_var_elem_p(struct1%reserve)
    CALL structOut%fruit%Init_elem_var_elem_p(struct1%fruit)
    CALL structOut%seed_bed%Init_elem_var_elem_p(struct1%seed_bed)
  
    !> 2.0 run calculation
    !! the differentiation between ON_CHUNK/ON_DOMAIN is done at the level of variables (i.e., the op_ov routines of the jsb_var)
    !> 2.1 element variables of the tissue pool itself
    structOut%carbon          = struct1%carbon      + scalar
    structOut%nitrogen        = struct1%nitrogen    + scalar
    structOut%phosphorus      = struct1%phosphorus  + scalar
    structOut%carbon13        = struct1%carbon13    + scalar
    structOut%carbon14        = struct1%carbon14    + scalar
    structOut%nitrogen15      = struct1%nitrogen15  + scalar
    ! adding the pointer is this needed here too (not only with the jsb_var op_ov)
    structOut%carbon%ptr      => structOut%carbon%ptr2d(:,:)
    structOut%nitrogen%ptr    => structOut%nitrogen%ptr2d(:,:)
    structOut%phosphorus%ptr  => structOut%phosphorus%ptr2d(:,:)
    structOut%carbon13%ptr    => structOut%carbon13%ptr2d(:,:)
    structOut%carbon14%ptr    => structOut%carbon14%ptr2d(:,:)
    structOut%nitrogen15%ptr  => structOut%nitrogen15%ptr2d(:,:)
    !> 2.2 element pools
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
    TYPE(t_veg_tissues_pool),       TARGET      :: structOut
  
    CHARACTER(len=*), PARAMETER :: routine = modname//':Add_scalar_and_veg_tiss_p'
  
    !> 1.0 init temporary pool structOut
    !! init tissue pool with element variables
    CALL structOut%Init_tiss_pool_for_op_ov(name='op_ov_help', shortname='op_ov_help', &
      & unit=TRIM(struct1%element_unit), pool_suffix='veg_p') ! l_elements=.TRUE. by default
    structOut%carbon%owner_model_id      = struct1%carbon%owner_model_id
    structOut%carbon%owner_proc_id       = struct1%carbon%owner_proc_id
    structOut%carbon%owner_tile_path     = struct1%carbon%owner_tile_path
    ALLOCATE(structOut%carbon%ptr2d,     source=struct1%carbon%ptr2d)
    structOut%nitrogen%owner_model_id    = struct1%nitrogen%owner_model_id
    structOut%nitrogen%owner_proc_id     = struct1%nitrogen%owner_proc_id
    structOut%nitrogen%owner_tile_path   = struct1%nitrogen%owner_tile_path
    ALLOCATE(structOut%nitrogen%ptr2d,   source=struct1%nitrogen%ptr2d)
    structOut%phosphorus%owner_model_id  = struct1%phosphorus%owner_model_id
    structOut%phosphorus%owner_proc_id   = struct1%phosphorus%owner_proc_id
    structOut%phosphorus%owner_tile_path = struct1%phosphorus%owner_tile_path
    ALLOCATE(structOut%phosphorus%ptr2d, source=struct1%phosphorus%ptr2d)
    structOut%carbon13%owner_model_id    = struct1%carbon13%owner_model_id
    structOut%carbon13%owner_proc_id     = struct1%carbon13%owner_proc_id
    structOut%carbon13%owner_tile_path   = struct1%carbon13%owner_tile_path
    ALLOCATE(structOut%carbon13%ptr2d,   source=struct1%carbon13%ptr2d)
    structOut%carbon14%owner_model_id    = struct1%carbon14%owner_model_id
    structOut%carbon14%owner_proc_id     = struct1%carbon14%owner_proc_id
    structOut%carbon14%owner_tile_path   = struct1%carbon14%owner_tile_path
    ALLOCATE(structOut%carbon14%ptr2d,   source=struct1%carbon14%ptr2d)
    structOut%nitrogen15%owner_model_id  = struct1%nitrogen15%owner_model_id
    structOut%nitrogen15%owner_proc_id   = struct1%nitrogen15%owner_proc_id
    structOut%nitrogen15%owner_tile_path = struct1%nitrogen15%owner_tile_path
    ALLOCATE(structOut%nitrogen15%ptr2d, source=struct1%nitrogen15%ptr2d)
    !! init element pools
    CALL structOut%leaf%Init_elem_var_elem_p(struct1%leaf)
    CALL structOut%fine_root%Init_elem_var_elem_p(struct1%fine_root)
    CALL structOut%coarse_root%Init_elem_var_elem_p(struct1%coarse_root)
    CALL structOut%sap_wood%Init_elem_var_elem_p(struct1%sap_wood)
    CALL structOut%heart_wood%Init_elem_var_elem_p(struct1%heart_wood)
    CALL structOut%labile%Init_elem_var_elem_p(struct1%labile)
    CALL structOut%reserve%Init_elem_var_elem_p(struct1%reserve)
    CALL structOut%fruit%Init_elem_var_elem_p(struct1%fruit)
    CALL structOut%seed_bed%Init_elem_var_elem_p(struct1%seed_bed)
  
    !> 2.0 run calculation
    !! the differentiation between ON_CHUNK/ON_DOMAIN is done at the level of variables (i.e., the op_ov routines of the jsb_var)
    !> 2.1 element variables of the tissue pool itself
    structOut%carbon          = scalar + struct1%carbon    
    structOut%nitrogen        = scalar + struct1%nitrogen  
    structOut%phosphorus      = scalar + struct1%phosphorus
    structOut%carbon13        = scalar + struct1%carbon13  
    structOut%carbon14        = scalar + struct1%carbon14  
    structOut%nitrogen15      = scalar + struct1%nitrogen15
    ! adding the pointer is this needed here too (not only with the jsb_var op_ov)
    structOut%carbon%ptr      => structOut%carbon%ptr2d(:,:)
    structOut%nitrogen%ptr    => structOut%nitrogen%ptr2d(:,:)
    structOut%phosphorus%ptr  => structOut%phosphorus%ptr2d(:,:)
    structOut%carbon13%ptr    => structOut%carbon13%ptr2d(:,:)
    structOut%carbon14%ptr    => structOut%carbon14%ptr2d(:,:)
    structOut%nitrogen15%ptr  => structOut%nitrogen15%ptr2d(:,:)
    !> 2.2 element pools
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
    TYPE(t_veg_tissues_pool),       TARGET      :: structOut
  
    CHARACTER(len=*), PARAMETER :: routine = modname//':Add_veg_tiss_p_and_field1d'
  
    !> 1.0 init temporary pool structOut
    !! init tissue pool with element variables
    CALL structOut%Init_tiss_pool_for_op_ov(name='op_ov_help', shortname='op_ov_help', &
      & unit=TRIM(struct1%element_unit), pool_suffix='veg_p') ! l_elements=.TRUE. by default
    structOut%carbon%owner_model_id      = struct1%carbon%owner_model_id
    structOut%carbon%owner_proc_id       = struct1%carbon%owner_proc_id
    structOut%carbon%owner_tile_path     = struct1%carbon%owner_tile_path
    ALLOCATE(structOut%carbon%ptr2d,     source=struct1%carbon%ptr2d)
    structOut%nitrogen%owner_model_id    = struct1%nitrogen%owner_model_id
    structOut%nitrogen%owner_proc_id     = struct1%nitrogen%owner_proc_id
    structOut%nitrogen%owner_tile_path   = struct1%nitrogen%owner_tile_path
    ALLOCATE(structOut%nitrogen%ptr2d,   source=struct1%nitrogen%ptr2d)
    structOut%phosphorus%owner_model_id  = struct1%phosphorus%owner_model_id
    structOut%phosphorus%owner_proc_id   = struct1%phosphorus%owner_proc_id
    structOut%phosphorus%owner_tile_path = struct1%phosphorus%owner_tile_path
    ALLOCATE(structOut%phosphorus%ptr2d, source=struct1%phosphorus%ptr2d)
    structOut%carbon13%owner_model_id    = struct1%carbon13%owner_model_id
    structOut%carbon13%owner_proc_id     = struct1%carbon13%owner_proc_id
    structOut%carbon13%owner_tile_path   = struct1%carbon13%owner_tile_path
    ALLOCATE(structOut%carbon13%ptr2d,   source=struct1%carbon13%ptr2d)
    structOut%carbon14%owner_model_id    = struct1%carbon14%owner_model_id
    structOut%carbon14%owner_proc_id     = struct1%carbon14%owner_proc_id
    structOut%carbon14%owner_tile_path   = struct1%carbon14%owner_tile_path
    ALLOCATE(structOut%carbon14%ptr2d,   source=struct1%carbon14%ptr2d)
    structOut%nitrogen15%owner_model_id  = struct1%nitrogen15%owner_model_id
    structOut%nitrogen15%owner_proc_id   = struct1%nitrogen15%owner_proc_id
    structOut%nitrogen15%owner_tile_path = struct1%nitrogen15%owner_tile_path
    ALLOCATE(structOut%nitrogen15%ptr2d, source=struct1%nitrogen15%ptr2d)
    !! init element pools
    CALL structOut%leaf%Init_elem_var_elem_p(struct1%leaf)
    CALL structOut%fine_root%Init_elem_var_elem_p(struct1%fine_root)
    CALL structOut%coarse_root%Init_elem_var_elem_p(struct1%coarse_root)
    CALL structOut%sap_wood%Init_elem_var_elem_p(struct1%sap_wood)
    CALL structOut%heart_wood%Init_elem_var_elem_p(struct1%heart_wood)
    CALL structOut%labile%Init_elem_var_elem_p(struct1%labile)
    CALL structOut%reserve%Init_elem_var_elem_p(struct1%reserve)
    CALL structOut%fruit%Init_elem_var_elem_p(struct1%fruit)
    CALL structOut%seed_bed%Init_elem_var_elem_p(struct1%seed_bed)
  
    !> 2.0 run calculation
    !! the differentiation between ON_CHUNK/ON_DOMAIN is done at the level of variables (i.e., the op_ov routines of the jsb_var)
    !> 2.1 element variables of the tissue pool itself
    structOut%carbon          = struct1%carbon      + field1d(:)
    structOut%nitrogen        = struct1%nitrogen    + field1d(:)
    structOut%phosphorus      = struct1%phosphorus  + field1d(:)
    structOut%carbon13        = struct1%carbon13    + field1d(:)
    structOut%carbon14        = struct1%carbon14    + field1d(:)
    structOut%nitrogen15      = struct1%nitrogen15  + field1d(:)
    ! adding the pointer is this needed here too (not only with the jsb_var op_ov)
    structOut%carbon%ptr      => structOut%carbon%ptr2d(:,:)
    structOut%nitrogen%ptr    => structOut%nitrogen%ptr2d(:,:)
    structOut%phosphorus%ptr  => structOut%phosphorus%ptr2d(:,:)
    structOut%carbon13%ptr    => structOut%carbon13%ptr2d(:,:)
    structOut%carbon14%ptr    => structOut%carbon14%ptr2d(:,:)
    structOut%nitrogen15%ptr  => structOut%nitrogen15%ptr2d(:,:)
    !> 2.2 element pools
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
    TYPE(t_veg_tissues_pool),       TARGET      :: structOut
  
    CHARACTER(len=*), PARAMETER :: routine = modname//':Add_field1d_and_veg_tiss_p'
  
    !> 1.0 init temporary pool structOut
    !! init tissue pool with element variables
    CALL structOut%Init_tiss_pool_for_op_ov(name='op_ov_help', shortname='op_ov_help', &
      & unit=TRIM(struct1%element_unit), pool_suffix='veg_p') ! l_elements=.TRUE. by default
    structOut%carbon%owner_model_id      = struct1%carbon%owner_model_id
    structOut%carbon%owner_proc_id       = struct1%carbon%owner_proc_id
    structOut%carbon%owner_tile_path     = struct1%carbon%owner_tile_path
    ALLOCATE(structOut%carbon%ptr2d,     source=struct1%carbon%ptr2d)
    structOut%nitrogen%owner_model_id    = struct1%nitrogen%owner_model_id
    structOut%nitrogen%owner_proc_id     = struct1%nitrogen%owner_proc_id
    structOut%nitrogen%owner_tile_path   = struct1%nitrogen%owner_tile_path
    ALLOCATE(structOut%nitrogen%ptr2d,   source=struct1%nitrogen%ptr2d)
    structOut%phosphorus%owner_model_id  = struct1%phosphorus%owner_model_id
    structOut%phosphorus%owner_proc_id   = struct1%phosphorus%owner_proc_id
    structOut%phosphorus%owner_tile_path = struct1%phosphorus%owner_tile_path
    ALLOCATE(structOut%phosphorus%ptr2d, source=struct1%phosphorus%ptr2d)
    structOut%carbon13%owner_model_id    = struct1%carbon13%owner_model_id
    structOut%carbon13%owner_proc_id     = struct1%carbon13%owner_proc_id
    structOut%carbon13%owner_tile_path   = struct1%carbon13%owner_tile_path
    ALLOCATE(structOut%carbon13%ptr2d,   source=struct1%carbon13%ptr2d)
    structOut%carbon14%owner_model_id    = struct1%carbon14%owner_model_id
    structOut%carbon14%owner_proc_id     = struct1%carbon14%owner_proc_id
    structOut%carbon14%owner_tile_path   = struct1%carbon14%owner_tile_path
    ALLOCATE(structOut%carbon14%ptr2d,   source=struct1%carbon14%ptr2d)
    structOut%nitrogen15%owner_model_id  = struct1%nitrogen15%owner_model_id
    structOut%nitrogen15%owner_proc_id   = struct1%nitrogen15%owner_proc_id
    structOut%nitrogen15%owner_tile_path = struct1%nitrogen15%owner_tile_path
    ALLOCATE(structOut%nitrogen15%ptr2d, source=struct1%nitrogen15%ptr2d)
    !! init element pools
    CALL structOut%leaf%Init_elem_var_elem_p(struct1%leaf)
    CALL structOut%fine_root%Init_elem_var_elem_p(struct1%fine_root)
    CALL structOut%coarse_root%Init_elem_var_elem_p(struct1%coarse_root)
    CALL structOut%sap_wood%Init_elem_var_elem_p(struct1%sap_wood)
    CALL structOut%heart_wood%Init_elem_var_elem_p(struct1%heart_wood)
    CALL structOut%labile%Init_elem_var_elem_p(struct1%labile)
    CALL structOut%reserve%Init_elem_var_elem_p(struct1%reserve)
    CALL structOut%fruit%Init_elem_var_elem_p(struct1%fruit)
    CALL structOut%seed_bed%Init_elem_var_elem_p(struct1%seed_bed)
  
    !> 2.0 run calculation
    !! the differentiation between ON_CHUNK/ON_DOMAIN is done at the level of variables (i.e., the op_ov routines of the jsb_var)
    !> 2.1 element variables of the tissue pool itself
    structOut%carbon          = field1d(:) + struct1%carbon    
    structOut%nitrogen        = field1d(:) + struct1%nitrogen  
    structOut%phosphorus      = field1d(:) + struct1%phosphorus
    structOut%carbon13        = field1d(:) + struct1%carbon13  
    structOut%carbon14        = field1d(:) + struct1%carbon14  
    structOut%nitrogen15      = field1d(:) + struct1%nitrogen15
    ! adding the pointer is this needed here too (not only with the jsb_var op_ov)
    structOut%carbon%ptr      => structOut%carbon%ptr2d(:,:)
    structOut%nitrogen%ptr    => structOut%nitrogen%ptr2d(:,:)
    structOut%phosphorus%ptr  => structOut%phosphorus%ptr2d(:,:)
    structOut%carbon13%ptr    => structOut%carbon13%ptr2d(:,:)
    structOut%carbon14%ptr    => structOut%carbon14%ptr2d(:,:)
    structOut%nitrogen15%ptr  => structOut%nitrogen15%ptr2d(:,:)
    !> 2.2 element pools
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
    TYPE(t_veg_tissues_pool),       TARGET      :: structOut
  
    CHARACTER(len=*), PARAMETER :: routine = modname//':Add_veg_tiss_p_and_field2d'
  
    !> 1.0 init temporary pool structOut
    !! init tissue pool with element variables
    CALL structOut%Init_tiss_pool_for_op_ov(name='op_ov_help', shortname='op_ov_help', &
      & unit=TRIM(struct1%element_unit), pool_suffix='veg_p') ! l_elements=.TRUE. by default
    structOut%carbon%owner_model_id      = struct1%carbon%owner_model_id
    structOut%carbon%owner_proc_id       = struct1%carbon%owner_proc_id
    structOut%carbon%owner_tile_path     = struct1%carbon%owner_tile_path
    ALLOCATE(structOut%carbon%ptr2d,     source=struct1%carbon%ptr2d)
    structOut%nitrogen%owner_model_id    = struct1%nitrogen%owner_model_id
    structOut%nitrogen%owner_proc_id     = struct1%nitrogen%owner_proc_id
    structOut%nitrogen%owner_tile_path   = struct1%nitrogen%owner_tile_path
    ALLOCATE(structOut%nitrogen%ptr2d,   source=struct1%nitrogen%ptr2d)
    structOut%phosphorus%owner_model_id  = struct1%phosphorus%owner_model_id
    structOut%phosphorus%owner_proc_id   = struct1%phosphorus%owner_proc_id
    structOut%phosphorus%owner_tile_path = struct1%phosphorus%owner_tile_path
    ALLOCATE(structOut%phosphorus%ptr2d, source=struct1%phosphorus%ptr2d)
    structOut%carbon13%owner_model_id    = struct1%carbon13%owner_model_id
    structOut%carbon13%owner_proc_id     = struct1%carbon13%owner_proc_id
    structOut%carbon13%owner_tile_path   = struct1%carbon13%owner_tile_path
    ALLOCATE(structOut%carbon13%ptr2d,   source=struct1%carbon13%ptr2d)
    structOut%carbon14%owner_model_id    = struct1%carbon14%owner_model_id
    structOut%carbon14%owner_proc_id     = struct1%carbon14%owner_proc_id
    structOut%carbon14%owner_tile_path   = struct1%carbon14%owner_tile_path
    ALLOCATE(structOut%carbon14%ptr2d,   source=struct1%carbon14%ptr2d)
    structOut%nitrogen15%owner_model_id  = struct1%nitrogen15%owner_model_id
    structOut%nitrogen15%owner_proc_id   = struct1%nitrogen15%owner_proc_id
    structOut%nitrogen15%owner_tile_path = struct1%nitrogen15%owner_tile_path
    ALLOCATE(structOut%nitrogen15%ptr2d, source=struct1%nitrogen15%ptr2d)
    !! init element pools
    CALL structOut%leaf%Init_elem_var_elem_p(struct1%leaf)
    CALL structOut%fine_root%Init_elem_var_elem_p(struct1%fine_root)
    CALL structOut%coarse_root%Init_elem_var_elem_p(struct1%coarse_root)
    CALL structOut%sap_wood%Init_elem_var_elem_p(struct1%sap_wood)
    CALL structOut%heart_wood%Init_elem_var_elem_p(struct1%heart_wood)
    CALL structOut%labile%Init_elem_var_elem_p(struct1%labile)
    CALL structOut%reserve%Init_elem_var_elem_p(struct1%reserve)
    CALL structOut%fruit%Init_elem_var_elem_p(struct1%fruit)
    CALL structOut%seed_bed%Init_elem_var_elem_p(struct1%seed_bed)
  
    !> 2.0 run calculation
    !! the differentiation between ON_CHUNK/ON_DOMAIN is done at the level of variables (i.e., the op_ov routines of the jsb_var)
    !> 2.1 element variables of the tissue pool itself
    structOut%carbon          = struct1%carbon      + field2d(:,:)
    structOut%nitrogen        = struct1%nitrogen    + field2d(:,:)
    structOut%phosphorus      = struct1%phosphorus  + field2d(:,:)
    structOut%carbon13        = struct1%carbon13    + field2d(:,:)
    structOut%carbon14        = struct1%carbon14    + field2d(:,:)
    structOut%nitrogen15      = struct1%nitrogen15  + field2d(:,:)
    ! adding the pointer is this needed here too (not only with the jsb_var op_ov)
    structOut%carbon%ptr      => structOut%carbon%ptr2d(:,:)
    structOut%nitrogen%ptr    => structOut%nitrogen%ptr2d(:,:)
    structOut%phosphorus%ptr  => structOut%phosphorus%ptr2d(:,:)
    structOut%carbon13%ptr    => structOut%carbon13%ptr2d(:,:)
    structOut%carbon14%ptr    => structOut%carbon14%ptr2d(:,:)
    structOut%nitrogen15%ptr  => structOut%nitrogen15%ptr2d(:,:)
    !> 2.2 element pools
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
    TYPE(t_veg_tissues_pool),       TARGET      :: structOut
  
    CHARACTER(len=*), PARAMETER :: routine = modname//':Add_veg_tiss_p_and_veg_tiss_p'
  
    !> 1.0 init temporary pool structOut
    !! init tissue pool with element variables
    CALL structOut%Init_tiss_pool_for_op_ov(name='op_ov_help', shortname='op_ov_help', &
      & unit=TRIM(struct1%element_unit), pool_suffix='veg_p') ! l_elements=.TRUE. by default
    structOut%carbon%owner_model_id      = struct1%carbon%owner_model_id
    structOut%carbon%owner_proc_id       = struct1%carbon%owner_proc_id
    structOut%carbon%owner_tile_path     = struct1%carbon%owner_tile_path
    ALLOCATE(structOut%carbon%ptr2d,     source=struct1%carbon%ptr2d)
    structOut%nitrogen%owner_model_id    = struct1%nitrogen%owner_model_id
    structOut%nitrogen%owner_proc_id     = struct1%nitrogen%owner_proc_id
    structOut%nitrogen%owner_tile_path   = struct1%nitrogen%owner_tile_path
    ALLOCATE(structOut%nitrogen%ptr2d,   source=struct1%nitrogen%ptr2d)
    structOut%phosphorus%owner_model_id  = struct1%phosphorus%owner_model_id
    structOut%phosphorus%owner_proc_id   = struct1%phosphorus%owner_proc_id
    structOut%phosphorus%owner_tile_path = struct1%phosphorus%owner_tile_path
    ALLOCATE(structOut%phosphorus%ptr2d, source=struct1%phosphorus%ptr2d)
    structOut%carbon13%owner_model_id    = struct1%carbon13%owner_model_id
    structOut%carbon13%owner_proc_id     = struct1%carbon13%owner_proc_id
    structOut%carbon13%owner_tile_path   = struct1%carbon13%owner_tile_path
    ALLOCATE(structOut%carbon13%ptr2d,   source=struct1%carbon13%ptr2d)
    structOut%carbon14%owner_model_id    = struct1%carbon14%owner_model_id
    structOut%carbon14%owner_proc_id     = struct1%carbon14%owner_proc_id
    structOut%carbon14%owner_tile_path   = struct1%carbon14%owner_tile_path
    ALLOCATE(structOut%carbon14%ptr2d,   source=struct1%carbon14%ptr2d)
    structOut%nitrogen15%owner_model_id  = struct1%nitrogen15%owner_model_id
    structOut%nitrogen15%owner_proc_id   = struct1%nitrogen15%owner_proc_id
    structOut%nitrogen15%owner_tile_path = struct1%nitrogen15%owner_tile_path
    ALLOCATE(structOut%nitrogen15%ptr2d, source=struct1%nitrogen15%ptr2d)
    !! init element pools
    CALL structOut%leaf%Init_elem_var_elem_p(struct1%leaf)
    CALL structOut%fine_root%Init_elem_var_elem_p(struct1%fine_root)
    CALL structOut%coarse_root%Init_elem_var_elem_p(struct1%coarse_root)
    CALL structOut%sap_wood%Init_elem_var_elem_p(struct1%sap_wood)
    CALL structOut%heart_wood%Init_elem_var_elem_p(struct1%heart_wood)
    CALL structOut%labile%Init_elem_var_elem_p(struct1%labile)
    CALL structOut%reserve%Init_elem_var_elem_p(struct1%reserve)
    CALL structOut%fruit%Init_elem_var_elem_p(struct1%fruit)
    CALL structOut%seed_bed%Init_elem_var_elem_p(struct1%seed_bed)
  
    !> 2.0 run calculation
    !! the differentiation between ON_CHUNK/ON_DOMAIN is done at the level of variables (i.e., the op_ov routines of the jsb_var)
    !> 2.1 element variables of the tissue pool itself
    structOut%carbon          = struct1%carbon      + struct2%carbon    
    structOut%nitrogen        = struct1%nitrogen    + struct2%nitrogen  
    structOut%phosphorus      = struct1%phosphorus  + struct2%phosphorus
    structOut%carbon13        = struct1%carbon13    + struct2%carbon13  
    structOut%carbon14        = struct1%carbon14    + struct2%carbon14  
    structOut%nitrogen15      = struct1%nitrogen15  + struct2%nitrogen15
    ! adding the pointer is this needed here too (not only with the jsb_var op_ov)
    structOut%carbon%ptr      => structOut%carbon%ptr2d(:,:)
    structOut%nitrogen%ptr    => structOut%nitrogen%ptr2d(:,:)
    structOut%phosphorus%ptr  => structOut%phosphorus%ptr2d(:,:)
    structOut%carbon13%ptr    => structOut%carbon13%ptr2d(:,:)
    structOut%carbon14%ptr    => structOut%carbon14%ptr2d(:,:)
    structOut%nitrogen15%ptr  => structOut%nitrogen15%ptr2d(:,:)
    !> 2.2 element pools
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
    TYPE(t_veg_tissues_pool),       TARGET      :: structOut
  
    CHARACTER(len=*), PARAMETER :: routine = modname//':Subtract_veg_tiss_p_minus_scalar'
  
    !> 1.0 init temporary pool structOut
    !! init tissue pool with element variables
    CALL structOut%Init_tiss_pool_for_op_ov(name='op_ov_help', shortname='op_ov_help', &
      & unit=TRIM(struct1%element_unit), pool_suffix='veg_p') ! l_elements=.TRUE. by default
    structOut%carbon%owner_model_id      = struct1%carbon%owner_model_id
    structOut%carbon%owner_proc_id       = struct1%carbon%owner_proc_id
    structOut%carbon%owner_tile_path     = struct1%carbon%owner_tile_path
    ALLOCATE(structOut%carbon%ptr2d,     source=struct1%carbon%ptr2d)
    structOut%nitrogen%owner_model_id    = struct1%nitrogen%owner_model_id
    structOut%nitrogen%owner_proc_id     = struct1%nitrogen%owner_proc_id
    structOut%nitrogen%owner_tile_path   = struct1%nitrogen%owner_tile_path
    ALLOCATE(structOut%nitrogen%ptr2d,   source=struct1%nitrogen%ptr2d)
    structOut%phosphorus%owner_model_id  = struct1%phosphorus%owner_model_id
    structOut%phosphorus%owner_proc_id   = struct1%phosphorus%owner_proc_id
    structOut%phosphorus%owner_tile_path = struct1%phosphorus%owner_tile_path
    ALLOCATE(structOut%phosphorus%ptr2d, source=struct1%phosphorus%ptr2d)
    structOut%carbon13%owner_model_id    = struct1%carbon13%owner_model_id
    structOut%carbon13%owner_proc_id     = struct1%carbon13%owner_proc_id
    structOut%carbon13%owner_tile_path   = struct1%carbon13%owner_tile_path
    ALLOCATE(structOut%carbon13%ptr2d,   source=struct1%carbon13%ptr2d)
    structOut%carbon14%owner_model_id    = struct1%carbon14%owner_model_id
    structOut%carbon14%owner_proc_id     = struct1%carbon14%owner_proc_id
    structOut%carbon14%owner_tile_path   = struct1%carbon14%owner_tile_path
    ALLOCATE(structOut%carbon14%ptr2d,   source=struct1%carbon14%ptr2d)
    structOut%nitrogen15%owner_model_id  = struct1%nitrogen15%owner_model_id
    structOut%nitrogen15%owner_proc_id   = struct1%nitrogen15%owner_proc_id
    structOut%nitrogen15%owner_tile_path = struct1%nitrogen15%owner_tile_path
    ALLOCATE(structOut%nitrogen15%ptr2d, source=struct1%nitrogen15%ptr2d)
    !! init element pools
    CALL structOut%leaf%Init_elem_var_elem_p(struct1%leaf)
    CALL structOut%fine_root%Init_elem_var_elem_p(struct1%fine_root)
    CALL structOut%coarse_root%Init_elem_var_elem_p(struct1%coarse_root)
    CALL structOut%sap_wood%Init_elem_var_elem_p(struct1%sap_wood)
    CALL structOut%heart_wood%Init_elem_var_elem_p(struct1%heart_wood)
    CALL structOut%labile%Init_elem_var_elem_p(struct1%labile)
    CALL structOut%reserve%Init_elem_var_elem_p(struct1%reserve)
    CALL structOut%fruit%Init_elem_var_elem_p(struct1%fruit)
    CALL structOut%seed_bed%Init_elem_var_elem_p(struct1%seed_bed)
  
    !> 2.0 run calculation
    !! the differentiation between ON_CHUNK/ON_DOMAIN is done at the level of variables (i.e., the op_ov routines of the jsb_var)
    !> 2.1 element variables of the tissue pool itself
    structOut%carbon          = struct1%carbon      - scalar
    structOut%nitrogen        = struct1%nitrogen    - scalar
    structOut%phosphorus      = struct1%phosphorus  - scalar
    structOut%carbon13        = struct1%carbon13    - scalar
    structOut%carbon14        = struct1%carbon14    - scalar
    structOut%nitrogen15      = struct1%nitrogen15  - scalar
    ! adding the pointer is this needed here too (not only with the jsb_var op_ov)
    structOut%carbon%ptr      => structOut%carbon%ptr2d(:,:)
    structOut%nitrogen%ptr    => structOut%nitrogen%ptr2d(:,:)
    structOut%phosphorus%ptr  => structOut%phosphorus%ptr2d(:,:)
    structOut%carbon13%ptr    => structOut%carbon13%ptr2d(:,:)
    structOut%carbon14%ptr    => structOut%carbon14%ptr2d(:,:)
    structOut%nitrogen15%ptr  => structOut%nitrogen15%ptr2d(:,:)
    !> 2.2 element pools
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
    TYPE(t_veg_tissues_pool),       TARGET      :: structOut
  
    CHARACTER(len=*), PARAMETER :: routine = modname//':Subtract_veg_tiss_p_minus_field1d'
  
    !> 1.0 init temporary pool structOut
    !! init tissue pool with element variables
    CALL structOut%Init_tiss_pool_for_op_ov(name='op_ov_help', shortname='op_ov_help', &
      & unit=TRIM(struct1%element_unit), pool_suffix='veg_p') ! l_elements=.TRUE. by default
    structOut%carbon%owner_model_id      = struct1%carbon%owner_model_id
    structOut%carbon%owner_proc_id       = struct1%carbon%owner_proc_id
    structOut%carbon%owner_tile_path     = struct1%carbon%owner_tile_path
    ALLOCATE(structOut%carbon%ptr2d,     source=struct1%carbon%ptr2d)
    structOut%nitrogen%owner_model_id    = struct1%nitrogen%owner_model_id
    structOut%nitrogen%owner_proc_id     = struct1%nitrogen%owner_proc_id
    structOut%nitrogen%owner_tile_path   = struct1%nitrogen%owner_tile_path
    ALLOCATE(structOut%nitrogen%ptr2d,   source=struct1%nitrogen%ptr2d)
    structOut%phosphorus%owner_model_id  = struct1%phosphorus%owner_model_id
    structOut%phosphorus%owner_proc_id   = struct1%phosphorus%owner_proc_id
    structOut%phosphorus%owner_tile_path = struct1%phosphorus%owner_tile_path
    ALLOCATE(structOut%phosphorus%ptr2d, source=struct1%phosphorus%ptr2d)
    structOut%carbon13%owner_model_id    = struct1%carbon13%owner_model_id
    structOut%carbon13%owner_proc_id     = struct1%carbon13%owner_proc_id
    structOut%carbon13%owner_tile_path   = struct1%carbon13%owner_tile_path
    ALLOCATE(structOut%carbon13%ptr2d,   source=struct1%carbon13%ptr2d)
    structOut%carbon14%owner_model_id    = struct1%carbon14%owner_model_id
    structOut%carbon14%owner_proc_id     = struct1%carbon14%owner_proc_id
    structOut%carbon14%owner_tile_path   = struct1%carbon14%owner_tile_path
    ALLOCATE(structOut%carbon14%ptr2d,   source=struct1%carbon14%ptr2d)
    structOut%nitrogen15%owner_model_id  = struct1%nitrogen15%owner_model_id
    structOut%nitrogen15%owner_proc_id   = struct1%nitrogen15%owner_proc_id
    structOut%nitrogen15%owner_tile_path = struct1%nitrogen15%owner_tile_path
    ALLOCATE(structOut%nitrogen15%ptr2d, source=struct1%nitrogen15%ptr2d)
    !! init element pools
    CALL structOut%leaf%Init_elem_var_elem_p(struct1%leaf)
    CALL structOut%fine_root%Init_elem_var_elem_p(struct1%fine_root)
    CALL structOut%coarse_root%Init_elem_var_elem_p(struct1%coarse_root)
    CALL structOut%sap_wood%Init_elem_var_elem_p(struct1%sap_wood)
    CALL structOut%heart_wood%Init_elem_var_elem_p(struct1%heart_wood)
    CALL structOut%labile%Init_elem_var_elem_p(struct1%labile)
    CALL structOut%reserve%Init_elem_var_elem_p(struct1%reserve)
    CALL structOut%fruit%Init_elem_var_elem_p(struct1%fruit)
    CALL structOut%seed_bed%Init_elem_var_elem_p(struct1%seed_bed)
  
    !> 2.0 run calculation
    !! the differentiation between ON_CHUNK/ON_DOMAIN is done at the level of variables (i.e., the op_ov routines of the jsb_var)
    !> 2.1 element variables of the tissue pool itself
    structOut%carbon          = struct1%carbon      - field1d(:)
    structOut%nitrogen        = struct1%nitrogen    - field1d(:)
    structOut%phosphorus      = struct1%phosphorus  - field1d(:)
    structOut%carbon13        = struct1%carbon13    - field1d(:)
    structOut%carbon14        = struct1%carbon14    - field1d(:)
    structOut%nitrogen15      = struct1%nitrogen15  - field1d(:)
    ! adding the pointer is this needed here too (not only with the jsb_var op_ov)
    structOut%carbon%ptr      => structOut%carbon%ptr2d(:,:)
    structOut%nitrogen%ptr    => structOut%nitrogen%ptr2d(:,:)
    structOut%phosphorus%ptr  => structOut%phosphorus%ptr2d(:,:)
    structOut%carbon13%ptr    => structOut%carbon13%ptr2d(:,:)
    structOut%carbon14%ptr    => structOut%carbon14%ptr2d(:,:)
    structOut%nitrogen15%ptr  => structOut%nitrogen15%ptr2d(:,:)
    !> 2.2 element pools
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
    TYPE(t_veg_tissues_pool),       TARGET      :: structOut
  
    CHARACTER(len=*), PARAMETER :: routine = modname//':Subtract_veg_tiss_p_minus_field2d'
  
    !> 1.0 init temporary pool structOut
    !! init tissue pool with element variables
    CALL structOut%Init_tiss_pool_for_op_ov(name='op_ov_help', shortname='op_ov_help', &
      & unit=TRIM(struct1%element_unit), pool_suffix='veg_p') ! l_elements=.TRUE. by default
    structOut%carbon%owner_model_id      = struct1%carbon%owner_model_id
    structOut%carbon%owner_proc_id       = struct1%carbon%owner_proc_id
    structOut%carbon%owner_tile_path     = struct1%carbon%owner_tile_path
    ALLOCATE(structOut%carbon%ptr2d,     source=struct1%carbon%ptr2d)
    structOut%nitrogen%owner_model_id    = struct1%nitrogen%owner_model_id
    structOut%nitrogen%owner_proc_id     = struct1%nitrogen%owner_proc_id
    structOut%nitrogen%owner_tile_path   = struct1%nitrogen%owner_tile_path
    ALLOCATE(structOut%nitrogen%ptr2d,   source=struct1%nitrogen%ptr2d)
    structOut%phosphorus%owner_model_id  = struct1%phosphorus%owner_model_id
    structOut%phosphorus%owner_proc_id   = struct1%phosphorus%owner_proc_id
    structOut%phosphorus%owner_tile_path = struct1%phosphorus%owner_tile_path
    ALLOCATE(structOut%phosphorus%ptr2d, source=struct1%phosphorus%ptr2d)
    structOut%carbon13%owner_model_id    = struct1%carbon13%owner_model_id
    structOut%carbon13%owner_proc_id     = struct1%carbon13%owner_proc_id
    structOut%carbon13%owner_tile_path   = struct1%carbon13%owner_tile_path
    ALLOCATE(structOut%carbon13%ptr2d,   source=struct1%carbon13%ptr2d)
    structOut%carbon14%owner_model_id    = struct1%carbon14%owner_model_id
    structOut%carbon14%owner_proc_id     = struct1%carbon14%owner_proc_id
    structOut%carbon14%owner_tile_path   = struct1%carbon14%owner_tile_path
    ALLOCATE(structOut%carbon14%ptr2d,   source=struct1%carbon14%ptr2d)
    structOut%nitrogen15%owner_model_id  = struct1%nitrogen15%owner_model_id
    structOut%nitrogen15%owner_proc_id   = struct1%nitrogen15%owner_proc_id
    structOut%nitrogen15%owner_tile_path = struct1%nitrogen15%owner_tile_path
    ALLOCATE(structOut%nitrogen15%ptr2d, source=struct1%nitrogen15%ptr2d)
    !! init element pools
    CALL structOut%leaf%Init_elem_var_elem_p(struct1%leaf)
    CALL structOut%fine_root%Init_elem_var_elem_p(struct1%fine_root)
    CALL structOut%coarse_root%Init_elem_var_elem_p(struct1%coarse_root)
    CALL structOut%sap_wood%Init_elem_var_elem_p(struct1%sap_wood)
    CALL structOut%heart_wood%Init_elem_var_elem_p(struct1%heart_wood)
    CALL structOut%labile%Init_elem_var_elem_p(struct1%labile)
    CALL structOut%reserve%Init_elem_var_elem_p(struct1%reserve)
    CALL structOut%fruit%Init_elem_var_elem_p(struct1%fruit)
    CALL structOut%seed_bed%Init_elem_var_elem_p(struct1%seed_bed)
  
    !> 2.0 run calculation
    !! the differentiation between ON_CHUNK/ON_DOMAIN is done at the level of variables (i.e., the op_ov routines of the jsb_var)
    !> 2.1 element variables of the tissue pool itself
    structOut%carbon          = struct1%carbon      - field2d(:,:)
    structOut%nitrogen        = struct1%nitrogen    - field2d(:,:)
    structOut%phosphorus      = struct1%phosphorus  - field2d(:,:)
    structOut%carbon13        = struct1%carbon13    - field2d(:,:)
    structOut%carbon14        = struct1%carbon14    - field2d(:,:)
    structOut%nitrogen15      = struct1%nitrogen15  - field2d(:,:)
    ! adding the pointer is this needed here too (not only with the jsb_var op_ov)
    structOut%carbon%ptr      => structOut%carbon%ptr2d(:,:)
    structOut%nitrogen%ptr    => structOut%nitrogen%ptr2d(:,:)
    structOut%phosphorus%ptr  => structOut%phosphorus%ptr2d(:,:)
    structOut%carbon13%ptr    => structOut%carbon13%ptr2d(:,:)
    structOut%carbon14%ptr    => structOut%carbon14%ptr2d(:,:)
    structOut%nitrogen15%ptr  => structOut%nitrogen15%ptr2d(:,:)
    !> 2.2 element pools
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
    TYPE(t_veg_tissues_pool),       TARGET      :: structOut
  
    CHARACTER(len=*), PARAMETER :: routine = modname//':Subtract_veg_tiss_p_minus_veg_tiss_p'
  
    !> 1.0 init temporary pool structOut
    !! init tissue pool with element variables
    CALL structOut%Init_tiss_pool_for_op_ov(name='op_ov_help', shortname='op_ov_help', &
      & unit=TRIM(struct1%element_unit), pool_suffix='veg_p') ! l_elements=.TRUE. by default
    structOut%carbon%owner_model_id      = struct1%carbon%owner_model_id
    structOut%carbon%owner_proc_id       = struct1%carbon%owner_proc_id
    structOut%carbon%owner_tile_path     = struct1%carbon%owner_tile_path
    ALLOCATE(structOut%carbon%ptr2d,     source=struct1%carbon%ptr2d)
    structOut%nitrogen%owner_model_id    = struct1%nitrogen%owner_model_id
    structOut%nitrogen%owner_proc_id     = struct1%nitrogen%owner_proc_id
    structOut%nitrogen%owner_tile_path   = struct1%nitrogen%owner_tile_path
    ALLOCATE(structOut%nitrogen%ptr2d,   source=struct1%nitrogen%ptr2d)
    structOut%phosphorus%owner_model_id  = struct1%phosphorus%owner_model_id
    structOut%phosphorus%owner_proc_id   = struct1%phosphorus%owner_proc_id
    structOut%phosphorus%owner_tile_path = struct1%phosphorus%owner_tile_path
    ALLOCATE(structOut%phosphorus%ptr2d, source=struct1%phosphorus%ptr2d)
    structOut%carbon13%owner_model_id    = struct1%carbon13%owner_model_id
    structOut%carbon13%owner_proc_id     = struct1%carbon13%owner_proc_id
    structOut%carbon13%owner_tile_path   = struct1%carbon13%owner_tile_path
    ALLOCATE(structOut%carbon13%ptr2d,   source=struct1%carbon13%ptr2d)
    structOut%carbon14%owner_model_id    = struct1%carbon14%owner_model_id
    structOut%carbon14%owner_proc_id     = struct1%carbon14%owner_proc_id
    structOut%carbon14%owner_tile_path   = struct1%carbon14%owner_tile_path
    ALLOCATE(structOut%carbon14%ptr2d,   source=struct1%carbon14%ptr2d)
    structOut%nitrogen15%owner_model_id  = struct1%nitrogen15%owner_model_id
    structOut%nitrogen15%owner_proc_id   = struct1%nitrogen15%owner_proc_id
    structOut%nitrogen15%owner_tile_path = struct1%nitrogen15%owner_tile_path
    ALLOCATE(structOut%nitrogen15%ptr2d, source=struct1%nitrogen15%ptr2d)
    !! init element pools
    CALL structOut%leaf%Init_elem_var_elem_p(struct1%leaf)
    CALL structOut%fine_root%Init_elem_var_elem_p(struct1%fine_root)
    CALL structOut%coarse_root%Init_elem_var_elem_p(struct1%coarse_root)
    CALL structOut%sap_wood%Init_elem_var_elem_p(struct1%sap_wood)
    CALL structOut%heart_wood%Init_elem_var_elem_p(struct1%heart_wood)
    CALL structOut%labile%Init_elem_var_elem_p(struct1%labile)
    CALL structOut%reserve%Init_elem_var_elem_p(struct1%reserve)
    CALL structOut%fruit%Init_elem_var_elem_p(struct1%fruit)
    CALL structOut%seed_bed%Init_elem_var_elem_p(struct1%seed_bed)
  
    !> 2.0 run calculation
    !! the differentiation between ON_CHUNK/ON_DOMAIN is done at the level of variables (i.e., the op_ov routines of the jsb_var)
    !> 2.1 element variables of the tissue pool itself
    structOut%carbon          = struct1%carbon      - struct2%carbon    
    structOut%nitrogen        = struct1%nitrogen    - struct2%nitrogen  
    structOut%phosphorus      = struct1%phosphorus  - struct2%phosphorus
    structOut%carbon13        = struct1%carbon13    - struct2%carbon13  
    structOut%carbon14        = struct1%carbon14    - struct2%carbon14  
    structOut%nitrogen15      = struct1%nitrogen15  - struct2%nitrogen15
    ! adding the pointer is this needed here too (not only with the jsb_var op_ov)
    structOut%carbon%ptr      => structOut%carbon%ptr2d(:,:)
    structOut%nitrogen%ptr    => structOut%nitrogen%ptr2d(:,:)
    structOut%phosphorus%ptr  => structOut%phosphorus%ptr2d(:,:)
    structOut%carbon13%ptr    => structOut%carbon13%ptr2d(:,:)
    structOut%carbon14%ptr    => structOut%carbon14%ptr2d(:,:)
    structOut%nitrogen15%ptr  => structOut%nitrogen15%ptr2d(:,:)
    !> 2.2 element pools
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
    TYPE(t_veg_tissues_pool),       TARGET      :: structOut
  
    CHARACTER(len=*), PARAMETER :: routine = modname//':Multiply_veg_tiss_p_with_scalar'
  
    !> 1.0 init temporary pool structOut
    !! init tissue pool with element variables
    CALL structOut%Init_tiss_pool_for_op_ov(name='op_ov_help', shortname='op_ov_help', &
      & unit=TRIM(struct1%element_unit), pool_suffix='veg_p') ! l_elements=.TRUE. by default
    structOut%carbon%owner_model_id      = struct1%carbon%owner_model_id
    structOut%carbon%owner_proc_id       = struct1%carbon%owner_proc_id
    structOut%carbon%owner_tile_path     = struct1%carbon%owner_tile_path
    ALLOCATE(structOut%carbon%ptr2d,     source=struct1%carbon%ptr2d)
    structOut%nitrogen%owner_model_id    = struct1%nitrogen%owner_model_id
    structOut%nitrogen%owner_proc_id     = struct1%nitrogen%owner_proc_id
    structOut%nitrogen%owner_tile_path   = struct1%nitrogen%owner_tile_path
    ALLOCATE(structOut%nitrogen%ptr2d,   source=struct1%nitrogen%ptr2d)
    structOut%phosphorus%owner_model_id  = struct1%phosphorus%owner_model_id
    structOut%phosphorus%owner_proc_id   = struct1%phosphorus%owner_proc_id
    structOut%phosphorus%owner_tile_path = struct1%phosphorus%owner_tile_path
    ALLOCATE(structOut%phosphorus%ptr2d, source=struct1%phosphorus%ptr2d)
    structOut%carbon13%owner_model_id    = struct1%carbon13%owner_model_id
    structOut%carbon13%owner_proc_id     = struct1%carbon13%owner_proc_id
    structOut%carbon13%owner_tile_path   = struct1%carbon13%owner_tile_path
    ALLOCATE(structOut%carbon13%ptr2d,   source=struct1%carbon13%ptr2d)
    structOut%carbon14%owner_model_id    = struct1%carbon14%owner_model_id
    structOut%carbon14%owner_proc_id     = struct1%carbon14%owner_proc_id
    structOut%carbon14%owner_tile_path   = struct1%carbon14%owner_tile_path
    ALLOCATE(structOut%carbon14%ptr2d,   source=struct1%carbon14%ptr2d)
    structOut%nitrogen15%owner_model_id  = struct1%nitrogen15%owner_model_id
    structOut%nitrogen15%owner_proc_id   = struct1%nitrogen15%owner_proc_id
    structOut%nitrogen15%owner_tile_path = struct1%nitrogen15%owner_tile_path
    ALLOCATE(structOut%nitrogen15%ptr2d, source=struct1%nitrogen15%ptr2d)
    !! init element pools
    CALL structOut%leaf%Init_elem_var_elem_p(struct1%leaf)
    CALL structOut%fine_root%Init_elem_var_elem_p(struct1%fine_root)
    CALL structOut%coarse_root%Init_elem_var_elem_p(struct1%coarse_root)
    CALL structOut%sap_wood%Init_elem_var_elem_p(struct1%sap_wood)
    CALL structOut%heart_wood%Init_elem_var_elem_p(struct1%heart_wood)
    CALL structOut%labile%Init_elem_var_elem_p(struct1%labile)
    CALL structOut%reserve%Init_elem_var_elem_p(struct1%reserve)
    CALL structOut%fruit%Init_elem_var_elem_p(struct1%fruit)
    CALL structOut%seed_bed%Init_elem_var_elem_p(struct1%seed_bed)
  
    !> 2.0 run calculation
    !! the differentiation between ON_CHUNK/ON_DOMAIN is done at the level of variables (i.e., the op_ov routines of the jsb_var)
    !> 2.1 element variables of the tissue pool itself
    structOut%carbon          = struct1%carbon      * scalar
    structOut%nitrogen        = struct1%nitrogen    * scalar
    structOut%phosphorus      = struct1%phosphorus  * scalar
    structOut%carbon13        = struct1%carbon13    * scalar
    structOut%carbon14        = struct1%carbon14    * scalar
    structOut%nitrogen15      = struct1%nitrogen15  * scalar
    ! adding the pointer is this needed here too (not only with the jsb_var op_ov)
    structOut%carbon%ptr      => structOut%carbon%ptr2d(:,:)
    structOut%nitrogen%ptr    => structOut%nitrogen%ptr2d(:,:)
    structOut%phosphorus%ptr  => structOut%phosphorus%ptr2d(:,:)
    structOut%carbon13%ptr    => structOut%carbon13%ptr2d(:,:)
    structOut%carbon14%ptr    => structOut%carbon14%ptr2d(:,:)
    structOut%nitrogen15%ptr  => structOut%nitrogen15%ptr2d(:,:)
    !> 2.2 element pools
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
    TYPE(t_veg_tissues_pool),       TARGET      :: structOut
  
    CHARACTER(len=*), PARAMETER :: routine = modname//':Multiply_scalar_with_veg_tiss_p'
  
    !> 1.0 init temporary pool structOut
    !! init tissue pool with element variables
    CALL structOut%Init_tiss_pool_for_op_ov(name='op_ov_help', shortname='op_ov_help', &
      & unit=TRIM(struct1%element_unit), pool_suffix='veg_p') ! l_elements=.TRUE. by default
    structOut%carbon%owner_model_id      = struct1%carbon%owner_model_id
    structOut%carbon%owner_proc_id       = struct1%carbon%owner_proc_id
    structOut%carbon%owner_tile_path     = struct1%carbon%owner_tile_path
    ALLOCATE(structOut%carbon%ptr2d,     source=struct1%carbon%ptr2d)
    structOut%nitrogen%owner_model_id    = struct1%nitrogen%owner_model_id
    structOut%nitrogen%owner_proc_id     = struct1%nitrogen%owner_proc_id
    structOut%nitrogen%owner_tile_path   = struct1%nitrogen%owner_tile_path
    ALLOCATE(structOut%nitrogen%ptr2d,   source=struct1%nitrogen%ptr2d)
    structOut%phosphorus%owner_model_id  = struct1%phosphorus%owner_model_id
    structOut%phosphorus%owner_proc_id   = struct1%phosphorus%owner_proc_id
    structOut%phosphorus%owner_tile_path = struct1%phosphorus%owner_tile_path
    ALLOCATE(structOut%phosphorus%ptr2d, source=struct1%phosphorus%ptr2d)
    structOut%carbon13%owner_model_id    = struct1%carbon13%owner_model_id
    structOut%carbon13%owner_proc_id     = struct1%carbon13%owner_proc_id
    structOut%carbon13%owner_tile_path   = struct1%carbon13%owner_tile_path
    ALLOCATE(structOut%carbon13%ptr2d,   source=struct1%carbon13%ptr2d)
    structOut%carbon14%owner_model_id    = struct1%carbon14%owner_model_id
    structOut%carbon14%owner_proc_id     = struct1%carbon14%owner_proc_id
    structOut%carbon14%owner_tile_path   = struct1%carbon14%owner_tile_path
    ALLOCATE(structOut%carbon14%ptr2d,   source=struct1%carbon14%ptr2d)
    structOut%nitrogen15%owner_model_id  = struct1%nitrogen15%owner_model_id
    structOut%nitrogen15%owner_proc_id   = struct1%nitrogen15%owner_proc_id
    structOut%nitrogen15%owner_tile_path = struct1%nitrogen15%owner_tile_path
    ALLOCATE(structOut%nitrogen15%ptr2d, source=struct1%nitrogen15%ptr2d)
    !! init element pools
    CALL structOut%leaf%Init_elem_var_elem_p(struct1%leaf)
    CALL structOut%fine_root%Init_elem_var_elem_p(struct1%fine_root)
    CALL structOut%coarse_root%Init_elem_var_elem_p(struct1%coarse_root)
    CALL structOut%sap_wood%Init_elem_var_elem_p(struct1%sap_wood)
    CALL structOut%heart_wood%Init_elem_var_elem_p(struct1%heart_wood)
    CALL structOut%labile%Init_elem_var_elem_p(struct1%labile)
    CALL structOut%reserve%Init_elem_var_elem_p(struct1%reserve)
    CALL structOut%fruit%Init_elem_var_elem_p(struct1%fruit)
    CALL structOut%seed_bed%Init_elem_var_elem_p(struct1%seed_bed)
  
    !> 2.0 run calculation
    !! the differentiation between ON_CHUNK/ON_DOMAIN is done at the level of variables (i.e., the op_ov routines of the jsb_var)
    !> 2.1 element variables of the tissue pool itself
    structOut%carbon          = scalar * struct1%carbon    
    structOut%nitrogen        = scalar * struct1%nitrogen  
    structOut%phosphorus      = scalar * struct1%phosphorus
    structOut%carbon13        = scalar * struct1%carbon13  
    structOut%carbon14        = scalar * struct1%carbon14  
    structOut%nitrogen15      = scalar * struct1%nitrogen15
    ! adding the pointer is this needed here too (not only with the jsb_var op_ov)
    structOut%carbon%ptr      => structOut%carbon%ptr2d(:,:)
    structOut%nitrogen%ptr    => structOut%nitrogen%ptr2d(:,:)
    structOut%phosphorus%ptr  => structOut%phosphorus%ptr2d(:,:)
    structOut%carbon13%ptr    => structOut%carbon13%ptr2d(:,:)
    structOut%carbon14%ptr    => structOut%carbon14%ptr2d(:,:)
    structOut%nitrogen15%ptr  => structOut%nitrogen15%ptr2d(:,:)
    !> 2.2 element pools
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
    TYPE(t_veg_tissues_pool),       TARGET      :: structOut
  
    CHARACTER(len=*), PARAMETER :: routine = modname//':Multiply_veg_tiss_p_with_field1d'
  
    !> 1.0 init temporary pool structOut
    !! init tissue pool with element variables
    CALL structOut%Init_tiss_pool_for_op_ov(name='op_ov_help', shortname='op_ov_help', &
      & unit=TRIM(struct1%element_unit), pool_suffix='veg_p') ! l_elements=.TRUE. by default
    structOut%carbon%owner_model_id      = struct1%carbon%owner_model_id
    structOut%carbon%owner_proc_id       = struct1%carbon%owner_proc_id
    structOut%carbon%owner_tile_path     = struct1%carbon%owner_tile_path
    ALLOCATE(structOut%carbon%ptr2d,     source=struct1%carbon%ptr2d)
    structOut%nitrogen%owner_model_id    = struct1%nitrogen%owner_model_id
    structOut%nitrogen%owner_proc_id     = struct1%nitrogen%owner_proc_id
    structOut%nitrogen%owner_tile_path   = struct1%nitrogen%owner_tile_path
    ALLOCATE(structOut%nitrogen%ptr2d,   source=struct1%nitrogen%ptr2d)
    structOut%phosphorus%owner_model_id  = struct1%phosphorus%owner_model_id
    structOut%phosphorus%owner_proc_id   = struct1%phosphorus%owner_proc_id
    structOut%phosphorus%owner_tile_path = struct1%phosphorus%owner_tile_path
    ALLOCATE(structOut%phosphorus%ptr2d, source=struct1%phosphorus%ptr2d)
    structOut%carbon13%owner_model_id    = struct1%carbon13%owner_model_id
    structOut%carbon13%owner_proc_id     = struct1%carbon13%owner_proc_id
    structOut%carbon13%owner_tile_path   = struct1%carbon13%owner_tile_path
    ALLOCATE(structOut%carbon13%ptr2d,   source=struct1%carbon13%ptr2d)
    structOut%carbon14%owner_model_id    = struct1%carbon14%owner_model_id
    structOut%carbon14%owner_proc_id     = struct1%carbon14%owner_proc_id
    structOut%carbon14%owner_tile_path   = struct1%carbon14%owner_tile_path
    ALLOCATE(structOut%carbon14%ptr2d,   source=struct1%carbon14%ptr2d)
    structOut%nitrogen15%owner_model_id  = struct1%nitrogen15%owner_model_id
    structOut%nitrogen15%owner_proc_id   = struct1%nitrogen15%owner_proc_id
    structOut%nitrogen15%owner_tile_path = struct1%nitrogen15%owner_tile_path
    ALLOCATE(structOut%nitrogen15%ptr2d, source=struct1%nitrogen15%ptr2d)
    !! init element pools
    CALL structOut%leaf%Init_elem_var_elem_p(struct1%leaf)
    CALL structOut%fine_root%Init_elem_var_elem_p(struct1%fine_root)
    CALL structOut%coarse_root%Init_elem_var_elem_p(struct1%coarse_root)
    CALL structOut%sap_wood%Init_elem_var_elem_p(struct1%sap_wood)
    CALL structOut%heart_wood%Init_elem_var_elem_p(struct1%heart_wood)
    CALL structOut%labile%Init_elem_var_elem_p(struct1%labile)
    CALL structOut%reserve%Init_elem_var_elem_p(struct1%reserve)
    CALL structOut%fruit%Init_elem_var_elem_p(struct1%fruit)
    CALL structOut%seed_bed%Init_elem_var_elem_p(struct1%seed_bed)
  
    !> 2.0 run calculation
    !! the differentiation between ON_CHUNK/ON_DOMAIN is done at the level of variables (i.e., the op_ov routines of the jsb_var)
    !> 2.1 element variables of the tissue pool itself
    structOut%carbon          = struct1%carbon      * field1d(:)
    structOut%nitrogen        = struct1%nitrogen    * field1d(:)
    structOut%phosphorus      = struct1%phosphorus  * field1d(:)
    structOut%carbon13        = struct1%carbon13    * field1d(:)
    structOut%carbon14        = struct1%carbon14    * field1d(:)
    structOut%nitrogen15      = struct1%nitrogen15  * field1d(:)
    ! adding the pointer is this needed here too (not only with the jsb_var op_ov)
    structOut%carbon%ptr      => structOut%carbon%ptr2d(:,:)
    structOut%nitrogen%ptr    => structOut%nitrogen%ptr2d(:,:)
    structOut%phosphorus%ptr  => structOut%phosphorus%ptr2d(:,:)
    structOut%carbon13%ptr    => structOut%carbon13%ptr2d(:,:)
    structOut%carbon14%ptr    => structOut%carbon14%ptr2d(:,:)
    structOut%nitrogen15%ptr  => structOut%nitrogen15%ptr2d(:,:)
    !> 2.2 element pools
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
    TYPE(t_veg_tissues_pool),       TARGET      :: structOut
  
    CHARACTER(len=*), PARAMETER :: routine = modname//':Multiply_field1d_with_veg_tiss_p'
  
    !> 1.0 init temporary pool structOut
    !! init tissue pool with element variables
    CALL structOut%Init_tiss_pool_for_op_ov(name='op_ov_help', shortname='op_ov_help', &
      & unit=TRIM(struct1%element_unit), pool_suffix='veg_p') ! l_elements=.TRUE. by default
    structOut%carbon%owner_model_id      = struct1%carbon%owner_model_id
    structOut%carbon%owner_proc_id       = struct1%carbon%owner_proc_id
    structOut%carbon%owner_tile_path     = struct1%carbon%owner_tile_path
    ALLOCATE(structOut%carbon%ptr2d,     source=struct1%carbon%ptr2d)
    structOut%nitrogen%owner_model_id    = struct1%nitrogen%owner_model_id
    structOut%nitrogen%owner_proc_id     = struct1%nitrogen%owner_proc_id
    structOut%nitrogen%owner_tile_path   = struct1%nitrogen%owner_tile_path
    ALLOCATE(structOut%nitrogen%ptr2d,   source=struct1%nitrogen%ptr2d)
    structOut%phosphorus%owner_model_id  = struct1%phosphorus%owner_model_id
    structOut%phosphorus%owner_proc_id   = struct1%phosphorus%owner_proc_id
    structOut%phosphorus%owner_tile_path = struct1%phosphorus%owner_tile_path
    ALLOCATE(structOut%phosphorus%ptr2d, source=struct1%phosphorus%ptr2d)
    structOut%carbon13%owner_model_id    = struct1%carbon13%owner_model_id
    structOut%carbon13%owner_proc_id     = struct1%carbon13%owner_proc_id
    structOut%carbon13%owner_tile_path   = struct1%carbon13%owner_tile_path
    ALLOCATE(structOut%carbon13%ptr2d,   source=struct1%carbon13%ptr2d)
    structOut%carbon14%owner_model_id    = struct1%carbon14%owner_model_id
    structOut%carbon14%owner_proc_id     = struct1%carbon14%owner_proc_id
    structOut%carbon14%owner_tile_path   = struct1%carbon14%owner_tile_path
    ALLOCATE(structOut%carbon14%ptr2d,   source=struct1%carbon14%ptr2d)
    structOut%nitrogen15%owner_model_id  = struct1%nitrogen15%owner_model_id
    structOut%nitrogen15%owner_proc_id   = struct1%nitrogen15%owner_proc_id
    structOut%nitrogen15%owner_tile_path = struct1%nitrogen15%owner_tile_path
    ALLOCATE(structOut%nitrogen15%ptr2d, source=struct1%nitrogen15%ptr2d)
    !! init element pools
    CALL structOut%leaf%Init_elem_var_elem_p(struct1%leaf)
    CALL structOut%fine_root%Init_elem_var_elem_p(struct1%fine_root)
    CALL structOut%coarse_root%Init_elem_var_elem_p(struct1%coarse_root)
    CALL structOut%sap_wood%Init_elem_var_elem_p(struct1%sap_wood)
    CALL structOut%heart_wood%Init_elem_var_elem_p(struct1%heart_wood)
    CALL structOut%labile%Init_elem_var_elem_p(struct1%labile)
    CALL structOut%reserve%Init_elem_var_elem_p(struct1%reserve)
    CALL structOut%fruit%Init_elem_var_elem_p(struct1%fruit)
    CALL structOut%seed_bed%Init_elem_var_elem_p(struct1%seed_bed)
  
    !> 2.0 run calculation
    !! the differentiation between ON_CHUNK/ON_DOMAIN is done at the level of variables (i.e., the op_ov routines of the jsb_var)
    !> 2.1 element variables of the tissue pool itself
    structOut%carbon          = field1d(:) * struct1%carbon    
    structOut%nitrogen        = field1d(:) * struct1%nitrogen  
    structOut%phosphorus      = field1d(:) * struct1%phosphorus
    structOut%carbon13        = field1d(:) * struct1%carbon13  
    structOut%carbon14        = field1d(:) * struct1%carbon14  
    structOut%nitrogen15      = field1d(:) * struct1%nitrogen15
    ! adding the pointer is this needed here too (not only with the jsb_var op_ov)
    structOut%carbon%ptr      => structOut%carbon%ptr2d(:,:)
    structOut%nitrogen%ptr    => structOut%nitrogen%ptr2d(:,:)
    structOut%phosphorus%ptr  => structOut%phosphorus%ptr2d(:,:)
    structOut%carbon13%ptr    => structOut%carbon13%ptr2d(:,:)
    structOut%carbon14%ptr    => structOut%carbon14%ptr2d(:,:)
    structOut%nitrogen15%ptr  => structOut%nitrogen15%ptr2d(:,:)
    !> 2.2 element pools
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
    TYPE(t_veg_tissues_pool),       TARGET      :: structOut
  
    CHARACTER(len=*), PARAMETER :: routine = modname//':Multiply_veg_tiss_p_with_field2d'
  
    !> 1.0 init temporary pool structOut
    !! init tissue pool with element variables
    CALL structOut%Init_tiss_pool_for_op_ov(name='op_ov_help', shortname='op_ov_help', &
      & unit=TRIM(struct1%element_unit), pool_suffix='veg_p') ! l_elements=.TRUE. by default
    structOut%carbon%owner_model_id      = struct1%carbon%owner_model_id
    structOut%carbon%owner_proc_id       = struct1%carbon%owner_proc_id
    structOut%carbon%owner_tile_path     = struct1%carbon%owner_tile_path
    ALLOCATE(structOut%carbon%ptr2d,     source=struct1%carbon%ptr2d)
    structOut%nitrogen%owner_model_id    = struct1%nitrogen%owner_model_id
    structOut%nitrogen%owner_proc_id     = struct1%nitrogen%owner_proc_id
    structOut%nitrogen%owner_tile_path   = struct1%nitrogen%owner_tile_path
    ALLOCATE(structOut%nitrogen%ptr2d,   source=struct1%nitrogen%ptr2d)
    structOut%phosphorus%owner_model_id  = struct1%phosphorus%owner_model_id
    structOut%phosphorus%owner_proc_id   = struct1%phosphorus%owner_proc_id
    structOut%phosphorus%owner_tile_path = struct1%phosphorus%owner_tile_path
    ALLOCATE(structOut%phosphorus%ptr2d, source=struct1%phosphorus%ptr2d)
    structOut%carbon13%owner_model_id    = struct1%carbon13%owner_model_id
    structOut%carbon13%owner_proc_id     = struct1%carbon13%owner_proc_id
    structOut%carbon13%owner_tile_path   = struct1%carbon13%owner_tile_path
    ALLOCATE(structOut%carbon13%ptr2d,   source=struct1%carbon13%ptr2d)
    structOut%carbon14%owner_model_id    = struct1%carbon14%owner_model_id
    structOut%carbon14%owner_proc_id     = struct1%carbon14%owner_proc_id
    structOut%carbon14%owner_tile_path   = struct1%carbon14%owner_tile_path
    ALLOCATE(structOut%carbon14%ptr2d,   source=struct1%carbon14%ptr2d)
    structOut%nitrogen15%owner_model_id  = struct1%nitrogen15%owner_model_id
    structOut%nitrogen15%owner_proc_id   = struct1%nitrogen15%owner_proc_id
    structOut%nitrogen15%owner_tile_path = struct1%nitrogen15%owner_tile_path
    ALLOCATE(structOut%nitrogen15%ptr2d, source=struct1%nitrogen15%ptr2d)
    !! init element pools
    CALL structOut%leaf%Init_elem_var_elem_p(struct1%leaf)
    CALL structOut%fine_root%Init_elem_var_elem_p(struct1%fine_root)
    CALL structOut%coarse_root%Init_elem_var_elem_p(struct1%coarse_root)
    CALL structOut%sap_wood%Init_elem_var_elem_p(struct1%sap_wood)
    CALL structOut%heart_wood%Init_elem_var_elem_p(struct1%heart_wood)
    CALL structOut%labile%Init_elem_var_elem_p(struct1%labile)
    CALL structOut%reserve%Init_elem_var_elem_p(struct1%reserve)
    CALL structOut%fruit%Init_elem_var_elem_p(struct1%fruit)
    CALL structOut%seed_bed%Init_elem_var_elem_p(struct1%seed_bed)
  
    !> 2.0 run calculation
    !! the differentiation between ON_CHUNK/ON_DOMAIN is done at the level of variables (i.e., the op_ov routines of the jsb_var)
    !> 2.1 element variables of the tissue pool itself
    structOut%carbon          = struct1%carbon      * field2d(:,:)
    structOut%nitrogen        = struct1%nitrogen    * field2d(:,:)
    structOut%phosphorus      = struct1%phosphorus  * field2d(:,:)
    structOut%carbon13        = struct1%carbon13    * field2d(:,:)
    structOut%carbon14        = struct1%carbon14    * field2d(:,:)
    structOut%nitrogen15      = struct1%nitrogen15  * field2d(:,:)
    ! adding the pointer is this needed here too (not only with the jsb_var op_ov)
    structOut%carbon%ptr      => structOut%carbon%ptr2d(:,:)
    structOut%nitrogen%ptr    => structOut%nitrogen%ptr2d(:,:)
    structOut%phosphorus%ptr  => structOut%phosphorus%ptr2d(:,:)
    structOut%carbon13%ptr    => structOut%carbon13%ptr2d(:,:)
    structOut%carbon14%ptr    => structOut%carbon14%ptr2d(:,:)
    structOut%nitrogen15%ptr  => structOut%nitrogen15%ptr2d(:,:)
    !> 2.2 element pools
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
    TYPE(t_veg_tissues_pool),       TARGET      :: structOut
  
    CHARACTER(len=*), PARAMETER :: routine = modname//':Multiply_field2d_with_veg_tiss_p'
  
    !> 1.0 init temporary pool structOut
    !! init tissue pool with element variables
    CALL structOut%Init_tiss_pool_for_op_ov(name='op_ov_help', shortname='op_ov_help', &
      & unit=TRIM(struct1%element_unit), pool_suffix='veg_p') ! l_elements=.TRUE. by default
    structOut%carbon%owner_model_id      = struct1%carbon%owner_model_id
    structOut%carbon%owner_proc_id       = struct1%carbon%owner_proc_id
    structOut%carbon%owner_tile_path     = struct1%carbon%owner_tile_path
    ALLOCATE(structOut%carbon%ptr2d,     source=struct1%carbon%ptr2d)
    structOut%nitrogen%owner_model_id    = struct1%nitrogen%owner_model_id
    structOut%nitrogen%owner_proc_id     = struct1%nitrogen%owner_proc_id
    structOut%nitrogen%owner_tile_path   = struct1%nitrogen%owner_tile_path
    ALLOCATE(structOut%nitrogen%ptr2d,   source=struct1%nitrogen%ptr2d)
    structOut%phosphorus%owner_model_id  = struct1%phosphorus%owner_model_id
    structOut%phosphorus%owner_proc_id   = struct1%phosphorus%owner_proc_id
    structOut%phosphorus%owner_tile_path = struct1%phosphorus%owner_tile_path
    ALLOCATE(structOut%phosphorus%ptr2d, source=struct1%phosphorus%ptr2d)
    structOut%carbon13%owner_model_id    = struct1%carbon13%owner_model_id
    structOut%carbon13%owner_proc_id     = struct1%carbon13%owner_proc_id
    structOut%carbon13%owner_tile_path   = struct1%carbon13%owner_tile_path
    ALLOCATE(structOut%carbon13%ptr2d,   source=struct1%carbon13%ptr2d)
    structOut%carbon14%owner_model_id    = struct1%carbon14%owner_model_id
    structOut%carbon14%owner_proc_id     = struct1%carbon14%owner_proc_id
    structOut%carbon14%owner_tile_path   = struct1%carbon14%owner_tile_path
    ALLOCATE(structOut%carbon14%ptr2d,   source=struct1%carbon14%ptr2d)
    structOut%nitrogen15%owner_model_id  = struct1%nitrogen15%owner_model_id
    structOut%nitrogen15%owner_proc_id   = struct1%nitrogen15%owner_proc_id
    structOut%nitrogen15%owner_tile_path = struct1%nitrogen15%owner_tile_path
    ALLOCATE(structOut%nitrogen15%ptr2d, source=struct1%nitrogen15%ptr2d)
    !! init element pools
    CALL structOut%leaf%Init_elem_var_elem_p(struct1%leaf)
    CALL structOut%fine_root%Init_elem_var_elem_p(struct1%fine_root)
    CALL structOut%coarse_root%Init_elem_var_elem_p(struct1%coarse_root)
    CALL structOut%sap_wood%Init_elem_var_elem_p(struct1%sap_wood)
    CALL structOut%heart_wood%Init_elem_var_elem_p(struct1%heart_wood)
    CALL structOut%labile%Init_elem_var_elem_p(struct1%labile)
    CALL structOut%reserve%Init_elem_var_elem_p(struct1%reserve)
    CALL structOut%fruit%Init_elem_var_elem_p(struct1%fruit)
    CALL structOut%seed_bed%Init_elem_var_elem_p(struct1%seed_bed)
  
    !> 2.0 run calculation
    !! the differentiation between ON_CHUNK/ON_DOMAIN is done at the level of variables (i.e., the op_ov routines of the jsb_var)
    !> 2.1 element variables of the tissue pool itself
    structOut%carbon          = field2d(:,:) * struct1%carbon    
    structOut%nitrogen        = field2d(:,:) * struct1%nitrogen  
    structOut%phosphorus      = field2d(:,:) * struct1%phosphorus
    structOut%carbon13        = field2d(:,:) * struct1%carbon13  
    structOut%carbon14        = field2d(:,:) * struct1%carbon14  
    structOut%nitrogen15      = field2d(:,:) * struct1%nitrogen15
    ! adding the pointer is this needed here too (not only with the jsb_var op_ov)
    structOut%carbon%ptr      => structOut%carbon%ptr2d(:,:)
    structOut%nitrogen%ptr    => structOut%nitrogen%ptr2d(:,:)
    structOut%phosphorus%ptr  => structOut%phosphorus%ptr2d(:,:)
    structOut%carbon13%ptr    => structOut%carbon13%ptr2d(:,:)
    structOut%carbon14%ptr    => structOut%carbon14%ptr2d(:,:)
    structOut%nitrogen15%ptr  => structOut%nitrogen15%ptr2d(:,:)
    !> 2.2 element pools
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
    TYPE(t_veg_tissues_pool),       TARGET      :: structOut
  
    CHARACTER(len=*), PARAMETER :: routine = modname//':Multiply_veg_tiss_p_with_veg_tiss_p'
  
    !> 1.0 init temporary pool structOut
    !! init tissue pool with element variables
    CALL structOut%Init_tiss_pool_for_op_ov(name='op_ov_help', shortname='op_ov_help', &
      & unit=TRIM(struct1%element_unit), pool_suffix='veg_p') ! l_elements=.TRUE. by default
    structOut%carbon%owner_model_id      = struct1%carbon%owner_model_id
    structOut%carbon%owner_proc_id       = struct1%carbon%owner_proc_id
    structOut%carbon%owner_tile_path     = struct1%carbon%owner_tile_path
    ALLOCATE(structOut%carbon%ptr2d,     source=struct1%carbon%ptr2d)
    structOut%nitrogen%owner_model_id    = struct1%nitrogen%owner_model_id
    structOut%nitrogen%owner_proc_id     = struct1%nitrogen%owner_proc_id
    structOut%nitrogen%owner_tile_path   = struct1%nitrogen%owner_tile_path
    ALLOCATE(structOut%nitrogen%ptr2d,   source=struct1%nitrogen%ptr2d)
    structOut%phosphorus%owner_model_id  = struct1%phosphorus%owner_model_id
    structOut%phosphorus%owner_proc_id   = struct1%phosphorus%owner_proc_id
    structOut%phosphorus%owner_tile_path = struct1%phosphorus%owner_tile_path
    ALLOCATE(structOut%phosphorus%ptr2d, source=struct1%phosphorus%ptr2d)
    structOut%carbon13%owner_model_id    = struct1%carbon13%owner_model_id
    structOut%carbon13%owner_proc_id     = struct1%carbon13%owner_proc_id
    structOut%carbon13%owner_tile_path   = struct1%carbon13%owner_tile_path
    ALLOCATE(structOut%carbon13%ptr2d,   source=struct1%carbon13%ptr2d)
    structOut%carbon14%owner_model_id    = struct1%carbon14%owner_model_id
    structOut%carbon14%owner_proc_id     = struct1%carbon14%owner_proc_id
    structOut%carbon14%owner_tile_path   = struct1%carbon14%owner_tile_path
    ALLOCATE(structOut%carbon14%ptr2d,   source=struct1%carbon14%ptr2d)
    structOut%nitrogen15%owner_model_id  = struct1%nitrogen15%owner_model_id
    structOut%nitrogen15%owner_proc_id   = struct1%nitrogen15%owner_proc_id
    structOut%nitrogen15%owner_tile_path = struct1%nitrogen15%owner_tile_path
    ALLOCATE(structOut%nitrogen15%ptr2d, source=struct1%nitrogen15%ptr2d)
    !! init element pools
    CALL structOut%leaf%Init_elem_var_elem_p(struct1%leaf)
    CALL structOut%fine_root%Init_elem_var_elem_p(struct1%fine_root)
    CALL structOut%coarse_root%Init_elem_var_elem_p(struct1%coarse_root)
    CALL structOut%sap_wood%Init_elem_var_elem_p(struct1%sap_wood)
    CALL structOut%heart_wood%Init_elem_var_elem_p(struct1%heart_wood)
    CALL structOut%labile%Init_elem_var_elem_p(struct1%labile)
    CALL structOut%reserve%Init_elem_var_elem_p(struct1%reserve)
    CALL structOut%fruit%Init_elem_var_elem_p(struct1%fruit)
    CALL structOut%seed_bed%Init_elem_var_elem_p(struct1%seed_bed)
  
    !> 2.0 run calculation
    !! the differentiation between ON_CHUNK/ON_DOMAIN is done at the level of variables (i.e., the op_ov routines of the jsb_var)
    !> 2.1 element variables of the tissue pool itself
    structOut%carbon          = struct1%carbon      * struct2%carbon    
    structOut%nitrogen        = struct1%nitrogen    * struct2%nitrogen  
    structOut%phosphorus      = struct1%phosphorus  * struct2%phosphorus
    structOut%carbon13        = struct1%carbon13    * struct2%carbon13  
    structOut%carbon14        = struct1%carbon14    * struct2%carbon14  
    structOut%nitrogen15      = struct1%nitrogen15  * struct2%nitrogen15
    ! adding the pointer is this needed here too (not only with the jsb_var op_ov)
    structOut%carbon%ptr      => structOut%carbon%ptr2d(:,:)
    structOut%nitrogen%ptr    => structOut%nitrogen%ptr2d(:,:)
    structOut%phosphorus%ptr  => structOut%phosphorus%ptr2d(:,:)
    structOut%carbon13%ptr    => structOut%carbon13%ptr2d(:,:)
    structOut%carbon14%ptr    => structOut%carbon14%ptr2d(:,:)
    structOut%nitrogen15%ptr  => structOut%nitrogen15%ptr2d(:,:)
    !> 2.2 element pools
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
    TYPE(t_veg_tissues_pool),       TARGET      :: structOut
  
    CHARACTER(len=*), PARAMETER :: routine = modname//':Div_veg_tiss_p_by_scalar'
  
    !> 1.0 init temporary pool structOut
    !! init tissue pool with element variables
    CALL structOut%Init_tiss_pool_for_op_ov(name='op_ov_help', shortname='op_ov_help', &
      & unit=TRIM(struct1%element_unit), pool_suffix='veg_p') ! l_elements=.TRUE. by default
    structOut%carbon%owner_model_id      = struct1%carbon%owner_model_id
    structOut%carbon%owner_proc_id       = struct1%carbon%owner_proc_id
    structOut%carbon%owner_tile_path     = struct1%carbon%owner_tile_path
    ALLOCATE(structOut%carbon%ptr2d,     source=struct1%carbon%ptr2d)
    structOut%nitrogen%owner_model_id    = struct1%nitrogen%owner_model_id
    structOut%nitrogen%owner_proc_id     = struct1%nitrogen%owner_proc_id
    structOut%nitrogen%owner_tile_path   = struct1%nitrogen%owner_tile_path
    ALLOCATE(structOut%nitrogen%ptr2d,   source=struct1%nitrogen%ptr2d)
    structOut%phosphorus%owner_model_id  = struct1%phosphorus%owner_model_id
    structOut%phosphorus%owner_proc_id   = struct1%phosphorus%owner_proc_id
    structOut%phosphorus%owner_tile_path = struct1%phosphorus%owner_tile_path
    ALLOCATE(structOut%phosphorus%ptr2d, source=struct1%phosphorus%ptr2d)
    structOut%carbon13%owner_model_id    = struct1%carbon13%owner_model_id
    structOut%carbon13%owner_proc_id     = struct1%carbon13%owner_proc_id
    structOut%carbon13%owner_tile_path   = struct1%carbon13%owner_tile_path
    ALLOCATE(structOut%carbon13%ptr2d,   source=struct1%carbon13%ptr2d)
    structOut%carbon14%owner_model_id    = struct1%carbon14%owner_model_id
    structOut%carbon14%owner_proc_id     = struct1%carbon14%owner_proc_id
    structOut%carbon14%owner_tile_path   = struct1%carbon14%owner_tile_path
    ALLOCATE(structOut%carbon14%ptr2d,   source=struct1%carbon14%ptr2d)
    structOut%nitrogen15%owner_model_id  = struct1%nitrogen15%owner_model_id
    structOut%nitrogen15%owner_proc_id   = struct1%nitrogen15%owner_proc_id
    structOut%nitrogen15%owner_tile_path = struct1%nitrogen15%owner_tile_path
    ALLOCATE(structOut%nitrogen15%ptr2d, source=struct1%nitrogen15%ptr2d)
    !! init element pools
    CALL structOut%leaf%Init_elem_var_elem_p(struct1%leaf)
    CALL structOut%fine_root%Init_elem_var_elem_p(struct1%fine_root)
    CALL structOut%coarse_root%Init_elem_var_elem_p(struct1%coarse_root)
    CALL structOut%sap_wood%Init_elem_var_elem_p(struct1%sap_wood)
    CALL structOut%heart_wood%Init_elem_var_elem_p(struct1%heart_wood)
    CALL structOut%labile%Init_elem_var_elem_p(struct1%labile)
    CALL structOut%reserve%Init_elem_var_elem_p(struct1%reserve)
    CALL structOut%fruit%Init_elem_var_elem_p(struct1%fruit)
    CALL structOut%seed_bed%Init_elem_var_elem_p(struct1%seed_bed)
  
    !> 2.0 run calculation
    !! the differentiation between ON_CHUNK/ON_DOMAIN is done at the level of variables (i.e., the op_ov routines of the jsb_var)
    !> 2.1 element variables of the tissue pool itself
    structOut%carbon          = struct1%carbon      / scalar
    structOut%nitrogen        = struct1%nitrogen    / scalar
    structOut%phosphorus      = struct1%phosphorus  / scalar
    structOut%carbon13        = struct1%carbon13    / scalar
    structOut%carbon14        = struct1%carbon14    / scalar
    structOut%nitrogen15      = struct1%nitrogen15  / scalar
    ! adding the pointer is this needed here too (not only with the jsb_var op_ov)
    structOut%carbon%ptr      => structOut%carbon%ptr2d(:,:)
    structOut%nitrogen%ptr    => structOut%nitrogen%ptr2d(:,:)
    structOut%phosphorus%ptr  => structOut%phosphorus%ptr2d(:,:)
    structOut%carbon13%ptr    => structOut%carbon13%ptr2d(:,:)
    structOut%carbon14%ptr    => structOut%carbon14%ptr2d(:,:)
    structOut%nitrogen15%ptr  => structOut%nitrogen15%ptr2d(:,:)
    !> 2.2 element pools
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
    TYPE(t_veg_tissues_pool),       TARGET      :: structOut
  
    CHARACTER(len=*), PARAMETER :: routine = modname//':Div_veg_tiss_p_by_field1d'
  
    !> 1.0 init temporary pool structOut
    !! init tissue pool with element variables
    CALL structOut%Init_tiss_pool_for_op_ov(name='op_ov_help', shortname='op_ov_help', &
      & unit=TRIM(struct1%element_unit), pool_suffix='veg_p') ! l_elements=.TRUE. by default
    structOut%carbon%owner_model_id      = struct1%carbon%owner_model_id
    structOut%carbon%owner_proc_id       = struct1%carbon%owner_proc_id
    structOut%carbon%owner_tile_path     = struct1%carbon%owner_tile_path
    ALLOCATE(structOut%carbon%ptr2d,     source=struct1%carbon%ptr2d)
    structOut%nitrogen%owner_model_id    = struct1%nitrogen%owner_model_id
    structOut%nitrogen%owner_proc_id     = struct1%nitrogen%owner_proc_id
    structOut%nitrogen%owner_tile_path   = struct1%nitrogen%owner_tile_path
    ALLOCATE(structOut%nitrogen%ptr2d,   source=struct1%nitrogen%ptr2d)
    structOut%phosphorus%owner_model_id  = struct1%phosphorus%owner_model_id
    structOut%phosphorus%owner_proc_id   = struct1%phosphorus%owner_proc_id
    structOut%phosphorus%owner_tile_path = struct1%phosphorus%owner_tile_path
    ALLOCATE(structOut%phosphorus%ptr2d, source=struct1%phosphorus%ptr2d)
    structOut%carbon13%owner_model_id    = struct1%carbon13%owner_model_id
    structOut%carbon13%owner_proc_id     = struct1%carbon13%owner_proc_id
    structOut%carbon13%owner_tile_path   = struct1%carbon13%owner_tile_path
    ALLOCATE(structOut%carbon13%ptr2d,   source=struct1%carbon13%ptr2d)
    structOut%carbon14%owner_model_id    = struct1%carbon14%owner_model_id
    structOut%carbon14%owner_proc_id     = struct1%carbon14%owner_proc_id
    structOut%carbon14%owner_tile_path   = struct1%carbon14%owner_tile_path
    ALLOCATE(structOut%carbon14%ptr2d,   source=struct1%carbon14%ptr2d)
    structOut%nitrogen15%owner_model_id  = struct1%nitrogen15%owner_model_id
    structOut%nitrogen15%owner_proc_id   = struct1%nitrogen15%owner_proc_id
    structOut%nitrogen15%owner_tile_path = struct1%nitrogen15%owner_tile_path
    ALLOCATE(structOut%nitrogen15%ptr2d, source=struct1%nitrogen15%ptr2d)
    !! init element pools
    CALL structOut%leaf%Init_elem_var_elem_p(struct1%leaf)
    CALL structOut%fine_root%Init_elem_var_elem_p(struct1%fine_root)
    CALL structOut%coarse_root%Init_elem_var_elem_p(struct1%coarse_root)
    CALL structOut%sap_wood%Init_elem_var_elem_p(struct1%sap_wood)
    CALL structOut%heart_wood%Init_elem_var_elem_p(struct1%heart_wood)
    CALL structOut%labile%Init_elem_var_elem_p(struct1%labile)
    CALL structOut%reserve%Init_elem_var_elem_p(struct1%reserve)
    CALL structOut%fruit%Init_elem_var_elem_p(struct1%fruit)
    CALL structOut%seed_bed%Init_elem_var_elem_p(struct1%seed_bed)
  
    !> 2.0 run calculation
    !! the differentiation between ON_CHUNK/ON_DOMAIN is done at the level of variables (i.e., the op_ov routines of the jsb_var)
    !> 2.1 element variables of the tissue pool itself
    structOut%carbon          = struct1%carbon      / field1d(:)
    structOut%nitrogen        = struct1%nitrogen    / field1d(:)
    structOut%phosphorus      = struct1%phosphorus  / field1d(:)
    structOut%carbon13        = struct1%carbon13    / field1d(:)
    structOut%carbon14        = struct1%carbon14    / field1d(:)
    structOut%nitrogen15      = struct1%nitrogen15  / field1d(:)
    ! adding the pointer is this needed here too (not only with the jsb_var op_ov)
    structOut%carbon%ptr      => structOut%carbon%ptr2d(:,:)
    structOut%nitrogen%ptr    => structOut%nitrogen%ptr2d(:,:)
    structOut%phosphorus%ptr  => structOut%phosphorus%ptr2d(:,:)
    structOut%carbon13%ptr    => structOut%carbon13%ptr2d(:,:)
    structOut%carbon14%ptr    => structOut%carbon14%ptr2d(:,:)
    structOut%nitrogen15%ptr  => structOut%nitrogen15%ptr2d(:,:)
    !> 2.2 element pools
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
    TYPE(t_veg_tissues_pool),       TARGET      :: structOut
  
    CHARACTER(len=*), PARAMETER :: routine = modname//':Div_veg_tiss_p_by_field2d'
  
    !> 1.0 init temporary pool structOut
    !! init tissue pool with element variables
    CALL structOut%Init_tiss_pool_for_op_ov(name='op_ov_help', shortname='op_ov_help', &
      & unit=TRIM(struct1%element_unit), pool_suffix='veg_p') ! l_elements=.TRUE. by default
    structOut%carbon%owner_model_id      = struct1%carbon%owner_model_id
    structOut%carbon%owner_proc_id       = struct1%carbon%owner_proc_id
    structOut%carbon%owner_tile_path     = struct1%carbon%owner_tile_path
    ALLOCATE(structOut%carbon%ptr2d,     source=struct1%carbon%ptr2d)
    structOut%nitrogen%owner_model_id    = struct1%nitrogen%owner_model_id
    structOut%nitrogen%owner_proc_id     = struct1%nitrogen%owner_proc_id
    structOut%nitrogen%owner_tile_path   = struct1%nitrogen%owner_tile_path
    ALLOCATE(structOut%nitrogen%ptr2d,   source=struct1%nitrogen%ptr2d)
    structOut%phosphorus%owner_model_id  = struct1%phosphorus%owner_model_id
    structOut%phosphorus%owner_proc_id   = struct1%phosphorus%owner_proc_id
    structOut%phosphorus%owner_tile_path = struct1%phosphorus%owner_tile_path
    ALLOCATE(structOut%phosphorus%ptr2d, source=struct1%phosphorus%ptr2d)
    structOut%carbon13%owner_model_id    = struct1%carbon13%owner_model_id
    structOut%carbon13%owner_proc_id     = struct1%carbon13%owner_proc_id
    structOut%carbon13%owner_tile_path   = struct1%carbon13%owner_tile_path
    ALLOCATE(structOut%carbon13%ptr2d,   source=struct1%carbon13%ptr2d)
    structOut%carbon14%owner_model_id    = struct1%carbon14%owner_model_id
    structOut%carbon14%owner_proc_id     = struct1%carbon14%owner_proc_id
    structOut%carbon14%owner_tile_path   = struct1%carbon14%owner_tile_path
    ALLOCATE(structOut%carbon14%ptr2d,   source=struct1%carbon14%ptr2d)
    structOut%nitrogen15%owner_model_id  = struct1%nitrogen15%owner_model_id
    structOut%nitrogen15%owner_proc_id   = struct1%nitrogen15%owner_proc_id
    structOut%nitrogen15%owner_tile_path = struct1%nitrogen15%owner_tile_path
    ALLOCATE(structOut%nitrogen15%ptr2d, source=struct1%nitrogen15%ptr2d)
    !! init element pools
    CALL structOut%leaf%Init_elem_var_elem_p(struct1%leaf)
    CALL structOut%fine_root%Init_elem_var_elem_p(struct1%fine_root)
    CALL structOut%coarse_root%Init_elem_var_elem_p(struct1%coarse_root)
    CALL structOut%sap_wood%Init_elem_var_elem_p(struct1%sap_wood)
    CALL structOut%heart_wood%Init_elem_var_elem_p(struct1%heart_wood)
    CALL structOut%labile%Init_elem_var_elem_p(struct1%labile)
    CALL structOut%reserve%Init_elem_var_elem_p(struct1%reserve)
    CALL structOut%fruit%Init_elem_var_elem_p(struct1%fruit)
    CALL structOut%seed_bed%Init_elem_var_elem_p(struct1%seed_bed)
  
    !> 2.0 run calculation
    !! the differentiation between ON_CHUNK/ON_DOMAIN is done at the level of variables (i.e., the op_ov routines of the jsb_var)
    !> 2.1 element variables of the tissue pool itself
    structOut%carbon          = struct1%carbon      / field2d(:,:)
    structOut%nitrogen        = struct1%nitrogen    / field2d(:,:)
    structOut%phosphorus      = struct1%phosphorus  / field2d(:,:)
    structOut%carbon13        = struct1%carbon13    / field2d(:,:)
    structOut%carbon14        = struct1%carbon14    / field2d(:,:)
    structOut%nitrogen15      = struct1%nitrogen15  / field2d(:,:)
    ! adding the pointer is this needed here too (not only with the jsb_var op_ov)
    structOut%carbon%ptr      => structOut%carbon%ptr2d(:,:)
    structOut%nitrogen%ptr    => structOut%nitrogen%ptr2d(:,:)
    structOut%phosphorus%ptr  => structOut%phosphorus%ptr2d(:,:)
    structOut%carbon13%ptr    => structOut%carbon13%ptr2d(:,:)
    structOut%carbon14%ptr    => structOut%carbon14%ptr2d(:,:)
    structOut%nitrogen15%ptr  => structOut%nitrogen15%ptr2d(:,:)
    !> 2.2 element pools
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
    TYPE(t_veg_tissues_pool),       TARGET      :: structOut
  
    CHARACTER(len=*), PARAMETER :: routine = modname//':Div_veg_tiss_p_by_veg_tiss_p'
  
    !> 1.0 init temporary pool structOut
    !! init tissue pool with element variables
    CALL structOut%Init_tiss_pool_for_op_ov(name='op_ov_help', shortname='op_ov_help', &
      & unit=TRIM(struct1%element_unit), pool_suffix='veg_p') ! l_elements=.TRUE. by default
    structOut%carbon%owner_model_id      = struct1%carbon%owner_model_id
    structOut%carbon%owner_proc_id       = struct1%carbon%owner_proc_id
    structOut%carbon%owner_tile_path     = struct1%carbon%owner_tile_path
    ALLOCATE(structOut%carbon%ptr2d,     source=struct1%carbon%ptr2d)
    structOut%nitrogen%owner_model_id    = struct1%nitrogen%owner_model_id
    structOut%nitrogen%owner_proc_id     = struct1%nitrogen%owner_proc_id
    structOut%nitrogen%owner_tile_path   = struct1%nitrogen%owner_tile_path
    ALLOCATE(structOut%nitrogen%ptr2d,   source=struct1%nitrogen%ptr2d)
    structOut%phosphorus%owner_model_id  = struct1%phosphorus%owner_model_id
    structOut%phosphorus%owner_proc_id   = struct1%phosphorus%owner_proc_id
    structOut%phosphorus%owner_tile_path = struct1%phosphorus%owner_tile_path
    ALLOCATE(structOut%phosphorus%ptr2d, source=struct1%phosphorus%ptr2d)
    structOut%carbon13%owner_model_id    = struct1%carbon13%owner_model_id
    structOut%carbon13%owner_proc_id     = struct1%carbon13%owner_proc_id
    structOut%carbon13%owner_tile_path   = struct1%carbon13%owner_tile_path
    ALLOCATE(structOut%carbon13%ptr2d,   source=struct1%carbon13%ptr2d)
    structOut%carbon14%owner_model_id    = struct1%carbon14%owner_model_id
    structOut%carbon14%owner_proc_id     = struct1%carbon14%owner_proc_id
    structOut%carbon14%owner_tile_path   = struct1%carbon14%owner_tile_path
    ALLOCATE(structOut%carbon14%ptr2d,   source=struct1%carbon14%ptr2d)
    structOut%nitrogen15%owner_model_id  = struct1%nitrogen15%owner_model_id
    structOut%nitrogen15%owner_proc_id   = struct1%nitrogen15%owner_proc_id
    structOut%nitrogen15%owner_tile_path = struct1%nitrogen15%owner_tile_path
    ALLOCATE(structOut%nitrogen15%ptr2d, source=struct1%nitrogen15%ptr2d)
    !! init element pools
    CALL structOut%leaf%Init_elem_var_elem_p(struct1%leaf)
    CALL structOut%fine_root%Init_elem_var_elem_p(struct1%fine_root)
    CALL structOut%coarse_root%Init_elem_var_elem_p(struct1%coarse_root)
    CALL structOut%sap_wood%Init_elem_var_elem_p(struct1%sap_wood)
    CALL structOut%heart_wood%Init_elem_var_elem_p(struct1%heart_wood)
    CALL structOut%labile%Init_elem_var_elem_p(struct1%labile)
    CALL structOut%reserve%Init_elem_var_elem_p(struct1%reserve)
    CALL structOut%fruit%Init_elem_var_elem_p(struct1%fruit)
    CALL structOut%seed_bed%Init_elem_var_elem_p(struct1%seed_bed)
  
    !> 2.0 run calculation
    !! the differentiation between ON_CHUNK/ON_DOMAIN is done at the level of variables (i.e., the op_ov routines of the jsb_var)
    !> 2.1 element variables of the tissue pool itself
    structOut%carbon          = struct1%carbon      / struct2%carbon    
    structOut%nitrogen        = struct1%nitrogen    / struct2%nitrogen  
    structOut%phosphorus      = struct1%phosphorus  / struct2%phosphorus
    structOut%carbon13        = struct1%carbon13    / struct2%carbon13  
    structOut%carbon14        = struct1%carbon14    / struct2%carbon14  
    structOut%nitrogen15      = struct1%nitrogen15  / struct2%nitrogen15
    ! adding the pointer is this needed here too (not only with the jsb_var op_ov)
    structOut%carbon%ptr      => structOut%carbon%ptr2d(:,:)
    structOut%nitrogen%ptr    => structOut%nitrogen%ptr2d(:,:)
    structOut%phosphorus%ptr  => structOut%phosphorus%ptr2d(:,:)
    structOut%carbon13%ptr    => structOut%carbon13%ptr2d(:,:)
    structOut%carbon14%ptr    => structOut%carbon14%ptr2d(:,:)
    structOut%nitrogen15%ptr  => structOut%nitrogen15%ptr2d(:,:)
    !> 2.2 element pools
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
