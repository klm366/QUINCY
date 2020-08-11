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
  USE mo_exception,              ONLY: finish
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
    TYPE(t_jsb_var_real2d), POINTER :: carbon
    TYPE(t_jsb_var_real2d), POINTER :: nitrogen
    TYPE(t_jsb_var_real2d), POINTER :: phosphorus
    TYPE(t_jsb_var_real2d), POINTER :: carbon13
    TYPE(t_jsb_var_real2d), POINTER :: carbon14
    TYPE(t_jsb_var_real2d), POINTER :: nitrogen15
  CONTAINS
    PRIVATE
    ! infrastructure
    PROCEDURE :: Init                   => t_veg_elements_pool_init
    PROCEDURE :: Get_element_name_by_id => t_veg_elements_pool_get_element_name_by_id
    ! operator overload
    ! + addition
    ! PROCEDURE         :: Add_element_a                      ! 
    ! PROCEDURE         :: Add_element_a_and_scalar           ! 
    ! PROCEDURE         :: Add_veg_elem_p_to_veg_elem_p       ! pool%elements + pool%elements
    ! GENERIC, PUBLIC   :: OPERATOR(+) => Add_veg_elem_p_to_veg_elem_p
    ! = assign
    PROCEDURE         :: Assign_scalar_to_veg_elem_p          ! pool%elements = scalar
    PROCEDURE         :: Assign_field1d_to_veg_elem_p         ! pool%elements = var(:)
    PROCEDURE         :: Assign_veg_elem_p_hlp_to_veg_elem_p  ! pool%elements = t_veg_elements_pool_hlp
    PROCEDURE         :: Assign_veg_elem_p_to_veg_elem_p      ! pool%elements = pool%elements  ! only numeric values of elements are assigned
    GENERIC, PUBLIC   :: ASSIGNMENT(=) => Assign_scalar_to_veg_elem_p, &
                                          Assign_field1d_to_veg_elem_p, &
                                          Assign_veg_elem_p_hlp_to_veg_elem_p, &
                                          Assign_veg_elem_p_to_veg_elem_p
  END TYPE t_veg_elements_pool

  !-----------------------------------------------------------------------------------------------------
  !> Simple veg_elements_pool TYPE for temporary objects with operator overload
  !! 
  !! 2D is needed in the context of ON_DOMAIN
  !-----------------------------------------------------------------------------------------------------
  TYPE t_veg_elements_pool_hlp
    REAL(wp), ALLOCATABLE :: carbon(:,:)      ! dimension of variable in memory is var(ics:ice,iblk)
    REAL(wp), ALLOCATABLE :: nitrogen(:,:)
    REAL(wp), ALLOCATABLE :: phosphorus(:,:)
    REAL(wp), ALLOCATABLE :: carbon13(:,:)
    REAL(wp), ALLOCATABLE :: carbon14(:,:)
    REAL(wp), ALLOCATABLE :: nitrogen15(:,:)
  END TYPE t_veg_elements_pool_hlp

  !-----------------------------------------------------------------------------------------------------
  !> Type definition for t_veg_tissues_pool
  !! vegetation pool of tissues
  !-----------------------------------------------------------------------------------------------------
  TYPE, EXTENDS(t_veg_elements_pool) :: t_veg_tissues_pool
    TYPE(t_veg_elements_pool), POINTER :: leaf             !< elements of leaf
    TYPE(t_veg_elements_pool), POINTER :: fine_root        !< elements of fine_root
    TYPE(t_veg_elements_pool), POINTER :: coarse_root      !< elements of coarse_root
    TYPE(t_veg_elements_pool), POINTER :: sap_wood         !< elements of sap_wood
    TYPE(t_veg_elements_pool), POINTER :: heart_wood       !< elements of heart_wood
    TYPE(t_veg_elements_pool), POINTER :: labile           !< elements of labile
    TYPE(t_veg_elements_pool), POINTER :: reserve          !< elements of reserve
    TYPE(t_veg_elements_pool), POINTER :: fruit            !< elements of fruit
    TYPE(t_veg_elements_pool), POINTER :: seed_bed         !< elements of seed_bed
  END TYPE t_veg_tissues_pool

  !-----------------------------------------------------------------------------------------------------
  !> Type definition for t_veg_main_pool
  !! pool-structure containing pools and fluxes
  !! 1st pool in line
  !! added directly into the t_veg_memory as "CLASS(t_jsb_pool), POINTER :: pools"
  !-----------------------------------------------------------------------------------------------------
  TYPE, EXTENDS(t_veg_elements_pool) :: t_veg_main_pool
    ! pools
    TYPE(t_veg_tissues_pool),   POINTER :: veg_pool
    TYPE(t_veg_elements_pool),  POINTER :: veg_total_biomass
    ! fluxes
    TYPE(t_veg_tissues_pool),   POINTER :: veg_growth    
    TYPE(t_veg_tissues_pool),   POINTER :: veg_litterfall
    TYPE(t_veg_tissues_pool),   POINTER :: veg_falloc    
    TYPE(t_veg_tissues_pool),   POINTER :: veg_total_flux
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
  !> create and init an object of the VEG_ tissues pool "t_veg_tissues_pool"
  !! init and add all its objects of the "t_veg_elements_pool" with "pool_object%Init()"
  !! "object => Create_veg_tissues_pool()"
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
  !> operator overload and assignment routines for pools
  !!  TYPE t_veg_elements_pool
  !! 
  !--------------------------------------------------------------------------------------------------------------------------------
  
  !--------------------------------------------------------------------------------------------------------------------------------
  ! assign
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

  SUBROUTINE Assign_veg_elem_p_hlp_to_veg_elem_p(lhs, rhs)
  
    CLASS(t_veg_elements_pool),     INTENT(inout) :: lhs
    TYPE(t_veg_elements_pool_hlp ), INTENT(in)    :: rhs

    CHARACTER(len=*), PARAMETER :: routine = modname//':Assign_veg_elem_p_hlp_to_veg_elem_p'

    lhs%carbon     = rhs%carbon(:,:)
    lhs%nitrogen   = rhs%nitrogen(:,:)
    lhs%phosphorus = rhs%phosphorus(:,:)
    lhs%carbon13   = rhs%carbon13(:,:)
    lhs%carbon14   = rhs%carbon14(:,:)
    lhs%nitrogen15 = rhs%nitrogen15(:,:)
  
  END SUBROUTINE Assign_veg_elem_p_hlp_to_veg_elem_p

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
  ! FUNCTION Add_veg_elem_p_to_veg_elem_p(struct1, struct2) RESULT(res)
  ! 
  !   USE mo_jsb_subset,        ONLY: t_subset, ON_DOMAIN, ON_CHUNK
  ! 
  !   CLASS(t_veg_elements_pool),     INTENT(in)  :: struct1, struct2
  !   TYPE(t_veg_elements_pool_hlp)               :: res
  ! 
  !   TYPE(t_subset) :: sub
  !   INTEGER        :: ics, ice, iblk, nb, nc
  ! 
  !   CHARACTER(len=*), PARAMETER :: routine = modname//':Add_veg_elem_p_to_veg_elem_p'
  ! 
  !   sub = struct1%carbon%Get_subset()
  !   ics = sub%ics ; ice = sub%ice ; iblk = sub%iblk ; nb = sub%nb ; nc = sub%nc
  ! 
  !   IF (sub%type == ON_CHUNK) THEN
  !     ALLOCATE(res%carbon(nc,1))
  !     ALLOCATE(res%nitrogen(nc,1))
  !     ALLOCATE(res%phosphorus(nc,1))
  !     ALLOCATE(res%carbon13(nc,1))
  !     ALLOCATE(res%carbon14(nc,1))
  !     ALLOCATE(res%nitrogen15(nc,1))
  !     res%carbon(:,:)     = struct1%carbon      + struct2%carbon
  !     res%nitrogen(:,:)   = struct1%nitrogen    + struct2%nitrogen  
  !     res%phosphorus(:,:) = struct1%phosphorus  + struct2%phosphorus
  !     res%carbon13(:,:)   = struct1%carbon13    + struct2%carbon13  
  !     res%carbon14(:,:)   = struct1%carbon14    + struct2%carbon14  
  !     res%nitrogen15(:,:) = struct1%nitrogen15  + struct2%nitrogen15
  !   ELSE
  !     ALLOCATE(res%carbon(nc,nb))
  !     ALLOCATE(res%nitrogen(nc,nb))
  !     ALLOCATE(res%phosphorus(nc,nb))
  !     ALLOCATE(res%carbon13(nc,nb))
  !     ALLOCATE(res%carbon14(nc,nb))
  !     ALLOCATE(res%nitrogen15(nc,nb))
  !     res%carbon(:,:)     = struct1%carbon      + struct2%carbon
  !     res%nitrogen(:,:)   = struct1%nitrogen    + struct2%nitrogen  
  !     res%phosphorus(:,:) = struct1%phosphorus  + struct2%phosphorus
  !     res%carbon13(:,:)   = struct1%carbon13    + struct2%carbon13  
  !     res%carbon14(:,:)   = struct1%carbon14    + struct2%carbon14  
  !     res%nitrogen15(:,:) = struct1%nitrogen15  + struct2%nitrogen15
  !   ENDIF
  ! 
  ! END FUNCTION Add_veg_elem_p_to_veg_elem_p  
    

  !--------------------------------------------------------------------------------------------------------------------------------
  !> operator overload and assignment routines for pools
  !!  TYPE t_veg_tissues_pool
  !! 
  !--------------------------------------------------------------------------------------------------------------------------------


END MODULE mo_veg_memory_class
