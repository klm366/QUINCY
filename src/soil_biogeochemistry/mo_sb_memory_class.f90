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
!! @brief Contains the memory class for the soil_biogeochemistry processes
!!
!! @author
!!  Jan Engel, MPI-BGC Jena (based on the JSBACH4 version from Rainer Schneck, MPI-M Hamburg)
!!
!! @par Revision History
!! First version by Jan Engel (20170103)
!!
MODULE mo_sb_memory_class

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
  PUBLIC :: t_sb_memory, max_no_of_vars
  PUBLIC :: t_sb_main_pool, t_sb_components_pool, t_sb_elements_pool

  INTEGER, PARAMETER :: max_no_of_vars      = 480  !> TYPE t_jsb_var | 'normal' var: ~ 200 | 4 pools: 240 var | 2 tissues: 12

  !-----------------------------------------------------------------------------------------------------
  !> Type definition for sb_elements_pool
  !! pool with only variables
  !! last pool in line
  !-----------------------------------------------------------------------------------------------------
  TYPE, EXTENDS(t_jsb_pool) :: t_sb_elements_pool
    TYPE(t_jsb_var_real3d), POINTER :: carbon
    TYPE(t_jsb_var_real3d), POINTER :: nitrogen
    TYPE(t_jsb_var_real3d), POINTER :: phosphorus
    TYPE(t_jsb_var_real3d), POINTER :: carbon13
    TYPE(t_jsb_var_real3d), POINTER :: carbon14
    TYPE(t_jsb_var_real3d), POINTER :: nitrogen15
  CONTAINS
    !! infrastructure
    PROCEDURE :: Init                   => t_sb_elements_pool_init
    PROCEDURE :: Get_element_name_by_id => t_sb_elements_pool_get_element_name_by_id
    PROCEDURE :: Init_elem_var_elem_p   => Init_element_var_sb_elements_pool            ! used in t_sb_components_pool operator overload routines to init temporary objects
    !! operator overload
    !! [def: field2d(nc,vgrid) | field3d(nc,vgrid,nblks)]
    ! +
    PROCEDURE                 :: Add_sb_elem_p_and_scalar            ! elements_pool + scalar        | ON_CHUNK & ON_DOMAIN
    PROCEDURE, PASS(struct1)  :: Add_scalar_and_sb_elem_p            ! scalar        + elements_pool | ON_CHUNK & ON_DOMAIN
    PROCEDURE                 :: Add_sb_elem_p_and_field2d           ! elements_pool + field2d       | ON_CHUNK
    PROCEDURE, PASS(struct1)  :: Add_field2d_and_sb_elem_p           ! field2d       + elements_pool | ON_CHUNK
    PROCEDURE                 :: Add_sb_elem_p_and_field3d           ! elements_pool + field3d       | ON_CHUNK & ON_DOMAIN
    PROCEDURE                 :: Add_sb_elem_p_and_sb_elem_p         ! elements_pool + elements_pool | ON_CHUNK & ON_DOMAIN
    GENERIC                   :: OPERATOR(+) => Add_sb_elem_p_and_scalar,  Add_scalar_and_sb_elem_p, &
      &                                         Add_sb_elem_p_and_field2d, Add_field2d_and_sb_elem_p, &
      &                                         Add_sb_elem_p_and_field3d, &
      &                                         Add_sb_elem_p_and_sb_elem_p
    ! -
    PROCEDURE                 :: Subtract_sb_elem_p_minus_scalar     ! elements_pool - scalar        | ON_CHUNK & ON_DOMAIN
    PROCEDURE                 :: Subtract_sb_elem_p_minus_field2d    ! elements_pool - field2d       | ON_CHUNK
    PROCEDURE                 :: Subtract_sb_elem_p_minus_field3d    ! elements_pool - field3d       | ON_CHUNK & ON_DOMAIN
    PROCEDURE                 :: Subtract_sb_elem_p_minus_sb_elem_p  ! elements_pool - elements_pool | ON_CHUNK & ON_DOMAIN
    GENERIC                   :: OPERATOR(-) => Subtract_sb_elem_p_minus_scalar, &
      &                                         Subtract_sb_elem_p_minus_field2d, &
      &                                         Subtract_sb_elem_p_minus_field3d, &
      &                                         Subtract_sb_elem_p_minus_sb_elem_p
    ! *
    PROCEDURE                 :: Multiply_sb_elem_p_with_scalar      ! elements_pool * scalar        | ON_CHUNK & ON_DOMAIN
    PROCEDURE, PASS(struct1)  :: Multiply_scalar_with_sb_elem_p      ! scalar        * elements_pool | ON_CHUNK & ON_DOMAIN
    PROCEDURE                 :: Multiply_sb_elem_p_with_field2d     ! elements_pool * field2d       | ON_CHUNK
    PROCEDURE, PASS(struct1)  :: Multiply_field2d_with_sb_elem_p     ! field2d       * elements_pool | ON_CHUNK
    PROCEDURE                 :: Multiply_sb_elem_p_with_field3d     ! elements_pool * field3d       | ON_CHUNK & ON_DOMAIN
    PROCEDURE, PASS(struct1)  :: Multiply_field3d_with_sb_elem_p     ! field3d       * elements_pool | ON_CHUNK & ON_DOMAIN
    PROCEDURE                 :: Multiply_sb_elem_p_with_sb_elem_p   ! elements_pool * elements_pool | ON_CHUNK & ON_DOMAIN
    GENERIC                   :: OPERATOR(*) => Multiply_sb_elem_p_with_scalar,  Multiply_scalar_with_sb_elem_p, &
      &                                         Multiply_sb_elem_p_with_field2d, Multiply_field2d_with_sb_elem_p, &
      &                                         Multiply_sb_elem_p_with_field3d, Multiply_field3d_with_sb_elem_p, &
      &                                         Multiply_sb_elem_p_with_sb_elem_p
    ! /
    PROCEDURE                 :: Div_sb_elem_p_by_scalar             ! elements_pool / scalar        | ON_CHUNK & ON_DOMAIN
    PROCEDURE                 :: Div_sb_elem_p_by_field2d            ! elements_pool / field2d       | ON_CHUNK
    PROCEDURE                 :: Div_sb_elem_p_by_field3d            ! elements_pool / field3d       | ON_CHUNK & ON_DOMAIN
    PROCEDURE                 :: Div_sb_elem_p_by_sb_elem_p          ! elements_pool / elements_pool | ON_CHUNK & ON_DOMAIN
    GENERIC                   :: OPERATOR(/) => Div_sb_elem_p_by_scalar, &
      &                                         Div_sb_elem_p_by_field2d, &
      &                                         Div_sb_elem_p_by_field3d, &
      &                                         Div_sb_elem_p_by_sb_elem_p
    ! =
    PROCEDURE                 :: Assign_scalar_to_sb_elem_p          ! elements_pool = scalar        | ON_CHUNK & ON_DOMAIN
    PROCEDURE                 :: Assign_field2d_to_sb_elem_p         ! elements_pool = field2d(:,:)  | ON_CHUNK
    PROCEDURE                 :: Assign_field3d_to_sb_elem_p         ! elements_pool = field3d(:,:,:)| ON_CHUNK & ON_DOMAIN
    PROCEDURE                 :: Assign_sb_elem_p_to_sb_elem_p       ! elements_pool = elements_pool | ON_CHUNK & ON_DOMAIN | only numeric values of elements are assigned
    GENERIC                   :: ASSIGNMENT(=) => Assign_scalar_to_sb_elem_p, &
      &                                           Assign_field2d_to_sb_elem_p, &
      &                                           Assign_field3d_to_sb_elem_p, &
      &                                           Assign_sb_elem_p_to_sb_elem_p
    !! finalize
    FINAL     :: finalize_sb_elements_pool
  END TYPE t_sb_elements_pool

  !-----------------------------------------------------------------------------------------------------
  !> Type definition for t_sb_components_pool
  !! soil_biogeochemistry pool of components
  !-----------------------------------------------------------------------------------------------------
  TYPE, EXTENDS(t_jsb_pool) :: t_sb_components_pool
    ! element variables
    TYPE(t_jsb_var_real3d), POINTER :: carbon
    TYPE(t_jsb_var_real3d), POINTER :: nitrogen
    TYPE(t_jsb_var_real3d), POINTER :: phosphorus
    TYPE(t_jsb_var_real3d), POINTER :: carbon13
    TYPE(t_jsb_var_real3d), POINTER :: carbon14
    TYPE(t_jsb_var_real3d), POINTER :: nitrogen15
    ! component pools
    TYPE(t_sb_elements_pool), POINTER :: dom                 !< elements of dissolved organic matter
    TYPE(t_sb_elements_pool), POINTER :: dom_assoc           !< elements of minerally associated disolved organic matter
    TYPE(t_sb_elements_pool), POINTER :: soluable_litter     !< elements of metabolic litter (sugars, cellulose etc.)
    TYPE(t_sb_elements_pool), POINTER :: polymeric_litter    !< elements of structural litter (lignified litter)
    TYPE(t_sb_elements_pool), POINTER :: woody_litter        !< elements of physically protected wood litter
    TYPE(t_sb_elements_pool), POINTER :: fungi               !< elements of saprophytic fungal biomass
    TYPE(t_sb_elements_pool), POINTER :: mycorrhiza          !< elements of mycorrhizal fungal biomass 
    TYPE(t_sb_elements_pool), POINTER :: microbial           !< elements of microbial biomass (bacteria!)
    TYPE(t_sb_elements_pool), POINTER :: residue             !< elements of necromass of fungi, microbial and mycorrhiza
    TYPE(t_sb_elements_pool), POINTER :: residue_assoc       !< elements of minerally associated microbial necromass 
  CONTAINS
    !! infrastructure
    PROCEDURE :: Init                     => t_sb_components_pool_init
    PROCEDURE :: Get_element_name_by_id   => t_sb_components_pool_get_element_name_by_id
    PROCEDURE :: Init_comp_pool_for_op_ov => Create_sb_components_pool_temporary_object_for_op_ov    ! create 
    !! operator overload
    !! [def: field2d(nc,vgrid) | field3d(nc,vgrid,nblks)]
    ! +
    PROCEDURE                 :: Add_sb_comp_p_and_scalar            ! elements_pool + scalar        | ON_CHUNK & ON_DOMAIN
    ! PROCEDURE, PASS(struct1)  :: Add_scalar_and_sb_comp_p            ! scalar        + elements_pool | ON_CHUNK & ON_DOMAIN
    ! PROCEDURE                 :: Add_sb_comp_p_and_field2d           ! elements_pool + field2d       | ON_CHUNK
    ! PROCEDURE, PASS(struct1)  :: Add_field2d_and_sb_comp_p           ! field2d       + elements_pool | ON_CHUNK
    ! PROCEDURE                 :: Add_sb_comp_p_and_field3d           ! elements_pool + field3d       | ON_CHUNK & ON_DOMAIN
    ! PROCEDURE                 :: Add_sb_comp_p_and_sb_comp_p         ! elements_pool + elements_pool | ON_CHUNK & ON_DOMAIN
    GENERIC                   :: OPERATOR(+) => Add_sb_comp_p_and_scalar !,  Add_scalar_and_sb_comp_p, &
    !   &                                         Add_sb_comp_p_and_field2d, Add_field2d_and_sb_comp_p, &
    !   &                                         Add_sb_comp_p_and_field3d, &
    !   &                                         Add_sb_comp_p_and_sb_comp_p
    ! ! -
    ! PROCEDURE                 :: Subtract_sb_comp_p_minus_scalar     ! elements_pool - scalar        | ON_CHUNK & ON_DOMAIN
    ! PROCEDURE                 :: Subtract_sb_comp_p_minus_field2d    ! elements_pool - field2d       | ON_CHUNK
    ! PROCEDURE                 :: Subtract_sb_comp_p_minus_field3d    ! elements_pool - field3d       | ON_CHUNK & ON_DOMAIN
    ! PROCEDURE                 :: Subtract_sb_comp_p_minus_sb_comp_p  ! elements_pool - elements_pool | ON_CHUNK & ON_DOMAIN
    ! GENERIC                   :: OPERATOR(-) => Subtract_sb_comp_p_minus_scalar, &
    !   &                                         Subtract_sb_comp_p_minus_field2d, &
    !   &                                         Subtract_sb_comp_p_minus_field3d, &
    !   &                                         Subtract_sb_comp_p_minus_sb_comp_p
    ! ! *
    ! PROCEDURE                 :: Multiply_sb_comp_p_with_scalar      ! elements_pool * scalar        | ON_CHUNK & ON_DOMAIN
    ! PROCEDURE, PASS(struct1)  :: Multiply_scalar_with_sb_comp_p      ! scalar        * elements_pool | ON_CHUNK & ON_DOMAIN
    ! PROCEDURE                 :: Multiply_sb_comp_p_with_field2d     ! elements_pool * field2d       | ON_CHUNK
    ! PROCEDURE, PASS(struct1)  :: Multiply_field2d_with_sb_comp_p     ! field2d       * elements_pool | ON_CHUNK
    ! PROCEDURE                 :: Multiply_sb_comp_p_with_field3d     ! elements_pool * field3d       | ON_CHUNK & ON_DOMAIN
    ! PROCEDURE, PASS(struct1)  :: Multiply_field3d_with_sb_comp_p     ! field3d       * elements_pool | ON_CHUNK & ON_DOMAIN
    ! PROCEDURE                 :: Multiply_sb_comp_p_with_sb_comp_p   ! elements_pool * elements_pool | ON_CHUNK & ON_DOMAIN
    ! GENERIC                   :: OPERATOR(*) => Multiply_sb_comp_p_with_scalar,  Multiply_scalar_with_sb_comp_p, &
    !   &                                         Multiply_sb_comp_p_with_field2d, Multiply_field2d_with_sb_comp_p, &
    !   &                                         Multiply_sb_comp_p_with_field3d, Multiply_field3d_with_sb_comp_p, &
    !   &                                         Multiply_sb_comp_p_with_sb_comp_p
    ! ! /
    ! PROCEDURE                 :: Div_sb_comp_p_by_scalar             ! elements_pool / scalar        | ON_CHUNK & ON_DOMAIN
    ! PROCEDURE                 :: Div_sb_comp_p_by_field2d            ! elements_pool / field2d       | ON_CHUNK
    ! PROCEDURE                 :: Div_sb_comp_p_by_field3d            ! elements_pool / field3d       | ON_CHUNK & ON_DOMAIN
    ! PROCEDURE                 :: Div_sb_comp_p_by_sb_comp_p          ! elements_pool / elements_pool | ON_CHUNK & ON_DOMAIN
    ! GENERIC                   :: OPERATOR(/) => Div_sb_comp_p_by_scalar, &
    !   &                                         Div_sb_comp_p_by_field2d, &
    !   &                                         Div_sb_comp_p_by_field3d, &
    !   &                                         Div_sb_comp_p_by_sb_comp_p
    ! =
    PROCEDURE                 :: Assign_scalar_to_sb_comp_p          ! elements_pool = scalar        | ON_CHUNK & ON_DOMAIN
    ! PROCEDURE                 :: Assign_field2d_to_sb_comp_p         ! elements_pool = field2d(:,:)  | ON_CHUNK
    ! PROCEDURE                 :: Assign_field3d_to_sb_comp_p         ! elements_pool = field3d(:,:,:)| ON_CHUNK & ON_DOMAIN
    PROCEDURE                 :: Assign_sb_comp_p_to_sb_comp_p       ! elements_pool = elements_pool | ON_CHUNK & ON_DOMAIN | only numeric values of elements are assigned
    GENERIC                   :: ASSIGNMENT(=) => Assign_scalar_to_sb_comp_p , &
    !   &                                           Assign_field2d_to_sb_comp_p, &
    !   &                                           Assign_field3d_to_sb_comp_p, &
      &                                           Assign_sb_comp_p_to_sb_comp_p
    !! finalize
    FINAL     :: finalize_sb_components_pool
  END TYPE t_sb_components_pool

  !-----------------------------------------------------------------------------------------------------
  !> Type definition for t_sb_main_pool
  !! pool-structure containing pools and fluxes
  !! 1st pool in line
  !! added directly into the t_sb_memory as "CLASS(t_jsb_pool), POINTER :: pools"
  !-----------------------------------------------------------------------------------------------------
  TYPE, EXTENDS(t_sb_elements_pool) :: t_sb_main_pool
    !> 1.0 element variables
    !! inherited from t_sb_elements_pool
    ! TYPE(t_jsb_var_real3d), POINTER :: carbon, nitrogen, phosphorus, carbon13, carbon14, nitrogen15
    !> 2.0 components_pools
    !! pools
    TYPE(t_sb_components_pool),   POINTER :: sb_pool               !< 
    !! fluxes
    TYPE(t_sb_components_pool),   POINTER :: sb_formation          !< 
    TYPE(t_sb_components_pool),   POINTER :: sb_loss               !< 
    TYPE(t_sb_components_pool),   POINTER :: sb_transport          !< 
    !> 3.0 element_pools
    !! pools
    TYPE(t_sb_elements_pool),     POINTER :: sb_delta_pools        !< 
    !! fluxes
    TYPE(t_sb_elements_pool),     POINTER :: sb_mycorrhiza_export  !< 
  END TYPE t_sb_main_pool

  !-----------------------------------------------------------------------------------------------------
  !> list with the element_id of all elements as defined above
  !! is allocated in Init_PROC_memory
  !-----------------------------------------------------------------------------------------------------
  INTEGER, ALLOCATABLE, SAVE :: all_elements(:)

  !-----------------------------------------------------------------------------------------------------
  !> Type definition for sb memory
  !!
  !! 
  !-----------------------------------------------------------------------------------------------------
  TYPE, EXTENDS(t_jsb_memory) :: t_sb_memory

    !> pools & fluxes
    !! 
    ! the sb_main_pool does contain all other pools & fluxes 
    TYPE(t_sb_main_pool), POINTER :: sb_main_pool

    ! 2D total fluxes for mass balance checks in soil_biogeochemistry model
    TYPE(t_jsb_var_real2d)          :: & 
                                    total_flux_carbon, &          !< carbon mass balance checks in soil_biogeochemistry model [mol C m-2 s-1]
                                    total_flux_nitrogen, &        !< nitrogen mass balance checks in soil_biogeochemistry model [mol N m-2 s-1]
                                    total_flux_phosphorus, &      !< phosphorus mass balance checks in soil_biogeochemistry model [mol P m-2 s-1]
                                    total_flux_carbon13, &        !< carbon13 mass balance checks in soil_biogeochemistry model [mol C13 m-2 s-1]
                                    total_flux_carbon14, &        !< carbon14 mass balance checks in soil_biogeochemistry model [mol C14 m-2 s-1]
                                    total_flux_nitrogen15         !< nitrogen15 mass balance checks in soil_biogeochemistry model [mol N15 m-2 s-1]

    ! 3D soil properties
    TYPE(t_jsb_var_real3d)          :: &
                                    ph_sl                         !< ph per soil layer

  CONTAINS
    PROCEDURE :: Init => Init_sb_memory
  END TYPE t_sb_memory

  CHARACTER(len=*), PARAMETER :: modname = 'mo_sb_memory_class'

CONTAINS

  !-----------------------------------------------------------------------------------------------------
  !> helper function for t_sb_elements_pool
  !! PROCEDURE :: Get_element_name_by_id => t_sb_elements_pool_get_element_name_by_id
  !-----------------------------------------------------------------------------------------------------
  FUNCTION t_sb_elements_pool_get_element_name_by_id(this, id) RESULT(name)

    CLASS(t_sb_elements_pool),  INTENT(in)  :: this
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

  END FUNCTION t_sb_elements_pool_get_element_name_by_id

  !-----------------------------------------------------------------------------------------------------
  !> initializing an object of the SB_ elements pool "t_sb_elements_pool"
  !! "pool_object%Init()"
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE t_sb_elements_pool_init(this, name, shortname, l_elements, element_list, element_unit)

    CLASS(t_sb_elements_pool),  INTENT(inout) :: this               ! object of TYPE t_sb_elements_pool
    CHARACTER(LEN=*),           INTENT(in)    :: name               ! pool name       (No spaces !)
    CHARACTER(LEN=*),           INTENT(in)    :: shortname          ! pool short name (No spaces !)
    LOGICAL,          OPTIONAL, INTENT(in)    :: l_elements         ! IF pool has elements 
    INTEGER,          OPTIONAL, INTENT(in)    :: element_list(:)    ! list of elements; use if only selected elements should be added
    CHARACTER(LEN=*), OPTIONAL, INTENT(in)    :: element_unit       ! unit of elements

    LOGICAL                     :: l_add_elements
    INTEGER, ALLOCATABLE        :: element_list_loc(:)
    INTEGER                     :: i
    CHARACTER(len=*), PARAMETER :: routine = modname//':t_sb_elements_pool_init'

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
      CALL this%Add_element('carbon', 'C', ELEM_C, dim=3)
      i = this%Find_element_position(ELEM_C)
      IF (i > 0) THEN
        SELECT TYPE (selector => this%element_list(i)%p)
        CLASS IS (t_jsb_var_real3d)
          this%carbon => selector
        END SELECT
      END IF
    END IF
    !> nitrogen
    IF (One_of(ELEM_N, element_list_loc) > 0) THEN
      CALL this%Add_element('nitrogen', 'N', ELEM_N, dim=3)
      i = this%Find_element_position(ELEM_N)
      IF (i > 0) THEN
        SELECT TYPE (selector => this%element_list(i)%p)
        CLASS IS (t_jsb_var_real3d)
          this%nitrogen => selector
        END SELECT
      END IF
    END IF
    !> phosphorus
    IF (One_of(ELEM_P, element_list_loc) > 0) THEN
      CALL this%Add_element('phosphorus', 'P', ELEM_P, dim=3)
      i = this%Find_element_position(ELEM_P)
      IF (i > 0) THEN
        SELECT TYPE (selector => this%element_list(i)%p)
        CLASS IS (t_jsb_var_real3d)
          this%phosphorus => selector
        END SELECT
      END IF
    END IF
    !> carbon13
    IF (One_of(ELEM_C13, element_list_loc) > 0) THEN
      CALL this%Add_element('carbon13', 'C13', ELEM_C13, dim=3)
      i = this%Find_element_position(ELEM_C13)
      IF (i > 0) THEN
        SELECT TYPE (selector => this%element_list(i)%p)
        CLASS IS (t_jsb_var_real3d)
          this%carbon13 => selector
        END SELECT
      END IF
    END IF
    !> carbon14
    IF (One_of(ELEM_C14, element_list_loc) > 0) THEN
      CALL this%Add_element('carbon14', 'C14', ELEM_C14, dim=3)
      i = this%Find_element_position(ELEM_C14)
      IF (i > 0) THEN
        SELECT TYPE (selector => this%element_list(i)%p)
        CLASS IS (t_jsb_var_real3d)
          this%carbon14 => selector
        END SELECT
      END IF
    END IF
    !> nitrogen15
    IF (One_of(ELEM_N15, element_list_loc) > 0) THEN
      CALL this%Add_element('nitrogen15', 'N15', ELEM_N15, dim=3)
      i = this%Find_element_position(ELEM_N15)
      IF (i > 0) THEN
        SELECT TYPE (selector => this%element_list(i)%p)
        CLASS IS (t_jsb_var_real3d)
          this%nitrogen15 => selector
        END SELECT
      END IF
    END IF

  END SUBROUTINE t_sb_elements_pool_init

  !-----------------------------------------------------------------------------------------------------
  !> helper function for t_sb_components_pool
  !! PROCEDURE :: Get_element_name_by_id => t_sb_components_pool_get_element_name_by_id
  !-----------------------------------------------------------------------------------------------------
  FUNCTION t_sb_components_pool_get_element_name_by_id(this, id) RESULT(name)

    CLASS(t_sb_components_pool),  INTENT(in)  :: this
    INTEGER,                      INTENT(in)  :: id
    CHARACTER(LEN=:),             ALLOCATABLE :: name

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

  END FUNCTION t_sb_components_pool_get_element_name_by_id

  !-----------------------------------------------------------------------------------------------------
  !> initializing the element-variables of an object of the SB_ components pool "t_sb_components_pool"
  !! "pool_object%Init()"
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE t_sb_components_pool_init(this, name, shortname, l_elements, element_list, element_unit)

    CLASS(t_sb_components_pool),  INTENT(inout) :: this               ! object of TYPE t_sb_components_pool
    CHARACTER(LEN=*),             INTENT(in)    :: name               ! pool name       (No spaces !)
    CHARACTER(LEN=*),             INTENT(in)    :: shortname          ! pool short name (No spaces !)
    LOGICAL,          OPTIONAL,   INTENT(in)    :: l_elements         ! IF pool has elements 
    INTEGER,          OPTIONAL,   INTENT(in)    :: element_list(:)    ! list of elements; use if only selected elements should be added
    CHARACTER(LEN=*), OPTIONAL,   INTENT(in)    :: element_unit       ! unit of elements

    LOGICAL                     :: l_add_elements
    INTEGER, ALLOCATABLE        :: element_list_loc(:)
    INTEGER                     :: i
    CHARACTER(len=*), PARAMETER :: routine = modname//':t_sb_components_pool_init'

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
      CALL this%Add_element('carbon', 'C', ELEM_C, dim=3)
      i = this%Find_element_position(ELEM_C)
      IF (i > 0) THEN
        SELECT TYPE (selector => this%element_list(i)%p)
        CLASS IS (t_jsb_var_real3d)
          this%carbon => selector
        END SELECT
      END IF
    END IF
    !> nitrogen
    IF (One_of(ELEM_N, element_list_loc) > 0) THEN
      CALL this%Add_element('nitrogen', 'N', ELEM_N, dim=3)
      i = this%Find_element_position(ELEM_N)
      IF (i > 0) THEN
        SELECT TYPE (selector => this%element_list(i)%p)
        CLASS IS (t_jsb_var_real3d)
          this%nitrogen => selector
        END SELECT
      END IF
    END IF
    !> phosphorus
    IF (One_of(ELEM_P, element_list_loc) > 0) THEN
      CALL this%Add_element('phosphorus', 'P', ELEM_P, dim=3)
      i = this%Find_element_position(ELEM_P)
      IF (i > 0) THEN
        SELECT TYPE (selector => this%element_list(i)%p)
        CLASS IS (t_jsb_var_real3d)
          this%phosphorus => selector
        END SELECT
      END IF
    END IF
    !> carbon13
    IF (One_of(ELEM_C13, element_list_loc) > 0) THEN
      CALL this%Add_element('carbon13', 'C13', ELEM_C13, dim=3)
      i = this%Find_element_position(ELEM_C13)
      IF (i > 0) THEN
        SELECT TYPE (selector => this%element_list(i)%p)
        CLASS IS (t_jsb_var_real3d)
          this%carbon13 => selector
        END SELECT
      END IF
    END IF
    !> carbon14
    IF (One_of(ELEM_C14, element_list_loc) > 0) THEN
      CALL this%Add_element('carbon14', 'C14', ELEM_C14, dim=3)
      i = this%Find_element_position(ELEM_C14)
      IF (i > 0) THEN
        SELECT TYPE (selector => this%element_list(i)%p)
        CLASS IS (t_jsb_var_real3d)
          this%carbon14 => selector
        END SELECT
      END IF
    END IF
    !> nitrogen15
    IF (One_of(ELEM_N15, element_list_loc) > 0) THEN
      CALL this%Add_element('nitrogen15', 'N15', ELEM_N15, dim=3)
      i = this%Find_element_position(ELEM_N15)
      IF (i > 0) THEN
        SELECT TYPE (selector => this%element_list(i)%p)
        CLASS IS (t_jsb_var_real3d)
          this%nitrogen15 => selector
        END SELECT
      END IF
    END IF

  END SUBROUTINE t_sb_components_pool_init

  !-----------------------------------------------------------------------------------------------------
  !> create and init an object of the SB_ tissues pool "t_sb_components_pool"
  !! init and add all its objects of the "t_sb_elements_pool" with "pool_object%Init()"
  !! "object => Create_sb_components_pool()"
  !-----------------------------------------------------------------------------------------------------
  FUNCTION Create_sb_components_pool(name, shortname, unit, lct_ids, pool_prefix, pool_suffix) RESULT(sb_components_pool)

    CHARACTER(LEN=*),           INTENT(in)  :: name
    CHARACTER(LEN=*),           INTENT(in)  :: shortname
    CHARACTER(LEN=*),           INTENT(in)  :: unit
    INTEGER,                    INTENT(in)  :: lct_ids(:)
    CHARACTER(LEN=*), OPTIONAL, INTENT(in)  :: pool_prefix
    CHARACTER(LEN=*), OPTIONAL, INTENT(in)  :: pool_suffix
    TYPE(t_sb_components_pool), POINTER     :: sb_components_pool

    TYPE(t_sb_elements_pool),   POINTER     :: sb_elements_pool
    CHARACTER(LEN=:),           ALLOCATABLE :: prefix_loc, suffix_loc

    prefix_loc = ''
    suffix_loc = ''
    CALL assign_if_present_allocatable(prefix_loc, pool_prefix)
    CALL assign_if_present_allocatable(suffix_loc, pool_suffix)
    IF (prefix_loc /= '') prefix_loc = prefix_loc // '_'
    IF (suffix_loc /= '') suffix_loc = '_' // suffix_loc

    ! tissue pool init
    ALLOCATE(sb_components_pool)
    CALL sb_components_pool%Init(name, shortname, l_elements=.TRUE., element_unit=unit)
    
    ! element pools allocate & init & add to tissue pool
    ! "shortname" is used for model output
    ! "element_list" is optional and can be used to select specific elements (not needed when all elements may be added)
    ! dom
    ALLOCATE(sb_elements_pool)
    CALL sb_elements_pool%Init(name='dom_pool'//suffix_loc, shortname='dom', l_elements=.TRUE., element_unit=unit)
    CALL sb_components_pool%Add_pool(sb_elements_pool)
    sb_components_pool%dom => sb_elements_pool
    NULLIFY(sb_elements_pool)
    ! dom_assoc
    ALLOCATE(sb_elements_pool)
    CALL sb_elements_pool%Init(name='dom_assoc_pool'//suffix_loc, shortname='dom_assoc', l_elements=.TRUE., element_unit=unit)
    CALL sb_components_pool%Add_pool(sb_elements_pool)
    sb_components_pool%dom_assoc => sb_elements_pool
    NULLIFY(sb_elements_pool)
    ! soluable_litter
    ALLOCATE(sb_elements_pool)
    CALL sb_elements_pool%Init(name='soluable_litter_pool'//suffix_loc, shortname='sol_litter', l_elements=.TRUE., &
                               element_unit=unit)
    CALL sb_components_pool%Add_pool(sb_elements_pool)
    sb_components_pool%soluable_litter => sb_elements_pool
    NULLIFY(sb_elements_pool)
    ! polymeric_litter
    ALLOCATE(sb_elements_pool)
    CALL sb_elements_pool%Init(name='polymeric_litter_pool'//suffix_loc, shortname='pol_litter', l_elements=.TRUE., &
                               element_unit=unit)
    CALL sb_components_pool%Add_pool(sb_elements_pool)
    sb_components_pool%polymeric_litter => sb_elements_pool
    NULLIFY(sb_elements_pool)
    ! woody_litter
    ALLOCATE(sb_elements_pool)
    CALL sb_elements_pool%Init(name='woody_litter_pool'//suffix_loc, shortname='woo_litter', l_elements=.TRUE., element_unit=unit)
    CALL sb_components_pool%Add_pool(sb_elements_pool)
    sb_components_pool%woody_litter => sb_elements_pool
    NULLIFY(sb_elements_pool)
    ! fungi
    ALLOCATE(sb_elements_pool)
    CALL sb_elements_pool%Init(name='fungi_pool'//suffix_loc, shortname='fungi', l_elements=.TRUE., element_unit=unit)
    CALL sb_components_pool%Add_pool(sb_elements_pool)
    sb_components_pool%fungi => sb_elements_pool
    NULLIFY(sb_elements_pool)
    ! mycorrhiza
    ALLOCATE(sb_elements_pool)
    CALL sb_elements_pool%Init(name='mycorrhiza_pool'//suffix_loc, shortname='mycorrhiza', l_elements=.TRUE., element_unit=unit)
    CALL sb_components_pool%Add_pool(sb_elements_pool)
    sb_components_pool%mycorrhiza => sb_elements_pool
    NULLIFY(sb_elements_pool)
    ! microbial
    ALLOCATE(sb_elements_pool)
    CALL sb_elements_pool%Init(name='microbial_pool'//suffix_loc, shortname='microbial', l_elements=.TRUE., element_unit=unit)
    CALL sb_components_pool%Add_pool(sb_elements_pool)
    sb_components_pool%microbial => sb_elements_pool
    NULLIFY(sb_elements_pool)
    ! residue
    ALLOCATE(sb_elements_pool)
    CALL sb_elements_pool%Init(name='residue_pool'//suffix_loc, shortname='residue', l_elements=.TRUE., element_unit=unit)
    CALL sb_components_pool%Add_pool(sb_elements_pool)
    sb_components_pool%residue => sb_elements_pool
    NULLIFY(sb_elements_pool)
    ! residue_assoc
    ALLOCATE(sb_elements_pool)
    CALL sb_elements_pool%Init(name='residue_assoc_pool'//suffix_loc, shortname='residue_assoc', l_elements=.TRUE., &
                               element_unit=unit)
    CALL sb_components_pool%Add_pool(sb_elements_pool)
    sb_components_pool%residue_assoc => sb_elements_pool
    NULLIFY(sb_elements_pool)

  END FUNCTION Create_sb_components_pool

  !-----------------------------------------------------------------------------------------------------
  !> create a components pool (t_sb_components_pool)
  !! used for operator overload routines of components pools
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE Create_sb_components_pool_temporary_object_for_op_ov(this, name, shortname, unit, pool_suffix)

    CLASS(t_sb_components_pool),  INTENT(inout) :: this
    CHARACTER(LEN=*),             INTENT(in)    :: name
    CHARACTER(LEN=*),             INTENT(in)    :: shortname
    CHARACTER(LEN=*),             INTENT(in)    :: unit
    CHARACTER(LEN=*), OPTIONAL,   INTENT(in)    :: pool_suffix

    TYPE(t_sb_elements_pool),     POINTER       :: sb_elements_pool
    CHARACTER(LEN=:),             ALLOCATABLE   :: suffix_loc

    ! pool name suffix
    suffix_loc = ''
    CALL assign_if_present_allocatable(suffix_loc, pool_suffix)
    IF (suffix_loc /= '') suffix_loc = '_' // suffix_loc

    ! components pool init element var
    CALL this%Init(name, shortname, l_elements=.TRUE., element_unit=unit)
    
    ! element pools allocate & init & add to components pool
    ! "shortname" is used for model output
    ! "element_list" is optional and can be used to select specific elements (not needed when all elements may be added)
    ! dom
    ALLOCATE(sb_elements_pool)
    CALL sb_elements_pool%Init(name='dom_pool'//suffix_loc, shortname='dom', l_elements=.TRUE., element_unit=unit)
    CALL this%Add_pool(sb_elements_pool)
    this%dom => sb_elements_pool
    NULLIFY(sb_elements_pool)
    ! dom_assoc
    ALLOCATE(sb_elements_pool)
    CALL sb_elements_pool%Init(name='dom_assoc_pool'//suffix_loc, shortname='dom_assoc', l_elements=.TRUE., element_unit=unit)
    CALL this%Add_pool(sb_elements_pool)
    this%dom_assoc => sb_elements_pool
    NULLIFY(sb_elements_pool)
    ! soluable_litter
    ALLOCATE(sb_elements_pool)
    CALL sb_elements_pool%Init(name='soluable_litter_pool'//suffix_loc, shortname='sol_litter', l_elements=.TRUE., &
                               element_unit=unit)
    CALL this%Add_pool(sb_elements_pool)
    this%soluable_litter => sb_elements_pool
    NULLIFY(sb_elements_pool)
    ! polymeric_litter
    ALLOCATE(sb_elements_pool)
    CALL sb_elements_pool%Init(name='polymeric_litter_pool'//suffix_loc, shortname='pol_litter', l_elements=.TRUE., &
                               element_unit=unit)
    CALL this%Add_pool(sb_elements_pool)
    this%polymeric_litter => sb_elements_pool
    NULLIFY(sb_elements_pool)
    ! woody_litter
    ALLOCATE(sb_elements_pool)
    CALL sb_elements_pool%Init(name='woody_litter_pool'//suffix_loc, shortname='woo_litter', l_elements=.TRUE., element_unit=unit)
    CALL this%Add_pool(sb_elements_pool)
    this%woody_litter => sb_elements_pool
    NULLIFY(sb_elements_pool)
    ! fungi
    ALLOCATE(sb_elements_pool)
    CALL sb_elements_pool%Init(name='fungi_pool'//suffix_loc, shortname='fungi', l_elements=.TRUE., element_unit=unit)
    CALL this%Add_pool(sb_elements_pool)
    this%fungi => sb_elements_pool
    NULLIFY(sb_elements_pool)
    ! mycorrhiza
    ALLOCATE(sb_elements_pool)
    CALL sb_elements_pool%Init(name='mycorrhiza_pool'//suffix_loc, shortname='mycorrhiza', l_elements=.TRUE., element_unit=unit)
    CALL this%Add_pool(sb_elements_pool)
    this%mycorrhiza => sb_elements_pool
    NULLIFY(sb_elements_pool)
    ! microbial
    ALLOCATE(sb_elements_pool)
    CALL sb_elements_pool%Init(name='microbial_pool'//suffix_loc, shortname='microbial', l_elements=.TRUE., element_unit=unit)
    CALL this%Add_pool(sb_elements_pool)
    this%microbial => sb_elements_pool
    NULLIFY(sb_elements_pool)
    ! residue
    ALLOCATE(sb_elements_pool)
    CALL sb_elements_pool%Init(name='residue_pool'//suffix_loc, shortname='residue', l_elements=.TRUE., element_unit=unit)
    CALL this%Add_pool(sb_elements_pool)
    this%residue => sb_elements_pool
    NULLIFY(sb_elements_pool)
    ! residue_assoc
    ALLOCATE(sb_elements_pool)
    CALL sb_elements_pool%Init(name='residue_assoc_pool'//suffix_loc, shortname='residue_assoc', l_elements=.TRUE., &
                               element_unit=unit)
    CALL this%Add_pool(sb_elements_pool)
    this%residue_assoc => sb_elements_pool
    NULLIFY(sb_elements_pool)

  END SUBROUTINE Create_sb_components_pool_temporary_object_for_op_ov

  !-----------------------------------------------------------------------------------------------------
  !> initialize an elements_pool that is part of a components pool
  !! used for operator overload routines of components pools
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE Init_element_var_sb_elements_pool(this, struct1)

    CLASS(t_sb_elements_pool),  INTENT(inout)  :: this
    CLASS(t_sb_elements_pool),  INTENT(in)     :: struct1

    this%carbon%owner_model_id      = struct1%carbon%owner_model_id
    this%carbon%owner_proc_id       = struct1%carbon%owner_proc_id
    this%carbon%owner_tile_path     = struct1%carbon%owner_tile_path
    ALLOCATE(this%carbon%ptr3d,     source=struct1%carbon%ptr3d)
    this%nitrogen%owner_model_id    = struct1%nitrogen%owner_model_id
    this%nitrogen%owner_proc_id     = struct1%nitrogen%owner_proc_id
    this%nitrogen%owner_tile_path   = struct1%nitrogen%owner_tile_path
    ALLOCATE(this%nitrogen%ptr3d,   source=struct1%nitrogen%ptr3d)
    this%phosphorus%owner_model_id  = struct1%phosphorus%owner_model_id
    this%phosphorus%owner_proc_id   = struct1%phosphorus%owner_proc_id
    this%phosphorus%owner_tile_path = struct1%phosphorus%owner_tile_path
    ALLOCATE(this%phosphorus%ptr3d, source=struct1%phosphorus%ptr3d)
    this%carbon13%owner_model_id    = struct1%carbon13%owner_model_id
    this%carbon13%owner_proc_id     = struct1%carbon13%owner_proc_id
    this%carbon13%owner_tile_path   = struct1%carbon13%owner_tile_path
    ALLOCATE(this%carbon13%ptr3d,   source=struct1%carbon13%ptr3d)
    this%carbon14%owner_model_id    = struct1%carbon14%owner_model_id
    this%carbon14%owner_proc_id     = struct1%carbon14%owner_proc_id
    this%carbon14%owner_tile_path   = struct1%carbon14%owner_tile_path
    ALLOCATE(this%carbon14%ptr3d,   source=struct1%carbon14%ptr3d)
    this%nitrogen15%owner_model_id  = struct1%nitrogen15%owner_model_id
    this%nitrogen15%owner_proc_id   = struct1%nitrogen15%owner_proc_id
    this%nitrogen15%owner_tile_path = struct1%nitrogen15%owner_tile_path
    ALLOCATE(this%nitrogen15%ptr3d, source=struct1%nitrogen15%ptr3d)

  END SUBROUTINE Init_element_var_sb_elements_pool

  !-----------------------------------------------------------------------------------------------------
  !> initialize memory for the SB_ process
  !! pools and 
  !! 2d/3d variables
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE Init_sb_memory(mem, prefix, suffix, lct_ids, model_id)

    USE mo_jsb_varlist,       ONLY: BASIC, MEDIUM, FULL
    USE mo_jsb_io,            ONLY: grib_bits, t_cf, t_grib1, t_grib2, tables
    USE mo_jsb_grid_class,    ONLY: t_jsb_grid, t_jsb_vgrid
    USE mo_jsb_grid,          ONLY: Get_grid, Get_vgrid
    USE mo_jsb_model_class,   ONLY: t_jsb_model
    USE mo_jsb_class,         ONLY: Get_model

    CLASS(t_sb_memory),   INTENT(inout), TARGET  :: mem             !> SB_ memory
    CHARACTER(len=*),     INTENT(in)             :: prefix          !> process name
    CHARACTER(len=*),     INTENT(in)             :: suffix          !> tile name
    INTEGER,              INTENT(in)             :: lct_ids(:)      !> vector of land cover type IDs
    INTEGER,              INTENT(in)             :: model_id        !> model ID model\%id
    ! ---------------------------
    TYPE(t_jsb_model),          POINTER :: model
    TYPE(t_jsb_grid),           POINTER :: hgrid                    ! Horizontal grid
    TYPE(t_jsb_vgrid),          POINTER :: vgrid_surface            ! Vertical grid
    TYPE(t_jsb_vgrid),          POINTER :: vgrid_soil_sb            ! Vertical grid
    TYPE(t_sb_elements_pool),   POINTER :: sb_elements_pool         ! pools/fluxes like sb_delta_pools & sb_mycorrhiza_export
    TYPE(t_sb_components_pool), POINTER :: sb_components_pool       ! pools/fluxes like sb_pool & sb_formation
    INTEGER                     :: table
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':Init_sb_memory'
    ! ---------------------------
    model         => Get_model(model_id)
    table         =  tables(1)
    hgrid         => Get_grid(model%grid_id)
    vgrid_surface => Get_vgrid('surface')
    vgrid_soil_sb => Get_vgrid('soil_layer_sb')


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
      ALLOCATE(t_sb_main_pool :: mem%pools)
      CALL mem%pools%Init(name='soil_biogeochemistry_main_pool', shortname='sb_main_pool', l_elements=.FALSE.)

      SELECT TYPE (selected => mem%pools)
      CLASS IS (t_sb_main_pool)
        mem%sb_main_pool => selected
      END SELECT

      !-----------------------------------------------------------------------------------------------------------------------------
      !> create the vegetation tissue pools & fluxes (sb_components_pool) 
      !! sb_pool / sb_formation / ...
      ! create pools & fluxes
      ! "shortname" is used for model output names
      ! sb_pool
      sb_components_pool => Create_sb_components_pool(name='soil_biogeochemistry_pool', shortname='sb_pool', unit='mol m-3', &
        & lct_ids=lct_ids, pool_suffix='sbpool') ! l_elements=.TRUE. by default
      CALL mem%sb_main_pool%Add_pool(sb_components_pool)
      mem%sb_main_pool%sb_pool => sb_components_pool
      ! sb_formation
      sb_components_pool => Create_sb_components_pool(name='sb_formation_flux', shortname='sb_formation', unit='mol m-3', &
        & lct_ids=lct_ids, pool_suffix='sbform') ! l_elements=.TRUE. by default
      CALL mem%sb_main_pool%Add_pool(sb_components_pool)
      mem%sb_main_pool%sb_formation => sb_components_pool
      ! sb_loss
      sb_components_pool => Create_sb_components_pool(name='sb_loss_flux', shortname='sb_loss', unit='mol m-3', &
        & lct_ids=lct_ids, pool_suffix='sbloss') ! l_elements=.TRUE. by default
      CALL mem%sb_main_pool%Add_pool(sb_components_pool)
      mem%sb_main_pool%sb_loss => sb_components_pool
      ! sb_transport
      sb_components_pool => Create_sb_components_pool(name='sb_transport_flux', shortname='sb_transport', unit='mol m-3', &
        & lct_ids=lct_ids, pool_suffix='sbtran') ! l_elements=.TRUE. by default
      CALL mem%sb_main_pool%Add_pool(sb_components_pool)
      mem%sb_main_pool%sb_transport => sb_components_pool
      NULLIFY(sb_components_pool)

      ! init and add element pools to the veg_main_pool
      ! (i.e., add them at the same level as the above veg_tissues_pool veg_pool / veg_growth / ...)
      ! "shortname" is used for model output names
      ! sb_delta_pools
      ALLOCATE(sb_elements_pool)
      CALL sb_elements_pool%Init(name='sb_delta_pools_elem_pool', shortname='sb_delta_pools', &
        &                        l_elements=.TRUE., element_unit='mol m-3')
      CALL mem%sb_main_pool%Add_pool(sb_elements_pool)
      mem%sb_main_pool%sb_delta_pools => sb_elements_pool
      NULLIFY(sb_elements_pool)
      ! sb_mycorrhiza_export
      ALLOCATE(sb_elements_pool)
      CALL sb_elements_pool%Init(name='sb_mycorrhiza_export_elem_flux', shortname='sb_mycorrhiza_export', &
        &                        l_elements=.TRUE., element_unit='mol m-3')
      CALL mem%sb_main_pool%Add_pool(sb_elements_pool)
      mem%sb_main_pool%sb_mycorrhiza_export => sb_elements_pool
      NULLIFY(sb_elements_pool)

      ! set paths of all above initialized & added pools (type-bound procedure of t_jsb_pool)
      CALL mem%pools%Set_paths()

      ! add all variables in the sb_components_pool to the SB_ memory
      CALL mem%Add_var(mem%sb_main_pool%sb_pool, &
        & hgrid, vgrid_surface, vgrid_soil_sb,                                                         &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                           &
        & prefix, suffix,                                                                              &
        & loutput=.TRUE.,                                                                              &
        & lrestart=.TRUE., initval_r=0.0_wp )
      CALL mem%Add_var(mem%sb_main_pool%sb_formation, &
        & hgrid, vgrid_surface, vgrid_soil_sb,                                                         &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                           &
        & prefix, suffix,                                                                              &
        & loutput=.TRUE.,                                                                              &
        & lrestart=.FALSE., initval_r=0.0_wp )                                  ! no restart files needed for fluxes
      CALL mem%Add_var(mem%sb_main_pool%sb_loss, &
        & hgrid, vgrid_surface, vgrid_soil_sb,                                                         &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                           &
        & prefix, suffix,                                                                              &
        & loutput=.TRUE.,                                                                              &
        & lrestart=.FALSE., initval_r=0.0_wp )                                  ! no restart files needed for fluxes
      CALL mem%Add_var(mem%sb_main_pool%sb_transport, &
        & hgrid, vgrid_surface, vgrid_soil_sb,                                                         &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                           &
        & prefix, suffix,                                                                              &
        & loutput=.TRUE.,                                                                              &
        & lrestart=.FALSE., initval_r=0.0_wp )                                  ! no restart files needed for fluxes

      ! add all variables in the sb_elements_pool to the SB_ memory
      CALL mem%Add_var(mem%sb_main_pool%sb_delta_pools, &
        & hgrid, vgrid_surface, vgrid_soil_sb,                                                         &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                           &
        & prefix, suffix,                                                                              &
        & loutput=.TRUE.,                                                                              &
        & lrestart=.TRUE., initval_r=0.0_wp )
      CALL mem%Add_var(mem%sb_main_pool%sb_mycorrhiza_export, &
        & hgrid, vgrid_surface, vgrid_soil_sb,                                                         &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                           &
        & prefix, suffix,                                                                              &
        & loutput=.TRUE.,                                                                              &
        & lrestart=.FALSE., initval_r=0.0_wp )                                  ! no restart files needed for fluxes

      ! output to logfile all pools & their variables at model-init 
      ! (maybe use for debugging; its quite a long list)
      !IF (.NOT. ASSOCIATED(mem%parent)) THEN
      !  CALL mem%sb_main_pool%Print()
      !  ! CALL mem%sb_main_pool%sb_pool%Print()
      !END IF

      !=============================================================================================================================
      !> 2D & 3D variables
      !=============================================================================================================================

      CALL mem%Add_var('total_flux_carbon', mem%total_flux_carbon, &
            hgrid, vgrid_surface, &
            t_cf('total_flux_carbon', 'mol C m-2 s-1', 'carbon mass balance checks in soil_biogeochemistry model'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            loutput=.TRUE., &
            initval_r=0.0_wp) ! initvalue

      CALL mem%Add_var('total_flux_nitrogen', mem%total_flux_nitrogen, &
            hgrid, vgrid_surface, &
            t_cf('total_flux_nitrogen', 'mol N m-2 s-1', 'nitrogen mass balance checks in soil_biogeochemistry model'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            loutput=.TRUE., &
            initval_r=0.0_wp) ! initvalue

      CALL mem%Add_var('total_flux_phosphorus', mem%total_flux_phosphorus, &
            hgrid, vgrid_surface, &
            t_cf('total_flux_phosphorus', 'mol P m-2 s-1', 'phosphorus mass balance checks in soil_biogeochemistry model'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            loutput=.TRUE., &
            initval_r=0.0_wp) ! initvalue

      CALL mem%Add_var('total_flux_carbon13', mem%total_flux_carbon13, &
            hgrid, vgrid_surface, &
            t_cf('total_flux_carbon13', 'mol C13 m-2 s-1', 'carbon13 mass balance checks in soil_biogeochemistry model'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            loutput=.TRUE., &
            initval_r=0.0_wp) ! initvalue

      CALL mem%Add_var('total_flux_carbon14', mem%total_flux_carbon14, &
            hgrid, vgrid_surface, &
            t_cf('total_flux_carbon14', 'mol C14 m-2 s-1', 'carbon14 mass balance checks in soil_biogeochemistry model'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            loutput=.TRUE., &
            initval_r=0.0_wp) ! initvalue

      CALL mem%Add_var('total_flux_nitrogen15', mem%total_flux_nitrogen15, &
            hgrid, vgrid_surface, &
            t_cf('total_flux_nitrogen15', 'mol N15 m-2 s-1', 'nitrogen15 mass balance checks in soil_biogeochemistry model'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            loutput=.TRUE., &
            initval_r=0.0_wp) ! initvalue
            
      CALL mem%Add_var('ph_sl', mem%ph_sl, &
            hgrid, vgrid_soil_sb, &
            t_cf('ph_sl', '', 'ph per soil layer'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            loutput=.TRUE., &
            initval_r=0.0_wp) ! initvalue
    ENDIF

  END SUBROUTINE Init_sb_memory

  !--------------------------------------------------------------------------------------------------------------------------------
  !> operator overload and assignment routines for element pools
  !!  TYPE t_sb_elements_pool
  !! 
  !--------------------------------------------------------------------------------------------------------------------------------

  !--------------------------------------------------------------------------------------------------------------------------------
  ! finalize after operator overload
  !--------------------------------------------------------------------------------------------------------------------------------
  SUBROUTINE finalize_sb_elements_pool(this)

    TYPE(t_sb_elements_pool) :: this

    CHARACTER(len=*), PARAMETER :: routine = modname//':finalize_sb_elements_pool'

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
   
  END SUBROUTINE finalize_sb_elements_pool

  !--------------------------------------------------------------------------------------------------------------------------------
  ! = assign
  !--------------------------------------------------------------------------------------------------------------------------------
  !--------------------------------------------------------------------------------------------------------------------------------
  ! t_sb_elements_pool = scalar | ON_CHUNK & ON_DOMAIN
  !--------------------------------------------------------------------------------------------------------------------------------
  SUBROUTINE Assign_scalar_to_sb_elem_p(lhs, rhs)

    CLASS(t_sb_elements_pool),  INTENT(inout) :: lhs
    REAL(wp),                   INTENT(in)    :: rhs

    CHARACTER(len=*), PARAMETER :: routine = modname//':Assign_scalar_to_sb_elem_p'

    lhs%carbon     = rhs
    lhs%nitrogen   = rhs
    lhs%phosphorus = rhs
    lhs%carbon13   = rhs
    lhs%carbon14   = rhs
    lhs%nitrogen15 = rhs

  END SUBROUTINE Assign_scalar_to_sb_elem_p
  !--------------------------------------------------------------------------------------------------------------------------------
  ! t_sb_elements_pool = field2d(:,:) | ON_CHUNK
  !--------------------------------------------------------------------------------------------------------------------------------
  SUBROUTINE Assign_field2d_to_sb_elem_p(lhs, rhs)
  
    CLASS(t_sb_elements_pool),  INTENT(inout) :: lhs
    REAL(wp),                   INTENT(in)    :: rhs(:,:)

    CHARACTER(len=*), PARAMETER :: routine = modname//':Assign_field2d_to_sb_elem_p'

    lhs%carbon     = rhs(:,:)
    lhs%nitrogen   = rhs(:,:)
    lhs%phosphorus = rhs(:,:)
    lhs%carbon13   = rhs(:,:)
    lhs%carbon14   = rhs(:,:)
    lhs%nitrogen15 = rhs(:,:)
  
  END SUBROUTINE Assign_field2d_to_sb_elem_p
  !--------------------------------------------------------------------------------------------------------------------------------
  ! t_sb_elements_pool = field3d(:,:,:) | ON_CHUNK & ON_DOMAIN | field3d(:,:,1) & field3d(:,:,:)
  !--------------------------------------------------------------------------------------------------------------------------------
  SUBROUTINE Assign_field3d_to_sb_elem_p(lhs, rhs)
  
    CLASS(t_sb_elements_pool),  INTENT(inout) :: lhs
    REAL(wp),                   INTENT(in)    :: rhs(:,:,:)

    CHARACTER(len=*), PARAMETER :: routine = modname//':Assign_field3d_to_sb_elem_p'

    lhs%carbon     = rhs(:,:,:)
    lhs%nitrogen   = rhs(:,:,:)
    lhs%phosphorus = rhs(:,:,:)
    lhs%carbon13   = rhs(:,:,:)
    lhs%carbon14   = rhs(:,:,:)
    lhs%nitrogen15 = rhs(:,:,:)
  
  END SUBROUTINE Assign_field3d_to_sb_elem_p
  !--------------------------------------------------------------------------------------------------------------------------------
  ! t_sb_elements_pool = t_sb_elements_pool | ON_CHUNK & ON_DOMAIN
  !--------------------------------------------------------------------------------------------------------------------------------
  SUBROUTINE Assign_sb_elem_p_to_sb_elem_p(lhs, rhs)

    CLASS(t_sb_elements_pool),  INTENT(inout) :: lhs
    CLASS(t_sb_elements_pool),  INTENT(in)    :: rhs

    CHARACTER(len=*), PARAMETER :: routine = modname//':Assign_sb_elem_p_to_sb_elem_p'

    lhs%carbon     = rhs%carbon    
    lhs%nitrogen   = rhs%nitrogen  
    lhs%phosphorus = rhs%phosphorus
    lhs%carbon13   = rhs%carbon13  
    lhs%carbon14   = rhs%carbon14  
    lhs%nitrogen15 = rhs%nitrogen15

  END SUBROUTINE Assign_sb_elem_p_to_sb_elem_p

  !--------------------------------------------------------------------------------------------------------------------------------
  ! + addition
  !--------------------------------------------------------------------------------------------------------------------------------
  !--------------------------------------------------------------------------------------------------------------------------------
  ! t_sb_elements_pool + scalar        | ON_CHUNK & ON_DOMAIN
  !--------------------------------------------------------------------------------------------------------------------------------
  FUNCTION Add_sb_elem_p_and_scalar(struct1, scalar) RESULT(structOut)
  
    CLASS(t_sb_elements_pool),      INTENT(in)  :: struct1
    REAL(wp),                       INTENT(in)  :: scalar
    TYPE(t_sb_elements_pool),       TARGET      :: structOut
  
    CHARACTER(len=*), PARAMETER :: routine = modname//':Add_sb_elem_p_and_scalar'
  
    ! init temporary pool structOut
    CALL structOut%Init(name='opov_help', shortname='opov_help', l_elements=.TRUE., element_unit=TRIM(struct1%element_unit))
    structOut%carbon%owner_model_id      = struct1%carbon%owner_model_id
    structOut%carbon%owner_proc_id       = struct1%carbon%owner_proc_id
    structOut%carbon%owner_tile_path     = struct1%carbon%owner_tile_path
    ALLOCATE(structOut%carbon%ptr3d,     source=struct1%carbon%ptr3d)
    structOut%nitrogen%owner_model_id    = struct1%nitrogen%owner_model_id
    structOut%nitrogen%owner_proc_id     = struct1%nitrogen%owner_proc_id
    structOut%nitrogen%owner_tile_path   = struct1%nitrogen%owner_tile_path
    ALLOCATE(structOut%nitrogen%ptr3d,   source=struct1%nitrogen%ptr3d)
    structOut%phosphorus%owner_model_id  = struct1%phosphorus%owner_model_id
    structOut%phosphorus%owner_proc_id   = struct1%phosphorus%owner_proc_id
    structOut%phosphorus%owner_tile_path = struct1%phosphorus%owner_tile_path
    ALLOCATE(structOut%phosphorus%ptr3d, source=struct1%phosphorus%ptr3d)
    structOut%carbon13%owner_model_id    = struct1%carbon13%owner_model_id
    structOut%carbon13%owner_proc_id     = struct1%carbon13%owner_proc_id
    structOut%carbon13%owner_tile_path   = struct1%carbon13%owner_tile_path
    ALLOCATE(structOut%carbon13%ptr3d,   source=struct1%carbon13%ptr3d)
    structOut%carbon14%owner_model_id    = struct1%carbon14%owner_model_id
    structOut%carbon14%owner_proc_id     = struct1%carbon14%owner_proc_id
    structOut%carbon14%owner_tile_path   = struct1%carbon14%owner_tile_path
    ALLOCATE(structOut%carbon14%ptr3d,   source=struct1%carbon14%ptr3d)
    structOut%nitrogen15%owner_model_id  = struct1%nitrogen15%owner_model_id
    structOut%nitrogen15%owner_proc_id   = struct1%nitrogen15%owner_proc_id
    structOut%nitrogen15%owner_tile_path = struct1%nitrogen15%owner_tile_path
    ALLOCATE(structOut%nitrogen15%ptr3d, source=struct1%nitrogen15%ptr3d)

    ! run calculation
    ! the differentiation between ON_CHUNK/ON_DOMAIN is done at the level of variables (i.e., the op-ov routines of the jsb_var)
    structOut%carbon          = struct1%carbon      + scalar
    structOut%nitrogen        = struct1%nitrogen    + scalar
    structOut%phosphorus      = struct1%phosphorus  + scalar
    structOut%carbon13        = struct1%carbon13    + scalar
    structOut%carbon14        = struct1%carbon14    + scalar
    structOut%nitrogen15      = struct1%nitrogen15  + scalar

    ! is this needed here too (not only with the jsb_var op-ov)
    structOut%carbon%ptr      => structOut%carbon%ptr3d(:,:,:)
    structOut%nitrogen%ptr    => structOut%nitrogen%ptr3d(:,:,:)
    structOut%phosphorus%ptr  => structOut%phosphorus%ptr3d(:,:,:)
    structOut%carbon13%ptr    => structOut%carbon13%ptr3d(:,:,:)
    structOut%carbon14%ptr    => structOut%carbon14%ptr3d(:,:,:)
    structOut%nitrogen15%ptr  => structOut%nitrogen15%ptr3d(:,:,:)
  
  END FUNCTION Add_sb_elem_p_and_scalar
  !--------------------------------------------------------------------------------------------------------------------------------
  ! scalar + t_sb_elements_pool        | ON_CHUNK & ON_DOMAIN
  !--------------------------------------------------------------------------------------------------------------------------------
  FUNCTION Add_scalar_and_sb_elem_p(scalar, struct1) RESULT(structOut)
  
    REAL(wp),                       INTENT(in)  :: scalar
    CLASS(t_sb_elements_pool),      INTENT(in)  :: struct1
    TYPE(t_sb_elements_pool),       TARGET      :: structOut
  
    CHARACTER(len=*), PARAMETER :: routine = modname//':Add_scalar_and_sb_elem_p'
  
    ! init temporary pool structOut
    CALL structOut%Init(name='opov_help', shortname='opov_help', l_elements=.TRUE., element_unit=TRIM(struct1%element_unit))
    structOut%carbon%owner_model_id      = struct1%carbon%owner_model_id
    structOut%carbon%owner_proc_id       = struct1%carbon%owner_proc_id
    structOut%carbon%owner_tile_path     = struct1%carbon%owner_tile_path
    ALLOCATE(structOut%carbon%ptr3d,     source=struct1%carbon%ptr3d)
    structOut%nitrogen%owner_model_id    = struct1%nitrogen%owner_model_id
    structOut%nitrogen%owner_proc_id     = struct1%nitrogen%owner_proc_id
    structOut%nitrogen%owner_tile_path   = struct1%nitrogen%owner_tile_path
    ALLOCATE(structOut%nitrogen%ptr3d,   source=struct1%nitrogen%ptr3d)
    structOut%phosphorus%owner_model_id  = struct1%phosphorus%owner_model_id
    structOut%phosphorus%owner_proc_id   = struct1%phosphorus%owner_proc_id
    structOut%phosphorus%owner_tile_path = struct1%phosphorus%owner_tile_path
    ALLOCATE(structOut%phosphorus%ptr3d, source=struct1%phosphorus%ptr3d)
    structOut%carbon13%owner_model_id    = struct1%carbon13%owner_model_id
    structOut%carbon13%owner_proc_id     = struct1%carbon13%owner_proc_id
    structOut%carbon13%owner_tile_path   = struct1%carbon13%owner_tile_path
    ALLOCATE(structOut%carbon13%ptr3d,   source=struct1%carbon13%ptr3d)
    structOut%carbon14%owner_model_id    = struct1%carbon14%owner_model_id
    structOut%carbon14%owner_proc_id     = struct1%carbon14%owner_proc_id
    structOut%carbon14%owner_tile_path   = struct1%carbon14%owner_tile_path
    ALLOCATE(structOut%carbon14%ptr3d,   source=struct1%carbon14%ptr3d)
    structOut%nitrogen15%owner_model_id  = struct1%nitrogen15%owner_model_id
    structOut%nitrogen15%owner_proc_id   = struct1%nitrogen15%owner_proc_id
    structOut%nitrogen15%owner_tile_path = struct1%nitrogen15%owner_tile_path
    ALLOCATE(structOut%nitrogen15%ptr3d, source=struct1%nitrogen15%ptr3d)

    ! run calculation
    ! the differentiation between ON_CHUNK/ON_DOMAIN is done at the level of variables (i.e., the op-ov routines of the jsb_var)
    structOut%carbon          = scalar + struct1%carbon    
    structOut%nitrogen        = scalar + struct1%nitrogen  
    structOut%phosphorus      = scalar + struct1%phosphorus
    structOut%carbon13        = scalar + struct1%carbon13  
    structOut%carbon14        = scalar + struct1%carbon14  
    structOut%nitrogen15      = scalar + struct1%nitrogen15

    ! is this needed here too (not only with the jsb_var op-ov)
    structOut%carbon%ptr      => structOut%carbon%ptr3d(:,:,:)
    structOut%nitrogen%ptr    => structOut%nitrogen%ptr3d(:,:,:)
    structOut%phosphorus%ptr  => structOut%phosphorus%ptr3d(:,:,:)
    structOut%carbon13%ptr    => structOut%carbon13%ptr3d(:,:,:)
    structOut%carbon14%ptr    => structOut%carbon14%ptr3d(:,:,:)
    structOut%nitrogen15%ptr  => structOut%nitrogen15%ptr3d(:,:,:)
  
  END FUNCTION Add_scalar_and_sb_elem_p
  !--------------------------------------------------------------------------------------------------------------------------------
  ! t_sb_elements_pool + field2d(:,:)        | ON_CHUNK
  !--------------------------------------------------------------------------------------------------------------------------------
  FUNCTION Add_sb_elem_p_and_field2d(struct1, field2d) RESULT(structOut)
  
    CLASS(t_sb_elements_pool),      INTENT(in)  :: struct1
    REAL(wp),                       INTENT(in)  :: field2d(:,:)
    TYPE(t_sb_elements_pool),       TARGET      :: structOut
  
    CHARACTER(len=*), PARAMETER :: routine = modname//':Add_sb_elem_p_and_field2d'
  
    ! init temporary pool structOut
    CALL structOut%Init(name='opov_help', shortname='opov_help', l_elements=.TRUE., element_unit=TRIM(struct1%element_unit))
    structOut%carbon%owner_model_id      = struct1%carbon%owner_model_id
    structOut%carbon%owner_proc_id       = struct1%carbon%owner_proc_id
    structOut%carbon%owner_tile_path     = struct1%carbon%owner_tile_path
    ALLOCATE(structOut%carbon%ptr3d,     source=struct1%carbon%ptr3d)
    structOut%nitrogen%owner_model_id    = struct1%nitrogen%owner_model_id
    structOut%nitrogen%owner_proc_id     = struct1%nitrogen%owner_proc_id
    structOut%nitrogen%owner_tile_path   = struct1%nitrogen%owner_tile_path
    ALLOCATE(structOut%nitrogen%ptr3d,   source=struct1%nitrogen%ptr3d)
    structOut%phosphorus%owner_model_id  = struct1%phosphorus%owner_model_id
    structOut%phosphorus%owner_proc_id   = struct1%phosphorus%owner_proc_id
    structOut%phosphorus%owner_tile_path = struct1%phosphorus%owner_tile_path
    ALLOCATE(structOut%phosphorus%ptr3d, source=struct1%phosphorus%ptr3d)
    structOut%carbon13%owner_model_id    = struct1%carbon13%owner_model_id
    structOut%carbon13%owner_proc_id     = struct1%carbon13%owner_proc_id
    structOut%carbon13%owner_tile_path   = struct1%carbon13%owner_tile_path
    ALLOCATE(structOut%carbon13%ptr3d,   source=struct1%carbon13%ptr3d)
    structOut%carbon14%owner_model_id    = struct1%carbon14%owner_model_id
    structOut%carbon14%owner_proc_id     = struct1%carbon14%owner_proc_id
    structOut%carbon14%owner_tile_path   = struct1%carbon14%owner_tile_path
    ALLOCATE(structOut%carbon14%ptr3d,   source=struct1%carbon14%ptr3d)
    structOut%nitrogen15%owner_model_id  = struct1%nitrogen15%owner_model_id
    structOut%nitrogen15%owner_proc_id   = struct1%nitrogen15%owner_proc_id
    structOut%nitrogen15%owner_tile_path = struct1%nitrogen15%owner_tile_path
    ALLOCATE(structOut%nitrogen15%ptr3d, source=struct1%nitrogen15%ptr3d)

    ! run calculation
    ! the differentiation between ON_CHUNK/ON_DOMAIN is done at the level of variables (i.e., the op-ov routines of the jsb_var)
    structOut%carbon          = struct1%carbon      + field2d(:,:)
    structOut%nitrogen        = struct1%nitrogen    + field2d(:,:)
    structOut%phosphorus      = struct1%phosphorus  + field2d(:,:)
    structOut%carbon13        = struct1%carbon13    + field2d(:,:)
    structOut%carbon14        = struct1%carbon14    + field2d(:,:)
    structOut%nitrogen15      = struct1%nitrogen15  + field2d(:,:)

    ! is this needed here too (not only with the jsb_var op-ov)
    structOut%carbon%ptr      => structOut%carbon%ptr3d(:,:,:)
    structOut%nitrogen%ptr    => structOut%nitrogen%ptr3d(:,:,:)
    structOut%phosphorus%ptr  => structOut%phosphorus%ptr3d(:,:,:)
    structOut%carbon13%ptr    => structOut%carbon13%ptr3d(:,:,:)
    structOut%carbon14%ptr    => structOut%carbon14%ptr3d(:,:,:)
    structOut%nitrogen15%ptr  => structOut%nitrogen15%ptr3d(:,:,:)
  
  END FUNCTION Add_sb_elem_p_and_field2d
  !--------------------------------------------------------------------------------------------------------------------------------
  ! field2d(:,:) + t_sb_elements_pool        | ON_CHUNK
  !--------------------------------------------------------------------------------------------------------------------------------
  FUNCTION Add_field2d_and_sb_elem_p(field2d, struct1) RESULT(structOut)
  
    REAL(wp),                       INTENT(in)  :: field2d(:,:)
    CLASS(t_sb_elements_pool),      INTENT(in)  :: struct1
    TYPE(t_sb_elements_pool),       TARGET      :: structOut
  
    CHARACTER(len=*), PARAMETER :: routine = modname//':Add_field2d_and_sb_elem_p'
  
    ! init temporary pool structOut
    CALL structOut%Init(name='opov_help', shortname='opov_help', l_elements=.TRUE., element_unit=TRIM(struct1%element_unit))
    structOut%carbon%owner_model_id      = struct1%carbon%owner_model_id
    structOut%carbon%owner_proc_id       = struct1%carbon%owner_proc_id
    structOut%carbon%owner_tile_path     = struct1%carbon%owner_tile_path
    ALLOCATE(structOut%carbon%ptr3d,     source=struct1%carbon%ptr3d)
    structOut%nitrogen%owner_model_id    = struct1%nitrogen%owner_model_id
    structOut%nitrogen%owner_proc_id     = struct1%nitrogen%owner_proc_id
    structOut%nitrogen%owner_tile_path   = struct1%nitrogen%owner_tile_path
    ALLOCATE(structOut%nitrogen%ptr3d,   source=struct1%nitrogen%ptr3d)
    structOut%phosphorus%owner_model_id  = struct1%phosphorus%owner_model_id
    structOut%phosphorus%owner_proc_id   = struct1%phosphorus%owner_proc_id
    structOut%phosphorus%owner_tile_path = struct1%phosphorus%owner_tile_path
    ALLOCATE(structOut%phosphorus%ptr3d, source=struct1%phosphorus%ptr3d)
    structOut%carbon13%owner_model_id    = struct1%carbon13%owner_model_id
    structOut%carbon13%owner_proc_id     = struct1%carbon13%owner_proc_id
    structOut%carbon13%owner_tile_path   = struct1%carbon13%owner_tile_path
    ALLOCATE(structOut%carbon13%ptr3d,   source=struct1%carbon13%ptr3d)
    structOut%carbon14%owner_model_id    = struct1%carbon14%owner_model_id
    structOut%carbon14%owner_proc_id     = struct1%carbon14%owner_proc_id
    structOut%carbon14%owner_tile_path   = struct1%carbon14%owner_tile_path
    ALLOCATE(structOut%carbon14%ptr3d,   source=struct1%carbon14%ptr3d)
    structOut%nitrogen15%owner_model_id  = struct1%nitrogen15%owner_model_id
    structOut%nitrogen15%owner_proc_id   = struct1%nitrogen15%owner_proc_id
    structOut%nitrogen15%owner_tile_path = struct1%nitrogen15%owner_tile_path
    ALLOCATE(structOut%nitrogen15%ptr3d, source=struct1%nitrogen15%ptr3d)

    ! run calculation
    ! the differentiation between ON_CHUNK/ON_DOMAIN is done at the level of variables (i.e., the op-ov routines of the jsb_var)
    structOut%carbon          = field2d(:,:) + struct1%carbon    
    structOut%nitrogen        = field2d(:,:) + struct1%nitrogen  
    structOut%phosphorus      = field2d(:,:) + struct1%phosphorus
    structOut%carbon13        = field2d(:,:) + struct1%carbon13  
    structOut%carbon14        = field2d(:,:) + struct1%carbon14  
    structOut%nitrogen15      = field2d(:,:) + struct1%nitrogen15

    ! is this needed here too (not only with the jsb_var op-ov)
    structOut%carbon%ptr      => structOut%carbon%ptr3d(:,:,:)
    structOut%nitrogen%ptr    => structOut%nitrogen%ptr3d(:,:,:)
    structOut%phosphorus%ptr  => structOut%phosphorus%ptr3d(:,:,:)
    structOut%carbon13%ptr    => structOut%carbon13%ptr3d(:,:,:)
    structOut%carbon14%ptr    => structOut%carbon14%ptr3d(:,:,:)
    structOut%nitrogen15%ptr  => structOut%nitrogen15%ptr3d(:,:,:)
  
  END FUNCTION Add_field2d_and_sb_elem_p
  !--------------------------------------------------------------------------------------------------------------------------------
  ! t_sb_elements_pool + field3d(:,:,:)        | ON_CHUNK & ON_DOMAIN | field3d(:,:,1) & field3d(:,:,:)
  !--------------------------------------------------------------------------------------------------------------------------------
  FUNCTION Add_sb_elem_p_and_field3d(struct1, field3d) RESULT(structOut)
  
    CLASS(t_sb_elements_pool),      INTENT(in)  :: struct1
    REAL(wp),                       INTENT(in)  :: field3d(:,:,:)
    TYPE(t_sb_elements_pool),       TARGET      :: structOut
  
    CHARACTER(len=*), PARAMETER :: routine = modname//':Add_sb_elem_p_and_field3d'
  
    ! init temporary pool structOut
    CALL structOut%Init(name='opov_help', shortname='opov_help', l_elements=.TRUE., element_unit=TRIM(struct1%element_unit))
    structOut%carbon%owner_model_id      = struct1%carbon%owner_model_id
    structOut%carbon%owner_proc_id       = struct1%carbon%owner_proc_id
    structOut%carbon%owner_tile_path     = struct1%carbon%owner_tile_path
    ALLOCATE(structOut%carbon%ptr3d,     source=struct1%carbon%ptr3d)
    structOut%nitrogen%owner_model_id    = struct1%nitrogen%owner_model_id
    structOut%nitrogen%owner_proc_id     = struct1%nitrogen%owner_proc_id
    structOut%nitrogen%owner_tile_path   = struct1%nitrogen%owner_tile_path
    ALLOCATE(structOut%nitrogen%ptr3d,   source=struct1%nitrogen%ptr3d)
    structOut%phosphorus%owner_model_id  = struct1%phosphorus%owner_model_id
    structOut%phosphorus%owner_proc_id   = struct1%phosphorus%owner_proc_id
    structOut%phosphorus%owner_tile_path = struct1%phosphorus%owner_tile_path
    ALLOCATE(structOut%phosphorus%ptr3d, source=struct1%phosphorus%ptr3d)
    structOut%carbon13%owner_model_id    = struct1%carbon13%owner_model_id
    structOut%carbon13%owner_proc_id     = struct1%carbon13%owner_proc_id
    structOut%carbon13%owner_tile_path   = struct1%carbon13%owner_tile_path
    ALLOCATE(structOut%carbon13%ptr3d,   source=struct1%carbon13%ptr3d)
    structOut%carbon14%owner_model_id    = struct1%carbon14%owner_model_id
    structOut%carbon14%owner_proc_id     = struct1%carbon14%owner_proc_id
    structOut%carbon14%owner_tile_path   = struct1%carbon14%owner_tile_path
    ALLOCATE(structOut%carbon14%ptr3d,   source=struct1%carbon14%ptr3d)
    structOut%nitrogen15%owner_model_id  = struct1%nitrogen15%owner_model_id
    structOut%nitrogen15%owner_proc_id   = struct1%nitrogen15%owner_proc_id
    structOut%nitrogen15%owner_tile_path = struct1%nitrogen15%owner_tile_path
    ALLOCATE(structOut%nitrogen15%ptr3d, source=struct1%nitrogen15%ptr3d)

    ! run calculation
    ! the differentiation between ON_CHUNK/ON_DOMAIN is done at the level of variables (i.e., the op-ov routines of the jsb_var)
    structOut%carbon          = struct1%carbon      + field3d(:,:,:)
    structOut%nitrogen        = struct1%nitrogen    + field3d(:,:,:)
    structOut%phosphorus      = struct1%phosphorus  + field3d(:,:,:)
    structOut%carbon13        = struct1%carbon13    + field3d(:,:,:)
    structOut%carbon14        = struct1%carbon14    + field3d(:,:,:)
    structOut%nitrogen15      = struct1%nitrogen15  + field3d(:,:,:)

    ! is this needed here too (not only with the jsb_var op-ov)
    structOut%carbon%ptr      => structOut%carbon%ptr3d(:,:,:)
    structOut%nitrogen%ptr    => structOut%nitrogen%ptr3d(:,:,:)
    structOut%phosphorus%ptr  => structOut%phosphorus%ptr3d(:,:,:)
    structOut%carbon13%ptr    => structOut%carbon13%ptr3d(:,:,:)
    structOut%carbon14%ptr    => structOut%carbon14%ptr3d(:,:,:)
    structOut%nitrogen15%ptr  => structOut%nitrogen15%ptr3d(:,:,:)
  
  END FUNCTION Add_sb_elem_p_and_field3d
  !--------------------------------------------------------------------------------------------------------------------------------
  ! t_sb_elements_pool + t_sb_elements_pool | ON_CHUNK & ON_DOMAIN
  !--------------------------------------------------------------------------------------------------------------------------------
  FUNCTION Add_sb_elem_p_and_sb_elem_p(struct1, struct2) RESULT(structOut)
  
    CLASS(t_sb_elements_pool),      INTENT(in)  :: struct1, struct2
    TYPE(t_sb_elements_pool),       TARGET      :: structOut
  
    CHARACTER(len=*), PARAMETER :: routine = modname//':Add_sb_elem_p_and_sb_elem_p'
  
    ! init temporary pool structOut
    CALL structOut%Init(name='opov_help', shortname='opov_help', l_elements=.TRUE., element_unit=TRIM(struct1%element_unit))
    structOut%carbon%owner_model_id      = struct1%carbon%owner_model_id
    structOut%carbon%owner_proc_id       = struct1%carbon%owner_proc_id
    structOut%carbon%owner_tile_path     = struct1%carbon%owner_tile_path
    ALLOCATE(structOut%carbon%ptr3d,     source=struct1%carbon%ptr3d)
    structOut%nitrogen%owner_model_id    = struct1%nitrogen%owner_model_id
    structOut%nitrogen%owner_proc_id     = struct1%nitrogen%owner_proc_id
    structOut%nitrogen%owner_tile_path   = struct1%nitrogen%owner_tile_path
    ALLOCATE(structOut%nitrogen%ptr3d,   source=struct1%nitrogen%ptr3d)
    structOut%phosphorus%owner_model_id  = struct1%phosphorus%owner_model_id
    structOut%phosphorus%owner_proc_id   = struct1%phosphorus%owner_proc_id
    structOut%phosphorus%owner_tile_path = struct1%phosphorus%owner_tile_path
    ALLOCATE(structOut%phosphorus%ptr3d, source=struct1%phosphorus%ptr3d)
    structOut%carbon13%owner_model_id    = struct1%carbon13%owner_model_id
    structOut%carbon13%owner_proc_id     = struct1%carbon13%owner_proc_id
    structOut%carbon13%owner_tile_path   = struct1%carbon13%owner_tile_path
    ALLOCATE(structOut%carbon13%ptr3d,   source=struct1%carbon13%ptr3d)
    structOut%carbon14%owner_model_id    = struct1%carbon14%owner_model_id
    structOut%carbon14%owner_proc_id     = struct1%carbon14%owner_proc_id
    structOut%carbon14%owner_tile_path   = struct1%carbon14%owner_tile_path
    ALLOCATE(structOut%carbon14%ptr3d,   source=struct1%carbon14%ptr3d)
    structOut%nitrogen15%owner_model_id  = struct1%nitrogen15%owner_model_id
    structOut%nitrogen15%owner_proc_id   = struct1%nitrogen15%owner_proc_id
    structOut%nitrogen15%owner_tile_path = struct1%nitrogen15%owner_tile_path
    ALLOCATE(structOut%nitrogen15%ptr3d, source=struct1%nitrogen15%ptr3d)

    ! run calculation
    ! the differentiation between ON_CHUNK/ON_DOMAIN is done at the level of variables (i.e., the op-ov routines of the jsb_var)
    structOut%carbon          = struct1%carbon      + struct2%carbon
    structOut%nitrogen        = struct1%nitrogen    + struct2%nitrogen  
    structOut%phosphorus      = struct1%phosphorus  + struct2%phosphorus
    structOut%carbon13        = struct1%carbon13    + struct2%carbon13  
    structOut%carbon14        = struct1%carbon14    + struct2%carbon14  
    structOut%nitrogen15      = struct1%nitrogen15  + struct2%nitrogen15

    ! is this needed here too (not only with the jsb_var op-ov)
    structOut%carbon%ptr      => structOut%carbon%ptr3d(:,:,:)
    structOut%nitrogen%ptr    => structOut%nitrogen%ptr3d(:,:,:)
    structOut%phosphorus%ptr  => structOut%phosphorus%ptr3d(:,:,:)
    structOut%carbon13%ptr    => structOut%carbon13%ptr3d(:,:,:)
    structOut%carbon14%ptr    => structOut%carbon14%ptr3d(:,:,:)
    structOut%nitrogen15%ptr  => structOut%nitrogen15%ptr3d(:,:,:)
  
  END FUNCTION Add_sb_elem_p_and_sb_elem_p

  !--------------------------------------------------------------------------------------------------------------------------------
  ! - minus
  !--------------------------------------------------------------------------------------------------------------------------------
  !--------------------------------------------------------------------------------------------------------------------------------
  ! t_sb_elements_pool - scalar        | ON_CHUNK & ON_DOMAIN
  !--------------------------------------------------------------------------------------------------------------------------------
  FUNCTION Subtract_sb_elem_p_minus_scalar(struct1, scalar) RESULT(structOut)
  
    CLASS(t_sb_elements_pool),      INTENT(in)  :: struct1
    REAL(wp),                       INTENT(in)  :: scalar
    TYPE(t_sb_elements_pool),       TARGET      :: structOut
  
    CHARACTER(len=*), PARAMETER :: routine = modname//':Subtract_sb_elem_p_minus_scalar'
  
    ! init temporary pool structOut
    CALL structOut%Init(name='opov_help', shortname='opov_help', l_elements=.TRUE., element_unit=TRIM(struct1%element_unit))
    structOut%carbon%owner_model_id      = struct1%carbon%owner_model_id
    structOut%carbon%owner_proc_id       = struct1%carbon%owner_proc_id
    structOut%carbon%owner_tile_path     = struct1%carbon%owner_tile_path
    ALLOCATE(structOut%carbon%ptr3d,     source=struct1%carbon%ptr3d)
    structOut%nitrogen%owner_model_id    = struct1%nitrogen%owner_model_id
    structOut%nitrogen%owner_proc_id     = struct1%nitrogen%owner_proc_id
    structOut%nitrogen%owner_tile_path   = struct1%nitrogen%owner_tile_path
    ALLOCATE(structOut%nitrogen%ptr3d,   source=struct1%nitrogen%ptr3d)
    structOut%phosphorus%owner_model_id  = struct1%phosphorus%owner_model_id
    structOut%phosphorus%owner_proc_id   = struct1%phosphorus%owner_proc_id
    structOut%phosphorus%owner_tile_path = struct1%phosphorus%owner_tile_path
    ALLOCATE(structOut%phosphorus%ptr3d, source=struct1%phosphorus%ptr3d)
    structOut%carbon13%owner_model_id    = struct1%carbon13%owner_model_id
    structOut%carbon13%owner_proc_id     = struct1%carbon13%owner_proc_id
    structOut%carbon13%owner_tile_path   = struct1%carbon13%owner_tile_path
    ALLOCATE(structOut%carbon13%ptr3d,   source=struct1%carbon13%ptr3d)
    structOut%carbon14%owner_model_id    = struct1%carbon14%owner_model_id
    structOut%carbon14%owner_proc_id     = struct1%carbon14%owner_proc_id
    structOut%carbon14%owner_tile_path   = struct1%carbon14%owner_tile_path
    ALLOCATE(structOut%carbon14%ptr3d,   source=struct1%carbon14%ptr3d)
    structOut%nitrogen15%owner_model_id  = struct1%nitrogen15%owner_model_id
    structOut%nitrogen15%owner_proc_id   = struct1%nitrogen15%owner_proc_id
    structOut%nitrogen15%owner_tile_path = struct1%nitrogen15%owner_tile_path
    ALLOCATE(structOut%nitrogen15%ptr3d, source=struct1%nitrogen15%ptr3d)

    ! run calculation
    ! the differentiation between ON_CHUNK/ON_DOMAIN is done at the level of variables (i.e., the op-ov routines of the jsb_var)
    structOut%carbon          = struct1%carbon      - scalar
    structOut%nitrogen        = struct1%nitrogen    - scalar
    structOut%phosphorus      = struct1%phosphorus  - scalar
    structOut%carbon13        = struct1%carbon13    - scalar
    structOut%carbon14        = struct1%carbon14    - scalar
    structOut%nitrogen15      = struct1%nitrogen15  - scalar

    ! is this needed here too (not only with the jsb_var op-ov)
    structOut%carbon%ptr      => structOut%carbon%ptr3d(:,:,:)
    structOut%nitrogen%ptr    => structOut%nitrogen%ptr3d(:,:,:)
    structOut%phosphorus%ptr  => structOut%phosphorus%ptr3d(:,:,:)
    structOut%carbon13%ptr    => structOut%carbon13%ptr3d(:,:,:)
    structOut%carbon14%ptr    => structOut%carbon14%ptr3d(:,:,:)
    structOut%nitrogen15%ptr  => structOut%nitrogen15%ptr3d(:,:,:)
  
  END FUNCTION Subtract_sb_elem_p_minus_scalar
  !--------------------------------------------------------------------------------------------------------------------------------
  ! t_sb_elements_pool - field2d(:,:)        | ON_CHUNK
  !--------------------------------------------------------------------------------------------------------------------------------
  FUNCTION Subtract_sb_elem_p_minus_field2d(struct1, field2d) RESULT(structOut)
  
    CLASS(t_sb_elements_pool),      INTENT(in)  :: struct1
    REAL(wp),                       INTENT(in)  :: field2d(:,:)
    TYPE(t_sb_elements_pool),       TARGET      :: structOut
  
    CHARACTER(len=*), PARAMETER :: routine = modname//':Subtract_sb_elem_p_minus_field2d'
  
    ! init temporary pool structOut
    CALL structOut%Init(name='opov_help', shortname='opov_help', l_elements=.TRUE., element_unit=TRIM(struct1%element_unit))
    structOut%carbon%owner_model_id      = struct1%carbon%owner_model_id
    structOut%carbon%owner_proc_id       = struct1%carbon%owner_proc_id
    structOut%carbon%owner_tile_path     = struct1%carbon%owner_tile_path
    ALLOCATE(structOut%carbon%ptr3d,     source=struct1%carbon%ptr3d)
    structOut%nitrogen%owner_model_id    = struct1%nitrogen%owner_model_id
    structOut%nitrogen%owner_proc_id     = struct1%nitrogen%owner_proc_id
    structOut%nitrogen%owner_tile_path   = struct1%nitrogen%owner_tile_path
    ALLOCATE(structOut%nitrogen%ptr3d,   source=struct1%nitrogen%ptr3d)
    structOut%phosphorus%owner_model_id  = struct1%phosphorus%owner_model_id
    structOut%phosphorus%owner_proc_id   = struct1%phosphorus%owner_proc_id
    structOut%phosphorus%owner_tile_path = struct1%phosphorus%owner_tile_path
    ALLOCATE(structOut%phosphorus%ptr3d, source=struct1%phosphorus%ptr3d)
    structOut%carbon13%owner_model_id    = struct1%carbon13%owner_model_id
    structOut%carbon13%owner_proc_id     = struct1%carbon13%owner_proc_id
    structOut%carbon13%owner_tile_path   = struct1%carbon13%owner_tile_path
    ALLOCATE(structOut%carbon13%ptr3d,   source=struct1%carbon13%ptr3d)
    structOut%carbon14%owner_model_id    = struct1%carbon14%owner_model_id
    structOut%carbon14%owner_proc_id     = struct1%carbon14%owner_proc_id
    structOut%carbon14%owner_tile_path   = struct1%carbon14%owner_tile_path
    ALLOCATE(structOut%carbon14%ptr3d,   source=struct1%carbon14%ptr3d)
    structOut%nitrogen15%owner_model_id  = struct1%nitrogen15%owner_model_id
    structOut%nitrogen15%owner_proc_id   = struct1%nitrogen15%owner_proc_id
    structOut%nitrogen15%owner_tile_path = struct1%nitrogen15%owner_tile_path
    ALLOCATE(structOut%nitrogen15%ptr3d, source=struct1%nitrogen15%ptr3d)

    ! run calculation
    ! the differentiation between ON_CHUNK/ON_DOMAIN is done at the level of variables (i.e., the op-ov routines of the jsb_var)
    structOut%carbon          = struct1%carbon      - field2d(:,:)
    structOut%nitrogen        = struct1%nitrogen    - field2d(:,:)
    structOut%phosphorus      = struct1%phosphorus  - field2d(:,:)
    structOut%carbon13        = struct1%carbon13    - field2d(:,:)
    structOut%carbon14        = struct1%carbon14    - field2d(:,:)
    structOut%nitrogen15      = struct1%nitrogen15  - field2d(:,:)

    ! is this needed here too (not only with the jsb_var op-ov)
    structOut%carbon%ptr      => structOut%carbon%ptr3d(:,:,:)
    structOut%nitrogen%ptr    => structOut%nitrogen%ptr3d(:,:,:)
    structOut%phosphorus%ptr  => structOut%phosphorus%ptr3d(:,:,:)
    structOut%carbon13%ptr    => structOut%carbon13%ptr3d(:,:,:)
    structOut%carbon14%ptr    => structOut%carbon14%ptr3d(:,:,:)
    structOut%nitrogen15%ptr  => structOut%nitrogen15%ptr3d(:,:,:)
  
  END FUNCTION Subtract_sb_elem_p_minus_field2d
  !--------------------------------------------------------------------------------------------------------------------------------
  ! t_sb_elements_pool - field3d(:,:,:)        | ON_CHUNK & ON_DOMAIN | field3d(:,:,1) & field3d(:,:,:)
  !--------------------------------------------------------------------------------------------------------------------------------
  FUNCTION Subtract_sb_elem_p_minus_field3d(struct1, field3d) RESULT(structOut)
  
    CLASS(t_sb_elements_pool),      INTENT(in)  :: struct1
    REAL(wp),                       INTENT(in)  :: field3d(:,:,:)
    TYPE(t_sb_elements_pool),       TARGET      :: structOut
  
    CHARACTER(len=*), PARAMETER :: routine = modname//':Subtract_sb_elem_p_minus_field3d'
  
    ! init temporary pool structOut
    CALL structOut%Init(name='opov_help', shortname='opov_help', l_elements=.TRUE., element_unit=TRIM(struct1%element_unit))
    structOut%carbon%owner_model_id      = struct1%carbon%owner_model_id
    structOut%carbon%owner_proc_id       = struct1%carbon%owner_proc_id
    structOut%carbon%owner_tile_path     = struct1%carbon%owner_tile_path
    ALLOCATE(structOut%carbon%ptr3d,     source=struct1%carbon%ptr3d)
    structOut%nitrogen%owner_model_id    = struct1%nitrogen%owner_model_id
    structOut%nitrogen%owner_proc_id     = struct1%nitrogen%owner_proc_id
    structOut%nitrogen%owner_tile_path   = struct1%nitrogen%owner_tile_path
    ALLOCATE(structOut%nitrogen%ptr3d,   source=struct1%nitrogen%ptr3d)
    structOut%phosphorus%owner_model_id  = struct1%phosphorus%owner_model_id
    structOut%phosphorus%owner_proc_id   = struct1%phosphorus%owner_proc_id
    structOut%phosphorus%owner_tile_path = struct1%phosphorus%owner_tile_path
    ALLOCATE(structOut%phosphorus%ptr3d, source=struct1%phosphorus%ptr3d)
    structOut%carbon13%owner_model_id    = struct1%carbon13%owner_model_id
    structOut%carbon13%owner_proc_id     = struct1%carbon13%owner_proc_id
    structOut%carbon13%owner_tile_path   = struct1%carbon13%owner_tile_path
    ALLOCATE(structOut%carbon13%ptr3d,   source=struct1%carbon13%ptr3d)
    structOut%carbon14%owner_model_id    = struct1%carbon14%owner_model_id
    structOut%carbon14%owner_proc_id     = struct1%carbon14%owner_proc_id
    structOut%carbon14%owner_tile_path   = struct1%carbon14%owner_tile_path
    ALLOCATE(structOut%carbon14%ptr3d,   source=struct1%carbon14%ptr3d)
    structOut%nitrogen15%owner_model_id  = struct1%nitrogen15%owner_model_id
    structOut%nitrogen15%owner_proc_id   = struct1%nitrogen15%owner_proc_id
    structOut%nitrogen15%owner_tile_path = struct1%nitrogen15%owner_tile_path
    ALLOCATE(structOut%nitrogen15%ptr3d, source=struct1%nitrogen15%ptr3d)

    ! run calculation
    ! the differentiation between ON_CHUNK/ON_DOMAIN is done at the level of variables (i.e., the op-ov routines of the jsb_var)
    structOut%carbon          = struct1%carbon      - field3d(:,:,:)
    structOut%nitrogen        = struct1%nitrogen    - field3d(:,:,:)
    structOut%phosphorus      = struct1%phosphorus  - field3d(:,:,:)
    structOut%carbon13        = struct1%carbon13    - field3d(:,:,:)
    structOut%carbon14        = struct1%carbon14    - field3d(:,:,:)
    structOut%nitrogen15      = struct1%nitrogen15  - field3d(:,:,:)

    ! is this needed here too (not only with the jsb_var op-ov)
    structOut%carbon%ptr      => structOut%carbon%ptr3d(:,:,:)
    structOut%nitrogen%ptr    => structOut%nitrogen%ptr3d(:,:,:)
    structOut%phosphorus%ptr  => structOut%phosphorus%ptr3d(:,:,:)
    structOut%carbon13%ptr    => structOut%carbon13%ptr3d(:,:,:)
    structOut%carbon14%ptr    => structOut%carbon14%ptr3d(:,:,:)
    structOut%nitrogen15%ptr  => structOut%nitrogen15%ptr3d(:,:,:)
  
  END FUNCTION Subtract_sb_elem_p_minus_field3d
  !--------------------------------------------------------------------------------------------------------------------------------
  ! t_sb_elements_pool - t_sb_elements_pool | ON_CHUNK & ON_DOMAIN
  !--------------------------------------------------------------------------------------------------------------------------------
  FUNCTION Subtract_sb_elem_p_minus_sb_elem_p(struct1, struct2) RESULT(structOut)
  
    CLASS(t_sb_elements_pool),      INTENT(in)  :: struct1, struct2
    TYPE(t_sb_elements_pool),       TARGET      :: structOut
  
    CHARACTER(len=*), PARAMETER :: routine = modname//':Subtract_sb_elem_p_minus_sb_elem_p'
  
    ! init temporary pool structOut
    CALL structOut%Init(name='opov_help', shortname='opov_help', l_elements=.TRUE., element_unit=TRIM(struct1%element_unit))
    structOut%carbon%owner_model_id      = struct1%carbon%owner_model_id
    structOut%carbon%owner_proc_id       = struct1%carbon%owner_proc_id
    structOut%carbon%owner_tile_path     = struct1%carbon%owner_tile_path
    ALLOCATE(structOut%carbon%ptr3d,     source=struct1%carbon%ptr3d)
    structOut%nitrogen%owner_model_id    = struct1%nitrogen%owner_model_id
    structOut%nitrogen%owner_proc_id     = struct1%nitrogen%owner_proc_id
    structOut%nitrogen%owner_tile_path   = struct1%nitrogen%owner_tile_path
    ALLOCATE(structOut%nitrogen%ptr3d,   source=struct1%nitrogen%ptr3d)
    structOut%phosphorus%owner_model_id  = struct1%phosphorus%owner_model_id
    structOut%phosphorus%owner_proc_id   = struct1%phosphorus%owner_proc_id
    structOut%phosphorus%owner_tile_path = struct1%phosphorus%owner_tile_path
    ALLOCATE(structOut%phosphorus%ptr3d, source=struct1%phosphorus%ptr3d)
    structOut%carbon13%owner_model_id    = struct1%carbon13%owner_model_id
    structOut%carbon13%owner_proc_id     = struct1%carbon13%owner_proc_id
    structOut%carbon13%owner_tile_path   = struct1%carbon13%owner_tile_path
    ALLOCATE(structOut%carbon13%ptr3d,   source=struct1%carbon13%ptr3d)
    structOut%carbon14%owner_model_id    = struct1%carbon14%owner_model_id
    structOut%carbon14%owner_proc_id     = struct1%carbon14%owner_proc_id
    structOut%carbon14%owner_tile_path   = struct1%carbon14%owner_tile_path
    ALLOCATE(structOut%carbon14%ptr3d,   source=struct1%carbon14%ptr3d)
    structOut%nitrogen15%owner_model_id  = struct1%nitrogen15%owner_model_id
    structOut%nitrogen15%owner_proc_id   = struct1%nitrogen15%owner_proc_id
    structOut%nitrogen15%owner_tile_path = struct1%nitrogen15%owner_tile_path
    ALLOCATE(structOut%nitrogen15%ptr3d, source=struct1%nitrogen15%ptr3d)

    ! run calculation
    ! the differentiation between ON_CHUNK/ON_DOMAIN is done at the level of variables (i.e., the op-ov routines of the jsb_var)
    structOut%carbon          = struct1%carbon      - struct2%carbon
    structOut%nitrogen        = struct1%nitrogen    - struct2%nitrogen  
    structOut%phosphorus      = struct1%phosphorus  - struct2%phosphorus
    structOut%carbon13        = struct1%carbon13    - struct2%carbon13  
    structOut%carbon14        = struct1%carbon14    - struct2%carbon14  
    structOut%nitrogen15      = struct1%nitrogen15  - struct2%nitrogen15

    ! is this needed here too (not only with the jsb_var op-ov)
    structOut%carbon%ptr      => structOut%carbon%ptr3d(:,:,:)
    structOut%nitrogen%ptr    => structOut%nitrogen%ptr3d(:,:,:)
    structOut%phosphorus%ptr  => structOut%phosphorus%ptr3d(:,:,:)
    structOut%carbon13%ptr    => structOut%carbon13%ptr3d(:,:,:)
    structOut%carbon14%ptr    => structOut%carbon14%ptr3d(:,:,:)
    structOut%nitrogen15%ptr  => structOut%nitrogen15%ptr3d(:,:,:)
  
  END FUNCTION Subtract_sb_elem_p_minus_sb_elem_p

  !--------------------------------------------------------------------------------------------------------------------------------
  ! * multiply
  !--------------------------------------------------------------------------------------------------------------------------------
  !--------------------------------------------------------------------------------------------------------------------------------
  ! t_sb_elements_pool * scalar        | ON_CHUNK & ON_DOMAIN
  !--------------------------------------------------------------------------------------------------------------------------------
  FUNCTION Multiply_sb_elem_p_with_scalar(struct1, scalar) RESULT(structOut)
  
    CLASS(t_sb_elements_pool),      INTENT(in)  :: struct1
    REAL(wp),                       INTENT(in)  :: scalar
    TYPE(t_sb_elements_pool),       TARGET      :: structOut
  
    CHARACTER(len=*), PARAMETER :: routine = modname//':Multiply_sb_elem_p_with_scalar'
  
    ! init temporary pool structOut
    CALL structOut%Init(name='opov_help', shortname='opov_help', l_elements=.TRUE., element_unit=TRIM(struct1%element_unit))
    structOut%carbon%owner_model_id      = struct1%carbon%owner_model_id
    structOut%carbon%owner_proc_id       = struct1%carbon%owner_proc_id
    structOut%carbon%owner_tile_path     = struct1%carbon%owner_tile_path
    ALLOCATE(structOut%carbon%ptr3d,     source=struct1%carbon%ptr3d)
    structOut%nitrogen%owner_model_id    = struct1%nitrogen%owner_model_id
    structOut%nitrogen%owner_proc_id     = struct1%nitrogen%owner_proc_id
    structOut%nitrogen%owner_tile_path   = struct1%nitrogen%owner_tile_path
    ALLOCATE(structOut%nitrogen%ptr3d,   source=struct1%nitrogen%ptr3d)
    structOut%phosphorus%owner_model_id  = struct1%phosphorus%owner_model_id
    structOut%phosphorus%owner_proc_id   = struct1%phosphorus%owner_proc_id
    structOut%phosphorus%owner_tile_path = struct1%phosphorus%owner_tile_path
    ALLOCATE(structOut%phosphorus%ptr3d, source=struct1%phosphorus%ptr3d)
    structOut%carbon13%owner_model_id    = struct1%carbon13%owner_model_id
    structOut%carbon13%owner_proc_id     = struct1%carbon13%owner_proc_id
    structOut%carbon13%owner_tile_path   = struct1%carbon13%owner_tile_path
    ALLOCATE(structOut%carbon13%ptr3d,   source=struct1%carbon13%ptr3d)
    structOut%carbon14%owner_model_id    = struct1%carbon14%owner_model_id
    structOut%carbon14%owner_proc_id     = struct1%carbon14%owner_proc_id
    structOut%carbon14%owner_tile_path   = struct1%carbon14%owner_tile_path
    ALLOCATE(structOut%carbon14%ptr3d,   source=struct1%carbon14%ptr3d)
    structOut%nitrogen15%owner_model_id  = struct1%nitrogen15%owner_model_id
    structOut%nitrogen15%owner_proc_id   = struct1%nitrogen15%owner_proc_id
    structOut%nitrogen15%owner_tile_path = struct1%nitrogen15%owner_tile_path
    ALLOCATE(structOut%nitrogen15%ptr3d, source=struct1%nitrogen15%ptr3d)

    ! run calculation
    ! the differentiation between ON_CHUNK/ON_DOMAIN is done at the level of variables (i.e., the op-ov routines of the jsb_var)
    structOut%carbon          = struct1%carbon      * scalar
    structOut%nitrogen        = struct1%nitrogen    * scalar
    structOut%phosphorus      = struct1%phosphorus  * scalar
    structOut%carbon13        = struct1%carbon13    * scalar
    structOut%carbon14        = struct1%carbon14    * scalar
    structOut%nitrogen15      = struct1%nitrogen15  * scalar

    ! is this needed here too (not only with the jsb_var op-ov)
    structOut%carbon%ptr      => structOut%carbon%ptr3d(:,:,:)
    structOut%nitrogen%ptr    => structOut%nitrogen%ptr3d(:,:,:)
    structOut%phosphorus%ptr  => structOut%phosphorus%ptr3d(:,:,:)
    structOut%carbon13%ptr    => structOut%carbon13%ptr3d(:,:,:)
    structOut%carbon14%ptr    => structOut%carbon14%ptr3d(:,:,:)
    structOut%nitrogen15%ptr  => structOut%nitrogen15%ptr3d(:,:,:)
  
  END FUNCTION Multiply_sb_elem_p_with_scalar
  !--------------------------------------------------------------------------------------------------------------------------------
  ! scalar * t_sb_elements_pool        | ON_CHUNK & ON_DOMAIN
  !--------------------------------------------------------------------------------------------------------------------------------
  FUNCTION Multiply_scalar_with_sb_elem_p(scalar, struct1) RESULT(structOut)
  
    REAL(wp),                       INTENT(in)  :: scalar
    CLASS(t_sb_elements_pool),      INTENT(in)  :: struct1
    TYPE(t_sb_elements_pool),       TARGET      :: structOut
  
    CHARACTER(len=*), PARAMETER :: routine = modname//':Multiply_scalar_with_sb_elem_p'
  
    ! init temporary pool structOut
    CALL structOut%Init(name='opov_help', shortname='opov_help', l_elements=.TRUE., element_unit=TRIM(struct1%element_unit))
    structOut%carbon%owner_model_id      = struct1%carbon%owner_model_id
    structOut%carbon%owner_proc_id       = struct1%carbon%owner_proc_id
    structOut%carbon%owner_tile_path     = struct1%carbon%owner_tile_path
    ALLOCATE(structOut%carbon%ptr3d,     source=struct1%carbon%ptr3d)
    structOut%nitrogen%owner_model_id    = struct1%nitrogen%owner_model_id
    structOut%nitrogen%owner_proc_id     = struct1%nitrogen%owner_proc_id
    structOut%nitrogen%owner_tile_path   = struct1%nitrogen%owner_tile_path
    ALLOCATE(structOut%nitrogen%ptr3d,   source=struct1%nitrogen%ptr3d)
    structOut%phosphorus%owner_model_id  = struct1%phosphorus%owner_model_id
    structOut%phosphorus%owner_proc_id   = struct1%phosphorus%owner_proc_id
    structOut%phosphorus%owner_tile_path = struct1%phosphorus%owner_tile_path
    ALLOCATE(structOut%phosphorus%ptr3d, source=struct1%phosphorus%ptr3d)
    structOut%carbon13%owner_model_id    = struct1%carbon13%owner_model_id
    structOut%carbon13%owner_proc_id     = struct1%carbon13%owner_proc_id
    structOut%carbon13%owner_tile_path   = struct1%carbon13%owner_tile_path
    ALLOCATE(structOut%carbon13%ptr3d,   source=struct1%carbon13%ptr3d)
    structOut%carbon14%owner_model_id    = struct1%carbon14%owner_model_id
    structOut%carbon14%owner_proc_id     = struct1%carbon14%owner_proc_id
    structOut%carbon14%owner_tile_path   = struct1%carbon14%owner_tile_path
    ALLOCATE(structOut%carbon14%ptr3d,   source=struct1%carbon14%ptr3d)
    structOut%nitrogen15%owner_model_id  = struct1%nitrogen15%owner_model_id
    structOut%nitrogen15%owner_proc_id   = struct1%nitrogen15%owner_proc_id
    structOut%nitrogen15%owner_tile_path = struct1%nitrogen15%owner_tile_path
    ALLOCATE(structOut%nitrogen15%ptr3d, source=struct1%nitrogen15%ptr3d)

    ! run calculation
    ! the differentiation between ON_CHUNK/ON_DOMAIN is done at the level of variables (i.e., the op-ov routines of the jsb_var)
    structOut%carbon          = scalar * struct1%carbon    
    structOut%nitrogen        = scalar * struct1%nitrogen  
    structOut%phosphorus      = scalar * struct1%phosphorus
    structOut%carbon13        = scalar * struct1%carbon13  
    structOut%carbon14        = scalar * struct1%carbon14  
    structOut%nitrogen15      = scalar * struct1%nitrogen15

    ! is this needed here too (not only with the jsb_var op-ov)
    structOut%carbon%ptr      => structOut%carbon%ptr3d(:,:,:)
    structOut%nitrogen%ptr    => structOut%nitrogen%ptr3d(:,:,:)
    structOut%phosphorus%ptr  => structOut%phosphorus%ptr3d(:,:,:)
    structOut%carbon13%ptr    => structOut%carbon13%ptr3d(:,:,:)
    structOut%carbon14%ptr    => structOut%carbon14%ptr3d(:,:,:)
    structOut%nitrogen15%ptr  => structOut%nitrogen15%ptr3d(:,:,:)
  
  END FUNCTION Multiply_scalar_with_sb_elem_p
  !--------------------------------------------------------------------------------------------------------------------------------
  ! t_sb_elements_pool * field2d(:,:)        | ON_CHUNK
  !--------------------------------------------------------------------------------------------------------------------------------
  FUNCTION Multiply_sb_elem_p_with_field2d(struct1, field2d) RESULT(structOut)
  
    CLASS(t_sb_elements_pool),      INTENT(in)  :: struct1
    REAL(wp),                       INTENT(in)  :: field2d(:,:)
    TYPE(t_sb_elements_pool),       TARGET      :: structOut
  
    CHARACTER(len=*), PARAMETER :: routine = modname//':Multiply_sb_elem_p_with_field2d'
  
    ! init temporary pool structOut
    CALL structOut%Init(name='opov_help', shortname='opov_help', l_elements=.TRUE., element_unit=TRIM(struct1%element_unit))
    structOut%carbon%owner_model_id      = struct1%carbon%owner_model_id
    structOut%carbon%owner_proc_id       = struct1%carbon%owner_proc_id
    structOut%carbon%owner_tile_path     = struct1%carbon%owner_tile_path
    ALLOCATE(structOut%carbon%ptr3d,     source=struct1%carbon%ptr3d)
    structOut%nitrogen%owner_model_id    = struct1%nitrogen%owner_model_id
    structOut%nitrogen%owner_proc_id     = struct1%nitrogen%owner_proc_id
    structOut%nitrogen%owner_tile_path   = struct1%nitrogen%owner_tile_path
    ALLOCATE(structOut%nitrogen%ptr3d,   source=struct1%nitrogen%ptr3d)
    structOut%phosphorus%owner_model_id  = struct1%phosphorus%owner_model_id
    structOut%phosphorus%owner_proc_id   = struct1%phosphorus%owner_proc_id
    structOut%phosphorus%owner_tile_path = struct1%phosphorus%owner_tile_path
    ALLOCATE(structOut%phosphorus%ptr3d, source=struct1%phosphorus%ptr3d)
    structOut%carbon13%owner_model_id    = struct1%carbon13%owner_model_id
    structOut%carbon13%owner_proc_id     = struct1%carbon13%owner_proc_id
    structOut%carbon13%owner_tile_path   = struct1%carbon13%owner_tile_path
    ALLOCATE(structOut%carbon13%ptr3d,   source=struct1%carbon13%ptr3d)
    structOut%carbon14%owner_model_id    = struct1%carbon14%owner_model_id
    structOut%carbon14%owner_proc_id     = struct1%carbon14%owner_proc_id
    structOut%carbon14%owner_tile_path   = struct1%carbon14%owner_tile_path
    ALLOCATE(structOut%carbon14%ptr3d,   source=struct1%carbon14%ptr3d)
    structOut%nitrogen15%owner_model_id  = struct1%nitrogen15%owner_model_id
    structOut%nitrogen15%owner_proc_id   = struct1%nitrogen15%owner_proc_id
    structOut%nitrogen15%owner_tile_path = struct1%nitrogen15%owner_tile_path
    ALLOCATE(structOut%nitrogen15%ptr3d, source=struct1%nitrogen15%ptr3d)

    ! run calculation
    ! the differentiation between ON_CHUNK/ON_DOMAIN is done at the level of variables (i.e., the op-ov routines of the jsb_var)
    structOut%carbon          = struct1%carbon      * field2d(:,:)
    structOut%nitrogen        = struct1%nitrogen    * field2d(:,:)
    structOut%phosphorus      = struct1%phosphorus  * field2d(:,:)
    structOut%carbon13        = struct1%carbon13    * field2d(:,:)
    structOut%carbon14        = struct1%carbon14    * field2d(:,:)
    structOut%nitrogen15      = struct1%nitrogen15  * field2d(:,:)

    ! is this needed here too (not only with the jsb_var op-ov)
    structOut%carbon%ptr      => structOut%carbon%ptr3d(:,:,:)
    structOut%nitrogen%ptr    => structOut%nitrogen%ptr3d(:,:,:)
    structOut%phosphorus%ptr  => structOut%phosphorus%ptr3d(:,:,:)
    structOut%carbon13%ptr    => structOut%carbon13%ptr3d(:,:,:)
    structOut%carbon14%ptr    => structOut%carbon14%ptr3d(:,:,:)
    structOut%nitrogen15%ptr  => structOut%nitrogen15%ptr3d(:,:,:)
  
  END FUNCTION Multiply_sb_elem_p_with_field2d
  !--------------------------------------------------------------------------------------------------------------------------------
  ! field2d(:,:) * t_sb_elements_pool        | ON_CHUNK
  !--------------------------------------------------------------------------------------------------------------------------------
  FUNCTION Multiply_field2d_with_sb_elem_p(field2d, struct1) RESULT(structOut)
  
    REAL(wp),                       INTENT(in)  :: field2d(:,:)
    CLASS(t_sb_elements_pool),      INTENT(in)  :: struct1
    TYPE(t_sb_elements_pool),       TARGET      :: structOut
  
    CHARACTER(len=*), PARAMETER :: routine = modname//':Multiply_field2d_with_sb_elem_p'
  
    ! init temporary pool structOut
    CALL structOut%Init(name='opov_help', shortname='opov_help', l_elements=.TRUE., element_unit=TRIM(struct1%element_unit))
    structOut%carbon%owner_model_id      = struct1%carbon%owner_model_id
    structOut%carbon%owner_proc_id       = struct1%carbon%owner_proc_id
    structOut%carbon%owner_tile_path     = struct1%carbon%owner_tile_path
    ALLOCATE(structOut%carbon%ptr3d,     source=struct1%carbon%ptr3d)
    structOut%nitrogen%owner_model_id    = struct1%nitrogen%owner_model_id
    structOut%nitrogen%owner_proc_id     = struct1%nitrogen%owner_proc_id
    structOut%nitrogen%owner_tile_path   = struct1%nitrogen%owner_tile_path
    ALLOCATE(structOut%nitrogen%ptr3d,   source=struct1%nitrogen%ptr3d)
    structOut%phosphorus%owner_model_id  = struct1%phosphorus%owner_model_id
    structOut%phosphorus%owner_proc_id   = struct1%phosphorus%owner_proc_id
    structOut%phosphorus%owner_tile_path = struct1%phosphorus%owner_tile_path
    ALLOCATE(structOut%phosphorus%ptr3d, source=struct1%phosphorus%ptr3d)
    structOut%carbon13%owner_model_id    = struct1%carbon13%owner_model_id
    structOut%carbon13%owner_proc_id     = struct1%carbon13%owner_proc_id
    structOut%carbon13%owner_tile_path   = struct1%carbon13%owner_tile_path
    ALLOCATE(structOut%carbon13%ptr3d,   source=struct1%carbon13%ptr3d)
    structOut%carbon14%owner_model_id    = struct1%carbon14%owner_model_id
    structOut%carbon14%owner_proc_id     = struct1%carbon14%owner_proc_id
    structOut%carbon14%owner_tile_path   = struct1%carbon14%owner_tile_path
    ALLOCATE(structOut%carbon14%ptr3d,   source=struct1%carbon14%ptr3d)
    structOut%nitrogen15%owner_model_id  = struct1%nitrogen15%owner_model_id
    structOut%nitrogen15%owner_proc_id   = struct1%nitrogen15%owner_proc_id
    structOut%nitrogen15%owner_tile_path = struct1%nitrogen15%owner_tile_path
    ALLOCATE(structOut%nitrogen15%ptr3d, source=struct1%nitrogen15%ptr3d)

    ! run calculation
    ! the differentiation between ON_CHUNK/ON_DOMAIN is done at the level of variables (i.e., the op-ov routines of the jsb_var)
    structOut%carbon          = field2d(:,:) * struct1%carbon    
    structOut%nitrogen        = field2d(:,:) * struct1%nitrogen  
    structOut%phosphorus      = field2d(:,:) * struct1%phosphorus
    structOut%carbon13        = field2d(:,:) * struct1%carbon13  
    structOut%carbon14        = field2d(:,:) * struct1%carbon14  
    structOut%nitrogen15      = field2d(:,:) * struct1%nitrogen15

    ! is this needed here too (not only with the jsb_var op-ov)
    structOut%carbon%ptr      => structOut%carbon%ptr3d(:,:,:)
    structOut%nitrogen%ptr    => structOut%nitrogen%ptr3d(:,:,:)
    structOut%phosphorus%ptr  => structOut%phosphorus%ptr3d(:,:,:)
    structOut%carbon13%ptr    => structOut%carbon13%ptr3d(:,:,:)
    structOut%carbon14%ptr    => structOut%carbon14%ptr3d(:,:,:)
    structOut%nitrogen15%ptr  => structOut%nitrogen15%ptr3d(:,:,:)
  
  END FUNCTION Multiply_field2d_with_sb_elem_p
  !--------------------------------------------------------------------------------------------------------------------------------
  ! t_sb_elements_pool * field3d(:,:,:)        | ON_CHUNK & ON_DOMAIN | field3d(:,:,1) & field3d(:,:,:)
  !--------------------------------------------------------------------------------------------------------------------------------
  FUNCTION Multiply_sb_elem_p_with_field3d(struct1, field3d) RESULT(structOut)
  
    CLASS(t_sb_elements_pool),      INTENT(in)  :: struct1
    REAL(wp),                       INTENT(in)  :: field3d(:,:,:)
    TYPE(t_sb_elements_pool),       TARGET      :: structOut
  
    CHARACTER(len=*), PARAMETER :: routine = modname//':Multiply_sb_elem_p_with_field3d'
  
    ! init temporary pool structOut
    CALL structOut%Init(name='opov_help', shortname='opov_help', l_elements=.TRUE., element_unit=TRIM(struct1%element_unit))
    structOut%carbon%owner_model_id      = struct1%carbon%owner_model_id
    structOut%carbon%owner_proc_id       = struct1%carbon%owner_proc_id
    structOut%carbon%owner_tile_path     = struct1%carbon%owner_tile_path
    ALLOCATE(structOut%carbon%ptr3d,     source=struct1%carbon%ptr3d)
    structOut%nitrogen%owner_model_id    = struct1%nitrogen%owner_model_id
    structOut%nitrogen%owner_proc_id     = struct1%nitrogen%owner_proc_id
    structOut%nitrogen%owner_tile_path   = struct1%nitrogen%owner_tile_path
    ALLOCATE(structOut%nitrogen%ptr3d,   source=struct1%nitrogen%ptr3d)
    structOut%phosphorus%owner_model_id  = struct1%phosphorus%owner_model_id
    structOut%phosphorus%owner_proc_id   = struct1%phosphorus%owner_proc_id
    structOut%phosphorus%owner_tile_path = struct1%phosphorus%owner_tile_path
    ALLOCATE(structOut%phosphorus%ptr3d, source=struct1%phosphorus%ptr3d)
    structOut%carbon13%owner_model_id    = struct1%carbon13%owner_model_id
    structOut%carbon13%owner_proc_id     = struct1%carbon13%owner_proc_id
    structOut%carbon13%owner_tile_path   = struct1%carbon13%owner_tile_path
    ALLOCATE(structOut%carbon13%ptr3d,   source=struct1%carbon13%ptr3d)
    structOut%carbon14%owner_model_id    = struct1%carbon14%owner_model_id
    structOut%carbon14%owner_proc_id     = struct1%carbon14%owner_proc_id
    structOut%carbon14%owner_tile_path   = struct1%carbon14%owner_tile_path
    ALLOCATE(structOut%carbon14%ptr3d,   source=struct1%carbon14%ptr3d)
    structOut%nitrogen15%owner_model_id  = struct1%nitrogen15%owner_model_id
    structOut%nitrogen15%owner_proc_id   = struct1%nitrogen15%owner_proc_id
    structOut%nitrogen15%owner_tile_path = struct1%nitrogen15%owner_tile_path
    ALLOCATE(structOut%nitrogen15%ptr3d, source=struct1%nitrogen15%ptr3d)

    ! run calculation
    ! the differentiation between ON_CHUNK/ON_DOMAIN is done at the level of variables (i.e., the op-ov routines of the jsb_var)
    structOut%carbon          = struct1%carbon      * field3d(:,:,:)
    structOut%nitrogen        = struct1%nitrogen    * field3d(:,:,:)
    structOut%phosphorus      = struct1%phosphorus  * field3d(:,:,:)
    structOut%carbon13        = struct1%carbon13    * field3d(:,:,:)
    structOut%carbon14        = struct1%carbon14    * field3d(:,:,:)
    structOut%nitrogen15      = struct1%nitrogen15  * field3d(:,:,:)

    ! is this needed here too (not only with the jsb_var op-ov)
    structOut%carbon%ptr      => structOut%carbon%ptr3d(:,:,:)
    structOut%nitrogen%ptr    => structOut%nitrogen%ptr3d(:,:,:)
    structOut%phosphorus%ptr  => structOut%phosphorus%ptr3d(:,:,:)
    structOut%carbon13%ptr    => structOut%carbon13%ptr3d(:,:,:)
    structOut%carbon14%ptr    => structOut%carbon14%ptr3d(:,:,:)
    structOut%nitrogen15%ptr  => structOut%nitrogen15%ptr3d(:,:,:)
  
  END FUNCTION Multiply_sb_elem_p_with_field3d
  !--------------------------------------------------------------------------------------------------------------------------------
  ! field3d(:,:,:) * t_sb_elements_pool        | ON_CHUNK & ON_DOMAIN | field3d(:,:,1) & field3d(:,:,:)
  !--------------------------------------------------------------------------------------------------------------------------------
  FUNCTION Multiply_field3d_with_sb_elem_p(field3d, struct1) RESULT(structOut)
  
    REAL(wp),                       INTENT(in)  :: field3d(:,:,:)
    CLASS(t_sb_elements_pool),      INTENT(in)  :: struct1
    TYPE(t_sb_elements_pool),       TARGET      :: structOut
  
    CHARACTER(len=*), PARAMETER :: routine = modname//':Multiply_field3d_with_sb_elem_p'
  
    ! init temporary pool structOut
    CALL structOut%Init(name='opov_help', shortname='opov_help', l_elements=.TRUE., element_unit=TRIM(struct1%element_unit))
    structOut%carbon%owner_model_id      = struct1%carbon%owner_model_id
    structOut%carbon%owner_proc_id       = struct1%carbon%owner_proc_id
    structOut%carbon%owner_tile_path     = struct1%carbon%owner_tile_path
    ALLOCATE(structOut%carbon%ptr3d,     source=struct1%carbon%ptr3d)
    structOut%nitrogen%owner_model_id    = struct1%nitrogen%owner_model_id
    structOut%nitrogen%owner_proc_id     = struct1%nitrogen%owner_proc_id
    structOut%nitrogen%owner_tile_path   = struct1%nitrogen%owner_tile_path
    ALLOCATE(structOut%nitrogen%ptr3d,   source=struct1%nitrogen%ptr3d)
    structOut%phosphorus%owner_model_id  = struct1%phosphorus%owner_model_id
    structOut%phosphorus%owner_proc_id   = struct1%phosphorus%owner_proc_id
    structOut%phosphorus%owner_tile_path = struct1%phosphorus%owner_tile_path
    ALLOCATE(structOut%phosphorus%ptr3d, source=struct1%phosphorus%ptr3d)
    structOut%carbon13%owner_model_id    = struct1%carbon13%owner_model_id
    structOut%carbon13%owner_proc_id     = struct1%carbon13%owner_proc_id
    structOut%carbon13%owner_tile_path   = struct1%carbon13%owner_tile_path
    ALLOCATE(structOut%carbon13%ptr3d,   source=struct1%carbon13%ptr3d)
    structOut%carbon14%owner_model_id    = struct1%carbon14%owner_model_id
    structOut%carbon14%owner_proc_id     = struct1%carbon14%owner_proc_id
    structOut%carbon14%owner_tile_path   = struct1%carbon14%owner_tile_path
    ALLOCATE(structOut%carbon14%ptr3d,   source=struct1%carbon14%ptr3d)
    structOut%nitrogen15%owner_model_id  = struct1%nitrogen15%owner_model_id
    structOut%nitrogen15%owner_proc_id   = struct1%nitrogen15%owner_proc_id
    structOut%nitrogen15%owner_tile_path = struct1%nitrogen15%owner_tile_path
    ALLOCATE(structOut%nitrogen15%ptr3d, source=struct1%nitrogen15%ptr3d)

    ! run calculation
    ! the differentiation between ON_CHUNK/ON_DOMAIN is done at the level of variables (i.e., the op-ov routines of the jsb_var)
    structOut%carbon          = field3d(:,:,:) * struct1%carbon    
    structOut%nitrogen        = field3d(:,:,:) * struct1%nitrogen  
    structOut%phosphorus      = field3d(:,:,:) * struct1%phosphorus
    structOut%carbon13        = field3d(:,:,:) * struct1%carbon13  
    structOut%carbon14        = field3d(:,:,:) * struct1%carbon14  
    structOut%nitrogen15      = field3d(:,:,:) * struct1%nitrogen15

    ! is this needed here too (not only with the jsb_var op-ov)
    structOut%carbon%ptr      => structOut%carbon%ptr3d(:,:,:)
    structOut%nitrogen%ptr    => structOut%nitrogen%ptr3d(:,:,:)
    structOut%phosphorus%ptr  => structOut%phosphorus%ptr3d(:,:,:)
    structOut%carbon13%ptr    => structOut%carbon13%ptr3d(:,:,:)
    structOut%carbon14%ptr    => structOut%carbon14%ptr3d(:,:,:)
    structOut%nitrogen15%ptr  => structOut%nitrogen15%ptr3d(:,:,:)
  
  END FUNCTION Multiply_field3d_with_sb_elem_p
  !--------------------------------------------------------------------------------------------------------------------------------
  ! t_sb_elements_pool * t_sb_elements_pool | ON_CHUNK & ON_DOMAIN
  !--------------------------------------------------------------------------------------------------------------------------------
  FUNCTION Multiply_sb_elem_p_with_sb_elem_p(struct1, struct2) RESULT(structOut)
  
    CLASS(t_sb_elements_pool),      INTENT(in)  :: struct1, struct2
    TYPE(t_sb_elements_pool),       TARGET      :: structOut
  
    CHARACTER(len=*), PARAMETER :: routine = modname//':Multiply_sb_elem_p_with_sb_elem_p'
  
    ! init temporary pool structOut
    CALL structOut%Init(name='opov_help', shortname='opov_help', l_elements=.TRUE., element_unit=TRIM(struct1%element_unit))
    structOut%carbon%owner_model_id      = struct1%carbon%owner_model_id
    structOut%carbon%owner_proc_id       = struct1%carbon%owner_proc_id
    structOut%carbon%owner_tile_path     = struct1%carbon%owner_tile_path
    ALLOCATE(structOut%carbon%ptr3d,     source=struct1%carbon%ptr3d)
    structOut%nitrogen%owner_model_id    = struct1%nitrogen%owner_model_id
    structOut%nitrogen%owner_proc_id     = struct1%nitrogen%owner_proc_id
    structOut%nitrogen%owner_tile_path   = struct1%nitrogen%owner_tile_path
    ALLOCATE(structOut%nitrogen%ptr3d,   source=struct1%nitrogen%ptr3d)
    structOut%phosphorus%owner_model_id  = struct1%phosphorus%owner_model_id
    structOut%phosphorus%owner_proc_id   = struct1%phosphorus%owner_proc_id
    structOut%phosphorus%owner_tile_path = struct1%phosphorus%owner_tile_path
    ALLOCATE(structOut%phosphorus%ptr3d, source=struct1%phosphorus%ptr3d)
    structOut%carbon13%owner_model_id    = struct1%carbon13%owner_model_id
    structOut%carbon13%owner_proc_id     = struct1%carbon13%owner_proc_id
    structOut%carbon13%owner_tile_path   = struct1%carbon13%owner_tile_path
    ALLOCATE(structOut%carbon13%ptr3d,   source=struct1%carbon13%ptr3d)
    structOut%carbon14%owner_model_id    = struct1%carbon14%owner_model_id
    structOut%carbon14%owner_proc_id     = struct1%carbon14%owner_proc_id
    structOut%carbon14%owner_tile_path   = struct1%carbon14%owner_tile_path
    ALLOCATE(structOut%carbon14%ptr3d,   source=struct1%carbon14%ptr3d)
    structOut%nitrogen15%owner_model_id  = struct1%nitrogen15%owner_model_id
    structOut%nitrogen15%owner_proc_id   = struct1%nitrogen15%owner_proc_id
    structOut%nitrogen15%owner_tile_path = struct1%nitrogen15%owner_tile_path
    ALLOCATE(structOut%nitrogen15%ptr3d, source=struct1%nitrogen15%ptr3d)

    ! run calculation
    ! the differentiation between ON_CHUNK/ON_DOMAIN is done at the level of variables (i.e., the op-ov routines of the jsb_var)
    structOut%carbon          = struct1%carbon      * struct2%carbon
    structOut%nitrogen        = struct1%nitrogen    * struct2%nitrogen  
    structOut%phosphorus      = struct1%phosphorus  * struct2%phosphorus
    structOut%carbon13        = struct1%carbon13    * struct2%carbon13  
    structOut%carbon14        = struct1%carbon14    * struct2%carbon14  
    structOut%nitrogen15      = struct1%nitrogen15  * struct2%nitrogen15

    ! is this needed here too (not only with the jsb_var op-ov)
    structOut%carbon%ptr      => structOut%carbon%ptr3d(:,:,:)
    structOut%nitrogen%ptr    => structOut%nitrogen%ptr3d(:,:,:)
    structOut%phosphorus%ptr  => structOut%phosphorus%ptr3d(:,:,:)
    structOut%carbon13%ptr    => structOut%carbon13%ptr3d(:,:,:)
    structOut%carbon14%ptr    => structOut%carbon14%ptr3d(:,:,:)
    structOut%nitrogen15%ptr  => structOut%nitrogen15%ptr3d(:,:,:)
  
  END FUNCTION Multiply_sb_elem_p_with_sb_elem_p

  !--------------------------------------------------------------------------------------------------------------------------------
  ! / divide
  !--------------------------------------------------------------------------------------------------------------------------------
  !--------------------------------------------------------------------------------------------------------------------------------
  ! t_sb_elements_pool / scalar        | ON_CHUNK & ON_DOMAIN
  !--------------------------------------------------------------------------------------------------------------------------------
  FUNCTION Div_sb_elem_p_by_scalar(struct1, scalar) RESULT(structOut)
  
    CLASS(t_sb_elements_pool),      INTENT(in)  :: struct1
    REAL(wp),                       INTENT(in)  :: scalar
    TYPE(t_sb_elements_pool),       TARGET      :: structOut
  
    CHARACTER(len=*), PARAMETER :: routine = modname//':Div_sb_elem_p_by_scalar'
  
    ! init temporary pool structOut
    CALL structOut%Init(name='opov_help', shortname='opov_help', l_elements=.TRUE., element_unit=TRIM(struct1%element_unit))
    structOut%carbon%owner_model_id      = struct1%carbon%owner_model_id
    structOut%carbon%owner_proc_id       = struct1%carbon%owner_proc_id
    structOut%carbon%owner_tile_path     = struct1%carbon%owner_tile_path
    ALLOCATE(structOut%carbon%ptr3d,     source=struct1%carbon%ptr3d)
    structOut%nitrogen%owner_model_id    = struct1%nitrogen%owner_model_id
    structOut%nitrogen%owner_proc_id     = struct1%nitrogen%owner_proc_id
    structOut%nitrogen%owner_tile_path   = struct1%nitrogen%owner_tile_path
    ALLOCATE(structOut%nitrogen%ptr3d,   source=struct1%nitrogen%ptr3d)
    structOut%phosphorus%owner_model_id  = struct1%phosphorus%owner_model_id
    structOut%phosphorus%owner_proc_id   = struct1%phosphorus%owner_proc_id
    structOut%phosphorus%owner_tile_path = struct1%phosphorus%owner_tile_path
    ALLOCATE(structOut%phosphorus%ptr3d, source=struct1%phosphorus%ptr3d)
    structOut%carbon13%owner_model_id    = struct1%carbon13%owner_model_id
    structOut%carbon13%owner_proc_id     = struct1%carbon13%owner_proc_id
    structOut%carbon13%owner_tile_path   = struct1%carbon13%owner_tile_path
    ALLOCATE(structOut%carbon13%ptr3d,   source=struct1%carbon13%ptr3d)
    structOut%carbon14%owner_model_id    = struct1%carbon14%owner_model_id
    structOut%carbon14%owner_proc_id     = struct1%carbon14%owner_proc_id
    structOut%carbon14%owner_tile_path   = struct1%carbon14%owner_tile_path
    ALLOCATE(structOut%carbon14%ptr3d,   source=struct1%carbon14%ptr3d)
    structOut%nitrogen15%owner_model_id  = struct1%nitrogen15%owner_model_id
    structOut%nitrogen15%owner_proc_id   = struct1%nitrogen15%owner_proc_id
    structOut%nitrogen15%owner_tile_path = struct1%nitrogen15%owner_tile_path
    ALLOCATE(structOut%nitrogen15%ptr3d, source=struct1%nitrogen15%ptr3d)

    ! run calculation
    ! the differentiation between ON_CHUNK/ON_DOMAIN is done at the level of variables (i.e., the op-ov routines of the jsb_var)
    structOut%carbon          = struct1%carbon      / scalar
    structOut%nitrogen        = struct1%nitrogen    / scalar
    structOut%phosphorus      = struct1%phosphorus  / scalar
    structOut%carbon13        = struct1%carbon13    / scalar
    structOut%carbon14        = struct1%carbon14    / scalar
    structOut%nitrogen15      = struct1%nitrogen15  / scalar

    ! is this needed here too (not only with the jsb_var op-ov)
    structOut%carbon%ptr      => structOut%carbon%ptr3d(:,:,:)
    structOut%nitrogen%ptr    => structOut%nitrogen%ptr3d(:,:,:)
    structOut%phosphorus%ptr  => structOut%phosphorus%ptr3d(:,:,:)
    structOut%carbon13%ptr    => structOut%carbon13%ptr3d(:,:,:)
    structOut%carbon14%ptr    => structOut%carbon14%ptr3d(:,:,:)
    structOut%nitrogen15%ptr  => structOut%nitrogen15%ptr3d(:,:,:)
  
  END FUNCTION Div_sb_elem_p_by_scalar
  !--------------------------------------------------------------------------------------------------------------------------------
  ! t_sb_elements_pool / field2d(:,:)        | ON_CHUNK
  !--------------------------------------------------------------------------------------------------------------------------------
  FUNCTION Div_sb_elem_p_by_field2d(struct1, field2d) RESULT(structOut)
  
    CLASS(t_sb_elements_pool),      INTENT(in)  :: struct1
    REAL(wp),                       INTENT(in)  :: field2d(:,:)
    TYPE(t_sb_elements_pool),       TARGET      :: structOut
  
    CHARACTER(len=*), PARAMETER :: routine = modname//':Div_sb_elem_p_by_field2d'
  
    ! init temporary pool structOut
    CALL structOut%Init(name='opov_help', shortname='opov_help', l_elements=.TRUE., element_unit=TRIM(struct1%element_unit))
    structOut%carbon%owner_model_id      = struct1%carbon%owner_model_id
    structOut%carbon%owner_proc_id       = struct1%carbon%owner_proc_id
    structOut%carbon%owner_tile_path     = struct1%carbon%owner_tile_path
    ALLOCATE(structOut%carbon%ptr3d,     source=struct1%carbon%ptr3d)
    structOut%nitrogen%owner_model_id    = struct1%nitrogen%owner_model_id
    structOut%nitrogen%owner_proc_id     = struct1%nitrogen%owner_proc_id
    structOut%nitrogen%owner_tile_path   = struct1%nitrogen%owner_tile_path
    ALLOCATE(structOut%nitrogen%ptr3d,   source=struct1%nitrogen%ptr3d)
    structOut%phosphorus%owner_model_id  = struct1%phosphorus%owner_model_id
    structOut%phosphorus%owner_proc_id   = struct1%phosphorus%owner_proc_id
    structOut%phosphorus%owner_tile_path = struct1%phosphorus%owner_tile_path
    ALLOCATE(structOut%phosphorus%ptr3d, source=struct1%phosphorus%ptr3d)
    structOut%carbon13%owner_model_id    = struct1%carbon13%owner_model_id
    structOut%carbon13%owner_proc_id     = struct1%carbon13%owner_proc_id
    structOut%carbon13%owner_tile_path   = struct1%carbon13%owner_tile_path
    ALLOCATE(structOut%carbon13%ptr3d,   source=struct1%carbon13%ptr3d)
    structOut%carbon14%owner_model_id    = struct1%carbon14%owner_model_id
    structOut%carbon14%owner_proc_id     = struct1%carbon14%owner_proc_id
    structOut%carbon14%owner_tile_path   = struct1%carbon14%owner_tile_path
    ALLOCATE(structOut%carbon14%ptr3d,   source=struct1%carbon14%ptr3d)
    structOut%nitrogen15%owner_model_id  = struct1%nitrogen15%owner_model_id
    structOut%nitrogen15%owner_proc_id   = struct1%nitrogen15%owner_proc_id
    structOut%nitrogen15%owner_tile_path = struct1%nitrogen15%owner_tile_path
    ALLOCATE(structOut%nitrogen15%ptr3d, source=struct1%nitrogen15%ptr3d)

    ! run calculation
    ! the differentiation between ON_CHUNK/ON_DOMAIN is done at the level of variables (i.e., the op-ov routines of the jsb_var)
    structOut%carbon          = struct1%carbon      / field2d(:,:)
    structOut%nitrogen        = struct1%nitrogen    / field2d(:,:)
    structOut%phosphorus      = struct1%phosphorus  / field2d(:,:)
    structOut%carbon13        = struct1%carbon13    / field2d(:,:)
    structOut%carbon14        = struct1%carbon14    / field2d(:,:)
    structOut%nitrogen15      = struct1%nitrogen15  / field2d(:,:)

    ! is this needed here too (not only with the jsb_var op-ov)
    structOut%carbon%ptr      => structOut%carbon%ptr3d(:,:,:)
    structOut%nitrogen%ptr    => structOut%nitrogen%ptr3d(:,:,:)
    structOut%phosphorus%ptr  => structOut%phosphorus%ptr3d(:,:,:)
    structOut%carbon13%ptr    => structOut%carbon13%ptr3d(:,:,:)
    structOut%carbon14%ptr    => structOut%carbon14%ptr3d(:,:,:)
    structOut%nitrogen15%ptr  => structOut%nitrogen15%ptr3d(:,:,:)
  
  END FUNCTION Div_sb_elem_p_by_field2d
  !--------------------------------------------------------------------------------------------------------------------------------
  ! t_sb_elements_pool / field3d(:,:,:)        | ON_CHUNK & ON_DOMAIN | field3d(:,:,1) & field3d(:,:,:)
  !--------------------------------------------------------------------------------------------------------------------------------
  FUNCTION Div_sb_elem_p_by_field3d(struct1, field3d) RESULT(structOut)
  
    CLASS(t_sb_elements_pool),      INTENT(in)  :: struct1
    REAL(wp),                       INTENT(in)  :: field3d(:,:,:)
    TYPE(t_sb_elements_pool),       TARGET      :: structOut
  
    CHARACTER(len=*), PARAMETER :: routine = modname//':Div_sb_elem_p_by_field3d'
  
    ! init temporary pool structOut
    CALL structOut%Init(name='opov_help', shortname='opov_help', l_elements=.TRUE., element_unit=TRIM(struct1%element_unit))
    structOut%carbon%owner_model_id      = struct1%carbon%owner_model_id
    structOut%carbon%owner_proc_id       = struct1%carbon%owner_proc_id
    structOut%carbon%owner_tile_path     = struct1%carbon%owner_tile_path
    ALLOCATE(structOut%carbon%ptr3d,     source=struct1%carbon%ptr3d)
    structOut%nitrogen%owner_model_id    = struct1%nitrogen%owner_model_id
    structOut%nitrogen%owner_proc_id     = struct1%nitrogen%owner_proc_id
    structOut%nitrogen%owner_tile_path   = struct1%nitrogen%owner_tile_path
    ALLOCATE(structOut%nitrogen%ptr3d,   source=struct1%nitrogen%ptr3d)
    structOut%phosphorus%owner_model_id  = struct1%phosphorus%owner_model_id
    structOut%phosphorus%owner_proc_id   = struct1%phosphorus%owner_proc_id
    structOut%phosphorus%owner_tile_path = struct1%phosphorus%owner_tile_path
    ALLOCATE(structOut%phosphorus%ptr3d, source=struct1%phosphorus%ptr3d)
    structOut%carbon13%owner_model_id    = struct1%carbon13%owner_model_id
    structOut%carbon13%owner_proc_id     = struct1%carbon13%owner_proc_id
    structOut%carbon13%owner_tile_path   = struct1%carbon13%owner_tile_path
    ALLOCATE(structOut%carbon13%ptr3d,   source=struct1%carbon13%ptr3d)
    structOut%carbon14%owner_model_id    = struct1%carbon14%owner_model_id
    structOut%carbon14%owner_proc_id     = struct1%carbon14%owner_proc_id
    structOut%carbon14%owner_tile_path   = struct1%carbon14%owner_tile_path
    ALLOCATE(structOut%carbon14%ptr3d,   source=struct1%carbon14%ptr3d)
    structOut%nitrogen15%owner_model_id  = struct1%nitrogen15%owner_model_id
    structOut%nitrogen15%owner_proc_id   = struct1%nitrogen15%owner_proc_id
    structOut%nitrogen15%owner_tile_path = struct1%nitrogen15%owner_tile_path
    ALLOCATE(structOut%nitrogen15%ptr3d, source=struct1%nitrogen15%ptr3d)

    ! run calculation
    ! the differentiation between ON_CHUNK/ON_DOMAIN is done at the level of variables (i.e., the op-ov routines of the jsb_var)
    structOut%carbon          = struct1%carbon      / field3d(:,:,:)
    structOut%nitrogen        = struct1%nitrogen    / field3d(:,:,:)
    structOut%phosphorus      = struct1%phosphorus  / field3d(:,:,:)
    structOut%carbon13        = struct1%carbon13    / field3d(:,:,:)
    structOut%carbon14        = struct1%carbon14    / field3d(:,:,:)
    structOut%nitrogen15      = struct1%nitrogen15  / field3d(:,:,:)

    ! is this needed here too (not only with the jsb_var op-ov)
    structOut%carbon%ptr      => structOut%carbon%ptr3d(:,:,:)
    structOut%nitrogen%ptr    => structOut%nitrogen%ptr3d(:,:,:)
    structOut%phosphorus%ptr  => structOut%phosphorus%ptr3d(:,:,:)
    structOut%carbon13%ptr    => structOut%carbon13%ptr3d(:,:,:)
    structOut%carbon14%ptr    => structOut%carbon14%ptr3d(:,:,:)
    structOut%nitrogen15%ptr  => structOut%nitrogen15%ptr3d(:,:,:)
  
  END FUNCTION Div_sb_elem_p_by_field3d
  !--------------------------------------------------------------------------------------------------------------------------------
  ! t_sb_elements_pool / t_sb_elements_pool | ON_CHUNK & ON_DOMAIN
  !--------------------------------------------------------------------------------------------------------------------------------
  FUNCTION Div_sb_elem_p_by_sb_elem_p(struct1, struct2) RESULT(structOut)
  
    CLASS(t_sb_elements_pool),      INTENT(in)  :: struct1, struct2
    TYPE(t_sb_elements_pool),       TARGET      :: structOut
  
    CHARACTER(len=*), PARAMETER :: routine = modname//':Div_sb_elem_p_by_sb_elem_p'
  
    ! init temporary pool structOut
    CALL structOut%Init(name='opov_help', shortname='opov_help', l_elements=.TRUE., element_unit=TRIM(struct1%element_unit))
    structOut%carbon%owner_model_id      = struct1%carbon%owner_model_id
    structOut%carbon%owner_proc_id       = struct1%carbon%owner_proc_id
    structOut%carbon%owner_tile_path     = struct1%carbon%owner_tile_path
    ALLOCATE(structOut%carbon%ptr3d,     source=struct1%carbon%ptr3d)
    structOut%nitrogen%owner_model_id    = struct1%nitrogen%owner_model_id
    structOut%nitrogen%owner_proc_id     = struct1%nitrogen%owner_proc_id
    structOut%nitrogen%owner_tile_path   = struct1%nitrogen%owner_tile_path
    ALLOCATE(structOut%nitrogen%ptr3d,   source=struct1%nitrogen%ptr3d)
    structOut%phosphorus%owner_model_id  = struct1%phosphorus%owner_model_id
    structOut%phosphorus%owner_proc_id   = struct1%phosphorus%owner_proc_id
    structOut%phosphorus%owner_tile_path = struct1%phosphorus%owner_tile_path
    ALLOCATE(structOut%phosphorus%ptr3d, source=struct1%phosphorus%ptr3d)
    structOut%carbon13%owner_model_id    = struct1%carbon13%owner_model_id
    structOut%carbon13%owner_proc_id     = struct1%carbon13%owner_proc_id
    structOut%carbon13%owner_tile_path   = struct1%carbon13%owner_tile_path
    ALLOCATE(structOut%carbon13%ptr3d,   source=struct1%carbon13%ptr3d)
    structOut%carbon14%owner_model_id    = struct1%carbon14%owner_model_id
    structOut%carbon14%owner_proc_id     = struct1%carbon14%owner_proc_id
    structOut%carbon14%owner_tile_path   = struct1%carbon14%owner_tile_path
    ALLOCATE(structOut%carbon14%ptr3d,   source=struct1%carbon14%ptr3d)
    structOut%nitrogen15%owner_model_id  = struct1%nitrogen15%owner_model_id
    structOut%nitrogen15%owner_proc_id   = struct1%nitrogen15%owner_proc_id
    structOut%nitrogen15%owner_tile_path = struct1%nitrogen15%owner_tile_path
    ALLOCATE(structOut%nitrogen15%ptr3d, source=struct1%nitrogen15%ptr3d)

    ! run calculation
    ! the differentiation between ON_CHUNK/ON_DOMAIN is done at the level of variables (i.e., the op-ov routines of the jsb_var)
    structOut%carbon          = struct1%carbon      / struct2%carbon
    structOut%nitrogen        = struct1%nitrogen    / struct2%nitrogen  
    structOut%phosphorus      = struct1%phosphorus  / struct2%phosphorus
    structOut%carbon13        = struct1%carbon13    / struct2%carbon13  
    structOut%carbon14        = struct1%carbon14    / struct2%carbon14  
    structOut%nitrogen15      = struct1%nitrogen15  / struct2%nitrogen15

    ! is this needed here too (not only with the jsb_var op-ov)
    structOut%carbon%ptr      => structOut%carbon%ptr3d(:,:,:)
    structOut%nitrogen%ptr    => structOut%nitrogen%ptr3d(:,:,:)
    structOut%phosphorus%ptr  => structOut%phosphorus%ptr3d(:,:,:)
    structOut%carbon13%ptr    => structOut%carbon13%ptr3d(:,:,:)
    structOut%carbon14%ptr    => structOut%carbon14%ptr3d(:,:,:)
    structOut%nitrogen15%ptr  => structOut%nitrogen15%ptr3d(:,:,:)
  
  END FUNCTION Div_sb_elem_p_by_sb_elem_p

  !--------------------------------------------------------------------------------------------------------------------------------
  !> operator overload and assignment routines for components pools
  !!  TYPE t_sb_components_pool
  !! 
  !--------------------------------------------------------------------------------------------------------------------------------

  !--------------------------------------------------------------------------------------------------------------------------------
  ! finalize after operator overload
  ! t_sb_components_pool
  !--------------------------------------------------------------------------------------------------------------------------------
  SUBROUTINE finalize_sb_components_pool(this)

    TYPE(t_sb_components_pool) :: this

    CHARACTER(len=*), PARAMETER :: routine = modname//':finalize_sb_components_pool'

    !! CALL message(routine, 'Finalizing '//TRIM(this%name))

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
    
    IF (ASSOCIATED(this%dom))               DEALLOCATE(this%dom)
    IF (ASSOCIATED(this%dom_assoc))         DEALLOCATE(this%dom_assoc)
    IF (ASSOCIATED(this%soluable_litter))   DEALLOCATE(this%soluable_litter)
    IF (ASSOCIATED(this%polymeric_litter))  DEALLOCATE(this%polymeric_litter)
    IF (ASSOCIATED(this%woody_litter))      DEALLOCATE(this%woody_litter)
    IF (ASSOCIATED(this%fungi))             DEALLOCATE(this%fungi)
    IF (ASSOCIATED(this%mycorrhiza))        DEALLOCATE(this%mycorrhiza)
    IF (ASSOCIATED(this%microbial))         DEALLOCATE(this%microbial)
    IF (ASSOCIATED(this%residue))           DEALLOCATE(this%residue)
    IF (ASSOCIATED(this%residue_assoc))     DEALLOCATE(this%residue_assoc)
   
  END SUBROUTINE finalize_sb_components_pool

  !--------------------------------------------------------------------------------------------------------------------------------
  ! = assign
  !--------------------------------------------------------------------------------------------------------------------------------
  !--------------------------------------------------------------------------------------------------------------------------------
  ! sb_components_pool = scalar | ON_CHUNK & ON_DOMAIN
  !--------------------------------------------------------------------------------------------------------------------------------
  RECURSIVE SUBROUTINE Assign_scalar_to_sb_comp_p(lhs, rhs)

    CLASS(t_sb_components_pool),  INTENT(inout) :: lhs
    REAL(wp),                     INTENT(in)    :: rhs

    CHARACTER(len=*), PARAMETER :: routine = modname//':Assign_scalar_to_sb_comp_p'

    lhs%carbon      = rhs
    lhs%nitrogen    = rhs
    lhs%phosphorus  = rhs
    lhs%carbon13    = rhs
    lhs%carbon14    = rhs
    lhs%nitrogen15  = rhs
    
    lhs%dom               = rhs
    lhs%dom_assoc         = rhs
    lhs%soluable_litter   = rhs
    lhs%polymeric_litter  = rhs
    lhs%woody_litter      = rhs
    lhs%fungi             = rhs
    lhs%mycorrhiza        = rhs
    lhs%microbial         = rhs
    lhs%residue           = rhs
    lhs%residue_assoc     = rhs
    
  END SUBROUTINE Assign_scalar_to_sb_comp_p

  !--------------------------------------------------------------------------------------------------------------------------------
  ! sb_components_pool = sb_components_pool | ON_CHUNK & ON_DOMAIN
  !--------------------------------------------------------------------------------------------------------------------------------
  SUBROUTINE Assign_sb_comp_p_to_sb_comp_p(lhs, rhs)

    CLASS(t_sb_components_pool), INTENT(inout) :: lhs
    CLASS(t_sb_components_pool), INTENT(in)    :: rhs

    CHARACTER(len=*), PARAMETER :: routine = modname//':Assign_sb_comp_p_to_sb_comp_p'

    lhs%carbon      = rhs%carbon    
    lhs%nitrogen    = rhs%nitrogen  
    lhs%phosphorus  = rhs%phosphorus
    lhs%carbon13    = rhs%carbon13  
    lhs%carbon14    = rhs%carbon14  
    lhs%nitrogen15  = rhs%nitrogen15
    
    lhs%dom               = rhs%dom             
    lhs%dom_assoc         = rhs%dom_assoc       
    lhs%soluable_litter   = rhs%soluable_litter 
    lhs%polymeric_litter  = rhs%polymeric_litter
    lhs%woody_litter      = rhs%woody_litter    
    lhs%fungi             = rhs%fungi           
    lhs%mycorrhiza        = rhs%mycorrhiza      
    lhs%microbial         = rhs%microbial       
    lhs%residue           = rhs%residue         
    lhs%residue_assoc     = rhs%residue_assoc   

  END SUBROUTINE Assign_sb_comp_p_to_sb_comp_p





  !--------------------------------------------------------------------------------------------------------------------------------
  ! + addition
  !--------------------------------------------------------------------------------------------------------------------------------
  !--------------------------------------------------------------------------------------------------------------------------------
  ! sb_components_pool + scalar        | ON_CHUNK & ON_DOMAIN
  !--------------------------------------------------------------------------------------------------------------------------------
  FUNCTION Add_sb_comp_p_and_scalar(struct1, scalar) RESULT(structOut)
  
    CLASS(t_sb_components_pool),      INTENT(in)  :: struct1
    REAL(wp),                         INTENT(in)  :: scalar
    TYPE(t_sb_components_pool),       TARGET      :: structOut
  
    CHARACTER(len=*), PARAMETER :: routine = modname//':Add_sb_comp_p_and_scalar'

    !> 1.0 init temporary pool structOut
    !! init components pool with element variables
    CALL structOut%Init_comp_pool_for_op_ov(name='comp_op_ov_help', shortname='c_op_ov_help', &
      & unit=TRIM(struct1%element_unit), pool_suffix='sb_p') ! l_elements=.TRUE. by default
    structOut%carbon%owner_model_id      = struct1%carbon%owner_model_id
    structOut%carbon%owner_proc_id       = struct1%carbon%owner_proc_id
    structOut%carbon%owner_tile_path     = struct1%carbon%owner_tile_path
    ALLOCATE(structOut%carbon%ptr3d,     source=struct1%carbon%ptr3d)
    structOut%nitrogen%owner_model_id    = struct1%nitrogen%owner_model_id
    structOut%nitrogen%owner_proc_id     = struct1%nitrogen%owner_proc_id
    structOut%nitrogen%owner_tile_path   = struct1%nitrogen%owner_tile_path
    ALLOCATE(structOut%nitrogen%ptr3d,   source=struct1%nitrogen%ptr3d)
    structOut%phosphorus%owner_model_id  = struct1%phosphorus%owner_model_id
    structOut%phosphorus%owner_proc_id   = struct1%phosphorus%owner_proc_id
    structOut%phosphorus%owner_tile_path = struct1%phosphorus%owner_tile_path
    ALLOCATE(structOut%phosphorus%ptr3d, source=struct1%phosphorus%ptr3d)
    structOut%carbon13%owner_model_id    = struct1%carbon13%owner_model_id
    structOut%carbon13%owner_proc_id     = struct1%carbon13%owner_proc_id
    structOut%carbon13%owner_tile_path   = struct1%carbon13%owner_tile_path
    ALLOCATE(structOut%carbon13%ptr3d,   source=struct1%carbon13%ptr3d)
    structOut%carbon14%owner_model_id    = struct1%carbon14%owner_model_id
    structOut%carbon14%owner_proc_id     = struct1%carbon14%owner_proc_id
    structOut%carbon14%owner_tile_path   = struct1%carbon14%owner_tile_path
    ALLOCATE(structOut%carbon14%ptr3d,   source=struct1%carbon14%ptr3d)
    structOut%nitrogen15%owner_model_id  = struct1%nitrogen15%owner_model_id
    structOut%nitrogen15%owner_proc_id   = struct1%nitrogen15%owner_proc_id
    structOut%nitrogen15%owner_tile_path = struct1%nitrogen15%owner_tile_path
    ALLOCATE(structOut%nitrogen15%ptr3d, source=struct1%nitrogen15%ptr3d)
    !! init element pools
    CALL structOut%dom%Init_elem_var_elem_p(struct1%dom)
    CALL structOut%dom_assoc%Init_elem_var_elem_p(struct1%dom_assoc)
    CALL structOut%soluable_litter%Init_elem_var_elem_p(struct1%soluable_litter)
    CALL structOut%polymeric_litter%Init_elem_var_elem_p(struct1%polymeric_litter)
    CALL structOut%woody_litter%Init_elem_var_elem_p(struct1%woody_litter)
    CALL structOut%fungi%Init_elem_var_elem_p(struct1%fungi)
    CALL structOut%mycorrhiza%Init_elem_var_elem_p(struct1%mycorrhiza)
    CALL structOut%microbial%Init_elem_var_elem_p(struct1%microbial)
    CALL structOut%residue%Init_elem_var_elem_p(struct1%residue)
    CALL structOut%residue_assoc%Init_elem_var_elem_p(struct1%residue_assoc)

    !> 2.0 run calculation
    !! the differentiation between ON_CHUNK/ON_DOMAIN is done at the level of variables (i.e., the op_ov routines of the jsb_var)
    !> 2.1 element variables of the components pool itself
    structOut%carbon          = struct1%carbon      + scalar
    structOut%nitrogen        = struct1%nitrogen    + scalar
    structOut%phosphorus      = struct1%phosphorus  + scalar
    structOut%carbon13        = struct1%carbon13    + scalar
    structOut%carbon14        = struct1%carbon14    + scalar
    structOut%nitrogen15      = struct1%nitrogen15  + scalar
    ! adding the pointer is this needed here too (not only with the jsb_var op_ov)
    structOut%carbon%ptr      => structOut%carbon%ptr3d(:,:,:)
    structOut%nitrogen%ptr    => structOut%nitrogen%ptr3d(:,:,:)
    structOut%phosphorus%ptr  => structOut%phosphorus%ptr3d(:,:,:)
    structOut%carbon13%ptr    => structOut%carbon13%ptr3d(:,:,:)
    structOut%carbon14%ptr    => structOut%carbon14%ptr3d(:,:,:)
    structOut%nitrogen15%ptr  => structOut%nitrogen15%ptr3d(:,:,:)
    !> 2.2 element pools
    structOut%dom               = struct1%dom               + scalar
    structOut%dom_assoc         = struct1%dom_assoc         + scalar
    structOut%soluable_litter   = struct1%soluable_litter   + scalar
    structOut%polymeric_litter  = struct1%polymeric_litter  + scalar
    structOut%woody_litter      = struct1%woody_litter      + scalar
    structOut%fungi             = struct1%fungi             + scalar
    structOut%mycorrhiza        = struct1%mycorrhiza        + scalar
    structOut%microbial         = struct1%microbial         + scalar
    structOut%residue           = struct1%residue           + scalar
    structOut%residue_assoc     = struct1%residue_assoc     + scalar

  END FUNCTION Add_sb_comp_p_and_scalar




END MODULE mo_sb_memory_class
