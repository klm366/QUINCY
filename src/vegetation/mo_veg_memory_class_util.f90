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
!! @brief Contains the major pool code used by the memory class for the vegetation processes
!!
!! @author
!!  Jan Engel, MPI-BGC Jena
!!
!!
MODULE mo_veg_memory_class_util

!! #################################################################################################################################
!! code for QUINCY model without JSBACH4 - ends with (searchterm) "quincy-no-jsb_end"
!! #################################################################################################################################
#ifdef __QUINCY_NO_JSB__

  USE mo_kind,                   ONLY: wp
  USE mo_jsb_var_class,          ONLY: t_jsb_var_real2d, &
                                       t_jsb_struct, t_jsb_n_struct      ! quincy specific

  IMPLICIT NONE

  ! all is public except the modename; code is used by the mo_veg_memory_class only
  PUBLIC
  PRIVATE :: modname

  !-----------------------------------------------------------------------------------------------------
  !> Type definition for helper elements with ALLOCATABLE variables for elements in veg memory
  !! 
  !! 
  !-----------------------------------------------------------------------------------------------------
  TYPE t_veg_element_a
    REAL(wp), ALLOCATABLE :: carbon(:)      ! dimension of variable in memory is var(ics:ice,iblk)
    REAL(wp), ALLOCATABLE :: nitrogen(:)
    REAL(wp), ALLOCATABLE :: phosphorus(:)
    REAL(wp), ALLOCATABLE :: carbon13(:)
    REAL(wp), ALLOCATABLE :: carbon14(:)
    REAL(wp), ALLOCATABLE :: nitrogen15(:)
  CONTAINS
    PRIVATE
    !! operator overload
    ! +
    PROCEDURE               :: Add_element_a                   ! Add two elements
    PROCEDURE               :: Add_element_a_and_scalar        ! Add scalar to element
    PROCEDURE               :: Add_element_a_and_vector       ! Add REAL 1D vector ( var(ics:ice) ) to element
    GENERIC, PUBLIC         :: OPERATOR(+) => Add_element_a, Add_element_a_and_scalar, Add_element_a_and_vector
    !! operator overload
    ! -
    PROCEDURE               :: Subtract_element_a                   ! Subtract one element from another element
    PROCEDURE               :: Subtract_element_a_and_scalar        ! Subtract scalar from element
    PROCEDURE               :: Subtract_element_a_and_vector        ! Subtract vector ( var(ics:ice) ) from element
    GENERIC, PUBLIC         :: OPERATOR(-) => Subtract_element_a, Subtract_element_a_and_scalar, Subtract_element_a_and_vector
    !! operator overload
    ! *
    PROCEDURE               :: Multiply_element_a                   ! Multiply two elements
    PROCEDURE               :: Multiply_element_a_and_scalar        ! Multiply element with 
    PROCEDURE, PASS(this2)  :: Multiply_scalar_and_element_a        ! Multiply scalar with element
    PROCEDURE               :: Multiply_element_a_and_vector        ! Multiply element with REAL 1D vector ( var(ics:ice) )
    PROCEDURE, PASS(this2)  :: Multiply_vector_and_element_a        ! Multiply REAL 1D vector ( var(ics:ice) ) with element
    GENERIC, PUBLIC         :: OPERATOR(*) => Multiply_element_a, Multiply_element_a_and_scalar, Multiply_scalar_and_element_a, &
                                              Multiply_element_a_and_vector, Multiply_vector_and_element_a
  END TYPE t_veg_element_a

  !-----------------------------------------------------------------------------------------------------
  !> Type definition for helper pools with ALLOCATABLE variables for pools in veg memory
  !! 
  !! 
  !-----------------------------------------------------------------------------------------------------
  TYPE t_veg_pool_a
    TYPE(t_veg_element_a)  :: leaf             !< elements of leaf
    TYPE(t_veg_element_a)  :: fine_root        !< elements of fine_root
    TYPE(t_veg_element_a)  :: coarse_root      !< elements of coarse_root
    TYPE(t_veg_element_a)  :: sap_wood         !< elements of sap_wood
    TYPE(t_veg_element_a)  :: heart_wood       !< elements of heart_wood
    TYPE(t_veg_element_a)  :: labile           !< elements of labile
    TYPE(t_veg_element_a)  :: reserve          !< elements of reserve
    TYPE(t_veg_element_a)  :: fruit            !< elements of fruit
    TYPE(t_veg_element_a)  :: seed_bed         !< elements of seed_bed
  CONTAINS
    PRIVATE
    !! operator overload
    ! +
    PROCEDURE               :: Add_pool_a                   ! Add two pools
    PROCEDURE               :: Add_pool_a_and_scalar        ! Add scalar to pool
    PROCEDURE               :: Add_pool_a_and_vector        ! Add REAL 1D vector ( var(ics:ice) ) to pool
    GENERIC, PUBLIC         :: OPERATOR(+) => Add_pool_a, Add_pool_a_and_scalar, Add_pool_a_and_vector
    !! operator overload
    ! -
    PROCEDURE               :: Subtract_pool_a                   ! Subtract one pool from another pool
    PROCEDURE               :: Subtract_pool_a_and_scalar        ! Subtract scalar from pool
    PROCEDURE               :: Subtract_pool_a_and_vector        ! Subtract vector ( var(ics:ice) ) from pool
    GENERIC, PUBLIC         :: OPERATOR(-) => Subtract_pool_a, Subtract_pool_a_and_scalar, Subtract_pool_a_and_vector
    !! operator overload
    ! *
    PROCEDURE               :: Multiply_pool_a                   ! Multiply two pools
    PROCEDURE               :: Multiply_pool_a_and_scalar        ! Multiply pool with scalar
    PROCEDURE, PASS(this2)  :: Multiply_scalar_and_pool_a        ! Multiply scalar with pool
    PROCEDURE               :: Multiply_pool_a_and_vector        ! Multiply pool with REAL 1D vector ( var(ics:ice) )
    PROCEDURE, PASS(this2)  :: Multiply_vector_and_pool_a        ! Multiply REAL 1D vector ( var(ics:ice) ) with pool
    GENERIC, PUBLIC         :: OPERATOR(*) => Multiply_pool_a, Multiply_pool_a_and_scalar, Multiply_scalar_and_pool_a, &
                                              Multiply_pool_a_and_vector, Multiply_vector_and_pool_a
  END TYPE t_veg_pool_a

  !-----------------------------------------------------------------------------------------------------
  !> Type definition for elements in veg pool
  !! 
  !! elements in veg would be 2D (while they would be 3D in soil)
  !! 
  !-----------------------------------------------------------------------------------------------------
  TYPE, EXTENDS(t_jsb_struct) :: t_quincy_struct_veg_element
    TYPE(t_jsb_var_real2d)    :: carbon       !< total carbon concentration of pool in mol / m2
    TYPE(t_jsb_var_real2d)    :: nitrogen     !< totol nitrogen concentration of pool in mol / m2
    TYPE(t_jsb_var_real2d)    :: phosphorus   !< phosphorus concentration of pool in mol / m2
    TYPE(t_jsb_var_real2d)    :: carbon13     !< 13carbon concentration of pool in mol / m2
    TYPE(t_jsb_var_real2d)    :: carbon14     !< 14carbon concentration of pool in pmol / m2
    TYPE(t_jsb_var_real2d)    :: nitrogen15   !< 15nitrogen concentration of pool in mol / m2
  CONTAINS
    PROCEDURE :: Init       => Init_t_quincy_struct_veg_element
    PROCEDURE :: Add_vars   => Add_element_variables_to_memory
  END TYPE t_quincy_struct_veg_element

  !-----------------------------------------------------------------------------------------------------
  !> Type definition for veg pool with nested TYPEs
  !!
  !! used for: pool, and fluxes: growth and litterfall
  !!  
  !-----------------------------------------------------------------------------------------------------
  TYPE, EXTENDS(t_jsb_n_struct) :: t_quincy_struct_veg_pool
    TYPE(t_quincy_struct_veg_element)  :: leaf             !< elements of leaf
    TYPE(t_quincy_struct_veg_element)  :: fine_root        !< elements of fine_root
    TYPE(t_quincy_struct_veg_element)  :: coarse_root      !< elements of coarse_root
    TYPE(t_quincy_struct_veg_element)  :: sap_wood         !< elements of sap_wood
    TYPE(t_quincy_struct_veg_element)  :: heart_wood       !< elements of heart_wood
    TYPE(t_quincy_struct_veg_element)  :: labile           !< elements of labile
    TYPE(t_quincy_struct_veg_element)  :: reserve          !< elements of reserve
    TYPE(t_quincy_struct_veg_element)  :: fruit            !< elements of fruit
    TYPE(t_quincy_struct_veg_element)  :: seed_bed         !< elements of seed_bed
  CONTAINS
    PROCEDURE :: Init           => Init_t_quincy_struct_veg_pool
    PROCEDURE :: Add_elements   => Add_pool_elements_to_memory
  END TYPE t_quincy_struct_veg_pool

  ! definition of veg_ pool and element structures with ALLOCATABLE variables
  ! allocation of variables done with Init_veg_pools_elements_a during model init in mo_jsb_model_init:quincy_init
  TYPE(t_veg_pool_a)     :: veg_pool           , &
                            veg_growth         , &
                            veg_litterfall     , &
                            veg_falloc         , &
                            veg_total_flux
  TYPE(t_veg_element_a)  :: veg_exudation      , &
                            veg_establishment  , &
                            veg_reserve_use    , &
                            veg_total_biomass

  CHARACTER(len=*), PARAMETER :: modname = 'mo_veg_memory_class_util'

CONTAINS

  !-----------------------------------------------------------------------------------------------------!
  ! t_veg_element_a operator overload                                                                   !
  !-----------------------------------------------------------------------------------------------------!
  !-----------------------------------------------------------------------------------------------------
  !> + (t_veg_element_a + t_veg_element_a)
  !!
  !-----------------------------------------------------------------------------------------------------
  PURE FUNCTION Add_element_a(this,element2) RESULT(element3)
    IMPLICIT NONE
    ! 0.1 InOut
    CLASS(t_veg_element_a), INTENT(in)  :: this, element2
    TYPE(t_veg_element_a)               :: element3
    ! 1.0 go
    ALLOCATE(element3%carbon,     source=this%carbon)       ! note: Function results are automatically deallocated after the result is "used" in the outer scope.
    ALLOCATE(element3%nitrogen,   source=this%nitrogen)
    ALLOCATE(element3%phosphorus, source=this%phosphorus)
    ALLOCATE(element3%carbon13,   source=this%carbon13)
    ALLOCATE(element3%carbon14,   source=this%carbon14)
    ALLOCATE(element3%nitrogen15, source=this%nitrogen15)
    element3%carbon(:)       = this%carbon(:)      + element2%carbon(:)
    element3%nitrogen(:)     = this%nitrogen(:)    + element2%nitrogen(:)
    element3%phosphorus(:)   = this%phosphorus(:)  + element2%phosphorus(:)
    element3%carbon13(:)     = this%carbon13(:)    + element2%carbon13(:)
    element3%carbon14(:)     = this%carbon14(:)    + element2%carbon14(:)
    element3%nitrogen15(:)   = this%nitrogen15(:)  + element2%nitrogen15(:)
  END FUNCTION Add_element_a
  !-----------------------------------------------------------------------------------------------------
  !> + (t_veg_element_a + scalar)
  !!
  !-----------------------------------------------------------------------------------------------------
  PURE FUNCTION Add_element_a_and_scalar(this,scalar) RESULT(element3)
    IMPLICIT NONE
    ! 0.1 InOut
    CLASS(t_veg_element_a), INTENT(in)  :: this
    REAL(wp),               INTENT(in)  :: scalar
    TYPE(t_veg_element_a)               :: element3
    ! 1.0 go
    ALLOCATE(element3%carbon,     source=this%carbon)
    ALLOCATE(element3%nitrogen,   source=this%nitrogen)
    ALLOCATE(element3%phosphorus, source=this%phosphorus)
    ALLOCATE(element3%carbon13,   source=this%carbon13)
    ALLOCATE(element3%carbon14,   source=this%carbon14)
    ALLOCATE(element3%nitrogen15, source=this%nitrogen15)
    element3%carbon(:)       = this%carbon(:)      + scalar
    element3%nitrogen(:)     = this%nitrogen(:)    + scalar
    element3%phosphorus(:)   = this%phosphorus(:)  + scalar
    element3%carbon13(:)     = this%carbon13(:)    + scalar
    element3%carbon14(:)     = this%carbon14(:)    + scalar
    element3%nitrogen15(:)   = this%nitrogen15(:)  + scalar
  END FUNCTION Add_element_a_and_scalar
  !-----------------------------------------------------------------------------------------------------
  !> + (t_veg_element_a + 1D vector)
  !!
  !-----------------------------------------------------------------------------------------------------
  PURE FUNCTION Add_element_a_and_vector(this,vector) RESULT(element3)
    IMPLICIT NONE
    ! 0.1 InOut
    CLASS(t_veg_element_a), INTENT(in)  :: this
    REAL(wp),               INTENT(in)  :: vector(:)
    TYPE(t_veg_element_a)               :: element3
    ! 1.0 go
    ALLOCATE(element3%carbon,     source=this%carbon)
    ALLOCATE(element3%nitrogen,   source=this%nitrogen)
    ALLOCATE(element3%phosphorus, source=this%phosphorus)
    ALLOCATE(element3%carbon13,   source=this%carbon13)
    ALLOCATE(element3%carbon14,   source=this%carbon14)
    ALLOCATE(element3%nitrogen15, source=this%nitrogen15)
    element3%carbon(:)       = this%carbon(:)      + vector(:)
    element3%nitrogen(:)     = this%nitrogen(:)    + vector(:)
    element3%phosphorus(:)   = this%phosphorus(:)  + vector(:)
    element3%carbon13(:)     = this%carbon13(:)    + vector(:)
    element3%carbon14(:)     = this%carbon14(:)    + vector(:)
    element3%nitrogen15(:)   = this%nitrogen15(:)  + vector(:)
  END FUNCTION Add_element_a_and_vector
  !-----------------------------------------------------------------------------------------------------
  !> - (t_veg_element_a - t_veg_element_a)
  !!
  !-----------------------------------------------------------------------------------------------------
  PURE FUNCTION Subtract_element_a(this,element2) RESULT(element3)
    IMPLICIT NONE
    ! 0.1 InOut
    CLASS(t_veg_element_a), INTENT(in)  :: this, element2
    TYPE(t_veg_element_a)               :: element3
    ! 1.0 go
    ALLOCATE(element3%carbon,     source=this%carbon)
    ALLOCATE(element3%nitrogen,   source=this%nitrogen)
    ALLOCATE(element3%phosphorus, source=this%phosphorus)
    ALLOCATE(element3%carbon13,   source=this%carbon13)
    ALLOCATE(element3%carbon14,   source=this%carbon14)
    ALLOCATE(element3%nitrogen15, source=this%nitrogen15)
    element3%carbon(:)       = this%carbon(:)      - element2%carbon(:)
    element3%nitrogen(:)     = this%nitrogen(:)    - element2%nitrogen(:)
    element3%phosphorus(:)   = this%phosphorus(:)  - element2%phosphorus(:)
    element3%carbon13(:)     = this%carbon13(:)    - element2%carbon13(:)
    element3%carbon14(:)     = this%carbon14(:)    - element2%carbon14(:)
    element3%nitrogen15(:)   = this%nitrogen15(:)  - element2%nitrogen15(:)
  END FUNCTION Subtract_element_a
  !-----------------------------------------------------------------------------------------------------
  !> - (t_veg_element_a - scalar)
  !!
  !-----------------------------------------------------------------------------------------------------
  PURE FUNCTION Subtract_element_a_and_scalar(this,scalar) RESULT(element3)
    IMPLICIT NONE
    ! 0.1 InOut
    CLASS(t_veg_element_a), INTENT(in)  :: this
    REAL(wp),               INTENT(in)  :: scalar
    TYPE(t_veg_element_a)               :: element3
    ! 1.0 go
    ALLOCATE(element3%carbon,     source=this%carbon)
    ALLOCATE(element3%nitrogen,   source=this%nitrogen)
    ALLOCATE(element3%phosphorus, source=this%phosphorus)
    ALLOCATE(element3%carbon13,   source=this%carbon13)
    ALLOCATE(element3%carbon14,   source=this%carbon14)
    ALLOCATE(element3%nitrogen15, source=this%nitrogen15)
    element3%carbon(:)       = this%carbon(:)      - scalar
    element3%nitrogen(:)     = this%nitrogen(:)    - scalar
    element3%phosphorus(:)   = this%phosphorus(:)  - scalar
    element3%carbon13(:)     = this%carbon13(:)    - scalar
    element3%carbon14(:)     = this%carbon14(:)    - scalar
    element3%nitrogen15(:)   = this%nitrogen15(:)  - scalar
  END FUNCTION Subtract_element_a_and_scalar
  !-----------------------------------------------------------------------------------------------------
  !> - (t_veg_element_a - 1D vector)
  !!
  !-----------------------------------------------------------------------------------------------------
  PURE FUNCTION Subtract_element_a_and_vector(this,vector) RESULT(element3)
    IMPLICIT NONE
    ! 0.1 InOut
    CLASS(t_veg_element_a), INTENT(in)  :: this
    REAL(wp),               INTENT(in)  :: vector(:)
    TYPE(t_veg_element_a)               :: element3
    ! 1.0 go
    ALLOCATE(element3%carbon,     source=this%carbon)
    ALLOCATE(element3%nitrogen,   source=this%nitrogen)
    ALLOCATE(element3%phosphorus, source=this%phosphorus)
    ALLOCATE(element3%carbon13,   source=this%carbon13)
    ALLOCATE(element3%carbon14,   source=this%carbon14)
    ALLOCATE(element3%nitrogen15, source=this%nitrogen15)
    element3%carbon(:)       = this%carbon(:)      - vector(:)
    element3%nitrogen(:)     = this%nitrogen(:)    - vector(:)
    element3%phosphorus(:)   = this%phosphorus(:)  - vector(:)
    element3%carbon13(:)     = this%carbon13(:)    - vector(:)
    element3%carbon14(:)     = this%carbon14(:)    - vector(:)
    element3%nitrogen15(:)   = this%nitrogen15(:)  - vector(:)
  END FUNCTION Subtract_element_a_and_vector
  !-----------------------------------------------------------------------------------------------------
  !> * (t_veg_element_a * t_veg_element_a)
  !!
  !-----------------------------------------------------------------------------------------------------
  PURE FUNCTION Multiply_element_a(this,element2) RESULT(element3)
    IMPLICIT NONE
    ! 0.1 InOut
    CLASS(t_veg_element_a), INTENT(in)  :: this, element2
    TYPE(t_veg_element_a)               :: element3
    ! 1.0 go
    ALLOCATE(element3%carbon,     source=this%carbon)
    ALLOCATE(element3%nitrogen,   source=this%nitrogen)
    ALLOCATE(element3%phosphorus, source=this%phosphorus)
    ALLOCATE(element3%carbon13,   source=this%carbon13)
    ALLOCATE(element3%carbon14,   source=this%carbon14)
    ALLOCATE(element3%nitrogen15, source=this%nitrogen15)
    element3%carbon(:)       = this%carbon(:)      * element2%carbon(:)
    element3%nitrogen(:)     = this%nitrogen(:)    * element2%nitrogen(:)
    element3%phosphorus(:)   = this%phosphorus(:)  * element2%phosphorus(:)
    element3%carbon13(:)     = this%carbon13(:)    * element2%carbon13(:)
    element3%carbon14(:)     = this%carbon14(:)    * element2%carbon14(:)
    element3%nitrogen15(:)   = this%nitrogen15(:)  * element2%nitrogen15(:)
  END FUNCTION Multiply_element_a
  !-----------------------------------------------------------------------------------------------------
  !> * (t_veg_element_a * scalar)
  !!
  !-----------------------------------------------------------------------------------------------------
  PURE FUNCTION Multiply_element_a_and_scalar(this,scalar) RESULT(element3)
    IMPLICIT NONE
    ! 0.1 InOut
    CLASS(t_veg_element_a), INTENT(in)  :: this
    REAL(wp),               INTENT(in)  :: scalar
    TYPE(t_veg_element_a)               :: element3
    ! 1.0 go
    ALLOCATE(element3%carbon,     source=this%carbon)
    ALLOCATE(element3%nitrogen,   source=this%nitrogen)
    ALLOCATE(element3%phosphorus, source=this%phosphorus)
    ALLOCATE(element3%carbon13,   source=this%carbon13)
    ALLOCATE(element3%carbon14,   source=this%carbon14)
    ALLOCATE(element3%nitrogen15, source=this%nitrogen15)
    element3%carbon(:)       = this%carbon(:)      * scalar
    element3%nitrogen(:)     = this%nitrogen(:)    * scalar
    element3%phosphorus(:)   = this%phosphorus(:)  * scalar
    element3%carbon13(:)     = this%carbon13(:)    * scalar
    element3%carbon14(:)     = this%carbon14(:)    * scalar
    element3%nitrogen15(:)   = this%nitrogen15(:)  * scalar
  END FUNCTION Multiply_element_a_and_scalar
  !-----------------------------------------------------------------------------------------------------
  !> * (scalar * t_veg_element_a)
  !!
  !-----------------------------------------------------------------------------------------------------
  PURE FUNCTION Multiply_scalar_and_element_a(scalar,this2) RESULT(element3)
    IMPLICIT NONE
    ! 0.1 InOut
    REAL(wp),               INTENT(in)  :: scalar
    CLASS(t_veg_element_a), INTENT(in)  :: this2
    TYPE(t_veg_element_a)               :: element3
    ! 1.0 go
    ALLOCATE(element3%carbon,     source=this2%carbon)
    ALLOCATE(element3%nitrogen,   source=this2%nitrogen)
    ALLOCATE(element3%phosphorus, source=this2%phosphorus)
    ALLOCATE(element3%carbon13,   source=this2%carbon13)
    ALLOCATE(element3%carbon14,   source=this2%carbon14)
    ALLOCATE(element3%nitrogen15, source=this2%nitrogen15)
    element3%carbon(:)       = scalar * this2%carbon(:)
    element3%nitrogen(:)     = scalar * this2%nitrogen(:)
    element3%phosphorus(:)   = scalar * this2%phosphorus(:)
    element3%carbon13(:)     = scalar * this2%carbon13(:)
    element3%carbon14(:)     = scalar * this2%carbon14(:)
    element3%nitrogen15(:)   = scalar * this2%nitrogen15(:)
  END FUNCTION Multiply_scalar_and_element_a
  !-----------------------------------------------------------------------------------------------------
  !> * (t_veg_element_a * 1D vector)
  !!
  !-----------------------------------------------------------------------------------------------------
  PURE FUNCTION Multiply_element_a_and_vector(this,vector) RESULT(element3)
    IMPLICIT NONE
    ! 0.1 InOut
    CLASS(t_veg_element_a), INTENT(in)  :: this
    REAL(wp),               INTENT(in)  :: vector(:)
    TYPE(t_veg_element_a)               :: element3
    ! 1.0 go
    ALLOCATE(element3%carbon,     source=this%carbon)
    ALLOCATE(element3%nitrogen,   source=this%nitrogen)
    ALLOCATE(element3%phosphorus, source=this%phosphorus)
    ALLOCATE(element3%carbon13,   source=this%carbon13)
    ALLOCATE(element3%carbon14,   source=this%carbon14)
    ALLOCATE(element3%nitrogen15, source=this%nitrogen15)
    element3%carbon(:)       = this%carbon(:)      * vector(:)
    element3%nitrogen(:)     = this%nitrogen(:)    * vector(:)
    element3%phosphorus(:)   = this%phosphorus(:)  * vector(:)
    element3%carbon13(:)     = this%carbon13(:)    * vector(:)
    element3%carbon14(:)     = this%carbon14(:)    * vector(:)
    element3%nitrogen15(:)   = this%nitrogen15(:)  * vector(:)
  END FUNCTION Multiply_element_a_and_vector
  !-----------------------------------------------------------------------------------------------------
  !> * (vector * t_veg_element_a)
  !!
  !-----------------------------------------------------------------------------------------------------
  PURE FUNCTION Multiply_vector_and_element_a(vector,this2) RESULT(element3)
    IMPLICIT NONE
    ! 0.1 InOut
    REAL(wp),               INTENT(in)  :: vector(:)
    CLASS(t_veg_element_a), INTENT(in)  :: this2
    TYPE(t_veg_element_a)               :: element3
    ! 1.0 go
    ALLOCATE(element3%carbon,     source=this2%carbon)
    ALLOCATE(element3%nitrogen,   source=this2%nitrogen)
    ALLOCATE(element3%phosphorus, source=this2%phosphorus)
    ALLOCATE(element3%carbon13,   source=this2%carbon13)
    ALLOCATE(element3%carbon14,   source=this2%carbon14)
    ALLOCATE(element3%nitrogen15, source=this2%nitrogen15)
    element3%carbon(:)       = vector(:) * this2%carbon(:)
    element3%nitrogen(:)     = vector(:) * this2%nitrogen(:)
    element3%phosphorus(:)   = vector(:) * this2%phosphorus(:)
    element3%carbon13(:)     = vector(:) * this2%carbon13(:)
    element3%carbon14(:)     = vector(:) * this2%carbon14(:)
    element3%nitrogen15(:)   = vector(:) * this2%nitrogen15(:)
  END FUNCTION Multiply_vector_and_element_a

  !-----------------------------------------------------------------------------------------------------!
  ! t_veg_pool_a operator overload                                                                      !
  !-----------------------------------------------------------------------------------------------------!
  !-----------------------------------------------------------------------------------------------------
  !> + (t_veg_pool_a + t_veg_pool_a)
  !!
  !-----------------------------------------------------------------------------------------------------
  PURE FUNCTION Add_pool_a(this,pool2) RESULT(pool3)
    IMPLICIT NONE
    ! 0.1 InOut
    CLASS(t_veg_pool_a), INTENT(in) :: this, pool2
    TYPE(t_veg_pool_a)              :: pool3
    ! 1.0 go
    pool3%leaf          = this%leaf         + pool2%leaf
    pool3%fine_root     = this%fine_root    + pool2%fine_root
    pool3%coarse_root   = this%coarse_root  + pool2%coarse_root
    pool3%sap_wood      = this%sap_wood     + pool2%sap_wood
    pool3%heart_wood    = this%heart_wood   + pool2%heart_wood
    pool3%labile        = this%labile       + pool2%labile
    pool3%reserve       = this%reserve      + pool2%reserve
    pool3%fruit         = this%fruit        + pool2%fruit
    pool3%seed_bed      = this%seed_bed     + pool2%seed_bed
  END FUNCTION Add_pool_a
  !-----------------------------------------------------------------------------------------------------
  !> + (t_veg_pool_a + scalar)
  !!
  !-----------------------------------------------------------------------------------------------------
  PURE FUNCTION Add_pool_a_and_scalar(this,scalar) RESULT(pool3)
    IMPLICIT NONE
    ! 0.1 InOut
    CLASS(t_veg_pool_a), INTENT(in) :: this
    REAL(wp),            INTENT(in) :: scalar
    TYPE(t_veg_pool_a)              :: pool3
    ! 1.0 go
    pool3%leaf          = this%leaf         + scalar
    pool3%fine_root     = this%fine_root    + scalar
    pool3%coarse_root   = this%coarse_root  + scalar
    pool3%sap_wood      = this%sap_wood     + scalar
    pool3%heart_wood    = this%heart_wood   + scalar
    pool3%labile        = this%labile       + scalar
    pool3%reserve       = this%reserve      + scalar
    pool3%fruit         = this%fruit        + scalar
    pool3%seed_bed      = this%seed_bed     + scalar
  END FUNCTION Add_pool_a_and_scalar
  !-----------------------------------------------------------------------------------------------------
  !> + (t_veg_pool_a + 1D vector)
  !!
  !-----------------------------------------------------------------------------------------------------
  PURE FUNCTION Add_pool_a_and_vector(this,vector) RESULT(pool3)
    IMPLICIT NONE
    ! 0.1 InOut
    CLASS(t_veg_pool_a), INTENT(in) :: this
    REAL(wp),            INTENT(in) :: vector(:)
    TYPE(t_veg_pool_a)              :: pool3
    ! 1.0 go
    pool3%leaf          = this%leaf         + vector(:)
    pool3%fine_root     = this%fine_root    + vector(:)
    pool3%coarse_root   = this%coarse_root  + vector(:)
    pool3%sap_wood      = this%sap_wood     + vector(:)
    pool3%heart_wood    = this%heart_wood   + vector(:)
    pool3%labile        = this%labile       + vector(:)
    pool3%reserve       = this%reserve      + vector(:)
    pool3%fruit         = this%fruit        + vector(:)
    pool3%seed_bed      = this%seed_bed     + vector(:)
  END FUNCTION Add_pool_a_and_vector
  !-----------------------------------------------------------------------------------------------------
  !> - (t_veg_pool_a - t_veg_pool_a)
  !!
  !-----------------------------------------------------------------------------------------------------
  PURE FUNCTION Subtract_pool_a(this,pool2) RESULT(pool3)
    IMPLICIT NONE
    ! 0.1 InOut
    CLASS(t_veg_pool_a), INTENT(in) :: this, pool2
    TYPE(t_veg_pool_a)              :: pool3
    ! 1.0 go
    pool3%leaf          = this%leaf         - pool2%leaf
    pool3%fine_root     = this%fine_root    - pool2%fine_root
    pool3%coarse_root   = this%coarse_root  - pool2%coarse_root
    pool3%sap_wood      = this%sap_wood     - pool2%sap_wood
    pool3%heart_wood    = this%heart_wood   - pool2%heart_wood
    pool3%labile        = this%labile       - pool2%labile
    pool3%reserve       = this%reserve      - pool2%reserve
    pool3%fruit         = this%fruit        - pool2%fruit
    pool3%seed_bed      = this%seed_bed     - pool2%seed_bed
  END FUNCTION Subtract_pool_a
  !-----------------------------------------------------------------------------------------------------
  !> - (t_veg_pool_a - scalar)
  !!
  !-----------------------------------------------------------------------------------------------------
  PURE FUNCTION Subtract_pool_a_and_scalar(this,scalar) RESULT(pool3)
    IMPLICIT NONE
    ! 0.1 InOut
    CLASS(t_veg_pool_a), INTENT(in) :: this
    REAL(wp),            INTENT(in) :: scalar
    TYPE(t_veg_pool_a)              :: pool3
    ! 1.0 go
    pool3%leaf          = this%leaf         - scalar
    pool3%fine_root     = this%fine_root    - scalar
    pool3%coarse_root   = this%coarse_root  - scalar
    pool3%sap_wood      = this%sap_wood     - scalar
    pool3%heart_wood    = this%heart_wood   - scalar
    pool3%labile        = this%labile       - scalar
    pool3%reserve       = this%reserve      - scalar
    pool3%fruit         = this%fruit        - scalar
    pool3%seed_bed      = this%seed_bed     - scalar
  END FUNCTION Subtract_pool_a_and_scalar
  !-----------------------------------------------------------------------------------------------------
  !> - (t_veg_pool_a - 1D vector)
  !!
  !-----------------------------------------------------------------------------------------------------
  PURE FUNCTION Subtract_pool_a_and_vector(this,vector) RESULT(pool3)
    IMPLICIT NONE
    ! 0.1 InOut
    CLASS(t_veg_pool_a), INTENT(in) :: this
    REAL(wp),            INTENT(in) :: vector(:)
    TYPE(t_veg_pool_a)              :: pool3
    ! 1.0 go
    pool3%leaf          = this%leaf         - vector(:)
    pool3%fine_root     = this%fine_root    - vector(:)
    pool3%coarse_root   = this%coarse_root  - vector(:)
    pool3%sap_wood      = this%sap_wood     - vector(:)
    pool3%heart_wood    = this%heart_wood   - vector(:)
    pool3%labile        = this%labile       - vector(:)
    pool3%reserve       = this%reserve      - vector(:)
    pool3%fruit         = this%fruit        - vector(:)
    pool3%seed_bed      = this%seed_bed     - vector(:)
  END FUNCTION Subtract_pool_a_and_vector
  !-----------------------------------------------------------------------------------------------------
  !> * (t_veg_pool_a * t_veg_pool_a)
  !!
  !-----------------------------------------------------------------------------------------------------
  PURE FUNCTION Multiply_pool_a(this,pool2) RESULT(pool3)
    IMPLICIT NONE
    ! 0.1 InOut
    CLASS(t_veg_pool_a), INTENT(in) :: this, pool2
    TYPE(t_veg_pool_a)              :: pool3
    ! 1.0 go
    pool3%leaf          = this%leaf         * pool2%leaf
    pool3%fine_root     = this%fine_root    * pool2%fine_root
    pool3%coarse_root   = this%coarse_root  * pool2%coarse_root
    pool3%sap_wood      = this%sap_wood     * pool2%sap_wood
    pool3%heart_wood    = this%heart_wood   * pool2%heart_wood
    pool3%labile        = this%labile       * pool2%labile
    pool3%reserve       = this%reserve      * pool2%reserve
    pool3%fruit         = this%fruit        * pool2%fruit
    pool3%seed_bed      = this%seed_bed     * pool2%seed_bed
  END FUNCTION Multiply_pool_a
  !-----------------------------------------------------------------------------------------------------
  !> * (t_veg_pool_a * scalar)
  !!
  !-----------------------------------------------------------------------------------------------------
  PURE FUNCTION Multiply_pool_a_and_scalar(this,scalar) RESULT(pool3)
    IMPLICIT NONE
    ! 0.1 InOut
    CLASS(t_veg_pool_a), INTENT(in) :: this
    REAL(wp),            INTENT(in) :: scalar
    TYPE(t_veg_pool_a)              :: pool3
    ! 1.0 go
    pool3%leaf          = this%leaf         * scalar
    pool3%fine_root     = this%fine_root    * scalar
    pool3%coarse_root   = this%coarse_root  * scalar
    pool3%sap_wood      = this%sap_wood     * scalar
    pool3%heart_wood    = this%heart_wood   * scalar
    pool3%labile        = this%labile       * scalar
    pool3%reserve       = this%reserve      * scalar
    pool3%fruit         = this%fruit        * scalar
    pool3%seed_bed      = this%seed_bed     * scalar
  END FUNCTION Multiply_pool_a_and_scalar
  !-----------------------------------------------------------------------------------------------------
  !> * (scalar * t_veg_pool_a)
  !!
  !-----------------------------------------------------------------------------------------------------
  PURE FUNCTION Multiply_scalar_and_pool_a(scalar,this2) RESULT(pool3)
    IMPLICIT NONE
    ! 0.1 InOut
    REAL(wp),            INTENT(in) :: scalar
    CLASS(t_veg_pool_a), INTENT(in) :: this2
    TYPE(t_veg_pool_a)              :: pool3
    ! 1.0 go
    pool3%leaf          = scalar * this2%leaf
    pool3%fine_root     = scalar * this2%fine_root
    pool3%coarse_root   = scalar * this2%coarse_root
    pool3%sap_wood      = scalar * this2%sap_wood
    pool3%heart_wood    = scalar * this2%heart_wood
    pool3%labile        = scalar * this2%labile
    pool3%reserve       = scalar * this2%reserve
    pool3%fruit         = scalar * this2%fruit
    pool3%seed_bed      = scalar * this2%seed_bed
  END FUNCTION Multiply_scalar_and_pool_a
  !-----------------------------------------------------------------------------------------------------
  !> * (t_veg_pool_a * 1D vector)
  !!
  !-----------------------------------------------------------------------------------------------------
  PURE FUNCTION Multiply_pool_a_and_vector(this,vector) RESULT(pool3)
    IMPLICIT NONE
    ! 0.1 InOut
    CLASS(t_veg_pool_a), INTENT(in) :: this
    REAL(wp),            INTENT(in) :: vector(:)
    TYPE(t_veg_pool_a)              :: pool3
    ! 1.0 go
    pool3%leaf          = this%leaf         * vector(:)
    pool3%fine_root     = this%fine_root    * vector(:)
    pool3%coarse_root   = this%coarse_root  * vector(:)
    pool3%sap_wood      = this%sap_wood     * vector(:)
    pool3%heart_wood    = this%heart_wood   * vector(:)
    pool3%labile        = this%labile       * vector(:)
    pool3%reserve       = this%reserve      * vector(:)
    pool3%fruit         = this%fruit        * vector(:)
    pool3%seed_bed      = this%seed_bed     * vector(:)
  END FUNCTION Multiply_pool_a_and_vector
  !-----------------------------------------------------------------------------------------------------
  !> * (1D vector * t_veg_pool_a)
  !!
  !-----------------------------------------------------------------------------------------------------
  PURE FUNCTION Multiply_vector_and_pool_a(vector,this2) RESULT(pool3)
    IMPLICIT NONE
    ! 0.1 InOut
    REAL(wp),            INTENT(in) :: vector(:)
    CLASS(t_veg_pool_a), INTENT(in) :: this2
    TYPE(t_veg_pool_a)              :: pool3
    ! 1.0 go
    pool3%leaf          = vector(:) * this2%leaf
    pool3%fine_root     = vector(:) * this2%fine_root
    pool3%coarse_root   = vector(:) * this2%coarse_root
    pool3%sap_wood      = vector(:) * this2%sap_wood
    pool3%heart_wood    = vector(:) * this2%heart_wood
    pool3%labile        = vector(:) * this2%labile
    pool3%reserve       = vector(:) * this2%reserve
    pool3%fruit         = vector(:) * this2%fruit
    pool3%seed_bed      = vector(:) * this2%seed_bed
  END FUNCTION Multiply_vector_and_pool_a

  !-----------------------------------------------------------------------------------------------------
  !> Initialisation t_quincy_struct_veg_element with a scalar
  !!
  !-----------------------------------------------------------------------------------------------------
  PURE ELEMENTAL SUBROUTINE Init_t_quincy_struct_veg_element(this,s)

    CLASS(t_quincy_struct_veg_element), INTENT(inout) :: this 
    REAL(wp),                           INTENT(in)    :: s

    this%carbon%ptr(:,:)     = s
    this%nitrogen%ptr(:,:)   = s
    this%phosphorus%ptr(:,:) = s
    this%carbon13%ptr(:,:)   = s
    this%carbon14%ptr(:,:)   = s
    this%nitrogen15%ptr(:,:) = s
  END SUBROUTINE Init_t_quincy_struct_veg_element 

  !-----------------------------------------------------------------------------------------------------
  !> Initialisation t_quincy_struct_veg_pool with a scalar
  !!
  !-----------------------------------------------------------------------------------------------------
  PURE ELEMENTAL SUBROUTINE Init_t_quincy_struct_veg_pool(this,s)

    CLASS(t_quincy_struct_veg_pool),  INTENT(inout) :: this 
    REAL(wp),                         INTENT(in)    :: s

    CALL this%leaf%init(s)
    CALL this%fine_root%init(s)
    CALL this%coarse_root%init(s)
    CALL this%sap_wood%init(s)
    CALL this%heart_wood%init(s)
    CALL this%labile%init(s)
    CALL this%reserve%init(s)
    CALL this%fruit%init(s)
    CALL this%seed_bed%init(s)
  END SUBROUTINE Init_t_quincy_struct_veg_pool 


  ! ======================================================================================================= !
  !>Add all variables of an element TYPE to the process memory
  !>pool%Add_vars()
  !> 
  SUBROUTINE Add_element_variables_to_memory(this, prefix, suffix, lct_ids, model_id, unit)

    USE mo_jsb_varlist,       ONLY: BASIC, MEDIUM, FULL
    USE mo_jsb_io,            ONLY: grib_bits, t_cf, t_grib1, t_grib2
    USE mo_jsb_grid_class,    ONLY: t_jsb_grid, t_jsb_vgrid
    USE mo_jsb_grid,          ONLY: Get_grid, Get_vgrid
    USE mo_jsb_class,         ONLY: Get_model
    USE mo_jsb_model_class,   ONLY: t_jsb_model
    ! ----------------------------------------------------------------------------------------------------- !
    CLASS(t_quincy_struct_veg_element), INTENT(inout) :: this 
    CHARACTER(len=*),                   INTENT(in)    :: prefix          !> process name
    CHARACTER(len=*),                   INTENT(in)    :: suffix          !> tile name
    INTEGER,                            INTENT(in)    :: lct_ids !(:)    !> in jsbach4: vector of land cover type IDs
    INTEGER,                            INTENT(in)    :: model_id        !> model ID model\%id
    CHARACTER(len=*),                   INTENT(in)    :: unit            !> unit of the variables in the elements
    ! ----------------------------------------------------------------------------------------------------- !
    TYPE(t_jsb_model), POINTER  :: model
    TYPE(t_jsb_grid),  POINTER  :: hgrid                        ! Horizontal grid
    TYPE(t_jsb_vgrid), POINTER  :: surface                      ! Vertical grid
    INTEGER                     :: table
    CHARACTER(len=50)           :: var_name
    CHARACTER(len=100)          :: var_name_long
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':Add_element_variables_to_memory'
    ! ----------------------------------------------------------------------------------------------------- !
    model         => Get_model(model_id)
    table         = 1     ! in jsbach4 this line is: table = tables(1) 
    hgrid         => Get_grid(model%grid_id)
    surface       => Get_vgrid('surface')

    ! add variables to pool 
    var_name      = TRIM(this%name)//'_carbon'
    var_name_long = 'Carbon in vegetation '//TRIM(this%name)
    CALL this%Add_var(TRIM(var_name), this%carbon, &
      & hgrid, surface, &
      & t_cf(TRIM(var_name), TRIM(unit), TRIM(var_name_long)), &
      & t_grib1(table, 255, grib_bits), &
      & t_grib2(255, 255, 255, grib_bits), &
      & prefix, suffix, &
      & loutput = .TRUE., &
      & lrestart = .TRUE., &
      & initval_r = 0.0_wp)
    var_name      = TRIM(this%name)//'_nitrogen'
    var_name_long = 'Nitrogen in vegetation '//TRIM(this%name)
    CALL this%Add_var(TRIM(var_name), this%nitrogen, &
      & hgrid, surface, &
      & t_cf(TRIM(var_name), TRIM(unit), TRIM(var_name_long)), &
      & t_grib1(table, 255, grib_bits), &
      & t_grib2(255, 255, 255, grib_bits), &
      & prefix, suffix, &
      & loutput = .TRUE., &
      & lrestart = .TRUE., &
      & initval_r = 0.0_wp)
    var_name      = TRIM(this%name)//'_phosphorus'
    var_name_long = 'Phosphorus in vegetation '//TRIM(this%name)
    CALL this%Add_var(TRIM(var_name), this%phosphorus, &
      & hgrid, surface, &
      & t_cf(TRIM(var_name), TRIM(unit), TRIM(var_name_long)), &
      & t_grib1(table, 255, grib_bits), &
      & t_grib2(255, 255, 255, grib_bits), &
      & prefix, suffix, &
      & loutput = .TRUE., &
      & lrestart = .TRUE., &
      & initval_r = 0.0_wp)
    var_name      = TRIM(this%name)//'_carbon13'
    var_name_long = 'Carbon13 in vegetation '//TRIM(this%name)
    CALL this%Add_var(TRIM(var_name), this%carbon13, &
      & hgrid, surface, &
      & t_cf(TRIM(var_name), TRIM(unit), TRIM(var_name_long)), &
      & t_grib1(table, 255, grib_bits), &
      & t_grib2(255, 255, 255, grib_bits), &
      & prefix, suffix, &
      & loutput = .TRUE., &
      & lrestart = .TRUE., &
      & initval_r = 0.0_wp)
    var_name      = TRIM(this%name)//'_carbon14'
    var_name_long = 'Carbon14 in vegetation '//TRIM(this%name)
    CALL this%Add_var(TRIM(var_name), this%carbon14, &
      & hgrid, surface, &
      & t_cf(TRIM(var_name), TRIM(unit), TRIM(var_name_long)), &
      & t_grib1(table, 255, grib_bits), &
      & t_grib2(255, 255, 255, grib_bits), &
      & prefix, suffix, &
      & loutput = .TRUE., &
      & lrestart = .TRUE., &
      & initval_r = 0.0_wp)
    var_name      = TRIM(this%name)//'_nitrogen15'
    var_name_long = 'Nitrogen15 in vegetation '//TRIM(this%name)
    CALL this%Add_var(TRIM(var_name), this%nitrogen15, &
      & hgrid, surface, &
      & t_cf(TRIM(var_name), TRIM(unit), TRIM(var_name_long)), &
      & t_grib1(table, 255, grib_bits), &
      & t_grib2(255, 255, 255, grib_bits), &
      & prefix, suffix, &
      & loutput = .TRUE., &
      & lrestart = .TRUE., &
      & initval_r = 0.0_wp)
      
  END SUBROUTINE Add_element_variables_to_memory


  ! ======================================================================================================= !
  !> Add all elements of a pool TYPE to the process memory
  !> pool%Add_elements()
  !>
  SUBROUTINE Add_pool_elements_to_memory(this, prefix, suffix, lct_ids, model_id, unit)

    USE mo_jsb_varlist,       ONLY: BASIC, MEDIUM, FULL
    USE mo_jsb_io,            ONLY: grib_bits, t_cf, t_grib1, t_grib2
    USE mo_jsb_grid_class,    ONLY: t_jsb_grid, t_jsb_vgrid
    USE mo_jsb_grid,          ONLY: Get_grid, Get_vgrid
    USE mo_jsb_class,         ONLY: Get_model
    USE mo_jsb_model_class,   ONLY: t_jsb_model
    ! ----------------------------------------------------------------------------------------------------- !
    CLASS(t_quincy_struct_veg_pool),  INTENT(inout) :: this 
    CHARACTER(len=*),                 INTENT(in)    :: prefix          !> ?
    CHARACTER(len=*),                 INTENT(in)    :: suffix          !> ?
    INTEGER,                          INTENT(in)    :: lct_ids !(:)    !> in jsbach4: vector of land cover type IDs
    INTEGER,                          INTENT(in)    :: model_id        !> model ID model\%id
    CHARACTER(len=*),                 INTENT(in)    :: unit            !> unit of the variables in the elements
    ! ----------------------------------------------------------------------------------------------------- !
    TYPE(t_jsb_model), POINTER  :: model
    TYPE(t_jsb_grid),  POINTER  :: hgrid                        ! Horizontal grid
    TYPE(t_jsb_vgrid), POINTER  :: surface                      ! Vertical grid
    INTEGER                     :: table
    CHARACTER(len=40)           :: element_name
    CHARACTER(len=80)           :: element_name_long    
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':Add_pool_elements_to_memory'
    ! ----------------------------------------------------------------------------------------------------- !
    model         => Get_model(model_id)
    table         = 1     ! in jsbach4 this line is: table = tables(1) 
    hgrid         => Get_grid(model%grid_id)
    surface       => Get_vgrid('surface')

    !> add element pools
    !> 
    element_name      = TRIM(this%name)//'_leaf'
    element_name_long = 'Leaf element of vegetation '//TRIM(this%name)
    CALL this%Add_pool(TRIM(element_name), this%leaf, &
      & hgrid, surface, &
      & t_cf(TRIM(element_name), TRIM(unit), TRIM(element_name_long)), &
      & t_grib1(table, 255, grib_bits), &
      & t_grib2(255, 255, 255, grib_bits), &
      & prefix, suffix, &
      & loutput = .TRUE., &
      & lrestart = .TRUE.)
    element_name      = TRIM(this%name)//'_fine_root'
    element_name_long = 'Fine_root element of vegetation '//TRIM(this%name)
    CALL this%Add_pool(TRIM(element_name), this%fine_root, &
      & hgrid, surface, &
      & t_cf(TRIM(element_name), TRIM(unit), TRIM(element_name_long)), &
      & t_grib1(table, 255, grib_bits), &
      & t_grib2(255, 255, 255, grib_bits), &
      & prefix, suffix, &
      & loutput = .TRUE., &
      & lrestart = .TRUE.)
    element_name      = TRIM(this%name)//'_coarse_root'
    element_name_long = 'Coarse_root element of vegetation '//TRIM(this%name)
    CALL this%Add_pool(TRIM(element_name), this%coarse_root, &
      & hgrid, surface, &
      & t_cf(TRIM(element_name), TRIM(unit), TRIM(element_name_long)), &
      & t_grib1(table, 255, grib_bits), &
      & t_grib2(255, 255, 255, grib_bits), &
      & prefix, suffix, &
      & loutput = .TRUE., &
      & lrestart = .TRUE.)
    element_name      = TRIM(this%name)//'_sap_wood'
    element_name_long = 'Sap_wood element of vegetation '//TRIM(this%name)
    CALL this%Add_pool(TRIM(element_name), this%sap_wood, &
      & hgrid, surface, &
      & t_cf(TRIM(element_name), TRIM(unit), TRIM(element_name_long)), &
      & t_grib1(table, 255, grib_bits), &
      & t_grib2(255, 255, 255, grib_bits), &
      & prefix, suffix, &
      & loutput = .TRUE., &
      & lrestart = .TRUE.)
    element_name      = TRIM(this%name)//'_heart_wood'
    element_name_long = 'Heart_wood element of vegetation '//TRIM(this%name)
    CALL this%Add_pool(TRIM(element_name), this%heart_wood, &
      & hgrid, surface, &
      & t_cf(TRIM(element_name), TRIM(unit), TRIM(element_name_long)), &
      & t_grib1(table, 255, grib_bits), &
      & t_grib2(255, 255, 255, grib_bits), &
      & prefix, suffix, &
      & loutput = .TRUE., &
      & lrestart = .TRUE.)
    element_name      = TRIM(this%name)//'_labile'
    element_name_long = 'Labile element of vegetation '//TRIM(this%name)
    CALL this%Add_pool(TRIM(element_name), this%labile, &
      & hgrid, surface, &
      & t_cf(TRIM(element_name), TRIM(unit), TRIM(element_name_long)), &
      & t_grib1(table, 255, grib_bits), &
      & t_grib2(255, 255, 255, grib_bits), &
      & prefix, suffix, &
      & loutput = .TRUE., &
      & lrestart = .TRUE.)
    element_name      = TRIM(this%name)//'_reserve'
    element_name_long = 'Reserve element of vegetation '//TRIM(this%name)
    CALL this%Add_pool(TRIM(element_name), this%reserve, &
      & hgrid, surface, &
      & t_cf(TRIM(element_name), TRIM(unit), TRIM(element_name_long)), &
      & t_grib1(table, 255, grib_bits), &
      & t_grib2(255, 255, 255, grib_bits), &
      & prefix, suffix, &
      & loutput = .TRUE., &
      & lrestart = .TRUE.)
    element_name      = TRIM(this%name)//'_fruit'
    element_name_long = 'Fruit element of vegetation '//TRIM(this%name)
    CALL this%Add_pool(TRIM(element_name), this%fruit, &
      & hgrid, surface, &
      & t_cf(TRIM(element_name), TRIM(unit), TRIM(element_name_long)), &
      & t_grib1(table, 255, grib_bits), &
      & t_grib2(255, 255, 255, grib_bits), &
      & prefix, suffix, &
      & loutput = .TRUE., &
      & lrestart = .TRUE.)
    element_name      = TRIM(this%name)//'_seed_bed'
    element_name_long = 'Seed_bed element of vegetation '//TRIM(this%name)
    CALL this%Add_pool(TRIM(element_name), this%seed_bed, &
      & hgrid, surface, &
      & t_cf(TRIM(element_name), TRIM(unit), TRIM(element_name_long)), &
      & t_grib1(table, 255, grib_bits), &
      & t_grib2(255, 255, 255, grib_bits), &
      & prefix, suffix, &
      & loutput = .TRUE., &
      & lrestart = .TRUE.)
  END SUBROUTINE Add_pool_elements_to_memory

  !-----------------------------------------------------------------------------------------------------
  !> ALLOCATE variables in veg_pool_a and veg_element_a structures
  !!
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE Init_veg_pools_elements_a(ics, ice, nc, iblk)

    INTEGER, INTENT(in)     :: ics, ice, nc, iblk             !< information about dimensions of variables in memory

    INTEGER                 :: npoints                        !< number of points in the chunk

    ! set dimension to number of points in this chunk
    npoints = nc

    !!  ALLOCATE variables
    ALLOCATE(veg_pool%leaf%carbon(npoints))
    ALLOCATE(veg_pool%leaf%nitrogen(npoints))
    ALLOCATE(veg_pool%leaf%phosphorus(npoints))
    ALLOCATE(veg_pool%leaf%carbon13(npoints))
    ALLOCATE(veg_pool%leaf%carbon14(npoints))
    ALLOCATE(veg_pool%leaf%nitrogen15(npoints))
    ALLOCATE(veg_pool%fine_root%carbon(npoints))
    ALLOCATE(veg_pool%fine_root%nitrogen(npoints))
    ALLOCATE(veg_pool%fine_root%phosphorus(npoints))
    ALLOCATE(veg_pool%fine_root%carbon13(npoints))
    ALLOCATE(veg_pool%fine_root%carbon14(npoints))
    ALLOCATE(veg_pool%fine_root%nitrogen15(npoints))
    ALLOCATE(veg_pool%coarse_root%carbon(npoints))
    ALLOCATE(veg_pool%coarse_root%nitrogen(npoints))
    ALLOCATE(veg_pool%coarse_root%phosphorus(npoints))
    ALLOCATE(veg_pool%coarse_root%carbon13(npoints))
    ALLOCATE(veg_pool%coarse_root%carbon14(npoints))
    ALLOCATE(veg_pool%coarse_root%nitrogen15(npoints))
    ALLOCATE(veg_pool%sap_wood%carbon(npoints))
    ALLOCATE(veg_pool%sap_wood%nitrogen(npoints))
    ALLOCATE(veg_pool%sap_wood%phosphorus(npoints))
    ALLOCATE(veg_pool%sap_wood%carbon13(npoints))
    ALLOCATE(veg_pool%sap_wood%carbon14(npoints))
    ALLOCATE(veg_pool%sap_wood%nitrogen15(npoints))
    ALLOCATE(veg_pool%heart_wood%carbon(npoints))
    ALLOCATE(veg_pool%heart_wood%nitrogen(npoints))
    ALLOCATE(veg_pool%heart_wood%phosphorus(npoints))
    ALLOCATE(veg_pool%heart_wood%carbon13(npoints))
    ALLOCATE(veg_pool%heart_wood%carbon14(npoints))
    ALLOCATE(veg_pool%heart_wood%nitrogen15(npoints))
    ALLOCATE(veg_pool%labile%carbon(npoints))
    ALLOCATE(veg_pool%labile%nitrogen(npoints))
    ALLOCATE(veg_pool%labile%phosphorus(npoints))
    ALLOCATE(veg_pool%labile%carbon13(npoints))
    ALLOCATE(veg_pool%labile%carbon14(npoints))
    ALLOCATE(veg_pool%labile%nitrogen15(npoints))
    ALLOCATE(veg_pool%reserve%carbon(npoints))
    ALLOCATE(veg_pool%reserve%nitrogen(npoints))
    ALLOCATE(veg_pool%reserve%phosphorus(npoints))
    ALLOCATE(veg_pool%reserve%carbon13(npoints))
    ALLOCATE(veg_pool%reserve%carbon14(npoints))
    ALLOCATE(veg_pool%reserve%nitrogen15(npoints))
    ALLOCATE(veg_pool%fruit%carbon(npoints))
    ALLOCATE(veg_pool%fruit%nitrogen(npoints))
    ALLOCATE(veg_pool%fruit%phosphorus(npoints))
    ALLOCATE(veg_pool%fruit%carbon13(npoints))
    ALLOCATE(veg_pool%fruit%carbon14(npoints))
    ALLOCATE(veg_pool%fruit%nitrogen15(npoints))
    ALLOCATE(veg_pool%seed_bed%carbon(npoints))
    ALLOCATE(veg_pool%seed_bed%nitrogen(npoints))
    ALLOCATE(veg_pool%seed_bed%phosphorus(npoints))
    ALLOCATE(veg_pool%seed_bed%carbon13(npoints))
    ALLOCATE(veg_pool%seed_bed%carbon14(npoints))
    ALLOCATE(veg_pool%seed_bed%nitrogen15(npoints))

    ALLOCATE(veg_growth%leaf%carbon(npoints))
    ALLOCATE(veg_growth%leaf%nitrogen(npoints))
    ALLOCATE(veg_growth%leaf%phosphorus(npoints))
    ALLOCATE(veg_growth%leaf%carbon13(npoints))
    ALLOCATE(veg_growth%leaf%carbon14(npoints))
    ALLOCATE(veg_growth%leaf%nitrogen15(npoints))
    ALLOCATE(veg_growth%fine_root%carbon(npoints))
    ALLOCATE(veg_growth%fine_root%nitrogen(npoints))
    ALLOCATE(veg_growth%fine_root%phosphorus(npoints))
    ALLOCATE(veg_growth%fine_root%carbon13(npoints))
    ALLOCATE(veg_growth%fine_root%carbon14(npoints))
    ALLOCATE(veg_growth%fine_root%nitrogen15(npoints))
    ALLOCATE(veg_growth%coarse_root%carbon(npoints))
    ALLOCATE(veg_growth%coarse_root%nitrogen(npoints))
    ALLOCATE(veg_growth%coarse_root%phosphorus(npoints))
    ALLOCATE(veg_growth%coarse_root%carbon13(npoints))
    ALLOCATE(veg_growth%coarse_root%carbon14(npoints))
    ALLOCATE(veg_growth%coarse_root%nitrogen15(npoints))
    ALLOCATE(veg_growth%sap_wood%carbon(npoints))
    ALLOCATE(veg_growth%sap_wood%nitrogen(npoints))
    ALLOCATE(veg_growth%sap_wood%phosphorus(npoints))
    ALLOCATE(veg_growth%sap_wood%carbon13(npoints))
    ALLOCATE(veg_growth%sap_wood%carbon14(npoints))
    ALLOCATE(veg_growth%sap_wood%nitrogen15(npoints))
    ALLOCATE(veg_growth%heart_wood%carbon(npoints))
    ALLOCATE(veg_growth%heart_wood%nitrogen(npoints))
    ALLOCATE(veg_growth%heart_wood%phosphorus(npoints))
    ALLOCATE(veg_growth%heart_wood%carbon13(npoints))
    ALLOCATE(veg_growth%heart_wood%carbon14(npoints))
    ALLOCATE(veg_growth%heart_wood%nitrogen15(npoints))
    ALLOCATE(veg_growth%labile%carbon(npoints))
    ALLOCATE(veg_growth%labile%nitrogen(npoints))
    ALLOCATE(veg_growth%labile%phosphorus(npoints))
    ALLOCATE(veg_growth%labile%carbon13(npoints))
    ALLOCATE(veg_growth%labile%carbon14(npoints))
    ALLOCATE(veg_growth%labile%nitrogen15(npoints))
    ALLOCATE(veg_growth%reserve%carbon(npoints))
    ALLOCATE(veg_growth%reserve%nitrogen(npoints))
    ALLOCATE(veg_growth%reserve%phosphorus(npoints))
    ALLOCATE(veg_growth%reserve%carbon13(npoints))
    ALLOCATE(veg_growth%reserve%carbon14(npoints))
    ALLOCATE(veg_growth%reserve%nitrogen15(npoints))
    ALLOCATE(veg_growth%fruit%carbon(npoints))
    ALLOCATE(veg_growth%fruit%nitrogen(npoints))
    ALLOCATE(veg_growth%fruit%phosphorus(npoints))
    ALLOCATE(veg_growth%fruit%carbon13(npoints))
    ALLOCATE(veg_growth%fruit%carbon14(npoints))
    ALLOCATE(veg_growth%fruit%nitrogen15(npoints))
    ALLOCATE(veg_growth%seed_bed%carbon(npoints))
    ALLOCATE(veg_growth%seed_bed%nitrogen(npoints))
    ALLOCATE(veg_growth%seed_bed%phosphorus(npoints))
    ALLOCATE(veg_growth%seed_bed%carbon13(npoints))
    ALLOCATE(veg_growth%seed_bed%carbon14(npoints))
    ALLOCATE(veg_growth%seed_bed%nitrogen15(npoints))

    ALLOCATE(veg_litterfall%leaf%carbon(npoints))
    ALLOCATE(veg_litterfall%leaf%nitrogen(npoints))
    ALLOCATE(veg_litterfall%leaf%phosphorus(npoints))
    ALLOCATE(veg_litterfall%leaf%carbon13(npoints))
    ALLOCATE(veg_litterfall%leaf%carbon14(npoints))
    ALLOCATE(veg_litterfall%leaf%nitrogen15(npoints))
    ALLOCATE(veg_litterfall%fine_root%carbon(npoints))
    ALLOCATE(veg_litterfall%fine_root%nitrogen(npoints))
    ALLOCATE(veg_litterfall%fine_root%phosphorus(npoints))
    ALLOCATE(veg_litterfall%fine_root%carbon13(npoints))
    ALLOCATE(veg_litterfall%fine_root%carbon14(npoints))
    ALLOCATE(veg_litterfall%fine_root%nitrogen15(npoints))
    ALLOCATE(veg_litterfall%coarse_root%carbon(npoints))
    ALLOCATE(veg_litterfall%coarse_root%nitrogen(npoints))
    ALLOCATE(veg_litterfall%coarse_root%phosphorus(npoints))
    ALLOCATE(veg_litterfall%coarse_root%carbon13(npoints))
    ALLOCATE(veg_litterfall%coarse_root%carbon14(npoints))
    ALLOCATE(veg_litterfall%coarse_root%nitrogen15(npoints))
    ALLOCATE(veg_litterfall%sap_wood%carbon(npoints))
    ALLOCATE(veg_litterfall%sap_wood%nitrogen(npoints))
    ALLOCATE(veg_litterfall%sap_wood%phosphorus(npoints))
    ALLOCATE(veg_litterfall%sap_wood%carbon13(npoints))
    ALLOCATE(veg_litterfall%sap_wood%carbon14(npoints))
    ALLOCATE(veg_litterfall%sap_wood%nitrogen15(npoints))
    ALLOCATE(veg_litterfall%heart_wood%carbon(npoints))
    ALLOCATE(veg_litterfall%heart_wood%nitrogen(npoints))
    ALLOCATE(veg_litterfall%heart_wood%phosphorus(npoints))
    ALLOCATE(veg_litterfall%heart_wood%carbon13(npoints))
    ALLOCATE(veg_litterfall%heart_wood%carbon14(npoints))
    ALLOCATE(veg_litterfall%heart_wood%nitrogen15(npoints))
    ALLOCATE(veg_litterfall%labile%carbon(npoints))
    ALLOCATE(veg_litterfall%labile%nitrogen(npoints))
    ALLOCATE(veg_litterfall%labile%phosphorus(npoints))
    ALLOCATE(veg_litterfall%labile%carbon13(npoints))
    ALLOCATE(veg_litterfall%labile%carbon14(npoints))
    ALLOCATE(veg_litterfall%labile%nitrogen15(npoints))
    ALLOCATE(veg_litterfall%reserve%carbon(npoints))
    ALLOCATE(veg_litterfall%reserve%nitrogen(npoints))
    ALLOCATE(veg_litterfall%reserve%phosphorus(npoints))
    ALLOCATE(veg_litterfall%reserve%carbon13(npoints))
    ALLOCATE(veg_litterfall%reserve%carbon14(npoints))
    ALLOCATE(veg_litterfall%reserve%nitrogen15(npoints))
    ALLOCATE(veg_litterfall%fruit%carbon(npoints))
    ALLOCATE(veg_litterfall%fruit%nitrogen(npoints))
    ALLOCATE(veg_litterfall%fruit%phosphorus(npoints))
    ALLOCATE(veg_litterfall%fruit%carbon13(npoints))
    ALLOCATE(veg_litterfall%fruit%carbon14(npoints))
    ALLOCATE(veg_litterfall%fruit%nitrogen15(npoints))
    ALLOCATE(veg_litterfall%seed_bed%carbon(npoints))
    ALLOCATE(veg_litterfall%seed_bed%nitrogen(npoints))
    ALLOCATE(veg_litterfall%seed_bed%phosphorus(npoints))
    ALLOCATE(veg_litterfall%seed_bed%carbon13(npoints))
    ALLOCATE(veg_litterfall%seed_bed%carbon14(npoints))
    ALLOCATE(veg_litterfall%seed_bed%nitrogen15(npoints))

    ALLOCATE(veg_falloc%leaf%carbon(npoints))
    ALLOCATE(veg_falloc%leaf%nitrogen(npoints))
    ALLOCATE(veg_falloc%leaf%phosphorus(npoints))
    ALLOCATE(veg_falloc%leaf%carbon13(npoints))
    ALLOCATE(veg_falloc%leaf%carbon14(npoints))
    ALLOCATE(veg_falloc%leaf%nitrogen15(npoints))
    ALLOCATE(veg_falloc%fine_root%carbon(npoints))
    ALLOCATE(veg_falloc%fine_root%nitrogen(npoints))
    ALLOCATE(veg_falloc%fine_root%phosphorus(npoints))
    ALLOCATE(veg_falloc%fine_root%carbon13(npoints))
    ALLOCATE(veg_falloc%fine_root%carbon14(npoints))
    ALLOCATE(veg_falloc%fine_root%nitrogen15(npoints))
    ALLOCATE(veg_falloc%coarse_root%carbon(npoints))
    ALLOCATE(veg_falloc%coarse_root%nitrogen(npoints))
    ALLOCATE(veg_falloc%coarse_root%phosphorus(npoints))
    ALLOCATE(veg_falloc%coarse_root%carbon13(npoints))
    ALLOCATE(veg_falloc%coarse_root%carbon14(npoints))
    ALLOCATE(veg_falloc%coarse_root%nitrogen15(npoints))
    ALLOCATE(veg_falloc%sap_wood%carbon(npoints))
    ALLOCATE(veg_falloc%sap_wood%nitrogen(npoints))
    ALLOCATE(veg_falloc%sap_wood%phosphorus(npoints))
    ALLOCATE(veg_falloc%sap_wood%carbon13(npoints))
    ALLOCATE(veg_falloc%sap_wood%carbon14(npoints))
    ALLOCATE(veg_falloc%sap_wood%nitrogen15(npoints))
    ALLOCATE(veg_falloc%heart_wood%carbon(npoints))
    ALLOCATE(veg_falloc%heart_wood%nitrogen(npoints))
    ALLOCATE(veg_falloc%heart_wood%phosphorus(npoints))
    ALLOCATE(veg_falloc%heart_wood%carbon13(npoints))
    ALLOCATE(veg_falloc%heart_wood%carbon14(npoints))
    ALLOCATE(veg_falloc%heart_wood%nitrogen15(npoints))
    ALLOCATE(veg_falloc%labile%carbon(npoints))
    ALLOCATE(veg_falloc%labile%nitrogen(npoints))
    ALLOCATE(veg_falloc%labile%phosphorus(npoints))
    ALLOCATE(veg_falloc%labile%carbon13(npoints))
    ALLOCATE(veg_falloc%labile%carbon14(npoints))
    ALLOCATE(veg_falloc%labile%nitrogen15(npoints))
    ALLOCATE(veg_falloc%reserve%carbon(npoints))
    ALLOCATE(veg_falloc%reserve%nitrogen(npoints))
    ALLOCATE(veg_falloc%reserve%phosphorus(npoints))
    ALLOCATE(veg_falloc%reserve%carbon13(npoints))
    ALLOCATE(veg_falloc%reserve%carbon14(npoints))
    ALLOCATE(veg_falloc%reserve%nitrogen15(npoints))
    ALLOCATE(veg_falloc%fruit%carbon(npoints))
    ALLOCATE(veg_falloc%fruit%nitrogen(npoints))
    ALLOCATE(veg_falloc%fruit%phosphorus(npoints))
    ALLOCATE(veg_falloc%fruit%carbon13(npoints))
    ALLOCATE(veg_falloc%fruit%carbon14(npoints))
    ALLOCATE(veg_falloc%fruit%nitrogen15(npoints))
    ALLOCATE(veg_falloc%seed_bed%carbon(npoints))
    ALLOCATE(veg_falloc%seed_bed%nitrogen(npoints))
    ALLOCATE(veg_falloc%seed_bed%phosphorus(npoints))
    ALLOCATE(veg_falloc%seed_bed%carbon13(npoints))
    ALLOCATE(veg_falloc%seed_bed%carbon14(npoints))
    ALLOCATE(veg_falloc%seed_bed%nitrogen15(npoints))

    ALLOCATE(veg_total_flux%leaf%carbon(npoints))
    ALLOCATE(veg_total_flux%leaf%nitrogen(npoints))
    ALLOCATE(veg_total_flux%leaf%phosphorus(npoints))
    ALLOCATE(veg_total_flux%leaf%carbon13(npoints))
    ALLOCATE(veg_total_flux%leaf%carbon14(npoints))
    ALLOCATE(veg_total_flux%leaf%nitrogen15(npoints))
    ALLOCATE(veg_total_flux%fine_root%carbon(npoints))
    ALLOCATE(veg_total_flux%fine_root%nitrogen(npoints))
    ALLOCATE(veg_total_flux%fine_root%phosphorus(npoints))
    ALLOCATE(veg_total_flux%fine_root%carbon13(npoints))
    ALLOCATE(veg_total_flux%fine_root%carbon14(npoints))
    ALLOCATE(veg_total_flux%fine_root%nitrogen15(npoints))
    ALLOCATE(veg_total_flux%coarse_root%carbon(npoints))
    ALLOCATE(veg_total_flux%coarse_root%nitrogen(npoints))
    ALLOCATE(veg_total_flux%coarse_root%phosphorus(npoints))
    ALLOCATE(veg_total_flux%coarse_root%carbon13(npoints))
    ALLOCATE(veg_total_flux%coarse_root%carbon14(npoints))
    ALLOCATE(veg_total_flux%coarse_root%nitrogen15(npoints))
    ALLOCATE(veg_total_flux%sap_wood%carbon(npoints))
    ALLOCATE(veg_total_flux%sap_wood%nitrogen(npoints))
    ALLOCATE(veg_total_flux%sap_wood%phosphorus(npoints))
    ALLOCATE(veg_total_flux%sap_wood%carbon13(npoints))
    ALLOCATE(veg_total_flux%sap_wood%carbon14(npoints))
    ALLOCATE(veg_total_flux%sap_wood%nitrogen15(npoints))
    ALLOCATE(veg_total_flux%heart_wood%carbon(npoints))
    ALLOCATE(veg_total_flux%heart_wood%nitrogen(npoints))
    ALLOCATE(veg_total_flux%heart_wood%phosphorus(npoints))
    ALLOCATE(veg_total_flux%heart_wood%carbon13(npoints))
    ALLOCATE(veg_total_flux%heart_wood%carbon14(npoints))
    ALLOCATE(veg_total_flux%heart_wood%nitrogen15(npoints))
    ALLOCATE(veg_total_flux%labile%carbon(npoints))
    ALLOCATE(veg_total_flux%labile%nitrogen(npoints))
    ALLOCATE(veg_total_flux%labile%phosphorus(npoints))
    ALLOCATE(veg_total_flux%labile%carbon13(npoints))
    ALLOCATE(veg_total_flux%labile%carbon14(npoints))
    ALLOCATE(veg_total_flux%labile%nitrogen15(npoints))
    ALLOCATE(veg_total_flux%reserve%carbon(npoints))
    ALLOCATE(veg_total_flux%reserve%nitrogen(npoints))
    ALLOCATE(veg_total_flux%reserve%phosphorus(npoints))
    ALLOCATE(veg_total_flux%reserve%carbon13(npoints))
    ALLOCATE(veg_total_flux%reserve%carbon14(npoints))
    ALLOCATE(veg_total_flux%reserve%nitrogen15(npoints))
    ALLOCATE(veg_total_flux%fruit%carbon(npoints))
    ALLOCATE(veg_total_flux%fruit%nitrogen(npoints))
    ALLOCATE(veg_total_flux%fruit%phosphorus(npoints))
    ALLOCATE(veg_total_flux%fruit%carbon13(npoints))
    ALLOCATE(veg_total_flux%fruit%carbon14(npoints))
    ALLOCATE(veg_total_flux%fruit%nitrogen15(npoints))
    ALLOCATE(veg_total_flux%seed_bed%carbon(npoints))
    ALLOCATE(veg_total_flux%seed_bed%nitrogen(npoints))
    ALLOCATE(veg_total_flux%seed_bed%phosphorus(npoints))
    ALLOCATE(veg_total_flux%seed_bed%carbon13(npoints))
    ALLOCATE(veg_total_flux%seed_bed%carbon14(npoints))
    ALLOCATE(veg_total_flux%seed_bed%nitrogen15(npoints))

    ALLOCATE(veg_exudation%carbon(npoints))
    ALLOCATE(veg_exudation%nitrogen(npoints))
    ALLOCATE(veg_exudation%phosphorus(npoints))
    ALLOCATE(veg_exudation%carbon13(npoints))
    ALLOCATE(veg_exudation%carbon14(npoints))
    ALLOCATE(veg_exudation%nitrogen15(npoints))

    ALLOCATE(veg_establishment%carbon(npoints))
    ALLOCATE(veg_establishment%nitrogen(npoints))
    ALLOCATE(veg_establishment%phosphorus(npoints))
    ALLOCATE(veg_establishment%carbon13(npoints))
    ALLOCATE(veg_establishment%carbon14(npoints))
    ALLOCATE(veg_establishment%nitrogen15(npoints))

    ALLOCATE(veg_reserve_use%carbon(npoints))
    ALLOCATE(veg_reserve_use%nitrogen(npoints))
    ALLOCATE(veg_reserve_use%phosphorus(npoints))
    ALLOCATE(veg_reserve_use%carbon13(npoints))
    ALLOCATE(veg_reserve_use%carbon14(npoints))
    ALLOCATE(veg_reserve_use%nitrogen15(npoints))

    ALLOCATE(veg_total_biomass%carbon(npoints))
    ALLOCATE(veg_total_biomass%nitrogen(npoints))
    ALLOCATE(veg_total_biomass%phosphorus(npoints))
    ALLOCATE(veg_total_biomass%carbon13(npoints))
    ALLOCATE(veg_total_biomass%carbon14(npoints))
    ALLOCATE(veg_total_biomass%nitrogen15(npoints))

  END SUBROUTINE Init_veg_pools_elements_a

  !-----------------------------------------------------------------------------------------------------
  !> DEALLOCATE variables from veg_pool_a and veg_element_a structures at the program end
  !!
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE Deallocate_var_veg_pools_elements_a

    !< 1.0  DEALLOCATE variables
    !! 
    DEALLOCATE(veg_pool%leaf%carbon)
    DEALLOCATE(veg_pool%leaf%nitrogen)
    DEALLOCATE(veg_pool%leaf%phosphorus)
    DEALLOCATE(veg_pool%leaf%carbon13)
    DEALLOCATE(veg_pool%leaf%carbon14)
    DEALLOCATE(veg_pool%leaf%nitrogen15)
    DEALLOCATE(veg_pool%fine_root%carbon)
    DEALLOCATE(veg_pool%fine_root%nitrogen)
    DEALLOCATE(veg_pool%fine_root%phosphorus)
    DEALLOCATE(veg_pool%fine_root%carbon13)
    DEALLOCATE(veg_pool%fine_root%carbon14)
    DEALLOCATE(veg_pool%fine_root%nitrogen15)
    DEALLOCATE(veg_pool%coarse_root%carbon)
    DEALLOCATE(veg_pool%coarse_root%nitrogen)
    DEALLOCATE(veg_pool%coarse_root%phosphorus)
    DEALLOCATE(veg_pool%coarse_root%carbon13)
    DEALLOCATE(veg_pool%coarse_root%carbon14)
    DEALLOCATE(veg_pool%coarse_root%nitrogen15)
    DEALLOCATE(veg_pool%sap_wood%carbon)
    DEALLOCATE(veg_pool%sap_wood%nitrogen)
    DEALLOCATE(veg_pool%sap_wood%phosphorus)
    DEALLOCATE(veg_pool%sap_wood%carbon13)
    DEALLOCATE(veg_pool%sap_wood%carbon14)
    DEALLOCATE(veg_pool%sap_wood%nitrogen15)
    DEALLOCATE(veg_pool%heart_wood%carbon)
    DEALLOCATE(veg_pool%heart_wood%nitrogen)
    DEALLOCATE(veg_pool%heart_wood%phosphorus)
    DEALLOCATE(veg_pool%heart_wood%carbon13)
    DEALLOCATE(veg_pool%heart_wood%carbon14)
    DEALLOCATE(veg_pool%heart_wood%nitrogen15)
    DEALLOCATE(veg_pool%labile%carbon)
    DEALLOCATE(veg_pool%labile%nitrogen)
    DEALLOCATE(veg_pool%labile%phosphorus)
    DEALLOCATE(veg_pool%labile%carbon13)
    DEALLOCATE(veg_pool%labile%carbon14)
    DEALLOCATE(veg_pool%labile%nitrogen15)
    DEALLOCATE(veg_pool%reserve%carbon)
    DEALLOCATE(veg_pool%reserve%nitrogen)
    DEALLOCATE(veg_pool%reserve%phosphorus)
    DEALLOCATE(veg_pool%reserve%carbon13)
    DEALLOCATE(veg_pool%reserve%carbon14)
    DEALLOCATE(veg_pool%reserve%nitrogen15)
    DEALLOCATE(veg_pool%fruit%carbon)
    DEALLOCATE(veg_pool%fruit%nitrogen)
    DEALLOCATE(veg_pool%fruit%phosphorus)
    DEALLOCATE(veg_pool%fruit%carbon13)
    DEALLOCATE(veg_pool%fruit%carbon14)
    DEALLOCATE(veg_pool%fruit%nitrogen15)
    DEALLOCATE(veg_pool%seed_bed%carbon)
    DEALLOCATE(veg_pool%seed_bed%nitrogen)
    DEALLOCATE(veg_pool%seed_bed%phosphorus)
    DEALLOCATE(veg_pool%seed_bed%carbon13)
    DEALLOCATE(veg_pool%seed_bed%carbon14)
    DEALLOCATE(veg_pool%seed_bed%nitrogen15)

    DEALLOCATE(veg_growth%leaf%carbon)
    DEALLOCATE(veg_growth%leaf%nitrogen)
    DEALLOCATE(veg_growth%leaf%phosphorus)
    DEALLOCATE(veg_growth%leaf%carbon13)
    DEALLOCATE(veg_growth%leaf%carbon14)
    DEALLOCATE(veg_growth%leaf%nitrogen15)
    DEALLOCATE(veg_growth%fine_root%carbon)
    DEALLOCATE(veg_growth%fine_root%nitrogen)
    DEALLOCATE(veg_growth%fine_root%phosphorus)
    DEALLOCATE(veg_growth%fine_root%carbon13)
    DEALLOCATE(veg_growth%fine_root%carbon14)
    DEALLOCATE(veg_growth%fine_root%nitrogen15)
    DEALLOCATE(veg_growth%coarse_root%carbon)
    DEALLOCATE(veg_growth%coarse_root%nitrogen)
    DEALLOCATE(veg_growth%coarse_root%phosphorus)
    DEALLOCATE(veg_growth%coarse_root%carbon13)
    DEALLOCATE(veg_growth%coarse_root%carbon14)
    DEALLOCATE(veg_growth%coarse_root%nitrogen15)
    DEALLOCATE(veg_growth%sap_wood%carbon)
    DEALLOCATE(veg_growth%sap_wood%nitrogen)
    DEALLOCATE(veg_growth%sap_wood%phosphorus)
    DEALLOCATE(veg_growth%sap_wood%carbon13)
    DEALLOCATE(veg_growth%sap_wood%carbon14)
    DEALLOCATE(veg_growth%sap_wood%nitrogen15)
    DEALLOCATE(veg_growth%heart_wood%carbon)
    DEALLOCATE(veg_growth%heart_wood%nitrogen)
    DEALLOCATE(veg_growth%heart_wood%phosphorus)
    DEALLOCATE(veg_growth%heart_wood%carbon13)
    DEALLOCATE(veg_growth%heart_wood%carbon14)
    DEALLOCATE(veg_growth%heart_wood%nitrogen15)
    DEALLOCATE(veg_growth%labile%carbon)
    DEALLOCATE(veg_growth%labile%nitrogen)
    DEALLOCATE(veg_growth%labile%phosphorus)
    DEALLOCATE(veg_growth%labile%carbon13)
    DEALLOCATE(veg_growth%labile%carbon14)
    DEALLOCATE(veg_growth%labile%nitrogen15)
    DEALLOCATE(veg_growth%reserve%carbon)
    DEALLOCATE(veg_growth%reserve%nitrogen)
    DEALLOCATE(veg_growth%reserve%phosphorus)
    DEALLOCATE(veg_growth%reserve%carbon13)
    DEALLOCATE(veg_growth%reserve%carbon14)
    DEALLOCATE(veg_growth%reserve%nitrogen15)
    DEALLOCATE(veg_growth%fruit%carbon)
    DEALLOCATE(veg_growth%fruit%nitrogen)
    DEALLOCATE(veg_growth%fruit%phosphorus)
    DEALLOCATE(veg_growth%fruit%carbon13)
    DEALLOCATE(veg_growth%fruit%carbon14)
    DEALLOCATE(veg_growth%fruit%nitrogen15)
    DEALLOCATE(veg_growth%seed_bed%carbon)
    DEALLOCATE(veg_growth%seed_bed%nitrogen)
    DEALLOCATE(veg_growth%seed_bed%phosphorus)
    DEALLOCATE(veg_growth%seed_bed%carbon13)
    DEALLOCATE(veg_growth%seed_bed%carbon14)
    DEALLOCATE(veg_growth%seed_bed%nitrogen15)

    DEALLOCATE(veg_litterfall%leaf%carbon)
    DEALLOCATE(veg_litterfall%leaf%nitrogen)
    DEALLOCATE(veg_litterfall%leaf%phosphorus)
    DEALLOCATE(veg_litterfall%leaf%carbon13)
    DEALLOCATE(veg_litterfall%leaf%carbon14)
    DEALLOCATE(veg_litterfall%leaf%nitrogen15)
    DEALLOCATE(veg_litterfall%fine_root%carbon)
    DEALLOCATE(veg_litterfall%fine_root%nitrogen)
    DEALLOCATE(veg_litterfall%fine_root%phosphorus)
    DEALLOCATE(veg_litterfall%fine_root%carbon13)
    DEALLOCATE(veg_litterfall%fine_root%carbon14)
    DEALLOCATE(veg_litterfall%fine_root%nitrogen15)
    DEALLOCATE(veg_litterfall%coarse_root%carbon)
    DEALLOCATE(veg_litterfall%coarse_root%nitrogen)
    DEALLOCATE(veg_litterfall%coarse_root%phosphorus)
    DEALLOCATE(veg_litterfall%coarse_root%carbon13)
    DEALLOCATE(veg_litterfall%coarse_root%carbon14)
    DEALLOCATE(veg_litterfall%coarse_root%nitrogen15)
    DEALLOCATE(veg_litterfall%sap_wood%carbon)
    DEALLOCATE(veg_litterfall%sap_wood%nitrogen)
    DEALLOCATE(veg_litterfall%sap_wood%phosphorus)
    DEALLOCATE(veg_litterfall%sap_wood%carbon13)
    DEALLOCATE(veg_litterfall%sap_wood%carbon14)
    DEALLOCATE(veg_litterfall%sap_wood%nitrogen15)
    DEALLOCATE(veg_litterfall%heart_wood%carbon)
    DEALLOCATE(veg_litterfall%heart_wood%nitrogen)
    DEALLOCATE(veg_litterfall%heart_wood%phosphorus)
    DEALLOCATE(veg_litterfall%heart_wood%carbon13)
    DEALLOCATE(veg_litterfall%heart_wood%carbon14)
    DEALLOCATE(veg_litterfall%heart_wood%nitrogen15)
    DEALLOCATE(veg_litterfall%labile%carbon)
    DEALLOCATE(veg_litterfall%labile%nitrogen)
    DEALLOCATE(veg_litterfall%labile%phosphorus)
    DEALLOCATE(veg_litterfall%labile%carbon13)
    DEALLOCATE(veg_litterfall%labile%carbon14)
    DEALLOCATE(veg_litterfall%labile%nitrogen15)
    DEALLOCATE(veg_litterfall%reserve%carbon)
    DEALLOCATE(veg_litterfall%reserve%nitrogen)
    DEALLOCATE(veg_litterfall%reserve%phosphorus)
    DEALLOCATE(veg_litterfall%reserve%carbon13)
    DEALLOCATE(veg_litterfall%reserve%carbon14)
    DEALLOCATE(veg_litterfall%reserve%nitrogen15)
    DEALLOCATE(veg_litterfall%fruit%carbon)
    DEALLOCATE(veg_litterfall%fruit%nitrogen)
    DEALLOCATE(veg_litterfall%fruit%phosphorus)
    DEALLOCATE(veg_litterfall%fruit%carbon13)
    DEALLOCATE(veg_litterfall%fruit%carbon14)
    DEALLOCATE(veg_litterfall%fruit%nitrogen15)
    DEALLOCATE(veg_litterfall%seed_bed%carbon)
    DEALLOCATE(veg_litterfall%seed_bed%nitrogen)
    DEALLOCATE(veg_litterfall%seed_bed%phosphorus)
    DEALLOCATE(veg_litterfall%seed_bed%carbon13)
    DEALLOCATE(veg_litterfall%seed_bed%carbon14)
    DEALLOCATE(veg_litterfall%seed_bed%nitrogen15)

    DEALLOCATE(veg_falloc%leaf%carbon)
    DEALLOCATE(veg_falloc%leaf%nitrogen)
    DEALLOCATE(veg_falloc%leaf%phosphorus)
    DEALLOCATE(veg_falloc%leaf%carbon13)
    DEALLOCATE(veg_falloc%leaf%carbon14)
    DEALLOCATE(veg_falloc%leaf%nitrogen15)
    DEALLOCATE(veg_falloc%fine_root%carbon)
    DEALLOCATE(veg_falloc%fine_root%nitrogen)
    DEALLOCATE(veg_falloc%fine_root%phosphorus)
    DEALLOCATE(veg_falloc%fine_root%carbon13)
    DEALLOCATE(veg_falloc%fine_root%carbon14)
    DEALLOCATE(veg_falloc%fine_root%nitrogen15)
    DEALLOCATE(veg_falloc%coarse_root%carbon)
    DEALLOCATE(veg_falloc%coarse_root%nitrogen)
    DEALLOCATE(veg_falloc%coarse_root%phosphorus)
    DEALLOCATE(veg_falloc%coarse_root%carbon13)
    DEALLOCATE(veg_falloc%coarse_root%carbon14)
    DEALLOCATE(veg_falloc%coarse_root%nitrogen15)
    DEALLOCATE(veg_falloc%sap_wood%carbon)
    DEALLOCATE(veg_falloc%sap_wood%nitrogen)
    DEALLOCATE(veg_falloc%sap_wood%phosphorus)
    DEALLOCATE(veg_falloc%sap_wood%carbon13)
    DEALLOCATE(veg_falloc%sap_wood%carbon14)
    DEALLOCATE(veg_falloc%sap_wood%nitrogen15)
    DEALLOCATE(veg_falloc%heart_wood%carbon)
    DEALLOCATE(veg_falloc%heart_wood%nitrogen)
    DEALLOCATE(veg_falloc%heart_wood%phosphorus)
    DEALLOCATE(veg_falloc%heart_wood%carbon13)
    DEALLOCATE(veg_falloc%heart_wood%carbon14)
    DEALLOCATE(veg_falloc%heart_wood%nitrogen15)
    DEALLOCATE(veg_falloc%labile%carbon)
    DEALLOCATE(veg_falloc%labile%nitrogen)
    DEALLOCATE(veg_falloc%labile%phosphorus)
    DEALLOCATE(veg_falloc%labile%carbon13)
    DEALLOCATE(veg_falloc%labile%carbon14)
    DEALLOCATE(veg_falloc%labile%nitrogen15)
    DEALLOCATE(veg_falloc%reserve%carbon)
    DEALLOCATE(veg_falloc%reserve%nitrogen)
    DEALLOCATE(veg_falloc%reserve%phosphorus)
    DEALLOCATE(veg_falloc%reserve%carbon13)
    DEALLOCATE(veg_falloc%reserve%carbon14)
    DEALLOCATE(veg_falloc%reserve%nitrogen15)
    DEALLOCATE(veg_falloc%fruit%carbon)
    DEALLOCATE(veg_falloc%fruit%nitrogen)
    DEALLOCATE(veg_falloc%fruit%phosphorus)
    DEALLOCATE(veg_falloc%fruit%carbon13)
    DEALLOCATE(veg_falloc%fruit%carbon14)
    DEALLOCATE(veg_falloc%fruit%nitrogen15)
    DEALLOCATE(veg_falloc%seed_bed%carbon)
    DEALLOCATE(veg_falloc%seed_bed%nitrogen)
    DEALLOCATE(veg_falloc%seed_bed%phosphorus)
    DEALLOCATE(veg_falloc%seed_bed%carbon13)
    DEALLOCATE(veg_falloc%seed_bed%carbon14)
    DEALLOCATE(veg_falloc%seed_bed%nitrogen15)

    DEALLOCATE(veg_total_flux%leaf%carbon)
    DEALLOCATE(veg_total_flux%leaf%nitrogen)
    DEALLOCATE(veg_total_flux%leaf%phosphorus)
    DEALLOCATE(veg_total_flux%leaf%carbon13)
    DEALLOCATE(veg_total_flux%leaf%carbon14)
    DEALLOCATE(veg_total_flux%leaf%nitrogen15)
    DEALLOCATE(veg_total_flux%fine_root%carbon)
    DEALLOCATE(veg_total_flux%fine_root%nitrogen)
    DEALLOCATE(veg_total_flux%fine_root%phosphorus)
    DEALLOCATE(veg_total_flux%fine_root%carbon13)
    DEALLOCATE(veg_total_flux%fine_root%carbon14)
    DEALLOCATE(veg_total_flux%fine_root%nitrogen15)
    DEALLOCATE(veg_total_flux%coarse_root%carbon)
    DEALLOCATE(veg_total_flux%coarse_root%nitrogen)
    DEALLOCATE(veg_total_flux%coarse_root%phosphorus)
    DEALLOCATE(veg_total_flux%coarse_root%carbon13)
    DEALLOCATE(veg_total_flux%coarse_root%carbon14)
    DEALLOCATE(veg_total_flux%coarse_root%nitrogen15)
    DEALLOCATE(veg_total_flux%sap_wood%carbon)
    DEALLOCATE(veg_total_flux%sap_wood%nitrogen)
    DEALLOCATE(veg_total_flux%sap_wood%phosphorus)
    DEALLOCATE(veg_total_flux%sap_wood%carbon13)
    DEALLOCATE(veg_total_flux%sap_wood%carbon14)
    DEALLOCATE(veg_total_flux%sap_wood%nitrogen15)
    DEALLOCATE(veg_total_flux%heart_wood%carbon)
    DEALLOCATE(veg_total_flux%heart_wood%nitrogen)
    DEALLOCATE(veg_total_flux%heart_wood%phosphorus)
    DEALLOCATE(veg_total_flux%heart_wood%carbon13)
    DEALLOCATE(veg_total_flux%heart_wood%carbon14)
    DEALLOCATE(veg_total_flux%heart_wood%nitrogen15)
    DEALLOCATE(veg_total_flux%labile%carbon)
    DEALLOCATE(veg_total_flux%labile%nitrogen)
    DEALLOCATE(veg_total_flux%labile%phosphorus)
    DEALLOCATE(veg_total_flux%labile%carbon13)
    DEALLOCATE(veg_total_flux%labile%carbon14)
    DEALLOCATE(veg_total_flux%labile%nitrogen15)
    DEALLOCATE(veg_total_flux%reserve%carbon)
    DEALLOCATE(veg_total_flux%reserve%nitrogen)
    DEALLOCATE(veg_total_flux%reserve%phosphorus)
    DEALLOCATE(veg_total_flux%reserve%carbon13)
    DEALLOCATE(veg_total_flux%reserve%carbon14)
    DEALLOCATE(veg_total_flux%reserve%nitrogen15)
    DEALLOCATE(veg_total_flux%fruit%carbon)
    DEALLOCATE(veg_total_flux%fruit%nitrogen)
    DEALLOCATE(veg_total_flux%fruit%phosphorus)
    DEALLOCATE(veg_total_flux%fruit%carbon13)
    DEALLOCATE(veg_total_flux%fruit%carbon14)
    DEALLOCATE(veg_total_flux%fruit%nitrogen15)
    DEALLOCATE(veg_total_flux%seed_bed%carbon)
    DEALLOCATE(veg_total_flux%seed_bed%nitrogen)
    DEALLOCATE(veg_total_flux%seed_bed%phosphorus)
    DEALLOCATE(veg_total_flux%seed_bed%carbon13)
    DEALLOCATE(veg_total_flux%seed_bed%carbon14)
    DEALLOCATE(veg_total_flux%seed_bed%nitrogen15)

    DEALLOCATE(veg_exudation%carbon)
    DEALLOCATE(veg_exudation%nitrogen)
    DEALLOCATE(veg_exudation%phosphorus)
    DEALLOCATE(veg_exudation%carbon13)
    DEALLOCATE(veg_exudation%carbon14)
    DEALLOCATE(veg_exudation%nitrogen15)

    DEALLOCATE(veg_establishment%carbon)
    DEALLOCATE(veg_establishment%nitrogen)
    DEALLOCATE(veg_establishment%phosphorus)
    DEALLOCATE(veg_establishment%carbon13)
    DEALLOCATE(veg_establishment%carbon14)
    DEALLOCATE(veg_establishment%nitrogen15)

    DEALLOCATE(veg_reserve_use%carbon)
    DEALLOCATE(veg_reserve_use%nitrogen)
    DEALLOCATE(veg_reserve_use%phosphorus)
    DEALLOCATE(veg_reserve_use%carbon13)
    DEALLOCATE(veg_reserve_use%carbon14)
    DEALLOCATE(veg_reserve_use%nitrogen15)

    DEALLOCATE(veg_total_biomass%carbon)
    DEALLOCATE(veg_total_biomass%nitrogen)
    DEALLOCATE(veg_total_biomass%phosphorus)
    DEALLOCATE(veg_total_biomass%carbon13)
    DEALLOCATE(veg_total_biomass%carbon14)
    DEALLOCATE(veg_total_biomass%nitrogen15)

  END SUBROUTINE Deallocate_var_veg_pools_elements_a

!! #################################################################################################################################
!! quincy-no-jsb_end
!! #################################################################################################################################
#else
!! #################################################################################################################################
!! code specific to QUINCY with JSBACH4 (until the end of this module)
!! #################################################################################################################################

  USE mo_kind,                   ONLY: wp
  USE mo_jsb_var_class,          ONLY: t_jsb_var_real2d

  IMPLICIT NONE

  ! all is public except the modename; code is used by the mo_veg_memory_class only
  PUBLIC
  PRIVATE :: modname


  CHARACTER(len=*), PARAMETER :: modname = 'mo_veg_memory_class_util'

!CONTAINS



!! #################################################################################################################################
!! code specific to QUINCY with JSBACH4 ends here
!! #################################################################################################################################
#endif

END MODULE mo_veg_memory_class_util
