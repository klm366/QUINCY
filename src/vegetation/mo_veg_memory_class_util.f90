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
  USE mo_exception,              ONLY: finish, message
  USE mo_util,                   ONLY: One_of
  USE mo_jsb_utils_iface,        ONLY: assign_if_present_allocatable
  USE mo_jsb_var_class,          ONLY: t_jsb_var_real2d
  USE mo_jsb_pool_class,         ONLY: t_jsb_pool_elem2d, ELEM_C, ELEM_N, ELEM_P, ELEM_C13, ELEM_C14, ELEM_N15
  USE mo_jsb_subset,             ONLY: t_subset, ON_DOMAIN, ON_CHUNK

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: all_elements, Create_veg_tissues_pool
  PUBLIC :: t_veg_main_pool, t_veg_tissues_pool, t_veg_elements_pool

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
    PROCEDURE :: Get_element_name_by_id     => t_veg_elements_pool_get_element_name_by_id
    PROCEDURE :: Init                       => t_veg_elements_pool_init
    PROCEDURE :: Init_op_ov_hlp             => t_veg_elements_pool_op_ov_hlp_init             ! init hlp object for operator overload
    PROCEDURE :: Associate_domain_var_ptr   => t_veg_elements_pool_associate_domain_variable_pointers
    PROCEDURE :: Set_to_scalar_on_chunk     => t_veg_elements_pool_set_all_elements_to_scalar_value_on_chunk  ! e.g., set all to zero
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
    PROCEDURE :: Get_element_name_by_id     => t_veg_tissues_pool_get_element_name_by_id
    PROCEDURE :: Init                       => t_veg_tissues_pool_init                      ! pool_object created with FUNCTION Create_veg_tissues_pool()
    PROCEDURE :: Init_op_ov_hlp             => t_veg_tissues_pool_op_ov_hlp_create_init     ! create and init hlp object for operator overload
    PROCEDURE :: Associate_domain_var_ptr   => t_veg_tissues_pool_associate_domain_variable_pointers
    PROCEDURE :: Set_to_scalar_on_chunk     => t_veg_tissues_pool_set_all_elements_to_scalar_value_on_chunk  ! e.g., set all to zero
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
    TYPE(t_veg_tissues_pool),   POINTER :: veg_total_flux   ! @TODO currently not used (Dec 2020)
    !! allocation factors: fractional allocation from growth%labile to growth%tissue_types
    TYPE(t_veg_tissues_pool),   POINTER :: veg_falloc
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

  CHARACTER(len=*), PARAMETER :: modname = 'mo_veg_memory_class_util'


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

  ! ======================================================================================================= !
  !>set all elements of the veg_elements_pool to the scalar value
  !> 
  SUBROUTINE t_veg_elements_pool_set_all_elements_to_scalar_value_on_chunk(this, scalar)

    CLASS(t_veg_elements_pool), INTENT(inout) :: this
    REAL(wp)                  , INTENT(in)    :: scalar
    
    this%carbon(:)       = scalar
    this%nitrogen(:)     = scalar
    this%phosphorus(:)   = scalar
    this%carbon13(:)     = scalar
    this%carbon14(:)     = scalar
    this%nitrogen15(:)   = scalar

  END SUBROUTINE t_veg_elements_pool_set_all_elements_to_scalar_value_on_chunk

  ! ======================================================================================================= !
  !>set all elements of the veg_elements_pool to the scalar value
  !> 
  SUBROUTINE t_veg_tissues_pool_set_all_elements_to_scalar_value_on_chunk(this, scalar)

    CLASS(t_veg_tissues_pool), INTENT(inout) :: this
    REAL(wp)                 , INTENT(in)    :: scalar
    
    ! elements
    this%carbon(:)       = scalar
    this%nitrogen(:)     = scalar
    this%phosphorus(:)   = scalar
    this%carbon13(:)     = scalar
    this%carbon14(:)     = scalar
    this%nitrogen15(:)   = scalar
    ! veg_elements_pool
    CALL this%leaf%Set_to_scalar_on_chunk(scalar)
    CALL this%fine_root%Set_to_scalar_on_chunk(scalar)
    CALL this%coarse_root%Set_to_scalar_on_chunk(scalar)
    CALL this%sap_wood%Set_to_scalar_on_chunk(scalar)
    CALL this%heart_wood%Set_to_scalar_on_chunk(scalar)
    CALL this%labile%Set_to_scalar_on_chunk(scalar)
    CALL this%reserve%Set_to_scalar_on_chunk(scalar)
    CALL this%fruit%Set_to_scalar_on_chunk(scalar)
    CALL this%seed_bed%Set_to_scalar_on_chunk(scalar)

  END SUBROUTINE t_veg_tissues_pool_set_all_elements_to_scalar_value_on_chunk  

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
    END IF
    !> nitrogen
    IF (One_of(ELEM_N, element_list_loc) > 0) THEN
      CALL this%Add_element('nitrogen', 'N', ELEM_N, dim=2)
    END IF
    !> phosphorus
    IF (One_of(ELEM_P, element_list_loc) > 0) THEN
      CALL this%Add_element('phosphorus', 'P', ELEM_P, dim=2)
    END IF
    !> carbon13
    IF (One_of(ELEM_C13, element_list_loc) > 0) THEN
      CALL this%Add_element('carbon13', 'C13', ELEM_C13, dim=2)
    END IF
    !> carbon14
    IF (One_of(ELEM_C14, element_list_loc) > 0) THEN
      CALL this%Add_element('carbon14', 'C14', ELEM_C14, dim=2)
    END IF
    !> nitrogen15
    IF (One_of(ELEM_N15, element_list_loc) > 0) THEN
      CALL this%Add_element('nitrogen15', 'N15', ELEM_N15, dim=2)
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
      IF (ASSOCIATED(this%carbon))     NULLIFY(this%carbon)
      this%carbon     => this%op_ov_hlp_1d_c
      ALLOCATE(this%op_ov_hlp_1d_n(nc))
      IF (ASSOCIATED(this%nitrogen))   NULLIFY(this%nitrogen)
      this%nitrogen   => this%op_ov_hlp_1d_n
      ALLOCATE(this%op_ov_hlp_1d_p(nc))
      IF (ASSOCIATED(this%phosphorus)) NULLIFY(this%phosphorus)
      this%phosphorus => this%op_ov_hlp_1d_p
      ALLOCATE(this%op_ov_hlp_1d_c13(nc))
      IF (ASSOCIATED(this%carbon13))   NULLIFY(this%carbon13)
      this%carbon13   => this%op_ov_hlp_1d_c13
      ALLOCATE(this%op_ov_hlp_1d_c14(nc))
      IF (ASSOCIATED(this%carbon14))   NULLIFY(this%carbon14)
      this%carbon14   => this%op_ov_hlp_1d_c14
      ALLOCATE(this%op_ov_hlp_1d_n15(nc))
      IF (ASSOCIATED(this%nitrogen15)) NULLIFY(this%nitrogen15)
      this%nitrogen15 => this%op_ov_hlp_1d_n15
    ELSE
      ALLOCATE(this%op_ov_hlp_2d_c(nc,nb))
      IF (ASSOCIATED(this%domain_c))   NULLIFY(this%domain_c)
      this%domain_c   => this%op_ov_hlp_2d_c
      ALLOCATE(this%op_ov_hlp_2d_n(nc,nb))
      IF (ASSOCIATED(this%domain_n))   NULLIFY(this%domain_n)
      this%domain_n   => this%op_ov_hlp_2d_n
      ALLOCATE(this%op_ov_hlp_2d_p(nc,nb))
      IF (ASSOCIATED(this%domain_p))   NULLIFY(this%domain_p)
      this%domain_p   => this%op_ov_hlp_2d_p
      ALLOCATE(this%op_ov_hlp_2d_c13(nc,nb))
      IF (ASSOCIATED(this%domain_c13)) NULLIFY(this%domain_c13)
      this%domain_c13 => this%op_ov_hlp_2d_c13
      ALLOCATE(this%op_ov_hlp_2d_c14(nc,nb))
      IF (ASSOCIATED(this%domain_c14)) NULLIFY(this%domain_c14)
      this%domain_c14 => this%op_ov_hlp_2d_c14
      ALLOCATE(this%op_ov_hlp_2d_n15(nc,nb))
      IF (ASSOCIATED(this%domain_n15)) NULLIFY(this%domain_n15)
      this%domain_n15 => this%op_ov_hlp_2d_n15
    ENDIF

  END SUBROUTINE t_veg_elements_pool_op_ov_hlp_init

  !-----------------------------------------------------------------------------------------------------
  !> set element-variables pointers to the domain pointer ptr2d(:,:) in the element list this%element_list(i)%p
  !! t_veg_elements_pool
  !! this need to be called for the VEG_ pools used in the veg_init() routine (or any other *_init() where these pools are used)
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE t_veg_elements_pool_associate_domain_variable_pointers(this)

    CLASS(t_veg_elements_pool),  INTENT(inout) :: this

    TYPE(t_subset)              :: sub    ! ON_CHUNK / ON_DOMAIN information
    INTEGER                     :: i      ! looping
    CHARACTER(len=*), PARAMETER :: routine = modname//':t_veg_elements_pool_associate_domain_variable_pointers'

    sub = this%Get_subset()

    ! test working ON_DOMAIN
    IF (sub%type == ON_CHUNK) THEN
      CALL finish(routine, ' May be called only in the process init rouines, i.e., when working ON_DOMAIN.')
    END IF

    ! loop over elements in element_list(:) of this element pool
    ! associate REAL element vectors with %ptr2d
    IF (ALLOCATED(this%element_list)) THEN
      DO i = 1,SIZE(this%element_list)
        SELECT TYPE (element => this%element_list(i)%p)
        CLASS IS (t_jsb_var_real2d)
          SELECT CASE(element%element_id)
          CASE (ELEM_C)
            IF (ASSOCIATED(this%domain_c)) NULLIFY(this%domain_c)
            this%domain_c    => element%ptr2d
          CASE (ELEM_N)
            IF (ASSOCIATED(this%domain_n)) NULLIFY(this%domain_n)
            this%domain_n    => element%ptr2d
          CASE (ELEM_P)
            IF (ASSOCIATED(this%domain_p))   NULLIFY(this%domain_p)
            this%domain_p    => element%ptr2d
          CASE (ELEM_C13)
            IF (ASSOCIATED(this%domain_c13)) NULLIFY(this%domain_c13)
            this%domain_c13  => element%ptr2d
          CASE (ELEM_C14)
            IF (ASSOCIATED(this%domain_c14)) NULLIFY(this%domain_c14)
            this%domain_c14  => element%ptr2d
          CASE (ELEM_N15)
            IF (ASSOCIATED(this%domain_n15)) NULLIFY(this%domain_n15)
            this%domain_n15  => element%ptr2d
          END SELECT
        END SELECT
      END DO
    END IF

  END SUBROUTINE t_veg_elements_pool_associate_domain_variable_pointers

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
    END IF
    !> nitrogen
    IF (One_of(ELEM_N, element_list_loc) > 0) THEN
      CALL this%Add_element('nitrogen', 'N', ELEM_N, dim=2)
    END IF
    !> phosphorus
    IF (One_of(ELEM_P, element_list_loc) > 0) THEN
      CALL this%Add_element('phosphorus', 'P', ELEM_P, dim=2)
    END IF
    !> carbon13
    IF (One_of(ELEM_C13, element_list_loc) > 0) THEN
      CALL this%Add_element('carbon13', 'C13', ELEM_C13, dim=2)
    END IF
    !> carbon14
    IF (One_of(ELEM_C14, element_list_loc) > 0) THEN
      CALL this%Add_element('carbon14', 'C14', ELEM_C14, dim=2)
    END IF
    !> nitrogen15
    IF (One_of(ELEM_N15, element_list_loc) > 0) THEN
      CALL this%Add_element('nitrogen15', 'N15', ELEM_N15, dim=2)
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
      IF (ASSOCIATED(this%carbon))     NULLIFY(this%carbon)
      this%carbon     => this%op_ov_hlp_1d_c
      ALLOCATE(this%op_ov_hlp_1d_n(nc))
      IF (ASSOCIATED(this%nitrogen))   NULLIFY(this%nitrogen)
      this%nitrogen   => this%op_ov_hlp_1d_n
      ALLOCATE(this%op_ov_hlp_1d_p(nc))
      IF (ASSOCIATED(this%phosphorus)) NULLIFY(this%phosphorus)
      this%phosphorus => this%op_ov_hlp_1d_p
      ALLOCATE(this%op_ov_hlp_1d_c13(nc))
      IF (ASSOCIATED(this%carbon13))   NULLIFY(this%carbon13)
      this%carbon13   => this%op_ov_hlp_1d_c13
      ALLOCATE(this%op_ov_hlp_1d_c14(nc))
      IF (ASSOCIATED(this%carbon14))   NULLIFY(this%carbon14)
      this%carbon14   => this%op_ov_hlp_1d_c14
      ALLOCATE(this%op_ov_hlp_1d_n15(nc))
      IF (ASSOCIATED(this%nitrogen15)) NULLIFY(this%nitrogen15)
      this%nitrogen15 => this%op_ov_hlp_1d_n15
    ELSE
      ALLOCATE(this%op_ov_hlp_2d_c(nc,nb))
      IF (ASSOCIATED(this%domain_c))   NULLIFY(this%domain_c)
      this%domain_c   => this%op_ov_hlp_2d_c
      ALLOCATE(this%op_ov_hlp_2d_n(nc,nb))
      IF (ASSOCIATED(this%domain_n))   NULLIFY(this%domain_n)
      this%domain_n   => this%op_ov_hlp_2d_n
      ALLOCATE(this%op_ov_hlp_2d_p(nc,nb))
      IF (ASSOCIATED(this%domain_p))   NULLIFY(this%domain_p)
      this%domain_p   => this%op_ov_hlp_2d_p
      ALLOCATE(this%op_ov_hlp_2d_c13(nc,nb))
      IF (ASSOCIATED(this%domain_c13)) NULLIFY(this%domain_c13)
      this%domain_c13 => this%op_ov_hlp_2d_c13
      ALLOCATE(this%op_ov_hlp_2d_c14(nc,nb))
      IF (ASSOCIATED(this%domain_c14)) NULLIFY(this%domain_c14)
      this%domain_c14 => this%op_ov_hlp_2d_c14
      ALLOCATE(this%op_ov_hlp_2d_n15(nc,nb))
      IF (ASSOCIATED(this%domain_n15)) NULLIFY(this%domain_n15)
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
  !> set element-variables pointers to the domain pointer ptr2d(:,:) in the element list this%element_list(i)%p
  !! t_veg_tissues_pool
  !! this need to be called for the VEG_ pools used in the veg_init() routine (or any other *_init() where these pools are used)
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE t_veg_tissues_pool_associate_domain_variable_pointers(this)

    CLASS(t_veg_tissues_pool),   INTENT(inout) :: this

    TYPE(t_subset)              :: sub    ! ON_CHUNK / ON_DOMAIN information
    INTEGER                     :: i      ! looping
    CHARACTER(len=*), PARAMETER :: routine = modname//':t_veg_tissues_pool_associate_domain_variable_pointers'

    sub = this%Get_subset()

    ! test working ON_DOMAIN
    IF (sub%type == ON_CHUNK) THEN
      CALL finish(routine, ' May be called only in the process init rouines, i.e., when working ON_DOMAIN.')
    END IF

    ! loop over elements in element_list(:) of this tissue pool
    ! associate REAL element vectors with %ptr2d
    IF (ALLOCATED(this%element_list)) THEN
      DO i = 1,SIZE(this%element_list)
        SELECT TYPE (element => this%element_list(i)%p)
        CLASS IS (t_jsb_var_real2d)
          SELECT CASE(element%element_id)
          CASE (ELEM_C)
            IF (ASSOCIATED(this%domain_c)) NULLIFY(this%domain_c)
            this%domain_c    => element%ptr2d
          CASE (ELEM_N)
            IF (ASSOCIATED(this%domain_n)) NULLIFY(this%domain_n)
            this%domain_n    => element%ptr2d
          CASE (ELEM_P)
            IF (ASSOCIATED(this%domain_p))   NULLIFY(this%domain_p)
            this%domain_p    => element%ptr2d
          CASE (ELEM_C13)
            IF (ASSOCIATED(this%domain_c13)) NULLIFY(this%domain_c13)
            this%domain_c13  => element%ptr2d
          CASE (ELEM_C14)
            IF (ASSOCIATED(this%domain_c14)) NULLIFY(this%domain_c14)
            this%domain_c14  => element%ptr2d
          CASE (ELEM_N15)
            IF (ASSOCIATED(this%domain_n15)) NULLIFY(this%domain_n15)
            this%domain_n15  => element%ptr2d
          END SELECT
        END SELECT
      END DO
    END IF
    
    ! associate pointers for all element pools of this tissue pool
    CALL this%leaf%Associate_domain_var_ptr()
    CALL this%fine_root%Associate_domain_var_ptr()
    CALL this%coarse_root%Associate_domain_var_ptr()
    CALL this%sap_wood%Associate_domain_var_ptr()
    CALL this%heart_wood%Associate_domain_var_ptr()
    CALL this%labile%Associate_domain_var_ptr()
    CALL this%reserve%Associate_domain_var_ptr()
    CALL this%fruit%Associate_domain_var_ptr()
    CALL this%seed_bed%Associate_domain_var_ptr()

  END SUBROUTINE t_veg_tissues_pool_associate_domain_variable_pointers

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
    IF (sub%type /= ON_CHUNK) THEN  ! get domain dimension not available via Get_subset() during process init
      sub%nc = SIZE(struct1%domain_c(:,:), 1)
      sub%nb = SIZE(struct1%domain_c(:,:), 2)
    ENDIF
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
    IF (sub%type /= ON_CHUNK) THEN  ! get domain dimension not available via Get_subset() during process init
      sub%nc = SIZE(struct1%domain_c(:,:), 1)
      sub%nb = SIZE(struct1%domain_c(:,:), 2)
    ENDIF
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
    IF (sub%type /= ON_CHUNK) THEN  ! get domain dimension not available via Get_subset() during process init
      sub%nc = SIZE(struct1%domain_c(:,:), 1)
      sub%nb = SIZE(struct1%domain_c(:,:), 2)
    ENDIF
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
    IF (sub%type /= ON_CHUNK) THEN  ! get domain dimension not available via Get_subset() during process init
      sub%nc = SIZE(struct1%domain_c(:,:), 1)
      sub%nb = SIZE(struct1%domain_c(:,:), 2)
    ENDIF
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
    IF (sub%type /= ON_CHUNK) THEN  ! get domain dimension not available via Get_subset() during process init
      sub%nc = SIZE(struct1%domain_c(:,:), 1)
      sub%nb = SIZE(struct1%domain_c(:,:), 2)
    ENDIF
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
    IF (sub%type /= ON_CHUNK) THEN  ! get domain dimension not available via Get_subset() during process init
      sub%nc = SIZE(struct1%domain_c(:,:), 1)
      sub%nb = SIZE(struct1%domain_c(:,:), 2)
    ENDIF
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
    IF (sub%type /= ON_CHUNK) THEN  ! get domain dimension not available via Get_subset() during process init
      sub%nc = SIZE(struct1%domain_c(:,:), 1)
      sub%nb = SIZE(struct1%domain_c(:,:), 2)
    ENDIF
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
    IF (sub%type /= ON_CHUNK) THEN  ! get domain dimension not available via Get_subset() during process init
      sub%nc = SIZE(struct1%domain_c(:,:), 1)
      sub%nb = SIZE(struct1%domain_c(:,:), 2)
    ENDIF
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
    IF (sub%type /= ON_CHUNK) THEN  ! get domain dimension not available via Get_subset() during process init
      sub%nc = SIZE(struct1%domain_c(:,:), 1)
      sub%nb = SIZE(struct1%domain_c(:,:), 2)
    ENDIF
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
    IF (sub%type /= ON_CHUNK) THEN  ! get domain dimension not available via Get_subset() during process init
      sub%nc = SIZE(struct1%domain_c(:,:), 1)
      sub%nb = SIZE(struct1%domain_c(:,:), 2)
    ENDIF
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
    IF (sub%type /= ON_CHUNK) THEN  ! get domain dimension not available via Get_subset() during process init
      sub%nc = SIZE(struct1%domain_c(:,:), 1)
      sub%nb = SIZE(struct1%domain_c(:,:), 2)
    ENDIF
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
    IF (sub%type /= ON_CHUNK) THEN  ! get domain dimension not available via Get_subset() during process init
      sub%nc = SIZE(struct1%domain_c(:,:), 1)
      sub%nb = SIZE(struct1%domain_c(:,:), 2)
    ENDIF
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
    IF (sub%type /= ON_CHUNK) THEN  ! get domain dimension not available via Get_subset() during process init
      sub%nc = SIZE(struct1%domain_c(:,:), 1)
      sub%nb = SIZE(struct1%domain_c(:,:), 2)
    ENDIF
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
    IF (sub%type /= ON_CHUNK) THEN  ! get domain dimension not available via Get_subset() during process init
      sub%nc = SIZE(struct1%domain_c(:,:), 1)
      sub%nb = SIZE(struct1%domain_c(:,:), 2)
    ENDIF
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
    IF (sub%type /= ON_CHUNK) THEN  ! get domain dimension not available via Get_subset() during process init
      sub%nc = SIZE(struct1%domain_c(:,:), 1)
      sub%nb = SIZE(struct1%domain_c(:,:), 2)
    ENDIF
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
    IF (sub%type /= ON_CHUNK) THEN  ! get domain dimension not available via Get_subset() during process init
      sub%nc = SIZE(struct1%domain_c(:,:), 1)
      sub%nb = SIZE(struct1%domain_c(:,:), 2)
    ENDIF
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
    IF (sub%type /= ON_CHUNK) THEN  ! get domain dimension not available via Get_subset() during process init
      sub%nc = SIZE(struct1%domain_c(:,:), 1)
      sub%nb = SIZE(struct1%domain_c(:,:), 2)
    ENDIF
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
    IF (sub%type /= ON_CHUNK) THEN  ! get domain dimension not available via Get_subset() during process init
      sub%nc = SIZE(struct1%domain_c(:,:), 1)
      sub%nb = SIZE(struct1%domain_c(:,:), 2)
    ENDIF
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
    IF (sub%type /= ON_CHUNK) THEN  ! get domain dimension not available via Get_subset() during process init
      sub%nc = SIZE(struct1%domain_c(:,:), 1)
      sub%nb = SIZE(struct1%domain_c(:,:), 2)
    ENDIF
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
    IF (sub%type /= ON_CHUNK) THEN  ! get domain dimension not available via Get_subset() during process init
      sub%nc = SIZE(struct1%domain_c(:,:), 1)
      sub%nb = SIZE(struct1%domain_c(:,:), 2)
    ENDIF
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
    IF (sub%type /= ON_CHUNK) THEN  ! get domain dimension not available via Get_subset() during process init
      sub%nc = SIZE(struct1%domain_c(:,:), 1)
      sub%nb = SIZE(struct1%domain_c(:,:), 2)
    ENDIF
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
    IF (sub%type /= ON_CHUNK) THEN  ! get domain dimension not available via Get_subset() during process init
      sub%nc = SIZE(struct1%domain_c(:,:), 1)
      sub%nb = SIZE(struct1%domain_c(:,:), 2)
    ENDIF
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
    IF (sub%type /= ON_CHUNK) THEN  ! get domain dimension not available via Get_subset() during process init
      sub%nc = SIZE(struct1%domain_c(:,:), 1)
      sub%nb = SIZE(struct1%domain_c(:,:), 2)
    ENDIF
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
    IF (sub%type /= ON_CHUNK) THEN  ! get domain dimension not available via Get_subset() during process init
      sub%nc = SIZE(struct1%domain_c(:,:), 1)
      sub%nb = SIZE(struct1%domain_c(:,:), 2)
    ENDIF
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
    IF (sub%type /= ON_CHUNK) THEN  ! get domain dimension not available via Get_subset() during process init
      sub%nc = SIZE(struct1%domain_c(:,:), 1)
      sub%nb = SIZE(struct1%domain_c(:,:), 2)
    ENDIF
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
    IF (sub%type /= ON_CHUNK) THEN  ! get domain dimension not available via Get_subset() during process init
      sub%nc = SIZE(struct1%domain_c(:,:), 1)
      sub%nb = SIZE(struct1%domain_c(:,:), 2)
    ENDIF
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
    IF (sub%type /= ON_CHUNK) THEN  ! get domain dimension not available via Get_subset() during process init
      sub%nc = SIZE(struct1%domain_c(:,:), 1)
      sub%nb = SIZE(struct1%domain_c(:,:), 2)
    ENDIF
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
    IF (sub%type /= ON_CHUNK) THEN  ! get domain dimension not available via Get_subset() during process init
      sub%nc = SIZE(struct1%domain_c(:,:), 1)
      sub%nb = SIZE(struct1%domain_c(:,:), 2)
    ENDIF
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
    IF (sub%type /= ON_CHUNK) THEN  ! get domain dimension not available via Get_subset() during process init
      sub%nc = SIZE(struct1%domain_c(:,:), 1)
      sub%nb = SIZE(struct1%domain_c(:,:), 2)
    ENDIF
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
    IF (sub%type /= ON_CHUNK) THEN  ! get domain dimension not available via Get_subset() during process init
      sub%nc = SIZE(struct1%domain_c(:,:), 1)
      sub%nb = SIZE(struct1%domain_c(:,:), 2)
    ENDIF
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
    IF (sub%type /= ON_CHUNK) THEN  ! get domain dimension not available via Get_subset() during process init
      sub%nc = SIZE(struct1%domain_c(:,:), 1)
      sub%nb = SIZE(struct1%domain_c(:,:), 2)
    ENDIF
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
    IF (sub%type /= ON_CHUNK) THEN  ! get domain dimension not available via Get_subset() during process init
      sub%nc = SIZE(struct1%domain_c(:,:), 1)
      sub%nb = SIZE(struct1%domain_c(:,:), 2)
    ENDIF
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
    IF (sub%type /= ON_CHUNK) THEN  ! get domain dimension not available via Get_subset() during process init
      sub%nc = SIZE(struct1%domain_c(:,:), 1)
      sub%nb = SIZE(struct1%domain_c(:,:), 2)
    ENDIF
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
    IF (sub%type /= ON_CHUNK) THEN  ! get domain dimension not available via Get_subset() during process init
      sub%nc = SIZE(struct1%domain_c(:,:), 1)
      sub%nb = SIZE(struct1%domain_c(:,:), 2)
    ENDIF
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
    IF (sub%type /= ON_CHUNK) THEN  ! get domain dimension not available via Get_subset() during process init
      sub%nc = SIZE(struct1%domain_c(:,:), 1)
      sub%nb = SIZE(struct1%domain_c(:,:), 2)
    ENDIF
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
    IF (sub%type /= ON_CHUNK) THEN  ! get domain dimension not available via Get_subset() during process init
      sub%nc = SIZE(struct1%domain_c(:,:), 1)
      sub%nb = SIZE(struct1%domain_c(:,:), 2)
    ENDIF
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
    IF (sub%type /= ON_CHUNK) THEN  ! get domain dimension not available via Get_subset() during process init
      sub%nc = SIZE(struct1%domain_c(:,:), 1)
      sub%nb = SIZE(struct1%domain_c(:,:), 2)
    ENDIF
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
    IF (sub%type /= ON_CHUNK) THEN  ! get domain dimension not available via Get_subset() during process init
      sub%nc = SIZE(struct1%domain_c(:,:), 1)
      sub%nb = SIZE(struct1%domain_c(:,:), 2)
    ENDIF
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
    IF (sub%type /= ON_CHUNK) THEN  ! get domain dimension not available via Get_subset() during process init
      sub%nc = SIZE(struct1%domain_c(:,:), 1)
      sub%nb = SIZE(struct1%domain_c(:,:), 2)
    ENDIF
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
    IF (sub%type /= ON_CHUNK) THEN  ! get domain dimension not available via Get_subset() during process init
      sub%nc = SIZE(struct1%domain_c(:,:), 1)
      sub%nb = SIZE(struct1%domain_c(:,:), 2)
    ENDIF
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
    IF (sub%type /= ON_CHUNK) THEN  ! get domain dimension not available via Get_subset() during process init
      sub%nc = SIZE(struct1%domain_c(:,:), 1)
      sub%nb = SIZE(struct1%domain_c(:,:), 2)
    ENDIF
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
    IF (sub%type /= ON_CHUNK) THEN  ! get domain dimension not available via Get_subset() during process init
      sub%nc = SIZE(struct1%domain_c(:,:), 1)
      sub%nb = SIZE(struct1%domain_c(:,:), 2)
    ENDIF
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

!! #################################################################################################################################
!! code specific to QUINCY with JSBACH4 ends here
!! #################################################################################################################################
#endif

END MODULE mo_veg_memory_class_util
