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
MODULE mo_sb_memory_class_util

!! #################################################################################################################################
!! code for QUINCY model without JSBACH4 - ends with (searchterm) "quincy-no-jsb_end"
!! #################################################################################################################################
#ifdef __QUINCY_NO_JSB__

  USE mo_kind,                   ONLY: wp
  USE mo_jsb_var_class,          ONLY: t_jsb_var_real3d, &
                                       t_jsb_struct, t_jsb_n_struct      ! quincy specific

  IMPLICIT NONE

  ! all is public except the modename; code is used by the mo_veg_memory_class only
  PUBLIC
  PRIVATE :: modname
  
  !-----------------------------------------------------------------------------------------------------
  !> Type definition for helper elements with ALLOCATABLE variables for elements in sb memory
  !! 
  !! 
  !-----------------------------------------------------------------------------------------------------
  TYPE t_sb_element_a
    REAL(wp), ALLOCATABLE :: carbon(:,:)      ! dimension of variable in memory is var(ics:ice,soil_layers,iblk)
    REAL(wp), ALLOCATABLE :: nitrogen(:,:)
    REAL(wp), ALLOCATABLE :: phosphorus(:,:)
    REAL(wp), ALLOCATABLE :: carbon13(:,:)
    REAL(wp), ALLOCATABLE :: carbon14(:,:)
    REAL(wp), ALLOCATABLE :: nitrogen15(:,:)
  CONTAINS
    PRIVATE
    !! operator overload
    ! +
    PROCEDURE               :: Add_element_a                   ! Add two elements
    PROCEDURE               :: Add_element_a_and_scalar        ! Add scalar to element
    PROCEDURE               :: Add_element_a_and_field2d       ! Add REAL field2d ( var(ics:ice,soil_layers) ) to element
    GENERIC, PUBLIC         :: OPERATOR(+) => Add_element_a, Add_element_a_and_scalar, Add_element_a_and_field2d
    !! operator overload
    ! -
    PROCEDURE               :: Subtract_element_a                   ! Subtract one element from another element
    PROCEDURE               :: Subtract_element_a_and_scalar        ! Subtract scalar from element
    PROCEDURE               :: Subtract_element_a_and_field2d       ! Subtract field2d ( var(ics:ice,soil_layers) ) from element
    GENERIC, PUBLIC         :: OPERATOR(-) => Subtract_element_a, Subtract_element_a_and_scalar, Subtract_element_a_and_field2d
    !! operator overload
    ! /
    PROCEDURE               :: Divide_element_a                   ! Divide two elements
    PROCEDURE               :: Divide_element_a_and_scalar        ! Divide scalar to element
    PROCEDURE               :: Divide_element_a_and_field2d       ! Divide REAL field2d ( var(ics:ice,soil_layers) ) to element
    GENERIC, PUBLIC         :: OPERATOR(/) => Divide_element_a, Divide_element_a_and_scalar, Divide_element_a_and_field2d
    !! operator overload
    ! *
    PROCEDURE               :: Multiply_element_a                   ! Multiply two elements
    PROCEDURE               :: Multiply_element_a_and_scalar        ! Multiply element with 
    PROCEDURE, PASS(this2)  :: Multiply_scalar_and_element_a        ! Multiply scalar with element
    PROCEDURE               :: Multiply_element_a_and_field2d       ! Multiply element with REAL field2d ( var(ics:ice,soil_layers) )
    PROCEDURE, PASS(this2)  :: Multiply_field2d_and_element_a       ! Multiply REAL field2d ( var(ics:ice,soil_layers) ) with element
    GENERIC, PUBLIC         :: OPERATOR(*) => Multiply_element_a, Multiply_element_a_and_scalar, Multiply_scalar_and_element_a, &
                                              Multiply_element_a_and_field2d, Multiply_field2d_and_element_a
  END TYPE t_sb_element_a

  !-----------------------------------------------------------------------------------------------------
  !> Type definition for helper pools with ALLOCATABLE variables for pools in sb memory
  !! 
  !! 
  !-----------------------------------------------------------------------------------------------------
  TYPE t_sb_pool_a
    TYPE(t_sb_element_a)  :: dom                 !< dissolved organic matter
    TYPE(t_sb_element_a)  :: dom_assoc           !< minerally associated disolved organic matter
    TYPE(t_sb_element_a)  :: soluable_litter     !< metabolic litter (sugars, cellulose etc.)
    TYPE(t_sb_element_a)  :: polymeric_litter    !< structural litter (lignified litter)
    TYPE(t_sb_element_a)  :: woody_litter        !< physically protected wood litter
    TYPE(t_sb_element_a)  :: fungi               !< saprophytic fungal biomass
    TYPE(t_sb_element_a)  :: mycorrhiza          !< mycorrhizal fungal biomass 
    TYPE(t_sb_element_a)  :: microbial           !< microbial biomass (bacteria!)
    TYPE(t_sb_element_a)  :: residue             !< necromass of fungi, microbial and mycorrhiza
    TYPE(t_sb_element_a)  :: residue_assoc       !< minerally associated microbial necromass 
  CONTAINS
    PRIVATE
    !! operator overload
    ! +
    PROCEDURE               :: Add_pool_a                   ! Add two pools
    PROCEDURE               :: Add_pool_a_and_scalar        ! Add scalar to pool
    PROCEDURE               :: Add_pool_a_and_field2d       ! Add REAL field2d ( var(ics:ice,soil_layers) ) to pool
    GENERIC, PUBLIC         :: OPERATOR(+) => Add_pool_a, Add_pool_a_and_scalar, Add_pool_a_and_field2d
    !! operator overload
    ! -
    PROCEDURE               :: Subtract_pool_a                   ! Subtract one pool from another pool
    PROCEDURE               :: Subtract_pool_a_and_scalar        ! Subtract scalar from pool
    PROCEDURE               :: Subtract_pool_a_and_field2d       ! Subtract field2d ( var(ics:ice,soil_layers) ) from pool
    GENERIC, PUBLIC         :: OPERATOR(-) => Subtract_pool_a, Subtract_pool_a_and_scalar, Subtract_pool_a_and_field2d
    !! operator overload
    ! /
    PROCEDURE               :: Divide_pool_a                   ! Divide two pools
    PROCEDURE               :: Divide_pool_a_and_scalar        ! Divide scalar to pool
    PROCEDURE               :: Divide_pool_a_and_field2d       ! Divide REAL field2d ( var(ics:ice,soil_layers) ) to pool
    GENERIC, PUBLIC         :: OPERATOR(/) => Divide_pool_a, Divide_pool_a_and_scalar, Divide_pool_a_and_field2d
    !! operator overload
    ! *
    PROCEDURE               :: Multiply_pool_a                   ! Multiply two pools
    PROCEDURE               :: Multiply_pool_a_and_scalar        ! Multiply pool with scalar
    PROCEDURE, PASS(this2)  :: Multiply_scalar_and_pool_a        ! Multiply scalar with pool
    PROCEDURE               :: Multiply_pool_a_and_field2d       ! Multiply pool with REAL 1D field2d ( var(ics:ice) )
    PROCEDURE, PASS(this2)  :: Multiply_field2d_and_pool_a       ! Multiply REAL 1D field2d ( var(ics:ice,soil_layers) ) with pool
    GENERIC, PUBLIC         :: OPERATOR(*) => Multiply_pool_a, Multiply_pool_a_and_scalar, Multiply_scalar_and_pool_a, &
                                              Multiply_pool_a_and_field2d, Multiply_field2d_and_pool_a
  END TYPE t_sb_pool_a

  !-----------------------------------------------------------------------------------------------------
  !> Type definition for elements in sb pool
  !! 
  !! elements in soil are 3D 
  !! 
  !-----------------------------------------------------------------------------------------------------
  TYPE, EXTENDS(t_jsb_struct) :: t_quincy_struct_sb_element
    TYPE(t_jsb_var_real3d)    :: carbon       !< total carbon concentration of pool [mol m-3]
    TYPE(t_jsb_var_real3d)    :: nitrogen     !< totol nitrogen concentration of pool [mol m-3]
    TYPE(t_jsb_var_real3d)    :: phosphorus   !< phosphorus concentration of pool [mol m-3]
    TYPE(t_jsb_var_real3d)    :: carbon13     !< 13carbon concentration of pool [mol m-3]
    TYPE(t_jsb_var_real3d)    :: carbon14     !< 14carbon concentration of pool [pmol m-3]
    TYPE(t_jsb_var_real3d)    :: nitrogen15   !< 15nitrogen concentration of pool [mol m-3]
  CONTAINS
    PROCEDURE :: Init       => Init_t_quincy_struct_sb_element
    PROCEDURE :: Add_vars   => Add_element_variables_to_memory
  END TYPE t_quincy_struct_sb_element

  !-----------------------------------------------------------------------------------------------------
  !> Type definition for sb pool with nested TYPEs
  !!
  !! used for: pool ...
  !!  
  !-----------------------------------------------------------------------------------------------------
  TYPE, EXTENDS(t_jsb_n_struct) :: t_quincy_struct_sb_pool
    TYPE(t_quincy_struct_sb_element)  :: dom                 !< dissolved organic matter
    TYPE(t_quincy_struct_sb_element)  :: dom_assoc           !< minerally associated disolved organic matter
    TYPE(t_quincy_struct_sb_element)  :: soluable_litter     !< metabolic litter (sugars, cellulose etc.)
    TYPE(t_quincy_struct_sb_element)  :: polymeric_litter    !< structural litter (lignified litter)
    TYPE(t_quincy_struct_sb_element)  :: woody_litter        !< physically protected wood litter
    TYPE(t_quincy_struct_sb_element)  :: fungi               !< saprophytic fungal biomass
    TYPE(t_quincy_struct_sb_element)  :: mycorrhiza          !< mycorrhizal fungal biomass 
    TYPE(t_quincy_struct_sb_element)  :: microbial           !< microbial biomass (bacteria!)
    TYPE(t_quincy_struct_sb_element)  :: residue             !< necromass of fungi, microbial and mycorrhiza
    TYPE(t_quincy_struct_sb_element)  :: residue_assoc       !< minerally associated microbial necromass 
  CONTAINS
    PROCEDURE :: Init           => Init_t_quincy_struct_sb_pool
    PROCEDURE :: Add_elements   => Add_pool_elements_to_memory
  END TYPE t_quincy_struct_sb_pool

  ! definition of sb_ pool and element structures with ALLOCATABLE variables
  ! allocation of variables done with Init_sb_pools_elements_a during model init in mo_jsb_model_init:quincy_init
  TYPE(t_sb_pool_a)       :: sb_pool             , &
                             sb_formation        , &
                             sb_loss             , &
                             sb_transport
  TYPE(t_sb_element_a)    :: sb_delta_pools      , &
                             sb_mycorrhiza_export

  CHARACTER(len=*), PARAMETER :: modname = 'mo_sb_memory_class_util'

CONTAINS

  !-----------------------------------------------------------------------------------------------------!
  ! t_sb_element_a operator overload                                                                   !
  !-----------------------------------------------------------------------------------------------------!
  !-----------------------------------------------------------------------------------------------------
  !> + (t_sb_element_a + t_sb_element_a)
  !!
  !-----------------------------------------------------------------------------------------------------
  PURE FUNCTION Add_element_a(this,element2) RESULT(element3)
    IMPLICIT NONE
    ! 0.1 InOut
    CLASS(t_sb_element_a), INTENT(in)  :: this, element2
    TYPE(t_sb_element_a)               :: element3
    ! 1.0 go
    ALLOCATE(element3%carbon,     source=this%carbon)       ! note: Function results are automatically deallocated after the result is "used" in the outer scope.
    ALLOCATE(element3%nitrogen,   source=this%nitrogen)
    ALLOCATE(element3%phosphorus, source=this%phosphorus)
    ALLOCATE(element3%carbon13,   source=this%carbon13)
    ALLOCATE(element3%carbon14,   source=this%carbon14)
    ALLOCATE(element3%nitrogen15, source=this%nitrogen15)
    element3%carbon(:,:)       = this%carbon(:,:)      + element2%carbon(:,:)
    element3%nitrogen(:,:)     = this%nitrogen(:,:)    + element2%nitrogen(:,:)
    element3%phosphorus(:,:)   = this%phosphorus(:,:)  + element2%phosphorus(:,:)
    element3%carbon13(:,:)     = this%carbon13(:,:)    + element2%carbon13(:,:)
    element3%carbon14(:,:)     = this%carbon14(:,:)    + element2%carbon14(:,:)
    element3%nitrogen15(:,:)   = this%nitrogen15(:,:)  + element2%nitrogen15(:,:)
  END FUNCTION Add_element_a
  !-----------------------------------------------------------------------------------------------------
  !> + (t_sb_element_a + scalar)
  !!
  !-----------------------------------------------------------------------------------------------------
  PURE FUNCTION Add_element_a_and_scalar(this,scalar) RESULT(element3)
    IMPLICIT NONE
    ! 0.1 InOut
    CLASS(t_sb_element_a), INTENT(in)  :: this
    REAL(wp),              INTENT(in)  :: scalar
    TYPE(t_sb_element_a)               :: element3
    ! 1.0 go
    ALLOCATE(element3%carbon,     source=this%carbon)
    ALLOCATE(element3%nitrogen,   source=this%nitrogen)
    ALLOCATE(element3%phosphorus, source=this%phosphorus)
    ALLOCATE(element3%carbon13,   source=this%carbon13)
    ALLOCATE(element3%carbon14,   source=this%carbon14)
    ALLOCATE(element3%nitrogen15, source=this%nitrogen15)
    element3%carbon(:,:)       = this%carbon(:,:)      + scalar
    element3%nitrogen(:,:)     = this%nitrogen(:,:)    + scalar
    element3%phosphorus(:,:)   = this%phosphorus(:,:)  + scalar
    element3%carbon13(:,:)     = this%carbon13(:,:)    + scalar
    element3%carbon14(:,:)     = this%carbon14(:,:)    + scalar
    element3%nitrogen15(:,:)   = this%nitrogen15(:,:)  + scalar
  END FUNCTION Add_element_a_and_scalar
  !-----------------------------------------------------------------------------------------------------
  !> + (t_sb_element_a + field2d)
  !!
  !-----------------------------------------------------------------------------------------------------
  PURE FUNCTION Add_element_a_and_field2d(this,field2d) RESULT(element3)
    IMPLICIT NONE
    ! 0.1 InOut
    CLASS(t_sb_element_a), INTENT(in)  :: this
    REAL(wp),              INTENT(in)  :: field2d(:,:)
    TYPE(t_sb_element_a)               :: element3
    ! 1.0 go
    ALLOCATE(element3%carbon,     source=this%carbon)
    ALLOCATE(element3%nitrogen,   source=this%nitrogen)
    ALLOCATE(element3%phosphorus, source=this%phosphorus)
    ALLOCATE(element3%carbon13,   source=this%carbon13)
    ALLOCATE(element3%carbon14,   source=this%carbon14)
    ALLOCATE(element3%nitrogen15, source=this%nitrogen15)
    element3%carbon(:,:)       = this%carbon(:,:)      + field2d(:,:)
    element3%nitrogen(:,:)     = this%nitrogen(:,:)    + field2d(:,:)
    element3%phosphorus(:,:)   = this%phosphorus(:,:)  + field2d(:,:)
    element3%carbon13(:,:)     = this%carbon13(:,:)    + field2d(:,:)
    element3%carbon14(:,:)     = this%carbon14(:,:)    + field2d(:,:)
    element3%nitrogen15(:,:)   = this%nitrogen15(:,:)  + field2d(:,:)
  END FUNCTION Add_element_a_and_field2d
  !-----------------------------------------------------------------------------------------------------
  !> - (t_sb_element_a - t_sb_element_a)
  !!
  !-----------------------------------------------------------------------------------------------------
  PURE FUNCTION Subtract_element_a(this,element2) RESULT(element3)
    IMPLICIT NONE
    ! 0.1 InOut
    CLASS(t_sb_element_a), INTENT(in)  :: this, element2
    TYPE(t_sb_element_a)               :: element3
    ! 1.0 go
    ALLOCATE(element3%carbon,     source=this%carbon)
    ALLOCATE(element3%nitrogen,   source=this%nitrogen)
    ALLOCATE(element3%phosphorus, source=this%phosphorus)
    ALLOCATE(element3%carbon13,   source=this%carbon13)
    ALLOCATE(element3%carbon14,   source=this%carbon14)
    ALLOCATE(element3%nitrogen15, source=this%nitrogen15)
    element3%carbon(:,:)       = this%carbon(:,:)      - element2%carbon(:,:)
    element3%nitrogen(:,:)     = this%nitrogen(:,:)    - element2%nitrogen(:,:)
    element3%phosphorus(:,:)   = this%phosphorus(:,:)  - element2%phosphorus(:,:)
    element3%carbon13(:,:)     = this%carbon13(:,:)    - element2%carbon13(:,:)
    element3%carbon14(:,:)     = this%carbon14(:,:)    - element2%carbon14(:,:)
    element3%nitrogen15(:,:)   = this%nitrogen15(:,:)  - element2%nitrogen15(:,:)
  END FUNCTION Subtract_element_a
  !-----------------------------------------------------------------------------------------------------
  !> - (t_sb_element_a - scalar)
  !!
  !-----------------------------------------------------------------------------------------------------
  PURE FUNCTION Subtract_element_a_and_scalar(this,scalar) RESULT(element3)
    IMPLICIT NONE
    ! 0.1 InOut
    CLASS(t_sb_element_a), INTENT(in)  :: this
    REAL(wp),              INTENT(in)  :: scalar
    TYPE(t_sb_element_a)               :: element3
    ! 1.0 go
    ALLOCATE(element3%carbon,     source=this%carbon)
    ALLOCATE(element3%nitrogen,   source=this%nitrogen)
    ALLOCATE(element3%phosphorus, source=this%phosphorus)
    ALLOCATE(element3%carbon13,   source=this%carbon13)
    ALLOCATE(element3%carbon14,   source=this%carbon14)
    ALLOCATE(element3%nitrogen15, source=this%nitrogen15)
    element3%carbon(:,:)       = this%carbon(:,:)      - scalar
    element3%nitrogen(:,:)     = this%nitrogen(:,:)    - scalar
    element3%phosphorus(:,:)   = this%phosphorus(:,:)  - scalar
    element3%carbon13(:,:)     = this%carbon13(:,:)    - scalar
    element3%carbon14(:,:)     = this%carbon14(:,:)    - scalar
    element3%nitrogen15(:,:)   = this%nitrogen15(:,:)  - scalar
  END FUNCTION Subtract_element_a_and_scalar
  !-----------------------------------------------------------------------------------------------------
  !> - (t_sb_element_a - field2d)
  !!
  !-----------------------------------------------------------------------------------------------------
  PURE FUNCTION Subtract_element_a_and_field2d(this,field2d) RESULT(element3)
    IMPLICIT NONE
    ! 0.1 InOut
    CLASS(t_sb_element_a), INTENT(in)  :: this
    REAL(wp),              INTENT(in)  :: field2d(:,:)
    TYPE(t_sb_element_a)               :: element3
    ! 1.0 go
    ALLOCATE(element3%carbon,     source=this%carbon)
    ALLOCATE(element3%nitrogen,   source=this%nitrogen)
    ALLOCATE(element3%phosphorus, source=this%phosphorus)
    ALLOCATE(element3%carbon13,   source=this%carbon13)
    ALLOCATE(element3%carbon14,   source=this%carbon14)
    ALLOCATE(element3%nitrogen15, source=this%nitrogen15)
    element3%carbon(:,:)       = this%carbon(:,:)      - field2d(:,:)
    element3%nitrogen(:,:)     = this%nitrogen(:,:)    - field2d(:,:)
    element3%phosphorus(:,:)   = this%phosphorus(:,:)  - field2d(:,:)
    element3%carbon13(:,:)     = this%carbon13(:,:)    - field2d(:,:)
    element3%carbon14(:,:)     = this%carbon14(:,:)    - field2d(:,:)
    element3%nitrogen15(:,:)   = this%nitrogen15(:,:)  - field2d(:,:)
  END FUNCTION Subtract_element_a_and_field2d
  !-----------------------------------------------------------------------------------------------------
  !> / (t_sb_element_a / t_sb_element_a)
  !!
  !-----------------------------------------------------------------------------------------------------
  PURE FUNCTION Divide_element_a(this,element2) RESULT(element3)
    IMPLICIT NONE
    ! 0.1 InOut
    CLASS(t_sb_element_a), INTENT(in)  :: this, element2
    TYPE(t_sb_element_a)               :: element3
    ! 1.0 go
    ALLOCATE(element3%carbon,     source=this%carbon)       ! note: Function results are automatically deallocated after the result is "used" in the outer scope.
    ALLOCATE(element3%nitrogen,   source=this%nitrogen)
    ALLOCATE(element3%phosphorus, source=this%phosphorus)
    ALLOCATE(element3%carbon13,   source=this%carbon13)
    ALLOCATE(element3%carbon14,   source=this%carbon14)
    ALLOCATE(element3%nitrogen15, source=this%nitrogen15)
    element3%carbon(:,:)       = this%carbon(:,:)      / element2%carbon(:,:)
    element3%nitrogen(:,:)     = this%nitrogen(:,:)    / element2%nitrogen(:,:)
    element3%phosphorus(:,:)   = this%phosphorus(:,:)  / element2%phosphorus(:,:)
    element3%carbon13(:,:)     = this%carbon13(:,:)    / element2%carbon13(:,:)
    element3%carbon14(:,:)     = this%carbon14(:,:)    / element2%carbon14(:,:)
    element3%nitrogen15(:,:)   = this%nitrogen15(:,:)  / element2%nitrogen15(:,:)
  END FUNCTION Divide_element_a
  !-----------------------------------------------------------------------------------------------------
  !> / (t_sb_element_a / scalar)
  !!
  !-----------------------------------------------------------------------------------------------------
  PURE FUNCTION Divide_element_a_and_scalar(this,scalar) RESULT(element3)
    IMPLICIT NONE
    ! 0.1 InOut
    CLASS(t_sb_element_a), INTENT(in)  :: this
    REAL(wp),              INTENT(in)  :: scalar
    TYPE(t_sb_element_a)               :: element3
    ! 1.0 go
    ALLOCATE(element3%carbon,     source=this%carbon)
    ALLOCATE(element3%nitrogen,   source=this%nitrogen)
    ALLOCATE(element3%phosphorus, source=this%phosphorus)
    ALLOCATE(element3%carbon13,   source=this%carbon13)
    ALLOCATE(element3%carbon14,   source=this%carbon14)
    ALLOCATE(element3%nitrogen15, source=this%nitrogen15)
    element3%carbon(:,:)       = this%carbon(:,:)      / scalar
    element3%nitrogen(:,:)     = this%nitrogen(:,:)    / scalar
    element3%phosphorus(:,:)   = this%phosphorus(:,:)  / scalar
    element3%carbon13(:,:)     = this%carbon13(:,:)    / scalar
    element3%carbon14(:,:)     = this%carbon14(:,:)    / scalar
    element3%nitrogen15(:,:)   = this%nitrogen15(:,:)  / scalar
  END FUNCTION Divide_element_a_and_scalar
  !-----------------------------------------------------------------------------------------------------
  !> / (t_sb_element_a / field2d)
  !!
  !-----------------------------------------------------------------------------------------------------
  PURE FUNCTION Divide_element_a_and_field2d(this,field2d) RESULT(element3)
    IMPLICIT NONE
    ! 0.1 InOut
    CLASS(t_sb_element_a), INTENT(in)  :: this
    REAL(wp),              INTENT(in)  :: field2d(:,:)
    TYPE(t_sb_element_a)               :: element3
    ! 1.0 go
    ALLOCATE(element3%carbon,     source=this%carbon)
    ALLOCATE(element3%nitrogen,   source=this%nitrogen)
    ALLOCATE(element3%phosphorus, source=this%phosphorus)
    ALLOCATE(element3%carbon13,   source=this%carbon13)
    ALLOCATE(element3%carbon14,   source=this%carbon14)
    ALLOCATE(element3%nitrogen15, source=this%nitrogen15)
    element3%carbon(:,:)       = this%carbon(:,:)      / field2d(:,:)
    element3%nitrogen(:,:)     = this%nitrogen(:,:)    / field2d(:,:)
    element3%phosphorus(:,:)   = this%phosphorus(:,:)  / field2d(:,:)
    element3%carbon13(:,:)     = this%carbon13(:,:)    / field2d(:,:)
    element3%carbon14(:,:)     = this%carbon14(:,:)    / field2d(:,:)
    element3%nitrogen15(:,:)   = this%nitrogen15(:,:)  / field2d(:,:)
  END FUNCTION Divide_element_a_and_field2d
  !-----------------------------------------------------------------------------------------------------
  !> * (t_sb_element_a * t_sb_element_a)
  !!
  !-----------------------------------------------------------------------------------------------------
  PURE FUNCTION Multiply_element_a(this,element2) RESULT(element3)
    IMPLICIT NONE
    ! 0.1 InOut
    CLASS(t_sb_element_a), INTENT(in)  :: this, element2
    TYPE(t_sb_element_a)               :: element3
    ! 1.0 go
    ALLOCATE(element3%carbon,     source=this%carbon)
    ALLOCATE(element3%nitrogen,   source=this%nitrogen)
    ALLOCATE(element3%phosphorus, source=this%phosphorus)
    ALLOCATE(element3%carbon13,   source=this%carbon13)
    ALLOCATE(element3%carbon14,   source=this%carbon14)
    ALLOCATE(element3%nitrogen15, source=this%nitrogen15)
    element3%carbon(:,:)       = this%carbon(:,:)      * element2%carbon(:,:)
    element3%nitrogen(:,:)     = this%nitrogen(:,:)    * element2%nitrogen(:,:)
    element3%phosphorus(:,:)   = this%phosphorus(:,:)  * element2%phosphorus(:,:)
    element3%carbon13(:,:)     = this%carbon13(:,:)    * element2%carbon13(:,:)
    element3%carbon14(:,:)     = this%carbon14(:,:)    * element2%carbon14(:,:)
    element3%nitrogen15(:,:)   = this%nitrogen15(:,:)  * element2%nitrogen15(:,:)
  END FUNCTION Multiply_element_a
  !-----------------------------------------------------------------------------------------------------
  !> * (t_sb_element_a * scalar)
  !!
  !-----------------------------------------------------------------------------------------------------
  PURE FUNCTION Multiply_element_a_and_scalar(this,scalar) RESULT(element3)
    IMPLICIT NONE
    ! 0.1 InOut
    CLASS(t_sb_element_a), INTENT(in)  :: this
    REAL(wp),              INTENT(in)  :: scalar
    TYPE(t_sb_element_a)               :: element3
    ! 1.0 go
    ALLOCATE(element3%carbon,     source=this%carbon)
    ALLOCATE(element3%nitrogen,   source=this%nitrogen)
    ALLOCATE(element3%phosphorus, source=this%phosphorus)
    ALLOCATE(element3%carbon13,   source=this%carbon13)
    ALLOCATE(element3%carbon14,   source=this%carbon14)
    ALLOCATE(element3%nitrogen15, source=this%nitrogen15)
    element3%carbon(:,:)       = this%carbon(:,:)      * scalar
    element3%nitrogen(:,:)     = this%nitrogen(:,:)    * scalar
    element3%phosphorus(:,:)   = this%phosphorus(:,:)  * scalar
    element3%carbon13(:,:)     = this%carbon13(:,:)    * scalar
    element3%carbon14(:,:)     = this%carbon14(:,:)    * scalar
    element3%nitrogen15(:,:)   = this%nitrogen15(:,:)  * scalar
  END FUNCTION Multiply_element_a_and_scalar
  !-----------------------------------------------------------------------------------------------------
  !> * (scalar * t_sb_element_a)
  !!
  !-----------------------------------------------------------------------------------------------------
  PURE FUNCTION Multiply_scalar_and_element_a(scalar,this2) RESULT(element3)
    IMPLICIT NONE
    ! 0.1 InOut
    REAL(wp),              INTENT(in)  :: scalar
    CLASS(t_sb_element_a), INTENT(in)  :: this2
    TYPE(t_sb_element_a)               :: element3
    ! 1.0 go
    ALLOCATE(element3%carbon,     source=this2%carbon)
    ALLOCATE(element3%nitrogen,   source=this2%nitrogen)
    ALLOCATE(element3%phosphorus, source=this2%phosphorus)
    ALLOCATE(element3%carbon13,   source=this2%carbon13)
    ALLOCATE(element3%carbon14,   source=this2%carbon14)
    ALLOCATE(element3%nitrogen15, source=this2%nitrogen15)
    element3%carbon(:,:)       = scalar * this2%carbon(:,:)
    element3%nitrogen(:,:)     = scalar * this2%nitrogen(:,:)
    element3%phosphorus(:,:)   = scalar * this2%phosphorus(:,:)
    element3%carbon13(:,:)     = scalar * this2%carbon13(:,:)
    element3%carbon14(:,:)     = scalar * this2%carbon14(:,:)
    element3%nitrogen15(:,:)   = scalar * this2%nitrogen15(:,:)
  END FUNCTION Multiply_scalar_and_element_a
  !-----------------------------------------------------------------------------------------------------
  !> * (t_sb_element_a * field2d)
  !!
  !-----------------------------------------------------------------------------------------------------
  PURE FUNCTION Multiply_element_a_and_field2d(this,field2d) RESULT(element3)
    IMPLICIT NONE
    ! 0.1 InOut
    CLASS(t_sb_element_a), INTENT(in)  :: this
    REAL(wp),              INTENT(in)  :: field2d(:,:)
    TYPE(t_sb_element_a)               :: element3
    ! 1.0 go
    ALLOCATE(element3%carbon,     source=this%carbon)
    ALLOCATE(element3%nitrogen,   source=this%nitrogen)
    ALLOCATE(element3%phosphorus, source=this%phosphorus)
    ALLOCATE(element3%carbon13,   source=this%carbon13)
    ALLOCATE(element3%carbon14,   source=this%carbon14)
    ALLOCATE(element3%nitrogen15, source=this%nitrogen15)
    element3%carbon(:,:)       = this%carbon(:,:)      * field2d(:,:)
    element3%nitrogen(:,:)     = this%nitrogen(:,:)    * field2d(:,:)
    element3%phosphorus(:,:)   = this%phosphorus(:,:)  * field2d(:,:)
    element3%carbon13(:,:)     = this%carbon13(:,:)    * field2d(:,:)
    element3%carbon14(:,:)     = this%carbon14(:,:)    * field2d(:,:)
    element3%nitrogen15(:,:)   = this%nitrogen15(:,:)  * field2d(:,:)
  END FUNCTION Multiply_element_a_and_field2d
  !-----------------------------------------------------------------------------------------------------
  !> * (field2d * t_sb_element_a)
  !!
  !-----------------------------------------------------------------------------------------------------
  PURE FUNCTION Multiply_field2d_and_element_a(field2d,this2) RESULT(element3)
    IMPLICIT NONE
    ! 0.1 InOut
    REAL(wp),              INTENT(in)  :: field2d(:,:)
    CLASS(t_sb_element_a), INTENT(in)  :: this2
    TYPE(t_sb_element_a)               :: element3
    ! 1.0 go
    ALLOCATE(element3%carbon,     source=this2%carbon)
    ALLOCATE(element3%nitrogen,   source=this2%nitrogen)
    ALLOCATE(element3%phosphorus, source=this2%phosphorus)
    ALLOCATE(element3%carbon13,   source=this2%carbon13)
    ALLOCATE(element3%carbon14,   source=this2%carbon14)
    ALLOCATE(element3%nitrogen15, source=this2%nitrogen15)
    element3%carbon(:,:)       = field2d(:,:) * this2%carbon(:,:)
    element3%nitrogen(:,:)     = field2d(:,:) * this2%nitrogen(:,:)
    element3%phosphorus(:,:)   = field2d(:,:) * this2%phosphorus(:,:)
    element3%carbon13(:,:)     = field2d(:,:) * this2%carbon13(:,:)
    element3%carbon14(:,:)     = field2d(:,:) * this2%carbon14(:,:)
    element3%nitrogen15(:,:)   = field2d(:,:) * this2%nitrogen15(:,:)
  END FUNCTION Multiply_field2d_and_element_a

  !-----------------------------------------------------------------------------------------------------!
  ! t_sb_pool_a operator overload                                                                       !
  !-----------------------------------------------------------------------------------------------------!
  !-----------------------------------------------------------------------------------------------------
  !> + (t_sb_pool_a + t_sb_pool_a)
  !!
  !-----------------------------------------------------------------------------------------------------
  PURE FUNCTION Add_pool_a(this,pool2) RESULT(pool3)
    IMPLICIT NONE
    ! 0.1 InOut
    CLASS(t_sb_pool_a), INTENT(in) :: this, pool2
    TYPE(t_sb_pool_a)              :: pool3
    ! 1.0 go
    pool3%dom               = this%dom              + pool2%dom
    pool3%dom_assoc         = this%dom_assoc        + pool2%dom_assoc
    pool3%soluable_litter   = this%soluable_litter  + pool2%soluable_litter
    pool3%polymeric_litter  = this%polymeric_litter + pool2%polymeric_litter
    pool3%woody_litter      = this%woody_litter     + pool2%woody_litter
    pool3%fungi             = this%fungi            + pool2%fungi
    pool3%mycorrhiza        = this%mycorrhiza       + pool2%mycorrhiza
    pool3%microbial         = this%microbial        + pool2%microbial
    pool3%residue           = this%residue          + pool2%residue
    pool3%residue_assoc     = this%residue_assoc    + pool2%residue_assoc
  END FUNCTION Add_pool_a
  !-----------------------------------------------------------------------------------------------------
  !> + (t_sb_pool_a + scalar)
  !!
  !-----------------------------------------------------------------------------------------------------
  PURE FUNCTION Add_pool_a_and_scalar(this,scalar) RESULT(pool3)
    IMPLICIT NONE
    ! 0.1 InOut
    CLASS(t_sb_pool_a), INTENT(in) :: this
    REAL(wp),           INTENT(in) :: scalar
    TYPE(t_sb_pool_a)              :: pool3
    ! 1.0 go
    pool3%dom               = this%dom              + scalar
    pool3%dom_assoc         = this%dom_assoc        + scalar
    pool3%soluable_litter   = this%soluable_litter  + scalar
    pool3%polymeric_litter  = this%polymeric_litter + scalar
    pool3%woody_litter      = this%woody_litter     + scalar
    pool3%fungi             = this%fungi            + scalar
    pool3%mycorrhiza        = this%mycorrhiza       + scalar
    pool3%microbial         = this%microbial        + scalar
    pool3%residue           = this%residue          + scalar
    pool3%residue_assoc     = this%residue_assoc    + scalar
  END FUNCTION Add_pool_a_and_scalar
  !-----------------------------------------------------------------------------------------------------
  !> + (t_sb_pool_a + field2d)
  !!
  !-----------------------------------------------------------------------------------------------------
  PURE FUNCTION Add_pool_a_and_field2d(this,field2d) RESULT(pool3)
    IMPLICIT NONE
    ! 0.1 InOut
    CLASS(t_sb_pool_a), INTENT(in) :: this
    REAL(wp),           INTENT(in) :: field2d(:,:)
    TYPE(t_sb_pool_a)              :: pool3
    ! 1.0 go
    pool3%dom               = this%dom              + field2d(:,:)
    pool3%dom_assoc         = this%dom_assoc        + field2d(:,:)
    pool3%soluable_litter   = this%soluable_litter  + field2d(:,:)
    pool3%polymeric_litter  = this%polymeric_litter + field2d(:,:)
    pool3%woody_litter      = this%woody_litter     + field2d(:,:)
    pool3%fungi             = this%fungi            + field2d(:,:)
    pool3%mycorrhiza        = this%mycorrhiza       + field2d(:,:)
    pool3%microbial         = this%microbial        + field2d(:,:)
    pool3%residue           = this%residue          + field2d(:,:)
    pool3%residue_assoc     = this%residue_assoc    + field2d(:,:)
  END FUNCTION Add_pool_a_and_field2d
  !-----------------------------------------------------------------------------------------------------
  !> - (t_sb_pool_a - t_sb_pool_a)
  !!
  !-----------------------------------------------------------------------------------------------------
  PURE FUNCTION Subtract_pool_a(this,pool2) RESULT(pool3)
    IMPLICIT NONE
    ! 0.1 InOut
    CLASS(t_sb_pool_a), INTENT(in) :: this, pool2
    TYPE(t_sb_pool_a)              :: pool3
    ! 1.0 go
    pool3%dom               = this%dom              - pool2%dom
    pool3%dom_assoc         = this%dom_assoc        - pool2%dom_assoc
    pool3%soluable_litter   = this%soluable_litter  - pool2%soluable_litter
    pool3%polymeric_litter  = this%polymeric_litter - pool2%polymeric_litter
    pool3%woody_litter      = this%woody_litter     - pool2%woody_litter
    pool3%fungi             = this%fungi            - pool2%fungi
    pool3%mycorrhiza        = this%mycorrhiza       - pool2%mycorrhiza
    pool3%microbial         = this%microbial        - pool2%microbial
    pool3%residue           = this%residue          - pool2%residue
    pool3%residue_assoc     = this%residue_assoc    - pool2%residue_assoc
  END FUNCTION Subtract_pool_a
  !-----------------------------------------------------------------------------------------------------
  !> - (t_sb_pool_a - scalar)
  !!
  !-----------------------------------------------------------------------------------------------------
  PURE FUNCTION Subtract_pool_a_and_scalar(this,scalar) RESULT(pool3)
    IMPLICIT NONE
    ! 0.1 InOut
    CLASS(t_sb_pool_a), INTENT(in) :: this
    REAL(wp),           INTENT(in) :: scalar
    TYPE(t_sb_pool_a)              :: pool3
    ! 1.0 go
    pool3%dom               = this%dom              - scalar
    pool3%dom_assoc         = this%dom_assoc        - scalar
    pool3%soluable_litter   = this%soluable_litter  - scalar
    pool3%polymeric_litter  = this%polymeric_litter - scalar
    pool3%woody_litter      = this%woody_litter     - scalar
    pool3%fungi             = this%fungi            - scalar
    pool3%mycorrhiza        = this%mycorrhiza       - scalar
    pool3%microbial         = this%microbial        - scalar
    pool3%residue           = this%residue          - scalar
    pool3%residue_assoc     = this%residue_assoc    - scalar
  END FUNCTION Subtract_pool_a_and_scalar
  !-----------------------------------------------------------------------------------------------------
  !> - (t_sb_pool_a - field2d)
  !!
  !-----------------------------------------------------------------------------------------------------
  PURE FUNCTION Subtract_pool_a_and_field2d(this,field2d) RESULT(pool3)
    IMPLICIT NONE
    ! 0.1 InOut
    CLASS(t_sb_pool_a), INTENT(in) :: this
    REAL(wp),           INTENT(in) :: field2d(:,:)
    TYPE(t_sb_pool_a)              :: pool3
    ! 1.0 go
    pool3%dom               = this%dom              - field2d(:,:)
    pool3%dom_assoc         = this%dom_assoc        - field2d(:,:)
    pool3%soluable_litter   = this%soluable_litter  - field2d(:,:)
    pool3%polymeric_litter  = this%polymeric_litter - field2d(:,:)
    pool3%woody_litter      = this%woody_litter     - field2d(:,:)
    pool3%fungi             = this%fungi            - field2d(:,:)
    pool3%mycorrhiza        = this%mycorrhiza       - field2d(:,:)
    pool3%microbial         = this%microbial        - field2d(:,:)
    pool3%residue           = this%residue          - field2d(:,:)
    pool3%residue_assoc     = this%residue_assoc    - field2d(:,:)
  END FUNCTION Subtract_pool_a_and_field2d
  !-----------------------------------------------------------------------------------------------------
  !> / (t_sb_pool_a / t_sb_pool_a)
  !!
  !-----------------------------------------------------------------------------------------------------
  PURE FUNCTION Divide_pool_a(this,pool2) RESULT(pool3)
    IMPLICIT NONE
    ! 0.1 InOut
    CLASS(t_sb_pool_a), INTENT(in) :: this, pool2
    TYPE(t_sb_pool_a)              :: pool3
    ! 1.0 go
    pool3%dom               = this%dom              / pool2%dom
    pool3%dom_assoc         = this%dom_assoc        / pool2%dom_assoc
    pool3%soluable_litter   = this%soluable_litter  / pool2%soluable_litter
    pool3%polymeric_litter  = this%polymeric_litter / pool2%polymeric_litter
    pool3%woody_litter      = this%woody_litter     / pool2%woody_litter
    pool3%fungi             = this%fungi            / pool2%fungi
    pool3%mycorrhiza        = this%mycorrhiza       / pool2%mycorrhiza
    pool3%microbial         = this%microbial        / pool2%microbial
    pool3%residue           = this%residue          / pool2%residue
    pool3%residue_assoc     = this%residue_assoc    / pool2%residue_assoc
  END FUNCTION Divide_pool_a
  !-----------------------------------------------------------------------------------------------------
  !> / (t_sb_pool_a / scalar)
  !!
  !-----------------------------------------------------------------------------------------------------
  PURE FUNCTION Divide_pool_a_and_scalar(this,scalar) RESULT(pool3)
    IMPLICIT NONE
    ! 0.1 InOut
    CLASS(t_sb_pool_a), INTENT(in) :: this
    REAL(wp),           INTENT(in) :: scalar
    TYPE(t_sb_pool_a)              :: pool3
    ! 1.0 go
    pool3%dom               = this%dom              / scalar
    pool3%dom_assoc         = this%dom_assoc        / scalar
    pool3%soluable_litter   = this%soluable_litter  / scalar
    pool3%polymeric_litter  = this%polymeric_litter / scalar
    pool3%woody_litter      = this%woody_litter     / scalar
    pool3%fungi             = this%fungi            / scalar
    pool3%mycorrhiza        = this%mycorrhiza       / scalar
    pool3%microbial         = this%microbial        / scalar
    pool3%residue           = this%residue          / scalar
    pool3%residue_assoc     = this%residue_assoc    / scalar
  END FUNCTION Divide_pool_a_and_scalar
  !-----------------------------------------------------------------------------------------------------
  !> / (t_sb_pool_a / field2d)
  !!
  !-----------------------------------------------------------------------------------------------------
  PURE FUNCTION Divide_pool_a_and_field2d(this,field2d) RESULT(pool3)
    IMPLICIT NONE
    ! 0.1 InOut
    CLASS(t_sb_pool_a), INTENT(in) :: this
    REAL(wp),           INTENT(in) :: field2d(:,:)
    TYPE(t_sb_pool_a)              :: pool3
    ! 1.0 go
    pool3%dom               = this%dom              / field2d(:,:)
    pool3%dom_assoc         = this%dom_assoc        / field2d(:,:)
    pool3%soluable_litter   = this%soluable_litter  / field2d(:,:)
    pool3%polymeric_litter  = this%polymeric_litter / field2d(:,:)
    pool3%woody_litter      = this%woody_litter     / field2d(:,:)
    pool3%fungi             = this%fungi            / field2d(:,:)
    pool3%mycorrhiza        = this%mycorrhiza       / field2d(:,:)
    pool3%microbial         = this%microbial        / field2d(:,:)
    pool3%residue           = this%residue          / field2d(:,:)
    pool3%residue_assoc     = this%residue_assoc    / field2d(:,:)
  END FUNCTION Divide_pool_a_and_field2d
  !-----------------------------------------------------------------------------------------------------
  !> * (t_sb_pool_a * t_sb_pool_a)
  !!
  !-----------------------------------------------------------------------------------------------------
  PURE FUNCTION Multiply_pool_a(this,pool2) RESULT(pool3)
    IMPLICIT NONE
    ! 0.1 InOut
    CLASS(t_sb_pool_a), INTENT(in) :: this, pool2
    TYPE(t_sb_pool_a)              :: pool3
    ! 1.0 go
    pool3%dom               = this%dom              * pool2%dom
    pool3%dom_assoc         = this%dom_assoc        * pool2%dom_assoc
    pool3%soluable_litter   = this%soluable_litter  * pool2%soluable_litter
    pool3%polymeric_litter  = this%polymeric_litter * pool2%polymeric_litter
    pool3%woody_litter      = this%woody_litter     * pool2%woody_litter
    pool3%fungi             = this%fungi            * pool2%fungi
    pool3%mycorrhiza        = this%mycorrhiza       * pool2%mycorrhiza
    pool3%microbial         = this%microbial        * pool2%microbial
    pool3%residue           = this%residue          * pool2%residue
    pool3%residue_assoc     = this%residue_assoc    * pool2%residue_assoc
  END FUNCTION Multiply_pool_a
  !-----------------------------------------------------------------------------------------------------
  !> * (t_sb_pool_a * scalar)
  !!
  !-----------------------------------------------------------------------------------------------------
  PURE FUNCTION Multiply_pool_a_and_scalar(this,scalar) RESULT(pool3)
    IMPLICIT NONE
    ! 0.1 InOut
    CLASS(t_sb_pool_a), INTENT(in) :: this
    REAL(wp),           INTENT(in) :: scalar
    TYPE(t_sb_pool_a)              :: pool3
    ! 1.0 go
    pool3%dom               = this%dom              * scalar
    pool3%dom_assoc         = this%dom_assoc        * scalar
    pool3%soluable_litter   = this%soluable_litter  * scalar
    pool3%polymeric_litter  = this%polymeric_litter * scalar
    pool3%woody_litter      = this%woody_litter     * scalar
    pool3%fungi             = this%fungi            * scalar
    pool3%mycorrhiza        = this%mycorrhiza       * scalar
    pool3%microbial         = this%microbial        * scalar
    pool3%residue           = this%residue          * scalar
    pool3%residue_assoc     = this%residue_assoc    * scalar
  END FUNCTION Multiply_pool_a_and_scalar
  !-----------------------------------------------------------------------------------------------------
  !> * (scalar * t_sb_pool_a)
  !!
  !-----------------------------------------------------------------------------------------------------
  PURE FUNCTION Multiply_scalar_and_pool_a(scalar,this2) RESULT(pool3)
    IMPLICIT NONE
    ! 0.1 InOut
    REAL(wp),           INTENT(in) :: scalar
    CLASS(t_sb_pool_a), INTENT(in) :: this2
    TYPE(t_sb_pool_a)              :: pool3
    ! 1.0 go
    pool3%dom               = scalar * this2%dom
    pool3%dom_assoc         = scalar * this2%dom_assoc
    pool3%soluable_litter   = scalar * this2%soluable_litter
    pool3%polymeric_litter  = scalar * this2%polymeric_litter
    pool3%woody_litter      = scalar * this2%woody_litter
    pool3%fungi             = scalar * this2%fungi
    pool3%mycorrhiza        = scalar * this2%mycorrhiza
    pool3%microbial         = scalar * this2%microbial
    pool3%residue           = scalar * this2%residue
    pool3%residue_assoc     = scalar * this2%residue_assoc
  END FUNCTION Multiply_scalar_and_pool_a
  !-----------------------------------------------------------------------------------------------------
  !> * (t_sb_pool_a * field2d)
  !!
  !-----------------------------------------------------------------------------------------------------
  PURE FUNCTION Multiply_pool_a_and_field2d(this,field2d) RESULT(pool3)
    IMPLICIT NONE
    ! 0.1 InOut
    CLASS(t_sb_pool_a), INTENT(in) :: this
    REAL(wp),           INTENT(in) :: field2d(:,:)
    TYPE(t_sb_pool_a)              :: pool3
    ! 1.0 go
    pool3%dom               = this%dom              * field2d(:,:)
    pool3%dom_assoc         = this%dom_assoc        * field2d(:,:)
    pool3%soluable_litter   = this%soluable_litter  * field2d(:,:)
    pool3%polymeric_litter  = this%polymeric_litter * field2d(:,:)
    pool3%woody_litter      = this%woody_litter     * field2d(:,:)
    pool3%fungi             = this%fungi            * field2d(:,:)
    pool3%mycorrhiza        = this%mycorrhiza       * field2d(:,:)
    pool3%microbial         = this%microbial        * field2d(:,:)
    pool3%residue           = this%residue          * field2d(:,:)
    pool3%residue_assoc     = this%residue_assoc    * field2d(:,:)
  END FUNCTION Multiply_pool_a_and_field2d
  !-----------------------------------------------------------------------------------------------------
  !> * (field2d * t_sb_pool_a)
  !!
  !-----------------------------------------------------------------------------------------------------
  PURE FUNCTION Multiply_field2d_and_pool_a(field2d,this2) RESULT(pool3)
    IMPLICIT NONE
    ! 0.1 InOut
    REAL(wp),           INTENT(in) :: field2d(:,:)
    CLASS(t_sb_pool_a), INTENT(in) :: this2
    TYPE(t_sb_pool_a)              :: pool3
    ! 1.0 go
    pool3%dom               = field2d(:,:) * this2%dom
    pool3%dom_assoc         = field2d(:,:) * this2%dom_assoc
    pool3%soluable_litter   = field2d(:,:) * this2%soluable_litter
    pool3%polymeric_litter  = field2d(:,:) * this2%polymeric_litter
    pool3%woody_litter      = field2d(:,:) * this2%woody_litter
    pool3%fungi             = field2d(:,:) * this2%fungi
    pool3%mycorrhiza        = field2d(:,:) * this2%mycorrhiza
    pool3%microbial         = field2d(:,:) * this2%microbial
    pool3%residue           = field2d(:,:) * this2%residue
    pool3%residue_assoc     = field2d(:,:) * this2%residue_assoc
  END FUNCTION Multiply_field2d_and_pool_a

  !-----------------------------------------------------------------------------------------------------
  !> Initialisation t_quincy_struct_sb_element with a scalar
  !!
  !-----------------------------------------------------------------------------------------------------
  PURE ELEMENTAL SUBROUTINE Init_t_quincy_struct_sb_element(this,s)

    CLASS(t_quincy_struct_sb_element),  INTENT(inout) :: this 
    REAL(wp),                           INTENT(in)    :: s

    this%carbon%ptr(:,:,:)     = s
    this%nitrogen%ptr(:,:,:)   = s
    this%phosphorus%ptr(:,:,:) = s
    this%carbon13%ptr(:,:,:)   = s
    this%carbon14%ptr(:,:,:)   = s
    this%nitrogen15%ptr(:,:,:) = s
  END SUBROUTINE Init_t_quincy_struct_sb_element 

  !-----------------------------------------------------------------------------------------------------
  !> Initialisation t_quincy_struct_sb_pool with a scalar
  !!
  !-----------------------------------------------------------------------------------------------------
  PURE ELEMENTAL SUBROUTINE Init_t_quincy_struct_sb_pool(this,s)

    CLASS(t_quincy_struct_sb_pool), INTENT(inout) :: this 
    REAL(wp),                       INTENT(in)    :: s

    CALL this%dom%init(s)
    CALL this%dom_assoc%init(s)
    CALL this%soluable_litter%init(s)
    CALL this%polymeric_litter%init(s)
    CALL this%woody_litter%init(s)
    CALL this%fungi%init(s)
    CALL this%mycorrhiza%init(s)
    CALL this%microbial%init(s)
    CALL this%residue%init(s)
    CALL this%residue_assoc%init(s)
  END SUBROUTINE Init_t_quincy_struct_sb_pool 

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
    CLASS(t_quincy_struct_sb_element),  INTENT(inout) :: this 
    CHARACTER(len=*),                   INTENT(in)    :: prefix          !> process name
    CHARACTER(len=*),                   INTENT(in)    :: suffix          !> tile name
    INTEGER,                            INTENT(in)    :: lct_ids !(:)    !> in jsbach4: vector of land cover type IDs
    INTEGER,                            INTENT(in)    :: model_id        !> model ID model\%id
    CHARACTER(len=*),                   INTENT(in)    :: unit            !> unit of the variables in the elements
    ! ----------------------------------------------------------------------------------------------------- !
    TYPE(t_jsb_model), POINTER  :: model
    TYPE(t_jsb_grid),  POINTER  :: hgrid                        ! Horizontal grid
    TYPE(t_jsb_vgrid), POINTER  :: vgrid_soil_sb                ! Vertical grid
    INTEGER                     :: table
    CHARACTER(len=50)           :: var_name
    CHARACTER(len=100)          :: var_name_long
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':Add_element_variables_to_memory'
    ! ----------------------------------------------------------------------------------------------------- !
    model         => Get_model(model_id)
    table         = 1     ! in jsbach4 this line is: table = tables(1)
    hgrid         => Get_grid(model%grid_id)
    vgrid_soil_sb => Get_vgrid('soil_layer_sb')

    ! add variables to pool 
    var_name      = TRIM(this%name)//'_carbon'
    var_name_long = 'Carbon in soil biogeochemistry '//TRIM(this%name)
    CALL this%Add_var(TRIM(var_name), this%carbon, &
      & hgrid, vgrid_soil_sb, &
      & t_cf(TRIM(var_name), TRIM(unit), TRIM(var_name_long)), &
      & t_grib1(table, 255, grib_bits), &
      & t_grib2(255, 255, 255, grib_bits), &
      & prefix, suffix, &
      & loutput = .TRUE., &
      & lrestart = .TRUE., &
      & initval_r = 0.0_wp)
    var_name      = TRIM(this%name)//'_nitrogen'
    var_name_long = 'Nitrogen in soil biogeochemistry '//TRIM(this%name)
    CALL this%Add_var(TRIM(var_name), this%nitrogen, &
      & hgrid, vgrid_soil_sb, &
      & t_cf(TRIM(var_name), TRIM(unit), TRIM(var_name_long)), &
      & t_grib1(table, 255, grib_bits), &
      & t_grib2(255, 255, 255, grib_bits), &
      & prefix, suffix, &
      & loutput = .TRUE., &
      & lrestart = .TRUE., &
      & initval_r = 0.0_wp)
    var_name      = TRIM(this%name)//'_phosphorus'
    var_name_long = 'Phosphorus in soil biogeochemistry '//TRIM(this%name)
    CALL this%Add_var(TRIM(var_name), this%phosphorus, &
      & hgrid, vgrid_soil_sb, &
      & t_cf(TRIM(var_name), TRIM(unit), TRIM(var_name_long)), &
      & t_grib1(table, 255, grib_bits), &
      & t_grib2(255, 255, 255, grib_bits), &
      & prefix, suffix, &
      & loutput = .TRUE., &
      & lrestart = .TRUE., &
      & initval_r = 0.0_wp)
    var_name      = TRIM(this%name)//'_carbon13'
    var_name_long = 'Carbon13 in soil biogeochemistry '//TRIM(this%name)
    CALL this%Add_var(TRIM(var_name), this%carbon13, &
      & hgrid, vgrid_soil_sb, &
      & t_cf(TRIM(var_name), TRIM(unit), TRIM(var_name_long)), &
      & t_grib1(table, 255, grib_bits), &
      & t_grib2(255, 255, 255, grib_bits), &
      & prefix, suffix, &
      & loutput = .TRUE., &
      & lrestart = .TRUE., &
      & initval_r = 0.0_wp)
    var_name      = TRIM(this%name)//'_carbon14'
    var_name_long = 'Carbon14 in soil biogeochemistry '//TRIM(this%name)
    CALL this%Add_var(TRIM(var_name), this%carbon14, &
      & hgrid, vgrid_soil_sb, &
      & t_cf(TRIM(var_name), TRIM(unit), TRIM(var_name_long)), &
      & t_grib1(table, 255, grib_bits), &
      & t_grib2(255, 255, 255, grib_bits), &
      & prefix, suffix, &
      & loutput = .TRUE., &
      & lrestart = .TRUE., &
      & initval_r = 0.0_wp)
    var_name      = TRIM(this%name)//'_nitrogen15'
    var_name_long = 'Nitrogen15 in soil biogeochemistry '//TRIM(this%name)
    CALL this%Add_var(TRIM(var_name), this%nitrogen15, &
      & hgrid, vgrid_soil_sb, &
      & t_cf(TRIM(var_name), TRIM(unit), TRIM(var_name_long)), &
      & t_grib1(table, 255, grib_bits), &
      & t_grib2(255, 255, 255, grib_bits), &
      & prefix, suffix, &
      & loutput = .TRUE., &
      & lrestart = .TRUE., &
      & initval_r = 0.0_wp)
  END SUBROUTINE Add_element_variables_to_memory

  !-----------------------------------------------------------------------------------------------------
  !> Add all elements of a pool TYPE to the process memory
  !!
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE Add_pool_elements_to_memory(this, prefix, suffix, lct_ids, model_id, unit)
  
    USE mo_jsb_varlist,       ONLY: BASIC, MEDIUM, FULL
    USE mo_jsb_io,            ONLY: grib_bits, t_cf, t_grib1, t_grib2
    USE mo_jsb_grid_class,    ONLY: t_jsb_grid, t_jsb_vgrid
    USE mo_jsb_grid,          ONLY: Get_grid, Get_vgrid
    USE mo_jsb_class,         ONLY: Get_model
    USE mo_jsb_model_class,   ONLY: t_jsb_model    
    ! ----------------------------------------------------------------------------------------------------- !
    CLASS(t_quincy_struct_sb_pool),  INTENT(inout) :: this 
    CHARACTER(len=*),                INTENT(in)    :: prefix          !> process name
    CHARACTER(len=*),                INTENT(in)    :: suffix          !> tile name
    INTEGER,                         INTENT(in)    :: lct_ids !(:)    !> in jsbach4: vector of land cover type IDs
    INTEGER,                         INTENT(in)    :: model_id        !> model ID model\%id
    CHARACTER(len=*),                INTENT(in)    :: unit            !> unit of the variables in the elements
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
    
    ! run add_var()
    element_name      = TRIM(this%name)//'_dom'
    element_name_long = 'Dissolved organic matter element of SB'//TRIM(this%name)
    CALL this%Add_pool(TRIM(element_name), this%dom, &
      & hgrid, surface, &
      & t_cf(TRIM(element_name), TRIM(unit), TRIM(element_name_long)), &
      & t_grib1(table, 255, grib_bits), &
      & t_grib2(255, 255, 255, grib_bits), &
      & prefix, suffix, &
      & loutput = .TRUE., &
      & lrestart = .TRUE.)
    element_name      = TRIM(this%name)//'_dom_assoc'
    element_name_long = 'Minerally associated disolved organic matter element of SB'//TRIM(this%name)
    CALL this%Add_pool(TRIM(element_name), this%dom_assoc, &
      & hgrid, surface, &
      & t_cf(TRIM(element_name), TRIM(unit), TRIM(element_name_long)), &
      & t_grib1(table, 255, grib_bits), &
      & t_grib2(255, 255, 255, grib_bits), &
      & prefix, suffix, &
      & loutput = .TRUE., &
      & lrestart = .TRUE.)
    element_name      = TRIM(this%name)//'_soluable_litter'
    element_name_long = 'Metabolic litter (sugars, cellulose etc.) element of SB'//TRIM(this%name)
    CALL this%Add_pool(TRIM(element_name), this%soluable_litter, &
      & hgrid, surface, &
      & t_cf(TRIM(element_name), TRIM(unit), TRIM(element_name_long)), &
      & t_grib1(table, 255, grib_bits), &
      & t_grib2(255, 255, 255, grib_bits), &
      & prefix, suffix, &
      & loutput = .TRUE., &
      & lrestart = .TRUE.)
    element_name      = TRIM(this%name)//'_polymeric_litter'
    element_name_long = 'Structural litter (lignified litter) element of SB'//TRIM(this%name)
    CALL this%Add_pool(TRIM(element_name), this%polymeric_litter, &
      & hgrid, surface, &
      & t_cf(TRIM(element_name), TRIM(unit), TRIM(element_name_long)), &
      & t_grib1(table, 255, grib_bits), &
      & t_grib2(255, 255, 255, grib_bits), &
      & prefix, suffix, &
      & loutput = .TRUE., &
      & lrestart = .TRUE.)
    element_name      = TRIM(this%name)//'_woody_litter'
    element_name_long = 'Physically protected wood litter element of SB'//TRIM(this%name)
    CALL this%Add_pool(TRIM(element_name), this%woody_litter, &
      & hgrid, surface, &
      & t_cf(TRIM(element_name), TRIM(unit), TRIM(element_name_long)), &
      & t_grib1(table, 255, grib_bits), &
      & t_grib2(255, 255, 255, grib_bits), &
      & prefix, suffix, &
      & loutput = .TRUE., &
      & lrestart = .TRUE.)
    element_name      = TRIM(this%name)//'_fungi'
    element_name_long = 'Saprophytic fungal biomass element of SB'//TRIM(this%name)
    CALL this%Add_pool(TRIM(element_name), this%fungi, &
      & hgrid, surface, &
      & t_cf(TRIM(element_name), TRIM(unit), TRIM(element_name_long)), &
      & t_grib1(table, 255, grib_bits), &
      & t_grib2(255, 255, 255, grib_bits), &
      & prefix, suffix, &
      & loutput = .TRUE., &
      & lrestart = .TRUE.)
    element_name      = TRIM(this%name)//'_mycorrhiza'
    element_name_long = 'Mycorrhizal fungal biomass element of SB'//TRIM(this%name)
    CALL this%Add_pool(TRIM(element_name), this%mycorrhiza, &
      & hgrid, surface, &
      & t_cf(TRIM(element_name), TRIM(unit), TRIM(element_name_long)), &
      & t_grib1(table, 255, grib_bits), &
      & t_grib2(255, 255, 255, grib_bits), &
      & prefix, suffix, &
      & loutput = .TRUE., &
      & lrestart = .TRUE.)
    element_name      = TRIM(this%name)//'_microbial'
    element_name_long = 'Microbial biomass (bacteria) element of SB'//TRIM(this%name)
    CALL this%Add_pool(TRIM(element_name), this%microbial, &
      & hgrid, surface, &
      & t_cf(TRIM(element_name), TRIM(unit), TRIM(element_name_long)), &
      & t_grib1(table, 255, grib_bits), &
      & t_grib2(255, 255, 255, grib_bits), &
      & prefix, suffix, &
      & loutput = .TRUE., &
      & lrestart = .TRUE.)
    element_name      = TRIM(this%name)//'_residue'
    element_name_long = 'Necromass of fungi, microbial and mycorrhiza element of SB'//TRIM(this%name)
    CALL this%Add_pool(TRIM(element_name), this%residue, &
      & hgrid, surface, &
      & t_cf(TRIM(element_name), TRIM(unit), TRIM(element_name_long)), &
      & t_grib1(table, 255, grib_bits), &
      & t_grib2(255, 255, 255, grib_bits), &
      & prefix, suffix, &
      & loutput = .TRUE., &
      & lrestart = .TRUE.)
    element_name      = TRIM(this%name)//'_residue_assoc'
    element_name_long = 'Minerally associated microbial necromass element of SB'//TRIM(this%name)
    CALL this%Add_pool(TRIM(element_name), this%residue_assoc, &
      & hgrid, surface, &
      & t_cf(TRIM(element_name), TRIM(unit), TRIM(element_name_long)), &
      & t_grib1(table, 255, grib_bits), &
      & t_grib2(255, 255, 255, grib_bits), &
      & prefix, suffix, &
      & loutput = .TRUE., &
      & lrestart = .TRUE.)  
  END SUBROUTINE Add_pool_elements_to_memory

  !-----------------------------------------------------------------------------------------------------
  !> ALLOCATE variables in sb_pool_a and sb_element_a structures
  !!
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE Init_sb_pools_elements_a(ics, ice, nc, iblk, nsoil_sb)

    INTEGER, INTENT(in)     :: ics, ice, nc, iblk, nsoil_sb   !< information about dimensions of variables in memory

    INTEGER                 :: npoints                        !< number of points in the chunk

    ! set dimension to number of points in this chunk
    npoints = nc

    !! 1.0  ALLOCATE variables
    ALLOCATE(sb_pool%dom%carbon(npoints,nsoil_sb))
    ALLOCATE(sb_pool%dom%nitrogen(npoints,nsoil_sb))
    ALLOCATE(sb_pool%dom%phosphorus(npoints,nsoil_sb))
    ALLOCATE(sb_pool%dom%carbon13(npoints,nsoil_sb))
    ALLOCATE(sb_pool%dom%carbon14(npoints,nsoil_sb))
    ALLOCATE(sb_pool%dom%nitrogen15(npoints,nsoil_sb))
    ALLOCATE(sb_pool%dom_assoc%carbon(npoints,nsoil_sb))
    ALLOCATE(sb_pool%dom_assoc%nitrogen(npoints,nsoil_sb))
    ALLOCATE(sb_pool%dom_assoc%phosphorus(npoints,nsoil_sb))
    ALLOCATE(sb_pool%dom_assoc%carbon13(npoints,nsoil_sb))
    ALLOCATE(sb_pool%dom_assoc%carbon14(npoints,nsoil_sb))
    ALLOCATE(sb_pool%dom_assoc%nitrogen15(npoints,nsoil_sb))
    ALLOCATE(sb_pool%soluable_litter%carbon(npoints,nsoil_sb))
    ALLOCATE(sb_pool%soluable_litter%nitrogen(npoints,nsoil_sb))
    ALLOCATE(sb_pool%soluable_litter%phosphorus(npoints,nsoil_sb))
    ALLOCATE(sb_pool%soluable_litter%carbon13(npoints,nsoil_sb))
    ALLOCATE(sb_pool%soluable_litter%carbon14(npoints,nsoil_sb))
    ALLOCATE(sb_pool%soluable_litter%nitrogen15(npoints,nsoil_sb))
    ALLOCATE(sb_pool%polymeric_litter%carbon(npoints,nsoil_sb))
    ALLOCATE(sb_pool%polymeric_litter%nitrogen(npoints,nsoil_sb))
    ALLOCATE(sb_pool%polymeric_litter%phosphorus(npoints,nsoil_sb))
    ALLOCATE(sb_pool%polymeric_litter%carbon13(npoints,nsoil_sb))
    ALLOCATE(sb_pool%polymeric_litter%carbon14(npoints,nsoil_sb))
    ALLOCATE(sb_pool%polymeric_litter%nitrogen15(npoints,nsoil_sb))
    ALLOCATE(sb_pool%woody_litter%carbon(npoints,nsoil_sb))
    ALLOCATE(sb_pool%woody_litter%nitrogen(npoints,nsoil_sb))
    ALLOCATE(sb_pool%woody_litter%phosphorus(npoints,nsoil_sb))
    ALLOCATE(sb_pool%woody_litter%carbon13(npoints,nsoil_sb))
    ALLOCATE(sb_pool%woody_litter%carbon14(npoints,nsoil_sb))
    ALLOCATE(sb_pool%woody_litter%nitrogen15(npoints,nsoil_sb))
    ALLOCATE(sb_pool%fungi%carbon(npoints,nsoil_sb))
    ALLOCATE(sb_pool%fungi%nitrogen(npoints,nsoil_sb))
    ALLOCATE(sb_pool%fungi%phosphorus(npoints,nsoil_sb))
    ALLOCATE(sb_pool%fungi%carbon13(npoints,nsoil_sb))
    ALLOCATE(sb_pool%fungi%carbon14(npoints,nsoil_sb))
    ALLOCATE(sb_pool%fungi%nitrogen15(npoints,nsoil_sb))
    ALLOCATE(sb_pool%mycorrhiza%carbon(npoints,nsoil_sb))
    ALLOCATE(sb_pool%mycorrhiza%nitrogen(npoints,nsoil_sb))
    ALLOCATE(sb_pool%mycorrhiza%phosphorus(npoints,nsoil_sb))
    ALLOCATE(sb_pool%mycorrhiza%carbon13(npoints,nsoil_sb))
    ALLOCATE(sb_pool%mycorrhiza%carbon14(npoints,nsoil_sb))
    ALLOCATE(sb_pool%mycorrhiza%nitrogen15(npoints,nsoil_sb))
    ALLOCATE(sb_pool%microbial%carbon(npoints,nsoil_sb))
    ALLOCATE(sb_pool%microbial%nitrogen(npoints,nsoil_sb))
    ALLOCATE(sb_pool%microbial%phosphorus(npoints,nsoil_sb))
    ALLOCATE(sb_pool%microbial%carbon13(npoints,nsoil_sb))
    ALLOCATE(sb_pool%microbial%carbon14(npoints,nsoil_sb))
    ALLOCATE(sb_pool%microbial%nitrogen15(npoints,nsoil_sb))
    ALLOCATE(sb_pool%residue%carbon(npoints,nsoil_sb))
    ALLOCATE(sb_pool%residue%nitrogen(npoints,nsoil_sb))
    ALLOCATE(sb_pool%residue%phosphorus(npoints,nsoil_sb))
    ALLOCATE(sb_pool%residue%carbon13(npoints,nsoil_sb))
    ALLOCATE(sb_pool%residue%carbon14(npoints,nsoil_sb))
    ALLOCATE(sb_pool%residue%nitrogen15(npoints,nsoil_sb))
    ALLOCATE(sb_pool%residue_assoc%carbon(npoints,nsoil_sb))
    ALLOCATE(sb_pool%residue_assoc%nitrogen(npoints,nsoil_sb))
    ALLOCATE(sb_pool%residue_assoc%phosphorus(npoints,nsoil_sb))
    ALLOCATE(sb_pool%residue_assoc%carbon13(npoints,nsoil_sb))
    ALLOCATE(sb_pool%residue_assoc%carbon14(npoints,nsoil_sb))
    ALLOCATE(sb_pool%residue_assoc%nitrogen15(npoints,nsoil_sb))

    ALLOCATE(sb_formation%dom%carbon(npoints,nsoil_sb))
    ALLOCATE(sb_formation%dom%nitrogen(npoints,nsoil_sb))
    ALLOCATE(sb_formation%dom%phosphorus(npoints,nsoil_sb))
    ALLOCATE(sb_formation%dom%carbon13(npoints,nsoil_sb))
    ALLOCATE(sb_formation%dom%carbon14(npoints,nsoil_sb))
    ALLOCATE(sb_formation%dom%nitrogen15(npoints,nsoil_sb))
    ALLOCATE(sb_formation%dom_assoc%carbon(npoints,nsoil_sb))
    ALLOCATE(sb_formation%dom_assoc%nitrogen(npoints,nsoil_sb))
    ALLOCATE(sb_formation%dom_assoc%phosphorus(npoints,nsoil_sb))
    ALLOCATE(sb_formation%dom_assoc%carbon13(npoints,nsoil_sb))
    ALLOCATE(sb_formation%dom_assoc%carbon14(npoints,nsoil_sb))
    ALLOCATE(sb_formation%dom_assoc%nitrogen15(npoints,nsoil_sb))
    ALLOCATE(sb_formation%soluable_litter%carbon(npoints,nsoil_sb))
    ALLOCATE(sb_formation%soluable_litter%nitrogen(npoints,nsoil_sb))
    ALLOCATE(sb_formation%soluable_litter%phosphorus(npoints,nsoil_sb))
    ALLOCATE(sb_formation%soluable_litter%carbon13(npoints,nsoil_sb))
    ALLOCATE(sb_formation%soluable_litter%carbon14(npoints,nsoil_sb))
    ALLOCATE(sb_formation%soluable_litter%nitrogen15(npoints,nsoil_sb))
    ALLOCATE(sb_formation%polymeric_litter%carbon(npoints,nsoil_sb))
    ALLOCATE(sb_formation%polymeric_litter%nitrogen(npoints,nsoil_sb))
    ALLOCATE(sb_formation%polymeric_litter%phosphorus(npoints,nsoil_sb))
    ALLOCATE(sb_formation%polymeric_litter%carbon13(npoints,nsoil_sb))
    ALLOCATE(sb_formation%polymeric_litter%carbon14(npoints,nsoil_sb))
    ALLOCATE(sb_formation%polymeric_litter%nitrogen15(npoints,nsoil_sb))
    ALLOCATE(sb_formation%woody_litter%carbon(npoints,nsoil_sb))
    ALLOCATE(sb_formation%woody_litter%nitrogen(npoints,nsoil_sb))
    ALLOCATE(sb_formation%woody_litter%phosphorus(npoints,nsoil_sb))
    ALLOCATE(sb_formation%woody_litter%carbon13(npoints,nsoil_sb))
    ALLOCATE(sb_formation%woody_litter%carbon14(npoints,nsoil_sb))
    ALLOCATE(sb_formation%woody_litter%nitrogen15(npoints,nsoil_sb))
    ALLOCATE(sb_formation%fungi%carbon(npoints,nsoil_sb))
    ALLOCATE(sb_formation%fungi%nitrogen(npoints,nsoil_sb))
    ALLOCATE(sb_formation%fungi%phosphorus(npoints,nsoil_sb))
    ALLOCATE(sb_formation%fungi%carbon13(npoints,nsoil_sb))
    ALLOCATE(sb_formation%fungi%carbon14(npoints,nsoil_sb))
    ALLOCATE(sb_formation%fungi%nitrogen15(npoints,nsoil_sb))
    ALLOCATE(sb_formation%mycorrhiza%carbon(npoints,nsoil_sb))
    ALLOCATE(sb_formation%mycorrhiza%nitrogen(npoints,nsoil_sb))
    ALLOCATE(sb_formation%mycorrhiza%phosphorus(npoints,nsoil_sb))
    ALLOCATE(sb_formation%mycorrhiza%carbon13(npoints,nsoil_sb))
    ALLOCATE(sb_formation%mycorrhiza%carbon14(npoints,nsoil_sb))
    ALLOCATE(sb_formation%mycorrhiza%nitrogen15(npoints,nsoil_sb))
    ALLOCATE(sb_formation%microbial%carbon(npoints,nsoil_sb))
    ALLOCATE(sb_formation%microbial%nitrogen(npoints,nsoil_sb))
    ALLOCATE(sb_formation%microbial%phosphorus(npoints,nsoil_sb))
    ALLOCATE(sb_formation%microbial%carbon13(npoints,nsoil_sb))
    ALLOCATE(sb_formation%microbial%carbon14(npoints,nsoil_sb))
    ALLOCATE(sb_formation%microbial%nitrogen15(npoints,nsoil_sb))
    ALLOCATE(sb_formation%residue%carbon(npoints,nsoil_sb))
    ALLOCATE(sb_formation%residue%nitrogen(npoints,nsoil_sb))
    ALLOCATE(sb_formation%residue%phosphorus(npoints,nsoil_sb))
    ALLOCATE(sb_formation%residue%carbon13(npoints,nsoil_sb))
    ALLOCATE(sb_formation%residue%carbon14(npoints,nsoil_sb))
    ALLOCATE(sb_formation%residue%nitrogen15(npoints,nsoil_sb))
    ALLOCATE(sb_formation%residue_assoc%carbon(npoints,nsoil_sb))
    ALLOCATE(sb_formation%residue_assoc%nitrogen(npoints,nsoil_sb))
    ALLOCATE(sb_formation%residue_assoc%phosphorus(npoints,nsoil_sb))
    ALLOCATE(sb_formation%residue_assoc%carbon13(npoints,nsoil_sb))
    ALLOCATE(sb_formation%residue_assoc%carbon14(npoints,nsoil_sb))
    ALLOCATE(sb_formation%residue_assoc%nitrogen15(npoints,nsoil_sb))

    ALLOCATE(sb_loss%dom%carbon(npoints,nsoil_sb))
    ALLOCATE(sb_loss%dom%nitrogen(npoints,nsoil_sb))
    ALLOCATE(sb_loss%dom%phosphorus(npoints,nsoil_sb))
    ALLOCATE(sb_loss%dom%carbon13(npoints,nsoil_sb))
    ALLOCATE(sb_loss%dom%carbon14(npoints,nsoil_sb))
    ALLOCATE(sb_loss%dom%nitrogen15(npoints,nsoil_sb))
    ALLOCATE(sb_loss%dom_assoc%carbon(npoints,nsoil_sb))
    ALLOCATE(sb_loss%dom_assoc%nitrogen(npoints,nsoil_sb))
    ALLOCATE(sb_loss%dom_assoc%phosphorus(npoints,nsoil_sb))
    ALLOCATE(sb_loss%dom_assoc%carbon13(npoints,nsoil_sb))
    ALLOCATE(sb_loss%dom_assoc%carbon14(npoints,nsoil_sb))
    ALLOCATE(sb_loss%dom_assoc%nitrogen15(npoints,nsoil_sb))
    ALLOCATE(sb_loss%soluable_litter%carbon(npoints,nsoil_sb))
    ALLOCATE(sb_loss%soluable_litter%nitrogen(npoints,nsoil_sb))
    ALLOCATE(sb_loss%soluable_litter%phosphorus(npoints,nsoil_sb))
    ALLOCATE(sb_loss%soluable_litter%carbon13(npoints,nsoil_sb))
    ALLOCATE(sb_loss%soluable_litter%carbon14(npoints,nsoil_sb))
    ALLOCATE(sb_loss%soluable_litter%nitrogen15(npoints,nsoil_sb))
    ALLOCATE(sb_loss%polymeric_litter%carbon(npoints,nsoil_sb))
    ALLOCATE(sb_loss%polymeric_litter%nitrogen(npoints,nsoil_sb))
    ALLOCATE(sb_loss%polymeric_litter%phosphorus(npoints,nsoil_sb))
    ALLOCATE(sb_loss%polymeric_litter%carbon13(npoints,nsoil_sb))
    ALLOCATE(sb_loss%polymeric_litter%carbon14(npoints,nsoil_sb))
    ALLOCATE(sb_loss%polymeric_litter%nitrogen15(npoints,nsoil_sb))
    ALLOCATE(sb_loss%woody_litter%carbon(npoints,nsoil_sb))
    ALLOCATE(sb_loss%woody_litter%nitrogen(npoints,nsoil_sb))
    ALLOCATE(sb_loss%woody_litter%phosphorus(npoints,nsoil_sb))
    ALLOCATE(sb_loss%woody_litter%carbon13(npoints,nsoil_sb))
    ALLOCATE(sb_loss%woody_litter%carbon14(npoints,nsoil_sb))
    ALLOCATE(sb_loss%woody_litter%nitrogen15(npoints,nsoil_sb))
    ALLOCATE(sb_loss%fungi%carbon(npoints,nsoil_sb))
    ALLOCATE(sb_loss%fungi%nitrogen(npoints,nsoil_sb))
    ALLOCATE(sb_loss%fungi%phosphorus(npoints,nsoil_sb))
    ALLOCATE(sb_loss%fungi%carbon13(npoints,nsoil_sb))
    ALLOCATE(sb_loss%fungi%carbon14(npoints,nsoil_sb))
    ALLOCATE(sb_loss%fungi%nitrogen15(npoints,nsoil_sb))
    ALLOCATE(sb_loss%mycorrhiza%carbon(npoints,nsoil_sb))
    ALLOCATE(sb_loss%mycorrhiza%nitrogen(npoints,nsoil_sb))
    ALLOCATE(sb_loss%mycorrhiza%phosphorus(npoints,nsoil_sb))
    ALLOCATE(sb_loss%mycorrhiza%carbon13(npoints,nsoil_sb))
    ALLOCATE(sb_loss%mycorrhiza%carbon14(npoints,nsoil_sb))
    ALLOCATE(sb_loss%mycorrhiza%nitrogen15(npoints,nsoil_sb))
    ALLOCATE(sb_loss%microbial%carbon(npoints,nsoil_sb))
    ALLOCATE(sb_loss%microbial%nitrogen(npoints,nsoil_sb))
    ALLOCATE(sb_loss%microbial%phosphorus(npoints,nsoil_sb))
    ALLOCATE(sb_loss%microbial%carbon13(npoints,nsoil_sb))
    ALLOCATE(sb_loss%microbial%carbon14(npoints,nsoil_sb))
    ALLOCATE(sb_loss%microbial%nitrogen15(npoints,nsoil_sb))
    ALLOCATE(sb_loss%residue%carbon(npoints,nsoil_sb))
    ALLOCATE(sb_loss%residue%nitrogen(npoints,nsoil_sb))
    ALLOCATE(sb_loss%residue%phosphorus(npoints,nsoil_sb))
    ALLOCATE(sb_loss%residue%carbon13(npoints,nsoil_sb))
    ALLOCATE(sb_loss%residue%carbon14(npoints,nsoil_sb))
    ALLOCATE(sb_loss%residue%nitrogen15(npoints,nsoil_sb))
    ALLOCATE(sb_loss%residue_assoc%carbon(npoints,nsoil_sb))
    ALLOCATE(sb_loss%residue_assoc%nitrogen(npoints,nsoil_sb))
    ALLOCATE(sb_loss%residue_assoc%phosphorus(npoints,nsoil_sb))
    ALLOCATE(sb_loss%residue_assoc%carbon13(npoints,nsoil_sb))
    ALLOCATE(sb_loss%residue_assoc%carbon14(npoints,nsoil_sb))
    ALLOCATE(sb_loss%residue_assoc%nitrogen15(npoints,nsoil_sb))

    ALLOCATE(sb_transport%dom%carbon(npoints,nsoil_sb))
    ALLOCATE(sb_transport%dom%nitrogen(npoints,nsoil_sb))
    ALLOCATE(sb_transport%dom%phosphorus(npoints,nsoil_sb))
    ALLOCATE(sb_transport%dom%carbon13(npoints,nsoil_sb))
    ALLOCATE(sb_transport%dom%carbon14(npoints,nsoil_sb))
    ALLOCATE(sb_transport%dom%nitrogen15(npoints,nsoil_sb))
    ALLOCATE(sb_transport%dom_assoc%carbon(npoints,nsoil_sb))
    ALLOCATE(sb_transport%dom_assoc%nitrogen(npoints,nsoil_sb))
    ALLOCATE(sb_transport%dom_assoc%phosphorus(npoints,nsoil_sb))
    ALLOCATE(sb_transport%dom_assoc%carbon13(npoints,nsoil_sb))
    ALLOCATE(sb_transport%dom_assoc%carbon14(npoints,nsoil_sb))
    ALLOCATE(sb_transport%dom_assoc%nitrogen15(npoints,nsoil_sb))
    ALLOCATE(sb_transport%soluable_litter%carbon(npoints,nsoil_sb))
    ALLOCATE(sb_transport%soluable_litter%nitrogen(npoints,nsoil_sb))
    ALLOCATE(sb_transport%soluable_litter%phosphorus(npoints,nsoil_sb))
    ALLOCATE(sb_transport%soluable_litter%carbon13(npoints,nsoil_sb))
    ALLOCATE(sb_transport%soluable_litter%carbon14(npoints,nsoil_sb))
    ALLOCATE(sb_transport%soluable_litter%nitrogen15(npoints,nsoil_sb))
    ALLOCATE(sb_transport%polymeric_litter%carbon(npoints,nsoil_sb))
    ALLOCATE(sb_transport%polymeric_litter%nitrogen(npoints,nsoil_sb))
    ALLOCATE(sb_transport%polymeric_litter%phosphorus(npoints,nsoil_sb))
    ALLOCATE(sb_transport%polymeric_litter%carbon13(npoints,nsoil_sb))
    ALLOCATE(sb_transport%polymeric_litter%carbon14(npoints,nsoil_sb))
    ALLOCATE(sb_transport%polymeric_litter%nitrogen15(npoints,nsoil_sb))
    ALLOCATE(sb_transport%woody_litter%carbon(npoints,nsoil_sb))
    ALLOCATE(sb_transport%woody_litter%nitrogen(npoints,nsoil_sb))
    ALLOCATE(sb_transport%woody_litter%phosphorus(npoints,nsoil_sb))
    ALLOCATE(sb_transport%woody_litter%carbon13(npoints,nsoil_sb))
    ALLOCATE(sb_transport%woody_litter%carbon14(npoints,nsoil_sb))
    ALLOCATE(sb_transport%woody_litter%nitrogen15(npoints,nsoil_sb))
    ALLOCATE(sb_transport%fungi%carbon(npoints,nsoil_sb))
    ALLOCATE(sb_transport%fungi%nitrogen(npoints,nsoil_sb))
    ALLOCATE(sb_transport%fungi%phosphorus(npoints,nsoil_sb))
    ALLOCATE(sb_transport%fungi%carbon13(npoints,nsoil_sb))
    ALLOCATE(sb_transport%fungi%carbon14(npoints,nsoil_sb))
    ALLOCATE(sb_transport%fungi%nitrogen15(npoints,nsoil_sb))
    ALLOCATE(sb_transport%mycorrhiza%carbon(npoints,nsoil_sb))
    ALLOCATE(sb_transport%mycorrhiza%nitrogen(npoints,nsoil_sb))
    ALLOCATE(sb_transport%mycorrhiza%phosphorus(npoints,nsoil_sb))
    ALLOCATE(sb_transport%mycorrhiza%carbon13(npoints,nsoil_sb))
    ALLOCATE(sb_transport%mycorrhiza%carbon14(npoints,nsoil_sb))
    ALLOCATE(sb_transport%mycorrhiza%nitrogen15(npoints,nsoil_sb))
    ALLOCATE(sb_transport%microbial%carbon(npoints,nsoil_sb))
    ALLOCATE(sb_transport%microbial%nitrogen(npoints,nsoil_sb))
    ALLOCATE(sb_transport%microbial%phosphorus(npoints,nsoil_sb))
    ALLOCATE(sb_transport%microbial%carbon13(npoints,nsoil_sb))
    ALLOCATE(sb_transport%microbial%carbon14(npoints,nsoil_sb))
    ALLOCATE(sb_transport%microbial%nitrogen15(npoints,nsoil_sb))
    ALLOCATE(sb_transport%residue%carbon(npoints,nsoil_sb))
    ALLOCATE(sb_transport%residue%nitrogen(npoints,nsoil_sb))
    ALLOCATE(sb_transport%residue%phosphorus(npoints,nsoil_sb))
    ALLOCATE(sb_transport%residue%carbon13(npoints,nsoil_sb))
    ALLOCATE(sb_transport%residue%carbon14(npoints,nsoil_sb))
    ALLOCATE(sb_transport%residue%nitrogen15(npoints,nsoil_sb))
    ALLOCATE(sb_transport%residue_assoc%carbon(npoints,nsoil_sb))
    ALLOCATE(sb_transport%residue_assoc%nitrogen(npoints,nsoil_sb))
    ALLOCATE(sb_transport%residue_assoc%phosphorus(npoints,nsoil_sb))
    ALLOCATE(sb_transport%residue_assoc%carbon13(npoints,nsoil_sb))
    ALLOCATE(sb_transport%residue_assoc%carbon14(npoints,nsoil_sb))
    ALLOCATE(sb_transport%residue_assoc%nitrogen15(npoints,nsoil_sb))

    ALLOCATE(sb_delta_pools%carbon(npoints,nsoil_sb))
    ALLOCATE(sb_delta_pools%nitrogen(npoints,nsoil_sb))
    ALLOCATE(sb_delta_pools%phosphorus(npoints,nsoil_sb))
    ALLOCATE(sb_delta_pools%carbon13(npoints,nsoil_sb))
    ALLOCATE(sb_delta_pools%carbon14(npoints,nsoil_sb))
    ALLOCATE(sb_delta_pools%nitrogen15(npoints,nsoil_sb))

    ALLOCATE(sb_mycorrhiza_export%carbon(npoints,nsoil_sb))
    ALLOCATE(sb_mycorrhiza_export%nitrogen(npoints,nsoil_sb))
    ALLOCATE(sb_mycorrhiza_export%phosphorus(npoints,nsoil_sb))
    ALLOCATE(sb_mycorrhiza_export%carbon13(npoints,nsoil_sb))
    ALLOCATE(sb_mycorrhiza_export%carbon14(npoints,nsoil_sb))
    ALLOCATE(sb_mycorrhiza_export%nitrogen15(npoints,nsoil_sb))
  END SUBROUTINE Init_sb_pools_elements_a

  !-----------------------------------------------------------------------------------------------------
  !> DEALLOCATE variables from sb_pool_a and sb_element_a structures at the program end
  !!
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE Deallocate_var_sb_pools_elements_a

    !< 1.0  DEALLOCATE variables
    !! 
    DEALLOCATE(sb_pool%dom%carbon)
    DEALLOCATE(sb_pool%dom%nitrogen)
    DEALLOCATE(sb_pool%dom%phosphorus)
    DEALLOCATE(sb_pool%dom%carbon13)
    DEALLOCATE(sb_pool%dom%carbon14)
    DEALLOCATE(sb_pool%dom%nitrogen15)
    DEALLOCATE(sb_pool%dom_assoc%carbon)
    DEALLOCATE(sb_pool%dom_assoc%nitrogen)
    DEALLOCATE(sb_pool%dom_assoc%phosphorus)
    DEALLOCATE(sb_pool%dom_assoc%carbon13)
    DEALLOCATE(sb_pool%dom_assoc%carbon14)
    DEALLOCATE(sb_pool%dom_assoc%nitrogen15)
    DEALLOCATE(sb_pool%soluable_litter%carbon)
    DEALLOCATE(sb_pool%soluable_litter%nitrogen)
    DEALLOCATE(sb_pool%soluable_litter%phosphorus)
    DEALLOCATE(sb_pool%soluable_litter%carbon13)
    DEALLOCATE(sb_pool%soluable_litter%carbon14)
    DEALLOCATE(sb_pool%soluable_litter%nitrogen15)
    DEALLOCATE(sb_pool%polymeric_litter%carbon)
    DEALLOCATE(sb_pool%polymeric_litter%nitrogen)
    DEALLOCATE(sb_pool%polymeric_litter%phosphorus)
    DEALLOCATE(sb_pool%polymeric_litter%carbon13)
    DEALLOCATE(sb_pool%polymeric_litter%carbon14)
    DEALLOCATE(sb_pool%polymeric_litter%nitrogen15)
    DEALLOCATE(sb_pool%woody_litter%carbon)
    DEALLOCATE(sb_pool%woody_litter%nitrogen)
    DEALLOCATE(sb_pool%woody_litter%phosphorus)
    DEALLOCATE(sb_pool%woody_litter%carbon13)
    DEALLOCATE(sb_pool%woody_litter%carbon14)
    DEALLOCATE(sb_pool%woody_litter%nitrogen15)
    DEALLOCATE(sb_pool%fungi%carbon)
    DEALLOCATE(sb_pool%fungi%nitrogen)
    DEALLOCATE(sb_pool%fungi%phosphorus)
    DEALLOCATE(sb_pool%fungi%carbon13)
    DEALLOCATE(sb_pool%fungi%carbon14)
    DEALLOCATE(sb_pool%fungi%nitrogen15)
    DEALLOCATE(sb_pool%mycorrhiza%carbon)
    DEALLOCATE(sb_pool%mycorrhiza%nitrogen)
    DEALLOCATE(sb_pool%mycorrhiza%phosphorus)
    DEALLOCATE(sb_pool%mycorrhiza%carbon13)
    DEALLOCATE(sb_pool%mycorrhiza%carbon14)
    DEALLOCATE(sb_pool%mycorrhiza%nitrogen15)
    DEALLOCATE(sb_pool%microbial%carbon)
    DEALLOCATE(sb_pool%microbial%nitrogen)
    DEALLOCATE(sb_pool%microbial%phosphorus)
    DEALLOCATE(sb_pool%microbial%carbon13)
    DEALLOCATE(sb_pool%microbial%carbon14)
    DEALLOCATE(sb_pool%microbial%nitrogen15)
    DEALLOCATE(sb_pool%residue%carbon)
    DEALLOCATE(sb_pool%residue%nitrogen)
    DEALLOCATE(sb_pool%residue%phosphorus)
    DEALLOCATE(sb_pool%residue%carbon13)
    DEALLOCATE(sb_pool%residue%carbon14)
    DEALLOCATE(sb_pool%residue%nitrogen15)
    DEALLOCATE(sb_pool%residue_assoc%carbon)
    DEALLOCATE(sb_pool%residue_assoc%nitrogen)
    DEALLOCATE(sb_pool%residue_assoc%phosphorus)
    DEALLOCATE(sb_pool%residue_assoc%carbon13)
    DEALLOCATE(sb_pool%residue_assoc%carbon14)
    DEALLOCATE(sb_pool%residue_assoc%nitrogen15)

    DEALLOCATE(sb_formation%dom%carbon)
    DEALLOCATE(sb_formation%dom%nitrogen)
    DEALLOCATE(sb_formation%dom%phosphorus)
    DEALLOCATE(sb_formation%dom%carbon13)
    DEALLOCATE(sb_formation%dom%carbon14)
    DEALLOCATE(sb_formation%dom%nitrogen15)
    DEALLOCATE(sb_formation%dom_assoc%carbon)
    DEALLOCATE(sb_formation%dom_assoc%nitrogen)
    DEALLOCATE(sb_formation%dom_assoc%phosphorus)
    DEALLOCATE(sb_formation%dom_assoc%carbon13)
    DEALLOCATE(sb_formation%dom_assoc%carbon14)
    DEALLOCATE(sb_formation%dom_assoc%nitrogen15)
    DEALLOCATE(sb_formation%soluable_litter%carbon)
    DEALLOCATE(sb_formation%soluable_litter%nitrogen)
    DEALLOCATE(sb_formation%soluable_litter%phosphorus)
    DEALLOCATE(sb_formation%soluable_litter%carbon13)
    DEALLOCATE(sb_formation%soluable_litter%carbon14)
    DEALLOCATE(sb_formation%soluable_litter%nitrogen15)
    DEALLOCATE(sb_formation%polymeric_litter%carbon)
    DEALLOCATE(sb_formation%polymeric_litter%nitrogen)
    DEALLOCATE(sb_formation%polymeric_litter%phosphorus)
    DEALLOCATE(sb_formation%polymeric_litter%carbon13)
    DEALLOCATE(sb_formation%polymeric_litter%carbon14)
    DEALLOCATE(sb_formation%polymeric_litter%nitrogen15)
    DEALLOCATE(sb_formation%woody_litter%carbon)
    DEALLOCATE(sb_formation%woody_litter%nitrogen)
    DEALLOCATE(sb_formation%woody_litter%phosphorus)
    DEALLOCATE(sb_formation%woody_litter%carbon13)
    DEALLOCATE(sb_formation%woody_litter%carbon14)
    DEALLOCATE(sb_formation%woody_litter%nitrogen15)
    DEALLOCATE(sb_formation%fungi%carbon)
    DEALLOCATE(sb_formation%fungi%nitrogen)
    DEALLOCATE(sb_formation%fungi%phosphorus)
    DEALLOCATE(sb_formation%fungi%carbon13)
    DEALLOCATE(sb_formation%fungi%carbon14)
    DEALLOCATE(sb_formation%fungi%nitrogen15)
    DEALLOCATE(sb_formation%mycorrhiza%carbon)
    DEALLOCATE(sb_formation%mycorrhiza%nitrogen)
    DEALLOCATE(sb_formation%mycorrhiza%phosphorus)
    DEALLOCATE(sb_formation%mycorrhiza%carbon13)
    DEALLOCATE(sb_formation%mycorrhiza%carbon14)
    DEALLOCATE(sb_formation%mycorrhiza%nitrogen15)
    DEALLOCATE(sb_formation%microbial%carbon)
    DEALLOCATE(sb_formation%microbial%nitrogen)
    DEALLOCATE(sb_formation%microbial%phosphorus)
    DEALLOCATE(sb_formation%microbial%carbon13)
    DEALLOCATE(sb_formation%microbial%carbon14)
    DEALLOCATE(sb_formation%microbial%nitrogen15)
    DEALLOCATE(sb_formation%residue%carbon)
    DEALLOCATE(sb_formation%residue%nitrogen)
    DEALLOCATE(sb_formation%residue%phosphorus)
    DEALLOCATE(sb_formation%residue%carbon13)
    DEALLOCATE(sb_formation%residue%carbon14)
    DEALLOCATE(sb_formation%residue%nitrogen15)
    DEALLOCATE(sb_formation%residue_assoc%carbon)
    DEALLOCATE(sb_formation%residue_assoc%nitrogen)
    DEALLOCATE(sb_formation%residue_assoc%phosphorus)
    DEALLOCATE(sb_formation%residue_assoc%carbon13)
    DEALLOCATE(sb_formation%residue_assoc%carbon14)
    DEALLOCATE(sb_formation%residue_assoc%nitrogen15)

    DEALLOCATE(sb_loss%dom%carbon)
    DEALLOCATE(sb_loss%dom%nitrogen)
    DEALLOCATE(sb_loss%dom%phosphorus)
    DEALLOCATE(sb_loss%dom%carbon13)
    DEALLOCATE(sb_loss%dom%carbon14)
    DEALLOCATE(sb_loss%dom%nitrogen15)
    DEALLOCATE(sb_loss%dom_assoc%carbon)
    DEALLOCATE(sb_loss%dom_assoc%nitrogen)
    DEALLOCATE(sb_loss%dom_assoc%phosphorus)
    DEALLOCATE(sb_loss%dom_assoc%carbon13)
    DEALLOCATE(sb_loss%dom_assoc%carbon14)
    DEALLOCATE(sb_loss%dom_assoc%nitrogen15)
    DEALLOCATE(sb_loss%soluable_litter%carbon)
    DEALLOCATE(sb_loss%soluable_litter%nitrogen)
    DEALLOCATE(sb_loss%soluable_litter%phosphorus)
    DEALLOCATE(sb_loss%soluable_litter%carbon13)
    DEALLOCATE(sb_loss%soluable_litter%carbon14)
    DEALLOCATE(sb_loss%soluable_litter%nitrogen15)
    DEALLOCATE(sb_loss%polymeric_litter%carbon)
    DEALLOCATE(sb_loss%polymeric_litter%nitrogen)
    DEALLOCATE(sb_loss%polymeric_litter%phosphorus)
    DEALLOCATE(sb_loss%polymeric_litter%carbon13)
    DEALLOCATE(sb_loss%polymeric_litter%carbon14)
    DEALLOCATE(sb_loss%polymeric_litter%nitrogen15)
    DEALLOCATE(sb_loss%woody_litter%carbon)
    DEALLOCATE(sb_loss%woody_litter%nitrogen)
    DEALLOCATE(sb_loss%woody_litter%phosphorus)
    DEALLOCATE(sb_loss%woody_litter%carbon13)
    DEALLOCATE(sb_loss%woody_litter%carbon14)
    DEALLOCATE(sb_loss%woody_litter%nitrogen15)
    DEALLOCATE(sb_loss%fungi%carbon)
    DEALLOCATE(sb_loss%fungi%nitrogen)
    DEALLOCATE(sb_loss%fungi%phosphorus)
    DEALLOCATE(sb_loss%fungi%carbon13)
    DEALLOCATE(sb_loss%fungi%carbon14)
    DEALLOCATE(sb_loss%fungi%nitrogen15)
    DEALLOCATE(sb_loss%mycorrhiza%carbon)
    DEALLOCATE(sb_loss%mycorrhiza%nitrogen)
    DEALLOCATE(sb_loss%mycorrhiza%phosphorus)
    DEALLOCATE(sb_loss%mycorrhiza%carbon13)
    DEALLOCATE(sb_loss%mycorrhiza%carbon14)
    DEALLOCATE(sb_loss%mycorrhiza%nitrogen15)
    DEALLOCATE(sb_loss%microbial%carbon)
    DEALLOCATE(sb_loss%microbial%nitrogen)
    DEALLOCATE(sb_loss%microbial%phosphorus)
    DEALLOCATE(sb_loss%microbial%carbon13)
    DEALLOCATE(sb_loss%microbial%carbon14)
    DEALLOCATE(sb_loss%microbial%nitrogen15)
    DEALLOCATE(sb_loss%residue%carbon)
    DEALLOCATE(sb_loss%residue%nitrogen)
    DEALLOCATE(sb_loss%residue%phosphorus)
    DEALLOCATE(sb_loss%residue%carbon13)
    DEALLOCATE(sb_loss%residue%carbon14)
    DEALLOCATE(sb_loss%residue%nitrogen15)
    DEALLOCATE(sb_loss%residue_assoc%carbon)
    DEALLOCATE(sb_loss%residue_assoc%nitrogen)
    DEALLOCATE(sb_loss%residue_assoc%phosphorus)
    DEALLOCATE(sb_loss%residue_assoc%carbon13)
    DEALLOCATE(sb_loss%residue_assoc%carbon14)
    DEALLOCATE(sb_loss%residue_assoc%nitrogen15)

    DEALLOCATE(sb_transport%dom%carbon)
    DEALLOCATE(sb_transport%dom%nitrogen)
    DEALLOCATE(sb_transport%dom%phosphorus)
    DEALLOCATE(sb_transport%dom%carbon13)
    DEALLOCATE(sb_transport%dom%carbon14)
    DEALLOCATE(sb_transport%dom%nitrogen15)
    DEALLOCATE(sb_transport%dom_assoc%carbon)
    DEALLOCATE(sb_transport%dom_assoc%nitrogen)
    DEALLOCATE(sb_transport%dom_assoc%phosphorus)
    DEALLOCATE(sb_transport%dom_assoc%carbon13)
    DEALLOCATE(sb_transport%dom_assoc%carbon14)
    DEALLOCATE(sb_transport%dom_assoc%nitrogen15)
    DEALLOCATE(sb_transport%soluable_litter%carbon)
    DEALLOCATE(sb_transport%soluable_litter%nitrogen)
    DEALLOCATE(sb_transport%soluable_litter%phosphorus)
    DEALLOCATE(sb_transport%soluable_litter%carbon13)
    DEALLOCATE(sb_transport%soluable_litter%carbon14)
    DEALLOCATE(sb_transport%soluable_litter%nitrogen15)
    DEALLOCATE(sb_transport%polymeric_litter%carbon)
    DEALLOCATE(sb_transport%polymeric_litter%nitrogen)
    DEALLOCATE(sb_transport%polymeric_litter%phosphorus)
    DEALLOCATE(sb_transport%polymeric_litter%carbon13)
    DEALLOCATE(sb_transport%polymeric_litter%carbon14)
    DEALLOCATE(sb_transport%polymeric_litter%nitrogen15)
    DEALLOCATE(sb_transport%woody_litter%carbon)
    DEALLOCATE(sb_transport%woody_litter%nitrogen)
    DEALLOCATE(sb_transport%woody_litter%phosphorus)
    DEALLOCATE(sb_transport%woody_litter%carbon13)
    DEALLOCATE(sb_transport%woody_litter%carbon14)
    DEALLOCATE(sb_transport%woody_litter%nitrogen15)
    DEALLOCATE(sb_transport%fungi%carbon)
    DEALLOCATE(sb_transport%fungi%nitrogen)
    DEALLOCATE(sb_transport%fungi%phosphorus)
    DEALLOCATE(sb_transport%fungi%carbon13)
    DEALLOCATE(sb_transport%fungi%carbon14)
    DEALLOCATE(sb_transport%fungi%nitrogen15)
    DEALLOCATE(sb_transport%mycorrhiza%carbon)
    DEALLOCATE(sb_transport%mycorrhiza%nitrogen)
    DEALLOCATE(sb_transport%mycorrhiza%phosphorus)
    DEALLOCATE(sb_transport%mycorrhiza%carbon13)
    DEALLOCATE(sb_transport%mycorrhiza%carbon14)
    DEALLOCATE(sb_transport%mycorrhiza%nitrogen15)
    DEALLOCATE(sb_transport%microbial%carbon)
    DEALLOCATE(sb_transport%microbial%nitrogen)
    DEALLOCATE(sb_transport%microbial%phosphorus)
    DEALLOCATE(sb_transport%microbial%carbon13)
    DEALLOCATE(sb_transport%microbial%carbon14)
    DEALLOCATE(sb_transport%microbial%nitrogen15)
    DEALLOCATE(sb_transport%residue%carbon)
    DEALLOCATE(sb_transport%residue%nitrogen)
    DEALLOCATE(sb_transport%residue%phosphorus)
    DEALLOCATE(sb_transport%residue%carbon13)
    DEALLOCATE(sb_transport%residue%carbon14)
    DEALLOCATE(sb_transport%residue%nitrogen15)
    DEALLOCATE(sb_transport%residue_assoc%carbon)
    DEALLOCATE(sb_transport%residue_assoc%nitrogen)
    DEALLOCATE(sb_transport%residue_assoc%phosphorus)
    DEALLOCATE(sb_transport%residue_assoc%carbon13)
    DEALLOCATE(sb_transport%residue_assoc%carbon14)
    DEALLOCATE(sb_transport%residue_assoc%nitrogen15)

    DEALLOCATE(sb_delta_pools%carbon)
    DEALLOCATE(sb_delta_pools%nitrogen)
    DEALLOCATE(sb_delta_pools%phosphorus)
    DEALLOCATE(sb_delta_pools%carbon13)
    DEALLOCATE(sb_delta_pools%carbon14)
    DEALLOCATE(sb_delta_pools%nitrogen15)

    DEALLOCATE(sb_mycorrhiza_export%carbon)
    DEALLOCATE(sb_mycorrhiza_export%nitrogen)
    DEALLOCATE(sb_mycorrhiza_export%phosphorus)
    DEALLOCATE(sb_mycorrhiza_export%carbon13)
    DEALLOCATE(sb_mycorrhiza_export%carbon14)
    DEALLOCATE(sb_mycorrhiza_export%nitrogen15)
  END SUBROUTINE Deallocate_var_sb_pools_elements_a

!! #################################################################################################################################
!! quincy-no-jsb_end
!! #################################################################################################################################
#else
!! #################################################################################################################################
!! code specific to QUINCY with JSBACH4 (until the end of this module)
!! #################################################################################################################################

  USE mo_kind,                   ONLY: wp
  USE mo_jsb_control,            ONLY: debug_on
  USE mo_exception,              ONLY: finish, message
  USE mo_util,                   ONLY: One_of
  USE mo_jsb_utils_iface,        ONLY: assign_if_present_allocatable
  USE mo_jsb_var_class,          ONLY: t_jsb_var_real3d
  USE mo_jsb_pool_class,         ONLY: t_jsb_pool_elem3d, ELEM_C, ELEM_N, ELEM_P, ELEM_C13, ELEM_C14, ELEM_N15
  USE mo_jsb_subset,             ONLY: t_subset, ON_DOMAIN, ON_CHUNK

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: all_elements, Create_sb_components_pool
  PUBLIC :: t_sb_main_pool, t_sb_components_pool, t_sb_elements_pool

  !-----------------------------------------------------------------------------------------------------
  !> Type definition for sb_elements_pool
  !! pool with only variables
  !! last pool in line
  !-----------------------------------------------------------------------------------------------------
  TYPE, EXTENDS(t_jsb_pool_elem3d) :: t_sb_elements_pool
    !
    ! element variable pointers are declared in mo_jsb_pool_class:t_jsb_pool_elem3d (see also mo_sb_memory_class:all_elements(:))
    !
    !! var for operator overload
    ! 2d ON_CHUNK
    REAL(wp), ALLOCATABLE   :: op_ov_hlp_2d_c(:,:)
    REAL(wp), ALLOCATABLE   :: op_ov_hlp_2d_n(:,:)
    REAL(wp), ALLOCATABLE   :: op_ov_hlp_2d_p(:,:)
    REAL(wp), ALLOCATABLE   :: op_ov_hlp_2d_c13(:,:)
    REAL(wp), ALLOCATABLE   :: op_ov_hlp_2d_c14(:,:)
    REAL(wp), ALLOCATABLE   :: op_ov_hlp_2d_n15(:,:)
    ! 3d ON_DOMAIN
    REAL(wp), ALLOCATABLE   :: op_ov_hlp_3d_c(:,:,:)
    REAL(wp), ALLOCATABLE   :: op_ov_hlp_3d_n(:,:,:)
    REAL(wp), ALLOCATABLE   :: op_ov_hlp_3d_p(:,:,:)
    REAL(wp), ALLOCATABLE   :: op_ov_hlp_3d_c13(:,:,:)
    REAL(wp), ALLOCATABLE   :: op_ov_hlp_3d_c14(:,:,:)
    REAL(wp), ALLOCATABLE   :: op_ov_hlp_3d_n15(:,:,:)
  CONTAINS
    !! infrastructure
    PROCEDURE :: Get_element_name_by_id     => t_sb_elements_pool_get_element_name_by_id
    PROCEDURE :: Init                       => t_sb_elements_pool_init
    PROCEDURE :: Init_op_ov_hlp             => t_sb_elements_pool_op_ov_hlp_init             ! init hlp object for operator overload
    PROCEDURE :: Associate_chunk_var_ptr    => t_sb_elements_pool_associate_chunk_variable_pointers
    PROCEDURE :: Associate_domain_var_ptr   => t_sb_elements_pool_associate_domain_variable_pointers
    PROCEDURE :: Set_to_scalar_on_chunk     => t_sb_elements_pool_set_all_elements_to_scalar_value_on_chunk  ! e.g., set all to zero
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
  TYPE, EXTENDS(t_jsb_pool_elem3d) :: t_sb_components_pool
    !
    ! element variable pointers are declared in mo_jsb_pool_class:t_jsb_pool_elem3d (see also mo_sb_memory_class:all_elements(:))
    !
    !! var for operator overload
    ! 2d ON_CHUNK
    REAL(wp), ALLOCATABLE   :: op_ov_hlp_2d_c(:,:)
    REAL(wp), ALLOCATABLE   :: op_ov_hlp_2d_n(:,:)
    REAL(wp), ALLOCATABLE   :: op_ov_hlp_2d_p(:,:)
    REAL(wp), ALLOCATABLE   :: op_ov_hlp_2d_c13(:,:)
    REAL(wp), ALLOCATABLE   :: op_ov_hlp_2d_c14(:,:)
    REAL(wp), ALLOCATABLE   :: op_ov_hlp_2d_n15(:,:)
    ! 3d ON_DOMAIN
    REAL(wp), ALLOCATABLE   :: op_ov_hlp_3d_c(:,:,:)
    REAL(wp), ALLOCATABLE   :: op_ov_hlp_3d_n(:,:,:)
    REAL(wp), ALLOCATABLE   :: op_ov_hlp_3d_p(:,:,:)
    REAL(wp), ALLOCATABLE   :: op_ov_hlp_3d_c13(:,:,:)
    REAL(wp), ALLOCATABLE   :: op_ov_hlp_3d_c14(:,:,:)
    REAL(wp), ALLOCATABLE   :: op_ov_hlp_3d_n15(:,:,:)
    !
    ! component pools
    !
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
    PROCEDURE :: Get_element_name_by_id     => t_sb_components_pool_get_element_name_by_id
    PROCEDURE :: Init                       => t_sb_components_pool_init
    PROCEDURE :: Init_op_ov_hlp             => t_sb_components_pool_op_ov_hlp_create_init     ! create and init hlp object for operator overload
    PROCEDURE :: Associate_chunk_var_ptr    => t_sb_components_pool_associate_chunk_variable_pointers
    PROCEDURE :: Associate_domain_var_ptr   => t_sb_components_pool_associate_domain_variable_pointers
    PROCEDURE :: Set_to_scalar_on_chunk     => t_sb_components_pool_set_all_elements_to_scalar_value_on_chunk  ! e.g., set all to zero
    !! operator overload
    !! [def: field2d(nc,vgrid) | field3d(nc,vgrid,nblks)]
    ! +
    PROCEDURE                 :: Add_sb_comp_p_and_scalar            ! elements_pool + scalar        | ON_CHUNK & ON_DOMAIN
    PROCEDURE, PASS(struct1)  :: Add_scalar_and_sb_comp_p            ! scalar        + elements_pool | ON_CHUNK & ON_DOMAIN
    PROCEDURE                 :: Add_sb_comp_p_and_field2d           ! elements_pool + field2d       | ON_CHUNK
    PROCEDURE, PASS(struct1)  :: Add_field2d_and_sb_comp_p           ! field2d       + elements_pool | ON_CHUNK
    PROCEDURE                 :: Add_sb_comp_p_and_field3d           ! elements_pool + field3d       | ON_CHUNK & ON_DOMAIN
    PROCEDURE                 :: Add_sb_comp_p_and_sb_comp_p         ! elements_pool + elements_pool | ON_CHUNK & ON_DOMAIN
    GENERIC                   :: OPERATOR(+) => Add_sb_comp_p_and_scalar , Add_scalar_and_sb_comp_p, &
      &                                         Add_sb_comp_p_and_field2d, Add_field2d_and_sb_comp_p, &
      &                                         Add_sb_comp_p_and_field3d, &
      &                                         Add_sb_comp_p_and_sb_comp_p
    ! -
    PROCEDURE                 :: Subtract_sb_comp_p_minus_scalar     ! elements_pool - scalar        | ON_CHUNK & ON_DOMAIN
    PROCEDURE                 :: Subtract_sb_comp_p_minus_field2d    ! elements_pool - field2d       | ON_CHUNK
    PROCEDURE                 :: Subtract_sb_comp_p_minus_field3d    ! elements_pool - field3d       | ON_CHUNK & ON_DOMAIN
    PROCEDURE                 :: Subtract_sb_comp_p_minus_sb_comp_p  ! elements_pool - elements_pool | ON_CHUNK & ON_DOMAIN
    GENERIC                   :: OPERATOR(-) => Subtract_sb_comp_p_minus_scalar, &
      &                                         Subtract_sb_comp_p_minus_field2d, &
      &                                         Subtract_sb_comp_p_minus_field3d, &
      &                                         Subtract_sb_comp_p_minus_sb_comp_p
    ! *
    PROCEDURE                 :: Multiply_sb_comp_p_with_scalar      ! elements_pool * scalar        | ON_CHUNK & ON_DOMAIN
    PROCEDURE, PASS(struct1)  :: Multiply_scalar_with_sb_comp_p      ! scalar        * elements_pool | ON_CHUNK & ON_DOMAIN
    PROCEDURE                 :: Multiply_sb_comp_p_with_field2d     ! elements_pool * field2d       | ON_CHUNK
    PROCEDURE, PASS(struct1)  :: Multiply_field2d_with_sb_comp_p     ! field2d       * elements_pool | ON_CHUNK
    PROCEDURE                 :: Multiply_sb_comp_p_with_field3d     ! elements_pool * field3d       | ON_CHUNK & ON_DOMAIN
    PROCEDURE, PASS(struct1)  :: Multiply_field3d_with_sb_comp_p     ! field3d       * elements_pool | ON_CHUNK & ON_DOMAIN
    PROCEDURE                 :: Multiply_sb_comp_p_with_sb_comp_p   ! elements_pool * elements_pool | ON_CHUNK & ON_DOMAIN
    GENERIC                   :: OPERATOR(*) => Multiply_sb_comp_p_with_scalar,  Multiply_scalar_with_sb_comp_p, &
      &                                         Multiply_sb_comp_p_with_field2d, Multiply_field2d_with_sb_comp_p, &
      &                                         Multiply_sb_comp_p_with_field3d, Multiply_field3d_with_sb_comp_p, &
      &                                         Multiply_sb_comp_p_with_sb_comp_p
    ! /
    PROCEDURE                 :: Div_sb_comp_p_by_scalar             ! elements_pool / scalar        | ON_CHUNK & ON_DOMAIN
    PROCEDURE                 :: Div_sb_comp_p_by_field2d            ! elements_pool / field2d       | ON_CHUNK
    PROCEDURE                 :: Div_sb_comp_p_by_field3d            ! elements_pool / field3d       | ON_CHUNK & ON_DOMAIN
    PROCEDURE                 :: Div_sb_comp_p_by_sb_comp_p          ! elements_pool / elements_pool | ON_CHUNK & ON_DOMAIN
    GENERIC                   :: OPERATOR(/) => Div_sb_comp_p_by_scalar, &
      &                                         Div_sb_comp_p_by_field2d, &
      &                                         Div_sb_comp_p_by_field3d, &
      &                                         Div_sb_comp_p_by_sb_comp_p
    ! =
    PROCEDURE                 :: Assign_scalar_to_sb_comp_p          ! elements_pool = scalar        | ON_CHUNK & ON_DOMAIN
    PROCEDURE                 :: Assign_field2d_to_sb_comp_p         ! elements_pool = field2d(:,:)  | ON_CHUNK
    PROCEDURE                 :: Assign_field3d_to_sb_comp_p         ! elements_pool = field3d(:,:,:)| ON_CHUNK & ON_DOMAIN
    PROCEDURE                 :: Assign_sb_comp_p_to_sb_comp_p       ! elements_pool = elements_pool | ON_CHUNK & ON_DOMAIN | only numeric values of elements are assigned
    GENERIC                   :: ASSIGNMENT(=) => Assign_scalar_to_sb_comp_p , &
      &                                           Assign_field2d_to_sb_comp_p, &
      &                                           Assign_field3d_to_sb_comp_p, &
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
    !! 
    !> 2.0 components_pools
    !! pools
    TYPE(t_sb_components_pool),   POINTER :: sb_pool
    !! fluxes
    TYPE(t_sb_components_pool),   POINTER :: sb_formation
    TYPE(t_sb_components_pool),   POINTER :: sb_loss
    TYPE(t_sb_components_pool),   POINTER :: sb_transport
    !> 3.0 element_pools
    !! pools
    TYPE(t_sb_elements_pool),     POINTER :: sb_delta_pools
    !! fluxes
    TYPE(t_sb_elements_pool),     POINTER :: sb_mycorrhiza_export
  END TYPE t_sb_main_pool

  !-----------------------------------------------------------------------------------------------------
  !> list with the element_id of all elements as defined above
  !! is allocated in Init_PROC_memory
  !-----------------------------------------------------------------------------------------------------
  INTEGER, ALLOCATABLE, SAVE :: all_elements(:)

  CHARACTER(len=*), PARAMETER :: modname = 'mo_sb_memory_class_util'


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

  ! ======================================================================================================= !
  !>set all elements of the sb_elements_pool to the scalar value
  !> 
  SUBROUTINE t_sb_elements_pool_set_all_elements_to_scalar_value_on_chunk(this, scalar)

    CLASS(t_sb_elements_pool), INTENT(inout) :: this
    REAL(wp)                 , INTENT(in)    :: scalar
    
    this%carbon(:,:)       = scalar
    this%nitrogen(:,:)     = scalar
    this%phosphorus(:,:)   = scalar
    this%carbon13(:,:)     = scalar
    this%carbon14(:,:)     = scalar
    this%nitrogen15(:,:)   = scalar

  END SUBROUTINE t_sb_elements_pool_set_all_elements_to_scalar_value_on_chunk

  ! ======================================================================================================= !
  !>set all elements of the sb_components_pool to the scalar value
  !> 
  SUBROUTINE t_sb_components_pool_set_all_elements_to_scalar_value_on_chunk(this, scalar)

    CLASS(t_sb_components_pool), INTENT(inout) :: this
    REAL(wp)                   , INTENT(in)    :: scalar
    
    ! elements
    this%carbon(:,:)       = scalar
    this%nitrogen(:,:)     = scalar
    this%phosphorus(:,:)   = scalar
    this%carbon13(:,:)     = scalar
    this%carbon14(:,:)     = scalar
    this%nitrogen15(:,:)   = scalar
    ! sb_elements_pool
    CALL this%dom%Set_to_scalar_on_chunk(scalar)
    CALL this%dom_assoc%Set_to_scalar_on_chunk(scalar)
    CALL this%soluable_litter%Set_to_scalar_on_chunk(scalar)
    CALL this%polymeric_litter%Set_to_scalar_on_chunk(scalar)
    CALL this%woody_litter%Set_to_scalar_on_chunk(scalar)
    CALL this%fungi%Set_to_scalar_on_chunk(scalar)
    CALL this%mycorrhiza%Set_to_scalar_on_chunk(scalar)
    CALL this%microbial%Set_to_scalar_on_chunk(scalar)
    CALL this%residue%Set_to_scalar_on_chunk(scalar)
    CALL this%residue_assoc%Set_to_scalar_on_chunk(scalar)

  END SUBROUTINE t_sb_components_pool_set_all_elements_to_scalar_value_on_chunk  

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
  
    this%element_unit                  = element_unit
    this%contains_elem_2d_ptr_to_chunk = .TRUE.
  
    !> carbon
    IF (One_of(ELEM_C, element_list_loc) > 0) THEN
      CALL this%Add_element('carbon', 'C', ELEM_C, dim=3)
    END IF
    !> nitrogen
    IF (One_of(ELEM_N, element_list_loc) > 0) THEN
      CALL this%Add_element('nitrogen', 'N', ELEM_N, dim=3)
    END IF
    !> phosphorus
    IF (One_of(ELEM_P, element_list_loc) > 0) THEN
      CALL this%Add_element('phosphorus', 'P', ELEM_P, dim=3)
    END IF
    !> carbon13
    IF (One_of(ELEM_C13, element_list_loc) > 0) THEN
      CALL this%Add_element('carbon13', 'C13', ELEM_C13, dim=3)
    END IF
    !> carbon14
    IF (One_of(ELEM_C14, element_list_loc) > 0) THEN
      CALL this%Add_element('carbon14', 'C14', ELEM_C14, dim=3)
    END IF
    !> nitrogen15
    IF (One_of(ELEM_N15, element_list_loc) > 0) THEN
      CALL this%Add_element('nitrogen15', 'N15', ELEM_N15, dim=3)
    END IF

  END SUBROUTINE t_sb_elements_pool_init

  !-----------------------------------------------------------------------------------------------------
  !> initializing a hlp object of the SB_ elements pool "t_sb_elements_pool" for operator overload
  !! "pool_object%Init_op_ov_hlp()"
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE t_sb_elements_pool_op_ov_hlp_init(this, name, shortname, owner_model_id, subset_type, &
    &                                          nc, nvgrid, nb, l_elements, element_unit)

    CLASS(t_sb_elements_pool),  INTENT(inout), TARGET :: this               ! object of TYPE t_sb_elements_pool
    CHARACTER(LEN=*),           INTENT(in)            :: name               ! pool name       (No spaces !)
    CHARACTER(LEN=*),           INTENT(in)            :: shortname          ! pool short name (No spaces !)
    INTEGER,                    INTENT(in)            :: owner_model_id     ! model ID
    INTEGER,                    INTENT(in)            :: subset_type        ! ON_DOMAIN=1 / ON_CHUNK=2
    INTEGER,                    INTENT(in)            :: nc                 ! number of points in chunk
    INTEGER,                    INTENT(in)            :: nvgrid             ! number of layers of the vertical grid
    INTEGER,                    INTENT(in)            :: nb                 ! number of blocks
    LOGICAL,                    INTENT(in)            :: l_elements         ! IF pool has elements 
    CHARACTER(LEN=*),           INTENT(in)            :: element_unit       ! unit of elements

    CHARACTER(len=*), PARAMETER :: routine = modname//':t_sb_elements_pool_op_ov_hlp_init'

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
      ALLOCATE(this%op_ov_hlp_2d_c(nc,nvgrid))
      IF (ASSOCIATED(this%carbon))     NULLIFY(this%carbon)
      this%carbon     => this%op_ov_hlp_2d_c
      ALLOCATE(this%op_ov_hlp_2d_n(nc,nvgrid))
      IF (ASSOCIATED(this%nitrogen))   NULLIFY(this%nitrogen)
      this%nitrogen   => this%op_ov_hlp_2d_n
      ALLOCATE(this%op_ov_hlp_2d_p(nc,nvgrid))
      IF (ASSOCIATED(this%phosphorus)) NULLIFY(this%phosphorus)
      this%phosphorus => this%op_ov_hlp_2d_p
      ALLOCATE(this%op_ov_hlp_2d_c13(nc,nvgrid))
      IF (ASSOCIATED(this%carbon13))   NULLIFY(this%carbon13)
      this%carbon13   => this%op_ov_hlp_2d_c13
      ALLOCATE(this%op_ov_hlp_2d_c14(nc,nvgrid))
      IF (ASSOCIATED(this%carbon14))   NULLIFY(this%carbon14)
      this%carbon14   => this%op_ov_hlp_2d_c14
      ALLOCATE(this%op_ov_hlp_2d_n15(nc,nvgrid))
      IF (ASSOCIATED(this%nitrogen15)) NULLIFY(this%nitrogen15)
      this%nitrogen15 => this%op_ov_hlp_2d_n15
    ELSE
      ALLOCATE(this%op_ov_hlp_3d_c(nc,nvgrid,nb))
      IF (ASSOCIATED(this%domain_c))   NULLIFY(this%domain_c)
      this%domain_c   => this%op_ov_hlp_3d_c
      ALLOCATE(this%op_ov_hlp_3d_n(nc,nvgrid,nb))
      IF (ASSOCIATED(this%domain_n))   NULLIFY(this%domain_n)
      this%domain_n   => this%op_ov_hlp_3d_n
      ALLOCATE(this%op_ov_hlp_3d_p(nc,nvgrid,nb))
      IF (ASSOCIATED(this%domain_p))   NULLIFY(this%domain_p)
      this%domain_p   => this%op_ov_hlp_3d_p
      ALLOCATE(this%op_ov_hlp_3d_c13(nc,nvgrid,nb))
      IF (ASSOCIATED(this%domain_c13)) NULLIFY(this%domain_c13)
      this%domain_c13 => this%op_ov_hlp_3d_c13
      ALLOCATE(this%op_ov_hlp_3d_c14(nc,nvgrid,nb))
      IF (ASSOCIATED(this%domain_c14)) NULLIFY(this%domain_c14)
      this%domain_c14 => this%op_ov_hlp_3d_c14
      ALLOCATE(this%op_ov_hlp_3d_n15(nc,nvgrid,nb))
      IF (ASSOCIATED(this%domain_n15)) NULLIFY(this%domain_n15)
      this%domain_n15 => this%op_ov_hlp_3d_n15
    ENDIF

  END SUBROUTINE t_sb_elements_pool_op_ov_hlp_init

  !-----------------------------------------------------------------------------------------------------
  !> set element-variables pointers to the chunk pointer ptr3d(ics:ice,:,iblk) in the element list this%element_list(i)%p
  !> t_sb_elements_pool
  !> this needs to be called for the SB_ pools used in the task routines (when working ON_CHUNK)
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE t_sb_elements_pool_associate_chunk_variable_pointers(this, ics, ice, iblk)

    CLASS(t_sb_elements_pool),  INTENT(inout)   :: this
    INTEGER,                    INTENT(in)      :: ics          !< chunk start
    INTEGER,                    INTENT(in)      :: ice          !< chunk end
    INTEGER,                    INTENT(in)      :: iblk         !< block number

    TYPE(t_subset)              :: sub    ! ON_CHUNK / ON_DOMAIN information
    INTEGER                     :: i      ! looping
    CHARACTER(len=*), PARAMETER :: routine = modname//':t_sb_elements_pool_associate_chunk_variable_pointers'

    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Set element-variables pointers to the chunk.')

    sub = this%Get_subset()

    ! test working ON_CHUNK
    IF (sub%type == ON_DOMAIN) THEN
      CALL finish(routine, ' May be called only in task rouines, i.e., when working ON_CHUNK.')
    END IF

    ! loop over elements in element_list(:) of this element pool
    ! associate REAL element vectors with %ptr2d(ics:ice,:,iblk)
    IF (ALLOCATED(this%element_list)) THEN
      DO i = 1,SIZE(this%element_list)
        SELECT TYPE (element => this%element_list(i)%p)
        CLASS IS (t_jsb_var_real3d)
          SELECT CASE(element%element_id)
          CASE (ELEM_C)
            IF (ASSOCIATED(this%carbon)) NULLIFY(this%carbon)
            this%carbon     => element%ptr3d(ics:ice, :, iblk)
          CASE (ELEM_N)
            IF (ASSOCIATED(this%nitrogen)) NULLIFY(this%nitrogen)
            this%nitrogen   => element%ptr3d(ics:ice, :, iblk)
          CASE (ELEM_P)
            IF (ASSOCIATED(this%phosphorus)) NULLIFY(this%phosphorus)
            this%phosphorus => element%ptr3d(ics:ice, :, iblk)
          CASE (ELEM_C13)
            IF (ASSOCIATED(this%carbon13)) NULLIFY(this%carbon13)
            this%carbon13   => element%ptr3d(ics:ice, :, iblk)
          CASE (ELEM_C14)
            IF (ASSOCIATED(this%carbon14)) NULLIFY(this%carbon14)
            this%carbon14   => element%ptr3d(ics:ice, :, iblk)
          CASE (ELEM_N15)
            IF (ASSOCIATED(this%nitrogen15)) NULLIFY(this%nitrogen15)
            this%nitrogen15 => element%ptr3d(ics:ice, :, iblk)
          END SELECT
        END SELECT
      END DO
    END IF

  END SUBROUTINE t_sb_elements_pool_associate_chunk_variable_pointers

  !-----------------------------------------------------------------------------------------------------
  !> set element-variables pointers to the domain pointer ptr3d(:,:,:) in the element list this%element_list(i)%p
  !! t_sb_elements_pool
  !! this need to be called for the SB_ pools used in the sb_init() routine (or any other *_init() where these pools are used)
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE t_sb_elements_pool_associate_domain_variable_pointers(this)

    CLASS(t_sb_elements_pool),  INTENT(inout) :: this

    TYPE(t_subset)              :: sub    ! ON_CHUNK / ON_DOMAIN information
    INTEGER                     :: i      ! looping
    CHARACTER(len=*), PARAMETER :: routine = modname//':t_sb_elements_pool_associate_domain_variable_pointers'

    sub = this%Get_subset()

    ! test working ON_DOMAIN
    IF (sub%type == ON_CHUNK) THEN
      CALL finish(routine, ' May be called only in the process init rouines, i.e., when working ON_DOMAIN.')
    END IF

    ! loop over elements in element_list(:) of this element pool
    ! associate REAL element vectors with %ptr3d
    IF (ALLOCATED(this%element_list)) THEN
      DO i = 1,SIZE(this%element_list)
        SELECT TYPE (element => this%element_list(i)%p)
        CLASS IS (t_jsb_var_real3d)
          SELECT CASE(element%element_id)
          CASE (ELEM_C)
            IF (ASSOCIATED(this%domain_c)) NULLIFY(this%domain_c)
            this%domain_c    => element%ptr3d
          CASE (ELEM_N)
            IF (ASSOCIATED(this%domain_n)) NULLIFY(this%domain_n)
            this%domain_n    => element%ptr3d
          CASE (ELEM_P)
            IF (ASSOCIATED(this%domain_p))   NULLIFY(this%domain_p)
            this%domain_p    => element%ptr3d
          CASE (ELEM_C13)
            IF (ASSOCIATED(this%domain_c13)) NULLIFY(this%domain_c13)
            this%domain_c13  => element%ptr3d
          CASE (ELEM_C14)
            IF (ASSOCIATED(this%domain_c14)) NULLIFY(this%domain_c14)
            this%domain_c14  => element%ptr3d
          CASE (ELEM_N15)
            IF (ASSOCIATED(this%domain_n15)) NULLIFY(this%domain_n15)
            this%domain_n15  => element%ptr3d
          END SELECT
        END SELECT
      END DO
    END IF

  END SUBROUTINE t_sb_elements_pool_associate_domain_variable_pointers

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
  
    this%element_unit                  = element_unit
    this%contains_elem_2d_ptr_to_chunk = .TRUE.
  
    !> carbon
    IF (One_of(ELEM_C, element_list_loc) > 0) THEN
      CALL this%Add_element('carbon', 'C', ELEM_C, dim=3)
    END IF
    !> nitrogen
    IF (One_of(ELEM_N, element_list_loc) > 0) THEN
      CALL this%Add_element('nitrogen', 'N', ELEM_N, dim=3)
    END IF
    !> phosphorus
    IF (One_of(ELEM_P, element_list_loc) > 0) THEN
      CALL this%Add_element('phosphorus', 'P', ELEM_P, dim=3)
    END IF
    !> carbon13
    IF (One_of(ELEM_C13, element_list_loc) > 0) THEN
      CALL this%Add_element('carbon13', 'C13', ELEM_C13, dim=3)
    END IF
    !> carbon14
    IF (One_of(ELEM_C14, element_list_loc) > 0) THEN
      CALL this%Add_element('carbon14', 'C14', ELEM_C14, dim=3)
    END IF
    !> nitrogen15
    IF (One_of(ELEM_N15, element_list_loc) > 0) THEN
      CALL this%Add_element('nitrogen15', 'N15', ELEM_N15, dim=3)
    END IF

  END SUBROUTINE t_sb_components_pool_init

  !-----------------------------------------------------------------------------------------------------
  !> create and init a hlp object of the SB_ components pool "t_sb_components_pool" for operator overload
  !! "pool_object%Init_op_ov_hlp()"
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE t_sb_components_pool_op_ov_hlp_create_init(this, name, shortname, owner_model_id, subset_type, &
    &                                                   nc, nvgrid, nb, l_elements, element_unit, pool_suffix)

    CLASS(t_sb_components_pool),  INTENT(inout), TARGET :: this
    CHARACTER(LEN=*),             INTENT(in)            :: name               ! pool name       (No spaces !)
    CHARACTER(LEN=*),             INTENT(in)            :: shortname          ! pool short name (No spaces !)
    INTEGER,                      INTENT(in)            :: owner_model_id     ! model ID
    INTEGER,                      INTENT(in)            :: subset_type        ! ON_DOMAIN=1 / ON_CHUNK=2
    INTEGER,                      INTENT(in)            :: nc                 ! number of points in chunk
    INTEGER,                      INTENT(in)            :: nvgrid             ! number of layers of the vertical grid
    INTEGER,                      INTENT(in)            :: nb                 ! number of blocks
    LOGICAL,                      INTENT(in)            :: l_elements         ! IF pool has elements 
    CHARACTER(LEN=*),             INTENT(in)            :: element_unit       ! unit of elements
    CHARACTER(LEN=*),             INTENT(in)            :: pool_suffix        ! add to the name of each elements_pool created here

    CHARACTER(len=*), PARAMETER :: routine = modname//':t_sb_components_pool_op_ov_hlp_create_init'

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
      ALLOCATE(this%op_ov_hlp_2d_c(nc,nvgrid))
      IF (ASSOCIATED(this%carbon))     NULLIFY(this%carbon)
      this%carbon     => this%op_ov_hlp_2d_c
      ALLOCATE(this%op_ov_hlp_2d_n(nc,nvgrid))
      IF (ASSOCIATED(this%nitrogen))   NULLIFY(this%nitrogen)
      this%nitrogen   => this%op_ov_hlp_2d_n
      ALLOCATE(this%op_ov_hlp_2d_p(nc,nvgrid))
      IF (ASSOCIATED(this%phosphorus)) NULLIFY(this%phosphorus)
      this%phosphorus => this%op_ov_hlp_2d_p
      ALLOCATE(this%op_ov_hlp_2d_c13(nc,nvgrid))
      IF (ASSOCIATED(this%carbon13))   NULLIFY(this%carbon13)
      this%carbon13   => this%op_ov_hlp_2d_c13
      ALLOCATE(this%op_ov_hlp_2d_c14(nc,nvgrid))
      IF (ASSOCIATED(this%carbon14))   NULLIFY(this%carbon14)
      this%carbon14   => this%op_ov_hlp_2d_c14
      ALLOCATE(this%op_ov_hlp_2d_n15(nc,nvgrid))
      IF (ASSOCIATED(this%nitrogen15)) NULLIFY(this%nitrogen15)
      this%nitrogen15 => this%op_ov_hlp_2d_n15
    ELSE
      ALLOCATE(this%op_ov_hlp_3d_c(nc,nvgrid,nb))
      IF (ASSOCIATED(this%domain_c))   NULLIFY(this%domain_c)
      this%domain_c   => this%op_ov_hlp_3d_c
      ALLOCATE(this%op_ov_hlp_3d_n(nc,nvgrid,nb))
      IF (ASSOCIATED(this%domain_n))   NULLIFY(this%domain_n)
      this%domain_n   => this%op_ov_hlp_3d_n
      ALLOCATE(this%op_ov_hlp_3d_p(nc,nvgrid,nb))
      IF (ASSOCIATED(this%domain_p))   NULLIFY(this%domain_p)
      this%domain_p   => this%op_ov_hlp_3d_p
      ALLOCATE(this%op_ov_hlp_3d_c13(nc,nvgrid,nb))
      IF (ASSOCIATED(this%domain_c13)) NULLIFY(this%domain_c13)
      this%domain_c13 => this%op_ov_hlp_3d_c13
      ALLOCATE(this%op_ov_hlp_3d_c14(nc,nvgrid,nb))
      IF (ASSOCIATED(this%domain_c14)) NULLIFY(this%domain_c14)
      this%domain_c14 => this%op_ov_hlp_3d_c14
      ALLOCATE(this%op_ov_hlp_3d_n15(nc,nvgrid,nb))
      IF (ASSOCIATED(this%domain_n15)) NULLIFY(this%domain_n15)
      this%domain_n15 => this%op_ov_hlp_3d_n15
    ENDIF

    !> 2.0 element pools: allocate & init & add to components pool
    !!     
    ALLOCATE(this%dom)
    CALL this%dom%Init_op_ov_hlp(               'op_ov_hlp_dom'//pool_suffix, 'op_ov_hlp_dom', &
      &                                         owner_model_id, subset_type, nc, nvgrid, nb, .TRUE., element_unit)
    ALLOCATE(this%dom_assoc)
    CALL this%dom_assoc%Init_op_ov_hlp(         'op_ov_hlp_domass'//pool_suffix, 'op_ov_hlp_domass', &
      &                                         owner_model_id, subset_type, nc, nvgrid, nb, .TRUE., element_unit)
    ALLOCATE(this%soluable_litter)
    CALL this%soluable_litter%Init_op_ov_hlp(   'op_ov_hlp_sollit'//pool_suffix, 'op_ov_hlp_sollit', &
      &                                         owner_model_id, subset_type, nc, nvgrid, nb, .TRUE., element_unit)
    ALLOCATE(this%polymeric_litter)
    CALL this%polymeric_litter%Init_op_ov_hlp(  'op_ov_hlp_pollit'//pool_suffix, 'op_ov_hlp_pollit', &
      &                                         owner_model_id, subset_type, nc, nvgrid, nb, .TRUE., element_unit)
    ALLOCATE(this%woody_litter)
    CALL this%woody_litter%Init_op_ov_hlp(      'op_ov_hlp_wodlit'//pool_suffix, 'op_ov_hlp_wodlit', &
      &                                         owner_model_id, subset_type, nc, nvgrid, nb, .TRUE., element_unit)
    ALLOCATE(this%fungi)
    CALL this%fungi%Init_op_ov_hlp(             'op_ov_hlp_fungi'//pool_suffix, 'op_ov_hlp_fungi', &
      &                                         owner_model_id, subset_type, nc, nvgrid, nb, .TRUE., element_unit)
    ALLOCATE(this%mycorrhiza)
    CALL this%mycorrhiza%Init_op_ov_hlp(        'op_ov_hlp_mycor'//pool_suffix, 'op_ov_hlp_mycor', &
      &                                         owner_model_id, subset_type, nc, nvgrid, nb, .TRUE., element_unit)
    ALLOCATE(this%microbial)
    CALL this%microbial%Init_op_ov_hlp(         'op_ov_hlp_micro'//pool_suffix, 'op_ov_hlp_micro', &
      &                                         owner_model_id, subset_type, nc, nvgrid, nb, .TRUE., element_unit)
    ALLOCATE(this%residue)
    CALL this%residue%Init_op_ov_hlp(           'op_ov_hlp_resid'//pool_suffix, 'op_ov_hlp_resid', &
      &                                         owner_model_id, subset_type, nc, nvgrid, nb, .TRUE., element_unit)
    ALLOCATE(this%residue_assoc)
    CALL this%residue_assoc%Init_op_ov_hlp(     'op_ov_hlp_residass'//pool_suffix, 'op_ov_hlp_residass', &
      &                                         owner_model_id, subset_type, nc, nvgrid, nb, .TRUE., element_unit)

  END SUBROUTINE t_sb_components_pool_op_ov_hlp_create_init

  !-----------------------------------------------------------------------------------------------------
  !> create and init an object of the SB_ tissues pool "t_sb_components_pool"
  !! used for SB_ memory init
  !!   init and add all its objects of the "t_sb_elements_pool" with "pool_object%Init()"
  !!   "object => Create_sb_components_pool()"
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
  !> set element-variables pointers to the chunk pointer ptr3d(ics:ice,:,iblk) in the element list this%element_list(i)%p
  !> t_sb_components_pool
  !> this needs to be called for the SB_ pools used in the task routines (when working ON_CHUNK)
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE t_sb_components_pool_associate_chunk_variable_pointers(this, ics, ice, iblk)

    CLASS(t_sb_components_pool),  INTENT(inout)   :: this
    INTEGER,                      INTENT(in)      :: ics          !< chunk start
    INTEGER,                      INTENT(in)      :: ice          !< chunk end
    INTEGER,                      INTENT(in)      :: iblk         !< block number

    TYPE(t_subset)              :: sub    ! ON_CHUNK / ON_DOMAIN information
    INTEGER                     :: i      ! looping
    CHARACTER(len=*), PARAMETER :: routine = modname//':t_sb_components_pool_associate_chunk_variable_pointers'

    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Set element-variables pointers to the chunk.')

    sub = this%Get_subset()

    ! test working ON_CHUNK
    IF (sub%type == ON_DOMAIN) THEN
      CALL finish(routine, ' May be called only in task rouines, i.e., when working ON_CHUNK.')
    END IF

    ! loop over elements in element_list(:) of this element pool
    ! associate REAL element vectors with %ptr3d(ics:ice,:,iblk)
    IF (ALLOCATED(this%element_list)) THEN
      DO i = 1,SIZE(this%element_list)
        SELECT TYPE (element => this%element_list(i)%p)
        CLASS IS (t_jsb_var_real3d)
          SELECT CASE(element%element_id)
          CASE (ELEM_C)
            IF (ASSOCIATED(this%carbon)) NULLIFY(this%carbon)
            this%carbon     => element%ptr3d(ics:ice, :, iblk)
          CASE (ELEM_N)
            IF (ASSOCIATED(this%nitrogen)) NULLIFY(this%nitrogen)
            this%nitrogen   => element%ptr3d(ics:ice, :, iblk)
          CASE (ELEM_P)
            IF (ASSOCIATED(this%phosphorus)) NULLIFY(this%phosphorus)
            this%phosphorus => element%ptr3d(ics:ice, :, iblk)
          CASE (ELEM_C13)
            IF (ASSOCIATED(this%carbon13)) NULLIFY(this%carbon13)
            this%carbon13   => element%ptr3d(ics:ice, :, iblk)
          CASE (ELEM_C14)
            IF (ASSOCIATED(this%carbon14)) NULLIFY(this%carbon14)
            this%carbon14   => element%ptr3d(ics:ice, :, iblk)
          CASE (ELEM_N15)
            IF (ASSOCIATED(this%nitrogen15)) NULLIFY(this%nitrogen15)
            this%nitrogen15 => element%ptr3d(ics:ice, :, iblk)
          END SELECT
        END SELECT
      END DO
    END IF

    ! associate pointers for all element pools of this components pool
    CALL this%dom%Associate_chunk_var_ptr(ics, ice, iblk)
    CALL this%dom_assoc%Associate_chunk_var_ptr(ics, ice, iblk)
    CALL this%soluable_litter%Associate_chunk_var_ptr(ics, ice, iblk)
    CALL this%polymeric_litter%Associate_chunk_var_ptr(ics, ice, iblk)
    CALL this%woody_litter%Associate_chunk_var_ptr(ics, ice, iblk)
    CALL this%fungi%Associate_chunk_var_ptr(ics, ice, iblk)
    CALL this%mycorrhiza%Associate_chunk_var_ptr(ics, ice, iblk)
    CALL this%microbial%Associate_chunk_var_ptr(ics, ice, iblk)
    CALL this%residue%Associate_chunk_var_ptr(ics, ice, iblk)
    CALL this%residue_assoc%Associate_chunk_var_ptr(ics, ice, iblk)    

  END SUBROUTINE t_sb_components_pool_associate_chunk_variable_pointers

  !-----------------------------------------------------------------------------------------------------
  !> set element-variables pointers to the domain pointer ptr3d(:,:) in the element list this%element_list(i)%p
  !! t_sb_components_pool
  !! this need to be called for the SB_ pools used in the sb_init() routine (or any other *_init() where these pools are used)
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE t_sb_components_pool_associate_domain_variable_pointers(this)

    CLASS(t_sb_components_pool),   INTENT(inout) :: this

    TYPE(t_subset)              :: sub    ! ON_CHUNK / ON_DOMAIN information
    INTEGER                     :: i      ! looping
    CHARACTER(len=*), PARAMETER :: routine = modname//':t_sb_components_pool_associate_domain_variable_pointers'

    sub = this%Get_subset()

    ! test working ON_DOMAIN
    IF (sub%type == ON_CHUNK) THEN
      CALL finish(routine, ' May be called only in the process init rouines, i.e., when working ON_DOMAIN.')
    END IF

    ! loop over elements in element_list(:) of this tissue pool
    ! associate REAL element vectors with %ptr3d
    IF (ALLOCATED(this%element_list)) THEN
      DO i = 1,SIZE(this%element_list)
        SELECT TYPE (element => this%element_list(i)%p)
        CLASS IS (t_jsb_var_real3d)
          SELECT CASE(element%element_id)
          CASE (ELEM_C)
            IF (ASSOCIATED(this%domain_c)) NULLIFY(this%domain_c)
            this%domain_c    => element%ptr3d
          CASE (ELEM_N)
            IF (ASSOCIATED(this%domain_n)) NULLIFY(this%domain_n)
            this%domain_n    => element%ptr3d
          CASE (ELEM_P)
            IF (ASSOCIATED(this%domain_p))   NULLIFY(this%domain_p)
            this%domain_p    => element%ptr3d
          CASE (ELEM_C13)
            IF (ASSOCIATED(this%domain_c13)) NULLIFY(this%domain_c13)
            this%domain_c13  => element%ptr3d
          CASE (ELEM_C14)
            IF (ASSOCIATED(this%domain_c14)) NULLIFY(this%domain_c14)
            this%domain_c14  => element%ptr3d
          CASE (ELEM_N15)
            IF (ASSOCIATED(this%domain_n15)) NULLIFY(this%domain_n15)
            this%domain_n15  => element%ptr3d
          END SELECT
        END SELECT
      END DO
    END IF
    
    ! associate pointers for all element pools of this components pool
    CALL this%dom%Associate_domain_var_ptr()
    CALL this%dom_assoc%Associate_domain_var_ptr()
    CALL this%soluable_litter%Associate_domain_var_ptr()
    CALL this%polymeric_litter%Associate_domain_var_ptr()
    CALL this%woody_litter%Associate_domain_var_ptr()
    CALL this%fungi%Associate_domain_var_ptr()
    CALL this%mycorrhiza%Associate_domain_var_ptr()
    CALL this%microbial%Associate_domain_var_ptr()
    CALL this%residue%Associate_domain_var_ptr()
    CALL this%residue_assoc%Associate_domain_var_ptr()

  END SUBROUTINE t_sb_components_pool_associate_domain_variable_pointers

  !--------------------------------------------------------------------------------------------------------------------------------
  !> operator overload and assignment routines for element pools
  !!  TYPE t_sb_elements_pool
  !! 
  !--------------------------------------------------------------------------------------------------------------------------------

  !--------------------------------------------------------------------------------------------------------------------------------
  ! finalize after operator overload
  ! t_sb_elements_pool
  !--------------------------------------------------------------------------------------------------------------------------------
  SUBROUTINE finalize_sb_elements_pool(this)

    TYPE(t_sb_elements_pool) :: this

    CHARACTER(len=*), PARAMETER :: routine = modname//':finalize_sb_elements_pool'

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

    IF (ALLOCATED(this%op_ov_hlp_2d_c))     DEALLOCATE(this%op_ov_hlp_2d_c)
    IF (ALLOCATED(this%op_ov_hlp_2d_n))     DEALLOCATE(this%op_ov_hlp_2d_n)
    IF (ALLOCATED(this%op_ov_hlp_2d_p))     DEALLOCATE(this%op_ov_hlp_2d_p)
    IF (ALLOCATED(this%op_ov_hlp_2d_c13))   DEALLOCATE(this%op_ov_hlp_2d_c13)
    IF (ALLOCATED(this%op_ov_hlp_2d_c14))   DEALLOCATE(this%op_ov_hlp_2d_c14)
    IF (ALLOCATED(this%op_ov_hlp_2d_n15))   DEALLOCATE(this%op_ov_hlp_2d_n15)
    IF (ALLOCATED(this%op_ov_hlp_3d_c))     DEALLOCATE(this%op_ov_hlp_3d_c)
    IF (ALLOCATED(this%op_ov_hlp_3d_n))     DEALLOCATE(this%op_ov_hlp_3d_n)
    IF (ALLOCATED(this%op_ov_hlp_3d_p))     DEALLOCATE(this%op_ov_hlp_3d_p)
    IF (ALLOCATED(this%op_ov_hlp_3d_c13))   DEALLOCATE(this%op_ov_hlp_3d_c13)
    IF (ALLOCATED(this%op_ov_hlp_3d_c14))   DEALLOCATE(this%op_ov_hlp_3d_c14)
    IF (ALLOCATED(this%op_ov_hlp_3d_n15))   DEALLOCATE(this%op_ov_hlp_3d_n15)
    
    IF (ALLOCATED(this%shortname))          DEALLOCATE(this%shortname)
    IF (ALLOCATED(this%name))               DEALLOCATE(this%name)
    ! not allocated: this%path
    IF (ALLOCATED(this%element_unit))       DEALLOCATE(this%element_unit)
    IF (ASSOCIATED(this%pool_list))         DEALLOCATE(this%pool_list)
    IF (ASSOCIATED(this%leaf_list))         DEALLOCATE(this%leaf_list)
    IF (ALLOCATED(this%element_list))       DEALLOCATE(this%element_list)
    IF (ALLOCATED(this%var_list))           DEALLOCATE(this%var_list)    
    ! this%temporary of t_hermetic is PRIVATE

  END SUBROUTINE finalize_sb_elements_pool
  
  !--------------------------------------------------------------------------------------------------------------------------------
  ! = assign
  !--------------------------------------------------------------------------------------------------------------------------------
  !--------------------------------------------------------------------------------------------------------------------------------
  ! t_sb_elements_pool = scalar | ON_CHUNK & ON_DOMAIN
  !--------------------------------------------------------------------------------------------------------------------------------
  SUBROUTINE Assign_scalar_to_sb_elem_p(lhs, rhs)

    CLASS(t_sb_elements_pool), INTENT(inout) :: lhs
    REAL(wp),                  INTENT(in)    :: rhs

    TYPE(t_subset)              :: sub
    CHARACTER(len=*), PARAMETER :: routine = modname//':Assign_scalar_to_sb_elem_p'

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

  END SUBROUTINE Assign_scalar_to_sb_elem_p
  !--------------------------------------------------------------------------------------------------------------------------------
  ! t_sb_elements_pool = field2d(:,:) | ON_CHUNK
  !--------------------------------------------------------------------------------------------------------------------------------
  SUBROUTINE Assign_field2d_to_sb_elem_p(lhs, rhs)
  
    CLASS(t_sb_elements_pool), INTENT(inout) :: lhs
    REAL(wp),                  INTENT(in)    :: rhs(:,:)

    TYPE(t_subset)              :: sub
    CHARACTER(len=*), PARAMETER :: routine = modname//':Assign_field2d_to_sb_elem_p'
    
    sub = lhs%Get_subset()
    
    IF (sub%type == ON_CHUNK) THEN
      lhs%carbon     = rhs(:,:)
      lhs%nitrogen   = rhs(:,:)
      lhs%phosphorus = rhs(:,:)
      lhs%carbon13   = rhs(:,:)
      lhs%carbon14   = rhs(:,:)
      lhs%nitrogen15 = rhs(:,:)
    ELSE
      CALL finish(routine, lhs%shortname//' Dimension mismatch, cannot assign 2d-field to jsb_var real3d on domain')
    END IF
  
  END SUBROUTINE Assign_field2d_to_sb_elem_p
  !--------------------------------------------------------------------------------------------------------------------------------
  ! t_sb_elements_pool = field3d(:,:,:) | ON_CHUNK & ON_DOMAIN | field3d(:,:,1) & field3d(:,:,:)
  !--------------------------------------------------------------------------------------------------------------------------------
  SUBROUTINE Assign_field3d_to_sb_elem_p(lhs, rhs)
  
    CLASS(t_sb_elements_pool), INTENT(inout) :: lhs
    REAL(wp),                  INTENT(in)    :: rhs(:,:,:)

    TYPE(t_subset)              :: sub
    CHARACTER(len=*), PARAMETER :: routine = modname//':Assign_field3d_to_sb_elem_p'

    sub = lhs%Get_subset()

    IF (sub%type == ON_CHUNK) THEN
      IF (SIZE(rhs,3) > 1) CALL finish(routine, 'Dimension mismatch, field3d with nblks > 1; assignment to REAL pool%element(:)')
      lhs%carbon     = rhs(:,:,1)
      lhs%nitrogen   = rhs(:,:,1)
      lhs%phosphorus = rhs(:,:,1)
      lhs%carbon13   = rhs(:,:,1)
      lhs%carbon14   = rhs(:,:,1)
      lhs%nitrogen15 = rhs(:,:,1)
    ELSE
      lhs%domain_c   = rhs(:,:,:)
      lhs%domain_n   = rhs(:,:,:)
      lhs%domain_p   = rhs(:,:,:)
      lhs%domain_c13 = rhs(:,:,:)
      lhs%domain_c14 = rhs(:,:,:)
      lhs%domain_n15 = rhs(:,:,:)
    END IF
  
  END SUBROUTINE Assign_field3d_to_sb_elem_p
  !--------------------------------------------------------------------------------------------------------------------------------
  ! t_sb_elements_pool = t_sb_elements_pool | ON_CHUNK & ON_DOMAIN
  !--------------------------------------------------------------------------------------------------------------------------------
  SUBROUTINE Assign_sb_elem_p_to_sb_elem_p(lhs, rhs)

    CLASS(t_sb_elements_pool), INTENT(inout) :: lhs
    CLASS(t_sb_elements_pool), INTENT(in)    :: rhs

    TYPE(t_subset)              :: sub
    CHARACTER(len=*), PARAMETER :: routine = modname//':Assign_sb_elem_p_to_sb_elem_p'
    
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

  END SUBROUTINE Assign_sb_elem_p_to_sb_elem_p

  !--------------------------------------------------------------------------------------------------------------------------------
  ! + addition
  !--------------------------------------------------------------------------------------------------------------------------------
  !--------------------------------------------------------------------------------------------------------------------------------
  ! t_sb_elements_pool + scalar        | ON_CHUNK & ON_DOMAIN
  !--------------------------------------------------------------------------------------------------------------------------------
  FUNCTION Add_sb_elem_p_and_scalar(struct1, scalar) RESULT(structOut)
  
    CLASS(t_sb_elements_pool),     INTENT(in)  :: struct1
    REAL(wp),                      INTENT(in)  :: scalar
    TYPE(t_sb_elements_pool)                   :: structOut

    TYPE(t_subset)              :: sub  
    INTEGER                     :: nvgrid
    CHARACTER(len=*), PARAMETER :: routine = modname//':Add_sb_elem_p_and_scalar'

    sub = struct1%Get_subset()
    IF (sub%type /= ON_CHUNK) THEN  ! get domain dimension not available via Get_subset() during process init
      sub%nc = SIZE(struct1%domain_c(:,:,:), 1)
      nvgrid = SIZE(struct1%domain_c(:,:,:), 2)
      sub%nb = SIZE(struct1%domain_c(:,:,:), 3)
    ELSE
      nvgrid = SIZE(struct1%carbon(:,:), 2)
    ENDIF
    CALL structOut%Init_op_ov_hlp('op_ov_hlp', 'op_ov_hlp', struct1%owner_model_id, &
      &                           sub%type, sub%nc, nvgrid, sub%nb, .TRUE., TRIM(struct1%element_unit))

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

  END FUNCTION Add_sb_elem_p_and_scalar
  !--------------------------------------------------------------------------------------------------------------------------------
  ! scalar + t_sb_elements_pool        | ON_CHUNK & ON_DOMAIN
  !--------------------------------------------------------------------------------------------------------------------------------
  FUNCTION Add_scalar_and_sb_elem_p(scalar, struct1) RESULT(structOut)
  
    REAL(wp),                      INTENT(in)  :: scalar
    CLASS(t_sb_elements_pool),     INTENT(in)  :: struct1
    TYPE(t_sb_elements_pool)                   :: structOut

    TYPE(t_subset)              :: sub  
    INTEGER                     :: nvgrid
    CHARACTER(len=*), PARAMETER :: routine = modname//':Add_scalar_and_sb_elem_p'

    sub = struct1%Get_subset()
    IF (sub%type /= ON_CHUNK) THEN  ! get domain dimension not available via Get_subset() during process init
      sub%nc = SIZE(struct1%domain_c(:,:,:), 1)
      nvgrid = SIZE(struct1%domain_c(:,:,:), 2)
      sub%nb = SIZE(struct1%domain_c(:,:,:), 3)
    ELSE
      nvgrid = SIZE(struct1%carbon(:,:), 2)
    ENDIF
    CALL structOut%Init_op_ov_hlp('op_ov_hlp', 'op_ov_hlp', struct1%owner_model_id, &
      &                           sub%type, sub%nc, nvgrid, sub%nb, .TRUE., TRIM(struct1%element_unit))

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

  END FUNCTION Add_scalar_and_sb_elem_p
  !--------------------------------------------------------------------------------------------------------------------------------
  ! t_sb_elements_pool + field2d(:,:)        | ON_CHUNK
  !--------------------------------------------------------------------------------------------------------------------------------
  FUNCTION Add_sb_elem_p_and_field2d(struct1, field2d) RESULT(structOut)
  
    CLASS(t_sb_elements_pool),     INTENT(in)  :: struct1
    REAL(wp),                      INTENT(in)  :: field2d(:,:)
    TYPE(t_sb_elements_pool)                   :: structOut

    TYPE(t_subset)              :: sub  
    INTEGER                     :: nvgrid
    CHARACTER(len=*), PARAMETER :: routine = modname//':Add_sb_elem_p_and_field2d'

    sub = struct1%Get_subset()
    IF (sub%type /= ON_CHUNK) THEN  ! get domain dimension not available via Get_subset() during process init
      sub%nc = SIZE(struct1%domain_c(:,:,:), 1)
      nvgrid = SIZE(struct1%domain_c(:,:,:), 2)
      sub%nb = SIZE(struct1%domain_c(:,:,:), 3)
    ELSE
      nvgrid = SIZE(struct1%carbon(:,:), 2)
    ENDIF
    CALL structOut%Init_op_ov_hlp('op_ov_hlp', 'op_ov_hlp', struct1%owner_model_id, &
      &                           sub%type, sub%nc, nvgrid, sub%nb, .TRUE., TRIM(struct1%element_unit))

    ! run calculation
    IF (sub%type == ON_CHUNK) THEN
      structOut%carbon          = struct1%carbon      + field2d(:,:)
      structOut%nitrogen        = struct1%nitrogen    + field2d(:,:)
      structOut%phosphorus      = struct1%phosphorus  + field2d(:,:)
      structOut%carbon13        = struct1%carbon13    + field2d(:,:)
      structOut%carbon14        = struct1%carbon14    + field2d(:,:)
      structOut%nitrogen15      = struct1%nitrogen15  + field2d(:,:)
    ELSE
      CALL finish(routine, 'Dimension mismatch, cannot calc with 2d-field and jsb_var real3d on domain')
    END IF

  END FUNCTION Add_sb_elem_p_and_field2d
  !--------------------------------------------------------------------------------------------------------------------------------
  ! field2d(:,:) + t_sb_elements_pool        | ON_CHUNK
  !--------------------------------------------------------------------------------------------------------------------------------
  FUNCTION Add_field2d_and_sb_elem_p(field2d, struct1) RESULT(structOut)
  
    REAL(wp),                      INTENT(in)  :: field2d(:,:)
    CLASS(t_sb_elements_pool),     INTENT(in)  :: struct1
    TYPE(t_sb_elements_pool)                   :: structOut

    TYPE(t_subset)              :: sub  
    INTEGER                     :: nvgrid
    CHARACTER(len=*), PARAMETER :: routine = modname//':Add_field2d_and_sb_elem_p'

    sub = struct1%Get_subset()
    IF (sub%type /= ON_CHUNK) THEN  ! get domain dimension not available via Get_subset() during process init
      sub%nc = SIZE(struct1%domain_c(:,:,:), 1)
      nvgrid = SIZE(struct1%domain_c(:,:,:), 2)
      sub%nb = SIZE(struct1%domain_c(:,:,:), 3)
    ELSE
      nvgrid = SIZE(struct1%carbon(:,:), 2)
    ENDIF
    CALL structOut%Init_op_ov_hlp('op_ov_hlp', 'op_ov_hlp', struct1%owner_model_id, &
      &                           sub%type, sub%nc, nvgrid, sub%nb, .TRUE., TRIM(struct1%element_unit))

    ! run calculation
    IF (sub%type == ON_CHUNK) THEN
      structOut%carbon          = field2d(:,:) + struct1%carbon    
      structOut%nitrogen        = field2d(:,:) + struct1%nitrogen  
      structOut%phosphorus      = field2d(:,:) + struct1%phosphorus
      structOut%carbon13        = field2d(:,:) + struct1%carbon13  
      structOut%carbon14        = field2d(:,:) + struct1%carbon14  
      structOut%nitrogen15      = field2d(:,:) + struct1%nitrogen15
    ELSE
      CALL finish(routine, 'Dimension mismatch, cannot calc with 2d-field and jsb_var real3d on domain')
    END IF

  END FUNCTION Add_field2d_and_sb_elem_p
  !--------------------------------------------------------------------------------------------------------------------------------
  ! t_sb_elements_pool + field3d(:,:,:)        | ON_CHUNK & ON_DOMAIN | field3d(:,:,1) & field3d(:,:,:)
  !--------------------------------------------------------------------------------------------------------------------------------
  FUNCTION Add_sb_elem_p_and_field3d(struct1, field3d) RESULT(structOut)
  
    CLASS(t_sb_elements_pool),     INTENT(in)  :: struct1
    REAL(wp),                      INTENT(in)  :: field3d(:,:,:)
    TYPE(t_sb_elements_pool)                   :: structOut

    TYPE(t_subset)              :: sub  
    INTEGER                     :: nvgrid
    CHARACTER(len=*), PARAMETER :: routine = modname//':Add_sb_elem_p_and_field3d'

    sub = struct1%Get_subset()
    IF (sub%type /= ON_CHUNK) THEN  ! get domain dimension not available via Get_subset() during process init
      sub%nc = SIZE(struct1%domain_c(:,:,:), 1)
      nvgrid = SIZE(struct1%domain_c(:,:,:), 2)
      sub%nb = SIZE(struct1%domain_c(:,:,:), 3)
    ELSE
      nvgrid = SIZE(struct1%carbon(:,:), 2)
    ENDIF
    CALL structOut%Init_op_ov_hlp('op_ov_hlp', 'op_ov_hlp', struct1%owner_model_id, &
      &                           sub%type, sub%nc, nvgrid, sub%nb, .TRUE., TRIM(struct1%element_unit))

    ! run calculation
    IF (sub%type == ON_CHUNK) THEN
      IF (SIZE(field3d,3) > 1) CALL finish(routine, 'Dimension mismatch, field3d with nblks > 1; calc with REAL pool%element(:)')
      structOut%carbon          = struct1%carbon      + field3d(:,:,1)
      structOut%nitrogen        = struct1%nitrogen    + field3d(:,:,1)
      structOut%phosphorus      = struct1%phosphorus  + field3d(:,:,1)
      structOut%carbon13        = struct1%carbon13    + field3d(:,:,1)
      structOut%carbon14        = struct1%carbon14    + field3d(:,:,1)
      structOut%nitrogen15      = struct1%nitrogen15  + field3d(:,:,1)
    ELSE
      structOut%domain_c        = struct1%domain_c    + field3d(:,:,:)
      structOut%domain_n        = struct1%domain_n    + field3d(:,:,:)
      structOut%domain_p        = struct1%domain_p    + field3d(:,:,:)
      structOut%domain_c13      = struct1%domain_c13  + field3d(:,:,:)
      structOut%domain_c14      = struct1%domain_c14  + field3d(:,:,:)
      structOut%domain_n15      = struct1%domain_n15  + field3d(:,:,:)
    END IF

  END FUNCTION Add_sb_elem_p_and_field3d
  !--------------------------------------------------------------------------------------------------------------------------------
  ! t_sb_elements_pool + t_sb_elements_pool | ON_CHUNK & ON_DOMAIN
  !--------------------------------------------------------------------------------------------------------------------------------
  FUNCTION Add_sb_elem_p_and_sb_elem_p(struct1, struct2) RESULT(structOut)
  
    CLASS(t_sb_elements_pool),     INTENT(in)  :: struct1, struct2
    TYPE(t_sb_elements_pool)                   :: structOut

    TYPE(t_subset)              :: sub  
    INTEGER                     :: nvgrid
    CHARACTER(len=*), PARAMETER :: routine = modname//':Add_sb_elem_p_and_sb_elem_p'

    sub = struct1%Get_subset()
    IF (sub%type /= ON_CHUNK) THEN  ! get domain dimension not available via Get_subset() during process init
      sub%nc = SIZE(struct1%domain_c(:,:,:), 1)
      nvgrid = SIZE(struct1%domain_c(:,:,:), 2)
      sub%nb = SIZE(struct1%domain_c(:,:,:), 3)
    ELSE
      nvgrid = SIZE(struct1%carbon(:,:), 2)
    ENDIF
    CALL structOut%Init_op_ov_hlp('op_ov_hlp', 'op_ov_hlp', struct1%owner_model_id, &
      &                           sub%type, sub%nc, nvgrid, sub%nb, .TRUE., TRIM(struct1%element_unit))

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

  END FUNCTION Add_sb_elem_p_and_sb_elem_p

  !--------------------------------------------------------------------------------------------------------------------------------
  ! - minus
  !--------------------------------------------------------------------------------------------------------------------------------
  !--------------------------------------------------------------------------------------------------------------------------------
  ! t_sb_elements_pool - scalar        | ON_CHUNK & ON_DOMAIN
  !--------------------------------------------------------------------------------------------------------------------------------
  FUNCTION Subtract_sb_elem_p_minus_scalar(struct1, scalar) RESULT(structOut)
  
    CLASS(t_sb_elements_pool),     INTENT(in)  :: struct1
    REAL(wp),                      INTENT(in)  :: scalar
    TYPE(t_sb_elements_pool)                   :: structOut

    TYPE(t_subset)              :: sub  
    INTEGER                     :: nvgrid
    CHARACTER(len=*), PARAMETER :: routine = modname//':Subtract_sb_elem_p_minus_scalar'

    sub = struct1%Get_subset()
    IF (sub%type /= ON_CHUNK) THEN  ! get domain dimension not available via Get_subset() during process init
      sub%nc = SIZE(struct1%domain_c(:,:,:), 1)
      nvgrid = SIZE(struct1%domain_c(:,:,:), 2)
      sub%nb = SIZE(struct1%domain_c(:,:,:), 3)
    ELSE
      nvgrid = SIZE(struct1%carbon(:,:), 2)
    ENDIF
    CALL structOut%Init_op_ov_hlp('op_ov_hlp', 'op_ov_hlp', struct1%owner_model_id, &
      &                           sub%type, sub%nc, nvgrid, sub%nb, .TRUE., TRIM(struct1%element_unit))

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

  END FUNCTION Subtract_sb_elem_p_minus_scalar
  !--------------------------------------------------------------------------------------------------------------------------------
  ! t_sb_elements_pool - field2d(:,:)        | ON_CHUNK
  !--------------------------------------------------------------------------------------------------------------------------------
  FUNCTION Subtract_sb_elem_p_minus_field2d(struct1, field2d) RESULT(structOut)
  
    CLASS(t_sb_elements_pool),     INTENT(in)  :: struct1
    REAL(wp),                      INTENT(in)  :: field2d(:,:)
    TYPE(t_sb_elements_pool)                   :: structOut

    TYPE(t_subset)              :: sub  
    INTEGER                     :: nvgrid
    CHARACTER(len=*), PARAMETER :: routine = modname//':Subtract_sb_elem_p_minus_field2d'

    sub = struct1%Get_subset()
    IF (sub%type /= ON_CHUNK) THEN  ! get domain dimension not available via Get_subset() during process init
      sub%nc = SIZE(struct1%domain_c(:,:,:), 1)
      nvgrid = SIZE(struct1%domain_c(:,:,:), 2)
      sub%nb = SIZE(struct1%domain_c(:,:,:), 3)
    ELSE
      nvgrid = SIZE(struct1%carbon(:,:), 2)
    ENDIF
    CALL structOut%Init_op_ov_hlp('op_ov_hlp', 'op_ov_hlp', struct1%owner_model_id, &
      &                           sub%type, sub%nc, nvgrid, sub%nb, .TRUE., TRIM(struct1%element_unit))

    ! run calculation
    IF (sub%type == ON_CHUNK) THEN
      structOut%carbon          = struct1%carbon      - field2d(:,:)
      structOut%nitrogen        = struct1%nitrogen    - field2d(:,:)
      structOut%phosphorus      = struct1%phosphorus  - field2d(:,:)
      structOut%carbon13        = struct1%carbon13    - field2d(:,:)
      structOut%carbon14        = struct1%carbon14    - field2d(:,:)
      structOut%nitrogen15      = struct1%nitrogen15  - field2d(:,:)
    ELSE
      CALL finish(routine, 'Dimension mismatch, cannot calc with 2d-field and jsb_var real3d on domain')
    END IF

  END FUNCTION Subtract_sb_elem_p_minus_field2d
  !--------------------------------------------------------------------------------------------------------------------------------
  ! t_sb_elements_pool - field3d(:,:,:)        | ON_CHUNK & ON_DOMAIN | field3d(:,:,1) & field3d(:,:,:)
  !--------------------------------------------------------------------------------------------------------------------------------
  FUNCTION Subtract_sb_elem_p_minus_field3d(struct1, field3d) RESULT(structOut)
  
    CLASS(t_sb_elements_pool),     INTENT(in)  :: struct1
    REAL(wp),                      INTENT(in)  :: field3d(:,:,:)
    TYPE(t_sb_elements_pool)                   :: structOut

    TYPE(t_subset)              :: sub  
    INTEGER                     :: nvgrid
    CHARACTER(len=*), PARAMETER :: routine = modname//':Subtract_sb_elem_p_minus_field3d'

    sub = struct1%Get_subset()
    IF (sub%type /= ON_CHUNK) THEN  ! get domain dimension not available via Get_subset() during process init
      sub%nc = SIZE(struct1%domain_c(:,:,:), 1)
      nvgrid = SIZE(struct1%domain_c(:,:,:), 2)
      sub%nb = SIZE(struct1%domain_c(:,:,:), 3)
    ELSE
      nvgrid = SIZE(struct1%carbon(:,:), 2)
    ENDIF
    CALL structOut%Init_op_ov_hlp('op_ov_hlp', 'op_ov_hlp', struct1%owner_model_id, &
      &                           sub%type, sub%nc, nvgrid, sub%nb, .TRUE., TRIM(struct1%element_unit))

    ! run calculation
    IF (sub%type == ON_CHUNK) THEN
      IF (SIZE(field3d,3) > 1) CALL finish(routine, 'Dimension mismatch, field3d with nblks > 1; calc with REAL pool%element(:)')
      structOut%carbon          = struct1%carbon      - field3d(:,:,1)
      structOut%nitrogen        = struct1%nitrogen    - field3d(:,:,1)
      structOut%phosphorus      = struct1%phosphorus  - field3d(:,:,1)
      structOut%carbon13        = struct1%carbon13    - field3d(:,:,1)
      structOut%carbon14        = struct1%carbon14    - field3d(:,:,1)
      structOut%nitrogen15      = struct1%nitrogen15  - field3d(:,:,1)
    ELSE
      structOut%domain_c        = struct1%domain_c    - field3d(:,:,:)
      structOut%domain_n        = struct1%domain_n    - field3d(:,:,:)
      structOut%domain_p        = struct1%domain_p    - field3d(:,:,:)
      structOut%domain_c13      = struct1%domain_c13  - field3d(:,:,:)
      structOut%domain_c14      = struct1%domain_c14  - field3d(:,:,:)
      structOut%domain_n15      = struct1%domain_n15  - field3d(:,:,:)
    END IF

  END FUNCTION Subtract_sb_elem_p_minus_field3d
  !--------------------------------------------------------------------------------------------------------------------------------
  ! t_sb_elements_pool - t_sb_elements_pool | ON_CHUNK & ON_DOMAIN
  !--------------------------------------------------------------------------------------------------------------------------------
  FUNCTION Subtract_sb_elem_p_minus_sb_elem_p(struct1, struct2) RESULT(structOut)
  
    CLASS(t_sb_elements_pool),     INTENT(in)  :: struct1, struct2
    TYPE(t_sb_elements_pool)                   :: structOut

    TYPE(t_subset)              :: sub  
    INTEGER                     :: nvgrid
    CHARACTER(len=*), PARAMETER :: routine = modname//':Subtract_sb_elem_p_minus_sb_elem_p'

    sub = struct1%Get_subset()
    IF (sub%type /= ON_CHUNK) THEN  ! get domain dimension not available via Get_subset() during process init
      sub%nc = SIZE(struct1%domain_c(:,:,:), 1)
      nvgrid = SIZE(struct1%domain_c(:,:,:), 2)
      sub%nb = SIZE(struct1%domain_c(:,:,:), 3)
    ELSE
      nvgrid = SIZE(struct1%carbon(:,:), 2)
    ENDIF
    CALL structOut%Init_op_ov_hlp('op_ov_hlp', 'op_ov_hlp', struct1%owner_model_id, &
      &                           sub%type, sub%nc, nvgrid, sub%nb, .TRUE., TRIM(struct1%element_unit))

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

  END FUNCTION Subtract_sb_elem_p_minus_sb_elem_p

  !--------------------------------------------------------------------------------------------------------------------------------
  ! * multiply
  !--------------------------------------------------------------------------------------------------------------------------------
  !--------------------------------------------------------------------------------------------------------------------------------
  ! t_sb_elements_pool * scalar        | ON_CHUNK & ON_DOMAIN
  !--------------------------------------------------------------------------------------------------------------------------------
  FUNCTION Multiply_sb_elem_p_with_scalar(struct1, scalar) RESULT(structOut)
  
    CLASS(t_sb_elements_pool),     INTENT(in)  :: struct1
    REAL(wp),                      INTENT(in)  :: scalar
    TYPE(t_sb_elements_pool)                   :: structOut

    TYPE(t_subset)              :: sub  
    INTEGER                     :: nvgrid
    CHARACTER(len=*), PARAMETER :: routine = modname//':Multiply_sb_elem_p_with_scalar'

    sub = struct1%Get_subset()
    IF (sub%type /= ON_CHUNK) THEN  ! get domain dimension not available via Get_subset() during process init
      sub%nc = SIZE(struct1%domain_c(:,:,:), 1)
      nvgrid = SIZE(struct1%domain_c(:,:,:), 2)
      sub%nb = SIZE(struct1%domain_c(:,:,:), 3)
    ELSE
      nvgrid = SIZE(struct1%carbon(:,:), 2)
    ENDIF
    CALL structOut%Init_op_ov_hlp('op_ov_hlp', 'op_ov_hlp', struct1%owner_model_id, &
      &                           sub%type, sub%nc, nvgrid, sub%nb, .TRUE., TRIM(struct1%element_unit))

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

  END FUNCTION Multiply_sb_elem_p_with_scalar
  !--------------------------------------------------------------------------------------------------------------------------------
  ! scalar * t_sb_elements_pool        | ON_CHUNK & ON_DOMAIN
  !--------------------------------------------------------------------------------------------------------------------------------
  FUNCTION Multiply_scalar_with_sb_elem_p(scalar, struct1) RESULT(structOut)
  
    REAL(wp),                      INTENT(in)  :: scalar
    CLASS(t_sb_elements_pool),     INTENT(in)  :: struct1
    TYPE(t_sb_elements_pool)                   :: structOut

    TYPE(t_subset)              :: sub  
    INTEGER                     :: nvgrid
    CHARACTER(len=*), PARAMETER :: routine = modname//':Multiply_scalar_with_sb_elem_p'

    sub = struct1%Get_subset()
    IF (sub%type /= ON_CHUNK) THEN  ! get domain dimension not available via Get_subset() during process init
      sub%nc = SIZE(struct1%domain_c(:,:,:), 1)
      nvgrid = SIZE(struct1%domain_c(:,:,:), 2)
      sub%nb = SIZE(struct1%domain_c(:,:,:), 3)
    ELSE
      nvgrid = SIZE(struct1%carbon(:,:), 2)
    ENDIF
    CALL structOut%Init_op_ov_hlp('op_ov_hlp', 'op_ov_hlp', struct1%owner_model_id, &
      &                           sub%type, sub%nc, nvgrid, sub%nb, .TRUE., TRIM(struct1%element_unit))

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

  END FUNCTION Multiply_scalar_with_sb_elem_p
  !--------------------------------------------------------------------------------------------------------------------------------
  ! t_sb_elements_pool * field2d(:,:)        | ON_CHUNK
  !--------------------------------------------------------------------------------------------------------------------------------
  FUNCTION Multiply_sb_elem_p_with_field2d(struct1, field2d) RESULT(structOut)
  
    CLASS(t_sb_elements_pool),     INTENT(in)  :: struct1
    REAL(wp),                      INTENT(in)  :: field2d(:,:)
    TYPE(t_sb_elements_pool)                   :: structOut

    TYPE(t_subset)              :: sub  
    INTEGER                     :: nvgrid
    CHARACTER(len=*), PARAMETER :: routine = modname//':Multiply_sb_elem_p_with_field2d'

    sub = struct1%Get_subset()
    IF (sub%type /= ON_CHUNK) THEN  ! get domain dimension not available via Get_subset() during process init
      sub%nc = SIZE(struct1%domain_c(:,:,:), 1)
      nvgrid = SIZE(struct1%domain_c(:,:,:), 2)
      sub%nb = SIZE(struct1%domain_c(:,:,:), 3)
    ELSE
      nvgrid = SIZE(struct1%carbon(:,:), 2)
    ENDIF
    CALL structOut%Init_op_ov_hlp('op_ov_hlp', 'op_ov_hlp', struct1%owner_model_id, &
      &                           sub%type, sub%nc, nvgrid, sub%nb, .TRUE., TRIM(struct1%element_unit))

    ! run calculation
    IF (sub%type == ON_CHUNK) THEN
      structOut%carbon          = struct1%carbon      * field2d(:,:)
      structOut%nitrogen        = struct1%nitrogen    * field2d(:,:)
      structOut%phosphorus      = struct1%phosphorus  * field2d(:,:)
      structOut%carbon13        = struct1%carbon13    * field2d(:,:)
      structOut%carbon14        = struct1%carbon14    * field2d(:,:)
      structOut%nitrogen15      = struct1%nitrogen15  * field2d(:,:)
    ELSE
      CALL finish(routine, 'Dimension mismatch, cannot calc with 2d-field and jsb_var real3d on domain')
    END IF

  END FUNCTION Multiply_sb_elem_p_with_field2d
  !--------------------------------------------------------------------------------------------------------------------------------
  ! field2d(:,:) * t_sb_elements_pool        | ON_CHUNK
  !--------------------------------------------------------------------------------------------------------------------------------
  FUNCTION Multiply_field2d_with_sb_elem_p(field2d, struct1) RESULT(structOut)
  
    REAL(wp),                      INTENT(in)  :: field2d(:,:)
    CLASS(t_sb_elements_pool),     INTENT(in)  :: struct1
    TYPE(t_sb_elements_pool)                   :: structOut

    TYPE(t_subset)              :: sub  
    INTEGER                     :: nvgrid
    CHARACTER(len=*), PARAMETER :: routine = modname//':Multiply_field2d_with_sb_elem_p'

    sub = struct1%Get_subset()
    IF (sub%type /= ON_CHUNK) THEN  ! get domain dimension not available via Get_subset() during process init
      sub%nc = SIZE(struct1%domain_c(:,:,:), 1)
      nvgrid = SIZE(struct1%domain_c(:,:,:), 2)
      sub%nb = SIZE(struct1%domain_c(:,:,:), 3)
    ELSE
      nvgrid = SIZE(struct1%carbon(:,:), 2)
    ENDIF
    CALL structOut%Init_op_ov_hlp('op_ov_hlp', 'op_ov_hlp', struct1%owner_model_id, &
      &                           sub%type, sub%nc, nvgrid, sub%nb, .TRUE., TRIM(struct1%element_unit))

    ! run calculation
    IF (sub%type == ON_CHUNK) THEN
      structOut%carbon          = field2d(:,:) * struct1%carbon    
      structOut%nitrogen        = field2d(:,:) * struct1%nitrogen  
      structOut%phosphorus      = field2d(:,:) * struct1%phosphorus
      structOut%carbon13        = field2d(:,:) * struct1%carbon13  
      structOut%carbon14        = field2d(:,:) * struct1%carbon14  
      structOut%nitrogen15      = field2d(:,:) * struct1%nitrogen15
    ELSE
      CALL finish(routine, 'Dimension mismatch, cannot calc with 2d-field and jsb_var real3d on domain')
    END IF

  END FUNCTION Multiply_field2d_with_sb_elem_p
  !--------------------------------------------------------------------------------------------------------------------------------
  ! t_sb_elements_pool * field3d(:,:,:)        | ON_CHUNK & ON_DOMAIN | field3d(:,:,1) & field3d(:,:,:)
  !--------------------------------------------------------------------------------------------------------------------------------
  FUNCTION Multiply_sb_elem_p_with_field3d(struct1, field3d) RESULT(structOut)
  
    CLASS(t_sb_elements_pool),     INTENT(in)  :: struct1
    REAL(wp),                      INTENT(in)  :: field3d(:,:,:)
    TYPE(t_sb_elements_pool)                   :: structOut

    TYPE(t_subset)              :: sub  
    INTEGER                     :: nvgrid
    CHARACTER(len=*), PARAMETER :: routine = modname//':Multiply_sb_elem_p_with_field3d'

    sub = struct1%Get_subset()
    IF (sub%type /= ON_CHUNK) THEN  ! get domain dimension not available via Get_subset() during process init
      sub%nc = SIZE(struct1%domain_c(:,:,:), 1)
      nvgrid = SIZE(struct1%domain_c(:,:,:), 2)
      sub%nb = SIZE(struct1%domain_c(:,:,:), 3)
    ELSE
      nvgrid = SIZE(struct1%carbon(:,:), 2)
    ENDIF
    CALL structOut%Init_op_ov_hlp('op_ov_hlp', 'op_ov_hlp', struct1%owner_model_id, &
      &                           sub%type, sub%nc, nvgrid, sub%nb, .TRUE., TRIM(struct1%element_unit))

    ! run calculation
    IF (sub%type == ON_CHUNK) THEN
      IF (SIZE(field3d,3) > 1) CALL finish(routine, 'Dimension mismatch, field3d with nblks > 1; calc with REAL pool%element(:)')
      structOut%carbon          = struct1%carbon      * field3d(:,:,1)
      structOut%nitrogen        = struct1%nitrogen    * field3d(:,:,1)
      structOut%phosphorus      = struct1%phosphorus  * field3d(:,:,1)
      structOut%carbon13        = struct1%carbon13    * field3d(:,:,1)
      structOut%carbon14        = struct1%carbon14    * field3d(:,:,1)
      structOut%nitrogen15      = struct1%nitrogen15  * field3d(:,:,1)
    ELSE
      structOut%domain_c        = struct1%domain_c    * field3d(:,:,:)
      structOut%domain_n        = struct1%domain_n    * field3d(:,:,:)
      structOut%domain_p        = struct1%domain_p    * field3d(:,:,:)
      structOut%domain_c13      = struct1%domain_c13  * field3d(:,:,:)
      structOut%domain_c14      = struct1%domain_c14  * field3d(:,:,:)
      structOut%domain_n15      = struct1%domain_n15  * field3d(:,:,:)
    END IF

  END FUNCTION Multiply_sb_elem_p_with_field3d
  !--------------------------------------------------------------------------------------------------------------------------------
  ! field3d(:,:,:) * t_sb_elements_pool        | ON_CHUNK & ON_DOMAIN | field3d(:,:,1) & field3d(:,:,:)
  !--------------------------------------------------------------------------------------------------------------------------------
  FUNCTION Multiply_field3d_with_sb_elem_p(field3d, struct1) RESULT(structOut)
  
    REAL(wp),                      INTENT(in)  :: field3d(:,:,:)
    CLASS(t_sb_elements_pool),     INTENT(in)  :: struct1
    TYPE(t_sb_elements_pool)                   :: structOut

    TYPE(t_subset)              :: sub  
    INTEGER                     :: nvgrid
    CHARACTER(len=*), PARAMETER :: routine = modname//':Multiply_field3d_with_sb_elem_p'

    sub = struct1%Get_subset()
    IF (sub%type /= ON_CHUNK) THEN  ! get domain dimension not available via Get_subset() during process init
      sub%nc = SIZE(struct1%domain_c(:,:,:), 1)
      nvgrid = SIZE(struct1%domain_c(:,:,:), 2)
      sub%nb = SIZE(struct1%domain_c(:,:,:), 3)
    ELSE
      nvgrid = SIZE(struct1%carbon(:,:), 2)
    ENDIF
    CALL structOut%Init_op_ov_hlp('op_ov_hlp', 'op_ov_hlp', struct1%owner_model_id, &
      &                           sub%type, sub%nc, nvgrid, sub%nb, .TRUE., TRIM(struct1%element_unit))

    ! run calculation
    IF (sub%type == ON_CHUNK) THEN
      IF (SIZE(field3d,3) > 1) CALL finish(routine, 'Dimension mismatch, field3d with nblks > 1; calc with REAL pool%element(:)')
      structOut%carbon          = field3d(:,:,1) * struct1%carbon    
      structOut%nitrogen        = field3d(:,:,1) * struct1%nitrogen  
      structOut%phosphorus      = field3d(:,:,1) * struct1%phosphorus
      structOut%carbon13        = field3d(:,:,1) * struct1%carbon13  
      structOut%carbon14        = field3d(:,:,1) * struct1%carbon14  
      structOut%nitrogen15      = field3d(:,:,1) * struct1%nitrogen15
    ELSE
      structOut%domain_c        = field3d(:,:,:) * struct1%domain_c  
      structOut%domain_n        = field3d(:,:,:) * struct1%domain_n  
      structOut%domain_p        = field3d(:,:,:) * struct1%domain_p  
      structOut%domain_c13      = field3d(:,:,:) * struct1%domain_c13
      structOut%domain_c14      = field3d(:,:,:) * struct1%domain_c14
      structOut%domain_n15      = field3d(:,:,:) * struct1%domain_n15
    END IF

  END FUNCTION Multiply_field3d_with_sb_elem_p
  !--------------------------------------------------------------------------------------------------------------------------------
  ! t_sb_elements_pool * t_sb_elements_pool | ON_CHUNK & ON_DOMAIN
  !--------------------------------------------------------------------------------------------------------------------------------
  FUNCTION Multiply_sb_elem_p_with_sb_elem_p(struct1, struct2) RESULT(structOut)
  
    CLASS(t_sb_elements_pool),     INTENT(in)  :: struct1, struct2
    TYPE(t_sb_elements_pool)                   :: structOut

    TYPE(t_subset)              :: sub  
    INTEGER                     :: nvgrid
    CHARACTER(len=*), PARAMETER :: routine = modname//':Multiply_sb_elem_p_with_sb_elem_p'

    sub = struct1%Get_subset()
    IF (sub%type /= ON_CHUNK) THEN  ! get domain dimension not available via Get_subset() during process init
      sub%nc = SIZE(struct1%domain_c(:,:,:), 1)
      nvgrid = SIZE(struct1%domain_c(:,:,:), 2)
      sub%nb = SIZE(struct1%domain_c(:,:,:), 3)
    ELSE
      nvgrid = SIZE(struct1%carbon(:,:), 2)
    ENDIF
    CALL structOut%Init_op_ov_hlp('op_ov_hlp', 'op_ov_hlp', struct1%owner_model_id, &
      &                           sub%type, sub%nc, nvgrid, sub%nb, .TRUE., TRIM(struct1%element_unit))

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

  END FUNCTION Multiply_sb_elem_p_with_sb_elem_p

  !--------------------------------------------------------------------------------------------------------------------------------
  ! / divide
  !--------------------------------------------------------------------------------------------------------------------------------
  !--------------------------------------------------------------------------------------------------------------------------------
  ! t_sb_elements_pool / scalar        | ON_CHUNK & ON_DOMAIN
  !--------------------------------------------------------------------------------------------------------------------------------
  FUNCTION Div_sb_elem_p_by_scalar(struct1, scalar) RESULT(structOut)
  
    CLASS(t_sb_elements_pool),     INTENT(in)  :: struct1
    REAL(wp),                      INTENT(in)  :: scalar
    TYPE(t_sb_elements_pool)                   :: structOut

    TYPE(t_subset)              :: sub  
    INTEGER                     :: nvgrid
    CHARACTER(len=*), PARAMETER :: routine = modname//':Div_sb_elem_p_by_scalar'

    sub = struct1%Get_subset()
    IF (sub%type /= ON_CHUNK) THEN  ! get domain dimension not available via Get_subset() during process init
      sub%nc = SIZE(struct1%domain_c(:,:,:), 1)
      nvgrid = SIZE(struct1%domain_c(:,:,:), 2)
      sub%nb = SIZE(struct1%domain_c(:,:,:), 3)
    ELSE
      nvgrid = SIZE(struct1%carbon(:,:), 2)
    ENDIF
    CALL structOut%Init_op_ov_hlp('op_ov_hlp', 'op_ov_hlp', struct1%owner_model_id, &
      &                           sub%type, sub%nc, nvgrid, sub%nb, .TRUE., TRIM(struct1%element_unit))

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

  END FUNCTION Div_sb_elem_p_by_scalar
  !--------------------------------------------------------------------------------------------------------------------------------
  ! t_sb_elements_pool / field2d(:,:)        | ON_CHUNK
  !--------------------------------------------------------------------------------------------------------------------------------
  FUNCTION Div_sb_elem_p_by_field2d(struct1, field2d) RESULT(structOut)
  
    CLASS(t_sb_elements_pool),     INTENT(in)  :: struct1
    REAL(wp),                      INTENT(in)  :: field2d(:,:)
    TYPE(t_sb_elements_pool)                   :: structOut

    TYPE(t_subset)              :: sub  
    INTEGER                     :: nvgrid
    CHARACTER(len=*), PARAMETER :: routine = modname//':Div_sb_elem_p_by_field2d'

    sub = struct1%Get_subset()
    IF (sub%type /= ON_CHUNK) THEN  ! get domain dimension not available via Get_subset() during process init
      sub%nc = SIZE(struct1%domain_c(:,:,:), 1)
      nvgrid = SIZE(struct1%domain_c(:,:,:), 2)
      sub%nb = SIZE(struct1%domain_c(:,:,:), 3)
    ELSE
      nvgrid = SIZE(struct1%carbon(:,:), 2)
    ENDIF
    CALL structOut%Init_op_ov_hlp('op_ov_hlp', 'op_ov_hlp', struct1%owner_model_id, &
      &                           sub%type, sub%nc, nvgrid, sub%nb, .TRUE., TRIM(struct1%element_unit))

    ! run calculation
    IF (sub%type == ON_CHUNK) THEN
      structOut%carbon          = struct1%carbon      / field2d(:,:)
      structOut%nitrogen        = struct1%nitrogen    / field2d(:,:)
      structOut%phosphorus      = struct1%phosphorus  / field2d(:,:)
      structOut%carbon13        = struct1%carbon13    / field2d(:,:)
      structOut%carbon14        = struct1%carbon14    / field2d(:,:)
      structOut%nitrogen15      = struct1%nitrogen15  / field2d(:,:)
    ELSE
      CALL finish(routine, 'Dimension mismatch, cannot calc with 2d-field and jsb_var real3d on domain')
    END IF

  END FUNCTION Div_sb_elem_p_by_field2d
  !--------------------------------------------------------------------------------------------------------------------------------
  ! t_sb_elements_pool / field3d(:,:,:)        | ON_CHUNK & ON_DOMAIN | field3d(:,:,1) & field3d(:,:,:)
  !--------------------------------------------------------------------------------------------------------------------------------
  FUNCTION Div_sb_elem_p_by_field3d(struct1, field3d) RESULT(structOut)
  
    CLASS(t_sb_elements_pool),     INTENT(in)  :: struct1
    REAL(wp),                      INTENT(in)  :: field3d(:,:,:)
    TYPE(t_sb_elements_pool)                   :: structOut

    TYPE(t_subset)              :: sub  
    INTEGER                     :: nvgrid
    CHARACTER(len=*), PARAMETER :: routine = modname//':Div_sb_elem_p_by_field3d'

    sub = struct1%Get_subset()
    IF (sub%type /= ON_CHUNK) THEN  ! get domain dimension not available via Get_subset() during process init
      sub%nc = SIZE(struct1%domain_c(:,:,:), 1)
      nvgrid = SIZE(struct1%domain_c(:,:,:), 2)
      sub%nb = SIZE(struct1%domain_c(:,:,:), 3)
    ELSE
      nvgrid = SIZE(struct1%carbon(:,:), 2)
    ENDIF
    CALL structOut%Init_op_ov_hlp('op_ov_hlp', 'op_ov_hlp', struct1%owner_model_id, &
      &                           sub%type, sub%nc, nvgrid, sub%nb, .TRUE., TRIM(struct1%element_unit))

    ! run calculation
    IF (sub%type == ON_CHUNK) THEN
      IF (SIZE(field3d,3) > 1) CALL finish(routine, 'Dimension mismatch, field3d with nblks > 1; calc with REAL pool%element(:)')
      structOut%carbon          = struct1%carbon      / field3d(:,:,1)
      structOut%nitrogen        = struct1%nitrogen    / field3d(:,:,1)
      structOut%phosphorus      = struct1%phosphorus  / field3d(:,:,1)
      structOut%carbon13        = struct1%carbon13    / field3d(:,:,1)
      structOut%carbon14        = struct1%carbon14    / field3d(:,:,1)
      structOut%nitrogen15      = struct1%nitrogen15  / field3d(:,:,1)
    ELSE
      structOut%domain_c        = struct1%domain_c    / field3d(:,:,:)
      structOut%domain_n        = struct1%domain_n    / field3d(:,:,:)
      structOut%domain_p        = struct1%domain_p    / field3d(:,:,:)
      structOut%domain_c13      = struct1%domain_c13  / field3d(:,:,:)
      structOut%domain_c14      = struct1%domain_c14  / field3d(:,:,:)
      structOut%domain_n15      = struct1%domain_n15  / field3d(:,:,:)
    END IF

  END FUNCTION Div_sb_elem_p_by_field3d
  !--------------------------------------------------------------------------------------------------------------------------------
  ! t_sb_elements_pool / t_sb_elements_pool | ON_CHUNK & ON_DOMAIN
  !--------------------------------------------------------------------------------------------------------------------------------
  FUNCTION Div_sb_elem_p_by_sb_elem_p(struct1, struct2) RESULT(structOut)
  
    CLASS(t_sb_elements_pool),     INTENT(in)  :: struct1, struct2
    TYPE(t_sb_elements_pool)                   :: structOut

    TYPE(t_subset)              :: sub  
    INTEGER                     :: nvgrid
    CHARACTER(len=*), PARAMETER :: routine = modname//':Add_sb_elem_p_and_sb_elem_p'

    sub = struct1%Get_subset()
    IF (sub%type /= ON_CHUNK) THEN  ! get domain dimension not available via Get_subset() during process init
      sub%nc = SIZE(struct1%domain_c(:,:,:), 1)
      nvgrid = SIZE(struct1%domain_c(:,:,:), 2)
      sub%nb = SIZE(struct1%domain_c(:,:,:), 3)
    ELSE
      nvgrid = SIZE(struct1%carbon(:,:), 2)
    ENDIF
    CALL structOut%Init_op_ov_hlp('op_ov_hlp', 'op_ov_hlp', struct1%owner_model_id, &
      &                           sub%type, sub%nc, nvgrid, sub%nb, .TRUE., TRIM(struct1%element_unit))

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

  END FUNCTION Div_sb_elem_p_by_sb_elem_p

  !--------------------------------------------------------------------------------------------------------------------------------
  !> operator overload and assignment routines for tissue pools
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

    IF (ALLOCATED(this%op_ov_hlp_2d_c))     DEALLOCATE(this%op_ov_hlp_2d_c)
    IF (ALLOCATED(this%op_ov_hlp_2d_n))     DEALLOCATE(this%op_ov_hlp_2d_n)
    IF (ALLOCATED(this%op_ov_hlp_2d_p))     DEALLOCATE(this%op_ov_hlp_2d_p)
    IF (ALLOCATED(this%op_ov_hlp_2d_c13))   DEALLOCATE(this%op_ov_hlp_2d_c13)
    IF (ALLOCATED(this%op_ov_hlp_2d_c14))   DEALLOCATE(this%op_ov_hlp_2d_c14)
    IF (ALLOCATED(this%op_ov_hlp_2d_n15))   DEALLOCATE(this%op_ov_hlp_2d_n15)
    IF (ALLOCATED(this%op_ov_hlp_3d_c))     DEALLOCATE(this%op_ov_hlp_3d_c)
    IF (ALLOCATED(this%op_ov_hlp_3d_n))     DEALLOCATE(this%op_ov_hlp_3d_n)
    IF (ALLOCATED(this%op_ov_hlp_3d_p))     DEALLOCATE(this%op_ov_hlp_3d_p)
    IF (ALLOCATED(this%op_ov_hlp_3d_c13))   DEALLOCATE(this%op_ov_hlp_3d_c13)
    IF (ALLOCATED(this%op_ov_hlp_3d_c14))   DEALLOCATE(this%op_ov_hlp_3d_c14)
    IF (ALLOCATED(this%op_ov_hlp_3d_n15))   DEALLOCATE(this%op_ov_hlp_3d_n15)
    
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
    
    IF (ALLOCATED(this%shortname))          DEALLOCATE(this%shortname)
    IF (ALLOCATED(this%name))               DEALLOCATE(this%name)
    ! not allocated: this%path
    IF (ALLOCATED(this%element_unit))       DEALLOCATE(this%element_unit)
    IF (ASSOCIATED(this%pool_list))         DEALLOCATE(this%pool_list)
    IF (ASSOCIATED(this%leaf_list))         DEALLOCATE(this%leaf_list)
    IF (ALLOCATED(this%element_list))       DEALLOCATE(this%element_list)
    IF (ALLOCATED(this%var_list))           DEALLOCATE(this%var_list)    
    ! this%temporary of t_hermetic is PRIVATE
   
  END SUBROUTINE finalize_sb_components_pool

  !--------------------------------------------------------------------------------------------------------------------------------
  ! = assign
  !--------------------------------------------------------------------------------------------------------------------------------
  !--------------------------------------------------------------------------------------------------------------------------------
  ! sb_components_pool = scalar | ON_CHUNK & ON_DOMAIN
  !--------------------------------------------------------------------------------------------------------------------------------
  SUBROUTINE Assign_scalar_to_sb_comp_p(lhs, rhs)

    CLASS(t_sb_components_pool),  INTENT(inout) :: lhs
    REAL(wp),                     INTENT(in)    :: rhs

    TYPE(t_subset)              :: sub
    CHARACTER(len=*), PARAMETER :: routine = modname//':Assign_scalar_to_sb_comp_p'
    
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
    
    lhs%dom              = rhs
    lhs%dom_assoc        = rhs
    lhs%soluable_litter  = rhs
    lhs%polymeric_litter = rhs
    lhs%woody_litter     = rhs
    lhs%fungi            = rhs
    lhs%mycorrhiza       = rhs
    lhs%microbial        = rhs
    lhs%residue          = rhs
    lhs%residue_assoc    = rhs
    
  END SUBROUTINE Assign_scalar_to_sb_comp_p
  !--------------------------------------------------------------------------------------------------------------------------------
  ! sb_components_pool = field2d(:,:) | ON_CHUNK
  !--------------------------------------------------------------------------------------------------------------------------------
  SUBROUTINE Assign_field2d_to_sb_comp_p(lhs, rhs)
  
    CLASS(t_sb_components_pool),  INTENT(inout) :: lhs
    REAL(wp),                     INTENT(in)    :: rhs(:,:)

    TYPE(t_subset)              :: sub
    CHARACTER(len=*), PARAMETER :: routine = modname//':Assign_field2d_to_sb_comp_p'
    
    sub = lhs%Get_subset()

    IF (sub%type == ON_CHUNK) THEN
      lhs%carbon     = rhs(:,:)
      lhs%nitrogen   = rhs(:,:)
      lhs%phosphorus = rhs(:,:)
      lhs%carbon13   = rhs(:,:)
      lhs%carbon14   = rhs(:,:)
      lhs%nitrogen15 = rhs(:,:)
    ELSE
      CALL finish(routine, lhs%shortname//' Dimension mismatch, cannot assign 2d-field to jsb_var real3d on domain')
    END IF

    lhs%dom              = rhs(:,:)
    lhs%dom_assoc        = rhs(:,:)
    lhs%soluable_litter  = rhs(:,:)
    lhs%polymeric_litter = rhs(:,:)
    lhs%woody_litter     = rhs(:,:)
    lhs%fungi            = rhs(:,:)
    lhs%mycorrhiza       = rhs(:,:)
    lhs%microbial        = rhs(:,:)
    lhs%residue          = rhs(:,:)
    lhs%residue_assoc    = rhs(:,:)
    
  END SUBROUTINE Assign_field2d_to_sb_comp_p
  !--------------------------------------------------------------------------------------------------------------------------------
  ! sb_components_pool = field3d(:,:,:) | ON_CHUNK & ON_DOMAIN | field3d(:,:,1) & field3d(:,:,:)
  !--------------------------------------------------------------------------------------------------------------------------------
  SUBROUTINE Assign_field3d_to_sb_comp_p(lhs, rhs)
  
    CLASS(t_sb_components_pool),  INTENT(inout) :: lhs
    REAL(wp),                     INTENT(in)    :: rhs(:,:,:)

    TYPE(t_subset)              :: sub
    CHARACTER(len=*), PARAMETER :: routine = modname//':Assign_field3d_to_sb_comp_p'

    sub = lhs%Get_subset()

    IF (sub%type == ON_CHUNK) THEN
      IF (SIZE(rhs,3) > 1) CALL finish(routine, 'Dimension mismatch, field3d with nblks > 1; assignment to REAL pool%element(:)')
      lhs%carbon     = rhs(:,:,1)
      lhs%nitrogen   = rhs(:,:,1)
      lhs%phosphorus = rhs(:,:,1)
      lhs%carbon13   = rhs(:,:,1)
      lhs%carbon14   = rhs(:,:,1)
      lhs%nitrogen15 = rhs(:,:,1)
    ELSE
      lhs%domain_c   = rhs(:,:,:)
      lhs%domain_n   = rhs(:,:,:)
      lhs%domain_p   = rhs(:,:,:)
      lhs%domain_c13 = rhs(:,:,:)
      lhs%domain_c14 = rhs(:,:,:)
      lhs%domain_n15 = rhs(:,:,:)
    END IF

    lhs%dom              = rhs(:,:,:)
    lhs%dom_assoc        = rhs(:,:,:)
    lhs%soluable_litter  = rhs(:,:,:)
    lhs%polymeric_litter = rhs(:,:,:)
    lhs%woody_litter     = rhs(:,:,:)
    lhs%fungi            = rhs(:,:,:)
    lhs%mycorrhiza       = rhs(:,:,:)
    lhs%microbial        = rhs(:,:,:)
    lhs%residue          = rhs(:,:,:)
    lhs%residue_assoc    = rhs(:,:,:)

  END SUBROUTINE Assign_field3d_to_sb_comp_p
  !--------------------------------------------------------------------------------------------------------------------------------
  ! sb_components_pool = sb_components_pool | ON_CHUNK & ON_DOMAIN
  !--------------------------------------------------------------------------------------------------------------------------------
  SUBROUTINE Assign_sb_comp_p_to_sb_comp_p(lhs, rhs)

    CLASS(t_sb_components_pool), INTENT(inout) :: lhs
    CLASS(t_sb_components_pool), INTENT(in)    :: rhs

    TYPE(t_subset)              :: sub
    CHARACTER(len=*), PARAMETER :: routine = modname//':Assign_sb_comp_p_to_sb_comp_p'
    
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
    
    lhs%dom              = rhs%dom             
    lhs%dom_assoc        = rhs%dom_assoc       
    lhs%soluable_litter  = rhs%soluable_litter 
    lhs%polymeric_litter = rhs%polymeric_litter
    lhs%woody_litter     = rhs%woody_litter    
    lhs%fungi            = rhs%fungi           
    lhs%mycorrhiza       = rhs%mycorrhiza      
    lhs%microbial        = rhs%microbial       
    lhs%residue          = rhs%residue         
    lhs%residue_assoc    = rhs%residue_assoc   

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
    TYPE(t_sb_components_pool)                    :: structOut

    TYPE(t_subset)              :: sub  
    INTEGER                     :: nvgrid
    CHARACTER(len=*), PARAMETER :: routine = modname//':Add_sb_comp_p_and_scalar'

    sub = struct1%Get_subset()
    IF (sub%type /= ON_CHUNK) THEN  ! get domain dimension not available via Get_subset() during process init
      sub%nc = SIZE(struct1%domain_c(:,:,:), 1)
      nvgrid = SIZE(struct1%domain_c(:,:,:), 2)
      sub%nb = SIZE(struct1%domain_c(:,:,:), 3)
    ELSE
      nvgrid = SIZE(struct1%carbon(:,:), 2)
    ENDIF
    CALL structOut%Init_op_ov_hlp('op_ov_hlp', 'op_ov_hlp', struct1%owner_model_id, &
      &                           sub%type, sub%nc, nvgrid, sub%nb, .TRUE., TRIM(struct1%element_unit), 'veg_p')

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

    structOut%dom              = struct1%dom              + scalar
    structOut%dom_assoc        = struct1%dom_assoc        + scalar
    structOut%soluable_litter  = struct1%soluable_litter  + scalar
    structOut%polymeric_litter = struct1%polymeric_litter + scalar
    structOut%woody_litter     = struct1%woody_litter     + scalar
    structOut%fungi            = struct1%fungi            + scalar
    structOut%mycorrhiza       = struct1%mycorrhiza       + scalar
    structOut%microbial        = struct1%microbial        + scalar
    structOut%residue          = struct1%residue          + scalar
    structOut%residue_assoc    = struct1%residue_assoc    + scalar

  END FUNCTION Add_sb_comp_p_and_scalar
  !--------------------------------------------------------------------------------------------------------------------------------
  ! scalar + sb_components_pool        | ON_CHUNK & ON_DOMAIN
  !--------------------------------------------------------------------------------------------------------------------------------
  FUNCTION Add_scalar_and_sb_comp_p(scalar, struct1) RESULT(structOut)
  
    REAL(wp),                         INTENT(in)  :: scalar
    CLASS(t_sb_components_pool),      INTENT(in)  :: struct1
    TYPE(t_sb_components_pool)                    :: structOut
  
    TYPE(t_subset)              :: sub  
    INTEGER                     :: nvgrid
    CHARACTER(len=*), PARAMETER :: routine = modname//':Add_scalar_and_sb_comp_p'
  
    sub = struct1%Get_subset()
    IF (sub%type /= ON_CHUNK) THEN  ! get domain dimension not available via Get_subset() during process init
      sub%nc = SIZE(struct1%domain_c(:,:,:), 1)
      nvgrid = SIZE(struct1%domain_c(:,:,:), 2)
      sub%nb = SIZE(struct1%domain_c(:,:,:), 3)
    ELSE
      nvgrid = SIZE(struct1%carbon(:,:), 2)
    ENDIF
    CALL structOut%Init_op_ov_hlp('op_ov_hlp', 'op_ov_hlp', struct1%owner_model_id, &
      &                           sub%type, sub%nc, nvgrid, sub%nb, .TRUE., TRIM(struct1%element_unit), 'veg_p')

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

    structOut%dom              = scalar + struct1%dom             
    structOut%dom_assoc        = scalar + struct1%dom_assoc       
    structOut%soluable_litter  = scalar + struct1%soluable_litter 
    structOut%polymeric_litter = scalar + struct1%polymeric_litter
    structOut%woody_litter     = scalar + struct1%woody_litter    
    structOut%fungi            = scalar + struct1%fungi           
    structOut%mycorrhiza       = scalar + struct1%mycorrhiza      
    structOut%microbial        = scalar + struct1%microbial       
    structOut%residue          = scalar + struct1%residue         
    structOut%residue_assoc    = scalar + struct1%residue_assoc   
  
  END FUNCTION Add_scalar_and_sb_comp_p
  !--------------------------------------------------------------------------------------------------------------------------------
  ! sb_components_pool + field2d(:,:)        | ON_CHUNK
  !--------------------------------------------------------------------------------------------------------------------------------
  FUNCTION Add_sb_comp_p_and_field2d(struct1, field2d) RESULT(structOut)
  
    CLASS(t_sb_components_pool),      INTENT(in)  :: struct1
    REAL(wp),                         INTENT(in)  :: field2d(:,:)
    TYPE(t_sb_components_pool)                    :: structOut
  
    TYPE(t_subset)              :: sub  
    INTEGER                     :: nvgrid
    CHARACTER(len=*), PARAMETER :: routine = modname//':Add_sb_comp_p_and_field2d'
  
    sub = struct1%Get_subset()
    IF (sub%type /= ON_CHUNK) THEN  ! get domain dimension not available via Get_subset() during process init
      sub%nc = SIZE(struct1%domain_c(:,:,:), 1)
      nvgrid = SIZE(struct1%domain_c(:,:,:), 2)
      sub%nb = SIZE(struct1%domain_c(:,:,:), 3)
    ELSE
      nvgrid = SIZE(struct1%carbon(:,:), 2)
    ENDIF
    CALL structOut%Init_op_ov_hlp('op_ov_hlp', 'op_ov_hlp', struct1%owner_model_id, &
      &                           sub%type, sub%nc, nvgrid, sub%nb, .TRUE., TRIM(struct1%element_unit), 'veg_p')
  
    ! run calculation
    IF (sub%type == ON_CHUNK) THEN
      structOut%carbon          = struct1%carbon      + field2d(:,:)
      structOut%nitrogen        = struct1%nitrogen    + field2d(:,:)
      structOut%phosphorus      = struct1%phosphorus  + field2d(:,:)
      structOut%carbon13        = struct1%carbon13    + field2d(:,:)
      structOut%carbon14        = struct1%carbon14    + field2d(:,:)
      structOut%nitrogen15      = struct1%nitrogen15  + field2d(:,:)
    ELSE
      CALL finish(routine, 'Dimension mismatch, cannot calc with 2d-field and jsb_var real3d on domain')
    END IF

    structOut%dom              = struct1%dom              + field2d(:,:)
    structOut%dom_assoc        = struct1%dom_assoc        + field2d(:,:)
    structOut%soluable_litter  = struct1%soluable_litter  + field2d(:,:)
    structOut%polymeric_litter = struct1%polymeric_litter + field2d(:,:)
    structOut%woody_litter     = struct1%woody_litter     + field2d(:,:)
    structOut%fungi            = struct1%fungi            + field2d(:,:)
    structOut%mycorrhiza       = struct1%mycorrhiza       + field2d(:,:)
    structOut%microbial        = struct1%microbial        + field2d(:,:)
    structOut%residue          = struct1%residue          + field2d(:,:)
    structOut%residue_assoc    = struct1%residue_assoc    + field2d(:,:)
  
  END FUNCTION Add_sb_comp_p_and_field2d
  !--------------------------------------------------------------------------------------------------------------------------------
  ! field2d(:,:) + sb_components_pool        | ON_CHUNK
  !--------------------------------------------------------------------------------------------------------------------------------
  FUNCTION Add_field2d_and_sb_comp_p(field2d, struct1) RESULT(structOut)
  
    REAL(wp),                         INTENT(in)  :: field2d(:,:)
    CLASS(t_sb_components_pool),      INTENT(in)  :: struct1
    TYPE(t_sb_components_pool)                    :: structOut
  
    TYPE(t_subset)              :: sub  
    INTEGER                     :: nvgrid
    CHARACTER(len=*), PARAMETER :: routine = modname//':Add_field2d_and_sb_comp_p'
  
    sub = struct1%Get_subset()
    IF (sub%type /= ON_CHUNK) THEN  ! get domain dimension not available via Get_subset() during process init
      sub%nc = SIZE(struct1%domain_c(:,:,:), 1)
      nvgrid = SIZE(struct1%domain_c(:,:,:), 2)
      sub%nb = SIZE(struct1%domain_c(:,:,:), 3)
    ELSE
      nvgrid = SIZE(struct1%carbon(:,:), 2)
    ENDIF
    CALL structOut%Init_op_ov_hlp('op_ov_hlp', 'op_ov_hlp', struct1%owner_model_id, &
      &                           sub%type, sub%nc, nvgrid, sub%nb, .TRUE., TRIM(struct1%element_unit), 'veg_p')
  
    ! run calculation
    IF (sub%type == ON_CHUNK) THEN
      structOut%carbon          = field2d(:,:) + struct1%carbon    
      structOut%nitrogen        = field2d(:,:) + struct1%nitrogen  
      structOut%phosphorus      = field2d(:,:) + struct1%phosphorus
      structOut%carbon13        = field2d(:,:) + struct1%carbon13  
      structOut%carbon14        = field2d(:,:) + struct1%carbon14  
      structOut%nitrogen15      = field2d(:,:) + struct1%nitrogen15
    ELSE
      CALL finish(routine, 'Dimension mismatch, cannot calc with 2d-field and jsb_var real3d on domain')
    END IF

    structOut%dom              = field2d(:,:) + struct1%dom             
    structOut%dom_assoc        = field2d(:,:) + struct1%dom_assoc       
    structOut%soluable_litter  = field2d(:,:) + struct1%soluable_litter 
    structOut%polymeric_litter = field2d(:,:) + struct1%polymeric_litter
    structOut%woody_litter     = field2d(:,:) + struct1%woody_litter    
    structOut%fungi            = field2d(:,:) + struct1%fungi           
    structOut%mycorrhiza       = field2d(:,:) + struct1%mycorrhiza      
    structOut%microbial        = field2d(:,:) + struct1%microbial       
    structOut%residue          = field2d(:,:) + struct1%residue         
    structOut%residue_assoc    = field2d(:,:) + struct1%residue_assoc   
  
  END FUNCTION Add_field2d_and_sb_comp_p
  !--------------------------------------------------------------------------------------------------------------------------------
  ! sb_components_pool + field3d(:,:,:)        | ON_CHUNK & ON_DOMAIN | field3d(:,:,1) & field3d(:,:,:)
  !--------------------------------------------------------------------------------------------------------------------------------
  FUNCTION Add_sb_comp_p_and_field3d(struct1, field3d) RESULT(structOut)
  
    CLASS(t_sb_components_pool),      INTENT(in)  :: struct1
    REAL(wp),                         INTENT(in)  :: field3d(:,:,:)
    TYPE(t_sb_components_pool)                    :: structOut
  
    TYPE(t_subset)              :: sub  
    INTEGER                     :: nvgrid
    CHARACTER(len=*), PARAMETER :: routine = modname//':Add_sb_comp_p_and_field3d'
  
    sub = struct1%Get_subset()
    IF (sub%type /= ON_CHUNK) THEN  ! get domain dimension not available via Get_subset() during process init
      sub%nc = SIZE(struct1%domain_c(:,:,:), 1)
      nvgrid = SIZE(struct1%domain_c(:,:,:), 2)
      sub%nb = SIZE(struct1%domain_c(:,:,:), 3)
    ELSE
      nvgrid = SIZE(struct1%carbon(:,:), 2)
    ENDIF
    CALL structOut%Init_op_ov_hlp('op_ov_hlp', 'op_ov_hlp', struct1%owner_model_id, &
      &                           sub%type, sub%nc, nvgrid, sub%nb, .TRUE., TRIM(struct1%element_unit), 'veg_p')
  
    ! run calculation
    IF (sub%type == ON_CHUNK) THEN
      IF (SIZE(field3d,3) > 1) CALL finish(routine, 'Dimension mismatch, field3d with nblks > 1; calc with REAL pool%element(:)')
      structOut%carbon          = struct1%carbon      + field3d(:,:,1)
      structOut%nitrogen        = struct1%nitrogen    + field3d(:,:,1)
      structOut%phosphorus      = struct1%phosphorus  + field3d(:,:,1)
      structOut%carbon13        = struct1%carbon13    + field3d(:,:,1)
      structOut%carbon14        = struct1%carbon14    + field3d(:,:,1)
      structOut%nitrogen15      = struct1%nitrogen15  + field3d(:,:,1)
    ELSE
      structOut%domain_c        = struct1%domain_c    + field3d(:,:,:)
      structOut%domain_n        = struct1%domain_n    + field3d(:,:,:)
      structOut%domain_p        = struct1%domain_p    + field3d(:,:,:)
      structOut%domain_c13      = struct1%domain_c13  + field3d(:,:,:)
      structOut%domain_c14      = struct1%domain_c14  + field3d(:,:,:)
      structOut%domain_n15      = struct1%domain_n15  + field3d(:,:,:)
    END IF

    structOut%dom              = struct1%dom              + field3d(:,:,:)
    structOut%dom_assoc        = struct1%dom_assoc        + field3d(:,:,:)
    structOut%soluable_litter  = struct1%soluable_litter  + field3d(:,:,:)
    structOut%polymeric_litter = struct1%polymeric_litter + field3d(:,:,:)
    structOut%woody_litter     = struct1%woody_litter     + field3d(:,:,:)
    structOut%fungi            = struct1%fungi            + field3d(:,:,:)
    structOut%mycorrhiza       = struct1%mycorrhiza       + field3d(:,:,:)
    structOut%microbial        = struct1%microbial        + field3d(:,:,:)
    structOut%residue          = struct1%residue          + field3d(:,:,:)
    structOut%residue_assoc    = struct1%residue_assoc    + field3d(:,:,:)
  
  END FUNCTION Add_sb_comp_p_and_field3d
  !--------------------------------------------------------------------------------------------------------------------------------
  ! sb_components_pool + sb_components_pool        | ON_CHUNK & ON_DOMAIN
  !--------------------------------------------------------------------------------------------------------------------------------
  FUNCTION Add_sb_comp_p_and_sb_comp_p(struct1, struct2) RESULT(structOut)
  
    CLASS(t_sb_components_pool),      INTENT(in)  :: struct1, struct2
    TYPE(t_sb_components_pool)                    :: structOut
  
    TYPE(t_subset)              :: sub  
    INTEGER                     :: nvgrid
    CHARACTER(len=*), PARAMETER :: routine = modname//':Add_sb_comp_p_and_sb_comp_p'
  
    sub = struct1%Get_subset()
    IF (sub%type /= ON_CHUNK) THEN  ! get domain dimension not available via Get_subset() during process init
      sub%nc = SIZE(struct1%domain_c(:,:,:), 1)
      nvgrid = SIZE(struct1%domain_c(:,:,:), 2)
      sub%nb = SIZE(struct1%domain_c(:,:,:), 3)
    ELSE
      nvgrid = SIZE(struct1%carbon(:,:), 2)
    ENDIF
    CALL structOut%Init_op_ov_hlp('op_ov_hlp', 'op_ov_hlp', struct1%owner_model_id, &
      &                           sub%type, sub%nc, nvgrid, sub%nb, .TRUE., TRIM(struct1%element_unit), 'veg_p')
  
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

    structOut%dom              = struct1%dom              + struct2%dom             
    structOut%dom_assoc        = struct1%dom_assoc        + struct2%dom_assoc       
    structOut%soluable_litter  = struct1%soluable_litter  + struct2%soluable_litter 
    structOut%polymeric_litter = struct1%polymeric_litter + struct2%polymeric_litter
    structOut%woody_litter     = struct1%woody_litter     + struct2%woody_litter    
    structOut%fungi            = struct1%fungi            + struct2%fungi           
    structOut%mycorrhiza       = struct1%mycorrhiza       + struct2%mycorrhiza      
    structOut%microbial        = struct1%microbial        + struct2%microbial       
    structOut%residue          = struct1%residue          + struct2%residue         
    structOut%residue_assoc    = struct1%residue_assoc    + struct2%residue_assoc   
  
  END FUNCTION Add_sb_comp_p_and_sb_comp_p
  
  !--------------------------------------------------------------------------------------------------------------------------------
  ! - minus
  !--------------------------------------------------------------------------------------------------------------------------------
  !--------------------------------------------------------------------------------------------------------------------------------
  ! sb_components_pool - scalar        | ON_CHUNK & ON_DOMAIN
  !--------------------------------------------------------------------------------------------------------------------------------
  FUNCTION Subtract_sb_comp_p_minus_scalar(struct1, scalar) RESULT(structOut)
  
    CLASS(t_sb_components_pool),      INTENT(in)  :: struct1
    REAL(wp),                         INTENT(in)  :: scalar
    TYPE(t_sb_components_pool)                    :: structOut

    TYPE(t_subset)              :: sub  
    INTEGER                     :: nvgrid
    CHARACTER(len=*), PARAMETER :: routine = modname//':Subtract_sb_comp_p_minus_scalar'

    sub = struct1%Get_subset()
    IF (sub%type /= ON_CHUNK) THEN  ! get domain dimension not available via Get_subset() during process init
      sub%nc = SIZE(struct1%domain_c(:,:,:), 1)
      nvgrid = SIZE(struct1%domain_c(:,:,:), 2)
      sub%nb = SIZE(struct1%domain_c(:,:,:), 3)
    ELSE
      nvgrid = SIZE(struct1%carbon(:,:), 2)
    ENDIF
    CALL structOut%Init_op_ov_hlp('op_ov_hlp', 'op_ov_hlp', struct1%owner_model_id, &
      &                           sub%type, sub%nc, nvgrid, sub%nb, .TRUE., TRIM(struct1%element_unit), 'veg_p')

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

    structOut%dom              = struct1%dom              - scalar
    structOut%dom_assoc        = struct1%dom_assoc        - scalar
    structOut%soluable_litter  = struct1%soluable_litter  - scalar
    structOut%polymeric_litter = struct1%polymeric_litter - scalar
    structOut%woody_litter     = struct1%woody_litter     - scalar
    structOut%fungi            = struct1%fungi            - scalar
    structOut%mycorrhiza       = struct1%mycorrhiza       - scalar
    structOut%microbial        = struct1%microbial        - scalar
    structOut%residue          = struct1%residue          - scalar
    structOut%residue_assoc    = struct1%residue_assoc    - scalar
  
  END FUNCTION Subtract_sb_comp_p_minus_scalar
  !--------------------------------------------------------------------------------------------------------------------------------
  ! sb_components_pool - field2d(:,:)        | ON_CHUNK
  !--------------------------------------------------------------------------------------------------------------------------------
  FUNCTION Subtract_sb_comp_p_minus_field2d(struct1, field2d) RESULT(structOut)
  
    CLASS(t_sb_components_pool),      INTENT(in)  :: struct1
    REAL(wp),                         INTENT(in)  :: field2d(:,:)
    TYPE(t_sb_components_pool)                    :: structOut
  
    TYPE(t_subset)              :: sub  
    INTEGER                     :: nvgrid
    CHARACTER(len=*), PARAMETER :: routine = modname//':Subtract_sb_comp_p_minus_field2d'
  
    sub = struct1%Get_subset()
    IF (sub%type /= ON_CHUNK) THEN  ! get domain dimension not available via Get_subset() during process init
      sub%nc = SIZE(struct1%domain_c(:,:,:), 1)
      nvgrid = SIZE(struct1%domain_c(:,:,:), 2)
      sub%nb = SIZE(struct1%domain_c(:,:,:), 3)
    ELSE
      nvgrid = SIZE(struct1%carbon(:,:), 2)
    ENDIF
    CALL structOut%Init_op_ov_hlp('op_ov_hlp', 'op_ov_hlp', struct1%owner_model_id, &
      &                           sub%type, sub%nc, nvgrid, sub%nb, .TRUE., TRIM(struct1%element_unit), 'veg_p')
  
    ! run calculation
    IF (sub%type == ON_CHUNK) THEN
      structOut%carbon          = struct1%carbon      - field2d(:,:)
      structOut%nitrogen        = struct1%nitrogen    - field2d(:,:)
      structOut%phosphorus      = struct1%phosphorus  - field2d(:,:)
      structOut%carbon13        = struct1%carbon13    - field2d(:,:)
      structOut%carbon14        = struct1%carbon14    - field2d(:,:)
      structOut%nitrogen15      = struct1%nitrogen15  - field2d(:,:)
    ELSE
      CALL finish(routine, 'Dimension mismatch, cannot calc with 2d-field and jsb_var real3d on domain')
    END IF

    structOut%dom              = struct1%dom              - field2d(:,:)
    structOut%dom_assoc        = struct1%dom_assoc        - field2d(:,:)
    structOut%soluable_litter  = struct1%soluable_litter  - field2d(:,:)
    structOut%polymeric_litter = struct1%polymeric_litter - field2d(:,:)
    structOut%woody_litter     = struct1%woody_litter     - field2d(:,:)
    structOut%fungi            = struct1%fungi            - field2d(:,:)
    structOut%mycorrhiza       = struct1%mycorrhiza       - field2d(:,:)
    structOut%microbial        = struct1%microbial        - field2d(:,:)
    structOut%residue          = struct1%residue          - field2d(:,:)
    structOut%residue_assoc    = struct1%residue_assoc    - field2d(:,:)
  
  END FUNCTION Subtract_sb_comp_p_minus_field2d
  !--------------------------------------------------------------------------------------------------------------------------------
  ! sb_components_pool - field3d(:,:,:)        | ON_CHUNK & ON_DOMAIN | field3d(:,:,1) & field3d(:,:,:)
  !--------------------------------------------------------------------------------------------------------------------------------
  FUNCTION Subtract_sb_comp_p_minus_field3d(struct1, field3d) RESULT(structOut)
  
    CLASS(t_sb_components_pool),      INTENT(in)  :: struct1
    REAL(wp),                         INTENT(in)  :: field3d(:,:,:)
    TYPE(t_sb_components_pool)                    :: structOut
  
    TYPE(t_subset)              :: sub  
    INTEGER                     :: nvgrid
    CHARACTER(len=*), PARAMETER :: routine = modname//':Subtract_sb_comp_p_minus_field3d'
  
    sub = struct1%Get_subset()
    IF (sub%type /= ON_CHUNK) THEN  ! get domain dimension not available via Get_subset() during process init
      sub%nc = SIZE(struct1%domain_c(:,:,:), 1)
      nvgrid = SIZE(struct1%domain_c(:,:,:), 2)
      sub%nb = SIZE(struct1%domain_c(:,:,:), 3)
    ELSE
      nvgrid = SIZE(struct1%carbon(:,:), 2)
    ENDIF
    CALL structOut%Init_op_ov_hlp('op_ov_hlp', 'op_ov_hlp', struct1%owner_model_id, &
      &                           sub%type, sub%nc, nvgrid, sub%nb, .TRUE., TRIM(struct1%element_unit), 'veg_p')
  
    ! run calculation
    IF (sub%type == ON_CHUNK) THEN
      IF (SIZE(field3d,3) > 1) CALL finish(routine, 'Dimension mismatch, field3d with nblks > 1; calc with REAL pool%element(:)')
      structOut%carbon          = struct1%carbon      - field3d(:,:,1)
      structOut%nitrogen        = struct1%nitrogen    - field3d(:,:,1)
      structOut%phosphorus      = struct1%phosphorus  - field3d(:,:,1)
      structOut%carbon13        = struct1%carbon13    - field3d(:,:,1)
      structOut%carbon14        = struct1%carbon14    - field3d(:,:,1)
      structOut%nitrogen15      = struct1%nitrogen15  - field3d(:,:,1)
    ELSE
      structOut%domain_c        = struct1%domain_c    - field3d(:,:,:)
      structOut%domain_n        = struct1%domain_n    - field3d(:,:,:)
      structOut%domain_p        = struct1%domain_p    - field3d(:,:,:)
      structOut%domain_c13      = struct1%domain_c13  - field3d(:,:,:)
      structOut%domain_c14      = struct1%domain_c14  - field3d(:,:,:)
      structOut%domain_n15      = struct1%domain_n15  - field3d(:,:,:)
    END IF

    structOut%dom              = struct1%dom              - field3d(:,:,:)
    structOut%dom_assoc        = struct1%dom_assoc        - field3d(:,:,:)
    structOut%soluable_litter  = struct1%soluable_litter  - field3d(:,:,:)
    structOut%polymeric_litter = struct1%polymeric_litter - field3d(:,:,:)
    structOut%woody_litter     = struct1%woody_litter     - field3d(:,:,:)
    structOut%fungi            = struct1%fungi            - field3d(:,:,:)
    structOut%mycorrhiza       = struct1%mycorrhiza       - field3d(:,:,:)
    structOut%microbial        = struct1%microbial        - field3d(:,:,:)
    structOut%residue          = struct1%residue          - field3d(:,:,:)
    structOut%residue_assoc    = struct1%residue_assoc    - field3d(:,:,:)
  
  END FUNCTION Subtract_sb_comp_p_minus_field3d
  !--------------------------------------------------------------------------------------------------------------------------------
  ! sb_components_pool - sb_components_pool        | ON_CHUNK & ON_DOMAIN
  !--------------------------------------------------------------------------------------------------------------------------------
  FUNCTION Subtract_sb_comp_p_minus_sb_comp_p(struct1, struct2) RESULT(structOut)
  
    CLASS(t_sb_components_pool),      INTENT(in)  :: struct1, struct2
    TYPE(t_sb_components_pool)                    :: structOut
  
    TYPE(t_subset)              :: sub  
    INTEGER                     :: nvgrid
    CHARACTER(len=*), PARAMETER :: routine = modname//':Subtract_sb_comp_p_minus_sb_comp_p'
  
    sub = struct1%Get_subset()
    IF (sub%type /= ON_CHUNK) THEN  ! get domain dimension not available via Get_subset() during process init
      sub%nc = SIZE(struct1%domain_c(:,:,:), 1)
      nvgrid = SIZE(struct1%domain_c(:,:,:), 2)
      sub%nb = SIZE(struct1%domain_c(:,:,:), 3)
    ELSE
      nvgrid = SIZE(struct1%carbon(:,:), 2)
    ENDIF
    CALL structOut%Init_op_ov_hlp('op_ov_hlp', 'op_ov_hlp', struct1%owner_model_id, &
      &                           sub%type, sub%nc, nvgrid, sub%nb, .TRUE., TRIM(struct1%element_unit), 'veg_p')
  
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

    structOut%dom              = struct1%dom              - struct2%dom             
    structOut%dom_assoc        = struct1%dom_assoc        - struct2%dom_assoc       
    structOut%soluable_litter  = struct1%soluable_litter  - struct2%soluable_litter 
    structOut%polymeric_litter = struct1%polymeric_litter - struct2%polymeric_litter
    structOut%woody_litter     = struct1%woody_litter     - struct2%woody_litter    
    structOut%fungi            = struct1%fungi            - struct2%fungi           
    structOut%mycorrhiza       = struct1%mycorrhiza       - struct2%mycorrhiza      
    structOut%microbial        = struct1%microbial        - struct2%microbial       
    structOut%residue          = struct1%residue          - struct2%residue         
    structOut%residue_assoc    = struct1%residue_assoc    - struct2%residue_assoc   
  
  END FUNCTION Subtract_sb_comp_p_minus_sb_comp_p
  
  !--------------------------------------------------------------------------------------------------------------------------------
  ! * multiply
  !--------------------------------------------------------------------------------------------------------------------------------
  !--------------------------------------------------------------------------------------------------------------------------------
  ! sb_components_pool * scalar        | ON_CHUNK & ON_DOMAIN
  !--------------------------------------------------------------------------------------------------------------------------------
  FUNCTION Multiply_sb_comp_p_with_scalar(struct1, scalar) RESULT(structOut)
  
    CLASS(t_sb_components_pool),      INTENT(in)  :: struct1
    REAL(wp),                         INTENT(in)  :: scalar
    TYPE(t_sb_components_pool)                    :: structOut

    TYPE(t_subset)              :: sub  
    INTEGER                     :: nvgrid
    CHARACTER(len=*), PARAMETER :: routine = modname//':Multiply_sb_comp_p_with_scalar'

    sub = struct1%Get_subset()
    IF (sub%type /= ON_CHUNK) THEN  ! get domain dimension not available via Get_subset() during process init
      sub%nc = SIZE(struct1%domain_c(:,:,:), 1)
      nvgrid = SIZE(struct1%domain_c(:,:,:), 2)
      sub%nb = SIZE(struct1%domain_c(:,:,:), 3)
    ELSE
      nvgrid = SIZE(struct1%carbon(:,:), 2)
    ENDIF
    CALL structOut%Init_op_ov_hlp('op_ov_hlp', 'op_ov_hlp', struct1%owner_model_id, &
      &                           sub%type, sub%nc, nvgrid, sub%nb, .TRUE., TRIM(struct1%element_unit), 'veg_p')

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

    structOut%dom              = struct1%dom              * scalar
    structOut%dom_assoc        = struct1%dom_assoc        * scalar
    structOut%soluable_litter  = struct1%soluable_litter  * scalar
    structOut%polymeric_litter = struct1%polymeric_litter * scalar
    structOut%woody_litter     = struct1%woody_litter     * scalar
    structOut%fungi            = struct1%fungi            * scalar
    structOut%mycorrhiza       = struct1%mycorrhiza       * scalar
    structOut%microbial        = struct1%microbial        * scalar
    structOut%residue          = struct1%residue          * scalar
    structOut%residue_assoc    = struct1%residue_assoc    * scalar
  
  END FUNCTION Multiply_sb_comp_p_with_scalar
  !--------------------------------------------------------------------------------------------------------------------------------
  ! scalar * sb_components_pool        | ON_CHUNK & ON_DOMAIN
  !--------------------------------------------------------------------------------------------------------------------------------
  FUNCTION Multiply_scalar_with_sb_comp_p(scalar, struct1) RESULT(structOut)
  
    REAL(wp),                         INTENT(in)  :: scalar
    CLASS(t_sb_components_pool),      INTENT(in)  :: struct1
    TYPE(t_sb_components_pool)                    :: structOut
  
    TYPE(t_subset)              :: sub  
    INTEGER                     :: nvgrid
    CHARACTER(len=*), PARAMETER :: routine = modname//':Multiply_scalar_with_sb_comp_p'
  
    sub = struct1%Get_subset()
    IF (sub%type /= ON_CHUNK) THEN  ! get domain dimension not available via Get_subset() during process init
      sub%nc = SIZE(struct1%domain_c(:,:,:), 1)
      nvgrid = SIZE(struct1%domain_c(:,:,:), 2)
      sub%nb = SIZE(struct1%domain_c(:,:,:), 3)
    ELSE
      nvgrid = SIZE(struct1%carbon(:,:), 2)
    ENDIF
    CALL structOut%Init_op_ov_hlp('op_ov_hlp', 'op_ov_hlp', struct1%owner_model_id, &
      &                           sub%type, sub%nc, nvgrid, sub%nb, .TRUE., TRIM(struct1%element_unit), 'veg_p')

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

    structOut%dom              = scalar * struct1%dom             
    structOut%dom_assoc        = scalar * struct1%dom_assoc       
    structOut%soluable_litter  = scalar * struct1%soluable_litter 
    structOut%polymeric_litter = scalar * struct1%polymeric_litter
    structOut%woody_litter     = scalar * struct1%woody_litter    
    structOut%fungi            = scalar * struct1%fungi           
    structOut%mycorrhiza       = scalar * struct1%mycorrhiza      
    structOut%microbial        = scalar * struct1%microbial       
    structOut%residue          = scalar * struct1%residue         
    structOut%residue_assoc    = scalar * struct1%residue_assoc   
  
  END FUNCTION Multiply_scalar_with_sb_comp_p
  !--------------------------------------------------------------------------------------------------------------------------------
  ! sb_components_pool * field2d(:,:)        | ON_CHUNK
  !--------------------------------------------------------------------------------------------------------------------------------
  FUNCTION Multiply_sb_comp_p_with_field2d(struct1, field2d) RESULT(structOut)
  
    CLASS(t_sb_components_pool),      INTENT(in)  :: struct1
    REAL(wp),                         INTENT(in)  :: field2d(:,:)
    TYPE(t_sb_components_pool)                    :: structOut
  
    TYPE(t_subset)              :: sub  
    INTEGER                     :: nvgrid
    CHARACTER(len=*), PARAMETER :: routine = modname//':Multiply_sb_comp_p_with_field2d'
  
    sub = struct1%Get_subset()
    IF (sub%type /= ON_CHUNK) THEN  ! get domain dimension not available via Get_subset() during process init
      sub%nc = SIZE(struct1%domain_c(:,:,:), 1)
      nvgrid = SIZE(struct1%domain_c(:,:,:), 2)
      sub%nb = SIZE(struct1%domain_c(:,:,:), 3)
    ELSE
      nvgrid = SIZE(struct1%carbon(:,:), 2)
    ENDIF
    CALL structOut%Init_op_ov_hlp('op_ov_hlp', 'op_ov_hlp', struct1%owner_model_id, &
      &                           sub%type, sub%nc, nvgrid, sub%nb, .TRUE., TRIM(struct1%element_unit), 'veg_p')
  
    ! run calculation
    IF (sub%type == ON_CHUNK) THEN
      structOut%carbon          = struct1%carbon      * field2d(:,:)
      structOut%nitrogen        = struct1%nitrogen    * field2d(:,:)
      structOut%phosphorus      = struct1%phosphorus  * field2d(:,:)
      structOut%carbon13        = struct1%carbon13    * field2d(:,:)
      structOut%carbon14        = struct1%carbon14    * field2d(:,:)
      structOut%nitrogen15      = struct1%nitrogen15  * field2d(:,:)
    ELSE
      CALL finish(routine, 'Dimension mismatch, cannot calc with 2d-field and jsb_var real3d on domain')
    END IF

    structOut%dom              = struct1%dom              * field2d(:,:)
    structOut%dom_assoc        = struct1%dom_assoc        * field2d(:,:)
    structOut%soluable_litter  = struct1%soluable_litter  * field2d(:,:)
    structOut%polymeric_litter = struct1%polymeric_litter * field2d(:,:)
    structOut%woody_litter     = struct1%woody_litter     * field2d(:,:)
    structOut%fungi            = struct1%fungi            * field2d(:,:)
    structOut%mycorrhiza       = struct1%mycorrhiza       * field2d(:,:)
    structOut%microbial        = struct1%microbial        * field2d(:,:)
    structOut%residue          = struct1%residue          * field2d(:,:)
    structOut%residue_assoc    = struct1%residue_assoc    * field2d(:,:)
  
  END FUNCTION Multiply_sb_comp_p_with_field2d
  !--------------------------------------------------------------------------------------------------------------------------------
  ! field2d(:,:) * sb_components_pool        | ON_CHUNK
  !--------------------------------------------------------------------------------------------------------------------------------
  FUNCTION Multiply_field2d_with_sb_comp_p(field2d, struct1) RESULT(structOut)
  
    REAL(wp),                         INTENT(in)  :: field2d(:,:)
    CLASS(t_sb_components_pool),      INTENT(in)  :: struct1
    TYPE(t_sb_components_pool)                    :: structOut
  
    TYPE(t_subset)              :: sub  
    INTEGER                     :: nvgrid
    CHARACTER(len=*), PARAMETER :: routine = modname//':Multiply_field2d_with_sb_comp_p'
  
    sub = struct1%Get_subset()
    IF (sub%type /= ON_CHUNK) THEN  ! get domain dimension not available via Get_subset() during process init
      sub%nc = SIZE(struct1%domain_c(:,:,:), 1)
      nvgrid = SIZE(struct1%domain_c(:,:,:), 2)
      sub%nb = SIZE(struct1%domain_c(:,:,:), 3)
    ELSE
      nvgrid = SIZE(struct1%carbon(:,:), 2)
    ENDIF
    CALL structOut%Init_op_ov_hlp('op_ov_hlp', 'op_ov_hlp', struct1%owner_model_id, &
      &                           sub%type, sub%nc, nvgrid, sub%nb, .TRUE., TRIM(struct1%element_unit), 'veg_p')
  
    ! run calculation
    IF (sub%type == ON_CHUNK) THEN
      structOut%carbon          = field2d(:,:) * struct1%carbon    
      structOut%nitrogen        = field2d(:,:) * struct1%nitrogen  
      structOut%phosphorus      = field2d(:,:) * struct1%phosphorus
      structOut%carbon13        = field2d(:,:) * struct1%carbon13  
      structOut%carbon14        = field2d(:,:) * struct1%carbon14  
      structOut%nitrogen15      = field2d(:,:) * struct1%nitrogen15
    ELSE
      CALL finish(routine, 'Dimension mismatch, cannot calc with 2d-field and jsb_var real3d on domain')
    END IF

    structOut%dom              = field2d(:,:) * struct1%dom             
    structOut%dom_assoc        = field2d(:,:) * struct1%dom_assoc       
    structOut%soluable_litter  = field2d(:,:) * struct1%soluable_litter 
    structOut%polymeric_litter = field2d(:,:) * struct1%polymeric_litter
    structOut%woody_litter     = field2d(:,:) * struct1%woody_litter    
    structOut%fungi            = field2d(:,:) * struct1%fungi           
    structOut%mycorrhiza       = field2d(:,:) * struct1%mycorrhiza      
    structOut%microbial        = field2d(:,:) * struct1%microbial       
    structOut%residue          = field2d(:,:) * struct1%residue         
    structOut%residue_assoc    = field2d(:,:) * struct1%residue_assoc   
  
  END FUNCTION Multiply_field2d_with_sb_comp_p
  !--------------------------------------------------------------------------------------------------------------------------------
  ! sb_components_pool * field3d(:,:,:)        | ON_CHUNK & ON_DOMAIN | field3d(:,:,1) & field3d(:,:,:)
  !--------------------------------------------------------------------------------------------------------------------------------
  FUNCTION Multiply_sb_comp_p_with_field3d(struct1, field3d) RESULT(structOut)
  
    CLASS(t_sb_components_pool),      INTENT(in)  :: struct1
    REAL(wp),                         INTENT(in)  :: field3d(:,:,:)
    TYPE(t_sb_components_pool)                    :: structOut
  
    TYPE(t_subset)              :: sub  
    INTEGER                     :: nvgrid
    CHARACTER(len=*), PARAMETER :: routine = modname//':Multiply_sb_comp_p_with_field3d'
  
    sub = struct1%Get_subset()
    IF (sub%type /= ON_CHUNK) THEN  ! get domain dimension not available via Get_subset() during process init
      sub%nc = SIZE(struct1%domain_c(:,:,:), 1)
      nvgrid = SIZE(struct1%domain_c(:,:,:), 2)
      sub%nb = SIZE(struct1%domain_c(:,:,:), 3)
    ELSE
      nvgrid = SIZE(struct1%carbon(:,:), 2)
    ENDIF
    CALL structOut%Init_op_ov_hlp('op_ov_hlp', 'op_ov_hlp', struct1%owner_model_id, &
      &                           sub%type, sub%nc, nvgrid, sub%nb, .TRUE., TRIM(struct1%element_unit), 'veg_p')
  
    ! run calculation
    IF (sub%type == ON_CHUNK) THEN
      IF (SIZE(field3d,3) > 1) CALL finish(routine, 'Dimension mismatch, field3d with nblks > 1; calc with REAL pool%element(:)')
      structOut%carbon          = struct1%carbon      * field3d(:,:,1)
      structOut%nitrogen        = struct1%nitrogen    * field3d(:,:,1)
      structOut%phosphorus      = struct1%phosphorus  * field3d(:,:,1)
      structOut%carbon13        = struct1%carbon13    * field3d(:,:,1)
      structOut%carbon14        = struct1%carbon14    * field3d(:,:,1)
      structOut%nitrogen15      = struct1%nitrogen15  * field3d(:,:,1)
    ELSE
      structOut%domain_c        = struct1%domain_c    * field3d(:,:,:)
      structOut%domain_n        = struct1%domain_n    * field3d(:,:,:)
      structOut%domain_p        = struct1%domain_p    * field3d(:,:,:)
      structOut%domain_c13      = struct1%domain_c13  * field3d(:,:,:)
      structOut%domain_c14      = struct1%domain_c14  * field3d(:,:,:)
      structOut%domain_n15      = struct1%domain_n15  * field3d(:,:,:)
    END IF

    structOut%dom              = struct1%dom              * field3d(:,:,:)
    structOut%dom_assoc        = struct1%dom_assoc        * field3d(:,:,:)
    structOut%soluable_litter  = struct1%soluable_litter  * field3d(:,:,:)
    structOut%polymeric_litter = struct1%polymeric_litter * field3d(:,:,:)
    structOut%woody_litter     = struct1%woody_litter     * field3d(:,:,:)
    structOut%fungi            = struct1%fungi            * field3d(:,:,:)
    structOut%mycorrhiza       = struct1%mycorrhiza       * field3d(:,:,:)
    structOut%microbial        = struct1%microbial        * field3d(:,:,:)
    structOut%residue          = struct1%residue          * field3d(:,:,:)
    structOut%residue_assoc    = struct1%residue_assoc    * field3d(:,:,:)
  
  END FUNCTION Multiply_sb_comp_p_with_field3d
  !--------------------------------------------------------------------------------------------------------------------------------
  ! field3d(:,:,:) * sb_components_pool        | ON_CHUNK & ON_DOMAIN | field3d(:,:,1) & field3d(:,:,:)
  !--------------------------------------------------------------------------------------------------------------------------------
  FUNCTION Multiply_field3d_with_sb_comp_p(field3d, struct1) RESULT(structOut)
  
    REAL(wp),                         INTENT(in)  :: field3d(:,:,:)
    CLASS(t_sb_components_pool),      INTENT(in)  :: struct1
    TYPE(t_sb_components_pool)                    :: structOut
  
    TYPE(t_subset)              :: sub
    INTEGER                     :: nvgrid
    CHARACTER(len=*), PARAMETER :: routine = modname//':Multiply_field3d_with_sb_comp_p'
  
    sub = struct1%Get_subset()
    IF (sub%type /= ON_CHUNK) THEN  ! get domain dimension not available via Get_subset() during process init
      sub%nc = SIZE(struct1%domain_c(:,:,:), 1)
      nvgrid = SIZE(struct1%domain_c(:,:,:), 2)
      sub%nb = SIZE(struct1%domain_c(:,:,:), 3)
    ELSE
      nvgrid = SIZE(struct1%carbon(:,:), 2)
    ENDIF
    CALL structOut%Init_op_ov_hlp('op_ov_hlp', 'op_ov_hlp', struct1%owner_model_id, &
      &                           sub%type, sub%nc, nvgrid, sub%nb, .TRUE., TRIM(struct1%element_unit), 'veg_p')
  
    ! run calculation
    IF (sub%type == ON_CHUNK) THEN
      IF (SIZE(field3d,3) > 1) CALL finish(routine, 'Dimension mismatch, field3d with nblks > 1; calc with REAL pool%element(:)')
      structOut%carbon          = field3d(:,:,1) * struct1%carbon    
      structOut%nitrogen        = field3d(:,:,1) * struct1%nitrogen  
      structOut%phosphorus      = field3d(:,:,1) * struct1%phosphorus
      structOut%carbon13        = field3d(:,:,1) * struct1%carbon13  
      structOut%carbon14        = field3d(:,:,1) * struct1%carbon14  
      structOut%nitrogen15      = field3d(:,:,1) * struct1%nitrogen15
    ELSE
      structOut%domain_c        = field3d(:,:,:) * struct1%domain_c  
      structOut%domain_n        = field3d(:,:,:) * struct1%domain_n  
      structOut%domain_p        = field3d(:,:,:) * struct1%domain_p  
      structOut%domain_c13      = field3d(:,:,:) * struct1%domain_c13
      structOut%domain_c14      = field3d(:,:,:) * struct1%domain_c14
      structOut%domain_n15      = field3d(:,:,:) * struct1%domain_n15
    END IF

    structOut%dom              = field3d(:,:,:) * struct1%dom             
    structOut%dom_assoc        = field3d(:,:,:) * struct1%dom_assoc       
    structOut%soluable_litter  = field3d(:,:,:) * struct1%soluable_litter 
    structOut%polymeric_litter = field3d(:,:,:) * struct1%polymeric_litter
    structOut%woody_litter     = field3d(:,:,:) * struct1%woody_litter    
    structOut%fungi            = field3d(:,:,:) * struct1%fungi           
    structOut%mycorrhiza       = field3d(:,:,:) * struct1%mycorrhiza      
    structOut%microbial        = field3d(:,:,:) * struct1%microbial       
    structOut%residue          = field3d(:,:,:) * struct1%residue         
    structOut%residue_assoc    = field3d(:,:,:) * struct1%residue_assoc   
  
  END FUNCTION Multiply_field3d_with_sb_comp_p
  !--------------------------------------------------------------------------------------------------------------------------------
  ! sb_components_pool * sb_components_pool        | ON_CHUNK & ON_DOMAIN
  !--------------------------------------------------------------------------------------------------------------------------------
  FUNCTION Multiply_sb_comp_p_with_sb_comp_p(struct1, struct2) RESULT(structOut)
  
    CLASS(t_sb_components_pool),      INTENT(in)  :: struct1, struct2
    TYPE(t_sb_components_pool)                    :: structOut
  
    TYPE(t_subset)              :: sub  
    INTEGER                     :: nvgrid
    CHARACTER(len=*), PARAMETER :: routine = modname//':Multiply_sb_comp_p_with_sb_comp_p'
  
    sub = struct1%Get_subset()
    IF (sub%type /= ON_CHUNK) THEN  ! get domain dimension not available via Get_subset() during process init
      sub%nc = SIZE(struct1%domain_c(:,:,:), 1)
      nvgrid = SIZE(struct1%domain_c(:,:,:), 2)
      sub%nb = SIZE(struct1%domain_c(:,:,:), 3)
    ELSE
      nvgrid = SIZE(struct1%carbon(:,:), 2)
    ENDIF
    CALL structOut%Init_op_ov_hlp('op_ov_hlp', 'op_ov_hlp', struct1%owner_model_id, &
      &                           sub%type, sub%nc, nvgrid, sub%nb, .TRUE., TRIM(struct1%element_unit), 'veg_p')
  
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

    structOut%dom              = struct1%dom              * struct2%dom             
    structOut%dom_assoc        = struct1%dom_assoc        * struct2%dom_assoc       
    structOut%soluable_litter  = struct1%soluable_litter  * struct2%soluable_litter 
    structOut%polymeric_litter = struct1%polymeric_litter * struct2%polymeric_litter
    structOut%woody_litter     = struct1%woody_litter     * struct2%woody_litter    
    structOut%fungi            = struct1%fungi            * struct2%fungi           
    structOut%mycorrhiza       = struct1%mycorrhiza       * struct2%mycorrhiza      
    structOut%microbial        = struct1%microbial        * struct2%microbial       
    structOut%residue          = struct1%residue          * struct2%residue         
    structOut%residue_assoc    = struct1%residue_assoc    * struct2%residue_assoc   
  
  END FUNCTION Multiply_sb_comp_p_with_sb_comp_p
  
  !--------------------------------------------------------------------------------------------------------------------------------
  ! / divide
  !--------------------------------------------------------------------------------------------------------------------------------
  !--------------------------------------------------------------------------------------------------------------------------------
  ! sb_components_pool / scalar        | ON_CHUNK & ON_DOMAIN
  !--------------------------------------------------------------------------------------------------------------------------------
  FUNCTION Div_sb_comp_p_by_scalar(struct1, scalar) RESULT(structOut)
  
    CLASS(t_sb_components_pool),      INTENT(in)  :: struct1
    REAL(wp),                         INTENT(in)  :: scalar
    TYPE(t_sb_components_pool)                    :: structOut

    TYPE(t_subset)              :: sub  
    INTEGER                     :: nvgrid
    CHARACTER(len=*), PARAMETER :: routine = modname//':Div_sb_comp_p_by_scalar'

    sub = struct1%Get_subset()
    IF (sub%type /= ON_CHUNK) THEN  ! get domain dimension not available via Get_subset() during process init
      sub%nc = SIZE(struct1%domain_c(:,:,:), 1)
      nvgrid = SIZE(struct1%domain_c(:,:,:), 2)
      sub%nb = SIZE(struct1%domain_c(:,:,:), 3)
    ELSE
      nvgrid = SIZE(struct1%carbon(:,:), 2)
    ENDIF
    CALL structOut%Init_op_ov_hlp('op_ov_hlp', 'op_ov_hlp', struct1%owner_model_id, &
      &                           sub%type, sub%nc, nvgrid, sub%nb, .TRUE., TRIM(struct1%element_unit), 'veg_p')

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

    structOut%dom              = struct1%dom              / scalar
    structOut%dom_assoc        = struct1%dom_assoc        / scalar
    structOut%soluable_litter  = struct1%soluable_litter  / scalar
    structOut%polymeric_litter = struct1%polymeric_litter / scalar
    structOut%woody_litter     = struct1%woody_litter     / scalar
    structOut%fungi            = struct1%fungi            / scalar
    structOut%mycorrhiza       = struct1%mycorrhiza       / scalar
    structOut%microbial        = struct1%microbial        / scalar
    structOut%residue          = struct1%residue          / scalar
    structOut%residue_assoc    = struct1%residue_assoc    / scalar
  
  END FUNCTION Div_sb_comp_p_by_scalar
  !--------------------------------------------------------------------------------------------------------------------------------
  ! sb_components_pool / field2d(:,:)        | ON_CHUNK
  !--------------------------------------------------------------------------------------------------------------------------------
  FUNCTION Div_sb_comp_p_by_field2d(struct1, field2d) RESULT(structOut)
  
    CLASS(t_sb_components_pool),      INTENT(in)  :: struct1
    REAL(wp),                         INTENT(in)  :: field2d(:,:)
    TYPE(t_sb_components_pool)                    :: structOut
  
    TYPE(t_subset)              :: sub  
    INTEGER                     :: nvgrid
    CHARACTER(len=*), PARAMETER :: routine = modname//':Div_sb_comp_p_by_field2d'
  
    sub = struct1%Get_subset()
    IF (sub%type /= ON_CHUNK) THEN  ! get domain dimension not available via Get_subset() during process init
      sub%nc = SIZE(struct1%domain_c(:,:,:), 1)
      nvgrid = SIZE(struct1%domain_c(:,:,:), 2)
      sub%nb = SIZE(struct1%domain_c(:,:,:), 3)
    ELSE
      nvgrid = SIZE(struct1%carbon(:,:), 2)
    ENDIF
    CALL structOut%Init_op_ov_hlp('op_ov_hlp', 'op_ov_hlp', struct1%owner_model_id, &
      &                           sub%type, sub%nc, nvgrid, sub%nb, .TRUE., TRIM(struct1%element_unit), 'veg_p')
  
    ! run calculation
    IF (sub%type == ON_CHUNK) THEN
      structOut%carbon          = struct1%carbon      / field2d(:,:)
      structOut%nitrogen        = struct1%nitrogen    / field2d(:,:)
      structOut%phosphorus      = struct1%phosphorus  / field2d(:,:)
      structOut%carbon13        = struct1%carbon13    / field2d(:,:)
      structOut%carbon14        = struct1%carbon14    / field2d(:,:)
      structOut%nitrogen15      = struct1%nitrogen15  / field2d(:,:)
    ELSE
      CALL finish(routine, 'Dimension mismatch, cannot calc with 2d-field and jsb_var real3d on domain')
    END IF

    structOut%dom              = struct1%dom              / field2d(:,:)
    structOut%dom_assoc        = struct1%dom_assoc        / field2d(:,:)
    structOut%soluable_litter  = struct1%soluable_litter  / field2d(:,:)
    structOut%polymeric_litter = struct1%polymeric_litter / field2d(:,:)
    structOut%woody_litter     = struct1%woody_litter     / field2d(:,:)
    structOut%fungi            = struct1%fungi            / field2d(:,:)
    structOut%mycorrhiza       = struct1%mycorrhiza       / field2d(:,:)
    structOut%microbial        = struct1%microbial        / field2d(:,:)
    structOut%residue          = struct1%residue          / field2d(:,:)
    structOut%residue_assoc    = struct1%residue_assoc    / field2d(:,:)
  
  END FUNCTION Div_sb_comp_p_by_field2d
  !--------------------------------------------------------------------------------------------------------------------------------
  ! sb_components_pool / field3d(:,:,:)        | ON_CHUNK & ON_DOMAIN | field3d(:,:,1) & field3d(:,:,:)
  !--------------------------------------------------------------------------------------------------------------------------------
  FUNCTION Div_sb_comp_p_by_field3d(struct1, field3d) RESULT(structOut)
  
    CLASS(t_sb_components_pool),      INTENT(in)  :: struct1
    REAL(wp),                         INTENT(in)  :: field3d(:,:,:)
    TYPE(t_sb_components_pool)                    :: structOut
  
    TYPE(t_subset)              :: sub  
    INTEGER                     :: nvgrid
    CHARACTER(len=*), PARAMETER :: routine = modname//':Div_sb_comp_p_by_field3d'
  
    sub = struct1%Get_subset()
    IF (sub%type /= ON_CHUNK) THEN  ! get domain dimension not available via Get_subset() during process init
      sub%nc = SIZE(struct1%domain_c(:,:,:), 1)
      nvgrid = SIZE(struct1%domain_c(:,:,:), 2)
      sub%nb = SIZE(struct1%domain_c(:,:,:), 3)
    ELSE
      nvgrid = SIZE(struct1%carbon(:,:), 2)
    ENDIF
    CALL structOut%Init_op_ov_hlp('op_ov_hlp', 'op_ov_hlp', struct1%owner_model_id, &
      &                           sub%type, sub%nc, nvgrid, sub%nb, .TRUE., TRIM(struct1%element_unit), 'veg_p')
  
    ! run calculation
    IF (sub%type == ON_CHUNK) THEN
      IF (SIZE(field3d,3) > 1) CALL finish(routine, 'Dimension mismatch, field3d with nblks > 1; calc with REAL pool%element(:)')
      structOut%carbon          = struct1%carbon      / field3d(:,:,1)
      structOut%nitrogen        = struct1%nitrogen    / field3d(:,:,1)
      structOut%phosphorus      = struct1%phosphorus  / field3d(:,:,1)
      structOut%carbon13        = struct1%carbon13    / field3d(:,:,1)
      structOut%carbon14        = struct1%carbon14    / field3d(:,:,1)
      structOut%nitrogen15      = struct1%nitrogen15  / field3d(:,:,1)
    ELSE
      structOut%domain_c        = struct1%domain_c    / field3d(:,:,:)
      structOut%domain_n        = struct1%domain_n    / field3d(:,:,:)
      structOut%domain_p        = struct1%domain_p    / field3d(:,:,:)
      structOut%domain_c13      = struct1%domain_c13  / field3d(:,:,:)
      structOut%domain_c14      = struct1%domain_c14  / field3d(:,:,:)
      structOut%domain_n15      = struct1%domain_n15  / field3d(:,:,:)
    END IF

    structOut%dom              = struct1%dom              / field3d(:,:,:)
    structOut%dom_assoc        = struct1%dom_assoc        / field3d(:,:,:)
    structOut%soluable_litter  = struct1%soluable_litter  / field3d(:,:,:)
    structOut%polymeric_litter = struct1%polymeric_litter / field3d(:,:,:)
    structOut%woody_litter     = struct1%woody_litter     / field3d(:,:,:)
    structOut%fungi            = struct1%fungi            / field3d(:,:,:)
    structOut%mycorrhiza       = struct1%mycorrhiza       / field3d(:,:,:)
    structOut%microbial        = struct1%microbial        / field3d(:,:,:)
    structOut%residue          = struct1%residue          / field3d(:,:,:)
    structOut%residue_assoc    = struct1%residue_assoc    / field3d(:,:,:)
  
  END FUNCTION Div_sb_comp_p_by_field3d
  !--------------------------------------------------------------------------------------------------------------------------------
  ! sb_components_pool / sb_components_pool        | ON_CHUNK & ON_DOMAIN
  !--------------------------------------------------------------------------------------------------------------------------------
  FUNCTION Div_sb_comp_p_by_sb_comp_p(struct1, struct2) RESULT(structOut)
  
    CLASS(t_sb_components_pool),      INTENT(in)  :: struct1, struct2
    TYPE(t_sb_components_pool)                    :: structOut
  
    TYPE(t_subset)              :: sub  
    INTEGER                     :: nvgrid
    CHARACTER(len=*), PARAMETER :: routine = modname//':Div_sb_comp_p_by_sb_comp_p'
  
    sub = struct1%Get_subset()
    IF (sub%type /= ON_CHUNK) THEN  ! get domain dimension not available via Get_subset() during process init
      sub%nc = SIZE(struct1%domain_c(:,:,:), 1)
      nvgrid = SIZE(struct1%domain_c(:,:,:), 2)
      sub%nb = SIZE(struct1%domain_c(:,:,:), 3)
    ELSE
      nvgrid = SIZE(struct1%carbon(:,:), 2)
    ENDIF
    CALL structOut%Init_op_ov_hlp('op_ov_hlp', 'op_ov_hlp', struct1%owner_model_id, &
      &                           sub%type, sub%nc, nvgrid, sub%nb, .TRUE., TRIM(struct1%element_unit), 'veg_p')
  
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

    structOut%dom              = struct1%dom              / struct2%dom             
    structOut%dom_assoc        = struct1%dom_assoc        / struct2%dom_assoc       
    structOut%soluable_litter  = struct1%soluable_litter  / struct2%soluable_litter 
    structOut%polymeric_litter = struct1%polymeric_litter / struct2%polymeric_litter
    structOut%woody_litter     = struct1%woody_litter     / struct2%woody_litter    
    structOut%fungi            = struct1%fungi            / struct2%fungi           
    structOut%mycorrhiza       = struct1%mycorrhiza       / struct2%mycorrhiza      
    structOut%microbial        = struct1%microbial        / struct2%microbial       
    structOut%residue          = struct1%residue          / struct2%residue         
    structOut%residue_assoc    = struct1%residue_assoc    / struct2%residue_assoc   
  
  END FUNCTION Div_sb_comp_p_by_sb_comp_p

!! #################################################################################################################################
!! code specific to QUINCY with JSBACH4 ends here
!! #################################################################################################################################
#endif

END MODULE mo_sb_memory_class_util
