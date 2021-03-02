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
  USE mo_jsb_var_class,          ONLY: t_jsb_var_real3d

  IMPLICIT NONE

  ! all is public except the modename; code is used by the mo_sb_memory_class only
  PUBLIC
  PRIVATE :: modname


  CHARACTER(len=*), PARAMETER :: modname = 'mo_sb_memory_class_util'

!CONTAINS



!! #################################################################################################################################
!! code specific to QUINCY with JSBACH4 ends here
!! #################################################################################################################################
#endif

END MODULE mo_sb_memory_class_util
