!>
!! @brief Contains definitions and methods for landcover types.
!!
!! @author
!!  Reiner Schnur, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Reiner Schnur (2016-02-01)
!!
MODULE mo_jsb_lct_class
#ifndef __NO_JSBACH__

  USE mo_exception,           ONLY: message_text, message, finish

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: BOX_TYPE, LAND_TYPE, PFT_TYPE, BARE_TYPE, GLACIER_TYPE, LAKE_TYPE
  PUBLIC :: max_no_of_lct, t_jsb_lct, Contains_lct

  INTEGER, PARAMETER :: max_no_of_lct = 10

  ENUM, BIND(C)
    ENUMERATOR :: BOX_TYPE=1,   & !< Grid box
      &           LAND_TYPE,    & !< Land except lakes
      &           PFT_TYPE,     & !< Vegetated with PFT
      &           BARE_TYPE,    & !< Bare soil
      &           GLACIER_TYPE, & !< Glacier
      &           LAKE_TYPE       !< Lake
  END ENUM

  TYPE t_jsb_lct
    INTEGER           :: id = 0           !< Id for lct, one of enumeration above
    CHARACTER(LEN=5)  :: name = ""        !< Short name for lct
    CHARACTER(LEN=30) :: long_name = ""   !< Long name for lct
  CONTAINS
    PROCEDURE :: Get_type_name     => Get_lct_type_name
    PROCEDURE :: Get_type_longname => Get_lct_type_longname
    PROCEDURE :: Contains
  END TYPE t_jsb_lct

  INTERFACE t_jsb_lct
    PROCEDURE Construct_jsb_lct
  END INTERFACE

  CHARACTER(len=*), PARAMETER :: modname = 'mo_jsb_lct_class'

CONTAINS

  FUNCTION Construct_jsb_lct(id, name, long_name) RESULT(return_value)

    INTEGER, INTENT(in)                    :: id
    CHARACTER(LEN=*), INTENT(in), OPTIONAL :: name
    CHARACTER(LEN=*), INTENT(in), OPTIONAL :: long_name
    TYPE(t_jsb_lct)                        :: return_value

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':Construct_jsb_lct'

    return_value%id = id
    IF (PRESENT(name)) THEN
      return_value%name = TRIM(name)
    ELSE
      return_value%name = TRIM(return_value%Get_type_name())
    END IF

    IF (PRESENT(long_name)) THEN
      return_value%long_name = TRIM(long_name)
    ELSE
      return_value%long_name = TRIM(return_value%Get_type_longname())
    END IF

  END FUNCTION Construct_jsb_lct

  FUNCTION Contains(this, lct_type) RESULT(return_value)

    CLASS(t_jsb_lct), INTENT(in) :: this
    INTEGER,         INTENT(in) :: lct_type
    LOGICAL                     :: return_value

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':Contains'

    return_value = Contains_lct(this%id, lct_type)

  END FUNCTION Contains

  FUNCTION Contains_lct(lct_type_current, lct_type) RESULT(return_value)

    INTEGER,         INTENT(in) :: lct_type_current
    INTEGER,         INTENT(in) :: lct_type
    LOGICAL                     :: return_value

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':Contains_lct'

    IF (lct_type_current == lct_type) THEN
      return_value = .TRUE.
      RETURN
    END IF

    SELECT CASE(lct_type_current)
    CASE (BOX_TYPE)
      return_value = .TRUE.
    CASE (LAND_TYPE)
      SELECT CASE(lct_type)
      CASE (PFT_TYPE, BARE_TYPE, GLACIER_TYPE)
        return_value = .TRUE.
      CASE DEFAULT
        return_value = .FALSE.
      END SELECT
    CASE DEFAULT
      return_value = .FALSE.
    END SELECT

  END FUNCTION Contains_lct

  FUNCTION Get_lct_type_name(this) RESULT(return_value)

    CLASS(t_jsb_lct),  INTENT(in)  :: this
    CHARACTER(len=:), ALLOCATABLE :: return_value

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':Get_lct_type_name'

    SELECT CASE(this%id)
    CASE (BOX_TYPE)
      return_value = 'box'
    CASE (LAND_TYPE)
      return_value = 'land'
    CASE (PFT_TYPE)
      return_value = 'pft'
    CASE (BARE_TYPE)
      return_value = 'bare'
    CASE (GLACIER_TYPE)
      return_value = 'glac'
    CASE (LAKE_TYPE)
      return_value = 'lake'
    END SELECT

  END FUNCTION Get_lct_type_name

  FUNCTION Get_lct_type_longname(this) RESULT(return_value)

    CLASS(t_jsb_lct),  INTENT(in)  :: this
    CHARACTER(len=:), ALLOCATABLE :: return_value

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':Get_lct_type_longname'

    SELECT CASE(this%id)
    CASE (BOX_TYPE)
      return_value = 'Grid box'
    CASE (LAND_TYPE)
      return_value = 'Land except lakes'
    CASE (PFT_TYPE)
      return_value = 'Vegetated with PFT'
    CASE (BARE_TYPE)
      return_value = 'Bare soil'
    CASE (GLACIER_TYPE)
      return_value = 'Glacier'
    CASE (LAKE_TYPE)
      return_value = 'Lake'
    END SELECT

  END FUNCTION Get_lct_type_longname

#endif
END MODULE mo_jsb_lct_class
