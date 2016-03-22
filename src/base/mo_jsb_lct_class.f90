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
  PUBLIC :: max_no_of_lct, Get_lct_name

  INTEGER, PARAMETER :: max_no_of_lct = 10

  ENUM, BIND(C)
    ENUMERATOR :: BOX_TYPE=1, LAND_TYPE=2, PFT_TYPE, BARE_TYPE, GLACIER_TYPE, LAKE_TYPE
  END ENUM

  CHARACTER(len=*), PARAMETER :: modname = 'mo_jsb_lct_class'

CONTAINS

  FUNCTION Get_lct_name(ilct) RESULT(return_value)

    INTEGER, INTENT(in) :: ilct
    CHARACTER(len=:), ALLOCATABLE :: return_value

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':Get_lct_name'

    SELECT CASE(ilct)
    CASE (BOX_TYPE)
      return_value = 'box'
    CASE (LAND_TYPE)
      return_value = 'land'
    CASE (PFT_TYPE)
      return_value = 'pft'
    CASE (BARE_TYPE)
      return_value = 'bare_soil'
    CASE (GLACIER_TYPE)
      return_value = 'glacier'
    CASE (LAKE_TYPE)
      return_value = 'lake'
    END SELECT

  END FUNCTION Get_lct_name

#endif
END MODULE mo_jsb_lct_class
