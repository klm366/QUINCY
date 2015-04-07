!>
!! @brief Contains the class and methods for the JSBACH memory structure
!!
!! @author
!!  Reiner Schnur, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Reiner Schnur (2013-04-04)
!!
MODULE mo_jsb_memory_class
#ifndef __NO_JSBACH__

  USE mo_atmland_memory_class, ONLY: t_atm2land_memory, t_land2atm_memory
  USE mo_srf_memory_class,     ONLY: t_srf_memory
  USE mo_soil_memory_class,    ONLY: t_soil_memory
  USE mo_veg_memory_class,     ONLY: t_veg_memory
  !USE mo_glacier_memory_class, ONLY: t_glacier
!!$  USE mo_cdi_constants,     ONLY: DATATYPE_FLT32, DATATYPE_FLT64, DATATYPE_PACK16, &
!!$                                  ZA_SURFACE, &
!!$                                  FILETYPE_NC2, &
!!$                                  GRID_CELL, GRID_REFERENCE, GRID_UNSTRUCTURED_CELL
  USE mo_kind,              ONLY: wp
!!$  USE mo_exception,         ONLY: message

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: t_jsb_memory!, new_jsbach_memory

  ! This type contains the memory state for each tile
  ! It includes all memory types from all sub-modules
  TYPE t_jsb_memory
    REAL(wp),                POINTER :: fract     (:,:)  !< Cover fraction of this tile (relative to parent tile)
    REAL(wp),                POINTER :: fract_max (:,:)  !< Maximum potential cover fraction for this tile (rel. to parent)
    TYPE(t_atm2land_memory), POINTER :: atm2land => NULL()
    TYPE(t_land2atm_memory), POINTER :: land2atm => NULL()
    TYPE(t_srf_memory),      POINTER :: srf      => NULL()
    TYPE(t_soil_memory),     POINTER :: soil     => NULL()
    TYPE(t_veg_memory),      POINTER :: veg      => NULL()
    TYPE(t_carbon_memory),   POINTER :: carbon   => NULL()
    !TYPE(t_glacier),         POINTER :: glacier
  END type t_jsb_memory

  CHARACTER(len=*), PARAMETER :: modname = 'mo_jsb_memory_class'

#endif
END MODULE mo_jsb_memory_class
