!>
!! @brief Contains process factory.
!!
!! @author
!!  Reiner Schnur, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Reiner Schnur (2015-06-29)
!!
MODULE mo_jsb_process_factory
#ifndef __NO_JSBACH__

  USE mo_jsb_memory_class,     ONLY: t_jsb_memory, t_jsb_memory_p
  USE mo_jsb_config_class,     ONLY: t_jsb_config
  USE mo_atmland_memory_class, ONLY: t_atm2land_memory, t_land2atm_memory
  USE mo_srf_memory_class,     ONLY: t_srf_memory
  USE mo_soil_memory_class,    ONLY: t_soil_memory
  USE mo_veg_memory_class,     ONLY: t_veg_memory
  USE mo_hd_memory_class,      ONLY: t_hd_memory
  USE mo_srf_config_class,     ONLY: t_srf_config
  USE mo_soil_config_class,    ONLY: t_soil_config
  USE mo_veg_config_class,     ONLY: t_veg_config
  USE mo_hd_config_class,      ONLY: t_hd_config

!!$  USE mo_jsb_process_class,    ONLY: t_jsb_process

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: max_no_of_processes, Create_process_config, &
            Create_process_memory, Get_process_name, Get_process_memory, Get_process_config
  PUBLIC :: SRF, ATM2LAND, LAND2ATM, SKIN, SOIL, VEG, SNOW, GLACIER, LAKE, HD

  INTEGER, PARAMETER :: max_no_of_processes = 20

  ENUM, BIND(C)
    ENUMERATOR :: SRF=1, ATM2LAND, LAND2ATM, SKIN, SOIL, VEG, SNOW, GLACIER, LAKE, HD
  END ENUM

  INTERFACE Get_process_memory
    MODULE PROCEDURE Get_atm2land_memory
    MODULE PROCEDURE Get_land2atm_memory
    MODULE PROCEDURE Get_srf_memory
    MODULE PROCEDURE Get_soil_memory
    MODULE PROCEDURE Get_veg_memory
    MODULE PROCEDURE Get_hd_memory
  END INTERFACE Get_process_memory

  INTERFACE Get_process_config
    MODULE PROCEDURE Get_srf_config
    MODULE PROCEDURE Get_soil_config
    MODULE PROCEDURE Get_veg_config
    MODULE PROCEDURE Get_hd_config
  END INTERFACE Get_process_config

CONTAINS

  FUNCTION Create_process_config(icomp, namelist_filename) RESULT(return_ptr)

    USE mo_jsb_config_class,   ONLY: t_jsb_config
    USE mo_srf_config_class,   ONLY: new_srf_config
    USE mo_soil_config_class,  ONLY: new_soil_config
    USE mo_veg_config_class,   ONLY: new_veg_config
#ifndef __NO_JSBACH_HD__
    USE mo_hd_config_class,    ONLY: new_hd_config
#endif
    USE mo_jsb_namelist_iface, ONLY: open_nml, close_nml

    INTEGER, INTENT(in) :: icomp
    CHARACTER(LEN=*), INTENT(in) :: namelist_filename
    CLASS(t_jsb_config), POINTER :: return_ptr

    INTEGER :: nml_handler

    nml_handler = open_nml(TRIM(namelist_filename))

    SELECT CASE(icomp)
    CASE (SRF)
      return_ptr  => new_srf_config(nml_handler)
    CASE (SOIL)
      return_ptr => new_soil_config(nml_handler)
    CASE (VEG)
      return_ptr  => new_veg_config(nml_handler)
#ifndef __NO_JSBACH_HD__
    CASE (HD)
      return_ptr => new_hd_config(nml_handler)
#endif
    CASE DEFAULT
      return_ptr => NULL()
    END SELECT

    CALL close_nml(nml_handler)

  END FUNCTION Create_process_config

  FUNCTION Create_process_memory(icomp) RESULT(return_ptr)

    USE mo_atmland_memory_class, ONLY: t_atm2land_memory, t_land2atm_memory, max_no_of_atmland_vars => max_no_of_vars
    USE mo_srf_memory_class,     ONLY: t_srf_memory, max_no_of_srf_vars => max_no_of_vars
    USE mo_soil_memory_class,    ONLY: t_soil_memory, max_no_of_soil_vars => max_no_of_vars
    USE mo_veg_memory_class,     ONLY: t_veg_memory, max_no_of_veg_vars => max_no_of_vars
#ifndef __NO_JSBACH_HD__
    USE mo_hd_memory_class,      ONLY: t_hd_memory, max_no_of_hd_vars => max_no_of_vars
#endif

    INTEGER, INTENT(in) :: icomp
    CLASS(t_jsb_memory), POINTER   :: return_ptr

    SELECT CASE(icomp)
    CASE (ATM2LAND)
      ALLOCATE(t_atm2land_memory::return_ptr)
      return_ptr%max_no_of_vars = max_no_of_atmland_vars
    CASE (LAND2ATM)
      ALLOCATE(t_land2atm_memory::return_ptr)
      return_ptr%max_no_of_vars = max_no_of_atmland_vars
    CASE (SRF)
      ALLOCATE(t_srf_memory::return_ptr)
      return_ptr%max_no_of_vars = max_no_of_srf_vars
    CASE (SOIL)
      ALLOCATE(t_soil_memory::return_ptr)
      return_ptr%max_no_of_vars = max_no_of_soil_vars
    CASE(VEG)
      ALLOCATE(t_veg_memory::return_ptr)
      return_ptr%max_no_of_vars = max_no_of_veg_vars
#ifndef __NO_JSBACH_HD__
    CASE (HD)
      ALLOCATE(t_hd_memory::return_ptr)
      return_ptr%max_no_of_vars = max_no_of_hd_vars
#endif
    END SELECT

    ALLOCATE(return_ptr%vars(return_ptr%max_no_of_vars))

  END FUNCTION Create_process_memory

  FUNCTION Get_process_name(icomp) RESULT(return_value)

    INTEGER, INTENT(in) :: icomp
    CHARACTER(len=:), ALLOCATABLE :: return_value

    SELECT CASE(icomp)
    CASE (ATM2LAND)
      return_value = 'atm2land'
    CASE (LAND2ATM)
      return_value = 'land2atm'
    CASE (SRF)
      return_value = 'srf'
    CASE (SOIL)
      return_value = 'soil'
    CASE (VEG)
      return_value = 'veg'
#ifndef __NO_JSBACH_HD__
    CASE (HD)
      return_value = 'hd'
#endif
    CASE (GLACIER)
      return_value = 'glac'
    CASE (LAKE)
      return_value = 'lake'
    END SELECT

  END FUNCTION Get_process_name

  SUBROUTINE Get_atm2land_memory(mem, iproc, return_ptr)

    TYPE(t_jsb_memory_p), INTENT(in) :: mem(:)
    INTEGER,              INTENT(in) :: iproc
    TYPE(t_atm2land_memory),  POINTER    :: return_ptr

    CLASS(*),  POINTER    :: ptr

    ptr => mem(iproc)%p
    SELECT TYPE(ptr)
    TYPE IS (t_atm2land_memory)
      return_ptr => ptr
    END SELECT
    NULLIFY(ptr)

  END SUBROUTINE Get_atm2land_memory

  SUBROUTINE Get_land2atm_memory(mem, iproc, return_ptr)

    TYPE(t_jsb_memory_p), INTENT(in) :: mem(:)
    INTEGER,              INTENT(in) :: iproc
    TYPE(t_land2atm_memory),  POINTER    :: return_ptr

    CLASS(*),  POINTER    :: ptr

    ptr => mem(iproc)%p
    SELECT TYPE(ptr)
    TYPE IS (t_land2atm_memory)
      return_ptr => ptr
    END SELECT
    NULLIFY(ptr)

  END SUBROUTINE Get_land2atm_memory

  SUBROUTINE Get_srf_memory(mem, iproc, return_ptr)

    TYPE(t_jsb_memory_p), INTENT(in) :: mem(:)
    INTEGER,              INTENT(in) :: iproc
    TYPE(t_srf_memory),  POINTER    :: return_ptr

    CLASS(*),  POINTER    :: ptr

    ptr => mem(iproc)%p
    SELECT TYPE(ptr)
    TYPE IS (t_srf_memory)
      return_ptr => ptr
    END SELECT
    NULLIFY(ptr)

  END SUBROUTINE Get_srf_memory

  SUBROUTINE Get_soil_memory(mem, iproc, return_ptr)

    TYPE(t_jsb_memory_p), INTENT(in) :: mem(:)
    INTEGER,              INTENT(in) :: iproc
    TYPE(t_soil_memory),  POINTER    :: return_ptr

    CLASS(*),  POINTER    :: ptr

    ptr => mem(iproc)%p
    SELECT TYPE(ptr)
    TYPE IS (t_soil_memory)
      return_ptr => ptr
    END SELECT
    NULLIFY(ptr)

  END SUBROUTINE Get_soil_memory

  SUBROUTINE Get_veg_memory(mem, iproc, return_ptr)

    TYPE(t_jsb_memory_p), INTENT(in) :: mem(:)
    INTEGER,              INTENT(in) :: iproc
    TYPE(t_veg_memory),  POINTER    :: return_ptr

    CLASS(*),  POINTER    :: ptr

    ptr => mem(iproc)%p
    SELECT TYPE(ptr)
    TYPE IS (t_veg_memory)
      return_ptr => ptr
    END SELECT
    NULLIFY(ptr)

  END SUBROUTINE Get_veg_memory

  SUBROUTINE Get_hd_memory(mem, iproc, return_ptr)

    TYPE(t_jsb_memory_p), INTENT(in) :: mem(:)
    INTEGER,              INTENT(in) :: iproc
    TYPE(t_hd_memory),  POINTER    :: return_ptr

    CLASS(*),  POINTER    :: ptr

    ptr => mem(iproc)%p
    SELECT TYPE(ptr)
    TYPE IS (t_hd_memory)
      return_ptr => ptr
    END SELECT
    NULLIFY(ptr)

  END SUBROUTINE Get_hd_memory

  SUBROUTINE Get_srf_config(config, return_ptr)

    CLASS(t_jsb_config), TARGET, INTENT(in) :: config
    TYPE(t_srf_config),  POINTER    :: return_ptr

    return_ptr => NULL()

    SELECT TYPE (config)
    TYPE IS (t_srf_config)
      return_ptr => config
    END SELECT

  END SUBROUTINE Get_srf_config

  SUBROUTINE Get_soil_config(config, return_ptr)

    CLASS(t_jsb_config), TARGET, INTENT(in) :: config
    TYPE(t_soil_config), POINTER    :: return_ptr

    return_ptr => NULL()

    SELECT TYPE (config)
    TYPE IS (t_soil_config)
      return_ptr => config
    END SELECT

  END SUBROUTINE Get_soil_config

  SUBROUTINE Get_veg_config(config, return_ptr)

    CLASS(t_jsb_config), TARGET, INTENT(in) :: config
    TYPE(t_veg_config), POINTER    :: return_ptr

    return_ptr => NULL()

    SELECT TYPE (config)
    TYPE IS (t_veg_config)
      return_ptr => config
    END SELECT

  END SUBROUTINE Get_veg_config

  SUBROUTINE Get_hd_config(config, return_ptr)

    CLASS(t_jsb_config), TARGET, INTENT(in) :: config
    TYPE(t_hd_config), POINTER    :: return_ptr

    return_ptr => NULL()

    SELECT TYPE (config)
    TYPE IS (t_hd_config)
      return_ptr => config
    END SELECT

  END SUBROUTINE Get_hd_config

#endif
END MODULE mo_jsb_process_factory
