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

  USE mo_exception,            ONLY: finish

  !USE mo_radiation_interface,  ONLY: tsk_Radiation_surface
  !
  USE mo_jsb_memory_class,     ONLY: t_jsb_memory, t_jsb_memory_p
  USE mo_jsb_config_class,     ONLY: t_jsb_config
  !
  USE mo_atmland_memory_class, ONLY: t_atm2land_memory, t_land2atm_memory
  USE mo_seb_memory_class,     ONLY: t_seb_memory
  USE mo_turb_memory_class,    ONLY: t_turb_memory
  USE mo_soil_memory_class,    ONLY: t_soil_memory
  USE mo_hydro_memory_class,   ONLY: t_hydro_memory
  USE mo_veg_memory_class,     ONLY: t_veg_memory
  USE mo_rad_memory_class,     ONLY: t_rad_memory
#ifndef __NO_JSBACH_HD__
  USE mo_hd_memory_class,      ONLY: t_hd_memory
#endif
  USE mo_seb_config_class,     ONLY: t_seb_config
  USE mo_turb_config_class,    ONLY: t_turb_config
  USE mo_soil_config_class,    ONLY: t_soil_config
  USE mo_hydro_config_class,   ONLY: t_hydro_config
  USE mo_veg_config_class,     ONLY: t_veg_config
  USE mo_rad_config_class,     ONLY: t_rad_config
#ifndef __NO_JSBACH_HD__
  USE mo_hd_config_class,      ONLY: t_hd_config
#endif

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: max_no_of_processes, & !Create_process, & !Register_process_tasks, &
            Create_process_memory, Get_process_name, Get_process_memory, Get_process_config
  PUBLIC :: ATM2LAND, LAND2ATM, SEB, TURB, SOIL, HYDRO, VEG, SNOW, HD, RAD

  INTEGER, PARAMETER :: max_no_of_processes = 20

  ENUM, BIND(C)
    ENUMERATOR :: ATM2LAND=1, LAND2ATM, SEB, TURB, SOIL, HYDRO, VEG, SNOW, HD, RAD
  END ENUM

  INTERFACE Get_process_memory
    MODULE PROCEDURE Get_atm2land_memory
    MODULE PROCEDURE Get_land2atm_memory
    MODULE PROCEDURE Get_seb_memory
    MODULE PROCEDURE Get_turb_memory
    MODULE PROCEDURE Get_soil_memory
    MODULE PROCEDURE Get_hydro_memory
    MODULE PROCEDURE Get_veg_memory
    MODULE PROCEDURE Get_rad_memory
#ifndef __NO_JSBACH_HD__
    MODULE PROCEDURE Get_hd_memory
#endif
  END INTERFACE Get_process_memory

  INTERFACE Get_process_config
    MODULE PROCEDURE Get_seb_config
    MODULE PROCEDURE Get_turb_config
    MODULE PROCEDURE Get_soil_config
    MODULE PROCEDURE Get_hydro_config
    MODULE PROCEDURE Get_veg_config
    MODULE PROCEDURE Get_rad_config
#ifndef __NO_JSBACH_HD__
    MODULE PROCEDURE Get_hd_config
#endif
  END INTERFACE Get_process_config

  CHARACTER(len=*), PARAMETER :: modname = 'mo_jsb_process_factory'

CONTAINS

  FUNCTION Create_process_memory(iproc) RESULT(return_ptr)

    USE mo_atmland_memory_class, ONLY: max_no_of_atmland_vars => max_no_of_vars
    USE mo_seb_memory_class,     ONLY: max_no_of_seb_vars => max_no_of_vars
    USE mo_turb_memory_class,    ONLY: max_no_of_turb_vars => max_no_of_vars
    USE mo_soil_memory_class,    ONLY: max_no_of_soil_vars => max_no_of_vars
    USE mo_hydro_memory_class,   ONLY: max_no_of_hydro_vars => max_no_of_vars
    USE mo_veg_memory_class,     ONLY: max_no_of_veg_vars => max_no_of_vars
    USE mo_rad_memory_class,     ONLY: max_no_of_rad_vars => max_no_of_vars
#ifndef __NO_JSBACH_HD__
    USE mo_hd_memory_class,      ONLY: max_no_of_hd_vars => max_no_of_vars
#endif

    INTEGER, INTENT(in) :: iproc
    CLASS(t_jsb_memory), POINTER   :: return_ptr

    SELECT CASE(iproc)
    CASE (ATM2LAND)
      ALLOCATE(t_atm2land_memory::return_ptr)
      return_ptr%max_no_of_vars = max_no_of_atmland_vars
    CASE (LAND2ATM)
      ALLOCATE(t_land2atm_memory::return_ptr)
      return_ptr%max_no_of_vars = max_no_of_atmland_vars
    CASE (SEB)
      ALLOCATE(t_seb_memory::return_ptr)
      return_ptr%max_no_of_vars = max_no_of_seb_vars
    CASE (TURB)
      ALLOCATE(t_turb_memory::return_ptr)
      return_ptr%max_no_of_vars = max_no_of_turb_vars
    CASE (SOIL)
      ALLOCATE(t_soil_memory::return_ptr)
      return_ptr%max_no_of_vars = max_no_of_soil_vars
    CASE (HYDRO)
      ALLOCATE(t_hydro_memory::return_ptr)
      return_ptr%max_no_of_vars = max_no_of_hydro_vars
    CASE(VEG)
      ALLOCATE(t_veg_memory::return_ptr)
      return_ptr%max_no_of_vars = max_no_of_veg_vars
    CASE(RAD)
      ALLOCATE(t_rad_memory::return_ptr)
      return_ptr%max_no_of_vars = max_no_of_rad_vars
#ifndef __NO_JSBACH_HD__
    CASE (HD)
      ALLOCATE(t_hd_memory::return_ptr)
      return_ptr%max_no_of_vars = max_no_of_hd_vars
#endif
    END SELECT

    ALLOCATE(return_ptr%vars(return_ptr%max_no_of_vars))

  END FUNCTION Create_process_memory

  FUNCTION Get_process_name(iproc) RESULT(return_value)

    INTEGER, INTENT(in) :: iproc
    CHARACTER(len=:), ALLOCATABLE :: return_value

    SELECT CASE(iproc)
    CASE (ATM2LAND)
      return_value = 'atm2land'
    CASE (LAND2ATM)
      return_value = 'land2atm'
    CASE (SEB)
      return_value = 'seb'
    CASE (TURB)
      return_value = 'turb'
    CASE (SOIL)
      return_value = 'soil'
    CASE (HYDRO)
      return_value = 'hydro'
    CASE (VEG)
      return_value = 'veg'
    CASE (RAD)
      return_value = 'rad'
#ifndef __NO_JSBACH_HD__
    CASE (HD)
      return_value = 'hd'
#endif
    CASE DEFAULT
      return_value = ''
    END SELECT

  END FUNCTION Get_process_name

  SUBROUTINE Get_atm2land_memory(mem, iproc, return_ptr)

    TYPE(t_jsb_memory_p), INTENT(in) :: mem(:)
    INTEGER,              INTENT(in) :: iproc
    TYPE(t_atm2land_memory),  POINTER    :: return_ptr

    CLASS(*),  POINTER    :: ptr

    ptr => mem(iproc)%p
    IF (.NOT. ASSOCIATED(ptr)) RETURN

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
    IF (.NOT. ASSOCIATED(ptr)) RETURN

    SELECT TYPE(ptr)
    TYPE IS (t_land2atm_memory)
      return_ptr => ptr
    END SELECT
    NULLIFY(ptr)

  END SUBROUTINE Get_land2atm_memory

  SUBROUTINE Get_seb_memory(mem, iproc, return_ptr)

    TYPE(t_jsb_memory_p), INTENT(in) :: mem(:)
    INTEGER,              INTENT(in) :: iproc
    TYPE(t_seb_memory),  POINTER    :: return_ptr

    CLASS(*),  POINTER    :: ptr

    ptr => mem(iproc)%p
    IF (.NOT. ASSOCIATED(ptr)) RETURN

    SELECT TYPE(ptr)
    TYPE IS (t_seb_memory)
      return_ptr => ptr
    END SELECT
    NULLIFY(ptr)

  END SUBROUTINE Get_seb_memory

  SUBROUTINE Get_turb_memory(mem, iproc, return_ptr)

    TYPE(t_jsb_memory_p), INTENT(in) :: mem(:)
    INTEGER,              INTENT(in) :: iproc
    TYPE(t_turb_memory),  POINTER    :: return_ptr

    CLASS(*),  POINTER    :: ptr

    ptr => mem(iproc)%p
    IF (.NOT. ASSOCIATED(ptr)) RETURN

    SELECT TYPE(ptr)
    TYPE IS (t_turb_memory)
      return_ptr => ptr
    END SELECT
    NULLIFY(ptr)

  END SUBROUTINE Get_turb_memory

  SUBROUTINE Get_soil_memory(mem, iproc, return_ptr)

    TYPE(t_jsb_memory_p), INTENT(in) :: mem(:)
    INTEGER,              INTENT(in) :: iproc
    TYPE(t_soil_memory),  POINTER    :: return_ptr

    CLASS(*),  POINTER    :: ptr

    ptr => mem(iproc)%p
    IF (.NOT. ASSOCIATED(ptr)) RETURN

    SELECT TYPE(ptr)
    TYPE IS (t_soil_memory)
      return_ptr => ptr
    END SELECT
    NULLIFY(ptr)

  END SUBROUTINE Get_soil_memory

  SUBROUTINE Get_hydro_memory(mem, iproc, return_ptr)

    TYPE(t_jsb_memory_p), INTENT(in) :: mem(:)
    INTEGER,              INTENT(in) :: iproc
    TYPE(t_hydro_memory), POINTER    :: return_ptr

    CLASS(*),  POINTER    :: ptr

    ptr => mem(iproc)%p
    IF (.NOT. ASSOCIATED(ptr)) RETURN

    SELECT TYPE(ptr)
    TYPE IS (t_hydro_memory)
      return_ptr => ptr
    END SELECT
    NULLIFY(ptr)

  END SUBROUTINE Get_hydro_memory

  SUBROUTINE Get_veg_memory(mem, iproc, return_ptr)

    TYPE(t_jsb_memory_p), INTENT(in) :: mem(:)
    INTEGER,              INTENT(in) :: iproc
    TYPE(t_veg_memory),  POINTER    :: return_ptr

    CLASS(*),  POINTER    :: ptr

    ptr => mem(iproc)%p
    IF (.NOT. ASSOCIATED(ptr)) RETURN

    SELECT TYPE(ptr)
    TYPE IS (t_veg_memory)
      return_ptr => ptr
    END SELECT
    NULLIFY(ptr)

  END SUBROUTINE Get_veg_memory

  SUBROUTINE Get_rad_memory(mem, iproc, return_ptr)

    TYPE(t_jsb_memory_p), INTENT(in) :: mem(:)
    INTEGER,              INTENT(in) :: iproc
    TYPE(t_rad_memory),  POINTER    :: return_ptr

    CLASS(*),  POINTER    :: ptr

    ptr => mem(iproc)%p
    IF (.NOT. ASSOCIATED(ptr)) RETURN

    SELECT TYPE(ptr)
    TYPE IS (t_rad_memory)
      return_ptr => ptr
    END SELECT
    NULLIFY(ptr)

  END SUBROUTINE Get_rad_memory

#ifndef __NO_JSBACH_HD__
  SUBROUTINE Get_hd_memory(mem, iproc, return_ptr)

    TYPE(t_jsb_memory_p), INTENT(in) :: mem(:)
    INTEGER,              INTENT(in) :: iproc
    TYPE(t_hd_memory),  POINTER    :: return_ptr

    CLASS(*),  POINTER    :: ptr

    ptr => mem(iproc)%p
    IF (.NOT. ASSOCIATED(ptr)) RETURN

    SELECT TYPE(ptr)
    TYPE IS (t_hd_memory)
      return_ptr => ptr
    END SELECT
    NULLIFY(ptr)

  END SUBROUTINE Get_hd_memory
#endif

  SUBROUTINE Get_seb_config(config, return_ptr)

    CLASS(t_jsb_config), POINTER :: config
    TYPE(t_seb_config),  POINTER :: return_ptr

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':Get_seb_config'

    return_ptr => NULL()

    IF (.NOT. ASSOCIATED(config)) RETURN

    SELECT TYPE (config)
    TYPE IS (t_seb_config)
      return_ptr => config
    END SELECT

    IF (.NOT. ASSOCIATED(return_ptr)) &
      CALL finish(TRIM(routine), 'Process has no configuration')

  END SUBROUTINE Get_seb_config

  SUBROUTINE Get_turb_config(config, return_ptr)

    CLASS(t_jsb_config), POINTER :: config
    TYPE(t_turb_config), POINTER :: return_ptr

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':Get_turb_config'

    return_ptr => NULL()

    IF (.NOT. ASSOCIATED(config)) RETURN

    SELECT TYPE (config)
    TYPE IS (t_turb_config)
      return_ptr => config
    END SELECT

    IF (.NOT. ASSOCIATED(return_ptr)) &
      CALL finish(TRIM(routine), 'Process has no configuration')

  END SUBROUTINE Get_turb_config

  SUBROUTINE Get_soil_config(config, return_ptr)

    CLASS(t_jsb_config), POINTER :: config
    TYPE(t_soil_config), POINTER :: return_ptr

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':Get_soil_config'

    return_ptr => NULL()

    IF (.NOT. ASSOCIATED(config)) RETURN

    SELECT TYPE (config)
    TYPE IS (t_soil_config)
      return_ptr => config
    END SELECT

    IF (.NOT. ASSOCIATED(return_ptr)) &
      CALL finish(TRIM(routine), 'Process has no configuration')

  END SUBROUTINE Get_soil_config

  SUBROUTINE Get_hydro_config(config, return_ptr)

    CLASS(t_jsb_config), POINTER :: config
    TYPE(t_hydro_config), POINTER :: return_ptr

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':Get_hydro_config'

    return_ptr => NULL()

    IF (.NOT. ASSOCIATED(config)) RETURN

    SELECT TYPE (config)
    TYPE IS (t_hydro_config)
      return_ptr => config
    END SELECT

    IF (.NOT. ASSOCIATED(return_ptr)) &
      CALL finish(TRIM(routine), 'Process has no configuration')

  END SUBROUTINE Get_hydro_config

  SUBROUTINE Get_veg_config(config, return_ptr)

    CLASS(t_jsb_config), POINTER :: config
    TYPE(t_veg_config),  POINTER :: return_ptr

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':Get_veg_config'

    return_ptr => NULL()

    IF (.NOT. ASSOCIATED(config)) RETURN

    SELECT TYPE (config)
    TYPE IS (t_veg_config)
      return_ptr => config
    END SELECT

    IF (.NOT. ASSOCIATED(return_ptr)) &
      CALL finish(TRIM(routine), 'Process has no configuration')

  END SUBROUTINE Get_veg_config

  SUBROUTINE Get_rad_config(config, return_ptr)

    CLASS(t_jsb_config), POINTER :: config
    TYPE(t_rad_config),  POINTER :: return_ptr

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':Get_rad_config'

    return_ptr => NULL()

    IF (.NOT. ASSOCIATED(config)) RETURN

    SELECT TYPE (config)
    TYPE IS (t_rad_config)
      return_ptr => config
    END SELECT

    IF (.NOT. ASSOCIATED(return_ptr)) &
      CALL finish(TRIM(routine), 'Process has no configuration')

  END SUBROUTINE Get_rad_config

#ifndef __NO_JSBACH_HD__
  SUBROUTINE Get_hd_config(config, return_ptr)

    CLASS(t_jsb_config), POINTER :: config
    TYPE(t_hd_config),   POINTER :: return_ptr

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':Get_hd_config'

    return_ptr => NULL()

    IF (.NOT. ASSOCIATED(config)) RETURN

    SELECT TYPE (config)
    TYPE IS (t_hd_config)
      return_ptr => config
    END SELECT

    IF (.NOT. ASSOCIATED(return_ptr)) &
      CALL finish(TRIM(routine), 'Process has no configuration')

  END SUBROUTINE Get_hd_config
#endif

#endif
END MODULE mo_jsb_process_factory
