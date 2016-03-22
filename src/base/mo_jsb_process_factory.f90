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
  USE mo_srf_memory_class,     ONLY: t_srf_memory
  USE mo_soil_memory_class,    ONLY: t_soil_memory
  USE mo_veg_memory_class,     ONLY: t_veg_memory
  USE mo_rad_memory_class,     ONLY: t_rad_memory
#ifndef __NO_JSBACH_HD__
  USE mo_hd_memory_class,      ONLY: t_hd_memory
#endif
  USE mo_glacier_memory_class, ONLY: t_glacier_memory
  USE mo_lake_memory_class,    ONLY: t_lake_memory
  USE mo_srf_config_class,     ONLY: t_srf_config
  USE mo_soil_config_class,    ONLY: t_soil_config
  USE mo_veg_config_class,     ONLY: t_veg_config
  USE mo_rad_config_class,     ONLY: t_rad_config
#ifndef __NO_JSBACH_HD__
  USE mo_hd_config_class,      ONLY: t_hd_config
#endif
  USE mo_glacier_config_class, ONLY: t_glacier_config
  USE mo_lake_config_class,    ONLY: t_lake_config

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: max_no_of_processes, & !Create_process, & !Register_process_tasks, &
            Create_process_memory, Get_process_name, Get_process_memory, Get_process_config
  PUBLIC :: SRF, ATM2LAND, LAND2ATM, SKIN, SOIL, VEG, SNOW, GLACIER, LAKE, HD, RAD

  INTEGER, PARAMETER :: max_no_of_processes = 20

  ENUM, BIND(C)
    ENUMERATOR :: SRF=1, ATM2LAND, LAND2ATM, SKIN, SOIL, VEG, SNOW, GLACIER, LAKE, HD, RAD
  END ENUM

  INTERFACE Get_process_memory
    MODULE PROCEDURE Get_atm2land_memory
    MODULE PROCEDURE Get_land2atm_memory
    MODULE PROCEDURE Get_srf_memory
    MODULE PROCEDURE Get_soil_memory
    MODULE PROCEDURE Get_veg_memory
    MODULE PROCEDURE Get_rad_memory
#ifndef __NO_JSBACH_HD__
    MODULE PROCEDURE Get_hd_memory
#endif
    MODULE PROCEDURE Get_glacier_memory
    MODULE PROCEDURE Get_lake_memory
  END INTERFACE Get_process_memory

  INTERFACE Get_process_config
    MODULE PROCEDURE Get_srf_config
    MODULE PROCEDURE Get_soil_config
    MODULE PROCEDURE Get_veg_config
    MODULE PROCEDURE Get_rad_config
#ifndef __NO_JSBACH_HD__
    MODULE PROCEDURE Get_hd_config
#endif
    MODULE PROCEDURE Get_glacier_config
    MODULE PROCEDURE Get_lake_config
  END INTERFACE Get_process_config

  CHARACTER(len=*), PARAMETER :: modname = 'mo_jsb_process_factory'

CONTAINS

!!$  FUNCTION Create_process(iproc, model_id) RESULT(return_ptr)
!!$
!!$    USE mo_jsb_process_class, ONLY: t_jsb_process
!!$
!!$    INTEGER, INTENT(in) :: iproc
!!$    INTEGER, INTENT(in) :: model_id
!!$    CLASS(t_jsb_process), POINTER :: return_ptr
!!$
!!$    SELECT CASE (iproc)
!!$    CASE (RAD)
!!$      ALLOCATE(t_jsb_process::return_ptr)
!!$      ALLOCATE(t_rad_config::return_ptr%config)
!!$      !CALL return_ptr%Register_task(tsk_Radiation_surface(model_id))
!!$
!!$    CASE DEFAULT
!!$      return_ptr => NULL()
!!$    END SELECT
!!$
!!$    IF (ASSOCIATED(return_ptr)) THEN
!!$      return_ptr%id             = iproc
!!$      return_ptr%owner_model_id = model_id
!!$
!!$      !CALL return_ptr%Register_tasks(model_id)
!!$    END IF
!!$
!!$  END FUNCTION Create_process

!!$  FUNCTION Create_process_config(iproc, namelist_filename) RESULT(return_ptr)
!!$
!!$    USE mo_jsb_config_class,   ONLY: t_jsb_config
!!$
!!$    INTEGER, INTENT(in) :: iproc
!!$    CHARACTER(LEN=*), INTENT(in) :: namelist_filename
!!$    CLASS(t_jsb_config), POINTER :: return_ptr
!!$
!!$    SELECT CASE(iproc)
 !!$    CASE (SRF)
!!$      ALLOCATE(t_srf_config::return_ptr)
!!$    CASE (SOIL)
!!$      ALLOCATE(t_soil_config::return_ptr)
!!$    CASE (VEG)
!!$      ALLOCATE(t_veg_config::return_ptr)
!!$    CASE (RAD)
!!$      ALLOCATE(t_rad_config::return_ptr)
!!$#ifndef __NO_JSBACH_HD__
!!$    CASE (HD)
!!$      ALLOCATE(t_hd_config::return_ptr)
!!$#endif
!!$    CASE (GLACIER)
!!$      ALLOCATE(t_glacier_config::return_ptr)
!!$    CASE (LAKE)
!!$      ALLOCATE(t_lake_config::return_ptr)
!!$    CASE DEFAULT
!!$      return_ptr => NULL()
!!$    END SELECT
!!$
!!$    IF (ASSOCIATED(return_ptr)) return_ptr%namelist_filename = namelist_filename
!!$
!!$  END FUNCTION Create_process_config

!!$  SUBROUTINE Register_process_tasks(iproc, model_id)
!!$
!!$    INTEGER, INTENT(in) :: iproc
!!$    INTEGER, INTENT(in) :: model_id
!!$
!!$    SELECT CASE (iproc)
!!$    CASE (RAD)
!!$      CALL Register_radiation_tasks(model_id)
!!$    END SELECT
!!$
!!$  END SUBROUTINE Register_process_tasks

  FUNCTION Create_process_memory(iproc) RESULT(return_ptr)

    USE mo_atmland_memory_class, ONLY: max_no_of_atmland_vars => max_no_of_vars
    USE mo_srf_memory_class,     ONLY: max_no_of_srf_vars => max_no_of_vars
    USE mo_soil_memory_class,    ONLY: max_no_of_soil_vars => max_no_of_vars
    USE mo_veg_memory_class,     ONLY: max_no_of_veg_vars => max_no_of_vars
    USE mo_rad_memory_class,     ONLY: max_no_of_rad_vars => max_no_of_vars
    USE mo_glacier_memory_class, ONLY: max_no_of_glacier_vars => max_no_of_vars
    USE mo_lake_memory_class,    ONLY: max_no_of_lake_vars => max_no_of_vars
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
    CASE (SRF)
      ALLOCATE(t_srf_memory::return_ptr)
      return_ptr%max_no_of_vars = max_no_of_srf_vars
    CASE (SOIL)
      ALLOCATE(t_soil_memory::return_ptr)
      return_ptr%max_no_of_vars = max_no_of_soil_vars
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
    CASE (GLACIER)
      ALLOCATE(t_glacier_memory::return_ptr)
      return_ptr%max_no_of_vars = max_no_of_glacier_vars
    CASE (LAKE)
      ALLOCATE(t_lake_memory::return_ptr)
      return_ptr%max_no_of_vars = max_no_of_lake_vars
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
    CASE (SRF)
      return_value = 'srf'
    CASE (SOIL)
      return_value = 'soil'
    CASE (VEG)
      return_value = 'veg'
    CASE (RAD)
      return_value = 'rad'
#ifndef __NO_JSBACH_HD__
    CASE (HD)
      return_value = 'hd'
#endif
    CASE (GLACIER)
      return_value = 'glacier'
    CASE (LAKE)
      return_value = 'lake'
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

  SUBROUTINE Get_srf_memory(mem, iproc, return_ptr)

    TYPE(t_jsb_memory_p), INTENT(in) :: mem(:)
    INTEGER,              INTENT(in) :: iproc
    TYPE(t_srf_memory),  POINTER    :: return_ptr

    CLASS(*),  POINTER    :: ptr

    ptr => mem(iproc)%p
    IF (.NOT. ASSOCIATED(ptr)) RETURN

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
    IF (.NOT. ASSOCIATED(ptr)) RETURN

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

  SUBROUTINE Get_glacier_memory(mem, iproc, return_ptr)

    TYPE(t_jsb_memory_p),   INTENT(in) :: mem(:)
    INTEGER,                INTENT(in) :: iproc
    TYPE(t_glacier_memory), POINTER    :: return_ptr

    CLASS(*),  POINTER    :: ptr

    ptr => mem(iproc)%p
    IF (.NOT. ASSOCIATED(ptr)) RETURN

    SELECT TYPE(ptr)
    TYPE IS (t_glacier_memory)
      return_ptr => ptr
    END SELECT
    NULLIFY(ptr)

  END SUBROUTINE Get_glacier_memory

  SUBROUTINE Get_lake_memory(mem, iproc, return_ptr)

    TYPE(t_jsb_memory_p),   INTENT(in) :: mem(:)
    INTEGER,                INTENT(in) :: iproc
    TYPE(t_lake_memory),    POINTER    :: return_ptr

    CLASS(*),  POINTER    :: ptr

    ptr => mem(iproc)%p
    IF (.NOT. ASSOCIATED(ptr)) RETURN

    SELECT TYPE(ptr)
    TYPE IS (t_lake_memory)
      return_ptr => ptr
    END SELECT
    NULLIFY(ptr)

  END SUBROUTINE Get_lake_memory

  SUBROUTINE Get_srf_config(config, return_ptr)

    CLASS(t_jsb_config), POINTER :: config
    TYPE(t_srf_config),  POINTER :: return_ptr

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':Get_srf_config'

    return_ptr => NULL()

    IF (.NOT. ASSOCIATED(config)) RETURN

    SELECT TYPE (config)
    TYPE IS (t_srf_config)
      return_ptr => config
    END SELECT

    IF (.NOT. ASSOCIATED(return_ptr)) &
      CALL finish(TRIM(routine), 'Process has no configuration')

  END SUBROUTINE Get_srf_config

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

  SUBROUTINE Get_glacier_config(config, return_ptr)

    CLASS(t_jsb_config), TARGET, INTENT(in) :: config
    TYPE(t_glacier_config), POINTER         :: return_ptr

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':Get_glacier_config'

    return_ptr => NULL()

    SELECT TYPE (config)
    TYPE IS (t_glacier_config)
      return_ptr => config
    END SELECT

    IF (.NOT. ASSOCIATED(return_ptr)) &
      CALL finish(TRIM(routine), 'Process has no configuration')

  END SUBROUTINE Get_glacier_config

  SUBROUTINE Get_lake_config(config, return_ptr)

    CLASS(t_jsb_config), TARGET, INTENT(in) :: config
    TYPE(t_lake_config), POINTER            :: return_ptr

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':Get_lake_config'

    return_ptr => NULL()

    SELECT TYPE (config)
    TYPE IS (t_lake_config)
      return_ptr => config
    END SELECT

    IF (.NOT. ASSOCIATED(return_ptr)) &
      CALL finish(TRIM(routine), 'Process has no configuration')

  END SUBROUTINE Get_lake_config

#endif
END MODULE mo_jsb_process_factory
