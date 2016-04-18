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

  !
  USE mo_jsb_process_class,    ONLY: t_jsb_process, Get_process_name
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
  USE mo_radiation_interface,  ONLY: Register_radiation_tasks
  USE mo_rad_config_class,     ONLY: t_rad_config
  USE mo_seb_interface,        ONLY: Register_seb_tasks
  USE mo_seb_config_class,     ONLY: t_seb_config
  USE mo_turb_interface,       ONLY: Register_turb_tasks
  USE mo_turb_config_class,    ONLY: t_turb_config
  USE mo_soil_interface,       ONLY: Register_soil_tasks
  USE mo_soil_config_class,    ONLY: t_soil_config
  USE mo_hydro_interface,      ONLY: Register_hydro_tasks
  USE mo_hydro_config_class,   ONLY: t_hydro_config
  USE mo_veg_interface,        ONLY: Register_veg_tasks
  USE mo_veg_config_class,     ONLY: t_veg_config
#ifndef __NO_JSBACH_HD__
  USE mo_hd_interface,         ONLY: Register_hd_tasks
  USE mo_hd_config_class,      ONLY: t_hd_config
#endif

  USE mo_jsb_process_class, ONLY: ATM2LAND, LAND2ATM, SEB, TURB, SOIL, HYDRO, VEG, HD, RAD

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: Create_process_memory
  PUBLIC :: Create_process, max_no_of_processes

  INTEGER, PARAMETER :: max_no_of_processes = 20

  CHARACTER(len=*), PARAMETER :: modname = 'mo_jsb_process_factory'

CONTAINS

  FUNCTION Create_process(iproc, model_id, namelist_filename) RESULT(return_ptr)

    INTEGER, INTENT(in) :: iproc
    INTEGER, INTENT(in) :: model_id
    CHARACTER(LEN=*), INTENT(in) :: namelist_filename
    CLASS(t_jsb_process), POINTER :: return_ptr

    ALLOCATE(t_jsb_process::return_ptr)

    return_ptr%id             = iproc
    return_ptr%name           = Get_process_name(iproc)
    return_ptr%owner_model_id = model_id

    SELECT CASE (iproc)
    CASE (RAD)
      CALL Register_radiation_tasks(return_ptr, model_id)
      ALLOCATE(t_rad_config::return_ptr%config)
    CASE (SEB)
      CALL Register_seb_tasks(return_ptr, model_id)
      ALLOCATE(t_seb_config::return_ptr%config)
    CASE (TURB)
      CALL Register_turb_tasks(return_ptr, model_id)
      ALLOCATE(t_turb_config::return_ptr%config)
    CASE (SOIL)
      CALL Register_soil_tasks(return_ptr, model_id)
      ALLOCATE(t_soil_config::return_ptr%config)
    CASE (HYDRO)
      CALL Register_hydro_tasks(return_ptr, model_id)
      ALLOCATE(t_hydro_config::return_ptr%config)
    CASE (VEG)
      CALL Register_veg_tasks(return_ptr, model_id)
      ALLOCATE(t_veg_config::return_ptr%config)
#ifndef __NO_JSBACH_HD__
    CASE (HD)
      CALL Register_hd_tasks(return_ptr, model_id)
      ALLOCATE(t_hd_config::return_ptr%config)
#endif
    CASE DEFAULT
      return_ptr => NULL()
    END SELECT

    IF (ASSOCIATED(return_ptr)) THEN
      return_ptr%config%namelist_filename = TRIM(namelist_filename)
    END IF

  END FUNCTION Create_process

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

#endif
END MODULE mo_jsb_process_factory
