!>
!! @brief Contains process class
!!
!! @author
!!  Reiner Schnur, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Reiner Schnur (2015-09-01)
!!
MODULE mo_jsb_process
#ifndef __NO_JSBACH__

  USE mo_exception, ONLY: finish

!!$  USE mo_jsb_tile_class, ONLY: t_jsb_tile
!!$  USE mo_jsb_model_class, ONLY: t_jsb_model
  USE mo_jsb_process_factory
  USE mo_jsb_process_class,    ONLY: t_jsb_process
  USE mo_radiation_interface,  ONLY: t_Radiation_process
  USE mo_rad_config_class,     ONLY: t_rad_config
  USE mo_srf_interface,        ONLY: t_srf_process
  USE mo_srf_config_class,     ONLY: t_srf_config
  USE mo_soil_interface,       ONLY: t_soil_process
  USE mo_soil_config_class,    ONLY: t_soil_config
  USE mo_veg_interface,        ONLY: t_veg_process
  USE mo_veg_config_class,     ONLY: t_veg_config
#ifndef __NO_JSBACH_HD__
  USE mo_hd_interface,         ONLY: t_hd_process
  USE mo_hd_config_class,      ONLY: t_hd_config
#endif
  USE mo_glacier_interface,    ONLY: t_glacier_process
  USE mo_glacier_config_class, ONLY: t_glacier_config
  USE mo_lake_interface,       ONLY: t_lake_process
  USE mo_lake_config_class,    ONLY: t_lake_config

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: Create_process

CONTAINS

  FUNCTION Create_process(iproc, model_id, namelist_filename) RESULT(return_ptr)

    USE mo_jsb_process_class, ONLY: t_jsb_process

    INTEGER, INTENT(in) :: iproc
    INTEGER, INTENT(in) :: model_id
    CHARACTER(LEN=*), INTENT(in) :: namelist_filename
    CLASS(t_jsb_process), POINTER :: return_ptr

    SELECT CASE (iproc)
    CASE (RAD)
      ALLOCATE(t_radiation_process::return_ptr)
      ALLOCATE(t_rad_config::return_ptr%config)
    CASE (SRF)
      ALLOCATE(t_srf_process::return_ptr)
      ALLOCATE(t_srf_config::return_ptr%config)
    CASE (SOIL)
      ALLOCATE(t_soil_process::return_ptr)
      ALLOCATE(t_soil_config::return_ptr%config)
    CASE (VEG)
      ALLOCATE(t_veg_process::return_ptr)
      ALLOCATE(t_veg_config::return_ptr%config)
#ifndef __NO_JSBACH_HD__
    CASE (HD)
      ALLOCATE(t_hd_process::return_ptr)
      ALLOCATE(t_hd_config::return_ptr%config)
#endif
    CASE (GLACIER)
      ALLOCATE(t_glacier_process::return_ptr)
      ALLOCATE(t_glacier_config::return_ptr%config)
    CASE (LAKE)
      ALLOCATE(t_lake_process::return_ptr)
      ALLOCATE(t_lake_config::return_ptr%config)
    CASE DEFAULT
      return_ptr => NULL()
    END SELECT

    IF (ASSOCIATED(return_ptr)) THEN
      return_ptr%id             = iproc
      return_ptr%owner_model_id = model_id
      return_ptr%config%namelist_filename = TRIM(namelist_filename)
      !ALLOCATE(return_ptr%task_list)

      CALL return_ptr%Register_tasks(model_id)
    END IF

  END FUNCTION Create_process
!!$  FUNCTION Create_process(iproc, model_id) RESULT(return_ptr)
!!$
!!$    INTEGER, INTENT(in) :: iproc
!!$    INTEGER, INTENT(in) :: model_id
!!$    CLASS(t_jsb_process), POINTER :: return_ptr
!!$
!!$    SELECT CASE (iproc)
!!$    CASE(RAD)
!!$      ALLOCATE(t_Radiation_process::return_ptr)
!!$    CASE DEFAULT
!!$      return_ptr => NULL()
!!$    END SELECT
!!$
!!$    IF (ASSOCIATED(return_ptr)) THEN
!!$      return_ptr%id             = iproc
!!$      return_ptr%owner_model_id = model_id
!!$
!!$      CALL return_ptr%Register_tasks(model_id)
!!$    END IF
!!$
!!$  END FUNCTION Create_process

#endif
END MODULE mo_jsb_process
