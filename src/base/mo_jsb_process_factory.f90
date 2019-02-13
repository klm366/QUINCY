!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/license .
!!
!! @brief Contains process factory.
!!
!! @author
!!  Reiner Schnur, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Reiner Schnur (2015-06-29)
!! Copied to QUINCY by Jan Engel MPI-BGC (2016-08-17)
!!
MODULE mo_jsb_process_factory
#ifndef __NO_JSBACH__

  USE mo_exception,            ONLY: finish

  ! 
  USE mo_jsb_process_class,    ONLY: t_jsb_process, Get_process_name
  USE mo_jsb_memory_class,     ONLY: t_jsb_memory, t_jsb_memory_p
  USE mo_jsb_config_class,     ONLY: t_jsb_config
  ! 
  USE mo_a2l_memory_class, ONLY: t_a2l_memory
  USE mo_l2a_memory_class, ONLY: t_l2a_memory
!!  USE mo_seb_memory_class,     ONLY: t_seb_memory
!!  USE mo_turb_memory_class,    ONLY: t_turb_memory
  USE mo_sse_memory_class,    ONLY: t_sse_memory
  USE mo_sb_memory_class,      ONLY: t_sb_memory
!!  USE mo_hydro_memory_class,   ONLY: t_hydro_memory
  USE mo_veg_memory_class,     ONLY: t_veg_memory
  USE mo_assimi_memory_class,  ONLY: t_assimi_memory
  USE mo_pheno_memory_class,   ONLY: t_pheno_memory
  USE mo_rad_memory_class,     ONLY: t_rad_memory
!!#ifndef __NO_JSBACH_HD__
!!  USE mo_hd_memory_class,      ONLY: t_hd_memory
!!#endif
  USE mo_rad_interface,        ONLY: Register_rad_tasks
  USE mo_rad_config_class,     ONLY: t_rad_config
!!  USE mo_seb_interface,        ONLY: Register_seb_tasks
!!  USE mo_seb_config_class,     ONLY: t_seb_config
!!  USE mo_turb_interface,       ONLY: Register_turb_tasks
!!  USE mo_turb_config_class,    ONLY: t_turb_config
  USE mo_sse_interface,       ONLY: Register_sse_tasks
  USE mo_sse_config_class,    ONLY: t_sse_config
!!  USE mo_hydro_interface,      ONLY: Register_hydro_tasks
!!  USE mo_hydro_config_class,   ONLY: t_hydro_config
  USE mo_veg_interface,        ONLY: Register_veg_tasks
  USE mo_veg_config_class,     ONLY: t_veg_config
  USE mo_assimi_interface,     ONLY: Register_assimi_tasks
  USE mo_assimi_config_class,  ONLY: t_assimi_config
  USE mo_pheno_interface,      ONLY: Register_pheno_tasks
  USE mo_pheno_config_class,   ONLY: t_pheno_config  
  USE mo_sb_interface,         ONLY: Register_sb_tasks
  USE mo_sb_config_class,      ONLY: t_sb_config  
!!#ifndef __NO_JSBACH_HD__
!!  USE mo_hd_interface,         ONLY: Register_hd_tasks
!!  USE mo_hd_config_class,      ONLY: t_hd_config
!!#endif

  USE mo_jsb_process_class, ONLY: A2L_, L2A_, SEB_, TURB_, SSE_, HYDRO_, VEG_, HD_, RAD_, ASSIMI_, PHENO_, SB_

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
   CASE (RAD_)
     CALL Register_rad_tasks(return_ptr, model_id)
     ALLOCATE(t_rad_config::return_ptr%config)
  !  CASE (SEB_)
  !    CALL Register_seb_tasks(return_ptr, model_id)
  !    ALLOCATE(t_seb_config::return_ptr%config)
  !  CASE (TURB_)
  !    CALL Register_turb_tasks(return_ptr, model_id)
  !    ALLOCATE(t_turb_config::return_ptr%config)
   CASE (SSE_)
     CALL Register_sse_tasks(return_ptr, model_id)
     ALLOCATE(t_sse_config::return_ptr%config)
  !  CASE (HYDRO_)
  !    CALL Register_hydro_tasks(return_ptr, model_id)
  !    ALLOCATE(t_hydro_config::return_ptr%config)
   CASE (VEG_)
     CALL Register_veg_tasks(return_ptr, model_id)
     ALLOCATE(t_veg_config::return_ptr%config)
   CASE (ASSIMI_)
     CALL Register_assimi_tasks(return_ptr, model_id)
     ALLOCATE(t_assimi_config::return_ptr%config)
   CASE (PHENO_)
     CALL Register_pheno_tasks(return_ptr, model_id)
     ALLOCATE(t_pheno_config::return_ptr%config)
   CASE (SB_)
     CALL Register_sb_tasks(return_ptr, model_id)
     ALLOCATE(t_sb_config::return_ptr%config)
! #ifndef __NO_JSBACH_HD__
!    CASE (HD_)
!      CALL Register_hd_tasks(return_ptr, model_id)
!      ALLOCATE(t_hd_config::return_ptr%config)
! #endif
   CASE DEFAULT
     return_ptr => NULL()
   END SELECT

   IF (ASSOCIATED(return_ptr)) THEN
     return_ptr%config%namelist_filename = TRIM(namelist_filename)
   END IF

 END FUNCTION Create_process

  FUNCTION Create_process_memory(iproc) RESULT(return_ptr)

    USE mo_a2l_memory_class, ONLY: max_no_of_a2l_vars => max_no_of_vars
    USE mo_l2a_memory_class, ONLY: max_no_of_l2a_vars => max_no_of_vars
!!    USE mo_seb_memory_class,     ONLY: max_no_of_seb_vars => max_no_of_vars
!!    USE mo_turb_memory_class,    ONLY: max_no_of_turb_vars => max_no_of_vars
    USE mo_sse_memory_class,    ONLY: max_no_of_sse_vars      => max_no_of_vars
    USE mo_sb_memory_class,      ONLY: max_no_of_sb_vars        => max_no_of_vars, &
                                       max_no_of_sb_structs     => max_no_of_structs, &    ! *slmjsb4
                                       max_no_of_sb_n_structs   => max_no_of_n_structs     ! *slmjsb4
!!    USE mo_hydro_memory_class,   ONLY: max_no_of_hydro_vars => max_no_of_vars
    USE mo_veg_memory_class,     ONLY: max_no_of_veg_vars      => max_no_of_vars, &
                                       max_no_of_veg_structs   => max_no_of_structs, &    ! *slmjsb4
                                       max_no_of_veg_n_structs => max_no_of_n_structs     ! *slmjsb4
    USE mo_assimi_memory_class,  ONLY: max_no_of_assimi_vars => max_no_of_vars
    USE mo_pheno_memory_class,   ONLY: max_no_of_pheno_vars => max_no_of_vars
    USE mo_rad_memory_class,     ONLY: max_no_of_rad_vars => max_no_of_vars
!!#ifndef __NO_JSBACH_HD__
!!    USE mo_hd_memory_class,      ONLY: max_no_of_hd_vars => max_no_of_vars
!!#endif

    IMPLICIT NONE
    ! 0.1 InOut
    INTEGER, INTENT(in) :: iproc
    CLASS(t_jsb_memory), POINTER   :: return_ptr
    ! 0.2 local
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':Create_process_memory'


    SELECT CASE(iproc)
    CASE (A2L_)
      ALLOCATE(t_a2l_memory::return_ptr)
      return_ptr%max_no_of_vars = max_no_of_a2l_vars
    CASE (L2A_)
      ALLOCATE(t_l2a_memory::return_ptr)
      return_ptr%max_no_of_vars = max_no_of_l2a_vars
!!    CASE (SEB_)
!!      ALLOCATE(t_seb_memory::return_ptr)
!!      return_ptr%max_no_of_vars = max_no_of_seb_vars
!!    CASE (TURB_)
!!      ALLOCATE(t_turb_memory::return_ptr)
!!      return_ptr%max_no_of_vars = max_no_of_turb_vars
    CASE (SSE_)
      ALLOCATE(t_sse_memory::return_ptr)
      return_ptr%max_no_of_vars      = max_no_of_sse_vars
    CASE (SB_)
      ALLOCATE(t_sb_memory::return_ptr)
      return_ptr%max_no_of_vars      = max_no_of_sb_vars
      return_ptr%max_no_of_structs   = max_no_of_sb_structs        ! *slmjsb4
      return_ptr%max_no_of_n_structs = max_no_of_sb_n_structs      ! *slmjsb4
!!    CASE (HYDRO_)
!!      ALLOCATE(t_hydro_memory::return_ptr)
!!      return_ptr%max_no_of_vars = max_no_of_hydro_vars
    CASE (VEG_)
      ALLOCATE(t_veg_memory::return_ptr)
      return_ptr%max_no_of_vars      = max_no_of_veg_vars
      return_ptr%max_no_of_structs   = max_no_of_veg_structs        ! *slmjsb4
      return_ptr%max_no_of_n_structs = max_no_of_veg_n_structs      ! *slmjsb4
    CASE (ASSIMI_)
      ALLOCATE(t_assimi_memory::return_ptr)
      return_ptr%max_no_of_vars = max_no_of_assimi_vars
    CASE (PHENO_)
      ALLOCATE(t_pheno_memory::return_ptr)
      return_ptr%max_no_of_vars = max_no_of_pheno_vars
    CASE (RAD_)
      ALLOCATE(t_rad_memory::return_ptr)
      return_ptr%max_no_of_vars = max_no_of_rad_vars
!!#ifndef __NO_JSBACH_HD__
!!    CASE (HD_)
!!      ALLOCATE(t_hd_memory::return_ptr)
!!      return_ptr%max_no_of_vars = max_no_of_hd_vars
!!#endif
    END SELECT

    ALLOCATE(return_ptr%vars(return_ptr%max_no_of_vars))
    ! allocate TYPEs and nested TYPEs only for certain processes   ! *slmjsb4 (VEG_ & SB_)
    SELECT CASE(iproc)
    CASE (VEG_, SB_)
      ALLOCATE(return_ptr%structs(return_ptr%max_no_of_structs))      ! *slmjsb4
      ALLOCATE(return_ptr%n_structs(return_ptr%max_no_of_n_structs))  ! *slmjsb4
    END SELECT

  END FUNCTION Create_process_memory

#endif
END MODULE mo_jsb_process_factory
