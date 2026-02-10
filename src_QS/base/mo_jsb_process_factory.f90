!> Contains process factory
!>
!> ICON-Land
!>
!> ---------------------------------------
!> Copyright (C) 2013-2026, MPI-M, MPI-BGC
!>
!> Contact: icon-model.org
!> Authors: AUTHORS.md
!> See LICENSES/ for license information
!> SPDX-License-Identifier: BSD-3-Clause
!> ---------------------------------------
!>
MODULE mo_jsb_process_factory

  USE mo_exception,            ONLY: finish

  ! 
  USE mo_jsb_process_class,    ONLY: t_jsb_process, Get_process_name
  USE mo_jsb_memory_class,     ONLY: t_jsb_memory, t_jsb_memory_p
  USE mo_jsb_config_class,     ONLY: t_jsb_config
  ! 
  USE mo_a2l_memory_class,     ONLY: t_a2l_memory
  USE mo_l2a_memory_class,     ONLY: t_l2a_memory
  USE mo_seb_memory_class,     ONLY: t_seb_memory
  USE mo_turb_memory_class,    ONLY: t_turb_memory
  USE mo_sse_memory_class,     ONLY: t_sse_memory
  USE mo_sb_memory_class,      ONLY: t_sb_memory
  USE mo_hydro_memory_class,   ONLY: t_hydro_memory
  USE mo_veg_memory_class,     ONLY: t_veg_memory
  USE mo_q_assimi_memory_class,ONLY: t_q_assimi_memory
  USE mo_q_pheno_memory_class, ONLY: t_q_pheno_memory
  USE mo_dfire_memory_class,   ONLY: t_dfire_memory
  USE mo_rad_memory_class,     ONLY: t_rad_memory
  USE mo_phyd_memory_class,    ONLY: t_phyd_memory
  USE mo_q_syl_memory_class,   ONLY: t_q_syl_memory
  USE mo_q_agr_memory_class,   ONLY: t_q_agr_memory
  USE mo_spq_memory_class,     ONLY: t_spq_memory

!!#ifndef __NO_JSBACH_HD__
!!  USE mo_hd_memory_class,      ONLY: t_hd_memory
!!#endif
  USE mo_rad_interface,        ONLY: Register_rad_tasks
  USE mo_rad_config_class,     ONLY: t_rad_config
!!  USE mo_seb_interface,        ONLY: Register_seb_tasks
  USE mo_seb_config_class,     ONLY: t_seb_config
  USE mo_turb_interface,       ONLY: Register_turb_tasks
  USE mo_turb_config_class,    ONLY: t_turb_config
!!  USE mo_sse_interface,        ONLY: Register_sse_tasks_quincy
  USE mo_sse_config_class,     ONLY: t_sse_config
!!  USE mo_hydro_interface,      ONLY: Register_hydro_tasks
  USE mo_hydro_config_class,   ONLY: t_hydro_config
  USE mo_q_veg_interface,      ONLY: Register_veg_tasks_quincy
  USE mo_veg_config_class,     ONLY: t_veg_config
  USE mo_q_assimi_interface,   ONLY: Register_q_assimi_tasks
  USE mo_q_assimi_config_class,ONLY: t_q_assimi_config
  USE mo_q_pheno_interface,    ONLY: Register_q_pheno_tasks
  USE mo_q_pheno_config_class, ONLY: t_q_pheno_config
  USE mo_dfire_interface,      ONLY: Register_dfire_tasks
  USE mo_dfire_config_class,   ONLY: t_dfire_config
  USE mo_q_sb_interface,       ONLY: Register_sb_tasks_quincy
  USE mo_sb_config_class,      ONLY: t_sb_config  
  USE mo_phyd_interface,       ONLY: Register_phyd_tasks
  USE mo_phyd_config_class,    ONLY: t_phyd_config
  USE mo_q_syl_interface,      ONLY: Register_q_syl_tasks
  USE mo_q_syl_config_class,   ONLY: t_q_syl_config
  USE mo_q_agr_interface,      ONLY: Register_q_agr_tasks
  USE mo_q_agr_config_class,   ONLY: t_q_agr_config
  USE mo_q_spq_interface,      ONLY: Register_spq_tasks_quincy
  USE mo_spq_config_class,     ONLY: t_spq_config
!!#ifndef __NO_JSBACH_HD__
!!  USE mo_hd_interface,         ONLY: Register_hd_tasks
!!  USE mo_hd_config_class,      ONLY: t_hd_config
!!#endif

  USE mo_jsb_process_class, ONLY: A2L_, L2A_, SEB_, TURB_, SSE_, HYDRO_, VEG_, HD_, Q_ASSIMI_, Q_PHENO_, &
    & SB_, DFIRE_, PHYD_, RAD_, Q_SYL_, Q_AGR_, SPQ_
  
  IMPLICIT NONE
  PRIVATE

  PUBLIC :: Create_process_memory
  PUBLIC :: Create_process, max_no_of_processes

  INTEGER, PARAMETER :: max_no_of_processes = 42

  CHARACTER(len=*), PARAMETER :: modname = 'mo_jsb_process_factory'

CONTAINS

  FUNCTION Create_process( &
    & iproc, &
    & model_id, &
    & namelist_filename, &
    & model_scheme) RESULT(return_ptr)

    INTEGER,          INTENT(in)  :: iproc
    INTEGER,          INTENT(in)  :: model_id
    CHARACTER(LEN=*), INTENT(in)  :: namelist_filename
    INTEGER,          INTENT(in)  :: model_scheme
    CLASS(t_jsb_process), POINTER :: return_ptr

   ALLOCATE(t_jsb_process::return_ptr)

   return_ptr%id             = iproc
   return_ptr%name           = Get_process_name(iproc)
   return_ptr%owner_model_id = model_id

   SELECT CASE (iproc)
   CASE (RAD_)
     CALL Register_rad_tasks(return_ptr, model_id)
     ALLOCATE(t_rad_config::return_ptr%config)
   CASE (SEB_)
  !   CALL Register_seb_tasks(return_ptr, model_id)
     ALLOCATE(t_seb_config::return_ptr%config)
   CASE (TURB_)
     CALL Register_turb_tasks(return_ptr, model_id)
     ALLOCATE(t_turb_config::return_ptr%config)
   CASE (SSE_)
  !   CALL Register_sse_tasks_quincy(return_ptr, model_id)
     ALLOCATE(t_sse_config::return_ptr%config)
   CASE (HYDRO_)
  !   CALL Register_hydro_tasks(return_ptr, model_id)
     ALLOCATE(t_hydro_config::return_ptr%config)
   CASE (VEG_)
     CALL Register_veg_tasks_quincy(return_ptr, model_id)
     ALLOCATE(t_veg_config::return_ptr%config)
   CASE (Q_ASSIMI_)
     CALL Register_q_assimi_tasks(return_ptr, model_id)
     ALLOCATE(t_q_assimi_config::return_ptr%config)
   CASE (Q_PHENO_)
     CALL Register_q_pheno_tasks(return_ptr, model_id)
     ALLOCATE(t_q_pheno_config::return_ptr%config)
   CASE (DFIRE_)
     CALL Register_dfire_tasks(return_ptr, model_id)
     ALLOCATE(t_dfire_config::return_ptr%config)
   CASE (SB_)
     CALL Register_sb_tasks_quincy(return_ptr, model_id)
     ALLOCATE(t_sb_config::return_ptr%config)
  CASE (PHYD_)
    CALL Register_phyd_tasks(return_ptr, model_id)
    ALLOCATE(t_phyd_config::return_ptr%config)
  CASE (Q_SYL_)
    CALL Register_q_syl_tasks(return_ptr, model_id)
    ALLOCATE(t_q_syl_config::return_ptr%config)
  CASE (Q_AGR_)
    CALL Register_q_agr_tasks(return_ptr, model_id)
    ALLOCATE(t_q_agr_config::return_ptr%config)
  CASE (SPQ_)
    CALL Register_spq_tasks_quincy(return_ptr, model_id)
    ALLOCATE(t_spq_config::return_ptr%config)
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

    USE mo_a2l_memory_class,     ONLY: max_no_of_a2l_vars       => max_no_of_vars
    USE mo_l2a_memory_class,     ONLY: max_no_of_l2a_vars       => max_no_of_vars
    USE mo_seb_memory_class,     ONLY: max_no_of_seb_vars       => max_no_of_vars
    USE mo_turb_memory_class,    ONLY: max_no_of_turb_vars      => max_no_of_vars
    USE mo_sse_memory_class,     ONLY: max_no_of_sse_vars       => max_no_of_vars
    USE mo_sb_memory_class,      ONLY: max_no_of_sb_vars        => max_no_of_vars
    USE mo_hydro_memory_class,   ONLY: max_no_of_hydro_vars     => max_no_of_vars
    USE mo_veg_memory_class,     ONLY: max_no_of_veg_vars       => max_no_of_vars
    USE mo_q_assimi_memory_class,ONLY: max_no_of_q_assimi_vars  => max_no_of_vars
    USE mo_q_pheno_memory_class, ONLY: max_no_of_q_pheno_vars   => max_no_of_vars
    USE mo_phyd_memory_class,    ONLY: max_no_of_phyd_vars      => max_no_of_vars
    USE mo_dfire_memory_class,   ONLY: max_no_of_dfire_vars     => max_no_of_vars
    USE mo_rad_memory_class,     ONLY: max_no_of_rad_vars       => max_no_of_vars
    USE mo_q_syl_memory_class,   ONLY: max_no_of_q_syl_vars     => max_no_of_vars
    USE mo_q_agr_memory_class,   ONLY: max_no_of_q_agr_vars     => max_no_of_vars
    USE mo_spq_memory_class,     ONLY: max_no_of_spq_vars       => max_no_of_vars
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
    CASE (SEB_)
      ALLOCATE(t_seb_memory::return_ptr)
      return_ptr%max_no_of_vars = max_no_of_seb_vars
    CASE (TURB_)
      ALLOCATE(t_turb_memory::return_ptr)
      return_ptr%max_no_of_vars = max_no_of_turb_vars
    CASE (SSE_)
      ALLOCATE(t_sse_memory::return_ptr)
      return_ptr%max_no_of_vars      = max_no_of_sse_vars
    CASE (SB_)
      ALLOCATE(t_sb_memory::return_ptr)
      return_ptr%max_no_of_vars      = max_no_of_sb_vars
    CASE (HYDRO_)
      ALLOCATE(t_hydro_memory::return_ptr)
      return_ptr%max_no_of_vars = max_no_of_hydro_vars
    CASE (VEG_)
      ALLOCATE(t_veg_memory::return_ptr)
      return_ptr%max_no_of_vars      = max_no_of_veg_vars
    CASE (Q_ASSIMI_)
      ALLOCATE(t_q_assimi_memory::return_ptr)
      return_ptr%max_no_of_vars = max_no_of_q_assimi_vars
    CASE (Q_PHENO_)
      ALLOCATE(t_q_pheno_memory::return_ptr)
      return_ptr%max_no_of_vars = max_no_of_q_pheno_vars
    CASE (DFIRE_)
      ALLOCATE(t_dfire_memory::return_ptr)
      return_ptr%max_no_of_vars = max_no_of_dfire_vars
    CASE (RAD_)
      ALLOCATE(t_rad_memory::return_ptr)
      return_ptr%max_no_of_vars = max_no_of_rad_vars
    CASE (PHYD_)
      ALLOCATE(t_phyd_memory::return_ptr)
      return_ptr%max_no_of_vars = max_no_of_phyd_vars
    CASE (Q_SYL_)
      ALLOCATE(t_q_syl_memory::return_ptr)
      return_ptr%max_no_of_vars = max_no_of_q_syl_vars
    CASE (Q_AGR_)
      ALLOCATE(t_q_agr_memory::return_ptr)
      return_ptr%max_no_of_vars = max_no_of_q_agr_vars
    CASE (SPQ_)
      ALLOCATE(t_spq_memory::return_ptr)
      return_ptr%max_no_of_vars = max_no_of_spq_vars
!!#ifndef __NO_JSBACH_HD__
!!    CASE (HD_)
!!      ALLOCATE(t_hd_memory::return_ptr)
!!      return_ptr%max_no_of_vars = max_no_of_hd_vars
!!#endif
    END SELECT

    ALLOCATE(return_ptr%vars(return_ptr%max_no_of_vars))
  END FUNCTION Create_process_memory

END MODULE mo_jsb_process_factory
