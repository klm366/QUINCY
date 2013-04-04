!>
!! @brief Contains interface to parallel scatter/gather routines
!!
!! @author
!!  Reiner Schnur, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Reiner Schnur (2012-10-09)
!!
MODULE mo_jsb_parallel_iface

#ifdef __ICON__
  USE mo_gather_scatter,  ONLY: gather_field => gather_cells, scatter_field => scatter_cells
  USE mo_mpi,             ONLY: p_comm_work_test, p_comm_work, my_process_is_stdio, my_process_is_mpi_parallel, &
                                get_my_global_mpi_id
  USE mo_parallel_config, ONLY: p_test_run
  USE mo_model_domain,    ONLY: t_patch, p_patch
#else
  USE mo_decomposition,   ONLY: local_decomposition, global_decomposition
  USE mo_tr_scatter,      ONLY: scatter_field
  USE mo_tr_gather,       ONLY: gather_field
  USE mo_mpi,             ONLY: p_all_comm
#endif
  USE mo_mpi, ONLY: p_io,    & !< Processor performing I/O
                    p_bcast    !< Broadcasting routine

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: init_parallel_iface, init_mpi_communicators, my_process_is_stdio, get_my_global_mpi_id
  PUBLIC :: gather_field, scatter_field
  PUBLIC :: p_io, mpi_comm, p_bcast

  INTEGER :: mpi_comm

CONTAINS

  SUBROUTINE init_parallel_iface

#ifdef __ICON__
    IF (p_test_run) THEN
      mpi_comm = p_comm_work_test
    ELSE
      mpi_comm = p_comm_work
    END IF
#else
    mpi_comm = p_all_comm
#endif

  END SUBROUTINE init_parallel_iface

!>
!! Only called for jsbach offline
  SUBROUTINE init_mpi_communicators

#ifdef __ICON__
    USE mo_mpi, ONLY: set_mpi_work_communicators
    USE mo_parallel_config, ONLY: p_test_run, l_test_openmp, num_io_procs
#else
    USE mo_mpi, ONLY: p_init_communicators
    USE mo_control, ONLY: nproca, nprocb, nprocio
#endif

#ifdef __ICON__
    CALL set_mpi_work_communicators(p_test_run, l_test_openmp, num_io_procs)
#else
    CALL p_init_communicators(nproca, nprocb, nprocio)
#endif

  END SUBROUTINE init_mpi_communicators

#ifndef __ICON__
  LOGICAL FUNCTION my_process_is_stdio()

  USE mo_mpi, ONLY: p_parallel_io

#ifdef NOMPI
    my_process_is_stdio = .TRUE.
#else
    my_process_is_stdio = p_parallel_io
#endif

  END FUNCTION my_process_is_stdio

  LOGICAL FUNCTION my_process_is_mpi_parallel()

    USE mo_mpi, ONLY: p_parallel

    my_process_is_mpi_parallel = p_parallel

  END FUNCTION my_process_is_mpi_parallel

  INTEGER FUNCTION get_my_global_mpi_id()

    USE mo_mpi, ONLY: p_pe

    get_my_global_mpi_id = p_pe

  END FUNCTION get_my_global_mpi_id
#endif

END MODULE mo_jsb_parallel_iface
