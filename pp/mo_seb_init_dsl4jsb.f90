!
! This file was automatically generated from:
!    /scratch/snx3000/rschnur/icon-rc-testing-plus-mc8.test/externals/jsbach/src/srf_energy_bal/mo_seb_init.f90
!    with md5 hash: 4e55a0fbf37ae06cabd1cd25f6f24c64
!
!>
!! @brief Initialization of the the seb memory.
!!
!! @author
!!  Reiner Schnur, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Reiner Schnur (2013-04-04)
!!
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/license .
!!

MODULE mo_seb_init
#ifndef __NO_JSBACH__

  USE mo_kind,              ONLY: wp
  USE mo_exception,         ONLY: message
  USE mo_jsb_control,       ONLY: debug_on

  !USE mo_jsb_model_class,   ONLY: t_jsb_model
  USE mo_jsb_tile_class,    ONLY: t_jsb_tile_abstract
  !USE mo_jsb_class,         ONLY: get_model
  !USE mo_jsb_io_netcdf,     ONLY: t_input_file !, jsb_netcdf_open_input

  USE mo_jsb_process_class, ONLY: seb_
  USE mo_seb_config_class, ONLY: t_seb_config
  USE mo_seb_memory_class, ONLY: t_seb_memory

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: seb_init

  CHARACTER(len=*), PARAMETER :: modname = 'mo_seb_init'

CONTAINS

  !
  !> Intialize soil process (after memory has been set up)
  !
  SUBROUTINE seb_init(tile)
    USE mo_jsb_time,        ONLY: timeStep_in_days
    USE mo_pheno_constants, ONLY: PhenoParam

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile

    CHARACTER(len=*), PARAMETER :: routine = modname//':seb_init'

    REAL(wp), POINTER, DIMENSION(:,:) :: &
      & F_pseudo_soil_temp,    &
      & N_pseudo_soil_temp

    TYPE(t_seb_memory), POINTER :: seb__mem

    CALL seb_init_bc(tile)
    !CALL seb_init_ic(tile)

    ! Initial values for the calculation of pseudo_soil_temp in SUBROUTINE calc_pseudo_soil_temp
    SELECT TYPE (mm => tile%mem(seb_)%p); TYPE IS (t_seb_memory); seb__mem => mm ; END SELECT


    IF (tile%contains_land) THEN
      n_pseudo_soil_temp => seb__mem%n_pseudo_soil_temp%ptr(:,:) ! OUT
      f_pseudo_soil_temp => seb__mem%f_pseudo_soil_temp%ptr(:,:) ! OUT

      F_pseudo_soil_temp = EXP(-timeStep_in_days()/PhenoParam%EG_SG%tau_pseudo_soil)
      N_pseudo_soil_temp = 1._wp / (1._wp - F_pseudo_soil_temp)

      !$acc update device(F_pseudo_soil_temp, N_pseudo_soil_temp)
    END IF

  END SUBROUTINE seb_init

  SUBROUTINE seb_init_bc(tile)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile

    !TYPE(t_jsb_model),       POINTER :: model

    CHARACTER(len=*), PARAMETER :: routine = modname//':seb_init_bc'

    IF (tile%level > 0) CONTINUE ! only to avoid compiler warnings
    !model => get_model(tile%owner_model_id)

  END SUBROUTINE seb_init_bc

  SUBROUTINE seb_init_ic(tile)

    !USE mo_jsb_time, ONLY: get_time_interpolation_weights

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile

    ! TYPE(t_jsb_model),       POINTER :: model

    ! TYPE(t_seb_config), POINTER :: seb__conf
    !TYPE(t_seb_memory), POINTER :: seb__mem

    !REAL(wp), POINTER :: &
    !  & ptr_2D(:,  :)       !< temporary pointer

    !TYPE(t_input_file) :: input_file

    !INTEGER  :: i

    CHARACTER(len=*), PARAMETER :: routine = modname//':seb_init_ic'

    !model => get_model(tile%owner_model_id)

    ! Get seb config
    ! SELECT TYPE (cc => model%processes(seb_)%p%config); TYPE IS (t_seb_config) ; seb__conf => cc ; END SELECT

    ! IF (.NOT. seb__conf%active) RETURN

    ! Get seb memory of the tile
    !SELECT TYPE (mm => tile%mem(seb_)%p); TYPE IS (t_seb_memory); seb__mem => mm ; END SELECT

    ! Initialize parameters
    !
    ! IF (debug_on()) CALL message(TRIM(routine), 'Setting initial state of seb memory for tile '// &
      ! &                          TRIM(tile%name)//' from '//TRIM(seb__conf%ic_filename))

    ! Initial surface temperature
    ! Is initialized in sse_init_ic together with soil temperatures

  END SUBROUTINE seb_init_ic

#endif
END MODULE mo_seb_init
