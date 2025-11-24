!>
!! Initialization of the topmodel wetland scheme.
!!
!>
!> ICON-Land
!>
!> ---------------------------------------
!> Copyright (C) 2013-2024, MPI-M, MPI-BGC
!>
!> Contact: icon-model.org
!> Authors: AUTHORS.md
!> See LICENSES/ for license information
!> SPDX-License-Identifier: BSD-3-Clause
!> ---------------------------------------
!>
!!
!! @author
!!  Stiig Wilkenskjeld, MPI-M Hamburg
!!
!! @par Revision History
!! Adopted from JSBACH3 code by Thomas Kleinen, 2022-12-06
!!
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/license .
!!

MODULE mo_top_init
#ifndef __NO_JSBACH__

  USE mo_kind,               ONLY: wp
  USE mo_exception,          ONLY: message, finish
  USE mo_jsb_control,        ONLY: debug_on

  USE mo_jsb_model_class,    ONLY: t_jsb_model
  USE mo_jsb_grid_class,     ONLY: t_jsb_grid, t_jsb_vgrid
  USE mo_jsb_grid,           ONLY: Get_grid, Get_vgrid
  USE mo_jsb_tile_class,     ONLY: t_jsb_tile_abstract
  USE mo_jsb_class,          ONLY: get_model
  USE mo_jsb_io_netcdf,      ONLY: t_input_file, jsb_netcdf_open_input
  USE mo_jsb_io,             ONLY: missval
  USE mo_jsb_impl_constants, ONLY: ifs_nsoil, ifs_soil_depth
  USE mo_util,               ONLY: i2s => int2string, l2s => logical2string, r2s => real2string
  USE mo_mpi,                ONLY: my_process_is_mpi_workroot
  USE mo_top_util,           ONLY: interpol_reset
  dsl4jsb_Use_processes TOP_, SSE_, SEB_, HYDRO_
  dsl4jsb_Use_config(TOP_)
  dsl4jsb_Use_memory(TOP_)

  IMPLICIT NONE
  PRIVATE
INCLUDE 'netcdf.inc'
  PUBLIC :: topmodel_init

  TYPE t_top_init_vars
    REAL(wp), POINTER ::                  &
      & cti_mean       (:,:)   => NULL(), &
      & cti_mu         (:,:)   => NULL(), &
      & cti_phi        (:,:)   => NULL(), &
      & cti_chi        (:,:)   => NULL(), &
      & cti_basin      (:,:)   => NULL(), &
      & cti_openwater  (:,:)   => NULL(), &
      & frac_inun_presc(:,:,:) => NULL(), &
      & wt_inun_presc  (:,:,:) => NULL()
  END TYPE t_top_init_vars

  TYPE(t_top_init_vars) :: top_init_vars

  CHARACTER(len=*), PARAMETER :: modname = 'mo_topmodel_init'

CONTAINS

  !
  !> Intialize topmodel wetland process (after memory has been set up)
  !
  SUBROUTINE topmodel_init(tile)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile

    TYPE(t_jsb_model), POINTER  :: model
    CHARACTER(len=*), PARAMETER :: routine = modname//':topmodel_init'

    model => Get_model(tile%owner_model_id)

    IF (.NOT. ASSOCIATED(tile%parent_tile)) THEN
      CALL topmodel_read_init_vars(tile)
    END IF

    CALL topmodel_init_bc(tile)
    CALL topmodel_init_ic(tile)

    IF (tile%Is_last_process_tile(TOP_)) THEN
      CALL topmodel_finalize_init_vars()
    END IF

  END SUBROUTINE topmodel_init

  SUBROUTINE topmodel_finalize_init_vars

    IF (ASSOCIATED(top_init_vars%cti_mean)) THEN
      DEALLOCATE(                          &
        & top_init_vars%cti_basin,         &
        & top_init_vars%cti_openwater,     &
        & top_init_vars%cti_mean,          &
        & top_init_vars%cti_mu,            &
        & top_init_vars%cti_phi,           &
        & top_init_vars%cti_chi            &
        & )
    ENDIF

    IF (ASSOCIATED(top_init_vars%frac_inun_presc)) DEALLOCATE(top_init_vars%frac_inun_presc)
    IF (ASSOCIATED(top_init_vars%wt_inun_presc))   DEALLOCATE(top_init_vars%wt_inun_presc)

  END SUBROUTINE topmodel_finalize_init_vars

  SUBROUTINE topmodel_read_init_vars(tile)
    USE mo_jsb_time,     ONLY: get_date_components
    USE mo_time_config,  ONLY: time_config
    USE mo_jsb_parallel, ONLY: my_process_is_stdio, p_bcast, p_io, mpi_comm

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile

    dsl4jsb_Def_config(TOP_)

    REAL(wp), POINTER :: ptr_2D(:,:), ptr_3D(:,:,:) !< temporary pointers

    TYPE(t_input_file) :: input_file

    TYPE(t_jsb_model), POINTER :: model
    TYPE(t_jsb_grid),  POINTER :: grid
    TYPE(t_jsb_vgrid), POINTER :: soil_e
    INTEGER :: nproma, nblks, isoil
    INTEGER :: ntime
    INTEGER :: syear, eyear, smonth, emonth, idst, irec, stat, varid
    INTEGER :: srec(2)
    REAL(wp), POINTER :: rtime(:)
    REAL(wp) :: ttime(2)

    CHARACTER(len=*), PARAMETER :: routine = modname//':topmodel_read_init_vars'

    model => get_model(tile%owner_model_id)
    grid  => get_grid(model%grid_id)
    nproma = grid%Get_nproma()
    nblks  = grid%Get_nblks()

    ! Get topmodel config
    dsl4jsb_Get_config(TOP_)

    IF (.NOT. tile%Is_process_active(TOP_)) RETURN

    IF (debug_on() .OR. dsl4jsb_Config(TOP_)%dbg_process) &
      CALL message(TRIM(routine), 'Reading topmodel boundary conditions for tile '//TRIM(tile%name))

    IF (dsl4jsb_Config(TOP_)%inun_mode == 0) THEN
      ALLOCATE(                                           &
        & top_init_vars%cti_basin    (nproma, nblks    ), &
        & top_init_vars%cti_openwater(nproma, nblks    ), &
        & top_init_vars%cti_mean     (nproma, nblks    ), &
        & top_init_vars%cti_mu       (nproma, nblks    ), &
        & top_init_vars%cti_phi      (nproma, nblks    ), &
        & top_init_vars%cti_chi      (nproma, nblks    )  &
        & )

      top_init_vars%cti_mean      (:,:) = missval
      top_init_vars%cti_mu        (:,:) = missval
      top_init_vars%cti_phi       (:,:) = missval
      top_init_vars%cti_chi       (:,:) = missval
      top_init_vars%cti_basin     (:,:) = missval
      top_init_vars%cti_openwater (:,:) = missval

      input_file = jsb_netcdf_open_input(TRIM(dsl4jsb_Config(TOP_)%ic_filename), model%grid_id)

      ptr_2D => input_file%Read_2d(variable_name='cti_mean'     , fill_array=top_init_vars%cti_mean)
      ptr_2D => input_file%Read_2d(variable_name='cti_mu'       , fill_array=top_init_vars%cti_mu)
      ptr_2D => input_file%Read_2d(variable_name='cti_phi'      , fill_array=top_init_vars%cti_phi)
      ptr_2D => input_file%Read_2d(variable_name='cti_chi'      , fill_array=top_init_vars%cti_chi)
      ptr_2D => input_file%Read_2d(variable_name='cti_basin'    , fill_array=top_init_vars%cti_basin)
      ptr_2D => input_file%Read_2d(variable_name='cti_openwater', fill_array=top_init_vars%cti_openwater)

      CALL input_file%Close()

    ELSE

      input_file = jsb_netcdf_open_input(TRIM(dsl4jsb_Config(TOP_)%bc_filename), model%grid_id)

      srec = (/1,12/) ! Default for climatological files
      ntime = -1
      IF (my_process_is_mpi_workroot()) THEN
        ntime = input_file%Get_dimlen("time")
      ENDIF
      CALL p_bcast(ntime,p_io,mpi_comm)
      IF (IAND(dsl4jsb_config(TOP_)%inun_mode,3)==2) THEN ! This is not a climatological file - calculate relvante times
        ! Unfortunately a mismatch between the parallel netcdf library and the mapped routines, it is not possible to
        ! read the time variable using the standard routines - thus we need to call the low-level nf_ routines on every PE
        ! and ignore errors from PEs which have not opened the file.
        ALLOCATE(rtime(ntime))
        stat = nf_inq_varid(input_file%file_id, "time", varid)
        stat = nf_get_var_double(input_file%file_id, varid, rtime(:))
        CALL p_bcast(rtime,p_io,mpi_comm)
        CALL get_date_components(time_config%tc_startdate,YEAR=syear, MONTH=smonth)
        CALL get_date_components(time_config%tc_stopdate,YEAR=eyear, MONTH=emonth)
        ttime = (/REAL(syear,wp)*10000._wp+REAL(smonth,wp)*100._wp+01._wp, & ! This is for CDO "Absolute dates" only
                  REAL(eyear,wp)*10000._wp+REAL(emonth,wp)*100._wp+01._wp/)
        CALL message(routine,'Times read: '//r2s(rtime(1))//'-'//r2s(rtime(ntime)))
        CALL message(routine,'Times modl: '//r2s(ttime(1))//'-'//r2s(ttime(2)))
        ! TODO: Need some test for time limits not overlapping
        iRec = 1;
        DO iDst = 1, 2
          DO WHILE ((iRec <= ntime) .AND. (ttime(iDst) > rtime(iRec)))
            iRec = iRec + 1
          ENDDO
          srec(iDst) = MIN(iRec,ntime)
        ENDDO
        IF (srec(1)>1) THEN
          dsl4jsb_config(TOP_)%inun_mode = IOR(dsl4jsb_config(TOP_)%inun_mode,8)
        ENDIF
        srec(1) = MAX(1,srec(1)-1) ! Better interpolation if we know the values of the month before
        dsl4jsb_config(TOP_)%step_offset = srec(1)-1

        DEALLOCATE(rtime)

        CALL message(routine,'Using inundation records: '//i2s(srec(1))//'-'//i2s(srec(2)))

      ENDIF

      ! If inundation info is only available on box level: read info now
      IF (IAND(dsl4jsb_config(TOP_)%inun_mode,4)==0) THEN

        ALLOCATE(top_init_vars%frac_inun_presc(nproma, nblks, srec(2)-srec(1)+1))
        top_init_vars%frac_inun_presc(:,:,:) = missval

        ptr_3D => input_file%Read_2d_time(variable_name='frac_inun', fill_array=top_init_vars%frac_inun_presc, &
          start_time_step=srec(1), end_time_step=srec(2))

        IF (input_file%Has_var("wt_inun")) THEN
          ALLOCATE(top_init_vars%wt_inun_presc(nproma, nblks, srec(2)-srec(1)+1))
          top_init_vars%wt_inun_presc(:,:,:) = missval

          ptr_3D => input_file%Read_2d_time(variable_name='wt_inun', fill_array=top_init_vars%wt_inun_presc, &
            start_time_step=srec(1), end_time_step=srec(2))
        ENDIF

      ENDIF

      CALL input_file%Close()

    ENDIF

  END SUBROUTINE topmodel_read_init_vars

  SUBROUTINE topmodel_init_bc(tile)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile

    dsl4jsb_Def_config(TOP_)
    dsl4jsb_Def_memory(TOP_)

    TYPE(t_jsb_model), POINTER :: model
    TYPE(t_jsb_grid),  POINTER :: hgrid
    INTEGER                    :: step, n, iblk

    CHARACTER(len=*), PARAMETER :: routine = modname//':topmodel_init_bc'

    model => get_model(tile%owner_model_id)

    ! Get soil config
    dsl4jsb_Get_config(TOP_)

    IF (.NOT. tile%Is_process_active(TOP_)) RETURN

    hgrid   => Get_grid(model%grid_id)

    IF (debug_on() .OR. dsl4jsb_Config(TOP_)%dbg_process) &
      CALL message(TRIM(routine), 'Setting top-model boundary conditions for tile '//TRIM(tile%name))

    IF (dsl4jsb_Config(TOP_)%inun_mode > 0) THEN

      ! Get soil memory of the tile
      dsl4jsb_Get_memory(TOP_)

      ALLOCATE(top__mem%frac_inun_presc(hgrid%nproma, hgrid%nblks, SIZE(top_init_vars%frac_inun_presc,3)))
      !$acc enter data pcreate(top__mem%frac_inun_presc)
      !$acc enter data attach(top__mem%frac_inun_presc)
      !
      dsl4jsb_memory(TOP_)%frac_inun_presc(:,:,:) = top_init_vars%frac_inun_presc(:,:,:)
      !$acc update device(top__mem%frac_inun_presc)

      n = hgrid%npromz
      DO iblk=hgrid%nblks,1,-1
        dsl4jsb_memory(TOP_)%interp_frac_inun(iblk)%p => interpol_reset(n)
        n = hgrid%nproma
      ENDDO

      IF (ASSOCIATED(top_init_vars%wt_inun_presc)) THEN
        ALLOCATE(top__mem%wt_inun_presc(hgrid%nproma, hgrid%nblks, SIZE(top_init_vars%wt_inun_presc,3)))
        !$acc enter data pcreate(top__mem%wt_inun_presc)
        !$acc enter data attach(top__mem%wt_inun_presc)
        !
        dsl4jsb_memory(TOP_)%wt_inun_presc(:,:,:) = top_init_vars%wt_inun_presc(:,:,:)
        !$acc update device(top__mem%wt_inun_presc)
        n = hgrid%npromz
        DO iblk=hgrid%nblks,1,-1
          dsl4jsb_memory(TOP_)%interp_wt_inun(iblk)%p => interpol_reset(n)
          n = hgrid%nproma
        ENDDO
      ENDIF

      step = 1
      IF (IAND(dsl4jsb_config(TOP_)%inun_mode,8) /= 0) step = 2
      dsl4jsb_memory(TOP_)%step(:) = step

    ENDIF

  END SUBROUTINE topmodel_init_bc

  SUBROUTINE topmodel_init_ic(tile)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile

    TYPE(t_jsb_model), POINTER :: model

    dsl4jsb_Def_config(TOP_)
    dsl4jsb_Def_memory(TOP_)

    INTEGER :: iblk

    CHARACTER(len=*), PARAMETER :: routine = modname//':topmodel_init_ic'

    model => get_model(tile%owner_model_id)

    ! Get soil memory of the tile
    dsl4jsb_Get_memory(TOP_)

    ! Get soil config
    dsl4jsb_Get_config(TOP_)

    IF (.NOT. tile%Is_process_active(TOP_)) RETURN

    IF (debug_on() .OR. dsl4jsb_Config(TOP_)%dbg_process) &
      CALL message(TRIM(routine), 'Initializing top-model memory for tile '//TRIM(tile%name))

    IF (dsl4jsb_Config(TOP_)%inun_mode == 0) THEN
      dsl4jsb_memory(TOP_)%cti_mean     (:,:) = top_init_vars%cti_mean(:,:)
      dsl4jsb_memory(TOP_)%cti_mu       (:,:) = top_init_vars%cti_mu(:,:)
      dsl4jsb_memory(TOP_)%cti_phi      (:,:) = top_init_vars%cti_phi(:,:)
      dsl4jsb_memory(TOP_)%cti_chi      (:,:) = top_init_vars%cti_chi(:,:)
      dsl4jsb_memory(TOP_)%cti_basin    (:,:) = top_init_vars%cti_basin(:,:)
      dsl4jsb_memory(TOP_)%cti_openwater(:,:) = top_init_vars%cti_openwater(:,:)

      !$acc update device(top__mem%cti_mean,top__mem%cti_basin,top__mem%cti_openwater)
      !$acc update device(top__mem%cti_mu,top__mem%cti_phi,top__mem%cti_chi)
    ENDIF

  END SUBROUTINE topmodel_init_ic

#endif
END MODULE mo_top_init
