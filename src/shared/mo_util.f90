!> Some utility functions
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
!>#### Various helper routines for ICON-Land
!>
MODULE mo_util
#ifndef __NO_JSBACH__

  USE mo_kind,            ONLY: wp
  USE mo_exception,       ONLY: message, finish
  USE mo_jsb_utils_iface, ONLY: toupper, tolower, int2string, real2string, logical2string, separator, one_of_str
  USE mo_math_utilities,  ONLY: tdma_solver_vec ! TODO from ICON, mo_math_utilities should be replaced once
                                               !      libmath-support is an external library

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: soil_depth_to_layers_2d, ifs2soil, one_of, &
    &       report_memory_usage, soil_init_from_texture, tdma_solver_vec, scaled_global_sum_2D
  PUBLIC :: toupper, tolower, int2string, real2string, logical2string, separator

#ifndef __QUINCY_STANDALONE__
  PUBLIC :: read_2D_var, read_3D_var
#endif

  INTERFACE one_of
    MODULE PROCEDURE one_of_str
    MODULE PROCEDURE one_of_int
  END INTERFACE one_of

  CHARACTER(len=*), PARAMETER :: modname = 'mo_util'

CONTAINS

  !
  !> Compute soil layer depths from (fixed) layer thicknesses and total depth until bedrock
  !>
  FUNCTION soil_depth_to_layers_2d(depth, dz) RESULT(depth_l)

    REAL(wp), INTENT(in) :: depth(:,:)
    REAL(wp), INTENT(in) :: dz(:)
    REAL(wp)             :: depth_l(SIZE(depth,1), SIZE(dz), SIZE(depth,2))

    REAL(wp) :: bottom, bottom_m1
    INTEGER  :: nsoil, i

    CHARACTER(len=*), PARAMETER :: routine = modname//':soil_depth_to_layer_2d'

    nsoil = SIZE(dz)

    depth_l(:,:,:) = 0._wp
    bottom = 0._wp
    bottom_m1 = 0._wp
    DO i=1,nsoil
      bottom = bottom + dz(i)
      WHERE (depth(:,:) >= bottom)
        depth_l(:,i,:) = dz(i)
      ELSEWHERE (depth(:,:) > bottom_m1)
        depth_l(:,i,:) = depth(:,:) - bottom_m1
      END WHERE
      bottom_m1 = bottom_m1 + dz(i)
    END DO

  END FUNCTION soil_depth_to_layers_2d

  SUBROUTINE ifs2soil(ifs_3d, ifs_depth, soil_3d, soil_depth, ifs_sfc)

    REAL(wp), INTENT(in) :: &
      & ifs_3d (:,:,:),     &
      & ifs_depth(:),       &
      & soil_depth(:)
    REAL(wp), OPTIONAL, INTENT(in) :: &
      & ifs_sfc(:,:)

    REAL(wp), INTENT(out) :: &
      & soil_3d(:,:,:)

    REAL(wp) :: ifs(SIZE(ifs_3d,1),0:SIZE(ifs_depth)+1,SIZE(ifs_3d,3)), zifs(0:SIZE(ifs_depth)+1), zsoil(0:SIZE(soil_depth))
    REAL(wp) :: zwt, zwt_sum
    INTEGER  :: nsoil_in, nsoil_out, isoil_in, isoil_out, ik1

    CHARACTER(len=*), PARAMETER :: routine = modname//':ifs2soil'

    nsoil_in  = SIZE(ifs_depth)
    nsoil_out = SIZE(soil_depth)

    zsoil(1:nsoil_out) = soil_depth(:)
    zsoil(0)           = 0._wp
    zifs(1:nsoil_in) = ifs_depth(:)
    zifs(0)          = 0._wp
    zifs(nsoil_in+1) = soil_depth(nsoil_out)
    ifs(:,1:nsoil_in,:) = ifs_3d(:,:,:)
    IF (PRESENT(ifs_sfc)) THEN
      ifs(:,0         ,:) = ifs_sfc(:,:)
    ELSE
      ifs(:,0         ,:) = ifs_3d(:,1,:)
    END IF
    ifs(:,nsoil_in+1,:) = ifs_3d(:,nsoil_in,:)

    soil_3d(:,:,:) = 0._wp

    IF (PRESENT(ifs_sfc)) THEN
      soil_3d(:,1,:) = ifs_sfc(:,:)
      ik1 = 2
    ELSE
      ik1 = 1
    END IF
    DO isoil_out = ik1, nsoil_out
      zwt_sum        = 0._wp
      DO isoil_in = 1,nsoil_in+1
        IF (zifs(isoil_in) <= zsoil(isoil_out-1)) THEN
          ! Do nothing
          zwt = 0._wp
        ELSE IF (zifs(isoil_in-1) > zsoil(isoil_out)) THEN
          ! Do nothing
          zwt = 0._wp
        ELSE IF (zifs(isoil_in-1) > zsoil(isoil_out-1) .AND. zifs(isoil_in) <= zsoil(isoil_out)) THEN
          zwt = 1._wp
          zwt_sum = zwt_sum + zwt
          soil_3d(:,isoil_out,:) = soil_3d(:,isoil_out,:) + zwt * ifs(:,isoil_in,:)
        ELSE IF (zifs(isoil_in) > zsoil(isoil_out-1) .AND. zifs(isoil_in) <= zsoil(isoil_out)) THEN
          zwt = (zifs(isoil_in) - zsoil(isoil_out-1)) / (zsoil(isoil_out) - zsoil(isoil_out-1))
          zwt_sum = zwt_sum + zwt
          soil_3d(:,isoil_out,:) = soil_3d(:,isoil_out,:) + zwt * ifs(:,isoil_in,:)
        ELSE IF (zifs(isoil_in) > zsoil(isoil_out) .AND. zifs(isoil_in-1) <= zsoil(isoil_out)) THEN
          zwt = (zsoil(isoil_out) - zifs(isoil_in-1)) / (zsoil(isoil_out) - zsoil(isoil_out-1))
          zwt_sum = zwt_sum + zwt
          soil_3d(:,isoil_out,:) = soil_3d(:,isoil_out,:) + zwt * ifs(:,isoil_in,:)
        END IF
        ! print*, 'AAA isoil_out, isoil_in, zwt, zwt_sum', isoil_out, isoil_in, zwt, zwt_sum
      END DO
      IF (zwt_sum < EPSILON(1._wp)) CALL finish(routine, 'zwt_sum = 0 (should not happen)!')
      soil_3d(:,isoil_out,:) = soil_3d(:,isoil_out,:) / zwt_sum
    END DO

  END SUBROUTINE ifs2soil

  !> Calculate mineral soil parameters fom soil texture.
  !!
  !! Input to this routine should conform:
  !! fraction_sand + fraction_silt + fraction_clay = 1
  !! fraction_sand_deep + fraction_silt_deep + fraction_clay_deep = 1
  !! The calculated mineral soil parameters do not have a vertical dimension,
  !! as calculations are based on averages of the upper and lower soil textures.
  !! TODO: Calculate layer dependant parameters
  PURE ELEMENTAL SUBROUTINE soil_init_from_texture(                 &
    & parameter_sand, parameter_silt, parameter_clay, parameter_oc, &
    & fraction_sand, fraction_silt, fraction_clay,                  &
    & fraction_sand_deep, fraction_silt_deep, fraction_clay_deep,   &
    & parameter_from_texture)

    REAL(wp), INTENT(in)  ::  parameter_sand, parameter_silt, parameter_clay, parameter_oc
    REAL(wp), INTENT(in)  ::                                      &
      fraction_sand, fraction_silt, fraction_clay,                &
      fraction_sand_deep, fraction_silt_deep, fraction_clay_deep
    REAL(wp), INTENT(out)  ::  parameter_from_texture

    ! Compute texture based soil parameter
    parameter_from_texture = ((fraction_sand + fraction_sand_deep)         &
      &                      * 0.5_wp * parameter_sand) +                  &
      &                      ((fraction_silt + fraction_silt_deep)         &
      &                      * 0.5_wp * parameter_silt) +                  &
      &                      ((fraction_clay + fraction_clay_deep)         &
      &                      * 0.5_wp * parameter_clay)

    ! Ensure the parameter is within the correct bounds
    parameter_from_texture = MAX(parameter_from_texture, MIN(parameter_sand,parameter_silt,parameter_clay))
    parameter_from_texture = MIN(parameter_from_texture, MAX(parameter_sand,parameter_silt,parameter_clay))

  END SUBROUTINE soil_init_from_texture

  FUNCTION one_of_int(in_int, arg)
    INTEGER :: one_of_int
    INTEGER, INTENT(IN)           :: in_int    ! input integer
    INTEGER, INTENT(IN)           :: arg(:)
    ! local variables:
    INTEGER :: i

    one_of_int = -1
    DO i=1,SIZE(arg)
      IF (in_int == arg(i)) THEN
        one_of_int=i
        EXIT
      END IF
    END DO
  END FUNCTION one_of_int

  SUBROUTINE report_memory_usage(msg)

    CHARACTER(len=*), INTENT(in), OPTIONAL:: msg     ! Message tag to be added to the origial OOM message

    INTEGER           :: size=0, resident=0, share=0, text=0, lib=0, data=0, dt=0, oom_score=0
    INTEGER           :: ierr=0

    OPEN(unit=1, file='/proc/self/statm', form='formatted', status='OLD', action='READ', iostat=ierr)
    READ(unit=1, fmt=*, iostat=ierr) size, resident, share, text, lib, data, dt
    CLOSE(unit=1)
    OPEN(unit=1, file='/proc/self/oom_score', form='formatted', status='OLD', action='READ', iostat=ierr)
    READ(unit=1, fmt=*, iostat=ierr) oom_score
    CLOSE(unit=1)

    IF (ierr < 0) THEN
      CALL message('', 'Problem reading /proc/self/statm')
    ELSE
      IF (PRESENT(msg)) THEN
        CALL message('INFO -- Memory usage in pages and OOM score', ' size='//TRIM(int2string(size))// &
          &                             ' resident='//TRIM(int2string(resident))// &
          &                             ' data='//TRIM(int2string(data))// &
          &                             ' oom_score='//TRIM(int2string(oom_score))//' - '//msg)
      ELSE
        CALL message('INFO -- Memory usage in pages and OOM score', ' size='//TRIM(int2string(size))// &
          &                             ' resident='//TRIM(int2string(resident))// &
          &                             ' data='//TRIM(int2string(data))// &
          &                             ' oom_score='//TRIM(int2string(oom_score)))
      END IF
    END IF

  END SUBROUTINE report_memory_usage


  ! ====================================================================================================== !
  !
  !> Calculate the scaled global sum of a given 2D variable using passed 2D and 1D scaling factors
  !
  ! ======================================================================================================= !
  SUBROUTINE scaled_global_sum_2D(scaled_global_sum, var_ptr, scaling_2D, scaling_factor, acc_stream)

    USE mo_sync,                    ONLY: global_sum_array
    ! ----------------------------------------------------------------------------------------------------- !
    REAL(wp),          INTENT(out) :: scaled_global_sum(:) !< Scaled variable
    REAL(wp), POINTER, INTENT(in)  :: var_ptr(:,:)         !< Variable pointer
    REAL(wp),          INTENT(in)  :: scaling_2D(:,:)      !< 2D scaling factor
    REAL(wp),          INTENT(in)  :: scaling_factor       !< Additional scaling factor
    INTEGER,           INTENT(in)  :: acc_stream           !< ACC stream - passed to avoid cyclic
                                                           !< dependencies
    ! ----------------------------------------------------------------------------------------------------- !
    REAL(wp) :: scaled_variable(SIZE(scaling_2D,1), SIZE(scaling_2D,2))
    INTEGER  :: blk, jc, nblks, nproma
    ! ----------------------------------------------------------------------------------------------------- !
    nblks  = SIZE(scaling_2D,2)
    nproma = SIZE(scaling_2D,1)

    !$ACC DATA ASYNC(acc_stream) CREATE(scaled_variable)

    !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(acc_stream) COLLAPSE(2)
    DO blk = 1, nblks
      DO jc = 1, nproma
        scaled_variable(jc,blk) = scaling_2D(jc,blk) * var_ptr(jc,blk)
      END DO
    END DO
    !$ACC END PARALLEL LOOP
    scaled_global_sum(:) = global_sum_array(scaled_variable(:,:), lacc=.TRUE.) * scaling_factor

    !$ACC END DATA

  END SUBROUTINE scaled_global_sum_2D


#ifndef __QUINCY_STANDALONE__
  ! ====================================================================================================== !
  !
  !> Small helper to read given 2D variable from given netcdf file
  !
  ! ======================================================================================================= !
  SUBROUTINE read_2D_var(input_file, ptr2D, var_name)
    USE mo_jsb_io_netcdf,           ONLY: t_input_file
    ! ----------------------------------------------------------------------------------------------------- !
    TYPE(t_input_file),      INTENT(INOUT) :: input_file !< file to read passed var from
    REAL(wp), POINTER,       INTENT(IN)    :: ptr2D(:,:) !< var pointer
    CHARACTER(*),            INTENT(IN)    :: var_name   !< name of the variable in the file
    ! ----------------------------------------------------------------------------------------------------- !
    REAL(wp), POINTER :: return_2D(:,:)   !< temporary pointer used to "call" the netcdf reading routine
    ! ----------------------------------------------------------------------------------------------------- !

    return_2D => input_file%Read_2d(variable_name = var_name)
    ptr2D(:,:) = return_2D(:,:)

  END SUBROUTINE read_2D_var

  ! ====================================================================================================== !
  !
  !> Small helper to read given 3D variable from given netcdf file
  !
  ! ======================================================================================================= !
  SUBROUTINE read_3D_var(input_file, ptr3D, var_name, name_extdim, size_extdim)
    USE mo_jsb_io_netcdf,           ONLY: t_input_file
    ! ----------------------------------------------------------------------------------------------------- !
    TYPE(t_input_file),      INTENT(INOUT) :: input_file   !< file to read passed var from
    REAL(wp), POINTER,       INTENT(IN)    :: ptr3D(:,:,:) !< var pointer
    CHARACTER(*),            INTENT(IN)    :: var_name     !< name of the variable in the file
    CHARACTER(*),            INTENT(IN)    :: name_extdim  !< name of the variable in the file
    INTEGER,                 INTENT(IN)    :: size_extdim  !< size of additional dim
    ! ----------------------------------------------------------------------------------------------------- !
    INTEGER           :: i_dim3            !< third dim index
    REAL(wp), POINTER :: return_3D(:,:, :) !< temporary pointer used to "call" the netcdf reading routine
    ! ----------------------------------------------------------------------------------------------------- !

    DO i_dim3 = 1,size_extdim
      return_3D => input_file%Read_2d_extdim(variable_name=var_name, &
        & start_extdim=i_dim3, end_extdim=i_dim3, extdim_name=name_extdim)
      ptr3D(:,i_dim3,:) = return_3D(:,:,1)
    END DO

  END SUBROUTINE read_3D_var
#endif

#endif
END MODULE mo_util
