!>
!! @brief Contains the utilities for the sse process
!!
!! @author
!!  Reiner Schnur, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Reiner Schnur (2013-04-04)
!! Second version by Jan Engel (2016-12-08 jsbach4:git) - integration into SLM, no code changes though
!!
MODULE mo_sse_util
#ifndef __NO_JSBACH__

#include "dsl_jsbach.inc"

  USE mo_kind,          ONLY: wp, zero
  USE mo_exception,     ONLY: finish, message, message_text
  USE mo_jsb_control,   ONLY: debug_on

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: soil_depth_to_layers_2d, calc_qmax_texture
  PUBLIC :: reset_sse_fluxes           ! called in the land model interface prior to all Tasks

  CHARACTER(len=*), PARAMETER :: modname = 'mo_sse_util'

CONTAINS

  !-----------------------------------------------------------------------------------------------------
  !> FUNCTION soil_depth_to_layers_2d
  !! 
  !! out: soil_depth_sl
  !-----------------------------------------------------------------------------------------------------
  FUNCTION soil_depth_to_layers_2d(soil_depth, dz) RESULT(soil_depth_sl)

    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    REAL(wp), INTENT(in) :: soil_depth(:,:)
    REAL(wp), INTENT(in) :: dz(:)
    REAL(wp)             :: soil_depth_sl(SIZE(soil_depth,1), SIZE(dz), SIZE(soil_depth,2))
    ! ---------------------------
    ! 0.2 Local
    REAL(wp) :: dz_cumsum(SIZE(dz)), dz_cumsum_m1(SIZE(dz))
    INTEGER  :: nsoil, isoil
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':soil_depth_to_layer_2d'


    nsoil = SIZE(dz)

    DO isoil = 2,nsoil
      dz_cumsum(isoil)    = SUM(dz(1:isoil))
      dz_cumsum_m1(isoil) = SUM(dz(1:isoil-1))
    END DO
    dz_cumsum(1)    = dz(1)
    dz_cumsum_m1(1) = 0._wp

    soil_depth_sl(:,:,:) = 0._wp

    DO isoil = 1,nsoil
      WHERE (soil_depth(:,:) > dz_cumsum(isoil))
        soil_depth_sl(:,isoil,:) = dz(isoil)
      ELSEWHERE (soil_depth(:,:) > dz_cumsum_m1(isoil))
        soil_depth_sl(:,isoil,:) = soil_depth(:,:) - dz_cumsum_m1(isoil)
      END WHERE
    END DO

  END FUNCTION soil_depth_to_layers_2d

  !-----------------------------------------------------------------------------------------------------
  !> FUNCTION calc_qmax_texture
  !!
  !! Output: sorption capacity for organic material 
  !-----------------------------------------------------------------------------------------------------
  ELEMENTAL FUNCTION calc_qmax_texture(qmax_fine_particle, silt_sl, clay_sl) RESULT(qmax_texture)


    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    REAL(wp), INTENT(in) :: qmax_fine_particle    !< maximum sorption capacity of fine soil particle [mol / kg fine particle (silt or clay)]
    REAL(wp), INTENT(in) :: silt_sl               !< silt content of the mineral soil [kg silt / kg mineral soil]
    REAL(wp), INTENT(in) :: clay_sl               !< clay content of the mineral soil [kg clay / kg mineral soil]
    REAL(wp)             :: qmax_texture          !< maximum sorption capacity of the soil [mol / kg]
    ! ---------------------------
    ! 0.2 Local
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_qmax_texture'

    ! maximum sorption capacity per kg mineral soil (mol/kg)
    qmax_texture = qmax_fine_particle * (silt_sl + clay_sl)

  END FUNCTION calc_qmax_texture

  !-----------------------------------------------------------------------------------------------------
  !> init/reset soil fluxes (with/to zero)
  !! 
  !! called in the land model interface prior to all Tasks
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE reset_sse_fluxes(tile, options)

    USE mo_jsb_class,             ONLY: Get_model
    USE mo_jsb_tile_class,        ONLY: t_jsb_tile_abstract
    USE mo_jsb_task_class,        ONLY: t_jsb_task_options
    USE mo_jsb_model_class,       ONLY: t_jsb_model
    USE mo_jsb_process_class,     ONLY: A2L_, L2A_, SEB_, TURB_, SSE_, HYDRO_, VEG_, HD_, RAD_, ASSIMI_, PHENO_, SB_

    ! Use of process memories
    dsl4jsb_Use_memory(SSE_)        !  USE mo_sse_memory_class,      ONLY: t_sse_memory

    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    CLASS(t_jsb_tile_abstract), INTENT(inout)     :: tile         !< one tile with data structure for one lct
    TYPE(t_jsb_task_options),   INTENT(in)        :: options      !< model options
    ! ---------------------------
    ! 0.2 Local 
    TYPE(t_jsb_model),      POINTER       :: model
    !TYPE(t_lctlib_element), POINTER       :: lctlib         !< land-cover-type library - parameter across pft's
    INTEGER                               :: iblk, ics, ice, nc
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':reset_sse_fluxes'
    ! ---------------------------
    ! 0.3 Declare Memory
    ! Declare process configuration and memory Pointers
    dsl4jsb_Def_memory(SSE_)                          ! TYPE(t_sse_memory),      POINTER :: sse__mem
    ! Declare pointers to variables in memory
    ! SSE_ 2D
    dsl4jsb_Real2D_onChunk      :: transpiration
    ! Get local variables from options argument
    iblk    = options%iblk
    ics     = options%ics
    ice     = options%ice
    nc      = options%nc
    ! ---------------------------
    ! 0.4 Process Activity, Debug Option
    IF (.NOT. tile%Is_process_active(SSE_)) RETURN
    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')
    ! ---------------------------
    ! 0.5 Get Memory
    ! Get process config
!    dsl4jsb_Get_config(SSE_)
    ! Get process memories
    dsl4jsb_Get_memory(SSE_)           ! sse__mem => sse_memory
    ! Get process variables (Set pointers to variables in memory)
    model  => Get_model(tile%owner_model_id)
    !lctlib => model%lctlib(tile%lcts%lib_id)
    ! ---------------------------
    dsl4jsb_Get_var2D_onChunk(SSE_, transpiration)



    ! ------------------------------------------------------------------------------------------------------------
    ! Go science
    ! ------------------------------------------------------------------------------------------------------------


    ! flux variables
    transpiration(:) = zero                  !< transpiration in kg m-2 s-1

    ! other variables
    ! ...


  END SUBROUTINE reset_sse_fluxes

#endif
END MODULE mo_sse_util
