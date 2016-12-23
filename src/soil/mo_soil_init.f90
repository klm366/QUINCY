!>
!! @brief Initialization of the the soil memory.
!!
!! @author
!!  Jan Engel, MPI-BGC Jena
!!
!! @par Revision History
!! First version for SLM by Jan Engel (20161208) based on the JBACH4 version by Reiner Schnur (2013-04-04)
!!
MODULE mo_soil_init

#include "dsl_jsbach.inc"

  USE mo_kind,                  ONLY: wp
  USE mo_exception,             ONLY: message


  PUBLIC :: soil_init

  CHARACTER(len=*), PARAMETER :: modname = 'mo_soil_init'


CONTAINS

  !-----------------------------------------------------------------------------------------------------
  !> Intialize soil process (after memory has been set up)
  !! 
  !! 
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE soil_init(tile)

    USE mo_jsb_class,             ONLY: Get_model
    USE mo_base_memory,           ONLY: t_slm_tile, t_model
    USE mo_jsb_lctlib_class,      ONLY: t_lctlib_element
    USE mo_jsb_process_class,     ONLY: ATM2LAND_, LAND2ATM_, SEB_, TURB_, SOIL_, HYDRO_, VEG_, HD_, RAD_, ASSIMI_, PHENO_
    USE mo_soil_config_class,     ONLY: soil_config
    USE mo_soil_constants
    USE mo_soil_util,             ONLY: soil_depth_to_layers_2d

    ! USE PROC_memory
    dsl4jsb_Use_memory(SOIL_)


  IMPLICIT NONE

    ! 0.1 InOut
    TYPE(t_slm_tile),       INTENT(inout)        :: tile       !< data structure for dynamic vegetation and soil variables
    ! 0.2 Local
    TYPE(t_model),          POINTER              :: model
    TYPE(t_lctlib_element), POINTER              :: lctlib     !< land-cover-type library - parameter across pft's
    ! 0.3 Memory allocation
    ! Declare process configuration and memory Pointers
    dsl4jsb_Def_memory(SOIL_)
    TYPE(t_soil_memory),    POINTER              :: soil       ! local pointer to soil__mem for SLM usability
    dsl4jsb_Get_memory(SOIL_)
    soil   => soil__mem 
    model  => Get_model(tile%owner_model_id)
    lctlib => model%lctlib(tile%lcts%lib_id)


!    ! Total soil depth until bedrock
!    return_pointer => input_file%Read_2d(         &
!      & variable_name='soil_depth',               &
!      & fill_array = dsl4jsb_var2D_ptr(SOIL_, soil_depth))
!    dsl4jsb_var2D_onDomain(SOIL_, soil_depth) = &
!      MAX(0.1_wp, MERGE(dsl4jsb_var2D_onDomain(SOIL_, soil_depth), 0._wp, dsl4jsb_var2D_onDomain(SOIL_, soil_depth) >= 0._wp))

    ! Compute actual soil layer depths based on the fixed layer thicknesses from namelist and
    ! the bedrock depth just read in.
!    dsl4jsb_var3D_onDomain(SOIL_, soil_depth_l) = soil_depth_to_layers_2d( &
!      & dsl4jsb_var2D_onDomain(SOIL_, soil_depth),                         & ! Total soil depth until bedrock (from textures)
!      & dsl4jsb_Config(SOIL_)%dz_water(1:dsl4jsb_Config(SOIL_)%nsoil_water))    ! Soil layer thicknesses from namelist
   

    soil%soil_depth%ptr = soil_config%soil_depth 
    soil%soil_depth_l%ptr(:,:,:) = soil_depth_to_layers_2d(soil%soil_depth%ptr(:,:), soil_config%dz_water(1:soil_config%nsoil_water))

    soil%bulk_dens_l%ptr  = rho_bulk_soil 
    soil%heat_capa_l%ptr  = soil_heat_cap
    soil%therm_cond_l%ptr = soil_therm_cond



  END SUBROUTINE soil_init

END MODULE mo_soil_init
