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

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: soil_init

  CHARACTER(len=*), PARAMETER :: modname = 'mo_soil_init'

CONTAINS

  !-----------------------------------------------------------------------------------------------------
  !> Intialize soil process (after memory has been set up)
  !! 
  !! 
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE soil_init(tile, options)

    USE mo_jsb_class,             ONLY: Get_model
    USE mo_base_memory,           ONLY: t_slm_tile, t_jsb_task_options
    USE mo_jsb_model_class,       ONLY: t_jsb_model
    USE mo_jsb_lctlib_class,      ONLY: t_lctlib_element
    USE mo_jsb_process_class,     ONLY: A2L_, L2A_, SEB_, TURB_, SOIL_, HYDRO_, VEG_, HD_, RAD_, ASSIMI_, PHENO_, SB_
    USE mo_soil_constants,        ONLY: theta_sat, theta_wilt, rho_bulk_soil, soil_heat_cap, soil_therm_cond
    USE mo_soil_util,             ONLY: soil_depth_to_layers_2d

    ! Use of process configurations (t_PROC_config)
    dsl4jsb_Use_config(SOIL_)

    ! USE PROC_memory
    dsl4jsb_Use_memory(SOIL_)


    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    TYPE(t_slm_tile),           INTENT(inout)     :: tile         !< data structure for dynamic phenoetation and soil variables
    TYPE(t_jsb_task_options),   INTENT(in)        :: options      !< model options
    ! ---------------------------
    ! 0.2 Local
    TYPE(t_jsb_model),      POINTER              :: model
    TYPE(t_lctlib_element), POINTER              :: lctlib     !< land-cover-type library - parameter across pft's
    INTEGER                                      :: isoil
    INTEGER                                      :: iblk, ics, ice, nc
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':soil_init'
    ! ---------------------------
    ! 0.3 Declare Memory
    ! Declare process configuration and memory Pointers
    dsl4jsb_Def_config(SOIL_)
    dsl4jsb_Def_memory(SOIL_)
    ! Declare pointers to variables in memory
    dsl4jsb_Real2D_onChunk    :: soil_depth
    dsl4jsb_Real2D_onChunk    :: root_depth
    dsl4jsb_Real2D_onChunk    :: w_skin
    dsl4jsb_Real2D_onChunk    :: w_soil_root
    dsl4jsb_Real2D_onChunk    :: w_soil_root_pwp
    dsl4jsb_Real2D_onChunk    :: w_soil_root_fc
    dsl4jsb_Real2D_onChunk    :: w_soil_root_theta
    dsl4jsb_Real3D_onChunk    :: soil_depth_sl
    dsl4jsb_Real3D_onChunk    :: bulk_dens_sl
    dsl4jsb_Real3D_onChunk    :: heat_capa_sl
    dsl4jsb_Real3D_onChunk    :: therm_cond_sl
    dsl4jsb_Real3D_onChunk    :: root_fraction_sl
    dsl4jsb_Real3D_onChunk    :: w_soil_pwp_sl
    dsl4jsb_Real3D_onChunk    :: w_soil_fc_sl 
    dsl4jsb_Real3D_onChunk    :: w_soil_sl
    ! Get local variables from options argument
    iblk    = options%iblk
    ics     = options%ics
    ice     = options%ice
    nc      = options%nc
    ! ---------------------------
    ! 0.4 Process Activity, Debug Option
!    IF (.NOT. tile%Is_process_active(SOIL_)) RETURN
!    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')
    ! ---------------------------
    ! 0.5 Get Memory
    model  => Get_model(tile%owner_model_id)
    lctlib => model%lctlib(tile%lcts%lib_id)
    ! Get process config
    dsl4jsb_Get_config(SOIL_)
    ! Get process memories
    dsl4jsb_Get_memory(SOIL_)
    ! Get process variables (Set pointers to variables in memory)
    dsl4jsb_Get_var2D_onChunk(SOIL_, soil_depth)
    dsl4jsb_Get_var2D_onChunk(SOIL_, root_depth)
    dsl4jsb_Get_var2D_onChunk(SOIL_, w_skin)
    dsl4jsb_Get_var2D_onChunk(SOIL_, w_soil_root)
    dsl4jsb_Get_var2D_onChunk(SOIL_, w_soil_root_pwp)
    dsl4jsb_Get_var2D_onChunk(SOIL_, w_soil_root_fc)
    dsl4jsb_Get_var2D_onChunk(SOIL_, w_soil_root_theta)
    dsl4jsb_Get_var3D_onChunk(SOIL_, soil_depth_sl)
    dsl4jsb_Get_var3D_onChunk(SOIL_, bulk_dens_sl)
    dsl4jsb_Get_var3D_onChunk(SOIL_, heat_capa_sl)
    dsl4jsb_Get_var3D_onChunk(SOIL_, therm_cond_sl)
    dsl4jsb_Get_var3D_onChunk(SOIL_, root_fraction_sl)
    dsl4jsb_Get_var3D_onChunk(SOIL_, w_soil_pwp_sl)
    dsl4jsb_Get_var3D_onChunk(SOIL_, w_soil_fc_sl)
    dsl4jsb_Get_var3D_onChunk(SOIL_, w_soil_sl)



    ! ------------------------------------------------------------------------------------------------------------
    ! Go science
    ! ------------------------------------------------------------------------------------------------------------


    ! Total soil depth until bedrock
!    return_pointer => input_file%Read_2d(         &
!      & variable_name='soil_depth',               &
!      & fill_array = dsl4jsb_var2D_ptr(SOIL_, soil_depth))
    soil_depth(:)     = dsl4jsb_Config(SOIL_)%soil_depth ! SLM code *slmtemp
    dsl4jsb_var2D_onDomain(SOIL_, soil_depth) = &
      MAX(0.1_wp, MERGE(dsl4jsb_var2D_onDomain(SOIL_, soil_depth), 0._wp, dsl4jsb_var2D_onDomain(SOIL_, soil_depth) >= 0._wp))

    ! Compute actual soil layer depths based on the fixed layer thicknesses from namelist and
    ! the bedrock depth just read in.
    dsl4jsb_var3D_onDomain(SOIL_, soil_depth_sl) = soil_depth_to_layers_2d( &
      & dsl4jsb_var2D_onDomain(SOIL_, soil_depth),                         & ! Total soil depth until bedrock (from textures)
      & dsl4jsb_Config(SOIL_)%dz_water(1:options%nsoil_sb))    ! Soil layer thicknesses from namelist

    bulk_dens_sl(:,:)  = rho_bulk_soil 
    heat_capa_sl(:,:)  = soil_heat_cap
    therm_cond_sl(:,:) = soil_therm_cond

    ! rooting depth given by prescribed AWC
    root_depth(:) = MAX(dsl4jsb_Config(SOIL_)%soil_depth_water(1), & 
                          ( dsl4jsb_Config(SOIL_)%soil_awc_prescribe / 1000._wp )/ (theta_sat - theta_wilt)) 

    ! water in the surface reservoir [m]
    w_skin(:) = 0.0_wp

    ! soil water in the rooting zone [m]
    w_soil_root_pwp(:)    = theta_wilt * root_depth(:) 
    w_soil_root_fc(:)     = dsl4jsb_Config(SOIL_)%soil_awc_prescribe / 1000._wp + w_soil_root_pwp(:)
    ! adjust initial water amount to prescribed theta 
    w_soil_root(:)        = w_soil_root_pwp(:) + & 
                            dsl4jsb_Config(SOIL_)%soil_theta_prescribe * (w_soil_root_fc(:) - w_soil_root_pwp(:)) 
    w_soil_root_theta(:)  = dsl4jsb_Config(SOIL_)%soil_theta_prescribe

    ! root fractions assuming an exponential distribution of roots
    ! Cr_l = Cr_0 * exp(-lctlib%k_root_dist*depth)
    ! integral to bottom of root zone: 1/lctlib%k_root_dist * (1 - exp(-lctlib%k_root_dist*depth) 
    ! thus for each layer
    ! layer 1
    root_fraction_sl(:,1) = &  
          (1._wp - exp(-lctlib%k_root_dist * dsl4jsb_Config(SOIL_)%soil_depth_water(1))) / &
          (1._wp - exp(-lctlib%k_root_dist * root_depth(:)))
    ! layer 2 -> options%nsoil_sb
    DO isoil=2,options%nsoil_sb
       WHERE(root_depth(:) > dsl4jsb_Config(SOIL_)%soil_depth_water(isoil-1))
          WHERE(root_depth(:) > dsl4jsb_Config(SOIL_)%soil_depth_water(isoil))
             root_fraction_sl(:,isoil) = & 
                (exp(-lctlib%k_root_dist * dsl4jsb_Config(SOIL_)%soil_depth_water(isoil-1)) - &
                 exp(-lctlib%k_root_dist * dsl4jsb_Config(SOIL_)%soil_depth_water(isoil))) / & 
                (1._wp - exp(-lctlib%k_root_dist * root_depth(:)))
          ELSEWHERE
             root_fraction_sl(:,isoil) = & 
                (exp(-lctlib%k_root_dist * dsl4jsb_Config(SOIL_)%soil_depth_water(isoil-1)) - &
                 exp(-lctlib%k_root_dist * root_depth(:))) / & 
                (1._wp - exp(-lctlib%k_root_dist * root_depth(:)))
          ENDWHERE
       ELSEWHERE
          root_fraction_sl(:,isoil) = 0.0_wp 
       ENDWHERE
    ENDDO


    ! soil water per layer depth
    DO isoil=1,options%nsoil_sb
       WHERE(soil_depth_sl(:,isoil) > 0.0_wp) 
          w_soil_pwp_sl(:,isoil) = theta_wilt * soil_depth_sl(:,isoil) 
          w_soil_fc_sl(:,isoil)  = theta_sat * soil_depth_sl(:,isoil) 
       ENDWHERE
    ENDDO
    ! adjust initial water amount to prescribed theta
    w_soil_sl(:,:) = w_soil_pwp_sl(:,:) + & 
                     dsl4jsb_Config(SOIL_)%soil_theta_prescribe * (w_soil_fc_sl(:,:) - w_soil_pwp_sl(:,:))

  END SUBROUTINE soil_init

END MODULE mo_soil_init
