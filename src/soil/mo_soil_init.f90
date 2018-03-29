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

  USE mo_kind,                  ONLY: wp, eps8
  USE mo_exception,             ONLY: message_text, finish, message
  USE mo_jsb_control,           ONLY: debug_on

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
    USE mo_soil_util,             ONLY: soil_depth_to_layers_2d
    USE mo_soil_constants,        ONLY: soil_heat_cap, soil_therm_cond, kdiff_sat_max, &
                                        k_pwp_s, k_pwp_c, k_pwp_sc, k_pwp_a, k_pwp_at, k_pwp_bt, &
                                        k_fc_s, k_fc_c, k_fc_sc, k_fc_a, k_fc_at, k_fc_bt, k_fc_ct, &
                                        k_sat_s, k_sat_c, k_sat_sc, k_sat_a, k_sat_at, k_sat_bt, k_sat_ct, k_sat_dt, &
                                        wpot_pwp, wpot_fc, wpot_fc

    ! Use of process configurations (t_PROC_config)
    dsl4jsb_Use_config(SOIL_)

    ! USE PROC_memory
    dsl4jsb_Use_memory(SOIL_)
    dsl4jsb_Use_memory(VEG_)

    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    TYPE(t_slm_tile),           INTENT(inout)        :: tile         !< data structure for dynamic phenoetation and soil variables
    TYPE(t_jsb_task_options),   INTENT(in)           :: options      !< model options
    ! ---------------------------
    ! 0.2 Local
    TYPE(t_jsb_model),      POINTER                  :: model
    TYPE(t_lctlib_element), POINTER                  :: lctlib     !< land-cover-type library - parameter across pft's
    REAL(wp), DIMENSION(options%nc,options%nsoil_sb) :: hlp1, &
                                                        theta_pwp_sl, theta_fc_sl, theta_sat_sl
    REAL(wp), DIMENSION(options%nc)                  :: soil_awc, &
                                                        hlp2
    INTEGER                                          :: nsoil_sb, isoil, ichunk
    INTEGER                                          :: iblk, ics, ice, nc
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':soil_init'
    ! ---------------------------
    ! 0.3 Declare Memory
    ! Declare process configuration and memory Pointers
    dsl4jsb_Def_config(SOIL_)
    dsl4jsb_Def_memory(SOIL_)
    dsl4jsb_Def_memory(VEG_)
    ! Declare pointers to variables in memory
    ! VEG_
    dsl4jsb_Real3D_onChunk    :: root_fraction_sl
    ! SOIL_
    dsl4jsb_Real2D_onChunk    :: soil_depth
    dsl4jsb_Real2D_onChunk    :: root_depth
    dsl4jsb_Real2D_onChunk    :: t_srf_old 
    dsl4jsb_Real2D_onChunk    :: w_skin
    dsl4jsb_Real2D_onChunk    :: w_soil_root
    dsl4jsb_Real2D_onChunk    :: w_soil_root_pwp
    dsl4jsb_Real2D_onChunk    :: w_soil_root_fc
    dsl4jsb_Real2D_onChunk    :: w_soil_root_theta
    dsl4jsb_Real2D_onChunk    :: w_soil_root_pot
    dsl4jsb_Real3D_onChunk    :: soil_depth_sl
    dsl4jsb_Real3D_onChunk    :: bulk_dens_sl
    dsl4jsb_Real3D_onChunk    :: heat_capa_sl
    dsl4jsb_Real3D_onChunk    :: therm_cond_sl
    dsl4jsb_Real3D_onChunk    :: w_soil_pwp_sl
    dsl4jsb_Real3D_onChunk    :: w_soil_fc_sl 
    dsl4jsb_Real3D_onChunk    :: w_soil_sat_sl 
    dsl4jsb_Real3D_onChunk    :: w_soil_sl
    dsl4jsb_Real3D_onChunk    :: w_soil_pot_sl
    dsl4jsb_Real3D_onChunk    :: saxtonA 
    dsl4jsb_Real3D_onChunk    :: saxtonB 
    dsl4jsb_Real3D_onChunk    :: saxtonC 
    dsl4jsb_Real3D_onChunk    :: kdiff_sat_sl 
    dsl4jsb_Real3D_onChunk    :: sand_sl
    dsl4jsb_Real3D_onChunk    :: silt_sl
    dsl4jsb_Real3D_onChunk    :: clay_sl
    ! Get local variables from options argument
    iblk      = options%iblk
    ics       = options%ics
    ice       = options%ice
    nc        = options%nc
    nsoil_sb  = options%nsoil_sb
    ! ---------------------------
    ! 0.4 Process Activity, Debug Option
    IF (.NOT. tile%Is_process_active(SOIL_)) RETURN
    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')
    ! ---------------------------
    ! 0.5 Get Memory
    model  => Get_model(tile%owner_model_id)
    lctlib => model%lctlib(tile%lcts%lib_id)
    ! Get process config
    dsl4jsb_Get_config(SOIL_)
    ! Get process memories
    dsl4jsb_Get_memory(SOIL_)
    dsl4jsb_Get_memory(VEG_)
    ! Get process variables (Set pointers to variables in memory)
    dsl4jsb_Get_var3D_onChunk(VEG_,  root_fraction_sl)
    ! ---------------------------
    dsl4jsb_Get_var2D_onChunk(SOIL_, soil_depth)
    dsl4jsb_Get_var2D_onChunk(SOIL_, root_depth)
    dsl4jsb_Get_var2D_onChunk(SOIL_, t_srf_old)
    dsl4jsb_Get_var2D_onChunk(SOIL_, w_skin)
    dsl4jsb_Get_var2D_onChunk(SOIL_, w_soil_root)
    dsl4jsb_Get_var2D_onChunk(SOIL_, w_soil_root_pwp)
    dsl4jsb_Get_var2D_onChunk(SOIL_, w_soil_root_fc)
    dsl4jsb_Get_var2D_onChunk(SOIL_, w_soil_root_theta)
    dsl4jsb_Get_var2D_onChunk(SOIL_, w_soil_root_pot)
    dsl4jsb_Get_var3D_onChunk(SOIL_, soil_depth_sl)
    dsl4jsb_Get_var3D_onChunk(SOIL_, bulk_dens_sl)
    dsl4jsb_Get_var3D_onChunk(SOIL_, heat_capa_sl)
    dsl4jsb_Get_var3D_onChunk(SOIL_, therm_cond_sl)
    dsl4jsb_Get_var3D_onChunk(SOIL_, w_soil_pwp_sl)
    dsl4jsb_Get_var3D_onChunk(SOIL_, w_soil_fc_sl)
    dsl4jsb_Get_var3D_onChunk(SOIL_, w_soil_sat_sl)
    dsl4jsb_Get_var3D_onChunk(SOIL_, w_soil_sl)
    dsl4jsb_Get_var3D_onChunk(SOIL_, w_soil_pot_sl)
    dsl4jsb_Get_var3D_onChunk(SOIL_, saxtonA)
    dsl4jsb_Get_var3D_onChunk(SOIL_, saxtonB)
    dsl4jsb_Get_var3D_onChunk(SOIL_, saxtonC)
    dsl4jsb_Get_var3D_onChunk(SOIL_, kdiff_sat_sl)
    dsl4jsb_Get_var3D_onChunk(SOIL_, sand_sl)
    dsl4jsb_Get_var3D_onChunk(SOIL_, silt_sl)
    dsl4jsb_Get_var3D_onChunk(SOIL_, clay_sl)


    ! ------------------------------------------------------------------------------------------------------------
    ! Go science
    ! ------------------------------------------------------------------------------------------------------------


    !> 1.0 Total soil depth until bedrock
    !! 
!    return_pointer => input_file%Read_2d(         &
!      & variable_name='soil_depth',               &
!      & fill_array = dsl4jsb_var2D_ptr(SOIL_, soil_depth))
    soil_depth(:)                             = dsl4jsb_Config(SOIL_)%soil_depth ! SLM code *slmtemp
    dsl4jsb_var2D_onDomain(SOIL_, soil_depth) = &
      MAX(0.1_wp, MERGE(dsl4jsb_var2D_onDomain(SOIL_, soil_depth), 0._wp, dsl4jsb_var2D_onDomain(SOIL_, soil_depth) >= 0._wp))

    !> 2.0 Compute actual soil layer depths based on the fixed layer thicknesses from namelist and
    !!  the bedrock depth just read in.
    dsl4jsb_var3D_onDomain(SOIL_, soil_depth_sl) = soil_depth_to_layers_2d( &
      & dsl4jsb_var2D_onDomain(SOIL_, soil_depth),                         & ! Total soil depth until bedrock (from textures)
      & dsl4jsb_Config(SOIL_)%dz_water(1:nsoil_sb))    ! Soil layer thicknesses from namelist

    !> 3.0 check if last soil layer has a thickness larger zero
    !! 
    DO ichunk  = 1,nc
      IF (soil_depth_sl(ichunk,nsoil_sb) < eps8) THEN
        WRITE(message_text,*) 'The last soil layer has a thickness < eps8. Increase soil_depth from ', soil_depth, &
                              ' to at least ', SUM(dsl4jsb_Config(SOIL_)%dz_water(1:nsoil_sb-1)), ' but better ', &
                              SUM(dsl4jsb_Config(SOIL_)%dz_water(1:nsoil_sb)), ' and re-start the simulation.'
        CALL finish(TRIM(routine), message_text)
      END IF  
    END DO
    
    !> 4.0 init soil-property variables
    !! 
    bulk_dens_sl(:,:)   = dsl4jsb_Config(SOIL_)%bulk_density 
    heat_capa_sl(:,:)   = soil_heat_cap
    therm_cond_sl(:,:)  = soil_therm_cond
    sand_sl(:,:)        = dsl4jsb_Config(SOIL_)%soil_sand
    silt_sl(:,:)        = dsl4jsb_Config(SOIL_)%soil_silt
    clay_sl(:,:)        = dsl4jsb_Config(SOIL_)%soil_clay
    t_srf_old(:)        = 273.15_wp

    !> 5.0 soil hydraulic properties, based on Saxton and Rawls 2006, SSSA
    !! 

    !> 5.1 fractional water holding capacity at permanent wilting point
    !! note that this assumes that the soil hydraulic properties are identical between the layers
    !!
    hlp1(:,:)         = k_pwp_s * sand_sl(:,:) + & 
                        k_pwp_c * clay_sl(:,:) + &
                        k_pwp_sc * sand_sl(:,:) * clay_sl(:,:) + &
                        k_pwp_a
    theta_pwp_sl(:,:) = k_pwp_at + (1.0_wp + k_pwp_bt) * hlp1(:,:) 

    !> 5.2 fractional water holding capacity at field capacity 
    !! 
    hlp1(:,:)         = k_fc_s * sand_sl(:,:) + &
                        k_fc_c * clay_sl(:,:) + &
                        k_fc_sc * sand_sl(:,:) * clay_sl(:,:) + &
                        k_fc_a
    theta_fc_sl(:,:)  = k_fc_at + (1.0_wp + k_fc_bt) * hlp1(:,:) + k_fc_ct * hlp1(:,:)**2._wp

    !> 5.3 fractional water holding capacity at saturation 
    !! 
    hlp1(:,:)         = k_sat_s * sand_sl(:,:) + & 
                        k_sat_c * clay_sl(:,:) + &
                        k_sat_sc * sand_sl(:,:) * clay_sl(:,:) + &
                        k_sat_a
    theta_sat_sl(:,:) = theta_fc_sl(:,:) + k_sat_at + (1.0_wp + k_sat_bt) * hlp1(:,:) - & 
                        k_sat_ct * sand_sl(:,:) + k_sat_dt 
    
    !> 5.4 parameters describing the moisture-tension curve                 
    !! note that for convenience, B has a negative sign and A is in MPa 
    saxtonB(:,:)      = - (log(wpot_pwp) - log(wpot_fc))/(log(theta_fc_sl(:,:)) - log(theta_pwp_sl(:,:)))
    saxtonA(:,:)      = - (exp(log(wpot_fc) - saxtonB(:,:) * log(theta_fc_sl(:,:))))/1000.0_wp

    !> 5.5 parameters describing the moisture-diffusivity curve
    !! 
    kdiff_sat_sl(:,:) = kdiff_sat_max * (theta_sat_sl(:,:) - theta_fc_sl(:,:)) ** (3.0_wp + 1._wp/saxtonB(:,:))
    saxtonC(:,:)      = 3.0_wp - 2._wp * saxtonB(:,:)

    !> 6.0 rooting depth given by prescribed AWC
    !! soil water in the rooting zone [m] \n
    !! the init of 'root_fraction_sl(:,:)' (sec 6.3) depends on this init of the root_depth(:)
    !! 
    soil_awc(:)         = 0._wp
    w_soil_root_pwp(:)  = 0._wp
    w_soil_root_fc(:)   = 0._wp
    root_depth(:)       = 0._wp
    DO isoil=1,nsoil_sb
      hlp1(:,isoil)     = dsl4jsb_Config(SOIL_)%dz_water(isoil) * (theta_fc_sl(:,isoil) - theta_pwp_sl(:,isoil))
      soil_awc(:)       = soil_awc(:) + hlp1(:,isoil)
      WHERE (soil_awc(:) < (dsl4jsb_Config(SOIL_)%soil_awc_prescribe / 1000._wp)) 
        root_depth(:)       = root_depth(:) + dsl4jsb_Config(SOIL_)%dz_water(isoil)
        w_soil_root_pwp(:)  = w_soil_root_pwp(:) + theta_pwp_sl(:,isoil) * dsl4jsb_Config(SOIL_)%dz_water(isoil)
        w_soil_root_fc(:)   = w_soil_root_fc(:) + theta_fc_sl(:,isoil) * dsl4jsb_Config(SOIL_)%dz_water(isoil)
      ELSEWHERE
        hlp2(:)             = MAX(0._wp, dsl4jsb_Config(SOIL_)%soil_awc_prescribe / 1000._wp - soil_awc(:) + hlp1(:,isoil))
        root_depth(:)       = root_depth(:) + hlp2(:) / (theta_fc_sl(:,isoil) - theta_pwp_sl(:,isoil))
        w_soil_root_pwp(:)  = w_soil_root_pwp(:) + theta_pwp_sl(:,isoil) * hlp2(:) / (theta_fc_sl(:,isoil) - theta_pwp_sl(:,isoil))
        w_soil_root_fc(:)   = w_soil_root_fc(:) + theta_fc_sl(:,isoil) * hlp2(:) / (theta_fc_sl(:,isoil) - theta_pwp_sl(:,isoil))
      END WHERE
    END DO

    !> 6.1 water in the surface reservoir [m]
    !! 
    w_skin(:) = 0.0_wp

    !> 6.2 adjust initial water amount to prescribed theta 
    !! 
    w_soil_root(:)        = w_soil_root_pwp(:) + & 
                            dsl4jsb_Config(SOIL_)%soil_theta_prescribe * (w_soil_root_fc(:) - w_soil_root_pwp(:)) 
    w_soil_root_theta(:)  = dsl4jsb_Config(SOIL_)%soil_theta_prescribe
    w_soil_root_pot(:)    = saxtonA(:,1) * (MIN(w_soil_root_fc,w_soil_root(:)) / root_depth(:)) ** saxtonB(:,1) 

    !> 6.3 root fraction 
    !! 
    !! root fractions assuming an exponential distribution of roots \n
    !!   Cr_l = Cr_0 * exp(-lctlib%k_root_dist*depth) \n
    !! integral to bottom of root zone: 1/lctlib%k_root_dist * (1 - exp(-lctlib%k_root_dist*depth) \n
    !!   thus for each layer
    !!
    !! root_fraction_sl(:,:) is a variable of veg_mem and depends on the init of root_depth(:), see section 6.0 and ticket #358
    !! 
    ! layer 1
    root_fraction_sl(:,1) = (1._wp - exp(-lctlib%k_root_dist * dsl4jsb_Config(SOIL_)%soil_depth_water(1))) / &
                            (1._wp - exp(-lctlib%k_root_dist * root_depth(:)))
    ! layer 2 -> nsoil_sb
    DO isoil=2,nsoil_sb
       WHERE(root_depth(:) > dsl4jsb_Config(SOIL_)%soil_depth_water(isoil-1))
          WHERE(root_depth(:) > dsl4jsb_Config(SOIL_)%soil_depth_water(isoil))
             root_fraction_sl(:,isoil) = (exp(-lctlib%k_root_dist * dsl4jsb_Config(SOIL_)%soil_depth_water(isoil-1)) - &
                                          exp(-lctlib%k_root_dist * dsl4jsb_Config(SOIL_)%soil_depth_water(isoil))) / & 
                                         (1._wp - exp(-lctlib%k_root_dist * root_depth(:)))
          ELSEWHERE
             root_fraction_sl(:,isoil) = (exp(-lctlib%k_root_dist * dsl4jsb_Config(SOIL_)%soil_depth_water(isoil-1)) - &
                                          exp(-lctlib%k_root_dist * root_depth(:))) / & 
                                         (1._wp - exp(-lctlib%k_root_dist * root_depth(:)))
          ENDWHERE
       ELSEWHERE
          root_fraction_sl(:,isoil) = 0.0_wp 
       ENDWHERE
    ENDDO

    !> 7.0 soil water per layer depth
    !! 
    w_soil_pwp_sl(:,:)  = 0.0_wp
    w_soil_fc_sl(:,:)   = 0.0_wp
    w_soil_sat_sl(:,:)  = 0.0_wp
    DO isoil=1,nsoil_sb
       WHERE(soil_depth_sl(:,isoil) > 0.0_wp) 
          w_soil_pwp_sl(:,isoil) = theta_pwp_sl(:,isoil) * soil_depth_sl(:,isoil) 
          w_soil_fc_sl(:,isoil)  = theta_fc_sl(:,isoil) * soil_depth_sl(:,isoil) 
          w_soil_sat_sl(:,isoil) = theta_sat_sl(:,isoil) * soil_depth_sl(:,isoil) 
       ENDWHERE
    ENDDO
    
    !> 7.1 adjust initial water amount to prescribed theta
    !! 
    w_soil_sl(:,:) = w_soil_pwp_sl(:,:) + dsl4jsb_Config(SOIL_)%soil_theta_prescribe * (w_soil_fc_sl(:,:) - w_soil_pwp_sl(:,:))

    !> 7.2 soil water potential per layer
    !! 
    WHERE(soil_depth_sl(:,:) > eps8)
       hlp1(:,:)          = MIN(w_soil_fc_sl(:,:),w_soil_sl(:,:)) / soil_depth_sl(:,:)
       w_soil_pot_sl(:,:) = saxtonA(:,:) * hlp1(:,:) ** saxtonB(:,:)
    ELSEWHERE
       w_soil_pot_sl(:,:) = 0.0_wp 
    ENDWHERE

  END SUBROUTINE soil_init

END MODULE mo_soil_init
