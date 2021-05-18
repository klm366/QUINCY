!>
!!  @par Copyright
!!    QUINCY terrestrial biosphere model with coupled carbon, nitrogen, and phosphorus cycles
!!    Copyright (C) 2019 Max Planck Institute for Biogeochemistry (MPI-BGC), Jena, Germany 
!!    Contributors: Caldararu S., Engel J., Kern M., Pallandt M., Schnur R., Thum T., Yu L., Zaehle S.
!!
!!    This file is part of the QUINCY model.
!!    The QUINCY model is free software: you can redistribute it and/or modify
!!    it under the terms of the GNU General Public License version 3,
!!    as published by the Free Software Foundation.
!!
!!    The QUINCY model is distributed in the hope that it will be useful,
!!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!!    GNU General Public License for more details.
!!
!!    You should have received a copy of the GNU General Public License
!!    along with the QUINCY model.  If not, see <https://www.gnu.org/licenses/>.
!! 
!! @brief Initialization of the the soil_biogeochemistry memory.
!!
!! @author
!!  Sönke Zaehle, MPI-BGC Jena
!!
MODULE mo_sb_init

#ifdef __QUINCY_NO_JSB__
#include "dsl_jsbach.inc"
#endif

  USE mo_kind,                    ONLY: wp
  USE mo_exception,               ONLY: message, finish, message_text
  USE mo_jsb_control,             ONLY: debug_on

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: sb_init

  CHARACTER(len=*), PARAMETER :: modname = 'mo_sb_init'

CONTAINS

  ! ======================================================================================================= !
  !>Intialize SB_ process memory for QUINCY
  !>
  SUBROUTINE sb_init(tile)

    USE mo_jsb_class,                   ONLY: Get_model
    USE mo_jsb_tile_class,              ONLY: t_jsb_tile_abstract
    USE mo_jsb_model_class,             ONLY: t_jsb_model
    USE mo_jsb_lctlib_class,            ONLY: t_lctlib_element
    USE mo_jsb_grid_class,              ONLY: t_jsb_grid, t_jsb_vgrid
    USE mo_jsb_grid,                    ONLY: Get_grid, Get_vgrid
    USE mo_jsb_process_class,           ONLY: SB_, SSE_

    USE mo_jsb_math_constants,          ONLY: eps4, eps8
    USE mo_sb_constants,                ONLY: nh4_solute_prescribe, no3_solute_prescribe
    dsl4jsb_Use_config(SB_)
    dsl4jsb_Use_memory(SB_)
    dsl4jsb_Use_memory(SSE_)
#ifdef __QUINCY_NO_JSB__
    dsl4slm_Use_memory_a_sb
    USE mo_sb_memory_class,       ONLY: sb_pool
#else
    dsl4jsb_Use_pools(SB_) t_sb_components_pool
#endif
    ! ----------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_tile_abstract), INTENT(inout)     :: tile          !< one tile with data structure for one lct
    ! ----------------------------------------------------------------------------------------------------- !
    TYPE(t_jsb_model),      POINTER               :: model
    TYPE(t_lctlib_element), POINTER               :: lctlib        !< land-cover-type library - parameter across pft's
    TYPE(t_jsb_grid),       POINTER               :: hgrid         !< Horizontal grid
    TYPE(t_jsb_vgrid),      POINTER               :: vgrid_soil_sb !< Vertical grid
    INTEGER                                       :: nsoil_sb      !< number of soil layers as used/defined by the SB_ process
    REAL(wp)                                      :: hlp1
    REAL(wp), ALLOCATABLE, DIMENSION(:,:)         :: hlp2
    REAL(wp), ALLOCATABLE, DIMENSION(:,:,:)       :: hlp3
    INTEGER                                       :: isoil, ic, nc, iblk, nblks
    INTEGER                                       :: nproma, npromz
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':sb_init'
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Def_config(SB_)
    dsl4jsb_Def_memory(SB_)
    dsl4jsb_Def_memory(SSE_)
#ifdef __QUINCY_NO_JSB__
#else
    dsl4jsb_Def_pool(t_sb_components_pool)  :: sb_pool
#endif
    ! ----------------------------------------------------------------------------------------------------- !
    ! SB_ 3D
    dsl4jsb_Real3D_onDomain :: nh4_solute
    dsl4jsb_Real3D_onDomain :: no3_solute
    ! SSE_ 3D  TODO soil_lay_* var may be moved to HYDRO_
    dsl4jsb_Real3D_onDomain :: soil_lay_width_sl
    dsl4jsb_Real3D_onDomain :: soil_lay_depth_center_sl
    dsl4jsb_Real3D_onDomain :: soil_lay_depth_ubound_sl
    dsl4jsb_Real3D_onDomain :: soil_lay_depth_lbound_sl
    ! ----------------------------------------------------------------------------------------------------- !
    IF (.NOT. tile%Is_process_active(SB_)) RETURN
    IF (debug_on()) CALL message(TRIM(routine), 'Initializing QUINCY SB_ memory for tile '//TRIM(tile%name))
    ! ----------------------------------------------------------------------------------------------------- !
    model         => Get_model(tile%owner_model_id)
    IF (.NOT. model%config%use_quincy) RETURN       ! only with QUINCY model
    IF (tile%lcts(1)%lib_id == 0) RETURN            ! run this init only if the present tile is a pft
    lctlib        => model%lctlib(tile%lcts(1)%lib_id)
    hgrid         => Get_grid(model%grid_id)
    vgrid_soil_sb => Get_vgrid('soil_layer_sb')
    nsoil_sb      =  vgrid_soil_sb%n_levels
    nproma        =  hgrid%nproma
    npromz        =  hgrid%npromz
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Get_config(SB_)
    dsl4jsb_Get_memory(SB_)
    dsl4jsb_Get_memory(SSE_)
#ifdef __QUINCY_NO_JSB__
#else
    ! Get pools
    dsl4jsb_Get_pool(SB_, sb_main_pool%sb_pool, sb_pool)
    ! associate domain variables
    CALL sb_pool%Associate_domain_var_ptr()
#endif
    ! ----------------------------------------------------------------------------------------------------- !
    ! SB_ 3D
    dsl4jsb_Get_var3D_onDomain(SB_, nh4_solute)
    dsl4jsb_Get_var3D_onDomain(SB_, no3_solute)
    ! SSE_ 3D  TODO soil_lay_* var might be moved to HYDRO_ ?
    dsl4jsb_Get_var3D_onDomain(SSE_, soil_lay_width_sl)       
    dsl4jsb_Get_var3D_onDomain(SSE_, soil_lay_depth_center_sl)
    dsl4jsb_Get_var3D_onDomain(SSE_, soil_lay_depth_ubound_sl)
    dsl4jsb_Get_var3D_onDomain(SSE_, soil_lay_depth_lbound_sl)
    ! ----------------------------------------------------------------------------------------------------- !
    nc    = SIZE(nh4_solute,1)
    nblks = SIZE(nh4_solute,3)
    ! ----------------------------------------------------------------------------------------------------- !
    ALLOCATE(hlp2(nc, nblks))
    ALLOCATE(hlp3(nc, nsoil_sb, nblks))
    ! ----------------------------------------------------------------------------------------------------- !
#ifdef __QUINCY_NO_JSB__
#define Get_sb_pool
#include "slm_pools_a_incl.f90"
#undef Get_sb_pool
#endif

    ! -------------------------------------------------------------------------------------------------------
    ! Go science
    ! -------------------------------------------------------------------------------------------------------

    ! TODO - remove the below test-code for pools when the QUINCY land model is implemented
    sb_pool                         = 1.2_wp
    sb_pool                         = sb_pool + sb_pool * 3.0_wp
    sb_pool%soluable_litter         = sb_pool%polymeric_litter + sb_pool%microbial
    sb_pool%residue%domain_c(:,:,:) = sb_pool%residue%domain_c(:,:,:) * sb_pool%residue%domain_n(:,:,:) + 2.0_wp
    !

    ! ----------------------------------------------------------------------------------------------------- !
    !>1.0 
    !>
    !>  1.1 simple_0d
    !>
    ! nh4_solute(:,1,:) = nh4_solute_prescribe / soil_lay_width_sl(:,1,:)
    ! no3_solute(:,1,:) = no3_solute_prescribe / soil_lay_width_sl(:,1,:)

    !>  1.2 simple_1d
    !>
    ! note: this differs on purpose from the calculation in mo_sb_update_pools for the prescribed option
    !       this is because here I'm trying to ensure that the models are initialised with the same amount of nutrient
    !       whereas in the prescribed option the target is unlimiting concentrations across the profile. The latter leads
    !       to more nutrient in the soil column in 1d compared to 0d.
    DO iblk = 1,nblks
      ! set nc correctly avoid division by zero (soil_lay_depth* at nblks might be initialized with zero for values larger npromz)
      IF (iblk < nblks) THEN
        nc = nproma
      ELSE
        nc = npromz
      ENDIF
      DO ic = 1,nc
        nh4_solute(ic,:,iblk) = nh4_solute_prescribe / soil_lay_depth_ubound_sl(ic,nsoil_sb,iblk)
        no3_solute(ic,:,iblk) = no3_solute_prescribe / soil_lay_depth_ubound_sl(ic,nsoil_sb,iblk)
      ENDDO
    ENDDO

    
  END SUBROUTINE sb_init

END MODULE mo_sb_init
