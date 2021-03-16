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
!! @brief Contains the utilities for the soil_biogeochemistry process
!!
!! @author
!!  Jan Engel, MPI-BGC Jena
!!
!! @par Revision History
!! First version by Jan Engel (20170104)
!!
MODULE mo_sb_util

  USE mo_kind,                        ONLY: wp
  USE mo_jsb_control,                 ONLY: debug_on
  USE mo_exception,                   ONLY: message, finish

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: reset_sb_fluxes, update_time_average_soilbiogeochemistry

  CHARACTER(len=*), PARAMETER, PRIVATE :: modname = 'mo_sb_util'


CONTAINS


  ! ======================================================================================================= !
  !>init/reset sb fluxes (with/to zero)
  !> 
  SUBROUTINE reset_sb_fluxes(tile, options)

    USE mo_jsb_class,             ONLY: Get_model
    USE mo_jsb_tile_class,        ONLY: t_jsb_tile_abstract
    USE mo_jsb_task_class,        ONLY: t_jsb_task_options
    USE mo_jsb_model_class,       ONLY: t_jsb_model
    USE mo_jsb_math_constants,    ONLY: zero
    dsl4jsb_Use_processes SB_
    dsl4jsb_Use_memory(SB_)

    CLASS(t_jsb_tile_abstract), INTENT(inout)     :: tile         !< one tile with data structure for one lct
    TYPE(t_jsb_task_options),   INTENT(in)        :: options      !< model options

    TYPE(t_jsb_model),      POINTER       :: model
    INTEGER                               :: iblk, ics, ice, nc
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':reset_sb_fluxes'

    dsl4jsb_Def_memory(SB_)                               ! TYPE(t_sb_memory),       POINTER :: sb__mem

    ! 2D SB_
    dsl4jsb_Real2D_onChunk :: total_flux_carbon
    dsl4jsb_Real2D_onChunk :: total_flux_nitrogen
    dsl4jsb_Real2D_onChunk :: total_flux_phosphorus
    dsl4jsb_Real2D_onChunk :: total_flux_carbon13
    dsl4jsb_Real2D_onChunk :: total_flux_carbon14
    dsl4jsb_Real2D_onChunk :: total_flux_nitrogen15

    iblk  = options%iblk
    ics   = options%ics
    ice   = options%ice
    nc    = options%nc

    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    dsl4jsb_Get_memory(SB_)

    model  => Get_model(tile%owner_model_id)

    dsl4jsb_Get_var2D_onChunk(SB_,       total_flux_carbon)           ! out
    dsl4jsb_Get_var2D_onChunk(SB_,       total_flux_nitrogen)         ! out
    dsl4jsb_Get_var2D_onChunk(SB_,       total_flux_phosphorus)       ! out
    dsl4jsb_Get_var2D_onChunk(SB_,       total_flux_carbon13)         ! out
    dsl4jsb_Get_var2D_onChunk(SB_,       total_flux_carbon14)         ! out
    dsl4jsb_Get_var2D_onChunk(SB_,       total_flux_nitrogen15)       ! out

    ! ------------------------------------------------------------------------------------------------------------
    ! Go science
    ! ------------------------------------------------------------------------------------------------------------

    !> 1.0 initialise/reset fluxes
    !!
    ! 2D fields  
    total_flux_carbon(:)         = zero
    total_flux_nitrogen(:)       = zero
    total_flux_phosphorus(:)     = zero
    total_flux_carbon13(:)       = zero
    total_flux_carbon14(:)       = zero
    total_flux_nitrogen15(:)     = zero

  END SUBROUTINE reset_sb_fluxes


  ! ======================================================================================================= !
  !> calculate time moving averages and daytime averages for SB_
  !>
  SUBROUTINE update_time_average_soilbiogeochemistry(tile, options)

    USE mo_jsb_control,             ONLY: debug_on
    USE mo_jsb_tile_class,          ONLY: t_jsb_tile_abstract
    USE mo_jsb_task_class,          ONLY: t_jsb_task_options
    USE mo_jsb_lctlib_class,        ONLY: t_lctlib_element
    USE mo_jsb_model_class,         ONLY: t_jsb_model
    USE mo_jsb_class,               ONLY: Get_model
    USE mo_time_averages            ! e.g. calc_time_mavg, mavg_period_tphen, mavg_period_weekly
    dsl4jsb_Use_processes SB_
    dsl4jsb_Use_memory(SB_)
#ifdef __QUINCY_NO_JSB__
    dsl4slm_Use_memory_a_sb
    USE mo_sb_memory_class,    ONLY: sb_mycorrhiza_export
#else
    dsl4jsb_Use_pools(SB_) t_sb_elements_pool
#endif
    ! ----------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options
    ! ----------------------------------------------------------------------------------------------------- !
    TYPE(t_jsb_model),      POINTER         :: model                !< instance of the model
    TYPE(t_lctlib_element), POINTER         :: lctlib               !< land-cover-type library - parameter across pft's
    INTEGER                                 :: iblk, ics, ice, nc   !< dimensions
    CHARACTER(len=*), PARAMETER             :: routine = modname//':update_time_average_soilbiogeochemistry'
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Def_memory(SB_)
#ifdef __QUINCY_NO_JSB__
#else
    dsl4jsb_Def_pool(t_sb_elements_pool)    :: sb_mycorrhiza_export
#endif
    ! ----------------------------------------------------------------------------------------------------- !
    ! 3D SB_
    dsl4jsb_Real3D_onChunk      :: microbial_cue_eff
    dsl4jsb_Real3D_onChunk      :: microbial_nue_eff
    dsl4jsb_Real3D_onChunk      :: microbial_pue_eff
    dsl4jsb_Real3D_onChunk      :: dom_cn
    dsl4jsb_Real3D_onChunk      :: dom_cp
    dsl4jsb_Real3D_onChunk      :: enzyme_frac_poly_c
    dsl4jsb_Real3D_onChunk      :: enzyme_frac_poly_n
    dsl4jsb_Real3D_onChunk      :: enzyme_frac_poly_p
    dsl4jsb_Real3D_onChunk      :: microbial_cue_eff_tmic_mavg
    dsl4jsb_Real3D_onChunk      :: microbial_nue_eff_tmic_mavg
    dsl4jsb_Real3D_onChunk      :: microbial_pue_eff_tmic_mavg
    dsl4jsb_Real3D_onChunk      :: myc_export_c_tmyc_mavg_sl
    dsl4jsb_Real3D_onChunk      :: myc_export_n_tmyc_mavg_sl
    dsl4jsb_Real3D_onChunk      :: myc_export_n_tlabile_mavg_sl
    dsl4jsb_Real3D_onChunk      :: myc_export_p_tlabile_mavg_sl
    dsl4jsb_Real3D_onChunk      :: dom_cn_mavg
    dsl4jsb_Real3D_onChunk      :: dom_cp_mavg
    dsl4jsb_Real3D_onChunk      :: enzyme_frac_poly_c_mavg
    dsl4jsb_Real3D_onChunk      :: enzyme_frac_poly_n_mavg
    dsl4jsb_Real3D_onChunk      :: enzyme_frac_poly_p_mavg
    ! ----------------------------------------------------------------------------------------------------- !
    iblk    = options%iblk
    ics     = options%ics
    ice     = options%ice
    nc      = options%nc
    ! ----------------------------------------------------------------------------------------------------- !
    IF (.NOT. tile%Is_process_active(SB_)) RETURN
    IF (tile%lcts(1)%lib_id == 0)          RETURN  ! only if the present tile is a pft
    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')
    ! ----------------------------------------------------------------------------------------------------- !
    model         => Get_model(tile%owner_model_id)
    lctlib        => model%lctlib(tile%lcts(1)%lib_id)
    dsl4jsb_Get_memory(SB_)
#ifdef __QUINCY_NO_JSB__
#else
    dsl4jsb_Get_pool(SB_, sb_main_pool%sb_mycorrhiza_export,    sb_mycorrhiza_export)
#endif
    ! ----------------------------------------------------------------------------------------------------- !
    ! 3D SB_
    dsl4jsb_Get_var3D_onChunk(SB_,  microbial_cue_eff)                  ! in
    dsl4jsb_Get_var3D_onChunk(SB_,  microbial_nue_eff)                  ! in
    dsl4jsb_Get_var3D_onChunk(SB_,  microbial_pue_eff)                  ! in
    dsl4jsb_Get_var3D_onChunk(SB_,  dom_cn)                             ! in
    dsl4jsb_Get_var3D_onChunk(SB_,  dom_cp)                             ! in
    dsl4jsb_Get_var3D_onChunk(SB_,  enzyme_frac_poly_c)                 ! in
    dsl4jsb_Get_var3D_onChunk(SB_,  enzyme_frac_poly_n)                 ! in
    dsl4jsb_Get_var3D_onChunk(SB_,  enzyme_frac_poly_p)                 ! in
    dsl4jsb_Get_var3D_onChunk(SB_,  microbial_cue_eff_tmic_mavg)        ! out
    dsl4jsb_Get_var3D_onChunk(SB_,  microbial_nue_eff_tmic_mavg)        ! out
    dsl4jsb_Get_var3D_onChunk(SB_,  microbial_pue_eff_tmic_mavg)        ! out
    dsl4jsb_Get_var3D_onChunk(SB_,  myc_export_c_tmyc_mavg_sl)          ! out
    dsl4jsb_Get_var3D_onChunk(SB_,  myc_export_n_tmyc_mavg_sl)          ! out
    dsl4jsb_Get_var3D_onChunk(SB_,  myc_export_n_tlabile_mavg_sl)       ! out
    dsl4jsb_Get_var3D_onChunk(SB_,  myc_export_p_tlabile_mavg_sl)       ! out
    dsl4jsb_Get_var3D_onChunk(SB_,  dom_cn_mavg)                        ! out
    dsl4jsb_Get_var3D_onChunk(SB_,  dom_cp_mavg)                        ! out
    dsl4jsb_Get_var3D_onChunk(SB_,  enzyme_frac_poly_c_mavg)            ! out
    dsl4jsb_Get_var3D_onChunk(SB_,  enzyme_frac_poly_n_mavg)            ! out
    dsl4jsb_Get_var3D_onChunk(SB_,  enzyme_frac_poly_p_mavg)            ! out
    ! ----------------------------------------------------------------------------------------------------- !
#ifdef __QUINCY_NO_JSB__
#define Get_sb_mycorrhiza_export
#include "slm_pools_a_incl.f90"
#undef Get_sb_mycorrhiza_export
#endif


    ! ----------------------------------------------------------------------------------------------------- !
    !> daytime averages - not used here
    !>

    ! ----------------------------------------------------------------------------------------------------- !
    !>1.0 moving averages
    !>
    !>

    ! docu:
    ! calc_time_mavg(current average, new value, length of avg_period,  !
    !                do_calc=LOGICAL, avg_period_unit='day')            ! OPTIONAL
    !                RETURN(new current average)
    ! the unit of the averaging period is 'day' by default, but can also be 'week' or 'year'
    
    !>  1.1 tlabile (averages at the timescale of the labile pool)
    !> 
    myc_export_n_tlabile_mavg_sl(:,:) = calc_mavg(myc_export_n_tlabile_mavg_sl(:,:), &
                                                  sb_mycorrhiza_export%nitrogen(:,:), mavg_period_tlabile) 
    myc_export_p_tlabile_mavg_sl(:,:) = calc_mavg(myc_export_p_tlabile_mavg_sl(:,:), &
                                                  sb_mycorrhiza_export%phosphorus(:,:), mavg_period_tlabile)

    !>  1.2 lctlib tau_mycorrhiza
    !>
    myc_export_c_tmyc_mavg_sl(:,:) = calc_mavg(myc_export_c_tmyc_mavg_sl(:,:), &
                                               sb_mycorrhiza_export%carbon(:,:), & 
 	                                             lctlib%tau_mycorrhiza, avg_period_unit='year')
    myc_export_n_tmyc_mavg_sl(:,:) = calc_mavg(myc_export_n_tmyc_mavg_sl(:,:), &
                                               sb_mycorrhiza_export%nitrogen(:,:), & 
 	                                             lctlib%tau_mycorrhiza, avg_period_unit='year')

    !>  1.3 tmic (microbial community acclimation)
    !> 
    microbial_cue_eff_tmic_mavg(:,:) = calc_mavg(microbial_cue_eff_tmic_mavg(:,:), microbial_cue_eff(:,:), &
                                                 mavg_period_tmic)
    microbial_nue_eff_tmic_mavg(:,:) = calc_mavg(microbial_nue_eff_tmic_mavg(:,:), microbial_nue_eff(:,:), &
                                                 mavg_period_tmic)
    microbial_pue_eff_tmic_mavg(:,:) = calc_mavg(microbial_pue_eff_tmic_mavg(:,:), microbial_pue_eff(:,:), &
                                                 mavg_period_tmic)

    !>  1.4 tenzyme (memory time-scale for enzyme allocation)
    !> 
    dom_cn_mavg(:,:)             = calc_mavg(dom_cn_mavg(:,:), dom_cn(:,:), &
                                             mavg_period_tenzyme)
    dom_cp_mavg(:,:)             = calc_mavg(dom_cp_mavg(:,:), dom_cp(:,:), &
                                             mavg_period_tenzyme)
    enzyme_frac_poly_c_mavg(:,:) = calc_mavg(enzyme_frac_poly_c_mavg(:,:), enzyme_frac_poly_c(:,:), &
                                             mavg_period_tenzyme)
    enzyme_frac_poly_n_mavg(:,:) = calc_mavg(enzyme_frac_poly_n_mavg(:,:), enzyme_frac_poly_n(:,:), &
                                             mavg_period_tenzyme)
    enzyme_frac_poly_p_mavg(:,:) = calc_mavg(enzyme_frac_poly_p_mavg(:,:), enzyme_frac_poly_p(:,:), &
                                             mavg_period_tenzyme)

! not needed, because veg_pool not modified:
! #ifdef __QUINCY_NO_JSB__
! #define Passback_sb_mycorrhiza_export
! #include "slm_pools_a_incl.f90"
! #undef Passback_sb_mycorrhiza_export
! #endif

  END SUBROUTINE update_time_average_soilbiogeochemistry


!#endif
END MODULE mo_sb_util
