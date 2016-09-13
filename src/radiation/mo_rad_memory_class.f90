!>
!! @brief Contains the memory class for the radiation process.
!!
!! @author
!!  Jan Engel, MPI-BGC Jena (based on the JSBACH4 version from Rainer Schneck, MPI-M Hamburg)
!!
!! @par Revision History
!! First version by Jan Engel (20160819)
!!
MODULE mo_rad_memory_class

  USE mo_kind, ONLY: wp

  USE mo_jsb_memory_class,       ONLY: t_jsb_memory
  USE mo_jsb_var_class,          ONLY: t_jsb_var_real2d, t_jsb_var_real3d, t_jsb_var_real1d

  USE mo_slm_varlist,            ONLY: slm_add_var   ! replace by mo_jsb_varlist later   *slmtemp



  IMPLICIT NONE
  PRIVATE
  PUBLIC :: t_rad_memory, max_no_of_vars

  INTEGER, PARAMETER :: max_no_of_vars = 40   ! default JSBACH4

  !-----------------------------------------------------------------------------------------------------
  !> Type definition for rad memory
  !! 
  !-----------------------------------------------------------------------------------------------------
  TYPE, EXTENDS(t_jsb_memory) :: t_rad_memory

! *slmtemp 1D -> 2D
    TYPE(t_jsb_var_real1d)  :: &
      & sw_srf_net,            & !< net shortwave radiation at surface
      & swvis_srf_net,         & !< net shortwave radiation at surface in visible range
      & swnir_srf_net,         & !< net shortwave radiation at surface in NIR range
      & lw_srf_net,            & !< net longwave radiation at surface
      & rad_srf_net              !< net radiation at surface

    TYPE(t_jsb_var_real1d)  :: &        ! REAL(wp) 1D ; implemented for SLM use only, see mo_jsb_var_class 
      & alb_vis,               & !< Albedo of the whole tile in the VIS (visible) range
      & alb_nir,               & !< Albedo of the whole tile in the NIR (near-infrared) range
      & alb_vis_snow,          & !< Albedo of snow in the VIS range
      & alb_nir_snow             !< Albedo of snow in the NIR range                               

    ! Additional variables for soil
    TYPE(t_jsb_var_real1d) :: &
      & alb_vis_soil,         & !< Albedo in VIS range of soil
                                ! Initialized from default:
      & alb_nir_soil,         & !< Albedo in VIS range of soil
      & arad_vis_soil,        & !< radiation in VIS range absorbed by the soil
      & arad_nir_soil           !< radiation in NIR range absorbed by the soil

    ! Additional variables for vegetation
    TYPE(t_jsb_var_real1d) :: &
      & alb_vis_can,          & !< Albedo of canopy (without snow cover) in VIS range.
      & alb_nir_can,          & !< Albedo of canopy (without snow cover) in NIR range. 
                                !< If use_alb_canopy=.TRUE. taken from bc_land_phys.nc file, else from lctlib-File (AlbedoCanopyNIR)
      & arad_vis_can,         & !< radiation in VIS range absorbed by the canopy      
      & arad_nir_can            !< radiation in NIR range absorbed by the canopy

    ! misses par,arad by canopy layer

  CONTAINS
    PROCEDURE :: Init => Init_rad_memory
  END TYPE t_rad_memory


  CHARACTER(len=*), PARAMETER :: modname = 'mo_rad_memory_class'




CONTAINS


  !-----------------------------------------------------------------------------------------------------
  !> initialize memory (variables) for the process: radiation
  !!
  !! this is a radically simplified version of the init routine when compared to JSBACH4 (20160819) \n
  !! please see Init_assimi_memory where the missign code is just commeted out (i.e. not deleted)
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE Init_rad_memory(mem, prefix, suffix, lct_ids, model_id, ics, ice, iblk)

    CLASS(t_rad_memory), INTENT(inout), TARGET :: mem   !> radiation memory
    CHARACTER(len=*),    INTENT(in)    :: prefix        !> ?
    CHARACTER(len=*),    INTENT(in)    :: suffix        !> ?
    INTEGER,             INTENT(in)    :: lct_ids!(:)   !> in jsbach4: vector of land cover type IDs
    INTEGER,             INTENT(in)    :: model_id      !> model ID
    INTEGER,             INTENT(in)    :: ics           !> model%options -> for SLM compatibility
    INTEGER,             INTENT(in)    :: ice           !> model%options -> for SLM compatibility
    INTEGER,             INTENT(in)    :: iblk          !> model%options -> for SLM compatibility


    IF (.NOT. ASSOCIATED(mem%sw_srf_net%ptr)) THEN

      CALL slm_add_var( mem, 'sw_srf_net', mem%sw_srf_net, &
            ics, ice, iblk, &   ! model%options -> for SLM compatibility
            initval_r=0.0_wp)   ! initval taken from JSBACH3

      CALL slm_add_var( mem, 'swvis_srf_net', mem%swvis_srf_net, &
            ics, ice, iblk, &   ! model%options -> for SLM compatibility
            initval_r=0.0_wp)   ! initval taken from JSBACH3

      CALL slm_add_var( mem, 'swnir_srf_net', mem%swnir_srf_net, &
            ics, ice, iblk, &   ! model%options -> for SLM compatibility
            initval_r=0.0_wp)   ! initval taken from JSBACH3

      CALL slm_add_var( mem, 'lw_srf_net', mem%lw_srf_net, &
            ics, ice, iblk, &   ! model%options -> for SLM compatibility
            initval_r=0.0_wp)   ! initval taken from JSBACH3

      CALL slm_add_var( mem, 'rad_srf_net', mem%rad_srf_net, &
            ics, ice, iblk, &   ! model%options -> for SLM compatibility
            initval_r=0.0_wp)   ! initval taken from JSBACH3

      CALL slm_add_var(mem, 'alb_vis', mem%alb_vis, &   
            ics, ice, iblk, &   ! model%options -> for SLM compatibility
            initval_r=0.0_wp)   ! initval taken from JSBACH3

      CALL slm_add_var(mem, 'alb_nir', mem%alb_nir, &   
            ics, ice, iblk, &   ! model%options -> for SLM compatibility
            initval_r=0.0_wp)   ! initval taken from JSBACH3

      CALL slm_add_var(mem, 'alb_vis_snow', mem%alb_vis_snow, &   
            ics, ice, iblk, &   ! model%options -> for SLM compatibility
            initval_r=0.0_wp)   ! initval taken from JSBACH3

      CALL slm_add_var(mem, 'alb_nir_snow', mem%alb_nir_snow, &   
            ics, ice, iblk, &   ! model%options -> for SLM compatibility
            initval_r=0.0_wp)   ! initval taken from JSBACH3

      CALL slm_add_var(mem, 'alb_vis_soil', mem%alb_vis_soil, &   
            ics, ice, iblk, &   ! model%options -> for SLM compatibility
            initval_r=0.0_wp)   ! initval taken from JSBACH3

      CALL slm_add_var(mem, 'alb_nir_soil', mem%alb_nir_soil, &   
            ics, ice, iblk, &   ! model%options -> for SLM compatibility
            initval_r=0.0_wp)   ! initval taken from JSBACH3

      CALL slm_add_var(mem, 'arad_vis_soil', mem%arad_vis_soil, &   
            ics, ice, iblk, &   ! model%options -> for SLM compatibility
            initval_r=0.0_wp)   ! initval taken from JSBACH3

      CALL slm_add_var(mem, 'arad_nir_soil', mem%arad_nir_soil, &   
            ics, ice, iblk, &   ! model%options -> for SLM compatibility
            initval_r=0.0_wp)   ! initval taken from JSBACH3

      CALL slm_add_var(mem, 'alb_vis_can', mem%alb_vis_can, &   
            ics, ice, iblk, &   ! model%options -> for SLM compatibility
            initval_r=0.0_wp)   ! initval taken from JSBACH3

      CALL slm_add_var(mem, 'alb_nir_can', mem%alb_nir_can, &   
            ics, ice, iblk, &   ! model%options -> for SLM compatibility
            initval_r=0.0_wp)   ! initval taken from JSBACH3

      CALL slm_add_var(mem, 'arad_vis_can', mem%arad_vis_can, &   
            ics, ice, iblk, &   ! model%options -> for SLM compatibility
            initval_r=0.0_wp)   ! initval taken from JSBACH3

      CALL slm_add_var(mem, 'arad_nir_can', mem%arad_nir_can, &   
            ics, ice, iblk, &   ! model%options -> for SLM compatibility
            initval_r=0.0_wp)   ! initval taken from JSBACH3

    END IF

  END SUBROUTINE Init_rad_memory

END MODULE mo_rad_memory_class
