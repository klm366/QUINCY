!>
!! @brief Contains the memory class for the assimilation process.
!!
!! @author
!!  Jan Engel, MPI-BGC Jena (based on the JSBACH4 version from Rainer Schneck, MPI-M Hamburg)
!!
!! @par Revision History
!! First version by Jan Engel (20160815)
!!
MODULE mo_assimi_memory_class

  USE mo_kind,                   ONLY: wp
  USE mo_exception,              ONLY: message
  USE mo_jsb_class,              ONLY: Get_model

  USE mo_jsb_memory_class,       ONLY: t_jsb_memory
  USE mo_jsb_var_class,          ONLY: t_jsb_var_real2d, t_jsb_var_real3d

  USE mo_slm_varlist,            ONLY: slm_add_var   ! replace by mo_jsb_varlist later  *slmtemp



  IMPLICIT NONE
  PRIVATE
  PUBLIC :: t_assimi_memory, max_no_of_vars

  INTEGER, PARAMETER :: max_no_of_vars = 40   ! default JSBACH4


  !-----------------------------------------------------------------------------------------------------
  !> Type definition for assimi memory
  !! 
  !-----------------------------------------------------------------------------------------------------
  TYPE, EXTENDS(t_jsb_memory) :: t_assimi_memory

  ! canopy integrated variables
    TYPE(t_jsb_var_real2d)    :: &
                              gross_assimilation    , &      !< canopy-integrated gross assimilation rate (µmol CO2/m2/s) 
                              gross_assimilation_C13, &      !< canopy-integrated gross assimilation rate (µmol 13CO2/m2/s) 
                              gross_assimilation_C14, &      !< canopy-integrated gross assimilation rate (µmol 14CO2/m2/s) 
                              net_assimilation      , &      !< canopy-integrated net assimilation rate (µmol CO2/m2/s) 
                              net_assimilation_boc  , &      !< net assimilation rate at the bottom of the canopy (µmol CO2/m2/s) 
                              maint_respiration_leaf, &      !< canopy-integrated foliar respiration rate (µmol CO2/m2/s)
                              canopy_cond           , &      !< canopy-integrated conductance for water (m/s)
                              co2_conc_leaf         , &      !< canopy-integrated internal leaf CO2 concentration (ppm) 
                              beta_air              , &      !< scaling factor to account for air humidity effects on conductance
                              beta_soil_ps          , &      !< scaling factor to account for soil moisture constraints on photosynthesis
                              beta_soil_gs          , &      !< scaling factor to account for soil moisture constraints on conductance
                              t_jmax_opt                     !< optimal temperature for electron transport 

    TYPE(t_jsb_var_real3d)    :: &
                              gross_assimilation_cl     , &   !< gross assimilation rate per canopy layer (µmol CO2/m2/s) 
                              net_assimilation_cl       , &   !< net assimilation rate per canopy layer (µmol CO2/m2/s)
                              maint_respiration_leaf_cl , &   !< foliar respiration rate per canopy layer (µmol CO2/m2/s)
                              canopy_cond_cl            , &   !< conductance for water per canopy layer (m/s)
                              co2_conc_leaf_cl                !< internal leaf CO2 concentration per canopy layer (ppm)

    TYPE(t_jsb_var_real2d)    :: &
                              beta_air_daytime          , &   !< average of the previous day
                              beta_air_daytime_dacc     , &   !< average of the previous day
                              beta_air_tfrac_mavg       , &   !< scaling factor to account for air humidity effects on conductance
                              beta_air_tcnl_mavg              !< scaling factor to account for air humidity effects on conductance

    TYPE(t_jsb_var_real2d)    :: &
                              beta_soil_ps_daytime      , &   !< average of the previous day
                              beta_soil_ps_daytime_dacc , &   !< average of the previous day
                              beta_soil_ps_tfrac_mavg   , &   !< scaling factor to account for soil moisture constraints on photosynthesis 
                              beta_soil_ps_tcnl_mavg    , &   !< scaling factor to account for soil moisture constraints on conductance
                              beta_soil_gs_daytime      , &   !< average of the previous day
                              beta_soil_gs_daytime_dacc , &   !< average of the previous day
                              beta_soil_gs_tphen_mavg   , &   !< scaling factor to account for soil moisture constraints on conductance
                              beta_soil_gs_tfrac_mavg   , &   !< scaling factor to account for soil moisture constraints on conductance
                              beta_soil_gs_tcnl_mavg          !< scaling factor to account for soil moisture constraints on conductance




  CONTAINS
    PROCEDURE :: Init => Init_assimi_memory
  END TYPE t_assimi_memory


  CHARACTER(len=*), PARAMETER :: modname = 'mo_assimi_memory_class'


CONTAINS


  !-----------------------------------------------------------------------------------------------------
  !> initialize memory (variables) for the process: assimilation
  !!
  !!  
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE Init_assimi_memory(mem, prefix, suffix, lct_ids, model_id)

!    USE mo_jsb_varlist,       ONLY: t_jsb_varlist, jsb_add_var
!    USE mo_jsb_io,            ONLY: grib_bits, t_cf, t_grib1, t_grib2, TSTEP_CONSTANT, tables
!    USE mo_jsb_grid_class,    ONLY: t_jsb_grid, t_jsb_vgrid
!    USE mo_jsb_grid,          ONLY: Get_grid, Get_vgrid
    USE mo_base_memory,         ONLY: t_model

    ! 0.1.InOut
    CLASS(t_assimi_memory), INTENT(inout), TARGET :: mem  !> assimi memory
    CHARACTER(len=*),    INTENT(in)    :: prefix          !> ?
    CHARACTER(len=*),    INTENT(in)    :: suffix          !> ?
    INTEGER,             INTENT(in)    :: lct_ids !(:)    !> in jsbach4: vector of land cover type IDs
    INTEGER,             INTENT(in)    :: model_id        !> model ID model\%id
    ! 0.2 Local
    TYPE(t_model), POINTER :: model
!    TYPE(t_jsb_model), POINTER :: model
!    TYPE(t_jsb_grid),  POINTER :: hgrid                        ! Horizontal grid
!    TYPE(t_jsb_vgrid), POINTER :: surface, canopy_layer        ! Vertical grid
!    INTEGER :: table

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':Init_assimi_memory'

    model => Get_model(model_id)

!    table = tables(1)

!    hgrid   => Get_grid(mem%grid_id)
!    surface => Get_vgrid('surface')
!    canopy_layer => Get_vgrid('canopy_layer')

    IF (.NOT. ALLOCATED(mem%gross_assimilation%ptr)) THEN

      !> add variables to assimi memory assimi\%VAR
      !! 
      !! slm_add_var -> slm_add_var_list_element_r1d    *slmtemp 1D -> 2D \n
      !! and \n
      !! slm_add_var -> slm_add_var_list_element_r1d    *slmtemp 1D -> 3D
      !!  
      CALL slm_add_var(mem, 'gross_assimilation', mem%gross_assimilation, &
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'gross_assimilation_C13', mem%gross_assimilation_C13, &
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'gross_assimilation_C14', mem%gross_assimilation_C14, &
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'net_assimilation', mem%net_assimilation, &
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'net_assimilation_boc', mem%net_assimilation_boc, &
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'maint_respiration_leaf', mem%maint_respiration_leaf, &
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'canopy_cond', mem%canopy_cond, &
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'co2_conc_leaf', mem%co2_conc_leaf, &
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'beta_air', mem%beta_air, &
            initval_r=1.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'beta_soil_ps', mem%beta_soil_ps, &
            initval_r=1.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'beta_soil_gs', mem%beta_soil_gs, &
            initval_r=1.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 't_jmax_opt', mem%t_jmax_opt, &
            initval_r=1.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'gross_assimilation_cl', mem%gross_assimilation_cl, &
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'net_assimilation_cl', mem%net_assimilation_cl, &
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'maint_respiration_leaf_cl', mem%maint_respiration_leaf_cl, &
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'canopy_cond_cl', mem%canopy_cond_cl, &
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'co2_conc_leaf_cl', mem%co2_conc_leaf_cl, &
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'beta_air_daytime', mem%beta_air_daytime, &
            initval_r=1.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'beta_air_daytime_dacc', mem%beta_air_daytime_dacc, &
            initval_r=1.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'beta_air_tfrac_mavg', mem%beta_air_tfrac_mavg, &
            initval_r=1.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'beta_air_tcnl_mavg', mem%beta_air_tcnl_mavg, &
            initval_r=1.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'beta_soil_ps_daytime', mem%beta_soil_ps_daytime, &
            initval_r=1.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'beta_soil_ps_daytime_dacc', mem%beta_soil_ps_daytime_dacc, &
            initval_r=1.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'beta_soil_ps_tfrac_mavg', mem%beta_soil_ps_tfrac_mavg, &
            initval_r=1.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'beta_soil_ps_tcnl_mavg', mem%beta_soil_ps_tcnl_mavg, &
            initval_r=1.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'beta_soil_gs_daytime', mem%beta_soil_gs_daytime, &
            initval_r=1.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'beta_soil_gs_daytime_dacc', mem%beta_soil_gs_daytime_dacc, &
            initval_r=1.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'beta_soil_gs_tphen_mavg', mem%beta_soil_gs_tphen_mavg, &
            initval_r=1.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'beta_soil_gs_tfrac_mavg', mem%beta_soil_gs_tfrac_mavg, &
            initval_r=1.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'beta_soil_gs_tcnl_mavg', mem%beta_soil_gs_tcnl_mavg, &
            initval_r=1.0_wp) ! initvalue

    END IF

  END SUBROUTINE Init_assimi_memory

END MODULE mo_assimi_memory_class
