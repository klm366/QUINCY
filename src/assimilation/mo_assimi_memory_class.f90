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
#ifndef __NO_JSBACH__

  USE mo_kind,                   ONLY: wp
  USE mo_exception,              ONLY: message
  USE mo_jsb_class,              ONLY: Get_model
  USE mo_jsb_memory_class,       ONLY: t_jsb_memory
  USE mo_jsb_var_class,          ONLY: t_jsb_var_real2d, t_jsb_var_real3d


  IMPLICIT NONE
  PRIVATE
  PUBLIC :: t_assimi_memory, max_no_of_vars

  INTEGER, PARAMETER :: max_no_of_vars = 60 

  !-----------------------------------------------------------------------------------------------------
  !> Type definition for assimi memory
  !! 
  !-----------------------------------------------------------------------------------------------------
  TYPE, EXTENDS(t_jsb_memory) :: t_assimi_memory

    ! fluxes
    TYPE(t_jsb_var_real2d)    :: &
                              gross_assimilation        , &  !< canopy-integrated gross assimilation rate [µmol CO2 m-2 s-1] 
                              gross_assimilation_C13    , &  !< canopy-integrated gross assimilation rate [µmol 13CO2 m-2 s-1]
                              gross_assimilation_C14    , &  !< canopy-integrated gross assimilation rate [µmol 14CO2 m-2 s-1]
                              net_assimilation          , &  !< canopy-integrated net assimilation rate [µmol CO2 m-2 s-1] 
                              net_assimilation_boc      , &  !< net assimilation rate at the bottom of the canopy [µmol CO2 m-2 s-1] 
                              maint_respiration_leaf         !< canopy-integrated foliar respiration rate [µmol CO2 m-2 s-1]

    TYPE(t_jsb_var_real3d)    :: &
                              gross_assimilation_cl     , &  !< gross assimilation rate per canopy layer [µmol CO2 m-2 s-1] 
                              net_assimilation_cl       , &  !< net assimilation rate per canopy layer [µmol CO2 m-2 s-1]
                              maint_respiration_leaf_cl      !< foliar respiration rate per canopy layer [µmol CO2 m-2 s-1]

    ! other 
    TYPE(t_jsb_var_real2d)    :: &
                              beta_air                  , &  !< scaling factor to account for air humidity effects on conductance [unitless]
                              beta_soa                  , &  !< scaling factor to account for spring in coniferous forests [unitless]
                              beta_soil_ps              , &  !< scaling factor to account for soil moisture constraints on photosynthesis [unitless]
                              beta_soil_gs              , &  !< scaling factor to account for soil moisture constraints on conductance [unitless]
                              canopy_cond               , &  !< canopy-integrated conductance for water [m s-1]
                              co2_conc_leaf             , &  !< canopy-integrated internal leaf CO2 concentration [ppm]
                              t_jmax_opt                     !< optimal temperature for electron transport of photosynthesis [°C]

    TYPE(t_jsb_var_real3d)    :: &
                              canopy_cond_cl            , &  !< conductance for water per canopy layer [m s-1]
                              co2_conc_leaf_cl               !< internal leaf CO2 concentration per canopy layer [ppm]

    ! diagnostics
    TYPE(t_jsb_var_real3d)    :: &
                              jmax_cl                   , &  !< maximum rate of electron transport [µmol CO2 m-2 s-1]
                              vcmax_cl                       !< maximum rate of carboxylation [µmol CO2 m-2 s-1]

    ! averaging
    TYPE(t_jsb_var_real2d)    :: &
                              soa_tsoa_mavg                  !< state of acclimation, used in calculation of beta_soa time-averaged

    TYPE(t_jsb_var_real2d)    :: &
                              beta_air_daytime          , &  !< scaling factor to account for air humidity effects on conductance average of the previous day [unitless]
                              beta_air_daytime_dacc     , &  !< scaling factor to account for air humidity effects on conductance average of the previous day [unitless]
                              beta_air_tfrac_mavg       , &  !< scaling factor to account for air humidity effects on conductance time-averaged [unitless]
                              beta_air_tcnl_mavg        , &  !< scaling factor to account for air humidity effects on conductance time-averaged [unitless]
                              beta_soa_tphen_mavg            !< scaling factor to account for spring in coniferous forests time-averaged [unitless]

    TYPE(t_jsb_var_real2d)    :: &
                              beta_soil_ps_daytime      , &  !< scaling factor to account for soil moisture constraints on photosynthesis average of the previous day [unitless]
                              beta_soil_ps_daytime_dacc , &  !< scaling factor to account for soil moisture constraints on photosynthesis average of the previous day [unitless]
                              beta_soil_ps_tfrac_mavg   , &  !< scaling factor to account for soil moisture constraints on photosynthesis time-averaged [unitless]
                              beta_soil_ps_tcnl_mavg    , &  !< scaling factor to account for soil moisture constraints on conductance time-averaged [unitless]
                              beta_soil_gs_daytime      , &  !< scaling factor to account for soil moisture constraints on conductance average of the previous day [unitless]
                              beta_soil_gs_daytime_dacc , &  !< scaling factor to account for soil moisture constraints on conductance average of the previous day [unitless]
                              beta_soil_gs_tphen_mavg   , &  !< scaling factor to account for soil moisture constraints on conductance time-averaged [unitless]
                              beta_soil_gs_tfrac_mavg   , &  !< scaling factor to account for soil moisture constraints on conductance time-averaged [unitless]
                              beta_soil_gs_tcnl_mavg         !< scaling factor to account for soil moisture constraints on conductance time-averaged [unitless]


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
    USE mo_slm_varlist,       ONLY: jsb_add_var   ! replace by mo_jsb_varlist later    *quincytemp
    USE mo_jsb_io,            ONLY: grib_bits, t_cf, t_grib1, t_grib2   ! , TSTEP_CONSTANT, tables
    USE mo_jsb_grid_class,    ONLY: t_jsb_grid, t_jsb_vgrid
    USE mo_jsb_grid,          ONLY: Get_grid, Get_vgrid
    USE mo_jsb_model_class,   ONLY: t_jsb_model

    IMPLICIT NONE
    ! ---------------------------
    ! 0.1.InOut
    CLASS(t_assimi_memory), INTENT(inout), TARGET  :: mem          !> assimi memory
    CHARACTER(len=*),     INTENT(in)            :: prefix          !> process name
    CHARACTER(len=*),     INTENT(in)            :: suffix          !> tile name
    INTEGER,              INTENT(in)            :: lct_ids !(:)    !> in jsbach4: vector of land cover type IDs
    INTEGER,              INTENT(in)            :: model_id        !> model ID model\%id
    ! ---------------------------
    ! 0.2 Local
    TYPE(t_jsb_model), POINTER :: model
    TYPE(t_jsb_grid),  POINTER :: hgrid                        ! Horizontal grid
    TYPE(t_jsb_vgrid), POINTER :: vgrid_surface                ! Vertical grid
    TYPE(t_jsb_vgrid), POINTER :: vgrid_canopy                 ! Vertical grid
    INTEGER                    :: table
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':Init_assimi_memory'
    ! ---------------------------
    ! 0.5 Get Memory
    model         => Get_model(model_id)
    table         = 1     ! in jsbach4 this line is: table = tables(1)     *quincytemp
    hgrid         => Get_grid(model%grid_id)
    vgrid_surface => Get_vgrid('surface')
    vgrid_canopy  => Get_vgrid('canopy_layer')



    !> 1.0  init memory
    !! 
    IF (.NOT. ASSOCIATED(mem%gross_assimilation%ptr)) THEN

      !> add variables to assimi memory assimi\%VAR
      !! 
      !! jsb_add_var -> jsb_add_var_list_element_r1d    *slmtemp 1D -> 2D \n
      !! and \n
      !! jsb_add_var -> jsb_add_var_list_element_r1d    *slmtemp 1D -> 3D
      !!  
      CALL jsb_add_var(mem, 'gross_assimilation', mem%gross_assimilation, &
            hgrid, vgrid_surface, &
            t_cf('gross_assimilation', 'µmol CO2 m-2 s-1', 'canopy-integrated gross assimilation rate'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'gross_assimilation_C13', mem%gross_assimilation_C13, &
            hgrid, vgrid_surface, &
            t_cf('gross_assimilation_C13', 'µmol 13CO2 m-2 s-1', 'canopy-integrated gross assimilation rate'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'gross_assimilation_C14', mem%gross_assimilation_C14, &
            hgrid, vgrid_surface, &
            t_cf('gross_assimilation_C14', 'µmol 14CO2 m-2 s-1', 'canopy-integrated gross assimilation rate'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'net_assimilation', mem%net_assimilation, &
            hgrid, vgrid_surface, &
            t_cf('net_assimilation', 'µmol CO2 m-2 s-1', 'canopy-integrated net assimilation rate'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'net_assimilation_boc', mem%net_assimilation_boc, &
            hgrid, vgrid_surface, &
            t_cf('net_assimilation_boc', 'µmol CO2 m-2 s-1', 'net assimilation rate at the bottom of the canopy'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'maint_respiration_leaf', mem%maint_respiration_leaf, &
            hgrid, vgrid_surface, &
            t_cf('maint_respiration_leaf', 'µmol CO2 m-2 s-1', 'canopy-integrated foliar respiration rate'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue

      CALL jsb_add_var(mem, 'gross_assimilation_cl', mem%gross_assimilation_cl, &
            hgrid, vgrid_canopy, &
            t_cf('gross_assimilation_cl', 'µmol CO2 m-2 s-1', 'gross assimilation rate per canopy layer'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'net_assimilation_cl', mem%net_assimilation_cl, &
            hgrid, vgrid_canopy, &
            t_cf('net_assimilation_cl', 'µmol CO2 m-2 s-1', 'net assimilation rate per canopy layer'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'maint_respiration_leaf_cl', mem%maint_respiration_leaf_cl, &
            hgrid, vgrid_canopy, &
            t_cf('maint_respiration_leaf_cl', 'µmol CO2 m-2 s-1', 'foliar respiration rate per canopy layer'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue

      CALL jsb_add_var(mem, 'beta_air', mem%beta_air, &
            hgrid, vgrid_surface, &
            t_cf('beta_air', 'unitless', 'scaling factor to account for air humidity effects on conductance'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
 
      CALL jsb_add_var(mem, 'beta_soa', mem%beta_soa, &
            hgrid, vgrid_surface, &
            t_cf('beta_soa', 'unitless', 'scaling factor to account for spring in coniferous forests'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue

      CALL jsb_add_var(mem, 'beta_soil_ps', mem%beta_soil_ps, &
            hgrid, vgrid_surface, &
            t_cf('beta_soil_ps', 'unitless', 'scaling factor to account for soil moisture constraints on photosynthesis'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'beta_soil_gs', mem%beta_soil_gs, &
            hgrid, vgrid_surface, &
            t_cf('beta_soil_gs', 'unitless', 'scaling factor to account for soil moisture constraints on conductance'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue

      CALL jsb_add_var(mem, 'canopy_cond', mem%canopy_cond, &
            hgrid, vgrid_surface, &
            t_cf('canopy_cond', 'm s-1', 'canopy-integrated conductance for water'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'co2_conc_leaf', mem%co2_conc_leaf, &
            hgrid, vgrid_surface, &
            t_cf('co2_conc_leaf', 'ppm', 'canopy-integrated internal leaf CO2 concentration'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue

      CALL jsb_add_var(mem, 't_jmax_opt', mem%t_jmax_opt, &
            hgrid, vgrid_surface, &
            t_cf('t_jmax_opt', 'deg C', 'optimal temperature for electron transport of photosynthesis'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
 
      CALL jsb_add_var(mem, 'canopy_cond_cl', mem%canopy_cond_cl, &
            hgrid, vgrid_canopy, &
            t_cf('canopy_cond_cl', 'm s-1', 'conductance for water per canopy layer'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
 
      CALL jsb_add_var(mem, 'co2_conc_leaf_cl', mem%co2_conc_leaf_cl, &
            hgrid, vgrid_canopy, &
            t_cf('co2_conc_leaf_cl', 'ppm', 'internal leaf CO2 concentration per canopy layer'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
 
      CALL jsb_add_var(mem, 'jmax_cl', mem%jmax_cl, &
            hgrid, vgrid_canopy, &
            t_cf('jmax_cl', 'µmol CO2 m-2 s-1', 'maximum rate of electron transport'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
 
      CALL jsb_add_var(mem, 'vcmax_cl', mem%vcmax_cl, &
            hgrid, vgrid_canopy, &
            t_cf('vcmax_cl', 'µmol CO2 m-2 s-1', 'maximum rate of carboxylation'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue

      CALL jsb_add_var(mem, 'soa_tsoa_mavg', mem%soa_tsoa_mavg, &
            hgrid, vgrid_surface, &
            t_cf('soa_tsoa_mavg', '', 'state of acclimation, used in calculation of beta_soa time-averaged'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue

      CALL jsb_add_var(mem, 'beta_air_daytime', mem%beta_air_daytime, &
            hgrid, vgrid_surface, &
            t_cf('beta_air_daytime', 'unitless', 'scaling factor to account for air humidity effects on conductance average of the previous day'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'beta_air_daytime_dacc', mem%beta_air_daytime_dacc, &
            hgrid, vgrid_surface, &
            t_cf('beta_air_daytime_dacc', 'unitless', 'scaling factor to account for air humidity effects on conductance average of the previous day'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'beta_air_tfrac_mavg', mem%beta_air_tfrac_mavg, &
            hgrid, vgrid_surface, &
            t_cf('beta_air_tfrac_mavg', 'unitless', 'scaling factor to account for air humidity effects on conductance time-averaged'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'beta_air_tcnl_mavg', mem%beta_air_tcnl_mavg, &
            hgrid, vgrid_surface, &
            t_cf('beta_air_tcnl_mavg', 'unitless', 'scaling factor to account for air humidity effects on conductance time-averaged'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue

      CALL jsb_add_var(mem, 'beta_soa_tphen_mavg', mem%beta_soa_tphen_mavg, &
            hgrid, vgrid_surface, &
            t_cf('beta_soa_tphen_mavg', 'unitless', 'scaling factor to account for spring in coniferous forests time-averaged'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue

      CALL jsb_add_var(mem, 'beta_soil_ps_daytime', mem%beta_soil_ps_daytime, &
            hgrid, vgrid_surface, &
            t_cf('beta_soil_ps_daytime', 'unitless', 'scaling factor to account for soil moisture constraints on photosynthesis average of the previous day'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'beta_soil_ps_daytime_dacc', mem%beta_soil_ps_daytime_dacc, &
            hgrid, vgrid_surface, &
            t_cf('beta_soil_ps_daytime_dacc', 'unitless', 'scaling factor to account for soil moisture constraints on photosynthesis average of the previous day'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'beta_soil_ps_tfrac_mavg', mem%beta_soil_ps_tfrac_mavg, &
            hgrid, vgrid_surface, &
            t_cf('beta_soil_ps_tfrac_mavg', 'unitless', 'scaling factor to account for soil moisture constraints on photosynthesis time-averaged'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'beta_soil_ps_tcnl_mavg', mem%beta_soil_ps_tcnl_mavg, &
            hgrid, vgrid_surface, &
            t_cf('beta_soil_ps_tcnl_mavg', 'unitless', 'scaling factor to account for soil moisture constraints on conductance time-averaged'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'beta_soil_gs_daytime', mem%beta_soil_gs_daytime, &
            hgrid, vgrid_surface, &
            t_cf('beta_soil_gs_daytime', 'unitless', 'scaling factor to account for soil moisture constraints on conductance average of the previous day'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'beta_soil_gs_daytime_dacc', mem%beta_soil_gs_daytime_dacc, &
            hgrid, vgrid_surface, &
            t_cf('beta_soil_gs_daytime_dacc', 'unitless', 'scaling factor to account for soil moisture constraints on conductance average of the previous day'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'beta_soil_gs_tphen_mavg', mem%beta_soil_gs_tphen_mavg, &
            hgrid, vgrid_surface, &
            t_cf('beta_soil_gs_tphen_mavg', 'unitless', 'scaling factor to account for soil moisture constraints on conductance time-averaged'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'beta_soil_gs_tfrac_mavg', mem%beta_soil_gs_tfrac_mavg, &
            hgrid, vgrid_surface, &
            t_cf('beta_soil_gs_tfrac_mavg', 'unitless', 'scaling factor to account for soil moisture constraints on conductance time-averaged'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'beta_soil_gs_tcnl_mavg', mem%beta_soil_gs_tcnl_mavg, &
            hgrid, vgrid_surface, &
            t_cf('beta_soil_gs_tcnl_mavg', 'unitless', 'scaling factor to account for soil moisture constraints on conductance time-averaged'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue

    END IF

  END SUBROUTINE Init_assimi_memory

#endif
END MODULE mo_assimi_memory_class
