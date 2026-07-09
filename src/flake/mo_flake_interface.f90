!>
!! @brief Contains the interfaces to the processes related to flake systems
!!
!! @author
!!  Zoé Rehder, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Zoé Rehder (2022-08-17)
!!
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/license .
!!

MODULE mo_flake_interface
#ifndef __NO_JSBACH__

  ! -------------------------------------------------------------------------------------------------------
  ! Used variables of module
  
  ! Use of basic structures
  USE mo_jsb_control,     ONLY: debug_on, jsbach_runs_standalone
  USE mo_kind,            ONLY: wp
  USE mo_exception,       ONLY: message

  USE mo_jsb_model_class,    ONLY: t_jsb_model
  USE mo_jsb_class,          ONLY: Get_model  
  USE mo_jsb_grid_class,     ONLY: t_jsb_grid
  USE mo_jsb_tile_class,     ONLY: t_jsb_tile_abstract, t_jsb_aggregator
  USE mo_jsb_config_class,   ONLY: t_jsb_config, t_jsb_config_p
  USE mo_jsb_process_class,  ONLY: t_jsb_process
  USE mo_jsb_task_class,     ONLY: t_jsb_process_task, t_jsb_task_options

  ! Use of processes in this module
  dsl4jsb_Use_processes FLAKE_, A2L_, SEB_, HYDRO_, RAD_

  ! Use of process configurations
  dsl4jsb_Use_config(FLAKE_)

  ! Use of process memories
  dsl4jsb_Use_memory(FLAKE_)
  dsl4jsb_Use_memory(A2L_)
  dsl4jsb_Use_memory(SEB_)
  dsl4jsb_Use_memory(HYDRO_)
  dsl4jsb_Use_memory(RAD_)

  ! -------------------------------------------------------------------------------------------------------
  ! Module variables 
  
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: Register_flake_tasks !, t_flake_process

  CHARACTER(len=*), PARAMETER :: modname = 'mo_flake_interface'
  
  !> Type definition for flake process
!   TYPE, EXTENDS(t_jsb_process) :: t_flake_process
!   CONTAINS
!     PROCEDURE :: Register_tasks => Register_flake_tasks        !< Registers tasks for flake process
!   END TYPE t_flake_process
  
  !> Type definition for flake task
  TYPE, EXTENDS(t_jsb_process_task) :: tsk_flake_wbphysics
  CONTAINS
    PROCEDURE, NOPASS :: Integrate => update_flake_wbphysics    !< Advances task computation for one timestep
    PROCEDURE, NOPASS :: Aggregate => aggregate_flake_wbphysics !< Aggregates computed task variables
  END TYPE tsk_flake_wbphysics
  
  !> Constructor interface for flake_wbphysics
  INTERFACE tsk_flake_wbphysics
    PROCEDURE Create_task_flake_wbphysics                        !< Constructor function for task
  END INTERFACE tsk_flake_wbphysics

  !X Repeat here for each additional task of the process the two blocks above ...
  
CONTAINS

  ! ================================================================================================================================
  !! Constructors for tasks
  
  ! -------------------------------------------------------------------------------------------------------
  !> Constructor for flake task
  !!
  !! @param[in]     model_id     Model id
  !! @return        return_ptr   Instance of process task "flake"
  !!
  FUNCTION Create_task_flake_wbphysics(model_id) RESULT(return_ptr)

    INTEGER,                   INTENT(in) :: model_id
    CLASS(t_jsb_process_task), POINTER    :: return_ptr

    ALLOCATE(tsk_flake_wbphysics::return_ptr)
    CALL return_ptr%Construct(name='flake_wbphysics', process_id=FLAKE_, owner_model_id=model_id)

  END FUNCTION Create_task_flake_wbphysics
  

  !X Repeat here for each additional task of the process the function above ...
    
  ! -------------------------------------------------------------------------------------------------------
  !> Register tasks for flake process
  !!
  !! @param[in,out] process      Instance of flake process class  
  !! @param[in]     model_id  Model id  
  !!
  SUBROUTINE Register_flake_tasks(process, model_id)

    CLASS(t_jsb_process),    INTENT(inout) :: process
    INTEGER,                 INTENT(in)    :: model_id

    CALL process%Register_task(tsk_flake_wbphysics(model_id))
    !X Repeat the line above here for each additional task of the process ...
    
  END SUBROUTINE Register_flake_tasks

  ! ================================================================================================================================
  !>
  !! @brief Implementation of "update" for task "flake_wbphysics"
  !! Task "flake_wbphysics" <EXPLANATIION_OF_PROCESS>.
  !!
  !! @param[in,out] tile    Tile for which routine is executed.
  !! @param[in]     options Additional run-time parameters.
  !!
  SUBROUTINE update_flake_wbphysics(tile, options)
    USE mo_util,                ONLY: real2string

    USE mo_jsb_grid,      ONLY: Get_grid  
    USE mo_phy_schemes,   ONLY: lwnet_from_lwdown, qsat_water, qsat_ice, surface_dry_static_energy


    USE mo_flake_process, ONLY: calc_surface_temperature     , &
                               & calc_radflux                 , &
                               & calc_ice_water_heatfluxes    , &
                               & calc_lake_ice                , &
                               & calc_lake_temp_stratification

    USE mo_flake_SfcFlx_process, ONLY: calc_momsenlat

    USE mo_flake_constants , ONLY : & 
      rot_earth             , & ! Rotation rate of the earth [rad/s]
      tpl_T_f                   ! Fresh water freezing point [K]  
    
    USE mo_flake_parameters, ONLY : &
      h_Ice_min             , & ! Minimum ice thickness
      h_Snow_min                ! Minimum snow thickness

    USE mo_jsb_physical_constants, ONLY: &
      & rhoh2o, & !> [kg/m3]  density of liquid water !! = tpl_rho_w_r
      & cpd, cvd

    ! Arguments
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    ! Local variables
    INTEGER                     :: iblk, ics, ice, nc
    REAL(wp), POINTER           :: lat(:)
    REAL(wp)                    :: dtime
    CHARACTER(len=*), PARAMETER :: routine = modname//':update_flake_wbphysics'
    REAL(wp), PARAMETER         :: h_atmvars = 2._wp !< height of atmospheric forcing variables (wind_air, t_air, q_air)
    TYPE(t_jsb_model), POINTER  :: model
    TYPE(t_jsb_grid),  POINTER  :: grid


    ! Declare pointers to process configs and memory
    dsl4jsb_Def_config(FLAKE_)
    dsl4jsb_Def_memory(A2L_)
    dsl4jsb_Def_memory(SEB_)
    dsl4jsb_Def_memory(FLAKE_) 
    dsl4jsb_Def_memory(HYDRO_) 
    dsl4jsb_Def_memory(RAD_) 


    ! Declare pointers to variables in memory
    dsl4jsb_Real2D_onChunk :: &
      & q_bot,                & !< Heat flux through the water-bottom sediment interface of lakes   [W m-2]
    ! atmospheric forcing
      & wind_air,             & !< 2m air wind speed                                                [m s-1]
      & t_air,                & !< 2m air temperature                                               [K]
      & q_air,                & !< 2m air specific humidity                                         [g g-1]
      & press_srf,            & !< surface air pressure                                             [Pa]
      & lw_srf_down,          & !< longwave surface downward flux                                   [W m-2]
      & swvis_srf_down,       & !< shortwave surface downward flux in the visible band              [W m-2]
      & swnir_srf_down,       & !< shortwave surface downward flux in the near infrared band        [W m-2]
    ! heat fluxes and atm exchange
      & albedo_lwtr,          & !< albedo of lake water
      & albedo_lice,          & !< albedo of ice on lake water
      & latent_hflx_wtr,      & !< latent heat flux over water                                      [W m-2]
      & sensible_hflx_wtr,    & !< sensible heat flux over water                                    [W m-2]
      & latent_hflx,          & !< latent heat flux                                                 [W m-2]
      & sensible_hflx,        & !< sensible heat flux                                    [W m-2]
      & qsat_lwtr,            & !< saturation specific humidity at surface over lake water          [kg kg-1]
      & s_lwtr,               & !< surface dry static energy over lake water                        [m2 s-2]
      & latent_hflx_ice,      & !< latent heat flux over ice                                        [W m-2]
      & sensible_hflx_ice,    & !< sensible heat flux over ice                                      [W m-2]
      & qsat_lice,            & !< saturation specific humidity at surface over lake ice            [kg kg-1]
      & s_lice,               & !< surface dry static energy over lake ice                          [m2 s-2]
      & evapo_wtr,            & !< evaporation over lake water
      & evapo_ice,            & !< evaporation over lake water
    ! topography and location
    !   & depth_w,              &
      & fetch,                &
      & par_Coriolis,         &
    !   & pondtype,             &
    ! FLake inout
    !   & t_snow,               &
      & t_lice,               &
      & t_mnw,                &
      & t_wml,                &
      & t_bot,                &
      & c_t,                  &
      & h_ice,                &
      & h_snow,               &
      & h_ml,                 &
      & T_lwtr,                &
      & z0t_sf,                &
      & z0u_sf,                &          
      & w_star

    REAL(wp) ::                   &
        ! & albedo_water,           & 
        ! & albedo_lice(options%nc), & 
        & albedo_snow(options%nc),& 
    !   & w_star(options%nc),     &
        & u_star(options%nc),     &
        & u_star_sf(options%nc),  &
        & I_atm(options%nc),      &
        & I_bot(options%nc),      &
        & I_h(options%nc),        &
        & I_ice(options%nc),      &
        & I_intm_0_h(options%nc), &
        & I_intm_h_D(options%nc), &
        & I_snow(options%nc),     & ! might not be needed
        & I_w(options%nc),        & 
        & Q_ice(options%nc),      & 
        & Q_momentum(options%nc), &
        & Q_sensible(options%nc), &
        & Q_latent(options%nc),   &
        & Q_watvap(options%nc),   &
        & Q_w(options%nc),        &
        & Q_star(options%nc),     &
        & Phi_T_pr0(options%nc)
    
    REAL(wp) :: cpd_or_cvd
    LOGICAL :: l_ice_create(options%nc), &
      & use_tmx, &
      & jsb_standalone

    ! Get local variables from options argument
    iblk    = options%iblk
    ics     = options%ics
    ice     = options%ice
    nc      = options%nc
    dtime   = options%dtime

    ! If process is not active on this tile, do nothing
    IF (.NOT. tile%Is_process_active(FLAKE_)) RETURN

    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    model => Get_model(tile%owner_model_id)
    grid => Get_grid(model%grid_id)
    lat  => grid%lat(ics:ice, iblk) 

    use_tmx = model%config%use_tmx
    IF (use_tmx) THEN
      cpd_or_cvd = cvd
    ELSE
      cpd_or_cvd = cpd
    END IF

    jsb_standalone = jsbach_runs_standalone()


    ! Get process memories
    dsl4jsb_Get_memory(FLAKE_)
    dsl4jsb_Get_memory(A2L_)
    dsl4jsb_Get_memory(SEB_)
    dsl4jsb_Get_memory(HYDRO_)
    dsl4jsb_Get_memory(RAD_)

    ! Get process variables
    dsl4jsb_Get_var2D_onChunk(A2L_, wind_air)           ! in
    dsl4jsb_Get_var2D_onChunk(A2L_, t_air)              ! in
    dsl4jsb_Get_var2D_onChunk(A2L_, q_air)              ! in
    dsl4jsb_Get_var2D_onChunk(A2L_, press_srf)          ! in
    dsl4jsb_Get_var2D_onChunk(A2L_, lw_srf_down)        ! in
    dsl4jsb_Get_var2D_onChunk(A2L_, swvis_srf_down)     ! in
    dsl4jsb_Get_var2D_onChunk(A2L_, swnir_srf_down)     ! in

    dsl4jsb_Get_var2D_onChunk(SEB_, sensible_hflx_wtr)  ! in
    dsl4jsb_Get_var2D_onChunk(SEB_, sensible_hflx_ice)  ! in
    dsl4jsb_Get_var2D_onChunk(SEB_, sensible_hflx)      ! in
    dsl4jsb_Get_var2D_onChunk(SEB_, latent_hflx_wtr)    ! in
    dsl4jsb_Get_var2D_onChunk(SEB_, latent_hflx_ice)    ! in
    dsl4jsb_Get_var2D_onChunk(SEB_, latent_hflx)        ! in
    dsl4jsb_Get_var2D_onChunk(HYDRO_, evapo_wtr)        ! in
    dsl4jsb_Get_var2D_onChunk(HYDRO_, evapo_ice)        ! in

    dsl4jsb_Get_var2D_onChunk(RAD_, albedo_lwtr)        ! in
    dsl4jsb_Get_var2D_onChunk(RAD_, albedo_lice)        ! in

    dsl4jsb_Get_var2D_onChunk(SEB_, qsat_lwtr)          ! out 
    dsl4jsb_Get_var2D_onChunk(SEB_, qsat_lice)          ! out 
    dsl4jsb_Get_var2D_onChunk(SEB_, s_lwtr)             ! out 
    dsl4jsb_Get_var2D_onChunk(SEB_, s_lice)             ! out 
    dsl4jsb_Get_var2D_onChunk(SEB_, T_lwtr)             ! inout

    dsl4jsb_Get_var2D_onChunk(FLAKE_, z0u_sf)           ! out
    dsl4jsb_Get_var2D_onChunk(FLAKE_, z0t_sf)           ! out

    ! dsl4jsb_Get_var2D_onChunk(FLAKE_, depth_w)         ! in
    dsl4jsb_Get_var2D_onChunk(FLAKE_, Q_bot)           ! in
    dsl4jsb_Get_var2D_onChunk(FLAKE_, fetch)           ! in
    ! dsl4jsb_Get_var2D_onChunk(FLAKE_, pondtype)        ! in

    ! dsl4jsb_Get_var2D_onChunk(FLAKE_, t_snow)        ! inout
    dsl4jsb_Get_var2D_onChunk(FLAKE_, t_lice)          ! inout
    dsl4jsb_Get_var2D_onChunk(FLAKE_, t_mnw)           ! inout
    dsl4jsb_Get_var2D_onChunk(FLAKE_, t_wml)           ! inout
    dsl4jsb_Get_var2D_onChunk(FLAKE_, t_bot)           ! inout
    dsl4jsb_Get_var2D_onChunk(FLAKE_, c_t)             ! inout
    dsl4jsb_Get_var2D_onChunk(FLAKE_, h_ice)           ! inout
    dsl4jsb_Get_var2D_onChunk(FLAKE_, h_snow)          ! inout
    dsl4jsb_Get_var2D_onChunk(FLAKE_, h_ml)            ! inout
    
    dsl4jsb_Get_var2D_onChunk(FLAKE_, w_star)          ! out


    ! ---------------------------
    ! Go    

    !  Compute coriolis parameter
    par_Coriolis => grid%coriolis(ics:ice, iblk)

    !  Set albedos of the lake water, lake ice and snow
    !  They are now computed in the RAD process

    ! Snow is not considered
    albedo_snow  = albedo_lice
    

    I_atm = swvis_srf_down + swnir_srf_down

    !  Compute solar radiation fluxes (positive downward)    (needed?)
    ! I_atm_flk = I_atm_in
    ! Q_bot_flk = Q_bot_in

    call calc_radflux ( &
                        !& depth_w,        & ! in
                        & albedo_lwtr,    & ! in
                        & albedo_lice,    & ! in
                        & albedo_snow,    & ! in
                        & h_ice,          & ! in 
                        & h_ML,           & ! in 
                        & h_snow,         & ! in
                        & I_atm,          & ! in
                        & I_bot,          & ! out
                        & I_h,            & ! out
                        & I_ice,          & ! out
                        & I_intm_0_h,     & ! out
                        & I_intm_h_D,     & ! out
                        & I_snow,         & ! out (but maybe not needed) 
                        & I_w             & ! out
                        )

    !  Compute long-wave radiation fluxes (positive downward)
    Q_w = lwnet_from_lwdown(lw_srf_down, T_lwtr)     ! Radiation of the atmosphere 

    !  Compute the surface friction velocity and fluxes of sensible and latent heat 
    ! CALL calc_momsenlat ( h_atmvars, h_atmvars, fetch, wind_air, t_air, q_air, T_lwtr, press_srf, h_ice,  &
    !                         Q_momentum, Q_sensible, Q_latent, Q_watvap, z0u_sf, z0t_sf, u_star ) 
    Q_momentum = 0._wp
    u_star_sf = SQRT(-Q_momentum/rhoh2o) !zr todo: Check differences between u_star and u_star_sf

    ! WHERE(h_ice.GE.h_Ice_min)  
    !   latent_hflx_wtr   = 0._wp
    !   sensible_hflx_wtr = 0._wp
    !   evapo_wtr         = 0._wp
    !   latent_hflx_ice   = Q_latent
    !   sensible_hflx_ice = Q_sensible
    !   evapo_ice        = Q_watvap
    ! ELSEWHERE
    !   latent_hflx_wtr   = Q_latent
    !   sensible_hflx_wtr = Q_sensible
    !   evapo_wtr         = Q_watvap
    !   latent_hflx_ice   = 0._wp
    !   sensible_hflx_ice = 0._wp
    !   evapo_ice        = 0._wp
    ! END WHERE

    !  Compute heat fluxes Q_snow, Q_ice, Q_w

    ! Q_w = Q_w - Q_sensible - Q_latent  ! Add sensible and latent heat fluxes (notice the signs)
    Q_w = Q_w + sensible_hflx + latent_hflx

    WHERE (h_ice >= h_Ice_min .AND. h_snow >= h_Snow_min)             ! Ice and snow exist
        !Q_snow = Q_w
        Q_ice  = 0._wp
        Q_w    = 0._wp
    ELSEWHERE (h_ice >= h_Ice_min .AND. h_snow < h_Snow_min)          ! Ice exists but no snow cover
        !Q_snow = 0._wp
        Q_ice  = Q_w
        Q_w    = 0._wp
    ELSEWHERE                                           ! No ice-snow cover
        !Q_snow = 0._wp
        Q_ice  = 0._wp
    END WHERE
    !Q_snow = 0._wp ! Snow is not considered for now

    !  Advance FLake variables 
    !  From here on the original file flake_driver is included but split into the routines

    call calc_ice_water_heatfluxes( h_ice,      & ! in
                                  & h_ML,       & ! in
                                  & T_bot,      & ! in
                                  & T_wML,      & ! in
                                  & C_T,        & ! in  ! & depth_w,    & ! in
                                  & I_w,        & ! in
                                  & I_h,        & ! in
                                  & I_intm_0_h, & ! in
                                  & I_bot,      & ! inout
                                  & Q_w,        & ! out
                                  & Q_star,     & ! out
                                  & Phi_T_pr0   & ! out
                                  )



    call calc_lake_ice( I_ice,        & ! in
                      & I_w,          & ! in 
                      & T_wML,        & ! in 
                      & Q_ice,        & ! in     
                      & dtime,        & ! in
                      & Q_w,          & ! inout
                      & h_ice,        & ! inout
                      & T_lice,       & ! inout
                      & l_ice_create  & ! out
                      )

    call calc_lake_temp_stratification( dtime,                & ! in & depth_w,  & extincoef_optic_water,   & ! in
                                      & l_ice_create,         & ! in
                                      & h_ice,                & ! in
                                      & I_bot,                & ! in
                                      & I_intm_0_h,           & ! in
                                      & I_intm_h_d,           & ! in
                                      & I_w,                  & ! in
                                      & I_h,                  & ! in
                                      & par_Coriolis,         & ! in
                                      & Q_bot,                & ! in
                                      & Q_star,               & ! in
                                      & u_star_sf,            & ! in 
                                      & C_T,                  & ! inout
                                      & Q_w,                  & ! inout
                                      & T_bot,                & ! inout
                                      & T_mnw,                & ! inout
                                      & T_wML,                & ! inout
                                      & h_ML,                 & ! inout
                                      & w_star                & ! out
                                      )    

    call calc_surface_temperature(h_snow, & ! in
                                & T_lice, & ! in /instead of t_snow
                                & h_ice,  & ! in
                                & T_lice, & ! in
                                & T_wML,  & ! in
                                & T_lwtr  ) ! out    

    CALL message(TRIM(routine), &
      & 'Lake temperatures, min: '&
      & //real2string(MINVAL(T_wML)) //'; max: '//real2string(MAXVAL(T_wML))  )    
    qsat_lwtr = qsat_water(T_wML, press_srf)
    s_lwtr(:)    = surface_dry_static_energy(T_wML(:), qsat_lwtr(:), cpd_or_cvd, jsb_standalone)
    WHERE(h_ice >= h_Ice_min)             ! Ice exists
      qsat_lice(:) = qsat_ice(t_lice(:), press_srf(:))
      s_lice(:)    = surface_dry_static_energy(t_lice(:), qsat_lice(:), cpd_or_cvd, jsb_standalone) 
    ELSEWHERE
      qsat_lice(:) = 0._wp
      s_lice(:)    = 0._wp
    END WHERE

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')
    
  END SUBROUTINE update_flake_wbphysics
  
  ! -------------------------------------------------------------------------------------------------------
  !>
  !! Implementation of "aggregate" for task "flake"
  !! 
  !! @param[in,out] tile    Tile for which routine is executed.  
  !! @param[in]     options Additional run-time parameters.  
  !!
  SUBROUTINE aggregate_flake_wbphysics(tile, options)

    ! Arguments
    TYPE(t_jsb_model), POINTER :: model
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    dsl4jsb_Def_memory(FLAKE_)

    CLASS(t_jsb_aggregator), POINTER :: weighted_by_fract

    ! Local variables
    INTEGER  :: iblk, ics, ice, nc

    CHARACTER(len=*), PARAMETER :: routine = modname//':aggregate_flake_wbphysics'
    
    ! Get local variables from options argument
    iblk = options%iblk
    ics  = options%ics
    ice  = options%ice
    nc   = options%nc
    
    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    dsl4jsb_Get_memory(FLAKE_)

    model             => Get_model(tile%owner_model_id)    
    weighted_by_fract => tile%Get_aggregator("weighted_by_fract")

    CALL weighted_by_fract%BeginAggregate()

    dsl4jsb_Aggregate_onChunk(FLAKE_, h_ice,        weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(FLAKE_, h_ml,         weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(FLAKE_, h_snow,       weighted_by_fract)
    ! dsl4jsb_Aggregate_onChunk(FLAKE_, T_sfc,        weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(FLAKE_, t_lice,       weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(FLAKE_, t_mnw,        weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(FLAKE_, t_wml,        weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(FLAKE_, t_bot,        weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(FLAKE_, c_t,          weighted_by_fract) 
    dsl4jsb_Aggregate_onChunk(FLAKE_, w_star,       weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(FLAKE_, z0t_sf,       weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(FLAKE_, z0u_sf,       weighted_by_fract)

    CALL weighted_by_fract%EndAggregate()

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE aggregate_flake_wbphysics

  
  ! ================================================================================================================================
  !>
  !! @brief Another task for "flake"
  !!
  !! @param[in]     options Additional run-time parameters.  
  !! @param[in,out] tile    Tile for which routine is executed.
  !  ....
  
#endif
END MODULE mo_flake_interface
