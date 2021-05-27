!>
!! @brief Contains the interfaces to the jedi process
!!
!! @author
!!  Rainer Schneck, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Rainer Schneck (2017-04-05)
!!
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/license .
!!

MODULE mo_jedi_interface
#ifndef __NO_JSBACH__

  ! -------------------------------------------------------------------------------------------------------
  ! Used variables of module

  ! Use of basic structures
  USE mo_jsb_control,     ONLY: debug_on
  USE mo_kind,            ONLY: wp
  USE mo_exception,       ONLY: message, finish

  USE mo_jsb_model_class,    ONLY: t_jsb_model
  USE mo_jsb_class,          ONLY: Get_model
  USE mo_jsb_grid_class,     ONLY: t_jsb_grid
  USE mo_jsb_tile_class,     ONLY: t_jsb_tile_abstract, t_jsb_aggregator
  USE mo_jsb_process_class,  ONLY: t_jsb_process
  USE mo_jsb_task_class,     ONLY: t_jsb_process_task, t_jsb_task_options

  USE mo_jedi_constants,     ONLY: molarMassCO2_kg, sec_per_day

  ! Use of processes in this module (Get_jedi_memory and Get_jedi_config)
  dsl4jsb_Use_processes JEDI_, ASSIMI_, PHENO_, DISTURB_, A2L_, L2A_,HYDRO_,RAD_

  ! Use of process configurations
  dsl4jsb_Use_config(JEDI_)
  dsl4jsb_Use_config(PHENO_)

  ! Use of process memories (t_jedi_memory)
  dsl4jsb_Use_memory(A2L_)
  dsl4jsb_Use_memory(ASSIMI_)
  dsl4jsb_Use_memory(PHENO_)
  dsl4jsb_Use_memory(JEDI_)
  dsl4jsb_Use_memory(L2A_)

  ! =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-  
  ! Use the ELEMENT FUNCTION from HYDRO for JeDi : get_water_stress_factor
  ! 
  dsl4jsb_Use_config(HYDRO_)
  dsl4jsb_Use_memory(HYDRO_)
  dsl4jsb_Use_memory(RAD_)
  ! =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-  
  

  ! -------------------------------------------------------------------------------------------------------
  ! Module variables

  IMPLICIT NONE
  PRIVATE

  PUBLIC ::  Register_jedi_tasks
  PUBLIC ::  recalc_per_tile_vars
  PUBLIC ::  jedi_biomass_to_fract

  CHARACTER(len=*), PARAMETER :: modname = 'mo_jedi_interface'

  !> Type definition for C_NPP_pot_allocation task
  TYPE, EXTENDS(t_jsb_process_task) ::   tsk_tradeoff
  CONTAINS
    PROCEDURE, NOPASS :: Integrate => update_tradeoff    !< Advances task computation for one timestep
    PROCEDURE, NOPASS :: Aggregate => aggregate_tradeoff !< Aggregates computed task variables
  END TYPE   tsk_tradeoff

  !> Constructor interface for tradeoff task
  INTERFACE tsk_tradeoff
    PROCEDURE Create_tradeoff        !< Constructor function for task
  END INTERFACE tsk_tradeoff
  
! =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- 
 !PHS 24.07.2020 Compute water stress in CARBON at PFT level. 
 !> Type definition for water_stress task
  TYPE, EXTENDS(t_jsb_process_task) :: tsk_water_stress
  CONTAINS
    PROCEDURE, NOPASS :: Integrate => update_water_stress
    PROCEDURE, NOPASS :: Aggregate => aggregate_water_stress
  END TYPE tsk_water_stress

  !> Constructor interface for water_stress task
  INTERFACE tsk_water_stress
    PROCEDURE Create_task_water_stress
  END INTERFACE tsk_water_stress
! =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- 

CONTAINS

  ! ================================================================================================================================
  !! Constructors for tasks

  ! -------------------------------------------------------------------------------------------------------
  !> Constructor for carbon task
  !!
  !! @param[in]     model_id     Model id
  !! @return        return_ptr   Instance of process task "carbon"
  !!
  FUNCTION Create_tradeoff(model_id) RESULT(return_ptr)

    INTEGER,                   INTENT(in) :: model_id
    CLASS(t_jsb_process_task), POINTER    :: return_ptr

    ALLOCATE(tsk_tradeoff::return_ptr)
    CALL return_ptr%Construct(name='tradeoff', process_id=JEDI_, owner_model_id=model_id)

  END FUNCTION Create_tradeoff


! =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-  
 !PHS 24.07.2020 Compute water stress in JEDI at PFT level. 
  ! -------------------------------------------------------------------------------------------------------
  !> Constructor for water_stress task
  !!
  !! @param[in]     model_id     Model id
  !! @return        return_ptr   Instance of process task "water_stress"
  !!
  FUNCTION Create_task_water_stress(model_id) RESULT(return_ptr)

    INTEGER,                   INTENT(in) :: model_id
    CLASS(t_jsb_process_task), POINTER    :: return_ptr

    ALLOCATE(tsk_water_stress::return_ptr)
    CALL return_ptr%Construct(name='water_stress', process_id=JEDI_, owner_model_id=model_id)

  END FUNCTION Create_task_water_stress
! =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- 


  ! -------------------------------------------------------------------------------------------------------
  !> Register tasks for jedi process
  !!
  !! @param[in,out] this      Instance of jedi process class
  !! @param[in]     model_id  Model id
  !!
  SUBROUTINE Register_jedi_tasks(this, model_id)

    USE mo_jsb_model_class,   ONLY: t_jsb_model
    USE mo_jsb_class,         ONLY: Get_model

    ! in/out
    CLASS(t_jsb_process), INTENT(inout) :: this
    INTEGER,                 INTENT(in) :: model_id

    ! local
    TYPE(t_jsb_model), POINTER :: model

    ! get var / objects 
    model   => Get_model(model_id)

    CALL this%Register_task(tsk_tradeoff(model_id))
    CALL this%Register_task(tsk_water_stress(model_id))

  END SUBROUTINE Register_jedi_tasks

  ! ================================================================================================================================
  !>
  !! @brief Implementation to calculate the allocation of NPP carbon to the plant and soil carbon pools
  !!        R: Corresponds to update_cbalance_bethy in JSBACH3.
  !!
  !! @param[in,out] tile    Tile for which routine is executed.
  !! @param[in]     options Additional run-time parameters.
  !!
  SUBROUTINE update_tradeoff(tile, options)

    ! Use declarations
    USE mo_jsb_grid,               ONLY: Get_grid, Get_vgrid
    USE mo_util,                   ONLY: soil_depth_to_layers_2d
    USE mo_jedi_process,           ONLY: calc_Ctradeoff,jedi_water_stress_factor,jedi_soil_depth_to_layers_1d!,drop_expensive_tissue
    USE mo_jsb_physical_constants, ONLY: tmelt  
    USE mo_jsb_time,               ONLY: t_datetime,is_newday, is_newyear, timesteps_per_day, &
      &  is_time_experiment_start, is_time_restart, get_year_at_model_start,     &
      &  get_time_experiment_start, get_year, deallocateDatetime
    USE mo_jedi_constants,         ONLY: max_pipe_transpir
    USE mo_jsb_grid_class,         ONLY:  t_jsb_vgrid !PHS: root_depth
    USE mo_hydro_util,             ONLY: get_water_in_root_zone !PHS: root_depth
    ! Arguments
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    ! Declare local variables
    TYPE(t_jsb_model), POINTER   :: model
    TYPE(t_jsb_vgrid), POINTER   :: soil_w
    INTEGER                      :: iblk, ics, ice, nc, i
    REAL(wp)                     :: dtime

    CHARACTER(len=*), PARAMETER  :: routine = modname//':update_tradeoff'

    ! =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-  
    ! Declare variavles for JeDi carbon allocation 
    ! 
    TYPE(t_datetime), POINTER    :: inidate
    INTEGER                      :: iniyear
    LOGICAL                      :: lspinup_jedi
  
    !REAL(wp), parameter  :: rootshoot_opt     = 3.58462104_wp ! median value of root-shoot ratio 50th percentile
    REAL(wp), parameter  :: rootshoot_opt     = 8.473583644190947_wp ! median value of root-shoot ratio 75th percentile
    REAL(wp), parameter  :: c_CN_woody       = 0.00303_wp * (12.01_wp / 14.0067_wp)  ! molC/molN Ratio for Woody Tissue
    REAL(wp), parameter  :: c_sapwood        = 0.05_wp
    REAL(wp), parameter  :: c_Res_N          = 0.218_wp * ( 14.0067_wp / 12.01_wp )  !! specific respiration in molC/molN/d
    REAL(wp), PARAMETER  :: Q10              = 1.8_wp          ! Empirical parameter for temperature dependence of
                                                               !   heterotrophic respiration rate. This value is 
                                                               !   taken from JB4.

    ! =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- 
    
    ! Declare process configuration and memory Pointers
    dsl4jsb_Def_config(PHENO_)
    dsl4jsb_Def_memory(PHENO_)
    dsl4jsb_Def_memory(ASSIMI_)
    dsl4jsb_Def_memory(JEDI_)
    dsl4jsb_Def_config(JEDI_)
    ! =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-  
    ! JeDi needs soil moisture from HYDRO process
    dsl4jsb_Def_memory(HYDRO_)
    dsl4jsb_Def_config(HYDRO_)
    dsl4jsb_Def_memory(RAD_)
    ! =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-  
    
    !dsl4jsb_Def_config(A2L_)
    dsl4jsb_Def_memory(A2L_)

    ! Declare pointers to variables in memory

    ! =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-  
    ! Variavles for JeDi carbon allocation 
    !
    dsl4jsb_Real2D_onChunk :: c_storage               !     Storage carbon pool 
    dsl4jsb_Real2D_onChunk :: c_leaves                !      Leaves carbon pool             
    dsl4jsb_Real2D_onChunk :: c_fnroot                !   Fine root carbon pool
    dsl4jsb_Real2D_onChunk :: c_csroot                ! Coarse root carbon pool
    dsl4jsb_Real2D_onChunk :: c_seed                  !        Seed carbon pool
    dsl4jsb_Real2D_onChunk :: c_stem                  !        Stem carbon pool
    dsl4jsb_Real2D_onChunk :: ff_grow_W    
    dsl4jsb_Real2D_onChunk :: ff_grow_T
    dsl4jsb_Real2D_onChunk :: ff_grow_G
    dsl4jsb_Real2D_onChunk :: Flg_grow
    dsl4jsb_Real2D_onChunk :: Flg_grow_G
    dsl4jsb_Real2D_onChunk :: ff_npp
    dsl4jsb_Real2D_onChunk :: flx_germ
    dsl4jsb_Real2D_onChunk :: adj_fraction
    dsl4jsb_Real2D_onChunk :: jedi_light_constraint
    dsl4jsb_Real2D_onChunk :: Flg_die
    dsl4jsb_Real2D_onChunk :: Flg_senesc
    dsl4jsb_Real2D_onChunk :: MaintRes
    dsl4jsb_Real2D_onChunk :: MaintRes_sum
    dsl4jsb_Real2D_onChunk :: MaintRes_yDayMean
    dsl4jsb_Real2D_onChunk :: temp_sum
    dsl4jsb_Real2D_onChunk :: temp_yDayMean
    dsl4jsb_Real2D_onChunk :: GrowthResp    
    dsl4jsb_Real2D_onChunk :: AutoRes
    dsl4jsb_Real2D_onChunk :: root_maxmoist
    dsl4jsb_Real2D_onChunk :: seed_maxmoist
    dsl4jsb_Real2D_onChunk :: jedi_pot_supply_transp
    dsl4jsb_Real2D_onChunk :: jedi_root_depth
    dsl4jsb_Real2D_onChunk :: soil_depth
    dsl4jsb_Real3D_onChunk :: jedi_root_depth_sl
    dsl4jsb_Real3D_onChunk :: vol_field_cap_sl
    dsl4jsb_Real3D_onChunk :: soil_depth_sl
    dsl4jsb_Real3D_onChunk :: w_soil_sl
    dsl4jsb_Real3D_onChunk :: seed_depth_sl
    dsl4jsb_Real2D_onChunk :: C_balan
    dsl4jsb_Real2D_onChunk :: C_delta
    dsl4jsb_Real2D_onChunk :: c_total
    dsl4jsb_Real2D_onChunk :: one_counter
    ! Taking relative soil moisture from the HYDRO process            
    dsl4jsb_Real2D_onChunk :: alpha_jedi
    dsl4jsb_Real2D_onChunk :: alpha_jedi_seed
    dsl4jsb_Real2D_onChunk :: max_moist
    dsl4jsb_Real2D_onChunk :: w_soil_column_seed
    dsl4jsb_Real2D_onChunk :: w_soil_column_jedi
    dsl4jsb_Real2D_onChunk :: ff_T_30d
    dsl4jsb_Real2D_onChunk :: ff_dTdt        
    dsl4jsb_Real2D_onChunk :: pseudo_c_total              
    dsl4jsb_Real2D_onChunk :: ff_T_30d_yDayMean
    dsl4jsb_Real2D_onChunk :: faPAR_prev_noon
    ! =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- 
    
    dsl4jsb_Real2D_onChunk :: gross_assimilation_ca  ! Gross assimilation of all canopy layers together
    dsl4jsb_Real2D_onChunk :: GPP_sum                ! used to accumulated GPP-Rate over a day
    dsl4jsb_Real2D_onChunk :: NPP_act_yDayMean       ! mean value of actual NPP-Rate yesterday, i.e. after N-limitation
                                                     ! Actual NPP after N-limitation and excess carbon drop.
                                                     ! [mol(CO2)/(m^2(canopy) s)]
    dsl4jsb_Real2D_onChunk :: GPP_yDayMean

    dsl4jsb_Real2D_onChunk :: veg_fract_correction
    dsl4jsb_Real2D_onChunk :: fract_fpc_max
    dsl4jsb_Real2D_onChunk :: t_air                  ! Atmosphere temperature (lowest layer) in Kelvin!
    dsl4jsb_Real2D_onChunk :: rain
    dsl4jsb_Real2D_onChunk :: snow


    !PHS: root depth
    soil_w => Get_vgrid('soil_depth_water')
    
    ! If process is not active on this tile, do nothing
    IF (.NOT. tile%Is_process_active(JEDI_)) RETURN

    iblk  = options%iblk
    ics   = options%ics
    ice   = options%ice
    nc    = options%nc
    dtime = options%dtime

    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    model => Get_model(tile%owner_model_id)

    dsl4jsb_Get_memory(ASSIMI_)
    dsl4jsb_Get_config(PHENO_)
    dsl4jsb_Get_memory(PHENO_)
    dsl4jsb_Get_memory(JEDI_)
    dsl4jsb_Get_config(JEDI_)
    dsl4jsb_Get_memory(A2L_)
    ! =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-  
    ! JeDi needs soil moisture from HYDRO process
    dsl4jsb_Get_memory(HYDRO_)
    dsl4jsb_Get_config(HYDRO_)
    dsl4jsb_Get_memory(RAD_)
    ! =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-  
    ! Variavles for JeDi carbon allocation 
    !
    dsl4jsb_Get_var2D_onChunk(JEDI_, c_storage         )     !     Storage carbon pool 
    dsl4jsb_Get_var2D_onChunk(JEDI_, c_leaves          )     !      Leaves carbon pool             
    dsl4jsb_Get_var2D_onChunk(JEDI_, c_fnroot          )     !   Fine root carbon pool
    dsl4jsb_Get_var2D_onChunk(JEDI_, c_csroot          )     ! Coarse root carbon pool
    dsl4jsb_Get_var2D_onChunk(JEDI_, c_seed            )     !        Seed carbon pool
    dsl4jsb_Get_var2D_onChunk(JEDI_, c_stem            )     !        Stem carbon pool
    dsl4jsb_Get_var2D_onChunk(JEDI_, ff_grow_W         )
    dsl4jsb_Get_var2D_onChunk(JEDI_, ff_grow_T         )
    dsl4jsb_Get_var2D_onChunk(JEDI_, ff_grow_G         )
    dsl4jsb_Get_var2D_onChunk(JEDI_, Flg_grow          )
    dsl4jsb_Get_var2D_onChunk(JEDI_, Flg_grow_G        )
    dsl4jsb_Get_var2D_onChunk(JEDI_, ff_npp            )
    dsl4jsb_Get_var2D_onChunk(JEDI_, flx_germ          )
    dsl4jsb_Get_var2D_onChunk(JEDI_, adj_fraction      )
    dsl4jsb_Get_var2D_onChunk(JEDI_, jedi_light_constraint )
    dsl4jsb_Get_var2D_onChunk(JEDI_, Flg_die           )
    dsl4jsb_Get_var2D_onChunk(JEDI_, Flg_senesc        )
    dsl4jsb_Get_var2D_onChunk(JEDI_, alpha_jedi        )
    
    dsl4jsb_Get_var2D_onChunk(JEDI_, alpha_jedi_seed)
    dsl4jsb_Get_var3D_onChunk(JEDI_, seed_depth_sl)
    dsl4jsb_Get_var2D_onChunk(JEDI_, seed_maxmoist)
    dsl4jsb_Get_var2D_onChunk(JEDI_, w_soil_column_seed)
    dsl4jsb_Get_var2D_onChunk(JEDI_, w_soil_column_jedi)
    dsl4jsb_Get_var2D_onChunk(JEDI_, root_maxmoist )
    dsl4jsb_Get_var2D_onChunk(JEDI_, jedi_pot_supply_transp )
    dsl4jsb_Get_var2D_onChunk(HYDRO_, max_moist         )

    dsl4jsb_Get_var2D_onChunk(JEDI_, jedi_root_depth  )
    dsl4jsb_Get_var2D_onChunk(HYDRO_,  soil_depth     )
    ! 3 Dimensional--------------------------------------------------------------------------
    dsl4jsb_Get_var3D_onChunk(HYDRO_,  w_soil_sl        )
    dsl4jsb_Get_var3D_onChunk(JEDI_, jedi_root_depth_sl ) 
    dsl4jsb_Get_var3D_onChunk(HYDRO_,  vol_field_cap_sl ) 
    dsl4jsb_Get_var3D_onChunk(HYDRO_,  soil_depth_sl    )
    ! 3 Dimensional--------------------------------------------------------------------------
    dsl4jsb_Get_var2D_onChunk(JEDI_, one_counter       )
    dsl4jsb_Get_var2D_onChunk(JEDI_, MaintRes          )
    dsl4jsb_Get_var2D_onChunk(JEDI_, MaintRes_sum      )
    dsl4jsb_Get_var2D_onChunk(JEDI_, MaintRes_yDayMean )
    dsl4jsb_Get_var2D_onChunk(JEDI_, temp_sum          )
    dsl4jsb_Get_var2D_onChunk(JEDI_, temp_yDayMean     )
    dsl4jsb_Get_var2D_onChunk(JEDI_, GrowthResp        )
    dsl4jsb_Get_var2D_onChunk(JEDI_, AutoRes           )
    dsl4jsb_Get_var2D_onChunk(JEDI_, C_balan           )
    dsl4jsb_Get_var2D_onChunk(JEDI_, C_delta           )
    dsl4jsb_Get_var2D_onChunk(JEDI_, c_total           )
    dsl4jsb_Get_var2D_onChunk(JEDI_, ff_dTdt           )
    dsl4jsb_Get_var2D_onChunk(JEDI_, ff_T_30d          )
    dsl4jsb_Get_var2D_onChunk(JEDI_, pseudo_c_total    )
    dsl4jsb_Get_var2D_onChunk(JEDI_, ff_T_30d_yDayMean ) 
    dsl4jsb_Get_var2D_onChunk(RAD_,  faPAR_prev_noon   )

    ! =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- 
   
    ! xx
    dsl4jsb_Get_var2D_onChunk(ASSIMI_,  gross_assimilation_ca )     ! in
    dsl4jsb_Get_var2D_onChunk(JEDI_,    GPP_sum )                   ! inout
    dsl4jsb_Get_var2D_onChunk(JEDI_,    NPP_act_yDayMean )          ! inout
    dsl4jsb_Get_var2D_onChunk(JEDI_,    GPP_yDayMean )              ! inout
    !
    dsl4jsb_Get_var2D_onChunk(PHENO_,   veg_fract_correction )      ! in
    dsl4jsb_Get_var2D_onChunk(PHENO_,   fract_fpc_max )             ! in
    dsl4jsb_Get_var2D_onChunk(A2L_,     t_air)                      ! in
    dsl4jsb_Get_var2D_onChunk(A2L_,     rain)                       ! in
    dsl4jsb_Get_var2D_onChunk(A2L_,     snow)                       ! in

    ! ---------------------------
    ! Go

    ! Put a comment on the screen...
    IF (debug_on() .AND. iblk==1 .AND. is_newday(options%current_datetime, dtime) ) CALL message(TRIM(routine), &
                                                           'update_tradeoff: First time step of new day')
    IF (debug_on() .AND. iblk==1 .AND. is_newyear(options%current_datetime, dtime) ) CALL message(TRIM(routine), &
                                                           'update_tradeoff: First time step of new year')
    
    inidate => get_time_experiment_start()
    iniyear = get_year(inidate)
    CALL deallocateDatetime(inidate)
    lspinup_jedi = dsl4jsb_Config(JEDI_)%l_spinup_jedi
    
    IF ((get_year(options%current_datetime) .LE. iniyear + 2 ) .AND. lspinup_jedi ) THEN    !  PHS: initial carbon
      ff_npp = 0._wp
      Flg_senesc = 0._wp
      Flg_die = 0._wp
      one_counter = 0._wp           
    END IF     !  PHS: initial carbon
    
    seed_depth_sl = jedi_soil_depth_to_layers_1d(  SPREAD(0.05_wp, NCOPIES=nc, DIM=1)  , soil_w%dz(:))         ! Jedi assumes germination of seed depends on the first 50mm of soil wetness
    CALL get_water_in_root_zone(vol_field_cap_sl(:,:) * soil_depth_sl(:,:),  soil_depth_sl(:,:), seed_depth_sl(:,:),     seed_maxmoist(:))    
    CALL get_water_in_root_zone(vol_field_cap_sl(:,:) * soil_depth_sl(:,:), soil_depth_sl(:,:), jedi_root_depth_sl(:,:), root_maxmoist(:))

    CALL get_water_in_root_zone(w_soil_sl(:,:), soil_depth_sl(:,:), seed_depth_sl(:,:), w_soil_column_seed(:))
    CALL get_water_in_root_zone(w_soil_sl(:,:), soil_depth_sl(:,:), jedi_root_depth_sl(:,:), w_soil_column_jedi(:))
      
    alpha_jedi(:)          = jedi_water_stress_factor( w_soil_column_jedi(:),      root_maxmoist(:),  dsl4jsb_Config(HYDRO_)%w_soil_wilt_fract)
    alpha_jedi_seed(:) = jedi_water_stress_factor( w_soil_column_seed(:), seed_maxmoist(:),  dsl4jsb_Config(HYDRO_)%w_soil_wilt_fract)            
    
    ! 2019.11.27 PHS: still at the stage of testing "water stress" 
    WHERE (c_leaves .gt. 0._wp)  
      ! The potential supply rate from soil depends on the root-shoot-ratio. The root-shoot-ratio is the ratio betweem fine-root-carbon and Leaf-area-index.
      !     a(r)=1-exp(- (-log(0.01))/r_optimal   * r  )     a : is a function of root-shoot-ratio  [0...1]
      ! r : root-shoot-ratio [10^-3  ~  10^3]
      ! r_optimal: the median value of the root-shoot-ratio among all PGS
      jedi_pot_supply_transp = (1._wp - exp( -(- log(0.01))/rootshoot_opt * &
        &                         (c_fnroot / (c_leaves*dsl4jsb_Lctlib_param(JeDiSLA)) ) )) * max_pipe_transpir * alpha_jedi
    ELSEWHERE 
      jedi_pot_supply_transp =  0._wp
    END WHERE
    
    ff_grow_T =  ( t_air - tmelt ) / dsl4jsb_Lctlib_param(jedi_N_pseudo_T)      + &
      &          ff_grow_T * dsl4jsb_Lctlib_param(jedi_F_pseudo_T)
    ff_grow_G =  alpha_jedi_seed   / dsl4jsb_Lctlib_param(jedi_N_pseudo_soilM)  + &
      &          ff_grow_G * dsl4jsb_Lctlib_param(jedi_F_pseudo_soilM)
    ff_grow_W =  alpha_jedi        / dsl4jsb_Lctlib_param (jedi_N_pseudo_soilM) + &
      &          ff_grow_W * dsl4jsb_Lctlib_param(jedi_F_pseudo_soilM)
    ff_npp    =  NPP_act_yDayMean  / dsl4jsb_Lctlib_param(jedi_N_pseudo_npp)    + &
      &          ff_npp    * dsl4jsb_Lctlib_param(jedi_F_pseudo_npp)
      
    ff_T_30d  =  ( t_air - tmelt ) / ( 1._wp / (1._wp-EXP(- dtime / (30._wp* 86400._wp )))) + &
      &          ff_T_30d * EXP(- dtime / (30._wp* 86400._wp ))

    !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=!     

    ! All per tile variables need to be re-calculated in case they are written to output
    CALL recalc_per_tile_vars(tile, options,                                      &
      & NPP_act_yDayMean = NPP_act_yDayMean(:), GPP_yDayMean = GPP_yDayMean(:))

    ! Perform daily sums if we ...
    IF ( (            is_time_experiment_start(options%current_datetime) & ! ... just started a totally new experiment
      &   .AND. .NOT. is_time_restart(options%current_datetime) )        & ! ... that is not a restart
      & .OR. .NOT. is_newday(options%current_datetime, dtime)            & ! ... or not starting a new day
      & ) THEN

      GPP_sum  = GPP_sum  + gross_assimilation_ca
      temp_sum = temp_sum + t_air

    ELSE                                    ! A NEW DAY BEGINS and we DID NOT START AN
                                            ! TOTALLY NEW EXPERIMENT ==> perform carbon balance

      ! First save as "previous previous days" and "previous days" means

      GPP_yDayMean = GPP_sum / timesteps_per_day(dtime)

      ! Then restart summing of this days values
      GPP_sum      = gross_assimilation_ca

      !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=!     
      ! PHS: Compute Maintenance respration for JeDi. 
      ! PHS try to compute maintenance respiration only once per day

      temp_yDayMean  = temp_sum / timesteps_per_day(dtime)
       
      MaintRes = 0._wp
      MaintRes = c_sapwood * c_CN_woody * ( c_stem + c_csroot ) + &
        &        dsl4jsb_Lctlib_param(NitroStatus) * dsl4jsb_Lctlib_param(JeDiSLA) * &
        &             ( c_leaves + c_fnroot )  !NNN
      MaintRes = MaintRes  * c_Res_N    * ( 1.6_wp ** ( ( temp_yDayMean - tmelt - 20._wp ) / 10._wp) ) ! JeDi-DGM uses 1.6 for Q10 value
    
      MaintRes = MaintRes / sec_per_day ! The maintanence respiration is daily computed in JeDi.
      MaintRes_yDayMean   = MaintRes
       
      temp_sum = t_air
       
      !PHS 05.08.2019 Tend
      !PHS  I think I made mistake for calculating ff_dTdt 
      !ff_dTdt            = ( ff_T_30d - ff_T_30d_yDayMean ) / ( 1._wp /  ( 1._wp - EXP(- dtime / (30._wp* 86400._wp ))  )) + &
      !  &                  EXP(- dtime / (30._wp* 86400._wp ))  * ff_dTdt  
      ff_dTdt            = ( ff_T_30d - ff_T_30d_yDayMean ) / ( 1._wp /  ( 1._wp - EXP(- 86400._wp / (30._wp* 86400._wp ))  )) + &
        &                  EXP(- 86400._wp / (30._wp* 86400._wp ))  * ff_dTdt   
      ff_T_30d_yDayMean = ff_T_30d

      IF ( .NOT. lspinup_jedi .OR. get_year(options%current_datetime) > iniyear + 2 ) THEN
        CALL calc_Ctradeoff( &
          & GPP_yDayMean,                                & ! in
          ! JeDi plant trait parameter comes in now <<<<<<<<<<<
          & dsl4jsb_Lctlib_param(GrowRspdMoist),         & ! in
          & dsl4jsb_Lctlib_param(GrowRspdTemp ),         & ! in
          & dsl4jsb_Lctlib_param(CritiGrowTemp),         & ! in
          & dsl4jsb_Lctlib_param(EndGrowTemp),           & ! in
          & dsl4jsb_Lctlib_param(GermFrac     ),         & ! in
          & dsl4jsb_Lctlib_param(Allo2Reproduc),         & ! in
          & dsl4jsb_Lctlib_param(Allo2Abground),         & ! in
          & dsl4jsb_Lctlib_param(Allo2Blground),         & ! in
          & dsl4jsb_Lctlib_param(Allo2Storage ),         & ! in
          & dsl4jsb_Lctlib_param(RelAllo2Abgrd),         & ! in
          & dsl4jsb_Lctlib_param(RelAllo2Blgrd),         & ! in
          & dsl4jsb_Lctlib_param(TurnoverStruc),         & ! in
          & dsl4jsb_Lctlib_param(TurnoverLeaf ),         & ! in
          & dsl4jsb_Lctlib_param(SenesRspdNPP ),         & ! in
          & dsl4jsb_Lctlib_param(NitroStatus  ),         & ! in
          & dsl4jsb_Lctlib_param(JeDiSLA),               & ! in
          & dsl4jsb_Lctlib_param(LandcoverClass),        & ! in
          ! JeDi plant trait parameter comes in now >>>>>>>>>>>
          & lspinup_jedi,                                & ! in 
          & iniyear,                                     & ! in
          & get_year(options%current_datetime),          & ! in
          & timesteps_per_day(dtime),                    & ! in
          & faPAR_prev_noon,                             & ! in
          & one_counter,                                 & ! inout
          & c_storage,                                   & ! inout
          & c_leaves,                                    & ! inout
          & c_fnroot,                                    & ! inout
          & c_csroot,                                    & ! inout
          & c_seed,                                      & ! inout
          & c_stem,                                      & ! inout
          & ff_grow_W,                                   & ! in
          & ff_grow_T,                                   & ! in
          & ff_grow_G,                                   & ! in
          & ff_dTdt,                                     & ! in
          & Flg_grow,                                    & ! out
          & Flg_grow_G ,                                 & ! out
          & ff_npp,                                      & ! inout
          & Flg_die,                                     & ! inout
          & Flg_senesc,                                  & ! out
          & flx_germ,                                    & ! out
          & adj_fraction,                                & ! out
          & jedi_light_constraint,                       & ! out
          & MaintRes_yDayMean,                           & ! in
          & GrowthResp,                                  & ! out
          & AutoRes,                                     & ! out
          & jedi_root_depth,                             & ! out 
          & C_balan,                                     & ! out
          & C_delta,                                     & ! out
          & c_total,                                     & !
          & NPP_act_yDayMean                             & ! out
          & )
      END IF  !  PHS: initial carbon

      jedi_root_depth = MAX(jedi_root_depth, 0.05_wp)
      
      jedi_root_depth = MIN(jedi_root_depth, soil_depth)
      
      jedi_root_depth_sl = jedi_soil_depth_to_layers_1d( jedi_root_depth, soil_w%dz(:) )   

      pseudo_c_total  =  c_total / ( 1._wp / (1._wp-EXP(- 86400._wp / (3650._wp* 86400._wp )))) + &
        &                            pseudo_c_total * EXP(- 86400._wp / (3650._wp* 86400._wp ))

      ! Compute carbon contents from canopy to tile box area by weighting pools with fractions of grid box covered by vegetation
      ! ------------------------------------------------------------------------------------------------------------------------

      CALL recalc_per_tile_vars(tile, options,                                      &
        ! =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-  
        ! Variavles for JeDi carbon allocation 
        & c_leaves = c_leaves(:) , c_storage = c_storage(:) , c_fnroot = c_fnroot(:) , &
        & c_csroot = c_csroot(:) , c_stem = c_stem(:)       , c_seed = c_seed(:)     , &
        & c_total = c_total(:),         &
        ! =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-  
        & NPP_act_yDayMean = NPP_act_yDayMean(:), GPP_yDayMean = GPP_yDayMean(:))

    END IF

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')    

  END SUBROUTINE update_tradeoff

! =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- 
      !PHS 24.07.2020 Compute water stress in CARBON at PFT level. 
  
! -------------------------------------------------------------------------------------------------------
  !>
  !! @brief Implementation of "update" for task "water_stress"
  !! Task "update_water_stress" calculates water stress
  !!
  !! @param[in,out] tile    Tile for which routine is executed.
  !! @param[in]     options Additional run-time parameters.
  !!
  SUBROUTINE update_water_stress(tile, options)
    ! ---------------------------
    ! Variables

    ! Used variables
    USE mo_jedi_constants,  ONLY: max_pipe_transpir

    ! Arguments
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    TYPE(t_jsb_model), POINTER :: model
    dsl4jsb_Def_config(HYDRO_)

    dsl4jsb_Def_memory(HYDRO_)
    dsl4jsb_Def_memory(JEDI_)   !PHS

 
    dsl4jsb_Real2D_onChunk :: water_stress
    dsl4jsb_Real2D_onChunk :: evapopot
    dsl4jsb_Real2D_onChunk :: jedi_pot_supply_transp  !PHS
    dsl4jsb_Real2D_onChunk :: alpha_jedi   !PHS	
   
    INTEGER :: iblk, ics, ice
    
    CHARACTER(len=*), PARAMETER :: routine = modname//':update_water_stress'

    iblk = options%iblk
    ics  = options%ics
    ice  = options%ice

    IF (.NOT. tile%Is_process_active(JEDI_) .OR. .NOT. tile%is_vegetation) RETURN

    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    model => Get_model(tile%owner_model_id)

    dsl4jsb_Get_config(HYDRO_)

    ! Get reference to variables for current block
    !
    dsl4jsb_Get_memory(HYDRO_)
    dsl4jsb_Get_memory(JEDI_)

    dsl4jsb_Get_var2D_onChunk(JEDI_,  water_stress)   ! out   !PHS
    
    dsl4jsb_Get_var2D_onChunk(HYDRO_, evapopot)               !PHS
    dsl4jsb_Get_var2D_onChunk(JEDI_,  jedi_pot_supply_transp) !PHS	
    dsl4jsb_Get_var2D_onChunk(JEDI_,  alpha_jedi)             !PHS

    WHERE ( evapopot(:) >= 0._wp )  
      water_stress(:) = 0._wp
    ELSEWHERE ( evapopot(:) < 0._wp )   
      water_stress(:) = MIN(1._wp , ABS(evapopot(:)) / max_pipe_transpir  * &
        &              (1._wp - exp( - jedi_pot_supply_transp(:) / ABS(evapopot(:)) ) ))
    END WHERE
      
    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE update_water_stress
! =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- 

  ! -------------------------------------------------------------------------------------------------------
  !>
  !! Implementation of "aggregate" for task "tradeoff"
  !!
  !! @param[in,out] tile    Tile for which aggregation of child tiles is executed.
  !! @param[in]     config  Vector of process configurations.
  !! @param[in]     options Additional run-time parameters.
  !!
  SUBROUTINE aggregate_tradeoff(tile, options)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    dsl4jsb_Def_memory(JEDI_)

    CLASS(t_jsb_aggregator), POINTER :: weighted_by_fract

    CHARACTER(len=*), PARAMETER :: routine = modname//':aggregate_tradeoff'

    INTEGER  :: iblk , ics, ice

    iblk = options%iblk
    ics  = options%ics
    ice  = options%ice

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    dsl4jsb_Get_memory(JEDI_)

    weighted_by_fract => tile%Get_aggregator("weighted_by_fract")

    !@todo Those commented out are currently not calculated

    ! Variables for JeDi carbon allocation 
    dsl4jsb_Aggregate_onChunk(JEDI_, c_leaves_ta,        weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(JEDI_, c_stem_ta,          weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(JEDI_, c_seed_ta,          weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(JEDI_, c_storage_ta,       weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(JEDI_, c_fnroot_ta,        weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(JEDI_, c_csroot_ta,        weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(JEDI_, c_total_ta,         weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(JEDI_, MaintRes_yDayMean,  weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(JEDI_, GrowthResp,         weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(JEDI_, AutoRes,            weighted_by_fract)
    !dsl4jsb_Aggregate_onChunk(JEDI_, root_maxmoist,      weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(JEDI_, jedi_root_depth,    weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(JEDI_, jedi_root_depth_sl, weighted_by_fract)
    
    ! fluxes
    dsl4jsb_Aggregate_onChunk(JEDI_, NPP_act_yDayMean_ta, weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(JEDI_, GPP_yDayMean_ta,     weighted_by_fract)

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE aggregate_tradeoff

  ! =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- 
      !PHS 24.07.2020 Compute water stress in CARBON at PFT level. 
  
  ! -------------------------------------------------------------------------------------------------------
  !>
  !! Implementation of "aggregate" for task "water_stress"
  !!
  !! @param[in,out] tile    Tile for which routine is executed.
  !! @param[in]     options Additional run-time parameters.
  !!
  SUBROUTINE aggregate_water_stress(tile, options)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    dsl4jsb_Def_memory(JEDI_)
    
    TYPE(t_jsb_model),       POINTER :: model
    CLASS(t_jsb_aggregator), POINTER :: weighted_by_fract

    CHARACTER(len=*), PARAMETER :: routine = modname//':aggregate_water_stress'

    INTEGER :: iblk, ics, ice

    iblk = options%iblk
    ics  = options%ics
    ice  = options%ice

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    model => Get_model(tile%owner_model_id)

    dsl4jsb_Get_memory(JEDI_)
  
    weighted_by_fract => tile%Get_aggregator("weighted_by_fract")

    dsl4jsb_Aggregate_onChunk(JEDI_, water_stress, weighted_by_fract)

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE aggregate_water_stress
  ! =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- 
  

  !##########################################################################################################################
  ! other subroutines (no further tasks)
  !##########################################################################################################################

  ! ================================================================================================================================
  !>
  !! @brief Calculates the per tile states from current canopy states
  !!
  !! @todo discuss how to make more general, e.g. get the info on which state variables to work on from somewhere else?
  !!
  !! @param[in,out] tile    Tile for which routine is executed.
  !! @param[in]     options Additional run-time parameters.
  !
  SUBROUTINE recalc_per_tile_vars(tile, options,                  &
    ! =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-  
    ! Variavles for JeDi carbon allocation 
    & c_leaves , c_storage , c_fnroot , c_csroot , c_stem , c_seed , c_total, &
    ! =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-  
    & NPP_act_yDayMean,          &
    & GPP_yDayMean)

    USE mo_jedi_process,  ONLY: get_per_tile

    ! Arguments
    CLASS(t_jsb_tile_abstract), INTENT(inout), TARGET :: tile
    TYPE(t_jsb_task_options),   INTENT(in)            :: options
    ! =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-  
    ! Variavles for JeDi carbon allocation 
    
    REAL(wp),                   INTENT(in), OPTIONAL  :: c_leaves(:)
    REAL(wp),                   INTENT(in), OPTIONAL  :: c_stem(:)
    REAL(wp),                   INTENT(in), OPTIONAL  :: c_fnroot(:)
    REAL(wp),                   INTENT(in), OPTIONAL  :: c_csroot(:)
    REAL(wp),                   INTENT(in), OPTIONAL  :: c_seed(:)
    REAL(wp),                   INTENT(in), OPTIONAL  :: c_storage(:)
    REAL(wp),                   INTENT(in), OPTIONAL  :: c_total(:)
    ! =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-  

    REAL(wp),                   INTENT(in), OPTIONAL  :: NPP_act_yDayMean(:)
    REAL(wp),                   INTENT(in), OPTIONAL  :: GPP_yDayMean(:)

    TYPE(t_jsb_model), POINTER :: model

    dsl4jsb_Def_memory(JEDI_)
    dsl4jsb_Def_memory(PHENO_)

    ! Local variables
    CHARACTER(len=*), PARAMETER :: routine = modname//':recalc_per_tile_vars'
    INTEGER                     :: ics, ice, iblk

    ! =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-  
    ! Variavles for JeDi carbon allocation 
    
    dsl4jsb_Real2D_onChunk :: c_storage_ta
    dsl4jsb_Real2D_onChunk :: c_seed_ta
    dsl4jsb_Real2D_onChunk :: c_stem_ta
    dsl4jsb_Real2D_onChunk :: c_leaves_ta
    dsl4jsb_Real2D_onChunk :: c_fnroot_ta
    dsl4jsb_Real2D_onChunk :: c_csroot_ta        
    dsl4jsb_Real2D_onChunk :: c_total_ta          
    ! =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-  

    dsl4jsb_Real2D_onChunk :: NPP_act_yDayMean_ta
    dsl4jsb_Real2D_onChunk :: GPP_yDayMean_ta

    dsl4jsb_Real2D_onChunk :: veg_fract_correction
    dsl4jsb_Real2D_onChunk :: fract_fpc_max

    model => Get_model(tile%owner_model_id)

    dsl4jsb_Get_memory(JEDI_)
    dsl4jsb_Get_memory(PHENO_)

    iblk  = options%iblk
    ics   = options%ics
    ice   = options%ice

    ! =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-  
    ! Variavles for JeDi carbon allocation 
    
    dsl4jsb_Get_var2D_onChunk(JEDI_,  c_leaves_ta  )
    dsl4jsb_Get_var2D_onChunk(JEDI_,  c_fnroot_ta  )
    dsl4jsb_Get_var2D_onChunk(JEDI_,  c_csroot_ta  )
    dsl4jsb_Get_var2D_onChunk(JEDI_,  c_stem_ta    )
    dsl4jsb_Get_var2D_onChunk(JEDI_,  c_storage_ta )
    dsl4jsb_Get_var2D_onChunk(JEDI_,  c_seed_ta    )
    dsl4jsb_Get_var2D_onChunk(JEDI_,  c_total_ta   )
    ! =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-  

    dsl4jsb_Get_var2D_onChunk(JEDI_,  NPP_act_yDayMean_ta )       ! opt out
    dsl4jsb_Get_var2D_onChunk(JEDI_,  GPP_yDayMean_ta )           ! opt out

    dsl4jsb_Get_var2D_onChunk(PHENO_, veg_fract_correction )      ! in
    dsl4jsb_Get_var2D_onChunk(PHENO_, fract_fpc_max )             ! in

    ! =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-  
    ! Variavles for JeDi carbon allocation 
    
    IF (PRESENT(c_leaves) )  CALL get_per_tile(c_leaves_ta(:) , c_leaves(:) , veg_fract_correction(:), fract_fpc_max(:))    
    IF (PRESENT(c_stem)   )  CALL get_per_tile(c_stem_ta(:)   , c_stem(:)   , veg_fract_correction(:), fract_fpc_max(:))    
    IF (PRESENT(c_fnroot) )  CALL get_per_tile(c_fnroot_ta(:) , c_fnroot(:) , veg_fract_correction(:), fract_fpc_max(:))    
    IF (PRESENT(c_csroot) )  CALL get_per_tile(c_csroot_ta(:) , c_csroot(:) , veg_fract_correction(:), fract_fpc_max(:))
    IF (PRESENT(c_seed)   )  CALL get_per_tile(c_seed_ta(:)   , c_seed(:)   , veg_fract_correction(:), fract_fpc_max(:))   
    IF (PRESENT(c_storage))  CALL get_per_tile(c_storage_ta(:), c_storage(:), veg_fract_correction(:), fract_fpc_max(:))
    IF (PRESENT(c_total)  )  CALL get_per_tile(c_total_ta(:)  , c_total(:)  , veg_fract_correction(:), fract_fpc_max(:))
    ! =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-  

    ! C fluxes
    IF (PRESENT(NPP_act_yDayMean)) CALL get_per_tile(NPP_act_yDayMean_ta(:), &
      & NPP_act_yDayMean(:), veg_fract_correction(:), fract_fpc_max(:))
    IF (PRESENT(GPP_yDayMean)) CALL get_per_tile(GPP_yDayMean_ta(:), &
      & GPP_yDayMean(:), veg_fract_correction(:), fract_fpc_max(:))

  END SUBROUTINE recalc_per_tile_vars

  SUBROUTINE jedi_biomass_to_fract(tile, options)

    USE mo_jsb_var_class,   ONLY:t_jsb_var_real2d
    USE mo_jsb_grid,           ONLY: Get_grid
    USE mo_jsb_tile_class,  ONLY: t_fract_memory

    CLASS(t_jsb_tile_abstract), INTENT(inout)         :: tile
    TYPE(t_jsb_task_options)  , INTENT(in)         :: options
    
    TYPE(t_jsb_var_real2d), POINTER :: var
    CLASS(t_jsb_tile_abstract), POINTER :: child

    dsl4jsb_Real2D_onDomain :: pseudo_c_total

    dsl4jsb_Def_memory(JEDI_)

    TYPE(t_jsb_model), POINTER      :: model
    TYPE(t_jsb_grid),  POINTER      :: grid

    REAL(wp), ALLOCATABLE :: c_total_sum(:,:), fract(:,:)
    INTEGER               :: nblks, nproma
    INTEGER               :: i, no_children

    CHARACTER(len=*), PARAMETER :: routine = modname//':jedi_biomass_to_fract'

    IF (TRIM(tile%name) /= 'veg') CALL finish(routine, 'should only be called on veg tile')

    dsl4jsb_Get_memory(JEDI_)

    model => Get_model(tile%owner_model_id)
    grid => Get_grid(model%grid_id)

    nblks = grid%nblks
    nproma = grid%nproma

    var => dsl4jsb_memory(JEDI_)%pseudo_c_total

    ALLOCATE(c_total_sum(nproma,nblks))
    c_total_sum(:,:) = 0._wp

    ALLOCATE(fract(nproma,nblks))
    CALL tile%Get_fraction(fract=fract(:,:)) ! Fractions of veg tile

    no_children = tile%Get_no_of_children()

    IF (ASSOCIATED(dsl4jsb_memory(JEDI_)%children(1)%p%vars(var%child_idx(1))%p%ptr2d)) THEN

      DO i=1,no_children
        IF (var%child_idx(i) > 0) THEN
          ASSOCIATE (data_child => dsl4jsb_memory(JEDI_)%children(i)%p%vars(var%child_idx(i))%p%ptr2d(:,:) )
            c_total_sum(:,:) = c_total_sum(:,:) + data_child(:,:)
          END ASSOCIATE
        END IF
      END DO
      c_total_sum(:,:) = c_total_sum(:,:) * fract(:,:) ! Fractions of children need to sum up to fraction of parent tile in block below

      child => tile%Get_first_child_tile()
      DO i=1,no_children
        SELECT TYPE (child__mem => child%mem(0)%p)
        TYPE IS (t_fract_memory)
          ASSOCIATE (data_child => dsl4jsb_memory(JEDI_)%children(i)%p%vars(var%child_idx(i))%p%ptr2d(:,:) )
            WHERE (c_total_sum(:,:) > 0._wp)
              child__mem%fract%ptr(:,:) = data_child(:,:) / c_total_sum(:,:)
            END WHERE
          END ASSOCIATE
          child => child%Get_next_sibling_tile()
        END SELECT
      END DO

    ELSE IF (ASSOCIATED(var%ptr2d_t)) THEN

      IF (SIZE(var%ptr_t, 2) /= no_children) CALL finish(routine, 'Dimension mismatch in c_total%ptr_t and no_children')

      c_total_sum(:,:) = SUM(var%ptr2d_t(:,:,:), DIM=2)
      child => tile%Get_first_child_tile()
      DO i=1,no_children
        SELECT TYPE (child__mem => child%mem(0)%p)
        TYPE IS (t_fract_memory)
          WHERE (c_total_sum(:,:) > 0._wp)
            child__mem%fract%ptr(:,:) = var%ptr2d_t(:,i,:) / c_total_sum(:,:) * fract(:,:) ! Fractions of children need to sum up to fraction of parent tile
          END WHERE
          child => child%Get_next_sibling_tile()
        END SELECT
      END DO

    ELSE

      CALL finish(routine, 'Could not get children of variable c_total')

    END IF
    
    tile%l_fract_children_changed(:,:) = .TRUE.

    DEALLOCATE(c_total_sum)
    DEALLOCATE(fract)

  END SUBROUTINE jedi_biomass_to_fract

#endif
END MODULE mo_jedi_interface
