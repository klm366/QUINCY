!>
!! @brief Contains the memory class for the jedi process.
!!
!! @author
!!  Rainer Schneck, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Rainer Schneck (2017-04-6)
!!
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/license .
!!

MODULE mo_jedi_memory_class
#ifndef __NO_JSBACH__

  USE mo_kind, ONLY: wp
  USE mo_util, ONLY: One_of

  USE mo_jsb_model_class,        ONLY: t_jsb_model
  USE mo_jsb_class,              ONLY: Get_model
  USE mo_jsb_lct_class,          ONLY: VEG_TYPE, LAND_TYPE
  USE mo_jsb_memory_class,       ONLY: t_jsb_memory
  USE mo_jsb_var_class,          ONLY: t_jsb_var_real2d, t_jsb_var_real3d
  USE mo_jsb_io,                 ONLY: TSTEP_CONSTANT
  ! Use of prcesses in this module
  dsl4jsb_Use_processes JEDI_

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: t_jedi_memory, max_no_of_vars

  INTEGER, PARAMETER :: max_no_of_vars = 170   !PHS changed from 110 to 120

  !> Type definition for memory
  TYPE, EXTENDS(t_jsb_memory) :: t_jedi_memory

    ! cbalance_type variables

    TYPE(t_jsb_var_real2d) ::       &
      !
      LAI_sum,                      & ! used to accumulate LAI over a day
      NPP_pot_sum,                  & ! used to accumulated NPP-Rate over a day
      GPP_sum,                      & ! used to accumulated GPP-Rate over a day
      NPP_act_yDayMean,             & ! mean value of actual NPP-Rate after N-limitation and excess carbon drop.
                                      ! "yDayMean" because the actual NPP depends on/is calculated from results of the
                                      ! previous day.
      GPP_yDayMean

    ! Output on tile areas:
    TYPE(t_jsb_var_real2d) ::       &
      !
      NPP_act_yDayMean_ta,          & ! mean value of actual NPP-Rate yesterday, i.e. after N-limitation.
                                      ! Actual NPP (after N-limitation and) excess carbon drop.
                                      ! [mol(CO2)/(m^2(tile area) s)]
      GPP_yDayMean_ta

    ! =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-  
    ! Variavles for JeDi carbon allocation 
   
    TYPE(t_jsb_var_real2d) ::       &
      c_storage,                    &
      c_storage_temp,               &
      c_seed,                       &
      c_stem,                       &
      c_leaves,                     &
      c_fnroot,                     &
      c_csroot,                     &
      c_total                       
   TYPE(t_jsb_var_real2d) ::        &
      c_storage_ta,                 &
      c_seed_ta,                    &
      c_stem_ta,                    &
      c_leaves_ta,                  &
      c_fnroot_ta,                  &
      c_csroot_ta,                  &
      c_total_ta
       
    TYPE(t_jsb_var_real2d) ::       &
      ff_grow_T,                    &
      ff_grow_G,                    &
      ff_grow_W,                    &
      Flg_grow,                     &
      Flg_grow_G,                   &
      ff_npp,                       &
      Flg_die,                      &
      Flg_senesc,                   &
      alpha_jedi,                   &
      MaintRes,                     &
      MaintRes_sum,                 &
      MaintRes_yDayMean,            &
      temp_sum,                     &
      temp_yDayMean,                &
      GrowthResp,                   &
      AutoRes ,                     &
      C_balan,                      &
      C_delta,                      &
      flx_germ,                     &
      one_counter,                  &
      ff_dTdt,                      &                     
      ff_T_30d_yDayMean,            &
      ff_T_30d,                     &
      pseudo_c_total,               &
      root_maxmoist,                &
      jedi_pot_supply_transp,       &
      jedi_root_depth,              &
      adj_fraction,                 &
      jedi_light_constraint,        &
      seed_maxmoist,                &
      alpha_jedi_seed,              &
      w_soil_column_jedi,           &
      w_soil_column_seed,           &
      water_stress
         
    TYPE(t_jsb_var_real3d) ::       &
      seed_depth_sl,                &
      jedi_root_depth_sl           
    ! =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-  


  CONTAINS
    PROCEDURE :: Init => Init_jedi_memory
  END TYPE t_jedi_memory


  CHARACTER(len=*), PARAMETER :: modname = 'mo_jedi_memory_class'

CONTAINS

  SUBROUTINE Init_jedi_memory(mem, prefix, suffix, lct_ids, model_id, l_fixed_fractions)

    USE mo_jsb_varlist,       ONLY: BASIC, MEDIUM, FULL
    USE mo_jsb_io,            ONLY: grib_bits, t_cf, t_grib1, t_grib2, tables
    USE mo_jsb_grid_class,    ONLY: t_jsb_grid, t_jsb_vgrid
    USE mo_jsb_grid,          ONLY: Get_grid, Get_vgrid

    CLASS(t_jedi_memory), INTENT(inout), TARGET :: mem
    CHARACTER(len=*),    INTENT(in)    :: prefix
    CHARACTER(len=*),    INTENT(in)    :: suffix
    INTEGER,             INTENT(in)    :: lct_ids(:)
    INTEGER,             INTENT(in)    :: model_id
    LOGICAL, OPTIONAL,   INTENT(in)    :: l_fixed_fractions

    TYPE(t_jsb_model), POINTER :: model
    TYPE(t_jsb_grid),  POINTER :: hgrid             ! Horizontal grid
    TYPE(t_jsb_vgrid), POINTER :: surface, soil_w   ! Vertical grids
    INTEGER :: table

    CHARACTER(len=*), PARAMETER :: routine = modname//':Init_jedi_memory'

    IF (lct_ids(1) > 0) CONTINUE ! avoid compiler warning about dummy argument not being used

    model => Get_model(model_id)

    table = tables(1)

    hgrid   => Get_grid(mem%grid_id)
    surface => Get_vgrid('surface')
    soil_w => Get_vgrid('soil_depth_water')

    ! =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-  
    ! Variavles for JeDi carbon allocation 
    ! =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-  
    CALL mem%Add_var( 'c_storage', mem%c_storage,                                     &
      & hgrid, surface,                                                                        &
      & t_cf('c_storage', 'mol(C) m-2(canopy)', 'C-Pool for green parts of vegetation.'),    &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                     &
      & prefix, suffix,                                                                        &
      & output_level=BASIC,loutput=.FALSE.,                                                      &
      & lrestart=.TRUE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) 
      
    CALL mem%Add_var( 'c_seed', mem%c_seed,                                     &
      & hgrid, surface,                                                                        &
      & t_cf('c_seed', 'mol(C) m-2(canopy)', 'C-Pool for green parts of vegetation.'),    &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                     &
      & prefix, suffix,                                                                        &
      & output_level=BASIC,loutput=.FALSE.,                                                      &
      & lrestart=.TRUE., initval_r=0.0_wp, l_aggregate_all=.TRUE. )
      
    CALL mem%Add_var( 'c_stem', mem%c_stem,                                     &
      & hgrid, surface,                                                                        &
      & t_cf('c_stem', 'mol(C) m-2(canopy)', 'C-Pool of stem.'),    &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                     &
      & prefix, suffix,                                                                        &
      & output_level=BASIC,loutput=.FALSE.,                                                      &
      & lrestart=.TRUE., initval_r=0.0_wp, l_aggregate_all=.TRUE. )

    CALL mem%Add_var( 'c_csroot', mem%c_csroot,                                     &
      & hgrid, surface,                                                                        &
      & t_cf('c_csroot', 'mol(C) m-2(canopy)', 'C-Pool of coarse root.'),    &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                     &
      & prefix, suffix,                                                                        &
      & output_level=BASIC,loutput=.FALSE.,                                                      &
      & lrestart=.TRUE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) 

    CALL mem%Add_var( 'c_fnroot', mem%c_fnroot,                                     &
      & hgrid, surface,                                                                        &
      & t_cf('c_fnroot', 'mol(C) m-2(canopy)', 'C-Pool of fine root.'),    &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                     &
      & prefix, suffix,                                                                        &
      & output_level=BASIC,loutput=.FALSE.,                                                      &
      & lrestart=.TRUE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) 

    CALL mem%Add_var( 'c_leaves', mem%c_leaves,                                     &
      & hgrid, surface,                                                                        &
      & t_cf('c_leaves', 'mol(C) m-2(canopy)', 'C-Pool of leaves'),    &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                     &
      & prefix, suffix,                                                                        &
      & output_level=BASIC,loutput=.FALSE.,                                                      &
      & lrestart=.TRUE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) 

    CALL mem%Add_var( 'c_total', mem%c_total,                              &
      & hgrid, surface,                                                                        &
      & t_cf('c_total','mol(C) m-2(canopy)' ,'Total carbon of six pools'),           &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                     &
      & prefix, suffix,                                                                        &
      & output_level=BASIC,loutput=.FALSE.,                                                      &
      & lrestart=.TRUE., initval_r=0.0_wp, l_aggregate_all=.TRUE. )   

    CALL mem%Add_var( 'pseudo_c_total', mem%pseudo_c_total,                                    &
      & hgrid, surface,                                                                        &
      & t_cf('pseudo_c_total','mol(C) m-2(canopy)' ,'Pseudo total carbon of six pools'),           &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                     &
      & prefix, suffix,                                                                        &
      & output_level=BASIC,loutput=.TRUE.,                                                      &
      & lrestart=.TRUE., initval_r=0.0_wp, l_aggregate_all=.TRUE. )   

    CALL mem%Add_var( 'c_storage_ta', mem%c_storage_ta,                                           &
      & hgrid, surface,                                                                        &
      & t_cf('c_storage_ta', 'mol(C) m-2', 'C-Pool for green parts of vegetation.'),    &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                     &
      & prefix, suffix,                                                                        &
      & output_level=BASIC,loutput=.TRUE.,                                                      &
      & lrestart=.TRUE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) 
      
    CALL mem%Add_var( 'c_seed_ta', mem%c_seed_ta,                                     &
      & hgrid, surface,                                                                        &
      & t_cf('c_seed_ta', 'mol(C) m-2', 'C-Pool of seed.'),    &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                     &
      & prefix, suffix,                                                                        &
      & output_level=BASIC,loutput=.TRUE.,                                                      &
      & lrestart=.TRUE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) 

    CALL mem%Add_var( 'c_stem_ta', mem%c_stem_ta,                                     &
      & hgrid, surface,                                                                        &
      & t_cf('c_stem_ta', 'mol(C) m-2', 'C-Pool of stem.'),    &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                     &
      & prefix, suffix,                                                                        &
      & output_level=BASIC,loutput=.TRUE.,                                                      &
      & lrestart=.TRUE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) 

    CALL mem%Add_var( 'c_csroot_ta', mem%c_csroot_ta,                                     &
      & hgrid, surface,                                                                        &
      & t_cf('c_csroot_ta', 'mol(C) m-2' ,'C-Pool of coarse root.'),    &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                     &
      & prefix, suffix,                                                                        &
      & output_level=BASIC,loutput=.TRUE.,                                                      &
      & lrestart=.TRUE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) 

    CALL mem%Add_var( 'c_fnroot_ta', mem%c_fnroot_ta,                                     &
      & hgrid, surface,                                                                        &
      & t_cf('c_fnroot_ta','mol(C) m-2', 'C-Pool of fine root.'),    &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                     &
      & prefix, suffix,                                                                        &
      & output_level=BASIC,loutput=.TRUE.,                                                      &
      & lrestart=.TRUE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) 
      
    CALL mem%Add_var( 'c_leaves_ta', mem%c_leaves_ta,                                     &
      & hgrid, surface,                                                                        &
      & t_cf('c_leaves_ta', 'mol(C) m-2', 'C-Pool od leaves.'),    &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                     &
      & prefix, suffix,                                                                        &
      & output_level=BASIC,loutput=.TRUE.,                                                      &
      & lrestart=.TRUE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) 
    
    CALL mem%Add_var( 'c_total_ta', mem%c_total_ta,                              &
      & hgrid, surface,                                                                        &
      & t_cf('c_total_ta','mol(C) m-2','Total carbon of six pools'),           &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                     &
      & prefix, suffix,                                                                        &
      & output_level=BASIC,loutput=.TRUE.,                                                      &
      & lrestart=.TRUE., initval_r=0.0_wp, l_aggregate_all=.TRUE. )   
      
    CALL mem%Add_var( 'ff_grow_T', mem%ff_grow_T,                                      &
      & hgrid, surface,                                                                        &
      & t_cf('ff_grow_T', 'deg C', 'Pseudo air temperature'),    &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                     &
      & prefix, suffix,                                                                        &
      & output_level=BASIC,loutput=.FALSE.,                                                      &
      & lrestart=.TRUE., initval_r=0.0_wp, l_aggregate_all=.FALSE. ) 
        
    CALL mem%Add_var( 'ff_grow_G', mem%ff_grow_G,                                      &
      & hgrid, surface,                                                                        &
      & t_cf('ff_grow_G', '-', 'Pseudo top 50mm soil moisture'),    &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                     &
      & prefix, suffix,                                                                        &
      & output_level=BASIC,loutput=.FALSE.,                                                      &
      & lrestart=.TRUE., initval_r=0.0_wp, l_aggregate_all=.FALSE. )  
      
    CALL mem%Add_var( 'ff_grow_W', mem%ff_grow_W,                                      &
      & hgrid, surface,                                                                        &
      & t_cf('ff_grow_W', '-', 'Pseudo soil moisture to root-depth'),    &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                     &
      & prefix, suffix,                                                                        &
      & output_level=BASIC,loutput=.FALSE.,                                                      &
      & lrestart=.TRUE., initval_r=0.0_wp, l_aggregate_all=.FALSE. )  
    
    CALL mem%Add_var( 'ff_npp', mem%ff_npp,                                              &
      & hgrid, surface,                                                                        &
      & t_cf('ff_npp', 'mol(C) m-2(canopy) s-1', 'Pseudo net primary production'),           &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                     &
      & prefix, suffix,                                                                        &
      & output_level=BASIC,loutput=.TRUE.,                                                      &
      & lrestart=.TRUE., initval_r=0.0_wp, l_aggregate_all=.FALSE. )  
    
    CALL mem%Add_var( 'flx_germ', mem%flx_germ,                              &
      & hgrid, surface,                                                                        &
      & t_cf('flx_germ', 'mol(C) m-2(canopy) s-1', 'carbon flux from germination '),           &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                     &
      & prefix, suffix,                                                                        &
      & output_level=BASIC,loutput=.FALSE.,                                                      &
      & lrestart=.TRUE., initval_r=0.0_wp, l_aggregate_all=.FALSE. )   
    
    CALL mem%Add_var( 'adj_fraction', mem%adj_fraction,                              &
      & hgrid, surface,                                                                        &
      & t_cf('adj_fraction', '-', 'Constraint factor for leaf growth '),           &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                     &
      & prefix, suffix,                                                                        &
      & output_level=BASIC,loutput=.FALSE.,                                                      &
      & lrestart=.FALSE., initval_r=1.0_wp, l_aggregate_all=.FALSE. )   

    CALL mem%Add_var( 'jedi_light_constraint', mem%jedi_light_constraint,                              &
      & hgrid, surface,                                                                        &
      & t_cf('jedi_light_constraint', '-', 'Light constraint factor for leaf growth '),          &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                     &
      & prefix, suffix,                                                                        &
      & output_level=BASIC,loutput=.FALSE.,                                                      &
      & lrestart=.FALSE., initval_r=1.0_wp, l_aggregate_all=.FALSE. )   

    
    CALL mem%Add_var( 'Flg_grow', mem%Flg_grow,                                              &
      & hgrid, surface,                                                                        &
      & t_cf('Flg_grow', '-', 'Flag of growth'),           &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                     &
      & prefix, suffix,                                                                        &
      & output_level=BASIC,loutput=.TRUE.,                                                      &
      & lrestart=.FALSE., initval_r=0.0_wp, l_aggregate_all=.FALSE. )    
    
    CALL mem%Add_var( 'Flg_grow_G', mem%Flg_grow_G,                                              &
      & hgrid, surface,                                                                        &
      & t_cf('Flg_grow_G', '-', 'Flag of germination'),           &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                     &
      & prefix, suffix,                                                                        &
      & output_level=BASIC,loutput=.TRUE.,                                                      &
      & lrestart=.FALSE., initval_r=0.0_wp, l_aggregate_all=.FALSE. )    
    
    CALL mem%Add_var( 'Flg_die', mem%Flg_die,                                              &
      & hgrid, surface,                                                                        &
      & t_cf('Flg_die', '-', 'Dead'),           &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                     &
      & prefix, suffix,                                                                        &
      & output_level=BASIC,loutput=.TRUE.,                                                      &
      & lrestart=.TRUE., initval_r=0.0_wp, l_aggregate_all=.FALSE. )    
    
    CALL mem%Add_var( 'Flg_senesc', mem%Flg_senesc,                                              &
      & hgrid, surface,                                                                        &
      & t_cf('Flg_senesc', '-', 'Senescence flag'),           &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                     &
      & prefix, suffix,                                                                        &
      & output_level=BASIC,loutput=.TRUE.,                                                      &
      & lrestart=.FALSE., initval_r=0.0_wp, l_aggregate_all=.FALSE. )    
      
    CALL mem%Add_var( 'alpha_jedi', mem%alpha_jedi,                              &
      & hgrid, surface,                                                                        &
      & t_cf('alpha_jedi', '-', 'relative soil wetness for JeDi'),           &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                     &
      & prefix, suffix,                                                                        &
      & output_level=BASIC,loutput=.TRUE.,                                                      &
      & lrestart=.TRUE., initval_r=0.0_wp, l_aggregate_all=.FALSE. )   
      
    CALL mem%Add_var( 'alpha_jedi_seed', mem%alpha_jedi_seed,                                  &
      & hgrid, surface,                                                                        &
      & t_cf('alpha_jedi_seed', '-', 'relative soil wetness for seed (upper 50mm of the soil)'),           &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                     &
      & prefix, suffix,                                                                        &
      & output_level=BASIC,loutput=.TRUE.,                                                      &
      & lrestart=.TRUE., initval_r=0.0_wp, l_aggregate_all=.FALSE. )   
      
    CALL mem%Add_var( 'jedi_root_depth', mem%jedi_root_depth,                   &
      &       hgrid, surface,                                                    &
      & t_cf('jedi_root_depth', 'm', 'Jedi root depth'),                                     &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits), &
      & prefix, suffix,&
      & output_level=BASIC, loutput=.TRUE.,                                                  &
      & lrestart=.TRUE.,initval_r=0.0_wp, l_aggregate_all=.FALSE.)!, isteptype=TSTEP_CONSTANT )
      
    CALL mem%Add_var( 'seed_maxmoist', mem%seed_maxmoist,                              &
            & hgrid, surface,                                                                        &
      & t_cf('seed_maxmoist','m' ,'Maximum moisture of the top 50 mm soil (used for seed germination)'),           &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                     &
      & prefix, suffix,                                                                        &
      & output_level=BASIC,loutput=.FALSE.,                                                      &
      & lrestart=.FALSE., initval_r=0.0_wp, l_aggregate_all=.FALSE.)!, isteptype=TSTEP_CONSTANT )
      
    CALL mem%Add_var( 'w_soil_column_seed', mem%w_soil_column_seed,                              &
      & hgrid, surface,                                                                        &
      & t_cf('w_soil_column_seed', 'm water equivalent' ,'cumulated soil water in the top 50mm'),           &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                     &
      & prefix, suffix,                                                                        &
      & output_level=BASIC,loutput=.FALSE.,                                                      &
      & lrestart=.FALSE., initval_r=0.0_wp, l_aggregate_all=.FALSE.)!, isteptype=TSTEP_CONSTANT ) 

    CALL mem%Add_var( 'w_soil_column_jedi', mem%w_soil_column_jedi,                                    &
      & hgrid, surface,                                                                        &
      & t_cf('w_soil_column_jedi', 'm water equivalent' ,'cumulated soil water in jedi root zone'),           &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                     &
      & prefix, suffix,                                                                        &
      & output_level=BASIC,loutput=.FALSE.,                                                      &
      & lrestart=.FALSE., initval_r=0.0_wp, l_aggregate_all=.FALSE.)!, isteptype=TSTEP_CONSTANT ) 
    
    CALL mem%Add_var( 'seed_depth_sl', mem%seed_depth_sl,                              &
      & hgrid, soil_w,                                                                        &
      & t_cf('seed_depth_sl','m', 'depth of the soil layer for the top 50 mm soil '),           &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                     &
      & prefix, suffix,                                                                        &
      & output_level=BASIC,loutput=.FALSE.,                                                      &
      & lrestart=.TRUE., initval_r=0.0_wp, l_aggregate_all=.FALSE.)!, isteptype=TSTEP_CONSTANT )
      
    CALL mem%Add_var( 'jedi_root_depth_sl', mem%jedi_root_depth_sl,                                    &
      & hgrid, soil_w,                                                                        &
      & t_cf('jedi_root_depth_sl','m', ' root depth in each soil layer '),           &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                     &
      & prefix, suffix,                                                                        &
      & output_level=BASIC,loutput=.FALSE.,                                                      &
      & lrestart=.TRUE., initval_r=0.0_wp, l_aggregate_all=.FALSE.)!, isteptype=TSTEP_CONSTANT )

    CALL mem%Add_var( 'MaintRes', mem%MaintRes,                              &
      & hgrid, surface,                                                                        &
      & t_cf('MaintRes', 'mol(C) m-2(canopy) s-1', 'Maintenance Respiration'),           &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                     &
      & prefix, suffix,                                                                        &
      & output_level=BASIC,loutput=.False.,                                                      &
      & lrestart=.FALSE., initval_r=0.0_wp, l_aggregate_all=.FALSE. )  
      
    CALL mem%Add_var( 'MaintRes_sum', mem%MaintRes_sum,                              &
      & hgrid, surface,                                                                        &
      & t_cf('MaintRes_sum', 'mol(C) m-2(canopy) s-1', 'Used to accumulate maintenance Respiration over a day '),           &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                     &
      & prefix, suffix,                                                                        &
      & output_level=BASIC,loutput=.False.,                                                      &
      & lrestart=.FALSE., initval_r=0.0_wp, l_aggregate_all=.FALSE. )      
        
    CALL mem%Add_var( 'MaintRes_ydaymean', mem%MaintRes_yDayMean,                              &
      & hgrid, surface,                                                                        &
      & t_cf('MaintRes_ydaymean', 'mol(C) m-2(canopy) s-1', 'Mean maintenance Respiration of previous day'),           &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                     &
      & prefix, suffix,                                                                        &
      & output_level=BASIC,loutput=.TRUE.,                                                      &
      & lrestart=.TRUE., initval_r=0.0_wp, l_aggregate_all=.FALSE. )    
    
    CALL mem%Add_var( 'temp_sum', mem%temp_sum,                              &
      & hgrid, surface,                                                                        &
      & t_cf('temp_sum', 'K', 'Used to accumulate temperature over a day. '),           &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                     &
      & prefix, suffix,                                                                        &
      & output_level=BASIC,loutput=.FALSE.,                                                      &
      & lrestart=.TRUE., initval_r=0.0_wp, l_aggregate_all=.FALSE. )    
      
    CALL mem%Add_var( 'temp_yDayMean', mem%temp_yDayMean,                              &
      & hgrid, surface,                                                                        &
      & t_cf('temp_yDayMean', 'K', 'Mean temperature of previous day'),           &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                     &
      & prefix, suffix,                                                                        &
      & output_level=BASIC,loutput=.FALSE.,                                                      &
      & lrestart=.FALSE., initval_r=15.0_wp, l_aggregate_all=.FALSE. )    
          
    CALL mem%Add_var( 'GrowthResp', mem%GrowthResp,                              &
      & hgrid, surface,                                                                        &
      & t_cf('GrowthResp', 'mol(C) m-2(canopy) s-1', 'Growth respiration'),           &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                     &
      & prefix, suffix,                                                                        &
      & output_level=BASIC,loutput=.TRUE.,                                                      &
      & lrestart=.FALSE., initval_r=0.0_wp, l_aggregate_all=.FALSE. )  
    
    
    CALL mem%Add_var( 'AutoRes', mem%AutoRes,                              &
      & hgrid, surface,                                                                        &
      & t_cf('AutoRes', 'mol(C) m-2(canopy) s-1', 'Autotrophic respiration'),           &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                     &
      & prefix, suffix,                                                                        &
      & output_level=BASIC,loutput=.TRUE.,                                                      &
      & lrestart=.FALSE., initval_r=0.0_wp, l_aggregate_all=.FALSE. )  
      
      
    CALL mem%Add_var( 'C_balan', mem%C_balan,                              &
      & hgrid, surface,                                                                        &
      & t_cf('C_balan', 'mol(C) m-2(canopy) s-1', 'Carbon balance check'),           &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                     &
      & prefix, suffix,                                                                        &
      & output_level=BASIC,loutput=.TRUE.,                                                      &
      & lrestart=.FALSE., initval_r=0.0_wp, l_aggregate_all=.FALSE. )   
      
    CALL mem%Add_var( 'C_delta', mem%C_delta,                              &
      & hgrid, surface,                                                                        &
      & t_cf('C_delta', 'mol(C) m-2(canopy) s-1', 'Carbon delta check'),           &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                     &
      & prefix, suffix,                                                                        &
      & output_level=BASIC,loutput=.TRUE.,                                                      &
      & lrestart=.FALSE., initval_r=0.0_wp, l_aggregate_all=.FALSE. )   
        
        
    CALL mem%Add_var( 'one_counter', mem%one_counter,                              &
      & hgrid, surface,                                                                        &
      & t_cf('one_counter', '-', 'Counter to contro initial phase of growth'),           &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                     &
      & prefix, suffix,                                                                        &
      & output_level=BASIC,loutput=.TRUE.,                                                      &
      & lrestart=.TRUE., initval_r=0.0_wp, l_aggregate_all=.FALSE. )   
      
      
    CALL mem%Add_var( 'ff_dTdt', mem%ff_dTdt,                              &
      & hgrid, surface,                                                                        &
      & t_cf('ff_dTdt', 'deg C', 'Pseudo air temperature change '),           &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                     &
      & prefix, suffix,                                                                        &
      & output_level=BASIC,loutput=.FALSE.,                                                      &
      & lrestart=.TRUE., initval_r=0.0_wp, l_aggregate_all=.FALSE. )   
      
    CALL mem%Add_var( 'ff_T_30d_yDayMean', mem%ff_T_30d_yDayMean,                              &
      & hgrid, surface,                                                                        &
      & t_cf('ff_T_30d_yDayMean', 'deg C', 'Yesterday mean of pseudo temperature for 30 days'),           &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                     &
      & prefix, suffix,                                                                        &
      & output_level=BASIC,loutput=.FALSE.,                                                      &
      & lrestart=.TRUE., initval_r=0.0_wp, l_aggregate_all=.FALSE. )   
        
    CALL mem%Add_var( 'ff_T_30d', mem%ff_T_30d,                              &
      & hgrid, surface,                                                                        &
      & t_cf('ff_T_30d', 'deg C', 'Pseudo temperature for 30 days '),           &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                     &
      & prefix, suffix,                                                                        &
      & output_level=BASIC,loutput=.FALSE.,                                                      &
      & lrestart=.TRUE., initval_r=0.0_wp, l_aggregate_all=.FALSE. )             
      
    CALL mem%Add_var( 'root_maxmoist', mem%root_maxmoist,                              &
      & hgrid, surface,                                                                        &
      & t_cf('root_maxmoist', 'm', 'jedi root depth max. moist'),           &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                     &
      & prefix, suffix,                                                                        &
      & output_level=BASIC,loutput=.FALSE.,                                                      &
      & lrestart=.FALSE., initval_r=0.0_wp, l_aggregate_all=.FALSE. )    
      
    CALL mem%Add_var( 'jedi_pot_supply_transp', mem%jedi_pot_supply_transp,                              &
      & hgrid, surface,                                                                        &
      & t_cf('jedi_pot_supply_transp', 'kg m-2 s-1', 'fine root potential supply rate for transpiration'),           &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                     &
      & prefix, suffix,                                                                        &
      & output_level=BASIC,loutput=.TRUE.,                                                      &
      & lrestart=.FALSE., initval_r=0.0_wp, l_aggregate_all=.FALSE. )     
           
    CALL mem%Add_var( 'water_stress', mem%water_stress,                &
      & hgrid, surface,                                                     &
      & t_cf('water_stress', '-', 'water stress for canopy'),                     &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),  &        
      & prefix, suffix,                                                      &
      & output_level=BASIC,loutput=.TRUE.,                                  &
      & lrestart=.FALSE., initval_r=0.0_wp )
        
    ! =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-  

    CALL mem%Add_var('gpp_sum', mem%GPP_sum,                                             &
      & hgrid, surface,                                                                  &
      & t_cf('GPP_sum', 'mol(C) m-2(canopy) s-1', 'Used to accumulate GPP over a day.'), &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),               &
      & prefix, suffix,                                                                  &
      & lrestart=.TRUE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

    CALL mem%Add_var('npp_act_ydaymean', mem%NPP_act_yDayMean,                                      &
      & hgrid, surface,                                                                             &
      & t_cf('NPP_act_yDayMean', 'mol(C) m-2(canopy) s-1', 'Actual NPP Rate of the previous day.'), &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                          &
      & prefix, suffix,                                                                             &
      & lrestart=.TRUE., initval_r=1.0E-13_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3. Initializing to zero would 
                                                  ! result in a zero divide in the first time step of an initialized run.

    CALL mem%Add_var('gpp_ydaymean', mem%GPP_yDayMean,                                        &
      & hgrid, surface,                                                                       &
      & t_cf('GPP_yDayMean', 'mol(C) m-2(canopy) s-1', 'Mean GPP Rate of the previous day.'), &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                    &
      & prefix, suffix,                                                                       &
      & lrestart=.TRUE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

    ! Carbon variables on tile area
    CALL mem%Add_var('npp_act_ydaymean_ta', mem%NPP_act_yDayMean_ta,                                &
      & hgrid, surface,                                                                             &
      & t_cf('NPP_act_yDayMean_ta', 'mol(C) m-2 s-1', 'Actual NPP Rate of the previous day.'), &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                          &
      & prefix, suffix,                                                                             &
      & output_level=BASIC,  loutput=.TRUE.,                                                        &
      & lrestart=.FALSE., initval_r=1.0E-13_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3. Initializing to zero would result
                                                  ! in a zero divide in the first time step of an initialized run.

    CALL mem%Add_var('gpp_ydaymean_ta', mem%GPP_yDayMean_ta,                             &
      & hgrid, surface,                                                                  &
      & t_cf('GPP_yDayMean_ta', 'mol(C) m-2 s-1', 'Mean GPP Rate of the previous day.'), &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),               &
      & prefix, suffix,                                                                  &
      & output_level=BASIC,  loutput=.TRUE.,                                             &
      & lrestart=.FALSE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

  END SUBROUTINE Init_jedi_memory

#endif
END MODULE mo_jedi_memory_class
