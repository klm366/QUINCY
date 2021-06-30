!>
!! @brief Contains the routines for the jedi processes
!!
!! @author
!!  Rainer Schneck, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Rainer Schneck (2015-07-29)
!! Julia Nabel (2018-01-10): l_forestRegrowth flag -> allow for changing max_LAI
!!
MODULE mo_jedi_process
#ifndef __NO_JSBACH__

  USE mo_kind,      ONLY: wp
  USE mo_exception, ONLY: message

  USE mo_jsb_lctlib_class, ONLY: t_lctlib_element

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: calc_Ctradeoff, get_per_tile,jedi_water_stress_factor,jedi_soil_depth_to_layers_1d!,drop_expensive_tissue

  CHARACTER(len=*), PARAMETER :: modname = 'mo_jedi_process'

CONTAINS

  FUNCTION jedi_soil_depth_to_layers_1d(depth, dz) RESULT(depth_l)

    REAL(wp), INTENT(in) :: depth(:)
    REAL(wp), INTENT(in) :: dz(:)
    REAL(wp)             :: depth_l(SIZE(depth,1), SIZE(dz))

    REAL(wp) :: bottom, bottom_m1
    INTEGER  :: nsoil, i

    CHARACTER(len=*), PARAMETER :: routine = modname//':jedi_soil_depth_to_layer_1d'

    nsoil = SIZE(dz)

    depth_l(:,:) = 0._wp
    bottom = 0._wp
    bottom_m1 = 0._wp
    DO i=1,nsoil
      bottom = bottom + dz(i)
      WHERE (depth(:) >= bottom)
        depth_l(:,i) = dz(i)
      ELSEWHERE (depth(:) > bottom_m1)
        depth_l(:,i) = depth(:) - bottom_m1
      END WHERE
      bottom_m1 = bottom_m1 + dz(i)
    END DO

  END FUNCTION jedi_soil_depth_to_layers_1d

  ELEMENTAL FUNCTION jedi_water_stress_factor ( &
    & w_soil, w_soil_max, w_soil_wilt_fract) RESULT(water_stress_factor)

    ! TODO
    ! TBD: Maybe change later so that it uses geographically dependend critical values (like wilting point)

    REAL(wp), INTENT(in) :: w_soil            !< Soil water content
    REAL(wp), INTENT(in) :: w_soil_max        !< Soil water content at field capacity
    
    REAL(wp), INTENT(in) :: w_soil_wilt_fract !< Fraction of max. soil water content at wilting point
    REAL(wp)             :: water_stress_factor

    REAL(wp) :: w_crit, w_wilt                       !< Soil water content at critical/wilting point

    
    w_wilt = w_soil_max * w_soil_wilt_fract

    IF (w_soil_max - w_wilt > 0._wp) THEN
      water_stress_factor = MAX (0._wp, &
        &                        MIN (1._wp, (w_soil - w_wilt) / (w_soil_max - w_wilt) ))
    ELSE
      water_stress_factor = 0._wp
    END IF

 END FUNCTION jedi_water_stress_factor



  ! Calculations to allocate NPP to the plant C pools.
  ELEMENTAL SUBROUTINE calc_Ctradeoff( &
    & GPP_yDayMean    ,           & ! in
    ! JeDi plant trait parameters <<<<<<<<<<
    & GrowRspdMoist,              & ! in
    & GrowRspdTemp ,              & ! in
    & CritiGrowTemp,              & ! in
    & EndGrowTemp,                & !in
    & GermFrac     ,              & ! in
    & Allo2Reproduc,              & ! in
    & Allo2Abground,              & ! in
    & Allo2Blground,              & ! in
    & Allo2Storage ,              & ! in
    & RelAllo2Abgrd,              & ! in
    & RelAllo2Blgrd,              & ! in
    & TurnoverStruc,              & ! in
    & TurnoverLeaf ,              & ! in
    & SenesRspdNPP ,              & ! in
    & NitroStatus  ,              & ! in
    & JeDiSLA,                    & ! in
    & LandcoverClass,             & !in
    ! JeDi plant trait parameters >>>>>>>>>>
    & lspinup_jedi,               & !in 
    & woody_root_coe,             & !in
    & fine_root_coe,              & !in
    & iniyear,                    & ! in
    & currentyear,                & !in
    & time_steps_per_day,         & ! in
    & faPAR_prev_noon,            & !in
    & one_counter,                & !inout
    & c_storage,                  & ! inout
    & c_leaves,                   & ! inout
    & c_fnroot,                   & ! inout
    & c_csroot,                   & ! inout
    & c_seed,                     & ! inout
    & c_stem,                     & ! inout
    & ff_grow_W,                  & ! in
    & ff_grow_T,                  & ! in
    & ff_grow_G,                  & ! in
    & ff_dTdt,                    & !in
    & Flg_grow,                   & ! out
    & Flg_grow_G ,                & ! out
    & ff_npp,                     & ! inout
    & Flg_die,                    & ! inout
    & Flg_senesc,                 & ! out
    & flx_germ,                   & ! out
    & adj_fraction,               &  !out
    & jedi_light_constraint,      &  !out
    & MaintRes_yDayMean,          & ! in
    & GrowthResp,                 & ! out
    & AutoRes,                    & ! out
    & root_depth,                 & ! out root_depth computed from coarse root carbon
    & C_balan,                    & ! out C_balan = GPP - Mainres - Growthres - litter
    & C_delta,                    & ! out C_delta = New total carbon - old total carbon
    & c_total,                    & ! out
    & NPP_act_yDayMean            & ! out
    & )

    USE mo_jedi_constants, ONLY : sec_per_day, days_per_year, sec_per_year
    USE mo_math_constants, ONLY: pi
    
    REAL(wp), INTENT(in) :: GPP_yDayMean           ! Yesterdays mean GPP rate [mol(C)/m^2 s]
    REAL(wp), INTENT(in) :: GrowRspdMoist         ! Jedi parameter -- P01 Growth response time to moisture conditions 
    REAL(wp), INTENT(in) :: GrowRspdTemp          ! Jedi parameter -- P02 Growth response time to temperature conditions 
    REAL(wp), INTENT(in) :: CritiGrowTemp         ! Jedi parameter -- P03 Critical temperature for growth        
    REAL(wp), INTENT(in) :: EndGrowTemp        
    REAL(wp), INTENT(in) :: GermFrac              ! Jedi parameter -- P04 Germination fraction                      
    REAL(wp), INTENT(in) :: Allo2Reproduc         ! Jedi parameter -- P05 Allocation to reproduction                
    REAL(wp), INTENT(in) :: Allo2Abground         ! Jedi parameter -- P06 Allocation to above ground growth           
    REAL(wp), INTENT(in) :: Allo2Blground         ! Jedi parameter -- P07 Allocation to below ground growth           
    REAL(wp), INTENT(in) :: Allo2Storage          ! Jedi parameter -- P08 Allocation to storage                       
    REAL(wp), INTENT(in) :: RelAllo2Abgrd         ! Jedi parameter -- P09 Relative allocation to aboveground growth   
    REAL(wp), INTENT(in) :: RelAllo2Blgrd         ! Jedi parameter -- P10 Relative alloction to belowground growth   
    REAL(wp), INTENT(in) :: TurnoverStruc         ! Jedi parameter -- P11 Turnover time of structural pools           
    REAL(wp), INTENT(in) :: TurnoverLeaf          ! Jedi parameter -- P12 Turnover time of leaf and fine root pools   
    REAL(wp), INTENT(in) :: SenesRspdNPP          ! Jedi parameter -- P13 Senescence response time to productivity conditions 
    REAL(wp), INTENT(in) :: NitroStatus           ! Jedi parameter -- P15 Plant nitrogen status    
    REAL(wp), INTENT(in) :: JeDiSLA      
    INTEGER,  INTENT(in) :: LandcoverClass             
    INTEGER,  INTENT(in) :: iniyear      
    INTEGER,  INTENT(in) :: currentyear  
    INTEGER,  INTENT(in) :: time_steps_per_day
    LOGICAL,  INTENT(IN) :: lspinup_jedi
    REAL(wp), INTENT(IN) :: woody_root_coe
    REAL(wp), INTENT(IN) :: fine_root_coe
    REAL(wp), INTENT(inout) :: one_counter
    !
    ! JeDi carbon allocation 
    REAL(wp), INTENT(inout) :: c_storage            ! Stored carbon pool: on input last value; updated on output
                                                     !    [mol(C)/m^2(canopy)]                     
    REAL(wp), INTENT(inout) :: c_leaves           ! Leaves carbon pool: on input last value; updated on output
                                                     !    [mol(C)/m^2(canopy)]
    REAL(wp), INTENT(inout) :: c_seed             ! Seed carbon pool: on input last value; updated on output
                                                     !    [mol(C)/m^2(canopy)]                                                                                                                                                   
    REAL(wp), INTENT(inout) :: c_fnroot           ! Fine root carbon pool: on input last value; updated on output
                                                     !    [mol(C)/m^2(canopy)]
    REAL(wp), INTENT(inout) :: c_csroot           ! Coarse root carbon pool: on input last value; updated on output
                                                     !    [mol(C)/m^2(canopy)]                                                     
    REAL(wp), INTENT(inout) :: c_stem             ! Stem carbon pool: on input last value; updated on output
                                                     !    [mol(C)/m^2(canopy)]        
    REAL(wp), INTENT(in)    :: ff_grow_T,ff_dTdt
    REAL(wp), INTENT(in)    :: ff_grow_G
    REAL(wp), INTENT(in)    :: ff_grow_W 
    REAL(wp), INTENT(out)   :: Flg_grow
    REAL(wp), INTENT(out)   :: Flg_grow_G
    REAL(wp), INTENT(inout) :: ff_npp
    REAL(wp), INTENT(out)   :: flx_germ  
    REAL(wp), INTENT(out)   :: Flg_senesc
    REAL(wp), INTENT(inout) :: Flg_die    
    REAL(wp), INTENT(out)   :: adj_fraction  
    REAL(wp), INTENT(in)    :: faPAR_prev_noon    
    REAL(wp), INTENT(out)   :: jedi_light_constraint                                                  
 
    REAL(wp), INTENT(out)   :: NPP_act_yDayMean       ! Actual NPP after N-limitation and excess carbon drop [mol(C)/m^2 s]
    REAL(wp), INTENT(in)    :: MaintRes_yDayMean
    REAL(wp), INTENT(out)   :: GrowthResp
    REAL(wp), INTENT(out)   :: AutoRes
    REAL(wp), INTENT(out)   :: C_balan
    REAL(wp), INTENT(out)   :: C_delta
    REAL(wp), INTENT(out)   :: c_total
    REAL(wp), INTENT(out)   :: root_depth
    !
    ! locals
    
    ! JeDi carbon allocation 

    REAL(wp), PARAMETER :: PA0_seed = 1.0_wp / (365._wp * 12.01_wp) ! Virtual carbon for seed pool
    
    REAL(wp) :: A_seed , A_storage
    REAL(wp) :: AG_struct, AG_growth
    REAL(wp) :: BG_struct, BG_growth
    REAL(wp) :: Flg_seed
    REAL(wp) :: Flg_grow_T, Flg_grow_W
    REAL(wp) :: tau_stem, tau_leaves, tau_fnroot, tau_csroot, tau_seed
    REAL(wp) :: ff_germ_Bsoil
    REAL(wp) :: T_degC, fract_veg
    
    REAL(wp), PARAMETER :: c_GRes_leaves    = 1.0_wp / (1.0_wp + 0.33_wp)
    REAL(wp), PARAMETER :: c_GRes_fnroot    = 1.0_wp / (1.0_wp + 0.33_wp)
    REAL(wp), PARAMETER :: c_GRes_wood      = 1.0_wp / (1.0_wp + 0.25_wp)
    REAL(wp), PARAMETER :: c_GRes_csroot    = 1.0_wp / (1.0_wp + 0.49_wp)
    REAL(wp), PARAMETER :: c_GRes_seed      = 1.0_wp / (1.0_wp + 0.66_wp)    
    REAL(wp), PARAMETER :: c_tau_seed       = 3._wp * 365._wp !  seed decomposition time scale (years), value is taken from JeDi 
    REAL(wp), PARAMETER :: c_tau_sen        = 10._wp 
    REAL(wp), PARAMETER :: c_tau_storage    = 10._wp * 365._wp ! decay of the storage pool (1-year), PHS decided arbitrarily
    REAL(wp), parameter :: cTiny            = 1.E-8_wp        
    
    REAL(wp), PARAMETER :: leaf_growth_fract = 0.03_wp   
    REAL(wp), PARAMETER :: fapar_max = 0.96_wp
    REAL(wp) :: allo_max_leaves,zAlloc_leaves_tmp
    
    REAL(wp) :: MaintRes
    REAL(wp) :: NPP
    REAL(WP) :: f_leafroot,f_finetissues,f_stemcsroot                              
    REAL(WP) :: zAlloc, zAlloc_leaves, zAlloc_stem, zAlloc_seed, zAlloc_fnroot, zAlloc_csroot
    REAL(WP) :: zLoss_leaves, zLoss_stem, zLoss_seed, zLoss_fnroot, zLoss_csroot, zLoss_storage 
    REAL(WP) :: c_total_temp
    REAL(wp) :: leaf_shedding 
    
    ! Note: All fields with intent(out) are non-defined on return, even if the calling routine has set them.

    NPP_act_yDayMean = 0.0_wp
      
    ! Initializations

    flx_germ   = 0._wp
    AutoRes    = 0._wp  
    
    GrowthResp = 0._wp
    zAlloc     = 0._wp
    zLoss_seed = 0._wp
    zLoss_storage = 0._wp
    
    AutoRes = MaintRes_yDayMean * sec_per_day
    NPP     = GPP_yDayMean * sec_per_day  - AutoRes
    
    c_total_temp = c_storage + c_seed + c_leaves + c_stem + c_csroot + c_fnroot 
   
    !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=!
    !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-! Growing condition !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-!
    !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=! 
    
    Flg_grow_T = 0._wp
    Flg_grow_W = 0._wp
    Flg_grow_G = 0._wp
    
  !  IF ( ff_grow_T >= CritiGrowTemp ) Flg_grow_T = 1.0_wp
    
    IF ( ff_grow_T >= EndGrowTemp   ) THEN
      Flg_grow_T = 1.0_wp
    ELSE IF ( ff_grow_T > CritiGrowTemp .AND. ff_grow_T <= EndGrowTemp ) THEN
      IF ( ff_dTdt .ge. 0._wp  ) THEN
        Flg_grow_T = 1.0_wp
      ELSE
        Flg_grow_T = 0.0_wp
      END IF
    ELSE 
      Flg_grow_T = 0.0_wp
    END IF
    
    IF ( ff_grow_W >= 0.5_wp ) Flg_grow_W = 1.0_wp
    IF ( ff_grow_G >= 0.5_wp ) Flg_grow_G = 1.0_wp  
    
    Flg_grow   = Flg_grow_T * Flg_grow_W 
    Flg_grow_G = Flg_grow_T * Flg_grow_G
    
    IF ( Flg_grow_G > 0._wp .AND. one_counter < 1._wp )  THEN
      IF ( lspinup_jedi ) THEN        
        c_storage    = 3000
        one_counter  = one_counter +1
      ELSE
        c_seed       = 0.5_wp
        one_counter  = one_counter +1
      END IF
    ENDIF
   
    !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=!
    !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-! Germination condition !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-!
    !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=! 
    
    IF ( Flg_die < 1._wp )  THEN  ! If the species is alive
  
      flx_germ  = Flg_grow_G * GermFrac * c_seed ! if germination is favorable then carbon flows from seed to storage pool 
      c_storage = c_storage + flx_germ
      c_seed    = c_seed    - flx_germ  
      
      zLoss_seed = c_seed / c_tau_seed  ! Why the decay of seed is only computed outside the loop?
      c_seed     = c_seed - zLoss_seed  !  need to think about when the seed will decay maybe only when its after 150 year of simulation or when its not spinup simulation??
      
      !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=!
      !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-! Carbon allocation factors !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-!
      !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=! 
      
      fract_veg = 1.0_wp  ! Dont know what this control, need to check later  ????????
      
      AG_growth = Flg_grow * Allo2Abground
      BG_growth = Flg_grow * Allo2Blground
      
      AG_struct = Flg_grow * RelAllo2Abgrd * fract_veg
      BG_struct = Flg_grow * RelAllo2Blgrd * fract_veg
      
      !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=!
      !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-! Carbon allocation !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-!
      !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=! 
      
      ! Compute JeDi carbon allocation of different C(tissue) pools
    
      IF ( c_storage > 0._wp ) THEN  
	      Flg_seed = 0._wp
    
        IF ( NPP > 0._wp ) Flg_seed = 1.0_wp ! OR NPP_yDayMean ???

        A_seed = Flg_grow * Flg_seed * Allo2Reproduc
        ! PHS:: 14.08.2020 try without turnover for storage ! 
        IF ( .NOT. lspinup_jedi .OR. currentyear >= iniyear + 50 ) THEN 
          IF ( c_storage >  0._wp ) zLoss_storage = c_storage / c_tau_storage 
        END IF      
        c_storage = c_storage + NPP - zLoss_storage
        
        zAlloc = MAX(0._wp,c_storage)  !  To make sure that no negative c_storage distributed to the pools
        
        !PHS: 10.10.2019  I introduce a break to stop rapid growth of leaf
        ! I define the initial amount of carbon that is added to the leaf pool to be 4g C.
        ! Afterwards, the maximal allocation to c_leaves can only be 3% c_leaves of the previous day. (leaf_growth_fract)
        ! Likewise, the other tissues' allocation are constrained by the proportional of the c_leaves allocation.  (adj_fraction)
        IF ( c_leaves <= 4._wp / 12.01_wp )  THEN
          allo_max_leaves = 4._wp / 12.01_wp
        ELSE
          allo_max_leaves = c_leaves * leaf_growth_fract
        END IF
        !PHS: 9.6.2020 FAPAR limitation on leaf growth

        !IF ( faPAR_prev_noon .LT. fapar_max ) THEN
        !PHS: 15.7.2020 try to  use different light constraint for grass 
        IF ( LandcoverClass == 3) THEN      
          IF ( faPAR_prev_noon < fapar_max ) THEN
            jedi_light_constraint = COS( 0.5_wp * pi * ( (faPAR_prev_noon / fapar_max )**30._wp) )
          ELSE
            jedi_light_constraint = 0._wp
          END IF
        ELSE IF (LandcoverClass == 4) THEN
          IF ( faPAR_prev_noon < 0.83_wp ) THEN 
            jedi_light_constraint = COS( 0.5_wp * pi * ( (faPAR_prev_noon / 0.83_wp )**30._wp) )
          ELSE
            jedi_light_constraint = 0._wp
          END IF
        END IF   
        !ELSE
            !jedi_light_constraint = 0._wp
        !END IF
        
        zAlloc_leaves_tmp = zAlloc * AG_growth * ( 1._wp - AG_struct ) * c_GRes_leaves  
        zAlloc_leaves = MIN( zAlloc_leaves_tmp, allo_max_leaves )  
        
        IF (zAlloc_leaves /= 0._wp  ) THEN   ! RS: this is dangerous! TODO
          adj_fraction = zAlloc_leaves / zAlloc_leaves_tmp
        ELSE
          adj_fraction  = 1._wp
        END IF
        
        zAlloc_leaves = zAlloc_leaves * jedi_light_constraint

        zAlloc_seed   = zAlloc * A_seed                            * c_GRes_seed   * adj_fraction  * jedi_light_constraint
        !zAlloc_leaves = zAlloc * AG_growth * ( 1._wp - AG_struct ) * c_GRes_leaves
        zAlloc_stem   = zAlloc * AG_growth * AG_struct             * c_GRes_wood   * adj_fraction  * jedi_light_constraint
        zAlloc_csroot = zAlloc * BG_growth * BG_struct             * c_GRes_csroot * adj_fraction  * jedi_light_constraint
        zAlloc_fnroot = zAlloc * BG_growth * ( 1._wp - BG_struct ) * c_GRes_fnroot * adj_fraction  * jedi_light_constraint
        
        c_storage = c_storage - ( A_seed + AG_growth + BG_growth ) * zAlloc        * adj_fraction  * jedi_light_constraint
        
        !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=!
        !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-! Growth respiration !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-!
        !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=! 

        GrowthResp =   ( zAlloc_leaves / c_GRes_leaves ) * ( 1.0_wp - c_GRes_leaves ) &
          &          + ( zAlloc_fnroot / c_GRes_fnroot ) * ( 1.0_wp - c_GRes_fnroot ) &
          &          + ( zAlloc_stem   / c_GRes_wood   ) * ( 1.0_wp - c_GRes_wood   ) &
          &          + ( zAlloc_csroot / c_GRes_csroot ) * ( 1.0_wp - c_GRes_csroot ) &
          &          + ( zAlloc_seed   / c_GRes_seed   ) * ( 1.0_wp - c_GRes_seed   ) 
        
        AutoRes = AutoRes + GrowthResp  
        
        !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=!
        !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-! Senescence condition !=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-!
        !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=! 

        Flg_senesc = 0._wp
        IF ( .NOT. lspinup_jedi .OR. currentyear >= iniyear + 50 ) THEN 
          IF ( ff_npp < 0._wp .AND. NPP < 0._wp ) Flg_senesc = 1.0_wp - Flg_grow   ! See remarks below.
         END IF
        ! No idea why the flag for senescence is using:  1 - flag_grow  
        ! I took the calculation from JeDi.
        ! I guess senescence happens, when: (1) Non-growing condition (2) NPP<0, which may suggest that plants spend many energy respiring
        !                                   (3) C_storage>0           (4) ff_npp<0 
        
        tau_stem   = TurnoverStruc
        tau_csroot = TurnoverStruc
        IF (LandcoverClass == 3) THEN     
          tau_leaves = 1._wp / ( (1._wp / TurnoverLeaf ) + ( 1._wp / TurnoverLeaf ) * Flg_senesc )  ! *             RelSenesAbgd )      !PHS 11.12.2019 play around with the senescence process
          tau_fnroot = 1._wp / ( (1._wp / TurnoverLeaf ) + ( 1._wp / TurnoverLeaf ) * Flg_senesc  ) ! * ( 1._wp -   RelSenesAbgd) )!PHS 11.12.2019 play around with the senescence process
        ELSE IF (LandcoverClass == 4) THEN  
          tau_leaves = 1._wp / ( (1._wp / TurnoverLeaf ) + ( 1._wp / (TurnoverLeaf/2._wp) ) * Flg_senesc )  ! *             RelSenesAbgd )      !PHS 11.12.2019 play around with the senescence process
          tau_fnroot = 1._wp / ( (1._wp / TurnoverLeaf ) + ( 1._wp / (TurnoverLeaf/2._wp) ) * Flg_senesc  ) ! * ( 1._wp -   RelSenesAbgd) )!PHS 11.12.2019 play around with the senescence process
        END IF

        !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=!
        !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-! Compute turnover Carbon !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-!
        !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=! 
        
        zLoss_leaves = c_leaves / tau_leaves
        zLoss_stem   = c_stem   / tau_stem
        zLoss_csroot = c_csroot / tau_csroot
        zLoss_fnroot = c_fnroot / tau_fnroot

        !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=!
        !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-! Update Carbon pools !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-!
        !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=! 

        c_seed   = c_seed   + zAlloc_seed        
        c_leaves = c_leaves + zAlloc_leaves - zLoss_leaves        
        c_stem   = c_stem   + zAlloc_stem   - zLoss_stem        
        c_csroot = c_csroot + zAlloc_csroot - zLoss_csroot       
        c_fnroot = c_fnroot + zAlloc_fnroot - zLoss_fnroot

        !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=!
        !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-! Death of plants     !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-!
        !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=! 
         
        NPP_act_yDayMean = ( GPP_yDayMean * sec_per_day  - AutoRes ) / (  sec_per_day )

        IF ( LandcoverClass == 3) THEN     
         ! root_depth = sqrt( c_csroot * 2250._wp * 12.01_wp ) * 0.001_wp
         root_depth = sqrt( c_csroot * woody_root_coe * 12.01_wp ) * 0.001_wp
        ELSE IF (LandcoverClass == 4) THEN
         ! root_depth = sqrt( c_fnroot * 550._wp * 12.01_wp ) * 0.001_wp
          root_depth = sqrt( c_fnroot * fine_root_coe * 12.01_wp ) * 0.001_wp
        END IF

        GrowthResp = GrowthResp / sec_per_day 
        AutoRes    = AutoRes    / sec_per_day 

        IF (  c_storage <= cTiny )  THEN 
          zLoss_leaves = zLoss_leaves + c_leaves
          zLoss_stem   = zLoss_stem   + c_stem    + c_storage ! I dont know where the dead storage carbon shall go so I first put it into stem letter pool
          zLoss_csroot = zLoss_csroot + c_csroot  + c_seed ! I dont know where should dead seed carbon go so i put into coarse root's litter pool first
          zLoss_fnroot = zLoss_fnroot + c_fnroot
            
          c_storage = 0._wp   ! Should the storage carbon go to some litter pools?
          c_seed    = 0._wp
          c_leaves  = 0._wp
          c_stem    = 0._wp
          c_csroot  = 0._wp
          c_fnroot  = 0._wp
          Flg_die   = 1._wp
        END IF

      ENDIF

    ELSE 
      GrowthResp                   = 0._wp
      AutoRes                      = 0._wp
      NPP_act_yDayMean             = 0._wp
      root_depth                   = 0._wp   
    ENDIF

    ! PHS: test biomass aggregtion without storage   
    c_total = c_leaves + c_stem + c_csroot + c_fnroot 
      !c_total = c_storage + c_leaves + c_stem + c_csroot + c_fnroot 

    C_delta =  ( c_storage + c_seed + c_leaves + c_stem + c_csroot + c_fnroot ) - c_total_temp  ! this is the difference of assimilated carbon added into c-pools in this timestep
   

    C_balan =   GPP_yDayMean   *  sec_per_day  &
      &       - MaintRes_yDayMean  *  sec_per_day  &
      &       - GrowthResp         *  sec_per_day  &
      &       - ( zLoss_leaves + zLoss_fnroot + zLoss_stem + zLoss_csroot + zLoss_seed + zLoss_storage) 
   
  END SUBROUTINE calc_Ctradeoff

  ! --- get_per_tile() ---------------------------------------------------------------------------------------------------
  !
  ! calculate per tile density from density per canopy m2
  !
  ELEMENTAL SUBROUTINE get_per_tile(  &
    & Cpool_per_tile,                 & ! out
    & Cpool_per_canopy,               & ! in
    & veg_fract_correction,           & ! in
    & fract_fpc_max)                    ! in

    REAL(wp), INTENT(IN) ::  fract_fpc_max, veg_fract_correction
    REAL(wp), INTENT(IN) ::  Cpool_per_canopy
    REAL(wp), INTENT(OUT)::  Cpool_per_tile

    Cpool_per_tile = Cpool_per_canopy * veg_fract_correction * fract_fpc_max

  END SUBROUTINE get_per_tile

#endif
END MODULE mo_jedi_process
