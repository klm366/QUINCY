!> Functions that define and create QUINCY model output as ascii text
!>
!> QUINCY standalone
!>
!> ---------------------------------------
!> Copyright (C) 2013-2026, MPI-M, MPI-BGC
!>
!> Contact: https://doi.org/10.17871/quincy-model-2019
!> Authors: AUTHORS.md
!> See LICENSES/ for license information
!> SPDX-License-Identifier: BSD-3-Clause
!> ---------------------------------------
!> 
!> three main sections
!>   - I   Pools 
!>   - II  Fluxes 
!>   - III Soil forcing
!> 
MODULE mo_qs_output_txt

  USE mo_kind,                  ONLY: wp
  USE mo_jsb_control,           ONLY: debug_on
  USE mo_exception,             ONLY: message_text, message, finish
  USE mo_jsb_lctlib_class,      ONLY: t_lctlib_element
  USE mo_qs_forcing,            ONLY: forcing_options
  USE mo_quincy_model_config,   ONLY: QLAND, QPLANT, QSOIL, QCANOPY, Q_TEST_CANOPY, Q_TEST_RADIATION
  USE mo_qs_output_txt_util

  USE mo_lnd_bgcm_idx
  USE mo_lnd_bgcm_store,          ONLY: t_lnd_bgcm_store
  USE mo_lnd_bgcm_store_class,    ONLY: VEG_BGCM_POOL_ID, VEG_BGCM_SEED_BED_POOL_ID, VEG_BGCM_GROWTH_ID, &
    &                                   VEG_BGCM_SEED_BED_GROWTH_ID, VEG_BGCM_LITTERFALL_ID, VEG_BGCM_FRAC_ALLOC_ID, &
    &                                   VEG_BGCM_SEED_BED_LITTERFALL_ID, VEG_BGCM_EXUDATION_ID, VEG_BGCM_ESTABLISHMENT_ID, &
    &                                   VEG_BGCM_RESERVE_USE_ID
  USE mo_lnd_bgcm_store_class,    ONLY: SB_BGCM_POOL_ID, SB_BGCM_DELTA_POOLS_ID, SB_BGCM_FORMATION_ID, SB_BGCM_LOSS_ID, &
    &                                   SB_BGCM_TRANSPORT_ID, SB_BGCM_MYCO_EXPORT_ID

  dsl4jsb_Use_config(SB_)

  IMPLICIT NONE

  PUBLIC  ! all is public

  INTEGER, PARAMETER, PRIVATE   :: ic = 1                             !< the point in ics:ice of which data will be written to file; per default = 1
  INTEGER, PARAMETER            :: lun_sim_runtime_estimate   = 90    !< used in mo_qs_time_util
  INTEGER, PRIVATE, PARAMETER   :: lun_ascii_soil_forc_start  = 200   !< start LUN for writing out forcing for the soil model-mode
  INTEGER, PRIVATE, PARAMETER   :: lun_ascii_soil_forc_max    = 204   !< max   LUN for writing out forcing for the soil model-mode

  CHARACTER(len=*), PARAMETER, PRIVATE :: modname = 'mo_qs_output_txt'

CONTAINS

  !-----------------------------------------------------------------------------------------------------!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! section I Pools                                                                                   !!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !-----------------------------------------------------------------------------------------------------!


  ! ======================================================================================================= !
  !>Pools: open files and write column headers
  !>
  SUBROUTINE open_pool_output_files(output_interval_pool, spinup_id)

    CHARACTER(len=*), INTENT(IN)      :: output_interval_pool   !< output interval for "pool" variables
    CHARACTER(len=*), INTENT(IN)      :: spinup_id              !< add ID to filename if file contains data from spinup simulation
    ! ----------------------------------------------------------------------------------------------------- !
    CHARACTER(100)              :: fname    !< filename helper
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':open_pool_output_files'
    ! ----------------------------------------------------------------------------------------------------- !

    !>1.0 define open() statements for active output files
    !>
    IF (is_active_txt_vegpool_c) THEN
      fname = name_vegpool_c//'_'//TRIM(spinup_id)//TRIM(output_interval_pool)//'.txt'
      OPEN (UNIT=lun_vegpool_c, FILE=fname, ACTION="write", STATUS="unknown")
    END IF
    IF (is_active_txt_vegpool_dc13dc14) THEN
      fname = name_vegpool_dc13dc14//'_'//TRIM(spinup_id)//TRIM(output_interval_pool)//'.txt'
      OPEN (UNIT=lun_vegpool_dc13dc14, FILE=fname, ACTION="write", STATUS="unknown")
    END IF
    IF (is_active_txt_vegpool_n) THEN
      fname = name_vegpool_n//'_'//TRIM(spinup_id)//TRIM(output_interval_pool)//'.txt'
      OPEN (UNIT=lun_vegpool_n, FILE=fname, ACTION="write", STATUS="unknown")
    END IF
    IF (is_active_txt_vegpool_p) THEN
      fname = name_vegpool_p//'_'//TRIM(spinup_id)//TRIM(output_interval_pool)//'.txt'
      OPEN (UNIT=lun_vegpool_p, FILE=fname, ACTION="write", STATUS="unknown")
    END IF
    IF (is_active_txt_vegpool_dn15) THEN
      fname = name_vegpool_dn15//'_'//TRIM(spinup_id)//TRIM(output_interval_pool)//'.txt'
      OPEN (UNIT=lun_vegpool_dn15, FILE=fname, ACTION="write", STATUS="unknown")
    END IF
    IF (is_active_txt_soil_physics) THEN
      fname = name_soil_physics//'_'//TRIM(spinup_id)//TRIM(output_interval_pool)//'.txt'
      OPEN (UNIT=lun_soil_physics, FILE=fname, ACTION="write", STATUS="unknown")
    END IF
    IF (is_active_txt_soil_physics_sl) THEN
      fname = name_soil_physics_sl//'_'//TRIM(spinup_id)//TRIM(output_interval_pool)//'.txt'
      OPEN (UNIT=lun_soil_physics_sl, FILE=fname, ACTION="write", STATUS="unknown")
    END IF
    IF (is_active_txt_soilpool_inorg) THEN
      fname = name_soilpool_inorg//'_'//TRIM(spinup_id)//TRIM(output_interval_pool)//'.txt'
      OPEN (UNIT=lun_soilpool_inorg, FILE=fname, ACTION="write", STATUS="unknown")
    END IF
    IF (is_active_txt_ecosyspool_cnp) THEN
      fname = name_ecosyspool_cnp//'_'//TRIM(spinup_id)//TRIM(output_interval_pool)//'.txt'
      OPEN (UNIT=lun_ecosyspool_cnp, FILE=fname, ACTION="write", STATUS="unknown")
    END IF
    IF (is_active_txt_soilpool_dc13dc14dn15) THEN
      fname = name_soilpool_dc13dc14dn15//'_'//TRIM(spinup_id)//TRIM(output_interval_pool)//'.txt'
      OPEN (UNIT=lun_soilpool_dc13dc14dn15, FILE=fname, ACTION="write", STATUS="unknown")
    END IF
    IF (is_active_txt_veg_diagnostics) THEN
      fname = name_veg_diagnostics//'_'//TRIM(spinup_id)//TRIM(output_interval_pool)//'.txt'
      OPEN (UNIT=lun_veg_diagnostics, FILE=fname, ACTION="write", STATUS="unknown")
    END IF
    IF (is_active_txt_veg_canopy_diagnostics_cl) THEN
      fname = name_veg_canopy_diagnostics_cl//'_'//TRIM(spinup_id)//TRIM(output_interval_pool)//'.txt'
      OPEN (UNIT=lun_veg_canopy_diagnostics_cl, FILE=fname, ACTION="write", STATUS="unknown")
    END IF
    IF (is_active_txt_soilpool_org) THEN
      fname = name_soilpool_org//'_'//TRIM(spinup_id)//TRIM(output_interval_pool)//'.txt'
      OPEN (UNIT=lun_soilpool_org, FILE=fname, ACTION="write", STATUS="unknown")
    END IF
    IF (is_active_txt_soilpool_c_org_sl) THEN
      fname = name_soilpool_c_org_sl//'_'//TRIM(spinup_id)//TRIM(output_interval_pool)//'.txt'
      OPEN (UNIT=lun_soilpool_c_org_sl, FILE=fname, ACTION="write", STATUS="unknown")
    END IF
    IF (is_active_txt_soilpool_n_org_sl) THEN
      fname = name_soilpool_n_org_sl//'_'//TRIM(spinup_id)//TRIM(output_interval_pool)//'.txt'
      OPEN (UNIT=lun_soilpool_n_org_sl, FILE=fname, ACTION="write", STATUS="unknown")
    END IF
    IF (is_active_txt_soilpool_p_org_sl) THEN
      fname = name_soilpool_p_org_sl//'_'//TRIM(spinup_id)//TRIM(output_interval_pool)//'.txt'
      OPEN (UNIT=lun_soilpool_p_org_sl, FILE=fname, ACTION="write", STATUS="unknown")
    END IF
    IF (is_active_txt_soilpool_ndn15_inorg_sl) THEN
      fname = name_soilpool_ndn15_inorg_sl//'_'//TRIM(spinup_id)//TRIM(output_interval_pool)//'.txt'
      OPEN (UNIT=lun_soilpool_ndn15_inorg_sl, FILE=fname, ACTION="write", STATUS="unknown")
    END IF
    IF (is_active_txt_soilpool_p_inorg_sl) THEN
      fname = name_soilpool_p_inorg_sl//'_'//TRIM(spinup_id)//TRIM(output_interval_pool)//'.txt'
      OPEN (UNIT=lun_soilpool_p_inorg_sl, FILE=fname, ACTION="write", STATUS="unknown")
    END IF
    IF (is_active_txt_soilpool_dc13_org_sl) THEN
      fname = name_soilpool_dc13_org_sl//'_'//TRIM(spinup_id)//TRIM(output_interval_pool)//'.txt'
      OPEN (UNIT=lun_soilpool_dc13_org_sl, FILE=fname, ACTION="write", STATUS="unknown")
    END IF
    IF (is_active_txt_soilpool_dc14_org_sl) THEN
      fname = name_soilpool_dc14_org_sl//'_'//TRIM(spinup_id)//TRIM(output_interval_pool)//'.txt'
      OPEN (UNIT=lun_soilpool_dc14_org_sl, FILE=fname, ACTION="write", STATUS="unknown")
    END IF
    IF (is_active_txt_soilpool_dn15_org_sl) THEN
      fname = name_soilpool_dn15_org_sl//'_'//TRIM(spinup_id)//TRIM(output_interval_pool)//'.txt'
      OPEN (UNIT=lun_soilpool_dn15_org_sl, FILE=fname, ACTION="write", STATUS="unknown")
    END IF
    IF (is_active_txt_soil_diagnostics_microbial_sl) THEN
      fname = name_soil_diagnostics_microbial_sl//'_'//TRIM(spinup_id)//TRIM(output_interval_pool)//'.txt'
      OPEN (UNIT=lun_soil_diagnostics_microbial_sl, FILE=fname, ACTION="write", STATUS="unknown")
    END IF
    IF (is_active_txt_soil_diagnostics_sorption_sl) THEN
      fname = name_soil_diagnostics_sorption_sl//'_'//TRIM(spinup_id)//TRIM(output_interval_pool)//'.txt'
      OPEN (UNIT=lun_soil_diagnostics_sorption_sl, FILE=fname, ACTION="write", STATUS="unknown")
    END IF

    !>2.0 write column headers into output files
    !>
    CALL init_output_writing_pool()
  END SUBROUTINE open_pool_output_files


  ! ======================================================================================================= !
  !>Pools: write column headers of output files
  !>
  SUBROUTINE init_output_writing_pool()

    USE mo_jsb_grid_class,        ONLY: t_jsb_vgrid
    USE mo_jsb_grid,              ONLY: Get_vgrid
    ! ----------------------------------------------------------------------------------------------------- !
    TYPE(t_jsb_vgrid),POINTER :: vgrid_soil_w                ! soil layer vertical grid
    TYPE(t_jsb_vgrid),POINTER :: vgrid_canopy_q_assimi       ! canopy layer vertical grid
    INTEGER                   :: nsoil_w                     ! number of soil layers as used/defined by the SB_ process
    INTEGER                   :: ncanopy                     ! number of canopy layers
    INTEGER                   :: isoil, icanopy
    CHARACTER(2)              :: layerstr
    CHARACTER(450)            :: p1str, p2str, p3str, p4str, p5str, p6str, p7str, p8str, p9str, &
      &                          p10str, p11str, p12str, p13str, p14str
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':init_output_writing_pool'
    ! ----------------------------------------------------------------------------------------------------- !
    vgrid_soil_w           => Get_vgrid('soil_depth_water')
    nsoil_w                =  vgrid_soil_w%n_levels
    vgrid_canopy_q_assimi  => Get_vgrid('q_canopy_layer')
    ncanopy                =  vgrid_canopy_q_assimi%n_levels
    ! ----------------------------------------------------------------------------------------------------- !

    !> 0.9 test max number of canopy layers (15) assumed for writing name_veg_canopy_diagnostics_cl
    !>
    IF (ncanopy > 15) THEN
      WRITE(message_text,'(a)') 'number of canopy layers (q_canopy_layer vgrid) is larger than 15'
      CALL finish(TRIM(routine), message_text)
    END IF

    !> 1.0 write headers of active output files
    !>
    IF (is_active_txt_vegpool_c) THEN
      ! name_vegpool_c [mol C / m2]
      WRITE(lun_vegpool_c,'(2a)') &
        &   'Time Leaf_C FineRoot_C SapWood_C HeartWood_C CoarseRoot_C ', &
        &   'Fruit_C Labile_C Reserve_C SeedBank_C Total_C Total_C13 Total_C14'
    END IF
    IF (is_active_txt_vegpool_dc13dc14) THEN
      ! name_vegpool_dc13dc14 [ ... ]
      WRITE(lun_vegpool_dc13dc14,'(4a)') &
        &   'Time Leaf_dC13 FineRoot_dC13 SapWood_dC13 HeartWood_dC13 CoarseRoot_dC13 ', &
        &   'Fruit_dC13 Labile_dC13 Reserve_dC13 SeedBank_dC13 Leaf_DC14 FineRoot_DC14 ', &
        &   'SapWood_DC14 HeartWood_DC14 CoarseRoot_DC14 Fruit_DC14 Labile_DC14 ', &
        &   'Reserve_DC14 SeedBank_DC14'
    END IF
    IF (is_active_txt_vegpool_n) THEN
      ! name_vegpool_n [mol N / m2]
      WRITE(lun_vegpool_n,'(2a)') &
        &   'Time Leaf_N FineRoot_N SapWood_N HeartWood_N CoarseRoot_N ', &
        &   'Fruit_N Labile_N Reserve_N SeedBank_N Total_N'
    END IF
    IF (is_active_txt_vegpool_p) THEN
      ! name_vegpool_p [mol P / m2]
      WRITE(lun_vegpool_p,'(2a)') &
        &   'Time Leaf_P FineRoot_P SapWood_P HeartWood_P CoarseRoot_P ', &
        &   'Fruit_P Labile_P Reserve_P SeedBank_P Total_P'
    END IF
    IF (is_active_txt_vegpool_dn15) THEN
      ! name_vegpool_dn15 [ ... ]
      WRITE(lun_vegpool_dn15,'(2a)') &
        &   'Time Leaf_dN15 FineRoot_dN15 SapWood_dN15 HeartWood_dN15 CoarseRoot_dN15 ', &
        &   'Fruit_dN15 Labile_dN15 Reserve_dN15 SeedBank_dN15 Total_dN15'
    END IF
    IF (is_active_txt_soil_physics) THEN
      ! name_soil_physics [m, m, m]
      WRITE(lun_soil_physics,'(a)') 'Time SnowHeight RootZoneSoilMoisture RootZoneSoilWaterPotential pH ph_decomp_rate_impact'
    END IF
    IF (is_active_txt_soil_physics_sl) THEN
      ! name_soil_physics_sl [m, K] (output per soil_layer)
      p1str=''
      p2str=''
      DO isoil = 1,nsoil_w
        IF (isoil < 10) THEN 
          WRITE(layerstr,'(I1)') isoil
        ELSE
          WRITE(layerstr,'(I2)') isoil
        END IF
        WRITE(p1str,'(a)') TRIM(p1str)//' WaterContentSoil_L'//TRIM(layerstr)
        WRITE(p2str,'(a)') TRIM(p2str)//' TemperatureSoil_L'//TRIM(layerstr)
      END DO
      WRITE(lun_soil_physics_sl,'(3a)') 'Time ',TRIM(p1str), TRIM(p2str)
    END IF
    IF (is_active_txt_soilpool_inorg) THEN
      ! name_soilpool_inorg [mol / m2]
      WRITE(lun_soilpool_inorg,'(a)') 'Time MineralAssocNH4 MineralNH4 MineralNO3 MineralP NOY N2O N2'
    END IF
    IF (is_active_txt_ecosyspool_cnp) THEN
      ! name_ecosyspool_cnp [mol / m2]
      WRITE(lun_ecosyspool_cnp,'(a)') 'Time TotalSoilOrgC TotalSoilOrgN TotalSoilInOrgN EcoC EcoN EcoP EcoN15'
    END IF
    IF (is_active_txt_soilpool_dc13dc14dn15) THEN
      ! name_soilpool_dc13dc14dn15 [ ... ]
      WRITE(lun_soilpool_dc13dc14dn15,'(4a)') &
        &   'Time MetLitterdC13 StrLitterdC13 WodLitterdC13 FastdC13 SlowdC13 MycdC13 MetLitterDC14 StrLitterDC14',& 
        &   ' WodLitterDC14 FastDC14 SlowDC14 MycDC14',&
        &   ' MetLitterdN15 StrLitterdN15 WodLitterN15', &
        &   ' FastdN15 SlowdN15 MineralAssocNH4_dN15 MineralNH4_dN15 MineralNO3_dN15'
    END IF
    IF (is_active_txt_veg_diagnostics) THEN
      ! name_veg_diagnostics [ ... ]
      WRITE(lun_veg_diagnostics,'(5a)') &
        &   'Time Height Diameter LAI ind_density LeafM2RootM', &
        &   ' growing_season Ndemand Pdemand sinklim', &
        &   ' growthReqN_uptMavg growthReqP_uptMavg', &
        &   ' BetaSoilGS BetaSoilPS ph_fine_root_factor'
    END IF
    IF (is_active_txt_veg_canopy_diagnostics_cl) THEN
      ! name_veg_canopy_diagnostics_cl [ ... ]
      p1str=''
      p2str=''
      p3str=''
      p4str=''
      p5str=''
      p6str=''
      p7str=''
      p8str=''
      p9str=''
      DO icanopy = 1, ncanopy
        IF (icanopy < 10) THEN 
          WRITE(layerstr,'(I1)') icanopy
        ELSE
          WRITE(layerstr,'(I2)') icanopy
        END IF
        WRITE(p1str,'(a)') TRIM(p1str)//' f_et_l'//TRIM(layerstr)
        WRITE(p2str,'(a)') TRIM(p2str)//' f_rub_l'//TRIM(layerstr)
        WRITE(p3str,'(a)') TRIM(p3str)//' f_chl_l'//TRIM(layerstr)
        WRITE(p4str,'(a)') TRIM(p4str)//' f_pepc_l'//TRIM(layerstr)
        WRITE(p5str,'(a)') TRIM(p5str)//' f_oth_l'//TRIM(layerstr)
        WRITE(p6str,'(a)') TRIM(p6str)//' leaf_n_l'//TRIM(layerstr)
        WRITE(p7str,'(a)') TRIM(p7str)//' jmax_l'//TRIM(layerstr)
        WRITE(p8str,'(a)') TRIM(p8str)//' vcmax_l'//TRIM(layerstr)
        WRITE(p9str,'(a)') TRIM(p9str)//' chlfl_yield_l'//TRIM(layerstr)
      END DO
      WRITE(lun_veg_canopy_diagnostics_cl,'(10a)') &
        &   'Time ', TRIM(p1str), TRIM(p2str), TRIM(p3str), TRIM(p4str), TRIM(p5str), &
        &   TRIM(p6str), TRIM(p7str), TRIM(p8str), TRIM(p9str)
    END IF
    IF (is_active_txt_soilpool_org) THEN
      ! name_soilpool_org [mol / m2]
      ! Note: with the ascii output the harvest litter flux and the product pools as well as product decay is not written out.
      !       If one wants to close the budget in the output when using sylviculture, one would need to add the litter flux,
      !       and when using product pools, in addition, also the product pools as well as product decay would need to be added.
      WRITE(lun_soilpool_org,'(7a)') &
        &   'Time ', &
        &   'poolC_dom poolC_domAs poolC_solLitter poolC_polLitter poolC_wooLitter ', &
        &   'poolC_myco poolC_mic poolC_res poolC_resAs ', &
        &   'poolN_dom poolN_domAs poolN_solLitter poolN_polLitter poolN_wooLitter ', &
        &   'poolN_myco poolN_mic poolN_res poolN_resAs ', &
        &   'poolP_dom poolP_domAs poolP_solLitter poolP_polLitter poolP_wooLitter ', &
        &   'poolP_myco poolP_mic poolP_res poolP_resAs'
    END IF
    IF (is_active_txt_soilpool_c_org_sl) THEN
      ! name_soilpool_c_org_sl [mol C / m3]
      p1str=''
      p2str=''
      p3str=''
      p4str=''
      p5str=''
      p6str=''
      p7str=''
      p8str=''
      p9str=''
      DO isoil = 1,nsoil_w
        IF (isoil < 10) THEN 
          WRITE(layerstr,'(I1)') isoil
        ELSE
          WRITE(layerstr,'(I2)') isoil
        END IF
        WRITE(p1str,'(a)') TRIM(p1str)//' SoluLit_L'//TRIM(layerstr)
        WRITE(p2str,'(a)') TRIM(p2str)//' PolyLit_L'//TRIM(layerstr)
        WRITE(p3str,'(a)') TRIM(p3str)//' WodLit_L'//TRIM(layerstr)
        WRITE(p4str,'(a)') TRIM(p4str)//' DOM_L'//TRIM(layerstr)
        WRITE(p5str,'(a)') TRIM(p5str)//' aDOM_L'//TRIM(layerstr)
        WRITE(p6str,'(a)') TRIM(p6str)//' Myc_L'//TRIM(layerstr)
        WRITE(p7str,'(a)') TRIM(p7str)//' Mic_L'//TRIM(layerstr)
        WRITE(p8str,'(a)') TRIM(p8str)//' Res_L'//TRIM(layerstr)
        WRITE(p9str,'(a)') TRIM(p9str)//' aRes_L'//TRIM(layerstr)
      END DO
      WRITE(lun_soilpool_c_org_sl,'(10a)') &
        &   'Time ',TRIM(p1str),TRIM(p2str),TRIM(p3str),TRIM(p4str),TRIM(p5str),TRIM(p6str),TRIM(p7str), &
        &   TRIM(p8str),TRIM(p9str)
    END IF
    IF (is_active_txt_soilpool_n_org_sl) THEN
      ! name_soilpool_n_org_sl [mol N / m3]
      WRITE(lun_soilpool_n_org_sl,'(10a)') &
        &   'Time ',TRIM(p1str),TRIM(p2str),TRIM(p3str),TRIM(p4str),TRIM(p5str),TRIM(p6str),TRIM(p7str), &
        &   TRIM(p8str),TRIM(p9str)
    END IF
    IF (is_active_txt_soilpool_p_org_sl) THEN
      ! name_soilpool_p_org_sl [mol P / m3]
      WRITE(lun_soilpool_p_org_sl,'(10a)') &
        &   'Time ',TRIM(p1str),TRIM(p2str),TRIM(p3str),TRIM(p4str),TRIM(p5str),TRIM(p6str),TRIM(p7str), &
        &   TRIM(p8str),TRIM(p9str)
    END IF
    IF (is_active_txt_soilpool_ndn15_inorg_sl) THEN
      ! name_soilpool_ndn15_inorg_sl [mol N / m3]
      p1str=''
      p2str=''
      p3str=''
      p4str=''
      p5str=''
      p6str=''
      p7str=''
      p8str=''
      p9str=''
      DO isoil = 1,nsoil_w
        IF (isoil < 10) THEN 
          WRITE(layerstr,'(I1)') isoil
        ELSE
          WRITE(layerstr,'(I2)') isoil
        END IF
        WRITE(p1str,'(a)') TRIM(p1str)//' NH4assoc_L'//TRIM(layerstr)
        WRITE(p2str,'(a)') TRIM(p2str)//' NH4solute_L'//TRIM(layerstr)
        WRITE(p3str,'(a)') TRIM(p3str)//' NO3solute_L'//TRIM(layerstr)
        WRITE(p4str,'(a)') TRIM(p4str)//' NOY_L'//TRIM(layerstr)
        WRITE(p5str,'(a)') TRIM(p5str)//' N2O_L'//TRIM(layerstr)
        WRITE(p6str,'(a)') TRIM(p6str)//' N2_L'//TRIM(layerstr)
        WRITE(p7str,'(a)') TRIM(p7str)//' dN15_NH4assoc_L'//TRIM(layerstr)
        WRITE(p8str,'(a)') TRIM(p8str)//' dN15_NH4solute_L'//TRIM(layerstr)
        WRITE(p9str,'(a)') TRIM(p9str)//' dN15_NO3solute_L'//TRIM(layerstr)
      END DO
      WRITE(lun_soilpool_ndn15_inorg_sl,'(10a)') &
        &   'Time ',TRIM(p1str),TRIM(p2str),TRIM(p3str),TRIM(p4str),TRIM(p5str),TRIM(p6str), &
        &   TRIM(p7str),TRIM(p8str),TRIM(p9str)
    END IF
    IF (is_active_txt_soilpool_p_inorg_sl) THEN
      ! name_soilpool_p_inorg_sl [mol P / m3]
      p1str=''
      p2str=''
      p3str=''
      p4str=''
      p5str=''
      DO isoil = 1,nsoil_w
        IF (isoil < 10) THEN 
          WRITE(layerstr,'(I1)') isoil
        ELSE
          WRITE(layerstr,'(I2)') isoil
        END IF
        WRITE(p1str,'(a)') TRIM(p1str)//' PO4assoc_fast_L'//TRIM(layerstr)
        WRITE(p2str,'(a)') TRIM(p2str)//' PO4assoc_slow_L'//TRIM(layerstr)
        WRITE(p3str,'(a)') TRIM(p3str)//' PO4assoc_occl_L'//TRIM(layerstr)
        WRITE(p4str,'(a)') TRIM(p4str)//' PO4solute_L'//TRIM(layerstr)
        WRITE(p5str,'(a)') TRIM(p5str)//' PO4primary_L'//TRIM(layerstr)
      END DO
      WRITE(lun_soilpool_p_inorg_sl,'(6a)') 'Time ',TRIM(p1str),TRIM(p2str),TRIM(p3str),TRIM(p4str),TRIM(p5str)
    END IF
    IF (is_active_txt_soil_diagnostics_microbial_sl) THEN
      ! name_soil_diagnostics_microbial_sl [mol C / m3, kg/m3, na, mol P / m3, mol P / m3, mol C / m3]
      p1str=''
      p2str=''
      p3str=''
      p4str=''
      p5str=''
      p6str=''
      p7str=''
      p8str=''
      p9str=''
      p10str=''
      p11str=''
      p12str=''
      p13str=''
      p14str=''
      DO isoil = 1,nsoil_w
        IF (isoil < 10) THEN 
          WRITE(layerstr,'(I1)') isoil
        ELSE
          WRITE(layerstr,'(I2)') isoil
        END IF
        WRITE(p1str,'(a)') TRIM(p1str)//' CUEmic_L'//TRIM(layerstr)
        WRITE(p2str,'(a)') TRIM(p2str)//' NUEmic_L'//TRIM(layerstr)
        WRITE(p3str,'(a)') TRIM(p3str)//' PUEmic_L'//TRIM(layerstr)
        WRITE(p4str,'(a)') TRIM(p4str)//' RootC_L'//TRIM(layerstr)
        WRITE(p5str,'(a)') TRIM(p5str)//' EnzPoly_L'//TRIM(layerstr) 
        WRITE(p6str,'(a)') TRIM(p6str)//' EnzRes_L'//TRIM(layerstr) 
        WRITE(p7str,'(a)') TRIM(p7str)//' EnzPolyC_L'//TRIM(layerstr) 
        WRITE(p8str,'(a)') TRIM(p8str)//' EnzPolyN_L'//TRIM(layerstr) 
        WRITE(p9str,'(a)') TRIM(p9str)//' EnzPolyP_L'//TRIM(layerstr) 
        WRITE(p10str,'(a)') TRIM(p10str)//' EnzPolyC_mavg_L'//TRIM(layerstr) 
        WRITE(p11str,'(a)') TRIM(p11str)//' EnzPolyN_mavg_L'//TRIM(layerstr) 
        WRITE(p12str,'(a)') TRIM(p12str)//' EnzPolyP_mavg_L'//TRIM(layerstr)
        WRITE(p13str,'(a)') TRIM(p13str)//' DOM_CN_mavg_L'//TRIM(layerstr)
        WRITE(p14str,'(a)') TRIM(p14str)//' DOM_CP_mavg_L'//TRIM(layerstr)
      END DO
      WRITE(lun_soil_diagnostics_microbial_sl,'(15a)') &
        &   'Time ',TRIM(p1str),TRIM(p2str),TRIM(p3str),TRIM(p4str),TRIM(p5str),TRIM(p6str),TRIM(p7str),&
        &   TRIM(p8str),TRIM(p9str),TRIM(p10str),TRIM(p11str),TRIM(p12str),TRIM(p13str),TRIM(p14str)
    END IF
    IF (is_active_txt_soil_diagnostics_sorption_sl) THEN
      ! name_soil_diagnostics_sorption_sl [ ... ]
      p1str=''
      p2str=''
      p3str=''
      p4str=''
      p5str=''
      p6str=''
      p7str=''
      p8str=''
      p9str=''
      DO isoil = 1,nsoil_w
        IF (isoil < 10) THEN 
          WRITE(layerstr,'(I1)') isoil
        ELSE
          WRITE(layerstr,'(I2)') isoil
        END IF
        WRITE(p1str,'(a)') TRIM(p1str)//' QmaxOrg_L'//TRIM(layerstr)
        WRITE(p2str,'(a)') TRIM(p2str)//' Rho_L'//TRIM(layerstr)
        WRITE(p3str,'(a)') TRIM(p3str)//' RhoCor_L'//TRIM(layerstr)
        WRITE(p4str,'(a)') TRIM(p4str)//' QmaxPo4_L'//TRIM(layerstr)
        WRITE(p5str,'(a)') TRIM(p5str)//' KmPo4_L'//TRIM(layerstr)
        WRITE(p9str,'(a)') TRIM(p9str)//' MinVol_L'//TRIM(layerstr)
      END DO
      WRITE(lun_soil_diagnostics_sorption_sl,'(10a)') &
        &   'Time ',TRIM(p1str),TRIM(p2str),TRIM(p3str),TRIM(p4str),TRIM(p5str),TRIM(p6str),TRIM(p7str),&
        &   TRIM(p8str),TRIM(p9str)
    END IF
    IF (is_active_txt_soilpool_dc13_org_sl .OR. &
      & is_active_txt_soilpool_dc14_org_sl .OR. &
      & is_active_txt_soilpool_dn15_org_sl) THEN
      ! prepare: name_soilpool_[dC13 DC14 dN15]_org_sl
      ! identical variable names across the below three output files
      p1str=''
      p2str=''
      p3str=''
      p4str=''
      p5str=''
      p6str=''
      p7str=''
      p8str=''
      p9str=''
      DO isoil = 1,nsoil_w
        IF (isoil < 10) THEN 
          WRITE(layerstr,'(I1)') isoil
        ELSE
          WRITE(layerstr,'(I2)') isoil
        END IF
        WRITE(p1str,'(a)') TRIM(p1str)//' SoluLit_L'//TRIM(layerstr)
        WRITE(p2str,'(a)') TRIM(p2str)//' PolyLit_L'//TRIM(layerstr)
        WRITE(p3str,'(a)') TRIM(p3str)//' WodLit_L'//TRIM(layerstr)
        WRITE(p4str,'(a)') TRIM(p4str)//' DOM_L'//TRIM(layerstr)
        WRITE(p5str,'(a)') TRIM(p5str)//' aDOM_L'//TRIM(layerstr)
        WRITE(p6str,'(a)') TRIM(p6str)//' Myc_L'//TRIM(layerstr)
        WRITE(p7str,'(a)') TRIM(p7str)//' Mic_L'//TRIM(layerstr)
        WRITE(p8str,'(a)') TRIM(p8str)//' Res_L'//TRIM(layerstr)
        WRITE(p9str,'(a)') TRIM(p9str)//' aRes_L'//TRIM(layerstr)
      END DO
    END IF
    IF (is_active_txt_soilpool_dc13_org_sl) THEN
      ! name_soilpool_dc13_org_sl [ ... ]
      WRITE(lun_soilpool_dc13_org_sl,'(11a)') &
        &   'Time ',TRIM(p1str),TRIM(p2str),TRIM(p3str),TRIM(p4str),TRIM(p5str),TRIM(p6str),TRIM(p7str),&
        &   TRIM(p8str),TRIM(p9str)
    END IF
    IF (is_active_txt_soilpool_dc14_org_sl) THEN
      ! name_soilpool_dc14_org_sl [ ... ]
      WRITE(lun_soilpool_dc14_org_sl,'(11a)') &
        &   'Time ',TRIM(p1str),TRIM(p2str),TRIM(p3str),TRIM(p4str),TRIM(p5str),TRIM(p6str),TRIM(p7str),&
        &   TRIM(p8str),TRIM(p9str)
    END IF
    IF (is_active_txt_soilpool_dn15_org_sl) THEN
      ! name_soilpool_dn15_org_sl [ ... ]
      WRITE(lun_soilpool_dn15_org_sl,'(11a)') &
        &   'Time ',TRIM(p1str),TRIM(p2str),TRIM(p3str),TRIM(p4str),TRIM(p5str),TRIM(p6str),TRIM(p7str),&
        &   TRIM(p8str),TRIM(p9str)
    END IF
  END SUBROUTINE init_output_writing_pool

  ! ======================================================================================================= !
  !>Pools: write data to output "pool" files
  !>
  !>  called by land_model_interface() at a specific "output interval"
  !>
  !>  NOTE: the max number of soil layers is assumed to be 20 (number of columns in output files are set to match this)
  !>
  SUBROUTINE write_output_pool(tile, options, time_d)

    USE mo_jsb_class,             ONLY: Get_model
    USE mo_jsb_tile_class,        ONLY: t_jsb_tile_abstract
    USE mo_jsb_task_class,        ONLY: t_jsb_task_options
    USE mo_jsb_model_class,       ONLY: t_jsb_model
    USE mo_veg_constants,         ONLY: itree
    USE mo_jsb_math_constants,    ONLY: eps8
    USE mo_isotope_util,          ONLY: calc_delta_C13, calc_delta_C14, calc_delta_N15
    USE mo_jsb_process_class,     ONLY: Q_ASSIMI_, VEG_, Q_PHENO_, SPQ_, SSE_, HYDRO_, SB_
    USE mo_jsb_grid_class,        ONLY: t_jsb_vgrid
    USE mo_jsb_grid,              ONLY: Get_vgrid
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Use_config(SB_)
    dsl4jsb_Use_memory(Q_ASSIMI_)
    dsl4jsb_Use_memory(VEG_)
    dsl4jsb_Use_memory(Q_PHENO_)
    dsl4jsb_Use_memory(SPQ_)
    dsl4jsb_Use_memory(SSE_)
    dsl4jsb_Use_memory(HYDRO_)
    dsl4jsb_Use_memory(SB_)
    ! ----------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_tile_abstract), INTENT(in) :: tile         !< one tile with data structure for one lct
    TYPE(t_jsb_task_options),   INTENT(in) :: options      !< model options
    REAL(wp),                   INTENT(in) :: time_d       !< time in days
    ! ----------------------------------------------------------------------------------------------------- !
    TYPE(t_jsb_model),      POINTER      :: model                     !< the model
    TYPE(t_lnd_bgcm_store), POINTER      :: bgcm_store                !< the bgcm store of this tile
    TYPE(t_lctlib_element), POINTER      :: lctlib                    !< land-cover-type library - parameter across pft's
    TYPE(t_jsb_vgrid),      POINTER      :: vgrid_canopy              !< Vertical grid
    TYPE(t_jsb_vgrid),      POINTER      :: vgrid_soil_w              !< Vertical grid
    INTEGER                              :: ncanopy                   !< number of canopy layers
    INTEGER                              :: nsoil_w                   !< number of soil layers as used/defined by the SB_ process
    INTEGER                              :: icanopy                   ! loop over ncanopy
    REAL(wp)                             :: leaf_mass_to_root_mass    ! name_veg_diagnostics output
    REAL(wp)                             :: dC13leaf,dC13fine_root, & !
      &                                     dC13sap_wood,dC13heart_wood,dC13coarse_root, & ! delta C13 in pools
      &                                     dC13fruit,dC13labile,dC13reserve,dC13seed_bed
    REAL(wp),ALLOCATABLE, DIMENSION(:,:) :: dC13soluable_litter,dC13polymeric_litter,dC13woody_litter, &
      &                                     dC13dom,dC13dom_assoc,dC13mycorrhiza,dC13microbial, &
      &                                     dC13residue,dC13residue_assoc
    REAL(wp),ALLOCATABLE, DIMENSION(:,:) :: hlp1_2d
    INTEGER                              :: iblk, ics, ice
    INTEGER                              :: nc, ic, is
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':write_output_pool'
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Def_mt2L2D :: veg_pool_mt
    dsl4jsb_Def_mt1L2D :: seed_bed_pool_mt
    dsl4jsb_Def_mt2L3D :: sb_pool_mt
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Def_config(SB_)
    dsl4jsb_Def_memory(Q_ASSIMI_)
    dsl4jsb_Def_memory(VEG_)
    dsl4jsb_Def_memory(Q_PHENO_)
    dsl4jsb_Def_memory(SPQ_)
    dsl4jsb_Def_memory(SSE_)
    dsl4jsb_Def_memory(HYDRO_)
    dsl4jsb_Def_memory(SB_)
    ! ----------------------------------------------------------------------------------------------------- !
    ! Q_PHENO_
    dsl4jsb_Real2D_onChunk      :: growing_season
    ! Q_ASSIMI_ 2D
    dsl4jsb_Real2D_onChunk      :: beta_soil_gs
    dsl4jsb_Real2D_onChunk      :: beta_soil_ps
    dsl4jsb_Real2D_onChunk      :: wtr_soil_root_pot
    ! Q_ASSIMI_ 2D
    dsl4jsb_Real3D_onChunk      :: jmax_cl
    dsl4jsb_Real3D_onChunk      :: vcmax_cl
    dsl4jsb_Real3D_onChunk      :: chlfl_yield_cl
    ! HYDRO_ 2D
    dsl4jsb_Real2D_onChunk      :: num_sl_above_bedrock
    dsl4jsb_Real2D_onChunk      :: wtr_rootzone
    dsl4jsb_Real2D_onChunk      :: root_depth
    dsl4jsb_Real2D_onChunk      :: weq_snow_soil
    ! HYDRO_ 3D
    dsl4jsb_Real3D_onChunk      :: soil_depth_sl
    dsl4jsb_Real3D_onChunk      :: wtr_soil_sl
    ! SPQ_ 2D
    dsl4jsb_Real2D_onChunk      :: snow_height      ! jsbach equivalent: H weq_snow_soil (Amount of snow on soil [m water equivalent])
    ! SSE_ 3D
    dsl4jsb_Real3D_onChunk      :: t_soil_sl
    dsl4jsb_Real3D_onChunk      :: bulk_dens_sl
    ! VEG_ 2D
    dsl4jsb_Real2D_onChunk      :: ph_fine_root_factor_acc
    dsl4jsb_Real2D_onChunk      :: height
    dsl4jsb_Real2D_onChunk      :: lai
    dsl4jsb_Real2D_onChunk      :: diameter
    dsl4jsb_Real2D_onChunk      :: dens_ind
    dsl4jsb_Real2D_onChunk      :: beta_sinklim_ps_tacclim_mavg
    dsl4jsb_Real2D_onChunk      :: demand_uptake_n_tuptake_mavg
    dsl4jsb_Real2D_onChunk      :: demand_uptake_p_tuptake_mavg
    dsl4jsb_Real2D_onChunk      :: growth_req_n_tuptake_mavg
    dsl4jsb_Real2D_onChunk      :: growth_req_p_tuptake_mavg
    dsl4jsb_Real2D_onChunk      :: veg_pool_total_c
    dsl4jsb_Real2D_onChunk      :: veg_pool_total_n
    dsl4jsb_Real2D_onChunk      :: veg_pool_total_p
    dsl4jsb_Real2D_onChunk      :: veg_pool_total_c13
    dsl4jsb_Real2D_onChunk      :: veg_pool_total_c14
    dsl4jsb_Real2D_onChunk      :: veg_pool_total_n15
    ! VEG_ 3D
    dsl4jsb_Real3D_onChunk      :: fn_et_cl
    dsl4jsb_Real3D_onChunk      :: fn_rub_cl
    dsl4jsb_Real3D_onChunk      :: fn_chl_cl
    dsl4jsb_Real3D_onChunk      :: fn_pepc_cl
    dsl4jsb_Real3D_onChunk      :: fn_oth_cl
    dsl4jsb_Real3D_onChunk      :: leaf_nitrogen_cl
    dsl4jsb_Real3D_onChunk      :: root_fraction_sl
    ! SB_ 3D
    dsl4jsb_Real3D_onChunk      :: rph_decomposition
    dsl4jsb_Real3D_onChunk      :: ph_sl
    dsl4jsb_Real3D_onChunk      :: nh4_solute
    dsl4jsb_Real3D_onChunk      :: nh4_n15_solute
    dsl4jsb_Real3D_onChunk      :: no3_solute
    dsl4jsb_Real3D_onChunk      :: no3_n15_solute
    dsl4jsb_Real3D_onChunk      :: po4_solute
    dsl4jsb_Real3D_onChunk      :: noy
    dsl4jsb_Real3D_onChunk      :: noy_n15
    dsl4jsb_Real3D_onChunk      :: n2o
    dsl4jsb_Real3D_onChunk      :: n2o_n15
    dsl4jsb_Real3D_onChunk      :: n2
    dsl4jsb_Real3D_onChunk      :: n2_n15
    dsl4jsb_Real3D_onChunk      :: nh4_assoc
    dsl4jsb_Real3D_onChunk      :: nh4_n15_assoc
    dsl4jsb_Real3D_onChunk      :: po4_assoc_fast
    dsl4jsb_Real3D_onChunk      :: po4_assoc_slow
    dsl4jsb_Real3D_onChunk      :: po4_occluded
    dsl4jsb_Real3D_onChunk      :: po4_primary
    dsl4jsb_Real3D_onChunk      :: qmax_org
    dsl4jsb_Real3D_onChunk      :: qmax_po4
    dsl4jsb_Real3D_onChunk      :: km_adsorpt_po4_sl
    dsl4jsb_Real3D_onChunk      :: bulk_dens_corr_sl
    dsl4jsb_Real3D_onChunk      :: microbial_cue_eff_tmic_mavg 
    dsl4jsb_Real3D_onChunk      :: microbial_nue_eff_tmic_mavg
    dsl4jsb_Real3D_onChunk      :: microbial_pue_eff_tmic_mavg
    dsl4jsb_Real3D_onChunk      :: dom_cn_mavg 
    dsl4jsb_Real3D_onChunk      :: dom_cp_mavg 
    dsl4jsb_Real3D_onChunk      :: enzyme_frac_poly
    dsl4jsb_Real3D_onChunk      :: enzyme_frac_residue
    dsl4jsb_Real3D_onChunk      :: enzyme_frac_poly_c
    dsl4jsb_Real3D_onChunk      :: enzyme_frac_poly_n
    dsl4jsb_Real3D_onChunk      :: enzyme_frac_poly_p
    dsl4jsb_Real3D_onChunk      :: enzyme_frac_poly_c_mavg
    dsl4jsb_Real3D_onChunk      :: enzyme_frac_poly_n_mavg
    dsl4jsb_Real3D_onChunk      :: enzyme_frac_poly_p_mavg
    dsl4jsb_Real3D_onChunk      :: sb_pool_total_c
    dsl4jsb_Real3D_onChunk      :: volume_min_sl
    ! ----------------------------------------------------------------------------------------------------- !
    iblk      = options%iblk
    ics       = options%ics
    ice       = options%ice
    nc        = options%nc
    ! ----------------------------------------------------------------------------------------------------- !
    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')
    ! ----------------------------------------------------------------------------------------------------- !
    model         => Get_model(tile%owner_model_id)
    lctlib        => model%lctlib(tile%lcts(1)%lib_id)
    vgrid_canopy  => Get_vgrid('q_canopy_layer')
    ncanopy       =  vgrid_canopy%n_levels
    vgrid_soil_w  => Get_vgrid('soil_depth_water')
    nsoil_w       =  vgrid_soil_w%n_levels
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Get_config(SB_)
    dsl4jsb_Get_memory(Q_ASSIMI_)
    dsl4jsb_Get_memory(VEG_)
    dsl4jsb_Get_memory(Q_PHENO_)
    dsl4jsb_Get_memory(SPQ_)
    dsl4jsb_Get_memory(SSE_)
    dsl4jsb_Get_memory(HYDRO_)
    dsl4jsb_Get_memory(SB_)
    ! ----------------------------------------------------------------------------------------------------- !
    ALLOCATE(dC13soluable_litter(nc, nsoil_w))
    ALLOCATE(dC13polymeric_litter(nc, nsoil_w))
    ALLOCATE(dC13woody_litter(nc, nsoil_w))
    ALLOCATE(dC13dom(nc, nsoil_w))
    ALLOCATE(dC13dom_assoc(nc, nsoil_w))
    ALLOCATE(dC13mycorrhiza(nc, nsoil_w))
    ALLOCATE(dC13microbial(nc, nsoil_w))
    ALLOCATE(dC13residue(nc, nsoil_w))
    ALLOCATE(dC13residue_assoc(nc, nsoil_w))
    ALLOCATE(hlp1_2d(nc, nsoil_w))
    ! ----------------------------------------------------------------------------------------------------- !
    bgcm_store => tile%bgcm_store
    dsl4jsb_Get_mt2L2D(VEG_BGCM_POOL_ID, veg_pool_mt)
    dsl4jsb_Get_mt1L2D(VEG_BGCM_POOL_ID, seed_bed_pool_mt)
    dsl4jsb_Get_mt2L3D(SB_BGCM_POOL_ID, sb_pool_mt)
    ! ----------------------------------------------------------------------------------------------------- !
    ! Q_PHENO_ 2D
    dsl4jsb_Get_var2D_onChunk(Q_PHENO_, growing_season)
    ! Q_ASSIMI_ 2D
    dsl4jsb_Get_var2D_onChunk(Q_ASSIMI_, beta_soil_gs)
    dsl4jsb_Get_var2D_onChunk(Q_ASSIMI_, beta_soil_ps)
    dsl4jsb_Get_var2D_onChunk(Q_ASSIMI_, wtr_soil_root_pot)
    ! Q_ASSIMI_ 3D
    dsl4jsb_Get_var3D_onChunk(Q_ASSIMI_, jmax_cl)
    dsl4jsb_Get_var3D_onChunk(Q_ASSIMI_, vcmax_cl)
    dsl4jsb_Get_var3D_onChunk(Q_ASSIMI_, chlfl_yield_cl)
    ! HYDRO_ 2D
    dsl4jsb_Get_var2D_onChunk(HYDRO_,    num_sl_above_bedrock)
    dsl4jsb_Get_var2D_onChunk(HYDRO_,    wtr_rootzone)
    dsl4jsb_Get_var2D_onChunk(HYDRO_,    root_depth)
    dsl4jsb_Get_var2D_onChunk(HYDRO_,    weq_snow_soil)
    ! HYDRO_ 3D
    dsl4jsb_Get_var3D_onChunk(HYDRO_,    soil_depth_sl)
    dsl4jsb_Get_var3D_onChunk(HYDRO_,    wtr_soil_sl)
    ! SPQ_ 2D
    dsl4jsb_Get_var2D_onChunk(SPQ_, snow_height)
    ! SSE_ 3D
    dsl4jsb_Get_var3D_onChunk(SSE_, t_soil_sl)
    dsl4jsb_Get_var3D_onChunk(SSE_, bulk_dens_sl)
    ! ---------------------------
    dsl4jsb_Get_var2D_onChunk(VEG_, ph_fine_root_factor_acc)
    dsl4jsb_Get_var2D_onChunk(VEG_, height)
    dsl4jsb_Get_var2D_onChunk(VEG_, lai)
    dsl4jsb_Get_var2D_onChunk(VEG_, diameter)
    dsl4jsb_Get_var2D_onChunk(VEG_, dens_ind)
    dsl4jsb_Get_var2D_onChunk(VEG_, beta_sinklim_ps_tacclim_mavg)
    dsl4jsb_Get_var2D_onChunk(VEG_, demand_uptake_n_tuptake_mavg)
    dsl4jsb_Get_var2D_onChunk(VEG_, demand_uptake_p_tuptake_mavg)
    dsl4jsb_Get_var2D_onChunk(VEG_, growth_req_n_tuptake_mavg)
    dsl4jsb_Get_var2D_onChunk(VEG_, growth_req_p_tuptake_mavg)
    dsl4jsb_Get_var2D_onChunk(VEG_, veg_pool_total_c)
    dsl4jsb_Get_var2D_onChunk(VEG_, veg_pool_total_n)
    dsl4jsb_Get_var2D_onChunk(VEG_, veg_pool_total_p)
    dsl4jsb_Get_var2D_onChunk(VEG_, veg_pool_total_c13)
    dsl4jsb_Get_var2D_onChunk(VEG_, veg_pool_total_c14)
    dsl4jsb_Get_var2D_onChunk(VEG_, veg_pool_total_n15)
    ! VEG_ 3D
    dsl4jsb_Get_var3D_onChunk(VEG_, fn_et_cl)
    dsl4jsb_Get_var3D_onChunk(VEG_, fn_rub_cl)
    dsl4jsb_Get_var3D_onChunk(VEG_, fn_chl_cl)
    dsl4jsb_Get_var3D_onChunk(VEG_, fn_pepc_cl)
    dsl4jsb_Get_var3D_onChunk(VEG_, fn_oth_cl)
    dsl4jsb_Get_var3D_onChunk(VEG_, leaf_nitrogen_cl)
    dsl4jsb_Get_var3D_onChunk(VEG_, root_fraction_sl)
    ! ---------------------------
    dsl4jsb_Get_var3D_onChunk(SB_,  rph_decomposition)
    dsl4jsb_Get_var3D_onChunk(SB_,  ph_sl)
    dsl4jsb_Get_var3D_onChunk(SB_,  nh4_solute)
    dsl4jsb_Get_var3D_onChunk(SB_,  nh4_n15_solute)
    dsl4jsb_Get_var3D_onChunk(SB_,  no3_solute)
    dsl4jsb_Get_var3D_onChunk(SB_,  no3_n15_solute)
    dsl4jsb_Get_var3D_onChunk(SB_,  po4_solute)
    dsl4jsb_Get_var3D_onChunk(SB_,  noy)
    dsl4jsb_Get_var3D_onChunk(SB_,  noy_n15)
    dsl4jsb_Get_var3D_onChunk(SB_,  n2o)
    dsl4jsb_Get_var3D_onChunk(SB_,  n2o_n15)
    dsl4jsb_Get_var3D_onChunk(SB_,  n2)
    dsl4jsb_Get_var3D_onChunk(SB_,  n2_n15)
    dsl4jsb_Get_var3D_onChunk(SB_,  nh4_assoc)
    dsl4jsb_Get_var3D_onChunk(SB_,  nh4_n15_assoc)
    dsl4jsb_Get_var3D_onChunk(SB_,  po4_assoc_fast)
    dsl4jsb_Get_var3D_onChunk(SB_,  po4_assoc_slow)
    dsl4jsb_Get_var3D_onChunk(SB_,  po4_occluded)
    dsl4jsb_Get_var3D_onChunk(SB_,  po4_primary)
    dsl4jsb_Get_var3D_onChunk(SB_,  qmax_org)
    dsl4jsb_Get_var3D_onChunk(SB_,  qmax_po4)
    dsl4jsb_Get_var3D_onChunk(SB_,  km_adsorpt_po4_sl)
    dsl4jsb_Get_var3D_onChunk(SB_,  bulk_dens_corr_sl)
    dsl4jsb_Get_var3D_onChunk(SB_,  microbial_cue_eff_tmic_mavg)
    dsl4jsb_Get_var3D_onChunk(SB_,  microbial_nue_eff_tmic_mavg)
    dsl4jsb_Get_var3D_onChunk(SB_,  microbial_pue_eff_tmic_mavg)
    dsl4jsb_Get_var3D_onChunk(SB_,  dom_cn_mavg)
    dsl4jsb_Get_var3D_onChunk(SB_,  dom_cp_mavg)
    dsl4jsb_Get_var3D_onChunk(SB_,  enzyme_frac_poly)
    dsl4jsb_Get_var3D_onChunk(SB_,  enzyme_frac_residue)
    dsl4jsb_Get_var3D_onChunk(SB_,  enzyme_frac_poly_c)
    dsl4jsb_Get_var3D_onChunk(SB_,  enzyme_frac_poly_n)
    dsl4jsb_Get_var3D_onChunk(SB_,  enzyme_frac_poly_p)
    dsl4jsb_Get_var3D_onChunk(SB_,  enzyme_frac_poly_c_mavg)
    dsl4jsb_Get_var3D_onChunk(SB_,  enzyme_frac_poly_n_mavg)
    dsl4jsb_Get_var3D_onChunk(SB_,  enzyme_frac_poly_p_mavg)
    dsl4jsb_Get_var3D_onChunk(SB_,  sb_pool_total_c)
    dsl4jsb_Get_var3D_onChunk(SB_,  volume_min_sl)
    ! ----------------------------------------------------------------------------------------------------- !

    !>0.9 the point in ics:ice of which data will be written to file; per default = 1
    !>
    !>  is actually defined and declared globally in his module above, but the ifort compiler complains anyway
    !>
    ic = 1

    !>1.0 write output
    !>
    IF (is_active_txt_vegpool_c) THEN
      ! name_vegpool_c
      write(lun_vegpool_c,'(13F12.6)') time_d, &
        &   veg_pool_mt(ix_leaf, ixC, ic),&
        &   veg_pool_mt(ix_fine_root, ixC, ic),veg_pool_mt(ix_sap_wood, ixC, ic),&
        &   veg_pool_mt(ix_heart_wood, ixC, ic),veg_pool_mt(ix_coarse_root, ixC, ic),&
        &   veg_pool_mt(ix_fruit, ixC, ic),veg_pool_mt(ix_labile, ixC, ic),&
        &   veg_pool_mt(ix_reserve, ixC, ic),seed_bed_pool_mt(ixC, ic),&
        &   veg_pool_total_c(ic), veg_pool_total_c13(ic), veg_pool_total_c14(ic)
    END IF
    IF (is_active_txt_vegpool_dc13dc14) THEN
      ! name_vegpool_dc13dc14
      dC13leaf        = calc_delta_C13(veg_pool_mt(ix_leaf, ixC, ic),veg_pool_mt(ix_leaf, ixC13, ic))
      dC13fine_root   = calc_delta_C13(veg_pool_mt(ix_fine_root, ixC, ic),veg_pool_mt(ix_fine_root, ixC13, ic))
      dC13sap_wood    = calc_delta_C13(veg_pool_mt(ix_sap_wood, ixC, ic),veg_pool_mt(ix_sap_wood, ixC13, ic))
      dC13heart_wood  = calc_delta_C13(veg_pool_mt(ix_heart_wood, ixC, ic),veg_pool_mt(ix_heart_wood, ixC13, ic))
      dC13coarse_root = calc_delta_C13(veg_pool_mt(ix_coarse_root, ixC, ic),veg_pool_mt(ix_coarse_root, ixC13, ic))
      dC13fruit       = calc_delta_C13(veg_pool_mt(ix_fruit, ixC, ic),veg_pool_mt(ix_fruit, ixC13, ic))
      dC13labile      = calc_delta_C13(veg_pool_mt(ix_labile, ixC, ic),veg_pool_mt(ix_labile, ixC13, ic))
      dC13reserve     = calc_delta_C13(veg_pool_mt(ix_reserve, ixC, ic),veg_pool_mt(ix_reserve, ixC13, ic))
      dC13seed_bed    = calc_delta_C13(seed_bed_pool_mt(ixC, ic),seed_bed_pool_mt(ixC13, ic))
      write(lun_vegpool_dc13dc14,'(19F18.6)') time_d, &
        &   dC13leaf,dC13fine_root,dC13sap_wood,dC13heart_wood,&
        &   dC13coarse_root,dC13fruit,dC13labile,dC13reserve,dC13seed_bed, &
        &   calc_delta_C14(veg_pool_mt(ix_leaf, ixC, ic),veg_pool_mt(ix_leaf, ixC14, ic),dC13leaf),&
        &   calc_delta_C14(veg_pool_mt(ix_fine_root, ixC, ic),veg_pool_mt(ix_fine_root, ixC14, ic),dC13fine_root),&
        &   calc_delta_C14(veg_pool_mt(ix_sap_wood, ixC, ic),veg_pool_mt(ix_sap_wood, ixC14, ic),dC13sap_wood),&
        &   calc_delta_C14(veg_pool_mt(ix_heart_wood, ixC, ic),veg_pool_mt(ix_heart_wood, ixC14, ic),dC13heart_wood),&
        &   calc_delta_C14(veg_pool_mt(ix_coarse_root, ixC, ic),veg_pool_mt(ix_coarse_root, ixC14, ic),dC13coarse_root),&
        &   calc_delta_C14(veg_pool_mt(ix_fruit, ixC, ic),veg_pool_mt(ix_fruit, ixC14, ic),dC13fruit),&
        &   calc_delta_C14(veg_pool_mt(ix_labile, ixC, ic),veg_pool_mt(ix_labile, ixC14, ic),dC13labile),&
        &   calc_delta_C14(veg_pool_mt(ix_reserve, ixC, ic),veg_pool_mt(ix_reserve, ixC14, ic),dC13reserve),&
        &   calc_delta_C14(seed_bed_pool_mt(ixC, ic),seed_bed_pool_mt(ixC14, ic),dC13seed_bed)
    END IF
    IF (is_active_txt_vegpool_n) THEN
      ! name_vegpool_n
      write(lun_vegpool_n,'(12F12.6)') time_d, &
        &   veg_pool_mt(ix_leaf, ixN, ic),&
        &   veg_pool_mt(ix_fine_root, ixN, ic),veg_pool_mt(ix_sap_wood, ixN, ic),&
        &   veg_pool_mt(ix_heart_wood, ixN, ic),veg_pool_mt(ix_coarse_root, ixN, ic),&
        &   veg_pool_mt(ix_fruit, ixN, ic),veg_pool_mt(ix_labile, ixN, ic),&
        &   veg_pool_mt(ix_reserve, ixN, ic),seed_bed_pool_mt(ixN, ic),&
        &   veg_pool_total_n(ic)
    END IF
    IF (is_active_txt_vegpool_p) THEN
      ! name_vegpool_p
      write(lun_vegpool_p,'(11F12.6)') time_d, &
        &   veg_pool_mt(ix_leaf, ixP, ic),&
        &   veg_pool_mt(ix_fine_root, ixP, ic),veg_pool_mt(ix_sap_wood, ixP, ic),&
        &   veg_pool_mt(ix_heart_wood, ixP, ic),veg_pool_mt(ix_coarse_root, ixP, ic),&
        &   veg_pool_mt(ix_fruit, ixP, ic),veg_pool_mt(ix_labile, ixP, ic),&
        &   veg_pool_mt(ix_reserve, ixP, ic),seed_bed_pool_mt(ixP, ic),&
        &   veg_pool_total_p(ic)
    END IF
    IF (is_active_txt_vegpool_dn15) THEN
      ! name_vegpool_dn15
      write(lun_vegpool_dn15,'(12F14.6)') time_d, &
        &   calc_delta_N15(veg_pool_mt(ix_leaf, ixN, ic),veg_pool_mt(ix_leaf, ixN15, ic)), &
        &   calc_delta_N15(veg_pool_mt(ix_fine_root, ixN, ic),veg_pool_mt(ix_fine_root, ixN15, ic)), &
        &   calc_delta_N15(veg_pool_mt(ix_sap_wood, ixN, ic),veg_pool_mt(ix_sap_wood, ixN15, ic)), &
        &   calc_delta_N15(veg_pool_mt(ix_heart_wood, ixN, ic),veg_pool_mt(ix_heart_wood, ixN15, ic)), &
        &   calc_delta_N15(veg_pool_mt(ix_coarse_root, ixN, ic),veg_pool_mt(ix_coarse_root, ixN15, ic)), &
        &   calc_delta_N15(veg_pool_mt(ix_fruit, ixN, ic),veg_pool_mt(ix_fruit, ixN15, ic)), &
        &   calc_delta_N15(veg_pool_mt(ix_labile, ixN, ic),veg_pool_mt(ix_labile, ixN15, ic)), &
        &   calc_delta_N15(veg_pool_mt(ix_reserve, ixN, ic),veg_pool_mt(ix_reserve, ixN15, ic)), &
        &   calc_delta_N15(seed_bed_pool_mt(ixN, ic),seed_bed_pool_mt(ixN15, ic)), &
        &   calc_delta_N15(veg_pool_total_n(ic), veg_pool_total_n15(ic))
    END IF
    IF (is_active_txt_soil_physics) THEN
      ! name_soil_physics
      ! output jsbach variable for snow height when running with "jsbach soil-physics"
      IF (model%config%use_soil_phys_jsbach) THEN
        snow_height(ic) = weq_snow_soil(ic)
      END IF
      write(lun_soil_physics,'(20F16.6)') time_d, &
        &   snow_height(ic), wtr_rootzone(ic), wtr_soil_root_pot(ic),ph_sl(ic,1), rph_decomposition(ic,1)
    END IF
    IF (is_active_txt_soil_physics_sl) THEN
      ! name_soil_physics_sl
      write(lun_soil_physics_sl,'(41F16.6)') time_d, wtr_soil_sl(ic,:), t_soil_sl(ic,:)
    END IF
    IF (is_active_txt_soilpool_inorg) THEN
      ! name_soilpool_inorg
      write(lun_soilpool_inorg,'(30F16.6)') time_d, &
        &   SUM(nh4_assoc(ic,:)*soil_depth_sl(ic,:))    , &
        &   SUM(nh4_solute(ic,:)*soil_depth_sl(ic,:))   , &
        &   SUM(no3_solute(ic,:)*soil_depth_sl(ic,:))   , &
        &   SUM(po4_solute(ic,:)*soil_depth_sl(ic,:))   , &
        &   SUM(noy(ic,:)*soil_depth_sl(ic,:))          , &
        &   SUM(n2o(ic,:)*soil_depth_sl(ic,:))          , &
        &   SUM(n2(ic,:)*soil_depth_sl(ic,:))
    END IF
    IF (is_active_txt_ecosyspool_cnp) THEN
      ! name_ecosyspool_cnp
      write(lun_ecosyspool_cnp,'(30F16.6)') time_d, &
        ! TotalSoilOrgC
        &   SUM(sb_pool_mt(ix_soluable_litter, ixC, ic, :)*soil_depth_sl(ic,:))     + & 
        &    SUM(sb_pool_mt(ix_polymeric_litter, ixC, ic, :)*soil_depth_sl(ic,:))     + & 
        &    SUM(sb_pool_mt(ix_woody_litter, ixC, ic, :)*soil_depth_sl(ic,:))         + &
        &    SUM(sb_pool_mt(ix_microbial, ixC, ic, :)*soil_depth_sl(ic,:))            + &
        &    SUM(sb_pool_mt(ix_residue, ixC, ic, :)*soil_depth_sl(ic,:))              + & 
        &    SUM(sb_pool_mt(ix_mycorrhiza, ixC, ic, :)*soil_depth_sl(ic,:))           , &
        ! TotalSoilOrgN
        &   SUM(sb_pool_mt(ix_soluable_litter, ixN, ic, :)*soil_depth_sl(ic,:))   + & 
        &    SUM(sb_pool_mt(ix_polymeric_litter, ixN, ic, :)*soil_depth_sl(ic,:))   + & 
        &    SUM(sb_pool_mt(ix_woody_litter, ixN, ic, :)*soil_depth_sl(ic,:))       + &
        &    SUM(sb_pool_mt(ix_microbial, ixN, ic, :)*soil_depth_sl(ic,:))          + &
        &    SUM(sb_pool_mt(ix_residue, ixN, ic, :)*soil_depth_sl(ic,:))            + &
        &    SUM(sb_pool_mt(ix_mycorrhiza, ixN, ic, :)*soil_depth_sl(ic,:))         , &
        ! TotalSoilInOrgN
        &   SUM(nh4_solute(ic,:)*soil_depth_sl(ic,:))                                + & 
        &    SUM(nh4_assoc(ic,:)*soil_depth_sl(ic,:))                                + & 
        &    SUM(no3_solute(ic,:)*soil_depth_sl(ic,:))                                + & 
        &    SUM(noy(ic,:)*soil_depth_sl(ic,:))                                       + &
        &    SUM(n2o(ic,:)*soil_depth_sl(ic,:))                                       + &
        &    SUM(n2(ic,:)*soil_depth_sl(ic,:))                                        , &
        ! EcoC
        &   SUM(sb_pool_total_c(ic,:) * soil_depth_sl(ic,:)) + veg_pool_total_c(ic), &
        ! EcoN
        ! (potentially missing some pools used with jsm)
        &   SUM(sb_pool_mt(ix_soluable_litter, ixN, ic, :)*soil_depth_sl(ic,:))   + & 
        &    SUM(sb_pool_mt(ix_polymeric_litter, ixN, ic, :)*soil_depth_sl(ic,:))   + & 
        &    SUM(sb_pool_mt(ix_woody_litter, ixN, ic, :)*soil_depth_sl(ic,:))       + &
        &    SUM(sb_pool_mt(ix_microbial, ixN, ic, :)*soil_depth_sl(ic,:))          + &
        &    SUM(sb_pool_mt(ix_residue, ixN, ic, :)*soil_depth_sl(ic,:))            + & 
        &    SUM(sb_pool_mt(ix_mycorrhiza, ixN, ic, :)*soil_depth_sl(ic,:))         + &
        &    SUM(nh4_solute(ic,:)*soil_depth_sl(ic,:))                                + & 
        &    SUM(nh4_assoc(ic,:)*soil_depth_sl(ic,:))                                + & 
        &    SUM(no3_solute(ic,:)*soil_depth_sl(ic,:))                                + &
        &    SUM(noy(ic,:)*soil_depth_sl(ic,:))                                       + &
        &    SUM(n2o(ic,:)*soil_depth_sl(ic,:))                                       + &
        &    SUM(n2(ic,:)*soil_depth_sl(ic,:))                                        + &     
        &    veg_pool_total_n(ic)                                       , &
        ! EcoP
        ! (potentially missing some pools used with jsm)
        &   SUM(sb_pool_mt(ix_soluable_litter, ixP, ic, :)*soil_depth_sl(ic,:))   + & 
        &    SUM(sb_pool_mt(ix_polymeric_litter, ixP, ic, :)*soil_depth_sl(ic,:))   + & 
        &    SUM(sb_pool_mt(ix_woody_litter, ixP, ic, :)*soil_depth_sl(ic,:))       + &
        &    SUM(sb_pool_mt(ix_microbial, ixP, ic, :)*soil_depth_sl(ic,:))          + &
        &    SUM(sb_pool_mt(ix_residue, ixP, ic, :)*soil_depth_sl(ic,:))            + & 
        &    SUM(sb_pool_mt(ix_mycorrhiza, ixP, ic, :)*soil_depth_sl(ic,:))         + &
        &    SUM(po4_solute(ic,:)*soil_depth_sl(ic,:))                            + & 
        &    SUM(po4_assoc_fast(ic,:)*soil_depth_sl(ic,:))                        + &
        &    SUM(po4_assoc_slow(ic,:)*soil_depth_sl(ic,:))                        + &
        &    SUM(po4_occluded(ic,:)*soil_depth_sl(ic,:))                          + &
        &    SUM(po4_primary(ic,:)*soil_depth_sl(ic,:))                           + &
        &    veg_pool_total_p(ic)                                     , &
        ! EcoN15
        ! (potentially missing some pools used with jsm)
        &   SUM(sb_pool_mt(ix_soluable_litter, ixN15, ic, :)*soil_depth_sl(ic,:))   + & 
        &    SUM(sb_pool_mt(ix_polymeric_litter, ixN15, ic, :)*soil_depth_sl(ic,:))   + & 
        &    SUM(sb_pool_mt(ix_woody_litter, ixN15, ic, :)*soil_depth_sl(ic,:))       + &
        &    SUM(sb_pool_mt(ix_microbial, ixN15, ic, :)*soil_depth_sl(ic,:))          + &
        &    SUM(sb_pool_mt(ix_residue, ixN15, ic, :)*soil_depth_sl(ic,:))            + & 
        &    SUM(sb_pool_mt(ix_mycorrhiza, ixN15, ic, :)*soil_depth_sl(ic,:))         + &
        &    SUM(nh4_n15_assoc(ic,:)*soil_depth_sl(ic,:))                                + & 
        &    SUM(nh4_n15_solute(ic,:)*soil_depth_sl(ic,:))                                + & 
        &    SUM(no3_n15_solute(ic,:)*soil_depth_sl(ic,:))                                + &
        &    SUM(noy_n15(ic,:)*soil_depth_sl(ic,:))                                       + &
        &    SUM(n2o_n15(ic,:)*soil_depth_sl(ic,:))                                       + &
        &    SUM(n2_n15(ic,:)*soil_depth_sl(ic,:))                                        + &     
        &    veg_pool_total_n15(ic)
    END IF
    IF (is_active_txt_soilpool_dc13dc14dn15) THEN
      ! name_soilpool_dc13dc14dn15
      dC13soluable_litter(ic,1)  = calc_delta_C13(SUM(sb_pool_mt(ix_soluable_litter, ixC, ic, :) &
        &   * soil_depth_sl(ic,:)), &
        &   SUM(sb_pool_mt(ix_soluable_litter, ixC13, ic, :)*soil_depth_sl(ic,:)))
      dC13polymeric_litter(ic,1) = calc_delta_C13(SUM(sb_pool_mt(ix_polymeric_litter, ixC, ic, :) &
        &   * soil_depth_sl(ic,:)), &
        &   SUM(sb_pool_mt(ix_polymeric_litter, ixC13, ic, :)*soil_depth_sl(ic,:)))
      dC13woody_litter(ic,1)     = calc_delta_C13(SUM(sb_pool_mt(ix_woody_litter, ixC, ic, :) &
        &   * soil_depth_sl(ic,:)), &
        &   SUM(sb_pool_mt(ix_woody_litter, ixC13, ic, :)*soil_depth_sl(ic,:)))
      dC13microbial(ic,1)        = calc_delta_C13(SUM(sb_pool_mt(ix_microbial, ixC, ic, :) &
        &   * soil_depth_sl(ic,:)), &
        &   SUM(sb_pool_mt(ix_microbial, ixC13, ic, :)*soil_depth_sl(ic,:)))
      dC13residue(ic,1)          = calc_delta_C13(SUM(sb_pool_mt(ix_residue, ixC, ic, :) &
        &   * soil_depth_sl(ic,:)), &
        &   SUM(sb_pool_mt(ix_residue, ixC13, ic, :)*soil_depth_sl(ic,:)))
      dC13mycorrhiza(ic,1)       = calc_delta_C13(SUM(sb_pool_mt(ix_mycorrhiza, ixC, ic, :) &
        &   * soil_depth_sl(ic,:)), &
        &   SUM(sb_pool_mt(ix_mycorrhiza, ixC13, ic, :)*soil_depth_sl(ic,:)))
      write(lun_soilpool_dc13dc14dn15,'(21F14.6)') time_d, &
        &   dC13soluable_litter(ic,1), &
        &   dC13polymeric_litter(ic,1), &
        &   dC13woody_litter(ic,1), &
        &   dC13microbial(ic,1), &
        &   dC13residue(ic,1), &
        &   dC13mycorrhiza(ic,1), &
        &   calc_delta_C14(SUM(sb_pool_mt(ix_soluable_litter, ixC, ic, :)*soil_depth_sl(ic,:)), &
        &                  SUM(sb_pool_mt(ix_soluable_litter, ixC14, ic, :)*soil_depth_sl(ic,:)), &
        &                       dC13soluable_litter(ic,1)) , &
        &   calc_delta_C14(SUM(sb_pool_mt(ix_polymeric_litter, ixC, ic, :)*soil_depth_sl(ic,:)), &
        &                  SUM(sb_pool_mt(ix_polymeric_litter, ixC14, ic, :)*soil_depth_sl(ic,:)), &
        &                       dC13polymeric_litter(ic,1)), &
        &   calc_delta_C14(SUM(sb_pool_mt(ix_woody_litter, ixC, ic, :)*soil_depth_sl(ic,:)), &
        &                  SUM(sb_pool_mt(ix_woody_litter, ixC14, ic, :)*soil_depth_sl(ic,:)), &
        &                       dC13woody_litter(ic,1))        , &
        &   calc_delta_C14(SUM(sb_pool_mt(ix_microbial, ixC, ic, :)*soil_depth_sl(ic,:)), &
        &                  SUM(sb_pool_mt(ix_microbial, ixC14, ic, :)*soil_depth_sl(ic,:)), dC13microbial(ic,1)), &
        &   calc_delta_C14(SUM(sb_pool_mt(ix_residue, ixC, ic, :)*soil_depth_sl(ic,:)), &
        &                  SUM(sb_pool_mt(ix_residue, ixC14, ic, :)*soil_depth_sl(ic,:)), dC13residue(ic,1)), &
        &   calc_delta_C14(SUM(sb_pool_mt(ix_mycorrhiza, ixC, ic, :)*soil_depth_sl(ic,:)), &
        &                  SUM(sb_pool_mt(ix_mycorrhiza, ixC14, ic, :)*soil_depth_sl(ic,:)),dC13mycorrhiza(ic,1)), &
        &   calc_delta_N15(SUM(sb_pool_mt(ix_soluable_litter, ixN, ic, :)*soil_depth_sl(ic,:)),&
        &                  SUM(sb_pool_mt(ix_soluable_litter, ixN15, ic, :)*soil_depth_sl(ic,:))), &
        &   calc_delta_N15(SUM(sb_pool_mt(ix_polymeric_litter, ixN, ic, :)*soil_depth_sl(ic,:)),&
        &                  SUM(sb_pool_mt(ix_polymeric_litter, ixN15, ic, :)*soil_depth_sl(ic,:))), &
        &   calc_delta_N15(SUM(sb_pool_mt(ix_woody_litter, ixN, ic, :)*soil_depth_sl(ic,:)), &
        &                  SUM(sb_pool_mt(ix_woody_litter, ixN15, ic, :)*soil_depth_sl(ic,:))), &
        &   calc_delta_N15(SUM(sb_pool_mt(ix_microbial, ixN, ic, :)*soil_depth_sl(ic,:)), &
        &                  SUM(sb_pool_mt(ix_microbial, ixN15, ic, :)*soil_depth_sl(ic,:))), &
        &   calc_delta_N15(SUM(sb_pool_mt(ix_residue, ixN, ic, :)*soil_depth_sl(ic,:)), &
        &                  SUM(sb_pool_mt(ix_residue, ixN15, ic, :)*soil_depth_sl(ic,:))), &
        &   calc_delta_N15(SUM(nh4_assoc(ic,:)*soil_depth_sl(ic,:)), &
        &                  SUM(nh4_n15_assoc(ic,:)*soil_depth_sl(ic,:))), &
        &   calc_delta_N15(SUM(nh4_solute(ic,:)*soil_depth_sl(ic,:)), &
        &                  SUM(nh4_n15_solute(ic,:)*soil_depth_sl(ic,:))), &
        &   calc_delta_N15(SUM(no3_solute(ic,:)*soil_depth_sl(ic,:)), &
        &                  SUM(no3_n15_solute(ic,:)*soil_depth_sl(ic,:)))
    END IF
    IF (is_active_txt_veg_diagnostics) THEN
      ! name_veg_diagnostics
      IF (veg_pool_mt(ix_fine_root, ixC, ic) > eps8) THEN
        leaf_mass_to_root_mass = veg_pool_mt(ix_leaf, ixC, ic) / veg_pool_mt(ix_fine_root, ixC, ic)
      ELSE
        leaf_mass_to_root_mass = 0._wp
      END IF
      write(lun_veg_diagnostics,'(15F24.8)') time_d, &
        &   height(ic), diameter(ic), lai(ic), dens_ind(ic), leaf_mass_to_root_mass, &
        &   growing_season(ic), demand_uptake_n_tuptake_mavg(ic),  demand_uptake_p_tuptake_mavg(ic), &
        &   beta_sinklim_ps_tacclim_mavg(ic), &
        &   growth_req_n_tuptake_mavg(ic), growth_req_p_tuptake_mavg(ic), &
        &   beta_soil_gs(ic), beta_soil_ps(ic),ph_fine_root_factor_acc(ic)
    END IF
    IF (is_active_txt_veg_canopy_diagnostics_cl) THEN
      ! name_veg_canopy_diagnostics_cl
      write(lun_veg_canopy_diagnostics_cl,'(136F18.8)') &
        &   time_d, &
        &   fn_et_cl(ic,:), &
        &   fn_rub_cl(ic,:), &
        &   fn_chl_cl(ic,:), &
        &   fn_pepc_cl(ic,:), &
        &   fn_oth_cl(ic,:), &
        &   leaf_nitrogen_cl(ic,:), &
        &   jmax_cl(ic,:), &
        &   vcmax_cl(ic,:), &
        &   chlfl_yield_cl(ic,:)
    END IF
    IF (is_active_txt_soilpool_org) THEN
      ! name_soilpool_org
      write(lun_soilpool_org,'(28F24.8)') time_d, &
        ! C
        &   SUM(sb_pool_mt(ix_dom, ixC, ic, :) * soil_depth_sl(ic,:)), &
        &   SUM(sb_pool_mt(ix_dom_assoc, ixC, ic, :) * soil_depth_sl(ic,:)), &
        &   SUM(sb_pool_mt(ix_soluable_litter, ixC, ic, :) * soil_depth_sl(ic,:)), &
        &   SUM(sb_pool_mt(ix_polymeric_litter, ixC, ic, :) * soil_depth_sl(ic,:)), &
        &   SUM(sb_pool_mt(ix_woody_litter, ixC, ic, :) * soil_depth_sl(ic,:)), &
        &   SUM(sb_pool_mt(ix_mycorrhiza, ixC, ic, :) * soil_depth_sl(ic,:)), &
        &   SUM(sb_pool_mt(ix_microbial, ixC, ic, :) * soil_depth_sl(ic,:)), &
        &   SUM(sb_pool_mt(ix_residue, ixC, ic, :) * soil_depth_sl(ic,:)), &
        &   SUM(sb_pool_mt(ix_residue_assoc, ixC, ic, :) * soil_depth_sl(ic,:)), &
        ! N
        &   SUM(sb_pool_mt(ix_dom, ixN, ic, :) * soil_depth_sl(ic,:)), &
        &   SUM(sb_pool_mt(ix_dom_assoc, ixN, ic, :) * soil_depth_sl(ic,:)), &
        &   SUM(sb_pool_mt(ix_soluable_litter, ixN, ic, :) * soil_depth_sl(ic,:)), &
        &   SUM(sb_pool_mt(ix_polymeric_litter, ixN, ic, :) * soil_depth_sl(ic,:)), &
        &   SUM(sb_pool_mt(ix_woody_litter, ixN, ic, :) * soil_depth_sl(ic,:)), &
        &   SUM(sb_pool_mt(ix_mycorrhiza, ixN, ic, :) * soil_depth_sl(ic,:)), &
        &   SUM(sb_pool_mt(ix_microbial, ixN, ic, :) * soil_depth_sl(ic,:)), &
        &   SUM(sb_pool_mt(ix_residue, ixN, ic, :) * soil_depth_sl(ic,:)), &
        &   SUM(sb_pool_mt(ix_residue_assoc, ixN, ic, :) * soil_depth_sl(ic,:)), &
        ! P
        &   SUM(sb_pool_mt(ix_dom, ixP, ic, :) * soil_depth_sl(ic,:)), &
        &   SUM(sb_pool_mt(ix_dom_assoc, ixP, ic, :) * soil_depth_sl(ic,:)), &
        &   SUM(sb_pool_mt(ix_soluable_litter, ixP, ic, :) * soil_depth_sl(ic,:)), &
        &   SUM(sb_pool_mt(ix_polymeric_litter, ixP, ic, :) * soil_depth_sl(ic,:)), &
        &   SUM(sb_pool_mt(ix_woody_litter, ixP, ic, :) * soil_depth_sl(ic,:)), &
        &   SUM(sb_pool_mt(ix_mycorrhiza, ixP, ic, :) * soil_depth_sl(ic,:)), &
        &   SUM(sb_pool_mt(ix_microbial, ixP, ic, :) * soil_depth_sl(ic,:)), &
        &   SUM(sb_pool_mt(ix_residue, ixP, ic, :) * soil_depth_sl(ic,:)), &
        &   SUM(sb_pool_mt(ix_residue_assoc, ixP, ic, :) * soil_depth_sl(ic,:))
    END IF
    IF (is_active_txt_soilpool_c_org_sl) THEN
      ! name_soilpool_c_org_sl
      write(lun_soilpool_c_org_sl,'(181F16.6)') &
        &   time_d,sb_pool_mt(ix_soluable_litter, ixC, ic, :),&
        &   sb_pool_mt(ix_polymeric_litter, ixC, ic, :),sb_pool_mt(ix_woody_litter, ixC, ic, :),&
        &   sb_pool_mt(ix_dom, ixC, ic, :),sb_pool_mt(ix_dom_assoc, ixC, ic, :),&
        &   sb_pool_mt(ix_mycorrhiza, ixC, ic, :),&
        &   sb_pool_mt(ix_microbial, ixC, ic, :),sb_pool_mt(ix_residue, ixC, ic, :),&
        &   sb_pool_mt(ix_residue_assoc, ixC, ic, :)
    END IF
    IF (is_active_txt_soilpool_n_org_sl) THEN
      ! name_soilpool_n_org_sl
      write(lun_soilpool_n_org_sl,'(181F16.6)') &
        &   time_d,sb_pool_mt(ix_soluable_litter, ixN, ic, :),&
        &   sb_pool_mt(ix_polymeric_litter, ixN, ic, :),sb_pool_mt(ix_woody_litter, ixN, ic, :),&
        &   sb_pool_mt(ix_dom, ixN, ic, :),sb_pool_mt(ix_dom_assoc, ixN, ic, :),&
        &   sb_pool_mt(ix_mycorrhiza, ixN, ic, :),&
        &   sb_pool_mt(ix_microbial, ixN, ic, :),sb_pool_mt(ix_residue, ixN, ic, :),&
        &   sb_pool_mt(ix_residue_assoc, ixN, ic, :)
    END IF
    IF (is_active_txt_soilpool_p_org_sl) THEN
      ! name_soilpool_p_org_sl
      write(lun_soilpool_p_org_sl,'(181F16.6)') &
        &   time_d,sb_pool_mt(ix_soluable_litter, ixP, ic, :),&
        &   sb_pool_mt(ix_polymeric_litter, ixP, ic, :),sb_pool_mt(ix_woody_litter, ixP, ic, :),&
        &   sb_pool_mt(ix_dom, ixP, ic, :),sb_pool_mt(ix_dom_assoc, ixP, ic, :),&
        &   sb_pool_mt(ix_mycorrhiza, ixP, ic, :),&
        &   sb_pool_mt(ix_microbial, ixP, ic, :),sb_pool_mt(ix_residue, ixP, ic, :),&
        &   sb_pool_mt(ix_residue_assoc, ixP, ic, :)
    END IF
    IF (is_active_txt_soilpool_ndn15_inorg_sl) THEN
      ! name_soilpool_ndn15_inorg_sl
      write(lun_soilpool_ndn15_inorg_sl,'(181F18.6)') time_d, &
        &   nh4_assoc(ic,:), nh4_solute(ic,:),no3_solute(ic,:)  , &
        &   noy(ic,:), n2o(ic,:), n2(ic,:), &
        &   calc_delta_N15(nh4_assoc(ic,:) ,nh4_n15_assoc(ic,:))  , &
        &   calc_delta_N15(nh4_solute(ic,:),nh4_n15_solute(ic,:)) , &
        &   calc_delta_N15(no3_solute(ic,:),no3_n15_solute(ic,:))
    END IF
    IF (is_active_txt_soilpool_p_inorg_sl) THEN
      ! name_soilpool_p_inorg_sl
      write(lun_soilpool_p_inorg_sl,'(101F16.6)') time_d, &
        &   po4_assoc_fast(ic,:), po4_assoc_slow(ic,:), po4_occluded(ic,:), po4_solute(ic,:), po4_primary(ic,:)
    END IF
    IF (is_active_txt_soilpool_dc13_org_sl) THEN
      ! name_soilpool_dc13_org_sl
      dC13soluable_litter(ic,:)   = calc_delta_C13(sb_pool_mt(ix_soluable_litter, ixC, ic, :), &
        &                           sb_pool_mt(ix_soluable_litter, ixC13, ic, :))
      dC13polymeric_litter(ic,:)  = calc_delta_C13(sb_pool_mt(ix_polymeric_litter, ixC, ic, :), &
        &                           sb_pool_mt(ix_polymeric_litter, ixC13, ic, :))
      dC13woody_litter(ic,:)      = calc_delta_C13(sb_pool_mt(ix_woody_litter, ixC, ic, :), &
        &                           sb_pool_mt(ix_woody_litter, ixC13, ic, :))
      dC13dom(ic,:)               = calc_delta_C13(sb_pool_mt(ix_dom, ixC, ic, :), &
        &                           sb_pool_mt(ix_dom, ixC13, ic, :))
      dC13dom_assoc(ic,:)         = calc_delta_C13(sb_pool_mt(ix_dom_assoc, ixC, ic, :), &
        &                           sb_pool_mt(ix_dom_assoc, ixC13, ic, :))
      dC13mycorrhiza(ic,:)        = calc_delta_C13(sb_pool_mt(ix_mycorrhiza, ixC, ic, :), &
        &                           sb_pool_mt(ix_mycorrhiza, ixC13, ic, :))
      dC13microbial(ic,:)         = calc_delta_C13(sb_pool_mt(ix_microbial, ixC, ic, :), &
        &                           sb_pool_mt(ix_microbial, ixC13, ic, :))
      dC13residue(ic,:)           = calc_delta_C13(sb_pool_mt(ix_residue, ixC, ic, :), &
        &                           sb_pool_mt(ix_residue, ixC13, ic, :))
      dC13residue_assoc(ic,:)     = calc_delta_C13(sb_pool_mt(ix_residue_assoc, ixC, ic, :), &
        &                           sb_pool_mt(ix_residue_assoc, ixC13, ic, :))
      write(lun_soilpool_dc13_org_sl,'(181F16.6)') time_d, &
        &   dC13soluable_litter(ic,:), dC13polymeric_litter(ic,:), &
        &   dC13woody_litter(ic,:), dC13dom(ic,:), &
        &   dC13dom_assoc(ic,:), &
        &   dC13mycorrhiza(ic,:), dC13microbial(ic,:), &
        &   dC13residue(ic,:), dC13residue_assoc(ic,:)
    END IF
    IF (is_active_txt_soilpool_dc14_org_sl) THEN
      ! name_soilpool_dc14_org_sl
      write(lun_soilpool_dc14_org_sl,'(181F16.6)') time_d, &
        &   calc_delta_C14(sb_pool_mt(ix_soluable_litter, ixC, ic, :), &
        &                  sb_pool_mt(ix_soluable_litter, ixC14, ic, :),dC13soluable_litter(ic,:)), &
        &   calc_delta_C14(sb_pool_mt(ix_polymeric_litter, ixC, ic, :), &
        &                  sb_pool_mt(ix_polymeric_litter, ixC14, ic, :),dC13polymeric_litter(ic,:)), &
        &   calc_delta_C14(sb_pool_mt(ix_woody_litter, ixC, ic, :), &
        &                  sb_pool_mt(ix_woody_litter, ixC14, ic, :),dC13woody_litter(ic,:)), &
        &   calc_delta_C14(sb_pool_mt(ix_dom, ixC, ic, :), &
        &                  sb_pool_mt(ix_dom, ixC14, ic, :),dC13dom(ic,:)), &
        &   calc_delta_C14(sb_pool_mt(ix_dom_assoc, ixC, ic, :), &
        &                  sb_pool_mt(ix_dom_assoc, ixC14, ic, :),dC13dom_assoc(ic,:)), &
        &   calc_delta_C14(sb_pool_mt(ix_mycorrhiza, ixC, ic, :), &
        &                  sb_pool_mt(ix_mycorrhiza, ixC14, ic, :),dC13mycorrhiza(ic,:)), &
        &   calc_delta_C14(sb_pool_mt(ix_microbial, ixC, ic, :), &
        &                  sb_pool_mt(ix_microbial, ixC14, ic, :),dC13microbial(ic,:)), &
        &   calc_delta_C14(sb_pool_mt(ix_residue, ixC, ic, :), &
        &                  sb_pool_mt(ix_residue, ixC14, ic, :),dC13residue(ic,:)), &
        &   calc_delta_C14(sb_pool_mt(ix_residue_assoc, ixC, ic, :), &
        &                  sb_pool_mt(ix_residue_assoc, ixC14, ic, :),dC13residue_assoc(ic,:))
    END IF
    IF (is_active_txt_soilpool_dn15_org_sl) THEN
      ! name_soilpool_dn15_org_sl
      write(lun_soilpool_dn15_org_sl,'(181F20.6)') time_d, &
        &   calc_delta_N15(sb_pool_mt(ix_soluable_litter, ixN, ic, :), &
        &                  sb_pool_mt(ix_soluable_litter, ixN15, ic, :)), &
        &   calc_delta_N15(sb_pool_mt(ix_polymeric_litter, ixN, ic, :), &
        &                  sb_pool_mt(ix_polymeric_litter, ixN15, ic, :)), &
        &   calc_delta_N15(sb_pool_mt(ix_woody_litter, ixN, ic, :), &
        &                  sb_pool_mt(ix_woody_litter, ixN15, ic, :)), &
        &   calc_delta_N15(sb_pool_mt(ix_dom, ixN, ic, :), &
        &                  sb_pool_mt(ix_dom, ixN15, ic, :)), &
        &   calc_delta_N15(sb_pool_mt(ix_dom_assoc, ixN, ic, :), &
        &                  sb_pool_mt(ix_dom_assoc, ixN15, ic, :)), &
        &   calc_delta_N15(sb_pool_mt(ix_mycorrhiza, ixN, ic, :), &
        &                  sb_pool_mt(ix_mycorrhiza, ixN15, ic, :)), &
        &   calc_delta_N15(sb_pool_mt(ix_microbial, ixN, ic, :), &
        &                  sb_pool_mt(ix_microbial, ixN15, ic, :)), &
        &   calc_delta_N15(sb_pool_mt(ix_residue, ixN, ic, :), &
        &                  sb_pool_mt(ix_residue, ixN15, ic, :)), &
        &   calc_delta_N15(sb_pool_mt(ix_residue_assoc, ixN, ic, :), &
        &                  sb_pool_mt(ix_residue_assoc, ixN15, ic, :))
    END IF
    IF (is_active_txt_soil_diagnostics_microbial_sl) THEN
      ! name_soil_diagnostics_microbial_sl
      hlp1_2d(:,:) = 0.0_wp
      DO is = 1, INT(num_sl_above_bedrock(ic))
        hlp1_2d(ic,is) = root_fraction_sl(ic,is) / soil_depth_sl(ic,is)
      END DO
      write(lun_soil_diagnostics_microbial_sl,'(281F16.6)') &
        &   time_d, microbial_cue_eff_tmic_mavg(ic,:), microbial_nue_eff_tmic_mavg(ic,:), &
        &   microbial_pue_eff_tmic_mavg(ic,:), &
        &   veg_pool_mt(ix_fine_root, ixC, ic) * hlp1_2d(ic,:), &
        &   enzyme_frac_poly(ic,:), enzyme_frac_residue(ic,:), enzyme_frac_poly_c(ic,:), &
        &   enzyme_frac_poly_n(ic,:), enzyme_frac_poly_p(ic,:), enzyme_frac_poly_c_mavg(ic,:), &
        &   enzyme_frac_poly_n_mavg(ic,:), enzyme_frac_poly_p_mavg(ic,:), &
        &   dom_cn_mavg(ic,:), dom_cp_mavg(ic,:)
    END IF
    IF (is_active_txt_soil_diagnostics_sorption_sl) THEN
      ! name_soil_diagnostics_sorption_sl
      write(lun_soil_diagnostics_sorption_sl,'(281F16.6)') &
        &   time_d, qmax_org(ic,:), bulk_dens_sl(ic,:), bulk_dens_corr_sl(ic,:), qmax_po4(ic,:), &
        &   km_adsorpt_po4_sl(ic,:), volume_min_sl(ic,:)
    END IF

    ! Deallocate local allocatable variables
    DEALLOCATE(dC13soluable_litter)
    DEALLOCATE(dC13polymeric_litter)
    DEALLOCATE(dC13woody_litter)
    DEALLOCATE(dC13dom)
    DEALLOCATE(dC13dom_assoc)
    DEALLOCATE(dC13mycorrhiza)
    DEALLOCATE(dC13microbial)
    DEALLOCATE(dC13residue)
    DEALLOCATE(dC13residue_assoc)
  END SUBROUTINE write_output_pool


  !-----------------------------------------------------------------------------------------------------!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! section II Fluxes                                                                                 !!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !-----------------------------------------------------------------------------------------------------!


  ! ======================================================================================================= !
  !>Fluxes: open files and write column headers
  !>
  SUBROUTINE open_flux_output_files(tile, options, &
    &                               output_interval_flux, spinup_id)

    USE mo_jsb_tile_class,         ONLY: t_jsb_tile_abstract
    USE mo_jsb_task_class,         ONLY: t_jsb_task_options
    ! ----------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_tile_abstract), INTENT(INOUT)   :: tile                 !< one tile with data structure for one lct
    TYPE(t_jsb_task_options),   INTENT(IN)      :: options              !< model%options
    CHARACTER(len=*),           INTENT(IN)      :: output_interval_flux !< output interval for "flux" variables
    CHARACTER(len=*),           INTENT(IN)      :: spinup_id            !< add ID to filename if file contains data from spinup simulation
    ! ----------------------------------------------------------------------------------------------------- !

    CHARACTER(100)              :: fname    !< filename helper
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':open_flux_output_files'
    ! ----------------------------------------------------------------------------------------------------- !

    !>1.0 define open() statements for active output files
    !>
    IF (is_active_txt_mycflux_cnp) THEN
      fname = name_mycflux_cnp//'_'//TRIM(spinup_id)//TRIM(output_interval_flux)//'.txt'
      OPEN (UNIT=lun_mycflux_cnp, FILE=fname, ACTION="write", STATUS="unknown")
    END IF
    IF (is_active_txt_myc_diagnostics) THEN
      fname = name_myc_diagnostics//'_'//TRIM(spinup_id)//TRIM(output_interval_flux)//'.txt'
      OPEN (UNIT=lun_myc_diagnostics, FILE=fname, ACTION="write", STATUS="unknown")
    END IF
    IF (is_active_txt_matter_conservation_test) THEN
      fname = name_matter_conservation_test//'_'//TRIM(spinup_id)//TRIM(output_interval_flux)//'.txt'
      OPEN (UNIT=lun_matter_conservation_test, FILE=fname, ACTION="write", STATUS="unknown")
    END IF
    IF (is_active_txt_vegflux_n) THEN
      fname = name_vegflux_n//'_'//TRIM(spinup_id)//TRIM(output_interval_flux)//'.txt'
      OPEN (UNIT=lun_vegflux_n, FILE=fname, ACTION="write", STATUS="unknown")
    END IF
    IF (is_active_txt_vegflux_p) THEN
      fname = name_vegflux_p//'_'//TRIM(spinup_id)//TRIM(output_interval_flux)//'.txt'
      OPEN (UNIT=lun_vegflux_p, FILE=fname, ACTION="write", STATUS="unknown")
    END IF
    IF (is_active_txt_vegflux_c) THEN
      fname = name_vegflux_c//'_'//TRIM(spinup_id)//TRIM(output_interval_flux)//'.txt'
      OPEN (UNIT=lun_vegflux_c, FILE=fname, ACTION="write", STATUS="unknown")
    END IF
    IF (is_active_txt_vegflux_dc13dc14) THEN
      fname = name_vegflux_dc13dc14//'_'//TRIM(spinup_id)//TRIM(output_interval_flux)//'.txt'
      OPEN (UNIT=lun_vegflux_dc13dc14, FILE=fname, ACTION="write", STATUS="unknown")
    END IF
    IF (is_active_txt_vegflux_h2o) THEN
      fname = name_vegflux_h2o//'_'//TRIM(spinup_id)//TRIM(output_interval_flux)//'.txt'
      OPEN (UNIT=lun_vegflux_h2o, FILE=fname, ACTION="write", STATUS="unknown")
    END IF
    IF (is_active_txt_vegflux_energy) THEN
      fname = name_vegflux_energy//'_'//TRIM(spinup_id)//TRIM(output_interval_flux)//'.txt'
      OPEN (UNIT=lun_vegflux_energy, FILE=fname, ACTION="write", STATUS="unknown")
    END IF
    IF (is_active_txt_atm_forcing) THEN
      fname = name_atm_forcing//'_'//TRIM(spinup_id)//TRIM(output_interval_flux)//'.txt'
      OPEN (UNIT=lun_atm_forcing, FILE=fname, ACTION="write", STATUS="unknown")
    END IF
    IF (is_active_txt_soilflux_cdc13dc14) THEN
      fname = name_soilflux_cdc13dc14//'_'//TRIM(spinup_id)//TRIM(output_interval_flux)//'.txt'
      OPEN (UNIT=lun_soilflux_cdc13dc14, FILE=fname, ACTION="write", STATUS="unknown")
    END IF
    IF (is_active_txt_soilflux_n) THEN
      fname = name_soilflux_n//'_'//TRIM(spinup_id)//TRIM(output_interval_flux)//'.txt'
      OPEN (UNIT=lun_soilflux_n, FILE=fname, ACTION="write", STATUS="unknown")
    END IF
    IF (is_active_txt_soilflux_p) THEN
      fname = name_soilflux_p//'_'//TRIM(spinup_id)//TRIM(output_interval_flux)//'.txt'
      OPEN (UNIT=lun_soilflux_p, FILE=fname, ACTION="write", STATUS="unknown")
    END IF
    IF (is_active_txt_soilflux_dn15) THEN
      fname = name_soilflux_dn15//'_'//TRIM(spinup_id)//TRIM(output_interval_flux)//'.txt'
      OPEN (UNIT=lun_soilflux_dn15, FILE=fname, ACTION="write", STATUS="unknown")
    END IF
    IF (is_active_txt_sb_mass_balance) THEN
      fname = name_sb_mass_balance//'_'//TRIM(spinup_id)//TRIM(output_interval_flux)//'.txt'
      OPEN (UNIT=lun_sb_mass_balance, FILE=fname, ACTION="write", STATUS="unknown")
    END IF
    IF (is_active_txt_soilflux_c_sl) THEN
      fname = name_soilflux_c_sl//'_'//TRIM(spinup_id)//TRIM(output_interval_flux)//'.txt'
      OPEN (UNIT=lun_soilflux_c_sl, FILE=fname, ACTION="write", STATUS="unknown")
    END IF
    IF (is_active_txt_soilflux_n_sl) THEN
      fname = name_soilflux_n_sl//'_'//TRIM(spinup_id)//TRIM(output_interval_flux)//'.txt'
      OPEN (UNIT=lun_soilflux_n_sl, FILE=fname, ACTION="write", STATUS="unknown")
    END IF
    IF (is_active_txt_soilflux_p_sl) THEN
      fname = name_soilflux_p_sl//'_'//TRIM(spinup_id)//TRIM(output_interval_flux)//'.txt'
      OPEN (UNIT=lun_soilflux_p_sl, FILE=fname, ACTION="write", STATUS="unknown")
    END IF

    !>2.0 write column headers into output files
    !>
    CALL init_output_writing_flux()

    !>3.0 set all flux-output *_acc variables to zero, but do not produce output
    !>  that is needed when, e.g., switching between diff. output periods (fluxnet diagnostics monthly -> timestep) 
    !>  to delete/set zero these *_acc variables that might contain data from previous output period & output file
    !>
    CALL write_output_flux(tile, &          ! tile
      &                    options, &       ! model&options
      &                    0.0_wp, &        ! time_d
      &                    0.0_wp, &        ! output_interval_flux 
      &                    .FALSE., &       ! write_flux_txt_output
      &                    .TRUE.)          ! set_zero_flux_acc_variables 
  END SUBROUTINE open_flux_output_files


  ! ======================================================================================================= !
  !>Fluxes: write column headers of output files
  !>
  SUBROUTINE init_output_writing_flux()

    USE mo_jsb_grid_class,        ONLY: t_jsb_vgrid
    USE mo_jsb_grid,              ONLY: Get_vgrid
    ! ----------------------------------------------------------------------------------------------------- !
    TYPE(t_jsb_vgrid),POINTER :: vgrid_soil_w                ! soil layer vertical grid
    INTEGER                   :: nsoil_w                     ! number of soil layers as used/defined by the SB_ process
    INTEGER                   :: isoil
    CHARACTER(2)              :: layerstr
    CHARACTER(450)            :: p1str,p2str,p3str,p4str,p5str,p6str,p7str,p8str,p9str,p10str, &
      &                          p11str,p12str,p13str,p14str,p15str,p16str,p17str,p18str,p19str,p20str 
    ! ----------------------------------------------------------------------------------------------------- !
    vgrid_soil_w => Get_vgrid('soil_depth_water')
    nsoil_w      =  vgrid_soil_w%n_levels
    ! ----------------------------------------------------------------------------------------------------- !

    !>1.0 write headers of active output files
    !>
    IF (is_active_txt_mycflux_cnp) THEN
      ! name_mycflux_cnp [mol / m2 / interval]
      WRITE(lun_mycflux_cnp,'(a)') 'Time mycExportC mycExportN mycExportP mycExportC13 mycExportC14 mycExportN15'
    END IF
    IF (is_active_txt_myc_diagnostics) THEN
      ! name_myc_diagnostics [...]
      WRITE(lun_myc_diagnostics,'(a)') 'Time fCm_resp fCm_SOM'
    END IF
    IF (is_active_txt_matter_conservation_test) THEN
      ! name_matter_conservation_test [mol m-2 timestep-1]
      WRITE(lun_matter_conservation_test,'(2a)') &
        & 'Time carbon_cons_test'
    END IF
    IF (is_active_txt_vegflux_n) THEN
      ! name_vegflux_n [mol N / m2 / interval]
      WRITE(lun_vegflux_n,'(7a)') &
        & 'Time Uptake_NH4 Uptake_NO3 RecyclingN ReserveUseN ', &
        & 'FixationN Fixation_dN15 Uptake_NH4_dN15 Uptake_NO3_dN15 GrowthN ', &
        & 'Growth_sapwood_dN15 ', &
        & 'GrowthN_leaf GrowthN_froot GrowthN_croot GrowthN_swood GrowthN_hwood ', &
        & 'GrowthN_labile GrowthN_reserve GrowthN_fruit GrowthN_sbed ', &
        & 'LitterfN_leaf LitterfN_froot LitterfN_croot LitterfN_swood LitterfN_hwood ', &
        & 'LitterfN_labile LitterfN_reserve LitterfN_fruit LitterfN_sbed'
    END IF
    IF (is_active_txt_vegflux_p) THEN
      ! name_vegflux_p [mol P / m2 / interval]
      WRITE(lun_vegflux_p,'(5a)') &
        & 'Time Uptake_PO4 RecyclingP ReserveUseP GrowthP ', &
        & 'GrowthP_leaf GrowthP_froot GrowthP_croot GrowthP_swood GrowthP_hwood ', &
        & 'GrowthP_labile GrowthP_reserve GrowthP_fruit GrowthP_sbed ', &
        & 'LitterfP_leaf LitterfP_froot LitterfP_croot LitterfP_swood LitterfP_hwood ', &
        & 'LitterfP_labile LitterfP_reserve LitterfP_fruit LitterfP_sbed'
    END IF
    IF (is_active_txt_vegflux_c) THEN
      ! name_vegflux_c [mol C / m2 / interval]
      WRITE(lun_vegflux_c,'(7a)') &
        & 'Time GPP GrowthResp MaintResp NTransformResp NFixationResp NPP ', &
        & 'leaf_herbivory fruit_herbivory ', &
        & 'ReserveUseC RootExudationC GrowthC ', &
        & 'GrowthC_leaf GrowthC_froot GrowthC_croot GrowthC_swood GrowthC_hwood ', &
        & 'GrowthC_labile GrowthC_reserve GrowthC_fruit GrowthC_sbed ', &
        & 'LitterfC_leaf LitterfC_froot LitterfC_croot LitterfC_swood LitterfC_hwood ', &
        & 'LitterfC_labile LitterfC_reserve LitterfC_fruit LitterfC_sbed'
    END IF
    IF (is_active_txt_vegflux_dc13dc14) THEN
      ! name_vegflux_dc13dc14 [...]
      WRITE(lun_vegflux_dc13dc14,'(2a)') 'Time GPP_dC13 GrowthResp_dC13 MaintResp_dC13 Growth_sapwood_dC13 ', &
        &   'GPP_DC14 GrowthResp_DC14 MaintResp_DC14 Growth_sapwood_DC14'
    END IF
    IF (is_active_txt_vegflux_h2o) THEN
      ! name_vegflux_h2o [kg/m2/interval, ... , m/interval, m/interval] (interval = output interval)
      WRITE(lun_vegflux_h2o,'(a)') &
        &   'Time Precip PotEvaporation Evaporation Interception Transpiration SrfRunoff Drainage GwRunoff GC GA'
    END IF
    IF (is_active_txt_vegflux_energy) THEN
      ! name_vegflux_energy [Wm-2, Wm-2, Wm-2, Wm-2, --, --, --, --]
      WRITE(lun_vegflux_energy,'(2a)') 'Time Rnet Qh Qle Qg Albedo_vis Albedo_nir Albedo_noon fPAR'
    END IF
    IF (is_active_txt_atm_forcing) THEN
      ! name_atm_forcing [...]
      WRITE(lun_atm_forcing,'(2a)') 'Time sum_Precip mean_AirTemperature mean_Wind mean_SWdown mean_LWdown mean_VPD ', &
        &   'atm_CO2 atm_CO2_dC13 atm_CO2_DC14'
    END IF
    IF (is_active_txt_soilflux_cdc13dc14) THEN
      ! name_soilflux_cdc13dc14 [mol C /m2 / interval]
      WRITE(lun_soilflux_cdc13dc14,'(4a)') 'Time Litterfall_C ', &
        &   'HetResp LeachingDOM_C LateralLossDOM_C Litterfall_dC13  ', &
        &   'HetResp_dC13 LeachingDOM_dC13 LateralLossDOM_dC13 ', &
        &   'Litterfall_DC14 HetResp_DC14 LeachingDOM_DC14 LateralLossDOM_DC14'
    END IF
    IF (is_active_txt_soilflux_n) THEN
      ! name_soilflux_n [mol N /m2 /interval]
      WRITE(lun_soilflux_n,'(5a)') 'Time Litterfall_N ', &
        &   'NDeposition Ndep_mavg Nfix_asymb NetMinNH4 NetMinNO3 Nitrification Denitrification Emission_NOY ', &
        &   'Emission_N2O Emission_N2 PlantUptakeNH4 PlantUptakeNO3 LeachingNH4 LeachingNO3 LeachingDOM_N ', &
        &   'LateralLoss_NH4 LateralLoss_NO3 LateralLossDOM_N ', &
        &   'Fast_exchange_nh4 MycUptake_NH4 MycUptake_NO3 MycUptake_NOrg' 
    END IF
    IF (is_active_txt_soilflux_p) THEN
      ! name_soilflux_p [mol P /m2 /interval]
      WRITE(lun_soilflux_p,'(4a)') 'Time Litterfall_P ', &
        &   'PDeposition Leaching_PO4 LeachingDOM_P LateralLoss_PO4 LateralLossDOM_P ', &
        &   'Net_Mineralisation_PO4 PlantUptake_PO4 Weathering_PO4 ', &
        &   'Occlusion_PO4 Slow_exchange_PO4 Fast_exchange_PO4 Biochem_mineralisation_PO4' 
    END IF
    IF (is_active_txt_soilflux_dn15) THEN
      ! name_soilflux_dn15 [...]
      WRITE(lun_soilflux_dn15,'(2a)') &
        &   'Time dN15_Deposition dN15_fix_asymb Emission_NOY_dN15 Emission_N2O_dN15 Emission_N2_dN15 ', &
        &   'Leaching_NH4_dN15 Leaching_NO3_dN15 LateralLoss_NH4_dN15 LateralLoss_NO3_dN15' 
    END IF
    IF (is_active_txt_sb_mass_balance) THEN
      ! name_sb_mass_balance [...]
      WRITE(lun_sb_mass_balance,'(2a)') &
        &   'Time sb_deltaPoolC sb_deltaPoolN sb_deltaPoolP sb_deltaPoolC13 sb_deltaPoolC14 sb_deltaPoolN15 ', &
        &   'sb_totalFluxC sb_totalFluxN sb_totalFluxP sb_totalFluxC13 sb_totalFluxC14 sb_totalFluxN15'
    END IF
    IF (is_active_txt_soilflux_c_sl) THEN
      ! name_soilflux_c_sl [mol C /m3 /interval]
      p1str=''
      p2str=''
      p3str=''
      p4str=''
      p5str=''
      p6str='' 
      DO isoil = 1, nsoil_w
        IF (isoil < 10) THEN 
          WRITE(layerstr,'(I1)') isoil
        ELSE
          WRITE(layerstr,'(I2)') isoil
        END IF
        WRITE(p1str,'(a)') TRIM(p1str)//' HetRespiration_L'//TRIM(layerstr)
        WRITE(p2str,'(a)') TRIM(p2str)//' HetRespirationC13_L'//TRIM(layerstr)
        WRITE(p3str,'(a)') TRIM(p3str)//' HetRespirationC14_L'//TRIM(layerstr)
        WRITE(p4str,'(a)') TRIM(p4str)//' Transport_DOC_L'//TRIM(layerstr)
        WRITE(p5str,'(a)') TRIM(p5str)//' Root_Respiration_L'//TRIM(layerstr)
        WRITE(p6str,'(a)') TRIM(p6str)//' LateralLoss_DOC_L'//TRIM(layerstr) 
      END DO
      WRITE(lun_soilflux_c_sl,'(7a)') 'Time ',TRIM(p1str),TRIM(p2str),TRIM(p3str),TRIM(p4str),TRIM(p5str),TRIM(p6str) 
    END IF
    IF (is_active_txt_soilflux_n_sl) THEN
      ! name_soilflux_n_sl [mol N /m3 /interval]
      p1str=''
      p2str=''
      p3str=''
      p4str=''
      p5str=''
      p6str=''
      p7str=''
      p8str=''
      p9str=''
      p10str=''
      p11str=''
      p12str=''
      p13str=''
      p14str=''
      p15str=''
      p16str=''
      p17str=''
      p18str='' 
      p19str='' 
      p20str='' 
      DO isoil = 1, nsoil_w
        IF (isoil < 10) THEN 
          WRITE(layerstr,'(I1)') isoil
        ELSE
          WRITE(layerstr,'(I2)') isoil
        END IF
        WRITE(p1str,'(a)') TRIM(p1str)//' NetMineral_nh4_L'//TRIM(layerstr)
        WRITE(p2str,'(a)') TRIM(p2str)//' PlantUptake_nh4_L'//TRIM(layerstr)
        WRITE(p3str,'(a)') TRIM(p3str)//' MicUptake_nh4_L'//TRIM(layerstr)
        WRITE(p4str,'(a)') TRIM(p4str)//' MycUptake_nh4_L'//TRIM(layerstr)
        WRITE(p5str,'(a)') TRIM(p5str)//' Transport_nh4_L'//TRIM(layerstr)
        WRITE(p6str,'(a)') TRIM(p6str)//' Transport_assoc_nh4_L'//TRIM(layerstr)
        WRITE(p7str,'(a)') TRIM(p7str)//' AsymbNfixation_L'//TRIM(layerstr)
        WRITE(p8str,'(a)') TRIM(p8str)//' NetMineral_no3_L'//TRIM(layerstr)
        WRITE(p9str,'(a)') TRIM(p9str)//' PlantUptake_no3_L'//TRIM(layerstr)
        WRITE(p10str,'(a)') TRIM(p10str)//' MicUptake_no3_L'//TRIM(layerstr)
        WRITE(p11str,'(a)') TRIM(p11str)//' MycUptake_no3_L'//TRIM(layerstr)
        WRITE(p12str,'(a)') TRIM(p12str)//' Transport_no3_L'//TRIM(layerstr)
        WRITE(p13str,'(a)') TRIM(p13str)//' Transport_DON_L'//TRIM(layerstr)
        WRITE(p14str,'(a)') TRIM(p14str)//' MycUptake_norg_L'//TRIM(layerstr)
        WRITE(p15str,'(a)') TRIM(p15str)//' Adsorpt_nh4_L'//TRIM(layerstr)
        WRITE(p16str,'(a)') TRIM(p16str)//' Nitrification_no3_L'//TRIM(layerstr)
        WRITE(p17str,'(a)') TRIM(p17str)//' Denitrification_n2_L'//TRIM(layerstr)
        WRITE(p18str,'(a)') TRIM(p18str)//' LateralLoss_nh4_L'//TRIM(layerstr) 
        WRITE(p19str,'(a)') TRIM(p19str)//' LateralLoss_no3_L'//TRIM(layerstr) 
        WRITE(p20str,'(a)') TRIM(p20str)//' LateralLoss_DON_L'//TRIM(layerstr) 
      END DO
      WRITE(lun_soilflux_n_sl,'(21a)') 'Time ', TRIM(p1str),TRIM(p2str),TRIM(p3str),TRIM(p4str),TRIM(p5str), &
        &   TRIM(p6str),TRIM(p7str),TRIM(p8str),TRIM(p9str),TRIM(p10str), &
        &   TRIM(p11str),TRIM(p12str),TRIM(p13str),TRIM(p14str),TRIM(p15str), &
        &   TRIM(p16str),TRIM(p17str), TRIM(p18str),TRIM(p19str),TRIM(p20str)
    END IF
    IF (is_active_txt_soilflux_p_sl) THEN
      ! name_soilflux_p_sl [mol P /m3 /interval]
      p1str=''
      p2str=''
      p3str=''
      p4str=''
      p6str=''
      p7str=''
      p8str=''
      p9str=''
      p10str=''
      p11str=''
      p12str=''
      p13str=''
      p14str=''
      p15str=''
      p16str=''
      p17str=''
      p18str='' 
      p19str='' 
      DO isoil = 1, nsoil_w
        IF (isoil < 10) THEN 
          WRITE(layerstr,'(I1)') isoil
        ELSE
          WRITE(layerstr,'(I2)') isoil
        END IF
        WRITE(p1str,'(a)') TRIM(p1str)//' Weath_L'//TRIM(layerstr)
        WRITE(p2str,'(a)') TRIM(p2str)//' Ocl_L'//TRIM(layerstr)
        WRITE(p3str,'(a)') TRIM(p3str)//' Exchange_fast_L'//TRIM(layerstr)
        WRITE(p4str,'(a)') TRIM(p4str)//' Exchange_slow_L'//TRIM(layerstr)
        WRITE(p6str,'(a)') TRIM(p6str)//' Adsorpt_fast_L'//TRIM(layerstr)
        WRITE(p7str,'(a)')  TRIM(p7str)//' NetMin_L'//TRIM(layerstr)
        WRITE(p8str,'(a)')  TRIM(p8str)//' Biochem_min_L'//TRIM(layerstr)
        WRITE(p9str,'(a)')  TRIM(p9str)//' Mic_upt_L'//TRIM(layerstr)
        WRITE(p10str,'(a)') TRIM(p10str)//' Plant_upt_L'//TRIM(layerstr)
        WRITE(p11str,'(a)') TRIM(p11str)//' Myc_upt_L'//TRIM(layerstr)
        WRITE(p12str,'(a)') TRIM(p12str)//' Transport_solute_L'//TRIM(layerstr)
        WRITE(p13str,'(a)') TRIM(p13str)//' Transport_assoc_fast_L'//TRIM(layerstr)
        WRITE(p14str,'(a)') TRIM(p14str)//' Transport_assoc_slow_L'//TRIM(layerstr)
        WRITE(p15str,'(a)') TRIM(p15str)//' Transport_occluded_L'//TRIM(layerstr)
        WRITE(p16str,'(a)') TRIM(p16str)//' Transport_primary_L'//TRIM(layerstr)
        WRITE(p17str,'(a)') TRIM(p17str)//' Transport_DOP_L'//TRIM(layerstr)
        WRITE(p18str,'(a)') TRIM(p18str)//' LateralLoss_solute_L'//TRIM(layerstr) 
        WRITE(p19str,'(a)') TRIM(p19str)//' LateralLoss_DOP_L'//TRIM(layerstr) 
      END DO
      WRITE(lun_soilflux_p_sl,'(19a)') 'Time ', TRIM(p1str),TRIM(p2str),TRIM(p3str),TRIM(p4str), &
        &   TRIM(p6str),TRIM(p7str),TRIM(p8str),TRIM(p9str),TRIM(p10str), &
        &   TRIM(p11str),TRIM(p12str),TRIM(p13str),TRIM(p14str),TRIM(p15str), &
        &   TRIM(p16str),TRIM(p17str),TRIM(p18str),TRIM(p19str) 
    END IF
  END SUBROUTINE init_output_writing_flux

  ! ======================================================================================================= !
  !>Fluxes: write data to files in a certain interval and accumulate fluxes
  !>
  SUBROUTINE write_output_flux(tile, options, time_d, output_interval_flux, &
    &                          write_flux_txt_output, &
    &                          set_zero_flux_acc_variables)

    USE mo_jsb_class,              ONLY: Get_model
    USE mo_jsb_model_class,        ONLY: t_jsb_model
    USE mo_jsb_tile_class,         ONLY: t_jsb_tile_abstract
    USE mo_jsb_task_class,         ONLY: t_jsb_task_options
    USE mo_jsb_process_class,      ONLY: SB_, A2L_, L2A_, Q_ASSIMI_, VEG_, RAD_, SPQ_, HYDRO_, SEB_, SSE_
    USE mo_jsb_grid_class,         ONLY: t_jsb_vgrid
    USE mo_jsb_grid,               ONLY: Get_vgrid
    USE mo_jsb_math_constants,     ONLY: eps8
    USE mo_jsb_physical_constants, ONLY: molar_mass_C, molar_mass_N
    USE mo_sb_constants,           ONLY: organic_cn
    USE mo_isotope_util,           ONLY: calc_delta_C13, calc_delta_C14, calc_delta_N15
    USE mo_spq_config_class,       ONLY: max_soil_layers
    USE mo_phy_schemes,            ONLY: qsat_mixed
    USE mo_atmland_constants,      ONLY: eps_vpd
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Use_config(SB_)
    dsl4jsb_Use_memory(A2L_)
    dsl4jsb_Use_memory(L2A_)
    dsl4jsb_Use_memory(Q_ASSIMI_)
    dsl4jsb_Use_memory(VEG_)
    dsl4jsb_Use_memory(RAD_)
    dsl4jsb_Use_memory(SPQ_)
    dsl4jsb_Use_memory(HYDRO_)
    dsl4jsb_Use_memory(SEB_)
    dsl4jsb_Use_memory(SSE_)
    dsl4jsb_Use_memory(SB_)
    ! ----------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile                     !< one tile with data structure for one lct
    TYPE(t_jsb_task_options),   INTENT(in)    :: options                  !< model options
    REAL(wp),                   INTENT(in)    :: time_d                   !< time in days
    REAL(wp),                   INTENT(in)    :: output_interval_flux     !< output interval in days (1.0 == one day exactly)
    LOGICAL,                    INTENT(in)    :: write_flux_txt_output        !< if TRUE: write output at current timestep TRUE/FALSE
    LOGICAL,                    INTENT(in)    :: set_zero_flux_acc_variables  !< if TRUE: run this routine but only set all *_acc variables to zero 
    ! ----------------------------------------------------------------------------------------------------- !
    TYPE(t_jsb_model),      POINTER   :: model                        !< the model
    TYPE(t_lnd_bgcm_store), POINTER   :: bgcm_store                   !< the bgcm store of this tile
    TYPE(t_lctlib_element), POINTER   :: lctlib                       !< land-cover-type library - parameter across pft's
    TYPE(t_jsb_vgrid), POINTER        :: vgrid_soil_w                 !< Vertical grid
    INTEGER                           :: nsoil_w                      !< number of soil layers as used/defined by the SB_ process
    LOGICAL                           :: set_zero_flux_acc_var_local  !< local variable of input variable: set_zero_flux_acc_variables
    ! local accumulation variables
      ! Note: the *_acc variables are initialized only the first time the subroutine is called. 
      ! On subsequent calls the old value of *_acc is retained.
      ! To improve code readability the SAVE attribute is added, even though this is redundant.
    REAL(wp), SAVE                        :: timestep_counter_flux_acc    = 0.0_wp
    REAL(wp), SAVE                        :: timestep_counter_albedo_acc  = 0.0_wp
    REAL(wp), SAVE                        :: gpp_acc             = 0.0_wp
    REAL(wp), SAVE                        :: growth_resp_acc     = 0.0_wp
    REAL(wp), SAVE                        :: maint_resp_acc      = 0.0_wp   
    REAL(wp), SAVE                        :: n_transform_resp_acc = 0.0_wp   
    REAL(wp), SAVE                        :: n_fixation_resp_acc = 0.0_wp 
    REAL(wp), SAVE                        :: exudation_c_acc     = 0.0_wp  
    REAL(wp), SAVE                        :: npp_acc             = 0.0_wp
    REAL(wp), SAVE                        :: leaf_herb_resp_acc  = 0.0_wp
    REAL(wp), SAVE                        :: fruit_herb_resp_acc = 0.0_wp
    REAL(wp), SAVE                        :: res_use_c_acc       = 0.0_wp   
    REAL(wp), SAVE                        :: gpp_C13_acc         = 0.0_wp
    REAL(wp), SAVE                        :: growth_resp_C13_acc = 0.0_wp
    REAL(wp), SAVE                        :: maint_resp_C13_acc  = 0.0_wp
    REAL(wp), SAVE                        :: gpp_C14_acc         = 0.0_wp
    REAL(wp), SAVE                        :: growth_resp_C14_acc = 0.0_wp
    REAL(wp), SAVE                        :: maint_resp_C14_acc  = 0.0_wp
    REAL(wp), SAVE                        :: uptake_nh4_acc      = 0.0_wp
    REAL(wp), SAVE                        :: plant_uptake_nh4_acc = 0.0_wp   
    REAL(wp), SAVE                        :: uptake_nh4_n15_acc  = 0.0_wp   
    REAL(wp), SAVE                        :: uptake_no3_acc      = 0.0_wp 
    REAL(wp), SAVE                        :: plant_uptake_no3_acc = 0.0_wp  
    REAL(wp), SAVE                        :: uptake_no3_n15_acc  = 0.0_wp   
    REAL(wp), SAVE                        :: recycling_n_acc     = 0.0_wp  
    REAL(wp), SAVE                        :: res_use_n_acc       = 0.0_wp  
    REAL(wp), SAVE                        :: n_fixation_acc      = 0.0_wp   
    REAL(wp), SAVE                        :: n_fixation_n15_acc  = 0.0_wp   
    REAL(wp), SAVE                        :: plant_uptake_po4_acc = 0.0_wp
    REAL(wp), SAVE                        :: uptake_po4_acc      = 0.0_wp
    REAL(wp), SAVE                        :: recycling_p_acc     = 0.0_wp  
    REAL(wp), SAVE                        :: res_use_p_acc       = 0.0_wp   
    REAL(wp), SAVE                        :: precip_acc          = 0.0_wp 
    REAL(wp), SAVE                        :: evaporation_acc     = 0.0_wp 
    REAL(wp), SAVE                        :: transpiration_acc   = 0.0_wp 
    REAL(wp), SAVE                        :: interception_acc    = 0.0_wp 
    REAL(wp), SAVE                        :: evapopot_acc        = 0.0_wp 
    REAL(wp), SAVE                        :: srf_runoff_acc      = 0.0_wp 
    REAL(wp), SAVE                        :: drainage_acc        = 0.0_wp 
    REAL(wp), SAVE                        :: gw_runoff_acc       = 0.0_wp
    REAL(wp), SAVE                        :: gc_acc              = 0.0_wp
    REAL(wp), SAVE                        :: ga_acc              = 0.0_wp
    REAL(wp), SAVE                        :: n_deposition_acc    = 0.0_wp
    REAL(wp), SAVE                        :: n_deposition_n15_acc = 0.0_wp
    REAL(wp), SAVE                        :: rnet_acc            = 0.0_wp
    REAL(wp), SAVE                        :: qh_acc              = 0.0_wp
    REAL(wp), SAVE                        :: qle_acc             = 0.0_wp
    REAL(wp), SAVE                        :: qg_acc              = 0.0_wp
    REAL(wp), SAVE                        :: alb_vis_acc         = 0.0_wp
    REAL(wp), SAVE                        :: alb_nir_acc         = 0.0_wp
    ! veg_growth veg_litterfall fluxes CNP
    REAL(wp), SAVE                        :: growth_c_acc        = 0.0_wp   ! sum of all tissues, but seed_bed
    REAL(wp), SAVE                        :: growth_n_acc        = 0.0_wp   ! sum of all tissues, but seed_bed
    REAL(wp), SAVE                        :: growth_p_acc        = 0.0_wp   ! sum of all tissues, but seed_bed
    REAL(wp), SAVE    :: veg_growth_leaf_c_acc                  = 0.0_wp
    REAL(wp), SAVE    :: veg_growth_fine_root_c_acc             = 0.0_wp
    REAL(wp), SAVE    :: veg_growth_coarse_root_c_acc           = 0.0_wp
    REAL(wp), SAVE    :: veg_growth_sap_wood_c_acc              = 0.0_wp
    REAL(wp), SAVE    :: veg_growth_heart_wood_c_acc            = 0.0_wp
    REAL(wp), SAVE    :: veg_growth_labile_c_acc                = 0.0_wp
    REAL(wp), SAVE    :: veg_growth_reserve_c_acc               = 0.0_wp
    REAL(wp), SAVE    :: veg_growth_fruit_c_acc                 = 0.0_wp
    REAL(wp), SAVE    :: veg_growth_seed_bed_c_acc              = 0.0_wp
    REAL(wp), SAVE    :: veg_growth_leaf_n_acc                  = 0.0_wp
    REAL(wp), SAVE    :: veg_growth_fine_root_n_acc             = 0.0_wp
    REAL(wp), SAVE    :: veg_growth_coarse_root_n_acc           = 0.0_wp
    REAL(wp), SAVE    :: veg_growth_sap_wood_n_acc              = 0.0_wp
    REAL(wp), SAVE    :: veg_growth_heart_wood_n_acc            = 0.0_wp
    REAL(wp), SAVE    :: veg_growth_labile_n_acc                = 0.0_wp
    REAL(wp), SAVE    :: veg_growth_reserve_n_acc               = 0.0_wp
    REAL(wp), SAVE    :: veg_growth_fruit_n_acc                 = 0.0_wp
    REAL(wp), SAVE    :: veg_growth_seed_bed_n_acc              = 0.0_wp
    REAL(wp), SAVE    :: veg_growth_leaf_p_acc                  = 0.0_wp
    REAL(wp), SAVE    :: veg_growth_fine_root_p_acc             = 0.0_wp
    REAL(wp), SAVE    :: veg_growth_coarse_root_p_acc           = 0.0_wp
    REAL(wp), SAVE    :: veg_growth_sap_wood_p_acc              = 0.0_wp
    REAL(wp), SAVE    :: veg_growth_heart_wood_p_acc            = 0.0_wp
    REAL(wp), SAVE    :: veg_growth_labile_p_acc                = 0.0_wp
    REAL(wp), SAVE    :: veg_growth_reserve_p_acc               = 0.0_wp
    REAL(wp), SAVE    :: veg_growth_fruit_p_acc                 = 0.0_wp
    REAL(wp), SAVE    :: veg_growth_seed_bed_p_acc              = 0.0_wp
    REAL(wp), SAVE    :: veg_litterfall_leaf_c_acc              = 0.0_wp
    REAL(wp), SAVE    :: veg_litterfall_fine_root_c_acc         = 0.0_wp
    REAL(wp), SAVE    :: veg_litterfall_coarse_root_c_acc       = 0.0_wp
    REAL(wp), SAVE    :: veg_litterfall_sap_wood_c_acc          = 0.0_wp
    REAL(wp), SAVE    :: veg_litterfall_heart_wood_c_acc        = 0.0_wp
    REAL(wp), SAVE    :: veg_litterfall_labile_c_acc            = 0.0_wp
    REAL(wp), SAVE    :: veg_litterfall_reserve_c_acc           = 0.0_wp
    REAL(wp), SAVE    :: veg_litterfall_fruit_c_acc             = 0.0_wp
    REAL(wp), SAVE    :: veg_litterfall_seed_bed_c_acc          = 0.0_wp
    REAL(wp), SAVE    :: veg_litterfall_leaf_n_acc              = 0.0_wp
    REAL(wp), SAVE    :: veg_litterfall_fine_root_n_acc         = 0.0_wp
    REAL(wp), SAVE    :: veg_litterfall_coarse_root_n_acc       = 0.0_wp
    REAL(wp), SAVE    :: veg_litterfall_sap_wood_n_acc          = 0.0_wp
    REAL(wp), SAVE    :: veg_litterfall_heart_wood_n_acc        = 0.0_wp
    REAL(wp), SAVE    :: veg_litterfall_labile_n_acc            = 0.0_wp
    REAL(wp), SAVE    :: veg_litterfall_reserve_n_acc           = 0.0_wp
    REAL(wp), SAVE    :: veg_litterfall_fruit_n_acc             = 0.0_wp
    REAL(wp), SAVE    :: veg_litterfall_seed_bed_n_acc          = 0.0_wp
    REAL(wp), SAVE    :: veg_litterfall_leaf_p_acc              = 0.0_wp
    REAL(wp), SAVE    :: veg_litterfall_fine_root_p_acc         = 0.0_wp
    REAL(wp), SAVE    :: veg_litterfall_coarse_root_p_acc       = 0.0_wp
    REAL(wp), SAVE    :: veg_litterfall_sap_wood_p_acc          = 0.0_wp
    REAL(wp), SAVE    :: veg_litterfall_heart_wood_p_acc        = 0.0_wp
    REAL(wp), SAVE    :: veg_litterfall_labile_p_acc            = 0.0_wp
    REAL(wp), SAVE    :: veg_litterfall_reserve_p_acc           = 0.0_wp
    REAL(wp), SAVE    :: veg_litterfall_fruit_p_acc             = 0.0_wp
    REAL(wp), SAVE    :: veg_litterfall_seed_bed_p_acc          = 0.0_wp
    REAL(wp), SAVE                        :: growth_swood_c13_acc = 0.0_wp
    REAL(wp), SAVE                        :: growth_swood_c14_acc = 0.0_wp
    REAL(wp), SAVE                        :: growth_swood_n15_acc = 0.0_wp
    ! soil - local accumulation variables (partly also used for land output)
    REAL(wp), SAVE                        :: litterfallC_acc     = 0.0_wp
    REAL(wp), SAVE                        :: litterfallC13_acc   = 0.0_wp
    REAL(wp), SAVE                        :: litterfallC14_acc   = 0.0_wp
    REAL(wp), SAVE                        :: litterfallN_acc     = 0.0_wp
    REAL(wp), SAVE                        :: litterfallP_acc     = 0.0_wp
    REAL(wp), SAVE                        :: het_respC_acc       = 0.0_wp
    REAL(wp), SAVE                        :: het_respC13_acc     = 0.0_wp
    REAL(wp), SAVE                        :: het_respC14_acc     = 0.0_wp
    REAL(wp), SAVE                        :: leachingCorg_acc    = 0.0_wp
    REAL(wp), SAVE                        :: leachingC13org_acc  = 0.0_wp
    REAL(wp), SAVE                        :: leachingC14org_acc  = 0.0_wp
    REAL(wp), SAVE                        :: leachingNorg_acc    = 0.0_wp
    REAL(wp), SAVE                        :: leachingPorg_acc    = 0.0_wp
    REAL(wp), SAVE                        :: lateral_loss_Corg_acc    = 0.0_wp 
    REAL(wp), SAVE                        :: lateral_loss_C13org_acc  = 0.0_wp 
    REAL(wp), SAVE                        :: lateral_loss_C14org_acc  = 0.0_wp 
    REAL(wp), SAVE                        :: lateral_loss_Norg_acc    = 0.0_wp 
    REAL(wp), SAVE                        :: lateral_loss_Porg_acc    = 0.0_wp 
    REAL(wp), SAVE                        :: net_min_nh4_acc     = 0.0_wp
    REAL(wp), SAVE                        :: net_min_no3_acc     = 0.0_wp
    REAL(wp), SAVE                        :: net_min_po4_acc     = 0.0_wp
    REAL(wp), SAVE                        :: nitrification_acc   = 0.0_wp
    REAL(wp), SAVE                        :: denitrification_acc = 0.0_wp
    REAL(wp), SAVE                        :: fast_exchange_nh4_acc  = 0.0_wp
    REAL(wp), SAVE                        :: emission_noy_n15_acc = 0.0_wp
    REAL(wp), SAVE                        :: emission_noy_acc    = 0.0_wp
    REAL(wp), SAVE                        :: emission_n2o_acc    = 0.0_wp
    REAL(wp), SAVE                        :: emission_n2o_n15_acc = 0.0_wp
    REAL(wp), SAVE                        :: emission_n2_acc     = 0.0_wp
    REAL(wp), SAVE                        :: emission_n2_n15_acc = 0.0_wp
    REAL(wp), SAVE                        :: leaching_nh4_acc    = 0.0_wp
    REAL(wp), SAVE                        :: leaching_nh4_n15_acc = 0.0_wp
    REAL(wp), SAVE                        :: leaching_no3_acc    = 0.0_wp
    REAL(wp), SAVE                        :: leaching_no3_n15_acc = 0.0_wp
    REAL(wp), SAVE                        :: p_deposition_acc = 0.0_wp
    REAL(wp), SAVE                        :: leaching_po4_acc    = 0.0_wp
    REAL(wp), SAVE                        :: lateral_loss_nh4_acc  = 0.0_wp 
    REAL(wp), SAVE                        :: lateral_loss_nh4_n15_acc = 0.0_wp
    REAL(wp), SAVE                        :: lateral_loss_no3_acc    = 0.0_wp
    REAL(wp), SAVE                        :: lateral_loss_no3_n15_acc = 0.0_wp 
    REAL(wp), SAVE                        :: lateral_loss_po4_acc    = 0.0_wp  
    REAL(wp), SAVE                        :: weathering_po4_acc    = 0.0_wp
    REAL(wp), SAVE                        :: occlusion_po4_acc    = 0.0_wp
    REAL(wp), SAVE                        :: slow_exchange_po4_acc    = 0.0_wp
    REAL(wp), SAVE                        :: fast_exchange_po4_acc  = 0.0_wp
    REAL(wp), SAVE                        :: biochem_mineralisation_po4_acc  = 0.0_wp
    REAL(wp), SAVE                        :: asymb_n_fixation_acc = 0.0_wp
    REAL(wp), SAVE                        :: asymb_n_fixation_n15_acc = 0.0_wp
    REAL(wp), SAVE                        :: sb_delta_pools_carbon_acc     = 0.0_wp
    REAL(wp), SAVE                        :: sb_delta_pools_nitrogen_acc   = 0.0_wp
    REAL(wp), SAVE                        :: sb_delta_pools_phosphorus_acc = 0.0_wp
    REAL(wp), SAVE                        :: sb_delta_pools_carbon13_acc   = 0.0_wp
    REAL(wp), SAVE                        :: sb_delta_pools_carbon14_acc   = 0.0_wp
    REAL(wp), SAVE                        :: sb_delta_pools_nitrogen15_acc = 0.0_wp
    REAL(wp), SAVE                        :: sb_total_flux_carbon_acc      = 0.0_wp
    REAL(wp), SAVE                        :: sb_total_flux_nitrogen_acc    = 0.0_wp
    REAL(wp), SAVE                        :: sb_total_flux_phosphorus_acc  = 0.0_wp
    REAL(wp), SAVE                        :: sb_total_flux_carbon13_acc    = 0.0_wp
    REAL(wp), SAVE                        :: sb_total_flux_carbon14_acc    = 0.0_wp
    REAL(wp), SAVE                        :: sb_total_flux_nitrogen15_acc  = 0.0_wp
    REAL(wp), SAVE                        :: myc_export_c_acc              = 0.0_wp
    REAL(wp), SAVE                        :: myc_export_n_acc              = 0.0_wp
    REAL(wp), SAVE                        :: myc_export_p_acc              = 0.0_wp
    REAL(wp), SAVE                        :: myc_export_c13_acc            = 0.0_wp
    REAL(wp), SAVE                        :: myc_export_c14_acc            = 0.0_wp
    REAL(wp), SAVE                        :: myc_export_n15_acc            = 0.0_wp
    REAL(wp), SAVE                        :: myc_uptake_nh4_acc            = 0.0_wp
    REAL(wp), SAVE                        :: myc_uptake_no3_acc            = 0.0_wp
    REAL(wp), SAVE                        :: myc_uptake_norg_acc           = 0.0_wp
    REAL(wp), SAVE                        :: myc_respiration_c_acc         = 0.0_wp
    REAL(wp), SAVE                        :: f_cm_som_acc                  = 0.0_wp
    ! atmospheric forcing   (note: 'precip_acc' is already calculated for name_vegflux_h2o output)
    REAL(wp), SAVE                        :: t_air_acc                     = 0.0_wp
    REAL(wp), SAVE                        :: wind_air_acc                  = 0.0_wp
    REAL(wp), SAVE                        :: swvis_srf_down_acc            = 0.0_wp
    REAL(wp), SAVE                        :: swnir_srf_down_acc            = 0.0_wp
    REAL(wp), SAVE                        :: lw_srf_down_acc               = 0.0_wp
    REAL(wp), SAVE                        :: vpd_acc                       = 0.0_wp
    REAL(wp), SAVE                        :: co2_mixing_ratio_acc          = 0.0_wp
    REAL(wp), SAVE                        :: co2_mixing_ratio_C13_acc      = 0.0_wp
    REAL(wp), SAVE                        :: co2_mixing_ratio_C14_acc      = 0.0_wp
    ! variables with accumulated information per soil layer
    ! hardcoded Dimension, using max soil layer number of 20 (which actually is the max soil layer number)
    REAL(wp), DIMENSION(max_soil_layers), SAVE :: &
      ! C
      &   het_respiration_lay_acc = 0.0_wp       , &
      &   het_respiration_c13_lay_acc = 0.0_wp   , &
      &   het_respiration_c14_lay_acc = 0.0_wp   , &
      &   transport_DOC_lay_acc = 0.0_wp         , &
      &   root_respiration_lay_acc = 0.0_wp      , &
      &   lateral_loss_Corg_lay_acc = 0.0_wp     , &
      ! N
      &   net_mineral_nh4_lay_acc = 0.0_wp       , &
      &   plant_uptake_nh4_lay_acc = 0.0_wp      , &
      &   microbial_upt_nh4_lay_acc = 0.0_wp     , &
      &   microbial_upt_no3_lay_acc = 0.0_wp     , &
      &   myc_uptake_nh4_lay_acc = 0.0_wp        , &
      &   asymb_n_fixation_lay_acc = 0.0_wp      , &
      &   net_mineral_no3_lay_acc = 0.0_wp       , &
      &   plant_uptake_no3_lay_acc = 0.0_wp      , &
      &   myc_uptake_no3_lay_acc = 0.0_wp        , &
      &   myc_uptake_norg_lay_acc = 0.0_wp       , &
      &   transport_solute_nh4_acc = 0.0_wp      , &
      &   transport_nh4_assoc_acc = 0.0_wp       , &
      &   transport_solute_no3_acc = 0.0_wp      , &
      &   transport_DON_lay_acc = 0.0_wp         , &
      &   f_adsorpt_nh4_lay_acc = 0.0_wp         , &
      &   nitrification_lay_acc = 0.0_wp         , &
      &   denitrification_lay_acc = 0.0_wp       , &
      &   lateral_loss_nh4_lay_acc = 0.0_wp      , &
      &   lateral_loss_no3_lay_acc = 0.0_wp      , &
      &   lateral_loss_Norg_lay_acc = 0.0_wp     , &
      ! P
      &   f_exchange_po4_acc = 0.0_wp        , &
      &   f_adsorpt_po4_acc = 0.0_wp         , &
      &   net_min_po4_lay_acc = 0.0_wp       , &
      &   ocl_po4_acc = 0.0_wp               , &
      &   weath_po4_acc = 0.0_wp             , &
      &   biochem_min_po4_acc = 0.0_wp       , &
      &   s_exchange_po4_acc = 0.0_wp        , &
      &   plant_upt_po4_acc = 0.0_wp         , &
      &   microbial_upt_po4_acc = 0.0_wp     , &
      &   myc_upt_po4_acc = 0.0_wp           , &
      &   transport_solute_po4_acc = 0.0_wp      , &
      &   transport_assoc_fast_po4_acc = 0.0_wp  , &
      &   transport_assoc_slow_po4_acc = 0.0_wp  , &
      &   transport_occluded_po4_acc = 0.0_wp    , &
      &   transport_primary_po4_acc = 0.0_wp     , &
      &   transport_DOP_lay_acc = 0.0_wp         , &
      &   lateral_loss_po4_lay_acc = 0.0_wp      , &
      &   lateral_loss_Porg_lay_acc = 0.0_wp     , &
      &   c_cons_test_acc = 0.0_wp
    REAL(wp), ALLOCATABLE, DIMENSION(:,:) :: hlp1_2d
    REAL(wp)                              :: hlp1
    REAL(wp)                              :: dC13gpp, dC13rg, dC13rm, dC13wgr
    REAL(wp)                              :: dC13litterfall, dC13het_resp, dC13leaching, dC13lateral_loss 
    REAL(wp)                              :: dC13atm
    REAL(wp)                              :: mean_sw_down
    REAL(wp)                              :: dtime_one
    REAL(wp)                              :: dtime                !< timestep lengt
    REAL(wp)                              :: f_veg_flux_output
    INTEGER                               :: iblk, ics, ice, nc, ic, is
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':write_output_flux'
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Def_mt2L2D :: veg_pool_mt
    dsl4jsb_Def_mt2L2D :: veg_litterfall_mt
    dsl4jsb_Def_mt1L2D :: seed_bed_growth_mt
    dsl4jsb_Def_mt1L2D :: seed_bed_litterfall_mt
    dsl4jsb_Def_mt2L2D :: veg_growth_mt
    dsl4jsb_Def_mt1L2D :: veg_exudation_mt
    dsl4jsb_Def_mt1L2D :: veg_reserve_use_mt
    dsl4jsb_Def_mt2L3D :: sb_transport_mt
    dsl4jsb_Def_mt1L3D :: sb_mycorrhiza_export_mt
    dsl4jsb_Def_mt1L3D :: sb_delta_pools_mt
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Def_config(SB_)
    dsl4jsb_Def_memory(A2L_)
    dsl4jsb_Def_memory(L2A_)
    dsl4jsb_Def_memory(Q_ASSIMI_)
    dsl4jsb_Def_memory(VEG_)
    dsl4jsb_Def_memory(RAD_)
    dsl4jsb_Def_memory(SPQ_)
    dsl4jsb_Def_memory(HYDRO_)
    dsl4jsb_Def_memory(SEB_)
    dsl4jsb_Def_memory(SSE_)
    dsl4jsb_Def_memory(SB_)
    ! ----------------------------------------------------------------------------------------------------- !
    ! A2L_
    dsl4jsb_Real2D_onChunk      :: swpar_srf_down
    dsl4jsb_Real2D_onChunk      :: nhx_deposition
    dsl4jsb_Real2D_onChunk      :: noy_deposition
    dsl4jsb_Real2D_onChunk      :: nhx_n15_deposition
    dsl4jsb_Real2D_onChunk      :: noy_n15_deposition
    dsl4jsb_Real2D_onChunk      :: p_deposition
    dsl4jsb_Real2D_onChunk      :: ndep_mavg
    dsl4jsb_Real2D_onChunk      :: rain 
    dsl4jsb_Real2D_onChunk      :: snow
    dsl4jsb_Real2D_onChunk      :: t_air
    dsl4jsb_Real2D_onChunk      :: wind_air
    dsl4jsb_Real2D_onChunk      :: swvis_srf_down
    dsl4jsb_Real2D_onChunk      :: swnir_srf_down
    dsl4jsb_Real2D_onChunk      :: lw_srf_down
    dsl4jsb_Real2D_onChunk      :: co2_mixing_ratio
    dsl4jsb_Real2D_onChunk      :: co2_mixing_ratio_C13
    dsl4jsb_Real2D_onChunk      :: co2_mixing_ratio_C14
    dsl4jsb_Real2D_onChunk      :: press_srf
    dsl4jsb_Real2D_onChunk      :: q_air
    ! L2A_
    dsl4jsb_Real2D_onChunk      :: q_c_conservation_test
    ! Q_ASSIMI_
    dsl4jsb_Real2D_onChunk      :: gross_assimilation
    dsl4jsb_Real2D_onChunk      :: gross_assimilation_C13
    dsl4jsb_Real2D_onChunk      :: gross_assimilation_C14
    dsl4jsb_Real2D_onChunk      :: canopy_cond
    dsl4jsb_Real2D_onChunk      :: aerodyn_cond
    ! HYDRO_ 2D
    dsl4jsb_Real2D_onChunk      :: num_sl_above_bedrock
    dsl4jsb_Real2D_onChunk      :: transpiration
    dsl4jsb_Real2D_onChunk      :: evapo_skin
    dsl4jsb_Real2D_onChunk      :: runoff
    ! HYDRO_ 3D
    dsl4jsb_Real3D_onChunk      :: soil_depth_sl
    ! HYDRO_ 2D / SPQ_ 2D
    dsl4jsb_Real2D_onChunk      :: evapopot
    dsl4jsb_Real2D_onChunk      :: drainage
    ! SPQ_ 2D
    dsl4jsb_Real2D_onChunk      :: evaporation
    dsl4jsb_Real2D_onChunk      :: interception
    dsl4jsb_Real2D_onChunk      :: srf_runoff 
    dsl4jsb_Real2D_onChunk      :: gw_runoff
    dsl4jsb_Real2D_onChunk      :: sensible_heat_flx
    dsl4jsb_Real2D_onChunk      :: latent_heat_flx
    dsl4jsb_Real2D_onChunk      :: ground_heat_flx
    ! SEB_ 2D
    dsl4jsb_Real2D_onChunk      :: sensible_hflx
    dsl4jsb_Real2D_onChunk      :: latent_hflx
    ! SSE_ 2D
    dsl4jsb_Real2D_onChunk      :: grnd_hflx
    ! RAD_
    dsl4jsb_Real2D_onChunk      :: albedo_noon
    dsl4jsb_Real2D_onChunk      :: rad_srf_net
    dsl4jsb_Real2D_onChunk      :: fapar_noon
    dsl4jsb_Real2D_onChunk      :: alb_vis_lnd
    dsl4jsb_Real2D_onChunk      :: alb_nir_lnd
    ! VEG_ 2D
    dsl4jsb_Real2D_onChunk      :: growth_respiration
    dsl4jsb_Real2D_onChunk      :: maint_respiration
    dsl4jsb_Real2D_onChunk      :: fmaint_rate_root
    dsl4jsb_Real2D_onChunk      :: n_transform_respiration
    dsl4jsb_Real2D_onChunk      :: n_fixation_respiration
    dsl4jsb_Real2D_onChunk      :: npp
    dsl4jsb_Real2D_onChunk      :: herbivory_leaf_resp
    dsl4jsb_Real2D_onChunk      :: herbivory_fruit_resp
    dsl4jsb_Real2D_onChunk      :: growth_respiration_c13
    dsl4jsb_Real2D_onChunk      :: maint_respiration_c13
    dsl4jsb_Real2D_onChunk      :: growth_respiration_c14
    dsl4jsb_Real2D_onChunk      :: maint_respiration_c14
    dsl4jsb_Real2D_onChunk      :: uptake_nh4
    dsl4jsb_Real2D_onChunk      :: uptake_nh4_n15
    dsl4jsb_Real2D_onChunk      :: uptake_no3
    dsl4jsb_Real2D_onChunk      :: uptake_no3_n15
    dsl4jsb_Real2D_onChunk      :: recycling_leaf_n
    dsl4jsb_Real2D_onChunk      :: recycling_fine_root_n
    dsl4jsb_Real2D_onChunk      :: recycling_heart_wood_n
    dsl4jsb_Real2D_onChunk      :: n_fixation
    dsl4jsb_Real2D_onChunk      :: n_fixation_n15
    dsl4jsb_Real2D_onChunk      :: uptake_po4
    dsl4jsb_Real2D_onChunk      :: recycling_leaf_p
    dsl4jsb_Real2D_onChunk      :: recycling_fine_root_p
    dsl4jsb_Real2D_onChunk      :: recycling_heart_wood_p
    ! VEG_ 3D
    dsl4jsb_Real3D_onChunk      :: root_fraction_sl
    ! SB_ 2D
    dsl4jsb_Real2D_onChunk      :: total_flux_carbon
    dsl4jsb_Real2D_onChunk      :: total_flux_nitrogen
    dsl4jsb_Real2D_onChunk      :: total_flux_phosphorus
    dsl4jsb_Real2D_onChunk      :: total_flux_carbon13
    dsl4jsb_Real2D_onChunk      :: total_flux_carbon14
    dsl4jsb_Real2D_onChunk      :: total_flux_nitrogen15
    dsl4jsb_Real2D_onChunk      :: emission_noy
    dsl4jsb_Real2D_onChunk      :: emission_noy_n15
    dsl4jsb_Real2D_onChunk      :: emission_n2o
    dsl4jsb_Real2D_onChunk      :: emission_n2o_n15
    dsl4jsb_Real2D_onChunk      :: emission_n2
    dsl4jsb_Real2D_onChunk      :: emission_n2_n15
    dsl4jsb_Real2D_onChunk      :: leaching_dom_carbon
    dsl4jsb_Real2D_onChunk      :: leaching_dom_nitrogen
    dsl4jsb_Real2D_onChunk      :: leaching_dom_phosphorus
    dsl4jsb_Real2D_onChunk      :: leaching_dom_carbon13
    dsl4jsb_Real2D_onChunk      :: leaching_dom_carbon14
    dsl4jsb_Real2D_onChunk      :: leaching_dom_nitrogen15
    dsl4jsb_Real2D_onChunk      :: leaching_nh4_solute
    dsl4jsb_Real2D_onChunk      :: leaching_nh4_n15_solute
    dsl4jsb_Real2D_onChunk      :: leaching_no3_solute
    dsl4jsb_Real2D_onChunk      :: leaching_no3_n15_solute
    dsl4jsb_Real2D_onChunk      :: leaching_po4_solute
    dsl4jsb_Real3D_onChunk      :: transport_nh4_solute 
    dsl4jsb_Real3D_onChunk      :: transport_nh4_assoc 
    dsl4jsb_Real3D_onChunk      :: transport_no3_solute
    ! SB_ 3D
    dsl4jsb_Real3D_onChunk      :: het_respiration
    dsl4jsb_Real3D_onChunk      :: het_respiration_c13
    dsl4jsb_Real3D_onChunk      :: het_respiration_c14
    dsl4jsb_Real3D_onChunk      :: net_mineralisation_nh4
    dsl4jsb_Real3D_onChunk      :: net_mineralisation_no3
    dsl4jsb_Real3D_onChunk      :: net_mineralisation_po4
    dsl4jsb_Real3D_onChunk      :: nitrification_no3
    dsl4jsb_Real3D_onChunk      :: nitrification_noy
    dsl4jsb_Real3D_onChunk      :: nitrification_n2o
    dsl4jsb_Real3D_onChunk      :: denitrification_n2
    dsl4jsb_Real3D_onChunk      :: denitrification_noy
    dsl4jsb_Real3D_onChunk      :: denitrification_n2o
    dsl4jsb_Real3D_onChunk      :: asymb_n_fixation 
    dsl4jsb_Real3D_onChunk      :: asymb_n_fixation_n15 
    dsl4jsb_Real3D_onChunk      :: weathering_po4
    dsl4jsb_Real3D_onChunk      :: occlusion_po4
    dsl4jsb_Real3D_onChunk      :: slow_exchange_po4
    dsl4jsb_Real3D_onChunk      :: fast_exchange_po4
    dsl4jsb_Real3D_onChunk      :: fast_adsorpt_po4
    dsl4jsb_Real3D_onChunk      :: fast_adsorpt_nh4
    dsl4jsb_Real3D_onChunk      :: fast_exchange_nh4
    dsl4jsb_Real3D_onChunk      :: biochem_mineralisation_po4
    dsl4jsb_Real3D_onChunk      :: plant_uptake_nh4_sl
    dsl4jsb_Real3D_onChunk      :: plant_uptake_no3_sl
    dsl4jsb_Real3D_onChunk      :: plant_uptake_po4_sl
    dsl4jsb_Real3D_onChunk      :: mycorrhiza_uptake_nh4_sl
    dsl4jsb_Real3D_onChunk      :: mycorrhiza_uptake_no3_sl
    dsl4jsb_Real3D_onChunk      :: mycorrhiza_uptake_norg_sl
    dsl4jsb_Real3D_onChunk      :: mycorrhiza_uptake_po4_sl
    dsl4jsb_Real3D_onChunk      :: microbial_uptake_nh4_sl
    dsl4jsb_Real3D_onChunk      :: microbial_uptake_no3_sl
    dsl4jsb_Real3D_onChunk      :: microbial_uptake_po4_sl
    dsl4jsb_Real3D_onChunk      :: transport_po4_solute
    dsl4jsb_Real3D_onChunk      :: transport_po4_assoc_fast
    dsl4jsb_Real3D_onChunk      :: transport_po4_assoc_slow
    dsl4jsb_Real3D_onChunk      :: transport_po4_occluded
    dsl4jsb_Real3D_onChunk      :: transport_po4_primary
    dsl4jsb_Real3D_onChunk      :: myc_respiration
    dsl4jsb_Real3D_onChunk      :: lateral_loss_dom_carbon_sl
    dsl4jsb_Real3D_onChunk      :: lateral_loss_dom_nitrogen_sl
    dsl4jsb_Real3D_onChunk      :: lateral_loss_dom_phosphorus_sl
    dsl4jsb_Real3D_onChunk      :: lateral_loss_dom_carbon13_sl
    dsl4jsb_Real3D_onChunk      :: lateral_loss_dom_carbon14_sl
    dsl4jsb_Real3D_onChunk      :: lateral_loss_dom_nitrogen15_sl
    dsl4jsb_Real3D_onChunk      :: lateral_loss_nh4_solute_sl
    dsl4jsb_Real3D_onChunk      :: lateral_loss_nh4_n15_solute_sl
    dsl4jsb_Real3D_onChunk      :: lateral_loss_no3_solute_sl
    dsl4jsb_Real3D_onChunk      :: lateral_loss_no3_n15_solute_sl
    dsl4jsb_Real3D_onChunk      :: lateral_loss_po4_solute_sl
    ! ----------------------------------------------------------------------------------------------------- !
    iblk     = options%iblk
    ics      = options%ics
    ice      = options%ice
    nc       = options%nc
    dtime    = options%dtime
    ! ----------------------------------------------------------------------------------------------------- !
    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')
    ! ----------------------------------------------------------------------------------------------------- !
    model         => Get_model(tile%owner_model_id)
    lctlib        => model%lctlib(tile%lcts(1)%lib_id)
    vgrid_soil_w  => Get_vgrid('soil_depth_water')
    nsoil_w       =  vgrid_soil_w%n_levels
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Get_config(SB_)
    dsl4jsb_Get_memory(A2L_) 
    dsl4jsb_Get_memory(L2A_) 
    dsl4jsb_Get_memory(Q_ASSIMI_)
    dsl4jsb_Get_memory(VEG_)
    dsl4jsb_Get_memory(RAD_)
    dsl4jsb_Get_memory(SPQ_)
    dsl4jsb_Get_memory(HYDRO_)
    dsl4jsb_Get_memory(SEB_)
    dsl4jsb_Get_memory(SSE_)
    dsl4jsb_Get_memory(SB_)
    ! ----------------------------------------------------------------------------------------------------- !
    ALLOCATE(hlp1_2d(nc, nsoil_w))
    ! ----------------------------------------------------------------------------------------------------- !
    bgcm_store => tile%bgcm_store
    dsl4jsb_Get_mt2L2D(VEG_BGCM_POOL_ID, veg_pool_mt)
    dsl4jsb_Get_mt2L2D(VEG_BGCM_LITTERFALL_ID, veg_litterfall_mt)
    dsl4jsb_Get_mt2L2D(VEG_BGCM_GROWTH_ID, veg_growth_mt)
    dsl4jsb_Get_mt1L2D(VEG_BGCM_SEED_BED_LITTERFALL_ID, seed_bed_litterfall_mt)
    dsl4jsb_Get_mt1L2D(VEG_BGCM_SEED_BED_GROWTH_ID, seed_bed_growth_mt)
    dsl4jsb_Get_mt1L2D(VEG_BGCM_EXUDATION_ID, veg_exudation_mt)
    dsl4jsb_Get_mt1L2D(VEG_BGCM_RESERVE_USE_ID, veg_reserve_use_mt)
    dsl4jsb_Get_mt2L3D(SB_BGCM_TRANSPORT_ID, sb_transport_mt)
    dsl4jsb_Get_mt1L3D(SB_BGCM_MYCO_EXPORT_ID, sb_mycorrhiza_export_mt)
    dsl4jsb_Get_mt1L3D(SB_BGCM_DELTA_POOLS_ID, sb_delta_pools_mt)
    ! ----------------------------------------------------------------------------------------------------- !
    ! A2L_
    dsl4jsb_Get_var2D_onChunk(A2L_, swpar_srf_down)
    dsl4jsb_Get_var2D_onChunk(A2L_, nhx_deposition)
    dsl4jsb_Get_var2D_onChunk(A2L_, noy_deposition)
    dsl4jsb_Get_var2D_onChunk(A2L_, nhx_n15_deposition)
    dsl4jsb_Get_var2D_onChunk(A2L_, noy_n15_deposition)
    dsl4jsb_Get_var2D_onChunk(A2L_, p_deposition)
    dsl4jsb_Get_var2D_onChunk(A2L_, ndep_mavg)
    dsl4jsb_Get_var2D_onChunk(A2L_, rain)
    dsl4jsb_Get_var2D_onChunk(A2L_, snow)
    dsl4jsb_Get_var2D_onChunk(A2L_, t_air)
    dsl4jsb_Get_var2D_onChunk(A2L_, wind_air)
    dsl4jsb_Get_var2D_onChunk(A2L_, swvis_srf_down)
    dsl4jsb_Get_var2D_onChunk(A2L_, swnir_srf_down)
    dsl4jsb_Get_var2D_onChunk(A2L_, lw_srf_down)
    dsl4jsb_Get_var2D_onChunk(A2L_, co2_mixing_ratio)
    dsl4jsb_Get_var2D_onChunk(A2L_, co2_mixing_ratio_C13)
    dsl4jsb_Get_var2D_onChunk(A2L_, co2_mixing_ratio_C14)
    dsl4jsb_Get_var2D_onChunk(A2L_, press_srf)
    dsl4jsb_Get_var2D_onChunk(A2L_, q_air)
    ! ---------------------------
    ! L2A_
    dsl4jsb_Get_var2D_onChunk(L2A_, q_c_conservation_test)
    ! ---------------------------
    ! Q_ASSIMI_
    dsl4jsb_Get_var2D_onChunk(Q_ASSIMI_, gross_assimilation)
    dsl4jsb_Get_var2D_onChunk(Q_ASSIMI_, gross_assimilation_C13)
    dsl4jsb_Get_var2D_onChunk(Q_ASSIMI_, gross_assimilation_C14)
    dsl4jsb_Get_var2D_onChunk(Q_ASSIMI_, canopy_cond)
    dsl4jsb_Get_var2D_onChunk(Q_ASSIMI_, aerodyn_cond)
    ! ---------------------------
    ! HYDRO_ 2D
    dsl4jsb_Get_var2D_onChunk(HYDRO_,    num_sl_above_bedrock)        ! in
    dsl4jsb_Get_var2D_onChunk(HYDRO_,    transpiration)
    dsl4jsb_Get_var2D_onChunk(HYDRO_,    evapo_skin)
    IF (model%config%use_soil_phys_jsbach) THEN
      dsl4jsb_Get_var2D_onChunk(HYDRO_,    evapopot)
      dsl4jsb_Get_var2D_onChunk(HYDRO_,    drainage)
    END IF
    dsl4jsb_Get_var2D_onChunk(HYDRO_,    runoff)
    ! HYDRO_ 3D
    dsl4jsb_Get_var3D_onChunk(HYDRO_,    soil_depth_sl)               ! in
    ! ---------------------------
    ! SEB_ 2D
    dsl4jsb_Get_var2D_onChunk(SEB_,      sensible_hflx)               ! in
    dsl4jsb_Get_var2D_onChunk(SEB_,      latent_hflx)                 ! in
    ! ---------------------------
    ! SSE_ 2D
    dsl4jsb_Get_var2D_onChunk(SSE_,      grnd_hflx)                   ! in
    ! ---------------------------
    ! SPQ_ 2D
    IF (.NOT. model%config%use_soil_phys_jsbach) THEN
      dsl4jsb_Get_var2D_onChunk(SPQ_, evaporation)           ! jsbach var: H evapo_skin (Evaporation from skin reservoir)
      dsl4jsb_Get_var2D_onChunk(SPQ_, interception)          ! jsbach var: local var 'rain_to_skinres' in calc_surface_hydrology_land() ?
      dsl4jsb_Get_var2D_onChunk(SPQ_, evapopot)              ! jsbach var: H evapopot (Potential evaporation)
      dsl4jsb_Get_var2D_onChunk(SPQ_, drainage)              ! jsbach var: H 
      dsl4jsb_Get_var2D_onChunk(SPQ_, srf_runoff)            ! jsbach var: H runoff (Surface runoff)
      dsl4jsb_Get_var2D_onChunk(SPQ_, gw_runoff)             ! jsbach var: H drainage (Subsurface drainage (horizontal) from soil layers)
      dsl4jsb_Get_var2D_onChunk(SPQ_, sensible_heat_flx)     ! jsbach var: SEB sensible_hflx (Sensible heat flux [W m-2])
      dsl4jsb_Get_var2D_onChunk(SPQ_, latent_heat_flx)       ! jsbach var: SEB latent_hflx   (Latent heat flux [W m-2])
      dsl4jsb_Get_var2D_onChunk(SPQ_, ground_heat_flx)       ! jsbach var: SSE grnd_hflx     (ground heat flux J m-2 s-1)
    END IF
    ! ---------------------------
    ! RAD_
    dsl4jsb_Get_var2D_onChunk(RAD_, albedo_noon)
    dsl4jsb_Get_var2D_onChunk(RAD_, rad_srf_net)
    dsl4jsb_Get_var2D_onChunk(RAD_, fapar_noon)
    dsl4jsb_Get_var2D_onChunk(RAD_, alb_vis_lnd)
    dsl4jsb_Get_var2D_onChunk(RAD_, alb_nir_lnd)
    ! ---------------------------
    ! VEG_ 2D
    dsl4jsb_Get_var2D_onChunk(VEG_, growth_respiration)
    dsl4jsb_Get_var2D_onChunk(VEG_, maint_respiration)
    dsl4jsb_Get_var2D_onChunk(VEG_, fmaint_rate_root)
    dsl4jsb_Get_var2D_onChunk(VEG_, n_transform_respiration)
    dsl4jsb_Get_var2D_onChunk(VEG_, n_fixation_respiration)
    dsl4jsb_Get_var2D_onChunk(VEG_, npp)
    dsl4jsb_Get_var2D_onChunk(VEG_, herbivory_leaf_resp)
    dsl4jsb_Get_var2D_onChunk(VEG_, herbivory_fruit_resp)
    dsl4jsb_Get_var2D_onChunk(VEG_, growth_respiration_c13)
    dsl4jsb_Get_var2D_onChunk(VEG_, maint_respiration_c13)
    dsl4jsb_Get_var2D_onChunk(VEG_, growth_respiration_c14)
    dsl4jsb_Get_var2D_onChunk(VEG_, maint_respiration_c14)
    dsl4jsb_Get_var2D_onChunk(VEG_, uptake_nh4)
    dsl4jsb_Get_var2D_onChunk(VEG_, uptake_nh4_n15)
    dsl4jsb_Get_var2D_onChunk(VEG_, uptake_no3)
    dsl4jsb_Get_var2D_onChunk(VEG_, uptake_no3_n15)
    dsl4jsb_Get_var2D_onChunk(VEG_, recycling_leaf_n)
    dsl4jsb_Get_var2D_onChunk(VEG_, recycling_fine_root_n)
    dsl4jsb_Get_var2D_onChunk(VEG_, recycling_heart_wood_n)
    dsl4jsb_Get_var2D_onChunk(VEG_, n_fixation)
    dsl4jsb_Get_var2D_onChunk(VEG_, n_fixation_n15)
    dsl4jsb_Get_var2D_onChunk(VEG_, uptake_po4)
    dsl4jsb_Get_var2D_onChunk(VEG_, recycling_leaf_p)
    dsl4jsb_Get_var2D_onChunk(VEG_, recycling_fine_root_p)
    dsl4jsb_Get_var2D_onChunk(VEG_, recycling_heart_wood_p)
    ! VEG_ 3D
    dsl4jsb_Get_var3D_onChunk(VEG_, root_fraction_sl)
    ! ---------------------------
    ! SB_ 2D
    dsl4jsb_Get_var2D_onChunk(SB_,  total_flux_carbon)    
    dsl4jsb_Get_var2D_onChunk(SB_,  total_flux_nitrogen)  
    dsl4jsb_Get_var2D_onChunk(SB_,  total_flux_phosphorus)
    dsl4jsb_Get_var2D_onChunk(SB_,  total_flux_carbon13)  
    dsl4jsb_Get_var2D_onChunk(SB_,  total_flux_carbon14)  
    dsl4jsb_Get_var2D_onChunk(SB_,  total_flux_nitrogen15)
    dsl4jsb_Get_var2D_onChunk(SB_,  emission_noy)
    dsl4jsb_Get_var2D_onChunk(SB_,  emission_noy_n15)
    dsl4jsb_Get_var2D_onChunk(SB_,  emission_n2o)
    dsl4jsb_Get_var2D_onChunk(SB_,  emission_n2o_n15)
    dsl4jsb_Get_var2D_onChunk(SB_,  emission_n2)
    dsl4jsb_Get_var2D_onChunk(SB_,  emission_n2_n15)
    dsl4jsb_Get_var2D_onChunk(SB_,  leaching_dom_carbon)
    dsl4jsb_Get_var2D_onChunk(SB_,  leaching_dom_nitrogen)
    dsl4jsb_Get_var2D_onChunk(SB_,  leaching_dom_phosphorus)
    dsl4jsb_Get_var2D_onChunk(SB_,  leaching_dom_carbon13)
    dsl4jsb_Get_var2D_onChunk(SB_,  leaching_dom_carbon14)
    dsl4jsb_Get_var2D_onChunk(SB_,  leaching_dom_nitrogen15)
    dsl4jsb_Get_var2D_onChunk(SB_,  leaching_nh4_solute)
    dsl4jsb_Get_var2D_onChunk(SB_,  leaching_nh4_n15_solute)
    dsl4jsb_Get_var2D_onChunk(SB_,  leaching_no3_solute)
    dsl4jsb_Get_var2D_onChunk(SB_,  leaching_no3_n15_solute)
    dsl4jsb_Get_var2D_onChunk(SB_,  leaching_po4_solute)
    ! SB_ 3D
    dsl4jsb_Get_var3D_onChunk(SB_,  transport_nh4_solute)
    dsl4jsb_Get_var3D_onChunk(SB_,  transport_nh4_assoc)
    dsl4jsb_Get_var3D_onChunk(SB_,  transport_no3_solute)
    dsl4jsb_Get_var3D_onChunk(SB_,  lateral_loss_dom_carbon_sl)
    dsl4jsb_Get_var3D_onChunk(SB_,  lateral_loss_dom_nitrogen_sl)
    dsl4jsb_Get_var3D_onChunk(SB_,  lateral_loss_dom_phosphorus_sl) 
    dsl4jsb_Get_var3D_onChunk(SB_,  lateral_loss_dom_carbon13_sl) 
    dsl4jsb_Get_var3D_onChunk(SB_,  lateral_loss_dom_carbon14_sl) 
    dsl4jsb_Get_var3D_onChunk(SB_,  lateral_loss_dom_nitrogen15_sl) 
    dsl4jsb_Get_var3D_onChunk(SB_,  lateral_loss_nh4_solute_sl) 
    dsl4jsb_Get_var3D_onChunk(SB_,  lateral_loss_nh4_n15_solute_sl) 
    dsl4jsb_Get_var3D_onChunk(SB_,  lateral_loss_no3_solute_sl) 
    dsl4jsb_Get_var3D_onChunk(SB_,  lateral_loss_no3_n15_solute_sl) 
    dsl4jsb_Get_var3D_onChunk(SB_,  lateral_loss_po4_solute_sl) 
    dsl4jsb_Get_var3D_onChunk(SB_,  het_respiration)
    dsl4jsb_Get_var3D_onChunk(SB_,  het_respiration_c13)
    dsl4jsb_Get_var3D_onChunk(SB_,  het_respiration_c14)
    dsl4jsb_Get_var3D_onChunk(SB_,  net_mineralisation_nh4)
    dsl4jsb_Get_var3D_onChunk(SB_,  net_mineralisation_no3)
    dsl4jsb_Get_var3D_onChunk(SB_,  net_mineralisation_po4)
    dsl4jsb_Get_var3D_onChunk(SB_,  nitrification_no3)
    dsl4jsb_Get_var3D_onChunk(SB_,  nitrification_noy)
    dsl4jsb_Get_var3D_onChunk(SB_,  nitrification_n2o)
    dsl4jsb_Get_var3D_onChunk(SB_,  denitrification_n2)
    dsl4jsb_Get_var3D_onChunk(SB_,  denitrification_noy)
    dsl4jsb_Get_var3D_onChunk(SB_,  denitrification_n2o)
    dsl4jsb_Get_var3D_onChunk(SB_,  asymb_n_fixation)
    dsl4jsb_Get_var3D_onChunk(SB_,  asymb_n_fixation_n15)
    dsl4jsb_Get_var3D_onChunk(SB_,  slow_exchange_po4)
    dsl4jsb_Get_var3D_onChunk(SB_,  fast_exchange_po4)
    dsl4jsb_Get_var3D_onChunk(SB_,  fast_adsorpt_po4)
    dsl4jsb_Get_var3D_onChunk(SB_,  fast_adsorpt_nh4)
    dsl4jsb_Get_var3D_onChunk(SB_,  fast_exchange_nh4)
    dsl4jsb_Get_var3D_onChunk(SB_,  weathering_po4)
    dsl4jsb_Get_var3D_onChunk(SB_,  occlusion_po4)
    dsl4jsb_Get_var3D_onChunk(SB_,  biochem_mineralisation_po4)
    dsl4jsb_Get_var3D_onChunk(SB_,  plant_uptake_nh4_sl)
    dsl4jsb_Get_var3D_onChunk(SB_,  plant_uptake_no3_sl)
    dsl4jsb_Get_var3D_onChunk(SB_,  plant_uptake_po4_sl)
    dsl4jsb_Get_var3D_onChunk(SB_,  mycorrhiza_uptake_nh4_sl)
    dsl4jsb_Get_var3D_onChunk(SB_,  mycorrhiza_uptake_no3_sl)
    dsl4jsb_Get_var3D_onChunk(SB_,  mycorrhiza_uptake_norg_sl)
    dsl4jsb_Get_var3D_onChunk(SB_,  mycorrhiza_uptake_po4_sl)
    dsl4jsb_Get_var3D_onChunk(SB_,  microbial_uptake_nh4_sl)
    dsl4jsb_Get_var3D_onChunk(SB_,  microbial_uptake_no3_sl)
    dsl4jsb_Get_var3D_onChunk(SB_,  microbial_uptake_po4_sl)
    dsl4jsb_Get_var3D_onChunk(SB_,  transport_po4_solute)
    dsl4jsb_Get_var3D_onChunk(SB_,  transport_po4_assoc_fast)
    dsl4jsb_Get_var3D_onChunk(SB_,  transport_po4_assoc_slow)
    dsl4jsb_Get_var3D_onChunk(SB_,  transport_po4_occluded)
    dsl4jsb_Get_var3D_onChunk(SB_,  transport_po4_primary)
    dsl4jsb_Get_var3D_onChunk(SB_,  myc_respiration)
    ! ----------------------------------------------------------------------------------------------------- !

    !>0.7 the point in ics:ice of which data will be written to file; per default = 1
    !>
    !> is actually defined and declared globally in his module above, but ifort complains anyway
    !>
    ic = 1

    !>0.8 modify 'dtime_one' factor accordingly to output interval
    !>
    ! accumulate fluxes across the timesteps of one interval:
    !   for timestep-interval leave it to seconds (which is the basic timestep for fluxes), i.e. multiply by 1
    !   for any interval larger timestep multiply by dtime_one
    IF (output_interval_flux < 1.0_wp) THEN
      dtime_one = 1.0_wp                 ! interval = one timestep
    ELSE
      dtime_one = dtime / 1000000.0_wp   ! interval > one timestep
    END IF
    ! veg_flux output factor [veg_flux * f_veg_flux_output]
    f_veg_flux_output = dtime_one * 1000000._wp / dtime

    !>0.9 pass value of "set *_acc variables to zero"
    !>
    set_zero_flux_acc_var_local = set_zero_flux_acc_variables


    !>1.0 accumulate fluxes across the timesteps of one interval
    !>
    !>  this is done per default
    !>    and is only not done in case 'set_zero_flux_acc_variables is TRUE'
    !>
    IF (.NOT. set_zero_flux_acc_var_local) THEN
      ! counter of flux accumulation timesteps
      timestep_counter_flux_acc = timestep_counter_flux_acc + 1.0_wp
      ! Vegetation
      gpp_acc              = gpp_acc               + gross_assimilation(ic)        * dtime_one
      growth_resp_acc      = growth_resp_acc       + growth_respiration(ic)        * dtime_one
      maint_resp_acc       = maint_resp_acc        + maint_respiration(ic)         * dtime_one
      n_transform_resp_acc = n_transform_resp_acc  + n_transform_respiration(ic)   * dtime_one
      n_fixation_resp_acc  = n_fixation_resp_acc   + n_fixation_respiration(ic)    * dtime_one
      exudation_c_acc      = exudation_c_acc       + veg_exudation_mt(ixC, ic)     * dtime_one *1000000._wp/dtime      
      npp_acc              = npp_acc               + npp(ic)                       * dtime_one
      leaf_herb_resp_acc   = leaf_herb_resp_acc    + herbivory_leaf_resp(ic)       * dtime_one
      fruit_herb_resp_acc  = fruit_herb_resp_acc   + herbivory_fruit_resp(ic)      * dtime_one
      res_use_c_acc        = res_use_c_acc         + veg_reserve_use_mt(ixC, ic)   * dtime_one *1000000._wp/dtime
      growth_c_acc         = growth_c_acc          + (veg_growth_mt(ix_leaf, ixC, ic)         &
        &                                           + veg_growth_mt(ix_fine_root, ixC, ic)    &
        &                                           + veg_growth_mt(ix_coarse_root, ixC, ic)  &
        &                                           + veg_growth_mt(ix_sap_wood, ixC, ic)     &
        &                                           + veg_growth_mt(ix_heart_wood, ixC, ic)   &
        &                                           + veg_growth_mt(ix_labile, ixC, ic)       &
        &                                           + veg_growth_mt(ix_reserve, ixC, ic)      &
        &                                           + veg_growth_mt(ix_fruit, ixC, ic))       &
        &                                           * dtime_one * 1000000._wp / dtime
      growth_swood_c13_acc = growth_swood_c13_acc  + veg_growth_mt(ix_sap_wood, ixC13, ic) * dtime_one *1000000._wp/dtime
      growth_swood_c14_acc = growth_swood_c14_acc  + veg_growth_mt(ix_sap_wood, ixC14, ic) * dtime_one *1000000._wp/dtime
      veg_growth_leaf_c_acc            = veg_growth_leaf_c_acc            &
        &                                + veg_growth_mt(ix_leaf, ixC, ic) * f_veg_flux_output
      veg_growth_fine_root_c_acc       = veg_growth_fine_root_c_acc       &
        &                                + veg_growth_mt(ix_fine_root, ixC, ic) * f_veg_flux_output
      veg_growth_coarse_root_c_acc     = veg_growth_coarse_root_c_acc     &
        &                                + veg_growth_mt(ix_coarse_root, ixC, ic) * f_veg_flux_output
      veg_growth_sap_wood_c_acc        = veg_growth_sap_wood_c_acc        &
        &                                + veg_growth_mt(ix_sap_wood, ixC, ic) * f_veg_flux_output
      veg_growth_heart_wood_c_acc      = veg_growth_heart_wood_c_acc      &
        &                                + veg_growth_mt(ix_heart_wood, ixC, ic) * f_veg_flux_output
      veg_growth_labile_c_acc          = veg_growth_labile_c_acc          &
        &                                + veg_growth_mt(ix_labile, ixC, ic) * f_veg_flux_output
      veg_growth_reserve_c_acc         = veg_growth_reserve_c_acc         &
        &                                + veg_growth_mt(ix_reserve, ixC, ic) * f_veg_flux_output
      veg_growth_fruit_c_acc           = veg_growth_fruit_c_acc           &
        &                                + veg_growth_mt(ix_fruit, ixC, ic) * f_veg_flux_output
      veg_growth_seed_bed_c_acc        = veg_growth_seed_bed_c_acc        &
        &                                + seed_bed_growth_mt(ixC, ic) * f_veg_flux_output
      veg_growth_leaf_n_acc            = veg_growth_leaf_n_acc            &
        &                                + veg_growth_mt(ix_leaf, ixN, ic) * f_veg_flux_output
      veg_growth_fine_root_n_acc       = veg_growth_fine_root_n_acc       &
        &                                + veg_growth_mt(ix_fine_root, ixN, ic) * f_veg_flux_output
      veg_growth_coarse_root_n_acc     = veg_growth_coarse_root_n_acc     &
        &                                + veg_growth_mt(ix_coarse_root, ixN, ic) * f_veg_flux_output
      veg_growth_sap_wood_n_acc        = veg_growth_sap_wood_n_acc        &
        &                                + veg_growth_mt(ix_sap_wood, ixN, ic) * f_veg_flux_output
      veg_growth_heart_wood_n_acc      = veg_growth_heart_wood_n_acc      &
        &                                + veg_growth_mt(ix_heart_wood, ixN, ic) * f_veg_flux_output
      veg_growth_labile_n_acc          = veg_growth_labile_n_acc          &
        &                                + veg_growth_mt(ix_labile, ixN, ic) * f_veg_flux_output
      veg_growth_reserve_n_acc         = veg_growth_reserve_n_acc         &
        &                                + veg_growth_mt(ix_reserve, ixN, ic) * f_veg_flux_output
      veg_growth_fruit_n_acc           = veg_growth_fruit_n_acc           &
        &                                + veg_growth_mt(ix_fruit, ixN, ic) * f_veg_flux_output
      veg_growth_seed_bed_n_acc        = veg_growth_seed_bed_n_acc        &
        &                                + seed_bed_growth_mt(ixN, ic) * f_veg_flux_output
      veg_growth_leaf_p_acc            = veg_growth_leaf_p_acc            &
        &                                + veg_growth_mt(ix_leaf, ixP, ic) * f_veg_flux_output
      veg_growth_fine_root_p_acc       = veg_growth_fine_root_p_acc       &
        &                                + veg_growth_mt(ix_fine_root, ixP, ic) * f_veg_flux_output
      veg_growth_coarse_root_p_acc     = veg_growth_coarse_root_p_acc     &
        &                                + veg_growth_mt(ix_coarse_root, ixP, ic) * f_veg_flux_output
      veg_growth_sap_wood_p_acc        = veg_growth_sap_wood_p_acc        &
        &                                + veg_growth_mt(ix_sap_wood, ixP, ic) * f_veg_flux_output
      veg_growth_heart_wood_p_acc      = veg_growth_heart_wood_p_acc      &
        &                                + veg_growth_mt(ix_heart_wood, ixP, ic) * f_veg_flux_output
      veg_growth_labile_p_acc          = veg_growth_labile_p_acc          &
        &                                + veg_growth_mt(ix_labile, ixP, ic) * f_veg_flux_output
      veg_growth_reserve_p_acc         = veg_growth_reserve_p_acc         &
        &                                + veg_growth_mt(ix_reserve, ixP, ic) * f_veg_flux_output
      veg_growth_fruit_p_acc           = veg_growth_fruit_p_acc           &
        &                                + veg_growth_mt(ix_fruit, ixP, ic) * f_veg_flux_output
      veg_growth_seed_bed_p_acc        = veg_growth_seed_bed_p_acc        &
        &                                + seed_bed_growth_mt(ixP, ic) * f_veg_flux_output
      veg_litterfall_leaf_c_acc        = veg_litterfall_leaf_c_acc        &
        &                                + veg_litterfall_mt(ix_leaf, ixC, ic) * f_veg_flux_output
      veg_litterfall_fine_root_c_acc   = veg_litterfall_fine_root_c_acc   &
        &                                + veg_litterfall_mt(ix_fine_root, ixC, ic) * f_veg_flux_output
      veg_litterfall_coarse_root_c_acc = veg_litterfall_coarse_root_c_acc &
        &                                + veg_litterfall_mt(ix_coarse_root, ixC, ic) * f_veg_flux_output
      veg_litterfall_sap_wood_c_acc    = veg_litterfall_sap_wood_c_acc    &
        &                                + veg_litterfall_mt(ix_sap_wood, ixC, ic) * f_veg_flux_output
      veg_litterfall_heart_wood_c_acc  = veg_litterfall_heart_wood_c_acc  &
        &                                + veg_litterfall_mt(ix_heart_wood, ixC, ic) * f_veg_flux_output
      veg_litterfall_labile_c_acc      = veg_litterfall_labile_c_acc      &
        &                                + veg_litterfall_mt(ix_labile, ixC, ic) * f_veg_flux_output
      veg_litterfall_reserve_c_acc     = veg_litterfall_reserve_c_acc     &
        &                                + veg_litterfall_mt(ix_reserve, ixC, ic) * f_veg_flux_output
      veg_litterfall_fruit_c_acc       = veg_litterfall_fruit_c_acc       &
        &                                + veg_litterfall_mt(ix_fruit, ixC, ic) * f_veg_flux_output
      veg_litterfall_seed_bed_c_acc    = veg_litterfall_seed_bed_c_acc    &
        &                                + seed_bed_litterfall_mt(ixC, ic) * f_veg_flux_output
      veg_litterfall_leaf_n_acc        = veg_litterfall_leaf_n_acc        &
        &                                + veg_litterfall_mt(ix_leaf, ixN, ic) * f_veg_flux_output
      veg_litterfall_fine_root_n_acc   = veg_litterfall_fine_root_n_acc   &
        &                                + veg_litterfall_mt(ix_fine_root, ixN, ic) * f_veg_flux_output
      veg_litterfall_coarse_root_n_acc = veg_litterfall_coarse_root_n_acc &
        &                                + veg_litterfall_mt(ix_coarse_root, ixN, ic) * f_veg_flux_output
      veg_litterfall_sap_wood_n_acc    = veg_litterfall_sap_wood_n_acc    &
        &                                + veg_litterfall_mt(ix_sap_wood, ixN, ic) * f_veg_flux_output
      veg_litterfall_heart_wood_n_acc  = veg_litterfall_heart_wood_n_acc  &
        &                                + veg_litterfall_mt(ix_heart_wood, ixN, ic) * f_veg_flux_output
      veg_litterfall_labile_n_acc      = veg_litterfall_labile_n_acc      &
        &                                + veg_litterfall_mt(ix_labile, ixN, ic) * f_veg_flux_output
      veg_litterfall_reserve_n_acc     = veg_litterfall_reserve_n_acc     &
        &                                + veg_litterfall_mt(ix_reserve, ixN, ic) * f_veg_flux_output
      veg_litterfall_fruit_n_acc       = veg_litterfall_fruit_n_acc       &
        &                                + veg_litterfall_mt(ix_fruit, ixN, ic) * f_veg_flux_output
      veg_litterfall_seed_bed_n_acc    = veg_litterfall_seed_bed_n_acc    &
        &                                + seed_bed_litterfall_mt(ixN, ic) * f_veg_flux_output
      veg_litterfall_leaf_p_acc        = veg_litterfall_leaf_p_acc        &
        &                                + veg_litterfall_mt(ix_leaf, ixP, ic) * f_veg_flux_output
      veg_litterfall_fine_root_p_acc   = veg_litterfall_fine_root_p_acc   &
        &                                + veg_litterfall_mt(ix_fine_root, ixP, ic) * f_veg_flux_output
      veg_litterfall_coarse_root_p_acc = veg_litterfall_coarse_root_p_acc &
        &                                + veg_litterfall_mt(ix_coarse_root, ixP, ic) * f_veg_flux_output
      veg_litterfall_sap_wood_p_acc    = veg_litterfall_sap_wood_p_acc    &
        &                                + veg_litterfall_mt(ix_sap_wood, ixP, ic) * f_veg_flux_output
      veg_litterfall_heart_wood_p_acc  = veg_litterfall_heart_wood_p_acc  &
        &                                + veg_litterfall_mt(ix_heart_wood, ixP, ic) * f_veg_flux_output
      veg_litterfall_labile_p_acc      = veg_litterfall_labile_p_acc      &
        &                                + veg_litterfall_mt(ix_labile, ixP, ic) * f_veg_flux_output
      veg_litterfall_reserve_p_acc     = veg_litterfall_reserve_p_acc     &
        &                                + veg_litterfall_mt(ix_reserve, ixP, ic) * f_veg_flux_output
      veg_litterfall_fruit_p_acc       = veg_litterfall_fruit_p_acc       &
        &                                + veg_litterfall_mt(ix_fruit, ixP, ic) * f_veg_flux_output
      veg_litterfall_seed_bed_p_acc    = veg_litterfall_seed_bed_p_acc    &
        &                                + seed_bed_litterfall_mt(ixP, ic) * f_veg_flux_output
      gpp_C13_acc          = gpp_C13_acc           + gross_assimilation_C13(ic)    * dtime_one
      growth_resp_C13_acc  = growth_resp_C13_acc   + growth_respiration_c13(ic)    * dtime_one
      maint_resp_C13_acc   = maint_resp_C13_acc    + maint_respiration_c13(ic)     * dtime_one
      gpp_C14_acc          = gpp_C14_acc           + gross_assimilation_C14(ic)    * dtime_one
      growth_resp_C14_acc  = growth_resp_C14_acc   + growth_respiration_c14(ic)    * dtime_one
      maint_resp_C14_acc   = maint_resp_C14_acc    + maint_respiration_c14(ic)     * dtime_one
      uptake_nh4_acc       = uptake_nh4_acc        + uptake_nh4(ic)                * dtime_one
      uptake_nh4_n15_acc   = uptake_nh4_n15_acc    + uptake_nh4_n15(ic)            * dtime_one
      uptake_no3_acc       = uptake_no3_acc        + uptake_no3(ic)                * dtime_one
      uptake_no3_n15_acc   = uptake_no3_n15_acc    + uptake_no3_n15(ic)            * dtime_one
      recycling_n_acc      = recycling_n_acc       + ( recycling_leaf_n(ic) + recycling_fine_root_n(ic) &
        &                                                                   + recycling_heart_wood_n(ic)) &
        &                                                                   * dtime_one * 1000000._wp / dtime
      res_use_n_acc        = res_use_n_acc         + veg_reserve_use_mt(ixN, ic)  * dtime_one *1000000._wp/dtime
      n_fixation_acc       = n_fixation_acc        + n_fixation(ic)                * dtime_one 
      n_fixation_n15_acc   = n_fixation_n15_acc    + n_fixation_n15(ic)            * dtime_one
      growth_n_acc         = growth_n_acc          + (veg_growth_mt(ix_leaf, ixN, ic)        + & 
                                                      veg_growth_mt(ix_fine_root, ixN, ic)   + &
                                                      veg_growth_mt(ix_coarse_root, ixN, ic) + & 
                                                      veg_growth_mt(ix_sap_wood, ixN, ic)    + & 
                                                      veg_growth_mt(ix_heart_wood, ixN, ic)  + & 
                                                      veg_growth_mt(ix_labile, ixN, ic)      + & 
                                                      veg_growth_mt(ix_reserve, ixN, ic)     + & 
                                                      veg_growth_mt(ix_fruit, ixN, ic)) * dtime_one *1000000._wp/dtime
      growth_swood_n15_acc = growth_swood_n15_acc  + veg_growth_mt(ix_sap_wood, ixN15, ic) * dtime_one *1000000._wp/dtime
      uptake_po4_acc       = uptake_po4_acc        + uptake_po4(ic)                * dtime_one
      recycling_p_acc      = recycling_p_acc       + (recycling_leaf_p(ic) + recycling_fine_root_p(ic) + &
                                                      recycling_heart_wood_p(ic))  * dtime_one *1000000._wp/dtime
      growth_p_acc         = growth_p_acc          + (veg_growth_mt(ix_leaf, ixP, ic) + & 
                                                      veg_growth_mt(ix_fine_root, ixP, ic)   + &
                                                      veg_growth_mt(ix_coarse_root, ixP, ic) + & 
                                                      veg_growth_mt(ix_sap_wood, ixP, ic)    + & 
                                                      veg_growth_mt(ix_heart_wood, ixP, ic)  + & 
                                                      veg_growth_mt(ix_labile, ixP, ic)      + & 
                                                      veg_growth_mt(ix_reserve, ixP, ic)     + & 
                                                      veg_growth_mt(ix_fruit, ixP, ic)) * dtime_one *1000000._wp/dtime     
      res_use_p_acc        = res_use_p_acc         + veg_reserve_use_mt(ixP, ic)  * dtime_one *1000000._wp/dtime
      ! use jsbach soil-physics or quincy SPQ_
      !
      ! quincy
      IF (.NOT. model%config%use_soil_phys_jsbach) THEN
        ! if interval == timestep (output_interval_flux is smaller 1.0 only for timestep interval)
        IF (output_interval_flux < 1.0_wp) THEN
          transpiration_acc = transpiration_acc + transpiration(ic)
          evaporation_acc   = evaporation_acc   + evaporation(ic)
          interception_acc  = interception_acc  + interception(ic)
          evapopot_acc      = evapopot_acc      + evapopot(ic)
          drainage_acc      = drainage_acc      + drainage(ic)
          srf_runoff_acc    = srf_runoff_acc    + srf_runoff(ic)
          gw_runoff_acc     = gw_runoff_acc     + gw_runoff(ic)
        ELSE
          transpiration_acc = transpiration_acc + transpiration(ic)                * dtime
          evaporation_acc   = evaporation_acc   + evaporation(ic)                  * dtime
          interception_acc  = interception_acc  + interception(ic)                 * dtime
          evapopot_acc      = evapopot_acc      + evapopot(ic)                     * dtime
          drainage_acc      = drainage_acc      + drainage(ic)                     * dtime
          srf_runoff_acc    = srf_runoff_acc    + srf_runoff(ic)                   * dtime
          gw_runoff_acc     = gw_runoff_acc     + gw_runoff(ic)                    * dtime
        END IF
        qh_acc              = qh_acc            - sensible_heat_flx(ic)            * dtime_one               
        qle_acc             = qle_acc           - latent_heat_flx(ic)              * dtime_one               
        qg_acc              = qg_acc            - ground_heat_flx(ic)              * dtime_one               
      ! jsbach
      ELSE
        ! if interval == timestep (output_interval_flux is smaller 1.0 only for timestep interval)
        IF (output_interval_flux < 1.0_wp) THEN 
          transpiration_acc = transpiration_acc + (transpiration(ic) * (-1.0_wp))
          evaporation_acc   = evaporation_acc   + (evapo_skin(ic) * (-1.0_wp))
          interception_acc  = 0.0_wp
          evapopot_acc      = evapopot_acc      + evapopot(ic)
          drainage_acc      = 0.0_wp
          srf_runoff_acc    = srf_runoff_acc    + runoff(ic)
          gw_runoff_acc     = gw_runoff_acc     + drainage(ic)                                  ! = SUM(drainage_sl) + drain_bot
        ELSE
          transpiration_acc = transpiration_acc + (transpiration(ic) * (-1.0_wp))  * dtime
          evaporation_acc   = evaporation_acc   + (evapo_skin(ic) * (-1.0_wp))     * dtime
          interception_acc  = 0.0_wp
          evapopot_acc      = evapopot_acc      + evapopot(ic)                     * dtime
          drainage_acc      = 0.0_wp
          srf_runoff_acc    = srf_runoff_acc    + runoff(ic)                       * dtime
          gw_runoff_acc     = gw_runoff_acc     + drainage(ic)                     * dtime      ! = SUM(drainage_sl) + drain_bot
        END IF
        qh_acc              = qh_acc            - sensible_hflx(ic)                * dtime_one               
        qle_acc             = qle_acc           - latent_hflx(ic)                  * dtime_one               
        qg_acc              = qg_acc            - grnd_hflx(ic)                    * dtime_one               
      END IF    ! IF (.NOT. model%config%use_soil_phys_jsbach)
      ! precipitation accumulation depends on output interval
      ! if interval == timestep (output_interval_flux is smaller 1.0 only for timestep interval)
      IF (output_interval_flux < 1.0_wp) THEN
        precip_acc = precip_acc + (rain(ic) + snow(ic))  
      ELSE
        precip_acc = precip_acc + (rain(ic) + snow(ic)) * dtime
      END IF 
      gc_acc = gc_acc + canopy_cond(ic)                                                      
      ga_acc = ga_acc + aerodyn_cond(ic) 
      ! albedo is accumulated only if there is sunlight
      IF (alb_vis_lnd(ic) > 0.0_wp) THEN
        alb_vis_acc                 = alb_vis_acc + alb_vis_lnd(ic)
        alb_nir_acc                 = alb_nir_acc + alb_nir_lnd(ic)
        timestep_counter_albedo_acc = timestep_counter_albedo_acc + 1._wp
      END IF
      rnet_acc         = rnet_acc          + rad_srf_net(ic) * dtime_one               
      ! Mycorrhiza
      myc_export_c_acc      = myc_export_c_acc      + SUM(sb_mycorrhiza_export_mt(ixC, ic, :) * &
                                                        soil_depth_sl(ic,:)) * dtime_one *1000000._wp/dtime
      myc_export_n_acc      = myc_export_n_acc      + SUM(sb_mycorrhiza_export_mt(ixN, ic, :) * &
                                                        soil_depth_sl(ic,:)) * dtime_one *1000000._wp/dtime
      myc_export_p_acc      = myc_export_p_acc      + SUM(sb_mycorrhiza_export_mt(ixP, ic, :) * &
                                                        soil_depth_sl(ic,:)) * dtime_one *1000000._wp/dtime
      myc_export_c13_acc    = myc_export_c13_acc    + SUM(sb_mycorrhiza_export_mt(ixC13, ic, :) * &
                                                        soil_depth_sl(ic,:)) * dtime_one *1000000._wp/dtime
      myc_export_c14_acc    = myc_export_c14_acc    + SUM(sb_mycorrhiza_export_mt(ixC14, ic, :) * &
                                                        soil_depth_sl(ic,:)) * dtime_one *1000000._wp/dtime
      myc_export_n15_acc    = myc_export_n15_acc    + SUM(sb_mycorrhiza_export_mt(ixN15, ic, :) * &
                                                        soil_depth_sl(ic,:)) * dtime_one *1000000._wp/dtime
      myc_respiration_c_acc = myc_respiration_c_acc + SUM(myc_respiration(ic,:) * soil_depth_sl(ic,:)) * &
                              molar_mass_C * dtime_one * 1000000._wp / dtime
      f_cm_som_acc          = f_cm_som_acc + organic_cn * molar_mass_C / molar_mass_N * &
                              SUM(mycorrhiza_uptake_norg_sl(ic,:) * soil_depth_sl(ic,:)) * &
                              dtime_one * 1000000._wp / dtime
      ! matter conservation test
      c_cons_test_acc       = c_cons_test_acc + q_c_conservation_test(ic) * dtime_one
      ! atmospheric forcing (A2L_)
      t_air_acc                = t_air_acc                + t_air(ic)
      wind_air_acc             = wind_air_acc             + wind_air(ic)
      swvis_srf_down_acc       = swvis_srf_down_acc       + swvis_srf_down(ic)
      swnir_srf_down_acc       = swnir_srf_down_acc       + swnir_srf_down(ic)
      lw_srf_down_acc          = lw_srf_down_acc          + lw_srf_down(ic)
      ! NOTE: in the forcing we read VPD air, therefore we should also only output VPD air
      vpd_acc                  = vpd_acc + ((qsat_mixed(t_air(ic), press_srf(ic)) - q_air(ic)) &
        &                                   / eps_vpd * press_srf(ic) / 1000.0_wp)
      co2_mixing_ratio_acc     = co2_mixing_ratio_acc     + co2_mixing_ratio(ic)
      co2_mixing_ratio_C13_acc = co2_mixing_ratio_C13_acc + co2_mixing_ratio_C13(ic)
      co2_mixing_ratio_C14_acc = co2_mixing_ratio_C14_acc + co2_mixing_ratio_C14(ic)
      ! Carbon
      litterfallC_acc = litterfallC_acc + (veg_litterfall_mt(ix_leaf, ixC, ic) + & 
                              veg_litterfall_mt(ix_fine_root, ixC, ic) + veg_litterfall_mt(ix_coarse_root, ixC, ic) + &
                              veg_litterfall_mt(ix_sap_wood, ixC, ic) + veg_litterfall_mt(ix_heart_wood, ixC, ic) + &
                              veg_litterfall_mt(ix_labile, ixC, ic) + veg_litterfall_mt(ix_reserve, ixC, ic) + &
                              veg_litterfall_mt(ix_fruit, ixC, ic) + seed_bed_litterfall_mt(ixC, ic)) * &
                              1000000.0_wp * dtime_one/dtime
      litterfallC13_acc = litterfallC13_acc + (veg_litterfall_mt(ix_leaf, ixC13, ic) + & 
                              veg_litterfall_mt(ix_fine_root, ixC13, ic) + veg_litterfall_mt(ix_coarse_root, ixC13, ic) + &
                              veg_litterfall_mt(ix_sap_wood, ixC13, ic) + veg_litterfall_mt(ix_heart_wood, ixC13, ic) + &
                              veg_litterfall_mt(ix_labile, ixC13, ic) + veg_litterfall_mt(ix_reserve, ixC13, ic) + &
                              veg_litterfall_mt(ix_fruit, ixC13, ic) + seed_bed_litterfall_mt(ixC13, ic)) * &
                              1000000.0_wp * dtime_one/dtime
      litterfallC14_acc = litterfallC14_acc + (veg_litterfall_mt(ix_leaf, ixC14, ic) + & 
                              veg_litterfall_mt(ix_fine_root, ixC14, ic) + veg_litterfall_mt(ix_coarse_root, ixC14, ic) + &
                              veg_litterfall_mt(ix_sap_wood, ixC14, ic) + veg_litterfall_mt(ix_heart_wood, ixC14, ic) + &
                              veg_litterfall_mt(ix_labile, ixC14, ic) + veg_litterfall_mt(ix_reserve, ixC14, ic) + &
                              veg_litterfall_mt(ix_fruit, ixC14, ic) + seed_bed_litterfall_mt(ixC14, ic)) * &
                              1000000.0_wp * dtime_one/dtime
      litterfallN_acc = litterfallN_acc + (veg_litterfall_mt(ix_leaf, ixN, ic) + & 
                              veg_litterfall_mt(ix_fine_root, ixN, ic) + veg_litterfall_mt(ix_coarse_root, ixN, ic) + &
                              veg_litterfall_mt(ix_sap_wood, ixN, ic) + veg_litterfall_mt(ix_heart_wood, ixN, ic) + &
                              veg_litterfall_mt(ix_labile, ixN, ic) + veg_litterfall_mt(ix_reserve, ixN, ic) + &
                              veg_litterfall_mt(ix_fruit, ixN, ic) + seed_bed_litterfall_mt(ixN, ic)) * &
                              1000000.0_wp * dtime_one/dtime
      litterfallP_acc = litterfallP_acc + (veg_litterfall_mt(ix_leaf, ixP, ic) + & 
                              veg_litterfall_mt(ix_fine_root, ixP, ic) + veg_litterfall_mt(ix_coarse_root, ixP, ic) + &
                              veg_litterfall_mt(ix_sap_wood, ixP, ic) + veg_litterfall_mt(ix_heart_wood, ixP, ic) + &
                              veg_litterfall_mt(ix_labile, ixP, ic) + veg_litterfall_mt(ix_reserve, ixP, ic) + &
                              veg_litterfall_mt(ix_fruit, ixP, ic) + seed_bed_litterfall_mt(ixP, ic)) * &
                              1000000.0_wp * dtime_one/dtime
      het_respC_acc      = het_respC_acc      + SUM(het_respiration(ic,:) * soil_depth_sl(ic,:)) * dtime_one
      het_respC13_acc    = het_respC13_acc    + SUM(het_respiration_c13(ic,:) * soil_depth_sl(ic,:)) * dtime_one
      het_respC14_acc    = het_respC14_acc    + SUM(het_respiration_c14(ic,:) * soil_depth_sl(ic,:)) * dtime_one
      leachingCorg_acc   = leachingCorg_acc   + leaching_dom_carbon(ic)   * dtime_one
      leachingC13org_acc = leachingC13org_acc + leaching_dom_carbon13(ic) * dtime_one
      leachingC14org_acc = leachingC14org_acc + leaching_dom_carbon14(ic) * dtime_one
      lateral_loss_Corg_acc   = lateral_loss_Corg_acc   + SUM(lateral_loss_dom_carbon_sl(ic,:) &
        &                       * soil_depth_sl(ic,:)) * dtime_one 
      lateral_loss_C13org_acc = lateral_loss_C13org_acc + SUM(lateral_loss_dom_carbon13_sl(ic,:) &
        &                       * soil_depth_sl(ic,:)) * dtime_one
      lateral_loss_C14org_acc = lateral_loss_C14org_acc + SUM(lateral_loss_dom_carbon14_sl(ic,:) &
        &                       * soil_depth_sl(ic,:)) * dtime_one 
      ! Nitrogen
      n_deposition_acc     = n_deposition_acc  + (nhx_deposition(ic) + noy_deposition(ic)) * dtime_one      ! unit: mumol/m2/sec * dtime_one
      n_deposition_n15_acc = n_deposition_n15_acc  + (nhx_n15_deposition(ic) + noy_n15_deposition(ic)) * dtime_one      ! unit: mumol/m2/sec * dtime_one
      net_min_nh4_acc      = net_min_nh4_acc + SUM(net_mineralisation_nh4(ic,:) * soil_depth_sl(ic,:)) * dtime_one
      net_min_no3_acc      = net_min_no3_acc  + SUM(net_mineralisation_no3(ic,:) &
        &                    * soil_depth_sl(ic,:)) * dtime_one                     
      plant_uptake_nh4_acc = plant_uptake_nh4_acc + &
                           SUM(plant_uptake_nh4_sl(ic,:) * soil_depth_sl(ic,:)) * dtime_one 
      plant_uptake_no3_acc = plant_uptake_no3_acc + &
                           SUM(plant_uptake_no3_sl(ic,:) * soil_depth_sl(ic,:)) * dtime_one
      nitrification_acc    = nitrification_acc + & 
                            SUM((nitrification_no3(ic,:) + nitrification_noy(ic,:) + nitrification_n2o(ic,:)) * &
                            soil_depth_sl(ic,:)) * dtime_one
      denitrification_acc  = denitrification_acc + & 
                            SUM((denitrification_n2(ic,:) + denitrification_noy(ic,:) + denitrification_n2o(ic,:))*&
                            soil_depth_sl(ic,:)) * dtime_one 
      emission_noy_acc = emission_noy_acc + emission_noy(ic) * dtime_one
      emission_noy_n15_acc = emission_noy_n15_acc + emission_noy_n15(ic) * dtime_one
      emission_n2o_acc = emission_n2o_acc + emission_n2o(ic) * dtime_one
      emission_n2o_n15_acc = emission_n2o_n15_acc + emission_n2o_n15(ic) * dtime_one
      emission_n2_acc  = emission_n2_acc  + emission_n2(ic) * dtime_one
      emission_n2_n15_acc  = emission_n2_n15_acc  + emission_n2_n15(ic) * dtime_one
      leaching_nh4_acc = leaching_nh4_acc + leaching_nh4_solute(ic) * dtime_one
      leaching_nh4_n15_acc = leaching_nh4_n15_acc + leaching_nh4_n15_solute(ic) * dtime_one
      leaching_no3_acc = leaching_no3_acc + leaching_no3_solute(ic) * dtime_one
      leaching_no3_n15_acc = leaching_no3_n15_acc + leaching_no3_n15_solute(ic) * dtime_one
      leachingNorg_acc = leachingNorg_acc + leaching_dom_nitrogen(ic) * dtime_one
      lateral_loss_nh4_acc     = lateral_loss_nh4_acc + SUM(lateral_loss_nh4_solute_sl(ic,:) &
        &                        * soil_depth_sl(ic,:)) * dtime_one 
      lateral_loss_nh4_n15_acc = lateral_loss_nh4_n15_acc + SUM(lateral_loss_nh4_n15_solute_sl(ic,:) &
        &                        * soil_depth_sl(ic,:)) * dtime_one 
      lateral_loss_no3_acc     = lateral_loss_no3_acc + SUM(lateral_loss_no3_solute_sl(ic,:) &
        &                        * soil_depth_sl(ic,:)) * dtime_one 
      lateral_loss_no3_n15_acc = lateral_loss_no3_n15_acc + SUM(lateral_loss_no3_n15_solute_sl(ic,:) &
        &                        * soil_depth_sl(ic,:)) * dtime_one 
      lateral_loss_Norg_acc    = lateral_loss_Norg_acc + SUM(lateral_loss_dom_nitrogen_sl(ic,:) &
        &                        * soil_depth_sl(ic,:)) * dtime_one 
      asymb_n_fixation_acc =  asymb_n_fixation_acc + SUM(asymb_n_fixation(ic,:)*soil_depth_sl(ic,:)) * dtime_one
      asymb_n_fixation_n15_acc =  asymb_n_fixation_n15_acc + SUM(asymb_n_fixation_n15(ic,:)*soil_depth_sl(ic,:)) * &
                                  dtime_one
      fast_exchange_nh4_acc  = fast_exchange_nh4_acc + &
                                 SUM(fast_exchange_nh4(ic,:) * soil_depth_sl(ic,:)) * dtime_one
      myc_uptake_nh4_acc  = myc_uptake_nh4_acc  + SUM(mycorrhiza_uptake_nh4_sl(ic,:)  * soil_depth_sl(ic,:)) * dtime_one
      myc_uptake_no3_acc  = myc_uptake_no3_acc  + SUM(mycorrhiza_uptake_no3_sl(ic,:)  * soil_depth_sl(ic,:)) * dtime_one
      myc_uptake_norg_acc = myc_uptake_norg_acc + SUM(mycorrhiza_uptake_norg_sl(ic,:) * soil_depth_sl(ic,:)) * dtime_one
      ! Phosphorus
      p_deposition_acc       = p_deposition_acc + p_deposition(ic) * dtime_one
      leaching_po4_acc       = leaching_po4_acc + leaching_po4_solute(ic) * dtime_one           
      leachingPorg_acc       = leachingPorg_acc + leaching_dom_phosphorus(ic) * dtime_one
      lateral_loss_po4_acc   = lateral_loss_po4_acc + SUM(lateral_loss_po4_solute_sl(ic,:) &
        &                      * soil_depth_sl(ic,:)) * dtime_one 
      lateral_loss_Porg_acc  = lateral_loss_Porg_acc + SUM(lateral_loss_dom_phosphorus_sl(ic,:) &
        &                      * soil_depth_sl(ic,:)) * dtime_one 
      net_min_po4_acc        = net_min_po4_acc + &
                                SUM(net_mineralisation_po4(ic,:) * soil_depth_sl(ic,:)) * dtime_one                
      plant_uptake_po4_acc   = plant_uptake_po4_acc + &
                              SUM(plant_uptake_po4_sl(ic,:) * soil_depth_sl(ic,:)) * dtime_one
      weathering_po4_acc     = weathering_po4_acc + &
                                SUM(weathering_po4(ic,:) * soil_depth_sl(ic,:)) * dtime_one
      occlusion_po4_acc      = occlusion_po4_acc + &
                                SUM(occlusion_po4(ic,:) * soil_depth_sl(ic,:)) * dtime_one
      slow_exchange_po4_acc  = slow_exchange_po4_acc + &
                                 SUM(slow_exchange_po4(ic,:) * soil_depth_sl(ic,:)) * dtime_one
      fast_exchange_po4_acc  = fast_exchange_po4_acc + &
                                 SUM(fast_exchange_po4(ic,:) * soil_depth_sl(ic,:)) * dtime_one
      biochem_mineralisation_po4_acc = biochem_mineralisation_po4_acc + &
                                      SUM(biochem_mineralisation_po4(ic,:) * soil_depth_sl(ic,:)) * dtime_one
      ! mass balance checks in SB_
      sb_delta_pools_carbon_acc     = sb_delta_pools_carbon_acc     + &
                                        SUM(sb_delta_pools_mt(ixC, ic, :)*soil_depth_sl(ic,:))
      sb_delta_pools_nitrogen_acc   = sb_delta_pools_nitrogen_acc   + &
                                        SUM(sb_delta_pools_mt(ixN, ic, :)*soil_depth_sl(ic,:))
      sb_delta_pools_phosphorus_acc = sb_delta_pools_phosphorus_acc + &
                                        SUM(sb_delta_pools_mt(ixP, ic, :)*soil_depth_sl(ic,:))
      sb_delta_pools_carbon13_acc   = sb_delta_pools_carbon13_acc   + &
                                        SUM(sb_delta_pools_mt(ixC13, ic, :)*soil_depth_sl(ic,:))
      sb_delta_pools_carbon14_acc   = sb_delta_pools_carbon14_acc   + &
                                        SUM(sb_delta_pools_mt(ixC14, ic, :)*soil_depth_sl(ic,:))
      sb_delta_pools_nitrogen15_acc = sb_delta_pools_nitrogen15_acc + &
                                        SUM(sb_delta_pools_mt(ixN15, ic, :)*soil_depth_sl(ic,:))
      sb_total_flux_carbon_acc      = sb_total_flux_carbon_acc      + total_flux_carbon(ic)
      sb_total_flux_nitrogen_acc    = sb_total_flux_nitrogen_acc    + total_flux_nitrogen(ic)
      sb_total_flux_phosphorus_acc  = sb_total_flux_phosphorus_acc  + total_flux_phosphorus(ic)
      sb_total_flux_carbon13_acc    = sb_total_flux_carbon13_acc    + total_flux_carbon13(ic)
      sb_total_flux_carbon14_acc    = sb_total_flux_carbon14_acc    + total_flux_carbon14(ic)
      sb_total_flux_nitrogen15_acc  = sb_total_flux_nitrogen15_acc  + total_flux_nitrogen15(ic)
      ! accumulate layer-variables
      ! C 
      het_respiration_lay_acc(1:nsoil_w)     = het_respiration_lay_acc(1:nsoil_w) + &
                                                het_respiration(ic,:) * dtime_one
      het_respiration_c13_lay_acc(1:nsoil_w) = het_respiration_c13_lay_acc(1:nsoil_w) + &
                                                het_respiration_c13(ic,:) * dtime_one
      het_respiration_c14_lay_acc(1:nsoil_w) = het_respiration_c14_lay_acc(1:nsoil_w) + &
                                                het_respiration_c14(ic,:) * dtime_one
      transport_DOC_lay_acc(1:nsoil_w)       = transport_DOC_lay_acc(1:nsoil_w) + &
                                                sb_transport_mt(ix_dom, ixC, ic, :) * &
                                                1000000.0_wp * dtime_one/dtime
      hlp1_2d(:,:) = 0.0_wp
      DO is = 1,INT(num_sl_above_bedrock(ic))
        hlp1_2d(ic,is) = root_fraction_sl(ic,is) / soil_depth_sl(ic,is)
      END DO
      root_respiration_lay_acc(1:nsoil_w)    = root_respiration_lay_acc(1:nsoil_w) + &
                                                fmaint_rate_root(ic) * &
                                                (veg_pool_mt(ix_fine_root, ixN, ic) &
                                                + veg_pool_mt(ix_coarse_root, ixN, ic)) * &
                                                hlp1_2d(ic,:) * dtime_one
                                                ! root_fraction_sl(ic,:) / soil_depth_sl(ic,:) * dtime_one
      lateral_loss_Corg_lay_acc(1:nsoil_w)   = lateral_loss_Corg_lay_acc(1:nsoil_w) + &
                                                lateral_loss_dom_carbon_sl(ic,:) * &
                                                dtime_one
      ! N
      net_mineral_nh4_lay_acc(1:nsoil_w)     = net_mineral_nh4_lay_acc(1:nsoil_w) + &
                                                net_mineralisation_nh4(ic,:) * dtime_one
      plant_uptake_nh4_lay_acc(1:nsoil_w)    = plant_uptake_nh4_lay_acc(1:nsoil_w) + &
                                                plant_uptake_nh4_sl(ic,:) * dtime_one
      microbial_upt_nh4_lay_acc(1:nsoil_w)   = microbial_upt_nh4_lay_acc(1:nsoil_w) + &
                                                microbial_uptake_nh4_sl(ic,:) * dtime_one
      microbial_upt_no3_lay_acc(1:nsoil_w)   = microbial_upt_no3_lay_acc(1:nsoil_w) + &
                                                microbial_uptake_no3_sl(ic,:) * dtime_one
      myc_uptake_nh4_lay_acc(1:nsoil_w)      = myc_uptake_nh4_lay_acc(1:nsoil_w) + &
                                                mycorrhiza_uptake_nh4_sl(ic,:) * dtime_one
      transport_solute_nh4_acc(1:nsoil_w)    = transport_solute_nh4_acc(1:nsoil_w) + &
                                                transport_nh4_solute(ic,:) * dtime_one
      transport_nh4_assoc_acc(1:nsoil_w)     = transport_nh4_assoc_acc(1:nsoil_w) + &
                                                transport_nh4_assoc(ic,:) * dtime_one
      asymb_n_fixation_lay_acc(1:nsoil_w)    = asymb_n_fixation_lay_acc(1:nsoil_w) + &
                                                asymb_n_fixation(ic,:) * dtime_one
      net_mineral_no3_lay_acc(1:nsoil_w)     = net_mineral_no3_lay_acc(1:nsoil_w) + &
                                                net_mineralisation_no3(ic,:) * dtime_one
      plant_uptake_no3_lay_acc(1:nsoil_w)    = plant_uptake_no3_lay_acc(1:nsoil_w) + &
                                                plant_uptake_no3_sl(ic,:) * dtime_one
      myc_uptake_no3_lay_acc(1:nsoil_w)      = myc_uptake_no3_lay_acc(1:nsoil_w) + &
                                                mycorrhiza_uptake_no3_sl(ic,:) * dtime_one
      transport_solute_no3_acc(1:nsoil_w)    = transport_solute_no3_acc(1:nsoil_w) + &
                                                transport_no3_solute(ic,:) * dtime_one
      transport_DON_lay_acc(1:nsoil_w)       = transport_DON_lay_acc(1:nsoil_w) + &
                                                sb_transport_mt(ix_dom, ixN, ic, :) * &
                                                1000000.0_wp * dtime_one/dtime
      myc_uptake_norg_lay_acc(1:nsoil_w)     = myc_uptake_norg_lay_acc(1:nsoil_w) + &
                                                mycorrhiza_uptake_norg_sl(ic,:) * dtime_one
      f_adsorpt_nh4_lay_acc(1:nsoil_w)       = f_adsorpt_nh4_lay_acc(1:nsoil_w) + &
                                                fast_adsorpt_nh4(ic,:) * dtime_one
      nitrification_lay_acc(1:nsoil_w)       = nitrification_lay_acc(1:nsoil_w) + &
                                                (nitrification_no3(ic,:) + nitrification_noy(ic,:) + &
                                                 nitrification_n2o(ic,:)) * dtime_one
      denitrification_lay_acc(1:nsoil_w)     = denitrification_lay_acc(1:nsoil_w) + &
                                                (denitrification_n2(ic,:) + denitrification_noy(ic,:) + &
                                                 denitrification_n2o(ic,:)) * dtime_one
      lateral_loss_nh4_lay_acc(1:nsoil_w)    = lateral_loss_nh4_lay_acc(1:nsoil_w) + &
                                                lateral_loss_nh4_solute_sl(ic,:) * dtime_one
      lateral_loss_no3_lay_acc(1:nsoil_w)    = lateral_loss_no3_lay_acc(1:nsoil_w) + &
                                                lateral_loss_no3_solute_sl(ic,:) * dtime_one
      lateral_loss_Norg_lay_acc(1:nsoil_w)   = lateral_loss_Norg_lay_acc(1:nsoil_w) + &
                                                lateral_loss_dom_nitrogen_sl(ic,:) * dtime_one
      ! P
      f_exchange_po4_acc(1:nsoil_w)       = f_exchange_po4_acc(1:nsoil_w) + &
                                             fast_exchange_po4(ic,:) * dtime_one
      f_adsorpt_po4_acc(1:nsoil_w)        = f_adsorpt_po4_acc(1:nsoil_w) + &
                                             fast_adsorpt_po4(ic,:) * dtime_one
      net_min_po4_lay_acc(1:nsoil_w)      = net_min_po4_lay_acc(1:nsoil_w) + &
                                             net_mineralisation_po4(ic,:) * dtime_one
      weath_po4_acc(1:nsoil_w)            = weath_po4_acc(1:nsoil_w) + &
                                             weathering_po4(ic,:) * dtime_one
      ocl_po4_acc(1:nsoil_w)              = ocl_po4_acc(1:nsoil_w) + &
                                             occlusion_po4(ic,:) * dtime_one
      s_exchange_po4_acc(1:nsoil_w)       = s_exchange_po4_acc(1:nsoil_w) + &
                                             slow_exchange_po4(ic,:) * dtime_one
      biochem_min_po4_acc(1:nsoil_w)      = biochem_min_po4_acc(1:nsoil_w) + &
                                             biochem_mineralisation_po4(ic,:) * dtime_one
      plant_upt_po4_acc(1:nsoil_w)        = plant_upt_po4_acc(1:nsoil_w) + &
                                             plant_uptake_po4_sl(ic,:) * dtime_one
      microbial_upt_po4_acc(1:nsoil_w)    = microbial_upt_po4_acc(1:nsoil_w) + &
                                             microbial_uptake_po4_sl(ic,:) * dtime_one
      myc_upt_po4_acc(1:nsoil_w)          = myc_upt_po4_acc(1:nsoil_w) + &
                                             mycorrhiza_uptake_po4_sl(ic,:) * dtime_one
      transport_solute_po4_acc(1:nsoil_w)     = transport_solute_po4_acc(1:nsoil_w) + &
                                                 transport_po4_solute(ic,:) * dtime_one
      transport_assoc_fast_po4_acc(1:nsoil_w) = transport_assoc_fast_po4_acc(1:nsoil_w) + &
                                                 transport_po4_assoc_fast(ic,:) * dtime_one
      transport_assoc_slow_po4_acc(1:nsoil_w) = transport_assoc_slow_po4_acc(1:nsoil_w) + &
                                                 transport_po4_assoc_slow(ic,:) * dtime_one
      transport_occluded_po4_acc(1:nsoil_w)   = transport_occluded_po4_acc(1:nsoil_w) + &
                                                 transport_po4_occluded(ic,:) * dtime_one
      transport_primary_po4_acc(1:nsoil_w)    = transport_primary_po4_acc(1:nsoil_w) + &
                                                 transport_po4_primary(ic,:) * dtime_one
      transport_DOP_lay_acc(1:nsoil_w)        = transport_DOP_lay_acc(1:nsoil_w) + &
                                                 sb_transport_mt(ix_dom, ixP, ic, :) * &
                                                 1000000.0_wp * dtime_one/dtime
      lateral_loss_po4_lay_acc(1:nsoil_w)     = lateral_loss_po4_lay_acc(1:nsoil_w) + &
                                                 lateral_loss_po4_solute_sl(ic,:) * dtime_one
      lateral_loss_Porg_lay_acc(1:nsoil_w)    = lateral_loss_Porg_lay_acc(1:nsoil_w) + &
                                                 lateral_loss_dom_phosphorus_sl(ic,:) * dtime_one
    END IF  ! IF (.NOT. set_zero_flux_acc_var_local)

    !>2.0 write flux output
    !> 
    !>  and set set_zero_flux_acc_var_local to TRUE
    !>
    IF (write_flux_txt_output .AND. .NOT. set_zero_flux_acc_var_local) THEN
      IF (is_active_txt_mycflux_cnp) THEN
        ! name_mycflux_cnp
        write(lun_mycflux_cnp,'(7F24.12)') time_d, &
          &   myc_export_c_acc, myc_export_n_acc, myc_export_p_acc, &
          &   myc_export_c13_acc, myc_export_c14_acc, myc_export_n15_acc
      END IF
      IF (is_active_txt_myc_diagnostics) THEN
        ! name_myc_diagnostics
        write(lun_myc_diagnostics,'(3F24.12)') time_d, &
          &   myc_respiration_c_acc, f_cm_som_acc
      END IF
      IF (is_active_txt_matter_conservation_test) THEN
        ! name_matter_conservation_test
        !   output values at ~ 14 digits or so might vary because of precision issues
        !   only write 10 digits to avoid "false not-identical output warning"
        write(lun_matter_conservation_test,'(2F20.10)') time_d, &
          &   q_c_conservation_test(ic)
      END IF
      IF (is_active_txt_vegflux_n) THEN
        ! name_vegflux_n
        write(lun_vegflux_n,'(30F24.8)') &
          &   time_d, uptake_nh4_acc, uptake_no3_acc, recycling_n_acc, &
          &   res_use_n_acc, n_fixation_acc, &
          &   calc_delta_N15(n_fixation_acc, n_fixation_n15_acc), &
          &   calc_delta_N15(uptake_nh4_acc, uptake_nh4_n15_acc), &
          &   calc_delta_N15(uptake_no3_acc, uptake_no3_n15_acc), &
          &   growth_n_acc, &
          &   calc_delta_N15(veg_growth_sap_wood_n_acc, growth_swood_n15_acc), &
          &   veg_growth_leaf_n_acc, veg_growth_fine_root_n_acc, &
          &   veg_growth_coarse_root_n_acc, veg_growth_sap_wood_n_acc, &
          &   veg_growth_heart_wood_n_acc, veg_growth_labile_n_acc, &
          &   veg_growth_reserve_n_acc, veg_growth_fruit_n_acc, &
          &   veg_growth_seed_bed_n_acc, &
          &   veg_litterfall_leaf_n_acc, veg_litterfall_fine_root_n_acc, &
          &   veg_litterfall_coarse_root_n_acc, veg_litterfall_sap_wood_n_acc, &
          &   veg_litterfall_heart_wood_n_acc, veg_litterfall_labile_n_acc, &
          &   veg_litterfall_reserve_n_acc, veg_litterfall_fruit_n_acc, &
          &   veg_litterfall_seed_bed_n_acc
      END IF
      IF (is_active_txt_vegflux_p) THEN
        ! name_vegflux_p
        write(lun_vegflux_p,'(25F24.8)') &
          &   time_d, uptake_po4_acc, recycling_p_acc, res_use_p_acc, growth_p_acc, &
          &   veg_growth_leaf_p_acc, veg_growth_fine_root_p_acc, &
          &   veg_growth_coarse_root_p_acc, veg_growth_sap_wood_p_acc, &
          &   veg_growth_heart_wood_p_acc, veg_growth_labile_p_acc, &
          &   veg_growth_reserve_p_acc, veg_growth_fruit_p_acc, &
          &   veg_growth_seed_bed_p_acc, &
          &   veg_litterfall_leaf_p_acc, veg_litterfall_fine_root_p_acc, &
          &   veg_litterfall_coarse_root_p_acc, veg_litterfall_sap_wood_p_acc, &
          &   veg_litterfall_heart_wood_p_acc, veg_litterfall_labile_p_acc, &
          &   veg_litterfall_reserve_p_acc, veg_litterfall_fruit_p_acc, &
          &   veg_litterfall_seed_bed_p_acc
      END IF
      IF (is_active_txt_vegflux_c) THEN
        ! name_vegflux_c
        write(lun_vegflux_c,'(30F24.8)') &
          &   time_d, gpp_acc, growth_resp_acc, maint_resp_acc , &
          &   n_transform_resp_acc, n_fixation_resp_acc, npp_acc, &
          &   leaf_herb_resp_acc, fruit_herb_resp_acc, &
          &   res_use_c_acc, exudation_c_acc, growth_c_acc, &
          &   veg_growth_leaf_c_acc, veg_growth_fine_root_c_acc, &
          &   veg_growth_coarse_root_c_acc, veg_growth_sap_wood_c_acc, &
          &   veg_growth_heart_wood_c_acc, veg_growth_labile_c_acc, &
          &   veg_growth_reserve_c_acc, veg_growth_fruit_c_acc, &
          &   veg_growth_seed_bed_c_acc, &
          &   veg_litterfall_leaf_c_acc, veg_litterfall_fine_root_c_acc, &
          &   veg_litterfall_coarse_root_c_acc, veg_litterfall_sap_wood_c_acc, &
          &   veg_litterfall_heart_wood_c_acc, veg_litterfall_labile_c_acc, &
          &   veg_litterfall_reserve_c_acc, veg_litterfall_fruit_c_acc, &
          &   veg_litterfall_seed_bed_c_acc
      END IF
      IF (is_active_txt_vegflux_dc13dc14) THEN
        ! name_vegflux_dc13dc14
        dC13gpp = calc_delta_C13(gpp_acc,gpp_C13_acc)
        dC13rg  = calc_delta_C13(growth_resp_acc,growth_resp_C13_acc)
        dC13rm  = calc_delta_C13(maint_resp_acc,maint_resp_C13_acc) 
        dC13wgr = calc_delta_C13(veg_growth_sap_wood_c_acc, growth_swood_c13_acc)
        write(lun_vegflux_dc13dc14,'(9F18.6)') time_d, dC13gpp, dC13rg, dC13rm, dC13wgr  , &
          &   calc_delta_C14(gpp_acc,gpp_C14_acc,dC13gpp)                , &
          &   calc_delta_C14(growth_resp_acc,growth_resp_C14_acc,dC13rg) , &
          &   calc_delta_C14(maint_resp_acc,maint_resp_C14_acc,dC13rm)   , &
          &   calc_delta_C14(veg_growth_sap_wood_c_acc, growth_swood_c14_acc, growth_swood_c13_acc)
      END IF
      IF (is_active_txt_vegflux_h2o) THEN
        ! name_vegflux_h2o
        ! GC and GA are written out as values 'per output interval'
        IF (timestep_counter_flux_acc > 1.0_wp) THEN
          gc_acc = gc_acc / timestep_counter_flux_acc
          ga_acc = ga_acc / timestep_counter_flux_acc
        END IF
        write(lun_vegflux_h2o,'(12F18.8)') time_d, &
          &   precip_acc, evapopot_acc, evaporation_acc, interception_acc, &
          &   transpiration_acc, srf_runoff_acc, drainage_acc, gw_runoff_acc, gc_acc, ga_acc
      END IF
      IF (is_active_txt_vegflux_energy) THEN
        ! name_vegflux_energy
        ! albedo is written out as average of 'time of day with sunlight, i.e., `IF (alb_vis_lnd(ic) > 0.0_wp)`'
        IF (timestep_counter_albedo_acc > 1.0_wp) THEN
          alb_vis_acc = alb_vis_acc / timestep_counter_albedo_acc
          alb_nir_acc = alb_nir_acc / timestep_counter_albedo_acc
        END IF
        write(lun_vegflux_energy,'(9F18.8)') time_d, &
          &   rnet_acc, qh_acc, qle_acc, -qg_acc, alb_vis_acc, alb_nir_acc, albedo_noon(ic), fapar_noon(ic)
      END IF
      IF (is_active_txt_atm_forcing) THEN
        ! name_atm_forcing
        IF(timestep_counter_flux_acc > 1.0_wp) THEN
          t_air_acc                = t_air_acc                / timestep_counter_flux_acc
          wind_air_acc             = wind_air_acc             / timestep_counter_flux_acc
          swvis_srf_down_acc       = swvis_srf_down_acc       / timestep_counter_flux_acc
          swnir_srf_down_acc       = swnir_srf_down_acc       / timestep_counter_flux_acc
          lw_srf_down_acc          = lw_srf_down_acc          / timestep_counter_flux_acc
          vpd_acc                  = vpd_acc                  / timestep_counter_flux_acc
          co2_mixing_ratio_acc     = co2_mixing_ratio_acc     / timestep_counter_flux_acc
          co2_mixing_ratio_C13_acc = co2_mixing_ratio_C13_acc / timestep_counter_flux_acc
          co2_mixing_ratio_C14_acc = co2_mixing_ratio_C14_acc / timestep_counter_flux_acc
        END IF
        mean_sw_down = swvis_srf_down_acc + swnir_srf_down_acc
        dC13atm      = calc_delta_C13(co2_mixing_ratio_acc, co2_mixing_ratio_C13_acc) 
        write(lun_atm_forcing,'(10F20.8)') &
          &   time_d, precip_acc, t_air_acc, wind_air_acc, mean_sw_down, lw_srf_down_acc, &
          &   vpd_acc, co2_mixing_ratio_acc, dC13atm, &
          &   calc_delta_C14(co2_mixing_ratio_acc, co2_mixing_ratio_C14_acc, dC13atm)
      END IF
      IF (is_active_txt_soilflux_cdc13dc14) THEN
        ! name_soilflux_cdc13dc14
        dC13litterfall    = calc_delta_C13(litterfallC_acc,litterfallC13_acc)
        dC13het_resp      = calc_delta_C13(het_respC_acc,het_respC13_acc)
        dC13leaching      = calc_delta_C13(leachingCorg_acc,leachingC13org_acc)
        dC13lateral_loss  = calc_delta_C13(lateral_loss_Corg_acc, lateral_loss_C13org_acc) 
        write(lun_soilflux_cdc13dc14,'(15F20.12)') &
          &   time_d, litterfallC_acc, het_respC_acc, leachingCorg_acc, lateral_loss_Corg_acc, &
          &   dC13litterfall, dC13het_resp, dC13leaching, dC13lateral_loss, &
          &   calc_delta_C14(litterfallC_acc,litterfallC14_acc,dC13litterfall), &
          &   calc_delta_C14(het_respC_acc,het_respC14_acc,dC13het_resp), &
          &   calc_delta_C14(leachingCorg_acc,leachingC14org_acc,dC13leaching), &
          &   calc_delta_C14(lateral_loss_Corg_acc,lateral_loss_C14org_acc,dC13lateral_loss) 
      END IF
      IF (is_active_txt_soilflux_n) THEN
        ! name_soilflux_n
        write(lun_soilflux_n,'(25F20.12)') &
          &   time_d,litterfallN_acc,n_deposition_acc,ndep_mavg(1),asymb_n_fixation_acc,net_min_nh4_acc, net_min_no3_acc, &
          &   nitrification_acc,denitrification_acc,emission_noy_acc,emission_n2o_acc,emission_n2_acc,& 
          &   plant_uptake_nh4_acc,plant_uptake_no3_acc,leaching_nh4_acc,leaching_no3_acc,leachingNorg_acc, &
          &   lateral_loss_nh4_acc,lateral_loss_no3_acc,lateral_loss_Norg_acc, &
          &   fast_exchange_nh4_acc, myc_uptake_nh4_acc, myc_uptake_no3_acc, myc_uptake_norg_acc
      END IF
      IF (is_active_txt_soilflux_p) THEN
        ! name_soilflux_p
        write(lun_soilflux_p,'(22F20.12)') &
          &   time_d, litterfallP_acc, p_deposition_acc, leaching_po4_acc, leachingPorg_acc, &
          &   lateral_loss_po4_acc, lateral_loss_Porg_acc, &
          &   net_min_po4_acc, plant_uptake_po4_acc, weathering_po4_acc, occlusion_po4_acc, &
          &   slow_exchange_po4_acc, fast_exchange_po4_acc, biochem_mineralisation_po4_acc
      END IF
      IF (is_active_txt_soilflux_dn15) THEN
        ! name_soilflux_dn15
        write(lun_soilflux_dn15,'(20F18.8)') &
          &   time_d, &
          &   calc_delta_N15(n_deposition_acc,n_deposition_n15_acc), &
          &   calc_delta_N15(asymb_n_fixation_acc,asymb_n_fixation_n15_acc), &
          &   calc_delta_N15(emission_noy_acc,emission_noy_n15_acc), &
          &   calc_delta_N15(emission_n2o_acc,emission_n2o_n15_acc), &
          &   calc_delta_N15(emission_n2_acc,emission_n2_n15_acc), &
          &   calc_delta_N15(leaching_nh4_acc,leaching_nh4_n15_acc), &
          &   calc_delta_N15(leaching_no3_acc,leaching_no3_n15_acc), &
          &   calc_delta_N15(lateral_loss_nh4_acc,lateral_loss_nh4_n15_acc), &
          &   calc_delta_N15(lateral_loss_no3_acc,lateral_loss_no3_n15_acc) 
      END IF
      IF (is_active_txt_sb_mass_balance) THEN
        ! name_sb_mass_balance
        write(lun_sb_mass_balance,'(13F24.12)') &
          &   time_d, &
          &   sb_delta_pools_carbon_acc, sb_delta_pools_nitrogen_acc, sb_delta_pools_phosphorus_acc, &
          &   sb_delta_pools_carbon13_acc, sb_delta_pools_carbon14_acc, sb_delta_pools_nitrogen15_acc, &
          &   sb_total_flux_carbon_acc, sb_total_flux_nitrogen_acc, sb_total_flux_phosphorus_acc, &
          &   sb_total_flux_carbon13_acc, sb_total_flux_carbon14_acc, sb_total_flux_nitrogen15_acc
      END IF
      IF (is_active_txt_soilflux_c_sl) THEN
        ! name_soilflux_c_sl
        write(lun_soilflux_c_sl,'(121F20.10)') time_d, &
          &   het_respiration_lay_acc(1:nsoil_w), &
          &   calc_delta_C13(het_respiration_lay_acc(1:nsoil_w), het_respiration_c13_lay_acc(1:nsoil_w)), &
          &   calc_delta_C14(het_respiration_lay_acc(1:nsoil_w), het_respiration_c14_lay_acc(1:nsoil_w), &
          &                  calc_delta_C13(het_respiration_lay_acc(1:nsoil_w), het_respiration_c13_lay_acc(1:nsoil_w))), &
          &   transport_DOC_lay_acc(1:nsoil_w), &
          &   root_respiration_lay_acc(1:nsoil_w), &
          &   lateral_loss_Corg_lay_acc(1:nsoil_w)
      END IF
      IF (is_active_txt_soilflux_n_sl) THEN
        ! name_soilflux_n_sl
        write(lun_soilflux_n_sl,'(401F20.10)') time_d, &
          &   net_mineral_nh4_lay_acc(1:nsoil_w),   plant_uptake_nh4_lay_acc(1:nsoil_w), &
          &   microbial_upt_nh4_lay_acc(1:nsoil_w), myc_uptake_nh4_lay_acc(1:nsoil_w), &
          &   transport_solute_nh4_acc(1:nsoil_w),  transport_nh4_assoc_acc(1:nsoil_w), &
          &   asymb_n_fixation_lay_acc(1:nsoil_w), &
          &   net_mineral_no3_lay_acc(1:nsoil_w),   plant_uptake_no3_lay_acc(1:nsoil_w), &
          &   microbial_upt_no3_lay_acc(1:nsoil_w), myc_uptake_no3_lay_acc(1:nsoil_w), &
          &   transport_solute_no3_acc(1:nsoil_w),  transport_DON_lay_acc(1:nsoil_w), &
          &   myc_uptake_norg_lay_acc(1:nsoil_w), f_adsorpt_nh4_lay_acc(1:nsoil_w), &
          &   nitrification_lay_acc(1:nsoil_w), denitrification_lay_acc(1:nsoil_w), &
          &   lateral_loss_nh4_lay_acc(1:nsoil_w), &
          &   lateral_loss_no3_lay_acc(1:nsoil_w), &
          &   lateral_loss_Norg_lay_acc(1:nsoil_w) 
      END IF
      IF (is_active_txt_soilflux_p_sl) THEN
        ! name_soilflux_p_sl
        write(lun_soilflux_p_sl,'(361F20.10)') time_d, &
          &   weath_po4_acc(1:nsoil_w), ocl_po4_acc(1:nsoil_w), &
          &   f_exchange_po4_acc(1:nsoil_w), s_exchange_po4_acc(1:nsoil_w), &
          &   f_adsorpt_po4_acc(1:nsoil_w), &
          &   net_min_po4_lay_acc(1:nsoil_w), biochem_min_po4_acc(1:nsoil_w), &
          &   microbial_upt_po4_acc(1:nsoil_w), plant_upt_po4_acc(1:nsoil_w), &
          &   myc_upt_po4_acc(1:nsoil_w), &
          &   transport_solute_po4_acc(1:nsoil_w), transport_assoc_fast_po4_acc(1:nsoil_w), &
          &   transport_assoc_slow_po4_acc(1:nsoil_w), transport_occluded_po4_acc(1:nsoil_w), &
          &   transport_primary_po4_acc(1:nsoil_w), transport_DOP_lay_acc(1:nsoil_w), &
          &   lateral_loss_po4_lay_acc(1:nsoil_w), &
          &   lateral_loss_Porg_lay_acc(1:nsoil_w) 
      END IF

      !
      ! switch on the "set to zero of *_acc variables" (after the *_acc variable values were written out)
      !
      set_zero_flux_acc_var_local = .TRUE.
    END IF  ! IF (write_flux_txt_output .AND. .NOT. set_zero_flux_acc_var_local)

    !>3.0 set all *_acc variables zero
    !>       
    !>  either if:
    !>    a) flux-output variables were written out, and thus need to be zero'ed (set_zero_flux_acc_var_local is set to TRUE above)
    !>    b) this routine is called only to zero the *_acc flux-output variables (set_zero_flux_acc_variables is TRUE)
    !>         because (the output interval is changed and) new output files are created deleting potential prev. accumulated values
    !>
     IF (set_zero_flux_acc_var_local) THEN
      timestep_counter_flux_acc          = 0.0_wp
      timestep_counter_albedo_acc        = 0.0_wp
      ! vegetation
      gpp_acc                            = 0.0_wp
      growth_resp_acc                    = 0.0_wp
      maint_resp_acc                     = 0.0_wp
      n_transform_resp_acc               = 0.0_wp
      n_fixation_resp_acc                = 0.0_wp
      exudation_c_acc                    = 0.0_wp
      npp_acc                            = 0.0_wp
      leaf_herb_resp_acc                 = 0.0_wp
      fruit_herb_resp_acc                = 0.0_wp
      res_use_c_acc                      = 0.0_wp
      growth_c_acc                       = 0.0_wp
      growth_swood_c13_acc               = 0.0_wp
      growth_swood_c14_acc               = 0.0_wp
      veg_growth_leaf_c_acc              = 0.0_wp
      veg_growth_fine_root_c_acc         = 0.0_wp
      veg_growth_coarse_root_c_acc       = 0.0_wp
      veg_growth_sap_wood_c_acc          = 0.0_wp
      veg_growth_heart_wood_c_acc        = 0.0_wp
      veg_growth_labile_c_acc            = 0.0_wp
      veg_growth_reserve_c_acc           = 0.0_wp
      veg_growth_fruit_c_acc             = 0.0_wp
      veg_growth_seed_bed_c_acc          = 0.0_wp
      veg_growth_leaf_n_acc              = 0.0_wp
      veg_growth_fine_root_n_acc         = 0.0_wp
      veg_growth_coarse_root_n_acc       = 0.0_wp
      veg_growth_sap_wood_n_acc          = 0.0_wp
      veg_growth_heart_wood_n_acc        = 0.0_wp
      veg_growth_labile_n_acc            = 0.0_wp
      veg_growth_reserve_n_acc           = 0.0_wp
      veg_growth_fruit_n_acc             = 0.0_wp
      veg_growth_seed_bed_n_acc          = 0.0_wp
      veg_growth_leaf_p_acc              = 0.0_wp
      veg_growth_fine_root_p_acc         = 0.0_wp
      veg_growth_coarse_root_p_acc       = 0.0_wp
      veg_growth_sap_wood_p_acc          = 0.0_wp
      veg_growth_heart_wood_p_acc        = 0.0_wp
      veg_growth_labile_p_acc            = 0.0_wp
      veg_growth_reserve_p_acc           = 0.0_wp
      veg_growth_fruit_p_acc             = 0.0_wp
      veg_growth_seed_bed_p_acc          = 0.0_wp
      veg_litterfall_leaf_c_acc          = 0.0_wp
      veg_litterfall_fine_root_c_acc     = 0.0_wp
      veg_litterfall_coarse_root_c_acc   = 0.0_wp
      veg_litterfall_sap_wood_c_acc      = 0.0_wp
      veg_litterfall_heart_wood_c_acc    = 0.0_wp
      veg_litterfall_labile_c_acc        = 0.0_wp
      veg_litterfall_reserve_c_acc       = 0.0_wp
      veg_litterfall_fruit_c_acc         = 0.0_wp
      veg_litterfall_seed_bed_c_acc      = 0.0_wp
      veg_litterfall_leaf_n_acc          = 0.0_wp
      veg_litterfall_fine_root_n_acc     = 0.0_wp
      veg_litterfall_coarse_root_n_acc   = 0.0_wp
      veg_litterfall_sap_wood_n_acc      = 0.0_wp
      veg_litterfall_heart_wood_n_acc    = 0.0_wp
      veg_litterfall_labile_n_acc        = 0.0_wp
      veg_litterfall_reserve_n_acc       = 0.0_wp
      veg_litterfall_fruit_n_acc         = 0.0_wp
      veg_litterfall_seed_bed_n_acc      = 0.0_wp
      veg_litterfall_leaf_p_acc          = 0.0_wp
      veg_litterfall_fine_root_p_acc     = 0.0_wp
      veg_litterfall_coarse_root_p_acc   = 0.0_wp
      veg_litterfall_sap_wood_p_acc      = 0.0_wp
      veg_litterfall_heart_wood_p_acc    = 0.0_wp
      veg_litterfall_labile_p_acc        = 0.0_wp
      veg_litterfall_reserve_p_acc       = 0.0_wp
      veg_litterfall_fruit_p_acc         = 0.0_wp
      veg_litterfall_seed_bed_p_acc      = 0.0_wp             
      gpp_C13_acc                        = 0.0_wp
      growth_resp_C13_acc                = 0.0_wp
      maint_resp_C13_acc                 = 0.0_wp
      gpp_C14_acc                        = 0.0_wp
      growth_resp_C14_acc                = 0.0_wp
      maint_resp_C14_acc                 = 0.0_wp
      recycling_n_acc                    = 0.0_wp
      res_use_n_acc                      = 0.0_wp
      growth_n_acc                       = 0.0_wp
      growth_swood_n15_acc               = 0.0_wp
      n_fixation_acc                     = 0.0_wp
      n_fixation_n15_acc                 = 0.0_wp
      uptake_po4_acc                     = 0.0_wp
      recycling_p_acc                    = 0.0_wp
      res_use_p_acc                      = 0.0_wp
      growth_p_acc                       = 0.0_wp
      ! various
      precip_acc                 = 0.0_wp
      evaporation_acc            = 0.0_wp
      interception_acc           = 0.0_wp
      evapopot_acc               = 0.0_wp
      transpiration_acc          = 0.0_wp
      drainage_acc               = 0.0_wp
      gw_runoff_acc              = 0.0_wp
      srf_runoff_acc             = 0.0_wp
      gc_acc                     = 0.0_wp
      ga_acc                     = 0.0_wp
      rnet_acc                   = 0.0_wp
      qh_acc                     = 0.0_wp
      qle_acc                    = 0.0_wp
      qg_acc                     = 0.0_wp
      alb_vis_acc                = 0.0_wp
      alb_nir_acc                = 0.0_wp
      uptake_no3_acc             = 0.0_wp
      uptake_no3_n15_acc         = 0.0_wp
      uptake_nh4_acc             = 0.0_wp
      uptake_nh4_n15_acc         = 0.0_wp
      myc_export_c_acc           = 0.0_wp
      myc_export_n_acc           = 0.0_wp
      myc_export_p_acc           = 0.0_wp
      myc_export_c13_acc         = 0.0_wp
      myc_export_c14_acc         = 0.0_wp
      myc_export_n15_acc         = 0.0_wp
      myc_respiration_c_acc      = 0.0_wp
      f_cm_som_acc               = 0.0_wp
      ! matter conservation test
      c_cons_test_acc            = 0.0_wp
      ! atmospheric forcing (precip_acc is already considered above, see name_vegflux_h2o)
      t_air_acc                  = 0.0_wp
      wind_air_acc               = 0.0_wp
      swvis_srf_down_acc         = 0.0_wp
      swnir_srf_down_acc         = 0.0_wp
      lw_srf_down_acc            = 0.0_wp
      vpd_acc                    = 0.0_wp
      co2_mixing_ratio_acc       = 0.0_wp
      co2_mixing_ratio_C13_acc   = 0.0_wp
      co2_mixing_ratio_C14_acc   = 0.0_wp
      ! land & soil model
      litterfallC_acc                   = 0.0_wp
      het_respC_acc                     = 0.0_wp
      leachingCorg_acc                  = 0.0_wp
      lateral_loss_Corg_acc             = 0.0_wp
      litterfallC13_acc                 = 0.0_wp
      het_respC13_acc                   = 0.0_wp
      leachingC13org_acc                = 0.0_wp
      lateral_loss_C13org_acc           = 0.0_wp 
      litterfallC14_acc                 = 0.0_wp
      het_respC14_acc                   = 0.0_wp
      leachingC14org_acc                = 0.0_wp
      lateral_loss_C14org_acc           = 0.0_wp
      nitrification_acc                 = 0.0_wp
      denitrification_acc               = 0.0_wp
      litterfallN_acc                   = 0.0_wp
      net_min_nh4_acc                   = 0.0_wp
      net_min_no3_acc                   = 0.0_wp
      plant_uptake_nh4_acc              = 0.0_wp
      plant_uptake_no3_acc              = 0.0_wp
      leachingNorg_acc                  = 0.0_wp
      lateral_loss_Norg_acc             = 0.0_wp 
      myc_uptake_nh4_acc                = 0.0_wp
      myc_uptake_no3_acc                = 0.0_wp
      myc_uptake_norg_acc               = 0.0_wp
      fast_exchange_nh4_acc             = 0.0_wp
      emission_noy_acc                  = 0.0_wp
      emission_noy_n15_acc              = 0.0_wp
      emission_n2o_acc                  = 0.0_wp
      emission_n2o_n15_acc              = 0.0_wp
      emission_n2_acc                   = 0.0_wp
      emission_n2_n15_acc               = 0.0_wp
      leaching_nh4_acc                  = 0.0_wp
      leaching_nh4_n15_acc              = 0.0_wp
      leaching_no3_acc                  = 0.0_wp
      leaching_no3_n15_acc              = 0.0_wp
      lateral_loss_nh4_acc              = 0.0_wp 
      lateral_loss_nh4_n15_acc          = 0.0_wp 
      lateral_loss_no3_acc              = 0.0_wp 
      lateral_loss_no3_n15_acc          = 0.0_wp 
      asymb_n_fixation_acc              = 0.0_wp
      asymb_n_fixation_n15_acc          = 0.0_wp
      litterfallP_acc                   = 0.0_wp
      leachingPorg_acc                  = 0.0_wp
      lateral_loss_Porg_acc             = 0.0_wp 
      p_deposition_acc                  = 0.0_wp
      n_deposition_acc                  = 0.0_wp
      n_deposition_n15_acc              = 0.0_wp
      leaching_po4_acc                  = 0.0_wp
      lateral_loss_po4_acc              = 0.0_wp 
      net_min_po4_acc                   = 0.0_wp
      plant_uptake_po4_acc              = 0.0_wp
      weathering_po4_acc                = 0.0_wp
      occlusion_po4_acc                 = 0.0_wp
      slow_exchange_po4_acc             = 0.0_wp
      fast_exchange_po4_acc             = 0.0_wp
      biochem_mineralisation_po4_acc    = 0.0_wp
      sb_delta_pools_carbon_acc         = 0.0_wp
      sb_delta_pools_nitrogen_acc       = 0.0_wp
      sb_delta_pools_phosphorus_acc     = 0.0_wp
      sb_delta_pools_carbon13_acc       = 0.0_wp
      sb_delta_pools_carbon14_acc       = 0.0_wp
      sb_delta_pools_nitrogen15_acc     = 0.0_wp
      sb_total_flux_carbon_acc          = 0.0_wp
      sb_total_flux_nitrogen_acc        = 0.0_wp
      sb_total_flux_phosphorus_acc      = 0.0_wp
      sb_total_flux_carbon13_acc        = 0.0_wp
      sb_total_flux_carbon14_acc        = 0.0_wp
      sb_total_flux_nitrogen15_acc      = 0.0_wp
      ! simple_1d & jsm soil models
      ! C
      het_respiration_lay_acc(1:nsoil_w)         = 0.0_wp
      het_respiration_c13_lay_acc(1:nsoil_w)     = 0.0_wp
      het_respiration_c14_lay_acc(1:nsoil_w)     = 0.0_wp
      transport_DOC_lay_acc(1:nsoil_w)           = 0.0_wp
      root_respiration_lay_acc(1:nsoil_w)        = 0.0_wp
      lateral_loss_Corg_lay_acc(1:nsoil_w)       = 0.0_wp 
      ! N
      net_mineral_nh4_lay_acc(1:nsoil_w)         = 0.0_wp
      plant_uptake_nh4_lay_acc(1:nsoil_w)        = 0.0_wp
      microbial_upt_nh4_lay_acc(1:nsoil_w)       = 0.0_wp
      microbial_upt_no3_lay_acc(1:nsoil_w)       = 0.0_wp
      myc_uptake_nh4_lay_acc(1:nsoil_w)          = 0.0_wp
      asymb_n_fixation_lay_acc(1:nsoil_w)        = 0.0_wp
      net_mineral_no3_lay_acc(1:nsoil_w)         = 0.0_wp
      plant_uptake_no3_lay_acc(1:nsoil_w)        = 0.0_wp
      myc_uptake_no3_lay_acc(1:nsoil_w)          = 0.0_wp
      transport_solute_nh4_acc(1:nsoil_w)        = 0.0_wp
      transport_nh4_assoc_acc(1:nsoil_w)         = 0.0_wp
      transport_solute_no3_acc(1:nsoil_w)        = 0.0_wp
      transport_DON_lay_acc(1:nsoil_w)           = 0.0_wp
      myc_uptake_norg_lay_acc(1:nsoil_w)         = 0.0_wp
      f_adsorpt_nh4_lay_acc(1:nsoil_w)           = 0.0_wp
      nitrification_lay_acc(1:nsoil_w)           = 0.0_wp
      denitrification_lay_acc(1:nsoil_w)         = 0.0_wp
      lateral_loss_nh4_lay_acc(1:nsoil_w)        = 0.0_wp 
      lateral_loss_no3_lay_acc(1:nsoil_w)        = 0.0_wp 
      lateral_loss_Norg_lay_acc(1:nsoil_w)       = 0.0_wp 
      ! P
      f_exchange_po4_acc(1:nsoil_w)              = 0.0_wp
      f_adsorpt_po4_acc(1:nsoil_w)               = 0.0_wp
      weath_po4_acc(1:nsoil_w)                   = 0.0_wp
      ocl_po4_acc(1:nsoil_w)                     = 0.0_wp
      s_exchange_po4_acc(1:nsoil_w)              = 0.0_wp
      net_min_po4_lay_acc(1:nsoil_w)             = 0.0_wp
      biochem_min_po4_acc(1:nsoil_w)             = 0.0_wp
      microbial_upt_po4_acc(1:nsoil_w)           = 0.0_wp
      plant_upt_po4_acc(1:nsoil_w)               = 0.0_wp 
      myc_upt_po4_acc(1:nsoil_w)                 = 0.0_wp  
      transport_solute_po4_acc(1:nsoil_w)        = 0.0_wp
      transport_assoc_fast_po4_acc(1:nsoil_w)    = 0.0_wp
      transport_assoc_slow_po4_acc(1:nsoil_w)    = 0.0_wp
      transport_occluded_po4_acc(1:nsoil_w)      = 0.0_wp
      transport_primary_po4_acc(1:nsoil_w)       = 0.0_wp
      transport_DOP_lay_acc(1:nsoil_w)           = 0.0_wp      
      lateral_loss_po4_lay_acc(1:nsoil_w)        = 0.0_wp 
      lateral_loss_Porg_lay_acc(1:nsoil_w)       = 0.0_wp 
    END IF    ! IF (set_zero_flux_acc_var_local)
  END SUBROUTINE write_output_flux


  !-----------------------------------------------------------------------------------------------------!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! section III Soil forcing                                                                          !!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !-----------------------------------------------------------------------------------------------------!


  ! ======================================================================================================= !
  !>Soil forcing file: open files and write column headers
  !>
  !>  open output files, define file names for pools, write column header names
  !>
  SUBROUTINE open_output_files_soil_forcing(tile, iblk, ics, ice, nc)

    USE mo_jsb_class,             ONLY: Get_model
    USE mo_jsb_tile_class,        ONLY: t_jsb_tile_abstract
    USE mo_jsb_model_class,       ONLY: t_jsb_model
    USE mo_jsb_process_class,     ONLY: VEG_, SB_
    
    ! Use of process configurations (t_PROC_config)
    dsl4jsb_Use_config(VEG_)
    dsl4jsb_Use_config(SB_)
    ! ---------------------------
    ! 0.1 InOut
    CLASS(t_jsb_tile_abstract), INTENT(in) :: tile         !< one tile with data structure for one lct
    INTEGER,                    INTENT(in) :: iblk, &      !< iBlock, start point of chunk, end point of chunk, number of points in the chunk
                                              ics, &
                                              ice, &
                                              nc
    ! ---------------------------
    ! 0.2 Local
    TYPE(t_jsb_model),POINTER             :: model
    INTEGER        :: lun           ! logical unit number
    CHARACTER(100) :: fname         ! filename
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':open_output_files_soil_forcing'
    ! ---------------------------
    ! 0.3 Declare Memory
    ! Declare process configuration and memory Pointers
    dsl4jsb_Def_config(VEG_)
    dsl4jsb_Def_config(SB_)
    ! ---------------------------
    ! 0.5 Get Memory
    model  => Get_model(tile%owner_model_id)
    ! Get process config
    dsl4jsb_Get_config(VEG_)
    dsl4jsb_Get_config(SB_)




    !! 
    !! output interval will be always timestep !
    !! 


    !> 1.0 define open() statements for file output
    !! 
    ! set initial lun to not conflict with the above pool & flux output
    lun = lun_ascii_soil_forc_start
    fname = 'soilphysics'//'_timestep.txt'
    OPEN (UNIT=lun, FILE=fname, ACTION="write", STATUS="unknown")
    lun=lun+1
    fname = 'veglitterfall'//'_timestep.txt'
    OPEN (UNIT=lun, FILE=fname, ACTION="write", STATUS="unknown")
    lun=lun+1

    !> 1.1  test if lun for soil forcing is too big
    !! 
    IF (lun > lun_ascii_soil_forc_max) THEN
      WRITE(message_text,'(a)') 'lun of soil forcing output files is too large'
      CALL finish(TRIM(routine),message_text)
    END IF

    !> 2.0 write column headers into output files
    !! 
    CALL init_output_files_soil_forcing()

  END SUBROUTINE open_output_files_soil_forcing


  ! ======================================================================================================= !
  !>Soil forcing files: write column headers of output files
  !>
  !>  complete litterfall output (used as input for the stand alone soil model)
  !>  and some soilphysics variables
  !>
  SUBROUTINE init_output_files_soil_forcing()

    USE mo_jsb_grid_class,        ONLY: t_jsb_vgrid
    USE mo_jsb_grid,              ONLY: Get_vgrid
    ! ---------------------------
    ! 0.1 InOut
    !..
    ! ---------------------------
    ! 0.2 Local
    TYPE(t_jsb_vgrid), POINTER    :: vgrid_soil_w                ! Vertical grid
    INTEGER                       :: nsoil_w                     ! number of soil layers as used/defined by the SB_ process
    INTEGER                       :: lun, isoil
    CHARACTER(2)                  :: layerstr
    CHARACTER(450)                :: p1str,p2str,p3str,p4str,p5str,p6str,p7str,p8str,p9str,p10str
    ! ---------------------------
    ! 0.5 Get Memory
    vgrid_soil_w => Get_vgrid('soil_depth_water')
    nsoil_w      =  vgrid_soil_w%n_levels



    !> 1.0 write headers
    !!
    ! set initial lun to not conflict with the above pool & flux output
    lun = lun_ascii_soil_forc_start
    
    ! soilphysics
    p1str=''
    p2str=''
    p3str=''
    p4str=''
    DO isoil = 1,nsoil_w
      IF (isoil < 10) THEN 
        WRITE(layerstr,'(I1)') isoil
      ELSE
        WRITE(layerstr,'(I2)') isoil
      END IF
      WRITE(p1str,'(a)') TRIM(p1str)//' TemperatureSoil_L'//TRIM(layerstr)
      WRITE(p2str,'(a)') TRIM(p2str)//' WaterContentSoil_L'//TRIM(layerstr)
      WRITE(p3str,'(a)') TRIM(p3str)//' PercolationSoil_L'//TRIM(layerstr)
      WRITE(p4str,'(a)') TRIM(p4str)//' RootFractionSoil_L'//TRIM(layerstr)
    END DO
    WRITE(lun,'(8a)') 'Time ', &
                      TRIM(p1str), TRIM(p2str), TRIM(p3str), TRIM(p4str), &
                      ' VegPoolFineRootCarbon', &
                      ' nhx_deposition noy_deposition nhx_dn15_deposition noy_dn15_deposition p_deposition', &
                      ' demandUptakeN demandUptakeP'
    lun=lun+1
    ! veglitterfall
    WRITE(lun,'(10a)')  'Time ', &
                        'leafC leafN leafP leafC13 leafC14 leafN15 ' , &
                        'fine_rootC fine_rootN fine_rootP fine_rootC13 fine_rootC14 fine_rootN15 ' , &
                        'coarse_rootC coarse_rootN coarse_rootP coarse_rootC13 coarse_rootC14 coarse_rootN15 ' , &
                        'sap_woodC sap_woodN sap_woodP sap_woodC13 sap_woodC14 sap_woodN15 ' , &
                        'heart_woodC heart_woodN heart_woodP heart_woodC13 heart_woodC14 heart_woodN15 ' , &
                        'labileC labileN labileP labileC13 labileC14 labileN15 ' , &
                        'reserveC reserveN reserveP reserveC13 reserveC14 reserveN15 ' , &
                        'fruitC fruitN fruitP fruitC13 fruitC14 fruitN15 ' , &
                        'seed_bedC seed_bedN seed_bedP seed_bedC13 seed_bedC14 seed_bedN15'
    lun=lun+1

  END SUBROUTINE init_output_files_soil_forcing

  ! ======================================================================================================= !
  !>Subroutine to generate the forcing file for the stand alone soil model from 
  !>  litterfall data and SB rate modifiers 
  !>
  SUBROUTINE write_soil_forcing(tile, options, time_d)

    USE mo_jsb_class,             ONLY: Get_model
    USE mo_jsb_tile_class,        ONLY: t_jsb_tile_abstract
    USE mo_jsb_task_class,        ONLY: t_jsb_task_options
    USE mo_jsb_model_class,       ONLY: t_jsb_model
    USE mo_jsb_process_class,     ONLY: A2L_, HYDRO_, SSE_, VEG_
    USE mo_isotope_util,          ONLY: calc_delta_C13,calc_delta_C14,calc_delta_N15
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Use_memory(A2L_)
    dsl4jsb_Use_memory(VEG_)
    dsl4jsb_Use_memory(HYDRO_)
    dsl4jsb_Use_memory(SSE_)
    ! ----------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_tile_abstract), INTENT(inout)     :: tile         !< one tile with data structure for one lct
    TYPE(t_jsb_task_options),   INTENT(in)        :: options      !< model options
    REAL(wp),                   INTENT(in)        :: time_d       !< time in days
    ! ----------------------------------------------------------------------------------------------------- !
    TYPE(t_jsb_model),      POINTER :: model          !< the model
    TYPE(t_lnd_bgcm_store), POINTER :: bgcm_store     !< the bgcm store of this tile
    TYPE(t_lctlib_element), POINTER :: lctlib         !< land-cover-type library - parameter across pft's
    REAL(wp)                        :: dC13_leaf,dC13_fine_root,dC13_coarse_root,&
      &                                dC13_sap_wood,dC13_heart_wood,dC13_labile, &
                                       dC13_reserve,dC13_fruit,dC13_seed_bed
    REAL(wp)                        :: DC14_leaf,DC14_fine_root,DC14_coarse_root,&
      &                                DC14_sap_wood,DC14_heart_wood,DC14_labile, &
                                       DC14_reserve,DC14_fruit,DC14_seed_bed
    REAL(wp)                        :: dN15_leaf,dN15_fine_root,dN15_coarse_root,& 
      &                                dN15_sap_wood,dN15_heart_wood,dN15_labile, &
                                       dN15_reserve,dN15_fruit,dN15_seed_bed
    INTEGER                         :: lun          ! logical unit number
    INTEGER                         :: iblk, ics, ice, nc
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':write_soil_forcing'
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Def_mt2L2D :: veg_pool_mt
    dsl4jsb_Def_mt2L2D :: veg_litterfall_mt
    dsl4jsb_Def_mt1L2D :: seed_bed_litterfall_mt
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Def_memory(A2L_)
    dsl4jsb_Def_memory(VEG_)
    dsl4jsb_Def_memory(HYDRO_)
    dsl4jsb_Def_memory(SSE_)
    ! ----------------------------------------------------------------------------------------------------- !
    ! A2L_
    dsl4jsb_Real2D_onChunk      :: nhx_deposition 
    dsl4jsb_Real2D_onChunk      :: noy_deposition 
    dsl4jsb_Real2D_onChunk      :: nhx_n15_deposition 
    dsl4jsb_Real2D_onChunk      :: noy_n15_deposition 
    dsl4jsb_Real2D_onChunk      :: p_deposition 
    ! VEG_ 2D
    dsl4jsb_Real2D_onChunk      :: demand_uptake_n_tuptake_mavg
    dsl4jsb_Real2D_onChunk      :: demand_uptake_p_tuptake_mavg
    ! VEG_ 3D
    dsl4jsb_Real3D_onChunk      :: root_fraction_sl
    ! HYDRO_ 3D
    dsl4jsb_Real3D_onChunk      :: frac_wtr_transp_down_sl
    dsl4jsb_Real3D_onChunk      :: wtr_soil_sl
    ! SSE_ 3D
    dsl4jsb_Real3D_onChunk      :: t_soil_sl
    ! ----------------------------------------------------------------------------------------------------- !
    iblk      = options%iblk
    ics       = options%ics
    ice       = options%ice
    nc        = options%nc
    ! ----------------------------------------------------------------------------------------------------- !
    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')
    ! ----------------------------------------------------------------------------------------------------- !
    model  => Get_model(tile%owner_model_id)
    lctlib => model%lctlib(tile%lcts(1)%lib_id) 
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Get_memory(A2L_)
    dsl4jsb_Get_memory(VEG_)
    dsl4jsb_Get_memory(HYDRO_)
    dsl4jsb_Get_memory(SSE_)
    ! ----------------------------------------------------------------------------------------------------- !
    bgcm_store => tile%bgcm_store
    dsl4jsb_Get_mt2L2D(VEG_BGCM_POOL_ID, veg_pool_mt)
    dsl4jsb_Get_mt2L2D(VEG_BGCM_LITTERFALL_ID, veg_litterfall_mt)
    dsl4jsb_Get_mt1L2D(VEG_BGCM_SEED_BED_LITTERFALL_ID, seed_bed_litterfall_mt)
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Get_var2D_onChunk(A2L_,  nhx_deposition)
    dsl4jsb_Get_var2D_onChunk(A2L_,  noy_deposition)
    dsl4jsb_Get_var2D_onChunk(A2L_,  nhx_n15_deposition)
    dsl4jsb_Get_var2D_onChunk(A2L_,  noy_n15_deposition)
    dsl4jsb_Get_var2D_onChunk(A2L_,  p_deposition)
    ! ---------------------------
    dsl4jsb_Get_var2D_onChunk(VEG_,  demand_uptake_n_tuptake_mavg)
    dsl4jsb_Get_var2D_onChunk(VEG_,  demand_uptake_p_tuptake_mavg)
    dsl4jsb_Get_var3D_onChunk(VEG_,  root_fraction_sl)
    ! HYDRO_ 3D
    dsl4jsb_Get_var3D_onChunk(HYDRO_, frac_wtr_transp_down_sl)
    dsl4jsb_Get_var3D_onChunk(HYDRO_, wtr_soil_sl)
    ! SSE_ 3D
    dsl4jsb_Get_var3D_onChunk(SSE_, t_soil_sl)                      ! in
    ! ----------------------------------------------------------------------------------------------------- !

    !>1.0 write output
    !>
    ! set initial lun to not conflict with the above pool & flux output
    lun = lun_ascii_soil_forc_start
    ! soilphysics
    write(lun,'(90F18.8)') time_d                                                       , &
                            ! variables with one value per soil layer (20 sl max)
                            t_soil_sl(ic,:)                                             , &
                            wtr_soil_sl(ic,:)                                           , &
                            frac_wtr_transp_down_sl(ic,:)                               , &
                            root_fraction_sl(ic,:)                                      , &
                            ! variables with one value only
                            veg_pool_mt(ix_fine_root, ixC, ic)                          , &
                            nhx_deposition(ic)                                          , &
                            noy_deposition(ic)                                          , &
                            calc_delta_N15(nhx_deposition(ic),nhx_n15_deposition(ic))   , &
                            calc_delta_N15(noy_deposition(ic),noy_n15_deposition(ic))   , &
                            p_deposition(ic)                                            , &
                            demand_uptake_n_tuptake_mavg(ic)                            , &
                            demand_uptake_p_tuptake_mavg(ic) 
                            
    lun=lun+1
    !> convert 13C and 14C in delta notation
    !! 
    dC13_leaf        = calc_delta_C13(veg_litterfall_mt(ix_leaf, ixC, ic),veg_litterfall_mt(ix_leaf, ixC13, ic))
    dC13_fine_root   = calc_delta_C13(veg_litterfall_mt(ix_fine_root, ixC, ic),veg_litterfall_mt(ix_fine_root, ixC13, ic))
    dC13_coarse_root = calc_delta_C13(veg_litterfall_mt(ix_coarse_root, ixC, ic),veg_litterfall_mt(ix_coarse_root, ixC13, ic))
    dC13_sap_wood    = calc_delta_C13(veg_litterfall_mt(ix_sap_wood, ixC, ic),veg_litterfall_mt(ix_sap_wood, ixC13, ic))
    dC13_heart_wood  = calc_delta_C13(veg_litterfall_mt(ix_heart_wood, ixC, ic),veg_litterfall_mt(ix_heart_wood, ixC13, ic))
    dC13_labile      = calc_delta_C13(veg_litterfall_mt(ix_labile, ixC, ic),veg_litterfall_mt(ix_labile, ixC13, ic))
    dC13_reserve     = calc_delta_C13(veg_litterfall_mt(ix_reserve, ixC, ic),veg_litterfall_mt(ix_reserve, ixC13, ic))
    dC13_fruit       = calc_delta_C13(veg_litterfall_mt(ix_fruit, ixC, ic),veg_litterfall_mt(ix_fruit, ixC13, ic))
    dC13_seed_bed    = calc_delta_C13(seed_bed_litterfall_mt(ixC, ic),seed_bed_litterfall_mt(ixC13, ic))
    DC14_leaf        = calc_delta_C14(veg_litterfall_mt(ix_leaf, ixC, ic), &
                                      veg_litterfall_mt(ix_leaf, ixC14, ic),dC13_leaf)
    DC14_fine_root   = calc_delta_C14(veg_litterfall_mt(ix_fine_root, ixC, ic), &
                                      veg_litterfall_mt(ix_fine_root, ixC14, ic),dC13_fine_root)
    DC14_coarse_root = calc_delta_C14(veg_litterfall_mt(ix_coarse_root, ixC, ic), &
                                      veg_litterfall_mt(ix_coarse_root, ixC14, ic),dC13_coarse_root)
    DC14_sap_wood    = calc_delta_C14(veg_litterfall_mt(ix_sap_wood, ixC, ic), &
                                      veg_litterfall_mt(ix_sap_wood, ixC14, ic),dC13_sap_wood)
    DC14_heart_wood  = calc_delta_C14(veg_litterfall_mt(ix_heart_wood, ixC, ic), &
                                      veg_litterfall_mt(ix_heart_wood, ixC14, ic),dC13_heart_wood)
    DC14_labile      = calc_delta_C14(veg_litterfall_mt(ix_labile, ixC, ic), &
                                      veg_litterfall_mt(ix_labile, ixC14, ic),dC13_labile)
    DC14_reserve     = calc_delta_C14(veg_litterfall_mt(ix_reserve, ixC, ic), &
                                      veg_litterfall_mt(ix_reserve, ixC14, ic),dC13_reserve)
    DC14_fruit       = calc_delta_C14(veg_litterfall_mt(ix_fruit, ixC, ic), &
                                      veg_litterfall_mt(ix_fruit, ixC14, ic),dC13_fruit)
    DC14_seed_bed    = calc_delta_C14(seed_bed_litterfall_mt(ixC, ic), &
                                      seed_bed_litterfall_mt(ixC14, ic),dC13_seed_bed)
    dN15_leaf        = calc_delta_N15(veg_litterfall_mt(ix_leaf, ixN, ic),veg_litterfall_mt(ix_leaf, ixN15, ic))
    dN15_fine_root   = calc_delta_N15(veg_litterfall_mt(ix_fine_root, ixN, ic),veg_litterfall_mt(ix_fine_root, ixN15, ic))
    dN15_coarse_root = calc_delta_N15(veg_litterfall_mt(ix_coarse_root, ixN, ic),veg_litterfall_mt(ix_coarse_root, ixN15, ic))
    dN15_sap_wood    = calc_delta_N15(veg_litterfall_mt(ix_sap_wood, ixN, ic),veg_litterfall_mt(ix_sap_wood, ixN15, ic))
    dN15_heart_wood  = calc_delta_N15(veg_litterfall_mt(ix_heart_wood, ixN, ic),veg_litterfall_mt(ix_heart_wood, ixN15, ic))
    dN15_labile      = calc_delta_N15(veg_litterfall_mt(ix_labile, ixN, ic),veg_litterfall_mt(ix_labile, ixN15, ic))
    dN15_reserve     = calc_delta_N15(veg_litterfall_mt(ix_reserve, ixN, ic),veg_litterfall_mt(ix_reserve, ixN15, ic))
    dN15_fruit       = calc_delta_N15(veg_litterfall_mt(ix_fruit, ixN, ic),veg_litterfall_mt(ix_fruit, ixN15, ic))
    dN15_seed_bed    = calc_delta_N15(seed_bed_litterfall_mt(ixN, ic),seed_bed_litterfall_mt(ixN15, ic))
    ! veglitterfall
    write(lun,'(60F18.8)') time_d, &   ! 6 variables for (currently) 9 elements + time_d = 55 columns
      ! leaf
      veg_litterfall_mt(ix_leaf, ixC, ic)*1000._wp, &
      veg_litterfall_mt(ix_leaf, ixN, ic)*1000000._wp, &
      veg_litterfall_mt(ix_leaf, ixP, ic)*1000000._wp    , &
      dC13_leaf, DC14_leaf, dN15_leaf, &
      ! fine_root
      veg_litterfall_mt(ix_fine_root, ixC, ic)*1000._wp, &
      veg_litterfall_mt(ix_fine_root, ixN, ic)*1000000._wp, &
      veg_litterfall_mt(ix_fine_root, ixP, ic)*1000000._wp   , &
      dC13_fine_root, DC14_fine_root, dN15_fine_root, &
      ! coarse_root
      veg_litterfall_mt(ix_coarse_root, ixC, ic)*1000._wp, &
      veg_litterfall_mt(ix_coarse_root, ixN, ic)*1000000._wp, &
      veg_litterfall_mt(ix_coarse_root, ixP, ic)*1000000._wp   , &
      dC13_coarse_root, DC14_coarse_root, dN15_coarse_root, &
      ! sap_wood
      veg_litterfall_mt(ix_sap_wood, ixC, ic)*1000._wp, &
      veg_litterfall_mt(ix_sap_wood, ixN, ic)*1000000._wp, &
      veg_litterfall_mt(ix_sap_wood, ixP, ic)*1000000._wp    , &
      dC13_sap_wood, DC14_sap_wood, dN15_sap_wood, &
      ! heart_wood
      veg_litterfall_mt(ix_heart_wood, ixC, ic)*1000._wp, &
      veg_litterfall_mt(ix_heart_wood, ixN, ic)*1000000._wp, &
      veg_litterfall_mt(ix_heart_wood, ixP, ic)*1000000._wp    , &
      dC13_heart_wood, DC14_heart_wood, dN15_heart_wood, &
      ! labile
      veg_litterfall_mt(ix_labile, ixC, ic)*1000._wp, &
      veg_litterfall_mt(ix_labile, ixN, ic)*1000000._wp, &
      veg_litterfall_mt(ix_labile, ixP, ic)*1000000._wp    , &
      dC13_labile, DC14_labile, dN15_labile, &
      ! reserve
      veg_litterfall_mt(ix_reserve, ixC, ic)*1000._wp, &
      veg_litterfall_mt(ix_reserve, ixN, ic)*1000000._wp, &
      veg_litterfall_mt(ix_reserve, ixP, ic)*1000000._wp   , &
      dC13_reserve, DC14_reserve, dN15_reserve, &
      ! fruit
      veg_litterfall_mt(ix_fruit, ixC, ic)*1000._wp, &
      veg_litterfall_mt(ix_fruit, ixN, ic)*1000000._wp, &
      veg_litterfall_mt(ix_fruit, ixP, ic)*1000000._wp   , &
      dC13_fruit, DC14_fruit, dN15_fruit, &
      ! seed_bed
      seed_bed_litterfall_mt(ixC, ic)*1000._wp, &
      seed_bed_litterfall_mt(ixN, ic)*1000000._wp, &
      seed_bed_litterfall_mt(ixP, ic)*1000000._wp  , &
      dC13_seed_bed, DC14_seed_bed, dN15_seed_bed 
    lun=lun+1

  END SUBROUTINE write_soil_forcing 


  ! ======================================================================================================= !
  !!
  !! Close output files
  !!
  ! ======================================================================================================= !


  ! ======================================================================================================= !
  !>Pools and Fluxes: close files
  !>
  !>  at this point of the simulation all data has been written
  !>
  SUBROUTINE close_output_files

    INTEGER :: ilun
    LOGICAL :: opened

    DO ilun = first_lun_ascii_txt_output_files, (last_lun_ascii_txt_output_files_plus_one - 1)
      INQUIRE(ilun, opened=opened)
      IF (opened) CLOSE(ilun) 
    END DO
  END SUBROUTINE close_output_files


  ! ======================================================================================================= !
  !>Soil forcing files: close files
  !> 
  !>  at this point of the simulation all data has been written
  !>
  SUBROUTINE close_output_files_soil_forcing

    INTEGER :: ilun
    LOGICAL :: opened

    DO ilun = lun_ascii_soil_forc_start, lun_ascii_soil_forc_max
      INQUIRE(ilun, opened=opened)
      IF (opened) CLOSE(ilun) 
    END DO
  END SUBROUTINE close_output_files_soil_forcing

END MODULE mo_qs_output_txt
