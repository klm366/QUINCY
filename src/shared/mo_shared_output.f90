!>
!! @brief functions / subroutines that help to create output from plant_model_driver
!!
!! @author
!!  Jan Engel, MPI-BGC Jena
!!
!! @par Revision History
!! First version by Jan Engel (2016-03-14)
!! Second version by Jan Engel (20170630): remove ALLOCATE (only one tile is used), and init *_acc VAR and pools within its routines
!!
MODULE mo_shared_output
  !! this module may be not necessary anymore when QUINCY is linked with JSBACH !!

#include "dsl_jsbach.inc"

  USE mo_kind,                  ONLY: wp, eps4, eps8
  USE mo_jsb_control,           ONLY: debug_on
  USE mo_exception,             ONLY: message_text, message, finish
  USE mo_jsb_lctlib_class,      ONLY: t_lctlib_element
  USE mo_jsb_model_init,        ONLY: fluxnet_output_status   ! used for 'special case: fluxnet output'
  USE mo_base_config_class,     ONLY: base_config
  USE mo_sb_config_class,       ONLY: t_sb_config
  USE mo_jsb4_forcing,          ONLY: forcing_options

  IMPLICIT NONE

  PUBLIC  ! all is public

  !-----------------------------------------------------------------------------------------------------
  !> Define local variables
  !! 
  !! output interval for pools and fluxes
  !-----------------------------------------------------------------------------------------------------
  REAL(wp), SAVE :: output_interval_pool   ! equals 'dtime/one_day' or is a mulitple of it
  REAL(wp), SAVE :: output_interval_flux   ! equals 'dtime/one_day' or is a mulitple of it
  INTEGER, PRIVATE :: ipoi = 1             ! the chunk in ics:ice of which data will be written to file; per default = 1

  CHARACTER(len=*), PARAMETER, PRIVATE :: modname = 'mo_shared_output'

CONTAINS

  !-----------------------------------------------------------------------------------------------------
  !> the input interval string is transformed into a number, which is used for timely model output
  !! 
  !! applied for pools and flux output
  !-----------------------------------------------------------------------------------------------------
  FUNCTION calc_output_time_interval() RESULT(interval_number)

    USE mo_base_constants,          ONLY: one_day, dtime
    USE mo_base_config_class,       ONLY: base_config

    IMPLICIT NONE
    ! 0.1 InOut
    REAL(wp)                 :: interval_number      !< number output

    !> 1.0 create interval_number depending on pool/flux
    !!
    SELECT CASE( TRIM(base_config%output_interval_pool) )
      CASE('timestep')
        interval_number = dtime/one_day  ! length of one time step as advanced in 'time' 
      CASE('daily')
        interval_number = 1.0_wp
      CASE('weekly')
        interval_number = 7.0_wp
      CASE('monthly')
        interval_number = 30.0_wp
      CASE('yearly')
        interval_number = 365.0_wp
    END SELECT
  END FUNCTION calc_output_time_interval


  !-----------------------------------------------------------------------------------------------------!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Pools                                                                                             !!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !-----------------------------------------------------------------------------------------------------!


  !-----------------------------------------------------------------------------------------------------
  !> Pools: open files and write column headers
  !! 
  !! open output files, define file names for pools, write column header names
  !! 
  !! lun: start at 100 with a current max of 129 as flux outputs starts with lun 130
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE open_pool_output_files(tile, iblk, ics, ice, nc, nsoil_sb)

    USE mo_jsb_class,             ONLY: Get_model
    USE mo_base_memory,           ONLY: t_slm_tile
    USE mo_jsb_model_class,       ONLY: t_jsb_model
    USE mo_jsb_process_class,     ONLY: A2L_, L2A_, SEB_, TURB_, SOIL_, HYDRO_, VEG_, HD_, RAD_, ASSIMI_, PHENO_, SB_
    
    ! Use of process configurations (t_PROC_config)
    dsl4jsb_Use_config(ASSIMI_)
    dsl4jsb_Use_config(VEG_)
    dsl4jsb_Use_config(SB_)

    IMPLICIT NONE
    ! 0.1 InOut
    TYPE(t_slm_tile),INTENT(in)                 :: tile           !< data structure for dynamic vegetation and soil variables
    INTEGER,INTENT(in)                          :: iblk, &        !< iBlock, start point of chunk, end point of chunk, number of points in the chunk
                                                   ics, &
                                                   ice, &
                                                   nc, &
                                                   nsoil_sb
    ! 0.2 Local
    TYPE(t_jsb_model),      POINTER       :: model
    INTEGER        :: lun           ! logical unit number
    INTEGER        :: max_lun_pool  ! currently 129
    CHARACTER(100) :: fname         ! filename
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':open_pool_output_files'
    ! ---------------------------
    ! 0.3 Declare Memory
    ! Declare process configuration and memory Pointers
    dsl4jsb_Def_config(ASSIMI_)
    dsl4jsb_Def_config(VEG_)
    dsl4jsb_Def_config(SB_)
    ! 0.5 Get Memory
    model  => Get_model(tile%owner_model_id)
    ! Get process config
    dsl4jsb_Get_config(ASSIMI_)
    dsl4jsb_Get_config(VEG_)
    dsl4jsb_Get_config(SB_)


    !> 1.0 init
    max_lun_pool          = 129
    output_interval_pool  = 0.0_wp

    !> 1.1 transform interval string to interval number (e.g. week -> 7)
    !! 
    output_interval_pool = calc_output_time_interval()


    !> 2.1 define open() statements for file output
    !! 
    !! 
    lun = 100 
    SELECT CASE(TRIM(base_config%slm_model_name))
      CASE('land','plant')
         fname = 'vegpoolC'//'_'//TRIM(base_config%output_interval_pool)//'.txt'
         OPEN (UNIT=lun, FILE=fname, ACTION="write", STATUS="unknown")
         lun=lun+1
         fname = 'vegpooldC13DC14'//'_'//TRIM(base_config%output_interval_pool)//'.txt'
         OPEN (UNIT=lun, FILE=fname, ACTION="write", STATUS="unknown")
         lun=lun+1
         fname = 'vegpoolN'//'_'//TRIM(base_config%output_interval_pool)//'.txt'
         OPEN (UNIT=lun, FILE=fname, ACTION="write", STATUS="unknown")
         lun=lun+1
         fname = 'vegpoolP'//'_'//TRIM(base_config%output_interval_pool)//'.txt'
         OPEN (UNIT=lun, FILE=fname, ACTION="write", STATUS="unknown")
         lun=lun+1
         fname = 'vegpooldN15'//'_'//TRIM(base_config%output_interval_pool)//'.txt'
         OPEN (UNIT=lun, FILE=fname, ACTION="write", STATUS="unknown")
         lun=lun+1
         fname = 'soilpoolPhysics'//'_'//TRIM(base_config%output_interval_pool)//'.txt'
         OPEN (UNIT=lun, FILE=fname, ACTION="write", STATUS="unknown")
         lun=lun+1
         fname = 'soilpoolPhysics_lay'//'_'//TRIM(base_config%output_interval_pool)//'.txt'
         OPEN (UNIT=lun, FILE=fname, ACTION="write", STATUS="unknown")
         lun=lun+1
         fname = 'soilpoolOrg'//'_'//TRIM(base_config%output_interval_pool)//'.txt'
         OPEN (UNIT=lun, FILE=fname, ACTION="write", STATUS="unknown")
         lun=lun+1
         fname = 'soilpoolInOrg'//'_'//TRIM(base_config%output_interval_pool)//'.txt'
         OPEN (UNIT=lun, FILE=fname, ACTION="write", STATUS="unknown")
         lun=lun+1
         fname = 'ecopoolCNP'//'_'//TRIM(base_config%output_interval_pool)//'.txt'
         OPEN (UNIT=lun, FILE=fname, ACTION="write", STATUS="unknown")
         lun=lun+1
         fname = 'soilpooldC13DC14dN15'//'_'//TRIM(base_config%output_interval_pool)//'.txt'
         OPEN (UNIT=lun, FILE=fname, ACTION="write", STATUS="unknown")
         lun=lun+1
         fname = 'veg_diagnostics'//'_'//TRIM(base_config%output_interval_pool)//'.txt'
         OPEN (UNIT=lun, FILE=fname, ACTION="write", STATUS="unknown")
         lun=lun+1
    END SELECT


    !! land & plant & canopy
    !! 
    SELECT CASE(TRIM(base_config%slm_model_name))
      CASE('land','plant','canopy')
        fname = 'veg_canopy'//'_'//TRIM(base_config%output_interval_pool)//'.txt'
        OPEN (UNIT=lun, FILE=fname, ACTION="write", STATUS="unknown")
        lun=lun+1
    END SELECT
    

    !! land & soil
    !! detailed data per layer only for simple_1d and jsm
    SELECT CASE(TRIM(base_config%slm_model_name))
      CASE('land','soil')
        SELECT CASE(TRIM(dsl4jsb_Config(SB_)%sb_model_scheme))
          CASE('simple_1d','jsm')
            fname = 'soilpoolCOrg_lay'//'_'//TRIM(base_config%output_interval_pool)//'.txt'
            OPEN (UNIT=lun, FILE=fname, ACTION="write", STATUS="unknown")
            lun=lun+1
            fname = 'soilpoolNOrg_lay'//'_'//TRIM(base_config%output_interval_pool)//'.txt'
            OPEN (UNIT=lun, FILE=fname, ACTION="write", STATUS="unknown")
            lun=lun+1
            fname = 'soilpoolPOrg_lay'//'_'//TRIM(base_config%output_interval_pool)//'.txt'
            OPEN (UNIT=lun, FILE=fname, ACTION="write", STATUS="unknown")
            lun=lun+1
            fname = 'soilpoolInNOrg_lay'//'_'//TRIM(base_config%output_interval_pool)//'.txt'
            OPEN (UNIT=lun, FILE=fname, ACTION="write", STATUS="unknown")
            lun=lun+1
            fname = 'soilpoolInPOrg_lay'//'_'//TRIM(base_config%output_interval_pool)//'.txt'
            OPEN (UNIT=lun, FILE=fname, ACTION="write", STATUS="unknown")
            lun=lun+1
            fname = 'soilpooldC13org_lay'//'_'//TRIM(base_config%output_interval_pool)//'.txt'
            OPEN (UNIT=lun, FILE=fname, ACTION="write", STATUS="unknown")
            lun=lun+1
            fname = 'soilpooldC14org_lay'//'_'//TRIM(base_config%output_interval_pool)//'.txt'
            OPEN (UNIT=lun, FILE=fname, ACTION="write", STATUS="unknown")
            lun=lun+1
            fname = 'soilpooldN15org_lay'//'_'//TRIM(base_config%output_interval_pool)//'.txt'
            OPEN (UNIT=lun, FILE=fname, ACTION="write", STATUS="unknown")
            lun=lun+1
            fname = 'soil_diagnostics_lay'//'_'//TRIM(base_config%output_interval_pool)//'.txt'
            OPEN (UNIT=lun, FILE=fname, ACTION="write", STATUS="unknown")
            lun=lun+1
        END SELECT
    END SELECT


    !> 2.2  test if lun for pools is too big
    !! 
    IF (lun > max_lun_pool) THEN
       WRITE(message_text,'(a)') 'lun of pools is too big, reaching values that are reserved for lun of flux output'
       CALL finish(TRIM(routine),message_text)
    END IF


    !> 3.0 write column headers into output files
    !! 
    !! 
    CALL init_output_writing_pool(nsoil_sb, dsl4jsb_Config(SB_)%sb_model_scheme)

  END SUBROUTINE open_pool_output_files


  !-----------------------------------------------------------------------------------------------------
  !> Pools: write column headers of output files
  !! 
  !! 
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE init_output_writing_pool(nsoil_sb, sb_model_scheme)

    IMPLICIT NONE
    ! 0.1 InOut
    INTEGER, INTENT(in)       :: nsoil_sb
    CHARACTER(15), INTENT(in) :: sb_model_scheme
    ! 0.2 Local
    INTEGER                   :: lun,isoil
    CHARACTER(2)              :: layerstr
    CHARACTER(450)            :: p1str,p2str,p3str,p4str,p5str,p6str,p7str,p8str,p9str,p10str

    !> 1.0 write headers
    !!
    !! lun of pools starts at 100 with a max of 30, as lun of flux output starts with 130 !
    lun = 100
    SELECT CASE(TRIM(base_config%slm_model_name))
      CASE('land','plant')

         WRITE(lun,'(2a)') 'Time Leaf_C FineRoot_C SapWood_C HeartWood_C CoarseRoot_C ', & 
                          'Fruit_C Labile_C Reserve_C SeedBank_C Total_C Total_C13 Total_C14'
         lun=lun+1
         WRITE(lun,'(4a)') 'Time Leaf_dC13 FineRoot_dC13 SapWood_dC13 HeartWood_dC13 CoarseRoot_dC13 ', & 
                          'Fruit_dC13 Labile_dC13 Reserve_dC13 SeedBank_dC13 Leaf_DC14 FineRoot_DC14 ', &
                          'SapWood_DC14 HeartWood_DC14 CoarseRoot_DC14 Fruit_DC14 Labile_DC14 ', &
                          'Reserve_DC14 SeedBank_DC14'
         lun=lun+1
         WRITE(lun,'(2a)') 'Time Leaf_N FineRoot_N SapWood_N HeartWood_N CoarseRoot_N ', & 
                          'Fruit_N Labile_N Reserve_N SeedBank_N Total_N'
         lun=lun+1
         WRITE(lun,'(2a)') 'Time Leaf_P FineRoot_P SapWood_P HeartWood_P CoarseRoot_P ', & 
                          'Fruit_P Labile_P Reserve_P SeedBank_P Total_P'
         lun=lun+1
         WRITE(lun,'(2a)') 'Time Leaf_d15N FineRoot_d15N SapWood_d15N HeartWood_d15N CoarseRoot_d15N ', & 
                          'Fruit_dN15 Labile_dN15 Reserve_dN15 SeedBank_dN15 Total_dN15'
         lun=lun+1
         ! soilpoolPhysics
         WRITE(lun,'(a)') 'Time RootZoneSoilMoisture RootZoneSoilWaterPotential BetaSoilGS BetaSoilPS '
         lun=lun+1
         ! soilpoolPhysics_lay output per soil_layer for two SOIL_ variables
         p1str=''
         p2str=''
         DO isoil = 1,nsoil_sb
            IF(isoil<10) THEN 
              WRITE(layerstr,'(I1)') isoil
            ELSE
              WRITE(layerstr,'(I2)') isoil
            ENDIF
            WRITE(p1str,'(a)') TRIM(p1str)//' WaterContentSoil_L'//TRIM(layerstr)
            WRITE(p2str,'(a)') TRIM(p2str)//' TemperatureSoil_L'//TRIM(layerstr)
         ENDDO
         WRITE(lun,'(3a)') 'Time ',TRIM(p1str), TRIM(p2str)
         lun=lun+1
         WRITE(lun,'(3a)') 'Time MetLitterC MetLitterN MetLitterP StrLitterC StrLitterN StrLitterP', &
                           ' WodLitterC WodLitterN WodLitterP',&
                           ' FastC FastN FastP SlowC SlowN SlowP MycC MycN MycP'
         lun=lun+1
         WRITE(lun,'(a)') 'Time MineralNH4 MineralNO3 MineralP NOY N2O N2'
         lun=lun+1
         WRITE(lun,'(a)') 'Time TotalSoilOrgC TotalSoilOrgN TotalSoilInOrgN EcoC EcoN EcoP EcoN15'
         lun=lun+1
         WRITE(lun,'(3a)') 'Time MetLitterdC13 StrLitterdC13 WodLitterdC13 FastdC13 SlowdC13 MycdC13 MetLitterDC14 StrLitterDC14',& 
                           ' WodLitterDC14 FastDC14 SlowDC14 MycDC14',&
                           ' MetLitterdN15 StrLitterdN15 WodLitterN15 FastdN15 SlowdN15 MineralNH4_dN15 MineralNO3_dN15'
         lun=lun+1
         WRITE(lun,'(2a)') 'Time Height Diameter LAI ind_density La2Sa L2R S2C BCF CN_target dNPP dNup cn_plant allometry_npp', & 
                           ' growing_season cost_nup cost_bnf unit_transp transp'
         lun=lun+1
    END SELECT


    !! land & plant & canopy
    !! 
    SELECT CASE(TRIM(base_config%slm_model_name))
      CASE('land','plant','canopy')
        ! veg_canopy
        WRITE(lun,'(a)') 'Time CanopyLayer f_et f_rub f_chl f_pepc f_oth leaf_nitrogen jmax vcmax'
        lun=lun+1
    END SELECT


    !! land & soil
    !! detailed data per layer only for simple_1d and jsm
    SELECT CASE(TRIM(base_config%slm_model_name))
      CASE('land','soil')
        SELECT CASE(TRIM(sb_model_scheme))
          CASE('simple_1d','jsm')
            ! soilpoolCOrg_lay
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
            DO isoil = 1,nsoil_sb
              IF(isoil<10) THEN 
                WRITE(layerstr,'(I1)') isoil
              ELSE
                WRITE(layerstr,'(I2)') isoil
              ENDIF
              WRITE(p1str,'(a)') TRIM(p1str)//' SoluLit_L'//TRIM(layerstr)
              WRITE(p2str,'(a)') TRIM(p2str)//' PolyLit_L'//TRIM(layerstr)
              WRITE(p3str,'(a)') TRIM(p3str)//' WodLit_L'//TRIM(layerstr)
              WRITE(p4str,'(a)') TRIM(p4str)//' DOM_L'//TRIM(layerstr)
              WRITE(p5str,'(a)') TRIM(p5str)//' aDOM_L'//TRIM(layerstr)
              WRITE(p6str,'(a)') TRIM(p6str)//' Fungi_L'//TRIM(layerstr)
              WRITE(p7str,'(a)') TRIM(p7str)//' Myc_L'//TRIM(layerstr)
              WRITE(p8str,'(a)') TRIM(p8str)//' Mic_L'//TRIM(layerstr)
              WRITE(p9str,'(a)') TRIM(p9str)//' Res_L'//TRIM(layerstr)
              WRITE(p10str,'(a)') TRIM(p10str)//' aRes_L'//TRIM(layerstr)
            ENDDO
            WRITE(lun,'(11a)') 'Time ',TRIM(p1str),TRIM(p2str),TRIM(p3str),TRIM(p4str),TRIM(p5str),TRIM(p6str),TRIM(p7str),&
              TRIM(p8str),TRIM(p9str),TRIM(p10str)
            lun=lun+1
            
            !soilpoolNOrg_lay
            WRITE(lun,'(11a)') 'Time ',TRIM(p1str),TRIM(p2str),TRIM(p3str),TRIM(p4str),TRIM(p5str),TRIM(p6str),TRIM(p7str),&
              TRIM(p8str),TRIM(p9str),TRIM(p10str)
            lun=lun+1
            
            !soilpoolPOrg_lay
            WRITE(lun,'(11a)') 'Time ',TRIM(p1str),TRIM(p2str),TRIM(p3str),TRIM(p4str),TRIM(p5str),TRIM(p6str),TRIM(p7str),&
              TRIM(p8str),TRIM(p9str),TRIM(p10str)
            lun=lun+1
            
            !soilpoolInNOrg_lay
            p1str=''
            p2str=''
            p3str=''
            p4str=''
            p5str=''
            p6str=''
            p7str=''
            p8str=''
            p9str=''
            DO isoil = 1,nsoil_sb
              IF(isoil<10) THEN 
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
            WRITE(lun,'(10a)') 'Time ',TRIM(p1str),TRIM(p2str),TRIM(p3str),TRIM(p4str),TRIM(p5str),TRIM(p6str), &
               TRIM(p7str),TRIM(p8str),TRIM(p9str)
            lun=lun+1
            
            !soilpoolInPOrg_lay
            p1str=''
            p2str=''
            p3str=''
            p4str=''
            p5str=''
            DO isoil = 1,nsoil_sb
              IF(isoil<10) THEN 
                WRITE(layerstr,'(I1)') isoil
              ELSE
                WRITE(layerstr,'(I2)') isoil
              ENDIF
              WRITE(p1str,'(a)') TRIM(p1str)//' PO4assoc_fast_L'//TRIM(layerstr)
              WRITE(p2str,'(a)') TRIM(p2str)//' PO4assoc_slow_L'//TRIM(layerstr)
              WRITE(p3str,'(a)') TRIM(p3str)//' PO4assoc_occl_L'//TRIM(layerstr)
              WRITE(p4str,'(a)') TRIM(p4str)//' PO4solute_L'//TRIM(layerstr)
              WRITE(p5str,'(a)') TRIM(p5str)//' PO4primary_L'//TRIM(layerstr)
            ENDDO
            WRITE(lun,'(6a)') 'Time ',TRIM(p1str),TRIM(p2str),TRIM(p3str),TRIM(p4str),TRIM(p5str)
            lun=lun+1
            
            ! dC13 dC14 dN15
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
            DO isoil = 1,nsoil_sb
              IF(isoil<10) THEN 
                WRITE(layerstr,'(I1)') isoil
              ELSE
                WRITE(layerstr,'(I2)') isoil
              ENDIF
              WRITE(p1str,'(a)') TRIM(p1str)//' SoluLit_L'//TRIM(layerstr)
              WRITE(p2str,'(a)') TRIM(p2str)//' PolyLit_L'//TRIM(layerstr)
              WRITE(p3str,'(a)') TRIM(p3str)//' WodLit_L'//TRIM(layerstr)
              WRITE(p4str,'(a)') TRIM(p4str)//' DOM_L'//TRIM(layerstr)
              WRITE(p5str,'(a)') TRIM(p5str)//' aDOM_L'//TRIM(layerstr)
              WRITE(p6str,'(a)') TRIM(p6str)//' Fungi_L'//TRIM(layerstr)
              WRITE(p7str,'(a)') TRIM(p7str)//' Myc_L'//TRIM(layerstr)
              WRITE(p8str,'(a)') TRIM(p8str)//' Mic_L'//TRIM(layerstr)
              WRITE(p9str,'(a)') TRIM(p9str)//' Res_L'//TRIM(layerstr)
              WRITE(p10str,'(a)') TRIM(p10str)//' aRes_L'//TRIM(layerstr)
            ENDDO
            ! soilpooldC13org_lay            
            WRITE(lun,'(11a)') 'Time ',TRIM(p1str),TRIM(p2str),TRIM(p3str),TRIM(p4str),TRIM(p5str),TRIM(p6str),TRIM(p7str),&
              TRIM(p8str),TRIM(p9str),TRIM(p10str)
            lun=lun+1          
            ! soilpooldC14org_lay       
            WRITE(lun,'(11a)') 'Time ',TRIM(p1str),TRIM(p2str),TRIM(p3str),TRIM(p4str),TRIM(p5str),TRIM(p6str),TRIM(p7str),&
              TRIM(p8str),TRIM(p9str),TRIM(p10str)
            lun=lun+1      
            ! soilpooldN15org_lay       
            WRITE(lun,'(11a)') 'Time ',TRIM(p1str),TRIM(p2str),TRIM(p3str),TRIM(p4str),TRIM(p5str),TRIM(p6str),TRIM(p7str),&
              TRIM(p8str),TRIM(p9str),TRIM(p10str)
            lun=lun+1 
            
            ! soil_diagnostics_lay
            p1str=''
            p2str=''
            p3str=''
            p4str=''
            p5str=''
            DO isoil = 1,nsoil_sb
              IF(isoil<10) THEN 
                WRITE(layerstr,'(I1)') isoil
              ELSE
                WRITE(layerstr,'(I2)') isoil
              ENDIF
              WRITE(p1str,'(a)') TRIM(p1str)//' QmaxOrg_L'//TRIM(layerstr)
              WRITE(p2str,'(a)') TRIM(p2str)//' Rho_L'//TRIM(layerstr)
              WRITE(p3str,'(a)') TRIM(p3str)//' CUEmic_L'//TRIM(layerstr)
              WRITE(p4str,'(a)') TRIM(p4str)//' QmaxPo4_L'//TRIM(layerstr)
              WRITE(p5str,'(a)') TRIM(p5str)//' KmPo4_L'//TRIM(layerstr)
            ENDDO
            WRITE(lun,'(6a)') 'Time ',TRIM(p1str),TRIM(p2str),TRIM(p3str),TRIM(p4str),TRIM(p5str)
            lun=lun+1                 
        END SELECT
    END SELECT

  END SUBROUTINE init_output_writing_pool


  !-----------------------------------------------------------------------------------------------------
  !> Pools: write data to files in a certain interval
  !! 
  !! write-interval accordingly to namelist
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE write_output_pool(tile, options, time_d, time_y, l_transient_spinup)

    USE mo_jsb_class,             ONLY: Get_model
    USE mo_base_memory,           ONLY: t_slm_tile, t_jsb_task_options
    USE mo_jsb_model_class,       ONLY: t_jsb_model
    USE mo_veg_constants,         ONLY: itree
    USE mo_shared_isotope,        ONLY: calc_delta_C13,calc_delta_C14,calc_delta_N15
    USE mo_jsb_process_class,     ONLY: A2L_, L2A_, SEB_, TURB_, SOIL_, HYDRO_, VEG_, HD_, RAD_, ASSIMI_, PHENO_, SB_

    ! Use of process configurations (t_PROC_config)
    dsl4jsb_Use_config(ASSIMI_)    
    dsl4jsb_Use_config(SB_)
    ! Use of process memories
    dsl4jsb_Use_memory(ASSIMI_)     !  USE mo_assimi_memory_class,   ONLY: t_assimi_memory
    dsl4jsb_Use_memory(VEG_)        !  USE mo_veg_memory_class,      ONLY: t_veg_memory
    dsl4slm_Use_memory_a_veg
    dsl4slm_Use_memory_a_sb
    USE mo_veg_memory_class,      ONLY: veg_pool, veg_total_biomass
    USE mo_sb_memory_class,       ONLY: sb_pool
    dsl4jsb_Use_memory(PHENO_)      !  USE mo_pheno_memory_class,    ONLY: t_pheno_memory
    dsl4jsb_Use_memory(SOIL_)       !  USE mo_soil_memory_class,     ONLY: t_soil_memory
    dsl4jsb_Use_memory(SB_)         !  USE mo_sb_memory_class,       ONLY: t_sb_memory

    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    TYPE(t_slm_tile), INTENT(in), TARGET :: tile        !< one tile
    TYPE(t_jsb_task_options), INTENT(in) :: options     !< model options
    REAL(wp), INTENT(in)                 :: time_d      !< time in days
    INTEGER,  INTENT(in)                 :: time_y      !< time in years; the 1st year of a simulation == year 1 (not zero)
    LOGICAL,  INTENT(in)                 :: l_transient_spinup
    ! ---------------------------
    ! 0.2 Local
    TYPE(t_jsb_model),      POINTER      :: model
    TYPE(t_lctlib_element), POINTER      :: lctlib                    !< land-cover-type library - parameter across pft's
    LOGICAL                              :: l_land_plant_canopy_out,& !< write land/plant/canopy output T/F
                                            l_soil_output             !< write soil       output T/F
    INTEGER                              :: lun                       ! logical unit number
    INTEGER                              :: icanopy                   ! loop over ncanopy
    REAL(wp)                             :: hlp1                      ! helper variable
    REAL(wp)                             :: l2r, s2c, bcf             ! helper ratio variables
    REAl(wp)                             :: growing_season_float      ! output variable for growing season
    REAL(wp)                             :: dC13leaf,dC13fine_root, & ! 
                                            dC13sap_wood,dC13heart_wood,dC13coarse_root, & ! delta C13 in pools 
                                            dC13fruit,dC13labile,dC13reserve,dC13seed_bed
    REAL(wp), DIMENSION(options%nc,options%nsoil_sb) & 
                                         :: dC13soluable_litter,dC13polymeric_litter,dC13woody_litter, & 
                                            dC13dom,dC13dom_assoc,dC13fungi,dC13mycorrhiza,dC13microbial, & 
                                            dC13residue,dC13residue_assoc
    INTEGER                              :: iblk, ics, ice, nc, nsoil_sb
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':write_output_pool'
    ! ---------------------------
    ! 0.3 Declare Memory
    ! Declare process configuration and memory Pointers
    dsl4jsb_Def_config(ASSIMI_)
    dsl4jsb_Def_config(SB_)
    dsl4jsb_Def_memory(ASSIMI_)                           ! TYPE(t_assimi_memory),   POINTER :: assimi__mem
    dsl4jsb_Def_memory(VEG_)                              ! TYPE(t_veg_memory),      POINTER :: veg__mem
    dsl4jsb_Def_memory(PHENO_)                            ! TYPE(t_pheno_memory),    POINTER :: pheno__mem
    dsl4jsb_Def_memory(SOIL_)                             ! TYPE(t_soil_memory),     POINTER :: soil__mem
    dsl4jsb_Def_memory(SB_)                               ! TYPE(t_sb_memory),       POINTER :: sb__mem
    ! Declare pointers to variables in memory
    ! PHENO_
    dsl4jsb_Logical2D_onChunk   :: growing_season
    ! ASSIMI_
    dsl4jsb_Real2D_onChunk      :: beta_soil_gs
    dsl4jsb_Real2D_onChunk      :: beta_soil_ps
    ! SOIL_ 2D
    dsl4jsb_Real2D_onChunk      :: w_soil_root 
    dsl4jsb_Real2D_onChunk      :: w_soil_root_pot
    ! SOIL_ 3D
    dsl4jsb_Real3D_onChunk      :: soil_depth_sl
    dsl4jsb_Real3D_onChunk      :: t_soil_sl
    dsl4jsb_Real3D_onChunk      :: w_soil_sl
    ! VEG_ 2D
    dsl4jsb_Real2D_onChunk      :: lai 
    dsl4jsb_Real2D_onChunk      :: height
    dsl4jsb_Real2D_onChunk      :: diameter
    dsl4jsb_Real2D_onChunk      :: dens_ind
    dsl4jsb_Real2D_onChunk      :: target_cn_leaf
    dsl4jsb_Real2D_onChunk      :: unit_npp_talloc_mavg
    dsl4jsb_Real2D_onChunk      :: unit_uptake_n_talloc_mavg
    dsl4jsb_Real2D_onChunk      :: growth_cn_talloc_mavg
    dsl4jsb_Real2D_onChunk      :: npp_talloc_mavg
    dsl4jsb_Real2D_onChunk      :: cost_n_uptake_root
    dsl4jsb_Real2D_onChunk      :: cost_n_fixation_tuptake_mavg
    dsl4jsb_Real2D_onChunk      :: unit_transpiration_talloc_mavg
    dsl4jsb_Real2D_onChunk      :: growth_req_n_tuptake_mavg
    ! VEG_ 3D
    dsl4jsb_Real3D_onChunk      :: fn_et_cl
    dsl4jsb_Real3D_onChunk      :: fn_rub_cl
    dsl4jsb_Real3D_onChunk      :: fn_chl_cl
    dsl4jsb_Real3D_onChunk      :: fn_pepc_cl
    dsl4jsb_Real3D_onChunk      :: fn_oth_cl
    dsl4jsb_Real3D_onChunk      :: leaf_nitrogen_cl
    dsl4jsb_Real3D_onChunk      :: jmax_cl
    dsl4jsb_Real3D_onChunk      :: vcmax_cl
    ! SB_ 3D
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
    ! Get local variables from options argument
    iblk      = options%iblk
    ics       = options%ics
    ice       = options%ice
    nc        = options%nc
    nsoil_sb  = options%nsoil_sb
    ! ---------------------------
    ! 0.4 Process Activity, Debug Option
!    IF (.NOT. tile%Is_process_active(xx_)) RETURN   ! ask for all processes considered here ?
    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')
    ! ---------------------------
    ! 0.5 Get Memory
    model  => Get_model(tile%owner_model_id)
    lctlib => model%lctlib(tile%lcts%lib_id)
    ! Get process config
    dsl4jsb_Get_config(ASSIMI_)
    dsl4jsb_Get_config(SB_)
    ! Get process memories
    dsl4jsb_Get_memory(ASSIMI_)
    dsl4jsb_Get_memory(VEG_)
    dsl4jsb_Get_memory(PHENO_)                            ! pheno__mem => pheno_memory
    dsl4jsb_Get_memory(SOIL_)                             ! soil__mem => soil_memory
    dsl4jsb_Get_memory(SB_)                               ! sb__mem => sb_memory
    ! Get process variables (Set pointers to variables in memory)
    dsl4jsb_Get_var2D_onChunk(PHENO_, growing_season)
    ! ---------------------------
    dsl4jsb_Get_var2D_onChunk(ASSIMI_, beta_soil_gs)
    dsl4jsb_Get_var2D_onChunk(ASSIMI_, beta_soil_ps)
    dsl4jsb_Get_var3D_onChunk(ASSIMI_, jmax_cl)
    dsl4jsb_Get_var3D_onChunk(ASSIMI_, vcmax_cl)
    ! ---------------------------
    dsl4jsb_Get_var2D_onChunk(SOIL_, w_soil_root)
    dsl4jsb_Get_var2D_onChunk(SOIL_, w_soil_root_pot)
    dsl4jsb_Get_var3D_onChunk(SOIL_, soil_depth_sl)
    dsl4jsb_Get_var3D_onChunk(SOIL_, t_soil_sl)
    dsl4jsb_Get_var3D_onChunk(SOIL_, w_soil_sl)
    ! ---------------------------
    dsl4jsb_Get_var2D_onChunk(VEG_, lai)
    dsl4jsb_Get_var2D_onChunk(VEG_, height)
    dsl4jsb_Get_var2D_onChunk(VEG_, diameter)
    dsl4jsb_Get_var2D_onChunk(VEG_, dens_ind)
    dsl4jsb_Get_var2D_onChunk(VEG_, target_cn_leaf)
    dsl4jsb_Get_var2D_onChunk(VEG_, unit_npp_talloc_mavg)
    dsl4jsb_Get_var2D_onChunk(VEG_, unit_uptake_n_talloc_mavg)
    dsl4jsb_Get_var2D_onChunk(VEG_, growth_cn_talloc_mavg)
    dsl4jsb_Get_var2D_onChunk(VEG_, npp_talloc_mavg)
    dsl4jsb_Get_var2D_onChunk(VEG_, cost_n_uptake_root)
    dsl4jsb_Get_var2D_onChunk(VEG_, cost_n_fixation_tuptake_mavg)
    dsl4jsb_Get_var2D_onChunk(VEG_, unit_transpiration_talloc_mavg)
    dsl4jsb_Get_var2D_onChunk(VEG_, growth_req_n_tuptake_mavg)
    dsl4jsb_Get_var3D_onChunk(VEG_, fn_et_cl)
    dsl4jsb_Get_var3D_onChunk(VEG_, fn_rub_cl)
    dsl4jsb_Get_var3D_onChunk(VEG_, fn_chl_cl)
    dsl4jsb_Get_var3D_onChunk(VEG_, fn_pepc_cl)
    dsl4jsb_Get_var3D_onChunk(VEG_, fn_oth_cl)
    dsl4jsb_Get_var3D_onChunk(VEG_, leaf_nitrogen_cl)
    ! ---------------------------
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
    ! Get process memory pool structures
#define Get_veg_pool
#define Get_veg_total_biomass
#define Get_sb_pool
#include "slm_pools_a_incl.f90"
#undef Get_veg_pool
#undef Get_veg_total_biomass
#undef Get_sb_pool



    ! ------------------------------------------------------------------------------------------------------------
    ! Go science
    ! ------------------------------------------------------------------------------------------------------------


    !< 0.9 init local variables
    !! 
    lun = 100

    !< 1.0 test if anything would be written out at this timestep
    !! 
    ! set to FALSE as default
    l_land_plant_canopy_out = .FALSE.
    l_soil_output           = .FALSE.
    ! test if TRUE
    SELECT CASE(TRIM(base_config%slm_model_name))
      CASE('land','plant','canopy')
        ! special fluxnet simulations; run only if 'special case: fluxnet output' is not active, or active & after spinup
        IF (fluxnet_output_status < 3) THEN 
          ! after spinup, when doing transient simulations
          ! output start & end
          IF(.NOT. l_transient_spinup .AND. &
             & time_y >= base_config%output_start_first_day_year .AND. time_y <= base_config%output_end_last_day_year) THEN
            ! matching output_interval
            IF(MODULO(time_d, output_interval_pool) < eps4 .OR. &
               & ABS(output_interval_pool - MODULO(time_d, output_interval_pool)) < eps4) THEN
              ! set TRUE
              l_land_plant_canopy_out = .TRUE.
            END IF
          END IF
        END IF
      CASE('soil')
        ! after spinup, when doing transient simulations
        ! output start & end
        IF(.NOT. l_transient_spinup .AND. &
           & time_y >= base_config%output_start_first_day_year .AND. time_y <= base_config%output_end_last_day_year) THEN
          ! matching output_interval
          IF(MODULO(time_d, output_interval_pool) < eps4 .OR. &
             & ABS(output_interval_pool - MODULO(time_d, output_interval_pool)) < eps4) THEN
            ! set TRUE
            l_soil_output = .TRUE.
          END IF
        END IF
    END SELECT


    !< 2.0  output
    !! 
    !! land/plant
    !! land/plant/canopy
    !! land/soil
    SELECT CASE(TRIM(base_config%slm_model_name))
      CASE('land','plant')
        ! output T/F
        IF (l_land_plant_canopy_out) THEN 

              write(lun,'(13F12.6)') time_d,veg_pool%leaf%carbon(ipoi),&
                veg_pool%fine_root%carbon(ipoi),veg_pool%sap_wood%carbon(ipoi),&
                veg_pool%heart_wood%carbon(ipoi),veg_pool%coarse_root%carbon(ipoi),&
                veg_pool%fruit%carbon(ipoi),veg_pool%labile%carbon(ipoi),&
                veg_pool%reserve%carbon(ipoi),veg_pool%seed_bed%carbon(ipoi),&
                veg_total_biomass%carbon(ipoi),veg_total_biomass%carbon13(ipoi),veg_total_biomass%carbon14(ipoi)
       
              dC13leaf        = calc_delta_C13(veg_pool%leaf%carbon(ipoi),veg_pool%leaf%carbon13(ipoi))
              dC13fine_root   = calc_delta_C13(veg_pool%fine_root%carbon(ipoi),veg_pool%fine_root%carbon13(ipoi))
              dC13sap_wood    = calc_delta_C13(veg_pool%sap_wood%carbon(ipoi),veg_pool%sap_wood%carbon13(ipoi))
              dC13heart_wood  = calc_delta_C13(veg_pool%heart_wood%carbon(ipoi),veg_pool%heart_wood%carbon13(ipoi))
              dC13coarse_root = calc_delta_C13(veg_pool%coarse_root%carbon(ipoi),veg_pool%coarse_root%carbon13(ipoi))
              dC13fruit       = calc_delta_C13(veg_pool%fruit%carbon(ipoi),veg_pool%fruit%carbon13(ipoi))
              dC13labile      = calc_delta_C13(veg_pool%labile%carbon(ipoi),veg_pool%labile%carbon13(ipoi))
              dC13reserve     = calc_delta_C13(veg_pool%reserve%carbon(ipoi),veg_pool%reserve%carbon13(ipoi))
              dC13seed_bed    = calc_delta_C13(veg_pool%seed_bed%carbon(ipoi),veg_pool%seed_bed%carbon13(ipoi))
              lun=lun+1   ! increase lun by 1 after each write statement
              write(lun,'(19F18.6)') time_d,dC13leaf,dC13fine_root,dC13sap_wood,dC13heart_wood,&
                dC13coarse_root,dC13fruit,dC13labile,dC13reserve,dC13seed_bed, & 
                calc_delta_C14(veg_pool%leaf%carbon(ipoi),veg_pool%leaf%carbon14(ipoi),dC13leaf),&
                calc_delta_C14(veg_pool%fine_root%carbon(ipoi),veg_pool%fine_root%carbon14(ipoi),dC13fine_root),&
                calc_delta_C14(veg_pool%sap_wood%carbon(ipoi),veg_pool%sap_wood%carbon14(ipoi),dC13sap_wood),&
                calc_delta_C14(veg_pool%heart_wood%carbon(ipoi),veg_pool%heart_wood%carbon14(ipoi),dC13heart_wood),&
                calc_delta_C14(veg_pool%coarse_root%carbon(ipoi),veg_pool%coarse_root%carbon14(ipoi),dC13coarse_root),&
                calc_delta_C14(veg_pool%fruit%carbon(ipoi),veg_pool%fruit%carbon14(ipoi),dC13fruit),&
                calc_delta_C14(veg_pool%labile%carbon(ipoi),veg_pool%labile%carbon14(ipoi),dC13labile),&
                calc_delta_C14(veg_pool%reserve%carbon(ipoi),veg_pool%reserve%carbon14(ipoi),dC13reserve),&
                calc_delta_C14(veg_pool%seed_bed%carbon(ipoi),veg_pool%seed_bed%carbon14(ipoi),dC13seed_bed)
              lun=lun+1
              write(lun,'(12F12.6)') time_d,veg_pool%leaf%nitrogen(ipoi),&
                veg_pool%fine_root%nitrogen(ipoi),veg_pool%sap_wood%nitrogen(ipoi),&
                veg_pool%heart_wood%nitrogen(ipoi),veg_pool%coarse_root%nitrogen(ipoi),&
                veg_pool%fruit%nitrogen(ipoi),veg_pool%labile%nitrogen(ipoi),&
                veg_pool%reserve%nitrogen(ipoi),veg_pool%seed_bed%nitrogen(ipoi),&
                veg_total_biomass%nitrogen(ipoi)
              lun=lun+1
              write(lun,'(11F12.6)') time_d,veg_pool%leaf%phosphorus(ipoi),&
                veg_pool%fine_root%phosphorus(ipoi),veg_pool%sap_wood%phosphorus(ipoi),&
                veg_pool%heart_wood%phosphorus(ipoi),veg_pool%coarse_root%phosphorus(ipoi),&
                veg_pool%fruit%phosphorus(ipoi),veg_pool%labile%phosphorus(ipoi),&
                veg_pool%reserve%phosphorus(ipoi),veg_pool%seed_bed%phosphorus(ipoi),& 
                veg_total_biomass%phosphorus(ipoi)
              lun=lun+1
                write(lun,'(12F12.6)') time_d, &
                calc_delta_N15(veg_pool%leaf%nitrogen(ipoi),veg_pool%leaf%nitrogen15(ipoi)), &
                calc_delta_N15(veg_pool%fine_root%nitrogen(ipoi),veg_pool%fine_root%nitrogen15(ipoi)), &
                calc_delta_N15(veg_pool%sap_wood%nitrogen(ipoi),veg_pool%sap_wood%nitrogen15(ipoi)), &
                calc_delta_N15(veg_pool%heart_wood%nitrogen(ipoi),veg_pool%heart_wood%nitrogen15(ipoi)), &
                calc_delta_N15(veg_pool%coarse_root%nitrogen(ipoi),veg_pool%coarse_root%nitrogen15(ipoi)), &
                calc_delta_N15(veg_pool%fruit%nitrogen(ipoi),veg_pool%fruit%nitrogen15(ipoi)), &
                calc_delta_N15(veg_pool%labile%nitrogen(ipoi),veg_pool%labile%nitrogen15(ipoi)), &
                calc_delta_N15(veg_pool%reserve%nitrogen(ipoi),veg_pool%reserve%nitrogen15(ipoi)), &
                calc_delta_N15(veg_pool%seed_bed%nitrogen(ipoi),veg_pool%seed_bed%nitrogen15(ipoi)), &
                calc_delta_N15(veg_total_biomass%nitrogen(ipoi),veg_total_biomass%nitrogen15(ipoi))
              lun=lun+1
              ! soilpoolPhysics
              write(lun,'(20F16.6)') time_d, &
                                     w_soil_root(ipoi),  w_soil_root_pot(ipoi), &
                                     beta_soil_gs(ipoi), beta_soil_ps(ipoi)
              lun=lun+1
              ! soilpoolPhysics_lay (max of 20 soil layers)
              write(lun,'(41F16.6)') time_d, w_soil_sl(ipoi,:), t_soil_sl(ipoi,:)
              lun=lun+1
              write(lun,'(30F16.6)') time_d, &
                                     SUM(sb_pool%soluable_litter%carbon(ipoi,:)*soil_depth_sl(ipoi,:))      , & 
                                     SUM(sb_pool%soluable_litter%nitrogen(ipoi,:)*soil_depth_sl(ipoi,:))    , & 
                                     SUM(sb_pool%soluable_litter%phosphorus(ipoi,:)*soil_depth_sl(ipoi,:))  , & 
                                     SUM(sb_pool%polymeric_litter%carbon(ipoi,:)*soil_depth_sl(ipoi,:))     , & 
                                     SUM(sb_pool%polymeric_litter%nitrogen(ipoi,:)*soil_depth_sl(ipoi,:))   , & 
                                     SUM(sb_pool%polymeric_litter%phosphorus(ipoi,:)*soil_depth_sl(ipoi,:)) , & 
                                     SUM(sb_pool%woody_litter%carbon(ipoi,:)*soil_depth_sl(ipoi,:))         , & 
                                     SUM(sb_pool%woody_litter%nitrogen(ipoi,:)*soil_depth_sl(ipoi,:))       , & 
                                     SUM(sb_pool%woody_litter%phosphorus(ipoi,:)*soil_depth_sl(ipoi,:))     , & 
                                     SUM(sb_pool%microbial%carbon(ipoi,:)*soil_depth_sl(ipoi,:))            , & 
                                     SUM(sb_pool%microbial%nitrogen(ipoi,:)*soil_depth_sl(ipoi,:))          , & 
                                     SUM(sb_pool%microbial%phosphorus(ipoi,:)*soil_depth_sl(ipoi,:))        , & 
                                     SUM(sb_pool%residue%carbon(ipoi,:)*soil_depth_sl(ipoi,:))              , & 
                                     SUM(sb_pool%residue%nitrogen(ipoi,:)*soil_depth_sl(ipoi,:))            , & 
                                     SUM(sb_pool%residue%phosphorus(ipoi,:)*soil_depth_sl(ipoi,:))          , & 
                                     SUM(sb_pool%mycorrhiza%carbon(ipoi,:)*soil_depth_sl(ipoi,:))           , &
                                     SUM(sb_pool%mycorrhiza%nitrogen(ipoi,:)*soil_depth_sl(ipoi,:))         , &
                                     SUM(sb_pool%mycorrhiza%phosphorus(ipoi,:)*soil_depth_sl(ipoi,:))       
              lun=lun+1
              write(lun,'(30F16.6)') time_d, &
                                     SUM(nh4_solute(ipoi,:)*soil_depth_sl(ipoi,:))                       , & 
                                     SUM(no3_solute(ipoi,:)*soil_depth_sl(ipoi,:))                       , & 
                                     SUM(po4_solute(ipoi,:)*soil_depth_sl(ipoi,:))                       , &
                                     SUM(noy(ipoi,:)*soil_depth_sl(ipoi,:))                              , &
                                     SUM(n2o(ipoi,:)*soil_depth_sl(ipoi,:))                              , &
                                     SUM(n2(ipoi,:)*soil_depth_sl(ipoi,:))  
              lun=lun+1
              write(lun,'(30F16.6)') time_d, &
                                     SUM(sb_pool%soluable_litter%carbon(ipoi,:)*soil_depth_sl(ipoi,:))     + & 
                                      SUM(sb_pool%polymeric_litter%carbon(ipoi,:)*soil_depth_sl(ipoi,:))     + & 
                                      SUM(sb_pool%woody_litter%carbon(ipoi,:)*soil_depth_sl(ipoi,:))         + &
                                      SUM(sb_pool%microbial%carbon(ipoi,:)*soil_depth_sl(ipoi,:))            + &
                                      SUM(sb_pool%residue%carbon(ipoi,:)*soil_depth_sl(ipoi,:))              + & 
                                      SUM(sb_pool%mycorrhiza%carbon(ipoi,:)*soil_depth_sl(ipoi,:))           , &
                                     SUM(sb_pool%soluable_litter%nitrogen(ipoi,:)*soil_depth_sl(ipoi,:))   + & 
                                      SUM(sb_pool%polymeric_litter%nitrogen(ipoi,:)*soil_depth_sl(ipoi,:))   + & 
                                      SUM(sb_pool%woody_litter%nitrogen(ipoi,:)*soil_depth_sl(ipoi,:))       + &
                                      SUM(sb_pool%microbial%nitrogen(ipoi,:)*soil_depth_sl(ipoi,:))          + &
                                      SUM(sb_pool%residue%nitrogen(ipoi,:)*soil_depth_sl(ipoi,:))            + & 
                                      SUM(sb_pool%mycorrhiza%nitrogen(ipoi,:)*soil_depth_sl(ipoi,:))         , &
                                     SUM(nh4_solute(ipoi,:)*soil_depth_sl(ipoi,:))                                + & 
                                      SUM(no3_solute(ipoi,:)*soil_depth_sl(ipoi,:))                                + & 
                                      SUM(noy(ipoi,:)*soil_depth_sl(ipoi,:))                                       + &
                                      SUM(n2o(ipoi,:)*soil_depth_sl(ipoi,:))                                       + &
                                      SUM(n2(ipoi,:)*soil_depth_sl(ipoi,:))                                        , &     
                                     SUM(sb_pool%soluable_litter%carbon(ipoi,:)*soil_depth_sl(ipoi,:))     + & 
                                      SUM(sb_pool%polymeric_litter%carbon(ipoi,:)*soil_depth_sl(ipoi,:))     + & 
                                      SUM(sb_pool%woody_litter%carbon(ipoi,:)*soil_depth_sl(ipoi,:))         + &
                                      SUM(sb_pool%microbial%carbon(ipoi,:)*soil_depth_sl(ipoi,:))            + &
                                      SUM(sb_pool%residue%carbon(ipoi,:)*soil_depth_sl(ipoi,:))              + & 
                                      SUM(sb_pool%mycorrhiza%carbon(ipoi,:)*soil_depth_sl(ipoi,:))           + &
                                      veg_total_biomass%carbon(ipoi)                                         , &
                                     SUM(sb_pool%soluable_litter%nitrogen(ipoi,:)*soil_depth_sl(ipoi,:))   + & 
                                      SUM(sb_pool%polymeric_litter%nitrogen(ipoi,:)*soil_depth_sl(ipoi,:))   + & 
                                      SUM(sb_pool%woody_litter%nitrogen(ipoi,:)*soil_depth_sl(ipoi,:))       + &
                                      SUM(sb_pool%microbial%nitrogen(ipoi,:)*soil_depth_sl(ipoi,:))          + &
                                      SUM(sb_pool%residue%nitrogen(ipoi,:)*soil_depth_sl(ipoi,:))            + & 
                                      SUM(sb_pool%mycorrhiza%nitrogen(ipoi,:)*soil_depth_sl(ipoi,:))         + &
                                      SUM(nh4_solute(ipoi,:)*soil_depth_sl(ipoi,:))                                + & 
                                      SUM(no3_solute(ipoi,:)*soil_depth_sl(ipoi,:))                                + &
                                      SUM(noy(ipoi,:)*soil_depth_sl(ipoi,:))                                       + &
                                      SUM(n2o(ipoi,:)*soil_depth_sl(ipoi,:))                                       + &
                                      SUM(n2(ipoi,:)*soil_depth_sl(ipoi,:))                                        + &     
                                      veg_total_biomass%nitrogen(ipoi)                                       , &
                                     SUM(sb_pool%soluable_litter%phosphorus(ipoi,:)*soil_depth_sl(ipoi,:))   + & 
                                      SUM(sb_pool%polymeric_litter%phosphorus(ipoi,:)*soil_depth_sl(ipoi,:))   + & 
                                      SUM(sb_pool%woody_litter%phosphorus(ipoi,:)*soil_depth_sl(ipoi,:))       + &
                                      SUM(sb_pool%microbial%phosphorus(ipoi,:)*soil_depth_sl(ipoi,:))          + &
                                      SUM(sb_pool%residue%phosphorus(ipoi,:)*soil_depth_sl(ipoi,:))            + & 
                                      SUM(sb_pool%mycorrhiza%phosphorus(ipoi,:)*soil_depth_sl(ipoi,:))         + &
                                      SUM(po4_solute(ipoi,:)*soil_depth_sl(ipoi,:))                            + & 
                                      SUM(po4_assoc_fast(ipoi,:)*soil_depth_sl(ipoi,:))                        + &
                                      SUM(po4_assoc_slow(ipoi,:)*soil_depth_sl(ipoi,:))                        + &
                                      SUM(po4_occluded(ipoi,:)*soil_depth_sl(ipoi,:))                          + &
                                      SUM(po4_primary(ipoi,:)*soil_depth_sl(ipoi,:))                           + &
                                      veg_total_biomass%phosphorus(ipoi)                                     , &
                                     SUM(sb_pool%soluable_litter%nitrogen15(ipoi,:)*soil_depth_sl(ipoi,:))   + & 
                                      SUM(sb_pool%polymeric_litter%nitrogen15(ipoi,:)*soil_depth_sl(ipoi,:))   + & 
                                      SUM(sb_pool%woody_litter%nitrogen15(ipoi,:)*soil_depth_sl(ipoi,:))       + &
                                      SUM(sb_pool%microbial%nitrogen15(ipoi,:)*soil_depth_sl(ipoi,:))          + &
                                      SUM(sb_pool%residue%nitrogen15(ipoi,:)*soil_depth_sl(ipoi,:))            + & 
                                      SUM(sb_pool%mycorrhiza%nitrogen15(ipoi,:)*soil_depth_sl(ipoi,:))         + &
                                      SUM(nh4_n15_solute(ipoi,:)*soil_depth_sl(ipoi,:))                                + & 
                                      SUM(no3_n15_solute(ipoi,:)*soil_depth_sl(ipoi,:))                                + &
                                      SUM(noy_n15(ipoi,:)*soil_depth_sl(ipoi,:))                                       + &
                                      SUM(n2o_n15(ipoi,:)*soil_depth_sl(ipoi,:))                                       + &
                                      SUM(n2_n15(ipoi,:)*soil_depth_sl(ipoi,:))                                        + &     
                                      veg_total_biomass%nitrogen15(ipoi)
              lun=lun+1 
              dC13soluable_litter(ipoi,1)  = calc_delta_C13(SUM(sb_pool%soluable_litter%carbon(ipoi,:)*&
                                              soil_depth_sl(ipoi,:)), &
                                              SUM(sb_pool%soluable_litter%carbon13(ipoi,:)*soil_depth_sl(ipoi,:)))
              dC13polymeric_litter(ipoi,1) = calc_delta_C13(SUM(sb_pool%polymeric_litter%carbon(ipoi,:)*&
                                              soil_depth_sl(ipoi,:)), &
                                              SUM(sb_pool%polymeric_litter%carbon13(ipoi,:)*soil_depth_sl(ipoi,:)))
              dC13woody_litter(ipoi,1)     = calc_delta_C13(SUM(sb_pool%woody_litter%carbon(ipoi,:)*&
                                              soil_depth_sl(ipoi,:)), &
                                              SUM(sb_pool%woody_litter%carbon13(ipoi,:)*soil_depth_sl(ipoi,:)))
              dC13microbial(ipoi,1)        = calc_delta_C13(SUM(sb_pool%microbial%carbon(ipoi,:)*&
                                              soil_depth_sl(ipoi,:)), &
                                              SUM(sb_pool%microbial%carbon13(ipoi,:)*soil_depth_sl(ipoi,:)))
              dC13residue(ipoi,1)          = calc_delta_C13(SUM(sb_pool%residue%carbon(ipoi,:)*&
                                              soil_depth_sl(ipoi,:)), &
                                              SUM(sb_pool%residue%carbon13(ipoi,:)*soil_depth_sl(ipoi,:)))
              dC13mycorrhiza(ipoi,1)       = calc_delta_C13(SUM(sb_pool%mycorrhiza%carbon(ipoi,:)*&
                                              soil_depth_sl(ipoi,:)), &
                                              SUM(sb_pool%mycorrhiza%carbon13(ipoi,:)*soil_depth_sl(ipoi,:)))
              write(lun,'(20F14.6)') time_d, &
                dC13soluable_litter(ipoi,1), &
                dC13polymeric_litter(ipoi,1), &
                dC13woody_litter(ipoi,1), &
                dC13microbial(ipoi,1), &
                dC13residue(ipoi,1), & 
                dC13mycorrhiza(ipoi,1), &
                calc_delta_C14(SUM(sb_pool%soluable_litter%carbon(ipoi,:)*soil_depth_sl(ipoi,:)), &
                               SUM(sb_pool%soluable_litter%carbon14(ipoi,:)*soil_depth_sl(ipoi,:)), &
                                    dC13soluable_litter(ipoi,1)) , &
                calc_delta_C14(SUM(sb_pool%polymeric_litter%carbon(ipoi,:)*soil_depth_sl(ipoi,:)), &
                               SUM(sb_pool%polymeric_litter%carbon14(ipoi,:)*soil_depth_sl(ipoi,:)), &
                                    dC13polymeric_litter(ipoi,1)), &
                calc_delta_C14(SUM(sb_pool%woody_litter%carbon(ipoi,:)*soil_depth_sl(ipoi,:)), &
                               SUM(sb_pool%woody_litter%carbon14(ipoi,:)*soil_depth_sl(ipoi,:)), &
                                    dC13woody_litter(ipoi,1))        , &
                calc_delta_C14(SUM(sb_pool%microbial%carbon(ipoi,:)*soil_depth_sl(ipoi,:)), &
                               SUM(sb_pool%microbial%carbon14(ipoi,:)*soil_depth_sl(ipoi,:)), dC13microbial(ipoi,1)), &
                calc_delta_C14(SUM(sb_pool%residue%carbon(ipoi,:)*soil_depth_sl(ipoi,:)), &
                               SUM(sb_pool%residue%carbon14(ipoi,:)*soil_depth_sl(ipoi,:)), dC13residue(ipoi,1)), &
                calc_delta_C14(SUM(sb_pool%mycorrhiza%carbon(ipoi,:)*soil_depth_sl(ipoi,:)), &
                               SUM(sb_pool%mycorrhiza%carbon14(ipoi,:)*soil_depth_sl(ipoi,:)),dC13mycorrhiza(ipoi,1)), &
                calc_delta_N15(sb_pool%soluable_litter%nitrogen(ipoi,1),&
                                sb_pool%soluable_litter%nitrogen15(ipoi,1)), &
                calc_delta_N15(sb_pool%polymeric_litter%nitrogen(ipoi,1),&
                                sb_pool%polymeric_litter%nitrogen15(ipoi,1)), &
                calc_delta_N15(sb_pool%woody_litter%nitrogen(ipoi,1),sb_pool%woody_litter%nitrogen15(ipoi,1)), &
                calc_delta_N15(sb_pool%microbial%nitrogen(ipoi,1),sb_pool%microbial%nitrogen15(ipoi,1)), &
                calc_delta_N15(sb_pool%residue%nitrogen(ipoi,1),sb_pool%residue%nitrogen15(ipoi,1)), &
                calc_delta_N15(nh4_solute(ipoi,1),nh4_n15_solute(ipoi,1)), &
                calc_delta_N15(no3_solute(ipoi,1),no3_n15_solute(ipoi,1))
              lun=lun+1

              IF(lctlib%growthform == itree) THEN 
                 IF (veg_pool%sap_wood%carbon(ipoi) > eps8) THEN
                    hlp1 = lai(ipoi)/(veg_pool%sap_wood%carbon(ipoi) / height(ipoi) / lctlib%wood_density)
                 ELSE
                    hlp1 = 0._wp
                 ENDIF
              ELSE
                 IF (veg_pool%leaf%carbon(ipoi)>eps8) THEN
                    hlp1 = veg_pool%sap_wood%carbon(ipoi)/veg_pool%leaf%carbon(ipoi)
                 ELSE
                    hlp1 = 0._wp
                 ENDIF
              ENDIF
              IF (veg_pool%fine_root%carbon(ipoi) > eps8) THEN
                 l2r = veg_pool%leaf%carbon(ipoi)/veg_pool%fine_root%carbon(ipoi)
              ELSE
                 l2r = 0._wp
              ENDIF
              IF (veg_pool%sap_wood%carbon(ipoi)>eps8) THEN
                 s2c =  veg_pool%coarse_root%carbon(ipoi)/veg_pool%sap_wood%carbon(ipoi)
              ELSE
                 s2c = 0._wp
              ENDIF
              IF ((veg_pool%leaf%carbon(ipoi)+veg_pool%sap_wood%carbon(ipoi)+veg_pool%heart_wood%carbon(ipoi)+&
                   veg_pool%fine_root%carbon(ipoi)+veg_pool%coarse_root%carbon(ipoi))>eps8) THEN
                 bcf =  (veg_pool%fine_root%carbon(ipoi)+veg_pool%coarse_root%carbon(ipoi))/ & 
                (veg_pool%leaf%carbon(ipoi)+veg_pool%sap_wood%carbon(ipoi)+veg_pool%heart_wood%carbon(ipoi)+&
                   veg_pool%fine_root%carbon(ipoi)+veg_pool%coarse_root%carbon(ipoi))
              ELSE
                 bcf = 0._wp
              ENDIF
              ! transition from bool to 0/1 for growing season
              IF (growing_season(ipoi)) THEN
                growing_season_float = 1.0_wp
                ELSE
                growing_season_float = 0.0_wp
              END IF 
              write(lun,'(22F20.6)') time_d, height(ipoi), diameter(ipoi), lai(ipoi), dens_ind(ipoi), &
                hlp1,l2r,s2c, bcf, &
                target_cn_leaf(ipoi), unit_npp_talloc_mavg(ipoi), unit_uptake_n_talloc_mavg(ipoi),&
                growth_cn_talloc_mavg(ipoi), npp_talloc_mavg(ipoi), &
                growing_season_float, & 
                cost_n_uptake_root(ipoi), cost_n_fixation_tuptake_mavg(ipoi), &
                unit_transpiration_talloc_mavg(ipoi) * 1.e10_wp, growth_req_n_tuptake_mavg(ipoi)
              lun=lun+1
           END IF  ! write output T/F (l_land_plant_canopy_out)
    END SELECT  ! select land/plant & canopy


    !! land & plant & canopy
    !! 
    SELECT CASE(TRIM(base_config%slm_model_name))
      CASE('land','plant','canopy')
        ! output T/F
        IF (l_land_plant_canopy_out) THEN 
          ! veg_canopy 
          DO icanopy=1, dsl4jsb_Config(ASSIMI_)%ncanopy
            write(lun,'(10F12.6)') time_d, REAL(icanopy), fn_et_cl(ipoi,icanopy), &
                                  fn_rub_cl(ipoi,icanopy), fn_chl_cl(ipoi,icanopy),&
                                  fn_pepc_cl(ipoi,icanopy), fn_oth_cl(ipoi,icanopy),leaf_nitrogen_cl(ipoi,icanopy), &
                                  jmax_cl(ipoi,icanopy), vcmax_cl(ipoi,icanopy)
          END DO
          lun=lun+1
        ENDIF
    END SELECT  ! select land/plant & canopy


    !! land & soil
    !! detailed data per layer only for simple_1d and jsm (max of 20 soil layers)
    SELECT CASE(TRIM(base_config%slm_model_name))
    CASE('land','soil')
      ! output T/F
      IF (l_land_plant_canopy_out .OR. l_soil_output) THEN 
        SELECT CASE(TRIM(dsl4jsb_Config(SB_)%sb_model_scheme))
        CASE('simple_1d','jsm')
          ! soilpoolCOrg_lay
          write(lun,'(201F16.6)') time_d,sb_pool%soluable_litter%carbon(ipoi,:),&
              sb_pool%polymeric_litter%carbon(ipoi,:),sb_pool%woody_litter%carbon(ipoi,:),&
              sb_pool%dom%carbon(ipoi,:),sb_pool%dom_assoc%carbon(ipoi,:),&
              sb_pool%fungi%carbon(ipoi,:),sb_pool%mycorrhiza%carbon(ipoi,:),&
              sb_pool%microbial%carbon(ipoi,:),sb_pool%residue%carbon(ipoi,:),&
              sb_pool%residue_assoc%carbon(ipoi,:)
          lun=lun+1
          ! soilpoolNOrg_lay
          write(lun,'(201F16.6)') time_d,sb_pool%soluable_litter%nitrogen(ipoi,:),&
              sb_pool%polymeric_litter%nitrogen(ipoi,:),sb_pool%woody_litter%nitrogen(ipoi,:),&
              sb_pool%dom%nitrogen(ipoi,:),sb_pool%dom_assoc%nitrogen(ipoi,:),&
              sb_pool%fungi%nitrogen(ipoi,:),sb_pool%mycorrhiza%nitrogen(ipoi,:),&
              sb_pool%microbial%nitrogen(ipoi,:),sb_pool%residue%nitrogen(ipoi,:),&
              sb_pool%residue_assoc%nitrogen(ipoi,:)
          lun=lun+1
          ! soilpoolPOrg_lay
          write(lun,'(201F16.6)') time_d,sb_pool%soluable_litter%phosphorus(ipoi,:),&
              sb_pool%polymeric_litter%phosphorus(ipoi,:),sb_pool%woody_litter%phosphorus(ipoi,:),&
              sb_pool%dom%phosphorus(ipoi,:),sb_pool%dom_assoc%phosphorus(ipoi,:),&
              sb_pool%fungi%phosphorus(ipoi,:),sb_pool%mycorrhiza%phosphorus(ipoi,:),&
              sb_pool%microbial%phosphorus(ipoi,:),sb_pool%residue%phosphorus(ipoi,:),&
              sb_pool%residue_assoc%phosphorus(ipoi,:)
          lun=lun+1
          ! soilpoolInNOrg_lay
          write(lun,'(181F18.6)') time_d, &
            nh4_assoc(ipoi,:), nh4_solute(ipoi,:),no3_solute(ipoi,:)  , &
            noy(ipoi,:), n2o(ipoi,:), n2(ipoi,:), &
            calc_delta_N15(nh4_assoc(ipoi,:) ,nh4_n15_assoc(ipoi,:))  , &
            calc_delta_N15(nh4_solute(ipoi,:),nh4_n15_solute(ipoi,:)) , &
            calc_delta_N15(no3_solute(ipoi,:),no3_n15_solute(ipoi,:))
          lun=lun+1
          ! soilpoolInPOrg_lay
          write(lun,'(101F16.6)') time_d, &
            po4_assoc_fast(ipoi,:), po4_assoc_slow(ipoi,:), po4_occluded(ipoi,:), po4_solute(ipoi,:), po4_primary(ipoi,:)
          lun=lun+1
          ! soilpooldC13org_lay
          dC13soluable_litter(ipoi,:)   = calc_delta_C13(sb_pool%soluable_litter%carbon(ipoi,:), &
                                                      sb_pool%soluable_litter%carbon13(ipoi,:))
          dC13polymeric_litter(ipoi,:)  = calc_delta_C13(sb_pool%polymeric_litter%carbon(ipoi,:), &
                                                      sb_pool%polymeric_litter%carbon13(ipoi,:))
          dC13woody_litter(ipoi,:)      = calc_delta_C13(sb_pool%woody_litter%carbon(ipoi,:), &
                                                      sb_pool%woody_litter%carbon13(ipoi,:))
          dC13dom(ipoi,:)               = calc_delta_C13(sb_pool%dom%carbon(ipoi,:), &
                                                      sb_pool%dom%carbon13(ipoi,:))
          dC13dom_assoc(ipoi,:)         = calc_delta_C13(sb_pool%dom_assoc%carbon(ipoi,:), &
                                                      sb_pool%dom_assoc%carbon13(ipoi,:))
          dC13fungi(ipoi,:)             = calc_delta_C13(sb_pool%fungi%carbon(ipoi,:), &
                                                      sb_pool%fungi%carbon13(ipoi,:))
          dC13mycorrhiza(ipoi,:)        = calc_delta_C13(sb_pool%mycorrhiza%carbon(ipoi,:), &
                                                      sb_pool%mycorrhiza%carbon13(ipoi,:))
          dC13microbial(ipoi,:)         = calc_delta_C13(sb_pool%microbial%carbon(ipoi,:), &
                                                      sb_pool%microbial%carbon13(ipoi,:))
          dC13residue(ipoi,:)           = calc_delta_C13(sb_pool%residue%carbon(ipoi,:), &
                                                      sb_pool%residue%carbon13(ipoi,:))
          dC13residue_assoc(ipoi,:)     = calc_delta_C13(sb_pool%residue_assoc%carbon(ipoi,:), &
                                                      sb_pool%residue_assoc%carbon13(ipoi,:))
          write(lun,'(201F16.6)') time_d, &
                                  dC13soluable_litter(ipoi,:), dC13polymeric_litter(ipoi,:), &
                                  dC13woody_litter(ipoi,:),    dC13dom(ipoi,:), &
                                  dC13dom_assoc(ipoi,:),       dC13fungi(ipoi,:), &
                                  dC13mycorrhiza(ipoi,:),      dC13microbial(ipoi,:), &
                                  dC13residue(ipoi,:),         dC13residue_assoc(ipoi,:)
          lun=lun+1
          ! soilpooldC14org_lay
          write(lun,'(201F16.6)') time_d,&
            calc_delta_C14(sb_pool%soluable_litter%carbon(ipoi,:), &
                            sb_pool%soluable_litter%carbon14(ipoi,:),dC13soluable_litter(ipoi,:)),&
            calc_delta_C14(sb_pool%polymeric_litter%carbon(ipoi,:), &
                            sb_pool%polymeric_litter%carbon14(ipoi,:),dC13polymeric_litter(ipoi,:)),&
            calc_delta_C14(sb_pool%woody_litter%carbon(ipoi,:), &
                            sb_pool%woody_litter%carbon14(ipoi,:),dC13woody_litter(ipoi,:)),&
            calc_delta_C14(sb_pool%dom%carbon(ipoi,:), &
                            sb_pool%dom%carbon14(ipoi,:),dC13dom(ipoi,:)),&
            calc_delta_C14(sb_pool%dom_assoc%carbon(ipoi,:), &
                            sb_pool%dom_assoc%carbon14(ipoi,:),dC13dom_assoc(ipoi,:)),&
            calc_delta_C14(sb_pool%fungi%carbon(ipoi,:), &
                            sb_pool%fungi%carbon14(ipoi,:),dC13fungi(ipoi,:)),&
            calc_delta_C14(sb_pool%mycorrhiza%carbon(ipoi,:), &
                            sb_pool%mycorrhiza%carbon14(ipoi,:),dC13mycorrhiza(ipoi,:)),&
            calc_delta_C14(sb_pool%microbial%carbon(ipoi,:), &
                            sb_pool%microbial%carbon14(ipoi,:),dC13microbial(ipoi,:)),&
            calc_delta_C14(sb_pool%residue%carbon(ipoi,:), &
                            sb_pool%residue%carbon14(ipoi,:),dC13residue(ipoi,:)),&
            calc_delta_C14(sb_pool%residue_assoc%carbon(ipoi,:), &
                            sb_pool%residue_assoc%carbon14(ipoi,:),dC13residue_assoc(ipoi,:))
          lun=lun+1
          ! soilpooldN15org_lay
          write(lun,'(201F16.6)') time_d, &
            calc_delta_N15(sb_pool%soluable_litter%nitrogen(ipoi,:), &
                            sb_pool%soluable_litter%nitrogen15(ipoi,:)), &
            calc_delta_N15(sb_pool%polymeric_litter%nitrogen(ipoi,:), &
                            sb_pool%polymeric_litter%nitrogen15(ipoi,:)), &
            calc_delta_N15(sb_pool%woody_litter%nitrogen(ipoi,:), &
                            sb_pool%woody_litter%nitrogen15(ipoi,:)), &
            calc_delta_N15(sb_pool%dom%nitrogen(ipoi,:), &
                            sb_pool%dom%nitrogen15(ipoi,:)), &
            calc_delta_N15(sb_pool%dom_assoc%nitrogen(ipoi,:), &
                            sb_pool%dom_assoc%nitrogen15(ipoi,:)), &
            calc_delta_N15(sb_pool%fungi%nitrogen(ipoi,:), &
                            sb_pool%fungi%nitrogen15(ipoi,:)), &
            calc_delta_N15(sb_pool%mycorrhiza%nitrogen(ipoi,:), &
                            sb_pool%mycorrhiza%nitrogen15(ipoi,:)), &
            calc_delta_N15(sb_pool%microbial%nitrogen(ipoi,:), &
                            sb_pool%microbial%nitrogen15(ipoi,:)), &
            calc_delta_N15(sb_pool%residue%nitrogen(ipoi,:), &
                            sb_pool%residue%nitrogen15(ipoi,:)), &
            calc_delta_N15(sb_pool%residue_assoc%nitrogen(ipoi,:), &
                            sb_pool%residue_assoc%nitrogen15(ipoi,:))
          lun=lun+1
          ! soil_diagnostics_lay
          write(lun,'(101F16.6)') time_d, qmax_org(ipoi,:), bulk_dens_corr_sl(ipoi,:), microbial_cue_eff_tmic_mavg(ipoi,:), &
                                 qmax_po4(ipoi,:), km_adsorpt_po4_sl(ipoi,:)
          lun=lun+1
        END SELECT  ! simple_1d/jsm
      END IF        ! l_land_plant_canopy_out .OR. l_soil_output
    END SELECT      ! land/soil


! Pass back values to process memory
#define Passback_veg_pool
#define Passback_veg_total_biomass
#define Passback_sb_pool
#include "slm_pools_a_incl.f90"
#undef Passback_veg_pool
#undef Passback_veg_total_biomass
#undef Passback_sb_pool

  END SUBROUTINE write_output_pool


  !-----------------------------------------------------------------------------------------------------!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Fluxes                                                                                            !!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !-----------------------------------------------------------------------------------------------------!


  !-----------------------------------------------------------------------------------------------------
  !> Fluxes: open files and write column headers
  !! 
  !! open output files, define file names for pools, write column header names
  !! 
  !! lun: start at 100 + 30
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE open_flux_output_files(tile, iblk, ics, ice, nc, nsoil_sb)

    USE mo_jsb_class,             ONLY: Get_model
    USE mo_base_memory,           ONLY: t_slm_tile
    USE mo_jsb_model_class,       ONLY: t_jsb_model
    USE mo_jsb_process_class,     ONLY: A2L_, L2A_, SEB_, TURB_, SOIL_, HYDRO_, VEG_, HD_, RAD_, ASSIMI_, PHENO_, SB_
    
    ! Use of process configurations (t_PROC_config)
    dsl4jsb_Use_config(ASSIMI_)
    dsl4jsb_Use_config(VEG_)
    dsl4jsb_Use_config(SB_)

    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    TYPE(t_slm_tile),INTENT(in)                 :: tile           !< data structure for dynamic vegetation and soil variables
    INTEGER,INTENT(in)                          :: iblk, &        !< iBlock, start point of chunk, end point of chunk, number of points in the chunk
                                                   ics, &
                                                   ice, &
                                                   nc, &
                                                   nsoil_sb 
    ! ---------------------------
    ! 0.2 Local
    TYPE(t_jsb_model),      POINTER       :: model
    INTEGER        :: lun           ! logical unit number
    CHARACTER(100) :: fname         ! filename
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':open_flux_output_files'
    ! ---------------------------
    ! 0.3 Declare Memory
    ! Declare process configuration and memory Pointers
    dsl4jsb_Def_config(ASSIMI_)
    dsl4jsb_Def_config(VEG_)
    dsl4jsb_Def_config(SB_)
    ! ---------------------------
    ! 0.5 Get Memory
    model  => Get_model(tile%owner_model_id)
    ! Get process config
    dsl4jsb_Get_config(ASSIMI_)
    dsl4jsb_Get_config(VEG_)
    dsl4jsb_Get_config(SB_)


    !> 1.0 init
    output_interval_flux = 0.0_wp

    !> 1.1 transform interval string to interval number (e.g. week -> 7)
    !! 
    output_interval_flux = calc_output_time_interval()


    !> 2.0 define open() statements for file output
    !! 
    !! 
    lun = 100 + 30    ! -> max of 30 files for pool output

    !! land
    SELECT CASE(TRIM(base_config%slm_model_name))
      CASE('land')  ! for "land & soil" see below
        ! mycofluxCNP
        fname = 'mycofluxCNP'//'_'//TRIM(base_config%output_interval_flux)//'.txt'
        OPEN (UNIT=lun, FILE=fname, ACTION="write", STATUS="unknown")
        lun=lun+1
    END SELECT

    !! land/plant & canopy  (land/soil below)
    SELECT CASE(TRIM(base_config%slm_model_name))
      CASE('land','plant') ! additional to CASE('land')
        fname = 'vegfluxC'//'_'//TRIM(base_config%output_interval_flux)//'.txt'
        OPEN (UNIT=lun, FILE=fname, ACTION="write", STATUS="unknown")
        lun=lun+1
        fname = 'vegfluxdC13DC14'//'_'//TRIM(base_config%output_interval_flux)//'.txt'
        OPEN (UNIT=lun, FILE=fname, ACTION="write", STATUS="unknown")
        lun=lun+1
        fname = 'vegfluxN'//'_'//TRIM(base_config%output_interval_flux)//'.txt'
        OPEN (UNIT=lun, FILE=fname, ACTION="write", STATUS="unknown")
        lun=lun+1
        fname = 'vegfluxP'//'_'//TRIM(base_config%output_interval_flux)//'.txt'
        OPEN (UNIT=lun, FILE=fname, ACTION="write", STATUS="unknown")
        lun=lun+1
        fname = 'vegfluxH2O'//'_'//TRIM(base_config%output_interval_flux)//'.txt'
        OPEN (UNIT=lun, FILE=fname, ACTION="write", STATUS="unknown")
        lun=lun+1
        fname = 'vegfluxEnergy'//'_'//TRIM(base_config%output_interval_flux)//'.txt'
        OPEN (UNIT=lun, FILE=fname, ACTION="write", STATUS="unknown")
        lun=lun+1
      CASE('canopy') ! subset of CASE('land','plant')
        fname = 'vegfluxC'//'_'//TRIM(base_config%output_interval_flux)//'.txt'
        OPEN (UNIT=lun, FILE=fname, ACTION="write", STATUS="unknown")
        lun=lun+1
        fname = 'vegfluxdC13dC14'//'_'//TRIM(base_config%output_interval_flux)//'.txt'
        OPEN (UNIT=lun, FILE=fname, ACTION="write", STATUS="unknown")
        lun=lun+1
        fname = 'vegfluxH2O'//'_'//TRIM(base_config%output_interval_flux)//'.txt'
        OPEN (UNIT=lun, FILE=fname, ACTION="write", STATUS="unknown")
        lun=lun+1
        fname = 'vegfluxEnergy'//'_'//TRIM(base_config%output_interval_flux)//'.txt'
        OPEN (UNIT=lun, FILE=fname, ACTION="write", STATUS="unknown")
        lun=lun+1
    END SELECT

    !! land & soil - additional to CASE('land')
    !!  simple_0d and simple_1d and jsm
    !!  layer-data for simple_1d and jsm
    SELECT CASE(TRIM(base_config%slm_model_name))
      CASE('land','soil')
        !! all soil models (simple_0d and simple_1d and jsm)
        ! soilfluxCdC13DC14
        fname = 'soilfluxCdC13DC14'//'_'//TRIM(base_config%output_interval_flux)//'.txt'
        OPEN (UNIT=lun, FILE=fname, ACTION="write", STATUS="unknown")
        lun=lun+1
        ! soilfluxN
        fname = 'soilfluxN'//'_'//TRIM(base_config%output_interval_flux)//'.txt'
        OPEN (UNIT=lun, FILE=fname, ACTION="write", STATUS="unknown")
        lun=lun+1
        ! soilfluxP
        fname = 'soilfluxP'//'_'//TRIM(base_config%output_interval_flux)//'.txt'
        OPEN (UNIT=lun, FILE=fname, ACTION="write", STATUS="unknown")
        lun=lun+1
        ! soilfluxdN
        fname = 'soilfluxdN'//'_'//TRIM(base_config%output_interval_flux)//'.txt'
        OPEN (UNIT=lun, FILE=fname, ACTION="write", STATUS="unknown")
        lun=lun+1
        ! sb_massBalance
        fname = 'sb_massBalance'//'_'//TRIM(base_config%output_interval_flux)//'.txt'
        OPEN (UNIT=lun, FILE=fname, ACTION="write", STATUS="unknown")
        lun=lun+1
        !! simple_1d & jsm layer-data
        SELECT CASE(TRIM(dsl4jsb_Config(SB_)%sb_model_scheme))
          CASE('simple_1d','jsm')
            ! soilfluxC_lay
            fname = 'soilfluxC_lay'//'_'//TRIM(base_config%output_interval_flux)//'.txt'
            OPEN (UNIT=lun, FILE=fname, ACTION="write", STATUS="unknown")
            lun=lun+1
            ! soilfluxN_lay
            fname = 'soilfluxN_lay'//'_'//TRIM(base_config%output_interval_flux)//'.txt'
            OPEN (UNIT=lun, FILE=fname, ACTION="write", STATUS="unknown")
            lun=lun+1
            ! soilfluxP_lay
            fname = 'soilfluxP_lay'//'_'//TRIM(base_config%output_interval_flux)//'.txt'
            OPEN (UNIT=lun, FILE=fname, ACTION="write", STATUS="unknown")
            lun=lun+1
        END SELECT
    END SELECT


    !> 3.0 write column headers into output files
    !! 
    !! 
    CALL init_output_writing_flux(nsoil_sb, dsl4jsb_Config(SB_)%sb_model_scheme)

  END SUBROUTINE open_flux_output_files


  !-----------------------------------------------------------------------------------------------------
  !> Fluxes: write column headers of output files
  !! 
  !! 
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE init_output_writing_flux(nsoil_sb, sb_model_scheme)

    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    ! 
    INTEGER,       INTENT(in) :: nsoil_sb
    CHARACTER(15), INTENT(in) :: sb_model_scheme
    ! ---------------------------
    ! 0.2 Local
    INTEGER         :: lun,isoil
    CHARACTER(2)    :: layerstr
    CHARACTER(450)  :: p1str,p2str,p3str,p4str,p5str,p6str,p7str,p8str,p9str,p10str,p11str,p12str,p13str,p14str,p15str


    !> 1.0 write headers
    !!
    !! lun of fluxes starts at a multiple of 100, plus 30 (any lun > 129 in pool output throws an exception)
    lun = 100 + 30
    
    !! land 
    SELECT CASE(TRIM(base_config%slm_model_name))
      CASE('land')
        ! mycofluxCNP
        WRITE(lun,'(a)') 'Time mycoExportC mycoExportN mycoExportP mycoExportC13 mycoExportC14 mycoExportN15'
        lun=lun+1
    END SELECT

    !! land/plant & canopy  (land/soil below)
    SELECT CASE(TRIM(base_config%slm_model_name))
      CASE('land','plant') ! additional to CASE('land')
        ! vegfluxC
        WRITE(lun,'(a)') 'Time GPP GrowthResp MaintResp NTransformResp NFixationResp NPP ReserveUseC RootExudationC'
        lun=lun+1 
        ! vegfluxdC13DC14
        WRITE(lun,'(a)') 'Time GPP_dC13 GrowthResp_dC13 MaintResp_dC13 GPP_DC14 GrowthResp_DC14 MaintResp_DC14'
        lun=lun+1
        ! vegfluxN
        WRITE(lun,'(a)') 'Time Uptake_NH4 Uptake_NO3 RecyclingN ReserveUseN FixationN Fixation_dN15 Uptake_d15NH4 Uptake_d15NO3'
        lun=lun+1
        ! vegfluxP
        WRITE(lun,'(a)') 'Time Uptake_PO4 RecyclingP ReserveUseP'
        lun=lun+1
        WRITE(lun,'(a)') 'Time Precip Evaporation Interception Transpiration Runoff Drainage GC GA'
        lun=lun+1
        WRITE(lun,'(a)') 'Time Rnet Qh Qle Qg Albedo PPFD aPPFD fAPPFD' 
        lun=lun+1
      CASE('canopy') ! subset of CASE('land','plant')
        ! files: vegfluxC, vegfluxdC13dC14, vegfluxH2O, vegfluxEnergy
        WRITE(lun,'(a)') 'Time GPP GrowthResp MaintResp NTransFResp NFixResp NPP ReserveUse'
        lun=lun+1   
        WRITE(lun,'(2a)') 'Time GPP_dC13 GrowthResp_dC13 MaintResp_dC13 ', &
                          'GPP_DC14 GrowthResp_DC14 MaintResp_DC14'
        lun=lun+1
        WRITE(lun,'(a)') 'Time Precip Evaporation Interception Transpiration Runoff Drainage GC GA'
        lun=lun+1
        WRITE(lun,'(a)') 'Time Rnet Qh Qle Qg Albedo PPFD aPPFD fAPPFD' 
        lun=lun+1
    END SELECT

    
    !! land & soil - additional to CASE('land')
    !!  simple_0d and simple_1d and jsm
    !!  layer-data for simple_1d and jsm
    SELECT CASE(TRIM(base_config%slm_model_name))
      CASE('land','soil')
        !! all soil models (simple_0d and simple_1d and jsm)
        ! soilfluxCdC13DC14
        WRITE(lun,'(2a)') 'Time Litterfall HetResp LeachingDOM Litterfall_dC13 HetResp_dC13 LeachingDOM_dC13 ', &
                          'Litterfall_dC14 HetResp_DC14 LeachingDOM_dC14'
        lun=lun+1
        ! soilfluxN
        WRITE(lun,'(2a)') 'Time Litterfall NDeposition Nfix_asymb NetMinNH4 NetMinNO3 Nitrification Denitrification Emission_NOY ', &
                          'Emission_N2O Emission_N2 PlantUptakeNH4 PlantUptakeNO3 LeachingNH4 LeachingNO3 LeachingDOM' 
        lun=lun+1
        ! soilfluxP
        WRITE(lun,'(3a)') 'Time Litterfall ', &
                          'PDeposition Leaching_PO4 Leaching_DOM Net_Mineralisation_PO4 PlantUptake_PO4 Weathering_PO4 ', &
                          'Occlusion_PO4 Slow_sorption_PO4 Fast_adsorption_PO4 Biochem_mineralisation_PO4' 
        lun=lun+1
        ! soilfluxdN
        WRITE(lun,'(2a)') 'Time d15NDeposition d15Nfix_asymb Emission_d15NOY Emission_d15N2O Emission_d15N2 ', &
                          'Leaching_d15NH4 Leaching_d15NO3'
        lun=lun+1
        ! sb_massBalance
        WRITE(lun,'(2a)') 'Time sb_deltaPoolC sb_deltaPoolN sb_deltaPoolP sb_deltaPoolC13 sb_deltaPoolC14 sb_deltaPoolN15 ', &
                          'sb_totalFluxC sb_totalFluxN sb_totalFluxP sb_totalFluxC13 sb_totalFluxC14 sb_totalFluxN15'
        lun=lun+1
        !! simple_1d & jsm layer-data
        SELECT CASE(TRIM(sb_model_scheme))
          CASE('simple_1d','jsm')
            ! soilfluxC_lay
            p1str=''
            p2str=''
            p3str=''
            p4str=''
            DO isoil = 1, nsoil_sb
              IF(isoil<10) THEN 
                WRITE(layerstr,'(I1)') isoil
              ELSE
                WRITE(layerstr,'(I2)') isoil
              ENDIF
              WRITE(p1str,'(a)') TRIM(p1str)//' HetRespiration_L'//TRIM(layerstr)
              WRITE(p2str,'(a)') TRIM(p2str)//' HetRespirationC13_L'//TRIM(layerstr)
              WRITE(p3str,'(a)') TRIM(p3str)//' HetRespirationC14_L'//TRIM(layerstr)
              WRITE(p4str,'(a)') TRIM(p4str)//' Transport_DOC_L'//TRIM(layerstr)
            END DO
            WRITE(lun,'(5a)') 'Time ',TRIM(p1str),TRIM(p2str),TRIM(p3str),TRIM(p4str)
            lun=lun+1
            ! soilfluxN_lay
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
            DO isoil = 1, nsoil_sb
              IF(isoil<10) THEN 
                WRITE(layerstr,'(I1)') isoil
              ELSE
                WRITE(layerstr,'(I2)') isoil
              ENDIF
              WRITE(p1str,'(a)') TRIM(p1str)//' NetMineral_nh4_L'//TRIM(layerstr)
              WRITE(p2str,'(a)') TRIM(p2str)//' PlantUptake_nh4_L'//TRIM(layerstr)
              WRITE(p3str,'(a)') TRIM(p3str)//' MicUptake_nh4_L'//TRIM(layerstr)
              WRITE(p4str,'(a)') TRIM(p4str)//' MycoUptake_nh4_L'//TRIM(layerstr)
              WRITE(p5str,'(a)') TRIM(p5str)//' AsymbNfixation_L'//TRIM(layerstr)
              WRITE(p6str,'(a)') TRIM(p6str)//' NetMineral_no3_L'//TRIM(layerstr)
              WRITE(p7str,'(a)') TRIM(p7str)//' PlantUptake_no3_L'//TRIM(layerstr)
              WRITE(p8str,'(a)') TRIM(p8str)//' MicUptake_no3_L'//TRIM(layerstr)
              WRITE(p9str,'(a)') TRIM(p9str)//' MycoUptake_no3_L'//TRIM(layerstr)
              WRITE(p10str,'(a)') TRIM(p10str)//' Transport_DON_L'//TRIM(layerstr)
            END DO
            WRITE(lun,'(11a)') 'Time ', &
                               TRIM(p1str),TRIM(p2str),TRIM(p3str),TRIM(p4str),TRIM(p5str), &
                               TRIM(p6str),TRIM(p7str),TRIM(p8str),TRIM(p9str),TRIM(p10str)
            lun=lun+1
            ! soilfluxP_lay
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
            DO isoil = 1, nsoil_sb
              IF(isoil<10) THEN 
                WRITE(layerstr,'(I1)') isoil
              ELSE
                WRITE(layerstr,'(I2)') isoil
              ENDIF
              WRITE(p1str,'(a)') TRIM(p1str)//' Adsorp_L'//TRIM(layerstr)
              WRITE(p2str,'(a)') TRIM(p2str)//' Weath_L'//TRIM(layerstr)
              WRITE(p3str,'(a)') TRIM(p3str)//' Ocl_L'//TRIM(layerstr)
              WRITE(p4str,'(a)') TRIM(p4str)//' Desorp_L'//TRIM(layerstr)

              WRITE(p5str,'(a)') TRIM(p5str)//' NetMin_L'//TRIM(layerstr)
              WRITE(p6str,'(a)') TRIM(p6str)//' Biochem_min_L'//TRIM(layerstr)
              WRITE(p7str,'(a)') TRIM(p7str)//' Mic_upt_L'//TRIM(layerstr)
              WRITE(p8str,'(a)') TRIM(p8str)//' Plant_upt_L'//TRIM(layerstr)
              WRITE(p9str,'(a)') TRIM(p9str)//' Myco_upt_L'//TRIM(layerstr)

              WRITE(p10str,'(a)') TRIM(p10str)//' Trans_solute_L'//TRIM(layerstr)
              WRITE(p11str,'(a)') TRIM(p11str)//' Trans_assoc_fast_L'//TRIM(layerstr)
              WRITE(p12str,'(a)') TRIM(p12str)//' Trans_assoc_slow_L'//TRIM(layerstr)
              WRITE(p13str,'(a)') TRIM(p13str)//' Trans_occluded_L'//TRIM(layerstr)
              WRITE(p14str,'(a)') TRIM(p14str)//' Trans_primary_L'//TRIM(layerstr)
              WRITE(p15str,'(a)') TRIM(p15str)//' Trans_DOP_L'//TRIM(layerstr)
            END DO
            WRITE(lun,'(16a)') 'Time ', &
                               TRIM(p1str),TRIM(p2str),TRIM(p3str),TRIM(p4str),TRIM(p5str),TRIM(p6str),TRIM(p7str),TRIM(p8str), &
                               TRIM(p9str),TRIM(p10str),TRIM(p11str),TRIM(p12str),TRIM(p13str),TRIM(p14str),TRIM(p15str)
            lun=lun+1
        END SELECT
    END SELECT
    
  END SUBROUTINE init_output_writing_flux


  !-----------------------------------------------------------------------------------------------------
  !> Fluxes: write data to files in a certain interval and accumulate fluxes
  !! 
  !! write-interval accordingly to namelist
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE write_output_flux(tile, options, time_d, time_y, l_transient_spinup)

    USE mo_jsb_class,             ONLY: Get_model
    USE mo_jsb_process_class,     ONLY: A2L_, L2A_, SEB_, TURB_, SOIL_, HYDRO_, VEG_, HD_, RAD_, ASSIMI_, PHENO_, SB_
    USE mo_base_memory,           ONLY: t_slm_tile, t_jsb_task_options
    USE mo_jsb_model_class,       ONLY: t_jsb_model
    USE mo_base_constants,        ONLY: dtime
    USE mo_rad_constants,         ONLY: rad2ppfd
    USE mo_shared_isotope,        ONLY: calc_delta_C13, calc_delta_C14, calc_delta_N15
    USE mo_soil_config_class,     ONLY: max_soil_layers

    ! Use of process configurations (t_PROC_config)
    dsl4jsb_Use_config(SB_)

    ! Use of process memories
    dsl4jsb_Use_memory(A2L_)        !  USE mo_a2l_memory_class, ONLY: t_a2l_memory
    dsl4jsb_Use_memory(L2A_)        !  USE mo_l2a_memory_class, ONLY: t_l2a_memory
    dsl4jsb_Use_memory(ASSIMI_)     !  USE mo_assimi_memory_class,   ONLY: t_assimi_memory
    dsl4jsb_Use_memory(VEG_)        !  USE mo_veg_memory_class,      ONLY: t_veg_memory
    dsl4slm_Use_memory_a_veg
    dsl4slm_Use_memory_a_sb
    USE mo_veg_memory_class,      ONLY: veg_litterfall, veg_exudation, veg_reserve_use
    USE mo_sb_memory_class,       ONLY: sb_transport, sb_mycorrhiza_export, sb_delta_pools
    dsl4jsb_Use_memory(RAD_)        !  USE mo_rad_memory_class,      ONLY: t_rad_memory
    dsl4jsb_Use_memory(SOIL_)       !  USE mo_soil_memory_class,     ONLY: t_soil_memory
    dsl4jsb_Use_memory(SB_)         !  USE mo_sb_memory_class,       ONLY: t_sb_memory


    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    TYPE(t_slm_tile),         INTENT(inout) :: tile           !< one tile
    TYPE(t_jsb_task_options), INTENT(in)    :: options        !< model options
    REAL(wp),                 INTENT(in)    :: time_d         !< time in days
    INTEGER,                  INTENT(in)    :: time_y         !< time in years; the 1st year of a simulation == year 1 (not zero)
    LOGICAL,                  INTENT(in)    :: l_transient_spinup
    ! ---------------------------
    ! 0.2 Local
    TYPE(t_jsb_model),      POINTER       :: model
    TYPE(t_lctlib_element), POINTER       :: lctlib             !< land-cover-type library - parameter across pft's
    INTEGER                               :: lun                ! locical unit number
    REAL(wp)                              :: dC13gpp,dC13rg,dC13rm,hlp1,dC13litterfall,dC13het_resp,dC13leaching
    ! local accumulation variables
    REAL(wp), SAVE                        :: gpp_acc   = 0.0_wp ! Note: gpp_acc is initialized only the first time the subroutine is called. 
                                                                ! On subsequent calls the old value of gpp_acc is retained.
                                                                ! To avoid confusion the SAVE attribute is added, even though this is redundant.
    REAL(wp), SAVE                        :: growth_resp_acc     = 0.0_wp
    REAL(wp), SAVE                        :: maint_resp_acc      = 0.0_wp   
    REAL(wp), SAVE                        :: n_transform_resp_acc = 0.0_wp   
    REAL(wp), SAVE                        :: n_fixation_resp_acc = 0.0_wp 
    REAL(wp), SAVE                        :: exudationC_acc      = 0.0_wp  
    REAL(wp), SAVE                        :: npp_acc             = 0.0_wp   
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
    REAL(wp), SAVE                        :: runoff_acc          = 0.0_wp 
    REAL(wp), SAVE                        :: drainage_acc        = 0.0_wp 
    REAL(wp), SAVE                        :: gc_acc              = 0.0_wp
    REAL(wp), SAVE                        :: ga_acc              = 0.0_wp
    REAL(wp), SAVE                        :: n_deposition_acc    = 0.0_wp
    REAL(wp), SAVE                        :: n_deposition_n15_acc = 0.0_wp
    REAL(wp), SAVE                        :: ppfd_acc            = 0.0_wp
    REAL(wp), SAVE                        :: rnet_acc            = 0.0_wp
    REAL(wp), SAVE                        :: qh_acc              = 0.0_wp
    REAL(wp), SAVE                        :: qle_acc             = 0.0_wp
    REAL(wp), SAVE                        :: qg_acc              = 0.0_wp
    REAL(wp), SAVE                        :: albedo_acc          = 0.0_wp
    REAL(wp), SAVE                        :: appfd_acc           = 0.0_wp
    INTEGER,  SAVE                        :: ndat_acc            = 0 
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
    REAL(wp), SAVE                        :: net_min_nh4_acc     = 0.0_wp
    REAL(wp), SAVE                        :: net_min_no3_acc     = 0.0_wp
    REAL(wp), SAVE                        :: net_min_po4_acc     = 0.0_wp
    REAL(wp), SAVE                        :: nitrification_acc   = 0.0_wp
    REAL(wp), SAVE                        :: denitrification_acc = 0.0_wp
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
    REAL(wp), SAVE                        :: weathering_po4_acc    = 0.0_wp
    REAL(wp), SAVE                        :: occlusion_po4_acc    = 0.0_wp
    REAL(wp), SAVE                        :: slow_sorption_po4_acc    = 0.0_wp
    REAL(wp), SAVE                        :: fast_adsorption_po4_acc  = 0.0_wp
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
    REAL(wp), SAVE                        :: myco_export_c_acc             = 0.0_wp
    REAL(wp), SAVE                        :: myco_export_n_acc             = 0.0_wp
    REAL(wp), SAVE                        :: myco_export_p_acc             = 0.0_wp
    REAL(wp), SAVE                        :: myco_export_c13_acc           = 0.0_wp
    REAL(wp), SAVE                        :: myco_export_c14_acc           = 0.0_wp
    REAL(wp), SAVE                        :: myco_export_n15_acc           = 0.0_wp
    ! variables with accumulated information per soil layer
    ! hardcoded max soil layer number (20)
    REAL(wp), DIMENSION(max_soil_layers), SAVE :: &
                                             ! C
                                             het_respiration_lay_acc = 0.0_wp       , &
                                             het_respiration_c13_lay_acc = 0.0_wp   , &
                                             het_respiration_c14_lay_acc = 0.0_wp   , &
                                             transport_DOC_lay_acc = 0.0_wp               , &
                                             ! N
                                             net_mineral_nh4_lay_acc = 0.0_wp       , &
                                             plant_uptake_nh4_lay_acc = 0.0_wp      , &
                                             microbial_upt_nh4_lay_acc = 0.0_wp     , &
                                             microbial_upt_no3_lay_acc = 0.0_wp     , &
                                             myco_uptake_nh4_lay_acc = 0.0_wp       , &
                                             asymb_n_fixation_lay_acc = 0.0_wp      , &
                                             net_mineral_no3_lay_acc = 0.0_wp       , &
                                             plant_uptake_no3_lay_acc = 0.0_wp      , &
                                             myco_uptake_no3_lay_acc = 0.0_wp       , &
                                             transport_DON_lay_acc = 0.0_wp               , &
                                             ! P
                                             adsorpt_po4_acc = 0.0_wp           , &  
                                             net_min_po4_lay_acc = 0.0_wp       , &
                                             ocl_po4_acc = 0.0_wp               , &
                                             weath_po4_acc = 0.0_wp             , &
                                             biochem_min_po4_acc = 0.0_wp       , &
                                             desorp_po4_acc = 0.0_wp            , &
                                             plant_upt_po4_acc = 0.0_wp         , &
                                             microbial_upt_po4_acc = 0.0_wp     , &
                                             myco_upt_po4_acc = 0.0_wp          , &
                                             trans_solute_po4_acc = 0.0_wp      , &
                                             trans_assoc_fast_po4_acc = 0.0_wp  , &
                                             trans_assoc_slow_po4_acc = 0.0_wp  , &
                                             trans_occluded_po4_acc = 0.0_wp    , &
                                             trans_primary_po4_acc = 0.0_wp     , &
                                             transport_DOP_lay_acc = 0.0_wp
    LOGICAL                               :: l_flux_calc_acc_output, &  !< do variable accumulation and enable output
                                             l_flux_output_interval     !< do output at this timestep and set accumulation variables to zero
    REAL(wp), SAVE                        :: dtime_one = 1.0_wp         !< just 1.0    
    INTEGER                               :: iblk, ics, ice, nc, nsoil_sb
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':write_output_flux'
    ! ---------------------------
    ! 0.3 Declare Memory
    ! Declare process configuration and memory Pointers
    dsl4jsb_Def_config(SB_)
    dsl4jsb_Def_memory(A2L_)                   !     TYPE(t_a2l_memory),    POINTER :: a2l__mem
    dsl4jsb_Def_memory(L2A_)                   !     TYPE(t_l2a_memory),    POINTER :: l2a__mem
    dsl4jsb_Def_memory(ASSIMI_)                !     TYPE(t_assimi_memory), POINTER :: assimi__mem
    dsl4jsb_Def_memory(VEG_)                   !     TYPE(t_veg_memory),    POINTER :: veg__mem
    dsl4jsb_Def_memory(RAD_)                   !     TYPE(t_rad_memory),    POINTER :: rad__mem
    dsl4jsb_Def_memory(SOIL_)                  !     TYPE(t_soil_memory),   POINTER :: soil__mem
    dsl4jsb_Def_memory(SB_)                    !     TYPE(t_sb_memory),     POINTER :: sb__mem
    ! Declare pointers to variables in memory
    ! A2L_
    dsl4jsb_Real2D_onChunk      :: ga
    dsl4jsb_Real2D_onChunk      :: swpar_srf_down
    dsl4jsb_Real2D_onChunk      :: nhx_deposition
    dsl4jsb_Real2D_onChunk      :: noy_deposition
    dsl4jsb_Real2D_onChunk      :: nhx_n15_deposition
    dsl4jsb_Real2D_onChunk      :: noy_n15_deposition
    dsl4jsb_Real2D_onChunk      :: p_deposition
    dsl4jsb_Real2D_onChunk      :: rain 
    dsl4jsb_Real2D_onChunk      :: snow
    ! L2A_
    dsl4jsb_Real2D_onChunk      :: albedo
    dsl4jsb_Real2D_onChunk      :: net_radiation
    dsl4jsb_Real2D_onChunk      :: sensible_heat_flx
    dsl4jsb_Real2D_onChunk      :: latent_heat_flx
    dsl4jsb_Real2D_onChunk      :: ground_heat_flx
    ! ASSIMI_
    dsl4jsb_Real2D_onChunk      :: gross_assimilation
    dsl4jsb_Real2D_onChunk      :: gross_assimilation_C13
    dsl4jsb_Real2D_onChunk      :: gross_assimilation_C14
    dsl4jsb_Real2D_onChunk      :: canopy_cond
    ! SOIL_ 2D
    dsl4jsb_Real2D_onChunk      :: drainage
    dsl4jsb_Real2D_onChunk      :: runoff 
    dsl4jsb_Real2D_onChunk      :: evaporation 
    dsl4jsb_Real2D_onChunk      :: transpiration
    ! SOIL_ 3D
    dsl4jsb_Real3D_onChunk      :: soil_depth_sl
    ! RAD_
    dsl4jsb_Real2D_onChunk      :: appfd
    ! VEG_
    dsl4jsb_Real2D_onChunk      :: growth_respiration
    dsl4jsb_Real2D_onChunk      :: maint_respiration
    dsl4jsb_Real2D_onChunk      :: n_transform_respiration
    dsl4jsb_Real2D_onChunk      :: n_fixation_respiration
    dsl4jsb_Real2D_onChunk      :: npp
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
    dsl4jsb_Real2D_onChunk      :: interception 
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
    dsl4jsb_Real3D_onChunk      :: slow_sorption_po4
    dsl4jsb_Real3D_onChunk      :: fast_adsorption_po4
    dsl4jsb_Real3D_onChunk      :: biochem_mineralisation_po4
    dsl4jsb_Real3D_onChunk      :: plant_uptake_nh4_sl
    dsl4jsb_Real3D_onChunk      :: plant_uptake_no3_sl
    dsl4jsb_Real3D_onChunk      :: plant_uptake_po4_sl
    dsl4jsb_Real3D_onChunk      :: mycorrhiza_uptake_nh4_sl
    dsl4jsb_Real3D_onChunk      :: mycorrhiza_uptake_no3_sl
    dsl4jsb_Real3D_onChunk      :: mycorrhiza_uptake_po4_sl
    dsl4jsb_Real3D_onChunk      :: microbial_uptake_nh4_sl
    dsl4jsb_Real3D_onChunk      :: microbial_uptake_no3_sl
    dsl4jsb_Real3D_onChunk      :: microbial_uptake_po4_sl
    dsl4jsb_Real3D_onChunk      :: transport_po4_solute
    dsl4jsb_Real3D_onChunk      :: transport_po4_assoc_fast
    dsl4jsb_Real3D_onChunk      :: transport_po4_assoc_slow
    dsl4jsb_Real3D_onChunk      :: transport_po4_occluded
    dsl4jsb_Real3D_onChunk      :: transport_po4_primary
    ! Get local variables from options argument
    iblk     = options%iblk
    ics      = options%ics
    ice      = options%ice
    nc       = options%nc
    nsoil_sb = options%nsoil_sb
    ! ---------------------------
    ! 0.4 Process Activity, Debug Option
!    IF (.NOT. tile%Is_process_active(xx_)) RETURN   ! ask for all processes considered here ?
    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')
    ! ---------------------------
    ! 0.5 Get Memory
    model  => Get_model(tile%owner_model_id)
    lctlib => model%lctlib(tile%lcts%lib_id)
    ! Get process config
    dsl4jsb_Get_config(SB_)
    ! Get process memories
    dsl4jsb_Get_memory(A2L_)
    dsl4jsb_Get_memory(L2A_)
    dsl4jsb_Get_memory(ASSIMI_)
    dsl4jsb_Get_memory(VEG_)
    dsl4jsb_Get_memory(RAD_)
    dsl4jsb_Get_memory(SOIL_)
    dsl4jsb_Get_memory(SB_)
    ! Get process variables (Set pointers to variables in memory)
    dsl4jsb_Get_var2D_onChunk(A2L_, swpar_srf_down)
    dsl4jsb_Get_var2D_onChunk(A2L_, ga)
    dsl4jsb_Get_var2D_onChunk(A2L_, nhx_deposition)
    dsl4jsb_Get_var2D_onChunk(A2L_, noy_deposition)
    dsl4jsb_Get_var2D_onChunk(A2L_, nhx_n15_deposition)
    dsl4jsb_Get_var2D_onChunk(A2L_, noy_n15_deposition)
    dsl4jsb_Get_var2D_onChunk(A2L_, p_deposition)
    dsl4jsb_Get_var2D_onChunk(A2L_, rain)
    dsl4jsb_Get_var2D_onChunk(A2L_, snow) 
    ! ---------------------------
    dsl4jsb_Get_var2D_onChunk(L2A_, albedo)
    dsl4jsb_Get_var2D_onChunk(L2A_, net_radiation)
    dsl4jsb_Get_var2D_onChunk(L2A_, sensible_heat_flx)
    dsl4jsb_Get_var2D_onChunk(L2A_, latent_heat_flx)
    dsl4jsb_Get_var2D_onChunk(L2A_, ground_heat_flx)
    ! ---------------------------
    dsl4jsb_Get_var2D_onChunk(ASSIMI_, gross_assimilation)
    dsl4jsb_Get_var2D_onChunk(ASSIMI_, gross_assimilation_C13)
    dsl4jsb_Get_var2D_onChunk(ASSIMI_, gross_assimilation_C14)
    dsl4jsb_Get_var2D_onChunk(ASSIMI_, canopy_cond)
    ! ---------------------------
    dsl4jsb_Get_var2D_onChunk(SOIL_, evaporation)
    dsl4jsb_Get_var2D_onChunk(SOIL_, interception)
    dsl4jsb_Get_var2D_onChunk(SOIL_, runoff)
    dsl4jsb_Get_var2D_onChunk(SOIL_, drainage)
    dsl4jsb_Get_var3D_onChunk(SOIL_, soil_depth_sl)
    dsl4jsb_Get_var2D_onChunk(SOIL_, transpiration)
    ! ---------------------------
    dsl4jsb_Get_var2D_onChunk(RAD_, appfd)
    ! ---------------------------
    dsl4jsb_Get_var2D_onChunk(VEG_, growth_respiration)
    dsl4jsb_Get_var2D_onChunk(VEG_, maint_respiration)
    dsl4jsb_Get_var2D_onChunk(VEG_, n_transform_respiration)
    dsl4jsb_Get_var2D_onChunk(VEG_, n_fixation_respiration)
    dsl4jsb_Get_var2D_onChunk(VEG_, npp)
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
    ! ---------------------------
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
    dsl4jsb_Get_var3D_onChunk(SB_,  slow_sorption_po4)
    dsl4jsb_Get_var3D_onChunk(SB_,  fast_adsorption_po4)
    dsl4jsb_Get_var3D_onChunk(SB_,  weathering_po4)
    dsl4jsb_Get_var3D_onChunk(SB_,  occlusion_po4)
    dsl4jsb_Get_var3D_onChunk(SB_,  biochem_mineralisation_po4)
    dsl4jsb_Get_var3D_onChunk(SB_,  plant_uptake_nh4_sl)
    dsl4jsb_Get_var3D_onChunk(SB_,  plant_uptake_no3_sl)
    dsl4jsb_Get_var3D_onChunk(SB_,  plant_uptake_po4_sl)
    dsl4jsb_Get_var3D_onChunk(SB_,  mycorrhiza_uptake_nh4_sl)
    dsl4jsb_Get_var3D_onChunk(SB_,  mycorrhiza_uptake_no3_sl)
    dsl4jsb_Get_var3D_onChunk(SB_,  mycorrhiza_uptake_po4_sl)
    dsl4jsb_Get_var3D_onChunk(SB_,  microbial_uptake_nh4_sl)
    dsl4jsb_Get_var3D_onChunk(SB_,  microbial_uptake_no3_sl)
    dsl4jsb_Get_var3D_onChunk(SB_,  microbial_uptake_po4_sl)
    dsl4jsb_Get_var3D_onChunk(SB_,  transport_po4_solute)
    dsl4jsb_Get_var3D_onChunk(SB_,  transport_po4_assoc_fast)
    dsl4jsb_Get_var3D_onChunk(SB_,  transport_po4_assoc_slow)
    dsl4jsb_Get_var3D_onChunk(SB_,  transport_po4_occluded)
    dsl4jsb_Get_var3D_onChunk(SB_,  transport_po4_primary)

    
    ! Get process memory pool structures
#define Get_veg_litterfall
#define Get_veg_exudation
#define Get_veg_reserve_use
#define Get_sb_transport
#define Get_sb_mycorrhiza_export
#define Get_sb_delta_pools
#include "slm_pools_a_incl.f90"
#undef Get_veg_litterfall
#undef Get_veg_exudation
#undef Get_veg_reserve_use
#undef Get_sb_transport
#undef Get_sb_mycorrhiza_export
#undef Get_sb_delta_pools

    


    ! ------------------------------------------------------------------------------------------------------------
    ! Go science
    ! ------------------------------------------------------------------------------------------------------------


    !< 1.0 test if anything would be calculated or/and written out at this timestep
    !! 
    ! set to FALSE as default
    l_flux_calc_acc_output = .FALSE.
    l_flux_output_interval = .FALSE.
    !! test if TRUE
    ! after spinup, when doing transient simulations
    ! output start & end
    IF(.NOT. l_transient_spinup .AND. &
       & time_y >= base_config%output_start_first_day_year .AND. time_y <= base_config%output_end_last_day_year) THEN
      ! set TRUE
      l_flux_calc_acc_output = .TRUE.
    END IF
    ! special fluxnet output
    SELECT CASE(TRIM(base_config%slm_model_name))
      CASE('land','plant','canopy')
        IF (l_flux_calc_acc_output .AND. fluxnet_output_status < 3) THEN  ! if l_flux_calc_acc_output==TRUE and fluxnet_output_status < 3
          l_flux_calc_acc_output = .TRUE.
        ELSE
          l_flux_calc_acc_output = .FALSE.
        END IF
    END SELECT
    ! matching output_interval
    IF(MODULO(time_d, output_interval_flux) < eps4 .OR. &
       & ABS(output_interval_flux - MODULO(time_d, output_interval_flux)) < eps4) THEN
      ! set TRUE
      l_flux_output_interval = .TRUE.
    END IF
    ! accumulate fluxes across the timesteps of one interval:
    ! for timestep-interval leave it to seconds (which is the basic timestep for fluxes), i.e. multiply by 1 \n
    ! for any interval larger timestep multiply by dtime
    IF(output_interval_flux < 1.0_wp) THEN
      dtime_one = 1.0_wp                 ! interval = timestep
    ELSE
      dtime_one = dtime / 1000000.0_wp   ! interval > timestep
    END IF


    !< 2.0  output
    !! 
    !! land/plant/canopy
    !! land/soil  (simple_0d/simple_1d/jsm  &  simple_1d/jsm)
    !! 
    ! set starting-lun for write
    lun = 100 + 30
    SELECT CASE(TRIM(base_config%slm_model_name))
      CASE('land','plant','canopy')

        !! run this routine ONLY if 'special case: fluxnet output' is not active, or active but not finished
        !! same for write_output_pool
        IF (fluxnet_output_status < 3) THEN   ! this switch includes the rest of this routine
           !! accumulate fluxes across the timesteps of one interval
           !! 
           IF(l_flux_calc_acc_output) THEN
             ! Vegetation
             gpp_acc              = gpp_acc               + gross_assimilation(ipoi)        * dtime_one
             growth_resp_acc      = growth_resp_acc       + growth_respiration(ipoi)        * dtime_one
             maint_resp_acc       = maint_resp_acc        + maint_respiration(ipoi)         * dtime_one
             n_transform_resp_acc = n_transform_resp_acc  + n_transform_respiration(ipoi)   * dtime_one
             n_fixation_resp_acc  = n_fixation_resp_acc   + n_fixation_respiration(ipoi)    * dtime_one
             exudationC_acc       = exudationC_acc        + veg_exudation%carbon(ipoi)      * dtime_one *1000000._wp/dtime      
             npp_acc              = npp_acc               + npp(ipoi)                       * dtime_one
             res_use_c_acc        = res_use_c_acc         + veg_reserve_use%carbon(ipoi)    * dtime_one *1000000._wp/dtime
             gpp_C13_acc          = gpp_C13_acc           + gross_assimilation_C13(ipoi)    * dtime_one
             growth_resp_C13_acc  = growth_resp_C13_acc   + growth_respiration_c13(ipoi)    * dtime_one
             maint_resp_C13_acc   = maint_resp_C13_acc    + maint_respiration_c13(ipoi)     * dtime_one
             gpp_C14_acc          = gpp_C14_acc           + gross_assimilation_C14(ipoi)    * dtime_one
             growth_resp_C14_acc  = growth_resp_C14_acc   + growth_respiration_c14(ipoi)    * dtime_one
             maint_resp_C14_acc   = maint_resp_C14_acc    + maint_respiration_c14(ipoi)     * dtime_one
             uptake_nh4_acc       = uptake_nh4_acc        + uptake_nh4(ipoi)                * dtime_one
             uptake_nh4_n15_acc   = uptake_nh4_n15_acc    + uptake_nh4_n15(ipoi)            * dtime_one
             uptake_no3_acc       = uptake_no3_acc        + uptake_no3(ipoi)                * dtime_one
             uptake_no3_n15_acc   = uptake_no3_n15_acc    + uptake_no3_n15(ipoi)            * dtime_one
             recycling_n_acc      = recycling_n_acc       + ( recycling_leaf_n(ipoi) + recycling_fine_root_n(ipoi) + &
                                                                                       recycling_heart_wood_n(ipoi)) * &
                                                                                              dtime_one *1000000._wp/dtime
             res_use_n_acc        = res_use_n_acc         + veg_reserve_use%nitrogen(ipoi)  * dtime_one *1000000._wp/dtime
             n_fixation_acc       = n_fixation_acc        + n_fixation(ipoi)                * dtime_one 
             n_fixation_n15_acc   = n_fixation_n15_acc    + n_fixation_n15(ipoi)            * dtime_one 
             uptake_po4_acc       = uptake_po4_acc        + uptake_po4(ipoi)                * dtime_one
             recycling_p_acc      = recycling_p_acc       + ( recycling_leaf_p(ipoi) + recycling_fine_root_p(ipoi) + &
                                                                                       recycling_heart_wood_p(ipoi)) * &
                                                                                                dtime_one *1000000._wp/dtime
             res_use_p_acc        = res_use_p_acc         + veg_reserve_use%phosphorus(ipoi)  * dtime_one *1000000._wp/dtime
             IF(output_interval_flux < 1.0_wp) THEN
                precip_acc       = precip_acc        + (rain(ipoi) + snow(ipoi))  
                interception_acc = interception_acc  + interception(ipoi)       
                evaporation_acc  = evaporation_acc   + evaporation(ipoi)       
                transpiration_acc= transpiration_acc + transpiration(ipoi)      
                runoff_acc       = runoff_acc        + runoff(ipoi)              
                drainage_acc     = drainage_acc      + drainage(ipoi)             
             ELSE
                precip_acc       = precip_acc        + (rain(ipoi) + snow(ipoi))          * dtime
                interception_acc = interception_acc  + interception(ipoi)                 * dtime 
                evaporation_acc  = evaporation_acc   + evaporation(ipoi)                  * dtime
                transpiration_acc= transpiration_acc + transpiration(ipoi)                * dtime
                runoff_acc       = runoff_acc        + runoff(ipoi)                       * dtime
                drainage_acc     = drainage_acc      + drainage(ipoi)                     * dtime
             ENDIF 
             gc_acc           = gc_acc            + canopy_cond(ipoi)                                                      
             ga_acc           = ga_acc            + ga(ipoi) 
             albedo_acc       = albedo_acc        + albedo(ipoi)
             ppfd_acc         = ppfd_acc          + swpar_srf_down(ipoi) * rad2ppfd  * dtime_one               
             rnet_acc         = rnet_acc          + net_radiation(ipoi)              * dtime_one               
             qh_acc           = qh_acc            - sensible_heat_flx(ipoi)          * dtime_one               
             qle_acc          = qle_acc           - latent_heat_flx(ipoi)            * dtime_one               
             qg_acc           = qg_acc            - ground_heat_flx(ipoi)            * dtime_one               
             appfd_acc        = appfd_acc         + appfd(ipoi)                      * dtime_one               
             ndat_acc         = ndat_acc          + 1 
             ! Mycorrhiza
             myco_export_c_acc   = myco_export_c_acc    + SUM(sb_mycorrhiza_export%carbon(ipoi,:) * &
                                                          soil_depth_sl(ipoi,:)) * dtime_one *1000000._wp/dtime
             myco_export_n_acc   = myco_export_n_acc    + SUM(sb_mycorrhiza_export%nitrogen(ipoi,:) * &
                                                          soil_depth_sl(ipoi,:)) * dtime_one *1000000._wp/dtime
             myco_export_p_acc   = myco_export_p_acc    + SUM(sb_mycorrhiza_export%phosphorus(ipoi,:) * &
                                                          soil_depth_sl(ipoi,:)) * dtime_one *1000000._wp/dtime
             myco_export_c13_acc = myco_export_c13_acc  + SUM(sb_mycorrhiza_export%carbon13(ipoi,:) * &
                                                          soil_depth_sl(ipoi,:)) * dtime_one *1000000._wp/dtime
             myco_export_c14_acc = myco_export_c14_acc  + SUM(sb_mycorrhiza_export%carbon14(ipoi,:) * &
                                                          soil_depth_sl(ipoi,:)) * dtime_one *1000000._wp/dtime
             myco_export_n15_acc = myco_export_n15_acc  + SUM(sb_mycorrhiza_export%nitrogen15(ipoi,:) * &
                                                          soil_depth_sl(ipoi,:)) * dtime_one *1000000._wp/dtime

             !> 1.2 write fluxes and set to zero
             !!
             IF(l_flux_output_interval) THEN

               IF(ndat_acc > 1) THEN
                 gc_acc = gc_acc / ndat_acc
                 ga_acc = ga_acc / ndat_acc
                 albedo_acc = albedo_acc / ndat_acc
               ENDIF
       
               ! land
               SELECT CASE(TRIM(base_config%slm_model_name))  ! cases: land & land/plant & canopy
                 CASE('land')
                   ! mycofluxCNP
                   write(lun,'(7F24.12)') time_d, & 
                                          myco_export_c_acc, myco_export_n_acc, myco_export_p_acc, &
                                          myco_export_c13_acc, myco_export_c14_acc, myco_export_n15_acc
                   lun=lun+1
               END SELECT

               ! land/plant & canopy & soil
               SELECT CASE(TRIM(base_config%slm_model_name))
                 CASE('land','plant') ! additional to CASE('land')
                   ! vegfluxC
                   write(lun,'(10F18.6)') time_d, gpp_acc, growth_resp_acc, maint_resp_acc , &
                                          n_transform_resp_acc,n_fixation_resp_acc         , & 
                                          npp_acc, res_use_c_acc, exudationC_acc
                   lun=lun+1
                   ! vegfluxdC13DC14
                   dC13gpp = calc_delta_C13(gpp_acc,gpp_C13_acc)
                   dC13rg  = calc_delta_C13(growth_resp_acc,growth_resp_C13_acc)
                   dC13rm  = calc_delta_C13(maint_resp_acc,maint_resp_C13_acc) 
                   write(lun,'(7F18.6)') time_d, dC13gpp, dC13rg, dC13rm                            , & 
                                         calc_delta_C14(gpp_acc,gpp_C14_acc,dC13gpp)                , & 
                                         calc_delta_C14(growth_resp_acc,growth_resp_C14_acc,dC13rg) , & 
                                         calc_delta_C14(maint_resp_acc,maint_resp_C14_acc,dC13rm)
                   lun=lun+1
                   ! vegfluxN
                   write(lun,'(10F18.6)') time_d, uptake_nh4_acc, uptake_no3_acc, recycling_n_acc   , & 
                                          res_use_n_acc, n_fixation_acc, & 
                                          calc_delta_N15(n_fixation_acc,n_fixation_n15_acc), &
                                          calc_delta_N15(uptake_nh4_acc,uptake_nh4_n15_acc), &
                                          calc_delta_N15(uptake_no3_acc,uptake_no3_n15_acc)

                   lun=lun+1
                   ! vegfluxP
                   write(lun,'(4F18.6)') time_d, uptake_po4_acc, recycling_p_acc, res_use_p_acc
                   lun=lun+1
                   write(lun,'(12F18.8)') time_d, precip_acc, evaporation_acc, interception_acc, & 
                                          transpiration_acc, runoff_acc, drainage_acc, gc_acc, ga_acc
                   lun=lun+1
                   IF(ppfd_acc>eps8) THEN
                      hlp1=appfd_acc/ppfd_acc
                   ELSE
                      hlp1=0.0_wp
                   ENDIF
                   write(lun,'(12F18.8)') time_d, rnet_acc, qh_acc, qle_acc, -qg_acc, & 
                                          albedo_acc, ppfd_acc, appfd_acc, hlp1
                   lun=lun+1
                 CASE('canopy') ! subset of the previous CASE('land','plant')
                   ! files: vegfluxC, vegfluxdC13dC14, vegfluxH2O, vegfluxEnergy
                   write(lun,'(8F18.6)') time_d, gpp_acc, growth_resp_acc, maint_resp_acc , &
                                         n_transform_resp_acc,n_fixation_resp_acc         , & 
                                         npp_acc, res_use_c_acc
                   lun=lun+1
                   dC13gpp = calc_delta_C13(gpp_acc,gpp_C13_acc)
                   dC13rg = calc_delta_C13(growth_resp_acc,growth_resp_C13_acc)
                   dC13rm = calc_delta_C13(maint_resp_acc,maint_resp_C13_acc) 
                   write(lun,'(7F18.6)') time_d, dC13gpp, dC13rg, dC13rm                            , & 
                                         calc_delta_C14(gpp_acc,gpp_C14_acc,dC13gpp)                , & 
                                         calc_delta_C14(growth_resp_acc,growth_resp_C14_acc,dC13rg) , & 
                                         calc_delta_C14(maint_resp_acc,maint_resp_C14_acc,dC13rm)
                   lun=lun+1
                   write(lun,'(12F18.8)') time_d, precip_acc, evaporation_acc, interception_acc, transpiration_acc, runoff_acc, & 
                                          drainage_acc, gc_acc, ga_acc
                   lun=lun+1
                   IF(ppfd_acc>eps8) THEN
                      hlp1=appfd_acc/ppfd_acc
                   ELSE
                      hlp1=0.0_wp
                   ENDIF
                   write(lun,'(12F18.8)') time_d, rnet_acc, qh_acc, qle_acc, -qg_acc, & 
                                          albedo_acc, ppfd_acc, appfd_acc, hlp1
                   lun=lun+1
               END SELECT ! land/plant & canopy

               !! set acc variables zero
               gpp_acc                    = 0.0_wp
               growth_resp_acc            = 0.0_wp
               maint_resp_acc             = 0.0_wp
               n_transform_resp_acc       = 0.0_wp
               n_fixation_resp_acc        = 0.0_wp
               exudationC_acc             = 0.0_wp
               npp_acc                    = 0.0_wp
               res_use_c_acc              = 0.0_wp
               gpp_C13_acc                = 0.0_wp
               growth_resp_C13_acc        = 0.0_wp
               maint_resp_C13_acc         = 0.0_wp
               gpp_C14_acc                = 0.0_wp
               growth_resp_C14_acc        = 0.0_wp
               maint_resp_C14_acc         = 0.0_wp
               recycling_n_acc            = 0.0_wp
               res_use_n_acc              = 0.0_wp
               n_fixation_acc             = 0.0_wp
               n_fixation_n15_acc         = 0.0_wp
               uptake_po4_acc             = 0.0_wp
               recycling_p_acc            = 0.0_wp
               res_use_p_acc              = 0.0_wp
               precip_acc                 = 0.0_wp
               evaporation_acc            = 0.0_wp
               interception_acc           = 0.0_wp
               transpiration_acc          = 0.0_wp
               drainage_acc               = 0.0_wp
               runoff_acc                 = 0.0_wp
               gc_acc                     = 0.0_wp
               ga_acc                     = 0.0_wp
               ppfd_acc                   = 0.0_wp
               rnet_acc                   = 0.0_wp
               qh_acc                     = 0.0_wp
               qle_acc                    = 0.0_wp
               qg_acc                     = 0.0_wp
               albedo_acc                 = 0.0_wp
               appfd_acc                  = 0.0_wp
               ndat_acc                   = 0        ! INTEGER 
               uptake_no3_acc             = 0.0_wp
               uptake_no3_n15_acc         = 0.0_wp
               uptake_nh4_acc             = 0.0_wp
               uptake_nh4_n15_acc         = 0.0_wp
               myco_export_c_acc          = 0.0_wp
               myco_export_n_acc          = 0.0_wp
               myco_export_p_acc          = 0.0_wp
               myco_export_c13_acc        = 0.0_wp
               myco_export_c14_acc        = 0.0_wp
               myco_export_n15_acc        = 0.0_wp
             END IF ! l_flux_output_interval
           END IF   ! l_flux_calc_acc_output
        END IF      ! on/off switch specific to 'special case: fluxnet output'
    END SELECT      ! land/plant/canopy


    !! land & soil - additional to CASE('land')
    !!  simple_0d and simple_1d and jsm
    !!  layer-data for simple_1d and jsm
    SELECT CASE(TRIM(base_config%slm_model_name))
    CASE('land','soil')
      IF (l_flux_calc_acc_output) THEN
        !! all soil models (simple_0d and simple_1d and jsm)
        ! accumulate variables
        ! Carbon
        litterfallC_acc = litterfallC_acc + (veg_litterfall%leaf%carbon(ipoi) + & 
                                veg_litterfall%fine_root%carbon(ipoi) + veg_litterfall%coarse_root%carbon(ipoi) + &
                                veg_litterfall%sap_wood%carbon(ipoi) + veg_litterfall%heart_wood%carbon(ipoi) + &
                                veg_litterfall%labile%carbon(ipoi) + veg_litterfall%reserve%carbon(ipoi) + &
                                veg_litterfall%fruit%carbon(ipoi) + veg_litterfall%seed_bed%carbon(ipoi)) * &
                                1000000.0_wp * dtime_one/dtime
        litterfallC13_acc = litterfallC13_acc + (veg_litterfall%leaf%carbon13(ipoi) + & 
                                veg_litterfall%fine_root%carbon13(ipoi) + veg_litterfall%coarse_root%carbon13(ipoi) + &
                                veg_litterfall%sap_wood%carbon13(ipoi) + veg_litterfall%heart_wood%carbon13(ipoi) + &
                                veg_litterfall%labile%carbon13(ipoi) + veg_litterfall%reserve%carbon13(ipoi) + &
                                veg_litterfall%fruit%carbon13(ipoi) + veg_litterfall%seed_bed%carbon13(ipoi)) * &
                                1000000.0_wp * dtime_one/dtime
        litterfallC14_acc = litterfallC14_acc + (veg_litterfall%leaf%carbon14(ipoi) + & 
                                veg_litterfall%fine_root%carbon14(ipoi) + veg_litterfall%coarse_root%carbon14(ipoi) + &
                                veg_litterfall%sap_wood%carbon14(ipoi) + veg_litterfall%heart_wood%carbon14(ipoi) + &
                                veg_litterfall%labile%carbon14(ipoi) + veg_litterfall%reserve%carbon14(ipoi) + &
                                veg_litterfall%fruit%carbon14(ipoi) + veg_litterfall%seed_bed%carbon14(ipoi)) * &
                                1000000.0_wp * dtime_one/dtime
        litterfallN_acc = litterfallN_acc + (veg_litterfall%leaf%nitrogen(ipoi) + & 
                                veg_litterfall%fine_root%nitrogen(ipoi) + veg_litterfall%coarse_root%nitrogen(ipoi) + &
                                veg_litterfall%sap_wood%nitrogen(ipoi) + veg_litterfall%heart_wood%nitrogen(ipoi) + &
                                veg_litterfall%labile%nitrogen(ipoi) + veg_litterfall%reserve%nitrogen(ipoi) + &
                                veg_litterfall%fruit%nitrogen(ipoi) + veg_litterfall%seed_bed%nitrogen(ipoi)) * &
                                1000000.0_wp * dtime_one/dtime
        litterfallP_acc = litterfallP_acc + (veg_litterfall%leaf%phosphorus(ipoi) + & 
                                veg_litterfall%fine_root%phosphorus(ipoi) + veg_litterfall%coarse_root%phosphorus(ipoi) + &
                                veg_litterfall%sap_wood%phosphorus(ipoi) + veg_litterfall%heart_wood%phosphorus(ipoi) + &
                                veg_litterfall%labile%phosphorus(ipoi) + veg_litterfall%reserve%phosphorus(ipoi) + &
                                veg_litterfall%fruit%phosphorus(ipoi) + veg_litterfall%seed_bed%phosphorus(ipoi)) * &
                                1000000.0_wp * dtime_one/dtime

        het_respC_acc      = het_respC_acc      + SUM(het_respiration(ipoi,:) * soil_depth_sl(ipoi,:)) * dtime_one
        het_respC13_acc    = het_respC13_acc    + SUM(het_respiration_c13(ipoi,:) * soil_depth_sl(ipoi,:)) * dtime_one
        het_respC14_acc    = het_respC14_acc    + SUM(het_respiration_c14(ipoi,:) * soil_depth_sl(ipoi,:)) * dtime_one
        leachingCorg_acc   = leachingCorg_acc   + leaching_dom_carbon(ipoi)   * dtime_one
        leachingC13org_acc = leachingC13org_acc + leaching_dom_carbon13(ipoi) * dtime_one
        leachingC14org_acc = leachingC14org_acc + leaching_dom_carbon14(ipoi) * dtime_one
        ! Nitrogen
        n_deposition_acc     = n_deposition_acc  + (nhx_deposition(ipoi) + noy_deposition(ipoi)) * dtime_one      ! unit: mumol/m2/sec * dtime_one
        n_deposition_n15_acc = n_deposition_n15_acc  + (nhx_n15_deposition(ipoi) + noy_n15_deposition(ipoi)) * dtime_one      ! unit: mumol/m2/sec * dtime_one
        net_min_nh4_acc      = net_min_nh4_acc + SUM(net_mineralisation_nh4(ipoi,:) * soil_depth_sl(ipoi,:)) * dtime_one
        net_min_no3_acc      = net_min_no3_acc  + SUM(net_mineralisation_no3(ipoi,:) * soil_depth_sl(ipoi,:)) * dtime_one                     
        plant_uptake_nh4_acc = plant_uptake_nh4_acc + &
                             SUM(plant_uptake_nh4_sl(ipoi,:) * soil_depth_sl(ipoi,:)) * dtime_one 
        plant_uptake_no3_acc = plant_uptake_no3_acc + &
                             SUM(plant_uptake_no3_sl(ipoi,:) * soil_depth_sl(ipoi,:)) * dtime_one
        nitrification_acc    = nitrification_acc + & 
                              SUM((nitrification_no3(ipoi,:) + nitrification_noy(ipoi,:) + nitrification_n2o(ipoi,:)) * &
                              soil_depth_sl(ipoi,:)) * dtime_one
        denitrification_acc  = denitrification_acc + & 
                              SUM((denitrification_n2(ipoi,:) + denitrification_noy(ipoi,:) + denitrification_n2o(ipoi,:))*&
                              soil_depth_sl(ipoi,:)) * dtime_one 
        emission_noy_acc = emission_noy_acc + emission_noy(ipoi) * dtime_one
        emission_noy_n15_acc = emission_noy_n15_acc + emission_noy_n15(ipoi) * dtime_one
        emission_n2o_acc = emission_n2o_acc + emission_n2o(ipoi) * dtime_one
        emission_n2o_n15_acc = emission_n2o_n15_acc + emission_n2o_n15(ipoi) * dtime_one
        emission_n2_acc  = emission_n2_acc  + emission_n2(ipoi) * dtime_one
        emission_n2_n15_acc  = emission_n2_n15_acc  + emission_n2_n15(ipoi) * dtime_one
        leaching_nh4_acc = leaching_nh4_acc + leaching_nh4_solute(ipoi) * dtime_one
        leaching_nh4_n15_acc = leaching_nh4_n15_acc + leaching_nh4_n15_solute(ipoi) * dtime_one
        leaching_no3_acc = leaching_no3_acc + leaching_no3_solute(ipoi) * dtime_one
        leaching_no3_n15_acc = leaching_no3_n15_acc + leaching_no3_n15_solute(ipoi) * dtime_one
        leachingNorg_acc = leachingNorg_acc + leaching_dom_nitrogen(ipoi) * dtime_one
        asymb_n_fixation_acc =  asymb_n_fixation_acc + SUM(asymb_n_fixation(ipoi,:)*soil_depth_sl(ipoi,:)) * dtime_one
        asymb_n_fixation_n15_acc =  asymb_n_fixation_n15_acc + SUM(asymb_n_fixation_n15(ipoi,:)*soil_depth_sl(ipoi,:)) * &
                                    dtime_one
        ! Phosphorus
        p_deposition_acc       = p_deposition_acc + p_deposition(ipoi) * dtime_one
        leaching_po4_acc       = leaching_po4_acc + leaching_po4_solute(ipoi) * dtime_one           
        leachingPorg_acc       = leachingPorg_acc + leaching_dom_phosphorus(ipoi) * dtime_one
        net_min_po4_acc        = net_min_po4_acc + &
                                  SUM(net_mineralisation_po4(ipoi,:) * soil_depth_sl(ipoi,:)) * dtime_one                
        plant_uptake_po4_acc   = plant_uptake_po4_acc + &
                                SUM(plant_uptake_po4_sl(ipoi,:) * soil_depth_sl(ipoi,:)) * dtime_one
        weathering_po4_acc     = weathering_po4_acc + &
                                  SUM(weathering_po4(ipoi,:) * soil_depth_sl(ipoi,:)) * dtime_one
        occlusion_po4_acc      = occlusion_po4_acc + &
                                  SUM(occlusion_po4(ipoi,:) * soil_depth_sl(ipoi,:)) * dtime_one
        slow_sorption_po4_acc  = slow_sorption_po4_acc + &
                                  SUM(slow_sorption_po4(ipoi,:) * soil_depth_sl(ipoi,:)) * dtime_one
        fast_adsorption_po4_acc        = fast_adsorption_po4_acc + &
                                          SUM(fast_adsorption_po4(ipoi,:) * soil_depth_sl(ipoi,:)) * dtime_one
        biochem_mineralisation_po4_acc = biochem_mineralisation_po4_acc + &
                                        SUM(biochem_mineralisation_po4(ipoi,:) * soil_depth_sl(ipoi,:)) * dtime_one
        ! mass balance checks in SB_
        sb_delta_pools_carbon_acc     = sb_delta_pools_carbon_acc     + &
                                          SUM(sb_delta_pools%carbon(ipoi,:)*soil_depth_sl(ipoi,:))
        sb_delta_pools_nitrogen_acc   = sb_delta_pools_nitrogen_acc   + &
                                          SUM(sb_delta_pools%nitrogen(ipoi,:)*soil_depth_sl(ipoi,:))
        sb_delta_pools_phosphorus_acc = sb_delta_pools_phosphorus_acc + &
                                          SUM(sb_delta_pools%phosphorus(ipoi,:)*soil_depth_sl(ipoi,:))
        sb_delta_pools_carbon13_acc   = sb_delta_pools_carbon13_acc   + &
                                          SUM(sb_delta_pools%carbon13(ipoi,:)*soil_depth_sl(ipoi,:))
        sb_delta_pools_carbon14_acc   = sb_delta_pools_carbon14_acc   + &
                                          SUM(sb_delta_pools%carbon14(ipoi,:)*soil_depth_sl(ipoi,:))
        sb_delta_pools_nitrogen15_acc = sb_delta_pools_nitrogen15_acc + &
                                          SUM(sb_delta_pools%nitrogen15(ipoi,:)*soil_depth_sl(ipoi,:))
        sb_total_flux_carbon_acc      = sb_total_flux_carbon_acc      + total_flux_carbon(ipoi)
        sb_total_flux_nitrogen_acc    = sb_total_flux_nitrogen_acc    + total_flux_nitrogen(ipoi)
        sb_total_flux_phosphorus_acc  = sb_total_flux_phosphorus_acc  + total_flux_phosphorus(ipoi)
        sb_total_flux_carbon13_acc    = sb_total_flux_carbon13_acc    + total_flux_carbon13(ipoi)
        sb_total_flux_carbon14_acc    = sb_total_flux_carbon14_acc    + total_flux_carbon14(ipoi)
        sb_total_flux_nitrogen15_acc  = sb_total_flux_nitrogen15_acc  + total_flux_nitrogen15(ipoi)

        ! write output & zero acc variables
        IF(l_flux_output_interval) THEN
          ! soilfluxCdC13DC14
          dC13litterfall = calc_delta_C13(litterfallC_acc,litterfallC13_acc)
          dC13het_resp   = calc_delta_C13(het_respC_acc,het_respC13_acc)
          dC13leaching   = calc_delta_C13(leachingCorg_acc,leachingC13org_acc)
          write(lun,'(12F20.12)') time_d, litterfallC_acc, het_respC_acc, leachingCorg_acc, &
             dC13litterfall, dC13het_resp, dC13leaching, &
             calc_delta_C14(litterfallC_acc,litterfallC14_acc,dC13litterfall), &
             calc_delta_C14(het_respC_acc,het_respC14_acc,dC13het_resp), &
             calc_delta_C14(leachingCorg_acc,leachingC14org_acc,dC13leaching)
          lun = lun + 1
          ! soilfluxN
          write(lun,'(20F20.12)') time_d,litterfallN_acc,n_deposition_acc,asymb_n_fixation_acc,net_min_nh4_acc, net_min_no3_acc, &
                                nitrification_acc,denitrification_acc,emission_noy_acc,emission_n2o_acc,emission_n2_acc,& 
                                plant_uptake_nh4_acc,plant_uptake_no3_acc,leaching_nh4_acc,leaching_no3_acc,leachingNorg_acc                    
          lun=lun+1
          ! soilfluxP
          write(lun,'(20F20.12)') time_d, litterfallP_acc, p_deposition_acc, leaching_po4_acc, leachingPorg_acc, &
                                 net_min_po4_acc, plant_uptake_po4_acc, weathering_po4_acc, occlusion_po4_acc, &
                                 slow_sorption_po4_acc, fast_adsorption_po4_acc, biochem_mineralisation_po4_acc
          lun = lun + 1 
          ! soilfluxdN
          write(lun,'(20F18.8)') time_d, & 
               calc_delta_N15(n_deposition_acc,n_deposition_n15_acc), & 
               calc_delta_N15(asymb_n_fixation_acc,asymb_n_fixation_n15_acc), & 
               calc_delta_N15(emission_noy_acc,emission_noy_n15_acc), &
               calc_delta_N15(emission_n2o_acc,emission_n2o_n15_acc), &
               calc_delta_N15(emission_n2_acc,emission_n2_n15_acc), &
               calc_delta_N15(leaching_nh4_acc,leaching_nh4_n15_acc), & 
               calc_delta_N15(leaching_no3_acc,leaching_no3_n15_acc)
          lun=lun+1
          ! sb_massBalance
          write(lun,'(13F24.12)') time_d, &
                                  sb_delta_pools_carbon_acc, sb_delta_pools_nitrogen_acc, sb_delta_pools_phosphorus_acc, &
                                  sb_delta_pools_carbon13_acc, sb_delta_pools_carbon14_acc, sb_delta_pools_nitrogen15_acc, &
                                  sb_total_flux_carbon_acc, sb_total_flux_nitrogen_acc, sb_total_flux_phosphorus_acc, &
                                  sb_total_flux_carbon13_acc, sb_total_flux_carbon14_acc, sb_total_flux_nitrogen15_acc
          lun=lun+1

          !! set acc variables to zero      
          litterfallC_acc                = 0.0_wp
          het_respC_acc                  = 0.0_wp
          leachingCorg_acc               = 0.0_wp
          litterfallC13_acc              = 0.0_wp
          het_respC13_acc                = 0.0_wp
          leachingC13org_acc             = 0.0_wp
          litterfallC14_acc              = 0.0_wp
          het_respC14_acc                = 0.0_wp
          leachingC14org_acc             = 0.0_wp
          nitrification_acc              = 0.0_wp
          denitrification_acc            = 0.0_wp
          litterfallN_acc                = 0.0_wp
          net_min_nh4_acc                = 0.0_wp
          net_min_no3_acc                = 0.0_wp
          plant_uptake_nh4_acc           = 0.0_wp
          plant_uptake_no3_acc           = 0.0_wp
          leachingNorg_acc               = 0.0_wp
          emission_noy_acc               = 0.0_wp
          emission_noy_n15_acc           = 0.0_wp
          emission_n2o_acc               = 0.0_wp
          emission_n2o_n15_acc           = 0.0_wp
          emission_n2_acc                = 0.0_wp
          emission_n2_n15_acc            = 0.0_wp
          leaching_nh4_acc               = 0.0_wp
          leaching_nh4_n15_acc           = 0.0_wp
          leaching_no3_acc               = 0.0_wp
          leaching_no3_n15_acc           = 0.0_wp
          asymb_n_fixation_acc           = 0.0_wp
          asymb_n_fixation_n15_acc       = 0.0_wp
          litterfallP_acc                = 0.0_wp
          leachingPorg_acc               = 0.0_wp
          p_deposition_acc               = 0.0_wp
          n_deposition_acc               = 0.0_wp
          n_deposition_n15_acc           = 0.0_wp
          leaching_po4_acc               = 0.0_wp
          net_min_po4_acc                = 0.0_wp
          plant_uptake_po4_acc           = 0.0_wp
          weathering_po4_acc             = 0.0_wp
          occlusion_po4_acc              = 0.0_wp
          slow_sorption_po4_acc          = 0.0_wp
          fast_adsorption_po4_acc        = 0.0_wp
          biochem_mineralisation_po4_acc = 0.0_wp
          sb_delta_pools_carbon_acc      = 0.0_wp
          sb_delta_pools_nitrogen_acc    = 0.0_wp
          sb_delta_pools_phosphorus_acc  = 0.0_wp
          sb_delta_pools_carbon13_acc    = 0.0_wp
          sb_delta_pools_carbon14_acc    = 0.0_wp
          sb_delta_pools_nitrogen15_acc  = 0.0_wp
          sb_total_flux_carbon_acc       = 0.0_wp
          sb_total_flux_nitrogen_acc     = 0.0_wp
          sb_total_flux_phosphorus_acc   = 0.0_wp
          sb_total_flux_carbon13_acc     = 0.0_wp
          sb_total_flux_carbon14_acc     = 0.0_wp
          sb_total_flux_nitrogen15_acc   = 0.0_wp
        END IF

        !! simple_1d & jsm layer-data
        SELECT CASE(TRIM(dsl4jsb_Config(SB_)%sb_model_scheme))
        CASE('simple_1d','jsm')
          ! accumulate variables
          ! C
          het_respiration_lay_acc(1:nsoil_sb)     = het_respiration_lay_acc(1:nsoil_sb) + &
                                                    het_respiration(ipoi,:) * dtime_one
          het_respiration_c13_lay_acc(1:nsoil_sb) = het_respiration_c13_lay_acc(1:nsoil_sb) + &
                                                    het_respiration_c13(ipoi,:) * dtime_one
          het_respiration_c14_lay_acc(1:nsoil_sb) = het_respiration_c14_lay_acc(1:nsoil_sb) + &
                                                    het_respiration_c14(ipoi,:) * dtime_one
          transport_DOC_lay_acc(1:nsoil_sb)       = transport_DOC_lay_acc(1:nsoil_sb) + &
                                                    sb_transport%dom%carbon(ipoi,:) * &
                                                    1000000.0_wp * dtime_one/dtime
          ! N
          net_mineral_nh4_lay_acc(1:nsoil_sb)     = net_mineral_nh4_lay_acc(1:nsoil_sb) + &
                                                    net_mineralisation_nh4(ipoi,:) * dtime_one
          plant_uptake_nh4_lay_acc(1:nsoil_sb)    = plant_uptake_nh4_lay_acc(1:nsoil_sb) + &
                                                    plant_uptake_nh4_sl(ipoi,:) * dtime_one
          microbial_upt_nh4_lay_acc(1:nsoil_sb)   = microbial_upt_nh4_lay_acc(1:nsoil_sb) + &
                                                    microbial_uptake_nh4_sl(ipoi,:) * dtime_one
          microbial_upt_no3_lay_acc(1:nsoil_sb)   = microbial_upt_no3_lay_acc(1:nsoil_sb) + &
                                                    microbial_uptake_no3_sl(ipoi,:) * dtime_one
          myco_uptake_nh4_lay_acc(1:nsoil_sb)     = myco_uptake_nh4_lay_acc(1:nsoil_sb) + &
                                                    mycorrhiza_uptake_nh4_sl(ipoi,:) * dtime_one
          asymb_n_fixation_lay_acc(1:nsoil_sb)    = asymb_n_fixation_lay_acc(1:nsoil_sb) + &
                                                    asymb_n_fixation(ipoi,:) * dtime_one
          net_mineral_no3_lay_acc(1:nsoil_sb)     = net_mineral_no3_lay_acc(1:nsoil_sb) + &
                                                    net_mineralisation_no3(ipoi,:) * dtime_one
          plant_uptake_no3_lay_acc(1:nsoil_sb)    = plant_uptake_no3_lay_acc(1:nsoil_sb) + &
                                                    plant_uptake_no3_sl(ipoi,:) * dtime_one
          myco_uptake_no3_lay_acc(1:nsoil_sb)     = myco_uptake_no3_lay_acc(1:nsoil_sb) + &
                                                    mycorrhiza_uptake_no3_sl(ipoi,:) * dtime_one
          transport_DON_lay_acc(1:nsoil_sb)       = transport_DON_lay_acc(1:nsoil_sb) + &
                                                    sb_transport%dom%nitrogen(ipoi,:) * &
                                                    1000000.0_wp * dtime_one/dtime
          ! P
          adsorpt_po4_acc(1:nsoil_sb)          = adsorpt_po4_acc(1:nsoil_sb) + &
                                                 fast_adsorption_po4(ipoi,:) * dtime_one
          net_min_po4_lay_acc(1:nsoil_sb)      = net_min_po4_lay_acc(1:nsoil_sb) + &
                                                 net_mineralisation_po4(ipoi,:) * dtime_one
          weath_po4_acc(1:nsoil_sb)            = weath_po4_acc(1:nsoil_sb) + &
                                                 weathering_po4(ipoi,:) * dtime_one
          ocl_po4_acc(1:nsoil_sb)              = ocl_po4_acc(1:nsoil_sb) + &
                                                 occlusion_po4(ipoi,:) * dtime_one
          desorp_po4_acc(1:nsoil_sb)           = desorp_po4_acc(1:nsoil_sb) + &
                                                 slow_sorption_po4(ipoi,:) * dtime_one
          biochem_min_po4_acc(1:nsoil_sb)      = biochem_min_po4_acc(1:nsoil_sb) + &
                                                 biochem_mineralisation_po4(ipoi,:) * dtime_one
          plant_upt_po4_acc(1:nsoil_sb)        = plant_upt_po4_acc(1:nsoil_sb) + &
                                                 plant_uptake_po4_sl(ipoi,:) * dtime_one
          microbial_upt_po4_acc(1:nsoil_sb)    = microbial_upt_po4_acc(1:nsoil_sb) + &
                                                 microbial_uptake_po4_sl(ipoi,:) * dtime_one
          myco_upt_po4_acc(1:nsoil_sb)         = myco_upt_po4_acc(1:nsoil_sb) + &
                                                 mycorrhiza_uptake_po4_sl(ipoi,:) * dtime_one
          trans_solute_po4_acc(1:nsoil_sb)     = trans_solute_po4_acc(1:nsoil_sb) + &
                                                 transport_po4_solute(ipoi,:) * dtime_one
          trans_assoc_fast_po4_acc(1:nsoil_sb) = trans_assoc_fast_po4_acc(1:nsoil_sb) + &
                                                 transport_po4_assoc_fast(ipoi,:) * dtime_one
          trans_assoc_slow_po4_acc(1:nsoil_sb) = trans_assoc_slow_po4_acc(1:nsoil_sb) + &
                                                 transport_po4_assoc_slow(ipoi,:) * dtime_one
          trans_occluded_po4_acc(1:nsoil_sb)   = trans_occluded_po4_acc(1:nsoil_sb) + &
                                                 transport_po4_occluded(ipoi,:) * dtime_one
          trans_primary_po4_acc(1:nsoil_sb)    = trans_primary_po4_acc(1:nsoil_sb) + &
                                                 transport_po4_primary(ipoi,:) * dtime_one
          transport_DOP_lay_acc(1:nsoil_sb)    = transport_DOP_lay_acc(1:nsoil_sb) + &
                                                 sb_transport%dom%phosphorus(ipoi,:) * &
                                                 1000000.0_wp * dtime_one/dtime
          ! write output & zero acc variables
          IF(l_flux_output_interval) THEN

            ! soilfluxC_lay
            write(lun,'(101F20.10)') time_d, &
              het_respiration_lay_acc(1:nsoil_sb), & 
              calc_delta_C13(het_respiration_lay_acc(1:nsoil_sb),het_respiration_c13_lay_acc(1:nsoil_sb)), & 
              calc_delta_C14(het_respiration_lay_acc(1:nsoil_sb),het_respiration_c14_lay_acc(1:nsoil_sb), &
              calc_delta_C13(het_respiration_lay_acc(1:nsoil_sb),het_respiration_c13_lay_acc(1:nsoil_sb))), &
              transport_DOC_lay_acc(1:nsoil_sb)
            lun = lun + 1
            ! soilfluxN_lay
            write(lun,'(201F20.10)') time_d, &
              net_mineral_nh4_lay_acc(1:nsoil_sb), plant_uptake_nh4_lay_acc(1:nsoil_sb), microbial_upt_nh4_lay_acc(1:nsoil_sb), &
              myco_uptake_nh4_lay_acc(1:nsoil_sb), asymb_n_fixation_lay_acc(1:nsoil_sb), &
              net_mineral_no3_lay_acc(1:nsoil_sb), plant_uptake_no3_lay_acc(1:nsoil_sb), microbial_upt_no3_lay_acc(1:nsoil_sb), &
              myco_uptake_no3_lay_acc(1:nsoil_sb), transport_DON_lay_acc(1:nsoil_sb)
            lun = lun + 1
            ! soilfluxP_lay
            write(lun,'(301F20.10)') time_d, &
              adsorpt_po4_acc(1:nsoil_sb), weath_po4_acc(1:nsoil_sb), &
              ocl_po4_acc(1:nsoil_sb), desorp_po4_acc(1:nsoil_sb), &
              net_min_po4_lay_acc(1:nsoil_sb), biochem_min_po4_acc(1:nsoil_sb), &
              microbial_upt_po4_acc(1:nsoil_sb), plant_upt_po4_acc(1:nsoil_sb), &
              myco_upt_po4_acc(1:nsoil_sb), trans_solute_po4_acc(1:nsoil_sb), &
              trans_assoc_fast_po4_acc(1:nsoil_sb) , trans_assoc_slow_po4_acc(1:nsoil_sb), &
              trans_occluded_po4_acc(1:nsoil_sb), trans_primary_po4_acc(1:nsoil_sb), &
              transport_DOP_lay_acc(1:nsoil_sb)
            lun = lun + 1

            !! set acc variables to zero
            ! C
            het_respiration_lay_acc(1:nsoil_sb)     = 0.0_wp
            het_respiration_c13_lay_acc(1:nsoil_sb) = 0.0_wp
            het_respiration_c14_lay_acc(1:nsoil_sb) = 0.0_wp
            transport_DOC_lay_acc(1:nsoil_sb)       = 0.0_wp
            ! N
            net_mineral_nh4_lay_acc(1:nsoil_sb)     = 0.0_wp
            plant_uptake_nh4_lay_acc(1:nsoil_sb)    = 0.0_wp
            microbial_upt_nh4_lay_acc(1:nsoil_sb)   = 0.0_wp
            microbial_upt_no3_lay_acc(1:nsoil_sb)   = 0.0_wp
            myco_uptake_nh4_lay_acc(1:nsoil_sb)     = 0.0_wp
            asymb_n_fixation_lay_acc(1:nsoil_sb)    = 0.0_wp
            net_mineral_no3_lay_acc(1:nsoil_sb)     = 0.0_wp
            plant_uptake_no3_lay_acc(1:nsoil_sb)    = 0.0_wp
            myco_uptake_no3_lay_acc(1:nsoil_sb)     = 0.0_wp
            transport_DON_lay_acc(1:nsoil_sb)       = 0.0_wp
            ! P
            adsorpt_po4_acc(1:nsoil_sb)             = 0.0_wp
            weath_po4_acc(1:nsoil_sb)               = 0.0_wp
            ocl_po4_acc(1:nsoil_sb)                 = 0.0_wp
            desorp_po4_acc(1:nsoil_sb)              = 0.0_wp
            net_min_po4_lay_acc(1:nsoil_sb)         = 0.0_wp
            biochem_min_po4_acc(1:nsoil_sb)         = 0.0_wp
            microbial_upt_po4_acc(1:nsoil_sb)       = 0.0_wp
            plant_upt_po4_acc(1:nsoil_sb)           = 0.0_wp 
            myco_upt_po4_acc(1:nsoil_sb)            = 0.0_wp  
            trans_solute_po4_acc(1:nsoil_sb)        = 0.0_wp
            trans_assoc_fast_po4_acc(1:nsoil_sb)    = 0.0_wp
            trans_assoc_slow_po4_acc(1:nsoil_sb)    = 0.0_wp
            trans_occluded_po4_acc(1:nsoil_sb)      = 0.0_wp
            trans_primary_po4_acc(1:nsoil_sb)       = 0.0_wp
            transport_DOP_lay_acc(1:nsoil_sb)       = 0.0_wp
          END IF      ! l_flux_output_interval
        END SELECT    ! simple_1d/jsm
      END IF          ! l_flux_calc_acc_output
    END SELECT        ! land/soil




! Pass back values to process memory
#define Passback_veg_litterfall
#define Passback_veg_exudation
#define Passback_veg_reserve_use
#define Passback_sb_transport
#define Passback_sb_mycorrhiza_export
#define Passback_sb_delta_pools
#include "slm_pools_a_incl.f90"
#undef Passback_veg_litterfall
#undef Passback_veg_exudation
#undef Passback_veg_reserve_use
#undef Passback_sb_transport
#undef Passback_sb_mycorrhiza_export
#undef Passback_sb_delta_pools

  END SUBROUTINE write_output_flux


  !-----------------------------------------------------------------------------------------------------!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Soil forcing                                                                                      !!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !-----------------------------------------------------------------------------------------------------!


  !-----------------------------------------------------------------------------------------------------
  !> Soil forcing file: open files and write column headers
  !! 
  !! open output files, define file names for pools, write column header names
  !! 
  !! lun: start at 100 + 60
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE open_output_files_soil_forcing(tile, iblk, ics, ice, nc, nsoil_sb)

    USE mo_jsb_class,             ONLY: Get_model
    USE mo_base_memory,           ONLY: t_slm_tile
    USE mo_jsb_model_class,       ONLY: t_jsb_model
    USE mo_jsb_process_class,     ONLY: A2L_, L2A_, SEB_, TURB_, SOIL_, HYDRO_, VEG_, HD_, RAD_, ASSIMI_, PHENO_, SB_
    
    ! Use of process configurations (t_PROC_config)
    dsl4jsb_Use_config(ASSIMI_)
    dsl4jsb_Use_config(VEG_)
    dsl4jsb_Use_config(SB_)

    IMPLICIT NONE
    ! 0.1 InOut
    TYPE(t_slm_tile), INTENT(in)          :: tile           !< data structure for dynamic vegetation and soil variables
    INTEGER,          INTENT(in)          :: iblk, &        !< iBlock, start point of chunk, end point of chunk, number of points in the chunk
                                             ics, &
                                             ice, &
                                             nc, &
                                             nsoil_sb 
    ! 0.2 Local
    TYPE(t_jsb_model),POINTER             :: model
    INTEGER        :: lun           ! logical unit number
    CHARACTER(100) :: fname         ! filename
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':open_output_files_soil_forcing'
    ! ---------------------------
    ! 0.3 Declare Memory
    ! Declare process configuration and memory Pointers
    dsl4jsb_Def_config(ASSIMI_)
    dsl4jsb_Def_config(VEG_)
    dsl4jsb_Def_config(SB_)
    ! 0.5 Get Memory
    model  => Get_model(tile%owner_model_id)
    ! Get process config
    dsl4jsb_Get_config(ASSIMI_)
    dsl4jsb_Get_config(VEG_)
    dsl4jsb_Get_config(SB_)

    !! 
    !! output interval will be always timestep !
    !! 

    !> 1.0 define open() statements for file output
    !! 
    lun = 100 + 60
    fname = 'soilphysics'//'_timestep.txt'
    OPEN (UNIT=lun, FILE=fname, ACTION="write", STATUS="unknown")
    lun=lun+1
    fname = 'veglitterfall'//'_timestep.txt'
    OPEN (UNIT=lun, FILE=fname, ACTION="write", STATUS="unknown")
    lun=lun+1

    !> 2.0 write column headers into output files
    !! 
    CALL init_output_files_soil_forcing(nsoil_sb)

  END SUBROUTINE open_output_files_soil_forcing


  !-----------------------------------------------------------------------------------------------------
  !> Soil forcing files: write column headers of output files
  !! 
  !! complete litterfall output (used as input for the stand alone soil model) \n
  !! and some soilphysics variables
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE init_output_files_soil_forcing(nsoil_sb)

    IMPLICIT NONE
    ! 0.1 InOut
    INTEGER, INTENT(in)       :: nsoil_sb
    ! 0.2 Local
    INTEGER                   :: lun, isoil
    CHARACTER(2)              :: layerstr
    CHARACTER(450)            :: p1str,p2str,p3str,p4str,p5str,p6str,p7str,p8str,p9str,p10str

    !> 1.0 write headers
    !!
    !! lun of fluxes starts at a multiple of 100, plus 60
    lun = 100 + 60
    
    ! soilphysics
    p1str=''
    p2str=''
    p3str=''
    p4str=''
    DO isoil = 1,nsoil_sb
      IF(isoil<10) THEN 
        WRITE(layerstr,'(I1)') isoil
      ELSE
        WRITE(layerstr,'(I2)') isoil
      END IF
      WRITE(p1str,'(a)') TRIM(p1str)//' TemperatureSoil_L'//TRIM(layerstr)
      WRITE(p2str,'(a)') TRIM(p2str)//' WaterContentSoil_L'//TRIM(layerstr)
      WRITE(p3str,'(a)') TRIM(p3str)//' PercolationSoil_L'//TRIM(layerstr)
      WRITE(p4str,'(a)') TRIM(p4str)//' RootFractionSoil_L'//TRIM(layerstr)
    END DO
    WRITE(lun,'(9a)') 'Time ', &
                      TRIM(p1str), TRIM(p2str), TRIM(p3str), TRIM(p4str), &
                      ' VegPoolFineRootCarbon ', &
                      'nhx_deposition noy_deposition ', &
                      'nhx_dn15_deposition noy_dn15_deposition ', &
                      'p_deposition '
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


  !-----------------------------------------------------------------------------------------------------
  !> Subroutine to generate the forcing file for the stand alone soil model from 
  !!  litterfall data and SB rate modifiers 
  !! 
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE write_soil_forcing(tile, options, time_d, time_y, l_transient_spinup)

    USE mo_jsb_class,             ONLY: Get_model
    USE mo_base_memory,           ONLY: t_slm_tile, t_jsb_task_options
    USE mo_jsb_model_class,       ONLY: t_jsb_model
    USE mo_jsb_process_class,     ONLY: A2L_, L2A_, SEB_, TURB_, SOIL_, HYDRO_, VEG_, HD_, RAD_, ASSIMI_, PHENO_, SB_
    USE mo_shared_isotope,        ONLY: calc_delta_C13,calc_delta_C14,calc_delta_N15

    ! Use of process memories
    dsl4jsb_Use_memory(A2L_)        !  USE mo_a2l_memory_class,      ONLY: t_a2l_memory
    dsl4jsb_Use_memory(VEG_)        !  USE mo_veg_memory_class,      ONLY: t_veg_memory
    dsl4slm_Use_memory_a_veg
    USE mo_veg_memory_class,        ONLY: veg_pool, veg_litterfall
    dsl4jsb_Use_memory(SOIL_)       !  USE mo_soil_memory_class,     ONLY: t_soil_memory

    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    TYPE(t_slm_tile),           INTENT(inout)     :: tile         !< one tile
    TYPE(t_jsb_task_options),   INTENT(in)        :: options      !< model options
    REAL(wp),                   INTENT(in)        :: time_d       !< time in days
    INTEGER,                    INTENT(in)        :: time_y       !< time in years; the 1st year of a simulation == year 1 (not zero)
    LOGICAL,                    INTENT(in)        :: l_transient_spinup
    ! ---------------------------
    ! 0.2 Local
    TYPE(t_jsb_model),      POINTER       :: model
    TYPE(t_lctlib_element), POINTER       :: lctlib     !< land-cover-type library - parameter across pft's
    REAL(wp)                              :: dC13_leaf,dC13_fine_root,dC13_coarse_root,dC13_sap_wood,dC13_heart_wood,dC13_labile, & 
                                             dC13_reserve,dC13_fruit,dC13_seed_bed
    REAL(wp)                              :: DC14_leaf,DC14_fine_root,DC14_coarse_root,DC14_sap_wood,DC14_heart_wood,DC14_labile, & 
                                             DC14_reserve,DC14_fruit,DC14_seed_bed
    REAL(wp)                              :: dN15_leaf,dN15_fine_root,dN15_coarse_root,dN15_sap_wood,dN15_heart_wood,dN15_labile, & 
                                             dN15_reserve,dN15_fruit,dN15_seed_bed
    INTEGER                               :: lun          ! logical unit number
    INTEGER                               :: iblk, ics, ice, nc, nsoil_sb
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':write_soil_forcing'
    ! ---------------------------
    ! 0.3 Declare Memory
    ! Declare process configuration and memory Pointers
    dsl4jsb_Def_memory(A2L_)                        ! TYPE(t_a2l_memory),     POINTER :: a2l__mem
    dsl4jsb_Def_memory(VEG_)                        ! TYPE(t_veg_memory),     POINTER :: veg__mem
    dsl4jsb_Def_memory(SOIL_)                       ! TYPE(t_soil_memory),    POINTER :: soil__mem
    ! Declare pointers to variables in memory
    ! A2L_
    dsl4jsb_Real2D_onChunk      :: nhx_deposition 
    dsl4jsb_Real2D_onChunk      :: noy_deposition 
    dsl4jsb_Real2D_onChunk      :: nhx_n15_deposition 
    dsl4jsb_Real2D_onChunk      :: noy_n15_deposition 
    dsl4jsb_Real2D_onChunk      :: p_deposition 
    ! VEG_
    dsl4jsb_Real3D_onChunk      :: root_fraction_sl
    ! SOIL_
    dsl4jsb_Real3D_onChunk      :: t_soil_sl
    dsl4jsb_Real3D_onChunk      :: w_soil_sl
    dsl4jsb_Real3D_onChunk      :: percolation_sl
    ! Get local variables from options argument
    iblk      = options%iblk
    ics       = options%ics
    ice       = options%ice
    nc        = options%nc
    nsoil_sb  = options%nsoil_sb
    ! ---------------------------
    ! 0.4 Process Activity, Debug Option
!    IF (.NOT. tile%Is_process_active(xx_)) RETURN
    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')
    ! ---------------------------
    ! 0.5 Get Memory
    model  => Get_model(tile%owner_model_id)
    lctlib => model%lctlib(tile%lcts%lib_id) 
    ! Get process config
!    dsl4jsb_Get_config(xx_)
    ! Get process memories
    dsl4jsb_Get_memory(A2L_)
    dsl4jsb_Get_memory(VEG_)
    dsl4jsb_Get_memory(SOIL_)                         
    ! Get process variables (Set pointers to variables in memory)
    dsl4jsb_Get_var2D_onChunk(A2L_,  nhx_deposition)
    dsl4jsb_Get_var2D_onChunk(A2L_,  noy_deposition)
    dsl4jsb_Get_var2D_onChunk(A2L_,  nhx_n15_deposition)
    dsl4jsb_Get_var2D_onChunk(A2L_,  noy_n15_deposition)
    dsl4jsb_Get_var2D_onChunk(A2L_,  p_deposition)
    ! ---------------------------
    dsl4jsb_Get_var3D_onChunk(VEG_,  root_fraction_sl)
    ! ---------------------------
    dsl4jsb_Get_var3D_onChunk(SOIL_, t_soil_sl)
    dsl4jsb_Get_var3D_onChunk(SOIL_, w_soil_sl)
    dsl4jsb_Get_var3D_onChunk(SOIL_, percolation_sl)
    ! Get process memory pool structures
#define Get_veg_pool
#define Get_veg_litterfall
#include "slm_pools_a_incl.f90"
#undef Get_veg_pool
#undef Get_veg_litterfall



    ! ------------------------------------------------------------------------------------------------------------
    ! Go science
    ! ------------------------------------------------------------------------------------------------------------


    !> 1.0 check output start & end, and write litterfall and soilphysics if appropriate
    !! 
    !! identical output period as pools & fluxes
    IF(.NOT. l_transient_spinup .AND. &
       time_y >= base_config%output_start_first_day_year .AND. time_y <= base_config%output_end_last_day_year) THEN

      lun = 100 + 60
      ! soilphysics
      write(lun,'(87F18.8)') time_d                             , &
                              ! variables with one values per soil layer (max of 20 soil layers)
                              t_soil_sl(ipoi,:)                 , &
                              w_soil_sl(ipoi,:)                 , &
                              percolation_sl(ipoi,:)            , &
                              root_fraction_sl(ipoi,:)          , &
                              ! variables with one value only
                              veg_pool%fine_root%carbon(ipoi)                               , &
                              nhx_deposition(ipoi)                                          , &
                              noy_deposition(ipoi)                                          , &
                              calc_delta_N15(nhx_deposition(ipoi),nhx_n15_deposition(ipoi)) , &
                              calc_delta_N15(noy_deposition(ipoi),noy_n15_deposition(ipoi)) , &
                              p_deposition(ipoi)
      lun=lun+1

      !> convert 13C and 14C in delta notation
      !! 
      dC13_leaf        = calc_delta_C13(veg_litterfall%leaf%carbon(ipoi),veg_litterfall%leaf%carbon13(ipoi))
      dC13_fine_root   = calc_delta_C13(veg_litterfall%fine_root%carbon(ipoi),veg_litterfall%fine_root%carbon13(ipoi))
      dC13_coarse_root = calc_delta_C13(veg_litterfall%coarse_root%carbon(ipoi),veg_litterfall%coarse_root%carbon13(ipoi))
      dC13_sap_wood    = calc_delta_C13(veg_litterfall%sap_wood%carbon(ipoi),veg_litterfall%sap_wood%carbon13(ipoi))
      dC13_heart_wood  = calc_delta_C13(veg_litterfall%heart_wood%carbon(ipoi),veg_litterfall%heart_wood%carbon13(ipoi))
      dC13_labile      = calc_delta_C13(veg_litterfall%labile%carbon(ipoi),veg_litterfall%labile%carbon13(ipoi))
      dC13_reserve     = calc_delta_C13(veg_litterfall%reserve%carbon(ipoi),veg_litterfall%reserve%carbon13(ipoi))
      dC13_fruit       = calc_delta_C13(veg_litterfall%fruit%carbon(ipoi),veg_litterfall%fruit%carbon13(ipoi))
      dC13_seed_bed    = calc_delta_C13(veg_litterfall%seed_bed%carbon(ipoi),veg_litterfall%seed_bed%carbon13(ipoi))

      DC14_leaf        = calc_delta_C14(veg_litterfall%leaf%carbon(ipoi), &
                                        veg_litterfall%leaf%carbon14(ipoi),dC13_leaf)
      DC14_fine_root   = calc_delta_C14(veg_litterfall%fine_root%carbon(ipoi), &
                                        veg_litterfall%fine_root%carbon14(ipoi),dC13_fine_root)
      DC14_coarse_root = calc_delta_C14(veg_litterfall%coarse_root%carbon(ipoi), &
                                        veg_litterfall%coarse_root%carbon14(ipoi),dC13_coarse_root)
      DC14_sap_wood    = calc_delta_C14(veg_litterfall%sap_wood%carbon(ipoi), &
                                        veg_litterfall%sap_wood%carbon14(ipoi),dC13_sap_wood)
      DC14_heart_wood  = calc_delta_C14(veg_litterfall%heart_wood%carbon(ipoi), &
                                        veg_litterfall%heart_wood%carbon14(ipoi),dC13_heart_wood)
      DC14_labile      = calc_delta_C14(veg_litterfall%labile%carbon(ipoi), &
                                        veg_litterfall%labile%carbon14(ipoi),dC13_labile)
      DC14_reserve     = calc_delta_C14(veg_litterfall%reserve%carbon(ipoi), &
                                        veg_litterfall%reserve%carbon14(ipoi),dC13_reserve)
      DC14_fruit       = calc_delta_C14(veg_litterfall%fruit%carbon(ipoi), &
                                        veg_litterfall%fruit%carbon14(ipoi),dC13_fruit)
      DC14_seed_bed    = calc_delta_C14(veg_litterfall%seed_bed%carbon(ipoi), &
                                        veg_litterfall%seed_bed%carbon14(ipoi),dC13_seed_bed)

      dN15_leaf        = calc_delta_N15(veg_litterfall%leaf%nitrogen(ipoi),veg_litterfall%leaf%nitrogen15(ipoi))
      dN15_fine_root   = calc_delta_N15(veg_litterfall%fine_root%nitrogen(ipoi),veg_litterfall%fine_root%nitrogen15(ipoi))
      dN15_coarse_root = calc_delta_N15(veg_litterfall%coarse_root%nitrogen(ipoi),veg_litterfall%coarse_root%nitrogen15(ipoi))
      dN15_sap_wood    = calc_delta_N15(veg_litterfall%sap_wood%nitrogen(ipoi),veg_litterfall%sap_wood%nitrogen15(ipoi))
      dN15_heart_wood  = calc_delta_N15(veg_litterfall%heart_wood%nitrogen(ipoi),veg_litterfall%heart_wood%nitrogen15(ipoi))
      dN15_labile      = calc_delta_N15(veg_litterfall%labile%nitrogen(ipoi),veg_litterfall%labile%nitrogen15(ipoi))
      dN15_reserve     = calc_delta_N15(veg_litterfall%reserve%nitrogen(ipoi),veg_litterfall%reserve%nitrogen15(ipoi))
      dN15_fruit       = calc_delta_N15(veg_litterfall%fruit%nitrogen(ipoi),veg_litterfall%fruit%nitrogen15(ipoi))
      dN15_seed_bed    = calc_delta_N15(veg_litterfall%seed_bed%nitrogen(ipoi),veg_litterfall%seed_bed%nitrogen15(ipoi))

      ! veglitterfall
      write(lun,'(60F18.8)') time_d, &   ! 6 variables for (currently) 9 elements + time_d = 55 columns
        ! leaf
        veg_litterfall%leaf%carbon(ipoi)*1000._wp, &
        veg_litterfall%leaf%nitrogen(ipoi)*1000000._wp, &
        veg_litterfall%leaf%phosphorus(ipoi)*1000000._wp    , &
        dC13_leaf, DC14_leaf, dN15_leaf, & 
        ! fine_root
        veg_litterfall%fine_root%carbon(ipoi)*1000._wp, &
        veg_litterfall%fine_root%nitrogen(ipoi)*1000000._wp, &
        veg_litterfall%fine_root%phosphorus(ipoi)*1000000._wp   , &
        dC13_fine_root, DC14_fine_root, dN15_fine_root, & 
        ! coarse_root
        veg_litterfall%coarse_root%carbon(ipoi)*1000._wp, &
        veg_litterfall%coarse_root%nitrogen(ipoi)*1000000._wp, &
        veg_litterfall%coarse_root%phosphorus(ipoi)*1000000._wp   , &
        dC13_coarse_root, DC14_coarse_root, dN15_coarse_root, & 
        ! sap_wood
        veg_litterfall%sap_wood%carbon(ipoi)*1000._wp, &
        veg_litterfall%sap_wood%nitrogen(ipoi)*1000000._wp, &
        veg_litterfall%sap_wood%phosphorus(ipoi)*1000000._wp    , &
        dC13_sap_wood, DC14_sap_wood, dN15_sap_wood, & 
        ! heart_wood
        veg_litterfall%heart_wood%carbon(ipoi)*1000._wp, &
        veg_litterfall%heart_wood%nitrogen(ipoi)*1000000._wp, &
        veg_litterfall%heart_wood%phosphorus(ipoi)*1000000._wp    , &
        dC13_heart_wood, DC14_heart_wood, dN15_heart_wood, & 
        ! labile
        veg_litterfall%labile%carbon(ipoi)*1000._wp, &
        veg_litterfall%labile%nitrogen(ipoi)*1000000._wp, &
        veg_litterfall%labile%phosphorus(ipoi)*1000000._wp    , &
        dC13_labile, DC14_labile, dN15_labile, & 
        ! reserve
        veg_litterfall%reserve%carbon(ipoi)*1000._wp, &
        veg_litterfall%reserve%nitrogen(ipoi)*1000000._wp, &
        veg_litterfall%reserve%phosphorus(ipoi)*1000000._wp   , &
        dC13_reserve, DC14_reserve, dN15_reserve, & 
        ! fruit
        veg_litterfall%fruit%carbon(ipoi)*1000._wp, &
        veg_litterfall%fruit%nitrogen(ipoi)*1000000._wp, &
        veg_litterfall%fruit%phosphorus(ipoi)*1000000._wp   , &
        dC13_fruit, DC14_fruit, dN15_fruit, & 
        ! seed_bed
        veg_litterfall%seed_bed%carbon(ipoi)*1000._wp, &
        veg_litterfall%seed_bed%nitrogen(ipoi)*1000000._wp, &
        veg_litterfall%seed_bed%phosphorus(ipoi)*1000000._wp  , &
        dC13_seed_bed, DC14_seed_bed, dN15_seed_bed 
      lun=lun+1

    END IF 

! Pass back values to process memory
#define Passback_veg_pool
#define Passback_veg_litterfall
#include "slm_pools_a_incl.f90"
#undef Passback_veg_pool
#undef Passback_veg_litterfall

  END SUBROUTINE write_soil_forcing 



  !-----------------------------------------------------------------------------------------------------!
  !!                                                                                                   !!
  !! Close output files                                                                                !!
  !!                                                                                                   !!
  !-----------------------------------------------------------------------------------------------------!


  !-----------------------------------------------------------------------------------------------------
  !> Pools and Fluxes: close files
  !! 
  !! at this point of the simulation all data has been written
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE close_output_files

    IMPLICIT NONE
    ! 0.1 InOut
    ! 
    ! 0.2 Local
    INTEGER :: lun, lun_start, lun_end
    LOGICAL :: opened

    ! 1.0 pools
    lun_start = 100
    lun_end = lun_start + 29
    DO lun = lun_start,lun_end
       INQUIRE(lun,opened=opened)
       IF(opened) CLOSE(lun) 
    ENDDO

    ! 2.0 fluxes
    lun_start = 100 + 30
    lun_end = lun_start + 29
    DO lun = lun_start,lun_end
       INQUIRE(lun,opened=opened)
       IF(opened) CLOSE(lun) 
    ENDDO

  END SUBROUTINE close_output_files


  !-----------------------------------------------------------------------------------------------------
  !> Soil forcing files: close files
  !! 
  !! at this point of the simulation all data has been written
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE close_output_files_soil_forcing

    IMPLICIT NONE
    ! 0.1 InOut
    ! 
    ! 0.2 Local
    INTEGER :: lun, lun_start, lun_end
    LOGICAL :: opened

    ! 1.0 close soil_forcing output files
    lun_start = 100 + 60
    lun_end = lun_start + 29
    DO lun = lun_start,lun_end
       INQUIRE(lun,opened=opened)
       IF(opened) CLOSE(lun) 
    ENDDO

  END SUBROUTINE close_output_files_soil_forcing

END MODULE mo_shared_output
