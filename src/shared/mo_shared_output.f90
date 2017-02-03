!>
!! @brief functions / subroutines that help to create output from plant_model_driver
!!
!! @author
!!  Jan Engel, MPI-BGC Jena
!!
!! @par Revision History
!! First version by Jan Engel (2016-03-14)
!!
MODULE mo_shared_output
  !! this module may be not necessary anymore when QUINCY is linked with JSBACH !!

#include "dsl_jsbach.inc"

  USE mo_kind,                  ONLY: wp, eps4, eps8
  USE mo_base_constants,        ONLY: ntile
  USE mo_jsb_lctlib_class,      ONLY: t_lctlib_element
  USE mo_jsb_model_init,        ONLY: fluxnet_output_status   ! used for 'special case: fluxnet output'
  USE mo_base_config_class,     ONLY: base_config
  USE mo_veg_config_class,      ONLY: veg_config
  USE mo_soil_config_class,     ONLY: soil_config
  USE mo_assimi_config_class,   ONLY: assimi_config


  IMPLICIT NONE

  REAL(wp), SAVE                   :: output_interval_pool   ! equals 'dtime/one_day' or is a mulitple of it
  REAL(wp), SAVE                   :: output_interval_flux   ! equals 'dtime/one_day' or is a mulitple of it
  ! flux accumulation, see write_output_flux
  ! these are set to zero frequently in section 1.2 of write_output_flux()
  REAL(wp), DIMENSION(ntile), SAVE :: gpp_acc             
  REAL(wp), DIMENSION(ntile), SAVE :: growth_resp_acc      
  REAL(wp), DIMENSION(ntile), SAVE :: maint_resp_acc       
  REAL(wp), DIMENSION(ntile), SAVE :: n_transform_resp_acc
  REAL(wp), DIMENSION(ntile), SAVE :: n_fixation_resp_acc
  REAL(wp), DIMENSION(ntile), SAVE :: npp_acc              
  REAL(wp), DIMENSION(ntile), SAVE :: res_use_c_acc        
  REAL(wp), DIMENSION(ntile), SAVE :: gpp_C13_acc             
  REAL(wp), DIMENSION(ntile), SAVE :: growth_resp_C13_acc      
  REAL(wp), DIMENSION(ntile), SAVE :: maint_resp_C13_acc       
  REAL(wp), DIMENSION(ntile), SAVE :: gpp_C14_acc             
  REAL(wp), DIMENSION(ntile), SAVE :: growth_resp_C14_acc      
  REAL(wp), DIMENSION(ntile), SAVE :: maint_resp_C14_acc       
  REAL(wp), DIMENSION(ntile), SAVE :: uptake_nh4_acc         
  REAL(wp), DIMENSION(ntile), SAVE :: uptake_no3_acc         
  REAL(wp), DIMENSION(ntile), SAVE :: n_fixation_acc         
  REAL(wp), DIMENSION(ntile), SAVE :: n_loss_acc         
  REAL(wp), DIMENSION(ntile), SAVE :: recycling_n_acc     
  REAL(wp), DIMENSION(ntile), SAVE :: res_use_n_acc        
  REAL(wp), DIMENSION(ntile), SAVE :: uptake_p_acc         
  REAL(wp), DIMENSION(ntile), SAVE :: recycling_p_acc     
  REAL(wp), DIMENSION(ntile), SAVE :: res_use_p_acc        
  REAL(wp), DIMENSION(ntile), SAVE :: transpiration_acc    
  REAL(wp), DIMENSION(ntile), SAVE :: drainage_acc    
  REAL(wp), DIMENSION(ntile), SAVE :: gc_acc               
  REAL(wp), DIMENSION(ntile), SAVE :: ga_acc               
  REAL(wp), DIMENSION(ntile), SAVE :: ppfd_acc 
  REAL(wp), DIMENSION(ntile), SAVE :: rnet_acc 
  REAL(wp), DIMENSION(ntile), SAVE :: qh_acc 
  REAL(wp), DIMENSION(ntile), SAVE :: qle_acc 
  REAL(wp), DIMENSION(ntile), SAVE :: qg_acc 
  REAL(wp), DIMENSION(ntile), SAVE :: albedo_acc 
  REAL(wp), DIMENSION(ntile), SAVE :: appfd_acc 
  INTEGER, DIMENSION(ntile), SAVE  :: ndat_acc ! number of accumulated data points

  CHARACTER(len=*), PARAMETER, PRIVATE :: modname = 'mo_shared_output'

CONTAINS

  !-----------------------------------------------------------------------------------------------------
  !> initialize variables for interval length and flux accumulation
  !! 
  !! 
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE init_veg_output_variables
    IMPLICIT NONE
    ! interval variables
    output_interval_pool = 0._wp  ! equals 'dtime/one_day' or is a mulitple of it
    output_interval_flux = 0._wp  ! equals 'dtime/one_day' or is a mulitple of it
    ! flux accumulation per tile [DIMENSION(ntile)], see write_output_flux
    ! these are set to zero frequently in section 1.2 of write_output_flux()
    gpp_acc           = 0._wp
    growth_resp_acc   = 0._wp
    maint_resp_acc    = 0._wp   
    n_transform_resp_acc    = 0._wp   
    n_fixation_resp_acc    = 0._wp   
    npp_acc           = 0._wp   
    res_use_c_acc     = 0._wp   
    gpp_C13_acc       = 0.0_wp
    growth_resp_C13_acc = 0.0_wp
    maint_resp_C13_acc = 0.0_wp
    gpp_C14_acc       = 0.0_wp
    growth_resp_C14_acc = 0.0_wp
    maint_resp_C14_acc = 0.0_wp
    uptake_nh4_acc      = 0._wp   
    uptake_no3_acc      = 0._wp   
    recycling_n_acc  = 0._wp  
    res_use_n_acc     = 0._wp  
    n_fixation_acc    = 0._wp   
    n_loss_acc        = 0._wp   
    uptake_p_acc      = 0._wp
    recycling_p_acc  = 0._wp  
    res_use_p_acc     = 0._wp   
    transpiration_acc = 0._wp 
    drainage_acc      = 0._wp 
    gc_acc            = 0._wp
    ga_acc            = 0._wp
    ppfd_acc          = 0._wp
    rnet_acc          = 0._wp
    qh_acc            = 0._wp
    qle_acc           = 0._wp
    qg_acc           = 0._wp
    albedo_acc        = 0._wp
    appfd_acc         = 0._wp
    ndat_acc          = 0
  END SUBROUTINE init_veg_output_variables


  !-----------------------------------------------------------------------------------------------------
  !> the input interval string is transformed into a number, which is used for timely model output
  !! 
  !! applied for pools and flux output
  !-----------------------------------------------------------------------------------------------------
  FUNCTION calc_output_time_interval(interval_string) RESULT(interval_number)

    USE mo_base_constants,     ONLY: one_day, dtime

    IMPLICIT NONE
    ! 0.1 InOut
    CHARACTER(8), INTENT(in) :: interval_string      !< string input
    REAL(wp)                 :: interval_number      !< number output

    !> 1.0 create interval_number depending on pool/flux
    !!
    SELECT CASE( TRIM(interval_string) )
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

  !-----------------------------------------------------------------------------------------------------
  !> Pools: write column headers of output files
  !! 
  !! 
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE init_output_writing_pool(itile)

    IMPLICIT NONE
    ! 0.1 InOut
    INTEGER, INTENT(in) :: itile    !< number of current tile
    ! 0.2 Local
    INTEGER :: lun,isoil
    CHARACTER(2) :: layerstr
    CHARACTER(100) :: p1str,p2str,p3str,p4str,p5str,p6str,p7str,p8str

    !> 1.0 write headers
    !!
    !! lun of pools starts at a multiple of 100
    lun = 100 * itile
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
                          'Fruit_N Labile_N Reserve_N SeedBank_N Total_N Total_N15'
         lun=lun+1
         WRITE(lun,'(2a)') 'Time Leaf_P FineRoot_P SapWood_P HeartWood_P CoarseRoot_P ', & 
                          'Fruit_P Labile_P Reserve_P SeedBank_P Total_P'
         lun=lun+1
         WRITE(lun,'(a)') 'Time RootZoneSoilMoisture RootZoneSoilWaterPotential BetaSoilGS BetaSoilPS'
         lun=lun+1
         WRITE(lun,'(a)') 'Time LitterN LitterP MineralNH4 MineralNO3 MineralP'
         lun=lun+1
         WRITE(lun,'(2a)') 'Time CanopyLayer f_et f_rub f_chl f_pepc f_oth leaf_nitrogen'
         lun=lun+1
         WRITE(lun,'(2a)') &
         & 'Time Height Diameter LAI ind_density La2Sa L2R S2C BCF CN_target dNPP dNup cn_plant allometry_npp ', & 
         & 'growing_season cost_nup cost_bnf unit_transp transp'

      CASE('soil')
         p1str=''
         p2str=''
         p3str=''
         p4str=''
         p5str=''
         p6str=''
         p7str=''
         p8str=''
         DO isoil = 1,soil_config%nsoil_water
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
            WRITE(p6str,'(a)') TRIM(p6str)//' Mic_L'//TRIM(layerstr)
            WRITE(p7str,'(a)') TRIM(p7str)//' Res_L'//TRIM(layerstr)
            WRITE(p8str,'(a)') TRIM(p8str)//' aRes_L'//TRIM(layerstr)
         ENDDO
         WRITE(lun,'(9a)') 'Time ',TRIM(p1str),TRIM(p2str),TRIM(p3str),TRIM(p4str),TRIM(p5str),TRIM(p6str),TRIM(p7str),TRIM(p8str)
         lun=lun+1                 
         WRITE(lun,'(9a)') 'Time ',TRIM(p1str),TRIM(p2str),TRIM(p3str),TRIM(p4str),TRIM(p5str),TRIM(p6str),TRIM(p7str),TRIM(p8str)
         lun=lun+1                 
         WRITE(lun,'(9a)') 'Time ',TRIM(p1str),TRIM(p2str),TRIM(p3str),TRIM(p4str),TRIM(p5str),TRIM(p6str),TRIM(p7str),TRIM(p8str)
         lun=lun+1                 
         WRITE(lun,'(9a)') 'Time ',TRIM(p1str),TRIM(p2str),TRIM(p3str),TRIM(p4str),TRIM(p5str),TRIM(p6str),TRIM(p7str),TRIM(p8str)
         lun=lun+1                 
         WRITE(lun,'(9a)') 'Time ',TRIM(p1str),TRIM(p2str),TRIM(p3str),TRIM(p4str),TRIM(p5str),TRIM(p6str),TRIM(p7str),TRIM(p8str)
         lun=lun+1                 
    END SELECT
     
  END SUBROUTINE init_output_writing_pool

  !-----------------------------------------------------------------------------------------------------
  !> Fluxes: write column headers of output files
  !! 
  !! 
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE init_output_writing_flux(itile)

    IMPLICIT NONE
    ! 0.1 InOut
    INTEGER, INTENT(in) :: itile    !< number of current tile
    ! 0.2 Local
    INTEGER :: lun

    !> 1.0 write headers
    !!
    !! lun of fluxes starts at a multiple of 100, plus 30
    lun = (100 * itile) + 30
    SELECT CASE(TRIM(base_config%slm_model_name))
      CASE('land','plant')
        WRITE(lun,'(a)') 'Time GPP GrowthResp MaintResp NTransFResp NFixResp NPP ReserveUse'
        lun=lun+1   
        WRITE(lun,'(2a)') 'Time GPP_dC13 GrowthResp_dC13 MaintResp_dC13 ', &
                          'GPP_DC14 GrowthResp_DC14 MaintResp_DC14'
        lun=lun+1
        WRITE(lun,'(a)') 'Time NUP-NH4 NUP-NO3 Resorption_N ReserveUse NFixation Nloss'
        lun=lun+1
        WRITE(lun,'(a)') 'Time PUP Resorption_P ReserveUse'
        lun=lun+1
        WRITE(lun,'(a)') 'Time Transpiration Drainage GC GA'
        lun=lun+1
        WRITE(lun,'(a)') 'Time Rnet Qh Qle Qg Albedo PPFD aPPFD fAPPFD' 

      CASE('soil')
    END SELECT
  END SUBROUTINE init_output_writing_flux

  !-----------------------------------------------------------------------------------------------------
  !> Soil forcing files: write column headers of output files
  !! 
  !! complete litterfall output (used as input for the stand alone soil model) \n
  !! and some soilphysics variables
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE init_output_files_soil_forcing(itile)

    IMPLICIT NONE
    ! 0.1 InOut
    INTEGER, INTENT(in) :: itile    !< number of current tile
    ! 0.2 Local
    INTEGER :: lun

    !> 1.0 write headers
    !!
    !! lun of fluxes starts at a multiple of 100, plus 60
    lun = (100 * itile) + 60
    WRITE(lun,'(4a)') 'Time stand alone RTM_decomposition*soil_layer RTM_mic_uptake*soil_layer RTM_sorption*soil_layer', &
                      'RTM_desorption*soil_layer RTM_hsc*soil_layer', & 
                      'RMM_decomposition*soil_layer RMM_mic_uptake*soil_layer RMM_sorption*soil_layer', &
                      'RMM_desorption*soil_layer RMM_hsc*soil_layer' 

    lun=lun+1   
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
  END SUBROUTINE init_output_files_soil_forcing


  !-----------------------------------------------------------------------------------------------------
  !> Pools: open files and write column headers
  !! 
  !! open output files, define file names for pools, write column header names
  !! 
  !! lun: start at 100 * itile
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE open_pool_output_files(itile)

    IMPLICIT NONE
    ! 0.1 InOut
    INTEGER, INTENT(in) :: itile    !< number of current tile
    ! 0.2 Local
    INTEGER        :: lun           ! logical unit number
    CHARACTER(100) :: fname         ! filename
    CHARACTER(30)  :: opt           ! option according to namlist parameters
    CHARACTER(5)   :: char_itile    ! itile transferred into a character -> 6 digits == max number of itile is 99999


    !> 1.0 init
    opt = " "

    !> 1.1 transform interval string to interval number (e.g. week -> 7)
    !! 
    !! it is not necessary to pass 'base_config%output_interval_pool' as argument to open_pool_output_files
    !! as 'base_config' is global & defined with SAVE and initialized before this function is called
    output_interval_pool = calc_output_time_interval(base_config%output_interval_pool)

    !> 1.2 transform itile to char_itile
    !! 
    IF(itile < 10) THEN
      WRITE(char_itile,'(I1)') itile   ! converting integer to string using a 'internal file', I1 == input is int with 2 digits
    ELSE IF(itile > 9 .AND. itile < 100) THEN
      WRITE(char_itile,'(I2)') itile   ! converting integer to string using a 'internal file', I2 == input is int with 2 digits
    ELSE IF(itile > 99 .AND. itile < 1000) THEN
      WRITE(char_itile,'(I3)') itile   ! converting integer to string using a 'internal file', I3 == input is int with 2 digits
    ELSE
      WRITE(char_itile,'(I4)') itile   ! converting integer to string using a 'internal file', I4 == input is int with 2 digits
    END IF

    !> 2.1 get options (opt) for filenames
    !! 
    !! 
    IF (assimi_config%flag_optimal_Nfraction) opt = trim(opt) //  "_optNfrac"
    IF (veg_config%flag_optimal_Nleaf) opt = trim(opt) // "_optNleaf"
    IF (veg_config%flag_optimal_allocation) opt = trim(opt) // "_optAlloc"
     

    !> 2.2 define open() statements for file output
    !! 
    !! 
    lun = 100 * itile
    SELECT CASE(TRIM(base_config%slm_model_name))
      CASE('land','plant')
         fname = 'vegpoolC_tile'//TRIM(char_itile)//'_'//TRIM(base_config%output_interval_pool)//TRIM(opt)//'.txt'
         OPEN (UNIT=lun, FILE=fname, ACTION="write", STATUS="unknown")
         lun=lun+1
         fname = 'vegpooldC13dC14_tile'//TRIM(char_itile)//'_'//TRIM(base_config%output_interval_pool)//TRIM(opt)//'.txt'
         OPEN (UNIT=lun, FILE=fname, ACTION="write", STATUS="unknown")
         lun=lun+1
         fname = 'vegpoolN_tile'//TRIM(char_itile)//'_'//TRIM(base_config%output_interval_pool)//TRIM(opt)//'.txt'
         OPEN (UNIT=lun, FILE=fname, ACTION="write", STATUS="unknown")
         lun=lun+1
         fname = 'vegpoolP_tile'//TRIM(char_itile)//'_'//TRIM(base_config%output_interval_pool)//TRIM(opt)//'.txt'
         OPEN (UNIT=lun, FILE=fname, ACTION="write", STATUS="unknown")
         lun=lun+1
         fname = 'vegpoolH2O_tile'//TRIM(char_itile)//'_'//TRIM(base_config%output_interval_pool)//TRIM(opt)//'.txt'
         OPEN (UNIT=lun, FILE=fname, ACTION="write", STATUS="unknown")
         lun=lun+1
         fname = 'soilpoolNP_tile'//TRIM(char_itile)//'_'//TRIM(base_config%output_interval_pool)//TRIM(opt)//'.txt'
         OPEN (UNIT=lun, FILE=fname, ACTION="write", STATUS="unknown")
         lun=lun+1
         fname = 'veg_canopy_tile'//TRIM(char_itile)//'_'//TRIM(base_config%output_interval_pool)//TRIM(opt)//'.txt'
         OPEN (UNIT=lun, FILE=fname, ACTION="write", STATUS="unknown")
         lun=lun+1
         fname = 'veg_diagnostics_tile'//TRIM(char_itile)//'_'//TRIM(base_config%output_interval_pool)//TRIM(opt)//'.txt'
         OPEN (UNIT=lun, FILE=fname, ACTION="write", STATUS="unknown")

      CASE('soil')
         fname = 'soilpoolCorg_tile'//TRIM(char_itile)//'_'//TRIM(base_config%output_interval_pool)//TRIM(opt)//'.txt'
         OPEN (UNIT=lun, FILE=fname, ACTION="write", STATUS="unknown")
         lun=lun+1
         fname = 'soilpooldC13org_tile'//TRIM(char_itile)//'_'//TRIM(base_config%output_interval_pool)//TRIM(opt)//'.txt'
         OPEN (UNIT=lun, FILE=fname, ACTION="write", STATUS="unknown")
         lun=lun+1
         fname = 'soilpoolDC14org_tile'//TRIM(char_itile)//'_'//TRIM(base_config%output_interval_pool)//TRIM(opt)//'.txt'
         OPEN (UNIT=lun, FILE=fname, ACTION="write", STATUS="unknown")
         lun=lun+1
         fname = 'soilpoolNorg_tile'//TRIM(char_itile)//'_'//TRIM(base_config%output_interval_pool)//TRIM(opt)//'.txt'
         OPEN (UNIT=lun, FILE=fname, ACTION="write", STATUS="unknown")
         lun=lun+1
         fname = 'soilpoolPorg_tile'//TRIM(char_itile)//'_'//TRIM(base_config%output_interval_pool)//TRIM(opt)//'.txt'
         OPEN (UNIT=lun, FILE=fname, ACTION="write", STATUS="unknown")
         lun=lun+1
    END SELECT

    !> 3.0 write column headers into output files
    !! 
    !! 
    CALL init_output_writing_pool(itile)

  END SUBROUTINE open_pool_output_files

  !-----------------------------------------------------------------------------------------------------
  !> Fluxes: open files and write column headers
  !! 
  !! open output files, define file names for pools, write column header names
  !! 
  !! lun: start at 100 * itile + 30
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE open_flux_output_files(itile)

    IMPLICIT NONE
    ! 0.1 InOut
    INTEGER, INTENT(in) :: itile    !< number of current tile
    ! 0.2 Local
    INTEGER        :: lun           ! logical unit number
    CHARACTER(100) :: fname         ! filename
    CHARACTER(30)  :: opt           ! option according to namlist parameters
    CHARACTER(5)   :: char_itile    ! itile transferred into a character -> 6 digits == max number of itile is 99999


    !> 1.0 init
    opt = " "

    !> 1.1 transform interval string to interval number (e.g. week -> 7)
    !! 
    !! it is not necessary to pass 'base_config%output_interval_pool' as argument to open_pool_output_files
    !! as 'base_config' is global & defined with SAVE and initialized before this function is called
    output_interval_flux = calc_output_time_interval(base_config%output_interval_flux)

    !> 1.2 transform itile to char_itile
    !! 
    IF(itile < 10) THEN
      WRITE(char_itile,'(I1)') itile   ! converting integer to string using a 'internal file', I1 == input is int with 2 digits
    ELSE IF(itile > 9 .AND. itile < 100) THEN
      WRITE(char_itile,'(I2)') itile   ! converting integer to string using a 'internal file', I2 == input is int with 2 digits
    ELSE IF(itile > 99 .AND. itile < 1000) THEN
      WRITE(char_itile,'(I3)') itile   ! converting integer to string using a 'internal file', I3 == input is int with 2 digits
    ELSE
      WRITE(char_itile,'(I4)') itile   ! converting integer to string using a 'internal file', I4 == input is int with 2 digits
    END IF

    !> 2.1 get options (opt) for filenames
    !! 
    !! 
    IF (assimi_config%flag_optimal_Nfraction) opt = trim(opt) //  "_optNfrac"
    IF (veg_config%flag_optimal_Nleaf) opt = trim(opt) // "_optNleaf"
    IF (veg_config%flag_optimal_allocation) opt = trim(opt) // "_optAlloc"

    !> 2.2 define open() statements for file output
    !! 
    !! 
    lun = (100 * itile) + 30
    SELECT CASE(TRIM(base_config%slm_model_name))
      CASE('land','plant')
        fname = 'vegfluxC_tile'//TRIM(char_itile)//'_'//TRIM(base_config%output_interval_flux)//TRIM(opt)//'.txt'
        OPEN (UNIT=lun, FILE=fname, ACTION="write", STATUS="unknown")
        lun=lun+1
        fname = 'vegfluxdC13dC14_tile'//TRIM(char_itile)//'_'//TRIM(base_config%output_interval_flux)//TRIM(opt)//'.txt'
        OPEN (UNIT=lun, FILE=fname, ACTION="write", STATUS="unknown")
        lun=lun+1
        fname = 'vegfluxN_tile'//TRIM(char_itile)//'_'//TRIM(base_config%output_interval_flux)//TRIM(opt)//'.txt'
        OPEN (UNIT=lun, FILE=fname, ACTION="write", STATUS="unknown")
        lun=lun+1
        fname = 'vegfluxP_tile'//TRIM(char_itile)//'_'//TRIM(base_config%output_interval_flux)//TRIM(opt)//'.txt'
        OPEN (UNIT=lun, FILE=fname, ACTION="write", STATUS="unknown")
        lun=lun+1
        fname = 'vegfluxH2O_tile'//TRIM(char_itile)//'_'//TRIM(base_config%output_interval_flux)//TRIM(opt)//'.txt'
        OPEN (UNIT=lun, FILE=fname, ACTION="write", STATUS="unknown")
        lun=lun+1
        fname = 'vegfluxEnergy_tile'//TRIM(char_itile)//'_'//TRIM(base_config%output_interval_flux)//TRIM(opt)//'.txt'
        OPEN (UNIT=lun, FILE=fname, ACTION="write", STATUS="unknown")

      CASE('soil')
    END SELECT

    !> 3.0 write column headers into output files
    !! 
    !! 
    CALL init_output_writing_flux(itile)

  END SUBROUTINE open_flux_output_files

  !-----------------------------------------------------------------------------------------------------
  !> Soil forcing file: open files and write column headers
  !! 
  !! open output files, define file names for pools, write column header names
  !! 
  !! lun: start at 100 * itile + 60
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE open_output_files_soil_forcing(itile)

    IMPLICIT NONE
    ! 0.1 InOut
    INTEGER, INTENT(in) :: itile    !< number of current tile
    ! 0.2 Local
    INTEGER        :: lun           ! logical unit number
    CHARACTER(100) :: fname         ! filename
    CHARACTER(30)  :: opt           ! option according to namlist parameters
    CHARACTER(5)   :: char_itile    ! itile transferred into a character -> 6 digits == max number of itile is 99999


    !> 1.0 init
    opt = " "

!    !> 1.1 transform interval string to interval number 
!    !! 
!    !! output interval will be always timestep
!    output_interval_litterfall_soilphysics = calc_output_time_interval('timestep')

    !> 1.1 transform itile to char_itile
    !! 
    IF(itile < 10) THEN
      WRITE(char_itile,'(I1)') itile   ! converting integer to string using a 'internal file', I1 == input is int with 2 digits
    ELSE IF(itile > 9 .AND. itile < 100) THEN
      WRITE(char_itile,'(I2)') itile   ! converting integer to string using a 'internal file', I2 == input is int with 2 digits
    ELSE IF(itile > 99 .AND. itile < 1000) THEN
      WRITE(char_itile,'(I3)') itile   ! converting integer to string using a 'internal file', I3 == input is int with 2 digits
    ELSE
      WRITE(char_itile,'(I4)') itile   ! converting integer to string using a 'internal file', I4 == input is int with 2 digits
    END IF

    !> 2.1 get options (opt) for filenames
    !! 
    !! 
    IF (assimi_config%flag_optimal_Nfraction) opt = trim(opt) //  "_optNfrac"
    IF (veg_config%flag_optimal_Nleaf) opt = trim(opt) // "_optNleaf"
    IF (veg_config%flag_optimal_allocation) opt = trim(opt) // "_optAlloc"

    !> 2.2 define open() statements for file output
    !! 
    !! 
    lun = (100 * itile) + 60
    fname = 'soilphysics_tile'//TRIM(char_itile)//'_timestep'//TRIM(opt)//'.txt'
    OPEN (UNIT=lun, FILE=fname, ACTION="write", STATUS="unknown")
    lun=lun+1
    fname = 'veglitterfall_tile'//TRIM(char_itile)//'_timestep'//TRIM(opt)//'.txt'
    OPEN (UNIT=lun, FILE=fname, ACTION="write", STATUS="unknown")

    !> 3.0 write column headers into output files
    !! 
    !! 
    CALL init_output_files_soil_forcing(itile)

  END SUBROUTINE open_output_files_soil_forcing

  !-----------------------------------------------------------------------------------------------------
  !> Pools: write data to files in a certain interval
  !! 
  !! write-interval accordingly to namelist
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE write_output_pool(tile, options, time_d, time_y, itile)

    USE mo_jsb_class,             ONLY: Get_model
    USE mo_base_memory,           ONLY: t_slm_tile, t_model, t_jsb_task_options
    USE mo_veg_constants,         ONLY: itree
    USE mo_assimi_config_class,   ONLY: assimi_config
    USE mo_shared_isotope,        ONLY: calc_delta_C13,calc_delta_C14
    USE mo_jsb_process_class,     ONLY: ATM2LAND_, LAND2ATM_, SEB_, TURB_, SOIL_, HYDRO_, VEG_, HD_, RAD_, ASSIMI_, PHENO_, SB_

    ! Use of process memories
    dsl4jsb_Use_memory(ASSIMI_)     !  USE mo_assimi_memory_class,   ONLY: t_assimi_memory
    dsl4jsb_Use_memory(VEG_)        !  USE mo_veg_memory_class,      ONLY: t_veg_memory
    dsl4jsb_Use_memory(PHENO_)      !  USE mo_pheno_memory_class,    ONLY: t_pheno_memory
    dsl4jsb_Use_memory(SOIL_)       !  USE mo_soil_memory_class,     ONLY: t_soil_memory
    dsl4jsb_Use_memory(SB_)         !  USE mo_sb_memory_class,       ONLY: t_sb_memory

    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    TYPE(t_slm_tile),         INTENT(in), TARGET  :: tile   !< one tile
    TYPE(t_jsb_task_options), INTENT(in) :: options         !< model options
    REAL(wp), INTENT(in)                 :: time_d      !< time in days
    INTEGER,  INTENT(in)                 :: time_y      !< time in years; the 1st year of a simulation == year 1 (not zero)
    INTEGER,  INTENT(in)                 :: itile       !< number of current tile
    ! ---------------------------
    ! 0.2 Local
    TYPE(t_model),          POINTER      :: model
    TYPE(t_lctlib_element), POINTER      :: lctlib     !< land-cover-type library - parameter across pft's
    INTEGER                          :: lun            ! locical unit number
    INTEGER                          :: icanopy         ! loop over ncanopy
    REAL(wp)                         :: hlp1           ! helper variable
    REAL(wp)                         :: l2r, s2c, bcf  !helper ratio variables
    REAl(wp)                         :: growing_season_float !output variable for growing season
    REAL(wp)                         :: dC13leaf,dC13fine_root,dC13sap_wood,dC13heart_wood,dC13coarse_root, & ! delta C13 in pools 
                                        dC13fruit,dC13labile,dC13reserve,dC13seed_bed
    INTEGER                          :: iblk, ics, ice, nc
    ! ---------------------------
    ! 0.3 Declare Memory
    ! Declare process configuration and memory Pointers
    dsl4jsb_Def_memory(ASSIMI_)                           ! TYPE(t_assimi_memory),   POINTER :: assimi__mem
    dsl4jsb_Def_memory(VEG_)                              ! TYPE(t_veg_memory),      POINTER :: veg__mem
    dsl4jsb_Def_memory(PHENO_)                            ! TYPE(t_pheno_memory),    POINTER :: pheno__mem
    dsl4jsb_Def_memory(SOIL_)                             ! TYPE(t_soil_memory),     POINTER :: soil__mem
    dsl4jsb_Def_memory(SB_)                               ! TYPE(t_sb_memory),       POINTER :: sb__mem
    ! Declare pointers to variables in memory
    dsl4jsb_Logical2D_onChunk                 :: growing_season
!! *slmtemp
    TYPE(t_assimi_memory),        POINTER     :: assimi   ! local pointer to veg__mem for SLM usability
    TYPE(t_veg_memory),           POINTER     :: veg      ! local pointer to veg__mem for SLM usability
    TYPE(t_soil_memory),          POINTER     :: soil     ! local pointer to soil__mem for SLM usability
    TYPE(t_sb_memory),            POINTER     :: sb       ! local pointer to sb__mem for SLM usability
    ! Get local variables from options argument
    iblk    = options%iblk
    ics     = options%ics
    ice     = options%ice
    nc      = options%nc
    ! ---------------------------
    ! 0.4 Process Activity, Debug Option
!    IF (.NOT. tile%Is_process_active(????)) RETURN   ! ask for all processes considered here ?
!    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')
    ! ---------------------------
    ! 0.5 Get Memory
    ! Get process config
!    dsl4jsb_Get_config(PHENO_)
    ! Get process memories
    dsl4jsb_Get_memory(ASSIMI_)
    dsl4jsb_Get_memory(VEG_)
    dsl4jsb_Get_memory(PHENO_)                            ! pheno__mem => pheno_memory
    dsl4jsb_Get_memory(SOIL_)                             ! soil__mem => soil_memory
    dsl4jsb_Get_memory(SB_)                               ! sb__mem => sb_memory
    ! Get process variables (Set pointers to variables in memory)
    model  => Get_model(tile%owner_model_id)
    lctlib => model%lctlib(tile%lcts%lib_id)
    ! ---------------------------
    dsl4jsb_Get_var2D_onChunk(PHENO_, growing_season)


    ! ------------------------------------------------------------------------------------------------------------
    ! Go science
    ! ------------------------------------------------------------------------------------------------------------

!! *slmtemp
    assimi  => assimi__mem 
    veg     => veg__mem 
    soil    => soil__mem  
    sb      => sb__mem  

    lun = 100 * itile
 
    SELECT CASE(TRIM(base_config%slm_model_name))
      CASE('land','plant')
        !> 0.9 run this routine ONLY if 'special case: fluxnet output' is not active, or active but not finished
        !! 
        !! same for write_output_flux
        IF (fluxnet_output_status < 3) THEN   ! this switch includes the rest of this routine
       
           !> 1.0 check output start & end, and write pools if appropriate with output_interval
           !! 
           !! 
           IF( time_y >= base_config%output_start_first_day_year .AND. time_y <= base_config%output_end_last_day_year) THEN
             IF( MODULO(time_d, output_interval_pool) < eps4 .OR. &
                & ABS(output_interval_pool - MODULO(time_d, output_interval_pool)) < eps4) THEN
              write(lun,'(13F12.6)') time_d,veg%pool%leaf%carbon%ptr,&
                veg%pool%fine_root%carbon%ptr,veg%pool%sap_wood%carbon%ptr,&
                veg%pool%heart_wood%carbon%ptr,veg%pool%coarse_root%carbon%ptr,&
                veg%pool%fruit%carbon%ptr,veg%pool%labile%carbon%ptr,&
                veg%pool%reserve%carbon%ptr,veg%pool%seed_bed%carbon%ptr,&
                veg%total_biomass%carbon%ptr,veg%total_biomass%carbon13%ptr,veg%total_biomass%carbon14%ptr
       
              dC13leaf        = calc_delta_C13(veg%pool%leaf%carbon%ptr,veg%pool%leaf%carbon13%ptr)
              dC13fine_root   = calc_delta_C13(veg%pool%fine_root%carbon%ptr,veg%pool%fine_root%carbon13%ptr)
              dC13sap_wood    = calc_delta_C13(veg%pool%sap_wood%carbon%ptr,veg%pool%sap_wood%carbon13%ptr)
              dC13heart_wood  = calc_delta_C13(veg%pool%heart_wood%carbon%ptr,veg%pool%heart_wood%carbon13%ptr)
              dC13coarse_root = calc_delta_C13(veg%pool%coarse_root%carbon%ptr,veg%pool%coarse_root%carbon13%ptr)
              dC13fruit       = calc_delta_C13(veg%pool%fruit%carbon%ptr,veg%pool%fruit%carbon13%ptr)
              dC13labile      = calc_delta_C13(veg%pool%labile%carbon%ptr,veg%pool%labile%carbon13%ptr)
              dC13reserve     = calc_delta_C13(veg%pool%reserve%carbon%ptr,veg%pool%reserve%carbon13%ptr)
              dC13seed_bed    = calc_delta_C13(veg%pool%seed_bed%carbon%ptr,veg%pool%seed_bed%carbon13%ptr)
              lun=lun+1
              write(lun,'(19F18.6)') time_d,dC13leaf,dC13fine_root,dC13sap_wood,dC13heart_wood,&
                dC13coarse_root,dC13fruit,dC13labile,dC13reserve,dC13seed_bed, & 
                calc_delta_C14(veg%pool%leaf%carbon%ptr,veg%pool%leaf%carbon14%ptr,dC13leaf),&
                calc_delta_C14(veg%pool%fine_root%carbon%ptr,veg%pool%fine_root%carbon14%ptr,dC13fine_root),&
                calc_delta_C14(veg%pool%sap_wood%carbon%ptr,veg%pool%sap_wood%carbon14%ptr,dC13sap_wood),&
                calc_delta_C14(veg%pool%heart_wood%carbon%ptr,veg%pool%heart_wood%carbon14%ptr,dC13heart_wood),&
                calc_delta_C14(veg%pool%coarse_root%carbon%ptr,veg%pool%coarse_root%carbon14%ptr,dC13coarse_root),&
                calc_delta_C14(veg%pool%fruit%carbon%ptr,veg%pool%fruit%carbon14%ptr,dC13fruit),&
                calc_delta_C14(veg%pool%labile%carbon%ptr,veg%pool%labile%carbon14%ptr,dC13labile),&
                calc_delta_C14(veg%pool%reserve%carbon%ptr,veg%pool%reserve%carbon14%ptr,dC13reserve),&
                calc_delta_C14(veg%pool%seed_bed%carbon%ptr,veg%pool%seed_bed%carbon14%ptr,dC13seed_bed)
              lun=lun+1
              write(lun,'(12F12.6)') time_d,veg%pool%leaf%nitrogen%ptr,&
                veg%pool%fine_root%nitrogen%ptr,veg%pool%sap_wood%nitrogen%ptr,&
                veg%pool%heart_wood%nitrogen%ptr,veg%pool%coarse_root%nitrogen%ptr,&
                veg%pool%fruit%nitrogen%ptr,veg%pool%labile%nitrogen%ptr,&
                veg%pool%reserve%nitrogen%ptr,veg%pool%seed_bed%nitrogen%ptr,&
                veg%total_biomass%nitrogen%ptr,veg%total_biomass%nitrogen15%ptr
              lun=lun+1
              write(lun,'(11F12.6)') time_d,veg%pool%leaf%phosphorus%ptr,&
                veg%pool%fine_root%phosphorus%ptr,veg%pool%sap_wood%phosphorus%ptr,&
                veg%pool%heart_wood%phosphorus%ptr,veg%pool%coarse_root%phosphorus%ptr,&
                veg%pool%fruit%phosphorus%ptr,veg%pool%labile%phosphorus%ptr,&
                veg%pool%reserve%phosphorus%ptr,veg%pool%seed_bed%phosphorus%ptr,veg%total_biomass%phosphorus%ptr
              lun=lun+1
              write(lun,'(9F12.6)') time_d, soil%w_soil_root%ptr(1,1), soil%w_soil_root_pot%ptr(1,1),assimi%beta_soil_gs%ptr,assimi%beta_soil_ps%ptr
              lun=lun+1
              write(lun,'(20F12.6)') time_d, sb%pool%polymeric_litter%nitrogen%ptr(1,1,1), sb%pool%polymeric_litter%phosphorus%ptr(1,1,1), & 
                   sb%nh4_solute%ptr(1,1,1), sb%no3_solute%ptr(1,1,1), sb%po4_solute%ptr(1,1,1)
              
              lun=lun+1
              DO icanopy=1, assimi_config%ncanopy
                write(lun,'(8F12.6)') time_d, REAL(icanopy), veg%fn_et_cl%ptr(icanopy,1), &
                                      veg%fn_rub_cl%ptr(icanopy,1), veg%fn_chl_cl%ptr(icanopy,1),&
                                      veg%fn_pepc_cl%ptr(icanopy,1), veg%fn_oth_cl%ptr(icanopy,1),veg%leaf_nitrogen_cl%ptr(icanopy,1)
              END DO
       
              IF(lctlib%growthform == itree) THEN 
                 IF (veg%pool%sap_wood%carbon%ptr > eps8) THEN
                    hlp1 = veg%lai%ptr/(veg%pool%sap_wood%carbon%ptr / veg%height%ptr / lctlib%wood_density)
                 ELSE
                    hlp1 = 0._wp
                 ENDIF
              ELSE
                 IF (veg%pool%leaf%carbon%ptr>eps8) THEN
                    hlp1 = veg%pool%sap_wood%carbon%ptr/veg%pool%leaf%carbon%ptr
                 ELSE
                    hlp1 = 0._wp
                 ENDIF
              ENDIF
              IF (veg%pool%fine_root%carbon%ptr > eps8) THEN
                 l2r = veg%pool%leaf%carbon%ptr/veg%pool%fine_root%carbon%ptr
              ELSE
                 l2r = 0._wp
              ENDIF
              IF (veg%pool%sap_wood%carbon%ptr>eps8) THEN
                 s2c =  veg%pool%coarse_root%carbon%ptr/veg%pool%sap_wood%carbon%ptr
              ELSE
                 s2c = 0._wp
              ENDIF
              IF ((veg%pool%leaf%carbon%ptr+veg%pool%sap_wood%carbon%ptr+veg%pool%heart_wood%carbon%ptr+&
                   veg%pool%fine_root%carbon%ptr+veg%pool%coarse_root%carbon%ptr)>eps8) THEN
                 bcf =  (veg%pool%fine_root%carbon%ptr+veg%pool%coarse_root%carbon%ptr)/ & 
                (veg%pool%leaf%carbon%ptr+veg%pool%sap_wood%carbon%ptr+veg%pool%heart_wood%carbon%ptr+&
                   veg%pool%fine_root%carbon%ptr+veg%pool%coarse_root%carbon%ptr)
              ELSE
                 bcf = 0._wp
              ENDIF
              ! transition from bool to 0/1 for growing season
              IF (growing_season(1)) THEN
                growing_season_float = 1.0_wp
                ELSE
                growing_season_float = 0.0_wp
              END IF 
              lun=lun+1
              write(lun,'(22F20.6)') time_d, veg%height%ptr, veg%diameter%ptr, veg%lai%ptr, veg%dens_ind%ptr, &
                hlp1,l2r,s2c, bcf, &
                veg%target_cn_leaf%ptr, veg%unit_npp_talloc_mavg%ptr, veg%unit_uptake_n_act_talloc_mavg%ptr,&
                veg%growth_cn_talloc_mavg%ptr, veg%npp_talloc_mavg%ptr, &
                growing_season_float, & 
                veg%cost_n_uptake_root%ptr, veg%cost_n_fixation_tuptake_mavg%ptr, &
                veg%unit_transpiration_talloc_mavg%ptr * 1.e10_wp, veg%growth_req_n_tuptake_mavg%ptr
             END IF
           END IF
        END IF  ! on/off switch specific to 'special case: fluxnet output'
      CASE('soil')
        IF( time_y >= base_config%output_start_first_day_year .AND. time_y <= base_config%output_end_last_day_year) THEN
          IF( MODULO(time_d, output_interval_pool) < eps4 .OR. &
             & ABS(output_interval_pool - MODULO(time_d, output_interval_pool)) < eps4) THEN

             write(lun,'(100F12.6)') time_d,sb%pool%soluable_litter%carbon%ptr(ics:ice,:,iblk),&
                                  sb%pool%polymeric_litter%carbon%ptr(ics:ice,:,iblk),sb%pool%woody_litter%carbon%ptr(ics:ice,:,iblk),&
                                  sb%pool%dom%carbon%ptr(ics:ice,:,iblk),sb%pool%dom_assoc%carbon%ptr(ics:ice,:,iblk),&
                                  sb%pool%microbial%carbon%ptr(ics:ice,:,iblk),sb%pool%residue%carbon%ptr(ics:ice,:,iblk),&
                                  sb%pool%residue_assoc%carbon%ptr(ics:ice,:,iblk)
             lun=lun+1
             lun=lun+1
             lun=lun+1
             write(lun,'(100F12.6)') time_d,sb%pool%soluable_litter%nitrogen%ptr(ics:ice,:,iblk),&
                                  sb%pool%polymeric_litter%nitrogen%ptr(ics:ice,:,iblk),sb%pool%woody_litter%nitrogen%ptr(ics:ice,:,iblk),&
                                  sb%pool%dom%nitrogen%ptr(ics:ice,:,iblk),sb%pool%dom_assoc%nitrogen%ptr(ics:ice,:,iblk),&
                                  sb%pool%microbial%nitrogen%ptr(ics:ice,:,iblk),sb%pool%residue%nitrogen%ptr(ics:ice,:,iblk),&
                                  sb%pool%residue_assoc%nitrogen%ptr(ics:ice,:,iblk)
             lun=lun+1
             write(lun,'(100F12.6)') time_d,sb%pool%soluable_litter%phosphorus%ptr(ics:ice,:,iblk),&
                                  sb%pool%polymeric_litter%phosphorus%ptr(ics:ice,:,iblk),sb%pool%woody_litter%phosphorus%ptr(ics:ice,:,iblk),&
                                  sb%pool%dom%phosphorus%ptr(ics:ice,:,iblk),sb%pool%dom_assoc%phosphorus%ptr(ics:ice,:,iblk),&
                                  sb%pool%microbial%phosphorus%ptr(ics:ice,:,iblk),sb%pool%residue%phosphorus%ptr(ics:ice,:,iblk),&
                                  sb%pool%residue_assoc%phosphorus%ptr(ics:ice,:,iblk)
             
          ENDIF
       ENDIF 
    END SELECT
  END SUBROUTINE write_output_pool

  !-----------------------------------------------------------------------------------------------------
  !> Fluxes: write data to files in a certain interval and accumulate fluxes
  !! 
  !! write-interval accordingly to namelist
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE write_output_flux(tile, time_d, time_y, itile)

    USE mo_jsb_process_class,     ONLY: ATM2LAND_, LAND2ATM_, SEB_, TURB_, SOIL_, HYDRO_, VEG_, HD_, RAD_, ASSIMI_, PHENO_, SB_
    USE mo_base_memory,           ONLY: t_slm_tile
    USE mo_base_constants,        ONLY: dtime
    USE mo_rad_constants,         ONLY: rad2ppfd
    USE mo_shared_isotope,        ONLY: calc_delta_C13,calc_delta_C14

    ! Use of process memories
    dsl4jsb_Use_memory(ATM2LAND_)   !  USE mo_atm2land_memory_class, ONLY: t_atm2land_memory
    dsl4jsb_Use_memory(LAND2ATM_)   !  USE mo_land2atm_memory_class, ONLY: t_land2atm_memory
    dsl4jsb_Use_memory(ASSIMI_)     !  USE mo_assimi_memory_class,   ONLY: t_assimi_memory
    dsl4jsb_Use_memory(VEG_)        !  USE mo_veg_memory_class,      ONLY: t_veg_memory
    dsl4jsb_Use_memory(RAD_)        !  USE mo_rad_memory_class,      ONLY: t_rad_memory
    dsl4jsb_Use_memory(SOIL_)       !  USE mo_soil_memory_class,     ONLY: t_soil_memory
    dsl4jsb_Use_memory(SB_)         !  USE mo_sb_memory_class,       ONLY: t_sb_memory


    IMPLICIT NONE
    ! 0.1 InOut
    TYPE(t_slm_tile),INTENT(in), TARGET  :: tile    !< one tile
    REAL(wp), INTENT(in)                 :: time_d  !< time in days
    INTEGER,  INTENT(in)                 :: time_y  !< time in years; the 1st year of a simulation == year 1 (not zero)
    INTEGER,  INTENT(in)                 :: itile   !< number of current tile
    ! 0.2 Local
    INTEGER                      :: lun             ! locical unit number
    REAL(wp), SAVE               :: dtime_one = 1.0_wp ! just 1.0
    REAL(wp)                     :: dC13gpp,dC13rg,dC13rm,hlp1
    ! 0.3 Memory allocation
    ! Declare process configuration and memory Pointers
    dsl4jsb_Def_memory(ATM2LAND_)                   !     TYPE(t_atm2land_memory), POINTER :: atm2land__mem
    dsl4jsb_Def_memory(LAND2ATM_)                   !     TYPE(t_land2atm_memory), POINTER :: land2atm__mem
    dsl4jsb_Def_memory(ASSIMI_)                     !     TYPE(t_assimi_memory),   POINTER :: assimi__mem
    dsl4jsb_Def_memory(VEG_)                        !     TYPE(t_veg_memory),      POINTER :: veg__mem
    dsl4jsb_Def_memory(RAD_)                        !     TYPE(t_rad_memory),      POINTER :: rad__mem
    dsl4jsb_Def_memory(SOIL_)                       !     TYPE(t_soil_memory),     POINTER :: soil__mem
    dsl4jsb_Def_memory(SB_)                         !     TYPE(t_sb_memory),       POINTER :: sb__mem
    TYPE(t_assimi_memory), POINTER   :: assimi      ! local pointer to assimi__mem for SLM usability
    TYPE(t_veg_memory), POINTER      :: veg         ! local pointer to veg__mem for SLM usability
    TYPE(t_rad_memory), POINTER      :: rad         ! local pointer to rad__mem for SLM usability
    TYPE(t_soil_memory), POINTER     :: soil        ! local pointer to soil__mem for SLM usability
    TYPE(t_sb_memory), POINTER       :: sb          ! local pointer to sb__mem for SLM usability
    REAL(wp), POINTER                :: ga 
    REAL(wp), POINTER                :: swpar_srf_down
    REAL(wp), POINTER                :: albedo
    REAL(wp), POINTER                :: net_radiation
    REAL(wp), POINTER                :: sensible_heat_flx
    REAL(wp), POINTER                :: latent_heat_flx
    REAL(wp), POINTER                :: ground_heat_flx
    dsl4jsb_Get_memory(ATM2LAND_)
    dsl4jsb_Get_memory(LAND2ATM_)
    dsl4jsb_Get_memory(ASSIMI_)
    dsl4jsb_Get_memory(VEG_)
    dsl4jsb_Get_memory(RAD_)
    dsl4jsb_Get_memory(SOIL_)
    dsl4jsb_Get_memory(SB_)

    assimi  => assimi__mem 
    veg     => veg__mem 
    rad     => rad__mem
    soil    => soil__mem
    sb      => sb__mem
    ga                => atm2land__mem%ga%ptr
    swpar_srf_down    => atm2land__mem%swpar_srf_down%ptr
    albedo            => land2atm__mem%albedo%ptr
    net_radiation     => land2atm__mem%net_radiation%ptr
    sensible_heat_flx => land2atm__mem%sensible_heat_flx%ptr
    latent_heat_flx   => land2atm__mem%latent_heat_flx%ptr
    ground_heat_flx   => land2atm__mem%ground_heat_flx%ptr

    SELECT CASE(TRIM(base_config%slm_model_name))
      CASE('land','plant')

        !> 0.9 run this routine ONLY if 'special case: fluxnet output' is not active, or active but not finished
        !! 
        !! same for write_output_pool
        IF (fluxnet_output_status < 3) THEN   ! this switch includes the rest of this routine
    
           !> 1.0 check output start & end, and write pools if appropriate with output_interval
           !! 
           !! 
           IF( time_y >= base_config%output_start_first_day_year .AND. time_y <= base_config%output_end_last_day_year) THEN
       
             !> 1.1 accumulate fluxes across the timesteps of one interval
             !! 
             !! for timestep-interval leave it to seconds (which is the basic timestep for fluxes), i.e. multiply by 1 \n
             !! for any interval larger timestep multiply by dtime
             IF(output_interval_flux < 1.0_wp) THEN
                 dtime_one = 1.0_wp                 ! interval = timestep
               ELSE
                 dtime_one = dtime / 1000000.0_wp      ! interval > timestep
             END IF
             gpp_acc(itile)          = gpp_acc(itile)           + assimi%gross_assimilation%ptr   * dtime_one
             growth_resp_acc(itile)  = growth_resp_acc(itile)   + veg%growth_respiration%ptr      * dtime_one
             maint_resp_acc(itile)   = maint_resp_acc(itile)    + veg%maint_respiration%ptr       * dtime_one
             n_transform_resp_acc(itile) = n_transform_resp_acc(itile) + veg%n_transform_respiration%ptr * dtime_one
             n_fixation_resp_acc(itile) = n_fixation_resp_acc(itile) + veg%n_fixation_respiration%ptr    * dtime_one
             npp_acc(itile)          = npp_acc(itile)           + veg%npp%ptr                 * dtime_one
             res_use_c_acc(itile)    = res_use_c_acc(itile)     + veg%reserve_use%carbon%ptr  * dtime_one *1000000._wp/dtime
             gpp_C13_acc(itile)          = gpp_C13_acc(itile)           + assimi%gross_assimilation_C13%ptr * dtime_one
             growth_resp_C13_acc(itile)  = growth_resp_C13_acc(itile)   + veg%growth_respiration_c13%ptr * dtime_one
             maint_resp_C13_acc(itile)   = maint_resp_C13_acc(itile)    + veg%maint_respiration_c13%ptr  * dtime_one
             gpp_C14_acc(itile)          = gpp_C14_acc(itile)           + assimi%gross_assimilation_C14%ptr * dtime_one
             growth_resp_C14_acc(itile)  = growth_resp_C14_acc(itile)   + veg%growth_respiration_c14%ptr * dtime_one
             maint_resp_C14_acc(itile)   = maint_resp_C14_acc(itile)    + veg%maint_respiration_c14%ptr  * dtime_one
             uptake_nh4_acc(itile)     = uptake_nh4_acc(itile)      + veg%uptake_nh4%ptr                * dtime_one
             uptake_no3_acc(itile)     = uptake_no3_acc(itile)      + veg%uptake_no3%ptr                * dtime_one
             recycling_n_acc(itile) = recycling_n_acc(itile)  + ( veg%recycling_leaf_n%ptr + &
                                       veg%recycling_fine_root_n%ptr + veg%recycling_heart_wood_n%ptr ) * dtime_one *1000000._wp/dtime
             res_use_n_acc(itile)    = res_use_n_acc(itile)     + veg%reserve_use%nitrogen%ptr    * dtime_one *1000000._wp/dtime
             n_fixation_acc(itile)   = n_fixation_acc(itile)    + veg%n_fixation%ptr              * dtime_one 
             n_loss_acc(itile)       = n_loss_acc(itile)        + sb%volatilisation_n%ptr(1,1,1)  * dtime_one 
             uptake_p_acc(itile)     = uptake_p_acc(itile)      + veg%uptake_p%ptr                * dtime_one
             recycling_p_acc(itile)  = recycling_p_acc(itile)  + ( veg%recycling_leaf_p%ptr + & 
                                       veg%recycling_fine_root_p%ptr + veg%recycling_heart_wood_p%ptr ) * dtime_one *1000000._wp/dtime
             res_use_p_acc(itile)    = res_use_p_acc(itile)     + veg%reserve_use%phosphorus%ptr  * dtime_one *1000000._wp/dtime
             transpiration_acc(itile)= transpiration_acc(itile) + veg%transpiration%ptr           * dtime_one 
             drainage_acc(itile)     = drainage_acc(itile)      + soil%drainage%ptr(1,1)          * dtime_one 
             gc_acc(itile)           = gc_acc(itile)            + assimi%canopy_cond%ptr                                                      
             ga_acc(itile)           = ga_acc(itile)            + ga                                                   
             albedo_acc(itile)       = albedo_acc(itile)        + albedo
             ppfd_acc(itile)         = ppfd_acc(itile)          + swpar_srf_down * rad2ppfd  * dtime_one               
             rnet_acc(itile)         = rnet_acc(itile)          + net_radiation              * dtime_one               
             qh_acc(itile)           = qh_acc(itile)            + sensible_heat_flx          * dtime_one               
             qle_acc(itile)          = qle_acc(itile)           + latent_heat_flx            * dtime_one               
             qg_acc(itile)           = qg_acc(itile)            + ground_heat_flx            * dtime_one               
             appfd_acc(itile)        = appfd_acc(itile)         + rad%appfd%ptr              * dtime_one               
             ndat_acc(itile)         = ndat_acc(itile)          + 1 
       
             !> 1.2 write fluxes and set to zero
             !!
             !!  
             IF( MODULO(time_d, output_interval_flux) < eps4 .OR. &
                & ABS(output_interval_flux - MODULO(time_d, output_interval_flux)) < eps4) THEN
       
               IF(ndat_acc(itile) > 1) THEN
                 gc_acc(itile) = gc_acc(itile) / ndat_acc(itile)
                 ga_acc(itile) = ga_acc(itile) / ndat_acc(itile)
                 albedo_acc(itile) = albedo_acc(itile) / ndat_acc(itile)
               ENDIF
       
               ! write
               lun = (100 * itile) + 30
               write(lun,'(8F12.6)') time_d, gpp_acc(itile), growth_resp_acc(itile), maint_resp_acc(itile), &
                                     n_transform_resp_acc(itile),n_fixation_resp_acc(itile), & 
                                     npp_acc(itile), res_use_c_acc(itile)
               lun=lun+1
               dC13gpp = calc_delta_C13(gpp_acc(itile),gpp_C13_acc(itile))
               dC13rg = calc_delta_C13(growth_resp_acc(itile),growth_resp_C13_acc(itile))
               dC13rm = calc_delta_C13(maint_resp_acc(itile),maint_resp_C13_acc(itile)) 
               write(lun,'(7F18.6)') time_d,dC13gpp,dC13rg,dC13rm, & 
                                     calc_delta_C14(gpp_acc(itile),gpp_C14_acc(itile),dC13gpp), & 
                                     calc_delta_C14(growth_resp_acc(itile),growth_resp_C14_acc(itile),dC13rg), & 
                                     calc_delta_C14(maint_resp_acc(itile),maint_resp_C14_acc(itile),dC13rm)
               lun=lun+1
               write(lun,'(20F12.6)') time_d, uptake_nh4_acc(itile), uptake_no3_acc(itile), recycling_n_acc(itile), & 
                                      res_use_n_acc(itile), n_fixation_acc(itile), n_loss_acc(itile)
               lun=lun+1
               write(lun,'(4F12.6)') time_d, uptake_p_acc(itile), recycling_p_acc(itile), res_use_p_acc(itile)
               lun=lun+1
               write(lun,'(12F16.8)') time_d, transpiration_acc(itile), drainage_acc(itile), gc_acc(itile), ga_acc(itile)
               lun=lun+1
               IF(ppfd_acc(itile)>eps8) THEN
                  hlp1=appfd_acc(itile)/ppfd_acc(itile)
               ELSE
                  hlp1=0.0_wp
               ENDIF
               write(lun,'(12F16.8)') time_d, rnet_acc(itile), qh_acc(itile), qle_acc(itile), -qg_acc(itile), & 
                                      albedo_acc(itile), ppfd_acc(itile), appfd_acc(itile), hlp1
               ! set zero
               gpp_acc(itile)           = 0.0_wp
               growth_resp_acc(itile)   = 0.0_wp
               maint_resp_acc(itile)    = 0.0_wp
               n_transform_resp_acc(itile) = 0.0_wp
               n_fixation_resp_acc(itile) = 0.0_wp
               npp_acc(itile)           = 0.0_wp
               res_use_c_acc(itile)     = 0.0_wp
               gpp_C13_acc(itile)       = 0.0_wp
               growth_resp_C13_acc(itile) = 0.0_wp
               maint_resp_C13_acc(itile) = 0.0_wp
               gpp_C14_acc(itile)       = 0.0_wp
               growth_resp_C14_acc(itile) = 0.0_wp
               maint_resp_C14_acc(itile) = 0.0_wp
               uptake_nh4_acc(itile)     = 0.0_wp
               uptake_no3_acc(itile)     = 0.0_wp
               recycling_n_acc(itile)   = 0.0_wp
               res_use_n_acc(itile)     = 0.0_wp
               n_fixation_acc(itile)    = 0.0_wp
               n_loss_acc(itile)        = 0.0_wp
               uptake_p_acc(itile)      = 0.0_wp
               recycling_p_acc(itile)   = 0.0_wp
               res_use_p_acc(itile)     = 0.0_wp
               transpiration_acc(itile) = 0.0_wp
               drainage_acc(itile)      = 0.0_wp
               gc_acc(itile)            = 0.0_wp
               ga_acc(itile)            = 0.0_wp
               ppfd_acc(itile)          = 0.0_wp
               rnet_acc(itile)          = 0.0_wp
               qh_acc(itile)            = 0.0_wp
               qle_acc(itile)           = 0.0_wp
               qg_acc(itile)            = 0.0_wp
               albedo_acc(itile)        = 0.0_wp
               appfd_acc(itile)          = 0.0_wp
               ndat_acc(itile)          = 0
             END IF
           END IF
        END IF  ! on/off switch specific to 'special case: fluxnet output'
      CASE('soil')
    END SELECT

  END SUBROUTINE write_output_flux

  !-----------------------------------------------------------------------------------------------------
  !> Subroutine to generate the forcing file for the stand alone soil model from 
  !!  litterfall data and SB rate modifiers 
  !! 
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE write_soil_forcing(tile, time_d, time_y, itile)

    USE mo_jsb_class,             ONLY: Get_model
    USE mo_base_memory,           ONLY: t_slm_tile, t_model
    USE mo_jsb_process_class,     ONLY: ATM2LAND_, LAND2ATM_, SEB_, TURB_, SOIL_, HYDRO_, VEG_, HD_, RAD_, ASSIMI_, PHENO_, SB_
    USE mo_veg_config_class,      ONLY: veg_config

    ! Use of process memories
    dsl4jsb_Use_memory(VEG_)        !  USE mo_veg_memory_class,      ONLY: t_veg_memory
    dsl4jsb_Use_memory(SB_)         !  USE mo_sb_memory_class,       ONLY: t_sb_memory
    dsl4jsb_Use_memory(SOIL_)       !  USE mo_soil_memory_class,     ONLY: t_soil_memory

    IMPLICIT NONE
    ! 0.1 InOut
    TYPE(t_slm_tile), INTENT(in), TARGET  :: tile         !< one tile
    REAL(wp),         INTENT(in)          :: time_d       !< time in days
    INTEGER,          INTENT(in)          :: time_y       !< time in years; the 1st year of a simulation == year 1 (not zero)
    INTEGER,          INTENT(in)          :: itile        !< number of current tile
    ! 0.2 Local
    TYPE(t_model),                POINTER :: model        ! slm model
    TYPE(t_lctlib_element),       POINTER :: lctlib       ! land-cover-type library - parameter across pft's
    INTEGER                               :: lun          ! logical unit number
    ! 0.3 Memory allocation
    ! Declare process configuration and memory Pointers
    dsl4jsb_Def_memory(VEG_)                        ! TYPE(t_veg_memory),     POINTER :: veg__mem
    dsl4jsb_Def_memory(SB_)                         ! TYPE(t_sb_memory),      POINTER :: sb__mem
    dsl4jsb_Def_memory(SOIL_)                       ! TYPE(t_soil_memory),    POINTER :: soil__mem
    TYPE(t_veg_memory),           POINTER     :: veg   ! local pointer to veg__mem for SLM usability
    TYPE(t_sb_memory),            POINTER     :: sb    ! local pointer to sb__mem for SLM usability
    TYPE(t_soil_memory),          POINTER     :: soil  ! local pointer to soil__mem for SLM usability
    dsl4jsb_Get_memory(VEG_)
    dsl4jsb_Get_memory(SB_)                         
    dsl4jsb_Get_memory(SOIL_)                         

    veg => veg__mem 
    sb  => sb__mem 
    soil  => soil__mem 
    model     => Get_model(tile%owner_model_id)
    lctlib    => model%lctlib(tile%lcts%lib_id)

    !> 1.0 check output start & end, and write litterfall and soilphysics if appropriate
    !! 
    !! identical output period as pools & fluxes
    IF( time_y >= base_config%output_start_first_day_year .AND. time_y <= base_config%output_end_last_day_year) THEN

      lun = 100 * itile + 60
      write(lun,'(56F12.6)') time_d, &
        sb%rtm_depolymerisation%ptr(1,:,1)  , & ! 1
        sb%rtm_mic_uptake%ptr(1,:,1)        , &
        sb%rtm_sorption%ptr(1,:,1)          , &
        sb%rtm_desorption%ptr(1,:,1)        , &
        sb%rtm_hsc%ptr(1,:,1)               , &
        sb%rmm_depolymerisation%ptr(1,:,1)  , & 
        sb%rmm_mic_uptake%ptr(1,:,1)        , & 
        sb%rmm_sorption%ptr(1,:,1)          , &
        sb%rmm_desorption%ptr(1,:,1)        , &
        sb%rmm_hsc%ptr(1,:,1)               , & ! 10
        soil%percolation_sl%ptr(1,:,1)

      lun=lun+1
      write(lun,'(55F12.6)') time_d, &
        ! leaf
        veg%litterfall%leaf%carbon%ptr, veg%litterfall%leaf%nitrogen%ptr, veg%litterfall%leaf%phosphorus%ptr    , &
        veg%litterfall%leaf%carbon13%ptr, veg%litterfall%leaf%carbon14%ptr, veg%litterfall%leaf%nitrogen15%ptr  , &
        ! fine_root
        veg%litterfall%fine_root%carbon%ptr, veg%litterfall%fine_root%nitrogen%ptr, veg%litterfall%fine_root%phosphorus%ptr   , &
        veg%litterfall%fine_root%carbon13%ptr, veg%litterfall%fine_root%carbon14%ptr, veg%litterfall%fine_root%nitrogen15%ptr , &
        ! coarse_root
        veg%litterfall%coarse_root%carbon%ptr, veg%litterfall%coarse_root%nitrogen%ptr, veg%litterfall%coarse_root%phosphorus%ptr   , &
        veg%litterfall%coarse_root%carbon13%ptr, veg%litterfall%coarse_root%carbon14%ptr, veg%litterfall%coarse_root%nitrogen15%ptr , &
        ! sap_wood
        veg%litterfall%sap_wood%carbon%ptr, veg%litterfall%sap_wood%nitrogen%ptr, veg%litterfall%sap_wood%phosphorus%ptr    , &
        veg%litterfall%sap_wood%carbon13%ptr, veg%litterfall%sap_wood%carbon14%ptr, veg%litterfall%sap_wood%nitrogen15%ptr  , &
        ! heart_wood
        veg%litterfall%heart_wood%carbon%ptr, veg%litterfall%heart_wood%nitrogen%ptr, veg%litterfall%heart_wood%phosphorus%ptr    , &
        veg%litterfall%heart_wood%carbon13%ptr, veg%litterfall%heart_wood%carbon14%ptr, veg%litterfall%heart_wood%nitrogen15%ptr  , &
        ! labile
        veg%litterfall%labile%carbon%ptr, veg%litterfall%labile%nitrogen%ptr, veg%litterfall%labile%phosphorus%ptr    , &
        veg%litterfall%labile%carbon13%ptr, veg%litterfall%labile%carbon14%ptr, veg%litterfall%labile%nitrogen15%ptr  , &
        ! reserve
        veg%litterfall%reserve%carbon%ptr, veg%litterfall%reserve%nitrogen%ptr, veg%litterfall%reserve%phosphorus%ptr   , &
        veg%litterfall%reserve%carbon13%ptr, veg%litterfall%reserve%carbon14%ptr, veg%litterfall%reserve%nitrogen15%ptr , &
        ! fruit
        veg%litterfall%fruit%carbon%ptr, veg%litterfall%fruit%nitrogen%ptr, veg%litterfall%fruit%phosphorus%ptr   , &
        veg%litterfall%fruit%carbon13%ptr, veg%litterfall%fruit%carbon14%ptr, veg%litterfall%fruit%nitrogen15%ptr , &
        ! seed_bed
        veg%litterfall%seed_bed%carbon%ptr, veg%litterfall%seed_bed%nitrogen%ptr, veg%litterfall%seed_bed%phosphorus%ptr  , &
        veg%litterfall%seed_bed%carbon13%ptr, veg%litterfall%seed_bed%carbon14%ptr, veg%litterfall%seed_bed%nitrogen15%ptr

    END IF 
  END SUBROUTINE write_soil_forcing 

  !-----------------------------------------------------------------------------------------------------
  !> Pools and Fluxes: close files
  !! 
  !! at this point of the simulation all data has been written
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE close_output_files(itile)

    IMPLICIT NONE
    ! 0.1 InOut
    INTEGER, INTENT(in)  :: itile   !< number of current tile 1-...
    ! 0.2 Local
    INTEGER :: lun, lun_start, lun_end
    LOGICAL :: opened

    ! 1.0 pools
    lun_start = 100 * itile
    lun_end = lun_start + 29
    DO lun = lun_start,lun_end
       INQUIRE(lun,opened=opened)
       IF(opened) CLOSE(lun) 
    ENDDO

    ! 2.0 fluxes
    lun_start = (100 * itile) + 30
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
  SUBROUTINE close_output_files_soil_forcing(itile)

    IMPLICIT NONE
    ! 0.1 InOut
    INTEGER, INTENT(in)  :: itile   !< number of current tile 1-...
    ! 0.2 Local
    INTEGER :: lun

    ! 1.0 close files
    lun = 100 * itile + 60
    CLOSE(lun) ! 0
    lun=lun+1
    CLOSE(lun) ! 1

  END SUBROUTINE close_output_files_soil_forcing

END MODULE mo_shared_output
