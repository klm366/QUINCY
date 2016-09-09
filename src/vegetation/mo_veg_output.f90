!>
!! @brief functions / subroutines that help to create output from plant_model_driver
!!
!! @author
!!  Jan Engel, MPI-BGC Jena
!!
!! @par Revision History
!! First version by Jan Engel (2016-03-14)
!!
MODULE mo_veg_output
  !! this module may be not necessary anymore when QUINCY is linked with JSBACH !!

  USE mo_kind,                  ONLY: wp, eps8
  USE mo_base_constants,        ONLY: ntile
  USE mo_jsb_lctlib_class,      ONLY: t_lctlib_element
  USE mo_jsb_model_init,        ONLY: fluxnet_output_status   ! used for 'special case: fluxnet output'
  USE mo_veg_config_class,      ONLY: veg_config
  USE mo_assimi_config_class,   ONLY: assimi_config


  IMPLICIT NONE

  REAL(wp), SAVE                   :: output_interval_pool   ! equals 'dt/one_day' or is a mulitple of it
  REAL(wp), SAVE                   :: output_interval_flux   ! equals 'dt/one_day' or is a mulitple of it
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
  REAL(wp), DIMENSION(ntile), SAVE :: uptake_n_acc         
  REAL(wp), DIMENSION(ntile), SAVE :: n_fixation_acc         
  REAL(wp), DIMENSION(ntile), SAVE :: recycling_n_acc     
  REAL(wp), DIMENSION(ntile), SAVE :: res_use_n_acc        
  REAL(wp), DIMENSION(ntile), SAVE :: uptake_p_acc         
  REAL(wp), DIMENSION(ntile), SAVE :: recycling_p_acc     
  REAL(wp), DIMENSION(ntile), SAVE :: res_use_p_acc        
  REAL(wp), DIMENSION(ntile), SAVE :: transpiration_acc    
  REAL(wp), DIMENSION(ntile), SAVE :: gc_acc               
  REAL(wp), DIMENSION(ntile), SAVE :: ga_acc               
  REAL(wp), DIMENSION(ntile), SAVE :: ppfd_acc 
  REAL(wp), DIMENSION(ntile), SAVE :: rnet_acc 
  REAL(wp), DIMENSION(ntile), SAVE :: qh_acc 
  REAL(wp), DIMENSION(ntile), SAVE :: qle_acc 
  REAL(wp), DIMENSION(ntile), SAVE :: albedo_acc 
  REAL(wp), DIMENSION(ntile), SAVE :: apar_acc 

CONTAINS

  !-----------------------------------------------------------------------------------------------------
  !> initialize variables for interval length and flux accumulation
  !! 
  !! 
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE init_veg_output_variables
    IMPLICIT NONE
    ! interval variables
    output_interval_pool = 0._wp  ! equals 'dt/one_day' or is a mulitple of it
    output_interval_flux = 0._wp  ! equals 'dt/one_day' or is a mulitple of it
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
    uptake_n_acc      = 0._wp   
    recycling_n_acc  = 0._wp  
    res_use_n_acc     = 0._wp  
    n_fixation_acc      = 0._wp   
    uptake_p_acc      = 0._wp
    recycling_p_acc  = 0._wp  
    res_use_p_acc     = 0._wp   
    transpiration_acc = 0._wp 
    gc_acc            = 0._wp
    ga_acc            = 0._wp
    ppfd_acc          = 0._wp
    rnet_acc          = 0._wp
    qh_acc            = 0._wp
    qle_acc           = 0._wp
    albedo_acc        = 0._wp
    apar_acc         = 0._wp
  END SUBROUTINE init_veg_output_variables


  !-----------------------------------------------------------------------------------------------------
  !> the input interval string is transformed into a number, which is used for timely model output
  !! 
  !! applied for pools and flux output
  !-----------------------------------------------------------------------------------------------------
  FUNCTION calc_output_time_interval(interval_string) RESULT(interval_number)

    USE mo_base_constants,     ONLY: one_day, dt

    IMPLICIT NONE
    ! 0.1 InOut
    CHARACTER(8), INTENT(in) :: interval_string      !< string input
    REAL(wp)                 :: interval_number      !< number output

    !> 1.0 create interval_number depending on pool/flux
    !!
    SELECT CASE( TRIM(interval_string) )
      CASE('timestep')
        interval_number = dt/one_day  ! length of one time step as advanced in 'time' 
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
    INTEGER :: lun

    !> 1.0 write headers
    !!
    !! lun of pools starts at a multiple of 100
    lun = 100 * itile
    WRITE(lun,'(2a)') 'Time Leaf_C FineRoot_C SapWood_C HeartWood_C CoarseRoot_C ', & 
                     'Fruit_C Labile_C Reserve_C SeedBank_C'
    lun=lun+1
    WRITE(lun,'(4a)') 'Time Leaf_dC13 FineRoot_dC13 SapWood_dC13 HeartWood_dC13 CoarseRoot_dC13 ', & 
                     'Fruit_dC13 Labile_dC13 Reserve_dC13 SeedBank_dC13 Leaf_DC14 FineRoot_DC14 ', &
                     'SapWood_DC14 HeartWood_DC14 CoarseRoot_DC14 Fruit_DC14 Labile_DC14 ', &
                     'Reserve_DC14 SeedBank_DC14'
    lun=lun+1
    WRITE(lun,'(2a)') 'Time Leaf_N FineRoot_N SapWood_N HeartWood_N CoarseRoot_N ', & 
                     'Fruit_N Labile_N Reserve_N SeedBank_N'
    lun=lun+1
    WRITE(lun,'(2a)') 'Time Leaf_P FineRoot_P SapWood_P HeartWood_P CoarseRoot_P ', & 
                     'Fruit_P Labile_P Reserve_P SeedBank_P'
    lun=lun+1
    WRITE(lun,'(a)') 'Time SoilMoisture'
    lun=lun+1
    WRITE(lun,'(a)') 'Time LitterN MineralN'
    lun=lun+1
    WRITE(lun,'(2a)') 'Time CanopyLayer f_et f_rub f_chl f_pepc f_oth nl_tot'
    lun=lun+1
    WRITE(lun,'(a)') &
    & 'Time Height Diameter LAI ind_density La2Sa L2R S2C BCF CN_target dNPP dNup cn_plant allometry_npp growing_season b_soil_gs cost_nup cost_bnf'

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
    !! lun of fluxes starts at a multiple of 100, plus 50
    lun = (100 * itile) + 50
    WRITE(lun,'(a)') 'Time GPP GrowthResp MaintResp NTransFResp NFixResp NPP ReserveUse'
    lun=lun+1   
    WRITE(lun,'(2a)') 'Time GPP_dC13 GrowthResp_dC13 MaintResp_dC13 ', &
                      'GPP_DC14 GrowthResp_DC14 MaintResp_DC14'
    lun=lun+1
    WRITE(lun,'(a)') 'Time NUP Resorption_N ReserveUse NFixation'
    lun=lun+1
    WRITE(lun,'(a)') 'Time PUP Resorption_P ReserveUse'
    lun=lun+1
    WRITE(lun,'(a)') 'Time Transpiration GC GA'
    lun=lun+1
    WRITE(lun,'(a)') 'Time PPFD Rnet Qh Qle Albedo APAR' 
  END SUBROUTINE init_output_writing_flux

  !-----------------------------------------------------------------------------------------------------
  !> Pools: open files and write column headers
  !! 
  !! open output files, define file names for pools, write column header names
  !! 
  !! lun: start at 10, increase per tile by 10, and increase per file by 1
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
    !! it is not necessary to pass 'veg_config%output_interval_pool' as argument to open_pool_output_files
    !! as 'veg_config' is global & defined with SAVE and initialized before this function is called
    output_interval_pool = calc_output_time_interval(veg_config%output_interval_pool)

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
    fname = 'poolC_tile'//TRIM(char_itile)//'_'//TRIM(veg_config%output_interval_pool)//TRIM(opt)//'.txt'
    OPEN (UNIT=lun, FILE=fname, ACTION="write", STATUS="unknown")
    lun=lun+1
    fname = 'pooldC13dC14_tile'//TRIM(char_itile)//'_'//TRIM(veg_config%output_interval_pool)//TRIM(opt)//'.txt'
    OPEN (UNIT=lun, FILE=fname, ACTION="write", STATUS="unknown")
    lun=lun+1
    fname = 'poolN_tile'//TRIM(char_itile)//'_'//TRIM(veg_config%output_interval_pool)//TRIM(opt)//'.txt'
    OPEN (UNIT=lun, FILE=fname, ACTION="write", STATUS="unknown")
    lun=lun+1
    fname = 'poolP_tile'//TRIM(char_itile)//'_'//TRIM(veg_config%output_interval_pool)//TRIM(opt)//'.txt'
    OPEN (UNIT=lun, FILE=fname, ACTION="write", STATUS="unknown")
    lun=lun+1
    fname = 'poolH2O_tile'//TRIM(char_itile)//'_'//TRIM(veg_config%output_interval_pool)//TRIM(opt)//'.txt'
    OPEN (UNIT=lun, FILE=fname, ACTION="write", STATUS="unknown")
    lun=lun+1
    fname = 'poolsoilN_tile'//TRIM(char_itile)//'_'//TRIM(veg_config%output_interval_pool)//TRIM(opt)//'.txt'
    OPEN (UNIT=lun, FILE=fname, ACTION="write", STATUS="unknown")
    lun=lun+1
    fname = 'canopy_tile'//TRIM(char_itile)//'_'//TRIM(veg_config%output_interval_pool)//TRIM(opt)//'.txt'
    OPEN (UNIT=lun, FILE=fname, ACTION="write", STATUS="unknown")
    lun=lun+1
    fname = 'veg_diagnostics_tile'//TRIM(char_itile)//'_'//TRIM(veg_config%output_interval_pool)//TRIM(opt)//'.txt'
    OPEN (UNIT=lun, FILE=fname, ACTION="write", STATUS="unknown")

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
  !! lun: start at 510, increase per tile by 10, and increase per file by 1
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
    !! it is not necessary to pass 'veg_config%output_interval_pool' as argument to open_pool_output_files
    !! as 'veg_config' is global & defined with SAVE and initialized before this function is called
    output_interval_flux = calc_output_time_interval(veg_config%output_interval_flux)

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
    lun = (100 * itile) + 50
    fname = 'fluxC_tile'//TRIM(char_itile)//'_'//TRIM(veg_config%output_interval_flux)//TRIM(opt)//'.txt'
    OPEN (UNIT=lun, FILE=fname, ACTION="write", STATUS="unknown")
    lun=lun+1
    fname = 'fluxdC13dC14_tile'//TRIM(char_itile)//'_'//TRIM(veg_config%output_interval_flux)//TRIM(opt)//'.txt'
    OPEN (UNIT=lun, FILE=fname, ACTION="write", STATUS="unknown")
    lun=lun+1
    fname = 'fluxN_tile'//TRIM(char_itile)//'_'//TRIM(veg_config%output_interval_flux)//TRIM(opt)//'.txt'
    OPEN (UNIT=lun, FILE=fname, ACTION="write", STATUS="unknown")
    lun=lun+1
    fname = 'fluxP_tile'//TRIM(char_itile)//'_'//TRIM(veg_config%output_interval_flux)//TRIM(opt)//'.txt'
    OPEN (UNIT=lun, FILE=fname, ACTION="write", STATUS="unknown")
    lun=lun+1
    fname = 'fluxH2O_tile'//TRIM(char_itile)//'_'//TRIM(veg_config%output_interval_flux)//TRIM(opt)//'.txt'
    OPEN (UNIT=lun, FILE=fname, ACTION="write", STATUS="unknown")
    lun=lun+1
    fname = 'fluxEnergy_tile'//TRIM(char_itile)//'_'//TRIM(veg_config%output_interval_flux)//TRIM(opt)//'.txt'
    OPEN (UNIT=lun, FILE=fname, ACTION="write", STATUS="unknown")

    !> 3.0 write column headers into output files
    !! 
    !! 
    CALL init_output_writing_flux(itile)
  END SUBROUTINE open_flux_output_files

  !-----------------------------------------------------------------------------------------------------
  !> Pools: write data to files in a certain interval
  !! 
  !! write-interval accordingly to namelist
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE write_output_pool(tile, lctlib, time_d, time_y, itile)

    USE mo_base_constants,   ONLY: ncanop
    USE mo_base_memory,      ONLY: t_slm_tile, t_slm_soil
    USE mo_veg_memory,       ONLY: t_slm_plant
    USE mo_veg_constants,    ONLY: itree
    USE mo_shared_isotope,   ONLY: calc_delta_C13,calc_delta_C14

    IMPLICIT NONE
    ! 0.1 InOut
    TYPE(t_slm_tile),INTENT(in), TARGET  :: tile        !< one tile
    TYPE(t_lctlib_element),INTENT(in)    :: lctlib      !< land-cover-type library - parameters for one pft
    REAL(wp), INTENT(in)                 :: time_d      !< time in days
    INTEGER,  INTENT(in)                 :: time_y      !< time in years; the 1st year of a simulation == year 1 (not zero)
    INTEGER,  INTENT(in)                 :: itile       !< number of current tile
    ! 0.2 Local
    TYPE(t_slm_plant),POINTER        :: plant          ! plant variables per tile
    TYPE(t_slm_soil),POINTER         :: soil           ! soil variables per tile
    INTEGER                          :: lun            ! locical unit number
    INTEGER                          :: icanop         ! loop over ncanop
    REAL(wp)                         :: hlp1           ! helper variable
    REAL(wp)                         :: l2r, s2c, bcf  !helper ratio variables
    REAl(wp)                         :: growing_season_float !output variable for growing season
    REAL(wp)                         :: dC13leaf,dC13fine_root,dC13sap_wood,dC13heart_wood,dC13coarse_root, & ! delta C13 in pools 
                                        dC13fruit,dC13labile,dC13reserve,dC13seedbank

    !  bi-directional pointer, though, bi-directionality is not necessary here
    plant         => tile%plant
    soil          => tile%soil


    !> 0.9 run this routine ONLY if 'special case: fluxnet output' is not active, or active but not finished
    !! 
    !! same for write_output_flux
    IF (fluxnet_output_status < 3) THEN   ! this switch includes the rest of this routine

    !> 1.0 check output start & end, and write pools if appropriate with output_interval
    !! 
    !! 
    IF( time_y >= veg_config%output_start_first_day_year .AND. time_y <= veg_config%output_end_last_day_year) THEN
      IF( MODULO(time_d, output_interval_pool) < eps8 .OR. &
         & ABS(output_interval_pool - MODULO(time_d, output_interval_pool)) < eps8) THEN
       lun = 100 * itile
       write(lun,'(10F12.6)') time_d,plant%pool%leaf%carbon,&
         plant%pool%fine_root%carbon,plant%pool%sap_wood%carbon,&
         plant%pool%heart_wood%carbon,plant%pool%coarse_root%carbon,&
         plant%pool%fruit%carbon,plant%pool%labile%carbon,&
         plant%pool%reserve%carbon,plant%pool%seed_bed%carbon
       lun=lun+1
       dC13leaf        = calc_delta_C13(plant%pool%leaf%carbon,plant%pool%leaf%carbon13)
       dC13fine_root   = calc_delta_C13(plant%pool%fine_root%carbon,plant%pool%fine_root%carbon13)
       dC13sap_wood    = calc_delta_C13(plant%pool%sap_wood%carbon,plant%pool%sap_wood%carbon13)
       dC13heart_wood  = calc_delta_C13(plant%pool%heart_wood%carbon,plant%pool%heart_wood%carbon13)
       dC13coarse_root = calc_delta_C13(plant%pool%coarse_root%carbon,plant%pool%coarse_root%carbon13)
       dC13fruit       = calc_delta_C13(plant%pool%fruit%carbon,plant%pool%fruit%carbon13)
       dC13labile      = calc_delta_C13(plant%pool%labile%carbon,plant%pool%labile%carbon13)
       dC13reserve     = calc_delta_C13(plant%pool%reserve%carbon,plant%pool%reserve%carbon13)
       dC13seedbank    = calc_delta_C13(plant%pool%seed_bed%carbon,plant%pool%seed_bed%carbon13)
       write(lun,'(19F18.6)') time_d,dC13leaf,dC13fine_root,dC13sap_wood,dC13heart_wood,&
         dC13coarse_root,dC13fruit,dC13labile,dC13reserve,dC13seedbank, & 
         calc_delta_C14(plant%pool%leaf%carbon,plant%pool%leaf%carbon14,dC13leaf),&
         calc_delta_C14(plant%pool%fine_root%carbon,plant%pool%fine_root%carbon14,dC13fine_root),&
         calc_delta_C14(plant%pool%sap_wood%carbon,plant%pool%sap_wood%carbon14,dC13sap_wood),&
         calc_delta_C14(plant%pool%heart_wood%carbon,plant%pool%heart_wood%carbon14,dC13heart_wood),&
         calc_delta_C14(plant%pool%coarse_root%carbon,plant%pool%coarse_root%carbon14,dC13coarse_root),&
         calc_delta_C14(plant%pool%fruit%carbon,plant%pool%fruit%carbon14,dC13fruit),&
         calc_delta_C14(plant%pool%labile%carbon,plant%pool%labile%carbon14,dC13labile),&
         calc_delta_C14(plant%pool%reserve%carbon,plant%pool%reserve%carbon14,dC13reserve),&
         calc_delta_C14(plant%pool%seed_bed%carbon,plant%pool%seed_bed%carbon14,dC13seedbank)
       lun=lun+1
       write(lun,'(10F12.6)') time_d,plant%pool%leaf%nitrogen,&
         plant%pool%fine_root%nitrogen,plant%pool%sap_wood%nitrogen,&
         plant%pool%heart_wood%nitrogen,plant%pool%coarse_root%nitrogen,&
         plant%pool%fruit%nitrogen,plant%pool%labile%nitrogen,&
         plant%pool%reserve%nitrogen,plant%pool%seed_bed%nitrogen
       lun=lun+1
       write(lun,'(10F12.6)') time_d,plant%pool%leaf%phosphorus,&
         plant%pool%fine_root%phosphorus,plant%pool%sap_wood%phosphorus,&
         plant%pool%heart_wood%phosphorus,plant%pool%coarse_root%phosphorus,&
         plant%pool%fruit%phosphorus,plant%pool%labile%phosphorus,&
         plant%pool%reserve%phosphorus,plant%pool%seed_bed%phosphorus
       lun=lun+1
       write(lun,'(9F12.6)') time_d, soil%moisture
       lun=lun+1
       write(lun,'(3F12.6)') time_d, soil%pool%litter%nitrogen, soil%ammonium
       lun=lun+1
       DO icanop=1, ncanop
         write(lun,'(8F12.6)') time_d, REAL(icanop), plant%fnl_et_cl(icanop), &
                               plant%fnl_rub_cl(icanop), plant%fnl_chl_cl(icanop),&
                               plant%fnl_pepc_cl(icanop), plant%fnl_oth_cl(icanop),plant%nl_tot_cl(icanop)
       END DO

       IF(lctlib%growthform == itree) THEN 
          IF (plant%pool%sap_wood%carbon > eps8) THEN
             hlp1 = plant%lai/(plant%pool%sap_wood%carbon / plant%height / lctlib%wood_density)
          ELSE
             hlp1 = 0._wp
          ENDIF
       ELSE
          IF (plant%pool%leaf%carbon>eps8) THEN
             hlp1 = plant%pool%sap_wood%carbon/plant%pool%leaf%carbon
          ELSE
             hlp1 = 0._wp
          ENDIF
       ENDIF
       IF (plant%pool%fine_root%carbon > eps8) THEN
          l2r = plant%pool%leaf%carbon/plant%pool%fine_root%carbon
       ELSE
          l2r = 0._wp
       ENDIF
       IF (plant%pool%sap_wood%carbon>eps8) THEN
          s2c =  plant%pool%coarse_root%carbon/plant%pool%sap_wood%carbon
       ELSE
          s2c = 0._wp
       ENDIF
       IF ((plant%pool%leaf%carbon+plant%pool%sap_wood%carbon+plant%pool%heart_wood%carbon+&
            plant%pool%fine_root%carbon+plant%pool%coarse_root%carbon)>eps8) THEN
          bcf =  (plant%pool%fine_root%carbon+plant%pool%coarse_root%carbon)/ & 
         (plant%pool%leaf%carbon+plant%pool%sap_wood%carbon+plant%pool%heart_wood%carbon+&
            plant%pool%fine_root%carbon+plant%pool%coarse_root%carbon)
       ELSE
          bcf = 0._wp
       ENDIF
       ! transition from bool to 0/1 for growing season
       IF (plant%growing_season) THEN
         growing_season_float = 1.0_wp
         ELSE
         growing_season_float = 0.0_wp
       END IF 
       lun=lun+1
       write(lun,'(20F20.6)') time_d, plant%height, plant%diameter, plant%lai, plant%dens_ind, &
         hlp1,l2r,s2c, bcf, &
         plant%target_cn_leaf, plant%unit_npp_mavg, plant%unit_uptake_n_talloc_mavg,&
         plant%growth_cn_plant_mavg, plant%allometry_npp_mavg, &
         growing_season_float, plant%beta_soil_gs_mavg, &
         plant%cost_n_uptake_root_mavg, plant%cost_n_fixation_mavg
      END IF
    END IF
    END IF  ! on/off switch specific to 'special case: fluxnet output'
  END SUBROUTINE write_output_pool

  !-----------------------------------------------------------------------------------------------------
  !> Fluxes: write data to files in a certain interval and accumulate fluxes
  !! 
  !! write-interval accordingly to namelist
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE write_output_flux(tile, time_d, time_y, itile)

    USE mo_jsb_process_class,     ONLY: ATM2LAND_, LAND2ATM_, SEB_, TURB_, SOIL_, HYDRO_, VEG_, HD_, RAD_, ASSIMI_, PHENO_
    USE mo_base_memory,           ONLY: t_slm_tile
    USE mo_veg_memory,            ONLY: t_slm_plant
    USE mo_atm2land_memory_class, ONLY: t_atm2land_memory
    USE mo_land2atm_memory_class, ONLY: t_land2atm_memory

    USE mo_base_constants,        ONLY: dt
    USE mo_atmland_constants,     ONLY: rad2par
    USE mo_shared_isotope,        ONLY: calc_delta_C13,calc_delta_C14

    IMPLICIT NONE
    ! 0.1 InOut
    TYPE(t_slm_tile),INTENT(in), TARGET  :: tile    !< one tile
    REAL(wp), INTENT(in)                 :: time_d  !< time in days
    INTEGER,  INTENT(in)                 :: time_y  !< time in years; the 1st year of a simulation == year 1 (not zero)
    INTEGER,  INTENT(in)                 :: itile   !< number of current tile
    ! 0.2 Local
    TYPE(t_slm_plant),POINTER    :: plant           ! plant variables per tile
    INTEGER                      :: lun             ! locical unit number
    REAL(wp), SAVE               :: dt_one = 1.0_wp ! just 1.0
    REAL(wp)                     :: dC13gpp,dC13rg,dC13rm
    ! 0.3 Memory allocation
    TYPE(t_atm2land_memory), POINTER :: atm2land__mem   ! dsl4jsb_Def_memory(ATM2LAND_)
    TYPE(t_land2atm_memory), POINTER :: land2atm__mem
    REAL(wp), POINTER                :: ga
    REAL(wp), POINTER                :: swpar_srf_down
    REAL(wp), POINTER                :: albedo
    REAL(wp), POINTER                :: net_radiation
    REAL(wp), POINTER                :: sensible_heat_flx
    REAL(wp), POINTER                :: latent_heat_flx
    SELECT TYPE (mm => tile%mem(ATM2LAND_)%p)            ! dsl4jsb_Get_memory(ATM2LAND_)
      TYPE IS(t_atm2land_memory)
        atm2land__mem => mm
    END SELECT
    SELECT TYPE (mm => tile%mem(LAND2ATM_)%p)
      TYPE IS(t_land2atm_memory)
        land2atm__mem => mm
    END SELECT
    ga              => atm2land__mem%ga%ptr
    swpar_srf_down  => atm2land__mem%swpar_srf_down%ptr
    albedo            => land2atm__mem%albedo%ptr
    net_radiation     => land2atm__mem%net_radiation%ptr
    sensible_heat_flx => land2atm__mem%sensible_heat_flx%ptr
    latent_heat_flx   => land2atm__mem%latent_heat_flx%ptr

    !  bi-directional pointer, though, bi-directionality is not necessary here
    plant => tile%plant


    !> 0.9 run this routine ONLY if 'special case: fluxnet output' is not active, or active but not finished
    !! 
    !! same for write_output_pool
    IF (fluxnet_output_status < 3) THEN   ! this switch includes the rest of this routine

    !> 1.0 check output start & end, and write pools if appropriate with output_interval
    !! 
    !! 
    IF( time_y >= veg_config%output_start_first_day_year .AND. time_y <= veg_config%output_end_last_day_year) THEN

      !> 1.1 accumulate fluxes across the timesteps of one interval
      !! 
      !! for timestep-interval leave it to seconds (which is the basic timestep for fluxes), i.e. multiply by 1 \n
      !! for any interval larger timestep multiply by dt
      IF(output_interval_flux < 1.0_wp) THEN
          dt_one = 1.0_wp                 ! interval = timestep
        ELSE
          dt_one = dt / 1000000.0_wp      ! interval > timestep
      END IF
      gpp_acc(itile)          = gpp_acc(itile)           + plant%gpp                     * dt_one
      growth_resp_acc(itile)  = growth_resp_acc(itile)   + plant%growth_respiration      * dt_one
      maint_resp_acc(itile)   = maint_resp_acc(itile)    + plant%maint_respiration       * dt_one
      n_transform_resp_acc(itile) = n_transform_resp_acc(itile) + plant%n_transform_respiration       * dt_one
      n_fixation_resp_acc(itile) = n_fixation_resp_acc(itile) + plant%n_fixation_respiration       * dt_one
      npp_acc(itile)          = npp_acc(itile)           + plant%npp                     * dt_one
      res_use_c_acc(itile)    = res_use_c_acc(itile)     + plant%reserve_use%carbon      * dt_one *1000000._wp/dt
      gpp_C13_acc(itile)          = gpp_C13_acc(itile)           + plant%gpp_c13                * dt_one
      growth_resp_C13_acc(itile)  = growth_resp_C13_acc(itile)   + plant%growth_respiration_c13 * dt_one
      maint_resp_C13_acc(itile)   = maint_resp_C13_acc(itile)    + plant%maint_respiration_c13  * dt_one
      gpp_C14_acc(itile)          = gpp_C14_acc(itile)           + plant%gpp_c14                * dt_one
      growth_resp_C14_acc(itile)  = growth_resp_C14_acc(itile)   + plant%growth_respiration_c14 * dt_one
      maint_resp_C14_acc(itile)   = maint_resp_C14_acc(itile)    + plant%maint_respiration_c14  * dt_one
      uptake_n_acc(itile)     = uptake_n_acc(itile)      + plant%uptake_n                * dt_one
      recycling_n_acc(itile) = recycling_n_acc(itile)  + ( plant%recycling_leaf_n + &
                                plant%recycling_fine_root_n + plant%recycling_sap_wood_n )  * dt_one *1000000._wp/dt
      res_use_n_acc(itile)    = res_use_n_acc(itile)     + plant%reserve_use%nitrogen    * dt_one *1000000._wp/dt
      n_fixation_acc(itile)   = n_fixation_acc(itile)    + plant%n_fixation              * dt_one 
      uptake_p_acc(itile)     = uptake_p_acc(itile)      + plant%uptake_p                * dt_one
      recycling_p_acc(itile)  = recycling_p_acc(itile)  + ( plant%recycling_leaf_p + & 
                                plant%recycling_fine_root_p + plant%recycling_sap_wood_p ) * dt_one *1000000._wp/dt
      res_use_p_acc(itile)    = res_use_p_acc(itile)     + plant%reserve_use%phosphorus  * dt_one *1000000._wp/dt
      transpiration_acc(itile)= transpiration_acc(itile) + plant%transpiration           * dt_one 
      gc_acc(itile)           = plant%gc                                                      
      ga_acc(itile)           = ga                                                   
      albedo_acc(itile)       = albedo
      ppfd_acc(itile)         = ppfd_acc(itile)          + swpar_srf_down * rad2par   * dt_one               
      rnet_acc(itile)         = rnet_acc(itile)          + net_radiation              * dt_one               
      qh_acc(itile)           = qh_acc(itile)            + sensible_heat_flx          * dt_one               
      qle_acc(itile)          = qle_acc(itile)           + latent_heat_flx            * dt_one               
      apar_acc(itile)         = apar_acc(itile)          + plant%apar                 * dt_one               

      !> 1.2 write fluxes and set to zero
      !!
      !!  
      IF( MODULO(time_d, output_interval_flux) < eps8 .OR. &
         & ABS(output_interval_flux - MODULO(time_d, output_interval_flux)) < eps8) THEN
        ! write
        lun = (100 * itile) + 50
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
        write(lun,'(5F12.6)') time_d, uptake_n_acc(itile), recycling_n_acc(itile), res_use_n_acc(itile), &
                              n_fixation_acc(itile)
        lun=lun+1
        write(lun,'(4F12.6)') time_d, uptake_p_acc(itile), recycling_p_acc(itile), res_use_p_acc(itile)
        lun=lun+1
        write(lun,'(12F16.8)') time_d, transpiration_acc(itile), gc_acc(itile), ga_acc(itile)
        lun=lun+1
        write(lun,'(12F16.8)') time_d, ppfd_acc(itile), rnet_acc(itile), qh_acc(itile), qle_acc(itile), & 
                               albedo_acc(itile), apar_acc(itile)
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
        uptake_n_acc(itile)      = 0.0_wp
        recycling_n_acc(itile)   = 0.0_wp
        res_use_n_acc(itile)     = 0.0_wp
        n_fixation_acc(itile)    = 0.0_wp
        uptake_p_acc(itile)      = 0.0_wp
        recycling_p_acc(itile)   = 0.0_wp
        res_use_p_acc(itile)     = 0.0_wp
        transpiration_acc(itile) = 0.0_wp
        gc_acc(itile)            = 0.0_wp
        ga_acc(itile)            = 0.0_wp
        ppfd_acc(itile)          = 0.0_wp
        rnet_acc(itile)          = 0.0_wp
        qh_acc(itile)            = 0.0_wp
        qle_acc(itile)           = 0.0_wp
        albedo_acc(itile)        = 0.0_wp
        apar_acc(itile)          = 0.0_wp
      END IF
    END IF
    END IF  ! on/off switch specific to 'special case: fluxnet output'
  END SUBROUTINE write_output_flux

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
    INTEGER :: lun

    ! 1.0 pools
    lun = 100 * itile
    CLOSE(lun) ! 0
    lun=lun+1
    CLOSE(lun) ! 1
    lun=lun+1
    CLOSE(lun) ! 2
    lun=lun+1
    CLOSE(lun) ! 3
    lun=lun+1
    CLOSE(lun) ! 4
    lun=lun+1
    CLOSE(lun) ! 5

    ! 2.0 fluxes
    lun = (100 * itile) + 50
    CLOSE(lun) ! 0
    lun=lun+1
    CLOSE(lun) ! 1
    lun=lun+1
    CLOSE(lun) ! 2
    lun=lun+1
    CLOSE(lun) ! 3
    lun=lun+1
    CLOSE(lun) ! 4
  END SUBROUTINE close_output_files

END MODULE mo_veg_output
