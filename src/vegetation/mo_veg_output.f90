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

  USE mo_kind,               ONLY: wp, eps8
  USE mo_base_constants,     ONLY: ntile
  USE mo_veg_config_class,   ONLY: veg_config

  IMPLICIT NONE

  REAL(wp), SAVE                   :: output_interval_pool   ! equals 'dt/one_day' or is a mulitple of it
  REAL(wp), SAVE                   :: output_interval_flux   ! equals 'dt/one_day' or is a mulitple of it
  ! flux accumulation, see write_output_flux
  ! these are set to zero frequently in section 1.2 of write_output_flux()
  REAL(wp), DIMENSION(ntile), SAVE :: gpp_acc             
  REAL(wp), DIMENSION(ntile), SAVE :: growth_resp_acc      
  REAL(wp), DIMENSION(ntile), SAVE :: maint_resp_acc       
  REAL(wp), DIMENSION(ntile), SAVE :: npp_acc              
  REAL(wp), DIMENSION(ntile), SAVE :: res_use_c_acc        
  REAL(wp), DIMENSION(ntile), SAVE :: uptake_n_acc         
  REAL(wp), DIMENSION(ntile), SAVE :: resorption_n_acc     
  REAL(wp), DIMENSION(ntile), SAVE :: res_use_n_acc        
  REAL(wp), DIMENSION(ntile), SAVE :: uptake_p_acc         
  REAL(wp), DIMENSION(ntile), SAVE :: resorption_p_acc     
  REAL(wp), DIMENSION(ntile), SAVE :: res_use_p_acc        
  REAL(wp), DIMENSION(ntile), SAVE :: transpiration_acc    
  REAL(wp), DIMENSION(ntile), SAVE :: gc_acc               


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
    npp_acc           = 0._wp   
    res_use_c_acc     = 0._wp   
    uptake_n_acc      = 0._wp   
    resorption_n_acc  = 0._wp  
    res_use_n_acc     = 0._wp  
    uptake_p_acc      = 0._wp
    resorption_p_acc  = 0._wp  
    res_use_p_acc     = 0._wp   
    transpiration_acc = 0._wp 
    gc_acc            = 0._wp
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
    lun = 10 * itile
    WRITE(lun,'(2a)'), 'Time Leaf_C FineRoot_C SapWood_C HeartWood_C CoarseRoot_C ', & 
                     'Fruit_C Labile_C Reserve_C SeedBank_C'
    lun=lun+1
    WRITE(lun,'(2a)'), 'Time Leaf_N FineRoot_N SapWood_N HeartWood_N CoarseRoot_N ', & 
                     'Fruit_N Labile_N Reserve_N SeedBank_N'
    lun=lun+1
    WRITE(lun,'(2a)'), 'Time Leaf_P FineRoot_P SapWood_P HeartWood_P CoarseRoot_P ', & 
                     'Fruit_P Labile_P Reserve_P SeedBank_P'
    lun=lun+1
    WRITE(lun,'(a)'), 'Time SoilMoisture'
    lun=lun+1
    WRITE(lun,'(a)'), 'Time Height Diameter LAI La2Sa L2R'
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
    lun = 500 + (10 * itile)
    WRITE(lun,'(a)'), 'Time GPP GrowthResp MaintResp NPP ReserveUse'
    lun=lun+1
    WRITE(lun,'(a)'), 'Time NUP Resorption_N ReserveUse'
    lun=lun+1
    WRITE(lun,'(a)'), 'Time PUP Resorption_P ReserveUse'
    lun=lun+1
    WRITE(lun,'(a)'), 'Time Transpiration GC'
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
    CHARACTER(20)  :: opt           ! option according to namlist parameters
    CHARACTER(5)   :: char_itile    ! itile transferred into a character -> 6 digits == max number of itile is 99999


    !> 1.0 init
    opt = ""

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
    IF(veg_config%flag_optimal_Nfraction .AND. veg_config%flag_optimal_Nleaf) THEN
      opt = "_optNfrac_optNleaf"
      ELSE IF(veg_config%flag_optimal_Nfraction) THEN
        opt = "_optNfrac"
      ELSE IF(veg_config%flag_optimal_Nleaf) THEN
        opt = "_optNleaf"
    END IF

    !> 2.2 define open() statements for file output
    !! 
    !! 
    lun = 10 * itile
    fname = 'poolC_tile'//TRIM(char_itile)//'_'//TRIM(veg_config%output_interval_pool)//TRIM(opt)//'.txt'
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
    CHARACTER(20)  :: opt           ! option according to namlist parameters
    CHARACTER(5)   :: char_itile    ! itile transferred into a character -> 6 digits == max number of itile is 99999


    !> 1.0 init
    opt = ""

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
    IF(veg_config%flag_optimal_Nfraction .AND. veg_config%flag_optimal_Nleaf) THEN
      opt = "_optNfrac_optNleaf"
      ELSE IF(veg_config%flag_optimal_Nfraction) THEN
        opt = "_optNfrac"
      ELSE IF(veg_config%flag_optimal_Nleaf) THEN
        opt = "_optNleaf"
    END IF

    !> 2.2 define open() statements for file output
    !! 
    !! 
    lun = 500 + (10 * itile)
    fname = 'fluxC_tile'//TRIM(char_itile)//'_'//TRIM(veg_config%output_interval_flux)//TRIM(opt)//'.txt'
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

    !> 3.0 write column headers into output files
    !! 
    !! 
    CALL init_output_writing_flux(itile)
  END SUBROUTINE open_flux_output_files

  !-----------------------------------------------------------------------------------------------------
  !> Pools: write data to files in a certain interval
  !! 
  !! write intervall accordingly to namelist
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE write_output_pool(time_d, time_y, plant, veg_param, soil_moisture, itile)

    USE mo_veg_memory,       ONLY: t_slm_plant
    USE mo_veg_parameters,   ONLY: itree, t_slm_veg_param

    IMPLICIT NONE
    ! 0.1 InOut
    REAL(wp), INTENT(in) :: time_d  !< time in days
    INTEGER,  INTENT(in) :: time_y  !< time in years; the 1st year of a simulation == year 1 (not zero)
    INTEGER,  INTENT(in) :: itile   !< number of current tile
    REAL(wp), INTENT(in) :: soil_moisture !< soil moisture; tile(itile)\%soil\%moisture
    TYPE(t_slm_plant), INTENT(in)      :: plant       !< plant parameters for one tile
    TYPE(t_slm_veg_param),INTENT(in)   :: veg_param   !< data structure for vegetation parameters
    ! 0.2 Local
    INTEGER  :: lun
    REAL(wp) :: hlp1

    !> 1.0 check output start & end, and write pools if appropriate with output_interval
    !! 
    !! 
    IF( time_y >= veg_config%output_start_first_day_year .AND. time_y <= veg_config%output_end_last_day_year) THEN
      IF( MODULO(time_d, output_interval_pool) < eps8 .OR. &
         & ABS(output_interval_pool - MODULO(time_d, output_interval_pool)) < eps8) THEN
       lun = 10 * itile
       write(lun,'(10F12.6)') time_d,plant%pool%leaf%carbon,&
         plant%pool%fine_root%carbon,plant%pool%sap_wood%carbon,&
         plant%pool%heart_wood%carbon,plant%pool%coarse_root%carbon,&
         plant%pool%fruit%carbon,plant%pool%labile%carbon,&
         plant%pool%reserve%carbon,plant%pool%seed_bed%carbon
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
       write(lun,'(9F12.6)') time_d, soil_moisture
       lun=lun+1
       IF(veg_param%growthform(plant%pftid) == itree) THEN 
          hlp1 = plant%lai/(plant%pool%sap_wood%carbon / plant%height / veg_param%wood_density(plant%pftid))
       ELSE
          hlp1 = plant%pool%sap_wood%carbon/plant%pool%leaf%carbon
       ENDIF
       write(lun,'(9F12.6)') time_d, plant%height, plant%diameter, plant%lai, &
         hlp1,plant%pool%leaf%carbon/plant%pool%fine_root%carbon
      END IF
    END IF
  END SUBROUTINE write_output_pool

  !-----------------------------------------------------------------------------------------------------
  !> Fluxes: write data to files in a certain interval and accumulate fluxes
  !! 
  !! write intervall accordingly to namelist
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE write_output_flux(time_d, time_y, plant, itile)

    USE mo_veg_memory,     ONLY: t_slm_plant
    USE mo_base_constants, ONLY: dt

    IMPLICIT NONE
    ! 0.1 InOut
    REAL(wp), INTENT(in) :: time_d  !< time in days
    INTEGER,  INTENT(in) :: time_y  !< time in years; the 1st year of a simulation == year 1 (not zero)
    INTEGER,  INTENT(in) :: itile   !< number of current tile
    TYPE(t_slm_plant), INTENT(in) :: plant  !< plant parameters for one tile
    ! 0.2 Local
    INTEGER :: lun
    REAL(wp) :: dt_one = 1.0_wp


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
      gpp_acc(itile)          = gpp_acc(itile)           + plant%flux%gpp                     * dt_one
      growth_resp_acc(itile)  = growth_resp_acc(itile)   + plant%flux%growth_respiration      * dt_one
      maint_resp_acc(itile)   = maint_resp_acc(itile)    + plant%flux%maint_respiration       * dt_one
      npp_acc(itile)          = npp_acc(itile)           + plant%flux%npp                     * dt_one
      res_use_c_acc(itile)    = res_use_c_acc(itile)     + plant%flux%reserve_use%carbon      * dt_one *1000000._wp/dt
      uptake_n_acc(itile)     = uptake_n_acc(itile)      + plant%flux%uptake_n                * dt_one
      resorption_n_acc(itile) = resorption_n_acc(itile)  + plant%flux%leaf_resorption_n       * dt_one *1000000._wp/dt
      res_use_n_acc(itile)    = res_use_n_acc(itile)     + plant%flux%reserve_use%nitrogen    * dt_one *1000000._wp/dt
      uptake_p_acc(itile)     = uptake_p_acc(itile)      + plant%flux%uptake_p                * dt_one
      resorption_p_acc(itile) = resorption_p_acc(itile)  + plant%flux%leaf_resorption_p       * dt_one *1000000._wp/dt
      res_use_p_acc(itile)    = res_use_p_acc(itile)     + plant%flux%reserve_use%phosphorus  * dt_one *1000000._wp/dt
      transpiration_acc(itile)= transpiration_acc(itile) + plant%flux%transpiration           * dt_one 
      gc_acc(itile)           = plant%gc                                                      * dt_one 

      !> 1.2 write fluxes and set to zero
      !!
      !!  
      IF( MODULO(time_d, output_interval_flux) < eps8 .OR. &
         & ABS(output_interval_flux - MODULO(time_d, output_interval_flux)) < eps8) THEN
        ! write
        lun = 500 + (10 * itile)
        write(lun,'(6F22.6)') time_d, gpp_acc(itile), growth_resp_acc(itile), maint_resp_acc(itile), &
                              npp_acc(itile), res_use_c_acc(itile)
        lun=lun+1
        write(lun,'(4F22.6)') time_d, uptake_n_acc(itile), resorption_n_acc(itile), res_use_n_acc(itile)
        lun=lun+1
        write(lun,'(4F22.6)') time_d, uptake_p_acc(itile), resorption_p_acc(itile), res_use_p_acc(itile)
        lun=lun+1
        write(lun,'(12F16.8)') time_d, transpiration_acc(itile), gc_acc(itile)
        ! set zero
        gpp_acc(itile)           = 0.0_wp
        growth_resp_acc(itile)   = 0.0_wp
        maint_resp_acc(itile)    = 0.0_wp
        npp_acc(itile)           = 0.0_wp
        res_use_c_acc(itile)     = 0.0_wp
        uptake_n_acc(itile)      = 0.0_wp
        resorption_n_acc(itile)  = 0.0_wp
        res_use_n_acc(itile)     = 0.0_wp
        uptake_p_acc(itile)      = 0.0_wp
        resorption_p_acc(itile)  = 0.0_wp
        res_use_p_acc(itile)     = 0.0_wp
        transpiration_acc(itile) = 0.0_wp
        gc_acc(itile)            = 0.0_wp
      END IF
    END IF
  END SUBROUTINE write_output_flux

  ! 
  ! close files
  ! 
  SUBROUTINE close_output_files(itile)
    ! close output files after all data have been written
    IMPLICIT NONE
    ! 0.1 Input
    INTEGER, INTENT(in)  :: itile   !< number of current tile 1-...
    ! 0.2 Local
    INTEGER :: lun

    ! 1.0 pools
    lun = 10 * itile
    CLOSE(lun) ! 0
    lun=lun+1
    CLOSE(lun) ! 1
    lun=lun+1
    CLOSE(lun) ! 2
    lun=lun+1
    CLOSE(lun) ! 3
    lun=lun+1
    CLOSE(lun) ! 4

    ! 2.0 fluxes
    lun = 500 + (10 * itile)
    CLOSE(lun) ! 0
    lun=lun+1
    CLOSE(lun) ! 1
    lun=lun+1
    CLOSE(lun) ! 2
    lun=lun+1
    CLOSE(lun) ! 3
  END SUBROUTINE close_output_files

END MODULE mo_veg_output
