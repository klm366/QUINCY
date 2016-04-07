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

  USE mo_kind, ONLY: wp, eps8
  USE mo_base_memory, ONLY: TileClass
  USE mo_base_constants, ONLY: ntile, one_day, dt
  USE mo_veg_config_class,   ONLY: veg_config

  IMPLICIT NONE

  REAL(wp), SAVE :: output_interval_pool = 0.  ! equals 'dt/one_day' or is a mulitple of it
  REAL(wp), SAVE :: output_interval_flux = 0.  ! equals 'dt/one_day' or is a mulitple of it

  REAL(wp), SAVE :: gpp_acc              = 0.  ! flux accumulation, see write_output_flux(time, tile)
  REAL(wp), SAVE :: growth_resp_acc      = 0.
  REAL(wp), SAVE :: maint_resp_acc       = 0.
  REAL(wp), SAVE :: npp_acc              = 0.
  REAL(wp), SAVE :: res_use_c_acc        = 0.
  REAL(wp), SAVE :: uptake_n_acc         = 0.
  REAL(wp), SAVE :: uptake_p_acc         = 0.

CONTAINS

  FUNCTION calc_output_time_interval(interval_string) RESULT(interval_number)
    ! the input interval string is transformed into a number, which is used for timely model output
    IMPLICIT NONE
    ! 0.1 Input
    CHARACTER(8), INTENT(in) :: interval_string
    ! 0.2 Output
    REAL(wp)                 :: interval_number 
    ! 0.3 Local
    ! 
    ! 1.0 create interval_number depending on pool/flux
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

  SUBROUTINE init_output_writing_pool
    ! write column headers of output files for pools
    IMPLICIT NONE
    ! 0.3 Local
    INTEGER :: lun = 10
    ! 1.0 write headers
    lun=10
    WRITE(lun,'(2a)'), 'Time Leaf_C FineRoot_C SapWood_C HeartWood_C CoarseRoot_C ', & 
                     'Fruit_C Labile_C Reserve_C SeedBank_C'
    lun=lun+1
    WRITE(lun,'(2a)'), 'Time Leaf_N FineRoot_N SapWood_N HeartWood_N CoarseRoot_N ', & 
                     'Fruit_N Labile_N Reserve_N SeedBank_N'
    lun=lun+1
    WRITE(lun,'(2a)'), 'Time Leaf_P FineRoot_P SapWood_P HeartWood_P CoarseRoot_P ', & 
                     'Fruit_P Labile_P Reserve_P SeedBank_P'
    !lun=lun+1
    !WRITE(lun,'(a)'), 'Time LAI'
  END SUBROUTINE init_output_writing_pool

  SUBROUTINE init_output_writing_flux
    ! write column headers of output files for fluxes
    IMPLICIT NONE
    ! 0.3 Local
    INTEGER :: lun = 30
    ! 1.0 write headers
    lun=30
    WRITE(lun,'(a)'), 'Time GPP GrowthResp MaintResp NPP ReserveUse'
    lun=lun+1
    WRITE(lun,'(a)'), 'Time NUP'
    lun=lun+1
    WRITE(lun,'(a)'), 'Time PUP'
  END SUBROUTINE init_output_writing_flux

  ! 
  ! open files and write column headers
  ! 
  SUBROUTINE open_pool_output_files
    ! open output files, define file names for pools, write column header names
    IMPLICIT NONE
    ! 0.3 Local
    INTEGER :: lun = 10
    CHARACTER(100) :: fname = "filename"
    CHARACTER(20)  :: opt = ""
    ! 
    ! 1.1 transform interval string to interval number (e.g. week -> 7)
    !   it is not necessary to pass 'veg_config%output_interval_pool' as argument to open_pool_output_files
    !   as 'veg_config' is global & defined with SAVE and initialized before this function is called
    ! 
    output_interval_pool = calc_output_time_interval(veg_config%output_interval_pool)
    ! 
    ! 2.1 get opt for filenames
    ! 
    IF(veg_config%flag_optimal_Nfraction .AND. veg_config%flag_optimal_Nleaf) THEN
      opt = "_optNfrac_optNleaf"
      ELSE IF(veg_config%flag_optimal_Nfraction) THEN
        opt = "_optNfrac"
      ELSE IF(veg_config%flag_optimal_Nleaf) THEN
        opt = "_optNleaf"
    END IF
    ! 
    ! 2.2 define open() statements for file output
    ! 
    lun=10
    fname = 'poolC_'//TRIM(veg_config%output_interval_pool)//TRIM(opt)//'.txt'
    OPEN (UNIT=lun, FILE=fname, ACTION="write", STATUS="unknown")
    lun=lun+1
    fname = 'poolN_'//TRIM(veg_config%output_interval_pool)//TRIM(opt)//'.txt'
    OPEN (UNIT=lun, FILE=fname, ACTION="write", STATUS="unknown")
    lun=lun+1
    fname = 'poolP_'//TRIM(veg_config%output_interval_pool)//TRIM(opt)//'.txt'
    OPEN (UNIT=lun, FILE=fname, ACTION="write", STATUS="unknown")
    !lun=lun+1
    !fname = 'LAI_'//TRIM(veg_config%output_interval_pool)//TRIM(opt)//'.txt'
    !OPEN (UNIT=lun, FILE=fname, ACTION="write", STATUS="unknown")
    ! 
    ! 3.0 write column headers into output files
    ! 
    CALL init_output_writing_pool
  END SUBROUTINE open_pool_output_files

  SUBROUTINE open_flux_output_files
    ! open output files, define file names for fluxes, write column header names
    IMPLICIT NONE
    ! 0.3 Local
    INTEGER :: lun = 10
    CHARACTER(100) :: fname = "filename"
    CHARACTER(20)  :: opt = ""
    ! 
    ! 1.1 transform interval string to interval number (e.g. week -> 7)
    ! 
    output_interval_flux = calc_output_time_interval(veg_config%output_interval_flux)
    ! 
    ! 2.1 get options for filenames
    ! 
    IF(veg_config%flag_optimal_Nfraction .AND. veg_config%flag_optimal_Nleaf) THEN
      opt = "_optNfrac_optNleaf"
      ELSE IF(veg_config%flag_optimal_Nfraction) THEN
        opt = "_optNfrac"
      ELSE IF(veg_config%flag_optimal_Nleaf) THEN
        opt = "_optNleaf"
    END IF
    ! 
    ! 2.2 define open() statements for file output
    ! 
    lun=30
    fname = 'fluxC_'//TRIM(veg_config%output_interval_flux)//TRIM(opt)//'.txt'
    OPEN (UNIT=lun, FILE=fname, ACTION="write", STATUS="unknown")
    lun=lun+1
    fname = 'fluxN_'//TRIM(veg_config%output_interval_flux)//TRIM(opt)//'.txt'
    OPEN (UNIT=lun, FILE=fname, ACTION="write", STATUS="unknown")
    lun=lun+1
    fname = 'fluxP_'//TRIM(veg_config%output_interval_flux)//TRIM(opt)//'.txt'
    OPEN (UNIT=lun, FILE=fname, ACTION="write", STATUS="unknown")
    ! 
    ! 3.0 write column headers into output files
    ! 
    CALL init_output_writing_flux
  END SUBROUTINE open_flux_output_files

  ! 
  ! write data to files in a certain interval and accumulate fluxes
  ! 
  SUBROUTINE write_output_pool(time, tile)
    ! write file output for pools at a certain interval accordingly to input
    IMPLICIT NONE
    ! 0.1 Input
    REAL(wp), INTENT(in) :: time
    !REAL(wp), INTENT(in) :: output_interval_pool
    TYPE(TileClass), DIMENSION(ntile), INTENT(in) :: tile
    ! 0.2 Output
    ! 
    ! 0.3 Local
    INTEGER :: lun = 10

    IF( MODULO(time, output_interval_pool) < eps8 .OR. ABS(output_interval_pool - MODULO(time, output_interval_pool)) < eps8) THEN
     lun=10
     write(lun,'(10F12.6)') time,tile(1)%plant%pool%leaf%carbon,&
         tile(1)%plant%pool%fine_root%carbon,tile(1)%plant%pool%sap_wood%carbon,&
         tile(1)%plant%pool%heart_wood%carbon,tile(1)%plant%pool%coarse_root%carbon,&
         tile(1)%plant%pool%fruit%carbon,tile(1)%plant%pool%labile%carbon,&
         tile(1)%plant%pool%reserve%carbon,tile(1)%plant%pool%seed_bed%carbon
     lun=lun+1
     write(lun,'(10F12.6)') time,tile(1)%plant%pool%leaf%nitrogen,&
         tile(1)%plant%pool%fine_root%nitrogen,tile(1)%plant%pool%sap_wood%nitrogen,&
         tile(1)%plant%pool%heart_wood%nitrogen,tile(1)%plant%pool%coarse_root%nitrogen,&
         tile(1)%plant%pool%fruit%nitrogen,tile(1)%plant%pool%labile%nitrogen,&
         tile(1)%plant%pool%reserve%nitrogen,tile(1)%plant%pool%seed_bed%nitrogen
     lun=lun+1
     write(lun,'(10F12.6)') time,tile(1)%plant%pool%leaf%phosphorus,&
         tile(1)%plant%pool%fine_root%phosphorus,tile(1)%plant%pool%sap_wood%phosphorus,&
         tile(1)%plant%pool%heart_wood%phosphorus,tile(1)%plant%pool%coarse_root%phosphorus,&
         tile(1)%plant%pool%fruit%phosphorus,tile(1)%plant%pool%labile%phosphorus,&
         tile(1)%plant%pool%reserve%phosphorus,tile(1)%plant%pool%seed_bed%carbon
     !lun=lun+1
     !write(lun,'(2F12.6)') time, tile(1)%plant%lai 
    END IF
  END SUBROUTINE write_output_pool

  SUBROUTINE write_output_flux(time, tile)
    ! write file output for fluxes at a certain interval accordingly to input
    IMPLICIT NONE
    ! 0.1 Input
    REAL(wp), INTENT(in) :: time
    !REAL(wp), INTENT(in) :: output_interval_flux
    TYPE(TileClass), DIMENSION(ntile), INTENT(in) :: tile
    ! 0.2 Output
    ! 
    ! 0.3 Local
    INTEGER :: lun = 30 
    REAL(wp) :: dt_one = 1.0

    ! 
    ! 1.0 accumulate fluxes across the timesteps of one interval
    !     for timestep-interval leave it to seconds (which is the basic timestep for fluxes), i.e. multiply by 1
    !     for any interval larger timestep multiply by dt
    ! 
    IF(output_interval_flux < 1.0_wp) THEN
        dt_one = 1.0_wp                 ! interval = timestep
      ELSE
        dt_one = dt / 1000000.0_wp      ! interval > timestep
    END IF
    gpp_acc          = gpp_acc         + tile(1)%plant%flux%gpp                * dt_one
    growth_resp_acc  = growth_resp_acc + tile(1)%plant%flux%growth_respiration * dt_one
    maint_resp_acc   = maint_resp_acc  + tile(1)%plant%flux%maint_respiration  * dt_one
    npp_acc          = npp_acc         + tile(1)%plant%flux%npp                * dt_one
    res_use_c_acc    = res_use_c_acc   + tile(1)%plant%flux%reserve_use%carbon * dt_one *1000000._wp/dt
    uptake_n_acc     = uptake_n_acc    + tile(1)%plant%flux%uptake_n           * dt_one
    uptake_p_acc     = uptake_p_acc    + tile(1)%plant%flux%uptake_p           * dt_one
    ! 
    ! 2.0 write fluxes and set to zero
    ! 
    IF( MODULO(time, output_interval_flux) < eps8 .OR. ABS(output_interval_flux - MODULO(time, output_interval_flux)) < eps8) THEN
     lun=30
     write(lun,'(6F22.6)') time, gpp_acc, growth_resp_acc, maint_resp_acc, npp_acc, res_use_c_acc
     lun=lun+1
     write(lun,'(2F22.6)') time, uptake_n_acc
     lun=lun+1
     write(lun,'(2F22.6)') time, uptake_p_acc
     gpp_acc         = 0._wp
     growth_resp_acc = 0._wp
     maint_resp_acc  = 0._wp
     npp_acc         = 0._wp
     res_use_c_acc   = 0._wp
     uptake_n_acc    = 0._wp
     uptake_p_acc    = 0._wp
    END IF
  END SUBROUTINE write_output_flux

  ! 
  ! close files
  ! 
  SUBROUTINE close_output_files
    ! close output files after all data have been written
    IMPLICIT NONE

    CLOSE(10)  ! please consider the unit names used in open_output_files
    CLOSE(11)
    CLOSE(12)
    CLOSE(30)
  END SUBROUTINE close_output_files

END MODULE mo_veg_output
