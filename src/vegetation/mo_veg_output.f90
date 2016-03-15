!!
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

  USE mo_kind, ONLY: wp, eps8
  USE mo_base_memory, ONLY: TileClass
  USE mo_base_constants, ONLY: ntile

  IMPLICIT NONE

CONTAINS

  FUNCTION calc_output_time_interval(interval_string, dt, one_day) RESULT(interval_number)
    ! the input interval string is transformed into a number, which is used for timely model output
    IMPLICIT NONE
    ! 0.1 Input
    CHARACTER(8), INTENT(in) :: interval_string
    REAL(wp), INTENT(in)     :: dt
    REAL(wp), INTENT(in)     :: one_day
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
    ! 1.0 Set write options and write headers
    lun=10
    WRITE(lun,'(2a)'), 'Time Leaf_C FineRoot_C SapWood_C HeartWood_C CoarseRoot_C ', & 
                     'Fruit_C Labile_C Reserve_C SeedBank_C'
    lun=lun+1
    WRITE(lun,'(2a)'), 'Time Leaf_N FineRoot_N SapWood_N HeartWood_N CoarseRoot_N ', & 
                     'Fruit_N Labile_N Reserve_N SeedBank_N'
    lun=lun+1
    WRITE(lun,'(2a)'), 'Time Leaf_P FineRoot_P SapWood_P HeartWood_P CoarseRoot_P ', & 
                     'Fruit_P Labile_P Reserve_P SeedBank_P'
  END SUBROUTINE init_output_writing_pool

  SUBROUTINE init_output_writing_flux
    ! write column headers of output files for fluxes
    IMPLICIT NONE
    ! 0.3 Local
    INTEGER :: lun = 30
    ! 1.0 Set write options
    lun=30
    WRITE(lun,'(a)'), 'Time GPP GrowthResp NPP NUP PUP LAI'
  END SUBROUTINE init_output_writing_flux

  SUBROUTINE write_output_pool(time, output_interval_pool, tile)
    ! write file output for pools at a certain interval accordingly to input
    IMPLICIT NONE
    ! 0.1 Input
    REAL(wp), INTENT(in) :: time
    REAL(wp), INTENT(in) :: output_interval_pool
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
    END IF
  END SUBROUTINE write_output_pool

  SUBROUTINE write_output_flux(time, output_interval_flux, tile)
    ! write file output for fluxes at a certain interval accordingly to input
    IMPLICIT NONE
    ! 0.1 Input
    REAL(wp), INTENT(in) :: time
    REAL(wp), INTENT(in) :: output_interval_flux
    TYPE(TileClass), DIMENSION(ntile), INTENT(in) :: tile
    ! 0.2 Output
    ! 
    ! 0.3 Local
    INTEGER :: lun = 30 

    IF( MODULO(time, output_interval_flux) < eps8 .OR. ABS(output_interval_flux - MODULO(time, output_interval_flux)) < eps8) THEN
     lun=30
     write(lun,'(10F12.6)') time,tile(1)%plant%flux%gpp,tile(1)%plant%flux%growth_respiration,& 
         tile(1)%plant%flux%npp,tile(1)%plant%flux%uptake_n,tile(1)%plant%flux%uptake_p,&
         tile(1)%plant%lai
    END IF
  END SUBROUTINE write_output_flux


END MODULE mo_veg_output
