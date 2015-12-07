!!
!! @brief driver for the stand-alone plant growth 
!!
!! @author
!!  Sönke Zaehle, MPI-BGC Jena 
!!
!! @par Revision History
!! First version by Sönke Zaehle (2015-12-04)
!!
PROGRAM plant_model_driver
  USE mo_kind, ONLY: wp
  USE mo_base_memory
  USE mo_base_constants
  USE mo_land_constants
  USE mo_veg_parameters
  USE mo_veg_interface, ONLY: update_vegetation

  IMPLICIT NONE

  TYPE(VegetationParameterClass),SAVE :: veg_param     ! data structure for vegetation parameters
  TYPE(TileClass),DIMENSION(ntile), SAVE :: tile       ! data structure for dynamic vegetation and soil variables

  REAL(wp) :: ref_t_air = 25.0_wp                      ! temperature in deg C 
  REAL(wp) :: ref_theta = 100._wp                      ! soil water availability in %

  REAL(wp) :: t_air                                    ! air temperature in °C         

  INTEGER :: i,lun
  !
  ! 1.0 Initialisation of vegetation parameters and constants
  ! 
  CALL init_pft_param_table(veg_param)
  CALL init_veg_parameters
  !
  ! 1.1 Initialise time step
  ! 
  dt = 1
  one_year = 365.25_wp
  !
  ! 1.2 Initialise tile
  ! 
  tile(:)%soil%theta = ref_theta
  !
  ! 1.3 Set conditions for run
  ! 
  t_air = 0.0_wp + Tzero
  ! set GPP to 12 µmol / m2 / s, but convert to mol 
  tile(1)%plant%flux%gpp = 12._wp * 86400._wp / 1000000._wp * dt 

  lun=10
  WRITE(lun,'(2a)'), 'Leaf_C FineRoot_C SapWood_C HeartWood_C CoarseRoot_C ', & 
                     'Fruit_C Labile_C Reserve_C SeedBank_C'
  lun=lun+1
  WRITE(lun,'(2a)'), 'Leaf_N FineRoot_N SapWood_N HeartWood_N CoarseRoot_N ', & 
                     'Fruit_N Labile_N Reserve_N SeedBank_N'
  lun=lun+1
  WRITE(lun,'(2a)'), 'Leaf_P FineRoot_P SapWood_P HeartWood_P CoarseRoot_P ', & 
                     'Fruit_P Labile_P Reserve_P SeedBank_P'
  lun=lun+1
  WRITE(lun,'(a)'), 'GPP GrowthResp NPP NUP PUP'

  DO i = 1,365*200
     CALL update_vegetation(tile,veg_param)

     lun=10
     write(lun,'(9F12.6)') tile(1)%plant%pool%leaf%carbon,&
         tile(1)%plant%pool%fine_root%carbon,tile(1)%plant%pool%sap_wood%carbon,&
         tile(1)%plant%pool%heart_wood%carbon,tile(1)%plant%pool%coarse_root%carbon,&
         tile(1)%plant%pool%fruit%carbon,tile(1)%plant%pool%labile%carbon,&
         tile(1)%plant%pool%reserve%carbon,tile(1)%plant%pool%seed_bed%carbon
     lun=lun+1
     write(lun,'(9F12.6)') tile(1)%plant%pool%leaf%nitrogen,&
         tile(1)%plant%pool%fine_root%nitrogen,tile(1)%plant%pool%sap_wood%nitrogen,&
         tile(1)%plant%pool%heart_wood%nitrogen,tile(1)%plant%pool%coarse_root%nitrogen,&
         tile(1)%plant%pool%fruit%nitrogen,tile(1)%plant%pool%labile%nitrogen,&
         tile(1)%plant%pool%reserve%nitrogen,tile(1)%plant%pool%seed_bed%nitrogen
     lun=lun+1
     write(lun,'(9F12.6)') tile(1)%plant%pool%leaf%phosphorus,&
         tile(1)%plant%pool%fine_root%phosphorus,tile(1)%plant%pool%sap_wood%phosphorus,&
         tile(1)%plant%pool%heart_wood%phosphorus,tile(1)%plant%pool%coarse_root%phosphorus,&
         tile(1)%plant%pool%fruit%phosphorus,tile(1)%plant%pool%labile%phosphorus,&
         tile(1)%plant%pool%reserve%phosphorus,tile(1)%plant%pool%seed_bed%carbon
     lun=lun+1
     write(lun,'(9F12.6)') tile(1)%plant%flux%gpp,tile(1)%plant%flux%growth_respiration,& 
         tile(1)%plant%flux%npp,tile(1)%plant%flux%uptake_n,tile(1)%plant%flux%uptake_p

  ENDDO

END PROGRAM plant_model_driver
