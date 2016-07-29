!>
!! @brief 
!!  driver for the stand-alone land model \n
!!  kind of replacement for JSBACH4
!!
!! @author
!!  Sönke Zaehle, MPI-BGC Jena 
!!
!! @par Revision History
!!  First version by Sönke Zaehle (2015-12-04) \n
!!  Reworked by Jan Engel (2016 July)
!!
PROGRAM land_model_driver

  ! echam
  USE mo_kind
  ! land model
  USE mo_land_model_interface,    ONLY: land_model_interface 
  ! base
  USE mo_base_constants,          ONLY: ntile, one_year
  USE mo_base_memory,             ONLY: t_slm_tile, t_model
  USE mo_base_config_class,       ONLY: base_config  !< TYPE, EXTENDS(t_jsb_config) :: t_base_config
  USE mo_jsb_model_init,          ONLY: init_model, &          ! 'init_model' is an interface
                                        time_day, time_year, ntimestep, itimestep, timestep_prop_day 
  ! atmland
  USE mo_atmland_forcing,         ONLY: update_forcing, open_forcing_input_file, close_forcing_input_file
  ! vegetation
  USE mo_veg_output,              ONLY: init_veg_output_variables, &
                                        open_pool_output_files, &
                                        open_flux_output_files, &
                                        close_output_files


  IMPLICIT NONE

  !! 1.0 -----------------------------------------------------------------------------------------------
  !! 
  !! Model init
  !! 
  TYPE(t_model), SAVE                      :: model         !< very basic version of one JSBACH4 model
  TYPE(t_slm_tile),DIMENSION(ntile), SAVE  :: tile          !< data structure for dynamic vegetation and soil variables

  CALL init_model(model, tile(1))                           !< model init : mo_jsb_model_init 


  !! 2.0 -----------------------------------------------------------------------------------------------
  !! 
  !! Forcing and Output files init
  !! 
  ! open input files (and read the header into dummy variables); see mo_atmland_forcing
  ! Initialise output and flux accumulation, and open output files; see mo_veg_output
  SELECT CASE(TRIM(base_config%slm_model_name))
    CASE('land')
      CALL open_forcing_input_file  ! mo_atmland_forcing; only if: option_forcing == 'readfile'
      CALL init_veg_output_variables     ! mo_veg_output : interval variables & flux accumulation per tile
      CALL open_pool_output_files(1)     ! 1 == itile 
      CALL open_flux_output_files(1)     ! 1 == itile 
    CASE('plant')
      CALL open_forcing_input_file  ! mo_atmland_forcing; only if: option_forcing == 'readfile'
      CALL init_veg_output_variables     ! mo_veg_output : interval variables & flux accumulation per tile
      CALL open_pool_output_files(1)     ! 1 == itile 
      CALL open_flux_output_files(1)     ! 1 == itile 
    CASE('canopy')
      CALL open_forcing_input_file  ! mo_atmland_forcing; only if: option_forcing == 'readfile'
    CASE('test_canopy')
      ntimestep = 1
   END SELECT


  !! 3.0 -----------------------------------------------------------------------------------------------
  !! 
  !! Time loop calling the land model interface
  !! 
  DO itimestep = 1, ntimestep

    !! 3.1 advance timer
    !! 
    time_day  = time_day  + timestep_prop_day         ! time in days

    !! 3.2 set orbital parameters and update time-variant forcing
    !!     only for models: land plant soil canopy
    SELECT CASE(TRIM(base_config%slm_model_name))
      CASE('land','plant','soil','canopy')
        CALL update_forcing(time_day)     !< mo_atmland_forcing
    END SELECT

    !! 3.3 land model interface
    !! 
    CALL land_model_interface(tile(1), &
                              model, &
                              time_day, time_year, timestep_prop_day &
                             )

    !! 3.4 check timer & adjust year counter
    !! 
    IF(time_day > one_year) THEN
      time_day  = time_day - one_year   ! largest value of time ~ 365; see mo_base_constants
      time_year = time_year + 1
    END IF

  END DO  ! end time loop 


  !! 4.0 -----------------------------------------------------------------------------------------------
  !! 
  !! Close output and forcing files
  !! 
  CALL close_output_files(1)     ! 1 == itile 
  CALL close_forcing_input_file  ! only if: option_forcing == 'readfile'

END PROGRAM land_model_driver
