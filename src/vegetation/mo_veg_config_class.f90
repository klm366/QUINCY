!!
!! @brief contains structures and methods for veg config, including namelist read
!!
!! @author
!!  Jan Engel, MPI-BGC Jena
!!
!! @par Revision History
!! First version by Jan Engel (2016-03-17)
!!
!! based on the JSBACH4 mo_veg_config_class written by Reiner Schnur, MPI-M Hamburg
!!
MODULE mo_veg_config_class
  USE mo_kind,                  ONLY: wp
  USE mo_base_config,           ONLY: options
  USE mo_namelist,              ONLY: POSITIONED, open_nml, position_nml  ! read model options from file

  IMPLICIT NONE

  ! variables for reading model-options from file (namelist.plant)
  INTEGER           :: f_unit      = 0
  INTEGER           :: read_status = 0
  INTEGER, SAVE     :: nml_unit    = 0     ! copied from JSBACH3 mo_jsbach.f90 | Unit of the jsbach namelist
  ! define variables, and set their default values, used for read model-options from file, and passed to options% ... then
  CHARACTER(8) :: read_output_interval_pool  = 'monthly'   ! used for read from file
  CHARACTER(8) :: read_output_interval_flux  = 'monthly'   ! 
  INTEGER      :: output_start_date          = 1
  INTEGER      :: output_end_date            = 100
  LOGICAL      :: flag_optimal_Nfraction     = .FALSE.  ! .TRUE. .FALSE.
  LOGICAL      :: flag_optimal_Nleaf         = .FALSE.  ! .TRUE. .FALSE.
  ! forcing options; define variables that contain values as read from file, passed to options% ... then
  CHARACTER(8) :: create_forcing             = 'generate'

  PRIVATE
  PUBLIC :: new_veg_config


CONTAINS

  SUBROUTINE new_veg_config
  
    IMPLICIT NONE

    ! 
    ! 1.3.1 &PLANT_MODEL_OPTIONS
    ! 
    NAMELIST /plant_model_options/  &    ! see file namelist & READ() function below
      read_output_interval_pool, &
      read_output_interval_flux, &
      output_start_date, &
      output_end_date, &
      flag_optimal_Nfraction, &
      flag_optimal_Nleaf
    ! read the namelist
    nml_unit = open_nml  ('namelist.plant')
    f_unit = position_nml ('PLANT_MODEL_OPTIONS', nml_unit, status=read_status)
      SELECT CASE (read_status)
        CASE (POSITIONED)
          READ (f_unit, plant_model_options)
      END SELECT
    ! pass values as read from file
    options%model%output_interval_pool   = read_output_interval_pool
    options%model%output_interval_flux   = read_output_interval_flux
    options%model%output_start_date      = output_start_date
    options%model%output_end_date        = output_end_date
    options%model%flag_optimal_Nfraction = flag_optimal_Nfraction
    options%model%flag_optimal_Nleaf     = flag_optimal_Nleaf
    ! 
    ! 1.3.2 &FORCING_CTL
    ! 
    NAMELIST /forcing_ctl/  &    ! see file namelist & READ() function below
      create_forcing
    ! read the namelist
    nml_unit = open_nml  ('namelist.plant')
    f_unit = position_nml ('FORCING_CTL', nml_unit, status=read_status)
      SELECT CASE (read_status)
        CASE (POSITIONED)
          READ (f_unit, forcing_ctl)
      END SELECT
    ! pass values as read from file to the variable set in mo_base_config
    options%forcing%option_forcing = create_forcing


  END SUBROUTINE new_veg_config

END MODULE mo_veg_config_class
