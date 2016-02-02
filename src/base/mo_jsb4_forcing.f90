!>
!! @brief Module to read and convert external driving/forcing data (expects annual files with year in filename)
!!
!! @author
!!  Julia Nabel, MPI-M Hamburg
!!
!! @par Revision History
!! Initial version: jsbach3 - mo_jsbalone_forcing.f90 -- current implementation only for daily, timestep or constant data

!TODO-check-with-RMS: Delete not required shortwave functions/ namelist flags and refactor whats left
! (The jsb3 functions heavily duplicate code!)  -- (check jsbach init header files in cosmos-landveg?)

MODULE mo_jsb4_forcing
#ifndef __NO_JSBACH__

#include "dsl_jsbach.inc"

  !>>> TODO: echam!
  USE mo_time_control,    ONLY: get_date_components, next_date, current_date, lstart, lresume, delta_time, time_step_len
  USE mo_netcdf,          ONLY: nf_check, nf_get_vara_double, nf_inq_dimid, nf_inq_dimlen, &
                                nf_inq_varid, NF_MAX_NAME, nf_get_att_text, NF_NOERR
  USE mo_time_conversion, ONLY: IMerge_HMS2Sec, day_in_year
  !<<<

  ! Apparently in echam and icon
  USE mo_math_constants,       ONLY: pi
  USE mo_kind,                 ONLY: wp
  USE mo_exception,            ONLY: message, message_text, finish
  USE mo_util_string,          ONLY: tolower, int2string
  USE mo_physical_constants,   ONLY: tmelt, p0sl_bg, grav, rd, amco2, amd, stbo, rdv, cpd, vtmpc1
  USE mo_echam_convect_tables, ONLY: tlucua, jptlucu1, jptlucu2

  ! jsb4
  USE mo_jsb_base,           ONLY: jsbach
  USE mo_jsb_model_class,    ONLY: t_jsb_model
  USE mo_jsb_tile_class,     ONLY: t_jsb_tile
  USE mo_jsb_parallel_iface, ONLY: my_process_is_stdio, my_process_is_mpi_parallel, p_io, p_bcast
  USE mo_jsb_io_netcdf,      ONLY: t_input_file, jsb_netcdf_open_input, jsb_netcdf_close_file, jsb_netcdf_read_2d_time

  dsl4jsb_Use_processes SRF
  dsl4jsb_Use_memory(SRF)

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: get_interface_variables_from_external_forcing, init_forcing

  ! --- CONSTANTS
  CHARACTER(len=*), PARAMETER :: modname = 'mo_jsbalone_forcing'

!>>>>> TODO-check-with-RMS: module for general constants / soubroutine initialising general constants?
  REAL(wp), PARAMETER :: gamma = 6.0E-3_wp  ! copied from jsbach3 - mo_jsbalone_forcing
                                            ! fixed temperature gradient of atmosphere [K/m] (see Knorr, p. 31). Note that
                                            ! for the standard ICAO-atmosphere this value is different, namely 6.5E-3 !!!
  REAL(wp), PARAMETER :: solar_const =  1361.371_wp               ! default solar constant (AMIP) [W/m2]
  REAL(wp), PARAMETER :: solar_const_vis = solar_const * 0.44_wp  ! solar constant corrected for visible rad. [W/m^2]
  REAL(wp), PARAMETER :: conv_prec_day = 1._wp/86400._wp ! copied from jsbach3 - mo_jsbalone_forcing
                                                         ! Conversion factor for precipitation in [mm/day] to [kg/m^2/s]:
                                                         ! 1 mm/day is equivalent to 1 kg /86400 s / m^2
  REAL(wp), PARAMETER :: molarMassDryAir_kg = amd * 1.e-3_wp  !  28.970e-3_wp   -  Mass of 1 mol of dry air in kg
  REAL(wp), PARAMETER :: molarMassCO2_kg = amco2 * 1.e-3_wp   ! 44.0095e-3_wp   -  Mass of 1 mol CO2 in kg

  !in jsb3: originally set in mo_jsbalone_forcing init_forcing to this constant
  REAL(wp), PARAMETER :: HeightWind = -1._wp
  REAL(wp), PARAMETER :: HeightHumidity = 10._wp

  !contained in jsb3 interface but (currently not in jsb4)
  REAL(wp), PARAMETER :: Veg_height = 10._wp

  !originally from echam mo_physc2.f90 and in jsb3 initialised with jsbalone_iniphy
  REAL(wp), PARAMETER :: cb = 5._wp
  REAL(wp), PARAMETER :: cc = 5._wp
  REAL(wp), PARAMETER :: cvdifts = 1.5_wp
  REAL(wp), PARAMETER :: ckap = 0.4_wp

  CHARACTER(len=*), PARAMETER :: name_of_unit_attribute = 'units'
!<<<<<<<

  ! expected varnames in forcing files
  CHARACTER(len=*), PARAMETER :: max_temp_varname = 'tmax'
  CHARACTER(len=*), PARAMETER :: min_temp_varname = 'tmin'
  CHARACTER(len=*), PARAMETER :: air_temp_varname = 'air_temp'
  CHARACTER(len=*), PARAMETER :: precipitation_varname = 'precip'
  CHARACTER(len=*), PARAMETER :: rel_humidity_varname = 'rel_humidity'
  CHARACTER(len=*), PARAMETER :: qair_varname = 'qair'
  CHARACTER(len=*), PARAMETER :: cloud_varname = 'cloud'
  CHARACTER(len=*), PARAMETER :: shortwave_varname = 'shortwave'
  CHARACTER(len=*), PARAMETER :: potential_shortwave_varname = 'mpot'
  CHARACTER(len=*), PARAMETER :: longwave_varname = 'longwave'
  CHARACTER(len=*), PARAMETER :: CO2_varname = 'CO2'
  CHARACTER(len=*), PARAMETER :: wspeed_varname = 'wspeed'

  ! expected units in forcing files (units of precipitation and co2 are specified in the namelist parameters)
  CHARACTER(len=*), PARAMETER :: exp_unit_temperature = 'degC'
  CHARACTER(len=*), PARAMETER :: exp_unit_radiation = 'W/m2'
  CHARACTER(len=*), PARAMETER :: exp_unit_qair = 'g/g'
  CHARACTER(len=*), PARAMETER :: exp_unit_wspeed = 'm/s'
  CHARACTER(len=*), PARAMETER :: unit_precip_mm_day = 'mm/day'
  CHARACTER(len=*), PARAMETER :: unit_precip_kg_m2_s = 'kg/m2/s'
  CHARACTER(len=*), PARAMETER :: unspecified_exp_unit = ''

  ! options controlling types of forcing data
  INTEGER,PARAMETER  :: DAILY_    = 1
  INTEGER,PARAMETER  :: CONST_    = 2
  INTEGER,PARAMETER  :: TIMESTEP_ = 3

  ! options controlling atmospheric humidity input data
  INTEGER,PARAMETER :: RH_   = 1         ! relative humidity as forcing data
  INTEGER,PARAMETER :: QAIR_ = 2         ! specific humidity as forcing data

  ! options controlling radiation input data (used for longwave and shortwave)
  INTEGER,PARAMETER :: CLOUD_    = 1     ! cloudcover as forcing data
  INTEGER,PARAMETER :: MEAN_RAD_ = 2     ! mean daily incoming radiation as input data

  ! options controlling the shortwave-radiation algorithms used
  INTEGER,PARAMETER :: SW_SCHEME_ORIGINAL_    = 1 ! Use shortwave_from_fpar_orig()  (scheme from Wolfgangs Diss.)
  INTEGER,PARAMETER :: SW_SCHEME_CCDAS_       = 2 ! Use shortwave_from_fpar_ccdas() (scheme from CCDAS)
  INTEGER,PARAMETER :: SW_SCHEME_NEW_         = 3 ! Use shortwave_from_fsw()  (Wolfgangs new scheme)
  INTEGER,PARAMETER :: SW_SCHEME_NEW_ORIG_CF_ = 4 ! Use shortwave_from_fsw()  (Wolfgangs new scheme with original conversion factor)
  INTEGER,PARAMETER :: SW_SCHEME_DIRECT_      = 5 ! Use shortwave_from_direct_shortwave() (scheme from shortwave record)

  ! --- Namelist parameters
  CHARACTER(NF_MAX_NAME) :: forcing_temp_file_prefix         !! Prefix for files with temperature forcing data [degC]
  CHARACTER(len=128)     :: forcing_temp_frequ               !! Frequency of temperature forcing data (DAILY/CONST/TIMESTEP)
  REAL(wp)               :: forcing_temp_const_tmin          !! Constant value for minimum daily temperature [degC]
  REAL(wp)               :: forcing_temp_const_tmax          !! Constant value for maximum daily temperature [degC]
  REAL(wp)               :: forcing_temp_ocean               !! Constant value for ocean cells [degC]

  CHARACTER(NF_MAX_NAME) :: forcing_precip_file_prefix       !! Prefix for files with precipitation forcing data (unit: mm/day
                                                             !! or kg/m^2/s, depending on 'forcing_precip_in_mm_per_day')
  CHARACTER(len=128)     :: forcing_precip_frequ             !! Frequency of precipitation forcing (DAILY/CONST/TIMESTEP)
  REAL(wp)               :: forcing_precip_const_precip      !! Constant value for precipitation
  REAL(wp)               :: forcing_precip_ocean             !! Constant value for ocean cells
  LOGICAL                :: forcing_precip_in_mm_per_day     !! Precipitation input in mm/day (true) or kg/m^2/s (false)

  CHARACTER(NF_MAX_NAME) :: forcing_sw_file_prefix           !! Prefix for files with shortwave radiation data
  CHARACTER(len=128)     :: forcing_sw_frequ                 !! Frequency of shortwave forcing data (DAILY/CONST/TIMESTEP)
  CHARACTER(len=128)     :: forcing_sw_type                  !! Type of input data used for forcing by shortwave radiation
  CHARACTER(len=128)     :: forcing_sw_scheme                !! Scheme to generate the  the components (PAR, NIR) of sw forcing
  REAL(wp)               :: forcing_sw_const_cloud           !! Constant value for cloudcover [percent]
  REAL(wp)               :: forcing_sw_const_shortwave       !! Constant value for downward shortwave radiation [W/m^2]
  REAL(wp)               :: forcing_sw_ocean                 !! Constant value for ocean cells

  CHARACTER(NF_MAX_NAME) :: forcing_table_sw_pot_file_prefix !! Prefix for files with potential short wave radiation data [W/m^2]
  CHARACTER(len=128)     :: forcing_table_sw_pot_frequ       !! Frequency of potential sw radiation (DAILY/CONST/TIMESTEP)
  REAL(wp)               :: forcing_table_sw_pot_const       !! Constant value for potential shortwave radiation [W/m^2]
  REAL(wp)               :: forcing_table_sw_pot_ocean       !! Constant value for ocean cells

  CHARACTER(NF_MAX_NAME) :: forcing_lw_file_prefix           !! Prefix for files with longwave radiation data [W/m^2]
  CHARACTER(len=128)     :: forcing_lw_frequ                 !! Frequency of longwave forcing data (DAILY/CONST/TIMESTEP)
  CHARACTER(len=128)     :: forcing_lw_type                  !! Type of forcing data for longwave forcing
  REAL(wp)               :: forcing_lw_const_cloud           !! Constant value for cloud cover [percent]
  REAL(wp)               :: forcing_lw_const_longwave        !! Constant value for downward longwave radiation [W/m^2]
  REAL(wp)               :: forcing_lw_ocean                 !! Constant value for ocean cells

  CHARACTER(NF_MAX_NAME) :: forcing_co2_file_prefix          !! Prefix for files with CO2 data (unit: mol(CO2)/mol(air)
                                                             !! or kg(CO2)/kg(air) or ppmv, depending on 'forcing_co2_unit'
  CHARACTER(len=128)     :: forcing_co2_frequ                !! Frequency of CO2 forcing data (DAILY/CONST)
  REAL(wp)               :: forcing_co2_const_co2            !! Constant value for CO2-concentration (unit: mol(CO2)/mol(air) or
                                                             !! kg(CO2)/kg(air) or ppmv, depending on 'forcing_co2_unit')
  REAL(wp)               :: forcing_co2_ocean                !! Constant value for ocean cells
  CHARACTER(len=128)     :: forcing_co2_unit                 !! Unit of CO2 input (PPMV/MOL_PER_MOL/KG_PER_KG)

  CHARACTER(NF_MAX_NAME) :: forcing_wind_file_prefix         !! Prefix for files with windspeed data [m/s]
  CHARACTER(len=128)     :: forcing_wind_frequ               !! Frequency of windspeed forcing data (DAILY/CONST/TIMESTEP)
  REAL(wp)               :: forcing_wind_const_wspeed        !! Constant value for windspeed [m/s]
  REAL(wp)               :: forcing_wind_ocean               !! Constant value for ocean cells

  CHARACTER(NF_MAX_NAME) :: forcing_qair_file_prefix         !! Prefix for files with qair data
  CHARACTER(len=128)     :: forcing_qair_frequ               !! Frequency of qair forcing data (DAILY/CONST)
  CHARACTER(len=128)     :: forcing_qair_type                !! Type of input data used for forcing by atm. humidity
  REAL(wp)               :: forcing_qair_const_rh            !! Constant value for relative humidity [%]
  REAL(wp)               :: forcing_qair_ocean               !! Constant value for ocean cells


  ! --- Structures for handling input data
  ! options relating forcing read from namelist
  TYPE forcing_options_type
    LOGICAL  :: air_temperature_as_timestep !! Specifies if the air temperatures are given per timestep
    INTEGER  :: type_of_qair_forcing        !! Possible values: see 'options controlling atmospheric humidity input data'
    INTEGER  :: type_of_shortwave_forcing   !! Possible values: see 'options controlling radiation input data'
    INTEGER  :: type_of_shortwave_scheme    !! Possible values: see 'options controlling the shortwave-radiation algorithms'
    INTEGER  :: type_of_longwave_forcing    !! Possible values: see 'options controlling radiation input data'
    LOGICAL  :: precip_in_mm_per_day        !! Specifies units of precipitation input data: true: mm/day; false: kg/(m^2*s)
    REAL(wp) :: conv_CO2_2_MassRatio        !! Factor for the conversion of CO2 input data to kg(CO2)/kg(dry air)
  END TYPE forcing_options_type
  TYPE(forcing_options_type), SAVE :: forcing_options

  ! information for single forcing variable and the file it is read from
  TYPE input_variable_type
    CHARACTER(NF_MAX_NAME)   :: variable_name=""       !! Name of the variable (name in netcdf-data file)
    CHARACTER(NF_MAX_NAME)   :: file_name_prefix = ""  !! Prefix of files in which the variable is located
    CHARACTER(NF_MAX_NAME)   :: expectedUnit = ""      !! Unit for the variable expected in the file
    INTEGER                  :: frequency = 0          !! Input data frequency (see 'options controlling types of forcing data')
    REAL(wp)                 :: ocean_value = 0._wp    !! Value used for ocean cells
    REAL(wp)                 :: constantValue = 0._wp  !! Value in case of frequency 'constant'
    INTEGER                  :: currentDayStart = 0    !! First time step in file corresponding to the current day
    INTEGER                  :: currentDayEnd = 0      !! Last time step in file corresponding to the current day
    INTEGER                  :: length_time_series     !! The number of time steps found in the data (only if my_process_is_stdio())
    REAL(wp), POINTER        :: data_read(:,:,:) => NULL()  !! Pointer to read data
    REAL(wp), POINTER        :: timevalues(:) => NULL()     !! The values of the time variable (only if my_process_is_stdio())
  END TYPE input_variable_type

  ! All theoretically possible input data, i.e. some are usually not used. (Input data is controlled by FORCING_CTL)
  TYPE forcing_input_type
     TYPE(input_variable_type) :: air_temp_min  !! daily minimum of air temperature at surface [Celsius]
     TYPE(input_variable_type) :: air_temp_max  !! daily maximum of air temperature at surface [Celsius]
     TYPE(input_variable_type) :: air_temp      !! air temperature [Celsius]
     TYPE(input_variable_type) :: precipitation !! daily precipitation (including snow); units may be
                                                !! [mm/day] or [kg/m^2/s] (see option  "precip_in_mm_per_day")
     TYPE(input_variable_type) :: cloud_cover   !! fraction of sky covered by clouds (from [0,1])
     TYPE(input_variable_type) :: shortwave     !! incident shortwave radiation [W/m^2]
     TYPE(input_variable_type) :: shortwave_pot !! potential (i.e. clear sky) shortwave radiation incident to surface [W/m^2]
     TYPE(input_variable_type) :: longwave      !! incident longwave radiation [W/m^2]
     TYPE(input_variable_type) :: CO2_concentr  !! CO2-concentration [mol(CO2)/mol(air)], [kg(CO2)/kg(dry Air)] or [ppmv]
     TYPE(input_variable_type) :: wind_speed    !! wind speed [m/s]
     TYPE(input_variable_type) :: rel_humidity  !! relative humidity [%]
     TYPE(input_variable_type) :: qair          !! specific humidity [kg/kg]
  END TYPE forcing_input_type
  TYPE(forcing_input_type),SAVE :: forcing_input

CONTAINS

  !>
  !! @brief provides the variables required for jsb4 interface calls in this timestep based on (namelist defined) external forcing
  !!
  !! @par History
  !! JN: used parts of jsbach3 mo_jsbalone_forcing: read_forcing / update_forcing and mo_jsbalone: update_driving
  !!
  SUBROUTINE get_interface_variables_from_external_forcing(model_id, ics, ice, nblks, sinlat, coslat, lon, &
    & t_srf_proc, fact_q_air_proc, fact_qsat_srf_proc, rough_h_srf_proc, cos_zenith_angle, CO2_concentration, &
    & t_air_K, q_air, rain, snow, wind_air, wind_10m, lw_srf_down, swvis_srf_down, &
    & swnir_srf_down, swpar_srf_down, press_srf, drag_srf, t_acoef, t_bcoef, q_acoef, q_bcoef, pch )

    INTEGER,  INTENT(IN)  :: model_id, ics, ice, nblks
    REAL(wp), INTENT(IN)  :: sinlat(:,:)
    REAL(wp), INTENT(IN)  :: coslat(:,:)
    REAL(wp), INTENT(IN)  :: lon(:,:)
    REAL(wp), INTENT(IN)  :: t_srf_proc(:,:)
    REAL(wp), INTENT(IN)  :: fact_q_air_proc(:,:)
    REAL(wp), INTENT(IN)  :: fact_qsat_srf_proc(:,:)
    REAL(wp), INTENT(IN)  :: rough_h_srf_proc(:,:)
    REAL(wp), INTENT(IN)  :: cos_zenith_angle(:,:)

    REAL(wp), INTENT(OUT) :: CO2_concentration(:,:)
    REAL(wp), INTENT(OUT) :: t_air_K(:,:)
    REAL(wp), INTENT(OUT) :: q_air(:,:)
    REAL(wp), INTENT(OUT) :: rain(:,:)
    REAL(wp), INTENT(OUT) :: snow(:,:)
    REAL(wp), INTENT(OUT) :: wind_air(:,:)
    REAL(wp), INTENT(OUT) :: wind_10m(:,:)
    REAL(wp), INTENT(OUT) :: lw_srf_down(:,:)
    REAL(wp), INTENT(OUT) :: swvis_srf_down(:,:)
    REAL(wp), INTENT(OUT) :: swnir_srf_down(:,:)
    REAL(wp), INTENT(OUT) :: swpar_srf_down(:,:)
    REAL(wp), INTENT(OUT) :: press_srf(:,:)
    REAL(wp), INTENT(OUT) :: drag_srf(:,:)
    REAL(wp), INTENT(OUT) :: t_acoef(:,:)
    REAL(wp), INTENT(OUT) :: t_bcoef(:,:)
    REAL(wp), INTENT(OUT) :: q_acoef(:,:)
    REAL(wp), INTENT(OUT) :: q_bcoef(:,:)
    REAL(wp), INTENT(OUT) :: pch(:,:)


    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':get_interface_variables_from_external_forcing'

    TYPE(t_jsb_model),  POINTER  :: model
    TYPE(t_jsb_tile), POINTER :: tile

    dsl4jsb_Def_memory(SRF)
    dsl4jsb_Real2D_onChunk :: elevation

    REAL(wp) :: hlp_r
    INTEGER  :: iblk
    INTEGER  :: day_of_year

    REAL(wp), ALLOCATABLE :: cloud_cover_data(:,:), rad_uv_down(:,:), vapor_pressure_data(:,:)
    REAL(wp), ALLOCATABLE :: precipitation_data(:,:), t_air(:,:), tmax_data(:,:), tmin_data(:,:), sat_spec_humidity(:,:)

    REAL(wp), ALLOCATABLE :: temporary_data1(:,:), temporary_data2(:,:)

    !TODO: refactor ... only required to handle unused frac_PAR_diffuse, rad_sw_down, rad_sw_down_pot_out
    REAL(wp), ALLOCATABLE :: help1(:,:), help2(:,:), help3(:,:)

    !---------------------------------------------------------------------------------
    ALLOCATE(help1(ice, nblks))
    ALLOCATE(help2(ice, nblks))
    ALLOCATE(help3(ice, nblks))

    ALLOCATE(temporary_data1(ice, nblks))
    ALLOCATE(temporary_data2(ice, nblks))

    ! --- temperature
    ! for temperature the variables used depend on the forcing frequency: tmin&tmax for daily data t_air for timestep
    ALLOCATE(t_air(ice, nblks))
    IF(forcing_options%air_temperature_as_timestep) THEN
      CALL get_data_and_update_field(nblks, ics, ice, model_id, forcing_input%air_temp, t_air)
    ELSE
      ALLOCATE(tmin_data(ice, nblks))
      ALLOCATE(tmax_data(ice, nblks))

      CALL get_data_and_update_field(nblks, ics, ice, model_id, forcing_input%air_temp_min, tmin_data)
      CALL get_data_and_update_field(nblks, ics, ice, model_id, forcing_input%air_temp_max, tmax_data)
      ! estimate temperature at time step from daily temperature maximum and temperature minimum
      CALL instantly_from_daily_temp_2(nblks, ice, lon, sinlat, coslat, tmin_data, tmax_data, next_date, t_air)
    ENDIF
    ! convert degC to Kelvin
    t_air_K = t_air + tmelt

    ! --- precipitation
    ALLOCATE(precipitation_data(ice, nblks))
    CALL get_data_and_update_field(nblks, ics, ice, model_id, forcing_input%precipitation, precipitation_data)
    ! convert if required
    IF(forcing_options%precip_in_mm_per_day) precipitation_data = precipitation_data * conv_prec_day
    ! compute rain and snow
    IF(forcing_options%air_temperature_as_timestep) THEN
      IF(forcing_input%precipitation%frequency .EQ. DAILY_) THEN
        ! In case of daily precipitation the minimum and maximum read for air_temp are used
        CALL snow_and_rain_from_precip(rain(:,:), snow(:,:), precipitation_data(:,:), &
          & 0.5_wp*(MAXVAL(forcing_input%air_temp%data_read(:,:,:),DIM=3) + MINVAL(forcing_input%air_temp%data_read(:,:,:),DIM=3)))
      ELSE
        ! else currentvalue
        CALL snow_and_rain_from_precip(rain(:,:), snow(:,:), precipitation_data(:,:), t_air(:,:))
      ENDIF
    ELSE
      CALL snow_and_rain_from_precip(rain(:,:), snow(:,:), precipitation_data(:,:), &
        & 0.5_wp*(tmax_data(:,:) + tmin_data(:,:)))
    ENDIF
    DEALLOCATE(precipitation_data)

    ! --- compute bottom pressure
    DO iblk = 1, nblks
      ! get elevation
      model => jsbach%models(model_id)%m
      tile => model%Get_top_tile()
      dsl4jsb_Get_memory(SRF)
      dsl4jsb_Get_var2D_onChunk(SRF, elevation, IN)

      press_srf(:,iblk) = bottom_pressure(elevation(:), t_air(:,iblk))
    ENDDO

    ! --- compute specific humidity
    ALLOCATE(sat_spec_humidity(ice, nblks))
    ALLOCATE(vapor_pressure_data(ice, nblks))
    SELECT CASE(forcing_options%type_of_qair_forcing)
    CASE(RH_)
      CALL get_data_and_update_field(nblks, ics, ice, model_id, forcing_input%rel_humidity, temporary_data1)
      ! check data
      hlp_r = MAXVAL(MAXVAL(temporary_data1,DIM=2))
      IF(hlp_r < 1._wp) THEN
        WRITE (message_text,*) 'Maximum value of relative humidity is ', hlp_r
        CALL message(TRIM(routine), message_text)
        CALL message(TRIM(routine),'relative humidity data probably represent fractions instead of percentages.')
      ENDIF
      ! Convert percentage data to fractional data:
      temporary_data1 = temporary_data1/100._wp
      ! Get specific humidity
      IF (forcing_input%rel_humidity%frequency .NE. TIMESTEP_) THEN
        IF (forcing_options%air_temperature_as_timestep) THEN
          DO iblk = 1, nblks
            CALL spec_humidity_from_rel_humidity(ice, sinlat(:,iblk), coslat(:,iblk), &
              & temporary_data1(:,iblk), MINVAL(forcing_input%air_temp%data_read(:,iblk,:),DIM=2), &
              & MAXVAL(forcing_input%air_temp%data_read(:,iblk,:),DIM=2), elevation(:), delta_time, q_air(:,iblk))
          ENDDO
        ELSE
          DO iblk = 1, nblks
            CALL spec_humidity_from_rel_humidity(ice, sinlat(:,iblk), coslat(:,iblk), &
              & temporary_data1(:,iblk), tmin_data(:,iblk), tmax_data(:,iblk), elevation(:), delta_time, q_air(:,iblk))
          ENDDO
        ENDIF
      ENDIF
      ! constrain daily specific humidity to be lower than or equal to the saturation spec. humidity
      DO iblk = 1, nblks
        CALL sat_specific_humidity(ice, t_air_K(:,iblk), press_srf(:,iblk), sat_spec_humidity(:,iblk))
      ENDDO
      IF (forcing_input%rel_humidity%frequency .EQ. TIMESTEP_) THEN
        q_air(:,:) = temporary_data1(:,:) * sat_spec_humidity(:,:)
      ELSE
        q_air(:,:) = MIN(q_air, sat_spec_humidity)
      ENDIF
      ! compute current vapour pressure
      vapor_pressure_data(:,:) = fun_vapour_pressure(q_air(:,:),press_srf(:,:))
    CASE(QAIR_)
      CALL get_data_and_update_field(nblks, ics, ice, model_id, forcing_input%qair, q_air(:,:))
      ! constrain daily specific humidity to be lower than or equal to the saturation spec. humidity
      DO iblk = 1, nblks
        CALL sat_specific_humidity(ice, t_air_K(:,iblk), press_srf(:,iblk), sat_spec_humidity(:,iblk))
      ENDDO
      q_air(:,:) = MIN(q_air(:,:),sat_spec_humidity(:,:))
      ! compute current vapour pressure
      vapor_pressure_data(:,:) = fun_vapour_pressure(q_air(:,:),press_srf(:,:))
    END SELECT
    DEALLOCATE(sat_spec_humidity)

    ! --- get cloud cover if required for short or longwave radiation
    IF ((forcing_options%type_of_shortwave_forcing .EQ. CLOUD_) &
        & .OR. (forcing_options%type_of_longwave_forcing .EQ. CLOUD_)) THEN

      ALLOCATE(cloud_cover_data(ice, nblks))

      IF(forcing_input%cloud_cover%frequency .EQ. TIMESTEP_) CALL finish(TRIM(routine), &
        & "Cloud cover not implemented as time-step forcing")
      CALL get_data_and_update_field(nblks, ics, ice, model_id, forcing_input%cloud_cover, cloud_cover_data)
      ! check data
      hlp_r = MAXVAL(MAXVAL(cloud_cover_data,DIM=2))
      IF(hlp_r < 1._wp) THEN
        WRITE (message_text,*) 'Maximum value of cloud cover is ', hlp_r
        CALL message(TRIM(routine), message_text)
        CALL message(TRIM(routine),'WARNING: cloud cover data probably represent fractions instead of percentages.')
      ENDIF
      ! Convert percentage data to fractional data:
      cloud_cover_data = cloud_cover_data/100._wp
    ENDIF

    ! --- compute shortwave radiation
    day_of_year=day_in_year(next_date)
    ALLOCATE(rad_uv_down(ice, nblks))
    SELECT CASE(forcing_options%type_of_shortwave_forcing)
    CASE(CLOUD_)   ! short wave radiation is computed from cloudcover
      ! Calculate shortwave
      CALL shortwave_from_cloudcover(day_of_year, cos_zenith_angle, press_srf, cloud_cover_data, &
        & rad_uv_down, swpar_srf_down, swnir_srf_down, help1, help2, help3)

    CASE(MEAN_RAD_)
      CALL get_data_and_update_field(nblks, ics, ice, model_id, forcing_input%shortwave, temporary_data1)
      IF(forcing_options%type_of_shortwave_scheme .NE. SW_SCHEME_DIRECT_) THEN
        CALL get_data_and_update_field(nblks, ics, ice, model_id, forcing_input%shortwave_pot, temporary_data2)
      ELSE
        temporary_data2 = 0._wp
      ENDIF
      CALL shortwave_from_shortwaveMean(day_of_year, cos_zenith_angle, coslat, sinlat, press_srf, &
        temporary_data1, temporary_data2, rad_uv_down, swpar_srf_down, swnir_srf_down, help1, help2, help3)

    CASE default
      CALL finish(TRIM(routine),'Please specifiy either cloud or mean_read for "FORCING_SW_TYPE" in forcing_ctl')
    END SELECT
    ! jsb3 standalone: swvis_srf_down = rad_uv_down + rad_par_down (& rad_par_down = swpar_srf_down)
    swvis_srf_down(:,:) = rad_uv_down(:,:) + swpar_srf_down(:,:)
    DEALLOCATE(rad_uv_down)

    ! --- compute longwave radiation
    SELECT CASE(forcing_options%type_of_longwave_forcing)
    CASE(CLOUD_)     ! longwave radiation is computed from cloudcover and vapour pressure
      lw_srf_down = longwave_rad_down(cloud_cover_data(:,:), vapor_pressure_data(:,:), t_air(:,:))
    CASE(MEAN_RAD_)
      IF (forcing_input%longwave%frequency .EQ. TIMESTEP_) THEN
        CALL get_data_and_update_field(nblks, ics, ice, model_id, forcing_input%longwave, lw_srf_down)
      ELSE
        ! jsbach3 driver: diurnal longwave downward radiation is generated such that the average radiation matches the
        ! observations but the diurnal cycle follows temperature to help calculation of a correct
        ! net radiation flux. Note that reduces biases resulting from assuming a constant flux, but
        ! this introduces uncertainty due to biases in the assumed dial course of temperature
        CALL get_data_and_update_field(nblks, ics, ice, model_id, forcing_input%longwave, temporary_data1)

        IF (forcing_options%air_temperature_as_timestep) THEN
          lw_srf_down = longwave_from_daily_longwave(temporary_data1(:,:), MINVAL(forcing_input%air_temp%data_read(:,:,:),DIM=3), &
            & MAXVAL(forcing_input%air_temp%data_read(:,:,:),DIM=3), t_air(:,:))
        ELSE
          lw_srf_down = longwave_from_daily_longwave(temporary_data1(:,:), tmin_data(:,:), tmax_data(:,:), t_air(:,:))
        ENDIF
      ENDIF
    CASE default
       CALL finish(TRIM(routine),'Please specifiy either cloud or mean_read for "FORCING_LW_TYPE" in forcing_ctl')
    END SELECT

    ! --- wind
    CALL get_data_and_update_field(nblks, ics, ice, model_id, forcing_input%wind_speed, wind_air)
    wind_10m = wind_air

    ! --- CO2 concentration
    CALL get_data_and_update_field(nblks, ics, ice, model_id, forcing_input%CO2_concentr, CO2_concentration)

    ! --- compute remaining required input for jsb4 interface
    CALL get_remaining_input_for_jsb4_interface(nblks, ice, t_air_K, press_srf, q_air, wind_air, &
      & t_srf_proc, fact_q_air_proc, fact_qsat_srf_proc, rough_h_srf_proc, &
      & drag_srf, t_acoef, t_bcoef, q_acoef, q_bcoef, pch)

    ! ---- clean up
    DEALLOCATE(t_air)
    IF(.NOT. forcing_options%air_temperature_as_timestep) THEN
      DEALLOCATE(tmin_data)
      DEALLOCATE(tmax_data)
    ENDIF

    IF ((forcing_options%type_of_shortwave_forcing .EQ. CLOUD_) &
        & .OR. (forcing_options%type_of_longwave_forcing .EQ. CLOUD_)) THEN
      DEALLOCATE(cloud_cover_data)
    ENDIF
    DEALLOCATE(vapor_pressure_data)

    DEALLOCATE(help1, help2, help3)
    DEALLOCATE(temporary_data1, temporary_data2)

  END SUBROUTINE get_interface_variables_from_external_forcing


  !>
  !! @brief init reading forcing data (read namelist, set options...)
  !!
  !! @par History
  !! JN: copied first version from jsbach3 mo_jsbalone_forcing
  !
  SUBROUTINE init_forcing()

    USE mo_jsb_namelist_iface, ONLY: POSITIONED, open_nml, position_nml

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':init_forcing'

    NAMELIST /forcing_ctl/          &
      & forcing_temp_file_prefix, forcing_temp_frequ, forcing_temp_const_tmin, forcing_temp_const_tmax, forcing_temp_ocean,   &
      & forcing_precip_file_prefix, forcing_precip_frequ, forcing_precip_const_precip,                                        &
        & forcing_precip_in_mm_per_day, forcing_precip_ocean,                                                                 &
      & forcing_sw_file_prefix, forcing_sw_frequ, forcing_sw_type, forcing_sw_scheme,                                         &
        & forcing_sw_const_cloud, forcing_sw_const_shortwave, forcing_sw_ocean,                                               &
      & forcing_table_sw_pot_file_prefix, forcing_table_sw_pot_frequ, forcing_table_sw_pot_const, forcing_table_sw_pot_ocean, &
      & forcing_lw_file_prefix, forcing_lw_frequ, forcing_lw_type,                                                            &
        & forcing_lw_const_cloud, forcing_lw_const_longwave, forcing_lw_ocean,                                                &
      & forcing_co2_file_prefix, forcing_co2_frequ, forcing_co2_const_co2, forcing_co2_unit, forcing_co2_ocean,               &
      & forcing_wind_file_prefix, forcing_wind_frequ, forcing_wind_const_wspeed, forcing_wind_ocean,                          &
      & forcing_qair_file_prefix, forcing_qair_frequ, forcing_qair_type, forcing_qair_const_rh, forcing_qair_ocean

    INTEGER:: nml_handler, read_status, f_unit

    !------------------------------------------------------------------------------------

    IF (my_process_is_stdio()) THEN
      ! ----- define default values for forcing_ctl namelist parameters
      forcing_temp_file_prefix = 'climate_'
      forcing_temp_frequ = 'DAILY'
      forcing_temp_const_tmin = 0.0_wp !jsb3: HUGE(0.0_wp)
      forcing_temp_const_tmax = 0.0_wp !jsb3: HUGE(0.0_wp)
      forcing_temp_ocean = 0._wp

      forcing_precip_file_prefix = 'climate_'
      forcing_precip_frequ = 'DAILY'
      forcing_precip_const_precip = 0.0_wp !jsb3: HUGE(0.0_wp)
      forcing_precip_in_mm_per_day = .TRUE.
      forcing_precip_ocean = 0._wp

      forcing_sw_file_prefix = 'climate_'
      forcing_sw_frequ = 'DAILY'
      forcing_sw_type = 'MEAN_RAD'
      forcing_sw_scheme = 'DIRECT'
      forcing_sw_const_cloud = 0.0_wp !jsb3: HUGE(0.0_wp)
      forcing_sw_const_shortwave = 0.0_wp !jsb3: HUGE(0.0_wp)
      forcing_sw_ocean = 0._wp
      forcing_table_sw_pot_file_prefix = 'climate_'
      forcing_table_sw_pot_frequ = 'DAILY'
      forcing_table_sw_pot_const = 0.0_wp !jsb3: HUGE(0.0_wp)
      forcing_table_sw_pot_ocean = 0._wp

      forcing_lw_file_prefix = 'climate_'
      forcing_lw_frequ = 'DAILY'
      forcing_lw_type = 'MEAN_RAD'
      forcing_lw_const_cloud = 0.0_wp !jsb3: HUGE(0.0_wp)
      forcing_lw_const_longwave = 0.0_wp !jsb3: HUGE(0.0_wp)
      forcing_lw_ocean = 100._wp

      forcing_co2_file_prefix = 'climate_'
      forcing_co2_frequ = 'DAILY'
      forcing_co2_const_co2 = 3.65e-4_wp
      forcing_co2_unit = ''
      forcing_co2_ocean = 0.00033114_wp

      forcing_wind_file_prefix = 'climate_'
      forcing_wind_frequ = 'DAILY'
      forcing_wind_const_wspeed = 0.0_wp !jsb3: HUGE(0.0_wp)
      forcing_wind_ocean = 0._wp

      forcing_qair_file_prefix = 'climate_'
      forcing_qair_frequ = 'DAILY'
      forcing_qair_type = 'QAIR'
      forcing_qair_const_rh = 100._wp
      forcing_qair_ocean = 1.0e-05

      ! ----- read the namelist
      nml_handler = open_nml ('namelist.jsbach')
      f_unit = position_nml ('forcing_ctl', nml_handler, status=read_status)
      IF (read_status .EQ. POSITIONED) READ(f_unit, forcing_ctl)

!TODO-JN: simplify and merge were possible
      ! ----- derive options
      forcing_options%air_temperature_as_timestep &
        & = (frequency_key_to_constant_val('temperature', forcing_temp_frequ) .EQ. TIMESTEP_)

      ! --- Precipitation
      forcing_options%precip_in_mm_per_day = forcing_precip_in_mm_per_day
      IF(forcing_options%precip_in_mm_per_day) THEN
        CALL message(TRIM(routine),'Precipitation forcing in mm/day')
      ELSE
        CALL message(TRIM(routine),'Precipitation forcing in kg/m^2/s')
      ENDIF
      ! --- QAIR
      SELECT CASE(TRIM(tolower(forcing_qair_type)))
      CASE("rh")
        forcing_options%type_of_qair_forcing = RH_
      CASE("qair")
        forcing_options%type_of_qair_forcing = QAIR_
      CASE default
        CALL finish(TRIM(routine),'Please specifiy either rh or qair for "FORCING_QAIR_TYPE" in forcing_ctl')
      END SELECT
      CALL message(TRIM(routine),"Type of air humidity forcing: "//TRIM(forcing_qair_type))
      ! --- SW
      SELECT CASE(TRIM(tolower(forcing_sw_type)))
      CASE("cloud")
        forcing_options%type_of_shortwave_forcing = CLOUD_
      CASE("mean_rad")
        forcing_options%type_of_shortwave_forcing = MEAN_RAD_
      CASE default
        CALL finish(TRIM(routine),'Please specifiy either cloud or mean_read for "FORCING_SW_TYPE" in forcing_ctl')
      END SELECT
      CALL message(TRIM(routine),"Type of shortwave radiation forcing: "//TRIM(forcing_sw_type))
      ! sw scheme
      SELECT CASE(TRIM(tolower(forcing_sw_scheme)))
      CASE("orig")
        forcing_options%type_of_shortwave_scheme= SW_SCHEME_ORIGINAL_
      CASE("ccdas")
        forcing_options%type_of_shortwave_scheme = SW_SCHEME_CCDAS_
      CASE("new")
        forcing_options%type_of_shortwave_scheme = SW_SCHEME_NEW_
      CASE("new_orig_cf")
        forcing_options%type_of_shortwave_scheme = SW_SCHEME_NEW_ORIG_CF_
      CASE("direct")
        forcing_options%type_of_shortwave_scheme = SW_SCHEME_DIRECT_
      CASE default
        CALL finish(TRIM(routine),'Please specifiy a value for "FORCING_SW_SCHEME" in forcing_ctl')
      END SELECT
      CALL message(TRIM(routine),"Radiation forcing scheme: "//TRIM(forcing_sw_scheme))
      !Assertion: sw type clouds only possible with the shortwave schemes original or ccdas
      IF(forcing_options%type_of_shortwave_forcing .EQ. CLOUD_) THEN
        IF(( forcing_options%type_of_shortwave_scheme .NE. SW_SCHEME_ORIGINAL_) &
            & .AND. ( forcing_options%type_of_shortwave_scheme .NE. SW_SCHEME_CCDAS_)) THEN
          CALL finish(TRIM(routine), &
            & "Calculating sw from cloud cover only implemented for sw scheme original or ccdas, please check selected sw scheme.")
        ENDIF
      ENDIF
      ! assert same frequency if not direct scheme
      IF (( forcing_options%type_of_shortwave_scheme .NE. SW_SCHEME_DIRECT_) &
          & .AND. (forcing_options%type_of_shortwave_forcing .NE. MEAN_RAD_ ) &
          & .AND. (forcing_sw_frequ .NE. forcing_table_sw_pot_frequ )) THEN
        CALL finish(TRIM(routine),"Invalid combination of shortwave ("//TRIM(forcing_sw_frequ)//") and mpot forcing frequency (" &
          & //TRIM(forcing_table_sw_pot_frequ)//"). [Need to be equal if not reading direct mean_rad forcing.]")
      ENDIF
      ! --- LW
      SELECT CASE(TRIM(tolower(forcing_lw_type)))
      CASE("cloud")
        forcing_options%type_of_longwave_forcing = CLOUD_
      CASE("mean_rad")
        forcing_options%type_of_longwave_forcing = MEAN_RAD_
      CASE default
        CALL finish(TRIM(routine),'Please specifiy either cloud or mean_read for "FORCING_LW_TYPE" in forcing_ctl')
      END SELECT
      CALL message(TRIM(routine),"Type of longwave radiation forcing: "//TRIM(forcing_lw_type))
      ! - CO2
      SELECT CASE(TRIM(tolower(forcing_co2_unit)))
      CASE("ppmv")
        CALL message(TRIM(routine),'CO2 forcing in ppmv')
        forcing_options%conv_CO2_2_MassRatio = molarMassCO2_kg / molarMassDryAir_kg * 0.000001_wp
      CASE("mol_per_mol")
        CALL message(TRIM(routine),'CO2 forcing in mol(CO2)/mol(Dry Air)')
        forcing_options%conv_CO2_2_MassRatio = molarMassCO2_kg / molarMassDryAir_kg
      CASE("kg_per_kg")
        CALL message(TRIM(routine),'CO2 forcing in kg(CO2)/kg(Dry Air)')
        forcing_options%conv_CO2_2_MassRatio = 1.0_wp
      CASE default
        CALL message(TRIM(routine),'WARNING: "FORCING_CO2_UNIT" missing in namelist forcing_ctl.')
        CALL message(TRIM(routine),'Assumption: unit is mol(CO2)/mol(DryAir). Verify your settings!')
        forcing_options%conv_CO2_2_MassRatio = molarMassCO2_kg / molarMassDryAir_kg
      END SELECT
      forcing_co2_const_co2 = forcing_co2_const_co2 * forcing_options%conv_CO2_2_MassRatio
    ENDIF

    ! ----- broadcast information
    IF (my_process_is_mpi_parallel()) THEN
      ! forcing options
      CALL p_bcast(forcing_options%precip_in_mm_per_day, p_io)
      CALL p_bcast(forcing_options%type_of_qair_forcing, p_io)
      CALL p_bcast(forcing_options%type_of_shortwave_forcing, p_io)
      CALL p_bcast(forcing_options%type_of_shortwave_scheme, p_io)
      CALL p_bcast(forcing_options%type_of_longwave_forcing, p_io)
      CALL p_bcast(forcing_options%conv_CO2_2_MassRatio, p_io)

      ! other namelist parameters
      CALL p_bcast(forcing_temp_file_prefix, p_io)
      CALL p_bcast(forcing_temp_frequ, p_io)
      CALL p_bcast(forcing_temp_const_tmin, p_io)
      CALL p_bcast(forcing_temp_const_tmax, p_io)
      CALL p_bcast(forcing_temp_ocean, p_io)

      CALL p_bcast(forcing_precip_file_prefix, p_io)
      CALL p_bcast(forcing_precip_frequ, p_io)
      CALL p_bcast(forcing_precip_const_precip, p_io)
      CALL p_bcast(forcing_precip_ocean, p_io)

      CALL p_bcast(forcing_sw_file_prefix, p_io)
      CALL p_bcast(forcing_sw_frequ, p_io)
      CALL p_bcast(forcing_sw_const_cloud, p_io)
      CALL p_bcast(forcing_sw_const_shortwave, p_io)
      CALL p_bcast(forcing_sw_ocean, p_io)
      CALL p_bcast(forcing_table_sw_pot_file_prefix, p_io)
      CALL p_bcast(forcing_table_sw_pot_frequ, p_io)
      CALL p_bcast(forcing_table_sw_pot_const, p_io)
      CALL p_bcast(forcing_table_sw_pot_ocean, p_io)

      CALL p_bcast(forcing_lw_file_prefix, p_io)
      CALL p_bcast(forcing_lw_frequ, p_io)
      CALL p_bcast(forcing_lw_const_cloud, p_io)
      CALL p_bcast(forcing_lw_const_longwave, p_io)
      CALL p_bcast(forcing_lw_ocean, p_io)

      CALL p_bcast(forcing_co2_file_prefix, p_io)
      CALL p_bcast(forcing_co2_frequ, p_io)
      CALL p_bcast(forcing_co2_const_co2, p_io)
      CALL p_bcast(forcing_co2_unit, p_io)
      CALL p_bcast(forcing_co2_ocean, p_io)

      CALL p_bcast(forcing_wind_file_prefix, p_io)
      CALL p_bcast(forcing_wind_frequ, p_io)
      CALL p_bcast(forcing_wind_const_wspeed, p_io)
      CALL p_bcast(forcing_wind_ocean, p_io)

      CALL p_bcast(forcing_qair_file_prefix, p_io)
      CALL p_bcast(forcing_qair_frequ, p_io)
      CALL p_bcast(forcing_qair_const_rh, p_io)
      CALL p_bcast(forcing_qair_ocean, p_io)
    ENDIF

    ! ----- init forcing_input
    ! --- temperature
    IF(tolower(forcing_temp_frequ) .NE. "timestep") THEN
      CALL init_forcing_input_variable(min_temp_varname, forcing_temp_frequ, forcing_temp_ocean, &
        & forcing_temp_file_prefix, exp_unit_temperature, forcing_temp_const_tmin, forcing_input%air_temp_min)
      CALL init_forcing_input_variable(max_temp_varname, forcing_temp_frequ, forcing_temp_ocean, &
        & forcing_temp_file_prefix, exp_unit_temperature, forcing_temp_const_tmax, forcing_input%air_temp_max)
    ELSE
      CALL init_forcing_input_variable(air_temp_varname, forcing_temp_frequ, forcing_temp_ocean, &
        & forcing_temp_file_prefix, exp_unit_temperature, forcing_temp_const_tmax, forcing_input%air_temp)
    ENDIF
    ! --- precipitation
    IF(forcing_options%precip_in_mm_per_day) THEN
      CALL init_forcing_input_variable(precipitation_varname, forcing_precip_frequ, forcing_precip_ocean, &
        & forcing_precip_file_prefix, unit_precip_mm_day, forcing_precip_const_precip, forcing_input%precipitation)
    ELSE
      CALL init_forcing_input_variable(precipitation_varname, forcing_precip_frequ, forcing_precip_ocean, &
        & forcing_precip_file_prefix, unit_precip_kg_m2_s, forcing_precip_const_precip, forcing_input%precipitation)
    ENDIF
    ! --- atmospheric humidity
    SELECT CASE(forcing_options%type_of_qair_forcing)
    CASE(RH_)
      CALL init_forcing_input_variable(rel_humidity_varname, forcing_qair_frequ, forcing_qair_ocean, &
        forcing_qair_file_prefix, unspecified_exp_unit ,forcing_qair_const_rh, forcing_input%rel_humidity)
    CASE(QAIR_)
      CALL init_forcing_input_variable(qair_varname, forcing_qair_frequ, forcing_qair_ocean, &
        forcing_qair_file_prefix, exp_unit_qair ,forcing_qair_const_rh, forcing_input%qair)
    CASE default
      CALL finish(TRIM(routine),'Please specifiy either qair or rh for "forcing_qair_type" in forcing_ctl')
    END SELECT
    ! --- shortwave radiation
    SELECT CASE(forcing_options%type_of_shortwave_forcing)
    CASE(CLOUD_)
      CALL init_forcing_input_variable(cloud_varname, forcing_sw_frequ, forcing_sw_ocean, &
        & forcing_sw_file_prefix, unspecified_exp_unit, forcing_sw_const_cloud, forcing_input%cloud_cover)
    CASE(MEAN_RAD_)
      CALL init_forcing_input_variable(shortwave_varname, forcing_sw_frequ, forcing_sw_ocean, &
        & forcing_sw_file_prefix, exp_unit_radiation, forcing_sw_const_shortwave, forcing_input%shortwave)
      ! jsbach3 text: "As a provisional solution we need here also a table of potential shortwave radiation"
      IF(forcing_options%type_of_shortwave_scheme .NE. SW_SCHEME_DIRECT_ ) THEN
        CALL init_forcing_input_variable(potential_shortwave_varname, forcing_table_sw_pot_frequ, forcing_table_sw_pot_ocean, &
          & forcing_table_sw_pot_file_prefix, exp_unit_radiation, forcing_table_sw_pot_const, forcing_input%shortwave_pot)
      ENDIF
    CASE default
      CALL finish(TRIM(routine),'Please specifiy either cloud or mean_read for "FORCING_SW_TYPE" in forcing_ctl')
    END SELECT
    ! --- longwave radiation
    SELECT CASE(forcing_options%type_of_longwave_forcing)
    CASE(CLOUD_)
      ! jsbach3 driver: Initialisation only if not already initialised for shortwave radiation
      IF (forcing_options%type_of_shortwave_forcing .NE. CLOUD_) THEN
        CALL init_forcing_input_variable(cloud_varname, forcing_lw_frequ, forcing_lw_ocean, &
          & forcing_lw_file_prefix, unspecified_exp_unit, forcing_lw_const_cloud, forcing_input%cloud_cover)
      ELSE !! Cloud cover is used for both shortwave and longwave
        ! In the jsbach3 driver the assumption in this case was accordingly that clouds were already read
        ! from the forcing sw file and that the constant is also the same!
        ! --> check compatibility
        IF (my_process_is_stdio()) THEN
          IF ((forcing_lw_frequ .NE. forcing_sw_frequ) .OR. (forcing_lw_const_cloud .NE. forcing_sw_const_cloud) &
              & .OR. (forcing_lw_file_prefix .NE. forcing_sw_file_prefix)) THEN
            CALL finish(TRIM(routine), &
              & "Inconsistent specifications for cloud data (FORCING_LW... and FORCING_SW..., please check)")
          ENDIF
        ENDIF
      ENDIF
    CASE(MEAN_RAD_)
      CALL init_forcing_input_variable(longwave_varname, forcing_lw_frequ, forcing_lw_ocean, &
        & forcing_lw_file_prefix, exp_unit_radiation, forcing_lw_const_longwave, forcing_input%longwave)
    CASE default
      CALL finish(TRIM(routine),'Please specifiy either cloud or mean_read for "FORCING_LW_TYPE" in forcing_ctl')
    END SELECT
    ! --- CO2
    CALL init_forcing_input_variable(CO2_varname, forcing_co2_frequ, forcing_co2_ocean, &
      & forcing_co2_file_prefix, unspecified_exp_unit, forcing_co2_const_co2, forcing_input%CO2_concentr)
    ! --- Wind
    CALL init_forcing_input_variable(wspeed_varname, forcing_wind_frequ, forcing_wind_ocean, &
      & forcing_wind_file_prefix, exp_unit_wspeed, forcing_wind_const_wspeed, forcing_input%wind_speed)

  END SUBROUTINE init_forcing

  !>
  !! @brief Provides forcing data for const, timestep or daily forcing (reads data from netcdf if required)
  !!
  !! @par History
  !! First version: JN (14.01.16)
  SUBROUTINE get_data_and_update_field(nblks, ics, ice, model_id, forcing_input_variable, current_forcing_array)

    INTEGER,                   INTENT(IN)    :: nblks, ics, ice, model_id
    TYPE(input_variable_type), INTENT(INOUT) :: forcing_input_variable
    REAL(wp),                  INTENT(OUT)   :: current_forcing_array(:,:)

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':get_data_and_update_field'

    LOGICAL ::  new_year_forcing     ! jsb3: defined comparing current_date with next_date, because with stepwise forcing
    LOGICAL ::  new_day_forcing      ! the dates in the forcing file shall correspond to the dates of stepwise output.

    INTEGER           :: year, month, day, hour, minute, second, ydate, time_steps_per_day, current_year, current_day
    INTEGER           :: frequency, currentTimeStep, thisDayStart, thisDayEnd, ith_timestep
    CHARACTER(len=8)  :: year_label
    CHARACTER(len=NF_MAX_NAME)  :: readUnit

    REAL(wp), POINTER :: return_pointer(:,:,:) !temporary pointer

    INTEGER :: nc_dim_id, nc_var_id, status, index_time
    TYPE(t_input_file) :: input_file

    !---------------------------
    frequency = forcing_input_variable%frequency

    IF(frequency .EQ. CONST_) THEN
      ! For constant forcing, only set correct constants and retun
      current_forcing_array = forcing_input_variable%constantValue
      CALL set_ocean_cells(nblks, ics, ice, model_id, forcing_input_variable%ocean_value, current_forcing_array)
      RETURN
    ENDIF

    !Assertion: Rest of the routine only implemented for daily or timestep frequency
    IF ( (frequency .NE. TIMESTEP_) .AND. (frequency .NE. DAILY_) ) CALL finish(TRIM(routine), &
      & 'Violation of assertion: routine only implemented for constant, daily or timestep frequency, please check program flow')
    !Assertion: delta_time needs to be smaller or equal 86400
    IF ((delta_time .GT. 86400._wp)) CALL finish(TRIM(routine), &
      & 'Violation of assertion: Reading forcing data not implemented for delta_time > day.')

    time_steps_per_day = 86400/INT(delta_time)
    CALL get_date_components(next_date, year, month, day, hour, minute, second)
    CALL get_date_components (current_date,  year=current_year, day=current_day)
    new_year_forcing = (year /= current_year)
    new_day_forcing = (day /= current_day)

    WRITE(year_label,'(i8.4)') year
    ydate = ISIGN(1,year)*(IABS(year)*10000+month*100+day)
    ith_timestep = IMerge_HMS2Sec(hour, minute, second)/int(delta_time)+1 !+1 because 00:00/X => 0; and 1 based vector
      !--> JN-TODO test in case of lresume

    ! --- if new current run or new file - get the timeseries of the forcing file and the current position in the timeseries
    IF (new_year_forcing .OR. lstart .OR. lresume) THEN
      IF(my_process_is_stdio()) THEN
        input_file = jsb_netcdf_open_input( &
          & trim(TRIM(ADJUSTL(forcing_input_variable%file_name_prefix)) // TRIM(ADJUSTL(year_label)) // '.nc') , model_id)

        ! TODO: via mo_jsb_io_netcdf_iface / mo_jsb_io_netcdf?
        ! get length and content of the time dimenson/var
        CALL nf_check(nf_inq_dimid(input_file%file_id,"time",nc_dim_id))
        CALL nf_check(nf_inq_dimlen(input_file%file_id, nc_dim_id, forcing_input_variable%length_time_series))
        IF(ASSOCIATED(forcing_input_variable%timevalues)) DEALLOCATE(forcing_input_variable%timevalues)
        ALLOCATE(forcing_input_variable%timevalues(forcing_input_variable%length_time_series))
        CALL nf_check(nf_inq_varid(input_file%file_id, 'time', nc_var_id))
        CALL nf_check(nf_get_vara_double(input_file%file_id, nc_var_id, (/ 1 /), &
          & (/ forcing_input_variable%length_time_series /), forcing_input_variable%timevalues))

        ! Assert: test that the netcdf file unit for this variable is the expected unit
        IF (TRIM(forcing_input_variable%expectedUnit) .NE. TRIM(unspecified_exp_unit)) THEN
          CALL nf_check(nf_inq_varid(input_file%file_id, TRIM(forcing_input_variable%variable_name), nc_var_id))
          status = nf_get_att_text(input_file%file_id, nc_var_id, name_of_unit_attribute, readUnit)
          IF (status /= NF_NOERR) THEN
            CALL message(TRIM(routine), &
                 & 'WARNING: netcdf error when accessing units attribute of variable '//TRIM(forcing_input_variable%variable_name))
          ELSE
            IF(TRIM(forcing_input_variable%expectedUnit) .NE. TRIM(readUnit)) CALL finish(TRIM(routine), &
              & 'Forcing file does not contain the expected unit "' //TRIM(forcing_input_variable%expectedUnit) &
              & //'" but: "'//TRIM(readUnit)//'" for the variable '//TRIM(forcing_input_variable%variable_name))
          ENDIF
        ENDIF

        ! TODO: further assertions? e.g. that lon lat are correct?

        currentTimeStep = 1
        thisDayStart = 0
        thisDayEnd = 0
        ! determine current date in the file
        DO WHILE (currentTimeStep .LE. forcing_input_variable%length_time_series)
          IF ( ydate .EQ. INT(forcing_input_variable%timevalues(currentTimeStep)) ) THEN
            IF (thisDayStart .EQ. 0) thisDayStart = currentTimeStep
            thisDayEnd = currentTimeStep
          ELSEIF ( ydate .LT. INT(forcing_input_variable%timevalues(currentTimeStep)) ) THEN
            EXIT
          ENDIF
          currentTimeStep = currentTimeStep + 1
        ENDDO

        CALL jsb_netcdf_close_file(input_file)
        forcing_input_variable%currentDayStart = thisDayStart
        forcing_input_variable%currentDayEnd = thisDayEnd

      ENDIF
        CALL p_bcast(forcing_input_variable%currentDayStart, p_io)
        CALL p_bcast(forcing_input_variable%currentDayEnd, p_io)

      IF(ASSOCIATED(forcing_input_variable%data_read)) DEALLOCATE(forcing_input_variable%data_read)
      ALLOCATE(forcing_input_variable%data_read(ice, nblks, time_steps_per_day))

    ENDIF

    thisDayStart = forcing_input_variable%currentDayStart
    thisDayEnd = forcing_input_variable%currentDayEnd

    ! --- if its a new day or new run - read the data
    IF (new_day_forcing .OR. lstart .OR. lresume) THEN
      ! check timeseries
      IF(my_process_is_stdio()) THEN
        !Assertion: currentTimeStep needs to be contained in the file
        IF(thisDayStart > forcing_input_variable%length_time_series) CALL finish(TRIM(routine), &
          & 'Day not found: ' //TRIM(int2string(ydate)))

        !Assertion: check for required number of timesteps
        IF( frequency .EQ. TIMESTEP_ ) THEN
          IF(thisDayEnd-thisDayStart+1 .NE. time_steps_per_day) CALL finish(TRIM(routine), &
             & 'Forcing file does not contain the demanded number of time steps per day: ' //TRIM(int2string(time_steps_per_day)) &
             & //' -- but: '//TRIM(int2string(thisDayEnd-thisDayStart+1)))
        ELSEIF( frequency .EQ. DAILY_ ) THEN
           IF(thisDayStart .NE. thisDayEnd) CALL finish(TRIM(routine), &
             & 'Forcing file contains more than one time step per day: ' //TRIM(int2string(thisDayEnd-thisDayStart+1)) &
             & //' -- daily reading only implemented for files with one time step per day.')
        ENDIF
      ENDIF
      ! get data
      input_file = jsb_netcdf_open_input( &
        & trim(TRIM(ADJUSTL(forcing_input_variable%file_name_prefix)) // TRIM(ADJUSTL(year_label)) // '.nc') , model_id)
      return_pointer => jsb_netcdf_read_2d_time(input_file, forcing_input_variable%variable_name, &
        & start_time_step=thisDayStart, end_time_step=thisDayEnd, fill_array=forcing_input_variable%data_read)
      CALL jsb_netcdf_close_file(input_file)

      ! Set the ocean cells of the read data to the constant ocean value
      IF(forcing_input_variable%frequency .EQ. TIMESTEP_) THEN
        DO index_time = 1,time_steps_per_day
          CALL set_ocean_cells(nblks, ics, ice, model_id, &
            & forcing_input_variable%ocean_value, forcing_input_variable%data_read(:,:,index_time))
        ENDDO
        ! Increase the start and end for the next day to be read
        forcing_input_variable%currentDayStart=forcing_input_variable%currentDayEnd+1
        forcing_input_variable%currentDayEnd=forcing_input_variable%currentDayEnd+time_steps_per_day
      ELSEIF(forcing_input_variable%frequency .EQ. DAILY_) THEN
        CALL set_ocean_cells(nblks, ics, ice, model_id, forcing_input_variable%ocean_value, forcing_input_variable%data_read(:,:,1))
        ! Increase the start and end for the next day to be read
        forcing_input_variable%currentDayStart=forcing_input_variable%currentDayStart+1
        forcing_input_variable%currentDayEnd=forcing_input_variable%currentDayEnd+1
      ENDIF
    ENDIF

    IF(forcing_input_variable%frequency .EQ. TIMESTEP_) THEN
      current_forcing_array = forcing_input_variable%data_read(:,:,ith_timestep)
    ELSEIF(forcing_input_variable%frequency .EQ. DAILY_) THEN
      current_forcing_array = forcing_input_variable%data_read(:,:,1)
    ENDIF

  END SUBROUTINE get_data_and_update_field

  !>
  !! @brief sets given values for input_variable_type field in forcing_input
  !!
  !! @par History
  !! JN: based on jsbach3 mo_jsbalone_forcing init_variable_access
  SUBROUTINE init_forcing_input_variable(variable_name, forcing_frequ_keyword, forcing_ocean_value, &
    & forcing_file_prefix, expected_unit, constantValue, forcing_input_variable)

    CHARACTER(len=*),INTENT(IN) :: variable_name
    CHARACTER(len=*),INTENT(IN) :: forcing_frequ_keyword
    REAL(wp),        INTENT(IN) :: forcing_ocean_value
    CHARACTER(len=*),INTENT(IN) :: forcing_file_prefix
    CHARACTER(len=*),INTENT(IN) :: expected_unit
    REAL(wp),        INTENT(IN) :: constantValue
    TYPE(input_variable_type),INTENT(INOUT) :: forcing_input_variable ! structure filled in this subroutine

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':init_forcing_input_variable'

    forcing_input_variable%frequency=frequency_key_to_constant_val(variable_name, forcing_frequ_keyword)
    forcing_input_variable%variable_name = variable_name
    forcing_input_variable%file_name_prefix = forcing_file_prefix
    forcing_input_variable%ocean_value = forcing_ocean_value
    forcing_input_variable%expectedUnit = expected_unit
    forcing_input_variable%constantValue = constantValue
    NULLIFY(forcing_input_variable%timevalues)

  END SUBROUTINE init_forcing_input_variable

  !>
  !! @brief assignes forcing frequency keyword to integer constant
  !!
  !! @par History
  !! JN (2016-01-14): based on jsbach3 mo_jsbalone_forcing init_variable_access
  INTEGER FUNCTION frequency_key_to_constant_val(variable_name, forcing_frequ_keyword)

    CHARACTER(len=*),INTENT(IN) :: variable_name
    CHARACTER(len=*),INTENT(IN) :: forcing_frequ_keyword

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':frequency_key_to_constant_val'

    SELECT CASE(TRIM(tolower(forcing_frequ_keyword)))
    CASE("monthly")
      CALL finish(TRIM(routine),'Reading files with monthly frequency is currently not supported, please implement.')
    CASE("daily")
      frequency_key_to_constant_val = DAILY_
    CASE("const")
      frequency_key_to_constant_val = CONST_
    CASE("timestep")
      frequency_key_to_constant_val = TIMESTEP_
    CASE default
      CALL finish(TRIM(routine), 'Please specify a valid reading frequency in your namelist for the variable '//TRIM(variable_name))
    END SELECT
    CALL message(TRIM(routine), 'Frequency for forcing variable "'//TRIM(variable_name)//'" is: '//TRIM(forcing_frequ_keyword))
  END FUNCTION frequency_key_to_constant_val


  !>
  !! @brief calculates the remaining drivers required for jsbach4, particularly drag coefficient
  !!
  !! @par History
  !! JN: copied required parts from jsbach3 mo_jsbalone and adapted for jsb4
  SUBROUTINE get_remaining_input_for_jsb4_interface(nblks, ice, temp_air, pressure, qair, wind, &
      & t_srf_proc, fact_q_air_proc, fact_qsat_srf_proc, rough_h_srf_proc, &
      & drag_srf, t_acoef, t_bcoef, q_acoef, q_bcoef, pch)

    INTEGER,  INTENT(IN)  :: nblks, ice
    REAL(wp), INTENT(IN)  :: temp_air(:,:)
    REAL(wp), INTENT(IN)  :: pressure(:,:)
    REAL(wp), INTENT(IN)  :: qair(:,:)
    REAL(wp), INTENT(IN)  :: wind(:,:)
    REAL(wp), INTENT(IN)  :: t_srf_proc(:,:)
    REAL(wp), INTENT(IN)  :: fact_q_air_proc(:,:)
    REAL(wp), INTENT(IN)  :: fact_qsat_srf_proc(:,:)
    REAL(wp), INTENT(IN)  :: rough_h_srf_proc(:,:)

    REAL(wp), INTENT(OUT) :: drag_srf(ice,nblks)
    REAL(wp), INTENT(OUT) :: t_acoef(ice,nblks)
    REAL(wp), INTENT(OUT) :: t_bcoef(ice,nblks)
    REAL(wp), INTENT(OUT) :: q_acoef(ice,nblks)
    REAL(wp), INTENT(OUT) :: q_bcoef(ice,nblks)
    REAL(wp), INTENT(OUT) :: pch(ice,nblks)

    INTEGER  :: iblk
    REAL(wp) :: zcons8, zcons9, zcons11, zcons12, zkappa
    REAL(wp), DIMENSION(ice):: zdu2, ztvd, qsat_surf, ztvs, zg, zril, zcons, zchnl, zcfnchl, zcfhl, zscfl, zucfhl

    ! in jsb3: rd = GasConstantDryAir; cpd = SpecificHeatDryAirConstPressure; grav = Gravity
    zcons8 = 2._wp * cb
    zcons9 = 3._wp * cb
    zcons11 = 3._wp * cb * cc
    zcons12 = cvdifts * time_step_len * grav / rd
    zkappa = rd / cpd

    DO iblk = 1,nblks

      !------------------------------------------------------------------------------------
      ! Approximation of cdrag
      !------------------------------------------------------------------------------------
      ! squared wind shear ! minimum wind speed square from echam/mo_surface_land.f90:precalc_land
      zdu2(:) = MAX(wind(:,iblk)**2, 1._wp)

      ! virtual potential air temperature (see mo_surface_boundary.f90)
      ! according to Saucier, WJ Principles of Meteoroligical Analyses
      ! tv = t * (1 + 0.61 * q * t)    ! virtual temperature
      ! td = t * ( 1000 / p_mb ) ^ R/cdp  ! potential temperature
      ! tvd = tair * (100000/p_pa)^zkappa * 1 + vtmpc1 * qair) ! virtual potential temperature
      ztvd(:) = ( temp_air(:,iblk) * ( 100000._wp/pressure(:,iblk))**zkappa ) * ( 1._wp + vtmpc1 * qair(:, iblk))

      ! virtual potential surface temperature
      CALL sat_specific_humidity(ice, t_srf_proc(:,iblk),pressure(:,iblk), qsat_surf(:))
      ztvs(:) = t_srf_proc(:,iblk) * ( 100000._wp/pressure(:,iblk))**zkappa * &
        & ( 1._wp + vtmpc1 * ( fact_qsat_srf_proc(:,iblk) * qsat_surf(:) + ( 1._wp - fact_q_air_proc(:,iblk) ) * qair(:,iblk)))

      ! geopotential of the surface layer (see echam's auxhybc.f90 & geopot.f90)
      ! If HeightWind is set, then the measurement height + an offset from the average vegetation height
      ! is used. Otherwise, the code defaults to the half-level of ECHAM's lowest atmospheric layer
      IF(HeightWind > 0._wp) THEN
        zg(:) = (0.75_wp * Veg_height + HeightWind) * grav
      ELSE
        zg(:) = ztvd(:) *  rd * LOG(1._wp / ( 1._wp - 0.0025_wp ))
      ENDIF

      ! Richardson number (dry, Brinkop & Roeckner 1995, Tellus)
      ! ztvd, ztvs are now virtual potential temperatures, changed by Thomas Raddatz 07.2014
      zril(:) = zg(:) * ( ztvd(:) - ztvs(:) ) / ( zdu2(:) * ztvd(:) )
      zril(:) = MAX(MIN(zril(:), 5._wp), -5._wp)

      ! Neutral drag coefficient for momentum and heat
      zchnl(:) = (ckap / LOG(1._wp + zg(:) / (grav * rough_h_srf_proc(:,iblk) )))**2

      ! account for stable/unstable case: helper variables
      zscfl(:) = SQRT ( 1._wp + ABS(zril(:)))
      zucfhl(:) = 1._wp / (1._wp + zcons11 * zchnl(:) * SQRT(ABS(zril(:)) * (1._wp &
        & + zg(:) / (grav * rough_h_srf_proc(:,iblk)))))

      ! ignoring cloud water correction (see mo_surface_land.f90)
      zcons(:) = zcons12 * pressure(:,iblk) / ( temp_air(:,iblk) * (1._wp + vtmpc1 * qair(:,iblk)))
      zcfnchl(:)  = zcons(:) * SQRT(zdu2(:)) * zchnl(:)

      ! Stable / Unstable case
      WHERE ( zril(:) .GT. 0._wp )
        zcfhl(:) = zcfnchl(:) / (1._wp + zcons8 * zril(:) * zscfl(:))
      ELSEWHERE
        zcfhl(:) = zcfnchl(:) * (1._wp - zcons9 * zril(:) * zucfhl(:))
      END WHERE

      drag_srf(:,iblk) = zcfhl(:)
      pch(:,iblk) = zcfhl(:) / zcfnchl(:) * zchnl(:)

    ENDDO

    !---------------------------------------------------------------------------------------------------------------
    ! Computation of Richtmeyr-Morton Coefficients
    ! This follows now Jan Polcher's explicit solution, i.e. atmospheric conditions at t+1 are assumed to be valid
    !---------------------------------------------------------------------------------------------------------------
    t_acoef = 0.0_wp
    t_bcoef = cpd * temp_air + HeightHumidity * grav
    q_acoef = 0.0_wp
    q_bcoef = qair

  END SUBROUTINE get_remaining_input_for_jsb4_interface

  !>
  !! @brief converts daily relative humidity to specific humidity considering the daily course of temperature
  !!
  !! @par History
  !! JN: copied from jsbach3 mo_jsbalone_forcing and adapted for jsb4 input
  SUBROUTINE spec_humidity_from_rel_humidity(ice, sinlat, coslat, rel_humidity, Tmin, Tmax, elevation, delta_time, qair)
    ! This subroutine converts daily relative humidity to specific humidity considering the daily course of temperature
    ! ops! Using this routine only makes sense together with the generated diurnal cycle of temperature as in
    ! SUBROUTINE instantly_from_daily_temp_2
    ! Otherwise the results will be faulty!

    !mo_orbit with the required routine seems to be available in ICON and ECHAM
    USE mo_orbit,      ONLY: inquire_declination

    IMPLICIT NONE

    INTEGER, INTENT(IN)  :: ice
    REAL(wp),INTENT(IN)  :: sinlat(1:ice)
    REAL(wp),INTENT(IN)  :: coslat(1:ice)
    REAL(wp),INTENT(IN)  :: rel_humidity(1:ice)  !! mean daily relative humidity
    REAL(wp),INTENT(IN)  :: Tmin(1:ice)          !! minimum temperature (= temperature at sunrise)
    REAL(wp),INTENT(IN)  :: Tmax(1:ice)          !! maximum temperature (= at 2pm)
    REAL(wp),INTENT(IN)  :: elevation(1:ice)     !! height above sea leavel
    REAL(wp),INTENT(IN)  :: delta_time           !! time step of model
    REAL(wp),INTENT(OUT) :: qair(1:ice)          !! the specific humidity of the day

    ! constants
    REAL(wp), PARAMETER :: pio180 = pi/180._wp     ! conversion factor from degrees to radians ("pi over 180")
    REAL(wp), PARAMETER :: rat1=14._wp/24._wp      ! day fraction at which maximum temperature shall occur
    REAL(wp), PARAMETER :: rat2=4._wp/24._wp       ! minimum length of normal day (as fraction of day)
    REAL(wp), PARAMETER :: rat3=20._wp/24._wp      ! maximum length of normal day (as fraction of day)
    REAL(wp), PARAMETER :: rat4=2._wp/24._wp       ! day fraction of dawn times (before sunrise + after sunset), used to set ...
                                                   ! ... temperature of a normal day to T_mean at sunrise-rat4/2 and sunset+rat4/2
    ! Temperature calculations
    REAL(wp) :: declination       ! solar declination
    REAL(wp) :: cpds, spds        ! cosine and sine contributions to zenith angle
    REAL(wp) :: frac_daylight     ! fraction of day between sunrise and sunset
    REAL(wp) :: frac_since_sunset ! time between current time step and sunset expressed as fraction of day
    REAL(wp) :: frac_sunrise      ! time of sunrise expressed as fraction of day since midnight
    REAL(wp) :: frac_sunset       ! time of sunset expressed as fraction of day since midnight
    REAL(wp) :: temp_sunset       ! temperature at sunset
    REAL(wp) :: dtmp              ! mean day temperature
    REAL(wp) :: dtran             ! day temperature range: maximum - minimum
    REAL(wp) :: sd                ! modulation factor for temperature
    REAL(wp) :: arg,nstep
    REAL(wp),ALLOCATABLE :: tmp(:)  ! diurnal temperature

    ! Qair calculations
    REAL(wp) :: rh_test,rh_last,qsat_last ! temporary RH estimates
    REAL(wp),ALLOCATABLE :: qsat(:)      ! diurnal saturated specific humidity
    REAL(wp),ALLOCATABLE :: inv_qsat(:)  ! inverse diurnal saturated specific humidity
    REAL(wp),ALLOCATABLE :: qair_test(:) ! diurnal actual specific humidity, temporary
    REAL(wp),ALLOCATABLE :: pressure(:)  ! diurnal air pressure
    REAL(wp) :: tstep        ! fractional time step of the day

    ! helpers
    INTEGER    :: i,j,jj,time_steps_per_day
    INTEGER    :: offset(1)

    CALL inquire_declination(declination)

    nstep = (86400.0_wp/delta_time)
    time_steps_per_day = 86400/int(delta_time)

    IF(.NOT.(ALLOCATED(tmp))) ALLOCATE(tmp(time_steps_per_day))
    IF(.NOT.(ALLOCATED(pressure))) ALLOCATE(pressure(time_steps_per_day))
    IF(.NOT.(ALLOCATED(qsat))) ALLOCATE(qsat(time_steps_per_day))
    IF(.NOT.(ALLOCATED(inv_qsat))) ALLOCATE(inv_qsat(time_steps_per_day))
    IF(.NOT.(ALLOCATED(qair_test))) ALLOCATE(qair_test(time_steps_per_day))

    DO i = 1,ice
      spds = sinlat(i)*SIN(declination)  ! sine contribution to solar zenith angle
      cpds = coslat(i)*COS(declination)  ! cosine contribution to solar zenith angle

      arg = -spds/cpds                        ! faulty?: for cpds=0 arg gets infinite
      IF (arg > 1._wp) THEN               ! polar night:
        frac_daylight = 0._wp            ! ... day length is zero
      ELSE IF (arg < -1._wp) THEN         ! polar day:
        frac_daylight = 1._wp            ! ... day length is the whole day
      ELSE                                ! normal day / night:
        frac_daylight = ACOS(arg)/pi    ! ... day length expressed as fraction of the whole day
      END IF
      dtmp = 0.5_wp*(Tmax(i)+Tmin(i)) ! The mean day temperature
      dtran = Tmax(i)- Tmin(i)        ! temperature range
      DO j = 1,time_steps_per_day
        tstep=j/nstep
        ! -------- normal day ---------------------------------------------
        IF (frac_daylight>=rat2 .AND. frac_daylight<=rat3) THEN
          frac_sunrise = 0.5_wp - frac_daylight/2._wp             ! time of sunrise expressed as fraction of day since midnight
          frac_sunset = frac_sunrise+frac_daylight                ! time of sunset expressed as fraction of day since midnight

          ! For time steps at daylight ...
          IF (tstep>frac_sunrise .AND. tstep<frac_sunset) THEN
            sd = COS(pi*(tstep-rat1)/(frac_daylight/2._wp+rat4)) ! modulation factor for day-temperature such that
                                                                 !   temperature is at maximum at rat1 (2 pm) and at average at
                                                                 !   sunrise-rat4/2 and sunset+rat4/2
            tmp(j) = dtmp + sd*dtran/2._wp                       ! temperature modulated by temperature range

          ! For time steps at night ...
          ELSE
            sd = COS(pi*(frac_sunset-rat1)/(frac_daylight/2._wp+rat4)) ! modulation factor for temperature at sunset
            temp_sunset = dtmp + dtran/2._wp*sd                        ! temperature at sunset
            frac_since_sunset = MOD(tstep-frac_sunset+1._wp,1._wp)     ! day fraction since sunset
            tmp(j) = Tmin(i)+(temp_sunset-Tmin(i))*(1._wp-frac_since_sunset/(1._wp-frac_daylight))! linear interpolation
                                                                     ! ... of temperature such that temperature is equal to the
                                                                     ! ... temperature at sunset and drops until the minimum ...
                                                                     ! ... temperature is met at sunrise.
          END IF
        ! --- long day (close to polar day) ------------------------------------------------------
        ELSEIF (frac_daylight>rat3) THEN
          sd = COS(pi*(tstep-rat1)/(frac_daylight/2._wp+rat4)) ! modulation factor for day-temperature such that ...
                                                               ! ... temperature is at maximum at rat1 (2 pm) and at average
                                                               ! ... at sunrise-rat4/2 and sunset+rat4/2
          tmp(j) = dtmp + dtran/2._wp*sd                       ! temperature modulated by temperature range
        ! ------------------------------- short day (close to polar night) ---------------------------------------------------
        ELSE
          tmp(j) = dtmp                                        ! temperature is set to mean temperture at all time steps
        END IF
      ENDDO

      ! from function bottom pressure; in jsb3: p0sl_bg = p_sealevel; grav = Gravity; rd = GasConstantDryAir
      pressure(:) = &
        & p0sl_bg * ( 1._wp - elevation(i)*gamma /(tmelt+tmp(:) + elevation(i)*gamma))**(grav/rd/gamma)

      CALL sat_specific_humidity(time_steps_per_day, tmp(:) + tmelt, pressure(:), qsat(:)) ! present saturated qair

      ! implied RH from minimum saturated specific humidity
      rh_test=SUM(MINVAL(qsat(:))/qsat(:))/nstep
      IF( rel_humidity(i) <= rh_test ) THEN
        ! Simply case, relative humitidy that assuming a constant qair there is never
        ! saturation during the day. In this case, qair can be inferred by inverting
        ! RH = 1/n * sum ( QAIR / QSAT ) = QAIR / n * sum ( 1 / QSAT )
        ! -> QAIR = n * RH / sum ( 1 / QSAT )
        inv_qsat(:) = 1._wp / qsat(:)
        qair(i) = nstep * rel_humidity(i) / SUM( inv_qsat(:) )
      ELSE
        ! The case where RH implies saturation at night. Under these conditions
        ! search for the qair value which is closed to giving the prescribed RH
        ! taking account of the implied period during the day where saturation occurs.
        ! Surely there's a better way to do this, but I cannot be bothered now.
        offset(:)=MAXLOC(qsat(:))
        IF( rel_humidity(i) >= 1._wp ) THEN
           qair(i)=qsat(offset(1))
        ELSE
          j=0
          jj=offset(1)
          rh_test=1._wp
          DO WHILE ( rh_test > rel_humidity(i) .AND. j < time_steps_per_day )
            ! save value from last time
            qsat_last=qsat(jj)
            rh_last=rh_test

            ! find tome step after hottest time of the day
            IF(j+offset(1)>time_steps_per_day) THEN
              jj=j+offset(1)-time_steps_per_day
            ELSE
              jj=j+offset(1)
            ENDIF

            ! create qair values constrained by qsat at low temperatures
            WHERE(qsat(:) <= qsat(jj))
              qair_test(:)=qsat(:)
            ELSEWHERE
              qair_test(:)=qsat(jj)
            END WHERE

            ! test whether target RH is reached
            rh_test = SUM(qair_test(:) / qsat(:)) /nstep
            j=j+1
          ENDDO
          ! The above solution is inprecise. Attempt to reconcile this by linear approximation
          qair(i)=qsat(jj)+(qsat_last-qsat(jj))*(rel_humidity(i)-rh_test)/(rh_last-rh_test)
        ENDIF
      ENDIF
    ENDDO

    IF((ALLOCATED(tmp))) DEALLOCATE(tmp)
    IF((ALLOCATED(pressure))) DEALLOCATE(pressure)
    IF((ALLOCATED(qsat))) DEALLOCATE(qsat)
    IF((ALLOCATED(inv_qsat))) DEALLOCATE(inv_qsat)
    IF((ALLOCATED(qair_test))) DEALLOCATE(qair_test)
  END SUBROUTINE spec_humidity_from_rel_humidity

  !>
  !! @brief Returns saturation specific humidity for given temperature and pressure
  !!
  !! @par History
  !! JN: copied from jsbach3 mo_atmosphere
  SUBROUTINE sat_specific_humidity(temp_dim, temp, pressure, qsat)
    ! Returns saturation specific humidity for given temperature and pressure (at some atmospheric level or at surface)
    ! Uses Eq. 2.27 of Fundamentals of Atmospheric Modeling for the saturated case, but the saturation vapor pressure
    ! of water over a liquid surface resp. ice (see pp 32-34 of Ref.) is computed as in ECHAM5 (mo_convect_tables)

    INTEGER,  INTENT(IN) :: temp_dim
    REAL(wp), INTENT(IN) :: temp(:)      ! Air temperature at level [K]
    REAL(wp), INTENT(IN) :: pressure(:)  ! Pressure at level [Pa]
    REAL(wp), INTENT(OUT):: qsat(:)

integer :: i
    REAL(wp) :: vtmpc1
    INTEGER  :: it(1:temp_dim)
    REAL(wp) :: tluc(1:temp_dim)

    vtmpc1 = 0.6077686814143877_wp
    it = NINT(temp*1000._wp)

    !Assertion: check tlucua bounds
    DO i = 1,temp_dim
      IF ((it(i) .LT. LBOUND(tlucua,dim=1)) .OR. (it(i) .GT. UBOUND(tlucua,dim=1))) THEN
        WRITE (message_text,*) &
          & 'Index out of tlucua bounds (would give a seg. fault!) - index was: ', it(i), ' i.e. temp was: ', temp(i)
        CALL finish('sat_specific_humidity', message_text)
      ENDIF
    ENDDO

    tluc = tlucua(it)

    WHERE (it >= jptlucu1 .AND. it <= jptlucu2)
      qsat = tluc / (pressure - vtmpc1*tluc)
    ELSEWHERE
      qsat = HUGE(1._wp)
    END WHERE

  END SUBROUTINE sat_specific_humidity

  !>
  !! @brief sets all not land (i.e. not notsea ;)) cells to the given constant
  !!
  !! @par History
  !! First version: JN (20.01.16)
  SUBROUTINE set_ocean_cells(nblks, ics, ice, model_id, constant, array)

    USE mo_jsb_tile_class,     ONLY: t_jsb_tile
    USE mo_jsb_model_class,    ONLY: t_jsb_model

    dsl4jsb_Def_memory(SRF)
    dsl4jsb_Real2D_onChunk :: land

    INTEGER, INTENT(IN) :: nblks, ics, ice, model_id
    REAL(wp), INTENT(IN) :: constant
    REAL(wp), INTENT(INOUT) :: array(:,:)
    INTEGER :: iblk

    TYPE(t_jsb_model),  POINTER  :: model
    TYPE(t_jsb_tile), POINTER :: tile

    DO iblk = 1, nblks
      model => jsbach%models(model_id)%m
      tile => model%Get_top_tile()

      dsl4jsb_Get_memory(SRF)
      dsl4jsb_Get_var2D_onChunk(SRF, land,    IN)

      WHERE (land .EQ. 0._wp)
        array(:,iblk) = constant
      END WHERE
    ENDDO
  END SUBROUTINE set_ocean_cells

  !>
  !! @brief estimate timestep temperature from daily min- and max temperature
  !!
  !! @par History
  !! JN: copied from jsbach3 mo_jsbalone_forcing and adapted for 2D arrays
  !
  ! --- instantly_from_daily_temp_2() --------------------------------------------------------------------------------------
  !
  ! The same as subroutine instantly_from_daily_temp_1(), but with temperature minimum and maximum as input instead of
  ! mean temperature and temperature range.
  ! In contrast to Marko's subroutine daytemp() (from which this code is derived) the routine
  ! correctly accounts for local time, i.e. temperatures are correct at all longitudes!
  !
  ! More precisely the temperature development is modelled as follows:
  ! Let: T_min = T_mean - 0.5*T_range,  and: T_max = T_mean + 0.5*T_range. Then for a
  !     (i) normal day (more than 4, but less than 20 hours of daylight):
  !         --> during daylight: temperature behaves like a positive sine-wave, such that T_max is reached at 14:00 and
  !             T_min at sunrise.
  !         --> during night: starting at sunset the temperature drops linearly until T_min at sunrise.
  !    (ii) short day (less than 4 hours of daylight):
  !         temperature is set to T_mean at all time steps, i.e. temperature is assumed to be constant.
  !   (iii) long day (more than 20 hours of daylight):
  !         temperature is modelled as in the daylight-case of (i), but for the whole day, even during night hours.
  !
  ! REMARKS:
  !      -- The time computations are based on the Julian/Gregorian calender. Hence there is a slight mismatch between the
  !         correct begin of a year, as defined by the orbit around the sun, and the begin of the year according to the calender.
  !         The mismatch is about 3 days in 10000 years. This affects the computation of daylength in the current routine, and
  !         thus the "phase" of the temperature curve, i.e. that e.g. the time when the temperature starts rising in the morning,
  !         may be somewhat shifted.
  SUBROUTINE instantly_from_daily_temp_2(nblks, ice, lon, sinlat, coslat, dtmp_min, dtmp_max, date, tmp)

    !>>> JN-TODO: echam!
    USE mo_time_conversion, ONLY: time_days
    USE mo_time_control,    ONLY: get_year_day
    USE mo_orbit,           ONLY: inquire_declination
    !<<<

    IMPLICIT NONE

    REAL(wp), PARAMETER :: pio180 = pi/180._wp ! conversion factor from degrees to radians ("pi over 180")
    REAL(wp), PARAMETER :: rat1=14._wp/24._wp  ! day fraction at which maximum temperature shall occur
    REAL(wp), PARAMETER :: rat2=4._wp/24._wp   ! minimum length of normal day (as fraction of day)
    REAL(wp), PARAMETER :: rat3=24._wp/24._wp  ! maximum length of normal day (as fraction of day)
    REAL(wp), PARAMETER :: rat4=2._wp/24._wp   ! day fraction of dawn times (before sunrise + after sunset), used to set ...
                                               ! ... temperature of a normal day to T_mean at sunrise-rat4/2 and sunset+rat4/2

    INTEGER, INTENT(IN) :: nblks, ice
    REAL(wp), INTENT (IN) :: lon(:,:)       ! minimal temperature at the considered day
    REAL(wp), INTENT (IN) :: sinlat(:,:)    ! maximal range at that day
    REAL(wp), INTENT (IN) :: coslat(:,:)    ! maximal range at that day
    REAL(wp), INTENT (IN) :: dtmp_min(:,:)  ! minimal temperature at the considered day
    REAL(wp), INTENT (IN) :: dtmp_max(:,:)  ! maximal range at that day
    TYPE (time_days), INTENT (IN)  :: date  ! the considered simulation instant at that day in time-seconds format
    REAL(wp), INTENT (OUT) :: tmp(:,:)      ! estimated temperature


    INTEGER    :: yearday           ! the number of the day of the current timestep in the current year
    REAL(wp)   :: year_instant      ! the time of the current time step in yearday-fraction-of-day format
    REAL(wp)   :: frac_timestep     ! the current time step expressed as the fraction of the current day
    REAL(wp)   :: frac_localtime    ! local time at a particular grid point expressed as the fraction of the current ...
                                    ! ... day since midnight (mod 1
    REAL(wp)   :: frac_UTC_offset   ! offset between UTC and local time (depends on longitude)
    REAL(wp)   :: declination       ! solar declination
    REAL(wp)   :: cpds, spds        ! cosine and sine contributions to zenith angle
    REAL(wp)   :: frac_daylight     ! fraction of day between sunrise and sunset
    REAL(wp)   :: frac_since_sunset ! time between current time step and sunset expressed as fraction of day
    REAL(wp)   :: frac_sunrise      ! time of sunrise expressed as fraction of day since midnight
    REAL(wp)   :: frac_sunset       ! time of sunset expressed as fraction of day since midnight
    REAL(wp)   :: temp_sunset       ! temperature at sunset
    REAL(wp)   :: dtmp              ! mean day temperature
    REAL(wp)   :: dtran             ! day temperature range: maximum - minimum

    REAL(wp)   :: sd                ! modulation factor for temperature
    REAL(wp)   :: arg
    INTEGER    :: i, j

    CALL inquire_declination(declination)

    year_instant=get_year_day(date)            ! conversion of current time step from day-seconds to yearday-dayfraction format
    yearday = INT(year_instant)                ! number of day of the current time step since the beginning of the current year
    frac_timestep=year_instant-REAL(yearday,wp)   ! fractional part of the day of the current time step since midnight in UTC

    DO j = 1,nblks
      DO i = 1,ice
        frac_UTC_offset = lon(i,j)/360._wp   ! difference between local time and UTC as fraction of a day
        frac_localtime = MOD(frac_timestep+frac_UTC_offset+2._wp,1._wp) ! local time as fraction of the local day (modulo 1)
        spds = sinlat(i,j)*SIN(declination)  ! sine contribution to solar zenith angle
        cpds = coslat(i,j)*COS(declination)  ! cosine contribution to solar zenith angle
        arg = -spds/cpds                        ! faulty?: for cpds=0 arg gets infinite
        IF (arg > 1._wp) THEN               ! polar night:
          frac_daylight = 0._wp            ! ... day length is zero
        ELSEIF (arg < -1._wp) THEN         ! polar day:
          frac_daylight = 1._wp            ! ... day length is the whole day
        ELSE                                ! normal day / night:
          frac_daylight = ACOS(arg)/pi    ! ... day length expressed as fraction of the whole day
        ENDIF

        !! BLARPP empirical correction to account for the fact that the scheme does actually not
        !! consered dtmp=0.5*(max+min). The night-day turn temperature is adjusted such that the mean
        !! temperature is within 0.05 K of the 0.5*(tmax+tmin) temperature
        IF (frac_daylight>=rat2 ) THEN ! -------- normal day ---------------------------------------------
          IF(frac_daylight < 0.36_wp)THEN
            arg=-8._wp * (frac_daylight-0.36_wp)**2 + 0.685_wp
          ELSE
            arg=-0.8_wp * (frac_daylight-0.36_wp)**2 + 0.685_wp
          ENDIF
        ELSE
          arg = 0.5_wp
        ENDIF

        dtmp = arg * dtmp_max(i,j) + (1._wp-arg) * dtmp_min(i,j) ! The dusk day temperature
        !dtmp = 0.5_wp*(dtmp_max(i,j)+dtmp_min(i,j)) ! The mean day temperature
        dtran = dtmp_max(i,j)- dtmp_min(i,j)        ! temperature range
        IF (frac_daylight>=rat2 .AND. frac_daylight<=rat3) THEN ! -------- normal day ---------------------------------------------
          frac_sunrise = 0.5_wp - frac_daylight/2._wp             ! time of sunrise expressed as fraction of day since midnight
          frac_sunset = frac_sunrise+frac_daylight                ! time of sunset expressed as fraction of day since midnight
          IF (frac_localtime > frac_sunrise .AND. frac_localtime <= frac_sunset) THEN  ! For time steps at daylight ...
            sd = COS(pi*(frac_localtime-rat1)/(frac_daylight/2._wp+rat4)) ! ... modulation factor for day-temperature such that
                                                                    ! .. temperature is at maximum at rat1 (2 pm) and at average at
                                                                    ! .. sunrise-rat4/2 and sunset+rat4/2
            IF (sd >= 0) THEN
              tmp(i,j) = dtmp + sd*(dtmp_max(i,j)-dtmp)                       ! temperature modulated by temperature range
            ELSE
              tmp(i,j) = dtmp + sd*(dtmp-dtmp_min(i,j))                       ! temperature modulated by temperature range
            ENDIF
            !tmp(i,j) = dtmp + sd*dtran/2._wp                       ! temperature modulated by temperature range
          ELSE                                                 ! For time steps at night ...
            sd = COS(pi*(frac_sunset-rat1)/(frac_daylight/2._wp+rat4)) ! modulation factor for temperature at sunset
            IF (sd >= 0) THEN
              temp_sunset = dtmp + sd*(dtmp_max(i,j)-dtmp)                       ! temperature modulated by temperature range
            ELSE
              temp_sunset = dtmp + sd*(dtmp-dtmp_min(i,j))                       ! temperature modulated by temperature range
            ENDIF
            !temp_sunset = dtmp + dtran/2._wp*sd                     ! temperature at sunset
            frac_since_sunset = MOD(frac_localtime-frac_sunset+1._wp,1._wp) ! day fraction since sunset
            tmp(i,j) = dtmp_min(i,j)+(temp_sunset-dtmp_min(i,j)) &
              & *(1._wp-frac_since_sunset/(1._wp-frac_daylight))! linear interpolation
                                                                ! ... of temperature such that temperature is equal to the
                                                                ! ... temperature at sunset and drops until the minimum ...
                                                                ! ... temperature is met at sunrise.
          ENDIF
        ELSEIF (frac_daylight>rat3) THEN ! --- long day (close to polar day) ------------------------------------------------------
          ! case is now obsolete
          sd = COS(pi*(frac_localtime-rat1)/(frac_daylight/2._wp+rat4)) ! modulation factor for day-temperature such that ...
                                                                       ! ... temperature is at maximum at rat1 (2 pm) and at average
                                                                       ! ... at sunrise-rat4/2 and sunset+rat4/2
          tmp(i,j) = dtmp + dtran/2._wp*sd                             ! temperature modulated by temperature range
        ELSE ! ------------------------------- short day (close to polar night) ---------------------------------------------------
          tmp(i,j) = dtmp                                        ! temperature is set to mean temperture at all time steps
        ENDIF
      ENDDO
    ENDDO

    RETURN
  END SUBROUTINE INSTANTLY_FROM_DAILY_TEMP_2

  !>
  !! @brief calculate bottom pressure from air temperature and elevation
  !!
  !! @par History
  !! JN: copied from jsbach3 mo_jsbalone_forcing and only slightly adapted for jsb4 context
  !
  ! === bottom pressure() =======================================================================================================
  !
  ! To estimate bottom pressure a so called "polytrope atmosphere" is assumed for the troposhere. Here a fixed temperature
  ! gradient "gamma" is prescribed, i.e. air temperature at height z is set to
  !     (1)     T(z)=T0-gamma*z,
  ! where T0 is the temperature at the bottom. Assuming further hydrostatic conditions the vertical pressure gradient is given by
  !              wp
  !     (2)     ---- = -g*rho,
  !              dz
  ! where rho is the air density and g the acceleration due to gravity. Pressure and Temperature are related by the state
  ! equation of an ideal gas
  !     (3)      p = rho*R*T,
  ! where R is the specific gas constant and T the temperature. By eliminating rho from (3) and (2) and entering (1) for T one
  ! obtains
  !              wp        p*g
  !     (4)     ---- = - --------------
  !              dz      R(T0-gamma*z)
  ! By integration one obtains the final result
  !                              gamma                              g
  !     (5)      p = p0 * (1 - z*-----)^k , with the exponent k= ------- ,
  !                               T0                             R*gamma
  ! where p0 and T0 are the bottom pressure and temperature. It should be noted that the expression in brackets can get negative
  ! only for extreme values of T0 and z. Assuming e.g. T0 > -80 Celsius (approx. 200 Kelvin), the bracket stays positive for
  ! z<30km, which is true even at the Mount Everest. (Eq. (5) coincides with Eq. (31) from Knorr for z << 1!)
  !
  ! If temperature T(z) at height z is known instead of T0, using Eq. (1) Eq. (5) translates to:
  !                                  gamma
  !     (6)      p = p0 * (1 - z*------------)^k
  !                              T(z)+z*gamma
  !
  ELEMENTAL FUNCTION bottom_pressure(elevation,air_temperature)
    REAL(wp)              :: bottom_pressure
    REAL(wp), INTENT(IN)  :: elevation
    REAL(wp), INTENT(IN)  :: air_temperature

    ! in jsb3: p0sl_bg = p_sealevel; grav = Gravity; rd = GasConstantDryAir
    bottom_pressure = &
      & p0sl_bg * ( 1._wp - elevation*gamma /(tmelt + air_temperature + elevation*gamma))**(grav/rd/gamma)
  END FUNCTION bottom_pressure

  !>
  !! @brief split precipitation in snow and rain depending on mean air temperature
  !!
  !! @par History
  !! JN: copied from jsbach3 mo_jsbalone_forcing
  !
  ! === snow_and_rain_from_precip() ==============================================================================================
  !
  ! Precipitation is divided into rain and snow according to a method by M.S. Wigmosta, L. Vail and D.P. Lettenmaier, "A 
  ! distributed hydrology-vegetation model for complex terrain", Water Resources Research 30 (1994) 1665-1679: depending on 
  ! average day temperature "T", the snow part of precipitation "P_sn" is computed from the total precipitation "P" as:
  !               / P               for T < -1.1 Celsius
  ! (1)   P_sn = |  P*(3.3 - T)/4.4 for -1.1 Celsius <= T <= 3.3 Celsius
  !               \ 0               for T > 3.3 Celsius
  ! Hence the rain part ist then
  ! (2)   P_rain = P-P_sn.
  !
  ELEMENTAL SUBROUTINE snow_and_rain_from_precip(rain,snow,precipitation,air_temp_daily)
    REAL(wp),INTENT(OUT) :: rain          !! rain [kg/m^2/s]
    REAL(wp),INTENT(OUT) :: snow          !! snow [kg/m^2/s]
    REAL(wp),INTENT(IN)  :: precipitation !! total precipitation (rain+snow) [kg/m^2/s]
    REAL(wp),INTENT(IN)  :: air_temp_daily!! daily mean air_temperature at surface [Celsius]

    IF(air_temp_daily < -1.1_wp) THEN
       snow = precipitation
    ELSEIF(air_temp_daily > 3.3_wp) THEN
       snow = 0.0_wp
    ELSE
       snow = precipitation*(3.3_wp-air_temp_daily)/4.4_wp
    ENDIF
    ! set the rain part according to Eq. (2):
    rain = precipitation - snow
  END SUBROUTINE snow_and_rain_from_precip

  !>
  !! @brief computes PAR and NIR, together with the diffuse part of PAR, from cloudcover
  !!
  !! @par History
  !! JN: copied from jsbach3 mo_jsbalone_forcing, removed not required output
  !
  ! === shortwave_from_cloudcover() ================================================================================================
  !
  ! Following the strategy described in [1], this routine computes PAR and NIR, together with the diffuse part of PAR, from
  ! cloudcover, by converting cloudcover "n" to "fpar" (fpar is the fraction of actual to potential radiation flux in the PAR
  ! band incident at the earth's surface):
  ! (1)           fpar = 0.5 + 0.4 * n.
  ! Then the main computations are done in the subroutine "shortwave_from_fpar()" (see there for a description of the method).
  !
  ! Literature:
  ! [1] W. Knorr, "Satellite remote sensing and modelling of the global CO2 exchange of land vegetation", Examensarbeit 49,
  ! (Max Planck Institut fuer Meteorologie, Hamburg, 1998).
  !
  ELEMENTAL SUBROUTINE shortwave_from_cloudcover(day_of_year, cos_zenith, air_pressure, cloud_cover, & !! input
    & rad_UV_down, rad_PAR_down, rad_NIR_down, frac_PAR_diffuse, rad_sw_down, rad_sw_down_pot) !! output

    INTEGER,INTENT(in)   :: day_of_year      !! day in year (from [1,365])
    REAL(wp),INTENT(in)  :: cos_zenith       !! cosine of zenith angle
    REAL(wp),INTENT(in)  :: air_pressure     !! air pressure at bottom (depending on local elevation) [N/m^2]
    REAL(wp),INTENT(in)  :: cloud_cover      !! fraction of sky covered by clouds (from [0,1])
    REAL(wp),INTENT(out) :: rad_UV_down      !! solar radiation flux (direct+diffuse) from the UV band 280-400 nm ...
                                             !! .. arriving at the bottom (includes cloud shading and scattering) [W/m^2]
    REAL(wp),INTENT(out) :: rad_PAR_down     !! solar radiation flux (direct+diffuse) from the visible band 400-700 nm ...
                                             !! .. arriving at the bottom (includes cloud shading and scattering) [W/m^2]
    REAL(wp),INTENT(out) :: rad_NIR_down     !! solar radiation flux (direct+diffuse) from the near infrared band 700-3000 nm ...
                                             !! .. arriving at the bottom (includes cloud shading and scattering) [W/m^2]
    REAL(wp),INTENT(out) :: frac_PAR_diffuse !! Fraction from  rad_PAR_down() that comes down as diffuse radiation ..
                                             !! .. (values in [0,1])
    REAL(wp),INTENT(out) :: rad_sw_down      !! total solar shortwave radiation : direct+diffuse, 280-3000 nm, (i.e. PAR+NIR), ..
                                             !! ... that actually arrives at the bottom [W/m^2]
    REAL(wp),INTENT(out) :: rad_sw_down_pot  !! total solar shortwave radiation : direct+diffuse, 280-3000 nm, (i.e. PAR+NIR), ..
                                             !! ... that actually arrives at the bottom [W/m^2]

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':shortwave_from_cloudcover'

    REAL(wp) :: fpar

    fpar = 0.5_wp + 0.4_wp * (1.0_wp - cloud_cover)
    SELECT CASE(forcing_options%type_of_shortwave_scheme)
    CASE(SW_SCHEME_ORIGINAL_)
      CALL shortwave_from_fpar_orig(day_of_year, cos_zenith, air_pressure, fpar,                            & !! inputs
        & rad_UV_down, rad_PAR_down, rad_NIR_down, frac_PAR_diffuse, rad_sw_down)   !! outputs
    CASE(SW_SCHEME_CCDAS_)
      CALL shortwave_from_fpar_ccdas(day_of_year, cos_zenith, air_pressure, fpar,                           & !! inputs
        & rad_UV_down, rad_PAR_down, rad_NIR_down, frac_PAR_diffuse, rad_sw_down)   !! outputs
    END SELECT
    rad_sw_down_pot = rad_sw_down / fpar

  END SUBROUTINE shortwave_from_cloudcover

  !>
  !! @brief computes PAR and NIR, together with the diffuse part of PAR, from fpar
  !!
  !! @par History
  !! JN: copied from jsbach3 mo_jsbalone_forcing
  !
  ! === shortwave_from_fpar_orig() ===============================================================================================
  !
  ! This routine follows exactly [1]. It computes PAR and NIR, together with the diffuse part of PAR, from "fpar", which is the
  ! quotient of the actual PAR flux "Rpar_act" (i.e. at cloudy sky) to the potential PAR flux "Rpar_pot" (i.e. at clear sky), both
  ! measured at the ground:
  !                       Rpar_act
  ! (1)           fpot := -------- ,
  !                       Rpar_pot
  ! Following [1] the various components of the total shortwave radiation are computed as follows: The total flux in the PAR band
  ! incident to the outer atmosphere "Rpar_top" is computed from
  ! (2)           Rpar_top = S0*cos(theta)*E,
  ! where S0=0.44*1360 W/m^2 =~= 600 W/m^2 is the PAR part of the solar constant (0.4-0.7 nm), "theta" is the zenith angle
  ! and
  ! (3)           E = (r0/r)^2 = 1.000110+0.034221 cos(psi)+0.001280*sin(psi)+0.000719*cos(2*psi)+0.000077*sin(2*psi)
  ! the inverse squared distance earth-sun "r" in units of the average distance "r0", where "psi" is the day-angle, i.e. the
  ! day of the year expressed as an angle (see e.g. [4], Eq. (3.3), [1] Eq. (24)).
  ! The direct part of potential PAR reaching ground "Rpar_pot_dir" is according to [2] approximately given by
  ! (4)          Rpar_pot_dir = Rpar_top*exp(-0.185*(p/p0)/cos(theta)),
  ! where "p" is the current and "p0" the sealevel pressure (see [1] Eqs. (27)-(29) and [2] Eqs. (1)-(3)), and the exponential
  ! term accounts for direct radiation reduced by Rayleigh scattering ( p/p0/cos(theta) is a measure for the mass of air along
  ! the beam ). According to [2] the diffuse part of potential PAR ("Rpar_pot_diff") reaching ground is approximately 40% of the
  ! difference between Rpar_top and the direct radiation (Equation (3) n [2] contains a printing error?):
  ! (5)           Rpar_pot_diff = 0.4*(Rpar_top - Rpar_pot_dir)
  ! so that potential PAR is
  ! (6)           Rpar_pot = Rpar_pot_dir + Rpar_pot_diff = Rpar_top*[0.4 + 0.6*exp(-0.185*(p/p0)/cos(theta))].
  ! By (1) the actual PAR then follows as
  ! (7)           Rpar_act = fpar * Rpar_pot.
  ! Next the near infrared (NIR) part of shortwave radiation has to be accounted for. Following [3] this is done by introducing
  ! a correction factor "F", such that the total downwards shortwave radiation "Rtot" reaching bottom is obtained from Rpar_act by
  ! (8)           Rtot = Rpar_act/F,
  ! where, using mfpar=1-fpar, the inverse factor is given by
  ! (9)           1/F = 1 + (1.185-0.437*mfpar-0.494*mfpar^2)* exp[(0.0305-0.208*mfpar+0.495*mfpar^2)/cos(theta)]
  ! Therefore the downward radiation reaching the surface in the near infrared band "Rnir" is
  ! (10)          Rnir = Rtot - Rpar_act
  ! Finally the fraction of direct radiation in PAR at the actual PAR
  !                           Rpar_act_dir
  ! (11)          fdirPar =   ------------
  !                             Rpar_act
  ! is estimated from the direct part in potential PAR
  ! (11)          Rpar_pot_dir = S0 * exp(-0.185*(p/p0)/cos(theta))
  ! total potential PAR "Rpar_pot" and "fpar" by (see Knorr Eq. (30))
  !                            /  0 for fpar < 0.2
  !                           /                   2/3
  !                          /      / 0.9 - fpar \      Rpar_pot_dir
  ! (12)          fdirPar = (  (1- (--------------)   ) ------------       for 0.2 <= fpar <= 0.9
  !                          \      \    0.7     /         Rpar_pot
  !                           \
  !                            \  1 for fpar > 0.9
  ! where the quotient showing up in (12) is
  !
  !               Rpar_pot_dir         exp(-0.185*(p/p0)/cos(theta))
  ! (13)          ------------ =   ---------------------------------------
  !                 Rpar_pot       0.4 + 0.6*exp(-0.185*(p/p0)/cos(theta))
  !
  ! It follows that the fraction of diffuse radiation "fdiffPar" in PAR is
  !
  ! (14)         fdiffPar = 1 - fdirPar.
  !
  ! REMARK: Eq. (12) is taken from [2] but contains an additional approximation: In [2] fdirPar does not depend on fpar, but on
  ! the transmission ratio for TOTAL SHORTWAVE flux, whereas fpar is defined as the transmission ratio for PAR only (see Eq. (1)).
  ! Hence, (12) includes the additional approximation of the shortwave transmission ratio by the PAR transmission ratio (fpar).
  !
  ! Literature:
  ! [1] W. Knorr, "Satellite remote sensing and modelling of the global CO2 exchange of land vegetation", Examensarbeit 49,
  ! (Max Planck Institut fuer Meteorologie, Hamburg, 1998).
  ! [2] A. Weiss & J.M. Norman, "Partitioning solar radiation into direct and diffuse visible and near-infrared components",
  ! Agricultural and Forest Meteorology, 34 (1985) 205-213.
  ! [3] R.T. Pinker & I. Laszlo, "Global distribution of photosynthetically active radiation as observed from satellites",
  ! J. Climate 5 (1992) 56-65.
  ! [4] G.W. Paltridge and C.M.R. Platt, "Radiative Processes in Meteorology and Climatology" (Elsevier, Amsterdam, 1976).
  ELEMENTAL SUBROUTINE shortwave_from_fpar_orig(day_of_year,cos_zenith,air_pressure,fpar, & !! inputs
    & rad_UV,rad_PAR,rad_NIR,frac_PAR_diffuse,rad_sw)  !! outputs

    INTEGER,INTENT(in)   :: day_of_year      !! day in year (from [1,365])
    REAL(wp),INTENT(in)  :: cos_zenith       !! cosine of zenith angle
    REAL(wp),INTENT(in)  :: air_pressure     !! air pressure at bottom (depending on local elevation) [N/m^2]
    REAL(wp),INTENT(in)  :: fpar             !! daily mean fraction of actual to potential PAR at bottom (from [0,1])
    REAL(wp),INTENT(out) :: rad_UV           !! "RUV": solar radiation flux (direct+diffuse) from the UV band 280-700 nm .
                                             !! .. arriving actually at the bottom (includes cloud shading and scattering) [W/m^2]
    REAL(wp),INTENT(out) :: rad_PAR          !! "Rpar_act": solar radiation flux (direct+diffuse) from the visible band 400-700 nm .
                                             !! .. arriving actually at the bottom (includes cloud shading and scattering) [W/m^2]

    REAL(wp),INTENT(out) :: rad_NIR          !! "Rnir": flux (direct+diffuse) from the near infrared band 700-3000 nm ..
                                             !! .. actually arriving at the bottom (includes cloud shading and scattering) [W/m^2]
    REAL(wp),INTENT(out) :: frac_PAR_diffuse !! "fdiffPar": Fraction from  rad_PAR() that comes down as diffuse radiation ..
                                             !! .. (values in [0,1])
    REAL(wp),INTENT(out) :: rad_sw           !! "Rtot": total solar shortwave radiation : direct+diffuse, 280-3000 nm, ,..
                                             !! ... (i.e. PAR+NIR) that actually arrives at the bottom [W/m^2]

    REAL(wp) :: day_angle     !! day of year expressed as angle in radians
    REAL(wp) :: r_sun_earth_2 !! inverse squared earth-sun-distance normalized to mean distance
    REAL(wp) :: Rpar_top      !! PAR flux incident to the outer atmosphere at particular zenith angle [W/m^2]
    REAL(wp) :: F_invers      !! factor 1/F to obtain total shortwave from visible radiation (see Eq. (9))
    REAL(wp) :: hlp_r,hlp2_r,hlp3_r,hlp4_r

    REAL(wp), PARAMETER  ::  frac_UV = 0._wp ! UV fraction of (UV + PAR)
                                             ! to be checked, if the ratio rad_PAR/rad_sw is within observed range (0.42-0.51) with this value of frac_UV

    day_angle = 2._wp*pi*(day_of_year-1._wp)/365._wp !! ATTENTION: this may be incorrect for paleo situations
    hlp_r=2._wp*day_angle
    r_sun_earth_2 = &         ! Eq. (3) from above
      & 1.000110_wp+3.4221E-2_wp*COS(day_angle)+1.280E-3_wp*SIN(day_angle)+7.19E-4_wp*COS(hlp_r)+7.7E-5_wp*SIN(hlp_r)

    IF(cos_zenith > 0.017452406_wp) THEN! zenith angle smaller than 89 degrees  (cos(89) = 0.017452406)
      Rpar_top   = solar_const_vis*cos_zenith*r_sun_earth_2       ! Eq. (2) from above
      hlp3_r = EXP(-0.185_wp*(air_pressure/p0sl_bg)/cos_zenith) ! in jsb3: p0sl_bg = p_sealevel
      hlp4_r = 0.4_wp+0.6_wp*hlp3_r
      rad_UV = fpar * Rpar_top * hlp4_r * frac_UV
      rad_PAR = fpar * Rpar_top * hlp4_r * (1._wp - frac_UV)   ! This is Rpar_act; Eq. (7) from above
      hlp_r=1._wp-fpar
      hlp2_r=hlp_r*hlp_r
      F_invers=1._wp + (1.185_wp-0.437_wp*hlp_r-0.494_wp*hlp2_r)* &
        & EXP(MIN(2.534_wp,(0.0305_wp-0.208_wp*hlp_r+0.495_wp*hlp2_r)/cos_zenith)) ! Eq. (9) from above
      rad_sw = (rad_PAR + rad_UV) *F_invers   ! This is "Rtot"; Eq. (8) from above
      rad_NIR = rad_sw - rad_PAR - rad_UV ! This is "Rnir"; Eq. (10) from above
      IF(fpar <= 0.2_wp) THEN
        frac_PAR_diffuse = 1.                                                     ! Eq. (12)-(15) from above
      ELSEIF(fpar >= 0.9_wp) THEN
        frac_PAR_diffuse = 0._wp                                                  ! Eq. (12)-(15) from above
      ELSE
        frac_PAR_diffuse = 1._wp - &
          & (1._wp- ((0.9_wp -fpar)/0.7_wp)**(2._wp/3._wp))*hlp3_r/hlp4_r ! Eq. (12)-(15) from above
      ENDIF
    ELSE ! zenith angle larger than 89 degrees
      rad_UV = 0._wp
      rad_PAR = 0._wp
      rad_NIR = 0._wp
      frac_PAR_diffuse = 1._wp ! (could also be another value in this case)
      rad_sw  = 0._wp
    ENDIF
  END SUBROUTINE shortwave_from_fpar_orig


  !>
  !! @brief computes shortwave parts (PAR,NIR,fraction of diffuse PAR)
  !!
  !! @par History
  !! JN: copied from jsbach3 mo_jsbalone_forcing
  !
  ELEMENTAL SUBROUTINE shortwave_from_direct_shortwave(day_of_year,cos_zenith,cos_lat,sin_lat,air_pressure, &
    & rad_sw_in, rad_UV,rad_PAR,rad_NIR,frac_PAR_diffuse,rad_sw,rad_sw_pot )
    ! similar concept as above, however, the observed radiation is used to derive the remaining
    ! parameters. In daily version, the instantaneous flux at noon required to calculate irradiation at any
    ! time of the day is approximated by solving the integral over the positive sun angle, following Prentice et al. Ecol. Mod. 1993

    INTEGER,INTENT(in)   :: day_of_year      !! day in year (from [1,365])
    REAL(wp),INTENT(in)  :: cos_zenith       !! cosine of zenith angle
    REAL(wp),INTENT(in)  :: cos_lat          !! cosine of latitude
    REAL(wp),INTENT(in)  :: sin_lat          !! sine of latitude
    REAL(wp),INTENT(in)  :: air_pressure     !! air pressure at bottom (depending on local elevation) [N/m^2]
    REAL(wp),INTENT(in)  :: rad_sw_in        !! "Rtot": daily average total solar shortwave radiation : direct+diffuse, 280-3000 nm,
                                             !! ... (i.e. PAR+NIR) that actually arrives at the bottom [W/m^2]
    REAL(wp),INTENT(out) :: rad_UV           !! "RUV": solar radiation flux (direct+diffuse) from the UV band 280-400 nm
                                             !! .. arriving actually at the bottom (includes cloud shading and scattering) [W/m^2]
    REAL(wp),INTENT(out) :: rad_PAR          !! "Rpar_act": solar radiation flux (direct+diffuse) from the visible band 400-700 nm .
                                             !! .. arriving actually at the bottom (includes cloud shading and scattering) [W/m^2]

    REAL(wp),INTENT(out) :: rad_NIR          !! "Rnir": flux (direct+diffuse) from the near infrared band 700-3000 nm ..
                                             !! .. actually arriving at the bottom (includes cloud shading and scattering) [W/m^2]
    REAL(wp),INTENT(out) :: frac_PAR_diffuse !! "fdiffPar": Fraction from  rad_PAR() that comes down as diffuse radiation ..
                                             !! .. (values in [0,1])
    REAL(wp),INTENT(out) :: rad_sw,rad_sw_pot

    REAL(wp) :: day_angle     !! day of year expressed as angle in radians
    REAL(wp) :: r_sun_earth_2 !! inverse squared earth-sun-distance normalized to mean distance
    REAL(wp) :: Rpar_top      !! PAR flux incident to the outer atmosphere at particular zenith angle [W/m^2]
    REAL(wp) :: fsw           !! "fsw": fraction of potential radiation solar shortwave radiation : direct+diffuse, 280-3000 nm, ,..
                              !! ... (i.e. PAR+NIR) that actually arrives at the bottom [W/m^2]
    REAL(wp) :: F_invers      !! factor 1/F to obtain total shortwave from visible radiation (see Eq. (9))
    REAL(wp) :: hlp_r,hlp2_r,hlp3_r,hlp4_r
    REAL(wp) :: rad_sw_in_corr
    REAL(wp) :: delta,v,u,hh,sinehh,K

    REAL(wp), PARAMETER  ::  frac_UV = 0.1_wp ! UV fraction of (UV + PAR)

    day_angle = 2._wp*pi*(day_of_year+10._wp)/365._wp !! ATTENTION: this may be incorrect for paleo situations

    hlp_r=2._wp*day_angle
    r_sun_earth_2 = &         ! Eq. (3) from above
      & 1.000110_wp+3.4221E-2_wp*COS(day_angle)+1.280E-3_wp*SIN(day_angle)+7.19E-4_wp*COS(hlp_r)+7.7E-5_wp*SIN(hlp_r)
    K=13750.98708_wp

    IF (forcing_input%shortwave%frequency .NE. TIMESTEP_) THEN
      ! today's sum of top of the atmosphere radiation
      delta=-23.4_wp*pi/180._wp*COS(day_angle)
      u=sin_lat*SIN(delta)
      v=cos_lat*COS(delta)
      IF (u.GE.v) THEN
         hh=pi
      ELSEIF (u.LE.(0._wp-v)) THEN
         hh=0.0_wp
      ELSE
         hh=ACOS(-u/v)
      ENDIF
      sinehh=SIN(hh)

      !estimate of noon radiation giving rise to the recorded radiation sum
      rad_sw_in_corr = 0._wp
      IF(u*hh+v*sinehh.GT.1e-5_wp)THEN
         rad_sw_in_corr = rad_sw_in / (2._wp*(u*hh+v*sinehh) * K) * 86400._wp
      ENDIF
      rad_sw = rad_sw_in_corr * cos_zenith

    ELSE
      ! time step forcing, use radiation directly
      rad_sw = rad_sw_in
    ENDIF

    IF(cos_zenith > 0.017452406_wp) THEN! zenith angle smaller than 89 degrees  (cos(89) = 0.017452406)
      Rpar_top   = solar_const*cos_zenith*r_sun_earth_2       ! Eq. (2) from above
      hlp3_r = EXP(-0.185_wp*(air_pressure/p0sl_bg)/cos_zenith) ! in jsb3: p0sl_bg = p_sealevel
      hlp4_r = 0.4_wp+0.6_wp*hlp3_r
      rad_sw_pot = Rpar_top * hlp4_r
      fsw = MAX(MIN(rad_sw/rad_sw_pot,1._wp),0._wp)
      hlp_r=1._wp-fsw
      hlp2_r=hlp_r*hlp_r
      F_invers=1._wp + (1.185_wp-0.437_wp*hlp_r-0.494_wp*hlp2_r)* &
        & EXP(MIN(2.534_wp,(0.0305_wp-0.208_wp*hlp_r+0.495_wp*hlp2_r)/cos_zenith)) ! Eq. (9) from above
      rad_UV = rad_sw * frac_UV / F_invers
      rad_PAR = rad_sw * (1._wp - frac_UV) / F_invers   ! This is "Rtot"; Eq. (8) from above
      rad_NIR = rad_sw - rad_PAR - rad_UV ! This is "Rnir"; Eq. (10) from above
      IF(fsw <= 0.2_wp) THEN
         frac_PAR_diffuse = 1._wp                                                  ! Eq. (12)-(15) from above
      ELSEIF(fsw >= 0.9_wp) THEN
         frac_PAR_diffuse = 0.08_wp                                                ! Eq. (12)-(15) from above
      ELSE
         frac_PAR_diffuse = 1._wp - &
           & (1._wp- ((0.9_wp -fsw)/0.7_wp)**(2._wp/3._wp))*hlp3_r/hlp4_r ! Eq. (12)-(15) from above
      ENDIF
    ELSE ! zenith angle larger than 89 degrees
      rad_UV = 0._wp
      rad_PAR = 0._wp
      rad_NIR = 0._wp
      frac_PAR_diffuse = 1._wp ! (could also be another value in this case)
      rad_sw  = 0._wp
      rad_sw_pot  = 0._wp
    ENDIF

  END SUBROUTINE shortwave_from_direct_shortwave


  !>
  !! @brief computes PAR and NIR, together with the diffuse part of PAR, from fpar (according to ccdas)
  !!
  !! @par History
  !! JN: copied from jsbach3 mo_jsbalone_forcing
  !
  ! === shortwave_from_fpar_ccdas() ================================================================================================
  !
  ! This routine implements the shortwave radiation scheme of ccdas. As far as I understand it is almost identical to
  ! to shortwave_from_fpar_orig(), with the  only exceptions that Equation (9) is replaced by the more simple relationsship  !
  !
  ! (9')        F= 0.43 +0.25*(1.-fpar)*cos(theta)
  ! (whereever this comes from).
  !
  ELEMENTAL SUBROUTINE shortwave_from_fpar_ccdas(day_of_year,cos_zenith,air_pressure,fpar, & !! inputs
    & rad_UV,rad_PAR,rad_NIR,frac_PAR_diffuse,rad_sw)  !! outputs

    INTEGER,INTENT(in)   :: day_of_year      !! day in year (from [1,365])
    REAL(wp),INTENT(in)  :: cos_zenith       !! cosine of zenith angle
    REAL(wp),INTENT(in)  :: air_pressure     !! air pressure at bottom (depending on local elevation) [N/m^2]
    REAL(wp),INTENT(in)  :: fpar             !! daily mean fraction of actual to potential PAR at bottom (from [0,1])
    REAL(wp),INTENT(out) :: rad_UV           !! "RUV_act": solar radiation flux (direct+diffuse) from the UV band 280-400 nm .
                                             !! .. arriving actually at the bottom (includes cloud shading and scattering) [W/m^2]
    REAL(wp),INTENT(out) :: rad_PAR          !! "Rpar_act": solar radiation flux (direct+diffuse) from the visible band 400-700 nm .
                                             !! .. arriving actually at the bottom (includes cloud shading and scattering) [W/m^2]

    REAL(wp),INTENT(out) :: rad_NIR          !! "Rnir": flux (direct+diffuse) from the near infrared band 700-3000 nm ..
                                             !! .. actually arriving at the bottom (includes cloud shading and scattering) [W/m^2]
    REAL(wp),INTENT(out) :: frac_PAR_diffuse !! "fdiffPar": Fraction from  rad_PAR() that comes down as diffuse radiation ..
                                             !! .. (values in [0,1])
    REAL(wp),INTENT(out) :: rad_sw           !! "Rtot": total solar shortwave radiation : direct+diffuse, 280-3000 nm, ,..
                                             !! ... (i.e. PAR+NIR) that actually arrives at the bottom [W/m^2]

    REAL(wp) :: day_angle     !! day of year expressed as angle in radians
    REAL(wp) :: r_sun_earth_2 !! inverse squared earth-sun-distance normalized to mean distance
    REAL(wp) :: Rpar_top      !! PAR flux incident to the outer atmosphere at particular zenith angle [W/m^2]
    REAL(wp) :: F             !! factor F to obtain total shortwave from visible radiation (see Eq. (9'))
    REAL(wp) :: hlp_r, hlp3_r,hlp4_r

    REAL(wp), PARAMETER  ::  frac_UV = 0._wp ! UV fraction of (UV + PAR)
                                             ! to be checked, if the ratio rad_PAR/rad_sw is within observed range (0.42-0.51) with this value of frac_UV

    day_angle = 2._wp*pi*(day_of_year-1._wp)/365._wp !! ATTENTION: this may be incorrect for paleo situations
    hlp_r=2._wp*day_angle
    r_sun_earth_2 = &         ! Eq. (3) from above
      & 1.000110_wp+3.4221E-2_wp*COS(day_angle)+1.280E-3_wp*SIN(day_angle)+7.19E-4_wp*COS(hlp_r)+7.7E-5_wp*SIN(hlp_r)
    IF(cos_zenith > 0.017452406_wp) THEN! zenith angle smaller than 89 degrees  (cos(89) = 0.017452406)
      Rpar_top   = solar_const_vis*cos_zenith*r_sun_earth_2       ! Eq. (2) from above
      hlp3_r = EXP(-0.185_wp*(air_pressure/p0sl_bg)/cos_zenith) ! in jsb3: p0sl_bg = p_sealevel
      hlp4_r = 0.4_wp+0.6_wp*hlp3_r
      rad_UV = fpar * Rpar_top * hlp4_r * frac_UV
      rad_PAR = fpar * Rpar_top * hlp4_r * (1._wp - frac_UV)    ! This is Rpar_act; Eq. (7) from above
      hlp_r=1.-fpar
      F = 0.43_wp +0.25_wp*(1.-fpar)*cos_zenith ! Eq. (9') from above
      rad_sw = (rad_PAR + rad_UV)/ F   ! This is "Rtot"; Eq. (8) from above
      rad_NIR = rad_sw - rad_PAR - rad_UV ! This is "Rnir"; Eq. (10) from above
      IF(fpar <= 0.2_wp) THEN
         frac_PAR_diffuse = 1._wp                                                     ! Eq. (12)-(14) from above
      ELSEIF(fpar >= 0.9_wp) THEN
         frac_PAR_diffuse = 0._wp                                                     ! Eq. (11)-(14) from above
      ELSE
         frac_PAR_diffuse = 1._wp - &
           & (1._wp- ((0.9_wp -fpar)/0.7_wp)**(2._wp/3._wp))*hlp3_r/hlp4_r ! Eq. (11)-(14) from above
      ENDIF
    ELSE ! zenith angle larger than 89 degrees
      rad_UV = 0._wp
      rad_PAR = 0._wp
      rad_NIR = 0._wp
      frac_PAR_diffuse = 1._wp ! (could also be another value in this case)
      rad_sw  = 0._wp
    ENDIF
  END SUBROUTINE shortwave_from_fpar_ccdas

  !>
  !! @brief computes PAR and NIR, together with the diffuse part of PAR, from fpar (according to ccdas but with Wolfgangs modification)
  !!
  !! @par History
  !! JN: copied from jsbach3 mo_jsbalone_forcing
  !
  ! === shortwave_from_fsw() =====================================================================================================
  !
  ! This routine is identical to the shortwave routine from CCDAS shortwave_from_fpar_ccdas(), but with a modification proposed by
  ! Wolfgang to use not fpar as input, but the transmission ratio fsw for TOTAL SHORTWAVE radiation instead.
  !
  ! Here
  !                     Rtot
  ! (A)        fsw := ----------
  !                     R_pot
  !
  ! the ratio between actual "Rtot" to potential "R_pot" shortwave radiation (daily means). "R_pot" is obtailed from the
  ! potential PAR radiation by
  !                      Rpar_pot
  ! (B)        R_pot = ---------
  !                       F_pot
  ! where the conversion factor "F_pot" is taken as
  ! (C)        F_pot = 0.43.
  ! With (1) then total incident shortwave radiation follows from
  !                                       Rpar_pot
  ! (D)        Rtot = fsw * R_pot = fsw * ---------
  !                                         F_pot
  ! The PAR part of total actual radiation is then obtained by
  ! (E)        Rpar_act = Rtot * F,
  ! where the conversion factor "F" is
  ! (F)        F= 0.43 +0.25*(1.-fpar)*cos(theta)
  ! Now by equations (10) to (14) from shortwave_from_fpar_orig() all output fields can be computed.
  !
  ELEMENTAL SUBROUTINE shortwave_from_fsw(day_of_year,cos_zenith,air_pressure,fsw, & !! inputs
    & rad_UV,rad_PAR,rad_NIR,frac_PAR_diffuse,rad_sw) !! outputs

    INTEGER,INTENT(in)   :: day_of_year      !! day in year (from [1,365])
    REAL(wp),INTENT(in)  :: cos_zenith       !! cosine of zenith angle
    REAL(wp),INTENT(in)  :: air_pressure     !! air pressure at bottom (depending on local elevation) [N/m^2]
    REAL(wp),INTENT(in)  :: fsw              !! daily mean fraction of actual to potential TOTAL SHORTWAVE flux at bottom ..
                                             !! .. (from [0,1])
    REAL(wp),INTENT(out) :: rad_UV           !! "RUV_act": solar radiation flux (direct+diffuse) from the UV band 280-400 nm
                                             !! .. arriving actually at the bottom (includes cloud shading and scattering) [W/m^2]
    REAL(wp),INTENT(out) :: rad_PAR          !! "Rpar_act": solar radiation flux (direct+diffuse) from the visible band 400-700 nm .
                                             !! .. arriving actually at the bottom (includes cloud shading and scattering) [W/m^2]

    REAL(wp),INTENT(out) :: rad_NIR          !! "Rnir": flux (direct+diffuse) from the near infrared band 700-3000 nm ..
                                             !! .. actually arriving at the bottom (includes cloud shading and scattering) [W/m^2]
    REAL(wp),INTENT(out) :: frac_PAR_diffuse !! "fdiffPar": Fraction from  rad_PAR() that comes down as diffuse radiation ..
                                             !! .. (values in [0,1])
    REAL(wp),INTENT(out) :: rad_sw           !! "Rtot": total solar shortwave radiation : direct+diffuse, 400-3000 nm, ,..
                                             !! ... (i.e. PAR+NIR) that actually arrives at the bottom [W/m^2]

    REAL(wp) :: day_angle     !! day of year expressed as angle in radians
    REAL(wp) :: r_sun_earth_2 !! inverse squared earth-sun-distance normalized to mean distance
    REAL(wp) :: Rpar_top      !! PAR flux incident to the outer atmosphere at particular zenith angle [W/m^2]
    REAL(wp) :: F             !! factor F to obtain total shortwave from visible radiation (see Eq. (9'))
    REAL(wp) :: hlp_r,hlp3_r,hlp4_r

    REAL(wp),PARAMETER :: F_pot = 0.43_wp !! See Eq. (C) from above
    REAL(wp),PARAMETER :: frac_UV = 0._wp !! UV fraction of (UV + PAR)
                                          !! to be checked, if the ratio rad_PAR/rad_sw is within observed range (0.42-0.51) with this value of frac_UV

    day_angle = 2._wp*pi*(day_of_year-1._wp)/365._wp !! ATTENTION: this may be incorrect for paleo situations
    hlp_r=2._wp*day_angle
    r_sun_earth_2 = &         ! Eq. (3) from above
      & 1.000110_wp+3.4221E-2_wp*COS(day_angle)+1.280E-3_wp*SIN(day_angle)+7.19E-4_wp*COS(hlp_r)+7.7E-5_wp*SIN(hlp_r)
    IF(cos_zenith > 0.017452406_wp) THEN! zenith angle smaller than 89 degrees  (cos(89) = 0.017452406)
       Rpar_top   = solar_const_vis*cos_zenith*r_sun_earth_2       ! Eq. (2) from above
       hlp3_r = EXP(-0.185_wp*(air_pressure/p0sl_bg)/cos_zenith) ! in jsb3: p0sl_bg = p_sealevel
       hlp4_r = 0.4_wp+0.6_wp*hlp3_r
       rad_sw = fsw*Rpar_top*hlp4_r/F_pot   ! This is "Rtot"; Eq. (D) from above
       hlp_r=1.-fsw
       F = 0.43_wp + 0.25_wp*(1.-fsw)*cos_zenith ! Eq. (9') from shortwave_from_fpar_ccdas(), but with the correct ratio
                                                 ! fsw instead fpar
       rad_UV = rad_sw * F * frac_UV
       rad_PAR = rad_sw * F * (1._wp - frac_UV)  ! This is Rpar_act; Eq. (E) from above
       rad_NIR = rad_sw - rad_PAR - rad_UV ! This is "Rnir"; Eq. (10) from above
       IF(fsw <= 0.2_wp) THEN
          frac_PAR_diffuse = 1._wp                                                     ! Eq. (12)-(14) from above
       ELSEIF(fsw >= 0.9_wp) THEN
          frac_PAR_diffuse = 0._wp                                                  ! Eq. (11)-(14) from above
       ELSE
          frac_PAR_diffuse = 1._wp - &
               (1._wp - ((0.9_wp - fsw)/0.7_wp)**(2._wp/3._wp))*hlp3_r/hlp4_r ! Eq. (11)-(14) from above
       ENDIF
    ELSE ! zenith angle larger than 89 degrees
       rad_UV = 0._wp
       rad_PAR = 0._wp
       rad_NIR = 0._wp
       frac_PAR_diffuse = 1._wp ! (could also be another value in this case)
       rad_sw  = 0._wp
    ENDIF
  END SUBROUTINE shortwave_from_fsw


  !>
  !! @brief computes PAR and NIR, together with the diffuse part of PAR, from fpar (according to ccdas but with some modifications)
  !!
  !! @par History
  !! JN: copied from jsbach3 mo_jsbalone_forcing
  !
  ! === shortwave_from_fsw_orig() ==================================================================================================
  !
  ! This routine is identical to the new shortwave routine from Wolfgang shortwave_from_fsw, except that the original conversion
  !
  ! factor
  ! (9)           1/F = 1 + (1.185-0.437*mfpar-0.494*mfpar^2)* exp[(0.0305-0.208*mfpar+0.495*mfpar^2)/cos(theta)]
  ! with mfpar=1-fpar is used (compare Eq. (9) in subroutine shortwave_from_fpar_orig()).
  ELEMENTAL SUBROUTINE shortwave_from_fsw_orig(day_of_year,cos_zenith,air_pressure,fsw, & !! inputs
      & rad_UV,rad_PAR,rad_NIR,frac_PAR_diffuse,rad_sw)  !! outputs

    INTEGER,INTENT(in)   :: day_of_year      !! day in year (from [1,365])
    REAL(wp),INTENT(in)  :: cos_zenith       !! cosine of zenith angle
    REAL(wp),INTENT(in)  :: air_pressure     !! air pressure at bottom (depending on local elevation) [N/m^2]
    REAL(wp),INTENT(in)  :: fsw              !! daily mean fraction of actual to potential TOTAL SHORTWAVE flux at bottom ..
                                             !! .. (from [0,1])
    REAL(wp),INTENT(out) :: rad_UV           !! "RUV_act": solar radiation flux (direct+diffuse) from the UV band 280-400 nm
                                             !! .. arriving actually at the bottom (includes cloud shading and scattering) [W/m^2]

    REAL(wp),INTENT(out) :: rad_PAR          !! "Rpar_act": solar radiation flux (direct+diffuse) from the visible band 400-700 nm .
                                             !! .. arriving actually at the bottom (includes cloud shading and scattering) [W/m^2]

    REAL(wp),INTENT(out) :: rad_NIR          !! "Rnir": flux (direct+diffuse) from the near infrared band 700-3000 nm ..
                                             !! .. actually arriving at the bottom (includes cloud shading and scattering) [W/m^2]
    REAL(wp),INTENT(out) :: frac_PAR_diffuse !! "fdiffPar": Fraction from  rad_PAR() that comes down as diffuse radiation ..
                                             !! .. (values in [0,1])
    REAL(wp),INTENT(out) :: rad_sw           !! "Rtot": total solar shortwave radiation : direct+diffuse, 280-3000 nm, ,..
                                             !! ... (i.e. PAR+NIR) that actually arrives at the bottom [W/m^2]

    REAL(wp) :: day_angle     !! day of year expressed as angle in radians
    REAL(wp) :: r_sun_earth_2 !! inverse squared earth-sun-distance normalized to mean distance
    REAL(wp) :: Rpar_top      !! PAR flux incident to the outer atmosphere at particular zenith angle [W/m^2]
    REAL(wp) :: F             !! factor F to obtain total shortwave from visible radiation (see Eq. (9'))
    REAL(wp) :: hlp_r,hlp2_r,hlp3_r,hlp4_r

    REAL(wp),PARAMETER :: F_pot = 0.43_wp !! See Eq. (C) from above
    REAL(wp),PARAMETER :: frac_UV = 0._wp !! UV fraction of (UV + PAR)
                                          !! to be checked, if the ratio rad_PAR/rad_sw is within observed range (0.42-0.51) with this value of frac_UV

    day_angle = 2._wp*pi*(day_of_year-1._wp)/365._wp !! ATTENTION: this may be incorrect for paleo situations
    hlp_r=2._wp*day_angle
    r_sun_earth_2 = &         ! Eq. (3) from above
      & 1.000110_wp+3.4221E-2_wp*COS(day_angle)+1.280E-3_wp*SIN(day_angle)+7.19E-4_wp*COS(hlp_r)+7.7E-5_wp*SIN(hlp_r)
    IF(cos_zenith > 0.017452406_wp) THEN! zenith angle smaller than 89 degrees  (cos(89) = 0.017452406)
       Rpar_top   = solar_const_vis*cos_zenith*r_sun_earth_2       ! Eq. (2) from above
       hlp3_r = EXP(-0.185_wp*(air_pressure/p0sl_bg)/cos_zenith) ! in jsb3: p0sl_bg = p_sealevel
       hlp4_r = 0.4_wp+0.6_wp*hlp3_r
       rad_sw = fsw*Rpar_top*hlp4_r/F_pot   ! This is "Rtot"; Eq. (D) from above
       hlp_r=1._wp-fsw
       hlp2_r=hlp_r*hlp_r
       F = 1._wp/(1._wp + (1.185_wp-0.437_wp*hlp_r-0.494_wp*hlp2_r) * &
         & EXP(MIN(2.534_wp,(0.0305_wp-0.208_wp*hlp_r+0.495_wp*hlp2_r)/cos_zenith))) ! Eq. (9) from above
       rad_UV = rad_sw * F * frac_UV
       rad_PAR = rad_sw * F * (1._wp - frac_UV)    ! This is Rpar_act; Eq. (E) from above
       rad_NIR = rad_sw - rad_PAR - rad_UV ! This is "Rnir"; Eq. (10) from above
       IF(fsw <= 0.2_wp) THEN
          frac_PAR_diffuse = 1._wp                                                     ! Eq. (12)-(14) from above
       ELSEIF(fsw >= 0.9_wp) THEN
          frac_PAR_diffuse = 0._wp                                                     ! Eq. (11)-(14) from above
       ELSE
          frac_PAR_diffuse = 1._wp - &
            & (1._wp - ((0.9_wp -fsw)/0.7_wp)**(2._wp/3._wp))*hlp3_r/hlp4_r ! Eq. (11)-(14) from above
       ENDIF
    ELSE ! zenith angle larger than 89 degrees
       rad_UV = 0._wp
       rad_PAR = 0._wp
       rad_NIR = 0._wp
       frac_PAR_diffuse = 1._wp ! (could also be another value in this case)
       rad_sw  = 0._wp
    ENDIF
  END SUBROUTINE shortwave_from_fsw_orig


  !>
  !! @brief computes shortwave parts (PAR,NIR,fraction of diffuse PAR) from the mean daily shortwave flux that arrives at the bottom
  !!
  !! @par History
  !! JN: copied from jsbach3 mo_jsbalone_forcing
  !
  ! === shortwave_from_shortwaveMean() =============================================================================================
  !
  ! This routine computes the various shortwave parts (PAR,NIR,fraction of diffuse PAR) from the mean daily shortwave flux that
  ! arrives at the bottom. The method is the same as that used in the subroutine "shortwave_from_cloudcover()", except that
  ! "fpar" (see description of  "shortwave_from_cloudcover()"), which is the fraction of actual to potential PAR at the earth's
  ! surface, is not computed from cloudcover, but estimated from the ratio of actual to potential mean daily shortwave flux
  !                       R0_act
  ! (1')           fpar = ------ ,
  !                       R0_pot
  ! where the actual mean daily shortwave flux "R0_act" is input to the routine and the potential mean daily shortwave flux
  ! is computed from the average shortwave flux reaching the bottom at a clear day:
  !
  !                                    1                                                             -0.185*p
  ! (A)           R0_pot = S * E *  ----- INTEGRAL[t_up:t_down] dt cos(theta(t)){ 0.4 + 0.6*exp(------------------) },
  !                                   24h                                                         p0*cos(theta(t))
  !
  ! where the term "0.4 + ..." accounts for Rayleigh scattering (which is insufficient!), "E" reflects the varying distance
  ! from the sun, "S=1360 W/m^2" is the solar constant and "t_up" and "t_down" are the times of sunrise and sunset
  ! (for other symbol names see subroutine "shortwave_from_fpar()").
  ! To compute the integral one needs the explicit representation of the cosine of the zenith-angle theta(t):
  ! (B)          cos(theta(t)) = sin(lat)*sin(dec) + cos(lat)*cos(dec)*cos(omega(t))
  ! where "lat" is the latitude, "dec" the declination and "omega(t)" is the hour-angle
  ! (C)          omega(t) = t*2*pi/24h ;
  ! here the time "t" is measured since/to the culmination of the sun at the particular day, i.e. in the morning "t" is positive
  ! and negative in the afternoon.
  ! The hour angles of sunrise and sunset are obtained from setting "theta(t) = pi/2" so that
  ! (D)          omega_sunrise = -omega_sunset = |arccos(tan(lat)*tan(dec))|.
  ! Since the cos(theta(t)) is symmetric the integral from (A) can be rewritten as
  !                                  1                                                       -0.185*p
  ! (A')          R0_pot = S * E *  -- INTEGRAL[0:omega_sunrise] dx cos(x){ 0.4 + 0.6*exp(--------------) },
  !                                 pi                                                       p0*cos(x)
  ! This value has to be provided when calling the subroutine (input variable "rad_sw_down_pot")
  !
  ELEMENTAL SUBROUTINE shortwave_from_shortwaveMean(day_of_year, cos_zenith, cos_lat, sin_lat, air_pressure,  &
    & rad_sw_down_act,rad_sw_down_pot, rad_UV_down, rad_PAR_down, rad_NIR_down, frac_PAR_diffuse, rad_sw_down, rad_sw_down_pot_out)

    INTEGER,INTENT(in)   :: day_of_year      !! day in year (from [1,365])
    REAL(wp),INTENT(in)  :: cos_zenith       !! cosine of zenith angle
    REAL(wp),INTENT(in)  :: cos_lat          !! cosine of latitude
    REAL(wp),INTENT(in)  :: sin_lat          !! sine of latitude
    REAL(wp),INTENT(in)  :: air_pressure     !! air pressure at bottom (depending on local elevation) [N/m^2]
    REAL(wp),INTENT(in)  :: rad_sw_down_act  !! daily mean of total solar shortwave radiation : direct+diffuse, 280-3000 nm, ..
                                             !! .. (i.e. PAR+NIR), that actually arrives at the bottom  [W/m^2]
    REAL(wp),INTENT(in)  :: rad_sw_down_pot  !! daily mean of total solar shortwave radiation : direct+diffuse, 280-3000 nm, ..
                                             !! ..(i.e. PAR+NIR), that potentially arrives at the bottom (i.e. at clear sky) [W/m^2]
    REAL(wp),INTENT(out) :: rad_UV_down      !! solar radiation flux (direct+diffuse) from the UV band 280-400 nm ...
                                             !! .. arriving at the bottom (includes cloud shading and scattering) [W/m^2]
    REAL(wp),INTENT(out) :: rad_PAR_down     !! solar radiation flux (direct+diffuse) from the visible band 400-700 nm ...
                                             !! .. arriving at the bottom (includes cloud shading and scattering) [W/m^2]
    REAL(wp),INTENT(out) :: rad_NIR_down     !! solar radiation flux (direct+diffuse) from the near infrared band 700-3000 nm ...
                                             !! .. arriving at the bottom (includes cloud shading and scattering) [W/m^2]
    REAL(wp),INTENT(out) :: frac_PAR_diffuse !! Fraction from  rad_PAR_down() that comes down as diffuse radiation ..
                                             !! .. (values in [0,1])
    REAL(wp),INTENT(out) :: rad_sw_down      !! total solar shortwave radiation flux: direct+diffuse, 280-3000 nm, (i.e. PAR+NIR),..
                                             !! ... that actually arrives at the bottom for the particular zenith angle [W/m^2]
    REAL(wp),INTENT(out) :: rad_sw_down_pot_out !! total solar shortwave radiation flux: direct+diffuse, 280-3000 nm, (i.e. PAR+NIR)
                                             !! ... that potentially arrives at the bottom for the particular zenith angle [W/m^2]

    REAL(wp),PARAMETER :: negligible_rad = 1e-8_wp !! Below this value (in [W/m^2]) radiation is assumed to be negligible
    REAL(wp) :: fsw         !! Fraction of actual out of potential radiation flux in the full shortwave band at surface

    IF(rad_sw_down_pot < negligible_rad ) THEN
       fsw = 0.5_wp
    ELSE
       fsw = MIN(rad_sw_down_act/rad_sw_down_pot,1._wp)
    ENDIF

    SELECT CASE(forcing_options%type_of_shortwave_scheme)
       CASE(SW_SCHEME_ORIGINAL_)
          CALL shortwave_from_fpar_orig(day_of_year,cos_zenith,air_pressure,fsw, &   !! inputs
            & rad_UV_down,rad_PAR_down,rad_NIR_down,frac_PAR_diffuse,rad_sw_down)   !! outputs
          rad_sw_down_pot_out = rad_sw_down_pot

       CASE(SW_SCHEME_CCDAS_)
          CALL shortwave_from_fpar_ccdas(day_of_year,cos_zenith,air_pressure,fsw, &   !! inputs
            & rad_UV_down,rad_PAR_down,rad_NIR_down,frac_PAR_diffuse,rad_sw_down)   !! outputs
          rad_sw_down_pot_out = rad_sw_down_pot

       CASE(SW_SCHEME_NEW_)
          CALL shortwave_from_fsw(day_of_year,cos_zenith,air_pressure,fsw, &   !! inputs
            & rad_UV_down,rad_PAR_down,rad_NIR_down,frac_PAR_diffuse,rad_sw_down)   !! outputs
          rad_sw_down_pot_out = rad_sw_down_pot

       CASE(SW_SCHEME_NEW_ORIG_CF_)
          CALL shortwave_from_fsw_orig(day_of_year,cos_zenith,air_pressure,fsw, &   !! inputs
            & rad_UV_down,rad_PAR_down,rad_NIR_down,frac_PAR_diffuse,rad_sw_down)   !! outputs
          rad_sw_down_pot_out = rad_sw_down_pot

       CASE(SW_SCHEME_DIRECT_)
          CALL shortwave_from_direct_shortwave(day_of_year,cos_zenith,cos_lat,sin_lat,air_pressure, &
            & rad_sw_down_act, rad_UV_down,rad_PAR_down,rad_NIR_down,frac_PAR_diffuse, rad_sw_down,rad_sw_down_pot_out)
   END SELECT
  END SUBROUTINE shortwave_from_shortwaveMean

  !>
  !! @brief computes longwave from observed longwave and diural course associated with air temperature
  !!
  !! @par History
  !! JN: copied from jsbach3 mo_jsbalone_forcing
  !
  ! === longwave_from_daily_longwave() ====================================================================================
  !
  ! Downward long wave radiation flux "R_d" [W/m^2] is according to [1],[2] computed by
  ! (1) R_d = r_cloud * epsA * sigma * T^4,
  ! where "T" is the air temperature [K] (in the bottom layer?), sigma = 5.6703e-8 W/(m^2 K^4) is the Stefan-Boltzmann constant,
  ! The mean daily longwave is known from observations, thus we are concerned only with the diurnal course associated with T
  ! [1] W. Knorr, "Satellite remote sensing and modelling of the global CO2 exchange of land vegetation", Examensarbeit 49,
  ! (Max Planck Institut fuer Meteorologie, Hamburg, 1998).
  ! [2] W. Brutsaert, "Evaporation into the Atmosphere", (Reidel, Dordrecht, 1982), pp. 299.
  !
  ELEMENTAL FUNCTION longwave_from_daily_longwave(longwave_down_daily, Tmin, Tmax, air_temp)
    REAL(wp)            :: longwave_from_daily_longwave !! hourly longwave downward radiation [W/m^2]
    REAL(wp),INTENT(in) :: longwave_down_daily          !! daily average downward longwave radiation [W/m^2]
    REAL(wp),INTENT(in) :: Tmin                         !! daily mimimum temperature [Celsius]
    REAL(wp),INTENT(in) :: Tmax                         !! daily maximum temperature [Celsius]
    REAL(wp),INTENT(in) :: air_temp                     !! air temperature at bottom [Celsius]

    REAL(wp) :: dtmp    ! daily average temperature
    REAL(wp) :: hlp_r   ! helper

    hlp_r = air_temp + tmelt                                ! Celsius --> Kelvin
    dtmp = 0.5_wp * (Tmax + Tmin) + tmelt
    longwave_from_daily_longwave = longwave_down_daily * (hlp_r**4) / (dtmp**4)      ! lw-radiation flux Eq. (1)

  END FUNCTION longwave_from_daily_longwave

  !>
  !! @brief computes longwave from cloud cover, vapor pressure and air temperature
  !!
  !! @par History
  !! JN: copied from jsbach3 mo_jsbalone_forcing
  !
  ! === longwave_rad_down() =======================================================================================================
  !
  ! Downward long wave radiation flux "R_d" [W/m^2] is according to [1],[2] computed by
  ! (1) R_d = r_cloud * epsA * sigma * T^4,
  ! where "T" is the air temperature [K] (in the bottom layer?), sigma = 5.6703e-8 W/(m^2 K^4) is the Stefan-Boltzmann constant,
  ! "epsA" is the emissivity of the cloudless atmosphere given by
  ! (2) epsA = epsA0 * (e_A/T)^(1/7).
  ! Here "e_vap" is the vapor pressure in [Pa] from above, and epsA0 = 0.64 (units: [(K/Pa)^(1/7)]) a constant (see [2]).
  ! The factor "r_cloud" in (1) corrects for clouds and is given by
  ! (3) r_cloud = 1+0.22*n^2,
  ! where "n" is the cloudcover.
  ! [1] W. Knorr, "Satellite remote sensing and modelling of the global CO2 exchange of land vegetation", Examensarbeit 49,
  ! (Max Planck Institut fuer Meteorologie, Hamburg, 1998).
  ! [2] W. Brutsaert, "Evaporation into the Atmosphere", (Reidel, Dordrecht, 1982), pp. 299.
  !
  ELEMENTAL FUNCTION longwave_rad_down(cloud_cover, vapor_pressure, air_temp)
    REAL(wp)            :: longwave_rad_down !! longwave downward radiation [W/m^2]
    REAL(wp),INTENT(in) :: cloud_cover       !! fraction of sky covered by clouds (from [0,1])
    REAL(wp),INTENT(in) :: vapor_pressure    !! vapor pressure [N/m^2]
    REAL(wp),INTENT(in) :: air_temp          !! air temperature at bottom [Celsius]

    REAL(wp) :: r_cloud ! factor correcting longwave radiation for clouds
    REAL(wp) :: epsA    ! emissivity of the cloudless atmosphere
    REAL(wp) :: hlp_r   ! helper

    r_cloud = 1._wp + 0.22_wp * cloud_cover * cloud_cover           ! cloud correction Eq. (3)
    hlp_r = air_temp + tmelt                                        ! Celsius --> Kelvin
    epsA = 0.64_wp * (vapor_pressure/hlp_r)**0.14285714_wp          ! emissivity Eq. (2)
    !jsbach3: stbo = StefanBoltzmann
    longwave_rad_down = r_cloud * epsA * stbo * hlp_r**4 ! lw-radiation flux Eq. (1)
  END FUNCTION longwave_rad_down


  !>
  !! @brief computes specific humidity
  !!
  !! @par History
  !! JN: copied from jsbach3 mo_jsbalone_forcing
  !
  ! === specific_humidity() =======================================================================================================
  !
  ! Specific humidity is computed here by assuming that each component of a mixture of dry air and vapor can be considered as an
  ! ideal gas. So, let m_d and m_v the masses of the two components in a fixed volume V. Then the specific humidity "q" is
  ! defined as
  ! (1) q = m_v/(m_v+m_d)
  ! Closely related is the moisture mixing ratio "r"
  ! (2) r = m_v/m_d
  ! Hence,
  ! (3) q = r/(1+r)
  ! Next it is assumed that both, dray air and vapor behave like ideal gases:
  ! (4) p_d*V = m_d*R_d*T and p_v*V = m_v*R_v*T,
  ! where p_d is the partial pressure of the dry air and p_v the vapor pressure, R_d and R_v the gas constants and T the
  ! temperature. Noting that air pressure "p" is the sum of the two partial pressures (p=p_d+p_v) the division of the two gas
  ! equation by each other gives by using (2):
  ! (5) r = eps*p_v/(p-p_v),
  ! where eps=R_d/R_v is the ratio of gas constants. Entering ths into (3) one finally obtains for the specific humidity:
  !                 p_v
  ! (6) q = eps -------------
  !             p-(1-eps)*p_v
  !
  ELEMENTAL FUNCTION specific_humidity(vapour_pressure, air_pressure)
    REAL(wp)            :: specific_humidity !! (from [0,1])
    REAL(wp),INTENT(IN) :: vapour_pressure   !! [N/m^2]
    REAL(wp),INTENT(IN) :: air_pressure      !! air pressure at bottom [N/m^2]

    !jsb3: rdv = eps
    specific_humidity = rdv * vapour_pressure/(air_pressure - (1._wp-rdv)*vapour_pressure)

  END FUNCTION specific_humidity

  !>
  !! @brief Inverse of specific_humidity
  !!
  !! @par History
  !! JN: copied from jsbach3 mo_jsbalone_forcing
  !
  ! === fun_vapour_pressure() =========================================================================================================
  !
  ! Inverse of specific_humidity()
  !
  ELEMENTAL FUNCTION fun_vapour_pressure(specific_humidity, air_pressure)
    REAL(wp) :: fun_vapour_pressure   !! [N/m^2]
    REAL(wp),INTENT(IN) :: specific_humidity !! (from [0,1])
    REAL(wp),INTENT(IN) :: air_pressure      !! air pressure at bottom [N/m^2]

    !jsb3: rdv = eps
    fun_vapour_pressure = air_pressure * specific_humidity / ( rdv + specific_humidity * ( 1._wp - rdv ) )

  END FUNCTION fun_vapour_pressure

#endif
END MODULE mo_jsb4_forcing
