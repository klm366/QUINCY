!>
!! @brief Contains the driver for jsbach standalone -> currently only and explicitly for use with Echam
!!        Initial version based on the jsbach3 driver 'jsbalone_driver'
!!
!! @author
!!  Julia Nabel, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Julia Nabel (2015-11-02), copied several parts from jsbach3 jsbalone_driver.f90
!!
PROGRAM jsb4_driver
#ifndef __ICON__

  USE mo_jsb_interface,      ONLY: jsbach_interface
  USE mo_jsb_base,           ONLY: init_base
  USE mo_jsb_model_init,     ONLY: init_model
  USE mo_jsb_version,        ONLY: jsbach_init_version, jsbach_label_run
  USE mo_jsb_control,        ONLY: jsbach_runs_standalone, jsbach_is_restarted
  USE mo_jsb_time,           ONLY: is_newday, is_newyear
  USE mo_jsb_class,          ONLY: get_model
  USE mo_jsb_tile_class,     ONLY: t_jsb_tile_abstract
  USE mo_jsb_model_class,    ONLY: t_jsb_model

  USE mo_carbon_interface,   ONLY: check_carbon_conservation_for_operation, calculate_current_c_state_sum
  USE mo_carbon_constants,   ONLY: molarMassCO2_kg, sec_per_day

  dsl4jsb_Use_processes SEB_, TURB_, HYDRO_, CARBON_, FM_
  dsl4jsb_Use_memory(SEB_)
  dsl4jsb_Use_memory(TURB_)
  dsl4jsb_Use_memory(HYDRO_)
  dsl4jsb_Use_memory(FM_)

  USE mo_jsb4_forcing,       ONLY: get_interface_variables_from_external_forcing, init_forcing &
    & , get_data_and_update_field, input_variable_type, ANNUAL_

  ! Currently from adapter
  USE mo_jsb_domain_iface,   ONLY: get_lat, get_lon
  USE mo_jsb_namelist_iface, ONLY: POSITIONED, open_nml, position_nml
  USE mo_jsb_parallel,       ONLY: my_process_is_stdio, p_io, p_bcast, my_process_is_mpi_parallel, Get_omp_thread
  USE mo_jsb_test,           ONLY: read_interface_vars, read_interface_variables

  ! Seems to be available in both ICON and ECHAM
  USE mo_exception,        ONLY: message, message_text
  USE mo_kind,             ONLY: wp
  USE mo_math_constants,   ONLY: pi
  USE mo_orbit,            ONLY: orbit_kepler, orbit_vsop87

!>>> TODO? echam:
  ! parallelisation
  USE mo_mpi,              ONLY: p_start, p_stop, p_init_communicators
  USE mo_machine,          ONLY: machine_setup

  ! io
  USE mo_io,               ONLY: write_streams, IO_read_streams
  USE mo_filename,         ONLY: out_filetype
  USE mo_netcdf,           ONLY: nf_check, nf_close
  USE mo_memory_base,      ONLY: default_output, NETCDF
  USE mo_output,           ONLY: open_output_streams, close_output_streams, out_streams, init_output
  USE mo_decomposition,    ONLY: local_decomposition

  ! time control and orbit calc
  USE mo_timer,            ONLY: init_timer
  USE mo_time_control,     ONLY: time_set, time_reset, lstop, l_trigfiles, l_putrerun, delta_time, &
                                    & time_step_len, get_orbit_times, l_orbvsop87, current_date, write_date

  USE mo_radiation_parameters, ONLY: lyr_perp, nmonth, yr_perp, solar_parameters
  ! ICON would be:
  ! USE mo_psrad_radiation_parameters, ONLY: lyr_perp, nmonth, yr_perp, solar_parameters
!<<<

  IMPLICIT NONE

  CHARACTER(len=*), PARAMETER :: progname = 'jsb4_driver'

  INTEGER,  EXTERNAL:: util_cputime
  REAL(wp), EXTERNAL:: util_walltime
  REAL(wp):: zwtime, zutime, zstime, zrtime

  INTEGER:: nml_handler, read_status, f_unit, status, iblk, nproma, nblks, ncid, model_id, ice, ics, no_omp_thread

  !! Variables of namelist jsb_parallel_nml
  INTEGER:: nproca  ! number of processors for jsbach
  INTEGER:: nprocb  ! has to be one
  INTEGER:: nprocio ! number of processors for I/O server
  INTEGER:: npedim  ! Working dimension for blocks in each domain.
                     ! Default (-1): each domain is processed in one call
  INTEGER:: buffer(4) ! buffer to broadcast the 4 above variables

  LOGICAL :: l_rad_call = .FALSE.  ! TODO? echam -- orbit calculation
  REAL(wp) :: time_of_day, orbit_date, rasc_sun, decl_sun, dist_sun, flx_ratio
  REAL(wp), ALLOCATABLE :: lon(:,:), lat(:,:), coslon(:,:), coslat(:,:), sinlon(:,:), sinlat(:,:)
  REAL(wp), ALLOCATABLE :: daylght_frc(:,:)

  TYPE(t_jsb_model), POINTER :: model
  CLASS(t_jsb_tile_abstract),  POINTER :: tile

  dsl4jsb_Def_memory(SEB_)
  dsl4jsb_Def_memory(TURB_)
  dsl4jsb_Def_memory(HYDRO_)
  dsl4jsb_Def_memory(FM_)

  !! JSB4 interface variables - input
  REAL(wp), ALLOCATABLE :: t_air_K(:,:)
  REAL(wp), ALLOCATABLE :: q_air(:,:)
  REAL(wp), ALLOCATABLE :: rain(:,:)
  REAL(wp), ALLOCATABLE :: snow(:,:)
  REAL(wp), ALLOCATABLE :: wind_air(:,:)
  REAL(wp), ALLOCATABLE :: wind_10m(:,:)
  REAL(wp), ALLOCATABLE :: lw_srf_down(:,:)
  REAL(wp), ALLOCATABLE :: swvis_srf_down(:,:)
  REAL(wp), ALLOCATABLE :: swnir_srf_down(:,:)
  REAL(wp), ALLOCATABLE :: swpar_srf_down(:,:)
  REAL(wp), ALLOCATABLE :: press_srf(:,:)
  REAL(wp), ALLOCATABLE :: drag_srf(:,:)
  REAL(wp), ALLOCATABLE :: t_acoef(:,:)
  REAL(wp), ALLOCATABLE :: t_bcoef(:,:)
  REAL(wp), ALLOCATABLE :: q_acoef(:,:)
  REAL(wp), ALLOCATABLE :: q_bcoef(:,:)
  REAL(wp), ALLOCATABLE :: pch(:,:)
  REAL(wp), ALLOCATABLE :: cos_zenith_angle(:,:)

  !! TODO: so far not used in jsb4 but already read from forcing file/ calculated in mo_jsb4_forcing
  REAL(wp), ALLOCATABLE :: CO2_concentration(:,:)
  REAL(wp), ALLOCATABLE :: frac_par_diffuse(:,:)

  !! JSB4 interface variables - output, required for next timesteps forcing (if not reading interface vars)
  REAL(wp), ALLOCATABLE :: t_srf_proc(:,:)
  REAL(wp), ALLOCATABLE :: fact_q_air_proc(:,:)
  REAL(wp), ALLOCATABLE :: fact_qsat_srf_proc(:,:)
  REAL(wp), ALLOCATABLE :: rough_h_srf_proc(:,:)
  REAL(wp), ALLOCATABLE :: rough_m_srf_proc(:,:)
  REAL(wp), ALLOCATABLE :: evap_act2pot_proc(:,:)
!  REAL(wp), ALLOCATABLE :: evap_act_sum_proc(:,:)
!  REAL(wp), ALLOCATABLE :: evap_pot_sum_proc(:,:)


  !! JSB4 interface variables - output
  REAL(wp), ALLOCATABLE :: t_eff_srf(:)
  REAL(wp), ALLOCATABLE :: qsat_srf(:)
  REAL(wp), ALLOCATABLE :: s_srf(:)
  REAL(wp), ALLOCATABLE :: evapotrans(:)
  REAL(wp), ALLOCATABLE :: evapopot(:)

  REAL(wp), ALLOCATABLE :: latent_hflx(:)
  REAL(wp), ALLOCATABLE :: sensible_hflx(:)
  REAL(wp), ALLOCATABLE :: grnd_hflx(:)
  REAL(wp), ALLOCATABLE :: grnd_hcap(:)
  REAL(wp), ALLOCATABLE :: q_snocpymlt(:)
  REAL(wp), ALLOCATABLE :: alb_vis_dir(:)
  REAL(wp), ALLOCATABLE :: alb_nir_dir(:)
  REAL(wp), ALLOCATABLE :: alb_vis_dif(:)
  REAL(wp), ALLOCATABLE :: alb_nir_dif(:)

  REAL(wp), ALLOCATABLE :: CO2_flux(:,:)

  REAL(wp), ALLOCATABLE :: old_c_state_sum(:,:), cconservation(:,:)

  !! only required for FM - temporal solution 
  REAL(wp), ALLOCATABLE :: area_to_harvest(:,:)
  TYPE(input_variable_type) :: harvest_area_input


    NAMELIST /jsb_parallel_nml/     &
      nproca, npedim,             &
      nprocb, nprocio

  !-------------------------------------------------------------------------------------------------

  model_id = 1
  default_output = .TRUE.
  out_filetype = NETCDF

  ! Start MPI
  CALL p_start

  ! Initialize wallclock timer
!$OMP PARALLEL
!$OMP MASTER
  zwtime = util_walltime(0)
!$OMP END MASTER
!$OMP END PARALLEL

  ! read processor decomposition
  IF (my_process_is_stdio()) THEN
    ! define default values
    nproca = 1
    nprocb = 1
    nprocio = 0
    npedim = -1

    ! read the namelist
    nml_handler = open_nml ('namelist.jsbach')
    f_unit = position_nml ('jsb_parallel_nml', nml_handler, status=read_status)
    IF (read_status == POSITIONED) READ(f_unit, jsb_parallel_nml)
  ENDIF

  IF (my_process_is_mpi_parallel()) THEN
    buffer(:) = (/ nproca, nprocb, nprocio, npedim /)
    CALL p_bcast (buffer, p_io)
    nproca  = buffer(1)
    nprocb  = buffer(2)
    nprocio = buffer(3)
    npedim  = buffer(4)
  END IF

  CALL p_init_communicators(nproca, nprocb, nprocio)
  CALL machine_setup

  CALL init_timer

  ! - jsbach_init
  CALL jsbach_init_version

  CALL init_base('namelist.jsbach')
  CALL init_model(nproca, nprocb, npedim, model_id)

  CALL jsbach_label_run(jsbach_runs_standalone())

  no_omp_thread = Get_omp_thread()

  ! Allocate arrays -> current assumption is for standalone nlon=nproma=npromz,
  !     therefore set (and allocated fields) here -- else need to strongly adapt
  !       -> i.e. set and allocate in the iblk loop for just one and change way of reading?!
  nproma = local_decomposition%nproma
  nblks = local_decomposition%ngpblks
  !TODO: get nproma/nblks with already existing function? 'get_nproma'
  !  -> hm, for that need to first get the model and then the grid

  ! forcing
  ALLOCATE(t_air_K(nproma, nblks))
  ALLOCATE(q_air(nproma, nblks))
  ALLOCATE(rain(nproma, nblks))
  ALLOCATE(snow(nproma, nblks))
  ALLOCATE(wind_air(nproma, nblks))
  ALLOCATE(wind_10m(nproma, nblks))
  ALLOCATE(lw_srf_down(nproma, nblks))
  ALLOCATE(swvis_srf_down(nproma, nblks))
  ALLOCATE(swnir_srf_down(nproma, nblks))
  ALLOCATE(swpar_srf_down(nproma, nblks))
  ALLOCATE(press_srf(nproma, nblks))
  ALLOCATE(drag_srf(nproma, nblks))
  ALLOCATE(t_acoef(nproma, nblks))
  ALLOCATE(t_bcoef(nproma, nblks))
  ALLOCATE(q_acoef(nproma, nblks))
  ALLOCATE(q_bcoef(nproma, nblks))
  ALLOCATE(pch(nproma, nblks))
  ALLOCATE(cos_zenith_angle(nproma, nblks))

  ALLOCATE(CO2_concentration(nproma, nblks))
  ALLOCATE(frac_par_diffuse(nproma, nblks))

  ! output required for next timesteps forcing (if not reading interface vars)
  ALLOCATE(t_srf_proc(nproma, nblks))
  ALLOCATE(fact_q_air_proc(nproma, nblks))
  ALLOCATE(fact_qsat_srf_proc(nproma, nblks))
  ALLOCATE(rough_h_srf_proc(nproma, nblks))
  ALLOCATE(rough_m_srf_proc(nproma, nblks))
  ALLOCATE(evap_act2pot_proc(nproma, nblks))
!  ALLOCATE(evap_act_sum_proc(nproma, nblks))
!  ALLOCATE(evap_pot_sum_proc(nproma, nblks))

  ! output
  ALLOCATE(t_eff_srf(nproma))
  ALLOCATE(qsat_srf(nproma))
  ALLOCATE(s_srf(nproma))
  ALLOCATE(evapotrans(nproma))
  ALLOCATE(evapopot(nproma))
  ALLOCATE(latent_hflx(nproma))
  ALLOCATE(sensible_hflx(nproma))
  ALLOCATE(grnd_hflx(nproma))
  ALLOCATE(grnd_hcap(nproma))
  ALLOCATE(q_snocpymlt(nproma))
  ALLOCATE(alb_vis_dir(nproma))
  ALLOCATE(alb_nir_dir(nproma))
  ALLOCATE(alb_vis_dif(nproma))
  ALLOCATE(alb_nir_dif(nproma))

  ALLOCATE(CO2_flux(nproma, nblks))

  ! Only required for c conservation test
  ALLOCATE(old_c_state_sum(nproma, nblks))
  old_c_state_sum = 0._wp
  CO2_flux = 0._wp
  ALLOCATE(cconservation(nproma, nblks))

  IF (.NOT. read_interface_vars) THEN
    CALL init_forcing(model_id)

    ALLOCATE(lon(nproma,nblks))
    ALLOCATE(lat(nproma,nblks))

    !TODO: get from grid? but how to name such a function get_lon_from_grid?
    !      and if, shouldn't then also get_nproma should be called get_nproma_from_grid?
    lon = get_lon(decomp = local_decomposition)
    lat = get_lat(decomp = local_decomposition)

    ALLOCATE(coslon(nproma,nblks))
    ALLOCATE(coslat(nproma,nblks))
    ALLOCATE(sinlon(nproma,nblks))
    ALLOCATE(sinlat(nproma,nblks))
    ALLOCATE(daylght_frc(nproma, nblks))

    !init for first timestep/day (according to jsb3: mo_jsbalone: init_driving)
    t_srf_proc = 273._wp
    fact_q_air_proc = 0.5_wp
    fact_qsat_srf_proc = 0.5_wp
    rough_h_srf_proc = 1._wp
    rough_m_srf_proc = 1._wp
    evap_act2pot_proc = 1._wp

!    evap_act_sum_proc = 0._wp
!    evap_pot_sum_proc = 0._wp
  ENDIF

  model => get_model(model_id)
  tile => model%Get_top_tile()
  ics = 1
  ice = nproma

!TODO-JN: habe ich das getestet?
  IF(jsbach_is_restarted()) THEN

    CALL IO_read_streams

    dsl4jsb_Get_memory(SEB_)
    dsl4jsb_Get_memory(TURB_)
    dsl4jsb_Get_memory(HYDRO_)
    DO iblk = 1, nblks
      t_srf_proc(:,iblk) = dsl4jsb_var2D_onChunk(SEB_, t)
      fact_q_air_proc(:,iblk) = dsl4jsb_var2D_onChunk(TURB_, fact_q_air)
      fact_qsat_srf_proc(:,iblk) = dsl4jsb_var2D_onChunk(TURB_, fact_qsat_srf)
      rough_h_srf_proc(:,iblk) = dsl4jsb_var2D_onChunk(TURB_, rough_h)
      rough_m_srf_proc(:,iblk) = dsl4jsb_var2D_onChunk(TURB_, rough_m)
      ! TODO: if really required calculate daily sums of evapotrans and evapopot in jsbach and get them here (see comment below)
      evap_act2pot_proc(:,iblk) = MIN(dsl4jsb_var2D_onChunk(HYDRO_, evapotranspiration),-EPSILON(1._wp)) &
        & / MIN(dsl4jsb_var2D_onChunk(HYDRO_, evapopot),  -EPSILON(1._wp))
    ENDDO
  ENDIF

  IF(tile%Is_process_active(FM_)) THEN
    dsl4jsb_Get_memory(FM_)
    ALLOCATE(area_to_harvest(nproma, nblks))
    area_to_harvest(:,:) = -1._wp
  ENDIF


  CALL init_output
  CALL open_output_streams

  ! - Main loop
  main_loop:DO 

    CALL time_set   !(echam!)
    IF(is_newday()) THEN
      CALL write_date(current_date,'Calculate for: ')
    ENDIF

    IF (l_trigfiles) THEN
      CALL close_output_streams
      CALL open_output_streams
    ENDIF

    IF (.NOT. read_interface_vars) THEN

      ! copied from jsbach3: mo_zenith.f90
      CALL get_orbit_times(l_rad_call, lyr_perp, nmonth, yr_perp, time_of_day, orbit_date)
      IF (l_orbvsop87) THEN
         CALL orbit_vsop87 (orbit_date, rasc_sun, decl_sun, dist_sun)
      ELSE
         CALL orbit_kepler (orbit_date, rasc_sun, decl_sun, dist_sun)
      END IF

      coslon = COS(pi * lon / 180._wp)
      sinlon = SIN(pi * lon / 180._wp)
      coslat = COS(pi * lat / 180._wp)
      sinlat = SIN(pi * lat / 180._wp)

      CALL solar_parameters(decl_sun, dist_sun, time_of_day, sinlon, sinlat, coslon, coslat, &
        & flx_ratio, cos_zenith_angle, daylght_frc)

      CALL get_interface_variables_from_external_forcing(1, 1, nproma, nblks, sinlat, coslat, lon, evap_act2pot_proc, &
        & t_srf_proc, fact_q_air_proc, fact_qsat_srf_proc, rough_h_srf_proc, rough_m_srf_proc, cos_zenith_angle,      &
        & CO2_concentration, t_air_K, q_air, rain, snow, wind_air, wind_10m,                                          &
        & lw_srf_down, swvis_srf_down, swnir_srf_down, swpar_srf_down, frac_par_diffuse,                              &
        & press_srf, drag_srf, t_acoef, t_bcoef, q_acoef, q_bcoef, pch)

      IF(tile%Is_process_active(FM_)) THEN
        IF(is_newyear() .AND. my_process_is_stdio()) THEN
          ! for forest management test, before there is a better way available: read 2 arrays: 
          ! one with harvest area and one with harvest age 
          harvest_area_input%file_name_prefix = "Harvest_"
          harvest_area_input%expectedUnit = ""
          harvest_area_input%ocean_value = 0._wp
          harvest_area_input%frequency = ANNUAL_
          NULLIFY(harvest_area_input%timevalues)

          harvest_area_input%variable_name = "HARVEST_AREA"
          CALL get_data_and_update_field(nblks, ics, ice, model_id, harvest_area_input, area_to_harvest)
        ENDIF

        DO iblk = 1, nblks
          dsl4jsb_var2D_onChunk(FM_, area_to_harvest) = area_to_harvest(:,iblk)
        ENDDO
      ENDIF

    ENDIF

    DO iblk = 1, nblks
      CALL model%Set_options(iblk=iblk)

      IF(is_newday()) THEN

        IF (model%processes(CARBON_)%p%config%active) THEN
          CALL check_carbon_conservation_for_operation(tile, model%options(no_omp_thread), old_c_state_sum(:,iblk), &
            & - CO2_flux(:,iblk)/molarMassCO2_kg * sec_per_day, cconservation(:,iblk))

          IF(.NOT. jsbach_is_restarted()) THEN
            !TODO-JN: where to put this cconservation test? 
            !IF (iblk .EQ. 1) THEN
            IF (ANY(cconservation(:,iblk) .GT. 1.e-12_wp)) THEN
              print *, '----------------'
              print *, 'iblk: ', iblk
              print *, 'current sum: ', old_c_state_sum(:,iblk)
              print *, 'flux: ', - CO2_flux(:,iblk)/molarMassCO2_kg * sec_per_day
              print *, 'error: ', cconservation(:,iblk)
            ENDIF
                      CALL calculate_current_c_state_sum(tile, model%options(no_omp_thread), old_c_state_sum(:,iblk))
            !IF (iblk .EQ. 1) THEN
            IF (ANY(cconservation(:,iblk) .GT. 1.e-12_wp)) THEN
              print *, '***********************'
              print *, 'iblk:', iblk
              print *, 'old sum:', old_c_state_sum(:,iblk)
            ENDIF
            IF (ANY(cconservation(:,iblk) .GT. 1.e-11_wp)) THEN
              time_loop:DO 
                CALL time_set   !(echam!)
                IF(is_newday()) THEN
                  CALL write_date(current_date,'Calculate for: ')
                ENDIF
                CALL out_streams

                CALL time_reset
                IF(is_newyear()) STOP
              ENDDO time_loop
            STOP
            ENDIF
          ENDIF
        ENDIF
      ENDIF
      
! current assumption for standalone: nlon=nproma=npromz, therefore set (and allocated fields) outside the loop
!   Else would need large refactoring and:
!          IF ( iblk == nblks ) THEN
!            nproma = local_decomposition%npromz
!          ELSE
!            nproma = local_decomposition%nproma
!          END IF

      IF (read_interface_vars) THEN
        CALL read_interface_variables( 1, nproma, iblk, ncid, t_air_K(:,iblk), q_air(:,iblk), rain(:,iblk), &
          & snow(:,iblk), wind_air(:,iblk), wind_10m(:,iblk), lw_srf_down(:,iblk), swvis_srf_down(:,iblk), &
          & swnir_srf_down(:,iblk), swpar_srf_down(:,iblk), press_srf(:,iblk), drag_srf(:,iblk), &
          & t_acoef(:,iblk), t_bcoef(:,iblk), q_acoef(:,iblk), q_bcoef(:,iblk), pch(:,iblk), cos_zenith_angle(:,iblk))
      ENDIF

      CALL jsbach_interface ( 1, iblk, 1, nproma, delta_time, time_step_len, t_air_K(:,iblk), &
        & q_air(:,iblk), rain(:,iblk), snow(:,iblk), wind_air(:,iblk), wind_10m(:,iblk), lw_srf_down(:,iblk), &
        & swvis_srf_down(:,iblk), swnir_srf_down(:,iblk), swpar_srf_down(:,iblk), frac_par_diffuse(:,iblk), &
        & press_srf(:,iblk), drag_srf(:,iblk), &
        & t_acoef(:,iblk), t_bcoef(:,iblk), q_acoef(:,iblk), q_bcoef(:,iblk), pch(:,iblk), cos_zenith_angle(:,iblk), &
        & CO2_concentration(:,iblk), &
        & t_srf_proc(:,iblk), t_eff_srf(:), qsat_srf(:), s_srf(:), fact_q_air_proc(:,iblk), fact_qsat_srf_proc(:,iblk), &
        & evapotrans(:), latent_hflx(:), sensible_hflx(:), grnd_hflx(:), grnd_hcap(:), &
        & rough_h_srf_proc(:,iblk), rough_m_srf_proc(:,iblk), q_snocpymlt(:), alb_vis_dir(:), alb_nir_dir(:), &
        & alb_vis_dif(:), alb_nir_dif(:), CO2_flux(:,iblk), evapopot=evapopot(:))

      IF (.NOT. read_interface_vars) THEN
        ! TODO: if really required (only ERA-interim forcing data):
        ! To ensure restart identity, the daily evap_act_sum and evap_pot_sum would need to be calculated in jsbach!
        ! => for now the instantaneous values are used instead of the daily sums, to ensure restart identity
        evap_act2pot_proc(:,iblk) = MIN(evapotrans(:),-EPSILON(1._wp)) / MIN(evapopot(:),  -EPSILON(1._wp))
!        evap_act_sum_proc(:,iblk) = evap_act_sum(:,iblk) + MIN(evapotrans(:),-EPSILON(1._wp))
!        evap_pot_sum_proc(:,iblk) = evap_pot_sum(:,iblk) + MIN(evapopot(:),  -EPSILON(1._wp))
      ENDIF

    ENDDO

!    IF ((.NOT. read_interface_vars) .AND. is_newday()) THEN
!      evap_act2pot_proc = evap_act_sum_proc / evap_pot_sum_proc
!      evap_act_sum_proc = 0._wp
!      evap_pot_sum_proc = 0._wp
!    ENDIF

    CALL out_streams

    IF (l_putrerun) CALL write_streams

    CALL time_reset
    IF (lstop) EXIT main_loop

  ENDDO main_loop

  IF (read_interface_vars) THEN
    CALL nf_check(nf_close(ncid))
  ENDIF

!$OMP PARALLEL
!$OMP MASTER
  status = util_cputime(zutime, zstime)
!$OMP END MASTER
!$OMP END PARALLEL
  IF (status == -1) THEN
    CALL message(TRIM(progname),'Cannot determine used CPU time')
  ELSE
!$OMP PARALLEL
!$OMP MASTER
    zwtime = util_walltime(0)
!$OMP END MASTER
!$OMP END PARALLEL
    zrtime = (zutime+zstime)/zwtime
    CALL message ('', '')
    WRITE (message_text,'(a,f10.2,a)') ' Wallclock        : ', zwtime, ' s'
    CALL message('',message_text)
    WRITE (message_text,'(a,f10.2,a)') ' CPU-time (user)  : ', zutime, ' s'
    CALL message('',message_text)
    WRITE (message_text,'(a,f10.2,a)') ' CPU-time (system): ', zstime, ' s'
    CALL message('',message_text)
    WRITE (message_text,'(a,f10.2,a)') ' Ratio            : ', 100*zrtime, ' %'
    CALL message('',message_text)
    CALL message ('', '')
  END IF

  ! deallocate in- and output arrays
  DEALLOCATE(t_air_K, q_air, rain, snow, wind_air, wind_10m)
  DEALLOCATE(lw_srf_down, swvis_srf_down, swnir_srf_down, swpar_srf_down, frac_par_diffuse)
  DEALLOCATE(press_srf, drag_srf, t_acoef, t_bcoef, q_acoef, q_bcoef, pch, cos_zenith_angle, CO2_concentration)

  DEALLOCATE(t_eff_srf, qsat_srf, s_srf)
  DEALLOCATE(evapotrans, latent_hflx, sensible_hflx, grnd_hflx, grnd_hcap)
  DEALLOCATE(q_snocpymlt, alb_vis_dir, alb_nir_dir, alb_vis_dif, alb_nir_dif)
  DEALLOCATE(CO2_flux)

  DEALLOCATE(old_c_state_sum)
  DEALLOCATE(cconservation)

  IF (.NOT. read_interface_vars) THEN
    DEALLOCATE(lon, lat, coslon, coslat, sinlon, sinlat, daylght_frc)
    DEALLOCATE(t_srf_proc, fact_q_air_proc, fact_qsat_srf_proc, rough_h_srf_proc, rough_m_srf_proc, evap_act2pot_proc)
    !DEALLOCATE(evap_act_sum_proc, evap_pot_sum_proc)
  ENDIF

  ! Stop MPI
  CALL p_stop
#endif
END PROGRAM jsb4_driver
