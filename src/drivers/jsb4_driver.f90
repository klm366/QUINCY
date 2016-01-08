!>
!! @brief Contains the driver for jsbach standalone -> currently only and explicitly for use with Echam
!!        Initial version based on the jsbach3 driver 'jsbalone_driver'
!!
!! @author
!!  Julia Nabel, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Julia Nabel (2015-11-02)
!!
PROGRAM jsb4_driver

  USE mo_jsb_interface,      ONLY: jsbach_interface
  USE mo_jsb_base,           ONLY: init_base
  USE mo_jsb_model_init,     ONLY: init_model
  USE mo_jsb_version,        ONLY: jsbach_init_version, jsbach_label_run
  USE mo_jsb_namelist_iface, ONLY: POSITIONED, open_nml, position_nml
  USE mo_jsb_control,        ONLY: jsbach_runs_standalone
  USE mo_jsb_parallel_iface, ONLY: my_process_is_stdio, p_io, p_bcast, my_process_is_mpi_parallel
  USE mo_jsb_test,           ONLY: read_interface_vars, read_interface_variables

! Seems to be available in both ICON and ECHAM
  USE mo_exception,        ONLY: message, message_text
  USE mo_kind,             ONLY: wp

! TODO-check-with-RMS: echam!
  ! parallelisation
  USE mo_mpi,              ONLY: p_start, p_stop, p_init_communicators, p_global_comm
  USE mo_machine,          ONLY: machine_setup

  ! io
  USE mo_netcdf,           ONLY: nf_check, nf_close
  USE mo_memory_base,      ONLY: default_output
  USE mo_output,           ONLY: open_output_streams, close_output_streams, out_streams, init_output
  USE mo_decomposition,    ONLY: local_decomposition

  ! time control
  USE mo_timer,            ONLY: init_timer
  USE mo_time_control,     ONLY: lresume, time_set, time_reset, lstop, l_trigfiles, delta_time, &
                                    & time_step_len


  implicit none

  INTEGER,  EXTERNAL:: util_cputime
  REAL(wp), EXTERNAL:: util_walltime
  REAL(wp):: zwtime, zutime, zstime, zrtime

  INTEGER:: nml_handler, read_status, f_unit, status, jrow, nproma, ncid

  !! Variables of namelist jsbalone_parctl                       ! TODO-check-with-RMS: echam!
  INTEGER:: nproca  ! number of processors for jsbach
  INTEGER:: nprocb  ! has to be one
  INTEGER:: nprocio ! number of processors for I/O server
  INTEGER:: npedim  ! Working dimension for blocks in each domain.
                     ! Default (-1): each domain is processed in one call
  INTEGER:: buffer(4) ! buffer to broadcast the 4 above variables


  !! JSB4 interface variables - input
  REAL(wp), ALLOCATABLE :: t_air(:)
  REAL(wp), ALLOCATABLE :: q_air(:)
  REAL(wp), ALLOCATABLE :: rain(:)
  REAL(wp), ALLOCATABLE :: snow(:)
  REAL(wp), ALLOCATABLE :: wind_air(:)
  REAL(wp), ALLOCATABLE :: wind_10m(:)
  REAL(wp), ALLOCATABLE :: lw_srf_down(:)
  REAL(wp), ALLOCATABLE :: swvis_srf_down(:)
  REAL(wp), ALLOCATABLE :: swnir_srf_down(:)
  REAL(wp), ALLOCATABLE :: swpar_srf_down(:)
  REAL(wp), ALLOCATABLE :: press_srf(:)
  REAL(wp), ALLOCATABLE :: drag_srf(:)
  REAL(wp), ALLOCATABLE :: t_acoef(:)
  REAL(wp), ALLOCATABLE :: t_bcoef(:)
  REAL(wp), ALLOCATABLE :: q_acoef(:)
  REAL(wp), ALLOCATABLE :: q_bcoef(:)
  REAL(wp), ALLOCATABLE :: pch(:)
  REAL(wp), ALLOCATABLE :: cos_zenith_angle(:)

  !! JSB4 interface variables - output
  REAL(wp), ALLOCATABLE :: t_srf(:)
  REAL(wp), ALLOCATABLE :: t_eff_srf(:)
  REAL(wp), ALLOCATABLE :: qsat_srf(:)
  REAL(wp), ALLOCATABLE :: s_srf(:)
  REAL(wp), ALLOCATABLE :: fact_q_air(:)
  REAL(wp), ALLOCATABLE :: fact_qsat_srf(:)
  REAL(wp), ALLOCATABLE :: evapotrans(:)
  REAL(wp), ALLOCATABLE :: latent_hflx(:)
  REAL(wp), ALLOCATABLE :: sensible_hflx(:)
  REAL(wp), ALLOCATABLE :: grnd_hflx(:)
  REAL(wp), ALLOCATABLE :: grnd_hcap(:)
  REAL(wp), ALLOCATABLE :: rough_h_srf(:)
  REAL(wp), ALLOCATABLE :: rough_m_srf(:)
  REAL(wp), ALLOCATABLE :: tte_corr(:)
  REAL(wp), ALLOCATABLE :: alb_vis_dir(:)
  REAL(wp), ALLOCATABLE :: alb_nir_dir(:)
  REAL(wp), ALLOCATABLE :: alb_vis_dif(:)
  REAL(wp), ALLOCATABLE :: alb_nir_dif(:)


  INCLUDE 'jsbalone_parctl.inc'  ! TODO-check-with-RMS: echam!

  default_output = .TRUE.

  !=================================================================================================
  ! Start MPI
  CALL p_start

  !=================================================================================================
  CALL jsbach_init_version

  !=================================================================================================
  ! Initialize wallclock timer
!$OMP PARALLEL
!$OMP MASTER
  zwtime = util_walltime(0)
!$OMP END MASTER
!$OMP END PARALLEL

  !=================================================================================================
  ! read processor decomposition
  IF (my_process_is_stdio()) THEN

      ! define default values
      nproca = 1
      nprocb = 1
      nprocio = 0
      npedim = -1

      ! read the namelist
      nml_handler = open_nml ('namelist.jsbach')
      f_unit = position_nml ('JSBALONE_PARCTL', nml_handler, status=read_status)
      SELECT CASE (read_status)
      CASE (POSITIONED)
          READ (f_unit, jsbalone_parctl)
      END SELECT
  ENDIF

  IF (my_process_is_mpi_parallel()) THEN
      buffer(:) = (/ nproca, nprocb, nprocio, npedim /)
      CALL p_bcast (buffer, p_io, comm=p_global_comm)
      nproca  = buffer(1)
      nprocb  = buffer(2)
      nprocio = buffer(3)
      npedim  = buffer(4)
  END IF

  CALL p_init_communicators(nproca, nprocb, nprocio)
  CALL machine_setup

  CALL init_timer

  ! - jsbach_init
  CALL init_base('namelist.jsbach')
  CALL init_model(nproca, nprocb, npedim, 1)

  CALL jsbach_label_run(jsbach_runs_standalone())

  !TODO-check-with-RMS: include in init_model?
  CALL init_output 
  CALL open_output_streams

  IF (read_interface_vars) THEN
      ncid = 0
      nproma = local_decomposition%nproma

      ! Allocate fields -> JN: TODO -- nicht interface test fall
      ALLOCATE(t_air(nproma))
      ALLOCATE(q_air(nproma))
      ALLOCATE(rain(nproma))
      ALLOCATE(snow(nproma))
      ALLOCATE(wind_air(nproma))
      ALLOCATE(wind_10m(nproma))
      ALLOCATE(lw_srf_down(nproma))
      ALLOCATE(swvis_srf_down(nproma))
      ALLOCATE(swnir_srf_down(nproma))
      ALLOCATE(swpar_srf_down(nproma))
      ALLOCATE(press_srf(nproma))
      ALLOCATE(drag_srf(nproma))
      ALLOCATE(t_acoef(nproma))
      ALLOCATE(t_bcoef(nproma))
      ALLOCATE(q_acoef(nproma))
      ALLOCATE(q_bcoef(nproma))
      ALLOCATE(pch(nproma))
      ALLOCATE(cos_zenith_angle(nproma))

      ALLOCATE(t_srf(nproma))
      ALLOCATE(t_eff_srf(nproma))
      ALLOCATE(qsat_srf(nproma))
      ALLOCATE(s_srf(nproma))
      ALLOCATE(fact_q_air(nproma))
      ALLOCATE(fact_qsat_srf(nproma))
      ALLOCATE(evapotrans(nproma))
      ALLOCATE(latent_hflx(nproma))
      ALLOCATE(sensible_hflx(nproma))
      ALLOCATE(grnd_hflx(nproma))
      ALLOCATE(grnd_hcap(nproma))
      ALLOCATE(rough_h_srf(nproma))
      ALLOCATE(rough_m_srf(nproma))
      ALLOCATE(tte_corr(nproma))
      ALLOCATE(alb_vis_dir(nproma))
      ALLOCATE(alb_nir_dir(nproma))
      ALLOCATE(alb_vis_dif(nproma))
      ALLOCATE(alb_nir_dif(nproma))

  ELSE
      ! JN-TODO
      print *, '...to be implemented...'
      stop
  ENDIF


! - Main loop
  main_loop:DO 
      CALL time_set   !(echam!)

      IF (l_trigfiles) THEN
          CALL close_output_streams
          CALL open_output_streams
      ENDIF

  !     CALL stepon_jsbach      -->     sowas gibts in jsbach4 gar nicht!
  !        --> wurde in jsb3 am anfang jedes zeitschritts gerufen 
  !            und macht sowas wie      
  !                 - wenn ganz zu anfang, cover fracts, glacier masks, potential nat veg, veg_ratio_max, ... initialisieren 
  !                 - jedes neue Jahr: jährliche karten einlesen (e.g. landcover/transitions/nitrogen)

!    IF (.NOT. read_interface_vars) THEN
    ! CALL compute_orbit_and_solar(domain)        -> in jsbach3 "mo_zenith.f90" ?? kopieren ?? wie ersetzen ??

    ! TODO: read forcing date from file
    ! CALL read_forcing und CALL update_forcing   -->     mo_jsbalone_frocing! (TODO... erstmal kopieren?)
    ! CALL update_driving                         -->     mo_jsbalone          (TODO... erstmal kopieren?)
!    ENDIF

      DO jrow = 1, local_decomposition%ngpblks 

! JN-TODO: hm, in standalone case always nproma=npromz???
!          IF ( jrow == local_decomposition%ngpblks ) THEN
!            nproma = local_decomposition%npromz
!          ELSE
!            nproma = local_decomposition%nproma
!          END IF

          IF (read_interface_vars) THEN
            CALL read_interface_variables( 1, nproma, jrow, ncid, t_air, q_air, rain, snow, wind_air, wind_10m, &
                        & lw_srf_down, swvis_srf_down, swnir_srf_down, swpar_srf_down, press_srf, drag_srf, &
                        & t_acoef, t_bcoef, q_acoef, q_bcoef, pch, cos_zenith_angle )
          ENDIF

          CALL jsbach_interface ( 1, jrow, 1, nproma, delta_time, time_step_len,               &
                        & t_air, q_air, rain, snow, wind_air, wind_10m, lw_srf_down,           &
                        & swvis_srf_down, swnir_srf_down, swpar_srf_down, press_srf, drag_srf, &
                        & t_acoef, t_bcoef, q_acoef, q_bcoef, pch, cos_zenith_angle,           &
                        & t_srf, t_eff_srf, qsat_srf, s_srf, fact_q_air, fact_qsat_srf,        &
                        & evapotrans, latent_hflx, sensible_hflx, grnd_hflx, grnd_hcap,        & 
                        & rough_h_srf, rough_m_srf, tte_corr, alb_vis_dir, alb_nir_dir,        & 
                        & alb_vis_dif, alb_nir_dif)

      ENDDO

      CALL out_streams

    ! CALL stepoff_jsbach in jsb3 (macht hd! und clean) 
    ! --> sowas gibts hier nicht --> hm, aber wann und wie dann hd (geht im moment ja aber eh nicht mit echam) (TODO-later)

    !IF (l_putrerun) THEN  CALL write_streams  -> echam

      CALL time_reset
      IF (lstop) EXIT main_loop

  ENDDO main_loop

!!!!!!!!!!!!!!!
print *, '...jsb4_driver in development...'

!TODO-JN: after main loop there is a lot of further stuff in the jsb3 driver!

  IF (read_interface_vars) THEN
      CALL nf_check(nf_close(ncid))
  ENDIF

!!!!!!!!!!!!!!!

!$OMP PARALLEL
!$OMP MASTER
  status = util_cputime(zutime, zstime)
!$OMP END MASTER
!$OMP END PARALLEL
  IF (status == -1) THEN
      CALL message('jsb4_driver','Cannot determine used CPU time')
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

  ! Stop MPI
  CALL p_stop


END PROGRAM jsb4_driver
