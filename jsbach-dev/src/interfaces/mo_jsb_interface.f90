!>
!! @brief Interface to run JSBACH for one time step
!!
!! @author
!!  Reiner Schnur, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Reiner Schnur (2013-04-04)
!!
MODULE mo_jsb_interface

  USE mo_jsb_base,           ONLY: get_model
  USE mo_jsb_model_class,    ONLY: t_jsb_model
  USE mo_jsb_tile_class,     ONLY: t_jsb_tile
  USE mo_atmland_interface,  ONLY: update_atm2land, update_land2atm
  USE mo_jsb_parallel_iface, ONLY: get_my_global_mpi_id
  USE mo_jsb_control,        ONLY: debug_on
  USE mo_timer,              ONLY: timer_start, timer_stop, timer_jsbach
  USE mo_exception,          ONLY: message_text, message, finish
  USE mo_util_string,        ONLY: int2string
  USE mo_kind,               ONLY: wp

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: jsbach_interface, interface_full

  INTERFACE jsbach_interface
    MODULE PROCEDURE interface_full
    MODULE PROCEDURE interface_inquire
  END INTERFACE jsbach_interface

  CHARACTER(len=*), PARAMETER :: modname = 'mo_jsb_interface'

CONTAINS

  SUBROUTINE interface_full( model_id, iblk, ics, ice, dtime, steplen,                    &
                           & t_air, q_air, rain, snow, wind_air, wind_10m,                &
                           & lwrad_srf_down, swrad_srf_down, press_srf, drag_srf,         &
                           & t_acoef, t_bcoef, q_acoef, q_bcoef,                          &
                           & pch, cos_zenith_angle,                                       &
                           & t_srf, t_rad_srf, qsat_srf, fact_q_air, fact_qsat_srf,       &
                           & evapotrans, latent_hflx, sensible_hflx,                      &
                           & swrad_net,                                                   &
!!$                           & zh_srf, zm_srf,                                              &
                           & alb_vis_dir, alb_nir_dir, alb_vis_dif, alb_nir_dif)

    USE mo_veg_interface,     ONLY: update_vegetation_1, update_vegetation_2
    USE mo_srf_interface,     ONLY: update_surface_radiation, update_snow_and_skin_fraction, &
                                    update_surface_energy, update_surface_fluxes, update_t_srf, &
                                    update_humidity_scaling
    USE mo_soil_interface,    ONLY: update_soil_1, update_soil_2
    USE mo_soil_memory_class, ONLY: copy_soil_mem

    INTEGER :: p_pe

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':interface_full'

    INTEGER, INTENT(in) :: &
      & model_id,          &
      & iblk,              &
      & ics,               &
      & ice

    REAL(wp), INTENT(in) :: dtime, steplen
    REAL(wp), INTENT(in) ::   &
      & t_air            (:), &
      & q_air            (:), &
      & rain             (:), &
      & snow             (:), &
      & wind_air         (:), &
      & wind_10m         (:), &
      & lwrad_srf_down   (:), &
      & swrad_srf_down   (:), &
      & press_srf        (:), &
      & drag_srf         (:), &
      & t_acoef          (:), &
      & t_bcoef          (:), &
      & q_acoef          (:), &
      & q_bcoef          (:), &
      & pch              (:), &
      & cos_zenith_angle (:)
    REAL(wp), INTENT(out) ::  &
      & t_srf            (:), &
      & t_rad_srf        (:), &
      & qsat_srf         (:), &
      & fact_q_air       (:), &
      & fact_qsat_srf    (:), &
      & evapotrans       (:), &
      & latent_hflx      (:), &
      & sensible_hflx    (:), &
      & swrad_net        (:), &
!!$      & zh_srf           (:), &
!!$      & zm_srf           (:), &
      & alb_vis_dir      (:), &
      & alb_nir_dir      (:), &
      & alb_vis_dif      (:), &
      & alb_nir_dif      (:)

    TYPE(t_jsb_model), POINTER :: model

    ! Local variables
    !
    TYPE(t_jsb_tile), POINTER :: tile
    INTEGER :: nc, i

    CALL timer_start(timer_jsbach)

    nc = ice - ics + 1

    model => get_model(model_id)

!!$    WRITE(message_text, '(A,I4,A,I6,A,I6)') 'Updating '//TRIM(model%name)//' for: block=', iblk, ', ics=', ics, ', ice=', ice
!!$    CALL message(TRIM(routine), message_text, all_print=.TRUE.)
    p_pe = get_my_global_mpi_id()
!!$    CALL message( TRIM(routine), 'Updating '//TRIM(model%name)//' on PE '//int2string(p_pe)// &
!!$                               & ', nblk='//int2string(iblk), all_print=.TRUE.)
    IF (debug_on() .AND. iblk == 1) CALL message( TRIM(routine), 'Updating '//TRIM(model%name))

    ! Put atmospheric forcing into the atm2land memory.
    ! Only has to be done for the root tile since, for all other sub-tiles, the atm2land structure points to this.
    CALL update_atm2land( model%land, iblk, ics, ice,                                                                         &
!!$      & t_air=t_air(:), q_air=q_air(:), rain=rain(:), snow=snow(:), wind_air=wind_air(:), wind_10m=wind_10m(:),               &
      & t_air=t_air(:), q_air=q_air(:), rain=rain(:), snow=snow(:), wind_air=wind_air(:), wind_10m=wind_air(:)*0.8_wp,        & ! Temporary
      & lwrad_srf_down=lwrad_srf_down(:), swrad_srf_down=swrad_srf_down(:), press_srf=press_srf(:), drag_srf=drag_srf(:),     &
      & t_acoef=t_acoef(:), t_bcoef=t_bcoef(:), q_acoef=q_acoef(:), q_bcoef=q_bcoef(:),                                       &
      & pch=pch(:), cos_zenith_angle=cos_zenith_angle(:))

    ! Now we have to loop through the tiles for which the calculations have to be performed.

    ! For the moment, we only have one root tile
    !
    tile => model%land

    ! Compute net surface radiation
    CALL update_surface_radiation(model%land, iblk, ics, ice, nc)

    ! Update wet skin and snow fractions
    CALL  update_snow_and_skin_fraction(model%land, iblk, ics, ice, nc, dtime, steplen)

    ! Compute new surface temperatures and saturation surface specific humidity
    CALL update_surface_energy(model%land, iblk, ics, ice, nc, dtime, steplen)

    ! Compute surface fluxes of heat and moisture based on the new surface temperature and qsat
    CALL update_surface_fluxes(model%land, iblk, ics, ice, nc, steplen)

    ! First part of vegetation that is applied before the soil is updated,
    ! in particular:
    ! 1. Compute canopy conductance not limited by water stress
    ! 2. Compute soil moisture in root zone
    ! 3. Compute water stress factor
    ! 4. Compute canopy conductance limited by water stress
    ! 5. Compute skin reservoir and snow fractions on canopy
    !    For the moment, this is done in update_soil_1 to be compatible with echam6
    CALL update_vegetation_1(tile, iblk, ics, ice, nc)

    ! Soil hydrology
    CALL update_soil_1(tile, iblk, ics, ice, nc, dtime, steplen)

    ! Second part of vegetation that is applied after the soil has been updated,
    ! in particular:
    ! 1. Compute canopy albedo
    CALL update_vegetation_2(tile, iblk, ics, ice, nc)

    ! Second part of soil that is applied after the soil energy and water balance has been updated,
    ! in particular:
    ! 1. Compute ground albedo
    CALL update_soil_2(tile, iblk, ics, ice, nc)

    ! Apply time filter for surface temperature
    CALL update_t_srf(tile, iblk, ics, ice, nc)

    ! Re-calculate wet skin and snow fractions with updated moisture parameters
    CALL  update_snow_and_skin_fraction(model%land, iblk, ics, ice, nc, dtime, steplen)

    ! Update effective evaporation weights for qsat and q_air
    CALL update_humidity_scaling(tile, iblk, ics, ice, nc)

    ! Finally, give the required variables back to the atmosphere
    CALL update_land2atm(tile, iblk, ics, ice,                                                                                  &
      & t_srf=t_srf(:), t_rad_srf=t_rad_srf(:), qsat_srf=qsat_srf(:), fact_q_air=fact_q_air(:), fact_qsat_srf=fact_qsat_srf(:), &
      & evapotrans=evapotrans(:), latent_hflx=latent_hflx(:), sensible_hflx=sensible_hflx(:),                                   &
      & swrad_net=swrad_net(:),                                                                                                 &
!!$      & zh_srf=zh_srf(:), zm_srf=zm_srf(:),                                                                                     &
      & alb_vis_dir=alb_vis_dir(:), alb_nir_dir=alb_nir_dir(:), alb_vis_dif=alb_vis_dif(:), alb_nir_dif=alb_nir_dif(:))

!!$    IF (iblk == 1) CALL message(TRIM(routine), 'Updating of '//TRIM(model%name)//' on PE '//int2string(p_pe)//' completed.', &
!!$                                all_print=.TRUE.)
    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Updating of '//TRIM(model%name)//' completed.')

    CALL timer_stop(timer_jsbach)

  END SUBROUTINE interface_full

  SUBROUTINE interface_inquire( model_id, iblk, ics, ice, &
                              & t_srf, t_rad_srf, qsat_srf, zh_srf, zm_srf, alb_vis, alb_nir)

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':interface_inquire'

    INTEGER, INTENT(in) :: &
      & model_id,          &
      & iblk,              &
      & ics,               &
      & ice

    REAL(wp), INTENT(inout) :: &
      & t_srf    (:), &
      & t_rad_srf(:), &
      & qsat_srf (:), &
      & zh_srf   (:), &
      & zm_srf   (:), &
      & alb_vis  (:), &
      & alb_nir  (:)

  END SUBROUTINE interface_inquire

END MODULE mo_jsb_interface
