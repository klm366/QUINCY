!>
!! @brief Interface to run JSBACH for one time step
!!
!! @author
!!  Rackhun Son, MPI-BGC Jena
!!
!! @par Revision History
!! First version by Rackhun Son (2022-12-29)
!!
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/license .
!!

MODULE mo_anthro2l_interface
#ifndef __NO_JSBACH__

  USE mo_exception,       ONLY: finish
  USE mo_kind,            ONLY: wp

  USE mo_jsb_control,        ONLY: jsbach_runs_standalone
  USE mo_jsb_model_class,    ONLY: t_jsb_model
  USE mo_jsb_class,          ONLY: Get_model
  USE mo_jsb_tile_class,     ONLY: t_jsb_tile_abstract
  USE mo_jsb_lct_class,      ONLY: LAKE_TYPE
  !USE mo_jsb_config_class,   ONLY: t_jsb_config
  !USE mo_jsb_process_class,  ONLY: t_jsb_process
  USE mo_jsb_task_class,     ONLY: t_jsb_task_options

  dsl4jsb_Use_processes ANTHRO2L_ 
  !dsl4jsb_Use_config(SEB_)
  dsl4jsb_Use_memory(ANTHRO2L_)

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: update_anthro2l

  CHARACTER(len=*), PARAMETER :: modname = 'mo_anthro2l_interface'


CONTAINS

  ! ======================================================================================================= !
  !> update anthropogenic2land
  !>
  SUBROUTINE update_anthro2l( &
    & tile, options,                    &
    & a_temperature,                    &
    & a_rel_humidity,                   &
    & fba,                              &
    & population,                       &
    & gdp,                              &
    & hdi,                              &
    & land_use_primf,                   &
    & land_use_c4per,                   &
    & land_use_urban,                   &
    & land_use_secdf,                   &
    & land_use_secma,                   &
    & land_use_primn,                   &
    & land_use_secdn,                   &
    & land_use_c4ann,                   &
    & land_use_pastr,                   &
    & land_use_range,                   &
    & land_use_c3ann,                   &
    & land_use_c3nfx,                   &
    & land_use_c3per,                   &
    & land_use_secmb,                   &
    & DEBUG_VAR                         &
    & )

    CLASS(t_jsb_tile_abstract), INTENT(inout)             :: tile
    TYPE(t_jsb_task_options),   INTENT(in)                :: options

    REAL(wp), OPTIONAL, DIMENSION(:), INTENT(in) ::  &
      & DEBUG_VAR,                                   &
      & a_temperature,                               &
      & a_rel_humidity,                              &
      & fba,                                         &
      & population,                                  &
      & gdp,                                         &
      & hdi,                                         &
      & land_use_primf,                              &
      & land_use_c4per,                              &
      & land_use_urban,                              &
      & land_use_secdf,                              &
      & land_use_secma,                              &
      & land_use_primn,                              &
      & land_use_secdn,                              &
      & land_use_c4ann,                              &
      & land_use_pastr,                              &
      & land_use_range,                              &
      & land_use_c3ann,                              &
      & land_use_c3nfx,                              &
      & land_use_c3per,                              &
      & land_use_secmb                   
    
    dsl4jsb_Def_memory(ANTHRO2L_)

    INTEGER ::  iblk, ics, ice, nc, i
    TYPE(t_jsb_model), POINTER :: model

    CHARACTER(len=*), PARAMETER :: routine = modname//':update_anthro2l'

    dsl4jsb_Real2D_onChunk ::  &
      & DEBUG_VAR_ptr,         &
      & a_temperature_ptr,     &
      & a_rel_humidity_ptr,    &
      & fba_ptr,               &
      & population_ptr,        &
      & gdp_ptr,               &
      & hdi_ptr,               &
      & land_use_primf_ptr,    &
      & land_use_c4per_ptr,    &
      & land_use_urban_ptr,    &
      & land_use_secdf_ptr,    &
      & land_use_secma_ptr,    &
      & land_use_primn_ptr,    &
      & land_use_secdn_ptr,    &
      & land_use_c4ann_ptr,    &
      & land_use_pastr_ptr,    &
      & land_use_range_ptr,    &
      & land_use_c3ann_ptr,    &
      & land_use_c3nfx_ptr,    &
      & land_use_c3per_ptr,    &
      & land_use_secmb_ptr

    iblk = options%iblk
    ics  = options%ics
    ice  = options%ice
    nc   = options%nc

    model => Get_model(tile%owner_model_id)

    IF (nc /= SIZE(a_rel_humidity,1)) CALL finish(TRIM(routine), 'Wrong dimensions')

    dsl4jsb_Get_memory(ANTHRO2L_)
    IF (PRESENT(DEBUG_VAR)) THEN
      DEBUG_VAR_ptr         => dsl4jsb_var2D_onChunk(ANTHRO2L_, DEBUG_VAR)
    END IF

    IF (PRESENT(a_temperature)) THEN
      a_temperature_ptr      => dsl4jsb_var2D_onChunk(ANTHRO2L_, a_temperature)
    END IF
    IF (PRESENT(a_rel_humidity)) THEN
      a_rel_humidity_ptr      => dsl4jsb_var2D_onChunk(ANTHRO2L_, a_rel_humidity)
    END IF
    IF (PRESENT(fba)) THEN
      fba_ptr      => anthro2l__mem%fba%ptr(ics:ice,iblk)
    END IF
    IF (PRESENT(population)) THEN
      population_ptr          => dsl4jsb_var2D_onChunk(ANTHRO2L_, population)
    END IF
    IF (PRESENT(gdp)) THEN
      gdp_ptr                 => dsl4jsb_var2D_onChunk(ANTHRO2L_, gdp)
    END IF
    IF (PRESENT(hdi)) THEN
      hdi_ptr                 => dsl4jsb_var2D_onChunk(ANTHRO2L_, hdi)
    END IF
    IF (PRESENT(land_use_primf)) THEN
      land_use_primf_ptr      => dsl4jsb_var2D_onChunk(ANTHRO2L_, land_use_primf)
    END IF
    IF (PRESENT(land_use_c4per)) THEN
      land_use_c4per_ptr      => dsl4jsb_var2D_onChunk(ANTHRO2L_, land_use_c4per)
    END IF
    IF (PRESENT(land_use_urban)) THEN
      land_use_urban_ptr      => dsl4jsb_var2D_onChunk(ANTHRO2L_, land_use_urban)
    END IF
    IF (PRESENT(land_use_secdf)) THEN
      land_use_secdf_ptr      => dsl4jsb_var2D_onChunk(ANTHRO2L_, land_use_secdf)
    END IF
    IF (PRESENT(land_use_secma)) THEN
      land_use_secma_ptr      => dsl4jsb_var2D_onChunk(ANTHRO2L_, land_use_secma)
    END IF
    IF (PRESENT(land_use_primn)) THEN
      land_use_primn_ptr      => dsl4jsb_var2D_onChunk(ANTHRO2L_, land_use_primn)
    END IF
    IF (PRESENT(land_use_secdn)) THEN
      land_use_secdn_ptr      => dsl4jsb_var2D_onChunk(ANTHRO2L_, land_use_secdn)
    END IF
    IF (PRESENT(land_use_c4ann)) THEN
      land_use_c4ann_ptr      => dsl4jsb_var2D_onChunk(ANTHRO2L_, land_use_c4ann)
    END IF
    IF (PRESENT(land_use_pastr)) THEN
      land_use_pastr_ptr      => dsl4jsb_var2D_onChunk(ANTHRO2L_, land_use_pastr)
    END IF
    IF (PRESENT(land_use_range)) THEN
      land_use_range_ptr      => dsl4jsb_var2D_onChunk(ANTHRO2L_, land_use_range)
    END IF
    IF (PRESENT(land_use_c3ann)) THEN
      land_use_c3ann_ptr      => dsl4jsb_var2D_onChunk(ANTHRO2L_, land_use_c3ann)
    END IF
    IF (PRESENT(land_use_c3nfx)) THEN
      land_use_c3nfx_ptr      => dsl4jsb_var2D_onChunk(ANTHRO2L_, land_use_c3nfx)
    END IF
    IF (PRESENT(land_use_c3per)) THEN
      land_use_c3per_ptr      => dsl4jsb_var2D_onChunk(ANTHRO2L_, land_use_c3per)
    END IF
    IF (PRESENT(land_use_secmb)) THEN
      land_use_secmb_ptr      => dsl4jsb_var2D_onChunk(ANTHRO2L_, land_use_secmb)
    END IF


    !$acc parallel default(present)

    IF (PRESENT(DEBUG_VAR)) THEN
      !$acc loop gang vector
      DO i=1,nc
        DEBUG_VAR_ptr(i) = DEBUG_VAR(i)
      END DO
    END IF

    IF (PRESENT(a_temperature)) THEN
      !$acc loop gang vector
      DO i=1,nc
        a_temperature_ptr(i) = a_temperature(i)
      END DO
    END IF

    IF (PRESENT(a_rel_humidity)) THEN
      !$acc loop gang vector
      DO i=1,nc
        a_rel_humidity_ptr(i) = a_rel_humidity(i)
      END DO
    END IF

    IF (PRESENT(fba)) THEN
      !$acc loop gang vector
      DO i=1,nc
        fba_ptr(i) = fba(i)
      END DO
    END IF

    IF (PRESENT(population)) THEN
      !$acc loop gang vector
      DO i=1,nc
        population_ptr(i) = population(i)
      END DO
    END IF

    IF (PRESENT(gdp)) THEN
      !$acc loop gang vector
      DO i=1,nc
        gdp_ptr(i) = gdp(i)
      END DO
    END IF

    IF (PRESENT(hdi)) THEN
      !$acc loop gang vector
      DO i=1,nc
        hdi_ptr(i) = hdi(i)
      END DO
    END IF

    IF (PRESENT(land_use_primf)) THEN
      !$acc loop gang vector
      DO i=1,nc
        land_use_primf_ptr(i) = land_use_primf(i)
      END DO
    END IF

    IF (PRESENT(land_use_c4per)) THEN
      !$acc loop gang vector
      DO i=1,nc
        land_use_c4per_ptr(i) = land_use_c4per(i)
      END DO
    END IF

    IF (PRESENT(land_use_urban)) THEN
      !$acc loop gang vector
      DO i=1,nc
        land_use_urban_ptr(i) = land_use_urban(i)
      END DO
    END IF

    IF (PRESENT(land_use_secdf)) THEN
      !$acc loop gang vector
      DO i=1,nc
        land_use_secdf_ptr(i) = land_use_secdf(i)
      END DO
    END IF

    IF (PRESENT(land_use_secma)) THEN
      !$acc loop gang vector
      DO i=1,nc
        land_use_secma_ptr(i) = land_use_secma(i)
      END DO
    END IF

    IF (PRESENT(land_use_primn)) THEN
      !$acc loop gang vector
      DO i=1,nc
        land_use_primn_ptr(i) = land_use_primn(i)
      END DO
    END IF

    IF (PRESENT(land_use_secdn)) THEN
      !$acc loop gang vector
      DO i=1,nc
        land_use_secdn_ptr(i) = land_use_secdn(i)
      END DO
    END IF

    IF (PRESENT(land_use_c4ann)) THEN
      !$acc loop gang vector
      DO i=1,nc
        land_use_c4ann_ptr(i) = land_use_c4ann(i)
      END DO
    END IF

    IF (PRESENT(land_use_pastr)) THEN
      !$acc loop gang vector
      DO i=1,nc
        land_use_pastr_ptr(i) = land_use_pastr(i)
      END DO
    END IF

    IF (PRESENT(land_use_range)) THEN
      !$acc loop gang vector
      DO i=1,nc
        land_use_range_ptr(i) = land_use_range(i)
      END DO
    END IF

    IF (PRESENT(land_use_c3ann)) THEN
      !$acc loop gang vector
      DO i=1,nc
        land_use_c3ann_ptr(i) = land_use_c3ann(i)
      END DO
    END IF

    IF (PRESENT(land_use_c3nfx)) THEN
      !$acc loop gang vector
      DO i=1,nc
        land_use_c3nfx_ptr(i) = land_use_c3nfx(i)
      END DO
    END IF

    IF (PRESENT(land_use_c3per)) THEN
      !$acc loop gang vector
      DO i=1,nc
        land_use_c3per_ptr(i) = land_use_c3per(i)
      END DO
    END IF

    IF (PRESENT(land_use_secmb)) THEN
      !$acc loop gang vector
      DO i=1,nc
        land_use_secmb_ptr(i) = land_use_secmb(i)
      END DO
    END IF
    !$acc end parallel

  END SUBROUTINE update_anthro2l

#endif
END MODULE mo_anthro2l_interface
