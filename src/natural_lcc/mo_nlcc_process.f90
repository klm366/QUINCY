!> Contains the routines for the nlcc processes
!>
!> ICON-Land
!>
!> ---------------------------------------
!> Copyright (C) 2013-2026, MPI-M, MPI-BGC
!>
!> Contact: icon-model.org
!> Authors: AUTHORS.md
!> See LICENSES/ for license information
!> SPDX-License-Identifier: BSD-3-Clause
!> ---------------------------------------
!>
MODULE mo_nlcc_process
#ifndef __NO_JSBACH__

  USE mo_jsb_control, ONLY: debug_on
  USE mo_kind,        ONLY: wp
  USE mo_exception,   ONLY: message, finish, message_text
  USE mo_jsb_physical_constants, ONLY: tmelt

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: calc_climbuf, bioclim_limits, scale_fpc, potential_tree_fpc, calc_cover_fract_inst
  PUBLIC :: fpc_to_cover_fract_pot, cover_fract_pot_to_cover_fract, desert_fraction, fpc_daily, fract_small

  CHARACTER(len=*), PARAMETER :: modname = 'mo_nlcc_process'
  REAL(wp), PARAMETER :: fract_small = 1.e-10_wp   !! smallest fraction of PFT tiles
  REAL(wp), PARAMETER :: act_fpc_min = 0.005_wp    !! minimum actual FPC to be considered to calculate potential FPC
  REAL(wp), PARAMETER :: sum_npp_min = 1.e-12_wp   !! minimum total npp for tree cover
  REAL(wp), PARAMETER :: tree_fpc_max = 1.0_wp     !! maximum FPC of woody PFTs
  REAL(wp), PARAMETER :: npp_nonlinearity = 1.5_wp !! parameter controlling the non-linearity in the dynamic equation with
                                                   !! respect to NPP
  REAL(wp), PARAMETER :: desert_extend = 0.65_wp   !! parameter controlling the extend of desert
                                                   !! (the lower the value the more desert)
  REAL(wp), PARAMETER :: desert_margin = 2._wp     !! parameter controlling the transition from vegetated land to desert
                                                   !! (the higher the value the sharper the transition)
  REAL(wp), PARAMETER :: tau_desert = 50._wp       !! time constant [year] by which the extension of deserts is adapted

CONTAINS
  !
  ! ===============================================================================================================================
  !>
  !! Compute running averages of climate variables for nlcc
  !!

#ifndef _OPENACC
  ELEMENTAL PURE &
#endif
  SUBROUTINE calc_climbuf (lstart, init_running_means, new_day, new_month, new_year, delta_time, t_air, &
    &                      seconds_day, seconds_month, temp_sum_day, temp_sum_month, &
    &                      min_mmtemp_of_yr, max_mmtemp_of_yr, min_mmtemp20, max_mmtemp20, &
    &                      gdd_sum_year, gdd_prev_year)

    !$ACC ROUTINE SEQ

    LOGICAL, INTENT(in) :: lstart ! first time step in a simulation
    LOGICAL, INTENT(in) :: init_running_means ! initialize the min/max temperature climatology
    LOGICAL, INTENT(in) :: new_day ! first time step in a day
    LOGICAL, INTENT(in) :: new_month ! firest time step in a month
    LOGICAL, INTENT(in) :: new_year ! first time step in a year
    REAL(wp), INTENT(in) :: delta_time ! time step length [s]
    REAL(wp), INTENT(in) :: seconds_day ! lenght of a day [s]
    REAL(wp), INTENT(in) :: seconds_month ! length of the previous month [s]
    REAL(wp), INTENT(in) :: t_air ! air temperature of lowest atmosphere level [K]
    REAL(wp), INTENT(inout) :: temp_sum_day
    REAL(wp), INTENT(inout) :: temp_sum_month
    REAL(wp), INTENT(inout) :: min_mmtemp_of_yr ! temperature of coldest month of previous year [K]
    REAL(wp), INTENT(inout) :: max_mmtemp_of_yr ! temperature of warmest month of previous year [K]
    REAL(wp), INTENT(inout) :: min_mmtemp20 ! temperature of coldest month climatology [K]
    REAL(wp), INTENT(inout) :: max_mmtemp20 ! temperature of warmest month climatology [K]
    REAL(wp), INTENT(inout) :: gdd_sum_year ! growing degree days of current year [C]
    REAL(wp), INTENT(inout) :: gdd_prev_year ! growing degree days of previous year [C]

!   local variables
    REAL(wp) :: prev_day_mean_temp   ! air temperature of previous day [K]
    REAL(wp) :: prev_month_mean_temp ! air temperature of previous month [K]

    IF (new_year) THEN
       gdd_prev_year = gdd_sum_year
       gdd_sum_year = 0._wp
    ENDIF

    IF (.NOT. new_day .OR. lstart) THEN
       temp_sum_day = temp_sum_day + t_air * delta_time
    ELSE
       prev_day_mean_temp = temp_sum_day / seconds_day
       temp_sum_day = t_air * delta_time
       gdd_sum_year = gdd_sum_year + MAX(0._wp,prev_day_mean_temp - tmelt - 5._wp)
    ENDIF

    IF (.NOT. new_month .OR. lstart) THEN
       temp_sum_month = temp_sum_month + t_air * delta_time
    ELSE
       prev_month_mean_temp = temp_sum_month / seconds_month
       temp_sum_month = t_air * delta_time
       ! --- minimum/maximum monthly mean temperature (20 year climatology)
       ! --- calculate coldest/warmest month of the year so far
       min_mmtemp_of_yr = min(min_mmtemp_of_yr, prev_month_mean_temp)
       max_mmtemp_of_yr = max(max_mmtemp_of_yr, prev_month_mean_temp)
       IF (new_year) THEN
          ! --- build 20yr running mean of coldest and warmest monthly temperature
          IF (init_running_means) THEN
             min_mmtemp20 = min_mmtemp_of_yr
             max_mmtemp20 = max_mmtemp_of_yr
          ELSE
             min_mmtemp20 = (min_mmtemp20 * 19._wp + min_mmtemp_of_yr) / 20._wp
             max_mmtemp20 = (max_mmtemp20 * 19._wp + max_mmtemp_of_yr) / 20._wp
          ENDIF
          ! --- reset for current year
          min_mmtemp_of_yr =  1000._wp
          max_mmtemp_of_yr = -1000._wp
       END IF ! new_year
    END IF ! new_month

  END SUBROUTINE calc_climbuf

! --- bioclim_limits ---------------------------------------------------------------------------------------

#ifndef _OPENACC
  ELEMENTAL PURE &
#endif
  SUBROUTINE bioclim_limits (is_dynamic, lct_tcmin, lct_tcmax, lct_twmax, lct_min_temprange, lct_gddmin, &
    &                        min_mmtemp20, max_mmtemp20, gdd_prev_year, bio_exist)

!$ACC ROUTINE SEQ

! calculation of bioclimatic limits for each PFT based on modified
! LPJ lookup-table
!

! input
    LOGICAL,  INTENT(in)  :: is_dynamic        ! mask indicating whether a PFT is subject to vegetation dynamics
    REAL(wp), INTENT(in)  :: lct_tcmin         ! PFT-specific minimum coldest-month temperature limit [C]
    REAL(wp), INTENT(in)  :: lct_tcmax         ! PFT-specific maximum coldest-month temperature limit [C]
    REAL(wp), INTENT(in)  :: lct_twmax         ! PFT-specific maximum warmest-month temperature limit [C]
    REAL(wp), INTENT(in)  :: lct_min_temprange ! PFT-specific minimum difference of 20-year average warme st
                                               ! minus coldest month temperature [C]
    REAL(wp), INTENT(in)  :: lct_gddmin        ! PFT-specific minimum GDD limit
    REAL(wp), INTENT(in)  :: min_mmtemp20      ! Minimum monthly mean temp. (20yr climatology) [K]
    REAL(wp), INTENT(in)  :: max_mmtemp20      ! Maximum monthly mean temp. (20yr climatology) [K]
    REAL(wp), INTENT(in)  :: gdd_prev_year     ! GDD of previous year

! output
    REAL(wp), INTENT(inout) :: bio_exist

!------------------------------------------------------------------------------

    bio_exist = 1._wp

    IF (is_dynamic .AND. (min_mmtemp20 - tmelt < lct_tcmin .OR. min_mmtemp20 - tmelt > lct_tcmax  &
       .OR. gdd_prev_year < lct_gddmin .OR. max_mmtemp20 - tmelt > lct_twmax                      &
       .OR. (max_mmtemp20 - min_mmtemp20) < lct_min_temprange)) bio_exist = 0._wp

  END SUBROUTINE bioclim_limits


! --- potential_tree_fpc ---------------------------------------------------------------------------

  SUBROUTINE potential_tree_fpc(nc, ntiles, is_dynamic, is_woody, bio_exist, &
                                npp_ave, act_fpc, pot_fpc)

!
! subroutine calculates potential FPC (in absence of disturbances) based on NPP for each PFT
! called once a year!

    INTEGER,  INTENT(in)  :: nc              ! vector length
    INTEGER,  INTENT(in)  :: ntiles          ! number of tiles (different land cover type mosaics per gridcell)
    LOGICAL,  INTENT(in)  :: is_dynamic(:)   ! flag indicating dynamic PFTs
    LOGICAL,  INTENT(in)  :: is_woody(:)     ! flag indicating woody PFTs
    REAL(wp), INTENT(in)  :: bio_exist(:,:)  ! bio_exist=1. indicates that all bio-climatic limits are fulfilled
    REAL(wp), INTENT(in)  :: npp_ave(:,:)    ! NPP averaged over 5 years
    REAL(wp), INTENT(in)  :: act_fpc(:,:)    ! actual foliage projective cover

    REAL(wp), INTENT(out) :: pot_fpc(:,:)

    ! Local variables
    INTEGER           :: ic, itile
    REAL(wp)          :: sum_npp(nc)

    CHARACTER(len=*), PARAMETER :: routine = modname//':potential_tree_fpc'

    !$ACC DATA ASYNC(acc_stream) CREATE(sum_npp)

    ! Summing up NPP
    !$ACC PARALLEL LOOP DEFAULT(PRESENT) ASYNC(acc_stream)
    DO ic = 1, nc
      sum_npp(ic) = 0._wp
    END DO
    !$ACC END PARALLEL LOOP

    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(acc_stream)
    !$ACC LOOP GANG VECTOR COLLAPSE(2)
    DO itile = 1, ntiles
      DO ic = 1, nc
        IF (is_woody(itile) .AND. is_dynamic(itile) .AND. &
            & npp_ave(ic,itile) > sum_npp_min .AND. bio_exist(ic,itile) > 0.5_wp) THEN
          sum_npp(ic) = sum_npp(ic) + npp_ave(ic,itile) ** npp_nonlinearity * MAX(act_fpc_min,act_fpc(ic,itile))
        END IF
      END DO
    END DO
    !$ACC END PARALLEL LOOP

    ! Determine potential FPC
    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(acc_stream)
    !$ACC LOOP GANG VECTOR COLLAPSE(2)
    DO itile = 1, ntiles
      DO ic = 1, nc
        IF (is_woody(itile) .AND. is_dynamic(itile) .AND. npp_ave(ic,itile) > sum_npp_min .AND. &
            & bio_exist(ic,itile) > 0.5_wp) THEN
          pot_fpc(ic,itile) = npp_ave(ic,itile) ** npp_nonlinearity / sum_npp(ic) * &
                             tree_fpc_max * MAX(act_fpc_min,act_fpc(ic,itile))
        ELSE
          pot_fpc(ic,itile) = 0._wp
        END IF
      END DO
    END DO
    !$ACC END PARALLEL LOOP

    !$ACC END DATA

  END SUBROUTINE potential_tree_fpc


! --- fpc_to_cover_fract_pot -------------------------------------------------------------------------

  SUBROUTINE fpc_to_cover_fract_pot (nc, ntiles, is_dynamic, is_woody, act_fpc, bare_fpc, cover_fract_pot)
!
! Convert fractional plant cover (act_fpc) calculated within nlcc to jsbach
! cover fractions reflecting natural vegetation, only (cover_fract_pot).
! The only difference between act_fpc and cover_fract_pot is caused by bare_fpc,
! which is added to the grass cover fractions.
!
! input
    INTEGER,  INTENT(in)  :: nc              ! vector length
    INTEGER,  INTENT(in)  :: ntiles          ! number of tiles (different land cover type mosaics per gridcell)
    LOGICAL,  INTENT(in)  :: is_dynamic(:)   !
    LOGICAL,  INTENT(in)  :: is_woody(:)     !
    REAL(wp), INTENT(in)  :: act_fpc(:,:)    ! actual fpc in nlcc
    REAL(wp), INTENT(in)  :: bare_fpc(:)     ! nlcc bare soil fraction
! output
    REAL(wp), INTENT(out) :: cover_fract_pot(:,:) ! cover fractions if there was only natural vegetation

! local variables
    INTEGER   :: ic, itile
    REAL(wp)  :: sum_grass_fpc(nc)     ! part of the vegetated area covered by grass

    CHARACTER(len=*), PARAMETER :: routine = modname//':fpc_to_cover_fract_pot'

    !$ACC DATA ASYNC(acc_stream) CREATE(sum_grass_fpc)

! -- initialization
    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(acc_stream)
    !$ACC LOOP GANG VECTOR COLLAPSE(2)
    DO itile = 1, ntiles
      DO ic = 1, nc
        cover_fract_pot(ic,itile) = 0._wp
      END DO
    END DO
    !$ACC END PARALLEL LOOP

! -- find out fpc of grass
    !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) ASYNC(acc_stream)
    DO ic = 1, nc
      sum_grass_fpc(ic) = 0._wp
    END DO
    !$ACC END PARALLEL LOOP

    DO itile = 1, ntiles
      IF (.NOT. is_woody(itile) .AND. is_dynamic(itile)) THEN
        !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) ASYNC(acc_stream)
        DO ic = 1, nc
          sum_grass_fpc(ic) = sum_grass_fpc(ic) + act_fpc(ic,itile)
        END DO
        !$ACC END PARALLEL LOOP
      END IF
    END DO

!-- calculate new cover_fractions
    DO itile = 1, ntiles
      IF (is_dynamic(itile) .AND. is_woody(itile)) THEN
        !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) ASYNC(acc_stream)
        DO ic = 1, nc
          cover_fract_pot(ic,itile) = act_fpc(ic,itile)
        END DO
        !$ACC END PARALLEL LOOP
      ELSE IF (is_dynamic(itile)) THEN
        !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) ASYNC(acc_stream)
        DO ic = 1, nc
          cover_fract_pot(ic,itile) = act_fpc(ic,itile) * (1._wp + bare_fpc(ic) / MAX(sum_grass_fpc(ic),EPSILON(1._wp)))
        END DO
        !$ACC END PARALLEL LOOP
      END IF
    END DO

    CALL scale_cover_fract (nc, ntiles, is_dynamic(:), cover_fract_pot(:,:))

    !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) ASYNC(acc_stream)
    DO ic = 1, nc
      IF (  SUM(cover_fract_pot(ic,:)) > 1._wp + REAL(ntiles,wp)*EPSILON(1._wp) .OR. &
          & SUM(cover_fract_pot(ic,:)) < 1._wp - REAL(ntiles,wp)*EPSILON(1._wp)) THEN
        WRITE (message_text,*) 'Sum of cover_fract_pot /= 1: ', SUM(cover_fract_pot(ic,:)), ' at ic: ', ic
        CALL finish (routine, message_text)
      END IF
    END DO
    !$ACC END PARALLEL LOOP

    DO itile = 1, ntiles
      !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) ASYNC(acc_stream)
      DO ic = 1, nc
        IF (is_dynamic(itile) .AND. cover_fract_pot(ic,itile) < fract_small) THEN
          WRITE (message_text,*) 'cover_fract_pot too small: ', cover_fract_pot(ic,itile), &
            & ' at ic: ', ic, ' itile: ', itile
          CALL finish (routine, message_text)
        END IF
      END DO
      !$ACC END PARALLEL LOOP
    END DO

    !$ACC END DATA

  END SUBROUTINE fpc_to_cover_fract_pot


! --- cover_fract_pot_to_cover_fract ---------------------------------------------------------------------------

  SUBROUTINE cover_fract_pot_to_cover_fract (nc, ntiles, is_dynamic, is_woody, is_pasture, is_crop, &
    &                                        cover_fract_pot, cover_fract_pot_previous_year, cover_fract)

!
! convert the fractional plant cover of potential natural vegetation to plant cover fractions including agricultural areas.
! Pasture is first established on grasslands, whereas crops are established at the expense of all natural cover types.
!
! In runs with land use transitions, the transitions from natural forests or grass lands to crops and pastures are predefined.
! It is not wanted, that the dynamic vegetation interferes with these transitions by e.g. re-establishing forests where they
! had just been removed. Thus the changes in forest fraction calculated by the dynamic vegetation are not instantanously
! applied. This results in an imbalance between actual cover_fractions and the ones the dynamic vegetation calculates for
! instantanous establishement of crops and pasture (cover_fract_inst).
! The final cover_fractions are based on cover_fract_inst with the exception, that the forest fraction is based on the forest
! fraction after land use transitions (cover_fract entering the routine), and only the natural change in cover_fractions in
! the direction to reduce the inbalance is taken into account (see below).
!
! input
    INTEGER,  INTENT(in)  :: nc                   ! vector length
    INTEGER,  INTENT(in)  :: ntiles               ! number of tiles (different land cover type mosaics per gridcell)
    LOGICAL,  INTENT(in)  :: is_dynamic(:)        ! PFTs taking part in competition
    LOGICAL,  INTENT(in)  :: is_woody(:)          ! PFTs that build wood
    LOGICAL,  INTENT(in)  :: is_pasture(:)        ! logical mask for pasture
    LOGICAL,  INTENT(in)  :: is_crop(:)           ! logical mask for crops
    REAL(wp), INTENT(in)  :: cover_fract_pot(:,:) ! cover fractions if there was only natural vegetation
    REAL(wp), INTENT(in)  :: cover_fract_pot_previous_year(:,:) ! last years cover_fract_pot

! inout
    REAL(wp), INTENT(inout) :: cover_fract(:,:)  ! cover fraction within jsbach

! local variables

    INTEGER   :: ic, itile
    INTEGER   :: nwoody                       ! number of woody PFTs
    INTEGER   :: nexcluded                    ! number of non-dynamic PFTs (crops, pastures)
    REAL(wp)  :: cover_fract_inst(nc,ntiles)  ! cover fractions that would result from an instantanous implementation of
                                                ! land use change and vegetation dynamics
    REAL(wp)  :: excluded_fract(nc)           ! fraction of the vegetated part of the grid box not considered for vegetation
                                                ! dynamics (e.g. agricultural areas)
    REAL(wp)  :: sum_woody_fract(nc)          ! fraction of the vegetated part of the grid box covered by woody plants
    REAL(wp)  :: sum_grass_fract(nc)          ! fraction of the vegetated part of the grid box covered by grass
    REAL(wp)  :: sum_woody_fract_inst(nc)     ! fraction of the vegetated part of the grid box covered by woody plants
                                                ! based on cover_fract_inst
    REAL(wp)  :: sum_woody_fract_pot(nc)      ! woody fraction of the vegetated part of the grid box considering natural
                                                ! vegetation, only
    REAL(wp)  :: sum_grass_fract_pot(nc)      ! grass fraction of the vegetated part of the grid box considering natural
                                                ! vegetation, only
    REAL(wp)  :: sum_woody_fract_old(nc)      ! fraction of the vegetated part of the grid box covered by woody plants
                                                ! (based on act_fpc_previous_year) that would have resulted from an
                                                ! instantaneous establishment of crop land and pasture
    REAL(wp)  :: sum_woody_fract_act(nc)      ! fraction of the vegetated part of the grid box covered by woody plants
    REAL(wp)  :: sum_grass_fract_act(nc)      ! fraction of the vegetated part of the grid box covered by grasses
    REAL(wp)  :: delta_woody(nc)              ! deviation of the fraction of the vegetated part of the grid box covered
                                                ! by woody plants from the fraction that would result from an instantaneous
                                                !  establishment of crop land and pasture

    CHARACTER(len=*), PARAMETER :: routine = modname//':cover_fract_pot_to_cover_fract'

    !$ACC DATA ASYNC(acc_stream) &
    !$ACC   CREATE(cover_fract_inst, excluded_fract, sum_woody_fract, sum_grass_fract, sum_woody_fract_inst) &
    !$ACC   CREATE(sum_woody_fract_pot, sum_grass_fract_pot, sum_woody_fract_old, sum_woody_fract_act) &
    !$ACC   CREATE(sum_grass_fract_act, delta_woody)

    !-- Find out number of woody types and number of exculded types
    nwoody = 0
    nexcluded = 0
    DO itile = 1, ntiles
      IF (is_dynamic(itile) .AND. is_woody(itile)) nwoody = nwoody + 1
      IF (.NOT. is_dynamic(itile))                 nexcluded = nexcluded + 1
    END DO

    !-- Sum up cover fractions of woody types, grasses and fractions not affected by NLCC.
    !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) ASYNC(acc_stream)
    DO ic = 1, nc
      sum_woody_fract(ic) = 0._wp
      sum_grass_fract(ic) = 0._wp
      excluded_fract(ic) = 0._wp
    END DO
    !$ACC END PARALLEL LOOP

    DO itile = 1, ntiles
      IF (is_woody(itile) .AND. is_dynamic(itile)) THEN
        !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) ASYNC(acc_stream)
        DO ic = 1, nc
          sum_woody_fract(ic) = sum_woody_fract(ic) + cover_fract(ic,itile)
        END DO
        !$ACC END PARALLEL LOOP
      ELSE IF (.NOT. is_woody(itile) .AND. is_dynamic(itile)) THEN
        !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) ASYNC(acc_stream)
        DO ic = 1, nc
          sum_grass_fract(ic) = sum_grass_fract(ic) + cover_fract(ic,itile)
        END DO
        !$ACC END PARALLEL LOOP
      ELSE IF (.NOT. is_dynamic(itile)) THEN
        !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) ASYNC(acc_stream)
        DO ic = 1, nc
          excluded_fract(ic) = excluded_fract(ic) + cover_fract(ic,itile)
        END DO
        !$ACC END PARALLEL LOOP
      END IF
    END DO

    !-- Sum up woody and grass fractions of the potential natural vegetation based on current vegetation dynamics
    !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) ASYNC(acc_stream)
    DO ic = 1, nc
      sum_woody_fract_pot(ic) = 0._wp
      sum_grass_fract_pot(ic) = 0._wp
    END DO
    !$ACC END PARALLEL LOOP

    DO itile = 1, ntiles
      IF (is_woody(itile) .AND. is_dynamic(itile)) THEN
        !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) ASYNC(acc_stream)
        DO ic = 1, nc
          sum_woody_fract_pot(ic) = sum_woody_fract_pot(ic) + cover_fract_pot(ic,itile)
        END DO
        !$ACC END PARALLEL LOOP
      ELSE IF (.NOT. is_woody(itile) .AND. is_dynamic(itile)) THEN
        !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) ASYNC(acc_stream)
        DO ic = 1, nc
          sum_grass_fract_pot(ic) = sum_grass_fract_pot(ic) + cover_fract_pot(ic,itile)
        END DO
        !$ACC END PARALLEL LOOP
      END IF
    END DO

    !-- Calculate the woody fraction resulting from last years potential vegetation assuming an instantaneous
    !   establishment of crops and pastures
    CALL calc_cover_fract_inst(nc, ntiles, is_woody(:), is_dynamic(:), is_pasture(:), is_crop(:),         &
                               cover_fract_pot_previous_year(:,:), cover_fract(:,:), cover_fract_inst(:,:))

    !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) ASYNC(acc_stream)
    DO ic = 1, nc
      sum_woody_fract_old(ic) = 0._wp
    END DO
    !$ACC END PARALLEL LOOP

    DO itile = 1, ntiles
      IF (is_woody(itile) .AND. is_dynamic(itile)) THEN
        !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) ASYNC(acc_stream)
        DO ic = 1, nc
          sum_woody_fract_old(ic) = sum_woody_fract_old(ic) + cover_fract_inst(ic,itile)
        END DO
        !$ACC END PARALLEL LOOP
      END IF
    END DO

    !-- Calculate the woody fraction resulting from current potential vegetation assuming an instantaneous
    !   establishment of crops and pastures
    CALL calc_cover_fract_inst(nc, ntiles, is_woody(:), is_dynamic(:), is_pasture(:), is_crop(:), &
                               cover_fract_pot(:,:), cover_fract(:,:), cover_fract_inst(:,:))

    !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) ASYNC(acc_stream)
    DO ic = 1, nc
      sum_woody_fract_inst(ic) = 0._wp
    END DO
    !$ACC END PARALLEL LOOP

    DO itile = 1, ntiles
      IF (is_woody(itile) .AND. is_dynamic(itile)) THEN
        !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) ASYNC(acc_stream)
        DO ic = 1, nc
          sum_woody_fract_inst(ic) = sum_woody_fract_inst(ic) + cover_fract_inst(ic,itile)
        END DO
        !$ACC END PARALLEL LOOP
      END IF
    END DO

    !-- Calculation of the final woody fraction. A difference between sum_woody_fract and sum_woody_fract_old indicates
    !   an imbalance of the state of the dynamic vegetation and the cover fractions calculated from landuse transitions.
    !   On the other hand, the difference between sum_woody_fract_inst and sum_woody_fract_old show the trend of the woody
    !   fraction caused by vegetation dynamics.
    !   The idea of this routine is, to allow shifts in forest fractions only if they reduce the imbalance between nlcc
    !   and landuse change.
    !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) ASYNC(acc_stream)
    DO ic = 1, nc
      IF (sum_woody_fract(ic) - sum_woody_fract_old(ic) > 0._wp) THEN
        ! -- woody fraction of cover_fract entering the routine (last updated by last years land use change) is greater
        !    than the the woody fraction resulting from last years potential vegetation and instantanous establishment of
        !    crops and pastures. Thus there is an inbalance, the woody fraction of cover_fract is too high.
        IF (sum_woody_fract_inst(ic) - sum_woody_fract_old(ic) > 0._wp) THEN
          ! -- vegetation dynamics even enlarge woody fraction
          !    ==> no change in actual woody fraction, unless the enlargement by vegetation dynamics exceeds the initial
          !        woody surplus
          delta_woody(ic) = MAX(0._wp, sum_woody_fract_inst(ic) - sum_woody_fract(ic))
        ELSE
          ! -- vegetation dynamics shrink woody fraction
          delta_woody(ic) = sum_woody_fract_inst(ic) - sum_woody_fract_old(ic)
        END IF
      ELSE
        ! -- The woody fraction of cover_fract is too small
        IF (sum_woody_fract_inst(ic) - sum_woody_fract_old(ic) > 0._wp) THEN
          ! -- vegetation dynamics enlarge woody fraction
          delta_woody(ic) = sum_woody_fract_inst(ic) - sum_woody_fract_old(ic)
        ELSE
          ! -- vegetation dynamics shrink woody fraction
          !    ==> no change in actual woody fraction, unless the reduction by vegetation dynamics exceeds the initial
          !        woody deficit
          delta_woody(ic) = MIN(0._wp, sum_woody_fract_inst(ic) - sum_woody_fract(ic))
        END IF
      END IF

      ! -- in some cases delta_woody can not be used to calculate the actual woody fraction. If there is no grass left
      !    after instantanous establishement of agricultural areas, the woody fractions of cover_fract_inst and
      !    cover_fract_inst_old are identical and delta_woody is zero, even if vegetation dynamics change considerably.
      !    The following formulation is needed to assure that sum_woody_fract_act is smaller than sum_woody_fract_pot and
      !    sum_grass_fract_act is smaller than sum_grass_fract_pot in all cases.
      sum_woody_fract_act(ic) = MIN(sum_woody_fract_pot(ic), &
        &                           MAX(REAL(nwoody,wp)*fract_small, &
        &                               sum_woody_fract_pot(ic) - excluded_fract(ic) + REAL(nexcluded,wp)*fract_small, &
        &                               sum_woody_fract(ic) + delta_woody(ic)))
    END DO
    !$ACC END PARALLEL LOOP

    !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) ASYNC(acc_stream)
    DO ic = 1, nc
      sum_grass_fract_act(ic) = 1._wp - sum_woody_fract_act(ic) - excluded_fract(ic)
    END DO
    !$ACC END PARALLEL LOOP

    !-- Use actual woody fraction to determine all actual fractions
    DO itile = 1,ntiles
      IF (is_dynamic(itile) .AND. is_woody(itile)) THEN
        !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) ASYNC(acc_stream)
        DO ic = 1, nc
          cover_fract(ic,itile) = cover_fract_pot(ic,itile) * sum_woody_fract_act(ic) / sum_woody_fract_pot(ic)
        END DO
        !$ACC END PARALLEL LOOP
      ELSE IF (is_dynamic(itile)) THEN
        !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) ASYNC(acc_stream)
        DO ic = 1, nc
          cover_fract(ic,itile) = cover_fract_pot(ic,itile) * sum_grass_fract_act(ic) / sum_grass_fract_pot(ic)
        END DO
        !$ACC END PARALLEL LOOP
      END IF
    END DO

    CALL scale_cover_fract(nc, ntiles, is_dynamic(:), cover_fract(:,:))

    !-- Stop simulation, if sum of cover fractions is not equal to 1 or a cover fraction is smaller than fract_small
    !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) ASYNC(acc_stream)
    DO ic = 1, nc
      IF (  SUM(cover_fract(ic,:)) > 1._wp + REAL(ntiles,wp)*EPSILON(1._wp) .OR. &
          & SUM(cover_fract(ic,:)) < 1._wp - REAL(ntiles,wp)*EPSILON(1._wp)) THEN
        WRITE (message_text,*) 'Sum of cover_fract /= 1: ', SUM(cover_fract(ic,:)), ' at ic: ', ic
        CALL finish (routine, message_text)
      END IF
    END DO
    !$ACC END PARALLEL LOOP
    DO itile = 1, ntiles
      !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) ASYNC(acc_stream)
      DO ic = 1, nc
        IF (cover_fract(ic,itile) < fract_small) THEN
          WRITE (message_text,*) 'cover_fract too small: ', cover_fract(ic,itile), ' at ic: ', ic
          CALL finish (routine, message_text)
        END IF
      END DO
      !$ACC END PARALLEL LOOP
    END DO

    !&ACC END DATA

  END SUBROUTINE cover_fract_pot_to_cover_fract


! --- desert_fraction --------------------------------------------------------

  SUBROUTINE desert_fraction(nc, ntiles, init_running_means, accelerate_nlcc,     &
    &                        is_woody, is_dynamic,                                &
    &                        act_fpc, bare_fpc,                                   &
    &                        sla, max_green_bio,                                  &
    &                        sum_green_bio_memory, desert_fpc)
! !DESCRIPTION:
!
! Calculation of desert fraction from green biomass
!
!------------------------------------------------------------------------------

! input
    INTEGER,  INTENT(in)  :: nc                ! vector length
    INTEGER,  INTENT(in)  :: ntiles            ! number of tiles (different land cover type mosaics per gridcell)
    LOGICAL,  INTENT(in)  :: init_running_means! initialisation of sum_green_bio_memory
    REAL(wp), INTENT(in)  :: accelerate_nlcc   ! acceleration factor for vegetation and desert dynamics
    LOGICAL,  INTENT(in)  :: is_woody(:)       ! flag to label a woody PFT
    LOGICAL,  INTENT(in)  :: is_dynamic(:)     ! flag to label PFTs that take part in the natural vegetation dynamics
    REAL(wp), INTENT(in)  :: act_fpc(:,:)      ! actual FPC
    REAL(wp), INTENT(in)  :: bare_fpc(:)       ! bare FPC
    REAL(wp), INTENT(in)  :: sla(:,:)          ! specific leaf area
    REAL(wp), INTENT(in)  :: max_green_bio(:,:)! maximum value of green biomass within a year

! input/output
    REAL(wp), INTENT(inout) :: sum_green_bio_memory(:) ! vegetated fraction calculated from green biomass

! output
    REAL(wp), INTENT(out) :: desert_fpc(:) ! desert FPC

! local variables
    INTEGER           :: ic, itile
    REAL(wp)          :: sum_green_bio(nc)
    REAL(wp)          :: sum_act_fpc(nc)
    REAL(wp)          :: sum_grass_fpc(nc)
    INTEGER           :: n_grass_pft

    CHARACTER(len=*), PARAMETER :: routine = modname//':desert_fraction'
!------------------------------------------------------------------------------
    !$ACC DATA ASYNC(acc_stream) &
    !$ACC   CREATE(sum_green_bio, sum_act_fpc, sum_grass_fpc)

    n_grass_pft = 0
    !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) ASYNC(acc_stream)
    DO ic = 1, nc
      sum_green_bio(ic) = 0._wp
      sum_act_fpc(ic)   = 0._wp
      sum_grass_fpc(ic) = 0._wp
    END DO
    !$ACC END PARALLEL LOOP

    DO itile = 1, ntiles
      IF (.NOT. is_woody(itile) .AND. is_dynamic(itile)) THEN
        n_grass_pft = n_grass_pft + 1
        !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) ASYNC(acc_stream)
        DO ic = 1, nc
          sum_grass_fpc(ic) = sum_grass_fpc(ic) + act_fpc(ic,itile)
        END DO
        !$ACC END PARALLEL LOOP
      END IF
    END DO

    DO itile = 1, ntiles
      IF (is_woody(itile) .AND. is_dynamic(itile)) THEN
        !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) ASYNC(acc_stream)
        DO ic = 1, nc
          sum_green_bio(ic) = sum_green_bio(ic) + MAX(0._wp, act_fpc(ic,itile) * (1._wp -          &
            &       exp(-desert_extend * ((max_green_bio(ic,itile) * sla(ic,itile)) ** desert_margin) /   &
            &       (3._wp ** (desert_margin - 1._wp)))))
          sum_act_fpc(ic) = sum_act_fpc(ic) + act_fpc(ic,itile)
        END DO
        !$ACC END PARALLEL LOOP
      ELSE IF (.NOT. is_woody(itile) .AND. is_dynamic(itile)) THEN
        !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) ASYNC(acc_stream)
        DO ic = 1, nc
          IF (sum_grass_fpc(ic) > EPSILON(1._wp)) THEN
            sum_green_bio(ic) = sum_green_bio(ic) +                                                    &
              &   MAX(0._wp, act_fpc(ic,itile) * (1._wp + bare_fpc(ic) / sum_grass_fpc(ic)) * (1._wp - &
              &   exp(-desert_extend * ((max_green_bio(ic,itile) * sla(ic,itile)) ** desert_margin) /  &
              &   (3._wp ** (desert_margin - 1._wp)))))
            sum_act_fpc(ic) = sum_act_fpc(ic) + act_fpc(ic,itile) * (1._wp + bare_fpc(ic) / sum_grass_fpc(ic))
          ELSE
            sum_green_bio(ic) = sum_green_bio(ic) +                                                   &
              &   MAX(0._wp, (bare_fpc(ic) / REAL(n_grass_pft,wp)) * (1._wp -                         &
              &   exp(-desert_extend * ((max_green_bio(ic,itile) * sla(ic,itile)) ** desert_margin) / &
              &   (3._wp ** (desert_margin - 1._wp)))))
            sum_act_fpc(ic) = sum_act_fpc(ic) + (bare_fpc(ic) / REAL(n_grass_pft,wp))
           END IF
        END DO
        !$ACC END PARALLEL LOOP
      END IF
    END DO

    !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) ASYNC(acc_stream)
    DO ic = 1, nc
      IF (sum_act_fpc(ic) > EPSILON(1._wp)) THEN
        sum_green_bio(ic) = sum_green_bio(ic) / sum_act_fpc(ic)
      END IF
    END DO
    !$ACC END PARALLEL LOOP
    IF (init_running_means) THEN
      !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) ASYNC(acc_stream)
      DO ic = 1, nc
        sum_green_bio_memory(ic) = sum_green_bio(ic)
      END DO
      !$ACC END PARALLEL LOOP
    ELSE
      !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) ASYNC(acc_stream)
      DO ic = 1, nc
        sum_green_bio_memory(ic) = (sum_green_bio_memory(ic) * (tau_desert/accelerate_nlcc - 1._wp) + sum_green_bio(ic)) &
          &                      / (tau_desert/accelerate_nlcc)
      END DO
      !$ACC END PARALLEL LOOP
    END IF

    !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) ASYNC(acc_stream)
    DO ic = 1, nc
      desert_fpc(ic) = MIN(1._wp - 2._wp * EPSILON(1._wp),MAX(0._wp,1._wp - sum_green_bio_memory(ic)))
    END DO
    !$ACC END PARALLEL LOOP

    !$ACC END DATA

  END SUBROUTINE desert_fraction

! --- fpc_daily ---------------------------------------------------------------------------

  SUBROUTINE fpc_daily(nc, ntiles, accelerate_nlcc, is_woody, is_dynamic, tau_Cpool_woods,    &
                       npp_ave, bio_exist, burned_fract, damaged_fract, pot_fpc, &
                       act_fpc, bare_fpc, fuel)

! Calculation of daily FPC for trees, shrubs, and herbaceous plants

! input
    INTEGER,         INTENT(in) :: nc              ! vector length
    INTEGER,         INTENT(in) :: ntiles          ! number of tiles (different land cover type mosaics per gridcell)
    REAL(wp),        INTENT(in) :: accelerate_nlcc ! factor to accellerate vegetation and desert dynamics
    LOGICAL,         INTENT(in) :: is_woody(:)
    LOGICAL,         INTENT(in) :: is_dynamic(:)
    REAL(wp),        INTENT(in) :: tau_Cpool_woods(:,:) ! lifetime of Cpool_woods [days]
    REAL(wp),        INTENT(in) :: npp_ave(:,:)
    REAL(wp),        INTENT(in) :: bio_exist(:,:)
    REAL(wp),        INTENT(in) :: burned_fract(:,:)
    REAL(wp),        INTENT(in) :: damaged_fract(:,:)
    REAL(wp),        INTENT(in) :: pot_fpc(:,:)

! inout
    REAL(wp), INTENT(inout) :: act_fpc(:,:)
    REAL(wp), INTENT(inout) :: bare_fpc(:)
    REAL(wp), INTENT(inout), optional :: fuel(:)

! local variables
    INTEGER     :: ic, itile
    REAL(wp)    :: act_fpc_previous_day(nc,ntiles)
    REAL(wp)    :: total_act_fpc(nc)
    REAL(wp)    :: grass_sum(nc)
    REAL(wp)    :: non_woody_fpc(nc)
    REAL(wp)    :: theta, veg_inc, tau_pft
    REAL(wp), PARAMETER  ::  min_dist_woody = 0.002_wp / 365._wp !! mortality of woody plants (not caused by fire/windthrow,
                                                                 !! but by other unresolved processes)
    REAL(wp), PARAMETER  ::  grass_mortality = 0.01_wp / 365._wp !! grass mortality (not caused by fire, but by other
                                                                 !! unresolved processes)
    REAL(wp), PARAMETER  ::  tau_Cpool_woods_to_tau_pft = 1._wp

    CHARACTER(len=*), PARAMETER :: routine = modname//':fpc_daily'

    !$ACC DATA ASYNC(acc_stream) &
    !$ACC   CREATE(act_fpc_previous_day, total_act_fpc, grass_sum, non_woody_fpc)

    ! initialisations

    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(acc_stream)
    !$ACC LOOP GANG VECTOR COLLAPSE(2)
    DO itile = 1, ntiles
      DO ic = 1, nc
        act_fpc_previous_day(:,:) = act_fpc(:,:)
      END DO
    END DO
    !$ACC END PARALLEL LOOP

    !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) ASYNC(acc_stream)
    DO ic = 1, nc
      non_woody_fpc(ic) = bare_fpc(ic)
    END DO
    !$ACC END PARALLEL LOOP

    DO itile = 1, ntiles
      IF (.NOT. is_dynamic(itile)) THEN
        !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) ASYNC(acc_stream)
        DO ic = 1, nc
          act_fpc(ic,itile) = 0._wp
        END DO
        !$ACC END PARALLEL LOOP
      END IF
      IF (.NOT. is_woody(itile)) THEN
        !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) ASYNC(acc_stream)
        DO ic = 1, nc
          non_woody_fpc(ic) = non_woody_fpc(ic) + act_fpc(ic,itile)
        END DO
        !$ACC END PARALLEL LOOP
      END IF
    END DO

    !-- adjust grass fractions to changes in woody type fpc

    !-- reset temporary arrays
    !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) ASYNC(acc_stream)
    DO ic = 1, nc
      grass_sum(ic) = 0._wp
    END DO
    !$ACC END PARALLEL LOOP

    !-- temporary sums for grass types
    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(acc_stream)
    !$ACC LOOP GANG VECTOR COLLAPSE(2)
    DO itile = 1, ntiles
      DO ic = 1, nc
        IF (.NOT. is_woody(itile) .AND. is_dynamic(itile) .AND. bio_exist(ic,itile) > 0.5_wp .AND. &
            &  npp_ave(ic,itile) > EPSILON(1._wp)) THEN
          grass_sum(ic) = grass_sum(ic) + npp_ave(ic,itile) * MAX(act_fpc_min,act_fpc(ic,itile))
        END IF
      END DO
    END DO
    !$ACC END PARALLEL LOOP

!!$    !-- area wind break of woody types (trees and shrubs) and burned
!!$    IF (PRESENT(burned_frac)) THEN
!!$      CALL disturbed_frac(lctlib, nidx,kidx0,kidx1,ntiles,DIST_FIRE, &
!!$                          dynamic_pft(:) .OR. (dist_opts%lburn_pasture .AND. lctlib%pasture_pft(:)), with_yasso, &
!!$                          act_fpc,veg_fract_correction,surf,cbal,climbuf,burned_frac,burned_frac_diag,fuel)
!!$      CALL disturbed_frac(lctlib, nidx,kidx0,kidx1,ntiles,DIST_WINDBREAK,woody_pft(:) .AND. dynamic_pft(:), with_yasso, &
!!$                          act_fpc,veg_fract_correction,surf,cbal,climbuf,damaged_frac)
!!$    END IF

    !-- dynamic equation for act_fpc, woody types


    DO ic = 1, nc
      veg_inc = 0._wp
      DO itile = 1, ntiles
        IF (is_dynamic(itile) .AND. is_woody(itile)) THEN
          veg_inc = veg_inc + MAX(fract_small,act_fpc(ic,itile) + accelerate_nlcc * &
          ((pot_fpc(ic,itile) - act_fpc(ic,itile)) / tau_pft &
          - (burned_fract(ic,itile) + damaged_fract(ic,itile) + min_dist_woody) * act_fpc_previous_day(ic,itile))) &
          - act_fpc(ic,itile)
        ELSEIF(is_dynamic(itile) .AND. .NOT. is_woody(itile)) THEN
          IF (bio_exist(ic,itile) > 0.5_wp .AND. npp_ave(ic,itile) > EPSILON(1._wp)) THEN
            veg_inc = veg_inc + MAX(fract_small,act_fpc(ic,itile) - accelerate_nlcc * &
            (burned_fract(ic,itile) * act_fpc_previous_day(ic,itile) + &
            grass_mortality * act_fpc(ic,itile) - &
            bare_fpc(ic) / grass_sum(ic) * npp_ave(ic,itile) * MAX(act_fpc_min,act_fpc(ic,itile)) / tau_pft)) &
            - act_fpc(ic,itile)
          ELSE
            veg_inc = veg_inc + MAX(fract_small,act_fpc(ic,itile) - accelerate_nlcc * &
            (burned_fract(ic,itile) * act_fpc_previous_day(ic,itile) + &
            act_fpc(ic,itile) / tau_pft)) &
            - act_fpc(ic,itile)
          ENDIF
        ENDIF
      END DO
      theta = 1._wp
      IF (bare_fpc(ic) < (veg_inc + fract_small)) theta = 0._wp
      DO itile = 1, ntiles
        tau_pft = tau_Cpool_woods(ic,itile) * tau_Cpool_woods_to_tau_pft
        IF (is_dynamic(itile) .AND. is_woody(itile)) THEN
          act_fpc(ic,itile) = MAX(fract_small,act_fpc(ic,itile) + accelerate_nlcc * &
            ((theta * pot_fpc(ic,itile) - act_fpc(ic,itile)) / tau_pft                &
            - (burned_fract(ic,itile) + damaged_fract(ic,itile) + min_dist_woody) * act_fpc_previous_day(ic,itile)))
        ENDIF
      ENDDO
    ENDDO

    !-- dynamic equation for act_fpc, grass types

    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(acc_stream)
    !$ACC LOOP GANG VECTOR COLLAPSE(2)
    DO itile = 1, ntiles
      DO ic = 1, nc
        tau_pft = tau_Cpool_woods(ic,itile) * tau_Cpool_woods_to_tau_pft
      END DO
    END DO
    !$ACC END PARALLEL LOOP
    DO itile = 1, ntiles
      IF (is_dynamic(itile) .AND. .NOT. is_woody(itile)) THEN
        IF (bio_exist(ic,itile) > 0.5_wp .AND. npp_ave(ic,itile) > EPSILON(1._wp)) THEN
          !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) ASYNC(acc_stream)
          DO ic = 1, nc
            act_fpc(ic,itile) = MAX(fract_small,act_fpc(ic,itile) - accelerate_nlcc * &
               (burned_fract(ic,itile) * act_fpc_previous_day(ic,itile) + &
               grass_mortality * act_fpc(ic,itile) - &
               bare_fpc(ic) / grass_sum(ic) * npp_ave(ic,itile) * MAX(act_fpc_min,act_fpc(ic,itile)) / tau_pft))
          END DO
          !$ACC END PARALLEL LOOP
        ELSE
          !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) ASYNC(acc_stream)
          DO ic = 1, nc
            act_fpc(ic,itile) = MAX(fract_small,act_fpc(ic,itile) - accelerate_nlcc * &
               (burned_fract(ic,itile) * act_fpc_previous_day(ic,itile) + &
               act_fpc(ic,itile) / tau_pft))
          END DO
          !$ACC END PARALLEL LOOP
        ENDIF
      ENDIF
    ENDDO

    !-- calculate total area fraction covered by dynamic vegetation

    !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) ASYNC(acc_stream)
    DO ic = 1, nc
      total_act_fpc(ic) = 0._wp
    END DO
    !$ACC END PARALLEL LOOP

    DO itile = 1, ntiles
      IF (is_dynamic(itile)) THEN
        !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) ASYNC(acc_stream)
        DO ic = 1, nc
          total_act_fpc(ic) = total_act_fpc(ic) + act_fpc(ic,itile)
        END DO
        !$ACC END PARALLEL LOOP
      END IF
    END DO

    !-- calculation of bare soil fraction

    !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) ASYNC(acc_stream)
    DO ic = 1, nc
      bare_fpc(ic) = MAX(0._wp, 1._wp - total_act_fpc(ic))
    END DO
    !$ACC END PARALLEL LOOP

    !$ACC END DATA

  END SUBROUTINE fpc_daily

! --- scale_fpc ---------------------------------------------------------------------------

  SUBROUTINE scale_fpc (calling_routine, nc, iblk, ntiles, is_dynamic, act_fpc, bare_fpc)
!
! Rescaling of fractional plant coverage to assure that
!  - the sum of act_fpc and bare soil is exactly one
!  - all tiles have a minimum vegetated fraction
!
! input
    CHARACTER, INTENT(in) :: calling_routine ! name of the calling routine
    INTEGER,  INTENT(in)  :: nc              ! vector length
    INTEGER,  INTENT(in)  :: iblk            ! number of block
    INTEGER,  INTENT(in)  :: ntiles          ! number of tiles (different land cover type mosaics per gridcell)
    LOGICAL,  INTENT(in)  :: is_dynamic(:)   ! PFTs taking part in competition

! inout
    REAL(wp), INTENT(inout) :: act_fpc(:,:) ! actual fpc in nlcc
    REAL(wp), INTENT(inout) :: bare_fpc(:)  ! fraction of bare ground

! local variables
    LOGICAL   :: just_return
    INTEGER   :: ic, itile      ! indices
    INTEGER   :: nsparce(nc)    ! number of PFTs with less then fract_small vegetation
    REAL(wp)  :: sum_fpc(nc)    ! sum of the different cover fractions
    REAL(wp)  :: excess(nc)     ! extra fraction, that needs to be redistibuted
    REAL(wp)  :: scalable(nc)   ! sum of fraction of tiles that can be scaled
    REAL(wp)  :: rescale(nc)    ! factor to rescale fraction of tiles

    CHARACTER(len=*), PARAMETER :: routine = modname//':scale_fpc'

    !$ACC DATA ASYNC(acc_stream) &
    !$ACC   CREATE(nsparce, sum_fpc, excess, scalable, rescale)

    ! test, if scaling is necessary
    just_return = .TRUE.

    !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) ASYNC(acc_stream)
    DO ic = 1, nc
      sum_fpc(ic) = bare_fpc(ic)
    END DO
    !$ACC END PARALLEL LOOP
    DO itile = 1, ntiles
      !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) ASYNC(acc_stream)
      DO ic = 1, nc
        sum_fpc(ic) = sum_fpc(ic) + act_fpc(ic,itile)
        IF (is_dynamic(itile) .AND. act_fpc(ic,itile) < fract_small) just_return = .FALSE.
      END DO
      !$ACC END PARALLEL LOOP
    END DO
    !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) ASYNC(acc_stream)
    DO ic = 1, nc
      IF (ABS(1._wp - sum_fpc(ic)) > 3._wp * EPSILON(1._wp)) just_return = .FALSE.
    END DO
    !$ACC END PARALLEL LOOP

    IF (just_return) RETURN

    IF (debug_on() .AND. iblk == 1) THEN
      WRITE (message_text,*) 'Attention: act_fpc is rescaled; called from ', calling_routine
      CALL message(routine, message_text)
    END IF

    ! sum act_fpc and bare_fpc
    !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) ASYNC(acc_stream)
    DO ic = 1, nc
      sum_fpc(ic) = bare_fpc(ic)
      nsparce(ic) = 0
    END DO
    !$ACC END PARALLEL LOOP

    DO itile = 1, ntiles
      IF (is_dynamic(itile)) THEN
        !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) ASYNC(acc_stream)
        DO ic = 1, nc
          IF (act_fpc(ic,itile) > fract_small) THEN
            sum_fpc(ic) = sum_fpc(ic) + act_fpc(ic,itile)
          ELSE
            nsparce(ic) = nsparce(ic) + 1
            act_fpc(ic,itile) = fract_small
          END IF
        END DO
        !$ACC END PARALLEL LOOP
      ELSE
        !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) ASYNC(acc_stream)
        DO ic = 1, nc
          act_fpc(ic,itile) = 0._wp
        END DO
        !$ACC END PARALLEL LOOP
      END IF
    END DO

    !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) ASYNC(acc_stream)
    DO ic = 1, nc
      sum_fpc(ic) = sum_fpc(ic) + REAL(nsparce(ic),wp) * fract_small
    END DO
    !$ACC END PARALLEL LOOP

    ! scaling of bare_fpc
    !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) ASYNC(acc_stream)
    DO ic = 1, nc
      IF (sum_fpc(ic) <= 1._wp - EPSILON(1._wp)) THEN
        bare_fpc(ic) = bare_fpc(ic) + (1._wp - sum_fpc(ic))
      ELSE IF (sum_fpc(ic) >= 1._wp + EPSILON(1._wp)) THEN
        bare_fpc(ic) = bare_fpc(ic) / sum_fpc(ic)
      END IF
    END DO
    !$ACC END PARALLEL LOOP

    ! scaling of act_fpc
    DO itile = 1, ntiles
      IF (is_dynamic(itile)) THEN
        !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) ASYNC(acc_stream)
        DO ic = 1, nc
          excess(ic) = 0._wp
          IF (sum_fpc(ic) >= 1._wp + EPSILON(1._wp) .AND. act_fpc(ic,itile) > fract_small * sum_fpc(ic)) THEN
            act_fpc(ic,itile) = act_fpc(ic,itile) / sum_fpc(ic)
          ELSE IF (sum_fpc(ic) >= 1._wp + EPSILON(1._wp)) THEN
            excess(ic) = excess(ic) + fract_small - act_fpc(ic,itile) / sum_fpc(ic)
            act_fpc(ic,itile) = fract_small
          END IF
        END DO
        !$ACC END PARALLEL LOOP
      END IF
    END DO

    !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) ASYNC(acc_stream)
    DO ic = 1, nc
      scalable(ic) = bare_fpc(ic)
    END DO
    !$ACC END PARALLEL LOOP

    DO itile = 1, ntiles
      IF (is_dynamic(itile)) THEN
        !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) ASYNC(acc_stream)
        DO ic = 1, nc
          IF (act_fpc(ic,itile) > 2._wp * fract_small) scalable = scalable + act_fpc(ic,itile)
        END DO
        !$ACC END PARALLEL LOOP
      END IF
    END DO
    !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) ASYNC(acc_stream)
    DO ic = 1, nc
      rescale(ic) = (scalable(ic) - excess(ic)) / scalable(ic)
      bare_fpc(ic) = bare_fpc(ic) * rescale(ic)
    END DO
    !$ACC END PARALLEL LOOP
    DO itile = 1, ntiles
      IF (is_dynamic(itile)) THEN
        !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) ASYNC(acc_stream)
        DO ic = 1, nc
          IF (act_fpc(ic,itile) > 2._wp * fract_small) act_fpc(ic,itile) = act_fpc(ic,itile) * rescale(ic)
        END DO
        !$ACC END PARALLEL LOOP
      END IF
    END DO

    !$ACC END DATA

  END SUBROUTINE scale_fpc


! --- calc_cover_fract_inst ---------------------------------------------------------------------------

  SUBROUTINE calc_cover_fract_inst(nc, ntiles, is_woody, is_dynamic, is_pasture, is_crop, &
                                   cover_fract_pot, cover_fract, cover_fract_inst)
!
! convert fractional plant cover calculated within nlcc to jsbach cover fractions
! that would result from an instantaneous establishment of crop land and pasture
! (pasture is preferentially established on grasslands)
!
! input
    INTEGER,  INTENT(in)  :: nc                ! vector length
    INTEGER,  INTENT(in)  :: ntiles            ! number of tiles (different land cover type mosaics per gridcell)
    LOGICAL,  INTENT(in)  :: is_woody(:)       ! PFTs that build wood
    LOGICAL,  INTENT(in)  :: is_dynamic(:)     ! PFTs taking part in competition
    LOGICAL,  INTENT(in)  :: is_pasture(:)     ! PFTs that are pasture
    LOGICAL,  INTENT(in)  :: is_crop(:)        ! PFTs that are crop
    REAL(wp), INTENT(in)  :: cover_fract_pot(:,:)  ! potential cover fractions (without land use)
    REAL(wp), INTENT(in)  :: cover_fract(:,:)  ! cover fraction within jsbach

! output
    REAL(wp), INTENT(out) :: cover_fract_inst(:,:) ! cover fraction within jsbach resulting from an
                                                   ! instantaneous establishment of crop land and pasture
! local variables
    INTEGER   :: ic, itile
    INTEGER   :: ngrass                  ! number of grass types
    INTEGER   :: ncrops                  ! number of crop types
    REAL(wp)  :: sum_grass_fract(nc)   ! potential fraction of grass if there was no land use
    REAL(wp)  :: sum_pasture_fract(nc) ! fraction of the vegetated part of the grid box covered with pasture
    REAL(wp)  :: sum_crop_fract(nc)    ! fraction of the vegetated part of the grid box covered with crops

    CHARACTER(len=*), PARAMETER :: routine = modname//':calc_cover_fract_inst'

    !$ACC DATA ASYNC(acc_stream) &
    !$ACC   CREATE(sum_grass_fract, sum_pasture_fract, sum_crop_fract)

    !-- sum up FPC of grasses, pastures and crops
    ngrass = 0
    ncrops = 0
    !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) ASYNC(acc_stream)
    DO ic = 1, nc
      sum_grass_fract(ic) = 0._wp
      sum_pasture_fract(ic) = 0._wp
      sum_crop_fract(ic) = 0._wp
    END DO
    !$ACC END PATALLEL LOOP
    DO itile = 1, ntiles
      IF (.NOT. is_woody(itile) .AND. is_dynamic(itile)) THEN
        ngrass = ngrass + 1
        !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) ASYNC(acc_stream)
        DO ic = 1, nc
          sum_grass_fract(ic) = sum_grass_fract(ic) + cover_fract_pot(ic,itile)
        END DO
        !$ACC END PATALLEL LOOP
      ELSE IF (is_pasture(itile)) THEN
        !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) ASYNC(acc_stream)
        DO ic = 1, nc
          sum_pasture_fract(ic) = sum_pasture_fract(ic) + cover_fract(ic,itile)
        END DO
        !$ACC END PATALLEL LOOP
      ELSE IF (is_crop(itile)) THEN
        ncrops = ncrops + 1
        !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) ASYNC(acc_stream)
        DO ic = 1, nc
          sum_crop_fract(ic) = sum_crop_fract(ic) + cover_fract(ic,itile)
        END DO
        !$ACC END PATALLEL LOOP
      END IF
    END DO

    !-- establishment of pastures (at the expense of grass land)
    DO itile = 1, ntiles
      IF (is_dynamic(itile)) THEN
        !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) ASYNC(acc_stream)
        DO ic = 1, nc
          IF (sum_pasture_fract(ic) <= sum_grass_fract(ic)) THEN
            ! There are enough grasslands
            IF (is_woody(itile)) THEN
               cover_fract_inst(ic,itile) = cover_fract_pot(ic,itile)
            ELSE
               cover_fract_inst(ic,itile) = cover_fract_pot(ic,itile) &
                 & * (sum_grass_fract(ic) - sum_pasture_fract(ic) + ngrass * fract_small) / sum_grass_fract(ic)
            END IF
          ELSE
            ! There are not enough grasslands for all pastures, so that woody types are also reduced
            IF (is_woody(itile)) THEN
              cover_fract_inst(ic,itile) = cover_fract_pot(ic,itile) &
                &  * (1._wp - sum_pasture_fract(ic)) / (1._wp - sum_grass_fract(ic))
            ELSE
              cover_fract_inst(ic,itile) = fract_small
            END IF
          END IF
        END DO
        !$ACC END PATALLEL LOOP
      ELSE IF (is_pasture(itile)) THEN
        !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) ASYNC(acc_stream)
        DO ic = 1, nc
          cover_fract_inst(ic,itile) = cover_fract(ic,itile)
        END DO
        !$ACC END PATALLEL LOOP
      END IF
    END DO

    !-- establishment of crops (at the expense of all natural pfts)
    DO itile = 1, ntiles
      IF (is_crop(itile)) THEN
        !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) ASYNC(acc_stream)
        DO ic = 1, nc
          cover_fract_inst(ic,itile) = cover_fract(ic,itile)
        END DO
        !$ACC END PATALLEL LOOP
      ELSE IF (.NOT. is_pasture(itile)) THEN
        !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) ASYNC(acc_stream)
        DO ic = 1, nc
          cover_fract_inst(ic,itile) = cover_fract_inst(ic,itile) &
            &  * (1._wp - sum_pasture_fract(ic) - sum_crop_fract(ic)) &
            &  / (1._wp - sum_pasture_fract(ic) - ncrops * fract_small)
        END DO
        !$ACC END PATALLEL LOOP
      END IF
    END DO

    CALL scale_cover_fract (nc, ntiles, is_dynamic(:), cover_fract_inst(:,:))

    !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) ASYNC(acc_stream)
    DO ic = 1, nc
      IF (  SUM(cover_fract_inst(ic,:)) > 1._wp + REAL(ntiles,wp)*EPSILON(1._wp) .OR. &
          & SUM(cover_fract_inst(ic,:)) < 1._wp - REAL(ntiles,wp)*EPSILON(1._wp)) THEN
        WRITE (message_text,*) 'Sum of cover_fract_inst /= 1: ',  &
          & SUM(cover_fract_inst(ic,:)), ' at ic: ', ic
        CALL finish (routine, message_text)
      END IF
    END DO
    !$ACC END PATALLEL LOOP
    DO itile = 1, ntiles
      !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) ASYNC(acc_stream)
      DO ic = 1, nc
        IF (cover_fract_inst(ic,itile) < fract_small) THEN
          WRITE (message_text,*) 'cover_fract_inst too small: ', cover_fract_inst(ic,itile), &
            &  ' at ic: ', ic
          CALL finish (routine, message_text)
        END IF
      END DO
      !$ACC END PATALLEL LOOP
    END DO

    !$ACC END DATA

  END SUBROUTINE calc_cover_fract_inst

!------------------------------------------------------------------------------

SUBROUTINE scale_cover_fract (nc, ntiles, is_dynamic, cover_fract)
!
! !DESCRIPTION:
!
! Rescaling of cover fractions to assure that
!  - the sum of cover fractions is one
!  - all tiles have at least a minimum vegetated fraction
!
!------------------------------------------------------------------------------
!
! !INPUT PARAMETERS:
!
    INTEGER,  INTENT(in)  :: nc                 ! vector length
    INTEGER,  INTENT(in)  :: ntiles             ! number of tiles
    LOGICAL,  INTENT(in)  :: is_dynamic(:)      ! mask for land surface treated by jsbach
!
! !IN- and OUTPUT PARAMETERS:
!
    REAL(wp), INTENT(inout) :: cover_fract(:,:) ! vegetated fraction

! !LOCAL VARIABLES:
!
    INTEGER   :: itile, ic, iter   ! looping indices
    INTEGER   :: niter             ! number of iterations needed
    INTEGER   :: nsparce(nc)     ! number of PFTs with a vegetated fraction of less then fract_small
    REAL(wp)  :: sum_fract(nc)   ! sum of all cover fractions
    REAL(wp)  :: excluded_fract(nc) ! sum of the PFT fractions that cannot be changed by natural vegetation
                                      ! dynamics: non-dynamic fractions and minimum fractions of dynamic PFTs
    REAL(wp)  :: anthro_fract(nc)   ! sum of non-dynamic PFT fractions
    REAL(wp)  :: considerable_nat_fract(nc) ! sum of natural PFTs with considerable fractions
    LOGICAL   :: corrected(nc)      ! flag indicating cells with corrected fractions

    CHARACTER(len=*), PARAMETER :: routine = modname//':scale_cover_fract'
!------------------------------------------------------------------------------

    !$ACC DATA ASYNC(acc_stream) &
    !$ACC   CREATE(nsparce, sum_fract, excluded_fract, anthro_fract, considerable_nat_fract, corrected)

    ! Make sure, crop and pasture have a cover fraction of at least fract_small
    DO itile = 1, ntiles
      IF (.NOT. is_dynamic(itile)) THEN
        !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) ASYNC(acc_stream)
        DO ic = 1, nc
          cover_fract(ic,itile) = MAX(fract_small,cover_fract(ic,itile))
        END DO
        !$ACC END PARALLEL LOOP
      END IF
    END DO

    ! Make sure, the sum of non-dynamic PFTs is smaller than (1 - fract_small*n_natural)
    !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) ASYNC(acc_stream)
    DO ic = 1, nc
      excluded_fract(ic) = 0._wp
    END DO
    !$ACC END PARALLEL LOOP

    DO itile = 1, ntiles
      IF (is_dynamic(itile)) THEN
        !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) ASYNC(acc_stream)
        DO ic = 1, nc
          excluded_fract(ic) = excluded_fract(ic) + fract_small
        END DO
        !$ACC END PARALLEL LOOP
      ELSE
        !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) ASYNC(acc_stream)
        DO ic = 1, nc
          excluded_fract(ic) = excluded_fract(ic) + cover_fract(ic,itile)
        END DO
        !$ACC END PARALLEL LOOP
      END IF
    END DO

    !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) ASYNC(acc_stream)
    DO ic = 1, nc
      IF (excluded_fract(ic) > 1._wp) THEN
        corrected(ic) = .FALSE.
        DO itile = 1, ntiles
          IF (excluded_fract(ic) > 1._wp .AND. .NOT. corrected(ic) .AND. .NOT. is_dynamic(itile) &
              &  .AND. cover_fract(ic,itile) > (fract_small + excluded_fract(ic)-1._wp)) THEN
            cover_fract(ic,itile) = cover_fract(ic,itile) - (excluded_fract(ic)-1._wp)
            corrected(ic) = .TRUE.
          END IF
        END DO
      END IF
    END DO
    !$ACC END PARALLEL LOOP

    ! The more tiles we have, the more iterations are needed. For binary identical results
    ! the number of iterations must not change.
    niter = ntiles
    DO iter = 1, niter

      !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) ASYNC(acc_stream)
      DO ic = 1, nc
        sum_fract(ic) = 0._wp
        considerable_nat_fract(ic) = 0._wp
        anthro_fract(ic) = 0._wp
        nsparce(ic) = 0
      END DO
      !$ACC END PARALLEL LOOP

      DO itile = 1, ntiles
        !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) ASYNC(acc_stream)
        DO ic = 1, nc
          sum_fract(ic) = sum_fract(ic) + cover_fract(ic,itile)
        END DO
        !$ACC END PARALLEL LOOP

        IF (is_dynamic(itile)) THEN
          !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) ASYNC(acc_stream)
          DO ic = 1, nc
            IF (cover_fract(ic,itile) > fract_small) THEN
              considerable_nat_fract(ic) = considerable_nat_fract(ic) + cover_fract(ic,itile)
            ELSE
              nsparce(ic) = nsparce(ic) + 1
            END IF
          END DO
          !$ACC END PARALLEL LOOP
        ELSE
          !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) ASYNC(acc_stream)
          DO ic = 1, nc
            anthro_fract(ic) = anthro_fract(ic) + cover_fract(ic,itile)
          END DO
          !$ACC END PARALLEL LOOP
        END IF
      END DO

      DO itile = 1, ntiles
        IF (is_dynamic(itile)) THEN
          !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) ASYNC(acc_stream)
          DO ic = 1, nc
            IF (cover_fract(ic,itile) > fract_small) THEN
              cover_fract(ic,itile) = cover_fract(ic,itile) * &
                &                    (1._wp - anthro_fract(ic) - REAL(nsparce(ic),wp)*fract_small) / considerable_nat_fract(ic)
            ELSE
              cover_fract(ic,itile) = fract_small
            END IF
          END DO
          !$ACC END PARALLEL LOOP
        ELSE
          !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) ASYNC(acc_stream)
          DO ic = 1, nc
            cover_fract(ic,itile) = MAX(fract_small, cover_fract(ic,itile))
          END DO
          !$ACC END PARALLEL LOOP
        END IF
      END DO
    END DO  ! iter

    ! Final corrections, to assure the sum of all PFT fractions equals 1. This sometimes might not be
    ! the case for numerical reasons.

    !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) ASYNC(acc_stream)
    DO ic = 1, nc
      sum_fract(ic) = SUM(cover_fract(ic,:))
      ! Exit if the cover fraction sum is considerably different from 1.
      IF (ABS(1._wp - sum_fract(ic)) > fract_small) THEN
        WRITE(message_text,*) 'Sum of PFT fractions differs from 1.: ', sum_fract(ic), ' at ic: ', ic
        CALL finish (routine, message_text)
      END IF
    END DO
    !$ACC END PARALLEL LOOP

    !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) ASYNC(acc_stream)
    DO ic = 1, nc
      corrected(ic) = .FALSE.
      DO itile = 1, ntiles
        IF (is_dynamic(itile) .AND. (sum_fract(ic) < 1._wp - EPSILON(1._wp)) &
            &  .AND. .NOT. corrected(ic) .AND. cover_fract(ic,itile) < sum_fract(ic)) THEN
          cover_fract(ic,itile) = cover_fract(ic,itile) + (1._wp - sum_fract(ic))
          corrected(ic) = .TRUE.
        END IF
      END DO
    END DO
    !$ACC END PARALLEL LOOP

    !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) ASYNC(acc_stream)
    DO ic = 1, nc
      corrected(ic) = .FALSE.
      DO itile = 1, ntiles
        IF (is_dynamic(itile) .AND. (sum_fract(ic) > 1._wp + EPSILON(1._wp)) &
            &  .AND. .NOT. corrected(ic) .AND. cover_fract(ic,itile) > 1._wp - sum_fract(ic) - fract_small) THEN
          cover_fract(ic,itile) = cover_fract(ic,itile) - (sum_fract(ic) - 1._wp)
          corrected(ic) = .TRUE.
        END IF
      END DO
    END DO
    !$ACC END PARALLEL LOOP

    DO itile = 1, ntiles
      !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) ASYNC(acc_stream)
      DO ic = 1, nc
        IF (cover_fract(ic,itile) < fract_small) THEN
          WRITE(message_text,*) 'cover_fract still smaller ', fract_small, ' after ', niter, ' iterations:', &
            &  ' cover_fract(ic,itile): ', cover_fract(ic,itile), ' at ic: ', ic, ' itile: ', itile
          CALL message (routine, message_text)

          ! Sometimes cover_fract remains slightly smaller than fract_small for numerical reasons
          IF (cover_fract(ic,itile) < fract_small &
              & .AND. cover_fract(ic,itile) >= fract_small - REAL(ntiles,wp)*EPSILON(1._wp)) THEN
            cover_fract(ic,itile) = fract_small
          END IF
          IF (cover_fract(ic,itile) < fract_small) THEN
            WRITE(message_text,*) 'cover_fract still smaller ', fract_small, ' after ', niter, ' iterations:', &
              &  ' cover_fract(ic,itile): ', cover_fract(ic,itile), ' at ic: ', ic, ' itile: ', itile
            CALL finish (routine, message_text)
          END IF
        END IF
      END DO
      !$ACC END PARALLEL LOOP
    END DO

    !$ACC END DATA

  END SUBROUTINE scale_cover_fract

#endif
END MODULE mo_nlcc_process
