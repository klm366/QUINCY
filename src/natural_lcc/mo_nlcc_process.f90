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
  USE mo_kind,        ONLY: jp
  USE mo_exception,   ONLY: message, finish, message_text
  USE mo_jsb_physical_constants, ONLY: tmelt

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: bioclim_limits, scale_fpc, potential_tree_fpc, calc_cover_fract_inst
  PUBLIC :: fpc_to_cover_fract_pot, cover_fract_pot_to_cover_fract, desert_fraction, fpc_daily, fract_small

  CHARACTER(len=*), PARAMETER :: modname = 'mo_nlcc_process'
  REAL(jp), PARAMETER :: fract_small = 1.e-10_jp   !! smallest fraction of PFT tiles
  REAL(jp), PARAMETER :: act_fpc_min = 0.005_jp    !! minimum actual FPC to be considered to calculate potential FPC
  REAL(jp), PARAMETER :: sum_npp_min = 1.e-12_jp   !! minimum total npp for tree cover
  REAL(jp), PARAMETER :: tree_fpc_max = 1.0_jp     !! maximum FPC of woody PFTs
  REAL(jp), PARAMETER :: npp_nonlinearity = 1.5_jp !! parameter controlling the non-linearity in the dynamic equation with
                                                   !! respect to NPP
  REAL(jp), PARAMETER :: desert_extend = 0.65_jp   !! parameter controlling the extend of desert
                                                   !! (the lower the value the more desert)
  REAL(jp), PARAMETER :: desert_margin = 2._jp     !! parameter controlling the transition from vegetated land to desert
                                                   !! (the higher the value the sharper the transition)
  REAL(jp), PARAMETER :: tau_desert = 50._jp       !! time constant [year] by which the extension of deserts is adapted

CONTAINS

! --- bioclim_limits ---------------------------------------------------------------------------------------

  ELEMENTAL PURE SUBROUTINE bioclim_limits (is_dynamic, lct_tcmin, lct_tcmax, lct_twmax, lct_min_temprange, lct_gddmin, &
                                            min_mmtemp20, max_mmtemp20, gdd_prev_year, bio_exist)

!
! calculation of bioclimatic limits for each PFT based on modified
! LPJ lookup-table
!

! input
    LOGICAL,  INTENT(in)  :: is_dynamic        ! mask indicating whether a PFT is subject to vegetation dynamics
    REAL(jp), INTENT(in)  :: lct_tcmin         ! PFT-specific minimum coldest-month temperature limit [C]
    REAL(jp), INTENT(in)  :: lct_tcmax         ! PFT-specific maximum coldest-month temperature limit [C]
    REAL(jp), INTENT(in)  :: lct_twmax         ! PFT-specific maximum warmest-month temperature limit [C]
    REAL(jp), INTENT(in)  :: lct_min_temprange ! PFT-specific minimum difference of 20-year average warme st
                                               ! minus coldest month temperature [C]
    REAL(jp), INTENT(in)  :: lct_gddmin        ! PFT-specific minimum GDD limit
    REAL(jp), INTENT(in)  :: min_mmtemp20      ! Minimum monthly mean temp. (20yr climatology) [K]
    REAL(jp), INTENT(in)  :: max_mmtemp20      ! Maximum monthly mean temp. (20yr climatology) [K]
    REAL(jp), INTENT(in)  :: gdd_prev_year     ! GDD of previous year

! output
    REAL(jp), INTENT(inout) :: bio_exist

!------------------------------------------------------------------------------

    bio_exist = 1._jp

    IF (is_dynamic .AND. (min_mmtemp20 - tmelt < lct_tcmin .OR. min_mmtemp20 - tmelt > lct_tcmax  &
       .OR. gdd_prev_year < lct_gddmin .OR. max_mmtemp20 - tmelt > lct_twmax                      &
       .OR. (max_mmtemp20 - min_mmtemp20) < lct_min_temprange)) bio_exist = 0._jp

  END SUBROUTINE bioclim_limits


! --- potential_tree_fpc ---------------------------------------------------------------------------

  SUBROUTINE potential_tree_fpc(nidx, ntiles, is_dynamic, is_woody, bio_exist, &
                                npp_ave, act_fpc, pot_fpc)

!
! subroutine calculates potential FPC (in absence of disturbances) based on NPP for each PFT
! called once a year!
!

! input
    INTEGER,  INTENT(in)  :: nidx            ! vector length
    INTEGER,  INTENT(in)  :: ntiles          ! number of tiles (different land cover type mosaics per gridcell)
    LOGICAL,  INTENT(in)  :: is_dynamic(:,:) ! flag indicating dynamic PFTs
    LOGICAL,  INTENT(in)  :: is_woody(:,:)   ! flag indicating woody PFTs
    REAL(jp), INTENT(in)  :: bio_exist(:,:)  ! bio_exist=1. indicates that all bio-climatic limits are fulfilled
    REAL(jp), INTENT(in)  :: npp_ave(:,:)    ! NPP averaged over 5 years
    REAL(jp), INTENT(in)  :: act_fpc(:,:)    ! actual foliage projective cover

! output
    REAL(jp), INTENT(out) :: pot_fpc(:,:)

! local variables
    INTEGER           :: i, itile
    REAL(jp)          :: sum_npp(nidx)

! summing up NPP
    sum_npp(:) = 0._jp

    DO i = 1,nidx
       DO itile = 1,ntiles
          IF (is_woody(i,itile) .AND. is_dynamic(i,itile) .AND. &
           npp_ave(i,itile) > sum_npp_min .AND. bio_exist(i,itile) > 0.5_jp) &
           sum_npp(i) = sum_npp(i) + npp_ave(i,itile) ** npp_nonlinearity * MAX(act_fpc_min,act_fpc(i,itile))
       END DO
    END DO

!-- determine potential FPC
    DO i = 1,nidx
       DO itile = 1,ntiles
          IF (is_woody(i,itile) .AND. is_dynamic(i,itile) .AND. npp_ave(i,itile) > sum_npp_min .AND. &
              bio_exist(i,itile) > 0.5_jp) THEN
             pot_fpc(i,itile) = npp_ave(i,itile) ** npp_nonlinearity / sum_npp(i) * &
                                tree_fpc_max * MAX(act_fpc_min,act_fpc(i,itile))
          ELSE
             pot_fpc(i,itile) = 0._jp
          END IF
       END DO
    END DO

  END SUBROUTINE potential_tree_fpc


! --- fpc_to_cover_fract_pot -------------------------------------------------------------------------

  SUBROUTINE fpc_to_cover_fract_pot (nidx, ntiles, is_dynamic, is_woody, act_fpc, bare_fpc, cover_fract_pot)
!
! Convert fractional plant cover (act_fpc) calculated within nlcc to jsbach
! cover fractions reflecting natural vegetation, only (cover_fract_pot).
! The only difference between act_fpc and cover_fract_pot is caused by bare_fpc,
! which is added to the grass cover fractions.
!

! input
    INTEGER,  INTENT(in)  :: nidx            ! vector length
    INTEGER,  INTENT(in)  :: ntiles          ! number of tiles (different land cover type mosaics per gridcell)
    LOGICAL,  INTENT(in)  :: is_dynamic(:,:) !
    LOGICAL,  INTENT(in)  :: is_woody(:,:)   !
    REAL(jp), INTENT(in)  :: act_fpc(:,:)    ! actual fpc in nlcc
    REAL(jp), INTENT(in)  :: bare_fpc(:)     ! nlcc bare soil fraction
! output
    REAL(jp), INTENT(out) :: cover_fract_pot(:,:) ! cover fractions if there was only natural vegetation

! local variables
    INTEGER   :: i, itile
    REAL(jp)  :: sum_grass_fpc(nidx)     ! part of the vegetated area covered by grass

! -- initialization
    cover_fract_pot(:,:) = 0._jp

! -- find out fpc of grass
    sum_grass_fpc(:) = 0._jp

    DO i = 1,nidx
      DO itile = 1,ntiles
         IF (.NOT. is_woody(i,itile) .AND. is_dynamic(i,itile)) &
            sum_grass_fpc(i) = sum_grass_fpc(i) + act_fpc(i,itile)
      END DO
    END DO

!-- calculate new cover_fractions
    DO i = 1,nidx
       DO itile = 1,ntiles
          IF (is_dynamic(i,itile)) THEN
             IF (is_woody(i,itile)) THEN
                cover_fract_pot(i,itile) = act_fpc(i,itile)
             ELSE
                cover_fract_pot(i,itile) = act_fpc(i,itile) * (1._jp + bare_fpc(i) / MAX(sum_grass_fpc(i),EPSILON(1._jp)))
             END IF
          END IF
       END DO
    END DO

    CALL scale_cover_fract (nidx, ntiles, is_dynamic(:,:), cover_fract_pot(:,:))

    IF (ANY(SUM(cover_fract_pot(:,:),DIM=2) > 1._jp + REAL(ntiles,jp)*EPSILON(1._jp)) .OR. &
        ANY(SUM(cover_fract_pot(:,:),DIM=2) < 1._jp - REAL(ntiles,jp)*EPSILON(1._jp))) THEN
       WRITE (message_text,*) 'sum of cover_fract_pot /= 1: ', &
            MINVAL(SUM(cover_fract_pot(:,:),DIM=2)), MAXVAL(SUM(cover_fract_pot(:,:),DIM=2)), &
            MAXLOC(SUM(cover_fract_pot(:,:),DIM=2))
       CALL finish ('fpc_to_cover_fract_pot', message_text)
    END IF
    DO itile = 1,ntiles
       IF (ANY(is_dynamic(:,itile) .AND. cover_fract_pot(:,itile) < fract_small)) THEN
          WRITE (message_text,*) 'cover_fract_pot too small: ', MINVAL(cover_fract_pot(:,itile)), itile
          CALL finish ('fpc_to_cover_fract_pot', message_text)
       END IF
    END DO

  END SUBROUTINE fpc_to_cover_fract_pot


! --- cover_fract_pot_to_cover_fract ---------------------------------------------------------------------------

  SUBROUTINE cover_fract_pot_to_cover_fract (nidx, ntiles, is_dynamic, is_woody, is_pasture, is_crop, &
                                             cover_fract_pot, cover_fract_pot_previous_year, cover_fract)

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
    INTEGER,  INTENT(in)  :: nidx                 ! vector length
    INTEGER,  INTENT(in)  :: ntiles               ! number of tiles (different land cover type mosaics per gridcell)
    LOGICAL,  INTENT(in)  :: is_dynamic(:,:)      ! PFTs taking part in competition
    LOGICAL,  INTENT(in)  :: is_woody(:,:)        ! PFTs that build wood
    LOGICAL,  INTENT(in)  :: is_pasture(:,:)      ! logical mask for pasture
    LOGICAL,  INTENT(in)  :: is_crop(:,:)         ! logical mask for crops
    REAL(jp), INTENT(in)  :: cover_fract_pot(:,:) ! cover fractions if there was only natural vegetation
    REAL(jp), INTENT(in)  :: cover_fract_pot_previous_year(:,:) ! last years cover_fract_pot

! inout
    REAL(jp), INTENT(inout) :: cover_fract(:,:)  ! cover fraction within jsbach

! local variables

    INTEGER   :: i, itile
    INTEGER   :: nwoody(nidx)                   ! number of woody PFTs
    REAL(jp)  :: nexcluded(nidx)                ! number of non-dynamic PFTs (crops, pastures)
    REAL(jp)  :: cover_fract_inst(nidx,ntiles)  ! cover fractions that would result from an instantanous implementation of
                                                ! land use change and vegetation dynamics
    REAL(jp)  :: excluded_fract(nidx)           ! fraction of the vegetated part of the grid box not considered for vegetation
                                                ! dynamics (e.g. agricultural areas)
    REAL(jp)  :: sum_woody_fract(nidx)          ! fraction of the vegetated part of the grid box covered by woody plants
    REAL(jp)  :: sum_grass_fract(nidx)          ! fraction of the vegetated part of the grid box covered by grass
    REAL(jp)  :: sum_woody_fract_inst(nidx)     ! fraction of the vegetated part of the grid box covered by woody plants
                                                ! based on cover_fract_inst
    REAL(jp)  :: sum_woody_fract_pot(nidx)      ! woody fraction of the vegetated part of the grid box considering natural
                                                ! vegetation, only
    REAL(jp)  :: sum_grass_fract_pot(nidx)      ! grass fraction of the vegetated part of the grid box considering natural
                                                ! vegetation, only
    REAL(jp)  :: sum_woody_fract_old(nidx)      ! fraction of the vegetated part of the grid box covered by woody plants
                                                ! (based on act_fpc_previous_year) that would have resulted from an
                                                ! instantaneous establishment of crop land and pasture
    REAL(jp)  :: sum_woody_fract_act(nidx)      ! fraction of the vegetated part of the grid box covered by woody plants
    REAL(jp)  :: sum_grass_fract_act(nidx)      ! fraction of the vegetated part of the grid box covered by grasses
    REAL(jp)  :: delta_woody(nidx)              ! deviation of the fraction of the vegetated part of the grid box covered
                                                ! by woody plants from the fraction that would result from an instantaneous
                                                !  establishment of crop land and pasture

    !-- find out number of woody types and number of exculded types
    nwoody(:) = 0
    nexcluded(:) = 0._jp
    DO i=1,nidx
       DO itile = 1,ntiles
          IF (is_dynamic(i,itile) .AND. is_woody(i,itile)) nwoody(i) = nwoody(i) + 1
          IF (.NOT. is_dynamic(i,itile)) nexcluded(i) = nexcluded(i) + 1._jp
       ENDDO
    ENDDO

    !-- Sum up cover fractions of woody types, grasses and fractions not affected by NLCC.
    sum_woody_fract(:) = 0._jp
    sum_grass_fract(:) = 0._jp
    excluded_fract(:) = 0._jp

    DO i = 1,nidx
       DO itile = 1,ntiles
          IF (is_woody(i,itile) .AND. is_dynamic(i,itile)) THEN
             sum_woody_fract(i) = sum_woody_fract(i) + cover_fract(i,itile)
          ELSE IF (.NOT. is_woody(i,itile) .AND. is_dynamic(i,itile)) THEN
             sum_grass_fract(i) = sum_grass_fract(i) + cover_fract(i,itile)
          ELSE IF (.NOT. is_dynamic(i,itile)) THEN
             excluded_fract(i) = excluded_fract(i) + cover_fract(i,itile)
          END IF
       END DO
    END DO

    !-- Sum up woody and grass fractions of the potential natural vegetation based on current vegetation dynamics
    sum_woody_fract_pot(:) = 0._jp
    sum_grass_fract_pot(:) = 0._jp

    DO i = 1,nidx
       DO itile = 1,ntiles
          IF (is_woody(i,itile) .AND. is_dynamic(i,itile)) THEN
             sum_woody_fract_pot(i) = sum_woody_fract_pot(i) + cover_fract_pot(i,itile)
          ELSE IF (.NOT. is_woody(i,itile) .AND. is_dynamic(i,itile)) THEN
             sum_grass_fract_pot(i) = sum_grass_fract_pot(i) + cover_fract_pot(i,itile)
          END IF
       END DO
    END DO

    !-- Calculate the woody fraction resulting from last years potential vegetation assuming an instantaneous
    !   establishment of crops and pastures
    CALL calc_cover_fract_inst(nidx, ntiles, is_woody(:,:), is_dynamic(:,:), is_pasture(:,:), is_crop(:,:), &
                               cover_fract_pot_previous_year(:,:), cover_fract(:,:), cover_fract_inst(:,:))

    sum_woody_fract_old(:) = 0._jp
    DO i = 1,nidx
       DO itile = 1,ntiles
          IF (is_woody(i,itile) .AND. is_dynamic(i,itile)) &
               sum_woody_fract_old(i) = sum_woody_fract_old(i) + cover_fract_inst(i,itile)
       END DO
    END DO

    !-- Calculate the woody fraction resulting from current potential vegetation assuming an instantaneous
    !   establishment of crops and pastures
    CALL calc_cover_fract_inst(nidx, ntiles, is_woody(:,:), is_dynamic(:,:), is_pasture(:,:), is_crop(:,:), &
                               cover_fract_pot(:,:), cover_fract(:,:), cover_fract_inst(:,:))

    sum_woody_fract_inst(:) = 0._jp
    DO i = 1,nidx
       DO itile = 1,ntiles
          IF (is_woody(i,itile) .AND. is_dynamic(i,itile)) THEN
             sum_woody_fract_inst(i) = sum_woody_fract_inst(i) + cover_fract_inst(i,itile)
          END IF
       END DO
    END DO

    !-- Calculation of the final woody fraction. A difference between sum_woody_fract and sum_woody_fract_old indicates
    !   an imbalance of the state of the dynamic vegetation and the cover fractions calculated from landuse transitions.
    !   On the other hand, the difference between sum_woody_fract_inst and sum_woody_fract_old show the trend of the woody
    !   fraction caused by vegetation dynamics.
    !   The idea of this routine is, to allow shifts in forest fractions only if they reduce the imbalance between nlcc
    !   and landuse change.
    DO i = 1,nidx
       IF (sum_woody_fract(i) - sum_woody_fract_old(i) > 0._jp) THEN
          ! -- woody fraction of cover_fract entering the routine (last updated by last years land use change) is greater
          !    than the the woody fraction resulting from last years potential vegetation and instantanous establishment of
          !    crops and pastures. Thus there is an inbalance, the woody fraction of cover_fract is too high.
          IF (sum_woody_fract_inst(i) - sum_woody_fract_old(i) > 0._jp) THEN
             ! -- vegetation dynamics even enlarge woody fraction
             !    ==> no change in actual woody fraction, unless the enlargement by vegetation dynamics exceeds the initial
             !        woody surplus
             delta_woody(i) = MAX(0._jp, sum_woody_fract_inst(i) - sum_woody_fract(i))
          ELSE
             ! -- vegetation dynamics shrink woody fraction
             delta_woody(i) = sum_woody_fract_inst(i) - sum_woody_fract_old(i)
          END IF
       ELSE
          ! -- The woody fraction of cover_fract is too small
          IF (sum_woody_fract_inst(i) - sum_woody_fract_old(i) > 0._jp) THEN
             ! -- vegetation dynamics enlarge woody fraction
             delta_woody(i) = sum_woody_fract_inst(i) - sum_woody_fract_old(i)
          ELSE
             ! -- vegetation dynamics shrink woody fraction
             !    ==> no change in actual woody fraction, unless the reduction by vegetation dynamics exceeds the initial
             !        woody deficit
             delta_woody(i) = MIN(0._jp, sum_woody_fract_inst(i) - sum_woody_fract(i))
          END IF
       END IF

       ! -- in some cases delta_woody can not be used to calculate the actual woody fraction. If there is no grass left
       !    after instantanous establishement of agricultural areas, the woody fractions of cover_fract_inst and
       !    cover_fract_inst_old are identical and delta_woody is zero, even if vegetation dynamics change considerably.
       !    The following formulation is needed to assure that sum_woody_fract_act is smaller than sum_woody_fract_pot and
       !    sum_grass_fract_act is smaller than sum_grass_fract_pot in all cases.
       sum_woody_fract_act(i) = MIN(sum_woody_fract_pot(i),&
                                MAX(nwoody(i)*fract_small, sum_woody_fract_pot(i)-excluded_fract(i)+nexcluded(i)*fract_small, &
                                    sum_woody_fract(i) + delta_woody(i)))
    END DO

    sum_grass_fract_act(:) = 1._jp - sum_woody_fract_act(:) - excluded_fract(:)

    !-- Use actual woody fraction to determine all actual fractions
    DO i = 1,nidx
       DO itile = 1,ntiles
          IF (is_dynamic(i,itile)) THEN
             IF (is_woody(i,itile)) THEN
                cover_fract(i,itile) = cover_fract_pot(i,itile) * sum_woody_fract_act(i) / sum_woody_fract_pot(i)
             ELSE
                cover_fract(i,itile) = cover_fract_pot(i,itile) * sum_grass_fract_act(i) / sum_grass_fract_pot(i)
             END IF
          END IF
       END DO
    END DO

    CALL scale_cover_fract(nidx, ntiles, is_dynamic(:,:), cover_fract(:,:))

    !-- Stop simulation, if sum of cover fractions is not equal to 1 or a cover fraction is smaller than fract_small
    IF (  ANY(SUM(cover_fract(:,:),DIM=2) > 1._jp + REAL(ntiles,jp)*EPSILON(1._jp)) .OR.   &
       &  ANY(SUM(cover_fract(:,:),DIM=2) < 1._jp - REAL(ntiles,jp)*EPSILON(1._jp))) THEN
       WRITE (message_text,*) 'sum of cover_fract /= 1: ', &
         &   MINVAL(SUM(cover_fract(:,:),DIM=2)), MAXVAL(SUM(cover_fract(:,:),DIM=2)), ', MINLOC: ', &
         &   MINLOC(SUM(cover_fract(:,:),DIM=2)), ', MAXLOC: ', MAXLOC(SUM(cover_fract(:,:),DIM=2))
       CALL finish ('cover_fract_pot_to_cover_fract', message_text)
    END IF
    DO itile = 1,ntiles
       IF (ANY(cover_fract(:,itile) < fract_small)) THEN
          WRITE (message_text,*) 'cover_fract too small: ', MINVAL(cover_fract(:,itile))
          CALL finish ('cover_fract_pot_to_cover_fract', message_text)
       END IF
    END DO

  END SUBROUTINE cover_fract_pot_to_cover_fract


! --- desert_fraction --------------------------------------------------------

  SUBROUTINE desert_fraction(nidx, ntiles, init_running_means, accelerate_nlcc,   &
                             is_woody, is_dynamic,                                &
                             act_fpc, bare_fpc,                                   &
                             sla, max_green_bio,                                  &
                             sum_green_bio_memory, desert_fpc)
! !DESCRIPTION:
!
! Calculation of desert fraction from green biomass
!
!------------------------------------------------------------------------------

! input
    INTEGER,  INTENT(in)  :: nidx              ! vector length
    INTEGER,  INTENT(in)  :: ntiles            ! number of tiles (different land cover type mosaics per gridcell)
    LOGICAL,  INTENT(in)  :: init_running_means! initialisation of sum_green_bio_memory
    REAL(jp), INTENT(in)  :: accelerate_nlcc   ! acceleration factor for vegetation and desert dynamics
    LOGICAL,  INTENT(in)  :: is_woody(:,:)     ! flag to label a woody PFT
    LOGICAL,  INTENT(in)  :: is_dynamic(:,:)   ! flag to label PFTs that take part in the natural vegetation dynamics
    REAL(jp), INTENT(in)  :: act_fpc(:,:)      ! actual FPC
    REAL(jp), INTENT(in)  :: bare_fpc(:)       ! bare FPC
    REAL(jp), INTENT(in)  :: sla(:,:)          ! specific leaf area
    REAL(jp), INTENT(in)  :: max_green_bio(:,:)! maximum value of green biomass within a year

! input/output
    REAL(jp), INTENT(inout) :: sum_green_bio_memory(:) ! vegetated fraction calculated from green biomass

! output
    REAL(jp), INTENT(out) :: desert_fpc(:) ! desert FPC

! local variables
    INTEGER           :: i, itile
    REAL(jp)          :: sum_green_bio(nidx)
    REAL(jp)          :: sum_act_fpc(nidx)
    REAL(jp)          :: sum_grass_fpc(nidx)
    INTEGER           :: n_grass_pft(nidx)

!------------------------------------------------------------------------------

    sum_green_bio(:) = 0._jp
    sum_act_fpc(:) = 0._jp
    sum_grass_fpc(:) = 0._jp
    n_grass_pft(:) = 0

    DO i = 1,nidx
       DO itile = 1,ntiles
          IF (.NOT. is_woody(i,itile) .AND. is_dynamic(i,itile)) THEN
             n_grass_pft(i) = n_grass_pft(i) + 1
             sum_grass_fpc(i) = sum_grass_fpc(i) + act_fpc(i,itile)
          END IF
       END DO
    END DO

    DO i = 1,nidx
       DO itile = 1,ntiles
          IF (is_woody(i,itile) .AND. is_dynamic(i,itile)) THEN
             sum_green_bio(i) = sum_green_bio(i) + MAX(0._jp, act_fpc(i,itile) * (1._jp -          &
               exp(-desert_extend * ((max_green_bio(i,itile) * sla(i,itile)) ** desert_margin) /   &
               (3._jp ** (desert_margin - 1._jp)))))
             sum_act_fpc(i) = sum_act_fpc(i) + act_fpc(i,itile)
          ELSE IF (.NOT. is_woody(i,itile) .AND. is_dynamic(i,itile)) THEN
             IF (sum_grass_fpc(i) > EPSILON(1._jp)) THEN
                sum_green_bio(i) = sum_green_bio(i) +                                              &
                 MAX(0._jp, act_fpc(i,itile) * (1._jp + bare_fpc(i) / sum_grass_fpc(i)) * (1._jp - &
                 exp(-desert_extend * ((max_green_bio(i,itile) * sla(i,itile)) ** desert_margin) / &
                (3._jp ** (desert_margin - 1._jp)))))
                sum_act_fpc(i) = sum_act_fpc(i) + act_fpc(i,itile) * (1._jp + bare_fpc(i) / sum_grass_fpc(i))
             ELSE
                sum_green_bio(i) = sum_green_bio(i) +                                              &
                 MAX(0._jp, (bare_fpc(i) / REAL(n_grass_pft(i),jp)) * (1._jp -                     &
                 exp(-desert_extend * ((max_green_bio(i,itile) * sla(i,itile)) ** desert_margin) / &
                 (3._jp ** (desert_margin - 1._jp)))))
                sum_act_fpc(i) = sum_act_fpc(i) + (bare_fpc(i) / REAL(n_grass_pft(i),jp))
             END IF
          END IF
       END DO
    END DO

    WHERE (sum_act_fpc(:) > EPSILON(1._jp))
       sum_green_bio(:) = sum_green_bio(:) / sum_act_fpc(:)
    END WHERE
    IF (init_running_means) THEN
       sum_green_bio_memory(:) = sum_green_bio(:)
    ELSE
       sum_green_bio_memory(:) = (sum_green_bio_memory(:) * (tau_desert/accelerate_nlcc - 1._jp) + sum_green_bio(:)) &
                                 / (tau_desert/accelerate_nlcc)
    END IF

    desert_fpc(:) = MIN(1._jp - 2._jp * EPSILON(1._jp),MAX(0._jp,1._jp - sum_green_bio_memory(:)))

  END SUBROUTINE desert_fraction

! --- fpc_daily ---------------------------------------------------------------------------

  SUBROUTINE fpc_daily(nidx, ntiles, accelerate_nlcc, is_woody, is_dynamic, tau_Cpool_woods,    &
                       npp_ave, bio_exist, burned_fract, damaged_fract, pot_fpc, &
                       act_fpc, bare_fpc, fuel)

! Calculation of daily FPC for trees, shrubs, and herbaceous plants

! input
    INTEGER,         INTENT(in) :: nidx            ! vector length
    INTEGER,         INTENT(in) :: ntiles          ! number of tiles (different land cover type mosaics per gridcell)
    REAL(jp),        INTENT(in) :: accelerate_nlcc ! factor to accellerate vegetation and desert dynamics
    LOGICAL,         INTENT(in) :: is_woody(:,:)
    LOGICAL,         INTENT(in) :: is_dynamic(:,:)
    REAL(jp),        INTENT(in) :: tau_Cpool_woods(:,:) ! lifetime of Cpool_woods [days]
    REAL(jp),        INTENT(in) :: npp_ave(:,:)
    REAL(jp),        INTENT(in) :: bio_exist(:,:)
    REAL(jp),        INTENT(in) :: burned_fract(:,:)
    REAL(jp),        INTENT(in) :: damaged_fract(:,:)
    REAL(jp),        INTENT(in) :: pot_fpc(:,:)

! inout
    REAL(jp), INTENT(inout) :: act_fpc(:,:)
    REAL(jp), INTENT(inout) :: bare_fpc(:)
    REAL(jp), INTENT(inout), optional :: fuel(:)

! local variables
    INTEGER     :: i, j, itile
    REAL(jp)    :: act_fpc_previous_day(nidx,ntiles)
    REAL(jp)    :: total_act_fpc(nidx)
    REAL(jp)    :: grass_sum(nidx)
    REAL(jp)    :: non_woody_fpc(nidx)
    REAL(jp)    :: woody_estab_fpc(nidx)
    REAL(jp)    :: theta, veg_inc, tau_pft
    REAL(jp), PARAMETER  ::  min_dist_woody = 0.002_jp / 365._jp !! mortality of woody plants (not caused by fire/windthrow,
                                                                 !! but by other unresolved processes)
    REAL(jp), PARAMETER  ::  grass_mortality = 0.01_jp / 365._jp !! grass mortality (not caused by fire, but by other
                                                                 !! unresolved processes)
    REAL(jp), PARAMETER  ::  tau_Cpool_woods_to_tau_pft = 1._jp

    ! initialisations

    act_fpc_previous_day(:,:) = act_fpc(:,:)
    non_woody_fpc(:) = bare_fpc(:)

    DO itile = 1,ntiles
       DO i = 1,nidx
!!$ TR          IF (is_dynamic(i,itile) .AND. is_woody(i,itile)) burned_local(i,itile) = 3.e-5_jp
!!$ TR          IF (is_dynamic(i,itile) .AND. .NOT. is_woody(i,itile)) burned_local(i,itile) = 1.e-4_jp
!!$ TR          IF (is_dynamic(i,itile) .AND. is_woody(i,itile)) damaged_local(i,itile) = 1.e-5_jp
          IF (.NOT. is_dynamic(i,itile)) act_fpc(i,itile) = 0._jp
          IF (.NOT. is_woody(i,itile)) non_woody_fpc(i) = non_woody_fpc(i) + act_fpc(i,itile)
       END DO
    END DO

    !-- adjust grass fractions to changes in woody type fpc

    !-- reset temporary arrays
    grass_sum(:)     = 0._jp

    !-- temporary sums for grass types
    DO itile = 1,ntiles
       DO i = 1,nidx
          IF (.NOT. is_woody(i,itile) .AND. is_dynamic(i,itile) .AND. bio_exist(i,itile) > 0.5_jp .AND. &
              npp_ave(i,itile) > EPSILON(1._jp)) THEN
            grass_sum(i) = grass_sum(i) + npp_ave(i,itile) * MAX(act_fpc_min,act_fpc(i,itile))
          END IF
       END DO
    END DO

!!$    !-- area wind break of woody types (trees and shrubs) and burned
!!$    IF (PRESENT(burned_frac)) THEN
!!$      CALL disturbed_frac(lctlib, nidx,kidx0,kidx1,ntiles,DIST_FIRE, &
!!$                          dynamic_pft(:) .OR. (dist_opts%lburn_pasture .AND. lctlib%pasture_pft(:)), with_yasso, &
!!$                          act_fpc,veg_fract_correction,surf,cbal,climbuf,burned_frac,burned_frac_diag,fuel)
!!$      CALL disturbed_frac(lctlib, nidx,kidx0,kidx1,ntiles,DIST_WINDBREAK,woody_pft(:) .AND. dynamic_pft(:), with_yasso, &
!!$                          act_fpc,veg_fract_correction,surf,cbal,climbuf,damaged_frac)
!!$    END IF

    !-- dynamic equation for act_fpc, woody types

    DO i = 1,nidx
      veg_inc = 0._jp
      DO itile = 1,ntiles
        tau_pft = tau_Cpool_woods(i,itile) * tau_Cpool_woods_to_tau_pft
        IF (is_dynamic(i,itile) .AND. is_woody(i,itile)) THEN
          veg_inc = veg_inc + MAX(fract_small,act_fpc(i,itile) + accelerate_nlcc * &
          ((pot_fpc(i,itile) - act_fpc(i,itile)) / tau_pft &
          - (burned_fract(i,itile) + damaged_fract(i,itile) + min_dist_woody) * act_fpc_previous_day(i,itile))) &
          - act_fpc(i,itile)
        ELSEIF(is_dynamic(i,itile) .AND. .NOT. is_woody(i,itile)) THEN
          IF (bio_exist(i,itile) > 0.5_jp .AND. npp_ave(i,itile) > EPSILON(1._jp)) THEN
            veg_inc = veg_inc + MAX(fract_small,act_fpc(i,itile) - accelerate_nlcc * &
            (burned_fract(i,itile) * act_fpc_previous_day(i,itile) + &
            grass_mortality * act_fpc(i,itile) - &
            bare_fpc(i) / grass_sum(i) * npp_ave(i,itile) * MAX(act_fpc_min,act_fpc(i,itile)) / tau_pft)) &
            - act_fpc(i,itile)
          ELSE
            veg_inc = veg_inc + MAX(fract_small,act_fpc(i,itile) - accelerate_nlcc * &
            (burned_fract(i,itile) * act_fpc_previous_day(i,itile) + &
            act_fpc(i,itile) / tau_pft)) &
            - act_fpc(i,itile)
          ENDIF
        ENDIF
      END DO
      theta = 1._jp
      IF (bare_fpc(i) < (veg_inc + fract_small)) theta = 0._jp
      DO itile = 1,ntiles
        tau_pft = tau_Cpool_woods(i,itile) * tau_Cpool_woods_to_tau_pft
        IF (is_dynamic(i,itile) .AND. is_woody(i,itile)) THEN
          act_fpc(i,itile) = MAX(fract_small,act_fpc(i,itile) + accelerate_nlcc * &
            ((theta * pot_fpc(i,itile) - act_fpc(i,itile)) / tau_pft                &
            - (burned_fract(i,itile) + damaged_fract(i,itile) + min_dist_woody) * act_fpc_previous_day(i,itile)))
        ENDIF
      ENDDO
    ENDDO

    !-- dynamic equation for act_fpc, grass types

    DO i = 1,nidx
      DO itile = 1,ntiles
        tau_pft = tau_Cpool_woods(i,itile) * tau_Cpool_woods_to_tau_pft
        IF (is_dynamic(i,itile) .AND. .NOT. is_woody(i,itile)) THEN
          IF (bio_exist(i,itile) > 0.5_jp .AND. npp_ave(i,itile) > EPSILON(1._jp)) THEN
            act_fpc(i,itile) = MAX(fract_small,act_fpc(i,itile) - accelerate_nlcc * &
               (burned_fract(i,itile) * act_fpc_previous_day(i,itile) + &
               grass_mortality * act_fpc(i,itile) - &
               bare_fpc(i) / grass_sum(i) * npp_ave(i,itile) * MAX(act_fpc_min,act_fpc(i,itile)) / tau_pft))
          ELSE
            act_fpc(i,itile) = MAX(fract_small,act_fpc(i,itile) - accelerate_nlcc * &
               (burned_fract(i,itile) * act_fpc_previous_day(i,itile) + &
               act_fpc(i,itile) / tau_pft))
          ENDIF
        ENDIF
      ENDDO
    ENDDO

    !-- calculate total area fraction covered by dynamic vegetation

    total_act_fpc(:) = 0._jp

    DO itile = 1,ntiles
       DO i = 1,nidx
          IF (is_dynamic(i,itile)) total_act_fpc(i) = total_act_fpc(i) + act_fpc(i,itile)
       END DO
    END DO

    !-- calculation of bare soil fraction

    bare_fpc(:) = MAX(0._jp,1._jp - total_act_fpc(:))

  END SUBROUTINE fpc_daily

! --- scale_fpc ---------------------------------------------------------------------------

  SUBROUTINE scale_fpc (routine, nidx, iblk, ntiles, is_dynamic, act_fpc, bare_fpc)
!
! Rescaling of fractional plant coverage to assure that
!  - the sum of act_fpc and bare soil is exactly one
!  - all tiles have a minimum vegetated fraction
!
! input
    CHARACTER, INTENT(in) :: routine         ! name of the calling routine
    INTEGER,  INTENT(in)  :: nidx            ! vector length
    INTEGER,  INTENT(in)  :: iblk            ! number of block
    INTEGER,  INTENT(in)  :: ntiles          ! number of tiles (different land cover type mosaics per gridcell)
    LOGICAL,  INTENT(in)  :: is_dynamic(:,:) ! PFTs taking part in competition

! inout
    REAL(jp), INTENT(inout) :: act_fpc(:,:) ! actual fpc in nlcc
    REAL(jp), INTENT(inout) :: bare_fpc(:)  ! fraction of bare ground

! local variables
    LOGICAL   :: just_return
    INTEGER   :: i, itile       ! indices
    INTEGER   :: nsparce(nidx)  ! number of PFTs with less then fract_small vegetation
    REAL(jp)  :: sum_fpc(nidx)  ! sum of the different cover fractions
    REAL(jp)  :: excess         ! extra fraction, that needs to be redistibuted
    REAL(jp)  :: scalable       ! sum of fraction of tiles that can be scaled
    REAL(jp)  :: rescale        ! factor to rescale fraction of tiles

    ! test, if scaling is necessary
    just_return = .true.
    sum_fpc(:) = bare_fpc(:)
    DO i = 1,nidx
       DO itile = 1,ntiles
          sum_fpc(i) = sum_fpc(i) + act_fpc(i,itile)
          IF (is_dynamic(i,itile) .AND. &
              act_fpc(i,itile) < fract_small) just_return = .false.
       END DO
       IF (ABS(1._jp - sum_fpc(i)) > 3._jp * EPSILON(1._jp)) just_return = .false.
    END DO

    IF (just_return) RETURN

    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine),'Attention: act_fpc is rescaled!')

    ! sum act_fpc and bare_fpc
    sum_fpc(:) = bare_fpc(:)
    nsparce(:) = 0

    DO i = 1,nidx
       DO itile = 1,ntiles
          IF (is_dynamic(i,itile)) THEN
             IF (act_fpc(i,itile) > fract_small) THEN
                sum_fpc(i) = sum_fpc(i) + act_fpc(i,itile)
             ELSE
                nsparce(i) = nsparce(i) + 1
                act_fpc(i,itile) = fract_small
             END IF
          ELSE
             act_fpc(i,itile) = 0._jp
          END IF
       END DO
    END DO

    sum_fpc(:) = sum_fpc(:) + REAL(nsparce(:),jp) * fract_small

    ! scaling of bare_fpc
    WHERE (sum_fpc(:) <= 1._jp - EPSILON(1._jp))
       bare_fpc(:) = bare_fpc(:) + (1._jp - sum_fpc(:))
    ELSEWHERE (sum_fpc(:) >= 1._jp + EPSILON(1._jp))
       bare_fpc(:) = bare_fpc(:) / sum_fpc(:)
    END WHERE

    ! scaling of act_fpc
    DO i = 1,nidx
       excess = 0._jp
       DO itile = 1,ntiles
          IF (is_dynamic(i,itile)) THEN
             IF (sum_fpc(i) >= 1._jp + EPSILON(1._jp) .AND. act_fpc(i,itile) > fract_small * sum_fpc(i)) THEN
                act_fpc(i,itile) = act_fpc(i,itile) / sum_fpc(i)
             ELSE IF (sum_fpc(i) >= 1._jp + EPSILON(1._jp)) THEN
                excess = excess + fract_small - act_fpc(i,itile) / sum_fpc(i)
                act_fpc(i,itile) = fract_small
             END IF
          END IF
       END DO
       scalable = bare_fpc(i)
       DO itile = 1,ntiles
          IF (is_dynamic(i,itile)) THEN
             IF (act_fpc(i,itile) > 2._jp * fract_small) scalable = scalable + act_fpc(i,itile)
          END IF
       END DO
       rescale = (scalable - excess) / scalable
       bare_fpc(i) = bare_fpc(i) * rescale
       DO itile = 1,ntiles
          IF (is_dynamic(i,itile)) THEN
             IF (act_fpc(i,itile) > 2._jp * fract_small) act_fpc(i,itile) = act_fpc(i,itile) * rescale
          END IF
       END DO
    END DO

  END SUBROUTINE scale_fpc


! --- calc_cover_fract_inst ---------------------------------------------------------------------------

  SUBROUTINE calc_cover_fract_inst(nidx, ntiles, is_woody, is_dynamic, is_pasture, is_crop, &
                                   cover_fract_pot, cover_fract, cover_fract_inst)
!
! convert fractional plant cover calculated within nlcc to jsbach cover fractions
! that would result from an instantaneous establishment of crop land and pasture
! (pasture is preferentially established on grasslands)
!
! input
    INTEGER,  INTENT(in)  :: nidx              ! vector length
    INTEGER,  INTENT(in)  :: ntiles            ! number of tiles (different land cover type mosaics per gridcell)
    LOGICAL,  INTENT(in)  :: is_woody(:,:)     ! PFTs that build wood
    LOGICAL,  INTENT(in)  :: is_dynamic(:,:)   ! PFTs taking part in competition
    LOGICAL,  INTENT(in)  :: is_pasture(:,:)   ! PFTs that are pasture
    LOGICAL,  INTENT(in)  :: is_crop(:,:)      ! PFTs that are crop
    REAL(jp), INTENT(in)  :: cover_fract_pot(:,:)  ! potential cover fractions (without land use)
    REAL(jp), INTENT(in)  :: cover_fract(:,:)  ! cover fraction within jsbach

! output
    REAL(jp), INTENT(out) :: cover_fract_inst(:,:) ! cover fraction within jsbach resulting from an
                                                   ! instantaneous establishment of crop land and pasture
! local variables
    INTEGER   :: i, itile
    INTEGER   :: ngrass                  ! number of grass types
    INTEGER   :: ncrops                  ! number of crop types
    REAL(jp)  :: sum_grass_fract(nidx)   ! potential fraction of grass if there was no land use
    REAL(jp)  :: sum_pasture_fract(nidx) ! fraction of the vegetated part of the grid box covered with pasture
    REAL(jp)  :: sum_crop_fract(nidx)    ! fraction of the vegetated part of the grid box covered with crops

    !-- sum up FPC of grasses, pastures and crops
    sum_grass_fract(:) = 0._jp
    sum_pasture_fract(:) = 0._jp
    sum_crop_fract(:) = 0._jp
    ngrass = 0
    ncrops = 0
    DO itile = 1, ntiles
       IF (.NOT. is_woody(1,itile) .AND. is_dynamic(1,itile)) THEN
          ngrass = ngrass + 1
          sum_grass_fract(:) = sum_grass_fract(:) + cover_fract_pot(:,itile)
       ELSE IF (is_pasture(1,itile)) THEN
          sum_pasture_fract(:) = sum_pasture_fract(:) + cover_fract(:,itile)
       ELSE IF (is_crop(1,itile)) THEN
          ncrops = ncrops + 1
          sum_crop_fract(:) = sum_crop_fract(:) + cover_fract(:,itile)
       END IF
    END DO

    !-- establishment of pastures (at the expense of grass land)
    DO i = 1,nidx
       DO itile = 1,ntiles
          IF (is_dynamic(i,itile)) THEN
             IF (sum_pasture_fract(i) <= sum_grass_fract(i)) THEN
                ! There are enough grasslands
                IF (is_woody(i,itile)) THEN
                   cover_fract_inst(i,itile) = cover_fract_pot(i,itile)
                ELSE
                   cover_fract_inst(i,itile) = cover_fract_pot(i,itile) &
                     & * (sum_grass_fract(i) - sum_pasture_fract(i) + ngrass * fract_small) / sum_grass_fract(i)
                END IF
             ELSE
                ! There are not enough grasslands for all pastures, so that woody types are also reduced
                IF (is_woody(i,itile)) THEN
                   cover_fract_inst(i,itile) = cover_fract_pot(i,itile) &
                     &  * (1._jp - sum_pasture_fract(i)) / (1._jp - sum_grass_fract(i))
                ELSE
                   cover_fract_inst(i,itile) = fract_small
                END IF
             END IF
          ELSE IF (is_pasture(i,itile)) THEN
             cover_fract_inst(i,itile) = cover_fract(i,itile)
          END IF
       END DO
    END DO

    !-- establishment of crops (at the expense of all natural pfts)
    DO itile = 1,ntiles
       WHERE (is_crop(:,itile))
          cover_fract_inst(:,itile) = cover_fract(:,itile)
       ELSEWHERE (.NOT. is_pasture(:,itile))
          cover_fract_inst(:,itile) = cover_fract_inst(:,itile) &
            &  * (1._jp - sum_pasture_fract(:) - sum_crop_fract(:)) &
            &  / (1._jp - sum_pasture_fract(:) - ncrops * fract_small)
       END WHERE
    END DO

    CALL scale_cover_fract (nidx, ntiles, is_dynamic(:,:), cover_fract_inst(:,:))

    IF (ANY(SUM(cover_fract_inst(:,:),DIM=2) > 1._jp + REAL(ntiles,jp)*EPSILON(1._jp)) .OR. &
      & ANY(SUM(cover_fract_inst(:,:),DIM=2) < 1._jp - REAL(ntiles,jp)*EPSILON(1._jp))) THEN
       WRITE (message_text,*) 'sum of cover_fract_inst /= 1: ', &
         &   MINVAL(SUM(cover_fract_inst(:,:),DIM=2)), MAXVAL(SUM(cover_fract_inst(:,:),DIM=2)), ', MINLOC: ', &
         &   MINLOC(SUM(cover_fract_inst(:,:),DIM=2)), ', MAXLOC: ', MAXLOC(SUM(cover_fract_inst(:,:),DIM=2))
       CALL finish ('calc_cover_fract_inst', message_text)
    END IF
    DO itile = 1,ntiles
       IF (ANY(cover_fract_inst(:,itile) < fract_small)) THEN
          WRITE (message_text,*) 'cover_fract_inst too small: ', MINVAL(cover_fract_inst(:,itile)), &
            &  MINLOC(cover_fract_inst(:,itile))
          CALL finish ('calc_cover_fract_inst', message_text)
       END IF
    END DO

  END SUBROUTINE calc_cover_fract_inst


!------------------------------------------------------------------------------

SUBROUTINE scale_cover_fract (nidx, ntiles, is_dynamic, cover_fract)
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
    INTEGER,  INTENT(in)  :: nidx               ! vector length
    INTEGER,  INTENT(in)  :: ntiles             ! number of tiles
    LOGICAL,  INTENT(in)  :: is_dynamic(:,:)    ! mask for land surface treated by jsbach
!
! !IN- and OUTPUT PARAMETERS:
!
    REAL(jp), INTENT(inout) :: cover_fract(:,:) ! vegetated fraction

! !LOCAL VARIABLES:
!
    INTEGER   :: itile, ic, iter   ! looping indices
    INTEGER   :: niter             ! number of iterations needed
    INTEGER   :: nsparce(nidx)     ! number of PFTs with a vegetated fraction of less then fract_small
    REAL(jp)  :: sum_fract(nidx)   ! sum of all cover fractions
    REAL(jp)  :: excluded_fract(nidx) ! sum of the PFT fractions that cannot be changed by natural vegetation
                                      ! dynamics: non-dynamic fractions and minimum fractions of dynamic PFTs
    REAL(jp)  :: anthro_fract(nidx)   ! sum of non-dynamic PFT fractions
    REAL(jp)  :: considerable_nat_fract(nidx) ! sum of natural PFTs with considerable fractions
    LOGICAL   :: corrected(nidx)      ! flag indicating cells with corrected fractions

    CHARACTER(len=*), PARAMETER :: routine = modname//':scale_cover_fract'
!------------------------------------------------------------------------------

    ! Make sure, crop and pasture have a cover fraction of at least fract_small
    WHERE (.NOT. is_dynamic(:,:))
       cover_fract(:,:) = MAX(fract_small,cover_fract(:,:))
    END WHERE

    ! Make sure, the sum of non-dynamic PFTs is smaller than (1 - fract_small*n_natural)
    excluded_fract(:) = 0._jp
    DO itile = 1,ntiles
       WHERE (is_dynamic(:,itile))
          excluded_fract(:) = excluded_fract(:) + fract_small
       ELSEWHERE
          excluded_fract(:) = excluded_fract(:) + cover_fract(:,itile)
       END WHERE
    END DO
    IF (ANY(excluded_fract(:) > 1._jp)) THEN
       corrected(:) = .FALSE.
       DO itile = 1, ntiles
          WHERE (excluded_fract(:) > 1._jp .AND. .NOT. corrected(:) .AND. .NOT. is_dynamic(:,itile) &
               &  .AND. cover_fract(:,itile) > (fract_small + excluded_fract(:)-1._jp))
             cover_fract(:,itile) = cover_fract(:,itile) - (excluded_fract(:)-1._jp)
             corrected(:) = .TRUE.
          END WHERE
       END DO
    END IF

    ! The more tiles we have, the more iterations are needed. For binary identical results
    ! the number of iterations must not change.
    niter = ntiles
    DO iter = 1, niter

       sum_fract(:) = 0._jp
       considerable_nat_fract(:) = 0._jp
       anthro_fract(:) = 0._jp
       nsparce(:) = 0

       DO itile = 1,ntiles
          sum_fract(:) = sum_fract(:) + cover_fract(:,itile)
          WHERE (cover_fract(:,itile) > fract_small .AND. is_dynamic(:,itile))
             considerable_nat_fract(:) = considerable_nat_fract(:) + cover_fract(:,itile)
          ELSEWHERE (is_dynamic(:,itile))
             nsparce(:) = nsparce(:) + 1
          ELSEWHERE
             anthro_fract(:) = anthro_fract(:) + cover_fract(:,itile)
          END WHERE
       END DO
       DO itile = 1,ntiles
          WHERE (cover_fract(:,itile) > fract_small .AND. is_dynamic(:,itile))
             cover_fract(:,itile) = cover_fract(:,itile) * &
               &                    (1._jp - anthro_fract(:) - REAL(nsparce,jp)*fract_small) / considerable_nat_fract(:)
          ELSEWHERE (is_dynamic(:,itile))
             cover_fract(:,itile) = fract_small
          ELSEWHERE
             cover_fract(:,itile) = MAX(fract_small,cover_fract(:,itile))
          END WHERE
       END DO
    END DO  ! iter

    ! Final corrections, to assure the sum of all PFT fractions equals 1. This sometimes might not be
    ! the case for numerical reasons.

    sum_fract(:) = SUM(cover_fract(:,:), DIM=2)
    ! Exit if the cover fraction sum is considerably different from 1.
    IF (ANY(ABS(1._jp - sum_fract(:)) > fract_small)) THEN
       WRITE(message_text,*) 'Sum of PFT fractions differs from 1.: ', &
         &  MINVAL(sum_fract(:)), MAXVAL(sum_fract(:)), ' MINLOC: ', &
         &  MINLOC(sum_fract(:)), ' MAXLOC: ', MAXLOC(sum_fract(:))
       CALL finish (routine, message_text)
    END IF

    corrected(:) = .FALSE.
    DO itile = 1,ntiles
       WHERE (is_dynamic(:,itile) .AND. (sum_fract(:) < 1._jp - EPSILON(1._jp)) &
            &  .AND. .NOT. corrected(:) .AND. cover_fract(:,itile) < sum_fract(:))
          cover_fract(:,itile) = cover_fract(:,itile) + (1._jp - sum_fract(:))
          corrected(:) = .TRUE.
       END WHERE
    END DO
    corrected(:) = .FALSE.
    DO itile = 1,ntiles
       WHERE (is_dynamic(:,itile) .AND. (sum_fract(:) > 1._jp + EPSILON(1._jp)) &
            &  .AND. .NOT. corrected(:) .AND. cover_fract(:,itile) > 1._jp - sum_fract(:) - fract_small)
          cover_fract(:,itile) = cover_fract(:,itile) - (sum_fract(:) - 1._jp)
          corrected(:) = .TRUE.
       END WHERE
    END DO

    IF (ANY(cover_fract(:,:) < fract_small)) THEN
       WRITE(message_text,*) 'cover_fract still smaller ', fract_small, ' after ', niter, ' iterations:', &
          &  MINVAL(MERGE(cover_fract(:,:), 1._jp, cover_fract > 0._jp)), 'MINLOC: ', &
          &  MINLOC(MERGE(cover_fract(:,:), 1._jp, cover_fract > 0._jp))
       CALL message (routine, message_text)

       ! Sometimes cover_fract remains slightly smaller than fract_small for numerical reasons
       WHERE (cover_fract(:,:) < fract_small .AND. cover_fract(:,:) >= fract_small - REAL(ntiles,jp)*EPSILON(1._jp))
          cover_fract(:,:) = fract_small
       END WHERE
       IF (ANY(cover_fract(:,:) < fract_small)) THEN
          WRITE(message_text,*) 'cover_fract still smaller ', fract_small, ' after ', niter, ' iterations:', &
             &  MINVAL(MERGE(cover_fract(:,:), 1._jp, cover_fract > 0._jp)),'MINLOC: ', &
             &  MINLOC(MERGE(cover_fract(:,:), 1._jp, cover_fract > 0._jp))
          CALL finish (routine, message_text)
       END IF
    END IF

  END SUBROUTINE scale_cover_fract

#endif
END MODULE mo_nlcc_process
