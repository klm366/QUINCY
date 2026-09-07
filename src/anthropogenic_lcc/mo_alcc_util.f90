!> Contains utilities for anthropogenic land cover change
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
MODULE mo_alcc_util
#ifndef __NO_JSBACH__

  USE mo_kind,               ONLY: jp
  USE mo_exception,          ONLY: finish, message, message_text
  USE mo_jsb_control,        ONLY: debug_on

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: scale_anthro_fractions_for_nlcc

  CHARACTER(len=*), PARAMETER :: modname = 'mo_alcc_util'

CONTAINS
  !================================================================================================
  !>
  !>#### Scaling of anthropogenic tile fractions
  !>
  !> In case anthropogenic and natural land cover change processes are active, additional scaling
  !> of the anthropogenic tile fractions might be necessary to assure all natural PFTs can have an
  !> absolute fraction of at least fract_small. This is assured for the relative tile fractions
  !> defined in the cover fraction boundary condition files, but it is also mandatory for absolute
  !> fractions.
  !>
  !------------------------------------------------------------------------------------------------
  SUBROUTINE scale_anthro_fractions_for_nlcc(nc, npft, nlcc_tile, cf_veg, cf_pft)

    USE mo_nlcc_process,      ONLY: fract_small

    INTEGER,  INTENT(in)    :: &
      & npft,                  & !< Number of PFT tiles
      & nc                       !< Vector length
    LOGICAL,  INTENT(in)    :: &
      & nlcc_tile(:)             !< The tile's fraction is calculated in NLCC
    REAL(jp), INTENT(in)    :: &
      & cf_veg(:)                !< Cover fraction of the veg tile
    REAL(jp), INTENT(inout) :: &
      & cf_pft(:,:)              !< Cover fractions of the PFT tiles
    !----------------------------------------------------------------------------------------------
    INTEGER ::                 &
      & ic,                    & !< Looping index for grid cells
      & itile,                 & !< Looping index for tiles
      & ndyn                     !< Number of PFTs participating in nlcc
    REAL(jp) ::                &
      & sum_anthro(nc),        & !< Total fraction of anthropogenic PFTs
      & diff,                  & !< Mismatch of anthropogenic PFT fractions
      & to_scale,              & !< Sum of anthropogenic fractions participating in scaling
      & test_sum                 !< Sum of PFT fractions for testing

    CHARACTER(len=*), PARAMETER :: routine = modname//':scale_anthro_fractions'
    !----------------------------------------------------------------------------------------------

    IF (debug_on()) CALL message(routine, 'Starting routine')


    ! Find out the total fraction of anthropogenic PFTs
    sum_anthro(:) = 0._jp
    DO itile = 1, npft
      DO ic = 1, nc
        IF (.NOT. nlcc_tile(itile) .AND. cf_veg(ic) > 0._jp) THEN
          sum_anthro(ic) = sum_anthro(ic) + MAX(cf_pft(ic,itile), fract_small)
        END IF
      END DO
    END DO

    ! Maximum allowed area fraction for anthropogenic PFTs is 1 - ndyn*fract_small
    ! Scaling is needed if the antropogenic PFT fractions exceed this area.
    ndyn = COUNT(nlcc_tile(:))
    DO ic = 1, nc
      IF (sum_anthro(ic) > 1._jp - REAL(ndyn,jp)*fract_small + EPSILON(1._jp)) THEN

        diff = sum_anthro(ic) - (1._jp - REAL(ndyn,jp)*fract_small)

        ! Scale anthropogenic PFTs with relevant fractions:
        ! We only consider for scaling the area fractions exceeding fract_small by
        ! first subtracting this fraction, an re-adding it after the scaling.
        to_scale = 0._jp
        DO itile = 1, npft
          IF (.NOT. nlcc_tile(itile) .AND. cf_pft(ic,itile) > fract_small) THEN
            to_scale = to_scale + cf_pft(ic,itile) - fract_small
          END IF
        END DO
        test_sum = 0._jp
        DO itile = 1, npft
          IF (.NOT. nlcc_tile(itile) .AND. cf_pft(ic,itile) > fract_small) THEN
            cf_pft(ic,itile) = (cf_pft(ic,itile) - fract_small) * (to_scale - diff) / to_scale + fract_small
            test_sum = test_sum + cf_pft(ic,itile)
          ELSE IF (.NOT. nlcc_tile(itile)) THEN
            test_sum = test_sum + fract_small
          END IF
        END DO

        ! Final test
        IF (test_sum > 1._jp - REAL(ndyn,jp)*fract_small + EPSILON(1._jp)) THEN
          WRITE (message_text,*) &
            & ' Sum of anthropogenic PFT fractions still greater than (1. - ndyn*fract_small): ', &
            & ' ic: ', ic, ' PFT tile fractions: ', cf_pft(ic,:), ' Sum anthro: ', test_sum
          CALL finish(routine, message_text)
        END IF

      END IF
    END DO

    IF (debug_on()) CALL message(routine, 'Finished routine')

  END SUBROUTINE scale_anthro_fractions_for_nlcc

#endif
END MODULE mo_alcc_util
