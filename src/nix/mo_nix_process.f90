!>
!> Contact: icon-model.org
!> Authors: AUTHORS.md
!> See LICENSES/ for license information
!> SPDX-License-Identifier: BSD-3-Clause
!> ---------------------------------------
!>
MODULE mo_nix_process
#ifndef __NO_JSBACH__

  USE mo_kind,      ONLY: wp
  USE mo_exception, ONLY: message, message_text


  IMPLICIT NONE
  PRIVATE

  PUBLIC ::                        &
    & calc_nix_layering

  CHARACTER(len=*), PARAMETER :: modname = 'mo_nix_process'

CONTAINS

  !--------------------------------------------------------------------
  ! Compute snow cover layering - NIX
  !--------------------------------------------------------------------  

  SUBROUTINE calc_nix_layering()

    ! ---------------------------
    ! Go

    CALL message (TRIM(modname), '  Hello NIX!')

!    WRITE(message_text, *) 'Hello NIX!'


  END SUBROUTINE calc_nix_layering        

#endif
END MODULE mo_nix_process

