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

  SUBROUTINE calc_nix_layering(  &
    & nc,            &
    & t_sn,          &           !< snow layer temperature - element      [K]     
    & dzm_sn,        &           !< snow layer thickness                  [m]
    & top_sn)                    !< index of top (upper) snow layer       [-]

    INTEGER, INTENT(in)  :: nc                      !< dimensions

    REAL(wp), INTENT(INOUT), DIMENSION(:,:) :: &
      & t_sn,           &
      & dzm_sn

    INTEGER, INTENT(INOUT), DIMENSION(:) :: &
      & top_sn

    ! Local
    INTEGER  :: ic                              ! iteration over grid points


    ! ---------------------------
    ! Go
    CALL message (TRIM(modname), ' Running NIX layering!')


!    DO ic = 1,nc
!
!      WRITE(*,*) 'NIX CALC LAYERING'
!      WRITE(*,"(A, 1I4)")     'top_sn:',     top_sn(ic)
!
!      WRITE (*,"(A, 10F8.3)") 't_sn:',       t_sn(ic,:)
!      WRITE (*,"(A, 10F8.3)") 'dzm_sn:',     dzm_sn(ic,:)
!
!    END DO


  END SUBROUTINE calc_nix_layering        

#endif
END MODULE mo_nix_process

