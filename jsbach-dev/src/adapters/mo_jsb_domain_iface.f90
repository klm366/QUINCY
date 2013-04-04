!>
!! @brief Contains interface to spatial domains from ICON/ECHAM6
!!
!! @author
!!  Reiner Schnur, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Reiner Schnur (2012-11-02)
!!
MODULE mo_jsb_domain_iface

#ifdef __ICON__
  USE mo_model_domain,       ONLY: t_patch, p_patch
  USE mo_impl_constants_grf, ONLY: grf_bdywidth_c
  USE mo_impl_constants,     ONLY: min_rlcell_int
  USE mo_loopindices,        ONLY: get_indices_c
#else
  USE mo_geoloc,          ONLY: philat_2d, philon_2d
#endif

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: t_patch
  PUBLIC :: get_blk_start, get_blk_end, get_col_start, get_col_end

#ifndef __ICON__
  TYPE t_patch
     INTEGER :: id = -1
  END type t_patch
#endif

CONTAINS

#ifdef __ICON__
  INTEGER FUNCTION get_blk_start(patch)

    TYPE(t_patch), INTENT(in) :: patch

    get_blk_start = patch%cells%start_blk(grf_bdywidth_c + 1, 1)

  END FUNCTION get_blk_start

  INTEGER FUNCTION get_blk_end(patch)

    TYPE(t_patch), INTENT(in) :: patch

    get_blk_end = patch%cells%end_blk(min_rlcell_int, MAX(1, patch%n_childdom))

  END FUNCTION get_blk_end

  INTEGER FUNCTION get_col_start(jb, patch)

    INTEGER,       INTENT(in) :: jb
    TYPE(t_patch), INTENT(in) :: patch

    INTEGER :: i_startblk, i_endblk, jcs, jce, rl_start, rl_end

    i_startblk = get_blk_start(patch)
    i_endblk   = get_blk_end(patch)
    rl_start   = grf_bdywidth_c + 1
    rl_end     = min_rlcell_int

    CALL get_indices_c(patch, jb, i_startblk, i_endblk, jcs, jce, rl_start, rl_end)

    get_col_start = jcs

  END FUNCTION get_col_start

  INTEGER FUNCTION get_col_end(jb, patch)

    INTEGER,       INTENT(in) :: jb
    TYPE(t_patch), INTENT(in) :: patch

    INTEGER :: i_startblk, i_endblk, jcs, jce, rl_start, rl_end

    i_startblk = get_blk_start(patch)
    i_endblk   = get_blk_end(patch)
    rl_start   = grf_bdywidth_c + 1
    rl_end     = min_rlcell_int

    CALL get_indices_c(patch, jb, i_startblk, i_endblk, jcs, jce, rl_start, rl_end)

    get_col_end = jce

  END FUNCTION get_col_end

#endif

END MODULE mo_jsb_domain_iface
