MODULE mo_jsb_subset

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: jsbach_subsets
  PUBLIC :: t_subset, ON_DOMAIN, ON_CHUNK

  ENUM, BIND(C)
    ENUMERATOR :: ON_DOMAIN=1, ON_CHUNK
  END ENUM

  TYPE t_subset
    INTEGER :: type
    INTEGER :: ics
    INTEGER :: ice
    INTEGER :: iblk
    INTEGER :: nb
    INTEGER :: nc
  END TYPE t_subset

  TYPE t_model_subsets
    TYPE(t_subset), ALLOCATABLE :: sub (:)
  END TYPE t_model_subsets

  TYPE(t_model_subsets), ALLOCATABLE, TARGET, SAVE :: jsbach_subsets(:)

END MODULE mo_jsb_subset