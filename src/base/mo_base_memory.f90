!>
!! @brief land model memory structure
!!
!! @author
!!  Sönke Zaehle, MPI-BGC Jena 
!!
!! @par Revision History
!! First version by Sönke Zaehle (2015-12-04)
!!
MODULE mo_base_memory

  USE mo_kind,                    ONLY: wp
  USE mo_jsb_lct_class
  USE mo_jsb_memory_class,        ONLY: t_jsb_memory_p

  IMPLICIT NONE

  !-----------------------------------------------------------------------------------------------------
  !> integrates all variables of one tile
  !!
  !! TYPE(t_slm_tile) :: tile --> declaration of tile in the driver \n
  !! 
  !-----------------------------------------------------------------------------------------------------
  TYPE t_slm_tile
    TYPE(t_jsb_lct)                   :: lcts      !< data on land cover type of the tile
    CLASS(t_jsb_memory_p), POINTER    :: mem(:)
    INTEGER                           :: owner_model_id = 0   ! necessary for using lctlib
  END TYPE t_slm_tile  

  !-----------------------------------------------------------------------------------------------------
  !> type definition of model\%options needed in t_jsb_model
  !!
  !! this is borrowed from JSBACH4 git mo_jsb_task_class [20160818 Jan Engel] \n
  !! enabling the use of ics ice iblk in the 2D and 3D variables
  !!
  !! initialization in JSBACH4 is done in the interface with CALL model%Set_options (... ), which is a type-bound procedure 
  !! from t_jsb_model
  !!
  !! initilization in SLM is (temporary) done in mo_jsb_model_init
  !-----------------------------------------------------------------------------------------------------
  TYPE t_jsb_task_options
    REAL(wp) :: dtime   = -1._wp  !< length of one time step
    REAL(wp) :: steplen = -1._wp  !< .. 
    INTEGER  :: iblk    = -1      !< Number of current block (chunk)
    INTEGER  :: ics     = -1      !< Index of chunk start
    INTEGER  :: ice     = -1      !< Index of chunk end
    INTEGER  :: nc      = -1      !< Length of chunk; nc = ice - ics + 1
  END TYPE t_jsb_task_options

  !-----------------------------------------------------------------------------------------------------
  !> longitude and latitude of a tile
  !!
  !! note: this TYPE is just a temporary solution \n
  !! migration to /base/mo_grid_class may follow, depending on structure in JSBACH4
  !-----------------------------------------------------------------------------------------------------
  TYPE t_slm_grid 
    REAL(wp) :: lon
    REAL(wp) :: lat
  END TYPE t_slm_grid

  TYPE(t_slm_grid)  :: grid           !< currently initialized with namelist parameter values

  CHARACTER(len=*), PARAMETER, PRIVATE :: modname = 'mo_base_memory'

  PRIVATE
  PUBLIC  :: t_slm_tile, t_jsb_task_options
  PUBLIC  :: grid
 
!CONTAINS

END MODULE mo_base_memory
