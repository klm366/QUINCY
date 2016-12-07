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
  USE mo_jsb_lctlib_class,        ONLY: t_lctlib_element
  USE mo_jsb_memory_class,        ONLY: t_jsb_memory_p, t_jsb_memory
  USE mo_jsb_process_class,       ONLY: t_jsb_process_p
  USE mo_shared_memory,           ONLY: t_element

  IMPLICIT NONE

  !-----------------------------------------------------------------------------------------------------
  !> soil pool
  !! 
  !-----------------------------------------------------------------------------------------------------
  TYPE t_slm_soil_pool
     TYPE(t_element)  :: litter         !< soil litter
     TYPE(t_element)  :: som            !< soil organic mater
  END TYPE t_slm_soil_pool

!  !-----------------------------------------------------------------------------------------------------
!  !> soil flux
!  !! 
!  !-----------------------------------------------------------------------------------------------------
!  TYPE t_slm_soil_flux
!    TYPE(t_slm_soil_pool) :: loss      !< decomposition flux from each pool
!  END TYPE t_slm_soil_flux

  !-----------------------------------------------------------------------------------------------------
  !> all (or most ?) variables of soil per tile
  !! 
  !! including all pools and all fluxes \n
  !! @todo this will be extended by an extra soil model
  !-----------------------------------------------------------------------------------------------------
  TYPE t_slm_soil
    TYPE(t_slm_soil_pool) :: pool     !< soil pool
    TYPE(t_slm_soil_pool) :: loss      !< decomposition flux from each pool
    !TYPE(t_slm_soil_flux) :: flux     !< soil flux
    REAL(wp) :: theta                 !< fractional water holding capacity
    REAL(wp) :: moisture              !< actual plant available water content (mm)
    REAL(wp) :: water_potential       !< actual soil water potential (MPa) 
    REAL(wp) :: ammonium              !< soil ammonium concentration (mmol/L)
    REAL(wp) :: phosphate             !< soil P concentration (mmol/L)
  END TYPE t_slm_soil

  !-----------------------------------------------------------------------------------------------------
  !> integrates all variables of one tile
  !!
  !! including soil memory \n
  !! TYPE(t_slm_tile) :: tile --> declaration of tile in the driver \n
  !! 
  !-----------------------------------------------------------------------------------------------------
  TYPE t_slm_tile
    TYPE(t_slm_soil)                  :: vegsoil      !< soil construct of tile  *slmdev

    TYPE(t_jsb_lct)                   :: lcts      !< data on land cover type of the tile
    CLASS(t_jsb_memory_p), POINTER    :: mem(:)
    INTEGER                           :: owner_model_id = 0   ! necessary for using lctlib
  END TYPE t_slm_tile  

  !-----------------------------------------------------------------------------------------------------
  !> type definition of model\%options needed in t_model
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
  !> very basic version of the JSBACH4: TYPE, EXTENDS(t_Hsm) :: t_jsb_model in mo_jsb_model_class
  !!
  !! 
  !! 
  !! 
  !-----------------------------------------------------------------------------------------------------
  TYPE t_model
    INTEGER                            :: id                !< Model id
    CHARACTER(len=10)                  :: shortname         !< Model short name
    !CHARACTER(len=132)                 :: description       !< Model description
    !CHARACTER(len=filename_max)        :: namelist_filename !< Model namelist file
    !TYPE(t_jsb_model_config),  POINTER :: config            !< Configuration
    TYPE(t_jsb_task_options)           :: options           !< Options needed by tasks
    !INTEGER                            :: grid_id
    TYPE(t_jsb_process_p),POINTER      :: processes(:)
    !TYPE(t_Task_queue),        POINTER :: queue        => NULL()
    !CLASS(t_jsb_process_task), POINTER :: current_task => NULL()
    TYPE(t_lctlib_element),POINTER     :: lctlib(:)         !< each lctlib(pft_id) contains all parameter values of pft_id \n
                                                            !! this structure makes parameters accessible throughout the model \n
                                                            !! consistently with jsbach4

  END TYPE t_model

  !-----------------------------------------------------------------------------------------------------
  !> Container CLASS for JSBACH model
  !! 
  !! Contains pointer to an instance of type t_jsb_model in order to allow arrays
  !! of pointers to models.
  !!
  !! This code is borrowed from jsbach4 mo_jsb_model_class
  !-----------------------------------------------------------------------------------------------------
!  TYPE t_jsb_model_m       ! jsbach4 code
!    TYPE(t_jsb_model), POINTER :: m
!  END TYPE t_jsb_model_m
  TYPE t_model_m            ! *slmtemp
    TYPE(t_model), POINTER :: m
  END TYPE t_model_m


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

 
CONTAINS

  !-----------------------------------------------------------------------------------------------------
  !> Constructor method for model instance
  !! 
  !!
  !! borrowed from jsbach4 mo_jsb_model_class and modified to fit SLM needs
  !!
  !-----------------------------------------------------------------------------------------------------
  FUNCTION new_model(id, name, shortname, description, namelist_filename) RESULT(return_ptr)

    IMPLICIT NONE

    ! 0.1 InOut
    INTEGER,                    INTENT(in) :: id
    CHARACTER(len=*),           INTENT(in) :: name
    CHARACTER(len=*),           INTENT(in) :: shortname
    CHARACTER(len=*),           INTENT(in) :: description
    CHARACTER(len=*),           INTENT(in) :: namelist_filename !< Model namelist file
    TYPE(t_model), POINTER                 :: return_ptr
    ! 0.2 Local
    !INTEGER :: iproc

    ALLOCATE(return_ptr)

!    CALL return_ptr%Init_hsm(name)
!    CALL return_ptr%Register_action(t_Action("INTEGRATE"))
!    CALL return_ptr%Register_action(t_Action("AGGREGATE"))

    return_ptr%id                = id
    return_ptr%shortname         = TRIM(shortname)
!    return_ptr%description       = TRIM(description)
!    return_ptr%namelist_filename = TRIM(namelist_filename)
!    return_ptr%config => NULL()

!!$    return_ptr%Hsm => return_ptr%t_Hsm

!!$    ALLOCATE(return_ptr%processes(max_no_of_processes))
!!$    DO iproc=1,max_no_of_processes
!!$      return_ptr%processes(iproc)%p => Create_process(iproc, return_ptr%id, TRIM(namelist_filename))
!!$    END DO

!    ALLOCATE(return_ptr%queue)

!    CALL message(TRIM(routine), 'construction of JSBACH model completed.')

  END FUNCTION new_model



END MODULE mo_base_memory
