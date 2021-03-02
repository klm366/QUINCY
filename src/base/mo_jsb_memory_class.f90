!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/license .
!!
!! @brief Contains the class and methods for the JSBACH memory structure
!!
!! @author
!!  Reiner Schnur, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Reiner Schnur (2013-04-04)
!! Copied to QUINCY by Jan Engel (2016-08-15) & updated 2017-05-16
!!
MODULE mo_jsb_memory_class
#ifndef __NO_JSBACH__

  USE mo_kind,                ONLY: wp
  USE mo_exception,           ONLY: finish
  USE mo_jsb_config_class,    ONLY: t_jsb_config
  USE mo_jsb_var_class,       ONLY: t_jsb_var_p, t_jsb_var, t_jsb_var_real2d, t_jsb_var_real3d, REAL2D, REAL3D, &
                                    t_jsb_struct_p, t_jsb_struct, t_jsb_n_struct_p, t_jsb_n_struct   ! quincy
  USE mo_jsb_grid_class,      ONLY: t_jsb_grid, t_jsb_vgrid
  USE mo_jsb_io,              ONLY: t_cf, t_grib1, t_grib2

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: t_jsb_memory_p, t_jsb_memory, t_jsb_param

  TYPE, ABSTRACT :: t_jsb_param
  END TYPE t_jsb_param

  ! This type contains the memory state for each tile for one process
  TYPE, ABSTRACT :: t_jsb_memory
    INTEGER :: id = 0
    INTEGER :: grid_id = 0
    INTEGER :: owner_model_id = 0
    INTEGER :: owner_proc_id = 0
    CHARACTER(LEN=:), ALLOCATABLE :: owner_proc_name
    INTEGER :: no_of_children = 0
    TYPE(t_jsb_memory_p), POINTER :: children(:) => NULL()
    CLASS(t_jsb_memory) , POINTER :: parent => NULL()
    INTEGER, ALLOCATABLE :: owner_tile_path(:)
    CHARACTER(len=10)    :: owner_tile_name
    INTEGER :: max_no_of_vars = 0
    INTEGER :: no_of_vars = 0
    INTEGER :: max_no_of_structs = 0                  ! *slmjsb4
    INTEGER :: no_of_structs = 0                      ! *slmjsb4
    INTEGER :: max_no_of_n_structs = 0                ! *slmjsb4    --> "n" is for "nested"
    INTEGER :: no_of_n_structs = 0                    ! *slmjsb4    --> "n" is for "nested"
    CHARACTER(LEN=50) :: varlist_name
    CHARACTER(LEN=50) :: structlist_name              ! *slmjsb4
    CHARACTER(LEN=50) :: n_structlist_name            ! *slmjsb4    --> "n" is for "nested"
    TYPE(t_jsb_var_p), ALLOCATABLE :: vars(:)
    TYPE(t_jsb_struct_p), ALLOCATABLE :: structs(:)       ! *slmjsb4
    TYPE(t_jsb_n_struct_p), ALLOCATABLE :: n_structs(:)   ! *slmjsb4  --> "n" is for "nested"
    CLASS(t_jsb_param), POINTER :: param => NULL()
    LOGICAL :: in_var_groups = .FALSE.
  CONTAINS
    PROCEDURE (Init_memory), DEFERRED, PASS(mem) :: Init
    !> QUINCY
    !>
#ifdef __QUINCY_NO_JSB__
    PROCEDURE                                    :: t_jsb_memory_quincy_add_var_r2d
    PROCEDURE                                    :: t_jsb_memory_quincy_add_var_r3d
    GENERIC                                      :: Add_var => t_jsb_memory_quincy_add_var_r2d, &
                                                               t_jsb_memory_quincy_add_var_r3d
    PROCEDURE                                    :: t_jsb_memory_quincy_add_simple_struct   ! mem%simple_struct 
    PROCEDURE                                    :: t_jsb_memory_quincy_add_nested_struct   ! mem%nested-struct
    GENERIC                                      :: Add_pool => t_jsb_memory_quincy_add_simple_struct, &
                                                                t_jsb_memory_quincy_add_nested_struct
    !> JSBACH4
    !>
#else
    ! PROCEDURE                                    :: Get_var
    ! PROCEDURE                                    :: t_jsb_memory_add_var_r2d
    ! PROCEDURE                                    :: t_jsb_memory_add_var_r3d
    ! PROCEDURE                                    :: t_jsb_memory_add_pool
    ! PROCEDURE                                    :: t_jsb_memory_add_pool_var_r2d
    ! GENERIC                                      :: Add_var => t_jsb_memory_add_var_r2d, t_jsb_memory_add_var_r3d, &
    !   &                                                        t_jsb_memory_add_pool, t_jsb_memory_add_pool_var_r2d
    ! PROCEDURE                                    :: Get_var_position
    ! PROCEDURE                                    :: Collect_var_real2d  ! Collect var from its child tiles
    !                                                                     ! R: Note, we do not want to keep a "var collection" in the memory
    !                                                                     !    because we already have all child values in the memory.
    ! PROCEDURE                                    :: Collect_var_real3d  ! Collect var from its child tiles
    ! PROCEDURE                                    :: HandDown_var_real2d ! Hand back var to child tiles
    ! PROCEDURE                                    :: Average_over_children_real2d ! Average var over the child tiles
#endif
  END TYPE t_jsb_memory

  TYPE t_jsb_memory_p
    CLASS(t_jsb_memory), POINTER :: p => NULL()
  END TYPE t_jsb_memory_p

  ! Note: argument "mem" needs to have the TARGET attribute so that we can point to individual variables
  !       (of class t_jsb_var) contained in the memory structure.
  ABSTRACT INTERFACE
    SUBROUTINE Init_memory(mem, prefix, suffix, lct_ids, model_id)
      IMPORT :: t_jsb_memory, t_jsb_config
      CLASS(t_jsb_memory), INTENT(inout), TARGET :: mem
      CHARACTER(LEN=*),    INTENT(in)    :: prefix, suffix
      INTEGER,             INTENT(in)    :: lct_ids!(:)   ! *slmtemp
      INTEGER,             INTENT(in)    :: model_id
    END SUBROUTINE Init_memory
  END INTERFACE

!  INTERFACE t_jsb_memory
!    PROCEDURE Create_process_memory
!  END INTERFACE t_jsb_memory

  CHARACTER(len=*), PARAMETER :: modname = 'mo_jsb_memory_class'


CONTAINS


  ! ======================================================================================================= !
  !> add 2D variable to process memory
  !>
  !> mem <- var
  !>
  SUBROUTINE t_jsb_memory_quincy_add_var_r2d( &
    mem, name, var, &
    hgrid, vgrid, cf, grib1, grib2, prefix, suffix, &
    loutput, lrestart, output_level, &
    initval_r)

    CLASS(t_jsb_memory),              INTENT(inout)   :: mem            ! memory structure
    CHARACTER(len=*),                 INTENT(in)      :: name           ! name of variable
    TYPE(t_jsb_var_real2d),   TARGET, INTENT(inout)   :: var            ! reference to jsbach var -> newly created plain scalar REAL
    TYPE(t_jsb_grid),        POINTER, INTENT(in)      :: hgrid          ! horizontal grid type
    TYPE(t_jsb_vgrid),       POINTER, INTENT(in)      :: vgrid          ! vertical grid type
    TYPE(t_cf),                       INTENT(in)      :: cf             ! CF related metadata
    TYPE(t_grib1),                    INTENT(in)      :: grib1          ! GRIB1 related metadata
    TYPE(t_grib2),                    INTENT(in)      :: grib2          ! GRIB2 related metadata
    CHARACTER(len=*),                 INTENT(in)      :: prefix         ! Prefix for variable name
    CHARACTER(len=*),                 INTENT(in)      :: suffix         ! Suffix for variable name
    LOGICAL,                OPTIONAL, INTENT(in)      :: loutput        ! output flag
    LOGICAL,                OPTIONAL, INTENT(in)      :: lrestart       ! restart flag
    INTEGER,                OPTIONAL, INTENT(in)      :: output_level   ! minimum level of output groups this var belongs to
    REAL(wp),               OPTIONAL, INTENT(in)      :: initval_r      ! value if var not available

    INTEGER                     :: pos_of_var
    INTEGER                     :: ice, ics, nproma, iblk
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':t_jsb_memory_quincy_add_var_r2d'

    ics     = 1              ! is always 1, see land_model_driver
    ice     = hgrid%nproma   ! number of points in the chunk
    nproma  = hgrid%nproma   ! number of points in the chunk
    iblk    = 1              ! is always 1, see land_model_driver

    ALLOCATE(var%field(nproma,iblk,1,1,1))
    var%ptr => var%field(:,:,1,1,1)

    var%name          = TRIM(name)
    var%owner_proc_id = mem%owner_proc_id

    pos_of_var = mem%no_of_vars + 1
    IF (pos_of_var > mem%max_no_of_vars) THEN
      CALL finish(TRIM(routine), 'Maximum number of vars exceeded: '//TRIM(mem%varlist_name)//', '//TRIM(name))
    ELSE
      mem%vars(pos_of_var)%p => var
      mem%no_of_vars = pos_of_var
    END IF

    ! added for QUINCY compatibility (in JSBACH4 the allocation is done differently)
    var%ptr(ics:ice,iblk) = initval_r

  END SUBROUTINE t_jsb_memory_quincy_add_var_r2d


  ! ======================================================================================================= !
  !> add 3D variable to process memory
  !>
  !> mem <- var
  !> 
  SUBROUTINE t_jsb_memory_quincy_add_var_r3d( &
    mem, name, var, &
    hgrid, vgrid, cf, grib1, grib2, prefix, suffix, &
    loutput, lrestart, output_level, &
    initval_r)
  
    CLASS(t_jsb_memory),              INTENT(inout)   :: mem            ! memory structure
    CHARACTER(len=*),                 INTENT(in)      :: name           ! name of variable
    TYPE(t_jsb_var_real3d),   TARGET, INTENT(inout)   :: var            ! reference to jsbach var -> newly created plain scalar REAL
    TYPE(t_jsb_grid),        POINTER, INTENT(in)      :: hgrid          ! horizontal grid type
    TYPE(t_jsb_vgrid),       POINTER, INTENT(in)      :: vgrid          ! vertical grid type
    TYPE(t_cf),                       INTENT(in)      :: cf             ! CF related metadata
    TYPE(t_grib1),                    INTENT(in)      :: grib1          ! GRIB1 related metadata
    TYPE(t_grib2),                    INTENT(in)      :: grib2          ! GRIB2 related metadata
    CHARACTER(len=*),                 INTENT(in)      :: prefix         ! Prefix for variable name
    CHARACTER(len=*),                 INTENT(in)      :: suffix         ! Suffix for variable name
    LOGICAL,                OPTIONAL, INTENT(in)      :: loutput        ! output flag
    LOGICAL,                OPTIONAL, INTENT(in)      :: lrestart       ! restart flag
    INTEGER,                OPTIONAL, INTENT(in)      :: output_level   ! minimum level of output groups this var belongs to
    REAL(wp),               OPTIONAL, INTENT(in)      :: initval_r      ! value if var not available

    INTEGER                     :: pos_of_var
    INTEGER                     :: ice, ics, nproma, iblk
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':t_jsb_memory_quincy_add_var_r3d'
  
    ics     = 1              ! is always 1, see land_model_driver
    ice     = hgrid%nproma   ! number of points in the chunk
    nproma  = hgrid%nproma   ! number of points in the chunk
    iblk    = 1              ! is always 1, see land_model_driver
  
    ALLOCATE(var%field(nproma, vgrid%n_levels, iblk, 1, 1))
    var%ptr => var%field(:,:,:,1,1)
  
    var%name          = TRIM(name)
    var%owner_proc_id = mem%owner_proc_id
  
    pos_of_var = mem%no_of_vars + 1
    IF (pos_of_var > mem%max_no_of_vars) THEN
      CALL finish(TRIM(routine), 'Maximum number of vars exceeded: '//TRIM(mem%varlist_name)//', '//TRIM(name))
    ELSE
      mem%vars(pos_of_var)%p => var
      mem%no_of_vars = pos_of_var
    END IF
  
    ! added for QUINCY compatibility (in JSBACH4 the allocation is done differently)
    var%ptr(ics:ice, 1:vgrid%n_levels, iblk) = initval_r
  
  END SUBROUTINE t_jsb_memory_quincy_add_var_r3d


  ! ======================================================================================================= !
  !> add simple pool structure to process memory
  !>
  !> mem <- pool  (a pool that contains only var, but no other pools)
  !>    
  SUBROUTINE t_jsb_memory_quincy_add_simple_struct( &
    mem, name, struct, &
    hgrid, vgrid, cf, grib1, grib2, prefix, suffix, &
    loutput, lrestart, output_level)

    CLASS(t_jsb_memory),           INTENT(inout) :: mem                 ! memory structure
    CHARACTER(len=*),              INTENT(in)    :: name                ! name of variable
    CLASS(t_jsb_struct),   TARGET, INTENT(inout) :: struct              ! struct
    TYPE(t_jsb_grid),     POINTER, INTENT(in)    :: hgrid               ! horizontal grid type
    TYPE(t_jsb_vgrid),    POINTER, INTENT(in)    :: vgrid               ! vertical grid type
    TYPE(t_cf),                    INTENT(in)    :: cf                  ! CF related metadata
    TYPE(t_grib1),                 INTENT(in)    :: grib1               ! GRIB1 related metadata
    TYPE(t_grib2),                 INTENT(in)    :: grib2               ! GRIB2 related metadata
    CHARACTER(len=*),              INTENT(in)    :: prefix              ! Prefix for pool name
    CHARACTER(len=*),              INTENT(in)    :: suffix              ! Suffix for pool name
    LOGICAL,             OPTIONAL, INTENT(in)    :: loutput             ! output flag
    LOGICAL,             OPTIONAL, INTENT(in)    :: lrestart            ! restart flag
    INTEGER,             OPTIONAL, INTENT(in)    :: output_level        ! minimum level of output groups this var belongs to

    INTEGER                     :: pos_of_struct
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':t_jsb_memory_quincy_add_simple_struct'

    struct%name          = TRIM(name)
    struct%owner_proc_id = mem%owner_proc_id

    pos_of_struct = mem%no_of_structs + 1
    IF (pos_of_struct > mem%max_no_of_structs) THEN
      CALL finish(TRIM(routine), 'Maximum number of structs exceeded: '//TRIM(mem%structlist_name)//', '//TRIM(name))
    ELSE
      mem%structs(pos_of_struct)%p => struct
      mem%no_of_structs = pos_of_struct
    END IF

    ALLOCATE(struct%vars(struct%max_no_of_vars))  !!! max_no_of_vars is hardcoded in mo_jsb_var_class:t_jsb_struct

  END SUBROUTINE t_jsb_memory_quincy_add_simple_struct


  ! ======================================================================================================= !
  !> add nested pool structure to process memory
  !>
  !> mem <- pool  (a pool that contains other pools, but no variables)
  !> 
  SUBROUTINE t_jsb_memory_quincy_add_nested_struct( &
    mem, name, struct, &
    hgrid, vgrid, cf, grib1, grib2, prefix, suffix, &
    loutput, lrestart, output_level)

    CLASS(t_jsb_memory),           INTENT(inout) :: mem                 ! memory structure
    CHARACTER(len=*),              INTENT(in)    :: name                ! name of variable
    CLASS(t_jsb_n_struct), TARGET, INTENT(inout) :: struct              ! nested struct
    TYPE(t_jsb_grid),     POINTER, INTENT(in)    :: hgrid               ! horizontal grid type
    TYPE(t_jsb_vgrid),    POINTER, INTENT(in)    :: vgrid               ! vertical grid type
    TYPE(t_cf),                    INTENT(in)    :: cf                  ! CF related metadata
    TYPE(t_grib1),                 INTENT(in)    :: grib1               ! GRIB1 related metadata
    TYPE(t_grib2),                 INTENT(in)    :: grib2               ! GRIB2 related metadata
    CHARACTER(len=*),              INTENT(in)    :: prefix              ! Prefix for pool name
    CHARACTER(len=*),              INTENT(in)    :: suffix              ! Suffix for pool name
    LOGICAL,             OPTIONAL, INTENT(in)    :: loutput             ! output flag
    LOGICAL,             OPTIONAL, INTENT(in)    :: lrestart            ! restart flag
    INTEGER,             OPTIONAL, INTENT(in)    :: output_level        ! minimum level of output groups this var belongs to

    INTEGER                     :: pos_of_n_struct
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':t_jsb_memory_quincy_add_nested_struct'

    struct%name          = TRIM(name)
    struct%owner_proc_id = mem%owner_proc_id

    pos_of_n_struct = mem%no_of_n_structs + 1
    IF (pos_of_n_struct > mem%max_no_of_n_structs) THEN
      CALL finish(TRIM(routine), 'Maximum number of nested structs exceeded: '//TRIM(mem%structlist_name)//', '//TRIM(name))
    ELSE
      mem%n_structs(pos_of_n_struct)%p => struct
      mem%no_of_n_structs = pos_of_n_struct
    END IF

    ALLOCATE(struct%l2_structs(struct%max_no_of_l2_structs))  !!! max_no_of_l2_structs is hardcoded in mo_jsb_var_class:t_jsb_n_struct

  END SUBROUTINE t_jsb_memory_quincy_add_nested_struct

#endif
END MODULE mo_jsb_memory_class
