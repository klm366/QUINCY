!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/license .
!!
!! @brief Contains the definition of variable types used across JSBACH4
!!
!! @par Revision History
!!  First version by Jan Engel (20160815); no changes compared to JSBACH4 except for not USEing mo_jsb_varlist_iface
!!  2nd   version by Jan Engel (20170927); remove operator overload functions from variables
!!  3rd   version by Jan Engel; add t_jsb_struct and t_jsb_n_struct for pool/fluxes infrastructure of QUINCY
!!
MODULE mo_jsb_var_class
#ifndef __NO_JSBACH__

  USE mo_kind,                ONLY: wp, dp
  USE mo_exception,           ONLY: finish, message
  USE mo_jsb_varlist_iface,   ONLY: VARNAME_LEN 
  USE mo_jsb_grid_class,      ONLY: t_jsb_grid, t_jsb_vgrid
  USE mo_jsb_io,              ONLY: t_cf, t_grib1, t_grib2

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: t_jsb_var_real2d, t_jsb_var_real3d, t_jsb_var, t_jsb_var_p
  PUBLIC :: t_jsb_struct, t_jsb_struct_p, t_jsb_n_struct, t_jsb_n_struct_p    ! quincy
  PUBLIC :: REAL2D, REAL3D


  ENUM, BIND(C)
    ENUMERATOR :: REAL2D=1, REAL3D
  END ENUM

  !-----------------------------------------------------------------------------------------------------
  !> Abstract Variable TYPE - t_jsb_var
  !! 
  !! 
  !-----------------------------------------------------------------------------------------------------
  TYPE, ABSTRACT :: t_jsb_var
    REAL(wp), ALLOCATABLE      :: field(:,:,:,:,:)
    CHARACTER(len=VARNAME_LEN) :: name            ! Name of variable
    CHARACTER(len=VARNAME_LEN) :: full_name       ! Full name of variable in varlist, including prefix and suffix
    INTEGER                    :: type            ! REAL2D / REAL3D
    INTEGER                    :: owner_proc_id = -1
    INTEGER, ALLOCATABLE       :: owner_tile_path(:)
    INTEGER, ALLOCATABLE       :: child_idx(:)    ! Position-numbers of same variable for for the child tiles (according to mem%vars(:)).
    REAL(dp)                   :: missval
    LOGICAL                    :: l_aggregate_all = .FALSE. ! When aggregating, consider fractions of all child tiles even
                                                  ! of those where variable doesn't exist
  END TYPE t_jsb_var

!#pragma GCC diagnostic push
!#pragma GCC diagnostic ignored "-Wsurprising"


  !-----------------------------------------------------------------------------------------------------
  !> 2D REAL Variable TYPE
  !! 
  !! ics:ice,iblk
  !-----------------------------------------------------------------------------------------------------
  TYPE, EXTENDS(t_jsb_var) :: t_jsb_var_real2d
    REAL(wp), POINTER :: ptr(:,:) => NULL()
  CONTAINS
    FINAL     :: finalize_jsb_var_real2d
  END TYPE t_jsb_var_real2d

  !-----------------------------------------------------------------------------------------------------
  !> 3D REAL Variable TYPE
  !! 
  !! ics:ice,:,iblk
  !-----------------------------------------------------------------------------------------------------
  TYPE, EXTENDS(t_jsb_var) :: t_jsb_var_real3d
    REAL(wp), POINTER :: ptr(:,:,:) => NULL()
  CONTAINS
    FINAL     :: finalize_jsb_var_real3d
  END TYPE t_jsb_var_real3d

  !-----------------------------------------------------------------------------------------------------
  !> Pointer to the abstract Variable TYPE
  !! 
  !! 
  !-----------------------------------------------------------------------------------------------------
  TYPE t_jsb_var_p
    CLASS(t_jsb_var), POINTER :: p
  END TYPE t_jsb_var_p

! quincy
  !-----------------------------------------------------------------------------------------------------
  !> TYPE that contains multiple variables
  !! 
  !! this type is aimed to be either included in: \n 
  !!    (1) the process_memory directly (mem\%TYPE) or \n
  !!    (2) a TYPE t_jsb_n_struct
  !! 
  !-----------------------------------------------------------------------------------------------------
  TYPE, ABSTRACT :: t_jsb_struct
    TYPE(t_jsb_var_p), ALLOCATABLE :: vars(:)
    CHARACTER(len=VARNAME_LEN)     :: name                      !< Name of variable
    CHARACTER(len=VARNAME_LEN)     :: full_name                 !< Full name of variable in varlist, including prefix and suffix
    INTEGER                        :: owner_proc_id      = -1   !< owner process id
    INTEGER, ALLOCATABLE           :: owner_tile_path(:)        !< owner tile path
    INTEGER                        :: no_of_vars         = 0    !< number of variables the TYPE contains in vars(:)
    INTEGER                        :: max_no_of_vars     = 20   !< ! hard coded limit of max var number - might be changed later !
    CHARACTER(LEN=50)              :: varlist_name              !< name of variable list
  CONTAINS
    PROCEDURE                      :: t_jsb_struct_quincy_add_struct_element_r2d
    PROCEDURE                      :: t_jsb_struct_quincy_add_struct_element_r3d
    GENERIC                        :: Add_var => t_jsb_struct_quincy_add_struct_element_r2d, &
                                                 t_jsb_struct_quincy_add_struct_element_r3d
  END TYPE t_jsb_struct

! quincy
  !-----------------------------------------------------------------------------------------------------
  !> POINTER to the t_jsb_struct
  !! 
  !-----------------------------------------------------------------------------------------------------
  TYPE t_jsb_struct_p
    CLASS(t_jsb_struct), POINTER :: p
  END TYPE t_jsb_struct_p

! quincy
  !-----------------------------------------------------------------------------------------------------
  !> TYPE that contains multiple TYPEs of t_jsb_struct
  !! 
  !! this type is aimed to be included in the process_memory directly (mem\%TYPE)
  !! 
  !-----------------------------------------------------------------------------------------------------
  TYPE, ABSTRACT :: t_jsb_n_struct
    TYPE(t_jsb_struct_p), ALLOCATABLE :: l2_structs(:)              !< 
    CHARACTER(len=VARNAME_LEN)        :: name                       !< Name of variable
    CHARACTER(len=VARNAME_LEN)        :: full_name                  !< Full name of variable in varlist, including prefix and suffix
    INTEGER                           :: owner_proc_id        = -1  !< 
    INTEGER, ALLOCATABLE              :: owner_tile_path(:)         !< 
    INTEGER                           :: no_of_l2_structs     = 0   !< 
    INTEGER                           :: max_no_of_l2_structs = 20  !<  ! hard coded limit - might be changed later !
    CHARACTER(LEN=50)                 :: l2_structlist_name         !< 
  CONTAINS
    PROCEDURE                         :: Add_pool => t_jsb_n_struct_quincy_add_level2_struct   ! mem%nested-struct%level2-struct
  END TYPE t_jsb_n_struct

! quincy
  !-----------------------------------------------------------------------------------------------------
  !> POINTER to the t_jsb_n_struct
  !! 
  !-----------------------------------------------------------------------------------------------------
  TYPE t_jsb_n_struct_p
    CLASS(t_jsb_n_struct), POINTER :: p
  END TYPE t_jsb_n_struct_p

!!$  INTERFACE t_jsb_var
!!$    PROCEDURE :: Create_jsb_var
!!$  END INTERFACE t_jsb_var

  CHARACTER(len=*), PARAMETER :: modname = 'mo_jsb_var_class'


CONTAINS


  ! ======================================================================================================= !
  !>
  !>
  SUBROUTINE finalize_jsb_var_real2d(this)

    TYPE(t_jsb_var_real2d) :: this

    CHARACTER(len=*), PARAMETER :: routine = modname//':finalize_jsb_var_real2d'

    !CALL message(routine, 'Finalizing '//TRIM(this%name))

    IF (ALLOCATED(this%owner_tile_path)) DEALLOCATE(this%owner_tile_path)
    IF (ALLOCATED(this%child_idx)) DEALLOCATE(this%child_idx)

    IF (ASSOCIATED(this%ptr)) NULLIFY(this%ptr)
    IF (ALLOCATED(this%field)) DEALLOCATE(this%field)

  END SUBROUTINE finalize_jsb_var_real2d


  ! ======================================================================================================= !
  !>
  !>
  SUBROUTINE finalize_jsb_var_real3d(this)

    TYPE(t_jsb_var_real3d) :: this

    CHARACTER(len=*), PARAMETER :: routine = modname//':finalize_jsb_var_real3d'

    !CALL message(routine, 'Finalizing '//TRIM(this%name))

    IF (ALLOCATED(this%owner_tile_path)) DEALLOCATE(this%owner_tile_path)
    IF (ALLOCATED(this%child_idx)) DEALLOCATE(this%child_idx)

    IF (ASSOCIATED(this%ptr)) NULLIFY(this%ptr)
    IF (ALLOCATED(this%field)) DEALLOCATE(this%field)

  END SUBROUTINE finalize_jsb_var_real3d

  
  ! ======================================================================================================= !
  !> add 2D variable to pool structure
  !>
  !> mem%pool <- var  or  mem%pool%pool <- var
  !>
  SUBROUTINE t_jsb_struct_quincy_add_struct_element_r2d( &
    struct, name, var, &
    hgrid, vgrid, cf, grib1, grib2, prefix, suffix, &
    loutput, lrestart, output_level, &
    initval_r)

    CLASS(t_jsb_struct),              INTENT(inout)   :: struct         ! pool structure
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
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':t_jsb_struct_quincy_add_struct_element_r2d'

    ics     = 1              ! is always 1, see land_model_driver
    ice     = hgrid%nproma   ! number of points in the chunk
    nproma  = hgrid%nproma   ! number of points in the chunk
    iblk    = 1              ! is always 1, see land_model_driver

    ALLOCATE(var%field(nproma,iblk,1,1,1))
    var%ptr => var%field(:,:,1,1,1)

    var%name          = TRIM(name)
    var%owner_proc_id = struct%owner_proc_id

    pos_of_var = struct%no_of_vars + 1
    IF (pos_of_var > struct%max_no_of_vars) THEN
      CALL finish(TRIM(routine), 'Maximum number of vars exceeded: '//TRIM(struct%varlist_name)//', '//TRIM(name))
    ELSE
      struct%vars(pos_of_var)%p => var
      struct%no_of_vars = pos_of_var
    END IF

    ! added for QUINCY compatibility (in JSBACH4 the allocation is done differently)
    var%ptr(ics:ice,iblk) = initval_r

  END SUBROUTINE t_jsb_struct_quincy_add_struct_element_r2d


  ! ======================================================================================================= !
  !> add 3D variable to pool structure
  !>
  !> mem%pool <- var  or  mem%pool%pool <- var
  !> 
  SUBROUTINE t_jsb_struct_quincy_add_struct_element_r3d( &
    struct, name, var, &
    hgrid, vgrid, cf, grib1, grib2, prefix, suffix, &
    loutput, lrestart, output_level, &
    initval_r)
  
    CLASS(t_jsb_struct),              INTENT(inout)   :: struct         ! pool structure
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
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':t_jsb_struct_quincy_add_struct_element_r3d'
  
    ics     = 1              ! is always 1, see land_model_driver
    ice     = hgrid%nproma   ! number of points in the chunk
    nproma  = hgrid%nproma   ! number of points in the chunk
    iblk    = 1              ! is always 1, see land_model_driver
  
    ALLOCATE(var%field(nproma, vgrid%n_levels, iblk, 1, 1))
    var%ptr => var%field(:,:,:,1,1)
  
    var%name          = TRIM(name)
    var%owner_proc_id = struct%owner_proc_id
  
    pos_of_var = struct%no_of_vars + 1
    IF (pos_of_var > struct%max_no_of_vars) THEN
      CALL finish(TRIM(routine), 'Maximum number of vars exceeded: '//TRIM(struct%varlist_name)//', '//TRIM(name))
    ELSE
      struct%vars(pos_of_var)%p => var
      struct%no_of_vars = pos_of_var
    END IF
  
    ! added for QUINCY compatibility (in JSBACH4 the allocation is done differently)
    var%ptr(ics:ice, 1:vgrid%n_levels, iblk) = initval_r
  
  END SUBROUTINE t_jsb_struct_quincy_add_struct_element_r3d


  ! ======================================================================================================= !
  !> add pool to a nested pool structure
  !>
  !> mem%pool <- pool  (add a pool that contains only var, but no other pools)
  !> 
  SUBROUTINE t_jsb_n_struct_quincy_add_level2_struct( &
    n_struct, name, l2_struct, &
    hgrid, vgrid, cf, grib1, grib2, prefix, suffix, &
    loutput, lrestart, output_level)

    CLASS(t_jsb_n_struct),         INTENT(inout) :: n_struct            ! memory structure
    CHARACTER(len=*),              INTENT(in)    :: name                ! name of variable
    CLASS(t_jsb_struct),   TARGET, INTENT(inout) :: l2_struct           ! l2_struct
    TYPE(t_jsb_grid),     POINTER, INTENT(in)    :: hgrid               ! horizontal grid type
    TYPE(t_jsb_vgrid),    POINTER, INTENT(in)    :: vgrid               ! vertical grid type
    TYPE(t_cf),                    INTENT(in)    :: cf                  ! CF related metadata
    TYPE(t_grib1),                 INTENT(in)    :: grib1               ! GRIB1 related metadata
    TYPE(t_grib2),                 INTENT(in)    :: grib2               ! GRIB2 related metadata
    CHARACTER(len=*),              INTENT(in)    :: prefix              ! Prefix for variable name
    CHARACTER(len=*),              INTENT(in)    :: suffix              ! Suffix for variable name
    LOGICAL,             OPTIONAL, INTENT(in)    :: loutput             ! output flag
    LOGICAL,             OPTIONAL, INTENT(in)    :: lrestart            ! restart flag
    INTEGER,             OPTIONAL, INTENT(in)    :: output_level        ! minimum level of output groups this var belongs to

    INTEGER                     :: pos_of_l2_struct
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':t_jsb_n_struct_quincy_add_level2_struct'

    l2_struct%name          = TRIM(name)
    l2_struct%owner_proc_id = n_struct%owner_proc_id

    pos_of_l2_struct = n_struct%no_of_l2_structs + 1
    IF (pos_of_l2_struct > n_struct%max_no_of_l2_structs) THEN
      CALL finish(TRIM(routine), 'Maximum number of l2_structs exceeded: '//TRIM(n_struct%l2_structlist_name)//', '//TRIM(name))
    ELSE
      n_struct%l2_structs(pos_of_l2_struct)%p => l2_struct
      n_struct%no_of_l2_structs = pos_of_l2_struct
    END IF

    ALLOCATE(l2_struct%vars(l2_struct%max_no_of_vars))  !!! max_no_of_vars is hardcoded in mo_jsb_var_class:t_jsb_struct

  END SUBROUTINE t_jsb_n_struct_quincy_add_level2_struct


#endif
END MODULE mo_jsb_var_class
