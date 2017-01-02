!>
!! @brief 
!!  In SLM: contains a small version of mo_jsb_varlist
!!  In JSBACH4: Contains routines for var_lists/streams
!!
!! @author
!!  Jan Engel, MPI-BGC Jena
!!
!! @par Revision History
!! First version by Jan Engel (20160816)
!!
MODULE mo_slm_varlist

  USE mo_kind,              ONLY: wp
  USE mo_exception,         ONLY: finish
  USE mo_jsb_class,         ONLY: Get_model   ! *slmtemp 
  USE mo_base_memory,       ONLY: t_model     ! *slmtemp  -> making available ice ics iblk, which is not needed in jsbach4

  USE mo_jsb_var_class,     ONLY: t_jsb_var_real2d, t_jsb_var_real3d, &
                                  t_jsb_var_logical2d                        , &                 ! *slmjsb4
                                  t_jsb_struct, t_jsb_n_struct               , &                 ! *slmjsb4
                                  t_jsb_var_real1d, t_jsb_var_logical1d                          ! *slmtemp
  USE mo_jsb_memory_class,  ONLY: t_jsb_memory

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: slm_add_var

  INTERFACE slm_add_var
  !! struct
    ! mem%struct 
    MODULE PROCEDURE slm_add_struct               ! *slmjsb4  a TYPE with only variables included
    ! mem%nested-struct
    MODULE PROCEDURE slm_add_n_struct             ! *slmjsb4  a TYPE with nested TYPEs
    ! mem%nested-struct%level2-struct
    MODULE PROCEDURE slm_add_l2_struct            ! *slmjsb4  a TYPE with only variables included within a TYPE
  !! variable -> struct (mem%struct or mem%nested-struct%level2-struct | works NOT for mem%nested-struct%VAR !)
    ! REAL 1D = scalar!
    MODULE PROCEDURE slm_add_struct_element_r1d   ! *slmjsb4 *slmtemp
    ! REAL 2D
    MODULE PROCEDURE slm_add_struct_element_r2d   ! *slmjsb4 *slmtemp
    ! REAL 3D
    MODULE PROCEDURE slm_add_struct_element_r3d   ! *slmjsb4 *slmtemp
  !! variable -> mem (process memory)
    ! REAL 1D = scalar!
    MODULE PROCEDURE slm_add_var_list_element_r1d ! *slmtemp
    ! LOGICAL 1D = scalar!
    MODULE PROCEDURE slm_add_var_list_element_l1d ! *slmjsb4 *slmtemp
    ! REAL 2D
    MODULE PROCEDURE slm_add_var_list_element_r2d ! *slmtemp
    ! REAL 3D
    MODULE PROCEDURE slm_add_var_list_element_r3d ! *slmtemp
  END INTERFACE slm_add_var

  CHARACTER(len=*), PARAMETER :: modname = 'mo_slm_varlist'

CONTAINS


  !-----------------------------------------------------------------------------------------------------!
  ! add struct                                                                                          !
  !-----------------------------------------------------------------------------------------------------!

  !-----------------------------------------------------------------------------------------------------
  !> add structure to process memory
  !!
  !!  t_jsb_memory:mem <- t_jsb_struct:struct
  !!    
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE slm_add_struct( &
    mem, name, struct) 
    IMPLICIT NONE
    ! 0.1 InOut
    CLASS(t_jsb_memory),  INTENT(inout)           :: mem                 ! memory structure
    CHARACTER(len=*),     INTENT(in)              :: name                ! name of variable
    CLASS(t_jsb_struct), TARGET, INTENT(inout)    :: struct              ! struct
    ! 0.2 Local
    INTEGER           :: pos_of_struct
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':slm_add_struct'

    struct%name = TRIM(name)
    struct%owner_proc_id = mem%owner_proc_id  ! ATM2LAND_=1, LAND2ATM_=2, SOIL_=5, VEG_=7, RAD_=9, ASSIMI_=10

    pos_of_struct = mem%no_of_structs + 1
    IF (pos_of_struct > mem%max_no_of_structs) THEN
      CALL finish(TRIM(routine), 'Maximum number of structs exceeded: '//TRIM(mem%structlist_name)//', '//TRIM(name))
    ELSE
      mem%structs(pos_of_struct)%p => struct
      mem%no_of_structs = pos_of_struct
    END IF

    ALLOCATE(struct%vars(struct%max_no_of_vars))  !!! hardcoded in   TYPE, ABSTRACT :: t_jsb_struct  ! *slmtemp

  END SUBROUTINE slm_add_struct


  !-----------------------------------------------------------------------------------------------------
  !> add nested structure to process memory
  !!
  !!  t_jsb_memory:mem <- t_jsb_n_struct:n_struct
  !!    
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE slm_add_n_struct( &
    mem, name, n_struct) 
    IMPLICIT NONE
    ! 0.1 InOut
    CLASS(t_jsb_memory),  INTENT(inout)           :: mem                 ! memory structure
    CHARACTER(len=*),     INTENT(in)              :: name                ! name of variable
    CLASS(t_jsb_n_struct), TARGET, INTENT(inout)  :: n_struct            ! n_struct
    ! 0.2 Local
    INTEGER           :: pos_of_n_struct
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':slm_add_n_struct'

    n_struct%name = TRIM(name)
    n_struct%owner_proc_id = mem%owner_proc_id  ! ATM2LAND_=1, LAND2ATM_=2, SOIL_=5, VEG_=7, RAD_=9, ASSIMI_=10

    pos_of_n_struct = mem%no_of_n_structs + 1
    IF (pos_of_n_struct > mem%max_no_of_n_structs) THEN
      CALL finish(TRIM(routine), 'Maximum number of n_structs exceeded: '//TRIM(mem%n_structlist_name)//', '//TRIM(name))
    ELSE
      mem%n_structs(pos_of_n_struct)%p => n_struct
      mem%no_of_n_structs = pos_of_n_struct
    END IF

    ALLOCATE(n_struct%l2_structs(n_struct%max_no_of_l2_structs))  !!! hardcoded in   TYPE, ABSTRACT :: t_jsb_n_struct  ! *slmtemp

  END SUBROUTINE slm_add_n_struct


  !-----------------------------------------------------------------------------------------------------
  !> add structure to a nested structure that is located directly in the process memory
  !!
  !!  t_jsb_n_struct:n_struct <- t_jsb_struct:l2_struct
  !!    
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE slm_add_l2_struct( &
    n_struct, name, l2_struct) 
    IMPLICIT NONE
    ! 0.1 InOut
    CLASS(t_jsb_n_struct),  INTENT(inout)       :: n_struct            ! memory structure
    CHARACTER(len=*),     INTENT(in)            :: name                ! name of variable
    CLASS(t_jsb_struct), TARGET, INTENT(inout)  :: l2_struct           ! l2_struct
    ! 0.2 Local
    INTEGER           :: pos_of_l2_struct
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':slm_add_l2_struct'

    l2_struct%name = TRIM(name)
    l2_struct%owner_proc_id = n_struct%owner_proc_id  ! ATM2LAND_=1, LAND2ATM_=2, SOIL_=5, VEG_=7, RAD_=9, ASSIMI_=10

    pos_of_l2_struct = n_struct%no_of_l2_structs + 1
    IF (pos_of_l2_struct > n_struct%max_no_of_l2_structs) THEN
      CALL finish(TRIM(routine), 'Maximum number of l2_structs exceeded: '//TRIM(n_struct%l2_structlist_name)//', '//TRIM(name))
    ELSE
      n_struct%l2_structs(pos_of_l2_struct)%p => l2_struct
      n_struct%no_of_l2_structs = pos_of_l2_struct
    END IF

    ALLOCATE(l2_struct%vars(l2_struct%max_no_of_vars))  !!! hardcoded in   TYPE, ABSTRACT :: t_jsb_struct  ! *slmtemp

  END SUBROUTINE slm_add_l2_struct



  !-----------------------------------------------------------------------------------------------------!
  ! add variable to struct                                                                              !
  !-----------------------------------------------------------------------------------------------------!

  !-----------------------------------------------------------------------------------------------------
  !> add 1D variable to struct or l2_struct
  !!
  !!  t_jsb_struct:struct <- t_jsb_var_real1d:var
  !! 
  !! SCALAR !
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE slm_add_struct_element_r1d( &
    struct, name, var, &
    initval_r) 
    IMPLICIT NONE
    ! 0.1 InOut
    CLASS(t_jsb_struct),              INTENT(inout)   :: struct              ! struct & l2_struct
    CHARACTER(len=*),                 INTENT(in)      :: name                ! name of variable
    TYPE(t_jsb_var_real1d), TARGET,   INTENT(inout)   :: var                 ! reference to jsbach var -> newly created plain scalar REAL
    REAL(wp),                         INTENT(in)      :: initval_r           ! value for initialization of the variable (usually 0.0_wp)
    ! 0.2 Local
    INTEGER                       :: pos_of_var
    TYPE(t_model), POINTER        :: model            ! *slmtemp
    INTEGER                       :: ice, ics, iblk   ! *slmtemp
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':slm_add_struct_element_real1d'

    model => Get_model(1)  !     model => Get_model(model_id)   *slmtemp hardcoded model ID !!!
    ice  = model%options%ice    ! *slmtemp
    ics  = model%options%ics    ! *slmtemp
    iblk = model%options%iblk   ! *slmtemp

    var%name = TRIM(name)
    var%owner_proc_id = struct%owner_proc_id  ! ATM2LAND_=1, LAND2ATM_=2, SOIL_=5, VEG_=7, RAD_=9, ASSIMI_=10

    pos_of_var = struct%no_of_vars + 1
    IF (pos_of_var > struct%max_no_of_vars) THEN
      CALL finish(TRIM(routine), 'Maximum number of vars exceeded: '//TRIM(struct%varlist_name)//', '//TRIM(name))
    ELSE
      struct%vars(pos_of_var)%p => var
      struct%no_of_vars = pos_of_var
    END IF

    ! added for SLM compatibility (in JSBACH4 the allocation is done differently)   *slmtemp
    ALLOCATE(var%ptr)
    var%ptr = initval_r

  END SUBROUTINE slm_add_struct_element_r1d


  !-----------------------------------------------------------------------------------------------------
  !> add 2D variable to struct or l2_struct
  !!
  !!  t_jsb_struct:struct <- t_jsb_var_real2d:var
  !!  
  !! ==> currectly used correctly for SOIL 2D
  !!  AND
  !! ==> fake/transitionally for VEG 3D variables
  !! 
  !!  see local variables levels:levele
  !! 
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE slm_add_struct_element_r2d( &
    struct, name, var, &
    initval_r) 

    USE mo_soil_config_class,               ONLY: soil_config
    USE mo_assimi_config_class,             ONLY: assimi_config

    IMPLICIT NONE
    ! 0.1 InOut
    CLASS(t_jsb_struct),              INTENT(inout)   :: struct              ! struct & l2_struct
    CHARACTER(len=*),                 INTENT(in)      :: name                ! name of variable
    TYPE(t_jsb_var_real2d), TARGET,   INTENT(inout)   :: var                 ! reference to jsbach var -> newly created plain scalar REAL
    REAL(wp),                         INTENT(in)      :: initval_r           ! value for initialization of the variable (usually 0.0_wp)
    ! 0.2 Local
    INTEGER                       :: pos_of_var
    TYPE(t_model), POINTER        :: model            ! *slmtemp
    INTEGER                       :: ice, ics, iblk   ! *slmtemp
    INTEGER                       :: levels, levele   ! level_start & level_end | SOIL ics:ice && VEG 1:ncanopy
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':slm_add_struct_element_real2d'

    model => Get_model(1)  !     model => Get_model(model_id)   *slmtemp hardcoded model ID !!!
    ice  = model%options%ice    ! *slmtemp
    ics  = model%options%ics    ! *slmtemp
    iblk = model%options%iblk   ! *slmtemp

    ! select nlevel according to process id 
    SELECT CASE(struct%owner_proc_id)
      CASE(5) ! soil
        levels = ics
        levele = ice
      CASE(7, 9, 10) ! veg, rad, assimi
        levels = 1
        levele = assimi_config%ncanopy
      CASE DEFAULT
        CALL finish(TRIM(routine), 'A variable that is not in memory of soil, veg or assimi aims beeing 2D, '//TRIM(name))
    END SELECT

    var%name = TRIM(name)
    var%owner_proc_id = struct%owner_proc_id    ! ATM2LAND_=1, LAND2ATM_=2, SOIL_=5, VEG_=7, RAD_=9, ASSIMI_=10

    pos_of_var = struct%no_of_vars + 1
    IF (pos_of_var > struct%max_no_of_vars) THEN
      CALL finish(TRIM(routine), 'Maximum number of vars exceeded: '//TRIM(struct%varlist_name)//', '//TRIM(name))
    ELSE
      struct%vars(pos_of_var)%p => var
      struct%no_of_vars = pos_of_var
    END IF

    ! added for SLM compatibility (in JSBACH4 the allocation is done differently)   *slmtemp
    ALLOCATE(var%ptr(levels:levele,iblk))
    var%ptr(levels:levele,iblk) = initval_r

  END SUBROUTINE slm_add_struct_element_r2d


  !-----------------------------------------------------------------------------------------------------
  !> add 3D variable to struct or l2_struct
  !!
  !!  t_jsb_struct:struct <- t_jsb_var_real2d:var
  !!    
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE slm_add_struct_element_r3d( &
    struct, name, var, &
    initval_r) 

    USE mo_soil_config_class,               ONLY: soil_config
    USE mo_assimi_config_class,             ONLY: assimi_config

    IMPLICIT NONE
    ! 0.1 InOut
    CLASS(t_jsb_struct),              INTENT(inout)   :: struct              ! struct & l2_struct
    CHARACTER(len=*),                 INTENT(in)      :: name                ! name of variable
    TYPE(t_jsb_var_real3d), TARGET,   INTENT(inout)   :: var                 ! reference to jsbach var -> newly created plain scalar REAL
    REAL(wp),                         INTENT(in)      :: initval_r           ! value for initialization of the variable (usually 0.0_wp)
    ! 0.2 Local
    INTEGER                       :: pos_of_var
    TYPE(t_model), POINTER        :: model            ! *slmtemp
    INTEGER                       :: ice, ics, iblk   ! *slmtemp
    INTEGER                       :: nlevel           ! number of levels at Z-Axis (canopy & soil)
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':slm_add_struct_element_real2d'

    model => Get_model(1)  !     model => Get_model(model_id)   *slmtemp hardcoded model ID !!!
    ice  = model%options%ice    ! *slmtemp
    ics  = model%options%ics    ! *slmtemp
    iblk = model%options%iblk   ! *slmtemp

    ! select nlevel according to process id 
    SELECT CASE(struct%owner_proc_id)
      CASE(5) ! soil
        nlevel = soil_config%nsoil_water
      CASE(7, 9, 10) ! veg, rad, assimi
        nlevel = assimi_config%ncanopy
      CASE DEFAULT
        CALL finish(TRIM(routine), 'A variable that is not in memory of soil, veg or assimi aims beeing 3D, '//TRIM(name))
    END SELECT

    var%name = TRIM(name)
    var%owner_proc_id = struct%owner_proc_id    ! ATM2LAND_=1, LAND2ATM_=2, SOIL_=5, VEG_=7, RAD_=9, ASSIMI_=10

    pos_of_var = struct%no_of_vars + 1
    IF (pos_of_var > struct%max_no_of_vars) THEN
      CALL finish(TRIM(routine), 'Maximum number of vars exceeded: '//TRIM(struct%varlist_name)//', '//TRIM(name))
    ELSE
      struct%vars(pos_of_var)%p => var
      struct%no_of_vars = pos_of_var
    END IF

    ! added for SLM compatibility (in JSBACH4 the allocation is done differently)   *slmtemp
    ALLOCATE(var%ptr(ics:ice,1:nlevel,iblk))
    var%ptr(ics:ice,1:nlevel,iblk) = initval_r

  END SUBROUTINE slm_add_struct_element_r3d



  !-----------------------------------------------------------------------------------------------------!
  ! add variable to process memory (mem%VAR)                                                            !
  !-----------------------------------------------------------------------------------------------------!

  !-----------------------------------------------------------------------------------------------------
  !> add 1D variable to to process memory
  !!
  !!  t_jsb_memory:mem <- t_jsb_var_r1d:var
  !! 
  !! SCALAR !
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE slm_add_var_list_element_r1d( &
    mem, name, var, &
    initval_r) 
    IMPLICIT NONE
    ! 0.1 InOut
    CLASS(t_jsb_memory),              INTENT(inout)   :: mem                 ! memory structure
    CHARACTER(len=*),                 INTENT(in)      :: name                ! name of variable
    TYPE(t_jsb_var_real1d), TARGET,   INTENT(inout)   :: var                 ! reference to jsbach var -> newly created plain scalar REAL
    REAL(wp),                         INTENT(in)      :: initval_r           ! value for initialization of the variable (usually 0.0_wp)
    ! 0.2 Local
    INTEGER                       :: pos_of_var
    TYPE(t_model), POINTER        :: model            ! *slmtemp
    INTEGER                       :: ice, ics, iblk   ! *slmtemp
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':slm_add_var_list_element_r1d'

    model => Get_model(1)  !     model => Get_model(model_id)   *slmtemp hardcoded model ID !!!
    ice  = model%options%ice    ! *slmtemp
    ics  = model%options%ics    ! *slmtemp
    iblk = model%options%iblk   ! *slmtemp

    var%name = TRIM(name)
    var%owner_proc_id = mem%owner_proc_id  ! ATM2LAND_=1, LAND2ATM_=2, SOIL_=5, VEG_=7, RAD_=9, ASSIMI_=10

    pos_of_var = mem%no_of_vars + 1
    IF (pos_of_var > mem%max_no_of_vars) THEN
      CALL finish(TRIM(routine), 'Maximum number of vars exceeded: '//TRIM(mem%varlist_name)//', '//TRIM(name))
    ELSE
      mem%vars(pos_of_var)%p => var
      mem%no_of_vars = pos_of_var
    END IF

    ! added for SLM compatibility (in JSBACH4 the allocation is done differently)   *slmtemp
    ALLOCATE(var%ptr)
    var%ptr = initval_r

  END SUBROUTINE slm_add_var_list_element_r1d


  !-----------------------------------------------------------------------------------------------------
  !> add 1D LOGICAL variable to process memory
  !!
  !!  t_jsb_memory:mem <- t_jsb_var_logical1d:var
  !!    
  !! SCALAR !
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE slm_add_var_list_element_l1d( &
    mem, name, var, &
    initval_r) 
    IMPLICIT NONE
    ! 0.1 InOut
    CLASS(t_jsb_memory),                INTENT(inout)   :: mem              ! memory structure
    CHARACTER(len=*),                   INTENT(in)      :: name             ! name of variable
    TYPE(t_jsb_var_logical1d), TARGET,  INTENT(inout)   :: var              ! reference to jsbach var -> newly created plain scalar REAL
    LOGICAL,                            INTENT(in)      :: initval_r        ! value for initialization of the variable (usually 0.0_wp)
    ! 0.2 Local
    INTEGER                   :: pos_of_var
    TYPE(t_model), POINTER    :: model            ! *slmtemp
    INTEGER                   :: ice, ics, iblk   ! *slmtemp
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':slm_add_var_list_element_l1d'

    model => Get_model(1)  !     model => Get_model(model_id)   *slmtemp hardcoded model ID !!!
    ice  = model%options%ice    ! *slmtemp
    ics  = model%options%ics    ! *slmtemp
    iblk = model%options%iblk   ! *slmtemp

    var%name = TRIM(name)
    var%owner_proc_id = mem%owner_proc_id

    pos_of_var = mem%no_of_vars + 1
    IF (pos_of_var > mem%max_no_of_vars) THEN
      CALL finish(TRIM(routine), 'Maximum number of vars exceeded: '//TRIM(mem%varlist_name)//', '//TRIM(name))
    ELSE
      mem%vars(pos_of_var)%p => var
      mem%no_of_vars = pos_of_var
    END IF

    ! added for SLM compatibility (in JSBACH4 the allocation is done differently)
    ALLOCATE(var%ptr)
    var%ptr = initval_r

  END SUBROUTINE slm_add_var_list_element_l1d


  !-----------------------------------------------------------------------------------------------------
  !> add 2D variable to process memory
  !!
  !!  t_jsb_memory:mem <- t_jsb_var_r2d:var
  !! 
  !!  
  !! ==> currectly used correctly for SOIL 2D
  !!  AND
  !! ==> fake/transitionally for VEG 3D variables
  !! 
  !!  see local variables levels:levele
  !! 
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE slm_add_var_list_element_r2d( &
    mem, name, var, &
    initval_r) 

    USE mo_soil_config_class,               ONLY: soil_config
    USE mo_assimi_config_class,             ONLY: assimi_config

    IMPLICIT NONE
    ! 0.1 InOut
    CLASS(t_jsb_memory),              INTENT(inout)   :: mem                 ! memory structure
    CHARACTER(len=*),                 INTENT(in)      :: name                ! name of variable
    TYPE(t_jsb_var_real2d), TARGET,   INTENT(inout)   :: var                 ! reference to jsbach var -> newly created plain scalar REAL
    REAL(wp),                         INTENT(in)      :: initval_r           ! value for initialization of the variable (usually 0.0_wp)
    ! 0.2 Local
    INTEGER                       :: pos_of_var
    TYPE(t_model), POINTER        :: model            ! *slmtemp
    INTEGER                       :: ice, ics, iblk   ! *slmtemp
    INTEGER                       :: levels, levele   ! level_start & level_end | SOIL ics:ice && VEG 1:ncanopy
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':slm_add_var_list_element_r2d'

    model => Get_model(1)  !     model => Get_model(model_id)   *slmtemp hardcoded model ID !!!
    ice  = model%options%ice    ! *slmtemp
    ics  = model%options%ics    ! *slmtemp
    iblk = model%options%iblk   ! *slmtemp

    ! select nlevel according to process id 
    SELECT CASE(mem%owner_proc_id)
      CASE(5) ! soil
        levels = ics
        levele = ice
      CASE(7, 9, 10) ! veg, rad, assimi
        levels = 1
        levele = assimi_config%ncanopy
      CASE DEFAULT
        CALL finish(TRIM(routine), 'A variable that is not in memory of soil, veg or assimi aims beeing 2D, '//TRIM(name))
    END SELECT

    var%name = TRIM(name)
    var%owner_proc_id = mem%owner_proc_id    ! ATM2LAND_=1, LAND2ATM_=2, SOIL_=5, VEG_=7, RAD_=9, ASSIMI_=10

    pos_of_var = mem%no_of_vars + 1
    IF (pos_of_var > mem%max_no_of_vars) THEN
      CALL finish(TRIM(routine), 'Maximum number of vars exceeded: '//TRIM(mem%varlist_name)//', '//TRIM(name))
    ELSE
      mem%vars(pos_of_var)%p => var
      mem%no_of_vars = pos_of_var
    END IF

    ! added for SLM compatibility (in JSBACH4 the allocation is done differently)   *slmtemp
    ALLOCATE(var%ptr(levels:levele,iblk))
    var%ptr(levels:levele,iblk) = initval_r

  END SUBROUTINE slm_add_var_list_element_r2d


  !-----------------------------------------------------------------------------------------------------
  !> add 3D variable to process memory
  !!
  !!  t_jsb_memory:mem <- t_jsb_var_r3d:var
  !! 
  !!    
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE slm_add_var_list_element_r3d( &
    mem, name, var, &
    initval_r) 

    USE mo_soil_config_class,               ONLY: soil_config
    USE mo_assimi_config_class,             ONLY: assimi_config

    IMPLICIT NONE
    ! 0.1 InOut
    CLASS(t_jsb_memory),              INTENT(inout)   :: mem                 ! memory structure
    CHARACTER(len=*),                 INTENT(in)      :: name                ! name of variable
    TYPE(t_jsb_var_real3d), TARGET,   INTENT(inout)   :: var                 ! reference to jsbach var -> newly created plain scalar REAL
    REAL(wp),                         INTENT(in)      :: initval_r           ! value for initialization of the variable (usually 0.0_wp)
    ! 0.2 Local
    INTEGER                       :: pos_of_var
    TYPE(t_model), POINTER        :: model            ! *slmtemp
    INTEGER                       :: ice, ics, iblk   ! *slmtemp
    INTEGER                       :: nlevel           ! number of levels at Z-Axis (canopy & soil)
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':slm_add_var_list_element_r3d'

    model => Get_model(1)  !     model => Get_model(model_id)   *slmtemp hardcoded model ID !!!
    ice  = model%options%ice    ! *slmtemp
    ics  = model%options%ics    ! *slmtemp
    iblk = model%options%iblk   ! *slmtemp

    ! select nlevel according to process id 
    SELECT CASE(mem%owner_proc_id)
      CASE(5) ! soil
        nlevel = soil_config%nsoil_water
      CASE(7, 9, 10) ! veg, rad, assimi
        nlevel = assimi_config%ncanopy
      CASE DEFAULT
        CALL finish(TRIM(routine), 'A variable that is not in memory of soil, veg or assimi aims beeing 3D, '//TRIM(name))
    END SELECT

    var%name = TRIM(name)
    var%owner_proc_id = mem%owner_proc_id    ! ATM2LAND_=1, LAND2ATM_=2, SOIL_=5, VEG_=7, RAD_=9, ASSIMI_=10

    pos_of_var = mem%no_of_vars + 1
    IF (pos_of_var > mem%max_no_of_vars) THEN
      CALL finish(TRIM(routine), 'Maximum number of vars exceeded: '//TRIM(mem%varlist_name)//', '//TRIM(name))
    ELSE
      mem%vars(pos_of_var)%p => var
      mem%no_of_vars = pos_of_var
    END IF

    ! added for SLM compatibility (in JSBACH4 the allocation is done differently)   *slmtemp
    ALLOCATE(var%ptr(ics:ice,1:nlevel,iblk))
    var%ptr(ics:ice,1:nlevel,iblk) = initval_r

  END SUBROUTINE slm_add_var_list_element_r3d



END MODULE mo_slm_varlist
