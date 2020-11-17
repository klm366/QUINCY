!>
!!  @par Copyright
!!    QUINCY terrestrial biosphere model with coupled carbon, nitrogen, and phosphorus cycles
!!    Copyright (C) 2019 Max Planck Institute for Biogeochemistry (MPI-BGC), Jena, Germany 
!!    Contributors: Caldararu S., Engel J., Kern M., Pallandt M., Schnur R., Thum T., Yu L., Zaehle S.
!!
!!    This file is part of the QUINCY model.
!!    The QUINCY model is free software: you can redistribute it and/or modify
!!    it under the terms of the GNU General Public License version 3,
!!    as published by the Free Software Foundation.
!!
!!    The QUINCY model is distributed in the hope that it will be useful,
!!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!!    GNU General Public License for more details.
!!
!!    You should have received a copy of the GNU General Public License
!!    along with the QUINCY model.  If not, see <https://www.gnu.org/licenses/>.
!! 
!! @brief 
!!  Contains routines for adding variables to var_lists
!!
!! @author
!!  Jan Engel, MPI-BGC Jena
!!
!! @par Revision History
!! First version by Jan Engel (20160816)
!!
MODULE mo_slm_varlist

#include "dsl_jsbach.inc"

  USE mo_kind,              ONLY: wp
  USE mo_exception,         ONLY: finish
  USE mo_jsb_class,         ONLY: Get_model   ! *slmtemp 
  USE mo_jsb_model_class,   ONLY: t_jsb_model
  USE mo_jsb_process_class, ONLY: A2L_, L2A_, SEB_, TURB_, SSE_, HYDRO_, VEG_, HD_, RAD_, ASSIMI_, PHENO_, SB_

  USE mo_jsb_var_class,     ONLY: t_jsb_var_real2d, t_jsb_var_real3d        , &                  ! *slmjsb4
                                  t_jsb_struct, t_jsb_n_struct                                   ! *slmjsb4
  USE mo_jsb_memory_class,  ONLY: t_jsb_memory
  USE mo_jsb_grid_class,    ONLY: t_jsb_grid, t_jsb_vgrid
  USE mo_jsb_io,            ONLY: t_cf, t_grib1, t_grib2

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: jsb_add_var

  INTERFACE jsb_add_var
  !! struct
    ! mem%struct 
    MODULE PROCEDURE slm_add_struct               ! *slmjsb4  a TYPE with only variables included
    ! mem%nested-struct
    MODULE PROCEDURE slm_add_n_struct             ! *slmjsb4  a TYPE with nested TYPEs
    ! mem%nested-struct%level2-struct
    MODULE PROCEDURE slm_add_l2_struct            ! *slmjsb4  a TYPE with only variables included within a TYPE

  !! variable -> struct (mem%struct or mem%nested-struct%level2-struct | works NOT for mem%nested-struct%VAR !)
    ! REAL 2D
    MODULE PROCEDURE slm_add_struct_element_r2d   ! *slmjsb4 *slmtemp
    ! REAL 3D
    MODULE PROCEDURE slm_add_struct_element_r3d   ! *slmjsb4 *slmtemp

  !! variable -> mem (process memory)
    ! REAL 2D
    MODULE PROCEDURE slm_add_var_list_element_r2d ! *slmtemp
    ! REAL 3D
    MODULE PROCEDURE slm_add_var_list_element_r3d ! *slmtemp
  END INTERFACE jsb_add_var

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
    mem, name, struct, &
    hgrid, vgrid, cf, grib1, grib2, prefix, suffix)
    IMPLICIT NONE

    ! 0.1 InOut
    CLASS(t_jsb_memory),  INTENT(inout)           :: mem                 ! memory structure
    CHARACTER(len=*),     INTENT(in)              :: name                ! name of variable
    CLASS(t_jsb_struct), TARGET, INTENT(inout)    :: struct              ! struct

    TYPE(t_jsb_grid),  POINTER, INTENT(in) :: hgrid     ! horizontal grid type
    TYPE(t_jsb_vgrid), POINTER, INTENT(in) :: vgrid     ! vertical grid type
    TYPE(t_cf),                 INTENT(in) :: cf        ! CF related metadata
    TYPE(t_grib1),              INTENT(in) :: grib1     ! GRIB1 related metadata
    TYPE(t_grib2),              INTENT(in) :: grib2     ! GRIB2 related metadata
    CHARACTER(len=*),           INTENT(in) :: prefix    ! Prefix for variable name
    CHARACTER(len=*),           INTENT(in) :: suffix    ! Suffix for variable name

    ! 0.2 Local
    INTEGER           :: pos_of_struct
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':slm_add_struct'

    struct%name = TRIM(name)
    struct%owner_proc_id = mem%owner_proc_id  ! A2L_=1, L2A_=2, SSE_=5, VEG_=7, RAD_=9, ASSIMI_=10, PHENO_=11, SB_=12

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
    mem, name, n_struct, &
    hgrid, vgrid, cf, grib1, grib2, prefix, suffix)
    IMPLICIT NONE

    ! 0.1 InOut
    CLASS(t_jsb_memory),  INTENT(inout)           :: mem                 ! memory structure
    CHARACTER(len=*),     INTENT(in)              :: name                ! name of variable
    CLASS(t_jsb_n_struct), TARGET, INTENT(inout)  :: n_struct            ! n_struct

    TYPE(t_jsb_grid),  POINTER, INTENT(in) :: hgrid     ! horizontal grid type
    TYPE(t_jsb_vgrid), POINTER, INTENT(in) :: vgrid     ! vertical grid type
    TYPE(t_cf),                 INTENT(in) :: cf        ! CF related metadata
    TYPE(t_grib1),              INTENT(in) :: grib1     ! GRIB1 related metadata
    TYPE(t_grib2),              INTENT(in) :: grib2     ! GRIB2 related metadata
    CHARACTER(len=*),           INTENT(in) :: prefix    ! Prefix for variable name
    CHARACTER(len=*),           INTENT(in) :: suffix    ! Suffix for variable name

    ! 0.2 Local
    INTEGER           :: pos_of_n_struct
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':slm_add_n_struct'

    n_struct%name = TRIM(name)
    n_struct%owner_proc_id = mem%owner_proc_id  ! A2L_=1, L2A_=2, SSE_=5, VEG_=7, RAD_=9, ASSIMI_=10, PHENO_=11, SB_=12

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
    n_struct, name, l2_struct, &
    hgrid, vgrid, cf, grib1, grib2, prefix, suffix)
    IMPLICIT NONE

    ! 0.1 InOut
    CLASS(t_jsb_n_struct),  INTENT(inout)       :: n_struct            ! memory structure
    CHARACTER(len=*),     INTENT(in)            :: name                ! name of variable
    CLASS(t_jsb_struct), TARGET, INTENT(inout)  :: l2_struct           ! l2_struct

    TYPE(t_jsb_grid),  POINTER, INTENT(in) :: hgrid     ! horizontal grid type
    TYPE(t_jsb_vgrid), POINTER, INTENT(in) :: vgrid     ! vertical grid type
    TYPE(t_cf),                 INTENT(in) :: cf        ! CF related metadata
    TYPE(t_grib1),              INTENT(in) :: grib1     ! GRIB1 related metadata
    TYPE(t_grib2),              INTENT(in) :: grib2     ! GRIB2 related metadata
    CHARACTER(len=*),           INTENT(in) :: prefix    ! Prefix for variable name
    CHARACTER(len=*),           INTENT(in) :: suffix    ! Suffix for variable name

    ! 0.2 Local
    INTEGER           :: pos_of_l2_struct
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':slm_add_l2_struct'

    l2_struct%name = TRIM(name)
    l2_struct%owner_proc_id = n_struct%owner_proc_id  ! A2L_=1, L2A_=2, SSE_=5, VEG_=7, RAD_=9, ASSIMI_=10, PHENO_=11, SB_=12

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
  !> add 2D variable to struct or l2_struct
  !!
  !!  t_jsb_struct:struct <- t_jsb_var_real2d:var
  !!  
  !! 
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE slm_add_struct_element_r2d( &
    struct, name, var, &
    hgrid, vgrid, cf, grib1, grib2, prefix, suffix, &
    initval_r)
    IMPLICIT NONE

    ! 0.1 InOut
    CLASS(t_jsb_struct),              INTENT(inout)   :: struct              ! struct & l2_struct
    CHARACTER(len=*),                 INTENT(in)      :: name                ! name of variable
    TYPE(t_jsb_var_real2d), TARGET,   INTENT(inout)   :: var                 ! reference to jsbach var -> newly created plain scalar REAL

    TYPE(t_jsb_grid),  POINTER, INTENT(in) :: hgrid     ! horizontal grid type
    TYPE(t_jsb_vgrid), POINTER, INTENT(in) :: vgrid     ! vertical grid type
    TYPE(t_cf),                 INTENT(in) :: cf        ! CF related metadata
    TYPE(t_grib1),              INTENT(in) :: grib1     ! GRIB1 related metadata
    TYPE(t_grib2),              INTENT(in) :: grib2     ! GRIB2 related metadata
    CHARACTER(len=*),           INTENT(in) :: prefix    ! Prefix for variable name
    CHARACTER(len=*),           INTENT(in) :: suffix    ! Suffix for variable name
    REAL(wp),         OPTIONAL, INTENT(in) :: initval_r ! value if var not available

    ! 0.2 Local
    INTEGER                       :: pos_of_var
    TYPE(t_jsb_model), POINTER    :: model
    INTEGER                       :: ice, ics, nproma, iblk   ! *slmtemp
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':slm_add_struct_element_real2d'

    model   => Get_model(1)        ! model => Get_model(model_id)   *slmtemp hardcoded model ID !!!
    ice     = model%options%ice    ! *slmtemp
    ics     = model%options%ics    ! *slmtemp
    nproma  = model%options%nproma ! *slmtemp
    iblk    = model%options%iblk   ! *slmtemp

    ALLOCATE(var%field(nproma,iblk,1,1,1))
    var%ptr => var%field(:,:,1,1,1)

    var%name = TRIM(name)
    var%owner_proc_id = struct%owner_proc_id    ! A2L_=1, L2A_=2, SSE_=5, VEG_=7, RAD_=9, ASSIMI_=10, PHENO_=11, SB_=12

    pos_of_var = struct%no_of_vars + 1
    IF (pos_of_var > struct%max_no_of_vars) THEN
      CALL finish(TRIM(routine), 'Maximum number of vars exceeded: '//TRIM(struct%varlist_name)//', '//TRIM(name))
    ELSE
      struct%vars(pos_of_var)%p => var
      struct%no_of_vars = pos_of_var
    END IF

    ! added for SLM compatibility (in JSBACH4 the allocation is done differently)   *slmtemp
    var%ptr(ics:ice,iblk) = initval_r

  END SUBROUTINE slm_add_struct_element_r2d


  !-----------------------------------------------------------------------------------------------------
  !> add 3D variable to struct or l2_struct
  !!
  !!  t_jsb_struct:struct <- t_jsb_var_real3d:var
  !!    
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE slm_add_struct_element_r3d( &
    struct, name, var, &
    hgrid, vgrid, cf, grib1, grib2, prefix, suffix, &
    initval_r)

    ! Use of process configurations (t_PROC_config)
    dsl4jsb_Use_config(ASSIMI_)
    dsl4jsb_Use_config(SSE_)

    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    CLASS(t_jsb_struct),              INTENT(inout)   :: struct              ! struct & l2_struct
    CHARACTER(len=*),                 INTENT(in)      :: name                ! name of variable
    TYPE(t_jsb_var_real3d), TARGET,   INTENT(inout)   :: var                 ! reference to jsbach var -> newly created plain scalar REAL

    TYPE(t_jsb_grid),  POINTER, INTENT(in) :: hgrid     ! horizontal grid type
    TYPE(t_jsb_vgrid), POINTER, INTENT(in) :: vgrid     ! vertical grid type
    TYPE(t_cf),                 INTENT(in) :: cf        ! CF related metadata
    TYPE(t_grib1),              INTENT(in) :: grib1     ! GRIB1 related metadata
    TYPE(t_grib2),              INTENT(in) :: grib2     ! GRIB2 related metadata
    CHARACTER(len=*),           INTENT(in) :: prefix    ! Prefix for variable name
    CHARACTER(len=*),           INTENT(in) :: suffix    ! Suffix for variable name
    REAL(wp),         OPTIONAL, INTENT(in) :: initval_r ! value if var not available
    ! ---------------------------
    ! 0.2 Local
    INTEGER                       :: pos_of_var
    TYPE(t_jsb_model), POINTER    :: model
    INTEGER                       :: ice, ics, nproma, iblk   ! *slmtemp
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':slm_add_struct_element_r3d'

    dsl4jsb_Def_config(ASSIMI_)
    dsl4jsb_Def_config(SSE_)

    model   => Get_model(1)        ! model => Get_model(model_id)   *slmtemp hardcoded model ID !!!
    ice     = model%options%ice    ! *slmtemp
    ics     = model%options%ics    ! *slmtemp
    nproma  = model%options%nproma ! *slmtemp
    iblk    = model%options%iblk   ! *slmtemp

    dsl4jsb_Get_config(ASSIMI_)
    dsl4jsb_Get_config(SSE_)


    ALLOCATE(var%field(nproma, vgrid%n_levels, iblk, 1, 1))
    var%ptr => var%field(:,:,:,1,1)
    
    var%name = TRIM(name)
    var%owner_proc_id = struct%owner_proc_id 

    pos_of_var = struct%no_of_vars + 1
    IF (pos_of_var > struct%max_no_of_vars) THEN
      CALL finish(TRIM(routine), 'Maximum number of vars exceeded: '//TRIM(struct%varlist_name)//', '//TRIM(name))
    ELSE
      struct%vars(pos_of_var)%p => var
      struct%no_of_vars = pos_of_var
    END IF

    ! added for SLM compatibility (in JSBACH4 the allocation is done differently)   *slmtemp
    var%ptr(ics:ice, 1:vgrid%n_levels, iblk) = initval_r

  END SUBROUTINE slm_add_struct_element_r3d



  !-----------------------------------------------------------------------------------------------------!
  ! add variable to process memory (mem%VAR)                                                            !
  !-----------------------------------------------------------------------------------------------------!

  !-----------------------------------------------------------------------------------------------------
  !> add 2D variable to process memory
  !!
  !!  t_jsb_memory:mem <- t_jsb_var_r2d:var
  !! 
  !! 
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE slm_add_var_list_element_r2d( &
    mem, name, var, &
    hgrid, vgrid, cf, grib1, grib2, prefix, suffix, &
    initval_r)
    IMPLICIT NONE

    ! 0.1 InOut
    CLASS(t_jsb_memory),              INTENT(inout)   :: mem                 ! memory structure
    CHARACTER(len=*),                 INTENT(in)      :: name                ! name of variable
    TYPE(t_jsb_var_real2d), TARGET,   INTENT(inout)   :: var                 ! reference to jsbach var -> newly created plain scalar REAL

    TYPE(t_jsb_grid),  POINTER, INTENT(in) :: hgrid     ! horizontal grid type
    TYPE(t_jsb_vgrid), POINTER, INTENT(in) :: vgrid     ! vertical grid type
    TYPE(t_cf),                 INTENT(in) :: cf        ! CF related metadata
    TYPE(t_grib1),              INTENT(in) :: grib1     ! GRIB1 related metadata
    TYPE(t_grib2),              INTENT(in) :: grib2     ! GRIB2 related metadata
    CHARACTER(len=*),           INTENT(in) :: prefix    ! Prefix for variable name
    CHARACTER(len=*),           INTENT(in) :: suffix    ! Suffix for variable name
    REAL(wp),         OPTIONAL, INTENT(in) :: initval_r ! value if var not available

    ! 0.2 Local
    INTEGER                       :: pos_of_var
    TYPE(t_jsb_model), POINTER    :: model
    INTEGER                       :: ice, ics, nproma, iblk   ! *slmtemp
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':slm_add_var_list_element_r2d'

    model   => Get_model(1)        ! model => Get_model(model_id)   *slmtemp hardcoded model ID !!!
    ice     = model%options%ice    ! *slmtemp
    ics     = model%options%ics    ! *slmtemp
    nproma  = model%options%nproma ! *slmtemp
    iblk    = model%options%iblk   ! *slmtemp

    ALLOCATE(var%field(nproma,iblk,1,1,1))
    var%ptr => var%field(:,:,1,1,1)

    var%name = TRIM(name)
    var%owner_proc_id = mem%owner_proc_id    ! A2L_=1, L2A_=2, SSE_=5, VEG_=7, RAD_=9, ASSIMI_=10, PHENO_=11, SB_=12

    pos_of_var = mem%no_of_vars + 1
    IF (pos_of_var > mem%max_no_of_vars) THEN
      CALL finish(TRIM(routine), 'Maximum number of vars exceeded: '//TRIM(mem%varlist_name)//', '//TRIM(name))
    ELSE
      mem%vars(pos_of_var)%p => var
      mem%no_of_vars = pos_of_var
    END IF

    ! added for SLM compatibility (in JSBACH4 the allocation is done differently)   *slmtemp
    var%ptr(ics:ice,iblk) = initval_r

  END SUBROUTINE slm_add_var_list_element_r2d


  !-----------------------------------------------------------------------------------------------------
  !> add 3D variable to process memory
  !!
  !!  t_jsb_memory:mem <- t_jsb_var_r3d:var
  !! 
  !!  Note ! the dimension of real3d variables [ncanopy / nsoil_water (nsoil_sb)] is not anymore exclusively determined by process_memory ! \n
  !!    with r1098 (ticket #342) a new CASE statement is introduced here, checking for some variables explicitly and set DIMENSION
  !! 
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE slm_add_var_list_element_r3d( &
    mem, name, var, &
    hgrid, vgrid, cf, grib1, grib2, prefix, suffix, &
    initval_r)

    ! Use of process configurations (t_PROC_config)
    dsl4jsb_Use_config(ASSIMI_)
    dsl4jsb_Use_config(SSE_)
    dsl4jsb_Use_config(RAD_)

    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    CLASS(t_jsb_memory),              INTENT(inout)   :: mem                 ! memory structure
    CHARACTER(len=*),                 INTENT(in)      :: name                ! name of variable
    TYPE(t_jsb_var_real3d), TARGET,   INTENT(inout)   :: var                 ! reference to jsbach var -> newly created plain scalar REAL

    TYPE(t_jsb_grid),  POINTER, INTENT(in) :: hgrid     ! horizontal grid type
    TYPE(t_jsb_vgrid), POINTER, INTENT(in) :: vgrid     ! vertical grid type
    TYPE(t_cf),                 INTENT(in) :: cf        ! CF related metadata
    TYPE(t_grib1),              INTENT(in) :: grib1     ! GRIB1 related metadata
    TYPE(t_grib2),              INTENT(in) :: grib2     ! GRIB2 related metadata
    CHARACTER(len=*),           INTENT(in) :: prefix    ! Prefix for variable name
    CHARACTER(len=*),           INTENT(in) :: suffix    ! Suffix for variable name
    REAL(wp),         OPTIONAL, INTENT(in) :: initval_r ! value if var not available
    ! ---------------------------
    ! 0.2 Local
    INTEGER                       :: pos_of_var
    TYPE(t_jsb_model), POINTER    :: model
    INTEGER                       :: ice, ics, nproma, iblk   ! *slmtemp
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':slm_add_var_list_element_r3d'

    dsl4jsb_Def_config(ASSIMI_)
    dsl4jsb_Def_config(SSE_)
    dsl4jsb_Def_config(RAD_)

    model   => Get_model(1)        ! model => Get_model(model_id)   *slmtemp hardcoded model ID !!!
    ice     = model%options%ice    ! *slmtemp
    ics     = model%options%ics    ! *slmtemp
    nproma  = model%options%nproma ! *slmtemp
    iblk    = model%options%iblk   ! *slmtemp

    dsl4jsb_Get_config(ASSIMI_)
    dsl4jsb_Get_config(SSE_)
    dsl4jsb_Get_config(RAD_)


    ALLOCATE(var%field(nproma, vgrid%n_levels, iblk, 1, 1))
    var%ptr => var%field(:,:,:,1,1)

    var%name = TRIM(name)
    var%owner_proc_id = mem%owner_proc_id

    pos_of_var = mem%no_of_vars + 1
    IF (pos_of_var > mem%max_no_of_vars) THEN
      CALL finish(TRIM(routine), 'Maximum number of vars exceeded: '//TRIM(mem%varlist_name)//', '//TRIM(name))
    ELSE
      mem%vars(pos_of_var)%p => var
      mem%no_of_vars = pos_of_var
    END IF

    ! added for SLM compatibility (in JSBACH4 the allocation is done differently)   *slmtemp
    var%ptr(ics:ice, 1:vgrid%n_levels, iblk) = initval_r

  END SUBROUTINE slm_add_var_list_element_r3d

END MODULE mo_slm_varlist
