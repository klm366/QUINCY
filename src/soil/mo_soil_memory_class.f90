!>
!! @brief Contains the memory class for the soil processes
!!
!! @author
!!  Jan Engel, MPI-BGC Jena (based on the JSBACH4 version from Rainer Schneck, MPI-M Hamburg)
!!
!! @par Revision History
!! First version by Jan Engel (20160913)
!!
MODULE mo_soil_memory_class

  USE mo_kind, ONLY: wp

  USE mo_jsb_memory_class,       ONLY: t_jsb_memory
  USE mo_jsb_var_class,          ONLY: t_jsb_var_real2d, t_jsb_var_real3d, &
                                       t_jsb_struct, t_jsb_n_struct             ! *slmjsb4

  USE mo_slm_varlist,            ONLY: slm_add_var   ! replace by mo_jsb_varlist later    *slmtemp



  IMPLICIT NONE
  PRIVATE
  PUBLIC :: t_soil_memory, max_no_of_vars
  PUBLIC :: max_no_of_structs, max_no_of_n_structs  ! *slmjsb4

  INTEGER, PARAMETER :: max_no_of_vars      = 40   !> TYPE t_jsb_var
  INTEGER, PARAMETER :: max_no_of_structs   = 20   !> TYPE t_jsb_struct     *slmjsb4
  INTEGER, PARAMETER :: max_no_of_n_structs = 20   !> TYPE t_jsb_n_struct   *slmjsb4



  !-----------------------------------------------------------------------------------------------------
  !> Type definition for elements in soil pool
  !! 
  !! elements in soil are 3D and would be 2D in veg
  !! 
  !-----------------------------------------------------------------------------------------------------
  TYPE, EXTENDS(t_jsb_struct) :: t_slm_struct_soil_element
    TYPE(t_jsb_var_real3d)    :: carbon       !< total carbon concentration of pool in mol / m2
    TYPE(t_jsb_var_real3d)    :: nitrogen     !< totol nitrogen concentration of pool in mol / m2
    TYPE(t_jsb_var_real3d)    :: phosphorus   !< phosphorus concentration of pool in mol / m2
    TYPE(t_jsb_var_real3d)    :: carbon13     !< 13carbon concentration of pool in mol / m2
    TYPE(t_jsb_var_real3d)    :: carbon14     !< 14carbon concentration of pool in pmol / m2
    TYPE(t_jsb_var_real3d)    :: nitrogen15   !< 15nitrogen concentration of pool in mol / m2
  CONTAINS
  END TYPE t_slm_struct_soil_element


  !-----------------------------------------------------------------------------------------------------
  !> Type definition for soil pool with nested TYPEs
  !!
  !! used for: pool ...
  !!  
  !-----------------------------------------------------------------------------------------------------
  TYPE, EXTENDS(t_jsb_n_struct) :: t_slm_struct_soil_pool
    TYPE(t_slm_struct_soil_element)  :: test1           !< elements of 
    TYPE(t_slm_struct_soil_element)  :: test2           !< elements of 
    TYPE(t_slm_struct_soil_element)  :: test3           !< elements of 
    TYPE(t_slm_struct_soil_element)  :: test4           !< elements of 
  END TYPE t_slm_struct_soil_pool


  !-----------------------------------------------------------------------------------------------------
  !> Type definition for soil memory
  !!
  !! @par includes: \n
  !!    soil variables soil\%VAR \n
  !!    pool  \n
  !!    ? flux \n
  !!    ? flux \n
  !!    moving averagages of variables VAR_mavg   \n
  !!    daytime averages of variables VAR_daytime \n
  !! 
  !-----------------------------------------------------------------------------------------------------
  TYPE, EXTENDS(t_jsb_memory) :: t_soil_memory


    TYPE(t_jsb_var_real3d)          :: &
                                    test          !< test


!    TYPE(t_slm_struct_soil_element)  :: &
!                                    allometry_labile_mavg   , &  !< moving average of labile pools neeeded for allometry calculations
!                                    establishment           , &  !< flux between seed bed and reserve (mol m-2 timestep-1)
!                                    reserve_use                  !< net flux between reserve and labile pool (mol m-2 timestep-1)

!    TYPE(t_slm_struct_soil_pool)     :: &             ! nested pool TYPE t_jsb_n_struct
!                                    pool       , &   !< pools of the plant
!                                    growth     , &   !< fluxes: growth
!                                    litterfall       !< fluxes: litterfall


  CONTAINS
    PROCEDURE :: Init => Init_soil_memory
  END TYPE t_soil_memory


  CHARACTER(len=*), PARAMETER :: modname = 'mo_soil_memory_class'


CONTAINS


  !-----------------------------------------------------------------------------------------------------
  !> initialize memory for the processes in: soil
  !!
  !!  
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE Init_soil_memory(mem, prefix, suffix, lct_ids, model_id, ics, ice, iblk)

!    USE mo_jsb_varlist,       ONLY: t_jsb_varlist, jsb_add_var
!    USE mo_jsb_io,            ONLY: grib_bits, t_cf, t_grib1, t_grib2, TSTEP_CONSTANT, tables
!    USE mo_jsb_grid_class,    ONLY: t_jsb_grid, t_jsb_vgrid
!    USE mo_jsb_grid,          ONLY: Get_grid, Get_vgrid

    CLASS(t_soil_memory), INTENT(inout), TARGET :: mem    !> soil memory
    CHARACTER(len=*),    INTENT(in)    :: prefix          !> ?
    CHARACTER(len=*),    INTENT(in)    :: suffix          !> ?
    INTEGER,             INTENT(in)    :: lct_ids!(:)     !> in jsbach4: vector of land cover type IDs
    INTEGER,             INTENT(in)    :: model_id        !> model ID
    INTEGER,             INTENT(in)    :: ics             !> model%options -> for SLM compatibility
    INTEGER,             INTENT(in)    :: ice             !> model%options -> for SLM compatibility
    INTEGER,             INTENT(in)    :: iblk            !> model%options -> for SLM compatibility

!    TYPE(t_jsb_model), POINTER :: model
!    TYPE(t_jsb_grid),  POINTER :: hgrid                        ! Horizontal grid
!    TYPE(t_jsb_vgrid), POINTER :: surface, canopy_layer        ! Vertical grid
!    INTEGER :: table

!    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':Init_assimi_memory'

!    model => Get_model(model_id)

!    table = tables(1)

!    hgrid   => Get_grid(mem%grid_id)
!    surface => Get_vgrid('surface')
!    canopy_layer => Get_vgrid('canopy_layer')

    IF (.NOT. ASSOCIATED(mem%test%ptr)) THEN

      !> add variables to soil memory soil\%VAR
      !! 
      !! 
      CALL slm_add_var(mem, 'test', mem%test, &
            ics, ice, iblk, &    ! model%options -> for SLM compatibility
            initval_r=0.0_wp)



    END IF

  END SUBROUTINE Init_soil_memory

END MODULE mo_soil_memory_class
