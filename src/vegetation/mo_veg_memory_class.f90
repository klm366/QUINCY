!>
!! @brief Contains the memory class for the vegetation processes
!!
!! @author
!!  Jan Engel, MPI-BGC Jena (based on the JSBACH4 version from Rainer Schneck, MPI-M Hamburg)
!!
!! @par Revision History
!! First version by Jan Engel (20160830)
!!
MODULE mo_veg_memory_class

  USE mo_kind, ONLY: wp

  USE mo_jsb_memory_class,       ONLY: t_jsb_memory
  USE mo_jsb_var_class,          ONLY: t_jsb_var_real2d, t_jsb_var_real3d, &
                                       t_jsb_var_real1d, &                      ! @slmtemp
                                       t_jsb_struct, t_jsb_n_struct             ! @slmjsb4

  USE mo_slm_varlist,            ONLY: slm_add_var   ! replace by mo_jsb_varlist later    @slmtemp



  IMPLICIT NONE
  PRIVATE
  PUBLIC :: t_veg_memory, max_no_of_vars
  PUBLIC :: max_no_of_structs, max_no_of_n_structs  ! @slmjsb4
  PUBLIC :: t_slm_struct_veg_pool                   ! needed in mo_veg_dynamics for passing littefall to the routine test_routine()

  INTEGER, PARAMETER :: max_no_of_vars      = 60
  INTEGER, PARAMETER :: max_no_of_structs   = 10   ! @slmjsb4
  INTEGER, PARAMETER :: max_no_of_n_structs = 10   ! @slmjsb4



  !-----------------------------------------------------------------------------------------------------
  !> Type definition for elements in veg pool
  !! 
  !! elements in veg would be 2D (while they would be 3D in soil)
  !-----------------------------------------------------------------------------------------------------
  TYPE, EXTENDS(t_jsb_struct) :: t_slm_struct_veg_element
  ! @slmtemp 1D -> 2D
    TYPE(t_jsb_var_real1d)    :: carbon
    TYPE(t_jsb_var_real1d)    :: nitrogen
    ! TYPE(t_jsb_var_real1d)    :: phosphorus ...
  END TYPE t_slm_struct_veg_element


  !-----------------------------------------------------------------------------------------------------
  !> Type definition for veg pool with nested TYPEs
  !! 
  !-----------------------------------------------------------------------------------------------------
  TYPE, EXTENDS(t_jsb_n_struct) :: t_slm_struct_veg_pool
    TYPE(t_slm_struct_veg_element)  :: leaf
    TYPE(t_slm_struct_veg_element)  :: labile
    !TYPE(t_slm_struct_veg_element)  :: fine_root ...
  END TYPE t_slm_struct_veg_pool


  !-----------------------------------------------------------------------------------------------------
  !> Type definition for veg memory
  !! 
  !-----------------------------------------------------------------------------------------------------
  TYPE, EXTENDS(t_jsb_memory) :: t_veg_memory


  ! @slmtemp  --> 2D
    TYPE(t_jsb_var_real1d)      :: &              ! REAL(wp) 1D ; implemented for SLM use only, see mo_jsb_var_class 
                                 height    , &    !< plant height in m
                                 diameter  , &    !< plant diameter in m
                                 dens_ind  , &    !< density of individuals (#/m2) 
                                 lai              !< actual leaf area index
                                 ! target_cn_leaf ...



    TYPE(t_slm_struct_veg_pool) :: &
                                   pool       , &   !< nested pool TYPE t_jsb_n_struct -> extented to t_struct_pool
                                   growth     , &   !< fluxes: growth
                                   litterfall       !< fluxes: litterfall


  CONTAINS
    PROCEDURE :: Init => Init_veg_memory
  END TYPE t_veg_memory


  CHARACTER(len=*), PARAMETER :: modname = 'mo_veg_memory_class'




CONTAINS


  !-----------------------------------------------------------------------------------------------------
  !> initialize memory for the processes in: vegetation
  !!
  !!  
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE Init_veg_memory(mem, prefix, suffix, lct_ids, model_id, ics, ice, iblk)

!    USE mo_jsb_varlist,       ONLY: t_jsb_varlist, jsb_add_var
!    USE mo_jsb_io,            ONLY: grib_bits, t_cf, t_grib1, t_grib2, TSTEP_CONSTANT, tables
!    USE mo_jsb_grid_class,    ONLY: t_jsb_grid, t_jsb_vgrid
!    USE mo_jsb_grid,          ONLY: Get_grid, Get_vgrid

    CLASS(t_veg_memory), INTENT(inout), TARGET :: mem
    CHARACTER(len=*),    INTENT(in)    :: prefix
    CHARACTER(len=*),    INTENT(in)    :: suffix
    INTEGER,             INTENT(in)    :: lct_ids !(:)
    INTEGER,             INTENT(in)    :: model_id
    INTEGER,             INTENT(in)    :: ics   ! model%options -> for SLM compatibility
    INTEGER,             INTENT(in)    :: ice   ! model%options -> for SLM compatibility
    INTEGER,             INTENT(in)    :: iblk  ! model%options -> for SLM compatibility

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

    IF (.NOT. ASSOCIATED(mem%height%ptr)) THEN

      CALL slm_add_var(mem, 'height', mem%height, &
            ics, ice, iblk, &    ! model%options -> for SLM compatibility
            initval_r=0.0_wp) ! initval taken from JSBACH3

      CALL slm_add_var(mem, 'diameter', mem%diameter, &
            ics, ice, iblk, &    ! model%options -> for SLM compatibility
            initval_r=0.0_wp) ! initval taken from JSBACH3

      CALL slm_add_var(mem, 'dens_ind', mem%dens_ind, &
            ics, ice, iblk, &    ! model%options -> for SLM compatibility
            initval_r=0.0_wp) ! initval taken from JSBACH3

      CALL slm_add_var(mem, 'lai', mem%lai, &
            ics, ice, iblk, &    ! model%options -> for SLM compatibility
            initval_r=0.0_wp) ! initval taken from JSBACH3


      ! nested TYPE (t_jsb_n_struct)
      CALL slm_add_var(mem, 'slm_veg_pool', mem%pool)
      CALL slm_add_var(mem, 'slm_veg_growth', mem%growth)
      CALL slm_add_var(mem, 'slm_veg_litterfall', mem%litterfall)
      ! level 2 TYPE (l2_struct) nested within the above mem%pool
      CALL slm_add_var(mem%pool, 'slm_veg_pool_leaf', mem%pool%leaf)                          ! pool
      CALL slm_add_var(mem%pool, 'slm_veg_pool_labile', mem%pool%labile)
      CALL slm_add_var(mem%growth, 'slm_veg_flux_growth_leaf', mem%growth%leaf)               ! flux growth
      CALL slm_add_var(mem%growth, 'slm_veg_flux_growth_labile', mem%growth%labile)
      CALL slm_add_var(mem%litterfall, 'slm_veg_flux_litterfall_leaf', mem%litterfall%leaf)   ! flux litterfall
      CALL slm_add_var(mem%litterfall, 'slm_veg_flux_litterfall_labile', mem%litterfall%labile)
      ! add variables to the l2 structs for each nested TYPE
      CALL slm_add_var(mem%pool%leaf, 'pool_leaf_carbon', mem%pool%leaf%carbon, &   
            ics, ice, iblk, &   ! model%options -> for SLM compatibility
            initval_r=0.0_wp)   ! initval taken from JSBACH3
      CALL slm_add_var(mem%pool%leaf, 'pool_leaf_nitrogen', mem%pool%leaf%nitrogen, &   
            ics, ice, iblk, &   ! model%options -> for SLM compatibility
            initval_r=0.0_wp)   ! initval taken from JSBACH3
      CALL slm_add_var(mem%pool%labile, 'pool_labile_carbon', mem%pool%labile%carbon, &   
            ics, ice, iblk, &   ! model%options -> for SLM compatibility
            initval_r=0.0_wp)   ! initval taken from JSBACH3
      CALL slm_add_var(mem%pool%labile, 'pool_labile_nitrogen', mem%pool%labile%nitrogen, &   
            ics, ice, iblk, &   ! model%options -> for SLM compatibility
            initval_r=0.0_wp)   ! initval taken from JSBACH3

      CALL slm_add_var(mem%growth%leaf, 'flux_growth_leaf_carbon', mem%growth%leaf%carbon, &   
            ics, ice, iblk, &   ! model%options -> for SLM compatibility
            initval_r=0.0_wp)   ! initval taken from JSBACH3
      CALL slm_add_var(mem%growth%leaf, 'flux_growth_leaf_nitrogen', mem%growth%leaf%nitrogen, &   
            ics, ice, iblk, &   ! model%options -> for SLM compatibility
            initval_r=0.0_wp)   ! initval taken from JSBACH3
      CALL slm_add_var(mem%growth%labile, 'flux_growth_labile_carbon', mem%growth%labile%carbon, &   
            ics, ice, iblk, &   ! model%options -> for SLM compatibility
            initval_r=0.0_wp)   ! initval taken from JSBACH3
      CALL slm_add_var(mem%growth%labile, 'flux_growth_labile_nitrogen', mem%growth%labile%nitrogen, &   
            ics, ice, iblk, &   ! model%options -> for SLM compatibility
            initval_r=0.0_wp)   ! initval taken from JSBACH3

      CALL slm_add_var(mem%litterfall%leaf, 'flux_litterfall_leaf_carbon', mem%litterfall%leaf%carbon, &   
            ics, ice, iblk, &   ! model%options -> for SLM compatibility
            initval_r=0.0_wp)   ! initval taken from JSBACH3
      CALL slm_add_var(mem%litterfall%leaf, 'flux_litterfall_leaf_nitrogen', mem%litterfall%leaf%nitrogen, &   
            ics, ice, iblk, &   ! model%options -> for SLM compatibility
            initval_r=0.0_wp)   ! initval taken from JSBACH3
      CALL slm_add_var(mem%litterfall%labile, 'flux_litterfall_labile_carbon', mem%litterfall%labile%carbon, &   
            ics, ice, iblk, &   ! model%options -> for SLM compatibility
            initval_r=0.0_wp)   ! initval taken from JSBACH3
      CALL slm_add_var(mem%litterfall%labile, 'flux_litterfall_labile_nitrogen', mem%litterfall%labile%nitrogen, &   
            ics, ice, iblk, &   ! model%options -> for SLM compatibility
            initval_r=0.0_wp)   ! initval taken from JSBACH3

    END IF

  END SUBROUTINE Init_veg_memory

END MODULE mo_veg_memory_class
