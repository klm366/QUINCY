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
                                       t_jsb_var_real1d, &                      ! *slmtemp
                                       t_jsb_struct, t_jsb_n_struct             ! *slmjsb4

  USE mo_slm_varlist,            ONLY: slm_add_var   ! replace by mo_jsb_varlist later    *slmtemp



  IMPLICIT NONE
  PRIVATE
  PUBLIC :: t_veg_memory, max_no_of_vars
  PUBLIC :: max_no_of_structs, max_no_of_n_structs  ! *slmjsb4
  PUBLIC :: t_slm_struct_veg_pool                   ! needed in mo_veg_dynamics for passing littefall to the routine test_routine()

  INTEGER, PARAMETER :: max_no_of_vars      = 80   !> TYPE t_jsb_var
  INTEGER, PARAMETER :: max_no_of_structs   = 20   !> TYPE t_jsb_struct     *slmjsb4
  INTEGER, PARAMETER :: max_no_of_n_structs = 20   !> TYPE t_jsb_n_struct   *slmjsb4



  !-----------------------------------------------------------------------------------------------------
  !> Type definition for elements in veg pool
  !! 
  !! elements in veg would be 2D (while they would be 3D in soil)
  !! 
  !! *slmtemp 1D -> 2D
  !-----------------------------------------------------------------------------------------------------
  TYPE, EXTENDS(t_jsb_struct) :: t_slm_struct_veg_element
    TYPE(t_jsb_var_real1d)    :: carbon       !< total carbon concentration of pool in mol / m2
    TYPE(t_jsb_var_real1d)    :: nitrogen     !< totol nitrogen concentration of pool in mol / m2
    TYPE(t_jsb_var_real1d)    :: phosphorus   !< phosphorus concentration of pool in mol / m2
    TYPE(t_jsb_var_real1d)    :: carbon13     !< 13carbon concentration of pool in mol / m2
    TYPE(t_jsb_var_real1d)    :: carbon14     !< 14carbon concentration of pool in pmol / m2
    TYPE(t_jsb_var_real1d)    :: nitrogen15   !< 15nitrogen concentration of pool in mol / m2
  CONTAINS
    PROCEDURE :: Add_struct                       ! Add two structs of type t_jsb_struct
    PROCEDURE :: Add_struct_and_scalar            ! Add struct of type t_jsb_struct with REAL scalar
    GENERIC   :: OPERATOR(+) => Add_struct, Add_struct_and_scalar

  END TYPE t_slm_struct_veg_element


  !-----------------------------------------------------------------------------------------------------
  !> Type definition for veg pool with nested TYPEs
  !!
  !! used for: pool, and fluxes: growth and litterfall
  !!  
  !-----------------------------------------------------------------------------------------------------
  TYPE, EXTENDS(t_jsb_n_struct) :: t_slm_struct_veg_pool
    TYPE(t_slm_struct_veg_element)  :: leaf             !< elements of leaf
    TYPE(t_slm_struct_veg_element)  :: fine_root        !< elements of fine_root
    TYPE(t_slm_struct_veg_element)  :: coarse_root      !< elements of coarse_root
    TYPE(t_slm_struct_veg_element)  :: sap_wood         !< elements of sap_wood
    TYPE(t_slm_struct_veg_element)  :: heart_wood       !< elements of heart_wood
    TYPE(t_slm_struct_veg_element)  :: labile           !< elements of labile
    TYPE(t_slm_struct_veg_element)  :: reserve          !< elements of reserve
    TYPE(t_slm_struct_veg_element)  :: fruit            !< elements of fruit
    TYPE(t_slm_struct_veg_element)  :: seed_bed         !< elements of seed_bed
  END TYPE t_slm_struct_veg_pool


  !-----------------------------------------------------------------------------------------------------
  !> Type definition for veg memory
  !!
  !! @par includes: \n
  !!    veg variables veg\%VAR \n
  !!    pool  \n
  !!    growth flux \n
  !!    litterfall flux \n
  !!    moving averagages of variables VAR_mavg   \n
  !!    daytime averages of variables VAR_daytime \n
  !! 
  !-----------------------------------------------------------------------------------------------------
  TYPE, EXTENDS(t_jsb_memory) :: t_veg_memory


  ! *slmtemp  --> 2D
    TYPE(t_jsb_var_real1d)          :: &              ! REAL(wp) 1D ; implemented for SLM use only, see mo_jsb_var_class 
                                    height      , &   !< plant height in m
                                    diameter    , &   !< plant diameter in m
                                    dens_ind    , &   !< density of individuals (#/m2) 
                                    lai               !< actual leaf area index
                                    ! target_cn_leaf ...

    TYPE(t_slm_struct_veg_element)  :: &
                                    allometry_labile_mavg   , &  !< moving average of labile pools neeeded for allometry calculations
                                    establishment           , &  !< flux between seed bed and reserve (mol m-2 timestep-1)
                                    reserve_use                  !< net flux between reserve and labile pool (mol m-2 timestep-1)

    TYPE(t_slm_struct_veg_pool)     :: &             ! nested pool TYPE t_jsb_n_struct
                                    pool       , &   !< pools of the plant
                                    growth     , &   !< fluxes: growth
                                    litterfall       !< fluxes: litterfall


  CONTAINS
    PROCEDURE :: Init => Init_veg_memory
  END TYPE t_veg_memory


  CHARACTER(len=*), PARAMETER :: modname = 'mo_veg_memory_class'


CONTAINS

  !-----------------------------------------------------------------------------------------------------
  !> operator overload +
  !!
  !! t_slm_struct_veg_element + t_slm_struct_veg_element
  !-----------------------------------------------------------------------------------------------------
  FUNCTION Add_struct(this,struct2) RESULT(struct3)

    CLASS(t_slm_struct_veg_element), INTENT(in) :: this, struct2
    TYPE(t_slm_struct_veg_element)              :: struct3

!    struct3%vars(1) = this%vars(1) + struct2%vars(1)
!    struct3%vars(2) = this%vars(2) + struct2%vars(2)
    struct3%carbon    = this%carbon    + struct2%carbon
    struct3%nitrogen  = this%nitrogen  + struct2%nitrogen
  END FUNCTION Add_struct

  !-----------------------------------------------------------------------------------------------------
  !> operator overload +
  !!
  !! t_slm_struct_veg_element + scalar
  !-----------------------------------------------------------------------------------------------------
  FUNCTION Add_struct_and_scalar(this,scalar) RESULT(struct3)

    CLASS(t_slm_struct_veg_element), INTENT(in) :: this
    REAL(wp),                        INTENT(in) :: scalar
    TYPE(t_slm_struct_veg_element)              :: struct3

!    struct3%vars(1) = this%vars(1) + scalar
!    struct3%vars(2) = this%vars(2) + scalar
    struct3%carbon    = this%carbon    + scalar
    struct3%nitrogen  = this%nitrogen  + scalar
  END FUNCTION Add_struct_and_scalar

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

      !> add variables to veg memory veg\%VAR
      !! 
      !! slm_add_var -> slm_add_var_list_element_r1d    *slmtemp 1D -> 2D
      !! 
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


      !> simple TYPE (t_jsb_struct)
      !! 
      !! slm_add_var -> slm_add_struct()
      !! 
      CALL slm_add_var(mem, 'veg_allometry_labile_mavg',  mem%allometry_labile_mavg)
      CALL slm_add_var(mem, 'veg_establishment',          mem%establishment)
      CALL slm_add_var(mem, 'veg_reserve_use',            mem%reserve_use)
      !> nested TYPE (t_jsb_n_struct)
      !! 
      !! slm_add_var -> slm_add_n_struct()
      !! 
      CALL slm_add_var(mem, 'veg_pool',       mem%pool)
      CALL slm_add_var(mem, 'veg_growth',     mem%growth)
      CALL slm_add_var(mem, 'veg_litterfall', mem%litterfall)
      !> level 2 TYPE (l2_struct) nested within the above nested TYPE structures (t_jsb_n_struct)
      !! 
      !! slm_add_var -> slm_add_l2_struct()
      !! 
      !! veg%pool \n
      !! veg%growth : fluxes  \n
      !! veg%litterfall : fluxes  \n
      CALL slm_add_var(mem%pool, 'veg_pool_leaf',         mem%pool%leaf)       
      CALL slm_add_var(mem%pool, 'veg_pool_fine_root',    mem%pool%fine_root)
      CALL slm_add_var(mem%pool, 'veg_pool_coarse_root',  mem%pool%coarse_root)
      CALL slm_add_var(mem%pool, 'veg_pool_sap_wood',     mem%pool%sap_wood)
      CALL slm_add_var(mem%pool, 'veg_pool_heart_wood',   mem%pool%heart_wood)
      CALL slm_add_var(mem%pool, 'veg_pool_labile',       mem%pool%labile)
      CALL slm_add_var(mem%pool, 'veg_pool_reserve',      mem%pool%reserve)
      CALL slm_add_var(mem%pool, 'veg_pool_fruit',        mem%pool%fruit)
      CALL slm_add_var(mem%pool, 'veg_pool_seed_bed',     mem%pool%seed_bed)
      CALL slm_add_var(mem%growth, 'veg_growth_leaf',         mem%growth%leaf)       
      CALL slm_add_var(mem%growth, 'veg_growth_fine_root',    mem%growth%fine_root)
      CALL slm_add_var(mem%growth, 'veg_growth_coarse_root',  mem%growth%coarse_root)
      CALL slm_add_var(mem%growth, 'veg_growth_sap_wood',     mem%growth%sap_wood)
      CALL slm_add_var(mem%growth, 'veg_growth_heart_wood',   mem%growth%heart_wood)
      CALL slm_add_var(mem%growth, 'veg_growth_labile',       mem%growth%labile)
      CALL slm_add_var(mem%growth, 'veg_growth_reserve',      mem%growth%reserve)
      CALL slm_add_var(mem%growth, 'veg_growth_fruit',        mem%growth%fruit)
      CALL slm_add_var(mem%growth, 'veg_growth_seed_bed',     mem%growth%seed_bed)
      CALL slm_add_var(mem%litterfall, 'veg_litterfall_leaf',         mem%litterfall%leaf)       
      CALL slm_add_var(mem%litterfall, 'veg_litterfall_fine_root',    mem%litterfall%fine_root)
      CALL slm_add_var(mem%litterfall, 'veg_litterfall_coarse_root',  mem%litterfall%coarse_root)
      CALL slm_add_var(mem%litterfall, 'veg_litterfall_sap_wood',     mem%litterfall%sap_wood)
      CALL slm_add_var(mem%litterfall, 'veg_litterfall_heart_wood',   mem%litterfall%heart_wood)
      CALL slm_add_var(mem%litterfall, 'veg_litterfall_labile',       mem%litterfall%labile)
      CALL slm_add_var(mem%litterfall, 'veg_litterfall_reserve',      mem%litterfall%reserve)
      CALL slm_add_var(mem%litterfall, 'veg_litterfall_fruit',        mem%litterfall%fruit)
      CALL slm_add_var(mem%litterfall, 'veg_litterfall_seed_bed',     mem%litterfall%seed_bed)
      !> add variables to the structs (veg\%t_jsb_struct\%VAR)
      !! 


      !> add variables to the l2 structs (veg\%t_jsb_n_struct\%t_jsb_struct\%VAR)
      !! 
      !! slm_add_var -> slm_add_struct_element_r1d    *slmtemp 1D -> 2D
      !! 
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
