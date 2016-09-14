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
  USE mo_jsb_var_class,          ONLY: t_jsb_var_real2d, t_jsb_var_real3d         , &
                                       t_jsb_var_integer2d, t_jsb_var_logical2d   , &                 ! *slmjsb4
                                       t_jsb_struct, t_jsb_n_struct               , &                 ! *slmjsb4
                                       t_jsb_var_real1d, t_jsb_var_integer1d, t_jsb_var_logical1d     ! *slmtemp

  USE mo_slm_varlist,            ONLY: slm_add_var   ! replace by mo_jsb_varlist later    *slmtemp



  IMPLICIT NONE
  PRIVATE
  PUBLIC :: t_veg_memory, max_no_of_vars
  PUBLIC :: max_no_of_structs, max_no_of_n_structs  ! *slmjsb4
  PUBLIC :: t_slm_struct_veg_pool                   ! needed in mo_veg_dynamics for passing littefall to the routine test_routine()

  INTEGER, PARAMETER :: max_no_of_vars      = 150  !> TYPE t_jsb_var
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
!    PROCEDURE :: Add_struct                       ! Add two structs of type t_jsb_struct
!    PROCEDURE :: Add_struct_and_scalar            ! Add struct of type t_jsb_struct with REAL scalar
!    GENERIC   :: OPERATOR(+) => Add_struct, Add_struct_and_scalar
! 
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


    TYPE(t_jsb_var_integer1d)      :: pftid             !< ID of the plant functional type; namelist parameter -> lctlib
    TYPE(t_jsb_var_logical1d)      :: growing_season    !< logical to identify whether the plant is in the growing season or not

  ! *slmtemp  --> 2D
    !> variables
    !! 
!    TYPE(t_jsb_var_real1d)          :: &   
!                                    height               , &             !< plant height in m
!                                    diameter             , &             !< plant diameter in m
!                                    dens_ind             , &             !< density of individuals (#/m2) 
!                                    lai                  , &             !< actual leaf area index
!                                    target_cn_leaf       , &             !< target foliar C:N ratio
!                                    target_np_leaf       , &             !< target foliar N:P ratio 
!                                    mean_leaf_age        , &             !< average foliage age in days 
!                                    gc                   , &             !< some kind of stomatal conductance ?
!                                    f_n_demand           , &             !< relative plant demand for N; see update_vegsoil_bgc
!                                    f_p_demand           , &             !< relative plant demand for P; see update_vegsoil_bgc
!                                    beta_air             , &             !< scaling factor to account for air humidity effects on conductance
!                                    beta_soil_ps         , &             !< scaling factor to account for soil moisture constraints on conductance
!                                    beta_soil_gs         , &             !< scaling factor to account for soil moisture constraints on conductance
!                                    dphi                 , &             !< plant-soil water potential gradient (MPa)
!                                    growth_req_n         , &             !< moles N required for a unit of C growth under current allocation fractions 
!                                    growth_req_p         , &             !< moles P required for a unit of N growth under current allocation fractions 
!                                    gdd                  , &             !< growing degree days above 5°C
!                                    nd_dormance          , &             !< number of days since last growing season 
!                                    t_jmax_opt           , &             !< current optimum temperature for eletron transport of photosynthesis
!                                    boc_an               , &             !< net photosynthesis at bottom of canopy (lowest canopy layer existing)
!                                    leaf_cn_direction    , &             !< indicates direction of change in leaf CN ratio to maximise NPP
!                                    apar                                 !< ..



    !> structures
    !! 
    TYPE(t_slm_struct_veg_element)  :: &
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

!  !-----------------------------------------------------------------------------------------------------
!  !> operator overload +
!  !!
!  !! t_slm_struct_veg_element + t_slm_struct_veg_element
!  !-----------------------------------------------------------------------------------------------------
!  FUNCTION Add_struct(this,struct2) RESULT(struct3)

!    CLASS(t_slm_struct_veg_element), INTENT(in) :: this, struct2
!    TYPE(t_slm_struct_veg_element)              :: struct3

!!    struct3%vars(1) = this%vars(1) + struct2%vars(1)
!!    struct3%vars(2) = this%vars(2) + struct2%vars(2)
!    ALLOCATE(struct3%vars(20))
!    struct3%carbon    = this%carbon    + struct2%carbon
!    struct3%nitrogen  = this%nitrogen  + struct2%nitrogen
!  END FUNCTION Add_struct

!  !-----------------------------------------------------------------------------------------------------
!  !> operator overload +
!  !!
!  !! t_slm_struct_veg_element + scalar
!  !-----------------------------------------------------------------------------------------------------
!  FUNCTION Add_struct_and_scalar(this,scalar) RESULT(struct3)

!    CLASS(t_slm_struct_veg_element), INTENT(in) :: this
!    REAL(wp),                        INTENT(in) :: scalar
!    TYPE(t_slm_struct_veg_element)              :: struct3

!!    struct3%vars(1) = this%vars(1) + scalar
!!    struct3%vars(2) = this%vars(2) + scalar
!    ALLOCATE(struct3%vars(20))
!    struct3%carbon    = this%carbon    + scalar
!    struct3%nitrogen  = this%nitrogen  + scalar
!  END FUNCTION Add_struct_and_scalar

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

    IF (.NOT. ASSOCIATED(mem%pftid%ptr)) THEN

      !> add variables to veg memory veg\%VAR
      !! 
      !! slm_add_var -> slm_add_var_list_element_r1d    *slmtemp 1D -> 2D
      !!                slm_add_var_list_element_i1d
      !!                slm_add_var_list_element_l1d
      !!  
      CALL slm_add_var(mem, 'veg_pftid', mem%pftid, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0)          ! initvalue 
      CALL slm_add_var(mem, 'veg_growing_season', mem%growing_season, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=.FALSE.)    ! initvalue 

      !> simple TYPE (t_jsb_struct)
      !! 
      !! slm_add_var -> slm_add_struct()
      !! 
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
      !! e.g. establishment
      !! slm_add_var -> slm_add_struct_element_r1d    *slmtemp 1D -> 2D
      !! 
      CALL slm_add_var(mem%establishment, 'veg_establishment_carbon', mem%establishment%carbon, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%establishment, 'veg_establishment_nitrogen', mem%establishment%nitrogen, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%establishment, 'veg_establishment_phosphorus', mem%establishment%phosphorus, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%establishment, 'veg_establishment_carbon13', mem%establishment%carbon13, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%establishment, 'veg_establishment_carbon14', mem%establishment%carbon14, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%establishment, 'veg_establishment_nitrogen15', mem%establishment%nitrogen15, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
 
      CALL slm_add_var(mem%reserve_use, 'veg_reserve_use_carbon', mem%reserve_use%carbon, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%reserve_use, 'veg_reserve_use_nitrogen', mem%reserve_use%nitrogen, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%reserve_use, 'veg_reserve_use_phosphorus', mem%reserve_use%phosphorus, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%reserve_use, 'veg_reserve_use_carbon13', mem%reserve_use%carbon13, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%reserve_use, 'veg_reserve_use_carbon14', mem%reserve_use%carbon14, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%reserve_use, 'veg_reserve_use_nitrogen15', mem%reserve_use%nitrogen15, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      !> add variables to the l2 structs (pool and flux variables, veg\%t_jsb_n_struct\%t_jsb_struct\%VAR)
      !! 
      !! slm_add_var -> slm_add_struct_element_r1d    *slmtemp 1D -> 2D
      !! 
      CALL slm_add_var(mem%pool%leaf, 'veg_pool_leaf_carbon', mem%pool%leaf%carbon, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%pool%leaf, 'veg_pool_leaf_nitrogen', mem%pool%leaf%nitrogen, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%pool%leaf, 'veg_pool_leaf_phosphorus', mem%pool%leaf%phosphorus, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%pool%leaf, 'veg_pool_leaf_carbon13', mem%pool%leaf%carbon13, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%pool%leaf, 'veg_pool_leaf_carbon14', mem%pool%leaf%carbon14, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%pool%leaf, 'veg_pool_leaf_nitrogen15', mem%pool%leaf%nitrogen15, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
 
      CALL slm_add_var(mem%pool%fine_root, 'veg_pool_fine_root_carbon', mem%pool%fine_root%carbon, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%pool%fine_root, 'veg_pool_fine_root_nitrogen', mem%pool%fine_root%nitrogen, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%pool%fine_root, 'veg_pool_fine_root_phosphorus', mem%pool%fine_root%phosphorus, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%pool%fine_root, 'veg_pool_fine_root_carbon13', mem%pool%fine_root%carbon13, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%pool%fine_root, 'veg_pool_fine_root_carbon14', mem%pool%fine_root%carbon14, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%pool%fine_root, 'veg_pool_fine_root_nitrogen15', mem%pool%fine_root%nitrogen15, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
 
      CALL slm_add_var(mem%pool%coarse_root, 'veg_pool_coarse_root_carbon', mem%pool%coarse_root%carbon, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%pool%coarse_root, 'veg_pool_coarse_root_nitrogen', mem%pool%coarse_root%nitrogen, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%pool%coarse_root, 'veg_pool_coarse_root_phosphorus', mem%pool%coarse_root%phosphorus, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%pool%coarse_root, 'veg_pool_coarse_root_carbon13', mem%pool%coarse_root%carbon13, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%pool%coarse_root, 'veg_pool_coarse_root_carbon14', mem%pool%coarse_root%carbon14, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%pool%coarse_root, 'veg_pool_coarse_root_nitrogen15', mem%pool%coarse_root%nitrogen15, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
 
      CALL slm_add_var(mem%pool%sap_wood, 'veg_pool_sap_wood_carbon', mem%pool%sap_wood%carbon, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%pool%sap_wood, 'veg_pool_sap_wood_nitrogen', mem%pool%sap_wood%nitrogen, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%pool%sap_wood, 'veg_pool_sap_wood_phosphorus', mem%pool%sap_wood%phosphorus, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%pool%sap_wood, 'veg_pool_sap_wood_carbon13', mem%pool%sap_wood%carbon13, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%pool%sap_wood, 'veg_pool_sap_wood_carbon14', mem%pool%sap_wood%carbon14, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%pool%sap_wood, 'veg_pool_sap_wood_nitrogen15', mem%pool%sap_wood%nitrogen15, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
 
      CALL slm_add_var(mem%pool%heart_wood, 'veg_pool_heart_wood_carbon', mem%pool%heart_wood%carbon, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%pool%heart_wood, 'veg_pool_heart_wood_nitrogen', mem%pool%heart_wood%nitrogen, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%pool%heart_wood, 'veg_pool_heart_wood_phosphorus', mem%pool%heart_wood%phosphorus, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%pool%heart_wood, 'veg_pool_heart_wood_carbon13', mem%pool%heart_wood%carbon13, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%pool%heart_wood, 'veg_pool_heart_wood_carbon14', mem%pool%heart_wood%carbon14, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%pool%heart_wood, 'veg_pool_heart_wood_nitrogen15', mem%pool%heart_wood%nitrogen15, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
 
      CALL slm_add_var(mem%pool%labile, 'veg_pool_labile_carbon', mem%pool%labile%carbon, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%pool%labile, 'veg_pool_labile_nitrogen', mem%pool%labile%nitrogen, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%pool%labile, 'veg_pool_labile_phosphorus', mem%pool%labile%phosphorus, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%pool%labile, 'veg_pool_labile_carbon13', mem%pool%labile%carbon13, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%pool%labile, 'veg_pool_labile_carbon14', mem%pool%labile%carbon14, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%pool%labile, 'veg_pool_labile_nitrogen15', mem%pool%labile%nitrogen15, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
 
      CALL slm_add_var(mem%pool%reserve, 'veg_pool_reserve_carbon', mem%pool%reserve%carbon, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%pool%reserve, 'veg_pool_reserve_nitrogen', mem%pool%reserve%nitrogen, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%pool%reserve, 'veg_pool_reserve_phosphorus', mem%pool%reserve%phosphorus, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%pool%reserve, 'veg_pool_reserve_carbon13', mem%pool%reserve%carbon13, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%pool%reserve, 'veg_pool_reserve_carbon14', mem%pool%reserve%carbon14, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%pool%reserve, 'veg_pool_reserve_nitrogen15', mem%pool%reserve%nitrogen15, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
 
      CALL slm_add_var(mem%pool%fruit, 'veg_pool_fruit_carbon', mem%pool%fruit%carbon, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%pool%fruit, 'veg_pool_fruit_nitrogen', mem%pool%fruit%nitrogen, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%pool%fruit, 'veg_pool_fruit_phosphorus', mem%pool%fruit%phosphorus, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%pool%fruit, 'veg_pool_fruit_carbon13', mem%pool%fruit%carbon13, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%pool%fruit, 'veg_pool_fruit_carbon14', mem%pool%fruit%carbon14, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%pool%fruit, 'veg_pool_fruit_nitrogen15', mem%pool%fruit%nitrogen15, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
 
      CALL slm_add_var(mem%pool%seed_bed, 'veg_pool_seed_bed_carbon', mem%pool%seed_bed%carbon, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%pool%seed_bed, 'veg_pool_seed_bed_nitrogen', mem%pool%seed_bed%nitrogen, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%pool%seed_bed, 'veg_pool_seed_bed_phosphorus', mem%pool%seed_bed%phosphorus, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%pool%seed_bed, 'veg_pool_seed_bed_carbon13', mem%pool%seed_bed%carbon13, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%pool%seed_bed, 'veg_pool_seed_bed_carbon14', mem%pool%seed_bed%carbon14, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%pool%seed_bed, 'veg_pool_seed_bed_nitrogen15', mem%pool%seed_bed%nitrogen15, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
 
      CALL slm_add_var(mem%growth%leaf, 'veg_growth_leaf_carbon', mem%growth%leaf%carbon, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%growth%leaf, 'veg_growth_leaf_nitrogen', mem%growth%leaf%nitrogen, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%growth%leaf, 'veg_growth_leaf_phosphorus', mem%growth%leaf%phosphorus, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%growth%leaf, 'veg_growth_leaf_carbon13', mem%growth%leaf%carbon13, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%growth%leaf, 'veg_growth_leaf_carbon14', mem%growth%leaf%carbon14, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%growth%leaf, 'veg_growth_leaf_nitrogen15', mem%growth%leaf%nitrogen15, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
 
      CALL slm_add_var(mem%growth%fine_root, 'veg_growth_fine_root_carbon', mem%growth%fine_root%carbon, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%growth%fine_root, 'veg_growth_fine_root_nitrogen', mem%growth%fine_root%nitrogen, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%growth%fine_root, 'veg_growth_fine_root_phosphorus', mem%growth%fine_root%phosphorus, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%growth%fine_root, 'veg_growth_fine_root_carbon13', mem%growth%fine_root%carbon13, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%growth%fine_root, 'veg_growth_fine_root_carbon14', mem%growth%fine_root%carbon14, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%growth%fine_root, 'veg_growth_fine_root_nitrogen15', mem%growth%fine_root%nitrogen15, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
 
      CALL slm_add_var(mem%growth%coarse_root, 'veg_growth_coarse_root_carbon', mem%growth%coarse_root%carbon, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%growth%coarse_root, 'veg_growth_coarse_root_nitrogen', mem%growth%coarse_root%nitrogen, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%growth%coarse_root, 'veg_growth_coarse_root_phosphorus', mem%growth%coarse_root%phosphorus, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%growth%coarse_root, 'veg_growth_coarse_root_carbon13', mem%growth%coarse_root%carbon13, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%growth%coarse_root, 'veg_growth_coarse_root_carbon14', mem%growth%coarse_root%carbon14, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%growth%coarse_root, 'veg_growth_coarse_root_nitrogen15', mem%growth%coarse_root%nitrogen15, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
 
      CALL slm_add_var(mem%growth%sap_wood, 'veg_growth_sap_wood_carbon', mem%growth%sap_wood%carbon, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%growth%sap_wood, 'veg_growth_sap_wood_nitrogen', mem%growth%sap_wood%nitrogen, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%growth%sap_wood, 'veg_growth_sap_wood_phosphorus', mem%growth%sap_wood%phosphorus, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%growth%sap_wood, 'veg_growth_sap_wood_carbon13', mem%growth%sap_wood%carbon13, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%growth%sap_wood, 'veg_growth_sap_wood_carbon14', mem%growth%sap_wood%carbon14, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%growth%sap_wood, 'veg_growth_sap_wood_nitrogen15', mem%growth%sap_wood%nitrogen15, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
 
      CALL slm_add_var(mem%growth%heart_wood, 'veg_growth_heart_wood_carbon', mem%growth%heart_wood%carbon, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%growth%heart_wood, 'veg_growth_heart_wood_nitrogen', mem%growth%heart_wood%nitrogen, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%growth%heart_wood, 'veg_growth_heart_wood_phosphorus', mem%growth%heart_wood%phosphorus, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%growth%heart_wood, 'veg_growth_heart_wood_carbon13', mem%growth%heart_wood%carbon13, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%growth%heart_wood, 'veg_growth_heart_wood_carbon14', mem%growth%heart_wood%carbon14, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%growth%heart_wood, 'veg_growth_heart_wood_nitrogen15', mem%growth%heart_wood%nitrogen15, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
 
      CALL slm_add_var(mem%growth%labile, 'veg_growth_labile_carbon', mem%growth%labile%carbon, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%growth%labile, 'veg_growth_labile_nitrogen', mem%growth%labile%nitrogen, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%growth%labile, 'veg_growth_labile_phosphorus', mem%growth%labile%phosphorus, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%growth%labile, 'veg_growth_labile_carbon13', mem%growth%labile%carbon13, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%growth%labile, 'veg_growth_labile_carbon14', mem%growth%labile%carbon14, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%growth%labile, 'veg_growth_labile_nitrogen15', mem%growth%labile%nitrogen15, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
 
      CALL slm_add_var(mem%growth%reserve, 'veg_growth_reserve_carbon', mem%growth%reserve%carbon, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%growth%reserve, 'veg_growth_reserve_nitrogen', mem%growth%reserve%nitrogen, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%growth%reserve, 'veg_growth_reserve_phosphorus', mem%growth%reserve%phosphorus, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%growth%reserve, 'veg_growth_reserve_carbon13', mem%growth%reserve%carbon13, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%growth%reserve, 'veg_growth_reserve_carbon14', mem%growth%reserve%carbon14, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%growth%reserve, 'veg_growth_reserve_nitrogen15', mem%growth%reserve%nitrogen15, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
 
      CALL slm_add_var(mem%growth%fruit, 'veg_growth_fruit_carbon', mem%growth%fruit%carbon, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%growth%fruit, 'veg_growth_fruit_nitrogen', mem%growth%fruit%nitrogen, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%growth%fruit, 'veg_growth_fruit_phosphorus', mem%growth%fruit%phosphorus, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%growth%fruit, 'veg_growth_fruit_carbon13', mem%growth%fruit%carbon13, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%growth%fruit, 'veg_growth_fruit_carbon14', mem%growth%fruit%carbon14, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%growth%fruit, 'veg_growth_fruit_nitrogen15', mem%growth%fruit%nitrogen15, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
 
      CALL slm_add_var(mem%growth%seed_bed, 'veg_growth_seed_bed_carbon', mem%growth%seed_bed%carbon, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%growth%seed_bed, 'veg_growth_seed_bed_nitrogen', mem%growth%seed_bed%nitrogen, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%growth%seed_bed, 'veg_growth_seed_bed_phosphorus', mem%growth%seed_bed%phosphorus, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%growth%seed_bed, 'veg_growth_seed_bed_carbon13', mem%growth%seed_bed%carbon13, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%growth%seed_bed, 'veg_growth_seed_bed_carbon14', mem%growth%seed_bed%carbon14, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%growth%seed_bed, 'veg_growth_seed_bed_nitrogen15', mem%growth%seed_bed%nitrogen15, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
 
      CALL slm_add_var(mem%litterfall%leaf, 'veg_litterfall_leaf_carbon', mem%litterfall%leaf%carbon, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%litterfall%leaf, 'veg_litterfall_leaf_nitrogen', mem%litterfall%leaf%nitrogen, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%litterfall%leaf, 'veg_litterfall_leaf_phosphorus', mem%litterfall%leaf%phosphorus, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%litterfall%leaf, 'veg_litterfall_leaf_carbon13', mem%litterfall%leaf%carbon13, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%litterfall%leaf, 'veg_litterfall_leaf_carbon14', mem%litterfall%leaf%carbon14, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%litterfall%leaf, 'veg_litterfall_leaf_nitrogen15', mem%litterfall%leaf%nitrogen15, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
 
      CALL slm_add_var(mem%litterfall%fine_root, 'veg_litterfall_fine_root_carbon', mem%litterfall%fine_root%carbon, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%litterfall%fine_root, 'veg_litterfall_fine_root_nitrogen', mem%litterfall%fine_root%nitrogen, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%litterfall%fine_root, 'veg_litterfall_fine_root_phosphorus', mem%litterfall%fine_root%phosphorus, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%litterfall%fine_root, 'veg_litterfall_fine_root_carbon13', mem%litterfall%fine_root%carbon13, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%litterfall%fine_root, 'veg_litterfall_fine_root_carbon14', mem%litterfall%fine_root%carbon14, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%litterfall%fine_root, 'veg_litterfall_fine_root_nitrogen15', mem%litterfall%fine_root%nitrogen15, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
 
      CALL slm_add_var(mem%litterfall%coarse_root, 'veg_litterfall_coarse_root_carbon', mem%litterfall%coarse_root%carbon, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%litterfall%coarse_root, 'veg_litterfall_coarse_root_nitrogen', mem%litterfall%coarse_root%nitrogen, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%litterfall%coarse_root, 'veg_litterfall_coarse_root_phosphorus', mem%litterfall%coarse_root%phosphorus, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%litterfall%coarse_root, 'veg_litterfall_coarse_root_carbon13', mem%litterfall%coarse_root%carbon13, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%litterfall%coarse_root, 'veg_litterfall_coarse_root_carbon14', mem%litterfall%coarse_root%carbon14, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%litterfall%coarse_root, 'veg_litterfall_coarse_root_nitrogen15', mem%litterfall%coarse_root%nitrogen15, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
 
      CALL slm_add_var(mem%litterfall%sap_wood, 'veg_litterfall_sap_wood_carbon', mem%litterfall%sap_wood%carbon, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%litterfall%sap_wood, 'veg_litterfall_sap_wood_nitrogen', mem%litterfall%sap_wood%nitrogen, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%litterfall%sap_wood, 'veg_litterfall_sap_wood_phosphorus', mem%litterfall%sap_wood%phosphorus, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%litterfall%sap_wood, 'veg_litterfall_sap_wood_carbon13', mem%litterfall%sap_wood%carbon13, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%litterfall%sap_wood, 'veg_litterfall_sap_wood_carbon14', mem%litterfall%sap_wood%carbon14, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%litterfall%sap_wood, 'veg_litterfall_sap_wood_nitrogen15', mem%litterfall%sap_wood%nitrogen15, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
 
      CALL slm_add_var(mem%litterfall%heart_wood, 'veg_litterfall_heart_wood_carbon', mem%litterfall%heart_wood%carbon, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%litterfall%heart_wood, 'veg_litterfall_heart_wood_nitrogen', mem%litterfall%heart_wood%nitrogen, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%litterfall%heart_wood, 'veg_litterfall_heart_wood_phosphorus', mem%litterfall%heart_wood%phosphorus, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%litterfall%heart_wood, 'veg_litterfall_heart_wood_carbon13', mem%litterfall%heart_wood%carbon13, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%litterfall%heart_wood, 'veg_litterfall_heart_wood_carbon14', mem%litterfall%heart_wood%carbon14, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%litterfall%heart_wood, 'veg_litterfall_heart_wood_nitrogen15', mem%litterfall%heart_wood%nitrogen15, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
 
      CALL slm_add_var(mem%litterfall%labile, 'veg_litterfall_labile_carbon', mem%litterfall%labile%carbon, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%litterfall%labile, 'veg_litterfall_labile_nitrogen', mem%litterfall%labile%nitrogen, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%litterfall%labile, 'veg_litterfall_labile_phosphorus', mem%litterfall%labile%phosphorus, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%litterfall%labile, 'veg_litterfall_labile_carbon13', mem%litterfall%labile%carbon13, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%litterfall%labile, 'veg_litterfall_labile_carbon14', mem%litterfall%labile%carbon14, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%litterfall%labile, 'veg_litterfall_labile_nitrogen15', mem%litterfall%labile%nitrogen15, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
 
      CALL slm_add_var(mem%litterfall%reserve, 'veg_litterfall_reserve_carbon', mem%litterfall%reserve%carbon, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%litterfall%reserve, 'veg_litterfall_reserve_nitrogen', mem%litterfall%reserve%nitrogen, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%litterfall%reserve, 'veg_litterfall_reserve_phosphorus', mem%litterfall%reserve%phosphorus, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%litterfall%reserve, 'veg_litterfall_reserve_carbon13', mem%litterfall%reserve%carbon13, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%litterfall%reserve, 'veg_litterfall_reserve_carbon14', mem%litterfall%reserve%carbon14, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%litterfall%reserve, 'veg_litterfall_reserve_nitrogen15', mem%litterfall%reserve%nitrogen15, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
 
      CALL slm_add_var(mem%litterfall%fruit, 'veg_litterfall_fruit_carbon', mem%litterfall%fruit%carbon, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%litterfall%fruit, 'veg_litterfall_fruit_nitrogen', mem%litterfall%fruit%nitrogen, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%litterfall%fruit, 'veg_litterfall_fruit_phosphorus', mem%litterfall%fruit%phosphorus, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%litterfall%fruit, 'veg_litterfall_fruit_carbon13', mem%litterfall%fruit%carbon13, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%litterfall%fruit, 'veg_litterfall_fruit_carbon14', mem%litterfall%fruit%carbon14, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%litterfall%fruit, 'veg_litterfall_fruit_nitrogen15', mem%litterfall%fruit%nitrogen15, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
 
      CALL slm_add_var(mem%litterfall%seed_bed, 'veg_litterfall_seed_bed_carbon', mem%litterfall%seed_bed%carbon, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%litterfall%seed_bed, 'veg_litterfall_seed_bed_nitrogen', mem%litterfall%seed_bed%nitrogen, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%litterfall%seed_bed, 'veg_litterfall_seed_bed_phosphorus', mem%litterfall%seed_bed%phosphorus, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%litterfall%seed_bed, 'veg_litterfall_seed_bed_carbon13', mem%litterfall%seed_bed%carbon13, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%litterfall%seed_bed, 'veg_litterfall_seed_bed_carbon14', mem%litterfall%seed_bed%carbon14, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%litterfall%seed_bed, 'veg_litterfall_seed_bed_nitrogen15', mem%litterfall%seed_bed%nitrogen15, &
            ics, ice, iblk, &     ! model%options -> for SLM compatibility
            initval_r=0.0_wp)     ! initvalue

    END IF

  END SUBROUTINE Init_veg_memory

END MODULE mo_veg_memory_class
