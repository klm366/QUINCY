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

  USE mo_kind,                   ONLY: wp
  USE mo_exception,              ONLY: message
  USE mo_jsb_class,              ONLY: Get_model

  USE mo_jsb_memory_class,       ONLY: t_jsb_memory
  USE mo_jsb_var_class,          ONLY: t_jsb_var_real2d, t_jsb_var_real3d, &
                                       t_jsb_struct, t_jsb_n_struct             ! *slmjsb4

  USE mo_slm_varlist,            ONLY: slm_add_var   ! replace by mo_jsb_varlist later    *slmtemp



  IMPLICIT NONE
  PRIVATE
  PUBLIC :: t_soil_memory, max_no_of_vars
  PUBLIC :: t_slm_struct_soil_element, t_slm_struct_soil_pool   ! *slmjsb4
  PUBLIC :: max_no_of_structs, max_no_of_n_structs              ! *slmjsb4

  INTEGER, PARAMETER :: max_no_of_vars      = 80   !> TYPE t_jsb_var
  INTEGER, PARAMETER :: max_no_of_structs   = 20   !> TYPE t_jsb_struct     *slmjsb4
  INTEGER, PARAMETER :: max_no_of_n_structs = 20   !> TYPE t_jsb_n_struct   *slmjsb4



  !-----------------------------------------------------------------------------------------------------
  !> Type definition for elements in soil pool
  !! 
  !! elements in soil are 3D (and would be 2D in veg)
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
    PRIVATE
    !! operator overload
    ! +
    PROCEDURE               :: Add_struct                   ! Add two variables of type t_slm_struct_soil_element
    PROCEDURE               :: Add_struct_and_scalar        ! Add variable of type t_slm_struct_soil_element with REAL scalar
    PROCEDURE               :: Add_struct_and_field3d       ! Add variable of type t_slm_struct_soil_element with REAL 2-d field2d
    GENERIC, PUBLIC         :: OPERATOR(+) => Add_struct, Add_struct_and_scalar, Add_struct_and_field3d
    ! -
    PROCEDURE               :: Subtract_struct              ! Subtract two variables of type t_slm_struct_soil_element
    PROCEDURE               :: Subtract_struct_and_scalar   ! Subtract variable of type t_slm_struct_soil_element with REAL scalar
    PROCEDURE               :: Subtract_struct_and_field3d  ! Subtract variable of type t_slm_struct_soil_element with REAL 2-d field2d
    GENERIC, PUBLIC         :: OPERATOR(-) => Subtract_struct, Subtract_struct_and_scalar, Subtract_struct_and_field3d
    ! *
    PROCEDURE               :: Multiply_struct              ! Multiply two variables of type t_slm_struct_soil_element
    PROCEDURE               :: Multiply_struct_and_scalar   ! Multiply variable of type t_slm_struct_soil_element with REAL scalar
    PROCEDURE               :: Multiply_struct_and_field3d  ! Multiply variable of type t_slm_struct_soil_element with REAL 2-d field2d
    GENERIC, PUBLIC         :: OPERATOR(*) => Multiply_struct, Multiply_struct_and_scalar, Multiply_struct_and_field3d
    ! /
    PROCEDURE               :: Divide_struct                ! Divide two variables of type t_slm_struct_soil_element
    PROCEDURE               :: Divide_struct_and_scalar     ! Divide variable of type t_slm_struct_soil_element by REAL scalar
    PROCEDURE               :: Divide_struct_and_field3d    ! Divide variable of type t_slm_struct_soil_element by REAL 2-d field2d
    GENERIC, PUBLIC         :: OPERATOR(/) => Divide_struct, Divide_struct_and_scalar, Divide_struct_and_field3d
    !! initialize with/ set to a certain REAL
    PROCEDURE, PUBLIC       :: init => init_t_slm_struct_soil_element
  END TYPE t_slm_struct_soil_element


  !-----------------------------------------------------------------------------------------------------
  !> Type definition for soil pool with nested TYPEs
  !!
  !! used for: pool ...
  !!  
  !-----------------------------------------------------------------------------------------------------
  TYPE, EXTENDS(t_jsb_n_struct) :: t_slm_struct_soil_pool
    TYPE(t_slm_struct_soil_element)  :: dom                 !< dissolved organic matter
    TYPE(t_slm_struct_soil_element)  :: dom_assoc           !< minerally associated disolved organic matter
    TYPE(t_slm_struct_soil_element)  :: metabolic_litter    !< metabolic litter (sugars, cellulose etc.)
    TYPE(t_slm_struct_soil_element)  :: structural_litter   !< structural litter (lignified litter)
    TYPE(t_slm_struct_soil_element)  :: woody_litter        !< physically protected wood litter
    TYPE(t_slm_struct_soil_element)  :: microbial           !< microbial biomass
    TYPE(t_slm_struct_soil_element)  :: residue             !< microbial necromass
    TYPE(t_slm_struct_soil_element)  :: residue_assoc       !< minerally associated microbial necromass 
  CONTAINS
    PRIVATE
    !! operator overload
    ! +
    PROCEDURE               :: Add_n_struct                   
    PROCEDURE               :: Add_n_struct_and_scalar        
    PROCEDURE               :: Add_n_struct_and_field3d       
    GENERIC, PUBLIC         :: OPERATOR(+) => Add_n_struct, Add_n_struct_and_scalar, Add_n_struct_and_field3d
    ! -
    PROCEDURE               :: Subtract_n_struct              
    PROCEDURE               :: Subtract_n_struct_and_scalar   
    PROCEDURE               :: Subtract_n_struct_and_field3d  
    GENERIC, PUBLIC         :: OPERATOR(-) => Subtract_n_struct, Subtract_n_struct_and_scalar, Subtract_n_struct_and_field3d
    ! *
    PROCEDURE               :: Multiply_n_struct              
    PROCEDURE               :: Multiply_n_struct_and_scalar   
    PROCEDURE               :: Multiply_n_struct_and_field3d  
    GENERIC, PUBLIC         :: OPERATOR(*) => Multiply_n_struct, Multiply_n_struct_and_scalar, Multiply_n_struct_and_field3d
    ! /
    PROCEDURE               :: Divide_n_struct                
    PROCEDURE               :: Divide_n_struct_and_scalar     
    PROCEDURE               :: Divide_n_struct_and_field3d    
    GENERIC, PUBLIC         :: OPERATOR(/) => Divide_n_struct, Divide_n_struct_and_scalar, Divide_n_struct_and_field3d
    !! initialize with/ set to a certain REAL
    PROCEDURE, PUBLIC       :: init => init_t_slm_struct_soil_pool
  END TYPE t_slm_struct_soil_pool


  !-----------------------------------------------------------------------------------------------------
  !> Type definition for soil memory
  !!
  !! @par includes: \n
  !!    soil variables soil\%VAR \n
  !!    pool  \n
  !!    ? flux \n
  !!    ? flux \n
  !!    ? moving averagages of variables VAR_mavg   \n
  !!    ? daytime averages of variables VAR_daytime \n
  !! 
  !-----------------------------------------------------------------------------------------------------
  TYPE, EXTENDS(t_jsb_memory) :: t_soil_memory


    TYPE(t_jsb_var_real2d)          :: &
                                    soil_depth          , & !< Soil depth derived from textures (bedrock) (for hydrology) [m]
                                    root_depth              !< actual rooting depth [m]


! JSBACH SOIL
    TYPE(t_jsb_var_real3d)          :: &
                                    soil_depth_l            !< Depth of each soil layer that can be saturated with water (until bedrock) [m]

! JSBACH HYDRO MEMORY (temporary here to get rid of vegsoil construct)
    TYPE(t_jsb_var_real2d) :: &
      & runoff,           & !< surface runoff
      & drainage            !< drainage

    ! Additional variables for soil
    TYPE(t_jsb_var_real2d) :: &
      & w_soil_root,      & !< Water content in root zone of the soil                                   [m]
      & w_soil_root_fc,   & !< Water content at field capacity in root zone of the soil                 [m]
      & w_soil_root_pwp     !< Water content at permanent wilting point in root zone of the soil        [m]

    TYPE(t_jsb_var_real3d) :: &
      & drainage_sl,      & !< drainage of soil layers                                                  [m]
      & w_soil_sl,        & !< Water content in soil layers                                             [m]
      & w_soil_fc_sl,     & !< Water content of soil layers at field capacity
                            !! (derived from volumetric field capacity)                                 [m]
      & w_soil_pwp_sl       !< Water content of soil layers at permanent wilting point
                            !! (derived from vol. perm wilt point)                                      [m]

!Soil Water addons for SLM in HYDRO
    TYPE(t_jsb_var_real2d) :: &
      & w_soil_root_theta, & !< root zone integrated plant available water (fraction of maximum) 
      & w_soil_root_pot      !< soil water potential (MPa)

    TYPE(t_jsb_var_real3d) :: &
      & w_soil_pot_sl        !< soil water potential per layer (MPa)

! should go to SLM VEG after change to 2D/3D
    TYPE(t_jsb_var_real3d)          :: &
                                    root_fraction_l         !< fractions of root per layer

! SLM SOILBGC
    TYPE(t_jsb_var_real3d)          :: &
                                    nh4_solute         , &  !< ammonium in soil solution
                                    nh4_assoc          , &  !< minerally associated ammonium
                                    no3_solute         , &  !< nitrate in soil solution
                                    noy                , &  !< NO_y
                                    n2o                , &  !< N2O
                                    n2                 , &  !< N2
                                    nh4_n15_solute     , &  !< N15-ammonium in soil solution
                                    nh4_n15_assoc      , &  !< minerally associated N15-ammonium
                                    no3_n15_solute     , &  !< N15-nitrate in soil solution
                                    noy_n15            , &  !< N15-NO_y
                                    n2o_n15            , &  !< N15-N2O
                                    n2_n15             , &  !< N15-N2
                                    po4_assoc          , &  !< phosphate in soil solution
                                    po4_solute              !< minerally associated phosphate 

    ! soil energy 
    TYPE(t_jsb_var_real3d)          :: &
                                    t_soil                  !< soil temperature  *slmdev currently we use t_air to get values here

    ! soil physical properties (not sure where they should belong to, but they are needed in the long-term
    TYPE(t_jsb_var_real3d)          :: & 
                                    heat_capa_l,         &  !< soil heat capacity (J / kg )
                                    bulk_dens_l,         &  !< soil bulk density (kg / m3)
                                    therm_cond_l            !< soil thermal conductivity (W / m2 / K)

!    TYPE(t_slm_struct_soil_element) :: &   ! simple structure for elements: t_jsb_struct
!                                    test_soil_element   !< test ...



    TYPE(t_slm_struct_soil_pool)    :: &   ! nested pool TYPE t_jsb_n_struct
                                    pool                !< pools of the soil and living fellows of it




  CONTAINS
    PROCEDURE :: Init => Init_soil_memory
  END TYPE t_soil_memory


  CHARACTER(len=*), PARAMETER :: modname = 'mo_soil_memory_class'


CONTAINS

  !-----------------------------------------------------------------------------------------------------!
  ! t_slm_struct_soil_element operator overload                                                         !
  !-----------------------------------------------------------------------------------------------------!

  !-----------------------------------------------------------------------------------------------------
  !> + (t_slm_struct_soil_element + t_slm_struct_soil_element)
  !!
  !-----------------------------------------------------------------------------------------------------
  FUNCTION Add_struct(this,struct2) RESULT(struct3)
    IMPLICIT NONE
    ! 0.1 InOut
    CLASS(t_slm_struct_soil_element), INTENT(in) :: this, struct2  !< t_slm_struct_soil_element in
    TYPE(t_slm_struct_soil_element), ALLOCATABLE :: struct3        !< t_slm_struct_soil_element out

    ALLOCATE(struct3, source=this)

    struct3%carbon       = this%carbon      + struct2%carbon
    struct3%nitrogen     = this%nitrogen    + struct2%nitrogen
    struct3%phosphorus   = this%phosphorus  + struct2%phosphorus
    struct3%carbon13     = this%carbon13    + struct2%carbon13
    struct3%carbon14     = this%carbon14    + struct2%carbon14
    struct3%nitrogen15   = this%nitrogen15  + struct2%nitrogen15
  END FUNCTION Add_struct

  !-----------------------------------------------------------------------------------------------------
  !> + (t_slm_struct_soil_element + scalar)
  !!
  !-----------------------------------------------------------------------------------------------------
  FUNCTION Add_struct_and_scalar(this,scalar) RESULT(struct3)
    IMPLICIT NONE
    ! 0.1 InOut
    CLASS(t_slm_struct_soil_element), INTENT(in) :: this
    REAL(wp),                        INTENT(in) :: scalar
    TYPE(t_slm_struct_soil_element), ALLOCATABLE :: struct3

    ALLOCATE(struct3, source=this)

    struct3%carbon       = this%carbon      + scalar
    struct3%nitrogen     = this%nitrogen    + scalar
    struct3%phosphorus   = this%phosphorus  + scalar
    struct3%carbon13     = this%carbon13    + scalar
    struct3%carbon14     = this%carbon14    + scalar
    struct3%nitrogen15   = this%nitrogen15  + scalar
  END FUNCTION Add_struct_and_scalar

  !-----------------------------------------------------------------------------------------------------
  !> + (t_slm_struct_soil_element + field3d)
  !!
  !-----------------------------------------------------------------------------------------------------
  FUNCTION Add_struct_and_field3d(this,field3d) RESULT(struct3)
    IMPLICIT NONE
    ! 0.1 InOut
    CLASS(t_slm_struct_soil_element), INTENT(in) :: this
    REAL(wp),                        INTENT(in) :: field3d(:,:,:)
    TYPE(t_slm_struct_soil_element), ALLOCATABLE :: struct3

    ALLOCATE(struct3, source=this)

    struct3%carbon       = this%carbon      + field3d
    struct3%nitrogen     = this%nitrogen    + field3d
    struct3%phosphorus   = this%phosphorus  + field3d
    struct3%carbon13     = this%carbon13    + field3d
    struct3%carbon14     = this%carbon14    + field3d
    struct3%nitrogen15   = this%nitrogen15  + field3d
  END FUNCTION Add_struct_and_field3d

  !-----------------------------------------------------------------------------------------------------
  !> - (t_slm_struct_soil_element - t_slm_struct_soil_element)
  !!
  !-----------------------------------------------------------------------------------------------------
  FUNCTION Subtract_struct(this,struct2) RESULT(struct3)
    IMPLICIT NONE
    ! 0.1 InOut
    CLASS(t_slm_struct_soil_element), INTENT(in) :: this, struct2  !< t_slm_struct_soil_element in
    TYPE(t_slm_struct_soil_element), ALLOCATABLE :: struct3        !< t_slm_struct_soil_element out

    ALLOCATE(struct3, source=this)

    struct3%carbon       = this%carbon      - struct2%carbon
    struct3%nitrogen     = this%nitrogen    - struct2%nitrogen
    struct3%phosphorus   = this%phosphorus  - struct2%phosphorus
    struct3%carbon13     = this%carbon13    - struct2%carbon13
    struct3%carbon14     = this%carbon14    - struct2%carbon14
    struct3%nitrogen15   = this%nitrogen15  - struct2%nitrogen15
  END FUNCTION Subtract_struct

  !-----------------------------------------------------------------------------------------------------
  !> - (t_slm_struct_soil_element - scalar)
  !!
  !-----------------------------------------------------------------------------------------------------
  FUNCTION Subtract_struct_and_scalar(this,scalar) RESULT(struct3)
    IMPLICIT NONE
    ! 0.1 InOut
    CLASS(t_slm_struct_soil_element), INTENT(in) :: this
    REAL(wp),                        INTENT(in) :: scalar
    TYPE(t_slm_struct_soil_element), ALLOCATABLE :: struct3

    ALLOCATE(struct3, source=this)

    struct3%carbon       = this%carbon      - scalar
    struct3%nitrogen     = this%nitrogen    - scalar
    struct3%phosphorus   = this%phosphorus  - scalar
    struct3%carbon13     = this%carbon13    - scalar
    struct3%carbon14     = this%carbon14    - scalar
    struct3%nitrogen15   = this%nitrogen15  - scalar
  END FUNCTION Subtract_struct_and_scalar

  !-----------------------------------------------------------------------------------------------------
  !> - (t_slm_struct_soil_element - field3d)
  !!
  !-----------------------------------------------------------------------------------------------------
  FUNCTION Subtract_struct_and_field3d(this,field3d) RESULT(struct3)
    IMPLICIT NONE
    ! 0.1 InOut
    CLASS(t_slm_struct_soil_element), INTENT(in) :: this
    REAL(wp),                        INTENT(in) :: field3d(:,:,:)
    TYPE(t_slm_struct_soil_element), ALLOCATABLE :: struct3

    ALLOCATE(struct3, source=this)

    struct3%carbon       = this%carbon      - field3d
    struct3%nitrogen     = this%nitrogen    - field3d
    struct3%phosphorus   = this%phosphorus  - field3d
    struct3%carbon13     = this%carbon13    - field3d
    struct3%carbon14     = this%carbon14    - field3d
    struct3%nitrogen15   = this%nitrogen15  - field3d
  END FUNCTION Subtract_struct_and_field3d

  !-----------------------------------------------------------------------------------------------------
  !> * (t_slm_struct_soil_element * t_slm_struct_soil_element)
  !!
  !-----------------------------------------------------------------------------------------------------
  FUNCTION Multiply_struct(this,struct2) RESULT(struct3)
    IMPLICIT NONE
    ! 0.1 InOut
    CLASS(t_slm_struct_soil_element), INTENT(in) :: this, struct2  !< t_slm_struct_soil_element in
    TYPE(t_slm_struct_soil_element), ALLOCATABLE :: struct3        !< t_slm_struct_soil_element out

    ALLOCATE(struct3, source=this)

    struct3%carbon       = this%carbon      * struct2%carbon
    struct3%nitrogen     = this%nitrogen    * struct2%nitrogen
    struct3%phosphorus   = this%phosphorus  * struct2%phosphorus
    struct3%carbon13     = this%carbon13    * struct2%carbon13
    struct3%carbon14     = this%carbon14    * struct2%carbon14
    struct3%nitrogen15   = this%nitrogen15  * struct2%nitrogen15
  END FUNCTION Multiply_struct

  !-----------------------------------------------------------------------------------------------------
  !> * (t_slm_struct_soil_element * scalar)
  !!
  !-----------------------------------------------------------------------------------------------------
  FUNCTION Multiply_struct_and_scalar(this,scalar) RESULT(struct3)
    IMPLICIT NONE
    ! 0.1 InOut
    CLASS(t_slm_struct_soil_element), INTENT(in) :: this
    REAL(wp),                        INTENT(in) :: scalar
    TYPE(t_slm_struct_soil_element), ALLOCATABLE :: struct3

    ALLOCATE(struct3, source=this)

    struct3%carbon       = this%carbon      * scalar
    struct3%nitrogen     = this%nitrogen    * scalar
    struct3%phosphorus   = this%phosphorus  * scalar
    struct3%carbon13     = this%carbon13    * scalar
    struct3%carbon14     = this%carbon14    * scalar
    struct3%nitrogen15   = this%nitrogen15  * scalar
  END FUNCTION Multiply_struct_and_scalar

  !-----------------------------------------------------------------------------------------------------
  !> * (t_slm_struct_soil_element * field3d)
  !!
  !-----------------------------------------------------------------------------------------------------
  FUNCTION Multiply_struct_and_field3d(this,field3d) RESULT(struct3)
    IMPLICIT NONE
    ! 0.1 InOut
    CLASS(t_slm_struct_soil_element), INTENT(in) :: this
    REAL(wp),                        INTENT(in) :: field3d(:,:,:)
    TYPE(t_slm_struct_soil_element), ALLOCATABLE :: struct3

    ALLOCATE(struct3, source=this)

    struct3%carbon       = this%carbon      * field3d
    struct3%nitrogen     = this%nitrogen    * field3d
    struct3%phosphorus   = this%phosphorus  * field3d
    struct3%carbon13     = this%carbon13    * field3d
    struct3%carbon14     = this%carbon14    * field3d
    struct3%nitrogen15   = this%nitrogen15  * field3d
  END FUNCTION Multiply_struct_and_field3d

  !-----------------------------------------------------------------------------------------------------
  !> / (t_slm_struct_soil_element / t_slm_struct_soil_element)
  !!
  !-----------------------------------------------------------------------------------------------------
  FUNCTION Divide_struct(this,struct2) RESULT(struct3)
    IMPLICIT NONE
    ! 0.1 InOut
    CLASS(t_slm_struct_soil_element), INTENT(in) :: this, struct2 
    TYPE(t_slm_struct_soil_element), ALLOCATABLE :: struct3       

    ALLOCATE(struct3, source=this)

    struct3%carbon       = this%carbon      / struct2%carbon
    struct3%nitrogen     = this%nitrogen    / struct2%nitrogen
    struct3%phosphorus   = this%phosphorus  / struct2%phosphorus
    struct3%carbon13     = this%carbon13    / struct2%carbon13
    struct3%carbon14     = this%carbon14    / struct2%carbon14
    struct3%nitrogen15   = this%nitrogen15  / struct2%nitrogen15
  END FUNCTION Divide_struct

  !-----------------------------------------------------------------------------------------------------
  !> / (t_slm_struct_soil_element / scalar)
  !!
  !-----------------------------------------------------------------------------------------------------
  FUNCTION Divide_struct_and_scalar(this,scalar) RESULT(struct3)
    IMPLICIT NONE
    ! 0.1 InOut
    CLASS(t_slm_struct_soil_element), INTENT(in) :: this
    REAL(wp),                        INTENT(in) :: scalar
    TYPE(t_slm_struct_soil_element), ALLOCATABLE :: struct3

    ALLOCATE(struct3, source=this)

    struct3%carbon       = this%carbon      / scalar
    struct3%nitrogen     = this%nitrogen    / scalar
    struct3%phosphorus   = this%phosphorus  / scalar
    struct3%carbon13     = this%carbon13    / scalar
    struct3%carbon14     = this%carbon14    / scalar
    struct3%nitrogen15   = this%nitrogen15  / scalar
  END FUNCTION Divide_struct_and_scalar

  !-----------------------------------------------------------------------------------------------------
  !> / (t_slm_struct_soil_element / field3d)
  !!
  !-----------------------------------------------------------------------------------------------------
  FUNCTION Divide_struct_and_field3d(this,field3d) RESULT(struct3)
    IMPLICIT NONE
    ! 0.1 InOut
    CLASS(t_slm_struct_soil_element), INTENT(in) :: this
    REAL(wp),                        INTENT(in) :: field3d(:,:,:)
    TYPE(t_slm_struct_soil_element), ALLOCATABLE :: struct3

    ALLOCATE(struct3, source=this)

    struct3%carbon       = this%carbon      / field3d
    struct3%nitrogen     = this%nitrogen    / field3d
    struct3%phosphorus   = this%phosphorus  / field3d
    struct3%carbon13     = this%carbon13    / field3d
    struct3%carbon14     = this%carbon14    / field3d
    struct3%nitrogen15   = this%nitrogen15  / field3d
  END FUNCTION Divide_struct_and_field3d

  !-----------------------------------------------------------------------------------------------------
  !> Initialisation t_slm_struct_soil_element with a scalar
  !!
  !-----------------------------------------------------------------------------------------------------
  PURE ELEMENTAL SUBROUTINE init_t_slm_struct_soil_element(this,s)
    IMPLICIT NONE
    ! 0.1 InOut
    CLASS(t_slm_struct_soil_element),INTENT(inout) :: this 
    REAL(wp),                       INTENT(in)    :: s

    CALL this%carbon%init(s) 
    CALL this%nitrogen%init(s) 
    CALL this%phosphorus%init(s) 
    CALL this%carbon13%init(s) 
    CALL this%carbon14%init(s) 
    CALL this%nitrogen15%init(s) 
  END SUBROUTINE init_t_slm_struct_soil_element 



  !-----------------------------------------------------------------------------------------------------!
  ! t_slm_struct_soil_pool operator overload                                                            !
  !-----------------------------------------------------------------------------------------------------!

  !-----------------------------------------------------------------------------------------------------
  !> + (t_slm_struct_soil_pool + t_slm_struct_soil_pool)
  !!
  !-----------------------------------------------------------------------------------------------------
  FUNCTION Add_n_struct(this,struct2) RESULT(struct3)
    IMPLICIT NONE
    ! 0.1 InOut
    CLASS(t_slm_struct_soil_pool), INTENT(in) :: this, struct2  !< t_slm_struct_soil_pool in
    TYPE(t_slm_struct_soil_pool), ALLOCATABLE :: struct3        !< t_slm_struct_soil_pool out

    ALLOCATE(struct3, source=this)

    struct3%dom       = this%dom      + struct2%dom
  END FUNCTION Add_n_struct

  !-----------------------------------------------------------------------------------------------------
  !> + (t_slm_struct_soil_pool + scalar)
  !!
  !-----------------------------------------------------------------------------------------------------
  FUNCTION Add_n_struct_and_scalar(this,scalar) RESULT(struct3)
    IMPLICIT NONE
    ! 0.1 InOut
    CLASS(t_slm_struct_soil_pool), INTENT(in) :: this
    REAL(wp),                        INTENT(in) :: scalar
    TYPE(t_slm_struct_soil_pool), ALLOCATABLE :: struct3

    ALLOCATE(struct3, source=this)

    struct3%dom       = this%dom      + scalar
  END FUNCTION Add_n_struct_and_scalar

  !-----------------------------------------------------------------------------------------------------
  !> + (t_slm_struct_soil_pool + field3d)
  !!
  !-----------------------------------------------------------------------------------------------------
  FUNCTION Add_n_struct_and_field3d(this,field3d) RESULT(struct3)
    IMPLICIT NONE
    ! 0.1 InOut
    CLASS(t_slm_struct_soil_pool), INTENT(in) :: this
    REAL(wp),                        INTENT(in) :: field3d(:,:,:)
    TYPE(t_slm_struct_soil_pool), ALLOCATABLE :: struct3

    ALLOCATE(struct3, source=this)

    struct3%dom       = this%dom      + field3d
  END FUNCTION Add_n_struct_and_field3d

  !-----------------------------------------------------------------------------------------------------
  !> - (t_slm_struct_soil_pool - t_slm_struct_soil_pool)
  !!
  !-----------------------------------------------------------------------------------------------------
  FUNCTION Subtract_n_struct(this,struct2) RESULT(struct3)
    IMPLICIT NONE
    ! 0.1 InOut
    CLASS(t_slm_struct_soil_pool), INTENT(in) :: this, struct2  !< t_slm_struct_soil_pool in
    TYPE(t_slm_struct_soil_pool), ALLOCATABLE :: struct3        !< t_slm_struct_soil_pool out

    ALLOCATE(struct3, source=this)

    struct3%dom       = this%dom      - struct2%dom
  END FUNCTION Subtract_n_struct

  !-----------------------------------------------------------------------------------------------------
  !> - (t_slm_struct_soil_pool - scalar)
  !!
  !-----------------------------------------------------------------------------------------------------
  FUNCTION Subtract_n_struct_and_scalar(this,scalar) RESULT(struct3)
    IMPLICIT NONE
    ! 0.1 InOut
    CLASS(t_slm_struct_soil_pool), INTENT(in) :: this
    REAL(wp),                        INTENT(in) :: scalar
    TYPE(t_slm_struct_soil_pool), ALLOCATABLE :: struct3

    ALLOCATE(struct3, source=this)

    struct3%dom       = this%dom      - scalar
  END FUNCTION Subtract_n_struct_and_scalar

  !-----------------------------------------------------------------------------------------------------
  !> - (t_slm_struct_soil_pool - field3d)
  !!
  !-----------------------------------------------------------------------------------------------------
  FUNCTION Subtract_n_struct_and_field3d(this,field3d) RESULT(struct3)
    IMPLICIT NONE
    ! 0.1 InOut
    CLASS(t_slm_struct_soil_pool), INTENT(in) :: this
    REAL(wp),                        INTENT(in) :: field3d(:,:,:)
    TYPE(t_slm_struct_soil_pool), ALLOCATABLE :: struct3

    ALLOCATE(struct3, source=this)

    struct3%dom       = this%dom      - field3d
  END FUNCTION Subtract_n_struct_and_field3d

  !-----------------------------------------------------------------------------------------------------
  !> * (t_slm_struct_soil_pool * t_slm_struct_soil_pool)
  !!
  !-----------------------------------------------------------------------------------------------------
  FUNCTION Multiply_n_struct(this,struct2) RESULT(struct3)
    IMPLICIT NONE
    ! 0.1 InOut
    CLASS(t_slm_struct_soil_pool), INTENT(in) :: this, struct2  !< t_slm_struct_soil_pool in
    TYPE(t_slm_struct_soil_pool), ALLOCATABLE :: struct3        !< t_slm_struct_soil_pool out

    ALLOCATE(struct3, source=this)

    struct3%dom       = this%dom      * struct2%dom
  END FUNCTION Multiply_n_struct

  !-----------------------------------------------------------------------------------------------------
  !> * (t_slm_struct_soil_pool * scalar)
  !!
  !-----------------------------------------------------------------------------------------------------
  FUNCTION Multiply_n_struct_and_scalar(this,scalar) RESULT(struct3)
    IMPLICIT NONE
    ! 0.1 InOut
    CLASS(t_slm_struct_soil_pool), INTENT(in) :: this
    REAL(wp),                        INTENT(in) :: scalar
    TYPE(t_slm_struct_soil_pool), ALLOCATABLE :: struct3

    ALLOCATE(struct3, source=this)

    struct3%dom       = this%dom      * scalar
  END FUNCTION Multiply_n_struct_and_scalar

  !-----------------------------------------------------------------------------------------------------
  !> * (t_slm_struct_soil_pool * field3d)
  !!
  !-----------------------------------------------------------------------------------------------------
  FUNCTION Multiply_n_struct_and_field3d(this,field3d) RESULT(struct3)
    IMPLICIT NONE
    ! 0.1 InOut
    CLASS(t_slm_struct_soil_pool), INTENT(in) :: this
    REAL(wp),                        INTENT(in) :: field3d(:,:,:)
    TYPE(t_slm_struct_soil_pool), ALLOCATABLE :: struct3

    ALLOCATE(struct3, source=this)

    struct3%dom       = this%dom      * field3d
  END FUNCTION Multiply_n_struct_and_field3d

  !-----------------------------------------------------------------------------------------------------
  !> / (t_slm_struct_soil_pool / t_slm_struct_soil_pool)
  !!
  !-----------------------------------------------------------------------------------------------------
  FUNCTION Divide_n_struct(this,struct2) RESULT(struct3)
    IMPLICIT NONE
    ! 0.1 InOut
    CLASS(t_slm_struct_soil_pool), INTENT(in) :: this, struct2 
    TYPE(t_slm_struct_soil_pool), ALLOCATABLE :: struct3       

    ALLOCATE(struct3, source=this)

    struct3%dom       = this%dom      / struct2%dom
  END FUNCTION Divide_n_struct

  !-----------------------------------------------------------------------------------------------------
  !> / (t_slm_struct_soil_pool / scalar)
  !!
  !-----------------------------------------------------------------------------------------------------
  FUNCTION Divide_n_struct_and_scalar(this,scalar) RESULT(struct3)
    IMPLICIT NONE
    ! 0.1 InOut
    CLASS(t_slm_struct_soil_pool), INTENT(in) :: this
    REAL(wp),                        INTENT(in) :: scalar
    TYPE(t_slm_struct_soil_pool), ALLOCATABLE :: struct3

    ALLOCATE(struct3, source=this)

    struct3%dom       = this%dom      / scalar
  END FUNCTION Divide_n_struct_and_scalar

  !-----------------------------------------------------------------------------------------------------
  !> / (t_slm_struct_soil_pool / field3d)
  !!
  !-----------------------------------------------------------------------------------------------------
  FUNCTION Divide_n_struct_and_field3d(this,field3d) RESULT(struct3)
    IMPLICIT NONE
    ! 0.1 InOut
    CLASS(t_slm_struct_soil_pool), INTENT(in) :: this
    REAL(wp),                        INTENT(in) :: field3d(:,:,:)
    TYPE(t_slm_struct_soil_pool), ALLOCATABLE :: struct3

    ALLOCATE(struct3, source=this)

    struct3%dom       = this%dom      / field3d
  END FUNCTION Divide_n_struct_and_field3d

  !-----------------------------------------------------------------------------------------------------
  !> Initialisation t_slm_struct_soil_pool with a scalar
  !!
  !-----------------------------------------------------------------------------------------------------
  PURE ELEMENTAL SUBROUTINE init_t_slm_struct_soil_pool(this,s)
    IMPLICIT NONE
    ! 0.1 InOut
    CLASS(t_slm_struct_soil_pool),INTENT(inout) :: this 
    REAL(wp),                       INTENT(in)    :: s

    CALL this%dom%init(s) 
  END SUBROUTINE init_t_slm_struct_soil_pool 



  !-----------------------------------------------------------------------------------------------------
  !> initialize memory for the processes in: soil
  !!
  !!  
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE Init_soil_memory(mem, prefix, suffix, lct_ids, model_id)

!    USE mo_jsb_varlist,       ONLY: t_jsb_varlist, jsb_add_var
!    USE mo_jsb_io,            ONLY: grib_bits, t_cf, t_grib1, t_grib2, TSTEP_CONSTANT, tables
!    USE mo_jsb_grid_class,    ONLY: t_jsb_grid, t_jsb_vgrid
!    USE mo_jsb_grid,          ONLY: Get_grid, Get_vgrid
    USE mo_base_memory,         ONLY: t_model

    ! 0.1.InOut
    CLASS(t_soil_memory), INTENT(inout), TARGET :: mem             !> soil memory
    CHARACTER(len=*),     INTENT(in)            :: prefix          !> ?
    CHARACTER(len=*),     INTENT(in)            :: suffix          !> ?
    INTEGER,              INTENT(in)            :: lct_ids !(:)    !> in jsbach4: vector of land cover type IDs
    INTEGER,              INTENT(in)            :: model_id        !> model ID model\%id
    ! 0.2 Local
    TYPE(t_model), POINTER :: model
!    TYPE(t_jsb_model), POINTER :: model
!    TYPE(t_jsb_grid),  POINTER :: hgrid                        ! Horizontal grid
!    TYPE(t_jsb_vgrid), POINTER :: surface, canopy_layer        ! Vertical grid
!    INTEGER :: table

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':Init_soil_memory'

    model => Get_model(model_id)

!    table = tables(1)

!    hgrid   => Get_grid(mem%grid_id)
!    surface => Get_vgrid('surface')
!    canopy_layer => Get_vgrid('canopy_layer')


    IF (.NOT. ALLOCATED(mem%soil_depth%ptr)) THEN

      !> add variables to soil memory soil\%VAR
      !! 
      !! soil\%VAR
      !!  
      CALL slm_add_var(mem, 'soil_depth', mem%soil_depth, &
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'soil_depth_l', mem%soil_depth_l, &
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'root_depth', mem%root_depth, &
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'root_fraction_l', mem%root_fraction_l, &
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'runoff', mem%runoff, &
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'drainage', mem%drainage, &
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'drainage_sl', mem%drainage_sl, &
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'w_soil_root', mem%w_soil_root, &
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'w_soil_root_fc', mem%w_soil_root_fc, &
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'w_soil_root_pwp', mem%w_soil_root_pwp, &
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'w_soil_root_pot', mem%w_soil_root_pot, &
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'w_soil_root_theta', mem%w_soil_root_theta, &
            initval_r=1.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'w_soil_sl', mem%w_soil_sl, &
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'w_soil_fc_sl', mem%w_soil_fc_sl, &
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'w_soil_pwp_sl', mem%w_soil_pwp_sl, &
            initval_r=0.0_wp) ! initvalue

      CALL slm_add_var(mem, 'w_soil_pot_sl', mem%w_soil_pot_sl, &
            initval_r=0.0_wp) ! initvalue

      CALL slm_add_var(mem, 'nh4_solute', mem%nh4_solute, &
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'nh4_assoc', mem%nh4_assoc, &
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'no3_solute', mem%no3_solute, &
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'noy', mem%noy, &
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'n2o', mem%n2o, &
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'n2', mem%n2, &
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'nh4_n15_solute', mem%nh4_n15_solute, &
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'nh4_n15_assoc', mem%nh4_n15_assoc, &
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'no3_n15_solute', mem%no3_n15_solute, &
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'noy_n15', mem%noy_n15, &
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'n2o_n15', mem%n2o_n15, &
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'n2_n15', mem%n2_n15, &
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'po4_assoc', mem%po4_assoc, &
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'po4_solute', mem%po4_solute, &
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 't_soil', mem%t_soil, &
            initval_r=273.15_wp) ! initvalue

      CALL slm_add_var(mem, 'heat_capa_l', mem%heat_capa_l, &
            initval_r=0.0_wp) ! initvalue

      CALL slm_add_var(mem, 'bulk_dens_l', mem%bulk_dens_l, &
            initval_r=0.0_wp) ! initvalue

      CALL slm_add_var(mem, 'therm_cond_l', mem%therm_cond_l, &
            initval_r=0.0_wp) ! initvalue

      !> simple TYPE (t_jsb_struct)
      !! 
      !! slm_add_var -> slm_add_struct()
      !! 
!      CALL slm_add_var(mem, 'soil_test_soil_element',          mem%test_soil_element)


      !> nested TYPE (t_jsb_n_struct)
      !! 
      !! slm_add_var -> slm_add_n_struct()
      !! 
      CALL slm_add_var(mem, 'soil_pool'        , mem%pool)

      !> level 2 TYPE (l2_struct) nested within the above nested TYPE structures (t_jsb_n_struct)
      !! 
      !! slm_add_var -> slm_add_l2_struct()
      !! 
      CALL slm_add_var(mem%pool,    'soil_pool_dom',                mem%pool%dom) 
      CALL slm_add_var(mem%pool,    'soil_pool_dom_assoc',          mem%pool%dom_assoc)
      CALL slm_add_var(mem%pool,    'soil_pool_metabolic_litter',   mem%pool%metabolic_litter) 
      CALL slm_add_var(mem%pool,    'soil_pool_structural_litter',  mem%pool%structural_litter) 
      CALL slm_add_var(mem%pool,    'soil_pool_woody_litter',       mem%pool%woody_litter) 
      CALL slm_add_var(mem%pool,    'soil_pool_microbial',          mem%pool%microbial) 
      CALL slm_add_var(mem%pool,    'soil_pool_residue',            mem%pool%residue) 
      CALL slm_add_var(mem%pool,    'soil_pool_residue_assoc',      mem%pool%residue_assoc)  

      !> add variables to the structs and the l2 structs
      !!
      !! soil\%t_jsb_struct\%VAR
      !! 
      !! soil\%t_jsb_n_struct\%t_jsb_struct\%VAR
      !! 
      CALL slm_add_var(mem%pool%dom, 'soil_pool_dom_carbon', mem%pool%dom%carbon, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%pool%dom, 'soil_pool_dom_nitrogen', mem%pool%dom%nitrogen, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%pool%dom, 'soil_pool_dom_phosphorus', mem%pool%dom%phosphorus, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%pool%dom, 'soil_pool_dom_carbon13', mem%pool%dom%carbon13, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%pool%dom, 'soil_pool_dom_carbon14', mem%pool%dom%carbon14, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%pool%dom, 'soil_pool_dom_nitrogen15', mem%pool%dom%nitrogen15, &
            initval_r=0.0_wp)     ! initvalue
 
      CALL slm_add_var(mem%pool%dom_assoc, 'soil_pool_dom_assoc_carbon', mem%pool%dom_assoc%carbon, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%pool%dom_assoc, 'soil_pool_dom_assoc_nitrogen', mem%pool%dom_assoc%nitrogen, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%pool%dom_assoc, 'soil_pool_dom_assoc_phosphorus', mem%pool%dom_assoc%phosphorus, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%pool%dom_assoc, 'soil_pool_dom_assoc_carbon13', mem%pool%dom_assoc%carbon13, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%pool%dom_assoc, 'soil_pool_dom_assoc_carbon14', mem%pool%dom_assoc%carbon14, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%pool%dom_assoc, 'soil_pool_dom_assoc_nitrogen15', mem%pool%dom_assoc%nitrogen15, &
            initval_r=0.0_wp)     ! initvalue
 
      CALL slm_add_var(mem%pool%metabolic_litter, 'soil_pool_metabolic_litter_carbon', mem%pool%metabolic_litter%carbon, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%pool%metabolic_litter, 'soil_pool_metabolic_litter_nitrogen', mem%pool%metabolic_litter%nitrogen, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%pool%metabolic_litter, 'soil_pool_metabolic_litter_phosphorus', mem%pool%metabolic_litter%phosphorus, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%pool%metabolic_litter, 'soil_pool_metabolic_litter_carbon13', mem%pool%metabolic_litter%carbon13, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%pool%metabolic_litter, 'soil_pool_metabolic_litter_carbon14', mem%pool%metabolic_litter%carbon14, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%pool%metabolic_litter, 'soil_pool_metabolic_litter_nitrogen15', mem%pool%metabolic_litter%nitrogen15, &
            initval_r=0.0_wp)     ! initvalue
 
      CALL slm_add_var(mem%pool%structural_litter, 'soil_pool_structural_litter_carbon', mem%pool%structural_litter%carbon, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%pool%structural_litter, 'soil_pool_structural_litter_nitrogen', mem%pool%structural_litter%nitrogen, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%pool%structural_litter, 'soil_pool_structural_litter_phosphorus', mem%pool%structural_litter%phosphorus, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%pool%structural_litter, 'soil_pool_structural_litter_carbon13', mem%pool%structural_litter%carbon13, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%pool%structural_litter, 'soil_pool_structural_litter_carbon14', mem%pool%structural_litter%carbon14, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%pool%structural_litter, 'soil_pool_structural_litter_nitrogen15', mem%pool%structural_litter%nitrogen15, &
            initval_r=0.0_wp)     ! initvalue
 
      CALL slm_add_var(mem%pool%woody_litter, 'soil_pool_woody_litter_carbon', mem%pool%woody_litter%carbon, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%pool%woody_litter, 'soil_pool_woody_litter_nitrogen', mem%pool%woody_litter%nitrogen, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%pool%woody_litter, 'soil_pool_woody_litter_phosphorus', mem%pool%woody_litter%phosphorus, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%pool%woody_litter, 'soil_pool_woody_litter_carbon13', mem%pool%woody_litter%carbon13, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%pool%woody_litter, 'soil_pool_woody_litter_carbon14', mem%pool%woody_litter%carbon14, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%pool%woody_litter, 'soil_pool_woody_litter_nitrogen15', mem%pool%woody_litter%nitrogen15, &
            initval_r=0.0_wp)     ! initvalue
 
      CALL slm_add_var(mem%pool%microbial, 'soil_pool_microbial_carbon', mem%pool%microbial%carbon, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%pool%microbial, 'soil_pool_microbial_nitrogen', mem%pool%microbial%nitrogen, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%pool%microbial, 'soil_pool_microbial_phosphorus', mem%pool%microbial%phosphorus, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%pool%microbial, 'soil_pool_microbial_carbon13', mem%pool%microbial%carbon13, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%pool%microbial, 'soil_pool_microbial_carbon14', mem%pool%microbial%carbon14, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%pool%microbial, 'soil_pool_microbial_nitrogen15', mem%pool%microbial%nitrogen15, &
            initval_r=0.0_wp)     ! initvalue
 
      CALL slm_add_var(mem%pool%residue, 'soil_pool_residue_carbon', mem%pool%residue%carbon, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%pool%residue, 'soil_pool_residue_nitrogen', mem%pool%residue%nitrogen, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%pool%residue, 'soil_pool_residue_phosphorus', mem%pool%residue%phosphorus, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%pool%residue, 'soil_pool_residue_carbon13', mem%pool%residue%carbon13, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%pool%residue, 'soil_pool_residue_carbon14', mem%pool%residue%carbon14, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%pool%residue, 'soil_pool_residue_nitrogen15', mem%pool%residue%nitrogen15, &
            initval_r=0.0_wp)     ! initvalue
 
      CALL slm_add_var(mem%pool%residue_assoc, 'soil_pool_residue_assoc_carbon', mem%pool%residue_assoc%carbon, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%pool%residue_assoc, 'soil_pool_residue_assoc_nitrogen', mem%pool%residue_assoc%nitrogen, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%pool%residue_assoc, 'soil_pool_residue_assoc_phosphorus', mem%pool%residue_assoc%phosphorus, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%pool%residue_assoc, 'soil_pool_residue_assoc_carbon13', mem%pool%residue_assoc%carbon13, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%pool%residue_assoc, 'soil_pool_residue_assoc_carbon14', mem%pool%residue_assoc%carbon14, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%pool%residue_assoc, 'soil_pool_residue_assoc_nitrogen15', mem%pool%residue_assoc%nitrogen15, &
            initval_r=0.0_wp)     ! initvalue

    END IF

  END SUBROUTINE Init_soil_memory

END MODULE mo_soil_memory_class
