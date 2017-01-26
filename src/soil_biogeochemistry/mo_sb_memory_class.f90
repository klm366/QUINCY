!>
!! @brief Contains the memory class for the soil_biogeochemistry processes
!!
!! @author
!!  Jan Engel, MPI-BGC Jena (based on the JSBACH4 version from Rainer Schneck, MPI-M Hamburg)
!!
!! @par Revision History
!! First version by Jan Engel (20170103)
!!
MODULE mo_sb_memory_class

  USE mo_kind,                   ONLY: wp
  USE mo_exception,              ONLY: message
  USE mo_jsb_class,              ONLY: Get_model

  USE mo_jsb_memory_class,       ONLY: t_jsb_memory
  USE mo_jsb_var_class,          ONLY: t_jsb_var_real2d, t_jsb_var_real3d, &
                                       t_jsb_struct, t_jsb_n_struct             ! *slmjsb4

  USE mo_slm_varlist,            ONLY: slm_add_var   ! replace by mo_jsb_varlist later    *slmtemp



  IMPLICIT NONE
  PRIVATE
  PUBLIC :: t_sb_memory, max_no_of_vars
  PUBLIC :: t_slm_struct_sb_element, t_slm_struct_sb_pool   ! *slmjsb4
  PUBLIC :: max_no_of_structs, max_no_of_n_structs                      ! *slmjsb4

  INTEGER, PARAMETER :: max_no_of_vars      = 80   !> TYPE t_jsb_var
  INTEGER, PARAMETER :: max_no_of_structs   = 20   !> TYPE t_jsb_struct     *slmjsb4
  INTEGER, PARAMETER :: max_no_of_n_structs = 20   !> TYPE t_jsb_n_struct   *slmjsb4



  !-----------------------------------------------------------------------------------------------------
  !> Type definition for elements in sb pool
  !! 
  !! elements in soil are 3D 
  !! 
  !-----------------------------------------------------------------------------------------------------
  TYPE, EXTENDS(t_jsb_struct) :: t_slm_struct_sb_element
    TYPE(t_jsb_var_real3d)    :: carbon       !< total carbon concentration of pool [mol m-3]
    TYPE(t_jsb_var_real3d)    :: nitrogen     !< totol nitrogen concentration of pool [mol m-3]
    TYPE(t_jsb_var_real3d)    :: phosphorus   !< phosphorus concentration of pool [mol m-3]
    TYPE(t_jsb_var_real3d)    :: carbon13     !< 13carbon concentration of pool [mol m-3]
    TYPE(t_jsb_var_real3d)    :: carbon14     !< 14carbon concentration of pool [pmol m-3]
    TYPE(t_jsb_var_real3d)    :: nitrogen15   !< 15nitrogen concentration of pool [mol m-3]
  CONTAINS
    PRIVATE
!!$    !! operator overload
!!$    ! +
!!$    PROCEDURE               :: Add_struct                   ! Add two variables of type t_slm_struct_sb_element
!!$    PROCEDURE               :: Add_struct_and_scalar        ! Add variable of type t_slm_struct_sb_element with REAL scalar
!!$    PROCEDURE               :: Add_struct_and_field3d       ! Add variable of type t_slm_struct_sb_element with REAL 2-d field2d
!!$    GENERIC, PUBLIC         :: OPERATOR(+) => Add_struct, Add_struct_and_scalar, Add_struct_and_field3d
!!$    ! -
!!$    PROCEDURE               :: Subtract_struct              ! Subtract two variables of type t_slm_struct_sb_element
!!$    PROCEDURE               :: Subtract_struct_and_scalar   ! Subtract variable of type t_slm_struct_sb_element with REAL scalar
!!$    PROCEDURE               :: Subtract_struct_and_field3d  ! Subtract variable of type t_slm_struct_sb_element with REAL 2-d field2d
!!$    GENERIC, PUBLIC         :: OPERATOR(-) => Subtract_struct, Subtract_struct_and_scalar, Subtract_struct_and_field3d
!!$    ! *
!!$    PROCEDURE               :: Multiply_struct              ! Multiply two variables of type t_slm_struct_sb_element
!!$    PROCEDURE               :: Multiply_struct_and_scalar   ! Multiply variable of type t_slm_struct_sb_element with REAL scalar
!!$    PROCEDURE               :: Multiply_struct_and_field3d  ! Multiply variable of type t_slm_struct_sb_element with REAL 2-d field2d
!!$    GENERIC, PUBLIC         :: OPERATOR(*) => Multiply_struct, Multiply_struct_and_scalar, Multiply_struct_and_field3d
!!$    ! /
!!$    PROCEDURE               :: Divide_struct                ! Divide two variables of type t_slm_struct_sb_element
!!$    PROCEDURE               :: Divide_struct_and_scalar     ! Divide variable of type t_slm_struct_sb_element by REAL scalar
!!$    PROCEDURE               :: Divide_struct_and_field3d    ! Divide variable of type t_slm_struct_sb_element by REAL 2-d field2d
!!$    GENERIC, PUBLIC         :: OPERATOR(/) => Divide_struct, Divide_struct_and_scalar, Divide_struct_and_field3d
    !! initialize with/ set to a certain REAL
    PROCEDURE, PUBLIC       :: init => init_t_slm_struct_sb_element
  END TYPE t_slm_struct_sb_element


  !-----------------------------------------------------------------------------------------------------
  !> Type definition for sb pool with nested TYPEs
  !!
  !! used for: pool ...
  !!  
  !-----------------------------------------------------------------------------------------------------
  TYPE, EXTENDS(t_jsb_n_struct) :: t_slm_struct_sb_pool
    TYPE(t_slm_struct_sb_element)  :: dom                 !< dissolved organic matter
    TYPE(t_slm_struct_sb_element)  :: dom_assoc           !< minerally associated disolved organic matter
    TYPE(t_slm_struct_sb_element)  :: soluable_litter    !< metabolic litter (sugars, cellulose etc.)
    TYPE(t_slm_struct_sb_element)  :: polymeric_litter   !< structural litter (lignified litter)
    TYPE(t_slm_struct_sb_element)  :: woody_litter        !< physically protected wood litter
    TYPE(t_slm_struct_sb_element)  :: microbial           !< microbial biomass
    TYPE(t_slm_struct_sb_element)  :: residue             !< microbial necromass
    TYPE(t_slm_struct_sb_element)  :: residue_assoc       !< minerally associated microbial necromass 
  CONTAINS
    PRIVATE
!!$    !! operator overload
!!$    ! +
!!$    PROCEDURE               :: Add_n_struct                   
!!$    PROCEDURE               :: Add_n_struct_and_scalar        
!!$    PROCEDURE               :: Add_n_struct_and_field3d       
!!$    GENERIC, PUBLIC         :: OPERATOR(+) => Add_n_struct, Add_n_struct_and_scalar, Add_n_struct_and_field3d
!!$    ! -
!!$    PROCEDURE               :: Subtract_n_struct              
!!$    PROCEDURE               :: Subtract_n_struct_and_scalar   
!!$    PROCEDURE               :: Subtract_n_struct_and_field3d  
!!$    GENERIC, PUBLIC         :: OPERATOR(-) => Subtract_n_struct, Subtract_n_struct_and_scalar, Subtract_n_struct_and_field3d
!!$    ! *
!!$    PROCEDURE               :: Multiply_n_struct              
!!$    PROCEDURE               :: Multiply_n_struct_and_scalar   
!!$    PROCEDURE               :: Multiply_n_struct_and_field3d  
!!$    GENERIC, PUBLIC         :: OPERATOR(*) => Multiply_n_struct, Multiply_n_struct_and_scalar, Multiply_n_struct_and_field3d
!!$    ! /
!!$    PROCEDURE               :: Divide_n_struct                
!!$    PROCEDURE               :: Divide_n_struct_and_scalar     
!!$    PROCEDURE               :: Divide_n_struct_and_field3d    
!!$    GENERIC, PUBLIC         :: OPERATOR(/) => Divide_n_struct, Divide_n_struct_and_scalar, Divide_n_struct_and_field3d
    !! initialize with/ set to a certain REAL
    PROCEDURE, PUBLIC       :: init => init_t_slm_struct_sb_pool
  END TYPE t_slm_struct_sb_pool


  !-----------------------------------------------------------------------------------------------------
  !> Type definition for sb memory
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
  TYPE, EXTENDS(t_jsb_memory) :: t_sb_memory

    ! pools
    TYPE(t_jsb_var_real3d)          :: &
                                    nh4_solute         , &  !< ammonium in soil solution [mol m-3]
                                    nh4_assoc          , &  !< minerally associated ammonium [mol m-3]
                                    no3_solute         , &  !< nitrate in soil solution [mol m-3]
                                    noy                , &  !< NO_y [mol m-3]
                                    n2o                , &  !< N2O [mol m-3]
                                    n2                 , &  !< N2 [mol m-3]
                                    nh4_n15_solute     , &  !< N15-ammonium in soil solution [mol m-3]
                                    nh4_n15_assoc      , &  !< minerally associated N15-ammonium [mol m-3]
                                    no3_n15_solute     , &  !< N15-nitrate in soil solution [mol m-3]
                                    noy_n15            , &  !< N15-NO_y [mol m-3]
                                    n2o_n15            , &  !< N15-N2O [mol m-3]
                                    n2_n15             , &  !< N15-N2 [mol m-3]
                                    po4_assoc          , &  !< phosphate in soil solution [mol m-3]
                                    po4_solute              !< minerally associated phosphate  [mol m-3]

    ! transport fluxes 
    TYPE(t_jsb_var_real3d)          :: &
                                    transport_nh4_solute         , &  !< transport of ammonium in soil solution [mol m-3 timestep-1]
                                    transport_nh4_assoc          , &  !< transport of minerally associated ammonium [mol m-3 timestep-1]
                                    transport_no3_solute         , &  !< transport of nitrate in soil solution [mol m-3 timestep-1]
                                    transport_noy                , &  !< transport of NO_y [mol m-3 timestep-1]
                                    transport_n2o                , &  !< transport of N2O [mol m-3 timestep-1]
                                    transport_n2                 , &  !< transport of N2 [mol m-3 timestep-1]
                                    transport_nh4_n15_solute     , &  !< transport of N15-ammonium in soil solution [mol m-3 timestep-1]
                                    transport_nh4_n15_assoc      , &  !< transport of minerally associated N15-ammonium [mol m-3 timestep-1]
                                    transport_no3_n15_solute     , &  !< transport of N15-nitrate in soil solution [mol m-3 timestep-1]
                                    transport_noy_n15            , &  !< transport of N15-NO_y [mol m-3 timestep-1]
                                    transport_n2o_n15            , &  !< transport of N15-N2O [mol m-3 timestep-1]
                                    transport_n2_n15             , &  !< transport of N15-N2 [mol m-3 timestep-1]
                                    transport_po4_assoc          , &  !< transport of phosphate in soil solution [mol m-3 timestep-1]
                                    transport_po4_solute              !< transport of minerally associated phosphate  [mol m-3 timestep-1]

    ! temperature response functions
    TYPE(t_jsb_var_real3d)          :: &
                                    rtm_depolymerisation, & !< rate modifier for depolymerisation processes [--]
                                    rtm_mic_uptake, &       !< rate modifier for microbial uptake [--] 
                                    rtm_sorption, &         !< rate modifier for sorption [--]
                                    rtm_desorption, &       !< rate modifier for desorption [--]
                                    rtm_hsc                 !< rate modifier for half-saturation constants [--] 

    ! moisture response functions
    TYPE(t_jsb_var_real3d)          :: &
                                    rmm_depolymerisation, & !< rate modifier for depolymerisation processes [--]
                                    rmm_mic_uptake, &       !< rate modifier for microbial uptake [--] 
                                    rmm_sorption, &         !< rate modifier for sorption [--]
                                    rmm_desorption, &       !< rate modifier for desorption [--]
                                    rmm_hsc                 !< rate modifier for half-saturation constants [--] 

    TYPE(t_jsb_var_real3d)          :: &
                                    het_respiration    ,    &  !< heterotrophic respiration  [mol m-3 s-1]
                                    het_respiration_c13,    &  !< c13 heterotrophic respiration   [mol m-3 s-1]
                                    het_respiration_c14,    &  !< c14 heterotrophic respiration   [mol m-3 s-1]
                                    net_mineralisation_n,   &  !< net mineralisation of N    [mol m-3 s-1]
                                    net_mineralisation_n15, &  !< net mineralisation of N15  [mol m-3 s-1] 
                                    net_mineralisation_p,   &  !< net mineralisation of P  [mol m-3 s-1]
                                    volatilisation_n,       &  !< volatilisation of N  [µmol m-3 s-1]
                                    nitrification              !< nitrification  [µmol m-3 s-1]
    
    TYPE(t_jsb_var_real3d)          :: &
      qmax_mineral_org, &        !< sorption capacity of mineral soil [kg C kg-1 mineral soil]
      k_bioturb               !< diffusion rate resulting from bioturbation 

!    TYPE(t_slm_struct_sb_element) :: &   ! simple structure for elements: t_jsb_struct
!                                    test_sb_element   !< test ...


    TYPE(t_slm_struct_sb_pool)    :: &   ! nested pool TYPE t_jsb_n_struct
                                    pool, &                !< pools of the sb
                                    formation, & 
                                    loss, & 
                                    transport



  CONTAINS
    PROCEDURE :: Init => Init_sb_memory
  END TYPE t_sb_memory


  CHARACTER(len=*), PARAMETER :: modname = 'mo_sb_memory_class'


CONTAINS

!!$  !-----------------------------------------------------------------------------------------------------!
!!$  ! t_slm_struct_sb_element operator overload                                                           !
!!$  !-----------------------------------------------------------------------------------------------------!
!!$
!!$  !-----------------------------------------------------------------------------------------------------
!!$  !> + (t_slm_struct_sb_element + t_slm_struct_sb_element)
!!$  !!
!!$  !-----------------------------------------------------------------------------------------------------
!!$  FUNCTION Add_struct(this,struct2) RESULT(struct3)
!!$    IMPLICIT NONE
!!$    ! 0.1 InOut
!!$    CLASS(t_slm_struct_sb_element), INTENT(in) :: this, struct2  !< t_slm_struct_sb_element in
!!$    TYPE(t_slm_struct_sb_element), ALLOCATABLE :: struct3        !< t_slm_struct_sb_element out
!!$
!!$    ALLOCATE(struct3, source=this)
!!$
!!$    struct3%carbon       = this%carbon      + struct2%carbon
!!$    struct3%nitrogen     = this%nitrogen    + struct2%nitrogen
!!$    struct3%phosphorus   = this%phosphorus  + struct2%phosphorus
!!$    struct3%carbon13     = this%carbon13    + struct2%carbon13
!!$    struct3%carbon14     = this%carbon14    + struct2%carbon14
!!$    struct3%nitrogen15   = this%nitrogen15  + struct2%nitrogen15
!!$  END FUNCTION Add_struct
!!$
!!$  !-----------------------------------------------------------------------------------------------------
!!$  !> + (t_slm_struct_sb_element + scalar)
!!$  !!
!!$  !-----------------------------------------------------------------------------------------------------
!!$  FUNCTION Add_struct_and_scalar(this,scalar) RESULT(struct3)
!!$    IMPLICIT NONE
!!$    ! 0.1 InOut
!!$    CLASS(t_slm_struct_sb_element), INTENT(in) :: this
!!$    REAL(wp),                        INTENT(in) :: scalar
!!$    TYPE(t_slm_struct_sb_element), ALLOCATABLE :: struct3
!!$
!!$    ALLOCATE(struct3, source=this)
!!$
!!$    struct3%carbon       = this%carbon      + scalar
!!$    struct3%nitrogen     = this%nitrogen    + scalar
!!$    struct3%phosphorus   = this%phosphorus  + scalar
!!$    struct3%carbon13     = this%carbon13    + scalar
!!$    struct3%carbon14     = this%carbon14    + scalar
!!$    struct3%nitrogen15   = this%nitrogen15  + scalar
!!$  END FUNCTION Add_struct_and_scalar
!!$
!!$  !-----------------------------------------------------------------------------------------------------
!!$  !> + (t_slm_struct_sb_element + field3d)
!!$  !!
!!$  !-----------------------------------------------------------------------------------------------------
!!$  FUNCTION Add_struct_and_field3d(this,field3d) RESULT(struct3)
!!$    IMPLICIT NONE
!!$    ! 0.1 InOut
!!$    CLASS(t_slm_struct_sb_element), INTENT(in) :: this
!!$    REAL(wp),                        INTENT(in) :: field3d(:,:,:)
!!$    TYPE(t_slm_struct_sb_element), ALLOCATABLE :: struct3
!!$
!!$    ALLOCATE(struct3, source=this)
!!$
!!$    struct3%carbon       = this%carbon      + field3d
!!$    struct3%nitrogen     = this%nitrogen    + field3d
!!$    struct3%phosphorus   = this%phosphorus  + field3d
!!$    struct3%carbon13     = this%carbon13    + field3d
!!$    struct3%carbon14     = this%carbon14    + field3d
!!$    struct3%nitrogen15   = this%nitrogen15  + field3d
!!$  END FUNCTION Add_struct_and_field3d
!!$
!!$  !-----------------------------------------------------------------------------------------------------
!!$  !> - (t_slm_struct_sb_element - t_slm_struct_sb_element)
!!$  !!
!!$  !-----------------------------------------------------------------------------------------------------
!!$  FUNCTION Subtract_struct(this,struct2) RESULT(struct3)
!!$    IMPLICIT NONE
!!$    ! 0.1 InOut
!!$    CLASS(t_slm_struct_sb_element), INTENT(in) :: this, struct2  !< t_slm_struct_sb_element in
!!$    TYPE(t_slm_struct_sb_element), ALLOCATABLE :: struct3        !< t_slm_struct_sb_element out
!!$
!!$    ALLOCATE(struct3, source=this)
!!$
!!$    struct3%carbon       = this%carbon      - struct2%carbon
!!$    struct3%nitrogen     = this%nitrogen    - struct2%nitrogen
!!$    struct3%phosphorus   = this%phosphorus  - struct2%phosphorus
!!$    struct3%carbon13     = this%carbon13    - struct2%carbon13
!!$    struct3%carbon14     = this%carbon14    - struct2%carbon14
!!$    struct3%nitrogen15   = this%nitrogen15  - struct2%nitrogen15
!!$  END FUNCTION Subtract_struct
!!$
!!$  !-----------------------------------------------------------------------------------------------------
!!$  !> - (t_slm_struct_sb_element - scalar)
!!$  !!
!!$  !-----------------------------------------------------------------------------------------------------
!!$  FUNCTION Subtract_struct_and_scalar(this,scalar) RESULT(struct3)
!!$    IMPLICIT NONE
!!$    ! 0.1 InOut
!!$    CLASS(t_slm_struct_sb_element), INTENT(in) :: this
!!$    REAL(wp),                        INTENT(in) :: scalar
!!$    TYPE(t_slm_struct_sb_element), ALLOCATABLE :: struct3
!!$
!!$    ALLOCATE(struct3, source=this)
!!$
!!$    struct3%carbon       = this%carbon      - scalar
!!$    struct3%nitrogen     = this%nitrogen    - scalar
!!$    struct3%phosphorus   = this%phosphorus  - scalar
!!$    struct3%carbon13     = this%carbon13    - scalar
!!$    struct3%carbon14     = this%carbon14    - scalar
!!$    struct3%nitrogen15   = this%nitrogen15  - scalar
!!$  END FUNCTION Subtract_struct_and_scalar
!!$
!!$  !-----------------------------------------------------------------------------------------------------
!!$  !> - (t_slm_struct_sb_element - field3d)
!!$  !!
!!$  !-----------------------------------------------------------------------------------------------------
!!$  FUNCTION Subtract_struct_and_field3d(this,field3d) RESULT(struct3)
!!$    IMPLICIT NONE
!!$    ! 0.1 InOut
!!$    CLASS(t_slm_struct_sb_element), INTENT(in) :: this
!!$    REAL(wp),                        INTENT(in) :: field3d(:,:,:)
!!$    TYPE(t_slm_struct_sb_element), ALLOCATABLE :: struct3
!!$
!!$    ALLOCATE(struct3, source=this)
!!$
!!$    struct3%carbon       = this%carbon      - field3d
!!$    struct3%nitrogen     = this%nitrogen    - field3d
!!$    struct3%phosphorus   = this%phosphorus  - field3d
!!$    struct3%carbon13     = this%carbon13    - field3d
!!$    struct3%carbon14     = this%carbon14    - field3d
!!$    struct3%nitrogen15   = this%nitrogen15  - field3d
!!$  END FUNCTION Subtract_struct_and_field3d
!!$
!!$  !-----------------------------------------------------------------------------------------------------
!!$  !> * (t_slm_struct_sb_element * t_slm_struct_sb_element)
!!$  !!
!!$  !-----------------------------------------------------------------------------------------------------
!!$  FUNCTION Multiply_struct(this,struct2) RESULT(struct3)
!!$    IMPLICIT NONE
!!$    ! 0.1 InOut
!!$    CLASS(t_slm_struct_sb_element), INTENT(in) :: this, struct2  !< t_slm_struct_sb_element in
!!$    TYPE(t_slm_struct_sb_element), ALLOCATABLE :: struct3        !< t_slm_struct_sb_element out
!!$
!!$    ALLOCATE(struct3, source=this)
!!$
!!$    struct3%carbon       = this%carbon      * struct2%carbon
!!$    struct3%nitrogen     = this%nitrogen    * struct2%nitrogen
!!$    struct3%phosphorus   = this%phosphorus  * struct2%phosphorus
!!$    struct3%carbon13     = this%carbon13    * struct2%carbon13
!!$    struct3%carbon14     = this%carbon14    * struct2%carbon14
!!$    struct3%nitrogen15   = this%nitrogen15  * struct2%nitrogen15
!!$  END FUNCTION Multiply_struct
!!$
!!$  !-----------------------------------------------------------------------------------------------------
!!$  !> * (t_slm_struct_sb_element * scalar)
!!$  !!
!!$  !-----------------------------------------------------------------------------------------------------
!!$  FUNCTION Multiply_struct_and_scalar(this,scalar) RESULT(struct3)
!!$    IMPLICIT NONE
!!$    ! 0.1 InOut
!!$    CLASS(t_slm_struct_sb_element), INTENT(in) :: this
!!$    REAL(wp),                        INTENT(in) :: scalar
!!$    TYPE(t_slm_struct_sb_element), ALLOCATABLE :: struct3
!!$
!!$    ALLOCATE(struct3, source=this)
!!$
!!$    struct3%carbon       = this%carbon      * scalar
!!$    struct3%nitrogen     = this%nitrogen    * scalar
!!$    struct3%phosphorus   = this%phosphorus  * scalar
!!$    struct3%carbon13     = this%carbon13    * scalar
!!$    struct3%carbon14     = this%carbon14    * scalar
!!$    struct3%nitrogen15   = this%nitrogen15  * scalar
!!$  END FUNCTION Multiply_struct_and_scalar
!!$
!!$  !-----------------------------------------------------------------------------------------------------
!!$  !> * (t_slm_struct_sb_element * field3d)
!!$  !!
!!$  !-----------------------------------------------------------------------------------------------------
!!$  FUNCTION Multiply_struct_and_field3d(this,field3d) RESULT(struct3)
!!$    IMPLICIT NONE
!!$    ! 0.1 InOut
!!$    CLASS(t_slm_struct_sb_element), INTENT(in) :: this
!!$    REAL(wp),                        INTENT(in) :: field3d(:,:,:)
!!$    TYPE(t_slm_struct_sb_element), ALLOCATABLE :: struct3
!!$
!!$    ALLOCATE(struct3, source=this)
!!$
!!$    struct3%carbon       = this%carbon      * field3d
!!$    struct3%nitrogen     = this%nitrogen    * field3d
!!$    struct3%phosphorus   = this%phosphorus  * field3d
!!$    struct3%carbon13     = this%carbon13    * field3d
!!$    struct3%carbon14     = this%carbon14    * field3d
!!$    struct3%nitrogen15   = this%nitrogen15  * field3d
!!$  END FUNCTION Multiply_struct_and_field3d
!!$
!!$  !-----------------------------------------------------------------------------------------------------
!!$  !> / (t_slm_struct_sb_element / t_slm_struct_sb_element)
!!$  !!
!!$  !-----------------------------------------------------------------------------------------------------
!!$  FUNCTION Divide_struct(this,struct2) RESULT(struct3)
!!$    IMPLICIT NONE
!!$    ! 0.1 InOut
!!$    CLASS(t_slm_struct_sb_element), INTENT(in) :: this, struct2 
!!$    TYPE(t_slm_struct_sb_element), ALLOCATABLE :: struct3       
!!$
!!$    ALLOCATE(struct3, source=this)
!!$
!!$    struct3%carbon       = this%carbon      / struct2%carbon
!!$    struct3%nitrogen     = this%nitrogen    / struct2%nitrogen
!!$    struct3%phosphorus   = this%phosphorus  / struct2%phosphorus
!!$    struct3%carbon13     = this%carbon13    / struct2%carbon13
!!$    struct3%carbon14     = this%carbon14    / struct2%carbon14
!!$    struct3%nitrogen15   = this%nitrogen15  / struct2%nitrogen15
!!$  END FUNCTION Divide_struct
!!$
!!$  !-----------------------------------------------------------------------------------------------------
!!$  !> / (t_slm_struct_sb_element / scalar)
!!$  !!
!!$  !-----------------------------------------------------------------------------------------------------
!!$  FUNCTION Divide_struct_and_scalar(this,scalar) RESULT(struct3)
!!$    IMPLICIT NONE
!!$    ! 0.1 InOut
!!$    CLASS(t_slm_struct_sb_element), INTENT(in) :: this
!!$    REAL(wp),                        INTENT(in) :: scalar
!!$    TYPE(t_slm_struct_sb_element), ALLOCATABLE :: struct3
!!$
!!$    ALLOCATE(struct3, source=this)
!!$
!!$    struct3%carbon       = this%carbon      / scalar
!!$    struct3%nitrogen     = this%nitrogen    / scalar
!!$    struct3%phosphorus   = this%phosphorus  / scalar
!!$    struct3%carbon13     = this%carbon13    / scalar
!!$    struct3%carbon14     = this%carbon14    / scalar
!!$    struct3%nitrogen15   = this%nitrogen15  / scalar
!!$  END FUNCTION Divide_struct_and_scalar
!!$
!!$  !-----------------------------------------------------------------------------------------------------
!!$  !> / (t_slm_struct_sb_element / field3d)
!!$  !!
!!$  !-----------------------------------------------------------------------------------------------------
!!$  FUNCTION Divide_struct_and_field3d(this,field3d) RESULT(struct3)
!!$    IMPLICIT NONE
!!$    ! 0.1 InOut
!!$    CLASS(t_slm_struct_sb_element), INTENT(in) :: this
!!$    REAL(wp),                        INTENT(in) :: field3d(:,:,:)
!!$    TYPE(t_slm_struct_sb_element), ALLOCATABLE :: struct3
!!$
!!$    ALLOCATE(struct3, source=this)
!!$
!!$    struct3%carbon       = this%carbon      / field3d
!!$    struct3%nitrogen     = this%nitrogen    / field3d
!!$    struct3%phosphorus   = this%phosphorus  / field3d
!!$    struct3%carbon13     = this%carbon13    / field3d
!!$    struct3%carbon14     = this%carbon14    / field3d
!!$    struct3%nitrogen15   = this%nitrogen15  / field3d
!!$  END FUNCTION Divide_struct_and_field3d
!!$
  !-----------------------------------------------------------------------------------------------------
  !> Initialisation t_slm_struct_sb_element with a scalar
  !!
  !-----------------------------------------------------------------------------------------------------
  PURE ELEMENTAL SUBROUTINE init_t_slm_struct_sb_element(this,s)
    IMPLICIT NONE
    ! 0.1 InOut
    CLASS(t_slm_struct_sb_element),INTENT(inout) :: this 
    REAL(wp),                       INTENT(in)    :: s

    CALL this%carbon%init(s) 
    CALL this%nitrogen%init(s) 
    CALL this%phosphorus%init(s) 
    CALL this%carbon13%init(s) 
    CALL this%carbon14%init(s) 
    CALL this%nitrogen15%init(s) 
  END SUBROUTINE init_t_slm_struct_sb_element 



!!$  !-----------------------------------------------------------------------------------------------------!
!!$  ! t_slm_struct_sb_pool operator overload                                                              !
!!$  !-----------------------------------------------------------------------------------------------------!
!!$
!!$  !-----------------------------------------------------------------------------------------------------
!!$  !> + (t_slm_struct_sb_pool + t_slm_struct_sb_pool)
!!$  !!
!!$  !-----------------------------------------------------------------------------------------------------
!!$  FUNCTION Add_n_struct(this,struct2) RESULT(struct3)
!!$    IMPLICIT NONE
!!$    ! 0.1 InOut
!!$    CLASS(t_slm_struct_sb_pool), INTENT(in) :: this, struct2  !< t_slm_struct_sb_pool in
!!$    TYPE(t_slm_struct_sb_pool), ALLOCATABLE :: struct3        !< t_slm_struct_sb_pool out
!!$
!!$    ALLOCATE(struct3, source=this)
!!$
!!$    struct3%dom               = this%dom                + struct2%dom
!!$    struct3%dom_assoc         = this%dom_assoc          + struct2%dom_assoc
!!$    struct3%soluable_litter  = this%soluable_litter   + struct2%soluable_litter
!!$    struct3%polymeric_litter = this%polymeric_litter  + struct2%polymeric_litter
!!$    struct3%woody_litter      = this%woody_litter       + struct2%woody_litter
!!$    struct3%microbial         = this%microbial          + struct2%microbial
!!$    struct3%residue           = this%residue            + struct2%residue
!!$    struct3%residue_assoc     = this%residue_assoc      + struct2%residue_assoc
!!$  END FUNCTION Add_n_struct
!!$
!!$  !-----------------------------------------------------------------------------------------------------
!!$  !> + (t_slm_struct_sb_pool + scalar)
!!$  !!
!!$  !-----------------------------------------------------------------------------------------------------
!!$  FUNCTION Add_n_struct_and_scalar(this,scalar) RESULT(struct3)
!!$    IMPLICIT NONE
!!$    ! 0.1 InOut
!!$    CLASS(t_slm_struct_sb_pool), INTENT(in) :: this
!!$    REAL(wp),                        INTENT(in) :: scalar
!!$    TYPE(t_slm_struct_sb_pool), ALLOCATABLE :: struct3
!!$
!!$    ALLOCATE(struct3, source=this)
!!$
!!$    struct3%dom               = this%dom                + scalar
!!$    struct3%dom_assoc         = this%dom_assoc          + scalar
!!$    struct3%soluable_litter  = this%soluable_litter   + scalar
!!$    struct3%polymeric_litter = this%polymeric_litter  + scalar
!!$    struct3%woody_litter      = this%woody_litter       + scalar
!!$    struct3%microbial         = this%microbial          + scalar
!!$    struct3%residue           = this%residue            + scalar
!!$    struct3%residue_assoc     = this%residue_assoc      + scalar
!!$  END FUNCTION Add_n_struct_and_scalar
!!$
!!$  !-----------------------------------------------------------------------------------------------------
!!$  !> + (t_slm_struct_sb_pool + field3d)
!!$  !!
!!$  !-----------------------------------------------------------------------------------------------------
!!$  FUNCTION Add_n_struct_and_field3d(this,field3d) RESULT(struct3)
!!$    IMPLICIT NONE
!!$    ! 0.1 InOut
!!$    CLASS(t_slm_struct_sb_pool), INTENT(in) :: this
!!$    REAL(wp),                        INTENT(in) :: field3d(:,:,:)
!!$    TYPE(t_slm_struct_sb_pool), ALLOCATABLE :: struct3
!!$
!!$    ALLOCATE(struct3, source=this)
!!$
!!$    struct3%dom               = this%dom                + field3d
!!$    struct3%dom_assoc         = this%dom_assoc          + field3d
!!$    struct3%soluable_litter  = this%soluable_litter   + field3d
!!$    struct3%polymeric_litter = this%polymeric_litter  + field3d
!!$    struct3%woody_litter      = this%woody_litter       + field3d
!!$    struct3%microbial         = this%microbial          + field3d
!!$    struct3%residue           = this%residue            + field3d
!!$    struct3%residue_assoc     = this%residue_assoc      + field3d
!!$  END FUNCTION Add_n_struct_and_field3d
!!$
!!$  !-----------------------------------------------------------------------------------------------------
!!$  !> - (t_slm_struct_sb_pool - t_slm_struct_sb_pool)
!!$  !!
!!$  !-----------------------------------------------------------------------------------------------------
!!$  FUNCTION Subtract_n_struct(this,struct2) RESULT(struct3)
!!$    IMPLICIT NONE
!!$    ! 0.1 InOut
!!$    CLASS(t_slm_struct_sb_pool), INTENT(in) :: this, struct2  !< t_slm_struct_sb_pool in
!!$    TYPE(t_slm_struct_sb_pool), ALLOCATABLE :: struct3        !< t_slm_struct_sb_pool out
!!$
!!$    ALLOCATE(struct3, source=this)
!!$
!!$    struct3%dom               = this%dom                - struct2%dom
!!$    struct3%dom_assoc         = this%dom_assoc          - struct2%dom_assoc
!!$    struct3%soluable_litter  = this%soluable_litter   - struct2%soluable_litter
!!$    struct3%polymeric_litter = this%polymeric_litter  - struct2%polymeric_litter
!!$    struct3%woody_litter      = this%woody_litter       - struct2%woody_litter
!!$    struct3%microbial         = this%microbial          - struct2%microbial
!!$    struct3%residue           = this%residue            - struct2%residue
!!$    struct3%residue_assoc     = this%residue_assoc      - struct2%residue_assoc
!!$  END FUNCTION Subtract_n_struct
!!$
!!$  !-----------------------------------------------------------------------------------------------------
!!$  !> - (t_slm_struct_sb_pool - scalar)
!!$  !!
!!$  !-----------------------------------------------------------------------------------------------------
!!$  FUNCTION Subtract_n_struct_and_scalar(this,scalar) RESULT(struct3)
!!$    IMPLICIT NONE
!!$    ! 0.1 InOut
!!$    CLASS(t_slm_struct_sb_pool), INTENT(in) :: this
!!$    REAL(wp),                        INTENT(in) :: scalar
!!$    TYPE(t_slm_struct_sb_pool), ALLOCATABLE :: struct3
!!$
!!$    ALLOCATE(struct3, source=this)
!!$
!!$    struct3%dom               = this%dom                - scalar
!!$    struct3%dom_assoc         = this%dom_assoc          - scalar
!!$    struct3%soluable_litter  = this%soluable_litter   - scalar
!!$    struct3%polymeric_litter = this%polymeric_litter  - scalar
!!$    struct3%woody_litter      = this%woody_litter       - scalar
!!$    struct3%microbial         = this%microbial          - scalar
!!$    struct3%residue           = this%residue            - scalar
!!$    struct3%residue_assoc     = this%residue_assoc      - scalar
!!$  END FUNCTION Subtract_n_struct_and_scalar
!!$
!!$  !-----------------------------------------------------------------------------------------------------
!!$  !> - (t_slm_struct_sb_pool - field3d)
!!$  !!
!!$  !-----------------------------------------------------------------------------------------------------
!!$  FUNCTION Subtract_n_struct_and_field3d(this,field3d) RESULT(struct3)
!!$    IMPLICIT NONE
!!$    ! 0.1 InOut
!!$    CLASS(t_slm_struct_sb_pool), INTENT(in) :: this
!!$    REAL(wp),                        INTENT(in) :: field3d(:,:,:)
!!$    TYPE(t_slm_struct_sb_pool), ALLOCATABLE :: struct3
!!$
!!$    ALLOCATE(struct3, source=this)
!!$
!!$    struct3%dom               = this%dom                - field3d
!!$    struct3%dom_assoc         = this%dom_assoc          - field3d
!!$    struct3%soluable_litter  = this%soluable_litter   - field3d
!!$    struct3%polymeric_litter = this%polymeric_litter  - field3d
!!$    struct3%woody_litter      = this%woody_litter       - field3d
!!$    struct3%microbial         = this%microbial          - field3d
!!$    struct3%residue           = this%residue            - field3d
!!$    struct3%residue_assoc     = this%residue_assoc      - field3d
!!$  END FUNCTION Subtract_n_struct_and_field3d
!!$
!!$  !-----------------------------------------------------------------------------------------------------
!!$  !> * (t_slm_struct_sb_pool * t_slm_struct_sb_pool)
!!$  !!
!!$  !-----------------------------------------------------------------------------------------------------
!!$  FUNCTION Multiply_n_struct(this,struct2) RESULT(struct3)
!!$    IMPLICIT NONE
!!$    ! 0.1 InOut
!!$    CLASS(t_slm_struct_sb_pool), INTENT(in) :: this, struct2  !< t_slm_struct_sb_pool in
!!$    TYPE(t_slm_struct_sb_pool), ALLOCATABLE :: struct3        !< t_slm_struct_sb_pool out
!!$
!!$    ALLOCATE(struct3, source=this)
!!$
!!$    struct3%dom               = this%dom                * struct2%dom
!!$    struct3%dom_assoc         = this%dom_assoc          * struct2%dom_assoc
!!$    struct3%soluable_litter  = this%soluable_litter   * struct2%soluable_litter
!!$    struct3%polymeric_litter = this%polymeric_litter  * struct2%polymeric_litter
!!$    struct3%woody_litter      = this%woody_litter       * struct2%woody_litter
!!$    struct3%microbial         = this%microbial          * struct2%microbial
!!$    struct3%residue           = this%residue            * struct2%residue
!!$    struct3%residue_assoc     = this%residue_assoc      * struct2%residue_assoc
!!$  END FUNCTION Multiply_n_struct
!!$
!!$  !-----------------------------------------------------------------------------------------------------
!!$  !> * (t_slm_struct_sb_pool * scalar)
!!$  !!
!!$  !-----------------------------------------------------------------------------------------------------
!!$  FUNCTION Multiply_n_struct_and_scalar(this,scalar) RESULT(struct3)
!!$    IMPLICIT NONE
!!$    ! 0.1 InOut
!!$    CLASS(t_slm_struct_sb_pool), INTENT(in) :: this
!!$    REAL(wp),                        INTENT(in) :: scalar
!!$    TYPE(t_slm_struct_sb_pool), ALLOCATABLE :: struct3
!!$
!!$    ALLOCATE(struct3, source=this)
!!$
!!$    struct3%dom               = this%dom                * scalar
!!$    struct3%dom_assoc         = this%dom_assoc          * scalar
!!$    struct3%soluable_litter  = this%soluable_litter   * scalar
!!$    struct3%polymeric_litter = this%polymeric_litter  * scalar
!!$    struct3%woody_litter      = this%woody_litter       * scalar
!!$    struct3%microbial         = this%microbial          * scalar
!!$    struct3%residue           = this%residue            * scalar
!!$    struct3%residue_assoc     = this%residue_assoc      * scalar
!!$  END FUNCTION Multiply_n_struct_and_scalar
!!$
!!$  !-----------------------------------------------------------------------------------------------------
!!$  !> * (t_slm_struct_sb_pool * field3d)
!!$  !!
!!$  !-----------------------------------------------------------------------------------------------------
!!$  FUNCTION Multiply_n_struct_and_field3d(this,field3d) RESULT(struct3)
!!$    IMPLICIT NONE
!!$    ! 0.1 InOut
!!$    CLASS(t_slm_struct_sb_pool), INTENT(in) :: this
!!$    REAL(wp),                        INTENT(in) :: field3d(:,:,:)
!!$    TYPE(t_slm_struct_sb_pool), ALLOCATABLE :: struct3
!!$
!!$    ALLOCATE(struct3, source=this)
!!$
!!$    struct3%dom               = this%dom                * field3d
!!$    struct3%dom_assoc         = this%dom_assoc          * field3d
!!$    struct3%soluable_litter  = this%soluable_litter   * field3d
!!$    struct3%polymeric_litter = this%polymeric_litter  * field3d
!!$    struct3%woody_litter      = this%woody_litter       * field3d
!!$    struct3%microbial         = this%microbial          * field3d
!!$    struct3%residue           = this%residue            * field3d
!!$    struct3%residue_assoc     = this%residue_assoc      * field3d
!!$  END FUNCTION Multiply_n_struct_and_field3d
!!$
!!$  !-----------------------------------------------------------------------------------------------------
!!$  !> / (t_slm_struct_sb_pool / t_slm_struct_sb_pool)
!!$  !!
!!$  !-----------------------------------------------------------------------------------------------------
!!$  FUNCTION Divide_n_struct(this,struct2) RESULT(struct3)
!!$    IMPLICIT NONE
!!$    ! 0.1 InOut
!!$    CLASS(t_slm_struct_sb_pool), INTENT(in) :: this, struct2 
!!$    TYPE(t_slm_struct_sb_pool), ALLOCATABLE :: struct3       
!!$
!!$    ALLOCATE(struct3, source=this)
!!$
!!$    struct3%dom               = this%dom                / struct2%dom
!!$    struct3%dom_assoc         = this%dom_assoc          / struct2%dom_assoc
!!$    struct3%soluable_litter  = this%soluable_litter   / struct2%soluable_litter
!!$    struct3%polymeric_litter = this%polymeric_litter  / struct2%polymeric_litter
!!$    struct3%woody_litter      = this%woody_litter       / struct2%woody_litter
!!$    struct3%microbial         = this%microbial          / struct2%microbial
!!$    struct3%residue           = this%residue            / struct2%residue
!!$    struct3%residue_assoc     = this%residue_assoc      / struct2%residue_assoc
!!$  END FUNCTION Divide_n_struct
!!$
!!$  !-----------------------------------------------------------------------------------------------------
!!$  !> / (t_slm_struct_sb_pool / scalar)
!!$  !!
!!$  !-----------------------------------------------------------------------------------------------------
!!$  FUNCTION Divide_n_struct_and_scalar(this,scalar) RESULT(struct3)
!!$    IMPLICIT NONE
!!$    ! 0.1 InOut
!!$    CLASS(t_slm_struct_sb_pool), INTENT(in) :: this
!!$    REAL(wp),                        INTENT(in) :: scalar
!!$    TYPE(t_slm_struct_sb_pool), ALLOCATABLE :: struct3
!!$
!!$    ALLOCATE(struct3, source=this)
!!$
!!$    struct3%dom               = this%dom                / scalar
!!$    struct3%dom_assoc         = this%dom_assoc          / scalar
!!$    struct3%soluable_litter  = this%soluable_litter   / scalar
!!$    struct3%polymeric_litter = this%polymeric_litter  / scalar
!!$    struct3%woody_litter      = this%woody_litter       / scalar
!!$    struct3%microbial         = this%microbial          / scalar
!!$    struct3%residue           = this%residue            / scalar
!!$    struct3%residue_assoc     = this%residue_assoc      / scalar
!!$  END FUNCTION Divide_n_struct_and_scalar
!!$
!!$  !-----------------------------------------------------------------------------------------------------
!!$  !> / (t_slm_struct_sb_pool / field3d)
!!$  !!
!!$  !-----------------------------------------------------------------------------------------------------
!!$  FUNCTION Divide_n_struct_and_field3d(this,field3d) RESULT(struct3)
!!$    IMPLICIT NONE
!!$    ! 0.1 InOut
!!$    CLASS(t_slm_struct_sb_pool), INTENT(in) :: this
!!$    REAL(wp),                        INTENT(in) :: field3d(:,:,:)
!!$    TYPE(t_slm_struct_sb_pool), ALLOCATABLE :: struct3
!!$
!!$    ALLOCATE(struct3, source=this)
!!$
!!$    struct3%dom               = this%dom                / field3d
!!$    struct3%dom_assoc         = this%dom_assoc          / field3d
!!$    struct3%soluable_litter  = this%soluable_litter   / field3d
!!$    struct3%polymeric_litter = this%polymeric_litter  / field3d
!!$    struct3%woody_litter      = this%woody_litter       / field3d
!!$    struct3%microbial         = this%microbial          / field3d
!!$    struct3%residue           = this%residue            / field3d
!!$    struct3%residue_assoc     = this%residue_assoc      / field3d
!!$  END FUNCTION Divide_n_struct_and_field3d
!!$
  !-----------------------------------------------------------------------------------------------------
  !> Initialisation t_slm_struct_sb_pool with a scalar
  !!
  !-----------------------------------------------------------------------------------------------------
  PURE ELEMENTAL SUBROUTINE init_t_slm_struct_sb_pool(this,s)
    IMPLICIT NONE
    ! 0.1 InOut
    CLASS(t_slm_struct_sb_pool),INTENT(inout) :: this 
    REAL(wp),                       INTENT(in)    :: s

    CALL this%dom%init(s)
    CALL this%dom_assoc%init(s)
    CALL this%soluable_litter%init(s)
    CALL this%polymeric_litter%init(s)
    CALL this%woody_litter%init(s)
    CALL this%microbial%init(s)
    CALL this%residue%init(s)
    CALL this%residue_assoc%init(s)
  END SUBROUTINE init_t_slm_struct_sb_pool 



  !-----------------------------------------------------------------------------------------------------
  !> initialize memory for the processes in: soil
  !!
  !!  
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE Init_sb_memory(mem, prefix, suffix, lct_ids, model_id)

!    USE mo_jsb_varlist,       ONLY: t_jsb_varlist, jsb_add_var
!    USE mo_jsb_io,            ONLY: grib_bits, t_cf, t_grib1, t_grib2, TSTEP_CONSTANT, tables
!    USE mo_jsb_grid_class,    ONLY: t_jsb_grid, t_jsb_vgrid
!    USE mo_jsb_grid,          ONLY: Get_grid, Get_vgrid
    USE mo_base_memory,         ONLY: t_model

    ! 0.1.InOut
    CLASS(t_sb_memory), INTENT(inout), TARGET :: mem             !> soil memory
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

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':Init_sb_memory'

    model => Get_model(model_id)

!    table = tables(1)

!    hgrid   => Get_grid(mem%grid_id)
!    surface => Get_vgrid('surface')
!    canopy_layer => Get_vgrid('canopy_layer')


    IF (.NOT. ALLOCATED(mem%nh4_solute%ptr)) THEN

      !> add variables to soil memory soil\%VAR
      !! 
      !! soil\%VAR
      !!  
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
  
      CALL slm_add_var(mem, 'transport_nh4_solute', mem%transport_nh4_solute, &
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'transport_nh4_assoc', mem%transport_nh4_assoc, &
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'transport_no3_solute', mem%transport_no3_solute, &
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'transport_noy', mem%transport_noy, &
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'transport_n2o', mem%transport_n2o, &
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'transport_n2', mem%transport_n2, &
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'transport_nh4_n15_solute', mem%transport_nh4_n15_solute, &
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'transport_nh4_n15_assoc', mem%transport_nh4_n15_assoc, &
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'transport_no3_n15_solute', mem%transport_no3_n15_solute, &
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'transport_noy_n15', mem%transport_noy_n15, &
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'transport_n2o_n15', mem%transport_n2o_n15, &
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'transport_n2_n15', mem%transport_n2_n15, &
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'transport_po4_assoc', mem%transport_po4_assoc, &
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'transport_po4_solute', mem%transport_po4_solute, &
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'rtm_depolymerisation', mem%rtm_depolymerisation, &
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'rtm_mic_uptake', mem%rtm_mic_uptake, &
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'rtm_sorption', mem%rtm_sorption, &
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'rtm_desorption', mem%rtm_desorption, &
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'rtm_hsc', mem%rtm_hsc, &
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'rmm_depolymerisation', mem%rmm_depolymerisation, &
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'rmm_mic_uptake', mem%rmm_mic_uptake, &
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'rmm_sorption', mem%rmm_sorption, &
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'rmm_desorption', mem%rmm_desorption, &
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'rmm_hsc', mem%rmm_hsc, &
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'het_respiration', mem%het_respiration, &
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'het_respiration_c13', mem%het_respiration_c13, &
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'het_respiration_c14', mem%het_respiration_c14, &
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'net_mineralisation_n', mem%net_mineralisation_n, &
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'net_mineralisation_n15', mem%net_mineralisation_n15, &
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'net_mineralisation_p', mem%net_mineralisation_p, &
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'volatilisation_n', mem%volatilisation_n, &
            initval_r=0.0_wp) ! initvalue
  
      CALL slm_add_var(mem, 'nitrification', mem%nitrification, &
            initval_r=0.0_wp) ! initvalue

      CALL slm_add_var(mem, 'qmax_mineral_org', mem%qmax_mineral_org, &
            initval_r=0.0_wp) ! initvalue

      CALL slm_add_var(mem, 'k_bioturb', mem%k_bioturb, &
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
      CALL slm_add_var(mem, 'sb_pool'           , mem%pool)
      CALL slm_add_var(mem, 'sb_formation'      , mem%formation)
      CALL slm_add_var(mem, 'sb_loss'  , mem%loss)
      CALL slm_add_var(mem, 'sb_transport'      , mem%transport)

      !> level 2 TYPE (l2_struct) nested within the above nested TYPE structures (t_jsb_n_struct)
      !! 
      !! slm_add_var -> slm_add_l2_struct()
      !! 
      CALL slm_add_var(mem%pool,    'sb_pool_dom',                mem%pool%dom) 
      CALL slm_add_var(mem%pool,    'sb_pool_dom_assoc',          mem%pool%dom_assoc)
      CALL slm_add_var(mem%pool,    'sb_pool_soluable_litter',   mem%pool%soluable_litter) 
      CALL slm_add_var(mem%pool,    'sb_pool_polymeric_litter',  mem%pool%polymeric_litter) 
      CALL slm_add_var(mem%pool,    'sb_pool_woody_litter',       mem%pool%woody_litter) 
      CALL slm_add_var(mem%pool,    'sb_pool_microbial',          mem%pool%microbial) 
      CALL slm_add_var(mem%pool,    'sb_pool_residue',            mem%pool%residue) 
      CALL slm_add_var(mem%pool,    'sb_pool_residue_assoc',      mem%pool%residue_assoc)  

      CALL slm_add_var(mem%formation,    'sb_formation_dom',                mem%formation%dom) 
      CALL slm_add_var(mem%formation,    'sb_formation_dom_assoc',          mem%formation%dom_assoc)
      CALL slm_add_var(mem%formation,    'sb_formation_soluable_litter',   mem%formation%soluable_litter) 
      CALL slm_add_var(mem%formation,    'sb_formation_polymeric_litter',  mem%formation%polymeric_litter) 
      CALL slm_add_var(mem%formation,    'sb_formation_woody_litter',       mem%formation%woody_litter) 
      CALL slm_add_var(mem%formation,    'sb_formation_microbial',          mem%formation%microbial) 
      CALL slm_add_var(mem%formation,    'sb_formation_residue',            mem%formation%residue) 
      CALL slm_add_var(mem%formation,    'sb_formation_residue_assoc',      mem%formation%residue_assoc)  

      CALL slm_add_var(mem%loss,    'sb_loss_dom',                mem%loss%dom) 
      CALL slm_add_var(mem%loss,    'sb_loss_dom_assoc',          mem%loss%dom_assoc)
      CALL slm_add_var(mem%loss,    'sb_loss_soluable_litter',   mem%loss%soluable_litter) 
      CALL slm_add_var(mem%loss,    'sb_loss_polymeric_litter',  mem%loss%polymeric_litter) 
      CALL slm_add_var(mem%loss,    'sb_loss_woody_litter',       mem%loss%woody_litter) 
      CALL slm_add_var(mem%loss,    'sb_loss_microbial',          mem%loss%microbial) 
      CALL slm_add_var(mem%loss,    'sb_loss_residue',            mem%loss%residue) 
      CALL slm_add_var(mem%loss,    'sb_loss_residue_assoc',      mem%loss%residue_assoc)  

      CALL slm_add_var(mem%transport,    'sb_transport_dom',                mem%transport%dom) 
      CALL slm_add_var(mem%transport,    'sb_transport_dom_assoc',          mem%transport%dom_assoc)
      CALL slm_add_var(mem%transport,    'sb_transport_soluable_litter',   mem%transport%soluable_litter) 
      CALL slm_add_var(mem%transport,    'sb_transport_polymeric_litter',  mem%transport%polymeric_litter) 
      CALL slm_add_var(mem%transport,    'sb_transport_woody_litter',       mem%transport%woody_litter) 
      CALL slm_add_var(mem%transport,    'sb_transport_microbial',          mem%transport%microbial) 
      CALL slm_add_var(mem%transport,    'sb_transport_residue',            mem%transport%residue) 
      CALL slm_add_var(mem%transport,    'sb_transport_residue_assoc',      mem%transport%residue_assoc)  

      !> add variables to the structs and the l2 structs
      !!
      !! soil\%t_jsb_struct\%VAR
      !! 
      !! soil\%t_jsb_n_struct\%t_jsb_struct\%VAR
      !! 
      CALL slm_add_var(mem%pool%dom, 'sb_pool_dom_carbon', mem%pool%dom%carbon, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%pool%dom, 'sb_pool_dom_nitrogen', mem%pool%dom%nitrogen, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%pool%dom, 'sb_pool_dom_phosphorus', mem%pool%dom%phosphorus, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%pool%dom, 'sb_pool_dom_carbon13', mem%pool%dom%carbon13, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%pool%dom, 'sb_pool_dom_carbon14', mem%pool%dom%carbon14, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%pool%dom, 'sb_pool_dom_nitrogen15', mem%pool%dom%nitrogen15, &
            initval_r=0.0_wp)     ! initvalue
 
      CALL slm_add_var(mem%pool%dom_assoc, 'sb_pool_dom_assoc_carbon', mem%pool%dom_assoc%carbon, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%pool%dom_assoc, 'sb_pool_dom_assoc_nitrogen', mem%pool%dom_assoc%nitrogen, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%pool%dom_assoc, 'sb_pool_dom_assoc_phosphorus', mem%pool%dom_assoc%phosphorus, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%pool%dom_assoc, 'sb_pool_dom_assoc_carbon13', mem%pool%dom_assoc%carbon13, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%pool%dom_assoc, 'sb_pool_dom_assoc_carbon14', mem%pool%dom_assoc%carbon14, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%pool%dom_assoc, 'sb_pool_dom_assoc_nitrogen15', mem%pool%dom_assoc%nitrogen15, &
            initval_r=0.0_wp)     ! initvalue
 
      CALL slm_add_var(mem%pool%soluable_litter, 'sb_pool_soluable_litter_carbon', mem%pool%soluable_litter%carbon, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%pool%soluable_litter, 'sb_pool_soluable_litter_nitrogen', mem%pool%soluable_litter%nitrogen, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%pool%soluable_litter, 'sb_pool_soluable_litter_phosphorus', mem%pool%soluable_litter%phosphorus, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%pool%soluable_litter, 'sb_pool_soluable_litter_carbon13', mem%pool%soluable_litter%carbon13, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%pool%soluable_litter, 'sb_pool_soluable_litter_carbon14', mem%pool%soluable_litter%carbon14, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%pool%soluable_litter, 'sb_pool_soluable_litter_nitrogen15', mem%pool%soluable_litter%nitrogen15, &
            initval_r=0.0_wp)     ! initvalue
 
      CALL slm_add_var(mem%pool%polymeric_litter, 'sb_pool_polymeric_litter_carbon', mem%pool%polymeric_litter%carbon, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%pool%polymeric_litter, 'sb_pool_polymeric_litter_nitrogen', mem%pool%polymeric_litter%nitrogen, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%pool%polymeric_litter, 'sb_pool_polymeric_litter_phosphorus', mem%pool%polymeric_litter%phosphorus, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%pool%polymeric_litter, 'sb_pool_polymeric_litter_carbon13', mem%pool%polymeric_litter%carbon13, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%pool%polymeric_litter, 'sb_pool_polymeric_litter_carbon14', mem%pool%polymeric_litter%carbon14, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%pool%polymeric_litter, 'sb_pool_polymeric_litter_nitrogen15', mem%pool%polymeric_litter%nitrogen15, &
            initval_r=0.0_wp)     ! initvalue
 
      CALL slm_add_var(mem%pool%woody_litter, 'sb_pool_woody_litter_carbon', mem%pool%woody_litter%carbon, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%pool%woody_litter, 'sb_pool_woody_litter_nitrogen', mem%pool%woody_litter%nitrogen, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%pool%woody_litter, 'sb_pool_woody_litter_phosphorus', mem%pool%woody_litter%phosphorus, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%pool%woody_litter, 'sb_pool_woody_litter_carbon13', mem%pool%woody_litter%carbon13, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%pool%woody_litter, 'sb_pool_woody_litter_carbon14', mem%pool%woody_litter%carbon14, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%pool%woody_litter, 'sb_pool_woody_litter_nitrogen15', mem%pool%woody_litter%nitrogen15, &
            initval_r=0.0_wp)     ! initvalue
 
      CALL slm_add_var(mem%pool%microbial, 'sb_pool_microbial_carbon', mem%pool%microbial%carbon, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%pool%microbial, 'sb_pool_microbial_nitrogen', mem%pool%microbial%nitrogen, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%pool%microbial, 'sb_pool_microbial_phosphorus', mem%pool%microbial%phosphorus, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%pool%microbial, 'sb_pool_microbial_carbon13', mem%pool%microbial%carbon13, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%pool%microbial, 'sb_pool_microbial_carbon14', mem%pool%microbial%carbon14, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%pool%microbial, 'sb_pool_microbial_nitrogen15', mem%pool%microbial%nitrogen15, &
            initval_r=0.0_wp)     ! initvalue
 
      CALL slm_add_var(mem%pool%residue, 'sb_pool_residue_carbon', mem%pool%residue%carbon, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%pool%residue, 'sb_pool_residue_nitrogen', mem%pool%residue%nitrogen, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%pool%residue, 'sb_pool_residue_phosphorus', mem%pool%residue%phosphorus, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%pool%residue, 'sb_pool_residue_carbon13', mem%pool%residue%carbon13, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%pool%residue, 'sb_pool_residue_carbon14', mem%pool%residue%carbon14, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%pool%residue, 'sb_pool_residue_nitrogen15', mem%pool%residue%nitrogen15, &
            initval_r=0.0_wp)     ! initvalue
 
      CALL slm_add_var(mem%pool%residue_assoc, 'sb_pool_residue_assoc_carbon', mem%pool%residue_assoc%carbon, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%pool%residue_assoc, 'sb_pool_residue_assoc_nitrogen', mem%pool%residue_assoc%nitrogen, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%pool%residue_assoc, 'sb_pool_residue_assoc_phosphorus', mem%pool%residue_assoc%phosphorus, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%pool%residue_assoc, 'sb_pool_residue_assoc_carbon13', mem%pool%residue_assoc%carbon13, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%pool%residue_assoc, 'sb_pool_residue_assoc_carbon14', mem%pool%residue_assoc%carbon14, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%pool%residue_assoc, 'sb_pool_residue_assoc_nitrogen15', mem%pool%residue_assoc%nitrogen15, &
            initval_r=0.0_wp)     ! initvalue

! formation
      CALL slm_add_var(mem%formation%dom, 'sb_formation_dom_carbon', mem%formation%dom%carbon, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%formation%dom, 'sb_formation_dom_nitrogen', mem%formation%dom%nitrogen, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%formation%dom, 'sb_formation_dom_phosphorus', mem%formation%dom%phosphorus, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%formation%dom, 'sb_formation_dom_carbon13', mem%formation%dom%carbon13, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%formation%dom, 'sb_formation_dom_carbon14', mem%formation%dom%carbon14, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%formation%dom, 'sb_formation_dom_nitrogen15', mem%formation%dom%nitrogen15, &
            initval_r=0.0_wp)     ! initvalue
 
      CALL slm_add_var(mem%formation%dom_assoc, 'sb_formation_dom_assoc_carbon', mem%formation%dom_assoc%carbon, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%formation%dom_assoc, 'sb_formation_dom_assoc_nitrogen', mem%formation%dom_assoc%nitrogen, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%formation%dom_assoc, 'sb_formation_dom_assoc_phosphorus', mem%formation%dom_assoc%phosphorus, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%formation%dom_assoc, 'sb_formation_dom_assoc_carbon13', mem%formation%dom_assoc%carbon13, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%formation%dom_assoc, 'sb_formation_dom_assoc_carbon14', mem%formation%dom_assoc%carbon14, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%formation%dom_assoc, 'sb_formation_dom_assoc_nitrogen15', mem%formation%dom_assoc%nitrogen15, &
            initval_r=0.0_wp)     ! initvalue
 
      CALL slm_add_var(mem%formation%soluable_litter, 'sb_formation_soluable_litter_carbon', mem%formation%soluable_litter%carbon, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%formation%soluable_litter, 'sb_formation_soluable_litter_nitrogen', mem%formation%soluable_litter%nitrogen, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%formation%soluable_litter, 'sb_formation_soluable_litter_phosphorus', mem%formation%soluable_litter%phosphorus, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%formation%soluable_litter, 'sb_formation_soluable_litter_carbon13', mem%formation%soluable_litter%carbon13, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%formation%soluable_litter, 'sb_formation_soluable_litter_carbon14', mem%formation%soluable_litter%carbon14, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%formation%soluable_litter, 'sb_formation_soluable_litter_nitrogen15', mem%formation%soluable_litter%nitrogen15, &
            initval_r=0.0_wp)     ! initvalue
 
      CALL slm_add_var(mem%formation%polymeric_litter, 'sb_formation_polymeric_litter_carbon', mem%formation%polymeric_litter%carbon, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%formation%polymeric_litter, 'sb_formation_polymeric_litter_nitrogen', mem%formation%polymeric_litter%nitrogen, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%formation%polymeric_litter, 'sb_formation_polymeric_litter_phosphorus', mem%formation%polymeric_litter%phosphorus, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%formation%polymeric_litter, 'sb_formation_polymeric_litter_carbon13', mem%formation%polymeric_litter%carbon13, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%formation%polymeric_litter, 'sb_formation_polymeric_litter_carbon14', mem%formation%polymeric_litter%carbon14, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%formation%polymeric_litter, 'sb_formation_polymeric_litter_nitrogen15', mem%formation%polymeric_litter%nitrogen15, &
            initval_r=0.0_wp)     ! initvalue
 
      CALL slm_add_var(mem%formation%woody_litter, 'sb_formation_woody_litter_carbon', mem%formation%woody_litter%carbon, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%formation%woody_litter, 'sb_formation_woody_litter_nitrogen', mem%formation%woody_litter%nitrogen, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%formation%woody_litter, 'sb_formation_woody_litter_phosphorus', mem%formation%woody_litter%phosphorus, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%formation%woody_litter, 'sb_formation_woody_litter_carbon13', mem%formation%woody_litter%carbon13, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%formation%woody_litter, 'sb_formation_woody_litter_carbon14', mem%formation%woody_litter%carbon14, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%formation%woody_litter, 'sb_formation_woody_litter_nitrogen15', mem%formation%woody_litter%nitrogen15, &
            initval_r=0.0_wp)     ! initvalue
 
      CALL slm_add_var(mem%formation%microbial, 'sb_formation_microbial_carbon', mem%formation%microbial%carbon, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%formation%microbial, 'sb_formation_microbial_nitrogen', mem%formation%microbial%nitrogen, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%formation%microbial, 'sb_formation_microbial_phosphorus', mem%formation%microbial%phosphorus, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%formation%microbial, 'sb_formation_microbial_carbon13', mem%formation%microbial%carbon13, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%formation%microbial, 'sb_formation_microbial_carbon14', mem%formation%microbial%carbon14, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%formation%microbial, 'sb_formation_microbial_nitrogen15', mem%formation%microbial%nitrogen15, &
            initval_r=0.0_wp)     ! initvalue
 
      CALL slm_add_var(mem%formation%residue, 'sb_formation_residue_carbon', mem%formation%residue%carbon, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%formation%residue, 'sb_formation_residue_nitrogen', mem%formation%residue%nitrogen, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%formation%residue, 'sb_formation_residue_phosphorus', mem%formation%residue%phosphorus, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%formation%residue, 'sb_formation_residue_carbon13', mem%formation%residue%carbon13, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%formation%residue, 'sb_formation_residue_carbon14', mem%formation%residue%carbon14, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%formation%residue, 'sb_formation_residue_nitrogen15', mem%formation%residue%nitrogen15, &
            initval_r=0.0_wp)     ! initvalue
 
      CALL slm_add_var(mem%formation%residue_assoc, 'sb_formation_residue_assoc_carbon', mem%formation%residue_assoc%carbon, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%formation%residue_assoc, 'sb_formation_residue_assoc_nitrogen', mem%formation%residue_assoc%nitrogen, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%formation%residue_assoc, 'sb_formation_residue_assoc_phosphorus', mem%formation%residue_assoc%phosphorus, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%formation%residue_assoc, 'sb_formation_residue_assoc_carbon13', mem%formation%residue_assoc%carbon13, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%formation%residue_assoc, 'sb_formation_residue_assoc_carbon14', mem%formation%residue_assoc%carbon14, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%formation%residue_assoc, 'sb_formation_residue_assoc_nitrogen15', mem%formation%residue_assoc%nitrogen15, &
            initval_r=0.0_wp)     ! initvalue

!loss
      CALL slm_add_var(mem%loss%dom, 'sb_loss_dom_carbon', mem%loss%dom%carbon, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%loss%dom, 'sb_loss_dom_nitrogen', mem%loss%dom%nitrogen, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%loss%dom, 'sb_loss_dom_phosphorus', mem%loss%dom%phosphorus, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%loss%dom, 'sb_loss_dom_carbon13', mem%loss%dom%carbon13, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%loss%dom, 'sb_loss_dom_carbon14', mem%loss%dom%carbon14, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%loss%dom, 'sb_loss_dom_nitrogen15', mem%loss%dom%nitrogen15, &
            initval_r=0.0_wp)     ! initvalue
 
      CALL slm_add_var(mem%loss%dom_assoc, 'sb_loss_dom_assoc_carbon', mem%loss%dom_assoc%carbon, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%loss%dom_assoc, 'sb_loss_dom_assoc_nitrogen', mem%loss%dom_assoc%nitrogen, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%loss%dom_assoc, 'sb_loss_dom_assoc_phosphorus', mem%loss%dom_assoc%phosphorus, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%loss%dom_assoc, 'sb_loss_dom_assoc_carbon13', mem%loss%dom_assoc%carbon13, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%loss%dom_assoc, 'sb_loss_dom_assoc_carbon14', mem%loss%dom_assoc%carbon14, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%loss%dom_assoc, 'sb_loss_dom_assoc_nitrogen15', mem%loss%dom_assoc%nitrogen15, &
            initval_r=0.0_wp)     ! initvalue
 
      CALL slm_add_var(mem%loss%soluable_litter, 'sb_loss_soluable_litter_carbon', mem%loss%soluable_litter%carbon, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%loss%soluable_litter, 'sb_loss_soluable_litter_nitrogen', mem%loss%soluable_litter%nitrogen, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%loss%soluable_litter, 'sb_loss_soluable_litter_phosphorus', mem%loss%soluable_litter%phosphorus, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%loss%soluable_litter, 'sb_loss_soluable_litter_carbon13', mem%loss%soluable_litter%carbon13, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%loss%soluable_litter, 'sb_loss_soluable_litter_carbon14', mem%loss%soluable_litter%carbon14, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%loss%soluable_litter, 'sb_loss_soluable_litter_nitrogen15', mem%loss%soluable_litter%nitrogen15, &
            initval_r=0.0_wp)     ! initvalue
 
      CALL slm_add_var(mem%loss%polymeric_litter, 'sb_loss_polymeric_litter_carbon', mem%loss%polymeric_litter%carbon, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%loss%polymeric_litter, 'sb_loss_polymeric_litter_nitrogen', mem%loss%polymeric_litter%nitrogen, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%loss%polymeric_litter, 'sb_loss_polymeric_litter_phosphorus', mem%loss%polymeric_litter%phosphorus, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%loss%polymeric_litter, 'sb_loss_polymeric_litter_carbon13', mem%loss%polymeric_litter%carbon13, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%loss%polymeric_litter, 'sb_loss_polymeric_litter_carbon14', mem%loss%polymeric_litter%carbon14, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%loss%polymeric_litter, 'sb_loss_polymeric_litter_nitrogen15', mem%loss%polymeric_litter%nitrogen15, &
            initval_r=0.0_wp)     ! initvalue
 
      CALL slm_add_var(mem%loss%woody_litter, 'sb_loss_woody_litter_carbon', mem%loss%woody_litter%carbon, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%loss%woody_litter, 'sb_loss_woody_litter_nitrogen', mem%loss%woody_litter%nitrogen, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%loss%woody_litter, 'sb_loss_woody_litter_phosphorus', mem%loss%woody_litter%phosphorus, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%loss%woody_litter, 'sb_loss_woody_litter_carbon13', mem%loss%woody_litter%carbon13, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%loss%woody_litter, 'sb_loss_woody_litter_carbon14', mem%loss%woody_litter%carbon14, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%loss%woody_litter, 'sb_loss_woody_litter_nitrogen15', mem%loss%woody_litter%nitrogen15, &
            initval_r=0.0_wp)     ! initvalue
 
      CALL slm_add_var(mem%loss%microbial, 'sb_loss_microbial_carbon', mem%loss%microbial%carbon, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%loss%microbial, 'sb_loss_microbial_nitrogen', mem%loss%microbial%nitrogen, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%loss%microbial, 'sb_loss_microbial_phosphorus', mem%loss%microbial%phosphorus, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%loss%microbial, 'sb_loss_microbial_carbon13', mem%loss%microbial%carbon13, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%loss%microbial, 'sb_loss_microbial_carbon14', mem%loss%microbial%carbon14, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%loss%microbial, 'sb_loss_microbial_nitrogen15', mem%loss%microbial%nitrogen15, &
            initval_r=0.0_wp)     ! initvalue
 
      CALL slm_add_var(mem%loss%residue, 'sb_loss_residue_carbon', mem%loss%residue%carbon, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%loss%residue, 'sb_loss_residue_nitrogen', mem%loss%residue%nitrogen, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%loss%residue, 'sb_loss_residue_phosphorus', mem%loss%residue%phosphorus, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%loss%residue, 'sb_loss_residue_carbon13', mem%loss%residue%carbon13, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%loss%residue, 'sb_loss_residue_carbon14', mem%loss%residue%carbon14, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%loss%residue, 'sb_loss_residue_nitrogen15', mem%loss%residue%nitrogen15, &
            initval_r=0.0_wp)     ! initvalue
 
      CALL slm_add_var(mem%loss%residue_assoc, 'sb_loss_residue_assoc_carbon', mem%loss%residue_assoc%carbon, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%loss%residue_assoc, 'sb_loss_residue_assoc_nitrogen', mem%loss%residue_assoc%nitrogen, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%loss%residue_assoc, 'sb_loss_residue_assoc_phosphorus', mem%loss%residue_assoc%phosphorus, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%loss%residue_assoc, 'sb_loss_residue_assoc_carbon13', mem%loss%residue_assoc%carbon13, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%loss%residue_assoc, 'sb_loss_residue_assoc_carbon14', mem%loss%residue_assoc%carbon14, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%loss%residue_assoc, 'sb_loss_residue_assoc_nitrogen15', mem%loss%residue_assoc%nitrogen15, &
            initval_r=0.0_wp)     ! initvalue

! Transport
      CALL slm_add_var(mem%transport%dom, 'sb_transport_dom_carbon', mem%transport%dom%carbon, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%transport%dom, 'sb_transport_dom_nitrogen', mem%transport%dom%nitrogen, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%transport%dom, 'sb_transport_dom_phosphorus', mem%transport%dom%phosphorus, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%transport%dom, 'sb_transport_dom_carbon13', mem%transport%dom%carbon13, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%transport%dom, 'sb_transport_dom_carbon14', mem%transport%dom%carbon14, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%transport%dom, 'sb_transport_dom_nitrogen15', mem%transport%dom%nitrogen15, &
            initval_r=0.0_wp)     ! initvalue
 
      CALL slm_add_var(mem%transport%dom_assoc, 'sb_transport_dom_assoc_carbon', mem%transport%dom_assoc%carbon, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%transport%dom_assoc, 'sb_transport_dom_assoc_nitrogen', mem%transport%dom_assoc%nitrogen, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%transport%dom_assoc, 'sb_transport_dom_assoc_phosphorus', mem%transport%dom_assoc%phosphorus, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%transport%dom_assoc, 'sb_transport_dom_assoc_carbon13', mem%transport%dom_assoc%carbon13, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%transport%dom_assoc, 'sb_transport_dom_assoc_carbon14', mem%transport%dom_assoc%carbon14, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%transport%dom_assoc, 'sb_transport_dom_assoc_nitrogen15', mem%transport%dom_assoc%nitrogen15, &
            initval_r=0.0_wp)     ! initvalue
 
      CALL slm_add_var(mem%transport%soluable_litter, 'sb_transport_soluable_litter_carbon', mem%transport%soluable_litter%carbon, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%transport%soluable_litter, 'sb_transport_soluable_litter_nitrogen', mem%transport%soluable_litter%nitrogen, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%transport%soluable_litter, 'sb_transport_soluable_litter_phosphorus', mem%transport%soluable_litter%phosphorus, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%transport%soluable_litter, 'sb_transport_soluable_litter_carbon13', mem%transport%soluable_litter%carbon13, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%transport%soluable_litter, 'sb_transport_soluable_litter_carbon14', mem%transport%soluable_litter%carbon14, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%transport%soluable_litter, 'sb_transport_soluable_litter_nitrogen15', mem%transport%soluable_litter%nitrogen15, &
            initval_r=0.0_wp)     ! initvalue
 
      CALL slm_add_var(mem%transport%polymeric_litter, 'sb_transport_polymeric_litter_carbon', mem%transport%polymeric_litter%carbon, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%transport%polymeric_litter, 'sb_transport_polymeric_litter_nitrogen', mem%transport%polymeric_litter%nitrogen, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%transport%polymeric_litter, 'sb_transport_polymeric_litter_phosphorus', mem%transport%polymeric_litter%phosphorus, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%transport%polymeric_litter, 'sb_transport_polymeric_litter_carbon13', mem%transport%polymeric_litter%carbon13, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%transport%polymeric_litter, 'sb_transport_polymeric_litter_carbon14', mem%transport%polymeric_litter%carbon14, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%transport%polymeric_litter, 'sb_transport_polymeric_litter_nitrogen15', mem%transport%polymeric_litter%nitrogen15, &
            initval_r=0.0_wp)     ! initvalue
 
      CALL slm_add_var(mem%transport%woody_litter, 'sb_transport_woody_litter_carbon', mem%transport%woody_litter%carbon, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%transport%woody_litter, 'sb_transport_woody_litter_nitrogen', mem%transport%woody_litter%nitrogen, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%transport%woody_litter, 'sb_transport_woody_litter_phosphorus', mem%transport%woody_litter%phosphorus, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%transport%woody_litter, 'sb_transport_woody_litter_carbon13', mem%transport%woody_litter%carbon13, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%transport%woody_litter, 'sb_transport_woody_litter_carbon14', mem%transport%woody_litter%carbon14, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%transport%woody_litter, 'sb_transport_woody_litter_nitrogen15', mem%transport%woody_litter%nitrogen15, &
            initval_r=0.0_wp)     ! initvalue
 
      CALL slm_add_var(mem%transport%microbial, 'sb_transport_microbial_carbon', mem%transport%microbial%carbon, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%transport%microbial, 'sb_transport_microbial_nitrogen', mem%transport%microbial%nitrogen, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%transport%microbial, 'sb_transport_microbial_phosphorus', mem%transport%microbial%phosphorus, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%transport%microbial, 'sb_transport_microbial_carbon13', mem%transport%microbial%carbon13, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%transport%microbial, 'sb_transport_microbial_carbon14', mem%transport%microbial%carbon14, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%transport%microbial, 'sb_transport_microbial_nitrogen15', mem%transport%microbial%nitrogen15, &
            initval_r=0.0_wp)     ! initvalue
 
      CALL slm_add_var(mem%transport%residue, 'sb_transport_residue_carbon', mem%transport%residue%carbon, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%transport%residue, 'sb_transport_residue_nitrogen', mem%transport%residue%nitrogen, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%transport%residue, 'sb_transport_residue_phosphorus', mem%transport%residue%phosphorus, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%transport%residue, 'sb_transport_residue_carbon13', mem%transport%residue%carbon13, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%transport%residue, 'sb_transport_residue_carbon14', mem%transport%residue%carbon14, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%transport%residue, 'sb_transport_residue_nitrogen15', mem%transport%residue%nitrogen15, &
            initval_r=0.0_wp)     ! initvalue
 
      CALL slm_add_var(mem%transport%residue_assoc, 'sb_transport_residue_assoc_carbon', mem%transport%residue_assoc%carbon, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%transport%residue_assoc, 'sb_transport_residue_assoc_nitrogen', mem%transport%residue_assoc%nitrogen, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%transport%residue_assoc, 'sb_transport_residue_assoc_phosphorus', mem%transport%residue_assoc%phosphorus, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%transport%residue_assoc, 'sb_transport_residue_assoc_carbon13', mem%transport%residue_assoc%carbon13, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%transport%residue_assoc, 'sb_transport_residue_assoc_carbon14', mem%transport%residue_assoc%carbon14, &
            initval_r=0.0_wp)     ! initvalue
      CALL slm_add_var(mem%transport%residue_assoc, 'sb_transport_residue_assoc_nitrogen15', mem%transport%residue_assoc%nitrogen15, &
            initval_r=0.0_wp)     ! initvalue

    END IF

  END SUBROUTINE Init_sb_memory

END MODULE mo_sb_memory_class
