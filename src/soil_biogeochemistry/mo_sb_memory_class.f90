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

  USE mo_slm_varlist,            ONLY: jsb_add_var   ! replace by mo_jsb_varlist later    *slmtemp



  IMPLICIT NONE
  PRIVATE
  PUBLIC :: t_sb_memory, max_no_of_vars
  PUBLIC :: max_no_of_structs, max_no_of_n_structs          ! *slmjsb4
  PUBLIC :: t_sb_element_a, t_sb_pool_a                     ! *slmtemp
  PUBLIC :: Init_sb_pools_elements_a, &                     ! *slmtemp
            Deallocate_var_sb_pools_elements_a
  PUBLIC :: sb_pool, sb_formation, sb_loss, sb_transport    ! *slmtemp
  PUBLIC :: sb_delta_pools, sb_mycorrhiza_export            ! *slmtemp

  INTEGER, PARAMETER :: max_no_of_vars      = 150   !> TYPE t_jsb_var
  INTEGER, PARAMETER :: max_no_of_structs   = 20   !> TYPE t_jsb_struct     *slmjsb4
  INTEGER, PARAMETER :: max_no_of_n_structs = 20   !> TYPE t_jsb_n_struct   *slmjsb4


  !-----------------------------------------------------------------------------------------------------
  !> Type definition for helper elements with ALLOCATABLE variables for elements in sb memory
  !! 
  !! 
  !-----------------------------------------------------------------------------------------------------
  TYPE t_sb_element_a
    REAL(wp), ALLOCATABLE :: carbon(:,:)      ! dimension of variable in memory is var(ics:ice,soil_layers,iblk)
    REAL(wp), ALLOCATABLE :: nitrogen(:,:)
    REAL(wp), ALLOCATABLE :: phosphorus(:,:)
    REAL(wp), ALLOCATABLE :: carbon13(:,:)
    REAL(wp), ALLOCATABLE :: carbon14(:,:)
    REAL(wp), ALLOCATABLE :: nitrogen15(:,:)
  CONTAINS
    PRIVATE
    !! operator overload
    ! +
    PROCEDURE               :: Add_element_a                   ! Add two elements
    PROCEDURE               :: Add_element_a_and_scalar        ! Add scalar to element
    PROCEDURE               :: Add_element_a_and_field2d       ! Add REAL field2d ( var(ics:ice,soil_layers) ) to element
    GENERIC, PUBLIC         :: OPERATOR(+) => Add_element_a, Add_element_a_and_scalar, Add_element_a_and_field2d
    !! operator overload
    ! -
    PROCEDURE               :: Subtract_element_a                   ! Subtract one element from another element
    PROCEDURE               :: Subtract_element_a_and_scalar        ! Subtract scalar from element
    PROCEDURE               :: Subtract_element_a_and_field2d       ! Subtract field2d ( var(ics:ice,soil_layers) ) from element
    GENERIC, PUBLIC         :: OPERATOR(-) => Subtract_element_a, Subtract_element_a_and_scalar, Subtract_element_a_and_field2d
    !! operator overload
    ! /
    PROCEDURE               :: Divide_element_a                   ! Divide two elements
    PROCEDURE               :: Divide_element_a_and_scalar        ! Divide scalar to element
    PROCEDURE               :: Divide_element_a_and_field2d       ! Divide REAL field2d ( var(ics:ice,soil_layers) ) to element
    GENERIC, PUBLIC         :: OPERATOR(/) => Divide_element_a, Divide_element_a_and_scalar, Divide_element_a_and_field2d
    !! operator overload
    ! *
    PROCEDURE               :: Multiply_element_a                   ! Multiply two elements
    PROCEDURE               :: Multiply_element_a_and_scalar        ! Multiply element with 
    PROCEDURE, PASS(this2)  :: Multiply_scalar_and_element_a        ! Multiply scalar with element
    PROCEDURE               :: Multiply_element_a_and_field2d       ! Multiply element with REAL field2d ( var(ics:ice,soil_layers) )
    PROCEDURE, PASS(this2)  :: Multiply_field2d_and_element_a       ! Multiply REAL field2d ( var(ics:ice,soil_layers) ) with element
    GENERIC, PUBLIC         :: OPERATOR(*) => Multiply_element_a, Multiply_element_a_and_scalar, Multiply_scalar_and_element_a, &
                                              Multiply_element_a_and_field2d, Multiply_field2d_and_element_a
  END TYPE t_sb_element_a


  !-----------------------------------------------------------------------------------------------------
  !> Type definition for helper pools with ALLOCATABLE variables for pools in sb memory
  !! 
  !! 
  !-----------------------------------------------------------------------------------------------------
  TYPE t_sb_pool_a
    TYPE(t_sb_element_a)  :: dom                 !< dissolved organic matter
    TYPE(t_sb_element_a)  :: dom_assoc           !< minerally associated disolved organic matter
    TYPE(t_sb_element_a)  :: soluable_litter     !< metabolic litter (sugars, cellulose etc.)
    TYPE(t_sb_element_a)  :: polymeric_litter    !< structural litter (lignified litter)
    TYPE(t_sb_element_a)  :: woody_litter        !< physically protected wood litter
    TYPE(t_sb_element_a)  :: fungi               !< saprophytic fungal biomass
    TYPE(t_sb_element_a)  :: mycorrhiza          !< mycorrhizal fungal biomass 
    TYPE(t_sb_element_a)  :: microbial           !< microbial biomass (bacteria!)
    TYPE(t_sb_element_a)  :: residue             !< necromass of fungi, microbial and mycorrhiza
    TYPE(t_sb_element_a)  :: residue_assoc       !< minerally associated microbial necromass 
  CONTAINS
    PRIVATE
    !! operator overload
    ! +
    PROCEDURE               :: Add_pool_a                   ! Add two pools
    PROCEDURE               :: Add_pool_a_and_scalar        ! Add scalar to pool
    PROCEDURE               :: Add_pool_a_and_field2d       ! Add REAL field2d ( var(ics:ice,soil_layers) ) to pool
    GENERIC, PUBLIC         :: OPERATOR(+) => Add_pool_a, Add_pool_a_and_scalar, Add_pool_a_and_field2d
    !! operator overload
    ! -
    PROCEDURE               :: Subtract_pool_a                   ! Subtract one pool from another pool
    PROCEDURE               :: Subtract_pool_a_and_scalar        ! Subtract scalar from pool
    PROCEDURE               :: Subtract_pool_a_and_field2d       ! Subtract field2d ( var(ics:ice,soil_layers) ) from pool
    GENERIC, PUBLIC         :: OPERATOR(-) => Subtract_pool_a, Subtract_pool_a_and_scalar, Subtract_pool_a_and_field2d
    !! operator overload
    ! /
    PROCEDURE               :: Divide_pool_a                   ! Divide two pools
    PROCEDURE               :: Divide_pool_a_and_scalar        ! Divide scalar to pool
    PROCEDURE               :: Divide_pool_a_and_field2d       ! Divide REAL field2d ( var(ics:ice,soil_layers) ) to pool
    GENERIC, PUBLIC         :: OPERATOR(/) => Divide_pool_a, Divide_pool_a_and_scalar, Divide_pool_a_and_field2d
    !! operator overload
    ! *
    PROCEDURE               :: Multiply_pool_a                   ! Multiply two pools
    PROCEDURE               :: Multiply_pool_a_and_scalar        ! Multiply pool with scalar
    PROCEDURE, PASS(this2)  :: Multiply_scalar_and_pool_a        ! Multiply scalar with pool
    PROCEDURE               :: Multiply_pool_a_and_field2d       ! Multiply pool with REAL 1D field2d ( var(ics:ice) )
    PROCEDURE, PASS(this2)  :: Multiply_field2d_and_pool_a       ! Multiply REAL 1D field2d ( var(ics:ice,soil_layers) ) with pool
    GENERIC, PUBLIC         :: OPERATOR(*) => Multiply_pool_a, Multiply_pool_a_and_scalar, Multiply_scalar_and_pool_a, &
                                              Multiply_pool_a_and_field2d, Multiply_field2d_and_pool_a
  END TYPE t_sb_pool_a


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
    PROCEDURE  :: init => init_t_slm_struct_sb_element
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
    TYPE(t_slm_struct_sb_element)  :: soluable_litter     !< metabolic litter (sugars, cellulose etc.)
    TYPE(t_slm_struct_sb_element)  :: polymeric_litter    !< structural litter (lignified litter)
    TYPE(t_slm_struct_sb_element)  :: woody_litter        !< physically protected wood litter
    TYPE(t_slm_struct_sb_element)  :: fungi               !< saprophytic fungal biomass
    TYPE(t_slm_struct_sb_element)  :: mycorrhiza          !< mycorrhizal fungal biomass 
    TYPE(t_slm_struct_sb_element)  :: microbial           !< microbial biomass (bacteria!)
    TYPE(t_slm_struct_sb_element)  :: residue             !< necromass of fungi, microbial and mycorrhiza
    TYPE(t_slm_struct_sb_element)  :: residue_assoc       !< minerally associated microbial necromass 
  CONTAINS
    PROCEDURE  :: init => init_t_slm_struct_sb_pool
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
                                    po4_assoc_fast     , &  !< fast minerally associated phosphate  [mol m-3]
                                    po4_assoc_slow     , &  !< slow minerally associated phosphate  [mol m-3] 
                                    po4_occluded       , &  !< occluded phosphate [mol m-3] 
                                    po4_primary        , &  !< primary phosphate mineral [mol m-3] 
                                    po4_solute              !< phophate in soil solution [mol m-3] 

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
                                    transport_po4_assoc_fast     , &  !< transport of fast minerally associated phosphate  [mol m-3 timestep-1]
                                    transport_po4_assoc_slow     , &  !< transport of slow minerally associated phosphate  [mol m-3 timestep-1] 
                                    transport_po4_occluded       , &  !< transport of occluded phosphate [mol m-3 timestep-1] 
                                    transport_po4_primary        , &  !< transport of primary phosphate [mol m-3 timestep-1] 
                                    transport_po4_solute              !< transport of minerally associated phosphate  [mol m-3 timestep-1]

    ! temperature response functions
    TYPE(t_jsb_var_real3d)          :: &
                                    rtm_decomposition, &    !< rate modifier for decomposition processes [--]
                                    rtm_depolymerisation, & !< rate modifier for depolymerisation processes [--]
                                    rtm_mic_uptake, &       !< rate modifier for microbial uptake [--] 
                                    rtm_sorption, &         !< rate modifier for sorption [--]
                                    rtm_desorption, &       !< rate modifier for desorption [--]
                                    rtm_hsc, &              !< rate modifier for half-saturation constants [--] 
                                    rtm_nitrification, &    !< rate modifier for nitrification [--] 
                                    rtm_denitrification, &  !< rate modifier for denitrification [--] 
                                    rtm_gasdiffusion        !< rate modifier for gasdiffusion [--] 

    ! moisture response functions
    TYPE(t_jsb_var_real3d)          :: &
                                    rmm_decomposition, &    !< rate modifier for decomposition processes [--]
                                    rmm_depolymerisation, & !< rate modifier for depolymerisation processes [--]
                                    rmm_mic_uptake, &       !< rate modifier for microbial uptake [--] 
                                    rmm_sorption, &         !< rate modifier for sorption [--]
                                    rmm_desorption, &       !< rate modifier for desorption [--]
                                    rmm_hsc, &              !< rate modifier for half-saturation constants [--] 
                                    rmm_nitrification, &    !< rate modifier for nitrification [--] 
                                    rmm_gasdiffusion        !< rate modifier for gasdiffusion [--] 

    TYPE(t_jsb_var_real3d)          :: &
                                    het_respiration    ,          &  !< heterotrophic respiration  [mol m-3 s-1]
                                    het_respiration_c13,          &  !< c13 heterotrophic respiration   [mol m-3 s-1]
                                    het_respiration_c14,          &  !< c14 heterotrophic respiration   [mol m-3 s-1]
                                    net_mineralisation_nh4,       &  !< net mineralisation of N    [mol m-3 s-1]
                                    net_mineralisation_nh4_n15,   &  !< net mineralisation of N15  [mol m-3 s-1] 
                                    net_mineralisation_no3,       &  !< net mineralisation of N    [mol m-3 s-1]
                                    net_mineralisation_no3_n15,   &  !< net mineralisation of N15  [mol m-3 s-1] 
                                    net_mineralisation_po4,       &  !< net mineralisation of P  [mol m-3 s-1]
                                    plant_uptake_nh4_sl,          &  !< plant uptake of N    [mol m-3 s-1]
                                    plant_uptake_nh4_n15_sl,      &  !< plant uptake of N15  [mol m-3 s-1] 
                                    plant_uptake_no3_sl,          &  !< plant uptake of N    [mol m-3 s-1]
                                    plant_uptake_no3_n15_sl,      &  !< plant uptake of N15  [mol m-3 s-1] 
                                    plant_uptake_po4_sl,          &  !< plant uptake of P  [mol m-3 s-1]
                                    microbial_uptake_nh4_sl,      &  !< microbial uptake of N    [mol m-3 s-1]
                                    microbial_uptake_nh4_n15_sl,  &  !< microbial uptake of N15  [mol m-3 s-1] 
                                    microbial_uptake_no3_sl,      &  !< microbial uptake of N    [mol m-3 s-1]
                                    microbial_uptake_no3_n15_sl,  &  !< microbial uptake of N15  [mol m-3 s-1] 
                                    microbial_uptake_po4_sl,      &  !< microbial uptake of P  [mol m-3 s-1]
                                    mycorrhiza_uptake_nh4_sl,     &  !< mycorrhiza uptake of N    [mol m-3 s-1]
                                    mycorrhiza_uptake_nh4_n15_sl, &  !< mycorrhiza uptake of N15  [mol m-3 s-1] 
                                    mycorrhiza_uptake_no3_sl,     &  !< mycorrhiza uptake of N    [mol m-3 s-1]
                                    mycorrhiza_uptake_no3_n15_sl, &  !< mycorrhiza uptake of N15  [mol m-3 s-1] 
                                    mycorrhiza_uptake_norg_sl,    &  !< mycorrhiza uptake of N    [mol m-3 s-1]
                                    mycorrhiza_uptake_norg_n15_sl,&  !< mycorrhiza uptake of N15  [mol m-3 s-1] 
                                    mycorrhiza_uptake_po4_sl,   &  !< mycorrhiza uptake of P  [mol m-3 s-1] 
                                    volatilisation_nh4,       &  !< volatilisation of NH4  [µmol m-3 s-1]
                                    volatilisation_nh4_n15, &    !< volatilisation of 15NH4  [µmol m-3 s-1]
                                    nitrification_no3, &         !< nitrification NO3 [µmol m-3 s-1]
                                    nitrification_no3_n15, &     !< nitrification 15NO3 [µmol m-3 s-1]
                                    nitrification_noy, &         !< nitrification to NOy  [µmol m-3 s-1]
                                    nitrification_noy_n15, &     !< nitrification to 15NOy  [µmol m-3 s-1]
                                    nitrification_n2o, &         !< nitrification to N2O [µmol m-3 s-1]
                                    nitrification_n2o_n15, &     !< nitrification to 15N2O [µmol m-3 s-1]
                                    denitrification_noy, &       !< nitrification NOy [µmol m-3 s-1]
                                    denitrification_n2o_n15, &   !< nitrification to 15N2O  [µmol m-3 s-1]
                                    denitrification_n2o, &       !< nitrification to N2O  [µmol m-3 s-1]
                                    denitrification_noy_n15, &   !< nitrification 15NOy [µmol m-3 s-1]
                                    denitrification_n2, &        !< denitrification to N2  [µmol m-3 s1]
                                    denitrification_n2_n15, &    !< denitrification to 15N2  [µmol m-3 s1]
                                    asymb_n_fixation, &          !< asymboitic N fixation [µmol m-3 s-1]
                                    asymb_n_fixation_n15, &      !< asymboitic 15N fixation [µmol m-3 s-1]
                                    weathering_po4, &            !< phosphorus weathering rate [µmol P/m3/s]
                                    occlusion_po4, &             !< phosphorus occlusion rate, assoc_slow --> occluded [µmol P/m3/s]
                                    slow_sorption_po4, &         !< phosphorus desorption rate, assoc_fast <--> assoc_slow [µmol P/m3/s]
                                    fast_adsorption_po4, &       !< phosphorus adsorption rate, solute <--> assoc_fast [µmol P/m3/s]
                                    biochem_mineralisation_po4, & !< phosphorus biochemical mineralization rate, [µmol P/m3/s]
                                    gross_mineralisation_po4     !< phosphorus gross (biological) mineralization rate, [µmol P/m3/s]

    ! 2D transport fluxes for communication with groundwater and atmosphere
    TYPE(t_jsb_var_real2d)          :: & 
                                    leaching_nh4_solute, &      !< land to groundwater flux of NH4 [µmol m-2 s-1]
                                    leaching_no3_solute, &      !< land to groundwater flux of NO3 [µmol m-2 s-1]
                                    leaching_po4_solute, &      !< land to groundwater flux of PO4 [µmol m-2 s-1]
                                    leaching_nh4_n15_solute, &  !< land to groundwater flux of 15NH4 [µmol m-2 s-1]
                                    leaching_no3_n15_solute, &  !< land to groundwater flux of 15NO3 [µmol m-2 s-1]
                                    emission_nh3, &             !< land to atmosphere flux of NH3 [µmol m-2 s-1] 
                                    emission_noy, &             !< land to atmosphere flux of NOy [µmol m-2 s-1] 
                                    emission_n2o, &             !< land to atmosphere flux of N2O [µmol m-2 s-1]
                                    emission_n2, &              !< land to atmosphere flux of N2 [µmol m-2 s-1]
                                    emission_nh3_n15, &         !< land to atmosphere flux of 15NH3 [µmol m-2 s-1]
                                    emission_noy_n15, &         !< land to atmosphere flux of 15NOy [µmol m-2 s-1]
                                    emission_n2o_n15, &         !< land to atmosphere flux of 15N2O [µmol m-2 s-1]
                                    emission_n2_n15             !< land to atmosphere flux of 15N2 [µmol m-2 s-1]
                       
    ! 3D 
    TYPE(t_jsb_var_real3d)          :: &
                                    qmax_org, &                   !< sorption capacity of dom of mineral soil [kg C kg-1 mineral soil]
                                    qmax_po4, &                   !< sorption capacity of po4 of mineral soil [kg C kg-1 mineral soil]
                                    k_bioturb, &                  !< diffusion rate resulting from bioturbation
                                    bulk_dens_corr_sl, &          !< organic fraction corrected bulk-density of the soil [kg m-3]
                                    bulk_soil_carbon_sl, &        !< bulk soil density of carbon [kg m-3] 
                                    particle_fluxrate, &          !< burial rate due to organic matter fall and decay [s-1?] 
                                    anaerobic_volume_fraction_sl  !< anaerobic volume fraction [--]

    ! 3D soil properties
    TYPE(t_jsb_var_real3d)          :: &
                                    ph_sl               , &       !< ph per soil layer
                                    lithology_class_sl            !< ..

    ! 2D total fluxes for mass balance checks in soil_biogeochemistry model
    TYPE(t_jsb_var_real2d)          :: & 
                                    total_flux_carbon, &          !< carbon mass balance checks in soil_biogeochemistry model
                                    total_flux_nitrogen, &        !< nitrogen mass balance checks in soil_biogeochemistry model
                                    total_flux_phosphorus, &      !< phosphorus mass balance checks in soil_biogeochemistry model
                                    total_flux_carbon13, &        !< carbon13 mass balance checks in soil_biogeochemistry model
                                    total_flux_carbon14, &        !< carbon14 mass balance checks in soil_biogeochemistry model
                                    total_flux_nitrogen15         !< nitrogen15 mass balance checks in soil_biogeochemistry model

   TYPE(t_slm_struct_sb_element)    :: &   ! simple structure for elements: t_jsb_struct
                                    delta_pools, &              !< for mass balance checks in soil_biogeochemistry model
                                    mycorrhiza_export           !< nutrient flux from mycorrhiza to host plant


    TYPE(t_slm_struct_sb_pool)      :: &   ! nested pool TYPE t_jsb_n_struct
                                    pool, &                     !< pools of the sb
                                    formation, &
                                    loss, &
                                    transport



  CONTAINS
    PROCEDURE :: Init => Init_sb_memory
  END TYPE t_sb_memory


  ! definition of sb_ pool and element structures with ALLOCATABLE variables
  ! allocation of variables done with Init_sb_pools_elements_a during model init in init_slm()
  TYPE(t_sb_pool_a)              :: sb_pool             , &
                                    sb_formation        , &
                                    sb_loss             , &
                                    sb_transport
  TYPE(t_sb_element_a)           :: sb_delta_pools      , &
                                    sb_mycorrhiza_export


  CHARACTER(len=*), PARAMETER :: modname = 'mo_sb_memory_class'


CONTAINS

  !-----------------------------------------------------------------------------------------------------!
  ! t_sb_element_a operator overload                                                                   !
  !-----------------------------------------------------------------------------------------------------!
  !-----------------------------------------------------------------------------------------------------
  !> + (t_sb_element_a + t_sb_element_a)
  !!
  !-----------------------------------------------------------------------------------------------------
  PURE FUNCTION Add_element_a(this,element2) RESULT(element3)
    IMPLICIT NONE
    ! 0.1 InOut
    CLASS(t_sb_element_a), INTENT(in)  :: this, element2
    TYPE(t_sb_element_a)               :: element3
    ! 1.0 go
    ALLOCATE(element3%carbon,     source=this%carbon)       ! note: Function results are automatically deallocated after the result is "used" in the outer scope.
    ALLOCATE(element3%nitrogen,   source=this%nitrogen)
    ALLOCATE(element3%phosphorus, source=this%phosphorus)
    ALLOCATE(element3%carbon13,   source=this%carbon13)
    ALLOCATE(element3%carbon14,   source=this%carbon14)
    ALLOCATE(element3%nitrogen15, source=this%nitrogen15)
    element3%carbon(:,:)       = this%carbon(:,:)      + element2%carbon(:,:)
    element3%nitrogen(:,:)     = this%nitrogen(:,:)    + element2%nitrogen(:,:)
    element3%phosphorus(:,:)   = this%phosphorus(:,:)  + element2%phosphorus(:,:)
    element3%carbon13(:,:)     = this%carbon13(:,:)    + element2%carbon13(:,:)
    element3%carbon14(:,:)     = this%carbon14(:,:)    + element2%carbon14(:,:)
    element3%nitrogen15(:,:)   = this%nitrogen15(:,:)  + element2%nitrogen15(:,:)
  END FUNCTION Add_element_a
  !-----------------------------------------------------------------------------------------------------
  !> + (t_sb_element_a + scalar)
  !!
  !-----------------------------------------------------------------------------------------------------
  PURE FUNCTION Add_element_a_and_scalar(this,scalar) RESULT(element3)
    IMPLICIT NONE
    ! 0.1 InOut
    CLASS(t_sb_element_a), INTENT(in)  :: this
    REAL(wp),              INTENT(in)  :: scalar
    TYPE(t_sb_element_a)               :: element3
    ! 1.0 go
    ALLOCATE(element3%carbon,     source=this%carbon)
    ALLOCATE(element3%nitrogen,   source=this%nitrogen)
    ALLOCATE(element3%phosphorus, source=this%phosphorus)
    ALLOCATE(element3%carbon13,   source=this%carbon13)
    ALLOCATE(element3%carbon14,   source=this%carbon14)
    ALLOCATE(element3%nitrogen15, source=this%nitrogen15)
    element3%carbon(:,:)       = this%carbon(:,:)      + scalar
    element3%nitrogen(:,:)     = this%nitrogen(:,:)    + scalar
    element3%phosphorus(:,:)   = this%phosphorus(:,:)  + scalar
    element3%carbon13(:,:)     = this%carbon13(:,:)    + scalar
    element3%carbon14(:,:)     = this%carbon14(:,:)    + scalar
    element3%nitrogen15(:,:)   = this%nitrogen15(:,:)  + scalar
  END FUNCTION Add_element_a_and_scalar
  !-----------------------------------------------------------------------------------------------------
  !> + (t_sb_element_a + field2d)
  !!
  !-----------------------------------------------------------------------------------------------------
  PURE FUNCTION Add_element_a_and_field2d(this,field2d) RESULT(element3)
    IMPLICIT NONE
    ! 0.1 InOut
    CLASS(t_sb_element_a), INTENT(in)  :: this
    REAL(wp),              INTENT(in)  :: field2d(:,:)
    TYPE(t_sb_element_a)               :: element3
    ! 1.0 go
    ALLOCATE(element3%carbon,     source=this%carbon)
    ALLOCATE(element3%nitrogen,   source=this%nitrogen)
    ALLOCATE(element3%phosphorus, source=this%phosphorus)
    ALLOCATE(element3%carbon13,   source=this%carbon13)
    ALLOCATE(element3%carbon14,   source=this%carbon14)
    ALLOCATE(element3%nitrogen15, source=this%nitrogen15)
    element3%carbon(:,:)       = this%carbon(:,:)      + field2d(:,:)
    element3%nitrogen(:,:)     = this%nitrogen(:,:)    + field2d(:,:)
    element3%phosphorus(:,:)   = this%phosphorus(:,:)  + field2d(:,:)
    element3%carbon13(:,:)     = this%carbon13(:,:)    + field2d(:,:)
    element3%carbon14(:,:)     = this%carbon14(:,:)    + field2d(:,:)
    element3%nitrogen15(:,:)   = this%nitrogen15(:,:)  + field2d(:,:)
  END FUNCTION Add_element_a_and_field2d
  !-----------------------------------------------------------------------------------------------------
  !> - (t_sb_element_a - t_sb_element_a)
  !!
  !-----------------------------------------------------------------------------------------------------
  PURE FUNCTION Subtract_element_a(this,element2) RESULT(element3)
    IMPLICIT NONE
    ! 0.1 InOut
    CLASS(t_sb_element_a), INTENT(in)  :: this, element2
    TYPE(t_sb_element_a)               :: element3
    ! 1.0 go
    ALLOCATE(element3%carbon,     source=this%carbon)
    ALLOCATE(element3%nitrogen,   source=this%nitrogen)
    ALLOCATE(element3%phosphorus, source=this%phosphorus)
    ALLOCATE(element3%carbon13,   source=this%carbon13)
    ALLOCATE(element3%carbon14,   source=this%carbon14)
    ALLOCATE(element3%nitrogen15, source=this%nitrogen15)
    element3%carbon(:,:)       = this%carbon(:,:)      - element2%carbon(:,:)
    element3%nitrogen(:,:)     = this%nitrogen(:,:)    - element2%nitrogen(:,:)
    element3%phosphorus(:,:)   = this%phosphorus(:,:)  - element2%phosphorus(:,:)
    element3%carbon13(:,:)     = this%carbon13(:,:)    - element2%carbon13(:,:)
    element3%carbon14(:,:)     = this%carbon14(:,:)    - element2%carbon14(:,:)
    element3%nitrogen15(:,:)   = this%nitrogen15(:,:)  - element2%nitrogen15(:,:)
  END FUNCTION Subtract_element_a
  !-----------------------------------------------------------------------------------------------------
  !> - (t_sb_element_a - scalar)
  !!
  !-----------------------------------------------------------------------------------------------------
  PURE FUNCTION Subtract_element_a_and_scalar(this,scalar) RESULT(element3)
    IMPLICIT NONE
    ! 0.1 InOut
    CLASS(t_sb_element_a), INTENT(in)  :: this
    REAL(wp),              INTENT(in)  :: scalar
    TYPE(t_sb_element_a)               :: element3
    ! 1.0 go
    ALLOCATE(element3%carbon,     source=this%carbon)
    ALLOCATE(element3%nitrogen,   source=this%nitrogen)
    ALLOCATE(element3%phosphorus, source=this%phosphorus)
    ALLOCATE(element3%carbon13,   source=this%carbon13)
    ALLOCATE(element3%carbon14,   source=this%carbon14)
    ALLOCATE(element3%nitrogen15, source=this%nitrogen15)
    element3%carbon(:,:)       = this%carbon(:,:)      - scalar
    element3%nitrogen(:,:)     = this%nitrogen(:,:)    - scalar
    element3%phosphorus(:,:)   = this%phosphorus(:,:)  - scalar
    element3%carbon13(:,:)     = this%carbon13(:,:)    - scalar
    element3%carbon14(:,:)     = this%carbon14(:,:)    - scalar
    element3%nitrogen15(:,:)   = this%nitrogen15(:,:)  - scalar
  END FUNCTION Subtract_element_a_and_scalar
  !-----------------------------------------------------------------------------------------------------
  !> - (t_sb_element_a - field2d)
  !!
  !-----------------------------------------------------------------------------------------------------
  PURE FUNCTION Subtract_element_a_and_field2d(this,field2d) RESULT(element3)
    IMPLICIT NONE
    ! 0.1 InOut
    CLASS(t_sb_element_a), INTENT(in)  :: this
    REAL(wp),              INTENT(in)  :: field2d(:,:)
    TYPE(t_sb_element_a)               :: element3
    ! 1.0 go
    ALLOCATE(element3%carbon,     source=this%carbon)
    ALLOCATE(element3%nitrogen,   source=this%nitrogen)
    ALLOCATE(element3%phosphorus, source=this%phosphorus)
    ALLOCATE(element3%carbon13,   source=this%carbon13)
    ALLOCATE(element3%carbon14,   source=this%carbon14)
    ALLOCATE(element3%nitrogen15, source=this%nitrogen15)
    element3%carbon(:,:)       = this%carbon(:,:)      - field2d(:,:)
    element3%nitrogen(:,:)     = this%nitrogen(:,:)    - field2d(:,:)
    element3%phosphorus(:,:)   = this%phosphorus(:,:)  - field2d(:,:)
    element3%carbon13(:,:)     = this%carbon13(:,:)    - field2d(:,:)
    element3%carbon14(:,:)     = this%carbon14(:,:)    - field2d(:,:)
    element3%nitrogen15(:,:)   = this%nitrogen15(:,:)  - field2d(:,:)
  END FUNCTION Subtract_element_a_and_field2d
  !-----------------------------------------------------------------------------------------------------
  !> / (t_sb_element_a / t_sb_element_a)
  !!
  !-----------------------------------------------------------------------------------------------------
  PURE FUNCTION Divide_element_a(this,element2) RESULT(element3)
    IMPLICIT NONE
    ! 0.1 InOut
    CLASS(t_sb_element_a), INTENT(in)  :: this, element2
    TYPE(t_sb_element_a)               :: element3
    ! 1.0 go
    ALLOCATE(element3%carbon,     source=this%carbon)       ! note: Function results are automatically deallocated after the result is "used" in the outer scope.
    ALLOCATE(element3%nitrogen,   source=this%nitrogen)
    ALLOCATE(element3%phosphorus, source=this%phosphorus)
    ALLOCATE(element3%carbon13,   source=this%carbon13)
    ALLOCATE(element3%carbon14,   source=this%carbon14)
    ALLOCATE(element3%nitrogen15, source=this%nitrogen15)
    element3%carbon(:,:)       = this%carbon(:,:)      / element2%carbon(:,:)
    element3%nitrogen(:,:)     = this%nitrogen(:,:)    / element2%nitrogen(:,:)
    element3%phosphorus(:,:)   = this%phosphorus(:,:)  / element2%phosphorus(:,:)
    element3%carbon13(:,:)     = this%carbon13(:,:)    / element2%carbon13(:,:)
    element3%carbon14(:,:)     = this%carbon14(:,:)    / element2%carbon14(:,:)
    element3%nitrogen15(:,:)   = this%nitrogen15(:,:)  / element2%nitrogen15(:,:)
  END FUNCTION Divide_element_a
  !-----------------------------------------------------------------------------------------------------
  !> / (t_sb_element_a / scalar)
  !!
  !-----------------------------------------------------------------------------------------------------
  PURE FUNCTION Divide_element_a_and_scalar(this,scalar) RESULT(element3)
    IMPLICIT NONE
    ! 0.1 InOut
    CLASS(t_sb_element_a), INTENT(in)  :: this
    REAL(wp),              INTENT(in)  :: scalar
    TYPE(t_sb_element_a)               :: element3
    ! 1.0 go
    ALLOCATE(element3%carbon,     source=this%carbon)
    ALLOCATE(element3%nitrogen,   source=this%nitrogen)
    ALLOCATE(element3%phosphorus, source=this%phosphorus)
    ALLOCATE(element3%carbon13,   source=this%carbon13)
    ALLOCATE(element3%carbon14,   source=this%carbon14)
    ALLOCATE(element3%nitrogen15, source=this%nitrogen15)
    element3%carbon(:,:)       = this%carbon(:,:)      / scalar
    element3%nitrogen(:,:)     = this%nitrogen(:,:)    / scalar
    element3%phosphorus(:,:)   = this%phosphorus(:,:)  / scalar
    element3%carbon13(:,:)     = this%carbon13(:,:)    / scalar
    element3%carbon14(:,:)     = this%carbon14(:,:)    / scalar
    element3%nitrogen15(:,:)   = this%nitrogen15(:,:)  / scalar
  END FUNCTION Divide_element_a_and_scalar
  !-----------------------------------------------------------------------------------------------------
  !> / (t_sb_element_a / field2d)
  !!
  !-----------------------------------------------------------------------------------------------------
  PURE FUNCTION Divide_element_a_and_field2d(this,field2d) RESULT(element3)
    IMPLICIT NONE
    ! 0.1 InOut
    CLASS(t_sb_element_a), INTENT(in)  :: this
    REAL(wp),              INTENT(in)  :: field2d(:,:)
    TYPE(t_sb_element_a)               :: element3
    ! 1.0 go
    ALLOCATE(element3%carbon,     source=this%carbon)
    ALLOCATE(element3%nitrogen,   source=this%nitrogen)
    ALLOCATE(element3%phosphorus, source=this%phosphorus)
    ALLOCATE(element3%carbon13,   source=this%carbon13)
    ALLOCATE(element3%carbon14,   source=this%carbon14)
    ALLOCATE(element3%nitrogen15, source=this%nitrogen15)
    element3%carbon(:,:)       = this%carbon(:,:)      / field2d(:,:)
    element3%nitrogen(:,:)     = this%nitrogen(:,:)    / field2d(:,:)
    element3%phosphorus(:,:)   = this%phosphorus(:,:)  / field2d(:,:)
    element3%carbon13(:,:)     = this%carbon13(:,:)    / field2d(:,:)
    element3%carbon14(:,:)     = this%carbon14(:,:)    / field2d(:,:)
    element3%nitrogen15(:,:)   = this%nitrogen15(:,:)  / field2d(:,:)
  END FUNCTION Divide_element_a_and_field2d
  !-----------------------------------------------------------------------------------------------------
  !> * (t_sb_element_a * t_sb_element_a)
  !!
  !-----------------------------------------------------------------------------------------------------
  PURE FUNCTION Multiply_element_a(this,element2) RESULT(element3)
    IMPLICIT NONE
    ! 0.1 InOut
    CLASS(t_sb_element_a), INTENT(in)  :: this, element2
    TYPE(t_sb_element_a)               :: element3
    ! 1.0 go
    ALLOCATE(element3%carbon,     source=this%carbon)
    ALLOCATE(element3%nitrogen,   source=this%nitrogen)
    ALLOCATE(element3%phosphorus, source=this%phosphorus)
    ALLOCATE(element3%carbon13,   source=this%carbon13)
    ALLOCATE(element3%carbon14,   source=this%carbon14)
    ALLOCATE(element3%nitrogen15, source=this%nitrogen15)
    element3%carbon(:,:)       = this%carbon(:,:)      * element2%carbon(:,:)
    element3%nitrogen(:,:)     = this%nitrogen(:,:)    * element2%nitrogen(:,:)
    element3%phosphorus(:,:)   = this%phosphorus(:,:)  * element2%phosphorus(:,:)
    element3%carbon13(:,:)     = this%carbon13(:,:)    * element2%carbon13(:,:)
    element3%carbon14(:,:)     = this%carbon14(:,:)    * element2%carbon14(:,:)
    element3%nitrogen15(:,:)   = this%nitrogen15(:,:)  * element2%nitrogen15(:,:)
  END FUNCTION Multiply_element_a
  !-----------------------------------------------------------------------------------------------------
  !> * (t_sb_element_a * scalar)
  !!
  !-----------------------------------------------------------------------------------------------------
  PURE FUNCTION Multiply_element_a_and_scalar(this,scalar) RESULT(element3)
    IMPLICIT NONE
    ! 0.1 InOut
    CLASS(t_sb_element_a), INTENT(in)  :: this
    REAL(wp),              INTENT(in)  :: scalar
    TYPE(t_sb_element_a)               :: element3
    ! 1.0 go
    ALLOCATE(element3%carbon,     source=this%carbon)
    ALLOCATE(element3%nitrogen,   source=this%nitrogen)
    ALLOCATE(element3%phosphorus, source=this%phosphorus)
    ALLOCATE(element3%carbon13,   source=this%carbon13)
    ALLOCATE(element3%carbon14,   source=this%carbon14)
    ALLOCATE(element3%nitrogen15, source=this%nitrogen15)
    element3%carbon(:,:)       = this%carbon(:,:)      * scalar
    element3%nitrogen(:,:)     = this%nitrogen(:,:)    * scalar
    element3%phosphorus(:,:)   = this%phosphorus(:,:)  * scalar
    element3%carbon13(:,:)     = this%carbon13(:,:)    * scalar
    element3%carbon14(:,:)     = this%carbon14(:,:)    * scalar
    element3%nitrogen15(:,:)   = this%nitrogen15(:,:)  * scalar
  END FUNCTION Multiply_element_a_and_scalar
  !-----------------------------------------------------------------------------------------------------
  !> * (scalar * t_sb_element_a)
  !!
  !-----------------------------------------------------------------------------------------------------
  PURE FUNCTION Multiply_scalar_and_element_a(scalar,this2) RESULT(element3)
    IMPLICIT NONE
    ! 0.1 InOut
    REAL(wp),              INTENT(in)  :: scalar
    CLASS(t_sb_element_a), INTENT(in)  :: this2
    TYPE(t_sb_element_a)               :: element3
    ! 1.0 go
    ALLOCATE(element3%carbon,     source=this2%carbon)
    ALLOCATE(element3%nitrogen,   source=this2%nitrogen)
    ALLOCATE(element3%phosphorus, source=this2%phosphorus)
    ALLOCATE(element3%carbon13,   source=this2%carbon13)
    ALLOCATE(element3%carbon14,   source=this2%carbon14)
    ALLOCATE(element3%nitrogen15, source=this2%nitrogen15)
    element3%carbon(:,:)       = scalar * this2%carbon(:,:)
    element3%nitrogen(:,:)     = scalar * this2%nitrogen(:,:)
    element3%phosphorus(:,:)   = scalar * this2%phosphorus(:,:)
    element3%carbon13(:,:)     = scalar * this2%carbon13(:,:)
    element3%carbon14(:,:)     = scalar * this2%carbon14(:,:)
    element3%nitrogen15(:,:)   = scalar * this2%nitrogen15(:,:)
  END FUNCTION Multiply_scalar_and_element_a
  !-----------------------------------------------------------------------------------------------------
  !> * (t_sb_element_a * field2d)
  !!
  !-----------------------------------------------------------------------------------------------------
  PURE FUNCTION Multiply_element_a_and_field2d(this,field2d) RESULT(element3)
    IMPLICIT NONE
    ! 0.1 InOut
    CLASS(t_sb_element_a), INTENT(in)  :: this
    REAL(wp),              INTENT(in)  :: field2d(:,:)
    TYPE(t_sb_element_a)               :: element3
    ! 1.0 go
    ALLOCATE(element3%carbon,     source=this%carbon)
    ALLOCATE(element3%nitrogen,   source=this%nitrogen)
    ALLOCATE(element3%phosphorus, source=this%phosphorus)
    ALLOCATE(element3%carbon13,   source=this%carbon13)
    ALLOCATE(element3%carbon14,   source=this%carbon14)
    ALLOCATE(element3%nitrogen15, source=this%nitrogen15)
    element3%carbon(:,:)       = this%carbon(:,:)      * field2d(:,:)
    element3%nitrogen(:,:)     = this%nitrogen(:,:)    * field2d(:,:)
    element3%phosphorus(:,:)   = this%phosphorus(:,:)  * field2d(:,:)
    element3%carbon13(:,:)     = this%carbon13(:,:)    * field2d(:,:)
    element3%carbon14(:,:)     = this%carbon14(:,:)    * field2d(:,:)
    element3%nitrogen15(:,:)   = this%nitrogen15(:,:)  * field2d(:,:)
  END FUNCTION Multiply_element_a_and_field2d
  !-----------------------------------------------------------------------------------------------------
  !> * (field2d * t_sb_element_a)
  !!
  !-----------------------------------------------------------------------------------------------------
  PURE FUNCTION Multiply_field2d_and_element_a(field2d,this2) RESULT(element3)
    IMPLICIT NONE
    ! 0.1 InOut
    REAL(wp),              INTENT(in)  :: field2d(:,:)
    CLASS(t_sb_element_a), INTENT(in)  :: this2
    TYPE(t_sb_element_a)               :: element3
    ! 1.0 go
    ALLOCATE(element3%carbon,     source=this2%carbon)
    ALLOCATE(element3%nitrogen,   source=this2%nitrogen)
    ALLOCATE(element3%phosphorus, source=this2%phosphorus)
    ALLOCATE(element3%carbon13,   source=this2%carbon13)
    ALLOCATE(element3%carbon14,   source=this2%carbon14)
    ALLOCATE(element3%nitrogen15, source=this2%nitrogen15)
    element3%carbon(:,:)       = field2d(:,:) * this2%carbon(:,:)
    element3%nitrogen(:,:)     = field2d(:,:) * this2%nitrogen(:,:)
    element3%phosphorus(:,:)   = field2d(:,:) * this2%phosphorus(:,:)
    element3%carbon13(:,:)     = field2d(:,:) * this2%carbon13(:,:)
    element3%carbon14(:,:)     = field2d(:,:) * this2%carbon14(:,:)
    element3%nitrogen15(:,:)   = field2d(:,:) * this2%nitrogen15(:,:)
  END FUNCTION Multiply_field2d_and_element_a
  

  !-----------------------------------------------------------------------------------------------------!
  ! t_sb_pool_a operator overload                                                                       !
  !-----------------------------------------------------------------------------------------------------!
  !-----------------------------------------------------------------------------------------------------
  !> + (t_sb_pool_a + t_sb_pool_a)
  !!
  !-----------------------------------------------------------------------------------------------------
  PURE FUNCTION Add_pool_a(this,pool2) RESULT(pool3)
    IMPLICIT NONE
    ! 0.1 InOut
    CLASS(t_sb_pool_a), INTENT(in) :: this, pool2
    TYPE(t_sb_pool_a)              :: pool3
    ! 1.0 go
    pool3%dom               = this%dom              + pool2%dom
    pool3%dom_assoc         = this%dom_assoc        + pool2%dom_assoc
    pool3%soluable_litter   = this%soluable_litter  + pool2%soluable_litter
    pool3%polymeric_litter  = this%polymeric_litter + pool2%polymeric_litter
    pool3%woody_litter      = this%woody_litter     + pool2%woody_litter
    pool3%fungi             = this%fungi            + pool2%fungi
    pool3%mycorrhiza        = this%mycorrhiza       + pool2%mycorrhiza
    pool3%microbial         = this%microbial        + pool2%microbial
    pool3%residue           = this%residue          + pool2%residue
    pool3%residue_assoc     = this%residue_assoc    + pool2%residue_assoc
  END FUNCTION Add_pool_a
  !-----------------------------------------------------------------------------------------------------
  !> + (t_sb_pool_a + scalar)
  !!
  !-----------------------------------------------------------------------------------------------------
  PURE FUNCTION Add_pool_a_and_scalar(this,scalar) RESULT(pool3)
    IMPLICIT NONE
    ! 0.1 InOut
    CLASS(t_sb_pool_a), INTENT(in) :: this
    REAL(wp),           INTENT(in) :: scalar
    TYPE(t_sb_pool_a)              :: pool3
    ! 1.0 go
    pool3%dom               = this%dom              + scalar
    pool3%dom_assoc         = this%dom_assoc        + scalar
    pool3%soluable_litter   = this%soluable_litter  + scalar
    pool3%polymeric_litter  = this%polymeric_litter + scalar
    pool3%woody_litter      = this%woody_litter     + scalar
    pool3%fungi             = this%fungi            + scalar
    pool3%mycorrhiza        = this%mycorrhiza       + scalar
    pool3%microbial         = this%microbial        + scalar
    pool3%residue           = this%residue          + scalar
    pool3%residue_assoc     = this%residue_assoc    + scalar
  END FUNCTION Add_pool_a_and_scalar
  !-----------------------------------------------------------------------------------------------------
  !> + (t_sb_pool_a + field2d)
  !!
  !-----------------------------------------------------------------------------------------------------
  PURE FUNCTION Add_pool_a_and_field2d(this,field2d) RESULT(pool3)
    IMPLICIT NONE
    ! 0.1 InOut
    CLASS(t_sb_pool_a), INTENT(in) :: this
    REAL(wp),           INTENT(in) :: field2d(:,:)
    TYPE(t_sb_pool_a)              :: pool3
    ! 1.0 go
    pool3%dom               = this%dom              + field2d(:,:)
    pool3%dom_assoc         = this%dom_assoc        + field2d(:,:)
    pool3%soluable_litter   = this%soluable_litter  + field2d(:,:)
    pool3%polymeric_litter  = this%polymeric_litter + field2d(:,:)
    pool3%woody_litter      = this%woody_litter     + field2d(:,:)
    pool3%fungi             = this%fungi            + field2d(:,:)
    pool3%mycorrhiza        = this%mycorrhiza       + field2d(:,:)
    pool3%microbial         = this%microbial        + field2d(:,:)
    pool3%residue           = this%residue          + field2d(:,:)
    pool3%residue_assoc     = this%residue_assoc    + field2d(:,:)
  END FUNCTION Add_pool_a_and_field2d
  !-----------------------------------------------------------------------------------------------------
  !> - (t_sb_pool_a - t_sb_pool_a)
  !!
  !-----------------------------------------------------------------------------------------------------
  PURE FUNCTION Subtract_pool_a(this,pool2) RESULT(pool3)
    IMPLICIT NONE
    ! 0.1 InOut
    CLASS(t_sb_pool_a), INTENT(in) :: this, pool2
    TYPE(t_sb_pool_a)              :: pool3
    ! 1.0 go
    pool3%dom               = this%dom              - pool2%dom
    pool3%dom_assoc         = this%dom_assoc        - pool2%dom_assoc
    pool3%soluable_litter   = this%soluable_litter  - pool2%soluable_litter
    pool3%polymeric_litter  = this%polymeric_litter - pool2%polymeric_litter
    pool3%woody_litter      = this%woody_litter     - pool2%woody_litter
    pool3%fungi             = this%fungi            - pool2%fungi
    pool3%mycorrhiza        = this%mycorrhiza       - pool2%mycorrhiza
    pool3%microbial         = this%microbial        - pool2%microbial
    pool3%residue           = this%residue          - pool2%residue
    pool3%residue_assoc     = this%residue_assoc    - pool2%residue_assoc
  END FUNCTION Subtract_pool_a
  !-----------------------------------------------------------------------------------------------------
  !> - (t_sb_pool_a - scalar)
  !!
  !-----------------------------------------------------------------------------------------------------
  PURE FUNCTION Subtract_pool_a_and_scalar(this,scalar) RESULT(pool3)
    IMPLICIT NONE
    ! 0.1 InOut
    CLASS(t_sb_pool_a), INTENT(in) :: this
    REAL(wp),           INTENT(in) :: scalar
    TYPE(t_sb_pool_a)              :: pool3
    ! 1.0 go
    pool3%dom               = this%dom              - scalar
    pool3%dom_assoc         = this%dom_assoc        - scalar
    pool3%soluable_litter   = this%soluable_litter  - scalar
    pool3%polymeric_litter  = this%polymeric_litter - scalar
    pool3%woody_litter      = this%woody_litter     - scalar
    pool3%fungi             = this%fungi            - scalar
    pool3%mycorrhiza        = this%mycorrhiza       - scalar
    pool3%microbial         = this%microbial        - scalar
    pool3%residue           = this%residue          - scalar
    pool3%residue_assoc     = this%residue_assoc    - scalar
  END FUNCTION Subtract_pool_a_and_scalar
  !-----------------------------------------------------------------------------------------------------
  !> - (t_sb_pool_a - field2d)
  !!
  !-----------------------------------------------------------------------------------------------------
  PURE FUNCTION Subtract_pool_a_and_field2d(this,field2d) RESULT(pool3)
    IMPLICIT NONE
    ! 0.1 InOut
    CLASS(t_sb_pool_a), INTENT(in) :: this
    REAL(wp),           INTENT(in) :: field2d(:,:)
    TYPE(t_sb_pool_a)              :: pool3
    ! 1.0 go
    pool3%dom               = this%dom              - field2d(:,:)
    pool3%dom_assoc         = this%dom_assoc        - field2d(:,:)
    pool3%soluable_litter   = this%soluable_litter  - field2d(:,:)
    pool3%polymeric_litter  = this%polymeric_litter - field2d(:,:)
    pool3%woody_litter      = this%woody_litter     - field2d(:,:)
    pool3%fungi             = this%fungi            - field2d(:,:)
    pool3%mycorrhiza        = this%mycorrhiza       - field2d(:,:)
    pool3%microbial         = this%microbial        - field2d(:,:)
    pool3%residue           = this%residue          - field2d(:,:)
    pool3%residue_assoc     = this%residue_assoc    - field2d(:,:)
  END FUNCTION Subtract_pool_a_and_field2d
  !-----------------------------------------------------------------------------------------------------
  !> / (t_sb_pool_a / t_sb_pool_a)
  !!
  !-----------------------------------------------------------------------------------------------------
  PURE FUNCTION Divide_pool_a(this,pool2) RESULT(pool3)
    IMPLICIT NONE
    ! 0.1 InOut
    CLASS(t_sb_pool_a), INTENT(in) :: this, pool2
    TYPE(t_sb_pool_a)              :: pool3
    ! 1.0 go
    pool3%dom               = this%dom              / pool2%dom
    pool3%dom_assoc         = this%dom_assoc        / pool2%dom_assoc
    pool3%soluable_litter   = this%soluable_litter  / pool2%soluable_litter
    pool3%polymeric_litter  = this%polymeric_litter / pool2%polymeric_litter
    pool3%woody_litter      = this%woody_litter     / pool2%woody_litter
    pool3%fungi             = this%fungi            / pool2%fungi
    pool3%mycorrhiza        = this%mycorrhiza       / pool2%mycorrhiza
    pool3%microbial         = this%microbial        / pool2%microbial
    pool3%residue           = this%residue          / pool2%residue
    pool3%residue_assoc     = this%residue_assoc    / pool2%residue_assoc
  END FUNCTION Divide_pool_a
  !-----------------------------------------------------------------------------------------------------
  !> / (t_sb_pool_a / scalar)
  !!
  !-----------------------------------------------------------------------------------------------------
  PURE FUNCTION Divide_pool_a_and_scalar(this,scalar) RESULT(pool3)
    IMPLICIT NONE
    ! 0.1 InOut
    CLASS(t_sb_pool_a), INTENT(in) :: this
    REAL(wp),           INTENT(in) :: scalar
    TYPE(t_sb_pool_a)              :: pool3
    ! 1.0 go
    pool3%dom               = this%dom              / scalar
    pool3%dom_assoc         = this%dom_assoc        / scalar
    pool3%soluable_litter   = this%soluable_litter  / scalar
    pool3%polymeric_litter  = this%polymeric_litter / scalar
    pool3%woody_litter      = this%woody_litter     / scalar
    pool3%fungi             = this%fungi            / scalar
    pool3%mycorrhiza        = this%mycorrhiza       / scalar
    pool3%microbial         = this%microbial        / scalar
    pool3%residue           = this%residue          / scalar
    pool3%residue_assoc     = this%residue_assoc    / scalar
  END FUNCTION Divide_pool_a_and_scalar
  !-----------------------------------------------------------------------------------------------------
  !> / (t_sb_pool_a / field2d)
  !!
  !-----------------------------------------------------------------------------------------------------
  PURE FUNCTION Divide_pool_a_and_field2d(this,field2d) RESULT(pool3)
    IMPLICIT NONE
    ! 0.1 InOut
    CLASS(t_sb_pool_a), INTENT(in) :: this
    REAL(wp),           INTENT(in) :: field2d(:,:)
    TYPE(t_sb_pool_a)              :: pool3
    ! 1.0 go
    pool3%dom               = this%dom              / field2d(:,:)
    pool3%dom_assoc         = this%dom_assoc        / field2d(:,:)
    pool3%soluable_litter   = this%soluable_litter  / field2d(:,:)
    pool3%polymeric_litter  = this%polymeric_litter / field2d(:,:)
    pool3%woody_litter      = this%woody_litter     / field2d(:,:)
    pool3%fungi             = this%fungi            / field2d(:,:)
    pool3%mycorrhiza        = this%mycorrhiza       / field2d(:,:)
    pool3%microbial         = this%microbial        / field2d(:,:)
    pool3%residue           = this%residue          / field2d(:,:)
    pool3%residue_assoc     = this%residue_assoc    / field2d(:,:)
  END FUNCTION Divide_pool_a_and_field2d
  !-----------------------------------------------------------------------------------------------------
  !> * (t_sb_pool_a * t_sb_pool_a)
  !!
  !-----------------------------------------------------------------------------------------------------
  PURE FUNCTION Multiply_pool_a(this,pool2) RESULT(pool3)
    IMPLICIT NONE
    ! 0.1 InOut
    CLASS(t_sb_pool_a), INTENT(in) :: this, pool2
    TYPE(t_sb_pool_a)              :: pool3
    ! 1.0 go
    pool3%dom               = this%dom              * pool2%dom
    pool3%dom_assoc         = this%dom_assoc        * pool2%dom_assoc
    pool3%soluable_litter   = this%soluable_litter  * pool2%soluable_litter
    pool3%polymeric_litter  = this%polymeric_litter * pool2%polymeric_litter
    pool3%woody_litter      = this%woody_litter     * pool2%woody_litter
    pool3%fungi             = this%fungi            * pool2%fungi
    pool3%mycorrhiza        = this%mycorrhiza       * pool2%mycorrhiza
    pool3%microbial         = this%microbial        * pool2%microbial
    pool3%residue           = this%residue          * pool2%residue
    pool3%residue_assoc     = this%residue_assoc    * pool2%residue_assoc
  END FUNCTION Multiply_pool_a
  !-----------------------------------------------------------------------------------------------------
  !> * (t_sb_pool_a * scalar)
  !!
  !-----------------------------------------------------------------------------------------------------
  PURE FUNCTION Multiply_pool_a_and_scalar(this,scalar) RESULT(pool3)
    IMPLICIT NONE
    ! 0.1 InOut
    CLASS(t_sb_pool_a), INTENT(in) :: this
    REAL(wp),           INTENT(in) :: scalar
    TYPE(t_sb_pool_a)              :: pool3
    ! 1.0 go
    pool3%dom               = this%dom              * scalar
    pool3%dom_assoc         = this%dom_assoc        * scalar
    pool3%soluable_litter   = this%soluable_litter  * scalar
    pool3%polymeric_litter  = this%polymeric_litter * scalar
    pool3%woody_litter      = this%woody_litter     * scalar
    pool3%fungi             = this%fungi            * scalar
    pool3%mycorrhiza        = this%mycorrhiza       * scalar
    pool3%microbial         = this%microbial        * scalar
    pool3%residue           = this%residue          * scalar
    pool3%residue_assoc     = this%residue_assoc    * scalar
  END FUNCTION Multiply_pool_a_and_scalar
  !-----------------------------------------------------------------------------------------------------
  !> * (scalar * t_sb_pool_a)
  !!
  !-----------------------------------------------------------------------------------------------------
  PURE FUNCTION Multiply_scalar_and_pool_a(scalar,this2) RESULT(pool3)
    IMPLICIT NONE
    ! 0.1 InOut
    REAL(wp),           INTENT(in) :: scalar
    CLASS(t_sb_pool_a), INTENT(in) :: this2
    TYPE(t_sb_pool_a)              :: pool3
    ! 1.0 go
    pool3%dom               = scalar * this2%dom
    pool3%dom_assoc         = scalar * this2%dom_assoc
    pool3%soluable_litter   = scalar * this2%soluable_litter
    pool3%polymeric_litter  = scalar * this2%polymeric_litter
    pool3%woody_litter      = scalar * this2%woody_litter
    pool3%fungi             = scalar * this2%fungi
    pool3%mycorrhiza        = scalar * this2%mycorrhiza
    pool3%microbial         = scalar * this2%microbial
    pool3%residue           = scalar * this2%residue
    pool3%residue_assoc     = scalar * this2%residue_assoc
  END FUNCTION Multiply_scalar_and_pool_a
  !-----------------------------------------------------------------------------------------------------
  !> * (t_sb_pool_a * field2d)
  !!
  !-----------------------------------------------------------------------------------------------------
  PURE FUNCTION Multiply_pool_a_and_field2d(this,field2d) RESULT(pool3)
    IMPLICIT NONE
    ! 0.1 InOut
    CLASS(t_sb_pool_a), INTENT(in) :: this
    REAL(wp),           INTENT(in) :: field2d(:,:)
    TYPE(t_sb_pool_a)              :: pool3
    ! 1.0 go
    pool3%dom               = this%dom              * field2d(:,:)
    pool3%dom_assoc         = this%dom_assoc        * field2d(:,:)
    pool3%soluable_litter   = this%soluable_litter  * field2d(:,:)
    pool3%polymeric_litter  = this%polymeric_litter * field2d(:,:)
    pool3%woody_litter      = this%woody_litter     * field2d(:,:)
    pool3%fungi             = this%fungi            * field2d(:,:)
    pool3%mycorrhiza        = this%mycorrhiza       * field2d(:,:)
    pool3%microbial         = this%microbial        * field2d(:,:)
    pool3%residue           = this%residue          * field2d(:,:)
    pool3%residue_assoc     = this%residue_assoc    * field2d(:,:)
  END FUNCTION Multiply_pool_a_and_field2d
  !-----------------------------------------------------------------------------------------------------
  !> * (field2d * t_sb_pool_a)
  !!
  !-----------------------------------------------------------------------------------------------------
  PURE FUNCTION Multiply_field2d_and_pool_a(field2d,this2) RESULT(pool3)
    IMPLICIT NONE
    ! 0.1 InOut
    REAL(wp),           INTENT(in) :: field2d(:,:)
    CLASS(t_sb_pool_a), INTENT(in) :: this2
    TYPE(t_sb_pool_a)              :: pool3
    ! 1.0 go
    pool3%dom               = field2d(:,:) * this2%dom
    pool3%dom_assoc         = field2d(:,:) * this2%dom_assoc
    pool3%soluable_litter   = field2d(:,:) * this2%soluable_litter
    pool3%polymeric_litter  = field2d(:,:) * this2%polymeric_litter
    pool3%woody_litter      = field2d(:,:) * this2%woody_litter
    pool3%fungi             = field2d(:,:) * this2%fungi
    pool3%mycorrhiza        = field2d(:,:) * this2%mycorrhiza
    pool3%microbial         = field2d(:,:) * this2%microbial
    pool3%residue           = field2d(:,:) * this2%residue
    pool3%residue_assoc     = field2d(:,:) * this2%residue_assoc
  END FUNCTION Multiply_field2d_and_pool_a




  !-----------------------------------------------------------------------------------------------------
  !> Initialisation t_slm_struct_sb_element with a scalar
  !!
  !-----------------------------------------------------------------------------------------------------
  PURE ELEMENTAL SUBROUTINE init_t_slm_struct_sb_element(this,s)
    IMPLICIT NONE
    ! 0.1 InOut
    CLASS(t_slm_struct_sb_element),INTENT(inout) :: this 
    REAL(wp),                      INTENT(in)    :: s
    ! 1.0 Go
    this%carbon%ptr(:,:,:)     = s
    this%nitrogen%ptr(:,:,:)   = s
    this%phosphorus%ptr(:,:,:) = s
    this%carbon13%ptr(:,:,:)   = s
    this%carbon14%ptr(:,:,:)   = s
    this%nitrogen15%ptr(:,:,:) = s
  END SUBROUTINE init_t_slm_struct_sb_element 

  !-----------------------------------------------------------------------------------------------------
  !> Initialisation t_slm_struct_sb_pool with a scalar
  !!
  !-----------------------------------------------------------------------------------------------------
  PURE ELEMENTAL SUBROUTINE init_t_slm_struct_sb_pool(this,s)
    IMPLICIT NONE
    ! 0.1 InOut
    CLASS(t_slm_struct_sb_pool),INTENT(inout) :: this 
    REAL(wp),                   INTENT(in)    :: s
    ! 1.0 Go
    CALL this%dom%init(s)
    CALL this%dom_assoc%init(s)
    CALL this%soluable_litter%init(s)
    CALL this%polymeric_litter%init(s)
    CALL this%woody_litter%init(s)
    CALL this%fungi%init(s)
    CALL this%mycorrhiza%init(s)
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
    USE mo_jsb_io,            ONLY: grib_bits
!    USE mo_jsb_io,            ONLY: TSTEP_CONSTANT, tables
    USE mo_jsb_io,            ONLY: t_cf, t_grib1, t_grib2
    USE mo_jsb_grid_class,    ONLY: t_jsb_grid, t_jsb_vgrid
!    USE mo_jsb_grid,          ONLY: Get_grid, Get_vgrid
    USE mo_jsb_model_class,   ONLY: t_jsb_model

    ! 0.1.InOut
    CLASS(t_sb_memory), INTENT(inout), TARGET :: mem             !> soil memory
    CHARACTER(len=*),     INTENT(in)            :: prefix          !> ?
    CHARACTER(len=*),     INTENT(in)            :: suffix          !> ?
    INTEGER,              INTENT(in)            :: lct_ids !(:)    !> in jsbach4: vector of land cover type IDs
    INTEGER,              INTENT(in)            :: model_id        !> model ID model\%id
    ! 0.2 Local
    TYPE(t_jsb_model), POINTER :: model
    TYPE(t_jsb_grid),  POINTER :: hgrid                        ! Horizontal grid
!    TYPE(t_jsb_vgrid), POINTER :: canopy_layer                 ! Vertical grid
    TYPE(t_jsb_vgrid), POINTER :: surface                      ! Vertical grid
    INTEGER :: table

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':Init_sb_memory'

    model => Get_model(model_id)

!    table = tables(1)
    table = 1

!    hgrid   => Get_grid(mem%grid_id)
!    surface => Get_vgrid('surface')
!    canopy_layer => Get_vgrid('canopy_layer')


    IF (.NOT. ALLOCATED(mem%nh4_solute%ptr)) THEN

      !> add variables to soil memory soil\%VAR
      !! 
      !! soil\%VAR
      !!  
      CALL jsb_add_var(mem, 'nh4_solute', mem%nh4_solute, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'nh4_assoc', mem%nh4_assoc, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'no3_solute', mem%no3_solute, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'noy', mem%noy, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'n2o', mem%n2o, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'n2', mem%n2, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'nh4_n15_solute', mem%nh4_n15_solute, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'nh4_n15_assoc', mem%nh4_n15_assoc, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'no3_n15_solute', mem%no3_n15_solute, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'noy_n15', mem%noy_n15, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'n2o_n15', mem%n2o_n15, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'n2_n15', mem%n2_n15, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'po4_assoc_fast', mem%po4_assoc_fast, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'po4_assoc_slow', mem%po4_assoc_slow, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'po4_occluded', mem%po4_occluded, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue

      CALL jsb_add_var(mem, 'po4_primary', mem%po4_primary, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'po4_solute', mem%po4_solute, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'transport_nh4_solute', mem%transport_nh4_solute, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'transport_nh4_assoc', mem%transport_nh4_assoc, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'transport_no3_solute', mem%transport_no3_solute, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'transport_noy', mem%transport_noy, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'transport_n2o', mem%transport_n2o, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'transport_n2', mem%transport_n2, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'transport_nh4_n15_solute', mem%transport_nh4_n15_solute, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'transport_nh4_n15_assoc', mem%transport_nh4_n15_assoc, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'transport_no3_n15_solute', mem%transport_no3_n15_solute, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'transport_noy_n15', mem%transport_noy_n15, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'transport_n2o_n15', mem%transport_n2o_n15, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'transport_n2_n15', mem%transport_n2_n15, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'transport_po4_assoc_fast', mem%transport_po4_assoc_fast, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'transport_po4_assoc_slow', mem%transport_po4_assoc_slow, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'transport_po4_occluded', mem%transport_po4_occluded, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'transport_po4_primary', mem%transport_po4_primary, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue  

      CALL jsb_add_var(mem, 'transport_po4_solute', mem%transport_po4_solute, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'rtm_decomposition', mem%rtm_decomposition, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'rtm_depolymerisation', mem%rtm_depolymerisation, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'rtm_mic_uptake', mem%rtm_mic_uptake, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'rtm_sorption', mem%rtm_sorption, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'rtm_desorption', mem%rtm_desorption, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'rtm_hsc', mem%rtm_hsc, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'rtm_nitrification', mem%rtm_nitrification, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'rtm_denitrification', mem%rtm_denitrification, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'rtm_gasdiffusion', mem%rtm_gasdiffusion, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'rmm_decomposition', mem%rmm_decomposition, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'rmm_depolymerisation', mem%rmm_depolymerisation, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'rmm_mic_uptake', mem%rmm_mic_uptake, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'rmm_sorption', mem%rmm_sorption, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'rmm_desorption', mem%rmm_desorption, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'rmm_hsc', mem%rmm_hsc, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'rmm_nitrification', mem%rmm_nitrification, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'rmm_gasdiffusion', mem%rmm_gasdiffusion, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'het_respiration', mem%het_respiration, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'het_respiration_c13', mem%het_respiration_c13, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'het_respiration_c14', mem%het_respiration_c14, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'net_mineralisation_nh4', mem%net_mineralisation_nh4, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'net_mineralisation_nh4_n15', mem%net_mineralisation_nh4_n15, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'net_mineralisation_no3', mem%net_mineralisation_no3, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'net_mineralisation_no3_n15', mem%net_mineralisation_no3_n15, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'net_mineralisation_po4', mem%net_mineralisation_po4, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'plant_uptake_nh4_sl', mem%plant_uptake_nh4_sl, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'plant_uptake_nh4_n15_sl', mem%plant_uptake_nh4_n15_sl, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'plant_uptake_no3_sl', mem%plant_uptake_no3_sl, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'plant_uptake_no3_n15_sl', mem%plant_uptake_no3_n15_sl, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue

      CALL jsb_add_var(mem, 'plant_uptake_po4_sl', mem%plant_uptake_po4_sl, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue

      CALL jsb_add_var(mem, 'microbial_uptake_nh4_sl', mem%microbial_uptake_nh4_sl, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'microbial_uptake_nh4_n15_sl', mem%microbial_uptake_nh4_n15_sl, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'microbial_uptake_no3_sl', mem%microbial_uptake_no3_sl, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'microbial_uptake_no3_n15_sl', mem%microbial_uptake_no3_n15_sl, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue

      CALL jsb_add_var(mem, 'microbial_uptake_po4_sl', mem%microbial_uptake_po4_sl, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue

      CALL jsb_add_var(mem, 'mycorrhiza_uptake_nh4_sl', mem%mycorrhiza_uptake_nh4_sl, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue

      CALL jsb_add_var(mem, 'mycorrhiza_uptake_nh4_n15_sl', mem%mycorrhiza_uptake_nh4_n15_sl, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'mycorrhiza_uptake_no3_sl', mem%mycorrhiza_uptake_no3_sl, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'mycorrhiza_uptake_no3_n15_sl', mem%mycorrhiza_uptake_no3_n15_sl, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue

      CALL jsb_add_var(mem, 'mycorrhiza_uptake_norg_sl', mem%mycorrhiza_uptake_norg_sl, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue

      CALL jsb_add_var(mem, 'mycorrhiza_uptake_norg_n15_sl', mem%mycorrhiza_uptake_norg_n15_sl, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue

      CALL jsb_add_var(mem, 'mycorrhiza_uptake_po4_sl', mem%mycorrhiza_uptake_po4_sl, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'volatilisation_nh4', mem%volatilisation_nh4, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'volatilisation_nh4_n15', mem%volatilisation_nh4_n15, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'nitrification_no3', mem%nitrification_no3, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue

      CALL jsb_add_var(mem, 'nitrification_no3_n15', mem%nitrification_no3_n15, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue

      CALL jsb_add_var(mem, 'nitrification_noy', mem%nitrification_noy, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue

      CALL jsb_add_var(mem, 'nitrification_noy_n15', mem%nitrification_noy_n15, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue

      CALL jsb_add_var(mem, 'nitrification_n2o', mem%nitrification_n2o, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue

      CALL jsb_add_var(mem, 'nitrification_n2o_n15', mem%nitrification_n2o_n15, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue

      CALL jsb_add_var(mem, 'denitrification_noy', mem%denitrification_noy, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue

      CALL jsb_add_var(mem, 'denitrification_noy_n15', mem%denitrification_noy_n15, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue

      CALL jsb_add_var(mem, 'denitrification_n2o', mem%denitrification_n2o, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue

      CALL jsb_add_var(mem, 'denitrification_n2o_n15', mem%denitrification_n2o_n15, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue

      CALL jsb_add_var(mem, 'denitrification_n2', mem%denitrification_n2, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue

      CALL jsb_add_var(mem, 'denitrification_n2_n15', mem%denitrification_n2_n15, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue

      CALL jsb_add_var(mem, 'asymb_n_fixation', mem%asymb_n_fixation, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue

      CALL jsb_add_var(mem, 'asymb_n_fixation_n15', mem%asymb_n_fixation_n15, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue

      CALL jsb_add_var(mem, 'weathering_po4', mem%weathering_po4, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue

      CALL jsb_add_var(mem, 'occlusion_po4', mem%occlusion_po4, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue

      CALL jsb_add_var(mem, 'slow_sorption_po4', mem%slow_sorption_po4, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue

      CALL jsb_add_var(mem, 'fast_adsorption_po4', mem%fast_adsorption_po4, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue

      CALL jsb_add_var(mem, 'biochem_mineralisation_po4', mem%biochem_mineralisation_po4, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue

      CALL jsb_add_var(mem, 'gross_mineralisation_po4', mem%gross_mineralisation_po4, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue

      CALL jsb_add_var(mem, 'qmax_po4', mem%qmax_po4, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue

      CALL jsb_add_var(mem, 'qmax_org', mem%qmax_org, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue

      CALL jsb_add_var(mem, 'k_bioturb', mem%k_bioturb, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue

      CALL jsb_add_var(mem, 'bulk_dens_corr_sl', mem%bulk_dens_corr_sl, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue

      CALL jsb_add_var(mem, 'bulk_soil_carbon_sl', mem%bulk_soil_carbon_sl, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue

      CALL jsb_add_var(mem, 'particle_fluxrate', mem%particle_fluxrate, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue

      CALL jsb_add_var(mem, 'anaerobic_volume_fraction_sl', mem%anaerobic_volume_fraction_sl, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue

      CALL jsb_add_var(mem, 'ph_sl', mem%ph_sl, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue

      CALL jsb_add_var(mem, 'lithology_class_sl', mem%lithology_class_sl, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue

      CALL jsb_add_var(mem, 'leaching_nh4_solute', mem%leaching_nh4_solute, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue

      CALL jsb_add_var(mem, 'leaching_nh4_n15_solute', mem%leaching_nh4_n15_solute, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue

      CALL jsb_add_var(mem, 'leaching_no3_solute', mem%leaching_no3_solute, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue

      CALL jsb_add_var(mem, 'leaching_no3_n15_solute', mem%leaching_no3_n15_solute, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue

      CALL jsb_add_var(mem, 'leaching_po4_solute', mem%leaching_po4_solute, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue

      CALL jsb_add_var(mem, 'emission_nh3', mem%emission_nh3, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue

      CALL jsb_add_var(mem, 'emission_nh3_n15', mem%emission_nh3_n15, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue

      CALL jsb_add_var(mem, 'emission_noy', mem%emission_noy, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue

      CALL jsb_add_var(mem, 'emission_noy_n15', mem%emission_noy_n15, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue

      CALL jsb_add_var(mem, 'emission_n2o', mem%emission_n2o, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue

      CALL jsb_add_var(mem, 'emission_n2o_n15', mem%emission_n2o_n15, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue

      CALL jsb_add_var(mem, 'emission_n2', mem%emission_n2, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue

      CALL jsb_add_var(mem, 'emission_n2_n15', mem%emission_n2_n15, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue

      CALL jsb_add_var(mem, 'total_flux_carbon', mem%total_flux_carbon, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue

      CALL jsb_add_var(mem, 'total_flux_nitrogen', mem%total_flux_nitrogen, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue

      CALL jsb_add_var(mem, 'total_flux_phosphorus', mem%total_flux_phosphorus, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue

      CALL jsb_add_var(mem, 'total_flux_carbon13', mem%total_flux_carbon13, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue

      CALL jsb_add_var(mem, 'total_flux_carbon14', mem%total_flux_carbon14, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue

      CALL jsb_add_var(mem, 'total_flux_nitrogen15', mem%total_flux_nitrogen15, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue


      !> simple TYPE (t_jsb_struct)
      !! 
      !! jsb_add_var -> slm_add_struct()
      !! 
      CALL jsb_add_var(mem, 'sb_delta_pools', mem%delta_pools, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix)

      CALL jsb_add_var(mem, 'sb_mycorrhiza_export', mem%mycorrhiza_export, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix)


      !> nested TYPE (t_jsb_n_struct)
      !! 
      !! jsb_add_var -> slm_add_n_struct()
      !! 
      CALL jsb_add_var(mem, 'sb_pool'           , mem%pool, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix)
      CALL jsb_add_var(mem, 'sb_formation'      , mem%formation, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix)
      CALL jsb_add_var(mem, 'sb_loss'  , mem%loss, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix)
      CALL jsb_add_var(mem, 'sb_transport'      , mem%transport, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix)

      !> level 2 TYPE (l2_struct) nested within the above nested TYPE structures (t_jsb_n_struct)
      !! 
      !! jsb_add_var -> slm_add_l2_struct()
      !! 
      CALL jsb_add_var(mem%pool,    'sb_pool_dom',                mem%pool%dom, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix) 
      CALL jsb_add_var(mem%pool,    'sb_pool_dom_assoc',          mem%pool%dom_assoc, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix)
      CALL jsb_add_var(mem%pool,    'sb_pool_soluable_litter',    mem%pool%soluable_litter, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix) 
      CALL jsb_add_var(mem%pool,    'sb_pool_polymeric_litter',   mem%pool%polymeric_litter, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix) 
      CALL jsb_add_var(mem%pool,    'sb_pool_woody_litter',       mem%pool%woody_litter, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix) 
      CALL jsb_add_var(mem%pool,    'sb_pool_fungi',              mem%pool%fungi, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix) 
      CALL jsb_add_var(mem%pool,    'sb_pool_mycorrhiza',         mem%pool%mycorrhiza, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix) 
      CALL jsb_add_var(mem%pool,    'sb_pool_microbial',          mem%pool%microbial, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix) 
      CALL jsb_add_var(mem%pool,    'sb_pool_residue',            mem%pool%residue, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix) 
      CALL jsb_add_var(mem%pool,    'sb_pool_residue_assoc',      mem%pool%residue_assoc, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix)  

      CALL jsb_add_var(mem%formation,    'sb_formation_dom',                mem%formation%dom, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix) 
      CALL jsb_add_var(mem%formation,    'sb_formation_dom_assoc',          mem%formation%dom_assoc, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix)
      CALL jsb_add_var(mem%formation,    'sb_formation_soluable_litter',    mem%formation%soluable_litter, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix) 
      CALL jsb_add_var(mem%formation,    'sb_formation_polymeric_litter',   mem%formation%polymeric_litter, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix) 
      CALL jsb_add_var(mem%formation,    'sb_formation_woody_litter',       mem%formation%woody_litter, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix) 
      CALL jsb_add_var(mem%formation,    'sb_formation_fungi',              mem%formation%fungi, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix) 
      CALL jsb_add_var(mem%formation,    'sb_formation_mycorrhiza',         mem%formation%mycorrhiza, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix) 
      CALL jsb_add_var(mem%formation,    'sb_formation_microbial',          mem%formation%microbial, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix) 
      CALL jsb_add_var(mem%formation,    'sb_formation_residue',            mem%formation%residue, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix) 
      CALL jsb_add_var(mem%formation,    'sb_formation_residue_assoc',      mem%formation%residue_assoc, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix)  

      CALL jsb_add_var(mem%loss,    'sb_loss_dom',                mem%loss%dom, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix) 
      CALL jsb_add_var(mem%loss,    'sb_loss_dom_assoc',          mem%loss%dom_assoc, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix)
      CALL jsb_add_var(mem%loss,    'sb_loss_soluable_litter',    mem%loss%soluable_litter, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix) 
      CALL jsb_add_var(mem%loss,    'sb_loss_polymeric_litter',   mem%loss%polymeric_litter, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix) 
      CALL jsb_add_var(mem%loss,    'sb_loss_woody_litter',       mem%loss%woody_litter, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix) 
      CALL jsb_add_var(mem%loss,    'sb_loss_fungi',              mem%loss%fungi, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix) 
      CALL jsb_add_var(mem%loss,    'sb_loss_mycorrhiza',         mem%loss%mycorrhiza, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix) 
      CALL jsb_add_var(mem%loss,    'sb_loss_microbial',          mem%loss%microbial, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix) 
      CALL jsb_add_var(mem%loss,    'sb_loss_residue',            mem%loss%residue, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix) 
      CALL jsb_add_var(mem%loss,    'sb_loss_residue_assoc',      mem%loss%residue_assoc, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix)  

      CALL jsb_add_var(mem%transport,    'sb_transport_dom',                mem%transport%dom, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix) 
      CALL jsb_add_var(mem%transport,    'sb_transport_dom_assoc',          mem%transport%dom_assoc, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix)
      CALL jsb_add_var(mem%transport,    'sb_transport_soluable_litter',    mem%transport%soluable_litter, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix) 
      CALL jsb_add_var(mem%transport,    'sb_transport_polymeric_litter',   mem%transport%polymeric_litter, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix) 
      CALL jsb_add_var(mem%transport,    'sb_transport_woody_litter',       mem%transport%woody_litter, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix) 
      CALL jsb_add_var(mem%transport,    'sb_transport_fungi',              mem%transport%fungi, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix) 
      CALL jsb_add_var(mem%transport,    'sb_transport_mycorrhiza',         mem%transport%mycorrhiza, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix) 
      CALL jsb_add_var(mem%transport,    'sb_transport_microbial',          mem%transport%microbial, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix) 
      CALL jsb_add_var(mem%transport,    'sb_transport_residue',            mem%transport%residue, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix) 
      CALL jsb_add_var(mem%transport,    'sb_transport_residue_assoc',      mem%transport%residue_assoc, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix)  

      !> add variables to the structs and the l2 structs
      !!
      !! soil\%t_jsb_struct\%VAR
      !! 
      !! soil\%t_jsb_n_struct\%t_jsb_struct\%VAR
      !! 
      CALL jsb_add_var(mem%pool%dom, 'sb_pool_dom_carbon', mem%pool%dom%carbon, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%pool%dom, 'sb_pool_dom_nitrogen', mem%pool%dom%nitrogen, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%pool%dom, 'sb_pool_dom_phosphorus', mem%pool%dom%phosphorus, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%pool%dom, 'sb_pool_dom_carbon13', mem%pool%dom%carbon13, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%pool%dom, 'sb_pool_dom_carbon14', mem%pool%dom%carbon14, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%pool%dom, 'sb_pool_dom_nitrogen15', mem%pool%dom%nitrogen15, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
 
      CALL jsb_add_var(mem%pool%dom_assoc, 'sb_pool_dom_assoc_carbon', mem%pool%dom_assoc%carbon, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%pool%dom_assoc, 'sb_pool_dom_assoc_nitrogen', mem%pool%dom_assoc%nitrogen, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%pool%dom_assoc, 'sb_pool_dom_assoc_phosphorus', mem%pool%dom_assoc%phosphorus, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%pool%dom_assoc, 'sb_pool_dom_assoc_carbon13', mem%pool%dom_assoc%carbon13, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%pool%dom_assoc, 'sb_pool_dom_assoc_carbon14', mem%pool%dom_assoc%carbon14, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%pool%dom_assoc, 'sb_pool_dom_assoc_nitrogen15', mem%pool%dom_assoc%nitrogen15, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
 
      CALL jsb_add_var(mem%pool%soluable_litter, 'sb_pool_soluable_litter_carbon', mem%pool%soluable_litter%carbon, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%pool%soluable_litter, 'sb_pool_soluable_litter_nitrogen', mem%pool%soluable_litter%nitrogen, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%pool%soluable_litter, 'sb_pool_soluable_litter_phosphorus', mem%pool%soluable_litter%phosphorus, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%pool%soluable_litter, 'sb_pool_soluable_litter_carbon13', mem%pool%soluable_litter%carbon13, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%pool%soluable_litter, 'sb_pool_soluable_litter_carbon14', mem%pool%soluable_litter%carbon14, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%pool%soluable_litter, 'sb_pool_soluable_litter_nitrogen15', mem%pool%soluable_litter%nitrogen15, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
 
      CALL jsb_add_var(mem%pool%polymeric_litter, 'sb_pool_polymeric_litter_carbon', mem%pool%polymeric_litter%carbon, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%pool%polymeric_litter, 'sb_pool_polymeric_litter_nitrogen', mem%pool%polymeric_litter%nitrogen, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%pool%polymeric_litter, 'sb_pool_polymeric_litter_phosphorus', mem%pool%polymeric_litter%phosphorus, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%pool%polymeric_litter, 'sb_pool_polymeric_litter_carbon13', mem%pool%polymeric_litter%carbon13, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%pool%polymeric_litter, 'sb_pool_polymeric_litter_carbon14', mem%pool%polymeric_litter%carbon14, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%pool%polymeric_litter, 'sb_pool_polymeric_litter_nitrogen15', mem%pool%polymeric_litter%nitrogen15, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
 
      CALL jsb_add_var(mem%pool%woody_litter, 'sb_pool_woody_litter_carbon', mem%pool%woody_litter%carbon, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%pool%woody_litter, 'sb_pool_woody_litter_nitrogen', mem%pool%woody_litter%nitrogen, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%pool%woody_litter, 'sb_pool_woody_litter_phosphorus', mem%pool%woody_litter%phosphorus, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%pool%woody_litter, 'sb_pool_woody_litter_carbon13', mem%pool%woody_litter%carbon13, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%pool%woody_litter, 'sb_pool_woody_litter_carbon14', mem%pool%woody_litter%carbon14, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%pool%woody_litter, 'sb_pool_woody_litter_nitrogen15', mem%pool%woody_litter%nitrogen15, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
 
      CALL jsb_add_var(mem%pool%fungi, 'sb_pool_fungi_carbon', mem%pool%fungi%carbon, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%pool%fungi, 'sb_pool_fungi_nitrogen', mem%pool%fungi%nitrogen, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%pool%fungi, 'sb_pool_fungi_phosphorus', mem%pool%fungi%phosphorus, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%pool%fungi, 'sb_pool_fungi_carbon13', mem%pool%fungi%carbon13, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%pool%fungi, 'sb_pool_fungi_carbon14', mem%pool%fungi%carbon14, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%pool%fungi, 'sb_pool_fungi_nitrogen15', mem%pool%fungi%nitrogen15, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
 
      CALL jsb_add_var(mem%pool%mycorrhiza, 'sb_pool_mycorrhiza_carbon', mem%pool%mycorrhiza%carbon, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%pool%mycorrhiza, 'sb_pool_mycorrhiza_nitrogen', mem%pool%mycorrhiza%nitrogen, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%pool%mycorrhiza, 'sb_pool_mycorrhiza_phosphorus', mem%pool%mycorrhiza%phosphorus, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%pool%mycorrhiza, 'sb_pool_mycorrhiza_carbon13', mem%pool%mycorrhiza%carbon13, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%pool%mycorrhiza, 'sb_pool_mycorrhiza_carbon14', mem%pool%mycorrhiza%carbon14, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%pool%mycorrhiza, 'sb_pool_mycorrhiza_nitrogen15', mem%pool%mycorrhiza%nitrogen15, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
 
      CALL jsb_add_var(mem%pool%microbial, 'sb_pool_microbial_carbon', mem%pool%microbial%carbon, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%pool%microbial, 'sb_pool_microbial_nitrogen', mem%pool%microbial%nitrogen, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%pool%microbial, 'sb_pool_microbial_phosphorus', mem%pool%microbial%phosphorus, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%pool%microbial, 'sb_pool_microbial_carbon13', mem%pool%microbial%carbon13, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%pool%microbial, 'sb_pool_microbial_carbon14', mem%pool%microbial%carbon14, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%pool%microbial, 'sb_pool_microbial_nitrogen15', mem%pool%microbial%nitrogen15, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
 
      CALL jsb_add_var(mem%pool%residue, 'sb_pool_residue_carbon', mem%pool%residue%carbon, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%pool%residue, 'sb_pool_residue_nitrogen', mem%pool%residue%nitrogen, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%pool%residue, 'sb_pool_residue_phosphorus', mem%pool%residue%phosphorus, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%pool%residue, 'sb_pool_residue_carbon13', mem%pool%residue%carbon13, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%pool%residue, 'sb_pool_residue_carbon14', mem%pool%residue%carbon14, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%pool%residue, 'sb_pool_residue_nitrogen15', mem%pool%residue%nitrogen15, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
 
      CALL jsb_add_var(mem%pool%residue_assoc, 'sb_pool_residue_assoc_carbon', mem%pool%residue_assoc%carbon, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%pool%residue_assoc, 'sb_pool_residue_assoc_nitrogen', mem%pool%residue_assoc%nitrogen, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%pool%residue_assoc, 'sb_pool_residue_assoc_phosphorus', mem%pool%residue_assoc%phosphorus, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%pool%residue_assoc, 'sb_pool_residue_assoc_carbon13', mem%pool%residue_assoc%carbon13, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%pool%residue_assoc, 'sb_pool_residue_assoc_carbon14', mem%pool%residue_assoc%carbon14, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%pool%residue_assoc, 'sb_pool_residue_assoc_nitrogen15', mem%pool%residue_assoc%nitrogen15, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue

! formation
      CALL jsb_add_var(mem%formation%dom, 'sb_formation_dom_carbon', mem%formation%dom%carbon, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%formation%dom, 'sb_formation_dom_nitrogen', mem%formation%dom%nitrogen, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%formation%dom, 'sb_formation_dom_phosphorus', mem%formation%dom%phosphorus, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%formation%dom, 'sb_formation_dom_carbon13', mem%formation%dom%carbon13, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%formation%dom, 'sb_formation_dom_carbon14', mem%formation%dom%carbon14, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%formation%dom, 'sb_formation_dom_nitrogen15', mem%formation%dom%nitrogen15, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
 
      CALL jsb_add_var(mem%formation%dom_assoc, 'sb_formation_dom_assoc_carbon', mem%formation%dom_assoc%carbon, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%formation%dom_assoc, 'sb_formation_dom_assoc_nitrogen', mem%formation%dom_assoc%nitrogen, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%formation%dom_assoc, 'sb_formation_dom_assoc_phosphorus', mem%formation%dom_assoc%phosphorus, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%formation%dom_assoc, 'sb_formation_dom_assoc_carbon13', mem%formation%dom_assoc%carbon13, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%formation%dom_assoc, 'sb_formation_dom_assoc_carbon14', mem%formation%dom_assoc%carbon14, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%formation%dom_assoc, 'sb_formation_dom_assoc_nitrogen15', mem%formation%dom_assoc%nitrogen15, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
 
      CALL jsb_add_var(mem%formation%soluable_litter, 'sb_formation_soluable_litter_carbon', mem%formation%soluable_litter%carbon, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%formation%soluable_litter, 'sb_formation_soluable_litter_nitrogen', mem%formation%soluable_litter%nitrogen, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%formation%soluable_litter, 'sb_formation_soluable_litter_phosphorus', mem%formation%soluable_litter%phosphorus, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%formation%soluable_litter, 'sb_formation_soluable_litter_carbon13', mem%formation%soluable_litter%carbon13, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%formation%soluable_litter, 'sb_formation_soluable_litter_carbon14', mem%formation%soluable_litter%carbon14, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%formation%soluable_litter, 'sb_formation_soluable_litter_nitrogen15', mem%formation%soluable_litter%nitrogen15, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
 
      CALL jsb_add_var(mem%formation%polymeric_litter, 'sb_formation_polymeric_litter_carbon', mem%formation%polymeric_litter%carbon, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%formation%polymeric_litter, 'sb_formation_polymeric_litter_nitrogen', mem%formation%polymeric_litter%nitrogen, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%formation%polymeric_litter, 'sb_formation_polymeric_litter_phosphorus', mem%formation%polymeric_litter%phosphorus, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%formation%polymeric_litter, 'sb_formation_polymeric_litter_carbon13', mem%formation%polymeric_litter%carbon13, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%formation%polymeric_litter, 'sb_formation_polymeric_litter_carbon14', mem%formation%polymeric_litter%carbon14, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%formation%polymeric_litter, 'sb_formation_polymeric_litter_nitrogen15', mem%formation%polymeric_litter%nitrogen15, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
 
      CALL jsb_add_var(mem%formation%woody_litter, 'sb_formation_woody_litter_carbon', mem%formation%woody_litter%carbon, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%formation%woody_litter, 'sb_formation_woody_litter_nitrogen', mem%formation%woody_litter%nitrogen, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%formation%woody_litter, 'sb_formation_woody_litter_phosphorus', mem%formation%woody_litter%phosphorus, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%formation%woody_litter, 'sb_formation_woody_litter_carbon13', mem%formation%woody_litter%carbon13, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%formation%woody_litter, 'sb_formation_woody_litter_carbon14', mem%formation%woody_litter%carbon14, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%formation%woody_litter, 'sb_formation_woody_litter_nitrogen15', mem%formation%woody_litter%nitrogen15, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
 
      CALL jsb_add_var(mem%formation%fungi, 'sb_formation_fungi_carbon', mem%formation%fungi%carbon, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%formation%fungi, 'sb_formation_fungi_nitrogen', mem%formation%fungi%nitrogen, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%formation%fungi, 'sb_formation_fungi_phosphorus', mem%formation%fungi%phosphorus, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%formation%fungi, 'sb_formation_fungi_carbon13', mem%formation%fungi%carbon13, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%formation%fungi, 'sb_formation_fungi_carbon14', mem%formation%fungi%carbon14, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%formation%fungi, 'sb_formation_fungi_nitrogen15', mem%formation%fungi%nitrogen15, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
 
      CALL jsb_add_var(mem%formation%mycorrhiza, 'sb_formation_mycorrhiza_carbon', mem%formation%mycorrhiza%carbon, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%formation%mycorrhiza, 'sb_formation_mycorrhiza_nitrogen', mem%formation%mycorrhiza%nitrogen, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%formation%mycorrhiza, 'sb_formation_mycorrhiza_phosphorus', mem%formation%mycorrhiza%phosphorus, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%formation%mycorrhiza, 'sb_formation_mycorrhiza_carbon13', mem%formation%mycorrhiza%carbon13, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%formation%mycorrhiza, 'sb_formation_mycorrhiza_carbon14', mem%formation%mycorrhiza%carbon14, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%formation%mycorrhiza, 'sb_formation_mycorrhiza_nitrogen15', mem%formation%mycorrhiza%nitrogen15, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
 
      CALL jsb_add_var(mem%formation%microbial, 'sb_formation_microbial_carbon', mem%formation%microbial%carbon, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%formation%microbial, 'sb_formation_microbial_nitrogen', mem%formation%microbial%nitrogen, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%formation%microbial, 'sb_formation_microbial_phosphorus', mem%formation%microbial%phosphorus, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%formation%microbial, 'sb_formation_microbial_carbon13', mem%formation%microbial%carbon13, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%formation%microbial, 'sb_formation_microbial_carbon14', mem%formation%microbial%carbon14, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%formation%microbial, 'sb_formation_microbial_nitrogen15', mem%formation%microbial%nitrogen15, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
 
      CALL jsb_add_var(mem%formation%residue, 'sb_formation_residue_carbon', mem%formation%residue%carbon, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%formation%residue, 'sb_formation_residue_nitrogen', mem%formation%residue%nitrogen, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%formation%residue, 'sb_formation_residue_phosphorus', mem%formation%residue%phosphorus, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%formation%residue, 'sb_formation_residue_carbon13', mem%formation%residue%carbon13, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%formation%residue, 'sb_formation_residue_carbon14', mem%formation%residue%carbon14, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%formation%residue, 'sb_formation_residue_nitrogen15', mem%formation%residue%nitrogen15, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
 
      CALL jsb_add_var(mem%formation%residue_assoc, 'sb_formation_residue_assoc_carbon', mem%formation%residue_assoc%carbon, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%formation%residue_assoc, 'sb_formation_residue_assoc_nitrogen', mem%formation%residue_assoc%nitrogen, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%formation%residue_assoc, 'sb_formation_residue_assoc_phosphorus', mem%formation%residue_assoc%phosphorus, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%formation%residue_assoc, 'sb_formation_residue_assoc_carbon13', mem%formation%residue_assoc%carbon13, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%formation%residue_assoc, 'sb_formation_residue_assoc_carbon14', mem%formation%residue_assoc%carbon14, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%formation%residue_assoc, 'sb_formation_residue_assoc_nitrogen15', mem%formation%residue_assoc%nitrogen15, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue

!loss
      CALL jsb_add_var(mem%loss%dom, 'sb_loss_dom_carbon', mem%loss%dom%carbon, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%loss%dom, 'sb_loss_dom_nitrogen', mem%loss%dom%nitrogen, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%loss%dom, 'sb_loss_dom_phosphorus', mem%loss%dom%phosphorus, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%loss%dom, 'sb_loss_dom_carbon13', mem%loss%dom%carbon13, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%loss%dom, 'sb_loss_dom_carbon14', mem%loss%dom%carbon14, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%loss%dom, 'sb_loss_dom_nitrogen15', mem%loss%dom%nitrogen15, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
 
      CALL jsb_add_var(mem%loss%dom_assoc, 'sb_loss_dom_assoc_carbon', mem%loss%dom_assoc%carbon, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%loss%dom_assoc, 'sb_loss_dom_assoc_nitrogen', mem%loss%dom_assoc%nitrogen, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%loss%dom_assoc, 'sb_loss_dom_assoc_phosphorus', mem%loss%dom_assoc%phosphorus, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%loss%dom_assoc, 'sb_loss_dom_assoc_carbon13', mem%loss%dom_assoc%carbon13, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%loss%dom_assoc, 'sb_loss_dom_assoc_carbon14', mem%loss%dom_assoc%carbon14, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%loss%dom_assoc, 'sb_loss_dom_assoc_nitrogen15', mem%loss%dom_assoc%nitrogen15, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
 
      CALL jsb_add_var(mem%loss%soluable_litter, 'sb_loss_soluable_litter_carbon', mem%loss%soluable_litter%carbon, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%loss%soluable_litter, 'sb_loss_soluable_litter_nitrogen', mem%loss%soluable_litter%nitrogen, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%loss%soluable_litter, 'sb_loss_soluable_litter_phosphorus', mem%loss%soluable_litter%phosphorus, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%loss%soluable_litter, 'sb_loss_soluable_litter_carbon13', mem%loss%soluable_litter%carbon13, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%loss%soluable_litter, 'sb_loss_soluable_litter_carbon14', mem%loss%soluable_litter%carbon14, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%loss%soluable_litter, 'sb_loss_soluable_litter_nitrogen15', mem%loss%soluable_litter%nitrogen15, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
 
      CALL jsb_add_var(mem%loss%polymeric_litter, 'sb_loss_polymeric_litter_carbon', mem%loss%polymeric_litter%carbon, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%loss%polymeric_litter, 'sb_loss_polymeric_litter_nitrogen', mem%loss%polymeric_litter%nitrogen, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%loss%polymeric_litter, 'sb_loss_polymeric_litter_phosphorus', mem%loss%polymeric_litter%phosphorus, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%loss%polymeric_litter, 'sb_loss_polymeric_litter_carbon13', mem%loss%polymeric_litter%carbon13, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%loss%polymeric_litter, 'sb_loss_polymeric_litter_carbon14', mem%loss%polymeric_litter%carbon14, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%loss%polymeric_litter, 'sb_loss_polymeric_litter_nitrogen15', mem%loss%polymeric_litter%nitrogen15, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
 
      CALL jsb_add_var(mem%loss%woody_litter, 'sb_loss_woody_litter_carbon', mem%loss%woody_litter%carbon, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%loss%woody_litter, 'sb_loss_woody_litter_nitrogen', mem%loss%woody_litter%nitrogen, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%loss%woody_litter, 'sb_loss_woody_litter_phosphorus', mem%loss%woody_litter%phosphorus, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%loss%woody_litter, 'sb_loss_woody_litter_carbon13', mem%loss%woody_litter%carbon13, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%loss%woody_litter, 'sb_loss_woody_litter_carbon14', mem%loss%woody_litter%carbon14, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%loss%woody_litter, 'sb_loss_woody_litter_nitrogen15', mem%loss%woody_litter%nitrogen15, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
 
      CALL jsb_add_var(mem%loss%fungi, 'sb_loss_fungi_carbon', mem%loss%fungi%carbon, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%loss%fungi, 'sb_loss_fungi_nitrogen', mem%loss%fungi%nitrogen, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%loss%fungi, 'sb_loss_fungi_phosphorus', mem%loss%fungi%phosphorus, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%loss%fungi, 'sb_loss_fungi_carbon13', mem%loss%fungi%carbon13, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%loss%fungi, 'sb_loss_fungi_carbon14', mem%loss%fungi%carbon14, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%loss%fungi, 'sb_loss_fungi_nitrogen15', mem%loss%fungi%nitrogen15, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
 
      CALL jsb_add_var(mem%loss%mycorrhiza, 'sb_loss_mycorrhiza_carbon', mem%loss%mycorrhiza%carbon, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%loss%mycorrhiza, 'sb_loss_mycorrhiza_nitrogen', mem%loss%mycorrhiza%nitrogen, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%loss%mycorrhiza, 'sb_loss_mycorrhiza_phosphorus', mem%loss%mycorrhiza%phosphorus, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%loss%mycorrhiza, 'sb_loss_mycorrhiza_carbon13', mem%loss%mycorrhiza%carbon13, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%loss%mycorrhiza, 'sb_loss_mycorrhiza_carbon14', mem%loss%mycorrhiza%carbon14, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%loss%mycorrhiza, 'sb_loss_mycorrhiza_nitrogen15', mem%loss%mycorrhiza%nitrogen15, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
 
      CALL jsb_add_var(mem%loss%microbial, 'sb_loss_microbial_carbon', mem%loss%microbial%carbon, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%loss%microbial, 'sb_loss_microbial_nitrogen', mem%loss%microbial%nitrogen, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%loss%microbial, 'sb_loss_microbial_phosphorus', mem%loss%microbial%phosphorus, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%loss%microbial, 'sb_loss_microbial_carbon13', mem%loss%microbial%carbon13, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%loss%microbial, 'sb_loss_microbial_carbon14', mem%loss%microbial%carbon14, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%loss%microbial, 'sb_loss_microbial_nitrogen15', mem%loss%microbial%nitrogen15, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
 
      CALL jsb_add_var(mem%loss%residue, 'sb_loss_residue_carbon', mem%loss%residue%carbon, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%loss%residue, 'sb_loss_residue_nitrogen', mem%loss%residue%nitrogen, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%loss%residue, 'sb_loss_residue_phosphorus', mem%loss%residue%phosphorus, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%loss%residue, 'sb_loss_residue_carbon13', mem%loss%residue%carbon13, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%loss%residue, 'sb_loss_residue_carbon14', mem%loss%residue%carbon14, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%loss%residue, 'sb_loss_residue_nitrogen15', mem%loss%residue%nitrogen15, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
 
      CALL jsb_add_var(mem%loss%residue_assoc, 'sb_loss_residue_assoc_carbon', mem%loss%residue_assoc%carbon, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%loss%residue_assoc, 'sb_loss_residue_assoc_nitrogen', mem%loss%residue_assoc%nitrogen, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%loss%residue_assoc, 'sb_loss_residue_assoc_phosphorus', mem%loss%residue_assoc%phosphorus, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%loss%residue_assoc, 'sb_loss_residue_assoc_carbon13', mem%loss%residue_assoc%carbon13, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%loss%residue_assoc, 'sb_loss_residue_assoc_carbon14', mem%loss%residue_assoc%carbon14, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%loss%residue_assoc, 'sb_loss_residue_assoc_nitrogen15', mem%loss%residue_assoc%nitrogen15, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue

! Transport
      CALL jsb_add_var(mem%transport%dom, 'sb_transport_dom_carbon', mem%transport%dom%carbon, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%transport%dom, 'sb_transport_dom_nitrogen', mem%transport%dom%nitrogen, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%transport%dom, 'sb_transport_dom_phosphorus', mem%transport%dom%phosphorus, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%transport%dom, 'sb_transport_dom_carbon13', mem%transport%dom%carbon13, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%transport%dom, 'sb_transport_dom_carbon14', mem%transport%dom%carbon14, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%transport%dom, 'sb_transport_dom_nitrogen15', mem%transport%dom%nitrogen15, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
 
      CALL jsb_add_var(mem%transport%dom_assoc, 'sb_transport_dom_assoc_carbon', mem%transport%dom_assoc%carbon, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%transport%dom_assoc, 'sb_transport_dom_assoc_nitrogen', mem%transport%dom_assoc%nitrogen, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%transport%dom_assoc, 'sb_transport_dom_assoc_phosphorus', mem%transport%dom_assoc%phosphorus, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%transport%dom_assoc, 'sb_transport_dom_assoc_carbon13', mem%transport%dom_assoc%carbon13, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%transport%dom_assoc, 'sb_transport_dom_assoc_carbon14', mem%transport%dom_assoc%carbon14, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%transport%dom_assoc, 'sb_transport_dom_assoc_nitrogen15', mem%transport%dom_assoc%nitrogen15, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
 
      CALL jsb_add_var(mem%transport%soluable_litter, 'sb_transport_soluable_litter_carbon', mem%transport%soluable_litter%carbon, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%transport%soluable_litter, 'sb_transport_soluable_litter_nitrogen', mem%transport%soluable_litter%nitrogen, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%transport%soluable_litter, 'sb_transport_soluable_litter_phosphorus', mem%transport%soluable_litter%phosphorus, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%transport%soluable_litter, 'sb_transport_soluable_litter_carbon13', mem%transport%soluable_litter%carbon13, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%transport%soluable_litter, 'sb_transport_soluable_litter_carbon14', mem%transport%soluable_litter%carbon14, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%transport%soluable_litter, 'sb_transport_soluable_litter_nitrogen15', mem%transport%soluable_litter%nitrogen15, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
 
      CALL jsb_add_var(mem%transport%polymeric_litter, 'sb_transport_polymeric_litter_carbon', mem%transport%polymeric_litter%carbon, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%transport%polymeric_litter, 'sb_transport_polymeric_litter_nitrogen', mem%transport%polymeric_litter%nitrogen, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%transport%polymeric_litter, 'sb_transport_polymeric_litter_phosphorus', mem%transport%polymeric_litter%phosphorus, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%transport%polymeric_litter, 'sb_transport_polymeric_litter_carbon13', mem%transport%polymeric_litter%carbon13, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%transport%polymeric_litter, 'sb_transport_polymeric_litter_carbon14', mem%transport%polymeric_litter%carbon14, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%transport%polymeric_litter, 'sb_transport_polymeric_litter_nitrogen15', mem%transport%polymeric_litter%nitrogen15, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
 
      CALL jsb_add_var(mem%transport%woody_litter, 'sb_transport_woody_litter_carbon', mem%transport%woody_litter%carbon, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%transport%woody_litter, 'sb_transport_woody_litter_nitrogen', mem%transport%woody_litter%nitrogen, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%transport%woody_litter, 'sb_transport_woody_litter_phosphorus', mem%transport%woody_litter%phosphorus, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%transport%woody_litter, 'sb_transport_woody_litter_carbon13', mem%transport%woody_litter%carbon13, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%transport%woody_litter, 'sb_transport_woody_litter_carbon14', mem%transport%woody_litter%carbon14, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%transport%woody_litter, 'sb_transport_woody_litter_nitrogen15', mem%transport%woody_litter%nitrogen15, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
 
      CALL jsb_add_var(mem%transport%fungi, 'sb_transport_fungi_carbon', mem%transport%fungi%carbon, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%transport%fungi, 'sb_transport_fungi_nitrogen', mem%transport%fungi%nitrogen, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%transport%fungi, 'sb_transport_fungi_phosphorus', mem%transport%fungi%phosphorus, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%transport%fungi, 'sb_transport_fungi_carbon13', mem%transport%fungi%carbon13, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%transport%fungi, 'sb_transport_fungi_carbon14', mem%transport%fungi%carbon14, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%transport%fungi, 'sb_transport_fungi_nitrogen15', mem%transport%fungi%nitrogen15, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
 
      CALL jsb_add_var(mem%transport%mycorrhiza, 'sb_transport_mycorrhiza_carbon', mem%transport%mycorrhiza%carbon, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%transport%mycorrhiza, 'sb_transport_mycorrhiza_nitrogen', mem%transport%mycorrhiza%nitrogen, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%transport%mycorrhiza, 'sb_transport_mycorrhiza_phosphorus', mem%transport%mycorrhiza%phosphorus, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%transport%mycorrhiza, 'sb_transport_mycorrhiza_carbon13', mem%transport%mycorrhiza%carbon13, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%transport%mycorrhiza, 'sb_transport_mycorrhiza_carbon14', mem%transport%mycorrhiza%carbon14, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%transport%mycorrhiza, 'sb_transport_mycorrhiza_nitrogen15', mem%transport%mycorrhiza%nitrogen15, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
 
      CALL jsb_add_var(mem%transport%microbial, 'sb_transport_microbial_carbon', mem%transport%microbial%carbon, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%transport%microbial, 'sb_transport_microbial_nitrogen', mem%transport%microbial%nitrogen, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%transport%microbial, 'sb_transport_microbial_phosphorus', mem%transport%microbial%phosphorus, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%transport%microbial, 'sb_transport_microbial_carbon13', mem%transport%microbial%carbon13, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%transport%microbial, 'sb_transport_microbial_carbon14', mem%transport%microbial%carbon14, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%transport%microbial, 'sb_transport_microbial_nitrogen15', mem%transport%microbial%nitrogen15, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
 
      CALL jsb_add_var(mem%transport%residue, 'sb_transport_residue_carbon', mem%transport%residue%carbon, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%transport%residue, 'sb_transport_residue_nitrogen', mem%transport%residue%nitrogen, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%transport%residue, 'sb_transport_residue_phosphorus', mem%transport%residue%phosphorus, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%transport%residue, 'sb_transport_residue_carbon13', mem%transport%residue%carbon13, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%transport%residue, 'sb_transport_residue_carbon14', mem%transport%residue%carbon14, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%transport%residue, 'sb_transport_residue_nitrogen15', mem%transport%residue%nitrogen15, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
 
      CALL jsb_add_var(mem%transport%residue_assoc, 'sb_transport_residue_assoc_carbon', mem%transport%residue_assoc%carbon, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%transport%residue_assoc, 'sb_transport_residue_assoc_nitrogen', mem%transport%residue_assoc%nitrogen, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%transport%residue_assoc, 'sb_transport_residue_assoc_phosphorus', mem%transport%residue_assoc%phosphorus, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%transport%residue_assoc, 'sb_transport_residue_assoc_carbon13', mem%transport%residue_assoc%carbon13, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%transport%residue_assoc, 'sb_transport_residue_assoc_carbon14', mem%transport%residue_assoc%carbon14, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%transport%residue_assoc, 'sb_transport_residue_assoc_nitrogen15', mem%transport%residue_assoc%nitrogen15, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue

      CALL jsb_add_var(mem%delta_pools, 'sb_delta_pools_carbon', mem%delta_pools%carbon, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%delta_pools, 'sb_delta_pools_nitrogen', mem%delta_pools%nitrogen, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%delta_pools, 'sb_delta_pools_phosphorus', mem%delta_pools%phosphorus, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%delta_pools, 'sb_delta_pools_carbon13', mem%delta_pools%carbon13, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%delta_pools, 'sb_delta_pools_carbon14', mem%delta_pools%carbon14, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%delta_pools, 'sb_delta_pools_nitrogen15', mem%delta_pools%nitrogen15, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue

      CALL jsb_add_var(mem%mycorrhiza_export, 'sb_mycorrhiza_export_carbon', mem%mycorrhiza_export%carbon, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%mycorrhiza_export, 'sb_mycorrhiza_export_nitrogen', mem%mycorrhiza_export%nitrogen, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%mycorrhiza_export, 'sb_mycorrhiza_export_phosphorus', mem%mycorrhiza_export%phosphorus, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%mycorrhiza_export, 'sb_mycorrhiza_export_carbon13', mem%mycorrhiza_export%carbon13, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%mycorrhiza_export, 'sb_mycorrhiza_export_carbon14', mem%mycorrhiza_export%carbon14, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue
      CALL jsb_add_var(mem%mycorrhiza_export, 'sb_mycorrhiza_export_nitrogen15', mem%mycorrhiza_export%nitrogen15, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp)     ! initvalue

    END IF

  END SUBROUTINE Init_sb_memory


  !-----------------------------------------------------------------------------------------------------
  !> ALLOCATE variables in sb_pool_a and sb_element_a structures
  !!
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE Init_sb_pools_elements_a(ics, ice, nc, iblk, nsoil_sb)

    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    INTEGER, INTENT(in)     :: ics, ice, nc, iblk, nsoil_sb   !< information about dimensions of variables in memory
    ! ---------------------------
    ! 0.2 Local
    INTEGER                 :: npoints                        !< number of point in the chunk
    
    
    ! temporary solution
    npoints = ice


    !< 1.0  ALLOCATE variables
    !! 
    ALLOCATE(sb_pool%dom%carbon(npoints,nsoil_sb))
    ALLOCATE(sb_pool%dom%nitrogen(npoints,nsoil_sb))
    ALLOCATE(sb_pool%dom%phosphorus(npoints,nsoil_sb))
    ALLOCATE(sb_pool%dom%carbon13(npoints,nsoil_sb))
    ALLOCATE(sb_pool%dom%carbon14(npoints,nsoil_sb))
    ALLOCATE(sb_pool%dom%nitrogen15(npoints,nsoil_sb))
    ALLOCATE(sb_pool%dom_assoc%carbon(npoints,nsoil_sb))
    ALLOCATE(sb_pool%dom_assoc%nitrogen(npoints,nsoil_sb))
    ALLOCATE(sb_pool%dom_assoc%phosphorus(npoints,nsoil_sb))
    ALLOCATE(sb_pool%dom_assoc%carbon13(npoints,nsoil_sb))
    ALLOCATE(sb_pool%dom_assoc%carbon14(npoints,nsoil_sb))
    ALLOCATE(sb_pool%dom_assoc%nitrogen15(npoints,nsoil_sb))
    ALLOCATE(sb_pool%soluable_litter%carbon(npoints,nsoil_sb))
    ALLOCATE(sb_pool%soluable_litter%nitrogen(npoints,nsoil_sb))
    ALLOCATE(sb_pool%soluable_litter%phosphorus(npoints,nsoil_sb))
    ALLOCATE(sb_pool%soluable_litter%carbon13(npoints,nsoil_sb))
    ALLOCATE(sb_pool%soluable_litter%carbon14(npoints,nsoil_sb))
    ALLOCATE(sb_pool%soluable_litter%nitrogen15(npoints,nsoil_sb))
    ALLOCATE(sb_pool%polymeric_litter%carbon(npoints,nsoil_sb))
    ALLOCATE(sb_pool%polymeric_litter%nitrogen(npoints,nsoil_sb))
    ALLOCATE(sb_pool%polymeric_litter%phosphorus(npoints,nsoil_sb))
    ALLOCATE(sb_pool%polymeric_litter%carbon13(npoints,nsoil_sb))
    ALLOCATE(sb_pool%polymeric_litter%carbon14(npoints,nsoil_sb))
    ALLOCATE(sb_pool%polymeric_litter%nitrogen15(npoints,nsoil_sb))
    ALLOCATE(sb_pool%woody_litter%carbon(npoints,nsoil_sb))
    ALLOCATE(sb_pool%woody_litter%nitrogen(npoints,nsoil_sb))
    ALLOCATE(sb_pool%woody_litter%phosphorus(npoints,nsoil_sb))
    ALLOCATE(sb_pool%woody_litter%carbon13(npoints,nsoil_sb))
    ALLOCATE(sb_pool%woody_litter%carbon14(npoints,nsoil_sb))
    ALLOCATE(sb_pool%woody_litter%nitrogen15(npoints,nsoil_sb))
    ALLOCATE(sb_pool%fungi%carbon(npoints,nsoil_sb))
    ALLOCATE(sb_pool%fungi%nitrogen(npoints,nsoil_sb))
    ALLOCATE(sb_pool%fungi%phosphorus(npoints,nsoil_sb))
    ALLOCATE(sb_pool%fungi%carbon13(npoints,nsoil_sb))
    ALLOCATE(sb_pool%fungi%carbon14(npoints,nsoil_sb))
    ALLOCATE(sb_pool%fungi%nitrogen15(npoints,nsoil_sb))
    ALLOCATE(sb_pool%mycorrhiza%carbon(npoints,nsoil_sb))
    ALLOCATE(sb_pool%mycorrhiza%nitrogen(npoints,nsoil_sb))
    ALLOCATE(sb_pool%mycorrhiza%phosphorus(npoints,nsoil_sb))
    ALLOCATE(sb_pool%mycorrhiza%carbon13(npoints,nsoil_sb))
    ALLOCATE(sb_pool%mycorrhiza%carbon14(npoints,nsoil_sb))
    ALLOCATE(sb_pool%mycorrhiza%nitrogen15(npoints,nsoil_sb))
    ALLOCATE(sb_pool%microbial%carbon(npoints,nsoil_sb))
    ALLOCATE(sb_pool%microbial%nitrogen(npoints,nsoil_sb))
    ALLOCATE(sb_pool%microbial%phosphorus(npoints,nsoil_sb))
    ALLOCATE(sb_pool%microbial%carbon13(npoints,nsoil_sb))
    ALLOCATE(sb_pool%microbial%carbon14(npoints,nsoil_sb))
    ALLOCATE(sb_pool%microbial%nitrogen15(npoints,nsoil_sb))
    ALLOCATE(sb_pool%residue%carbon(npoints,nsoil_sb))
    ALLOCATE(sb_pool%residue%nitrogen(npoints,nsoil_sb))
    ALLOCATE(sb_pool%residue%phosphorus(npoints,nsoil_sb))
    ALLOCATE(sb_pool%residue%carbon13(npoints,nsoil_sb))
    ALLOCATE(sb_pool%residue%carbon14(npoints,nsoil_sb))
    ALLOCATE(sb_pool%residue%nitrogen15(npoints,nsoil_sb))
    ALLOCATE(sb_pool%residue_assoc%carbon(npoints,nsoil_sb))
    ALLOCATE(sb_pool%residue_assoc%nitrogen(npoints,nsoil_sb))
    ALLOCATE(sb_pool%residue_assoc%phosphorus(npoints,nsoil_sb))
    ALLOCATE(sb_pool%residue_assoc%carbon13(npoints,nsoil_sb))
    ALLOCATE(sb_pool%residue_assoc%carbon14(npoints,nsoil_sb))
    ALLOCATE(sb_pool%residue_assoc%nitrogen15(npoints,nsoil_sb))

    ALLOCATE(sb_formation%dom%carbon(npoints,nsoil_sb))
    ALLOCATE(sb_formation%dom%nitrogen(npoints,nsoil_sb))
    ALLOCATE(sb_formation%dom%phosphorus(npoints,nsoil_sb))
    ALLOCATE(sb_formation%dom%carbon13(npoints,nsoil_sb))
    ALLOCATE(sb_formation%dom%carbon14(npoints,nsoil_sb))
    ALLOCATE(sb_formation%dom%nitrogen15(npoints,nsoil_sb))
    ALLOCATE(sb_formation%dom_assoc%carbon(npoints,nsoil_sb))
    ALLOCATE(sb_formation%dom_assoc%nitrogen(npoints,nsoil_sb))
    ALLOCATE(sb_formation%dom_assoc%phosphorus(npoints,nsoil_sb))
    ALLOCATE(sb_formation%dom_assoc%carbon13(npoints,nsoil_sb))
    ALLOCATE(sb_formation%dom_assoc%carbon14(npoints,nsoil_sb))
    ALLOCATE(sb_formation%dom_assoc%nitrogen15(npoints,nsoil_sb))
    ALLOCATE(sb_formation%soluable_litter%carbon(npoints,nsoil_sb))
    ALLOCATE(sb_formation%soluable_litter%nitrogen(npoints,nsoil_sb))
    ALLOCATE(sb_formation%soluable_litter%phosphorus(npoints,nsoil_sb))
    ALLOCATE(sb_formation%soluable_litter%carbon13(npoints,nsoil_sb))
    ALLOCATE(sb_formation%soluable_litter%carbon14(npoints,nsoil_sb))
    ALLOCATE(sb_formation%soluable_litter%nitrogen15(npoints,nsoil_sb))
    ALLOCATE(sb_formation%polymeric_litter%carbon(npoints,nsoil_sb))
    ALLOCATE(sb_formation%polymeric_litter%nitrogen(npoints,nsoil_sb))
    ALLOCATE(sb_formation%polymeric_litter%phosphorus(npoints,nsoil_sb))
    ALLOCATE(sb_formation%polymeric_litter%carbon13(npoints,nsoil_sb))
    ALLOCATE(sb_formation%polymeric_litter%carbon14(npoints,nsoil_sb))
    ALLOCATE(sb_formation%polymeric_litter%nitrogen15(npoints,nsoil_sb))
    ALLOCATE(sb_formation%woody_litter%carbon(npoints,nsoil_sb))
    ALLOCATE(sb_formation%woody_litter%nitrogen(npoints,nsoil_sb))
    ALLOCATE(sb_formation%woody_litter%phosphorus(npoints,nsoil_sb))
    ALLOCATE(sb_formation%woody_litter%carbon13(npoints,nsoil_sb))
    ALLOCATE(sb_formation%woody_litter%carbon14(npoints,nsoil_sb))
    ALLOCATE(sb_formation%woody_litter%nitrogen15(npoints,nsoil_sb))
    ALLOCATE(sb_formation%fungi%carbon(npoints,nsoil_sb))
    ALLOCATE(sb_formation%fungi%nitrogen(npoints,nsoil_sb))
    ALLOCATE(sb_formation%fungi%phosphorus(npoints,nsoil_sb))
    ALLOCATE(sb_formation%fungi%carbon13(npoints,nsoil_sb))
    ALLOCATE(sb_formation%fungi%carbon14(npoints,nsoil_sb))
    ALLOCATE(sb_formation%fungi%nitrogen15(npoints,nsoil_sb))
    ALLOCATE(sb_formation%mycorrhiza%carbon(npoints,nsoil_sb))
    ALLOCATE(sb_formation%mycorrhiza%nitrogen(npoints,nsoil_sb))
    ALLOCATE(sb_formation%mycorrhiza%phosphorus(npoints,nsoil_sb))
    ALLOCATE(sb_formation%mycorrhiza%carbon13(npoints,nsoil_sb))
    ALLOCATE(sb_formation%mycorrhiza%carbon14(npoints,nsoil_sb))
    ALLOCATE(sb_formation%mycorrhiza%nitrogen15(npoints,nsoil_sb))
    ALLOCATE(sb_formation%microbial%carbon(npoints,nsoil_sb))
    ALLOCATE(sb_formation%microbial%nitrogen(npoints,nsoil_sb))
    ALLOCATE(sb_formation%microbial%phosphorus(npoints,nsoil_sb))
    ALLOCATE(sb_formation%microbial%carbon13(npoints,nsoil_sb))
    ALLOCATE(sb_formation%microbial%carbon14(npoints,nsoil_sb))
    ALLOCATE(sb_formation%microbial%nitrogen15(npoints,nsoil_sb))
    ALLOCATE(sb_formation%residue%carbon(npoints,nsoil_sb))
    ALLOCATE(sb_formation%residue%nitrogen(npoints,nsoil_sb))
    ALLOCATE(sb_formation%residue%phosphorus(npoints,nsoil_sb))
    ALLOCATE(sb_formation%residue%carbon13(npoints,nsoil_sb))
    ALLOCATE(sb_formation%residue%carbon14(npoints,nsoil_sb))
    ALLOCATE(sb_formation%residue%nitrogen15(npoints,nsoil_sb))
    ALLOCATE(sb_formation%residue_assoc%carbon(npoints,nsoil_sb))
    ALLOCATE(sb_formation%residue_assoc%nitrogen(npoints,nsoil_sb))
    ALLOCATE(sb_formation%residue_assoc%phosphorus(npoints,nsoil_sb))
    ALLOCATE(sb_formation%residue_assoc%carbon13(npoints,nsoil_sb))
    ALLOCATE(sb_formation%residue_assoc%carbon14(npoints,nsoil_sb))
    ALLOCATE(sb_formation%residue_assoc%nitrogen15(npoints,nsoil_sb))

    ALLOCATE(sb_loss%dom%carbon(npoints,nsoil_sb))
    ALLOCATE(sb_loss%dom%nitrogen(npoints,nsoil_sb))
    ALLOCATE(sb_loss%dom%phosphorus(npoints,nsoil_sb))
    ALLOCATE(sb_loss%dom%carbon13(npoints,nsoil_sb))
    ALLOCATE(sb_loss%dom%carbon14(npoints,nsoil_sb))
    ALLOCATE(sb_loss%dom%nitrogen15(npoints,nsoil_sb))
    ALLOCATE(sb_loss%dom_assoc%carbon(npoints,nsoil_sb))
    ALLOCATE(sb_loss%dom_assoc%nitrogen(npoints,nsoil_sb))
    ALLOCATE(sb_loss%dom_assoc%phosphorus(npoints,nsoil_sb))
    ALLOCATE(sb_loss%dom_assoc%carbon13(npoints,nsoil_sb))
    ALLOCATE(sb_loss%dom_assoc%carbon14(npoints,nsoil_sb))
    ALLOCATE(sb_loss%dom_assoc%nitrogen15(npoints,nsoil_sb))
    ALLOCATE(sb_loss%soluable_litter%carbon(npoints,nsoil_sb))
    ALLOCATE(sb_loss%soluable_litter%nitrogen(npoints,nsoil_sb))
    ALLOCATE(sb_loss%soluable_litter%phosphorus(npoints,nsoil_sb))
    ALLOCATE(sb_loss%soluable_litter%carbon13(npoints,nsoil_sb))
    ALLOCATE(sb_loss%soluable_litter%carbon14(npoints,nsoil_sb))
    ALLOCATE(sb_loss%soluable_litter%nitrogen15(npoints,nsoil_sb))
    ALLOCATE(sb_loss%polymeric_litter%carbon(npoints,nsoil_sb))
    ALLOCATE(sb_loss%polymeric_litter%nitrogen(npoints,nsoil_sb))
    ALLOCATE(sb_loss%polymeric_litter%phosphorus(npoints,nsoil_sb))
    ALLOCATE(sb_loss%polymeric_litter%carbon13(npoints,nsoil_sb))
    ALLOCATE(sb_loss%polymeric_litter%carbon14(npoints,nsoil_sb))
    ALLOCATE(sb_loss%polymeric_litter%nitrogen15(npoints,nsoil_sb))
    ALLOCATE(sb_loss%woody_litter%carbon(npoints,nsoil_sb))
    ALLOCATE(sb_loss%woody_litter%nitrogen(npoints,nsoil_sb))
    ALLOCATE(sb_loss%woody_litter%phosphorus(npoints,nsoil_sb))
    ALLOCATE(sb_loss%woody_litter%carbon13(npoints,nsoil_sb))
    ALLOCATE(sb_loss%woody_litter%carbon14(npoints,nsoil_sb))
    ALLOCATE(sb_loss%woody_litter%nitrogen15(npoints,nsoil_sb))
    ALLOCATE(sb_loss%fungi%carbon(npoints,nsoil_sb))
    ALLOCATE(sb_loss%fungi%nitrogen(npoints,nsoil_sb))
    ALLOCATE(sb_loss%fungi%phosphorus(npoints,nsoil_sb))
    ALLOCATE(sb_loss%fungi%carbon13(npoints,nsoil_sb))
    ALLOCATE(sb_loss%fungi%carbon14(npoints,nsoil_sb))
    ALLOCATE(sb_loss%fungi%nitrogen15(npoints,nsoil_sb))
    ALLOCATE(sb_loss%mycorrhiza%carbon(npoints,nsoil_sb))
    ALLOCATE(sb_loss%mycorrhiza%nitrogen(npoints,nsoil_sb))
    ALLOCATE(sb_loss%mycorrhiza%phosphorus(npoints,nsoil_sb))
    ALLOCATE(sb_loss%mycorrhiza%carbon13(npoints,nsoil_sb))
    ALLOCATE(sb_loss%mycorrhiza%carbon14(npoints,nsoil_sb))
    ALLOCATE(sb_loss%mycorrhiza%nitrogen15(npoints,nsoil_sb))
    ALLOCATE(sb_loss%microbial%carbon(npoints,nsoil_sb))
    ALLOCATE(sb_loss%microbial%nitrogen(npoints,nsoil_sb))
    ALLOCATE(sb_loss%microbial%phosphorus(npoints,nsoil_sb))
    ALLOCATE(sb_loss%microbial%carbon13(npoints,nsoil_sb))
    ALLOCATE(sb_loss%microbial%carbon14(npoints,nsoil_sb))
    ALLOCATE(sb_loss%microbial%nitrogen15(npoints,nsoil_sb))
    ALLOCATE(sb_loss%residue%carbon(npoints,nsoil_sb))
    ALLOCATE(sb_loss%residue%nitrogen(npoints,nsoil_sb))
    ALLOCATE(sb_loss%residue%phosphorus(npoints,nsoil_sb))
    ALLOCATE(sb_loss%residue%carbon13(npoints,nsoil_sb))
    ALLOCATE(sb_loss%residue%carbon14(npoints,nsoil_sb))
    ALLOCATE(sb_loss%residue%nitrogen15(npoints,nsoil_sb))
    ALLOCATE(sb_loss%residue_assoc%carbon(npoints,nsoil_sb))
    ALLOCATE(sb_loss%residue_assoc%nitrogen(npoints,nsoil_sb))
    ALLOCATE(sb_loss%residue_assoc%phosphorus(npoints,nsoil_sb))
    ALLOCATE(sb_loss%residue_assoc%carbon13(npoints,nsoil_sb))
    ALLOCATE(sb_loss%residue_assoc%carbon14(npoints,nsoil_sb))
    ALLOCATE(sb_loss%residue_assoc%nitrogen15(npoints,nsoil_sb))

    ALLOCATE(sb_transport%dom%carbon(npoints,nsoil_sb))
    ALLOCATE(sb_transport%dom%nitrogen(npoints,nsoil_sb))
    ALLOCATE(sb_transport%dom%phosphorus(npoints,nsoil_sb))
    ALLOCATE(sb_transport%dom%carbon13(npoints,nsoil_sb))
    ALLOCATE(sb_transport%dom%carbon14(npoints,nsoil_sb))
    ALLOCATE(sb_transport%dom%nitrogen15(npoints,nsoil_sb))
    ALLOCATE(sb_transport%dom_assoc%carbon(npoints,nsoil_sb))
    ALLOCATE(sb_transport%dom_assoc%nitrogen(npoints,nsoil_sb))
    ALLOCATE(sb_transport%dom_assoc%phosphorus(npoints,nsoil_sb))
    ALLOCATE(sb_transport%dom_assoc%carbon13(npoints,nsoil_sb))
    ALLOCATE(sb_transport%dom_assoc%carbon14(npoints,nsoil_sb))
    ALLOCATE(sb_transport%dom_assoc%nitrogen15(npoints,nsoil_sb))
    ALLOCATE(sb_transport%soluable_litter%carbon(npoints,nsoil_sb))
    ALLOCATE(sb_transport%soluable_litter%nitrogen(npoints,nsoil_sb))
    ALLOCATE(sb_transport%soluable_litter%phosphorus(npoints,nsoil_sb))
    ALLOCATE(sb_transport%soluable_litter%carbon13(npoints,nsoil_sb))
    ALLOCATE(sb_transport%soluable_litter%carbon14(npoints,nsoil_sb))
    ALLOCATE(sb_transport%soluable_litter%nitrogen15(npoints,nsoil_sb))
    ALLOCATE(sb_transport%polymeric_litter%carbon(npoints,nsoil_sb))
    ALLOCATE(sb_transport%polymeric_litter%nitrogen(npoints,nsoil_sb))
    ALLOCATE(sb_transport%polymeric_litter%phosphorus(npoints,nsoil_sb))
    ALLOCATE(sb_transport%polymeric_litter%carbon13(npoints,nsoil_sb))
    ALLOCATE(sb_transport%polymeric_litter%carbon14(npoints,nsoil_sb))
    ALLOCATE(sb_transport%polymeric_litter%nitrogen15(npoints,nsoil_sb))
    ALLOCATE(sb_transport%woody_litter%carbon(npoints,nsoil_sb))
    ALLOCATE(sb_transport%woody_litter%nitrogen(npoints,nsoil_sb))
    ALLOCATE(sb_transport%woody_litter%phosphorus(npoints,nsoil_sb))
    ALLOCATE(sb_transport%woody_litter%carbon13(npoints,nsoil_sb))
    ALLOCATE(sb_transport%woody_litter%carbon14(npoints,nsoil_sb))
    ALLOCATE(sb_transport%woody_litter%nitrogen15(npoints,nsoil_sb))
    ALLOCATE(sb_transport%fungi%carbon(npoints,nsoil_sb))
    ALLOCATE(sb_transport%fungi%nitrogen(npoints,nsoil_sb))
    ALLOCATE(sb_transport%fungi%phosphorus(npoints,nsoil_sb))
    ALLOCATE(sb_transport%fungi%carbon13(npoints,nsoil_sb))
    ALLOCATE(sb_transport%fungi%carbon14(npoints,nsoil_sb))
    ALLOCATE(sb_transport%fungi%nitrogen15(npoints,nsoil_sb))
    ALLOCATE(sb_transport%mycorrhiza%carbon(npoints,nsoil_sb))
    ALLOCATE(sb_transport%mycorrhiza%nitrogen(npoints,nsoil_sb))
    ALLOCATE(sb_transport%mycorrhiza%phosphorus(npoints,nsoil_sb))
    ALLOCATE(sb_transport%mycorrhiza%carbon13(npoints,nsoil_sb))
    ALLOCATE(sb_transport%mycorrhiza%carbon14(npoints,nsoil_sb))
    ALLOCATE(sb_transport%mycorrhiza%nitrogen15(npoints,nsoil_sb))
    ALLOCATE(sb_transport%microbial%carbon(npoints,nsoil_sb))
    ALLOCATE(sb_transport%microbial%nitrogen(npoints,nsoil_sb))
    ALLOCATE(sb_transport%microbial%phosphorus(npoints,nsoil_sb))
    ALLOCATE(sb_transport%microbial%carbon13(npoints,nsoil_sb))
    ALLOCATE(sb_transport%microbial%carbon14(npoints,nsoil_sb))
    ALLOCATE(sb_transport%microbial%nitrogen15(npoints,nsoil_sb))
    ALLOCATE(sb_transport%residue%carbon(npoints,nsoil_sb))
    ALLOCATE(sb_transport%residue%nitrogen(npoints,nsoil_sb))
    ALLOCATE(sb_transport%residue%phosphorus(npoints,nsoil_sb))
    ALLOCATE(sb_transport%residue%carbon13(npoints,nsoil_sb))
    ALLOCATE(sb_transport%residue%carbon14(npoints,nsoil_sb))
    ALLOCATE(sb_transport%residue%nitrogen15(npoints,nsoil_sb))
    ALLOCATE(sb_transport%residue_assoc%carbon(npoints,nsoil_sb))
    ALLOCATE(sb_transport%residue_assoc%nitrogen(npoints,nsoil_sb))
    ALLOCATE(sb_transport%residue_assoc%phosphorus(npoints,nsoil_sb))
    ALLOCATE(sb_transport%residue_assoc%carbon13(npoints,nsoil_sb))
    ALLOCATE(sb_transport%residue_assoc%carbon14(npoints,nsoil_sb))
    ALLOCATE(sb_transport%residue_assoc%nitrogen15(npoints,nsoil_sb))

    ALLOCATE(sb_delta_pools%carbon(npoints,nsoil_sb))
    ALLOCATE(sb_delta_pools%nitrogen(npoints,nsoil_sb))
    ALLOCATE(sb_delta_pools%phosphorus(npoints,nsoil_sb))
    ALLOCATE(sb_delta_pools%carbon13(npoints,nsoil_sb))
    ALLOCATE(sb_delta_pools%carbon14(npoints,nsoil_sb))
    ALLOCATE(sb_delta_pools%nitrogen15(npoints,nsoil_sb))

    ALLOCATE(sb_mycorrhiza_export%carbon(npoints,nsoil_sb))
    ALLOCATE(sb_mycorrhiza_export%nitrogen(npoints,nsoil_sb))
    ALLOCATE(sb_mycorrhiza_export%phosphorus(npoints,nsoil_sb))
    ALLOCATE(sb_mycorrhiza_export%carbon13(npoints,nsoil_sb))
    ALLOCATE(sb_mycorrhiza_export%carbon14(npoints,nsoil_sb))
    ALLOCATE(sb_mycorrhiza_export%nitrogen15(npoints,nsoil_sb))
  END SUBROUTINE Init_sb_pools_elements_a


  !-----------------------------------------------------------------------------------------------------
  !> DEALLOCATE variables from sb_pool_a and sb_element_a structures at the program end
  !!
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE Deallocate_var_sb_pools_elements_a

    IMPLICIT NONE

    !< 1.0  DEALLOCATE variables
    !! 
    DEALLOCATE(sb_pool%dom%carbon)
    DEALLOCATE(sb_pool%dom%nitrogen)
    DEALLOCATE(sb_pool%dom%phosphorus)
    DEALLOCATE(sb_pool%dom%carbon13)
    DEALLOCATE(sb_pool%dom%carbon14)
    DEALLOCATE(sb_pool%dom%nitrogen15)
    DEALLOCATE(sb_pool%dom_assoc%carbon)
    DEALLOCATE(sb_pool%dom_assoc%nitrogen)
    DEALLOCATE(sb_pool%dom_assoc%phosphorus)
    DEALLOCATE(sb_pool%dom_assoc%carbon13)
    DEALLOCATE(sb_pool%dom_assoc%carbon14)
    DEALLOCATE(sb_pool%dom_assoc%nitrogen15)
    DEALLOCATE(sb_pool%soluable_litter%carbon)
    DEALLOCATE(sb_pool%soluable_litter%nitrogen)
    DEALLOCATE(sb_pool%soluable_litter%phosphorus)
    DEALLOCATE(sb_pool%soluable_litter%carbon13)
    DEALLOCATE(sb_pool%soluable_litter%carbon14)
    DEALLOCATE(sb_pool%soluable_litter%nitrogen15)
    DEALLOCATE(sb_pool%polymeric_litter%carbon)
    DEALLOCATE(sb_pool%polymeric_litter%nitrogen)
    DEALLOCATE(sb_pool%polymeric_litter%phosphorus)
    DEALLOCATE(sb_pool%polymeric_litter%carbon13)
    DEALLOCATE(sb_pool%polymeric_litter%carbon14)
    DEALLOCATE(sb_pool%polymeric_litter%nitrogen15)
    DEALLOCATE(sb_pool%woody_litter%carbon)
    DEALLOCATE(sb_pool%woody_litter%nitrogen)
    DEALLOCATE(sb_pool%woody_litter%phosphorus)
    DEALLOCATE(sb_pool%woody_litter%carbon13)
    DEALLOCATE(sb_pool%woody_litter%carbon14)
    DEALLOCATE(sb_pool%woody_litter%nitrogen15)
    DEALLOCATE(sb_pool%fungi%carbon)
    DEALLOCATE(sb_pool%fungi%nitrogen)
    DEALLOCATE(sb_pool%fungi%phosphorus)
    DEALLOCATE(sb_pool%fungi%carbon13)
    DEALLOCATE(sb_pool%fungi%carbon14)
    DEALLOCATE(sb_pool%fungi%nitrogen15)
    DEALLOCATE(sb_pool%mycorrhiza%carbon)
    DEALLOCATE(sb_pool%mycorrhiza%nitrogen)
    DEALLOCATE(sb_pool%mycorrhiza%phosphorus)
    DEALLOCATE(sb_pool%mycorrhiza%carbon13)
    DEALLOCATE(sb_pool%mycorrhiza%carbon14)
    DEALLOCATE(sb_pool%mycorrhiza%nitrogen15)
    DEALLOCATE(sb_pool%microbial%carbon)
    DEALLOCATE(sb_pool%microbial%nitrogen)
    DEALLOCATE(sb_pool%microbial%phosphorus)
    DEALLOCATE(sb_pool%microbial%carbon13)
    DEALLOCATE(sb_pool%microbial%carbon14)
    DEALLOCATE(sb_pool%microbial%nitrogen15)
    DEALLOCATE(sb_pool%residue%carbon)
    DEALLOCATE(sb_pool%residue%nitrogen)
    DEALLOCATE(sb_pool%residue%phosphorus)
    DEALLOCATE(sb_pool%residue%carbon13)
    DEALLOCATE(sb_pool%residue%carbon14)
    DEALLOCATE(sb_pool%residue%nitrogen15)
    DEALLOCATE(sb_pool%residue_assoc%carbon)
    DEALLOCATE(sb_pool%residue_assoc%nitrogen)
    DEALLOCATE(sb_pool%residue_assoc%phosphorus)
    DEALLOCATE(sb_pool%residue_assoc%carbon13)
    DEALLOCATE(sb_pool%residue_assoc%carbon14)
    DEALLOCATE(sb_pool%residue_assoc%nitrogen15)

    DEALLOCATE(sb_formation%dom%carbon)
    DEALLOCATE(sb_formation%dom%nitrogen)
    DEALLOCATE(sb_formation%dom%phosphorus)
    DEALLOCATE(sb_formation%dom%carbon13)
    DEALLOCATE(sb_formation%dom%carbon14)
    DEALLOCATE(sb_formation%dom%nitrogen15)
    DEALLOCATE(sb_formation%dom_assoc%carbon)
    DEALLOCATE(sb_formation%dom_assoc%nitrogen)
    DEALLOCATE(sb_formation%dom_assoc%phosphorus)
    DEALLOCATE(sb_formation%dom_assoc%carbon13)
    DEALLOCATE(sb_formation%dom_assoc%carbon14)
    DEALLOCATE(sb_formation%dom_assoc%nitrogen15)
    DEALLOCATE(sb_formation%soluable_litter%carbon)
    DEALLOCATE(sb_formation%soluable_litter%nitrogen)
    DEALLOCATE(sb_formation%soluable_litter%phosphorus)
    DEALLOCATE(sb_formation%soluable_litter%carbon13)
    DEALLOCATE(sb_formation%soluable_litter%carbon14)
    DEALLOCATE(sb_formation%soluable_litter%nitrogen15)
    DEALLOCATE(sb_formation%polymeric_litter%carbon)
    DEALLOCATE(sb_formation%polymeric_litter%nitrogen)
    DEALLOCATE(sb_formation%polymeric_litter%phosphorus)
    DEALLOCATE(sb_formation%polymeric_litter%carbon13)
    DEALLOCATE(sb_formation%polymeric_litter%carbon14)
    DEALLOCATE(sb_formation%polymeric_litter%nitrogen15)
    DEALLOCATE(sb_formation%woody_litter%carbon)
    DEALLOCATE(sb_formation%woody_litter%nitrogen)
    DEALLOCATE(sb_formation%woody_litter%phosphorus)
    DEALLOCATE(sb_formation%woody_litter%carbon13)
    DEALLOCATE(sb_formation%woody_litter%carbon14)
    DEALLOCATE(sb_formation%woody_litter%nitrogen15)
    DEALLOCATE(sb_formation%fungi%carbon)
    DEALLOCATE(sb_formation%fungi%nitrogen)
    DEALLOCATE(sb_formation%fungi%phosphorus)
    DEALLOCATE(sb_formation%fungi%carbon13)
    DEALLOCATE(sb_formation%fungi%carbon14)
    DEALLOCATE(sb_formation%fungi%nitrogen15)
    DEALLOCATE(sb_formation%mycorrhiza%carbon)
    DEALLOCATE(sb_formation%mycorrhiza%nitrogen)
    DEALLOCATE(sb_formation%mycorrhiza%phosphorus)
    DEALLOCATE(sb_formation%mycorrhiza%carbon13)
    DEALLOCATE(sb_formation%mycorrhiza%carbon14)
    DEALLOCATE(sb_formation%mycorrhiza%nitrogen15)
    DEALLOCATE(sb_formation%microbial%carbon)
    DEALLOCATE(sb_formation%microbial%nitrogen)
    DEALLOCATE(sb_formation%microbial%phosphorus)
    DEALLOCATE(sb_formation%microbial%carbon13)
    DEALLOCATE(sb_formation%microbial%carbon14)
    DEALLOCATE(sb_formation%microbial%nitrogen15)
    DEALLOCATE(sb_formation%residue%carbon)
    DEALLOCATE(sb_formation%residue%nitrogen)
    DEALLOCATE(sb_formation%residue%phosphorus)
    DEALLOCATE(sb_formation%residue%carbon13)
    DEALLOCATE(sb_formation%residue%carbon14)
    DEALLOCATE(sb_formation%residue%nitrogen15)
    DEALLOCATE(sb_formation%residue_assoc%carbon)
    DEALLOCATE(sb_formation%residue_assoc%nitrogen)
    DEALLOCATE(sb_formation%residue_assoc%phosphorus)
    DEALLOCATE(sb_formation%residue_assoc%carbon13)
    DEALLOCATE(sb_formation%residue_assoc%carbon14)
    DEALLOCATE(sb_formation%residue_assoc%nitrogen15)

    DEALLOCATE(sb_loss%dom%carbon)
    DEALLOCATE(sb_loss%dom%nitrogen)
    DEALLOCATE(sb_loss%dom%phosphorus)
    DEALLOCATE(sb_loss%dom%carbon13)
    DEALLOCATE(sb_loss%dom%carbon14)
    DEALLOCATE(sb_loss%dom%nitrogen15)
    DEALLOCATE(sb_loss%dom_assoc%carbon)
    DEALLOCATE(sb_loss%dom_assoc%nitrogen)
    DEALLOCATE(sb_loss%dom_assoc%phosphorus)
    DEALLOCATE(sb_loss%dom_assoc%carbon13)
    DEALLOCATE(sb_loss%dom_assoc%carbon14)
    DEALLOCATE(sb_loss%dom_assoc%nitrogen15)
    DEALLOCATE(sb_loss%soluable_litter%carbon)
    DEALLOCATE(sb_loss%soluable_litter%nitrogen)
    DEALLOCATE(sb_loss%soluable_litter%phosphorus)
    DEALLOCATE(sb_loss%soluable_litter%carbon13)
    DEALLOCATE(sb_loss%soluable_litter%carbon14)
    DEALLOCATE(sb_loss%soluable_litter%nitrogen15)
    DEALLOCATE(sb_loss%polymeric_litter%carbon)
    DEALLOCATE(sb_loss%polymeric_litter%nitrogen)
    DEALLOCATE(sb_loss%polymeric_litter%phosphorus)
    DEALLOCATE(sb_loss%polymeric_litter%carbon13)
    DEALLOCATE(sb_loss%polymeric_litter%carbon14)
    DEALLOCATE(sb_loss%polymeric_litter%nitrogen15)
    DEALLOCATE(sb_loss%woody_litter%carbon)
    DEALLOCATE(sb_loss%woody_litter%nitrogen)
    DEALLOCATE(sb_loss%woody_litter%phosphorus)
    DEALLOCATE(sb_loss%woody_litter%carbon13)
    DEALLOCATE(sb_loss%woody_litter%carbon14)
    DEALLOCATE(sb_loss%woody_litter%nitrogen15)
    DEALLOCATE(sb_loss%fungi%carbon)
    DEALLOCATE(sb_loss%fungi%nitrogen)
    DEALLOCATE(sb_loss%fungi%phosphorus)
    DEALLOCATE(sb_loss%fungi%carbon13)
    DEALLOCATE(sb_loss%fungi%carbon14)
    DEALLOCATE(sb_loss%fungi%nitrogen15)
    DEALLOCATE(sb_loss%mycorrhiza%carbon)
    DEALLOCATE(sb_loss%mycorrhiza%nitrogen)
    DEALLOCATE(sb_loss%mycorrhiza%phosphorus)
    DEALLOCATE(sb_loss%mycorrhiza%carbon13)
    DEALLOCATE(sb_loss%mycorrhiza%carbon14)
    DEALLOCATE(sb_loss%mycorrhiza%nitrogen15)
    DEALLOCATE(sb_loss%microbial%carbon)
    DEALLOCATE(sb_loss%microbial%nitrogen)
    DEALLOCATE(sb_loss%microbial%phosphorus)
    DEALLOCATE(sb_loss%microbial%carbon13)
    DEALLOCATE(sb_loss%microbial%carbon14)
    DEALLOCATE(sb_loss%microbial%nitrogen15)
    DEALLOCATE(sb_loss%residue%carbon)
    DEALLOCATE(sb_loss%residue%nitrogen)
    DEALLOCATE(sb_loss%residue%phosphorus)
    DEALLOCATE(sb_loss%residue%carbon13)
    DEALLOCATE(sb_loss%residue%carbon14)
    DEALLOCATE(sb_loss%residue%nitrogen15)
    DEALLOCATE(sb_loss%residue_assoc%carbon)
    DEALLOCATE(sb_loss%residue_assoc%nitrogen)
    DEALLOCATE(sb_loss%residue_assoc%phosphorus)
    DEALLOCATE(sb_loss%residue_assoc%carbon13)
    DEALLOCATE(sb_loss%residue_assoc%carbon14)
    DEALLOCATE(sb_loss%residue_assoc%nitrogen15)

    DEALLOCATE(sb_transport%dom%carbon)
    DEALLOCATE(sb_transport%dom%nitrogen)
    DEALLOCATE(sb_transport%dom%phosphorus)
    DEALLOCATE(sb_transport%dom%carbon13)
    DEALLOCATE(sb_transport%dom%carbon14)
    DEALLOCATE(sb_transport%dom%nitrogen15)
    DEALLOCATE(sb_transport%dom_assoc%carbon)
    DEALLOCATE(sb_transport%dom_assoc%nitrogen)
    DEALLOCATE(sb_transport%dom_assoc%phosphorus)
    DEALLOCATE(sb_transport%dom_assoc%carbon13)
    DEALLOCATE(sb_transport%dom_assoc%carbon14)
    DEALLOCATE(sb_transport%dom_assoc%nitrogen15)
    DEALLOCATE(sb_transport%soluable_litter%carbon)
    DEALLOCATE(sb_transport%soluable_litter%nitrogen)
    DEALLOCATE(sb_transport%soluable_litter%phosphorus)
    DEALLOCATE(sb_transport%soluable_litter%carbon13)
    DEALLOCATE(sb_transport%soluable_litter%carbon14)
    DEALLOCATE(sb_transport%soluable_litter%nitrogen15)
    DEALLOCATE(sb_transport%polymeric_litter%carbon)
    DEALLOCATE(sb_transport%polymeric_litter%nitrogen)
    DEALLOCATE(sb_transport%polymeric_litter%phosphorus)
    DEALLOCATE(sb_transport%polymeric_litter%carbon13)
    DEALLOCATE(sb_transport%polymeric_litter%carbon14)
    DEALLOCATE(sb_transport%polymeric_litter%nitrogen15)
    DEALLOCATE(sb_transport%woody_litter%carbon)
    DEALLOCATE(sb_transport%woody_litter%nitrogen)
    DEALLOCATE(sb_transport%woody_litter%phosphorus)
    DEALLOCATE(sb_transport%woody_litter%carbon13)
    DEALLOCATE(sb_transport%woody_litter%carbon14)
    DEALLOCATE(sb_transport%woody_litter%nitrogen15)
    DEALLOCATE(sb_transport%fungi%carbon)
    DEALLOCATE(sb_transport%fungi%nitrogen)
    DEALLOCATE(sb_transport%fungi%phosphorus)
    DEALLOCATE(sb_transport%fungi%carbon13)
    DEALLOCATE(sb_transport%fungi%carbon14)
    DEALLOCATE(sb_transport%fungi%nitrogen15)
    DEALLOCATE(sb_transport%mycorrhiza%carbon)
    DEALLOCATE(sb_transport%mycorrhiza%nitrogen)
    DEALLOCATE(sb_transport%mycorrhiza%phosphorus)
    DEALLOCATE(sb_transport%mycorrhiza%carbon13)
    DEALLOCATE(sb_transport%mycorrhiza%carbon14)
    DEALLOCATE(sb_transport%mycorrhiza%nitrogen15)
    DEALLOCATE(sb_transport%microbial%carbon)
    DEALLOCATE(sb_transport%microbial%nitrogen)
    DEALLOCATE(sb_transport%microbial%phosphorus)
    DEALLOCATE(sb_transport%microbial%carbon13)
    DEALLOCATE(sb_transport%microbial%carbon14)
    DEALLOCATE(sb_transport%microbial%nitrogen15)
    DEALLOCATE(sb_transport%residue%carbon)
    DEALLOCATE(sb_transport%residue%nitrogen)
    DEALLOCATE(sb_transport%residue%phosphorus)
    DEALLOCATE(sb_transport%residue%carbon13)
    DEALLOCATE(sb_transport%residue%carbon14)
    DEALLOCATE(sb_transport%residue%nitrogen15)
    DEALLOCATE(sb_transport%residue_assoc%carbon)
    DEALLOCATE(sb_transport%residue_assoc%nitrogen)
    DEALLOCATE(sb_transport%residue_assoc%phosphorus)
    DEALLOCATE(sb_transport%residue_assoc%carbon13)
    DEALLOCATE(sb_transport%residue_assoc%carbon14)
    DEALLOCATE(sb_transport%residue_assoc%nitrogen15)

    DEALLOCATE(sb_delta_pools%carbon)
    DEALLOCATE(sb_delta_pools%nitrogen)
    DEALLOCATE(sb_delta_pools%phosphorus)
    DEALLOCATE(sb_delta_pools%carbon13)
    DEALLOCATE(sb_delta_pools%carbon14)
    DEALLOCATE(sb_delta_pools%nitrogen15)

    DEALLOCATE(sb_mycorrhiza_export%carbon)
    DEALLOCATE(sb_mycorrhiza_export%nitrogen)
    DEALLOCATE(sb_mycorrhiza_export%phosphorus)
    DEALLOCATE(sb_mycorrhiza_export%carbon13)
    DEALLOCATE(sb_mycorrhiza_export%carbon14)
    DEALLOCATE(sb_mycorrhiza_export%nitrogen15)
  END SUBROUTINE Deallocate_var_sb_pools_elements_a

END MODULE mo_sb_memory_class
