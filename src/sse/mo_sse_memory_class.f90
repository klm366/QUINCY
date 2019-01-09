!>
!! @brief Contains the memory class for the sse processes
!!
!! @author
!!  Jan Engel, MPI-BGC Jena (based on the JSBACH4 version from Rainer Schneck, MPI-M Hamburg)
!!
!! @par Revision History
!! First version by Jan Engel (20160913)
!!
MODULE mo_sse_memory_class
#ifndef __NO_JSBACH__

  USE mo_kind,                   ONLY: wp
  USE mo_exception,              ONLY: message
  USE mo_jsb_class,              ONLY: Get_model
  USE mo_jsb_memory_class,       ONLY: t_jsb_memory
  USE mo_jsb_var_class,          ONLY: t_jsb_var_real2d, t_jsb_var_real3d


  IMPLICIT NONE
  PRIVATE
  PUBLIC :: t_sse_memory, max_no_of_vars

  INTEGER, PARAMETER :: max_no_of_vars      = 60 


  !-----------------------------------------------------------------------------------------------------
  !> Type definition for sse memory
  !!
  !! @par includes: \n
  !!    sse variables sse\%VAR
  !! 
  !-----------------------------------------------------------------------------------------------------
  TYPE, EXTENDS(t_jsb_memory) :: t_sse_memory


    TYPE(t_jsb_var_real2d)          :: &
                                    soil_depth      , & !< Soil depth derived from textures (bedrock) (for hydrology) [m]
                                    root_depth      , & !< actual rooting depth [m]
                                    t_srf_old       , & !< previous surface temperature [K]
                                    zril_old            !< previous timestep's Reynold's number

    ! JSBACH SOIL
    TYPE(t_jsb_var_real3d)          :: &
                                    soil_depth_sl       !< Depth of each soil layer that can be saturated with water (until bedrock) [m]

    ! JSBACH HYDRO MEMORY (temporary here to get rid of vegsoil construct)
    TYPE(t_jsb_var_real2d)          :: &
                                    interception,     & !< surface interception evaporation  [kg m-2 s-1]
                                    pet,              & !< potential evaporation  [kg m-2 s-1]
                                    evaporation,      & !< surface evaporation  [kg m-2 s-1]
                                    runoff,           & !< surface runoff  [kg m-2 s-1]
                                    drainage,         & !< drainage        [kg m-2 s-1]
                                    drainage_fraction   !< water mass fraction lost from the soil column by drainage [1 s-1] 

    ! Additional variables for sse
    TYPE(t_jsb_var_real2d)          :: &
                                    w_skin,           & !< Water content in skin reservoir                                          [m]
                                    w_soil_root,      & !< Water content in root zone of the soil                                   [m]
                                    w_soil_root_fc,   & !< Water content at field capacity in root zone of the soil                 [m]
                                    w_soil_root_pwp,  & !< Water content at permanent wilting point in root zone of the soil        [m]
                                    t_soil_root         !< Temperature of fine roots in the soil                                    [K]

    TYPE(t_jsb_var_real3d)          :: &
                                    percolation_sl,   & !< fraction of mass transferred to below layer                              [fraction s-1]
                                    drainage_sl,      & !< drainage of soil layers                                                  [kg m-2 s-1]
                                    w_soil_sl,        & !< Water content in soil layers                                             [m]
                                    w_soil_sat_sl,    & !< Water content of soil layers at saturation 
                                                        !! (derived from volumetric capacity at saturation)                         [m]
                                    w_soil_fc_sl,     & !< Water content of soil layers at field capacity
                                                        !! (derived from volumetric field capacity)                                 [m]
                                    w_soil_pwp_sl,    & !< Water content of soil layers at permanent wilting point
                                                        !! (derived from vol. perm wilt point)                                      [m]
                                    saxtonA,          & !< coefficient in moisture-tension relationship
                                    saxtonB,          & !< coefficient in moisture-tension relationship
                                    saxtonC,          & !< exponent of moisture-conductivity relationship
                                    kdiff_sat_sl        !< saturated hydraulic conductivity                                         [m s-1]

    ! transpiration
    TYPE(t_jsb_var_real2d)          :: &
                                    transpiration       !< transpiration in [kg m-2 s-1]

    !Soil Water addons for SLM in HYDRO
    TYPE(t_jsb_var_real2d)          :: &
                                    w_soil_root_theta, & !< root zone integrated plant available water, fraction of maximum 
                                    w_soil_root_pot      !< soil water potential [MPa]

    TYPE(t_jsb_var_real3d)          :: &
                                    w_soil_pot_sl        !< soil water potential per layer [MPa]

    ! soil energy 
    TYPE(t_jsb_var_real3d)          :: &
                                    t_soil_sl            !< soil temperature [K] 

    ! soil physical properties (not sure where they should belong to, but they are needed in the long-term)
    TYPE(t_jsb_var_real3d)          :: & 
                                    heat_capa_sl,         &  !< soil heat capacity [J kg-1]
                                    bulk_dens_sl,         &  !< soil bulk density [kg m-3]
                                    therm_cond_sl            !< soil thermal conductivity [W m-2 K-1]
    TYPE(t_jsb_var_real3d)          :: & 
                                    sand_sl             , &  !< soil sand proportion
                                    silt_sl             , &  !< soil silt proportion
                                    clay_sl             , &  !< soil clay proportion (re-calculated as: clay = 1.0 -sand -silt)
                                    qmax_org_min_sl     , &  !< maximum OM sorption capacity of mineral soil, [mol C kg-1 mineral soil]
                                    qmax_po4_min_sl     , &  !< maximum po4 sorption capacity of SOM, [mol P kg-1 SOM]
                                    qmax_po4_om_sl           !< maximum po4 sorption capacity of mineral soil, [mol P kg-1 mineral soil]

                                    


  CONTAINS
    PROCEDURE :: Init => Init_sse_memory
  END TYPE t_sse_memory


  CHARACTER(len=*), PARAMETER :: modname = 'mo_sse_memory_class'


CONTAINS

  !-----------------------------------------------------------------------------------------------------
  !> initialize memory for the SSE_ process
  !!
  !!  
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE Init_sse_memory(mem, prefix, suffix, lct_ids, model_id)

!    USE mo_jsb_varlist,       ONLY: t_jsb_varlist, jsb_add_var
    USE mo_slm_varlist,       ONLY: jsb_add_var   ! replace by mo_jsb_varlist later    *slmtemp
    USE mo_jsb_io,            ONLY: grib_bits, t_cf, t_grib1, t_grib2   ! , TSTEP_CONSTANT, tables
    USE mo_jsb_grid_class,    ONLY: t_jsb_grid, t_jsb_vgrid
    USE mo_jsb_grid,          ONLY: Get_grid, Get_vgrid
    USE mo_jsb_model_class,   ONLY: t_jsb_model

    IMPLICIT NONE
    ! ---------------------------
    ! 0.1.InOut
    CLASS(t_sse_memory), INTENT(inout), TARGET  :: mem             !> sse memory
    CHARACTER(len=*),     INTENT(in)            :: prefix          !> process name
    CHARACTER(len=*),     INTENT(in)            :: suffix          !> tile name
    INTEGER,              INTENT(in)            :: lct_ids !(:)    !> in jsbach4: vector of land cover type IDs
    INTEGER,              INTENT(in)            :: model_id        !> model ID model\%id
    ! ---------------------------
    ! 0.2 Local
    TYPE(t_jsb_model), POINTER :: model
    TYPE(t_jsb_grid),  POINTER :: hgrid                        ! Horizontal grid
    TYPE(t_jsb_vgrid), POINTER :: vgrid_surface                ! Vertical grid
    TYPE(t_jsb_vgrid), POINTER :: vgrid_soil_sb                ! Vertical grid
    INTEGER                    :: table
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':Init_sse_memory'
    ! ---------------------------
    ! 0.5 Get Memory
    model         => Get_model(model_id)
    table         = 1     ! in jsbach4 this line is: table = tables(1)     *slmtemp
    hgrid         => Get_grid(model%grid_id)
    vgrid_surface => Get_vgrid('surface')
    vgrid_soil_sb => Get_vgrid('soil_layer_sb')



    !> 1.0  init memory
    !! 
    IF (.NOT. ASSOCIATED(mem%soil_depth%ptr)) THEN

      !> add variables to sse memory sse\%VAR
      !! 
      !! sse\%VAR
      !!  
      CALL jsb_add_var(mem, 'soil_depth', mem%soil_depth, &
            hgrid, vgrid_surface, &
            t_cf('soil_depth', 'm', 'Soil depth derived from textures (bedrock), for hydrology'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue

      CALL jsb_add_var(mem, 'root_depth', mem%root_depth, &
            hgrid, vgrid_surface, &
            t_cf('root_depth', 'm', 'actual rooting depth'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue

      CALL jsb_add_var(mem, 't_srf_old', mem%t_srf_old, &
            hgrid, vgrid_surface, &
            t_cf('t_srf_old', 'K', 'previous surface temperature'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'zril_old', mem%zril_old, &
            hgrid, vgrid_surface, &
            t_cf('zril_old', '', 'previous timesteps Reynolds number'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue

      CALL jsb_add_var(mem, 'soil_depth_sl', mem%soil_depth_sl, &
            hgrid, vgrid_soil_sb, &
            t_cf('soil_depth_sl', 'm', 'Depth of each soil layer that can be saturated with water, until bedrock'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'interception', mem%interception, &
            hgrid, vgrid_surface, &
            t_cf('interception', 'kg m-2 s-1', 'surface interception evaporation'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'pet', mem%pet, &
            hgrid, vgrid_surface, &
            t_cf('pet', 'kg m-2 s-1', 'potential evaporation'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'evaporation', mem%evaporation, &
            hgrid, vgrid_surface, &
            t_cf('evaporation', 'kg m-2 s-1', 'surface evaporation'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'runoff', mem%runoff, &
            hgrid, vgrid_surface, &
            t_cf('runoff', 'kg m-2 s-1', 'surface runoff'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'drainage', mem%drainage, &
            hgrid, vgrid_surface, &
            t_cf('drainage', 'kg m-2 s-1', 'drainage'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'drainage_fraction', mem%drainage_fraction, &
            hgrid, vgrid_surface, &
            t_cf('drainage_fraction', '1 s-1', 'water mass fraction lost from the soil column by drainage'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'w_skin', mem%w_skin, &
            hgrid, vgrid_surface, &
            t_cf('w_skin', 'm', 'Water content in skin reservoir'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'w_soil_root', mem%w_soil_root, &
            hgrid, vgrid_surface, &
            t_cf('w_soil_root', 'm', 'Water content in root zone of the soil'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'w_soil_root_fc', mem%w_soil_root_fc, &
            hgrid, vgrid_surface, &
            t_cf('w_soil_root_fc', 'm', 'Water content at field capacity in root zone of the soil'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'w_soil_root_pwp', mem%w_soil_root_pwp, &
            hgrid, vgrid_surface, &
            t_cf('w_soil_root_pwp', 'm', 'Water content at permanent wilting point in root zone of the soil'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 't_soil_root', mem%t_soil_root, &
            hgrid, vgrid_surface, &
            t_cf('t_soil_root', 'K', 'Temperature of fine roots in the soil'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue

      CALL jsb_add_var(mem, 'percolation_sl', mem%percolation_sl, &
            hgrid, vgrid_soil_sb, &
            t_cf('percolation_sl', 'fraction s-1', 'fraction of mass transferred to below layer'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue

      CALL jsb_add_var(mem, 'drainage_sl', mem%drainage_sl, &
            hgrid, vgrid_soil_sb, &
            t_cf('drainage_sl', 'kg m-2 s-1', 'drainage of soil layers'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'w_soil_sl', mem%w_soil_sl, &
            hgrid, vgrid_soil_sb, &
            t_cf('w_soil_sl', 'm', 'Water content in soil layers'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'w_soil_sat_sl', mem%w_soil_sat_sl, &
            hgrid, vgrid_soil_sb, &
            t_cf('w_soil_sat_sl', 'm', 'Water content of soil layers at saturation'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'w_soil_fc_sl', mem%w_soil_fc_sl, &
            hgrid, vgrid_soil_sb, &
            t_cf('w_soil_fc_sl', 'm', 'Water content of soil layers at field capacity'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'w_soil_pwp_sl', mem%w_soil_pwp_sl, &
            hgrid, vgrid_soil_sb, &
            t_cf('w_soil_pwp_sl', 'm', 'Water content of soil layers at permanent wilting point'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue

      CALL jsb_add_var(mem, 'saxtonA', mem%saxtonA, &
            hgrid, vgrid_soil_sb, &
            t_cf('saxtonA', '', 'coefficient in moisture-tension relationship'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue

      CALL jsb_add_var(mem, 'saxtonB', mem%saxtonB, &
            hgrid, vgrid_soil_sb, &
            t_cf('saxtonB', '', 'coefficient in moisture-tension relationship'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue

      CALL jsb_add_var(mem, 'saxtonC', mem%saxtonC, &
            hgrid, vgrid_soil_sb, &
            t_cf('saxtonC', '', 'exponent of moisture-conductivity relationship'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue

      CALL jsb_add_var(mem, 'kdiff_sat_sl', mem%kdiff_sat_sl, &
            hgrid, vgrid_soil_sb, &
            t_cf('kdiff_sat_sl', 'm s-1', 'saturated hydraulic conductivity'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue

      CALL jsb_add_var(mem, 'transpiration', mem%transpiration, &
            hgrid, vgrid_surface, &
            t_cf('transpiration', 'kg m-2 s-1', 'transpiration in'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue

      CALL jsb_add_var(mem, 'w_soil_root_theta', mem%w_soil_root_theta, &
            hgrid, vgrid_surface, &
            t_cf('w_soil_root_theta', '', 'root zone integrated plant available water, fraction of maximum '), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=1.0_wp) ! initvalue

      CALL jsb_add_var(mem, 'w_soil_root_pot', mem%w_soil_root_pot, &
            hgrid, vgrid_surface, &
            t_cf('w_soil_root_pot', 'MPa', 'soil water potential'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue

      CALL jsb_add_var(mem, 'w_soil_pot_sl', mem%w_soil_pot_sl, &
            hgrid, vgrid_soil_sb, &
            t_cf('w_soil_pot_sl', 'MPa', 'soil water potential per layer'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue

      CALL jsb_add_var(mem, 't_soil_sl', mem%t_soil_sl, &
            hgrid, vgrid_soil_sb, &
            t_cf('t_soil_sl', 'K', 'soil temperature'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=273.15_wp) ! initvalue

      CALL jsb_add_var(mem, 'heat_capa_sl', mem%heat_capa_sl, &
            hgrid, vgrid_soil_sb, &
            t_cf('heat_capa_sl', 'J kg-1', 'soil heat capacity'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue

      CALL jsb_add_var(mem, 'bulk_dens_sl', mem%bulk_dens_sl, &
            hgrid, vgrid_soil_sb, &
            t_cf('bulk_dens_sl', 'kg m-3', 'soil bulk density'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue

      CALL jsb_add_var(mem, 'therm_cond_sl', mem%therm_cond_sl, &
            hgrid, vgrid_soil_sb, &
            t_cf('therm_cond_sl', 'W m-2 K-1', 'soil thermal conductivity'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue

      CALL jsb_add_var(mem, 'sand_sl', mem%sand_sl, &
            hgrid, vgrid_soil_sb, &
            t_cf('sand_sl', '', 'soil sand proportion'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
            
      CALL jsb_add_var(mem, 'silt_sl', mem%silt_sl, &
            hgrid, vgrid_soil_sb, &
            t_cf('silt_sl', '', 'soil silt proportion'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
            
      CALL jsb_add_var(mem, 'clay_sl', mem%clay_sl, &
            hgrid, vgrid_soil_sb, &
            t_cf('clay_sl', '', 'soil clay proportion'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue

      CALL jsb_add_var(mem, 'qmax_org_min_sl', mem%qmax_org_min_sl, &
            hgrid, vgrid_soil_sb, &
            t_cf('qmax_org_min_sl', 'mol C kg-1 mineral soil', 'maximum OM sorption capacity of mineral soil'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue

      CALL jsb_add_var(mem, 'qmax_po4_min_sl', mem%qmax_po4_min_sl, &
            hgrid, vgrid_soil_sb, &
            t_cf('qmax_po4_min_sl', 'mol P kg-1 SOM', 'maximum po4 sorption capacity of SOM'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue

      CALL jsb_add_var(mem, 'qmax_po4_om_sl', mem%qmax_po4_om_sl, &
            hgrid, vgrid_soil_sb, &
            t_cf('qmax_po4_om_sl', 'mol P kg-1 mineral soil', 'maximum po4 sorption capacity of mineral soil'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
    END IF

  END SUBROUTINE Init_sse_memory

#endif
END MODULE mo_sse_memory_class
