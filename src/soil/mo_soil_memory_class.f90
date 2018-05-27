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
  USE mo_jsb_var_class,          ONLY: t_jsb_var_real2d, t_jsb_var_real3d
  USE mo_slm_varlist,            ONLY: jsb_add_var   ! replace by mo_jsb_varlist later    *slmtemp


  IMPLICIT NONE
  PRIVATE
  PUBLIC :: t_soil_memory, max_no_of_vars

  INTEGER, PARAMETER :: max_no_of_vars      = 50 


  !-----------------------------------------------------------------------------------------------------
  !> Type definition for soil memory
  !!
  !! @par includes: \n
  !!    soil variables soil\%VAR
  !! 
  !-----------------------------------------------------------------------------------------------------
  TYPE, EXTENDS(t_jsb_memory) :: t_soil_memory


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
                                    interception,     & !< surface interception evaporation  [kg / m2 / s]
                                    evaporation,      & !< surface evaporation  [kg / m2 / s]
                                    runoff,           & !< surface runoff  [kg / m2 / s]
                                    drainage,         & !< drainage        [kg / m2 / s]
                                    drainage_fraction   !< water mass fraction lost from the soil column by drainage [1 / s] 

    ! Additional variables for soil
    TYPE(t_jsb_var_real2d)          :: &
                                    w_skin,           & !< Water content in skin reservoir                                          [m]
                                    w_soil_root,      & !< Water content in root zone of the soil                                   [m]
                                    w_soil_root_fc,   & !< Water content at field capacity in root zone of the soil                 [m]
                                    w_soil_root_pwp,  & !< Water content at permanent wilting point in root zone of the soil        [m]
                                    t_soil_root         !< Temperature of fine roots in the soil                                    [K]

    TYPE(t_jsb_var_real3d)          :: &
                                    percolation_sl,   & !< fraction of mass transferred to below layer                              [fraction/s]
                                    drainage_sl,      & !< drainage of soil layers                                                  [kg / m2 / s]
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
                                    kdiff_sat_sl        !< saturated hydraulic conductivity                                         [m / s]

    ! transpiration
    TYPE(t_jsb_var_real2d)          :: &
                                    transpiration       !< transpiration in kg m-2 s-1

    !Soil Water addons for SLM in HYDRO
    TYPE(t_jsb_var_real2d)          :: &
                                    w_soil_root_theta, & !< root zone integrated plant available water (fraction of maximum) 
                                    w_soil_root_pot      !< soil water potential [MPa]

    TYPE(t_jsb_var_real3d)          :: &
                                    w_soil_pot_sl        !< soil water potential per layer [MPa]

    ! soil energy 
    TYPE(t_jsb_var_real3d)          :: &
                                    t_soil_sl            !< soil temperature [K] 

    ! soil physical properties (not sure where they should belong to, but they are needed in the long-term)
    TYPE(t_jsb_var_real3d)          :: & 
                                    heat_capa_sl,         &  !< soil heat capacity (J / kg )
                                    bulk_dens_sl,         &  !< soil bulk density (kg / m3)
                                    therm_cond_sl            !< soil thermal conductivity (W / m2 / K)
    TYPE(t_jsb_var_real3d)          :: & 
                                    sand_sl             , &  !< soil sand proportion
                                    silt_sl             , &  !< soil silt proportion
                                    clay_sl             , &  !< soil clay proportion - recalculated from: clay = 1.0 -sand -silt
                                    qmax_org_min_sl     , &  !< maximum OM sorption capacity of mineral soil, [mol C kg-1 mineral soil]
                                    qmax_po4_min_sl     , &  !< maximum po4 sorption capacity of SOM, [mol P kg-1 SOM]
                                    qmax_po4_om_sl           !< maximum po4 sorption capacity of mineral soil, [mol P kg-1 mineral soil]

                                    


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
  SUBROUTINE Init_soil_memory(mem, prefix, suffix, lct_ids, model_id)

!    USE mo_jsb_varlist,       ONLY: t_jsb_varlist, jsb_add_var
    USE mo_jsb_io,            ONLY: grib_bits
!    USE mo_jsb_io,            ONLY: TSTEP_CONSTANT, tables
    USE mo_jsb_io,            ONLY: t_cf, t_grib1, t_grib2
    USE mo_jsb_grid_class,    ONLY: t_jsb_grid, t_jsb_vgrid
!    USE mo_jsb_grid,          ONLY: Get_grid, Get_vgrid
    USE mo_jsb_model_class,   ONLY: t_jsb_model

    IMPLICIT NONE
    ! ---------------------------
    ! 0.1.InOut
    CLASS(t_soil_memory), INTENT(inout), TARGET :: mem             !> soil memory
    CHARACTER(len=*),     INTENT(in)            :: prefix          !> ?
    CHARACTER(len=*),     INTENT(in)            :: suffix          !> ?
    INTEGER,              INTENT(in)            :: lct_ids !(:)    !> in jsbach4: vector of land cover type IDs
    INTEGER,              INTENT(in)            :: model_id        !> model ID model\%id
    ! ---------------------------
    ! 0.2 Local
    TYPE(t_jsb_model), POINTER :: model
    TYPE(t_jsb_grid),  POINTER :: hgrid                        ! Horizontal grid
!    TYPE(t_jsb_vgrid), POINTER :: canopy_layer                 ! Vertical grid
    TYPE(t_jsb_vgrid), POINTER :: surface                      ! Vertical grid
    INTEGER                    :: table

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':Init_soil_memory'

    model => Get_model(model_id)

!    table = tables(1)
    table = 1

!    hgrid   => Get_grid(mem%grid_id)
!    surface => Get_vgrid('surface')
!    canopy_layer => Get_vgrid('canopy_layer')


    IF (.NOT. ALLOCATED(mem%soil_depth%ptr)) THEN

      !> add variables to soil memory soil\%VAR
      !! 
      !! soil\%VAR
      !!  
      CALL jsb_add_var(mem, 'soil_depth', mem%soil_depth, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'soil_depth_sl', mem%soil_depth_sl, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'zril_old', mem%zril_old, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 't_srf_old', mem%t_srf_old, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'root_depth', mem%root_depth, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'interception', mem%interception, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'evaporation', mem%evaporation, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'runoff', mem%runoff, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'drainage', mem%drainage, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'drainage_fraction', mem%drainage_fraction, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'drainage_sl', mem%drainage_sl, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'percolation_sl', mem%percolation_sl, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'w_soil_root', mem%w_soil_root, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'w_skin', mem%w_skin, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'w_soil_root_fc', mem%w_soil_root_fc, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'w_soil_root_pwp', mem%w_soil_root_pwp, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 't_soil_root', mem%t_soil_root, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'w_soil_root_pot', mem%w_soil_root_pot, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'w_soil_root_theta', mem%w_soil_root_theta, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=1.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'w_soil_sl', mem%w_soil_sl, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'w_soil_sat_sl', mem%w_soil_sat_sl, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'w_soil_fc_sl', mem%w_soil_fc_sl, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
  
      CALL jsb_add_var(mem, 'w_soil_pwp_sl', mem%w_soil_pwp_sl, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue

      CALL jsb_add_var(mem, 'saxtonA', mem%saxtonA, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue

      CALL jsb_add_var(mem, 'saxtonB', mem%saxtonB, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue

      CALL jsb_add_var(mem, 'saxtonC', mem%saxtonC, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue

      CALL jsb_add_var(mem, 'kdiff_sat_sl', mem%kdiff_sat_sl, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue

      CALL jsb_add_var(mem, 'transpiration', mem%transpiration, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue

      CALL jsb_add_var(mem, 'w_soil_pot_sl', mem%w_soil_pot_sl, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue

      CALL jsb_add_var(mem, 't_soil_sl', mem%t_soil_sl, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=273.15_wp) ! initvalue

      CALL jsb_add_var(mem, 'heat_capa_sl', mem%heat_capa_sl, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue

      CALL jsb_add_var(mem, 'bulk_dens_sl', mem%bulk_dens_sl, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue

      CALL jsb_add_var(mem, 'therm_cond_sl', mem%therm_cond_sl, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue

      CALL jsb_add_var(mem, 'sand_sl', mem%sand_sl, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
            
      CALL jsb_add_var(mem, 'silt_sl', mem%silt_sl, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue
            
      CALL jsb_add_var(mem, 'clay_sl', mem%clay_sl, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue

      CALL jsb_add_var(mem, 'qmax_org_min_sl', mem%qmax_org_min_sl, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue

      CALL jsb_add_var(mem, 'qmax_po4_min_sl', mem%qmax_po4_min_sl, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue

      CALL jsb_add_var(mem, 'qmax_po4_om_sl', mem%qmax_po4_om_sl, &
            hgrid, surface, &
            t_cf('', '', ''), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            initval_r=0.0_wp) ! initvalue

    END IF

  END SUBROUTINE Init_soil_memory

END MODULE mo_soil_memory_class
