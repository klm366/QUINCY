!>
!! @brief Contains methods for the carbon memory.
!!
!! @author
!!  Rainer Schneck, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Rainer Schneck (2015-04-02)
!!
MODULE mo_carbon_memory
#ifndef __NO_JSBACH__

  USE mo_carbon_memory_class, ONLY: t_carbon_memory, t_soil_pools, t_veg_pools, t_photosyn, t_canopy_layer
  USE mo_kind,                ONLY: wp
  USE mo_exception,           ONLY: message

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: new_carbon_memory

  CHARACTER(len=*), PARAMETER :: modname = 'mo_carbon_memory'

CONTAINS

  FUNCTION new_carbon_memory(varlist_name, prefix, suffix, tile) RESULT(mem)

    USE mo_jsb_varlist,       ONLY: t_var_list, jsb_get_var_list, jsb_new_var_list, jsb_add_var
    USE mo_jsb_io,            ONLY: grib_bits, nc_bits, t_cf, t_grib1, t_grib2, &
                                    FILETYPE_NC2, TSTEP_CONSTANT, ZAXIS_SURFACE, ZAXIS_DEPTH_BELOW_LAND, tables
    USE mo_jsb_grid_class,    ONLY: t_jsb_grid, t_jsb_vgrid, new_vgrid
    USE mo_jsb_tile_class,    ONLY: t_jsb_tile
    USE mo_jsb_model_class,   ONLY: t_jsb_model, get_vgrid, register_vgrid
    USE mo_jsb_base,          ONLY: get_model

    CHARACTER(len=*),           INTENT(in) :: varlist_name
    CHARACTER(len=*),           INTENT(in) :: prefix
    CHARACTER(len=*),           INTENT(in) :: suffix
#ifdef HAVE_F2003
    TYPE(t_jsb_tile), POINTER , INTENT(in) :: tile
#else
    TYPE(t_jsb_tile), POINTER              :: tile
#endif
    TYPE(t_carbon_memory), POINTER :: mem
    TYPE(t_carbon_memory), POINTER :: new
    TYPE(t_var_list),   POINTER :: this_list
    TYPE(t_jsb_model), POINTER  :: model

    INTEGER                    :: nsoil_water
    TYPE(t_jsb_grid),  POINTER :: hgrid                 ! Horizontal grid
    TYPE(t_jsb_vgrid), POINTER :: surface               ! Vertical grids
    INTEGER :: table

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':new_carbon_memory'

    model => get_model(tile%owner_model_id)

    ALLOCATE(new)

    table = tables(1)

    !this_list => get_var_list(TRIM(varlist_name))
    CALL jsb_get_var_list(this_list, TRIM(varlist_name))
    IF (.NOT. ASSOCIATED(this_list)) THEN
      !this_list => new_var_list(TRIM(varlist_name), patch_id=patch_id)
      ALLOCATE(this_list)
      CALL jsb_new_var_list(this_list, TRIM(varlist_name),              &
                            lrestart=.TRUE., restart_type=FILETYPE_NC2, &
                            patch_id=model%grid%patch%id, table=table)
    END IF

    hgrid        => model%grid

    surface => get_vgrid(model, 'surface')
    IF (.NOT. ASSOCIATED(surface)) THEN
      surface => new_vgrid('surface', ZAXIS_SURFACE, 1)
      CALL register_vgrid(model, surface)
    END IF
    
    new%t_soil_pools   => new_soil_pools()
    new%t_veg_pools    => new_veg_pools()
    new%t_t_photosyn   => new_t_photosyn()
    new%t_canopy_layer => new_canopy_layer()
    
    mem => new

  CONTAINS

    FUNCTION new_soil_pools() RESULT(soilstate)

      TYPE(t_soil_pools), POINTER :: soilstate

      ALLOCATE(soilstate)

! R: nur zur Übersicht über call Optionen:      
!     TYPE(t_var_list),     INTENT(inout)        :: this_list           ! list
!     CHARACTER(len=*),     INTENT(in)           :: name                ! name of variable
!     REAL(wp),             POINTER              :: ptr(:,:)            ! reference to field
!     TYPE(t_jsb_grid),     POINTER, INTENT(in)  :: hgrid               ! horizontal grid type
!     TYPE(t_jsb_vgrid),    POINTER, INTENT(in)  :: vgrid               ! vertical grid type
!     TYPE(t_cf),           INTENT(in)           :: cf                  ! CF related metadata
!     TYPE(t_grib1),        INTENT(in)           :: grib1               ! GRIB1 related metadata
!     TYPE(t_grib2),        INTENT(in)           :: grib2               ! GRIB2 related metadata
!     LOGICAL,              INTENT(in), OPTIONAL :: loutput             ! output flag
!     LOGICAL,              INTENT(in), OPTIONAL :: lcontainer          ! container flag
!     LOGICAL,              INTENT(in), OPTIONAL :: lrestart            ! restart flag
!     LOGICAL,              INTENT(in), OPTIONAL :: lrestart_cont       ! continue restart if var not available
!     REAL(wp),             INTENT(in), OPTIONAL :: initval_r           ! value if var not available
!     INTEGER,              INTENT(in), OPTIONAL :: isteptype           ! type of statistical processing
!     REAL(wp),             INTENT(in), OPTIONAL :: resetval_r          ! reset value (after accumulation)
!     LOGICAL,              INTENT(in), OPTIONAL :: lmiss               ! missing value flag
!     REAL(wp),             INTENT(in), OPTIONAL :: missval_r           ! missing value
!     INTEGER,              INTENT(in), OPTIONAL :: tlev_source         ! actual TL for TL dependent vars
!     TYPE(t_var_metadata), POINTER,    OPTIONAL :: info                ! returns reference to metadata
!     LOGICAL,              INTENT(in), OPTIONAL :: verbose             ! print information
!     TYPE(t_list_element), POINTER, OPTIONAL  :: new_element           ! pointer to new var list element
!     INTEGER,              INTENT(in), OPTIONAL :: cdiTimeID           ! CDI time mode (TIME_VARIABLE/TIME_CONSTANT)
      

                                                                 
      CALL jsb_add_var( this_list, prefix//'YCpool_acid_ag1'//suffix, soilstate%YCpool_acid_ag1,  & ! R: Alle Variablen REAL(wp), POINTER, DIMENSION(:,:)
        & hgrid, surface,                                                                         &
        & t_cf('Yasso above ground litter-pool for acid soluble litter', '', '', nc_bits),        &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                      & ! R: übergibt den TYPE t_grids und keine Instanz davon.  255 =table number?
        & initval_r=0.0_wp )                                                                        ! R: initval_r  =value if var not available
     

      CALL jsb_add_var( this_list, prefix//'YCpool_water_ag1'//suffix, soilstate%YCpool_water_ag1,  & 
        & hgrid, surface,                                                                         &
        & t_cf('Yasso above ground litter-pool for water soluble litter', '', '', nc_bits),        &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                      & 
        & initval_r=0.0_wp )

      CALL jsb_add_var( this_list, prefix//'YCpool_ethanol_ag1'//suffix, soilstate%YCpool_ethanol_ag1,  & 
        & hgrid, surface,                                                                         &
        & t_cf('Yasso above ground litter-pool for ethanol soluble litter', '', '', nc_bits),        &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                      & 
        & initval_r=0.0_wp )

      CALL jsb_add_var( this_list, prefix//'YCpool_nonsoluble_ag1'//suffix, soilstate%YCpool_nonsoluble_ag1,  & 
        & hgrid, surface,                                                                         &
        & t_cf('Yasso above ground litter-pool for non-soluble litter', '', '', nc_bits),        &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                      & 
        & initval_r=0.0_wp )
 
      CALL jsb_add_var( this_list, prefix//'YCpool_acid_bg1'//suffix, soilstate%YCpool_acid_bg1,  & 
        & hgrid, surface,                                                                         &
        & t_cf('Yasso below ground litter-pool for acid soluble litter', '', '', nc_bits),        &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                      & 
        & initval_r=0.0_wp )

      CALL jsb_add_var( this_list, prefix//'YCpool_water_bg1'//suffix, soilstate%YCpool_water_bg1,  & 
        & hgrid, surface,                                                                         &
        & t_cf('Yasso below ground litter-pool for water soluble litter', '', '', nc_bits),        &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                      & 
        & initval_r=0.0_wp )

      CALL jsb_add_var( this_list, prefix//'YCpool_ethanol_bg1'//suffix, soilstate%YCpool_ethanol_bg1,  & 
        & hgrid, surface,                                                                         &
        & t_cf('Yasso below ground litter-pool for ethanol soluble litter', '', '', nc_bits),        &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                      & 
        & initval_r=0.0_wp )

      CALL jsb_add_var( this_list, prefix//'YCpool_nonsoluble_bg1'//suffix, soilstate%YCpool_nonsoluble_bg1,  & 
        & hgrid, surface,                                                                         &
        & t_cf('Yasso below ground litter-pool for non-soluble litter', '', '', nc_bits),        &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                      & 
        & initval_r=0.0_wp )

      CALL jsb_add_var( this_list, prefix//'YCpool_humus_1'//suffix, soilstate%YCpool_humus_1,  & 
        & hgrid, surface,                                                                         &
        & t_cf('Yasso below ground litter-pool for slow C compartment', '', '', nc_bits),        &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                      & 
        & initval_r=0.0_wp )

      CALL jsb_add_var( this_list, prefix//'YCpool_acid_ag2'//suffix, soilstate%YCpool_acid_ag2,  & 
        & hgrid, surface,                                                                         &
        & t_cf('Yasso above ground litter-pool for acid soluble litter', '', '', nc_bits),        &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                      & 
        & initval_r=0.0_wp )

      CALL jsb_add_var( this_list, prefix//'YCpool_water_ag2'//suffix, soilstate%YCpool_water_ag2,  & 
        & hgrid, surface,                                                                         &
        & t_cf('Yasso above ground litter-pool for water soluble litter', '', '', nc_bits),        &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                      & 
        & initval_r=0.0_wp )

      CALL jsb_add_var( this_list, prefix//'YCpool_ethanol_ag2'//suffix, soilstate%YCpool_ethanol_ag2,  & 
        & hgrid, surface,                                                                         &
        & t_cf('Yasso above ground litter-pool for ethanol soluble litter', '', '', nc_bits),        &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                      & 
        & initval_r=0.0_wp )

      CALL jsb_add_var( this_list, prefix//'YCpool_nonsoluble_ag2'//suffix, soilstate%YCpool_nonsoluble_ag2,  & 
        & hgrid, surface,                                                                         &
        & t_cf('Yasso above ground litter-pool for non-soluble litter', '', '', nc_bits),        &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                      & 
        & initval_r=0.0_wp )

      CALL jsb_add_var( this_list, prefix//'YCpool_acid_bg2'//suffix, soilstate%YCpool_acid_bg2,  & 
        & hgrid, surface,                                                                         &
        & t_cf('Yasso below ground litter-pool for acid soluble litter', '', '', nc_bits),        &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                      & 
        & initval_r=0.0_wp )

      CALL jsb_add_var( this_list, prefix//'YCpool_water_bg2'//suffix, soilstate%YCpool_water_bg2,  & 
        & hgrid, surface,                                                                         &
        & t_cf('Yasso below ground litter-pool for water soluble litter', '', '', nc_bits),        &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                      & 
        & initval_r=0.0_wp )

      CALL jsb_add_var( this_list, prefix//'YCpool_ethanol_bg2'//suffix, soilstate%YCpool_ethanol_bg2,  & 
        & hgrid, surface,                                                                         &
        & t_cf('Yasso below ground litter-pool for ethanol soluble litter', '', '', nc_bits),        &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                      & 
        & initval_r=0.0_wp )

      CALL jsb_add_var( this_list, prefix//'YCpool_nonsoluble_bg2'//suffix, soilstate%YCpool_nonsoluble_bg2,  & 
        & hgrid, surface,                                                                         &
        & t_cf('Yasso below ground litter-pool for non-soluble litter', '', '', nc_bits),        &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                      & 
        & initval_r=0.0_wp )

      CALL jsb_add_var( this_list, prefix//'YCpool_humus_2'//suffix, soilstate%YCpool_humus_2,  & 
        & hgrid, surface,                                                                         &
        & t_cf('Yasso below ground litter-pool for slow C compartment', '', '', nc_bits),        &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                      & 
        & initval_r=0.0_wp )

    END FUNCTION new_soil_pools

    FUNCTION new_water_state() RESULT(wstate)

      TYPE(t_carbon_water_state), POINTER :: wstate

      ALLOCATE(wstate)

      CALL jsb_add_var( this_list, prefix//'wfract_can'//suffix, wstate%wfract_can,                         &
        & hgrid, surface,                                                           &
        & t_cf('wet_fraction_canopy', '', '', nc_bits),                                      &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                       &
        & initval_r=0.0_wp )

      CALL jsb_add_var( this_list, prefix//'wsr_can'//suffix, wstate%wsr_can,                               &
        & hgrid, surface,                                                           &
        & t_cf('skin_reservoir_canopy', 'm', '', nc_bits),                                   &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                       &
        & initval_r=0.0_wp )

      CALL jsb_add_var( this_list, prefix//'sfract_can'//suffix, wstate%sfract_can,                         &
        & hgrid, surface,                                                           &
        & t_cf('snow_fraction_canopy', '', '', nc_bits),                                     &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                       &
        & initval_r=0.0_wp )

      CALL jsb_add_var( this_list, prefix//'wsn_can'//suffix, wstate%wsn_can,                               &
        & hgrid, surface,                                                           &
        & t_cf('snow_water_equivalent_canopy', 'm', '', nc_bits),                            &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                       &
        & initval_r=0.0_wp )

      CALL jsb_add_var( this_list, prefix//'ws_root'//suffix, wstate%ws_root,                               &
        & hgrid, surface,                                                           &
        & t_cf('soil_moisture_root', 'm', '', nc_bits),                                      &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                       &
        & initval_r=0.0_wp )

      CALL jsb_add_var( this_list, prefix//'ws_fc_root'//suffix, wstate%ws_fc_root,                         &
        & hgrid, surface,                                                           &
        & t_cf('field_capacity_root', 'm', '', nc_bits),                                     &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                       &
        & initval_r=0.0_wp )

      CALL jsb_add_var( this_list, prefix//'ws_pwp_root'//suffix, wstate%ws_pwp_root,                       &
        & hgrid, surface,                                                           &
        & t_cf('wilting_root', 'm', '', nc_bits),                                            &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                       &
        & initval_r=0.0_wp )

      CALL jsb_add_var( this_list, prefix//'water_stress'//suffix, wstate%water_stress,                     &
        & hgrid, surface,                                                           &
        & t_cf('water_stress_factor', '', '', nc_bits),                                      &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                       &
        & initval_r=0.0_wp )

    END FUNCTION new_water_state

  END FUNCTION new_carbon_memory

#endif
END MODULE mo_carbon_memory
