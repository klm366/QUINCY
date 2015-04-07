!>
!! @brief Initialization of the the carbon memory.
!!
!! @author
!!  Rainer Schneck, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Rainer Schneck (2015-04-29)
!!


Die SUBROUTINE jsbach_init wird nur von jsbalone_driver.f90:127: aufgerufen und nur diese Subroutine ruft   
      IF (theOptions%Usecarbon) THEN
       CALL init_cbalance_carbon(theLand%Grid, theLand%Domain,theLand%Cbalance, theLand%Nbalance, &
                                theOptions%WithNitrogen,theOptions%UseDynveg, theOptions%lcc_scheme, &
                                theOptions%UseLanduseTransitions, &
                                theOptions%FileType, theOptions%FileZtype, theLand%IO_veg, theLand%IO_nitro, &
                                theLand%IO_yasso, theOptions%WithYasso)
       CALL init_carbon(theLand%Grid, theLand%Domain, theLand%carbon, theOptions%FileType, theOptions%FileZtype, &
                       theLand%IO_diag_stream, theLand%IO_stream)
    ENDIF
diese beiden Routinen auf.
=> da driver erst später eingebaut werden sind diese Routinen unnötig.


MODULE mo_carbon_init
#ifndef __NO_JSBACH__

  USE mo_jsb_model_class,     ONLY: t_jsb_model
  USE mo_jsb_tile_class,      ONLY: t_jsb_tile
  USE mo_carbon_memory_class, ONLY: t_carbon_memory, t_soil_pools, t_veg_pools, t_photosyn, t_canopy_layer
  USE mo_jsb_base,            ONLY: get_model
  USE mo_jsb_io_netcdf,       ONLY: t_input_file, jsb_netcdf_close_file, jsb_netcdf_open_input, &
                                    jsb_netcdf_read_2d_time, jsb_netcdf_read_2d
  USE mo_kind,                ONLY: wp
  USE mo_exception,           ONLY: message

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: carbon_init_bc, carbon_init_ic

  CHARACTER(len=*), PARAMETER :: modname = 'mo_carbon_init'

CONTAINS



   SUBROUTINE carbon_init_bc(tile)
!vergleiche SUBROUTINE init_bethy
    TYPE(t_jsb_tile), INTENT(inout) :: tile

    TYPE(t_jsb_model),     POINTER  :: model
    TYPE(t_carbon_memory), POINTER  :: carbon

    REAL(wp), POINTER :: lai_mon(:,:,:)

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':carbon_init_bc'

    model => get_model(tile%owner_model_id)

    ! Get carbon memory of the tile
    carbon => tile%mem%carbon

    ! ...

    ! Propagate boundary conditions down the tile structure
    ! TBD
    
    
    
    carbon%e_transport_rate = 0._wp
    carbon%canopy_conductance = 0._wp
    carbon%canopy_conductance_unlimited = 0._wp
    carbon%CO2_conc_leaf_acc = 0._wp
    carbon%CO2_conc_leaf_unlimited_acc = 0._wp
    carbon%net_assimilation = 0._wp
    carbon%apar_acc = 0._wp
    carbon%par_acc = 0._wp
    carbon%gross_assimilation_acc = 0._wp
    carbon%dark_respiration_acc = 0._wp
    carbon%CO2_flux_net_acc = 0._wp
    carbon%CO2_flux_herbivory_acc = 0._wp
    carbon%CO2_emission_landcover_change_acc = 0._wp
    carbon%CO2_emission_harvest_acc = 0._wp
    carbon%CO2_flux_dynveg_acc = 0._wp

    carbon_diag%net_assimilation = 0._wp
    carbon_diag%faPAR = 0._wp
    carbon_diag%apar_acc = 0._wp
    carbon_diag%par_acc = 0._wp
    carbon_diag%gross_assimilation = 0._wp
    carbon_diag%dark_respiration = 0._wp


  END SUBROUTINE carbon_init_bc

   RECURSIVE SUBROUTINE carbon_init_ic(tile)
! 
!     USE mo_veg_process, ONLY: foliage_projected_cover
!     USE mo_veg_util,    ONLY: water_in_root_zone
!     USE mo_soil_util,   ONLY: soil_depth_to_layers_2d
     USE mo_jsb_time,    ONLY: get_time_interpolation_weights

    TYPE(t_jsb_tile),  INTENT(inout) :: tile

    TYPE(t_jsb_model),        POINTER :: model
    TYPE(t_carbon_memory),    POINTER :: carbon
    TYPE(t_soil_pools),       POINTER :: soil_pools
    TYPE(t_veg_pools),        POINTER :: veg_pools
    INTEGER :: i
! 
!     REAL(wp), POINTER ::        &
!       & lai_mon(:,:,:),         &
!       & return_pointer (:,:)       !< temporary pointer
! 
     TYPE(t_input_file) :: input_file
! 
!     REAL(wp) :: wgt1, wgt2
!     INTEGER  :: nmw1,nmw2
! 
     CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':carbon_init_ic'
 
     model => get_model(tile%owner_model_id)
 
!     ! Initialize tile fractions
!     ! TBD
! 
!     ! Get veg memory of the tile
     carbon => tile%mem%carbon
! 
!     IF (.NOT. ASSOCIATED(veg)) RETURN
! 
!     CALL message(TRIM(routine), 'Initializing veg memory for tile '//TRIM(tile%name))
! 
!     ! Initialize physical state
!     !
!     phys => veg%phys_state
! 
!     input_file = jsb_netcdf_open_input(TRIM(model%config%veg%bc_filename), model%id)
! 
!     ! Leaf area index
!     ! Get monthly LAI climatology
!     lai_mon => jsb_netcdf_read_2d_time(    &
!       & input_file,                        &
!       & variable_name='lai_clim')
!     lai_mon(:,:,:) = MERGE(lai_mon(:,:,:), 0._wp, lai_mon(:,:,:) > 0._wp)
!     phys%lai_mon(:,:,2:13) = lai_mon(:,:,1:12)
!     phys%lai_mon(:,:,1   ) = lai_mon(:,:,12  )
!     phys%lai_mon(:,:,14  ) = lai_mon(:,:,1   )
!     DEALLOCATE(lai_mon)
! 
!     CALL get_time_interpolation_weights(wgt1, wgt2, nmw1, nmw2)
!     phys%lai(:,:) = wgt1 * phys%lai_mon(:,:,nmw1+1) + wgt2 * phys%lai_mon(:,:,nmw2+1)
! 
!     ! Forest fraction
!     return_pointer => jsb_netcdf_read_2d(        &
!       & input_file,                              &
!       & variable_name='forest_fract',            &
!       & fill_array = phys%fract_forest)
!     phys%fract_forest(:,:) = MERGE(phys%fract_forest(:,:), 0._wp, phys%fract_forest(:,:) >= 0._wp)
! 
!     phys%canopy_cond(:,:) = 0._wp
! 
!     phys%skin_can_depth(:,:) = model%config%veg%wsr_capacity * phys%lai(:,:)
! 
!     return_pointer => jsb_netcdf_read_2d(         &
!       & input_file,                               &
!       & variable_name='albedo_veg_vis',           &
!       & fill_array = phys%alb_vis_can)
!     phys%alb_vis_can(:,:) = MERGE(phys%alb_vis_can(:,:), 0._wp, phys%alb_vis_can(:,:) >= 0._wp)
! 
!     return_pointer => jsb_netcdf_read_2d(         &
!       & input_file,                               &
!       & variable_name='albedo_veg_nir',           &
!       & fill_array = phys%alb_nir_can)
!     phys%alb_nir_can(:,:) = MERGE(phys%alb_nir_can(:,:), 0._wp, phys%alb_nir_can(:,:) >= 0._wp)
! 
!     CALL jsb_netcdf_close_file(input_file)
!     input_file = jsb_netcdf_open_input(TRIM(model%config%soil%bc_filename), model%id)
! 
!     return_pointer => jsb_netcdf_read_2d(         &
!       & input_file,                               &
!       & variable_name='root_depth',               &
!       & fill_array = phys%root_depth)
!     phys%root_depth(:,:) = MERGE(phys%root_depth(:,:), 0._wp, phys%root_depth(:,:) >= 0._wp)
!     phys%root_depth(:,:) = MIN(phys%root_depth(:,:), tile%mem%soil%phys_state%soil_depth(:,:))  ! Temporary fix (root depth sometimes larger then soil depth)
! !!$    phys%root_depth(:,:) = tile%mem%soil%phys_state%soil_depth(:,:) * 0.75_wp  ! Temporary
! 
!     phys%root_depth_l(:,:,:) = soil_depth_to_layers_2d(               &
!       & phys%root_depth(:,:),                                         & ! Total rooting depth
!       & model%config%soil%dz_water(1:model%config%soil%nsoil_water))    ! Soil layer thicknesses from namelist
! 
!     CALL jsb_netcdf_close_file(input_file)                              ! Temporary
!     input_file = jsb_netcdf_open_input('bc_land_frac.nc', model%id)     ! Temporary
! 
!     ! Maximum foliage projected cover
!     return_pointer => jsb_netcdf_read_2d(        &
!       & input_file,                              &
!       & variable_name='veg_ratio_max',           &
!       & fill_array = phys%fract_fpc_max)
!     phys%fract_fpc_max(:,:) = MERGE(phys%fract_fpc_max(:,:), 0._wp, phys%fract_fpc_max(:,:) >= 0._wp)
! 
! !!$    phys%fract_fpc(:,:) = foliage_projected_cover(phys%fract_fpc_max(:,:), phys%lai(:,:), phys%clump(:,:))
!     phys%fract_fpc(:,:) = foliage_projected_cover(phys%fract_fpc_max(:,:), phys%lai(:,:), &
!                           SPREAD(SPREAD(2._wp, DIM=1, NCOPIES=SIZE(phys%lai,1)), DIM=2, NCOPIES=SIZE(phys%lai,2))) ! Temporary
! 
!     CALL jsb_netcdf_close_file(input_file)
! 
!     ! Initialize water state
!     !
!     water => veg%water_state
! 
!     water%wfract_can  (:,:) = 0._wp
!     water%wsr_can     (:,:) = 0._wp
!     water%sfract_can  (:,:) = 0._wp
!     water%wsn_can     (:,:) = 0._wp
! 
!     ! Calculate water content and field capacity in the root zone
!     soil_water => tile%mem%soil%water_state
!     soil_phys => tile%mem%soil%phys_state
!     DO i=1,SIZE(water%ws_root,2)
!       water%ws_root     (:,i) = water_in_root_zone( soil_water%ws_l     (:,:,i), &
!                                                   & soil_phys%soil_depth_l(:,:,i), phys%root_depth_l(:,:,i))
!       water%ws_fc_root  (:,i) = water_in_root_zone( soil_water%ws_fc_l  (:,:,i), &
!                                                   & soil_phys%soil_depth_l(:,:,i), phys%root_depth_l(:,:,i))
!       water%ws_pwp_root (:,i) = water_in_root_zone( soil_water%ws_pwp_l (:,:,i), &
!                                                   & soil_phys%soil_depth_l(:,:,i), phys%root_depth_l(:,:,i))
!     END DO
! !!$    water%ws_root    (:,:) = 0._wp
! !!$    water%ws_fc_root (:,:) = 0._wp
! !!$    water%ws_pwp_root(:,:) = 0._wp
! 
!     water%water_stress(:,:) = 0._wp
! 
!     ! Propagate initial conditions down the tile structure, starting with the tile with which this subroutine was called initially
!     ! TBD: sub-tiles should not read the netcdf files again
!     IF (tile%no_of_sub_tiles > 0) THEN
!       DO i=1,tile%no_of_sub_tiles
!         CALL veg_init_ic(tile%sub_tiles(i))
!       END DO
!     END IF
! 
   END SUBROUTINE carbon_init_ic
! 
! #endif
 END MODULE mo_carbon_init
