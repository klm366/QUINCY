!>
!!  @par Copyright
!!    QUINCY terrestrial biosphere model with coupled carbon, nitrogen, and phosphorus cycles
!!    Copyright (C) 2019 Max Planck Institute for Biogeochemistry (MPI-BGC), Jena, Germany 
!!    Contributors: Caldararu S., Engel J., Kern M., Pallandt M., Schnur R., Thum T., Yu L., Zaehle S.
!!
!!    This file is part of the QUINCY model.
!!    The QUINCY model is free software: you can redistribute it and/or modify
!!    it under the terms of the GNU General Public License version 3,
!!    as published by the Free Software Foundation.
!!
!!    The QUINCY model is distributed in the hope that it will be useful,
!!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!!    GNU General Public License for more details.
!!
!!    You should have received a copy of the GNU General Public License
!!    along with the QUINCY model.  If not, see <https://www.gnu.org/licenses/>.
!! 
!! @brief Contains structures and methods for config of the assimilation process
!!
!! @author
!!  Reiner Schnur, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Reiner Schnur (2013-04-04);
!! integration into SLM by Jan Engel, MPI BGC (20160722)
!!
MODULE mo_assimi_config_class
#ifndef __NO_JSBACH__
  ! -------------------------------------------------------------------------------------------------------
  ! Used variables of module

  USE mo_exception,         ONLY: message_text, message, finish
  USE mo_io_units,          ONLY: filename_max
  USE mo_kind,              ONLY: wp
  USE mo_jsb_math_constants,ONLY: eps4
  USE mo_jsb_control,       ONLY: debug_on
  USE mo_jsb_config_class,  ONLY: t_jsb_config

  ! -------------------------------------------------------------------------------------------------------
  ! Module variables
  IMPLICIT NONE
  
  PRIVATE
  PUBLIC :: t_assimi_config

  !-----------------------------------------------------------------------------------------------------
  !> configuration of the assimilation process, derived from t_jsb_config
  !! 
  !! currently it does mainly: reading parameters from namelist
  !-----------------------------------------------------------------------------------------------------
  TYPE, EXTENDS(t_jsb_config) :: t_assimi_config
    LOGICAL                 :: flag_optimal_Nfraction     !< optimise leaf internal N allocation
    LOGICAL                 :: flag_t_resp_acclimation    !< whether or not to use the respiration temperature acclimation factor
    LOGICAL                 :: flag_t_jmax_acclimation    !< acclimation of optimum temperature for Jmax 
    CHARACTER(15)           :: canopy_layer_scheme        !< select scheme to calculate canopy layer thicknesses [standard|fapar]
    CHARACTER(15)           :: canopy_conductance_scheme  !< select scheme to calculate canopy conductance [medlyn|ballberry]
  CONTAINS
    PROCEDURE :: Init => Init_assimi_config
  END type t_assimi_config


  CHARACTER(len=*), PARAMETER :: modname = 'mo_assimi_config_class'

CONTAINS

  !-----------------------------------------------------------------------------------------------------
  !> configuration routine of t_assimi_config
  !! 
  !! currently it does only: read parameters from namelist
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE Init_assimi_config(config)   ! 'config' == t_assimi_config, see ABSTRACT INTERFACE -> SUBROUTINE Init_config(config)

!    USE mo_jsb_namelist_iface, ONLY: open_nml, POSITIONED, position_nml, close_nml
!    using the old JSBACH3 functions instead of wrappers implemented in JSBACH4 -> mo_jsb_namelist_iface
    USE mo_namelist,           ONLY: open_nml, POSITIONED, position_nml, close_nml  ! read model options from file
    USE mo_jsb_grid_class,     ONLY: t_jsb_vgrid, new_vgrid
    USE mo_jsb_grid,           ONLY: Register_vgrid
    USE mo_jsb_grid_class,     ONLY: ZAXIS_GENERIC   ! in JSBACH4 used from mo_jsb_io    *slmtemp

    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    CLASS(t_assimi_config), INTENT(inout) :: config    !< config type for assimi
    ! ---------------------------
    ! 0.2 Local
    ! variables for reading from namlist, identical to variable-name in namelist
    LOGICAL                     :: active
    CHARACTER(len=filename_max) :: ic_filename, bc_filename
    REAL(wp), ALLOCATABLE       :: canopy_layer_thickness_profile(:)  !< use LAI per layer [m2 m-2], indirectly defining the thickness (i.e., vertical height) of each single canopy layer 
    LOGICAL                     :: flag_optimal_Nfraction     , &
                                   flag_t_resp_acclimation    , &
                                   flag_t_jmax_acclimation
    CHARACTER(15)               :: canopy_layer_scheme
    CHARACTER(15)               :: canopy_conductance_scheme
    INTEGER                     :: ncanopy    , &
                                   icanopy
    REAL(wp)                    :: dAPAR
    REAL(wp)                    :: k_canopy_layer                !< extinction coefficient to calculate canopy layer depth
    TYPE(t_jsb_vgrid), POINTER  :: vgrid_canopy
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':Init_assimi_config'

    NAMELIST /assimilation_ctl/    &
        active,                    &
        ic_filename,               &
        bc_filename,               &
        flag_optimal_Nfraction,    &
        flag_t_resp_acclimation,   &
        flag_t_jmax_acclimation,   &
        ncanopy,                   &
        canopy_layer_scheme,       &
        canopy_conductance_scheme

    ! variables for reading model-options from namelist
    INTEGER  :: nml_handler, nml_unit, istat

    ! extinction coefficient for canopy layer depth calculation
    k_canopy_layer = 0.7_wp

    IF (debug_on()) CALL message(TRIM(routine), 'Starting assimi configuration')

    ! Set defaults
    active                         = .TRUE.
    bc_filename                    = 'bc_land_assimi.nc'
    ic_filename                    = 'ic_land_assimi.nc'
    flag_optimal_Nfraction         = .FALSE.      ! .TRUE. .FALSE.
    flag_t_resp_acclimation        = .FALSE.      ! .TRUE. .FALSE.
    flag_t_jmax_acclimation        = .FALSE.      ! .TRUE. .FALSE.
    ncanopy                        = 20           ! (default for FAPAR based canopy layers would be 10)
    canopy_layer_scheme            = "standard"   ! 
    canopy_conductance_scheme      = "medlyn"

    ! read the namelist
    ! &ASSIMILATION_CTL
    nml_handler = open_nml(TRIM(config%namelist_filename))

    nml_unit = position_nml('ASSIMILATION_CTL', nml_handler, STATUS=istat)
    IF (istat == POSITIONED) READ(nml_unit, assimilation_ctl)

    CALL close_nml(nml_handler)

    ! pass values as read from file (ncanopy is read from namelist and saved in the canopy_layer vgrid)
    config%active                        = active
    config%ic_filename                   = ic_filename
    config%bc_filename                   = bc_filename
    config%flag_optimal_Nfraction        = flag_optimal_Nfraction
    config%flag_t_resp_acclimation       = flag_t_resp_acclimation
    config%flag_t_jmax_acclimation       = flag_t_jmax_acclimation
    config%canopy_layer_scheme           = TRIM(canopy_layer_scheme)
    config%canopy_conductance_scheme     = TRIM(canopy_conductance_scheme)


    !< canopy layer depth calculation
    !! 
    !! here the thickness of each canopy layer is calculated (saved in the vgrid as dz(:)) \n
    !! this thickness is not directly given (by e.g. using the unit meter), but represented by the leaf mass in this layer \n
    !! this leaf mass is measured as LAI [m2 m-2], that is, the actual thickness (in m) of canopy layer 2 depends on how much 
    !!   vertical space (hight) needs to be covered to get the LAI specific to canopy layer 2
    !! 
    !! layer 1 is the top canopy layer (reversed order of vgrid axis) \n
    !! and the value of the upper bound of a layer is larger than the value of the lower bound of the same layer
    ALLOCATE(canopy_layer_thickness_profile(ncanopy))
    SELECT CASE(TRIM(config%canopy_layer_scheme))
      CASE ("standard")
        canopy_layer_thickness_profile(:) = 0.5_wp

      CASE ("fapar")
        ! alternative canopy layer depth calculation, approximately proportional to FAPAR to  
        ! reduce impact of actual canopy layering on PS calculation
        dAPAR = 1._wp/REAL(ncanopy,KIND=wp)
        DO icanopy=1,ncanopy
          canopy_layer_thickness_profile(icanopy) = -1.0_wp/k_canopy_layer * LOG(1._wp-(dAPAR*icanopy)+eps4)
        END DO

      CASE DEFAULT
        WRITE(message_text,'(2a)') 'invalid canopy_layer_scheme',canopy_layer_scheme
        CALL finish(routine, message_text) 
    END SELECT

    !< check canopy-conductance scheme
    !! 
    SELECT CASE(TRIM(config%canopy_conductance_scheme))
      CASE ("medlyn")
        print*, TRIM(routine)
        print*, '  using the canopy_conductance_scheme: ', TRIM(config%canopy_conductance_scheme), &
                '  the local lctlib_g1 variable is set to lctlib%g1_medlyn'
      CASE ("ballberry")
        print*, TRIM(routine)
        print*, '  using the canopy_conductance_scheme: ', TRIM(config%canopy_conductance_scheme), &
                '  the local lctlib_g1 variable is set to lctlib%g1_bberry'
      CASE DEFAULT
        WRITE(message_text,'(2a)') 'invalid canopy_conductance_scheme: ', TRIM(config%canopy_conductance_scheme)
        CALL finish(routine, message_text) 
    END SELECT

    !< Create vertical axis
    !! 
    !! these are "infrastructure" values, constant across PFTs and sites, and they do not interfere with tree height \n
    !! vgrid_canopy% lbounds, ubounds, and levels is calculated by 'new_vgrid()'
    !! 
    !! note: the axis of the canopy_layer vgrid is used as "reversed", i.e., layer 1 is the top canopy layer \n
    !! levels is defined as the depth at the center of the layer: levels(:) = 0.5_wp * (lbounds(:) + ubounds(:))
    vgrid_canopy  => new_vgrid('canopy_layer', ZAXIS_GENERIC, ncanopy, &
                      longname='Canopy layers (leaf mass per layer as LAI) from ASSIMI_ QUINCY', &
                      units='m2 m-2 LAI', &
                      ! levels=   , &
                      ! lbounds=  , & 
                      ! ubounds=  , &
                      dz=canopy_layer_thickness_profile(:))
    CALL register_vgrid(vgrid_canopy)

  END SUBROUTINE Init_assimi_config

#endif
END MODULE mo_assimi_config_class
