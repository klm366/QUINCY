!>
!! @brief Contains structures and methods for radiation config
!!
!! @author
!!  Rainer Schneck, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Rainer Schneck (2015-04-27)
!!
MODULE mo_rad_config_class
#ifndef __NO_JSBACH__
  ! -------------------------------------------------------------------------------------------------------
  ! Used variables of module

  USE mo_exception,         ONLY: message_text, message, finish
  USE mo_io_units,          ONLY: filename_max
  USE mo_kind,              ONLY: wp
  USE mo_jsb_config_class,  ONLY: t_jsb_config

  ! -------------------------------------------------------------------------------------------------------
  ! Module variables
  IMPLICIT NONE
  PRIVATE

  PUBLIC :: t_rad_config !, Get_rad_config

  TYPE, EXTENDS(t_jsb_config) :: t_rad_config
    LOGICAL           :: use_alb_veg_simple        !< Use the simplified albedo calculations for vegetation tiles
    LOGICAL           :: use_alb_soil_scheme       !< For vegetation and bare soil tiles:
                                                   !!   Use the complex calculations of the soil albedo from carbon in and
                                                   !!   on the soil (accunts for the soil under vegetation and for bare soil)
                                                   !!   if vegetation tile: only used if use_alb_veg_simple=.FALSE.
    LOGICAL           :: use_alb_mineralsoil_const !< For vegetation and bare soil tiles:
                                                   !!   if TRUE constant parameter values are used for albedo of the mineral soil
                                                   !!   if FALSE values for albedo of the mineral soil are taken from
                                                   !!   bc_land_phys.nc file.
    CHARACTER(len=10) :: use_alb_soil_organic_C    !< For vegetation and bare soil tiles:
                                                   !!   if linear: the albedo of soil without covering litter but with carbon in it
                                                   !!     is calculated from the mineral soil albedo minus a reduction calculated
                                                   !!     linearily from the carbon in the soil
                                                   !!   if log: the albedo of soil without covering litter but with carbon in it
                                                   !!     is calculated from a logarithmically decreasing function shaped by
                                                   !!     the mineral soil albedo and the carbon content within the soil
                                                   !!     linearily from the carbon in the soil
    LOGICAL           :: use_alb_soil_litter       !< For vegetation and bare soil tiles:
                                                   !!   TRUE: soil albedo includes the albedo of the litter cover above it
                                                   !!   FALSE soil albedo ignores the albedo of the litter cover above it

    LOGICAL           :: use_alb_canopy            !< Only for vegetation tiles:
                                                   !!   TRUE:  canopy albedo is taken from bc_land_phys.nc file
                                                   !!   FALSE: canopy albedo is taken from lctlib file
    REAL(wp)          :: albedo_age_weight         !< 0: snow age scheme with zero weight: only temperature scheme for snow albedo
                                                   !!      is used
                                                   !! 1: only snow age scheme is used, temperature scheme with zero weight
                                                   !! 0 < albedo_age_weight < 1: snow albedo is calculated by linearly weighting
                                                   !!                            the snow albedo resulting from both schemes
   CONTAINS
     PROCEDURE :: Init => Init_rad_config
  END type t_rad_config

  CHARACTER(len=*), PARAMETER :: modname = 'mo_rad_config_class'

CONTAINS
  ! -------------------------------------------------------------------------------------------------------
  !> Initialize radiation process
  !!
  !! @param[inout]     config     Configuration type of process (t_rad_config)
  !!
  SUBROUTINE Init_rad_config(config)

    USE mo_jsb_namelist_iface, ONLY: open_nml, POSITIONED, position_nml, close_nml

    CLASS(t_rad_config), INTENT(inout) :: config

    LOGICAL                     :: active
    CHARACTER(len=filename_max) :: ic_filename, bc_filename
    LOGICAL                     :: use_alb_veg_simple
    LOGICAL                     :: use_alb_mineralsoil_const
    CHARACTER(len=10)           :: use_alb_soil_organic_C
    LOGICAL                     :: use_alb_soil_litter
    LOGICAL                     :: use_alb_soil_scheme
    LOGICAL                     :: use_alb_canopy
    REAL(wp)                    :: albedo_age_weight

    NAMELIST /jsb_radiation_nml/      &
         active,                      &
         ic_filename,                 &
         bc_filename,                 &
         use_alb_veg_simple,          &
         use_alb_mineralsoil_const,   &
         use_alb_soil_organic_C,      &
         use_alb_soil_litter,         &
         use_alb_soil_scheme,         &
         use_alb_canopy,              &
         albedo_age_weight

    INTEGER :: nml_handler, nml_unit, istat

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':Init_rad_config'

    CALL message(TRIM(routine), 'Starting radiation configuration')

    ! Set defaults
    active              = .TRUE.
    bc_filename         = 'bc_land_radiation.nc'  !R: gibts im Moment noch nicht
    ic_filename         = 'ic_land_radiation.nc'  !R: gibts im Moment noch nicht
    use_alb_veg_simple        = .FALSE.
    use_alb_mineralsoil_const = .FALSE.
    use_alb_soil_organic_C    = 'linear'
    use_alb_soil_litter       = .TRUE.
    use_alb_soil_scheme       = .FALSE.
    use_alb_canopy            = .TRUE.
    albedo_age_weight         = 0.5_wp

    nml_handler = open_nml(TRIM(config%namelist_filename))

    nml_unit = position_nml('jsb_radiation_nml', nml_handler, STATUS=istat)
    IF (istat == POSITIONED) READ(nml_unit, jsb_radiation_nml)

    CALL close_nml(nml_handler)

    config%active              = active
    config%ic_filename         = ic_filename
    config%bc_filename         = bc_filename
    config%use_alb_veg_simple        = use_alb_veg_simple
    config%use_alb_mineralsoil_const = use_alb_mineralsoil_const
    config%use_alb_soil_organic_C    = use_alb_soil_organic_C
    config%use_alb_soil_litter       = use_alb_soil_litter
    config%use_alb_soil_scheme       = use_alb_soil_scheme
    config%use_alb_canopy            = use_alb_canopy
    config%albedo_age_weight         = albedo_age_weight

  END SUBROUTINE Init_rad_config

!!$  SUBROUTINE Get_rad_config(config, return_ptr)
!!$
!!$    CLASS(t_jsb_config), POINTER :: config
!!$    TYPE(t_rad_config),  POINTER :: return_ptr
!!$
!!$    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':Get_rad_config'
!!$
!!$    return_ptr => NULL()
!!$
!!$    IF (.NOT. ASSOCIATED(config)) RETURN
!!$
!!$    SELECT TYPE (config)
!!$    TYPE IS (t_rad_config)
!!$      return_ptr => config
!!$    END SELECT
!!$
!!$    IF (.NOT. ASSOCIATED(return_ptr)) &
!!$      CALL finish(TRIM(routine), 'Process has no configuration')
!!$
!!$  END SUBROUTINE Get_rad_config

#endif
END MODULE mo_rad_config_class
