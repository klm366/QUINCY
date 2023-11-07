!>
!! @brief Contains structures and methods for methane config
!!
!! @author
!!  Stiig Wilkenskjeld, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Stiig Wilkenskjeld (06.03.2023)
!!
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/license .
!!

MODULE mo_ch4_config_class
#ifndef __NO_JSBACH__
  ! -------------------------------------------------------------------------------------------------------
  ! Used variables of module

  USE mo_exception,         ONLY: message
  USE mo_io_units,          ONLY: filename_max
  USE mo_kind,              ONLY: wp
  USE mo_jsb_config_class,  ONLY: t_jsb_config

  ! -------------------------------------------------------------------------------------------------------
  ! Module variables
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: t_ch4_config

  TYPE, EXTENDS(t_jsb_config) :: t_ch4_config
    LOGICAL           :: dbg_process
    LOGICAL           :: has_peatland_tile             ! Is one of the tiles specific for peatlands?
    REAL(wp)          :: frac_CH4_anox_mineral         ! [-] fraction of anerobically decomposed carbon that becomes CH4 for mineral soils
    REAL(wp)          :: frac_CH4_anox_organic         ! [-] fraction of anerobically decomposed carbon that becomes CH4 for organic soils
    REAL(wp)          :: tile_threshold                ! tile size cutoff for CH4 calculations (0 - 0.01 resonable)
    REAL(wp)          :: meth_q10                      ! Q10 for CH4:CO2 ratio
    REAL(wp)          :: co2_conc                      ! Atmospheric concentration of CO2 [mol m-3]
    REAL(wp)          :: ch4_conc                      ! Atmospheric concentration of CH4 [mol m-3]
    REAL(wp)          ::  o2_conc                      ! Atmospheric concentration of  O2 [mol m-3]
    REAL(wp)          :: scale_termite                 ! scale termite CH4 output
    REAL(wp)          :: scale_herbivore               ! scale herbivore CH4 output
    REAL(wp)          :: scale_oxid_dry                ! scale CH4 oxidation dry part
    REAL(wp)          :: scale_oxid_wet                ! scale CH4 oxidation wet part
    REAL(wp)          :: implicit_contrib              ! contribution from implicit part of Crank-Nicolson diffusion equation
    REAL(wp)          :: frac_inun_min                 ! Threshold for tile fraction to be considered for methane processes
    LOGICAL           :: GHG_from_atm                  ! get GHG concentrations from ECHAM
    LOGICAL           :: do_fire                       ! determine fire CH4 emissions
    LOGICAL           :: do_termite                    ! determine termite CH4 emissions
    LOGICAL           :: do_herbivore                  ! determine herbivory CH4 emissions
    INTEGER           :: nstep_diffusion               ! number of diffusion steps per timestep
    INTEGER           :: nstep_per_day                 ! number of methane process time steps per day
    INTEGER           :: minimum_no_layers             ! Minimum number of soil layers to consider for methane
  CONTAINS
     PROCEDURE :: Init => Init_methane_config
  END type t_ch4_config

  CHARACTER(len=*), PARAMETER :: modname = 'mo_ch4_config_class'

CONTAINS
  ! -------------------------------------------------------------------------------------------------------
  !> Initialize methane process
  !!
  !! @param[inout]     config     Configuration type of process (t_carbon_config)
  !!
  SUBROUTINE Init_methane_config(config)

    USE mo_methane_util,       ONLY: convert_2_mol_m3, concentration_mol_m3_2_vmr
    USE mo_methane_constants,  ONLY: CH4_consts, CO2_consts, O2_consts
    USE mo_jsb_namelist_iface, ONLY: open_nml, POSITIONED, position_nml, close_nml
    USE mo_util,               ONLY: real2string, int2string, logical2string

    CLASS(t_ch4_config), INTENT(inout) :: config

    ! Local namelist parameters
    LOGICAL                     :: active
    CHARACTER(len=filename_max) :: ic_filename, bc_filename
    LOGICAL                     :: dbg_process
    LOGICAL                     :: has_peatland_tile         ! Is one of the tiles specific for peatlands?
    REAL(wp)                    :: frac_CH4_anox_mineral     ! [-] fraction of anerobically decomposed carbon that becomes CH4 for mineral soils
    REAL(wp)                    :: frac_CH4_anox_organic     ! [-] fraction of anerobically decomposed carbon that becomes CH4 for organic soils
    REAL(wp)                    :: tile_threshold            ! [-] tile size cutoff for CH4 calculations (0 - 0.01 resonable)
    REAL(wp)                    :: meth_q10                  ! Q10 for CH4:CO2 ratio
    REAL(wp)                    :: co2_vmr                   ! VMR of CO2 [m3 m-3]
    REAL(wp)                    :: ch4_vmr                   ! VMR of CH4 [m3 m-3]
    REAL(wp)                    ::  O2_vmr                   ! VMR of O2  [m3 m-3]
    REAL(wp)                    :: co2_mmr                   ! MMR of CO2 [kg kg-1]
    REAL(wp)                    :: ch4_mmr                   ! MMR of CH4 [kg kg-1]
    REAL(wp)                    ::  O2_mmr                   ! VMR of O2  [m3 m-3]
    REAL(wp)                    :: co2_conc                  ! Concentration of CO2 [mol m-3]
    REAL(wp)                    :: ch4_conc                  ! Concentration of CH4 [mol m-3]
    REAL(wp)                    ::  O2_conc                  ! VMR of O2 [m3 m-3]
    REAL(wp)                    :: scale_termite             ! scale termite CH4 output
    REAL(wp)                    :: scale_herbivore           ! scale herbivore CH4 output
    REAL(wp)                    :: scale_oxid_dry            ! scale CH4 oxidation dry part
    REAL(wp)                    :: scale_oxid_wet            ! scale CH4 oxidation wet part
    REAL(wp)                    :: implicit_contrib          ! contribution from implicit part of Crank-Nicolson diffusion equation
    REAL(wp)                    :: frac_inun_min             ! Threshold for tile fraction to be considered for methane processes
    LOGICAL                     :: GHG_from_atm              ! get GHG concentration from ECHAM
    LOGICAL                     :: do_fire                   ! determine fire CH4 emissions
    LOGICAL                     :: do_termite                ! determine termite CH4 emissions
    LOGICAL                     :: do_herbivore              ! determine herbivory CH4 emissions
    INTEGER                     :: minimum_no_layers         ! Minimum number of soil layers to consider for methane
    INTEGER                     :: nstep_per_day             ! Number of methane process time steps per day
    INTEGER                     :: nstep_diffusion           ! Number of diffusion time steps per methane process time step

    NAMELIST /jsb_methane_nml/     &
         active,                   &
         ic_filename,              &
         bc_filename,              &
         dbg_process,              &
         has_peatland_tile,        &
         frac_CH4_anox_mineral,    &
         frac_CH4_anox_organic,    &
         tile_threshold,           &
         meth_q10,                 &
         co2_vmr,                  &
         ch4_vmr,                  &
         co2_mmr,                  &
         ch4_mmr,                  &
         co2_conc,                 &
         ch4_conc,                 &
         scale_termite,            &
         scale_herbivore,          &
         scale_oxid_dry,           &
         scale_oxid_wet,           &
         implicit_contrib,         &
         GHG_from_atm,             &
         do_fire,                  &
         do_termite,               &
         do_herbivore,             &
         minimum_no_layers,        &
         frac_inun_min,            &
         nstep_per_day,            &
         nstep_diffusion

    INTEGER :: nml_handler, nml_unit, istat

    CHARACTER(len=*), PARAMETER :: routine = modname//':Init_methane_config'

    CALL message(TRIM(routine), 'Starting methane configuration')

    ! Set defaults
    active                 = .FALSE.
    bc_filename            = 'bc_land_carbon.nc'
    ic_filename            = 'ic_land_carbon.nc'
    dbg_process            = .FALSE.
    has_peatland_tile      = .FALSE.
    frac_CH4_anox_mineral  =    0.1_wp
    frac_CH4_anox_organic  =    0.1_wp
    tile_threshold         =    0.0_wp
    meth_q10               =    1.0_wp
    CO2_vmr                =  400.e-6_wp
    CH4_vmr                = 1885.e-9_wp
     O2_vmr                =    2.0946e-1_wp
    CO2_mmr                =   -1.0_wp ! Impossible value to determine specification
    CH4_mmr                =   -1.0_wp ! Impossible value to determine specification
     O2_mmr                =   -1.0_wp ! Impossible value to determine specification
    CO2_conc               =   -1.0_wp ! Impossible value to determine specification
    CH4_conc               =   -1.0_wp ! Impossible value to determine specification
     O2_conc               =   -1.0_wp ! Impossible value to determine specification
    scale_termite          =    1.0_wp
    scale_herbivore        =    1.0_wp
    scale_oxid_dry         =    1.0_wp
    scale_oxid_wet         =    1.0_wp
    implicit_contrib       =    0.5_wp
    frac_inun_min          =    1.e-3_wp
    minimum_no_layers      =    3
    nstep_per_day          =    2
    nstep_diffusion        =    1
    GHG_from_atm           = .FALSE.
    do_fire                = .TRUE.
    do_termite             = .TRUE.
    do_herbivore           = .TRUE.

    ! Read namelist
    nml_handler = open_nml(TRIM(config%namelist_filename))

    nml_unit = position_nml('jsb_methane_nml', nml_handler, STATUS=istat)
    IF (istat == POSITIONED) READ(nml_unit, jsb_methane_nml)

    CALL close_nml(nml_handler)

    ! Process units of CH4, CO2 and O2
    IF (.NOT. GHG_from_atm) THEN
      CALL convert_2_mol_m3(CO2_consts,CO2_mmr,CO2_vmr,CO2_conc)
!      CALL convert_2_mol_m3(CH4_consts,CH4_mmr,CH4_vmr,CH4_conc)
!      CALL convert_2_mol_m3( O2_consts, O2_mmr, O2_vmr, O2_conc)
    ENDIF
    ! Presently only CO2 concentration can be provided by the atmosphere
    CALL convert_2_mol_m3(CH4_consts,CH4_mmr,CH4_vmr,CH4_conc)
    CALL convert_2_mol_m3( O2_consts, O2_mmr, O2_vmr, O2_conc)

    ! Write resulting values in config memory
    config%active                 = active
    config%ic_filename            = ic_filename
    config%bc_filename            = bc_filename
    config%dbg_process            = dbg_process
    config%run_case               = 1
    IF (has_peatland_tile) THEN
      config%run_case               = 2
    ENDIF
    config%has_peatland_tile      = has_peatland_tile
    config%frac_CH4_anox_mineral  = frac_CH4_anox_mineral
    config%frac_CH4_anox_organic  = frac_CH4_anox_organic
    config%tile_threshold         = tile_threshold
    config%meth_q10               = meth_q10
    config%co2_conc               = co2_conc
    config%ch4_conc               = ch4_conc
    config%o2_conc                = o2_conc
    config%scale_termite          = scale_termite
    config%scale_herbivore        = scale_herbivore
    config%scale_oxid_dry         = scale_oxid_dry
    config%scale_oxid_wet         = scale_oxid_wet
    config%implicit_contrib       = implicit_contrib
    config%GHG_from_atm           = GHG_from_atm
    config%do_fire                = do_fire
    config%do_termite             = do_termite
    config%do_herbivore           = do_herbivore
    config%nstep_diffusion        = nstep_diffusion
    config%nstep_per_day          = nstep_per_day
    config%minimum_no_layers      = minimum_no_layers
    config%frac_inun_min          = frac_inun_min

    IF (active) THEN
      CALL message(routine,'Running CH4 soil transport model with parameters:')
      CALL message(routine,'CH4 production fraction organic soil '//real2string(frac_CH4_anox_organic))
      CALL message(routine,'tile size threshold '//real2string(tile_threshold))
      CALL message(routine,'Q10 for methane production fraction '//real2string(meth_q10))
      CALL message(routine,'implicit contribution in Crank-Nicolson diffusion '//real2string(implicit_contrib))
      CALL message(routine,'methane time steps per day '//int2string(nstep_per_day))
      CALL message(routine,'diffusion steps per timestep '//int2string(nstep_diffusion))
      CALL message(routine,'fire methane emissions '//logical2string(do_fire))
      CALL message(routine,'termite methane emissions '//logical2string(do_termite))
      CALL message(routine,'herbivore methane emissions '//logical2string(do_herbivore))
      IF (GHG_from_atm) THEN
        CALL message(routine,'Using ECHAM values for GHG mmr')
      ELSE
        CALL message(routine,'Using namelist values for GHG vmr')
        CALL message(routine,'CO2 VMR: '//real2string(concentration_mol_m3_2_vmr(CO2_consts,CO2_conc)))
        CALL message(routine,'CH4 VMR: '//real2string(concentration_mol_m3_2_vmr(CH4_consts,CH4_conc)))
      END IF
      CALL message(routine,'O2 VMR (namelist): '//real2string(concentration_mol_m3_2_vmr(O2_consts,O2_conc)))
      IF (scale_termite /= 1.0_wp) THEN
        CALL message(routine,'Scaling termite CH4 production by '//real2string(scale_termite))
      END IF
      IF (scale_herbivore /= 1.0_wp) THEN
        CALL message(routine,'Scaling herbivore CH4 production by '//real2string(scale_herbivore))
      END IF
      IF (scale_oxid_dry /= 1.0_wp) THEN
        CALL message(routine,'Scaling dry CH4 oxidation by '//real2string(scale_oxid_dry))
      END IF
      IF (scale_oxid_wet /= 1.0_wp) THEN
        CALL message(routine,'Scaling wet CH4 oxidation by '//real2string(scale_oxid_wet))
      END IF
    END IF

  END SUBROUTINE Init_methane_config

#endif
END MODULE mo_ch4_config_class
