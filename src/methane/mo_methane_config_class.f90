!> Contains structures and methods for methane config
!>
!> ICON-Land
!>
!> ---------------------------------------
!> Copyright (C) 2013-2026, MPI-M, MPI-BGC
!>
!> Contact: icon-model.org
!> Authors: AUTHORS.md
!> See LICENSES/ for license information
!> SPDX-License-Identifier: BSD-3-Clause
!> ---------------------------------------
!>
MODULE mo_ch4_config_class
#ifndef __NO_JSBACH__
  ! -------------------------------------------------------------------------------------------------------
  ! Used variables of module

  USE mo_exception,         ONLY: message, finish
  USE mo_io_units,          ONLY: filename_max
  USE mo_kind,              ONLY: wp
  USE mo_jsb_config_class,  ONLY: t_jsb_config

  ! -------------------------------------------------------------------------------------------------------
  ! Module variables
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: t_ch4_config

  TYPE, EXTENDS(t_jsb_config) :: t_ch4_config
    CHARACTER(len=32) :: CH4_varname                   ! Variable name of CH4 in netcdf file for CH4_source == 2
    CHARACTER(len=32) :: CO2_varname                   ! Variable name of CO2 in netcdf file for CO2_source == 2
    CHARACTER(len=32) ::  O2_varname                   ! Variable name of  O2 in netcdf file for  O2_source == 2
    LOGICAL           :: dbg_process
    LOGICAL           :: has_peatland_tile             ! Is one of the tiles specific for peatlands?
    REAL(wp)          :: frac_CH4_anox_mineral         ! [-] fraction of anerobically decomposed carbon that becomes CH4 for mineral soils
    REAL(wp)          :: frac_CH4_anox_organic         ! [-] fraction of anerobically decomposed carbon that becomes CH4 for organic soils
    REAL(wp)          :: tile_threshold                ! tile size cutoff for CH4 calculations (0 - 0.01 resonable)
    REAL(wp)          :: meth_q10                      ! Q10 for CH4:CO2 ratio
    REAL(wp)          :: oxid_q10                      ! Q10 for CH4 oxidation
    REAL(wp)          :: co2_conc                      ! Atmospheric concentration of CO2 [mol m-3]
    REAL(wp)          :: ch4_conc                      ! Atmospheric concentration of CH4 [mol m-3]
    REAL(wp)          ::  o2_conc                      ! Atmospheric concentration of  O2 [mol m-3]
    REAL(wp)          :: co2_scale                     ! Manual unit rescale factor for read timeseries of co2
    REAL(wp)          :: ch4_scale                     ! Manual unit rescale factor for read timeseries of ch4
    REAL(wp)          ::  o2_scale                     ! Manual unit rescale factor for read timeseries of o2
    REAL(wp)          :: scale_termite                 ! scale termite CH4 output
    REAL(wp)          :: scale_herbivore               ! scale herbivore CH4 output
    REAL(wp)          :: scale_oxid_dry                ! scale CH4 oxidation dry part
    REAL(wp)          :: scale_oxid_wet                ! scale CH4 oxidation wet part
    REAL(wp)          :: implicit_contrib              ! contribution from implicit part of Crank-Nicolson diffusion equation
    REAL(wp)          :: frac_inun_min                 ! Threshold for tile fraction to be considered for methane processes
    REAL(wp)          :: polar_anoxic_fac              ! Multiply polar anoxic respiration by this factor
    LOGICAL           :: do_fire                       ! determine fire CH4 emissions
    LOGICAL           :: do_termite                    ! determine termite CH4 emissions
    LOGICAL           :: do_herbivore                  ! determine herbivory CH4 emissions
    INTEGER           :: nstep_diffusion               ! number of diffusion steps per timestep
    INTEGER           :: nstep_per_day                 ! number of methane process time steps per day
    INTEGER           :: minimum_no_layers             ! Minimum number of soil layers to consider for methane
    INTEGER           :: CH4_source                    ! Source of atmospheric CH4 (opt 2 not available yet)
    INTEGER           :: CO2_source                    ! Source of atmospheric CO2
    INTEGER           ::  O2_source                    ! Source of atmospheric O2 (opt 2 not available yet)
    ! Source type for GHG are:
    ! 1 : Constant in time and space from this namelist
    ! 2 : Temporarily and spatially varying field from atmospheric host model (not available for all GHGs yet)
    ! 3 : Timeseries with yearly vmr values from file
    ! 4 : Timeseries with yearly mmr values from file
    ! 5 : Timeseries with yearly conc values from file
    ! Presently all timeseries of GHGs need to be in the same file (bc_filename)
  CONTAINS
     PROCEDURE :: Init => Init_methane_config
  END type t_ch4_config

  CHARACTER(len=*), PARAMETER :: modname = 'mo_ch4_config_class'

CONTAINS
  SUBROUTINE methane_GHG_atm_proc_and_msg(specie,source,mul,mmr,vmr,conc)
    USE mo_methane_constants,  ONLY: t_specie_const
    USE mo_methane_util,       ONLY: convert_2_mol_m3, concentration_mol_m3_2_vmr
    USE mo_util,               ONLY: real2string

    TYPE(t_specie_const), INTENT(in)    :: specie
    INTEGER,              INTENT(in)    :: source
    REAL(wp),             INTENT(in)    :: mul
    REAL(wp),             INTENT(in)    :: mmr
    REAL(wp),             INTENT(in)    :: vmr
    REAL(wp),             INTENT(inout) :: conc

    CHARACTER(len=*), PARAMETER :: routine = modname//':methane_GHG_atm_proc_and_msg'

    ! Process units of CH4, CO2 and O2
    SELECT CASE (source)
      CASE (1) ! Constant from namelist: Convert from whatever specified
        IF (conc >= 0._wp) conc = conc * mul
        CALL convert_2_mol_m3(specie,mmr*mul,vmr*mul,conc)
        CALL message(routine,'Atmospheric value of '//TRIM(specie%specie_name)//' is constant at vmr: '// &
                             real2string(concentration_mol_m3_2_vmr(specie,conc)))
      CASE (2) ! Take time and spatially dependent input from atm (presently only available for CO2)
        CALL message(routine,'Using atmospheric values of '//TRIM(specie%specie_name)//' from host model')
      CASE (3) ! Time dependent input from file
        CALL message(routine,'Taking atmospheric vmr values of '//TRIM(specie%specie_name)//' from time series file')
      CASE (4) ! Time dependent input from file
        CALL message(routine,'Taking atmospheric mmr values of '//TRIM(specie%specie_name)//' from time series file')
      CASE (5) ! Time dependent input from file
        CALL message(routine,'Taking atmospheric concentration values of '//TRIM(specie%specie_name)//' from time series file')
    END SELECT

  END SUBROUTINE methane_GHG_atm_proc_and_msg

  ! -------------------------------------------------------------------------------------------------------
  !> Initialize methane process
  !!
  !! @param[inout]     config     Configuration type of process (t_carbon_config)
  !!
  SUBROUTINE Init_methane_config(config)

    USE mo_methane_constants,  ONLY: CH4_consts, CO2_consts, O2_consts
    USE mo_jsb_namelist_iface, ONLY: open_nml, POSITIONED, position_nml, close_nml
    USE mo_util,               ONLY: real2string, int2string, logical2string

    CLASS(t_ch4_config), INTENT(inout) :: config

    ! Local namelist parameters
    LOGICAL                     :: active
    CHARACTER(len=filename_max) :: ic_filename, bc_filename
    CHARACTER(len=32)           :: CH4_varname
    CHARACTER(len=32)           :: CO2_varname
    CHARACTER(len=32)           ::  O2_varname
    LOGICAL                     :: dbg_process
    LOGICAL                     :: has_peatland_tile         ! Is one of the tiles specific for peatlands?
    REAL(wp)                    :: frac_CH4_anox_mineral     ! [-] fraction of anerobically decomposed carbon that becomes CH4 for mineral soils
    REAL(wp)                    :: frac_CH4_anox_organic     ! [-] fraction of anerobically decomposed carbon that becomes CH4 for organic soils
    REAL(wp)                    :: tile_threshold            ! [-] tile size cutoff for CH4 calculations (0 - 0.01 resonable)
    REAL(wp)                    :: meth_q10                  ! Q10 for CH4:CO2 ratio
    REAL(wp)                    :: oxid_q10                  ! Q10 for CH4 oxidation
    REAL(wp)                    :: co2_vmr                   ! VMR of CO2 [m3 m-3]
    REAL(wp)                    :: ch4_vmr                   ! VMR of CH4 [m3 m-3]
    REAL(wp)                    ::  O2_vmr                   ! VMR of O2  [m3 m-3]
    REAL(wp)                    :: co2_mmr                   ! MMR of CO2 [kg kg-1]
    REAL(wp)                    :: ch4_mmr                   ! MMR of CH4 [kg kg-1]
    REAL(wp)                    ::  O2_mmr                   ! VMR of O2  [m3 m-3]
    REAL(wp)                    :: co2_conc                  ! Concentration of CO2 [mol m-3]
    REAL(wp)                    :: ch4_conc                  ! Concentration of CH4 [mol m-3]
    REAL(wp)                    ::  O2_conc                  ! Concentration of  O2 [m3 m-3]
    REAL(wp)                    :: co2_scale                 ! Manual unit rescale factor for read timeseries of co2
    REAL(wp)                    :: ch4_scale                 ! Manual unit rescale factor for read timeseries of ch4
    REAL(wp)                    ::  o2_scale                 ! Manual unit rescale factor for read timeseries of o2
    REAL(wp)                    :: scale_termite             ! scale termite CH4 output
    REAL(wp)                    :: scale_herbivore           ! scale herbivore CH4 output
    REAL(wp)                    :: scale_oxid_dry            ! scale CH4 oxidation dry part
    REAL(wp)                    :: scale_oxid_wet            ! scale CH4 oxidation wet part
    REAL(wp)                    :: implicit_contrib          ! contr from implicit part of Crank-Nicolson diffusion equation
    REAL(wp)                    :: frac_inun_min             ! Threshold for tile fraction to be considered for methane processes
    REAL(wp)                    :: polar_anoxic_fac          ! Multiply pole-near anoxic respiration by this factor
    LOGICAL                     :: do_fire                   ! determine fire CH4 emissions
    LOGICAL                     :: do_termite                ! determine termite CH4 emissions
    LOGICAL                     :: do_herbivore              ! determine herbivory CH4 emissions
    INTEGER                     :: CH4_source                ! Source for atmospheric CH4 (see above)
    INTEGER                     :: CO2_source                ! Source for atmospheric CO2 (see above)
    INTEGER                     ::  O2_source                ! Source for atmospheric O2 (see above)
    INTEGER                     :: minimum_no_layers         ! Minimum number of soil layers to consider for methane
    INTEGER                     :: nstep_per_day             ! Number of methane process time steps per day
    INTEGER                     :: nstep_diffusion           ! Number of diffusion time steps per methane process time step

    NAMELIST /jsb_methane_nml/     &
         active,                   &
         ic_filename,              &
         bc_filename,              &
         CH4_varname,              &
         dbg_process,              &
         has_peatland_tile,        &
         frac_CH4_anox_mineral,    &
         frac_CH4_anox_organic,    &
         tile_threshold,           &
         meth_q10,                 &
         oxid_q10,                 &
         co2_vmr,                  &
         ch4_vmr,                  &
         co2_mmr,                  &
         ch4_mmr,                  &
         co2_conc,                 &
          o2_conc,                 &
         ch4_conc,                 &
         co2_scale,                &
          o2_scale,                &
         ch4_scale,                &
         ch4_source,               &
         co2_source,               &
          o2_source,               &
         scale_termite,            &
         scale_herbivore,          &
         scale_oxid_dry,           &
         scale_oxid_wet,           &
         implicit_contrib,         &
         do_fire,                  &
         do_termite,               &
         do_herbivore,             &
         minimum_no_layers,        &
         frac_inun_min,            &
         polar_anoxic_fac,         &
         nstep_per_day,            &
         nstep_diffusion

    INTEGER :: nml_handler, nml_unit, istat

    CHARACTER(len=*), PARAMETER :: routine = modname//':Init_methane_config'

    CALL message(TRIM(routine), 'Starting methane configuration')

    ! Set defaults
    active                 = .FALSE.
    bc_filename            = 'bc_land_methane.nc'
    ic_filename            = 'ic_land_methane.nc'
    CH4_varname            = 'CH4'
    CO2_varname            = 'CO2'
     O2_varname            = 'O2'
    dbg_process            = .FALSE.
    has_peatland_tile      = .FALSE.
    frac_CH4_anox_mineral  =    0.1_wp
    frac_CH4_anox_organic  =    0.1_wp
    tile_threshold         =    0.0_wp
    meth_q10               =    1.0_wp
    oxid_q10               =    1.9_wp
    CO2_vmr                =  400._wp
    CH4_vmr                = 1885._wp
     O2_vmr                =   20.946_wp
    CO2_mmr                =   -1.0_wp    ! Impossible value to determine specification
    CH4_mmr                =   -1.0_wp    ! Impossible value to determine specification
     O2_mmr                =   -1.0_wp    ! Impossible value to determine specification
    CO2_conc               =   -1.0_wp    ! Impossible value to determine specification
    CH4_conc               =   -1.0_wp    ! Impossible value to determine specification
     O2_conc               =   -1.0_wp    ! Impossible value to determine specification
    CH4_scale              =    1.e-9_wp  ! Default unit is ppbv
    CO2_scale              =    1.e-6_wp  ! Default unit is ppmv
     O2_scale              =    1.e-2_wp  ! Default unit is vol%
    CH4_source             =    1
    CO2_source             =    1
     O2_source             =    1
    scale_termite          =    1.0_wp
    scale_herbivore        =    1.0_wp
    scale_oxid_dry         =    1.0_wp
    scale_oxid_wet         =    1.0_wp
    implicit_contrib       =    0.5_wp
    frac_inun_min          =    1.e-3_wp
    polar_anoxic_fac       =    1.0_wp
    minimum_no_layers      =    3
    nstep_per_day          =    2
    nstep_diffusion        =    1
    do_fire                = .TRUE.
    do_termite             = .TRUE.
    do_herbivore           = .TRUE.

    ! Read namelist
    nml_handler = open_nml(TRIM(config%namelist_filename))

    nml_unit = position_nml('jsb_methane_nml', nml_handler, STATUS=istat)
    IF (istat == POSITIONED) READ(nml_unit, jsb_methane_nml)

    CALL close_nml(nml_handler)

    ! Consistency checks
    IF (CH4_source < 1 .OR. CH4_source > 3) &
      CALL finish(routine,'Bad specification of atmospheric CH4 source')

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
    config%oxid_q10               = oxid_q10
    config%scale_termite          = scale_termite
    config%scale_herbivore        = scale_herbivore
    config%scale_oxid_dry         = scale_oxid_dry
    config%scale_oxid_wet         = scale_oxid_wet
    config%implicit_contrib       = implicit_contrib
    config%CH4_source             = CH4_source
    config%CO2_source             = CO2_source
    config%O2_source              =  O2_source
    config%CH4_scale             = CH4_scale
    config%CO2_scale             = CO2_scale
    config%O2_scale              =  O2_scale
    config%CH4_varname            = CH4_varname
    config%CO2_varname            = CO2_varname
    config%O2_varname             =  O2_varname
    config%do_fire                = do_fire
    config%do_termite             = do_termite
    config%do_herbivore           = do_herbivore
    config%nstep_diffusion        = nstep_diffusion
    config%nstep_per_day          = nstep_per_day
    config%minimum_no_layers      = minimum_no_layers
    config%frac_inun_min          = frac_inun_min
    config%polar_anoxic_fac       = polar_anoxic_fac

    IF (active) THEN
      CALL message(routine,'Running CH4 soil transport model with parameters:')
      CALL message(routine,'CH4 production fraction organic soil '//real2string(frac_CH4_anox_organic))
      CALL message(routine,'tile size threshold '//real2string(tile_threshold))
      CALL message(routine,'Q10 for methane production fraction '//real2string(meth_q10))
      CALL message(routine,'Q10 for methane oxidation           '//real2string(oxid_q10))
      CALL message(routine,'implicit contribution in Crank-Nicolson diffusion '//real2string(implicit_contrib))
      CALL message(routine,'methane time steps per day '//int2string(nstep_per_day))
      CALL message(routine,'diffusion steps per timestep '//int2string(nstep_diffusion))
      CALL message(routine,'fire methane emissions '//logical2string(do_fire))
      CALL message(routine,'termite methane emissions '//logical2string(do_termite))
      CALL message(routine,'herbivore methane emissions '//logical2string(do_herbivore))
      IF (CH4_source==2) CALL finish(routine,'Atmospheric CH4 from atm model not yet available')
      IF (O2_source ==2) CALL finish(routine,'Atmospheric O2 from atm model not yet available')
      CALL methane_GHG_atm_proc_and_msg(CH4_consts,CH4_source,CH4_scale,CH4_mmr,CH4_vmr,CH4_conc)
      CALL methane_GHG_atm_proc_and_msg(CO2_consts,CO2_source,CO2_scale,CO2_mmr,CO2_vmr,CO2_conc)
      CALL methane_GHG_atm_proc_and_msg(O2_consts, O2_source,  O2_scale,O2_mmr, O2_vmr,  O2_conc)
      config%ch4_conc = CH4_conc
      config%co2_conc = CO2_conc
      config%o2_conc  = O2_conc
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
