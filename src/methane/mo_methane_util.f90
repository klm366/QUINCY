!>
!! @brief Contains utitilties for the methane production and transport
!!
!! @author
!!  Stiig Wilkenskjeld, MPI-M Hamburg
!!
!! @par Revision History
!! Adopted from JSBACH3 code by Thomas Kleinen, 2023-21-03
!!
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/license .
!!

MODULE mo_methane_util
#ifndef __NO_JSBACH__

  USE mo_kind,                   ONLY: wp
  USE mo_methane_constants,      ONLY: T_ref, t_specie_const
  USE mo_jsb_physical_constants, ONLY: p_sealevel => p0sl_bg, argas, molarMassDryAir, tmelt
  USE mo_exception,              ONLY: finish, message, message_text !, em_warn

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: find_layer, find_bedrock_layer, add_surface_layer
  PUBLIC :: dissolved_gas, determine_henry_cc
  PUBLIC :: concentration_vmr_2_mol_m3, concentration_mmr_2_mol_m3
  PUBLIC :: concentration_mol_m3_2_mmr, concentration_vmr_2_mmr
  PUBLIC :: concentration_mol_m3_2_vmr, concentration_mmr_2_vmr
  PUBLIC :: convert_2_vmr, convert_2_mmr, convert_2_mol_m3
#ifndef _OPENACC
  PUBLIC :: calc_diffusivity, diffusion_cn
#endif

  CHARACTER(len=*), PARAMETER :: modname = 'mo_methane_util'

CONTAINS

  !!
  ! @brief      Find layer of a given depth (eg. roots end, water table or acrotelm)
  !
  !!
!#ifndef _OPENACC
  PURE &
!#endif
  INTEGER FUNCTION find_layer(depth, nsoil, midpoints)
  !$acc routine seq
    INTEGER,  INTENT(in)  :: nsoil
    REAL(wp), INTENT(in)  :: depth, midpoints(nsoil)

    ! locals
    CHARACTER(len=*), PARAMETER :: routine = modname//':find_layer'
    INTEGER :: iLayer                        ! layer counter

    find_layer = 0                ! default value -> at or above surface

    IF (depth >= midpoints(1)) THEN
      find_layer = nsoil
      DO iLayer = 1, nsoil-1
        IF ((depth >= midpoints(iLayer)) .AND. (depth < midpoints(iLayer+1))) THEN
          find_layer = iLayer
        END IF
      END DO
    END IF

  END FUNCTION find_layer

  !!
  ! @brief      Find layer of a given depth
  !
  !!
!#ifndef _OPENACC
  PURE &
!#endif
  INTEGER FUNCTION find_bedrock_layer(nsoil, layer_bound, soil_depth)
  !$acc routine seq
    INTEGER,                    INTENT(in)  :: nsoil
    REAL(wp), DIMENSION(nsoil), INTENT(in)  :: layer_bound
    REAL(wp),                   INTENT(in)  :: soil_depth

    ! locals
    CHARACTER(len=*), PARAMETER :: routine = modname//':find_bedrock_layer'
    INTEGER :: bedrock_layer
    INTEGER :: iLayer                        ! layer counter

    bedrock_layer = nsoil
    DO iLayer = nsoil,2,-1
      IF(layer_bound(iLayer-1) >= soil_depth) THEN
        bedrock_layer = iLayer - 1
      ELSE
        IF ((soil_depth - layer_bound(iLayer-1)) < 1.e-3_wp) THEN
          bedrock_layer = bedrock_layer - 1
        END IF
      END IF
    END DO

    find_bedrock_layer = bedrock_layer

  END FUNCTION find_bedrock_layer

  !!
  ! @brief  determine gas concentration in water [mol m-2]
  ! based on Henry's law, as described by Wania
  !!
#ifndef _OPENACC
  PURE &
#endif
  REAL(wp) FUNCTION dissolved_gas(GAS, Tsoil, GAS_dens, dens2conc)
  !$acc routine seq
    TYPE (t_specie_const),  INTENT(in)  :: GAS       ! Type with constants for specific gas
    REAL(wp),               INTENT(in)  :: Tsoil     ! Soil temperature in layer [K]
    REAL(wp),               INTENT(in)  :: GAS_dens  ! GAS density in ice free pore space of layer [mol]
    REAL(wp),               INTENT(in)  :: dens2conc ! conversion factor from density to concentration []

    ! locals
    CHARACTER(len=*), PARAMETER :: routine = modname//':dissolved_gas'
    REAL(wp) :: GAS_concentration ! concentration of GAS in layer [mol m-3]
    REAL(wp) :: kHCC              ! Henry coefficient relating aqueous to gaseous concentration []
    REAL(wp) :: GAS_conc_diss     ! concentration of GAS dissolved in water [mol m-3]

    GAS_concentration = GAS_dens * dens2conc

    kHCC = determine_Henry_cc(GAS, Tsoil)

    ! Henry's law: kHCC = ca/cg -> ca = kHCC * cg
    ! ca: aqueous concentration cg: gaseous concentration

    GAS_conc_diss = kHCC * GAS_concentration
    dissolved_gas = GAS_conc_diss / dens2conc

  END FUNCTION dissolved_gas

  !!
  ! @brief  determine gas concentration in water based on Henry's law
  ! result is dimensionless Henry coefficient relating aqueous to gaseous concentration
  !!
!#ifndef _OPENACC
  PURE &
!#endif
  REAL(wp) FUNCTION determine_henry_cc(GAS, temp)
  !$acc routine seq
    TYPE (t_specie_const), INTENT(in) :: GAS            ! Structure with constanst for a specific gass
    REAL(wp),              INTENT(in) :: temp           ! Soil temperature in layer [K]

    ! Parameters
    CHARACTER(len=*), PARAMETER :: routine = modname//':determine_Henry_cc'
    REAL(wp), PARAMETER :: atm = 1._wp / p_sealevel  ! [atm/Pa] conversion factor for converting Pa in atm
    REAL(wp), PARAMETER :: r_gas_atm_l = argas * atm * 1000._wp  ! gas constant in [atm*l/K/mol]

    ! locals
    REAL(wp)                          :: henry_coeff_inv ! Henry's T dependent coefficient, inverted
    REAL(wp)                          :: henry_coeff     ! non-inverted Henry coefficient

    ! Calculate concentration of GAS in equilibrium with the
    ! atmosphere using Henry's law
    ! Henry's constant/coefficient depends on T
    ! UNITS: henry_k [L atm mol-1], const_GAS [K], T [K]
    ! henry_coeff_inv [L atm mol-1]
    ! -------------------------------------------------------

    henry_coeff_inv = GAS%henry_k_inv * EXP(-GAS%henry_const * &
                      (1.0_wp / temp - 1.0_wp / (tmelt + 25._wp)))
    henry_coeff = 1._wp / henry_coeff_inv

    ! Henry's law: kHCC = ca/cg -> ca = kHCC * cg
    ! ca: aqueous concentration cg: gaseous concentration
    determine_Henry_cc = henry_coeff * temp * r_gas_atm_l

  END FUNCTION determine_henry_cc

  !!
  ! @brief Converts a gas concentration [mol m-3] to volume mixing ration at reference temperature and pressure
  !
  !!
#ifndef _OPENACC
  PURE &
#endif
  REAL(wp) FUNCTION concentration_mol_m3_2_vmr(GAS,GAS_conc)
  !$acc routine seq
    TYPE (t_specie_const), INTENT(in) :: GAS      ! Type with constants specific for the gas
    REAL(wp),              INTENT(in) :: GAS_conc ! Gas concentration [mol m-3]
    CHARACTER(len=*), PARAMETER :: routine = modname//':concentration_mol_m3_2_vmr'
    concentration_mol_m3_2_vmr = GAS_conc * (argas * T_ref / p_sealevel)
  END FUNCTION concentration_mol_m3_2_vmr

  !!
  ! @brief Converts a gas volume mixing ration to concentration [mol m-3] at reference temperature and pressure
  !
  !!
#ifndef _OPENACC
  PURE &
#endif
  REAL(wp) FUNCTION concentration_vmr_2_mol_m3(GAS,GAS_vmr)
  !$acc routine seq
    TYPE (t_specie_const), INTENT(in) :: GAS     ! Type with constants specific for the gas
    REAL(wp),              INTENT(in) :: GAS_vmr ! Volume mixing ratio of gas [m3 m-3]
    CHARACTER(len=*), PARAMETER :: routine = modname//':concentration_vmr_2_mol_m3'
    ! from ideal gas law: PV = nRT <=> V/n = RT/P <=> n = V/(RT/P)
    ! volume mixing ratio is molar mixing ratio (1 mol gas = constant volume, independent of gas)
    ! => determine # of moles in volume
    concentration_vmr_2_mol_m3 = GAS_vmr / (argas * T_ref / p_sealevel)

  END FUNCTION concentration_vmr_2_mol_m3

  !!
  ! @brief Converts a gas concentration to mass mixing ratio at reference temperature and pressure
  !
  !!
#ifndef _OPENACC
  PURE &
#endif
  REAL(wp) FUNCTION concentration_mol_m3_2_mmr(GAS, GAS_conc)
  !$acc routine seq
    TYPE (t_specie_const), INTENT(in) :: GAS      ! Type with constants specific for the gas
    REAL(wp),              INTENT(in) :: GAS_conc ! Concentration of gas [mol m-3]
    CHARACTER(len=*), PARAMETER :: routine = modname//':concentration_mol_m3_2_mmr'
    concentration_mol_m3_2_mmr = GAS_conc / (molarMassDryAir/GAS%molarMass_kg) * (argas * T_ref / p_sealevel)
  END FUNCTION concentration_mol_m3_2_mmr

  !!
  ! @brief Converts a gas concentration to volume mixing ratio at reference temperature and pressure
  !
  !!
#ifndef _OPENACC
  PURE &
#endif
  REAL(wp) FUNCTION concentration_mmr_2_mol_m3(GAS, GAS_mmr)
  !$acc routine seq
    TYPE (t_specie_const), INTENT(in) :: GAS     ! Type with constants specific for the gas
    REAL(wp),              INTENT(in) :: GAS_mmr ! Mass mixing ratio of gas
    CHARACTER(len=*), PARAMETER :: routine = modname//':concentration_mmr_2_mol_m3'
    ! from ideal gas law: PV = nRT <=> V/n = RT/P
    ! volume mixing ratio is molar mixing ratio (1 mol gas = constant volume, independent of gas)
    ! => determine # of moles in volume
    concentration_mmr_2_mol_m3 = GAS_mmr * (molarMassDryAir/GAS%molarMass_kg) / (argas * T_ref / p_sealevel)
  END FUNCTION concentration_mmr_2_mol_m3

  !!
  ! @brief Converts a gas volume mixing ratio to mass mixing ratio at reference temperature and pressure
  !
  !!
#ifndef _OPENACC
  PURE &
#endif
  REAL(wp) FUNCTION concentration_vmr_2_mmr(GAS, GAS_vmr)
  !$acc routine seq
    TYPE (t_specie_const), INTENT(in) :: GAS     ! Type with constants specific for the gas
    REAL(wp),              INTENT(in) :: GAS_vmr ! Volume mixing ratio of gas
    CHARACTER(len=*), PARAMETER :: routine = modname//':concentration_vmr_2_mmr'
    concentration_vmr_2_mmr = concentration_mol_m3_2_mmr(GAS,concentration_vmr_2_mol_m3(GAS,GAS_vmr))
  END FUNCTION concentration_vmr_2_mmr

  !!
  ! @brief Converts a gas mass mixing ratio to volume mixing ratio at reference temperature and pressure
  !
  !!
#ifndef _OPENACC
  PURE &
#endif
  REAL(wp) FUNCTION concentration_mmr_2_vmr(GAS, GAS_mmr)
  !$acc routine seq
    TYPE (t_specie_const), INTENT(in) :: GAS     ! Type with constants specific for the gas
    REAL(wp),              INTENT(in) :: GAS_mmr ! Mass mixing ratio of gas [kg kg-1]
    CHARACTER(len=*), PARAMETER :: routine = modname//':concentration_mmr_2_vmr'
    concentration_mmr_2_vmr = concentration_mol_m3_2_vmr(GAS,concentration_mmr_2_mol_m3(GAS,GAS_mmr))
  END FUNCTION concentration_mmr_2_vmr

  !!
  ! @brief Determine which unit of a gas is properly defined
  !
  !!
#ifndef _OPENACC
  PURE &
#endif
  INTEGER FUNCTION which_valid(GAS,mmr,vmr,conc)
  !$acc routine seq
    TYPE (t_specie_const), INTENT(in)  :: GAS
    REAL(wp), OPTIONAL,    INTENT(in)  :: mmr
    REAL(wp), OPTIONAL,    INTENT(in)  :: vmr
    REAL(wp), OPTIONAL,    INTENT(in)  :: conc

    CHARACTER(len=*), PARAMETER :: routine = modname//':which_valid'
    INTEGER n, which

    ! Check presence of exactly 1 unit
    n = 0
    which = 0
    IF (PRESENT(mmr)) THEN
      IF (mmr  >= 0._wp) THEN
        n = n + 1
        which = 1
      ENDIF
    ENDIF
    IF (PRESENT(vmr)) THEN
      IF (vmr  >= 0._wp) THEN
        n = n + 1
        which = 2
      ENDIF
    ENDIF
    IF (PRESENT(conc)) THEN
      IF (conc >= 0._wp) THEN
        n = n + 1
        which = 3
      ENDIF
    ENDIF

    which_valid = which

  END FUNCTION which_valid

  !!
  ! @brief Converts gas concentration or volume mixing ratio to mass mixing ratio
  ! Invalid units are marked with negative values or not present
  ! Exactly 1 unit must be specified
#ifndef _OPENACC
  PURE &
#endif
  SUBROUTINE convert_2_mmr(GAS,mmr,vmr,conc)
  !$acc routine seq
    TYPE (t_specie_const),  INTENT(in)    :: GAS
    REAL(wp),     OPTIONAL, INTENT(inout) :: mmr
    REAL(wp),     OPTIONAL, INTENT(in)    :: vmr
    REAL(wp),     OPTIONAL, INTENT(in)    :: conc

    CHARACTER(len=*), PARAMETER :: routine = modname//':convert_2_mmr'
    INTEGER :: which

    IF (.NOT. PRESENT(mmr)) RETURN
    which = which_valid(GAS,mmr,vmr,conc)

    SELECT CASE (which)
    CASE (1) ! mmr given => nothing to do
    CASE (2) ! vmr given
      mmr = concentration_vmr_2_mmr(GAS,vmr)
    CASE (3) ! concentration given
      mmr = concentration_mol_m3_2_mmr(GAS,conc)
    END SELECT

  END SUBROUTINE convert_2_mmr

  !!
  ! @brief Converts gas concentration or mass mixing ratio to volume mixing ratio
  ! Invalid units are marked with negative values or not present
  ! Exactly 1 unit must be specified
#ifndef _OPENACC
  PURE &
#endif
  SUBROUTINE convert_2_vmr(GAS,mmr,vmr,conc)
  !$acc routine seq
    TYPE (t_specie_const), INTENT(in)    :: GAS
    REAL(wp), OPTIONAL,    INTENT(in)    :: mmr
    REAL(wp), OPTIONAL,    INTENT(inout) :: vmr
    REAL(wp), OPTIONAL,    INTENT(in)    :: conc

    CHARACTER(len=*), PARAMETER :: routine = modname//':convert_2_vmr'
    INTEGER :: which

    IF (.NOT. PRESENT(vmr)) RETURN
    which = which_valid(GAS,mmr,vmr,conc)

    SELECT CASE (which)
    CASE (1) ! mmr given
      vmr = concentration_mmr_2_vmr(GAS,mmr)
    CASE (2) ! vmr given => nothing to do
    CASE (3) ! concentration given
      vmr = concentration_mol_m3_2_vmr(GAS,conc)
    END SELECT

  END SUBROUTINE convert_2_vmr

  !!
  ! @brief Converts volume or mass mixing ratio to concentration
  ! Invalid units are marked with negative values or not present
  ! Exactly 1 unit must be specified
#ifndef _OPENACC
  PURE &
#endif
  SUBROUTINE convert_2_mol_m3(GAS,mmr,vmr,conc)
  !$acc routine seq
    TYPE (t_specie_const), INTENT(in)    :: GAS
    REAL(wp), OPTIONAL,    INTENT(in)    :: mmr
    REAL(wp), OPTIONAL,    INTENT(in)    :: vmr
    REAL(wp), OPTIONAL,    INTENT(inout) :: conc

    CHARACTER(len=*), PARAMETER :: routine = modname//':convert_2_mol_m3'
    INTEGER :: which

    IF (.NOT. PRESENT(conc)) RETURN
    which = which_valid(GAS,mmr,vmr,conc)

    SELECT CASE (which)
    CASE (1) ! mmr given
      conc = concentration_mmr_2_mol_m3(GAS,mmr)
    CASE (2) ! vmr given
      conc = concentration_vmr_2_mol_m3(GAS,vmr)
    CASE (3) ! concentration given => nothing to do
    END SELECT

  END SUBROUTINE convert_2_mol_m3

  !!
  ! @brief Calculates diffusivities for a gas in air and water
  !
  !!
#ifndef _OPENACC
  PURE &
#endif
  SUBROUTINE calc_diffusivity(GAS,TK,D_air,D_water)
  !$acc routine seq

    USE mo_jsb_physical_constants, ONLY: tmelt
    USE mo_methane_constants,      ONLY: t_specie_const

    TYPE (t_specie_const), INTENT(in)  :: GAS     ! Constants for relevante GAS
    REAL(wp),              INTENT(in)  :: TK      ! Temperature to calculate diffusivities for [K]
    REAL(wp), OPTIONAL,    INTENT(out) :: D_air   ! Diffusivity of gas in air [m2 s-1]
    REAL(wp), OPTIONAL,    INTENT(out) :: D_water ! Diffusivity of gas in air [m2 s-1]

    REAL(wp) :: T

    T = TK - tmelt ! Convert to degC

    ! Constants are in [cm2 s-1]. convert to [m2 s-1] on the fly
    IF (PRESENT(D_air)) &
      D_air   = (GAS%D_air_0   + GAS%D_air_1   * T)                                             * 1.e-4_wp
    IF (PRESENT(D_water)) &
      D_water = (GAS%D_water_0 + GAS%D_water_1 * T + GAS%D_water_2 * T * T) * GAS%D_water_scale * 1.e-4_wp

  END SUBROUTINE calc_diffusivity

  !!
  ! @brief Adds a surface layer to a layer vector. Vectors must have the right size upon call
  !
  !!
  ! If surf not specified, the "upper layer" of src is copied to the new upper layer
#ifndef _OPENACC
  PURE &
#endif
  SUBROUTINE add_surface_layer(src,dst,surf)
  !$acc routine seq
    REAL(wp),           INTENT(in)  :: src(:)
    REAL(wp),           INTENT(out) :: dst(:)
    REAL(wp), OPTIONAL, INTENT(in)  :: surf

    CHARACTER(len=*), PARAMETER :: routine = modname//':add_surface_layer'
    INTEGER :: is

    DO is=2,SIZE(dst)
      dst(is) = src(is-1)
    ENDDO
    IF (PRESENT(surf)) THEN
      dst(1) = surf
    ELSE
      dst(1) = src(1)
    ENDIF

  END SUBROUTINE add_surface_layer

#ifndef _OPENACC
  !!
  ! @brief   Setup and calc vertical diffusion of a gas using a Crank-Nicolson integration
  !
  !!
  SUBROUTINE diffusion_cn(nstep, dt, mxlay_above_rock, mxlay_above_wt,                   &
                          implicit_contrib, surf_GAS_conc, GAS_source, kHCC, diff_coeff, &
                          diff_coeff_surf, layer_thick, dens2conc, conc2dens,            &
                          GAS_dens, surf_GAS_flux, lissue)

  !$acc routine seq
    INTEGER,  INTENT(in)    :: nstep              ! #diffusion timesteps per model timestep
    REAL(wp), INTENT(in)    :: dt                 ! Model time step [s]
    INTEGER,  INTENT(in)    :: mxlay_above_rock   ! layer number above bedrock
    INTEGER,  INTENT(in)    :: mxlay_above_wt     ! layer number above water table
    REAL(wp), INTENT(in)    :: implicit_contrib   ! Implicit contribution for Crank-Nicholson
    REAL(wp), INTENT(in)    :: surf_GAS_conc      ! GAS concentration at surface [mol m-3]
    REAL(wp), INTENT(in)    :: GAS_source(:)      ! GAS source / sink in layer [mol m-2 ts-1]
    REAL(wp), INTENT(in)    :: kHCC(:)            ! Henry coefficient relating aqueous to gaseous concentration
    REAL(wp), INTENT(in)    :: diff_coeff(:)      ! Layer diffusion coefficients  [m2 s-1]
    REAL(wp), INTENT(in)    :: diff_coeff_surf    ! surface diffusion coefficient [m2 s-1]
    REAL(wp), INTENT(in)    :: layer_thick(:)     ! Layer thickness [m]
    REAL(wp), INTENT(in)    :: dens2conc(:)       ! conversion factor from density to concentration []
    REAL(wp), INTENT(in)    :: conc2dens(:)       ! conversion factor from concentration to density []
    REAL(wp), INTENT(inout) :: GAS_dens(:)        ! GAS density in layer [mol m-2]
    REAL(wp), INTENT(out)   :: surf_GAS_flux      ! GAS flux through surface [mol m-2 ts-1]
    LOGICAL,  INTENT(out), OPTIONAL :: lissue     ! Error message issued from matrix solver

    ! locals
    CHARACTER(len=*), PARAMETER :: routine = modname//':diffusion_cn'

    ! For diffusion calculation, we use a different layer counting scheme:
    ! We add a surface layer at atmospheric CH4 concentration, therefore
    ! all soil layers are moved down by one
    ! surface layer has same thickness as first belowsurface layer

    ! layer counters
    INTEGER  :: outside_layer, iLayer, num_iLayers
    INTEGER  :: iStep
    REAL(wp) :: nStep_fac
    LOGICAL  :: lissue_tst, lissue_loc

    ! local vectors of concentrations, cofficients, depths -- one element longer to consider surface state as well
    REAL(wp) :: GAS_conc        (mxlay_above_rock+1)            ! GAS concentration [mol m-3]
    REAL(wp) :: src_GAS_conc_dts(mxlay_above_rock+1)            ! source GAS concentration [mol m-3]
    REAL(wp) :: kHCC_loc        (mxlay_above_rock+1)            ! Henry coefficient [-]
    REAL(wp) :: diff_coeff_loc  (mxlay_above_rock+1)            ! diffusivity [m2 s-1]
    REAL(wp) :: dz              (mxlay_above_rock+1)            ! layer thickness [m]

    ! bookkeeping (all [mol m-2])
    REAL(wp) :: tot_dens_initial      ! total GAS density in soil column before diffusion
    REAL(wp) :: tot_dens_final        ! total GAS density in soil column after diffusion
    REAL(wp) :: tot_dens_source       ! total GAS density in soil column from sources and sinks
    REAL(wp) :: delta_time

    nStep_fac = REAL(nstep,wp)
    ! initialisation
    num_iLayers = mxlay_above_rock + 1

    ! concentration / density conversions
    CALL add_surface_layer(GAS_dens(:)   * dens2conc(:),     GAS_conc,     surf_GAS_conc)
    CALL add_surface_layer(GAS_source(:) * dens2conc(:), src_GAS_conc_dts, 0._wp)

    ! totals for balance checks
    tot_dens_initial = SUM(GAS_dens  (1:mxlay_above_rock))
    tot_dens_source  = SUM(GAS_source(1:mxlay_above_rock))

    CALL add_surface_layer(kHCC,kHCC_loc)
    CALL add_surface_layer(layer_thick,dz)
    CALL add_surface_layer(diff_coeff,diff_coeff_loc,diff_coeff_surf)

    delta_time       = dt               / nStep_fac
    src_GAS_conc_dts = src_GAS_conc_dts / nStep_fac

    lissue_loc = .FALSE.
    DO iStep = 1, nStep
      CALL diffusion_step_variable(num_iLayers, mxlay_above_wt, implicit_contrib, src_GAS_conc_dts, &
                                   kHCC_loc, diff_coeff_loc, dz, delta_time, GAS_conc, lissue_tst)
      lissue_loc = lissue_loc .OR. lissue_tst
    END DO

    ! final balance and output flux
    ! negative densities might be possible --> reset to zero, difference will be added to surface flux
    GAS_dens(1:mxlay_above_rock) = MAX(0._wp,GAS_conc(2:mxlay_above_rock+1) * conc2dens(1:mxlay_above_rock))

    tot_dens_final = SUM(GAS_dens(1:mxlay_above_rock))
    surf_GAS_flux = tot_dens_final - (tot_dens_initial + tot_dens_source)

    IF (PRESENT(lissue)) lissue = lissue_loc

  END SUBROUTINE diffusion_cn

  !!
  ! @brief   Calc vertical diffusion of a gas
  !
  !!
  !-----------------------------------------------------------------------------
  ! diffusion_step_variable
  ! should allow variable set up of implicit / explicit
  ! contributions to diffision problem
  !

  SUBROUTINE diffusion_step_variable(num_iLayers, mxlay_above_wt, implicit_contrib, src_GAS_conc, &
                                     kHCC, diff_coeff, dz, delta_time, GAS_conc, lissue)
  !$acc routine seq
    INTEGER , INTENT(in)    :: num_iLayers      ! no. of internal layers
    INTEGER,  INTENT(in)    :: mxlay_above_wt   ! (outside) layer number above water table
    REAL(wp), INTENT(in)    :: implicit_contrib ! implicit contribution to Crank-Nicolson integration
    REAL(wp), INTENT(in)    :: src_GAS_conc(:)  ! source GAS concentration in layer [mol m-3]
    REAL(wp), INTENT(in)    :: kHCC(:)          ! dimensionless Henry coefficient relating aqueous to gaseous concentration
    REAL(wp), INTENT(in)    :: diff_coeff(:)    ! Layer diffusion coefficients [m2 s-1]
    REAL(wp), INTENT(in)    :: dz(:)            ! Layer thickness [m]
    REAL(wp), INTENT(in)    :: delta_time       ! time step length [s]
    REAL(wp), INTENT(inout) :: GAS_conc(:)      ! GAS concentration in layer [mol m-3]
    LOGICAL, INTENT(out), OPTIONAL :: lissue ! Error message from matrix solver issued

    ! locals

    CHARACTER(len=*), PARAMETER :: routine = modname//':diffusion_step_variable'

    ! For diffusion calculation, we use a different layer counting scheme:
    ! We add a surface layer at atmospheric CH4 concentration, therefore
    ! all soil layers are moved down by one
    ! surface layer has same thickness as first belowsurface layer

    ! layer counters
    INTEGER :: outside_layer, iLayer

    ! Diffusion constants averaged over adjacent layers.
    REAL(wp) :: D_below, D_above

    ! Layer-dependent weighting factors in Crank-Nicholson scheme.
    REAL(wp) :: dz_factor, W_below, W_above
    REAL(wp) :: dz_here, dz_above, dz_below
    REAL(wp) :: C_here, C_above, C_below

    ! diagonal, left and right subdiagonals for Crank-Nicholson matrix
    REAL(wp), DIMENSION(num_iLayers) :: left
    REAL(wp), DIMENSION(num_iLayers) :: diag
    REAL(wp), DIMENSION(num_iLayers) :: right

    ! Right hand side vector for Crank-Nicholson scheme equations
    ! This is overwritten with solution by dgtsv tridiagonal solver
    REAL(wp), DIMENSION(num_iLayers) :: rhs

    ! bookkeeping
    REAL(wp)                                    :: correction       ! correction from inaccuracies of CN-scheme [mol m-2]
    INTEGER                                     :: info             ! dgtsv status variable
    LOGICAL                                     :: consider_aqueous ! switch to consider aqueous concentrations at interface

    consider_aqueous = .TRUE.

    !------------------------------------------------------------
    ! build tridiagonal matrix and known right hand side

    ! Process the local layers
    DO iLayer = 1, num_iLayers
      ! outside layers are shifted by one
      outside_layer = iLayer-1

      ! Extract sensible values to use for concentration and
      ! thickness of the current layer, the layer above and below
      dz_here = dz(iLayer)
      C_here = GAS_conc(iLayer)
      IF (iLayer == 1) THEN
        dz_above = 0._wp
        C_above = 0._wp
      ELSE
        dz_above = dz(iLayer-1)
        C_above = GAS_conc(iLayer-1)
        ! just below water table: use aqueous concentration of gas
        IF ((outside_layer == mxlay_above_wt+1) .AND. consider_aqueous) THEN
          C_above = C_above*kHCC(iLayer-1)
        END IF
      END IF
      IF (iLayer == num_iLayers) THEN
        dz_below = 0._wp
        C_below = 0._wp
      ELSE
        dz_below = dz(iLayer+1)
        C_below = GAS_conc(iLayer+1)
      END IF

      ! Determine diffusion coefficients to the layers above and below the current layer.
      IF (iLayer == 1) THEN
        ! upwards diffusivity in uppermost layer is zero
        D_above = 0._wp
      ELSE
        ! for the layers below, upwards diffusivity is the weighted mean of the layer
        ! diffusivity in layer and layer above
        D_above = (dz(iLayer-1)*diff_coeff(iLayer-1) + &
                   dz(iLayer  )*diff_coeff(iLayer))  / &
                  (dz(iLayer-1) + dz(iLayer))
        IF ((outside_layer == mxlay_above_wt+1) .AND. consider_aqueous) THEN
          D_above = (dz(iLayer-1) * (diff_coeff(iLayer-1)  / kHCC(iLayer-1)) + &
                     dz(iLayer  ) *  diff_coeff(iLayer  )) / (dz(iLayer-1) + dz(iLayer))
        END IF
      END IF
      IF (iLayer == num_iLayers) THEN
        ! downwards diffusivity in the bottom layer is set to zero
        ! to impose no flux boundary condition
        D_below = 0._wp
      ELSE
        ! for the layers above, downwards diffusivity is the weighted mean of the
        ! diffusivity in the layer and the layer below
        D_below = (dz(iLayer  )*diff_coeff(iLayer  )  + &
                   dz(iLayer+1)*diff_coeff(iLayer+1)) / &
                  (dz(iLayer) + dz(iLayer+1))
        IF ((outside_layer == mxlay_above_wt) .AND. consider_aqueous) THEN
          D_below = (dz(iLayer  ) * (diff_coeff(iLayer  )  / kHCC(iLayer)) + &
                     dz(iLayer+1) *  diff_coeff(iLayer+1)) / (dz(iLayer) + dz(iLayer+1))
        END IF
      END IF

      ! Calculate layer-dependent weighting factors.
      dz_factor = 0.25_wp * (dz_below + 2._wp * dz_here + dz_above)
      W_below = D_below * (delta_time / dz_factor) * (2._wp / (dz_below + dz_here))
      W_above = D_above * (delta_time / dz_factor) * (2._wp / (dz_here + dz_above))

      ! boundary conditions for surface and bottom layers: no flux through top / bottom
      IF (iLayer == 1) THEN
        W_above = 0._wp
      END IF
      IF (iLayer == num_iLayers) THEN
        W_below = 0._wp
      END IF

      ! Fill in matrix diagonal and off-diagonal elements.
      diag(iLayer) = 1.0_wp + (implicit_contrib * (W_below + W_above))
      IF ((outside_layer == mxlay_above_wt) .AND. consider_aqueous) THEN
        diag(iLayer) = 1.0_wp + (implicit_contrib * (W_below*kHCC(iLayer) + W_above))
      END IF
      left(iLayer) = -implicit_contrib * W_above
      IF ((outside_layer == mxlay_above_wt + 1) .AND. consider_aqueous) THEN
        left(iLayer) = -implicit_contrib * W_above*kHCC(iLayer-1)
      END IF
      right(iLayer) = -implicit_contrib * W_below

      ! Calculate right hand side vector values.
      rhs(iLayer) = (1.0_wp - (1.0_wp - implicit_contrib)*(W_below + W_above)) * C_here + &
                    (1.0_wp - implicit_contrib)*W_below*C_below +                         &
                    (1.0_wp - implicit_contrib)*W_above*C_above +                         &
                     src_GAS_conc(iLayer)
      IF ((outside_layer == mxlay_above_wt) .AND. consider_aqueous) THEN
        rhs(iLayer) = (1.0_wp - (1.0_wp - implicit_contrib)*((W_below*kHCC(iLayer)) + W_above))*C_here + &
                      (1.0_wp - implicit_contrib)*W_below*C_below +                       &
                      (1.0_wp - implicit_contrib)*W_above*C_above +                       &
                       src_GAS_conc(iLayer)
      END IF
    END DO

    ! Semi-implicit case -- whenever implicit_contrib > 0.1,
    ! fully explicit case is solced already ( = rhs)
    IF (implicit_contrib >= 0.1_wp) THEN
      !------------------------------------------------------------
      ! solve tridiagonal system using dgtsv from NETLIB
      CALL dgtsv(num_iLayers, 1, left(2:num_iLayers), diag, &
                 right(1:num_iLayers-1), rhs, num_iLayers, info)
      IF (PRESENT(lissue)) lissue = info /= 0
!      IF (info /= 0) THEN
!        WRITE(message_text,*) 'dgtsv info ', info
!        CALL message(routine, message_text, all_print=.TRUE., level=em_warn)
!      END IF
    END IF

    ! solution is ALWAYS written into the rhs vector, also in the fully explicit case
    GAS_conc = rhs

  END SUBROUTINE diffusion_step_variable
#endif
#endif
END MODULE mo_methane_util
