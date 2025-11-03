!>
!! @brief Contains the routines for the methane production and transport
!!
!! @author
!!  Stiig Wilkenskjeld, MPI-M Hamburg
!!
!! @par Revision History
!! Adopted from JSBACH3 code by Thomas Kleinen, 2023-03-21
!!
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/license .
!!

MODULE mo_methane_process
#ifndef __NO_JSBACH__

  USE mo_kind,     ONLY: wp

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: soil_column_setup, diffusion_setup, production, oxidation_mm, ebullition
  PUBLIC :: plant_setup, plant_transport, init_gas_pools, reaction_diffusion, determine_root_dist
  PUBLIC :: calc_methane_fire, calc_methane_termites, calc_methane_herbivores

  CHARACTER(len=*), PARAMETER :: modname = 'mo_methane_process'

CONTAINS

  !!
  ! @brief     determine soil column setup
  !
  !!
  PURE SUBROUTINE soil_column_setup(permafrost_flag,is_peatland,is_wet, mxlay_above_rock,                      &
                                    layer_thick, layer_mid, layer_bound, soil_depth, soil_porosity,            &
                                    water_table, acrotelm_height,                                              &
                                    layer_moisture, ice_content,                                               & 
                                    relative_moist, relative_ice, relative_porevol, relative_pore_air,         &
                                    relative_pore_moist, dens2conc, conc2dens, mxlay_above_wt, mxlay_acrotelm)
    USE mo_methane_constants, ONLY: minimum_porevol
    USE mo_methane_util,      ONLY: find_layer
    LOGICAL,  INTENT(in)  :: permafrost_flag        ! permafrost flag
    LOGICAL,  INTENT(in)  :: is_peatland            ! peatland flag
    LOGICAL,  INTENT(in)  :: is_wet                 ! indicate dry / wet part of tile: .FALSE. = dry, .TRUE. = wet
    REAL(wp), INTENT(in)  :: layer_thick(:)         ! soil layer thickness [m]
    REAL(wp), INTENT(in)  :: layer_mid(:)           ! soil layer mid depth [m]
    REAL(wp), INTENT(in)  :: layer_bound(:)         ! soil layer boundary [m] TODO: Upper/lower?
    REAL(wp), INTENT(in)  :: soil_depth             ! grid cell soil depth [m]
    REAL(wp), INTENT(in)  :: soil_porosity          ! grid cell soil porosity []
    INTEGER,  INTENT(in)  :: mxlay_above_rock       ! #lowest layer above bedrock
    REAL(wp), INTENT(in)  :: water_table            ! water table position [m]
    REAL(wp), INTENT(in)  :: acrotelm_height        ! height of acrotelm [m]
    REAL(wp), INTENT(in)  :: layer_moisture(:)      ! soil layer moisture content [?]
    REAL(wp), INTENT(in)  :: ice_content(:)         ! soil layer ice content [?]
    REAL(wp), INTENT(out) :: relative_moist(:)      ! soil layer relative moisture content []
    REAL(wp), INTENT(out) :: relative_ice(:)        ! soil layer relative ice content []
    REAL(wp), INTENT(out) :: relative_porevol(:)    ! relative pore volume in layer []
    REAL(wp), INTENT(out) :: relative_pore_air(:)   ! relative pore air content in layer []
    REAL(wp), INTENT(out) :: relative_pore_moist(:) ! relative pore moisture in layer []
    REAL(wp), INTENT(out) :: dens2conc(:)           ! conversion factor from density to concentration [m-1]
    REAL(wp), INTENT(out) :: conc2dens(:)           ! conversion factor from concentration to density [m]
    INTEGER,  INTENT(out) :: mxlay_above_wt         ! layer number above water table
    INTEGER,  INTENT(out) :: mxlay_acrotelm         ! layer number for acrotelm end

    CHARACTER(len=*), PARAMETER :: routine = modname//':soil_colum_setup'
    REAL(wp),         PARAMETER :: relative_moist_saturated = 0.95_wp
    REAL(wp),         PARAMETER :: relative_porevol_min     = minimum_porevol/2._wp

    ! locals
    REAL(wp), ALLOCATABLE :: effective_layer_volume(:)
    REAL(wp) :: current_thickness
    INTEGER  :: nsoil, iLayer

    relative_moist         = 0._wp
    relative_ice           = 0._wp
    relative_porevol       = 0._wp
    relative_pore_moist    = 0._wp
    relative_pore_air      = 0._wp
    effective_layer_volume = 0._wp
    conc2dens              = 0._wp
    dens2conc              = 0._wp

    nsoil = SIZE(layer_thick)
    ALLOCATE(effective_layer_volume(nsoil))

    !  determine relative moisture and ice content in soil layers ->  divide by layer thickness
    DO iLayer = 1,mxlay_above_rock
      current_thickness = layer_thick(iLayer)
      IF ((iLayer > 1) .AND. (layer_bound(iLayer) > soil_depth)) THEN
        current_thickness = soil_depth - layer_bound(iLayer-1)
      END IF
      relative_moist(iLayer) = layer_moisture(iLayer) / current_thickness
      relative_ice  (iLayer) = ice_content(iLayer)    / current_thickness
    END DO
    relative_moist(mxlay_above_rock+1:nsoil) = 0._wp
    relative_ice  (mxlay_above_rock+1:nsoil) = 1._wp

    IF (.NOT. permafrost_flag) THEN
      relative_ice = 0._wp
    END IF

    IF (is_wet) THEN
      ! determine, in which layer water table is (remember: water table below surface is negative)
      mxlay_above_wt = find_layer(-water_table, mxlay_above_rock, layer_mid)

      ! determine, in which layer acrotelm ends
      mxlay_acrotelm = find_layer(acrotelm_height, mxlay_above_rock, layer_mid)
    ELSE
      mxlay_above_wt = mxlay_above_rock
      mxlay_acrotelm = 0
    END IF

    ! determine relative pore volume, corrected by ice
    relative_porevol(1:mxlay_above_rock) = MAX(relative_porevol_min, soil_porosity - relative_ice(1:mxlay_above_rock))

    ! calculate the relative moisture content of ice free pores
    relative_pore_moist(1:mxlay_above_rock) = MIN(relative_moist_saturated, &
                                                  relative_moist(1:mxlay_above_rock) / relative_porevol(1:mxlay_above_rock))

    ! correct relative moisture and relative moisture content for layers below 
    ! water table -> set to relative_moist_saturated 
    relative_pore_moist((mxlay_above_wt+1):mxlay_above_rock) = relative_moist_saturated
    relative_moist     ((mxlay_above_wt+1):mxlay_above_rock) = relative_pore_moist((mxlay_above_wt+1):mxlay_above_rock) * &
                                                               relative_porevol   ((mxlay_above_wt+1):mxlay_above_rock)

    relative_pore_air(1:mxlay_above_rock) = 1._wp - relative_pore_moist(1:mxlay_above_rock)

    ! helper for transferring between gas densities and concentrations
    effective_layer_volume(1:mxlay_above_rock) = layer_thick(1:mxlay_above_rock) * relative_porevol(1:mxlay_above_rock)

    conc2dens(1:mxlay_above_rock) =         effective_layer_volume(1:mxlay_above_rock)
    dens2conc(1:mxlay_above_rock) = 1._wp / effective_layer_volume(1:mxlay_above_rock)

    DEALLOCATE(effective_layer_volume)

  END SUBROUTINE soil_column_setup

  !!
  ! @brief     Calculate the diffusivities for each gas in both water
  !            and air once, to be allocated as needed
  !
  !!
  PURE SUBROUTINE diffusion_setup(is_peatland, mxlay_above_rock, mxlay_above_wt, bclapp, &
                                  relative_porevol, relative_pore_air, Tsoil, &
                                  D_CH4, D_CO2, D_O2, D_CH4_surface, D_CO2_surface, D_O2_surface)
    USE mo_methane_constants, ONLY: CO2_consts, CH4_consts, O2_consts
    USE mo_methane_util,      ONLY: calc_diffusivity

    LOGICAL,  INTENT(in)  :: is_peatland          ! peatland flag (TRUE = peatland)
    INTEGER,  INTENT(in)  :: mxlay_above_rock     ! layer number above bedrock
    INTEGER,  INTENT(in)  :: mxlay_above_wt       ! layer number above water table
    REAL(wp), INTENT(in)  :: bclapp               ! exponent b in Clapp and Hornberger
    REAL(wp), INTENT(in)  :: relative_pore_air(:) ! relative pore air content in layer []
    REAL(wp), INTENT(in)  :: relative_porevol(:)  ! relative pore volume in layer []
    REAL(wp), INTENT(in)  :: Tsoil(:)             ! Soil temperature [K]
    REAL(wp), INTENT(out) :: D_CH4(:)             ! diffusivity of CH4 on layers [m2 s-1]
    REAL(wp), INTENT(out) :: D_CO2(:)             ! diffusivity of CO2 on layers [m2 s-1]
    REAL(wp), INTENT(out) :: D_O2(:)              ! diffusivity of O2 on layers [m2 s-1]
    REAL(wp), INTENT(out) :: D_CH4_surface        ! diffusivity of CH4 in surface air [m2 s-1]
    REAL(wp), INTENT(out) :: D_CO2_surface        ! diffusivity of CO2 in surface air [m2 s-1]
    REAL(wp), INTENT(out) :: D_O2_surface         ! diffusivity of O2 in surface air [m2 s-1]

    CHARACTER(len=*), PARAMETER :: routine = modname//':diffusion_setup'

    ! locals
    REAL(wp) :: D_CH4_air  (mxlay_above_rock)        ! diffusivity of CH4 in air [m2 s-1]
    REAL(wp) :: D_CO2_air  (mxlay_above_rock)        ! diffusivity of CO2 in air [m2 s-1]
    REAL(wp) :: D_O2_air   (mxlay_above_rock)        ! diffusivity of O2 in air [m2 s-1]
    REAL(wp) :: D_CH4_water(mxlay_above_rock)        ! diffusivity of CH4 in water [m2 s-1]
    REAL(wp) :: D_CO2_water(mxlay_above_rock)        ! diffusivity of CO2 in water [m2 s-1]
    REAL(wp) :: D_O2_water (mxlay_above_rock)        ! diffusivity of O2 in water [m2 s-1]

    REAL(wp) :: conv                                            ! conversion factor
    INTEGER  :: is                                              ! layer counter

    DO is = 1, mxlay_above_rock
      CALL calc_diffusivity(CH4_consts,Tsoil(is),D_CH4_air(is),D_CH4_water(is))
      CALL calc_diffusivity(CO2_consts,Tsoil(is),D_CO2_air(is),D_CO2_water(is))
      CALL calc_diffusivity( O2_consts,Tsoil(is),D_O2_air (is),D_O2_water (is))
    ENDDO

    D_CH4 = 0._wp
    D_CO2 = 0._wp
    D_O2  = 0._wp

    ! Use Millington and Quirk 1961 (MQ) Model tested in
    ! Iiyama and Hasegawa 2005.
    ! When Fair is less then 0.05, then use the diffusivity
    ! of water.  The cut-off of 0.05 is where Dwater is
    ! about equal D using the MQ model.

    ! Original Wania code assigns D_GAS_air, if relative_pore_air(iLayer) > 0.05
    ! We assume that this is the case ABOVE the water table

    DO is = 1, mxlay_above_wt
      IF (.NOT. is_peatland) THEN
        ! Moldrup et al. diffusivity modification for soil pore air in mineral soils
        conv = ((relative_pore_air(is) * relative_porevol(is)) ** 2._wp) * (relative_pore_air(is) ** (3._wp/bclapp))
      ELSE
        ! Millington & Quirk diffusivity modification for soil pore air in organic soils
        conv = ((relative_pore_air(is) * relative_porevol(is)) ** (10._wp/3._wp)) / (relative_porevol(is) ** 2._wp)
      END IF
      D_CH4(is) = D_CH4_air(is) * conv
      D_CO2(is) = D_CO2_air(is) * conv
      D_O2 (is) = D_O2_air (is) * conv
    END DO
    DO is = (mxlay_above_wt + 1), mxlay_above_rock
      ! Moldrup et al. diffusivity modification for soil pore water
      conv = relative_porevol(is) ** 2._wp
      D_CH4(is) = D_CH4_water(is) * conv
      D_CO2(is) = D_CO2_water(is) * conv
      D_O2 (is) = D_O2_water (is) * conv
    END DO

    ! Ensure that layer diffusivity is not larger than for layer(s) above
    DO is = 2, mxlay_above_rock
      IF (D_CH4(is) > D_CH4(is-1)) &
          D_CH4(is) = D_CH4(is-1)
      IF (D_CO2(is) > D_CO2(is-1)) &
          D_CO2(is) = D_CO2(is-1)
      IF (D_O2 (is) > D_O2 (is-1)) &
          D_O2 (is) = D_O2 (is-1)
    END DO

    ! Diffusivity to the air is assumed to be the same as for the upper soil layer
    D_CH4_surface = D_CH4_air(1)
    D_CO2_surface = D_CO2_air(1)
    D_O2_surface  = D_O2_air (1)

  END SUBROUTINE diffusion_setup

  !!
  ! @brief     set up plant transport parameters for all plant transports on tile
  !            determines tiller_area and tiller_surface
  !
  !!
  PURE SUBROUTINE plant_setup(is_wet,is_woody,LAI_max,mxlay_above_rock,LAI,Cpool_green,root_dist, &
                              tiller_area,tillers)
    USE mo_jsb_math_constants, ONLY: pi
    LOGICAL,            INTENT(in)  :: is_wet           ! wetness state: .FALSE. = dry, .TRUE. = wet
    LOGICAL,            INTENT(in)  :: is_woody         ! woody landcover type
    INTEGER,            INTENT(in)  :: mxlay_above_rock ! layer number above bedrock
    REAL(wp),           INTENT(in)  :: LAI_max          ! maximum leaf area index []
    REAL(wp),           INTENT(in)  :: LAI              ! actual, leaf area index []
    REAL(wp),           INTENT(in)  :: Cpool_green      ! C-pool for leaves, fine roots, vegetative organs and other living parts 
                                                        ! of vegetation [mol m-2(canopy)]
    REAL(wp),           INTENT(in)  :: root_dist(:)     ! root distribution in each layer
    REAL(wp),           INTENT(out) :: tiller_area(:)   ! area of tillers corresponding to the root fraction in each layer [m2?]
    REAL(wp), OPTIONAL, INTENT(out) :: tillers          ! number of tillers

    ! Parameters
    CHARACTER(len=*), PARAMETER :: routine = modname//':plant_setup'
    REAL(wp), PARAMETER   :: tiller_weight    = 0.0183_wp ! Individual tiller weight, Wania original is 0.22 gC -> conv to [mol]
    REAL(wp), PARAMETER   :: tiller_radius    = 0.003_wp  ! Radius of an average tiller [m]
    REAL(wp), PARAMETER   :: tiller_por_grass = 0.7_wp    ! Tiller porosity grass (Wania)
    REAL(wp), PARAMETER   :: tiller_por_tree  = 0.3_wp    ! Tiller porosity tree (arbitrary)
    REAL(wp), PARAMETER   :: tiller_por_min   = 0.05_wp   ! Tiller porosity min (Riley)
    REAL(wp), PARAMETER   :: dry_ratio        = 0.3_wp    ! Tiller porosity tree (arbitrary)

    ! locals
    REAL(wp)              :: n                            ! #tillers
    REAL(wp)              :: tiller_frac                  ! fraction of tillers weighted by root biomass in each layer
    REAL(wp)              :: tiller_por                   ! Tiller porosity
    INTEGER               :: iLayer                       ! layer counter

    tiller_area = 0._wp

    ! Root biomass distribution of all plants, modelled on Wania's 
    ! parameterisation for flood-tolerant graminoids
    ! ------------------------------------------------

    ! Wania uses dphen(d,pft) to have a phenology in the roots as
    ! they'll have to be built up in spring and plant transport
    ! will depend on existence of stems and leaves.
    ! We use ratio of LAI to LAI_max instead, which should give similar results

    ! Wania thinks that the limiting factor to plant transport are
    ! the number of tillers as they represent the
    ! cross-sectional area available to gas transport. She relates
    ! the biomass to tiller density based on Schimel 1995. The
    ! average biomass was 185g/m2 and the average number of
    ! tillers was 380/m2. This gives us 2.05 tillers per g
    ! biomass and 0.48 g per tiller, which corresponds to 0.22 g
    ! C per tiller

    ! TILLER SETUP for all gases
    ! --------------------------
    ! tiller # = phenology state x root mass / tiller_weight
    !
    ! We assume that 1/2 of Cpool_green is fine root mass
    ! we scale this with tiller_weight in mol C / m2
    n = (LAI / LAI_max) * (Cpool_green * 0.5_wp / tiller_weight)

    IF (is_woody) THEN    ! tree and shrub PFTs
      tiller_por = tiller_por_tree
    ELSE
      tiller_por = tiller_por_grass
    END IF
    IF (.NOT. is_wet) THEN
      tiller_por = tiller_por * dry_ratio
    END IF
    tiller_por = MAX(tiller_por,tiller_por_min)

    DO iLayer = 1, mxlay_above_rock
      tiller_frac = n * root_dist(iLayer)
      tiller_area(iLayer) = MAX(0.0001_wp, tiller_radius ** 2._wp * pi * tiller_frac * tiller_por)
    END DO

    IF (PRESENT(tillers)) tillers = n

  END SUBROUTINE plant_setup

  !!
  ! @brief     Setup gas distribution in soil
  !
  !!
  PURE SUBROUTINE init_gas_pools(is_peat, mxlay,                             &
                                 frac_inun_min, frac_inun,    frac_inun_pre, &
                                 CH4_dens_dry,  CO2_dens_dry, O2_dens_dry,   &
                                 CH4_dens_wet,  CO2_dens_wet, O2_dens_wet,   &
                                 CH4_dry,       CO2_dry,      O2_dry,        &
                                 CH4_wet,       CO2_wet,      O2_wet)
    LOGICAL,  INTENT(in)    :: is_peat            ! peatland flag
    INTEGER,  INTENT(in)    :: mxlay              ! layer number above bedrock
    REAL(wp), INTENT(in)    :: frac_inun_min      ! minimum inundated fraction of grid cell to regard for anoxic processes
    REAL(wp), INTENT(inout) :: frac_inun          ! inundated fraction of grid cell
    REAL(wp), INTENT(in)    :: frac_inun_pre      ! inundated fraction of grid cell in previous timestep
    ! These are the global values to be scaled on change of frac_inun - otherwise unchanged
    REAL(wp), INTENT(inout) :: CH4_dens_dry(:)    ! [mol m-2(tile)] CH4 density in layer, dry fraction
    REAL(wp), INTENT(inout) :: CO2_dens_dry(:)    ! [mol m-2(tile)] CO2 density in layer, dry fraction
    REAL(wp), INTENT(inout) :: O2_dens_dry(:)     ! [mol m-2(tile)]  O2 density in layer, dry fraction
    REAL(wp), INTENT(inout) :: CH4_dens_wet(:)    ! [mol m-2(tile)] CH4 density in layer, wet fraction
    REAL(wp), INTENT(inout) :: CO2_dens_wet(:)    ! [mol m-2(tile)] CO2 density in layer, wet fraction
    REAL(wp), INTENT(inout) :: O2_dens_wet(:)     ! [mol m-2(tile)]  O2 density in layer, wet fraction
    ! These are the local values scaled with the frac_inun
    REAL(wp), INTENT(out)   :: CH4_dry(:)         ! [mol m-2(dry frac)] CH4 density in layer, dry fraction
    REAL(wp), INTENT(out)   :: CO2_dry(:)         ! [mol m-2(dry frac)] CO2 density in layer, dry fraction
    REAL(wp), INTENT(out)   :: O2_dry (:)         ! [mol m-2(dry frac)]  O2 density in layer, dry fraction
    REAL(wp), INTENT(out)   :: CH4_wet(:)         ! [mol m-2(wet frac)] CH4 density in layer, wet fraction
    REAL(wp), INTENT(out)   :: CO2_wet(:)         ! [mol m-2(wet frac)] CO2 density in layer, wet fraction
    REAL(wp), INTENT(out)   :: O2_wet (:)         ! [mol m-2(wet frac)]  O2 density in layer, wet fraction

    CHARACTER(len=*), PARAMETER :: routine = modname//':init_gas_pools'

    ! locals
    REAL(wp) :: delta_inun                        ! change in inundated fraction
    REAL(wp) :: frac_transfer                     ! fraction to be transfered
    REAL(wp) :: tmp_GAS(mxlay)                    ! temporary GAS storage

    IF (.NOT. is_peat) THEN
      ! for mineral soils, gas densities in methane%GAS_WETNESS are relative to tile / gridcell 
      ! -> convert to densities relative to inundated / non-inundated area of grid cell
      ! changes in inundation require a gas transfer from wet to dry part or vice versa

      IF (frac_inun <  frac_inun_min) THEN   ! impose minimum inundation size
        frac_inun = 0._wp
      END IF
      frac_inun = MIN(frac_inun, 1._wp)   ! insure inundated fraction <= 1

      ! transfer amount and direction depends on change in inundation
      delta_inun = frac_inun - frac_inun_pre
      IF (delta_inun > EPSILON(1._wp)) THEN

        ! inundated area grows -> transfer from dry to wet part of tile
        frac_transfer = MIN(1._wp,delta_inun / (1._wp - frac_inun_pre))

        ! make sure we don't get incredible values from restart
        IF (frac_inun_pre < frac_inun_min) THEN
          CH4_dens_wet(1:mxlay) = 0._wp
          CO2_dens_wet(1:mxlay) = 0._wp
           O2_dens_wet(1:mxlay) = 0._wp
        END IF

        ! CH4
        tmp_GAS(1:mxlay) = CH4_dens_dry(1:mxlay) * frac_transfer
        CH4_dens_dry(1:mxlay) = CH4_dens_dry(1:mxlay) - tmp_GAS(1:mxlay)
        CH4_dens_wet(1:mxlay) = CH4_dens_wet(1:mxlay) + tmp_GAS(1:mxlay)

        ! CO2
        tmp_GAS(1:mxlay) = CO2_dens_dry(1:mxlay) * frac_transfer
        CO2_dens_dry(1:mxlay) = CO2_dens_dry(1:mxlay) - tmp_GAS(1:mxlay)
        CO2_dens_wet(1:mxlay) = CO2_dens_wet(1:mxlay) + tmp_GAS(1:mxlay)

        ! O2
        tmp_GAS(1:mxlay) = O2_dens_dry(1:mxlay) * frac_transfer
        O2_dens_dry (1:mxlay) = O2_dens_dry (1:mxlay) - tmp_GAS(1:mxlay)
        O2_dens_wet (1:mxlay) = O2_dens_wet (1:mxlay) + tmp_GAS(1:mxlay)

      ELSEIF (delta_inun < -EPSILON(1._wp)) THEN

        ! inundated area shrinks -> transfer from wet to dry part of tile
        frac_transfer = MIN(1._wp,-delta_inun / frac_inun_pre)

        ! CH4
        tmp_GAS(1:mxlay) = CH4_dens_wet(1:mxlay) * frac_transfer
        CH4_dens_dry(1:mxlay) = CH4_dens_dry(1:mxlay) + tmp_GAS(1:mxlay)
        CH4_dens_wet(1:mxlay) = CH4_dens_wet(1:mxlay) - tmp_GAS(1:mxlay)

        ! CO2
        tmp_GAS(1:mxlay) = CO2_dens_wet(1:mxlay) * frac_transfer
        CO2_dens_dry(1:mxlay) = CO2_dens_dry(1:mxlay) + tmp_GAS(1:mxlay)
        CO2_dens_wet(1:mxlay) = CO2_dens_wet(1:mxlay) - tmp_GAS(1:mxlay)

        ! O2
        tmp_GAS(1:mxlay) = O2_dens_wet(1:mxlay) * frac_transfer
        O2_dens_dry (1:mxlay) = O2_dens_dry (1:mxlay) + tmp_GAS(1:mxlay)
        O2_dens_wet (1:mxlay) = O2_dens_wet (1:mxlay) - tmp_GAS(1:mxlay)

        ! make sure we don't get incredible values from restart
        IF (frac_inun < frac_inun_min) THEN
          CH4_dens_wet(1:mxlay) = 0._wp
          CO2_dens_wet(1:mxlay) = 0._wp
          O2_dens_wet (1:mxlay) = 0._wp
        END IF

      ENDIF

    ELSE ! is_peat

      frac_inun = 1._wp

    ENDIF

    ! Define local arrays scaled to respective (non-)inundated parts of tile
    CH4_wet = CH4_dens_wet
    CO2_wet = CO2_dens_wet
     O2_wet =  O2_dens_wet
    CH4_dry = CH4_dens_dry
    CO2_dry = CO2_dens_dry
     O2_dry =  O2_dens_dry

    IF (frac_inun > 0._wp) THEN
      CH4_wet = CH4_dens_wet / frac_inun
      CO2_wet = CO2_dens_wet / frac_inun
       O2_wet =  O2_dens_wet / frac_inun
    ENDIF

    IF (frac_inun < 1._wp) THEN
      CH4_dry = CH4_dens_dry / (1._wp - frac_inun)
      CO2_dry = CO2_dens_dry / (1._wp - frac_inun)
       O2_dry =  O2_dens_dry / (1._wp - frac_inun)
    ENDIF

  END SUBROUTINE init_gas_pools

  !!
  ! @brief     Determine the CH4 / CO2 production and vertical distribution
  !
  !!
  PURE SUBROUTINE production(is_peatland, is_wet, mxlay_above_rock, mxlay_above_wt,        &
                        frac_CH4_anox_mineral, frac_CH4_anox_organic, meth_q10, root_dist, &
                        inun_frac, Tsoil, OxicDecompCarbon, AnoxDecompCarbon,              &
                        CH4_Prod, CO2_Prod, O2_Uptake, CH4_prod_ratio, resp_sl)
    LOGICAL,  INTENT(in)  :: is_peatland      ! peatland flag (.TRUE. if peatland)
    LOGICAL,  INTENT(in)  :: is_wet           ! wetness state, .FALSE. = dry, .TRUE. = wet
    INTEGER,  INTENT(in)  :: mxlay_above_rock ! layer number above bedrock
    INTEGER,  INTENT(in)  :: mxlay_above_wt   ! layer number above water table
    REAL(wp), INTENT(in)  :: frac_CH4_anox_mineral
    REAL(wp), INTENT(in)  :: frac_CH4_anox_organic
    REAL(wp), INTENT(in)  :: meth_q10
    REAL(wp), INTENT(in)  :: root_dist(:)     ! root distribution in each layer
    REAL(wp), INTENT(in)  :: inun_frac        ! grid cell inundated fraction []
    REAL(wp), INTENT(in)  :: Tsoil(:)         ! soil temperature [K]
    REAL(wp), INTENT(in)  :: OxicDecompCarbon ! Carbon production from oxic decomposition in tile [mol m-2 ts-1]
    REAL(wp), INTENT(in)  :: AnoxDecompCarbon ! Carbon production from anaerobic decomposition in tile [mol m-2 ts-1]
    REAL(wp), INTENT(out) :: CH4_Prod(:)      ! CH4 production per layer [mol m-2 ts-1]
    REAL(wp), INTENT(out) :: CO2_Prod(:)      ! CO2 production per layer [mol m-2 ts-1]
    REAL(wp), INTENT(out) :: O2_Uptake(:)     ! O2 uptake during oxic CO2 production per layer [mol m-2 ts-1]
    REAL(wp), INTENT(out) :: CH4_Prod_ratio   ! ratio of CH4 to total anoxic production []
    REAL(wp), OPTIONAL, INTENT(in) :: resp_sl(:) ! Distribution of respiration on individual soil layers if avail [mol m-2 ts-1]

    ! Parameters
    CHARACTER(len=*), PARAMETER :: routine = modname//':production'
    REAL(wp), PARAMETER :: Q10basetemp = 295._wp

    ! locals
    REAL(wp) :: carbon_dist   (mxlay_above_rock)     ! vertical distribution of carbon [mol m-2 ts-1]
    REAL(wp) :: CO2_Prod_oxic (mxlay_above_rock)     ! CO2 production per layer, oxic [mol m-2 ts-1]
    REAL(wp) :: CO2_Prod_anox (mxlay_above_rock)     ! CO2 production per layer, anoxic [mol m-2 ts-1]
    REAL(wp) :: root_dist_peat(mxlay_above_rock)     ! root distribution in peatland
    REAL(wp) :: f_ch4         (mxlay_above_rock)     ! fraction of C decomposed into CH4 [-]
    REAL(wp) :: carbon_oxic                          ! carbon decomposed oxically [-]
    REAL(wp) :: carbon_anox                          ! carbon decomposed anaerobically [-]
    REAL(wp) :: carbon_dist_sum_oxic                 ! sum of carbon in oxic layers [-]
    REAL(wp) :: carbon_dist_sum_anox                 ! sum of carbon in anoxic layer [-]
    REAL(wp) :: carbon_dist_sum_acro                 ! sum of carbon in anoxic acrotelm [-]
    REAL(wp) :: carbon_dist_sum_cato                 ! sum of carbon in catotelm [-]
    REAL(wp) :: carbon_dist_sum_oxic_inv             ! inverse of the sum of carbon in oxic layers [-]
    REAL(wp) :: carbon_dist_sum_anox_inv             ! inverse of the sum of carbon in anoxic layers [-]
    REAL(wp) :: carbon_dist_sum_acro_inv             ! inverse of the sum of carbon in anoxic acrotelm [-]
    REAL(wp) :: carbon_dist_sum_cato_inv             ! inverse of the sum of carbon in catotelm [-]
    REAL(wp) :: resp_sum

    CH4_Prod       = 0._wp
    CO2_Prod       = 0._wp
    CO2_Prod_oxic  = 0._wp
    CO2_Prod_anox  = 0._wp
    O2_Uptake      = 0._wp
    carbon_dist    = 0._wp
    root_dist_peat = 1._wp                ! we assume equally distributed C emissions per layer in peatland
    CH4_Prod_ratio = 0._wp

    ! production is different for peatlands and inundated areas:
    ! peatlands have vertically stacked oxic / anoxic layers, while
    ! inundated areas are completely anoxic in the inundated fraction
    ! and completely oxic in the non-inundated fraction
    IF (.NOT. is_peatland) THEN
      ! convert fluxes in non-peatlands to values per m2 in inundated / non-inundated fraction
      IF (inun_frac < 1._wp) THEN
        carbon_oxic = OxicDecompCarbon / (1._wp - inun_frac)
      ELSE
        carbon_oxic = 0._wp
      END IF

      IF (inun_frac > 0._wp) THEN
        carbon_anox = AnoxDecompCarbon / inun_frac
      ELSE
        carbon_anox = 0._wp
      END IF
    ELSE
      carbon_oxic = OxicDecompCarbon
      carbon_anox = AnoxDecompCarbon
    END IF ! .NOT. is_peatland

    ! In non-peatlands, inundated fraction is completely anoxic. This simplifies things substantially.
    ! we assume mineral soils in this case, using frac_CH4_anox_mineral
    IF (.NOT. is_peatland) THEN
      ! vertical distribution of anoxic carbon - entire soil column is either saturated or not
      IF (.NOT. is_wet) THEN

        IF (PRESENT(resp_sl)) THEN
          resp_sum = SUM(resp_sl(1:mxlay_above_rock))
          IF (resp_sum /= 0._wp) THEN ! If not, carbon_oxic should also be 0
            ! Distribute according to respiration
            carbon_dist(1:mxlay_above_rock) = carbon_oxic * resp_sl(1:mxlay_above_rock) / resp_sum
          ENDIF
        ELSE
          ! carbon is input according to root distribution
          carbon_dist(1:mxlay_above_rock) = root_dist(1:mxlay_above_rock) * carbon_oxic
        ENDIF

        ! dry part of grid cell is completely oxic, becoming CO2 
        CO2_Prod_oxic(1:mxlay_above_rock) = carbon_dist(1:mxlay_above_rock)
      ELSE
        f_ch4(:) = MIN(0.5_wp,frac_CH4_anox_mineral * meth_q10**((Tsoil(1:mxlay_above_rock)-Q10basetemp)/10._wp))
        IF (PRESENT(resp_sl)) THEN
          resp_sum = SUM(resp_sl(1:mxlay_above_rock))
          IF (resp_sum /= 0._wp) THEN ! If not, carbon_oxic should also be 0
            ! Distribute according to respiration
            carbon_dist(1:mxlay_above_rock) = carbon_oxic * resp_sl(1:mxlay_above_rock) / resp_sum
          ENDIF
        ELSE
          ! carbon is input according to root distribution
          carbon_dist(1:mxlay_above_rock) = root_dist(1:mxlay_above_rock) * carbon_anox
        ENDIF

        ! C decomposed under anaerobic conditions becomes both CH4 and CO2, with
        ! fractions determined by frac_CH4_anox_mineral
        ! determine CO2 densities and concentrations produced anaerobically
        CO2_Prod_anox(1:mxlay_above_rock) = (1._wp - f_ch4(1:mxlay_above_rock)) * carbon_dist(1:mxlay_above_rock)

        ! determine CH4 densities and concentrations produced
        CH4_Prod(1:mxlay_above_rock) = f_ch4(1:mxlay_above_rock) * carbon_dist(1:mxlay_above_rock)
      END IF
    ELSE   ! is_peatland
      f_ch4(:) = MIN(0.5_wp,frac_CH4_anox_organic * meth_q10**((Tsoil(1:mxlay_above_rock)-Q10basetemp)/10._wp))

      ! we assume organic soils in this case, using frac_CH4_anox_organic

      ! sum of carbon in oxic, anoxic acrotelm and catotelm
      ! NEED TO MODIFY THIS ONCE WE HAVE CATOTELM and also for the availability of layered respiration
      carbon_dist_sum_oxic = SUM(root_dist_peat(1:mxlay_above_wt))
      carbon_dist_sum_anox = SUM(root_dist_peat((mxlay_above_wt+1):mxlay_above_rock))

      ! inverse of the carbon distribution sum
      IF (carbon_dist_sum_oxic > 0._wp) THEN
        carbon_dist_sum_oxic_inv = 1._wp / carbon_dist_sum_oxic
      ELSE
        carbon_dist_sum_oxic_inv = 0._wp
      END IF
      IF (carbon_dist_sum_anox > 0._wp) THEN
        carbon_dist_sum_anox_inv = 1._wp / carbon_dist_sum_anox
      ELSE
        carbon_dist_sum_anox_inv = 0._wp
      END IF

      ! calculate the concentration of CO2 out of the
      ! density of decomposed carbon, assuming, that all
      ! carbon above the water table and half of the
      ! carbon below the water table becomes CO2

      ! determine CO2 densities produced above the water table
      ! consider very high water table as well to conserve carbon
      IF ((mxlay_above_wt <= 1._wp) .OR. (carbon_dist_sum_oxic_inv < 1.e-12_wp)) THEN
        CO2_Prod_oxic(1) = carbon_oxic
      ELSE
        CO2_Prod_oxic(1:mxlay_above_wt) = &
                root_dist_peat(1:mxlay_above_wt) * carbon_dist_sum_oxic_inv * carbon_oxic
      END IF

      ! vertical distribution of anoxic carbon in saturated part
      ! consider very low water table as well to conserve carbon
      IF ((mxlay_above_wt < mxlay_above_rock) .AND. (carbon_dist_sum_anox_inv > 1.e-12_wp)) THEN
        carbon_dist((mxlay_above_wt+1):mxlay_above_rock) = &
                root_dist_peat((mxlay_above_wt+1):mxlay_above_rock) * carbon_dist_sum_anox_inv * carbon_anox
      ELSE
        carbon_dist(mxlay_above_rock) = carbon_anox
      END IF

      ! determine CO2 densities produced below the water table
      CO2_Prod_anox(1:mxlay_above_rock) = (1._wp - f_ch4(1:mxlay_above_rock)) * carbon_dist(1:mxlay_above_rock)
      ! determine CH4 densities produced
      CH4_Prod(1:mxlay_above_rock) = f_ch4(1:mxlay_above_rock) * carbon_dist(1:mxlay_above_rock)
    END IF ! NOT is_peatland

    CO2_Prod (1:mxlay_above_rock) =  CO2_Prod_oxic(1:mxlay_above_rock) + CO2_Prod_anox(1:mxlay_above_rock)
    O2_Uptake(1:mxlay_above_rock) = -CO2_Prod(1:mxlay_above_rock) * 2._wp   ! 2 mol of O2 for 1 mol of CO2
    ! StW: The factor of 2 is a special case, assuming that the carbon is in the form of CH4 (other species may
    !      also fulfill this criteria), but for general hydrocarbons it should be x + y/4, where x is the 
    !      #C-atoms in the compound and y is the #H-atoms in the compound. E.g. For propane the factor is 5/3.
    !      Further modifications apply if the compound also contains oxygen.

    IF (SUM(CH4_Prod(1:mxlay_above_rock)+CO2_Prod_anox(1:mxlay_above_rock)) > 0._wp) &
      CH4_Prod_ratio = SUM(CH4_prod(1:mxlay_above_rock))/SUM(CH4_Prod(1:mxlay_above_rock)+CO2_Prod_anox(1:mxlay_above_rock))

  END SUBROUTINE production

  !!
  ! @brief     oxidation following Michaelis-Menten kinetics
  !
  !!
  PURE SUBROUTINE oxidation_mm(is_wet, dt, oxid_scale_dry, oxid_scale_wet, oxid_q10, Tsoil, rel_soil_moist, &
                                         dens2conc, conc2dens, CH4_dens, CO2_dens, O2_dens, O2_Uptake_Prod, &
                                         CH4_Oxid, CO2_Oxid, O2_Oxid, O2_Other)

    USE mo_jsb_physical_constants, ONLY: tmelt

    LOGICAL,  INTENT(in)  :: is_wet            ! wetness state: .FALSE. = dry, .TRUE. = wet
    REAL(wp), INTENT(in)  :: dt                ! time step [s]
    REAL(wp), INTENT(in)  :: oxid_scale_dry    ! 
    REAL(wp), INTENT(in)  :: oxid_scale_wet    !
    REAL(wp), INTENT(in)  :: oxid_q10
    REAL(wp), INTENT(in)  :: Tsoil             ! Soil temperature [K]
    REAL(wp), INTENT(in)  :: rel_soil_moist    ! relative soil moisture
    REAL(wp), INTENT(in)  :: dens2conc         ! conversion factor from density to concentration []
    REAL(wp), INTENT(in)  :: conc2dens         ! conversion factor from concentration to density []
    REAL(wp), INTENT(in)  :: CH4_dens          ! CH4 density in layer [mol m-2]
    REAL(wp), INTENT(in)  :: CO2_dens          ! CO2 density in layer [mol m-2]
    REAL(wp), INTENT(in)  :: O2_dens           ! O2  density in layer [mol m-2]
    REAL(wp), INTENT(in)  :: O2_Uptake_Prod    ! O2  uptake during CO2 production [mol m-2 ts-1]
    REAL(wp), INTENT(out) :: CH4_Oxid          ! CH4 density oxidised in layer [mol m-2 ts-1]
    REAL(wp), INTENT(out) :: CO2_Oxid          ! CO2 density oxidised in layer [mol m-2 ts-1]
    REAL(wp), INTENT(out) :: O2_Oxid           ! O2  density oxidised in layer [mol m-2 ts-1]
    REAL(wp), INTENT(out) :: O2_Other          ! O2  density utilised by other processes [mol m-2 ts-1]

    ! Parameters Riley
    CHARACTER(len=*), PARAMETER :: routine = modname//':oxidation_mm'
    REAL(wp), PARAMETER :: v_Max_wet = 1.25e-5_wp           ! max. reaction velocity (based on W&H, 2000) [mol m-3 s-1]
    REAL(wp), PARAMETER :: K_M_CH4_wet = 5.e-3_wp           ! Michaelis-Menten constant [mol m-3]
                                                            !   CH4: concentration where v = 0.5 * v_Max (W&H, 2000)
    REAL(wp), PARAMETER :: v_Max_dry = v_Max_wet/10._wp     ! max. reaction velocity (reduced for upland) [mol m-3 s-1]
    REAL(wp), PARAMETER :: K_M_CH4_dry = K_M_CH4_wet/10._wp ! Michaelis-Menten constant [mol m-3]
                                                            !   CH4: concentration where v = 0.5 * v_Max (reduced for upland)
    REAL(wp), PARAMETER :: K_M_O2  = 2.e-2_wp               ! Michaelis-Menten constant [mol m-3]
                                                            !   O2: concentration where v = 0.5 * v_Max
    REAL(wp), PARAMETER :: ten = 10._wp                     ! divisor for Q10 exponent []
    REAL(wp), PARAMETER :: stoich = 2._wp                   ! stoichiometric CH4 oxidation: 2mol O2 for 1mol CH4 []
    REAL(wp), PARAMETER :: frac_O2_Ox  = 0.4_wp             ! fraction of O2 needed for other oxidation processes []
    REAL(wp), PARAMETER :: frac_O2_Ph  = 0.1_wp             ! fraction of O2 needed for other physical processes []
    REAL(wp), PARAMETER :: frac_O2_Other = frac_O2_Ox + frac_O2_Ph ! fraction of O2 needed for other processes []

    ! local variable
    INTEGER :: iLayer
    ! REAL(wp) :: O2_avail
    REAL(wp) :: MMkin
    REAL(wp) :: moist_lim                                   ! moisture limitation
    REAL(wp) :: v_Max_eff
    REAL(wp) :: K_M_CH4_eff

    REAL(wp) :: CH4_conc                                    ! CH4 concentration [mol m-3]
    REAL(wp) :: CH4_conc_oxid                               ! CH4 concentration oxidized [mol m-3]
    REAL(wp) :: O2_conc                                     ! O2 concentration [mol m-3]
    REAL(wp) :: O2_dens_loc  

    CH4_conc_oxid = 0._wp

    IF (.NOT. is_wet) THEN
      v_Max_eff = v_Max_dry*oxid_scale_dry
      K_M_CH4_eff = K_M_CH4_dry
    ELSE
      v_Max_eff = v_Max_wet*oxid_scale_wet
      K_M_CH4_eff = K_M_CH4_wet
    END IF

    ! other O2 requirements: uptake during CO2 production, physical processes, other oxidation processes
    ! maximum requirement is half availability
    O2_Other = -MIN(-O2_Uptake_Prod + frac_O2_Ph * O2_dens, frac_O2_Other * O2_dens)
    O2_dens_loc = O2_dens + O2_Other

    CH4_conc = CH4_dens    * dens2conc
    O2_conc  = O2_dens_loc * dens2conc

    moist_lim = 1.0_wp
    IF (CH4_conc > 0._wp) THEN
      ! O2 availability - StW: Not used?!
      ! O2_avail = O2_conc / stoich
      ! Michaelis-Menten kinetics
      MMkin = v_Max_eff * &
              (CH4_conc / (K_M_CH4_eff + CH4_conc)) * &
              (O2_conc  / (K_M_O2      + O2_conc )) * &
              oxid_q10**((Tsoil - (tmelt + ten)) / ten) * &
              moist_lim * dt
      CH4_conc_oxid = -MMkin
    END IF

    CH4_Oxid = CH4_conc_oxid * conc2dens
    CO2_Oxid = -CH4_Oxid
    O2_Oxid  = CH4_Oxid * stoich

  END SUBROUTINE oxidation_mm

  !!
  ! @brief      plant transport for GAS
  !             determines GAS transported through plants
  !
  !!
  PURE SUBROUTINE plant_transport(GAS, layer_thick, layer_depth, mxlay_above_rock, mxlay_above_wt, &
                                  dt, Tsoil, relative_porevol, dens2conc, conc2dens, tiller_area,  &
                                  D_GAS, GAS_conc_atm, GAS_dens, GAS_Plant, GAS_excess)

    USE mo_methane_constants, ONLY: t_specie_const, minimum_porevol
    USE mo_methane_util,      ONLY: determine_Henry_cc

    TYPE (t_specie_const), INTENT(in)  :: GAS                 ! Constants for specific gas
    REAL(wp),              INTENT(in)  :: layer_thick(:)      ! soil layer thickness [m]
    REAL(wp),              INTENT(in)  :: layer_depth(:)      ! soil layer depth [m]
    INTEGER,               INTENT(in)  :: mxlay_above_rock    ! layer number above bedrock
    INTEGER,               INTENT(in)  :: mxlay_above_wt      ! layer number above water table
    REAL(wp),              INTENT(in)  :: dt                  ! time step [s]
    REAL(wp),              INTENT(in)  :: Tsoil(:)            ! Soil temperature [K]
    REAL(wp),              INTENT(in)  :: relative_porevol(:) ! relative pore volume in layer []
    REAL(wp),              INTENT(in)  :: dens2conc(:)        ! conversion factor from density to concentration [m-1]
    REAL(wp),              INTENT(in)  :: conc2dens(:)        ! conversion factor from concentration to density [m]
    REAL(wp),              INTENT(in)  :: tiller_area(:)      ! area of tillers corresponding to the root fraction in each layer
    REAL(wp),              INTENT(in)  :: D_GAS               ! Diffusivity of GAS (currently = SURFACE AIR diffusivity) [m2 s-1]
    REAL(wp),              INTENT(in)  :: GAS_conc_atm        ! GAS concentration in atmosphere [mol m-3]
    REAL(wp),              INTENT(in)  :: GAS_dens(:)         ! GAS density in ice free pore space of layer [mol m-2]
    REAL(wp),              INTENT(out) :: GAS_Plant(:)        ! GAS density transported by plants per layer [mol m-2 ts-1]
    REAL(wp), OPTIONAL,    INTENT(out) :: GAS_excess(:)       ! max(0,GAS_Plant_uncorr+GAS_dens) [mol m-2 ts-1]

    ! Parameters
    CHARACTER(len=*), PARAMETER :: routine = modname//':plant_transport'
    REAL(wp),         PARAMETER :: root_length_ratio = 3.0_wp  ! [m m-1] ratio of root length to depth

    ! locals
    REAL(wp) :: GAS_conc              (mxlay_above_rock) ! GAS concentration in ice free pore space of layer [mol m-3]
    REAL(wp) :: GAS_conc_gaseous      (mxlay_above_rock) ! GAS concentration in gaseous form [mol m-3]
    REAL(wp) :: kHCC                  (mxlay_above_rock) ! Henry factor relating gaseous to aqueous concentration []
    REAL(wp) :: effective_layer_volume(mxlay_above_rock) ! effective volume of layer [m]
    REAL(wp) :: J                                        ! GAS flux in layer [mol m-2 s-1]
    INTEGER  :: is                                       ! layer counter

    GAS_conc   = 0._wp
    GAS_Plant  = 0._wp
    IF (PRESENT(GAS_excess)) GAS_excess = 0._wp

    DO is = 1, mxlay_above_rock
      effective_layer_volume(is) = conc2dens(is)

      kHCC(is) = determine_Henry_cc(GAS, Tsoil(is))

      ! Convert GAS densities [mol layer-1] to concentrations [mol m-3]
      GAS_conc(is) = GAS_dens(is) * dens2conc(is)
    ENDDO

    ! determine concentration of gas in gaseous form: Below water table this is GAS_conc/kHCC
    DO is = 1, mxlay_above_wt
      GAS_conc_gaseous(is) = GAS_conc(is)
    ENDDO
    DO is = mxlay_above_wt+1, mxlay_above_rock
      GAS_conc_gaseous(is) = GAS_conc(is) / kHCC(is)
    ENDDO

    DO is = 1, mxlay_above_rock
      IF (relative_porevol(is) > minimum_porevol) THEN
        J = -D_GAS * tiller_area(is) * (GAS_conc_gaseous(is) - GAS_conc_atm) / (layer_depth(is) * root_length_ratio)
        GAS_Plant(is) = J * dt * (effective_layer_volume(is) / layer_thick(is))
        IF (PRESENT(GAS_excess)) THEN
          IF (GAS_Plant(is) + GAS_dens(is) < 0._wp) THEN
            GAS_excess(is) =  GAS_dens(is) + GAS_Plant(is)
            GAS_Plant (is) = -GAS_dens(is)
          ENDIF
        ENDIF
      END IF
    END DO

  END SUBROUTINE plant_transport

  !!
  ! @brief      calculate ebulition transport of CH4
  !
  !!
  !  PURE 
  SUBROUTINE ebullition(is_peatland, mxlay_above_rock, mxlay_above_wt, Tsoil, surface_air_pressure, &
                        Dz, relative_porevol, water_table, &
                        dens2conc, conc2dens, CH4_dens, CH4_ebul)

    USE mo_jsb_physical_constants, ONLY: UniversalGasConstant => argas, tmelt, rhoh2o, grav
    USE mo_methane_util,           ONLY: dissolved_gas
    USE mo_methane_constants,      ONLY: CH4_consts

    LOGICAL,  INTENT(in)  :: is_peatland          ! peatland flag
    INTEGER,  INTENT(in)  :: mxlay_above_rock     ! layer number above bedrock
    INTEGER,  INTENT(in)  :: mxlay_above_wt       ! layer number above water table
    REAL(wp), INTENT(in)  :: Tsoil(:)             ! Soil temperature [K]
    REAL(wp), INTENT(in)  :: surface_air_pressure ! surface air pressure [Pa]
    REAL(wp), INTENT(in)  :: Dz(:)                ! distance from soil layer above [m]
    REAL(wp), INTENT(in)  :: relative_porevol(:)  ! relative pore volume in layer []
    REAL(wp), INTENT(in)  :: water_table          ! water table position [m]
    REAL(wp), INTENT(in)  :: dens2conc(:)         ! conversion factor from density to concentration []
    REAL(wp), INTENT(in)  :: conc2dens(:)         ! conversion factor from concentration to density []
    REAL(wp), INTENT(in)  :: CH4_dens(:)          ! CH4 density in ice free pore space of layer [mol m-2]
    REAL(wp), INTENT(out) :: CH4_ebul(:)          ! CH4 density per layer transported by ebullition [mol m-2 ts-1]

    ! Parameters
    CHARACTER(len=*), PARAMETER :: routine = modname//':ebullition'
    REAL(wp), PARAMETER   :: vgc_low         = 0.145_wp
    REAL(wp), PARAMETER   :: vgc_high        = 0.15_wp
    REAL(wp), PARAMETER   :: bubble_CH4_frac = 0.57_wp

    ! locals
    REAL(wp)              :: CH4_diss_max               ! maximum concentration of dissolved CH4 [ml CH4 mlH2O]
    REAL(wp)              :: CH4_diss_max_m3            ! [mol m-3] maximum dissolved CH4 concentration
    REAL(wp)              :: CH4_diss_max_mol           ! [mol layer-1] maximum dissolved CH4 density
    REAL(wp)              :: CH4_gas (mxlay_above_rock) ! [mol] CH4 density as gas in layer
    REAL(wp)              :: CH4_diss(mxlay_above_rock) ! [mol] CH4 density dissolved in layer
    REAL(wp)              :: effective_layer_volume(mxlay_above_rock) ! effective pore volume of layer [m]
    REAL(wp)              :: water_height               ! height of water at current layer
    REAL(wp)              :: hydro_press                ! [Pa] hydrostatic pressure caused by water above layer
    INTEGER               :: is                         ! layer counter
    LOGICAL               :: newebul_flag               ! swich to Wania's new ebullotion scheme
    REAL(wp)              :: CH4_gas_vol
    REAL(wp)              :: CH4_ebul_vol
    REAL(wp)              :: CH4_vgc

    newebul_flag = .FALSE. ! StW: TODO: should this be a configurable parameter?

    DO is = 1, mxlay_above_rock
      effective_layer_volume(is) = conc2dens(is)
      CH4_Ebul(is) = 0._wp
      CH4_diss(is) = 0._wp
      CH4_gas (is) = 0._wp
    ENDDO
  
    ! water height is >= 0
    ! water table should be > 0 in inundated areas, <= 0 in peatlands
    ! water height is height of wt above surface
    water_height = MAX(water_table, 0._wp)

    DO is = 1, mxlay_above_rock
      CH4_diss(is) = dissolved_gas(CH4_consts, Tsoil(is), CH4_dens(is), dens2conc(is))
      CH4_gas (is) = CH4_dens(is) - CH4_diss(is)
    ENDDO

    ! SOLUBILITY OF CH4
    ! -----------------
    ! Estimate how much methane can be dissolved at current T
    ! from Yamamoto et al. 1976. Yamamoto's relationship
    ! holds only up to about 37 degrees, after which the
    ! dissolubility increases again when using the eqn below.
    ! UNITS: [ml CH4 per ml H20]

    DO is = mxlay_above_wt+1, mxlay_above_rock

      CH4_diss_max = 0.05708_wp - 0.001545_wp * (MAX(tmelt, Tsoil(is)) - tmelt) + &
                     0.00002069_wp * ((MAX(tmelt, Tsoil(is)) - tmelt) ** 2._wp) ![ml CH4 ml-1 H2O]

      ! below water table: increase water_height by layer thickness
      IF (is > mxlay_above_wt) THEN ! TODO: Isn't this by definition fullfilled?
        water_height = water_height + Dz(is)  ![m]
      END IF

      ! Calculate hydrostatic pressure
      hydro_press = rhoh2o * grav * water_height ![kg m s-2 or N]

      ! [mL CH4 mL-1 H2O] is equivalent to [m3 CH4 m-3 H2O in soil] 
      ! To calculate the volume of CH4 per layer, multiply by the
      ! the layer's water content V_water.
      ! We use effective_layer_volume
      ! -----------------------------------------------------------
      CH4_diss_max_m3 = CH4_diss_max * effective_layer_volume(is) ![m3 CH4 layer-1]

      ! Use the Ideal Gas Law to convert [m3 CH4 layer-1] to [mol CH4 layer-1]
      ! Add hydrostatic pressure to atmospheric pressure.
      ! pressure is in [Pa], which is [N m-2], so hydro_press
      ! has to be divided by 1m2 to have the units [N m-2] or [Pa]
      ! ----------------------------------------------------------
      CH4_diss_max_mol = CH4_diss_max_m3 * (surface_air_pressure + hydro_press) / (UniversalGasConstant * Tsoil(is))

      ! If the total density of CH4 is higher than what can be
      ! dissolved in the given water volume times the
      ! threshold for bubble initiation, then increase
      ! the gaseous CH4 pool.
      ! Units for concentrations are in [g CH4-C layer-1]
      ! Units for volume is in [m3]
      ! ------------------------------------------------------
      ! determine dissolved CH4
      IF (newebul_flag) THEN
        CH4_diss(is) = MIN(CH4_dens(is), CH4_diss_max_mol) ![mol CH4 layer-1]

        ! determine CH4 as gas in layer
        CH4_gas(is) = CH4_dens(is) - CH4_diss(is) ![mol CH4 m-2]

        ! ebullition flux is gas density in layer (but negative as flux removes from soil)
        CH4_Ebul(is) = -CH4_gas(is) ! [mol CH4 layer-1 timestep-1] 
        CH4_gas (is) = 0.0_wp
      ELSE

        CH4_gas_vol = CH4_gas(is) * ((UniversalGasConstant * Tsoil(is)) / (surface_air_pressure + hydro_press))
        CH4_vgc     = CH4_gas_vol / effective_layer_volume(is)

        IF (CH4_vgc > vgc_high*bubble_CH4_frac) THEN
          CH4_ebul_vol = (CH4_vgc - (vgc_low*bubble_CH4_frac)) * effective_layer_volume(is)
          CH4_vgc      = vgc_low*bubble_CH4_frac
          CH4_gas_vol  = CH4_vgc * effective_layer_volume(is)
          CH4_gas (is) =  CH4_gas_vol  / ((UniversalGasConstant * Tsoil(is)) / (surface_air_pressure + hydro_press))
          CH4_Ebul(is) = -CH4_ebul_vol / ((UniversalGasConstant * Tsoil(is)) / (surface_air_pressure + hydro_press))
        END IF

      END IF

    END DO
    ! ebullition flux SHOULD be added to soil just above water table, if wt is below surface
    ! done in reaction-diffusion

  END SUBROUTINE ebullition

  !!
  ! @brief      reaction-diffusion for all soil layers and all gases
  !
  !!
  SUBROUTINE reaction_diffusion(nstep, implicit_contribution,                  &
                   dt, mxlay_above_rock, mxlay_above_wt,                       &
                   layer_thick, Tsoil, relative_porevol, dens2conc, conc2dens, &
                   CH4_Prod, CO2_Prod, CH4_Plant, CO2_Plant, O2_Plant,         &
                   D_CH4, D_CO2, D_O2, D_CH4_surf, D_CO2_surf, D_O2_surf,      &
                   surf_CH4_conc, surf_CO2_conc, surf_O2_conc,                 &
                   CH4_dens, CO2_dens, O2_dens,                                &
                   CH4_Oxid, CO2_Oxid, O2_Oxid, CH4_Ebul, O2_Other,            &
                   surf_CH4_flux, surf_CO2_flux, surf_O2_flux, D_used, lissue, nissue)

    USE mo_methane_util,      ONLY: determine_Henry_cc, diffusion_cn
    USE mo_methane_constants, ONLY: CH4_consts, CO2_consts, O2_consts

    INTEGER,  INTENT(in)    :: nstep                  ! #diffusion time steps per model time step
    REAL(wp), INTENT(in)    :: implicit_contribution  ! Implicit contribution for Crank-Nicolson
    REAL(wp), INTENT(in)    :: dt                     ! Time step [s]
    INTEGER,  INTENT(in)    :: mxlay_above_rock       ! layer number above bedrock
    INTEGER,  INTENT(in)    :: mxlay_above_wt         ! layer number above water table
    REAL(wp), INTENT(in)    :: Tsoil(:)               ! [K] Soil temperature [K]
    REAL(wp), INTENT(in)    :: layer_thick(:)         ! soil layer thickness [m]
    REAL(wp), INTENT(in)    :: relative_porevol(:)    ! relative pore volume in layer []
    REAL(wp), INTENT(in)    :: dens2conc(:)           ! conversion factor from density to concentration []
    REAL(wp), INTENT(in)    :: conc2dens(:)           ! conversion factor from concentration to density []
    REAL(wp), INTENT(in)    :: CH4_Prod(:)            ! CH4 production in layer [mol m-2 ts-1]
    REAL(wp), INTENT(in)    :: CO2_Prod(:)            ! CO2 production in layer [mol m-2 ts-1]
    REAL(wp), INTENT(in)    :: CO2_Plant(:)           ! CO2 plant flux on layer [mol m-2 ts-1]
    REAL(wp), INTENT(in)    :: O2_Plant(:)            ! O2 plant  flux on layer [mol m-2 ts-1]
    REAL(wp), INTENT(in)    :: D_CH4(:)               ! diffusivity of CH4 on layers  [m2 s-1]
    REAL(wp), INTENT(in)    :: D_CO2(:)               ! diffusivity of CO2 on layers  [m2 s-1]
    REAL(wp), INTENT(in)    :: D_O2(:)                ! diffusivity of O2  on layers  [m2 s-1]
    REAL(wp), INTENT(in)    :: D_CH4_surf             ! diffusivity of CH4 on surface [m2 s-1]
    REAL(wp), INTENT(in)    :: D_CO2_surf             ! diffusivity of CO2 on surface [m2 s-1]
    REAL(wp), INTENT(in)    :: D_O2_surf              ! diffusivity of O2  on surface [m2 s-1]
    REAL(wp), INTENT(in)    :: surf_CH4_conc          ! concentration of CH4 on surface [mol m-3]
    REAL(wp), INTENT(in)    :: surf_CO2_conc          ! concentration of CO2 on surface [mol m-3]
    REAL(wp), INTENT(in)    :: surf_O2_conc           ! concentration of O2 on surface  [mol m-3]
    REAL(wp), INTENT(inout) :: CH4_dens(:)            ! CH4 density in ice free pore space of layer [mol m-2]
    REAL(wp), INTENT(inout) :: CO2_dens(:)            ! CO2 density in ice free pore space of layer [mol m-2]
    REAL(wp), INTENT(inout) :: O2_dens(:)             ! O2  density in ice free pore space of layer [mol m-2]
    REAL(wp), INTENT(inout) :: CH4_Oxid(:)            ! CH4 removal through oxidation in layer  [mol m-2 ts-1]
    REAL(wp), INTENT(inout) :: CO2_Oxid(:)            ! CO2 production from oxidation in layer  [mol m-2 ts-1]
    REAL(wp), INTENT(inout) :: O2_Other(:)            ! O2 consumption through other processes  [mol m-2 ts-1]
    REAL(wp), INTENT(inout) :: O2_Oxid(:)             ! O2 consumption through oxidation of CH4 [mol m-2 ts-1]
    REAL(wp), INTENT(inout) :: CH4_Ebul(:)            ! CH4 flux by ebullition on layert [mol m-2 ts-1]
    REAL(wp), INTENT(inout) :: CH4_Plant(:)           ! CH4 plant flux on layer [mol m-2 ts-1]
    REAL(wp), INTENT(out)   :: surf_CH4_flux          ! CH4 flux on surface [mol m-2 ts-1]
    REAL(wp), INTENT(out)   :: surf_CO2_flux          ! CO2 flux on surface [mol m-2 ts-1]
    REAL(wp), INTENT(out)   :: surf_O2_flux           ! O2  flux on surface [mol m-2 ts-1]
    REAL(wp), INTENT(out)   :: D_used(:)
    LOGICAL,  INTENT(out), OPTIONAL :: lissue
    INTEGER,  INTENT(out), OPTIONAL :: nissue

    ! locals
    CHARACTER(len=*), PARAMETER :: routine = modname//':reaction_diffusion'
    REAL(wp) :: O2_demand                             ! O2  demand in layer [mol m-2 ts-1]
    REAL(wp) :: CH4_demand                            ! CH4 demand in layer [mol m-2 ts-1]
    REAL(wp) :: O2_stress                             ! O2 stress factor []
    REAL(wp) :: CH4_stress                            ! CH4 stress factor []
    REAL(wp) :: CH4_Ebul_sum                          ! total Ch4 transported by ebullition [mol m-2 ts-1]
    REAL(wp) :: Di         (mxlay_above_rock)         ! diffusivity of GAS on layers used [m2 s-1]
    REAL(wp) :: CH4_source (mxlay_above_rock)         ! CH4 source / sink on layer [mol m-2 ts-1]
    REAL(wp) :: CO2_source (mxlay_above_rock)         ! CO2 source / sink on layer [mol m-2 ts-1]
    REAL(wp) :: O2_source  (mxlay_above_rock)         ! O2  source / sink on layer [mol m-2 ts-1]
    REAL(wp) :: kHCC       (mxlay_above_rock)         ! Henry coefficient relating aqueous to gaseous concentration []
    INTEGER :: is, nissue_loc
    LOGICAL :: lissue_loc, lissue_tst

    DO is = 1, mxlay_above_rock
      CH4_source(is) = 0._wp
      CO2_source(is) = 0._wp
      O2_source (is) = 0._wp
    ENDDO
    nissue_loc = 0

    DO is = 1,mxlay_above_rock
      O2_demand = -(O2_Oxid(is) + O2_Other(is))   ! O2_Oxid and O2_Other are defined negative
      O2_stress = 1._wp
      IF (O2_demand > 0._wp) THEN
        O2_stress = MIN((O2_dens(is) + O2_Plant(is)) / O2_demand, 1._wp)
        IF (O2_stress < 0._wp) THEN
          nissue_loc = nissue_loc + 1
          O2_stress = 0._wp
        ENDIF
      END IF

      CH4_demand = -(CH4_Oxid(is) + CH4_Plant(is) + CH4_Ebul(is))  ! demand variables are defined negative
      CH4_stress = 1._wp
      IF (CH4_demand > 0._wp) THEN
        CH4_stress = MIN((CH4_dens(is) + CH4_Prod(is)) / CH4_demand, 1._wp)
      END IF

      ! resolve oxidation stress
      IF (CH4_stress <= O2_stress) THEN   ! CH4 limited
        IF (O2_stress < 1) THEN
          O2_demand = -O2_Other(is)   ! O2_Oxid and O2_Other are defined negative
          O2_stress = 1._wp
          IF (O2_demand > 0._wp) THEN
            O2_stress = MIN((O2_dens(is) + O2_Plant(is) + &   ! should be -, but O2_oxid defined negative
                        (CH4_stress * O2_Oxid(is))) / O2_demand, 1._wp)
          END IF
        END IF
        CH4_Oxid(is) = CH4_Oxid(is) * CH4_stress
         O2_Oxid(is) =  O2_Oxid(is) * CH4_stress
      ELSE    ! O2 limited
        IF (CH4_stress < 1._wp) THEN
          CH4_demand = -(CH4_Ebul(is) + CH4_Plant(is))
          CH4_stress = 1._wp
          IF (CH4_demand > 0._wp) THEN
            CH4_stress = MIN((CH4_dens(is) + CH4_Prod(is) + &   ! should be -, but CH4_oxid defined negative
                          (O2_stress * CH4_Oxid(is))) / CH4_demand, 1._wp)
          END IF
        END IF
        CH4_Oxid(is) = CH4_Oxid(is) * O2_stress
         O2_Oxid(is) =  O2_Oxid(is) * O2_stress
      END IF
      CO2_Oxid (is) = -CH4_Oxid (is)
      CH4_Ebul (is) =  CH4_Ebul (is)  * CH4_stress
      CH4_Plant(is) =  CH4_Plant(is)  * CH4_stress
       O2_Other(is) =   O2_Other(is)  * O2_stress

    END DO   ! iLayer

    ! add density from CH4 ebullition to layer above water table, if water table is below surface
    CH4_Ebul_sum = SUM(CH4_Ebul(1:mxlay_above_rock))   ! this is NEGATIVE, since CH4_Ebul is defined negative
    IF (mxlay_above_wt > 0) THEN
      CH4_Ebul(mxlay_above_wt) = CH4_Ebul(mxlay_above_wt)-CH4_Ebul_sum
      CH4_Ebul_sum = 0._wp
    END IF

    DO is = 1, mxlay_above_rock
      CH4_source(is) = (CH4_Prod (is) + CH4_Oxid(is) + CH4_Plant(is) + CH4_Ebul(is))
      CO2_source(is) = (CO2_Prod (is) + CO2_Oxid(is) + CO2_Plant(is))
      O2_source (is) = (O2_Other (is) +  O2_Oxid(is) +  O2_Plant(is))
    ENDDO

    lissue_loc = .FALSE.
    ! diffusion calculation: CH4

    ! set up diffusion coefficients
    DO is = 1, mxlay_above_rock
      Di    (is) = D_CH4(is)
      kHCC  (is) = determine_Henry_cc(CH4_consts, Tsoil(is))
      D_used(is) = surf_CH4_flux
    ENDDO

    CALL diffusion_cn(nstep, dt, mxlay_above_rock, mxlay_above_wt, implicit_contribution, surf_CH4_conc, &
                      CH4_source, kHCC, Di, D_CH4_surf, layer_thick, dens2conc, conc2dens, CH4_dens, surf_CH4_flux,lissue_tst)
    lissue_loc = lissue_loc .OR. lissue_tst

    ! diffusion calculation: CO2

    ! set up diffusion coefficients
    DO is = 1, mxlay_above_rock
      Di  (is) = D_CO2(is)
      kHCC(is) = determine_Henry_cc(CO2_consts, Tsoil(is))
    ENDDO

    CALL diffusion_cn(nstep, dt, mxlay_above_rock, mxlay_above_wt, implicit_contribution, surf_CO2_conc, &
                      CO2_source, kHCC, Di, D_CO2_surf, layer_thick, dens2conc, conc2dens, CO2_dens, surf_CO2_flux,lissue_tst)
    lissue_loc = lissue_loc .OR. lissue_tst

    ! diffusion calculation: O2

    ! set up diffusion coefficients
    DO is = 1, mxlay_above_rock
      Di  (is) = D_O2(is)
      kHCC(is) = determine_Henry_cc(O2_consts, Tsoil(is))
    ENDDO

    CALL diffusion_cn(nstep, dt, mxlay_above_rock, mxlay_above_wt, implicit_contribution, surf_O2_conc, &
                      O2_source, kHCC, Di, D_O2_surf, layer_thick, dens2conc, conc2dens, O2_dens, surf_O2_flux,lissue_tst)
    lissue_loc = lissue_loc .OR. lissue_tst

    IF (PRESENT(lissue)) lissue = lissue_loc
    IF (PRESENT(nissue)) nissue = nissue_loc

  END SUBROUTINE reaction_diffusion

  !!
  ! @brief      Calculates the flux of methane produced by wildfires
  !
  !!
  !!
  PURE SUBROUTINE calc_methane_fire(CO2_fire_emission, landscape_type,                &
                                    CH4_out_fire, fire_em_NO, fire_em_CO, fire_em_iso)

    USE mo_methane_constants,      ONLY: max_specietypes, ef => emission_factors_for_methane, &
                                         SPECIE_CH4, SPECIE_NOx, SPECIE_CO, SPECIE_ISOPRENE,  &
                                         molarMassCO2_kg, molarMassCH4_kg, molarMassC_kg, C_2_biomass_ratio

    REAL(wp), INTENT(in)  :: CO2_fire_emission                ! CO2 emissions from fire [kg(CO2) m-2(grid cell) s-1]
    INTEGER,  INTENT(in)  :: landscape_type                   ! Emission factors for tile veg type and various chem species
                                                              ! [g(specie) kg-1(dry matter)]
    REAL(wp), INTENT(out) :: CH4_out_fire                     ! Negative methane emission from fire [mol m-2 s-1]
    REAL(wp), INTENT(out) :: fire_em_NO                       ! NOx emission      from fire [kg m-2 s-1]
    REAL(wp), INTENT(out) :: fire_em_CO                       ! CO emission       from fire [kg m-2 s-1]
    REAL(wp), INTENT(out) :: fire_em_iso                      ! Isoprene emission from fire [kg m-2 s-1]

    CHARACTER(len=*), PARAMETER :: routine = modname//':calc_methane_fire'

    REAL(wp) :: matter_burned_kg

    ! Very likely CO2_fire_emission should not actually be interpreted as a CO2 flux, but rather as a C loss flux.
    ! Therefore convert CO2_fire_emission to C (dry biomass)
    matter_burned_kg = ((CO2_fire_emission / molarMassCO2_kg) * molarMassC_kg) / C_2_biomass_ratio
    ! this SHOULD now be kg(dry matter)/m2(grid cell)/s

    ! flip sign for JSBACH sign convention -> loss of CH4 is negative
    CH4_out_fire = ef(SPECIE_CH4,     landscape_type) * matter_burned_kg / 1000._wp / molarMassCH4_kg ! [mol(CH4)]
    fire_em_NO   = ef(SPECIE_NOx,     landscape_type) * matter_burned_kg / 1000._wp                   ! [kg(NO)]
    fire_em_CO   = ef(SPECIE_CO,      landscape_type) * matter_burned_kg / 1000._wp                   ! [kg(CO)]
    fire_em_iso  = ef(SPECIE_ISOPRENE,landscape_type) * matter_burned_kg / 1000._wp                   ! [kg(C5H8)]

  END SUBROUTINE calc_methane_fire

  ! @brief      Calculates biomass of termites and the flux of methane that the produce
  !
  !!
! GPP dep for:
!    (cover_type < 17): NaturalVegFlag .OR. PastureFlag (no crops)
!      ! tropical termite biomass / emission factor assumed for cover types 2 (tropical broadleaf evergreen),
!      ! 3 (tropical broadleaf deciduous), 10 (raingreen shrubs), 13 (C4 grass), 16 (C4 pasture)
!      IF ((cover_type== 2) .OR. (cover_type== 3) .OR. &: 
!          (cover_type==10) .OR. &                      : PhenologyType == 2 (raingreen)
!          (cover_type==13) .OR. (cover_type==16))      : C4Flag
  PURE SUBROUTINE calc_methane_termites(is_crop, frac_tile, GPP, GDD, min_mmtemp20, GPP_dep_termites, &
                                        termite_biomass, CH4_out_termite, frac_non_crop, sum_non_crop)

    USE mo_methane_constants,      ONLY: molarMassCO2_kg, molarMassCH4_kg, molarMassC_kg
    USE mo_jsb_physical_constants, ONLY: tmelt

    LOGICAL,            INTENT(in)    :: is_crop            ! Is this a crop type
    REAL(wp),           INTENT(in)    :: frac_tile          ! Fraction of tile
    REAL(wp),           INTENT(in)    :: GPP                ! GPP of last year [mol(C) m-2 s-1]
    REAL(wp),           INTENT(in)    :: GDD                ! Growing degree days (of last year) [degC d]
    REAL(wp),           INTENT(in)    :: min_mmtemp20       ! Long term minimum temperature [degC]
    LOGICAL,            INTENT(in)    :: GPP_dep_termites   ! Does termite biomass depend on GPP for this tile?
    REAL(wp),           INTENT(out)   :: termite_biomass    ! Biomass of termites [g m-2]
    REAL(wp),           INTENT(out)   :: CH4_out_termite    ! Methane production flux by termites [mol m-2 s-1]
    REAL(wp),           INTENT(inout) :: frac_non_crop      ! Fraction of area over which sum_non_crop has been summed
    REAL(wp),           INTENT(inout) :: sum_non_crop       ! CH4 production from all non-crop tiles

    CHARACTER(len=*), PARAMETER :: routine = modname//':calc_methane_termites'
    REAL(wp), PARAMETER  :: efactor_CH4_TRO = 2.8e-6_wp     ! [g CH4 / g termite / h]
    REAL(wp), PARAMETER  :: efactor_CH4_TEM = 1.7e-6_wp     ! [g CH4 / g termite / h]
    REAL(wp), PARAMETER  :: seconds_per_hour = 3600._wp

    ! locals
    REAL(wp) :: termite_CH4_g_h

    termite_biomass = 0._wp
    termite_CH4_g_h = 0._wp
    IF (.NOT. is_crop) THEN
      ! tropical termite biomass / emission factor assumed for cover types 2 (tropical broadleaf evergreen),
      ! 3 (tropical broadleaf deciduous), 10 (raingreen shrubs), 13 (C4 grass), 16 (C4 pasture)
      IF (GPP_dep_termites) THEN
        termite_biomass = 1.21_wp * EXP(0.0008_wp*GPP*molarMassC_kg*1000._wp)
        termite_CH4_g_h = termite_biomass * efactor_CH4_TRO
      ELSE
        ! non-tropical environment fixed termite biomass
        ! only if gddsum > 1200 and tmin > 3
        IF ((min_mmtemp20-tmelt > 3._wp) .AND. (GDD > 1200._wp)) THEN
          termite_biomass = 3._wp
          termite_CH4_g_h = termite_biomass * efactor_CH4_TEM
        END IF
      END IF
    END IF
    IF (is_crop) THEN
      IF (frac_non_crop > 0._wp) THEN
        !  Crops are 40% of non-crop mean
        termite_CH4_g_h = sum_non_crop / frac_non_crop * 0.4_wp
      ENDIF
    ELSE
      frac_non_crop = frac_non_crop + frac_tile
      sum_non_crop = sum_non_crop + termite_CH4_g_h * frac_tile
    ENDIF

    CH4_out_termite = -termite_CH4_g_h / 1000._wp / molarMassCH4_kg / seconds_per_hour ! flip sign -> loss of CH4 is negative

  END SUBROUTINE calc_methane_termites

  PURE SUBROUTINE calc_methane_herbivores(herbivory2atm, iswoody, CH4_out_herbivore)

    USE mo_methane_constants, ONLY: molarMassCH4_kg, molarMassC_kg, C_2_biomass_ratio

    REAL(wp), INTENT(in)  :: herbivory2atm                 ! herbivory C-flux to atmosphere [mol(C) m-2 s-1]
    LOGICAL,  INTENT(in)  :: iswoody                       ! Is type woody (forest/shrub)
    REAL(wp), INTENT(out) :: CH4_out_herbivore             ! methane production from herbivores

    ! Parameters
    CHARACTER(len=*), PARAMETER :: routine = modname//':calc_methane_herbivores'
    REAL(wp), PARAMETER  :: CH4prod_mammal    = 14.9_wp   ! mean species value by Clauss et al., 2020
                                                          ! [g(CH4) / kg (DM)] We assume: DM = biomass
    REAL(wp), PARAMETER  :: mammalfrac_forest =  0.02_wp  ! herbivory fraction by mammals in forest
    REAL(wp), PARAMETER  :: mammalfrac_grass  =  0.4_wp   ! herbivory fraction by mammals in grassland

    ! locals
    REAL(wp)            :: herbivore_methane              ! methane flux from herbivores
    REAL(wp)            :: mammalfrac

    IF (iswoody) THEN
      mammalfrac = mammalfrac_forest
    ELSE
      mammalfrac = mammalfrac_grass
    ENDIF

    herbivore_methane = herbivory2atm * molarMassC_kg / C_2_biomass_ratio * mammalfrac * CH4prod_mammal
    CH4_out_herbivore = -herbivore_methane / 1000._wp / molarMassCH4_kg   ! convert to mol and switch sign -> loss of CH4 is < 0

  END SUBROUTINE calc_methane_herbivores

  !!
  ! @brief     determine distribution of plant roots in model soil layers (based on Jackson 1996)
  !
  !!
  PURE SUBROUTINE determine_root_dist(mxlay_above_rock, layer_thick, layer_depth, root_dist_beta, root_dist)
    INTEGER,  INTENT(in)  :: mxlay_above_rock     ! layer number above bedrock
    REAL(wp), INTENT(in)  :: layer_thick(:)       ! soil layer thickness
    REAL(wp), INTENT(in)  :: layer_depth(:)       ! soil layer depth
    REAL(wp), INTENT(in)  :: root_dist_beta       ! exponential factor for root distribution
    REAL(wp), INTENT(out) :: root_dist(:)         ! root distribution in each layer

    CHARACTER(len=*), PARAMETER :: routine = modname//':determine_root_dist'

    ! locals
    REAL(wp)  :: root_dist_cum(mxlay_above_rock) ! cumulative root distribution for layers
    INTEGER  :: iLayer                                              ! layer counter

    root_dist(:) = 0._wp

    IF (root_dist_beta < 0._wp) THEN
      root_dist(1:mxlay_above_rock) = layer_thick(1:mxlay_above_rock)
    ELSE
      DO iLayer = 1, mxlay_above_rock
        root_dist_cum(iLayer) = 1._wp - (root_dist_beta ** (layer_depth(iLayer) * 100._wp))
        IF (iLayer == 1) THEN
          root_dist(iLayer) = root_dist_cum(iLayer)
        ELSE
          root_dist(iLayer) = root_dist_cum(iLayer) - root_dist_cum(iLayer-1)
        END IF
      END DO
    END IF
    IF (SUM(root_dist) > 0._wp) THEN
      root_dist(:) = root_dist(:) / SUM(root_dist(:))
    ELSE
      root_dist(:) = 0._wp
    END IF

  END SUBROUTINE determine_root_dist

#endif
END MODULE mo_methane_process
