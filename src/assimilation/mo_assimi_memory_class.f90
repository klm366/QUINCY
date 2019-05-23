!>
!! @brief Contains the memory class for the assimilation process.
!!
!! @author
!!  Rainer Schneck, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Rainer Schneck (2015-05-8)
!!
MODULE mo_assimi_memory_class
#ifndef __NO_JSBACH__

  USE mo_kind, ONLY: wp
  !USE mo_util, ONLY: One_of

  !USE mo_jsb_model_class,        ONLY: t_jsb_model
  !USE mo_jsb_class,              ONLY: Get_model
  USE mo_jsb_memory_class,       ONLY: t_jsb_memory
  !USE mo_jsb_lct_class,          ONLY: VEG_TYPE
  USE mo_jsb_var_class,          ONLY: t_jsb_var_real2d, t_jsb_var_real3d

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: t_assimi_memory, max_no_of_vars

  INTEGER, PARAMETER :: max_no_of_vars = 40

  !> Type definition for memory
  TYPE, EXTENDS(t_jsb_memory) :: t_assimi_memory
     ! t_canopy_sum variables
    TYPE(t_jsb_var_real2d) ::     &
      gross_assimilation_ca,      & ! Gross photosynthesis at each time step
                                    ! [mol(co2) / ( m^2(canopy ground) * s) ]  (mean value)
                                    ! Note, this is the gross assimilation on the ground below canopy area.
                                    ! This area is independent of the LAI. So even when the LAI is below 1 this area remains
                                    ! the same (it is still "covered").
      gross_assimilation,         & ! ...on tile area
      dark_respiration_ca,        & ! Dark respiration of leaf at each time step
                                    ! [mol(co2) / (m^2(canopy ground) * s) ]  (mean value)
      dark_respiration,           & ! see above but here on tile area
      NPP_pot_rate_ca,            & ! Net primary production rate [mol(CO2)/(m^2 (canopy ground) * s)]
      canopy_cond_over_time

!                                         CO2_conc_leaf,                  & ! CO2 concentration inside leaf [MOL(CO2)/MOL(AIR)]
!                                         CO2_conc_leaf_acc,              & ! CO2 concentration inside leaf [MOL(CO2)/MOL(AIR)]
!                                                                           ! (mean value)
!                                         CO2_conc_leaf_unlimited_acc,    & ! CO2 concentration inside leaf for no water stress
!                                                                           ! [MOL(CO2)/MOL(AIR)] (mean value)
!
!                                         carbox_rate_max,                & ! Maximum carboxylation rate (=VCmax)
!                                                                           ! [MOL(CO2)/M^2 S]
!                                         e_transport_rate_max,           & ! Maximum electron transport rate (=Jmax)
!                                                                           ! (only C3 plants)[MOL(CO2)/M^2 S]
!                                         carbox_rate,                    & ! Actual carboxylation rate (=JC)
!                                                                           ! [MOL(CO2)/M^2 S]
!                                         e_transport_rate,               & ! Actual electron transport rate (=JE)
!                                                                           ! [MOL(CO2)/M^2 S]
!                                         net_assimilation,               & ! Canopy net assimilation assimilation
!                                                                           ! [MOL(CO2)/m^2 s] (mean value)

!
!      REAL(wp), POINTER, DIMENSION(:,:) ::                               &
!                                         CO2_flux_net_acc,               & ! net CO2 flux to the atmosphere
!                                         CO2_flux_herbivory_acc,         & ! CO2 flux to the atmosphere from grazing
!                                         CO2_emission_landcover_change_acc, & ! CO2 flux to the atmosphere from land cover changes
!                                         CO2_emission_harvest_acc,       & ! CO2 flux to the atmosphere from harvest
!                                         CO2_flux_dynveg_acc               ! CO2 flux to the atmosphere from vegetation dynamics
!                                                                           ! (disturbance?)
!


    ! t_canopy variables
    TYPE(t_jsb_var_real3d) ::     &
      & assimi_scaling_fact_cl,   &
      & canopy_cond_cl,           & ! Canopy conductance per layer and per leaf area
      & gross_assimilation_cl,    & ! Gross photosynthesis [ mol(co2) / (m^2(leaf area) * s) ]
                                    ! Note, this is the gross assimilation on the leaf area that actually exists.
                                    ! The influence of the lai on the existing leaf area is aready included.
                                    ! It is not the ground area below canopy!
      & dark_respiration_cl,      & ! Dark respiration of leaf [ mol(co2) / (m^2(leaf area) * s) ]
      & carbox_rate_max_cl,       & ! [ mol(co2)/(m^2(leaf area) * s) ]
      & e_transport_rate_max_cl,  &
      & carbox_rate_cl,           &
      & e_transport_rate_cl,      &
      & CO2_conc_leaf_cl            ! CO2 concentration inside leaf [ mol(CO2)/ mol(AIR) ]

    TYPE(t_jsb_var_real2d) :: &
      & canopy_cond                 ! Canopy conductance without water limit  [m/s] (mean value)

    ! R: This variable is needed for process DISTURB and should be later in the process vegdyn
    !    However it has to exist on each pft, therefore I put it to this process for now!
    TYPE(t_jsb_var_real2d) :: &
      & cover_fract_pot             ! Potential Natural Land Cover Fraction

  CONTAINS
    PROCEDURE :: Init => Init_assimi_memory
  END TYPE t_assimi_memory


  CHARACTER(len=*), PARAMETER :: modname = 'mo_assimi_memory_class'

CONTAINS

  SUBROUTINE Init_assimi_memory(mem, prefix, suffix, lct_ids, model_id)

    USE mo_jsb_varlist,       ONLY: jsb_add_var, BASIC !, MEDIUM, FULL
    USE mo_jsb_io,            ONLY: grib_bits, t_cf, t_grib1, t_grib2, tables !, TSTEP_CONSTANT
    USE mo_jsb_grid_class,    ONLY: t_jsb_grid, t_jsb_vgrid
    USE mo_jsb_grid,          ONLY: Get_grid, Get_vgrid

    CLASS(t_assimi_memory), INTENT(inout), TARGET :: mem
    CHARACTER(len=*),    INTENT(in)    :: prefix
    CHARACTER(len=*),    INTENT(in)    :: suffix
    INTEGER,             INTENT(in)    :: lct_ids(:)
    INTEGER,             INTENT(in)    :: model_id

    TYPE(t_jsb_grid),  POINTER :: hgrid                        ! Horizontal grid
    TYPE(t_jsb_vgrid), POINTER :: surface, canopy_layer        ! Vertical grid
    INTEGER :: table

    CHARACTER(len=*), PARAMETER :: routine = modname//':Init_assimi_memory'

    IF (lct_ids(1) > 0)  CONTINUE ! avoid compiler warning about dummy argument not being used
    IF (model_id > 0)    CONTINUE ! avoid compiler warning about dummy argument not being used

    table = tables(1)

    hgrid   => Get_grid(mem%grid_id)
    surface => Get_vgrid('surface')
    canopy_layer => Get_vgrid('canopy_layer')

    !$acc enter data pcreate(mem)

    IF (.NOT. ASSOCIATED(mem%gross_assimilation_ca%ptr)) THEN

      CALL jsb_add_var( mem, 'gross_assimilation_ca', mem%gross_assimilation_ca,     &
        & hgrid, surface,                                                            &
        & t_cf('gross_assimilation_ca', 'mol(CO2)/m^2(canopy ground) / s',           &
        &      'Gross photosynthesis on canopy ground area'),                        &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),         &
        & prefix, suffix,                                                            &
        & output_level=BASIC,                                                        &
        & lrestart=.FALSE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

      CALL jsb_add_var( mem, 'gross_assimilation', mem%gross_assimilation,                                          &
        & hgrid, surface,                                                                                           &
        & t_cf('gross_assimilation', 'mol(CO2)/m^2 (tile area) / s', 'Gross photosynthesis on tile area'),          &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                        &
        & prefix, suffix,                                                                                           &
        & output_level=BASIC,                                                                                       &
        & lrestart=.FALSE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

      CALL jsb_add_var( mem, 'dark_respiration_ca', mem%dark_respiration_ca,                                        &
        & hgrid, surface,                                                                                           &
        & t_cf('dark_respiration_ca', 'mol(CO2)/m^2(canopy ground) / s', 'Dark_respiration on canopy ground area'), &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                        &
        & prefix, suffix,                                                                                           &
        & output_level=BASIC,                                                                                       &
        & lrestart=.FALSE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

      CALL jsb_add_var( mem, 'dark_respiration', mem%dark_respiration,                                              &
        & hgrid, surface,                                                                                           &
        & t_cf('dark_respiration', 'mol(CO2)/m^2 (tile area) / s', 'Dark_respiration on tile area'),                &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                        &
        & prefix, suffix,                                                                                           &
        & output_level=BASIC,                                                                                       &
        & lrestart=.FALSE., initval_r=0.0_wp, l_aggregate_all=.TRUE. ) ! initval taken from JSBACH3

      CALL jsb_add_var( mem, 'npp_pot_rate_ca', mem%NPP_pot_rate_ca,                        &
        & hgrid, surface,                                                                   &
        & t_cf('NPP_pot_rate_ca', 'mol(CO2)/(m^2 (canopy ground) / s)', 'NPP_pot_rate_ca'), &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                &
        & prefix, suffix,                                                                   &
        & output_level=BASIC,                                                               &
        & lrestart=.FALSE., initval_r=0.0_wp, l_aggregate_all=.TRUE. )

      CALL jsb_add_var( mem, 'canopy_cond_over_time', mem%canopy_cond_over_time, &
        & hgrid, surface,                                                        &
        & t_cf('canopy_cond_over_time', '-',                                     &
        'Canopy cond. integrated over time'),                                    &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),     &
        & prefix, suffix,                                                        &
        & lrestart=.FALSE., initval_r=0.0_wp )

      ! REAL(wp), POINTER, DIMENSION(:,:,:) ::
      CALL jsb_add_var( mem, 'scaling_fact_cl', mem%assimi_scaling_fact_cl,  &
        & hgrid, canopy_layer,                                               &
        & t_cf('assimi_scaling_fact_cl', '-',                                &
        'assimilation scaling factors per canopy layer'),                    &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix,                                                    &
        & lrestart=.FALSE., initval_r=0.0_wp ) ! initval taken from JSBACH3

      CALL jsb_add_var( mem, 'canopy_cond_cl', mem%canopy_cond_cl,           &
        & hgrid, canopy_layer,                                               &
        & t_cf('canopy_cond_cl', '-',                                        &
        'Canopy cond. per canopy layer and per leaf area'),                  &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix,                                                    &
        & lrestart=.FALSE., initval_r=0.0_wp ) ! initval taken from JSBACH3

      CALL jsb_add_var( mem, 'gross_assimilation_cl', mem%gross_assimilation_cl,                                              &
        & hgrid, canopy_layer,                                                                                                &
        & t_cf('gross_assimilation_cl', 'mol(CO2)/m^2(leaf area) / s', 'Gross photosynthesis per canopy layer on leaf area'), &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                                  &
        & prefix, suffix,                                                                                                     &
        & lrestart=.FALSE., initval_r=0.0_wp, l_aggregate_all=.TRUE.  )

      CALL jsb_add_var( mem, 'dark_respiration_cl', mem%dark_respiration_cl,                                                    &
        & hgrid, canopy_layer,                                                                                                  &
        & t_cf('dark_respiration_cl', 'mol(CO2)/m^2(leaf area) / s', 'Dark respiration of leaf per canopy layer on leaf area'), &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                                    &
        & prefix, suffix,                                                                                                       &
        & lrestart=.FALSE., initval_r=0.0_wp, l_aggregate_all=.TRUE.  )

      CALL jsb_add_var( mem, 'carbox_rate_max_cl', mem%carbox_rate_max_cl,                                          &
        & hgrid, canopy_layer,                                                                                      &
        & t_cf('carbox_rate_max_cl', 'mol(CO2)/m^2(leaf area) / s', 'Maximum carboxylation rate per canopy layer'), &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                        &
        & prefix, suffix,                                                                                           &
        & lrestart=.FALSE., initval_r=0.0_wp, l_aggregate_all=.TRUE.  )

      CALL jsb_add_var( mem, 'e_transport_rate_max_cl', mem%e_transport_rate_max_cl, &
        & hgrid, canopy_layer,                                                       &
        & t_cf('e_transport_rate_max_cl', 'dimensionen',                             &
        'Maximum electron transport rate per canopy layer'),                         &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),         &
        & prefix, suffix,                                                            &
        & lrestart=.FALSE., initval_r=0.0_wp )

      CALL jsb_add_var( mem, 'carbox_rate_cl', mem%carbox_rate_cl,                                             &
        & hgrid, canopy_layer,                                                                                 &
        & t_cf('carbox_rate_cl', 'mol(CO2)/m^2(leaf area) / s', 'Actual carboxylation rate per canopy layer'), &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                   &
        & prefix, suffix,                                                                                      &
        & lrestart=.FALSE., initval_r=0.0_wp, l_aggregate_all=.TRUE.  )

      CALL jsb_add_var( mem, 'e_transport_rate_cl', mem%e_transport_rate_cl,                                             &
        & hgrid, canopy_layer,                                                                                           &
        & t_cf('e_transport_rate_cl', 'mol(CO2)/m^2(leaf area) / s', 'Actual electron transport rate per canopy layer'), &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                             &
        & prefix, suffix,                                                                                                &
        & lrestart=.FALSE., initval_r=0.0_wp, l_aggregate_all=.TRUE.  )! initval taken from JSBACH3

      CALL jsb_add_var( mem, 'co2_conc_leaf_cl', mem%CO2_conc_leaf_cl,                               &
        & hgrid, canopy_layer,                                                                       &
        & t_cf('CO2_conc_leaf_cl', 'dimensionen', 'CO2 concentration inside leaf per canopy layer'), &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                         &
        & prefix, suffix,                                                                            &
        & lrestart=.FALSE., initval_r=0.0_wp ) ! initval taken from JSBACH3 CO2_conc_leaf_acc

      CALL jsb_add_var( mem, 'canopy_cond', mem%canopy_cond,                 &
        & hgrid, surface,                                                    &
        & t_cf('canopy_conductance', '', ''),                                &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix,                                                    &
        & initval_r=0.0_wp )

      CALL jsb_add_var( mem, 'cover_fract_pot', mem%cover_fract_pot,         &
        & hgrid, surface,                                                    &
        & t_cf('cover_fract_pot', 'fract', ''),                              &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix,                                                    &
        & initval_r=0.1_wp ) ! initval choosen freely


    END IF


  END SUBROUTINE Init_assimi_memory

#endif
END MODULE mo_assimi_memory_class
