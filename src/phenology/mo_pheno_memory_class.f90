!>
!! @brief Contains the memory class for phenology.
!!
!! @author
!!  Rainer Schneck, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Rainer Schneck (2015-09-09)
!!

MODULE mo_pheno_memory_class
#ifndef __NO_JSBACH__

  USE mo_kind, ONLY: wp
  USE mo_util, ONLY: One_of

  USE mo_jsb_model_class,        ONLY: t_jsb_model
  USE mo_jsb_class,              ONLY: Get_model
  USE mo_jsb_memory_class,       ONLY: t_jsb_memory
  USE mo_jsb_lct_class,          ONLY: VEG_TYPE, LAND_TYPE
  USE mo_jsb_var_class,          ONLY: t_jsb_var_real2d, t_jsb_var_real3d
  USE mo_jsb_physical_constants, ONLY: tmelt

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: t_pheno_memory, max_no_of_vars

  INTEGER, PARAMETER :: max_no_of_vars = 40

  TYPE, EXTENDS(t_jsb_memory) :: t_pheno_memory

     TYPE(t_jsb_var_real2d)           ::                   &
                                         w_soil_rel,       & ! Filling of the soil water bucket relative to MaxMoisture (=field capacity?)
                                         growth_phase_SG,  & ! For summergreen only: values are:
                                                             !   -1.0  .. during rest phase (i.e. from autumn event to spring event)
                                                             !    0.0  .. during growth phase (i.e. some weeks after spring event)
                                                             !    1.0  .. in vegetative phase (i.e. from end of growth phase to autum event)
                                                             ! Updated every day by the subroutine "update_growth_phase()".
                                         growth_phase_EG,  & ! For evergreen only: values are:
                                                             !    0.0  .. during growth phase (i.e. some weeks after spring event)
                                                             !    1.0  .. in vegetative phase (i.e. from ent of growth phase to spring event)
                                                             ! Updated every day by the subroutine "update_growth_phase()".
                                         growth_phase_CRP, & ! For extra-tropical crops only: values are:
                                                             !     -2.0  .. during rest phase of winter crops (i.e. from spring to autumn)
                                                             !     -1.0  .. during rest phase of summer crops (i.e. from autumn to spring)
                                                             !      0.0  .. during growth phase of summer crops (i.e. from spring to autumn)
                                                             !      1.0  .. during growth phase of second crop (i.e. from summer to autumn)
                                                             !      2.0  .. during growth phase of winter crops (i.e. from autumn to spring)
                                                             ! Updated every day by the subroutine "update_growth_phase()".
                                         days_since_growth_begin_EG, & ! For summergreen and evergreen only: serves to remember the number of
                                                                       ! days elapsed since the growth phase began 
                                         days_since_growth_begin_SG    ! (as real numbers to save it in a stream).

    TYPE(t_jsb_var_real2d)   ::                                 &  
                                         previous_day_NPP_rate, & ! The mean primary production rate [mol(CO2)/(m^2 s)] at each grid point
                                                                  ! at previous day, because not known for current day; updated at ..
                                                                  ! the beginning of each new day (stream) 
                                         laiMax_dyn,            & ! Dynamical maximum LAI [m^2(leaf)/m^2(ground)]-- not larger than the maximum  lai
                                                                  ! coming from JSBACH, but smaller than that value when NPP is negative for more
                                                                  ! than 1 year. (stream) 
                                         year_NPP_sum,          & ! This field sums the NPP-rate during the whole year, starting on Jan,1, in the
                                                                  ! northern hemisphere and at July,1, in the southern hemisphere. A negative value
                                                                  ! indicates that summation has not yet started.
                                         day_NPP_sum   ! sum of NPP since midnight (needed to compute previous_day_NPP_rate(:,:))
                                                       ! [mol(CO2)/m^2]


    TYPE(t_jsb_var_real2d)   ::                        & ! Fields used to calculate the growth phases
                                      heat_sum_EG,     & ! Heat sums --- needed for the phenology of summergreens, evergreens and crops.
                                      heat_sum_SG,     & ! Where heat summation has not yet started, "heat_sum" is set to a value < -1.
                                      heat_sum_CRP,    & ! For summer crops heat summation starts with "0".
                                      heat_sum_winter, & ! For winter crops heat summation starts with EPSILON(1.) and is positive in the
                                                         ! growing period (autumn to spring) and negative during summer. (stream)
                                      chill_days_EG,   & ! Number of chill days --- needed for the phenology of summer- and evergreens 
                                      chill_days_SG      ! (stream)


  CONTAINS
    PROCEDURE :: Init => Init_pheno_memory
  END TYPE t_pheno_memory

  CHARACTER(len=*), PARAMETER :: modname = 'mo_pheno_memory_class'

CONTAINS


  SUBROUTINE Init_pheno_memory(mem, prefix, suffix, lct_ids, model_id)

    ! Declarations
    USE mo_jsb_varlist,       ONLY: t_jsb_varlist, jsb_add_var
    USE mo_jsb_io,            ONLY: grib_bits, t_cf, t_grib1, t_grib2, &
                                    TSTEP_CONSTANT, tables
    USE mo_jsb_grid_class,    ONLY: t_jsb_grid, t_jsb_vgrid
    USE mo_jsb_grid,          ONLY: Get_grid, Get_vgrid

    ! Arguments
    CLASS(t_pheno_memory), INTENT(inout), TARGET :: mem
    CHARACTER(len=*),      INTENT(in)            :: prefix
    CHARACTER(len=*),      INTENT(in)            :: suffix
    INTEGER,               INTENT(in)            :: lct_ids(:)
    INTEGER,               INTENT(in)            :: model_id
    
    ! Local variables
    TYPE(t_jsb_model), POINTER :: model
    TYPE(t_jsb_grid),  POINTER :: hgrid                        ! Horizontal grid
    TYPE(t_jsb_vgrid), POINTER :: surface, canopy_layer        ! Vertical grid
    INTEGER :: table

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':Init_pheno_memory'

    model => Get_model(model_id)

    table = tables(1)

    ! Set pointers to variables in the "grid"
    hgrid        => Get_grid(mem%grid_id)
    surface      => Get_vgrid('surface')
    canopy_layer => Get_vgrid('canopy_layer')

    ! ---------------------------
    ! Go

    IF ( (     One_of(LAND_TYPE, lct_ids(:)) > 0 &
      &   .OR. One_of(VEG_TYPE,  lct_ids(:)) > 0 &
      &  )                                       &
      & .AND. .NOT. ASSOCIATED(mem%w_soil_rel%ptr)) THEN

      CALL jsb_add_var( mem,'w_soil_rel', mem%w_soil_rel,                     &
        &  hgrid, surface,                                                    &
        &  t_cf('water_soil_relative', '', ''),                               &
        &  t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix,                                                     &
        & initval_r=0.0_wp )

      CALL jsb_add_var( mem,'growth_phase_SG', mem%growth_phase_SG,            &
          & hgrid, surface,                                                    &
          & t_cf('growth_phase_SummerGreen', '', ''),                          &
          & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits), &
          & prefix, suffix,                                                    &
          initval_r=1.0_wp) ! R: changed from -1.0

      CALL jsb_add_var( mem,'growth_phase_EG', mem%growth_phase_EG,            &
          & hgrid, surface,                                                    &
          & t_cf('growth_phase_EverGreen', '', ''),                            &
          & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits), &
          & prefix, suffix,                                                    &
          initval_r=1.0_wp ) ! R: changed from -1.0 as the number also changed in the code...

      CALL jsb_add_var( mem,'growth_phase_CRP', mem%growth_phase_CRP,        &
        & hgrid, surface,                                                    &
        & t_cf('growth_phase_Crop', '', ''),                                 &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix,                                                    &
        & initval_r=-1.0_wp )

      CALL jsb_add_var( mem,'days_since_growth_begin_EG', mem%days_since_growth_begin_EG, &
        & hgrid, surface,                                                                 &
        & t_cf('days_since_growth_begin_EverGreen', '', ''),                              &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),              &
        & prefix, suffix,                                                                 &
        & initval_r=0.0_wp )

      CALL jsb_add_var( mem,'days_since_growth_begin_SG', mem%days_since_growth_begin_SG, &
        & hgrid, surface,                                                                 &
        & t_cf('days_since_growth_begin_SummerGreen', '', ''),                            &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),              &
        & prefix, suffix,                                                                 &
        & loutput=.TRUE., lrestart=.TRUE., initval_r=0.0_wp )


        
      ! 3D:
      CALL jsb_add_var( mem,'previous_day_NPP_rate', mem%previous_day_NPP_rate, &
        & hgrid, canopy_layer,                                                  &
        & t_cf('previous_day_NPP_rate', '', ''),                                &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),    &
        & prefix, suffix,                                                       &
        & loutput=.TRUE., lrestart=.TRUE.)

      CALL jsb_add_var( mem,'laiMax_dyn', mem%laiMax_dyn,                    &
        & hgrid, canopy_layer,                                               &
        & t_cf('laiMax_dyn', '', ''),                                        &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix,                                                    &
        & loutput=.TRUE., lrestart=.TRUE., initval_r=-1.0_wp)

      CALL jsb_add_var( mem,'year_NPP_sum', mem%year_NPP_sum,                 &
        & hgrid, canopy_layer,                                                &
        & t_cf('year_NPP_sum', '', ''),                                       &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),  &
        & prefix, suffix,                                                     &
        & loutput=.TRUE., lrestart=.TRUE., initval_r=2.0_wp)
        
      CALL jsb_add_var( mem,'day_NPP_sum', mem%day_NPP_sum,                  &
        & hgrid, canopy_layer,                                               &
        & t_cf('day_NPP_sum', '', ''),                                       &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix,                                                    &
        & loutput=.TRUE., lrestart=.TRUE.)

        
      ! 2D:
      CALL jsb_add_var( mem,'heat_sum_EG', mem%heat_sum_EG,                   &
        & hgrid, surface,                                                    &
        & t_cf('heat_sum_EG', '', ''),                                       &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix,                                                    &
        & loutput=.TRUE., lrestart=.TRUE., initval_r=-99.0_wp)
        
      CALL jsb_add_var( mem,'heat_sum_SG', mem%heat_sum_SG,                   &
        & hgrid, surface,                                                    &
        & t_cf('heat_sum_SG', '', ''),                                       &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix,                                                    &
        & loutput=.TRUE., lrestart=.TRUE., initval_r=-99.0_wp) ! has to be lower than -1
        
      CALL jsb_add_var( mem,'heat_sum_CRP', mem%heat_sum_CRP,                 &
        & hgrid, surface,                                                    &
        & t_cf('heat_sum_CRP', '', ''),                                      &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix,                                                    &
        & loutput=.TRUE., lrestart=.TRUE., initval_r=9999.0_wp)
        
      CALL jsb_add_var( mem,'heat_sum_winter', mem%heat_sum_winter,           &
        & hgrid, surface,                                                    &
        & t_cf('heat_sum_winter', '', ''),                                   &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix,                                                    &
        & loutput=.TRUE., lrestart=.TRUE., initval_r=-EPSILON(1._wp))
        
      CALL jsb_add_var( mem,'chill_days_EG', mem%chill_days_EG,               &
        & hgrid, surface,                                                    &
        & t_cf('chill_days_EG', '', ''),                                     &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix,                                                    &
        & loutput=.TRUE., lrestart=.TRUE.)
        
      CALL jsb_add_var( mem,'chill_days_SG', mem%chill_days_SG,               &
        & hgrid, surface,                                                    &
        & t_cf('chill_days_SG', '', ''),                                     &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix,                                                    &
        & loutput=.TRUE., lrestart=.TRUE.)
        
    END IF

  END SUBROUTINE Init_pheno_memory


#endif
END MODULE mo_pheno_memory_class
