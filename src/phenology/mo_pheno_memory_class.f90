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
!! @brief Contains the memory class for the phenology process.
!!
!! @author
!!  Jan Engel, MPI-BGC Jena (based on the JSBACH4 version from Rainer Schneck, MPI-M Hamburg)
!!
!! @par Revision History
!! First version by Jan Engel (20161214)
!!
MODULE mo_pheno_memory_class
#ifndef __NO_JSBACH__

  USE mo_kind,                   ONLY: wp
  USE mo_exception,              ONLY: message
  USE mo_jsb_class,              ONLY: Get_model
  USE mo_jsb_memory_class,       ONLY: t_jsb_memory
  USE mo_jsb_var_class,          ONLY: t_jsb_var_real2d, t_jsb_var_real3d


  IMPLICIT NONE
  PRIVATE
  PUBLIC :: t_pheno_memory, max_no_of_vars

  INTEGER, PARAMETER :: max_no_of_vars = 20   ! default JSBACH4


  !-----------------------------------------------------------------------------------------------------
  !> Type definition for pheno memory
  !! 
  !-----------------------------------------------------------------------------------------------------
  TYPE, EXTENDS(t_jsb_memory) :: t_pheno_memory


    ! used as logical: 0 = FALSE / 1 = TRUE
    TYPE(t_jsb_var_real2d)          :: &
                                    growing_season          !< logical to identify whether the plant is in the growing season or not


    TYPE(t_jsb_var_real2d)          :: &
                                    gdd                 , & !< growing degree days above 5°C [K day]
                                    nd_dormance         , & !< number of days since last growing season [# days]
                                    lai_max                 !< maximum lai - from site specific lai data [m2 m-2]


  CONTAINS
    PROCEDURE :: Init => Init_pheno_memory
  END TYPE t_pheno_memory


  CHARACTER(len=*), PARAMETER :: modname = 'mo_pheno_memory_class'


CONTAINS


  !-----------------------------------------------------------------------------------------------------
  !> initialize memory (variables) for the process: phenology
  !!
  !!  
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE Init_pheno_memory(mem, prefix, suffix, lct_ids, model_id)

    USE mo_jsb_varlist,       ONLY: BASIC !, MEDIUM, FULL
    USE mo_jsb_io,            ONLY: grib_bits, t_cf, t_grib1, t_grib2   ! , TSTEP_CONSTANT, tables
    USE mo_jsb_grid_class,    ONLY: t_jsb_grid, t_jsb_vgrid
    USE mo_jsb_grid,          ONLY: Get_grid, Get_vgrid
    USE mo_jsb_model_class,   ONLY: t_jsb_model
    ! ----------------------------------------------------------------------------------------------------- !
    CLASS(t_pheno_memory), INTENT(inout), TARGET  :: mem           !> pheno memory
    CHARACTER(len=*),     INTENT(in)            :: prefix          !> process name
    CHARACTER(len=*),     INTENT(in)            :: suffix          !> tile name
    INTEGER,              INTENT(in)            :: lct_ids !(:)    !> in jsbach4: vector of land cover type IDs
    INTEGER,              INTENT(in)            :: model_id        !> model ID model\%id
    ! ----------------------------------------------------------------------------------------------------- !
    TYPE(t_jsb_model), POINTER :: model
    TYPE(t_jsb_grid),  POINTER :: hgrid                        ! Horizontal grid
    TYPE(t_jsb_vgrid), POINTER :: surface                      ! Vertical grid
    INTEGER                    :: table
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':Init_pheno_memory'
    ! ----------------------------------------------------------------------------------------------------- !
    model         => Get_model(model_id)
    table         = 1     ! in jsbach4 this line is: table = tables(1)     *slmtemp
    hgrid         => Get_grid(model%grid_id)
    surface       => Get_vgrid('surface')

    
    ! ----------------------------------------------------------------------------------------------------- !
    ! quincy
    ! IF (model%config%use_quincy) THEN
    !   IF ( One_of(LAND_TYPE, lct_ids(:)) > 0 .OR. &
    !      & One_of(VEG_TYPE,  lct_ids(:)) > 0) THEN
    
        CALL mem%Add_var('growing_season', mem%growing_season, &
          & hgrid, surface, &
          & t_cf('growing_season', 'logical', 'logical to identify whether the plant is in the growing season'), &
          & t_grib1(table, 255, grib_bits), &
          & t_grib2(255, 255, 255, grib_bits), &
          & prefix, suffix, &
          & loutput = .TRUE., &
          & lrestart = .TRUE., &
          & initval_r = 0.0_wp)
  
        CALL mem%Add_var('gdd', mem%gdd, &
          & hgrid, surface, &
          & t_cf('gdd', 'K days', 'growing degree days above 5 degC'), &
          & t_grib1(table, 255, grib_bits), &
          & t_grib2(255, 255, 255, grib_bits), &
          & prefix, suffix, &
          & loutput = .TRUE., &
          & lrestart = .TRUE., &
          & initval_r = 0.0_wp)
    
        CALL mem%Add_var('nd_dormance', mem%nd_dormance, &
          & hgrid, surface, &
          & t_cf('nd_dormance', '# days', 'number of days since last growing season'), &
          & t_grib1(table, 255, grib_bits), &
          & t_grib2(255, 255, 255, grib_bits), &
          & prefix, suffix, &
          & loutput = .TRUE., &
          & lrestart = .TRUE., &
          & initval_r = 0.0_wp)
  
        CALL mem%Add_var('lai_max', mem%lai_max, &
          & hgrid, surface, &
          & t_cf('lai_max', 'm2 m-2', 'maximum lai (site specific)'), &
          & t_grib1(table, 255, grib_bits), &
          & t_grib2(255, 255, 255, grib_bits), &
          & prefix, suffix, &
          & loutput = .TRUE., &
          & lrestart = .TRUE., &
          & initval_r = 0.0_wp)

    !   END IF
    ! END IF

  END SUBROUTINE Init_pheno_memory

#endif
END MODULE mo_pheno_memory_class
