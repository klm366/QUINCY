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
!! @brief Contains the memory class for the soil_biogeochemistry processes
!!
!! @author
!!  Jan Engel, MPI-BGC Jena (based on the JSBACH4 version from Rainer Schneck, MPI-M Hamburg)
!!
!! @par Revision History
!! First version by Jan Engel (20170103)
!!
MODULE mo_sb_memory_class

  USE mo_kind,                   ONLY: wp
  USE mo_exception,              ONLY: finish
  USE mo_util,                   ONLY: One_of
  USE mo_jsb_utils_iface,        ONLY: assign_if_present_allocatable
  USE mo_jsb_impl_constants,     ONLY: false
  
  USE mo_jsb_model_class,        ONLY: t_jsb_model
  USE mo_jsb_class,              ONLY: Get_model
  USE mo_jsb_memory_class,       ONLY: t_jsb_memory
  USE mo_jsb_var_class,          ONLY: t_jsb_var, t_jsb_var_real2d, t_jsb_var_real3d
  USE mo_jsb_pool_class,         ONLY: t_jsb_pool, ELEM_C, ELEM_N, ELEM_P, ELEM_C13, ELEM_C14, ELEM_N15
  USE mo_jsb_lct_class,          ONLY: LAND_TYPE, VEG_TYPE
  

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: t_sb_memory, max_no_of_vars

  INTEGER, PARAMETER :: max_no_of_vars      = 200  !> TYPE t_jsb_var 


  !-----------------------------------------------------------------------------------------------------
  !> Type definition for sb memory
  !!
  !! 
  !-----------------------------------------------------------------------------------------------------
  TYPE, EXTENDS(t_jsb_memory) :: t_sb_memory


    ! 2D total fluxes for mass balance checks in soil_biogeochemistry model
    TYPE(t_jsb_var_real2d)          :: & 
                                    total_flux_carbon, &          !< carbon mass balance checks in soil_biogeochemistry model [mol C m-2 s-1]
                                    total_flux_nitrogen, &        !< nitrogen mass balance checks in soil_biogeochemistry model [mol N m-2 s-1]
                                    total_flux_phosphorus, &      !< phosphorus mass balance checks in soil_biogeochemistry model [mol P m-2 s-1]
                                    total_flux_carbon13, &        !< carbon13 mass balance checks in soil_biogeochemistry model [mol C13 m-2 s-1]
                                    total_flux_carbon14, &        !< carbon14 mass balance checks in soil_biogeochemistry model [mol C14 m-2 s-1]
                                    total_flux_nitrogen15         !< nitrogen15 mass balance checks in soil_biogeochemistry model [mol N15 m-2 s-1]


  CONTAINS
    PROCEDURE :: Init => Init_sb_memory
  END TYPE t_sb_memory

  CHARACTER(len=*), PARAMETER :: modname = 'mo_sb_memory_class'

CONTAINS

  !-----------------------------------------------------------------------------------------------------
  !> initialize memory for the SB_ process
  !!
  !!  
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE Init_sb_memory(mem, prefix, suffix, lct_ids, model_id)

    USE mo_jsb_varlist,       ONLY: BASIC, MEDIUM, FULL
    USE mo_jsb_io,            ONLY: grib_bits, t_cf, t_grib1, t_grib2, tables
    USE mo_jsb_grid_class,    ONLY: t_jsb_grid, t_jsb_vgrid
    USE mo_jsb_grid,          ONLY: Get_grid, Get_vgrid
    USE mo_jsb_model_class,   ONLY: t_jsb_model
    USE mo_jsb_class,         ONLY: Get_model

    CLASS(t_sb_memory),   INTENT(inout), TARGET  :: mem             !> SB_ memory
    CHARACTER(len=*),     INTENT(in)             :: prefix          !> process name
    CHARACTER(len=*),     INTENT(in)             :: suffix          !> tile name
    INTEGER,              INTENT(in)             :: lct_ids(:)      !> vector of land cover type IDs
    INTEGER,              INTENT(in)             :: model_id        !> model ID model\%id
    ! ---------------------------
    TYPE(t_jsb_model), POINTER  :: model
    TYPE(t_jsb_grid),  POINTER  :: hgrid                        ! Horizontal grid
    TYPE(t_jsb_vgrid), POINTER  :: vgrid_surface                ! Vertical grid
!    TYPE(t_jsb_vgrid), POINTER  :: vgrid_soil_sb                ! Vertical grid
    INTEGER                     :: table
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':Init_sb_memory'
    ! ---------------------------
    model         => Get_model(model_id)
    table         =  tables(1)
    hgrid         => Get_grid(model%grid_id)
    vgrid_surface => Get_vgrid('surface')
!    vgrid_soil_sb => Get_vgrid('soil_layer_sb')


    !> add memory only for LAND & PFT tiles
    !!
    IF ( One_of(LAND_TYPE, lct_ids(:)) > 0 .OR. &
       & One_of(VEG_TYPE,  lct_ids(:)) > 0) THEN

      !=============================================================================================================================
      !> pools & fluxes
      !=============================================================================================================================





      !=============================================================================================================================
      !> 2D & 3D variables
      !=============================================================================================================================

      CALL mem%Add_var('total_flux_carbon', mem%total_flux_carbon, &
            hgrid, vgrid_surface, &
            t_cf('total_flux_carbon', 'mol C m-2 s-1', 'carbon mass balance checks in soil_biogeochemistry model'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            loutput=.TRUE., &
            initval_r=0.0_wp) ! initvalue

      CALL mem%Add_var('total_flux_nitrogen', mem%total_flux_nitrogen, &
            hgrid, vgrid_surface, &
            t_cf('total_flux_nitrogen', 'mol N m-2 s-1', 'nitrogen mass balance checks in soil_biogeochemistry model'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            loutput=.TRUE., &
            initval_r=0.0_wp) ! initvalue

      CALL mem%Add_var('total_flux_phosphorus', mem%total_flux_phosphorus, &
            hgrid, vgrid_surface, &
            t_cf('total_flux_phosphorus', 'mol P m-2 s-1', 'phosphorus mass balance checks in soil_biogeochemistry model'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            loutput=.TRUE., &
            initval_r=0.0_wp) ! initvalue

      CALL mem%Add_var('total_flux_carbon13', mem%total_flux_carbon13, &
            hgrid, vgrid_surface, &
            t_cf('total_flux_carbon13', 'mol C13 m-2 s-1', 'carbon13 mass balance checks in soil_biogeochemistry model'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            loutput=.TRUE., &
            initval_r=0.0_wp) ! initvalue

      CALL mem%Add_var('total_flux_carbon14', mem%total_flux_carbon14, &
            hgrid, vgrid_surface, &
            t_cf('total_flux_carbon14', 'mol C14 m-2 s-1', 'carbon14 mass balance checks in soil_biogeochemistry model'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            loutput=.TRUE., &
            initval_r=0.0_wp) ! initvalue

      CALL mem%Add_var('total_flux_nitrogen15', mem%total_flux_nitrogen15, &
            hgrid, vgrid_surface, &
            t_cf('total_flux_nitrogen15', 'mol N15 m-2 s-1', 'nitrogen15 mass balance checks in soil_biogeochemistry model'), &
            t_grib1(table, 255, grib_bits), &
            t_grib2(255, 255, 255, grib_bits), &
            prefix, suffix, &
            loutput=.TRUE., &
            initval_r=0.0_wp) ! initvalue
    ENDIF

  END SUBROUTINE Init_sb_memory

END MODULE mo_sb_memory_class
