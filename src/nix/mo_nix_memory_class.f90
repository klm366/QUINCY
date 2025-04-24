!> Contains the memory class for nix.
!>
!> ICON-Land
!>
!> ---------------------------------------
!> Copyright (C) 2013-2024, MPI-M, MPI-BGC
!>
!> Contact: icon-model.org
!> Authors: AUTHORS.md
!> See LICENSES/ for license information
!> SPDX-License-Identifier: BSD-3-Clause
!> ---------------------------------------
!>
MODULE mo_nix_memory_class
#ifndef __NO_JSBACH__

  USE mo_kind,             ONLY: wp
  USE mo_util,             ONLY: One_of
  USE mo_jsb_memory_class, ONLY: t_jsb_memory
  USE mo_jsb_var_class,    ONLY: t_jsb_var_real2d, t_jsb_var_real3d
  USE mo_jsb_lct_class,    ONLY: BARE_TYPE, VEG_TYPE, LAND_TYPE, GLACIER_TYPE, LAKE_TYPE

  ! Use of prcesses in this module
  dsl4jsb_Use_processes NIX_

  ! Use process configurations
  dsl4jsb_Use_config(NIX_)

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: t_nix_memory, max_no_of_vars

  INTEGER, PARAMETER :: max_no_of_vars = 6

  TYPE, EXTENDS(t_jsb_memory) :: t_nix_memory

 ! Main prognostic variables describing the state of the snow cover - Note: These are also input/output therefore  mandatory
    TYPE(t_jsb_var_real3d)  :: &
      &  t_sn                , & ! snow element temperature             [ K ]
      &  t_sn_n              , & ! snow nodal temperature               [ K ]
      &  theta_i_sn          , & ! volumetric ice content               [ - ]
      &  theta_w_sn          , & ! volumetric water content             [ - ]
      &  theta_a_sn          , & ! volumetric air content               [ - ]
      &  dzm_sn                  ! snow layer thickness                 [ m ]


  CONTAINS
    PROCEDURE :: Init => Init_nix_memory
  END TYPE t_nix_memory

  CHARACTER(len=*), PARAMETER :: modname = 'mo_nix_memory_class'

CONTAINS

  ! ======================================================================================================= !
  !> initialize memory (variables) for the process: nix
  !>
  SUBROUTINE Init_nix_memory(mem, prefix, suffix, lct_ids, lib_id, model_id)

    USE mo_jsb_varlist,       ONLY: BASIC, MEDIUM, FULL
    USE mo_jsb_io,            ONLY: grib_bits, t_cf, t_grib1, t_grib2, tables, TSTEP_CONSTANT
    USE mo_jsb_grid_class,    ONLY: t_jsb_grid, t_jsb_vgrid
    USE mo_jsb_grid,          ONLY: Get_grid, Get_vgrid
    USE mo_jsb_model_class,   ONLY: t_jsb_model, MODEL_QUINCY, MODEL_JSBACH
    USE mo_jsb_class,         ONLY: Get_model
!   USE mo_nix_constants

    CLASS(t_nix_memory), INTENT(inout), TARGET :: mem
    CHARACTER(len=*),    INTENT(in)            :: prefix
    CHARACTER(len=*),    INTENT(in)            :: suffix
    INTEGER,             INTENT(in)            :: lct_ids(:)
    INTEGER,             INTENT(in)            :: lib_id     !< id of primary lct in lctlib
    INTEGER,             INTENT(in)            :: model_id

    dsl4jsb_Def_config(NIX_)

    TYPE(t_jsb_model), POINTER  :: model                        !< model
    TYPE(t_jsb_grid),  POINTER  :: hgrid                        !< Horizontal grid
    TYPE(t_jsb_vgrid), POINTER  :: surface, snow_layer          !< Vertical grids
    INTEGER                     :: table                        !< ...

    CHARACTER(len=*), PARAMETER :: routine = modname//':Init_nix_memory'

    IF (model_id > 0) CONTINUE ! avoid compiler warning about dummy argument not being used

    model        => Get_model(model_id)
    table        =  tables(1)
    hgrid        => Get_grid(mem%grid_id)
    surface      => Get_vgrid('surface')
    snow_layer   => Get_vgrid('nix_vgrid')

      ! Prognostic variables  
      CALL mem%Add_var('t_sn', mem%t_sn,                                &
        & hgrid, snow_layer,                                            &
        & t_cf('weighted snow layer (element) temperature', '', ''),              &
        & t_grib1(table, 255, grib_bits), t_grib2(2, 3, 28, grib_bits), &
        & prefix, suffix,                                               &
        & lrestart=.FALSE., initval_r=0.0_wp )


      CALL mem%Add_var('t_sn_n', mem%t_sn_n,                                &
        & hgrid, snow_layer,                                            &
        & t_cf('weighted snow layer (nodal) temperature', '', ''),              &
        & t_grib1(table, 255, grib_bits), t_grib2(2, 3, 208, grib_bits), &
        & prefix, suffix,                                               &
        & lrestart=.FALSE., initval_r=0.0_wp )


      CALL mem%Add_var('theta_i_sn', mem%theta_i_sn,                                &
        & hgrid, snow_layer,                                            &
        & t_cf('weighted volumetric ice content', '', ''),              &
        & t_grib1(table, 255, grib_bits), t_grib2(2, 3, 209, grib_bits), &
        & prefix, suffix,                                               &
        & lrestart=.FALSE., initval_r=0.0_wp )


      CALL mem%Add_var('theta_w_sn', mem%theta_w_sn,                                &
        & hgrid, snow_layer,                                            &
        & t_cf('weighted volumetric water content', '', ''),              &
        & t_grib1(table, 255, grib_bits), t_grib2(2, 3, 210, grib_bits), &
        & prefix, suffix,                                               &
        & lrestart=.FALSE., initval_r=0.0_wp )


      CALL mem%Add_var('theta_a_sn', mem%theta_a_sn,                                &
        & hgrid, snow_layer,                                            &
        & t_cf('weighted volumetric air content', '', ''),              &
        & t_grib1(table, 255, grib_bits), t_grib2(2, 3, 211, grib_bits), &
        & prefix, suffix,                                               &
        & lrestart=.FALSE., initval_r=0.0_wp )


      CALL mem%Add_var('dzm_sn', mem%dzm_sn,                                &
        & hgrid, snow_layer,                                            &
        & t_cf('weighted snow layer thickness', '', ''),              &
        & t_grib1(table, 255, grib_bits), t_grib2(0, 1, 11, grib_bits), &
        & prefix, suffix,                                               &
        & lrestart=.FALSE., initval_r=0.0_wp )

















  END SUBROUTINE Init_nix_memory







#endif
END MODULE mo_nix_memory_class



