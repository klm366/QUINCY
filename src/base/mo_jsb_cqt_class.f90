!> Contains conserved quantity types, definitions and methods for to-be-conserved quantities
!>
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
!>#### Defines types dealing with conserved quantity types (CQTs) and methods to conserve quantities
!>
!> **Note:**
!>
!> - A CQT is not a fortran derived TYPE. It is rather a keyword specifying a certain way of
!>   handling conserved quantities (CT).
!> - All variables sharing a CQT need to have the same unit and number of dimensions (2D and 3D
!>   variables cannot be of the same CQT).
!> - All variables sharing a CQT will be treated in the same way by a certain lcc process
!>   (scaling, active or passive matter relocation).
!>
MODULE mo_jsb_cqt_class
#ifndef __NO_JSBACH__

  USE mo_util,                ONLY: int2string
  USE mo_exception,           ONLY: finish
  USE mo_jsb_var_class,       ONLY: t_jsb_var_p

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: WATER_CQ_TYPE
  PUBLIC :: FLUX_C_CQ_TYPE, LIVE_CARBON_CQ_TYPE, AG_DEAD_C_CQ_TYPE, BG_DEAD_C_CQ_TYPE, PRODUCT_CARBON_CQ_TYPE
  PUBLIC :: IQ_2L2D_POOL_CQ_TYPE, IQ_1L2D_POOL_CQ_TYPE, IQ_FLUX_CQ_TYPE, IQ_SL_POOL_CQ_TYPE, IQ_SL_FLUX_CQ_TYPE
  PUBLIC :: Get_number_of_types, Get_cqt_name, t_jsb_consQuan, t_jsb_consQuan_p, max_cqt_name_length

  ENUM, BIND(C)
    ENUMERATOR ::                 &
      & WATER_CQ_TYPE = 0,        & !< Currently still unused conserved quantity type for water pools
      & LIVE_CARBON_CQ_TYPE,      & !< Jsbach: CQT for living carbon pools
      & AG_DEAD_C_CQ_TYPE,        & !< Jsbach: CQT for above ground dead carbon pools
      & BG_DEAD_C_CQ_TYPE,        & !< Jsbach: CQT for below ground dead carbon pools
      & PRODUCT_CARBON_CQ_TYPE,   & !< Jsbach: CQT for carbon products
      & FLUX_C_CQ_TYPE,           & !< Jsbach: CQT for carbon fluxes
      & IQ_2L2D_POOL_CQ_TYPE,     & !< Quincy: CQT for two layer 2D pool type
      & IQ_1L2D_POOL_CQ_TYPE,     & !< Quincy: CQT for one layer 2D pool type
      & IQ_FLUX_CQ_TYPE,          & !< Quincy: CQT for fluxes
      & IQ_SL_POOL_CQ_TYPE,       & !< Quincy: CQT for soil layered pools
      & IQ_SL_FLUX_CQ_TYPE,       & !< Quincy: CQT for soil layered fluxes
      & LAST_CQ_TYPE !< needs to always be the last -- it is used to determine the max number of CQTs
  END ENUM

  !> Type used on tiles to collect conserved quantities from the process memories of the tile
  !
  TYPE :: t_jsb_consQuan
    INTEGER :: type_id = -1           !< ID of the CQT (see above)
    INTEGER :: var_type = 0           !< Variable type (REAL2D or REAL3D) of the CQ variables of this CQT.
    LOGICAL :: transfer_all = .FALSE. !< Flag indicating whether matter of the whole tile is collected
                                      !< (true) or if only a fraction proportional to the lost area fraction will
                                      !< be relocated (false; default); same for all variables of the CQT
    INTEGER :: no_of_vars = 0         !< Number of CQ variables in this CQT
    INTEGER :: last_index_used = 0    !< Last index used - required for filling of cq_vars
    INTEGER, ALLOCATABLE :: associated_process(:) !< Processes to which the variables of this CQT belong
    TYPE(t_jsb_var_p), ALLOCATABLE :: cq_vars(:)  !< Collection of all variables of the CQT
  END TYPE t_jsb_consQuan

  !> Type to hold vectors of pointers to t_jsb_consQuan objects
  TYPE :: t_jsb_consQuan_p
    TYPE(t_jsb_consQuan), POINTER :: p => NULL()
  END TYPE t_jsb_consQuan_p

  CHARACTER(len=*), PARAMETER :: modname = 'mo_jsb_cqt_class'

  INTEGER, PARAMETER :: max_cqt_name_length = 25 !< Maximum number of characters allowed for the name of a CQT

CONTAINS

  ! ====================================================================================================== !
  !
  !> Returns the number of CQTs defined in the enumerator
  !
  FUNCTION Get_number_of_types() RESULT(last_type_id)
    ! -------------------------------------------------------------------------------------------------- !
    INTEGER :: last_type_id  !< Number of CQTs defined in the enumerator
    ! -------------------------------------------------------------------------------------------------- !
    last_type_id = LAST_CQ_TYPE

  END FUNCTION Get_number_of_types

  ! ====================================================================================================== !
  !
  !> Returns the name of the CQT with the given id - restricted length (max_cqt_name_length)
  !
  FUNCTION Get_cqt_name(id) RESULT(return_value)
    ! -------------------------------------------------------------------------------------------------- !
    INTEGER,  INTENT(in)  :: id                     !< Id of the CQT whose name is requested
    CHARACTER(len=:), ALLOCATABLE :: return_value   !< Name of the CQT with the given id
    ! -------------------------------------------------------------------------------------------------- !
    CHARACTER(len=*), PARAMETER :: routine = modname//':Get_cqt_name'
    ! -------------------------------------------------------------------------------------------------- !

    SELECT CASE(id)
      CASE (WATER_CQ_TYPE)
        return_value = 'water'
      CASE (LIVE_CARBON_CQ_TYPE)
        return_value = 'living_carbon'
      CASE (AG_DEAD_C_CQ_TYPE)
        return_value = 'ag_dead_carbon'
      CASE (BG_DEAD_C_CQ_TYPE)
        return_value = 'bg_dead_carbon'
      CASE (PRODUCT_CARBON_CQ_TYPE)
        return_value = 'product_carbon'
      CASE (FLUX_C_CQ_TYPE)
        return_value = 'carbon_flux'
      CASE (IQ_2L2D_POOL_CQ_TYPE)
        return_value = 'quincy_bgcm_2L2D_pool'
      CASE (IQ_1L2D_POOL_CQ_TYPE)
        return_value = 'quincy_bgcm_1L2D_pool'
      CASE (IQ_FLUX_CQ_TYPE)
        return_value = 'quincy_bgcm_flux'
      CASE (IQ_SL_POOL_CQ_TYPE)
        return_value = 'quincy_sl_bgcm_pool'
      CASE (IQ_SL_FLUX_CQ_TYPE)
        return_value = 'quincy_sl_bgcm_flux'
      CASE DEFAULT
        CALL finish(TRIM(routine), 'No name specified for CQT of id '//int2string(id))
    END SELECT

  END FUNCTION Get_cqt_name

#endif
END MODULE mo_jsb_cqt_class
