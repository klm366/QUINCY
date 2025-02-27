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
!> For more information on the QUINCY model see: <https://doi.org/10.17871/quincy-model-2019>
!>
!> @author
!>  Albert Jornet Puig, University of Leipzig / MPI-BGC Jena
!>  Ana Bastos, University of Leipzig / MPI-BGC Jena
!>
MODULE mo_dfire_sf_adaptor
#ifndef __NO_JSBACH__

  USE mo_kind,                ONLY: wp
  USE mo_jsb_control,         ONLY: debug_on
  USE mo_exception,           ONLY: finish, message_text
  USE mo_dfire_constants,     ONLY: df_sf_carbon_to_biomass_factor, df_sf_units_g_to_kg

  USE mo_dfire_constants,     ONLY: df_sk_units_kgm2s_to_mmday, df_sf_units_k_to_c, &
                                    df_sf_units_sec_to_min, df_sf_units_m_to_cm, &
                                    df_sf_units_m2_to_ha
  USE mo_jsb_physical_constants, ONLY: molar_mass_C

  IMPLICIT NONE

  PUBLIC :: sf_convert_units!, sf_convert_units_mapped_wp_r2, convert_units_real_wp_r3

  PRIVATE

  CHARACTER(len=*), PARAMETER :: modname = 'mo_dfire_sf_adaptor'

  INTERFACE sf_convert_units
    MODULE PROCEDURE :: convert_units_real_wp_r1, convert_units_real_wp_r2, &
                        convert_units_real_wp_scalar, convert_units_real_wp_r3, &
                        convert_units_real_wp_r4, sf_convert_units_mapped_wp_r2
  END INTERFACE sf_convert_units

CONTAINS

  ! ======================================================================================================= !
  !>
  !> Units conversion from 'src' to 'dst' from the given variable 'var' for scalar wp values
  !>
  RECURSIVE FUNCTION convert_units_real_wp_scalar(src, dst, varIn, aux) RESULT (varOut)
    USE mo_dfire_constants, ONLY: df_sf_units_km2_to_ha
    CHARACTER(*), INTENT(in)        :: src, dst
    REAL(wp), INTENT(in)            :: varIn
    REAL(wp), INTENT(in), OPTIONAL  :: aux        !< auxiliary values for conversion
    REAL(wp)                        :: varOut     !< output

    ! ----------------------------------------------------------------------------------------------------- !
    IF (TRIM(src) == "km2" .AND. TRIM(dst) == "ha") THEN
      varOut = varIn * df_sf_units_km2_to_ha
    ELSE IF (TRIM(src) == "g" .AND. TRIM(dst) == "kg") THEN
      varOut = varIn * 0.001
    ELSE
      message_text = "Units conversion from "// src //" to "// dst &
          & // " not yet implemented"
      call finish("convert_units_real_wp_scalar", message_text)
    END IF
  END FUNCTION convert_units_real_wp_scalar


  ! ======================================================================================================= !
  !>
  !> Units conversion from 'src' to 'dst' from the given variable 'var'
  !>
  RECURSIVE FUNCTION convert_units_real_wp_r1(src, dst, varIn, aux) RESULT (varOut)
    CHARACTER(*), INTENT(in)        :: src, dst
    REAL(wp), INTENT(in)            :: varIn(:)
    REAL(wp), INTENT(in), OPTIONAL  :: aux(:)    !< auxiliary values for conversion
    REAL(wp), ALLOCATABLE           :: var_m2(:), var_g_m2(:) !< tmp variables
    REAL(wp), ALLOCATABLE           :: varOut(:) !< output
    ! ----------------------------------------------------------------------------------------------------- !
    IF (.NOT. ALLOCATED(varOut)) ALLOCATE(varOut, mold = varIn)

    IF (TRIM(src) == "Kelvin" .AND. TRIM(dst) == "Celsius") THEN
      varOut = varIn + df_sf_units_k_to_c
    ELSE IF (TRIM(src) == "m2" .AND. TRIM(dst) == "ha") THEN
      varOut = varIn * df_sf_units_m2_to_ha
    ELSE IF (TRIM(src) == "g" .AND. TRIM(dst) == "kg") THEN
      varOut = varIn * df_sf_units_g_to_kg
    ELSE IF (TRIM(src) == "kg m-2 s-1" .AND. TRIM(dst) == "mm/day") THEN
      varOut = varIn * df_sk_units_kgm2s_to_mmday
    ELSE IF (TRIM(src) == "m/s" .AND. TRIM(dst) == "m/min") THEN
      varOut = varIn * df_sf_units_sec_to_min
    ELSE IF (TRIM(src) == "m" .AND. TRIM(dst) == "cm") THEN
      varOut = varIn * df_sf_units_m_to_cm
    ELSE IF (TRIM(src) == "-/m3" .AND. TRIM(dst) == "-/m2") THEN
      IF (.NOT. PRESENT(aux))  THEN
        CALL finish("convert_units_real_wp_r1", &
                "When converting from x/m3 to x/m2 it is required to provide " // &
                        "the conversion factor for each layer in [m]")
      END IF

      varOut = varIn * aux
    ELSE IF (TRIM(src) == "-/m2" .AND. TRIM(dst) == "-/m3") THEN
      IF (.NOT. PRESENT(aux))  THEN
        CALL finish("convert_units_real_wp_r1", &
                "When converting from x/m2 to x/m3 it is required to provide " // &
                        "the conversion factor for each layer in [m]")
      END IF

      varOut = varIn / aux
    ELSE IF (TRIM(src) == "molC" .AND. TRIM(dst) == "gC") THEN
      varOut = varIn * molar_mass_C
    ELSE IF (TRIM(src) == "gC" .AND. TRIM(dst) == "molC") THEN
      varOut = varIn / molar_mass_C
    ELSE IF (TRIM(src) == "C" .AND. TRIM(dst) == "DM") THEN
        varOut = varIn * df_sf_carbon_to_biomass_factor
    ELSE IF (TRIM(src) == "DM" .AND. TRIM(dst) == "C") THEN
        varOut = varIn / df_sf_carbon_to_biomass_factor
    ELSE IF (TRIM(src) == "-/m3" .AND. TRIM(dst) == "-/m2") THEN
      IF (.NOT. PRESENT(aux))  THEN
        CALL finish("convert_units_real_wp_r1", &
                "When converting from x/m3 to x/m2 it is required to provide " // &
                        "the conversion factor for each layer in [m]")
      END IF
      varOut = varIn * aux
    ELSE IF (TRIM(src) == "molC/m3" .AND. TRIM(dst) == "gDM/m2") THEN !< higher level functions
        ALLOCATE(var_m2, var_g_m2, mold = varIn)
        var_m2 = sf_convert_units("-/m3", "-/m2", varIn, aux)
        var_g_m2 = sf_convert_units("molC", "gC", var_m2)
        varOut = sf_convert_units("C", "DM", var_g_m2)
    ELSE IF (TRIM(src) == "gDM/m2" .AND. TRIM(dst) == "molC/m3") THEN !< higher level functions
        ALLOCATE(var_m2, var_g_m2, mold = varIn)
        var_m2 = sf_convert_units("DM", "C", varIn, aux)
        var_g_m2 = sf_convert_units("gC", "molC", var_m2)
        varOut = sf_convert_units("-/m2", "-/m3", var_g_m2, aux)
    ELSE IF (TRIM(src) == "molC/m2" .AND. TRIM(dst) == "gDM/m2") THEN !< higher level functions
        ALLOCATE(var_g_m2, mold = varIn)
        var_g_m2 = sf_convert_units("molC", "gC", varIn)
        varOut = sf_convert_units("C", "DM", var_g_m2)
    ELSE IF (TRIM(src) == "gDM/m2" .AND. TRIM(dst) == "molC/m2") THEN !< higher level functions
        ALLOCATE(var_g_m2, mold = varIn)
        var_g_m2 = sf_convert_units("DM", "C", varIn)
        varOut = sf_convert_units("gC", "molC", var_g_m2)
    ELSE
        message_text = "Units conversion from "// src //" to "// dst &
        & // " not yet implemented"
        call finish("convert_units_real_wp_r1", message_text)
    END IF
  END FUNCTION convert_units_real_wp_r1


  RECURSIVE FUNCTION convert_units_real_wp_r2(src, dst, varIn, aux, nc_dim_opt) RESULT (varOut)
    CHARACTER(*), INTENT(in)        :: src, dst
    REAL(wp), INTENT(in)            :: varIn(:,:)   !< (#nc, #ndim)
    REAL(wp), INTENT(in), OPTIONAL  :: aux(:)       !< auxiliary values for conversion
    REAL(wp), ALLOCATABLE           :: varOut(:,:)  !< output
    INTEGER, OPTIONAL, INTENT(in)   :: nc_dim_opt
    ! ----------------------------------------------------------------------------------------------------- !
    REAL(wp), ALLOCATABLE           :: var_m2(:,:), var_g_m2(:,:) !< tmp variables
    INTEGER :: inc, nc_dim
    ! ----------------------------------------------------------------------------------------------------- !
    nc_dim = 1 !< default, first dimension of array
    IF (PRESENT(nc_dim_opt)) THEN
      nc_dim = nc_dim_opt
    END IF

    IF (.NOT. ALLOCATED(varOut)) ALLOCATE(varOut, mold = varIn)

    IF (TRIM(src) == "Kelvin" .AND. TRIM(dst) == "Celsius") THEN
      varOut = varIn + df_sf_units_k_to_c
    ELSE IF (TRIM(src) == "g" .AND. TRIM(dst) == "kg") THEN
      varOut = varIn * df_sf_units_g_to_kg
    ELSE IF (TRIM(src) == "kg m-2 s-1" .AND. TRIM(dst) == "mm/day") THEN
      varOut = varIn * df_sk_units_kgm2s_to_mmday
    ELSE IF (TRIM(src) == "molC" .AND. TRIM(dst) == "gC") THEN
      varOut = varIn * molar_mass_C
    ELSE IF (TRIM(src) == "gC" .AND. TRIM(dst) == "molC") THEN
      varOut = varIn / molar_mass_C
    ELSE IF (TRIM(src) == "C" .AND. TRIM(dst) == "DM") THEN
      varOut = varIn * df_sf_carbon_to_biomass_factor
    ELSE IF (TRIM(src) == "DM" .AND. TRIM(dst) == "C") THEN
      varOut = varIn / df_sf_carbon_to_biomass_factor
    ELSE IF (TRIM(src) == "-/m3" .AND. TRIM(dst) == "-/m2") THEN
      IF (.NOT. PRESENT(aux))  THEN
        CALL finish("convert_units_real_wp_r2", &
                "When converting from x/m3 to x/m2 it is required to provide " // &
                        "the conversion factor for each layer in [m]")
      END IF
      varOut = multiply_by_ncdim_wp_r2(varIn, aux, nc_dim)

    ELSE IF (TRIM(src) == "-/m2" .AND. TRIM(dst) == "-/m3") THEN
      IF (.NOT. PRESENT(aux))  THEN
        CALL finish("convert_units_real_wp_r2", &
                "When converting from x/m2 to x/m3 it is required to provide " // &
                        "the conversion factor for each layer in [m]")
      END IF
      varOut = division_by_ncdim_wp_r2(varIn, aux, nc_dim)

    ELSE IF (TRIM(src) == "molC/m3" .AND. TRIM(dst) == "gC/m2") THEN !< higher level functions
      ALLOCATE(var_m2, mold = varIn)
      var_m2 = sf_convert_units("-/m3", "-/m2", varIn, aux)
      varOut = sf_convert_units("molC", "gC", var_m2)
    ELSE IF (TRIM(src) == "molC/m3" .AND. TRIM(dst) == "gDM/m2") THEN !< higher level functions
        ALLOCATE(var_m2, mold = varIn)
        ALLOCATE(var_g_m2, mold = varIn)
        var_m2 = sf_convert_units("-/m3", "-/m2", varIn, aux, nc_dim)
        var_g_m2 = sf_convert_units("molC", "gC", var_m2)
        varOut = sf_convert_units("C", "DM", var_g_m2)
    ELSE IF (TRIM(src) == "gDM/m2" .AND. TRIM(dst) == "molC/m3") THEN !< higher level functions
        ALLOCATE(var_m2, mold = varIn)
        ALLOCATE(var_g_m2, mold = varIn)
        var_m2 = sf_convert_units("DM", "C", varIn, aux)
        var_g_m2 = sf_convert_units("gC", "molC", var_m2)
        varOut = sf_convert_units("-/m2", "-/m3", var_g_m2, aux)
    ELSE IF (TRIM(src) == "molC/m2" .AND. TRIM(dst) == "gDM/m2") THEN !< higher level functions
        ALLOCATE(var_g_m2, mold = varIn)
        var_g_m2 = sf_convert_units("molC", "gC", varIn)
        varOut = sf_convert_units("C", "DM", var_g_m2)
    ELSE IF (TRIM(src) == "gDM/m2" .AND. TRIM(dst) == "molC/m2") THEN !< higher level functions
        ALLOCATE(var_g_m2, mold = varIN)
        var_g_m2 = sf_convert_units("DM", "C", varIn)
        varOut = sf_convert_units("gC", "molC", var_g_m2)
    ELSE
      message_text = "Units conversion from "// src //" to "// dst &
              & // " not yet implemented"
      call finish("convert_units_real_wp_r2", message_text)
    END IF

  CONTAINS
    !< Account for nc dimension (gridcell) as aux must match this dimensions
    FUNCTION multiply_by_ncdim_wp_r2(varIn, aux, nc_dim) RESULT (varOut)
      REAL(wp), INTENT(in) :: varIn(:,:)
      REAL(wp), INTENT(in) :: aux(:)            !< (#nc)
      INTEGER, INTENT(in)  :: nc_dim
      REAL(wp), ALLOCATABLE :: varOut(:,:)

      INTEGER :: inc

      IF (.NOT. ALLOCATED(varOut)) ALLOCATE(varOut, mold=varIn)

      SELECT CASE (nc_dim)
      CASE (1)
        DO inc = 1, size(varIn, dim=NC_DIM)
          varOut(inc,:) = varIn(inc,:) * aux(inc)
        END DO
      CASE(2)
        DO inc = 1, size(varIn, dim=NC_DIM)
          varOut(:,inc) = varIn(:,inc) * aux(inc)
        END DO
      CASE DEFAULT
        CALL finish("multiply_by_ncdim_wp_r2", "Undefined nc_dim_opt dimension, " // &
            "only allowed  1(default),2,3")
      END SELECT
    END FUNCTION multiply_by_ncdim_wp_r2

    !< Same as multiply but for division
    FUNCTION division_by_ncdim_wp_r2(varIn, aux, nc_dim) RESULT (varOut)
      REAL(wp), INTENT(in) :: varIn(:,:)
      REAL(wp), INTENT(in) :: aux(:)
      INTEGER, INTENT(in)  :: nc_dim
      REAL(wp), ALLOCATABLE :: varOut(:,:)

      INTEGER :: inc

      IF (.NOT. ALLOCATED(varOut)) ALLOCATE(varOut, mold=varIn)

      SELECT CASE (nc_dim)
      CASE (1)
        DO inc = 1, size(varIn, dim=NC_DIM)
          varOut(inc,:) = varIn(inc,:) / aux(inc)
        END DO
      CASE(2)
        DO inc = 1, size(varIn, dim=NC_DIM)
          varOut(:,inc) = varIn(:,inc) / aux(inc)
        END DO
      CASE DEFAULT
        CALL finish("division_by_ncdim_wp_r2", "Undefined nc_dim_opt dimension, " // &
            "only allowed  1(default),2,3")
      END SELECT
    END FUNCTION division_by_ncdim_wp_r2

  END FUNCTION convert_units_real_wp_r2


  RECURSIVE FUNCTION convert_units_real_wp_r3(src, dst, varIn, aux, nc_dim_opt) RESULT (varOut)
    CHARACTER(*), INTENT(in)        :: src, dst
    REAL(wp), INTENT(in)            :: varIn(:,:,:)
    REAL(wp), INTENT(in), OPTIONAL  :: aux(:)          !< auxiliary values for conversion
    INTEGER, INTENT(in), OPTIONAL   :: nc_dim_opt      !< Define which dimension is nc. It works along with 'aux' as it has to match the nc dimension
    REAL(wp), ALLOCATABLE           :: varOut(:,:,:)   !< output
    ! ----------------------------------------------------------------------------------------------------- !
    REAL(wp), ALLOCATABLE           :: var_m2(:,:,:), var_g_m2(:,:,:) !< tmp variables
    INTEGER, PARAMETER              :: nc_dim_def = 1
    INTEGER                         :: nc_dim
    ! ----------------------------------------------------------------------------------------------------- !
    IF (.NOT. ALLOCATED(varOut)) ALLOCATE(varOut, mold=varIn)
    nc_dim = nc_dim_def
    IF (PRESENT(nc_dim_opt)) THEN
      nc_dim = nc_dim_opt
    ENDIF

    IF (TRIM(src) == "Kelvin" .AND. TRIM(dst) == "Celsius") THEN
      varOut = varIn + df_sf_units_k_to_c
    ELSE IF (TRIM(src) == "g" .AND. TRIM(dst) == "kg") THEN
      varOut = varIn * df_sf_units_g_to_kg
    ELSE IF (TRIM(src) == "kg m-2 s-1" .AND. TRIM(dst) == "mm/day") THEN
      varOut = varIn * df_sk_units_kgm2s_to_mmday
    ELSE IF (TRIM(src) == "molC" .AND. TRIM(dst) == "gC") THEN
      varOut = varIn * molar_mass_C
    ELSE IF (TRIM(src) == "gC" .AND. TRIM(dst) == "molC") THEN
      varOut = varIn / molar_mass_C
    ELSE IF (TRIM(src) == "C" .AND. TRIM(dst) == "DM") THEN
      varOut = varIn * df_sf_carbon_to_biomass_factor
    ELSE IF (TRIM(src) == "DM" .AND. TRIM(dst) == "C") THEN
      varOut = varIn / df_sf_carbon_to_biomass_factor
    ELSE IF (TRIM(src) == "-/m3" .AND. TRIM(dst) == "-/m2") THEN
      IF (.NOT. PRESENT(aux))  THEN
        CALL finish("convert_units_real_wp_r3", &
                "When converting from x/m3 to x/m2 it is required to provide " // &
                        "the conversion factor for each layer in [m]")
      END IF
      varOut = multiply_by_ncdim_wp_r3(varIn, aux, nc_dim)

    ELSE IF (TRIM(src) == "-/m2" .AND. TRIM(dst) == "-/m3") THEN
      IF (.NOT. PRESENT(aux))  THEN
        CALL finish("convert_units_real_wp_r3", &
                "When converting from x/m2 to x/m3 it is required to provide " // &
                        "the conversion factor for each layer in [m]")
      END IF
      varOut = division_by_ncdim_wp_r3(varIn, aux, nc_dim)

    ELSE IF (TRIM(src) == "molC/m3" .AND. TRIM(dst) == "gC/m2") THEN !< higher level functions
      ALLOCATE(var_m2, mold = varIn)
      var_m2 = sf_convert_units("-/m3", "-/m2", varIn, aux, nc_dim)
      varOut = sf_convert_units("molC", "gC", var_m2)
    ELSE IF (TRIM(src) == "molC/m3" .AND. TRIM(dst) == "gDM/m2") THEN !< higher level functions
        ALLOCATE(var_m2, mold = varIn)
        ALLOCATE(var_g_m2, mold = varIn)
        var_m2 = sf_convert_units("-/m3", "-/m2", varIn, aux, nc_dim)
        var_g_m2 = sf_convert_units("molC", "gC", var_m2)
        varOut = sf_convert_units("C", "DM", var_g_m2)
    ELSE IF (TRIM(src) == "gDM/m2" .AND. TRIM(dst) == "molC/m3") THEN !< higher level functions
        ALLOCATE(var_m2, mold = varIn)
        ALLOCATE(var_g_m2, mold = varIn)
        var_m2 = sf_convert_units("DM", "C", varIn, aux, nc_dim)
        var_g_m2 = sf_convert_units("gC", "molC", var_m2)
        varOut = sf_convert_units("-/m2", "-/m3", var_g_m2, aux, nc_dim)
    ELSE IF (TRIM(src) == "molC/m2" .AND. TRIM(dst) == "gDM/m2") THEN !< higher level functions
        ALLOCATE(var_g_m2, mold = varIn)
        var_g_m2 = sf_convert_units("molC", "gC", varIn)
        varOut = sf_convert_units("C", "DM", var_g_m2)
    ELSE IF (TRIM(src) == "gDM/m2" .AND. TRIM(dst) == "molC/m2") THEN !< higher level functions
        ALLOCATE(var_g_m2, mold = varIN)
        var_g_m2 = sf_convert_units("DM", "C", varIn)
        varOut = sf_convert_units("gC", "molC", var_g_m2)
    ELSE
      message_text = "Units conversion from "// src //" to "// dst &
              & // " not yet implemented"
      call finish("convert_units_real_wp_r3", message_text)
    END IF

  CONTAINS
    !< Account for nc dimension (gridcell) as aux must match this dimensions
    FUNCTION multiply_by_ncdim_wp_r3(varIn, aux, nc_dim) RESULT (varOut)
      REAL(wp), INTENT(in) :: varIn(:,:,:)
      REAL(wp), INTENT(in) :: aux(:)            !< (#nc)
      INTEGER, INTENT(in)  :: nc_dim
      REAL(wp), ALLOCATABLE :: varOut(:,:,:)

      INTEGER :: inc

      IF (.NOT. ALLOCATED(varOut)) ALLOCATE(varOut, mold=varIn)

      SELECT CASE (nc_dim)
      CASE (1)
        DO inc = 1, size(varIn, dim=NC_DIM)
          varOut(inc,:,:) = varIn(inc,:,:) * aux(inc)
        END DO
      CASE(2)
        DO inc = 1, size(varIn, dim=NC_DIM)
          varOut(:,inc,:) = varIn(:,inc,:) * aux(inc)
        END DO
      CASE(3)
        DO inc = 1, size(varIn, dim=NC_DIM)
          varOut(:,:,inc) = varIn(:,:,inc) * aux(inc)
        END DO
      CASE DEFAULT
        CALL finish("multiply_by_ncdim_wp_r3", "Undefined nc_dim_opt dimension, " // &
                        "only allowed  1(default),2,3")
      END SELECT
    END FUNCTION multiply_by_ncdim_wp_r3

    !< Same as multiply but for division
    FUNCTION division_by_ncdim_wp_r3(varIn, aux, nc_dim) RESULT (varOut)
      REAL(wp), INTENT(in) :: varIn(:,:,:)
      REAL(wp), INTENT(in) :: aux(:)
      INTEGER, INTENT(in)  :: nc_dim
      REAL(wp), ALLOCATABLE :: varOut(:,:,:)

      INTEGER :: inc

      IF (.NOT. ALLOCATED(varOut)) ALLOCATE(varOut, mold=varIn)

      SELECT CASE (nc_dim)
      CASE (1)
        DO inc = 1, size(varIn, dim=NC_DIM)
          varOut(inc,:,:) = varIn(inc,:,:) / aux(inc)
        END DO
      CASE(2)
        DO inc = 1, size(varIn, dim=NC_DIM)
          varOut(:,inc,:) = varIn(:,inc,:) / aux(inc)
        END DO
      CASE(3)
        DO inc = 1, size(varIn, dim=NC_DIM)
          varOut(:,:,inc) = varIn(:,:,inc) / aux(inc)
        END DO
      CASE DEFAULT
        CALL finish("division_by_ncdim_wp_r3", "Undefined nc_dim_opt dimension, " // &
                        "only allowed  1(default),2,3")
      END SELECT
    END FUNCTION division_by_ncdim_wp_r3

  END FUNCTION convert_units_real_wp_r3



  FUNCTION convert_units_real_wp_r4(src, dst, varIn, aux, nc_dim_opt) RESULT (varOut)
    CHARACTER(*), INTENT(in)        :: src, dst
    REAL(wp), INTENT(in)            :: varIn(:,:,:,:)
    REAL(wp), INTENT(in), OPTIONAL  :: aux(:)          !< auxiliary values for conversion
    INTEGER, INTENT(in), OPTIONAL   :: nc_dim_opt      !< Define which dimension is nc. It works along with 'aux' as it has to match the nc dimension
    REAL(wp), ALLOCATABLE           :: varOut(:,:,:,:)   !< output
    ! ----------------------------------------------------------------------------------------------------- !
    INTEGER, PARAMETER              :: nc_dim_def = 1
    INTEGER                         :: nc_dim
    ! ----------------------------------------------------------------------------------------------------- !
    IF (.NOT. ALLOCATED(varOut)) ALLOCATE(varOut, mold=varIn)

    nc_dim = nc_dim_def
    IF (PRESENT(nc_dim_opt)) THEN
      nc_dim = nc_dim_opt
    ENDIF

    IF (TRIM(src) == "molC/m3" .AND. TRIM(dst) == "gDM/m2") THEN !< higher level functions
        !< molC / m3 -> molC/m2 -> gC / m2 -> gDM / m2
        varOut = multiply_by_ncdim_wp_r4(varIn, aux, nc_dim) * molar_mass_C * df_sf_carbon_to_biomass_factor
    ELSE
      message_text = "Units conversion from "// src //" to "// dst &
              & // " not yet implemented"
      call finish("convert_units_real_wp_r4", message_text)
    END IF

  CONTAINS
    !< Account for nc dimension (gridcell) as aux must match this dimensions
    FUNCTION multiply_by_ncdim_wp_r4(varIn, aux, nc_dim) RESULT (varOut)
      REAL(wp), INTENT(in) :: varIn(:,:,:,:)
      REAL(wp), INTENT(in) :: aux(:)            !< (#nc)
      INTEGER, INTENT(in)  :: nc_dim
      REAL(wp), ALLOCATABLE :: varOut(:,:,:,:)

      INTEGER :: inc

      IF (.NOT. ALLOCATED(varOut)) ALLOCATE(varOut, mold=varIn)

      SELECT CASE (nc_dim)
      CASE (1)
        DO inc = 1, size(varIn, dim=NC_DIM)
          varOut(inc,:,:,:) = varIn(inc,:,:,:) * aux(inc)
        END DO
      CASE(2)
        DO inc = 1, size(varIn, dim=NC_DIM)
          varOut(:,inc,:,:) = varIn(:,inc,:,:) * aux(inc)
        END DO
      CASE(3)
        DO inc = 1, size(varIn, dim=NC_DIM)
          varOut(:,:,inc,:) = varIn(:,:,inc,:) * aux(inc)
        END DO
      CASE(4)
        DO inc = 1, size(varIn, dim=NC_DIM)
          varOut(:,:,:,inc) = varIn(:,:,:,inc) * aux(inc)
        END DO
      CASE DEFAULT
        CALL finish("multiply_by_ncdim_wp_r4", "Undefined nc_dim_opt dimension, " // &
                        "only allowed  1(default),2,3,4")
      END SELECT
    END FUNCTION multiply_by_ncdim_wp_r4

  END FUNCTION convert_units_real_wp_r4


  ! ======================================================================================================= !
  !>
  !> Convert to different units according to dsts.
  !>
  FUNCTION sf_convert_units_mapped_wp_r2 (src, dsts, varIn, dim, aux) result (varOut)
    CHARACTER(len=*), INTENT(in) :: src       !<
    CHARACTER(len=*), INTENT(in) :: dsts(:)   !< dsts size == as SIZE(varIn, dim="dim")
    REAL(wp), INTENT(in) :: varIn(:,:)        !< (#dim1, #dim2), It assumes nc=dim1.
    INTEGER, INTENT(in) :: dim                !< Dimension from varIn must have the same size as dsts
    REAL(wp), INTENT(in), OPTIONAL :: aux(:)  !< Auxiliary data for units conversion (e.g: soil_depth_sl)
    REAL(wp), ALLOCATABLE :: varOut(:,:)      !< (#dim1, #dim2), It assumes nc=dim1

    INTEGER :: iunit

    IF (.NOT. ALLOCATED(varOut)) THEN
      ALLOCATE(varOut, mold=varIn)
      varOut = 0._wp
    END IF

    IF (SIZE(dsts) /= SIZE(varIn, dim=2)) THEN
      WRITE(*,*) "sf_convert_units_mapped_wp_r2:: size dsts=", SIZE(dsts)
      WRITE(*,*) "sf_convert_units_mapped_wp_r2:: size varIn(.,X)=", SIZE(varIn,dim=2)
      CALL finish("sf_convert_units_mapped_wp_r2", "Differents dimensions on 'dsts' and 'varIn(:,X)'")
    ENDIF

    DO iunit=1, SIZE(varOut, dim=dim)
      varOut(:,iunit) = convert_units_real_wp_r1(src, dsts(iunit), varIn(:,iunit), aux=aux)
    END DO
  END FUNCTION sf_convert_units_mapped_wp_r2

#endif
END MODULE mo_dfire_sf_adaptor
