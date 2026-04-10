!> QUINCY model output parameters
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
!> For more information on the QUINCY model see: <https://doi.org/10.17871/quincy-model-2019>
!>
!>#### defines QUINCY model parameters, units and helpers for model output
!>
MODULE mo_quincy_output_class
#ifndef __NO_QUINCY__

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: unit_flux, unit_pool, unit_soillayer_flux, unit_soillayer_pool, unit_cn, unit_nc, unit_np, unit_pn, unit_cp
  PUBLIC :: unitless, unit_hydro_pool, unit_hydro_flux, unit_water_potential, unit_veg_pool_flux, unit_assimi_rates
  PUBLIC :: unit_fraction, unit_energy, unit_temperature, unit_speed, unit_size, unit_bulk_dens, unit_rain
  PUBLIC :: unit_conductivity, unit_sb_pool, unit_sb_flux

  CHARACTER(len=*), PARAMETER, PRIVATE :: modname = 'mo_quincy_output_class'

  ! units
  CHARACTER(len=*), PARAMETER :: unit_flux            = 'micro mol m-2 s-1'
  CHARACTER(len=*), PARAMETER :: unit_pool            = 'mol m-2'
  CHARACTER(len=*), PARAMETER :: unit_soillayer_flux  = 'micro mol m-3 s-1'
  CHARACTER(len=*), PARAMETER :: unit_soillayer_pool  = 'mol m-3'
  !@TODO: unit? mol/mol?
  CHARACTER(len=*), PARAMETER :: unit_cn              = 'mol C mol-1 N'
  CHARACTER(len=*), PARAMETER :: unit_cp              = 'mol C mol-1 P'
  CHARACTER(len=*), PARAMETER :: unit_nc              = 'mol N mol-1 C'
  CHARACTER(len=*), PARAMETER :: unit_np              = 'mol N mol-1 P'
  CHARACTER(len=*), PARAMETER :: unit_pn              = 'mol P mol-1 N'

  CHARACTER(len=*), PARAMETER :: unit_assimi_rates    = 'micro mol CO2 m-2 s-1'

  CHARACTER(len=*), PARAMETER :: unit_bulk_dens       = 'kg m-3'

  CHARACTER(len=*), PARAMETER :: unit_rain            = 'mm s-1'
  CHARACTER(len=*), PARAMETER :: unit_hydro_flux      = 'kg m-2 s-1'
  CHARACTER(len=*), PARAMETER :: unit_hydro_pool      = 'm'
  CHARACTER(len=*), PARAMETER :: unit_water_potential = 'MPa'

  CHARACTER(len=*), PARAMETER :: unit_size            = "m"
  CHARACTER(len=*), PARAMETER :: unit_speed           = "m s-1"
  CHARACTER(len=*), PARAMETER :: unit_temperature     = "K"
  CHARACTER(len=*), PARAMETER :: unit_fraction        = "ppm"
  CHARACTER(len=*), PARAMETER :: unit_energy          = "W m-2"
  CHARACTER(len=*), PARAMETER :: unitless             = '--'
  CHARACTER(len=*), PARAMETER :: unit_conductivity    = "m s-1"

  CHARACTER(len=*), PARAMETER :: unit_veg_pool_flux   = 'mol m-2 timestep-1'

  CHARACTER(len=*), PARAMETER :: unit_sb_pool   = 'mol m-3'
  CHARACTER(len=*), PARAMETER :: unit_sb_flux   = 'mol m-3 timestep-1'

#endif
END MODULE mo_quincy_output_class
