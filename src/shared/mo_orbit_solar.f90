!> Summary: Functions for orbital and solar parameters from host model
!>
!>### Provides functions from the host file adapters for orbital and solar computations
!>
!>#### Revision History
!> First version by Reiner Schnur (2022-08-26)
!>
!>#### Copyright and License
!> Copyright 2022 Max Planck Institute for Meteorology
!>
!> Redistribution and use in source and binary forms, with or
!> without modification, are permitted provided that the following
!> conditions are met:
!>
!> 1. Redistributions of source code must retain the above copyright
!>    notice, this list of conditions and the following disclaimer.
!>
!> 2. Redistributions in binary form must reproduce the above
!>    copyright notice, this list of conditions and the following
!>    disclaimer in the documentation and/or other materials
!>    provided with the distribution.
!>
!> 3. Neither the name of the copyright holder nor the names of its
!>    contributors may be used to endorse or promote products
!>    derived from this software without specific prior written
!>    permission.
!>
!> THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
!> "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
!> LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
!> A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
!> HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
!> SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
!> LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
!> DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
!> THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
!> (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
!> OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
!>

MODULE mo_orbit_solar
#ifndef __NO_JSBACH__

  USE mo_jsb_orbit_solar_iface, ONLY: inquire_declination, compute_cos_zenith_angle

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: inquire_declination, compute_cos_zenith_angle

  CHARACTER(len=*), PARAMETER :: modname = 'mo_orbit_solar'

! CONTAINS

#endif
END MODULE mo_orbit_solar
