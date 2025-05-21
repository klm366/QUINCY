!> ---------------------------------------
!> Copyright (C) 2013-2024, MPI-M, MPI-BGC
!>
!> Contact: icon-model.org
!> Authors: AUTHORS.md
!> See LICENSES/ for license information
!> SPDX-License-Identifier: BSD-3-Clause
!> ---------------------------------------
!>
MODULE mo_nix_init
#ifndef __NO_JSBACH__

  USE mo_kind,                   ONLY: wp
  USE mo_exception,              ONLY: message
  USE mo_jsb_control,            ONLY: debug_on

  USE mo_jsb_model_class,        ONLY: t_jsb_model, MODEL_QUINCY, MODEL_JSBACH
  USE mo_jsb_grid_class,         ONLY: t_jsb_grid
  USE mo_jsb_grid,               ONLY: Get_grid
  USE mo_jsb_class,              ONLY: Get_model
  USE mo_jsb_tile_class,         ONLY: t_jsb_tile_abstract
  USE mo_jsb_io_netcdf,          ONLY: t_input_file, jsb_netcdf_open_input
  USE mo_jsb_io,                 ONLY: missval

  USE mo_jsb_physical_constants, ONLY: rhoh2o

  ! NIX
  dsl4jsb_Use_processes NIX_
  dsl4jsb_Use_config(NIX_)
  dsl4jsb_Use_memory(NIX_)

  ! HYDRO
  dsl4jsb_Use_processes HYDRO_
  dsl4jsb_Use_config(HYDRO_)
  dsl4jsb_Use_memory(HYDRO_)


  IMPLICIT NONE
  PRIVATE

  PUBLIC :: nix_init

  CHARACTER(len=*), PARAMETER :: modname = 'mo_nix_init'

CONTAINS

  !
  !> Intialize nix process (after memory has been set up)
  !


  SUBROUTINE nix_init(tile)

    USE mo_jsb_class,             ONLY: Get_model
    USE mo_jsb_tile_class,        ONLY: t_jsb_tile_abstract

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_model), POINTER :: model
 
    ! NIX_
    dsl4jsb_Def_config(NIX_)
    dsl4jsb_Def_memory(NIX_)

    ! Pointers to variables in memory - NIX 
    dsl4jsb_Real3D_onDomain :: &
      & t_sn,                  &
      & dzm_sn          
    dsl4jsb_Int2D_onDomain ::  &
      & top_sn

    ! HYDRO_
    dsl4jsb_Def_config(HYDRO_)
    dsl4jsb_Def_memory(HYDRO_)

    ! Pointers to variables in memory - HYDRO
    dsl4jsb_Real2D_onDomain :: &
      & weq_snow,              &
      & snow_soil_dens               

    ! Local variables
    INTEGER :: nblks
    INTEGER :: ib        ! loop indices

    CHARACTER(len=*), PARAMETER :: routine = modname//':nix_init'

    model => get_model(tile%owner_model_id)

 
    ! Get required configs
    dsl4jsb_Get_config(NIX_)
    dsl4jsb_Get_config(HYDRO_)


    ! Get references to variables for domain - NIX_
    dsl4jsb_Get_memory(NIX_)

    dsl4jsb_Get_var3D_onDomain(NIX_, t_sn)
    dsl4jsb_Get_var3D_onDomain(NIX_, dzm_sn)

    dsl4jsb_Get_var2D_onDomain(NIX_, top_sn)


    ! Get references to variables for domain - HYDRO_
    dsl4jsb_Get_memory(HYDRO_)

    dsl4jsb_Get_var2D_onDomain(HYDRO_, weq_snow)
    dsl4jsb_Get_var2D_onDomain(HYDRO_, snow_soil_dens)


    !-----------------------------------------------------------------------
    ! Go   

    ! Get some array sizes for loops    
    nblks  = SIZE(t_sn,3)

    ! Initiate stratigraphy
    DO ib = 1, nblks 

      CALL nix_init_stratigraphy(                &
        & t_sn(:,:,ib),                          &
        & dzm_sn(:,:,ib),                        &
        & weq_snow(:,ib),                        &
        & snow_soil_dens(:,ib),                  &
        & dsl4jsb_Config(NIX_)%itype_nix_start   &
        & ) 

    ENDDO


  END SUBROUTINE nix_init

! =============================================================================
! + BEGIN subroutine: nix_init_stratigraphy
! =============================================================================

  SUBROUTINE nix_init_stratigraphy(    &
    & t_sn,                            &
    & dzm_sn,                          &
    & weq_snow,                        &
    & snow_soil_dens,                  &
    & itype_nix_start                  &
    & ) 

    !-----------------------------------------------------------------------
    !  Arguments

    REAL(wp), INTENT(INOUT) ::  t_sn(:,:)           !! snow layer (element) temperature [K]
    REAL(wp), INTENT(INOUT) ::  dzm_sn(:,:)         !! snow layer thickness             [m]

    REAL(wp), INTENT(IN)    ::  weq_snow(:)         !! snow water equivialent           [m]
    REAL(wp), INTENT(IN)    ::  snow_soil_dens(:)   !! snow density                     [kg/m3]    

    INTEGER,  INTENT(IN)    ::  itype_nix_start     !! type of initialization           [-]


    !-----------------------------------------------------------------------
    !  Local Variables

    INTEGER  :: ic, is
    INTEGER  :: nproma, nsnow

    REAL(KIND = wp) :: &
     & h_snow_nix
 
    CHARACTER(len=*), PARAMETER ::  &
      &  routine = 'mo_nix_init:nix_init_stratigtaphy'

    ! Assign local variables    
    nproma = SIZE(t_sn,1)
    nsnow  = SIZE(t_sn,2)


    SELECT CASE (itype_nix_start)
 
    !------------------------------------------------------------------------
    !  GO

    !----------------------------
    ! Case 1 - Hard cold start - All snow is wiped out.

    CASE (1)

    !----------------------------
    ! Case 2 - Soft cold start - Use mandatory analysis fields to initiate

    CASE (2)

      CALL message(TRIM(routine), "itype_nix_start = 2 was chosen, soft cold start")

      ! --------------------------
      ! Compute snow height and top layer index
      ! --------------------------

      DO ic = 1,nproma

        h_snow_nix = weq_snow(ic) * rhoh2o / snow_soil_dens(ic) 

       ENDDO


      DO is = 1, nsnow
        DO ic = 1, nproma

          t_sn(ic,is)   = 273.15_wp
          dzm_sn(ic,is) = h_snow_nix    
        
        END DO
      END DO

    END SELECT

  END SUBROUTINE nix_init_stratigraphy

! =============================================================================
! - END subroutine: nix_init_stratigraphy
! =============================================================================


#endif
END MODULE mo_nix_init








