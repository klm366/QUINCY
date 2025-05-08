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

      CALL nix_init_stratigraphy(      &
        & t_sn(:,:,ib),                &
        & dzm_sn(:,:,ib),              &
        & weq_snow(:,ib),              &
        & snow_soil_dens(:,ib)         &
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
    & snow_soil_dens                   &
    ) 

    !-----------------------------------------------------------------------
    !  Arguments

    REAL(wp), INTENT(INOUT) ::  t_sn(:,:)           !! snow layer (element) temperature [K]
    REAL(wp), INTENT(INOUT) ::  dzm_sn(:,:)         !! snow layer thickness             [m]

    REAL(wp), INTENT(IN)    ::  weq_snow(:)         !! snow water equivialent           [m]
    REAL(wp), INTENT(IN)    ::  snow_soil_dens(:)   !! snow density                     [kg/m3]    


    !-----------------------------------------------------------------------
    !  Local Variables

    INTEGER  :: ic, is
    INTEGER  :: nproma, nsnow

    REAL(KIND = wp) :: &
     & h_snow_nix


    ! Assign local variables    
    nproma = SIZE(t_sn,1)
    nsnow  = SIZE(t_sn,2)

 
    !------------------------------
    !  GO

    DO ic = 1,nproma

      h_snow_nix = weq_snow(ic) * rhoh2o / snow_soil_dens(ic) 

!IF(h_snow_nix .GT. 0.25_wp) THEN
!
!  WRITE(*,*) 'H_SNOW:', h_snow_nix
!
!ENDIF

    ENDDO


    DO is = 1, nsnow
      DO ic = 1, nproma

        t_sn(ic,is)   = 273.15_wp
        dzm_sn(ic,is) = 0.0_wp    
        

      END DO
    END DO


  END SUBROUTINE nix_init_stratigraphy

! =============================================================================
! - END subroutine: nix_init_stratigraphy
! =============================================================================


#endif
END MODULE mo_nix_init








