!>
!! @brief Initialization of the the phenology memory.
!!
!! @author
!!  Rainer Schneck, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Rainer Schneck (2015-09-29)
!!



MODULE mo_pheno_init
#ifndef __NO_JSBACH__
#include "dsl_jsbach.inc"

  USE mo_kind,                ONLY: wp
  USE mo_exception,           ONLY: message

  USE mo_jsb_model_class,     ONLY: t_jsb_model  
  USE mo_jsb_class,           ONLY: get_model
  USE mo_jsb_tile_class,      ONLY: t_jsb_tile_abstract
  USE mo_jsb_io_netcdf,       ONLY: t_input_file, jsb_netcdf_open_input


  dsl4jsb_Use_processes PHENO_
  dsl4jsb_Use_config(PHENO_)
  dsl4jsb_Use_memory(PHENO_)
    
  IMPLICIT NONE
  PRIVATE
  
  PUBLIC :: pheno_init

  CHARACTER(len=*), PARAMETER :: modname = 'mo_pheno_init'

CONTAINS

  SUBROUTINE pheno_init(tile)
    
    ! Arguments
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    
    ! Local variables
    TYPE(t_jsb_model), POINTER  :: model

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':pheno_init'

    ! Declare pointers for process configuration and memory

    ! Declare pointers for variables in memory

    ! ---------------------------
    ! Go
    
    model => get_model(tile%owner_model_id)

    CALL pheno_init_bc(tile)
    CALL pheno_init_ic(tile)

  END SUBROUTINE pheno_init

  SUBROUTINE pheno_init_bc(tile)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile

    TYPE(t_jsb_model),       POINTER :: model

    dsl4jsb_Def_config(PHENO_)
    dsl4jsb_Def_memory(PHENO_)

    REAL(wp), POINTER ::        &
      & lai_mon(:,:,:)

    TYPE(t_input_file) :: input_file

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':pheno_init_bc'

    model => get_model(tile%owner_model_id)

    dsl4jsb_Get_config(PHENO_)
    dsl4jsb_Get_memory(PHENO_)

    CALL message(TRIM(routine), 'Setting boundary conditions of pheno memory for tile '// &
      &                          TRIM(tile%name)//' from '//TRIM(dsl4jsb_Config(PHENO_)%bc_filename))

    input_file = jsb_netcdf_open_input(TRIM(dsl4jsb_Config(PHENO_)%bc_filename), model%grid_id)

    ! Leaf area index
    ! Get monthly LAI climatology
    lai_mon => input_file%Read_2d_time(    &
      & variable_name='lai_clim')
    lai_mon(:,:,:) = MERGE(lai_mon(:,:,:), 0._wp, lai_mon(:,:,:) > 0._wp)
    pheno__mem%lai_mon(:,:,1:12) = lai_mon(:,:,1:12)
    pheno__mem%lai_mon(:,:,0   ) = lai_mon(:,:,12  )
    pheno__mem%lai_mon(:,:,13  ) = lai_mon(:,:,1   )
    DEALLOCATE(lai_mon)

    CALL input_file%Close()

  END SUBROUTINE pheno_init_bc

  SUBROUTINE pheno_init_ic(tile)

    USE mo_jsb_time,    ONLY: get_time_interpolation_weights

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile

    TYPE(t_jsb_model),       POINTER :: model

    dsl4jsb_Def_config(PHENO_)
    dsl4jsb_Def_memory(PHENO_)

    REAL(wp) :: wgt1, wgt2
    INTEGER  :: nmw1,nmw2

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':pheno_init_ic'

    model => get_model(tile%owner_model_id)

    dsl4jsb_Get_config(PHENO_)
    dsl4jsb_Get_memory(PHENO_)

    CALL message(TRIM(routine), 'Setting initial conditions of pheno memory (lai) for tile '// &
      &                          TRIM(tile%name)//' from '//TRIM(dsl4jsb_Config(PHENO_)%bc_filename))

    ! Set initial leaf area index from climatology (read in pheno_init_bc)
    CALL get_time_interpolation_weights(wgt1, wgt2, nmw1, nmw2)
    dsl4jsb_var2D_onDomain(PHENO_, lai) = wgt1 * pheno__mem%lai_mon(:,:,nmw1) + wgt2 * pheno__mem%lai_mon(:,:,nmw2)

    !dsl4jsb_var2D_onDomain(HYDRO_, w_skin_can) = dsl4jsb_Config(VEG_)%wsr_capacity * dsl4jsb_var2D_onDomain(PHENO_, lai) 


    !CALL message(TRIM(routine), 'Setting initial conditions of pheno memory for tile '// &
    !  &                          TRIM(tile%name)//' from '//TRIM(dsl4jsb_Config(PHENO_)%ic_filename))

  END SUBROUTINE pheno_init_ic
  
#endif
END MODULE mo_pheno_init
