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
  USE mo_pheno_memory_class,  ONLY: t_pheno_memory
  USE mo_jsb_io_netcdf,       ONLY: t_input_file


  dsl4jsb_Use_processes PHENO
  dsl4jsb_Use_memory(PHENO)
    
  IMPLICIT NONE
  PRIVATE
  
  PUBLIC :: pheno_init

  CHARACTER(len=*), PARAMETER :: modname = 'mo_pheno_init'

CONTAINS

  SUBROUTINE pheno_init(tile)
    
    ! Declarations
    USE mo_pheno_constants, ONLY: PhenoParam
    USE mo_jsb_time,        ONLY: timeStep_in_days

    ! Arguments
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    
    ! Local variables
    TYPE(t_jsb_model), POINTER  :: model
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':pheno_init'

    ! Declare pointers for process configuration and memory
    dsl4jsb_Def_memory(PHENO)

    ! Declare pointers for variables in memory
    REAL(wp), POINTER, DIMENSION(:,:) :: F_pseudo_soil_temp, &
                                         N_pseudo_soil_temp

    ! ---------------------------
    ! Go
    
    dsl4jsb_Get_memory(PHENO)

    model => get_model(tile%owner_model_id)

    ! R: Jetzt in seb_init
    !dsl4jsb_Get_var2D_onDomain(PHENO, N_pseudo_soil_temp) ! OUT
    !dsl4jsb_Get_var2D_onDomain(PHENO, F_pseudo_soil_temp) ! OUT
    ! For computation of the pseudo soil temperature
    !F_pseudo_soil_temp = EXP(-timeStep_in_days()/PhenoParam%EG_SG%tau_pseudo_soil)
    !N_pseudo_soil_temp = 1._wp / (1._wp - F_pseudo_soil_temp)

    CALL message('init_phenology()','Initialization of PHENOLOGY finished.')

  END SUBROUTINE pheno_init

  
  
#endif
END MODULE mo_pheno_init
