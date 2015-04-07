!>
!! @brief Contains the memory class for carbon.
!!
!! @author
!!  Rainer Schneck, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Rainer Schneck (2015-04-30)
!!
MODULE mo_carbon_memory_class
#ifndef __NO_JSBACH__

  USE mo_kind, ONLY: wp

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: t_carbon_memory, t_soil_pools, t_veg_pools, t_canopy_layer,   &
            t_photosyn,  t_carbon_water_flux, copy_carbon_mem, add_carbon_mem

  TYPE t_soil_pools
    REAL(wp), POINTER, DIMENSION(:,:) ::                       &
                                        ! SIZE CLASS 1
                                        YCpool_acid_ag1,       & ! Yasso above ground litter-pool for acid soluble litter [mol(C)/m^2(canopy)]
                                        YCpool_water_ag1,      & ! Yasso above ground litter-pool for water soluble litter [mol(C)/m^2(canopy)]
                                        YCpool_ethanol_ag1,    & ! Yasso above ground litter-pool for ethanol soluble litter [mol(C)/m^2(canopy)]
                                        YCpool_nonsoluble_ag1, & ! Yasso above ground litter-pool for non-soluble litter [mol(C)/m^2(canopy)]
                                        YCpool_acid_bg1,       & ! Yasso below ground litter-pool for acid soluble litter [mol(C)/m^2(canopy)]
                                        YCpool_water_bg1,      & ! Yasso below ground litter-pool for water soluble litter [mol(C)/m^2(canopy)]
                                        YCpool_ethanol_bg1,    & ! Yasso below ground litter-pool for ethanol soluble litter [mol(C)/m^2(canopy)]
                                        YCpool_nonsoluble_bg1, & ! Yasso below ground litter-pool for non-soluble litter [mol(C)/m^2(canopy)]
                                        YCpool_humus_1,        & ! Yasso below ground litter-pool for slow C compartment [mol(C)^2/m(canopy)]
                                        ! SIZE CLASS 2
                                        YCpool_acid_ag2,       & ! Yasso above ground litter-pool for acid soluble litter
                                                                 ! [mol(C)/m^2(canopy)]
                                        YCpool_water_ag2,      & ! Yasso above ground litter-pool for water soluble litter
                                                                 ! [mol(C)/m^2(canopy)]
                                        YCpool_ethanol_ag2,    & ! Yasso above ground litter-pool for ethanol soluble litter
                                                                 ! [mol(C)/m^2(canopy)]
                                        YCpool_nonsoluble_ag2, & ! Yasso above ground litter-pool for non-soluble litter
                                                                 ! [mol(C)/m^2(canopy)]
                                        YCpool_acid_bg2,       & ! Yasso below ground litter-pool for acid soluble litter
                                                                 ! [mol(C)/m^2(canopy)]
                                        YCpool_water_bg2,      & ! Yasso below ground litter-pool for water soluble litter
                                                                 ! [mol(C)/m^2(canopy)]
                                        YCpool_ethanol_bg2,    & ! Yasso below ground litter-pool for ethanol soluble litter
                                                                 ! [mol(C)/m^2(canopy)]
                                        YCpool_nonsoluble_bg2, & ! Yasso below ground litter-pool for non-soluble litter
                                                                 ! [mol(C)/m^2(canopy)]
                                        YCpool_humus_2           ! Yasso below ground litter-pool for slow C compartment
                                                                 ! [mol(C)^2/m(canopy)]
  END TYPE t_soil_pools

  TYPE t_veg_pools
    REAL(wp), POINTER, DIMENSION(:,:) ::                       &
                                        Cpool_green,           & ! C-pool for leaves, fine roots, vegetative organs and other
                                                                 ! green (living) parts of vegetation [mol(C)/m^2(canopy)]
                                        Cpool_reserve,         & ! C-pool for carbohydrate reserve (sugars, starches) that 
                                                                 ! allows plants to survive bad times[mol(C)/m^2(canopy)]
                                        Cpool_woods,           & ! C-pool for stems, thick roots and other (dead) structural 
                                                                 ! material of living plants [mol(C)/m^2(canopy)]
                                        Cpool_litter_green_ag, & ! C-pool for above ground non-woody litter especially from 
                                                                 ! leaves [mol(C)/m^2(canopy)]
                                        Cpool_litter_green_bg, & ! C-pool for below ground non-woody litter especially from 
                                                                 ! fine roots [mol(C)/m^2(canopy)]
                                        Cpool_litter_wood_ag,  & ! C-pool for litter originating from above ground woody 
                                                                 ! parts of the plants [mol(C)/m^2(canopy)]
                                        Cpool_litter_wood_bg,  & ! C-pool for litter originating from below ground woody 
                                                                 ! parts of the plants [mol(C)/m^2(canopy)]
                                        Cpool_slow               ! C-pool for below ground organic material (coming from 
                                                                 ! the litter pools) [mol(C)/m^2(canopy)]
  END TYPE t_veg_pools

  
  TYPE t_photosyn
    REAL(wp), POINTER ::                                                &
                                        CO2_conc_leaf,                  & ! CO2 concentration inside leaf [MOL(CO2)/MOL(AIR)]
                                        CO2_conc_leaf_acc,              & ! CO2 concentration inside leaf [MOL(CO2)/MOL(AIR)]
                                                                          ! (mean value)
                                        CO2_conc_leaf_unlimited_acc,    & ! CO2 concentration inside leaf for no water stress
                                                                          ! [MOL(CO2)/MOL(AIR)] (mean value)
                                        canopy_conductance,             & ! Canopy conductance for water vapour 
                                                                          ! [m/s] (mean value)
                                        canopy_conductance_unlimited,   & ! Unlimited canopy conductance for water vapour
                                                                          ! [m/s] (mean value)
                                        gross_assimilation_acc,         & ! GROSS PHOTOSYNTHESIS [MOL(CO2) / M^2(ground) S] 
                                                                          ! (mean value)
                                        dark_respiration_acc,           & ! DARK RESPIRATION OF LEAF 
                                                                          ! [MOL(CO2) / M^2(ground) S] (mean value)
                                        max_carbox_rate,                & ! Maximum carboxylation rate (=VCmax)
                                                                          ! [MOL(CO2)/M^2 S]
                                        max_e_transport_rate,           & ! Maximum electron transport rate (=Jmax)
                                                                          ! (only C3 plants)[MOL(CO2)/M^2 S]
                                        carbox_rate,                    & ! Actual carboxylation rate (=JC)
                                                                          ! [MOL(CO2)/M^2 S]
                                        e_transport_rate,               & ! Actual electron transport rate (=JE)
                                                                          ! [MOL(CO2)/M^2 S]
                                        net_assimilation,               & ! Canopy net carbon assimilation
                                                                          ! [MOL(CO2)/m^2 s] (mean value)
                                        gross_assimilation,             & ! GROSS PHOTOSYNTHESIS at each time step
                                                                          ! [MOL(CO2) / M^2(ground) S]
                                        dark_respiration,               & ! DARK RESPIRATION OF LEAF at each time step
                                                                          ! [MOL(CO2) / M^2(ground) S]
                                        faPAR,                          & ! fraction of PAR absorbed by the canopy
                                        apar_acc,                       & ! PAR absorbed by the canopy
                                                                          ! [MOL PHOTONS/M^2 S] (mean value)
                                        par_acc                           ! incoming PAR [MOL PHOTONS/M^2 S] (mean value)
     REAL(wp), POINTER, DIMENSION(:) ::                                 &
                                        CO2_flux_net_acc,               & ! net CO2 flux to the atmosphere
                                        CO2_flux_herbivory_acc,         & ! CO2 flux to the atmosphere from grazing
                                        CO2_emission_landcover_change_acc, & ! CO2 flux to the atmosphere from land cover changes
                                        CO2_emission_harvest_acc,       & ! CO2 flux to the atmosphere from harvest
                                        CO2_flux_dynveg_acc               ! CO2 flux to the atmosphere from vegetation dynamics (disturbance?)    
  END TYPE t_photosyn

  
  TYPE t_canopy_layer
     INTEGER                             :: ncanopy                 ! Number of canopy layers
     INTEGER                             :: ntiles                  ! Number of tiles
     REAL(wp), POINTER, DIMENSION(:,:,:) ::                       & ! (nland, ncanopy, ntiles)
                                            lai,                  & ! LAI per canopy layer
                                            canopy_conductance,   & ! Canopy conductance per layer and per leaf area
                                            gross_assimilation,   & ! GROSS PHOTOSYNTHESIS [MOL(CO2) / M^2(leaf area) S]
                                            dark_respiration,     & ! DARK RESPIRATION OF LEAF [MOL(CO2) / M^2(leaf area) S]
                                            max_carbox_rate,      &
                                            max_e_transport_rate, &
                                            carbox_rate,          &
                                            e_transport_rate,     &
                                            faPAR                   ! Absorbed PAR per canopy layer [fraction of irradiance]
  END TYPE t_canopy_layer


  TYPE t_carbon_water_flux
    REAL(wp), POINTER :: &
    
  END TYPE t_carbon_water_flux

  TYPE t_carbon_memory
    TYPE(t_soil_pools),   POINTER :: soil_pools   => NULL()
    TYPE(t_veg_pools),    POINTER :: veg_pools    => NULL()
    TYPE(t_photosyn),     POINTER :: photosyn     => NULL()
    TYPE(t_canopy_layer), POINTER :: canopy_layer => NULL()
  END TYPE t_carbon_memory

! ! !   ! Note: Would be nicer to be able to overload the "=" operator, but that would be more
! ! !   ! complicated because we only want to copy parts of the arrays (ics:ice,iblk)
! ! !   INTERFACE copy_carbon_mem
! ! !     MODULE PROCEDURE copy_t_veg_pools
! ! !   END INTERFACE copy_carbon_mem
! ! ! 
! ! !   INTERFACE add_carbon_mem
! ! !     MODULE PROCEDURE add_t_veg_pools
! ! !   END INTERFACE add_carbon_mem

  CHARACTER(len=*), PARAMETER :: modname = 'mo_carbon_memory_class'

CONTAINS

! ! !   SUBROUTINE copy_t_veg_pools(in, out, iblk, ics, ice)
! ! ! 
! ! !     TYPE(t_veg_pools), INTENT(in)  :: in
! ! !     TYPE(t_veg_pools), INTENT(out) :: out
! ! !     INTEGER,                   INTENT(in)  :: iblk, ics, ice
! ! ! 
! ! !     out%t_srf         (ics:ice,    iblk) = in%t_srf         (ics:ice,    iblk)
! ! !     out%t_carbon_old     (ics:ice,    iblk) = in%t_carbon_old     (ics:ice,    iblk)
! ! !     out%t_carbon_unfilt  (ics:ice,    iblk) = in%t_carbon_unfilt  (ics:ice,    iblk)
! ! !     out%t_eff_srf     (ics:ice,    iblk) = in%t_eff_srf     (ics:ice,    iblk)
! ! !     out%s_srf         (ics:ice,    iblk) = in%s_srf         (ics:ice,    iblk)
! ! ! !!$    out%s_carbon_old     (ics:ice,    iblk) = in%s_carbon_old     (ics:ice,    iblk)
! ! !     out%qsat_srf      (ics:ice,    iblk) = in%qsat_srf      (ics:ice,    iblk)
! ! ! !!$    out%qsat_carbon_old  (ics:ice,    iblk) = in%qsat_carbon_old  (ics:ice,    iblk)
! ! !     out%tte_corr      (ics:ice,    iblk) = in%tte_corr      (ics:ice,    iblk)
! ! !     out%fact_q_air    (ics:ice,    iblk) = in%fact_q_air    (ics:ice,    iblk)
! ! !     out%fact_qsat_srf (ics:ice,    iblk) = in%fact_qsat_srf (ics:ice,    iblk)
! ! !     out%fact_qsat_trans_srf (ics:ice,    iblk) = in%fact_qsat_trans_srf (ics:ice,    iblk)
! ! ! 
! ! !   END SUBROUTINE copy_t_veg_pools
! ! ! 
! ! !   SUBROUTINE add_t_veg_pools(in1, in2, out, iblk, ics, ice)
! ! ! 
! ! !     TYPE(t_veg_pools), INTENT(in)  :: in1, in2
! ! !     TYPE(t_veg_pools), INTENT(out) :: out
! ! !     INTEGER,                   INTENT(in)  :: iblk, ics, ice
! ! ! 
! ! !     out%t_srf         (ics:ice,    iblk) = in1%t_srf         (ics:ice,    iblk) + in2%t_srf         (ics:ice,    iblk)
! ! !     out%t_carbon_old     (ics:ice,    iblk) = in1%t_carbon_old     (ics:ice,    iblk) + in2%t_carbon_old     (ics:ice,    iblk)
! ! !     out%t_carbon_unfilt  (ics:ice,    iblk) = in1%t_carbon_unfilt  (ics:ice,    iblk) + in2%t_carbon_unfilt  (ics:ice,    iblk)
! ! !     out%t_eff_srf     (ics:ice,    iblk) = in1%t_eff_srf     (ics:ice,    iblk) + in2%t_eff_srf     (ics:ice,    iblk)
! ! !     out%s_srf         (ics:ice,    iblk) = in1%s_srf         (ics:ice,    iblk) + in2%s_srf         (ics:ice,    iblk)
! ! ! !!$    out%s_carbon_old     (ics:ice,    iblk) = in1%s_carbon_old     (ics:ice,    iblk) + in2%s_carbon_old     (ics:ice,    iblk)
! ! !     out%qsat_srf      (ics:ice,    iblk) = in1%qsat_srf      (ics:ice,    iblk) + in2%qsat_srf      (ics:ice,    iblk)
! ! ! !!$    out%qsat_carbon_old  (ics:ice,    iblk) = in1%qsat_carbon_old  (ics:ice,    iblk) + in2%qsat_carbon_old  (ics:ice,    iblk)
! ! !     out%tte_corr      (ics:ice,    iblk) = in1%tte_corr      (ics:ice,    iblk) + in2%tte_corr      (ics:ice,    iblk)
! ! !     out%fact_q_air    (ics:ice,    iblk) = in1%fact_q_air    (ics:ice,    iblk) + in2%fact_q_air    (ics:ice,    iblk)
! ! !     out%fact_qsat_srf (ics:ice,    iblk) = in1%fact_qsat_srf (ics:ice,    iblk) + in2%fact_qsat_srf (ics:ice,    iblk)
! ! !     out%fact_qsat_trans_srf (ics:ice,  iblk) = in1%fact_qsat_trans_srf (ics:ice,  iblk) + in2%fact_qsat_trans_srf (ics:ice,  iblk)
! ! ! 
! ! !   END SUBROUTINE add_t_veg_pools

#endif
END MODULE mo_carbon_memory_class
