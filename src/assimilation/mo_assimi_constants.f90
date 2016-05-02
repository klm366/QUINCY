!>
!! @brief Contains constants for the assimilation processes
!!
!! @author
!!  Rainer Schneck, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Rainer Schneck (2015-05-11)
!!
MODULE mo_assimi_constants
#ifndef __NO_JSBACH__

  USE mo_kind, ONLY: wp
  USE mo_physical_constants,       ONLY: amd, amco2

  IMPLICIT NONE
  PUBLIC

  REAL(wp), PARAMETER ::                  &
         ! C3 PLANTS: FARQHUAR, G.D., S. VON CAEMMERER AND J.A. BERRY, 1980.           
         !            A BIOCHEMICAL MODEL OF PHOTOYNTHESIS IN LEAVES OF C3 SPECIES.    
         !            PLANTA 149, 78-90.                                               
                          ALPHA = 0.28,    & ! EFFICIENCY OF OF PHOTON CAPTURE
                          OX    = 0.21,    & ! OXYGEN CONCENTRATION [MOL(O2) / MOL(AIR)] 
                          KC0   = 460.E-6, & ! MICHAELIS-MENTEN CONSTANT FOR CO2 AT 25C [MOL(CO2) / MOL(AIR)] 
                          KO0   = 330.E-3, & ! MICHAELIS-MENTEN CONSTANT FOR O2 AT 25C [MOL(O2) / MOL(AIR)] 
                          EC    = 59356.,  & ! ACTIVATION ENERGY FOR KC [J / MOL]
                          EO    = 35948.,  & ! ACTIVATION ENERGY FOR KO [J / MOL]
                          EV    = 58520.,  & ! ACTIVATION ENERGY FOR VCMAX [J / MOL]
                          ER    = 45000.,  & ! ACTIVATION ENERGY FOR DARK RESPIRATION [J / MOL]
                          EK    = 50967.,  & !  = Q10=2 (Collatz et al. 1992)
                          FRDC3 = 0.011,   & ! RATIO OF DARK RESPIRATION TO "PVM" AT 25C for C3
         ! C4 PLANTS: COLLATZ, G.J., M. RIBAS-CARBO AND J.A. BERRY, 1992.              
         !            COUPLED PHOTOSYNTHESIS-STOMATAL CONDUCTANCE MODEL FOR LEAVES     
         !            OF C4 PLANTS. AUST. J. PLANT PHYSIOL. 19, 519-538.               
                          FRDC4 = 0.031,   & !RATIO OF DARK RESPIRATION TO "PVM" AT 25C for C4
                          ALC4  = 0.04,    & !EFFECTIVE QUANTUM EFFICIENCY
                          THETA = 0.83,    & !CURVATURE PARAMETER
                          minOfMaxCarboxrate = 1.0e-12, & ! Minimum of maximum carboxylation rate [10^(-6) mol/(m^2 s)].
         ! both                          
                          minStomaConductance = 0.0, & ! Minimum stomatal conductance [mol H2O /(m^2 s) ??]
         ! Factors that relates leaf internal CO2-concentration to CO2-concentration of ambient air: 
                          FCI1C3        = 0.87, & ! For C3 plants
                          FCI1C4        = 0.67, & ! For C4 plants
         ! Factors needed to calculate NPP_rate:                  
                          f_aut_leaf    = 0.40, & ! leaf fraction of plant-total (autotrophic) respiration
                          cCost         = 1.25    ! relative costs (measured in carbon) to produce 1 unit of carbon                          

  REAL(wp), PARAMETER ::                                     &
                         molarMassDryAir = amd   * 1.e-3_wp, &  ! 28.970e-3_wp. Mass [kg] per amount of substance [mol] of dry air
                         molarMassCO2    = amco2 * 1.e-3_wp    ! 44.0095e-3_wp. Mass [kg] per amount of substance [mol] of CO2                          
#endif
END MODULE mo_assimi_constants
