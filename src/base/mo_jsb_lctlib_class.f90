!>
!! This module provides the data from the landcover library file, which contains additional information on the
!! landcover data used in a particular run of JSBACH.
!! The name of the landcover library file is specified in the JSBACH configuration file (run.def) under the
!! keyword "LCT_FILE".
!!
!! Authors: Reiner Schnur, 7/2003
!!          Christian H. Reick, 8/2003
!!
!! @par Revision History
!! Julia Nabel (2018-01-10): added params for optional functionality in PHENO_: forest regrowth (l_forestRegrowth)
!!                           (derive maxLAI from biomass via allometric relationships)
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/license .
!!

MODULE mo_jsb_lctlib_class

  USE mo_jsb_parallel, ONLY: my_process_is_stdio, my_process_is_mpi_parallel, p_io, mpi_comm, p_bcast
  USE mo_kind,         ONLY: wp
  USE mo_exception,    ONLY: message_text, finish
  USE mo_util,         ONLY: int2string

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: t_lctlib_element, Read_lctlib

  TYPE t_lctlib_element
     ! --- general parameters ------------------------------------------------------------------------------------------------------
    INTEGER  :: LctNumber            !! Unique number for each landcover type (1, 2, 3, ..)
    CHARACTER(LEN=16)  :: LctName    !! Names of landcover types
    INTEGER  :: LandcoverClass       !! Landcover classes (not to be confused with landcover types!):
                                                 !! .. 0: Bare soil; 1:Glacier; 2: Lake;
                                                 !! .. 3: Natural forest; 4: Natural Grassland;
                                                 !! .. 5: Other natural vegetation; 6: crops; 7: pastures
    LOGICAL  :: NaturalVegFlag       !! "true" for natural vegetation landcover types, i.e. non-agricultural vegetation
    LOGICAL  :: ForestFlag           !! "true" for forest landcover types
    LOGICAL  :: GrassFlag            !! "True" for natural graslands
    LOGICAL  :: CropFlag             !! "true" for croplands
    LOGICAL  :: PastureFlag          !! "true" for landcover of type "pasture"
    LOGICAL  :: LakeFlag             !! "true" for a land surface class that signifies lake or ocean
    LOGICAL  :: GlacierFlag          !! "true" for a land surface class that signifies glaciers
    LOGICAL  :: BareSoilFlag         !! "true" if none of the foregoing landcover classes

    ! --- Albedo ------------------------------------------------------------------------------------------------------------------
    REAL(wp) :: AlbedoSnowVisMin      !! Minimum snow albedo in the visible range
    REAL(wp) :: AlbedoSnowVisMax      !! Maximum snow albedo in the visible range
    REAL(wp) :: AlbedoSnowNirMin      !! Minimum snow albedo in the NIR range
    REAL(wp) :: AlbedoSnowNirMax      !! Maximum snow albedo in the NIR range
    REAL(wp) :: AlbedoSnowMin         !! Minimum snow albedo
    REAL(wp) :: AlbedoSnowMax         !! Maximum snow albedo
    REAL(wp) :: AlbedoCanopyVIS       !! Albedo of the canopy (vegetation) in the visible range
    REAL(wp) :: AlbedoCanopyNIR       !! Albedo of the canopy (vegetation) in the NIR range
    REAL(wp) :: AlbedoLitterVis       !! Albedo of the leaf litter in the visible range
    REAL(wp) :: AlbedoLitterNir       !! Albedo of the leaf litter in the NIR range

    ! --- parameters used by bethy ------------------------------------------------------------------------------------------------
    LOGICAL  :: NitrogenScalingFlag   !! Indicates, whether nitrogen scaling shall be applied to that vegetation type
    LOGICAL  :: C4flag                !! Photosynthetic pathway: C4=.true or C3=.false.
    REAL(wp) :: CarboxRate            !! Maximum carboxilation rate at 25 degrees Celsius [1.E-6 * Mol(CO2)/m^2/s] ...
                                                  !! ... (Table 2.6 in Knorr)
    REAL(wp) :: ETransport            !! Maximum electron transport rate at 25 degrees Celsius [1.E-6 * Mol/m^2/s] ...
                                                  !! ... (Table 2.6 in Knorr)
    REAL(wp) :: VegHeight             !! Typical height of the vegetation classes [m]
    REAL(wp) :: VegRoughness          !! Roughness length of the vegetation classes [m]
    REAL(wp) :: MinVegRoughness       !! Minimal roughness length of the vegetation classes (LAI = 0) [m]
    REAL(wp) :: MaxVegRoughness       !! Maximal roughness length of the vegetation classes (LAI = LAI_max) [m]
    ! --- Parameters used in LoGroP Phenology scheme --------------------------------------------------------------------------
    INTEGER  :: PhenologyType         !! Phenology type (only for natural vegetation):
                                                  !! ... none: 0; evergreen: 1; summergreen: 2; raingreen: 3; grasses: 4; crops: 5
    REAL(wp) :: MaxLAI                !! Upper LAI boundary for LoGoP-Scheme (phenology) in [m2/m2]
    REAL(wp) :: specificLeafArea_C    !! Carbon content per leaf area in [mol(Carbon)/m^2(leaf)]

    ! ++ additional parameters required when l_forestRegrowth = true  --------------------------------------------------------
    REAL(wp) :: alpha_nr_ind          !! Parameter of self-thinning relationship
    REAL(wp) :: beta_nr_ind           !! Parameter of self-thinning relationship
    REAL(wp) :: alpha_leaf            !! Parameter of relationship between total biomass per ind and leaf biomass
    REAL(wp) :: beta_leaf             !! Parameter of relationship between total biomass per ind and leaf biomass

     ! --- Parameters used in Knorr Phenology scheme ---------------------------------------------------------------------------
    REAL(wp) :: knorr_Tau_w           !! Time before leaf shedding  [days]
    REAL(wp) :: knorr_T_phi           !! Temperature trigger for leaf growth [deg C]
    REAL(wp) :: knorr_T_r             !! "Spread" (sigma) of T_phi [deg C]
    REAL(wp) :: knorr_Day_c           !! Day-length at leaf shedding [hours]
    REAL(wp) :: knorr_Day_r           !! "Spread" (sigma) of Day_c [hours]
    REAL(wp) :: knorr_k_l             !! Inverse of leaf longevity [days-1]
    REAL(wp) :: knorr_leaf_growth_rate !! Initial leaf growth rate [days-1]
    REAL(wp) :: knorr_max_lai         !! Maximum LAI

    ! --- Parameters used in the carbon balance model -----------------------------------------------------------------------------
    REAL(wp) :: reserveC2leafC        !! Ratio of opt. C content of the reserve pool to the opt. C content of leaves
    REAL(wp) :: fract_npp_2_woodPool   !! Maximum fraction of NPP put into the wood pool of the carbon balance model
    REAL(wp) :: fract_npp_2_reservePool!! Optimal fraction of NPP put into the reserve pool of the carbon balance model
    REAL(wp) :: fract_npp_2_exudates   !! Optimal fraction of NPP put into the root exudates of the carbon balance model
    REAL(wp) :: fract_green_2_herbivory!! Fraction of NPP that is grazed by herbivores
    REAL(wp) :: tau_c_litter_leaf !! Time constant by which leaf litter pool is depreciated  [days]
    REAL(wp) :: tau_c_litter_wood !! Time constant by which woody litter pool is depreciated  [days]
    REAL(wp) :: LAI_shed_constant     !!
    REAL(wp) :: fract_C_litter_green2atmos !! Fraction of heterotrophic loss of green litter pools emitted to the
                                                      !! .. atmosphere (rest enters slow pool)
    REAL(wp) :: Max_C_content_woods   !! Maximum carbon content in the woody parts of plants [mol(C)/m^2(canopy)] (for
                                                  !! .. forests this is closely related to the yield, usually measured in
                                                  !! .. m^3(wood)/hectar)
    REAL(wp) :: ClumpinessFactor      !! Factor to calculate vegetation clumpiness:
                                                  !! Former in JSBACH3:   veg_ratio=veg_ratio_max*(1-exp(-LAI_max/ClumpinessFactor))
                                                  !! converting natural PFTs into agricultural PFTs
    REAL(wp) :: fract_wood_2_onSite    !! Fraction of wood pool to anthropogenically controlled onSite pool wenn
                                                  !! converting natural PFTs into agricultural PFTs
    REAL(wp) :: fract_wood_2_paper     !! Fraction of wood pool to paper (intermediately longlived) pool wenn...
    REAL(wp) :: fract_wood_2_construction !! Fract of wood pool to constructions (e.g. houses, furniture = longlived)
                                                  !!   pool when converting natural PFTs into agricultural PFTs

    ! --- Parameters used in the carbon balance model and the dynamic vegetation --------------------------------------------------
    REAL(wp) :: tau_c_woods           !! PFT-specific time scale for woody (lignified) plant tissue

    ! --- Parameters used by the dynamic vegetation -------------------------------------------------------------------------------
    LOGICAL  :: dynamic_PFT               !! indicates those PFTs which shall take part in the vegetation dynamics
    LOGICAL  :: woody_PFT                 !! indicates those PFTs which are of woody type (in contrast to grasses)
    LOGICAL  :: pasture_PFT               !! indicates those PFTs which are pasture
    REAL(wp) :: bclimit_min_cold_mmtemp   !! PFT-specific minimum coldest monthly mean temperature
    REAL(wp) :: bclimit_max_cold_mmtemp   !! PFT-specific maximum coldest monthly mean temperature
    REAL(wp) :: bclimit_max_warm_mmtemp   !! PFT-specific upper limit of warmest-month temperature
    REAL(wp) :: bclimit_min_temprange     !! PFT-specific 20-year average min warmest - coldest month temperature range
    REAL(wp) :: bclimit_min_gdd           !! PFT-specific minimum growing degree days (at or above 5 deg C)

    ! --- Parameters used in SPITFIRE (mo_disturbance_thonicke.f90) ---------------------------------------------------------------
    REAL(wp) :: moist_extinction          !! PFT specific moisture of extinction
    REAL(wp) :: fuel_dens                 !! fuel bulk density
    REAL(wp) :: flame_length_f            !! f parameter for flame length(scorch height) parameter 45 in lpj
    REAL(wp) :: crown_length              !! crown length parameter see table 1 in thonicke et al. 2010
    REAL(wp) :: bark_par1                 !! bark thickness parameter 1 see table 1 and eq. 21 in thonicke et al. 2010
    REAL(wp) :: bark_par2                 !! bark thickness parameter 2 see table 1 and eq. 21 in thonicke et al. 2010
    REAL(wp) :: RCK                       !! resistance factor to crown damage, tab1 and eq22 in thonicke et al. 2010
    REAL(wp) :: mort_prob                 !! mortality probability  see table 1 and eq. 22 in thonicke et al. 2010

    ! --- Parameters for the climate buffer ---------------------------------------------------------------------------------------
    REAL(wp) :: gdd_base                  !! PFT-specific GDD base
    REAL(wp) :: upper_tlim                !! PFT-specific base to calculate GDD_upper_tlim

    ! --- Parameters used by the yasso soil carbon model -------------------------------------------------------------------------
    REAL(wp) :: LitVeg_coef(5)       !! Yasso coefficient for separating litter into chemical pools
    REAL(wp) :: LeafLit_coef(5)      !! Yasso coefficient for separating leaf litter into chemical pools
    REAL(wp) :: WoodLit_coef(5)      !! Yasso coefficient for separating woody litter into chemical pools
    REAL(wp) :: WoodLitterSize       !! Yasso: size of the litter

    ! --- Other parameters --------------------------------------------------------------------------------------------------------
    REAL(wp) :: CanopyResistanceMin  !! Minimum canopy resistance (optional, only used in VIC scheme
                                                 !! and if BETHY is not used)
    REAL(wp) :: StemArea             !! Area of stems and branches of woody plants
  END TYPE t_lctlib_element

CONTAINS

  FUNCTION Read_lctlib(lctlib_file_name) RESULT(return_value)
    !
    ! Get landcover library, i.e. lookup table for landcover types
    !
    ! The structure of the lookup table in the input file is as follows. It can contain only three types of lines in arbitrary
    ! order:
    !
    !     Comment lines: These contain before a '#' only blanks
    !       Blank lines: These contain only blanks
    !        Data lines: The first nonblank characters form a keyword which is follwed by the data (separated by blanks)
    !
    ! Keywords are the same as the names of the lctlib components, with one exception:
    !
    !        NLCT : This is the keyword after which the number of landcover types, i.e. the number of data columns
    !               in the file is expected. This number has to be identical with the number of landcover types
    !               from the landcover data file (keyword: LCTLIB_FILE in the JSBACH configuration file).
    !
    ! After all other keywords NLCT columns of data are expected.
    !
    ! Example:
    !            # ---- LANDCOVER LIBRARY FILE -------
    !            NLCT 3
    !            LctNumber               3 5 9
    !            # Phenology types: 0= none 1=summergreen, 2=evergreen, 3=raingreen, 4=grasses
    !            PhenologyType          2 2 4
    !            # C4flag: 0=C3, 1=C4
    !            C4flag                 0 1 1
    !
    ! The first string on each line (except first line and comment lines) must correspond to the name of the lctlib component
    ! (case doesn't matter)
    ! For lctlib components of type LOGICAL use 0/1 in the file to indicate .FALSE./.TRUE.
    ! For lctlib components of type INTEGER the numbers on the line must be integer values
    ! For lctlib components of type REAL the numbers on the line can be either integer or floating point
    ! Comments can appear on any line, everything to the right of and including "#" is disregarded
    !
    ! The file can contain more keywords than needed --- therefore the same file can be used by several model components.
    !
    USE mo_util_string, ONLY: tolower

    CHARACTER(len=*), INTENT(in) :: lctlib_file_name
    TYPE(t_lctlib_element), POINTER    :: return_value(:)

    TYPE(t_lctlib_element), ALLOCATABLE :: lctlib(:)

    INTEGER, PARAMETER            :: lctlib_file_unit = 66

    CHARACTER(len=30)  :: key
    CHARACTER(len=256) :: line
    INTEGER            :: pos_comment, read_status
    INTEGER            :: pos,length
    CHARACTER(len=2)   :: blank_set = " "//achar(9) !! the blank characters: BLANK and TAB
    INTEGER,ALLOCATABLE:: itmp(:)                   !! temporary array used for input of logicals
    INTEGER            :: apu  ! THT 30.10.2012

    LOGICAL            :: exist_LctNumber              = .FALSE.
    LOGICAL            :: exist_LctName                = .FALSE.
    LOGICAL            :: exist_LandcoverClass         = .FALSE.
    LOGICAL            :: exist_AlbedoSnowVisMin       = .FALSE.
    LOGICAL            :: exist_AlbedoSnowVisMax       = .FALSE.
    LOGICAL            :: exist_AlbedoSnowNirMin       = .FALSE.
    LOGICAL            :: exist_AlbedoSnowNirMax       = .FALSE.
    LOGICAL            :: exist_AlbedoSnowMin          = .FALSE.
    LOGICAL            :: exist_AlbedoSnowMax          = .FALSE.
    LOGICAL            :: exist_AlbedoCanopyVIS        = .FALSE.
    LOGICAL            :: exist_AlbedoCanopyNIR        = .FALSE.
    LOGICAL            :: exist_AlbedoLitterVis        = .FALSE.
    LOGICAL            :: exist_AlbedoLitterNir        = .FALSE.
    LOGICAL            :: exist_nitrogenScalingFlag    = .FALSE.
    LOGICAL            :: exist_C4flag                 = .FALSE.
    LOGICAL            :: exist_CarboxRate             = .FALSE.
    LOGICAL            :: exist_ETransport             = .FALSE.
    LOGICAL            :: exist_VegHeight              = .FALSE.
    LOGICAL            :: exist_VegRoughness           = .FALSE.
    LOGICAL            :: exist_MinVegRoughness        = .FALSE.
    LOGICAL            :: exist_MaxVegRoughness        = .FALSE.
    LOGICAL            :: exist_PhenologyType          = .FALSE.
    LOGICAL            :: exist_CanopyResistanceMin    = .FALSE.
    LOGICAL            :: exist_MaxLAI                 = .FALSE.
    LOGICAL            :: exist_StemArea               = .FALSE.
    LOGICAL            :: exist_specificLeafArea_C     = .FALSE.
    LOGICAL            :: exist_alpha_nr_ind           = .FALSE.
    LOGICAL            :: exist_beta_nr_ind            = .FALSE.
    LOGICAL            :: exist_alpha_leaf             = .FALSE.
    LOGICAL            :: exist_beta_leaf              = .FALSE.
    LOGICAL            :: exist_knorr_Tau_w            = .FALSE.
    LOGICAL            :: exist_knorr_T_phi            = .FALSE.
    LOGICAL            :: exist_knorr_T_r              = .FALSE.
    LOGICAL            :: exist_knorr_Day_c            = .FALSE.
    LOGICAL            :: exist_knorr_Day_r            = .FALSE.
    LOGICAL            :: exist_knorr_k_l              = .FALSE.
    LOGICAL            :: exist_knorr_leaf_growth_rate = .FALSE.
    LOGICAL            :: exist_knorr_max_lai          = .FALSE.
    LOGICAL            :: exist_reserveC2leafC         = .FALSE.
    LOGICAL            :: exist_fract_npp_2_woodPool    = .FALSE.
    LOGICAL            :: exist_fract_npp_2_reservePool = .FALSE.
    LOGICAL            :: exist_fract_npp_2_exudates    = .FALSE.
    LOGICAL            :: exist_fract_green_2_herbivory = .FALSE.
    LOGICAL            :: exist_tau_c_litter_leaf  = .FALSE.
    LOGICAL            :: exist_tau_c_litter_wood  = .FALSE.
    LOGICAL            :: exist_tau_c_woods        = .FALSE.
    LOGICAL            :: exist_LAI_shed_constant      = .FALSE.
    LOGICAL            :: exist_fract_C_litter_green2atmos = .FALSE.
    LOGICAL            :: exist_Max_C_content_woods       = .FALSE.
    LOGICAL            :: exist_ClumpinessFactor          = .FALSE.
    LOGICAL            :: exists_dynamic_PFT              = .FALSE.
    LOGICAL            :: exists_woody_PFT                = .FALSE.
    LOGICAL            :: exists_pasture_PFT              = .FALSE.
    LOGICAL            :: exists_bclimit_min_cold_mmtemp  = .FALSE.
    LOGICAL            :: exists_bclimit_max_cold_mmtemp  = .FALSE.
    LOGICAL            :: exists_bclimit_max_warm_mmtemp  = .FALSE.
    LOGICAL            :: exists_bclimit_min_temprange    = .FALSE.
    LOGICAL            :: exists_bclimit_min_gdd          = .FALSE.
    LOGICAL            :: exists_gdd_base                 = .FALSE.
    LOGICAL            :: exists_upper_tlim               = .FALSE.
    LOGICAL            :: exists_fract_wood_2_onSite       = .FALSE.
    LOGICAL            :: exists_fract_wood_2_paper        = .FALSE.
    LOGICAL            :: exists_fract_wood_2_construction = .FALSE.
    LOGICAL            :: exists_moist_extinction      = .FALSE.
    LOGICAL            :: exists_fuel_dens             = .FALSE.
    LOGICAL            :: exists_flame_length_f        = .FALSE.
    LOGICAL            :: exists_crown_length          = .FALSE.
    LOGICAL            :: exists_bark_par1             = .FALSE.
    LOGICAL            :: exists_bark_par2             = .FALSE.
    LOGICAL            :: exists_RCK                   = .FALSE.
    LOGICAL            :: exists_mort_prob             = .FALSE.
    LOGICAL            :: exists_woodlittersize_coef      = .FALSE.
    LOGICAL            :: exists_woodlit_coef             = .FALSE.
    LOGICAL            :: exists_leaflit_coef             = .FALSE.
    LOGICAL            :: exists_litveg_coef              = .FALSE.
    INTEGER            :: i

    INTEGER :: nlct, npft

    ! Determine name of landcover library file

    ! Read lctlib data and eventually allocate also memory for lctlib data on io-processor

    IF (my_process_is_stdio()) THEN

       OPEN(unit=lctlib_file_unit, file=lctlib_file_name, form='FORMATTED', status='OLD', iostat=read_status)
       IF (read_status /= 0) CALL finish('init_lctlib','Error opening landcover library file')

       ! --- find keyword NLCT in landcover library file

       DO
          READ(lctlib_file_unit,'(A256)', IOSTAT=read_status) line
          IF (read_status /= 0) THEN
             CALL finish('init_lctlib','No keyword NLCT found in land cover library file '//TRIM(lctlib_file_name))
          END IF
          ! Look for comment
          pos_comment = SCAN(line,"#")                 ! Start position of comment, 0 if none present
          IF (pos_comment > 1) THEN
             line = TRIM(ADJUSTL(line(1:pos_comment))) ! Disregard everything to the right of "#", ..
             ! .. adjust to the left and trim to the right
          ELSE IF (pos_comment == 1) THEN
             CYCLE                                     ! Disregard whole line
          ELSE
             line = TRIM(ADJUSTL(line))                ! Only disregard blanks to the left and right
          ENDIF
          length=LEN_TRIM(line)
          IF(length== 0) CYCLE                 ! Line is empty

          pos = SCAN(line,blank_set)                   ! Position of first blank character
          IF(pos == 0)   CALL finish('init_lctlib',"Wrong syntax in lctlib definitions")

          READ(line(1:pos-1),'(A)') key

          IF(tolower(TRIM(key)) == 'nlct') THEN
             READ(line(pos:length),*,IOSTAT=read_status) nlct
             IF (read_status /= 0) THEN
                CALL finish('init_lctlib','Could not read number of landcover types (keyword: NLCT) from '//TRIM(lctlib_file_name))
             END IF
             EXIT ! found number of landcover types in landcover library file --- continue after loop
          END IF
       END DO

    END IF

    IF (my_process_is_mpi_parallel()) CALL p_bcast(nlct, p_io, mpi_comm)

    ALLOCATE(lctlib(nlct))

    IF (my_process_is_stdio()) THEN

       REWIND(unit=lctlib_file_unit) ! Go back to beginning of landcover library file

       ALLOCATE(itmp(nlct))

       DO
          READ(lctlib_file_unit,'(A256)', IOSTAT=read_status) line
          IF (read_status /= 0) EXIT                   ! Finished reading

          ! Look for comment
          pos_comment = SCAN(line,"#")                 ! Start position of comment, 0 if none present
          IF (pos_comment > 1) THEN
             line = TRIM(ADJUSTL(line(1:pos_comment))) ! Disregard everything to the right of "#", ..
             ! .. adjust to the left and trim to the right
          ELSE IF (pos_comment == 1) THEN
             CYCLE                                     ! Disregard whole line
          ELSE
             line = TRIM(ADJUSTL(line))                ! Only disregard blanks to the left and right
          ENDIF
          length=LEN_TRIM(line)
          IF(length== 0) CYCLE                 ! Line is empty

          pos = SCAN(line,blank_set)                   ! Position of first blank character
          IF(pos == 0)   CALL finish('init_lctlib',"Wrong syntax in lctlib definitions")

          READ(line(1:pos-1),'(A)') key
          key = tolower(key)

          IF(TRIM(key) == "nlct") CYCLE ! nlct already read above

          SELECT CASE (TRIM(key))
          CASE ('lctnumber')
            READ(line(pos:length),*) lctlib(1:nlct)%LctNumber
            exist_LctNumber = .TRUE.
          CASE ('lctname')                                           ! Name of landcover type
            READ(line(pos:length),*) lctlib(1:nlct)%LctName
            exist_LctName = .TRUE.
          CASE ('landcoverclass')                                    ! Landcover class
            READ(line(pos:length),*) lctlib(1:nlct)%LandcoverClass
            exist_LandcoverClass = .TRUE.
          CASE ('albedosnowvismin')
            READ(line(pos:length),*) lctlib(1:nlct)%AlbedoSnowVisMin
            exist_AlbedoSnowVisMin = .TRUE.
          CASE ('albedosnowvismax')
            READ(line(pos:length),*) lctlib(1:nlct)%AlbedoSnowVisMax
            exist_AlbedoSnowVisMax = .TRUE.
          CASE ('albedosnownirmin')
            READ(line(pos:length),*) lctlib(1:nlct)%AlbedoSnowNirMin
            exist_AlbedoSnowNirMin = .TRUE.
          CASE ('albedosnownirmax')
            READ(line(pos:length),*) lctlib(1:nlct)%AlbedoSnowNirMax
            exist_AlbedoSnowNirMax = .TRUE.
          CASE ('albedosnowmin')
            READ(line(pos:length),*) lctlib(1:nlct)%AlbedoSnowMin
            exist_AlbedoSnowMin = .TRUE.
          CASE ('albedosnowmax')
            READ(line(pos:length),*) lctlib(1:nlct)%AlbedoSnowMax
            exist_AlbedoSnowMax = .TRUE.
          CASE ('albedocanopyvis')
            READ(line(pos:length),*) lctlib(1:nlct)%AlbedoCanopyVIS
            exist_AlbedoCanopyVIS = .TRUE.
          CASE ('albedocanopynir')
            READ(line(pos:length),*) lctlib(1:nlct)%AlbedoCanopyNIR
            exist_AlbedoCanopyNIR = .TRUE.
          CASE ('albedolittervis')
            READ(line(pos:length),*) lctlib(1:nlct)%AlbedoLitterVis
            exist_AlbedoLitterVis = .TRUE.
          CASE ('albedolitternir')
            READ(line(pos:length),*) lctlib(1:nlct)%AlbedoLitterNir
            exist_AlbedoLitterNir = .TRUE.
          CASE ('nitrogenscalingflag')        !Whether nitrogen scaling should be accounted for (.false. input as 0 and .true. as 1)
            READ(line(pos:length),*) itmp(1:nlct)
            lctlib(1:nlct)%nitrogenscalingflag = .FALSE.
            DO i=1,nlct
              IF (itmp(i) /= 0) lctlib(i)%nitrogenscalingflag = .TRUE.
            END DO
            exist_nitrogenScalingFlag  = .TRUE.
          CASE ('c4flag')                              ! Photosynthetic pathway (C4: .true.; C3: .false.)
            READ(line(pos:length),*) itmp(1:nlct)
            lctlib(1:nlct)%c4flag = .FALSE.
            DO i=1,nlct
              IF (itmp(i) == 1) lctlib(i)%c4flag = .TRUE.
            END DO
            exist_C4flag = .TRUE.
          CASE ('carboxrate')                                         ! Carbox rate
            READ(line(pos:length),*) lctlib(1:nlct)%CarboxRate
            exist_CarboxRate = .TRUE.
          CASE ('etransport')                                         ! E-transport
             READ(line(pos:length),*) lctLib(1:nlct)%ETransport
           exist_ETransport = .TRUE.
          CASE ('vegheight')                                          ! typical vegetation height
             READ(line(pos:length),*) lctLib(1:nlct)%VegHeight
             exist_VegHeight = .TRUE.
          CASE ('vegroughness')                                       ! typical vegetation roughness length
             READ(line(pos:length),*) lctLib(1:nlct)%VegRoughness
             exist_VegRoughness = .TRUE.
          CASE ('minvegroughness')                                    ! typical vegetation roughness length at LAI=0
             READ(line(pos:length),*) lctLib(1:nlct)%MinVegRoughness
             exist_MinVegRoughness = .TRUE.
          CASE ('maxvegroughness')                                    ! typical vegetation roughness length at LAI-->inf
             READ(line(pos:length),*) lctLib(1:nlct)%MaxVegRoughness
             exist_MaxVegRoughness = .TRUE.
          CASE ('phenologytype')                                      ! Phenology type
             READ(line(pos:length),*) lctlib(1:nlct)%PhenologyType
             exist_PhenologyType = .TRUE.
          CASE ('canopyresistancemin')                                ! Minimum canopy resistance
             READ(line(pos:length),*) lctlib(1:nlct)%CanopyResistanceMin
             exist_CanopyResistanceMin = .TRUE.
          CASE ('maxlai')                                              ! Maximum LAI for LoGro-P
             READ(line(pos:length),*) lctlib(1:nlct)%MaxLAI
             exist_MaxLAI = .TRUE.
          CASE ('stemarea')                                            ! Area of stems and branches
             READ(line(pos:length),*) lctlib(1:nlct)%StemArea
             exist_StemArea = .TRUE.
          CASE ('specificleafarea_c')                                  ! Carbon content per leaf area for LoGro-P
             READ(line(pos:length),*) lctlib(1:nlct)%specificLeafArea_C
             exist_specificLeafArea_C = .TRUE.
          CASE ('alpha_nr_ind')                                              ! Maximum LAI for LoGro-P
             READ(line(pos:length),*) lctlib(1:nlct)%alpha_nr_ind
             exist_alpha_nr_ind = .TRUE.
          CASE ('beta_nr_ind')                                              ! Maximum LAI for LoGro-P
             READ(line(pos:length),*) lctlib(1:nlct)%beta_nr_ind
             exist_beta_nr_ind = .TRUE.
          CASE ('alpha_leaf')                                              ! Maximum LAI for LoGro-P
             READ(line(pos:length),*) lctlib(1:nlct)%alpha_leaf
             exist_alpha_leaf = .TRUE.
          CASE ('beta_leaf')                                              ! Maximum LAI for LoGro-P
             READ(line(pos:length),*) lctlib(1:nlct)%beta_leaf
             exist_beta_leaf = .TRUE.
          CASE ('knorr_tau_w')                                         ! time before leaf shedding
             READ(line(pos:length),*) lctlib(1:nlct)%knorr_Tau_w
             exist_knorr_tau_w = .TRUE.
          CASE ('knorr_t_phi')                                         ! Temperature trigger for leaf growth
             READ(line(pos:length),*) lctlib(1:nlct)%knorr_T_phi
             exist_knorr_t_phi = .TRUE.
          CASE ('knorr_t_r')                                           ! "Spread" (sigma) of T_phi
             READ(line(pos:length),*) lctlib(1:nlct)%knorr_T_r
             exist_knorr_t_r = .TRUE.
          CASE ('knorr_day_c')                                         ! Day-length at leaf shedding
             READ(line(pos:length),*) lctlib(1:nlct)%knorr_Day_c
             exist_knorr_day_c = .TRUE.
          CASE ('knorr_day_r')                                         ! "Spread" (sigma) of Day_c
             READ(line(pos:length),*) lctlib(1:nlct)%knorr_Day_r
             exist_knorr_day_r = .TRUE.
          CASE ('knorr_k_l')                                           ! Inverse of leaf longevity
             READ(line(pos:length),*) lctlib(1:nlct)%knorr_k_l
             exist_knorr_k_l = .TRUE.
          CASE ('knorr_max_lai')                                       ! maximum LAI
             READ(line(pos:length),*) lctlib(1:nlct)%knorr_max_lai
             exist_knorr_max_lai = .TRUE.
          CASE ('knorr_leaf_growth_rate')                              ! Initial leaf growth rate
             READ(line(pos:length),*) lctlib(1:nlct)%knorr_leaf_growth_rate
             exist_knorr_leaf_growth_rate = .TRUE.
          CASE ('reservec2leafc')                                      ! ratio of C in reserve pool to C in leaves
             READ(line(pos:length),*) lctlib(1:nlct)%reserveC2leafC
             exist_reserveC2leafC = .TRUE.
          CASE ('fract_npp_2_woodpool')                                 ! fraction of NPP directed to wood pool (see cbalance)
             READ(line(pos:length),*) lctlib(1:nlct)%fract_npp_2_woodPool
             exist_fract_npp_2_woodPool = .TRUE.
          CASE ('fract_npp_2_reservepool')                              ! fraction of NPP directed to reserve pool (see cbalance)
             READ(line(pos:length),*) lctlib(1:nlct)%fract_npp_2_reservePool
             exist_fract_npp_2_reservePool = .TRUE.
          CASE ('fract_npp_2_exudates')                                 ! fraction of NPP directed to root exudates (see cbalance)
             READ(line(pos:length),*) lctlib(1:nlct)%fract_npp_2_exudates
             exist_fract_npp_2_exudates = .TRUE.
          CASE ('fract_green_2_herbivory')                              ! fraction of Green directed to herbivory (see cbalance)
             READ(line(pos:length),*) lctlib(1:nlct)%fract_green_2_herbivory
             exist_fract_green_2_herbivory = .TRUE.
          CASE ('tau_c_litter_leaf')                           ! Time constant by which leaf litter pool is depreciated  [days]
             READ(line(pos:length),*) lctlib(1:nlct)%tau_c_litter_leaf
             exist_tau_c_litter_leaf = .TRUE.
          CASE ('tau_c_litter_wood')                           ! Time constant by which woody litter pool is depreciated  [days]
             READ(line(pos:length),*) lctlib(1:nlct)%tau_c_litter_wood
             exist_tau_c_litter_wood = .TRUE.
          CASE ('tau_c_woods')                              ! PFT-specific time_scale for c_woods (and vegetation dynamics)
             READ(line(pos:length),*) lctlib(1:nlct)%tau_c_woods
             ! Conversion from years to days
             lctlib(1:nlct)%tau_c_woods = lctlib(1:nlct)%tau_c_woods * 365._wp
             exist_tau_c_woods = .TRUE.
          CASE ('lai_shed_constant')                               ! Time constant by which leaves are constantly shedded [days-1]
             READ(line(pos:length),*) lctlib(1:nlct)%LAI_shed_constant
             exist_LAI_shed_constant = .TRUE.
          CASE ('fract_c_litter_green2atmos')                    ! fraction of loss of fast pool emitted to atmosphere (see cbalance)
             READ(line(pos:length),*) lctlib(1:nlct)%fract_C_litter_green2atmos
             exist_fract_C_litter_green2atmos = .TRUE.
          CASE ('max_c_content_woods')                          ! Maximum carbon content of wood pool (see cbalance)
             READ(line(pos:length),*) lctlib(1:nlct)%Max_C_content_woods
             exist_Max_C_content_woods = .TRUE.
          CASE ('clumpinessfactor')                             ! Factor to calculate clumpiness of vegetation
             READ(line(pos:length),*) lctlib(1:nlct)%ClumpinessFactor
             exist_ClumpinessFactor = .TRUE.
          CASE ('dynamic_pft')                                  ! Indicates those PFTs that shall take part in vegetation dynamics
            READ(line(pos:length),*) itmp(1:nlct)
            lctlib(1:nlct)%dynamic_PFT = .FALSE.
            DO i=1,nlct
              IF (itmp(i) /= 0) lctlib(i)%dynamic_PFT = .TRUE.
            END DO
            exists_dynamic_PFT = .TRUE.
          CASE ('woody_pft')                                    ! Indicates woody PFTs
            READ(line(pos:length),*) itmp(1:nlct)
            lctlib(1:nlct)%woody_PFT = .FALSE.
            DO i=1,nlct
              IF (itmp(i) /= 0) lctlib(i)%woody_PFT = .TRUE.
            END DO
            exists_woody_PFT = .TRUE.
          CASE ('pasture_pft')                                  ! Indicates pasture PFTs
            READ(line(pos:length),*) itmp(1:nlct)
            lctlib(1:nlct)%pasture_PFT = .FALSE.
            DO i=1,nlct
              IF (itmp(i) /= 0) lctlib(i)%pasture_PFT = .TRUE.
            END DO
            exists_pasture_PFT = .TRUE.
          CASE ('bclimit_min_cold_mmtemp')                      ! PFT-specific minimum coldest monthly mean temperature
             READ(line(pos:length),*) lctlib(1:nlct)%bclimit_min_cold_mmtemp
             exists_bclimit_min_cold_mmtemp = .TRUE.
          CASE ('bclimit_max_cold_mmtemp')                      ! PFT-specific maximum coldest monthly mean temperature
             READ(line(pos:length),*) lctlib(1:nlct)%bclimit_max_cold_mmtemp
             exists_bclimit_max_cold_mmtemp = .TRUE.
          CASE ('bclimit_max_warm_mmtemp')                      ! PFT-specific maximum warmest monthly mean temperature
             READ(line(pos:length),*) lctlib(1:nlct)%bclimit_max_warm_mmtemp
             exists_bclimit_max_warm_mmtemp = .TRUE.
          CASE ('bclimit_min_temprange')                        ! PFT-specific minimum temperature range between warmest and coldest
                                                                !     month (20-year average)
             READ(line(pos:length),*) lctlib(1:nlct)%bclimit_min_temprange
             exists_bclimit_min_temprange = .TRUE.
          CASE ('bclimit_min_gdd')                              ! PFT-specific minimum growing degree days (at or above 5 deg C)
             READ(line(pos:length),*) lctlib(1:nlct)%bclimit_min_gdd
             exists_bclimit_min_gdd = .TRUE.
          CASE ('gdd_base')                                     ! PFT-specific GDD base
             READ(line(pos:length),*) lctlib(1:nlct)%gdd_base
             exists_gdd_base = .TRUE.
          CASE ('upper_tlim')                                   ! PFT-specific upper temperature limit (used for gdd_upper_tlim)
             READ(line(pos:length),*) lctlib(1:nlct)%upper_tlim
             exists_upper_tlim = .TRUE.
          CASE ('fract_wood_2_onsite')                           ! PFT-specific fraction for wood to onSite conversion
             READ(line(pos:length),*) lctlib(1:nlct)%fract_wood_2_onSite
             exists_fract_wood_2_onSite = .TRUE.
          CASE ('fract_wood_2_paper')                            ! PFT-specific fraction for wood to paper conversion
             READ(line(pos:length),*) lctlib(1:nlct)%fract_wood_2_paper
             exists_fract_wood_2_paper = .TRUE.
          CASE ('fract_wood_2_construction')                     ! PFT-specific fraction for wood to construction conversion
             READ(line(pos:length),*) lctlib(1:nlct)%fract_wood_2_construction
             exists_fract_wood_2_construction = .TRUE.
          CASE ('moist_extinction')
            READ(line(pos:length),*) lctlib(1:nlct)%moist_extinction    ! PFT-specific moisture of extinction
             exists_moist_extinction = .TRUE.
          CASE ('fuel_dens')
            READ(line(pos:length),*) lctlib(1:nlct)%fuel_dens           ! PFT-specific fuel density
             exists_fuel_dens = .TRUE.
          CASE ('flame_length_f')
            READ(line(pos:length),*) lctlib(1:nlct)%flame_length_f      ! PFT-specific flame length parameter
             exists_flame_length_f = .TRUE.
          CASE ('crown_length')
            READ(line(pos:length),*) lctlib(1:nlct)%crown_length        ! PFT-specific crown length parameter
             exists_crown_length = .TRUE.
          CASE ('bark_par1')
            READ(line(pos:length),*) lctlib(1:nlct)%bark_par1           ! PFT-specific bark thickness parameter 1
             exists_bark_par1 = .TRUE.
          CASE ('bark_par2')
            READ(line(pos:length),*) lctlib(1:nlct)%bark_par2           ! PFT-specific bark thickness parameter 2
             exists_bark_par2 = .TRUE.
          CASE ('rck')
            READ(line(pos:length),*) lctlib(1:nlct)%RCK                 ! PFT-specific resistance to crown damage
             exists_RCK = .TRUE.
          CASE ('mort_prob')
            READ(line(pos:length),*) lctlib(1:nlct)%mort_prob           ! PFT-specific parameter for mortality
             exists_mort_prob = .TRUE.
          CASE ('litveg_coef')                                  ! PFT-specific fractions to seperate total litter into yasso pools
                                                                        !! (NOT NEEDED!)
             READ(line(pos:length),*) lctlib(1:nlct)%LitVeg_coef(1)
             DO apu = 2, 5
                READ(lctlib_file_unit,'(A256)', IOSTAT=read_status) line
                READ(line(1:length),*) lctlib(1:nlct)%LitVeg_coef(apu)
             END DO
             DO i = 1, nlct
                ! DSG: this check is important, if values don't add up to 1 carbon is not conserved
                IF (SUM(lctlib(i)%LitVeg_coef(:)) .NE. 1.) &
                     CALL finish('init_lctlib','invalid value for LitVeg_coef in '//TRIM(lctlib_file_name))
             END DO
             exists_litveg_coef = .TRUE.
          CASE ('leaflit_coef')                                 ! PFT-specific fractions to seperate leaf litter into yasso pools
             READ(line(pos:length),*) lctlib(1:nlct)%LeafLit_coef(1)
             DO apu = 2, 5
                READ(lctlib_file_unit,'(A256)', IOSTAT=read_status) line
                READ(line(1:length),*) lctlib(1:nlct)%LeafLit_coef(apu)
             END DO
             DO i = 1, nlct
                ! DSG: this check is important, if values dont add up to 1 carbon is not conserved
                IF (SUM(lctlib(i)%LeafLit_coef(:)) .NE. 1.) &
                     CALL finish('init_lctlib','invalid value for LeafLit_coef in '//TRIM(lctlib_file_name))
             END DO
             exists_leaflit_coef = .TRUE.
          CASE ('woodlit_coef')                                 ! PFT-specific  fractions to seperate wood litter into yasso pools
             READ(line(pos:length),*) lctlib(1:nlct)%WoodLit_coef(1)
             DO apu = 2, 5
                READ(lctlib_file_unit,'(A256)', IOSTAT=read_status) line
                READ(line(1:length),*) lctlib(1:nlct)%WoodLit_coef(apu)
             END DO
             DO i = 1, nlct
                ! DSG: this check is important, if values dont add up to 1 carbon is not conserved
                IF (SUM(lctlib(i)%WoodLit_coef(:)) .NE. 1.) &
                     CALL finish('init_lctlib','invalid value for WoodLit_coef_coef in '//TRIM(lctlib_file_name))
             END DO
             exists_woodlit_coef = .TRUE.
          CASE ('woodlittersize')                                ! PFT-specific wood liter size (needed for wood decomposition)
             READ(line(pos:length),*) lctlib(1:nlct)%WoodLitterSize
             exists_woodlittersize_coef = .TRUE.
          CASE default
             ! nothing to do
          END SELECT
       END DO
       DEALLOCATE(itmp)

       IF(.NOT. exist_LctNumber) &
            CALL finish('init_lctlib','No data for LctNumber found in '//TRIM(lctlib_file_name))
       IF(.NOT. exist_LctName) &
!            CALL finish('init_lctlib','No data for LctName found in '//TRIM(lctlib_file_name))
            ! Set default
            lctlib%LctName = "undefined"
       IF(.NOT. exist_LandcoverClass) &
            CALL finish('init_lctlib','No data for LandcoverClass found in '//TRIM(lctlib_file_name))
       IF(.NOT. exist_AlbedoSnowVisMin) &
            CALL finish('init_lctlib','No data for AlbedoSnowVisMin found in '//TRIM(lctlib_file_name))
       IF(.NOT. exist_AlbedoSnowVisMax) &
            CALL finish('init_lctlib','No data for AlbedoSnowVisMax found in '//TRIM(lctlib_file_name))
       IF(.NOT. exist_AlbedoSnowNirMin) &
            CALL finish('init_lctlib','No data for AlbedoSnowNirMin found in '//TRIM(lctlib_file_name))
       IF(.NOT. exist_AlbedoSnowNirMax) &
            CALL finish('init_lctlib','No data for AlbedoSnowNirMax found in '//TRIM(lctlib_file_name))
       IF(.NOT. exist_AlbedoSnowMin) &
            CALL finish('init_lctlib','No data for AlbedoSnowMin found in '//TRIM(lctlib_file_name))
       IF(.NOT. exist_AlbedoSnowMax) &
            CALL finish('init_lctlib','No data for AlbedoSnowMax found in '//TRIM(lctlib_file_name))
       IF(.NOT. exist_AlbedoCanopyVIS) &
            CALL finish('init_lctlib','No data for AlbedoCanopyVIS found in '//TRIM(lctlib_file_name))
       IF(.NOT. exist_AlbedoCanopyNIR) &
            CALL finish('init_lctlib','No data for albedoCanopyNIR found in '//TRIM(lctlib_file_name))
       IF(.NOT. exist_AlbedoLitterVis) &
            CALL finish('init_lctlib','No data for AlbedoLitterVis found in '//TRIM(lctlib_file_name))
       IF(.NOT. exist_AlbedoLitterNir) &
            CALL finish('init_lctlib','No data for AlbedoLitterNir found in '//TRIM(lctlib_file_name))
       IF(.NOT. exist_nitrogenScalingFlag) &
            CALL finish('init_lctLib','No data for NitrogenScalingFlag found in '//TRIM(lctLib_file_name))
       IF(.NOT. exist_C4flag) &
            CALL finish('init_lctLib','No data for C4flag found in '//TRIM(lctLib_file_name))
       IF(.NOT. exist_CarboxRate) &
            CALL finish('init_lctLib','No data for CarboxRate found in '//TRIM(lctLib_file_name))
       IF(.NOT. exist_ETransport) &
            CALL finish('init_lctLib','No data for ETransport found in '//TRIM(lctLib_file_name))
       IF(.NOT. exist_VegHeight) &
            CALL finish('init_lctLib','No data for VegHeight found in '//TRIM(lctLib_file_name))
       IF(.NOT. exist_VegRoughness) &
            CALL finish('init_lctLib','No data for VegRoughness found in '//TRIM(lctLib_file_name))
       IF(.NOT. exist_MinVegRoughness) &
            CALL finish('init_lctLib','No data for MinVegRoughness found in '//TRIM(lctLib_file_name))
       IF(.NOT. exist_MaxVegRoughness) &
            CALL finish('init_lctLib','No data for MaxVegRoughness found in '//TRIM(lctLib_file_name))
       IF(.NOT. exist_PhenologyType) &
            CALL finish('init_lctlib','No data for PhenologyType found in '//TRIM(lctlib_file_name))
       IF (.NOT. exist_CanopyResistanceMin) &
            lctlib%CanopyResistanceMin = -1.0_wp
       IF(.NOT. exist_MaxLAI) &
            CALL finish('init_lctlib','No data for MaxLAI found in '//TRIM(lctlib_file_name))
       IF(.NOT. exist_StemArea) &
            CALL finish('init_lctlib','No data for StemArea found in '//TRIM(lctlib_file_name))
       IF(.NOT. exist_specificLeafArea_C) &
            CALL finish('init_lctlib','No data for specificLeafArea_C found in '//TRIM(lctlib_file_name))

       IF(.NOT. exist_alpha_nr_ind) &
            CALL finish('init_lctlib','No data for alpha_nr_ind found in '//TRIM(lctlib_file_name))
       IF(.NOT. exist_beta_nr_ind) &
            CALL finish('init_lctlib','No data for beta_nr_ind found in '//TRIM(lctlib_file_name))
       IF(.NOT. exist_alpha_leaf) &
            CALL finish('init_lctlib','No data for alpha_leaf found in '//TRIM(lctlib_file_name))
       IF(.NOT. exist_beta_leaf) &
            CALL finish('init_lctlib','No data for beta_leaf found in '//TRIM(lctlib_file_name))

      IF(.NOT. exist_knorr_tau_w) &
        CALL finish('init_lctlib','No data for knorr_Tau_w found in '//TRIM(lctlib_file_name))
      IF(.NOT. exist_knorr_t_phi) &
        CALL finish('init_lctlib','No data for knorr_T_phi found in '//TRIM(lctlib_file_name))
      IF(.NOT. exist_knorr_t_r) &
        CALL finish('init_lctlib','No data for knorr_T_rfound in '//TRIM(lctlib_file_name))
      IF(.NOT. exist_knorr_day_c) &
        CALL finish('init_lctlib','No data for knorr_Day_c found in '//TRIM(lctlib_file_name))
      IF(.NOT. exist_knorr_day_r) &
        CALL finish('init_lctlib','No data for knorr_Day_r found in '//TRIM(lctlib_file_name))
      IF(.NOT. exist_knorr_k_l) &
        CALL finish('init_lctlib','No data for knorr_k_l found in '//TRIM(lctlib_file_name))
      IF(.NOT. exist_knorr_max_lai) &
        CALL finish('init_lctlib','No data for knorr_max_lai found in '//TRIM(lctlib_file_name))
      IF(.NOT. exist_knorr_leaf_growth_rate) &
        CALL finish('init_lctlib','No data for leaf_growth_rate found in '//TRIM(lctlib_file_name))

       IF(.NOT. exist_reserveC2leafC) &
            CALL finish('init_lctlib','No data for reserveC2leafC found in '//TRIM(lctlib_file_name))
       IF(.NOT. exist_fract_npp_2_woodPool) &
            CALL finish('init_lctlib','No data for fract_npp_2_woodPool found in '//TRIM(lctlib_file_name))
       IF(.NOT. exist_fract_npp_2_reservePool) &
            CALL finish('init_lctlib','No data for fract_npp_2_reservePool found in '//TRIM(lctlib_file_name))
       IF(.NOT. exist_fract_npp_2_exudates) &
            CALL finish('init_lctlib','No data for fract_npp_2_exudates found in '//TRIM(lctlib_file_name))
       IF(.NOT. exist_fract_green_2_herbivory) &
            CALL finish('init_lctlib','No data for fract_green_2_herbivory found in '//TRIM(lctlib_file_name))
       IF(.NOT. exist_tau_c_litter_leaf) &
            CALL finish('init_lctlib','No data for tau_c_litter_leaf found in '//TRIM(lctlib_file_name))
       IF(.NOT. exist_tau_c_litter_wood) &
            CALL finish('init_lctlib','No data for tau_c_litter_wood found in '//TRIM(lctlib_file_name))
       IF(.NOT. exist_tau_c_woods) &
            CALL finish('init_lctlib','No data for tau_c_woods found in '//TRIM(lctlib_file_name))
       IF(.NOT. exist_LAI_shed_constant) &
            CALL finish('init_lctlib','No data for LAI_shed_constant found in '//TRIM(lctlib_file_name))
       IF(.NOT. exist_fract_C_litter_green2atmos) &
            CALL finish('init_lctlib','No data for fract_C_litter_green2atmos found in '//TRIM(lctlib_file_name))
       IF(.NOT. exist_Max_C_content_woods) &
            CALL finish('init_lctlib','No data for Max_C_content_woods found in '//TRIM(lctlib_file_name))
       IF(.NOT. exist_ClumpinessFactor) &
            CALL finish('init_lctlib','No data for ClupinessFactor found in '//TRIM(lctlib_file_name))
       IF(.NOT. exists_dynamic_PFT) &
            CALL finish('init_lctlib','No data for DYNAMIC_PFT found in '//TRIM(lctlib_file_name))
       IF(.NOT. exists_woody_PFT) &
            CALL finish('init_lctlib','No data for WOODY_PFT found in '//TRIM(lctlib_file_name))
       IF(.NOT. exists_pasture_PFT) &
            CALL finish('init_lctlib','No data for PASTURE_PFT found in '//TRIM(lctlib_file_name))
       IF(.NOT. exists_bclimit_min_cold_mmtemp) &
            CALL finish('init_lctlib','No data for BCLIMIT_MIN_COLD_MMTEMP found in '//TRIM(lctlib_file_name))
       IF(.NOT. exists_bclimit_max_cold_mmtemp) &
            CALL finish('init_lctlib','No data for BCLIMIT_MAX_COLD_MMTEMP found in '//TRIM(lctlib_file_name))
       IF(.NOT. exists_bclimit_max_warm_mmtemp) &
            CALL finish('init_lctlib','No data for BCLIMIT_MAX_WARM_MMTEMP found in '//TRIM(lctlib_file_name))
       IF(.NOT. exists_bclimit_min_temprange ) &
            CALL finish('init_lctlib','No data for BCLIMIT_MIN_TEMPRANGE found in '//TRIM(lctlib_file_name))
       IF(.NOT.exists_bclimit_min_gdd ) &
            CALL finish('init_lctlib','No data for BCLIMIT_MIN_GDD found in '//TRIM(lctlib_file_name))
       IF(.NOT. exists_gdd_base) &
            CALL finish('init_lctlib','No data for GDD_BASE found in '//TRIM(lctlib_file_name))
       IF(.NOT. exists_upper_tlim) &
            CALL finish('init_lctlib','No data for UPPER_TLIM found in '//TRIM(lctlib_file_name))
       IF(.NOT. exists_fract_wood_2_onSite) &
            CALL finish('init_lctlib','No data for FRACT_WOOD_2_ONSITE found in '//TRIM(lctlib_file_name))
       IF(.NOT. exists_fract_wood_2_paper) &
            CALL finish('init_lctlib','No data for FRACT_WOOD_2_PAPER found in '//TRIM(lctlib_file_name))
       IF(.NOT. exists_fract_wood_2_construction) &
            CALL finish('init_lctlib','No data for FRACT_WOOD_2_CONSTRUCTION found in '//TRIM(lctlib_file_name))
       IF(.NOT. exists_moist_extinction) &
            CALL finish('init_lctlib','No data for moisture of extinction found in '//TRIM(lctlib_file_name))
       IF(.NOT. exists_fuel_dens) &
            CALL finish('init_lctlib','No data for bulk fuel density found in '//TRIM(lctlib_file_name))
       IF(.NOT. exists_flame_length_f) &
            CALL finish('init_lctlib','No data for flame length f found in '//TRIM(lctlib_file_name))
       IF(.NOT. exists_crown_length) &
            CALL finish('init_lctlib','No data for crown length found in '//TRIM(lctlib_file_name))
       IF(.NOT. exists_bark_par1) &
            CALL finish('init_lctlib','No data for bark_par1 found in '//TRIM(lctlib_file_name))
       IF(.NOT. exists_bark_par2) &
            CALL finish('init_lctlib','No data for bark_par2 found in '//TRIM(lctlib_file_name))
       IF(.NOT. exists_RCK) &
            CALL finish('init_lctlib','No data for RCK (resistance to crown damage) found in '//TRIM(lctlib_file_name))
       IF(.NOT. exists_mort_prob) &
            CALL finish('init_lctlib','No data for mortality probability found in '//TRIM(lctlib_file_name))
       IF(.NOT. exists_woodlittersize_coef) &
            CALL finish('init_lctlib','No data for WOODLITTERSIZE found in '//TRIM(lctlib_file_name))
       IF(.NOT. exists_woodlit_coef) &
            CALL finish('init_lctlib','No data for WOODLIT_COEF found in '//TRIM(lctlib_file_name))
       IF(.NOT. exists_leaflit_coef) &
            CALL finish('init_lctlib','No data for LEAFLIT_COEF found in '//TRIM(lctlib_file_name))
       IF(.NOT. exists_litveg_coef) &
            CALL finish('init_lctlib','No data for LITVEG_COEF found in '//TRIM(lctlib_file_name))

       CLOSE(unit=lctlib_file_unit) !! Closing access to lctlib file

       !------------------------------------------------------------------------------------------
       ! kalle, 151004
       ! Set carbox rate to  [Mol(CO2)/m^2/s]
       ! SET e-transport rate to [Mol(CO2)/m^2/s]
       DO i=1,nlct
         lctLib(i)%ETransport=lctLib(i)%ETransport * 1.e-06_wp
         lctlib(i)%CarboxRate=lctlib(i)%CarboxRate * 1.e-06_wp
       END DO
       !------------------------------------------------------------------------------------------

       !! --- translate landcover class information into landcover class flags

       lctlib(1:nlct)%BareSoilFlag   = .FALSE.
       lctlib(1:nlct)%GlacierFlag    = .FALSE.
       lctlib(1:nlct)%LakeFlag       = .FALSE.
       lctlib(1:nlct)%ForestFlag     = .FALSE.
       lctlib(1:nlct)%GrassFlag      = .FALSE.
       lctlib(1:nlct)%CropFlag       = .FALSE.
       lctlib(1:nlct)%PastureFlag    = .FALSE.
       lctlib(1:nlct)%NaturalVegFlag = .FALSE.

       DO i=1,nlct
          SELECT CASE (lctlib(i)%LandcoverClass)
          CASE (0)
             lctlib(i)%BareSoilFlag = .TRUE.
          CASE (1)
             lctlib(i)%GlacierFlag = .TRUE.
          CASE (2)
             lctlib(i)%LakeFlag = .TRUE.
          CASE (3)
             lctlib(i)%ForestFlag = .TRUE.
             lctlib(i)%NaturalVegFlag = .TRUE.
          CASE (4)
             lctlib(i)%GrassFlag = .TRUE.
             lctlib(i)%NaturalVegFlag = .TRUE.
          CASE (5)
             lctlib(i)%NaturalVegFlag = .TRUE.
          CASE (6)
             lctlib(i)%CropFlag = .TRUE.
          CASE (7)
             lctlib(i)%PastureFlag = .TRUE.
          CASE default
             WRITE(message_text,*) 'LandcoverClass ', lctlib(i)%LandcoverClass,' not allowed in ', TRIM(lctlib_file_name)
             CALL finish('init_lctlib', message_text)
          END SELECT
       END DO

       npft = 0
       DO i=1,nlct
          IF(lctlib(i)%naturalVegFlag .OR. lctlib(i)%cropFlag .OR. lctlib(i)%PastureFlag) THEN

             npft = npft + 1 ! count PFTs


             ! Check consistency with PhenologyTypes
             IF(lctlib(i)%naturalVegFlag .AND. &
                  (lctlib(i)%PhenologyType .LE.0 .OR. lctlib(i)%PhenologyType > 4) ) THEN
                CALL finish('init_lctlib',&
                     'Inconsistent entries for one LandcoverClass and PhenologyType in '//TRIM(lctlib_file_name))
             END IF
             IF(lctlib(i)%CropFlag .AND. .NOT. lctlib(i)%PhenologyType==5 ) THEN
                CALL finish('init_lctlib',&
                     'Inconsistent entries for LandcoverClass "crops" and PhenologyType in '//TRIM(lctlib_file_name))
             END IF

          END IF

          !! --- check consistency between landcover classes and specifications for the dynamic vegetation

          IF (lctlib(i)%dynamic_pft  .AND. .NOT. lctlib(i)%NaturalVegFlag) THEN
             CALL finish('init_lctlib',&
                     'Inconsistent entries: Every DYNAMIC_PFT must belong to LandcoverClass "natural" in '//TRIM(lctlib_file_name))
          END IF

          IF(lctlib(i)%woody_pft .AND. .NOT. lctlib(i)%NaturalVegFlag) THEN
             CALL finish('init_lctlib',&
                     'Inconsistent entries: Every WOODY_PFT must belong to LandcoverClass "natural" in '//TRIM(lctlib_file_name))
          END IF

          IF( lctlib(i)%GrassFlag .NEQV. (lctlib(i)%dynamic_pft .AND. .NOT. lctlib(i)%woody_PFT )) THEN
             CALL finish('init_lctlib',&
                  'Combination of DYNAMIC_PFT and WOODY_PFT is not consistent with LandcoverClass "grasses" in '//&
                                                                                                     TRIM(lctlib_file_name))
          END IF
       END DO

    END IF !p_parallel_io


    IF(my_process_is_mpi_parallel()) THEN
       CALL p_bcast(npft, p_io, mpi_comm)
       CALL p_bcast(lctlib(:)%LctNumber, p_io, mpi_comm)
       DO i=1,nlct
          CALL p_bcast(lctlib(i)%LctName, p_io, mpi_comm)
       ENDDO
       CALL p_bcast(lctlib(:)%LandcoverClass, p_io, mpi_comm)
       CALL p_bcast(lctlib(:)%NaturalVegFlag,p_io, mpi_comm)
       CALL p_bcast(lctlib(:)%ForestFlag,p_io, mpi_comm)
       CALL p_bcast(lctlib(:)%GrassFlag,p_io, mpi_comm)
       CALL p_bcast(lctlib(:)%CropFlag,p_io, mpi_comm)
       CALL p_bcast(lctlib(:)%PastureFlag,p_io, mpi_comm)
       CALL p_bcast(lctlib(:)%LakeFlag, p_io, mpi_comm)
       CALL p_bcast(lctlib(:)%GlacierFlag, p_io, mpi_comm)
       CALL p_bcast(lctlib(:)%BareSoilFlag, p_io, mpi_comm)
       CALL p_bcast(lctlib(:)%AlbedoSnowVisMin, p_io, mpi_comm)
       CALL p_bcast(lctlib(:)%AlbedoSnowVisMax, p_io, mpi_comm)
       CALL p_bcast(lctlib(:)%AlbedoSnowNirMin, p_io, mpi_comm)
       CALL p_bcast(lctlib(:)%AlbedoSnowNirMax, p_io, mpi_comm)
       CALL p_bcast(lctlib(:)%AlbedoSnowMin, p_io, mpi_comm)
       CALL p_bcast(lctlib(:)%AlbedoSnowMax, p_io, mpi_comm)
       CALL p_bcast(lctlib(:)%AlbedoCanopyVIS, p_io, mpi_comm)
       CALL p_bcast(lctlib(:)%AlbedoCanopyNIR, p_io, mpi_comm)
       CALL p_bcast(lctlib(:)%AlbedoLitterVis, p_io, mpi_comm)
       CALL p_bcast(lctlib(:)%AlbedoLitterNir, p_io, mpi_comm)
       CALL p_bcast(lctlib(:)%NitrogenScalingFlag,p_io, mpi_comm)
       CALL p_bcast(lctlib(:)%C4flag, p_io, mpi_comm)
       CALL p_bcast(lctlib(:)%CarboxRate, p_io, mpi_comm)
       CALL p_bcast(lctlib(:)%ETransport, p_io, mpi_comm)
       CALL p_bcast(lctlib(:)%VegHeight, p_io, mpi_comm)
       CALL p_bcast(lctlib(:)%VegRoughness, p_io, mpi_comm)
       CALL p_bcast(lctlib(:)%MinVegRoughness, p_io, mpi_comm)
       CALL p_bcast(lctlib(:)%MaxVegRoughness, p_io, mpi_comm)
       CALL p_bcast(lctlib(:)%PhenologyType, p_io, mpi_comm)
       CALL p_bcast(lctlib(:)%CanopyResistanceMin, p_io, mpi_comm)
       CALL p_bcast(lctlib(:)%MaxLAI, p_io, mpi_comm)
       CALL p_bcast(lctlib(:)%StemArea, p_io, mpi_comm)
       CALL p_bcast(lctlib(:)%specificLeafArea_C, p_io, mpi_comm)
       CALL p_bcast(lctlib(:)%alpha_nr_ind, p_io, mpi_comm)
       CALL p_bcast(lctlib(:)%beta_nr_ind, p_io, mpi_comm)
       CALL p_bcast(lctlib(:)%alpha_leaf, p_io, mpi_comm)
       CALL p_bcast(lctlib(:)%beta_leaf, p_io, mpi_comm)
       CALL p_bcast(lctlib(:)%knorr_Tau_w, p_io, mpi_comm)
       CALL p_bcast(lctlib(:)%knorr_T_phi, p_io, mpi_comm)
       CALL p_bcast(lctlib(:)%knorr_T_r, p_io, mpi_comm)
       CALL p_bcast(lctlib(:)%knorr_Day_c, p_io, mpi_comm)
       CALL p_bcast(lctlib(:)%knorr_Day_r, p_io, mpi_comm)
       CALL p_bcast(lctlib(:)%knorr_k_l, p_io, mpi_comm)
       CALL p_bcast(lctlib(:)%knorr_leaf_growth_rate, p_io, mpi_comm)
       CALL p_bcast(lctlib(:)%knorr_max_lai, p_io, mpi_comm)
       CALL p_bcast(lctlib(:)%reserveC2leafC, p_io, mpi_comm)
       CALL p_bcast(lctlib(:)%fract_npp_2_woodPool, p_io, mpi_comm)
       CALL p_bcast(lctlib(:)%fract_npp_2_reservePool, p_io, mpi_comm)
       CALL p_bcast(lctlib(:)%fract_npp_2_exudates, p_io, mpi_comm)
       CALL p_bcast(lctlib(:)%fract_green_2_herbivory, p_io, mpi_comm)
       CALL p_bcast(lctlib(:)%tau_c_litter_leaf, p_io, mpi_comm)
       CALL p_bcast(lctlib(:)%tau_c_litter_wood, p_io, mpi_comm)
       CALL p_bcast(lctlib(:)%tau_c_woods, p_io, mpi_comm)
       CALL p_bcast(lctlib(:)%LAI_shed_constant, p_io, mpi_comm)
       CALL p_bcast(lctlib(:)%fract_C_litter_green2atmos, p_io, mpi_comm)
       CALL p_bcast(lctlib(:)%Max_C_content_woods, p_io, mpi_comm)
       CALL p_bcast(lctlib(:)%ClumpinessFactor, p_io, mpi_comm)
       CALL p_bcast(lctlib(:)%dynamic_pft, p_io, mpi_comm)
       CALL p_bcast(lctlib(:)%woody_pft, p_io, mpi_comm)
       CALL p_bcast(lctlib(:)%pasture_pft, p_io, mpi_comm)
       CALL p_bcast(lctlib(:)%bclimit_min_cold_mmtemp, p_io, mpi_comm)
       CALL p_bcast(lctlib(:)%bclimit_max_cold_mmtemp, p_io, mpi_comm)
       CALL p_bcast(lctlib(:)%bclimit_max_warm_mmtemp, p_io, mpi_comm)
       CALL p_bcast(lctlib(:)%bclimit_min_temprange, p_io, mpi_comm)
       CALL p_bcast(lctlib(:)%bclimit_min_gdd, p_io, mpi_comm)
       CALL p_bcast(lctlib(:)%gdd_base, p_io, mpi_comm)
       CALL p_bcast(lctlib(:)%upper_tlim, p_io, mpi_comm)
       CALL p_bcast(lctlib(:)%fract_wood_2_onSite, p_io, mpi_comm)
       CALL p_bcast(lctlib(:)%fract_wood_2_paper, p_io, mpi_comm)
       CALL p_bcast(lctlib(:)%fract_wood_2_construction, p_io, mpi_comm)
       CALL p_bcast(lctlib(:)%moist_extinction, p_io, mpi_comm)
       CALL p_bcast(lctlib(:)%fuel_dens, p_io, mpi_comm)
       CALL p_bcast(lctlib(:)%flame_length_f, p_io, mpi_comm)
       CALL p_bcast(lctlib(:)%crown_length, p_io, mpi_comm)
       CALL p_bcast(lctlib(:)%bark_par1, p_io, mpi_comm)
       CALL p_bcast(lctlib(:)%bark_par2, p_io, mpi_comm)
       CALL p_bcast(lctlib(:)%RCK, p_io, mpi_comm)
       CALL p_bcast(lctlib(:)%mort_prob, p_io, mpi_comm)
       DO i=1,5
         CALL p_bcast(lctlib(:)%LitVeg_coef(i), p_io, mpi_comm)
         CALL p_bcast(lctlib(:)%LeafLit_coef(i), p_io, mpi_comm)
         CALL p_bcast(lctlib(:)%WoodLit_coef(i), p_io, mpi_comm)
       END DO
       CALL p_bcast(lctlib(:)%WoodLitterSize, p_io, mpi_comm)
    END IF

    IF(npft>nlct) CALL finish('init_lctlib','npft larger than nlct: '//int2string(npft)//' > '//int2string(nlct))

    ! Throw away values for glacier (first column)
    nlct = nlct - 1
    ALLOCATE(return_value(nlct))
    return_value(1:nlct) = lctlib(2:nlct+1)

  END FUNCTION Read_lctlib

END MODULE mo_jsb_lctlib_class
