!>
!! @brief Contains methods for initialization of JSBACH usecases.
!!
!! @author
!!  Reiner Schnur, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Reiner Schnur (2013-04-04)
!!
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/license .
!!

MODULE mo_jsb_model_usecases
#ifndef __NO_JSBACH__

  USE mo_exception,           ONLY: message, finish

  USE mo_jsb_model_class,     ONLY: t_jsb_model
  USE mo_jsb_process_class,   ONLY: A2L_, L2A_, SEB_, TURB_, SSE_, HYDRO_, RAD_, HD_, ASSIMI_, PHENO_, CARBON_, DISTURB_, &
                                    FUEL_, JEDI_, ALCC_, VEG_
  USE mo_jsb_process_class,   ONLY: ON_LEAFS_, ON_TILE_, AGGREGATE_ !, ON_SUBTREE_
  USE mo_jsb_tile,            ONLY: t_jsb_tile
  USE mo_jsb_lct_class,       ONLY: t_jsb_lct, &
    &                               LAND_TYPE, VEG_TYPE, GLACIER_TYPE, LAKE_TYPE !, BARE_TYPE

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: init_usecase

  CHARACTER(len=*), PARAMETER :: modname = 'mo_jsb_model_usecases'

CONTAINS

  SUBROUTINE init_usecase(model)

    TYPE(t_jsb_model), POINTER, INTENT(inout) :: model

    CHARACTER(len=*), PARAMETER :: routine = modname//':init_usecase'

    IF (model%config%use_lakes) THEN
      CALL message(TRIM(routine), 'Initializing JSBACH usecase "'//TRIM(model%config%usecase)//'" with lakes')
    ELSE
      CALL message(TRIM(routine), 'Initializing JSBACH usecase "'//TRIM(model%config%usecase)//'"')
    END IF


    SELECT CASE (TRIM(model%config%usecase))

    CASE ('jsbach_lite')
      CALL init_usecase_lite(model)

    CASE ('jsbach_pfts')
      CALL init_usecase_pfts(model)

    CASE ('jsbach_jedi')           !PHS added
      CALL init_usecase_jedi(model)  !PHS added

    CASE DEFAULT
      CALL finish(TRIM(routine), 'JSBACH usecase not defined.')

    END SELECT

  END SUBROUTINE init_usecase

  SUBROUTINE init_usecase_lite(model)

    TYPE(t_jsb_model), POINTER, INTENT(inout) :: model

    CLASS(t_jsb_tile), POINTER :: box_tile, lake_tile, land_tile, glacier_tile, veg_tile

    CHARACTER(len=*), PARAMETER :: routine = modname//':init_usecase_lite'

    ! Create root of the tile structure representing the whole grid box
    box_tile => t_jsb_tile('box', 'Top tile', &
      & processes      =(/A2L_,     L2A_,     RAD_,      HYDRO_,    TURB_    /), &
      & process_actions=(/ON_TILE_, ON_TILE_, ON_LEAFS_, ON_LEAFS_, ON_LEAFS_/), &
      & model_id=model%id, fract_varname='notsea')
    model%top => box_tile

#ifndef __NO_JSBACH_HD__
    IF (model%processes(HD_)%p%config%active) THEN
      CALL box_tile%Set_process_action(HD_, ON_TILE_)
    END IF
#endif

      ! Create sub-tiles
      IF (model%config%use_lakes) THEN
        lake_tile => t_jsb_tile('lake', 'Lake tile', &
          & processes=(/SEB_/), process_actions=(/ON_TILE_/), &
          & lct=t_jsb_lct(LAKE_TYPE), parent=box_tile, fract_varname='fract_lake')
        land_tile => t_jsb_tile('land', 'Land tile', &
          & processes=(/SEB_/), process_actions=(/ON_TILE_/), &
          & lct=t_jsb_lct(LAND_TYPE), parent=box_tile)
      ELSE
        land_tile => t_jsb_tile('land', 'Land tile', &
          & processes=(/SEB_/), process_actions=(/ON_TILE_/), &
          & lct=t_jsb_lct(LAND_TYPE), &
          & parent=box_tile, &
          & fract_varname='equal')   ! Fraction = 1.
      END IF

        IF (model%config%use_glacier) THEN
          glacier_tile => t_jsb_tile('glac', 'Glacier tile', processes=(/SSE_/), &
            & lct=t_jsb_lct(GLACIER_TYPE), parent=land_tile, &
            & fract_varname='fract_glac')
        END IF

          veg_tile => t_jsb_tile('veg', 'Veg tile', processes=(/SSE_, PHENO_/), &
            & lct=t_jsb_lct(VEG_TYPE), parent=land_tile,  &
            & fract_varname='fract_veg')

    ! Next 3 lines only to avoid compiler warnings about pointers not being dereferenced
    IF (ASSOCIATED(glacier_tile)) CONTINUE
    IF (ASSOCIATED(lake_tile)) CONTINUE
    IF (ASSOCIATED(veg_Tile)) CONTINUE

    NULLIFY(box_tile, lake_tile, glacier_tile, land_tile, veg_tile)

  END SUBROUTINE init_usecase_lite

  SUBROUTINE init_usecase_pfts(model)

    TYPE(t_jsb_model), POINTER, INTENT(inout) :: model

    CLASS(t_jsb_tile), POINTER :: box_tile, lake_tile, land_tile, glacier_tile, veg_tile, pft_tile

    INTEGER :: i

    INTEGER, PARAMETER :: npft = 11

    CHARACTER(len=5), PARAMETER :: pft_shortnames(npft) = [character(len=5) ::  &
      & 'pft01', 'pft02', 'pft03', 'pft04', 'pft05', 'pft06', &
      & 'pft07', 'pft08', 'pft09', 'pft10', 'pft11' ]

    CHARACTER(len=30), PARAMETER :: pft_longnames(npft) = [character(len=30) :: &
      & 'Tropical broadleaf evergreen', &
      & 'Tropical broadleaf deciduous', &
      & 'Extra-tropical evergreen',     &
      & 'Extra-tropical deciduous',     &
      & 'Raingreen shrubs',             &
      & 'Deciduous shrubs',             &
      & 'C3 grass',                     &
      & 'C4 grass',                     &
      & 'C3 pasture',                   &
      & 'C4 pasture',                   &
      & 'C3 crops'                    ]

    ! Index of PFT into lctlib columns
    ! id is one less than in lctlib file because the first column (glacier) has been thrown away
    INTEGER, PARAMETER :: pft_ids(npft) = (/ 1,2,3,4,9,10,11,12,14,15,19 /)

    CHARACTER(len=*), PARAMETER :: routine = modname//':init_usecase_pfts'


    ! Create root of the tile structure representing the whole grid box
    box_tile => t_jsb_tile('box', 'Top tile', &
      & processes      =(/A2L_,     L2A_,     RAD_  ,    TURB_    /), &
      & process_actions=(/ON_TILE_, ON_TILE_, ON_LEAFS_, ON_LEAFS_/), &
      & model_id=model%id, fract_varname='notsea')
    model%top => box_tile

#ifndef __NO_JSBACH_HD__
    IF (model%processes(HD_)%p%config%active) THEN
      CALL box_tile%Set_process_action(HD_, ON_TILE_)
    END IF
#endif

      IF (model%config%use_lakes) THEN
        lake_tile => t_jsb_tile('lake', 'Lake tile', &
          & processes=(/SEB_, HYDRO_/), process_actions=(/ON_TILE_, ON_TILE_/), &
          & lct=t_jsb_lct(LAKE_TYPE), parent=box_tile)
        land_tile => t_jsb_tile('land', 'Land tile', &
          & processes=(/SEB_/), process_actions=(/ON_TILE_/), &
          & lct=t_jsb_lct(LAND_TYPE), parent=box_tile)
      ELSE
        land_tile => t_jsb_tile('land', 'Land tile', &
          & processes=(/SEB_/), process_actions=(/ON_TILE_/), &
          & lct=t_jsb_lct(LAND_TYPE), parent=box_tile, &
          & fract_varname='equal')   ! Fraction = 1.
      END IF

        IF (model%config%use_glacier) THEN
          glacier_tile => t_jsb_tile('glac', 'Glacier tile', &
            & processes = (/HYDRO_, SSE_/), &
            & process_actions = (/ON_TILE_, ON_TILE_/), &
            & lct=t_jsb_lct(GLACIER_TYPE), parent=land_tile)
            ! & lct=t_jsb_lct(GLACIER_TYPE, lib_id=1), parent=land_tile)
        END IF

        ! default jsb4: HYDRO_ is active at veg tile only to save runtime
        IF (.NOT. model%config%use_quincy) THEN
          veg_tile => t_jsb_tile('veg', 'Veg tile', &
            & processes       = (/HYDRO_,   SSE_,     ASSIMI_,   PHENO_/),    &
            & process_actions = (/ON_TILE_, ON_TILE_, ON_LEAFS_, ON_LEAFS_/), &
            & lct=t_jsb_lct(VEG_TYPE), parent=land_tile,  &
            & fract_varname='fract_veg')
          IF (model%processes(ALCC_)%p%config%active) THEN
            CALL veg_tile%Set_process_action(ALCC_, ON_TILE_, AGGREGATE_)
          END IF
          !
          ! use_quincy: 
        !  HYDRO_ is active at PFT tiles to make all the HYDRO_ variables available to the PFT tiles
        !    note: HYDRO_ results differ only marginally between ON_TILE_ & ON_LEAFS_ (tested by Reiner)
        !  VEG_ is active
        ELSE
          veg_tile => t_jsb_tile('veg', 'Veg tile', &
            & processes       = (/HYDRO_,    SSE_,     ASSIMI_,   PHENO_,    VEG_   /), &
            & process_actions = (/ON_LEAFS_, ON_TILE_, ON_LEAFS_, ON_LEAFS_, ON_LEAFS_/), &
            & lct=t_jsb_lct(VEG_TYPE), parent=land_tile,  &
            & fract_varname='fract_veg')
        ENDIF

          ! R: TBD: DISTURB_ ist kein Task von CARBON_,
          !         da es die cover fractions der pfts verändert und nur als Nebeneffekt auch die C-Pools
          !         verändert. Er würde im Prinzip auch ablaufen ohne den Prozess CARBON_ und würde auch
          !         auf alle anderen Stoffpools wirken. Außerdem soll CARBON_ später wie NITROGEN_
          !         oder andere Stoffpools gehandelt werden. Diese Prozesse haben bis
          !         auf ihre Fixierung kaum eigene Tasks sondern werden von anderen Prozessen verändert.
          !         => die Stoffpools sollten langfristig evtl. gar keine Prozesse sein
          !
          !         DISTURB_ braucht CARBON,
          !         da es C-Pools als Input für die Feuerberechnung braucht und sollte diese auch verändern
          !         wenn der CARBON Prozess existiert. Dies ist ein grundsätzliches Problem: Prozesse benötigen
          !         - oder sogar verändern - evtl. Variablen, die eigentlich zu einem anderen Prozess gehören.
          !         Das bedeutet der Speicher des anderen Prozesses muß angelegt werden auch wenn dieser
          !         Prozess gar nicht abläuft. Dazu müsste wohl am besten eine Prozedur die die Namelist-
          !         Einträge entwprechend umändert geschrieben werden oder wenn Usecases benutzt werden
          !         (wie hier) der entsprechende Prozess dennoch definiert werden, auch wenn in der Queue-
          !         Liste nachher kein Task dieses Prozesses auftaucht.
          !
          !         DISTURB_ braucht PHENO wg. veg_fract_correct s.o.
          !
          !         CARBON_ braucht DISTURB_ nicht!
          !
          !         Manche Prozeduren sollen nur ausgeführt werden wenn beide Prozesse aktiv sind.
          !         In DISTURB_ z.B. die dazugehörige C-Umverteilung bei Feuer und Winbruch.
          !         Dies habe ich mit "IF ( tile%Is_process_active(CARBON_) ) THEN" gemacht.
          !
          IF (model%processes(CARBON_)%p%config%active) THEN
            CALL veg_tile%Set_process_action(CARBON_, ON_LEAFS_, AGGREGATE_)
          END IF
          IF (model%processes(FUEL_)%p%config%active) THEN
            CALL veg_tile%Set_process_action(FUEL_, ON_TILE_, AGGREGATE_)
          END IF
          IF (model%processes(DISTURB_)%p%config%active) THEN
            CALL veg_tile%Set_process_action(DISTURB_, ON_LEAFS_, AGGREGATE_)
          END IF

            DO i=1,npft
              pft_tile => t_jsb_tile(pft_shortnames(i), TRIM(pft_longnames(i)), &
               & lct=t_jsb_lct(VEG_TYPE, lib_id=pft_ids(i)), parent=veg_tile)
            END DO

    ! ! Next 3 lines only to avoid compiler warnings about pointers not being dereferenced
    ! IF (ASSOCIATED(glacier_tile)) CONTINUE
    ! IF (ASSOCIATED(lake_tile)) CONTINUE
    ! IF (ASSOCIATED(pft_tile)) CONTINUE
        
    NULLIFY(box_tile, land_tile, veg_tile, pft_tile)
    IF (model%config%use_glacier) NULLIFY(glacier_tile)
    IF (model%config%use_lakes) NULLIFY(lake_tile)

  END SUBROUTINE init_usecase_pfts

  SUBROUTINE init_usecase_jedi(model)

    TYPE(t_jsb_model), POINTER, INTENT(inout) :: model

    CLASS(t_jsb_tile), POINTER :: box_tile, lake_tile, land_tile, glacier_tile, veg_tile, pft_tile

    INTEGER :: i

    INTEGER :: npft

    CHARACTER(len=6), PARAMETER :: pft_shortnames(600) = [character(len=6) ::  &  
      &'pft001','pft002','pft003','pft004','pft005',&
      &'pft006','pft007','pft008','pft009','pft010',&
      &'pft011','pft012','pft013','pft014','pft015',&
      &'pft016','pft017','pft018','pft019','pft020',&
      &'pft021','pft022','pft023','pft024','pft025',&
      &'pft026','pft027','pft028','pft029','pft030',&
      &'pft031','pft032','pft033','pft034','pft035',&
      &'pft036','pft037','pft038','pft039','pft040',&
      &'pft041','pft042','pft043','pft044','pft045',&
      &'pft046','pft047','pft048','pft049','pft050',&
      &'pft051','pft052','pft053','pft054','pft055',&
      &'pft056','pft057','pft058','pft059','pft060',&
      &'pft061','pft062','pft063','pft064','pft065',&
      &'pft066','pft067','pft068','pft069','pft070',&
      &'pft071','pft072','pft073','pft074','pft075',&
      &'pft076','pft077','pft078','pft079','pft080',&
      &'pft081','pft082','pft083','pft084','pft085',&
      &'pft086','pft087','pft088','pft089','pft090',&
      &'pft091','pft092','pft093','pft094','pft095',&
      &'pft096','pft097','pft098','pft099','pft100',&
      &'pft101','pft102','pft103','pft104','pft105',&
      &'pft106','pft107','pft108','pft109','pft110',&
      &'pft111','pft112','pft113','pft114','pft115',&
      &'pft116','pft117','pft118','pft119','pft120',&
      &'pft121','pft122','pft123','pft124','pft125',&
      &'pft126','pft127','pft128','pft129','pft130',&
      &'pft131','pft132','pft133','pft134','pft135',&
      &'pft136','pft137','pft138','pft139','pft140',&
      &'pft141','pft142','pft143','pft144','pft145',&
      &'pft146','pft147','pft148','pft149','pft150',&
      &'pft151','pft152','pft153','pft154','pft155',&
      &'pft156','pft157','pft158','pft159','pft160',&
      &'pft161','pft162','pft163','pft164','pft165',&
      &'pft166','pft167','pft168','pft169','pft170',&
      &'pft171','pft172','pft173','pft174','pft175',&
      &'pft176','pft177','pft178','pft179','pft180',&
      &'pft181','pft182','pft183','pft184','pft185',&
      &'pft186','pft187','pft188','pft189','pft190',&
      &'pft191','pft192','pft193','pft194','pft195',&
      &'pft196','pft197','pft198','pft199','pft200',&
      &'pft201','pft202','pft203','pft204','pft205',&
      &'pft206','pft207','pft208','pft209','pft210',&
      &'pft211','pft212','pft213','pft214','pft215',&
      &'pft216','pft217','pft218','pft219','pft220',&
      &'pft221','pft222','pft223','pft224','pft225',&
      &'pft226','pft227','pft228','pft229','pft230',&
      &'pft231','pft232','pft233','pft234','pft235',&
      &'pft236','pft237','pft238','pft239','pft240',&
      &'pft241','pft242','pft243','pft244','pft245',&
      &'pft246','pft247','pft248','pft249','pft250',&      
      &'pft251','pft252','pft253','pft254','pft255',&
      &'pft256','pft257','pft258','pft259','pft260',&
       &'pft261','pft262','pft263','pft264','pft265',&
       &'pft266','pft267','pft268','pft269','pft270',&
       &'pft271','pft272','pft273','pft274','pft275',&
       &'pft276','pft277','pft278','pft279','pft280',&
       &'pft281','pft282','pft283','pft284','pft285',&
       &'pft286','pft287','pft288','pft289','pft290',&
       &'pft291','pft292','pft293','pft294','pft295',&
       &'pft296','pft297','pft298','pft299','pft300',&
       &'pft301','pft302','pft303','pft304','pft305',&
       &'pft306','pft307','pft308','pft309','pft310',&
       &'pft311','pft312','pft313','pft314','pft315',&
       &'pft316','pft317','pft318','pft319','pft320',&
       &'pft321','pft322','pft323','pft324','pft325',&
       &'pft326','pft327','pft328','pft329','pft330',&
       &'pft331','pft332','pft333','pft334','pft335',&
       &'pft336','pft337','pft338','pft339','pft340',&
       &'pft341','pft342','pft343','pft344','pft345',&
       &'pft346','pft347','pft348','pft349','pft350',&
       &'pft351','pft352','pft353','pft354','pft355',&
       &'pft356','pft357','pft358','pft359','pft360',&
       &'pft361','pft362','pft363','pft364','pft365',&
       &'pft366','pft367','pft368','pft369','pft370',&
       &'pft371','pft372','pft373','pft374','pft375',&
       &'pft376','pft377','pft378','pft379','pft380',&
       &'pft381','pft382','pft383','pft384','pft385',&
       &'pft386','pft387','pft388','pft389','pft390',&
       &'pft391','pft392','pft393','pft394','pft395',&
       &'pft396','pft397','pft398','pft399','pft400',&
       &'pft401','pft402','pft403','pft404','pft405',&
       &'pft406','pft407','pft408','pft409','pft410',&
       &'pft411','pft412','pft413','pft414','pft415',&
       &'pft416','pft417','pft418','pft419','pft420',&
       &'pft421','pft422','pft423','pft424','pft425',&
       &'pft426','pft427','pft428','pft429','pft430',&
       &'pft431','pft432','pft433','pft434','pft435',&
       &'pft436','pft437','pft438','pft439','pft440',&
       &'pft441','pft442','pft443','pft444','pft445',&
       &'pft446','pft447','pft448','pft449','pft450',&
       &'pft451','pft452','pft453','pft454','pft455',&
       &'pft456','pft457','pft458','pft459','pft460',&
       &'pft461','pft462','pft463','pft464','pft465',&
       &'pft466','pft467','pft468','pft469','pft470',&
       &'pft471','pft472','pft473','pft474','pft475',&
       &'pft476','pft477','pft478','pft479','pft480',&
       &'pft481','pft482','pft483','pft484','pft485',&
       &'pft486','pft487','pft488','pft489','pft490',&
       &'pft491','pft492','pft493','pft494','pft495',&
       &'pft496','pft497','pft498','pft499','pft500',&!]
       &'pft501','pft502','pft503','pft504','pft505',&
       &'pft506','pft507','pft508','pft509','pft510',&
       &'pft511','pft512','pft513','pft514','pft515',&
       &'pft516','pft517','pft518','pft519','pft520',&
       &'pft521','pft522','pft523','pft524','pft525',&
       &'pft526','pft527','pft528','pft529','pft530',&
       &'pft531','pft532','pft533','pft534','pft535',&
       &'pft536','pft537','pft538','pft539','pft540',&
       &'pft541','pft542','pft543','pft544','pft545',&
       &'pft546','pft547','pft548','pft549','pft550',&
       &'pft551','pft552','pft553','pft554','pft555',&
       &'pft556','pft557','pft558','pft559','pft560',&
       &'pft561','pft562','pft563','pft564','pft565',&
       &'pft566','pft567','pft568','pft569','pft570',&
       &'pft571','pft572','pft573','pft574','pft575',&
       &'pft576','pft577','pft578','pft579','pft580',&
       &'pft581','pft582','pft583','pft584','pft585',&
       &'pft586','pft587','pft588','pft589','pft590',&
       &'pft591','pft592','pft593','pft594','pft595',&
       &'pft596','pft597','pft598','pft599','pft600']

!       &'pft01', 'pft02', 'pft03', 'pft04', 'pft05',&
!      &'pft06', 'pft07', 'pft08', 'pft09', 'pft10',&
!      &'pft11', 'pft12', 'pft13', 'pft14', 'pft15',&
!      &'pft16', 'pft17', 'pft18', 'pft19', 'pft20',&
!      &'pft21', 'pft22', 'pft23', 'pft24', 'pft25',&
!      &'pft26', 'pft27', 'pft28', 'pft29', 'pft30',&
!      &'pft31', 'pft32', 'pft33', 'pft34', 'pft35',&
!      &'pft36', 'pft37', 'pft38', 'pft39', 'pft40']
                    
      
      

    CHARACTER(len=30), PARAMETER :: pft_longnames(600) = [character(len=30) :: &
      & 'Plant Growth Strategy 01',                                             &
      & 'Plant Growth Strategy 02',                                             &
      & 'Plant Growth Strategy 03',                                             &
      & 'Plant Growth Strategy 04',                                             &
      & 'Plant Growth Strategy 05',                                             &
      & 'Plant Growth Strategy 06',                                             &
      & 'Plant Growth Strategy 07',                                             &
      & 'Plant Growth Strategy 08',                                             &
      & 'Plant Growth Strategy 09',                                             &
      & 'Plant Growth Strategy 10',                                             &
      & 'Plant Growth Strategy 11',                                             &
      & 'Plant Growth Strategy 12',                                             &
      & 'Plant Growth Strategy 13',                                             &
      & 'Plant Growth Strategy 14',                                             &
      & 'Plant Growth Strategy 15',                                             &
      & 'Plant Growth Strategy 16',                                             &
      & 'Plant Growth Strategy 17',                                             &
      & 'Plant Growth Strategy 18',                                             &
      & 'Plant Growth Strategy 19',                                             &
      & 'Plant Growth Strategy 20',&
      & 'Plant Growth Strategy 21',                                             &
      & 'Plant Growth Strategy 22',                                             &
      & 'Plant Growth Strategy 23',                                             &
      & 'Plant Growth Strategy 24',                                             &
      & 'Plant Growth Strategy 25',                                             &
      & 'Plant Growth Strategy 26',                                             &
      & 'Plant Growth Strategy 27',                                             &
      & 'Plant Growth Strategy 28',                                             &
      & 'Plant Growth Strategy 29',                                             &
      & 'Plant Growth Strategy 30',  &
      & 'Plant Growth Strategy 31',                                             &
      & 'Plant Growth Strategy 32',                                             &
      & 'Plant Growth Strategy 33',                                             &
      & 'Plant Growth Strategy 34',                                             &
      & 'Plant Growth Strategy 35',                                             &
      & 'Plant Growth Strategy 36',                                             &
      & 'Plant Growth Strategy 37',                                             &
      & 'Plant Growth Strategy 38',                                             &
      & 'Plant Growth Strategy 39',                                             &
      & 'Plant Growth Strategy 40',&
      &'Plant Growth Strategy 041',&
      &'Plant Growth Strategy 042',&
      &'Plant Growth Strategy 043',&
      &'Plant Growth Strategy 044',&
      &'Plant Growth Strategy 045',&
      &'Plant Growth Strategy 046',&
      &'Plant Growth Strategy 047',&
      &'Plant Growth Strategy 048',&
      &'Plant Growth Strategy 049',&
      &'Plant Growth Strategy 050',&
      &'Plant Growth Strategy 051',&
      &'Plant Growth Strategy 052',&
      &'Plant Growth Strategy 053',&
      &'Plant Growth Strategy 054',&
      &'Plant Growth Strategy 055',&
      &'Plant Growth Strategy 056',&
      &'Plant Growth Strategy 057',&
      &'Plant Growth Strategy 058',&
      &'Plant Growth Strategy 059',&
      &'Plant Growth Strategy 060',&
      &'Plant Growth Strategy 061',&
      &'Plant Growth Strategy 062',&
      &'Plant Growth Strategy 063',&
      &'Plant Growth Strategy 064',&
      &'Plant Growth Strategy 065',&
      &'Plant Growth Strategy 066',&
      &'Plant Growth Strategy 067',&
      &'Plant Growth Strategy 068',&
      &'Plant Growth Strategy 069',&
      &'Plant Growth Strategy 070',&
      &'Plant Growth Strategy 071',&
      &'Plant Growth Strategy 072',&
      &'Plant Growth Strategy 073',&
      &'Plant Growth Strategy 074',&
      &'Plant Growth Strategy 075',&
      &'Plant Growth Strategy 076',&
      &'Plant Growth Strategy 077',&
      &'Plant Growth Strategy 078',&
      &'Plant Growth Strategy 079',&
      &'Plant Growth Strategy 080',&
      &'Plant Growth Strategy 081',&
      &'Plant Growth Strategy 082',&
      &'Plant Growth Strategy 083',&
      &'Plant Growth Strategy 084',&
      &'Plant Growth Strategy 085',&
      &'Plant Growth Strategy 086',&
      &'Plant Growth Strategy 087',&
      &'Plant Growth Strategy 088',&
      &'Plant Growth Strategy 089',&
      &'Plant Growth Strategy 090',&
      &'Plant Growth Strategy 091',&
      &'Plant Growth Strategy 092',&
      &'Plant Growth Strategy 093',&
      &'Plant Growth Strategy 094',&
      &'Plant Growth Strategy 095',&
      &'Plant Growth Strategy 096',&
      &'Plant Growth Strategy 097',&
      &'Plant Growth Strategy 098',&
      &'Plant Growth Strategy 099',&
      &'Plant Growth Strategy 100',&
      &'Plant Growth Strategy 101',&
      &'Plant Growth Strategy 102',&
      &'Plant Growth Strategy 103',&
      &'Plant Growth Strategy 104',&
      &'Plant Growth Strategy 105',&
      &'Plant Growth Strategy 106',&
      &'Plant Growth Strategy 107',&
      &'Plant Growth Strategy 108',&
      &'Plant Growth Strategy 109',&
      &'Plant Growth Strategy 110',&
      &'Plant Growth Strategy 111',&
      &'Plant Growth Strategy 112',&
      &'Plant Growth Strategy 113',&
      &'Plant Growth Strategy 114',&
      &'Plant Growth Strategy 115',&
      &'Plant Growth Strategy 116',&
      &'Plant Growth Strategy 117',&
      &'Plant Growth Strategy 118',&
      &'Plant Growth Strategy 119',&
      &'Plant Growth Strategy 120',&
      &'Plant Growth Strategy 121',&
      &'Plant Growth Strategy 122',&
      &'Plant Growth Strategy 123',&
      &'Plant Growth Strategy 124',&
      &'Plant Growth Strategy 125',&
      &'Plant Growth Strategy 126',&
      &'Plant Growth Strategy 127',&
      &'Plant Growth Strategy 128',&
      &'Plant Growth Strategy 129',&
      &'Plant Growth Strategy 130',&
      &'Plant Growth Strategy 131',&
      &'Plant Growth Strategy 132',&
      &'Plant Growth Strategy 133',&
      &'Plant Growth Strategy 134',&
      &'Plant Growth Strategy 135',&
      &'Plant Growth Strategy 136',&
      &'Plant Growth Strategy 137',&
      &'Plant Growth Strategy 138',&
      &'Plant Growth Strategy 139',&
      &'Plant Growth Strategy 140',&
      &'Plant Growth Strategy 141',&
      &'Plant Growth Strategy 142',&
      &'Plant Growth Strategy 143',&
      &'Plant Growth Strategy 144',&
      &'Plant Growth Strategy 145',&
      &'Plant Growth Strategy 146',&
      &'Plant Growth Strategy 147',&
      &'Plant Growth Strategy 148',&
      &'Plant Growth Strategy 149',&
      &'Plant Growth Strategy 150',&
      &'Plant Growth Strategy 151',&
      &'Plant Growth Strategy 152',&
      &'Plant Growth Strategy 153',&
      &'Plant Growth Strategy 154',&
      &'Plant Growth Strategy 155',&
      &'Plant Growth Strategy 156',&
      &'Plant Growth Strategy 157',&
      &'Plant Growth Strategy 158',&
      &'Plant Growth Strategy 159',&
      &'Plant Growth Strategy 160',&
      &'Plant Growth Strategy 161',&
      &'Plant Growth Strategy 162',&
      &'Plant Growth Strategy 163',&
      &'Plant Growth Strategy 164',&
      &'Plant Growth Strategy 165',&
      &'Plant Growth Strategy 166',&
      &'Plant Growth Strategy 167',&
      &'Plant Growth Strategy 168',&
      &'Plant Growth Strategy 169',&
      &'Plant Growth Strategy 170',&
      &'Plant Growth Strategy 171',&
      &'Plant Growth Strategy 172',&
      &'Plant Growth Strategy 173',&
      &'Plant Growth Strategy 174',&
      &'Plant Growth Strategy 175',&
      &'Plant Growth Strategy 176',&
      &'Plant Growth Strategy 177',&
      &'Plant Growth Strategy 178',&
      &'Plant Growth Strategy 179',&
      &'Plant Growth Strategy 180',&
      &'Plant Growth Strategy 181',&
      &'Plant Growth Strategy 182',&
      &'Plant Growth Strategy 183',&
      &'Plant Growth Strategy 184',&
      &'Plant Growth Strategy 185',&
      &'Plant Growth Strategy 186',&
      &'Plant Growth Strategy 187',&
      &'Plant Growth Strategy 188',&
      &'Plant Growth Strategy 189',&
      &'Plant Growth Strategy 190',&
      &'Plant Growth Strategy 191',&
      &'Plant Growth Strategy 192',&
      &'Plant Growth Strategy 193',&
      &'Plant Growth Strategy 194',&
      &'Plant Growth Strategy 195',&
      &'Plant Growth Strategy 196',&
      &'Plant Growth Strategy 197',&
      &'Plant Growth Strategy 198',&
      &'Plant Growth Strategy 199',&
      &'Plant Growth Strategy 200',&
      &'Plant Growth Strategy 201',&
      &'Plant Growth Strategy 202',&
      &'Plant Growth Strategy 203',&
      &'Plant Growth Strategy 204',&
      &'Plant Growth Strategy 205',&
      &'Plant Growth Strategy 206',&
      &'Plant Growth Strategy 207',&
      &'Plant Growth Strategy 208',&
      &'Plant Growth Strategy 209',&
      &'Plant Growth Strategy 210',&
      &'Plant Growth Strategy 211',&
      &'Plant Growth Strategy 212',&
      &'Plant Growth Strategy 213',&
      &'Plant Growth Strategy 214',&
      &'Plant Growth Strategy 215',&
      &'Plant Growth Strategy 216',&
      &'Plant Growth Strategy 217',&
      &'Plant Growth Strategy 218',&
      &'Plant Growth Strategy 219',&
      &'Plant Growth Strategy 220',&
      &'Plant Growth Strategy 221',&
      &'Plant Growth Strategy 222',&
      &'Plant Growth Strategy 223',&
      &'Plant Growth Strategy 224',&
      &'Plant Growth Strategy 225',&
      &'Plant Growth Strategy 226',&
      &'Plant Growth Strategy 227',&
      &'Plant Growth Strategy 228',&
      &'Plant Growth Strategy 229',&
      &'Plant Growth Strategy 230',&
      &'Plant Growth Strategy 231',&
      &'Plant Growth Strategy 232',&
      &'Plant Growth Strategy 233',&
      &'Plant Growth Strategy 234',&
      &'Plant Growth Strategy 235',&
      &'Plant Growth Strategy 236',&
      &'Plant Growth Strategy 237',&
      &'Plant Growth Strategy 238',&
      &'Plant Growth Strategy 239',&
      &'Plant Growth Strategy 240',&
      &'Plant Growth Strategy 241',&
      &'Plant Growth Strategy 242',&
      &'Plant Growth Strategy 243',&
      &'Plant Growth Strategy 244',&
      &'Plant Growth Strategy 245',&
      &'Plant Growth Strategy 246',&
      &'Plant Growth Strategy 247',&
      &'Plant Growth Strategy 248',&
      &'Plant Growth Strategy 249',&
      &'Plant Growth Strategy 250',&!]      
      &'Plant Growth Strategy 251',&
      &'Plant Growth Strategy 252',&
      &'Plant Growth Strategy 253',&
      &'Plant Growth Strategy 254',&
      &'Plant Growth Strategy 255',&
      &'Plant Growth Strategy 256',&
      &'Plant Growth Strategy 257',&
      &'Plant Growth Strategy 258',&
      &'Plant Growth Strategy 259',&
      &'Plant Growth Strategy 260',&
      &'Plant Growth Strategy 261',&
      &'Plant Growth Strategy 262',&
      &'Plant Growth Strategy 263',&
      &'Plant Growth Strategy 264',&
      &'Plant Growth Strategy 265',&
      &'Plant Growth Strategy 266',&
      &'Plant Growth Strategy 267',&
      &'Plant Growth Strategy 268',&
      &'Plant Growth Strategy 269',&
      &'Plant Growth Strategy 270',&
      &'Plant Growth Strategy 271',&
      &'Plant Growth Strategy 272',&
      &'Plant Growth Strategy 273',&
      &'Plant Growth Strategy 274',&
      &'Plant Growth Strategy 275',&
      &'Plant Growth Strategy 276',&
      &'Plant Growth Strategy 277',&
      &'Plant Growth Strategy 278',&
      &'Plant Growth Strategy 279',&
      &'Plant Growth Strategy 280',&
      &'Plant Growth Strategy 281',&
      &'Plant Growth Strategy 282',&
      &'Plant Growth Strategy 283',&
      &'Plant Growth Strategy 284',&
      &'Plant Growth Strategy 285',&
      &'Plant Growth Strategy 286',&
      &'Plant Growth Strategy 287',&
      &'Plant Growth Strategy 288',&
      &'Plant Growth Strategy 289',&
      &'Plant Growth Strategy 290',&
      &'Plant Growth Strategy 291',&
      &'Plant Growth Strategy 292',&
      &'Plant Growth Strategy 293',&
      &'Plant Growth Strategy 294',&
      &'Plant Growth Strategy 295',&
      &'Plant Growth Strategy 296',&
      &'Plant Growth Strategy 297',&
      &'Plant Growth Strategy 298',&
      &'Plant Growth Strategy 299',&
      &'Plant Growth Strategy 300',&
      &'Plant Growth Strategy 301',&
      &'Plant Growth Strategy 302',&
      &'Plant Growth Strategy 303',&
      &'Plant Growth Strategy 304',&
      &'Plant Growth Strategy 305',&
      &'Plant Growth Strategy 306',&
      &'Plant Growth Strategy 307',&
      &'Plant Growth Strategy 308',&
      &'Plant Growth Strategy 309',&
      &'Plant Growth Strategy 310',&
      &'Plant Growth Strategy 311',&
      &'Plant Growth Strategy 312',&
      &'Plant Growth Strategy 313',&
      &'Plant Growth Strategy 314',&
      &'Plant Growth Strategy 315',&
      &'Plant Growth Strategy 316',&
      &'Plant Growth Strategy 317',&
      &'Plant Growth Strategy 318',&
      &'Plant Growth Strategy 319',&
      &'Plant Growth Strategy 320',&
      &'Plant Growth Strategy 321',&
      &'Plant Growth Strategy 322',&
      &'Plant Growth Strategy 323',&
      &'Plant Growth Strategy 324',&
      &'Plant Growth Strategy 325',&
      &'Plant Growth Strategy 326',&
      &'Plant Growth Strategy 327',&
      &'Plant Growth Strategy 328',&
      &'Plant Growth Strategy 329',&
      &'Plant Growth Strategy 330',&
      &'Plant Growth Strategy 331',&
      &'Plant Growth Strategy 332',&
      &'Plant Growth Strategy 333',&
      &'Plant Growth Strategy 334',&
      &'Plant Growth Strategy 335',&
      &'Plant Growth Strategy 336',&
      &'Plant Growth Strategy 337',&
      &'Plant Growth Strategy 338',&
      &'Plant Growth Strategy 339',&
      &'Plant Growth Strategy 340',&
      &'Plant Growth Strategy 341',&
      &'Plant Growth Strategy 342',&
      &'Plant Growth Strategy 343',&
      &'Plant Growth Strategy 344',&
      &'Plant Growth Strategy 345',&
      &'Plant Growth Strategy 346',&
      &'Plant Growth Strategy 347',&
      &'Plant Growth Strategy 348',&
      &'Plant Growth Strategy 349',&
      &'Plant Growth Strategy 350',&
      &'Plant Growth Strategy 351',&
      &'Plant Growth Strategy 352',&
      &'Plant Growth Strategy 353',&
      &'Plant Growth Strategy 354',&
      &'Plant Growth Strategy 355',&
      &'Plant Growth Strategy 356',&
      &'Plant Growth Strategy 357',&
      &'Plant Growth Strategy 358',&
      &'Plant Growth Strategy 359',&
      &'Plant Growth Strategy 360',&
      &'Plant Growth Strategy 361',&
      &'Plant Growth Strategy 362',&
      &'Plant Growth Strategy 363',&
      &'Plant Growth Strategy 364',&
      &'Plant Growth Strategy 365',&
      &'Plant Growth Strategy 366',&
      &'Plant Growth Strategy 367',&
      &'Plant Growth Strategy 368',&
      &'Plant Growth Strategy 369',&
      &'Plant Growth Strategy 370',&
      &'Plant Growth Strategy 371',&
      &'Plant Growth Strategy 372',&
      &'Plant Growth Strategy 373',&
      &'Plant Growth Strategy 374',&
      &'Plant Growth Strategy 375',&
      &'Plant Growth Strategy 376',&
      &'Plant Growth Strategy 377',&
      &'Plant Growth Strategy 378',&
      &'Plant Growth Strategy 379',&
      &'Plant Growth Strategy 380',&
      &'Plant Growth Strategy 381',&
      &'Plant Growth Strategy 382',&
      &'Plant Growth Strategy 383',&
      &'Plant Growth Strategy 384',&
      &'Plant Growth Strategy 385',&
      &'Plant Growth Strategy 386',&
      &'Plant Growth Strategy 387',&
      &'Plant Growth Strategy 388',&
      &'Plant Growth Strategy 389',&
      &'Plant Growth Strategy 390',&
      &'Plant Growth Strategy 391',&
      &'Plant Growth Strategy 392',&
      &'Plant Growth Strategy 393',&
      &'Plant Growth Strategy 394',&
      &'Plant Growth Strategy 395',&
      &'Plant Growth Strategy 396',&
      &'Plant Growth Strategy 397',&
      &'Plant Growth Strategy 398',&
      &'Plant Growth Strategy 399',&
      &'Plant Growth Strategy 400',&
      &'Plant Growth Strategy 401',&
      &'Plant Growth Strategy 402',&
      &'Plant Growth Strategy 403',&
      &'Plant Growth Strategy 404',&
      &'Plant Growth Strategy 405',&
      &'Plant Growth Strategy 406',&
      &'Plant Growth Strategy 407',&
      &'Plant Growth Strategy 408',&
      &'Plant Growth Strategy 409',&
      &'Plant Growth Strategy 410',&
      &'Plant Growth Strategy 411',&
      &'Plant Growth Strategy 412',&
      &'Plant Growth Strategy 413',&
      &'Plant Growth Strategy 414',&
      &'Plant Growth Strategy 415',&
      &'Plant Growth Strategy 416',&
      &'Plant Growth Strategy 417',&
      &'Plant Growth Strategy 418',&
      &'Plant Growth Strategy 419',&
      &'Plant Growth Strategy 420',&
      &'Plant Growth Strategy 421',&
      &'Plant Growth Strategy 422',&
      &'Plant Growth Strategy 423',&
      &'Plant Growth Strategy 424',&
      &'Plant Growth Strategy 425',&
      &'Plant Growth Strategy 426',&
      &'Plant Growth Strategy 427',&
      &'Plant Growth Strategy 428',&
      &'Plant Growth Strategy 429',&
      &'Plant Growth Strategy 430',&
      &'Plant Growth Strategy 431',&
      &'Plant Growth Strategy 432',&
      &'Plant Growth Strategy 433',&
      &'Plant Growth Strategy 434',&
      &'Plant Growth Strategy 435',&
      &'Plant Growth Strategy 436',&
      &'Plant Growth Strategy 437',&
      &'Plant Growth Strategy 438',&
      &'Plant Growth Strategy 439',&
      &'Plant Growth Strategy 440',&
      &'Plant Growth Strategy 441',&
      &'Plant Growth Strategy 442',&
      &'Plant Growth Strategy 443',&
      &'Plant Growth Strategy 444',&
      &'Plant Growth Strategy 445',&
      &'Plant Growth Strategy 446',&
      &'Plant Growth Strategy 447',&
      &'Plant Growth Strategy 448',&
      &'Plant Growth Strategy 449',&
      &'Plant Growth Strategy 450',&
      &'Plant Growth Strategy 451',&
      &'Plant Growth Strategy 452',&
      &'Plant Growth Strategy 453',&
      &'Plant Growth Strategy 454',&
      &'Plant Growth Strategy 455',&
      &'Plant Growth Strategy 456',&
      &'Plant Growth Strategy 457',&
      &'Plant Growth Strategy 458',&
      &'Plant Growth Strategy 459',&
      &'Plant Growth Strategy 460',&
      &'Plant Growth Strategy 461',&
      &'Plant Growth Strategy 462',&
      &'Plant Growth Strategy 463',&
      &'Plant Growth Strategy 464',&
      &'Plant Growth Strategy 465',&
      &'Plant Growth Strategy 466',&
      &'Plant Growth Strategy 467',&
      &'Plant Growth Strategy 468',&
      &'Plant Growth Strategy 469',&
      &'Plant Growth Strategy 470',&
      &'Plant Growth Strategy 471',&
      &'Plant Growth Strategy 472',&
      &'Plant Growth Strategy 473',&
      &'Plant Growth Strategy 474',&
      &'Plant Growth Strategy 475',&
      &'Plant Growth Strategy 476',&
      &'Plant Growth Strategy 477',&
      &'Plant Growth Strategy 478',&
      &'Plant Growth Strategy 479',&
      &'Plant Growth Strategy 480',&
      &'Plant Growth Strategy 481',&
      &'Plant Growth Strategy 482',&
      &'Plant Growth Strategy 483',&
      &'Plant Growth Strategy 484',&
      &'Plant Growth Strategy 485',&
      &'Plant Growth Strategy 486',&
      &'Plant Growth Strategy 487',&
      &'Plant Growth Strategy 488',&
      &'Plant Growth Strategy 489',&
      &'Plant Growth Strategy 490',&
      &'Plant Growth Strategy 491',&
      &'Plant Growth Strategy 492',&
      &'Plant Growth Strategy 493',&
      &'Plant Growth Strategy 494',&
      &'Plant Growth Strategy 495',&
      &'Plant Growth Strategy 496',&
      &'Plant Growth Strategy 497',&
      &'Plant Growth Strategy 498',&
      &'Plant Growth Strategy 499',&
      &'Plant Growth Strategy 500',&!]
      &'Plant Growth Strategy 501',&
      &'Plant Growth Strategy 502',&
      &'Plant Growth Strategy 503',&
      &'Plant Growth Strategy 504',&
      &'Plant Growth Strategy 505',&
      &'Plant Growth Strategy 506',&
      &'Plant Growth Strategy 507',&
      &'Plant Growth Strategy 508',&
      &'Plant Growth Strategy 509',&
      &'Plant Growth Strategy 510',&
      &'Plant Growth Strategy 511',&
      &'Plant Growth Strategy 512',&
      &'Plant Growth Strategy 513',&
      &'Plant Growth Strategy 514',&
      &'Plant Growth Strategy 515',&
      &'Plant Growth Strategy 516',&
      &'Plant Growth Strategy 517',&
      &'Plant Growth Strategy 518',&
      &'Plant Growth Strategy 519',&
      &'Plant Growth Strategy 520',&
      &'Plant Growth Strategy 521',&
      &'Plant Growth Strategy 522',&
      &'Plant Growth Strategy 523',&
      &'Plant Growth Strategy 524',&
      &'Plant Growth Strategy 525',&
      &'Plant Growth Strategy 526',&
      &'Plant Growth Strategy 527',&
      &'Plant Growth Strategy 528',&
      &'Plant Growth Strategy 529',&
      &'Plant Growth Strategy 530',&
      &'Plant Growth Strategy 531',&
      &'Plant Growth Strategy 532',&
      &'Plant Growth Strategy 533',&
      &'Plant Growth Strategy 534',&
      &'Plant Growth Strategy 535',&
      &'Plant Growth Strategy 536',&
      &'Plant Growth Strategy 537',&
      &'Plant Growth Strategy 538',&
      &'Plant Growth Strategy 539',&
      &'Plant Growth Strategy 540',&
      &'Plant Growth Strategy 541',&
      &'Plant Growth Strategy 542',&
      &'Plant Growth Strategy 543',&
      &'Plant Growth Strategy 544',&
      &'Plant Growth Strategy 545',&
      &'Plant Growth Strategy 546',&
      &'Plant Growth Strategy 547',&
      &'Plant Growth Strategy 548',&
      &'Plant Growth Strategy 549',&
      &'Plant Growth Strategy 550',&
      &'Plant Growth Strategy 551',&
      &'Plant Growth Strategy 552',&
      &'Plant Growth Strategy 553',&
      &'Plant Growth Strategy 554',&
      &'Plant Growth Strategy 555',&
      &'Plant Growth Strategy 556',&
      &'Plant Growth Strategy 557',&
      &'Plant Growth Strategy 558',&
      &'Plant Growth Strategy 559',&
      &'Plant Growth Strategy 560',&
      &'Plant Growth Strategy 561',&
      &'Plant Growth Strategy 562',&
      &'Plant Growth Strategy 563',&
      &'Plant Growth Strategy 564',&
      &'Plant Growth Strategy 565',&
      &'Plant Growth Strategy 566',&
      &'Plant Growth Strategy 567',&
      &'Plant Growth Strategy 568',&
      &'Plant Growth Strategy 569',&
      &'Plant Growth Strategy 570',&
      &'Plant Growth Strategy 571',&
      &'Plant Growth Strategy 572',&
      &'Plant Growth Strategy 573',&
      &'Plant Growth Strategy 574',&
      &'Plant Growth Strategy 575',&
      &'Plant Growth Strategy 576',&
      &'Plant Growth Strategy 577',&
      &'Plant Growth Strategy 578',&
      &'Plant Growth Strategy 579',&
      &'Plant Growth Strategy 580',&
      &'Plant Growth Strategy 581',&
      &'Plant Growth Strategy 582',&
      &'Plant Growth Strategy 583',&
      &'Plant Growth Strategy 584',&
      &'Plant Growth Strategy 585',&
      &'Plant Growth Strategy 586',&
      &'Plant Growth Strategy 587',&
      &'Plant Growth Strategy 588',&
      &'Plant Growth Strategy 589',&
      &'Plant Growth Strategy 590',&
      &'Plant Growth Strategy 591',&
      &'Plant Growth Strategy 592',&
      &'Plant Growth Strategy 593',&
      &'Plant Growth Strategy 594',&
      &'Plant Growth Strategy 595',&
      &'Plant Growth Strategy 596',&
      &'Plant Growth Strategy 597',&
      &'Plant Growth Strategy 598',&
      &'Plant Growth Strategy 599',&
      &'Plant Growth Strategy 600']!,&
   
      

    ! Index of PFT into lctlib columns
    ! id is one less than in lctlib file because the first column (glacier) has been thrown away
    INTEGER, ALLOCATABLE :: pft_ids(:)

    CHARACTER(len=*), PARAMETER :: routine = modname//':init_usecase_pfts'

    npft = model%config%npft
    IF (npft > 600) CALL finish(routine, 'You have to edit mo_jsb_model_usecases.f90 to use npft > 600')
    ALLOCATE(pft_ids(npft))
    DO i=1,npft
      pft_ids(i) = i
    END DO

    ! Create root of the tile structure representing the whole grid box
    box_tile => t_jsb_tile('box', 'Top tile', &
      & processes      =(/A2L_,     L2A_,     RAD_  ,    TURB_    /), &
      & process_actions=(/ON_TILE_, ON_TILE_, ON_LEAFS_, ON_LEAFS_/), &
      & model_id=model%id, fract_varname='notsea')
    model%top => box_tile

#ifndef __NO_JSBACH_HD__
    IF (model%processes(HD_)%p%config%active) THEN
      CALL box_tile%Set_process_action(HD_, ON_TILE_)
    END IF
#endif

      IF (model%config%use_lakes) THEN
        lake_tile => t_jsb_tile('lake', 'Lake tile', &
          & processes=(/SEB_, HYDRO_/), process_actions=(/ON_TILE_, ON_TILE_/), &
          & lct=t_jsb_lct(LAKE_TYPE), parent=box_tile)
        land_tile => t_jsb_tile('land', 'Land tile', &
          & processes=(/SEB_/), process_actions=(/ON_TILE_/), &
          & lct=t_jsb_lct(LAND_TYPE), parent=box_tile)
      ELSE
        land_tile => t_jsb_tile('land', 'Land tile', &
          & processes=(/SEB_/), process_actions=(/ON_TILE_/), &
          & lct=t_jsb_lct(LAND_TYPE), parent=box_tile, &
          & fract_varname='equal')   ! Fraction = 1.
      END IF

        IF (model%config%use_glacier) THEN
          glacier_tile => t_jsb_tile('glac', 'Glacier tile', &
            & processes = (/HYDRO_, SSE_/), &
            & process_actions = (/ON_TILE_, ON_TILE_/), &
            & lct=t_jsb_lct(GLACIER_TYPE), parent=land_tile)
            ! & lct=t_jsb_lct(GLACIER_TYPE, lib_id=1), parent=land_tile)
        END IF
    !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-    
    ! March 22 2019
    !PHS: THERE IS A BUG IN THE WATER LIMITED ASSIMILATION PROCESSES WHEN TAKING THE WATER STRESS FROM THR HYDRO PROCESS
!   !PHS:    SO RAINER CHANGE THE COMPUTAION OF THE HYDRO PROCESS FROM TILE LEVEL (ON_TILE_) TO PFT (ON_LEAFS_)
!           veg_tile => t_jsb_tile('veg', 'Veg tile', &
!            & processes       = (/HYDRO_,   SSE_,     ASSIMI_,   PHENO_   /), &
!            & process_actions = (/ON_TILE_, ON_TILE_, ON_LEAFS_, ON_LEAFS_/), &
!            & lct=t_jsb_lct(VEG_TYPE), parent=land_tile,  &
!            & fract_varname='fract_veg')
    !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

          veg_tile => t_jsb_tile('veg', 'Veg tile', &
            & processes       = (/HYDRO_,   SSE_,     ASSIMI_,   PHENO_   /), &
            & process_actions = (/ON_TILE_, ON_TILE_, ON_LEAFS_, ON_LEAFS_/), &
            & lct=t_jsb_lct(VEG_TYPE), parent=land_tile,  &
            & fract_varname='fract_veg')

          IF (model%processes(JEDI_)%p%config%active) THEN
            CALL veg_tile%Set_process_action(JEDI_, ON_LEAFS_, AGGREGATE_)
          END IF
       !   IF (model%processes(VC_)%p%config%active) THEN
       !     CALL veg_tile%Set_process_action(VC_, ON_LEAFS_, AGGREGATE_)
       !   END IF

            DO i=1,npft
              pft_tile => t_jsb_tile(pft_shortnames(i), TRIM(pft_longnames(i)), &
               & lct=t_jsb_lct(VEG_TYPE, lib_id=pft_ids(i)), parent=veg_tile,fract_varname='equal') 
            END DO

    NULLIFY(box_tile)
    NULLIFY(lake_tile)
    NULLIFY(land_tile)
    NULLIFY(glacier_tile)
    NULLIFY(veg_tile)
    NULLIFY(pft_tile)

  END SUBROUTINE init_usecase_jedi

#endif
END MODULE mo_jsb_model_usecases
