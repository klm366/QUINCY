! -------------------------------------------------------
! definition of structures for vegetation
! you can define/use pools either individually or all
! -------------------------------------------------------

#if defined(Def_veg_pool)
  TYPE(t_veg_pool_a)              :: veg_pool
#endif

#if defined(Def_veg_growth)
  TYPE(t_veg_pool_a)              :: veg_growth
#endif

#if defined(Def_veg_litterfall)
  TYPE(t_veg_pool_a)              :: veg_litterfall
#endif

#if defined(Def_veg_falloc)
  TYPE(t_veg_pool_a)              :: veg_falloc
#endif

#if defined(Def_veg_total_flux)
  TYPE(t_veg_pool_a)              :: veg_total_flux
#endif


#ifdef Def_veg_exudation
  TYPE(t_veg_element_a)           :: veg_exudation
#endif

#ifdef Def_veg_establishment
  TYPE(t_veg_element_a)           :: veg_establishment
#endif

#ifdef Def_veg_reserve_use
  TYPE(t_veg_element_a)           :: veg_reserve_use
#endif

#ifdef Def_veg_total_biomass
  TYPE(t_veg_element_a)           :: veg_total_biomass
#endif

! -------------------------------------------------------
! definition of structures for sb
! you can define/use pools either individually or all
! -------------------------------------------------------

#if defined(Def_sb_pool)
  TYPE(t_sb_pool_a)              :: sb_pool
#endif

#if defined(Def_sb_formation)
  TYPE(t_sb_pool_a)              :: sb_formation
#endif

#if defined(Def_sb_loss)
  TYPE(t_sb_pool_a)              :: sb_loss
#endif

#if defined(Def_sb_transport)
  TYPE(t_sb_pool_a)              :: sb_transport
#endif


#ifdef Def_sb_delta_pools
  TYPE(t_sb_element_a)           :: sb_delta_pools
#endif

#ifdef Def_sb_mycorrhiza_export
  TYPE(t_sb_element_a)           :: sb_mycorrhiza_export
#endif


! -------------------------------------------------------
! allocate variables in pools of vegetation memory
! -------------------------------------------------------

#if defined(Get_veg_pool)
  allocate(veg_pool%leaf%carbon              , source=veg__mem%pool%leaf%carbon%ptr)
  allocate(veg_pool%leaf%nitrogen            , source=veg__mem%pool%leaf%nitrogen%ptr)
  allocate(veg_pool%leaf%phosphorus          , source=veg__mem%pool%leaf%phosphorus%ptr)
  allocate(veg_pool%leaf%carbon13            , source=veg__mem%pool%leaf%carbon13%ptr)
  allocate(veg_pool%leaf%carbon14            , source=veg__mem%pool%leaf%carbon14%ptr)
  allocate(veg_pool%leaf%nitrogen15          , source=veg__mem%pool%leaf%nitrogen15%ptr)
  allocate(veg_pool%fine_root%carbon         , source=veg__mem%pool%fine_root%carbon%ptr)
  allocate(veg_pool%fine_root%nitrogen       , source=veg__mem%pool%fine_root%nitrogen%ptr)
  allocate(veg_pool%fine_root%phosphorus     , source=veg__mem%pool%fine_root%phosphorus%ptr)
  allocate(veg_pool%fine_root%carbon13       , source=veg__mem%pool%fine_root%carbon13%ptr)
  allocate(veg_pool%fine_root%carbon14       , source=veg__mem%pool%fine_root%carbon14%ptr)
  allocate(veg_pool%fine_root%nitrogen15     , source=veg__mem%pool%fine_root%nitrogen15%ptr)
  allocate(veg_pool%coarse_root%carbon       , source=veg__mem%pool%coarse_root%carbon%ptr)
  allocate(veg_pool%coarse_root%nitrogen     , source=veg__mem%pool%coarse_root%nitrogen%ptr)
  allocate(veg_pool%coarse_root%phosphorus   , source=veg__mem%pool%coarse_root%phosphorus%ptr)
  allocate(veg_pool%coarse_root%carbon13     , source=veg__mem%pool%coarse_root%carbon13%ptr)
  allocate(veg_pool%coarse_root%carbon14     , source=veg__mem%pool%coarse_root%carbon14%ptr)
  allocate(veg_pool%coarse_root%nitrogen15   , source=veg__mem%pool%coarse_root%nitrogen15%ptr)
  allocate(veg_pool%sap_wood%carbon          , source=veg__mem%pool%sap_wood%carbon%ptr)
  allocate(veg_pool%sap_wood%nitrogen        , source=veg__mem%pool%sap_wood%nitrogen%ptr)
  allocate(veg_pool%sap_wood%phosphorus      , source=veg__mem%pool%sap_wood%phosphorus%ptr)
  allocate(veg_pool%sap_wood%carbon13        , source=veg__mem%pool%sap_wood%carbon13%ptr)
  allocate(veg_pool%sap_wood%carbon14        , source=veg__mem%pool%sap_wood%carbon14%ptr)
  allocate(veg_pool%sap_wood%nitrogen15      , source=veg__mem%pool%sap_wood%nitrogen15%ptr)
  allocate(veg_pool%heart_wood%carbon        , source=veg__mem%pool%heart_wood%carbon%ptr)
  allocate(veg_pool%heart_wood%nitrogen      , source=veg__mem%pool%heart_wood%nitrogen%ptr)
  allocate(veg_pool%heart_wood%phosphorus    , source=veg__mem%pool%heart_wood%phosphorus%ptr)
  allocate(veg_pool%heart_wood%carbon13      , source=veg__mem%pool%heart_wood%carbon13%ptr)
  allocate(veg_pool%heart_wood%carbon14      , source=veg__mem%pool%heart_wood%carbon14%ptr)
  allocate(veg_pool%heart_wood%nitrogen15    , source=veg__mem%pool%heart_wood%nitrogen15%ptr)
  allocate(veg_pool%labile%carbon            , source=veg__mem%pool%labile%carbon%ptr)
  allocate(veg_pool%labile%nitrogen          , source=veg__mem%pool%labile%nitrogen%ptr)
  allocate(veg_pool%labile%phosphorus        , source=veg__mem%pool%labile%phosphorus%ptr)
  allocate(veg_pool%labile%carbon13          , source=veg__mem%pool%labile%carbon13%ptr)
  allocate(veg_pool%labile%carbon14          , source=veg__mem%pool%labile%carbon14%ptr)
  allocate(veg_pool%labile%nitrogen15        , source=veg__mem%pool%labile%nitrogen15%ptr)
  allocate(veg_pool%reserve%carbon           , source=veg__mem%pool%reserve%carbon%ptr)
  allocate(veg_pool%reserve%nitrogen         , source=veg__mem%pool%reserve%nitrogen%ptr)
  allocate(veg_pool%reserve%phosphorus       , source=veg__mem%pool%reserve%phosphorus%ptr)
  allocate(veg_pool%reserve%carbon13         , source=veg__mem%pool%reserve%carbon13%ptr)
  allocate(veg_pool%reserve%carbon14         , source=veg__mem%pool%reserve%carbon14%ptr)
  allocate(veg_pool%reserve%nitrogen15       , source=veg__mem%pool%reserve%nitrogen15%ptr)
  allocate(veg_pool%fruit%carbon             , source=veg__mem%pool%fruit%carbon%ptr)
  allocate(veg_pool%fruit%nitrogen           , source=veg__mem%pool%fruit%nitrogen%ptr)
  allocate(veg_pool%fruit%phosphorus         , source=veg__mem%pool%fruit%phosphorus%ptr)
  allocate(veg_pool%fruit%carbon13           , source=veg__mem%pool%fruit%carbon13%ptr)
  allocate(veg_pool%fruit%carbon14           , source=veg__mem%pool%fruit%carbon14%ptr)
  allocate(veg_pool%fruit%nitrogen15         , source=veg__mem%pool%fruit%nitrogen15%ptr)
  allocate(veg_pool%seed_bed%carbon          , source=veg__mem%pool%seed_bed%carbon%ptr)
  allocate(veg_pool%seed_bed%nitrogen        , source=veg__mem%pool%seed_bed%nitrogen%ptr)
  allocate(veg_pool%seed_bed%phosphorus      , source=veg__mem%pool%seed_bed%phosphorus%ptr)
  allocate(veg_pool%seed_bed%carbon13        , source=veg__mem%pool%seed_bed%carbon13%ptr)
  allocate(veg_pool%seed_bed%carbon14        , source=veg__mem%pool%seed_bed%carbon14%ptr)
  allocate(veg_pool%seed_bed%nitrogen15      , source=veg__mem%pool%seed_bed%nitrogen15%ptr)
#endif

#if defined(Get_veg_growth)
  allocate(veg_growth%leaf%carbon              , source=veg__mem%growth%leaf%carbon%ptr)
  allocate(veg_growth%leaf%nitrogen            , source=veg__mem%growth%leaf%nitrogen%ptr)
  allocate(veg_growth%leaf%phosphorus          , source=veg__mem%growth%leaf%phosphorus%ptr)
  allocate(veg_growth%leaf%carbon13            , source=veg__mem%growth%leaf%carbon13%ptr)
  allocate(veg_growth%leaf%carbon14            , source=veg__mem%growth%leaf%carbon14%ptr)
  allocate(veg_growth%leaf%nitrogen15          , source=veg__mem%growth%leaf%nitrogen15%ptr)
  allocate(veg_growth%fine_root%carbon         , source=veg__mem%growth%fine_root%carbon%ptr)
  allocate(veg_growth%fine_root%nitrogen       , source=veg__mem%growth%fine_root%nitrogen%ptr)
  allocate(veg_growth%fine_root%phosphorus     , source=veg__mem%growth%fine_root%phosphorus%ptr)
  allocate(veg_growth%fine_root%carbon13       , source=veg__mem%growth%fine_root%carbon13%ptr)
  allocate(veg_growth%fine_root%carbon14       , source=veg__mem%growth%fine_root%carbon14%ptr)
  allocate(veg_growth%fine_root%nitrogen15     , source=veg__mem%growth%fine_root%nitrogen15%ptr)
  allocate(veg_growth%coarse_root%carbon       , source=veg__mem%growth%coarse_root%carbon%ptr)
  allocate(veg_growth%coarse_root%nitrogen     , source=veg__mem%growth%coarse_root%nitrogen%ptr)
  allocate(veg_growth%coarse_root%phosphorus   , source=veg__mem%growth%coarse_root%phosphorus%ptr)
  allocate(veg_growth%coarse_root%carbon13     , source=veg__mem%growth%coarse_root%carbon13%ptr)
  allocate(veg_growth%coarse_root%carbon14     , source=veg__mem%growth%coarse_root%carbon14%ptr)
  allocate(veg_growth%coarse_root%nitrogen15   , source=veg__mem%growth%coarse_root%nitrogen15%ptr)
  allocate(veg_growth%sap_wood%carbon          , source=veg__mem%growth%sap_wood%carbon%ptr)
  allocate(veg_growth%sap_wood%nitrogen        , source=veg__mem%growth%sap_wood%nitrogen%ptr)
  allocate(veg_growth%sap_wood%phosphorus      , source=veg__mem%growth%sap_wood%phosphorus%ptr)
  allocate(veg_growth%sap_wood%carbon13        , source=veg__mem%growth%sap_wood%carbon13%ptr)
  allocate(veg_growth%sap_wood%carbon14        , source=veg__mem%growth%sap_wood%carbon14%ptr)
  allocate(veg_growth%sap_wood%nitrogen15      , source=veg__mem%growth%sap_wood%nitrogen15%ptr)
  allocate(veg_growth%heart_wood%carbon        , source=veg__mem%growth%heart_wood%carbon%ptr)
  allocate(veg_growth%heart_wood%nitrogen      , source=veg__mem%growth%heart_wood%nitrogen%ptr)
  allocate(veg_growth%heart_wood%phosphorus    , source=veg__mem%growth%heart_wood%phosphorus%ptr)
  allocate(veg_growth%heart_wood%carbon13      , source=veg__mem%growth%heart_wood%carbon13%ptr)
  allocate(veg_growth%heart_wood%carbon14      , source=veg__mem%growth%heart_wood%carbon14%ptr)
  allocate(veg_growth%heart_wood%nitrogen15    , source=veg__mem%growth%heart_wood%nitrogen15%ptr)
  allocate(veg_growth%labile%carbon            , source=veg__mem%growth%labile%carbon%ptr)
  allocate(veg_growth%labile%nitrogen          , source=veg__mem%growth%labile%nitrogen%ptr)
  allocate(veg_growth%labile%phosphorus        , source=veg__mem%growth%labile%phosphorus%ptr)
  allocate(veg_growth%labile%carbon13          , source=veg__mem%growth%labile%carbon13%ptr)
  allocate(veg_growth%labile%carbon14          , source=veg__mem%growth%labile%carbon14%ptr)
  allocate(veg_growth%labile%nitrogen15        , source=veg__mem%growth%labile%nitrogen15%ptr)
  allocate(veg_growth%reserve%carbon           , source=veg__mem%growth%reserve%carbon%ptr)
  allocate(veg_growth%reserve%nitrogen         , source=veg__mem%growth%reserve%nitrogen%ptr)
  allocate(veg_growth%reserve%phosphorus       , source=veg__mem%growth%reserve%phosphorus%ptr)
  allocate(veg_growth%reserve%carbon13         , source=veg__mem%growth%reserve%carbon13%ptr)
  allocate(veg_growth%reserve%carbon14         , source=veg__mem%growth%reserve%carbon14%ptr)
  allocate(veg_growth%reserve%nitrogen15       , source=veg__mem%growth%reserve%nitrogen15%ptr)
  allocate(veg_growth%fruit%carbon             , source=veg__mem%growth%fruit%carbon%ptr)
  allocate(veg_growth%fruit%nitrogen           , source=veg__mem%growth%fruit%nitrogen%ptr)
  allocate(veg_growth%fruit%phosphorus         , source=veg__mem%growth%fruit%phosphorus%ptr)
  allocate(veg_growth%fruit%carbon13           , source=veg__mem%growth%fruit%carbon13%ptr)
  allocate(veg_growth%fruit%carbon14           , source=veg__mem%growth%fruit%carbon14%ptr)
  allocate(veg_growth%fruit%nitrogen15         , source=veg__mem%growth%fruit%nitrogen15%ptr)
  allocate(veg_growth%seed_bed%carbon          , source=veg__mem%growth%seed_bed%carbon%ptr)
  allocate(veg_growth%seed_bed%nitrogen        , source=veg__mem%growth%seed_bed%nitrogen%ptr)
  allocate(veg_growth%seed_bed%phosphorus      , source=veg__mem%growth%seed_bed%phosphorus%ptr)
  allocate(veg_growth%seed_bed%carbon13        , source=veg__mem%growth%seed_bed%carbon13%ptr)
  allocate(veg_growth%seed_bed%carbon14        , source=veg__mem%growth%seed_bed%carbon14%ptr)
  allocate(veg_growth%seed_bed%nitrogen15      , source=veg__mem%growth%seed_bed%nitrogen15%ptr)
#endif

#if defined(Get_veg_litterfall)
  allocate(veg_litterfall%leaf%carbon              , source=veg__mem%litterfall%leaf%carbon%ptr)
  allocate(veg_litterfall%leaf%nitrogen            , source=veg__mem%litterfall%leaf%nitrogen%ptr)
  allocate(veg_litterfall%leaf%phosphorus          , source=veg__mem%litterfall%leaf%phosphorus%ptr)
  allocate(veg_litterfall%leaf%carbon13            , source=veg__mem%litterfall%leaf%carbon13%ptr)
  allocate(veg_litterfall%leaf%carbon14            , source=veg__mem%litterfall%leaf%carbon14%ptr)
  allocate(veg_litterfall%leaf%nitrogen15          , source=veg__mem%litterfall%leaf%nitrogen15%ptr)
  allocate(veg_litterfall%fine_root%carbon         , source=veg__mem%litterfall%fine_root%carbon%ptr)
  allocate(veg_litterfall%fine_root%nitrogen       , source=veg__mem%litterfall%fine_root%nitrogen%ptr)
  allocate(veg_litterfall%fine_root%phosphorus     , source=veg__mem%litterfall%fine_root%phosphorus%ptr)
  allocate(veg_litterfall%fine_root%carbon13       , source=veg__mem%litterfall%fine_root%carbon13%ptr)
  allocate(veg_litterfall%fine_root%carbon14       , source=veg__mem%litterfall%fine_root%carbon14%ptr)
  allocate(veg_litterfall%fine_root%nitrogen15     , source=veg__mem%litterfall%fine_root%nitrogen15%ptr)
  allocate(veg_litterfall%coarse_root%carbon       , source=veg__mem%litterfall%coarse_root%carbon%ptr)
  allocate(veg_litterfall%coarse_root%nitrogen     , source=veg__mem%litterfall%coarse_root%nitrogen%ptr)
  allocate(veg_litterfall%coarse_root%phosphorus   , source=veg__mem%litterfall%coarse_root%phosphorus%ptr)
  allocate(veg_litterfall%coarse_root%carbon13     , source=veg__mem%litterfall%coarse_root%carbon13%ptr)
  allocate(veg_litterfall%coarse_root%carbon14     , source=veg__mem%litterfall%coarse_root%carbon14%ptr)
  allocate(veg_litterfall%coarse_root%nitrogen15   , source=veg__mem%litterfall%coarse_root%nitrogen15%ptr)
  allocate(veg_litterfall%sap_wood%carbon          , source=veg__mem%litterfall%sap_wood%carbon%ptr)
  allocate(veg_litterfall%sap_wood%nitrogen        , source=veg__mem%litterfall%sap_wood%nitrogen%ptr)
  allocate(veg_litterfall%sap_wood%phosphorus      , source=veg__mem%litterfall%sap_wood%phosphorus%ptr)
  allocate(veg_litterfall%sap_wood%carbon13        , source=veg__mem%litterfall%sap_wood%carbon13%ptr)
  allocate(veg_litterfall%sap_wood%carbon14        , source=veg__mem%litterfall%sap_wood%carbon14%ptr)
  allocate(veg_litterfall%sap_wood%nitrogen15      , source=veg__mem%litterfall%sap_wood%nitrogen15%ptr)
  allocate(veg_litterfall%heart_wood%carbon        , source=veg__mem%litterfall%heart_wood%carbon%ptr)
  allocate(veg_litterfall%heart_wood%nitrogen      , source=veg__mem%litterfall%heart_wood%nitrogen%ptr)
  allocate(veg_litterfall%heart_wood%phosphorus    , source=veg__mem%litterfall%heart_wood%phosphorus%ptr)
  allocate(veg_litterfall%heart_wood%carbon13      , source=veg__mem%litterfall%heart_wood%carbon13%ptr)
  allocate(veg_litterfall%heart_wood%carbon14      , source=veg__mem%litterfall%heart_wood%carbon14%ptr)
  allocate(veg_litterfall%heart_wood%nitrogen15    , source=veg__mem%litterfall%heart_wood%nitrogen15%ptr)
  allocate(veg_litterfall%labile%carbon            , source=veg__mem%litterfall%labile%carbon%ptr)
  allocate(veg_litterfall%labile%nitrogen          , source=veg__mem%litterfall%labile%nitrogen%ptr)
  allocate(veg_litterfall%labile%phosphorus        , source=veg__mem%litterfall%labile%phosphorus%ptr)
  allocate(veg_litterfall%labile%carbon13          , source=veg__mem%litterfall%labile%carbon13%ptr)
  allocate(veg_litterfall%labile%carbon14          , source=veg__mem%litterfall%labile%carbon14%ptr)
  allocate(veg_litterfall%labile%nitrogen15        , source=veg__mem%litterfall%labile%nitrogen15%ptr)
  allocate(veg_litterfall%reserve%carbon           , source=veg__mem%litterfall%reserve%carbon%ptr)
  allocate(veg_litterfall%reserve%nitrogen         , source=veg__mem%litterfall%reserve%nitrogen%ptr)
  allocate(veg_litterfall%reserve%phosphorus       , source=veg__mem%litterfall%reserve%phosphorus%ptr)
  allocate(veg_litterfall%reserve%carbon13         , source=veg__mem%litterfall%reserve%carbon13%ptr)
  allocate(veg_litterfall%reserve%carbon14         , source=veg__mem%litterfall%reserve%carbon14%ptr)
  allocate(veg_litterfall%reserve%nitrogen15       , source=veg__mem%litterfall%reserve%nitrogen15%ptr)
  allocate(veg_litterfall%fruit%carbon             , source=veg__mem%litterfall%fruit%carbon%ptr)
  allocate(veg_litterfall%fruit%nitrogen           , source=veg__mem%litterfall%fruit%nitrogen%ptr)
  allocate(veg_litterfall%fruit%phosphorus         , source=veg__mem%litterfall%fruit%phosphorus%ptr)
  allocate(veg_litterfall%fruit%carbon13           , source=veg__mem%litterfall%fruit%carbon13%ptr)
  allocate(veg_litterfall%fruit%carbon14           , source=veg__mem%litterfall%fruit%carbon14%ptr)
  allocate(veg_litterfall%fruit%nitrogen15         , source=veg__mem%litterfall%fruit%nitrogen15%ptr)
  allocate(veg_litterfall%seed_bed%carbon          , source=veg__mem%litterfall%seed_bed%carbon%ptr)
  allocate(veg_litterfall%seed_bed%nitrogen        , source=veg__mem%litterfall%seed_bed%nitrogen%ptr)
  allocate(veg_litterfall%seed_bed%phosphorus      , source=veg__mem%litterfall%seed_bed%phosphorus%ptr)
  allocate(veg_litterfall%seed_bed%carbon13        , source=veg__mem%litterfall%seed_bed%carbon13%ptr)
  allocate(veg_litterfall%seed_bed%carbon14        , source=veg__mem%litterfall%seed_bed%carbon14%ptr)
  allocate(veg_litterfall%seed_bed%nitrogen15      , source=veg__mem%litterfall%seed_bed%nitrogen15%ptr)
#endif

#if defined(Get_veg_falloc)
  allocate(veg_falloc%leaf%carbon              , source=veg__mem%falloc%leaf%carbon%ptr)
  allocate(veg_falloc%leaf%nitrogen            , source=veg__mem%falloc%leaf%nitrogen%ptr)
  allocate(veg_falloc%leaf%phosphorus          , source=veg__mem%falloc%leaf%phosphorus%ptr)
  allocate(veg_falloc%leaf%carbon13            , source=veg__mem%falloc%leaf%carbon13%ptr)
  allocate(veg_falloc%leaf%carbon14            , source=veg__mem%falloc%leaf%carbon14%ptr)
  allocate(veg_falloc%leaf%nitrogen15          , source=veg__mem%falloc%leaf%nitrogen15%ptr)
  allocate(veg_falloc%fine_root%carbon         , source=veg__mem%falloc%fine_root%carbon%ptr)
  allocate(veg_falloc%fine_root%nitrogen       , source=veg__mem%falloc%fine_root%nitrogen%ptr)
  allocate(veg_falloc%fine_root%phosphorus     , source=veg__mem%falloc%fine_root%phosphorus%ptr)
  allocate(veg_falloc%fine_root%carbon13       , source=veg__mem%falloc%fine_root%carbon13%ptr)
  allocate(veg_falloc%fine_root%carbon14       , source=veg__mem%falloc%fine_root%carbon14%ptr)
  allocate(veg_falloc%fine_root%nitrogen15     , source=veg__mem%falloc%fine_root%nitrogen15%ptr)
  allocate(veg_falloc%coarse_root%carbon       , source=veg__mem%falloc%coarse_root%carbon%ptr)
  allocate(veg_falloc%coarse_root%nitrogen     , source=veg__mem%falloc%coarse_root%nitrogen%ptr)
  allocate(veg_falloc%coarse_root%phosphorus   , source=veg__mem%falloc%coarse_root%phosphorus%ptr)
  allocate(veg_falloc%coarse_root%carbon13     , source=veg__mem%falloc%coarse_root%carbon13%ptr)
  allocate(veg_falloc%coarse_root%carbon14     , source=veg__mem%falloc%coarse_root%carbon14%ptr)
  allocate(veg_falloc%coarse_root%nitrogen15   , source=veg__mem%falloc%coarse_root%nitrogen15%ptr)
  allocate(veg_falloc%sap_wood%carbon          , source=veg__mem%falloc%sap_wood%carbon%ptr)
  allocate(veg_falloc%sap_wood%nitrogen        , source=veg__mem%falloc%sap_wood%nitrogen%ptr)
  allocate(veg_falloc%sap_wood%phosphorus      , source=veg__mem%falloc%sap_wood%phosphorus%ptr)
  allocate(veg_falloc%sap_wood%carbon13        , source=veg__mem%falloc%sap_wood%carbon13%ptr)
  allocate(veg_falloc%sap_wood%carbon14        , source=veg__mem%falloc%sap_wood%carbon14%ptr)
  allocate(veg_falloc%sap_wood%nitrogen15      , source=veg__mem%falloc%sap_wood%nitrogen15%ptr)
  allocate(veg_falloc%heart_wood%carbon        , source=veg__mem%falloc%heart_wood%carbon%ptr)
  allocate(veg_falloc%heart_wood%nitrogen      , source=veg__mem%falloc%heart_wood%nitrogen%ptr)
  allocate(veg_falloc%heart_wood%phosphorus    , source=veg__mem%falloc%heart_wood%phosphorus%ptr)
  allocate(veg_falloc%heart_wood%carbon13      , source=veg__mem%falloc%heart_wood%carbon13%ptr)
  allocate(veg_falloc%heart_wood%carbon14      , source=veg__mem%falloc%heart_wood%carbon14%ptr)
  allocate(veg_falloc%heart_wood%nitrogen15    , source=veg__mem%falloc%heart_wood%nitrogen15%ptr)
  allocate(veg_falloc%labile%carbon            , source=veg__mem%falloc%labile%carbon%ptr)
  allocate(veg_falloc%labile%nitrogen          , source=veg__mem%falloc%labile%nitrogen%ptr)
  allocate(veg_falloc%labile%phosphorus        , source=veg__mem%falloc%labile%phosphorus%ptr)
  allocate(veg_falloc%labile%carbon13          , source=veg__mem%falloc%labile%carbon13%ptr)
  allocate(veg_falloc%labile%carbon14          , source=veg__mem%falloc%labile%carbon14%ptr)
  allocate(veg_falloc%labile%nitrogen15        , source=veg__mem%falloc%labile%nitrogen15%ptr)
  allocate(veg_falloc%reserve%carbon           , source=veg__mem%falloc%reserve%carbon%ptr)
  allocate(veg_falloc%reserve%nitrogen         , source=veg__mem%falloc%reserve%nitrogen%ptr)
  allocate(veg_falloc%reserve%phosphorus       , source=veg__mem%falloc%reserve%phosphorus%ptr)
  allocate(veg_falloc%reserve%carbon13         , source=veg__mem%falloc%reserve%carbon13%ptr)
  allocate(veg_falloc%reserve%carbon14         , source=veg__mem%falloc%reserve%carbon14%ptr)
  allocate(veg_falloc%reserve%nitrogen15       , source=veg__mem%falloc%reserve%nitrogen15%ptr)
  allocate(veg_falloc%fruit%carbon             , source=veg__mem%falloc%fruit%carbon%ptr)
  allocate(veg_falloc%fruit%nitrogen           , source=veg__mem%falloc%fruit%nitrogen%ptr)
  allocate(veg_falloc%fruit%phosphorus         , source=veg__mem%falloc%fruit%phosphorus%ptr)
  allocate(veg_falloc%fruit%carbon13           , source=veg__mem%falloc%fruit%carbon13%ptr)
  allocate(veg_falloc%fruit%carbon14           , source=veg__mem%falloc%fruit%carbon14%ptr)
  allocate(veg_falloc%fruit%nitrogen15         , source=veg__mem%falloc%fruit%nitrogen15%ptr)
  allocate(veg_falloc%seed_bed%carbon          , source=veg__mem%falloc%seed_bed%carbon%ptr)
  allocate(veg_falloc%seed_bed%nitrogen        , source=veg__mem%falloc%seed_bed%nitrogen%ptr)
  allocate(veg_falloc%seed_bed%phosphorus      , source=veg__mem%falloc%seed_bed%phosphorus%ptr)
  allocate(veg_falloc%seed_bed%carbon13        , source=veg__mem%falloc%seed_bed%carbon13%ptr)
  allocate(veg_falloc%seed_bed%carbon14        , source=veg__mem%falloc%seed_bed%carbon14%ptr)
  allocate(veg_falloc%seed_bed%nitrogen15      , source=veg__mem%falloc%seed_bed%nitrogen15%ptr)
#endif

#if defined(Get_veg_total_flux)
  allocate(veg_total_flux%leaf%carbon              , source=veg__mem%total_flux%leaf%carbon%ptr)
  allocate(veg_total_flux%leaf%nitrogen            , source=veg__mem%total_flux%leaf%nitrogen%ptr)
  allocate(veg_total_flux%leaf%phosphorus          , source=veg__mem%total_flux%leaf%phosphorus%ptr)
  allocate(veg_total_flux%leaf%carbon13            , source=veg__mem%total_flux%leaf%carbon13%ptr)
  allocate(veg_total_flux%leaf%carbon14            , source=veg__mem%total_flux%leaf%carbon14%ptr)
  allocate(veg_total_flux%leaf%nitrogen15          , source=veg__mem%total_flux%leaf%nitrogen15%ptr)
  allocate(veg_total_flux%fine_root%carbon         , source=veg__mem%total_flux%fine_root%carbon%ptr)
  allocate(veg_total_flux%fine_root%nitrogen       , source=veg__mem%total_flux%fine_root%nitrogen%ptr)
  allocate(veg_total_flux%fine_root%phosphorus     , source=veg__mem%total_flux%fine_root%phosphorus%ptr)
  allocate(veg_total_flux%fine_root%carbon13       , source=veg__mem%total_flux%fine_root%carbon13%ptr)
  allocate(veg_total_flux%fine_root%carbon14       , source=veg__mem%total_flux%fine_root%carbon14%ptr)
  allocate(veg_total_flux%fine_root%nitrogen15     , source=veg__mem%total_flux%fine_root%nitrogen15%ptr)
  allocate(veg_total_flux%coarse_root%carbon       , source=veg__mem%total_flux%coarse_root%carbon%ptr)
  allocate(veg_total_flux%coarse_root%nitrogen     , source=veg__mem%total_flux%coarse_root%nitrogen%ptr)
  allocate(veg_total_flux%coarse_root%phosphorus   , source=veg__mem%total_flux%coarse_root%phosphorus%ptr)
  allocate(veg_total_flux%coarse_root%carbon13     , source=veg__mem%total_flux%coarse_root%carbon13%ptr)
  allocate(veg_total_flux%coarse_root%carbon14     , source=veg__mem%total_flux%coarse_root%carbon14%ptr)
  allocate(veg_total_flux%coarse_root%nitrogen15   , source=veg__mem%total_flux%coarse_root%nitrogen15%ptr)
  allocate(veg_total_flux%sap_wood%carbon          , source=veg__mem%total_flux%sap_wood%carbon%ptr)
  allocate(veg_total_flux%sap_wood%nitrogen        , source=veg__mem%total_flux%sap_wood%nitrogen%ptr)
  allocate(veg_total_flux%sap_wood%phosphorus      , source=veg__mem%total_flux%sap_wood%phosphorus%ptr)
  allocate(veg_total_flux%sap_wood%carbon13        , source=veg__mem%total_flux%sap_wood%carbon13%ptr)
  allocate(veg_total_flux%sap_wood%carbon14        , source=veg__mem%total_flux%sap_wood%carbon14%ptr)
  allocate(veg_total_flux%sap_wood%nitrogen15      , source=veg__mem%total_flux%sap_wood%nitrogen15%ptr)
  allocate(veg_total_flux%heart_wood%carbon        , source=veg__mem%total_flux%heart_wood%carbon%ptr)
  allocate(veg_total_flux%heart_wood%nitrogen      , source=veg__mem%total_flux%heart_wood%nitrogen%ptr)
  allocate(veg_total_flux%heart_wood%phosphorus    , source=veg__mem%total_flux%heart_wood%phosphorus%ptr)
  allocate(veg_total_flux%heart_wood%carbon13      , source=veg__mem%total_flux%heart_wood%carbon13%ptr)
  allocate(veg_total_flux%heart_wood%carbon14      , source=veg__mem%total_flux%heart_wood%carbon14%ptr)
  allocate(veg_total_flux%heart_wood%nitrogen15    , source=veg__mem%total_flux%heart_wood%nitrogen15%ptr)
  allocate(veg_total_flux%labile%carbon            , source=veg__mem%total_flux%labile%carbon%ptr)
  allocate(veg_total_flux%labile%nitrogen          , source=veg__mem%total_flux%labile%nitrogen%ptr)
  allocate(veg_total_flux%labile%phosphorus        , source=veg__mem%total_flux%labile%phosphorus%ptr)
  allocate(veg_total_flux%labile%carbon13          , source=veg__mem%total_flux%labile%carbon13%ptr)
  allocate(veg_total_flux%labile%carbon14          , source=veg__mem%total_flux%labile%carbon14%ptr)
  allocate(veg_total_flux%labile%nitrogen15        , source=veg__mem%total_flux%labile%nitrogen15%ptr)
  allocate(veg_total_flux%reserve%carbon           , source=veg__mem%total_flux%reserve%carbon%ptr)
  allocate(veg_total_flux%reserve%nitrogen         , source=veg__mem%total_flux%reserve%nitrogen%ptr)
  allocate(veg_total_flux%reserve%phosphorus       , source=veg__mem%total_flux%reserve%phosphorus%ptr)
  allocate(veg_total_flux%reserve%carbon13         , source=veg__mem%total_flux%reserve%carbon13%ptr)
  allocate(veg_total_flux%reserve%carbon14         , source=veg__mem%total_flux%reserve%carbon14%ptr)
  allocate(veg_total_flux%reserve%nitrogen15       , source=veg__mem%total_flux%reserve%nitrogen15%ptr)
  allocate(veg_total_flux%fruit%carbon             , source=veg__mem%total_flux%fruit%carbon%ptr)
  allocate(veg_total_flux%fruit%nitrogen           , source=veg__mem%total_flux%fruit%nitrogen%ptr)
  allocate(veg_total_flux%fruit%phosphorus         , source=veg__mem%total_flux%fruit%phosphorus%ptr)
  allocate(veg_total_flux%fruit%carbon13           , source=veg__mem%total_flux%fruit%carbon13%ptr)
  allocate(veg_total_flux%fruit%carbon14           , source=veg__mem%total_flux%fruit%carbon14%ptr)
  allocate(veg_total_flux%fruit%nitrogen15         , source=veg__mem%total_flux%fruit%nitrogen15%ptr)
  allocate(veg_total_flux%seed_bed%carbon          , source=veg__mem%total_flux%seed_bed%carbon%ptr)
  allocate(veg_total_flux%seed_bed%nitrogen        , source=veg__mem%total_flux%seed_bed%nitrogen%ptr)
  allocate(veg_total_flux%seed_bed%phosphorus      , source=veg__mem%total_flux%seed_bed%phosphorus%ptr)
  allocate(veg_total_flux%seed_bed%carbon13        , source=veg__mem%total_flux%seed_bed%carbon13%ptr)
  allocate(veg_total_flux%seed_bed%carbon14        , source=veg__mem%total_flux%seed_bed%carbon14%ptr)
  allocate(veg_total_flux%seed_bed%nitrogen15      , source=veg__mem%total_flux%seed_bed%nitrogen15%ptr)
#endif

#ifdef Get_veg_exudation
  allocate(veg_exudation%carbon       , source=veg__mem%exudation%carbon%ptr)
  allocate(veg_exudation%nitrogen     , source=veg__mem%exudation%nitrogen%ptr)
  allocate(veg_exudation%phosphorus   , source=veg__mem%exudation%phosphorus%ptr)
  allocate(veg_exudation%carbon13     , source=veg__mem%exudation%carbon13%ptr)
  allocate(veg_exudation%carbon14     , source=veg__mem%exudation%carbon14%ptr)
  allocate(veg_exudation%nitrogen15   , source=veg__mem%exudation%nitrogen15%ptr)
#endif

#ifdef Get_veg_establishment
  allocate(veg_establishment%carbon       , source=veg__mem%establishment%carbon%ptr)
  allocate(veg_establishment%nitrogen     , source=veg__mem%establishment%nitrogen%ptr)
  allocate(veg_establishment%phosphorus   , source=veg__mem%establishment%phosphorus%ptr)
  allocate(veg_establishment%carbon13     , source=veg__mem%establishment%carbon13%ptr)
  allocate(veg_establishment%carbon14     , source=veg__mem%establishment%carbon14%ptr)
  allocate(veg_establishment%nitrogen15   , source=veg__mem%establishment%nitrogen15%ptr)
#endif

#ifdef Get_veg_reserve_use
  allocate(veg_reserve_use%carbon       , source=veg__mem%reserve_use%carbon%ptr)
  allocate(veg_reserve_use%nitrogen     , source=veg__mem%reserve_use%nitrogen%ptr)
  allocate(veg_reserve_use%phosphorus   , source=veg__mem%reserve_use%phosphorus%ptr)
  allocate(veg_reserve_use%carbon13     , source=veg__mem%reserve_use%carbon13%ptr)
  allocate(veg_reserve_use%carbon14     , source=veg__mem%reserve_use%carbon14%ptr)
  allocate(veg_reserve_use%nitrogen15   , source=veg__mem%reserve_use%nitrogen15%ptr)
#endif

#ifdef Get_veg_total_biomass
  allocate(veg_total_biomass%carbon       , source=veg__mem%total_biomass%carbon%ptr)
  allocate(veg_total_biomass%nitrogen     , source=veg__mem%total_biomass%nitrogen%ptr)
  allocate(veg_total_biomass%phosphorus   , source=veg__mem%total_biomass%phosphorus%ptr)
  allocate(veg_total_biomass%carbon13     , source=veg__mem%total_biomass%carbon13%ptr)
  allocate(veg_total_biomass%carbon14     , source=veg__mem%total_biomass%carbon14%ptr)
  allocate(veg_total_biomass%nitrogen15   , source=veg__mem%total_biomass%nitrogen15%ptr)
#endif

! -------------------------------------------------------
! deallocate helper variables for vegetation pools
! and pass values to pools of vegetation memory
! -------------------------------------------------------

#if defined(Dealloc_veg_pool)
  ! pass the values (back) to process memory
  veg__mem%pool%leaf%carbon%ptr             = veg_pool%leaf%carbon
  veg__mem%pool%leaf%nitrogen%ptr           = veg_pool%leaf%nitrogen
  veg__mem%pool%leaf%phosphorus%ptr         = veg_pool%leaf%phosphorus
  veg__mem%pool%leaf%carbon13%ptr           = veg_pool%leaf%carbon13
  veg__mem%pool%leaf%carbon14%ptr           = veg_pool%leaf%carbon14
  veg__mem%pool%leaf%nitrogen15%ptr         = veg_pool%leaf%nitrogen15
  veg__mem%pool%fine_root%carbon%ptr        = veg_pool%fine_root%carbon
  veg__mem%pool%fine_root%nitrogen%ptr      = veg_pool%fine_root%nitrogen
  veg__mem%pool%fine_root%phosphorus%ptr    = veg_pool%fine_root%phosphorus
  veg__mem%pool%fine_root%carbon13%ptr      = veg_pool%fine_root%carbon13
  veg__mem%pool%fine_root%carbon14%ptr      = veg_pool%fine_root%carbon14
  veg__mem%pool%fine_root%nitrogen15%ptr    = veg_pool%fine_root%nitrogen15
  veg__mem%pool%coarse_root%carbon%ptr      = veg_pool%coarse_root%carbon
  veg__mem%pool%coarse_root%nitrogen%ptr    = veg_pool%coarse_root%nitrogen
  veg__mem%pool%coarse_root%phosphorus%ptr  = veg_pool%coarse_root%phosphorus
  veg__mem%pool%coarse_root%carbon13%ptr    = veg_pool%coarse_root%carbon13
  veg__mem%pool%coarse_root%carbon14%ptr    = veg_pool%coarse_root%carbon14
  veg__mem%pool%coarse_root%nitrogen15%ptr  = veg_pool%coarse_root%nitrogen15
  veg__mem%pool%sap_wood%carbon%ptr         = veg_pool%sap_wood%carbon
  veg__mem%pool%sap_wood%nitrogen%ptr       = veg_pool%sap_wood%nitrogen
  veg__mem%pool%sap_wood%phosphorus%ptr     = veg_pool%sap_wood%phosphorus
  veg__mem%pool%sap_wood%carbon13%ptr       = veg_pool%sap_wood%carbon13
  veg__mem%pool%sap_wood%carbon14%ptr       = veg_pool%sap_wood%carbon14
  veg__mem%pool%sap_wood%nitrogen15%ptr     = veg_pool%sap_wood%nitrogen15
  veg__mem%pool%heart_wood%carbon%ptr       = veg_pool%heart_wood%carbon
  veg__mem%pool%heart_wood%nitrogen%ptr     = veg_pool%heart_wood%nitrogen
  veg__mem%pool%heart_wood%phosphorus%ptr   = veg_pool%heart_wood%phosphorus
  veg__mem%pool%heart_wood%carbon13%ptr     = veg_pool%heart_wood%carbon13
  veg__mem%pool%heart_wood%carbon14%ptr     = veg_pool%heart_wood%carbon14
  veg__mem%pool%heart_wood%nitrogen15%ptr   = veg_pool%heart_wood%nitrogen15
  veg__mem%pool%labile%carbon%ptr           = veg_pool%labile%carbon
  veg__mem%pool%labile%nitrogen%ptr         = veg_pool%labile%nitrogen
  veg__mem%pool%labile%phosphorus%ptr       = veg_pool%labile%phosphorus
  veg__mem%pool%labile%carbon13%ptr         = veg_pool%labile%carbon13
  veg__mem%pool%labile%carbon14%ptr         = veg_pool%labile%carbon14
  veg__mem%pool%labile%nitrogen15%ptr       = veg_pool%labile%nitrogen15
  veg__mem%pool%reserve%carbon%ptr          = veg_pool%reserve%carbon
  veg__mem%pool%reserve%nitrogen%ptr        = veg_pool%reserve%nitrogen
  veg__mem%pool%reserve%phosphorus%ptr      = veg_pool%reserve%phosphorus
  veg__mem%pool%reserve%carbon13%ptr        = veg_pool%reserve%carbon13
  veg__mem%pool%reserve%carbon14%ptr        = veg_pool%reserve%carbon14
  veg__mem%pool%reserve%nitrogen15%ptr      = veg_pool%reserve%nitrogen15
  veg__mem%pool%fruit%carbon%ptr            = veg_pool%fruit%carbon
  veg__mem%pool%fruit%nitrogen%ptr          = veg_pool%fruit%nitrogen
  veg__mem%pool%fruit%phosphorus%ptr        = veg_pool%fruit%phosphorus
  veg__mem%pool%fruit%carbon13%ptr          = veg_pool%fruit%carbon13
  veg__mem%pool%fruit%carbon14%ptr          = veg_pool%fruit%carbon14
  veg__mem%pool%fruit%nitrogen15%ptr        = veg_pool%fruit%nitrogen15
  veg__mem%pool%seed_bed%carbon%ptr         = veg_pool%seed_bed%carbon
  veg__mem%pool%seed_bed%nitrogen%ptr       = veg_pool%seed_bed%nitrogen
  veg__mem%pool%seed_bed%phosphorus%ptr     = veg_pool%seed_bed%phosphorus
  veg__mem%pool%seed_bed%carbon13%ptr       = veg_pool%seed_bed%carbon13
  veg__mem%pool%seed_bed%carbon14%ptr       = veg_pool%seed_bed%carbon14
  veg__mem%pool%seed_bed%nitrogen15%ptr     = veg_pool%seed_bed%nitrogen15
  ! deallocate variables in helper pool structures
  deallocate(veg_pool%leaf%carbon)
  deallocate(veg_pool%leaf%nitrogen)
  deallocate(veg_pool%leaf%phosphorus)
  deallocate(veg_pool%leaf%carbon13)
  deallocate(veg_pool%leaf%carbon14)
  deallocate(veg_pool%leaf%nitrogen15)
  deallocate(veg_pool%fine_root%carbon)
  deallocate(veg_pool%fine_root%nitrogen)
  deallocate(veg_pool%fine_root%phosphorus)
  deallocate(veg_pool%fine_root%carbon13)
  deallocate(veg_pool%fine_root%carbon14)
  deallocate(veg_pool%fine_root%nitrogen15)
  deallocate(veg_pool%coarse_root%carbon)
  deallocate(veg_pool%coarse_root%nitrogen)
  deallocate(veg_pool%coarse_root%phosphorus)
  deallocate(veg_pool%coarse_root%carbon13)
  deallocate(veg_pool%coarse_root%carbon14)
  deallocate(veg_pool%coarse_root%nitrogen15)
  deallocate(veg_pool%sap_wood%carbon)
  deallocate(veg_pool%sap_wood%nitrogen)
  deallocate(veg_pool%sap_wood%phosphorus)
  deallocate(veg_pool%sap_wood%carbon13)
  deallocate(veg_pool%sap_wood%carbon14)
  deallocate(veg_pool%sap_wood%nitrogen15)
  deallocate(veg_pool%heart_wood%carbon)
  deallocate(veg_pool%heart_wood%nitrogen)
  deallocate(veg_pool%heart_wood%phosphorus)
  deallocate(veg_pool%heart_wood%carbon13)
  deallocate(veg_pool%heart_wood%carbon14)
  deallocate(veg_pool%heart_wood%nitrogen15)
  deallocate(veg_pool%labile%carbon)
  deallocate(veg_pool%labile%nitrogen)
  deallocate(veg_pool%labile%phosphorus)
  deallocate(veg_pool%labile%carbon13)
  deallocate(veg_pool%labile%carbon14)
  deallocate(veg_pool%labile%nitrogen15)
  deallocate(veg_pool%reserve%carbon)
  deallocate(veg_pool%reserve%nitrogen)
  deallocate(veg_pool%reserve%phosphorus)
  deallocate(veg_pool%reserve%carbon13)
  deallocate(veg_pool%reserve%carbon14)
  deallocate(veg_pool%reserve%nitrogen15)
  deallocate(veg_pool%fruit%carbon)
  deallocate(veg_pool%fruit%nitrogen)
  deallocate(veg_pool%fruit%phosphorus)
  deallocate(veg_pool%fruit%carbon13)
  deallocate(veg_pool%fruit%carbon14)
  deallocate(veg_pool%fruit%nitrogen15)
  deallocate(veg_pool%seed_bed%carbon)
  deallocate(veg_pool%seed_bed%nitrogen)
  deallocate(veg_pool%seed_bed%phosphorus)
  deallocate(veg_pool%seed_bed%carbon13)
  deallocate(veg_pool%seed_bed%carbon14)
  deallocate(veg_pool%seed_bed%nitrogen15)
#endif

#if defined(Dealloc_veg_growth)
  ! pass the values (back) to process memory
  veg__mem%growth%leaf%carbon%ptr             = veg_growth%leaf%carbon
  veg__mem%growth%leaf%nitrogen%ptr           = veg_growth%leaf%nitrogen
  veg__mem%growth%leaf%phosphorus%ptr         = veg_growth%leaf%phosphorus
  veg__mem%growth%leaf%carbon13%ptr           = veg_growth%leaf%carbon13
  veg__mem%growth%leaf%carbon14%ptr           = veg_growth%leaf%carbon14
  veg__mem%growth%leaf%nitrogen15%ptr         = veg_growth%leaf%nitrogen15
  veg__mem%growth%fine_root%carbon%ptr        = veg_growth%fine_root%carbon
  veg__mem%growth%fine_root%nitrogen%ptr      = veg_growth%fine_root%nitrogen
  veg__mem%growth%fine_root%phosphorus%ptr    = veg_growth%fine_root%phosphorus
  veg__mem%growth%fine_root%carbon13%ptr      = veg_growth%fine_root%carbon13
  veg__mem%growth%fine_root%carbon14%ptr      = veg_growth%fine_root%carbon14
  veg__mem%growth%fine_root%nitrogen15%ptr    = veg_growth%fine_root%nitrogen15
  veg__mem%growth%coarse_root%carbon%ptr      = veg_growth%coarse_root%carbon
  veg__mem%growth%coarse_root%nitrogen%ptr    = veg_growth%coarse_root%nitrogen
  veg__mem%growth%coarse_root%phosphorus%ptr  = veg_growth%coarse_root%phosphorus
  veg__mem%growth%coarse_root%carbon13%ptr    = veg_growth%coarse_root%carbon13
  veg__mem%growth%coarse_root%carbon14%ptr    = veg_growth%coarse_root%carbon14
  veg__mem%growth%coarse_root%nitrogen15%ptr  = veg_growth%coarse_root%nitrogen15
  veg__mem%growth%sap_wood%carbon%ptr         = veg_growth%sap_wood%carbon
  veg__mem%growth%sap_wood%nitrogen%ptr       = veg_growth%sap_wood%nitrogen
  veg__mem%growth%sap_wood%phosphorus%ptr     = veg_growth%sap_wood%phosphorus
  veg__mem%growth%sap_wood%carbon13%ptr       = veg_growth%sap_wood%carbon13
  veg__mem%growth%sap_wood%carbon14%ptr       = veg_growth%sap_wood%carbon14
  veg__mem%growth%sap_wood%nitrogen15%ptr     = veg_growth%sap_wood%nitrogen15
  veg__mem%growth%heart_wood%carbon%ptr       = veg_growth%heart_wood%carbon
  veg__mem%growth%heart_wood%nitrogen%ptr     = veg_growth%heart_wood%nitrogen
  veg__mem%growth%heart_wood%phosphorus%ptr   = veg_growth%heart_wood%phosphorus
  veg__mem%growth%heart_wood%carbon13%ptr     = veg_growth%heart_wood%carbon13
  veg__mem%growth%heart_wood%carbon14%ptr     = veg_growth%heart_wood%carbon14
  veg__mem%growth%heart_wood%nitrogen15%ptr   = veg_growth%heart_wood%nitrogen15
  veg__mem%growth%labile%carbon%ptr           = veg_growth%labile%carbon
  veg__mem%growth%labile%nitrogen%ptr         = veg_growth%labile%nitrogen
  veg__mem%growth%labile%phosphorus%ptr       = veg_growth%labile%phosphorus
  veg__mem%growth%labile%carbon13%ptr         = veg_growth%labile%carbon13
  veg__mem%growth%labile%carbon14%ptr         = veg_growth%labile%carbon14
  veg__mem%growth%labile%nitrogen15%ptr       = veg_growth%labile%nitrogen15
  veg__mem%growth%reserve%carbon%ptr          = veg_growth%reserve%carbon
  veg__mem%growth%reserve%nitrogen%ptr        = veg_growth%reserve%nitrogen
  veg__mem%growth%reserve%phosphorus%ptr      = veg_growth%reserve%phosphorus
  veg__mem%growth%reserve%carbon13%ptr        = veg_growth%reserve%carbon13
  veg__mem%growth%reserve%carbon14%ptr        = veg_growth%reserve%carbon14
  veg__mem%growth%reserve%nitrogen15%ptr      = veg_growth%reserve%nitrogen15
  veg__mem%growth%fruit%carbon%ptr            = veg_growth%fruit%carbon
  veg__mem%growth%fruit%nitrogen%ptr          = veg_growth%fruit%nitrogen
  veg__mem%growth%fruit%phosphorus%ptr        = veg_growth%fruit%phosphorus
  veg__mem%growth%fruit%carbon13%ptr          = veg_growth%fruit%carbon13
  veg__mem%growth%fruit%carbon14%ptr          = veg_growth%fruit%carbon14
  veg__mem%growth%fruit%nitrogen15%ptr        = veg_growth%fruit%nitrogen15
  veg__mem%growth%seed_bed%carbon%ptr         = veg_growth%seed_bed%carbon
  veg__mem%growth%seed_bed%nitrogen%ptr       = veg_growth%seed_bed%nitrogen
  veg__mem%growth%seed_bed%phosphorus%ptr     = veg_growth%seed_bed%phosphorus
  veg__mem%growth%seed_bed%carbon13%ptr       = veg_growth%seed_bed%carbon13
  veg__mem%growth%seed_bed%carbon14%ptr       = veg_growth%seed_bed%carbon14
  veg__mem%growth%seed_bed%nitrogen15%ptr     = veg_growth%seed_bed%nitrogen15
  ! deallocate variables in helper pool structures
  deallocate(veg_growth%leaf%carbon)
  deallocate(veg_growth%leaf%nitrogen)
  deallocate(veg_growth%leaf%phosphorus)
  deallocate(veg_growth%leaf%carbon13)
  deallocate(veg_growth%leaf%carbon14)
  deallocate(veg_growth%leaf%nitrogen15)
  deallocate(veg_growth%fine_root%carbon)
  deallocate(veg_growth%fine_root%nitrogen)
  deallocate(veg_growth%fine_root%phosphorus)
  deallocate(veg_growth%fine_root%carbon13)
  deallocate(veg_growth%fine_root%carbon14)
  deallocate(veg_growth%fine_root%nitrogen15)
  deallocate(veg_growth%coarse_root%carbon)
  deallocate(veg_growth%coarse_root%nitrogen)
  deallocate(veg_growth%coarse_root%phosphorus)
  deallocate(veg_growth%coarse_root%carbon13)
  deallocate(veg_growth%coarse_root%carbon14)
  deallocate(veg_growth%coarse_root%nitrogen15)
  deallocate(veg_growth%sap_wood%carbon)
  deallocate(veg_growth%sap_wood%nitrogen)
  deallocate(veg_growth%sap_wood%phosphorus)
  deallocate(veg_growth%sap_wood%carbon13)
  deallocate(veg_growth%sap_wood%carbon14)
  deallocate(veg_growth%sap_wood%nitrogen15)
  deallocate(veg_growth%heart_wood%carbon)
  deallocate(veg_growth%heart_wood%nitrogen)
  deallocate(veg_growth%heart_wood%phosphorus)
  deallocate(veg_growth%heart_wood%carbon13)
  deallocate(veg_growth%heart_wood%carbon14)
  deallocate(veg_growth%heart_wood%nitrogen15)
  deallocate(veg_growth%labile%carbon)
  deallocate(veg_growth%labile%nitrogen)
  deallocate(veg_growth%labile%phosphorus)
  deallocate(veg_growth%labile%carbon13)
  deallocate(veg_growth%labile%carbon14)
  deallocate(veg_growth%labile%nitrogen15)
  deallocate(veg_growth%reserve%carbon)
  deallocate(veg_growth%reserve%nitrogen)
  deallocate(veg_growth%reserve%phosphorus)
  deallocate(veg_growth%reserve%carbon13)
  deallocate(veg_growth%reserve%carbon14)
  deallocate(veg_growth%reserve%nitrogen15)
  deallocate(veg_growth%fruit%carbon)
  deallocate(veg_growth%fruit%nitrogen)
  deallocate(veg_growth%fruit%phosphorus)
  deallocate(veg_growth%fruit%carbon13)
  deallocate(veg_growth%fruit%carbon14)
  deallocate(veg_growth%fruit%nitrogen15)
  deallocate(veg_growth%seed_bed%carbon)
  deallocate(veg_growth%seed_bed%nitrogen)
  deallocate(veg_growth%seed_bed%phosphorus)
  deallocate(veg_growth%seed_bed%carbon13)
  deallocate(veg_growth%seed_bed%carbon14)
  deallocate(veg_growth%seed_bed%nitrogen15)
#endif

#if defined(Dealloc_veg_litterfall)
  ! pass the values (back) to process memory
  veg__mem%litterfall%leaf%carbon%ptr             = veg_litterfall%leaf%carbon
  veg__mem%litterfall%leaf%nitrogen%ptr           = veg_litterfall%leaf%nitrogen
  veg__mem%litterfall%leaf%phosphorus%ptr         = veg_litterfall%leaf%phosphorus
  veg__mem%litterfall%leaf%carbon13%ptr           = veg_litterfall%leaf%carbon13
  veg__mem%litterfall%leaf%carbon14%ptr           = veg_litterfall%leaf%carbon14
  veg__mem%litterfall%leaf%nitrogen15%ptr         = veg_litterfall%leaf%nitrogen15
  veg__mem%litterfall%fine_root%carbon%ptr        = veg_litterfall%fine_root%carbon
  veg__mem%litterfall%fine_root%nitrogen%ptr      = veg_litterfall%fine_root%nitrogen
  veg__mem%litterfall%fine_root%phosphorus%ptr    = veg_litterfall%fine_root%phosphorus
  veg__mem%litterfall%fine_root%carbon13%ptr      = veg_litterfall%fine_root%carbon13
  veg__mem%litterfall%fine_root%carbon14%ptr      = veg_litterfall%fine_root%carbon14
  veg__mem%litterfall%fine_root%nitrogen15%ptr    = veg_litterfall%fine_root%nitrogen15
  veg__mem%litterfall%coarse_root%carbon%ptr      = veg_litterfall%coarse_root%carbon
  veg__mem%litterfall%coarse_root%nitrogen%ptr    = veg_litterfall%coarse_root%nitrogen
  veg__mem%litterfall%coarse_root%phosphorus%ptr  = veg_litterfall%coarse_root%phosphorus
  veg__mem%litterfall%coarse_root%carbon13%ptr    = veg_litterfall%coarse_root%carbon13
  veg__mem%litterfall%coarse_root%carbon14%ptr    = veg_litterfall%coarse_root%carbon14
  veg__mem%litterfall%coarse_root%nitrogen15%ptr  = veg_litterfall%coarse_root%nitrogen15
  veg__mem%litterfall%sap_wood%carbon%ptr         = veg_litterfall%sap_wood%carbon
  veg__mem%litterfall%sap_wood%nitrogen%ptr       = veg_litterfall%sap_wood%nitrogen
  veg__mem%litterfall%sap_wood%phosphorus%ptr     = veg_litterfall%sap_wood%phosphorus
  veg__mem%litterfall%sap_wood%carbon13%ptr       = veg_litterfall%sap_wood%carbon13
  veg__mem%litterfall%sap_wood%carbon14%ptr       = veg_litterfall%sap_wood%carbon14
  veg__mem%litterfall%sap_wood%nitrogen15%ptr     = veg_litterfall%sap_wood%nitrogen15
  veg__mem%litterfall%heart_wood%carbon%ptr       = veg_litterfall%heart_wood%carbon
  veg__mem%litterfall%heart_wood%nitrogen%ptr     = veg_litterfall%heart_wood%nitrogen
  veg__mem%litterfall%heart_wood%phosphorus%ptr   = veg_litterfall%heart_wood%phosphorus
  veg__mem%litterfall%heart_wood%carbon13%ptr     = veg_litterfall%heart_wood%carbon13
  veg__mem%litterfall%heart_wood%carbon14%ptr     = veg_litterfall%heart_wood%carbon14
  veg__mem%litterfall%heart_wood%nitrogen15%ptr   = veg_litterfall%heart_wood%nitrogen15
  veg__mem%litterfall%labile%carbon%ptr           = veg_litterfall%labile%carbon
  veg__mem%litterfall%labile%nitrogen%ptr         = veg_litterfall%labile%nitrogen
  veg__mem%litterfall%labile%phosphorus%ptr       = veg_litterfall%labile%phosphorus
  veg__mem%litterfall%labile%carbon13%ptr         = veg_litterfall%labile%carbon13
  veg__mem%litterfall%labile%carbon14%ptr         = veg_litterfall%labile%carbon14
  veg__mem%litterfall%labile%nitrogen15%ptr       = veg_litterfall%labile%nitrogen15
  veg__mem%litterfall%reserve%carbon%ptr          = veg_litterfall%reserve%carbon
  veg__mem%litterfall%reserve%nitrogen%ptr        = veg_litterfall%reserve%nitrogen
  veg__mem%litterfall%reserve%phosphorus%ptr      = veg_litterfall%reserve%phosphorus
  veg__mem%litterfall%reserve%carbon13%ptr        = veg_litterfall%reserve%carbon13
  veg__mem%litterfall%reserve%carbon14%ptr        = veg_litterfall%reserve%carbon14
  veg__mem%litterfall%reserve%nitrogen15%ptr      = veg_litterfall%reserve%nitrogen15
  veg__mem%litterfall%fruit%carbon%ptr            = veg_litterfall%fruit%carbon
  veg__mem%litterfall%fruit%nitrogen%ptr          = veg_litterfall%fruit%nitrogen
  veg__mem%litterfall%fruit%phosphorus%ptr        = veg_litterfall%fruit%phosphorus
  veg__mem%litterfall%fruit%carbon13%ptr          = veg_litterfall%fruit%carbon13
  veg__mem%litterfall%fruit%carbon14%ptr          = veg_litterfall%fruit%carbon14
  veg__mem%litterfall%fruit%nitrogen15%ptr        = veg_litterfall%fruit%nitrogen15
  veg__mem%litterfall%seed_bed%carbon%ptr         = veg_litterfall%seed_bed%carbon
  veg__mem%litterfall%seed_bed%nitrogen%ptr       = veg_litterfall%seed_bed%nitrogen
  veg__mem%litterfall%seed_bed%phosphorus%ptr     = veg_litterfall%seed_bed%phosphorus
  veg__mem%litterfall%seed_bed%carbon13%ptr       = veg_litterfall%seed_bed%carbon13
  veg__mem%litterfall%seed_bed%carbon14%ptr       = veg_litterfall%seed_bed%carbon14
  veg__mem%litterfall%seed_bed%nitrogen15%ptr     = veg_litterfall%seed_bed%nitrogen15
  ! deallocate variables in helper pool structures
  deallocate(veg_litterfall%leaf%carbon)
  deallocate(veg_litterfall%leaf%nitrogen)
  deallocate(veg_litterfall%leaf%phosphorus)
  deallocate(veg_litterfall%leaf%carbon13)
  deallocate(veg_litterfall%leaf%carbon14)
  deallocate(veg_litterfall%leaf%nitrogen15)
  deallocate(veg_litterfall%fine_root%carbon)
  deallocate(veg_litterfall%fine_root%nitrogen)
  deallocate(veg_litterfall%fine_root%phosphorus)
  deallocate(veg_litterfall%fine_root%carbon13)
  deallocate(veg_litterfall%fine_root%carbon14)
  deallocate(veg_litterfall%fine_root%nitrogen15)
  deallocate(veg_litterfall%coarse_root%carbon)
  deallocate(veg_litterfall%coarse_root%nitrogen)
  deallocate(veg_litterfall%coarse_root%phosphorus)
  deallocate(veg_litterfall%coarse_root%carbon13)
  deallocate(veg_litterfall%coarse_root%carbon14)
  deallocate(veg_litterfall%coarse_root%nitrogen15)
  deallocate(veg_litterfall%sap_wood%carbon)
  deallocate(veg_litterfall%sap_wood%nitrogen)
  deallocate(veg_litterfall%sap_wood%phosphorus)
  deallocate(veg_litterfall%sap_wood%carbon13)
  deallocate(veg_litterfall%sap_wood%carbon14)
  deallocate(veg_litterfall%sap_wood%nitrogen15)
  deallocate(veg_litterfall%heart_wood%carbon)
  deallocate(veg_litterfall%heart_wood%nitrogen)
  deallocate(veg_litterfall%heart_wood%phosphorus)
  deallocate(veg_litterfall%heart_wood%carbon13)
  deallocate(veg_litterfall%heart_wood%carbon14)
  deallocate(veg_litterfall%heart_wood%nitrogen15)
  deallocate(veg_litterfall%labile%carbon)
  deallocate(veg_litterfall%labile%nitrogen)
  deallocate(veg_litterfall%labile%phosphorus)
  deallocate(veg_litterfall%labile%carbon13)
  deallocate(veg_litterfall%labile%carbon14)
  deallocate(veg_litterfall%labile%nitrogen15)
  deallocate(veg_litterfall%reserve%carbon)
  deallocate(veg_litterfall%reserve%nitrogen)
  deallocate(veg_litterfall%reserve%phosphorus)
  deallocate(veg_litterfall%reserve%carbon13)
  deallocate(veg_litterfall%reserve%carbon14)
  deallocate(veg_litterfall%reserve%nitrogen15)
  deallocate(veg_litterfall%fruit%carbon)
  deallocate(veg_litterfall%fruit%nitrogen)
  deallocate(veg_litterfall%fruit%phosphorus)
  deallocate(veg_litterfall%fruit%carbon13)
  deallocate(veg_litterfall%fruit%carbon14)
  deallocate(veg_litterfall%fruit%nitrogen15)
  deallocate(veg_litterfall%seed_bed%carbon)
  deallocate(veg_litterfall%seed_bed%nitrogen)
  deallocate(veg_litterfall%seed_bed%phosphorus)
  deallocate(veg_litterfall%seed_bed%carbon13)
  deallocate(veg_litterfall%seed_bed%carbon14)
  deallocate(veg_litterfall%seed_bed%nitrogen15)
#endif

#if defined(Dealloc_veg_falloc)
  ! pass the values (back) to process memory
  veg__mem%falloc%leaf%carbon%ptr             = veg_falloc%leaf%carbon
  veg__mem%falloc%leaf%nitrogen%ptr           = veg_falloc%leaf%nitrogen
  veg__mem%falloc%leaf%phosphorus%ptr         = veg_falloc%leaf%phosphorus
  veg__mem%falloc%leaf%carbon13%ptr           = veg_falloc%leaf%carbon13
  veg__mem%falloc%leaf%carbon14%ptr           = veg_falloc%leaf%carbon14
  veg__mem%falloc%leaf%nitrogen15%ptr         = veg_falloc%leaf%nitrogen15
  veg__mem%falloc%fine_root%carbon%ptr        = veg_falloc%fine_root%carbon
  veg__mem%falloc%fine_root%nitrogen%ptr      = veg_falloc%fine_root%nitrogen
  veg__mem%falloc%fine_root%phosphorus%ptr    = veg_falloc%fine_root%phosphorus
  veg__mem%falloc%fine_root%carbon13%ptr      = veg_falloc%fine_root%carbon13
  veg__mem%falloc%fine_root%carbon14%ptr      = veg_falloc%fine_root%carbon14
  veg__mem%falloc%fine_root%nitrogen15%ptr    = veg_falloc%fine_root%nitrogen15
  veg__mem%falloc%coarse_root%carbon%ptr      = veg_falloc%coarse_root%carbon
  veg__mem%falloc%coarse_root%nitrogen%ptr    = veg_falloc%coarse_root%nitrogen
  veg__mem%falloc%coarse_root%phosphorus%ptr  = veg_falloc%coarse_root%phosphorus
  veg__mem%falloc%coarse_root%carbon13%ptr    = veg_falloc%coarse_root%carbon13
  veg__mem%falloc%coarse_root%carbon14%ptr    = veg_falloc%coarse_root%carbon14
  veg__mem%falloc%coarse_root%nitrogen15%ptr  = veg_falloc%coarse_root%nitrogen15
  veg__mem%falloc%sap_wood%carbon%ptr         = veg_falloc%sap_wood%carbon
  veg__mem%falloc%sap_wood%nitrogen%ptr       = veg_falloc%sap_wood%nitrogen
  veg__mem%falloc%sap_wood%phosphorus%ptr     = veg_falloc%sap_wood%phosphorus
  veg__mem%falloc%sap_wood%carbon13%ptr       = veg_falloc%sap_wood%carbon13
  veg__mem%falloc%sap_wood%carbon14%ptr       = veg_falloc%sap_wood%carbon14
  veg__mem%falloc%sap_wood%nitrogen15%ptr     = veg_falloc%sap_wood%nitrogen15
  veg__mem%falloc%heart_wood%carbon%ptr       = veg_falloc%heart_wood%carbon
  veg__mem%falloc%heart_wood%nitrogen%ptr     = veg_falloc%heart_wood%nitrogen
  veg__mem%falloc%heart_wood%phosphorus%ptr   = veg_falloc%heart_wood%phosphorus
  veg__mem%falloc%heart_wood%carbon13%ptr     = veg_falloc%heart_wood%carbon13
  veg__mem%falloc%heart_wood%carbon14%ptr     = veg_falloc%heart_wood%carbon14
  veg__mem%falloc%heart_wood%nitrogen15%ptr   = veg_falloc%heart_wood%nitrogen15
  veg__mem%falloc%labile%carbon%ptr           = veg_falloc%labile%carbon
  veg__mem%falloc%labile%nitrogen%ptr         = veg_falloc%labile%nitrogen
  veg__mem%falloc%labile%phosphorus%ptr       = veg_falloc%labile%phosphorus
  veg__mem%falloc%labile%carbon13%ptr         = veg_falloc%labile%carbon13
  veg__mem%falloc%labile%carbon14%ptr         = veg_falloc%labile%carbon14
  veg__mem%falloc%labile%nitrogen15%ptr       = veg_falloc%labile%nitrogen15
  veg__mem%falloc%reserve%carbon%ptr          = veg_falloc%reserve%carbon
  veg__mem%falloc%reserve%nitrogen%ptr        = veg_falloc%reserve%nitrogen
  veg__mem%falloc%reserve%phosphorus%ptr      = veg_falloc%reserve%phosphorus
  veg__mem%falloc%reserve%carbon13%ptr        = veg_falloc%reserve%carbon13
  veg__mem%falloc%reserve%carbon14%ptr        = veg_falloc%reserve%carbon14
  veg__mem%falloc%reserve%nitrogen15%ptr      = veg_falloc%reserve%nitrogen15
  veg__mem%falloc%fruit%carbon%ptr            = veg_falloc%fruit%carbon
  veg__mem%falloc%fruit%nitrogen%ptr          = veg_falloc%fruit%nitrogen
  veg__mem%falloc%fruit%phosphorus%ptr        = veg_falloc%fruit%phosphorus
  veg__mem%falloc%fruit%carbon13%ptr          = veg_falloc%fruit%carbon13
  veg__mem%falloc%fruit%carbon14%ptr          = veg_falloc%fruit%carbon14
  veg__mem%falloc%fruit%nitrogen15%ptr        = veg_falloc%fruit%nitrogen15
  veg__mem%falloc%seed_bed%carbon%ptr         = veg_falloc%seed_bed%carbon
  veg__mem%falloc%seed_bed%nitrogen%ptr       = veg_falloc%seed_bed%nitrogen
  veg__mem%falloc%seed_bed%phosphorus%ptr     = veg_falloc%seed_bed%phosphorus
  veg__mem%falloc%seed_bed%carbon13%ptr       = veg_falloc%seed_bed%carbon13
  veg__mem%falloc%seed_bed%carbon14%ptr       = veg_falloc%seed_bed%carbon14
  veg__mem%falloc%seed_bed%nitrogen15%ptr     = veg_falloc%seed_bed%nitrogen15
  ! deallocate variables in helper pool structures
  deallocate(veg_falloc%leaf%carbon)
  deallocate(veg_falloc%leaf%nitrogen)
  deallocate(veg_falloc%leaf%phosphorus)
  deallocate(veg_falloc%leaf%carbon13)
  deallocate(veg_falloc%leaf%carbon14)
  deallocate(veg_falloc%leaf%nitrogen15)
  deallocate(veg_falloc%fine_root%carbon)
  deallocate(veg_falloc%fine_root%nitrogen)
  deallocate(veg_falloc%fine_root%phosphorus)
  deallocate(veg_falloc%fine_root%carbon13)
  deallocate(veg_falloc%fine_root%carbon14)
  deallocate(veg_falloc%fine_root%nitrogen15)
  deallocate(veg_falloc%coarse_root%carbon)
  deallocate(veg_falloc%coarse_root%nitrogen)
  deallocate(veg_falloc%coarse_root%phosphorus)
  deallocate(veg_falloc%coarse_root%carbon13)
  deallocate(veg_falloc%coarse_root%carbon14)
  deallocate(veg_falloc%coarse_root%nitrogen15)
  deallocate(veg_falloc%sap_wood%carbon)
  deallocate(veg_falloc%sap_wood%nitrogen)
  deallocate(veg_falloc%sap_wood%phosphorus)
  deallocate(veg_falloc%sap_wood%carbon13)
  deallocate(veg_falloc%sap_wood%carbon14)
  deallocate(veg_falloc%sap_wood%nitrogen15)
  deallocate(veg_falloc%heart_wood%carbon)
  deallocate(veg_falloc%heart_wood%nitrogen)
  deallocate(veg_falloc%heart_wood%phosphorus)
  deallocate(veg_falloc%heart_wood%carbon13)
  deallocate(veg_falloc%heart_wood%carbon14)
  deallocate(veg_falloc%heart_wood%nitrogen15)
  deallocate(veg_falloc%labile%carbon)
  deallocate(veg_falloc%labile%nitrogen)
  deallocate(veg_falloc%labile%phosphorus)
  deallocate(veg_falloc%labile%carbon13)
  deallocate(veg_falloc%labile%carbon14)
  deallocate(veg_falloc%labile%nitrogen15)
  deallocate(veg_falloc%reserve%carbon)
  deallocate(veg_falloc%reserve%nitrogen)
  deallocate(veg_falloc%reserve%phosphorus)
  deallocate(veg_falloc%reserve%carbon13)
  deallocate(veg_falloc%reserve%carbon14)
  deallocate(veg_falloc%reserve%nitrogen15)
  deallocate(veg_falloc%fruit%carbon)
  deallocate(veg_falloc%fruit%nitrogen)
  deallocate(veg_falloc%fruit%phosphorus)
  deallocate(veg_falloc%fruit%carbon13)
  deallocate(veg_falloc%fruit%carbon14)
  deallocate(veg_falloc%fruit%nitrogen15)
  deallocate(veg_falloc%seed_bed%carbon)
  deallocate(veg_falloc%seed_bed%nitrogen)
  deallocate(veg_falloc%seed_bed%phosphorus)
  deallocate(veg_falloc%seed_bed%carbon13)
  deallocate(veg_falloc%seed_bed%carbon14)
  deallocate(veg_falloc%seed_bed%nitrogen15)
#endif

#if defined(Dealloc_veg_total_flux)
  ! pass the values (back) to process memory
  veg__mem%total_flux%leaf%carbon%ptr             = veg_total_flux%leaf%carbon
  veg__mem%total_flux%leaf%nitrogen%ptr           = veg_total_flux%leaf%nitrogen
  veg__mem%total_flux%leaf%phosphorus%ptr         = veg_total_flux%leaf%phosphorus
  veg__mem%total_flux%leaf%carbon13%ptr           = veg_total_flux%leaf%carbon13
  veg__mem%total_flux%leaf%carbon14%ptr           = veg_total_flux%leaf%carbon14
  veg__mem%total_flux%leaf%nitrogen15%ptr         = veg_total_flux%leaf%nitrogen15
  veg__mem%total_flux%fine_root%carbon%ptr        = veg_total_flux%fine_root%carbon
  veg__mem%total_flux%fine_root%nitrogen%ptr      = veg_total_flux%fine_root%nitrogen
  veg__mem%total_flux%fine_root%phosphorus%ptr    = veg_total_flux%fine_root%phosphorus
  veg__mem%total_flux%fine_root%carbon13%ptr      = veg_total_flux%fine_root%carbon13
  veg__mem%total_flux%fine_root%carbon14%ptr      = veg_total_flux%fine_root%carbon14
  veg__mem%total_flux%fine_root%nitrogen15%ptr    = veg_total_flux%fine_root%nitrogen15
  veg__mem%total_flux%coarse_root%carbon%ptr      = veg_total_flux%coarse_root%carbon
  veg__mem%total_flux%coarse_root%nitrogen%ptr    = veg_total_flux%coarse_root%nitrogen
  veg__mem%total_flux%coarse_root%phosphorus%ptr  = veg_total_flux%coarse_root%phosphorus
  veg__mem%total_flux%coarse_root%carbon13%ptr    = veg_total_flux%coarse_root%carbon13
  veg__mem%total_flux%coarse_root%carbon14%ptr    = veg_total_flux%coarse_root%carbon14
  veg__mem%total_flux%coarse_root%nitrogen15%ptr  = veg_total_flux%coarse_root%nitrogen15
  veg__mem%total_flux%sap_wood%carbon%ptr         = veg_total_flux%sap_wood%carbon
  veg__mem%total_flux%sap_wood%nitrogen%ptr       = veg_total_flux%sap_wood%nitrogen
  veg__mem%total_flux%sap_wood%phosphorus%ptr     = veg_total_flux%sap_wood%phosphorus
  veg__mem%total_flux%sap_wood%carbon13%ptr       = veg_total_flux%sap_wood%carbon13
  veg__mem%total_flux%sap_wood%carbon14%ptr       = veg_total_flux%sap_wood%carbon14
  veg__mem%total_flux%sap_wood%nitrogen15%ptr     = veg_total_flux%sap_wood%nitrogen15
  veg__mem%total_flux%heart_wood%carbon%ptr       = veg_total_flux%heart_wood%carbon
  veg__mem%total_flux%heart_wood%nitrogen%ptr     = veg_total_flux%heart_wood%nitrogen
  veg__mem%total_flux%heart_wood%phosphorus%ptr   = veg_total_flux%heart_wood%phosphorus
  veg__mem%total_flux%heart_wood%carbon13%ptr     = veg_total_flux%heart_wood%carbon13
  veg__mem%total_flux%heart_wood%carbon14%ptr     = veg_total_flux%heart_wood%carbon14
  veg__mem%total_flux%heart_wood%nitrogen15%ptr   = veg_total_flux%heart_wood%nitrogen15
  veg__mem%total_flux%labile%carbon%ptr           = veg_total_flux%labile%carbon
  veg__mem%total_flux%labile%nitrogen%ptr         = veg_total_flux%labile%nitrogen
  veg__mem%total_flux%labile%phosphorus%ptr       = veg_total_flux%labile%phosphorus
  veg__mem%total_flux%labile%carbon13%ptr         = veg_total_flux%labile%carbon13
  veg__mem%total_flux%labile%carbon14%ptr         = veg_total_flux%labile%carbon14
  veg__mem%total_flux%labile%nitrogen15%ptr       = veg_total_flux%labile%nitrogen15
  veg__mem%total_flux%reserve%carbon%ptr          = veg_total_flux%reserve%carbon
  veg__mem%total_flux%reserve%nitrogen%ptr        = veg_total_flux%reserve%nitrogen
  veg__mem%total_flux%reserve%phosphorus%ptr      = veg_total_flux%reserve%phosphorus
  veg__mem%total_flux%reserve%carbon13%ptr        = veg_total_flux%reserve%carbon13
  veg__mem%total_flux%reserve%carbon14%ptr        = veg_total_flux%reserve%carbon14
  veg__mem%total_flux%reserve%nitrogen15%ptr      = veg_total_flux%reserve%nitrogen15
  veg__mem%total_flux%fruit%carbon%ptr            = veg_total_flux%fruit%carbon
  veg__mem%total_flux%fruit%nitrogen%ptr          = veg_total_flux%fruit%nitrogen
  veg__mem%total_flux%fruit%phosphorus%ptr        = veg_total_flux%fruit%phosphorus
  veg__mem%total_flux%fruit%carbon13%ptr          = veg_total_flux%fruit%carbon13
  veg__mem%total_flux%fruit%carbon14%ptr          = veg_total_flux%fruit%carbon14
  veg__mem%total_flux%fruit%nitrogen15%ptr        = veg_total_flux%fruit%nitrogen15
  veg__mem%total_flux%seed_bed%carbon%ptr         = veg_total_flux%seed_bed%carbon
  veg__mem%total_flux%seed_bed%nitrogen%ptr       = veg_total_flux%seed_bed%nitrogen
  veg__mem%total_flux%seed_bed%phosphorus%ptr     = veg_total_flux%seed_bed%phosphorus
  veg__mem%total_flux%seed_bed%carbon13%ptr       = veg_total_flux%seed_bed%carbon13
  veg__mem%total_flux%seed_bed%carbon14%ptr       = veg_total_flux%seed_bed%carbon14
  veg__mem%total_flux%seed_bed%nitrogen15%ptr     = veg_total_flux%seed_bed%nitrogen15
  ! deallocate variables in helper pool structures
  deallocate(veg_total_flux%leaf%carbon)
  deallocate(veg_total_flux%leaf%nitrogen)
  deallocate(veg_total_flux%leaf%phosphorus)
  deallocate(veg_total_flux%leaf%carbon13)
  deallocate(veg_total_flux%leaf%carbon14)
  deallocate(veg_total_flux%leaf%nitrogen15)
  deallocate(veg_total_flux%fine_root%carbon)
  deallocate(veg_total_flux%fine_root%nitrogen)
  deallocate(veg_total_flux%fine_root%phosphorus)
  deallocate(veg_total_flux%fine_root%carbon13)
  deallocate(veg_total_flux%fine_root%carbon14)
  deallocate(veg_total_flux%fine_root%nitrogen15)
  deallocate(veg_total_flux%coarse_root%carbon)
  deallocate(veg_total_flux%coarse_root%nitrogen)
  deallocate(veg_total_flux%coarse_root%phosphorus)
  deallocate(veg_total_flux%coarse_root%carbon13)
  deallocate(veg_total_flux%coarse_root%carbon14)
  deallocate(veg_total_flux%coarse_root%nitrogen15)
  deallocate(veg_total_flux%sap_wood%carbon)
  deallocate(veg_total_flux%sap_wood%nitrogen)
  deallocate(veg_total_flux%sap_wood%phosphorus)
  deallocate(veg_total_flux%sap_wood%carbon13)
  deallocate(veg_total_flux%sap_wood%carbon14)
  deallocate(veg_total_flux%sap_wood%nitrogen15)
  deallocate(veg_total_flux%heart_wood%carbon)
  deallocate(veg_total_flux%heart_wood%nitrogen)
  deallocate(veg_total_flux%heart_wood%phosphorus)
  deallocate(veg_total_flux%heart_wood%carbon13)
  deallocate(veg_total_flux%heart_wood%carbon14)
  deallocate(veg_total_flux%heart_wood%nitrogen15)
  deallocate(veg_total_flux%labile%carbon)
  deallocate(veg_total_flux%labile%nitrogen)
  deallocate(veg_total_flux%labile%phosphorus)
  deallocate(veg_total_flux%labile%carbon13)
  deallocate(veg_total_flux%labile%carbon14)
  deallocate(veg_total_flux%labile%nitrogen15)
  deallocate(veg_total_flux%reserve%carbon)
  deallocate(veg_total_flux%reserve%nitrogen)
  deallocate(veg_total_flux%reserve%phosphorus)
  deallocate(veg_total_flux%reserve%carbon13)
  deallocate(veg_total_flux%reserve%carbon14)
  deallocate(veg_total_flux%reserve%nitrogen15)
  deallocate(veg_total_flux%fruit%carbon)
  deallocate(veg_total_flux%fruit%nitrogen)
  deallocate(veg_total_flux%fruit%phosphorus)
  deallocate(veg_total_flux%fruit%carbon13)
  deallocate(veg_total_flux%fruit%carbon14)
  deallocate(veg_total_flux%fruit%nitrogen15)
  deallocate(veg_total_flux%seed_bed%carbon)
  deallocate(veg_total_flux%seed_bed%nitrogen)
  deallocate(veg_total_flux%seed_bed%phosphorus)
  deallocate(veg_total_flux%seed_bed%carbon13)
  deallocate(veg_total_flux%seed_bed%carbon14)
  deallocate(veg_total_flux%seed_bed%nitrogen15)
#endif

#ifdef Dealloc_veg_exudation
  ! pass the values (back) to process memory
  veg__mem%exudation%carbon%ptr             = veg_exudation%carbon
  veg__mem%exudation%nitrogen%ptr           = veg_exudation%nitrogen
  veg__mem%exudation%phosphorus%ptr         = veg_exudation%phosphorus
  veg__mem%exudation%carbon13%ptr           = veg_exudation%carbon13
  veg__mem%exudation%carbon14%ptr           = veg_exudation%carbon14
  veg__mem%exudation%nitrogen15%ptr         = veg_exudation%nitrogen15
  ! deallocate variables in helper pool structures
  deallocate(veg_exudation%carbon)
  deallocate(veg_exudation%nitrogen)
  deallocate(veg_exudation%phosphorus)
  deallocate(veg_exudation%carbon13)
  deallocate(veg_exudation%carbon14)
  deallocate(veg_exudation%nitrogen15)
#endif

#ifdef Dealloc_veg_establishment
  ! pass the values (back) to process memory
  veg__mem%establishment%carbon%ptr             = veg_establishment%carbon
  veg__mem%establishment%nitrogen%ptr           = veg_establishment%nitrogen
  veg__mem%establishment%phosphorus%ptr         = veg_establishment%phosphorus
  veg__mem%establishment%carbon13%ptr           = veg_establishment%carbon13
  veg__mem%establishment%carbon14%ptr           = veg_establishment%carbon14
  veg__mem%establishment%nitrogen15%ptr         = veg_establishment%nitrogen15
  ! deallocate variables in helper pool structures
  deallocate(veg_establishment%carbon)
  deallocate(veg_establishment%nitrogen)
  deallocate(veg_establishment%phosphorus)
  deallocate(veg_establishment%carbon13)
  deallocate(veg_establishment%carbon14)
  deallocate(veg_establishment%nitrogen15)
#endif

#ifdef Dealloc_veg_reserve_use
  ! pass the values (back) to process memory
  veg__mem%reserve_use%carbon%ptr             = veg_reserve_use%carbon
  veg__mem%reserve_use%nitrogen%ptr           = veg_reserve_use%nitrogen
  veg__mem%reserve_use%phosphorus%ptr         = veg_reserve_use%phosphorus
  veg__mem%reserve_use%carbon13%ptr           = veg_reserve_use%carbon13
  veg__mem%reserve_use%carbon14%ptr           = veg_reserve_use%carbon14
  veg__mem%reserve_use%nitrogen15%ptr         = veg_reserve_use%nitrogen15
  ! deallocate variables in helper pool structures
  deallocate(veg_reserve_use%carbon)
  deallocate(veg_reserve_use%nitrogen)
  deallocate(veg_reserve_use%phosphorus)
  deallocate(veg_reserve_use%carbon13)
  deallocate(veg_reserve_use%carbon14)
  deallocate(veg_reserve_use%nitrogen15)
#endif

#ifdef Dealloc_veg_total_biomass
  ! pass the values (back) to process memory
  veg__mem%total_biomass%carbon%ptr             = veg_total_biomass%carbon
  veg__mem%total_biomass%nitrogen%ptr           = veg_total_biomass%nitrogen
  veg__mem%total_biomass%phosphorus%ptr         = veg_total_biomass%phosphorus
  veg__mem%total_biomass%carbon13%ptr           = veg_total_biomass%carbon13
  veg__mem%total_biomass%carbon14%ptr           = veg_total_biomass%carbon14
  veg__mem%total_biomass%nitrogen15%ptr         = veg_total_biomass%nitrogen15
  ! deallocate variables in helper pool structures
  deallocate(veg_total_biomass%carbon)
  deallocate(veg_total_biomass%nitrogen)
  deallocate(veg_total_biomass%phosphorus)
  deallocate(veg_total_biomass%carbon13)
  deallocate(veg_total_biomass%carbon14)
  deallocate(veg_total_biomass%nitrogen15)
#endif
