! ----------------------------------------------
! definition of structures for vegetation
! ----------------------------------------------

#ifdef Def_veg_pool
  TYPE(t_veg_pool_p)              :: veg_pool
#endif

#ifdef Def_veg_growth
  TYPE(t_veg_pool_p)              :: veg_growth
#endif

#ifdef Def_veg_litterfall
  TYPE(t_veg_pool_p)              :: veg_litterfall
#endif

#ifdef Def_veg_falloc
  TYPE(t_veg_pool_p)              :: veg_falloc
#endif

#ifdef Def_veg_total_flux
  TYPE(t_veg_pool_p)              :: veg_total_flux
#endif

#ifdef Def_veg_exudation
  TYPE(t_veg_element_p)           :: veg_exudation
#endif

#ifdef Def_veg_establishment
  TYPE(t_veg_element_p)           :: veg_establishment
#endif

#ifdef Def_veg_reserve_use
  TYPE(t_veg_element_p)           :: veg_reserve_use
#endif

#ifdef Def_veg_total_biomass
  TYPE(t_veg_element_p)           :: veg_total_biomass
#endif


! ----------------------------------------------
! set pointers to structures for vegetation
! ----------------------------------------------

#ifdef Get_veg_pool
  veg_pool%leaf%carbon              => veg__mem%pool%leaf%carbon%ptr
  veg_pool%leaf%nitrogen            => veg__mem%pool%leaf%nitrogen%ptr
  veg_pool%leaf%phosphorus          => veg__mem%pool%leaf%phosphorus%ptr
  veg_pool%leaf%carbon13            => veg__mem%pool%leaf%carbon13%ptr
  veg_pool%leaf%carbon14            => veg__mem%pool%leaf%carbon14%ptr
  veg_pool%leaf%nitrogen15          => veg__mem%pool%leaf%nitrogen15%ptr
  veg_pool%fine_root%carbon         => veg__mem%pool%fine_root%carbon%ptr
  veg_pool%fine_root%nitrogen       => veg__mem%pool%fine_root%nitrogen%ptr
  veg_pool%fine_root%phosphorus     => veg__mem%pool%fine_root%phosphorus%ptr
  veg_pool%fine_root%carbon13       => veg__mem%pool%fine_root%carbon13%ptr
  veg_pool%fine_root%carbon14       => veg__mem%pool%fine_root%carbon14%ptr
  veg_pool%fine_root%nitrogen15     => veg__mem%pool%fine_root%nitrogen15%ptr
  veg_pool%coarse_root%carbon       => veg__mem%pool%coarse_root%carbon%ptr
  veg_pool%coarse_root%nitrogen     => veg__mem%pool%coarse_root%nitrogen%ptr
  veg_pool%coarse_root%phosphorus   => veg__mem%pool%coarse_root%phosphorus%ptr
  veg_pool%coarse_root%carbon13     => veg__mem%pool%coarse_root%carbon13%ptr
  veg_pool%coarse_root%carbon14     => veg__mem%pool%coarse_root%carbon14%ptr
  veg_pool%coarse_root%nitrogen15   => veg__mem%pool%coarse_root%nitrogen15%ptr
  veg_pool%sap_wood%carbon          => veg__mem%pool%sap_wood%carbon%ptr
  veg_pool%sap_wood%nitrogen        => veg__mem%pool%sap_wood%nitrogen%ptr
  veg_pool%sap_wood%phosphorus      => veg__mem%pool%sap_wood%phosphorus%ptr
  veg_pool%sap_wood%carbon13        => veg__mem%pool%sap_wood%carbon13%ptr
  veg_pool%sap_wood%carbon14        => veg__mem%pool%sap_wood%carbon14%ptr
  veg_pool%sap_wood%nitrogen15      => veg__mem%pool%sap_wood%nitrogen15%ptr
  veg_pool%heart_wood%carbon        => veg__mem%pool%heart_wood%carbon%ptr
  veg_pool%heart_wood%nitrogen      => veg__mem%pool%heart_wood%nitrogen%ptr
  veg_pool%heart_wood%phosphorus    => veg__mem%pool%heart_wood%phosphorus%ptr
  veg_pool%heart_wood%carbon13      => veg__mem%pool%heart_wood%carbon13%ptr
  veg_pool%heart_wood%carbon14      => veg__mem%pool%heart_wood%carbon14%ptr
  veg_pool%heart_wood%nitrogen15    => veg__mem%pool%heart_wood%nitrogen15%ptr
  veg_pool%labile%carbon            => veg__mem%pool%labile%carbon%ptr
  veg_pool%labile%nitrogen          => veg__mem%pool%labile%nitrogen%ptr
  veg_pool%labile%phosphorus        => veg__mem%pool%labile%phosphorus%ptr
  veg_pool%labile%carbon13          => veg__mem%pool%labile%carbon13%ptr
  veg_pool%labile%carbon14          => veg__mem%pool%labile%carbon14%ptr
  veg_pool%labile%nitrogen15        => veg__mem%pool%labile%nitrogen15%ptr
  veg_pool%reserve%carbon           => veg__mem%pool%reserve%carbon%ptr
  veg_pool%reserve%nitrogen         => veg__mem%pool%reserve%nitrogen%ptr
  veg_pool%reserve%phosphorus       => veg__mem%pool%reserve%phosphorus%ptr
  veg_pool%reserve%carbon13         => veg__mem%pool%reserve%carbon13%ptr
  veg_pool%reserve%carbon14         => veg__mem%pool%reserve%carbon14%ptr
  veg_pool%reserve%nitrogen15       => veg__mem%pool%reserve%nitrogen15%ptr
  veg_pool%fruit%carbon             => veg__mem%pool%fruit%carbon%ptr
  veg_pool%fruit%nitrogen           => veg__mem%pool%fruit%nitrogen%ptr
  veg_pool%fruit%phosphorus         => veg__mem%pool%fruit%phosphorus%ptr
  veg_pool%fruit%carbon13           => veg__mem%pool%fruit%carbon13%ptr
  veg_pool%fruit%carbon14           => veg__mem%pool%fruit%carbon14%ptr
  veg_pool%fruit%nitrogen15         => veg__mem%pool%fruit%nitrogen15%ptr
  veg_pool%seed_bed%carbon          => veg__mem%pool%seed_bed%carbon%ptr
  veg_pool%seed_bed%nitrogen        => veg__mem%pool%seed_bed%nitrogen%ptr
  veg_pool%seed_bed%phosphorus      => veg__mem%pool%seed_bed%phosphorus%ptr
  veg_pool%seed_bed%carbon13        => veg__mem%pool%seed_bed%carbon13%ptr
  veg_pool%seed_bed%carbon14        => veg__mem%pool%seed_bed%carbon14%ptr
  veg_pool%seed_bed%nitrogen15      => veg__mem%pool%seed_bed%nitrogen15%ptr
#endif

#ifdef Get_veg_growth
  veg_growth%leaf%carbon              => veg__mem%growth%leaf%carbon%ptr
  veg_growth%leaf%nitrogen            => veg__mem%growth%leaf%nitrogen%ptr
  veg_growth%leaf%phosphorus          => veg__mem%growth%leaf%phosphorus%ptr
  veg_growth%leaf%carbon13            => veg__mem%growth%leaf%carbon13%ptr
  veg_growth%leaf%carbon14            => veg__mem%growth%leaf%carbon14%ptr
  veg_growth%leaf%nitrogen15          => veg__mem%growth%leaf%nitrogen15%ptr
  veg_growth%fine_root%carbon         => veg__mem%growth%fine_root%carbon%ptr
  veg_growth%fine_root%nitrogen       => veg__mem%growth%fine_root%nitrogen%ptr
  veg_growth%fine_root%phosphorus     => veg__mem%growth%fine_root%phosphorus%ptr
  veg_growth%fine_root%carbon13       => veg__mem%growth%fine_root%carbon13%ptr
  veg_growth%fine_root%carbon14       => veg__mem%growth%fine_root%carbon14%ptr
  veg_growth%fine_root%nitrogen15     => veg__mem%growth%fine_root%nitrogen15%ptr
  veg_growth%coarse_root%carbon       => veg__mem%growth%coarse_root%carbon%ptr
  veg_growth%coarse_root%nitrogen     => veg__mem%growth%coarse_root%nitrogen%ptr
  veg_growth%coarse_root%phosphorus   => veg__mem%growth%coarse_root%phosphorus%ptr
  veg_growth%coarse_root%carbon13     => veg__mem%growth%coarse_root%carbon13%ptr
  veg_growth%coarse_root%carbon14     => veg__mem%growth%coarse_root%carbon14%ptr
  veg_growth%coarse_root%nitrogen15   => veg__mem%growth%coarse_root%nitrogen15%ptr
  veg_growth%sap_wood%carbon          => veg__mem%growth%sap_wood%carbon%ptr
  veg_growth%sap_wood%nitrogen        => veg__mem%growth%sap_wood%nitrogen%ptr
  veg_growth%sap_wood%phosphorus      => veg__mem%growth%sap_wood%phosphorus%ptr
  veg_growth%sap_wood%carbon13        => veg__mem%growth%sap_wood%carbon13%ptr
  veg_growth%sap_wood%carbon14        => veg__mem%growth%sap_wood%carbon14%ptr
  veg_growth%sap_wood%nitrogen15      => veg__mem%growth%sap_wood%nitrogen15%ptr
  veg_growth%heart_wood%carbon        => veg__mem%growth%heart_wood%carbon%ptr
  veg_growth%heart_wood%nitrogen      => veg__mem%growth%heart_wood%nitrogen%ptr
  veg_growth%heart_wood%phosphorus    => veg__mem%growth%heart_wood%phosphorus%ptr
  veg_growth%heart_wood%carbon13      => veg__mem%growth%heart_wood%carbon13%ptr
  veg_growth%heart_wood%carbon14      => veg__mem%growth%heart_wood%carbon14%ptr
  veg_growth%heart_wood%nitrogen15    => veg__mem%growth%heart_wood%nitrogen15%ptr
  veg_growth%labile%carbon            => veg__mem%growth%labile%carbon%ptr
  veg_growth%labile%nitrogen          => veg__mem%growth%labile%nitrogen%ptr
  veg_growth%labile%phosphorus        => veg__mem%growth%labile%phosphorus%ptr
  veg_growth%labile%carbon13          => veg__mem%growth%labile%carbon13%ptr
  veg_growth%labile%carbon14          => veg__mem%growth%labile%carbon14%ptr
  veg_growth%labile%nitrogen15        => veg__mem%growth%labile%nitrogen15%ptr
  veg_growth%reserve%carbon           => veg__mem%growth%reserve%carbon%ptr
  veg_growth%reserve%nitrogen         => veg__mem%growth%reserve%nitrogen%ptr
  veg_growth%reserve%phosphorus       => veg__mem%growth%reserve%phosphorus%ptr
  veg_growth%reserve%carbon13         => veg__mem%growth%reserve%carbon13%ptr
  veg_growth%reserve%carbon14         => veg__mem%growth%reserve%carbon14%ptr
  veg_growth%reserve%nitrogen15       => veg__mem%growth%reserve%nitrogen15%ptr
  veg_growth%fruit%carbon             => veg__mem%growth%fruit%carbon%ptr
  veg_growth%fruit%nitrogen           => veg__mem%growth%fruit%nitrogen%ptr
  veg_growth%fruit%phosphorus         => veg__mem%growth%fruit%phosphorus%ptr
  veg_growth%fruit%carbon13           => veg__mem%growth%fruit%carbon13%ptr
  veg_growth%fruit%carbon14           => veg__mem%growth%fruit%carbon14%ptr
  veg_growth%fruit%nitrogen15         => veg__mem%growth%fruit%nitrogen15%ptr
  veg_growth%seed_bed%carbon          => veg__mem%growth%seed_bed%carbon%ptr
  veg_growth%seed_bed%nitrogen        => veg__mem%growth%seed_bed%nitrogen%ptr
  veg_growth%seed_bed%phosphorus      => veg__mem%growth%seed_bed%phosphorus%ptr
  veg_growth%seed_bed%carbon13        => veg__mem%growth%seed_bed%carbon13%ptr
  veg_growth%seed_bed%carbon14        => veg__mem%growth%seed_bed%carbon14%ptr
  veg_growth%seed_bed%nitrogen15      => veg__mem%growth%seed_bed%nitrogen15%ptr
#endif

#ifdef Get_veg_litterfall
  veg_litterfall%leaf%carbon              => veg__mem%litterfall%leaf%carbon%ptr
  veg_litterfall%leaf%nitrogen            => veg__mem%litterfall%leaf%nitrogen%ptr
  veg_litterfall%leaf%phosphorus          => veg__mem%litterfall%leaf%phosphorus%ptr
  veg_litterfall%leaf%carbon13            => veg__mem%litterfall%leaf%carbon13%ptr
  veg_litterfall%leaf%carbon14            => veg__mem%litterfall%leaf%carbon14%ptr
  veg_litterfall%leaf%nitrogen15          => veg__mem%litterfall%leaf%nitrogen15%ptr
  veg_litterfall%fine_root%carbon         => veg__mem%litterfall%fine_root%carbon%ptr
  veg_litterfall%fine_root%nitrogen       => veg__mem%litterfall%fine_root%nitrogen%ptr
  veg_litterfall%fine_root%phosphorus     => veg__mem%litterfall%fine_root%phosphorus%ptr
  veg_litterfall%fine_root%carbon13       => veg__mem%litterfall%fine_root%carbon13%ptr
  veg_litterfall%fine_root%carbon14       => veg__mem%litterfall%fine_root%carbon14%ptr
  veg_litterfall%fine_root%nitrogen15     => veg__mem%litterfall%fine_root%nitrogen15%ptr
  veg_litterfall%coarse_root%carbon       => veg__mem%litterfall%coarse_root%carbon%ptr
  veg_litterfall%coarse_root%nitrogen     => veg__mem%litterfall%coarse_root%nitrogen%ptr
  veg_litterfall%coarse_root%phosphorus   => veg__mem%litterfall%coarse_root%phosphorus%ptr
  veg_litterfall%coarse_root%carbon13     => veg__mem%litterfall%coarse_root%carbon13%ptr
  veg_litterfall%coarse_root%carbon14     => veg__mem%litterfall%coarse_root%carbon14%ptr
  veg_litterfall%coarse_root%nitrogen15   => veg__mem%litterfall%coarse_root%nitrogen15%ptr
  veg_litterfall%sap_wood%carbon          => veg__mem%litterfall%sap_wood%carbon%ptr
  veg_litterfall%sap_wood%nitrogen        => veg__mem%litterfall%sap_wood%nitrogen%ptr
  veg_litterfall%sap_wood%phosphorus      => veg__mem%litterfall%sap_wood%phosphorus%ptr
  veg_litterfall%sap_wood%carbon13        => veg__mem%litterfall%sap_wood%carbon13%ptr
  veg_litterfall%sap_wood%carbon14        => veg__mem%litterfall%sap_wood%carbon14%ptr
  veg_litterfall%sap_wood%nitrogen15      => veg__mem%litterfall%sap_wood%nitrogen15%ptr
  veg_litterfall%heart_wood%carbon        => veg__mem%litterfall%heart_wood%carbon%ptr
  veg_litterfall%heart_wood%nitrogen      => veg__mem%litterfall%heart_wood%nitrogen%ptr
  veg_litterfall%heart_wood%phosphorus    => veg__mem%litterfall%heart_wood%phosphorus%ptr
  veg_litterfall%heart_wood%carbon13      => veg__mem%litterfall%heart_wood%carbon13%ptr
  veg_litterfall%heart_wood%carbon14      => veg__mem%litterfall%heart_wood%carbon14%ptr
  veg_litterfall%heart_wood%nitrogen15    => veg__mem%litterfall%heart_wood%nitrogen15%ptr
  veg_litterfall%labile%carbon            => veg__mem%litterfall%labile%carbon%ptr
  veg_litterfall%labile%nitrogen          => veg__mem%litterfall%labile%nitrogen%ptr
  veg_litterfall%labile%phosphorus        => veg__mem%litterfall%labile%phosphorus%ptr
  veg_litterfall%labile%carbon13          => veg__mem%litterfall%labile%carbon13%ptr
  veg_litterfall%labile%carbon14          => veg__mem%litterfall%labile%carbon14%ptr
  veg_litterfall%labile%nitrogen15        => veg__mem%litterfall%labile%nitrogen15%ptr
  veg_litterfall%reserve%carbon           => veg__mem%litterfall%reserve%carbon%ptr
  veg_litterfall%reserve%nitrogen         => veg__mem%litterfall%reserve%nitrogen%ptr
  veg_litterfall%reserve%phosphorus       => veg__mem%litterfall%reserve%phosphorus%ptr
  veg_litterfall%reserve%carbon13         => veg__mem%litterfall%reserve%carbon13%ptr
  veg_litterfall%reserve%carbon14         => veg__mem%litterfall%reserve%carbon14%ptr
  veg_litterfall%reserve%nitrogen15       => veg__mem%litterfall%reserve%nitrogen15%ptr
  veg_litterfall%fruit%carbon             => veg__mem%litterfall%fruit%carbon%ptr
  veg_litterfall%fruit%nitrogen           => veg__mem%litterfall%fruit%nitrogen%ptr
  veg_litterfall%fruit%phosphorus         => veg__mem%litterfall%fruit%phosphorus%ptr
  veg_litterfall%fruit%carbon13           => veg__mem%litterfall%fruit%carbon13%ptr
  veg_litterfall%fruit%carbon14           => veg__mem%litterfall%fruit%carbon14%ptr
  veg_litterfall%fruit%nitrogen15         => veg__mem%litterfall%fruit%nitrogen15%ptr
  veg_litterfall%seed_bed%carbon          => veg__mem%litterfall%seed_bed%carbon%ptr
  veg_litterfall%seed_bed%nitrogen        => veg__mem%litterfall%seed_bed%nitrogen%ptr
  veg_litterfall%seed_bed%phosphorus      => veg__mem%litterfall%seed_bed%phosphorus%ptr
  veg_litterfall%seed_bed%carbon13        => veg__mem%litterfall%seed_bed%carbon13%ptr
  veg_litterfall%seed_bed%carbon14        => veg__mem%litterfall%seed_bed%carbon14%ptr
  veg_litterfall%seed_bed%nitrogen15      => veg__mem%litterfall%seed_bed%nitrogen15%ptr
#endif

#ifdef Get_veg_falloc
  veg_falloc%leaf%carbon              => veg__mem%falloc%leaf%carbon%ptr
  veg_falloc%leaf%nitrogen            => veg__mem%falloc%leaf%nitrogen%ptr
  veg_falloc%leaf%phosphorus          => veg__mem%falloc%leaf%phosphorus%ptr
  veg_falloc%leaf%carbon13            => veg__mem%falloc%leaf%carbon13%ptr
  veg_falloc%leaf%carbon14            => veg__mem%falloc%leaf%carbon14%ptr
  veg_falloc%leaf%nitrogen15          => veg__mem%falloc%leaf%nitrogen15%ptr
  veg_falloc%fine_root%carbon         => veg__mem%falloc%fine_root%carbon%ptr
  veg_falloc%fine_root%nitrogen       => veg__mem%falloc%fine_root%nitrogen%ptr
  veg_falloc%fine_root%phosphorus     => veg__mem%falloc%fine_root%phosphorus%ptr
  veg_falloc%fine_root%carbon13       => veg__mem%falloc%fine_root%carbon13%ptr
  veg_falloc%fine_root%carbon14       => veg__mem%falloc%fine_root%carbon14%ptr
  veg_falloc%fine_root%nitrogen15     => veg__mem%falloc%fine_root%nitrogen15%ptr
  veg_falloc%coarse_root%carbon       => veg__mem%falloc%coarse_root%carbon%ptr
  veg_falloc%coarse_root%nitrogen     => veg__mem%falloc%coarse_root%nitrogen%ptr
  veg_falloc%coarse_root%phosphorus   => veg__mem%falloc%coarse_root%phosphorus%ptr
  veg_falloc%coarse_root%carbon13     => veg__mem%falloc%coarse_root%carbon13%ptr
  veg_falloc%coarse_root%carbon14     => veg__mem%falloc%coarse_root%carbon14%ptr
  veg_falloc%coarse_root%nitrogen15   => veg__mem%falloc%coarse_root%nitrogen15%ptr
  veg_falloc%sap_wood%carbon          => veg__mem%falloc%sap_wood%carbon%ptr
  veg_falloc%sap_wood%nitrogen        => veg__mem%falloc%sap_wood%nitrogen%ptr
  veg_falloc%sap_wood%phosphorus      => veg__mem%falloc%sap_wood%phosphorus%ptr
  veg_falloc%sap_wood%carbon13        => veg__mem%falloc%sap_wood%carbon13%ptr
  veg_falloc%sap_wood%carbon14        => veg__mem%falloc%sap_wood%carbon14%ptr
  veg_falloc%sap_wood%nitrogen15      => veg__mem%falloc%sap_wood%nitrogen15%ptr
  veg_falloc%heart_wood%carbon        => veg__mem%falloc%heart_wood%carbon%ptr
  veg_falloc%heart_wood%nitrogen      => veg__mem%falloc%heart_wood%nitrogen%ptr
  veg_falloc%heart_wood%phosphorus    => veg__mem%falloc%heart_wood%phosphorus%ptr
  veg_falloc%heart_wood%carbon13      => veg__mem%falloc%heart_wood%carbon13%ptr
  veg_falloc%heart_wood%carbon14      => veg__mem%falloc%heart_wood%carbon14%ptr
  veg_falloc%heart_wood%nitrogen15    => veg__mem%falloc%heart_wood%nitrogen15%ptr
  veg_falloc%labile%carbon            => veg__mem%falloc%labile%carbon%ptr
  veg_falloc%labile%nitrogen          => veg__mem%falloc%labile%nitrogen%ptr
  veg_falloc%labile%phosphorus        => veg__mem%falloc%labile%phosphorus%ptr
  veg_falloc%labile%carbon13          => veg__mem%falloc%labile%carbon13%ptr
  veg_falloc%labile%carbon14          => veg__mem%falloc%labile%carbon14%ptr
  veg_falloc%labile%nitrogen15        => veg__mem%falloc%labile%nitrogen15%ptr
  veg_falloc%reserve%carbon           => veg__mem%falloc%reserve%carbon%ptr
  veg_falloc%reserve%nitrogen         => veg__mem%falloc%reserve%nitrogen%ptr
  veg_falloc%reserve%phosphorus       => veg__mem%falloc%reserve%phosphorus%ptr
  veg_falloc%reserve%carbon13         => veg__mem%falloc%reserve%carbon13%ptr
  veg_falloc%reserve%carbon14         => veg__mem%falloc%reserve%carbon14%ptr
  veg_falloc%reserve%nitrogen15       => veg__mem%falloc%reserve%nitrogen15%ptr
  veg_falloc%fruit%carbon             => veg__mem%falloc%fruit%carbon%ptr
  veg_falloc%fruit%nitrogen           => veg__mem%falloc%fruit%nitrogen%ptr
  veg_falloc%fruit%phosphorus         => veg__mem%falloc%fruit%phosphorus%ptr
  veg_falloc%fruit%carbon13           => veg__mem%falloc%fruit%carbon13%ptr
  veg_falloc%fruit%carbon14           => veg__mem%falloc%fruit%carbon14%ptr
  veg_falloc%fruit%nitrogen15         => veg__mem%falloc%fruit%nitrogen15%ptr
  veg_falloc%seed_bed%carbon          => veg__mem%falloc%seed_bed%carbon%ptr
  veg_falloc%seed_bed%nitrogen        => veg__mem%falloc%seed_bed%nitrogen%ptr
  veg_falloc%seed_bed%phosphorus      => veg__mem%falloc%seed_bed%phosphorus%ptr
  veg_falloc%seed_bed%carbon13        => veg__mem%falloc%seed_bed%carbon13%ptr
  veg_falloc%seed_bed%carbon14        => veg__mem%falloc%seed_bed%carbon14%ptr
  veg_falloc%seed_bed%nitrogen15      => veg__mem%falloc%seed_bed%nitrogen15%ptr
#endif

#ifdef Get_veg_total_flux
  veg_total_flux%leaf%carbon              => veg__mem%total_flux%leaf%carbon%ptr
  veg_total_flux%leaf%nitrogen            => veg__mem%total_flux%leaf%nitrogen%ptr
  veg_total_flux%leaf%phosphorus          => veg__mem%total_flux%leaf%phosphorus%ptr
  veg_total_flux%leaf%carbon13            => veg__mem%total_flux%leaf%carbon13%ptr
  veg_total_flux%leaf%carbon14            => veg__mem%total_flux%leaf%carbon14%ptr
  veg_total_flux%leaf%nitrogen15          => veg__mem%total_flux%leaf%nitrogen15%ptr
  veg_total_flux%fine_root%carbon         => veg__mem%total_flux%fine_root%carbon%ptr
  veg_total_flux%fine_root%nitrogen       => veg__mem%total_flux%fine_root%nitrogen%ptr
  veg_total_flux%fine_root%phosphorus     => veg__mem%total_flux%fine_root%phosphorus%ptr
  veg_total_flux%fine_root%carbon13       => veg__mem%total_flux%fine_root%carbon13%ptr
  veg_total_flux%fine_root%carbon14       => veg__mem%total_flux%fine_root%carbon14%ptr
  veg_total_flux%fine_root%nitrogen15     => veg__mem%total_flux%fine_root%nitrogen15%ptr
  veg_total_flux%coarse_root%carbon       => veg__mem%total_flux%coarse_root%carbon%ptr
  veg_total_flux%coarse_root%nitrogen     => veg__mem%total_flux%coarse_root%nitrogen%ptr
  veg_total_flux%coarse_root%phosphorus   => veg__mem%total_flux%coarse_root%phosphorus%ptr
  veg_total_flux%coarse_root%carbon13     => veg__mem%total_flux%coarse_root%carbon13%ptr
  veg_total_flux%coarse_root%carbon14     => veg__mem%total_flux%coarse_root%carbon14%ptr
  veg_total_flux%coarse_root%nitrogen15   => veg__mem%total_flux%coarse_root%nitrogen15%ptr
  veg_total_flux%sap_wood%carbon          => veg__mem%total_flux%sap_wood%carbon%ptr
  veg_total_flux%sap_wood%nitrogen        => veg__mem%total_flux%sap_wood%nitrogen%ptr
  veg_total_flux%sap_wood%phosphorus      => veg__mem%total_flux%sap_wood%phosphorus%ptr
  veg_total_flux%sap_wood%carbon13        => veg__mem%total_flux%sap_wood%carbon13%ptr
  veg_total_flux%sap_wood%carbon14        => veg__mem%total_flux%sap_wood%carbon14%ptr
  veg_total_flux%sap_wood%nitrogen15      => veg__mem%total_flux%sap_wood%nitrogen15%ptr
  veg_total_flux%heart_wood%carbon        => veg__mem%total_flux%heart_wood%carbon%ptr
  veg_total_flux%heart_wood%nitrogen      => veg__mem%total_flux%heart_wood%nitrogen%ptr
  veg_total_flux%heart_wood%phosphorus    => veg__mem%total_flux%heart_wood%phosphorus%ptr
  veg_total_flux%heart_wood%carbon13      => veg__mem%total_flux%heart_wood%carbon13%ptr
  veg_total_flux%heart_wood%carbon14      => veg__mem%total_flux%heart_wood%carbon14%ptr
  veg_total_flux%heart_wood%nitrogen15    => veg__mem%total_flux%heart_wood%nitrogen15%ptr
  veg_total_flux%labile%carbon            => veg__mem%total_flux%labile%carbon%ptr
  veg_total_flux%labile%nitrogen          => veg__mem%total_flux%labile%nitrogen%ptr
  veg_total_flux%labile%phosphorus        => veg__mem%total_flux%labile%phosphorus%ptr
  veg_total_flux%labile%carbon13          => veg__mem%total_flux%labile%carbon13%ptr
  veg_total_flux%labile%carbon14          => veg__mem%total_flux%labile%carbon14%ptr
  veg_total_flux%labile%nitrogen15        => veg__mem%total_flux%labile%nitrogen15%ptr
  veg_total_flux%reserve%carbon           => veg__mem%total_flux%reserve%carbon%ptr
  veg_total_flux%reserve%nitrogen         => veg__mem%total_flux%reserve%nitrogen%ptr
  veg_total_flux%reserve%phosphorus       => veg__mem%total_flux%reserve%phosphorus%ptr
  veg_total_flux%reserve%carbon13         => veg__mem%total_flux%reserve%carbon13%ptr
  veg_total_flux%reserve%carbon14         => veg__mem%total_flux%reserve%carbon14%ptr
  veg_total_flux%reserve%nitrogen15       => veg__mem%total_flux%reserve%nitrogen15%ptr
  veg_total_flux%fruit%carbon             => veg__mem%total_flux%fruit%carbon%ptr
  veg_total_flux%fruit%nitrogen           => veg__mem%total_flux%fruit%nitrogen%ptr
  veg_total_flux%fruit%phosphorus         => veg__mem%total_flux%fruit%phosphorus%ptr
  veg_total_flux%fruit%carbon13           => veg__mem%total_flux%fruit%carbon13%ptr
  veg_total_flux%fruit%carbon14           => veg__mem%total_flux%fruit%carbon14%ptr
  veg_total_flux%fruit%nitrogen15         => veg__mem%total_flux%fruit%nitrogen15%ptr
  veg_total_flux%seed_bed%carbon          => veg__mem%total_flux%seed_bed%carbon%ptr
  veg_total_flux%seed_bed%nitrogen        => veg__mem%total_flux%seed_bed%nitrogen%ptr
  veg_total_flux%seed_bed%phosphorus      => veg__mem%total_flux%seed_bed%phosphorus%ptr
  veg_total_flux%seed_bed%carbon13        => veg__mem%total_flux%seed_bed%carbon13%ptr
  veg_total_flux%seed_bed%carbon14        => veg__mem%total_flux%seed_bed%carbon14%ptr
  veg_total_flux%seed_bed%nitrogen15      => veg__mem%total_flux%seed_bed%nitrogen15%ptr
#endif

#ifdef Get_veg_exudation
  veg_exudation%carbon       => veg__mem%exudation%carbon%ptr
  veg_exudation%nitrogen     => veg__mem%exudation%nitrogen%ptr
  veg_exudation%phosphorus   => veg__mem%exudation%phosphorus%ptr
  veg_exudation%carbon13     => veg__mem%exudation%carbon13%ptr
  veg_exudation%carbon14     => veg__mem%exudation%carbon14%ptr
  veg_exudation%nitrogen15   => veg__mem%exudation%nitrogen15%ptr
  
  ! veg_exudation%carbon       => veg__mem%exudation%carbon%ptr(ics:ice,iblk)
  ! veg_exudation%nitrogen     => veg__mem%exudation%nitrogen%ptr(ics:ice,iblk)
  ! veg_exudation%phosphorus   => veg__mem%exudation%phosphorus%ptr(ics:ice,iblk)
  ! veg_exudation%carbon13     => veg__mem%exudation%carbon13%ptr(ics:ice,iblk)
  ! veg_exudation%carbon14     => veg__mem%exudation%carbon14%ptr(ics:ice,iblk)
  ! veg_exudation%nitrogen15   => veg__mem%exudation%nitrogen15%ptr(ics:ice,iblk)
#endif

#ifdef Get_veg_establishment
  veg_establishment%carbon       => veg__mem%establishment%carbon%ptr
  veg_establishment%nitrogen     => veg__mem%establishment%nitrogen%ptr
  veg_establishment%phosphorus   => veg__mem%establishment%phosphorus%ptr
  veg_establishment%carbon13     => veg__mem%establishment%carbon13%ptr
  veg_establishment%carbon14     => veg__mem%establishment%carbon14%ptr
  veg_establishment%nitrogen15   => veg__mem%establishment%nitrogen15%ptr
  
  ! veg_establishment%carbon       => veg__mem%establishment%carbon%ptr(ics:ice,iblk)
  ! veg_establishment%nitrogen     => veg__mem%establishment%nitrogen%ptr(ics:ice,iblk)
  ! veg_establishment%phosphorus   => veg__mem%establishment%phosphorus%ptr(ics:ice,iblk)
  ! veg_establishment%carbon13     => veg__mem%establishment%carbon13%ptr(ics:ice,iblk)
  ! veg_establishment%carbon14     => veg__mem%establishment%carbon14%ptr(ics:ice,iblk)
  ! veg_establishment%nitrogen15   => veg__mem%establishment%nitrogen15%ptr(ics:ice,iblk)
#endif

#ifdef Get_veg_reserve_use
  veg_reserve_use%carbon       => veg__mem%reserve_use%carbon%ptr
  veg_reserve_use%nitrogen     => veg__mem%reserve_use%nitrogen%ptr
  veg_reserve_use%phosphorus   => veg__mem%reserve_use%phosphorus%ptr
  veg_reserve_use%carbon13     => veg__mem%reserve_use%carbon13%ptr
  veg_reserve_use%carbon14     => veg__mem%reserve_use%carbon14%ptr
  veg_reserve_use%nitrogen15   => veg__mem%reserve_use%nitrogen15%ptr
  
  ! veg_reserve_use%carbon       => veg__mem%reserve_use%carbon%ptr(ics:ice,iblk)
  ! veg_reserve_use%nitrogen     => veg__mem%reserve_use%nitrogen%ptr(ics:ice,iblk)
  ! veg_reserve_use%phosphorus   => veg__mem%reserve_use%phosphorus%ptr(ics:ice,iblk)
  ! veg_reserve_use%carbon13     => veg__mem%reserve_use%carbon13%ptr(ics:ice,iblk)
  ! veg_reserve_use%carbon14     => veg__mem%reserve_use%carbon14%ptr(ics:ice,iblk)
  ! veg_reserve_use%nitrogen15   => veg__mem%reserve_use%nitrogen15%ptr(ics:ice,iblk)
#endif

#ifdef Get_veg_total_biomass
  veg_total_biomass%carbon       => veg__mem%total_biomass%carbon%ptr
  veg_total_biomass%nitrogen     => veg__mem%total_biomass%nitrogen%ptr
  veg_total_biomass%phosphorus   => veg__mem%total_biomass%phosphorus%ptr
  veg_total_biomass%carbon13     => veg__mem%total_biomass%carbon13%ptr
  veg_total_biomass%carbon14     => veg__mem%total_biomass%carbon14%ptr
  veg_total_biomass%nitrogen15   => veg__mem%total_biomass%nitrogen15%ptr
  
  ! veg_total_biomass%carbon       => veg__mem%total_biomass%carbon%ptr(ics:ice,iblk)
  ! veg_total_biomass%nitrogen     => veg__mem%total_biomass%nitrogen%ptr(ics:ice,iblk)
  ! veg_total_biomass%phosphorus   => veg__mem%total_biomass%phosphorus%ptr(ics:ice,iblk)
  ! veg_total_biomass%carbon13     => veg__mem%total_biomass%carbon13%ptr(ics:ice,iblk)
  ! veg_total_biomass%carbon14     => veg__mem%total_biomass%carbon14%ptr(ics:ice,iblk)
  ! veg_total_biomass%nitrogen15   => veg__mem%total_biomass%nitrogen15%ptr(ics:ice,iblk)
#endif
