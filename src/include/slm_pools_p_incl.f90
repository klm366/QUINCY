#ifdef Def_veg_pool
  TYPE(t_veg_element_p), TARGET   :: veg_pool_leaf    
  TYPE(t_veg_element_p), TARGET   :: veg_pool_fine_root
  TYPE(t_veg_element_p), TARGET   :: veg_pool_coarse_root
  TYPE(t_veg_element_p), TARGET   :: veg_pool_sap_wood
  TYPE(t_veg_element_p), TARGET   :: veg_pool_heart_wood
  TYPE(t_veg_element_p), TARGET   :: veg_pool_labile
  TYPE(t_veg_element_p), TARGET   :: veg_pool_reserve
  TYPE(t_veg_element_p), TARGET   :: veg_pool_fruit
  TYPE(t_veg_element_p), TARGET   :: veg_pool_seed_bed
  TYPE(t_veg_pool_p)              :: veg_pool
#endif

#ifdef Def_veg_growth
  TYPE(t_veg_element_p), TARGET   :: veg_growth_leaf   
  TYPE(t_veg_element_p), TARGET   :: veg_growth_fine_root
  TYPE(t_veg_element_p), TARGET   :: veg_growth_coarse_root
  TYPE(t_veg_element_p), TARGET   :: veg_growth_sap_wood
  TYPE(t_veg_element_p), TARGET   :: veg_growth_heart_wood
  TYPE(t_veg_element_p), TARGET   :: veg_growth_labile
  TYPE(t_veg_element_p), TARGET   :: veg_growth_reserve
  TYPE(t_veg_element_p), TARGET   :: veg_growth_fruit
  TYPE(t_veg_element_p), TARGET   :: veg_growth_seed_bed
  TYPE(t_veg_pool_p)              :: veg_growth
#endif

#ifdef Def_veg_litterfall
  TYPE(t_veg_element_p), TARGET   :: veg_litterfall_leaf   
  TYPE(t_veg_element_p), TARGET   :: veg_litterfall_fine_root
  TYPE(t_veg_element_p), TARGET   :: veg_litterfall_coarse_root
  TYPE(t_veg_element_p), TARGET   :: veg_litterfall_sap_wood
  TYPE(t_veg_element_p), TARGET   :: veg_litterfall_heart_wood
  TYPE(t_veg_element_p), TARGET   :: veg_litterfall_labile
  TYPE(t_veg_element_p), TARGET   :: veg_litterfall_reserve
  TYPE(t_veg_element_p), TARGET   :: veg_litterfall_fruit
  TYPE(t_veg_element_p), TARGET   :: veg_litterfall_seed_bed
  TYPE(t_veg_pool_p)              :: veg_litterfall
#endif

#ifdef Def_veg_falloc
  TYPE(t_veg_element_p), TARGET   :: veg_falloc_leaf   
  TYPE(t_veg_element_p), TARGET   :: veg_falloc_fine_root
  TYPE(t_veg_element_p), TARGET   :: veg_falloc_coarse_root
  TYPE(t_veg_element_p), TARGET   :: veg_falloc_sap_wood
  TYPE(t_veg_element_p), TARGET   :: veg_falloc_heart_wood
  TYPE(t_veg_element_p), TARGET   :: veg_falloc_labile
  TYPE(t_veg_element_p), TARGET   :: veg_falloc_reserve
  TYPE(t_veg_element_p), TARGET   :: veg_falloc_fruit
  TYPE(t_veg_element_p), TARGET   :: veg_falloc_seed_bed
  TYPE(t_veg_pool_p)              :: veg_falloc
#endif

#ifdef Def_veg_total_flux
  TYPE(t_veg_element_p), TARGET   :: veg_total_flux_leaf   
  TYPE(t_veg_element_p), TARGET   :: veg_total_flux_fine_root
  TYPE(t_veg_element_p), TARGET   :: veg_total_flux_coarse_root
  TYPE(t_veg_element_p), TARGET   :: veg_total_flux_sap_wood
  TYPE(t_veg_element_p), TARGET   :: veg_total_flux_heart_wood
  TYPE(t_veg_element_p), TARGET   :: veg_total_flux_labile
  TYPE(t_veg_element_p), TARGET   :: veg_total_flux_reserve
  TYPE(t_veg_element_p), TARGET   :: veg_total_flux_fruit
  TYPE(t_veg_element_p), TARGET   :: veg_total_flux_seed_bed
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


#ifdef Get_veg_pool
  nullify(veg_pool%leaf)
  nullify(veg_pool%fine_root)
  nullify(veg_pool%coarse_root)
  nullify(veg_pool%sap_wood)
  nullify(veg_pool%heart_wood)
  nullify(veg_pool%labile)
  nullify(veg_pool%reserve)
  nullify(veg_pool%fruit)
  nullify(veg_pool%seed_bed)

  nullify(veg_pool_leaf%carbon)
  nullify(veg_pool_leaf%nitrogen)
  nullify(veg_pool_leaf%phosphorus)
  nullify(veg_pool_leaf%carbon13)
  nullify(veg_pool_leaf%carbon14)
  nullify(veg_pool_leaf%nitrogen15)
  nullify(veg_pool_fine_root%carbon)
  nullify(veg_pool_fine_root%nitrogen)
  nullify(veg_pool_fine_root%phosphorus)
  nullify(veg_pool_fine_root%carbon13)
  nullify(veg_pool_fine_root%carbon14)
  nullify(veg_pool_fine_root%nitrogen15)
  nullify(veg_pool_coarse_root%carbon)
  nullify(veg_pool_coarse_root%nitrogen)
  nullify(veg_pool_coarse_root%phosphorus)
  nullify(veg_pool_coarse_root%carbon13)
  nullify(veg_pool_coarse_root%carbon14)
  nullify(veg_pool_coarse_root%nitrogen15)
  nullify(veg_pool_sap_wood%carbon)
  nullify(veg_pool_sap_wood%nitrogen)
  nullify(veg_pool_sap_wood%phosphorus)
  nullify(veg_pool_sap_wood%carbon13)
  nullify(veg_pool_sap_wood%carbon14)
  nullify(veg_pool_sap_wood%nitrogen15)
  nullify(veg_pool_heart_wood%carbon)
  nullify(veg_pool_heart_wood%nitrogen)
  nullify(veg_pool_heart_wood%phosphorus)
  nullify(veg_pool_heart_wood%carbon13)
  nullify(veg_pool_heart_wood%carbon14)
  nullify(veg_pool_heart_wood%nitrogen15)
  nullify(veg_pool_labile%carbon)
  nullify(veg_pool_labile%nitrogen)
  nullify(veg_pool_labile%phosphorus)
  nullify(veg_pool_labile%carbon13)
  nullify(veg_pool_labile%carbon14)
  nullify(veg_pool_labile%nitrogen15)
  nullify(veg_pool_reserve%carbon)
  nullify(veg_pool_reserve%nitrogen)
  nullify(veg_pool_reserve%phosphorus)
  nullify(veg_pool_reserve%carbon13)
  nullify(veg_pool_reserve%carbon14)
  nullify(veg_pool_reserve%nitrogen15)
  nullify(veg_pool_fruit%carbon)
  nullify(veg_pool_fruit%nitrogen)
  nullify(veg_pool_fruit%phosphorus)
  nullify(veg_pool_fruit%carbon13)
  nullify(veg_pool_fruit%carbon14)
  nullify(veg_pool_fruit%nitrogen15)
  nullify(veg_pool_seed_bed%carbon)
  nullify(veg_pool_seed_bed%nitrogen)
  nullify(veg_pool_seed_bed%phosphorus)
  nullify(veg_pool_seed_bed%carbon13)
  nullify(veg_pool_seed_bed%carbon14)
  nullify(veg_pool_seed_bed%nitrogen15)

  veg_pool%leaf         => veg_pool_leaf
  veg_pool%fine_root    => veg_pool_fine_root
  veg_pool%coarse_root  => veg_pool_coarse_root
  veg_pool%sap_wood     => veg_pool_sap_wood
  veg_pool%heart_wood   => veg_pool_heart_wood
  veg_pool%labile       => veg_pool_labile
  veg_pool%reserve      => veg_pool_reserve
  veg_pool%fruit        => veg_pool_fruit
  veg_pool%seed_bed     => veg_pool_seed_bed

  veg_pool_leaf%carbon              => veg__mem%pool%leaf%carbon%ptr
  veg_pool_leaf%nitrogen            => veg__mem%pool%leaf%nitrogen%ptr
  veg_pool_leaf%phosphorus          => veg__mem%pool%leaf%phosphorus%ptr
  veg_pool_leaf%carbon13            => veg__mem%pool%leaf%carbon13%ptr
  veg_pool_leaf%carbon14            => veg__mem%pool%leaf%carbon14%ptr
  veg_pool_leaf%nitrogen15          => veg__mem%pool%leaf%nitrogen15%ptr
  veg_pool_fine_root%carbon         => veg__mem%pool%fine_root%carbon%ptr
  veg_pool_fine_root%nitrogen       => veg__mem%pool%fine_root%nitrogen%ptr
  veg_pool_fine_root%phosphorus     => veg__mem%pool%fine_root%phosphorus%ptr
  veg_pool_fine_root%carbon13       => veg__mem%pool%fine_root%carbon13%ptr
  veg_pool_fine_root%carbon14       => veg__mem%pool%fine_root%carbon14%ptr
  veg_pool_fine_root%nitrogen15     => veg__mem%pool%fine_root%nitrogen15%ptr
  veg_pool_coarse_root%carbon       => veg__mem%pool%coarse_root%carbon%ptr
  veg_pool_coarse_root%nitrogen     => veg__mem%pool%coarse_root%nitrogen%ptr
  veg_pool_coarse_root%phosphorus   => veg__mem%pool%coarse_root%phosphorus%ptr
  veg_pool_coarse_root%carbon13     => veg__mem%pool%coarse_root%carbon13%ptr
  veg_pool_coarse_root%carbon14     => veg__mem%pool%coarse_root%carbon14%ptr
  veg_pool_coarse_root%nitrogen15   => veg__mem%pool%coarse_root%nitrogen15%ptr
  veg_pool_sap_wood%carbon          => veg__mem%pool%sap_wood%carbon%ptr
  veg_pool_sap_wood%nitrogen        => veg__mem%pool%sap_wood%nitrogen%ptr
  veg_pool_sap_wood%phosphorus      => veg__mem%pool%sap_wood%phosphorus%ptr
  veg_pool_sap_wood%carbon13        => veg__mem%pool%sap_wood%carbon13%ptr
  veg_pool_sap_wood%carbon14        => veg__mem%pool%sap_wood%carbon14%ptr
  veg_pool_sap_wood%nitrogen15      => veg__mem%pool%sap_wood%nitrogen15%ptr
  veg_pool_heart_wood%carbon        => veg__mem%pool%heart_wood%carbon%ptr
  veg_pool_heart_wood%nitrogen      => veg__mem%pool%heart_wood%nitrogen%ptr
  veg_pool_heart_wood%phosphorus    => veg__mem%pool%heart_wood%phosphorus%ptr
  veg_pool_heart_wood%carbon13      => veg__mem%pool%heart_wood%carbon13%ptr
  veg_pool_heart_wood%carbon14      => veg__mem%pool%heart_wood%carbon14%ptr
  veg_pool_heart_wood%nitrogen15    => veg__mem%pool%heart_wood%nitrogen15%ptr
  veg_pool_labile%carbon            => veg__mem%pool%labile%carbon%ptr
  veg_pool_labile%nitrogen          => veg__mem%pool%labile%nitrogen%ptr
  veg_pool_labile%phosphorus        => veg__mem%pool%labile%phosphorus%ptr
  veg_pool_labile%carbon13          => veg__mem%pool%labile%carbon13%ptr
  veg_pool_labile%carbon14          => veg__mem%pool%labile%carbon14%ptr
  veg_pool_labile%nitrogen15        => veg__mem%pool%labile%nitrogen15%ptr
  veg_pool_reserve%carbon           => veg__mem%pool%reserve%carbon%ptr
  veg_pool_reserve%nitrogen         => veg__mem%pool%reserve%nitrogen%ptr
  veg_pool_reserve%phosphorus       => veg__mem%pool%reserve%phosphorus%ptr
  veg_pool_reserve%carbon13         => veg__mem%pool%reserve%carbon13%ptr
  veg_pool_reserve%carbon14         => veg__mem%pool%reserve%carbon14%ptr
  veg_pool_reserve%nitrogen15       => veg__mem%pool%reserve%nitrogen15%ptr
  veg_pool_fruit%carbon             => veg__mem%pool%fruit%carbon%ptr
  veg_pool_fruit%nitrogen           => veg__mem%pool%fruit%nitrogen%ptr
  veg_pool_fruit%phosphorus         => veg__mem%pool%fruit%phosphorus%ptr
  veg_pool_fruit%carbon13           => veg__mem%pool%fruit%carbon13%ptr
  veg_pool_fruit%carbon14           => veg__mem%pool%fruit%carbon14%ptr
  veg_pool_fruit%nitrogen15         => veg__mem%pool%fruit%nitrogen15%ptr
  veg_pool_seed_bed%carbon          => veg__mem%pool%seed_bed%carbon%ptr
  veg_pool_seed_bed%nitrogen        => veg__mem%pool%seed_bed%nitrogen%ptr
  veg_pool_seed_bed%phosphorus      => veg__mem%pool%seed_bed%phosphorus%ptr
  veg_pool_seed_bed%carbon13        => veg__mem%pool%seed_bed%carbon13%ptr
  veg_pool_seed_bed%carbon14        => veg__mem%pool%seed_bed%carbon14%ptr
  veg_pool_seed_bed%nitrogen15      => veg__mem%pool%seed_bed%nitrogen15%ptr
#endif

#ifdef Get_veg_growth
  nullify(veg_growth%leaf)
  nullify(veg_growth%fine_root)
  nullify(veg_growth%coarse_root)
  nullify(veg_growth%sap_wood)
  nullify(veg_growth%heart_wood)
  nullify(veg_growth%labile)
  nullify(veg_growth%reserve)
  nullify(veg_growth%fruit)
  nullify(veg_growth%seed_bed)

  nullify(veg_growth_leaf%carbon)
  nullify(veg_growth_leaf%nitrogen)
  nullify(veg_growth_leaf%phosphorus)
  nullify(veg_growth_leaf%carbon13)
  nullify(veg_growth_leaf%carbon14)
  nullify(veg_growth_leaf%nitrogen15)
  nullify(veg_growth_fine_root%carbon)
  nullify(veg_growth_fine_root%nitrogen)
  nullify(veg_growth_fine_root%phosphorus)
  nullify(veg_growth_fine_root%carbon13)
  nullify(veg_growth_fine_root%carbon14)
  nullify(veg_growth_fine_root%nitrogen15)
  nullify(veg_growth_coarse_root%carbon)
  nullify(veg_growth_coarse_root%nitrogen)
  nullify(veg_growth_coarse_root%phosphorus)
  nullify(veg_growth_coarse_root%carbon13)
  nullify(veg_growth_coarse_root%carbon14)
  nullify(veg_growth_coarse_root%nitrogen15)
  nullify(veg_growth_sap_wood%carbon)
  nullify(veg_growth_sap_wood%nitrogen)
  nullify(veg_growth_sap_wood%phosphorus)
  nullify(veg_growth_sap_wood%carbon13)
  nullify(veg_growth_sap_wood%carbon14)
  nullify(veg_growth_sap_wood%nitrogen15)
  nullify(veg_growth_heart_wood%carbon)
  nullify(veg_growth_heart_wood%nitrogen)
  nullify(veg_growth_heart_wood%phosphorus)
  nullify(veg_growth_heart_wood%carbon13)
  nullify(veg_growth_heart_wood%carbon14)
  nullify(veg_growth_heart_wood%nitrogen15)
  nullify(veg_growth_labile%carbon)
  nullify(veg_growth_labile%nitrogen)
  nullify(veg_growth_labile%phosphorus)
  nullify(veg_growth_labile%carbon13)
  nullify(veg_growth_labile%carbon14)
  nullify(veg_growth_labile%nitrogen15)
  nullify(veg_growth_reserve%carbon)
  nullify(veg_growth_reserve%nitrogen)
  nullify(veg_growth_reserve%phosphorus)
  nullify(veg_growth_reserve%carbon13)
  nullify(veg_growth_reserve%carbon14)
  nullify(veg_growth_reserve%nitrogen15)
  nullify(veg_growth_fruit%carbon)
  nullify(veg_growth_fruit%nitrogen)
  nullify(veg_growth_fruit%phosphorus)
  nullify(veg_growth_fruit%carbon13)
  nullify(veg_growth_fruit%carbon14)
  nullify(veg_growth_fruit%nitrogen15)
  nullify(veg_growth_seed_bed%carbon)
  nullify(veg_growth_seed_bed%nitrogen)
  nullify(veg_growth_seed_bed%phosphorus)
  nullify(veg_growth_seed_bed%carbon13)
  nullify(veg_growth_seed_bed%carbon14)
  nullify(veg_growth_seed_bed%nitrogen15)

  veg_growth%leaf         => veg_growth_leaf
  veg_growth%fine_root    => veg_growth_fine_root
  veg_growth%coarse_root  => veg_growth_coarse_root
  veg_growth%sap_wood     => veg_growth_sap_wood
  veg_growth%heart_wood   => veg_growth_heart_wood
  veg_growth%labile       => veg_growth_labile
  veg_growth%reserve      => veg_growth_reserve
  veg_growth%fruit        => veg_growth_fruit
  veg_growth%seed_bed     => veg_growth_seed_bed

  veg_growth_leaf%carbon              => veg__mem%growth%leaf%carbon%ptr
  veg_growth_leaf%nitrogen            => veg__mem%growth%leaf%nitrogen%ptr
  veg_growth_leaf%phosphorus          => veg__mem%growth%leaf%phosphorus%ptr
  veg_growth_leaf%carbon13            => veg__mem%growth%leaf%carbon13%ptr
  veg_growth_leaf%carbon14            => veg__mem%growth%leaf%carbon14%ptr
  veg_growth_leaf%nitrogen15          => veg__mem%growth%leaf%nitrogen15%ptr
  veg_growth_fine_root%carbon         => veg__mem%growth%fine_root%carbon%ptr
  veg_growth_fine_root%nitrogen       => veg__mem%growth%fine_root%nitrogen%ptr
  veg_growth_fine_root%phosphorus     => veg__mem%growth%fine_root%phosphorus%ptr
  veg_growth_fine_root%carbon13       => veg__mem%growth%fine_root%carbon13%ptr
  veg_growth_fine_root%carbon14       => veg__mem%growth%fine_root%carbon14%ptr
  veg_growth_fine_root%nitrogen15     => veg__mem%growth%fine_root%nitrogen15%ptr
  veg_growth_coarse_root%carbon       => veg__mem%growth%coarse_root%carbon%ptr
  veg_growth_coarse_root%nitrogen     => veg__mem%growth%coarse_root%nitrogen%ptr
  veg_growth_coarse_root%phosphorus   => veg__mem%growth%coarse_root%phosphorus%ptr
  veg_growth_coarse_root%carbon13     => veg__mem%growth%coarse_root%carbon13%ptr
  veg_growth_coarse_root%carbon14     => veg__mem%growth%coarse_root%carbon14%ptr
  veg_growth_coarse_root%nitrogen15   => veg__mem%growth%coarse_root%nitrogen15%ptr
  veg_growth_sap_wood%carbon          => veg__mem%growth%sap_wood%carbon%ptr
  veg_growth_sap_wood%nitrogen        => veg__mem%growth%sap_wood%nitrogen%ptr
  veg_growth_sap_wood%phosphorus      => veg__mem%growth%sap_wood%phosphorus%ptr
  veg_growth_sap_wood%carbon13        => veg__mem%growth%sap_wood%carbon13%ptr
  veg_growth_sap_wood%carbon14        => veg__mem%growth%sap_wood%carbon14%ptr
  veg_growth_sap_wood%nitrogen15      => veg__mem%growth%sap_wood%nitrogen15%ptr
  veg_growth_heart_wood%carbon        => veg__mem%growth%heart_wood%carbon%ptr
  veg_growth_heart_wood%nitrogen      => veg__mem%growth%heart_wood%nitrogen%ptr
  veg_growth_heart_wood%phosphorus    => veg__mem%growth%heart_wood%phosphorus%ptr
  veg_growth_heart_wood%carbon13      => veg__mem%growth%heart_wood%carbon13%ptr
  veg_growth_heart_wood%carbon14      => veg__mem%growth%heart_wood%carbon14%ptr
  veg_growth_heart_wood%nitrogen15    => veg__mem%growth%heart_wood%nitrogen15%ptr
  veg_growth_labile%carbon            => veg__mem%growth%labile%carbon%ptr
  veg_growth_labile%nitrogen          => veg__mem%growth%labile%nitrogen%ptr
  veg_growth_labile%phosphorus        => veg__mem%growth%labile%phosphorus%ptr
  veg_growth_labile%carbon13          => veg__mem%growth%labile%carbon13%ptr
  veg_growth_labile%carbon14          => veg__mem%growth%labile%carbon14%ptr
  veg_growth_labile%nitrogen15        => veg__mem%growth%labile%nitrogen15%ptr
  veg_growth_reserve%carbon           => veg__mem%growth%reserve%carbon%ptr
  veg_growth_reserve%nitrogen         => veg__mem%growth%reserve%nitrogen%ptr
  veg_growth_reserve%phosphorus       => veg__mem%growth%reserve%phosphorus%ptr
  veg_growth_reserve%carbon13         => veg__mem%growth%reserve%carbon13%ptr
  veg_growth_reserve%carbon14         => veg__mem%growth%reserve%carbon14%ptr
  veg_growth_reserve%nitrogen15       => veg__mem%growth%reserve%nitrogen15%ptr
  veg_growth_fruit%carbon             => veg__mem%growth%fruit%carbon%ptr
  veg_growth_fruit%nitrogen           => veg__mem%growth%fruit%nitrogen%ptr
  veg_growth_fruit%phosphorus         => veg__mem%growth%fruit%phosphorus%ptr
  veg_growth_fruit%carbon13           => veg__mem%growth%fruit%carbon13%ptr
  veg_growth_fruit%carbon14           => veg__mem%growth%fruit%carbon14%ptr
  veg_growth_fruit%nitrogen15         => veg__mem%growth%fruit%nitrogen15%ptr
  veg_growth_seed_bed%carbon          => veg__mem%growth%seed_bed%carbon%ptr
  veg_growth_seed_bed%nitrogen        => veg__mem%growth%seed_bed%nitrogen%ptr
  veg_growth_seed_bed%phosphorus      => veg__mem%growth%seed_bed%phosphorus%ptr
  veg_growth_seed_bed%carbon13        => veg__mem%growth%seed_bed%carbon13%ptr
  veg_growth_seed_bed%carbon14        => veg__mem%growth%seed_bed%carbon14%ptr
  veg_growth_seed_bed%nitrogen15      => veg__mem%growth%seed_bed%nitrogen15%ptr
#endif

#ifdef Get_veg_litterfall
  nullify(veg_litterfall%leaf)
  nullify(veg_litterfall%fine_root)
  nullify(veg_litterfall%coarse_root)
  nullify(veg_litterfall%sap_wood)
  nullify(veg_litterfall%heart_wood)
  nullify(veg_litterfall%labile)
  nullify(veg_litterfall%reserve)
  nullify(veg_litterfall%fruit)
  nullify(veg_litterfall%seed_bed)

  nullify(veg_litterfall_leaf%carbon)
  nullify(veg_litterfall_leaf%nitrogen)
  nullify(veg_litterfall_leaf%phosphorus)
  nullify(veg_litterfall_leaf%carbon13)
  nullify(veg_litterfall_leaf%carbon14)
  nullify(veg_litterfall_leaf%nitrogen15)
  nullify(veg_litterfall_fine_root%carbon)
  nullify(veg_litterfall_fine_root%nitrogen)
  nullify(veg_litterfall_fine_root%phosphorus)
  nullify(veg_litterfall_fine_root%carbon13)
  nullify(veg_litterfall_fine_root%carbon14)
  nullify(veg_litterfall_fine_root%nitrogen15)
  nullify(veg_litterfall_coarse_root%carbon)
  nullify(veg_litterfall_coarse_root%nitrogen)
  nullify(veg_litterfall_coarse_root%phosphorus)
  nullify(veg_litterfall_coarse_root%carbon13)
  nullify(veg_litterfall_coarse_root%carbon14)
  nullify(veg_litterfall_coarse_root%nitrogen15)
  nullify(veg_litterfall_sap_wood%carbon)
  nullify(veg_litterfall_sap_wood%nitrogen)
  nullify(veg_litterfall_sap_wood%phosphorus)
  nullify(veg_litterfall_sap_wood%carbon13)
  nullify(veg_litterfall_sap_wood%carbon14)
  nullify(veg_litterfall_sap_wood%nitrogen15)
  nullify(veg_litterfall_heart_wood%carbon)
  nullify(veg_litterfall_heart_wood%nitrogen)
  nullify(veg_litterfall_heart_wood%phosphorus)
  nullify(veg_litterfall_heart_wood%carbon13)
  nullify(veg_litterfall_heart_wood%carbon14)
  nullify(veg_litterfall_heart_wood%nitrogen15)
  nullify(veg_litterfall_labile%carbon)
  nullify(veg_litterfall_labile%nitrogen)
  nullify(veg_litterfall_labile%phosphorus)
  nullify(veg_litterfall_labile%carbon13)
  nullify(veg_litterfall_labile%carbon14)
  nullify(veg_litterfall_labile%nitrogen15)
  nullify(veg_litterfall_reserve%carbon)
  nullify(veg_litterfall_reserve%nitrogen)
  nullify(veg_litterfall_reserve%phosphorus)
  nullify(veg_litterfall_reserve%carbon13)
  nullify(veg_litterfall_reserve%carbon14)
  nullify(veg_litterfall_reserve%nitrogen15)
  nullify(veg_litterfall_fruit%carbon)
  nullify(veg_litterfall_fruit%nitrogen)
  nullify(veg_litterfall_fruit%phosphorus)
  nullify(veg_litterfall_fruit%carbon13)
  nullify(veg_litterfall_fruit%carbon14)
  nullify(veg_litterfall_fruit%nitrogen15)
  nullify(veg_litterfall_seed_bed%carbon)
  nullify(veg_litterfall_seed_bed%nitrogen)
  nullify(veg_litterfall_seed_bed%phosphorus)
  nullify(veg_litterfall_seed_bed%carbon13)
  nullify(veg_litterfall_seed_bed%carbon14)
  nullify(veg_litterfall_seed_bed%nitrogen15)

  veg_litterfall%leaf         => veg_litterfall_leaf
  veg_litterfall%fine_root    => veg_litterfall_fine_root
  veg_litterfall%coarse_root  => veg_litterfall_coarse_root
  veg_litterfall%sap_wood     => veg_litterfall_sap_wood
  veg_litterfall%heart_wood   => veg_litterfall_heart_wood
  veg_litterfall%labile       => veg_litterfall_labile
  veg_litterfall%reserve      => veg_litterfall_reserve
  veg_litterfall%fruit        => veg_litterfall_fruit
  veg_litterfall%seed_bed     => veg_litterfall_seed_bed

  veg_litterfall_leaf%carbon              => veg__mem%litterfall%leaf%carbon%ptr
  veg_litterfall_leaf%nitrogen            => veg__mem%litterfall%leaf%nitrogen%ptr
  veg_litterfall_leaf%phosphorus          => veg__mem%litterfall%leaf%phosphorus%ptr
  veg_litterfall_leaf%carbon13            => veg__mem%litterfall%leaf%carbon13%ptr
  veg_litterfall_leaf%carbon14            => veg__mem%litterfall%leaf%carbon14%ptr
  veg_litterfall_leaf%nitrogen15          => veg__mem%litterfall%leaf%nitrogen15%ptr
  veg_litterfall_fine_root%carbon         => veg__mem%litterfall%fine_root%carbon%ptr
  veg_litterfall_fine_root%nitrogen       => veg__mem%litterfall%fine_root%nitrogen%ptr
  veg_litterfall_fine_root%phosphorus     => veg__mem%litterfall%fine_root%phosphorus%ptr
  veg_litterfall_fine_root%carbon13       => veg__mem%litterfall%fine_root%carbon13%ptr
  veg_litterfall_fine_root%carbon14       => veg__mem%litterfall%fine_root%carbon14%ptr
  veg_litterfall_fine_root%nitrogen15     => veg__mem%litterfall%fine_root%nitrogen15%ptr
  veg_litterfall_coarse_root%carbon       => veg__mem%litterfall%coarse_root%carbon%ptr
  veg_litterfall_coarse_root%nitrogen     => veg__mem%litterfall%coarse_root%nitrogen%ptr
  veg_litterfall_coarse_root%phosphorus   => veg__mem%litterfall%coarse_root%phosphorus%ptr
  veg_litterfall_coarse_root%carbon13     => veg__mem%litterfall%coarse_root%carbon13%ptr
  veg_litterfall_coarse_root%carbon14     => veg__mem%litterfall%coarse_root%carbon14%ptr
  veg_litterfall_coarse_root%nitrogen15   => veg__mem%litterfall%coarse_root%nitrogen15%ptr
  veg_litterfall_sap_wood%carbon          => veg__mem%litterfall%sap_wood%carbon%ptr
  veg_litterfall_sap_wood%nitrogen        => veg__mem%litterfall%sap_wood%nitrogen%ptr
  veg_litterfall_sap_wood%phosphorus      => veg__mem%litterfall%sap_wood%phosphorus%ptr
  veg_litterfall_sap_wood%carbon13        => veg__mem%litterfall%sap_wood%carbon13%ptr
  veg_litterfall_sap_wood%carbon14        => veg__mem%litterfall%sap_wood%carbon14%ptr
  veg_litterfall_sap_wood%nitrogen15      => veg__mem%litterfall%sap_wood%nitrogen15%ptr
  veg_litterfall_heart_wood%carbon        => veg__mem%litterfall%heart_wood%carbon%ptr
  veg_litterfall_heart_wood%nitrogen      => veg__mem%litterfall%heart_wood%nitrogen%ptr
  veg_litterfall_heart_wood%phosphorus    => veg__mem%litterfall%heart_wood%phosphorus%ptr
  veg_litterfall_heart_wood%carbon13      => veg__mem%litterfall%heart_wood%carbon13%ptr
  veg_litterfall_heart_wood%carbon14      => veg__mem%litterfall%heart_wood%carbon14%ptr
  veg_litterfall_heart_wood%nitrogen15    => veg__mem%litterfall%heart_wood%nitrogen15%ptr
  veg_litterfall_labile%carbon            => veg__mem%litterfall%labile%carbon%ptr
  veg_litterfall_labile%nitrogen          => veg__mem%litterfall%labile%nitrogen%ptr
  veg_litterfall_labile%phosphorus        => veg__mem%litterfall%labile%phosphorus%ptr
  veg_litterfall_labile%carbon13          => veg__mem%litterfall%labile%carbon13%ptr
  veg_litterfall_labile%carbon14          => veg__mem%litterfall%labile%carbon14%ptr
  veg_litterfall_labile%nitrogen15        => veg__mem%litterfall%labile%nitrogen15%ptr
  veg_litterfall_reserve%carbon           => veg__mem%litterfall%reserve%carbon%ptr
  veg_litterfall_reserve%nitrogen         => veg__mem%litterfall%reserve%nitrogen%ptr
  veg_litterfall_reserve%phosphorus       => veg__mem%litterfall%reserve%phosphorus%ptr
  veg_litterfall_reserve%carbon13         => veg__mem%litterfall%reserve%carbon13%ptr
  veg_litterfall_reserve%carbon14         => veg__mem%litterfall%reserve%carbon14%ptr
  veg_litterfall_reserve%nitrogen15       => veg__mem%litterfall%reserve%nitrogen15%ptr
  veg_litterfall_fruit%carbon             => veg__mem%litterfall%fruit%carbon%ptr
  veg_litterfall_fruit%nitrogen           => veg__mem%litterfall%fruit%nitrogen%ptr
  veg_litterfall_fruit%phosphorus         => veg__mem%litterfall%fruit%phosphorus%ptr
  veg_litterfall_fruit%carbon13           => veg__mem%litterfall%fruit%carbon13%ptr
  veg_litterfall_fruit%carbon14           => veg__mem%litterfall%fruit%carbon14%ptr
  veg_litterfall_fruit%nitrogen15         => veg__mem%litterfall%fruit%nitrogen15%ptr
  veg_litterfall_seed_bed%carbon          => veg__mem%litterfall%seed_bed%carbon%ptr
  veg_litterfall_seed_bed%nitrogen        => veg__mem%litterfall%seed_bed%nitrogen%ptr
  veg_litterfall_seed_bed%phosphorus      => veg__mem%litterfall%seed_bed%phosphorus%ptr
  veg_litterfall_seed_bed%carbon13        => veg__mem%litterfall%seed_bed%carbon13%ptr
  veg_litterfall_seed_bed%carbon14        => veg__mem%litterfall%seed_bed%carbon14%ptr
  veg_litterfall_seed_bed%nitrogen15      => veg__mem%litterfall%seed_bed%nitrogen15%ptr
#endif

#ifdef Get_veg_falloc
  nullify(veg_falloc%leaf)
  nullify(veg_falloc%fine_root)
  nullify(veg_falloc%coarse_root)
  nullify(veg_falloc%sap_wood)
  nullify(veg_falloc%heart_wood)
  nullify(veg_falloc%labile)
  nullify(veg_falloc%reserve)
  nullify(veg_falloc%fruit)
  nullify(veg_falloc%seed_bed)

  nullify(veg_falloc_leaf%carbon)
  nullify(veg_falloc_leaf%nitrogen)
  nullify(veg_falloc_leaf%phosphorus)
  nullify(veg_falloc_leaf%carbon13)
  nullify(veg_falloc_leaf%carbon14)
  nullify(veg_falloc_leaf%nitrogen15)
  nullify(veg_falloc_fine_root%carbon)
  nullify(veg_falloc_fine_root%nitrogen)
  nullify(veg_falloc_fine_root%phosphorus)
  nullify(veg_falloc_fine_root%carbon13)
  nullify(veg_falloc_fine_root%carbon14)
  nullify(veg_falloc_fine_root%nitrogen15)
  nullify(veg_falloc_coarse_root%carbon)
  nullify(veg_falloc_coarse_root%nitrogen)
  nullify(veg_falloc_coarse_root%phosphorus)
  nullify(veg_falloc_coarse_root%carbon13)
  nullify(veg_falloc_coarse_root%carbon14)
  nullify(veg_falloc_coarse_root%nitrogen15)
  nullify(veg_falloc_sap_wood%carbon)
  nullify(veg_falloc_sap_wood%nitrogen)
  nullify(veg_falloc_sap_wood%phosphorus)
  nullify(veg_falloc_sap_wood%carbon13)
  nullify(veg_falloc_sap_wood%carbon14)
  nullify(veg_falloc_sap_wood%nitrogen15)
  nullify(veg_falloc_heart_wood%carbon)
  nullify(veg_falloc_heart_wood%nitrogen)
  nullify(veg_falloc_heart_wood%phosphorus)
  nullify(veg_falloc_heart_wood%carbon13)
  nullify(veg_falloc_heart_wood%carbon14)
  nullify(veg_falloc_heart_wood%nitrogen15)
  nullify(veg_falloc_labile%carbon)
  nullify(veg_falloc_labile%nitrogen)
  nullify(veg_falloc_labile%phosphorus)
  nullify(veg_falloc_labile%carbon13)
  nullify(veg_falloc_labile%carbon14)
  nullify(veg_falloc_labile%nitrogen15)
  nullify(veg_falloc_reserve%carbon)
  nullify(veg_falloc_reserve%nitrogen)
  nullify(veg_falloc_reserve%phosphorus)
  nullify(veg_falloc_reserve%carbon13)
  nullify(veg_falloc_reserve%carbon14)
  nullify(veg_falloc_reserve%nitrogen15)
  nullify(veg_falloc_fruit%carbon)
  nullify(veg_falloc_fruit%nitrogen)
  nullify(veg_falloc_fruit%phosphorus)
  nullify(veg_falloc_fruit%carbon13)
  nullify(veg_falloc_fruit%carbon14)
  nullify(veg_falloc_fruit%nitrogen15)
  nullify(veg_falloc_seed_bed%carbon)
  nullify(veg_falloc_seed_bed%nitrogen)
  nullify(veg_falloc_seed_bed%phosphorus)
  nullify(veg_falloc_seed_bed%carbon13)
  nullify(veg_falloc_seed_bed%carbon14)
  nullify(veg_falloc_seed_bed%nitrogen15)

  veg_falloc%leaf         => veg_falloc_leaf
  veg_falloc%fine_root    => veg_falloc_fine_root
  veg_falloc%coarse_root  => veg_falloc_coarse_root
  veg_falloc%sap_wood     => veg_falloc_sap_wood
  veg_falloc%heart_wood   => veg_falloc_heart_wood
  veg_falloc%labile       => veg_falloc_labile
  veg_falloc%reserve      => veg_falloc_reserve
  veg_falloc%fruit        => veg_falloc_fruit
  veg_falloc%seed_bed     => veg_falloc_seed_bed

  veg_falloc_leaf%carbon              => veg__mem%falloc%leaf%carbon%ptr
  veg_falloc_leaf%nitrogen            => veg__mem%falloc%leaf%nitrogen%ptr
  veg_falloc_leaf%phosphorus          => veg__mem%falloc%leaf%phosphorus%ptr
  veg_falloc_leaf%carbon13            => veg__mem%falloc%leaf%carbon13%ptr
  veg_falloc_leaf%carbon14            => veg__mem%falloc%leaf%carbon14%ptr
  veg_falloc_leaf%nitrogen15          => veg__mem%falloc%leaf%nitrogen15%ptr
  veg_falloc_fine_root%carbon         => veg__mem%falloc%fine_root%carbon%ptr
  veg_falloc_fine_root%nitrogen       => veg__mem%falloc%fine_root%nitrogen%ptr
  veg_falloc_fine_root%phosphorus     => veg__mem%falloc%fine_root%phosphorus%ptr
  veg_falloc_fine_root%carbon13       => veg__mem%falloc%fine_root%carbon13%ptr
  veg_falloc_fine_root%carbon14       => veg__mem%falloc%fine_root%carbon14%ptr
  veg_falloc_fine_root%nitrogen15     => veg__mem%falloc%fine_root%nitrogen15%ptr
  veg_falloc_coarse_root%carbon       => veg__mem%falloc%coarse_root%carbon%ptr
  veg_falloc_coarse_root%nitrogen     => veg__mem%falloc%coarse_root%nitrogen%ptr
  veg_falloc_coarse_root%phosphorus   => veg__mem%falloc%coarse_root%phosphorus%ptr
  veg_falloc_coarse_root%carbon13     => veg__mem%falloc%coarse_root%carbon13%ptr
  veg_falloc_coarse_root%carbon14     => veg__mem%falloc%coarse_root%carbon14%ptr
  veg_falloc_coarse_root%nitrogen15   => veg__mem%falloc%coarse_root%nitrogen15%ptr
  veg_falloc_sap_wood%carbon          => veg__mem%falloc%sap_wood%carbon%ptr
  veg_falloc_sap_wood%nitrogen        => veg__mem%falloc%sap_wood%nitrogen%ptr
  veg_falloc_sap_wood%phosphorus      => veg__mem%falloc%sap_wood%phosphorus%ptr
  veg_falloc_sap_wood%carbon13        => veg__mem%falloc%sap_wood%carbon13%ptr
  veg_falloc_sap_wood%carbon14        => veg__mem%falloc%sap_wood%carbon14%ptr
  veg_falloc_sap_wood%nitrogen15      => veg__mem%falloc%sap_wood%nitrogen15%ptr
  veg_falloc_heart_wood%carbon        => veg__mem%falloc%heart_wood%carbon%ptr
  veg_falloc_heart_wood%nitrogen      => veg__mem%falloc%heart_wood%nitrogen%ptr
  veg_falloc_heart_wood%phosphorus    => veg__mem%falloc%heart_wood%phosphorus%ptr
  veg_falloc_heart_wood%carbon13      => veg__mem%falloc%heart_wood%carbon13%ptr
  veg_falloc_heart_wood%carbon14      => veg__mem%falloc%heart_wood%carbon14%ptr
  veg_falloc_heart_wood%nitrogen15    => veg__mem%falloc%heart_wood%nitrogen15%ptr
  veg_falloc_labile%carbon            => veg__mem%falloc%labile%carbon%ptr
  veg_falloc_labile%nitrogen          => veg__mem%falloc%labile%nitrogen%ptr
  veg_falloc_labile%phosphorus        => veg__mem%falloc%labile%phosphorus%ptr
  veg_falloc_labile%carbon13          => veg__mem%falloc%labile%carbon13%ptr
  veg_falloc_labile%carbon14          => veg__mem%falloc%labile%carbon14%ptr
  veg_falloc_labile%nitrogen15        => veg__mem%falloc%labile%nitrogen15%ptr
  veg_falloc_reserve%carbon           => veg__mem%falloc%reserve%carbon%ptr
  veg_falloc_reserve%nitrogen         => veg__mem%falloc%reserve%nitrogen%ptr
  veg_falloc_reserve%phosphorus       => veg__mem%falloc%reserve%phosphorus%ptr
  veg_falloc_reserve%carbon13         => veg__mem%falloc%reserve%carbon13%ptr
  veg_falloc_reserve%carbon14         => veg__mem%falloc%reserve%carbon14%ptr
  veg_falloc_reserve%nitrogen15       => veg__mem%falloc%reserve%nitrogen15%ptr
  veg_falloc_fruit%carbon             => veg__mem%falloc%fruit%carbon%ptr
  veg_falloc_fruit%nitrogen           => veg__mem%falloc%fruit%nitrogen%ptr
  veg_falloc_fruit%phosphorus         => veg__mem%falloc%fruit%phosphorus%ptr
  veg_falloc_fruit%carbon13           => veg__mem%falloc%fruit%carbon13%ptr
  veg_falloc_fruit%carbon14           => veg__mem%falloc%fruit%carbon14%ptr
  veg_falloc_fruit%nitrogen15         => veg__mem%falloc%fruit%nitrogen15%ptr
  veg_falloc_seed_bed%carbon          => veg__mem%falloc%seed_bed%carbon%ptr
  veg_falloc_seed_bed%nitrogen        => veg__mem%falloc%seed_bed%nitrogen%ptr
  veg_falloc_seed_bed%phosphorus      => veg__mem%falloc%seed_bed%phosphorus%ptr
  veg_falloc_seed_bed%carbon13        => veg__mem%falloc%seed_bed%carbon13%ptr
  veg_falloc_seed_bed%carbon14        => veg__mem%falloc%seed_bed%carbon14%ptr
  veg_falloc_seed_bed%nitrogen15      => veg__mem%falloc%seed_bed%nitrogen15%ptr
#endif

#ifdef Get_veg_total_flux
  nullify(veg_total_flux%leaf)
  nullify(veg_total_flux%fine_root)
  nullify(veg_total_flux%coarse_root)
  nullify(veg_total_flux%sap_wood)
  nullify(veg_total_flux%heart_wood)
  nullify(veg_total_flux%labile)
  nullify(veg_total_flux%reserve)
  nullify(veg_total_flux%fruit)
  nullify(veg_total_flux%seed_bed)

  nullify(veg_total_flux_leaf%carbon)
  nullify(veg_total_flux_leaf%nitrogen)
  nullify(veg_total_flux_leaf%phosphorus)
  nullify(veg_total_flux_leaf%carbon13)
  nullify(veg_total_flux_leaf%carbon14)
  nullify(veg_total_flux_leaf%nitrogen15)
  nullify(veg_total_flux_fine_root%carbon)
  nullify(veg_total_flux_fine_root%nitrogen)
  nullify(veg_total_flux_fine_root%phosphorus)
  nullify(veg_total_flux_fine_root%carbon13)
  nullify(veg_total_flux_fine_root%carbon14)
  nullify(veg_total_flux_fine_root%nitrogen15)
  nullify(veg_total_flux_coarse_root%carbon)
  nullify(veg_total_flux_coarse_root%nitrogen)
  nullify(veg_total_flux_coarse_root%phosphorus)
  nullify(veg_total_flux_coarse_root%carbon13)
  nullify(veg_total_flux_coarse_root%carbon14)
  nullify(veg_total_flux_coarse_root%nitrogen15)
  nullify(veg_total_flux_sap_wood%carbon)
  nullify(veg_total_flux_sap_wood%nitrogen)
  nullify(veg_total_flux_sap_wood%phosphorus)
  nullify(veg_total_flux_sap_wood%carbon13)
  nullify(veg_total_flux_sap_wood%carbon14)
  nullify(veg_total_flux_sap_wood%nitrogen15)
  nullify(veg_total_flux_heart_wood%carbon)
  nullify(veg_total_flux_heart_wood%nitrogen)
  nullify(veg_total_flux_heart_wood%phosphorus)
  nullify(veg_total_flux_heart_wood%carbon13)
  nullify(veg_total_flux_heart_wood%carbon14)
  nullify(veg_total_flux_heart_wood%nitrogen15)
  nullify(veg_total_flux_labile%carbon)
  nullify(veg_total_flux_labile%nitrogen)
  nullify(veg_total_flux_labile%phosphorus)
  nullify(veg_total_flux_labile%carbon13)
  nullify(veg_total_flux_labile%carbon14)
  nullify(veg_total_flux_labile%nitrogen15)
  nullify(veg_total_flux_reserve%carbon)
  nullify(veg_total_flux_reserve%nitrogen)
  nullify(veg_total_flux_reserve%phosphorus)
  nullify(veg_total_flux_reserve%carbon13)
  nullify(veg_total_flux_reserve%carbon14)
  nullify(veg_total_flux_reserve%nitrogen15)
  nullify(veg_total_flux_fruit%carbon)
  nullify(veg_total_flux_fruit%nitrogen)
  nullify(veg_total_flux_fruit%phosphorus)
  nullify(veg_total_flux_fruit%carbon13)
  nullify(veg_total_flux_fruit%carbon14)
  nullify(veg_total_flux_fruit%nitrogen15)
  nullify(veg_total_flux_seed_bed%carbon)
  nullify(veg_total_flux_seed_bed%nitrogen)
  nullify(veg_total_flux_seed_bed%phosphorus)
  nullify(veg_total_flux_seed_bed%carbon13)
  nullify(veg_total_flux_seed_bed%carbon14)
  nullify(veg_total_flux_seed_bed%nitrogen15)

  veg_total_flux%leaf         => veg_total_flux_leaf
  veg_total_flux%fine_root    => veg_total_flux_fine_root
  veg_total_flux%coarse_root  => veg_total_flux_coarse_root
  veg_total_flux%sap_wood     => veg_total_flux_sap_wood
  veg_total_flux%heart_wood   => veg_total_flux_heart_wood
  veg_total_flux%labile       => veg_total_flux_labile
  veg_total_flux%reserve      => veg_total_flux_reserve
  veg_total_flux%fruit        => veg_total_flux_fruit
  veg_total_flux%seed_bed     => veg_total_flux_seed_bed

  veg_total_flux_leaf%carbon              => veg__mem%total_flux%leaf%carbon%ptr
  veg_total_flux_leaf%nitrogen            => veg__mem%total_flux%leaf%nitrogen%ptr
  veg_total_flux_leaf%phosphorus          => veg__mem%total_flux%leaf%phosphorus%ptr
  veg_total_flux_leaf%carbon13            => veg__mem%total_flux%leaf%carbon13%ptr
  veg_total_flux_leaf%carbon14            => veg__mem%total_flux%leaf%carbon14%ptr
  veg_total_flux_leaf%nitrogen15          => veg__mem%total_flux%leaf%nitrogen15%ptr
  veg_total_flux_fine_root%carbon         => veg__mem%total_flux%fine_root%carbon%ptr
  veg_total_flux_fine_root%nitrogen       => veg__mem%total_flux%fine_root%nitrogen%ptr
  veg_total_flux_fine_root%phosphorus     => veg__mem%total_flux%fine_root%phosphorus%ptr
  veg_total_flux_fine_root%carbon13       => veg__mem%total_flux%fine_root%carbon13%ptr
  veg_total_flux_fine_root%carbon14       => veg__mem%total_flux%fine_root%carbon14%ptr
  veg_total_flux_fine_root%nitrogen15     => veg__mem%total_flux%fine_root%nitrogen15%ptr
  veg_total_flux_coarse_root%carbon       => veg__mem%total_flux%coarse_root%carbon%ptr
  veg_total_flux_coarse_root%nitrogen     => veg__mem%total_flux%coarse_root%nitrogen%ptr
  veg_total_flux_coarse_root%phosphorus   => veg__mem%total_flux%coarse_root%phosphorus%ptr
  veg_total_flux_coarse_root%carbon13     => veg__mem%total_flux%coarse_root%carbon13%ptr
  veg_total_flux_coarse_root%carbon14     => veg__mem%total_flux%coarse_root%carbon14%ptr
  veg_total_flux_coarse_root%nitrogen15   => veg__mem%total_flux%coarse_root%nitrogen15%ptr
  veg_total_flux_sap_wood%carbon          => veg__mem%total_flux%sap_wood%carbon%ptr
  veg_total_flux_sap_wood%nitrogen        => veg__mem%total_flux%sap_wood%nitrogen%ptr
  veg_total_flux_sap_wood%phosphorus      => veg__mem%total_flux%sap_wood%phosphorus%ptr
  veg_total_flux_sap_wood%carbon13        => veg__mem%total_flux%sap_wood%carbon13%ptr
  veg_total_flux_sap_wood%carbon14        => veg__mem%total_flux%sap_wood%carbon14%ptr
  veg_total_flux_sap_wood%nitrogen15      => veg__mem%total_flux%sap_wood%nitrogen15%ptr
  veg_total_flux_heart_wood%carbon        => veg__mem%total_flux%heart_wood%carbon%ptr
  veg_total_flux_heart_wood%nitrogen      => veg__mem%total_flux%heart_wood%nitrogen%ptr
  veg_total_flux_heart_wood%phosphorus    => veg__mem%total_flux%heart_wood%phosphorus%ptr
  veg_total_flux_heart_wood%carbon13      => veg__mem%total_flux%heart_wood%carbon13%ptr
  veg_total_flux_heart_wood%carbon14      => veg__mem%total_flux%heart_wood%carbon14%ptr
  veg_total_flux_heart_wood%nitrogen15    => veg__mem%total_flux%heart_wood%nitrogen15%ptr
  veg_total_flux_labile%carbon            => veg__mem%total_flux%labile%carbon%ptr
  veg_total_flux_labile%nitrogen          => veg__mem%total_flux%labile%nitrogen%ptr
  veg_total_flux_labile%phosphorus        => veg__mem%total_flux%labile%phosphorus%ptr
  veg_total_flux_labile%carbon13          => veg__mem%total_flux%labile%carbon13%ptr
  veg_total_flux_labile%carbon14          => veg__mem%total_flux%labile%carbon14%ptr
  veg_total_flux_labile%nitrogen15        => veg__mem%total_flux%labile%nitrogen15%ptr
  veg_total_flux_reserve%carbon           => veg__mem%total_flux%reserve%carbon%ptr
  veg_total_flux_reserve%nitrogen         => veg__mem%total_flux%reserve%nitrogen%ptr
  veg_total_flux_reserve%phosphorus       => veg__mem%total_flux%reserve%phosphorus%ptr
  veg_total_flux_reserve%carbon13         => veg__mem%total_flux%reserve%carbon13%ptr
  veg_total_flux_reserve%carbon14         => veg__mem%total_flux%reserve%carbon14%ptr
  veg_total_flux_reserve%nitrogen15       => veg__mem%total_flux%reserve%nitrogen15%ptr
  veg_total_flux_fruit%carbon             => veg__mem%total_flux%fruit%carbon%ptr
  veg_total_flux_fruit%nitrogen           => veg__mem%total_flux%fruit%nitrogen%ptr
  veg_total_flux_fruit%phosphorus         => veg__mem%total_flux%fruit%phosphorus%ptr
  veg_total_flux_fruit%carbon13           => veg__mem%total_flux%fruit%carbon13%ptr
  veg_total_flux_fruit%carbon14           => veg__mem%total_flux%fruit%carbon14%ptr
  veg_total_flux_fruit%nitrogen15         => veg__mem%total_flux%fruit%nitrogen15%ptr
  veg_total_flux_seed_bed%carbon          => veg__mem%total_flux%seed_bed%carbon%ptr
  veg_total_flux_seed_bed%nitrogen        => veg__mem%total_flux%seed_bed%nitrogen%ptr
  veg_total_flux_seed_bed%phosphorus      => veg__mem%total_flux%seed_bed%phosphorus%ptr
  veg_total_flux_seed_bed%carbon13        => veg__mem%total_flux%seed_bed%carbon13%ptr
  veg_total_flux_seed_bed%carbon14        => veg__mem%total_flux%seed_bed%carbon14%ptr
  veg_total_flux_seed_bed%nitrogen15      => veg__mem%total_flux%seed_bed%nitrogen15%ptr
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
