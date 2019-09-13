!>
!!  @par Copyright
!!    QUINCY terrestrial biosphere model with coupled carbon, nitrogen, and phosphorus cycles
!!    Copyright (C) 2019 Max Planck Institute for Biogeochemistry (MPI-BGC), Jena, Germany 
!!    Contributors: Caldararu S., Engel J., Kern M., Pallandt M., Schnur R., Thum T., Yu L., Zaehle S.
!!
!!    This file is part of the QUINCY model.
!!    The QUINCY model is free software: you can redistribute it and/or modify
!!    it under the terms of the GNU General Public License version 3,
!!    as published by the Free Software Foundation.
!!
!!    The QUINCY model is distributed in the hope that it will be useful,
!!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!!    GNU General Public License for more details.
!!
!!    You should have received a copy of the GNU General Public License
!!    along with the QUINCY model.  If not, see <https://www.gnu.org/licenses/>.
!! 

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
allocate(veg_pool%leaf%carbon            ,source=reshape(veg__mem%pool%leaf%carbon%ptr(ics:ice,iblk),[nc]))
allocate(veg_pool%leaf%nitrogen          ,source=reshape(veg__mem%pool%leaf%nitrogen%ptr(ics:ice,iblk),[nc]))
allocate(veg_pool%leaf%phosphorus        ,source=reshape(veg__mem%pool%leaf%phosphorus%ptr(ics:ice,iblk),[nc]))
allocate(veg_pool%leaf%carbon13          ,source=reshape(veg__mem%pool%leaf%carbon13%ptr(ics:ice,iblk),[nc]))
allocate(veg_pool%leaf%carbon14          ,source=reshape(veg__mem%pool%leaf%carbon14%ptr(ics:ice,iblk),[nc]))
allocate(veg_pool%leaf%nitrogen15        ,source=reshape(veg__mem%pool%leaf%nitrogen15%ptr(ics:ice,iblk),[nc]))
allocate(veg_pool%fine_root%carbon       ,source=reshape(veg__mem%pool%fine_root%carbon%ptr(ics:ice,iblk),[nc]))
allocate(veg_pool%fine_root%nitrogen     ,source=reshape(veg__mem%pool%fine_root%nitrogen%ptr(ics:ice,iblk),[nc]))
allocate(veg_pool%fine_root%phosphorus   ,source=reshape(veg__mem%pool%fine_root%phosphorus%ptr(ics:ice,iblk),[nc]))
allocate(veg_pool%fine_root%carbon13     ,source=reshape(veg__mem%pool%fine_root%carbon13%ptr(ics:ice,iblk),[nc]))
allocate(veg_pool%fine_root%carbon14     ,source=reshape(veg__mem%pool%fine_root%carbon14%ptr(ics:ice,iblk),[nc]))
allocate(veg_pool%fine_root%nitrogen15   ,source=reshape(veg__mem%pool%fine_root%nitrogen15%ptr(ics:ice,iblk),[nc]))
allocate(veg_pool%coarse_root%carbon     ,source=reshape(veg__mem%pool%coarse_root%carbon%ptr(ics:ice,iblk),[nc]))
allocate(veg_pool%coarse_root%nitrogen   ,source=reshape(veg__mem%pool%coarse_root%nitrogen%ptr(ics:ice,iblk),[nc]))
allocate(veg_pool%coarse_root%phosphorus ,source=reshape(veg__mem%pool%coarse_root%phosphorus%ptr(ics:ice,iblk),[nc]))
allocate(veg_pool%coarse_root%carbon13   ,source=reshape(veg__mem%pool%coarse_root%carbon13%ptr(ics:ice,iblk),[nc]))
allocate(veg_pool%coarse_root%carbon14   ,source=reshape(veg__mem%pool%coarse_root%carbon14%ptr(ics:ice,iblk),[nc]))
allocate(veg_pool%coarse_root%nitrogen15 ,source=reshape(veg__mem%pool%coarse_root%nitrogen15%ptr(ics:ice,iblk),[nc]))
allocate(veg_pool%sap_wood%carbon        ,source=reshape(veg__mem%pool%sap_wood%carbon%ptr(ics:ice,iblk),[nc]))
allocate(veg_pool%sap_wood%nitrogen      ,source=reshape(veg__mem%pool%sap_wood%nitrogen%ptr(ics:ice,iblk),[nc]))
allocate(veg_pool%sap_wood%phosphorus    ,source=reshape(veg__mem%pool%sap_wood%phosphorus%ptr(ics:ice,iblk),[nc]))
allocate(veg_pool%sap_wood%carbon13      ,source=reshape(veg__mem%pool%sap_wood%carbon13%ptr(ics:ice,iblk),[nc]))
allocate(veg_pool%sap_wood%carbon14      ,source=reshape(veg__mem%pool%sap_wood%carbon14%ptr(ics:ice,iblk),[nc]))
allocate(veg_pool%sap_wood%nitrogen15    ,source=reshape(veg__mem%pool%sap_wood%nitrogen15%ptr(ics:ice,iblk),[nc]))
allocate(veg_pool%heart_wood%carbon      ,source=reshape(veg__mem%pool%heart_wood%carbon%ptr(ics:ice,iblk),[nc]))
allocate(veg_pool%heart_wood%nitrogen    ,source=reshape(veg__mem%pool%heart_wood%nitrogen%ptr(ics:ice,iblk),[nc]))
allocate(veg_pool%heart_wood%phosphorus  ,source=reshape(veg__mem%pool%heart_wood%phosphorus%ptr(ics:ice,iblk),[nc]))
allocate(veg_pool%heart_wood%carbon13    ,source=reshape(veg__mem%pool%heart_wood%carbon13%ptr(ics:ice,iblk),[nc]))
allocate(veg_pool%heart_wood%carbon14    ,source=reshape(veg__mem%pool%heart_wood%carbon14%ptr(ics:ice,iblk),[nc]))
allocate(veg_pool%heart_wood%nitrogen15  ,source=reshape(veg__mem%pool%heart_wood%nitrogen15%ptr(ics:ice,iblk),[nc]))
allocate(veg_pool%labile%carbon          ,source=reshape(veg__mem%pool%labile%carbon%ptr(ics:ice,iblk),[nc]))
allocate(veg_pool%labile%nitrogen        ,source=reshape(veg__mem%pool%labile%nitrogen%ptr(ics:ice,iblk),[nc]))
allocate(veg_pool%labile%phosphorus      ,source=reshape(veg__mem%pool%labile%phosphorus%ptr(ics:ice,iblk),[nc]))
allocate(veg_pool%labile%carbon13        ,source=reshape(veg__mem%pool%labile%carbon13%ptr(ics:ice,iblk),[nc]))
allocate(veg_pool%labile%carbon14        ,source=reshape(veg__mem%pool%labile%carbon14%ptr(ics:ice,iblk),[nc]))
allocate(veg_pool%labile%nitrogen15      ,source=reshape(veg__mem%pool%labile%nitrogen15%ptr(ics:ice,iblk),[nc]))
allocate(veg_pool%reserve%carbon         ,source=reshape(veg__mem%pool%reserve%carbon%ptr(ics:ice,iblk),[nc]))
allocate(veg_pool%reserve%nitrogen       ,source=reshape(veg__mem%pool%reserve%nitrogen%ptr(ics:ice,iblk),[nc]))
allocate(veg_pool%reserve%phosphorus     ,source=reshape(veg__mem%pool%reserve%phosphorus%ptr(ics:ice,iblk),[nc]))
allocate(veg_pool%reserve%carbon13       ,source=reshape(veg__mem%pool%reserve%carbon13%ptr(ics:ice,iblk),[nc]))
allocate(veg_pool%reserve%carbon14       ,source=reshape(veg__mem%pool%reserve%carbon14%ptr(ics:ice,iblk),[nc]))
allocate(veg_pool%reserve%nitrogen15     ,source=reshape(veg__mem%pool%reserve%nitrogen15%ptr(ics:ice,iblk),[nc]))
allocate(veg_pool%fruit%carbon           ,source=reshape(veg__mem%pool%fruit%carbon%ptr(ics:ice,iblk),[nc]))
allocate(veg_pool%fruit%nitrogen         ,source=reshape(veg__mem%pool%fruit%nitrogen%ptr(ics:ice,iblk),[nc]))
allocate(veg_pool%fruit%phosphorus       ,source=reshape(veg__mem%pool%fruit%phosphorus%ptr(ics:ice,iblk),[nc]))
allocate(veg_pool%fruit%carbon13         ,source=reshape(veg__mem%pool%fruit%carbon13%ptr(ics:ice,iblk),[nc]))
allocate(veg_pool%fruit%carbon14         ,source=reshape(veg__mem%pool%fruit%carbon14%ptr(ics:ice,iblk),[nc]))
allocate(veg_pool%fruit%nitrogen15       ,source=reshape(veg__mem%pool%fruit%nitrogen15%ptr(ics:ice,iblk),[nc]))
allocate(veg_pool%seed_bed%carbon        ,source=reshape(veg__mem%pool%seed_bed%carbon%ptr(ics:ice,iblk),[nc]))
allocate(veg_pool%seed_bed%nitrogen      ,source=reshape(veg__mem%pool%seed_bed%nitrogen%ptr(ics:ice,iblk),[nc]))
allocate(veg_pool%seed_bed%phosphorus    ,source=reshape(veg__mem%pool%seed_bed%phosphorus%ptr(ics:ice,iblk),[nc]))
allocate(veg_pool%seed_bed%carbon13      ,source=reshape(veg__mem%pool%seed_bed%carbon13%ptr(ics:ice,iblk),[nc]))
allocate(veg_pool%seed_bed%carbon14      ,source=reshape(veg__mem%pool%seed_bed%carbon14%ptr(ics:ice,iblk),[nc]))
allocate(veg_pool%seed_bed%nitrogen15    ,source=reshape(veg__mem%pool%seed_bed%nitrogen15%ptr(ics:ice,iblk),[nc]))

#endif

#if defined(Get_veg_growth)
allocate(veg_growth%leaf%carbon            ,source=reshape(veg__mem%growth%leaf%carbon%ptr(ics:ice,iblk),[nc]))
allocate(veg_growth%leaf%nitrogen          ,source=reshape(veg__mem%growth%leaf%nitrogen%ptr(ics:ice,iblk),[nc]))
allocate(veg_growth%leaf%phosphorus        ,source=reshape(veg__mem%growth%leaf%phosphorus%ptr(ics:ice,iblk),[nc]))
allocate(veg_growth%leaf%carbon13          ,source=reshape(veg__mem%growth%leaf%carbon13%ptr(ics:ice,iblk),[nc]))
allocate(veg_growth%leaf%carbon14          ,source=reshape(veg__mem%growth%leaf%carbon14%ptr(ics:ice,iblk),[nc]))
allocate(veg_growth%leaf%nitrogen15        ,source=reshape(veg__mem%growth%leaf%nitrogen15%ptr(ics:ice,iblk),[nc]))
allocate(veg_growth%fine_root%carbon       ,source=reshape(veg__mem%growth%fine_root%carbon%ptr(ics:ice,iblk),[nc]))
allocate(veg_growth%fine_root%nitrogen     ,source=reshape(veg__mem%growth%fine_root%nitrogen%ptr(ics:ice,iblk),[nc]))
allocate(veg_growth%fine_root%phosphorus   ,source=reshape(veg__mem%growth%fine_root%phosphorus%ptr(ics:ice,iblk),[nc]))
allocate(veg_growth%fine_root%carbon13     ,source=reshape(veg__mem%growth%fine_root%carbon13%ptr(ics:ice,iblk),[nc]))
allocate(veg_growth%fine_root%carbon14     ,source=reshape(veg__mem%growth%fine_root%carbon14%ptr(ics:ice,iblk),[nc]))
allocate(veg_growth%fine_root%nitrogen15   ,source=reshape(veg__mem%growth%fine_root%nitrogen15%ptr(ics:ice,iblk),[nc]))
allocate(veg_growth%coarse_root%carbon     ,source=reshape(veg__mem%growth%coarse_root%carbon%ptr(ics:ice,iblk),[nc]))
allocate(veg_growth%coarse_root%nitrogen   ,source=reshape(veg__mem%growth%coarse_root%nitrogen%ptr(ics:ice,iblk),[nc]))
allocate(veg_growth%coarse_root%phosphorus ,source=reshape(veg__mem%growth%coarse_root%phosphorus%ptr(ics:ice,iblk),[nc]))
allocate(veg_growth%coarse_root%carbon13   ,source=reshape(veg__mem%growth%coarse_root%carbon13%ptr(ics:ice,iblk),[nc]))
allocate(veg_growth%coarse_root%carbon14   ,source=reshape(veg__mem%growth%coarse_root%carbon14%ptr(ics:ice,iblk),[nc]))
allocate(veg_growth%coarse_root%nitrogen15 ,source=reshape(veg__mem%growth%coarse_root%nitrogen15%ptr(ics:ice,iblk),[nc]))
allocate(veg_growth%sap_wood%carbon        ,source=reshape(veg__mem%growth%sap_wood%carbon%ptr(ics:ice,iblk),[nc]))
allocate(veg_growth%sap_wood%nitrogen      ,source=reshape(veg__mem%growth%sap_wood%nitrogen%ptr(ics:ice,iblk),[nc]))
allocate(veg_growth%sap_wood%phosphorus    ,source=reshape(veg__mem%growth%sap_wood%phosphorus%ptr(ics:ice,iblk),[nc]))
allocate(veg_growth%sap_wood%carbon13      ,source=reshape(veg__mem%growth%sap_wood%carbon13%ptr(ics:ice,iblk),[nc]))
allocate(veg_growth%sap_wood%carbon14      ,source=reshape(veg__mem%growth%sap_wood%carbon14%ptr(ics:ice,iblk),[nc]))
allocate(veg_growth%sap_wood%nitrogen15    ,source=reshape(veg__mem%growth%sap_wood%nitrogen15%ptr(ics:ice,iblk),[nc]))
allocate(veg_growth%heart_wood%carbon      ,source=reshape(veg__mem%growth%heart_wood%carbon%ptr(ics:ice,iblk),[nc]))
allocate(veg_growth%heart_wood%nitrogen    ,source=reshape(veg__mem%growth%heart_wood%nitrogen%ptr(ics:ice,iblk),[nc]))
allocate(veg_growth%heart_wood%phosphorus  ,source=reshape(veg__mem%growth%heart_wood%phosphorus%ptr(ics:ice,iblk),[nc]))
allocate(veg_growth%heart_wood%carbon13    ,source=reshape(veg__mem%growth%heart_wood%carbon13%ptr(ics:ice,iblk),[nc]))
allocate(veg_growth%heart_wood%carbon14    ,source=reshape(veg__mem%growth%heart_wood%carbon14%ptr(ics:ice,iblk),[nc]))
allocate(veg_growth%heart_wood%nitrogen15  ,source=reshape(veg__mem%growth%heart_wood%nitrogen15%ptr(ics:ice,iblk),[nc]))
allocate(veg_growth%labile%carbon          ,source=reshape(veg__mem%growth%labile%carbon%ptr(ics:ice,iblk),[nc]))
allocate(veg_growth%labile%nitrogen        ,source=reshape(veg__mem%growth%labile%nitrogen%ptr(ics:ice,iblk),[nc]))
allocate(veg_growth%labile%phosphorus      ,source=reshape(veg__mem%growth%labile%phosphorus%ptr(ics:ice,iblk),[nc]))
allocate(veg_growth%labile%carbon13        ,source=reshape(veg__mem%growth%labile%carbon13%ptr(ics:ice,iblk),[nc]))
allocate(veg_growth%labile%carbon14        ,source=reshape(veg__mem%growth%labile%carbon14%ptr(ics:ice,iblk),[nc]))
allocate(veg_growth%labile%nitrogen15      ,source=reshape(veg__mem%growth%labile%nitrogen15%ptr(ics:ice,iblk),[nc]))
allocate(veg_growth%reserve%carbon         ,source=reshape(veg__mem%growth%reserve%carbon%ptr(ics:ice,iblk),[nc]))
allocate(veg_growth%reserve%nitrogen       ,source=reshape(veg__mem%growth%reserve%nitrogen%ptr(ics:ice,iblk),[nc]))
allocate(veg_growth%reserve%phosphorus     ,source=reshape(veg__mem%growth%reserve%phosphorus%ptr(ics:ice,iblk),[nc]))
allocate(veg_growth%reserve%carbon13       ,source=reshape(veg__mem%growth%reserve%carbon13%ptr(ics:ice,iblk),[nc]))
allocate(veg_growth%reserve%carbon14       ,source=reshape(veg__mem%growth%reserve%carbon14%ptr(ics:ice,iblk),[nc]))
allocate(veg_growth%reserve%nitrogen15     ,source=reshape(veg__mem%growth%reserve%nitrogen15%ptr(ics:ice,iblk),[nc]))
allocate(veg_growth%fruit%carbon           ,source=reshape(veg__mem%growth%fruit%carbon%ptr(ics:ice,iblk),[nc]))
allocate(veg_growth%fruit%nitrogen         ,source=reshape(veg__mem%growth%fruit%nitrogen%ptr(ics:ice,iblk),[nc]))
allocate(veg_growth%fruit%phosphorus       ,source=reshape(veg__mem%growth%fruit%phosphorus%ptr(ics:ice,iblk),[nc]))
allocate(veg_growth%fruit%carbon13         ,source=reshape(veg__mem%growth%fruit%carbon13%ptr(ics:ice,iblk),[nc]))
allocate(veg_growth%fruit%carbon14         ,source=reshape(veg__mem%growth%fruit%carbon14%ptr(ics:ice,iblk),[nc]))
allocate(veg_growth%fruit%nitrogen15       ,source=reshape(veg__mem%growth%fruit%nitrogen15%ptr(ics:ice,iblk),[nc]))
allocate(veg_growth%seed_bed%carbon        ,source=reshape(veg__mem%growth%seed_bed%carbon%ptr(ics:ice,iblk),[nc]))
allocate(veg_growth%seed_bed%nitrogen      ,source=reshape(veg__mem%growth%seed_bed%nitrogen%ptr(ics:ice,iblk),[nc]))
allocate(veg_growth%seed_bed%phosphorus    ,source=reshape(veg__mem%growth%seed_bed%phosphorus%ptr(ics:ice,iblk),[nc]))
allocate(veg_growth%seed_bed%carbon13      ,source=reshape(veg__mem%growth%seed_bed%carbon13%ptr(ics:ice,iblk),[nc]))
allocate(veg_growth%seed_bed%carbon14      ,source=reshape(veg__mem%growth%seed_bed%carbon14%ptr(ics:ice,iblk),[nc]))
allocate(veg_growth%seed_bed%nitrogen15    ,source=reshape(veg__mem%growth%seed_bed%nitrogen15%ptr(ics:ice,iblk),[nc]))

#endif

#if defined(Get_veg_litterfall)
allocate(veg_litterfall%leaf%carbon            ,source=reshape(veg__mem%litterfall%leaf%carbon%ptr(ics:ice,iblk),[nc]))
allocate(veg_litterfall%leaf%nitrogen          ,source=reshape(veg__mem%litterfall%leaf%nitrogen%ptr(ics:ice,iblk),[nc]))
allocate(veg_litterfall%leaf%phosphorus        ,source=reshape(veg__mem%litterfall%leaf%phosphorus%ptr(ics:ice,iblk),[nc]))
allocate(veg_litterfall%leaf%carbon13          ,source=reshape(veg__mem%litterfall%leaf%carbon13%ptr(ics:ice,iblk),[nc]))
allocate(veg_litterfall%leaf%carbon14          ,source=reshape(veg__mem%litterfall%leaf%carbon14%ptr(ics:ice,iblk),[nc]))
allocate(veg_litterfall%leaf%nitrogen15        ,source=reshape(veg__mem%litterfall%leaf%nitrogen15%ptr(ics:ice,iblk),[nc]))
allocate(veg_litterfall%fine_root%carbon       ,source=reshape(veg__mem%litterfall%fine_root%carbon%ptr(ics:ice,iblk),[nc]))
allocate(veg_litterfall%fine_root%nitrogen     ,source=reshape(veg__mem%litterfall%fine_root%nitrogen%ptr(ics:ice,iblk),[nc]))
allocate(veg_litterfall%fine_root%phosphorus   ,source=reshape(veg__mem%litterfall%fine_root%phosphorus%ptr(ics:ice,iblk),[nc]))
allocate(veg_litterfall%fine_root%carbon13     ,source=reshape(veg__mem%litterfall%fine_root%carbon13%ptr(ics:ice,iblk),[nc]))
allocate(veg_litterfall%fine_root%carbon14     ,source=reshape(veg__mem%litterfall%fine_root%carbon14%ptr(ics:ice,iblk),[nc]))
allocate(veg_litterfall%fine_root%nitrogen15   ,source=reshape(veg__mem%litterfall%fine_root%nitrogen15%ptr(ics:ice,iblk),[nc]))
allocate(veg_litterfall%coarse_root%carbon     ,source=reshape(veg__mem%litterfall%coarse_root%carbon%ptr(ics:ice,iblk),[nc]))
allocate(veg_litterfall%coarse_root%nitrogen   ,source=reshape(veg__mem%litterfall%coarse_root%nitrogen%ptr(ics:ice,iblk),[nc]))
allocate(veg_litterfall%coarse_root%phosphorus ,source=reshape(veg__mem%litterfall%coarse_root%phosphorus%ptr(ics:ice,iblk),[nc]))
allocate(veg_litterfall%coarse_root%carbon13   ,source=reshape(veg__mem%litterfall%coarse_root%carbon13%ptr(ics:ice,iblk),[nc]))
allocate(veg_litterfall%coarse_root%carbon14   ,source=reshape(veg__mem%litterfall%coarse_root%carbon14%ptr(ics:ice,iblk),[nc]))
allocate(veg_litterfall%coarse_root%nitrogen15 ,source=reshape(veg__mem%litterfall%coarse_root%nitrogen15%ptr(ics:ice,iblk),[nc]))
allocate(veg_litterfall%sap_wood%carbon        ,source=reshape(veg__mem%litterfall%sap_wood%carbon%ptr(ics:ice,iblk),[nc]))
allocate(veg_litterfall%sap_wood%nitrogen      ,source=reshape(veg__mem%litterfall%sap_wood%nitrogen%ptr(ics:ice,iblk),[nc]))
allocate(veg_litterfall%sap_wood%phosphorus    ,source=reshape(veg__mem%litterfall%sap_wood%phosphorus%ptr(ics:ice,iblk),[nc]))
allocate(veg_litterfall%sap_wood%carbon13      ,source=reshape(veg__mem%litterfall%sap_wood%carbon13%ptr(ics:ice,iblk),[nc]))
allocate(veg_litterfall%sap_wood%carbon14      ,source=reshape(veg__mem%litterfall%sap_wood%carbon14%ptr(ics:ice,iblk),[nc]))
allocate(veg_litterfall%sap_wood%nitrogen15    ,source=reshape(veg__mem%litterfall%sap_wood%nitrogen15%ptr(ics:ice,iblk),[nc]))
allocate(veg_litterfall%heart_wood%carbon      ,source=reshape(veg__mem%litterfall%heart_wood%carbon%ptr(ics:ice,iblk),[nc]))
allocate(veg_litterfall%heart_wood%nitrogen    ,source=reshape(veg__mem%litterfall%heart_wood%nitrogen%ptr(ics:ice,iblk),[nc]))
allocate(veg_litterfall%heart_wood%phosphorus  ,source=reshape(veg__mem%litterfall%heart_wood%phosphorus%ptr(ics:ice,iblk),[nc]))
allocate(veg_litterfall%heart_wood%carbon13    ,source=reshape(veg__mem%litterfall%heart_wood%carbon13%ptr(ics:ice,iblk),[nc]))
allocate(veg_litterfall%heart_wood%carbon14    ,source=reshape(veg__mem%litterfall%heart_wood%carbon14%ptr(ics:ice,iblk),[nc]))
allocate(veg_litterfall%heart_wood%nitrogen15  ,source=reshape(veg__mem%litterfall%heart_wood%nitrogen15%ptr(ics:ice,iblk),[nc]))
allocate(veg_litterfall%labile%carbon          ,source=reshape(veg__mem%litterfall%labile%carbon%ptr(ics:ice,iblk),[nc]))
allocate(veg_litterfall%labile%nitrogen        ,source=reshape(veg__mem%litterfall%labile%nitrogen%ptr(ics:ice,iblk),[nc]))
allocate(veg_litterfall%labile%phosphorus      ,source=reshape(veg__mem%litterfall%labile%phosphorus%ptr(ics:ice,iblk),[nc]))
allocate(veg_litterfall%labile%carbon13        ,source=reshape(veg__mem%litterfall%labile%carbon13%ptr(ics:ice,iblk),[nc]))
allocate(veg_litterfall%labile%carbon14        ,source=reshape(veg__mem%litterfall%labile%carbon14%ptr(ics:ice,iblk),[nc]))
allocate(veg_litterfall%labile%nitrogen15      ,source=reshape(veg__mem%litterfall%labile%nitrogen15%ptr(ics:ice,iblk),[nc]))
allocate(veg_litterfall%reserve%carbon         ,source=reshape(veg__mem%litterfall%reserve%carbon%ptr(ics:ice,iblk),[nc]))
allocate(veg_litterfall%reserve%nitrogen       ,source=reshape(veg__mem%litterfall%reserve%nitrogen%ptr(ics:ice,iblk),[nc]))
allocate(veg_litterfall%reserve%phosphorus     ,source=reshape(veg__mem%litterfall%reserve%phosphorus%ptr(ics:ice,iblk),[nc]))
allocate(veg_litterfall%reserve%carbon13       ,source=reshape(veg__mem%litterfall%reserve%carbon13%ptr(ics:ice,iblk),[nc]))
allocate(veg_litterfall%reserve%carbon14       ,source=reshape(veg__mem%litterfall%reserve%carbon14%ptr(ics:ice,iblk),[nc]))
allocate(veg_litterfall%reserve%nitrogen15     ,source=reshape(veg__mem%litterfall%reserve%nitrogen15%ptr(ics:ice,iblk),[nc]))
allocate(veg_litterfall%fruit%carbon           ,source=reshape(veg__mem%litterfall%fruit%carbon%ptr(ics:ice,iblk),[nc]))
allocate(veg_litterfall%fruit%nitrogen         ,source=reshape(veg__mem%litterfall%fruit%nitrogen%ptr(ics:ice,iblk),[nc]))
allocate(veg_litterfall%fruit%phosphorus       ,source=reshape(veg__mem%litterfall%fruit%phosphorus%ptr(ics:ice,iblk),[nc]))
allocate(veg_litterfall%fruit%carbon13         ,source=reshape(veg__mem%litterfall%fruit%carbon13%ptr(ics:ice,iblk),[nc]))
allocate(veg_litterfall%fruit%carbon14         ,source=reshape(veg__mem%litterfall%fruit%carbon14%ptr(ics:ice,iblk),[nc]))
allocate(veg_litterfall%fruit%nitrogen15       ,source=reshape(veg__mem%litterfall%fruit%nitrogen15%ptr(ics:ice,iblk),[nc]))
allocate(veg_litterfall%seed_bed%carbon        ,source=reshape(veg__mem%litterfall%seed_bed%carbon%ptr(ics:ice,iblk),[nc]))
allocate(veg_litterfall%seed_bed%nitrogen      ,source=reshape(veg__mem%litterfall%seed_bed%nitrogen%ptr(ics:ice,iblk),[nc]))
allocate(veg_litterfall%seed_bed%phosphorus    ,source=reshape(veg__mem%litterfall%seed_bed%phosphorus%ptr(ics:ice,iblk),[nc]))
allocate(veg_litterfall%seed_bed%carbon13      ,source=reshape(veg__mem%litterfall%seed_bed%carbon13%ptr(ics:ice,iblk),[nc]))
allocate(veg_litterfall%seed_bed%carbon14      ,source=reshape(veg__mem%litterfall%seed_bed%carbon14%ptr(ics:ice,iblk),[nc]))
allocate(veg_litterfall%seed_bed%nitrogen15    ,source=reshape(veg__mem%litterfall%seed_bed%nitrogen15%ptr(ics:ice,iblk),[nc]))

#endif

#if defined(Get_veg_falloc)
allocate(veg_falloc%leaf%carbon            ,source=reshape(veg__mem%falloc%leaf%carbon%ptr(ics:ice,iblk),[nc]))
allocate(veg_falloc%leaf%nitrogen          ,source=reshape(veg__mem%falloc%leaf%nitrogen%ptr(ics:ice,iblk),[nc]))
allocate(veg_falloc%leaf%phosphorus        ,source=reshape(veg__mem%falloc%leaf%phosphorus%ptr(ics:ice,iblk),[nc]))
allocate(veg_falloc%leaf%carbon13          ,source=reshape(veg__mem%falloc%leaf%carbon13%ptr(ics:ice,iblk),[nc]))
allocate(veg_falloc%leaf%carbon14          ,source=reshape(veg__mem%falloc%leaf%carbon14%ptr(ics:ice,iblk),[nc]))
allocate(veg_falloc%leaf%nitrogen15        ,source=reshape(veg__mem%falloc%leaf%nitrogen15%ptr(ics:ice,iblk),[nc]))
allocate(veg_falloc%fine_root%carbon       ,source=reshape(veg__mem%falloc%fine_root%carbon%ptr(ics:ice,iblk),[nc]))
allocate(veg_falloc%fine_root%nitrogen     ,source=reshape(veg__mem%falloc%fine_root%nitrogen%ptr(ics:ice,iblk),[nc]))
allocate(veg_falloc%fine_root%phosphorus   ,source=reshape(veg__mem%falloc%fine_root%phosphorus%ptr(ics:ice,iblk),[nc]))
allocate(veg_falloc%fine_root%carbon13     ,source=reshape(veg__mem%falloc%fine_root%carbon13%ptr(ics:ice,iblk),[nc]))
allocate(veg_falloc%fine_root%carbon14     ,source=reshape(veg__mem%falloc%fine_root%carbon14%ptr(ics:ice,iblk),[nc]))
allocate(veg_falloc%fine_root%nitrogen15   ,source=reshape(veg__mem%falloc%fine_root%nitrogen15%ptr(ics:ice,iblk),[nc]))
allocate(veg_falloc%coarse_root%carbon     ,source=reshape(veg__mem%falloc%coarse_root%carbon%ptr(ics:ice,iblk),[nc]))
allocate(veg_falloc%coarse_root%nitrogen   ,source=reshape(veg__mem%falloc%coarse_root%nitrogen%ptr(ics:ice,iblk),[nc]))
allocate(veg_falloc%coarse_root%phosphorus ,source=reshape(veg__mem%falloc%coarse_root%phosphorus%ptr(ics:ice,iblk),[nc]))
allocate(veg_falloc%coarse_root%carbon13   ,source=reshape(veg__mem%falloc%coarse_root%carbon13%ptr(ics:ice,iblk),[nc]))
allocate(veg_falloc%coarse_root%carbon14   ,source=reshape(veg__mem%falloc%coarse_root%carbon14%ptr(ics:ice,iblk),[nc]))
allocate(veg_falloc%coarse_root%nitrogen15 ,source=reshape(veg__mem%falloc%coarse_root%nitrogen15%ptr(ics:ice,iblk),[nc]))
allocate(veg_falloc%sap_wood%carbon        ,source=reshape(veg__mem%falloc%sap_wood%carbon%ptr(ics:ice,iblk),[nc]))
allocate(veg_falloc%sap_wood%nitrogen      ,source=reshape(veg__mem%falloc%sap_wood%nitrogen%ptr(ics:ice,iblk),[nc]))
allocate(veg_falloc%sap_wood%phosphorus    ,source=reshape(veg__mem%falloc%sap_wood%phosphorus%ptr(ics:ice,iblk),[nc]))
allocate(veg_falloc%sap_wood%carbon13      ,source=reshape(veg__mem%falloc%sap_wood%carbon13%ptr(ics:ice,iblk),[nc]))
allocate(veg_falloc%sap_wood%carbon14      ,source=reshape(veg__mem%falloc%sap_wood%carbon14%ptr(ics:ice,iblk),[nc]))
allocate(veg_falloc%sap_wood%nitrogen15    ,source=reshape(veg__mem%falloc%sap_wood%nitrogen15%ptr(ics:ice,iblk),[nc]))
allocate(veg_falloc%heart_wood%carbon      ,source=reshape(veg__mem%falloc%heart_wood%carbon%ptr(ics:ice,iblk),[nc]))
allocate(veg_falloc%heart_wood%nitrogen    ,source=reshape(veg__mem%falloc%heart_wood%nitrogen%ptr(ics:ice,iblk),[nc]))
allocate(veg_falloc%heart_wood%phosphorus  ,source=reshape(veg__mem%falloc%heart_wood%phosphorus%ptr(ics:ice,iblk),[nc]))
allocate(veg_falloc%heart_wood%carbon13    ,source=reshape(veg__mem%falloc%heart_wood%carbon13%ptr(ics:ice,iblk),[nc]))
allocate(veg_falloc%heart_wood%carbon14    ,source=reshape(veg__mem%falloc%heart_wood%carbon14%ptr(ics:ice,iblk),[nc]))
allocate(veg_falloc%heart_wood%nitrogen15  ,source=reshape(veg__mem%falloc%heart_wood%nitrogen15%ptr(ics:ice,iblk),[nc]))
allocate(veg_falloc%labile%carbon          ,source=reshape(veg__mem%falloc%labile%carbon%ptr(ics:ice,iblk),[nc]))
allocate(veg_falloc%labile%nitrogen        ,source=reshape(veg__mem%falloc%labile%nitrogen%ptr(ics:ice,iblk),[nc]))
allocate(veg_falloc%labile%phosphorus      ,source=reshape(veg__mem%falloc%labile%phosphorus%ptr(ics:ice,iblk),[nc]))
allocate(veg_falloc%labile%carbon13        ,source=reshape(veg__mem%falloc%labile%carbon13%ptr(ics:ice,iblk),[nc]))
allocate(veg_falloc%labile%carbon14        ,source=reshape(veg__mem%falloc%labile%carbon14%ptr(ics:ice,iblk),[nc]))
allocate(veg_falloc%labile%nitrogen15      ,source=reshape(veg__mem%falloc%labile%nitrogen15%ptr(ics:ice,iblk),[nc]))
allocate(veg_falloc%reserve%carbon         ,source=reshape(veg__mem%falloc%reserve%carbon%ptr(ics:ice,iblk),[nc]))
allocate(veg_falloc%reserve%nitrogen       ,source=reshape(veg__mem%falloc%reserve%nitrogen%ptr(ics:ice,iblk),[nc]))
allocate(veg_falloc%reserve%phosphorus     ,source=reshape(veg__mem%falloc%reserve%phosphorus%ptr(ics:ice,iblk),[nc]))
allocate(veg_falloc%reserve%carbon13       ,source=reshape(veg__mem%falloc%reserve%carbon13%ptr(ics:ice,iblk),[nc]))
allocate(veg_falloc%reserve%carbon14       ,source=reshape(veg__mem%falloc%reserve%carbon14%ptr(ics:ice,iblk),[nc]))
allocate(veg_falloc%reserve%nitrogen15     ,source=reshape(veg__mem%falloc%reserve%nitrogen15%ptr(ics:ice,iblk),[nc]))
allocate(veg_falloc%fruit%carbon           ,source=reshape(veg__mem%falloc%fruit%carbon%ptr(ics:ice,iblk),[nc]))
allocate(veg_falloc%fruit%nitrogen         ,source=reshape(veg__mem%falloc%fruit%nitrogen%ptr(ics:ice,iblk),[nc]))
allocate(veg_falloc%fruit%phosphorus       ,source=reshape(veg__mem%falloc%fruit%phosphorus%ptr(ics:ice,iblk),[nc]))
allocate(veg_falloc%fruit%carbon13         ,source=reshape(veg__mem%falloc%fruit%carbon13%ptr(ics:ice,iblk),[nc]))
allocate(veg_falloc%fruit%carbon14         ,source=reshape(veg__mem%falloc%fruit%carbon14%ptr(ics:ice,iblk),[nc]))
allocate(veg_falloc%fruit%nitrogen15       ,source=reshape(veg__mem%falloc%fruit%nitrogen15%ptr(ics:ice,iblk),[nc]))
allocate(veg_falloc%seed_bed%carbon        ,source=reshape(veg__mem%falloc%seed_bed%carbon%ptr(ics:ice,iblk),[nc]))
allocate(veg_falloc%seed_bed%nitrogen      ,source=reshape(veg__mem%falloc%seed_bed%nitrogen%ptr(ics:ice,iblk),[nc]))
allocate(veg_falloc%seed_bed%phosphorus    ,source=reshape(veg__mem%falloc%seed_bed%phosphorus%ptr(ics:ice,iblk),[nc]))
allocate(veg_falloc%seed_bed%carbon13      ,source=reshape(veg__mem%falloc%seed_bed%carbon13%ptr(ics:ice,iblk),[nc]))
allocate(veg_falloc%seed_bed%carbon14      ,source=reshape(veg__mem%falloc%seed_bed%carbon14%ptr(ics:ice,iblk),[nc]))
allocate(veg_falloc%seed_bed%nitrogen15    ,source=reshape(veg__mem%falloc%seed_bed%nitrogen15%ptr(ics:ice,iblk),[nc]))

#endif

#if defined(Get_veg_total_flux)
allocate(veg_total_flux%leaf%carbon            ,source=reshape(veg__mem%total_flux%leaf%carbon%ptr(ics:ice,iblk),[nc]))
allocate(veg_total_flux%leaf%nitrogen          ,source=reshape(veg__mem%total_flux%leaf%nitrogen%ptr(ics:ice,iblk),[nc]))
allocate(veg_total_flux%leaf%phosphorus        ,source=reshape(veg__mem%total_flux%leaf%phosphorus%ptr(ics:ice,iblk),[nc]))
allocate(veg_total_flux%leaf%carbon13          ,source=reshape(veg__mem%total_flux%leaf%carbon13%ptr(ics:ice,iblk),[nc]))
allocate(veg_total_flux%leaf%carbon14          ,source=reshape(veg__mem%total_flux%leaf%carbon14%ptr(ics:ice,iblk),[nc]))
allocate(veg_total_flux%leaf%nitrogen15        ,source=reshape(veg__mem%total_flux%leaf%nitrogen15%ptr(ics:ice,iblk),[nc]))
allocate(veg_total_flux%fine_root%carbon       ,source=reshape(veg__mem%total_flux%fine_root%carbon%ptr(ics:ice,iblk),[nc]))
allocate(veg_total_flux%fine_root%nitrogen     ,source=reshape(veg__mem%total_flux%fine_root%nitrogen%ptr(ics:ice,iblk),[nc]))
allocate(veg_total_flux%fine_root%phosphorus   ,source=reshape(veg__mem%total_flux%fine_root%phosphorus%ptr(ics:ice,iblk),[nc]))
allocate(veg_total_flux%fine_root%carbon13     ,source=reshape(veg__mem%total_flux%fine_root%carbon13%ptr(ics:ice,iblk),[nc]))
allocate(veg_total_flux%fine_root%carbon14     ,source=reshape(veg__mem%total_flux%fine_root%carbon14%ptr(ics:ice,iblk),[nc]))
allocate(veg_total_flux%fine_root%nitrogen15   ,source=reshape(veg__mem%total_flux%fine_root%nitrogen15%ptr(ics:ice,iblk),[nc]))
allocate(veg_total_flux%coarse_root%carbon     ,source=reshape(veg__mem%total_flux%coarse_root%carbon%ptr(ics:ice,iblk),[nc]))
allocate(veg_total_flux%coarse_root%nitrogen   ,source=reshape(veg__mem%total_flux%coarse_root%nitrogen%ptr(ics:ice,iblk),[nc]))
allocate(veg_total_flux%coarse_root%phosphorus ,source=reshape(veg__mem%total_flux%coarse_root%phosphorus%ptr(ics:ice,iblk),[nc]))
allocate(veg_total_flux%coarse_root%carbon13   ,source=reshape(veg__mem%total_flux%coarse_root%carbon13%ptr(ics:ice,iblk),[nc]))
allocate(veg_total_flux%coarse_root%carbon14   ,source=reshape(veg__mem%total_flux%coarse_root%carbon14%ptr(ics:ice,iblk),[nc]))
allocate(veg_total_flux%coarse_root%nitrogen15 ,source=reshape(veg__mem%total_flux%coarse_root%nitrogen15%ptr(ics:ice,iblk),[nc]))
allocate(veg_total_flux%sap_wood%carbon        ,source=reshape(veg__mem%total_flux%sap_wood%carbon%ptr(ics:ice,iblk),[nc]))
allocate(veg_total_flux%sap_wood%nitrogen      ,source=reshape(veg__mem%total_flux%sap_wood%nitrogen%ptr(ics:ice,iblk),[nc]))
allocate(veg_total_flux%sap_wood%phosphorus    ,source=reshape(veg__mem%total_flux%sap_wood%phosphorus%ptr(ics:ice,iblk),[nc]))
allocate(veg_total_flux%sap_wood%carbon13      ,source=reshape(veg__mem%total_flux%sap_wood%carbon13%ptr(ics:ice,iblk),[nc]))
allocate(veg_total_flux%sap_wood%carbon14      ,source=reshape(veg__mem%total_flux%sap_wood%carbon14%ptr(ics:ice,iblk),[nc]))
allocate(veg_total_flux%sap_wood%nitrogen15    ,source=reshape(veg__mem%total_flux%sap_wood%nitrogen15%ptr(ics:ice,iblk),[nc]))
allocate(veg_total_flux%heart_wood%carbon      ,source=reshape(veg__mem%total_flux%heart_wood%carbon%ptr(ics:ice,iblk),[nc]))
allocate(veg_total_flux%heart_wood%nitrogen    ,source=reshape(veg__mem%total_flux%heart_wood%nitrogen%ptr(ics:ice,iblk),[nc]))
allocate(veg_total_flux%heart_wood%phosphorus  ,source=reshape(veg__mem%total_flux%heart_wood%phosphorus%ptr(ics:ice,iblk),[nc]))
allocate(veg_total_flux%heart_wood%carbon13    ,source=reshape(veg__mem%total_flux%heart_wood%carbon13%ptr(ics:ice,iblk),[nc]))
allocate(veg_total_flux%heart_wood%carbon14    ,source=reshape(veg__mem%total_flux%heart_wood%carbon14%ptr(ics:ice,iblk),[nc]))
allocate(veg_total_flux%heart_wood%nitrogen15  ,source=reshape(veg__mem%total_flux%heart_wood%nitrogen15%ptr(ics:ice,iblk),[nc]))
allocate(veg_total_flux%labile%carbon          ,source=reshape(veg__mem%total_flux%labile%carbon%ptr(ics:ice,iblk),[nc]))
allocate(veg_total_flux%labile%nitrogen        ,source=reshape(veg__mem%total_flux%labile%nitrogen%ptr(ics:ice,iblk),[nc]))
allocate(veg_total_flux%labile%phosphorus      ,source=reshape(veg__mem%total_flux%labile%phosphorus%ptr(ics:ice,iblk),[nc]))
allocate(veg_total_flux%labile%carbon13        ,source=reshape(veg__mem%total_flux%labile%carbon13%ptr(ics:ice,iblk),[nc]))
allocate(veg_total_flux%labile%carbon14        ,source=reshape(veg__mem%total_flux%labile%carbon14%ptr(ics:ice,iblk),[nc]))
allocate(veg_total_flux%labile%nitrogen15      ,source=reshape(veg__mem%total_flux%labile%nitrogen15%ptr(ics:ice,iblk),[nc]))
allocate(veg_total_flux%reserve%carbon         ,source=reshape(veg__mem%total_flux%reserve%carbon%ptr(ics:ice,iblk),[nc]))
allocate(veg_total_flux%reserve%nitrogen       ,source=reshape(veg__mem%total_flux%reserve%nitrogen%ptr(ics:ice,iblk),[nc]))
allocate(veg_total_flux%reserve%phosphorus     ,source=reshape(veg__mem%total_flux%reserve%phosphorus%ptr(ics:ice,iblk),[nc]))
allocate(veg_total_flux%reserve%carbon13       ,source=reshape(veg__mem%total_flux%reserve%carbon13%ptr(ics:ice,iblk),[nc]))
allocate(veg_total_flux%reserve%carbon14       ,source=reshape(veg__mem%total_flux%reserve%carbon14%ptr(ics:ice,iblk),[nc]))
allocate(veg_total_flux%reserve%nitrogen15     ,source=reshape(veg__mem%total_flux%reserve%nitrogen15%ptr(ics:ice,iblk),[nc]))
allocate(veg_total_flux%fruit%carbon           ,source=reshape(veg__mem%total_flux%fruit%carbon%ptr(ics:ice,iblk),[nc]))
allocate(veg_total_flux%fruit%nitrogen         ,source=reshape(veg__mem%total_flux%fruit%nitrogen%ptr(ics:ice,iblk),[nc]))
allocate(veg_total_flux%fruit%phosphorus       ,source=reshape(veg__mem%total_flux%fruit%phosphorus%ptr(ics:ice,iblk),[nc]))
allocate(veg_total_flux%fruit%carbon13         ,source=reshape(veg__mem%total_flux%fruit%carbon13%ptr(ics:ice,iblk),[nc]))
allocate(veg_total_flux%fruit%carbon14         ,source=reshape(veg__mem%total_flux%fruit%carbon14%ptr(ics:ice,iblk),[nc]))
allocate(veg_total_flux%fruit%nitrogen15       ,source=reshape(veg__mem%total_flux%fruit%nitrogen15%ptr(ics:ice,iblk),[nc]))
allocate(veg_total_flux%seed_bed%carbon        ,source=reshape(veg__mem%total_flux%seed_bed%carbon%ptr(ics:ice,iblk),[nc]))
allocate(veg_total_flux%seed_bed%nitrogen      ,source=reshape(veg__mem%total_flux%seed_bed%nitrogen%ptr(ics:ice,iblk),[nc]))
allocate(veg_total_flux%seed_bed%phosphorus    ,source=reshape(veg__mem%total_flux%seed_bed%phosphorus%ptr(ics:ice,iblk),[nc]))
allocate(veg_total_flux%seed_bed%carbon13      ,source=reshape(veg__mem%total_flux%seed_bed%carbon13%ptr(ics:ice,iblk),[nc]))
allocate(veg_total_flux%seed_bed%carbon14      ,source=reshape(veg__mem%total_flux%seed_bed%carbon14%ptr(ics:ice,iblk),[nc]))
allocate(veg_total_flux%seed_bed%nitrogen15    ,source=reshape(veg__mem%total_flux%seed_bed%nitrogen15%ptr(ics:ice,iblk),[nc]))

#endif

#ifdef Get_veg_exudation
allocate(veg_exudation%carbon     ,source=reshape(veg__mem%exudation%carbon%ptr(ics:ice,iblk),[nc]))
allocate(veg_exudation%nitrogen   ,source=reshape(veg__mem%exudation%nitrogen%ptr(ics:ice,iblk),[nc]))
allocate(veg_exudation%phosphorus ,source=reshape(veg__mem%exudation%phosphorus%ptr(ics:ice,iblk),[nc]))
allocate(veg_exudation%carbon13   ,source=reshape(veg__mem%exudation%carbon13%ptr(ics:ice,iblk),[nc]))
allocate(veg_exudation%carbon14   ,source=reshape(veg__mem%exudation%carbon14%ptr(ics:ice,iblk),[nc]))
allocate(veg_exudation%nitrogen15 ,source=reshape(veg__mem%exudation%nitrogen15%ptr(ics:ice,iblk),[nc]))

#endif

#ifdef Get_veg_establishment
allocate(veg_establishment%carbon     ,source=reshape(veg__mem%establishment%carbon%ptr(ics:ice,iblk),[nc]))
allocate(veg_establishment%nitrogen   ,source=reshape(veg__mem%establishment%nitrogen%ptr(ics:ice,iblk),[nc]))
allocate(veg_establishment%phosphorus ,source=reshape(veg__mem%establishment%phosphorus%ptr(ics:ice,iblk),[nc]))
allocate(veg_establishment%carbon13   ,source=reshape(veg__mem%establishment%carbon13%ptr(ics:ice,iblk),[nc]))
allocate(veg_establishment%carbon14   ,source=reshape(veg__mem%establishment%carbon14%ptr(ics:ice,iblk),[nc]))
allocate(veg_establishment%nitrogen15 ,source=reshape(veg__mem%establishment%nitrogen15%ptr(ics:ice,iblk),[nc]))

#endif

#ifdef Get_veg_reserve_use
allocate(veg_reserve_use%carbon     ,source=reshape(veg__mem%reserve_use%carbon%ptr(ics:ice,iblk),[nc]))
allocate(veg_reserve_use%nitrogen   ,source=reshape(veg__mem%reserve_use%nitrogen%ptr(ics:ice,iblk),[nc]))
allocate(veg_reserve_use%phosphorus ,source=reshape(veg__mem%reserve_use%phosphorus%ptr(ics:ice,iblk),[nc]))
allocate(veg_reserve_use%carbon13   ,source=reshape(veg__mem%reserve_use%carbon13%ptr(ics:ice,iblk),[nc]))
allocate(veg_reserve_use%carbon14   ,source=reshape(veg__mem%reserve_use%carbon14%ptr(ics:ice,iblk),[nc]))
allocate(veg_reserve_use%nitrogen15 ,source=reshape(veg__mem%reserve_use%nitrogen15%ptr(ics:ice,iblk),[nc]))

#endif

#ifdef Get_veg_total_biomass
allocate(veg_total_biomass%carbon     ,source=reshape(veg__mem%total_biomass%carbon%ptr(ics:ice,iblk),[nc]))
allocate(veg_total_biomass%nitrogen   ,source=reshape(veg__mem%total_biomass%nitrogen%ptr(ics:ice,iblk),[nc]))
allocate(veg_total_biomass%phosphorus ,source=reshape(veg__mem%total_biomass%phosphorus%ptr(ics:ice,iblk),[nc]))
allocate(veg_total_biomass%carbon13   ,source=reshape(veg__mem%total_biomass%carbon13%ptr(ics:ice,iblk),[nc]))
allocate(veg_total_biomass%carbon14   ,source=reshape(veg__mem%total_biomass%carbon14%ptr(ics:ice,iblk),[nc]))
allocate(veg_total_biomass%nitrogen15 ,source=reshape(veg__mem%total_biomass%nitrogen15%ptr(ics:ice,iblk),[nc]))

#endif


! -------------------------------------------------------
! allocate variables in pools of sb memory
! -------------------------------------------------------

#if defined(Get_sb_pool)
allocate(sb_pool%dom%carbon                   ,& 
                    source=reshape(sb__mem%pool%dom%carbon%ptr(ics:ice,1:nsoil_sb,iblk),[nc,nsoil_sb]))
allocate(sb_pool%dom%nitrogen                 ,& 
                    source=reshape(sb__mem%pool%dom%nitrogen%ptr(ics:ice,1:nsoil_sb,iblk),[nc,nsoil_sb]))
allocate(sb_pool%dom%phosphorus               ,& 
                    source=reshape(sb__mem%pool%dom%phosphorus%ptr(ics:ice,1:nsoil_sb,iblk),[nc,nsoil_sb]))
allocate(sb_pool%dom%carbon13                 ,& 
                    source=reshape(sb__mem%pool%dom%carbon13%ptr(ics:ice,1:nsoil_sb,iblk),[nc,nsoil_sb]))
allocate(sb_pool%dom%carbon14                 ,& 
                    source=reshape(sb__mem%pool%dom%carbon14%ptr(ics:ice,1:nsoil_sb,iblk),[nc,nsoil_sb]))
allocate(sb_pool%dom%nitrogen15               ,& 
                    source=reshape(sb__mem%pool%dom%nitrogen15%ptr(ics:ice,1:nsoil_sb,iblk),[nc,nsoil_sb]))
allocate(sb_pool%dom_assoc%carbon             ,& 
                    source=reshape(sb__mem%pool%dom_assoc%carbon%ptr(ics:ice,1:nsoil_sb,iblk),[nc,nsoil_sb]))
allocate(sb_pool%dom_assoc%nitrogen           ,& 
                    source=reshape(sb__mem%pool%dom_assoc%nitrogen%ptr(ics:ice,1:nsoil_sb,iblk),[nc,nsoil_sb]))
allocate(sb_pool%dom_assoc%phosphorus         ,& 
                    source=reshape(sb__mem%pool%dom_assoc%phosphorus%ptr(ics:ice,1:nsoil_sb,iblk),[nc,nsoil_sb]))
allocate(sb_pool%dom_assoc%carbon13           ,& 
                    source=reshape(sb__mem%pool%dom_assoc%carbon13%ptr(ics:ice,1:nsoil_sb,iblk),[nc,nsoil_sb]))
allocate(sb_pool%dom_assoc%carbon14           ,& 
                    source=reshape(sb__mem%pool%dom_assoc%carbon14%ptr(ics:ice,1:nsoil_sb,iblk),[nc,nsoil_sb]))
allocate(sb_pool%dom_assoc%nitrogen15         ,& 
                    source=reshape(sb__mem%pool%dom_assoc%nitrogen15%ptr(ics:ice,1:nsoil_sb,iblk),[nc,nsoil_sb]))
allocate(sb_pool%soluable_litter%carbon       ,& 
                    source=reshape(sb__mem%pool%soluable_litter%carbon%ptr(ics:ice,1:nsoil_sb,iblk),[nc,nsoil_sb]))
allocate(sb_pool%soluable_litter%nitrogen     ,& 
                    source=reshape(sb__mem%pool%soluable_litter%nitrogen%ptr(ics:ice,1:nsoil_sb,iblk),[nc,nsoil_sb]))
allocate(sb_pool%soluable_litter%phosphorus   ,& 
                    source=reshape(sb__mem%pool%soluable_litter%phosphorus%ptr(ics:ice,1:nsoil_sb,iblk),[nc,nsoil_sb]))
allocate(sb_pool%soluable_litter%carbon13     ,& 
                    source=reshape(sb__mem%pool%soluable_litter%carbon13%ptr(ics:ice,1:nsoil_sb,iblk),[nc,nsoil_sb]))
allocate(sb_pool%soluable_litter%carbon14     ,& 
                    source=reshape(sb__mem%pool%soluable_litter%carbon14%ptr(ics:ice,1:nsoil_sb,iblk),[nc,nsoil_sb]))
allocate(sb_pool%soluable_litter%nitrogen15   ,& 
                    source=reshape(sb__mem%pool%soluable_litter%nitrogen15%ptr(ics:ice,1:nsoil_sb,iblk),[nc,nsoil_sb]))
allocate(sb_pool%polymeric_litter%carbon      ,& 
                    source=reshape(sb__mem%pool%polymeric_litter%carbon%ptr(ics:ice,1:nsoil_sb,iblk),[nc,nsoil_sb]))
allocate(sb_pool%polymeric_litter%nitrogen    ,& 
                    source=reshape(sb__mem%pool%polymeric_litter%nitrogen%ptr(ics:ice,1:nsoil_sb,iblk),[nc,nsoil_sb]))
allocate(sb_pool%polymeric_litter%phosphorus  ,& 
                    source=reshape(sb__mem%pool%polymeric_litter%phosphorus%ptr(ics:ice,1:nsoil_sb,iblk),[nc,nsoil_sb]))
allocate(sb_pool%polymeric_litter%carbon13    ,& 
                    source=reshape(sb__mem%pool%polymeric_litter%carbon13%ptr(ics:ice,1:nsoil_sb,iblk),[nc,nsoil_sb]))
allocate(sb_pool%polymeric_litter%carbon14    ,& 
                    source=reshape(sb__mem%pool%polymeric_litter%carbon14%ptr(ics:ice,1:nsoil_sb,iblk),[nc,nsoil_sb]))
allocate(sb_pool%polymeric_litter%nitrogen15  ,& 
                    source=reshape(sb__mem%pool%polymeric_litter%nitrogen15%ptr(ics:ice,1:nsoil_sb,iblk),[nc,nsoil_sb]))
allocate(sb_pool%woody_litter%carbon          ,& 
                    source=reshape(sb__mem%pool%woody_litter%carbon%ptr(ics:ice,1:nsoil_sb,iblk),[nc,nsoil_sb]))
allocate(sb_pool%woody_litter%nitrogen        ,& 
                    source=reshape(sb__mem%pool%woody_litter%nitrogen%ptr(ics:ice,1:nsoil_sb,iblk),[nc,nsoil_sb]))
allocate(sb_pool%woody_litter%phosphorus      ,& 
                    source=reshape(sb__mem%pool%woody_litter%phosphorus%ptr(ics:ice,1:nsoil_sb,iblk),[nc,nsoil_sb]))
allocate(sb_pool%woody_litter%carbon13        ,& 
                    source=reshape(sb__mem%pool%woody_litter%carbon13%ptr(ics:ice,1:nsoil_sb,iblk),[nc,nsoil_sb]))
allocate(sb_pool%woody_litter%carbon14        ,& 
                    source=reshape(sb__mem%pool%woody_litter%carbon14%ptr(ics:ice,1:nsoil_sb,iblk),[nc,nsoil_sb]))
allocate(sb_pool%woody_litter%nitrogen15      ,& 
                    source=reshape(sb__mem%pool%woody_litter%nitrogen15%ptr(ics:ice,1:nsoil_sb,iblk),[nc,nsoil_sb]))
allocate(sb_pool%fungi%carbon                 ,& 
                    source=reshape(sb__mem%pool%fungi%carbon%ptr(ics:ice,1:nsoil_sb,iblk),[nc,nsoil_sb]))
allocate(sb_pool%fungi%nitrogen               ,& 
                    source=reshape(sb__mem%pool%fungi%nitrogen%ptr(ics:ice,1:nsoil_sb,iblk),[nc,nsoil_sb]))
allocate(sb_pool%fungi%phosphorus             ,& 
                    source=reshape(sb__mem%pool%fungi%phosphorus%ptr(ics:ice,1:nsoil_sb,iblk),[nc,nsoil_sb]))
allocate(sb_pool%fungi%carbon13               ,& 
                    source=reshape(sb__mem%pool%fungi%carbon13%ptr(ics:ice,1:nsoil_sb,iblk),[nc,nsoil_sb]))
allocate(sb_pool%fungi%carbon14               ,& 
                    source=reshape(sb__mem%pool%fungi%carbon14%ptr(ics:ice,1:nsoil_sb,iblk),[nc,nsoil_sb]))
allocate(sb_pool%fungi%nitrogen15             ,& 
                    source=reshape(sb__mem%pool%fungi%nitrogen15%ptr(ics:ice,1:nsoil_sb,iblk),[nc,nsoil_sb]))
allocate(sb_pool%mycorrhiza%carbon            ,& 
                    source=reshape(sb__mem%pool%mycorrhiza%carbon%ptr(ics:ice,1:nsoil_sb,iblk),[nc,nsoil_sb]))
allocate(sb_pool%mycorrhiza%nitrogen          ,& 
                    source=reshape(sb__mem%pool%mycorrhiza%nitrogen%ptr(ics:ice,1:nsoil_sb,iblk),[nc,nsoil_sb]))
allocate(sb_pool%mycorrhiza%phosphorus        ,& 
                    source=reshape(sb__mem%pool%mycorrhiza%phosphorus%ptr(ics:ice,1:nsoil_sb,iblk),[nc,nsoil_sb]))
allocate(sb_pool%mycorrhiza%carbon13          ,& 
                    source=reshape(sb__mem%pool%mycorrhiza%carbon13%ptr(ics:ice,1:nsoil_sb,iblk),[nc,nsoil_sb]))
allocate(sb_pool%mycorrhiza%carbon14          ,& 
                    source=reshape(sb__mem%pool%mycorrhiza%carbon14%ptr(ics:ice,1:nsoil_sb,iblk),[nc,nsoil_sb]))
allocate(sb_pool%mycorrhiza%nitrogen15        ,& 
                    source=reshape(sb__mem%pool%mycorrhiza%nitrogen15%ptr(ics:ice,1:nsoil_sb,iblk),[nc,nsoil_sb]))
allocate(sb_pool%microbial%carbon             ,& 
                    source=reshape(sb__mem%pool%microbial%carbon%ptr(ics:ice,1:nsoil_sb,iblk),[nc,nsoil_sb]))
allocate(sb_pool%microbial%nitrogen           ,& 
                    source=reshape(sb__mem%pool%microbial%nitrogen%ptr(ics:ice,1:nsoil_sb,iblk),[nc,nsoil_sb]))
allocate(sb_pool%microbial%phosphorus         ,& 
                    source=reshape(sb__mem%pool%microbial%phosphorus%ptr(ics:ice,1:nsoil_sb,iblk),[nc,nsoil_sb]))
allocate(sb_pool%microbial%carbon13           ,& 
                    source=reshape(sb__mem%pool%microbial%carbon13%ptr(ics:ice,1:nsoil_sb,iblk),[nc,nsoil_sb]))
allocate(sb_pool%microbial%carbon14           ,& 
                    source=reshape(sb__mem%pool%microbial%carbon14%ptr(ics:ice,1:nsoil_sb,iblk),[nc,nsoil_sb]))
allocate(sb_pool%microbial%nitrogen15         ,& 
                    source=reshape(sb__mem%pool%microbial%nitrogen15%ptr(ics:ice,1:nsoil_sb,iblk),[nc,nsoil_sb]))
allocate(sb_pool%residue%carbon               ,& 
                    source=reshape(sb__mem%pool%residue%carbon%ptr(ics:ice,1:nsoil_sb,iblk),[nc,nsoil_sb]))
allocate(sb_pool%residue%nitrogen             ,& 
                    source=reshape(sb__mem%pool%residue%nitrogen%ptr(ics:ice,1:nsoil_sb,iblk),[nc,nsoil_sb]))
allocate(sb_pool%residue%phosphorus           ,& 
                    source=reshape(sb__mem%pool%residue%phosphorus%ptr(ics:ice,1:nsoil_sb,iblk),[nc,nsoil_sb]))
allocate(sb_pool%residue%carbon13             ,& 
                    source=reshape(sb__mem%pool%residue%carbon13%ptr(ics:ice,1:nsoil_sb,iblk),[nc,nsoil_sb]))
allocate(sb_pool%residue%carbon14             ,& 
                    source=reshape(sb__mem%pool%residue%carbon14%ptr(ics:ice,1:nsoil_sb,iblk),[nc,nsoil_sb]))
allocate(sb_pool%residue%nitrogen15           ,& 
                    source=reshape(sb__mem%pool%residue%nitrogen15%ptr(ics:ice,1:nsoil_sb,iblk),[nc,nsoil_sb]))
allocate(sb_pool%residue_assoc%carbon         ,& 
                    source=reshape(sb__mem%pool%residue_assoc%carbon%ptr(ics:ice,1:nsoil_sb,iblk),[nc,nsoil_sb]))
allocate(sb_pool%residue_assoc%nitrogen       ,& 
                    source=reshape(sb__mem%pool%residue_assoc%nitrogen%ptr(ics:ice,1:nsoil_sb,iblk),[nc,nsoil_sb]))
allocate(sb_pool%residue_assoc%phosphorus     ,& 
                    source=reshape(sb__mem%pool%residue_assoc%phosphorus%ptr(ics:ice,1:nsoil_sb,iblk),[nc,nsoil_sb]))
allocate(sb_pool%residue_assoc%carbon13       ,& 
                    source=reshape(sb__mem%pool%residue_assoc%carbon13%ptr(ics:ice,1:nsoil_sb,iblk),[nc,nsoil_sb]))
allocate(sb_pool%residue_assoc%carbon14       ,& 
                    source=reshape(sb__mem%pool%residue_assoc%carbon14%ptr(ics:ice,1:nsoil_sb,iblk),[nc,nsoil_sb]))
allocate(sb_pool%residue_assoc%nitrogen15     ,& 
                    source=reshape(sb__mem%pool%residue_assoc%nitrogen15%ptr(ics:ice,1:nsoil_sb,iblk),[nc,nsoil_sb]))

#endif

#if defined(Get_sb_formation)
allocate(sb_formation%dom%carbon                   ,& 
                    source=reshape(sb__mem%formation%dom%carbon%ptr(ics:ice,1:nsoil_sb,iblk),[nc,nsoil_sb]))
allocate(sb_formation%dom%nitrogen                 ,& 
                    source=reshape(sb__mem%formation%dom%nitrogen%ptr(ics:ice,1:nsoil_sb,iblk),[nc,nsoil_sb]))
allocate(sb_formation%dom%phosphorus               ,& 
                    source=reshape(sb__mem%formation%dom%phosphorus%ptr(ics:ice,1:nsoil_sb,iblk),[nc,nsoil_sb]))
allocate(sb_formation%dom%carbon13                 ,& 
                    source=reshape(sb__mem%formation%dom%carbon13%ptr(ics:ice,1:nsoil_sb,iblk),[nc,nsoil_sb]))
allocate(sb_formation%dom%carbon14                 ,& 
                    source=reshape(sb__mem%formation%dom%carbon14%ptr(ics:ice,1:nsoil_sb,iblk),[nc,nsoil_sb]))
allocate(sb_formation%dom%nitrogen15               ,& 
                    source=reshape(sb__mem%formation%dom%nitrogen15%ptr(ics:ice,1:nsoil_sb,iblk),[nc,nsoil_sb]))
allocate(sb_formation%dom_assoc%carbon             ,& 
                    source=reshape(sb__mem%formation%dom_assoc%carbon%ptr(ics:ice,1:nsoil_sb,iblk),[nc,nsoil_sb]))
allocate(sb_formation%dom_assoc%nitrogen           ,& 
                    source=reshape(sb__mem%formation%dom_assoc%nitrogen%ptr(ics:ice,1:nsoil_sb,iblk),[nc,nsoil_sb]))
allocate(sb_formation%dom_assoc%phosphorus         ,& 
                    source=reshape(sb__mem%formation%dom_assoc%phosphorus%ptr(ics:ice,1:nsoil_sb,iblk),[nc,nsoil_sb]))
allocate(sb_formation%dom_assoc%carbon13           ,& 
                    source=reshape(sb__mem%formation%dom_assoc%carbon13%ptr(ics:ice,1:nsoil_sb,iblk),[nc,nsoil_sb]))
allocate(sb_formation%dom_assoc%carbon14           ,& 
                    source=reshape(sb__mem%formation%dom_assoc%carbon14%ptr(ics:ice,1:nsoil_sb,iblk),[nc,nsoil_sb]))
allocate(sb_formation%dom_assoc%nitrogen15         ,& 
                    source=reshape(sb__mem%formation%dom_assoc%nitrogen15%ptr(ics:ice,1:nsoil_sb,iblk),[nc,nsoil_sb]))
allocate(sb_formation%soluable_litter%carbon       ,& 
                    source=reshape(sb__mem%formation%soluable_litter%carbon%ptr(ics:ice,1:nsoil_sb,iblk),[nc,nsoil_sb]))
allocate(sb_formation%soluable_litter%nitrogen     ,& 
                    source=reshape(sb__mem%formation%soluable_litter%nitrogen%ptr(ics:ice,1:nsoil_sb,iblk),[nc,nsoil_sb]))
allocate(sb_formation%soluable_litter%phosphorus   ,& 
                    source=reshape(sb__mem%formation%soluable_litter%phosphorus%ptr(ics:ice,1:nsoil_sb,iblk),[nc,nsoil_sb]))
allocate(sb_formation%soluable_litter%carbon13     ,& 
                    source=reshape(sb__mem%formation%soluable_litter%carbon13%ptr(ics:ice,1:nsoil_sb,iblk),[nc,nsoil_sb]))
allocate(sb_formation%soluable_litter%carbon14     ,& 
                    source=reshape(sb__mem%formation%soluable_litter%carbon14%ptr(ics:ice,1:nsoil_sb,iblk),[nc,nsoil_sb]))
allocate(sb_formation%soluable_litter%nitrogen15   ,& 
                    source=reshape(sb__mem%formation%soluable_litter%nitrogen15%ptr(ics:ice,1:nsoil_sb,iblk),[nc,nsoil_sb]))
allocate(sb_formation%polymeric_litter%carbon      ,& 
                    source=reshape(sb__mem%formation%polymeric_litter%carbon%ptr(ics:ice,1:nsoil_sb,iblk),[nc,nsoil_sb]))
allocate(sb_formation%polymeric_litter%nitrogen    ,& 
                    source=reshape(sb__mem%formation%polymeric_litter%nitrogen%ptr(ics:ice,1:nsoil_sb,iblk),[nc,nsoil_sb]))
allocate(sb_formation%polymeric_litter%phosphorus  ,& 
                    source=reshape(sb__mem%formation%polymeric_litter%phosphorus%ptr(ics:ice,1:nsoil_sb,iblk),[nc,nsoil_sb]))
allocate(sb_formation%polymeric_litter%carbon13    ,& 
                    source=reshape(sb__mem%formation%polymeric_litter%carbon13%ptr(ics:ice,1:nsoil_sb,iblk),[nc,nsoil_sb]))
allocate(sb_formation%polymeric_litter%carbon14    ,& 
                    source=reshape(sb__mem%formation%polymeric_litter%carbon14%ptr(ics:ice,1:nsoil_sb,iblk),[nc,nsoil_sb]))
allocate(sb_formation%polymeric_litter%nitrogen15  ,& 
                    source=reshape(sb__mem%formation%polymeric_litter%nitrogen15%ptr(ics:ice,1:nsoil_sb,iblk),[nc,nsoil_sb]))
allocate(sb_formation%woody_litter%carbon          ,& 
                    source=reshape(sb__mem%formation%woody_litter%carbon%ptr(ics:ice,1:nsoil_sb,iblk),[nc,nsoil_sb]))
allocate(sb_formation%woody_litter%nitrogen        ,& 
                    source=reshape(sb__mem%formation%woody_litter%nitrogen%ptr(ics:ice,1:nsoil_sb,iblk),[nc,nsoil_sb]))
allocate(sb_formation%woody_litter%phosphorus      ,& 
                    source=reshape(sb__mem%formation%woody_litter%phosphorus%ptr(ics:ice,1:nsoil_sb,iblk),[nc,nsoil_sb]))
allocate(sb_formation%woody_litter%carbon13        ,& 
                    source=reshape(sb__mem%formation%woody_litter%carbon13%ptr(ics:ice,1:nsoil_sb,iblk),[nc,nsoil_sb]))
allocate(sb_formation%woody_litter%carbon14        ,& 
                    source=reshape(sb__mem%formation%woody_litter%carbon14%ptr(ics:ice,1:nsoil_sb,iblk),[nc,nsoil_sb]))
allocate(sb_formation%woody_litter%nitrogen15      ,& 
                    source=reshape(sb__mem%formation%woody_litter%nitrogen15%ptr(ics:ice,1:nsoil_sb,iblk),[nc,nsoil_sb]))
allocate(sb_formation%fungi%carbon                 ,& 
                    source=reshape(sb__mem%formation%fungi%carbon%ptr(ics:ice,1:nsoil_sb,iblk),[nc,nsoil_sb]))
allocate(sb_formation%fungi%nitrogen               ,& 
                    source=reshape(sb__mem%formation%fungi%nitrogen%ptr(ics:ice,1:nsoil_sb,iblk),[nc,nsoil_sb]))
allocate(sb_formation%fungi%phosphorus             ,& 
                    source=reshape(sb__mem%formation%fungi%phosphorus%ptr(ics:ice,1:nsoil_sb,iblk),[nc,nsoil_sb]))
allocate(sb_formation%fungi%carbon13               ,& 
                    source=reshape(sb__mem%formation%fungi%carbon13%ptr(ics:ice,1:nsoil_sb,iblk),[nc,nsoil_sb]))
allocate(sb_formation%fungi%carbon14               ,& 
                    source=reshape(sb__mem%formation%fungi%carbon14%ptr(ics:ice,1:nsoil_sb,iblk),[nc,nsoil_sb]))
allocate(sb_formation%fungi%nitrogen15             ,& 
                    source=reshape(sb__mem%formation%fungi%nitrogen15%ptr(ics:ice,1:nsoil_sb,iblk),[nc,nsoil_sb]))
allocate(sb_formation%mycorrhiza%carbon            ,& 
                    source=reshape(sb__mem%formation%mycorrhiza%carbon%ptr(ics:ice,1:nsoil_sb,iblk),[nc,nsoil_sb]))
allocate(sb_formation%mycorrhiza%nitrogen          ,& 
                    source=reshape(sb__mem%formation%mycorrhiza%nitrogen%ptr(ics:ice,1:nsoil_sb,iblk),[nc,nsoil_sb]))
allocate(sb_formation%mycorrhiza%phosphorus        ,& 
                    source=reshape(sb__mem%formation%mycorrhiza%phosphorus%ptr(ics:ice,1:nsoil_sb,iblk),[nc,nsoil_sb]))
allocate(sb_formation%mycorrhiza%carbon13          ,& 
                    source=reshape(sb__mem%formation%mycorrhiza%carbon13%ptr(ics:ice,1:nsoil_sb,iblk),[nc,nsoil_sb]))
allocate(sb_formation%mycorrhiza%carbon14          ,& 
                    source=reshape(sb__mem%formation%mycorrhiza%carbon14%ptr(ics:ice,1:nsoil_sb,iblk),[nc,nsoil_sb]))
allocate(sb_formation%mycorrhiza%nitrogen15        ,& 
                    source=reshape(sb__mem%formation%mycorrhiza%nitrogen15%ptr(ics:ice,1:nsoil_sb,iblk),[nc,nsoil_sb]))
allocate(sb_formation%microbial%carbon             ,& 
                    source=reshape(sb__mem%formation%microbial%carbon%ptr(ics:ice,1:nsoil_sb,iblk),[nc,nsoil_sb]))
allocate(sb_formation%microbial%nitrogen           ,& 
                    source=reshape(sb__mem%formation%microbial%nitrogen%ptr(ics:ice,1:nsoil_sb,iblk),[nc,nsoil_sb]))
allocate(sb_formation%microbial%phosphorus         ,& 
                    source=reshape(sb__mem%formation%microbial%phosphorus%ptr(ics:ice,1:nsoil_sb,iblk),[nc,nsoil_sb]))
allocate(sb_formation%microbial%carbon13           ,& 
                    source=reshape(sb__mem%formation%microbial%carbon13%ptr(ics:ice,1:nsoil_sb,iblk),[nc,nsoil_sb]))
allocate(sb_formation%microbial%carbon14           ,& 
                    source=reshape(sb__mem%formation%microbial%carbon14%ptr(ics:ice,1:nsoil_sb,iblk),[nc,nsoil_sb]))
allocate(sb_formation%microbial%nitrogen15         ,& 
                    source=reshape(sb__mem%formation%microbial%nitrogen15%ptr(ics:ice,1:nsoil_sb,iblk),[nc,nsoil_sb]))
allocate(sb_formation%residue%carbon               ,& 
                    source=reshape(sb__mem%formation%residue%carbon%ptr(ics:ice,1:nsoil_sb,iblk),[nc,nsoil_sb]))
allocate(sb_formation%residue%nitrogen             ,& 
                    source=reshape(sb__mem%formation%residue%nitrogen%ptr(ics:ice,1:nsoil_sb,iblk),[nc,nsoil_sb]))
allocate(sb_formation%residue%phosphorus           ,& 
                    source=reshape(sb__mem%formation%residue%phosphorus%ptr(ics:ice,1:nsoil_sb,iblk),[nc,nsoil_sb]))
allocate(sb_formation%residue%carbon13             ,& 
                    source=reshape(sb__mem%formation%residue%carbon13%ptr(ics:ice,1:nsoil_sb,iblk),[nc,nsoil_sb]))
allocate(sb_formation%residue%carbon14             ,& 
                    source=reshape(sb__mem%formation%residue%carbon14%ptr(ics:ice,1:nsoil_sb,iblk),[nc,nsoil_sb]))
allocate(sb_formation%residue%nitrogen15           ,& 
                    source=reshape(sb__mem%formation%residue%nitrogen15%ptr(ics:ice,1:nsoil_sb,iblk),[nc,nsoil_sb]))
allocate(sb_formation%residue_assoc%carbon         ,& 
                    source=reshape(sb__mem%formation%residue_assoc%carbon%ptr(ics:ice,1:nsoil_sb,iblk),[nc,nsoil_sb]))
allocate(sb_formation%residue_assoc%nitrogen       ,& 
                    source=reshape(sb__mem%formation%residue_assoc%nitrogen%ptr(ics:ice,1:nsoil_sb,iblk),[nc,nsoil_sb]))
allocate(sb_formation%residue_assoc%phosphorus     ,& 
                    source=reshape(sb__mem%formation%residue_assoc%phosphorus%ptr(ics:ice,1:nsoil_sb,iblk),[nc,nsoil_sb]))
allocate(sb_formation%residue_assoc%carbon13       ,& 
                    source=reshape(sb__mem%formation%residue_assoc%carbon13%ptr(ics:ice,1:nsoil_sb,iblk),[nc,nsoil_sb]))
allocate(sb_formation%residue_assoc%carbon14       ,& 
                    source=reshape(sb__mem%formation%residue_assoc%carbon14%ptr(ics:ice,1:nsoil_sb,iblk),[nc,nsoil_sb]))
allocate(sb_formation%residue_assoc%nitrogen15     ,& 
                    source=reshape(sb__mem%formation%residue_assoc%nitrogen15%ptr(ics:ice,1:nsoil_sb,iblk),[nc,nsoil_sb]))

#endif

#if defined(Get_sb_loss)
allocate(sb_loss%dom%carbon                   ,& 
                    source=reshape(sb__mem%loss%dom%carbon%ptr(ics:ice,1:nsoil_sb,iblk),[nc,nsoil_sb]))
allocate(sb_loss%dom%nitrogen                 ,& 
                    source=reshape(sb__mem%loss%dom%nitrogen%ptr(ics:ice,1:nsoil_sb,iblk),[nc,nsoil_sb]))
allocate(sb_loss%dom%phosphorus               ,& 
                    source=reshape(sb__mem%loss%dom%phosphorus%ptr(ics:ice,1:nsoil_sb,iblk),[nc,nsoil_sb]))
allocate(sb_loss%dom%carbon13                 ,& 
                    source=reshape(sb__mem%loss%dom%carbon13%ptr(ics:ice,1:nsoil_sb,iblk),[nc,nsoil_sb]))
allocate(sb_loss%dom%carbon14                 ,& 
                    source=reshape(sb__mem%loss%dom%carbon14%ptr(ics:ice,1:nsoil_sb,iblk),[nc,nsoil_sb]))
allocate(sb_loss%dom%nitrogen15               ,& 
                    source=reshape(sb__mem%loss%dom%nitrogen15%ptr(ics:ice,1:nsoil_sb,iblk),[nc,nsoil_sb]))
allocate(sb_loss%dom_assoc%carbon             ,& 
                    source=reshape(sb__mem%loss%dom_assoc%carbon%ptr(ics:ice,1:nsoil_sb,iblk),[nc,nsoil_sb]))
allocate(sb_loss%dom_assoc%nitrogen           ,& 
                    source=reshape(sb__mem%loss%dom_assoc%nitrogen%ptr(ics:ice,1:nsoil_sb,iblk),[nc,nsoil_sb]))
allocate(sb_loss%dom_assoc%phosphorus         ,& 
                    source=reshape(sb__mem%loss%dom_assoc%phosphorus%ptr(ics:ice,1:nsoil_sb,iblk),[nc,nsoil_sb]))
allocate(sb_loss%dom_assoc%carbon13           ,& 
                    source=reshape(sb__mem%loss%dom_assoc%carbon13%ptr(ics:ice,1:nsoil_sb,iblk),[nc,nsoil_sb]))
allocate(sb_loss%dom_assoc%carbon14           ,& 
                    source=reshape(sb__mem%loss%dom_assoc%carbon14%ptr(ics:ice,1:nsoil_sb,iblk),[nc,nsoil_sb]))
allocate(sb_loss%dom_assoc%nitrogen15         ,& 
                    source=reshape(sb__mem%loss%dom_assoc%nitrogen15%ptr(ics:ice,1:nsoil_sb,iblk),[nc,nsoil_sb]))
allocate(sb_loss%soluable_litter%carbon       ,& 
                    source=reshape(sb__mem%loss%soluable_litter%carbon%ptr(ics:ice,1:nsoil_sb,iblk),[nc,nsoil_sb]))
allocate(sb_loss%soluable_litter%nitrogen     ,& 
                    source=reshape(sb__mem%loss%soluable_litter%nitrogen%ptr(ics:ice,1:nsoil_sb,iblk),[nc,nsoil_sb]))
allocate(sb_loss%soluable_litter%phosphorus   ,& 
                    source=reshape(sb__mem%loss%soluable_litter%phosphorus%ptr(ics:ice,1:nsoil_sb,iblk),[nc,nsoil_sb]))
allocate(sb_loss%soluable_litter%carbon13     ,& 
                    source=reshape(sb__mem%loss%soluable_litter%carbon13%ptr(ics:ice,1:nsoil_sb,iblk),[nc,nsoil_sb]))
allocate(sb_loss%soluable_litter%carbon14     ,& 
                    source=reshape(sb__mem%loss%soluable_litter%carbon14%ptr(ics:ice,1:nsoil_sb,iblk),[nc,nsoil_sb]))
allocate(sb_loss%soluable_litter%nitrogen15   ,& 
                    source=reshape(sb__mem%loss%soluable_litter%nitrogen15%ptr(ics:ice,1:nsoil_sb,iblk),[nc,nsoil_sb]))
allocate(sb_loss%polymeric_litter%carbon      ,& 
                    source=reshape(sb__mem%loss%polymeric_litter%carbon%ptr(ics:ice,1:nsoil_sb,iblk),[nc,nsoil_sb]))
allocate(sb_loss%polymeric_litter%nitrogen    ,& 
                    source=reshape(sb__mem%loss%polymeric_litter%nitrogen%ptr(ics:ice,1:nsoil_sb,iblk),[nc,nsoil_sb]))
allocate(sb_loss%polymeric_litter%phosphorus  ,& 
                    source=reshape(sb__mem%loss%polymeric_litter%phosphorus%ptr(ics:ice,1:nsoil_sb,iblk),[nc,nsoil_sb]))
allocate(sb_loss%polymeric_litter%carbon13    ,& 
                    source=reshape(sb__mem%loss%polymeric_litter%carbon13%ptr(ics:ice,1:nsoil_sb,iblk),[nc,nsoil_sb]))
allocate(sb_loss%polymeric_litter%carbon14    ,& 
                    source=reshape(sb__mem%loss%polymeric_litter%carbon14%ptr(ics:ice,1:nsoil_sb,iblk),[nc,nsoil_sb]))
allocate(sb_loss%polymeric_litter%nitrogen15  ,& 
                    source=reshape(sb__mem%loss%polymeric_litter%nitrogen15%ptr(ics:ice,1:nsoil_sb,iblk),[nc,nsoil_sb]))
allocate(sb_loss%woody_litter%carbon          ,& 
                    source=reshape(sb__mem%loss%woody_litter%carbon%ptr(ics:ice,1:nsoil_sb,iblk),[nc,nsoil_sb]))
allocate(sb_loss%woody_litter%nitrogen        ,& 
                    source=reshape(sb__mem%loss%woody_litter%nitrogen%ptr(ics:ice,1:nsoil_sb,iblk),[nc,nsoil_sb]))
allocate(sb_loss%woody_litter%phosphorus      ,& 
                    source=reshape(sb__mem%loss%woody_litter%phosphorus%ptr(ics:ice,1:nsoil_sb,iblk),[nc,nsoil_sb]))
allocate(sb_loss%woody_litter%carbon13        ,& 
                    source=reshape(sb__mem%loss%woody_litter%carbon13%ptr(ics:ice,1:nsoil_sb,iblk),[nc,nsoil_sb]))
allocate(sb_loss%woody_litter%carbon14        ,& 
                    source=reshape(sb__mem%loss%woody_litter%carbon14%ptr(ics:ice,1:nsoil_sb,iblk),[nc,nsoil_sb]))
allocate(sb_loss%woody_litter%nitrogen15      ,& 
                    source=reshape(sb__mem%loss%woody_litter%nitrogen15%ptr(ics:ice,1:nsoil_sb,iblk),[nc,nsoil_sb]))
allocate(sb_loss%fungi%carbon                 ,& 
                    source=reshape(sb__mem%loss%fungi%carbon%ptr(ics:ice,1:nsoil_sb,iblk),[nc,nsoil_sb]))
allocate(sb_loss%fungi%nitrogen               ,& 
                    source=reshape(sb__mem%loss%fungi%nitrogen%ptr(ics:ice,1:nsoil_sb,iblk),[nc,nsoil_sb]))
allocate(sb_loss%fungi%phosphorus             ,& 
                    source=reshape(sb__mem%loss%fungi%phosphorus%ptr(ics:ice,1:nsoil_sb,iblk),[nc,nsoil_sb]))
allocate(sb_loss%fungi%carbon13               ,& 
                    source=reshape(sb__mem%loss%fungi%carbon13%ptr(ics:ice,1:nsoil_sb,iblk),[nc,nsoil_sb]))
allocate(sb_loss%fungi%carbon14               ,& 
                    source=reshape(sb__mem%loss%fungi%carbon14%ptr(ics:ice,1:nsoil_sb,iblk),[nc,nsoil_sb]))
allocate(sb_loss%fungi%nitrogen15             ,& 
                    source=reshape(sb__mem%loss%fungi%nitrogen15%ptr(ics:ice,1:nsoil_sb,iblk),[nc,nsoil_sb]))
allocate(sb_loss%mycorrhiza%carbon            ,& 
                    source=reshape(sb__mem%loss%mycorrhiza%carbon%ptr(ics:ice,1:nsoil_sb,iblk),[nc,nsoil_sb]))
allocate(sb_loss%mycorrhiza%nitrogen          ,& 
                    source=reshape(sb__mem%loss%mycorrhiza%nitrogen%ptr(ics:ice,1:nsoil_sb,iblk),[nc,nsoil_sb]))
allocate(sb_loss%mycorrhiza%phosphorus        ,& 
                    source=reshape(sb__mem%loss%mycorrhiza%phosphorus%ptr(ics:ice,1:nsoil_sb,iblk),[nc,nsoil_sb]))
allocate(sb_loss%mycorrhiza%carbon13          ,& 
                    source=reshape(sb__mem%loss%mycorrhiza%carbon13%ptr(ics:ice,1:nsoil_sb,iblk),[nc,nsoil_sb]))
allocate(sb_loss%mycorrhiza%carbon14          ,& 
                    source=reshape(sb__mem%loss%mycorrhiza%carbon14%ptr(ics:ice,1:nsoil_sb,iblk),[nc,nsoil_sb]))
allocate(sb_loss%mycorrhiza%nitrogen15        ,& 
                    source=reshape(sb__mem%loss%mycorrhiza%nitrogen15%ptr(ics:ice,1:nsoil_sb,iblk),[nc,nsoil_sb]))
allocate(sb_loss%microbial%carbon             ,& 
                    source=reshape(sb__mem%loss%microbial%carbon%ptr(ics:ice,1:nsoil_sb,iblk),[nc,nsoil_sb]))
allocate(sb_loss%microbial%nitrogen           ,& 
                    source=reshape(sb__mem%loss%microbial%nitrogen%ptr(ics:ice,1:nsoil_sb,iblk),[nc,nsoil_sb]))
allocate(sb_loss%microbial%phosphorus         ,& 
                    source=reshape(sb__mem%loss%microbial%phosphorus%ptr(ics:ice,1:nsoil_sb,iblk),[nc,nsoil_sb]))
allocate(sb_loss%microbial%carbon13           ,& 
                    source=reshape(sb__mem%loss%microbial%carbon13%ptr(ics:ice,1:nsoil_sb,iblk),[nc,nsoil_sb]))
allocate(sb_loss%microbial%carbon14           ,& 
                    source=reshape(sb__mem%loss%microbial%carbon14%ptr(ics:ice,1:nsoil_sb,iblk),[nc,nsoil_sb]))
allocate(sb_loss%microbial%nitrogen15         ,& 
                    source=reshape(sb__mem%loss%microbial%nitrogen15%ptr(ics:ice,1:nsoil_sb,iblk),[nc,nsoil_sb]))
allocate(sb_loss%residue%carbon               ,& 
                    source=reshape(sb__mem%loss%residue%carbon%ptr(ics:ice,1:nsoil_sb,iblk),[nc,nsoil_sb]))
allocate(sb_loss%residue%nitrogen             ,& 
                    source=reshape(sb__mem%loss%residue%nitrogen%ptr(ics:ice,1:nsoil_sb,iblk),[nc,nsoil_sb]))
allocate(sb_loss%residue%phosphorus           ,& 
                    source=reshape(sb__mem%loss%residue%phosphorus%ptr(ics:ice,1:nsoil_sb,iblk),[nc,nsoil_sb]))
allocate(sb_loss%residue%carbon13             ,& 
                    source=reshape(sb__mem%loss%residue%carbon13%ptr(ics:ice,1:nsoil_sb,iblk),[nc,nsoil_sb]))
allocate(sb_loss%residue%carbon14             ,& 
                    source=reshape(sb__mem%loss%residue%carbon14%ptr(ics:ice,1:nsoil_sb,iblk),[nc,nsoil_sb]))
allocate(sb_loss%residue%nitrogen15           ,& 
                    source=reshape(sb__mem%loss%residue%nitrogen15%ptr(ics:ice,1:nsoil_sb,iblk),[nc,nsoil_sb]))
allocate(sb_loss%residue_assoc%carbon         ,& 
                    source=reshape(sb__mem%loss%residue_assoc%carbon%ptr(ics:ice,1:nsoil_sb,iblk),[nc,nsoil_sb]))
allocate(sb_loss%residue_assoc%nitrogen       ,& 
                    source=reshape(sb__mem%loss%residue_assoc%nitrogen%ptr(ics:ice,1:nsoil_sb,iblk),[nc,nsoil_sb]))
allocate(sb_loss%residue_assoc%phosphorus     ,& 
                    source=reshape(sb__mem%loss%residue_assoc%phosphorus%ptr(ics:ice,1:nsoil_sb,iblk),[nc,nsoil_sb]))
allocate(sb_loss%residue_assoc%carbon13       ,& 
                    source=reshape(sb__mem%loss%residue_assoc%carbon13%ptr(ics:ice,1:nsoil_sb,iblk),[nc,nsoil_sb]))
allocate(sb_loss%residue_assoc%carbon14       ,& 
                    source=reshape(sb__mem%loss%residue_assoc%carbon14%ptr(ics:ice,1:nsoil_sb,iblk),[nc,nsoil_sb]))
allocate(sb_loss%residue_assoc%nitrogen15     ,& 
                    source=reshape(sb__mem%loss%residue_assoc%nitrogen15%ptr(ics:ice,1:nsoil_sb,iblk),[nc,nsoil_sb]))

#endif

#if defined(Get_sb_transport)
allocate(sb_transport%dom%carbon                   ,& 
                    source=reshape(sb__mem%transport%dom%carbon%ptr(ics:ice,1:nsoil_sb,iblk),[nc,nsoil_sb]))
allocate(sb_transport%dom%nitrogen                 ,& 
                    source=reshape(sb__mem%transport%dom%nitrogen%ptr(ics:ice,1:nsoil_sb,iblk),[nc,nsoil_sb]))
allocate(sb_transport%dom%phosphorus               ,& 
                    source=reshape(sb__mem%transport%dom%phosphorus%ptr(ics:ice,1:nsoil_sb,iblk),[nc,nsoil_sb]))
allocate(sb_transport%dom%carbon13                 ,& 
                    source=reshape(sb__mem%transport%dom%carbon13%ptr(ics:ice,1:nsoil_sb,iblk),[nc,nsoil_sb]))
allocate(sb_transport%dom%carbon14                 ,& 
                    source=reshape(sb__mem%transport%dom%carbon14%ptr(ics:ice,1:nsoil_sb,iblk),[nc,nsoil_sb]))
allocate(sb_transport%dom%nitrogen15               ,& 
                    source=reshape(sb__mem%transport%dom%nitrogen15%ptr(ics:ice,1:nsoil_sb,iblk),[nc,nsoil_sb]))
allocate(sb_transport%dom_assoc%carbon             ,& 
                    source=reshape(sb__mem%transport%dom_assoc%carbon%ptr(ics:ice,1:nsoil_sb,iblk),[nc,nsoil_sb]))
allocate(sb_transport%dom_assoc%nitrogen           ,& 
                    source=reshape(sb__mem%transport%dom_assoc%nitrogen%ptr(ics:ice,1:nsoil_sb,iblk),[nc,nsoil_sb]))
allocate(sb_transport%dom_assoc%phosphorus         ,& 
                    source=reshape(sb__mem%transport%dom_assoc%phosphorus%ptr(ics:ice,1:nsoil_sb,iblk),[nc,nsoil_sb]))
allocate(sb_transport%dom_assoc%carbon13           ,& 
                    source=reshape(sb__mem%transport%dom_assoc%carbon13%ptr(ics:ice,1:nsoil_sb,iblk),[nc,nsoil_sb]))
allocate(sb_transport%dom_assoc%carbon14           ,& 
                    source=reshape(sb__mem%transport%dom_assoc%carbon14%ptr(ics:ice,1:nsoil_sb,iblk),[nc,nsoil_sb]))
allocate(sb_transport%dom_assoc%nitrogen15         ,& 
                    source=reshape(sb__mem%transport%dom_assoc%nitrogen15%ptr(ics:ice,1:nsoil_sb,iblk),[nc,nsoil_sb]))
allocate(sb_transport%soluable_litter%carbon       ,& 
                    source=reshape(sb__mem%transport%soluable_litter%carbon%ptr(ics:ice,1:nsoil_sb,iblk),[nc,nsoil_sb]))
allocate(sb_transport%soluable_litter%nitrogen     ,& 
                    source=reshape(sb__mem%transport%soluable_litter%nitrogen%ptr(ics:ice,1:nsoil_sb,iblk),[nc,nsoil_sb]))
allocate(sb_transport%soluable_litter%phosphorus   ,& 
                    source=reshape(sb__mem%transport%soluable_litter%phosphorus%ptr(ics:ice,1:nsoil_sb,iblk),[nc,nsoil_sb]))
allocate(sb_transport%soluable_litter%carbon13     ,& 
                    source=reshape(sb__mem%transport%soluable_litter%carbon13%ptr(ics:ice,1:nsoil_sb,iblk),[nc,nsoil_sb]))
allocate(sb_transport%soluable_litter%carbon14     ,& 
                    source=reshape(sb__mem%transport%soluable_litter%carbon14%ptr(ics:ice,1:nsoil_sb,iblk),[nc,nsoil_sb]))
allocate(sb_transport%soluable_litter%nitrogen15   ,& 
                    source=reshape(sb__mem%transport%soluable_litter%nitrogen15%ptr(ics:ice,1:nsoil_sb,iblk),[nc,nsoil_sb]))
allocate(sb_transport%polymeric_litter%carbon      ,& 
                    source=reshape(sb__mem%transport%polymeric_litter%carbon%ptr(ics:ice,1:nsoil_sb,iblk),[nc,nsoil_sb]))
allocate(sb_transport%polymeric_litter%nitrogen    ,& 
                    source=reshape(sb__mem%transport%polymeric_litter%nitrogen%ptr(ics:ice,1:nsoil_sb,iblk),[nc,nsoil_sb]))
allocate(sb_transport%polymeric_litter%phosphorus  ,& 
                    source=reshape(sb__mem%transport%polymeric_litter%phosphorus%ptr(ics:ice,1:nsoil_sb,iblk),[nc,nsoil_sb]))
allocate(sb_transport%polymeric_litter%carbon13    ,& 
                    source=reshape(sb__mem%transport%polymeric_litter%carbon13%ptr(ics:ice,1:nsoil_sb,iblk),[nc,nsoil_sb]))
allocate(sb_transport%polymeric_litter%carbon14    ,& 
                    source=reshape(sb__mem%transport%polymeric_litter%carbon14%ptr(ics:ice,1:nsoil_sb,iblk),[nc,nsoil_sb]))
allocate(sb_transport%polymeric_litter%nitrogen15  ,& 
                    source=reshape(sb__mem%transport%polymeric_litter%nitrogen15%ptr(ics:ice,1:nsoil_sb,iblk),[nc,nsoil_sb]))
allocate(sb_transport%woody_litter%carbon          ,& 
                    source=reshape(sb__mem%transport%woody_litter%carbon%ptr(ics:ice,1:nsoil_sb,iblk),[nc,nsoil_sb]))
allocate(sb_transport%woody_litter%nitrogen        ,& 
                    source=reshape(sb__mem%transport%woody_litter%nitrogen%ptr(ics:ice,1:nsoil_sb,iblk),[nc,nsoil_sb]))
allocate(sb_transport%woody_litter%phosphorus      ,& 
                    source=reshape(sb__mem%transport%woody_litter%phosphorus%ptr(ics:ice,1:nsoil_sb,iblk),[nc,nsoil_sb]))
allocate(sb_transport%woody_litter%carbon13        ,& 
                    source=reshape(sb__mem%transport%woody_litter%carbon13%ptr(ics:ice,1:nsoil_sb,iblk),[nc,nsoil_sb]))
allocate(sb_transport%woody_litter%carbon14        ,& 
                    source=reshape(sb__mem%transport%woody_litter%carbon14%ptr(ics:ice,1:nsoil_sb,iblk),[nc,nsoil_sb]))
allocate(sb_transport%woody_litter%nitrogen15      ,& 
                    source=reshape(sb__mem%transport%woody_litter%nitrogen15%ptr(ics:ice,1:nsoil_sb,iblk),[nc,nsoil_sb]))
allocate(sb_transport%fungi%carbon                 ,& 
                    source=reshape(sb__mem%transport%fungi%carbon%ptr(ics:ice,1:nsoil_sb,iblk),[nc,nsoil_sb]))
allocate(sb_transport%fungi%nitrogen               ,& 
                    source=reshape(sb__mem%transport%fungi%nitrogen%ptr(ics:ice,1:nsoil_sb,iblk),[nc,nsoil_sb]))
allocate(sb_transport%fungi%phosphorus             ,& 
                    source=reshape(sb__mem%transport%fungi%phosphorus%ptr(ics:ice,1:nsoil_sb,iblk),[nc,nsoil_sb]))
allocate(sb_transport%fungi%carbon13               ,& 
                    source=reshape(sb__mem%transport%fungi%carbon13%ptr(ics:ice,1:nsoil_sb,iblk),[nc,nsoil_sb]))
allocate(sb_transport%fungi%carbon14               ,& 
                    source=reshape(sb__mem%transport%fungi%carbon14%ptr(ics:ice,1:nsoil_sb,iblk),[nc,nsoil_sb]))
allocate(sb_transport%fungi%nitrogen15             ,& 
                    source=reshape(sb__mem%transport%fungi%nitrogen15%ptr(ics:ice,1:nsoil_sb,iblk),[nc,nsoil_sb]))
allocate(sb_transport%mycorrhiza%carbon            ,& 
                    source=reshape(sb__mem%transport%mycorrhiza%carbon%ptr(ics:ice,1:nsoil_sb,iblk),[nc,nsoil_sb]))
allocate(sb_transport%mycorrhiza%nitrogen          ,& 
                    source=reshape(sb__mem%transport%mycorrhiza%nitrogen%ptr(ics:ice,1:nsoil_sb,iblk),[nc,nsoil_sb]))
allocate(sb_transport%mycorrhiza%phosphorus        ,& 
                    source=reshape(sb__mem%transport%mycorrhiza%phosphorus%ptr(ics:ice,1:nsoil_sb,iblk),[nc,nsoil_sb]))
allocate(sb_transport%mycorrhiza%carbon13          ,& 
                    source=reshape(sb__mem%transport%mycorrhiza%carbon13%ptr(ics:ice,1:nsoil_sb,iblk),[nc,nsoil_sb]))
allocate(sb_transport%mycorrhiza%carbon14          ,& 
                    source=reshape(sb__mem%transport%mycorrhiza%carbon14%ptr(ics:ice,1:nsoil_sb,iblk),[nc,nsoil_sb]))
allocate(sb_transport%mycorrhiza%nitrogen15        ,& 
                    source=reshape(sb__mem%transport%mycorrhiza%nitrogen15%ptr(ics:ice,1:nsoil_sb,iblk),[nc,nsoil_sb]))
allocate(sb_transport%microbial%carbon             ,& 
                    source=reshape(sb__mem%transport%microbial%carbon%ptr(ics:ice,1:nsoil_sb,iblk),[nc,nsoil_sb]))
allocate(sb_transport%microbial%nitrogen           ,& 
                    source=reshape(sb__mem%transport%microbial%nitrogen%ptr(ics:ice,1:nsoil_sb,iblk),[nc,nsoil_sb]))
allocate(sb_transport%microbial%phosphorus         ,& 
                    source=reshape(sb__mem%transport%microbial%phosphorus%ptr(ics:ice,1:nsoil_sb,iblk),[nc,nsoil_sb]))
allocate(sb_transport%microbial%carbon13           ,& 
                    source=reshape(sb__mem%transport%microbial%carbon13%ptr(ics:ice,1:nsoil_sb,iblk),[nc,nsoil_sb]))
allocate(sb_transport%microbial%carbon14           ,& 
                    source=reshape(sb__mem%transport%microbial%carbon14%ptr(ics:ice,1:nsoil_sb,iblk),[nc,nsoil_sb]))
allocate(sb_transport%microbial%nitrogen15         ,& 
                    source=reshape(sb__mem%transport%microbial%nitrogen15%ptr(ics:ice,1:nsoil_sb,iblk),[nc,nsoil_sb]))
allocate(sb_transport%residue%carbon               ,& 
                    source=reshape(sb__mem%transport%residue%carbon%ptr(ics:ice,1:nsoil_sb,iblk),[nc,nsoil_sb]))
allocate(sb_transport%residue%nitrogen             ,& 
                    source=reshape(sb__mem%transport%residue%nitrogen%ptr(ics:ice,1:nsoil_sb,iblk),[nc,nsoil_sb]))
allocate(sb_transport%residue%phosphorus           ,& 
                    source=reshape(sb__mem%transport%residue%phosphorus%ptr(ics:ice,1:nsoil_sb,iblk),[nc,nsoil_sb]))
allocate(sb_transport%residue%carbon13             ,& 
                    source=reshape(sb__mem%transport%residue%carbon13%ptr(ics:ice,1:nsoil_sb,iblk),[nc,nsoil_sb]))
allocate(sb_transport%residue%carbon14             ,& 
                    source=reshape(sb__mem%transport%residue%carbon14%ptr(ics:ice,1:nsoil_sb,iblk),[nc,nsoil_sb]))
allocate(sb_transport%residue%nitrogen15           ,& 
                    source=reshape(sb__mem%transport%residue%nitrogen15%ptr(ics:ice,1:nsoil_sb,iblk),[nc,nsoil_sb]))
allocate(sb_transport%residue_assoc%carbon         ,& 
                    source=reshape(sb__mem%transport%residue_assoc%carbon%ptr(ics:ice,1:nsoil_sb,iblk),[nc,nsoil_sb]))
allocate(sb_transport%residue_assoc%nitrogen       ,& 
                    source=reshape(sb__mem%transport%residue_assoc%nitrogen%ptr(ics:ice,1:nsoil_sb,iblk),[nc,nsoil_sb]))
allocate(sb_transport%residue_assoc%phosphorus     ,& 
                    source=reshape(sb__mem%transport%residue_assoc%phosphorus%ptr(ics:ice,1:nsoil_sb,iblk),[nc,nsoil_sb]))
allocate(sb_transport%residue_assoc%carbon13       ,& 
                    source=reshape(sb__mem%transport%residue_assoc%carbon13%ptr(ics:ice,1:nsoil_sb,iblk),[nc,nsoil_sb]))
allocate(sb_transport%residue_assoc%carbon14       ,& 
                    source=reshape(sb__mem%transport%residue_assoc%carbon14%ptr(ics:ice,1:nsoil_sb,iblk),[nc,nsoil_sb]))
allocate(sb_transport%residue_assoc%nitrogen15     ,& 
                    source=reshape(sb__mem%transport%residue_assoc%nitrogen15%ptr(ics:ice,1:nsoil_sb,iblk),[nc,nsoil_sb]))

#endif

#if defined(Get_sb_delta_pools)
allocate(sb_delta_pools%carbon      ,& 
                    source=reshape(sb__mem%delta_pools%carbon%ptr(ics:ice,1:nsoil_sb,iblk),[nc,nsoil_sb]))
allocate(sb_delta_pools%nitrogen    ,& 
                    source=reshape(sb__mem%delta_pools%nitrogen%ptr(ics:ice,1:nsoil_sb,iblk),[nc,nsoil_sb]))
allocate(sb_delta_pools%phosphorus  ,& 
                    source=reshape(sb__mem%delta_pools%phosphorus%ptr(ics:ice,1:nsoil_sb,iblk),[nc,nsoil_sb]))
allocate(sb_delta_pools%carbon13    ,& 
                    source=reshape(sb__mem%delta_pools%carbon13%ptr(ics:ice,1:nsoil_sb,iblk),[nc,nsoil_sb]))
allocate(sb_delta_pools%carbon14    ,& 
                    source=reshape(sb__mem%delta_pools%carbon14%ptr(ics:ice,1:nsoil_sb,iblk),[nc,nsoil_sb]))
allocate(sb_delta_pools%nitrogen15  ,& 
                    source=reshape(sb__mem%delta_pools%nitrogen15%ptr(ics:ice,1:nsoil_sb,iblk),[nc,nsoil_sb]))

#endif

#if defined(Get_sb_mycorrhiza_export)
allocate(sb_mycorrhiza_export%carbon      ,& 
                    source=reshape(sb__mem%mycorrhiza_export%carbon%ptr(ics:ice,1:nsoil_sb,iblk),[nc,nsoil_sb]))
allocate(sb_mycorrhiza_export%nitrogen    ,& 
                    source=reshape(sb__mem%mycorrhiza_export%nitrogen%ptr(ics:ice,1:nsoil_sb,iblk),[nc,nsoil_sb]))
allocate(sb_mycorrhiza_export%phosphorus  ,& 
                    source=reshape(sb__mem%mycorrhiza_export%phosphorus%ptr(ics:ice,1:nsoil_sb,iblk),[nc,nsoil_sb]))
allocate(sb_mycorrhiza_export%carbon13    ,& 
                    source=reshape(sb__mem%mycorrhiza_export%carbon13%ptr(ics:ice,1:nsoil_sb,iblk),[nc,nsoil_sb]))
allocate(sb_mycorrhiza_export%carbon14    ,& 
                    source=reshape(sb__mem%mycorrhiza_export%carbon14%ptr(ics:ice,1:nsoil_sb,iblk),[nc,nsoil_sb]))
allocate(sb_mycorrhiza_export%nitrogen15  ,& 
                    source=reshape(sb__mem%mycorrhiza_export%nitrogen15%ptr(ics:ice,1:nsoil_sb,iblk),[nc,nsoil_sb]))

#endif


! -------------------------------------------------------
! deallocate helper variables for vegetation pools
! and pass values to pools of vegetation memory
! -------------------------------------------------------

#if defined(Passback_veg_pool)
  ! pass the values (back) to process memory
  veg__mem%pool%leaf%carbon%ptr(ics:ice,iblk)             = veg_pool%leaf%carbon(:)
  veg__mem%pool%leaf%nitrogen%ptr(ics:ice,iblk)           = veg_pool%leaf%nitrogen(:)
  veg__mem%pool%leaf%phosphorus%ptr(ics:ice,iblk)         = veg_pool%leaf%phosphorus(:)
  veg__mem%pool%leaf%carbon13%ptr(ics:ice,iblk)           = veg_pool%leaf%carbon13(:)
  veg__mem%pool%leaf%carbon14%ptr(ics:ice,iblk)           = veg_pool%leaf%carbon14(:)
  veg__mem%pool%leaf%nitrogen15%ptr(ics:ice,iblk)         = veg_pool%leaf%nitrogen15(:)
  veg__mem%pool%fine_root%carbon%ptr(ics:ice,iblk)        = veg_pool%fine_root%carbon(:)
  veg__mem%pool%fine_root%nitrogen%ptr(ics:ice,iblk)      = veg_pool%fine_root%nitrogen(:)
  veg__mem%pool%fine_root%phosphorus%ptr(ics:ice,iblk)    = veg_pool%fine_root%phosphorus(:)
  veg__mem%pool%fine_root%carbon13%ptr(ics:ice,iblk)      = veg_pool%fine_root%carbon13(:)
  veg__mem%pool%fine_root%carbon14%ptr(ics:ice,iblk)      = veg_pool%fine_root%carbon14(:)
  veg__mem%pool%fine_root%nitrogen15%ptr(ics:ice,iblk)    = veg_pool%fine_root%nitrogen15(:)
  veg__mem%pool%coarse_root%carbon%ptr(ics:ice,iblk)      = veg_pool%coarse_root%carbon(:)
  veg__mem%pool%coarse_root%nitrogen%ptr(ics:ice,iblk)    = veg_pool%coarse_root%nitrogen(:)
  veg__mem%pool%coarse_root%phosphorus%ptr(ics:ice,iblk)  = veg_pool%coarse_root%phosphorus(:)
  veg__mem%pool%coarse_root%carbon13%ptr(ics:ice,iblk)    = veg_pool%coarse_root%carbon13(:)
  veg__mem%pool%coarse_root%carbon14%ptr(ics:ice,iblk)    = veg_pool%coarse_root%carbon14(:)
  veg__mem%pool%coarse_root%nitrogen15%ptr(ics:ice,iblk)  = veg_pool%coarse_root%nitrogen15(:)
  veg__mem%pool%sap_wood%carbon%ptr(ics:ice,iblk)         = veg_pool%sap_wood%carbon(:)
  veg__mem%pool%sap_wood%nitrogen%ptr(ics:ice,iblk)       = veg_pool%sap_wood%nitrogen(:)
  veg__mem%pool%sap_wood%phosphorus%ptr(ics:ice,iblk)     = veg_pool%sap_wood%phosphorus(:)
  veg__mem%pool%sap_wood%carbon13%ptr(ics:ice,iblk)       = veg_pool%sap_wood%carbon13(:)
  veg__mem%pool%sap_wood%carbon14%ptr(ics:ice,iblk)       = veg_pool%sap_wood%carbon14(:)
  veg__mem%pool%sap_wood%nitrogen15%ptr(ics:ice,iblk)     = veg_pool%sap_wood%nitrogen15(:)
  veg__mem%pool%heart_wood%carbon%ptr(ics:ice,iblk)       = veg_pool%heart_wood%carbon(:)
  veg__mem%pool%heart_wood%nitrogen%ptr(ics:ice,iblk)     = veg_pool%heart_wood%nitrogen(:)
  veg__mem%pool%heart_wood%phosphorus%ptr(ics:ice,iblk)   = veg_pool%heart_wood%phosphorus(:)
  veg__mem%pool%heart_wood%carbon13%ptr(ics:ice,iblk)     = veg_pool%heart_wood%carbon13(:)
  veg__mem%pool%heart_wood%carbon14%ptr(ics:ice,iblk)     = veg_pool%heart_wood%carbon14(:)
  veg__mem%pool%heart_wood%nitrogen15%ptr(ics:ice,iblk)   = veg_pool%heart_wood%nitrogen15(:)
  veg__mem%pool%labile%carbon%ptr(ics:ice,iblk)           = veg_pool%labile%carbon(:)
  veg__mem%pool%labile%nitrogen%ptr(ics:ice,iblk)         = veg_pool%labile%nitrogen(:)
  veg__mem%pool%labile%phosphorus%ptr(ics:ice,iblk)       = veg_pool%labile%phosphorus(:)
  veg__mem%pool%labile%carbon13%ptr(ics:ice,iblk)         = veg_pool%labile%carbon13(:)
  veg__mem%pool%labile%carbon14%ptr(ics:ice,iblk)         = veg_pool%labile%carbon14(:)
  veg__mem%pool%labile%nitrogen15%ptr(ics:ice,iblk)       = veg_pool%labile%nitrogen15(:)
  veg__mem%pool%reserve%carbon%ptr(ics:ice,iblk)          = veg_pool%reserve%carbon(:)
  veg__mem%pool%reserve%nitrogen%ptr(ics:ice,iblk)        = veg_pool%reserve%nitrogen(:)
  veg__mem%pool%reserve%phosphorus%ptr(ics:ice,iblk)      = veg_pool%reserve%phosphorus(:)
  veg__mem%pool%reserve%carbon13%ptr(ics:ice,iblk)        = veg_pool%reserve%carbon13(:)
  veg__mem%pool%reserve%carbon14%ptr(ics:ice,iblk)        = veg_pool%reserve%carbon14(:)
  veg__mem%pool%reserve%nitrogen15%ptr(ics:ice,iblk)      = veg_pool%reserve%nitrogen15(:)
  veg__mem%pool%fruit%carbon%ptr(ics:ice,iblk)            = veg_pool%fruit%carbon(:)
  veg__mem%pool%fruit%nitrogen%ptr(ics:ice,iblk)          = veg_pool%fruit%nitrogen(:)
  veg__mem%pool%fruit%phosphorus%ptr(ics:ice,iblk)        = veg_pool%fruit%phosphorus(:)
  veg__mem%pool%fruit%carbon13%ptr(ics:ice,iblk)          = veg_pool%fruit%carbon13(:)
  veg__mem%pool%fruit%carbon14%ptr(ics:ice,iblk)          = veg_pool%fruit%carbon14(:)
  veg__mem%pool%fruit%nitrogen15%ptr(ics:ice,iblk)        = veg_pool%fruit%nitrogen15(:)
  veg__mem%pool%seed_bed%carbon%ptr(ics:ice,iblk)         = veg_pool%seed_bed%carbon(:)
  veg__mem%pool%seed_bed%nitrogen%ptr(ics:ice,iblk)       = veg_pool%seed_bed%nitrogen(:)
  veg__mem%pool%seed_bed%phosphorus%ptr(ics:ice,iblk)     = veg_pool%seed_bed%phosphorus(:)
  veg__mem%pool%seed_bed%carbon13%ptr(ics:ice,iblk)       = veg_pool%seed_bed%carbon13(:)
  veg__mem%pool%seed_bed%carbon14%ptr(ics:ice,iblk)       = veg_pool%seed_bed%carbon14(:)
  veg__mem%pool%seed_bed%nitrogen15%ptr(ics:ice,iblk)     = veg_pool%seed_bed%nitrogen15(:)
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

#if defined(Passback_veg_growth)
  ! pass the values (back) to process memory
  veg__mem%growth%leaf%carbon%ptr(ics:ice,iblk)             = veg_growth%leaf%carbon(:)
  veg__mem%growth%leaf%nitrogen%ptr(ics:ice,iblk)           = veg_growth%leaf%nitrogen(:)
  veg__mem%growth%leaf%phosphorus%ptr(ics:ice,iblk)         = veg_growth%leaf%phosphorus(:)
  veg__mem%growth%leaf%carbon13%ptr(ics:ice,iblk)           = veg_growth%leaf%carbon13(:)
  veg__mem%growth%leaf%carbon14%ptr(ics:ice,iblk)           = veg_growth%leaf%carbon14(:)
  veg__mem%growth%leaf%nitrogen15%ptr(ics:ice,iblk)         = veg_growth%leaf%nitrogen15(:)
  veg__mem%growth%fine_root%carbon%ptr(ics:ice,iblk)        = veg_growth%fine_root%carbon(:)
  veg__mem%growth%fine_root%nitrogen%ptr(ics:ice,iblk)      = veg_growth%fine_root%nitrogen(:)
  veg__mem%growth%fine_root%phosphorus%ptr(ics:ice,iblk)    = veg_growth%fine_root%phosphorus(:)
  veg__mem%growth%fine_root%carbon13%ptr(ics:ice,iblk)      = veg_growth%fine_root%carbon13(:)
  veg__mem%growth%fine_root%carbon14%ptr(ics:ice,iblk)      = veg_growth%fine_root%carbon14(:)
  veg__mem%growth%fine_root%nitrogen15%ptr(ics:ice,iblk)    = veg_growth%fine_root%nitrogen15(:)
  veg__mem%growth%coarse_root%carbon%ptr(ics:ice,iblk)      = veg_growth%coarse_root%carbon(:)
  veg__mem%growth%coarse_root%nitrogen%ptr(ics:ice,iblk)    = veg_growth%coarse_root%nitrogen(:)
  veg__mem%growth%coarse_root%phosphorus%ptr(ics:ice,iblk)  = veg_growth%coarse_root%phosphorus(:)
  veg__mem%growth%coarse_root%carbon13%ptr(ics:ice,iblk)    = veg_growth%coarse_root%carbon13(:)
  veg__mem%growth%coarse_root%carbon14%ptr(ics:ice,iblk)    = veg_growth%coarse_root%carbon14(:)
  veg__mem%growth%coarse_root%nitrogen15%ptr(ics:ice,iblk)  = veg_growth%coarse_root%nitrogen15(:)
  veg__mem%growth%sap_wood%carbon%ptr(ics:ice,iblk)         = veg_growth%sap_wood%carbon(:)
  veg__mem%growth%sap_wood%nitrogen%ptr(ics:ice,iblk)       = veg_growth%sap_wood%nitrogen(:)
  veg__mem%growth%sap_wood%phosphorus%ptr(ics:ice,iblk)     = veg_growth%sap_wood%phosphorus(:)
  veg__mem%growth%sap_wood%carbon13%ptr(ics:ice,iblk)       = veg_growth%sap_wood%carbon13(:)
  veg__mem%growth%sap_wood%carbon14%ptr(ics:ice,iblk)       = veg_growth%sap_wood%carbon14(:)
  veg__mem%growth%sap_wood%nitrogen15%ptr(ics:ice,iblk)     = veg_growth%sap_wood%nitrogen15(:)
  veg__mem%growth%heart_wood%carbon%ptr(ics:ice,iblk)       = veg_growth%heart_wood%carbon(:)
  veg__mem%growth%heart_wood%nitrogen%ptr(ics:ice,iblk)     = veg_growth%heart_wood%nitrogen(:)
  veg__mem%growth%heart_wood%phosphorus%ptr(ics:ice,iblk)   = veg_growth%heart_wood%phosphorus(:)
  veg__mem%growth%heart_wood%carbon13%ptr(ics:ice,iblk)     = veg_growth%heart_wood%carbon13(:)
  veg__mem%growth%heart_wood%carbon14%ptr(ics:ice,iblk)     = veg_growth%heart_wood%carbon14(:)
  veg__mem%growth%heart_wood%nitrogen15%ptr(ics:ice,iblk)   = veg_growth%heart_wood%nitrogen15(:)
  veg__mem%growth%labile%carbon%ptr(ics:ice,iblk)           = veg_growth%labile%carbon(:)
  veg__mem%growth%labile%nitrogen%ptr(ics:ice,iblk)         = veg_growth%labile%nitrogen(:)
  veg__mem%growth%labile%phosphorus%ptr(ics:ice,iblk)       = veg_growth%labile%phosphorus(:)
  veg__mem%growth%labile%carbon13%ptr(ics:ice,iblk)         = veg_growth%labile%carbon13(:)
  veg__mem%growth%labile%carbon14%ptr(ics:ice,iblk)         = veg_growth%labile%carbon14(:)
  veg__mem%growth%labile%nitrogen15%ptr(ics:ice,iblk)       = veg_growth%labile%nitrogen15(:)
  veg__mem%growth%reserve%carbon%ptr(ics:ice,iblk)          = veg_growth%reserve%carbon(:)
  veg__mem%growth%reserve%nitrogen%ptr(ics:ice,iblk)        = veg_growth%reserve%nitrogen(:)
  veg__mem%growth%reserve%phosphorus%ptr(ics:ice,iblk)      = veg_growth%reserve%phosphorus(:)
  veg__mem%growth%reserve%carbon13%ptr(ics:ice,iblk)        = veg_growth%reserve%carbon13(:)
  veg__mem%growth%reserve%carbon14%ptr(ics:ice,iblk)        = veg_growth%reserve%carbon14(:)
  veg__mem%growth%reserve%nitrogen15%ptr(ics:ice,iblk)      = veg_growth%reserve%nitrogen15(:)
  veg__mem%growth%fruit%carbon%ptr(ics:ice,iblk)            = veg_growth%fruit%carbon(:)
  veg__mem%growth%fruit%nitrogen%ptr(ics:ice,iblk)          = veg_growth%fruit%nitrogen(:)
  veg__mem%growth%fruit%phosphorus%ptr(ics:ice,iblk)        = veg_growth%fruit%phosphorus(:)
  veg__mem%growth%fruit%carbon13%ptr(ics:ice,iblk)          = veg_growth%fruit%carbon13(:)
  veg__mem%growth%fruit%carbon14%ptr(ics:ice,iblk)          = veg_growth%fruit%carbon14(:)
  veg__mem%growth%fruit%nitrogen15%ptr(ics:ice,iblk)        = veg_growth%fruit%nitrogen15(:)
  veg__mem%growth%seed_bed%carbon%ptr(ics:ice,iblk)         = veg_growth%seed_bed%carbon(:)
  veg__mem%growth%seed_bed%nitrogen%ptr(ics:ice,iblk)       = veg_growth%seed_bed%nitrogen(:)
  veg__mem%growth%seed_bed%phosphorus%ptr(ics:ice,iblk)     = veg_growth%seed_bed%phosphorus(:)
  veg__mem%growth%seed_bed%carbon13%ptr(ics:ice,iblk)       = veg_growth%seed_bed%carbon13(:)
  veg__mem%growth%seed_bed%carbon14%ptr(ics:ice,iblk)       = veg_growth%seed_bed%carbon14(:)
  veg__mem%growth%seed_bed%nitrogen15%ptr(ics:ice,iblk)     = veg_growth%seed_bed%nitrogen15(:)
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

#if defined(Passback_veg_litterfall)
  ! pass the values (back) to process memory
  veg__mem%litterfall%leaf%carbon%ptr(ics:ice,iblk)             = veg_litterfall%leaf%carbon(:)
  veg__mem%litterfall%leaf%nitrogen%ptr(ics:ice,iblk)           = veg_litterfall%leaf%nitrogen(:)
  veg__mem%litterfall%leaf%phosphorus%ptr(ics:ice,iblk)         = veg_litterfall%leaf%phosphorus(:)
  veg__mem%litterfall%leaf%carbon13%ptr(ics:ice,iblk)           = veg_litterfall%leaf%carbon13(:)
  veg__mem%litterfall%leaf%carbon14%ptr(ics:ice,iblk)           = veg_litterfall%leaf%carbon14(:)
  veg__mem%litterfall%leaf%nitrogen15%ptr(ics:ice,iblk)         = veg_litterfall%leaf%nitrogen15(:)
  veg__mem%litterfall%fine_root%carbon%ptr(ics:ice,iblk)        = veg_litterfall%fine_root%carbon(:)
  veg__mem%litterfall%fine_root%nitrogen%ptr(ics:ice,iblk)      = veg_litterfall%fine_root%nitrogen(:)
  veg__mem%litterfall%fine_root%phosphorus%ptr(ics:ice,iblk)    = veg_litterfall%fine_root%phosphorus(:)
  veg__mem%litterfall%fine_root%carbon13%ptr(ics:ice,iblk)      = veg_litterfall%fine_root%carbon13(:)
  veg__mem%litterfall%fine_root%carbon14%ptr(ics:ice,iblk)      = veg_litterfall%fine_root%carbon14(:)
  veg__mem%litterfall%fine_root%nitrogen15%ptr(ics:ice,iblk)    = veg_litterfall%fine_root%nitrogen15(:)
  veg__mem%litterfall%coarse_root%carbon%ptr(ics:ice,iblk)      = veg_litterfall%coarse_root%carbon(:)
  veg__mem%litterfall%coarse_root%nitrogen%ptr(ics:ice,iblk)    = veg_litterfall%coarse_root%nitrogen(:)
  veg__mem%litterfall%coarse_root%phosphorus%ptr(ics:ice,iblk)  = veg_litterfall%coarse_root%phosphorus(:)
  veg__mem%litterfall%coarse_root%carbon13%ptr(ics:ice,iblk)    = veg_litterfall%coarse_root%carbon13(:)
  veg__mem%litterfall%coarse_root%carbon14%ptr(ics:ice,iblk)    = veg_litterfall%coarse_root%carbon14(:)
  veg__mem%litterfall%coarse_root%nitrogen15%ptr(ics:ice,iblk)  = veg_litterfall%coarse_root%nitrogen15(:)
  veg__mem%litterfall%sap_wood%carbon%ptr(ics:ice,iblk)         = veg_litterfall%sap_wood%carbon(:)
  veg__mem%litterfall%sap_wood%nitrogen%ptr(ics:ice,iblk)       = veg_litterfall%sap_wood%nitrogen(:)
  veg__mem%litterfall%sap_wood%phosphorus%ptr(ics:ice,iblk)     = veg_litterfall%sap_wood%phosphorus(:)
  veg__mem%litterfall%sap_wood%carbon13%ptr(ics:ice,iblk)       = veg_litterfall%sap_wood%carbon13(:)
  veg__mem%litterfall%sap_wood%carbon14%ptr(ics:ice,iblk)       = veg_litterfall%sap_wood%carbon14(:)
  veg__mem%litterfall%sap_wood%nitrogen15%ptr(ics:ice,iblk)     = veg_litterfall%sap_wood%nitrogen15(:)
  veg__mem%litterfall%heart_wood%carbon%ptr(ics:ice,iblk)       = veg_litterfall%heart_wood%carbon(:)
  veg__mem%litterfall%heart_wood%nitrogen%ptr(ics:ice,iblk)     = veg_litterfall%heart_wood%nitrogen(:)
  veg__mem%litterfall%heart_wood%phosphorus%ptr(ics:ice,iblk)   = veg_litterfall%heart_wood%phosphorus(:)
  veg__mem%litterfall%heart_wood%carbon13%ptr(ics:ice,iblk)     = veg_litterfall%heart_wood%carbon13(:)
  veg__mem%litterfall%heart_wood%carbon14%ptr(ics:ice,iblk)     = veg_litterfall%heart_wood%carbon14(:)
  veg__mem%litterfall%heart_wood%nitrogen15%ptr(ics:ice,iblk)   = veg_litterfall%heart_wood%nitrogen15(:)
  veg__mem%litterfall%labile%carbon%ptr(ics:ice,iblk)           = veg_litterfall%labile%carbon(:)
  veg__mem%litterfall%labile%nitrogen%ptr(ics:ice,iblk)         = veg_litterfall%labile%nitrogen(:)
  veg__mem%litterfall%labile%phosphorus%ptr(ics:ice,iblk)       = veg_litterfall%labile%phosphorus(:)
  veg__mem%litterfall%labile%carbon13%ptr(ics:ice,iblk)         = veg_litterfall%labile%carbon13(:)
  veg__mem%litterfall%labile%carbon14%ptr(ics:ice,iblk)         = veg_litterfall%labile%carbon14(:)
  veg__mem%litterfall%labile%nitrogen15%ptr(ics:ice,iblk)       = veg_litterfall%labile%nitrogen15(:)
  veg__mem%litterfall%reserve%carbon%ptr(ics:ice,iblk)          = veg_litterfall%reserve%carbon(:)
  veg__mem%litterfall%reserve%nitrogen%ptr(ics:ice,iblk)        = veg_litterfall%reserve%nitrogen(:)
  veg__mem%litterfall%reserve%phosphorus%ptr(ics:ice,iblk)      = veg_litterfall%reserve%phosphorus(:)
  veg__mem%litterfall%reserve%carbon13%ptr(ics:ice,iblk)        = veg_litterfall%reserve%carbon13(:)
  veg__mem%litterfall%reserve%carbon14%ptr(ics:ice,iblk)        = veg_litterfall%reserve%carbon14(:)
  veg__mem%litterfall%reserve%nitrogen15%ptr(ics:ice,iblk)      = veg_litterfall%reserve%nitrogen15(:)
  veg__mem%litterfall%fruit%carbon%ptr(ics:ice,iblk)            = veg_litterfall%fruit%carbon(:)
  veg__mem%litterfall%fruit%nitrogen%ptr(ics:ice,iblk)          = veg_litterfall%fruit%nitrogen(:)
  veg__mem%litterfall%fruit%phosphorus%ptr(ics:ice,iblk)        = veg_litterfall%fruit%phosphorus(:)
  veg__mem%litterfall%fruit%carbon13%ptr(ics:ice,iblk)          = veg_litterfall%fruit%carbon13(:)
  veg__mem%litterfall%fruit%carbon14%ptr(ics:ice,iblk)          = veg_litterfall%fruit%carbon14(:)
  veg__mem%litterfall%fruit%nitrogen15%ptr(ics:ice,iblk)        = veg_litterfall%fruit%nitrogen15(:)
  veg__mem%litterfall%seed_bed%carbon%ptr(ics:ice,iblk)         = veg_litterfall%seed_bed%carbon(:)
  veg__mem%litterfall%seed_bed%nitrogen%ptr(ics:ice,iblk)       = veg_litterfall%seed_bed%nitrogen(:)
  veg__mem%litterfall%seed_bed%phosphorus%ptr(ics:ice,iblk)     = veg_litterfall%seed_bed%phosphorus(:)
  veg__mem%litterfall%seed_bed%carbon13%ptr(ics:ice,iblk)       = veg_litterfall%seed_bed%carbon13(:)
  veg__mem%litterfall%seed_bed%carbon14%ptr(ics:ice,iblk)       = veg_litterfall%seed_bed%carbon14(:)
  veg__mem%litterfall%seed_bed%nitrogen15%ptr(ics:ice,iblk)     = veg_litterfall%seed_bed%nitrogen15(:)
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

#if defined(Passback_veg_falloc)
  ! pass the values (back) to process memory
  veg__mem%falloc%leaf%carbon%ptr(ics:ice,iblk)             = veg_falloc%leaf%carbon(:)
  veg__mem%falloc%leaf%nitrogen%ptr(ics:ice,iblk)           = veg_falloc%leaf%nitrogen(:)
  veg__mem%falloc%leaf%phosphorus%ptr(ics:ice,iblk)         = veg_falloc%leaf%phosphorus(:)
  veg__mem%falloc%leaf%carbon13%ptr(ics:ice,iblk)           = veg_falloc%leaf%carbon13(:)
  veg__mem%falloc%leaf%carbon14%ptr(ics:ice,iblk)           = veg_falloc%leaf%carbon14(:)
  veg__mem%falloc%leaf%nitrogen15%ptr(ics:ice,iblk)         = veg_falloc%leaf%nitrogen15(:)
  veg__mem%falloc%fine_root%carbon%ptr(ics:ice,iblk)        = veg_falloc%fine_root%carbon(:)
  veg__mem%falloc%fine_root%nitrogen%ptr(ics:ice,iblk)      = veg_falloc%fine_root%nitrogen(:)
  veg__mem%falloc%fine_root%phosphorus%ptr(ics:ice,iblk)    = veg_falloc%fine_root%phosphorus(:)
  veg__mem%falloc%fine_root%carbon13%ptr(ics:ice,iblk)      = veg_falloc%fine_root%carbon13(:)
  veg__mem%falloc%fine_root%carbon14%ptr(ics:ice,iblk)      = veg_falloc%fine_root%carbon14(:)
  veg__mem%falloc%fine_root%nitrogen15%ptr(ics:ice,iblk)    = veg_falloc%fine_root%nitrogen15(:)
  veg__mem%falloc%coarse_root%carbon%ptr(ics:ice,iblk)      = veg_falloc%coarse_root%carbon(:)
  veg__mem%falloc%coarse_root%nitrogen%ptr(ics:ice,iblk)    = veg_falloc%coarse_root%nitrogen(:)
  veg__mem%falloc%coarse_root%phosphorus%ptr(ics:ice,iblk)  = veg_falloc%coarse_root%phosphorus(:)
  veg__mem%falloc%coarse_root%carbon13%ptr(ics:ice,iblk)    = veg_falloc%coarse_root%carbon13(:)
  veg__mem%falloc%coarse_root%carbon14%ptr(ics:ice,iblk)    = veg_falloc%coarse_root%carbon14(:)
  veg__mem%falloc%coarse_root%nitrogen15%ptr(ics:ice,iblk)  = veg_falloc%coarse_root%nitrogen15(:)
  veg__mem%falloc%sap_wood%carbon%ptr(ics:ice,iblk)         = veg_falloc%sap_wood%carbon(:)
  veg__mem%falloc%sap_wood%nitrogen%ptr(ics:ice,iblk)       = veg_falloc%sap_wood%nitrogen(:)
  veg__mem%falloc%sap_wood%phosphorus%ptr(ics:ice,iblk)     = veg_falloc%sap_wood%phosphorus(:)
  veg__mem%falloc%sap_wood%carbon13%ptr(ics:ice,iblk)       = veg_falloc%sap_wood%carbon13(:)
  veg__mem%falloc%sap_wood%carbon14%ptr(ics:ice,iblk)       = veg_falloc%sap_wood%carbon14(:)
  veg__mem%falloc%sap_wood%nitrogen15%ptr(ics:ice,iblk)     = veg_falloc%sap_wood%nitrogen15(:)
  veg__mem%falloc%heart_wood%carbon%ptr(ics:ice,iblk)       = veg_falloc%heart_wood%carbon(:)
  veg__mem%falloc%heart_wood%nitrogen%ptr(ics:ice,iblk)     = veg_falloc%heart_wood%nitrogen(:)
  veg__mem%falloc%heart_wood%phosphorus%ptr(ics:ice,iblk)   = veg_falloc%heart_wood%phosphorus(:)
  veg__mem%falloc%heart_wood%carbon13%ptr(ics:ice,iblk)     = veg_falloc%heart_wood%carbon13(:)
  veg__mem%falloc%heart_wood%carbon14%ptr(ics:ice,iblk)     = veg_falloc%heart_wood%carbon14(:)
  veg__mem%falloc%heart_wood%nitrogen15%ptr(ics:ice,iblk)   = veg_falloc%heart_wood%nitrogen15(:)
  veg__mem%falloc%labile%carbon%ptr(ics:ice,iblk)           = veg_falloc%labile%carbon(:)
  veg__mem%falloc%labile%nitrogen%ptr(ics:ice,iblk)         = veg_falloc%labile%nitrogen(:)
  veg__mem%falloc%labile%phosphorus%ptr(ics:ice,iblk)       = veg_falloc%labile%phosphorus(:)
  veg__mem%falloc%labile%carbon13%ptr(ics:ice,iblk)         = veg_falloc%labile%carbon13(:)
  veg__mem%falloc%labile%carbon14%ptr(ics:ice,iblk)         = veg_falloc%labile%carbon14(:)
  veg__mem%falloc%labile%nitrogen15%ptr(ics:ice,iblk)       = veg_falloc%labile%nitrogen15(:)
  veg__mem%falloc%reserve%carbon%ptr(ics:ice,iblk)          = veg_falloc%reserve%carbon(:)
  veg__mem%falloc%reserve%nitrogen%ptr(ics:ice,iblk)        = veg_falloc%reserve%nitrogen(:)
  veg__mem%falloc%reserve%phosphorus%ptr(ics:ice,iblk)      = veg_falloc%reserve%phosphorus(:)
  veg__mem%falloc%reserve%carbon13%ptr(ics:ice,iblk)        = veg_falloc%reserve%carbon13(:)
  veg__mem%falloc%reserve%carbon14%ptr(ics:ice,iblk)        = veg_falloc%reserve%carbon14(:)
  veg__mem%falloc%reserve%nitrogen15%ptr(ics:ice,iblk)      = veg_falloc%reserve%nitrogen15(:)
  veg__mem%falloc%fruit%carbon%ptr(ics:ice,iblk)            = veg_falloc%fruit%carbon(:)
  veg__mem%falloc%fruit%nitrogen%ptr(ics:ice,iblk)          = veg_falloc%fruit%nitrogen(:)
  veg__mem%falloc%fruit%phosphorus%ptr(ics:ice,iblk)        = veg_falloc%fruit%phosphorus(:)
  veg__mem%falloc%fruit%carbon13%ptr(ics:ice,iblk)          = veg_falloc%fruit%carbon13(:)
  veg__mem%falloc%fruit%carbon14%ptr(ics:ice,iblk)          = veg_falloc%fruit%carbon14(:)
  veg__mem%falloc%fruit%nitrogen15%ptr(ics:ice,iblk)        = veg_falloc%fruit%nitrogen15(:)
  veg__mem%falloc%seed_bed%carbon%ptr(ics:ice,iblk)         = veg_falloc%seed_bed%carbon(:)
  veg__mem%falloc%seed_bed%nitrogen%ptr(ics:ice,iblk)       = veg_falloc%seed_bed%nitrogen(:)
  veg__mem%falloc%seed_bed%phosphorus%ptr(ics:ice,iblk)     = veg_falloc%seed_bed%phosphorus(:)
  veg__mem%falloc%seed_bed%carbon13%ptr(ics:ice,iblk)       = veg_falloc%seed_bed%carbon13(:)
  veg__mem%falloc%seed_bed%carbon14%ptr(ics:ice,iblk)       = veg_falloc%seed_bed%carbon14(:)
  veg__mem%falloc%seed_bed%nitrogen15%ptr(ics:ice,iblk)     = veg_falloc%seed_bed%nitrogen15(:)
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

#if defined(Passback_veg_total_flux)
  ! pass the values (back) to process memory
  veg__mem%total_flux%leaf%carbon%ptr(ics:ice,iblk)             = veg_total_flux%leaf%carbon(:)
  veg__mem%total_flux%leaf%nitrogen%ptr(ics:ice,iblk)           = veg_total_flux%leaf%nitrogen(:)
  veg__mem%total_flux%leaf%phosphorus%ptr(ics:ice,iblk)         = veg_total_flux%leaf%phosphorus(:)
  veg__mem%total_flux%leaf%carbon13%ptr(ics:ice,iblk)           = veg_total_flux%leaf%carbon13(:)
  veg__mem%total_flux%leaf%carbon14%ptr(ics:ice,iblk)           = veg_total_flux%leaf%carbon14(:)
  veg__mem%total_flux%leaf%nitrogen15%ptr(ics:ice,iblk)         = veg_total_flux%leaf%nitrogen15(:)
  veg__mem%total_flux%fine_root%carbon%ptr(ics:ice,iblk)        = veg_total_flux%fine_root%carbon(:)
  veg__mem%total_flux%fine_root%nitrogen%ptr(ics:ice,iblk)      = veg_total_flux%fine_root%nitrogen(:)
  veg__mem%total_flux%fine_root%phosphorus%ptr(ics:ice,iblk)    = veg_total_flux%fine_root%phosphorus(:)
  veg__mem%total_flux%fine_root%carbon13%ptr(ics:ice,iblk)      = veg_total_flux%fine_root%carbon13(:)
  veg__mem%total_flux%fine_root%carbon14%ptr(ics:ice,iblk)      = veg_total_flux%fine_root%carbon14(:)
  veg__mem%total_flux%fine_root%nitrogen15%ptr(ics:ice,iblk)    = veg_total_flux%fine_root%nitrogen15(:)
  veg__mem%total_flux%coarse_root%carbon%ptr(ics:ice,iblk)      = veg_total_flux%coarse_root%carbon(:)
  veg__mem%total_flux%coarse_root%nitrogen%ptr(ics:ice,iblk)    = veg_total_flux%coarse_root%nitrogen(:)
  veg__mem%total_flux%coarse_root%phosphorus%ptr(ics:ice,iblk)  = veg_total_flux%coarse_root%phosphorus(:)
  veg__mem%total_flux%coarse_root%carbon13%ptr(ics:ice,iblk)    = veg_total_flux%coarse_root%carbon13(:)
  veg__mem%total_flux%coarse_root%carbon14%ptr(ics:ice,iblk)    = veg_total_flux%coarse_root%carbon14(:)
  veg__mem%total_flux%coarse_root%nitrogen15%ptr(ics:ice,iblk)  = veg_total_flux%coarse_root%nitrogen15(:)
  veg__mem%total_flux%sap_wood%carbon%ptr(ics:ice,iblk)         = veg_total_flux%sap_wood%carbon(:)
  veg__mem%total_flux%sap_wood%nitrogen%ptr(ics:ice,iblk)       = veg_total_flux%sap_wood%nitrogen(:)
  veg__mem%total_flux%sap_wood%phosphorus%ptr(ics:ice,iblk)     = veg_total_flux%sap_wood%phosphorus(:)
  veg__mem%total_flux%sap_wood%carbon13%ptr(ics:ice,iblk)       = veg_total_flux%sap_wood%carbon13(:)
  veg__mem%total_flux%sap_wood%carbon14%ptr(ics:ice,iblk)       = veg_total_flux%sap_wood%carbon14(:)
  veg__mem%total_flux%sap_wood%nitrogen15%ptr(ics:ice,iblk)     = veg_total_flux%sap_wood%nitrogen15(:)
  veg__mem%total_flux%heart_wood%carbon%ptr(ics:ice,iblk)       = veg_total_flux%heart_wood%carbon(:)
  veg__mem%total_flux%heart_wood%nitrogen%ptr(ics:ice,iblk)     = veg_total_flux%heart_wood%nitrogen(:)
  veg__mem%total_flux%heart_wood%phosphorus%ptr(ics:ice,iblk)   = veg_total_flux%heart_wood%phosphorus(:)
  veg__mem%total_flux%heart_wood%carbon13%ptr(ics:ice,iblk)     = veg_total_flux%heart_wood%carbon13(:)
  veg__mem%total_flux%heart_wood%carbon14%ptr(ics:ice,iblk)     = veg_total_flux%heart_wood%carbon14(:)
  veg__mem%total_flux%heart_wood%nitrogen15%ptr(ics:ice,iblk)   = veg_total_flux%heart_wood%nitrogen15(:)
  veg__mem%total_flux%labile%carbon%ptr(ics:ice,iblk)           = veg_total_flux%labile%carbon(:)
  veg__mem%total_flux%labile%nitrogen%ptr(ics:ice,iblk)         = veg_total_flux%labile%nitrogen(:)
  veg__mem%total_flux%labile%phosphorus%ptr(ics:ice,iblk)       = veg_total_flux%labile%phosphorus(:)
  veg__mem%total_flux%labile%carbon13%ptr(ics:ice,iblk)         = veg_total_flux%labile%carbon13(:)
  veg__mem%total_flux%labile%carbon14%ptr(ics:ice,iblk)         = veg_total_flux%labile%carbon14(:)
  veg__mem%total_flux%labile%nitrogen15%ptr(ics:ice,iblk)       = veg_total_flux%labile%nitrogen15(:)
  veg__mem%total_flux%reserve%carbon%ptr(ics:ice,iblk)          = veg_total_flux%reserve%carbon(:)
  veg__mem%total_flux%reserve%nitrogen%ptr(ics:ice,iblk)        = veg_total_flux%reserve%nitrogen(:)
  veg__mem%total_flux%reserve%phosphorus%ptr(ics:ice,iblk)      = veg_total_flux%reserve%phosphorus(:)
  veg__mem%total_flux%reserve%carbon13%ptr(ics:ice,iblk)        = veg_total_flux%reserve%carbon13(:)
  veg__mem%total_flux%reserve%carbon14%ptr(ics:ice,iblk)        = veg_total_flux%reserve%carbon14(:)
  veg__mem%total_flux%reserve%nitrogen15%ptr(ics:ice,iblk)      = veg_total_flux%reserve%nitrogen15(:)
  veg__mem%total_flux%fruit%carbon%ptr(ics:ice,iblk)            = veg_total_flux%fruit%carbon(:)
  veg__mem%total_flux%fruit%nitrogen%ptr(ics:ice,iblk)          = veg_total_flux%fruit%nitrogen(:)
  veg__mem%total_flux%fruit%phosphorus%ptr(ics:ice,iblk)        = veg_total_flux%fruit%phosphorus(:)
  veg__mem%total_flux%fruit%carbon13%ptr(ics:ice,iblk)          = veg_total_flux%fruit%carbon13(:)
  veg__mem%total_flux%fruit%carbon14%ptr(ics:ice,iblk)          = veg_total_flux%fruit%carbon14(:)
  veg__mem%total_flux%fruit%nitrogen15%ptr(ics:ice,iblk)        = veg_total_flux%fruit%nitrogen15(:)
  veg__mem%total_flux%seed_bed%carbon%ptr(ics:ice,iblk)         = veg_total_flux%seed_bed%carbon(:)
  veg__mem%total_flux%seed_bed%nitrogen%ptr(ics:ice,iblk)       = veg_total_flux%seed_bed%nitrogen(:)
  veg__mem%total_flux%seed_bed%phosphorus%ptr(ics:ice,iblk)     = veg_total_flux%seed_bed%phosphorus(:)
  veg__mem%total_flux%seed_bed%carbon13%ptr(ics:ice,iblk)       = veg_total_flux%seed_bed%carbon13(:)
  veg__mem%total_flux%seed_bed%carbon14%ptr(ics:ice,iblk)       = veg_total_flux%seed_bed%carbon14(:)
  veg__mem%total_flux%seed_bed%nitrogen15%ptr(ics:ice,iblk)     = veg_total_flux%seed_bed%nitrogen15(:)
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

#ifdef Passback_veg_exudation
  ! pass the values (back) to process memory
  veg__mem%exudation%carbon%ptr(ics:ice,iblk)             = veg_exudation%carbon(:)
  veg__mem%exudation%nitrogen%ptr(ics:ice,iblk)           = veg_exudation%nitrogen(:)
  veg__mem%exudation%phosphorus%ptr(ics:ice,iblk)         = veg_exudation%phosphorus(:)
  veg__mem%exudation%carbon13%ptr(ics:ice,iblk)           = veg_exudation%carbon13(:)
  veg__mem%exudation%carbon14%ptr(ics:ice,iblk)           = veg_exudation%carbon14(:)
  veg__mem%exudation%nitrogen15%ptr(ics:ice,iblk)         = veg_exudation%nitrogen15(:)
  ! deallocate variables in helper pool structures
  deallocate(veg_exudation%carbon)
  deallocate(veg_exudation%nitrogen)
  deallocate(veg_exudation%phosphorus)
  deallocate(veg_exudation%carbon13)
  deallocate(veg_exudation%carbon14)
  deallocate(veg_exudation%nitrogen15)

#endif

#ifdef Passback_veg_establishment
  ! pass the values (back) to process memory
  veg__mem%establishment%carbon%ptr(ics:ice,iblk)             = veg_establishment%carbon(:)
  veg__mem%establishment%nitrogen%ptr(ics:ice,iblk)           = veg_establishment%nitrogen(:)
  veg__mem%establishment%phosphorus%ptr(ics:ice,iblk)         = veg_establishment%phosphorus(:)
  veg__mem%establishment%carbon13%ptr(ics:ice,iblk)           = veg_establishment%carbon13(:)
  veg__mem%establishment%carbon14%ptr(ics:ice,iblk)           = veg_establishment%carbon14(:)
  veg__mem%establishment%nitrogen15%ptr(ics:ice,iblk)         = veg_establishment%nitrogen15(:)
  ! deallocate variables in helper pool structures
  deallocate(veg_establishment%carbon)
  deallocate(veg_establishment%nitrogen)
  deallocate(veg_establishment%phosphorus)
  deallocate(veg_establishment%carbon13)
  deallocate(veg_establishment%carbon14)
  deallocate(veg_establishment%nitrogen15)

#endif

#ifdef Passback_veg_reserve_use
  ! pass the values (back) to process memory
  veg__mem%reserve_use%carbon%ptr(ics:ice,iblk)             = veg_reserve_use%carbon(:)
  veg__mem%reserve_use%nitrogen%ptr(ics:ice,iblk)           = veg_reserve_use%nitrogen(:)
  veg__mem%reserve_use%phosphorus%ptr(ics:ice,iblk)         = veg_reserve_use%phosphorus(:)
  veg__mem%reserve_use%carbon13%ptr(ics:ice,iblk)           = veg_reserve_use%carbon13(:)
  veg__mem%reserve_use%carbon14%ptr(ics:ice,iblk)           = veg_reserve_use%carbon14(:)
  veg__mem%reserve_use%nitrogen15%ptr(ics:ice,iblk)         = veg_reserve_use%nitrogen15(:)
  ! deallocate variables in helper pool structures
  deallocate(veg_reserve_use%carbon)
  deallocate(veg_reserve_use%nitrogen)
  deallocate(veg_reserve_use%phosphorus)
  deallocate(veg_reserve_use%carbon13)
  deallocate(veg_reserve_use%carbon14)
  deallocate(veg_reserve_use%nitrogen15)

#endif

#ifdef Passback_veg_total_biomass
  ! pass the values (back) to process memory
  veg__mem%total_biomass%carbon%ptr(ics:ice,iblk)             = veg_total_biomass%carbon(:)
  veg__mem%total_biomass%nitrogen%ptr(ics:ice,iblk)           = veg_total_biomass%nitrogen(:)
  veg__mem%total_biomass%phosphorus%ptr(ics:ice,iblk)         = veg_total_biomass%phosphorus(:)
  veg__mem%total_biomass%carbon13%ptr(ics:ice,iblk)           = veg_total_biomass%carbon13(:)
  veg__mem%total_biomass%carbon14%ptr(ics:ice,iblk)           = veg_total_biomass%carbon14(:)
  veg__mem%total_biomass%nitrogen15%ptr(ics:ice,iblk)         = veg_total_biomass%nitrogen15(:)
  ! deallocate variables in helper pool structures
  deallocate(veg_total_biomass%carbon)
  deallocate(veg_total_biomass%nitrogen)
  deallocate(veg_total_biomass%phosphorus)
  deallocate(veg_total_biomass%carbon13)
  deallocate(veg_total_biomass%carbon14)
  deallocate(veg_total_biomass%nitrogen15)

#endif


! -------------------------------------------------------
! deallocate helper variables for sb pools
! and pass values to pools of sb memory
! -------------------------------------------------------

#if defined(Passback_sb_pool)
  ! pass the values (back) to process memory
  sb__mem%pool%dom%carbon%ptr(ics:ice,1:nsoil_sb,iblk)                  = sb_pool%dom%carbon(:,:)
  sb__mem%pool%dom%nitrogen%ptr(ics:ice,1:nsoil_sb,iblk)                = sb_pool%dom%nitrogen(:,:)
  sb__mem%pool%dom%phosphorus%ptr(ics:ice,1:nsoil_sb,iblk)              = sb_pool%dom%phosphorus(:,:)
  sb__mem%pool%dom%carbon13%ptr(ics:ice,1:nsoil_sb,iblk)                = sb_pool%dom%carbon13(:,:)
  sb__mem%pool%dom%carbon14%ptr(ics:ice,1:nsoil_sb,iblk)                = sb_pool%dom%carbon14(:,:)
  sb__mem%pool%dom%nitrogen15%ptr(ics:ice,1:nsoil_sb,iblk)              = sb_pool%dom%nitrogen15(:,:)
  sb__mem%pool%dom_assoc%carbon%ptr(ics:ice,1:nsoil_sb,iblk)            = sb_pool%dom_assoc%carbon(:,:)
  sb__mem%pool%dom_assoc%nitrogen%ptr(ics:ice,1:nsoil_sb,iblk)          = sb_pool%dom_assoc%nitrogen(:,:)
  sb__mem%pool%dom_assoc%phosphorus%ptr(ics:ice,1:nsoil_sb,iblk)        = sb_pool%dom_assoc%phosphorus(:,:)
  sb__mem%pool%dom_assoc%carbon13%ptr(ics:ice,1:nsoil_sb,iblk)          = sb_pool%dom_assoc%carbon13(:,:)
  sb__mem%pool%dom_assoc%carbon14%ptr(ics:ice,1:nsoil_sb,iblk)          = sb_pool%dom_assoc%carbon14(:,:)
  sb__mem%pool%dom_assoc%nitrogen15%ptr(ics:ice,1:nsoil_sb,iblk)        = sb_pool%dom_assoc%nitrogen15(:,:)
  sb__mem%pool%soluable_litter%carbon%ptr(ics:ice,1:nsoil_sb,iblk)      = sb_pool%soluable_litter%carbon(:,:)
  sb__mem%pool%soluable_litter%nitrogen%ptr(ics:ice,1:nsoil_sb,iblk)    = sb_pool%soluable_litter%nitrogen(:,:)
  sb__mem%pool%soluable_litter%phosphorus%ptr(ics:ice,1:nsoil_sb,iblk)  = sb_pool%soluable_litter%phosphorus(:,:)
  sb__mem%pool%soluable_litter%carbon13%ptr(ics:ice,1:nsoil_sb,iblk)    = sb_pool%soluable_litter%carbon13(:,:)
  sb__mem%pool%soluable_litter%carbon14%ptr(ics:ice,1:nsoil_sb,iblk)    = sb_pool%soluable_litter%carbon14(:,:)
  sb__mem%pool%soluable_litter%nitrogen15%ptr(ics:ice,1:nsoil_sb,iblk)  = sb_pool%soluable_litter%nitrogen15(:,:)
  sb__mem%pool%polymeric_litter%carbon%ptr(ics:ice,1:nsoil_sb,iblk)     = sb_pool%polymeric_litter%carbon(:,:)
  sb__mem%pool%polymeric_litter%nitrogen%ptr(ics:ice,1:nsoil_sb,iblk)   = sb_pool%polymeric_litter%nitrogen(:,:)
  sb__mem%pool%polymeric_litter%phosphorus%ptr(ics:ice,1:nsoil_sb,iblk) = sb_pool%polymeric_litter%phosphorus(:,:)
  sb__mem%pool%polymeric_litter%carbon13%ptr(ics:ice,1:nsoil_sb,iblk)   = sb_pool%polymeric_litter%carbon13(:,:)
  sb__mem%pool%polymeric_litter%carbon14%ptr(ics:ice,1:nsoil_sb,iblk)   = sb_pool%polymeric_litter%carbon14(:,:)
  sb__mem%pool%polymeric_litter%nitrogen15%ptr(ics:ice,1:nsoil_sb,iblk) = sb_pool%polymeric_litter%nitrogen15(:,:)
  sb__mem%pool%woody_litter%carbon%ptr(ics:ice,1:nsoil_sb,iblk)         = sb_pool%woody_litter%carbon(:,:)
  sb__mem%pool%woody_litter%nitrogen%ptr(ics:ice,1:nsoil_sb,iblk)       = sb_pool%woody_litter%nitrogen(:,:)
  sb__mem%pool%woody_litter%phosphorus%ptr(ics:ice,1:nsoil_sb,iblk)     = sb_pool%woody_litter%phosphorus(:,:)
  sb__mem%pool%woody_litter%carbon13%ptr(ics:ice,1:nsoil_sb,iblk)       = sb_pool%woody_litter%carbon13(:,:)
  sb__mem%pool%woody_litter%carbon14%ptr(ics:ice,1:nsoil_sb,iblk)       = sb_pool%woody_litter%carbon14(:,:)
  sb__mem%pool%woody_litter%nitrogen15%ptr(ics:ice,1:nsoil_sb,iblk)     = sb_pool%woody_litter%nitrogen15(:,:)
  sb__mem%pool%fungi%carbon%ptr(ics:ice,1:nsoil_sb,iblk)                = sb_pool%fungi%carbon(:,:)
  sb__mem%pool%fungi%nitrogen%ptr(ics:ice,1:nsoil_sb,iblk)              = sb_pool%fungi%nitrogen(:,:)
  sb__mem%pool%fungi%phosphorus%ptr(ics:ice,1:nsoil_sb,iblk)            = sb_pool%fungi%phosphorus(:,:)
  sb__mem%pool%fungi%carbon13%ptr(ics:ice,1:nsoil_sb,iblk)              = sb_pool%fungi%carbon13(:,:)
  sb__mem%pool%fungi%carbon14%ptr(ics:ice,1:nsoil_sb,iblk)              = sb_pool%fungi%carbon14(:,:)
  sb__mem%pool%fungi%nitrogen15%ptr(ics:ice,1:nsoil_sb,iblk)            = sb_pool%fungi%nitrogen15(:,:)
  sb__mem%pool%mycorrhiza%carbon%ptr(ics:ice,1:nsoil_sb,iblk)           = sb_pool%mycorrhiza%carbon(:,:)
  sb__mem%pool%mycorrhiza%nitrogen%ptr(ics:ice,1:nsoil_sb,iblk)         = sb_pool%mycorrhiza%nitrogen(:,:)
  sb__mem%pool%mycorrhiza%phosphorus%ptr(ics:ice,1:nsoil_sb,iblk)       = sb_pool%mycorrhiza%phosphorus(:,:)
  sb__mem%pool%mycorrhiza%carbon13%ptr(ics:ice,1:nsoil_sb,iblk)         = sb_pool%mycorrhiza%carbon13(:,:)
  sb__mem%pool%mycorrhiza%carbon14%ptr(ics:ice,1:nsoil_sb,iblk)         = sb_pool%mycorrhiza%carbon14(:,:)
  sb__mem%pool%mycorrhiza%nitrogen15%ptr(ics:ice,1:nsoil_sb,iblk)       = sb_pool%mycorrhiza%nitrogen15(:,:)
  sb__mem%pool%microbial%carbon%ptr(ics:ice,1:nsoil_sb,iblk)            = sb_pool%microbial%carbon(:,:)
  sb__mem%pool%microbial%nitrogen%ptr(ics:ice,1:nsoil_sb,iblk)          = sb_pool%microbial%nitrogen(:,:)
  sb__mem%pool%microbial%phosphorus%ptr(ics:ice,1:nsoil_sb,iblk)        = sb_pool%microbial%phosphorus(:,:)
  sb__mem%pool%microbial%carbon13%ptr(ics:ice,1:nsoil_sb,iblk)          = sb_pool%microbial%carbon13(:,:)
  sb__mem%pool%microbial%carbon14%ptr(ics:ice,1:nsoil_sb,iblk)          = sb_pool%microbial%carbon14(:,:)
  sb__mem%pool%microbial%nitrogen15%ptr(ics:ice,1:nsoil_sb,iblk)        = sb_pool%microbial%nitrogen15(:,:)
  sb__mem%pool%residue%carbon%ptr(ics:ice,1:nsoil_sb,iblk)              = sb_pool%residue%carbon(:,:)
  sb__mem%pool%residue%nitrogen%ptr(ics:ice,1:nsoil_sb,iblk)            = sb_pool%residue%nitrogen(:,:)
  sb__mem%pool%residue%phosphorus%ptr(ics:ice,1:nsoil_sb,iblk)          = sb_pool%residue%phosphorus(:,:)
  sb__mem%pool%residue%carbon13%ptr(ics:ice,1:nsoil_sb,iblk)            = sb_pool%residue%carbon13(:,:)
  sb__mem%pool%residue%carbon14%ptr(ics:ice,1:nsoil_sb,iblk)            = sb_pool%residue%carbon14(:,:)
  sb__mem%pool%residue%nitrogen15%ptr(ics:ice,1:nsoil_sb,iblk)          = sb_pool%residue%nitrogen15(:,:)
  sb__mem%pool%residue_assoc%carbon%ptr(ics:ice,1:nsoil_sb,iblk)        = sb_pool%residue_assoc%carbon(:,:)
  sb__mem%pool%residue_assoc%nitrogen%ptr(ics:ice,1:nsoil_sb,iblk)      = sb_pool%residue_assoc%nitrogen(:,:)
  sb__mem%pool%residue_assoc%phosphorus%ptr(ics:ice,1:nsoil_sb,iblk)    = sb_pool%residue_assoc%phosphorus(:,:)
  sb__mem%pool%residue_assoc%carbon13%ptr(ics:ice,1:nsoil_sb,iblk)      = sb_pool%residue_assoc%carbon13(:,:)
  sb__mem%pool%residue_assoc%carbon14%ptr(ics:ice,1:nsoil_sb,iblk)      = sb_pool%residue_assoc%carbon14(:,:)
  sb__mem%pool%residue_assoc%nitrogen15%ptr(ics:ice,1:nsoil_sb,iblk)    = sb_pool%residue_assoc%nitrogen15(:,:)
  ! deallocate variables in helper pool structures
  deallocate(sb_pool%dom%carbon)
  deallocate(sb_pool%dom%nitrogen)
  deallocate(sb_pool%dom%phosphorus)
  deallocate(sb_pool%dom%carbon13)
  deallocate(sb_pool%dom%carbon14)
  deallocate(sb_pool%dom%nitrogen15)
  deallocate(sb_pool%dom_assoc%carbon)
  deallocate(sb_pool%dom_assoc%nitrogen)
  deallocate(sb_pool%dom_assoc%phosphorus)
  deallocate(sb_pool%dom_assoc%carbon13)
  deallocate(sb_pool%dom_assoc%carbon14)
  deallocate(sb_pool%dom_assoc%nitrogen15)
  deallocate(sb_pool%soluable_litter%carbon)
  deallocate(sb_pool%soluable_litter%nitrogen)
  deallocate(sb_pool%soluable_litter%phosphorus)
  deallocate(sb_pool%soluable_litter%carbon13)
  deallocate(sb_pool%soluable_litter%carbon14)
  deallocate(sb_pool%soluable_litter%nitrogen15)
  deallocate(sb_pool%polymeric_litter%carbon)
  deallocate(sb_pool%polymeric_litter%nitrogen)
  deallocate(sb_pool%polymeric_litter%phosphorus)
  deallocate(sb_pool%polymeric_litter%carbon13)
  deallocate(sb_pool%polymeric_litter%carbon14)
  deallocate(sb_pool%polymeric_litter%nitrogen15)
  deallocate(sb_pool%woody_litter%carbon)
  deallocate(sb_pool%woody_litter%nitrogen)
  deallocate(sb_pool%woody_litter%phosphorus)
  deallocate(sb_pool%woody_litter%carbon13)
  deallocate(sb_pool%woody_litter%carbon14)
  deallocate(sb_pool%woody_litter%nitrogen15)
  deallocate(sb_pool%fungi%carbon)
  deallocate(sb_pool%fungi%nitrogen)
  deallocate(sb_pool%fungi%phosphorus)
  deallocate(sb_pool%fungi%carbon13)
  deallocate(sb_pool%fungi%carbon14)
  deallocate(sb_pool%fungi%nitrogen15)
  deallocate(sb_pool%mycorrhiza%carbon)
  deallocate(sb_pool%mycorrhiza%nitrogen)
  deallocate(sb_pool%mycorrhiza%phosphorus)
  deallocate(sb_pool%mycorrhiza%carbon13)
  deallocate(sb_pool%mycorrhiza%carbon14)
  deallocate(sb_pool%mycorrhiza%nitrogen15)
  deallocate(sb_pool%microbial%carbon)
  deallocate(sb_pool%microbial%nitrogen)
  deallocate(sb_pool%microbial%phosphorus)
  deallocate(sb_pool%microbial%carbon13)
  deallocate(sb_pool%microbial%carbon14)
  deallocate(sb_pool%microbial%nitrogen15)
  deallocate(sb_pool%residue%carbon)
  deallocate(sb_pool%residue%nitrogen)
  deallocate(sb_pool%residue%phosphorus)
  deallocate(sb_pool%residue%carbon13)
  deallocate(sb_pool%residue%carbon14)
  deallocate(sb_pool%residue%nitrogen15)
  deallocate(sb_pool%residue_assoc%carbon)
  deallocate(sb_pool%residue_assoc%nitrogen)
  deallocate(sb_pool%residue_assoc%phosphorus)
  deallocate(sb_pool%residue_assoc%carbon13)
  deallocate(sb_pool%residue_assoc%carbon14)
  deallocate(sb_pool%residue_assoc%nitrogen15)

#endif

#if defined(Passback_sb_formation)
  ! pass the values (back) to process memory
  sb__mem%formation%dom%carbon%ptr(ics:ice,1:nsoil_sb,iblk)                  = sb_formation%dom%carbon(:,:)
  sb__mem%formation%dom%nitrogen%ptr(ics:ice,1:nsoil_sb,iblk)                = sb_formation%dom%nitrogen(:,:)
  sb__mem%formation%dom%phosphorus%ptr(ics:ice,1:nsoil_sb,iblk)              = sb_formation%dom%phosphorus(:,:)
  sb__mem%formation%dom%carbon13%ptr(ics:ice,1:nsoil_sb,iblk)                = sb_formation%dom%carbon13(:,:)
  sb__mem%formation%dom%carbon14%ptr(ics:ice,1:nsoil_sb,iblk)                = sb_formation%dom%carbon14(:,:)
  sb__mem%formation%dom%nitrogen15%ptr(ics:ice,1:nsoil_sb,iblk)              = sb_formation%dom%nitrogen15(:,:)
  sb__mem%formation%dom_assoc%carbon%ptr(ics:ice,1:nsoil_sb,iblk)            = sb_formation%dom_assoc%carbon(:,:)
  sb__mem%formation%dom_assoc%nitrogen%ptr(ics:ice,1:nsoil_sb,iblk)          = sb_formation%dom_assoc%nitrogen(:,:)
  sb__mem%formation%dom_assoc%phosphorus%ptr(ics:ice,1:nsoil_sb,iblk)        = sb_formation%dom_assoc%phosphorus(:,:)
  sb__mem%formation%dom_assoc%carbon13%ptr(ics:ice,1:nsoil_sb,iblk)          = sb_formation%dom_assoc%carbon13(:,:)
  sb__mem%formation%dom_assoc%carbon14%ptr(ics:ice,1:nsoil_sb,iblk)          = sb_formation%dom_assoc%carbon14(:,:)
  sb__mem%formation%dom_assoc%nitrogen15%ptr(ics:ice,1:nsoil_sb,iblk)        = sb_formation%dom_assoc%nitrogen15(:,:)
  sb__mem%formation%soluable_litter%carbon%ptr(ics:ice,1:nsoil_sb,iblk)      = sb_formation%soluable_litter%carbon(:,:)
  sb__mem%formation%soluable_litter%nitrogen%ptr(ics:ice,1:nsoil_sb,iblk)    = sb_formation%soluable_litter%nitrogen(:,:)
  sb__mem%formation%soluable_litter%phosphorus%ptr(ics:ice,1:nsoil_sb,iblk)  = sb_formation%soluable_litter%phosphorus(:,:)
  sb__mem%formation%soluable_litter%carbon13%ptr(ics:ice,1:nsoil_sb,iblk)    = sb_formation%soluable_litter%carbon13(:,:)
  sb__mem%formation%soluable_litter%carbon14%ptr(ics:ice,1:nsoil_sb,iblk)    = sb_formation%soluable_litter%carbon14(:,:)
  sb__mem%formation%soluable_litter%nitrogen15%ptr(ics:ice,1:nsoil_sb,iblk)  = sb_formation%soluable_litter%nitrogen15(:,:)
  sb__mem%formation%polymeric_litter%carbon%ptr(ics:ice,1:nsoil_sb,iblk)     = sb_formation%polymeric_litter%carbon(:,:)
  sb__mem%formation%polymeric_litter%nitrogen%ptr(ics:ice,1:nsoil_sb,iblk)   = sb_formation%polymeric_litter%nitrogen(:,:)
  sb__mem%formation%polymeric_litter%phosphorus%ptr(ics:ice,1:nsoil_sb,iblk) = sb_formation%polymeric_litter%phosphorus(:,:)
  sb__mem%formation%polymeric_litter%carbon13%ptr(ics:ice,1:nsoil_sb,iblk)   = sb_formation%polymeric_litter%carbon13(:,:)
  sb__mem%formation%polymeric_litter%carbon14%ptr(ics:ice,1:nsoil_sb,iblk)   = sb_formation%polymeric_litter%carbon14(:,:)
  sb__mem%formation%polymeric_litter%nitrogen15%ptr(ics:ice,1:nsoil_sb,iblk) = sb_formation%polymeric_litter%nitrogen15(:,:)
  sb__mem%formation%woody_litter%carbon%ptr(ics:ice,1:nsoil_sb,iblk)         = sb_formation%woody_litter%carbon(:,:)
  sb__mem%formation%woody_litter%nitrogen%ptr(ics:ice,1:nsoil_sb,iblk)       = sb_formation%woody_litter%nitrogen(:,:)
  sb__mem%formation%woody_litter%phosphorus%ptr(ics:ice,1:nsoil_sb,iblk)     = sb_formation%woody_litter%phosphorus(:,:)
  sb__mem%formation%woody_litter%carbon13%ptr(ics:ice,1:nsoil_sb,iblk)       = sb_formation%woody_litter%carbon13(:,:)
  sb__mem%formation%woody_litter%carbon14%ptr(ics:ice,1:nsoil_sb,iblk)       = sb_formation%woody_litter%carbon14(:,:)
  sb__mem%formation%woody_litter%nitrogen15%ptr(ics:ice,1:nsoil_sb,iblk)     = sb_formation%woody_litter%nitrogen15(:,:)
  sb__mem%formation%fungi%carbon%ptr(ics:ice,1:nsoil_sb,iblk)                = sb_formation%fungi%carbon(:,:)
  sb__mem%formation%fungi%nitrogen%ptr(ics:ice,1:nsoil_sb,iblk)              = sb_formation%fungi%nitrogen(:,:)
  sb__mem%formation%fungi%phosphorus%ptr(ics:ice,1:nsoil_sb,iblk)            = sb_formation%fungi%phosphorus(:,:)
  sb__mem%formation%fungi%carbon13%ptr(ics:ice,1:nsoil_sb,iblk)              = sb_formation%fungi%carbon13(:,:)
  sb__mem%formation%fungi%carbon14%ptr(ics:ice,1:nsoil_sb,iblk)              = sb_formation%fungi%carbon14(:,:)
  sb__mem%formation%fungi%nitrogen15%ptr(ics:ice,1:nsoil_sb,iblk)            = sb_formation%fungi%nitrogen15(:,:)
  sb__mem%formation%mycorrhiza%carbon%ptr(ics:ice,1:nsoil_sb,iblk)           = sb_formation%mycorrhiza%carbon(:,:)
  sb__mem%formation%mycorrhiza%nitrogen%ptr(ics:ice,1:nsoil_sb,iblk)         = sb_formation%mycorrhiza%nitrogen(:,:)
  sb__mem%formation%mycorrhiza%phosphorus%ptr(ics:ice,1:nsoil_sb,iblk)       = sb_formation%mycorrhiza%phosphorus(:,:)
  sb__mem%formation%mycorrhiza%carbon13%ptr(ics:ice,1:nsoil_sb,iblk)         = sb_formation%mycorrhiza%carbon13(:,:)
  sb__mem%formation%mycorrhiza%carbon14%ptr(ics:ice,1:nsoil_sb,iblk)         = sb_formation%mycorrhiza%carbon14(:,:)
  sb__mem%formation%mycorrhiza%nitrogen15%ptr(ics:ice,1:nsoil_sb,iblk)       = sb_formation%mycorrhiza%nitrogen15(:,:)
  sb__mem%formation%microbial%carbon%ptr(ics:ice,1:nsoil_sb,iblk)            = sb_formation%microbial%carbon(:,:)
  sb__mem%formation%microbial%nitrogen%ptr(ics:ice,1:nsoil_sb,iblk)          = sb_formation%microbial%nitrogen(:,:)
  sb__mem%formation%microbial%phosphorus%ptr(ics:ice,1:nsoil_sb,iblk)        = sb_formation%microbial%phosphorus(:,:)
  sb__mem%formation%microbial%carbon13%ptr(ics:ice,1:nsoil_sb,iblk)          = sb_formation%microbial%carbon13(:,:)
  sb__mem%formation%microbial%carbon14%ptr(ics:ice,1:nsoil_sb,iblk)          = sb_formation%microbial%carbon14(:,:)
  sb__mem%formation%microbial%nitrogen15%ptr(ics:ice,1:nsoil_sb,iblk)        = sb_formation%microbial%nitrogen15(:,:)
  sb__mem%formation%residue%carbon%ptr(ics:ice,1:nsoil_sb,iblk)              = sb_formation%residue%carbon(:,:)
  sb__mem%formation%residue%nitrogen%ptr(ics:ice,1:nsoil_sb,iblk)            = sb_formation%residue%nitrogen(:,:)
  sb__mem%formation%residue%phosphorus%ptr(ics:ice,1:nsoil_sb,iblk)          = sb_formation%residue%phosphorus(:,:)
  sb__mem%formation%residue%carbon13%ptr(ics:ice,1:nsoil_sb,iblk)            = sb_formation%residue%carbon13(:,:)
  sb__mem%formation%residue%carbon14%ptr(ics:ice,1:nsoil_sb,iblk)            = sb_formation%residue%carbon14(:,:)
  sb__mem%formation%residue%nitrogen15%ptr(ics:ice,1:nsoil_sb,iblk)          = sb_formation%residue%nitrogen15(:,:)
  sb__mem%formation%residue_assoc%carbon%ptr(ics:ice,1:nsoil_sb,iblk)        = sb_formation%residue_assoc%carbon(:,:)
  sb__mem%formation%residue_assoc%nitrogen%ptr(ics:ice,1:nsoil_sb,iblk)      = sb_formation%residue_assoc%nitrogen(:,:)
  sb__mem%formation%residue_assoc%phosphorus%ptr(ics:ice,1:nsoil_sb,iblk)    = sb_formation%residue_assoc%phosphorus(:,:)
  sb__mem%formation%residue_assoc%carbon13%ptr(ics:ice,1:nsoil_sb,iblk)      = sb_formation%residue_assoc%carbon13(:,:)
  sb__mem%formation%residue_assoc%carbon14%ptr(ics:ice,1:nsoil_sb,iblk)      = sb_formation%residue_assoc%carbon14(:,:)
  sb__mem%formation%residue_assoc%nitrogen15%ptr(ics:ice,1:nsoil_sb,iblk)    = sb_formation%residue_assoc%nitrogen15(:,:)
  ! deallocate variables in helper pool structures
  deallocate(sb_formation%dom%carbon)
  deallocate(sb_formation%dom%nitrogen)
  deallocate(sb_formation%dom%phosphorus)
  deallocate(sb_formation%dom%carbon13)
  deallocate(sb_formation%dom%carbon14)
  deallocate(sb_formation%dom%nitrogen15)
  deallocate(sb_formation%dom_assoc%carbon)
  deallocate(sb_formation%dom_assoc%nitrogen)
  deallocate(sb_formation%dom_assoc%phosphorus)
  deallocate(sb_formation%dom_assoc%carbon13)
  deallocate(sb_formation%dom_assoc%carbon14)
  deallocate(sb_formation%dom_assoc%nitrogen15)
  deallocate(sb_formation%soluable_litter%carbon)
  deallocate(sb_formation%soluable_litter%nitrogen)
  deallocate(sb_formation%soluable_litter%phosphorus)
  deallocate(sb_formation%soluable_litter%carbon13)
  deallocate(sb_formation%soluable_litter%carbon14)
  deallocate(sb_formation%soluable_litter%nitrogen15)
  deallocate(sb_formation%polymeric_litter%carbon)
  deallocate(sb_formation%polymeric_litter%nitrogen)
  deallocate(sb_formation%polymeric_litter%phosphorus)
  deallocate(sb_formation%polymeric_litter%carbon13)
  deallocate(sb_formation%polymeric_litter%carbon14)
  deallocate(sb_formation%polymeric_litter%nitrogen15)
  deallocate(sb_formation%woody_litter%carbon)
  deallocate(sb_formation%woody_litter%nitrogen)
  deallocate(sb_formation%woody_litter%phosphorus)
  deallocate(sb_formation%woody_litter%carbon13)
  deallocate(sb_formation%woody_litter%carbon14)
  deallocate(sb_formation%woody_litter%nitrogen15)
  deallocate(sb_formation%fungi%carbon)
  deallocate(sb_formation%fungi%nitrogen)
  deallocate(sb_formation%fungi%phosphorus)
  deallocate(sb_formation%fungi%carbon13)
  deallocate(sb_formation%fungi%carbon14)
  deallocate(sb_formation%fungi%nitrogen15)
  deallocate(sb_formation%mycorrhiza%carbon)
  deallocate(sb_formation%mycorrhiza%nitrogen)
  deallocate(sb_formation%mycorrhiza%phosphorus)
  deallocate(sb_formation%mycorrhiza%carbon13)
  deallocate(sb_formation%mycorrhiza%carbon14)
  deallocate(sb_formation%mycorrhiza%nitrogen15)
  deallocate(sb_formation%microbial%carbon)
  deallocate(sb_formation%microbial%nitrogen)
  deallocate(sb_formation%microbial%phosphorus)
  deallocate(sb_formation%microbial%carbon13)
  deallocate(sb_formation%microbial%carbon14)
  deallocate(sb_formation%microbial%nitrogen15)
  deallocate(sb_formation%residue%carbon)
  deallocate(sb_formation%residue%nitrogen)
  deallocate(sb_formation%residue%phosphorus)
  deallocate(sb_formation%residue%carbon13)
  deallocate(sb_formation%residue%carbon14)
  deallocate(sb_formation%residue%nitrogen15)
  deallocate(sb_formation%residue_assoc%carbon)
  deallocate(sb_formation%residue_assoc%nitrogen)
  deallocate(sb_formation%residue_assoc%phosphorus)
  deallocate(sb_formation%residue_assoc%carbon13)
  deallocate(sb_formation%residue_assoc%carbon14)
  deallocate(sb_formation%residue_assoc%nitrogen15)

#endif

#if defined(Passback_sb_loss)
  ! pass the values (back) to process memory
  sb__mem%loss%dom%carbon%ptr(ics:ice,1:nsoil_sb,iblk)                  = sb_loss%dom%carbon(:,:)
  sb__mem%loss%dom%nitrogen%ptr(ics:ice,1:nsoil_sb,iblk)                = sb_loss%dom%nitrogen(:,:)
  sb__mem%loss%dom%phosphorus%ptr(ics:ice,1:nsoil_sb,iblk)              = sb_loss%dom%phosphorus(:,:)
  sb__mem%loss%dom%carbon13%ptr(ics:ice,1:nsoil_sb,iblk)                = sb_loss%dom%carbon13(:,:)
  sb__mem%loss%dom%carbon14%ptr(ics:ice,1:nsoil_sb,iblk)                = sb_loss%dom%carbon14(:,:)
  sb__mem%loss%dom%nitrogen15%ptr(ics:ice,1:nsoil_sb,iblk)              = sb_loss%dom%nitrogen15(:,:)
  sb__mem%loss%dom_assoc%carbon%ptr(ics:ice,1:nsoil_sb,iblk)            = sb_loss%dom_assoc%carbon(:,:)
  sb__mem%loss%dom_assoc%nitrogen%ptr(ics:ice,1:nsoil_sb,iblk)          = sb_loss%dom_assoc%nitrogen(:,:)
  sb__mem%loss%dom_assoc%phosphorus%ptr(ics:ice,1:nsoil_sb,iblk)        = sb_loss%dom_assoc%phosphorus(:,:)
  sb__mem%loss%dom_assoc%carbon13%ptr(ics:ice,1:nsoil_sb,iblk)          = sb_loss%dom_assoc%carbon13(:,:)
  sb__mem%loss%dom_assoc%carbon14%ptr(ics:ice,1:nsoil_sb,iblk)          = sb_loss%dom_assoc%carbon14(:,:)
  sb__mem%loss%dom_assoc%nitrogen15%ptr(ics:ice,1:nsoil_sb,iblk)        = sb_loss%dom_assoc%nitrogen15(:,:)
  sb__mem%loss%soluable_litter%carbon%ptr(ics:ice,1:nsoil_sb,iblk)      = sb_loss%soluable_litter%carbon(:,:)
  sb__mem%loss%soluable_litter%nitrogen%ptr(ics:ice,1:nsoil_sb,iblk)    = sb_loss%soluable_litter%nitrogen(:,:)
  sb__mem%loss%soluable_litter%phosphorus%ptr(ics:ice,1:nsoil_sb,iblk)  = sb_loss%soluable_litter%phosphorus(:,:)
  sb__mem%loss%soluable_litter%carbon13%ptr(ics:ice,1:nsoil_sb,iblk)    = sb_loss%soluable_litter%carbon13(:,:)
  sb__mem%loss%soluable_litter%carbon14%ptr(ics:ice,1:nsoil_sb,iblk)    = sb_loss%soluable_litter%carbon14(:,:)
  sb__mem%loss%soluable_litter%nitrogen15%ptr(ics:ice,1:nsoil_sb,iblk)  = sb_loss%soluable_litter%nitrogen15(:,:)
  sb__mem%loss%polymeric_litter%carbon%ptr(ics:ice,1:nsoil_sb,iblk)     = sb_loss%polymeric_litter%carbon(:,:)
  sb__mem%loss%polymeric_litter%nitrogen%ptr(ics:ice,1:nsoil_sb,iblk)   = sb_loss%polymeric_litter%nitrogen(:,:)
  sb__mem%loss%polymeric_litter%phosphorus%ptr(ics:ice,1:nsoil_sb,iblk) = sb_loss%polymeric_litter%phosphorus(:,:)
  sb__mem%loss%polymeric_litter%carbon13%ptr(ics:ice,1:nsoil_sb,iblk)   = sb_loss%polymeric_litter%carbon13(:,:)
  sb__mem%loss%polymeric_litter%carbon14%ptr(ics:ice,1:nsoil_sb,iblk)   = sb_loss%polymeric_litter%carbon14(:,:)
  sb__mem%loss%polymeric_litter%nitrogen15%ptr(ics:ice,1:nsoil_sb,iblk) = sb_loss%polymeric_litter%nitrogen15(:,:)
  sb__mem%loss%woody_litter%carbon%ptr(ics:ice,1:nsoil_sb,iblk)         = sb_loss%woody_litter%carbon(:,:)
  sb__mem%loss%woody_litter%nitrogen%ptr(ics:ice,1:nsoil_sb,iblk)       = sb_loss%woody_litter%nitrogen(:,:)
  sb__mem%loss%woody_litter%phosphorus%ptr(ics:ice,1:nsoil_sb,iblk)     = sb_loss%woody_litter%phosphorus(:,:)
  sb__mem%loss%woody_litter%carbon13%ptr(ics:ice,1:nsoil_sb,iblk)       = sb_loss%woody_litter%carbon13(:,:)
  sb__mem%loss%woody_litter%carbon14%ptr(ics:ice,1:nsoil_sb,iblk)       = sb_loss%woody_litter%carbon14(:,:)
  sb__mem%loss%woody_litter%nitrogen15%ptr(ics:ice,1:nsoil_sb,iblk)     = sb_loss%woody_litter%nitrogen15(:,:)
  sb__mem%loss%fungi%carbon%ptr(ics:ice,1:nsoil_sb,iblk)                = sb_loss%fungi%carbon(:,:)
  sb__mem%loss%fungi%nitrogen%ptr(ics:ice,1:nsoil_sb,iblk)              = sb_loss%fungi%nitrogen(:,:)
  sb__mem%loss%fungi%phosphorus%ptr(ics:ice,1:nsoil_sb,iblk)            = sb_loss%fungi%phosphorus(:,:)
  sb__mem%loss%fungi%carbon13%ptr(ics:ice,1:nsoil_sb,iblk)              = sb_loss%fungi%carbon13(:,:)
  sb__mem%loss%fungi%carbon14%ptr(ics:ice,1:nsoil_sb,iblk)              = sb_loss%fungi%carbon14(:,:)
  sb__mem%loss%fungi%nitrogen15%ptr(ics:ice,1:nsoil_sb,iblk)            = sb_loss%fungi%nitrogen15(:,:)
  sb__mem%loss%mycorrhiza%carbon%ptr(ics:ice,1:nsoil_sb,iblk)           = sb_loss%mycorrhiza%carbon(:,:)
  sb__mem%loss%mycorrhiza%nitrogen%ptr(ics:ice,1:nsoil_sb,iblk)         = sb_loss%mycorrhiza%nitrogen(:,:)
  sb__mem%loss%mycorrhiza%phosphorus%ptr(ics:ice,1:nsoil_sb,iblk)       = sb_loss%mycorrhiza%phosphorus(:,:)
  sb__mem%loss%mycorrhiza%carbon13%ptr(ics:ice,1:nsoil_sb,iblk)         = sb_loss%mycorrhiza%carbon13(:,:)
  sb__mem%loss%mycorrhiza%carbon14%ptr(ics:ice,1:nsoil_sb,iblk)         = sb_loss%mycorrhiza%carbon14(:,:)
  sb__mem%loss%mycorrhiza%nitrogen15%ptr(ics:ice,1:nsoil_sb,iblk)       = sb_loss%mycorrhiza%nitrogen15(:,:)
  sb__mem%loss%microbial%carbon%ptr(ics:ice,1:nsoil_sb,iblk)            = sb_loss%microbial%carbon(:,:)
  sb__mem%loss%microbial%nitrogen%ptr(ics:ice,1:nsoil_sb,iblk)          = sb_loss%microbial%nitrogen(:,:)
  sb__mem%loss%microbial%phosphorus%ptr(ics:ice,1:nsoil_sb,iblk)        = sb_loss%microbial%phosphorus(:,:)
  sb__mem%loss%microbial%carbon13%ptr(ics:ice,1:nsoil_sb,iblk)          = sb_loss%microbial%carbon13(:,:)
  sb__mem%loss%microbial%carbon14%ptr(ics:ice,1:nsoil_sb,iblk)          = sb_loss%microbial%carbon14(:,:)
  sb__mem%loss%microbial%nitrogen15%ptr(ics:ice,1:nsoil_sb,iblk)        = sb_loss%microbial%nitrogen15(:,:)
  sb__mem%loss%residue%carbon%ptr(ics:ice,1:nsoil_sb,iblk)              = sb_loss%residue%carbon(:,:)
  sb__mem%loss%residue%nitrogen%ptr(ics:ice,1:nsoil_sb,iblk)            = sb_loss%residue%nitrogen(:,:)
  sb__mem%loss%residue%phosphorus%ptr(ics:ice,1:nsoil_sb,iblk)          = sb_loss%residue%phosphorus(:,:)
  sb__mem%loss%residue%carbon13%ptr(ics:ice,1:nsoil_sb,iblk)            = sb_loss%residue%carbon13(:,:)
  sb__mem%loss%residue%carbon14%ptr(ics:ice,1:nsoil_sb,iblk)            = sb_loss%residue%carbon14(:,:)
  sb__mem%loss%residue%nitrogen15%ptr(ics:ice,1:nsoil_sb,iblk)          = sb_loss%residue%nitrogen15(:,:)
  sb__mem%loss%residue_assoc%carbon%ptr(ics:ice,1:nsoil_sb,iblk)        = sb_loss%residue_assoc%carbon(:,:)
  sb__mem%loss%residue_assoc%nitrogen%ptr(ics:ice,1:nsoil_sb,iblk)      = sb_loss%residue_assoc%nitrogen(:,:)
  sb__mem%loss%residue_assoc%phosphorus%ptr(ics:ice,1:nsoil_sb,iblk)    = sb_loss%residue_assoc%phosphorus(:,:)
  sb__mem%loss%residue_assoc%carbon13%ptr(ics:ice,1:nsoil_sb,iblk)      = sb_loss%residue_assoc%carbon13(:,:)
  sb__mem%loss%residue_assoc%carbon14%ptr(ics:ice,1:nsoil_sb,iblk)      = sb_loss%residue_assoc%carbon14(:,:)
  sb__mem%loss%residue_assoc%nitrogen15%ptr(ics:ice,1:nsoil_sb,iblk)    = sb_loss%residue_assoc%nitrogen15(:,:)
  ! deallocate variables in helper pool structures
  deallocate(sb_loss%dom%carbon)
  deallocate(sb_loss%dom%nitrogen)
  deallocate(sb_loss%dom%phosphorus)
  deallocate(sb_loss%dom%carbon13)
  deallocate(sb_loss%dom%carbon14)
  deallocate(sb_loss%dom%nitrogen15)
  deallocate(sb_loss%dom_assoc%carbon)
  deallocate(sb_loss%dom_assoc%nitrogen)
  deallocate(sb_loss%dom_assoc%phosphorus)
  deallocate(sb_loss%dom_assoc%carbon13)
  deallocate(sb_loss%dom_assoc%carbon14)
  deallocate(sb_loss%dom_assoc%nitrogen15)
  deallocate(sb_loss%soluable_litter%carbon)
  deallocate(sb_loss%soluable_litter%nitrogen)
  deallocate(sb_loss%soluable_litter%phosphorus)
  deallocate(sb_loss%soluable_litter%carbon13)
  deallocate(sb_loss%soluable_litter%carbon14)
  deallocate(sb_loss%soluable_litter%nitrogen15)
  deallocate(sb_loss%polymeric_litter%carbon)
  deallocate(sb_loss%polymeric_litter%nitrogen)
  deallocate(sb_loss%polymeric_litter%phosphorus)
  deallocate(sb_loss%polymeric_litter%carbon13)
  deallocate(sb_loss%polymeric_litter%carbon14)
  deallocate(sb_loss%polymeric_litter%nitrogen15)
  deallocate(sb_loss%woody_litter%carbon)
  deallocate(sb_loss%woody_litter%nitrogen)
  deallocate(sb_loss%woody_litter%phosphorus)
  deallocate(sb_loss%woody_litter%carbon13)
  deallocate(sb_loss%woody_litter%carbon14)
  deallocate(sb_loss%woody_litter%nitrogen15)
  deallocate(sb_loss%fungi%carbon)
  deallocate(sb_loss%fungi%nitrogen)
  deallocate(sb_loss%fungi%phosphorus)
  deallocate(sb_loss%fungi%carbon13)
  deallocate(sb_loss%fungi%carbon14)
  deallocate(sb_loss%fungi%nitrogen15)
  deallocate(sb_loss%mycorrhiza%carbon)
  deallocate(sb_loss%mycorrhiza%nitrogen)
  deallocate(sb_loss%mycorrhiza%phosphorus)
  deallocate(sb_loss%mycorrhiza%carbon13)
  deallocate(sb_loss%mycorrhiza%carbon14)
  deallocate(sb_loss%mycorrhiza%nitrogen15)
  deallocate(sb_loss%microbial%carbon)
  deallocate(sb_loss%microbial%nitrogen)
  deallocate(sb_loss%microbial%phosphorus)
  deallocate(sb_loss%microbial%carbon13)
  deallocate(sb_loss%microbial%carbon14)
  deallocate(sb_loss%microbial%nitrogen15)
  deallocate(sb_loss%residue%carbon)
  deallocate(sb_loss%residue%nitrogen)
  deallocate(sb_loss%residue%phosphorus)
  deallocate(sb_loss%residue%carbon13)
  deallocate(sb_loss%residue%carbon14)
  deallocate(sb_loss%residue%nitrogen15)
  deallocate(sb_loss%residue_assoc%carbon)
  deallocate(sb_loss%residue_assoc%nitrogen)
  deallocate(sb_loss%residue_assoc%phosphorus)
  deallocate(sb_loss%residue_assoc%carbon13)
  deallocate(sb_loss%residue_assoc%carbon14)
  deallocate(sb_loss%residue_assoc%nitrogen15)

#endif

#if defined(Passback_sb_transport)
  ! pass the values (back) to process memory
  sb__mem%transport%dom%carbon%ptr(ics:ice,1:nsoil_sb,iblk)                  = sb_transport%dom%carbon(:,:)
  sb__mem%transport%dom%nitrogen%ptr(ics:ice,1:nsoil_sb,iblk)                = sb_transport%dom%nitrogen(:,:)
  sb__mem%transport%dom%phosphorus%ptr(ics:ice,1:nsoil_sb,iblk)              = sb_transport%dom%phosphorus(:,:)
  sb__mem%transport%dom%carbon13%ptr(ics:ice,1:nsoil_sb,iblk)                = sb_transport%dom%carbon13(:,:)
  sb__mem%transport%dom%carbon14%ptr(ics:ice,1:nsoil_sb,iblk)                = sb_transport%dom%carbon14(:,:)
  sb__mem%transport%dom%nitrogen15%ptr(ics:ice,1:nsoil_sb,iblk)              = sb_transport%dom%nitrogen15(:,:)
  sb__mem%transport%dom_assoc%carbon%ptr(ics:ice,1:nsoil_sb,iblk)            = sb_transport%dom_assoc%carbon(:,:)
  sb__mem%transport%dom_assoc%nitrogen%ptr(ics:ice,1:nsoil_sb,iblk)          = sb_transport%dom_assoc%nitrogen(:,:)
  sb__mem%transport%dom_assoc%phosphorus%ptr(ics:ice,1:nsoil_sb,iblk)        = sb_transport%dom_assoc%phosphorus(:,:)
  sb__mem%transport%dom_assoc%carbon13%ptr(ics:ice,1:nsoil_sb,iblk)          = sb_transport%dom_assoc%carbon13(:,:)
  sb__mem%transport%dom_assoc%carbon14%ptr(ics:ice,1:nsoil_sb,iblk)          = sb_transport%dom_assoc%carbon14(:,:)
  sb__mem%transport%dom_assoc%nitrogen15%ptr(ics:ice,1:nsoil_sb,iblk)        = sb_transport%dom_assoc%nitrogen15(:,:)
  sb__mem%transport%soluable_litter%carbon%ptr(ics:ice,1:nsoil_sb,iblk)      = sb_transport%soluable_litter%carbon(:,:)
  sb__mem%transport%soluable_litter%nitrogen%ptr(ics:ice,1:nsoil_sb,iblk)    = sb_transport%soluable_litter%nitrogen(:,:)
  sb__mem%transport%soluable_litter%phosphorus%ptr(ics:ice,1:nsoil_sb,iblk)  = sb_transport%soluable_litter%phosphorus(:,:)
  sb__mem%transport%soluable_litter%carbon13%ptr(ics:ice,1:nsoil_sb,iblk)    = sb_transport%soluable_litter%carbon13(:,:)
  sb__mem%transport%soluable_litter%carbon14%ptr(ics:ice,1:nsoil_sb,iblk)    = sb_transport%soluable_litter%carbon14(:,:)
  sb__mem%transport%soluable_litter%nitrogen15%ptr(ics:ice,1:nsoil_sb,iblk)  = sb_transport%soluable_litter%nitrogen15(:,:)
  sb__mem%transport%polymeric_litter%carbon%ptr(ics:ice,1:nsoil_sb,iblk)     = sb_transport%polymeric_litter%carbon(:,:)
  sb__mem%transport%polymeric_litter%nitrogen%ptr(ics:ice,1:nsoil_sb,iblk)   = sb_transport%polymeric_litter%nitrogen(:,:)
  sb__mem%transport%polymeric_litter%phosphorus%ptr(ics:ice,1:nsoil_sb,iblk) = sb_transport%polymeric_litter%phosphorus(:,:)
  sb__mem%transport%polymeric_litter%carbon13%ptr(ics:ice,1:nsoil_sb,iblk)   = sb_transport%polymeric_litter%carbon13(:,:)
  sb__mem%transport%polymeric_litter%carbon14%ptr(ics:ice,1:nsoil_sb,iblk)   = sb_transport%polymeric_litter%carbon14(:,:)
  sb__mem%transport%polymeric_litter%nitrogen15%ptr(ics:ice,1:nsoil_sb,iblk) = sb_transport%polymeric_litter%nitrogen15(:,:)
  sb__mem%transport%woody_litter%carbon%ptr(ics:ice,1:nsoil_sb,iblk)         = sb_transport%woody_litter%carbon(:,:)
  sb__mem%transport%woody_litter%nitrogen%ptr(ics:ice,1:nsoil_sb,iblk)       = sb_transport%woody_litter%nitrogen(:,:)
  sb__mem%transport%woody_litter%phosphorus%ptr(ics:ice,1:nsoil_sb,iblk)     = sb_transport%woody_litter%phosphorus(:,:)
  sb__mem%transport%woody_litter%carbon13%ptr(ics:ice,1:nsoil_sb,iblk)       = sb_transport%woody_litter%carbon13(:,:)
  sb__mem%transport%woody_litter%carbon14%ptr(ics:ice,1:nsoil_sb,iblk)       = sb_transport%woody_litter%carbon14(:,:)
  sb__mem%transport%woody_litter%nitrogen15%ptr(ics:ice,1:nsoil_sb,iblk)     = sb_transport%woody_litter%nitrogen15(:,:)
  sb__mem%transport%fungi%carbon%ptr(ics:ice,1:nsoil_sb,iblk)                = sb_transport%fungi%carbon(:,:)
  sb__mem%transport%fungi%nitrogen%ptr(ics:ice,1:nsoil_sb,iblk)              = sb_transport%fungi%nitrogen(:,:)
  sb__mem%transport%fungi%phosphorus%ptr(ics:ice,1:nsoil_sb,iblk)            = sb_transport%fungi%phosphorus(:,:)
  sb__mem%transport%fungi%carbon13%ptr(ics:ice,1:nsoil_sb,iblk)              = sb_transport%fungi%carbon13(:,:)
  sb__mem%transport%fungi%carbon14%ptr(ics:ice,1:nsoil_sb,iblk)              = sb_transport%fungi%carbon14(:,:)
  sb__mem%transport%fungi%nitrogen15%ptr(ics:ice,1:nsoil_sb,iblk)            = sb_transport%fungi%nitrogen15(:,:)
  sb__mem%transport%mycorrhiza%carbon%ptr(ics:ice,1:nsoil_sb,iblk)           = sb_transport%mycorrhiza%carbon(:,:)
  sb__mem%transport%mycorrhiza%nitrogen%ptr(ics:ice,1:nsoil_sb,iblk)         = sb_transport%mycorrhiza%nitrogen(:,:)
  sb__mem%transport%mycorrhiza%phosphorus%ptr(ics:ice,1:nsoil_sb,iblk)       = sb_transport%mycorrhiza%phosphorus(:,:)
  sb__mem%transport%mycorrhiza%carbon13%ptr(ics:ice,1:nsoil_sb,iblk)         = sb_transport%mycorrhiza%carbon13(:,:)
  sb__mem%transport%mycorrhiza%carbon14%ptr(ics:ice,1:nsoil_sb,iblk)         = sb_transport%mycorrhiza%carbon14(:,:)
  sb__mem%transport%mycorrhiza%nitrogen15%ptr(ics:ice,1:nsoil_sb,iblk)       = sb_transport%mycorrhiza%nitrogen15(:,:)
  sb__mem%transport%microbial%carbon%ptr(ics:ice,1:nsoil_sb,iblk)            = sb_transport%microbial%carbon(:,:)
  sb__mem%transport%microbial%nitrogen%ptr(ics:ice,1:nsoil_sb,iblk)          = sb_transport%microbial%nitrogen(:,:)
  sb__mem%transport%microbial%phosphorus%ptr(ics:ice,1:nsoil_sb,iblk)        = sb_transport%microbial%phosphorus(:,:)
  sb__mem%transport%microbial%carbon13%ptr(ics:ice,1:nsoil_sb,iblk)          = sb_transport%microbial%carbon13(:,:)
  sb__mem%transport%microbial%carbon14%ptr(ics:ice,1:nsoil_sb,iblk)          = sb_transport%microbial%carbon14(:,:)
  sb__mem%transport%microbial%nitrogen15%ptr(ics:ice,1:nsoil_sb,iblk)        = sb_transport%microbial%nitrogen15(:,:)
  sb__mem%transport%residue%carbon%ptr(ics:ice,1:nsoil_sb,iblk)              = sb_transport%residue%carbon(:,:)
  sb__mem%transport%residue%nitrogen%ptr(ics:ice,1:nsoil_sb,iblk)            = sb_transport%residue%nitrogen(:,:)
  sb__mem%transport%residue%phosphorus%ptr(ics:ice,1:nsoil_sb,iblk)          = sb_transport%residue%phosphorus(:,:)
  sb__mem%transport%residue%carbon13%ptr(ics:ice,1:nsoil_sb,iblk)            = sb_transport%residue%carbon13(:,:)
  sb__mem%transport%residue%carbon14%ptr(ics:ice,1:nsoil_sb,iblk)            = sb_transport%residue%carbon14(:,:)
  sb__mem%transport%residue%nitrogen15%ptr(ics:ice,1:nsoil_sb,iblk)          = sb_transport%residue%nitrogen15(:,:)
  sb__mem%transport%residue_assoc%carbon%ptr(ics:ice,1:nsoil_sb,iblk)        = sb_transport%residue_assoc%carbon(:,:)
  sb__mem%transport%residue_assoc%nitrogen%ptr(ics:ice,1:nsoil_sb,iblk)      = sb_transport%residue_assoc%nitrogen(:,:)
  sb__mem%transport%residue_assoc%phosphorus%ptr(ics:ice,1:nsoil_sb,iblk)    = sb_transport%residue_assoc%phosphorus(:,:)
  sb__mem%transport%residue_assoc%carbon13%ptr(ics:ice,1:nsoil_sb,iblk)      = sb_transport%residue_assoc%carbon13(:,:)
  sb__mem%transport%residue_assoc%carbon14%ptr(ics:ice,1:nsoil_sb,iblk)      = sb_transport%residue_assoc%carbon14(:,:)
  sb__mem%transport%residue_assoc%nitrogen15%ptr(ics:ice,1:nsoil_sb,iblk)    = sb_transport%residue_assoc%nitrogen15(:,:)
  ! deallocate variables in helper pool structures
  deallocate(sb_transport%dom%carbon)
  deallocate(sb_transport%dom%nitrogen)
  deallocate(sb_transport%dom%phosphorus)
  deallocate(sb_transport%dom%carbon13)
  deallocate(sb_transport%dom%carbon14)
  deallocate(sb_transport%dom%nitrogen15)
  deallocate(sb_transport%dom_assoc%carbon)
  deallocate(sb_transport%dom_assoc%nitrogen)
  deallocate(sb_transport%dom_assoc%phosphorus)
  deallocate(sb_transport%dom_assoc%carbon13)
  deallocate(sb_transport%dom_assoc%carbon14)
  deallocate(sb_transport%dom_assoc%nitrogen15)
  deallocate(sb_transport%soluable_litter%carbon)
  deallocate(sb_transport%soluable_litter%nitrogen)
  deallocate(sb_transport%soluable_litter%phosphorus)
  deallocate(sb_transport%soluable_litter%carbon13)
  deallocate(sb_transport%soluable_litter%carbon14)
  deallocate(sb_transport%soluable_litter%nitrogen15)
  deallocate(sb_transport%polymeric_litter%carbon)
  deallocate(sb_transport%polymeric_litter%nitrogen)
  deallocate(sb_transport%polymeric_litter%phosphorus)
  deallocate(sb_transport%polymeric_litter%carbon13)
  deallocate(sb_transport%polymeric_litter%carbon14)
  deallocate(sb_transport%polymeric_litter%nitrogen15)
  deallocate(sb_transport%woody_litter%carbon)
  deallocate(sb_transport%woody_litter%nitrogen)
  deallocate(sb_transport%woody_litter%phosphorus)
  deallocate(sb_transport%woody_litter%carbon13)
  deallocate(sb_transport%woody_litter%carbon14)
  deallocate(sb_transport%woody_litter%nitrogen15)
  deallocate(sb_transport%fungi%carbon)
  deallocate(sb_transport%fungi%nitrogen)
  deallocate(sb_transport%fungi%phosphorus)
  deallocate(sb_transport%fungi%carbon13)
  deallocate(sb_transport%fungi%carbon14)
  deallocate(sb_transport%fungi%nitrogen15)
  deallocate(sb_transport%mycorrhiza%carbon)
  deallocate(sb_transport%mycorrhiza%nitrogen)
  deallocate(sb_transport%mycorrhiza%phosphorus)
  deallocate(sb_transport%mycorrhiza%carbon13)
  deallocate(sb_transport%mycorrhiza%carbon14)
  deallocate(sb_transport%mycorrhiza%nitrogen15)
  deallocate(sb_transport%microbial%carbon)
  deallocate(sb_transport%microbial%nitrogen)
  deallocate(sb_transport%microbial%phosphorus)
  deallocate(sb_transport%microbial%carbon13)
  deallocate(sb_transport%microbial%carbon14)
  deallocate(sb_transport%microbial%nitrogen15)
  deallocate(sb_transport%residue%carbon)
  deallocate(sb_transport%residue%nitrogen)
  deallocate(sb_transport%residue%phosphorus)
  deallocate(sb_transport%residue%carbon13)
  deallocate(sb_transport%residue%carbon14)
  deallocate(sb_transport%residue%nitrogen15)
  deallocate(sb_transport%residue_assoc%carbon)
  deallocate(sb_transport%residue_assoc%nitrogen)
  deallocate(sb_transport%residue_assoc%phosphorus)
  deallocate(sb_transport%residue_assoc%carbon13)
  deallocate(sb_transport%residue_assoc%carbon14)
  deallocate(sb_transport%residue_assoc%nitrogen15)

#endif

#if defined(Passback_sb_delta_pools)
  ! pass the values (back) to process memory
  sb__mem%delta_pools%carbon%ptr(ics:ice,1:nsoil_sb,iblk)       = sb_delta_pools%carbon(:,:)
  sb__mem%delta_pools%nitrogen%ptr(ics:ice,1:nsoil_sb,iblk)     = sb_delta_pools%nitrogen(:,:)
  sb__mem%delta_pools%phosphorus%ptr(ics:ice,1:nsoil_sb,iblk)   = sb_delta_pools%phosphorus(:,:)
  sb__mem%delta_pools%carbon13%ptr(ics:ice,1:nsoil_sb,iblk)     = sb_delta_pools%carbon13(:,:)
  sb__mem%delta_pools%carbon14%ptr(ics:ice,1:nsoil_sb,iblk)     = sb_delta_pools%carbon14(:,:)
  sb__mem%delta_pools%nitrogen15%ptr(ics:ice,1:nsoil_sb,iblk)   = sb_delta_pools%nitrogen15(:,:)
  ! deallocate variables in helper pool structures
  deallocate(sb_delta_pools%carbon)
  deallocate(sb_delta_pools%nitrogen)
  deallocate(sb_delta_pools%phosphorus)
  deallocate(sb_delta_pools%carbon13)
  deallocate(sb_delta_pools%carbon14)
  deallocate(sb_delta_pools%nitrogen15)

#endif

#if defined(Passback_sb_mycorrhiza_export)
  ! pass the values (back) to process memory
  sb__mem%mycorrhiza_export%carbon%ptr(ics:ice,1:nsoil_sb,iblk)       = sb_mycorrhiza_export%carbon(:,:)
  sb__mem%mycorrhiza_export%nitrogen%ptr(ics:ice,1:nsoil_sb,iblk)     = sb_mycorrhiza_export%nitrogen(:,:)
  sb__mem%mycorrhiza_export%phosphorus%ptr(ics:ice,1:nsoil_sb,iblk)   = sb_mycorrhiza_export%phosphorus(:,:)
  sb__mem%mycorrhiza_export%carbon13%ptr(ics:ice,1:nsoil_sb,iblk)     = sb_mycorrhiza_export%carbon13(:,:)
  sb__mem%mycorrhiza_export%carbon14%ptr(ics:ice,1:nsoil_sb,iblk)     = sb_mycorrhiza_export%carbon14(:,:)
  sb__mem%mycorrhiza_export%nitrogen15%ptr(ics:ice,1:nsoil_sb,iblk)   = sb_mycorrhiza_export%nitrogen15(:,:)
  ! deallocate variables in helper pool structures
  deallocate(sb_mycorrhiza_export%carbon)
  deallocate(sb_mycorrhiza_export%nitrogen)
  deallocate(sb_mycorrhiza_export%phosphorus)
  deallocate(sb_mycorrhiza_export%carbon13)
  deallocate(sb_mycorrhiza_export%carbon14)
  deallocate(sb_mycorrhiza_export%nitrogen15)

#endif
