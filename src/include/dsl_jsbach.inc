#if defined(__GFORTRAN__)
#define STRINGIFY(a) "a"
#else
#if defined(NAGFOR)
#define STRINGIFY(a)
#else
#define STRINGIFY(a) #a
#endif
#endif
#define TOSTRING(x) STRINGIFY(x)

#define IDENTITY(x) x
#define CONCATENATE2(a,b) IDENTITY(a)IDENTITY(b)
#define CONCATENATE3(a,b,c) IDENTITY(a)IDENTITY(b)IDENTITY(c)

#define dsl4jsb_Use_processes USE mo_jsb_process_class, ONLY: 

#define dsl4jsb_Use_memory(PROC) USE CONCATENATE3(mo_,PROC,memory_class), ONLY: CONCATENATE3(t_,PROC,memory)
#define dsl4jsb_Def_memory(PROC) TYPE(CONCATENATE3(t_,PROC,memory)), POINTER :: CONCATENATE2(PROC,_mem)
#define dsl4jsb_Get_memory(PROC) SELECT TYPE (mm => tile%mem(PROC)%p) ; TYPE IS( CONCATENATE3(t_,PROC,memor)y) ; CONCATENATE2(PROC,_mem => mm ; END SELECT)

#define dsl4jsb_Use_config(PROC) USE CONCATENATE3(mo_,PROC,config_class), ONLY: CONCATENATE3(t_,PROC,config)
#define dsl4jsb_Def_config(PROC) TYPE(CONCATENATE3(t_,PROC,config)), POINTER :: CONCATENATE2(PROC,_conf)
#define dsl4jsb_Get_config(PROC) SELECT TYPE (cc => model%processes(PROC)%p%config) ; TYPE IS( CONCATENATE3(t_,PROC,confi)g) ; CONCATENATE2(PROC,_conf => cc ; END SELECT)

#define dsl4jsb_Config(PROC) CONCATENATE2(PROC,_conf)

#define dsl4jsb_Lctlib_param(PARAM) CONCATENATE2(model%lctlib,(tile%lcts(1)%lib_id)%PARAM)
#define dsl4jsb_Lctlib CONCATENATE2(model%lctlib,(tile%lcts(1)%lib_id))

#define dsl4jsb_Real2D_onChunk REAL(wp), POINTER, DIMENSION(:)
#define dsl4jsb_Real3D_onChunk REAL(wp), POINTER, DIMENSION(:,:)

#define dsl4jsb_Get_var2D_onChunk(PROC,var) var => CONCATENATE2(PROC,_mem)%var%ptr(ics:ice,iblk)

#define dsl4jsb_Get_var3D_onChunk(PROC,var) var => CONCATENATE2(PROC,_mem)%var%ptr(ics:ice,:,iblk)

#define dsl4jsb_var2D_onChunk(PROC,var) CONCATENATE2(PROC,_mem)%var%ptr(ics:ice,iblk)
#define dsl4jsb_var3D_onChunk(PROC,var) CONCATENATE2(PROC,_mem)%var%ptr(ics:ice,:,iblk)

#define dsl4jsb_Get_var2d_onChunk dsl4jsb_Get_var2D_onChunk
#define dsl4jsb_Get_var3d_onChunk dsl4jsb_Get_var3D_onChunk
#define dsl4jsb_var2d_onChunk dsl4jsb_var2D_onChunk
#define dsl4jsb_var3d_onChunk dsl4jsb_var3D_onChunk

#define dsl4jsb_Real2D_onDomain REAL(wp), POINTER, DIMENSION(:,:)
#define dsl4jsb_Real3D_onDomain REAL(wp), POINTER, DIMENSION(:,:,:)

#define dsl4jsb_Get_var2D_onDomain(PROC,var) var => CONCATENATE2(PROC,_mem)%var%ptr(:,:)

#define dsl4jsb_Get_var3D_onDomain(PROC,var) var => CONCATENATE2(PROC,_mem)%var%ptr(:,:,:)

#define dsl4jsb_var2D_onDomain(PROC,var) CONCATENATE2(PROC,_mem)%var%ptr(:,:)
#define dsl4jsb_var3D_onDomain(PROC,var) CONCATENATE2(PROC,_mem)%var%ptr(:,:,:)

#define dsl4jsb_Get_var2d_onDomain dsl4jsb_Get_var2D_onDomain
#define dsl4jsb_Get_var3d_onDomain dsl4jsb_Get_var3D_onDomain
#define dsl4jsb_var2d_onDomain dsl4jsb_var2D_onDomain
#define dsl4jsb_var3d_onDomain dsl4jsb_var3D_onDomain

#define dsl4jsb_var2D_ptr(PROC,var) CONCATENATE2(PROC,_mem)%var%ptr
#define dsl4jsb_var3D_ptr(PROC,var) CONCATENATE2(PROC,_mem)%var%ptr
#define dsl4jsb_var2d_ptr dsl4jsb_var2D_ptr
#define dsl4jsb_var3d_ptr dsl4jsb_var3D_ptr

#define dsl4jsb_Aggregate_onChunk(PROC,var,agg) CALL agg%Aggregate(tile, CONCATENATE2(PROC,_mem)%var,ics,ice,iblk,routine)

#define dsl4slm_Use_memory_a_veg USE mo_veg_memory_class, ONLY: t_veg_element_a, t_veg_pool_a
#define dsl4slm_Use_memory_a_sb USE mo_sb_memory_class, ONLY: t_sb_element_a, t_sb_pool_a

#define dsl4slm_get_vegpool_var(poolN,elemN,varN) allocate(CONCATENATE2(veg_,poolN)%elemN%varN, source=reshape(veg__mem%poolN%elemN%varN%ptr(ics:ice,iblk),[nc]))
#define dsl4slm_passback_vegpool_var(poolN,elemN,varN) veg__mem%poolN%elemN%varN%ptr(ics:ice,iblk) = CONCATENATE2(veg_,poolN)%elemN%CONCATENATE2(varN,(:)) ; deallocate(CONCATENATE2(veg_,poolN)%elemN%varN)

#define dsl4slm_get_sbpool_var(poolN,elemN,varN) allocate(CONCATENATE2(sb_,poolN)%elemN%varN, source=reshape(sb__mem%poolN%elemN%varN%ptr(ics:ice,1:nsoil_sb,iblk),[nc]))
#define dsl4slm_passback_sbpool_var(poolN,elemN,varN) sb__mem%poolN%elemN%varN%ptr(ics:ice,1:nsoil_sb,iblk) = CONCATENATE2(sb_,poolN)%elemN%CONCATENATE2(varN,(:,:)) ; deallocate(CONCATENATE2(sb_,poolN)%elemN%varN)
