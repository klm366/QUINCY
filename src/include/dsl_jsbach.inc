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

#define dsl4jsb_Use_memory(PROC) USE CONCATENATE3(mo_,PROC,_memory_class), ONLY: CONCATENATE3(t_,PROC,_memory)
#define dsl4jsb_Def_memory(PROC) TYPE(CONCATENATE3(t_,PROC,_memory)), POINTER :: CONCATENATE2(PROC,__mem)
#define dsl4jsb_Get_memory(PROC) SELECT TYPE (mm => tile%mem(PROC)%p) ; TYPE IS( CONCATENATE3(t_,PROC,_memor)y) ; CONCATENATE2(PROC,__mem => mm ; END SELECT)

#define dsl4jsb_Use_config(PROC) USE CONCATENATE3(mo_,PROC,_config_class), ONLY: CONCATENATE3(t_,PROC,_config)
#define dsl4jsb_Def_config(PROC) TYPE(CONCATENATE3(t_,PROC,_config)), POINTER :: CONCATENATE2(PROC,__conf)
#define dsl4jsb_Get_config(PROC) SELECT TYPE (cc => model%processes(PROC)%p%config) ; TYPE IS( CONCATENATE3(t_,PROC,_confi)g) ; CONCATENATE2(PROC,__conf => cc ; END SELECT)

#define dsl4jsb_Config(PROC) CONCATENATE2(PROC,__conf)

#define dsl4jsb_Lctlib(PARAM) CONCATENATE2(model%lctlib%PARAM,(tile%lcts(1)%lib_id))

#define dsl4jsb_Real2D_onChunk REAL(wp), POINTER, DIMENSION(:)
#define dsl4jsb_Real3D_onChunk REAL(wp), POINTER, DIMENSION(:,:)

#define dsl4jsb_Get_var2D_onChunk(PROC,var) var => CONCATENATE2(PROC,__mem)%var%ptr(ics:ice,iblk)

#define dsl4jsb_Get_var3D_onChunk(PROC,var) var => CONCATENATE2(PROC,__mem)%var%ptr(ics:ice,:,iblk)

#define dsl4jsb_var2D_onChunk(PROC,var) CONCATENATE2(PROC,__mem)%var%ptr(ics:ice,iblk)
#define dsl4jsb_var3D_onChunk(PROC,var) CONCATENATE2(PROC,__mem)%var%ptr(ics:ice,:,iblk)

#define dsl4jsb_Get_var2d_onChunk dsl4jsb_Get_var2D_onChunk
#define dsl4jsb_Get_var3d_onChunk dsl4jsb_Get_var3D_onChunk
#define dsl4jsb_var2d_onChunk dsl4jsb_var2D_onChunk
#define dsl4jsb_var3d_onChunk dsl4jsb_var3D_onChunk

#define dsl4jsb_Real2D_onDomain REAL(wp), POINTER, DIMENSION(:,:)
#define dsl4jsb_Real3D_onDomain REAL(wp), POINTER, DIMENSION(:,:,:)

#define dsl4jsb_Get_var2D_onDomain(PROC,var) var => CONCATENATE2(PROC,__mem)%var%ptr(:,:)

#define dsl4jsb_Get_var3D_onDomain(PROC,var) var => CONCATENATE2(PROC,__mem)%var%ptr(:,:,:)

#define dsl4jsb_var2D_onDomain(PROC,var) CONCATENATE2(PROC,__mem)%var%ptr(:,:)
#define dsl4jsb_var3D_onDomain(PROC,var) CONCATENATE2(PROC,__mem)%var%ptr(:,:,:)

#define dsl4jsb_Get_var2d_onDomain dsl4jsb_Get_var2D_onDomain
#define dsl4jsb_Get_var3d_onDomain dsl4jsb_Get_var3D_onDomain
#define dsl4jsb_var2d_onDomain dsl4jsb_var2D_onDomain
#define dsl4jsb_var3d_onDomain dsl4jsb_var3D_onDomain

#define dsl4jsb_var2D_ptr(PROC,var) CONCATENATE2(PROC,__mem)%var%ptr
#define dsl4jsb_var3D_ptr(PROC,var) CONCATENATE2(PROC,__mem)%var%ptr
#define dsl4jsb_var2d_ptr dsl4jsb_var2D_ptr
#define dsl4jsb_var3d_ptr dsl4jsb_var3D_ptr

#define dsl4jsb_Aggregate_onChunk(PROC,var,agg) CALL agg%Aggregate(tile, CONCATENATE2(PROC,__mem)%var,ics,ice,iblk)
