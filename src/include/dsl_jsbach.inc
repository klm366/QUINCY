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

#define dsl4jsb_Use_processes USE mo_jsb_process_factory, ONLY: Get_process_memory, Get_process_config,

#define dsl4jsb_Use_memory(PROC) USE CONCATENATE3(mo_,PROC,_memory_class), ONLY: CONCATENATE3(t_,PROC,_memory)
#define dsl4jsb_Def_memory(PROC) TYPE(CONCATENATE3(t_,PROC,_memory)), POINTER :: CONCATENATE2(PROC,__mem)
#define dsl4jsb_Get_memory(PROC) CALL Get_process_memory(tile%mem, PROC, CONCATENATE2(PROC,__mem))

#define dsl4jsb_Use_config(PROC) USE CONCATENATE3(mo_,PROC,_config_class), ONLY: CONCATENATE3(t_,PROC,_config)
#define dsl4jsb_Def_config(PROC) TYPE(CONCATENATE3(t_,PROC,_config)), POINTER :: CONCATENATE2(PROC,__conf)
#define dsl4jsb_Get_config(PROC) CALL Get_process_config(model%processes(PROC)%p%config, CONCATENATE2(PROC,__conf))
#define dsl4jsb_Config(PROC) CONCATENATE2(PROC,__conf)

#define dsl4jsb_Real2D_onChunk REAL(wp), POINTER, DIMENSION(:)
#define dsl4jsb_Real3D_onChunk REAL(wp), POINTER, DIMENSION(:,:)

#define dsl4jsb_Get_var2D_onChunk(PROC,var,intent) var => CONCATENATE2(PROC,__mem)%var%ptr(ics:ice,iblk); CALL CONCATENATE2(PROC,__mem)%var%Set_inout(TOSTRING(intent))
#define dsl4jsb_Get_var3D_onChunk(PROC,var,intent) var => CONCATENATE2(PROC,__mem)%var%ptr(ics:ice,:,iblk); CALL CONCATENATE2(PROC,__mem)%var%Set_inout(TOSTRING(intent))
#define dsl4jsb_var2D_onChunk(PROC,var) CONCATENATE2(PROC,__mem)%var%ptr(ics:ice,iblk)
#define dsl4jsb_var3D_onChunk(PROC,var) CONCATENATE2(PROC,__mem)%var%ptr(ics:ice,:,iblk)

#define dsl4jsb_Get_var2d_onChunk dsl4jsb_Get_var2D_onChunk
#define dsl4jsb_Get_var3d_onChunk dsl4jsb_Get_var3D_onChunk
#define dsl4jsb_var2d_onChunk dsl4jsb_var2D_onChunk
#define dsl4jsb_var3d_onChunk dsl4jsb_var3D_onChunk

#define dsl4jsb_Real2D_onDomain REAL(wp), POINTER, DIMENSION(:,:)
#define dsl4jsb_Real3D_onDomain REAL(wp), POINTER, DIMENSION(:,:,:)

#define dsl4jsb_Get_var2D_onDomain(PROC,var,intent) var => CONCATENATE2(PROC,__mem)%var%ptr(:,:); CALL CONCATENATE2(PROC,__mem)%var%Set_inout(TOSTRING(intent))
#define dsl4jsb_Get_var3D_onDomain(PROC,var,intent) var => CONCATENATE2(PROC,__mem)%var%ptr(:,:,:); CALL CONCATENATE2(PROC,__mem)%var%Set_inout(TOSTRING(intent))
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

/*
#define VA_NUM_ARGS(...) VA_NUM_ARGS_IMPL(__VA_ARGS__, 5,4,3,2,1)
#define VA_NUM_ARGS_IMPL(_1,_2,_3,_4,_5,N,...) N

#define macro_dispatcher(func, ...) macro_dispatcher_(func, VA_NUM_ARGS(__VA_ARGS__))
#define macro_dispatcher_(func, nargs) macro_dispatcher__(func, nargs)
#define macro_dispatcher__(func, nargs) func ## nargs

#define Get_var2D_onChunk(...) macro_dispatcher(Get_var2D_onChunk, __VA_ARGS__)(__VA_ARGS__)

#define Get_var2D_onChunk2(mem,var) var => mem%var%ptr(ics:ice,iblk)
#define Get_var2D_onChunk3(mem,var,intent) var => mem%var%ptr(ics:ice,iblk); mem%var%inout=#intent
*/

