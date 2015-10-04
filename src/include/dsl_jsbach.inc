#define dsl$jsb

#define IDENTITY(x) x
#define CONCATENATE2(a,b) IDENTITY(a)IDENTITY(b)
#define CONCATENATE3(a,b,c) IDENTITY(a)IDENTITY(b)IDENTITY(c)

#define Use_processes USE mo_jsb_process_factory, ONLY: Get_process_memory, Get_process_config,

#define Use_memory(PROC) USE CONCATENATE3(mo_,PROC,_memory_class), ONLY: CONCATENATE3(t_,PROC,_memory)
#define Def_memory(PTR, PROC) TYPE(CONCATENATE3(t_,PROC,_memory)), POINTER :: PTR => NULL()
#define Get_memory(PTR, PROC) CALL Get_process_memory(tile%mem, PROC, PTR)

#define Use_config(PROC) USE CONCATENATE3(mo_,PROC,_config_class), ONLY: CONCATENATE3(t_,PROC,_config)
#define Def_config(PTR, PROC) TYPE(CONCATENATE3(t_,PROC,_config)), POINTER :: PTR => NULL()
#define Get_config(PTR, PROC) CALL Get_process_config(model%processes(PROC)%config, PTR)

#define onChunk_Real2D REAL(wp), POINTER, DIMENSION(:)
#define onChunk_Real3D REAL(wp), POINTER, DIMENSION(:,:)

#define Get_var2D_onChunk(mem,var,intent) var => mem%var%ptr(ics:ice,iblk); CONCATENATE3(CALL mem%var%Set_to_,intent,())
#define Get_var3D_onChunk(mem,var,intent) var => mem%var%ptr(ics:ice,:,iblk); CONCATENATE3(CALL mem%var%Set_to_,intent,())
#define Get_var2d_onChunk Get_var2D_onChunk
#define Get_var3d_onChunk Get_var3D_onChunk

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

