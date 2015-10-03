#define IDENTITY(x) x
#define CONCATENATE(a,b,c) IDENTITY(a)IDENTITY(b)IDENTITY(c)

#define Use_processes USE mo_jsb_process_factory, ONLY: Get_process_memory, Get_process_config,

#define Use_memory(PROC) USE CONCATENATE(mo_,PROC,_memory_class), ONLY: CONCATENATE(t_,PROC,_memory)
#define Def_memory(PTR, PROC) TYPE(CONCATENATE(t_,PROC,_memory)), POINTER :: PTR => NULL()
         !!!#define Def_memory(PTR, PROC) TYPE(t_##PROC##_memory), POINTER :: PTR => NULL()
#define Get_memory(PTR, PROC) CALL Get_process_memory(tile%mem, PROC, PTR)

#define Use_config(PROC) USE CONCATENATE(mo_,PROC,_config_class), ONLY: CONCATENATE(t_,PROC,_config)
#define Def_config(PTR, PROC) TYPE(CONCATENATE(t_,PROC,_config)), POINTER :: PTR => NULL()
#define Get_config(PTR, PROC) CALL Get_process_config(model%processes(PROC)%config, PTR)

#define Def_real2d_onChunk(var) REAL(wp), POINTER :: var(:)
#define Def_real3d_onChunk(var) REAL(wp), POINTER :: var(:,:)
