PROGRAM example

  dsl4jsb_Use_processes SOIL, &
    & VEG
  dsl4jsb_Use_memory(SoIL)
  dsl4jsb_Def_memory(SOIL)
  dsl4jsb_get_memory(SOIL) ! comment
  dsl4jsb_Use_config(SOIL)
  dsl4jsb_Def_config(SOIL)
  dsl4jsb_Get_config(SOIL)
  x = dsl4jsb_Config(SOIL)%par1 * dsl4jsb_Config(SOIL)%par2 * ws(:) ! comment
  dsl4jsb_Real2D_onChunk :: ws
  dsl4jsb_Get_var2D_onChunk(SOIL,ws,out)
  dsl4jsb_var2d_onChunk(SOIL,ws)
  dsl4jsb_Real2D_onDomain :: ws
  dsl4jsb_Get_var2D_onDomain(SOIL,ws,out)
  dsl4jsb_var2D_onDomain(SOIL,ws)
  dsl4jsb_var2D_ptr(SOIL,ws)
  dsl4jsb_Real3D_onChunk :: ws_l
  dsl4jsb_Get_var3D_onChunk(SOIL,ws_l,out)
  dsl4jsb_var3D_onChunk(SOIL,ws_l)
  dsl4jsb_Real3D_onDomain :: ws_l
  dsl4jsb_Get_var3D_onDomain(SOIL,ws_l,out)
  dsl4jsb_var3D_onDomain(SOIL,ws_l)
  dsl4jsb_var3D_ptr(SOIL,ws_l)

  ! Test abort on error
  dsl4jsb_Def_menory(SOIL)

END PROGRAM example
