PROGRAM example

  dsl4jsb_Use_processes SOIL_, &
    & VEG_
  dsl4jsb_Use_memory(SoIL_)
  dsl4jsb_Def_memory(SOIL_)
  dsl4jsb_get_memory(SOIL_) ! comment
  dsl4jsb_Use_config(SOIL_)
  dsl4jsb_Def_config(SOIL_)
  dsl4jsb_Get_config(SOIL_)
  x = dsl4jsb_Config(SOIL_)%par1 * dsl4jsb_Config(SOIL_)%par2 * ws(:) ! comment
  CALL sub(dsl4jsb_Lctlib , arg2, dsl4jsb_Lctlib_param(abc)) ! This doesn't work!!!
  CALL sub(dsl4jsb_Lctlib, arg2)
  CALL sub(dsl4jsb_Lctlib)
  CALL sub(dsl4jsb_Lctlib_param(abc), dsl4jsb_Lctlib_param(def))
  dsl4jsb_Real2D_onChunk :: ws
  dsl4jsb_Get_var2D_onChunk(SOIL_,ws)
  dsl4jsb_Get_var2D_onChunk(SOIL_, ws)
  dsl4jsb_var2d_onChunk(SOIL_,ws)
  dsl4jsb_Real2D_onDomain :: ws
  dsl4jsb_Get_var2D_onDomain(SOIL_,ws)
  dsl4jsb_var2D_onDomain(SOIL_,ws)
  dsl4jsb_var2D_ptr(SOIL_,ws)
  x = dsl4jsb_Config(SOIL_)%par1 * dsl4jsb_var2D_ptr(SOIL_, ws)(:) ! comment
  dsl4jsb_Real3D_onChunk :: ws_l
  dsl4jsb_Get_var3D_onChunk(SOIL_,ws_l)
  dsl4jsb_var3D_onChunk(SOIL_,ws_l)
  dsl4jsb_Real3D_onDomain :: ws_l
  dsl4jsb_Get_var3D_onDomain(SOIL_,ws_l)
  dsl4jsb_var3D_onDomain(SOIL_,ws_l)
  dsl4jsb_var3D_ptr(SOIL_,ws_l)
  dsl4jsb_Aggregate_onChunk(SOIL_,ws, aggregator)

  ! Test abort on error
  dsl4jsb_Def_menory(SOIL_)

END PROGRAM example
