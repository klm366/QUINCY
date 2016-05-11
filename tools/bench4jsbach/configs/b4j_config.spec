bench_name=string()
author=string()
reporting=boolean()
land_sea_mask=string()

[experiments]
[[__many__]]
name=string()
type=option('CMIP5','MPIESM1')
path=string()
iscan=boolean()
isref=boolean()

[fields]
[[__many__]]
config_fname=string()
