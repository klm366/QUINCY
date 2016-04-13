bench_name=string()
bench_configspec_fname=string()
author=string(default=None)
reporting=boolean(default=False)
land_sea_mask=string(default=None)

[experiments]
[[__many__]]
name=string()
type=option('CMIP5','MPIESM1',default=None)
path=string(default=None)
iscan=boolean(default=False)
isref=boolean(default=False)

[fields]
field_configspec_fname=string()
[[__many__]]
config_fname=string()
