field_name=string()
suffix=string()
frequency=option(yearly,monthly,6-hourly,default=monthly)
documentation=string()

[fetch]
[[__many__]]
longname=string()
unit=string()
[[[__many__]]]
name=string()
code=string()
submodel=string(default=None)
stream=string(default=None)

[preprocessing]
exec=string()
input=string_list
[[__many__]]

[diagnose]
exec=string()
input=string_list
[[__many__]]

[plotting]
suffix=string()
[[__many__]]
exec=string()
input=string_list
nref=option(single,multiple)
plot_title=string_list(default=None)
[[[__many__]]]

[report]
title=string(default=None)
abstract=string(default=None)
technical_teXed=string(default=None)
[[__many__]]
title=string(default=None)
name=string()
type=option(map,metric)
isbench=boolean()
benchvar=string(default=None)
benchunit=string(default=None)
ismatrix=boolean()
isvector=boolean()
abstract=string(default=None)
technical_teXed=string(default=None)
