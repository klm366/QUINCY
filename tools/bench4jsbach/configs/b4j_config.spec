bench_name=string()
author=string()
institute=string()
department=string()
git_repository=string_list(default=None)#string(default=None)
git_branch=string(default=None)
git_sha=string(default=None)
git_status=option('clean','modified',None)
hostname=string()
ip_address=string(default=None)
reporting=boolean()
land_sea_mask=string()

[experiments]
[[__many__]]
name=string()
type=option('CMIP5','MPIESM1')
path=string(default=None)
iscan=boolean()
isref=boolean()   

[fields]
[[__many__]]
config_filename=string()
