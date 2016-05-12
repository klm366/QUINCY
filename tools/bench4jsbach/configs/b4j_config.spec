bench_name=string()
author=string()
git_repository=string_list(default=None)
git_branch=string(default=None)
git_sha=string(default=None)
git_status=option('clean','modified',None)
redmine_url=string()
redmine_project=string()
redmine_wikipage=sting()
host_name=string(default=None)
host_ipaddress=string(default=None)
b4j_base=string()
timestamp=string()

redmine_user=sting()
redmine_key=string()
institute=string()
department=string()
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
