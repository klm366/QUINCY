bench_name=string()                         # AUTO: name of b4j run
author=string()                             # AUTO: username
git_repository=string_list(default=None)    # AUTO: b4j repository
git_branch=string(default=None)             # AUTO: 
git_sha=string(default=None)		    # AUTO: 
git_status=option('clean','modified',None)  # AUTO: 
redmine_url=string()			    # AUTO: b4j redmine site 
redmine_project=string()		    # AUTO: 
redmine_wikipage=string()		    # AUTO: 
host_name=string(default=None)		    # AUTO: host of b4j run
host_ipaddress=string(default=None)	    # AUTO: 
b4j_base=string()			    # AUTO: 
timestamp=string()			    # AUTO: b4j run config file creation time   

redmine_user=string()			    # USER: user redmine login
redmine_key=string()			    # USER: to find in redmine, "My account", right panel, "Api access key"
institute=string()			    # USER: user institute
department=string()			    # USER: user department
reporting=boolean()			    # USER: write report for b4j run, default "True"

[experiments]				    # USER: section experiments for b4j run
[[__many__]]				    # USER: subsection <A>, <A>:= id for one experiment
name_orig=string()			    # USER: <A>: original infix/id
type=option('CMIP5','MPIESM1')		    # USER: <A>: experiment type
path=string(default=None)		    # USER: <A>: location
host=string(default=None)		    # USER: <A>: 
iscan=boolean()				    # USER: <A>: candidate role
isref=boolean() 			    # USER: <A>: reference role
startyear=string()			    # USER: <A>: startyear of <A> in b4j run
endyear=string()			    # USER: <A>: endyear of <A> in b4j run
initialyear=string(default=None)	    # USER: <A>: initialyear of experiment 
finalyear=string(default=None)  	    # USER: <A>: finalyear of experiment 

[fields]				    # USER: section diagnostic fields for b4j run
[[__many__]]				    # USER: subsection <F>, <F>:= id for one diagnostic field
config_filename=string()                    # USER: <F>: configuration file's relative path from b4j_base or filename in <b4j_base>/prod/scripts                     
