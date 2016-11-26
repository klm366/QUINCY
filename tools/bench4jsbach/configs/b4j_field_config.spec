field_name=string()			    # USER: field id conform to id/name of diagnostic field <F> in main configuration 
suffix=string()				    # USER: output format suffix for processed data of b4j run
frequency=option('yearly','monthly','6-hourly',default='monthly')	  # USER: output frequence to use in b4j run
documentation=string()			    # USER: description

[fetch]					    # USER: section fetch for diagnostic field
[[__many__]]				    # USER: subsection <V>, <V>:= id/infix for one variable file to fetch
longname=string(default=None)		    # USER: <V>: descriptive name
unit=string()				    # USER: <V>: physical unit
[[[__many__]]]				    # USER: <V>: experiment type conform to type of experiment in main configuration
name=string()				    # USER: <V>: name
code=string()				    # USER: <V>: code 
submodel=string(default=None)		    # USER: <V>: submodel
stream=string(default=None)		    # USER: <V>: stream		  

[preprocessing]				    # USER: section preprocessing for diagnostic field
exec=string()				    # USER: executable script to call in b4j run
input=string_list
[[__many__]]				    # USER: id/infix for output file written by the script  

[diagnose]				    # USER: section diagnose for diagnostic field 
exec=string()				    # USER: executable script to call in b4j run       
input=string_list
[[__many__]]				    # USER: id/infix <D> for output file written by the script	       

[plotting]				    # USER: section plotting for diagnostic field
suffix=string()				    # USER: output format's suffix
[[__many__]]				    # USER: subsection <O>, <O>:= id of plot script
exec=string()				    # USER: <O>: executable script to call in b4j run
input=string_list
nref=option('single','multiple','None')	    # USER: <O>: as indicated by plot: comparison with single or all reference experiments 
plot_title=string_list(default=None)	    # USER: <O>: title string for plot
[[[__many__]]]				    # USER: <O>: id/infix for output file written by the script	    

[report]				    # USER: section report for diagnostic field
title=string(default=None)		    # USER: title for the diagnostic field's section in the report
abstract=string(default=None)		    # USER: scientific description of the diagnostic field
technical_teXed=string(default=None)	    # USER: technical description of the diagnostic field 
[[__many__]]				    # USER: subsection <SD> for subsection for a diagnostic method, key <SD> conform to key <D>
title=string(default=None)		    # USER: <SD>: title of report subsection
name=string()				    # USER: <SD>: name of diagnostic method
type=option('map','metric','barchart')	    # USER: <SD>: type of plot
isbench=boolean()			    # USER: <SD>: diagnostic method is used as benchmark
benchvar=string(default=None)		    # USER: <SD>: name of displayed variable in plotfile
benchunit=string(default=None)		    # USER: <SD>: physical unit
ismatrix=boolean()			    # USER: <SD>: comparison with all references
isvector=boolean()			    # USER: <SD>: comparison pairwise
abstract=string(default=None)		    # USER: <SD>: scientific description of diagnostic method
technical_teXed=string(default=None)	    # USER: <SD>: technical description of diagnostic method 


