#!/usr/bin/env Rscript

# run the SLM canopy model locally, i.e. not at the HPC
# do NOT run this Script in contrib/r_script/ BUT cd into contrib/test_bed/canopy/ and run from there !!
# call postprocessing_canopy.R to create figures

# prevent running this script in contrib/r_script/
cwd<-getwd()
cwd<-strsplit(cwd, "/")
wrongFolder<-grep("r_script",cwd)
if (length(wrongFolder) > 0) {
  stop("\n >> You may not run this script from contrib/r_script/ BUT cd into contrib/test_bed/canopy/ and run from there <<") }

# message
print(" running canopy-model simulations (using land.x) !")

make=read.delim('../../../util/Make.inc.local',header=F,sep="=",stringsAsFactors=F)
comp=gsub(" $","",gsub("^ ","",make[gsub(" $","",make[,1])=="COMP",2]))
system(paste("../../../",comp,"/bin/land.x",sep=''))

# call postprocessing_canopy.R 
postprocessing<-"no"  # no/yes
if (postprocessing == "yes")
{
  system("Rscript ./postprocessing_canopy.R")
} else print(" => postprocessing_canopy.R is off by default; you can switch it on in contrib/r_script/run_canopy.R")

