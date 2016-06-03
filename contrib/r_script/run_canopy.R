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
print(" running model simulations !")

make=read.delim('../../../util/Make.inc.local',header=F,sep="=",stringsAsFactors=F)
comp=gsub(" $","",gsub("^ ","",make[gsub(" $","",make[,1])=="COMP",2]))
system(paste("../../../",comp,"/bin/canopy.x",sep=''))

# call postprocessing_canopy.R 
system("Rscript ./postprocessing_canopy.R")
