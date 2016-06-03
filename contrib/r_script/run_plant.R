#!/usr/bin/env Rscript

# program to run SLM locally, i.e. not at the HPC
# do NOT run this Script in contrib/r_script/ BUT cd into contrib/test_bed/plant/ and run from there !!
# postprocessign of output data (e.g. plot plant model output): postprocessing_plant.R 

# prevent running this script in contrib/r_script/
cwd<-getwd()
cwd<-strsplit(cwd, "/")
wrongFolder<-grep("r_script",cwd)
if (length(wrongFolder) > 0) {
  stop("\n >> You may not run this script from contrib/r_script/ BUT cd into contrib/test_bed/plant/ and run from there <<") }

# string function to remove trailing/leading blanks
trim <- function (x) gsub("^\\s+|\\s+$", "", x)

# message
print(" running model simulations !")

# run plant.x
make=read.delim('../../../util/Make.inc.local',header=F,sep="=",stringsAsFactors=F)
comp=trim(make[trim(make[,1])=="COMP",2])
system(paste("../../../",comp,"/bin/plant.x",sep=''))

# call postprocessing_plant.R 
system("Rscript ./postprocessing_plant.R")
