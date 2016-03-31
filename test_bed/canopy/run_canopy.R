#!/usr/bin/env Rscript
#program to plot cexchange output 
make=read.delim('../../util/Make.inc.local',header=F,sep="=",stringsAsFactors=F)
comp=gsub(" $","",gsub("^ ","",make[gsub(" $","",make[,1])=="COMP",2]))
system(paste("../../",comp,"/bin/transient_canopy.x",sep=''))

flist = list.files(path= ".",patter="fort*")
nplot=length(flist)

pdf("test_canopy_plots.pdf", paper="a4", height=10, width=8, pointsize=12, bg="white", family="Helvetica", version="1.6")
mat=matrix(seq(1,3*nplot,1),nrow=nplot,ncol=3,byrow=T)
layout(mat)
par(mai=c(0.5,0.5,0,0))
for (i in 1:nplot) {
  dat=read.table(flist[i],header=T)
  plot(dat[,1],dat[,2],xlab=colnames(dat)[1],ylab=colnames(dat)[2],type='l')
  plot(dat[,1],dat[,3],xlab=colnames(dat)[1],ylab=colnames(dat)[3],type='l')
  plot(dat[,1],dat[,4],xlab=colnames(dat)[1],ylab=colnames(dat)[4],type='l')
}
dev.off()
