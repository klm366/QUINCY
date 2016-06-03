#!/usr/bin/env Rscript

# postprocessing of output data (e.g. plot canopy model output)
# called from contrib/site_level/run_slm.sh (HPC simulations) or run_canopy.R (local SLM runs)

# message
print(" postprocessing SLM output")

flist = list.files(path= ".",patter="fort*")
nplot=length(flist)

pdf("canopy_plots.pdf", paper="a4", height=10, width=8, pointsize=12, bg="white", family="Helvetica", version="1.6")
mat=matrix(seq(1,3*nplot,1),nrow=nplot,ncol=3,byrow=T)
layout(mat)
par(mai=c(0.5,0.5,0,0))
for (i in 1:nplot) {
  dat=read.table(flist[i],header=T)
  plot(dat[,1],dat[,2],xlab=colnames(dat)[1],ylab=colnames(dat)[2],type='l')
  plot(dat[,1],dat[,3],xlab=colnames(dat)[1],ylab=colnames(dat)[3],type='l')
  plot(dat[,1],dat[,4],xlab=colnames(dat)[1],ylab=colnames(dat)[4],type='l')
}
garbage <- dev.off()  # run dev.off() but suppress the notification it usually provides

# message
print(" finished postprocessing SLM output !")
