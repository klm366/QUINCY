#!/usr/bin/env Rscript

# postprocessing of output data from QUINCY test_canopy model
# called from contrib/test_bed/run-slm_test-bed.sh

# create overall file list
flist = list.files(path= ".",patter="fort*")
nplot=length(flist)

# plot data
pdf("test_canopy_plots.pdf", paper="a4", height=10, width=8, pointsize=12, bg="white", family="Helvetica", version="1.6")
  mat=matrix(seq(1,3*nplot,1),nrow=nplot,ncol=3,byrow=T)
  layout(mat)
  par(mai=c(0.3,0.3,0.08,0.06))
  for (i in 1:nplot) {
    dat=read.table(flist[i],header=T)
    plot(dat[,1],dat[,2], xlab=colnames(dat)[1], ylab=colnames(dat)[2], type='l', tck=0.02, mgp=c(1.2, 0.14, 0))
    plot(dat[,1],dat[,3], xlab=colnames(dat)[1], ylab=colnames(dat)[3], type='l', tck=0.02, mgp=c(1.2, 0.14, 0))
    plot(dat[,1],dat[,4], xlab=colnames(dat)[1], ylab=colnames(dat)[4], type='l', tck=0.02, mgp=c(1.2, 0.14, 0))
  }

garbage <- dev.off()  # run dev.off() but suppress the notification it usually provides
