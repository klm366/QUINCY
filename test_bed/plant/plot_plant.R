# program to plot plant model output 
#  changed pattern to "*.txt" in line 7 [Jan Engel 20160315]
make=read.delim('../../util/Make.inc.local',header=F,sep="=",stringsAsFactors=F)
comp=gsub(" $","",gsub("^ ","",make[gsub(" $","",make[,1])=="COMP",2]))
system(paste("../../",comp,"/bin/plant.x",sep=''))

flist = list.files(path= ".",patter="*.txt")
nplot=length(flist)

pdf('plant_plots.pdf')
for (i in 1:nplot) {
  dat=read.table(flist[i],header=T)
  nvar=dim(dat)[2]
  mat=matrix(seq(1,2*round(nvar/2)+2,1),nrow=round(nvar/2)+1,ncol=2,byrow=T) 
  layout(mat)
  par(mai=c(0.5,1,0,0))
  for (v in 1:nvar) plot(dat[,v],ylab=colnames(dat)[v],type='l')
}
dev.off()
