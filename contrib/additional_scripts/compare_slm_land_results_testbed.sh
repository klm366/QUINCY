#!/bin/bash

## 
## this script does compare the results of two simulations limited to vegpoolC_*.txt and ecopoolCNP_*.txt
##  by:
##  1) doing a md5 checksum comparison of these output files
##  2) creating a comparison figure in the test_bed/ folder
## 

# Your options:
# give the folder names of the simulation results in the test_bed/
landOld="land_old"
landNew="land"



########################################################################################################
####### do NOT change any below this line unless you know what you are doing! ##########################
########################################################################################################


## change into the directory where this script is located; thus, you can start this script from anywhere (no need to cd into a specific directory) ##
cd "$( dirname "$0" )"

## change then to the test_bed folder using a relative path and work there
cd ../test_bed

## create a temporary folder to work in
tempFolder="./compare_slm_land_results_automatic_simple_script/"
mkdir -p "$tempFolder"

## filenames to be compared
vegPoolCFilename="vegpoolC"
ecopoolCNPFilename="ecopoolCNP"

## the names of the copy of the files to be compared
cVegPoolCOld="$tempFolder"vegpoolc_old.txt
cVegPoolCNew="$tempFolder"vegpoolc_new.txt
cEcoPoolOld="$tempFolder"ecopoolcnp_old.txt
cEcoPoolNew="$tempFolder"ecopoolcnp_new.txt

## copy the files there
cp ./"$landOld"/"$vegPoolCFilename"_*     "$cVegPoolCOld"
cp ./"$landNew"/"$vegPoolCFilename"_*     "$cVegPoolCNew"
cp ./"$landOld"/"$ecopoolCNPFilename"_*   "$cEcoPoolOld"
cp ./"$landNew"/"$ecopoolCNPFilename"_*   "$cEcoPoolNew"

## md5sum comparison
md5VegPoolCOld=$(md5sum "$cVegPoolCOld" | awk '{print $1}')
md5VegPoolCNew=$(md5sum "$cVegPoolCNew" | awk '{print $1}')
md5EcoPoolOld=$(md5sum "$cEcoPoolOld" | awk '{print $1}')
md5EcoPoolNew=$(md5sum "$cEcoPoolNew" | awk '{print $1}')

echo " "
echo " Checksum comparison of files: "

if [ "$md5VegPoolCOld" == "$md5VegPoolCNew" ]
 then
  echo "  ""$vegPoolCFilename""  are identical"
else
  echo "  ""$vegPoolCFilename""  are NOT identical"
fi

if [ "$md5EcoPoolOld" == "$md5EcoPoolNew" ]
 then
  echo "  ""$ecopoolCNPFilename""  are identical"
else
  echo "  ""$ecopoolCNPFilename""  are NOT identical"
fi

## create comparison figures using R
echo " "
echo " Creating comparison figures anyway ..."

# R script
R --slave --silent --no-save --vanilla << EOF

rVegPoolCOld<-read.table("$cVegPoolCOld", header=T, sep="")
rVegPoolCNew<-read.table("$cVegPoolCNew", header=T, sep="")
rEcoPoolOld<-read.table("$cEcoPoolOld", header=T, sep="")
rEcoPoolNew<-read.table("$cEcoPoolNew", header=T, sep="")

pdf(paste0("land_comparison_plots_","$landOld","_","$landNew",".pdf"), 
    paper="a4", height=10, width=8, pointsize=12, 
    bg="white", family="Helvetica", version="1.6")

  dat    <- rVegPoolCOld
  datNew <- rVegPoolCNew
  nvar=dim(dat)[2]
  effColNum <- ncol(dat) -1 
  mat=matrix(seq(1,2*round(effColNum/2)+2,1),nrow=round(effColNum/2)+1,ncol=2,byrow=T) 
  layout(mat)
  par(mai=c(0.4,0.5,0.14,0),oma=c(0.1,0.1,0.1,0.1)) # oma = outer margins
  for (v in 2:nvar) {
    if (v == 2) header <- paste0("$landOld",": blue  ||  ","$landNew",": red")
    else header <- NA
    if(!all(is.na(dat[,v]))){
      plot(dat[,v],ylab=colnames(dat)[v],main=header,type='l',col=NA,tck=0.03, mgp=c(1.5, 0.2, 0))
      lines(datNew[,v],lwd=0.4,col="red")
      lines(dat[,v],lwd=0.4,col="blue2")
    }
  }

  dat    <- rEcoPoolOld
  datNew <- rEcoPoolNew
  nvar=dim(dat)[2]
  effColNum <- ncol(dat) -1 
  mat=matrix(seq(1,2*round(effColNum/2)+2,1),nrow=round(effColNum/2)+1,ncol=2,byrow=T) 
  layout(mat)
  par(mai=c(0.4,0.5,0.14,0),oma=c(0.1,0.1,0.1,0.1)) # oma = outer margins
  for (v in 2:nvar) {
    if (v == 2) header <- paste0("$landOld",": blue  ||  ","$landNew",": red")
    else header <- NA
    if(!all(is.na(dat[,v]))){
      plot(dat[,v],ylab=colnames(dat)[v],main=header,type='l',col=NA,tck=0.03, mgp=c(1.5, 0.2, 0))
      lines(datNew[,v],lwd=0.4,col="red")
      lines(dat[,v],lwd=0.4,col="blue2")
    }
  }

garbage <- dev.off()

EOF
# R script end

## remove the temporary folder created above
rm -rf "$tempFolder"

## done message
echo " "
echo "  ... done!"

## fin ##
