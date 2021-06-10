#!/usr/bin/env Rscript

#
# this script is called by ./compare_two_quincy_experiments.sh in the same directory
# it creates comparison plots for results of two simulation runs into the ./figs directory next to this script
#

## get arguments
args    <- commandArgs(TRUE)
workDir <- args[1]
oldID   <- args[2]
newID   <- args[3]
outFPer <- args[4]

## set working dir
setwd(workDir)


###########################################################
## read simulation output files
###########################################################
# new simulation (black lines in the figs)
flist_new = Sys.glob("new-output/vegpool*.txt") 
flist_new = c(flist_new,Sys.glob("new-output/soilpool*.txt"))
flist_new = c(flist_new,Sys.glob("new-output/veg_diagnostics*.txt"))
flist_new = c(flist_new,Sys.glob("new-output/soil_diagnostics*.txt"))
flist_new = c(flist_new,Sys.glob("new-output/ecopool*.txt"))
flist_new = c(flist_new,Sys.glob("new-output/vegflux*.txt"))
flist_new = c(flist_new,Sys.glob("new-output/soilflux*.txt"))
flist_new = c(flist_new,Sys.glob("new-output/sb_massBalance*.txt"))
flist_new = c(flist_new,Sys.glob("new-output/atm_forcing*.txt"))
# old simulation (orange lines in the figs)
flist_old = Sys.glob("old-output/vegpool*.txt") 
flist_old = c(flist_old,Sys.glob("old-output/soilpool*.txt"))
flist_old = c(flist_old,Sys.glob("old-output/veg_diagnostics*.txt"))
flist_old = c(flist_old,Sys.glob("old-output/soil_diagnostics*.txt"))
flist_old = c(flist_old,Sys.glob("old-output/ecopool*.txt"))
flist_old = c(flist_old,Sys.glob("old-output/vegflux*.txt"))
flist_old = c(flist_old,Sys.glob("old-output/soilflux*.txt"))
flist_old = c(flist_old,Sys.glob("old-output/sb_massBalance*.txt"))
flist_old = c(flist_old,Sys.glob("old-output/atm_forcing*.txt"))

## get & check number output files
nFilesO <- length(flist_old)
nFilesN <- length(flist_new)
if (nFilesO != nFilesN)
{
  print("EXIT! The number of output files considered, does differ between simulations.")
  print("       You could change the files used by the results_comparison_two_exp.r script")
  print(" ")
  print(paste0("old: ", flist_old))
  print(paste0("new: ", flist_new))
  quit(save = "no", status = 1)
}
nFiles <- nFilesN


###########################################################
## plot data of old and new
###########################################################
print("<< creating comparison plots >>")
pdf("figs/comp_results_two_exp.pdf", paper="a4", height=12.4, width=7.6, pointsize=12, 
    bg="white", family="Helvetica", version="1.6")
  for (i in 1:nFiles)
  {
    # read files and print filename
    print(flist_new[i])
    dat_new <- read.table(flist_new[i], header=T, sep="")
    print(flist_old[i])
    dat_old <- read.table(flist_old[i], header=T, sep="")
    # number of columns
    colNumN  <- ncol(dat_new)
    colNumO  <- ncol(dat_old)

    # get the file name without directory
    fNameN  <- gsub("new-output/", "", flist_new[i])
    fNameN  <- gsub(".txt", "", fNameN)
    
    # check whether column-names match (i.e., names & sequence of variables)
    # skip this file if column names of 'old' and 'new' are no match (using the R statement 'next')
    if (colNumN != colNumO)
    {
      print(paste0("   Number of variables differ between old and new output file: ", fNameN, " - skipping this file."))
      next
    }
    if ( all(colnames(dat_new) != colnames(dat_old)) )
    {
      print(paste0("   Names of variables differ between old and new output file: ", fNameN, " - skipping this file."))
      next
    }
    
    # collect information to distribute figs over multiple pages if more than 14
    # the number of variables is the 'number of columns - 1', as the 1st column is always 'time'
    varNum   <- colNumN - 1
    maxFigPg <- 14
    if (varNum < (maxFigPg + 1))
    {
      # page numbers
      pagesNum <- 1
      # column numbers of output file per page 
      colNumStarts <- 2
      colNumEnds   <- colNumN
    } else{
      # page numbers
      pagesNum <- as.integer(varNum / maxFigPg)  # number of pages minus 1 | == pages with 14 figs each
                                                 # last page may has less figs
      # test if there is a last page with /= 14 figs, and increase pagesNum if true 
      if ((varNum / maxFigPg) > pagesNum)
      {
        pagesNum <- pagesNum + 1
      }
      # column numbers of output file per page 
      colNumLast   <- colNumN
      colNumStarts <- seq(2, colNumLast, maxFigPg)    # e.g., seq(2,37,14) = 2 16 30
      colNumEnds   <- colNumStarts[2:length(colNumStarts)] - 1
      colNumEnds   <- c(colNumEnds, colNumLast)      
    }

    # create figures per page
    for (p in 1:pagesNum)
    {
      # get the effective (actual) number of variables at this page
      effVarNum <- colNumEnds[p] - colNumStarts[p] + 1

      # page layout
      # see also layout information in below 'plot()' command
      # set the minimum of effVarNum to 4 (so that figures get not too big if their number is)
      mat  <- matrix(seq(1,2*round(max(effVarNum,4)/2)+2,1), nrow=round(max(effVarNum,4)/2)+1, ncol=2, byrow=T) 
      layout(mat)
      par(mai=c(0.32, 0.4, 0.14, 0.02), oma=c(0.1, 0.1, 2.0, 0.1)) # oma = outer margins  || c(bottom, left, top, right)

      # loop over columns (variables) and create figs   
      idInfo1 <- paste0("orange -> old: ", oldID)
      idInfo2 <- paste0("black -> new: ", newID)  
      for (v in colNumStarts[p]:colNumEnds[p])
      {
        # plot simulation ID info at uppermost figures per page
        idInfo <- ""
        if(v == colNumStarts[p])     idInfo <- idInfo1
        if(v == colNumStarts[p] + 1) idInfo <- idInfo2
        # set y-axis range to default values
        yMin <- -1
        yMax <-  1
        # test if all values are NA
        if(all(is.na(dat_new[,v])) | all(is.na(dat_old[,v]))){
          print(paste0("  creating a green zero-line diff-figure for ", colnames(dat_new)[v], " because all values of either old or new are NA"))
          plot(1:nrow(dat_new), rep(0, nrow(dat_new)),
               type='l',lwd=0.4, col="forestgreen", ylim=c(yMin,yMax), ylab=colnames(dat_new)[v], xlab=outFPer, main=idInfo, 
               tck=0.03, mgp=c(1.5, 0.2, 0))
          if(v == colNumStarts[p]) mtext(fNameN, side=3, srt = 90, cex=1.0, line=1.8, col="red3")
        } else if(all(is.na(dat_new[,v])) & all(is.na(dat_old[,v]))){
          print(paste0("  creating a green zero-line diff-figure for ", colnames(dat_new)[v], " because all values of either old and new are NA"))
          plot(1:nrow(dat_new), rep(0, nrow(dat_new)),
               type='l',lwd=0.4, col="forestgreen", ylim=c(yMin,yMax), ylab=colnames(dat_new)[v], xlab=outFPer, main=idInfo, 
               tck=0.03, mgp=c(1.5, 0.2, 0))
          if(v == colNumStarts[p]) mtext(fNameN, side=3, srt = 90, cex=1.0, line=1.8, col="red3")
        } else{
          # y-axis min & max
          yMin <- min( c(min(dat_new[,v], na.rm=T), min(dat_old[,v], na.rm=T)) )
          yMax <- max( c(max(dat_new[,v], na.rm=T), max(dat_old[,v], na.rm=T)) )
          # plot
          plot(dat_new[,v], 
               type='l',lwd=0.4, col="black", ylim=c(yMin,yMax), ylab=colnames(dat_new)[v], xlab=outFPer, main=idInfo, 
               tck=0.03, mgp=c(1.5, 0.2, 0))
          if(v == colNumStarts[p]) mtext(fNameN, side=3, srt = 90, cex=1.0, line=1.8, col="red3")
          lines(dat_old[,v], lwd=0.4, col="orange")
        }
      }
    } # p in 1:pagesNum
  }   # for i in 1:nFiles
garbage <- dev.off()  # run dev.off() but suppress the notification it usually provides


###########################################################
## plot diff-data 'old - new'
###########################################################
print("<< creating diff plots >>")
pdf("figs/diff-prop_two-exp_new-minus-old.pdf", paper="a4", height=12.4, width=7.6, pointsize=12, 
    bg="white", family="Helvetica", version="1.6")
  for (i in 1:nFiles)
  {
    # read files and print filename
    print(flist_new[i])
    dat_new <- read.table(flist_new[i], header=T, sep="")
    print(flist_old[i])
    dat_old <- read.table(flist_old[i], header=T, sep="")
    # number of columns
    colNumN  <- ncol(dat_new)
    colNumO  <- ncol(dat_old)
    
    # get the file name without directory
    fNameN  <- gsub("new-output/", "", flist_new[i])
    fNameN  <- gsub(".txt", "", fNameN)

    # check whether column-names match (i.e., names & sequence of variables)
    # skip this file if column names of 'old' and 'new' are no match (using the R statement 'next')
    if (colNumN != colNumO)
    {
      print(paste0("   Number of variables differ between old and new output file: ", fNameN, " - skipping this file."))
      next
    }
    if ( all(colnames(dat_new) != colnames(dat_old)) )
    {
      print(paste0("   Names of variables differ between old and new output file: ", fNameN, " - skipping this file."))
      next
    }

    # calc proportional difference
    # show the proportional difference of the new to the old simulation, based on the old simulation data
    # calc absolute diff
    dat_diff <- dat_new - dat_old
    # set all 0 in dat_old == NA to avoid division by zero
    # info: in plot() the NA will be simply left out
    dat_old[dat_old == 0] <- NA
    # calc proportional diff
    dat_diff <- (dat_diff / dat_old) * 100

    # collect information to distribute figs over multiple pages if more than 14
    # the number of variables is the 'number of columns - 1', as the 1st column is always 'time'
    varNum   <- colNumN - 1
    maxFigPg <- 14
    if (varNum < (maxFigPg + 1))
    {
      # page numbers
      pagesNum <- 1
      # column numbers of output file per page 
      colNumStarts <- 2
      colNumEnds   <- colNumN
    } else{
      # page numbers
      pagesNum <- as.integer(varNum / maxFigPg)  # number of pages minus 1 | == pages with 14 figs each
                                                 # last page may has less figs
      # test if there is a last page with /= 14 figs, and increase pagesNum if true 
      if ((varNum / maxFigPg) > pagesNum)
      {
        pagesNum <- pagesNum + 1
      }
      # column numbers of output file per page 
      colNumLast   <- colNumN
      colNumStarts <- seq(2, colNumLast, maxFigPg)    # e.g., seq(2,37,14) = 2 16 30
      colNumEnds   <- colNumStarts[2:length(colNumStarts)] - 1
      colNumEnds   <- c(colNumEnds, colNumLast)      
    }

    # create figures per page
    for (p in 1:pagesNum)
    {
      # get the effective (actual) number of variables at this page
      effVarNum <- colNumEnds[p] - colNumStarts[p] + 1

      # page layout
      # see also layout information in below 'plot()' command
      # set the minimum of effVarNum to 4 (so that figures get not too big if their number is)
      mat  <- matrix(seq(1,2*round(max(effVarNum,4)/2)+2,1), nrow=round(max(effVarNum,4)/2)+1, ncol=2, byrow=T) 
      layout(mat)
      par(mai=c(0.32, 0.4, 0.14, 0.02), oma=c(0.1, 0.1, 2.0, 0.1)) # oma = outer margins  || c(bottom, left, top, right)

      # loop over columns (variables) and create figs   
      idInfo1 <- "diff: (new - old) / old * 100"
      for (v in colNumStarts[p]:colNumEnds[p])
      {
        # plot simulation ID info at uppermost figures per page
        idInfo <- ""
        if(v == colNumStarts[p])     idInfo <- idInfo1
        # set y-axis range to default values
        yMin <- -1
        yMax <-  1
        # test if all values are NA
        if(all(is.na(dat_diff[,v])))
        {
          print(paste0("  creating a green zero-line diff-figure for ", colnames(dat_diff)[v], " because all values are NA"))
          plot(1:nrow(dat_diff), rep(0, nrow(dat_diff)), 
               type='l',lwd=0.4, col="forestgreen", ylim=c(yMin,yMax), ylab=colnames(dat_diff)[v], xlab=outFPer, main=idInfo, 
               tck=0.03, mgp=c(1.5, 0.2, 0))
          if(v == colNumStarts[p]) mtext(fNameN, side=3, srt = 90, cex=1.0, line=1.8, col="red3")
        } else{
          # get useful min/max of y-axis
          yMin <- min(dat_diff[,v], na.rm = T)
          yMax <- max(dat_diff[,v], na.rm = T)
          if(!is.na(yMin)){ 
            if(yMin < -120.0) yMin <- -120.0
          } else yMin <- -1
          if(!is.na(yMax)){ 
            if(yMax > 120.0) yMax <-  120.0
          } else yMax <- 1
          # plot
          plot(dat_diff[,v], 
               type='l',lwd=0.4, col="black", ylim=c(yMin,yMax), ylab=colnames(dat_diff)[v], xlab=outFPer, main=idInfo, 
               tck=0.03, mgp=c(1.5, 0.2, 0))
          if(v == colNumStarts[p]) mtext(fNameN, side=3, srt = 90, cex=1.0, line=1.8, col="red3")
          abline(h=c(-4,4), lwd=0.8, col="red")
        }
      }
    } # p in 1:pagesNum
  }   # for i in 1:nFiles
garbage <- dev.off()  # run dev.off() but suppress the notification it usually provides

# message
print(" finished postprocessing comparison of two QUINCY experiments!")

