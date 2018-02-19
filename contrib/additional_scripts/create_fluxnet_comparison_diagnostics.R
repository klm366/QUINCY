#!/usr/bin/env Rscript
# 
# This script generates comparison diagnostics and plot for MODEL fluxnet simulations
# originally written by Lihui Luo and Soenke Zaehle
# initial version of this script is based on svn version r996 of generate_FLUXNET_diagnostics.R
# unused functions and code are removed

# contains functions:
# 
# MODEL_FLUXNET_GRAPHIC          - ordinary fluxnet plots
#
# MODEL_FLUXNET_COMPUTE_GROUPING - aggregate to daily, monthly, annual and diurnal cycle
# MODEL_FLUXNET_EVALUATION       - compute Root Mean Square Error, R sqaure, NSE, PBIAS, RSR


#################################################################################################################
# (1) Init                                                                                                      #
#################################################################################################################

# Libraries
library(ncdf4)
library(lattice)


# Set parameters
Args <- commandArgs(trailingOnly=TRUE);

startyear <- as.numeric(Args[1])
starttime <- paste(Args[1],"-01-01",sep='')
endyear   <- as.numeric(Args[2])

baselineFolder   <- Args[3]
simulationFolder <- Args[4]

siteName <- Args[5]

setwd(simulationFolder)


#################################################################################################################
# (2) Definition the graphic function                                                                           #
#################################################################################################################

# Graphic Function (1) FLUXNET STUFF-----------------------------------------------------------------
MODEL_FLUXNET_GRAPHIC <- function(MODEL_V, MODEL_2, FLUXNET_V, TIME, 
                                  mean, gtype, labs, variable_info, 
                                  RMSE, NRMSE, RR, NSE, PBIAS, RSR){

opar<-par( cex=0.6, pin=c(6.0, 2.0))

split.screen(c(2,1))
split.screen(c(1,2),screen=2)
screen(1)
  plot(TIME, MODEL_V*NA, type='n', 
       ylim=c(min(MODEL_V, MODEL_2, FLUXNET_V,na.rm=TRUE), max(MODEL_V, MODEL_2, FLUXNET_V,na.rm=TRUE)), 
       xlab="Time", ylab=paste(variable_info[1],"(", variable_info[2], ")", sep=""), 
       main=paste0(typestr, " ", variable_info[1], " (", variable_info[3],"-",variable_info[4], ")"),
       las=1, axes=F)
  axis(1, at=1:length(TIME), lab=labs,las=1)
  axis(2,las=1)
  box() 
  if (!all(is.na( FLUXNET_V ))) {
    lines(TIME, FLUXNET_V, type=gtype, col="black",pch=21,cex=0.8,lwd=0.9)
  } else {
    print("  All FLUXNET_V is NaN. Skip plotting.")
  }
  if (!all(is.na( MODEL_2 ))) {
    lines(TIME, MODEL_2, type=gtype, col="dodgerblue3",pch=4,cex=0.8,lwd=0.9)
  } else {
    print("  All MODEL_2 is NaN. Skip plotting.")
  }
  if (!all(is.na( MODEL_V ))) {
    lines(TIME, MODEL_V, type=gtype, col="red" ,pch=21,cex=0.8,lwd=0.9)
  } else {
    print("  All MODEL_V is NaN. Skip plotting.")
  }
  legend("topleft", legend = c("baseline","model","fluxnet"),  col=c("dodgerblue3","red","black"), pch=c(4,21,21), bty="n", ncol=1)
  legend("topright", legend = c(paste("NSE = ", NSE), paste("PBIAS = ", PBIAS, "%", sep=""), paste("RSR = ", RSR)), bty="n")

screen(3)
  if(!all(is.na(FLUXNET_V))) {
  if(!all(is.na(MODEL_V))) {
  if(!all(is.na(MODEL_2))) {
    plot(FLUXNET_V,MODEL_2, 
         xlab='fluxnet',ylab='model',las=1,
         col="dodgerblue3",pch=4,cex=0.8,
         ylim=c(min(MODEL_V, MODEL_2, FLUXNET_V,na.rm=TRUE), max(MODEL_V, MODEL_2, FLUXNET_V,na.rm=TRUE)),
         xlim=c(min(MODEL_V, MODEL_2, FLUXNET_V,na.rm=TRUE), max(MODEL_V, MODEL_2, FLUXNET_V,na.rm=TRUE)),
         main=paste0(variable_info[1]," (",variable_info[2],")"))
    points(FLUXNET_V, MODEL_V, col="red",pch=21,cex=0.8)
    line_FJ2 <- line(FLUXNET_V,MODEL_2)
    line_FJ  <- line(FLUXNET_V,MODEL_V)
    if (is.finite(coef(line_FJ2)[1]) & is.finite(coef(line_FJ2)[2])) {
      a=NA;b=NA
      if(!any(is.na(coef(line_FJ2)))) {
        a=sprintf("%.2f",coef(line_FJ2,use.na=false)[1])
        b=sprintf("%.2f",coef(line_FJ2,use.na=false)[2])
        abline(line_FJ2, col="royalblue4",lwd=0.9)
      }
    } else {
      print("  Coefficients of Baseline-Fluxnet have Inf values. Skip plotting regression line.")
    }
    if (is.finite(coef(line_FJ)[1]) & is.finite(coef(line_FJ)[2])) {
      a=NA;b=NA
      if(!any(is.na(coef(line_FJ)))) {
        a=sprintf("%.2f",coef(line_FJ,use.na=false)[1])
        b=sprintf("%.2f",coef(line_FJ,use.na=false)[2])
        abline(line_FJ, col="red3",lwd=0.9)
      }
      abline(0,1,col="gray50",lty=2)
      legend("bottomright", legend = c("model-fluxnet:",paste0("Int = ", a), paste0("slope = ", b), paste0("RMSE = ", RMSE), 
        paste0("NRMSE = ", NRMSE), as.expression(bquote(R^2 == .(RR)))), bty="n") 
    } else {
      print("  Coefficients of Model-Fluxnet have Inf values. Skip plotting regression line.")
    }
  }}}

screen(4)
  if(!all(is.na(FLUXNET_V))) {
  if(!all(is.na(MODEL_V))) {
  if(!all(is.na(MODEL_2))) {
    boxplot(MODEL_2, MODEL_V, FLUXNET_V, 
            names=c("baseline","model", "fluxnet"),las=1,
            ylab=variable_info[1],
            col=c("dodgerblue2","red",NA),
            main=paste0(variable_info[1]," (", variable_info[2], ")"))
  }}}
close.screen(1:4)
}


#################################################################################################################
# (3) Definition of an "aggregation by time" function                                                           #
#################################################################################################################

MODEL_FLUXNET_COMPUTE_GROUPING <- function(Value, TIME, period, method){

  if (period=="perhour"){
    time <- format(as.POSIXct(as.Date(as.numeric(TIME-1/24), origin=starttime)),"%m-%H-%M")
  }

  if (period=="perhour_15"){
    time <- format(as.POSIXct(as.Date(as.numeric(TIME-1/24), origin=starttime)),"%m-%d-%H-%M")
  }

  if (period=="daily"){
    time <- format(as.Date(as.numeric(TIME),origin=starttime), "%Y-%m-%d")
  }

  if (period=="halfmonth"){
    time <- format(as.Date(as.numeric(TIME),origin=starttime), "%Y-%m-%d")
  }

  if (period=="monthly"){
    time <- format(as.Date(as.numeric(TIME),origin=starttime), "%Y-%m")
  }
 
  if (period=="annual"){
    time <- format(as.Date(as.numeric(TIME),origin=starttime), "%Y")
  }

  JF_V <- data.frame(time, Value)
  JF_AGG <- aggregate(JF_V$Value, list(JF_V$time), method, na.rm=TRUE)

  if (period=="perhour_15"){
    time <- JF_AGG$Group.1
    time <- paste(substr(time,1,3), ifelse(as.integer(substr(time,4,5)) %in% c(1:15), 15, 30), substr(time,6,11),sep="")
    JF_V <- data.frame(time, Value=JF_AGG$x)
    JF_AGG <- aggregate(JF_V$Value, list(JF_V$time), method, na.rm=TRUE)
  }
  if (period=="halfmonth"){
    time <- JF_AGG$Group.1
    time <- paste(substr(time,1,8), ifelse(as.integer(substr(time,9,10)) %in% c(1:15), 15, 30),sep="")
    JF_V <- data.frame(time, Value=JF_AGG$x)
    JF_AGG <- aggregate(JF_V$Value, list(JF_V$time), method, na.rm=TRUE)
  }
  return (JF_AGG)

}


#################################################################################################################
# (4) Definition of a "model evaluation" function                                                               #
#################################################################################################################

MODEL_FLUXNET_EVALUATION <- function(method, MODEL_V, FLUXNET_V){
  
  d1=FLUXNET_V[!is.na(FLUXNET_V)&!is.na(MODEL_V)]
  d2=MODEL_V[!is.na(FLUXNET_V)&!is.na(MODEL_V)]

  if(length(d1)>0) {

   # Root Mean Square Error
   if (method=="rmse") r_value <- sqrt(mean((d1-d2)^2, na.rm=TRUE))

   # Normalized RMSE
   if (method=="nrmse") r_value <- (sqrt(mean((d1-d2)^2, na.rm=TRUE)))/mean(d1)

   # Coefficient of determination (R Squared)
   if (method=="RR") r_value <- summary(lm(d2 ~ d1, na.action=na.omit))[c("r.squared")]
  
   # Pearson's correlation coefficient
   if (method=="cor_test") r_value <- cor.test(d1, d2, na.action=na.omit )$p.value
     #cor.test(FLUXNET_V, MODEL_V, method="spearman")
     #cor.test(FLUXNET_V, MODEL_V, method="kendall")

   # Nash-Sutchliffe efficiency (NSE)
   if (method=="nse") r_value <- 1- sum((d1 - d2)^2)/sum((d1 - mean(d2))^2)

   # RMSE-observations standard deviation ratio (RSR)
   if (method=="rsr") r_value <- sqrt(sum((d1 - d2)^2))/sqrt(sum((d1 - mean(d2))^2)) 

   # Percent bias (PBIAS) 
   if (method=="pbias") r_value <- sum(d2 - d1)*100/sum(d1)

   # Modelling efficiency
   if (method=="mef") r_value <- (sum((d1-mean(d1))^2)-sum((d2-d1)^2))/sum((d1-mean(d1))^2)

   # Normalized average error
   if (method=="nae") r_value <- (mean(d2)-mean(d1))/mean(d1)
  
   # Variance Ratio
   if (method=="vr") r_value <- var(d2)/var(d1) 
    
  } else {
   r_value <- NA 
  }

  return(sprintf("%.2f", r_value))
}


#################################################################################################################
# (5) Here comes the main code ----------------------------------------------------------------------------------
#################################################################################################################

# define required compared variables, can be added later.
nvar <- 11
aggregation_level <- c("perhour", "daily", "monthly", "annual")   # possible options: perhour, perhour_15, daily, halfmonth, monthly, annual
stats_time        <- c("halfhourly", "daily", "monthly")
annual_values     <- array(NA,c(3,nvar+3))
variable_names <- array(c("PPFD", "PPFD", "PPFD", "PPFD",
                          "Rn", "Rn", "Net Radiation", "Rn",
                          "Qh", "H", "Sensible Heat Flux", "Qh" , 
                          "Qle", "LE", "Latent Heat Flux", "Qle",  
                          "Ga", "Ga", "Aerodynamic Conductance", "Ga" ,
                          "Gc", "Gc", "Canopy Conductance", "Gc" ,
                          "GPP", "GPP", "GPP", "GPP",
                          "Reco", "Reco", "Reco", "Reco",
                          "NEE", "NEE", "NEE", "NEE",
                          "WUE", "WUE", "Water-use Efficiency", "WUE",
                          "FAPAR", "fapar", "FAPAR", "fapar" )
                          , dim=c(4,nvar))
variable_units <- array(c("mol/m^2s", "umol/m^2s", "umol/m^2s",
                          "W/m^2", "W/m^2", "W/m^2", 
                          "W/m^2", "W/m^2", "W/m^2", 
                          "W/m^2", "W/m^2", "W/m^2", 
                          "m/s","m/s","m/s",
                          "m/s","m/s","m/s",
                          "mol/m^2s", "umol/m^2s", "umol/m^2s",
                          "mol/m^2s", "umol/m^2s", "umol/m^2s",
                          "mol/m^2s", "umol/m^2s", "umol/m^2s",
                          "g/kg", "g/kg", "g/kg",
                          "-","-","-")
                          , dim=c(3,nvar))
variable_units_conversion <- array(c(1,1,
                                     1,1,
                                     1,1,
                                     1,1,
                                     1,1,
                                     1,1,
                                     1,1,
                                     1,1,
                                     1,1,
                                     1,1,
                                     1,1)
                                     , dim=c(2,nvar))

# open/read input files
datafile   <- nc_open("./output/obs.nc")                               ## measured fluxnet data
modelfile  <- nc_open("./output/slm4fluxnet.nc")                       ## new/testing simulation results
modelfile2 <- nc_open(paste0(baselineFolder,"/output/slm4fluxnet.nc")) ## baseline simulation results for comparision

TIME   <- ncvar_get( modelfile, "time")
TIME_F <- ncvar_get( datafile, "time" )
TIME   <- seq(0,length(TIME)-1,1)/48

# diagnostic output is saved in the file: exp_*/site_folder/std_out/out.*.txt
print("                                                                 ")
print(" Create comparison Fluxnet diagnostics output.                   ")
print(" Output file is saved in:                                        ")
print(paste0(" ",simulationFolder,"/postprocessing"))
print("                                                                 ")

# create .pdf output file
pdf(paste0("./postprocessing/fluxnet_diagnostics_comp_",siteName,".pdf"))
## loop over variables to create figures from
for (V in 1:nvar){
     print(variable_names[3,V])
     MODEL_V <- ncvar_get( modelfile, variable_names[1,V] )
       MODEL_V <- MODEL_V*(variable_units_conversion[1,V])       
     MODEL_V[MODEL_V <= -9999.] <- NA
     MODEL_2 <- ncvar_get( modelfile2, variable_names[1,V] )
     MODEL_2 <- MODEL_2*(variable_units_conversion[1,V])       
     MODEL_2[MODEL_2 <= -9999.] <- NA

     # fapar
     if(variable_names[1,V]=='FAPAR') {
       FLUXNET_V <- ncvar_get( datafile, 'fapar_seawifs' )[1:length(TIME)]
       FLUXNET_V[FLUXNET_V <= -9999.] <- NA
       
       FLUXNET_V_FLAG <- ncvar_get( datafile, 'fapar_seawifs_flag' )[1:length(TIME)]
       FLUXNET_V_FILTERED <- FLUXNET_V
       FLUXNET_V_FILTERED[FLUXNET_V_FLAG > 0.5] <- NA
       
       FLUXNET_V2 <- ncvar_get( datafile, 'fapar_modis' )[1:length(TIME)]
       FLUXNET_V2[FLUXNET_V2 <= -9999.] <- NA
       
       FLUXNET_V2_FLAG <- ncvar_get( datafile, 'fapar_modis_flag' )[1:length(TIME)]
       FLUXNET_V2_FILTERED <- FLUXNET_V2
       FLUXNET_V2_FILTERED[FLUXNET_V2_FLAG > 0.5] <- NA

     # not fapar
     } else {
       FLUXNET_V <- ncvar_get( datafile, variable_names[2,V] )[1:length(TIME)]
       FLUXNET_V[FLUXNET_V <= -9999.] <- NA
       FLUXNET_V <- FLUXNET_V*(variable_units_conversion[2,V]) 
       FLUXNET_V_FILTERED <- FLUXNET_V
       if(!variable_names[2,V]=='Rn'&!variable_names[2,V]=='Gc'&!variable_names[2,V]=='Ga'&!variable_names[2,V]=='WUE') {
          FLUXNET_V_FLAG <- ncvar_get( datafile, paste(variable_names[2,V],'_flag',sep='') )[1:length(TIME)]
          FLUXNET_V_FILTERED[FLUXNET_V_FLAG > 1] <- NA
       }
     }

     # loop over "perhour", "daily", "monthly", "annual"
     for (time in aggregation_level ){ 
       if ( time != "annual" ) {
         
         # run this for all but "annual"
         variable_info       <- c(variable_names[3,V],variable_units[3,V],startyear,endyear)

         MODEL_AGG           <- MODEL_FLUXNET_COMPUTE_GROUPING(MODEL_V,TIME, time,'mean')
         MODEL_AGG2          <- MODEL_FLUXNET_COMPUTE_GROUPING(MODEL_2,TIME, time,'mean')
         MODEL_AGG           <- cbind(MODEL_AGG, MODEL_AGG2$x)
         colnames(MODEL_AGG) <- c("Group.1","x","x2")
         remove(MODEL_AGG2)
         
         FLUXNET_AGG   <- MODEL_FLUXNET_COMPUTE_GROUPING(FLUXNET_V,TIME, time,'mean')

         RMSE          <- MODEL_FLUXNET_EVALUATION("rmse", MODEL_AGG$x, FLUXNET_AGG$x)
         NRMSE         <- MODEL_FLUXNET_EVALUATION("nrmse", MODEL_AGG$x, FLUXNET_AGG$x)
         RR            <- MODEL_FLUXNET_EVALUATION("RR", MODEL_AGG$x, FLUXNET_AGG$x)
         MODEL_SD      <- (MODEL_FLUXNET_COMPUTE_GROUPING(MODEL_V,TIME, time,'sd'))$x
         MODEL_SD2     <- (MODEL_FLUXNET_COMPUTE_GROUPING(MODEL_2,TIME, time,'sd'))$x
         FLUXNET_SD    <- (MODEL_FLUXNET_COMPUTE_GROUPING(FLUXNET_V,TIME, time,'sd'))$x
         
         NSE           <- MODEL_FLUXNET_EVALUATION("nse", MODEL_AGG$x, FLUXNET_AGG$x)
         PBIAS         <- MODEL_FLUXNET_EVALUATION("pbias", MODEL_AGG$x, FLUXNET_AGG$x)
         RSR           <- MODEL_FLUXNET_EVALUATION("rsr", MODEL_AGG$x, FLUXNET_AGG$x)
         COR_TEST      <- MODEL_FLUXNET_EVALUATION("cor_test", MODEL_AGG$x, FLUXNET_AGG$x)

         ## additionally run for: daily, halfmonth, monthly, perhour, perhour_15
         if (time=="perhour") {
               typestr <- c("Per Hour Monthly")
               gtype   <- c("b")
	             labs    <- month.abb[as.numeric(gsub("0","",substr(MODEL_AGG$Group.1, 1, 2)))]
         }
         else if (time=="perhour_15") {
               typestr <- c("Per Hour 15 days")
               gtype   <- c("p")
	             labs    <- month.abb[as.numeric(gsub("0","",substr(MODEL_AGG$Group.1, 1, 2)))]
         }
         else if (time=="daily"){
               typestr <- c("Daily Mean")
               gtype   <- c("p")
	             labs    <- substr(MODEL_AGG$Group.1, 1, 4) 
         }
	       else if (time=="halfmonth"){
               typestr <- c("Half Month")
               gtype   <- c("b")
	             labs    <- substr(MODEL_AGG$Group.1, 1, 4)
         }
         else if (time=="monthly"){           
               typestr <- c("Monthly Mean")
               gtype   <- c("b")
               labs    <- substr(MODEL_AGG$Group.1, 1, 4)
         }


         ## create figures for each of the variable and time using the same function
         ## previously the below outcommented code was used to run diff graphics functions depending on fapar/other variable
         MODEL_FLUXNET_GRAPHIC(MODEL_AGG$x, MODEL_AGG$x2, FLUXNET_AGG$x, MODEL_AGG$Group.1, 
                               typestr, gtype, labs, variable_info, 
                               RMSE, NRMSE, RR, NSE, PBIAS, RSR)

       # end of if ( time != "annual" )
       # else here leads to time == "annual"
       } else {   
         annual_values[1,V]=mean(MODEL_V,na.rm=TRUE)
         annual_values[2,V]=mean(FLUXNET_V,na.rm=TRUE)
       } # end if/else (time != "annual")
     }   # end loop over aggregation_level ("perhour", "daily", "monthly", "annual")

} ## end loop over variables to create figures from: for (V in 1:nvar)

garbage <- dev.off()
