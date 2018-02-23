#!/usr/bin/env Rscript
# 
# This script generates some simple diagnostics and plot for MODEL fluxnet simulations
# originally written by Lihui Luo and Soenke Zaehle

# contains functions:
# 
# MODEL_FLUXNET_GRAPHIC          - ordinary fluxnet plots
#
# MODEL_FLUXNET_COMPUTE_GROUPING - aggregate to daily, monthly, annual and diurnal cycle
# MODEL_FLUXNET_EVALUATION       - compute Root Mean Square Error, R sqaure, NSE, PBIAS, RSR

# functions currently not used and commented out:
# 
# MODEL_FLUXNET_GRAPHIC_FAPAR    - satellite comparison plots
# MODEL_FLUXNET_GRAPHIC_FORCING  - compare forcing data with observation 
# MODEL_FLUXNET_GRAPHIC_ERR      - error bars 
# MODEL_FLUXNET_STATS_GRAPHIC    

#################################################################################################################
# (1) Init                                                                                                      #
#################################################################################################################

# Set parameters
Args <- commandArgs(trailingOnly=TRUE);
if(length(Args) != 2) {
 message("requires startyear endyear as input (e.g. generate_FLUXNET_diagnostics.R 2000 2006). Terminating");quit()
}


#################################################################################################################
# (2) Definition of five graphic functions                                                                      #
#################################################################################################################

# Graphic Function (1) FLUXNET STUFF-----------------------------------------------------------------
MODEL_FLUXNET_GRAPHIC <- function(MODEL_V, FLUXNET_V, TIME, 
                                  mean, gtype, labs, variable_info, 
                                  RMSE, NRMSE, RR, NSE, PBIAS, RSR){

opar<-par( cex=0.6, pin=c(6.0, 2.0))

split.screen(c(2,1))
split.screen(c(1,2),screen=2)
screen(1)
  plot(TIME, MODEL_V*NA, type='n',
    ylim=c(min(MODEL_V, FLUXNET_V,na.rm=TRUE), max(MODEL_V, FLUXNET_V,na.rm=TRUE)), 
    xlab="Time", ylab=paste(variable_info[1],"(", variable_info[2], ")", sep=""), 
    main=paste0(typestr, " ", variable_info[1], " (", variable_info[3],"-",variable_info[4], ")"),
    las=1, axes=F)
  axis(1, at=1:length(TIME), lab=labs,las=1)
  axis(2,las=1)
  box() 
  if (!all(is.na( FLUXNET_V ))) {
    lines(TIME, FLUXNET_V, type=gtype, col="black",pch=21,cex=0.8)
  } else {
    print("  All FLUXNET_V is NaN. Skip plotting.")
  }
  if (!all(is.na( MODEL_V ))) {
    lines(TIME, MODEL_V, type=gtype, col="red" ,pch=21,cex=0.8)
  } else {
    print("  All MODEL_V is NaN. Skip plotting.")
  }
  legend("topleft", legend = c("model","fluxnet"),  col=c("red","black"), pch=c(21,21), bty="n", ncol=1)
  legend("topright", legend = c(paste("NSE = ", NSE), paste("PBIAS = ", PBIAS, "%", sep=""), paste("RSR = ", RSR)), bty="n")

screen(3)
  if(!all(is.na(FLUXNET_V))) {
  if(!all(is.na(MODEL_V))) {
    plot(FLUXNET_V,MODEL_V, xlab='fluxnet',ylab='model',las=1,
         col="red",pch=21,cex=0.8,
         ylim=c(min(MODEL_V, FLUXNET_V,na.rm=TRUE), max(MODEL_V, FLUXNET_V,na.rm=TRUE)),
         xlim=c(min(MODEL_V, FLUXNET_V,na.rm=TRUE), max(MODEL_V, FLUXNET_V,na.rm=TRUE)),
         main=paste0(variable_info[1]," (",variable_info[2],")"))
    line_FJ <- line(FLUXNET_V,MODEL_V) 
    if (is.finite(coef(line_FJ)[1]) & is.finite(coef(line_FJ)[2])) {
      a=NA;b=NA
      if(!any(is.na(coef(line_FJ)))) {
        a=sprintf("%.2f",coef(line_FJ,use.na=false)[1])
        b=sprintf("%.2f",coef(line_FJ,use.na=false)[2])
        abline(line_FJ, col="black",lwd=0.9)
      }
      abline(0,1,col="gray50",lty=2)
      legend("bottomright", legend = c(paste("Int = ", a), paste("slope = ", b), paste("RMSE = ", RMSE), 
        paste("NRMSE = ", NRMSE), as.expression(bquote(R^2 == .(RR)))), bty="n") 
    } else {
      print("  Coefficients have Inf values. Skip plotting regression line.")
    }
  }
}
screen(4)
  if(!all(is.na(FLUXNET_V))) {
  if(!all(is.na(MODEL_V))) {
    boxplot(MODEL_V, FLUXNET_V, names=c("model", "fluxnet"),las=1,
            ylab=variable_info[1], 
            main=paste0(variable_info[1]," (", variable_info[2], ")"))
  }}
close.screen(1:4)
}

# Graphic Function (2) SATELLITE STUFF-----------------------------------------------------------------
# MODEL_FLUXNET_GRAPHIC_FAPAR <- function(MODEL_V, FLUXNET_V, FLUXNET_V2, FLUXNET_V3, TIME, typestr, gtype, variable_info, 
#   RMSE_seawifs, RMSE_modis, RR_seawifs, RR_modis, NSE, PBIAS, RSR){
# 
# opar<-par( cex=0.6, pin=c(6.0, 2.0))
# 
# split.screen(c(2,1))
# split.screen(c(1,2),screen=2)
# screen(1)
#   plot(TIME, MODEL_V, type=gtype,  col="red", pch=19,
#     ylim=c(min(MODEL_V, c(FLUXNET_V,FLUXNET_V2,FLUXNET_V3),na.rm=TRUE), 
#            max(MODEL_V, c(FLUXNET_V,FLUXNET_V2,FLUXNET_V3),na.rm=TRUE)), 
#     xlab="Time", ylab=paste(variable_info[1],"(", variable_info[2], ")", sep=""), 
#     main=paste(typestr, " ", variable_info[1], " (", variable_info[3],"-",variable_info[4], ")", sep=""),
#     las=1)
#   lines(TIME, FLUXNET_V, type=gtype, col="black" )
#   lines(TIME, FLUXNET_V2, type=gtype, col="blue" )
#   lines(TIME, MODEL_V, type=gtype, col="red" ,pch=19)
#   legend("bottomright", legend = c("MODEL","SeaWifs","Modis","Cyclopes"),  
#     col=c("red","black","blue","grey"), pch=c(19,21,21,21), bty="n", ncol=4)
#   legend("topright", legend = c(paste("NSE = ", NSE), paste("PBIAS = ", PBIAS, "%", sep=""), paste("RSR = ", RSR)), bty="n")
# 
# screen(3)
#   plot(FLUXNET_V,MODEL_V, xlab='SeaWifs',ylab='MODEL',las=1,xlim=c(0,1),ylim=c(0,1),
#       main= paste(variable_info[1], "(" ,  variable_info[2], ")", sep= "" ))
#   legend("topleft", legend = c(paste("RMSE (SeaWifs) = ", RMSE_seawifs), paste("RMSE (Modis) = ", RMSE_modis),
#     as.expression(bquote(R^2 (SeaWifs) == .(RR_seawifs))), as.expression(bquote(R^2 (Modis) == .(RR_modis)))), bty="n")
#   points(FLUXNET_V2, MODEL_V,col='blue',pch=1)
# 
# screen(4)
#   boxplot(MODEL_V, FLUXNET_V, names=c("MODEL", "FLUXNET"), col = "light grey", las=1,
#     ylab=paste(variable_info[1],"(", variable_info[2], ")", sep=""), 
#     main = paste("Boxplot of ", variable_info[1]))
# 
# close.screen(1:4)
# }


# Graphic Function (3) compare Forcing data with simulation and evaluation------------------------
# MODEL_FLUXNET_GRAPHIC_FORCING <- function(MODEL_V, FLUXNET_V, FORCING_V, Parameter_info, TIME, typestr, gtype, variable_info){
# 
#   opar<-par( cex=0.6, pin=c(6.0, 2.0))
#   plot(TIME, MODEL_V, type="n",  col="green", 
#     ylim=c(min(MODEL_V, FLUXNET_V,na.rm=TRUE), max(MODEL_V, FLUXNET_V,na.rm=TRUE)), 
#     xlab="Time", ylab=paste(variable_info[1],"(", variable_info[2], ")", sep=""), 
#     main=paste(typestr, " ", variable_info[1], " (", variable_info[3],"-",variable_info[4], ")", sep=""),
#     las=1)
#   lines(TIME, MODEL_V, lwd=2, xlab='Time', ylab=paste(variable_info[1],"(", variable_info[2], ")", sep=""), ylim=c(min(MODEL_V, FLUXNET_V,na.rm=TRUE), max(MODEL_V, FLUXNET_V,na.rm=TRUE)), col="green")
#   lines(TIME, FLUXNET_V, type=gtype, col="blue")
#  
#   par(new=T)
#   plot(TIME, FORCING_V, col='black', lwd=2, ann=F, axes=F)
#   mtext(Parameter_info, side=4, line=3, col='black')
#   axis(4, col.axis='black', col='black')
#   legend("topleft", legend = c("MODEL","FLUXNET","FORCING"),  col=c("green","blue","black"), pch=19, bty="n", ncol=2)
# }


# Graphic Function (4) Error Bars----------------------------------------------------------------------
# MODEL_FLUXNET_GRAPHIC_ERR <- function(MODEL_V, FLUXNET_V, TIME, MODEL_SD, FLUXNET_SD, variable_info, typestr){
# 
# opar<-par( cex=0.6, pin=c(6.0, 2.0))
# 
# split.screen(c(2,1))
# 
# screen(1)
#   plotCI(1:length(TIME), MODEL_V, MODEL_SD, pch=21, col="red",sfrac=0.001,gap=0,
#     ylim=c(min(MODEL_V-MODEL_SD, FLUXNET_V-FLUXNET_SD, na.rm=TRUE),
#            max(MODEL_V+MODEL_SD, FLUXNET_V+FLUXNET_SD, na.rm=TRUE)),
#     xlab="Time", ylab=paste(variable_info[1],"(", variable_info[2], ")", sep=""),
#     main=paste("Standard Deviation of ", typestr, " ", variable_info[1], " (", variable_info[3],"-",variable_info[4], ")", sep=""))
# 
# 
#   plotCI(1:length(TIME),FLUXNET_V, FLUXNET_SD, pch=21,col="black",sfrac=0.001,gap=0, add=TRUE,
#     ylim=c(min(MODEL_V-MODEL_SD, FLUXNET_V-FLUXNET_SD, na.rm=TRUE),
#            max(MODEL_V+MODEL_SD, FLUXNET_V+FLUXNET_SD, na.rm=TRUE)),
#     xlab="Time", ylab=paste(variable_info[1],"(", variable_info[2], ")", sep=""),
#     main=paste("Standard Deviation of ", typestr, " ", variable_info[1], " (", variable_info[3],"-",variable_info[4], ")", sep=""))
#     #lines(1:length(TIME), FLUXNET_V, type="l", col="light grey")
#     legend("topleft", legend = c("MODEL","FLUXNET"),
#     col=c("red","black"), pch=c(21,21), bty="n", ncol=2)
# 
# screen(2)
#   if (typestr=="Per Hour Monthly"){
#      TIME <- substr(TIME, 4, 5)
#      at <- 1:24
#      labs <- 1:24
#   }
#   if (typestr=="Half Month"){
#      TIME <- substr(TIME, 6, 10)
#      at <- 1:24
#      labs <-NA
#      for (V in 1:12) labs<-c(labs,rep(month.abb[V],2))
#      labs <- labs[2:25] 
#   }
#   if (typestr=="Monthly Mean"){
#      TIME <- substr(TIME, 6, 7)
#      at <- 1:12
#      labs <- month.abb
#   }
#   if (typestr=="Per Hour Monthly" || typestr=="Half Month" || typestr=="Monthly Mean"){
#     MODEL_V <- data.frame(TIME,MODEL_V)
#     FLUXNET_V <- data.frame(TIME,FLUXNET_V)
#     boxplot(MODEL_V ~ TIME, MODEL_V, col = "tomato", boxwex = 0.25, at = at - 0.2, 
#       ylab=paste(variable_info[1],"(", variable_info[2], ")", sep=""), 
#       main = paste(typestr, variable_info[1]),axes=F)
#     axis(1, at=at, lab=labs)
#     axis(2)
#     box()   
#     boxplot(FLUXNET_V ~ TIME, FLUXNET_V, add= TRUE, col = "orange", boxwex = 0.25, at = at + 0.2, axes=F,
#       ylab=paste(variable_info[1],"(", variable_info[2], ")", sep=""), 
#       main = paste(typestr, variable_info[1])) 
#     smartlegend(x="left",y="top", inset = 0, c("MODEL", "FLUXNET"), fill = c("tomato", "orange"), bty="n")
#   }
# close.screen(all = TRUE)
# }


# Graphic Function (5) Variables Statistics----------------------------------------------------------------------
# MODEL_FLUXNET_STATS_GRAPHIC <- function(stat_v, stat_name){
#   def.par <- par(no.readonly = TRUE)
#   nf <- layout(matrix(c(1,2,3,4,5,6),3,2,byrow=TRUE), heights=c(1,1,1), widths=c(3,1))
#   
#   for (V in 1:3){
# 	
#     ymin <- min(as.numeric(as.vector(eval_stats[stat_v[V],,])), na.rm=TRUE)
#     ymax <- max(as.numeric(as.vector(eval_stats[stat_v[V],,])), na.rm=TRUE)
#     plot(1:nvar, eval_stats[stat_v[V],,1], ylim=c(ymin, ymax), xlab="", ylab=stat_name[V], main="Variables Statistics", type="o", col="red", xaxt = "n", yaxt="n")
#     axis(1, at=1:nvar, lab=variable_names[4,1:nvar])
#     axis(2)
#     box()
#     lines(1:nvar, eval_stats[stat_v[V],,2], type="o", col="blue")
#     lines(1:nvar, eval_stats[stat_v[V],,3], type="o", col="green")
#     lines(1:nvar, eval_stats[stat_v[V],,4], type="o", col="orange")
#     legend("topright", legend = stat_name, col=c("red", "blue", "green", "orange"), pch=21, bty="n")
#     boxplot(as.numeric(as.vector(eval_stats[stat_v[V],,])))	
#   }
#   par(def.par)
# }


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

startyear <- as.numeric(Args[1])
starttime <- paste(Args[1],"-01-01",sep='')
endyear   <- as.numeric(Args[2])
work_dir  <- try(system("pwd",intern=TRUE))

setwd(work_dir)
library(ncdf4)
library(lattice)

# define required compared variables, can be added later.
nvar <- 11
aggregation_level <- c("perhour", "daily", "monthly", "annual")   # possible options: perhour, perhour_15, daily, halfmonth, monthly, annual
stats_time        <- c("halfhourly", "daily", "monthly")
annual_values     <- array(NA,c(3,nvar+3))
eval_stats        <- array(NA, dim=c(14,nvar,4))
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

# open input files
datafile  <- nc_open("obs.nc")            ## measured fluxnet data
modelfile <- nc_open("slm4fluxnet.nc")    ## simulation results

TIME   <- ncvar_get( modelfile, "time")
TIME_F <- ncvar_get( datafile, "time" )
TIME   <- seq(0,length(TIME)-1,1)/48

# diagnostic output is saved in the file: exp_*/site_folder/std_out/out.*.txt
print("                                                                 ")
print("-----------------------------------------------------------------")
print(" Diagnostic output of the script: generate_FLUXNET_diagnostics.R ")
print("                                                                 ")

# create .pdf output file
pdf('FLUXNET_diagnostics.pdf')
## loop over variables to create figures from
for (V in 1:nvar){
     print(variable_names[3,V])
     MODEL_V <- ncvar_get( modelfile, variable_names[1,V] )
     MODEL_V <- MODEL_V*(variable_units_conversion[1,V])       
     MODEL_V[MODEL_V <= -9999.] <- NA

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
       
       if( any(!is.na(FLUXNET_V_FILTERED)) & any(!is.na(MODEL_V)) ) {

         dummy <- MODEL_V
         dummy[is.na(FLUXNET_V_FILTERED)] <- NA

         eval_stats[1,V,1]=coef(line(FLUXNET_V_FILTERED[!is.na(FLUXNET_V_FILTERED)], MODEL_V[!is.na(FLUXNET_V_FILTERED)]), use.na=false)[1]
         eval_stats[2,V,1]=coef(line(FLUXNET_V_FILTERED[!is.na(FLUXNET_V_FILTERED)], MODEL_V[!is.na(FLUXNET_V_FILTERED)]), use.na=false)[2]
         eval_stats[3,V,1]=cor(MODEL_V, FLUXNET_V_FILTERED,use="na.or.complete")
         eval_stats[4,V,1]=MODEL_FLUXNET_EVALUATION("RR", MODEL_V, FLUXNET_V_FILTERED)
         eval_stats[5,V,1]=MODEL_FLUXNET_EVALUATION("nse", MODEL_V, FLUXNET_V_FILTERED)
         eval_stats[6,V,1]=MODEL_FLUXNET_EVALUATION("rmse", MODEL_V, FLUXNET_V_FILTERED)
         eval_stats[7,V,1]=MODEL_FLUXNET_EVALUATION("nrmse", MODEL_V, FLUXNET_V_FILTERED)
         eval_stats[8,V,1]=MODEL_FLUXNET_EVALUATION("pbias", MODEL_V, FLUXNET_V_FILTERED)
         eval_stats[9,V,1]=sd(dummy,na.rm=TRUE)
         eval_stats[10,V,1]=sd(FLUXNET_V_FILTERED,na.rm=TRUE)
         eval_stats[11,V,1]=MODEL_FLUXNET_EVALUATION("rsr", MODEL_V, FLUXNET_V_FILTERED)
         eval_stats[12,V,1]=MODEL_FLUXNET_EVALUATION("mef", MODEL_V, FLUXNET_V_FILTERED)
         eval_stats[13,V,1]=MODEL_FLUXNET_EVALUATION("nae", MODEL_V, FLUXNET_V_FILTERED)
         eval_stats[14,V,1]=MODEL_FLUXNET_EVALUATION("vr", MODEL_V, FLUXNET_V_FILTERED) 
       } 
     }

     # loop over "perhour", "daily", "monthly", "annual"
     si <- 1  # addressing the 3rd dimension of eval_stats[] starting with 2, as the '1' was already used above
     for (time in aggregation_level ){ 
       if ( time != "annual" ) {
         
         # run this for all but "annual"
         variable_info <- c(variable_names[3,V],variable_units[3,V],startyear,endyear)
         MODEL_AGG     <- MODEL_FLUXNET_COMPUTE_GROUPING(MODEL_V,TIME, time,'mean')
         FLUXNET_AGG   <- MODEL_FLUXNET_COMPUTE_GROUPING(FLUXNET_V,TIME, time,'mean')

         RMSE          <- MODEL_FLUXNET_EVALUATION("rmse", MODEL_AGG$x, FLUXNET_AGG$x)
         NRMSE         <- MODEL_FLUXNET_EVALUATION("nrmse", MODEL_AGG$x, FLUXNET_AGG$x)
         RR            <- MODEL_FLUXNET_EVALUATION("RR", MODEL_AGG$x, FLUXNET_AGG$x)
         MODEL_SD      <- (MODEL_FLUXNET_COMPUTE_GROUPING(MODEL_V,TIME, time,'sd'))$x
         FLUXNET_SD    <- (MODEL_FLUXNET_COMPUTE_GROUPING(FLUXNET_V,TIME, time,'sd'))$x
         
         NSE           <- MODEL_FLUXNET_EVALUATION("nse", MODEL_AGG$x, FLUXNET_AGG$x)
         PBIAS         <- MODEL_FLUXNET_EVALUATION("pbias", MODEL_AGG$x, FLUXNET_AGG$x)
         RSR           <- MODEL_FLUXNET_EVALUATION("rsr", MODEL_AGG$x, FLUXNET_AGG$x)
         COR_TEST      <- MODEL_FLUXNET_EVALUATION("cor_test", MODEL_AGG$x, FLUXNET_AGG$x)

         ## run for: daily, halfmonth, monthly
         if (time == "daily" || time == "halfmonth" || time == "monthly"){

          if ( !all(is.na(FLUXNET_AGG$x)) ) {
          if ( !all(is.na(MODEL_AGG$x)) ) {
   
             si <- si+1
             dummy <- MODEL_AGG$x
             dummy[is.na(FLUXNET_AGG$x)] <- NA

             eval_stats[1,V,si]=coef(line(FLUXNET_AGG$x, MODEL_AGG$x),use.na=false)[1]
             eval_stats[2,V,si]=coef(line(FLUXNET_AGG$x, MODEL_AGG$x),use.na=false)[2]
             eval_stats[3,V,si]=cor(MODEL_AGG$x, FLUXNET_AGG$x,use="na.or.complete")
             eval_stats[4,V,si]=MODEL_FLUXNET_EVALUATION("RR", MODEL_AGG$x, FLUXNET_AGG$x)
             eval_stats[5,V,si]=MODEL_FLUXNET_EVALUATION("nse", MODEL_AGG$x, FLUXNET_AGG$x)
             eval_stats[6,V,si]=MODEL_FLUXNET_EVALUATION("rmse",MODEL_AGG$x, FLUXNET_AGG$x )
  	         eval_stats[7,V,si]=MODEL_FLUXNET_EVALUATION("nrmse",MODEL_AGG$x, FLUXNET_AGG$x) 
             eval_stats[8,V,si]=MODEL_FLUXNET_EVALUATION("pbias", MODEL_AGG$x, FLUXNET_AGG$x)
             eval_stats[9,V,si]=sd(dummy,na.rm=TRUE)
             eval_stats[10,V,si]=sd(FLUXNET_AGG$x,na.rm=TRUE)
             eval_stats[11,V,si]=MODEL_FLUXNET_EVALUATION("rsr", MODEL_AGG$x, FLUXNET_AGG$x)
             eval_stats[12,V,si]=MODEL_FLUXNET_EVALUATION("mef", MODEL_AGG$x, FLUXNET_AGG$x)
             eval_stats[13,V,si]=MODEL_FLUXNET_EVALUATION("nae", MODEL_AGG$x, FLUXNET_AGG$x)
             eval_stats[14,V,si]=MODEL_FLUXNET_EVALUATION("vr", MODEL_AGG$x, FLUXNET_AGG$x) 
          }} # test NA 
         }
         
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
         MODEL_FLUXNET_GRAPHIC(MODEL_AGG$x, FLUXNET_AGG$x, MODEL_AGG$Group.1, 
                               typestr, gtype, labs, variable_info, 
                               RMSE, NRMSE, RR, NSE, PBIAS, RSR)

          #  ## run calc & graphics function MODEL_FLUXNET_GRAPHIC_FAPAR -> comparison against MODIS data
          #  if(variable_names[1,V]=="FAPAR"){
          #    if (time != "perhour" ){
          #      FLUXNET_AGG2 <- MODEL_FLUXNET_COMPUTE_GROUPING(FLUXNET_V2,TIME, time,'mean')
          # 
          #      RMSE_seawifs <- MODEL_FLUXNET_EVALUATION("rmse", MODEL_AGG$x, FLUXNET_AGG$x)
          #      RMSE_modis   <- MODEL_FLUXNET_EVALUATION("rmse", MODEL_AGG$x, FLUXNET_AGG2$x)
          #      RR_seawifs   <- MODEL_FLUXNET_EVALUATION("RR", MODEL_AGG$x, FLUXNET_AGG$x)
          #      RR_modis     <- MODEL_FLUXNET_EVALUATION("RR", MODEL_AGG$x, FLUXNET_AGG2$x)
          # 
          #      MODEL_FLUXNET_GRAPHIC_FAPAR(MODEL_AGG$x, FLUXNET_AGG$x, FLUXNET_AGG2$x, FLUXNET_AGG3$x, MODEL_AGG$Group.1, 
          #                                  typestr, gtype, variable_info, 
          #                                  RMSE_seawifs, RMSE_modis, RR_seawifs, RR_modis, NSE, PBIAS, RSR)
          #    }
          #  ## run graphics function MODEL_FLUXNET_GRAPHIC if not variable fapar
          #  } else{
          #      MODEL_FLUXNET_GRAPHIC(MODEL_AGG$x, FLUXNET_AGG$x, MODEL_AGG$Group.1, 
          #                            typestr, gtype, labs, variable_info, 
          #                            RMSE, NRMSE, RR, NSE, PBIAS, RSR)
          #  }

       # end of if ( time != "annual" )
       # else here leads to time == "annual"
       } else {   
         annual_values[1,V]=mean(MODEL_V,na.rm=TRUE)
         annual_values[2,V]=mean(FLUXNET_V,na.rm=TRUE)
       } # end if/else (time != "annual")
     }   # end loop over aggregation_level ("perhour", "daily", "monthly", "annual")

} ## end loop over variables to create figures from: for (V in 1:nvar)


## write tables with 'anual values' and 'eval_stats'
rownames                    <- c(variable_names[4,],"vegC","litterC","soilC")
annual_values[3,1:nvar]     <- variable_units[3,1:nvar]
annual_values[3,nvar+(1:3)] <- "gC/m2"

write.table(t(annual_values),file="summary.annual.txt",
            col.names=c('MODEL','Observed','Units'),
            row.names=rownames) 

for (s in 1:length(stats_time)){
  write.table(t(eval_stats[,,s]), file=paste0("statistics.",stats_time[s],".txt"),
              col.names=c("Intercept","Slope", "R", "R^2", "NSE", "RMSE", "NRMSE", "PBIAS", "SD_MOD", "SD_FLX", "RSR", "MEF", "NAE", "VR" ), 
              row.names=variable_names[4,])
}

# MODEL_FLUXNET_STATS_GRAPHIC(c(5,8,11),c("NSE","PBIAS","RSR"))

garbage <- dev.off()
