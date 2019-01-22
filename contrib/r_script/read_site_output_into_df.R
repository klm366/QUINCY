#################################################################################################################################
# this file defines a function that reads all/some of the default output files into one data.frame and returns this data.frame
# 
# you can define the set of files you want to read at line 80  (where the variable 'files' is defined)
# 
# in your R script simply source this file:
#   source("/path/repo/contrib/r_script/read_site_output_into_df.R")
# and execute this function:
#   myDF <- read_site_output_into_df()
#################################################################################################################################


## function definition ##
read_site_output_into_df <- function(){

  ## set experiment directory with site name
  ## you can also select a test_bed folder, in which case you may simply set siteName to ''
  expDir      <- '/Net/Groups/BGI/work_2/QUINCY/SLM/ReferenceRuns/fluxnet_diagnostics_default/'   # must end with a slash '.../'  !
  siteName    <- 'DE-Hai'   # for test_bed set to: ''
  
  ## set output interval
  outInterval <- 'monthly'  # timestep daily weekly monthly yearly spinup_yearly
  
  ## sets output files you could read
  f_atmforcing_ecopool <- (c(
                    'atm_forcing_',
                    'ecopoolCNP_'))
  f_myco <- (c(
                    'mycofluxCNP_',
                    'myc_diagnostics'))
  f_soil_mass_balance <- (c(   
                    'sb_massBalance_'))
  f_soil_flux <- (c(   
                    'soilfluxCdC13DC14_',
                    'soilfluxdN_',
                    'soilfluxN_',
                    'soilfluxP_'))
  f_soil_pool <- (c(   
                    'soilpooldC13DC14dN15_',
                    'soilpoolInOrg_',
                    'soilpoolOrg_',
                    'soilpoolPhysics_'))
  f_soil_lay <- (c(
                    'soil_diagnostics_lay_',
                    'soilfluxC_lay_',
                    'soilfluxN_lay_',
                    'soilfluxP_lay_',
                    'soilpoolCOrg_lay_',
                    'soilpooldC13org_lay_',
                    'soilpoolDC14org_lay_',
                    'soilpooldN15org_lay_',
                    'soilpoolInNOrg_lay_',
                    'soilpoolInPOrg_lay_',
                    'soilpoolNOrg_lay_',
                    'soilpoolPhysics_lay_',
                    'soilpoolPOrg_lay_'))
  f_veg_canopy_diagnostics <- (c(
                    'veg_canopy_',
                    'veg_diagnostics_'))
  f_veg_flux <- (c(
                    'vegfluxC_',
                    'vegfluxdC13DC14_',
                    'vegfluxEnergy_',
                    'vegfluxH2O_',
                    'vegfluxN_',
                    'vegfluxP_'))
  f_veg_pool <- (c(
                    'vegpoolC_',
                    'vegpooldC13DC14_',
                    'vegpooldN15_',
                    'vegpoolN_',
                    'vegpoolP_'))
  
  
  #########################################################
  ## compile the set of files you need in the df         ##
  ## simply comment in/out the ones you do/not need      ##
  files <- c(   f_atmforcing_ecopool,
                #f_myco,
                #f_soil_mass_balance,
                f_soil_flux,
                f_soil_pool,
                #f_soil_lay,
                #f_veg_canopy_diagnostics,
                f_veg_flux,
                f_veg_pool)
  ##                                                     ##
  #########################################################


  
  
  
  #########################################################
  ##  code you might not want to change                  ##
  ##                                                     ##
  
  ## default file type of output files
  fileType <- '.txt'
  
  # differentiate beetween test_bed and HPC simulations
  outDir <- ''
  if (nchar(siteName, type="chars") > 1){
    outDir <- 'output/'  
  }
  
  # read 1st file
  print(paste0('reading... ', files[1], outInterval, fileType))
  dafr   <- read.table(paste0(expDir, '/', siteName, '/', outDir, files[1], outInterval, fileType), sep="", header=T, stringsAsFactors=F)
  # read all other files
  for(f in 2:length(files)) {
    print(paste0('reading... ', files[f], outInterval, fileType))
    df2  <- read.table(paste0(expDir, '/', siteName, '/', outDir, files[f], outInterval, fileType), sep="", header=T, stringsAsFactors=F)
    df2  <- df2[,2:ncol(df2)]
    dafr <- cbind(dafr,df2)
  }
  # rm temp df
  remove(df2)
  
  # all good
  print("              ")
  print(" ok files read")
  print(" returning data.frame")

  # return the dafr
  return(dafr)
  
} # end function defintion #
