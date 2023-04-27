#summarise gridded ClimBRA daily files to monthly
#for gridded data from https://doi.org/10.57760/sciencedb.02316

library(terra)
library(logger)

path <- "/home/james/Documents/data/"
logfile <- "ClimBRA-summary.log"
log_appender(appender_tee(logfile))

#set model, scenarios to summarise
model <- c("EC-EARTH3")
scenario <- c("ssp585")
log_info(model)
log_info(scenario)

for(m in model){
  for(s in scenario){
    
    #precipitation with sum
    fpin <- paste0(path,model,"-pr-",scenario,".nc")
    
    log_info(paste0("reading ",fpin))
    dat <- rast(fpin)
    
    dat <- subset(dat,1:365)
    
    log_info("starting summary (year-months sum)")
    datym <- tapp(dat, "yearmonths", sum)

    fpout <- paste0(path,model,"-pr-",scenario,"-months.nc" )
    log_info(paste0("writing ", fpout))
    writeCDF(datym, fpout)
    
    #temperature variables with mean
    for(tas in c("tasmax", "tasmin")){
      
      fpin <- paste0(path,model,"-",tas,"-",scenario,".nc")
      log_info(paste0("reading ",fpin))
      dat <- rast(fpin)
      
      dat <- subset(dat,1:365)
      
      log_info("starting summary (year-months mean)")
      datym <- tapp(dat, "yearmonths", mean)

      fpout <- paste0(path,model,"-",tas,"-",scenario,"-months.nc" )
      log_info(paste0("writing ", fpout))
      writeCDF(datym, fpout)
    }
  }
}

