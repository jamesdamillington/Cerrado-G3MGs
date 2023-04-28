#summarise gridded ClimBRA daily files to monthly
#for gridded data from https://doi.org/10.57760/sciencedb.02316

library(terra)
library(lubridate)  #for ym() and ceiling_date()
library(logger)

path <- "/home/james/Documents/data/"
logfile <- "ClimBRA-summary.log"
log_appender(appender_tee(logfile))

#set model, scenarios, start, end to summarise
model <- c("EC-EARTH3")
scenario <- c("ssp585")
startym <- "2020-01"
endym <- "2020-12"

for(m in model){
  for(s in scenario){
    
    #precipitation with sum
    fpin <- paste0(path,model,"-pr-",scenario,".nc")
    
    log_info(paste0("reading ",fpin))
    dat <- rast(fpin)
    
    log_info("subset")
    dat <- subset(dat,time(dat) >= ym(startym) & time(dat) < ceiling_date(ym(endym), unit="months"))
    
    log_info("starting summary (year-months sum)")
    datym <- tapp(dat, "yearmonths", sum, na.rm=T) #na.rm to handle 29 Feb missing data

    #terra cannot write year-month to nc as it uses GDAL engine
    #terra also does not preserve layer names (unless they are defined as sub-datasets)
    #also zvals must be values (cannot be times)
    #so best we can do is save as much info in globals (and filename) and 
    #then set time when reading in elsewhere (e.g. moistureMap_future_cerrado.r) 
    fpout <- paste0(path,model,"-pr-",scenario,"-totalmonth-",startym,"-to-",endym,".nc")
    log_info(paste0("writing ", fpout))
    writeCDF(datym, fpout, overwrite=T,
             varname="pr", unit="mm", 
             longname=paste0("total monthly precipitation ",startym," to ",endym))
    
    #temperature variables with mean
    for(tas in c("tasmax", "tasmin")){
      
      fpin <- paste0(path,model,"-",tas,"-",scenario,".nc")
      log_info(paste0("reading ",fpin))
      dat <- rast(fpin)
      
      log_info("subset")
      dat <- subset(dat, time(dat) >= ym(startym) & time(dat) < ceiling_date(ym(endym), unit="months"))

      log_info("starting summary (year-months mean)")
      datym <- tapp(dat, "yearmonths", mean, na.rm=T)  #na.rm to handle 29 Feb missing data

      fpout <- paste0(path,model,"-",tas,"-",scenario,"-meanmonth-",startym,"-to-",endym,".nc")
      log_info(paste0("writing ", fpout))
      writeCDF(datym, fpout, overwrite=T,
               varname=tas, unit="degC", 
               longname=paste0("mean monthly ", tas,", ",startym," to ",endym))
    }
  }
}

